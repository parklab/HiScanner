"""
Helpers for multi-cell BAM input (e.g. 10x Genomics), where all cells share one
BAM and each read carries a cell barcode tag (CB by default).

The BAM is streamed one chromosome at a time and never split into per-cell BAMs.
Outputs use the same per-cell file layout as the per-cell-BAM mode, so all
downstream steps are unchanged.
"""
import gzip
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, Union

import pandas as pd
import pysam

from ..logger import logger

# Same filters as the per-cell `samtools view -q 30 -F 1284` (unmapped, secondary, duplicate)
READPOS_MIN_MAPQ = 30
READPOS_EXCLUDE_FLAGS = 1284

# Allele counting mirrors SCAN2's hc_raw.mmq60 input: MAPQ >= 60, and additionally
# skips supplementary (2048) and QC-fail (512) reads
ALLELE_MIN_MAPQ = 60
ALLELE_EXCLUDE_FLAGS = 1284 | 512 | 2048
ALLELE_MIN_BASEQ = 20

# Number of buffered read positions (over all cells) before flushing to disk
FLUSH_EVERY = 1_000_000


def is_multi_cell(config: Dict) -> bool:
    """Return True if the pipeline is configured for a shared multi-cell BAM."""
    return config.get('input_mode', 'per_cell_bam') == 'multi_cell_bam'


def find_bam_index(bam_path: Union[str, Path]) -> Optional[Path]:
    """Return the path of an existing .bai/.csi index for a BAM, or None."""
    bam_path = Path(bam_path)
    candidates = [
        Path(f"{bam_path}.bai"),
        Path(f"{bam_path}.csi"),
        bam_path.with_suffix('.bai'),
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def load_barcode_map(metadata_path: Union[str, Path]) -> Dict[str, str]:
    """
    Map barcode -> bamID for single cells in a multi-cell metadata file.

    The metadata must contain the columns bamID, barcode and singlecell.
    """
    metadata = pd.read_csv(metadata_path, sep='\t', dtype=str)
    missing = {'bamID', 'barcode', 'singlecell'} - set(metadata.columns)
    if missing:
        raise ValueError(
            f"Metadata {metadata_path} is missing column(s) required for "
            f"multi_cell_bam mode: {', '.join(sorted(missing))}"
        )
    cells = metadata[metadata.singlecell == 'Y']
    duplicated = cells.barcode[cells.barcode.duplicated()].unique()
    if len(duplicated) > 0:
        raise ValueError(
            f"Duplicated barcodes in metadata: {', '.join(duplicated[:10])}"
        )
    duplicated = cells.bamID[cells.bamID.duplicated()].unique()
    if len(duplicated) > 0:
        raise ValueError(
            f"Duplicated bamIDs in metadata: {', '.join(duplicated[:10])}"
        )
    return dict(zip(cells.barcode, cells.bamID))


def read_passes(read: pysam.AlignedSegment, min_mapq: int, exclude_flags: int) -> bool:
    """Apply MAPQ and flag filters equivalent to `samtools view -q <min_mapq> -F <exclude_flags>`."""
    return read.mapping_quality >= min_mapq and not (read.flag & exclude_flags)


def _cell_of(read: pysam.AlignedSegment, cell_tag: str, barcode_map: Dict[str, str]) -> Optional[str]:
    if not read.has_tag(cell_tag):
        return None
    return barcode_map.get(read.get_tag(cell_tag))


def readpos_path(outdir: Union[str, Path], cell: str, chrom: str) -> Path:
    return Path(outdir) / 'readpos' / cell / f'{chrom}.readpos.seq'


def extract_readpos_by_cell(
    bam_path: Union[str, Path],
    chrom: str,
    barcode_map: Dict[str, str],
    outdir: Union[str, Path],
    cell_tag: str = 'CB',
    min_mapq: int = READPOS_MIN_MAPQ,
    exclude_flags: int = READPOS_EXCLUDE_FLAGS,
    flush_every: int = FLUSH_EVERY,
) -> Dict[str, int]:
    """
    Stream one chromosome of a multi-cell BAM and write per-cell read start positions.

    Writes readpos/{cell}/{chrom}.readpos.seq (one 1-based POS per line, identical to
    `samtools view ... | perl -ane 'print $F[3]'`) for every cell in barcode_map,
    including empty files for cells without reads.

    Returns the number of reads written per cell.
    """
    cells = sorted(set(barcode_map.values()))

    # Create/truncate every output so all expected files exist
    for cell in cells:
        path = readpos_path(outdir, cell, chrom)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('')

    buffers: Dict[str, List[int]] = defaultdict(list)
    counts: Dict[str, int] = defaultdict(int)
    buffered = 0

    def flush() -> None:
        for cell, positions in buffers.items():
            if positions:
                with open(readpos_path(outdir, cell, chrom), 'a') as f:
                    f.write('\n'.join(map(str, positions)))
                    f.write('\n')
        buffers.clear()

    with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
        if chrom not in bam.references:
            logger.warning(f"{chrom} is not in the BAM header of {bam_path}; writing empty files")
            return {}
        for read in bam.fetch(chrom):
            if not read_passes(read, min_mapq, exclude_flags):
                continue
            cell = _cell_of(read, cell_tag, barcode_map)
            if cell is None:
                continue
            buffers[cell].append(read.reference_start + 1)
            counts[cell] += 1
            buffered += 1
            if buffered >= flush_every:
                flush()
                buffered = 0
    flush()

    logger.info(
        f"{chrom}: wrote {sum(counts.values())} read positions for "
        f"{len(counts)}/{len(cells)} cells"
    )
    return dict(counts)


def load_phased_sites(phased_vcf: Union[str, Path], chroms: Optional[Iterable[str]] = None) -> pd.DataFrame:
    """
    Read biallelic SNV sites from a phased hets VCF.

    Returns a DataFrame with columns CHROM, POS (1-based), REF, ALT.
    """
    skiprows = 0
    with gzip.open(phased_vcf, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                break
            skiprows += 1

    sites = pd.read_csv(
        phased_vcf,
        compression='gzip',
        sep='\t',
        skiprows=skiprows,
        usecols=['#CHROM', 'POS', 'REF', 'ALT'],
        dtype={'#CHROM': str, 'POS': int, 'REF': str, 'ALT': str},
        low_memory=False,
    ).rename(columns={'#CHROM': 'CHROM'})

    sites = sites[(sites.REF.str.len() == 1) & (sites.ALT.str.len() == 1)]
    if chroms is not None:
        sites = sites[sites.CHROM.isin(set(chroms))]
    return sites.drop_duplicates(['CHROM', 'POS']).sort_values(['CHROM', 'POS'])


def count_alleles_by_cell(
    bam_path: Union[str, Path],
    chrom: str,
    sites: pd.DataFrame,
    barcode_map: Dict[str, str],
    output_path: Union[str, Path],
    cell_tag: str = 'CB',
    min_mapq: int = ALLELE_MIN_MAPQ,
    exclude_flags: int = ALLELE_EXCLUDE_FLAGS,
    min_baseq: int = ALLELE_MIN_BASEQ,
) -> int:
    """
    Stream one chromosome of a multi-cell BAM and count REF/ALT reads per cell at het sites.

    `sites` holds the CHROM/POS/REF/ALT rows for this chromosome. Only cell x site pairs
    with at least one REF or ALT read are kept (sparse), and are written to output_path as
    a tab-separated table: cell, CHROM, POS, REF, ALT, AD_REF, AD_ALT.

    Returns the number of rows written.
    """
    sites = sites[sites.CHROM == chrom]
    # 0-based positions, sorted, for a quick overlap check per read
    site_pos0 = (sites.POS.values - 1).tolist()
    alleles = dict(zip(site_pos0, zip(sites.REF.str.upper(), sites.ALT.str.upper())))

    # (cell, pos0) -> [ref_count, alt_count]
    counts: Dict[Tuple[str, int], List[int]] = defaultdict(lambda: [0, 0])

    with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
        if chrom not in bam.references:
            logger.warning(f"{chrom} is not in the BAM header of {bam_path}; no allele counts")
            site_pos0 = []
        if site_pos0:
            for read in bam.fetch(chrom):
                if not read_passes(read, min_mapq, exclude_flags):
                    continue
                # Most reads overlap no het site; skip them before the per-base work
                i = bisect_left(site_pos0, read.reference_start)
                if i == len(site_pos0) or site_pos0[i] >= read.reference_end:
                    continue
                cell = _cell_of(read, cell_tag, barcode_map)
                if cell is None:
                    continue
                seq = read.query_sequence
                if seq is None:
                    continue
                quals = read.query_qualities
                for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                    if rpos not in alleles:
                        continue
                    if quals is not None and quals[qpos] < min_baseq:
                        continue
                    base = seq[qpos].upper()
                    ref, alt = alleles[rpos]
                    if base == ref:
                        counts[(cell, rpos)][0] += 1
                    elif base == alt:
                        counts[(cell, rpos)][1] += 1

    rows = [
        (cell, chrom, pos0 + 1, alleles[pos0][0], alleles[pos0][1], ref_n, alt_n)
        for (cell, pos0), (ref_n, alt_n) in counts.items()
    ]
    pd.DataFrame(
        rows,
        columns=['cell', 'CHROM', 'POS', 'REF', 'ALT', 'AD_REF', 'AD_ALT'],
    ).to_csv(output_path, sep='\t', index=False)

    logger.info(f"{chrom}: counted alleles at {len(site_pos0)} sites, {len(rows)} cell-site pairs")
    return len(rows)


def count_reads_by_tag(
    bam_path: Union[str, Path],
    cell_tag: str = 'CB',
    min_mapq: int = READPOS_MIN_MAPQ,
    exclude_flags: int = READPOS_EXCLUDE_FLAGS,
) -> Dict[str, int]:
    """Count passing reads per tag value over the whole BAM (one streaming pass)."""
    counts: Dict[str, int] = defaultdict(int)
    with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
        for read in bam.fetch(until_eof=True):
            if not read_passes(read, min_mapq, exclude_flags) or not read.has_tag(cell_tag):
                continue
            counts[read.get_tag(cell_tag)] += 1
    return dict(counts)


def read_barcode_list(path: Union[str, Path]) -> List[str]:
    """
    Read barcodes from a plain list (e.g. Cell Ranger barcodes.tsv[.gz], one per line)
    or from a CSV/TSV with a 'barcode' column (e.g. per_cell_summary_metrics.csv).
    """
    path = Path(path)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as f:
        first = f.readline().strip()
    sep = ',' if ',' in first else '\t'
    header = [h.strip() for h in first.split(sep)]
    if 'barcode' in header:
        table = pd.read_csv(path, sep=sep, dtype=str)
        barcodes = table['barcode'].dropna().tolist()
    else:
        table = pd.read_csv(path, sep=sep, header=None, dtype=str)
        barcodes = table[0].dropna().tolist()
    # Keep order, drop duplicates
    return list(dict.fromkeys(b.strip() for b in barcodes if b.strip()))
