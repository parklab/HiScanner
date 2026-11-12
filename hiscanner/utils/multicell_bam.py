"""
Support for multi-cell BAM input (e.g. 10x Genomics), where all cells share one
BAM and each read carries a cell barcode tag (CB by default).

The BAM is read one chromosome at a time and is never split into per-cell BAMs.
Outputs follow the per-cell file layout, so the downstream steps don't change.
"""
import gzip
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, Union

import pandas as pd
import pysam

from ..logger import logger

PathLike = Union[str, Path]

# Same as the per-cell `samtools view -q 30 -F 1284` (unmapped, secondary, duplicate)
READPOS_MIN_MAPQ = 30
READPOS_EXCLUDE_FLAGS = 4 | 256 | 1024

# Close to SCAN2's hc_raw.mmq60 input; also drops QC-fail and supplementary reads
ALLELE_MIN_MAPQ = 60
ALLELE_EXCLUDE_FLAGS = READPOS_EXCLUDE_FLAGS | 512 | 2048
ALLELE_MIN_BASEQ = 20

# Read positions kept in memory (all cells together) before writing to disk
FLUSH_EVERY = 1_000_000


def is_multi_cell(config: Dict) -> bool:
    return config.get('input_mode', 'per_cell_bam') == 'multi_cell_bam'


def find_bam_index(bam_path: PathLike) -> Optional[Path]:
    """Return the .bai/.csi index of a BAM, or None if there is none."""
    bam_path = Path(bam_path)
    for index in (Path(f'{bam_path}.bai'), Path(f'{bam_path}.csi'), bam_path.with_suffix('.bai')):
        if index.exists():
            return index
    return None


def load_barcode_map(metadata_path: PathLike) -> Dict[str, str]:
    """
    Map barcode -> bamID for the single cells in a multi-cell metadata file
    (columns: bamID, barcode, singlecell).
    """
    metadata = pd.read_csv(metadata_path, sep='\t', dtype=str)
    missing = {'bamID', 'barcode', 'singlecell'} - set(metadata.columns)
    if missing:
        raise ValueError(
            f"Metadata {metadata_path} is missing column(s) required for "
            f"multi_cell_bam mode: {', '.join(sorted(missing))}"
        )

    cells = metadata[metadata.singlecell == 'Y']
    for column in ('barcode', 'bamID'):
        dups = cells[column][cells[column].duplicated()].unique()
        if len(dups):
            raise ValueError(f"Duplicated {column}s in metadata: {', '.join(dups[:10])}")

    return dict(zip(cells.barcode, cells.bamID))


def read_passes(read: pysam.AlignedSegment, min_mapq: int, exclude_flags: int) -> bool:
    return read.mapping_quality >= min_mapq and not read.flag & exclude_flags


def _cell_of(read: pysam.AlignedSegment, cell_tag: str, barcode_map: Dict[str, str]) -> Optional[str]:
    if not read.has_tag(cell_tag):
        return None
    return barcode_map.get(read.get_tag(cell_tag))


def _check_chrom(bam: pysam.AlignmentFile, chrom: str, bam_path: PathLike) -> None:
    # samtools fails on an unknown region in per-cell mode; do the same here instead of
    # silently writing empty results (usually a 'chr' prefix mismatch with chrom_list)
    if chrom not in bam.references:
        raise ValueError(
            f"Chromosome {chrom} is not in the header of {bam_path} "
            f"(BAM has e.g. {', '.join(bam.references[:3])})"
        )


def readpos_path(outdir: PathLike, cell: str, chrom: str) -> Path:
    return Path(outdir) / 'readpos' / cell / f'{chrom}.readpos.seq'


def extract_readpos_by_cell(
    bam_path: PathLike,
    chrom: str,
    barcode_map: Dict[str, str],
    outdir: PathLike,
    cell_tag: str = 'CB',
    min_mapq: int = READPOS_MIN_MAPQ,
    exclude_flags: int = READPOS_EXCLUDE_FLAGS,
    flush_every: int = FLUSH_EVERY,
) -> Dict[str, int]:
    """
    Write readpos/{cell}/{chrom}.readpos.seq for every cell from one pass over a chromosome.

    Each file holds one 1-based read start per line, the same as the per-cell
    `samtools view | perl -ane 'print $F[3]'`. Cells without reads get an empty file.
    Returns the number of reads written per cell.
    """
    cells = sorted(set(barcode_map.values()))
    for cell in cells:
        path = readpos_path(outdir, cell, chrom)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('')

    buffers: Dict[str, List[int]] = defaultdict(list)
    counts: Dict[str, int] = defaultdict(int)
    buffered = 0

    def flush() -> None:
        for cell, positions in buffers.items():
            with open(readpos_path(outdir, cell, chrom), 'a') as f:
                f.write('\n'.join(map(str, positions)) + '\n')
        buffers.clear()

    with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
        _check_chrom(bam, chrom, bam_path)
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

    logger.info(f"{chrom}: wrote {sum(counts.values())} read positions for {len(counts)}/{len(cells)} cells")
    return dict(counts)


def _skip_vcf_header(vcf_path: PathLike) -> int:
    """Number of lines before the #CHROM line of a gzipped VCF."""
    with gzip.open(vcf_path, 'rt') as f:
        for n, line in enumerate(f):
            if line.startswith('#CHROM'):
                return n
    raise ValueError(f"No #CHROM header line in {vcf_path}")


def read_vcf_columns(vcf_path: PathLike, columns: Dict[str, type]) -> pd.DataFrame:
    """Read selected columns of a gzipped VCF; '#CHROM' is renamed to 'CHROM'."""
    return pd.read_csv(
        vcf_path,
        compression='gzip',
        sep='\t',
        skiprows=_skip_vcf_header(vcf_path),
        usecols=list(columns),
        dtype=columns,
        low_memory=False,
    ).rename(columns={'#CHROM': 'CHROM'})


def load_phased_sites(phased_vcf: PathLike, chroms: Optional[Iterable[str]] = None) -> pd.DataFrame:
    """Biallelic SNVs (CHROM, 1-based POS, REF, ALT) from the phased hets VCF."""
    sites = read_vcf_columns(phased_vcf, {'#CHROM': str, 'POS': int, 'REF': str, 'ALT': str})
    sites = sites[(sites.REF.str.len() == 1) & (sites.ALT.str.len() == 1)]
    if chroms is not None:
        sites = sites[sites.CHROM.isin(set(chroms))]
    return sites.drop_duplicates(['CHROM', 'POS']).sort_values(['CHROM', 'POS'])


def count_alleles_by_cell(
    bam_path: PathLike,
    chrom: str,
    sites: pd.DataFrame,
    barcode_map: Dict[str, str],
    output_path: PathLike,
    cell_tag: str = 'CB',
    min_mapq: int = ALLELE_MIN_MAPQ,
    exclude_flags: int = ALLELE_EXCLUDE_FLAGS,
    min_baseq: int = ALLELE_MIN_BASEQ,
) -> int:
    """
    Count REF/ALT reads per cell at the het sites of one chromosome.

    Writes only cell/site pairs with at least one REF or ALT read, as a TSV with
    columns cell, CHROM, POS, REF, ALT, AD_REF, AD_ALT. Overlapping mates of a
    read pair are counted once, as GATK does. Returns the number of rows written.
    """
    sites = sites[sites.CHROM == chrom]
    site_pos0 = (sites.POS.values - 1).tolist()
    alleles = dict(zip(site_pos0, zip(sites.REF.str.upper(), sites.ALT.str.upper())))

    # (cell, pos0) -> [ref, alt]
    counts: Dict[Tuple[str, int], List[int]] = defaultdict(lambda: [0, 0])
    # Sites already counted for the first mate of a pair whose second mate overlaps it.
    # read name -> (mate start, sites)
    pending: Dict[str, Tuple[int, set]] = {}

    with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
        _check_chrom(bam, chrom, bam_path)
        for read in bam.fetch(chrom) if site_pos0 else ():
            if not read_passes(read, min_mapq, exclude_flags):
                continue

            # Most reads cover no het site; skip them before looking at bases
            i = bisect_left(site_pos0, read.reference_start)
            if i == len(site_pos0) or site_pos0[i] >= read.reference_end:
                continue
            cell = _cell_of(read, cell_tag, barcode_map)
            if cell is None or read.query_sequence is None:
                continue

            already_counted = pending.pop(read.query_name, (None, ()))[1]
            mate_overlaps = (
                read.is_paired
                and not read.mate_is_unmapped
                and read.next_reference_id == read.reference_id
                and read.reference_start <= read.next_reference_start < read.reference_end
                and not already_counted
            )

            seq = read.query_sequence
            quals = read.query_qualities
            counted = set()
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if rpos not in alleles or rpos in already_counted:
                    continue
                if quals is not None and quals[qpos] < min_baseq:
                    continue
                ref, alt = alleles[rpos]
                base = seq[qpos].upper()
                if base == ref:
                    counts[(cell, rpos)][0] += 1
                elif base == alt:
                    counts[(cell, rpos)][1] += 1
                else:
                    continue
                counted.add(rpos)

            if mate_overlaps and counted:
                pending[read.query_name] = (read.next_reference_start, counted)

            # Mates that never showed up (filtered out) would otherwise pile up
            if len(pending) > 100_000:
                pending = {name: v for name, v in pending.items() if v[0] >= read.reference_start}

    rows = [
        (cell, chrom, pos0 + 1, *alleles[pos0], n_ref, n_alt)
        for (cell, pos0), (n_ref, n_alt) in counts.items()
    ]
    pd.DataFrame(
        rows, columns=['cell', 'CHROM', 'POS', 'REF', 'ALT', 'AD_REF', 'AD_ALT']
    ).to_csv(output_path, sep='\t', index=False)

    logger.info(f"{chrom}: {len(site_pos0)} het sites, {len(rows)} cell/site pairs with reads")
    return len(rows)


def count_reads_by_tag(
    bam_path: PathLike,
    cell_tag: str = 'CB',
    min_mapq: int = READPOS_MIN_MAPQ,
    exclude_flags: int = READPOS_EXCLUDE_FLAGS,
) -> Dict[str, int]:
    """Number of passing reads per barcode over the whole BAM."""
    counts: Dict[str, int] = defaultdict(int)
    with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
        for read in bam.fetch(until_eof=True):
            if read_passes(read, min_mapq, exclude_flags) and read.has_tag(cell_tag):
                counts[read.get_tag(cell_tag)] += 1
    return dict(counts)


def read_barcode_list(path: PathLike) -> List[str]:
    """
    Barcodes from a plain list (Cell Ranger barcodes.tsv[.gz], one per line) or from
    a CSV/TSV with a 'barcode' column (e.g. per_cell_summary_metrics.csv).
    """
    path = Path(path)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as f:
        first = f.readline().strip()

    sep = ',' if ',' in first else '\t'
    if 'barcode' in [h.strip() for h in first.split(sep)]:
        barcodes = pd.read_csv(path, sep=sep, dtype=str)['barcode']
    else:
        barcodes = pd.read_csv(path, sep=sep, header=None, dtype=str)[0]

    # keep the file order, drop blanks and duplicates
    return list(dict.fromkeys(b.strip() for b in barcodes.dropna() if b.strip()))
