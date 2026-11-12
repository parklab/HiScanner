import os
import subprocess
import pandas as pd
from contextlib import suppress
from pathlib import Path
from typing import Dict, Any, List, Optional
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from ..logger import logger
from ..utils.multicell_bam import is_multi_cell, read_vcf_columns
from .snp_calling import SNPCallingError
import numpy as np

class PhasingError(Exception):
    pass

def validate_inputs(config: Dict[str, Any]) -> Dict[str, Path]:
    """
    Validate input files and paths required for phasing.

    Parameters
    ----------
    config : Dict[str, Any]
        Pipeline configuration

    Returns
    -------
    Dict[str, Path]
        Dictionary of validated paths
    """
    try:
        # Get paths from config
        output_dir = Path(config['outdir'])
        phased_hets_dir = output_dir / 'phased_hets'

        # Create output directory
        phased_hets_dir.mkdir(parents=True, exist_ok=True)

        # Validate metadata file
        metadata_path = Path(config['metadata_path'])
        if not metadata_path.exists():
            raise PhasingError(f"Metadata file not found: {metadata_path}")

        if is_multi_cell(config):
            # allele counts come straight from the shared BAM, so only the phased hets are needed
            from .snp_calling import get_phased_vcf
            return {
                'output_dir': output_dir,
                'phased_hets_dir': phased_hets_dir,
                'metadata_path': metadata_path,
                'multi_cell_bam': Path(config['multi_cell_bam']),
                'phased_vcf': get_phased_vcf(config)
            }

        # Validate SCAN2 outputs through snp_calling module
        from .snp_calling import prepare_scan2_results
        scan2_files = prepare_scan2_results(config)

        return {
            'output_dir': output_dir,
            'phased_hets_dir': phased_hets_dir,
            'metadata_path': metadata_path,
            'raw_vcf': scan2_files['raw_vcf'],
            'phased_vcf': scan2_files['phased_vcf']
        }
    except Exception as e:
        raise PhasingError(f"Input validation failed: {e}")




def run_phasing(config: Dict[str, Any]) -> None:
    """
    Run the phasing step of the pipeline.

    Parameters
    ----------
    config : Dict[str, Any]
        Pipeline configuration
    """
    logger.info("Starting phasing step")

    try:
        # Validate inputs
        paths = validate_inputs(config)

        if is_multi_cell(config):
            run_phasing_multicell(config, paths)
            logger.info("Phasing step completed successfully")
            return

        # Read metadata
        metadata = pd.read_csv(paths['metadata_path'], sep='\t')
        single_cells = metadata[metadata.singlecell == 'Y']['bamID'].tolist()

        if not single_cells:
            raise PhasingError("No single cells found in metadata")

        logger.info(f"Processing {len(single_cells)} cells")

        # Process cells in parallel
        with ThreadPoolExecutor(max_workers=config.get('threads', 1)) as executor:
            futures = []
            for cell in single_cells:
                futures.append(
                    executor.submit(
                        process_cell_phasing,
                        cell,
                        paths['raw_vcf'],
                        paths['phased_vcf'],
                        paths['phased_hets_dir'],
                        config.get('rerun', False),
                        config  # Pass the full config dictionary
                    )
                )

            # Wait for all processes and check for errors
            for future in futures:
                future.result()

        logger.info("Phasing step completed successfully")

    except Exception as e:
        raise PhasingError(f"Error in phasing step: {e}")

def run_phasing_multicell(config: Dict[str, Any], paths: Dict[str, Path]) -> None:
    """
    Phasing for a shared multi-cell BAM (e.g. 10x, cells identified by the CB tag).

    REF/ALT counts at the phased het sites are taken from the BAM (one pass per
    chromosome) and written as {cell}.raw.txt, the same format the bcftools path
    produces, so apply_phasing can be reused as is.
    """
    from ..utils.multicell_bam import load_barcode_map, load_phased_sites, count_alleles_by_cell

    output_dir = paths['phased_hets_dir']
    threads = config.get('threads', 1)
    barcode_map = load_barcode_map(paths['metadata_path'])
    if not barcode_map:
        raise PhasingError("No single cells found in metadata")

    todo = sorted(
        cell for cell in set(barcode_map.values())
        if config.get('rerun', False) or not (output_dir / f'{cell}.hetsnp.txt').exists()
    )
    if not todo:
        logger.info("Skipping existing results for all cells")
        return
    logger.info(f"Processing {len(todo)} cells from {paths['multi_cell_bam']}")

    # no need to count reads for cells that are already done
    todo_set = set(todo)
    barcode_map = {bc: cell for bc, cell in barcode_map.items() if cell in todo_set}

    sites = load_phased_sites(paths['phased_vcf'], config.get('chrom_list'))
    chroms = list(dict.fromkeys(sites.CHROM))
    if not chroms:
        raise PhasingError(
            f"No biallelic SNVs in {paths['phased_vcf']} on chromosomes {config.get('chrom_list')}"
        )
    logger.info(f"Loaded {len(sites)} phased het SNVs on {len(chroms)} chromosomes")

    counts_dir = output_dir / 'multicell_counts'
    counts_dir.mkdir(parents=True, exist_ok=True)
    count_files = [counts_dir / f'{chrom}.counts.tsv' for chrom in chroms]

    # pysam + per-base work is CPU bound, so use processes (one chromosome each)
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(
                count_alleles_by_cell,
                paths['multi_cell_bam'],
                chrom,
                sites[sites.CHROM == chrom],
                barcode_map,
                count_file,
                config.get('cell_tag', 'CB'),
            )
            for chrom, count_file in zip(chroms, count_files)
        ]
        for future in futures:
            future.result()

    counts = pd.concat(
        [pd.read_csv(f, sep='\t', dtype={'cell': str, 'CHROM': str}) for f in count_files],
        ignore_index=True,
    )
    by_cell = dict(tuple(counts.groupby('cell')))
    phase_data = load_phase_data(paths['phased_vcf'])
    keep_raw_files = config.get('keep_raw_files', False)

    def finish_cell(cell: str) -> None:
        raw_output = output_dir / f'{cell}.raw.txt'
        cell_output = output_dir / f'{cell}.hetsnp.txt'
        cell_counts = by_cell.get(cell)
        if cell_counts is None:
            logger.warning(f"No reads at phased het sites for {cell}")
            pd.DataFrame(columns=['CHROM', 'POS', 'A', 'B', 'TOTAL', 'BAF']).to_csv(
                cell_output, sep='\t', index=False
            )
            return
        cell_counts[['CHROM', 'POS', 'REF', 'ALT', 'AD_REF', 'AD_ALT']].to_csv(
            raw_output, sep='\t', header=False, index=False
        )
        apply_phasing(cell, raw_output, phase_data, cell_output)
        if not keep_raw_files:
            raw_output.unlink(missing_ok=True)

    with ThreadPoolExecutor(max_workers=threads) as executor:
        for future in [executor.submit(finish_cell, cell) for cell in todo]:
            future.result()

    if not keep_raw_files:
        for count_file in count_files:
            count_file.unlink(missing_ok=True)
        with suppress(OSError):
            counts_dir.rmdir()

def load_phase_data(phased_vcf: Path) -> pd.DataFrame:
    """CHROM, POS and phased genotype (phasedgt column) from the phased hets VCF."""
    phase_data = read_vcf_columns(phased_vcf, {'#CHROM': str, 'POS': int, 'phasedgt': str})
    return phase_data.rename(columns={'phasedgt': 'GT'})

def apply_phasing(
    cell: str,
    raw_output: Path,
    phase_data: pd.DataFrame,
    cell_output: Path
) -> None:
    """
    Combine a cell's REF/ALT counts (raw_output) with the phased genotypes and
    write haplotype counts A/B and BAF to cell_output.
    """
    raw_data = pd.read_csv(
        raw_output,
        sep='\t',
        header=None,
        names=['CHROM', 'POS', 'REF', 'ALT', 'AD_REF', 'AD_ALT'],
        low_memory=False
    )

    # bcftools writes '.' for missing AD
    raw_data = raw_data.replace('.', np.nan).dropna()
    raw_data = raw_data.astype({
        'CHROM': str,
        'POS': int,
        'REF': str,
        'ALT': str,
        'AD_REF': int,
        'AD_ALT': int
    })

    logger.info(f"Applying phasing for {cell}")
    merged = raw_data.merge(phase_data, on=['CHROM', 'POS'], how='inner')

    # haplotype A carries ALT when the phased genotype is 1|0, REF otherwise
    alt_on_a = merged.GT == '1|0'
    merged['A'] = np.where(alt_on_a, merged.AD_ALT, merged.AD_REF)
    merged['B'] = np.where(alt_on_a, merged.AD_REF, merged.AD_ALT)
    merged['TOTAL'] = merged.A + merged.B
    merged = merged[merged.TOTAL > 0].copy()
    merged['BAF'] = np.minimum(merged.A, merged.B) / merged.TOTAL

    merged[['CHROM', 'POS', 'A', 'B', 'TOTAL', 'BAF']].to_csv(cell_output, sep='\t', index=False)
    logger.info(f"Results saved for {cell}")

def process_cell_phasing(
    cell: str,
    raw_vcf: Path,
    phased_vcf: Path,
    output_dir: Path,
    rerun: bool = False,
    config: Optional[Dict[str, Any]] = None
) -> None:
    """
    Process phasing for a single cell.

    Parameters
    ----------
    cell : str
        Cell identifier
    raw_vcf : Path
        Path to raw VCF file
    phased_vcf : Path
        Path to phased VCF file
    output_dir : Path
        Output directory
    rerun : bool, optional
        Whether to rerun existing results
    config : Dict[str, Any], optional
        Configuration dictionary
    """
    cell_output = output_dir / f'{cell}.hetsnp.txt'
    raw_output = output_dir / f'{cell}.raw.txt'

    if cell_output.exists() and not rerun:
        logger.info(f"Skipping existing results for {cell}")
        return

    try:
        if not raw_output.exists() or rerun:
            # Extract heterozygous SNPs using bcftools
            logger.info(f"Running bcftools query for {cell}")
            cmd = [
                'bcftools', 'query',
                '-s', cell,
                '-f', '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD{0}\\t%AD{1}]\\n',
                str(raw_vcf)
            ]

            with open(raw_output, 'w') as f:
                subprocess.run(cmd, stdout=f, check=True)
        else:
            logger.info(f"Using existing raw file for {cell}")

        logger.info(f"Processing phasing information for {cell}")
        phase_data = load_phase_data(phased_vcf)
        apply_phasing(cell, raw_output, phase_data, cell_output)

        # Clean up intermediate file if needed
        # Make this operation more defensive
        keep_raw_files = config.get('keep_raw_files', False) if config else False
        if not keep_raw_files and raw_output.exists():
            try:
                raw_output.unlink()
            except Exception as e:
                logger.warning(f"Could not remove raw output file for {cell}: {e}")

    except Exception as e:
        logger.error(f"Detailed error processing cell {cell}: {str(e)}")
        raise PhasingError(f"Error processing cell {cell}: {str(e)}")


if __name__ == "__main__":
    # Example usage
    example_config = {
        'outdir': './output',
        'metadata_path': './metadata.txt',
        'scan2_output': './scan2_out',
        'threads': 4
    }

    try:
        run_phasing(example_config)
    except PhasingError as e:
        print(f"Error: {e}")
