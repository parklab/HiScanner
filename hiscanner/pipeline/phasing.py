import os
import subprocess
import pandas as pd
from pathlib import Path
from typing import Dict, Any, List, Optional
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from ..logger import logger
from .snp_calling import SNPCallingError
import numpy as np
import gzip
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

        if config.get('input_mode', 'per_cell_bam') == 'multi_cell_bam':
            # Allele counts come from the multi-cell BAM; only phased hets are needed
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

        if config.get('input_mode', 'per_cell_bam') == 'multi_cell_bam':
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

    REF/ALT read counts at phased het sites are taken directly from the BAM, one
    streaming pass per chromosome, and written as {cell}.raw.txt in the same format
    as the bcftools path. Phasing and BAF are then computed by apply_phasing.
    """
    from ..utils.multicell_bam import (
        load_barcode_map, load_phased_sites, count_alleles_by_cell
    )

    output_dir = paths['phased_hets_dir']
    rerun = config.get('rerun', False)
    barcode_map = load_barcode_map(paths['metadata_path'])
    single_cells = sorted(set(barcode_map.values()))
    if not single_cells:
        raise PhasingError("No single cells found in metadata")

    todo = [c for c in single_cells
            if rerun or not (output_dir / f'{c}.hetsnp.txt').exists()]
    if not todo:
        logger.info("Skipping existing results for all cells")
        return
    logger.info(f"Processing {len(todo)} cells from {paths['multi_cell_bam']}")

    # Only count for cells that still need results
    barcode_map = {bc: cell for bc, cell in barcode_map.items() if cell in set(todo)}

    sites = load_phased_sites(paths['phased_vcf'], config.get('chrom_list'))
    chroms = list(dict.fromkeys(sites.CHROM))
    if not chroms:
        raise PhasingError(
            f"No biallelic SNVs in {paths['phased_vcf']} on chromosomes {config.get('chrom_list')}"
        )
    logger.info(f"Loaded {len(sites)} phased het SNVs on {len(chroms)} chromosomes")

    counts_dir = output_dir / 'multicell_counts'
    counts_dir.mkdir(parents=True, exist_ok=True)

    # One process per chromosome; each streams its region of the BAM independently
    with ProcessPoolExecutor(max_workers=config.get('threads', 1)) as executor:
        futures = {
            chrom: executor.submit(
                count_alleles_by_cell,
                paths['multi_cell_bam'],
                chrom,
                sites[sites.CHROM == chrom],
                barcode_map,
                counts_dir / f'{chrom}.counts.tsv',
                config.get('cell_tag', 'CB')
            )
            for chrom in chroms
        }
        for chrom, future in futures.items():
            future.result()

    counts = pd.concat(
        [pd.read_csv(counts_dir / f'{chrom}.counts.tsv', sep='\t',
                     dtype={'cell': str, 'CHROM': str})
         for chrom in chroms],
        ignore_index=True
    )
    by_cell = dict(tuple(counts.groupby('cell')))
    empty = counts.iloc[0:0]

    phase_data = load_phase_data(paths['phased_vcf'])
    keep_raw_files = config.get('keep_raw_files', False)

    def finish_cell(cell: str) -> None:
        raw_output = output_dir / f'{cell}.raw.txt'
        cell_counts = by_cell.get(cell, empty)
        if cell_counts.empty:
            logger.warning(f"No reads at phased het sites for {cell}")
            pd.DataFrame(columns=['CHROM', 'POS', 'A', 'B', 'TOTAL', 'BAF']).to_csv(
                output_dir / f'{cell}.hetsnp.txt', sep='\t', index=False
            )
            return
        cell_counts[['CHROM', 'POS', 'REF', 'ALT', 'AD_REF', 'AD_ALT']].to_csv(
            raw_output, sep='\t', header=False, index=False
        )
        apply_phasing(cell, raw_output, phase_data, output_dir / f'{cell}.hetsnp.txt')
        if not keep_raw_files:
            raw_output.unlink(missing_ok=True)

    with ThreadPoolExecutor(max_workers=config.get('threads', 1)) as executor:
        for future in [executor.submit(finish_cell, cell) for cell in todo]:
            future.result()

    if not keep_raw_files:
        for chrom in chroms:
            (counts_dir / f'{chrom}.counts.tsv').unlink(missing_ok=True)
        try:
            counts_dir.rmdir()
        except OSError:
            pass

def load_phase_data(phased_vcf: Path) -> pd.DataFrame:
    """
    Read CHROM, POS and phased genotype (phasedgt column) from the phased hets VCF.
    """
    # Count header lines in phased VCF
    skiprows = 0
    with gzip.open(phased_vcf, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                break
            skiprows += 1
    logger.info(f"Skipping {skiprows} rows in phased VCF (header)")

    # Read phasing information using gzip
    phase_data = pd.read_csv(
        phased_vcf,
        compression='gzip',
        sep='\t',
        skiprows=skiprows,
        usecols=['#CHROM', 'POS', 'phasedgt'],
        low_memory=False,
        dtype={
            '#CHROM': str,
            'POS': int,
            'phasedgt': str
        }
    )
    phase_data.columns = ['CHROM', 'POS', 'GT']
    return phase_data

def apply_phasing(
    cell: str,
    raw_output: Path,
    phase_data: pd.DataFrame,
    cell_output: Path
) -> None:
    """
    Merge per-cell REF/ALT counts (raw_output) with phased genotypes and write
    haplotype-specific counts and BAF to cell_output.
    """
    raw_data = pd.read_csv(
        raw_output,
        sep='\t',
        header=None,
        names=['CHROM', 'POS', 'REF', 'ALT', 'AD_REF', 'AD_ALT'],
        low_memory=False
    )

    # Replace '.' with NaN and drop those rows
    raw_data = raw_data.replace('.', np.nan)
    raw_data = raw_data.dropna()

    # Convert columns to appropriate types after cleaning
    raw_data = raw_data.astype({
        'CHROM': str,
        'POS': int,
        'REF': str,
        'ALT': str,
        'AD_REF': int,
        'AD_ALT': int
    })

    logger.info(f"Merging raw and phased data for {cell}")
    # Merge and process
    merged = pd.merge(
        raw_data,
        phase_data,
        on=['CHROM', 'POS'],
        how='inner'
    )

    logger.info(f"Applying phasing for {cell}")
    # Apply phasing
    merged['A'] = merged.apply(
        lambda x: x.AD_ALT if x.GT == '1|0' else x.AD_REF,
        axis=1
    )
    merged['B'] = merged.apply(
        lambda x: x.AD_REF if x.GT == '1|0' else x.AD_ALT,
        axis=1
    )

    logger.info(f"Calculating BAF for {cell}")
    # Calculate BAF
    merged['TOTAL'] = merged['A'] + merged['B']
    merged = merged[merged.TOTAL > 0]
    merged['BAF'] = merged.apply(
        lambda x: min(x.A, x.B) / x.TOTAL,
        axis=1
    )

    # Save results
    result = merged[['CHROM', 'POS', 'A', 'B', 'TOTAL', 'BAF']]
    result.to_csv(cell_output, sep='\t', index=False)
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
