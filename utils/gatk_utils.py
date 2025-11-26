"""
GATK utilities for variant calling.

Handles GATK HaplotypeCaller execution and reference genome indexing.
"""

import subprocess
from typing import Optional
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


def ensure_reference_indexes(reference_path: Path) -> bool:
    """
    Ensure reference genome has required indexes (dict and fai).
    Uses the same Docker commands as test_step5_miner.py.

    Args:
        reference_path: Path to reference FASTA

    Returns:
        True if all indexes exist or were created, False otherwise
    """
    reference_path = Path(reference_path).resolve()

    # Check/create reference dict (required by GATK)
    # Handle both .fa and .fasta extensions
    if reference_path.suffix == '.fa':
        ref_dict = reference_path.with_suffix('.dict')
    elif reference_path.suffix == '.fasta':
        ref_dict = reference_path.with_suffix('.dict')
    else:
        # Handle cases like .fa.gz or other compound extensions
        ref_dict = Path(str(reference_path).replace('.fa', '.dict').replace('.fasta', '.dict'))

    if not ref_dict.exists():
        logger.info("Creating reference dictionary (required by GATK)...")
        dict_cmd = f"""docker run --rm \
  -v {reference_path.parent}:/data/reference \
  broadinstitute/gatk:4.5.0.0 \
  gatk CreateSequenceDictionary \
    -R /data/reference/{reference_path.name} \
    -O /data/reference/{ref_dict.name}"""

        result = subprocess.run(dict_cmd, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            logger.info(f"Reference dict created: {ref_dict}")
        else:
            logger.error(f"Failed to create dict: {result.stderr[:200] if result.stderr else 'Unknown error'}")
            return False
    else:
        logger.debug(f"Reference dict exists: {ref_dict}")

    # Check/create reference index (required by GATK)
    ref_fai = Path(str(reference_path) + ".fai")
    if not ref_fai.exists():
        logger.info("Creating reference index (required by GATK)...")
        fai_cmd = f"""docker run --rm \
  -v {reference_path.parent}:/data/reference \
  biocontainers/samtools:v1.9-4-deb_cv1 \
  samtools faidx /data/reference/{reference_path.name}"""

        result = subprocess.run(fai_cmd, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            logger.info(f"Reference index created: {ref_fai}")
        else:
            logger.error(f"Failed to create index: {result.stderr[:200] if result.stderr else 'Unknown error'}")
            return False
    else:
        logger.debug(f"Reference index exists: {ref_fai}")

    return True


def run_gatk_docker(bam_path: Path, reference_path: Path, output_vcf: Path,
                    interval: Optional[str] = None, threads: int = 4) -> bool:
    """
    Run GATK HaplotypeCaller using Docker.
    Uses the exact same Docker command format as test_step5_miner.py.

    Args:
        bam_path: Path to BAM file
        reference_path: Path to reference FASTA
        output_vcf: Path for output VCF
        interval: Genomic interval (e.g., "chr20:10000000-15000000")
        threads: Number of threads for native pair-hmm

    Returns:
        True if successful, False otherwise
    """
    # Ensure paths are absolute
    bam_path = Path(bam_path).resolve()
    reference_path = Path(reference_path).resolve()
    output_vcf = Path(output_vcf).resolve()

    # Create output directory
    output_vcf.parent.mkdir(parents=True, exist_ok=True)

    # Ensure reference indexes exist
    if not ensure_reference_indexes(reference_path):
        logger.error("Failed to create reference indexes")
        return False

    # Build GATK Docker command (same format as test_step5_miner.py)
    gatk_cmd = f"""docker run --rm \
  -v {bam_path.parent}:/data/bams \
  -v {reference_path.parent}:/data/reference \
  -v {output_vcf.parent}:/data/output \
  broadinstitute/gatk:4.5.0.0 \
  gatk HaplotypeCaller \
    -R /data/reference/{reference_path.name} \
    -I /data/bams/{bam_path.name} \
    -O /data/output/{output_vcf.name} \
    --native-pair-hmm-threads {threads}"""

    if interval:
        gatk_cmd += f" \\\n    -L {interval}"

    try:
        logger.info(f"Running GATK HaplotypeCaller on region: {interval or 'full genome'}...")
        logger.debug(f"Command: {gatk_cmd}")

        result = subprocess.run(gatk_cmd, shell=True, capture_output=True, text=True, timeout=3600)

        if result.returncode == 0:
            # Verify output VCF was created
            if output_vcf.exists() and output_vcf.stat().st_size > 0:
                # Count variants
                with open(output_vcf, 'r') as f:
                    variants_called = sum(1 for line in f if not line.startswith('#'))
                logger.info(f"GATK completed successfully. Variants called: {variants_called}")
                return True
            else:
                logger.error(f"GATK completed but VCF not created: {output_vcf}")
                return False
        else:
            logger.error(f"GATK failed with return code {result.returncode}")
            logger.error(f"stderr: {result.stderr[:500] if result.stderr else 'None'}")
            return False

    except subprocess.TimeoutExpired:
        logger.error("GATK timed out after 1 hour")
        return False
    except Exception as e:
        logger.error(f"Error running GATK: {e}")
        return False
