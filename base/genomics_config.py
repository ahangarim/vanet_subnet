"""
Genomics Configuration Module

Configuration settings for the genomics variant-calling subnet.
"""

import os
from pathlib import Path

# Base paths
BASE_DIR = Path(__file__).parent.parent  # vanet/ directory (parent of base/)
DATA_DIR = BASE_DIR / "datasets"
MANIFEST_DIR = BASE_DIR / "base"  # Manifest in vanet/base/ folder

# Genomics task configuration
GENOMICS_CONFIG = {
    # Task generation
    "window_size": 5_000_000,  # 5Mb windows
    "chromosome": "chr20",
    "num_synthetic_mutations": 10,  # Synthetic SNPs for anti-cheating - reads actual ref bases from FASTA
    "reference_build": "GRCh38",

    # Scoring parameters
    "snp_weight": 0.7,
    "indel_weight": 0.3,
    "innovation_threshold": 0.05,  # 5% improvement for boost
    "boost_factor": 0.15,  # 15% boost multiplier
    "boost_half_life": 3600,  # 1 hour in seconds
    "ema_alpha": 0.1,  # Smoothing factor

    # Emphasis gammas
    "core_gamma": 4.0,
    "completeness_gamma": 3.0,
    "quality_gamma": 2.0,

    # Docker settings
    "use_docker": True,
    "use_bamsurgeon": True,  # Enable BamSurgeon for mutation injection
    "happy_docker_image": "mgibio/hap.py:v0.3.12",
    "gatk_docker_image": "broadinstitute/gatk:4.5.0.0",
    "deepvariant_docker_image": "google/deepvariant:1.5.0",
    "bamsurgeon_docker_image": "quay.io/biocontainers/bamsurgeon:1.4.1--pyhdfd78af_0",

    # Timeout settings
    "variant_calling_timeout": 3600,  # 1 hour - time for miner to complete GATK analysis
    "scoring_timeout": 600,  # 10 minutes - time for hap.py validation
    "task_interval": 14400,  # 4 hours - time between sending tasks to miners

    # S3 Configuration
    "s3_bucket": "vanetdata",  # Your S3 bucket name
    "s3_base_url": "https://vanetdata.s3.us-east-1.amazonaws.com",  # Public S3 URL base
    "reference_url": "https://vanetdata.s3.us-east-1.amazonaws.com/reference/chr20.fa",  # Reference genome on S3
    "bam_base_url": "https://vanetdata.s3.us-east-1.amazonaws.com/bams",   # S3 URL prefix for BAM files
    "truth_base_url": "https://vanetdata.s3.us-east-1.amazonaws.com/truth",  # S3 URL prefix for truth files

    # File paths (will be auto-downloaded from S3 if not present locally)
    "reference_fasta": "datasets/reference/chr20.fa",
    "reference_sdf": "datasets/reference/chr20.sdf",
    "giab_confident_bed": "datasets/truth/sample_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.chr20.bed",

    # Manifest file (contains reference data paths and region bounds)
    "manifests": [
        "s3_manifest.json"
    ],

    # Miner tool options
    "supported_callers": ["gatk", "deepvariant", "custom"],

    # Validation settings
    "min_variants_required": 10,  # Minimum variants in VCF
    "max_vcf_size_mb": 100,  # Maximum VCF file size

    # Local testing paths (when not using real genomic data)
    "use_synthetic_data": False,  # Set to True for testing without real BAMs
    "synthetic_bam_path": "/tmp/synthetic_test.bam",
    "synthetic_truth_vcf": "/tmp/synthetic_truth.vcf",
}

# Validator specific settings
VALIDATOR_CONFIG = {
    "query_batch_size": 10,  # Number of miners to query per round
    "scoring_interval": 100,  # Steps between weight updates
    "min_stake_requirement": 0.0,  # Minimum stake to participate
    "task_rotation_interval": 50,  # Steps before rotating to new tasks
}

# Miner specific settings
MINER_CONFIG = {
    "default_caller": "gatk",  # Default variant caller - use real GATK for production
    "cache_results": True,  # Cache results for repeated tasks
    "max_cache_size": 100,  # Maximum cached results
    "parallel_processing": False,  # Run multiple tasks in parallel
    "num_threads": 4,  # Threads per task
}


def get_data_path(relative_path: str) -> str:
    """Get absolute path for data files."""
    if GENOMICS_CONFIG["use_synthetic_data"]:
        # Return synthetic paths for testing
        if "bam" in relative_path.lower():
            return GENOMICS_CONFIG["synthetic_bam_path"]
        elif "vcf" in relative_path.lower():
            return GENOMICS_CONFIG["synthetic_truth_vcf"]

    # Check if path is already absolute
    if os.path.isabs(relative_path):
        return relative_path

    # Try relative to BASE_DIR (vanet directory)
    local_path = BASE_DIR / relative_path
    if local_path.exists():
        return str(local_path.resolve())

    # Return as-is if not found (might be S3 URI or other remote path)
    return relative_path


def get_manifest_path(manifest_name: str = None) -> str:
    """Get path to manifest file."""
    if manifest_name:
        return str(MANIFEST_DIR / manifest_name)

    # Return first available manifest
    for manifest in GENOMICS_CONFIG["manifests"]:
        path = MANIFEST_DIR / manifest
        if path.exists():
            return str(path)

    # Return empty if no manifest found
    return None


def is_docker_available() -> bool:
    """Check if Docker is available on the system."""
    try:
        import subprocess
        result = subprocess.run(["docker", "--version"], capture_output=True, timeout=5)
        return result.returncode == 0
    except:
        return False