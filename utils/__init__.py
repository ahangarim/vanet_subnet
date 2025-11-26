"""
VANET Genomics Utilities

Modular utilities for genomics task generation, mutation injection,
scoring, and weight tracking in the VANET subnet.
"""

# Task generation
from .task_generation import GenomicsTaskGenerator

# Mutation injection
from .mutation_injection import MutationInjector

# Scoring
from .scoring import HappyScorer, AdvancedScorer, subset_bed, slice_truth_vcf

# Weight tracking
from .weight_tracking import ScoreTracker

# File utilities
from .file_utils import download_file, ensure_s3_file

# GATK utilities
from .gatk_utils import run_gatk_docker, ensure_reference_indexes

__all__ = [
    # Task generation
    'GenomicsTaskGenerator',

    # Mutation injection
    'MutationInjector',

    # Scoring
    'HappyScorer',
    'AdvancedScorer',
    'subset_bed',
    'slice_truth_vcf',

    # Weight tracking
    'ScoreTracker',

    # File utilities
    'download_file',
    'ensure_s3_file',

    # GATK utilities
    'run_gatk_docker',
    'ensure_reference_indexes',
]
