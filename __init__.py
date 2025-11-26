"""
VANET - Genomics Variant Calling Subnet

A Bittensor subnet for distributed genomic variant calling
using GATK HaplotypeCaller with hap.py validation.
"""

__version__ = "0.1.0"

# Import base components for library usage
from .base import TaskSynapse, ComputeSynapse, GenomicsTaskSynapse
from .utils import GenomicsTaskGenerator, MutationInjector, HappyScorer

__all__ = [
    "TaskSynapse",
    "ComputeSynapse",
    "GenomicsTaskSynapse",
    "GenomicsTaskGenerator",
    "MutationInjector",
    "HappyScorer",
]
