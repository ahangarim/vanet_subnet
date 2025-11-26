"""
Task generation utilities for VANET subnet.

Generates random genomic tasks for validators with unique IDs and random regions.
"""

import os
import json
import random
import time
from typing import Dict, List, Optional, Any
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


class GenomicsTaskGenerator:
    """Generate random genomics tasks for validators."""

    def __init__(self, manifest_path: Optional[str] = None, seed: int = 42):
        self.manifest_path = manifest_path
        self.random_gen = random.Random(seed)
        self.manifest_data = {}

        self.chromosome = "chr20"
        self.contig_length = 64444167
        self.min_start = 10_000_000
        self.max_start = 55_000_000
        self.window_size = 5_000_000

        if manifest_path and os.path.exists(manifest_path):
            self.load_manifest()

    def load_manifest(self):
        """Load reference data manifest from JSON file."""
        try:
            with open(self.manifest_path, 'r') as f:
                self.manifest_data = json.load(f)

            # Load chromosome and contig length from manifest
            self.chromosome = self.manifest_data.get("chromosome", self.chromosome)
            self.contig_length = self.manifest_data.get("contig_length", self.contig_length)

            # Load region bounds from manifest if available
            bounds = self.manifest_data.get("region_bounds", {})
            self.min_start = bounds.get("min_start", self.min_start)
            self.max_start = bounds.get("max_start", self.max_start)
            self.window_size = bounds.get("window_size", self.window_size)

            logger.info(f"Loaded manifest: {self.manifest_data.get('description', 'reference data')}")
        except Exception as e:
            logger.error(f"Failed to load manifest: {e}")

    def generate_task(self, use_manifest: bool = True) -> Dict[str, Any]:
        """
        Generate a genomic task with unique ID and random region.

        Args:
            use_manifest: Ignored - always generates fresh random tasks

        Returns:
            Task dictionary with unique task_id and random region
        """
        # Generate random region within bounds
        start = self.random_gen.randint(self.min_start, self.max_start)
        end = start + self.window_size

        # Unique task ID using timestamp + random suffix
        timestamp = int(time.time())
        rand_suffix = self.random_gen.randint(1000, 9999)
        task_id = f"task_{self.chromosome}_{start}_{end}_{timestamp}_{rand_suffix}"

        return {
            "task_id": task_id,
            "task_type": "variant_calling",
            "data_locator": {},  # Filled by validator with actual BAM path
            "regions": [f"{self.chromosome}:{start}-{end}"],
            "ref_build": self.manifest_data.get("ref_build", "GRCh38"),
            "input_type": "BAM",
            "output_format": "VCF"
        }

    def get_reference_paths(self) -> Dict[str, Dict[str, str]]:
        """Get reference data paths from manifest (local or S3)."""
        data = self.manifest_data.get("data", {})
        return {
            "bam": data.get("bam", {}),
            "bam_index": data.get("bam_index", {}),
            "reference": data.get("reference", {}),
            "reference_index": data.get("reference_index", {}),
            "truth_vcf": data.get("truth_vcf", {}),
            "truth_bed": data.get("truth_bed", {})
        }

    def get_truth_for_task(self, task_id: str) -> Optional[Dict[str, Any]]:
        """Get truth data paths. Returns same paths for all tasks since truth is global."""
        data = self.manifest_data.get("data", {})
        if not data:
            return None

        return {
            "truth_vcf": data.get("truth_vcf", {}).get("local"),
            "truth_vcf_s3": data.get("truth_vcf", {}).get("s3"),
            "truth_bed": data.get("truth_bed", {}).get("local"),
            "truth_bed_s3": data.get("truth_bed", {}).get("s3"),
            "reference_fasta": data.get("reference", {}).get("local"),
            "reference_fasta_s3": data.get("reference", {}).get("s3")
        }
