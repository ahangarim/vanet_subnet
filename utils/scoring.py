"""
VCF scoring utilities using hap.py.

Validates variant calls against ground truth and computes accuracy metrics.
"""

import os
import csv
import subprocess
from typing import Dict, Optional
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


def subset_bed(source_bed: Path, target_bed: Path, region: str) -> bool:
    """Filter BED file to entries overlapping region (matching GenomeNet step 1)."""
    try:
        chrom, coords = region.split(":")
        start_str, end_str = coords.split("-")
        start = int(start_str.replace(",", ""))
        end = int(end_str.replace(",", ""))

        target_bed.parent.mkdir(parents=True, exist_ok=True)

        import gzip
        open_func = gzip.open if str(source_bed).endswith('.gz') else open

        entries_written = 0
        with open_func(source_bed, 'rt', encoding='utf-8') as src, \
             target_bed.open('w', encoding='utf-8') as dst:
            for line in src:
                if not line.strip() or line.startswith('#'):
                    continue

                parts = line.rstrip('\n').split('\t')
                if len(parts) < 3:
                    continue

                if parts[0] != chrom:
                    continue

                entry_start = int(parts[1])
                entry_end = int(parts[2])

                if entry_end <= start or entry_start >= end:
                    continue

                dst.write(line)
                entries_written += 1

        logger.info(f"Subset BED: {entries_written} entries in {region}")
        return True

    except Exception as e:
        logger.error(f"BED subset failed: {e}")
        return False


def slice_truth_vcf(source_vcf: Path, target_vcf: Path, region: str) -> bool:
    """Slice truth VCF to region using bcftools (matching GenomeNet step 1)."""
    try:
        source_vcf = Path(source_vcf).resolve()
        target_vcf = Path(target_vcf).resolve()
        target_vcf.parent.mkdir(parents=True, exist_ok=True)

        if not source_vcf.exists():
            logger.error(f"Source VCF not found: {source_vcf}")
            return False

        index_csi = source_vcf.with_suffix('.vcf.gz.csi')
        index_tbi = source_vcf.with_suffix('.vcf.gz.tbi')
        if not index_csi.exists() and not index_tbi.exists():
            logger.warning(f"VCF not indexed, extraction will be slow")

        source_dir = source_vcf.parent
        target_dir = target_vcf.parent

        if source_dir == target_dir:
            slice_cmd = f"""docker run --rm \
  -v {source_dir}:/data \
  biocontainers/bcftools:v1.9-1-deb_cv1 \
  bcftools view -r {region} /data/{source_vcf.name} -Oz -o /data/{target_vcf.name}"""
        else:
            slice_cmd = f"""docker run --rm \
  -v {source_dir}:/source \
  -v {target_dir}:/target \
  biocontainers/bcftools:v1.9-1-deb_cv1 \
  bcftools view -r {region} /source/{source_vcf.name} -Oz -o /target/{target_vcf.name}"""

        logger.info(f"Slicing truth VCF: {region}")
        result = subprocess.run(slice_cmd, shell=True, capture_output=True, text=True, timeout=300)

        if result.returncode != 0:
            logger.error(f"bcftools failed: {result.stderr}")
            return False

        if not target_vcf.exists():
            logger.error(f"Output not created: {target_vcf}")
            return False

        index_cmd = f"""docker run --rm \
  -v {target_dir}:/data \
  biocontainers/bcftools:v1.9-1-deb_cv1 \
  bcftools index /data/{target_vcf.name}"""

        result = subprocess.run(index_cmd, shell=True, capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            logger.warning(f"Index failed: {result.stderr}")

        logger.info(f"Created sliced VCF: {target_vcf.name}")
        return True

    except subprocess.TimeoutExpired:
        logger.error("VCF slicing timed out")
        return False
    except Exception as e:
        logger.error(f"VCF slicing failed: {e}")
        return False


class HappyScorer:
    """Score VCF outputs using hap.py (matching GenomeNet step 6)."""

    def __init__(self, use_docker: bool = True):
        self.use_docker = use_docker

    def score_vcf(self, truth_vcf: str, query_vcf: str,
                  reference_fasta: str = None, confident_bed: str = None) -> Dict[str, float]:
        """Score miner VCF against ground truth."""
        if self.use_docker and os.path.exists('/usr/bin/docker'):
            return self._score_with_happy(truth_vcf, query_vcf, reference_fasta, confident_bed)
        else:
            return self._score_simple(truth_vcf, query_vcf)

    def _score_with_happy(self, truth_vcf: str, query_vcf: str,
                         reference_fasta: str, confident_bed: str,
                         region: str = None) -> Dict[str, float]:
        """Run hap.py validation via Docker."""
        try:
            truth_vcf = Path(truth_vcf).resolve()
            query_vcf = Path(query_vcf).resolve()
            ref_path = Path(reference_fasta).resolve() if reference_fasta else None
            bed_path = Path(confident_bed).resolve() if confident_bed else None

            output_dir = query_vcf.parent
            output_prefix = output_dir / f"happy_{query_vcf.stem}"

            if not truth_vcf.exists():
                logger.error(f"Truth VCF not found: {truth_vcf}")
                return self._get_zero_scores()
            if not query_vcf.exists():
                logger.error(f"Query VCF not found: {query_vcf}")
                return self._get_zero_scores()
            if ref_path and not ref_path.exists():
                logger.error(f"Reference not found: {ref_path}")
                return self._get_zero_scores()

            # Validate BED file - if not found, run hap.py without confident regions filter
            use_bed = False
            subset_bed_path = None

            if bed_path and bed_path.exists():
                # Subset BED file to only include regions overlapping with task region
                # This improves performance and accuracy
                if region:
                    subset_bed_path = output_dir / f"confident_{region.replace(':', '_').replace('-', '_')}.bed"
                    logger.info(f"Creating subset BED for region {region}...")
                    if subset_bed(bed_path, subset_bed_path, region):
                        bed_path = subset_bed_path
                        use_bed = True
                        logger.info(f"Using subset confident regions BED: {bed_path}")
                    else:
                        logger.warning(f"Failed to subset BED, using full BED")
                        use_bed = True
                else:
                    use_bed = True
                    logger.info(f"Using full confident regions BED: {bed_path}")
            elif bed_path:
                logger.warning(f"Confident BED not found: {bed_path} - scoring without region filter")
            else:
                logger.info("No confident regions BED provided - scoring whole region")

            # Slice truth VCF to region for performance (like GenomeNet does)
            sliced_truth_vcf = None
            if region and truth_vcf.exists():
                sliced_truth_vcf = output_dir / f"truth_{region.replace(':', '_').replace('-', '_')}.vcf.gz"
                logger.info(f"Creating sliced truth VCF for region {region}...")
                if slice_truth_vcf(truth_vcf, sliced_truth_vcf, region):
                    truth_vcf = sliced_truth_vcf
                    logger.info(f"Using sliced truth VCF: {truth_vcf.name}")
                else:
                    logger.warning(f"Failed to slice truth VCF, using full chromosome VCF")

            # Build Docker command (matching GenomeNet exactly with vcfeval engine)
            # Mount BED file directory separately since merged truth VCF may be in different location
            happy_cmd = f"""docker run --rm \
  -e HGREF=/data/reference/{ref_path.name if ref_path else 'ref.fa'} \
  -v {truth_vcf.parent}:/data/truth \
  -v {query_vcf.parent}:/data/query \
  -v {ref_path.parent if ref_path else '/tmp'}:/data/reference \
  -v {output_dir}:/data/output"""

            # Add BED mount only if using BED and it's in a different directory than truth VCF
            if use_bed and bed_path.parent != truth_vcf.parent:
                happy_cmd += f" \\\n  -v {bed_path.parent}:/data/bed"
                bed_mount_path = f"/data/bed/{bed_path.name}"
            elif use_bed:
                bed_mount_path = f"/data/truth/{bed_path.name}"
            else:
                bed_mount_path = None

            # Check if SDF reference exists for vcfeval engine
            sdf_path = ref_path.parent / f"{ref_path.stem}.sdf" if ref_path else None
            use_vcfeval = sdf_path and sdf_path.exists() and sdf_path.is_dir()

            happy_cmd += f""" \\
  mgibio/hap.py:v0.3.12 \\
  /opt/hap.py/bin/hap.py \\
    /data/truth/{truth_vcf.name} \\
    /data/query/{query_vcf.name} \\
    -r /data/reference/{ref_path.name if ref_path else 'ref.fa'} \\
    -o /data/output/{output_prefix.name} \\
    --threads 4"""

            # Add vcfeval engine with SDF template (more accurate than default xcmp)
            if use_vcfeval:
                happy_cmd += f""" \\
    --engine vcfeval \\
    --engine-vcfeval-template /data/reference/{sdf_path.name}"""
                logger.info(f"Using vcfeval engine with SDF: {sdf_path.name}")
            else:
                logger.warning("SDF not found, using default hap.py engine (less accurate than vcfeval)")

            # Add confident regions filter only if BED file exists
            if use_bed and bed_mount_path:
                happy_cmd += f" \\\n    -f {bed_mount_path}"

            # Add region if specified (like test script does)
            if region:
                happy_cmd += f" \\\n    -l {region}"

            logger.info(f"Running hap.py validation on region: {region}")
            logger.debug(f"Command: {happy_cmd}")

            result = subprocess.run(happy_cmd, shell=True, capture_output=True, text=True, timeout=600)

            # Check if summary CSV was created (hap.py may return code 1 with warnings but still produce output)
            summary_csv = Path(f"{output_prefix}.summary.csv")

            if summary_csv.exists():
                logger.info(f"hap.py completed, parsing results from {summary_csv}")
                if result.returncode != 0:
                    logger.debug(f"hap.py returned code {result.returncode} (warnings only, output is valid)")

                # Parse results (same logic as test_step6_7_validator.py)
                happy_results = {}
                with open(summary_csv, 'r') as f:
                    reader = csv.DictReader(f)
                    # Debug: log CSV headers to catch format changes
                    if reader.fieldnames:
                        logger.debug(f"hap.py CSV columns: {reader.fieldnames}")
                    else:
                        logger.warning("hap.py CSV has no headers!")

                    rows_parsed = 0
                    for row in reader:
                        variant_type = row.get('Type', '')
                        logger.debug(f"Parsing row: Type={variant_type}, Filter={row.get('Filter', '')}")

                        if variant_type in ['INDEL', 'SNP']:
                            rows_parsed += 1

                            def safe_float(val):
                                try:
                                    return float(val) if val and val != 'nan' else 0.0
                                except:
                                    return 0.0

                            if variant_type == 'SNP':
                                happy_results['precision_snp'] = safe_float(row.get('METRIC.Precision', 0))
                                happy_results['recall_snp'] = safe_float(row.get('METRIC.Recall', 0))
                                happy_results['f1_snp'] = safe_float(row.get('METRIC.F1_Score', 0))
                                logger.debug(f"SNP metrics: P={happy_results['precision_snp']:.3f}, R={happy_results['recall_snp']:.3f}, F1={happy_results['f1_snp']:.3f}")
                            elif variant_type == 'INDEL':
                                happy_results['precision_indel'] = safe_float(row.get('METRIC.Precision', 0))
                                happy_results['recall_indel'] = safe_float(row.get('METRIC.Recall', 0))
                                happy_results['f1_indel'] = safe_float(row.get('METRIC.F1_Score', 0))
                                logger.debug(f"INDEL metrics: P={happy_results['precision_indel']:.3f}, R={happy_results['recall_indel']:.3f}, F1={happy_results['f1_indel']:.3f}")

                    if rows_parsed == 0:
                        logger.warning("hap.py CSV had no SNP/INDEL rows - check CSV format")

                # Fill in missing keys with defaults
                for key in ['f1_snp', 'f1_indel', 'precision_snp', 'recall_snp', 'precision_indel', 'recall_indel']:
                    if key not in happy_results:
                        happy_results[key] = 0.0

                # Calculate weighted F1
                happy_results['weighted_f1'] = 0.7 * happy_results['f1_snp'] + 0.3 * happy_results['f1_indel']

                logger.info(f"hap.py results: SNP F1={happy_results['f1_snp']:.3f}, INDEL F1={happy_results['f1_indel']:.3f}")
                return happy_results
            else:
                logger.error(f"hap.py failed - no summary CSV created. Return code: {result.returncode}")
                logger.error(f"stderr: {result.stderr[:500] if result.stderr else 'None'}")
                logger.error(f"stdout: {result.stdout[:500] if result.stdout else 'None'}")
                return self._get_zero_scores()

        except subprocess.TimeoutExpired:
            logger.error("hap.py execution timed out after 5 minutes")
            return self._get_zero_scores()
        except Exception as e:
            logger.error(f"Error running hap.py: {e}")
            import traceback
            logger.debug(traceback.format_exc())
            return self._get_zero_scores()

    def _score_simple(self, truth_vcf: str, query_vcf: str) -> Dict[str, float]:
        """
        Simple VCF scoring without hap.py (for testing).
        This is a simplified scoring that just compares variant positions.
        """
        try:
            truth_variants = self._load_vcf_positions(truth_vcf)
            query_variants = self._load_vcf_positions(query_vcf)

            if not truth_variants:
                return self._get_zero_scores()

            # Calculate basic metrics
            tp = len(truth_variants & query_variants)
            fp = len(query_variants - truth_variants)
            fn = len(truth_variants - query_variants)

            precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
            f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

            return {
                'f1_snp': f1,
                'f1_indel': f1 * 0.8,  # Simulate lower INDEL performance
                'precision_snp': precision,
                'recall_snp': recall,
                'precision_indel': precision * 0.8,
                'recall_indel': recall * 0.8,
                'weighted_f1': 0.7 * f1 + 0.3 * (f1 * 0.8)
            }
        except Exception as e:
            logger.error(f"Error in simple scoring: {e}")
            return self._get_zero_scores()

    def _load_vcf_positions(self, vcf_path: str) -> set:
        """Load variant positions from VCF file."""
        positions = set()

        if not os.path.exists(vcf_path):
            return positions

        with open(vcf_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    chrom, pos = parts[0], parts[1]
                    positions.add(f"{chrom}:{pos}")

        return positions

    def _get_zero_scores(self) -> Dict[str, float]:
        """Return zero scores structure."""
        return {
            'f1_snp': 0.0,
            'f1_indel': 0.0,
            'precision_snp': 0.0,
            'recall_snp': 0.0,
            'precision_indel': 0.0,
            'recall_indel': 0.0,
            'weighted_f1': 0.0
        }


class AdvancedScorer:
    """Advanced scoring with emphasis functions from GenomeNet."""

    @staticmethod
    def emphasis(metric: float, gamma: float = 3.0) -> float:
        """
        Apply nonlinear emphasis to push scores toward extremes.

        Args:
            metric: Raw metric (0-1)
            gamma: Emphasis power (higher = more extreme)

        Returns:
            Emphasized metric
        """
        return 1.0 - (1.0 - metric) ** gamma

    @staticmethod
    def compute_advanced_score(metrics: Dict[str, float]) -> float:
        """
        Compute advanced score with three components.

        Args:
            metrics: Dictionary with f1_snp, f1_indel, precision, recall, etc.

        Returns:
            Final score (0-100)
        """
        # Weighted F1
        weighted_f1 = 0.7 * metrics.get('f1_snp', 0) + 0.3 * metrics.get('f1_indel', 0)

        # Component 1: Core F1 with emphasis (50% weight)
        core_score = AdvancedScorer.emphasis(weighted_f1, gamma=4.0)

        # Component 2: Completeness (30% weight)
        avg_recall = (metrics.get('recall_snp', 0) + metrics.get('recall_indel', 0)) / 2
        completeness_score = AdvancedScorer.emphasis(avg_recall, gamma=3.0)

        # Component 3: Quality (20% weight)
        avg_precision = (metrics.get('precision_snp', 0) + metrics.get('precision_indel', 0)) / 2
        quality_score = AdvancedScorer.emphasis(avg_precision, gamma=2.0)

        # Final weighted score
        final_score = (
            0.5 * core_score +
            0.3 * completeness_score +
            0.2 * quality_score
        ) * 100

        return final_score
