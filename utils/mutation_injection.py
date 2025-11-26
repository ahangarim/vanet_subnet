"""
Mutation injection utilities using BamSurgeon.

Generates synthetic mutations and injects them into BAM files for anti-cheating.
"""

import os
import random
import subprocess
import shutil
import time
import platform
from typing import Dict, List, Optional
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


class MutationInjector:
    """Inject synthetic mutations using BamSurgeon (matching GenomeNet step 2)."""

    def __init__(self, seed: int = None):
        self.random_gen = random.Random(seed)
        self.current_seed = seed

    def reseed(self, seed: int = None):
        """Reseed RNG for fresh mutations each round (anti-cheating)."""
        if seed is None:
            seed = int(time.time() * 1000000) % (2**32)
        self.current_seed = seed
        self.random_gen = random.Random(seed)
        return seed

    def inject_mutations(self, task_id: str, region: str, num_mutations: int = 200,
                         reference_path: Optional[Path] = None) -> List[Dict]:
        """Generate synthetic SNP mutations for region. Reads actual ref bases from FASTA if provided."""
        chrom, coords = region.split(':')
        start, end = map(int, coords.split('-'))

        mutations = []
        nucleotides = ['A', 'C', 'G', 'T']

        positions = set()
        attempts = 0
        max_attempts = num_mutations * 10
        while len(positions) < num_mutations and attempts < max_attempts:
            pos = self.random_gen.randint(start, end)
            positions.add(pos)
            attempts += 1

        # If reference path is provided, read actual ref bases from FASTA
        # This matches GenomeNet's approach in make_window_mutations.py
        # IMPORTANT: pysam is REQUIRED for proper mutation injection!
        # Without pysam, mutations will use random ref bases which may fail BamSurgeon validation
        fasta_handle = None
        if reference_path:
            try:
                import pysam
                reference_path = Path(reference_path).resolve()
                if reference_path.exists():
                    fasta_handle = pysam.FastaFile(str(reference_path))
                    logger.info(f"Reading ref bases from FASTA: {reference_path}")
                else:
                    logger.warning(f"Reference FASTA not found: {reference_path}, using random ref bases")
            except ImportError:
                logger.error("pysam not installed! Install with: pip install pysam")
                logger.error("Without pysam, BamSurgeon mutations may fail validation")
                print("   ⚠️  WARNING: pysam not installed - mutations may fail!", flush=True)
                print("   Install pysam: pip install pysam", flush=True)
            except Exception as e:
                logger.warning(f"Failed to open FASTA: {e}, using random ref bases")

        for pos in sorted(positions):
            if fasta_handle:
                # Read actual ref base from FASTA (like GenomeNet does)
                try:
                    ref = fasta_handle.fetch(chrom, pos - 1, pos).upper()
                    if ref not in nucleotides:
                        # Skip positions with non-standard bases (N, etc)
                        continue
                except Exception as e:
                    logger.debug(f"Could not fetch ref at {chrom}:{pos}: {e}")
                    continue
            else:
                # Fallback: random reference (may fail BamSurgeon validation)
                ref = self.random_gen.choice(nucleotides)

            # Pick random alternate allele different from ref
            alt_choices = [n for n in nucleotides if n != ref]
            alt = self.random_gen.choice(alt_choices)

            mutations.append({
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'task_id': task_id,
                'synthetic': True
            })

        # Close FASTA handle
        if fasta_handle:
            fasta_handle.close()

        logger.info(f"Generated {len(mutations)} mutations for {region}")
        return mutations

    def inject_mutations_into_bam(self, bam_path: Path, reference_path: Path,
                                  mutations: List[Dict], output_bam: Path,
                                  region: str = None) -> bool:
        """
        Use BamSurgeon Docker to inject mutations into BAM file.

        Like GenomeNet, extracts a window BAM first, then injects mutations.
        This is much faster than processing the full BAM.

        Args:
            bam_path: Input BAM file path
            reference_path: Reference genome FASTA
            mutations: List of mutations to inject
            output_bam: Output mutated BAM path
            region: Genomic region (e.g., "chr20:20000000-25000000")

        Returns:
            True if successful, False otherwise
        """
        # Ensure paths are absolute
        bam_path = Path(bam_path).resolve()
        reference_path = Path(reference_path).resolve()
        output_bam = Path(output_bam).resolve()

        # Create output directory
        output_bam.parent.mkdir(parents=True, exist_ok=True)

        # Check if BAM index exists, create if not
        bam_index = bam_path.with_suffix('.bam.bai')
        if not bam_index.exists():
            bam_index_alt = Path(str(bam_path) + '.bai')
            if not bam_index_alt.exists():
                logger.info("Creating BAM index (required by BamSurgeon)...")
                index_cmd = f"""docker run --rm \
  -v {bam_path.parent}:/data/bams \
  biocontainers/samtools:v1.9-4-deb_cv1 \
  samtools index /data/bams/{bam_path.name}"""
                result = subprocess.run(index_cmd, shell=True, capture_output=True, text=True)
                if result.returncode != 0:
                    logger.error(f"Failed to create BAM index: {result.stderr}")
                    return False

        # Step 1: Extract window BAM from full BAM (like GenomeNet)
        # This makes BamSurgeon much faster and produces correct output
        window_bam = output_bam.parent / f"{output_bam.stem}_window.bam"
        if region:
            logger.info(f"Extracting window BAM for region {region}...")
            print(f"   Extracting window BAM for region {region}...", flush=True)

            # Verify BAM index exists and is accessible - CRITICAL for fast extraction
            if not bam_index.exists() and not bam_index_alt.exists():
                logger.error(f"BAM index not found at {bam_index} or {bam_index_alt}")
                print(f"   ERROR: BAM index missing! Extraction will be extremely slow.", flush=True)
                return False

            index_to_use = bam_index if bam_index.exists() else bam_index_alt
            index_size = index_to_use.stat().st_size
            logger.info(f"Using BAM index: {index_to_use} ({index_size} bytes)")
            print(f"   Using BAM index: {index_to_use.name} ({index_size / 1024:.1f} KB)", flush=True)

            # Warn if index seems too small (possibly corrupted)
            if index_size < 1000:
                logger.warning(f"BAM index is very small ({index_size} bytes) - may be corrupted")
                print(f"   WARNING: BAM index may be corrupted (only {index_size} bytes)", flush=True)

            # Use 4 threads for faster extraction of large BAM files
            extract_cmd = f"""docker run --rm \
  -v {bam_path.parent}:/data/bams \
  -v {output_bam.parent}:/data/output \
  biocontainers/samtools:v1.9-4-deb_cv1 \
  samtools view -@ 4 -b /data/bams/{bam_path.name} {region} -o /data/output/{window_bam.name}"""

            start_time = time.time()
            print(f"   Starting samtools extraction (timeout: 1800s / 30min)...", flush=True)
            result = subprocess.run(extract_cmd, shell=True, capture_output=True, text=True, timeout=1800)
            elapsed = time.time() - start_time
            print(f"   Extraction completed in {elapsed:.1f}s", flush=True)
            if result.returncode != 0:
                logger.error(f"Failed to extract window BAM: {result.stderr}")
                print(f"   Window extraction failed: {result.stderr}", flush=True)
                return False

            # Index the window BAM
            index_window_cmd = f"""docker run --rm \
  -v {output_bam.parent}:/data/output \
  biocontainers/samtools:v1.9-4-deb_cv1 \
  samtools index /data/output/{window_bam.name}"""

            result = subprocess.run(index_window_cmd, shell=True, capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(f"Failed to index window BAM: {result.stderr}")
                return False

            # Use window BAM for BamSurgeon
            bam_for_mutation = window_bam
            bam_mount_dir = output_bam.parent
            window_size_mb = window_bam.stat().st_size / (1024**2) if window_bam.exists() else 0
            print(f"   Window BAM extracted: {window_size_mb:.1f} MB", flush=True)
        else:
            bam_for_mutation = bam_path
            bam_mount_dir = bam_path.parent

        # Check if BWA index exists, create if not (required by BamSurgeon's internal bwa aln)
        # BWA index files: .amb, .ann, .bwt, .pac, .sa
        bwa_index_file = reference_path.with_suffix('.fa.bwt')
        if not bwa_index_file.exists():
            # Also check without .fa extension
            bwa_index_file_alt = Path(str(reference_path) + '.bwt')
            if not bwa_index_file_alt.exists():
                logger.info("Creating BWA index (required by BamSurgeon)...")
                print("   Creating BWA index (this may take a few minutes)...", flush=True)
                bwa_cmd = f"""docker run --rm \
  -v {reference_path.parent}:/data/reference \
  biocontainers/bwa:v0.7.17_cv1 \
  bwa index /data/reference/{reference_path.name}"""
                result = subprocess.run(bwa_cmd, shell=True, capture_output=True, text=True, timeout=600)
                if result.returncode != 0:
                    logger.error(f"Failed to create BWA index: {result.stderr}")
                    print(f"   BWA index failed: {result.stderr}", flush=True)
                    return False
                print("   BWA index created successfully", flush=True)

        # Create mutations file for BamSurgeon addsnv.py
        # BamSurgeon addsnv.py expects: chr start end VAF [alt]
        # - 4 columns: chr start(1-based) end(1-based) VAF (auto-selects alt base)
        # - 5 columns: chr start(1-based) end(1-based) VAF alt (specifies alt base)
        # NOTE: Do NOT include ref base - BamSurgeon reads it from reference FASTA
        mutations_bed = output_bam.parent / "mutations.bed"
        snp_count = 0
        with open(mutations_bed, 'w') as f:
            for mut in mutations:
                # Only include SNPs (single base ref and alt)
                if mut['ref'] in 'ACGT' and mut['alt'] in 'ACGT' and len(mut['ref']) == 1 and len(mut['alt']) == 1:
                    # BamSurgeon format: chr start(1-based) end(1-based) VAF alt
                    # Position is 1-based for BamSurgeon (not 0-based BED)
                    f.write(f"{mut['chrom']}\t{mut['pos']}\t{mut['pos']}\t0.5\t{mut['alt']}\n")
                    snp_count += 1

        logger.info(f"Created mutations BED with {snp_count} SNPs (from {len(mutations)} total mutations)")

        try:
            # Build BamSurgeon Docker command (matches GenomeNet step3.sh exactly)
            # CRITICAL: GenomeNet mounts everything under /work and uses -w /work
            # This ensures BamSurgeon's temp files and output are in the same location
            #
            # Key insight: BamSurgeon expects all files (BAM, reference, mutations, output)
            # to be accessible from its working directory. When paths are split across
            # multiple mounts, the final merge step can fail silently.
            #
            # Solution: Copy window BAM to output directory so everything is in /work

            # Copy window BAM to output directory if it's in a different location
            work_dir = output_bam.parent
            if bam_for_mutation.parent != work_dir:
                work_bam = work_dir / bam_for_mutation.name
                work_bai = work_dir / (bam_for_mutation.name + '.bai')
                shutil.copy2(bam_for_mutation, work_bam)
                # Also copy index if it exists
                src_bai = Path(str(bam_for_mutation) + '.bai')
                if src_bai.exists():
                    shutil.copy2(src_bai, work_bai)
                bam_for_mutation = work_bam
                logger.info(f"Copied window BAM to work directory: {work_bam}")

            # Copy reference to work directory (BamSurgeon needs BWA index files too)
            work_ref = work_dir / reference_path.name
            if not work_ref.exists():
                shutil.copy2(reference_path, work_ref)
                # Copy all index files
                for ext in ['.fai', '.amb', '.ann', '.bwt', '.pac', '.sa']:
                    idx_src = Path(str(reference_path) + ext)
                    if idx_src.exists():
                        shutil.copy2(idx_src, work_dir / (reference_path.name + ext))
                # Also copy .dict if exists
                dict_src = reference_path.with_suffix('.dict')
                if dict_src.exists():
                    shutil.copy2(dict_src, work_dir / dict_src.name)
                logger.info(f"Copied reference and indexes to work directory")

            # Now run BamSurgeon with everything in one mount (like GenomeNet)
            #
            # BamSurgeon addsnv.py flow:
            # 1. For each mutation, creates a .muts.bam with just the mutated reads
            # 2. At the end, merges all .muts.bam files with original BAM to create output
            #
            # The merge step uses: samtools merge -f output.bam original.bam *.muts.bam
            # This can fail silently if samtools has issues
            #
            # Detect platform for Docker compatibility
            platform_flag = "--platform linux/amd64" if platform.system() == "Darwin" else ""

            # STRATEGY: Use --skipmerge to let BamSurgeon create .muts.bam files,
            # then we do our own merge with samtools. This is more reliable than
            # BamSurgeon's internal merge which can fail silently.
            #
            # Available options from help: --tmpdir, --seed, --force, --ignorepileup, --skipmerge
            bamsurgeon_cmd = f"""docker run --rm {platform_flag} \
  -v {work_dir}:/work \
  -w /work \
  quay.io/biocontainers/bamsurgeon:1.4.1--pyhdfd78af_0 \
  addsnv.py \
    -v /work/{mutations_bed.name} \
    -f /work/{bam_for_mutation.name} \
    -r /work/{work_ref.name} \
    -o /work/{output_bam.name} \
    --mindepth 10 --maxdepth 2000 --minmutreads 3 \
    --tmpdir /work \
    --seed {self.current_seed or 42} \
    --force \
    --skipmerge"""

            logger.info(f"Running BamSurgeon to inject {len(mutations)} mutations (with --skipmerge)...")
            logger.debug(f"Command: {bamsurgeon_cmd}")
            print(f"   Running BamSurgeon (timeout: 1800s / 30min)...", flush=True)

            result = subprocess.run(bamsurgeon_cmd, shell=True, capture_output=True, text=True, timeout=1800)

            # Print any errors from BamSurgeon
            if result.returncode != 0:
                print(f"   BamSurgeon returned code {result.returncode}", flush=True)
                if result.stderr:
                    for line in result.stderr.split('\n')[:5]:
                        if line.strip():
                            print(f"   stderr: {line}", flush=True)

            # With --skipmerge, BamSurgeon creates .muts.bam files but doesn't merge
            # We need to merge them ourselves using samtools
            muts_bams = list(work_dir.glob("addsnv.*.muts.bam"))
            valid_muts = [m for m in muts_bams if m.stat().st_size > 10000]  # >10KB

            print(f"   BamSurgeon created {len(muts_bams)} .muts.bam files, {len(valid_muts)} valid", flush=True)

            if valid_muts and bam_for_mutation.exists():
                # Build list of BAM files to merge: original window BAM + all valid .muts.bam
                bams_to_merge = [f"/work/{bam_for_mutation.name}"] + [f"/work/{m.name}" for m in valid_muts]
                bam_list = " ".join(bams_to_merge)

                print(f"   Merging {len(bams_to_merge)} BAM files with samtools...", flush=True)

                # Use samtools merge to combine all BAMs (with 4 threads)
                merge_cmd = f"""docker run --rm {platform_flag} \
  -v {work_dir}:/work \
  -w /work \
  biocontainers/samtools:v1.9-4-deb_cv1 \
  samtools merge -@ 4 -f /work/{output_bam.name} {bam_list}"""

                merge_result = subprocess.run(merge_cmd, shell=True, capture_output=True, text=True, timeout=1800)

                if merge_result.returncode == 0 and output_bam.exists():
                    output_size = output_bam.stat().st_size
                    min_bam_size = 1024 * 1024  # At least 1MB

                    if output_size > min_bam_size:
                        output_size_mb = output_size / (1024**2)
                        logger.info(f"Successfully created mutated BAM: {output_bam} ({output_size_mb:.1f} MB)")
                        print(f"   Mutated BAM created: {output_size_mb:.1f} MB", flush=True)

                        # Index the merged BAM
                        print(f"   Indexing merged BAM...", flush=True)
                        index_cmd = f"""docker run --rm {platform_flag} \
  -v {work_dir}:/work \
  biocontainers/samtools:v1.9-4-deb_cv1 \
  samtools index /work/{output_bam.name}"""
                        subprocess.run(index_cmd, shell=True, capture_output=True, text=True)

                        # Sort the BAM with 4 threads (like GenomeNet step4.sh does)
                        sorted_bam = work_dir / f"{output_bam.stem}.sorted.bam"
                        print(f"   Sorting merged BAM...", flush=True)
                        sort_cmd = f"""docker run --rm {platform_flag} \
  -v {work_dir}:/work \
  biocontainers/samtools:v1.9-4-deb_cv1 \
  samtools sort -@ 4 /work/{output_bam.name} -o /work/{sorted_bam.name}"""
                        sort_result = subprocess.run(sort_cmd, shell=True, capture_output=True, text=True, timeout=1800)

                        if sort_result.returncode == 0 and sorted_bam.exists():
                            # Replace unsorted with sorted
                            os.unlink(output_bam)
                            shutil.move(sorted_bam, output_bam)
                            # Re-index
                            subprocess.run(index_cmd, shell=True, capture_output=True, text=True)
                            print(f"   BAM sorted and indexed", flush=True)

                        self._cleanup_bamsurgeon_temp_files(work_dir, output_bam.name, keep_output=True)
                        return True
                    else:
                        print(f"   Merge produced small file ({output_size} bytes)", flush=True)
                else:
                    print(f"   samtools merge failed: {merge_result.stderr[:200] if merge_result.stderr else 'unknown'}", flush=True)

            # If we get here, something went wrong
            print(f"   Files in work directory:", flush=True)
            for f in sorted(work_dir.iterdir(), key=lambda x: x.stat().st_size, reverse=True)[:15]:
                print(f"      {f.name}: {f.stat().st_size} bytes", flush=True)
            logger.error(f"Mutations file kept for inspection: {mutations_bed}")
            return False

        except subprocess.TimeoutExpired:
            logger.error("BamSurgeon timed out after 30 minutes")
            if mutations_bed.exists():
                os.unlink(mutations_bed)
            return False
        except Exception as e:
            logger.error(f"Error running BamSurgeon: {e}")
            if mutations_bed.exists():
                os.unlink(mutations_bed)
            return False

    def _cleanup_bamsurgeon_temp_files(self, work_dir: Path, output_name: str, keep_output: bool = True):
        """
        Clean up temporary files created by BamSurgeon.

        BamSurgeon creates many temp files: .muts.bam, logs, vcfs, etc.
        We want to keep only the final output BAM.
        """
        work_dir = Path(work_dir)

        # Patterns to clean up
        patterns_to_remove = [
            'mutations.bed',       # Our mutations input file
            '*_window.bam',        # Window BAM we extracted
            '*_window.bam.bai',    # Window BAM index
            '*.muts.bam',          # BamSurgeon per-mutation temp files
            'addsnv_logs_*',       # BamSurgeon log directories
            'addsnv.tmp*',         # BamSurgeon temp directories
            '*.addsnv.*.vcf',      # BamSurgeon output VCFs
            'chr20.fa',            # Copied reference
            'chr20.fa.*',          # Copied reference indexes
            'chr20.dict',          # Copied reference dict
        ]

        for pattern in patterns_to_remove:
            for match in work_dir.glob(pattern):
                try:
                    if match.is_dir():
                        shutil.rmtree(match)
                    else:
                        # Don't delete the output file
                        if keep_output and match.name == output_name:
                            continue
                        os.unlink(match)
                except Exception as e:
                    logger.debug(f"Could not remove {match}: {e}")

    def create_merged_truth_vcf(self, original_truth_vcf: Path, mutations: List[Dict],
                                 output_vcf: Path, reference_path: Path = None) -> bool:
        """
        Create a merged truth VCF containing original truth + synthetic mutations.

        Uses Docker for all bcftools/bgzip/tabix operations to avoid host dependency issues.
        Matches GenomeNet's step4.sh approach exactly.

        Args:
            original_truth_vcf: Path to original GIAB truth VCF
            mutations: List of synthetic mutations (from inject_mutations)
            output_vcf: Path for merged output VCF
            reference_path: Path to reference FASTA (unused, kept for API compat)

        Returns:
            True if successful, False otherwise
        """
        try:
            original_truth_vcf = Path(original_truth_vcf).resolve()
            output_vcf = Path(output_vcf).resolve()
            output_vcf.parent.mkdir(parents=True, exist_ok=True)

            # Create VCF from synthetic mutations (matching GenomeNet format)
            # Get chromosome and contig length from first mutation
            chrom = mutations[0].get('chrom', 'chr20') if mutations else 'chr20'
            # Default contig length for chr20 - could be passed as parameter if needed
            contig_length = 64444167

            synthetic_vcf = output_vcf.parent / f"{output_vcf.stem}_synthetic.vcf"
            with open(synthetic_vcf, 'w') as f:
                # VCF header - must include contig definition for bcftools
                f.write("##fileformat=VCFv4.2\n")
                f.write("##source=MutationInjector\n")
                f.write(f"##contig=<ID={chrom},length={contig_length}>\n")
                f.write("##INFO=<ID=SYNTHETIC,Number=0,Type=Flag,Description=\"Synthetic mutation injected by BamSurgeon\">\n")
                f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
                f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n")

                # Write mutations sorted by position
                for mut in sorted(mutations, key=lambda x: x['pos']):
                    chrom = mut['chrom']
                    pos = mut['pos']
                    ref = mut['ref']
                    alt = mut['alt']
                    # Write as heterozygous (0/1) with PASS filter
                    f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\tSYNTHETIC\tGT\t0/1\n")

            logger.info(f"Created synthetic mutations VCF with {len(mutations)} variants: {synthetic_vcf}")

            # Use Docker bcftools for all VCF operations (avoids host dependency issues)
            # Docker image: biocontainers/bcftools
            work_dir = output_vcf.parent
            truth_dir = original_truth_vcf.parent

            # Step 1: Sort and compress synthetic VCF
            sort_synthetic_cmd = f"""docker run --rm \
  -v {work_dir}:/work \
  biocontainers/bcftools:v1.9-1-deb_cv1 \
  bcftools sort /work/{synthetic_vcf.name} -Oz -o /work/{synthetic_vcf.name}.gz"""

            result = subprocess.run(sort_synthetic_cmd, shell=True, capture_output=True, text=True, timeout=120)
            if result.returncode != 0:
                logger.error(f"bcftools sort synthetic failed: {result.stderr}")
                return False

            synthetic_vcf_gz = synthetic_vcf.with_suffix('.vcf.gz')
            if not synthetic_vcf_gz.exists():
                # Docker might create .vcf.gz instead of adding .gz
                synthetic_vcf_gz = Path(str(synthetic_vcf) + '.gz')

            # Step 2: Index synthetic VCF
            index_synthetic_cmd = f"""docker run --rm \
  -v {work_dir}:/work \
  biocontainers/bcftools:v1.9-1-deb_cv1 \
  bcftools index /work/{synthetic_vcf_gz.name}"""

            result = subprocess.run(index_synthetic_cmd, shell=True, capture_output=True, text=True, timeout=60)
            if result.returncode != 0:
                logger.error(f"bcftools index synthetic failed: {result.stderr}")
                return False

            # Step 3: Concat original truth + synthetic (use -a for allow-overlaps)
            # Mount both directories if they're different
            tmp_merged = work_dir / "merged_tmp.vcf.gz"

            if truth_dir == work_dir:
                concat_cmd = f"""docker run --rm \
  -v {work_dir}:/work \
  biocontainers/bcftools:v1.9-1-deb_cv1 \
  bcftools concat -a /work/{original_truth_vcf.name} /work/{synthetic_vcf_gz.name} -Oz -o /work/merged_tmp.vcf.gz"""
            else:
                concat_cmd = f"""docker run --rm \
  -v {truth_dir}:/truth \
  -v {work_dir}:/work \
  biocontainers/bcftools:v1.9-1-deb_cv1 \
  bcftools concat -a /truth/{original_truth_vcf.name} /work/{synthetic_vcf_gz.name} -Oz -o /work/merged_tmp.vcf.gz"""

            result = subprocess.run(concat_cmd, shell=True, capture_output=True, text=True, timeout=120)
            if result.returncode != 0:
                logger.error(f"bcftools concat failed: {result.stderr}")
                return False

            # Step 4: Sort merged VCF
            sort_merged_cmd = f"""docker run --rm \
  -v {work_dir}:/work \
  biocontainers/bcftools:v1.9-1-deb_cv1 \
  bcftools sort /work/merged_tmp.vcf.gz -Oz -o /work/{output_vcf.name}"""

            result = subprocess.run(sort_merged_cmd, shell=True, capture_output=True, text=True, timeout=120)
            if result.returncode != 0:
                logger.error(f"bcftools sort merged failed: {result.stderr}")
                return False

            # Step 5: Index final merged VCF
            index_final_cmd = f"""docker run --rm \
  -v {work_dir}:/work \
  biocontainers/bcftools:v1.9-1-deb_cv1 \
  bcftools index /work/{output_vcf.name}"""

            result = subprocess.run(index_final_cmd, shell=True, capture_output=True, text=True, timeout=60)
            if result.returncode != 0:
                logger.error(f"bcftools index final failed: {result.stderr}")
                return False

            # Clean up temp files
            for temp_file in [synthetic_vcf, synthetic_vcf_gz,
                              Path(str(synthetic_vcf_gz) + ".csi"),
                              tmp_merged]:
                if temp_file.exists():
                    try:
                        os.unlink(temp_file)
                    except:
                        pass

            logger.info(f"Created merged truth VCF: {output_vcf}")
            return True

        except subprocess.TimeoutExpired:
            logger.error("VCF merging timed out")
            return False
        except Exception as e:
            logger.error(f"Error creating merged truth VCF: {e}")
            return False
