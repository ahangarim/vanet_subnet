# VANET Utils

This folder contains all the genomics processing tools used by miners and validators. Think of these as the toolbox that makes everything work.

---

## Overview

The utils folder is organized into 6 focused modules, each handling a specific part of the genomics pipeline:

```
utils/
├── task_generation.py       # Creates random genomic challenges
├── mutation_injection.py    # Adds synthetic mutations to BAM files
├── scoring.py               # Validates and scores VCF results
├── weight_tracking.py       # Tracks miner performance over time
├── file_utils.py            # Downloads and manages data files
└── gatk_utils.py            # Runs GATK variant calling
```

---

## Module Details

### 1. Task Generation ([task_generation.py](task_generation.py))

**What it does:** Creates random genomic tasks for validators to send to miners.

**Main Class: `GenomicsTaskGenerator`**

Think of this as a random challenge generator. It:
- Picks a random 5MB region on chromosome 20
- Creates a unique task ID with timestamp
- Provides paths to reference data (BAM files, truth VCFs, etc.)

**Example:**
```python
from utils import GenomicsTaskGenerator

generator = GenomicsTaskGenerator(manifest_path="base/s3_manifest.json")
task = generator.generate_task()

# Returns something like:
# {
#   "task_id": "task_chr20_25000000_30000000_1234567890_4567",
#   "regions": ["chr20:25000000-30000000"],
#   "ref_build": "GRCh38",
#   "input_type": "BAM",
#   "output_format": "VCF"
# }
```

**Key Features:**
- Random region generation (between 10MB and 55MB on chr20)
- Always 5MB windows to keep tasks consistent
- Unique task IDs to prevent duplicates
- Loads paths from manifest (local or S3)

---

### 2. Mutation Injection ([mutation_injection.py](mutation_injection.py))

**What it does:** Adds synthetic mutations to BAM files using BamSurgeon.

**Main Class: `MutationInjector`**

This is the anti-cheating mechanism. It:
- Creates random SNP mutations (A→T, C→G, etc.)
- Reads the actual reference genome to get correct base pairs
- Injects mutations into BAM files using Docker
- Creates a "truth VCF" showing what mutations were added

**Example:**
```python
from utils import MutationInjector

injector = MutationInjector(seed=12345)

# Step 1: Generate mutations
mutations = injector.inject_mutations(
    task_id="task_123",
    region="chr20:20000000-25000000",
    num_mutations=50,
    reference_path="datasets/reference/chr20.fa"
)
# Returns: [{'chrom': 'chr20', 'pos': 20123456, 'ref': 'A', 'alt': 'T'}, ...]

# Step 2: Inject into BAM file
success = injector.inject_mutations_into_bam(
    bam_path="input.bam",
    reference_path="chr20.fa",
    mutations=mutations,
    output_bam="mutated.bam",
    region="chr20:20000000-25000000"
)

# Step 3: Create truth VCF
injector.create_merged_truth_vcf(
    original_truth_vcf="original_truth.vcf.gz",
    mutations=mutations,
    output_vcf="merged_truth.vcf.gz"
)
```

**How It Works:**
1. Generates random positions in the region
2. Reads actual reference bases from FASTA file (requires pysam)
3. Picks random alternate alleles (different from reference)
4. Uses BamSurgeon Docker to modify the BAM file
5. Merges synthetic + original variants into truth VCF

**Important:** Reseeds the random number generator each round so miners can't memorize mutations!

---

### 3. Scoring ([scoring.py](scoring.py))

**What it does:** Validates miner VCF files against ground truth using hap.py.

**Main Classes:**
- `HappyScorer` - Runs hap.py validation via Docker
- `AdvancedScorer` - Computes advanced scores with multiple components

**HappyScorer Example:**
```python
from utils import HappyScorer

scorer = HappyScorer()

results = scorer.score_vcf(
    truth_vcf="merged_truth.vcf.gz",
    query_vcf="miner_output.vcf",
    reference_fasta="chr20.fa",
    confident_bed="confident_regions.bed"
)

# Returns:
# {
#   'f1_snp': 0.95,           # SNP accuracy
#   'f1_indel': 0.88,         # INDEL accuracy
#   'precision_snp': 0.96,    # SNP precision
#   'recall_snp': 0.94,       # SNP recall
#   'weighted_f1': 0.933      # Combined score (70% SNP + 30% INDEL)
# }
```

**AdvancedScorer Example:**
```python
from utils import AdvancedScorer

# Takes hap.py metrics and computes advanced score
score = AdvancedScorer.compute_advanced_score(results)
# Returns: 87.5 (on 0-100 scale)
```

**Scoring Components:**
1. **Core Score (50%)**: Weighted F1 (70% SNP + 30% INDEL)
2. **Completeness (30%)**: Recall and coverage
3. **Quality (20%)**: Precision and biological realism

**Helper Functions:**
- `subset_bed()` - Filters BED file to specific region
- `slice_truth_vcf()` - Extracts VCF region using bcftools

---

### 4. Weight Tracking ([weight_tracking.py](weight_tracking.py))

**What it does:** Tracks miner performance over time with EMA and innovation rewards.

**Main Class: `ScoreTracker`**

This manages how miners are rewarded:
- Smooths scores over time (prevents lucky one-time scores)
- Rewards breakthrough performance with temporary boosts
- Normalizes scores into weights for the blockchain

**Example:**
```python
from utils import ScoreTracker

tracker = ScoreTracker(
    num_miners=10,
    alpha=0.1,          # EMA smoothing (10% weight on new scores)
    boost_factor=0.15   # 15% boost for innovation
)

# Update a miner's score
import time
timestamp = time.time()
effective_score = tracker.update(
    miner_id=3,
    new_score=0.92,
    timestamp=timestamp
)

# Get weights for blockchain
weights = tracker.get_normalized_weights()
# Returns: [0.05, 0.08, 0.12, 0.25, ...] (sums to 1.0)
```

**How It Works:**
1. **EMA Smoothing**: New score averaged with historical scores
   - `new_ema = 0.9 * old_ema + 0.1 * new_score`
2. **Innovation Detection**: If score beats global best by 5%, activate boost
3. **Boost Decay**: Boost decreases over 1 hour (half-life)
4. **Normalization**: Convert to weights summing to 1.0

**Why This Matters:**
- Prevents miners from gaming the system with one lucky score
- Rewards sustained good performance
- Encourages innovation and improvement

---

### 5. File Utils ([file_utils.py](file_utils.py))

**What it does:** Downloads and manages genomic data files.

**Main Functions:**
- `download_file()` - Downloads from HTTP/S3 with caching
- `ensure_s3_file()` - Ensures file exists locally, downloads if needed

**Example:**
```python
from utils import download_file, ensure_s3_file

# Download a file
download_file(
    url="https://example.com/data.bam",
    local_path="datasets/bams/data.bam",
    use_cache=True  # Skip if already exists
)

# Auto-download from S3 if needed
bam_path = ensure_s3_file(
    file_path="datasets/bams/sample.bam",
    file_type="bam"
)
# If file doesn't exist locally, downloads from S3
# Also downloads .bai index automatically
```

**Features:**
- Automatic caching (won't re-download existing files)
- Auto-downloads index files (.bai for BAM, .fai/.dict for reference)
- Supports both local paths and S3 URLs
- Creates directories as needed

---

### 6. GATK Utils ([gatk_utils.py](gatk_utils.py))

**What it does:** Runs GATK HaplotypeCaller via Docker.

**Main Functions:**
- `run_gatk_docker()` - Runs variant calling
- `ensure_reference_indexes()` - Creates required .dict and .fai files

**Example:**
```python
from utils import run_gatk_docker, ensure_reference_indexes

# Make sure reference is indexed
ensure_reference_indexes("datasets/reference/chr20.fa")

# Run GATK HaplotypeCaller
success = run_gatk_docker(
    bam_path="input.bam",
    reference_path="chr20.fa",
    output_vcf="output.vcf",
    interval="chr20:10000000-15000000",  # Optional region
    threads=4
)

if success:
    print("Variant calling complete!")
```

**What It Does:**
1. Creates reference dictionary (.dict) if missing
2. Creates reference index (.fai) if missing
3. Runs GATK in Docker container
4. Returns VCF with all variants found
5. Counts variants and logs results

**Docker Command Used:**
```bash
docker run --rm \
  -v /path/to/bams:/data/bams \
  -v /path/to/reference:/data/reference \
  -v /path/to/output:/data/output \
  broadinstitute/gatk:4.5.0.0 \
  gatk HaplotypeCaller \
    -R /data/reference/chr20.fa \
    -I /data/bams/input.bam \
    -O /data/output/output.vcf \
    -L chr20:10000000-15000000
```

---

## How They Work Together

Here's how all these modules collaborate in a typical validator round:

```
1. GenomicsTaskGenerator
   ↓ (creates random task)

2. MutationInjector
   ↓ (adds synthetic mutations)

3. Send to Miners → Miners use gatk_utils.run_gatk_docker()
   ↓ (miners return VCF)

4. HappyScorer
   ↓ (validates VCF against truth)

5. AdvancedScorer
   ↓ (computes final score)

6. ScoreTracker
   ↓ (updates EMA + boost)

7. Blockchain
   (weights submitted)
```

---

## Dependencies

Each module has specific dependencies:

| Module | Requires |
|--------|----------|
| task_generation.py | ✅ None (pure Python) |
| mutation_injection.py | 🐳 Docker (BamSurgeon, samtools, bcftools)<br>📦 pysam (for reading FASTA) |
| scoring.py | 🐳 Docker (hap.py, bcftools) |
| weight_tracking.py | 📦 numpy |
| file_utils.py | 📦 boto3 (optional, for S3)<br>📦 urllib (standard lib) |
| gatk_utils.py | 🐳 Docker (GATK, samtools) |

---

## Import Guide

All modules are exported from the utils package:

```python
# Task generation
from utils import GenomicsTaskGenerator

# Mutation injection
from utils import MutationInjector

# Scoring
from utils import HappyScorer, AdvancedScorer, subset_bed, slice_truth_vcf

# Weight tracking
from utils import ScoreTracker

# File utilities
from utils import download_file, ensure_s3_file

# GATK utilities
from utils import run_gatk_docker, ensure_reference_indexes
```

---

## Common Workflows

### Validator Workflow:
```python
from utils import (
    GenomicsTaskGenerator,
    MutationInjector,
    HappyScorer,
    ScoreTracker
)

# 1. Generate task
generator = GenomicsTaskGenerator(manifest_path="base/s3_manifest.json")
task = generator.generate_task()

# 2. Inject mutations
injector = MutationInjector()
mutations = injector.inject_mutations(...)
injector.inject_mutations_into_bam(...)
injector.create_merged_truth_vcf(...)

# 3. Send to miners (via dendrite)
# ...miners process and return VCFs...

# 4. Score results
scorer = HappyScorer()
results = scorer.score_vcf(truth_vcf, miner_vcf, reference, bed)

# 5. Update weights
tracker = ScoreTracker(num_miners=100)
tracker.update(miner_id=0, new_score=results['weighted_f1'], timestamp=time.time())
weights = tracker.get_normalized_weights()
```

### Miner Workflow:
```python
from utils import run_gatk_docker, download_file

# 1. Receive task from validator
# task = synapse.task_id, synapse.data_locator, synapse.regions

# 2. Download BAM if needed
download_file(task['data_locator']['uri'], local_bam_path)

# 3. Run GATK
success = run_gatk_docker(
    bam_path=local_bam_path,
    reference_path="chr20.fa",
    output_vcf="output.vcf",
    interval=task['regions'][0]
)

# 4. Return VCF content to validator
```

---

## Performance Tips

1. **Use caching**: Enable `use_cache=True` in download functions
2. **Limit mutations**: More mutations = slower BamSurgeon processing
3. **Use regions**: Specify genomic intervals to avoid processing whole genome
4. **Docker resources**: Allocate enough memory for GATK (4GB minimum)
5. **Parallel processing**: BamSurgeon and GATK support multi-threading

---

## Troubleshooting

**BamSurgeon Issues:**
- "pysam not installed" → `pip install pysam`
- "Mutations may fail" → Make sure pysam is reading actual ref bases
- "Merge produced small file" → Check that .muts.bam files were created

**hap.py Issues:**
- "No summary CSV created" → Check that VCF files are valid
- "SDF not found" → Create SDF directory with rtg format command
- "Return code 1" → Often just warnings, check if CSV was still created

**GATK Issues:**
- "Reference dict not found" → Run `ensure_reference_indexes()` first
- "Timeout" → Increase timeout or reduce region size
- "Out of memory" → Increase Docker memory allocation

---

## Learn More

- See `neurons/` folder for how these modules are used by miner and validator
- See `base/` folder for configuration options
- Check individual module files for detailed docstrings
