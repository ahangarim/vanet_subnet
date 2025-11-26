# VANET – Genomics Variant-Calling Subnet

VANET is the production-ready Bittensor subnet that rewards miners for accurately genotyping hidden 5 Mb windows on chr20. Validators synthesize benchmarks from GIAB truth data, inject biologically realistic mutations with BamSurgeon, distribute mutated BAMs to miners, and evaluate the returned VCFs with hap.py/rtg tools. This README walks through the entire pipeline, from environment provisioning to end-to-end validator/miner operation.

## Table of Contents
- [Repository Layout](#repository-layout)
- [System Prerequisites](#system-prerequisites)
- [Local Setup](#local-setup)
- [Reference Data Download](#reference-data-download)
- [Validator Execution Pipeline](#validator-execution-pipeline)
- [Miner Execution Pipeline](#miner-execution-pipeline)
- [Automation via Bittensor Neurons](#automation-via-bittensor-neurons)
- [Scoring, Leaderboards, and Anti-Cheating](#scoring-leaderboards-and-anti-cheating)
- [Monitoring & Troubleshooting](#monitoring--troubleshooting)
- [Additional Documentation](#additional-documentation)

---

## Repository Layout

```
vanet/
├── base/                     # Core subnet config + protocol definitions
│   ├── __init__.py           # Module exports
│   ├── config.py             # Standardized config dataclasses (net, data, miner, validator)
│   ├── genomics_config.py    # Window size, mutation counts, Docker image names, scoring weights
│   ├── protocol.py           # Task/response synapse dataclasses shared by neurons
│   └── s3_manifest.json      # Reference data paths (local + S3 URLs)
├── neurons/                  # Bittensor neuron entrypoints
│   ├── __init__.py           # Module exports
│   ├── validator.py          # Full validator orchestration loop
│   └── miner.py              # Reference miner (GATK by default)
├── utils/                    # Genomics utility modules
│   ├── __init__.py           # Module exports
│   ├── task_generation.py   # Random genomic task generator
│   ├── mutation_injection.py # BamSurgeon mutation injection
│   ├── scoring.py            # hap.py validation and advanced scoring
│   ├── weight_tracking.py    # EMA and innovation boost tracking
│   ├── gatk_utils.py         # GATK HaplotypeCaller utilities
│   └── file_utils.py         # S3 download and file management
├── benchmarking_honest_miners/  # Mutation-space analysis proving anti-cheat guarantees
│   ├── analyze_mutation_space.py
│   └── mutation_space_scenarios.py
├── docs/                     # Architecture notes, hap.py Docker instructions, tooling background
├── datasets/                 # (gitignored) GIAB assets downloaded locally
│   ├── bams/                 # 300x chr20 BAM files
│   ├── reference/            # GRCh38 chr20 FASTA + indexes
│   └── truth/                # GIAB truth VCF + BED masks
├── __init__.py               # Package initialization
├── requirements.txt          # Python dependencies for both miners and validators
├── .env.example              # Example environment configuration
├── .gitignore                # Git ignore rules
└── README.md                 # This document
```

Each validator/miner launch uses `base/config.py` to pull environment variables (or CLI args) into a consistent config, then hands those settings to `neurons/*.py`.

---

## System Prerequisites

| Component | Requirement | Notes |
|-----------|-------------|-------|
| OS | Linux (Ubuntu 20.04+), macOS 13+, or WSL2 | Docker + Bittensor run best on Linux; macOS requires Rosetta for amd64 containers |
| CPU/RAM (Validator) | ≥4 cores / 16 GB RAM | Hap.py + BamSurgeon benefit from more cores; allocate ≥100 GB disk for datasets |
| CPU/RAM (Miner) | ≥4 cores / 8 GB RAM | Variant calling scales with threads; SSD recommended |
| GPU | Not required | DeepVariant will use CPU unless you supply CUDA image |
| Docker | 24.0+ | Needed for GATK, BamSurgeon, bcftools, hap.py |
| Python | 3.10+ (we test on 3.12) | Use virtualenv per instructions below |
| Bittensor | Latest pip install of `bittensor` | Provides wallet/subtensor/axon/dendrite APIs |

Optional developer tooling: `ruff`, `pytest`, `boto3` (for private S3 buckets).

---

## Local Setup

1. **Clone & virtual environment**
   ```bash
   git clone https://github.com/<org>/vanet.git
   cd vanet
   python3 -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip
   pip install -r requirements.txt   # install bittensor + genomics deps
   ```

2. **Install Bittensor CLI (optional but recommended)**
   ```bash
   pip install bittensor  # already included in requirements but safe to confirm
   btcli --version
   ```

3. **Install bioinformatics utilities (host-level)**
   ```bash
   # Ubuntu/Debian
   sudo apt-get update && sudo apt-get install -y samtools bcftools tabix parallel

   # macOS (Homebrew)
   brew install samtools bcftools htslib gnu-parallel
   ```
   > BamSurgeon, GATK, and hap.py are executed via Docker images, so no additional local installs are required beyond Docker.

4. **Docker & Rosetta validation**
   ```bash
   docker run --rm hello-world
   docker run --rm --platform linux/amd64 alpine uname -m  # expect x86_64 even on Apple Silicon
   docker pull --platform linux/amd64 mgibio/hap.py:v0.3.12
   docker pull broadinstitute/gatk:4.5.0.0
   docker pull quay.io/biocontainers/bamsurgeon:1.4.1--pyhdfd78af_0
   ```
   Apple Silicon: install Rosetta (`softwareupdate --install-rosetta`) and enable “Use Rosetta for x86/amd64 emulation” in Docker Desktop.

5. **Configure environment variables**
   ```bash
   cp .env.example .env
   # Edit using your favorite editor
   nano .env
   ```
   Key fields for validators:
   ```bash
   NETUID=2
   SUBTENSOR_NETWORK=test
   WALLET_NAME=my_validator
   WALLET_HOTKEY=default
   VALIDATOR_INTERVAL=240        # minutes between new tasks
   VARIANT_TIMEOUT=3600          # miner SLA in seconds
   USE_MUTATIONS=true            # enable BamSurgeon injection
   NUM_MUTATIONS=10              # synthetic SNPs per window
   ```
   Key fields for miners:
   ```bash
   NETUID=2
   SUBTENSOR_NETWORK=test
   WALLET_NAME=my_miner
   WALLET_HOTKEY=default
   MINER_TOOL=gatk               # or deepvariant/custom
   MINER_THREADS=4
   MINER_TIMEOUT=3600
   ```

6. **(Optional) Developer extras**
   ```bash
   pip install -e .[dev]
   ruff check .
   pytest
   ```

---

## Reference Data Download

Validators and miners must share the same chr20 datasets. Pre-trimmed GIAB assets live in our public S3 bucket (`base/s3_manifest.json` lists canonical URLs). Run the following once per host:

```bash
mkdir -p datasets/{bams,reference,truth}

# 300x chr20 BAM + index
ewget -O datasets/bams/sample.GRCh38.300x_chr20.bam \
  https://vanetdata.s3.us-east-1.amazonaws.com/bams/sample.GRCh38.300x_chr20.bam
wget -O datasets/bams/sample.GRCh38.300x_chr20.bam.bai \
  https://vanetdata.s3.us-east-1.amazonaws.com/bams/sample.GRCh38.300x_chr20.bam.bai

# GRCh38 chr20 FASTA + metadata
wget -O datasets/reference/chr20.fa \
  https://vanetdata.s3.us-east-1.amazonaws.com/reference/chr20.fa
wget -O datasets/reference/chr20.fa.fai \
  https://vanetdata.s3.us-east-1.amazonaws.com/reference/chr20.fa.fai
wget -O datasets/reference/chr20.dict \
  https://vanetdata.s3.us-east-1.amazonaws.com/reference/chr20.dict

# SDF template for rtg vcfeval (hap.py engine)
wget -O datasets/reference/chr20.sdf.tar.gz \
  https://vanetdata.s3.us-east-1.amazonaws.com/reference/chr20.sdf.tar.gz
tar -xzf datasets/reference/chr20.sdf.tar.gz -C datasets/reference/
rm datasets/reference/chr20.sdf.tar.gz

# Truth VCF + index + confident regions BED
wget -O datasets/truth/sample_GRCh38_1_22_v4.2.1_benchmark.chr20.vcf.gz \
  https://vanetdata.s3.us-east-1.amazonaws.com/truth/sample_GRCh38_1_22_v4.2.1_benchmark.chr20.vcf.gz
wget -O datasets/truth/sample_GRCh38_1_22_v4.2.1_benchmark.chr20.vcf.gz.csi \
  https://vanetdata.s3.us-east-1.amazonaws.com/truth/sample_GRCh38_1_22_v4.2.1_benchmark.chr20.vcf.gz.csi
wget -O datasets/truth/sample_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.chr20.bed \
  https://vanetdata.s3.us-east-1.amazonaws.com/truth/sample_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.chr20.bed
```

The validator caches per-task artifacts under `output/` (`mutated_bams/`, `merged_truth/`). Old files are automatically pruned after 24 h.

---

## Validator Execution Pipeline

The validator neuron mirrors the seven-stage manual flow described in `docs/architecture.md` and the original `step1.sh`–`step7`. Below is the production variant (all steps orchestrated inside `neurons/validator.py`).

1. **Window sampling** – `GenomicsTaskGenerator` (utils/genomics.py) samples a random 5 Mb interval within bounds defined in `base/s3_manifest.json`. BAMs are sliced via `samtools view` (Dockerized) and stored under `output/mutated_bams/<task_id>.bam` with clean `@RG` headers. BED masks and truth VCFs are subsetted using bcftools (see `subset_bed()` / `slice_truth_vcf()` for parity with manual `generate_chr20_windows.py`).

2. **Synthetic mutation selection** – `MutationInjector` calls BamSurgeon’s `randomsites.py` to pick biologically plausible SNPs inside the sampled window. Counts and seed values are recorded inside the task metadata to ensure reproducibility.

3. **Mutation application** – BamSurgeon `addsnv.py` injects the selected SNPs into the per-window BAM. The mutated BAM is sorted/indexed, and the intermediate VCF capturing inserted alleles is produced. Temporary directories (`addsnv.tmp`, logs) are cleaned automatically.

4. **Truth packaging** – The injected VCF is merged with the sliced GIAB truth (bcftools concat → sort) to create `*_truth_merged.vcf.gz` plus an updated BED mask capturing confident loci that intersect the window. These files become the ground truth bundle stored in `output/merged_truth/`.

5. **Task distribution** – The validator writes a `GenomicsTaskSynapse` containing paths (or presigned URLs) to the mutated BAM, reference FASTA/SDF, and mask. It pushes the task to a subset of miners (batch size defined in `VALIDATOR_CONFIG[
query_batch_size]`). Weighting/prioritization leverages metagraph stake information.

6. **Miner SLA & result ingestion** – Each miner receives the BAM/region, runs its caller, and returns a VCF. The validator enforces `variant_calling_timeout` (default 3600 s). Late or malformed submissions incur penalties.

7. **hap.py scoring & weight updates** – `HappyScorer` wraps the exact Docker command from `scripts/run_hap_validation.py`, translating host paths into `/work` bindings and calling:
   ```bash
   docker run --rm --platform linux/amd64 \
     -v $REPO:/work \
     mgibio/hap.py:v0.3.12 \
     /opt/hap.py/bin/hap.py TRUTH.vcf.gz SAMPLE.vcf \
       -r /work/datasets/reference/chr20.fa \
       -o /work/output/happy/<task_id> \
       --engine vcfeval --engine-vcfeval-template /work/datasets/reference/chr20.sdf \
       --target-regions /work/output/masks/<task_id>.bed
   ```
   The resulting `benchmark.summary.csv` feeds into `AdvancedScorer`, which replicates the `score_hap_results.py` formula (weighted SNP/INDEL F1, completeness, FP penalties, Ti/Tv + het/hom sanity). Scores are fed into an exponential moving average with an innovation boost. Updated weights are written back to the subnet at the configured interval.

8. **Cleanup** – Any mutated BAMs/truth bundles older than 24 h are deleted to reclaim disk (`Validator._cleanup_old_files`).

> **Manual dry-runs:** You can still run the original scripts (`scripts/generate_chr20_windows.py`, `run_hap_validation.py`, etc.) for debugging. They now live in the `docs` folder and match the automated pipeline.

### Running the validator neuron

```bash
source .venv/bin/activate
python -m vanet.neurons.validator \
  --netuid 2 \
  --subtensor.network test \
  --wallet.name my_validator \
  --wallet.hotkey default \
  --logging.debug
```

The neuron prints a 10-step initialization checklist (config, wallet, subtensor, metagraph, dendrite, scoring objects, etc.) before entering the task loop.

---

## Miner Execution Pipeline

Miners implement the inverse workflow: accept BAMs, run variant callers, and return VCFs.

1. **Task reception** – The validator contacts miners via dendrite using the `GenomicsTaskSynapse`. Each miner exposes `forward_genomics`, `forward_compute`, and `forward_task` on its axon; tasks are routed to `forward_genomics`.

2. **Variant calling** – `neurons/miner.py` reads the requested tool from config (`MINER_TOOL`). The default is Dockerized GATK HaplotypeCaller:
   ```bash
   docker run --rm \
     -v $PWD/datasets:/datasets \
     -v /tmp:/tmp \
     broadinstitute/gatk:4.5.0.0 \
     gatk HaplotypeCaller \
       -I /datasets/bams/<mutated>.bam \
       -R /datasets/reference/chr20.fa \
       -L chr20:start-end \
       -O /tmp/output.vcf.gz \
       -ERC GVCF
   ```
   The reference miner converts GVCF → VCF (bcftools) before responding. Miners can swap in DeepVariant or proprietary callers by extending `setup_variant_caller()`.

3. **Result caching** – To conserve bandwidth, miners hash the task payload (BAM path + region). Repeated tasks return cached VCFs instantly if `MINER_CONFIG["cache_results"]` is enabled.

4. **Response** – Successful calls return the VCF bytes plus metadata (runtime, caller version). Failures send structured error messages so validators can penalize appropriately.

### Running the miner neuron

```bash
source .venv/bin/activate
python -m vanet.neurons.miner \
  --netuid 2 \
  --subtensor.network test \
  --wallet.name my_miner \
  --wallet.hotkey default \
  --axon.port 8091 \
  --logging.debug
```

The miner automatically registers if needed, spins up an axon, and advertises supported variant callers. Ensure Docker has permission to mount the `datasets/` directory containing the validator-provided BAMs.

---

## Automation via Bittensor Neurons

- **Validator neuron** handles: task scheduling (`VALIDATOR_CONFIG['query_batch_size']`), miner selection, SLA tracking, hap.py scoring, EMA/boost updates, and weight writes.
- **Miner neuron** exposes axon endpoints, prioritizes requests (stake-aware priority functions), validates payloads (blacklist functions), and streams artifacts back to validators.

Both neurons rely on environment variables or CLI flags for wallet, subtensor endpoint, axon IP/port, logging verbosity, and Docker toggles. Consult `bt.wallet.add_args`, `bt.subtensor.add_args`, and `bt.axon.add_args` inside each neuron for advanced overrides (coldkey path, external IP, etc.).

---

## Scoring, Leaderboards, and Anti-Cheating

### Scoring formula

`utils/genomics.AdvancedScorer` mirrors `scripts/score_hap_results.py` (carried over from `genomenet`). Key rules:
- **Core component (60 %)** – Weighted harmonic mean of SNP/INDEL F1 from hap.py PASS rows (weights based on truth counts) transformed with `gamma=4` to emphasize deltas near 1.0.
- **Completeness (15 %)** – Average recall (`METRIC.Recall`) and PASS coverage (`1 - METRIC.Frac_NA`).
- **False-positive penalty (15 %)** – Exponential penalty if FP rate exceeds 0.2% or if `QUERY.TOTAL / TRUTH.TOTAL` deviates beyond ±5%.
- **Quality sanity (10 %)** – Ti/Tv and het/hom ratios must stay close to truth slices (ratios deviating >0.1 decay the score).

Scores feed an EMA (`alpha=0.1`) plus an innovation boost (15 % for miners beating their previous best by ≥5 %). Validators can export per-round summaries using the standalone `scripts/scoreboard.py` for public leaderboards.

### Anti-cheating via mutation-space explosion

Every validation round uses a new random window and a fresh BamSurgeon mutation set. `benchmarking_honest_miners/README.md` details why brute-force memorization is infeasible: with only ten injected SNVs per 5 Mb window, the variant space already exceeds 10^61 combinations; doubling to twenty SNVs pushes beyond 10^124. Coupled with rotating windows and GIAB masks, miners must genuinely run variant callers to succeed.

---

## Monitoring & Troubleshooting

| Symptom | Likely Cause | Fix |
|---------|--------------|-----|
| `docker: permission denied` | User not in docker group | `sudo usermod -aG docker $USER && newgrp docker` |
| Validator logs `No manifest found` | `base/s3_manifest.json` missing | Ensure repo copied `docs/` & `benchmarking_honest_miners/`; manifest ships in base/ |
| Hap.py exits with `No sample name provided` | Miner returned multi-sample VCF | Ensure miners emit single-sample VCF (`bcftools reheader -s <(echo sample) ...`) |
| GATK OOM | Not enough RAM / threads too high | Reduce `MINER_THREADS`, increase Docker memory limit |
| Slow hap.py scoring | SDF missing | Confirm `datasets/reference/chr20.sdf` exists; rerun download step |
| Validator not registered | Wallet missing stake/registration | Run `btcli subnets list` & `btcli subnets register` before launching |

Logs are emitted via `bt.logging`. Use `--logging.trace` for verbose debugging.

---

## Additional Documentation

- `docs/architecture.md` – Deep dive into the miner/validator incentive loop, anti-cheat logic, and subnet economics.
- `docs/hap_py_docker.md` – Instructions for rebuilding a modern hap.py Docker image (`genomenet/hap-py:0.3.15`) if `mgibio/hap.py` ever disappears.
- `docs/tools_and_data.md` – Inventory of scripts (`scripts/generate_chr20_windows.py`, `run_hap_validation.py`, etc.) and how they map to the automated pipeline.
- `benchmarking_honest_miners/README.md` – Formal combinatorial analysis of mutation search space.
  - Generated figures (`mutation_space_windows.png`, `mutation_space_heatmap.png`, `mutation_space_entropy.png`) quantify cheating odds across 1/5/10 Mb windows and SNV/INDEL mixes using `mutation_space_scenarios.py`.

Need a quick standalone dry-run? Use the transferred scripts in `docs/tools_and_data.md` to reproduce any stage outside the Bittensor network (ideal for integration tests and CI).

---

## License & Support

- License: [Add license here]
- Issues: submit to this repository with logs + software versions
- Community: join the project’s Discord or Bittensor forums for validator/miner coordination

Happy mining and validating! 🚀
