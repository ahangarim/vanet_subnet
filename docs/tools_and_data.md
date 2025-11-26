# Tooling & Data Assets

This document enumerates the external tools and datasets required to stand up the validator/miner prototype. Every asset is tracked with a version, canonical download location, approximate size, and local target path so the environment can be reproduced deterministically.

---

## 1. Tools (Docker-first)

| Name | Version | Docker Image | Purpose |
| --- | --- | --- | --- |
| **hap.py** | 0.3.12 | `mgibio/hap.py:v0.3.12` | Validator scoring loop (hap.py + vcfeval); same image used in `nist/run_snake_sr.sh` (see `docs/hap_py_docker.md`). |
| **RTG Tools** | 3.12.1 | `realtimegenomics/rtg-tools:3.12.1` | Direct vcfeval CLI / supplementary utilities. |
| **BamSurgeon** | latest | `fac2003/bamsurgeon:latest` | Inject synthetic variants for hidden eval pool. |
| **GATK** | 4.5.0.0 | `broadinstitute/gatk:4.5.0.0` | Reference miner pipeline (HaplotypeCaller). |
| **DeepVariant** | 1.6.1 | `google/deepvariant:1.6.1` | GPU/CPU miner baseline. |
| **samtools / htslib suite** | 1.22.1 | host install (`brew install samtools`) | Slicing BAMs, indexing FASTA/VCFs. |

We standardize on Docker images for every heavy tool so validator/miner nodes can be reproduced exactly, audited easily, and managed by Watchtower (or similar) for automated updates. For production you can run Watchtower against these images to ensure they stay patched without manual intervention.

### 1.1 Pulling the images

Use the helper script to pull everything at once:

```bash
python scripts/fetch_assets.py --manifest assets/tools.json --include-docker
```

Or pull manually:

```bash
docker pull pkrusche/hap.py:0.3.15
docker pull realtimegenomics/rtg-tools:3.12.1
docker pull fac2003/bamsurgeon:latest
docker pull broadinstitute/gatk:4.5.0.0
docker pull google/deepvariant:1.6.1
```

To enable automated image refreshes, attach Watchtower to your Docker daemon and whitelist these repository names so validators/miners receive updates without downtime. We rely on `mgibio/hap.py:v0.3.12`, which is actively maintained and already powers the `nist/run_snake_sr.sh` workflows—no schema-v1 issues like the legacy `pkrusche/hap.py`.

---

## 2. Datasets (Initial 10 Mb Working Set)

| Asset | Description | Source | Target |
| --- | --- | --- | --- |
| **GRCh38 chr20 FASTA** | Reference sequence (chr20 only) for tractable development. | `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz` | `datasets/reference/GRCh38_chr20.fa` (+ `.fai`) |
| **sample truth VCF (chr20:10–20 Mb)** | Ground-truth variants limited to development region. | `tabix -h ftp://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/.../sample_GRCh38_1_22_v4.2.1_benchmark.vcf.gz chr20:10000000-20000000` | `datasets/giab/truth/sample_chr20_10M_20M_truth.vcf.gz` (+ `.tbi`) |
| **sample high-confidence BED** | Region mask for scoring (full bed + optional regional slice). | `https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/.../sample_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed` | `datasets/giab/truth/sample_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed` |
| **sample BAM (chr20:10–20 Mb)** | Miner input reads (10 Mb window). Extracted from novoalign 300× chr20 BAM. | `samtools view -b ftp://.../sample.GRCh38.300x_chr20.bam chr20:10000000-20000000 -o datasets/giab/bams/sample_chr20_10M_20M.bam` | `datasets/giab/bams/sample_chr20_10M_20M.bam` (+ `.bai`) |

Future phases will expand the hidden evaluation pool via BamSurgeon-derived BAMs, additional chromosomes, and multiple genomes (HG003/4, HG005/6/7).

---

## 3. Tracking & Automation

- Structured manifests (`assets/tools.json`, `assets/data.json`) enumerate every artifact with URLs, hashes (when available), and target paths.
- `scripts/fetch_assets.py` consumes those manifests to download/copy assets idempotently.
- Large assets (Docker images, >10 GB BAMs) are referenced but not stored in git; `.gitignore` excludes `tools/` and `datasets/`.
- Each download script logs checksum verification so validators/miners can trust their environment.

Run order for a fresh machine:

```bash
# 1. Install native dependencies (samtools, curl, docker, python >=3.10)
brew install samtools

# 2. Fetch tools archives/repos
python scripts/fetch_assets.py --manifest assets/tools.json --category tools

# 3. Fetch reference + truth + sample BAM slices
python scripts/fetch_assets.py --manifest assets/data.json --category data
```

Once these steps complete you have:

- `tools/` with hap.py, RTG Tools, GATK zip, BamSurgeon repo ready for setup.
- `datasets/` with chr20 FASTA/index, sample truth subset, high-confidence BED, and an aligned BAM slice miners can consume immediately.

This provides the concrete working example required to exercise the validator loop end-to-end (task generation → miner run → hap.py scoring → weight update).
