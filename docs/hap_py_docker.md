# hap.py Docker Setup (mirroring `nist/` workflow)

This guide captures exactly how the legacy `nist/` automation invoked hap.py so we can remove that directory without losing institutional knowledge. The goal is to stand up the same Dockerized tooling (hap.py + rtg-tools) that powered the buccal-swab validation pipeline.

---

## 1. Why `mgibio/hap.py:v0.3.12`?

The NIST scripts (`run_snake.sh`, `run_snake_sr.sh`) pin the Docker image to `mgibio/hap.py:v0.3.12`. That image bundles:

- hap.py 0.3.12 binaries under `/opt/hap.py/bin/hap.py`
- RTG Tools (vcfeval) templates
- Python + bcftools/tabix dependencies

It uses the newer Docker manifest format, so Docker ≥28/containerd ≥2.1 can pull it (unlike the deprecated `pkrusche/hap.py` images).

---

## 2. Pull the Image

```bash
docker pull mgibio/hap.py:v0.3.12
docker run --rm mgibio/hap.py:v0.3.12 --version    # sanity check
```

If you run Watchtower, add `mgibio/hap.py` to your allowlist so the image stays patched.

For automated validator runs use `scripts/run_hap_validation.py`, which wraps the same Docker invocation and produces the combined TSV outputs (`combined_results_longread*.tsv`). Example:

```bash
python scripts/run_hap_validation.py \
  /abs/path/sample.vcf.gz \
  /abs/path/truth.vcf.gz \
  /abs/path/reference.fa \
  /abs/path/reference.sdf/ \
  /abs/path/confident_regions.bed.gz \
  --bind /abs/path:/abs/path \
  --output-dir /tmp/happy_out
```

--- 

## 3. Mounting Data (match `run_snake*.sh`)

In `nist/run_snake.sh` and `nist/run_snake_sr.sh`, the host paths were:

- `/Users/mahangari/Desktop/work/nist/nist-analysis-buccal-swab-capaccred-102025/…` for reference FASTA + SDF
- `/Users/mahangari/Desktop/work/nist/GRCh38_notinalldifficultregions.bed.gz` for masks
- Sample / truth VCFs inside the same tree

Docker is invoked with:

```bash
docker_binds="/Users/mahangari:/Users/mahangari"
docker_platform=linux/amd64
container_hap_py_bin=/opt/hap.py/bin/hap.py
```

To reproduce this without Snakemake, mount the same directory and run:

```bash
HOST_ROOT=/Users/mahangari/Desktop/work/nist
docker run --rm \
  --platform linux/amd64 \
  -v ${HOST_ROOT}:${HOST_ROOT} \
  mgibio/hap.py:v0.3.12 \
  /opt/hap.py/bin/hap.py \
    ${HOST_ROOT}/nist-analysis-buccal-swab-capaccred-102025/sample.vcf.gz \
    ${HOST_ROOT}/nist-analysis-buccal-swab-capaccred-102025/truth.vcf.gz \
    -r ${HOST_ROOT}/nist-analysis-buccal-swab-capaccred-102025/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
    --engine=vcfeval \
    --engine-vcfeval-template ${HOST_ROOT}/nist-analysis-buccal-swab-capaccred-102025/GRCh38_no_alt_analysis_set_GCA_000001405.15.sdf/ \
    --target-regions ${HOST_ROOT}/GRCh38_notinalldifficultregions.bed.gz \
    -o /tmp/happy_out
```

Replace the sample/truth paths with whatever inputs you are validating.

---

## 4. Snakemake Automation Reference

The `run_snake.sh` wrapper did nothing more than:

```bash
snakemake --snakefile nist-analysis-buccal-swab-capaccred-102025/stats/Snakefile_long_read \
  --config \
    sample_vcf=/abs/path/to/sample.vcf.gz \
    truth_vcf=/abs/path/to/truth.vcf.gz \
    reference=/abs/path/to/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
    reference_sdf=/abs/path/to/GRCh38_no_alt_analysis_set_GCA_000001405.15.sdf/ \
    mask_regions=/abs/path/to/GRCh38_notinalldifficultregions.bed.gz \
    use_docker=true \
    docker_use_sudo=false \
    docker_binds="/Users/mahangari:/Users/mahangari" \
    docker_image=mgibio/hap.py:v0.3.12 \
    docker_platform=linux/amd64 \
    container_hap_py_bin=/opt/hap.py/bin/hap.py \
  --cores 8
```

The Snakefile simply:
1. Runs hap.py with those paths.
2. Aggregates the resulting CSV summaries into `combined_results_longread*.tsv`.

You can copy the same config pattern for any future Snakemake-based validation, or re-use the logic in our Python-based validator loop by invoking hap.py through Docker with identical binds.

---

## 5. Summary Checklist

1. `docker pull mgibio/hap.py:v0.3.12`
2. Ensure host directories containing FASTA/SDF/VCFs are accessible (and optionally mirror `docker_binds="/Users/mahangari:/Users/mahangari"`).
3. Run hap.py via Docker (either manually or through Snakemake) using `/opt/hap.py/bin/hap.py` inside the container.
4. Optional: Configure Watchtower to auto-update the image.

With this captured in the repo, we can delete `nist/` once its data assets are archived elsewhere—the operational knowledge now lives here.
