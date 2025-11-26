# Genomics Variant-Calling Subnet – Architecture Blueprint

This document captures the concrete architecture we will implement for the genomics variant-calling subnet. It codifies the miner/validator game, data assets, incentive design, anti-cheating strategy, and implementation responsibilities so that subsequent engineering work stays aligned.

---

## 1. Problem Statement

Variant-calling accuracy remains the bottleneck for real-world genomics, yet labs benchmark privately and repeatedly. This subnet converts that fragmented workflow into a public, trusted market:

- **Miners**: run any variant-calling pipeline they like and emit VCFs.
- **Validators**: secretly benchmark those VCFs against private hold-outs and emit accuracy-weighted rewards / certificates.
- **Users**: receive validator-signed genotypes with known accuracy.

The subnet monetizes *trust* rather than raw compute, which makes it an ideal Bittensor workload.

---

## 2. Design Constraints (Jacob Steeves’ Requirements)

1. **Never score miners on public truth sets directly.** Public GIAB remains for human validation/marketing only.
2. **Validator loop must be simple** (≈100 lines) and continuously optimize a clear metric.
3. **Miners must retain wide degrees of freedom** (pipelines, ensembles, tuning), constrained only by the task interface.
4. **Copy resistance relies on incentives + variance**, not secrecy alone.
5. **Hold-out datasets are the “light in the dark.”** Validators guard them; miners never see which sample/region is being scored.

---

## 3. Data Asset Strategy

| Pool | Contents | Visibility | Purpose |
| --- | --- | --- | --- |
| **Public/Tuning** | Raw GIAB truth, 1KG benchmarks, FDA challenge data | Public | Miners tune pipelines; no emissions tied. |
| **Hidden Evaluation** | Derived truth genomes (mutated via BamSurgeon, downsampled BAM/FASTQ, region slices) | Validator-only | Drives emissions. Provides unpredictable tasks with known truth. |
| **User Pool** | Real uploads (FASTQ/BAM) | Miners + validators | Production jobs; scored when possible via consensus or later truth. |

### 3.1 Hidden Evaluation Construction

1. Start from GIAB samples (sample–HG007).
2. Apply BamSurgeon to insert synthetic SNPs/indels/SVs in high-confidence regions to destroy recognizable fingerprints.
3. Downsample to multiple coverages (6×, 15×, 30×) and slice into regions (10 Mb windows or stratified tracks).
4. Generate manifests linking each task ID to its truth VCF/BED and data locator (private to validators).
5. Optional simulated FASTQs (wgsim/ART) extend diversity further.

The validator thus knows the precise merged truth (`truth_GIAB ∪ inserted_variants`) for every derived dataset.

---

## 4. Task Interface (Miner-Facing)

Validators dispatch JSON tasks with the following schema:

```json
{
  "task_id": "UUID",
  "mode": "benchmark" | "user",
  "ref_build": "GRCh38",
  "input_type": "BAM" | "FASTQ",
  "data_locator": "s3://bucket/path" | {"r1": "...", "r2": "..."},
  "regions": ["chr20:10000000-12000000"],
  "output_format": "VCF",
  "constraints": {
    "include_genotypes": true,
    "include_filters": true
  }
}
```

Miners:

- Download the data.
- Run any internal pipeline(s).
- Return a VCF (and optionally upload to a validator-provided bucket).

Degrees of freedom remain unconstrained so miners can innovate.

---

## 5. Validator Loop

The validator process is intentionally simple and repeats indefinitely:

```python
while True:
    tasks = make_tasks()              # sample from hidden pool
    responses = query_miners(tasks)   # call miners via dendrite/HTTP
    scores = score_responses(responses)  # hap.py/vcfeval + truth lookup
    weights = update_weights(scores)  # EMA smoothing + boosts
    submit_weights_to_chain(weights)
```

### 5.1 Task Generation

- Randomly sample `K` entries per round from the hidden evaluation manifest.
- Mix genomes, coverages, and region difficulties to maintain variance.

### 5.2 Miner Querying

- Select a subset of miners for each task (cost control).
- Publish a pool of claimable mutated windows (e.g., 200 per round). Each entry tracks `available → claimed → scored`.
- When a miner requests work, assign a **unique** mutated BAM + reference slice pulled from the same distribution (window size, mutation count, difficulty). Once claimed, no other miner sees that BAM.
- Send the serialized task payload (data locator, regions, reference build), await the VCF (with embedded hotkey ID in the header).

### 5.3 Scoring

1. Recover the private truth `(truth_vcf, region_mask)` via `task_id`.
2. Run `hap.py` or `vcfeval` to compute SNP/INDEL precision & recall.
3. Combine into a weighted F1: `score = 0.7 * F1_snp + 0.3 * F1_indel`.

### 5.4 Weight Updates

- Maintain an exponential moving average (EMA) per miner to smooth noise.
- Apply an *innovation boost* when a miner exceeds the historical best by a statistically significant margin; decay the boost over time.
- Normalize effective scores into final weights for on-chain submission.

---

## 6. Anti-Cheating & Copy Resistance

1. **Hold-out secrecy**: tasks never expose sample identity; derived data is unrecognizable thanks to BamSurgeon + slicing.
2. **Hotkey-linked outputs**: require `##miner_hotkey=<hotkey>` in VCF headers to spot identical outputs.
3. **Variance-based detection**: copying another miner rarely matches their performance across many samples due to inherent variance.
4. **Boost with decay**: innovators get a temporary reward uplift; copycats matching later get only the baseline EMA.
5. **Task diversity**: substantial evaluation pool prevents memorization and ensures incremental improvements remain measurable.

---

## 7. Implementation Plan (Phase 1)

1. **Data preparation scripts**
   - Automate BamSurgeon workflows for each GIAB genome.
   - Downsample and slice BAM/FASTQ assets.
   - Produce a validator-only manifest (JSON/Parquet) linking `task_id → truth metadata`.

2. **Reference miner implementations**
   - Containerized DeepVariant/GATK pipelines for local simulation/testing.

3. **Offline validator harness**
   - Python module implementing task sampling, miner invocation (local call), hap.py scoring, EMA/boost logic, and weight logs.
   - CLI/pytest coverage to validate scoring behavior and anti-cheat heuristics.

4. **Integration hooks**
   - Prepare Docker images and configuration for Bittensor validator/miner templates.
   - Abstract network I/O so the same core logic can run locally or via TAO dendrite.

---

## 8. Deliverables Checklist

- `docs/architecture.md` (this document).
- `docs/data_generation.md` – detailed BamSurgeon + manifest instructions (to be written).
- `scripts/` – reproducible data prep tooling.
- `validator/` – Python package with task sampling, scoring, EMA, weight updates, CLI entrypoint.
- `miners/` – reference miner wrappers.
- `tests/` – automated checks on scoring logic and anti-cheat heuristics.

Maintaining this blueprint as the single source of truth ensures we keep the subnet aligned with the required incentives and can reason rigorously about each subsequent engineering step.
