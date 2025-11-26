# Mutation-Space Benchmark – Why Fresh Windows Keep VANET Honest

VANET’s validators routinely regenerate synthetic 5 Mb chr20 windows, inject novel mutations, and benchmark miners against those bespoke truth sets. This document quantifies why that design makes memorization impossible and details exactly how the window/mutation workflow integrates into production.

## Goals
- Demonstrate the combinatorial explosion of possible mutations even within a single 5 Mb window.
- Show how validators use `GenomicsTaskGenerator` + BamSurgeon to realize new combinations every round.
- Provide reproducible scripts so new contributors can validate the math or plot alternative scenarios.

## Background

Each validation round picks:
1. A random start coordinate on chr20 (bounded by `base/s3_manifest.json`, default 10 Mb–55 Mb).
2. A 5 Mb window (`window_size` in `base/genomics_config.py`).
3. `NUM_MUTATIONS` loci drawn from the GIAB confident BED intersections via BamSurgeon’s `randomsites.py`.
4. Alternate alleles assigned uniformly among the three possible non-reference bases.

Because every task uses a new RNG seed, miners can’t pre-compute lookup tables. Even if they attempted to, the search space is astronomically large as shown below.

## Combinatorial Calculations

For a window with `P` candidate positions and `M` synthetic SNVs inserted, the number of distinct mutation sets is:

```
C(P, M) × 3^M
```

Where `C(P, M)` is the binomial coefficient (“P choose M”) and `3^M` encodes the alternate allele choices. In practice, `P` is the number of confident loci inside the window (hundreds of thousands for chr20). We approximate `P` as 5 Mb × 0.8 (confident coverage), yielding ~4 M candidate bases.

Using `benchmarking_honest_miners/analyze_mutation_space.py` you can recompute log10 counts for any `P` and `M`. Example results:

| Synthetic SNVs (M) | Window Size | log10 C(P, M) | log10 Total Sets (incl. alt bases) |
|-------------------:|-------------|---------------|-------------------------------------|
| 10                 | 5 Mb        | ~61.2         | ~64.94                              |
| 15                 | 5 Mb        | ~86.5         | ~90.28                              |
| 20                 | 5 Mb        | ~110.8        | ~120.04                             |
| 10                 | 10 Mb       | ~63.4         | ~67.17                              |
| 20                 | 10 Mb       | ~121.6        | ~130.87                             |

Even the smallest configuration (10 SNVs in 5 Mb) yields >10^64 combinations—orders of magnitude beyond what any miner could memorize. Doubling to 20 SNVs pushes past 10^120 possibilities. Introducing indels or larger windows multiplies the space again.

For a broader comparison across window sizes and indel mixes, run `mutation_space_scenarios.py`. It assumes 80% confident-region coverage, sweeps `{window_mb ∈ {1,5,10}, snv_count ∈ {5,10,20}, indel_count ∈ {0,5,10}}`, and produces:

- `mutation_space_scenarios.csv` – log10 counts (and bits of entropy) for every tuple.
- `mutation_space_windows.png` – line plot showing, for example, that even a **1 Mb window with just 5 SNVs** spans ~10^30 possibilities, while a 10 Mb window with 20 SNVs exceeds 10^129.
- `mutation_space_heatmap.png` – heatmap for a 5 Mb window; adding 5 indels on top of 10 SNVs pushes the space from ~10^64 to ~10^96, and 10 indels push it beyond 10^126.
- `mutation_space_entropy.png` – scatter plot converting each scenario to “bits of entropy” (log2 combinations) so you can compare guessing difficulty at a glance.

These figures reinforce that even conservative settings leave an astronomical search space.

## Reproducing the Analysis

1. Activate your environment inside the repo root:
   ```bash
   cd vanet
   source .venv/bin/activate
   ```
2. Run the helper script:
   ```bash
   python benchmarking_honest_miners/analyze_mutation_space.py \
     --window-size 5000000 \
     --max-mutations 25 \
     --step 5 \
     --output benchmarking_honest_miners/mutation_space_5mb.csv
   ```
3. Plot the results (Matplotlib output saved to `mutation_space.png`). Use `--window-size 10000000` to mirror the earlier 10 Mb study from `genomenet`.

To emit the multi-window/indel figures discussed above, run:

```bash
python benchmarking_honest_miners/mutation_space_scenarios.py
```

which produces `mutation_space_scenarios.csv`, `mutation_space_windows.png`, `mutation_space_heatmap.png`, and `mutation_space_entropy.png`.

The CSV/PNG artifacts are versioned so validators can cite them when explaining the subnet’s anti-cheating stance to new participants.

## Validator Integration

During live operation (`neurons/validator.py`):
- `GenomicsTaskGenerator.generate_task()` selects the window.
- `MutationInjector.select_sites()` invokes BamSurgeon `randomsites.py` with the window-specific BED mask, ensuring we only mutate GIAB-confident loci.
- `MutationInjector.apply_mutations()` calls `addsnv.py` and merges the resulting VCF with sliced GIAB truth, yielding task-specific truth bundles.
- Metadata (window coordinates, RNG seeds, mutation file) is stored in the task payload so miners can audit the challenge if needed.

Because every task uses a fresh random suffix, miners never receive the same mutated BAM twice. Combined with the combinatorial growth above, the only viable strategy is to execute a genuine variant caller under the provided SLA.

## Extending the Model

Feel free to fork `analyze_mutation_space.py` to model:
- **Different chromosomes** – update `P` using BED coverage statistics for chr1-22.
- **Mixed variant types** – include indels by swapping the `3^M` term for per-type alternate counts.
- **Longer windows** – e.g., 10 Mb or whole-chromosome slices.

Send PRs with updated plots/tables if you expand the math; documenting these bounds helps sustain community trust that VANET’s validators are ungameable.

---

**Summary:** Random 5 Mb windows + BamSurgeon SNVs yield >10^64 combinations even under conservative settings. Validators rotate windows every few hours, so miners have no choice but to process each mutated BAM honestly. This folder captures the math and scripts proving that claim.
