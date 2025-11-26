#!/usr/bin/env python3
"""Generate comparative mutation-space figures across windows and variant types."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np

# Assumptions based on GIAB chr20 confident regions and BamSurgeon settings
WINDOW_SIZES_MB = [1, 5, 10]
SNV_COUNTS = [5, 10, 20]
INDEL_COUNTS = [0, 5, 10]
CONFIDENT_FRACTION = 0.8  # fraction of window covered by GIAB confident regions
SNV_ALLELE_CHOICES = 3  # A/C/G/T minus reference
INDEL_ALLELE_CHOICES = 6  # approx. insertion/deletion directions × common lengths


def compute_combinations(window_mb: int, snv_count: int, indel_count: int) -> Dict[str, float]:
    total_variants = snv_count + indel_count
    confident_bases = int(window_mb * 1_000_000 * CONFIDENT_FRACTION)
    if total_variants > confident_bases:
        raise ValueError("Requested variants exceed available confident positions")

    combos_positions = math.comb(confident_bases, total_variants)
    combos_alt = (SNV_ALLELE_CHOICES ** snv_count) * (INDEL_ALLELE_CHOICES ** max(indel_count, 0))
    total_combos = combos_positions * combos_alt
    log10 = math.log10(total_combos)
    return {
        "window_mb": window_mb,
        "snv_count": snv_count,
        "indel_count": indel_count,
        "total_variants": total_variants,
        "combos": total_combos,
        "log10_combos": log10,
        "bits": log10 / math.log10(2),
    }


def main() -> None:
    output_dir = Path(__file__).resolve().parent
    records: List[Dict[str, float]] = []
    for window in WINDOW_SIZES_MB:
        for snv in SNV_COUNTS:
            for indel in INDEL_COUNTS:
                records.append(compute_combinations(window, snv, indel))

    csv_path = output_dir / "mutation_space_scenarios.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "window_mb",
                "snv_count",
                "indel_count",
                "total_variants",
                "combos",
                "log10_combos",
                "bits",
            ],
        )
        writer.writeheader()
        writer.writerows(records)

    # Figure 1: log10 combos for zero-indel case across windows
    plt.figure(figsize=(7, 4))
    for window in WINDOW_SIZES_MB:
        subset = [r for r in records if r["window_mb"] == window and r["indel_count"] == 0]
        subset.sort(key=lambda r: r["snv_count"])
        plt.plot(
            [r["snv_count"] for r in subset],
            [r["log10_combos"] for r in subset],
            marker="o",
            label=f"{window} Mb window",
        )
    plt.xlabel("Synthetic SNVs per window (indels = 0)")
    plt.ylabel("log10(# mutation sets)")
    plt.title("Window size vs combinatorial space")
    plt.grid(True, linestyle=":")
    plt.legend()
    window_fig = output_dir / "mutation_space_windows.png"
    plt.tight_layout()
    plt.savefig(window_fig, dpi=200)

    # Figure 2: heatmap for 5 Mb window with varying SNVs/indels
    base_window = 5
    heat_data = np.zeros((len(SNV_COUNTS), len(INDEL_COUNTS)))
    for i, snv in enumerate(SNV_COUNTS):
        for j, indel in enumerate(INDEL_COUNTS):
            rec = compute_combinations(base_window, snv, indel)
            heat_data[i, j] = rec["log10_combos"]

    plt.figure(figsize=(6, 4))
    im = plt.imshow(heat_data, cmap="viridis", origin="lower")
    plt.colorbar(im, label="log10(# mutation sets)")
    plt.xticks(range(len(INDEL_COUNTS)), INDEL_COUNTS)
    plt.yticks(range(len(SNV_COUNTS)), SNV_COUNTS)
    plt.xlabel("Injected indels per 5 Mb window")
    plt.ylabel("Injected SNVs per 5 Mb window")
    plt.title("Mutation-space heatmap (5 Mb window)")
    for i, snv in enumerate(SNV_COUNTS):
        for j, indel in enumerate(INDEL_COUNTS):
            plt.text(j, i, f"{heat_data[i, j]:.1f}", ha="center", va="center", color="white")
    heatmap_fig = output_dir / "mutation_space_heatmap.png"
    plt.tight_layout()
    plt.savefig(heatmap_fig, dpi=200)
    # Figure 3: entropy scatter vs total variants
    plt.figure(figsize=(7, 4))
    markers = {0: "o", 5: "s", 10: "^"}
    colors = {1: "#1f77b4", 5: "#ff7f0e", 10: "#2ca02c"}
    for rec in records:
        plt.scatter(
            rec["total_variants"],
            rec["bits"],
            color=colors[rec["window_mb"]],
            marker=markers.get(rec["indel_count"], "o"),
            s=60,
            alpha=0.8,
        )

    from matplotlib.lines import Line2D

    color_handles = [
        Line2D([0], [0], marker="o", color="w", label=f"{w} Mb window", markerfacecolor=c, markersize=8)
        for w, c in colors.items()
    ]
    marker_handles = [
        Line2D([0], [0], marker=m, color="w", label=f"{indel} indels", markerfacecolor="gray", markersize=8)
        for indel, m in markers.items()
    ]

    plt.xlabel("Total injected variants (SNVs + indels)")
    plt.ylabel("Bits of entropy (log2 combinations)")
    plt.title("Guessing difficulty by window size and indels")
    plt.grid(True, linestyle=":")
    plt.legend(handles=color_handles + marker_handles, loc="upper left", fontsize=8)
    entropy_fig = output_dir / "mutation_space_entropy.png"
    plt.tight_layout()
    plt.savefig(entropy_fig, dpi=200)

    print(f"Wrote {csv_path}")
    print(f"Wrote {window_fig}")
    print(f"Wrote {heatmap_fig}")
    print(f"Wrote {entropy_fig}")


if __name__ == "__main__":
    main()
