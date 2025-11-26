import csv
import math
from pathlib import Path
import matplotlib.pyplot as plt

WINDOW_SIZE = 10_000_000  # 10 Mb
variants_to_test = [5, 10, 15, 20, 30, 40, 50]

records = []
for k in variants_to_test:
    comb = math.comb(WINDOW_SIZE, k)
    log10 = math.log10(comb)
    records.append({
        "variants": k,
        "combinations": comb,
        "log10_combinations": log10,
    })

output_dir = Path(__file__).resolve().parent
csv_path = output_dir / "mutation_space.csv"
with csv_path.open("w", newline="") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=["variants", "combinations", "log10_combinations"])
    writer.writeheader()
    for row in records:
        writer.writerow(row)

# Plot log10 combinations vs variant count
plt.figure(figsize=(6, 4))
plt.plot([r["variants"] for r in records], [r["log10_combinations"] for r in records], marker="o")
plt.xlabel("Number of synthetic SNVs in 10 Mb window")
plt.ylabel("log10(# of possible mutation sets)")
plt.title("Combinatorial explosion of synthetic SNV sets")
plt.grid(True, linestyle=":")
figure_path = output_dir / "mutation_space.png"
plt.tight_layout()
plt.savefig(figure_path, dpi=200)

print(f"Wrote {csv_path}")
print(f"Wrote {figure_path}")
