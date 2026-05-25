"""
Build the headline three-stratum figure for the paper.

Shows the imm:cold gating ratio in three nested strata
(full corpus -> ICB-only -> no-ICB -> no-ICB no-PARPi) as a
side-by-side dodged bar chart, with Wilson 95% CIs.
"""

import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SPIKE = Path(__file__).resolve().parents[1]
STATS = SPIKE / "data" / "headline_stats.json"
OUT = SPIKE / "figures" / "headline_three_strata.png"


def main():
    with STATS.open() as f:
        s = json.load(f)

    strata = [
        ("Full corpus", s["full_corpus"]),
        ("ICB-only", s["icb_stratum"]),
        ("Non-ICB", s["non_icb_stratum"]),
        ("Non-ICB,\nnon-PARPi", s["no_icb_no_parpi_stratum"]),
    ]

    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    x = np.arange(len(strata))
    width = 0.32

    imm_rates = [st["rates"]["immunogenic_gated"] for _, st in strata]
    cold_rates = [st["rates"]["cold_gated"] for _, st in strata]
    imm_ci = [st["ci"]["immunogenic_gated"] for _, st in strata]
    cold_ci = [st["ci"]["cold_gated"] for _, st in strata]

    imm_errl = [max(0, r - c[0]) for r, c in zip(imm_rates, imm_ci)]
    imm_errh = [max(0, c[1] - r) for r, c in zip(imm_rates, imm_ci)]
    cold_errl = [max(0, r - c[0]) for r, c in zip(cold_rates, cold_ci)]
    cold_errh = [max(0, c[1] - r) for r, c in zip(cold_rates, cold_ci)]

    bars1 = ax.bar(x - width / 2, imm_rates, width,
                   yerr=[imm_errl, imm_errh], capsize=4,
                   color="#c44e52", alpha=0.85, label="immunogenic-gated")
    bars2 = ax.bar(x + width / 2, cold_rates, width,
                   yerr=[cold_errl, cold_errh], capsize=4,
                   color="#4c72b0", alpha=0.85, label="cold-gated")

    for xi, (r, st) in zip(x, strata):
        n = st["n"]
        ax.text(xi, -0.012, f"n={n}", ha="center", va="top", fontsize=9, color="#444")

    for b, r in zip(bars1, imm_rates):
        ax.text(b.get_x() + b.get_width() / 2, r + 0.003,
                f"{r:.1%}", ha="center", va="bottom", fontsize=8.5)
    for b, r in zip(bars2, cold_rates):
        ax.text(b.get_x() + b.get_width() / 2, r + 0.003,
                f"{r:.1%}", ha="center", va="bottom", fontsize=8.5)

    # imm:cold ratio annotation above each stratum
    ax.set_xticks(x)
    ax.set_xticklabels([name for name, _ in strata], fontsize=10)
    for xi, (_, st) in zip(x, strata):
        ratio = st["imm_to_cold_ratio"]
        if ratio == float("inf"):
            label = "∞"
        else:
            label = f"{ratio:.2f}×"
        top = max(st["rates"]["immunogenic_gated"], st["rates"]["cold_gated"]) + 0.03
        ax.text(xi, top + 0.018, f"imm:cold\n= {label}", ha="center",
                va="bottom", fontsize=8.5,
                bbox=dict(boxstyle="round,pad=0.25", facecolor="#f4f4f4",
                          edgecolor="#bbb", linewidth=0.5))

    ax.set_ylim(0, 0.21)
    ax.set_ylabel("Fraction of trials in stratum")
    ax.set_title("Immune-phenotype gating in TNBC trial eligibility, across nested strata\n"
                 "(Wilson 95% CIs; n labels below the x-axis)")
    ax.legend(loc="upper right", framealpha=0.95)
    ax.grid(axis="y", alpha=0.3)
    ax.set_axisbelow(True)

    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"Written: {OUT}")


if __name__ == "__main__":
    main()
