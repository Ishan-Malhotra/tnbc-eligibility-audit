"""
Stratified analysis of subtype-compatibility classifications.

Headline claim: trials enrich for BLIA over BLIS, and this skew persists
in the non-ICB subset.

We bucket each trial into one of three exclusive eligibility-gate categories:

  immunogenic_gated : eligibility requires immune-hot markers (PD-L1+, CPS+,
                      TIL+) -- favors BLIA, disadvantages BLIS
  cold_gated        : eligibility requires immune-cold-associated markers
                      (gBRCA-mut, HRD+, PD-L1-negative) -- favors BLIS
  permissive        : no immune-axis gating

When a trial hits both immunogenic and cold gates (rare), we call it
'mixed'.  Trials are not double-counted.

Stratification axes:
  - is_immunotherapy (ICB intervention present)
  - is_parpi (PARP inhibitor intervention present)
  - setting: neoadjuvant / metastatic / refractory / other

Outputs:
  - figures/coverage_full.png         : full-corpus gate-category bar chart
  - figures/coverage_icb_split.png    : same, stratified by is_immunotherapy
  - figures/coverage_parpi_excluded.png : non-ICB AND non-PARPi (the strongest
                                          confound-controlled cut)
  - figures/subtype_enrichment.png    : raw subtype-enrichment counts
  - data/headline_stats.json          : numeric results
"""

from __future__ import annotations

import json
import math
import os
from collections import Counter
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SPIKE = Path(__file__).resolve().parents[1]
COMPAT = SPIKE / "data" / "subtype_compat.json"
STATS_OUT = SPIKE / "data" / "headline_stats.json"
FIG = SPIKE / "figures"


def gate_category(r):
    immuno = r["blia_enriched"]
    cold = r["blis_enriched"]
    if immuno and cold:
        return "mixed"
    if immuno:
        return "immunogenic_gated"
    if cold:
        return "cold_gated"
    return "permissive"


def wilson_ci(k, n, z=1.96):
    """Wilson score interval for a binomial proportion."""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    halfw = (z / denom) * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, centre - halfw), min(1.0, centre + halfw))


def two_proportion_z(k1, n1, k2, n2):
    """Two-proportion z-test (independence of gating rates between two strata).
    Returns (z, two-sided p) via the normal approximation."""
    if n1 == 0 or n2 == 0:
        return float("nan"), float("nan")
    p1, p2 = k1 / n1, k2 / n2
    pooled = (k1 + k2) / (n1 + n2)
    se = math.sqrt(pooled * (1 - pooled) * (1 / n1 + 1 / n2))
    if se == 0:
        return float("nan"), float("nan")
    z = (p1 - p2) / se
    # two-sided p via erfc
    p = math.erfc(abs(z) / math.sqrt(2))
    return z, p


def stratum_stats(records):
    cats = Counter(gate_category(r) for r in records)
    n = len(records)
    out = {
        "n": n,
        "counts": dict(cats),
        "rates": {k: cats[k] / n if n else 0.0 for k in
                  ("immunogenic_gated", "cold_gated", "permissive", "mixed")},
        "ci": {k: wilson_ci(cats[k], n) for k in
               ("immunogenic_gated", "cold_gated", "permissive", "mixed")},
        "imm_to_cold_ratio": (cats["immunogenic_gated"] / cats["cold_gated"])
        if cats["cold_gated"] else float("inf"),
    }
    return out


def bar_with_ci(ax, labels, rates, cis, color, title):
    x = range(len(labels))
    errs_low = [max(0, rate - ci[0]) for rate, ci in zip(rates, cis)]
    errs_high = [max(0, ci[1] - rate) for rate, ci in zip(rates, cis)]
    ax.bar(x, rates, color=color, yerr=[errs_low, errs_high], capsize=4, alpha=0.85)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=15)
    ax.set_ylabel("Fraction of trials")
    ax.set_ylim(0, max(0.85, max(rates) + 0.15 if rates else 0.5))
    ax.set_title(title)
    for xi, r in zip(x, rates):
        ax.text(xi, r + 0.01, f"{r:.1%}", ha="center", va="bottom", fontsize=9)


def make_full_corpus_figure(stats_full):
    cats = ["immunogenic_gated", "cold_gated", "mixed", "permissive"]
    rates = [stats_full["rates"][c] for c in cats]
    cis = [stats_full["ci"][c] for c in cats]
    fig, ax = plt.subplots(figsize=(6, 4.2))
    bar_with_ci(ax, [c.replace("_", "\n") for c in cats], rates, cis,
                color="#5778a4",
                title=f"Eligibility-gate category, full corpus (n={stats_full['n']})")
    fig.tight_layout()
    out = FIG / "coverage_full.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return out


def make_icb_split_figure(stats_icb, stats_nonicb):
    cats = ["immunogenic_gated", "cold_gated", "mixed", "permissive"]
    icb_rates = [stats_icb["rates"][c] for c in cats]
    icb_ci = [stats_icb["ci"][c] for c in cats]
    nonicb_rates = [stats_nonicb["rates"][c] for c in cats]
    nonicb_ci = [stats_nonicb["ci"][c] for c in cats]

    fig, (ax_icb, ax_non) = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    bar_with_ci(ax_icb, [c.replace("_", "\n") for c in cats], icb_rates, icb_ci,
                color="#e49444",
                title=f"ICB trials (n={stats_icb['n']})")
    bar_with_ci(ax_non, [c.replace("_", "\n") for c in cats], nonicb_rates, nonicb_ci,
                color="#5778a4",
                title=f"Non-ICB trials (n={stats_nonicb['n']})")
    ax_icb.set_ylabel("Fraction of trials")
    ax_non.set_ylabel("")
    fig.suptitle("Eligibility-gate category by ICB stratum  (Wilson 95% CI)",
                 fontsize=12)
    fig.tight_layout()
    out = FIG / "coverage_icb_split.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return out


def make_no_icb_no_parpi_figure(stats):
    cats = ["immunogenic_gated", "cold_gated", "mixed", "permissive"]
    rates = [stats["rates"][c] for c in cats]
    cis = [stats["ci"][c] for c in cats]
    fig, ax = plt.subplots(figsize=(6, 4.2))
    bar_with_ci(ax, [c.replace("_", "\n") for c in cats], rates, cis,
                color="#85b6b2",
                title=f"Non-ICB, non-PARPi trials (n={stats['n']})\n"
                      "(strongest confound-controlled stratum)")
    fig.tight_layout()
    out = FIG / "coverage_no_icb_no_parpi.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return out


def make_subtype_counts_figure(records):
    flags = ["blia_enriched", "blia_excluded",
             "blis_enriched", "blis_excluded",
             "lar_enriched", "lar_excluded",
             "mes_enriched", "mes_excluded"]
    counts = {f: sum(1 for r in records if r[f]) for f in flags}
    n = len(records)
    rates = [counts[f] / n for f in flags]
    labels = [f.replace("_", "\n") for f in flags]
    colors = []
    for f in flags:
        if f.endswith("_enriched"):
            colors.append("#5778a4")
        else:
            colors.append("#e49444")
    fig, ax = plt.subplots(figsize=(10, 4.2))
    ax.bar(range(len(flags)), rates, color=colors, alpha=0.85)
    ax.set_xticks(range(len(flags)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Fraction of trials")
    for i, r in enumerate(rates):
        ax.text(i, r + 0.003, f"{r:.1%}", ha="center", va="bottom", fontsize=8)
    ax.set_title(f"Per-subtype enrichment vs. exclusion flags (n={n})")
    # legend
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#5778a4", label="enriched"),
                       Patch(color="#e49444", label="excluded")],
              loc="upper right")
    fig.tight_layout()
    out = FIG / "subtype_enrichment.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return out


def setting_breakdown(records):
    rows = {}
    for setting in ("neoadjuvant", "metastatic", "refractory"):
        sub = [r for r in records if r["settings"].get(setting)]
        rows[setting] = stratum_stats(sub)
    return rows


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    with COMPAT.open() as f:
        records = json.load(f)

    icb = [r for r in records if r["is_immunotherapy"]]
    nonicb = [r for r in records if not r["is_immunotherapy"]]
    no_icb_no_parpi = [r for r in records
                       if not r["is_immunotherapy"] and not r["is_parpi"]]

    stats_full = stratum_stats(records)
    stats_icb = stratum_stats(icb)
    stats_non = stratum_stats(nonicb)
    stats_clean = stratum_stats(no_icb_no_parpi)

    # Headline statistical test: immunogenic_gated rate in non-ICB vs ICB
    z_imm, p_imm = two_proportion_z(
        stats_icb["counts"].get("immunogenic_gated", 0), stats_icb["n"],
        stats_non["counts"].get("immunogenic_gated", 0), stats_non["n"],
    )
    # And: immunogenic_gated rate vs cold_gated rate within non-ICB subset
    # (one-sample test against null of equality)
    k_im_non = stats_non["counts"].get("immunogenic_gated", 0)
    k_co_non = stats_non["counts"].get("cold_gated", 0)
    z_imc_non, p_imc_non = two_proportion_z(
        k_im_non, stats_non["n"], k_co_non, stats_non["n"],
    )
    z_imc_clean, p_imc_clean = two_proportion_z(
        stats_clean["counts"].get("immunogenic_gated", 0), stats_clean["n"],
        stats_clean["counts"].get("cold_gated", 0), stats_clean["n"],
    )

    headline = {
        "full_corpus": stats_full,
        "icb_stratum": stats_icb,
        "non_icb_stratum": stats_non,
        "no_icb_no_parpi_stratum": stats_clean,
        "tests": {
            "immunogenic_gated_rate_icb_vs_nonicb": {"z": z_imm, "p_two_sided": p_imm},
            "imm_vs_cold_within_nonicb": {"z": z_imc_non, "p_two_sided": p_imc_non},
            "imm_vs_cold_within_no_icb_no_parpi": {"z": z_imc_clean,
                                                   "p_two_sided": p_imc_clean},
        },
        "settings_breakdown_full": setting_breakdown(records),
    }

    STATS_OUT.parent.mkdir(parents=True, exist_ok=True)
    with STATS_OUT.open("w") as f:
        json.dump(headline, f, indent=2)

    # Figures
    make_full_corpus_figure(stats_full)
    make_icb_split_figure(stats_icb, stats_non)
    make_no_icb_no_parpi_figure(stats_clean)
    make_subtype_counts_figure(records)

    # Print headline to stdout
    def fmt(stats, label):
        c = stats["counts"]
        n = stats["n"]
        print(f"\n[{label}]  n={n}")
        for k in ("immunogenic_gated", "cold_gated", "mixed", "permissive"):
            ci = stats["ci"][k]
            print(f"  {k:<22} {c.get(k, 0):>5}  ({100 * c.get(k, 0) / n:>5.1f}%, "
                  f"95% CI {ci[0]*100:>5.1f}–{ci[1]*100:>5.1f}%)")
        ratio = stats["imm_to_cold_ratio"]
        print(f"  imm:cold ratio       {ratio:.2f}")

    fmt(stats_full, "Full corpus")
    fmt(stats_icb, "ICB stratum")
    fmt(stats_non, "Non-ICB stratum")
    fmt(stats_clean, "No-ICB, no-PARPi stratum")

    print("\n--- statistical tests ---")
    print(f"immunogenic_gated rate ICB vs non-ICB:        "
          f"z={z_imm:.2f}, p={p_imm:.2e}")
    print(f"imm_gated vs cold_gated within non-ICB:       "
          f"z={z_imc_non:.2f}, p={p_imc_non:.2e}")
    print(f"imm_gated vs cold_gated within no-ICB no-PARPi: "
          f"z={z_imc_clean:.2f}, p={p_imc_clean:.2e}")
    print(f"\nWritten: {STATS_OUT}")
    print(f"Figures in: {FIG}")


if __name__ == "__main__":
    main()
