"""
Lexical 4-subtype compatibility classifier for TNBC trial eligibility text.

Burstein taxonomy (LAR / MES / BLIS / BLIA).  For each trial we produce two
booleans per subtype:

  <subtype>_enriched : trial actively selects FOR this subtype's phenotype
  <subtype>_excluded : trial actively selects AGAINST this subtype's phenotype

These are not mutually exclusive; a trial can be enriched for BLIA AND
excluded for BLIS (e.g. PD-L1 CPS>=10 enrichment).  Most "any TNBC" trials
will be neither enriched nor excluded for most subtypes -- that's the
default-permissive case.

We also flag is_immunotherapy (ICB intervention present) and is_parpi
(PARP inhibitor intervention present) for downstream stratification.

Output: data/subtype_compat.json -- list of records keyed by nctId.
"""

from __future__ import annotations
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
TRIALS = ROOT / "data" / "trials.json"
OUT = Path(__file__).resolve().parents[1] / "data" / "subtype_compat.json"


# ---------------------------------------------------------------------------
# Drug lexicons
# ---------------------------------------------------------------------------
ICB_DRUGS = {
    # PD-1
    "pembrolizumab", "keytruda", "mk-3475",
    "nivolumab", "opdivo",
    "cemiplimab", "libtayo",
    "dostarlimab", "jemperli",
    "toripalimab", "tislelizumab", "camrelizumab", "sintilimab",
    "retifanlimab", "spartalizumab", "balstilimab", "zimberelimab",
    "serplulimab", "penpulimab",
    # PD-L1
    "atezolizumab", "tecentriq",
    "durvalumab", "imfinzi",
    "avelumab", "bavencio",
    # CTLA-4
    "ipilimumab", "yervoy",
    "tremelimumab",
    # LAG-3 / TIM-3 etc. (broader checkpoint family)
    "relatlimab", "lag525", "imp701",
    "tsr-022",
}

PARPI_DRUGS = {
    "olaparib", "lynparza",
    "talazoparib", "talzenna",
    "niraparib", "zejula",
    "rucaparib", "rubraca",
    "veliparib",
    "pamiparib",
    "fluzoparib",
}

AR_ANTAGONIST_DRUGS = {
    "bicalutamide",
    "enzalutamide",
    "darolutamide",
    "abiraterone",
    "apalutamide",
    "seviteronel",
    "orteronel",
}

MES_PATHWAY_DRUGS = {
    # mTOR / Notch / EMT-axis drugs whose presence signals mesenchymal-biology targeting
    "everolimus", "afinitor",
    "sirolimus", "rapamycin",
    "temsirolimus",
    "cb-103",
    "ro4929097",
    "mk-0752",
    "pf-03084014", "nirogacestat",
    "ipatasertib",
    "capivasertib",
    "buparlisib", "bkm120",
    "alpelisib",
}


# ---------------------------------------------------------------------------
# Eligibility text patterns
# ---------------------------------------------------------------------------
# BLIA enrichment = trial requires immunogenic-tumor markers
BLIA_ENRICH_PATTERNS = [
    re.compile(r"\bPD[-\s]?L1\b[^.]{0,80}\b(?:positive|expression|express|\+|positivity|≥|>=|of\s+at\s+least|combined\s+positive|cps)\b", re.I),
    re.compile(r"\bcps\b[^.]{0,40}\b(?:≥|>=|of\s+at\s+least|greater\s+than|≥\s*1|≥\s*10|>=\s*1|>=\s*10)", re.I),
    re.compile(r"\bcombined\s+positive\s+score\b", re.I),
    re.compile(r"\btil(?:s)?\b[^.]{0,60}\b(?:≥|>=|positive|greater\s+than|>\s*\d+\s*%)", re.I),
    re.compile(r"\btumou?r[-\s]infiltrating\s+lymphocytes?\b[^.]{0,80}(?:≥|>=|positive|greater\s+than|enriched)", re.I),
    re.compile(r"\bstil\b[^.]{0,40}(?:≥|>=|positive|greater\s+than)", re.I),
    re.compile(r"\bimmune[-\s]activated\b", re.I),
    re.compile(r"\b(?:hot|inflamed)\s+tumou?rs?\b", re.I),
]

# BLIS exclusion = trial requires PD-L1+/CPS+ enrollment, excluding immune-cold tumors
BLIS_EXCLUDE_PATTERNS = BLIA_ENRICH_PATTERNS  # same gate excludes immune-cold

# BLIA exclusion = trial requires PD-L1-negative (rare) -- explicit anti-immune-hot gate
BLIA_EXCLUDE_PATTERNS = [
    re.compile(r"\bPD[-\s]?L1\b[^.]{0,60}\b(?:negative|≤\s*0|<\s*1\s*%|absent|no\s+expression)\b", re.I),
    re.compile(r"\bimmune[-\s]suppressed\b", re.I),
]

# BLIS enrichment = trial actively enriches for immune-cold / BRCA-mut / p53-mut / etc.
# Rare in practice; mostly captured by PARPi trials enriching for gBRCA-mut.
BLIS_ENRICH_PATTERNS = [
    re.compile(r"\b(?:germline\s+)?BRCA\s*(?:1|2|1/2)?\b[^.]{0,40}\b(?:positive|mutated|mutation|carrier|alteration|deleterious)\b", re.I),
    re.compile(r"\bhomologous\s+recombination\s+deficien", re.I),
    re.compile(r"\bHRD\b\s*(?:positive|\+|high)", re.I),
    re.compile(r"\bPD[-\s]?L1\b[^.]{0,60}\b(?:negative|low|absent|<\s*1\s*%)", re.I),
    re.compile(r"\bimmune[-\s](?:cold|desert)\b", re.I),
    re.compile(r"\bBRCA\s*(?:1|2|1/2)?[-\s]?(?:mutant|mut|deficient)\b", re.I),
]

# LAR enrichment
LAR_ENRICH_PATTERNS = [
    re.compile(r"\bandrogen\s+receptor\b[^.]{0,40}\b(?:positive|expression|≥|>=|\+)", re.I),
    re.compile(r"\bAR[-\s]?(?:positive|\+)\b", re.I),
    re.compile(r"\bAR\s+ihc\b[^.]{0,30}\b(?:≥|>=|positive|greater)", re.I),
    re.compile(r"\bluminal\s+androgen\s+receptor\b", re.I),
    re.compile(r"\bLAR\b\s+(?:subtype|tnbc|positive)", re.I),
]

# LAR exclusion (rare)
LAR_EXCLUDE_PATTERNS = [
    re.compile(r"\bAR[-\s]?negative\b", re.I),
    re.compile(r"\bandrogen\s+receptor\b[^.]{0,40}\bnegative\b", re.I),
]

# MES enrichment - hard to detect lexically; use intervention class as proxy + EMT terms
MES_ENRICH_PATTERNS = [
    re.compile(r"\bmesenchymal\b[^.]{0,40}\b(?:subtype|like|tnbc)", re.I),
    re.compile(r"\bepithelial[-\s]?mesenchymal\s+transition\b", re.I),
    re.compile(r"\bEMT\b\s+(?:signature|positive|high)", re.I),
]

# MES exclusion (very rare)
MES_EXCLUDE_PATTERNS = []


# ---------------------------------------------------------------------------
# Line-of-treatment / setting tags (since phase field is broken)
# ---------------------------------------------------------------------------
SETTING_PATTERNS = {
    "neoadjuvant": re.compile(r"\bneo[-\s]?adjuvant\b", re.I),
    "adjuvant": re.compile(r"\badjuvant\b", re.I),
    "metastatic": re.compile(r"\b(?:metastatic|stage\s*IV|advanced)\b", re.I),
    "refractory": re.compile(r"\b(?:refractory|previously\s+treated|prior\s+(?:therapy|treatment|line)|relapsed)\b", re.I),
}


# ---------------------------------------------------------------------------
# Per-trial classifier
# ---------------------------------------------------------------------------
def intervention_names(trial):
    names = set()
    for iv in trial.get("interventions") or []:
        if not isinstance(iv, dict):
            continue
        n = (iv.get("name") or "").lower().strip()
        if n:
            names.add(n)
        for alt in iv.get("otherNames") or []:
            a = (alt or "").lower().strip()
            if a:
                names.add(a)
    return names


def any_match(patterns, text):
    return any(p.search(text) for p in patterns)


def name_overlap(names, lexicon):
    """Token-level membership: a drug name is a 'hit' if any of its
    whitespace-separated tokens (or the full string) appears in the lexicon."""
    for n in names:
        if n in lexicon:
            return True
        # also check token-by-token to catch e.g. "nab-paclitaxel" vs "paclitaxel"
        # and combo names like "ATRA+Toripalimab+chemo"
        toks = re.split(r"[\s\-/,()+:]+", n)
        for t in toks:
            if t and t in lexicon:
                return True
    return False


def classify_trial(trial):
    text = trial.get("eligibilityCriteria") or ""
    names = intervention_names(trial)

    is_immunotherapy = name_overlap(names, ICB_DRUGS)
    is_parpi = name_overlap(names, PARPI_DRUGS)
    has_ar_antagonist = name_overlap(names, AR_ANTAGONIST_DRUGS)
    has_mes_pathway = name_overlap(names, MES_PATHWAY_DRUGS)

    blia_enriched = any_match(BLIA_ENRICH_PATTERNS, text)
    blia_excluded = any_match(BLIA_EXCLUDE_PATTERNS, text)
    # NB: is_parpi is NOT used as a BLIS enrichment signal -- spot-checking
    # showed many umbrella/basket trials have a PARPi arm without actually
    # enriching for BRCA-mut patients (e.g. I-SPY 2).  We rely on
    # eligibility-text BRCA/HRD language only.
    blis_enriched = any_match(BLIS_ENRICH_PATTERNS, text)
    blis_excluded = any_match(BLIS_EXCLUDE_PATTERNS, text)
    # Note: BLIS_EXCLUDE_PATTERNS == BLIA_ENRICH_PATTERNS by design (PD-L1+
    # enrichment simultaneously selects in BLIA and selects out BLIS).

    lar_enriched = any_match(LAR_ENRICH_PATTERNS, text) or has_ar_antagonist
    lar_excluded = any_match(LAR_EXCLUDE_PATTERNS, text)

    mes_enriched = any_match(MES_ENRICH_PATTERNS, text) or has_mes_pathway
    mes_excluded = any_match(MES_EXCLUDE_PATTERNS, text)

    settings = {k: bool(p.search(text)) for k, p in SETTING_PATTERNS.items()}

    return {
        "nctId": trial.get("nctId"),
        "is_immunotherapy": is_immunotherapy,
        "is_parpi": is_parpi,
        "has_ar_antagonist": has_ar_antagonist,
        "has_mes_pathway": has_mes_pathway,
        "blia_enriched": blia_enriched,
        "blia_excluded": blia_excluded,
        "blis_enriched": blis_enriched,
        "blis_excluded": blis_excluded,
        "lar_enriched": lar_enriched,
        "lar_excluded": lar_excluded,
        "mes_enriched": mes_enriched,
        "mes_excluded": mes_excluded,
        "settings": settings,
    }


def main():
    with TRIALS.open() as f:
        trials = json.load(f)

    out = [classify_trial(t) for t in trials]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as f:
        json.dump(out, f, indent=2)

    # quick headline tally
    n = len(out)
    def pct(c):
        return f"{c} ({100 * c / n:.1f}%)"

    icb = sum(1 for r in out if r["is_immunotherapy"])
    parpi = sum(1 for r in out if r["is_parpi"])
    blia_e = sum(1 for r in out if r["blia_enriched"])
    blis_e = sum(1 for r in out if r["blis_enriched"])
    blis_x = sum(1 for r in out if r["blis_excluded"])
    lar_e = sum(1 for r in out if r["lar_enriched"])
    mes_e = sum(1 for r in out if r["mes_enriched"])

    print(f"Trials classified: {n}")
    print(f"  is_immunotherapy:     {pct(icb)}")
    print(f"  is_parpi:             {pct(parpi)}")
    print(f"  BLIA-enriched:        {pct(blia_e)}")
    print(f"  BLIS-enriched:        {pct(blis_e)}")
    print(f"  BLIS-excluded:        {pct(blis_x)}")
    print(f"  LAR-enriched:         {pct(lar_e)}")
    print(f"  MES-enriched:         {pct(mes_e)}")
    print(f"\nWritten: {OUT}")


if __name__ == "__main__":
    main()
