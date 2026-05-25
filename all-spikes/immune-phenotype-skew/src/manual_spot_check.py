"""
Manual spot-check substitute for the LLM validation slice when no API
key is available.  Prints 2 sampled trials from each gate bucket along
with their classifier flags and eligibility-text excerpt, for human
review.
"""

import json
import random
from pathlib import Path

SPIKE = Path(__file__).resolve().parents[1]
ROOT = SPIKE.parents[1]
COMPAT = SPIKE / "data" / "subtype_compat.json"
TRIALS = ROOT / "data" / "trials.json"

SEED = 17
N = 2


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


def main():
    with COMPAT.open() as f:
        compat = json.load(f)
    with TRIALS.open() as f:
        trials = {t["nctId"]: t for t in json.load(f)}

    rnd = random.Random(SEED)
    for bucket in ("immunogenic_gated", "cold_gated", "permissive", "mixed"):
        pool = [r for r in compat if gate_category(r) == bucket]
        if not pool:
            continue
        sample = rnd.sample(pool, min(N, len(pool)))
        print("=" * 80)
        print(f"BUCKET: {bucket}  (pool size {len(pool)})")
        print("=" * 80)
        for r in sample:
            t = trials[r["nctId"]]
            print(f"\n--- {r['nctId']} ---")
            print(f"Title: {t['title'][:200]}")
            print(f"  is_immunotherapy={r['is_immunotherapy']}  "
                  f"is_parpi={r['is_parpi']}  blia_enriched={r['blia_enriched']}  "
                  f"blis_enriched={r['blis_enriched']}")
            interv = [(iv.get("type"), iv.get("name"))
                      for iv in t.get("interventions") or [] if isinstance(iv, dict)]
            print(f"  Interventions: {interv[:8]}")
            elig = t.get("eligibilityCriteria") or ""
            # show the first eligibility-text region that mentioned any
            # gate term, plus the first 600 chars as fallback
            print(f"  Eligibility excerpt (first 800 chars):")
            print("    " + elig[:800].replace("\n", "\n    "))


if __name__ == "__main__":
    main()
