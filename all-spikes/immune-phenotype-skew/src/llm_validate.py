"""
LLM-as-judge validation for the lexical 4-subtype compatibility classifier.

Takes a stratified random sample of 100 trials (25 per lexical bucket:
immunogenic_gated, cold_gated, mixed, permissive) and asks Claude Sonnet
to independently bucket each one.  Compares against the lexical labels
and reports per-bucket precision plus the full 4x4 confusion matrix.

Uses prompt caching on the system instruction so all 100 calls share
the cached prefix.

Run only if ANTHROPIC_API_KEY is set in env.  Cost estimate:
~100 calls x ~2K input tokens (eligibility text capped at 6K chars)
-> ~$0.70 with prompt caching.
"""

from __future__ import annotations
import json
import os
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

SPIKE = Path(__file__).resolve().parents[1]
ROOT = SPIKE.parents[1]
COMPAT = SPIKE / "data" / "subtype_compat.json"
TRIALS = ROOT / "data" / "trials.json"
OUT = SPIKE / "data" / "llm_validation.json"

MODEL = "claude-sonnet-4-6"
N_PER_BUCKET = 25
BUCKETS = ("immunogenic_gated", "cold_gated", "mixed", "permissive")
SEED = 17

JUDGE_SYSTEM = """\
You are a molecular oncologist auditing a TNBC clinical trial's eligibility \
criteria for the immune-phenotype prior it encodes.  You will be given a \
trial's title, interventions, and eligibility text.  Decide which ONE of the \
following four eligibility-gate categories best describes how the trial \
selects (or fails to select) on immune-phenotype markers:

  immunogenic_gated : eligibility REQUIRES immune-hot markers
                      (PD-L1+, CPS>=1 or >=10, sTIL+, "inflamed tumor")
                      and would therefore exclude immune-cold (BLIS-like)
                      patients
  cold_gated        : eligibility REQUIRES immune-cold-associated markers
                      (gBRCA1/2 mutation carriers, HRD+, PD-L1 negative,
                      "immune-cold" tumors) and would therefore exclude
                      immune-hot (BLIA-like) patients
  permissive        : no immune-axis enrichment criterion -- allows any
                      TNBC patient regardless of PD-L1 / TIL / BRCA status
  mixed             : the trial has multiple arms or strata that gate on
                      BOTH immunogenic and cold markers simultaneously

Return STRICTLY a JSON object on a single line, no prose:
{"category": "immunogenic_gated" | "cold_gated" | "permissive" | "mixed",
 "rationale": "<=30 word string"}
"""


def stratified_sample(compat_records):
    by_bucket = defaultdict(list)
    for r in compat_records:
        b = gate_category(r)
        by_bucket[b].append(r)
    rnd = random.Random(SEED)
    sample = []
    for bucket in BUCKETS:
        pool = by_bucket.get(bucket, [])
        n_take = min(len(pool), N_PER_BUCKET)
        sample.extend(rnd.sample(pool, n_take))
    return sample


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


def user_message(trial):
    text = (trial.get("eligibilityCriteria") or "")[:6000]
    interv = [(iv.get("type"), iv.get("name")) for iv in trial.get("interventions") or []
              if isinstance(iv, dict)]
    return (
        f"Trial: {trial.get('nctId')}\n"
        f"Title: {trial.get('title')}\n"
        f"Interventions: {json.dumps(interv)}\n\n"
        f"Eligibility criteria:\n{text}"
    )


def main():
    key = os.environ.get("ANTHROPIC_API_KEY")
    if not key:
        print("ANTHROPIC_API_KEY not set -- skipping LLM validation. "
              "Set it and re-run to populate llm_validation.json.")
        sys.exit(0)

    import anthropic
    client = anthropic.Anthropic()

    with COMPAT.open() as f:
        compat = json.load(f)
    by_id = {r["nctId"]: r for r in compat}
    with TRIALS.open() as f:
        trials = {t["nctId"]: t for t in json.load(f)}

    sample = stratified_sample(compat)
    print(f"Validating {len(sample)} trials across {len(BUCKETS)} buckets "
          f"({N_PER_BUCKET}/bucket nominal)...")

    results = []
    for i, r in enumerate(sample, 1):
        trial = trials[r["nctId"]]
        try:
            resp = client.messages.create(
                model=MODEL,
                max_tokens=300,
                system=[{"type": "text", "text": JUDGE_SYSTEM,
                         "cache_control": {"type": "ephemeral"}}],
                messages=[{"role": "user", "content": user_message(trial)}],
            )
            raw = resp.content[0].text.strip()
            # extract JSON object
            try:
                start = raw.index("{")
                end = raw.rindex("}") + 1
                judge = json.loads(raw[start:end])
            except Exception:
                judge = {"category": "PARSE_ERROR", "rationale": raw[:200]}
        except Exception as e:
            judge = {"category": "API_ERROR", "rationale": str(e)[:200]}

        lexical = gate_category(r)
        agree = judge.get("category") == lexical
        results.append({
            "nctId": r["nctId"],
            "lexical": lexical,
            "judge": judge.get("category"),
            "agree": agree,
            "rationale": judge.get("rationale"),
        })
        print(f"  [{i:>2}/{len(sample)}] {r['nctId']:<14} lex={lexical:<18} "
              f"judge={judge.get('category','?'):<18} agree={agree}")

    # confusion + agreement
    by_lex = defaultdict(Counter)
    for r in results:
        by_lex[r["lexical"]][r["judge"]] += 1
    overall = sum(1 for r in results if r["agree"]) / len(results) if results else 0.0
    summary = {
        "n_total": len(results),
        "overall_agreement": overall,
        "per_lexical_bucket": {k: dict(v) for k, v in by_lex.items()},
        "results": results,
    }

    with OUT.open("w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nOverall agreement: {overall:.1%}")
    for lex_bucket, conf in by_lex.items():
        agree = conf[lex_bucket]
        total = sum(conf.values())
        print(f"  {lex_bucket:<22} precision: {agree}/{total} = {agree/total:.1%}")
    print(f"Written: {OUT}")


if __name__ == "__main__":
    main()
