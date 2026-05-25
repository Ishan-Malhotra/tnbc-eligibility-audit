# Immune-Phenotype Skew in TNBC Trial Eligibility (ICB-Controlled)

## Research question

Across the 1,452 TNBC trials registered on ClinicalTrials.gov (as of 2026-05-24), does eligibility-criterion design structurally enrich for the BLIA (basal-like immune-activated) subtype over the BLIS (basal-like immune-suppressed) subtype, **and does that skew persist after excluding immune-checkpoint-blocker (ICB) trials whose enrichment is mandated by FDA labels?**

## Claim under test

Eligibility design encodes an immunogenic-tumor prior that disadvantages BLIS patients, independent of the obvious ICB-labeling confound.

## Operationalization

1. **Subtype-compatibility scoring.** For each trial, derive a 4-element compatibility vector `(LAR, MES, BLIS, BLIA)` from eligibility text + interventions via lexical rules (regex over biomarker tokens) plus an LLM validation slice (~50 trials) to bound classifier error.
2. **ICB stratification.** Tag each trial `is_immunotherapy = True` if any intervention name matches a curated ICB list. Run the BLIA/BLIS enrichment test on the full corpus and on the `is_immunotherapy = False` subset separately.
3. **Foreclosure tests.** Confounders to address explicitly: trial line of treatment (neoadjuvant / adjuvant / metastatic — recovered from eligibility text since the `phase` field is broken in the fetched data); presence of any biomarker gating at all; PARP-inhibitor trials which require BRCA enrichment and have their own selection logic.

## Why this passes the philosophy rubric

- **Surprising:** The naive expectation is "PD-L1 gating tracks ICB drugs." The conditional claim — *that the skew persists in non-ICB trials* — is non-obvious.
- **Fruitful:** Motivates concrete trial-design recommendations (BLIS-enriched arms; design-of-control language audits).
- **Foreclosing alternatives:** The within-non-ICB stratification is the methodological backbone; secondary stratifications (line of treatment, PARPi-only) close the remaining obvious alternatives.
- **Feasible at 2h:** ~30m classifier rules, ~30m LLM slice, ~30m analysis+figures, ~30m TeX draft.

## File layout

- `src/classify_subtypes.py` — lexical 4-subtype compatibility classifier
- `src/llm_validate.py` — Anthropic API call for ~50-trial validation slice
- `src/analyze.py` — stratified statistics, figure generation
- `data/subtype_compat.json` — per-trial compatibility vectors + ICB flag
- `data/llm_validation.json` — LLM-judged labels for the validation slice
- `figures/` — coverage plots
- `paper/main.tex` — paper following `draft-format/caisc_2026.sty`
- `paper/main.pdf` — compiled output

## Progress log (append-only)

- 2026-05-25 T+00:00 — Spike scaffolded. Framing chosen after presenting 3 candidates rated against research_philosophy.md. Budget 2h, Anthropic API used sparingly (validation slice only).
- 2026-05-25 T+00:25 — First-pass lexical classifier shipped. 1,452 trials labeled. Headline ran. Initial result: imm:cold ratio 0.74 full corpus, 1.94 in no-ICB no-PARPi residual (p=0.018).
- 2026-05-25 T+00:35 — Manual spot-check (2/bucket) surfaced two classifier bugs: `name_overlap` not splitting on `+` (under-counted ICB trials by 4); `is_parpi → blis_enriched` over-counted umbrella trials (I-SPY 2, NCT06245889 pembro+olaparib pilot). Both fixed. The PARPi-conflation fix *strengthened* the headline.
- 2026-05-25 T+00:45 — Re-ran classifier + analysis. Final headline: no-ICB no-PARPi imm:cold = 2.00×, z=2.41, p=0.016. ICB stratum 3.12×; non-ICB w/PARPi 0.72×.
- 2026-05-25 T+00:55 — Built the three-stratum headline figure.
- 2026-05-25 T+01:05 — LLM validation harness (`src/llm_validate.py`) written but deferred: no API credits available in env. Manual spot-check (8 trials) substitutes; 8/8 consistent with assigned bucket after the PARPi-conflation fix.
- 2026-05-25 T+01:10 — Reviewer pass (see below).
- 2026-05-25 T+01:15 — Paper draft (`paper/main.tex`) written. 8-page body + appendix + both mandatory CAISc 2026 checklists.
- 2026-05-25 T+01:30 — PDF compilation blocked: no LaTeX on host; `brew install --cask basictex` failed because `/opt/homebrew` required `sudo chown` (user ran the chown in external Terminal).
- 2026-05-25 T+01:45 — User installed BasicTeX externally + `tlmgr install environ lineno units microtype collection-fontsrecommended`. natbib defaulted to author-year and rejected the numeric `thebibliography`; fixed with `\PassOptionsToPackage{numbers,compress}{natbib}` before the style file. Removed duplicate `\section*{References}` that conflicted with the env's own heading.
- 2026-05-25 T+01:55 — **First PDF compiled: `paper/main.pdf` (10 pages, 350 KB).** Headline 2.00× imm:cold finding rendered correctly in abstract.
- 2026-05-25 T+02:30 — LLM validation harness executed against funded API (100-trial stratified sample, Sonnet-class judge with prompt caching, ~$0.70 spend). Headline metrics: overall agreement 55%, Cohen's κ=0.40. Per-bucket precision: permissive 96%, cold-gated 76%, immunogenic-gated 36%, mixed 12%.
- 2026-05-25 T+02:35 — **Residual-stratum-specific re-analysis: immunogenic-bucket precision drops to 25% (3/12); cold-bucket rises to 100% (5/5).** Applying these as multiplicative corrections to the headline 34/17 counts: estimated true imm ≈ 8.5, true cold ≈ 17, LLM-grounded gate-level ratio ≈ **0.50× (inverts the headline direction).** Surfaced the inversion to user before writing A.6.
- 2026-05-25 T+02:40 — User chose honest reframe over "conservative floor" framing. Rewrote abstract, intro contributions, discussion §6.1, conclusion to lead with prose-level template diffusion (real, 2:1, lexical-measurable) as the substantive claim; gate-level inversion documented in §3.5 and Appendix A.6. Lexical 2.00× retained as a measurement but with its referent precisely identified.
- 2026-05-25 T+02:50 — Added `\usepackage{amsmath}` (needed for align*/text in A.6 equation block). Final PDF: **13 pages, 367 KB.** Body pp.1–7 (within 8-page CAISc limit), refs p.7, appendix incl. A.6 pp.8–10, AI Involvement Checklist pp.11–12, Reproducibility Checklist p.13. Spike complete; framing is now defensible at reviewer level.

## Reviewer pass (voila step 14)

I read the draft as a skeptical AI4Science workshop reviewer.

**Strengths:**
- Falsifiable conditional claim with a clean test (nested-stratum design)
- Confounders (ICB labels, PARPi labels) explicitly stripped; not hand-waved
- Wilson CIs + z-test reported; Bonferroni mentioned and survived
- Lexical classifier brittleness flagged honestly; the manual spot-check found a real bug (the PARPi conflation) and the fix *strengthened* the result, which is the right direction
- Code + data + figures released; reproducible end-to-end
- Three competing readings of the residual signal are negotiated rather than picked

**Weaknesses a real reviewer would flag:**
1. **LLM-judge validation deferred.** Single biggest gap. The harness exists, costs ~$0.20 to run, and would either confirm or break the 2× residual. Without it the headline rests on a hand-curated regex bank. **Workshop OK; main-conference requires this run.**
2. **Effect size is modest.** 34 vs 17 phenotype-gated trials. CIs are wide and Bonferroni-corrected p (~0.048) is borderline. Reviewer will (correctly) push on power.
3. **LAR/MES axes are weak.** The four-subtype framing is honest but the BLIA/BLIS axis is doing all the analytical work. Could be reframed as a 2-axis study without loss; the broader framing was a deliberate choice for fruitfulness over narrowness.
4. **No temporal slice.** Single 2026-05-24 snapshot. Whether the immunogenic prior is increasing or decreasing over time is the obvious follow-up.
5. **"Implications for trial design" are not operationalized.** A concrete checklist of "audit these PD-L1 gates" patterns would close this.

**Verdict:** Workshop-acceptable as-is. Main-conference bar requires the LLM-judge validation pass + at minimum a temporal-slice ablation.

