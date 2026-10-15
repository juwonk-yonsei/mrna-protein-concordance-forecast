<!--
B_revision_B14_claude.md — Claude's manuscript corrections after B 보고서 13, applied by the agent in B 지시서 14 §3.
Same conventions as the earlier revision files.
  [OPTION W/X] choose by B 지시서 14 §2.4: W if the two fits are shown to use different training rows, X if the cause is still not identified.
All tags and comments must be gone from the submitted files, except the [AUTHOR …] items, which the researcher fills.
-->

# D1. Supplementary Note 2 — the reproduction paragraph (factual error, must fix)

The present text says "The LightGBM version used for the earlier fit was not recorded". Report 13 shows it was recorded: 4.7.0 in `reports/B_보고서_2.md` (1 October 2026, 21:57 +0900, before the fit) and in `reports/B_보고서_3.md` (2 October 2026, 00:40 +0900, nine minutes after it), with `lightgbm-4.7.0.dist-info` installed on 1 October at 19:53 and unchanged. The sentence must go.

**Find** the sentence beginning "Refits with 1, 16 and 72 threads" and everything to the end of that paragraph, **and replace with:**

Refits with 1, 16 and 72 threads and with automatic, row-wise and column-wise histogram construction all produced the same predictions, and every setting was deterministic across repeats, so neither thread count nor histogram construction explains the difference. The script that wrote the frozen predictions is byte-identical to the released one, the training rows come from the same file, and LightGBM 4.7.0 was installed on 1 October 2026 and recorded both before and after the two fits. [OPTION W: The two fits differ in the rows they were given: the model frozen for RPPA and ProCan was fitted on [FILL: stage-3 training set in one clause, with its row count], whereas the APOLLO fit and the released model used [FILL: one clause, with its row count]. Refitting on the former reproduced the frozen RPPA and ProCan predictions (maximum absolute differences [FILL] and [FILL]), and that model is released alongside the other.] [OPTION X: We could not identify the cause.] The preregistered statistics recomputed with the released model differed from the frozen values by −2.9 × 10^−5 on RPPA and −3.3 × 10^−4 on ProCan. All results in the paper use the frozen predictions, which are released.

<!-- src: results/stage12/repro_grid.tsv; reports/B_보고서_13.md §3; results/stage14/trainset_compare.tsv -->

---

# D2. Methods, "Forecast models and baselines" — the two-fit paragraph

[OPTION X] Leave the paragraph as it stands.

[OPTION W] **Replace its last two sentences with:**

The two fits were given different training rows ([FILL: the two row counts]); each reproduces its own frozen predictions exactly, and both models are released. All reported results use the frozen predictions.

---

# D3. Discussion, limitations paragraph (number must match the released table)

Report 13 §9 shows that 7,449 genes in the released table have stored CPTAC correlations: the 7,432 of the P set plus 17 ADC targets outside it. The present sentence says 7,432.

**Find:** "covered only the 7,432 genes whose CPTAC correlations were stored at development, so for most genes in the released table the forecast rests on query-side RNA alone;"

**Replace with:**

covered only the 7,449 genes with stored CPTAC correlations, the 7,432 of the P set and 17 further antibody–drug conjugate targets, so for most of the 19,211 genes in the released table the forecast rests on query-side RNA alone;

<!-- src: results/stage12/forecast_v3_check.tsv -->
