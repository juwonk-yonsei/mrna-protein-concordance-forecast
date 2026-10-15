<!--
B_revision_B13_claude.md — Claude's manuscript corrections after B 보고서 12, applied by the agent in B 지시서 13 §4.
Same conventions as B_revision_B12_claude.md.
  [OPTION V/N] choose by B 지시서 13 §2.4: V if the LightGBM version explains the RPPA/ProCan difference, N if the cause is still not identified.
All tags and comments must be gone from the submitted files.
-->

# C1. Results, ADC subsection — leftover sentence (correctness, must fix)

**Find:** "For the other 23 targets, the FDA novel-approval pages did not open, so whether an approved ADC exists was not established."

**Replace with:**

A search of US Food and Drug Administration labels found no approved antibody–drug conjugate acting on the other 23 targets. The same search found seven approved conjugates against targets outside this set; forecasts for those targets are in the released genome-wide table.

<!-- src: results/stage12/adc_gene_calls.tsv; results/stage12/adc_approved_v2.tsv (ADCETRIS/TNFRSF8, BESPONSA/CD22, MYLOTARG/CD33, POLIVY/CD79B, ZYNLONTA/CD19, BLENREP/TNFRSF17, DECNUPAZ/IL3RA). [CHECK: that all seven genes are rows of forecast_all_genes_v3; if any is absent, give the number present and name the absent ones.] -->

---

# C2. Methods, "ADC target atlas" — first sentence (overclaim, must fix)

The sentence now says the targets were listed "before any data were analysed". The resource log shows model fitting on 1 October 2026 at 14:55 UTC, before `config/adc_targets.tsv` was committed at 15:10:33 +0900. The claim as written is false.

**Replace the first sentence with:**

The 29 targets were fixed in the project's first task instruction and committed on 1 October 2026, before any analysis of these targets was run; no formal selection rule was recorded. They do not include every target of an approved antibody–drug conjugate.

<!-- [CHECK: the first commit time of any ADC-specific output (results/stage5/atlas/ or results/stage7/adc_*), and confirm it is later than the commit of config/adc_targets.tsv. Report both times. If the ADC output is earlier, tell Claude and do not apply this block.] -->

---

# C3. Supplementary Note 2 — training-row file (wrong hash, must fix)

Block R11 of the previous revision cited SHA-256 `9a00d5dd…`, which is `results/reliability_training_rows.tsv`, a different file. The F1s training rows are derived from `results/f1s_predictions.tsv`.

**Find:** "identical training rows (SHA-256 `9a00d5dd5b304770929711b407b1c9203072411ad33c465be16e69c1559c4976`)"

**Replace with:**

identical training rows, the 72,880 rows of `results/f1s_predictions.tsv` (SHA-256 `3ab5b9378020ef2e3710ba70e6e817416d9a20e055ca8c8086c64b02071334d6`) that are in the P gene set with a finite observed correlation

---

# C4. Supplementary Note 2 — the reproduction paragraph, last sentences

**Find** the sentence that lists the settings tested ("Refits with 1, 16 and 72 threads…") **and replace it with** one of the following.

[OPTION V]
Refits with 1, 16 and 72 threads and with row-wise, column-wise and automatic histogram construction all gave the same result, and each setting was deterministic across repeats. The difference is explained by the LightGBM version: the RPPA and ProCan model was fitted with version [FILL], and refitting with that version reproduced the frozen RPPA and ProCan predictions (maximum absolute differences [FILL] and [FILL]) while failing to reproduce the APOLLO predictions ([FILL]), which version 4.7.0 reproduces exactly. Both versions are recorded in the release.

[OPTION N]
Refits with 1, 16 and 72 threads and with row-wise, column-wise and automatic histogram construction all gave the same result, and each setting was deterministic across repeats, so the difference is not caused by thread count or histogram construction. The LightGBM version used for the earlier fit was not recorded and we could not identify the cause.

<!-- src: results/stage12/repro_grid.tsv; results/stage13/version_refit.tsv -->

---

# C5. Methods, "Forecast models and baselines" — the two-fit paragraph

Keep OPTION N text from the previous revision unless §2.4 gives OPTION V, in which case **replace its last two sentences with:**

The two fits used different LightGBM versions ([FILL] and 4.7.0); each version reproduces its own frozen predictions exactly. All reported results use the frozen predictions, which are released together with both model files.

---

# C6. Data availability — one added sentence

**Append after the sentence listing the forecast tables:**

For the 11,762 of 19,211 genes without stored CPTAC correlations, the gene-level features take their default values, so the forecast rests on the gene's RNA distribution in the query cohort; the released table marks these genes.

<!-- src: results/stage12/forecast_all_genes_v3.tsv.gz; results/stage12/forecast_v3_check.tsv -->

---

# C7. Discussion, limitations paragraph — one added clause

**Find:** "The forecast's gene-level features covered only genes whose CPTAC correlations were stored at development;"

**Replace with:**

The forecast's gene-level features covered only the 7,432 genes whose CPTAC correlations were stored at development, so for most genes in the released table the forecast rests on query-side RNA alone;

---

# C8. Supplementary Note 7 — two added sentences

**Append at the end of the first paragraph:**

The 29 targets were fixed before any analysis of them and are not the set of approved conjugate targets: the label search found approved conjugates against seven further targets (TNFRSF8, CD22, CD33, CD79B, CD19, TNFRSF17 and IL3RA), whose forecasts are in the released genome-wide table.

<!-- [CHECK: that each of the seven is a row of forecast_all_genes_v3; drop any that is not, and say so.] -->
