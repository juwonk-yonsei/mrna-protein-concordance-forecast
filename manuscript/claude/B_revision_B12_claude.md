<!--
B_revision_B12_claude.md — Claude's manuscript revisions after B 보고서 10–11, applied by the agent in B 지시서 12 §5.
How to use (agent):
  - Each block R# names the file, the text to find, and the replacement. Replace whole paragraphs or sentences exactly as stated.
  - The working files are manuscript/srep/B_main_srep.md, manuscript/srep/B_SI_srep.md and the cover letter source.
  - B 보고서 10 already corrected some sentences (ProCan ADC 0.829 with 6 pairs, interval half-width sentence, haematopoietic lineage sentence, Supplementary Note 2 value 2.4 × 10^−16). Keep those corrections; the blocks below do not touch them unless stated.
  - [FILL: ...]  value from the named file. Never guess.
  - [CHECK: ...] verify against code, config or result files; fix the text if they disagree and report the change.
  - [OPTION R/N] choose by B 지시서 12 §2.4: R if a setting reproduced the frozen RPPA and ProCan predictions, N if not.
  - [OPTION F/U] choose by B 지시서 12 §4: F if the openFDA search ran, U if it could not be run.
  - <!-- src --> comments give the source for the ledger.
  - All tags and comments must be gone from the submitted files.
-->

# R1. Abstract (B_main_srep.md)

**Find:** "In post hoc comparisons, it was at least as accurate as looking up observed CPTAC correlations, even from the matched cancer type."

**Replace with:**

Post hoc, it was at least as accurate as looking up observed CPTAC correlations, pooled or cancer-matched, for genes whose CPTAC data it had been given.

**Find:** "Adjusting prediction intervals with 50 locally measured genes brought 80% interval coverage to 0.78–0.80 on all three platforms."

**Replace with:**

Recalibrating intervals with 50 locally measured genes brought 80% coverage to 0.78–0.80 on all three platforms.

<!-- src: results/stage11/l1_split_diff.tsv (layer a); results/stage10/lpool_diff.tsv; results/stage7/calibration_table.tsv. Recount: abstract ≤ 200 words. -->

---

# R2. Results, "Comparison with looking up observed correlations" (B_main_srep.md)

Keep the first two paragraphs of this subsection (the definition of the cascade lookup and the paragraph beginning "The forecast was more accurate than the cascade lookup on all three platforms").

**Find and delete:** the paragraph inserted in B10 that begins "Averaging a gene's correlations over all CPTAC cancer types was itself more accurate…", and the paragraph that begins "For genes outside the P set, for which CPTAC information was sparse or absent…".

**Insert in their place these two paragraphs:**

Because a matched-cancer correlation is estimated in a single cohort, we also compared the forecast with the mean of the gene's CPTAC correlations across cancer types, the pooled lookup, on the pairs that had a matched-cancer correlation (L1 pairs). The pooled lookup was more accurate than the matched-cancer correlation on APOLLO (0.028, 0.019–0.037) and not clearly different from it on RPPA (−0.032, −0.085 to 0.012) and ProCan (0.005, −0.017 to 0.026). The forecast was more accurate than the pooled lookup on RPPA (0.058, 0.019–0.105) and ProCan (0.021, 0.009–0.032) and less accurate on APOLLO (−0.010, −0.016 to −0.003) (Table 3).
<!-- src: results/stage10/lpool_diff.tsv -->

The APOLLO deficit came from 600 of the 6,982 L1 pairs, all for genes outside the P set. No CPTAC correlation had been stored for these genes at development, so the forecast's gene-level features took their default values, whereas the lookups used correlations recomputed later for all genes quantified in each CPTAC cohort (Methods). On these pairs, the forecast was less accurate than the pooled lookup (−0.122, −0.188 to −0.055). On the other 6,382 pairs, the forecast and the pooled lookup were not clearly different (0.001, −0.004 to 0.006), and both were more accurate than the matched-cancer correlation (forecast 0.031, 0.021–0.040; pooled lookup 0.030, 0.021–0.039). On RPPA and ProCan, all L1 pairs were of the second kind. This split was defined after the pooled-lookup results had been seen. Against the CPTAC mean lookup, which falls back to a global mean for genes without stored correlations, the forecast gained most outside the P set: 0.139 (0.070–0.209) on ProCan and 0.276 (0.220–0.328) on APOLLO, against 0.002 (−0.003 to 0.007) inside the P set on APOLLO (Supplementary Table S2).
<!-- src: results/stage11/l1_split_diff.tsv; results/stage11/l1_split_counts.tsv; results/stage7/f0_diff.tsv -->

---

# R3. Results, ADC sentence (B_main_srep.md)

**Find:** the sentence "Eight approved ADCs act on six of the 29 targets (Supplementary Table S7)."

**Replace with** [OPTION F]:

[FILL: number] approved ADCs act on [FILL: number] of the 29 targets (Supplementary Table S7).

[OPTION U]: keep the sentence and add: "Approval status could not be verified for [FILL: number] targets."

<!-- src: results/stage12/adc_approved_v2.tsv -->

---

# R4. Discussion (B_main_srep.md)

**Replace the first paragraph** (beginning "The central result is…") **with:**

The central result is that the concordance between a gene's mRNA and its protein can be forecast in cohorts where protein was not measured, and that the forecast transferred to three platforms that differ from the training data in technology and sample type. The preregistered tests compared the forecast with an RNA-variance baseline. The post hoc comparisons added a more demanding benchmark, looking up observed correlations. Where the forecast had been given a gene's CPTAC data, it was at least as accurate as these lookups on all three platforms; where a lookup drew on CPTAC correlations that the forecast's features had not used, the lookup was better.

**Replace the second paragraph** (beginning "Why can a forecast compete with a measurement?") **with:**

Why can a forecast compete with a measurement? A matched-cancer correlation is a measurement in a different cohort. CPTAC cohorts have 82 to 121 tumors per cancer type, so each correlation carries sampling error, and cohorts differ in sample handling and assays. Averaging a gene's correlations over cancer types reduces this error, and on APOLLO the pooled lookup was more accurate than the matched correlation. The forecast adds the gene's RNA distribution in the query cohort. This matters because a correlation is limited by how much its variables vary in the sample: a gene whose RNA varies little across the tumors of a cohort cannot have high concordance there, whatever its concordance elsewhere. On APOLLO, a TMT-labelled tumor cohort like the CPTAC cohorts, the forecast and the pooled lookup were equally accurate for the genes the forecast had CPTAC data on. On RPPA and on cell-line proteomes, the forecast was more accurate than the pooled lookup, which suggests that the query-side information helps most when the platform or sample type differs from the training data.
<!-- src: results/cptac_sample_counts.tsv (OV 82, BRCA 121); results/stage10/lpool_diff.tsv; results/stage11/l1_split_diff.tsv -->

**Replace the third paragraph** (beginning "These results suggest rules of practice.") **with:**

These results suggest rules of practice. They rest on post hoc analyses and should be read with that status. First, where CPTAC has measured a gene, its correlations pooled across cancer types are a better reference than the matched cancer type alone. If the forecast was built with the gene's CPTAC data, it was at least as accurate as the pooled value; if not, the pooled value should be preferred. The released forecast table gives the pooled correlation and marks the genes whose CPTAC data the forecast used. Second, for genes without CPTAC correlations, the forecast with its interval is the available estimate, and it gained most over a lookup that falls back to a global average. Third, on a new platform, intervals should be recalibrated with at least 50 genes measured on that platform; uncalibrated intervals undercovered on RPPA and ProCan. Antibody-based assays deserve particular caution: on RPPA lung adenocarcinoma, observed concordance was on average 0.184 lower than forecast, in line with the known variability of antibody quality[@Upadhya2023].
<!-- src: results/stage6/posthoc/posthoc_02_bias.tsv -->

**In the limitations paragraph, replace the first sentence** ("The comparisons with lookups and the calibration recommendation are post hoc.") **with:**

The comparisons with lookups, the split of the APOLLO lookup comparison by the forecast's CPTAC coverage and the calibration recommendation are post hoc. The forecast's gene-level features covered only genes whose CPTAC correlations were stored at development; retraining with all CPTAC-quantified genes may close the gap seen for the other genes, but we did not test this.

**Replace the summary paragraph** (beginning "In summary,") **with:**

In summary, a concordance forecast trained on CPTAC tumors ranked mRNA–protein concordance on three external platforms better than an RNA-variance baseline in preregistered tests, and in post hoc comparisons it was at least as accurate as looking up observed CPTAC correlations for genes whose CPTAC data it had been given. With a short local calibration, it gives a practical estimate of how far mRNA can stand in for protein in cohorts without proteomes.

---

# R5. Methods, "Forecast models and baselines" (B_main_srep.md)

**Append at the end of the first paragraph** (after the interval half-width sentence):

[OPTION R]
The forecast model was fitted twice with the same training rows, code and parameters: before the RPPA and ProCan predictions were frozen and before the APOLLO predictions were frozen. The two fits differed in [FILL: the setting identified in B12 §2, e.g. "the number of threads (16 and 64)"], which changes LightGBM's results through [FILL: one clause on the mechanism, e.g. "the order of floating-point summation"]. Both fitted models are released, and each reproduces its frozen predictions exactly (Supplementary Note 2).

[OPTION N]
The forecast model was fitted twice with the same training rows, code and parameters: before the RPPA and ProCan predictions were frozen and before the APOLLO predictions were frozen. The released model reproduces the APOLLO predictions exactly but not the RPPA and ProCan predictions (maximum absolute differences 0.031 and 0.044; Spearman correlations 0.9991 and 0.9992), and the cause was not identified (Supplementary Note 2). All reported results use the frozen predictions, which are released.

<!-- src: results/stage12/repro_grid.tsv; results/stage9/refit_compare.tsv; results/stage11/b4_refit_compare.tsv -->

---

# R6. Methods, "Post hoc analyses" (B_main_srep.md)

**Find** the sentence that defines the pooled lookup (inserted in B10 for "[FILL: definition of the pooled lookup used in B10 §2, one sentence.]") **and replace it with:**

The pooled lookup was the mean of the gene's correlations across CPTAC cancer types in this recomputed set (level L2 of the cascade) [CHECK: whether the mean includes the matched cancer type's own correlation; if it does, add "including the matched cancer type"; if not, add "excluding the matched cancer type"], and it was compared with the matched correlation and with the forecast on pairs that had a matched correlation (L1 pairs). After these comparisons had been seen, the APOLLO L1 pairs were split by whether F0 was gene-specific or equal to the global mean, that is, whether the forecast's gene-level features had been computed from stored CPTAC correlations; this split was committed before calculation (Supplementary Note 3).

<!-- src: config/analysis_plan_B10.yaml (or the B10 plan file); config/analysis_plan_B11.yaml -->

---

# R7. Methods, "ADC target atlas" (B_main_srep.md)

**Replace the first sentence** ("The 29 targets were [CHECK: …].") **with:**

The 29 targets were listed in the first analysis plan of the study, before any data were analysed; no formal selection rule was recorded. [CHECK: the first commit of config/adc_targets.tsv precedes every analysis output; give its date in the report.]

**Replace the last sentence** ("Approved ADCs were identified from US Food and Drug Administration prescribing information as of 2 October 2026 (Supplementary Table S7).") **with** [OPTION F]:

Approved ADCs were identified as of 2 October 2026 by searching US Food and Drug Administration labels through the openFDA interface for antibody–drug conjugates and reading the mechanism-of-action section of each label (Supplementary Note 7, Supplementary Table S7).

[OPTION U]: keep the sentence as it is.

---

# R8. Data availability (B_main_srep.md)

- DepMap 24Q4 must read `https://doi.org/10.25452/figshare.plus.27993248.v1` here and in the reference list (with `.v1`).
- [OPTION R]: change "the trained forecast model" to "the two fitted forecast models".

---

# R9. Table 3 (B_main_srep.md) — replace the whole table and its notes

**Table 3. Post hoc comparisons with lookups of observed CPTAC correlations.**

| Contrast | Pairs | RPPA | ProCan | APOLLO |
|---|---|---|---|---|
| Forecast − cascade lookup | All | 0.069 (0.021 to 0.129) | 0.029 (0.011 to 0.044) | 0.021 (0.011 to 0.032) |
| Forecast − CPTAC mean lookup | All | 0.064 (0.021 to 0.112) | 0.032 (0.019 to 0.044) | 0.028 (0.021 to 0.035) |
| Hybrid rule − forecast | All | −0.008 (−0.024 to 0.009) | −0.009 (−0.022 to 0.004) | −0.016 (−0.026 to −0.007) |
| Forecast − matched lookup | L1 | 0.026 (−0.025 to 0.074) | 0.026 (0.002 to 0.050) | 0.018 (0.008 to 0.028) |
| Pooled lookup − matched lookup | L1 | −0.032 (−0.085 to 0.012) | 0.005 (−0.017 to 0.026) | 0.028 (0.019 to 0.037) |
| Forecast − pooled lookup | L1 | 0.058 (0.019 to 0.105) | 0.021 (0.009 to 0.032) | −0.010 (−0.016 to −0.003) |
| Forecast − matched lookup | L1, covered | same as L1 | same as L1 | 0.031 (0.021 to 0.040) |
| Pooled lookup − matched lookup | L1, covered | same as L1 | same as L1 | 0.030 (0.021 to 0.039) |
| Forecast − pooled lookup | L1, covered | same as L1 | same as L1 | 0.001 (−0.004 to 0.006) |
| Forecast − pooled lookup | L1, not covered | no pairs | no pairs | −0.122 (−0.188 to −0.055) |

<!-- src: results/stage8/spearman_diff.tsv; results/stage9/f0_lstar_diff.tsv (or the file behind the current "CPTAC mean lookup" row); results/stage10/lpool_diff.tsv; results/stage11/l1_split_diff.tsv. Keep the current CPTAC mean lookup row values if the working table differs from the numbers above, and report the difference. -->

Notes: Differences in Spearman correlation with observed concordance; 95% intervals from 2,000 resamples of antibodies (RPPA), proteins (ProCan) or genes (APOLLO). All pairs: 4,045 (134 antibodies), 6,118 (2,746 proteins) and 7,196 genes. L1, pairs with a matched-cancer CPTAC correlation: 1,206 (127 antibodies), 3,356 (2,091 proteins) and 6,982 genes. On L1 pairs, the cascade lookup equals the matched lookup. Covered, pairs whose gene had stored CPTAC correlations, so that the forecast's gene-level features were gene-specific (6,382 APOLLO genes; all RPPA and ProCan L1 pairs); not covered, pairs whose gene had none (600 APOLLO genes, all outside the P set). The split was defined after the L1 results had been seen. The CPTAC mean lookup used the correlations stored at development; the cascade, matched and pooled lookups used a recomputed set (Methods).

---

# R10. Cover letter

**Find:** "In post hoc comparisons, it was at least as accurate as looking up observed CPTAC correlations, even from the matched cancer type, and a short calibration with 50 locally measured genes made its 80% intervals hold on all three platforms."

**Replace with:**

In post hoc comparisons, it was at least as accurate as looking up observed CPTAC correlations, pooled across cancer types or taken from the matched type, for genes whose CPTAC data it had been given; for genes whose CPTAC correlations the model had not used, a lookup was better, and we report this openly. A short calibration with 50 locally measured genes made its 80% intervals hold on all three platforms.

Also fill "[FILL: "and for N genes" if OPTION A]" with "and for 19,211 genes" [CHECK: forecast_all_genes_v3 gene count].

---

# R11. Supplementary Note 2, "Reproducibility" paragraph (B_SI_srep.md) — replace the whole paragraph

[OPTION R]
**Reproducibility.** For public release, the forecast model was refitted with the frozen data and parameters, because the fitted model objects had not been stored. The model had been fitted twice with identical training rows (SHA-256 `9a00d5dd5b304770929711b407b1c9203072411ad33c465be16e69c1559c4976`), code and parameters: at commit `224779714bffd58b02059f51de48d92d1324a24f` before the RPPA and ProCan predictions were frozen, and at commit `ce7949a323791175ee4d60c0b7d79d8e8c5bda02` before the APOLLO predictions were frozen. The two fits differed in [FILL: setting and values]. Refitting with the setting of each fit reproduced its frozen predictions (maximum absolute differences [FILL] for RPPA, [FILL] for ProCan and [FILL] for APOLLO). Both models are released ([FILL: file names]). Between the two models, RPPA and ProCan predictions differed by at most 0.031 and 0.044 (Spearman correlations 0.9991 and 0.9992), and the preregistered statistics recomputed with the second model differed from the frozen values by −2.9 × 10^−5 and −3.3 × 10^−4. All results in the paper use the frozen predictions.

[OPTION N]
**Reproducibility.** For public release, the forecast model was refitted with the frozen data and parameters, because the fitted model objects had not been stored. The model had been fitted twice with identical training rows (SHA-256 `9a00d5dd5b304770929711b407b1c9203072411ad33c465be16e69c1559c4976`), code and parameters: at commit `224779714bffd58b02059f51de48d92d1324a24f` before the RPPA and ProCan predictions were frozen, and at commit `ce7949a323791175ee4d60c0b7d79d8e8c5bda02` before the APOLLO predictions were frozen. The released model reproduced the frozen APOLLO predictions exactly (maximum absolute difference 2.4 × 10^−16 over 7,318 genes) but not the RPPA and ProCan predictions (0.031 and 0.044; Spearman correlations 0.9991 and 0.9992). The input features were identical (maximum absolute difference 1.4 × 10^−14). [FILL: the settings tested in B12 §2 in one sentence, e.g. "Refits with 1, 16 and 64 threads and with row-wise and column-wise histogram construction did not reproduce them."] The preregistered statistics recomputed with the refitted predictions differed from the frozen values by −2.9 × 10^−5 and −3.3 × 10^−4. All results in the paper use the frozen predictions, which are released.

<!-- src: results/stage11/b4_refit_compare.tsv; results/stage10/feature_diagnosis_summary.tsv; results/stage9/refit_compare.tsv; results/stage12/repro_grid.tsv -->

---

# R12. Supplementary Note 3, post hoc plans (B_SI_srep.md)

Make the last sentence of the post hoc paragraph list, in time order, every post hoc plan commit with its output commit. Add, after the B10 pooled-lookup plan:

…; the split of the APOLLO L1 pairs by the forecast's CPTAC coverage, planned after the pooled-lookup results had been seen, at `773bd579f9c2bb4eb65e1e0f0cf955fbbe71a6eb` (15 October 2026, 11:47:50 +0900), with output at `8720727cadc81ea4fcf488ab1ac536614accda8c`; and the reproduction analysis of the forecast model at [FILL: B12 plan commit and time], with output at [FILL].

---

# R13. Supplementary Note 6 (B_SI_srep.md)

**Replace the sentence** that filled "[FILL: from B10 §2, the number of evaluated pairs per platform for which F0 was the global mean but L* had a CPTAC value.]" **with:**

Among the evaluated pairs with a matched-cancer correlation (L1 pairs), F0 was the global mean while the lookups had a CPTAC value for 600 APOLLO pairs, all for genes outside the P set, and for no RPPA or ProCan pairs; for all 600, the forecast's gene-level features were at their default values. The pooled lookup is level L2 of the cascade, evaluated on L1 pairs.

<!-- src: results/stage11/l1_split_counts.tsv -->

---

# R14. Supplementary Note 7 (B_SI_srep.md)

**Replace** "[FILL: the search method of B10 §4 in two sentences.]" (or the sentences that replaced it) **with** [OPTION F]:

Labels were retrieved through the openFDA drug label and Drugs@FDA interfaces on [FILL: date] with the queries listed in Supplementary Table S7, covering the term "antibody-drug conjugate" and the payload stems of approved and investigational ADC names. For each product approved by 2 October 2026, the target was taken from the mechanism-of-action section of the label; "none found" means that no product in this search acted on the target.

[OPTION U]: keep the current sentences.
