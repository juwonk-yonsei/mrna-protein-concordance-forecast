<!--
B_SI_srep_claude.md — Supplementary Information for the Scientific Reports manuscript, written by Claude for B 지시서 10.
Same tags as the main text: [FILL], [CHECK], [@Key], <!-- src -->. All tags must be gone from the submitted PDF.
Scientific Reports rules: one PDF; items named "Supplementary Note", "Supplementary Table S1", "Supplementary Figure S1"; notes are not numbered by the journal but we number them for cross-reference, which the journal allows as titles. Do not refer to panels of supplementary figures.
-->

# Supplementary Information

**A preregistered model forecasts mRNA–protein concordance in unmeasured cancer contexts across three proteomic platforms**

[AUTHOR TO FILL: authors]

## Supplementary Note 1. Data sources and processing

**CPTAC.** [FILL: freeze name, download date, file names; RNA and protein units; how paired samples were matched; definition of the P gene set (7,432 genes).] The numbers of tumors with both RNA and protein are given in Supplementary Table S1.

**TCGA RPPA.** [FILL: TCPA file and level; RNA matrix used for RPPA samples; sample matching; antibody-to-gene mapping and handling of phospho-antibodies and multi-gene antibodies; minimum samples per cancer type.] The matrix had 7,754 samples and 258 antibodies, and the confirmation file had 4,154 antibody–cancer pairs, of which 4,045 had a defined observed concordance.

**ProCan.** [FILL: protein matrix file; RNA source and matching of cell lines; lineage labels.] The structure table had 8,498 protein columns, and the confirmation file had 12,159 protein–lineage pairs, of which 6,118 had a defined observed concordance.

**APOLLO.** Protein data were the Genomic Data Commons open supplement `APOLLO1_Level3_Gprot_v061418.n87.csv` (SHA-256 `621c04d7ee91d4a8295891c03513c872d731aeae772744045d3748040715f2f4`, 9,060,133 bytes) for 87 patients of the US Lung Cancer Biospecimen Resource Network, 2012–2018, with resected primary tumors[@Soltis2022]. Protein values are log2 TMT ratios to a common reference in channel m/z 126. Quantification used the median of unique peptide-spectrum matches, with at least two unique matches. Proteins quantified in at least 50% of samples had missing values filled by k-nearest neighbors, and proteins quantified in fewer than 50% were excluded (ProteomeXchange PXD036025). RNA files were processed with MapSplice 2.2.2, hg38, HTSeq 0.9.1 and GENCODE v28 basic annotation; TPM used the mean transcript length, with the per-million factor computed over protein-coding genes. We used log2(TPM + 1) without restricting the matrix to coding genes. All 87 samples were kept and no further quality filter was applied.

The processed RNA matrix `data/processed/apollo_luad/rna_log2tpm.parquet` has SHA-256 `24b4aeb2a9bced86a53645a28879f18125f69f24a9e9a443df1ebfd26c3cce96` (20,714,409 bytes), and the sample map `config/apollo_sample_map.tsv` (87 one-to-one identifiers) has SHA-256 `4b9045ad3f9faee678cc85e8769391da939bdf9675d1e72f70476b8064cdb6ef`. The RNA files had 56,832 symbol rows, no duplicate symbols and no Ensembl identifiers. Of the 19,211 genes used as RNA inputs of the anchor-residual predictor, 1,287 were absent from the APOLLO RNA matrix and were given z = 0. Per-sample TPM sums ranged from 1.87 × 10^6 to 2.90 × 10^6 (median 2.18 × 10^6) because non-coding genes remain in the file; the sums were not used as a filter.

Protein rows whose symbol contained `;`, `|` or `,` were dropped, and all rows of a symbol that occupied more than one row were dropped. From 7,631 rows, no multi-symbol rows and 25 rows of 9 duplicated symbols (CUX1, HLA-A, HLA-B, HLA-C, HLA-DRB1, MOCS2, NA, RABGAP1L and TMPO) were removed, leaving 7,606 symbols. The upper bounds of units before finite-value filters were 7,148 for P1, 7,318 for P2 and 6,323 for Q1 and Q2. Of the 7,318 P2 genes, 120 had no RNA row, so their observed concordance was undefined; with two further genes without a defined concordance, 7,196 genes entered P2. Public identifiers of APOLLO and CPTAC lung adenocarcinoma samples did not overlap; enrollment of the same patient under different identifiers cannot be excluded.
<!-- src: reports/B_methods_facts.md; results/stage6/confirm/descriptive_predictors.tsv (122 excluded) [CHECK: "two further genes" = 122 − 120] -->

**Other data.** CCLE lung (74 cell lines; 6,911 genes with an own-mRNA correlation) and skin (33 cell lines; 5,947 genes) lineages were used during development [CHECK: purpose]. File names, sources and SHA-256 hashes of all downloaded files are in Supplementary Table S9.
<!-- src: results/figures/source/table1_data.tsv -->

## Supplementary Note 2. Forecast model and reproducibility of frozen predictions

**Features.** [FILL: table of the nine F1s features with definitions, marking which belong to F1.]

**Model.** LightGBM used 300 trees, learning rate 0.05, 31 leaves, a minimum of 20 samples per leaf, row subsample 0.8 at every iteration, column subsample 0.8 and random state 20261001, with all available threads. Missing feature values were replaced by the column median of the training matrix. The 80% and 90% interval half-widths were 0.202 and 0.259. [CHECK: derivation of the half-widths.]

**Reproducibility.** For public release, the forecast model was refitted with the frozen data and parameters, because the fitted model object had not been stored. The refitted model reproduced the frozen APOLLO predictions exactly (maximum absolute difference 3.1 × 10^−16 over 7,318 genes). For RPPA and ProCan, the maximum absolute differences were 0.031 and 0.044, the Spearman correlations between refitted and frozen predictions were 0.9991 and 0.9992, and the preregistered statistics recomputed with the refitted predictions differed from the frozen values by −2.9 × 10^−5 and −3.3 × 10^−4. [FILL: cause of the RPPA and ProCan differences from B10 §3.] All results in the paper use the frozen predictions.
<!-- src: results/stage9/refit_compare.tsv -->

## Supplementary Note 3. Preregistration history and departures

The RPPA and ProCan tests (H1–H4) were preregistered at commit `f7befd8324bcef63d89cebb57de36c0232b5b18f` (2 October 2026, 09:57:46 +0900), and their confirmation output was committed at `a4764533d08bd09f404f827f7e01726c6a897cff` (10:29:53). For APOLLO, the RNA processing was committed at `ef3ba3ff12bdae349666b1fc32d783ecfa4f673b` (15:48:40), the preregistration at `87967f098ea1434f39c3cd3109d2edeea1c1a0ef` (15:50:05), the frozen predictions at `ce7949a323791175ee4d60c0b7d79d8e8c5bda02` (16:16:58), the code freeze at `759b3b1fc6f49db8398ba17f9769b6818bd67703` (16:58:18), the unlock at `f1c1c6370187456b831abef46078e509014e61fc` (16:59:01) and the confirmation output at `57a3a3816bec154ae0df4f38dfdae5155eedc464` (16:59:47). Each confirmation command was run once.

Departures accepted before preregistration and before unlocking: the first power calculation was discarded, and the rank definition was edited before the counted tests on simulated data.

Holm's procedure was applied to one family of four hypotheses (H1–H4) and to two families of two hypotheses for APOLLO (P1 and P2, primary; Q1 and Q2, secondary) [CHECK: thresholds]. Because the two APOLLO families were controlled separately, the familywise error rate across them was not controlled; the worst case stated in the preregistration is about 0.10. Primary claims rest on the primary family.

Post hoc analysis plans were committed before calculation: the lookup plan at `703c94594efb81f2764ec1ec5c442398ef52b0c9` (13 October 2026, 16:44:38 +0900; an earlier commit, `c33e5aaf54d4608b763d4b10c4d752183aed5045` at 16:41:11, did not parse and was corrected before any calculation, without changes to the context map or cascade), with output at `9984b2dc569d3f1a0f6a6b2b8b5a76852f1db86e` (16:57:21); the comparison of the CPTAC mean lookup with the cascade lookup at `a77b269950af12c1dc061c7cdf8f2fd2caa51775` (18:05:17); and [FILL: B10 §2 plan commit and time].

Commit times were recorded by git on the analysis workstation and were not certified by a third party. The published history fixes the order of commits through their hashes but does not independently prove their times.

## Supplementary Note 4. Power analysis and tests on simulated data

The power file used in the APOLLO preregistration was the second run. All nine scenarios had a rejection rate of 1.0 and no Wilcoxon errors. The P1 source had n = 7,418, mean 0.0969 and median 0.0624, with draw size 7,148. The Q1 source had n = 7,403, mean 0.0319 and draw size 6,323, and the Q2 source n = 7,403, mean 0.112 and draw size 6,323. Fifteen anchor-residual genes lacked leave-one-cancer-out R1 or R2 predictions and were excluded from Q only. The CPTAC lung adenocarcinoma reference for P2 had n = 7,418, a forecast Spearman correlation of 0.845, an RNA-variance baseline correlation of 0.386 and a difference of 0.459.
<!-- src: results/stage6/prep/power.tsv; results/stage6/prep/p2_cptac_luad_reference.tsv -->

On planted data, the primary tests rejected: P1 had n = 7,029, mean difference 0.0643 and stored P = 0, and P2 had n = 7,196, a forecast correlation of 0.928 against 0.544 for the baseline (difference 0.384, P = 0.0005). On null data, neither test rejected; the absolute mean difference was 0.000105, and the absolute Spearman correlations were 0.0229 for the forecast and 0.0162 for the baseline. The planted and null data sets were local checks and were not committed.
<!-- src: results/stage6/fake_planted/; results/stage6/fake_null/ -->

## Supplementary Note 5. Sensitivity analyses

**Non-imputed APOLLO matrix (post hoc).** The non-imputed protein file (GDC publication page APOLLO-LUAD-2022, file `48ef3b6e-4546-4833-ab18-debc2fa81a3c`, SHA-256 `71210627045eb233613ce7e3035f44e5aecb2d24c7da5685dcb551bb0bb4f549`, 10,954,279 bytes, downloaded 2 October 2026) was analysed after excluding 13 sample columns outside the sample map, blank gene-symbol rows and duplicated symbols. The P1-style mean difference was 0.0847 (n = 7,027; stored P = 0) and the P2-style difference was 0.231 (forecast 0.752 against baseline 0.521; n = 7,194; P = 0.0005). These are not preregistered tests.
<!-- src: results/stage6/posthoc/posthoc_05_raw.tsv -->

**Anchor-free baselines on RPPA and ProCan (post hoc).** Mean ρ_AR − ρ_R1 was 0.0250 on RPPA (n = 100) and −0.0489 on ProCan (n = 2,399); mean ρ_AR − ρ_R2 was −0.0343 on RPPA and −0.103 on ProCan. These stage-5 summaries have no intervals. Refitting R1 and R2 with the stage-5 function reproduced the stored coefficients (n = 4,998; maximum absolute difference 9.97 × 10^−17).
<!-- src: results/stage5/posthoc/posthoc_23_summary.tsv; reports/B_methods_facts.md -->

**Pre-specified descriptive summaries of the APOLLO anchor-residual differences.** Median 0.0556, first quartile 0.0106, third quartile 0.131, fraction above 0.05 0.529 and fraction below −0.05 0.0414 (n = 7,029). The median was 0.0151 outside the P set and 0.0623 inside it. By the selected anchor weight λ, the median difference was 0 at λ = 0, 0.00779 at 0.25, 0.0261 at 0.5, 0.0520 at 0.75 and 0.0753 at 1.
<!-- src: results/stage6/confirm/descriptive_ar.tsv -->

**Similarity of APOLLO to CPTAC.** The preregistered descriptive similarity s was 0.925, inside the CPTAC holdout range of 0.883–0.979.
<!-- src: results/stage6/confirm/s.tsv -->

## Supplementary Note 6. Lookup definitions and coverage

The CPTAC mean lookup F0 used the CPTAC correlations stored during model development. Among RPPA antibodies, 100 genes were in the P set, 28 outside it had a gene-specific F0 and 6 received the global mean; among ProCan proteins, the counts were 3,052, 819 and 182; and among APOLLO genes, 6,323, 201 and 794. These counts include genes whose observed concordance is undefined.
<!-- src: results/stage7/f0_fill_counts.tsv -->

The cascade lookup L* recomputed CPTAC correlations for all genes quantified in each cohort, which added 9,318 gene–cohort correlations. [FILL: from B10 §2, the number of evaluated pairs per platform for which F0 was the global mean but L* had a CPTAC value.] The context map is in Supplementary Table S10, and the numbers of pairs in each cascade level are in Supplementary Table S8. With only two RPPA antibodies, comparisons in level L3 were unstable (1,034 of 2,000 resamples undefined) and are not interpreted; comparisons in level L4, where the lookup is constant, are undefined.
<!-- src: results/stage8/lookup_build.tsv; results/stage8/spearman_diff.tsv -->

## Supplementary Note 7. Antibody–drug conjugate atlas

The atlas has 957 rows, 29 targets in each of 33 TCGA cancer types. Each row gives the frozen forecast and its 80% interval (the forecast plus or minus the stored half-width, not clipped). Every CLDN18 row states that gene-level assays cannot separate isoforms and that the ADC target is CLDN18.2. Column definitions are in the released file [FILL: path in the release].
<!-- src: results/stage7/adc_atlas_final.tsv; results/stage7/adc_columns.tsv -->

In a descriptive comparison with the Human Protein Atlas[@Uhlen2017], 24 of the 29 targets had a defined correlation, and the Spearman correlation between that concordance and the mean forecast was 0.0443. The two quantities differ in kind: the Human Protein Atlas compares staining levels between cancers, whereas the forecast is a within-cancer correlation across patients.
<!-- src: results/stage5/atlas/hpa_summary.tsv -->

Approved ADCs were identified from US Food and Drug Administration prescribing information as of 2 October 2026. [FILL: the search method of B10 §4 in two sentences.] The result is in Supplementary Table S7.

## Supplementary tables

**Supplementary Table S1. CPTAC cohorts.**

| Cancer type | Tumors with RNA and protein | Cohort reference |
|---|---|---|
| BRCA | 121 | [@Krug2020] |
| CCRCC | 103 | [@Clark2019] |
| COAD | 96 | [@Vasaikar2019] |
| GBM | 99 | [@WangLB2021] |
| HNSCC | 108 | [@Huang2021] |
| LSCC | 108 | [@Satpathy2021] |
| LUAD | 110 | [@Gillette2020] |
| OV | 82 | [@McDermott2020] |
| PDAC | 105 | [@Cao2021] |
| UCEC | 95 | [@Dou2020] |
| Total | 1,027 | |

<!-- src: results/cptac_sample_counts.tsv -->

**Supplementary Table S2. Forecast minus CPTAC mean lookup by membership in the P gene set (post hoc).**

| Model | Genes | RPPA | ProCan | APOLLO |
|---|---|---|---|---|
| F1s-LightGBM | All | 0.0643 (0.0213 to 0.112); 4,045 pairs | 0.0317 (0.0195 to 0.0440); 6,118 pairs | 0.0276 (0.0206 to 0.0348); 7,196 genes |
| F1s-LightGBM | Inside P | 0.0718 (0.0349 to 0.117); 3,051 | 0.0208 (0.0106 to 0.0308); 5,547 | 0.00187 (−0.00323 to 0.00678); 6,221 |
| F1s-LightGBM | Outside P | 0.0405 (−0.0912 to 0.204); 994 (34 antibodies) | 0.139 (0.0702 to 0.209); 571 (347 proteins) | 0.276 (0.220 to 0.328); 975 |
| F1s-ridge | All | 0.0573 (0.0266 to 0.0914) | 0.0386 (0.0295 to 0.0470) | 0.00462 (−0.00238 to 0.0117) |
| F1s-ridge | Inside P | [FILL] | [FILL] | [FILL] |
| F1s-ridge | Outside P | [FILL] | [FILL] | [FILL] |

<!-- src: results/stage7/f0_diff.tsv -->

Notes: Differences in Spearman correlation with observed concordance (95% interval from 2,000 unit resamples).

**Supplementary Table S3. Discrimination of genes with observed concordance of at least 0.5 (post hoc).**

| | RPPA | ProCan | APOLLO |
|---|---|---|---|
| Pairs (positives) | 4,045 (951) | 6,118 (1,993) | 7,196 (3,188) |
| AUROC, forecast | 0.746 (0.701–0.787) | 0.738 (0.723–0.753) | 0.868 (0.860–0.876) |
| AUROC, CPTAC mean lookup | 0.696 (0.636–0.748) | 0.725 (0.709–0.739) | 0.857 (0.848–0.865) |
| AUROC, RNA-variance baseline | 0.677 (0.630–0.723) | 0.620 (0.603–0.637) | 0.753 (0.742–0.764) |
| Forecast − CPTAC mean lookup, all pairs | 0.0498 (0.0264 to 0.0774) | 0.0136 (0.00739 to 0.0195) | 0.0112 (0.00674 to 0.0157) |
| Forecast − cascade lookup, all pairs | 0.0457 (0.0200 to 0.0756) | 0.0140 (0.00352 to 0.0239) | 0.0110 (0.00396 to 0.0180) |
| Forecast − cascade lookup, matched cancer (L1) | 0.00798 (−0.0238 to 0.0399) | 0.0173 (0.00245 to 0.0329) | 0.0100 (0.00274 to 0.0173) |
| Forecast − RNA-variance baseline, all pairs | 0.0683 (0.0157 to 0.120) | 0.119 (0.103 to 0.135) | 0.115 (0.104 to 0.125) |

<!-- src: results/stage7/f0_auroc.tsv; results/stage8/auroc_diff.tsv -->

**Supplementary Table S4. Within-context and within-antibody comparisons (post hoc).**

| Comparison | Value |
|---|---|
| Median within-context Spearman correlation, forecast: RPPA (31 cancer types), ProCan (3 lineages), APOLLO | 0.439, 0.440, 0.747 |
| Same, CPTAC mean lookup | 0.398, 0.404, 0.719 |
| Same, RNA-variance baseline | 0.246, 0.135, 0.492 |
| RPPA, forecast − CPTAC mean lookup within cancer type | 30 of 31 positive; median 0.0600 (0.0498 to 0.0773); sign test P = 3.0 × 10^−8 |
| RPPA, forecast − cascade lookup within cancer type | 28 of 31 positive; median 0.0799 (0.0574 to 0.0982); P = 4.6 × 10^−6 |
| RPPA within antibody: median Spearman correlation across cancer types, forecast and RNA-variance baseline (130 antibodies) | 0.255 and 0.203; median difference 0.0242 (−0.0161 to 0.0959) |

<!-- src: results/stage7/f0_context_median.tsv; results/stage8/within_context_sign.tsv; results/stage7/f0_within_antibody_summary.tsv -->

**Supplementary Table S5. Pre-specified descriptive metrics of all predictors on APOLLO (7,196 genes).**

| Predictor | Spearman | AUROC (ρ ≥ 0.5) | Pearson | Mean absolute error | Bias |
|---|---|---|---|---|---|
| F1s-LightGBM (forecast) | 0.747 | 0.868 | 0.750 | 0.116 | −0.0227 |
| F1s-ridge | 0.724 | [FILL] | [FILL] | [FILL] | −0.0242 |
| F1-LightGBM | 0.627 | [FILL] | [FILL] | [FILL] | 0.189 |
| F1-ridge | 0.719 | [FILL] | [FILL] | [FILL] | 0.155 |
| CPTAC mean lookup (F0) | 0.719 | 0.857 | 0.726 | 0.121 | −0.0124 |
| RNA-variance baseline | 0.492 | 0.753 | n/a | n/a | n/a |

<!-- src: results/stage6/confirm/descriptive_predictors.tsv. Bias = mean of predicted minus observed or observed minus predicted: [CHECK sign convention and state it in the note]. -->

Notes: [FILL: bias definition]. One hundred and twenty-two genes of the 7,318 P2 genes had no defined observed concordance and were excluded.

**Supplementary Table S6. Coverage of 80% intervals by calibration method and number of calibration genes K.**

[FILL: from results/stage7/calibration_table.tsv, F1s-LightGBM rows only: platform, method, K, mean coverage (2.5th–97.5th percentile of draws), mean width. Also list the F1s-ridge rows in a second block, labelled as not the preregistered predictor.]

**Supplementary Table S7. FDA-approved ADCs acting on atlas targets (as of 2 October 2026).**

[FILL: 29 rows from B10 §4: gene, protein, approved ADC (name, BLA, label URL) or "none found" with the evidence, search date.]

**Supplementary Table S8. Pairs in each level of the cascade lookup.**

| Level | RPPA pairs (antibodies) | ProCan pairs (proteins) | APOLLO genes |
|---|---|---|---|
| L1, matched cancer | 1,206 (127) | 3,356 (2,091) | 6,982 |
| L2, CPTAC mean | 2,653 (127) | 2,616 (2,570) | 93 |
| L3, ProCan mean | 62 (2) | not used | 26 |
| L4, global mean | 124 (4) | 146 (79) | 95 |
| Total | 4,045 (134) | 6,118 (2,746) | 7,196 |

<!-- src: results/stage9/f0_lstar_counts.tsv -->

**Supplementary Table S9. Downloaded files.**

[FILL: file name, source URL, release or version, download date, size, SHA-256, for every downloaded file used in the paper, from logs/download_log.tsv.]

**Supplementary Table S10. Context map of the matched lookup (L1).**

| Evaluated data set | Context | Matched CPTAC cohort |
|---|---|---|
| RPPA | breast invasive carcinoma | BRCA |
| RPPA | kidney clear cell carcinoma | CCRCC |
| RPPA | colon adenocarcinoma; rectum adenocarcinoma | COAD |
| RPPA | glioblastoma multiforme | GBM |
| RPPA | head and neck squamous cell carcinoma | HNSCC |
| RPPA | lung adenocarcinoma | LUAD |
| RPPA | lung squamous cell carcinoma | LSCC |
| RPPA | ovarian serous cystadenocarcinoma | OV |
| RPPA | pancreatic adenocarcinoma | PDAC |
| RPPA | uterine corpus endometrioid carcinoma | UCEC |
| RPPA | other cancer types | none |
| ProCan | Large Intestine | COAD |
| ProCan | Lung | LUAD and LSCC (mean of finite values) |
| ProCan | Haematopoietic and Lymphoid | none |
| APOLLO | lung adenocarcinoma | LUAD |

<!-- src: config/analysis_plan_B8.yaml -->

## Supplementary figures

**Supplementary Figure S1. Forecast minus CPTAC mean lookup by membership in the P gene set (post hoc).** Differences in Spearman correlation with 95% intervals for F1s-LightGBM and F1s-ridge on all pairs, pairs inside P and pairs outside P, for RPPA, ProCan and APOLLO. Values are in Supplementary Table S2.

**Supplementary Figure S2. Within-cancer-type comparison on RPPA (post hoc).** Spearman correlation of the forecast against that of the cascade lookup with observed concordance in each of the 31 RPPA cancer types; the line marks equality. [FILL: source file results/stage8/within_context.tsv]

## Supplementary references

[Converted from citation keys; numbered separately from the main text, or merged into the main reference list if the journal prefers. CHECK the Scientific Reports rule and follow it.]
