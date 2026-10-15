# Supplementary Information

**A preregistered model forecasts mRNA–protein concordance in unmeasured cancer contexts across three proteomic platforms**

Juwon Kang^1 and Junjeong Choi^1*<br><br>^1 College of Pharmacy, Yonsei Institute of Pharmaceutical Sciences, Yonsei University, Incheon 21983, Republic of Korea

## Supplementary Note 1. Data sources and processing

**CPTAC.** Paired tumor RNA and protein were the CPTAC pan-cancer data freeze v1.2, downloaded on 1 October 2026. RNA files are RSEM coding upper-quartile log2 abundances and protein files are log2 ratios to a reference intensity. A sample entered the paired set when it was listed as a tumor with both an RNA column and a protein column. The P gene set contains the 7,432 genes measured in at least 8 of the 10 cancer types at a tumor measurement rate of at least 0.70 and present in the RNA gene universe. The numbers of tumors with both RNA and protein are given in Supplementary Table S1.

**TCGA RPPA.** The protein matrix was the PanCanAtlas Xena file TCGA-RPPA-pancan-clean, downloaded on 1 October 2026. A TCPA data level is not recorded for this file. RNA was the PanCanAtlas batch-corrected RNA-seq matrix. Samples were the intersection of the two matrices after duplicate RPPA sample identifiers were dropped. Antibody names containing `_p`, `pS`, `pT`, `pY`, `cleaved` or `_c` were excluded, and an antibody mapped to more than one gene symbol was dropped. A cancer type was evaluated when at least 30 samples had both measurements. The matrix had 7,754 samples and 258 antibodies, and the confirmation file had 4,154 antibody–cancer pairs, of which 4,045 had a defined observed concordance.

**ProCan.** The protein matrix was ProCan-DepMapSanger_protein_matrix_8498_averaged.txt (figshare file 34411172). RNA was DepMap 24Q4 log2(TPM + 1), joined by Sanger model identifier. Lineage labels are the Tissue_type field of the ProCan mapping file. The structure table had 8,498 protein columns, and the confirmation file had 12,159 protein–lineage pairs, of which 6,118 had a defined observed concordance.

**APOLLO.** Protein data were the Genomic Data Commons open supplement `APOLLO1_Level3_Gprot_v061418.n87.csv` (SHA-256 `621c04d7ee91d4a8295891c03513c872d731aeae772744045d3748040715f2f4`, 9,060,133 bytes) for 87 patients of the US Lung Cancer Biospecimen Resource Network, 2012–2018, with resected primary tumors<sup>1</sup>. Protein values are log2 TMT ratios to a common reference in channel m/z 126. Quantification used the median of unique peptide-spectrum matches, with at least two unique matches. Proteins quantified in at least 50% of samples had missing values filled by k-nearest neighbors, and proteins quantified in fewer than 50% were excluded (ProteomeXchange PXD036025). RNA files were processed with MapSplice 2.2.2, hg38, HTSeq 0.9.1 and GENCODE v28 basic annotation; TPM used the mean transcript length, with the per-million factor computed over protein-coding genes. We used log2(TPM + 1) without restricting the matrix to coding genes. All 87 samples were kept and no further quality filter was applied.

The processed RNA matrix `data/processed/apollo_luad/rna_log2tpm.parquet` has SHA-256 `24b4aeb2a9bced86a53645a28879f18125f69f24a9e9a443df1ebfd26c3cce96` (20,714,409 bytes), and the sample map `config/apollo_sample_map.tsv` (87 one-to-one identifiers) has SHA-256 `4b9045ad3f9faee678cc85e8769391da939bdf9675d1e72f70476b8064cdb6ef`. The RNA files had 56,832 symbol rows, no duplicate symbols and no Ensembl identifiers. Of the 19,211 genes used as RNA inputs of the anchor-residual predictor, 1,287 were absent from the APOLLO RNA matrix and were given z = 0. Per-sample TPM sums ranged from 1.87 × 10^6 to 2.90 × 10^6 (median 2.18 × 10^6) because non-coding genes remain in the file; the sums were not used as a filter.

Protein rows whose symbol contained `;`, `|` or `,` were dropped, and all rows of a symbol that occupied more than one row were dropped. From 7,631 rows, no multi-symbol rows and 25 rows of 9 duplicated symbols (CUX1, HLA-A, HLA-B, HLA-C, HLA-DRB1, MOCS2, NA, RABGAP1L and TMPO) were removed, leaving 7,606 symbols. The upper bounds of units before finite-value filters were 7,148 for P1, 7,318 for P2 and 6,323 for Q1 and Q2. Of the 7,318 P2 genes, 120 had no RNA row, so their observed concordance was undefined; with two further genes without a defined concordance, 7,196 genes entered P2. Public identifiers of APOLLO and CPTAC lung adenocarcinoma samples did not overlap; enrollment of the same patient under different identifiers cannot be excluded.


**Other data.** CCLE lung (74 cell lines; 6,911 genes with an own-mRNA correlation) and skin (33 cell lines; 5,947 genes) lineages were used during development to compare the anchor-residual predictor with the gene's own mRNA in cell lines. File names, sources and SHA-256 hashes of all downloaded files are in Supplementary Table S9.


## Supplementary Note 2. Forecast model and reproducibility of frozen predictions

**Features.** F1s uses nine features. F1, fitted in the preceding stage, uses the five correlation summaries, the query RNA mean and standard deviation and the percentile of that standard deviation, the mean of the training RNA means, and the training protein measurement rate.

| Feature | Definition | Model |
|---|---|---|
| rho_mean | Mean of the gene's defined CPTAC correlations in the training cancers | F1, F1s |
| rho_sd | Standard deviation of those correlations (divisor n; 0 if fewer than two) | F1, F1s |
| rho_min | Minimum of those correlations | F1, F1s |
| rho_max | Maximum of those correlations | F1, F1s |
| n_defined_cancers | Number of training cancers with a defined correlation | F1, F1s |
| train_protein_rate_mean | Mean fraction of paired tumors in which the protein was measured | F1, F1s |
| holdout_rna_sd_percentile | Percentile of the gene's RNA standard deviation in the query context | F1, F1s |
| train_rna_mean | Mean of the gene's per-cancer RNA means in the training cancers | F1 |
| holdout_rna_mean | RNA mean in the query context | F1 |
| holdout_rna_sd | RNA standard deviation in the query context | F1 |
| train_rna_mean_percentile_mean | Mean, across training cancers, of the percentile of the gene's RNA mean | F1s |
| holdout_rna_mean_percentile | Percentile of the gene's RNA mean in the query context | F1s |

**Model.** LightGBM used 300 trees, learning rate 0.05, 31 leaves, a minimum of 20 samples per leaf, row subsample 0.8 at every iteration, column subsample 0.8 and random state 20261001, with all available threads. Missing feature values were replaced by the column median of the training matrix. The 80% and 90% interval half-widths were 0.202 and 0.259, the conformal quantiles of absolute residuals on the stacked leave-one-cancer-out rows.

**Reproducibility.** For public release, the forecast model was refitted with the frozen data and parameters, because the fitted model objects had not been stored. The model had been fitted twice with identical training rows, the 72,880 rows of `results/f1s_predictions.tsv` (SHA-256 `3ab5b9378020ef2e3710ba70e6e817416d9a20e055ca8c8086c64b02071334d6`) that are in the P gene set with a finite observed correlation, code and parameters: at commit `224779714bffd58b02059f51de48d92d1324a24f` before the RPPA and ProCan predictions were frozen, and at commit `ce7949a323791175ee4d60c0b7d79d8e8c5bda02` before the APOLLO predictions were frozen. The released model reproduced the frozen APOLLO predictions exactly (maximum absolute difference 2.4 × 10^−16 over 7,318 genes) but not the RPPA and ProCan predictions (0.031 and 0.044; Spearman correlations 0.9991 and 0.9992). The input features were identical (maximum absolute difference 1.4 × 10^−14). Refits with 1, 16 and 72 threads and with automatic, row-wise and column-wise histogram construction all produced the same predictions, and every setting was deterministic across repeats, so neither thread count nor histogram construction explains the difference. The script that wrote the frozen predictions is byte-identical to the released one, the training rows come from the same file, and LightGBM 4.7.0 was installed on 1 October 2026 and recorded both before and after the two fits. We could not identify the cause. The preregistered statistics recomputed with the released model differed from the frozen values by −2.9 × 10^−5 on RPPA and −3.3 × 10^−4 on ProCan. All results in the paper use the frozen predictions, which are released. Code that fits a model on a table in memory and later reads the same table from disk should treat the two as different inputs.


## Supplementary Note 3. Preregistration history and departures

The RPPA and ProCan tests (H1–H4) were preregistered at commit `f7befd8324bcef63d89cebb57de36c0232b5b18f` (2 October 2026, 09:57:46 +0900), and their confirmation output was committed at `a4764533d08bd09f404f827f7e01726c6a897cff` (10:29:53). For APOLLO, the RNA processing was committed at `ef3ba3ff12bdae349666b1fc32d783ecfa4f673b` (15:48:40), the preregistration at `87967f098ea1434f39c3cd3109d2edeea1c1a0ef` (15:50:05), the frozen predictions at `ce7949a323791175ee4d60c0b7d79d8e8c5bda02` (16:16:58), the code freeze at `759b3b1fc6f49db8398ba17f9769b6818bd67703` (16:58:18), the unlock at `f1c1c6370187456b831abef46078e509014e61fc` (16:59:01) and the confirmation output at `57a3a3816bec154ae0df4f38dfdae5155eedc464` (16:59:47). Each confirmation command was run once.

Departures accepted before preregistration and before unlocking: the first power calculation was discarded, and the rank definition was edited before the counted tests on simulated data.

Holm's procedure was applied to one family of four hypotheses (H1–H4) and to two families of two hypotheses for APOLLO (P1 and P2, primary; Q1 and Q2, secondary). Ordered thresholds are 0.0125, 0.05/3, 0.025 and 0.05 for H1–H4, and 0.025 then 0.05 within each APOLLO family. Because the two APOLLO families were controlled separately, the familywise error rate across them was not controlled; the worst case stated in the preregistration is about 0.10. Primary claims rest on the primary family.

Post hoc analysis plans were committed before calculation: the lookup plan at `703c94594efb81f2764ec1ec5c442398ef52b0c9` (13 October 2026, 16:44:38 +0900; an earlier commit, `c33e5aaf54d4608b763d4b10c4d752183aed5045` at 16:41:11, did not parse and was corrected before any calculation, without changes to the context map or cascade), with output at `9984b2dc569d3f1a0f6a6b2b8b5a76852f1db86e` (16:57:21); the comparison of the CPTAC mean lookup with the cascade lookup at `a77b269950af12c1dc061c7cdf8f2fd2caa51775` (18:05:17); and the pooled-lookup plan at `68a32e5974eb0426032c8310d13ec7d716fc27c7` (14 October 2026, 10:38:30 +0900); the split of the APOLLO L1 pairs by the forecast's CPTAC coverage, planned after the pooled-lookup results had been seen, at `773bd579f9c2bb4eb65e1e0f0cf955fbbe71a6eb` (15 October 2026, 11:47:50 +0900), with output at `8720727cadc81ea4fcf488ab1ac536614accda8c`; and the reproduction analysis of the forecast model at `32d4e7191a140668a27cd716d728e63daa73aa71` (15 October 2026, 14:07:53 +0900), with output at `1c02e3776e87334edf9e2b71690717b958b1eb31`.

Commit times were recorded by git on the analysis workstation and were not certified by a third party. The published history fixes the order of commits through their hashes but does not independently prove their times.

## Supplementary Note 4. Power analysis and tests on simulated data

The power file used in the APOLLO preregistration was the second run. All nine scenarios had a rejection rate of 1.0 and no Wilcoxon errors. The P1 source had n = 7,418, mean 0.0969 and median 0.0624, with draw size 7,148. The Q1 source had n = 7,403, mean 0.0319 and draw size 6,323, and the Q2 source n = 7,403, mean 0.112 and draw size 6,323. Fifteen anchor-residual genes lacked leave-one-cancer-out R1 or R2 predictions and were excluded from Q only. The CPTAC lung adenocarcinoma reference for P2 had n = 7,418, a forecast Spearman correlation of 0.845, an RNA-variance baseline correlation of 0.386 and a difference of 0.459.


On planted data, the primary tests rejected: P1 had n = 7,029, mean difference 0.0643 and stored P = 0, and P2 had n = 7,196, a forecast correlation of 0.928 against 0.544 for the baseline (difference 0.384, P = 0.0005). On null data, neither test rejected; the absolute mean difference was 0.000105, and the absolute Spearman correlations were 0.0229 for the forecast and 0.0162 for the baseline. The planted and null data sets were local checks and were not committed.


## Supplementary Note 5. Sensitivity analyses

**Non-imputed APOLLO matrix (post hoc).** The non-imputed protein file (GDC publication page APOLLO-LUAD-2022, file `48ef3b6e-4546-4833-ab18-debc2fa81a3c`, SHA-256 `71210627045eb233613ce7e3035f44e5aecb2d24c7da5685dcb551bb0bb4f549`, 10,954,279 bytes, downloaded 2 October 2026) was analysed after excluding 13 sample columns outside the sample map, blank gene-symbol rows and duplicated symbols. The P1-style mean difference was 0.0847 (n = 7,027; stored P = 0) and the P2-style difference was 0.231 (forecast 0.752 against baseline 0.521; n = 7,194; P = 0.0005). These are not preregistered tests.


**Anchor-free baselines on RPPA and ProCan (post hoc).** Mean ρ_AR − ρ_R1 was 0.0250 on RPPA (n = 100) and −0.0489 on ProCan (n = 2,399); mean ρ_AR − ρ_R2 was −0.0343 on RPPA and −0.103 on ProCan. These stage-5 summaries have no intervals. Refitting R1 and R2 with the stage-5 function reproduced the stored coefficients (n = 4,998; maximum absolute difference 9.97 × 10^−17).


**Pre-specified descriptive summaries of the APOLLO anchor-residual differences.** Median 0.0556, first quartile 0.0106, third quartile 0.131, fraction above 0.05 0.529 and fraction below −0.05 0.0414 (n = 7,029). The median was 0.0151 outside the P set and 0.0623 inside it. By the selected anchor weight λ, the median difference was 0 at λ = 0, 0.00779 at 0.25, 0.0261 at 0.5, 0.0520 at 0.75 and 0.0753 at 1.


**Similarity of APOLLO to CPTAC.** The preregistered descriptive similarity s was 0.925, inside the CPTAC holdout range of 0.883–0.979.


## Supplementary Note 6. Lookup definitions and coverage

The CPTAC mean lookup F0 used the CPTAC correlations stored during model development. Among RPPA antibodies, 100 genes were in the P set, 28 outside it had a gene-specific F0 and 6 received the global mean; among ProCan proteins, the counts were 3,052, 819 and 182; and among APOLLO genes, 6,323, 201 and 794. These counts include genes whose observed concordance is undefined.


The cascade lookup L* recomputed CPTAC correlations for all genes quantified in each cohort, which added 9,318 gene–cohort correlations. Among pairs with a finite observed concordance, F0 was within 10^−8 of the global mean while L* was a CPTAC value (levels L1–L3) for 681 APOLLO genes, 62 RPPA pairs (2 antibodies) and 0 ProCan pairs. Among the evaluated pairs with a matched-cancer correlation (L1 pairs), F0 was the global mean while the lookups had a CPTAC value for 600 APOLLO pairs, all for genes outside the P set, and for no RPPA or ProCan pairs; for all 600, the forecast's gene-level features were at their default values. The pooled lookup is level L2 of the cascade, evaluated on L1 pairs. The context map is in Supplementary Table S10, and the numbers of pairs in each cascade level are in Supplementary Table S8. With only two RPPA antibodies, comparisons in level L3 were unstable (1,034 of 2,000 resamples undefined) and are not interpreted; comparisons in level L4, where the lookup is constant, are undefined.


## Supplementary Note 7. Antibody–drug conjugate atlas

The atlas has 957 rows, 29 targets in each of 33 TCGA cancer types. Each row gives the frozen forecast and its 80% interval (the forecast plus or minus the stored half-width, not clipped). Every CLDN18 row states that gene-level assays cannot separate isoforms and that the ADC target is CLDN18.2. Column definitions accompany `results/stage7/adc_atlas_final.tsv` in the release. The 29 targets were fixed before any analysis of them and are not the set of approved conjugate targets: the label search found approved conjugates against seven further targets (TNFRSF8, CD22, CD33, CD79B, CD19, TNFRSF17 and IL3RA), whose forecasts are in the released genome-wide table.


In a descriptive comparison with the Human Protein Atlas<sup>2</sup>, 24 of the 29 targets had a defined correlation, and the Spearman correlation between that concordance and the mean forecast was 0.0443. The two quantities differ in kind: the Human Protein Atlas compares staining levels between cancers, whereas the forecast is a within-cancer correlation across patients.


Labels were retrieved through the openFDA drug label and Drugs@FDA interfaces on 15 October 2026 with the queries listed in Supplementary Table S7, covering the term "antibody-drug conjugate" and the payload stems of approved and investigational ADC names. For each product approved by 2 October 2026, the target was taken from the mechanism-of-action section of the label; "none found" means that no product in this search acted on the target.

## Supplementary tables

**Supplementary Table S1. CPTAC cohorts.**

| Cancer type | Tumors with RNA and protein | Cohort reference |
|---|---|---|
| BRCA | 121 | <sup>3</sup> |
| CCRCC | 103 | <sup>4</sup> |
| COAD | 96 | <sup>5</sup> |
| GBM | 99 | <sup>6</sup> |
| HNSCC | 108 | <sup>7</sup> |
| LSCC | 108 | <sup>8</sup> |
| LUAD | 110 | <sup>9</sup> |
| OV | 82 | <sup>10</sup> |
| PDAC | 105 | <sup>11</sup> |
| UCEC | 95 | <sup>12</sup> |
| Total | 1,027 | |



**Supplementary Table S2. Forecast minus CPTAC mean lookup by membership in the P gene set (post hoc).**

| Model | Genes | RPPA | ProCan | APOLLO |
|---|---|---|---|---|
| F1s-LightGBM | All | 0.0643 (0.0213 to 0.112); 4,045 pairs | 0.0317 (0.0195 to 0.0440); 6,118 pairs | 0.0276 (0.0206 to 0.0348); 7,196 genes |
| F1s-LightGBM | Inside P | 0.0718 (0.0349 to 0.117); 3,051 | 0.0208 (0.0106 to 0.0308); 5,547 | 0.00187 (−0.00323 to 0.00678); 6,221 |
| F1s-LightGBM | Outside P | 0.0405 (−0.0912 to 0.204); 994 (34 antibodies) | 0.139 (0.0702 to 0.209); 571 (347 proteins) | 0.276 (0.220 to 0.328); 975 |
| F1s-ridge | All | 0.0573 (0.0266 to 0.0914) | 0.0386 (0.0295 to 0.0470) | 0.00462 (−0.00238 to 0.0117) |
| F1s-ridge | Inside P | 0.0573 (0.02932 to 0.09418); 3,051 | 0.02943 (0.02216 to 0.03676); 5,547 | 0.007169 (0.003163 to 0.01106); 6,221 |
| F1s-ridge | Outside P | 0.04224 (-0.0519 to 0.1495); 994 (34 antibodies) | 0.08092 (0.02566 to 0.1354); 571 (347 proteins) | 0.1577 (0.1048 to 0.2098); 975 |



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



**Supplementary Table S4. Within-context and within-antibody comparisons (post hoc).**

| Comparison | Value |
|---|---|
| Median within-context Spearman correlation, forecast: RPPA (31 cancer types), ProCan (3 lineages), APOLLO | 0.439, 0.440, 0.747 |
| Same, CPTAC mean lookup | 0.398, 0.404, 0.719 |
| Same, RNA-variance baseline | 0.246, 0.135, 0.492 |
| RPPA, forecast − CPTAC mean lookup within cancer type | 30 of 31 positive; median 0.0600 (0.0498 to 0.0773); sign test P = 3.0 × 10^−8 |
| RPPA, forecast − cascade lookup within cancer type | 28 of 31 positive; median 0.0799 (0.0574 to 0.0982); P = 4.6 × 10^−6 |
| RPPA within antibody: median Spearman correlation across cancer types, forecast and RNA-variance baseline (130 antibodies) | 0.255 and 0.203; median difference 0.0242 (−0.0161 to 0.0959) |



**Supplementary Table S5. Pre-specified descriptive metrics of all predictors on APOLLO (7,196 genes).**

| Predictor | Spearman | AUROC (ρ ≥ 0.5) | Pearson | Mean absolute error | Bias |
|---|---|---|---|---|---|
| F1s-LightGBM (forecast) | 0.747 | 0.868 | 0.750 | 0.116 | −0.0227 |
| F1s-ridge | 0.724 | 0.853 | 0.730 | 0.120 | −0.0242 |
| F1-LightGBM | 0.627 | 0.808 | 0.620 | 0.222 | 0.189 |
| F1-ridge | 0.719 | 0.852 | 0.722 | 0.186 | 0.155 |
| CPTAC mean lookup (F0) | 0.719 | 0.857 | 0.726 | 0.121 | −0.0124 |
| RNA-variance baseline | 0.492 | 0.753 | n/a | n/a | n/a |



Notes: Bias is the mean of observed minus predicted concordance. One hundred and twenty-two genes of the 7,318 P2 genes had no defined observed concordance and were excluded.

**Supplementary Table S6. Coverage of 80% intervals by calibration method and number of calibration genes K.**

| Platform | Method | K | Mean 80% coverage (2.5th–97.5th) | Mean width |
|---|---|---|---|---|
| F1s-LightGBM, the preregistered forecast | | | | |
| RPPA | No adjustment | 10 | 0.511 (0.491–0.539) | 0.404 |
| RPPA | No adjustment | 20 | 0.511 (0.491–0.539) | 0.404 |
| RPPA | No adjustment | 50 | 0.511 (0.491–0.539) | 0.404 |
| RPPA | No adjustment | 100 | 0.511 (0.491–0.539) | 0.404 |
| RPPA | Shift (δ) | 10 | 0.615 (0.540–0.671) | 0.404 |
| RPPA | Shift (δ) | 20 | 0.618 (0.542–0.664) | 0.404 |
| RPPA | Shift (δ) | 50 | 0.621 (0.579–0.659) | 0.404 |
| RPPA | Shift (δ) | 100 | 0.623 (0.595–0.659) | 0.404 |
| RPPA | Shift + quantile | 10 | 0.775 (0.633–0.924) | 0.566 |
| RPPA | Shift + quantile | 20 | 0.796 (0.693–0.890) | 0.581 |
| RPPA | Shift + quantile | 50 | 0.798 (0.714–0.879) | 0.579 |
| RPPA | Shift + quantile | 100 | 0.801 (0.734–0.853) | 0.582 |
| ProCan | No adjustment | 10 | 0.628 (0.617–0.640) | 0.404 |
| ProCan | No adjustment | 20 | 0.628 (0.617–0.640) | 0.404 |
| ProCan | No adjustment | 50 | 0.628 (0.617–0.640) | 0.404 |
| ProCan | No adjustment | 100 | 0.628 (0.617–0.640) | 0.404 |
| ProCan | Shift (δ) | 10 | 0.604 (0.501–0.646) | 0.404 |
| ProCan | Shift (δ) | 20 | 0.609 (0.548–0.649) | 0.404 |
| ProCan | Shift (δ) | 50 | 0.616 (0.580–0.649) | 0.404 |
| ProCan | Shift (δ) | 100 | 0.620 (0.592–0.647) | 0.404 |
| ProCan | Shift + quantile | 10 | 0.735 (0.587–0.889) | 0.545 |
| ProCan | Shift + quantile | 20 | 0.767 (0.581–0.892) | 0.574 |
| ProCan | Shift + quantile | 50 | 0.789 (0.701–0.850) | 0.586 |
| ProCan | Shift + quantile | 100 | 0.792 (0.745–0.851) | 0.587 |
| APOLLO | No adjustment | 10 | 0.833 (0.824–0.846) | 0.404 |
| APOLLO | No adjustment | 20 | 0.833 (0.824–0.846) | 0.404 |
| APOLLO | No adjustment | 50 | 0.833 (0.824–0.846) | 0.404 |
| APOLLO | No adjustment | 100 | 0.833 (0.824–0.846) | 0.404 |
| APOLLO | Shift (δ) | 10 | 0.816 (0.739–0.845) | 0.404 |
| APOLLO | Shift (δ) | 20 | 0.824 (0.737–0.843) | 0.404 |
| APOLLO | Shift (δ) | 50 | 0.830 (0.810–0.843) | 0.404 |
| APOLLO | Shift (δ) | 100 | 0.832 (0.819–0.844) | 0.404 |
| APOLLO | Shift + quantile | 10 | 0.719 (0.499–0.947) | 0.342 |
| APOLLO | Shift + quantile | 20 | 0.742 (0.537–0.906) | 0.342 |
| APOLLO | Shift + quantile | 50 | 0.779 (0.673–0.893) | 0.360 |
| APOLLO | Shift + quantile | 100 | 0.791 (0.722–0.861) | 0.364 |
| F1s-ridge, not the preregistered predictor | | | | |
| RPPA | No adjustment | 10 | 0.472 (0.425–0.520) | 0.404 |
| RPPA | No adjustment | 20 | 0.472 (0.425–0.520) | 0.404 |
| RPPA | No adjustment | 50 | 0.472 (0.425–0.520) | 0.404 |
| RPPA | No adjustment | 100 | 0.472 (0.425–0.520) | 0.404 |
| RPPA | Shift (δ) | 10 | 0.610 (0.526–0.695) | 0.404 |
| RPPA | Shift (δ) | 20 | 0.614 (0.530–0.689) | 0.404 |
| RPPA | Shift (δ) | 50 | 0.617 (0.530–0.692) | 0.404 |
| RPPA | Shift (δ) | 100 | 0.618 (0.532–0.689) | 0.404 |
| RPPA | Shift + quantile | 10 | 0.761 (0.583–0.904) | 0.555 |
| RPPA | Shift + quantile | 20 | 0.792 (0.640–0.904) | 0.585 |
| RPPA | Shift + quantile | 50 | 0.794 (0.670–0.876) | 0.581 |
| RPPA | Shift + quantile | 100 | 0.800 (0.686–0.863) | 0.585 |
| ProCan | No adjustment | 10 | 0.616 (0.591–0.637) | 0.404 |
| ProCan | No adjustment | 20 | 0.616 (0.591–0.637) | 0.404 |
| ProCan | No adjustment | 50 | 0.616 (0.591–0.637) | 0.404 |
| ProCan | No adjustment | 100 | 0.616 (0.591–0.637) | 0.404 |
| ProCan | Shift (δ) | 10 | 0.603 (0.479–0.645) | 0.404 |
| ProCan | Shift (δ) | 20 | 0.615 (0.557–0.641) | 0.404 |
| ProCan | Shift (δ) | 50 | 0.622 (0.581–0.645) | 0.404 |
| ProCan | Shift (δ) | 100 | 0.619 (0.583–0.641) | 0.404 |
| ProCan | Shift + quantile | 10 | 0.722 (0.551–0.904) | 0.529 |
| ProCan | Shift + quantile | 20 | 0.778 (0.661–0.878) | 0.578 |
| ProCan | Shift + quantile | 50 | 0.780 (0.682–0.858) | 0.568 |
| ProCan | Shift + quantile | 100 | 0.795 (0.738–0.846) | 0.590 |
| APOLLO | No adjustment | 10 | 0.824 (0.817–0.831) | 0.404 |
| APOLLO | No adjustment | 20 | 0.824 (0.817–0.831) | 0.404 |
| APOLLO | No adjustment | 50 | 0.824 (0.817–0.831) | 0.404 |
| APOLLO | No adjustment | 100 | 0.824 (0.817–0.831) | 0.404 |
| APOLLO | Shift (δ) | 10 | 0.807 (0.712–0.834) | 0.404 |
| APOLLO | Shift (δ) | 20 | 0.817 (0.754–0.834) | 0.404 |
| APOLLO | Shift (δ) | 50 | 0.826 (0.812–0.833) | 0.404 |
| APOLLO | Shift (δ) | 100 | 0.827 (0.814–0.835) | 0.404 |
| APOLLO | Shift + quantile | 10 | 0.705 (0.435–0.924) | 0.341 |
| APOLLO | Shift + quantile | 20 | 0.765 (0.551–0.913) | 0.369 |
| APOLLO | Shift + quantile | 50 | 0.775 (0.664–0.857) | 0.358 |
| APOLLO | Shift + quantile | 100 | 0.789 (0.725–0.861) | 0.367 |

**Supplementary Table S7. openFDA and Drugs@FDA search for antibody–drug conjugates (search date 15 October 2026; approval cutoff 2 October 2026).**

Queries and the number of records returned. A count of 0 is an empty response, not a failed request. The active-ingredient queries were added after `openfda.generic_name:mafodotin` returned 0 and missed Blenrep; the belantamab label query was added because application BLA761440 has no openFDA label record; the sunirine query was added for a phrase hit whose payload stem was not in the planned list.

| Endpoint | Query | Hits | Note |
|---|---|---|---|
| label | `"antibody-drug conjugate"` | 20 | label text contains the phrase |
| drugsfda | `openfda.generic_name:vedotin` | 5 | planned generic-name stem query |
| drugsfda | `openfda.generic_name:mafodotin` | 0 | planned generic-name stem query |
| drugsfda | `openfda.generic_name:deruxtecan` | 2 | planned generic-name stem query |
| drugsfda | `openfda.generic_name:govitecan` | 1 | planned generic-name stem query |
| drugsfda | `openfda.generic_name:emtansine` | 1 | planned generic-name stem query |
| drugsfda | `openfda.generic_name:soravtansine` | 1 | planned generic-name stem query |
| drugsfda | `openfda.generic_name:tesirine` | 1 | planned generic-name stem query |
| drugsfda | `openfda.generic_name:ozogamicin` | 2 | planned generic-name stem query |
| drugsfda | `openfda.generic_name:duocarmazine` | 0 | planned generic-name stem query |
| drugsfda | `openfda.generic_name:tirumotecan` | 0 | planned generic-name stem query |
| drugsfda | `products.active_ingredients.name:vedotin` | 6 | added after openfda.generic_name:mafodotin returned 0 and missed Blenrep |
| drugsfda | `products.active_ingredients.name:mafodotin` | 1 | added after openfda.generic_name:mafodotin returned 0 and missed Blenrep |
| drugsfda | `products.active_ingredients.name:deruxtecan` | 3 | added after openfda.generic_name:mafodotin returned 0 and missed Blenrep |
| drugsfda | `products.active_ingredients.name:govitecan` | 1 | added after openfda.generic_name:mafodotin returned 0 and missed Blenrep |
| drugsfda | `products.active_ingredients.name:emtansine` | 1 | added after openfda.generic_name:mafodotin returned 0 and missed Blenrep |
| drugsfda | `products.active_ingredients.name:soravtansine` | 1 | added after openfda.generic_name:mafodotin returned 0 and missed Blenrep |
| drugsfda | `products.active_ingredients.name:tesirine` | 1 | added after openfda.generic_name:mafodotin returned 0 and missed Blenrep |
| drugsfda | `products.active_ingredients.name:ozogamicin` | 2 | added after openfda.generic_name:mafodotin returned 0 and missed Blenrep |
| drugsfda | `products.active_ingredients.name:duocarmazine` | 0 | added after openfda.generic_name:mafodotin returned 0 and missed Blenrep |
| drugsfda | `products.active_ingredients.name:tirumotecan` | 0 | added after openfda.generic_name:mafodotin returned 0 and missed Blenrep |
| label | `"belantamab"` | 2 | added because the Blenrep application number has no openFDA label record |
| drugsfda | `products.active_ingredients.name:sunirine` | 1 | phrase hit whose payload stem is not in the planned list |

Products whose active-ingredient name contained a payload stem, plus the phrase hit whose mechanism says the product is an antibody–drug conjugate. The HGNC symbol is the symbol whose name appears in the quoted mechanism text. Two applications of Adcetris and two of Datroway are listed separately. Of the 29 atlas targets, this search found an approved product for ERBB2, TACSTD2, NECTIN4, F3, FOLR1 and MET (eight products). For the other 23 targets, no product in this search acted on the gene.

| Brand | Generic name | Application | First approval | Label revision | Mechanism of action (quoted) | HGNC symbol |
|---|---|---|---|---|---|---|
| ADCETRIS | BRENTUXIMAB VEDOTIN | BLA125388 | 2011-08-19 | 2025-11-11 | CD30 is a member of the tumor necrosis factor receptor family and is expressed on the surface of sALCL cells and on Hodgkin Reed-Sternberg (HRS) cells in cHL. Brentuximab vedotin is an antibody-drug conjugate (ADC). | TNFRSF8 |
| ADCETRIS | BRENTUXIMAB VEDOTIN | BLA125399 | 2011-08-19 | 2025-11-11 | CD30 is a member of the tumor necrosis factor receptor family and is expressed on the surface of sALCL cells and on Hodgkin Reed-Sternberg (HRS) cells in cHL. Brentuximab vedotin is an antibody-drug conjugate (ADC). no label record for BLA125399; mechanism quoted from BLA125388, the other application of the same brand | TNFRSF8 |
| KADCYLA | ADO-TRASTUZUMAB EMTANSINE | BLA125427 | 2013-02-22 | 2026-04-01 | Ado-trastuzumab emtansine is a HER2-targeted antibody-drug conjugate. | ERBB2 |
| BESPONSA | INOTUZUMAB OZOGAMICIN | BLA761040 | 2017-08-17 | 2026-02-06 | Inotuzumab ozogamicin is a CD22-directed antibody drug conjugate (ADC). | CD22 |
| MYLOTARG | GEMTUZUMAB OZOGAMICIN | BLA761060 | 2017-09-01 | 2025-12-19 | Gemtuzumab ozogamicin is a CD33-directed antibody-drug conjugate (ADC). | CD33 |
| TRODELVY | SACITUZUMAB GOVITECAN-HZIY | BLA761115 | 2020-04-22 | 2026-06-24 | Sacituzumab govitecan-hziy is a Trop-2-directed antibody-drug conjugate. | TACSTD2 |
| POLIVY | POLATUZUMAB VEDOTIN-PIIQ | BLA761121 | 2019-06-10 | 2026-08-24 | Polatuzumab vedotin-piiq is a CD79b-directed antibody-drug conjugate with activity against dividing B cells. | CD79B |
| PADCEV | ENFORTUMAB VEDOTIN-EJFV | BLA761137 | 2019-12-18 | 2026-08-04 | The antibody is a human IgG1 kappa directed against Nectin-4, an adhesion protein located on the surface of cells. Enfortumab vedotin-ejfv is an ADC. | NECTIN4 |
| ENHERTU | FAM-TRASTUZUMAB DERUXTECAN-NXKI | BLA761139 | 2019-12-20 | 2026-05-15 | Fam-trastuzumab deruxtecan-nxki is a HER2-directed antibody-drug conjugate. | ERBB2 |
| ZYNLONTA | LONCASTUXIMAB TESIRINE-LPYL | BLA761196 | 2021-04-23 | 2026-06-12 | Loncastuximab tesirine-lpyl is an antibody-drug conjugate (ADC) targeting CD19. | CD19 |
| TIVDAK | TISOTUMAB VEDOTIN-TFTV | BLA761208 | 2021-09-20 | 2025-11-20 | Tisotumab vedotin-tftv is a tissue factor (TF)-directed antibody drug conjugate (ADC). | F3 |
| ELAHERE | MIRVETUXIMAB SORAVTANSINE-GYNX | BLA761310 | 2022-11-14 | 2025-07-14 | The antibody is a chimeric IgG1 directed against folate receptor alpha (FRα). Mirvetuximab soravtansine-gynx is an antibody-drug conjugate (ADC). | FOLR1 |
| EMRELIS | TELISOTUZUMAB VEDOTIN-TLLV | BLA761384 | 2025-05-14 | 2025-05-25 | Telisotuzumab vedotin-tllv is a c-Met-directed antibody drug conjugate (ADC). | MET |
| DATROWAY | DATOPOTAMAB DERUXTECAN-DLNK | BLA761394 | 2025-01-17 | 2026-05-27 | Datopotamab deruxtecan-dlnk, is a Trop-2-directed antibody-drug conjugate. | TACSTD2 |
| BLENREP | BELANTAMAB MAFODOTIN-BLMF | BLA761440 | 2025-10-23 | 2025-10-23 | The antibody component is an afucosylated IgG1 directed against BCMA, a protein expressed on normal B lymphocytes and multiple myeloma cells. Belantamab mafodotin‑blmf is an antibody‑drug conjugate (ADC). | TNFRSF17 |
| DECNUPAZ | PIVEKIMAB SUNIRINE-PVZY | BLA761460 | 2026-05-27 | 2026-05-27 | Pivekimab sunirine-pvzy is a CD123 (alpha-subunit of the interleukin-3 receptor)-directed antibody-drug conjugate (ADC). | IL3RA |
| DATROWAY | DATOPOTAMAB DERUXTECAN-DLNK | BLA761464 | 2025-06-23 | 2026-05-27 | Datopotamab deruxtecan-dlnk, is a Trop-2-directed antibody-drug conjugate. no label record for BLA761464; mechanism quoted from BLA761394, the other application of the same brand | TACSTD2 |

**Supplementary Table S8. Pairs in each level of the cascade lookup.**

| Level | RPPA pairs (antibodies) | ProCan pairs (proteins) | APOLLO genes |
|---|---|---|---|
| L1, matched cancer | 1,206 (127) | 3,356 (2,091) | 6,982 |
| L2, CPTAC mean | 2,653 (127) | 2,616 (2,570) | 93 |
| L3, ProCan mean | 62 (2) | not used | 26 |
| L4, global mean | 124 (4) | 146 (79) | 95 |
| Total | 4,045 (134) | 6,118 (2,746) | 7,196 |



**Supplementary Table S9. Downloaded files.**

| File | URL | Bytes | SHA-256 | Downloaded |
|---|---|---|---|---|
| `data/raw/cptac/BRCA/BRCA_meta.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/BRCA/BRCA_meta.txt | 7369 | `16d00a3022fbafa9d29e74c0ee689ae63279911f04aeef2d21344ed51a0a5cb7` | 2026-10-01T06:03:52Z |
| `data/raw/cptac/BRCA/BRCA_Tumor_CaseList.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/BRCA/BRCA_Tumor_CaseList.txt | 1094 | `27c5bbd3b25f2421db1fcc174a1d0cde9835a191eb8b8d957241ef0b531b9993` | 2026-10-01T06:03:53Z |
| `data/raw/cptac/BRCA/BRCA_phenotype.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/BRCA/BRCA_phenotype.txt | 227024 | `7f9fa329e1352695ea898de6caecde43dfb45445c63911d9b649747517dd1091` | 2026-10-01T06:03:54Z |
| `data/raw/cptac/BRCA/BRCA_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/BRCA/BRCA_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt | 33078031 | `1d7a10c7ed2df611ec807e08e0f8e0694689a8582f9c3a98fbc34235b002d819` | 2026-10-01T06:04:12Z |
| `data/raw/cptac/BRCA/BRCA_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/BRCA/BRCA_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt | 21398351 | `e584331bc77a440af50692c7874b736b3d135d01d96ce95d66af18cc64009b21` | 2026-10-01T06:04:21Z |
| `data/raw/cptac/CCRCC/CCRCC_meta.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/CCRCC/CCRCC_meta.txt | 10979 | `c332231f3cbd204094236dfd302977e70bbdddcf0049e862b60dd10244d04dae` | 2026-10-01T06:04:22Z |
| `data/raw/cptac/CCRCC/CCRCC_Tumor_CaseList.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/CCRCC/CCRCC_Tumor_CaseList.txt | 1030 | `ee692071c563ae797e5d7670e2b8698049e2467eb30f5916f7707e798d1441be` | 2026-10-01T06:04:23Z |
| `data/raw/cptac/CCRCC/CCRCC_Normal_CaseList.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/CCRCC/CCRCC_Normal_CaseList.txt | 800 | `4e631863d653d09b11e249176f3318c96928240f02f631725b2053be7d92f6ea` | 2026-10-01T06:04:24Z |
| `data/raw/cptac/CCRCC/CCRCC_phenotype.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/CCRCC/CCRCC_phenotype.txt | 191249 | `929ce841da91caf2a03a58b1561929648f4be50b8b081f7182e2ad95ee936d0a` | 2026-10-01T06:04:25Z |
| `data/raw/cptac/CCRCC/CCRCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/CCRCC/CCRCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt | 33113844 | `bcef3edfb1fcc738859369d78f2b7edb3f0972d438b4b1f50a114519f2ab9a50` | 2026-10-01T06:05:17Z |
| `data/raw/cptac/CCRCC/CCRCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/CCRCC/CCRCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt | 23810094 | `e4505de10482e5738937ebe7c9b4239ffdc86e51cb092db7f2473e92040aa727` | 2026-10-01T06:06:15Z |
| `data/raw/cptac/CCRCC/CCRCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/CCRCC/CCRCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt | 16629110 | `090077c25f0eeb59b737912f39dec98901db949e1f143bbc481d5cd731550f60` | 2026-10-01T06:06:20Z |
| `data/raw/cptac/CCRCC/CCRCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/CCRCC/CCRCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt | 13026061 | `1b00d684cd99d870bee199128485c22d6cb079e0a01de3c044b7a4fdb1735618` | 2026-10-01T06:06:25Z |
| `data/raw/cptac/COAD/COAD_meta.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/COAD/COAD_meta.txt | 11672 | `5242103f57fb5d919ab7e228dd36c97b669915d6bee88ae1e5a43c3aa270295e` | 2026-10-01T06:06:26Z |
| `data/raw/cptac/COAD/COAD_Tumor_CaseList.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/COAD/COAD_Tumor_CaseList.txt | 988 | `1495c045402199d98670c8862ac01035b5521b3b8538b38c117fb80340009ac8` | 2026-10-01T06:06:27Z |
| `data/raw/cptac/COAD/COAD_Normal_CaseList.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/COAD/COAD_Normal_CaseList.txt | 898 | `7e43c6619464ae23a7c1162e9266f0c71175946fa290b37e0405a2e28a4f390c` | 2026-10-01T06:06:28Z |
| `data/raw/cptac/COAD/COAD_phenotype.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/COAD/COAD_phenotype.txt | 204627 | `959ec27c2b668bfd2b0cd9ff811319efe62d3078c37c6ca590aa8be8b3b225d9` | 2026-10-01T06:06:29Z |
| `data/raw/cptac/COAD/COAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/COAD/COAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt | 28863817 | `e7c5108dcea2aa615b0c1288d7a8b43bdc1c2469f1a9543a4078f7ee048ef253` | 2026-10-01T06:06:37Z |
| `data/raw/cptac/COAD/COAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/COAD/COAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt | 12293004 | `fa2dde5f9165acdf6320d28e01dbad24b735ede320f5fab4a6ceab6c15aa0580` | 2026-10-01T06:06:42Z |
| `data/raw/cptac/COAD/COAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/COAD/COAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt | 12515956 | `a76a0cc0ab3986635e8745028eb0bcc64470035ab737cc0937ec7be114a09177` | 2026-10-01T06:06:47Z |
| `data/raw/cptac/GBM/GBM_meta.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/GBM/GBM_meta.txt | 6809 | `6617ef2f9010acec1eba09a765ca3e35e4ba2ea7ba75ba9e8734ed9ad94feecc` | 2026-10-01T06:06:48Z |
| `data/raw/cptac/GBM/GBM_Tumor_CaseList.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/GBM/GBM_Tumor_CaseList.txt | 1087 | `51dc5cda8b5974690e9065751c5c3437b04ea210dfc895bcaee2e960e40ff2e1` | 2026-10-01T06:06:48Z |
| `data/raw/cptac/GBM/GBM_phenotype.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/GBM/GBM_phenotype.txt | 209875 | `807c6b83236dfb96c10b8ad224f68b21fb356bca4bd14514095cb1c5256ea575` | 2026-10-01T06:06:50Z |
| `data/raw/cptac/GBM/GBM_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/GBM/GBM_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt | 32926067 | `2c3c0e0e8f7f5466a3d078cf868e1f7f5e3529b5d7ee04eabec97e9ade9a1ea1` | 2026-10-01T06:06:58Z |
| `data/raw/cptac/GBM/GBM_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/GBM/GBM_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt | 18033471 | `d7393915eb8ec10c1c48bea03d78cd9ee2a3477f7bb4b93f808ffade03c2f01b` | 2026-10-01T06:07:06Z |
| `data/raw/cptac/HNSCC/HNSCC_meta.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/HNSCC/HNSCC_meta.txt | 12087 | `3c9e32ff8b51b4a22f967ca753185e01d79f69d860f9530f5587314bd06c80ac` | 2026-10-01T06:07:07Z |
| `data/raw/cptac/HNSCC/HNSCC_Tumor_CaseList.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/HNSCC/HNSCC_Tumor_CaseList.txt | 1080 | `18037cf6f7322421d06094aab41cff95ad1ea0d3cbf8298013a1d315feb3c4e4` | 2026-10-01T06:07:08Z |
| `data/raw/cptac/HNSCC/HNSCC_Normal_CaseList.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/HNSCC/HNSCC_Normal_CaseList.txt | 620 | `7c943da5e8f6152a399433a2f166294a2de327d7d0e6693dd75896747f32cec1` | 2026-10-01T06:07:09Z |
| `data/raw/cptac/HNSCC/HNSCC_phenotype.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/HNSCC/HNSCC_phenotype.txt | 214773 | `579f51424bfbff5f9baf17553c95c12dfd94fc4a51b2f3562768d7c255adbd10` | 2026-10-01T06:07:10Z |
| `data/raw/cptac/HNSCC/HNSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/HNSCC/HNSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt | 35019254 | `a308c4410eb347e9036dff918935804981f5c67d7a2f4f9003c47d41de8cf195` | 2026-10-01T06:07:20Z |
| `data/raw/cptac/HNSCC/HNSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/HNSCC/HNSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt | 15523255 | `bb190f7493f83d693e8232f43ebf6a3e09b12ea23f84528fa04f66e56d2638db` | 2026-10-01T06:07:25Z |
| `data/raw/cptac/HNSCC/HNSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/HNSCC/HNSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt | 18184418 | `6d52dac50c55a4aa30e0aabaf3047a3a099875ea075b65d8f0bc83ac45a5c4b3` | 2026-10-01T06:07:31Z |
| `data/raw/cptac/HNSCC/HNSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/HNSCC/HNSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt | 10542735 | `0a067886f9b553d3de24d2dbab21e18037c44a50a52b9d60004d22dba3442adc` | 2026-10-01T06:07:36Z |
| `data/raw/cptac/LSCC/LSCC_meta.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LSCC/LSCC_meta.txt | 12051 | `6c334bbd8c9c4acfc94709a7fea5fd0181efa0881cee8ce1793085a5281a9d99` | 2026-10-01T06:07:37Z |
| `data/raw/cptac/LSCC/LSCC_Tumor_CaseList.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LSCC/LSCC_Tumor_CaseList.txt | 1186 | `8198fcc0e10e6bfaaabbb9886437b8e805af11e8f581fbff5ba63ec562c5c0eb` | 2026-10-01T06:07:38Z |
| `data/raw/cptac/LSCC/LSCC_Normal_CaseList.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LSCC/LSCC_Normal_CaseList.txt | 1087 | `5fdf25da74eba7a2f1cfeb3b7424b1cf6cc45cfc8163d721785abf908f424a27` | 2026-10-01T06:07:39Z |
| `data/raw/cptac/LSCC/LSCC_phenotype.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LSCC/LSCC_phenotype.txt | 212329 | `34f17a653dd5ebce0cdf1257c4498fe6579e314d18d435d333e4bcc4db22ef8d` | 2026-10-01T06:07:40Z |
| `data/raw/cptac/LSCC/LSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LSCC/LSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt | 35179517 | `0da58ae002ccc65cc83f13038f7f8439ad672898822e545e475f99540f5b7377` | 2026-10-01T06:07:48Z |
| `data/raw/cptac/LSCC/LSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LSCC/LSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt | 30624657 | `61822d27dda6ef9cf67743e8b9511f84620042d1e703fac8c5b25cfdc2a27f07` | 2026-10-01T06:07:54Z |
| `data/raw/cptac/LSCC/LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LSCC/LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt | 20631271 | `fd9f9cc0aed40942fd33d750d1d207bcdf9ea231fcddeff3cc8ac8fdf3e1384d` | 2026-10-01T06:08:00Z |
| `data/raw/cptac/LSCC/LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LSCC/LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt | 18920447 | `ea3f8e5a1e42f708af87f6cdfc86d3fa63a96bfdc4f6f3316ba54380ea579aaf` | 2026-10-01T06:08:05Z |
| `data/raw/cptac/LUAD/LUAD_meta.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LUAD/LUAD_meta.txt | 12232 | `f9facf713af130ad097c56794f1fbc9107334e02ad9bf0df3958b40d68b6af10` | 2026-10-01T06:08:06Z |
| `data/raw/cptac/LUAD/LUAD_Tumor_CaseList.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LUAD/LUAD_Tumor_CaseList.txt | 1092 | `92aba1a2d482feccb2f3f1873e147758674ee7472519340d3859e0c43f983ddf` | 2026-10-01T06:08:07Z |
| `data/raw/cptac/LUAD/LUAD_Normal_CaseList.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LUAD/LUAD_Normal_CaseList.txt | 1010 | `c6425a338baa30cdf33fcb617946ba6693d90100e729bbcc4ece7e9eea0ce623` | 2026-10-01T06:08:08Z |
| `data/raw/cptac/LUAD/LUAD_phenotype.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LUAD/LUAD_phenotype.txt | 210137 | `9f777402f7d7d53856a5adfc73c8291136430ddfad85641584426921609b349a` | 2026-10-01T06:08:09Z |
| `data/raw/cptac/LUAD/LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LUAD/LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt | 35757367 | `c2c5e97fd469591d7aec13bd28ba68aa19fe7ced161592f69924ecc388513e05` | 2026-10-01T06:08:19Z |
| `data/raw/cptac/LUAD/LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LUAD/LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt | 32856789 | `a54d7023a79ccd5159c01f70dbc6cec472e02b0ed9142d052c816cee51220b70` | 2026-10-01T06:08:28Z |
| `data/raw/cptac/LUAD/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LUAD/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt | 20132337 | `746900f54072198481207ea13f8a7ba14dcabd52345bace0bef7438120c03845` | 2026-10-01T06:08:36Z |
| `data/raw/cptac/LUAD/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LUAD/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt | 18491662 | `3cfc62b159276837dea644f457ebf6b212964320e0b3d8ac87168d82eb73747c` | 2026-10-01T06:08:43Z |
| `data/raw/cptac/OV/OV_meta.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/OV/OV_meta.txt | 6029 | `ab6b3bb70460c682492a98f7620cbbeda61baab6e4d8bc1e25d19c0580390123` | 2026-10-01T06:08:44Z |
| `data/raw/cptac/OV/OV_Tumor_CaseList.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/OV/OV_Tumor_CaseList.txt | 664 | `2c122ca238e35876f1166d7edef43ab87bd03b1c6baadcd97b1d9a9715d46312` | 2026-10-01T06:08:45Z |
| `data/raw/cptac/OV/OV_Normal_CaseList.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/OV/OV_Normal_CaseList.txt | 160 | `318a82c237b00e8f47b86a6cd90eae685aff01b8d02a83b49b355608e1c09d6e` | 2026-10-01T06:08:46Z |
| `data/raw/cptac/OV/OV_phenotype.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/OV/OV_phenotype.txt | 157772 | `58c75e054cb98a551e6a06e21784eb70f428da14779ace530c5549bef9352e8b` | 2026-10-01T06:08:47Z |
| `data/raw/cptac/OV/OV_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/OV/OV_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt | 22239335 | `5b25a3cbea0bb041f40fccb8b4120cf2b8135f16f102aacb1e4cb39923509d97` | 2026-10-01T06:08:53Z |
| `data/raw/cptac/OV/OV_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/OV/OV_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt | 12506387 | `274c22512ab70e2a07197ca6e51d268038c90e903bac823c1b27b7950af4ae6f` | 2026-10-01T06:08:57Z |
| `data/raw/cptac/OV/OV_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/OV/OV_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt | 3007610 | `6fc249b27c9284b598a94860c1e36f4fd45843f65b14c4ccdfbc0153725bb07a` | 2026-10-01T06:09:00Z |
| `data/raw/cptac/PDAC/PDAC_meta.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/PDAC/PDAC_meta.txt | 11358 | `ce53f0acf92cbbbf165c8c27425c27b7b6de301c2e95279eaf124c37a2316d33` | 2026-10-01T06:09:01Z |
| `data/raw/cptac/PDAC/PDAC_Tumor_CaseList.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/PDAC/PDAC_Tumor_CaseList.txt | 1155 | `089f7a166f4b90e888c367aa1582d6bc675b572c7207ee4d7670daa325b1e36b` | 2026-10-01T06:09:02Z |
| `data/raw/cptac/PDAC/PDAC_Normal_CaseList.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/PDAC/PDAC_Normal_CaseList.txt | 484 | `9f1dcb19bfb6a4630901edbc7b4e1ad77d70fb2d47eb6c8459620a524a3dd9c4` | 2026-10-01T06:09:03Z |
| `data/raw/cptac/PDAC/PDAC_phenotype.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/PDAC/PDAC_phenotype.txt | 206286 | `91c4511765e6c99e3e9b1a3711f31d313a1dee408072dafe375fcfc6f1af4a52` | 2026-10-01T06:09:04Z |
| `data/raw/cptac/PDAC/PDAC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/PDAC/PDAC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt | 34095120 | `3f00449e7112f8d15478ceb0563df1c0ad633823756638fdca66bef4986b4bad` | 2026-10-01T06:09:58Z |
| `data/raw/cptac/PDAC/PDAC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/PDAC/PDAC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt | 5852445 | `14c9e9b708763678a9fed459a81dc017564aa38afad9e667ab714f58594a8774` | 2026-10-01T06:10:19Z |
| `data/raw/cptac/PDAC/PDAC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/PDAC/PDAC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt | 16800925 | `e1be4bed67bfdac6fff7b12487b8af9e0e2a0e93901039f5f9124e180f24f6bd` | 2026-10-01T06:11:42Z |
| `data/raw/cptac/PDAC/PDAC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/PDAC/PDAC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt | 7177097 | `beb1bc33ee52d10f0f7a6fa0bfc0de97916471c37d80c5cb655ee0d0a2452ce9` | 2026-10-01T06:12:05Z |
| `data/raw/cptac/UCEC/UCEC_meta.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/UCEC/UCEC_meta.txt | 13972 | `4058a8420129dd37f4d9035e1d87beb70a800ce9cb42754f74e7ad3d35f3020b` | 2026-10-01T06:12:06Z |
| `data/raw/cptac/UCEC/UCEC_Tumor_CaseList.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/UCEC/UCEC_Tumor_CaseList.txt | 949 | `9f8db0ff8d82a54f9fc962f53ca4300107d08bc644b734912eccc3147a647545` | 2026-10-01T06:12:07Z |
| `data/raw/cptac/UCEC/UCEC_Normal_CaseList.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/UCEC/UCEC_Normal_CaseList.txt | 97 | `79ee61ccfa92e1e8bb1f17b94a7176a062aa40e1c0acb46294af99c7a84fce9e` | 2026-10-01T06:12:08Z |
| `data/raw/cptac/UCEC/UCEC_phenotype.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/UCEC/UCEC_phenotype.txt | 197288 | `89ddc8c4953184d32a592cf0240ba97658fc2bb78e7caa37ab90dd74dd5e87ad` | 2026-10-01T06:12:10Z |
| `data/raw/cptac/UCEC/UCEC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/UCEC/UCEC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt | 30753085 | `3d708a346477a51f2aa1fc80479cd60b1888644145b69b1082c31c6405770882` | 2026-10-01T06:14:39Z |
| `data/raw/cptac/UCEC/UCEC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/UCEC/UCEC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt | 16560392 | `aea5d956b3496530faf86052886de19e857e96294a2dbcb3adb9702457865d36` | 2026-10-01T06:15:52Z |
| `data/raw/cptac/UCEC/UCEC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt` | https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/UCEC/UCEC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt | 3354332 | `225df3e7827fbeec499deba3d773828935dbe722559730a902c66f9dabe1d229` | 2026-10-01T06:16:07Z |
| `data/raw/hgnc/hgnc_complete_set.txt` | https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt | 16953545 | `91d0ad20c34f26fb12c8cce3f80a65a6fc694f3a508040b19f6bfa29ae0d889e` | 2026-10-01T06:16:19Z |
| `data/raw/ccle/protein_quant_current_normalized.csv.gz` | https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz | 30916232 | `b72a9ff39b123c801249ea212844ac1f500c3f38fa0f399ec6ad08ee1be2a7ed` | 2026-10-01T06:20:10Z |
| `data/raw/depmap/download_files_nocaptcha.csv` | https://depmap.org/portal/api/no-captcha/download/files | 163720 | `efe6f26a1109880c37e8bc98c1b24d4bf217dd4b3520d78a7b2a15efae24f538` | 2026-10-01T06:20:12Z |
| `data/raw/depmap/Model.csv` | https://ndownloader.figshare.com/files/51065297 | 645696 | `b7a0c1385e6cef30132b56aff61f1261d11e3f490490b355c430d32ee0dbdcfa` | 2026-10-01T06:20:18Z |
| `data/raw/depmap/OmicsExpressionProteinCodingGenesTPMLogp1.csv` | https://ndownloader.figshare.com/files/51065489 | 506628654 | `2a71dc94110efcc0221eae821bb93a9f03b54bea16f005818911a09d33383d56` | 2026-10-01T07:36:23Z |
| `data/raw/procan/ProCan-DepMapSanger_protein_matrix_8498_averaged.txt` | https://ndownloader.figshare.com/files/34411172 | 94222540 | `5e7adf7e8217e5a02389c299dce1266e74c7a3b786b8c177e56d12402fe783ac` | 2026-10-01T07:36:23Z |
| `data/raw/procan/ProCan-DepMapSanger_mapping_file_averaged.txt` | https://ndownloader.figshare.com/files/34411133 | 82593 | `68ff17c98d3f163992418f66e4ea776d484eb11475b4a8e1d079767eaa703b7e` | 2026-10-01T07:36:23Z |
| `data/raw/procan/figshare_README.md` | https://ndownloader.figshare.com/files/36155439 | 3218 | `aaa5ca2539d67f12dcbbd893faff523e328cc69a994e81c2eeff8de0d93102ca` | 2026-10-01T07:36:23Z |
| `data/raw/depmap/README_24Q4.txt` | https://ndownloader.figshare.com/files/51065795 | 43103 | `f744a5cf5c112e01276e772679157aac8e85cdc4a8f6d3b86d808e3bfd332891` | 2026-10-01T07:36:23Z |
| `data/raw/tcga/TCGA-RPPA-pancan-clean.xena.gz` | https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/TCGA-RPPA-pancan-clean.xena.gz | 5602638 | `2dca46e09280d2c26f9f83e82ee98c8c41d1b896f2966918529e913bfad1868d` | 2026-10-01T07:36:23Z |
| `data/raw/tcga/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz` | https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz | 331000731 | `a00532ec86af8c07630c618f10f6277f09c484d0a9c17db5901edf95c7714b38` | 2026-10-01T07:36:24Z |
| `data/raw/tcga/TCGA_phenotype_denseDataOnlyDownload.tsv.gz` | https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz | 61165 | `ac3923bb155ada3391b653a811f6c8c9c4734593b0c2b0a1555728928d83e5a2` | 2026-10-01T07:36:24Z |
| `data/raw/hpa/normal_ihc_data.tsv.zip` | https://www.proteinatlas.org/download/tsv/normal_ihc_data.tsv.zip | 5732831 | `fa2315e5a53ecffd9427a5bfe219c8cc41b54089ef6f6be7b5d79b5455e0677f` | 2026-10-01T07:36:24Z |
| `data/raw/hpa/cancer_data.tsv.zip` | https://www.proteinatlas.org/download/tsv/cancer_data.tsv.zip | 1724018 | `6e771693dc1019a3536d8834175ee97dcfde4f835d5cb3bed7c5a48803668135` | 2026-10-01T07:36:24Z |
| `data/raw/cmp/model_list_20260921.csv` | https://cog.sanger.ac.uk/cmp/download/model_list_20260921.csv | 957555 | `419f65d1d6534454194f1a744550fab632cd1d7411f6b13411a885735aa37006` | 2026-10-01T07:36:24Z |
| `data/raw/docs/sanger_expression.html` | https://depmap.sanger.ac.uk/documentation/datasets/expression/ | 30458 | `ff89f0a5bb3bef4e588e6cac4a1433229a807c236390bdc4eb12d0c7012686bc` | 2026-10-01T07:36:24Z |
| `data/raw/docs/sanger_proteomics.html` | https://depmap.sanger.ac.uk/documentation/datasets/proteomics/ | 27484 | `a55e5ba6405f4cd92b61d5abc6252c49b38bb68fc4a85096d65b9f8e781c1c21` | 2026-10-01T07:36:24Z |
| `data/raw/cmp/rnaseq_merged_20260323.zip` | https://cog.sanger.ac.uk/cmp/download/rnaseq_merged_20260323.zip | 1291769946 | `d2a6acd70647cc4098a990cdedaca564c74411a87385b6f0d33e8081ce27b4c9` | 2026-10-01T08:03:27Z |
| `data/raw/rppa_map/acc_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/acc_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/blca_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/blca_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/brca_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/brca_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/cesc_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/cesc_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/chol_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/chol_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/coadread_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/coadread_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/dlbc_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/dlbc_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/esca_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/esca_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/gbm_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/gbm_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/hnsc_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/hnsc_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/kich_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/kich_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/kirc_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/kirc_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/kirp_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/kirp_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/lgg_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/lgg_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/lihc_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/lihc_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/luad_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/luad_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/lusc_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/lusc_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/meso_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/meso_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/ov_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/ov_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/paad_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/paad_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/pcpg_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/pcpg_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/prad_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/prad_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/sarc_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/sarc_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/skcm_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/skcm_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/stad_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/stad_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/tgct_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/tgct_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/thca_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/thca_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/thym_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/thym_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/ucec_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/ucec_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/ucs_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/ucs_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |
| `data/raw/rppa_map/uvm_tcga_pan_can_atlas_2018.antibody_col1.txt` | https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/uvm_tcga_pan_can_atlas_2018/data_rppa.txt | 2953 | `5160f8582dc0e8d549ec1fce53c94b9e075bc5f8e0588114c309c19bab905773` | 2026-10-01T11:08:10Z |

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



## Supplementary figures

**Supplementary Figure S1. Forecast minus CPTAC mean lookup by membership in the P gene set (post hoc).** Differences in Spearman correlation with 95% intervals for F1s-LightGBM and F1s-ridge on all pairs, pairs inside P and pairs outside P, for RPPA, ProCan and APOLLO. Pale bands and horizontal lines separate the three platforms. Values are in Supplementary Table S2.

**Supplementary Figure S2. Within-cancer-type comparison on RPPA (post hoc).** Spearman correlation of the forecast against that of the cascade lookup with observed concordance in each of the 31 RPPA cancer types; the line marks equality. Cancer types below the line are labelled. 

## Supplementary references

[Converted from citation keys; numbered separately from the main text, or merged into the main reference list if the journal prefers. CHECK the Scientific Reports rule and follow it.]
