<!--
B_main_srep_claude.md — Scientific Reports manuscript written by Claude for B 지시서 10.
How to use (agent):
  [FILL: ...]   value to fill from the named file. Never guess.
  [CHECK: ...]  statement to verify against code, config or result files. Fix the text if the file says otherwise and list the change in the report.
  [@Key]        citation key. The key table is in B 지시서 10 §6. Convert to numbered Scientific Reports references in order of first citation.
  <!-- src -->  source of the numbers in the preceding sentence, for the ledger.
  [OPTION A/B]  choose by the rule stated in B 지시서 10.
All tags, comments and keys must be gone from the submitted files.
-->

# A preregistered model forecasts mRNA–protein concordance in unmeasured cancer contexts across three proteomic platforms

<!-- Title: 15 words. Alternative for the author: "Forecasting how far mRNA can stand in for protein in tumor cohorts without proteomes: a preregistered three-platform evaluation" -->

[AUTHOR TO FILL: authors, affiliations, corresponding author and email]

## Abstract

Many tumor cohorts have transcriptomes but no matched proteomes, so mRNA is used as a stand-in for protein. Whether this is reasonable depends on each gene's mRNA–protein concordance in each context, known only where both were measured. We trained a model that forecasts this concordance from gene-level information and the gene's RNA distribution in the query cohort, using 1,027 CPTAC tumors of 10 cancer types, and preregistered tests on three external platforms. The forecast ranked observed concordance better than an RNA-variance baseline on reverse-phase protein arrays of TCGA tumors (Spearman correlation 0.43 versus 0.23), on DIA proteomes of cancer cell lines (0.42 versus 0.14) and on a TMT lung adenocarcinoma cohort (0.75 versus 0.49). In post hoc comparisons, it was at least as accurate as looking up observed CPTAC correlations, even from the matched cancer type. Adjusting prediction intervals with 50 locally measured genes brought 80% interval coverage to 0.78–0.80 on all three platforms. An anchor-residual protein predictor tested alongside improved on the gene's own mRNA in both tumor cohorts but not in cell lines. [OPTION A: We provide forecasts for [FILL: number of genes] genes across 33 cancer types, including 29 antibody–drug conjugate targets.] [OPTION B: We provide forecasts for 29 antibody–drug conjugate targets across 33 cancer types.]
<!-- src: results/cptac_sample_counts.tsv; results/stage4/confirm/hypotheses.tsv; results/stage6/confirm/hypotheses.tsv; results/stage8/spearman_diff.tsv; results/stage7/calibration_table.tsv (K = 50: 0.798, 0.789, 0.779); results/stage7/adc_atlas_final.tsv. Word limit 200; check after OPTION choice. -->

## Introduction

Large tumor cohorts such as The Cancer Genome Atlas (TCGA) profiled RNA in thousands of tumors[@Hoadley2018], whereas matched proteomes exist for far fewer[@LiY2023]. As a result, mRNA is routinely used in place of protein when a question is really about protein: whether a candidate drug target is expressed, which tumors carry more of it, or which cancer types to prioritize for an antibody–drug conjugate (ADC)[@Drago2021; @Bosi2023; @Savage2024].

This substitution is only as good as the concordance between a gene's mRNA and its protein across tumors. Concordance varies widely between genes[@Vogel2012; @Liu2016], differs between tumor types[@Savage2024] [CHECK: keep this citation only if the paper reports mRNA–protein correlations that differ between tumor types; otherwise use LiY2023 if it does, or delete "differs between tumor types"], is bounded by the reproducibility of both measurements[@Upadhya2022] and is lowered further when protein is measured with antibodies of uneven quality[@Upadhya2023]. For a gene and tumor type that have been profiled by proteogenomics, the observed correlation can be looked up in resources such as LinkedOmicsKB[@Liao2023] or OnCorr[@Nawaz2026]. For other genes, for tumor types without a proteome and for new cohorts, no lookup exists. Where one exists, it was estimated in a different cohort of about a hundred tumors, with its own sampling error.

Predicting protein abundance from transcriptomes has been studied extensively, including in the NCI-CPTAC DREAM challenge[@Yang2020] and in models that share information across genes[@LiH2019; @PCB2022]. Here we ask a different question: not what a protein's abundance is, but how far a gene's own mRNA can be trusted as a stand-in for its protein in a cohort where protein was not measured. We call a model that answers this question a concordance forecast. A useful forecast should rank genes by their concordance in a new cohort, should be at least as accurate as the lookups that users would otherwise consult, and should come with prediction intervals that hold on the platform at hand.

We trained concordance forecasts on correlations in 1,027 tumors of 10 cancer types from the Clinical Proteomic Tumor Analysis Consortium (CPTAC) with a leave-one-cancer-out design. We then preregistered tests on three external platforms that differ from CPTAC in technology and sample type: reverse-phase protein arrays (RPPA) of TCGA tumors[@LiJ2013], data-independent acquisition (DIA) mass spectrometry of cancer cell lines (ProCan)[@Goncalves2022] and tandem mass tag (TMT) proteomics of lung adenocarcinoma (APOLLO)[@Soltis2022]. The same preregistrations tested an existing type of protein predictor, an anchor-residual model. After the confirmatory tests, we compared the forecast with lookups of observed correlations, studied how to calibrate prediction intervals on a new platform with a few local measurements, and computed forecasts for ADC targets across 33 TCGA cancer types.

## Results

### Concordance and its forecast

We define the concordance of a gene in a context, such as a cancer type, a cell-line lineage or a cohort, as the Spearman correlation across samples between the gene's mRNA level and its protein level (ρ; at least 30 samples with both values). In CPTAC, concordance was computed in each of 10 cancer types for 7,432 genes, the P gene set (Table 1).
<!-- src: results/cptac_sample_counts.tsv; results/P_genes.txt; table1_data.tsv -->

The forecast, F1s, is a gradient-boosted tree model (LightGBM[@Ke2017]) that predicts ρ for a gene in a query context from nine rank- or percentile-based features. The features summarize the gene's concordance across CPTAC cancer types and describe the gene's RNA distribution in the query context (Methods). [CHECK: confirm this two-part description against the nine F1s features in config/analysis_plan_B3.yaml; if any feature is of another kind, rewrite these two sentences to match.] Because the query-side features use only RNA, the forecast can be computed for any cohort with RNA data.

We compared the forecast with two simpler predictors. The RNA-variance baseline ranks genes by the percentile of their RNA standard deviation in the query context, on the reasoning that a gene whose mRNA barely varies across samples cannot correlate strongly with its protein. The CPTAC mean lookup, F0, is the mean of the gene's observed CPTAC correlations, which is what a user consulting a proteogenomic resource would read off; it falls back to the global mean of 0.456 for genes without any CPTAC correlation.
<!-- src: results/stage7/f0_fill_counts.tsv (global mean 0.456 from 72,880 pairs) -->

In leave-one-cancer-out development, the Spearman correlation between the forecast and observed concordance in the held-out cancer type ranged from [FILL: min] to [FILL: max] across the 10 cancer types, against [FILL: min] to [FILL: max] for the RNA-variance baseline (Fig. 1a).
<!-- src: results/figures/source/fig1_loco_reliability.tsv -->

### Preregistered external tests

The external tests were preregistered in two stages, each before the protein data of the platforms concerned were used: RPPA and ProCan (hypotheses H1–H4) and APOLLO (P1, P2, Q1 and Q2). Predictions were frozen and committed before the protein data were unlocked, and each confirmation was run once (Methods; Supplementary Note 3). Table 2 lists all eight tests.

The forecast ranked observed concordance better than the RNA-variance baseline on all three platforms (Table 2, Fig. 2a). On TCGA RPPA (4,045 antibody–cancer pairs; 134 antibodies), the Spearman correlation between forecast and observed concordance was 0.432 (95% confidence interval 0.340–0.517), against 0.227 (0.119–0.329) for the baseline (H3; difference 0.205, 0.110–0.295; P = 0.0005). On ProCan cell lines (6,118 protein–lineage pairs; 2,746 proteins), it was 0.417 (0.390–0.444) against 0.141 (0.109–0.174) (H4; difference 0.276, 0.245–0.304; P = 0.0005). On APOLLO (7,196 genes), it was 0.747 (0.734–0.758) against 0.492 (0.473–0.510) (P2; difference 0.255, 0.238–0.273; P = 0.0005, the smallest value attainable with 2,000 resamples).
<!-- src: results/stage4/confirm/hypotheses.tsv (H3, H4); results/stage6/confirm/hypotheses.tsv (P2) -->

In a pre-specified descriptive analysis of APOLLO, the forecast identified genes with an observed concordance of at least 0.5 with an area under the receiver operating characteristic curve (AUROC) of 0.868, against 0.753 for the RNA-variance baseline and 0.857 for the CPTAC mean lookup (Supplementary Table S5).
<!-- src: results/stage6/confirm/descriptive_predictors.tsv -->

### Comparison with looking up observed correlations

The preregistered baseline was deliberately simple. A user with access to CPTAC would more likely look up the observed correlation, preferably from the same cancer type. We therefore compared the forecast with lookups in post hoc analyses whose definitions were committed before any lookup value was computed (Methods). The cascade lookup, L*, takes for each gene and context the first available of four values: the CPTAC correlation in the matched cancer type (L1), the mean of the gene's CPTAC correlations (L2), the mean ProCan correlation (L3; not used when ProCan itself was evaluated) and the global mean (L4).
<!-- src: config/analysis_plan_B8.yaml -->

The forecast was more accurate than the cascade lookup on all three platforms (Table 3, Fig. 2b). The difference in Spearman correlation with observed concordance was 0.069 (0.021–0.129) on RPPA, 0.029 (0.011–0.044) on ProCan and 0.021 (0.011–0.032) on APOLLO. Among pairs with a matched-cancer CPTAC correlation (L1), the forecast was more accurate on ProCan (0.026, 0.002–0.050) and APOLLO (0.018, 0.008–0.028), and the interval included zero on RPPA (0.026, −0.025 to 0.074). Using the matched correlation wherever it existed and the forecast elsewhere, a hybrid rule, did not improve on the forecast alone: the differences were −0.008 (−0.024 to 0.009) on RPPA, −0.009 (−0.022 to 0.004) on ProCan and −0.016 (−0.026 to −0.007) on APOLLO. Within RPPA cancer types, the forecast ranked antibodies better than the cascade lookup in 28 of 31 cancer types (median difference 0.080, 0.057–0.098; two-sided sign test P = 4.6 × 10^−6). AUROC differences pointed in the same direction (Supplementary Table S3).
<!-- src: results/stage8/spearman_diff.tsv; results/stage8/within_context_sign.tsv; results/stage8/auroc_diff.tsv -->

[FILL: one or two sentences from the B10 §2 analysis (pooled lookup versus matched lookup on L1 pairs), written by the rule in B 지시서 10 §2.4.]

For genes outside the P set, for which CPTAC information was sparse or absent, the forecast gained most over the CPTAC mean lookup: 0.139 (0.070–0.209) on ProCan and 0.276 (0.220–0.328) on APOLLO. Inside the P set, the gain on APOLLO was close to zero (0.002, −0.003 to 0.007) (Supplementary Table S2).
<!-- src: results/stage7/f0_diff.tsv -->

### Calibrating intervals on a new platform

The forecast comes with 80% and 90% intervals whose half-widths were set from CPTAC holdout errors (Methods). [CHECK: how the half-widths 0.202 and 0.259 were derived.] Forecast errors shifted between platforms: the mean of observed minus predicted concordance was −0.184 on RPPA lung adenocarcinoma, −0.026 on the ProCan lung lineage and −0.023 on APOLLO, against 0.047 on the CPTAC lung adenocarcinoma holdout. Without adjustment, the 80% intervals covered 0.511 of RPPA pairs, 0.628 of ProCan pairs and 0.833 of APOLLO genes (Fig. 3).
<!-- src: results/stage6/posthoc/posthoc_02_bias.tsv; results/stage7/calibration_table.tsv -->

We therefore evaluated two adjustments that use K genes measured on the new platform. The first shifts every forecast by δ, the mean of observed minus predicted concordance in the K genes. The second applies the same shift and replaces the interval half-width with the 80% quantile of absolute errors in the K genes. The shift alone did not bring RPPA or ProCan coverage into the range 0.75–0.85. With the shift and quantile, mean coverage was inside this range on all three platforms at K = 50 (RPPA 0.798, ProCan 0.789, APOLLO 0.779) and K = 100 (0.801, 0.792, 0.791), but not at K = 20 on APOLLO (0.742) (Fig. 3, Supplementary Table S6). The RPPA and ProCan draws were post hoc and the APOLLO draws pre-specified descriptive, so the recommendation of at least 50 calibration genes is post hoc.
<!-- src: results/stage7/calibration_table.tsv; results/stage7/calibration_recommend.tsv -->

### An anchor-residual protein predictor

The preregistrations also tested an anchor-residual (AR) predictor of protein abundance, which combines a gene's own mRNA, the anchor, with a ridge model of the remaining protein variation built from the mRNA of 50 other genes (Methods). [CHECK: this one-sentence description against the stage-5 AR code and config; correct it if it differs.] Predictors of this kind are not new[@Yang2020; @LiH2019]; we tested this one externally because its accuracy determines whether a protein estimate can improve on the gene's own mRNA in cohorts without proteomes. The comparison was the correlation of the AR prediction with protein (ρ_AR) against the correlation of the gene's own mRNA with protein (ρ_S0). In CPTAC development, the mean difference ρ_AR − ρ_S0 ranged from 0.085 to 0.121 across held-out cancer types (Fig. 1b).
<!-- src: results/figures/source/fig1_loco_delta.tsv (CCRCC 0.0854 to LSCC 0.121) -->

Externally, AR improved on the gene's own mRNA on RPPA (H1; mean difference 0.034, 0.024–0.043; 123 antibodies; P = 1.5 × 10^−13) and on APOLLO (P1; 0.085, 0.082–0.088; 7,029 genes), and it outperformed two anchor-free ridge models on APOLLO (Q1, 0.023; Q2, 0.063; 6,322 genes) (Table 2, Fig. 4a). It did not improve on the gene's own mRNA in ProCan cell lines (H2; 0.001, −0.002 to 0.004; 2,623 proteins; P = 0.10). In a post hoc description across 34 contexts (31 RPPA cancer types and 3 ProCan lineages), the median AR gain in a context increased with the similarity of its RNA profile to CPTAC (Spearman correlation 0.521; Fig. 4b). APOLLO's similarity, 0.925, lay inside the range of the CPTAC holdout cancer types (0.883–0.979). [CHECK: which contexts have the lowest similarity in results/figures/source/fig4_scope.tsv; if the lowest is the ProCan haematopoietic and lymphoid lineage, add "The ProCan haematopoietic and lymphoid lineage had the lowest similarity ([FILL: s]) and a median gain of [FILL: value]."]
<!-- src: results/stage4/confirm/hypotheses.tsv (H1, H2); results/stage6/confirm/hypotheses.tsv (P1, Q1, Q2); results/stage5/posthoc/posthoc_22_similarity.tsv; results/stage6/confirm/s.tsv -->

### Concordance forecasts for ADC targets

We computed forecasts for 29 ADC targets in each of 33 TCGA cancer types from TCGA RNA (957 forecasts; post hoc; Fig. 5a). [OPTION A: The same procedure gave forecasts for [FILL: number of genes] genes, with the RNA-profile similarity of each cancer type to CPTAC as a scope indicator (Data availability).] Where the target protein had been measured externally, forecasts agreed with observed concordance with Spearman correlations of 0.308 on RPPA (124 target–cancer pairs), 0.584 on APOLLO (21 targets) and 0.793 on ProCan (12 target–lineage pairs) [CHECK: which ProCan series in results/stage7/adc_concordance.tsv is the external comparison, and what the "ProCan frozen forecast" series is; keep one series in the text and explain the other in the Fig. 5 legend]. Twenty of the 21 APOLLO targets fell inside their 80% interval (Fig. 5b). The exception was CLDN18 (observed concordance 0.388, forecast 0.803). The clinically targeted isoform, CLDN18.2, cannot be separated from other CLDN18 isoforms by gene-level RNA or protein measurements. Eight approved ADCs act on six of the 29 targets (Supplementary Table S7). [CHECK: update both counts from B10 §4.]
<!-- src: results/stage7/adc_atlas_final.tsv; results/stage7/adc_concordance.tsv; results/stage6/confirm/descriptive_adc.tsv; results/stage9/adc_approved_full.tsv -->

## Discussion

The central result is that the concordance between a gene's mRNA and its protein can be forecast in cohorts where protein was not measured, and that the forecast transferred to three platforms that differ from the training data in technology and sample type. The preregistered tests compared the forecast with an RNA-variance baseline. The post hoc comparisons added a more demanding benchmark, looking up observed correlations, and the forecast was at least as accurate as the lookup even when the lookup came from the matched cancer type.

Why can a forecast compete with a measurement? A matched-cancer correlation is a measurement in a different cohort. CPTAC cohorts have 82 to 121 tumors per cancer type, so each correlation carries sampling error, and cohorts differ in sample handling and assays. The forecast instead combines what is known about the gene across all CPTAC cancer types with the gene's RNA distribution in the query cohort itself. The second part matters because a correlation is limited by how much its variables vary in the sample: a gene whose RNA varies little across the tumors of a cohort cannot have high concordance there, whatever its concordance elsewhere. The RNA-variance baseline alone carried this information (Fig. 2a), and the forecast combines it with gene-level information. [FILL: one sentence from B10 §2 by the rule in B 지시서 10 §2.4.]
<!-- src: results/cptac_sample_counts.tsv (OV 82, BRCA 121) -->

These results suggest rules of practice. They rest on post hoc analyses and should be read with that status. First, for a gene and context without a proteome, the forecast with its interval is a better guide than either a cross-cancer average or the RNA variance alone. Second, a matched-cancer correlation from another cohort need not be preferred over the forecast. Third, on a new platform, intervals should be recalibrated with at least 50 genes measured on that platform; uncalibrated intervals undercovered on RPPA and ProCan. Antibody-based assays deserve particular caution: on RPPA lung adenocarcinoma, observed concordance was on average 0.184 lower than forecast, in line with the known variability of antibody quality[@Upadhya2023].
<!-- src: results/stage6/posthoc/posthoc_02_bias.tsv -->

Concordance describes ranking across tumors. High concordance means that RNA orders tumors much as protein would, which is what matters when RNA is used to select tumors or patients with high target expression. It does not say whether a protein is abundant in absolute terms or present at the cell surface, both of which also matter for ADC targets. For CLDN18, the gene-level forecast was high, but the observed concordance in APOLLO was lower, and the clinically relevant isoform cannot be separated by gene-level measurements. Isoform-specific targets need isoform-level data.

The AR predictor was a secondary component. It improved on the gene's own mRNA in the two tumor cohorts but not in ProCan cell lines, and a post hoc description associated smaller gains with lower similarity to CPTAC. The failed test H2 is a preregistered result; the similarity description is post hoc and does not retest it. Protein predictors trained on tumors should not be assumed to transfer to cell lines.

The study has limitations. The comparisons with lookups and the calibration recommendation are post hoc. The preregistration was internal: plans and frozen predictions were committed to version control before unlocking, but the commit times were recorded by the analysis workstation and are not certified by a third party (Supplementary Note 3). The two APOLLO families were controlled separately, so the familywise error rate across them could reach about 0.10. Within an antibody, ranking cancer types by the forecast was not clearly better than ranking them by RNA variance (median difference 0.024, −0.016 to 0.096); the forecast is better suited to comparing genes than to comparing cancer types for one gene. Training used 10 CPTAC cancer types, and contexts whose RNA profiles differ from CPTAC may lie outside the forecast's scope. APOLLO protein values were imputed for proteins quantified in at least half of the samples; a post hoc analysis of the non-imputed matrix gave similar results (Supplementary Note 5). Patient overlap between APOLLO and CPTAC lung adenocarcinoma was excluded only by public identifiers.
<!-- src: results/stage7/f0_within_antibody_summary.tsv; results/stage6/posthoc/posthoc_05_raw.tsv -->

In summary, a concordance forecast trained on CPTAC tumors ranked mRNA–protein concordance on three external platforms better than an RNA-variance baseline in preregistered tests, and in post hoc comparisons it was at least as accurate as looking up observed correlations. With a short local calibration, it gives a practical estimate of how far mRNA can stand in for protein in cohorts without proteomes.

## Methods

### Data

**CPTAC.** Paired tumor RNA and protein for 10 cancer types were obtained from the LinkedOmicsKB data freeze [FILL: freeze name or version from logs/download_log.tsv][@Liao2023], which is built on the harmonized CPTAC pan-cancer data[@LiY2023] [CHECK: this relation from the LinkedOmicsKB page or Liao2023]: breast cancer (BRCA; 121 tumors)[@Krug2020], clear cell renal cell carcinoma (CCRCC; 103)[@Clark2019], colon adenocarcinoma (COAD; 96)[@Vasaikar2019], glioblastoma (GBM; 99)[@WangLB2021], head and neck squamous cell carcinoma (HNSCC; 108)[@Huang2021], lung squamous cell carcinoma (LSCC; 108)[@Satpathy2021], lung adenocarcinoma (LUAD; 110)[@Gillette2020], ovarian cancer (OV; 82)[@McDermott2020], pancreatic ductal adenocarcinoma (PDAC; 105)[@Cao2021] and uterine corpus endometrial carcinoma (UCEC; 95)[@Dou2020], 1,027 tumors in total. [CHECK: RNA and protein units in the freeze, and the definition of the P gene set (7,432 genes).]
<!-- src: results/cptac_sample_counts.tsv; results/P_genes.txt -->

**TCGA RPPA.** RPPA data for TCGA tumors were obtained from The Cancer Proteome Atlas (TCPA)[@LiJ2013] [FILL: data level and release], with 7,754 samples and 258 antibodies. [CHECK: the RNA matrix paired with RPPA samples, the antibody-to-gene mapping rule, and the minimum sample count per cancer type.] Concordance was evaluated for 4,045 antibody–cancer pairs (134 antibodies, 31 cancer types) with defined correlations.
<!-- src: table1_data.tsv; results/stage4/confirm/hypotheses.tsv -->

**ProCan.** DIA proteomes of cancer cell lines were those of the ProCan-DepMapSanger resource[@Goncalves2022] [FILL: file and download page], with 8,498 proteins. RNA was taken from DepMap 24Q4[@DepMap24Q4; @Ghandi2019]. [CHECK: that DepMap 24Q4 RNA is the RNA used for ProCan contexts.] Three lineages were evaluated: large intestine, lung, and haematopoietic and lymphoid (6,118 protein–lineage pairs; 2,746 proteins).

**APOLLO.** Protein (TMT log2 ratios to a common reference, proteins quantified in at least 50% of samples with k-nearest-neighbor imputation) and RNA (TPM) for 87 resected lung adenocarcinomas[@Soltis2022] were obtained from the open files of the Genomic Data Commons (ProteomeXchange PXD036025). File names, processing rules and hashes are given in Supplementary Note 1.

**Other data.** Cancer Cell Line Encyclopedia (CCLE) RNA[@Ghandi2019; @DepMap24Q4] and proteomes[@Nusinow2020] of lung (74 cell lines) and skin (33 cell lines) lineages were used during development only [CHECK: for what; write the purpose in one clause]. TCGA PanCanAtlas RNA[@Hoadley2018] was used to compute forecasts in 33 TCGA cancer types. Human Protein Atlas staining summaries[@Uhlen2017] were used for a descriptive comparison (Supplementary Note 7). Sources, file names and SHA-256 hashes of all downloaded files are listed in Supplementary Table S9.
<!-- src: table1_data.tsv -->

### Concordance

For a gene in a context, concordance ρ was the Spearman correlation between log2(TPM + 1) RNA and protein across samples in which both were finite, with average ranks for ties. It was computed when at least 30 such samples were available and was undefined when either vector was constant. [CHECK: RPPA and ProCan RNA units, if not TPM.]

### Forecast models and baselines

The training data were CPTAC leave-one-cancer-out rows: for each held-out cancer type, the target was ρ in that cancer type, gene-level features were computed from the other nine cancer types and query features from the RNA of the held-out cancer type. [CHECK against the stage-3 code.] The nine features of F1s were [FILL: list each feature with a one-line definition from config/analysis_plan_B3.yaml]. F1 used the gene-level features only. Each model was fitted as a ridge regression and as LightGBM[@Ke2017] (300 trees, learning rate 0.05, 31 leaves, minimum 20 samples per leaf, row subsample 0.8 at every iteration, column subsample 0.8, random state 20261001); F1s-LightGBM was the preregistered forecast. Missing feature values were replaced by the column median of the training matrix. [CHECK: how the final model applied to external data was fitted (all leave-one-cancer-out rows or all CPTAC rows).] The 80% and 90% interval half-widths were 0.202 and 0.259 [CHECK: derivation], and intervals were not clipped to [−1, 1].
<!-- src: reports/B_methods_facts.md; results/stage9/refit_compare.tsv -->

The RNA-variance baseline was the percentile of the gene's RNA standard deviation among genes in the query context. The CPTAC mean lookup F0 was the mean of the gene's defined CPTAC correlations stored at development; genes without any stored CPTAC correlation received the global mean of P-gene correlations, 0.456 (72,880 gene–cancer pairs).
<!-- src: results/stage7/f0_fill_counts.tsv -->

### Anchor-residual predictor

[FILL: three to five sentences from the stage-5 code and config: the anchor; how the 50 residual genes are chosen; the residual ridge model and its penalty; the λ grid and how λ is selected; the anchor-free baselines R1 and R2 (R2 ridge α = 1000).] The own-mRNA baseline S0 is the gene's own RNA, so ρ_S0 equals the observed concordance [CHECK]. For APOLLO, predictions used the stored final model without refitting.

### Preregistration and confirmatory tests

The tests on RPPA and ProCan (H1–H4) and on APOLLO (P1, P2, Q1, Q2) were preregistered in two analysis plans committed to version control on 2 October 2026, each before the protein data of the platforms concerned were used. Predictions were frozen and committed before an unlock commit, and the confirmation scripts were run once. Commit hashes and times, power analyses, tests on simulated planted and null data, and the departures accepted before unlocking are reported in Supplementary Notes 3 and 4.

H1, H2, P1, Q1 and Q2 compared unit-level differences in correlation with protein (ρ_AR − ρ_S0, ρ_AR − ρ_R1 or ρ_AR − ρ_R2) with a one-sided Wilcoxon signed-rank test[@Wilcoxon1945]; units were antibodies for RPPA, proteins for ProCan and genes for APOLLO. [CHECK: for H1 and H2, how differences were summarized across contexts within a unit.] H3, H4 and P2 compared the Spearman correlation of the forecast with observed concordance against that of the RNA-variance baseline. The P value was one plus the number of 2,000 unit resamples in which the difference was at or below zero, divided by 2,001. Confidence intervals were the 2.5th and 97.5th percentiles of 2,000 unit resamples. Holm's procedure[@Holm1979] controlled the familywise error rate at 0.05 within each of three families: H1–H4; P1 and P2 (primary); and Q1 and Q2 (secondary). [CHECK: thresholds per family against config/prereg_B4 and prereg_B6.] Each test restarted the generator with seed 20261001.

### Post hoc analyses

All analyses in this section were specified after the confirmations. The lookup comparisons (cascade lookup, hybrid rule and pooled lookup) were defined in plans committed before any lookup value was computed (Supplementary Note 3). [CHECK: whether the stage-5 and stage-7 post hoc analyses also had plan commits before calculation; if so, say so here and list them in Supplementary Note 3; if not, leave this sentence as it is.]

The cascade lookup L* and the hybrid rule were defined as in the Results, with the context map in Supplementary Table S10. L1 and L2 were computed from CPTAC correlations for all genes quantified in each cohort, which added 9,318 gene–cohort correlations to those stored at development [CHECK: wording against results/stage8/lookup_build.tsv]; F0 used only the stored correlations. A difference involving a constant predictor was undefined. [FILL: definition of the pooled lookup used in B10 §2, one sentence.]

Differences in Spearman correlation and AUROC (positive class: observed ρ of at least 0.5) were estimated with 2,000 resamples of units (antibodies, proteins or genes) with the generator restarted at seed 20261001 for every contrast; a resample kept all pairs of each sampled unit, and intervals were the 2.5th and 97.5th percentiles of resamples with a defined difference. Within-context comparisons used a two-sided sign test across RPPA cancer types.

### Calibration

For each platform, units were split into five folds, and calibration was repeated 20 times with seed 20261014, giving 100 draws for each K of 10, 20, 50 and 100. [CHECK: how the K calibration units and the evaluation units were drawn from the folds.] The shift δ was the mean of observed minus predicted concordance in the K units. The shift-and-quantile method also replaced the half-width by the 0.80 (or 0.90) quantile of absolute errors in the K units, an approach related to split-conformal calibration[@Lei2018]. Coverage was the fraction of evaluation units whose observed concordance lay inside the interval, averaged over draws. A method was recommended only if its mean 80% coverage lay within 0.75–0.85, inclusive, on all three platforms.

### RNA-profile similarity

[FILL: definition of the similarity s between a context and CPTAC, from the stage-5 code.] The reference range was the range of s over the CPTAC holdout cancer types (0.883–0.979).
<!-- src: results/stage6/confirm/s.tsv -->

### ADC target atlas

The 29 targets were [CHECK: the selection rule and source behind config/adc_targets.tsv]. For each target and TCGA cancer type, the forecast and its 80% interval were computed from TCGA PanCanAtlas RNA. Observed concordance was compared where the target protein had been measured. Approved ADCs were identified from US Food and Drug Administration prescribing information as of 2 October 2026 (Supplementary Table S7).

### Software

Analyses used Python 3.12.4 with NumPy 1.26.4[@Harris2020], pandas 2.2.2[@McKinney2010], SciPy 1.13.1[@Virtanen2020], scikit-learn 1.4.2[@Pedregosa2011], LightGBM 4.7.0[@Ke2017] and Matplotlib 3.8.4[@Hunter2007] on a CPU server without a GPU. [CHECK: that these versions were the versions at the time of the frozen predictions.]

### Use of large language models

Large language models were used throughout this study under the author's direction. Claude (Anthropic; [AUTHOR TO FILL: model and version]) drafted the analysis plans, the stage-by-stage task instructions and the text of this manuscript. A coding agent in the Cursor editor ([AUTHOR TO FILL: model and version]) wrote and ran the analysis code, produced the figures and earlier drafts of the text, and ran automated checks of the references against Crossref and of every number in the manuscript against the result files. The author reviewed each analysis plan before it was committed and executed; preregistered predictions were frozen before the external protein data were unlocked. The author reviewed all plans, code outputs and text and takes full responsibility for the content. Language models are not authors of this work. [AUTHOR TO CONFIRM: this paragraph]

## Data availability

All data analysed in this study are public. CPTAC data were obtained from LinkedOmicsKB (https://kb.linkedomics.org/), TCGA RPPA data from TCPA ([FILL: URL]), ProCan proteomes from [FILL: page] (PRIDE PXD030304), APOLLO data from the Genomic Data Commons ([FILL: publication page URL]; ProteomeXchange PXD036025), CCLE proteomes from [FILL: URL], DepMap 24Q4 data from https://doi.org/10.25452/figshare.plus.27993248.v1 and TCGA PanCanAtlas RNA from https://gdc.cancer.gov/about-data/publications/pancanatlas. File names, sources and SHA-256 hashes are listed in Supplementary Table S9. Code, frozen predictions, the trained forecast model, forecast tables [OPTION A: for all genes and] for the ADC targets, and figure source data are available at [FILL: GitHub URL] and archived at Zenodo ([FILL: DOI]).

## References

[Converted from citation keys by the agent; Scientific Reports style.]

## Acknowledgements

The results shown here are in whole or part based upon data generated by the TCGA Research Network: https://www.cancer.gov/tcga. Data used in this publication were generated by the Clinical Proteomic Tumor Analysis Consortium (NCI/NIH). [AUTHOR TO FILL: funding, with grant numbers.]

## Author contributions

[AUTHOR TO FILL]

## Competing interests

The author(s) declare no competing interests. [AUTHOR TO CONFIRM]

## Figure legends

**Figure 1. Development of the concordance forecast and the anchor-residual predictor on CPTAC tumors.** (a) Spearman correlation between predicted and observed concordance in each held-out cancer type in leave-one-cancer-out analysis, for the forecast (F1s-LightGBM) and the RNA-variance baseline. (b) Mean difference between the correlation of the anchor-residual prediction with protein and that of the gene's own mRNA (ρ_AR − ρ_S0) in each held-out cancer type; genes per cancer type: 6,778–7,426. Development analyses, not preregistered tests.
<!-- src: results/figures/source/fig1_loco_reliability.tsv; results/figures/source/fig1_loco_delta.tsv -->

**Figure 2. External evaluation of the concordance forecast.** (a) Spearman correlation between predicted and observed concordance on TCGA RPPA (4,045 antibody–cancer pairs; 134 antibodies), ProCan cell lines (6,118 protein–lineage pairs; 2,746 proteins) and APOLLO lung adenocarcinoma (7,196 genes). Circles: forecast and RNA-variance baseline with preregistered 95% intervals. Squares: CPTAC mean lookup and cascade lookup (post hoc). (b) Differences in Spearman correlation with 95% intervals: forecast minus RNA-variance baseline (preregistered tests H3, H4 and P2; circles), and, post hoc (diamonds), forecast minus cascade lookup on all pairs and on pairs with a matched-cancer CPTAC correlation, and hybrid rule minus forecast. Intervals are 2.5th–97.5th percentiles of 2,000 resamples of antibodies, proteins or genes.

**Figure 3. Calibration of 80% forecast intervals with locally measured genes.** Mean coverage of nominal 80% intervals over 100 draws (five folds by 20 repeats) with K = 10, 20, 50 or 100 calibration genes, for no adjustment, a shift by the mean error δ, and the shift plus an empirical quantile, on RPPA and ProCan (post hoc) and APOLLO (pre-specified descriptive). Error bars, 2.5th–97.5th percentiles of draws; grey band, coverage 0.75–0.85; dashed line, 0.80.

**Figure 4. External tests of the anchor-residual predictor.** (a) Mean unit-level differences with 95% intervals for the preregistered tests H1 (RPPA, 123 antibodies), H2 (ProCan, 2,623 proteins), P1 (APOLLO, 7,029 genes), and Q1 and Q2 (APOLLO, 6,322 genes). Grey squares: post hoc comparisons with the anchor-free baselines on RPPA (100 antibodies) and ProCan (2,399 proteins), without intervals. (b) Median difference ρ_AR − ρ_S0 per context against RNA-profile similarity to CPTAC, for 31 RPPA cancer types and 3 ProCan lineages (post hoc; circles) and APOLLO (preregistered; diamond). Dashed lines mark the range of similarity across CPTAC holdout cancer types (0.883–0.979).
<!-- src: results/figures/source/fig2_effects.tsv; results/figures/source/fig4_scope.tsv -->

**Figure 5. Concordance forecasts for antibody–drug conjugate targets.** (a) Forecast concordance for 29 ADC targets in 33 TCGA cancer types (post hoc). CLDN18 is marked because gene-level assays cannot separate the targeted isoform CLDN18.2. (b) Forecast against observed concordance for targets whose protein was measured: RPPA (124 target–cancer pairs; Spearman correlation 0.308), ProCan ([CHECK: series], 12 pairs; 0.793) and APOLLO (21 targets; 0.584), with CPTAC development data (228 pairs; 0.554) for reference [CHECK: whether the CPTAC points use leave-one-cancer-out predictions; say so, or say that they are in-sample]. Open circle, CLDN18 in APOLLO.
<!-- src: results/stage7/adc_concordance.tsv -->

## Tables

**Table 1. Data sets.**

| Data set | Platform | Samples | Units analysed | Role |
|---|---|---|---|---|
| CPTAC, 10 cancer types | RNA-seq; mass spectrometry | 1,027 tumors | 7,432 genes (P set) | Training and development |
| CCLE lung and skin lineages | RNA-seq; TMT mass spectrometry | 74 and 33 cell lines | [FILL] genes | Development |
| TCGA RPPA, 31 cancer types | RNA-seq; reverse-phase protein array | 7,754 samples in the matrix | 134 antibodies (4,045 antibody–cancer pairs) | Preregistered tests H1, H3 |
| ProCan, 3 lineages | RNA-seq; DIA mass spectrometry | [FILL] cell lines in the 3 lineages | 2,746 proteins (6,118 protein–lineage pairs) | Preregistered tests H2, H4 |
| APOLLO lung adenocarcinoma | RNA-seq; TMT mass spectrometry | 87 tumors | 7,196 genes | Preregistered tests P1, P2, Q1, Q2 |
| TCGA PanCanAtlas, 33 cancer types | RNA-seq | [FILL] tumors | 29 ADC targets [OPTION A: and [FILL] genes] | ADC atlas (post hoc) |

<!-- src: results/figures/source/table1_data.tsv; results/cptac_sample_counts.tsv. Units for H1 and H2 differ from those for H3 and H4 (123 antibodies and 2,623 proteins); give both in a table note. -->

Notes: Units analysed are those of the concordance tests (H3, H4, P2). The anchor-residual tests used 123 RPPA antibodies (H1), 2,623 ProCan proteins (H2), 7,029 APOLLO genes (P1) and 6,322 APOLLO genes (Q1, Q2).

**Table 2. Preregistered tests.**

| Test | Data set | Comparison | Units | Estimate (95% CI) | P | Holm decision |
|---|---|---|---|---|---|---|
| H1 | TCGA RPPA | ρ_AR − ρ_S0, mean | 123 antibodies | 0.0339 (0.0242–0.0432) | 1.5 × 10^−13 | Rejected |
| H2 | ProCan | ρ_AR − ρ_S0, mean | 2,623 proteins | 0.000932 (−0.00236 to 0.00423) | 0.10 | Not rejected |
| H3 | TCGA RPPA | Forecast − RNA-variance baseline, Spearman | 134 antibodies (4,045 pairs) | 0.205 (0.110–0.295) | 0.0005 | Rejected |
| H4 | ProCan | Forecast − RNA-variance baseline, Spearman | 2,746 proteins (6,118 pairs) | 0.276 (0.245–0.304) | 0.0005 | Rejected |
| P1 | APOLLO | ρ_AR − ρ_S0, mean | 7,029 genes | 0.0850 (0.0822–0.0878) | [FILL] | Rejected |
| P2 | APOLLO | Forecast − RNA-variance baseline, Spearman | 7,196 genes | 0.255 (0.238–0.273) | 0.0005 | Rejected |
| Q1 | APOLLO | ρ_AR − ρ_R1, mean | 6,322 genes | 0.0228 ([FILL]) | 4.5 × 10^−203 | Rejected |
| Q2 | APOLLO | ρ_AR − ρ_R2, mean | 6,322 genes | 0.0627 ([FILL]) | [FILL] | Rejected |

<!-- src: results/stage4/confirm/hypotheses.tsv; results/stage4/confirm/holm.tsv; results/stage6/confirm/hypotheses.tsv; results/stage6/confirm/holm.tsv -->

Notes: Wilcoxon signed-rank tests (one-sided) for H1, H2, P1, Q1 and Q2; resampling tests of the difference in Spearman correlation for H3, H4 and P2 (2,000 resamples; smallest attainable P = 0.0005). Intervals are 2.5th–97.5th percentiles of 2,000 unit resamples. Holm families: H1–H4; P1 and P2; Q1 and Q2. [FILL: how P values that underflow to 0 are reported, e.g. "P < 10^−300", after checking the test statistic.]

**Table 3. Post hoc comparisons with lookups of observed CPTAC correlations.**

| Contrast | Pairs | RPPA | ProCan | APOLLO |
|---|---|---|---|---|
| Forecast − cascade lookup | All | 0.069 (0.021 to 0.129) | 0.029 (0.011 to 0.044) | 0.021 (0.011 to 0.032) |
| Forecast − cascade lookup | Matched cancer (L1) | 0.026 (−0.025 to 0.074) | 0.026 (0.002 to 0.050) | 0.018 (0.008 to 0.028) |
| Forecast − CPTAC mean lookup | All | 0.064 (0.021 to 0.112) | 0.032 (0.019 to 0.044) | 0.028 (0.021 to 0.035) |
| Hybrid rule − forecast | All | −0.008 (−0.024 to 0.009) | −0.009 (−0.022 to 0.004) | −0.016 (−0.026 to −0.007) |
| [FILL: B10 §2 rows] | Matched cancer (L1) | [FILL] | [FILL] | [FILL] |

<!-- src: results/stage8/spearman_diff.tsv; results/stage9/f0_lstar_diff.tsv; results/stage10/ (B10 §2) -->

Notes: Differences in Spearman correlation with observed concordance; 95% intervals from 2,000 resamples of antibodies (RPPA), proteins (ProCan) or genes (APOLLO). All pairs: 4,045 (134 antibodies), 6,118 (2,746 proteins) and 7,196 genes. Matched cancer (L1): 1,206 (127 antibodies), 3,356 (2,091 proteins) and 6,982 genes. The CPTAC mean lookup used the CPTAC correlations stored at development, and the cascade lookup an expanded set (Methods).
