# A preregistered model forecasts mRNA–protein concordance in unmeasured cancer contexts across three proteomic platforms



Juwon Kang^1 (ORCID 0009-0009-2186-0038) and Junjeong Choi^1* (ORCID 0000-0003-1339-593X)<br><br>^1 College of Pharmacy, Yonsei Institute of Pharmaceutical Sciences, Yonsei University, Incheon 21983, Republic of Korea<br><br>* Corresponding author: Junjeong Choi, junjeong@yonsei.ac.kr<br><br>E-mail addresses: Juwon Kang, juwonk@yonsei.ac.kr; Junjeong Choi, junjeong@yonsei.ac.kr

## Abstract

Many tumor cohorts have transcriptomes but no matched proteomes, so mRNA is used as a stand-in for protein. Whether this is reasonable depends on each gene's mRNA–protein concordance in each context, known only where both were measured. We trained a model that forecasts this concordance from gene-level information and the gene's RNA distribution in the query cohort, using 1,027 CPTAC tumors of 10 cancer types, and preregistered tests on three external platforms. The forecast ranked observed concordance better than an RNA-variance baseline on reverse-phase protein arrays of TCGA tumors (Spearman correlation 0.43 versus 0.23), on DIA proteomes of cancer cell lines (0.42 versus 0.14) and on a TMT lung adenocarcinoma cohort (0.75 versus 0.49). Post hoc, it was at least as accurate as looking up observed CPTAC correlations, pooled or cancer-matched, for genes whose CPTAC data it had been given. Recalibrating intervals with 50 locally measured genes brought 80% coverage to 0.78–0.80 on all three platforms. An anchor-residual protein predictor tested alongside improved on the gene's own mRNA in both tumor cohorts but not in cell lines. We provide forecasts for 19,211 genes across 33 cancer types, including 29 antibody–drug conjugate targets.


## Introduction

Large tumor cohorts such as The Cancer Genome Atlas (TCGA) profiled RNA in thousands of tumors<sup>1</sup>, whereas matched proteomes exist for far fewer<sup>2</sup>. As a result, mRNA is routinely used in place of protein when a question is really about protein: whether a candidate drug target is expressed, which tumors carry more of it, or which cancer types to prioritize for an antibody–drug conjugate (ADC)<sup>3</sup><sup>4</sup><sup>5</sup>.

This substitution is only as good as the concordance between a gene's mRNA and its protein across tumors. Concordance varies widely between genes<sup>6</sup><sup>7</sup>, is bounded by the reproducibility of both measurements<sup>8</sup> and is lowered further when protein is measured with antibodies of uneven quality<sup>9</sup>. For a gene and tumor type that have been profiled by proteogenomics, the observed correlation can be looked up in resources such as LinkedOmicsKB<sup>10</sup> or OnCorr<sup>11</sup>. For other genes, for tumor types without a proteome and for new cohorts, no lookup exists. Where one exists, it was estimated in a different cohort of about a hundred tumors, with its own sampling error.

Predicting protein abundance from transcriptomes has been studied extensively, including in the NCI-CPTAC DREAM challenge<sup>12</sup> and in models that share information across genes<sup>13</sup><sup>14</sup>. Here we ask a different question: not what a protein's abundance is, but how far a gene's own mRNA can be trusted as a stand-in for its protein in a cohort where protein was not measured. We call a model that answers this question a concordance forecast. A useful forecast should rank genes by their concordance in a new cohort, should be at least as accurate as the lookups that users would otherwise consult, and should come with prediction intervals that hold on the platform at hand.

We trained concordance forecasts on correlations in 1,027 tumors of 10 cancer types from the Clinical Proteomic Tumor Analysis Consortium (CPTAC) with a leave-one-cancer-out design. We then preregistered tests on three external platforms that differ from CPTAC in technology and sample type: reverse-phase protein arrays (RPPA) of TCGA tumors<sup>15</sup>, data-independent acquisition (DIA) mass spectrometry of cancer cell lines (ProCan)<sup>16</sup> and tandem mass tag (TMT) proteomics of lung adenocarcinoma (APOLLO)<sup>17</sup>. The same preregistrations tested an existing type of protein predictor, an anchor-residual model. After the confirmatory tests, we compared the forecast with lookups of observed correlations, studied how to calibrate prediction intervals on a new platform with a few local measurements, and computed forecasts for ADC targets across 33 TCGA cancer types.

## Results

### Concordance and its forecast

We define the concordance of a gene in a context, such as a cancer type, a cell-line lineage or a cohort, as the Spearman correlation across samples between the gene's mRNA level and its protein level (ρ; at least 30 samples with both values). In CPTAC, concordance was computed in each of 10 cancer types for 7,432 genes, the P gene set (Table 1).


The forecast, F1s, is a gradient-boosted tree model (LightGBM<sup>18</sup>) that predicts ρ for a gene in a query context from nine rank- or percentile-based features. Five features summarize the gene's concordance across CPTAC cancer types. Two describe the gene's RNA distribution in the query context, one is the mean percentile of its RNA mean across the training cancer types, and one is the mean fraction of CPTAC tumors in which its protein was measured (Methods). Because the query-side features use only RNA, the forecast can be computed for any cohort with RNA data.

We compared the forecast with two simpler predictors. The RNA-variance baseline ranks genes by the percentile of their RNA standard deviation in the query context, on the reasoning that a gene whose mRNA barely varies across samples cannot correlate strongly with its protein. The CPTAC mean lookup, F0, is the mean of the gene's observed CPTAC correlations, which is what a user consulting a proteogenomic resource would read off; it falls back to the global mean of 0.456 for genes without any CPTAC correlation.


In leave-one-cancer-out development, the Spearman correlation between the forecast and observed concordance in the held-out cancer type ranged from 0.696 to 0.868 across the 10 cancer types, against 0.180 to 0.530 for the RNA-variance baseline (Fig. 1a).


### Preregistered external tests

The external tests were preregistered in two stages, each before the protein data of the platforms concerned were used: RPPA and ProCan (hypotheses H1–H4) and APOLLO (P1, P2, Q1 and Q2). Predictions were frozen and committed before the protein data were unlocked, and each confirmation was run once (Methods; Supplementary Note 3). Table 2 lists all eight tests.

The forecast ranked observed concordance better than the RNA-variance baseline on all three platforms (Table 2, Fig. 2a). On TCGA RPPA (4,045 antibody–cancer pairs; 134 antibodies), the Spearman correlation between forecast and observed concordance was 0.432 (95% confidence interval 0.340–0.517), against 0.227 (0.119–0.329) for the baseline (H3; difference 0.205, 0.110–0.295; P = 0.0005). On ProCan cell lines (6,118 protein–lineage pairs; 2,746 proteins), it was 0.417 (0.390–0.444) against 0.141 (0.109–0.174) (H4; difference 0.276, 0.245–0.304; P = 0.0005). On APOLLO (7,196 genes), it was 0.747 (0.734–0.758) against 0.492 (0.473–0.510) (P2; difference 0.255, 0.238–0.273; P = 0.0005, the smallest value attainable with 2,000 resamples).


In a pre-specified descriptive analysis of APOLLO, the forecast identified genes with an observed concordance of at least 0.5 with an area under the receiver operating characteristic curve (AUROC) of 0.868, against 0.753 for the RNA-variance baseline and 0.857 for the CPTAC mean lookup (Supplementary Table S5).


### Comparison with looking up observed correlations

The preregistered baseline was deliberately simple. A user with access to CPTAC would more likely look up the observed correlation, preferably from the same cancer type. We therefore compared the forecast with lookups in post hoc analyses whose definitions were committed before any lookup value was computed (Methods). The cascade lookup, L*, takes for each gene and context the first available of four values: the CPTAC correlation in the matched cancer type (L1), the mean of the gene's CPTAC correlations (L2), the mean ProCan correlation (L3; not used when ProCan itself was evaluated) and the global mean (L4).


The forecast was more accurate than the cascade lookup on all three platforms (Table 3, Fig. 2b). The difference in Spearman correlation with observed concordance was 0.069 (0.021–0.129) on RPPA, 0.029 (0.011–0.044) on ProCan and 0.021 (0.011–0.032) on APOLLO. Among pairs with a matched-cancer CPTAC correlation (L1), the forecast was more accurate on ProCan (0.026, 0.002–0.050) and APOLLO (0.018, 0.008–0.028), and the interval included zero on RPPA (0.026, −0.025 to 0.074). Using the matched correlation wherever it existed and the forecast elsewhere, a hybrid rule, did not improve on the forecast alone: the differences were −0.008 (−0.024 to 0.009) on RPPA, −0.009 (−0.022 to 0.004) on ProCan and −0.016 (−0.026 to −0.007) on APOLLO. Within RPPA cancer types, the forecast ranked antibodies better than the cascade lookup in 28 of 31 cancer types (median difference 0.080, 0.057–0.098; two-sided sign test P = 4.6 × 10^−6). AUROC differences pointed in the same direction (Supplementary Table S3).


Because a matched-cancer correlation is estimated in a single cohort, we also compared the forecast with the mean of the gene's CPTAC correlations across cancer types, the pooled lookup, on the pairs that had a matched-cancer correlation (L1 pairs). The pooled lookup was more accurate than the matched-cancer correlation on APOLLO (0.028, 0.019–0.037) and not clearly different from it on RPPA (−0.032, −0.085 to 0.012) and ProCan (0.005, −0.017 to 0.026). The forecast was more accurate than the pooled lookup on RPPA (0.058, 0.019–0.105) and ProCan (0.021, 0.009–0.032) and less accurate on APOLLO (−0.010, −0.016 to −0.003) (Table 3).

The APOLLO deficit came from 600 of the 6,982 L1 pairs, all for genes outside the P set. No CPTAC correlation had been stored for these genes at development, so the forecast's gene-level features took their default values, whereas the lookups used correlations recomputed later for all genes quantified in each CPTAC cohort (Methods). On these pairs, the forecast was less accurate than the pooled lookup (−0.122, −0.188 to −0.055). On the other 6,382 pairs, the forecast and the pooled lookup were not clearly different (0.001, −0.004 to 0.006), and both were more accurate than the matched-cancer correlation (forecast 0.031, 0.021–0.040; pooled lookup 0.030, 0.021–0.039). On RPPA and ProCan, all L1 pairs were of the second kind. This split was defined after the pooled-lookup results had been seen. Against the CPTAC mean lookup, which falls back to a global mean for genes without stored correlations, the forecast gained most outside the P set: 0.139 (0.070–0.209) on ProCan and 0.276 (0.220–0.328) on APOLLO, against 0.002 (−0.003 to 0.007) inside the P set on APOLLO (Supplementary Table S2).


### Calibrating intervals on a new platform

The forecast comes with 80% and 90% intervals whose half-widths, 0.202 and 0.259, are conformal quantiles of absolute residuals on the pooled CPTAC holdout rows (Methods). Forecast errors shifted between platforms: the mean of observed minus predicted concordance was −0.184 on RPPA lung adenocarcinoma, −0.026 on the ProCan lung lineage and −0.023 on APOLLO, against 0.047 on the CPTAC lung adenocarcinoma holdout. Without adjustment, the 80% intervals covered 0.511 of RPPA pairs, 0.628 of ProCan pairs and 0.833 of APOLLO genes (Fig. 3).


We therefore evaluated two adjustments that use K genes measured on the new platform. The first shifts every forecast by δ, the mean of observed minus predicted concordance in the K genes. The second applies the same shift and replaces the interval half-width with the 80% quantile of absolute errors in the K genes. The shift alone did not bring RPPA or ProCan coverage into the range 0.75–0.85. With the shift and quantile, mean coverage was inside this range on all three platforms at K = 50 (RPPA 0.798, ProCan 0.789, APOLLO 0.779) and K = 100 (0.801, 0.792, 0.791), but not at K = 20 on APOLLO (0.742) (Fig. 3, Supplementary Table S6). The RPPA and ProCan draws were post hoc and the APOLLO draws pre-specified descriptive, so the recommendation of at least 50 calibration genes is post hoc.


### An anchor-residual protein predictor

The preregistrations also tested an anchor-residual (AR) predictor of protein abundance, which combines a gene's own mRNA, the anchor, with a ridge model of the remaining protein variation built from the mRNA of 50 other genes (Methods). Predictors of this kind are not new<sup>12</sup><sup>13</sup>; we tested this one externally because its accuracy determines whether a protein estimate can improve on the gene's own mRNA in cohorts without proteomes. The comparison was the correlation of the AR prediction with protein (ρ_AR) against the correlation of the gene's own mRNA with protein (ρ_S0). In CPTAC development, the mean difference ρ_AR − ρ_S0 ranged from 0.085 to 0.121 across held-out cancer types (Fig. 1b).


Externally, AR improved on the gene's own mRNA on RPPA (H1; mean difference 0.034, 0.024–0.043; 123 antibodies; P = 1.5 × 10^−13) and on APOLLO (P1; 0.085, 0.082–0.088; 7,029 genes), and it outperformed two anchor-free ridge models on APOLLO (Q1, 0.023; Q2, 0.063; 6,322 genes) (Table 2, Fig. 4a). It did not improve on the gene's own mRNA in ProCan cell lines (H2; 0.001, −0.002 to 0.004; 2,623 proteins; P = 0.10). In a post hoc description across 34 contexts (31 RPPA cancer types and 3 ProCan lineages), the median AR gain in a context increased with the similarity of its RNA profile to CPTAC (Spearman correlation 0.521; Fig. 4b). APOLLO's similarity, 0.925, lay inside the range of the CPTAC holdout cancer types (0.883–0.979). The ProCan haematopoietic and lymphoid lineage had the lowest similarity (0.808) and a median gain of −0.003.


### Concordance forecasts for ADC targets

We computed forecasts for 29 ADC targets in each of 33 TCGA cancer types from TCGA RNA (957 forecasts; post hoc; Fig. 5a). The same procedure gave forecasts for 19,211 genes, with the RNA-profile similarity of each cancer type to CPTAC as a scope indicator (Data availability). Where the target protein had been measured externally, forecasts agreed with observed concordance with Spearman correlations of 0.308 on RPPA (124 target–cancer pairs), 0.584 on APOLLO (21 targets) and 0.829 on ProCan (6 target–lineage pairs; the frozen lineage forecasts). Twenty of the 21 APOLLO targets fell inside their 80% interval (Fig. 5b). The exception was CLDN18 (observed concordance 0.388, forecast 0.803). The clinically targeted isoform, CLDN18.2, cannot be separated from other CLDN18 isoforms by gene-level RNA or protein measurements. Eight approved ADCs act on six of the 29 targets (Supplementary Table S7). A search of US Food and Drug Administration labels found no approved antibody–drug conjugate acting on the other 23 targets. The same search found seven approved conjugates against targets outside this set; forecasts for those targets are in the released genome-wide table.


## Discussion

The central result is that the concordance between a gene's mRNA and its protein can be forecast in cohorts where protein was not measured, and that the forecast transferred to three platforms that differ from the training data in technology and sample type. The preregistered tests compared the forecast with an RNA-variance baseline. The post hoc comparisons added a more demanding benchmark, looking up observed correlations. Where the forecast had been given a gene's CPTAC data, it was at least as accurate as these lookups on all three platforms; where a lookup drew on CPTAC correlations that the forecast's features had not used, the lookup was better.

Why can a forecast compete with a measurement? A matched-cancer correlation is a measurement in a different cohort. CPTAC cohorts have 82 to 121 tumors per cancer type, so each correlation carries sampling error, and cohorts differ in sample handling and assays. Averaging a gene's correlations over cancer types reduces this error, and on APOLLO the pooled lookup was more accurate than the matched correlation. The forecast adds the gene's RNA distribution in the query cohort. This matters because a correlation is limited by how much its variables vary in the sample: a gene whose RNA varies little across the tumors of a cohort cannot have high concordance there, whatever its concordance elsewhere. On APOLLO, a TMT-labelled tumor cohort like the CPTAC cohorts, the forecast and the pooled lookup were equally accurate for the genes the forecast had CPTAC data on. On RPPA and on cell-line proteomes, the forecast was more accurate than the pooled lookup, which suggests that the query-side information helps most when the platform or sample type differs from the training data.


These results suggest rules of practice. They rest on post hoc analyses and should be read with that status. First, where CPTAC has measured a gene, its correlations pooled across cancer types are a better reference than the matched cancer type alone. If the forecast was built with the gene's CPTAC data, it was at least as accurate as the pooled value; if not, the pooled value should be preferred. The released forecast table gives the pooled correlation and marks the genes whose CPTAC data the forecast used. Second, for genes without CPTAC correlations, the forecast with its interval is the available estimate, and it gained most over a lookup that falls back to a global average. Third, on a new platform, intervals should be recalibrated with at least 50 genes measured on that platform; uncalibrated intervals undercovered on RPPA and ProCan. Antibody-based assays deserve particular caution: on RPPA lung adenocarcinoma, observed concordance was on average 0.184 lower than forecast, in line with the known variability of antibody quality<sup>9</sup>.


Concordance describes ranking across tumors. High concordance means that RNA orders tumors much as protein would, which is what matters when RNA is used to select tumors or patients with high target expression. It does not say whether a protein is abundant in absolute terms or present at the cell surface, both of which also matter for ADC targets. For CLDN18, the gene-level forecast was high, but the observed concordance in APOLLO was lower, and the clinically relevant isoform cannot be separated by gene-level measurements. Isoform-specific targets need isoform-level data.

The AR predictor was a secondary component. It improved on the gene's own mRNA in the two tumor cohorts but not in ProCan cell lines, and a post hoc description associated smaller gains with lower similarity to CPTAC. The failed test H2 is a preregistered result; the similarity description is post hoc and does not retest it. Protein predictors trained on tumors should not be assumed to transfer to cell lines.

The study has limitations. The comparisons with lookups, the split of the APOLLO lookup comparison by the forecast's CPTAC coverage and the calibration recommendation are post hoc. The forecast's gene-level features covered only the 7,449 genes with stored CPTAC correlations, the 7,432 of the P set and 17 further antibody–drug conjugate targets, so for most of the 19,211 genes in the released table the forecast rests on query-side RNA alone; retraining with all CPTAC-quantified genes may close the gap seen for the other genes, but we did not test this. The preregistration was internal: plans and frozen predictions were committed to version control before unlocking, but the commit times were recorded by the analysis workstation and are not certified by a third party (Supplementary Note 3). The two APOLLO families were controlled separately, so the familywise error rate across them could reach about 0.10. Within an antibody, ranking cancer types by the forecast was not clearly better than ranking them by RNA variance (median difference 0.024, −0.016 to 0.096); the forecast is better suited to comparing genes than to comparing cancer types for one gene. Training used 10 CPTAC cancer types, and contexts whose RNA profiles differ from CPTAC may lie outside the forecast's scope. APOLLO protein values were imputed for proteins quantified in at least half of the samples; a post hoc analysis of the non-imputed matrix gave similar results (Supplementary Note 5). Patient overlap between APOLLO and CPTAC lung adenocarcinoma was excluded only by public identifiers.


In summary, a concordance forecast trained on CPTAC tumors ranked mRNA–protein concordance on three external platforms better than an RNA-variance baseline in preregistered tests, and in post hoc comparisons it was at least as accurate as looking up observed CPTAC correlations for genes whose CPTAC data it had been given. With a short local calibration, it gives a practical estimate of how far mRNA can stand in for protein in cohorts without proteomes.

## Methods

### Data

**CPTAC.** Paired tumor RNA and protein for 10 cancer types were obtained from the CPTAC pan-cancer data freeze v1.2<sup>2</sup> as distributed for the LinkedOmicsKB resource<sup>10</sup>: breast cancer (BRCA; 121 tumors)<sup>19</sup>, clear cell renal cell carcinoma (CCRCC; 103)<sup>20</sup>, colon adenocarcinoma (COAD; 96)<sup>21</sup>, glioblastoma (GBM; 99)<sup>22</sup>, head and neck squamous cell carcinoma (HNSCC; 108)<sup>23</sup>, lung squamous cell carcinoma (LSCC; 108)<sup>24</sup>, lung adenocarcinoma (LUAD; 110)<sup>25</sup>, ovarian cancer (OV; 82)<sup>26</sup>, pancreatic ductal adenocarcinoma (PDAC; 105)<sup>27</sup> and uterine corpus endometrial carcinoma (UCEC; 95)<sup>28</sup>, 1,027 tumors in total. RNA values are RSEM coding upper-quartile log2 abundances, and protein values are log2 ratios to a reference intensity. The P gene set contains the 7,432 genes measured in at least 8 of the 10 cancer types at a tumor measurement rate of at least 0.70 and present in the RNA gene universe.


**TCGA RPPA.** RPPA data for TCGA tumors were the PanCanAtlas Xena matrix TCGA-RPPA-pancan-clean<sup>15</sup>, with 7,754 samples and 258 antibodies. A TCPA data level is not recorded for this file. RNA was the PanCanAtlas batch-corrected RNA-seq matrix, restricted to samples present in both matrices after duplicate RPPA sample identifiers were dropped. Antibody names containing `_p`, `pS`, `pT`, `pY`, `cleaved` or `_c` were excluded, and an antibody mapped to more than one gene symbol was dropped. Concordance was evaluated for 4,045 antibody–cancer pairs (134 antibodies, 31 cancer types) with at least 30 samples and a defined correlation.


**ProCan.** DIA proteomes of cancer cell lines were the averaged ProCan-DepMapSanger protein matrix of 8,498 proteins<sup>16</sup> (figshare file 34411172; PRIDE PXD030304). RNA was DepMap 24Q4 log2(TPM + 1)<sup>29</sup><sup>30</sup>, joined to Sanger cell lines by model identifier. Three lineages with at least 30 cell lines were evaluated: large intestine (31), lung (41), and haematopoietic and lymphoid (74), giving 6,118 protein–lineage pairs and 2,746 proteins.

**APOLLO.** Protein (TMT log2 ratios to a common reference, proteins quantified in at least 50% of samples with k-nearest-neighbor imputation) and RNA (TPM) for 87 resected lung adenocarcinomas<sup>17</sup> were obtained from the open files of the Genomic Data Commons (ProteomeXchange PXD036025). File names, processing rules and hashes are given in Supplementary Note 1.

**Other data.** Cancer Cell Line Encyclopedia (CCLE) RNA<sup>30</sup><sup>29</sup> and proteomes<sup>31</sup> of lung (74 cell lines) and skin (33 cell lines) lineages were used during development only, to compare the anchor-residual predictor with the gene's own mRNA in cell lines. TCGA PanCanAtlas RNA<sup>1</sup> was used to compute forecasts in 33 TCGA cancer types. Human Protein Atlas staining summaries<sup>32</sup> were used for a descriptive comparison (Supplementary Note 7). Sources, file names and SHA-256 hashes of all downloaded files are listed in Supplementary Table S9.


### Concordance

For a gene in a context, concordance ρ was the Spearman correlation between RNA and protein across samples in which both were finite, with average ranks for ties. CPTAC RNA was the stored RSEM log2 abundance and CPTAC protein the stored log2 ratio. APOLLO, ProCan and CCLE RNA were log2(TPM + 1). TCGA RNA was the stored PanCanAtlas matrix. Concordance was computed when at least 30 such samples were available and was undefined when either vector was constant.

### Forecast models and baselines

The training data were CPTAC leave-one-cancer-out rows: for each held-out cancer type, the target was ρ in that cancer type, gene-level features were computed from the other nine cancer types and query features from the RNA of the held-out cancer type. The nine features of F1s were the mean, standard deviation, minimum and maximum of the gene's CPTAC correlations in the training cancers; the number of those cancers with a defined correlation; the mean, across training cancers, of the percentile of the gene's RNA mean; the mean fraction of paired tumors in which the protein was measured; and the percentiles of the gene's RNA mean and RNA standard deviation in the query context. F1, fitted in the preceding stage, used the five correlation summaries, the query RNA mean and standard deviation and the percentile of that standard deviation, the mean of the training RNA means, and the training protein measurement rate. Each model was fitted as a ridge regression and as LightGBM<sup>18</sup> (300 trees, learning rate 0.05, 31 leaves, minimum 20 samples per leaf, row subsample 0.8 at every iteration, column subsample 0.8, random state 20261001); F1s-LightGBM was the preregistered forecast. Missing feature values were replaced by the column median of the training matrix. The model applied to external data was fitted to the stacked leave-one-cancer-out rows for P genes with a finite observed concordance, not to a single matrix that used all ten cancer types in every feature. The 80% and 90% interval half-widths, 0.202 and 0.259, are the conformal quantiles of absolute residuals on those rows: the residual at position ceil((n+1)(1−α)) in the ascending absolute residuals, with α = 0.20 and 0.10. Intervals were not clipped to [−1, 1]. The forecast model was fitted twice with the same training rows, code and parameters: before the RPPA and ProCan predictions were frozen and before the APOLLO predictions were frozen. The released model reproduces the APOLLO predictions exactly but not the RPPA and ProCan predictions (maximum absolute differences 0.031 and 0.044; Spearman correlations 0.9991 and 0.9992), and the cause was not identified (Supplementary Note 2). All reported results use the frozen predictions, which are released.


The RNA-variance baseline was the percentile of the gene's RNA standard deviation among genes in the query context. The CPTAC mean lookup F0 was the mean of the gene's defined CPTAC correlations stored at development; genes without any stored CPTAC correlation received the global mean of P-gene correlations, 0.456 (72,880 gene–cancer pairs).


### Anchor-residual predictor

The anchor is the gene's own mRNA. Its slope on protein is the least-squares slope without an intercept on measured samples, set to zero when negative or undefined. The residual is protein minus that slope times the gene's own RNA, after within-context standardization. The 50 residual genes are the other RNA genes with the largest absolute Pearson correlation to that residual; the gene itself is excluded, and ties break toward the smaller gene index. A ridge regression without an intercept is fitted to the residual on those 50 genes. The penalty α is chosen from 1, 10, 100, 1,000 and 10,000, and the residual weight λ from 0, 0.25, 0.5, 0.75 and 1, by the mean Spearman correlation in inner leave-one-cancer-out folds; ties take the smaller λ and then the larger α. The prediction is the anchor term plus λ times the residual prediction. The anchor-free baselines are R1, a ridge model on the gene's own RNA plus 50 other genes selected by absolute correlation with protein, and R2, a ridge model on all RNA genes. The R2 penalty used for the external comparison was 1,000. The own-mRNA baseline S0 is the gene's own RNA, so ρ_S0 equals the observed concordance. For APOLLO, predictions used the stored final model without refitting.

### Preregistration and confirmatory tests

The tests on RPPA and ProCan (H1–H4) and on APOLLO (P1, P2, Q1, Q2) were preregistered in two analysis plans committed to version control on 2 October 2026, each before the protein data of the platforms concerned were used. Predictions were frozen and committed before an unlock commit, and the confirmation scripts were run once. Commit hashes and times, power analyses, tests on simulated planted and null data, and the departures accepted before unlocking are reported in Supplementary Notes 3 and 4.

H1, H2, P1, Q1 and Q2 compared unit-level differences in correlation with protein (ρ_AR − ρ_S0, ρ_AR − ρ_R1 or ρ_AR − ρ_R2) with a one-sided Wilcoxon signed-rank test<sup>33</sup>; units were antibodies for RPPA, proteins for ProCan and genes for APOLLO. For H1 and H2 the tested value was the mean of the context-level differences within each antibody or protein. H3, H4 and P2 compared the Spearman correlation of the forecast with observed concordance against that of the RNA-variance baseline. The P value was one plus the number of 2,000 unit resamples in which the difference was at or below zero, divided by 2,001. Confidence intervals were the 2.5th and 97.5th percentiles of 2,000 unit resamples. Holm's procedure<sup>40</sup> controlled the familywise error rate at 0.05 within each of three families: H1–H4, with ordered thresholds 0.0125, 0.05/3, 0.025 and 0.05; and P1 and P2, and Q1 and Q2, each with ordered thresholds 0.025 and 0.05. Each test restarted the generator with seed 20261001.

### Post hoc analyses

All analyses in this section were specified after the confirmations. The lookup comparisons (cascade lookup, hybrid rule and pooled lookup) were defined in plans committed before any lookup value was computed (Supplementary Note 3).

The cascade lookup L* and the hybrid rule were defined as in the Results, with the context map in Supplementary Table S10. L1 and L2 were computed from CPTAC correlations for all genes quantified in each cohort, which added 9,318 gene–cohort correlations to those stored at development; F0 used only the stored correlations. A difference involving a constant predictor was undefined. The pooled lookup was the mean of the gene's correlations across CPTAC cancer types in this recomputed set (level L2 of the cascade), including the matched cancer type, and it was compared with the matched correlation and with the forecast on pairs that had a matched correlation (L1 pairs). After these comparisons had been seen, the APOLLO L1 pairs were split by whether F0 was gene-specific or equal to the global mean, that is, whether the forecast's gene-level features had been computed from stored CPTAC correlations; this split was committed before calculation (Supplementary Note 3).

Differences in Spearman correlation and AUROC (positive class: observed ρ of at least 0.5) were estimated with 2,000 resamples of units (antibodies, proteins or genes) with the generator restarted at seed 20261001 for every contrast; a resample kept all pairs of each sampled unit, and intervals were the 2.5th and 97.5th percentiles of resamples with a defined difference. Within-context comparisons used a two-sided sign test across RPPA cancer types.

### Calibration

For each platform, units with finite observed and predicted concordance were shuffled with seed 20261014 and split into five folds. Each fold was evaluated in turn. From the units in the other folds, K units were drawn without replacement, and the draw was repeated 20 times, giving 100 draws for each K of 10, 20, 50 and 100. The shift δ was the mean of observed minus predicted concordance in the K units. The shift-and-quantile method also replaced the half-width by the 0.80 (or 0.90) quantile of absolute errors in the K units, an approach related to split-conformal calibration<sup>34</sup>. Coverage was the fraction of evaluation units whose observed concordance lay inside the interval, averaged over draws. A method was recommended only if its mean 80% coverage lay within 0.75–0.85, inclusive, on all three platforms.

### RNA-profile similarity

Similarity s is the maximum Spearman correlation, over the 10 CPTAC cancer types, between the context's gene-wise RNA percentile profile and that cancer type's profile. Both profiles are percentiles of gene-wise RNA medians on the same gene universe. The reference range was the range of s over the CPTAC holdout cancer types (0.883–0.979).


### ADC target atlas

The 29 targets were fixed in the project's first task instruction and committed on 1 October 2026, before any analysis of these targets was run; no formal selection rule was recorded. For each target and TCGA cancer type, the forecast and its 80% interval were computed from TCGA PanCanAtlas RNA. Observed concordance was compared where the target protein had been measured. Approved ADCs were identified as of 2 October 2026 by searching US Food and Drug Administration labels through the openFDA interface for antibody–drug conjugates and reading the mechanism-of-action section of each label (Supplementary Note 7, Supplementary Table S7).

### Software

Analyses used Python 3.12.4 with NumPy 1.26.4<sup>35</sup>, pandas 2.2.2<sup>36</sup>, SciPy 1.13.1<sup>37</sup>, scikit-learn 1.4.2<sup>38</sup>, LightGBM 4.7.0<sup>18</sup> and Matplotlib 3.8.4<sup>39</sup> on a CPU server without a GPU. These versions are those recorded for this virtual environment before the 2 October 2026 freeze.

### Use of large language models

Large language models were used throughout this study under the author's direction. Claude Opus 5.5 (Anthropic) drafted the analysis plans, the stage-by-stage task instructions and the text of this manuscript. A coding agent in the Cursor editor running Grok 4.7 (xAI) wrote and ran the analysis code, produced the figures and earlier drafts of the text, and ran automated checks of the references against Crossref and of every number in the manuscript against the result files. The author reviewed each analysis plan before it was committed and executed; preregistered predictions were frozen before the external protein data were unlocked. The author reviewed all plans, code outputs and text and takes full responsibility for the content. Language models are not authors of this work.

## Data availability

All data analysed in this study are public. CPTAC data were obtained from LinkedOmicsKB (https://kb.linkedomics.org/), TCGA RPPA data from the PanCanAtlas Xena matrix (https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/TCGA-RPPA-pancan-clean.xena.gz), ProCan proteomes from figshare file 34411172 (PRIDE PXD030304), APOLLO data from the Genomic Data Commons (https://gdc.cancer.gov/about-data/publications/APOLLO-LUAD-2022; ProteomeXchange PXD036025), CCLE proteomes from https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz, DepMap 24Q4 data from https://doi.org/10.25452/figshare.plus.27993248.v1 and TCGA PanCanAtlas RNA from https://gdc.cancer.gov/about-data/publications/pancanatlas. File names, sources and SHA-256 hashes are listed in Supplementary Table S9. Code, frozen predictions, the trained forecast model, forecast tables for all genes and for the ADC targets, and figure source data are available at [AUTHOR TO FILL: GitHub URL] and archived at Zenodo ([AUTHOR TO FILL: DOI]). For the 11,762 of 19,211 genes without stored CPTAC correlations, the gene-level features take their default values, so the forecast rests on the gene's RNA distribution in the query cohort; the released table marks these genes.

## References

1. Hoadley, K. A. et al. Cell-of-Origin patterns dominate the molecular classification of 10,000 tumors from 33 types of cancer. *Cell* **173**, 291–304.e6 (2018).
2. Li, Y. et al. Proteogenomic data and resources for pan-cancer analysis. *Cancer Cell* **41**, 1397–1406 (2023).
3. Drago, J. Z., Modi, S. & Chandarlapaty, S. Unlocking the potential of antibody–drug conjugates for cancer therapy. *Nat. Rev. Clin. Oncol.* **18**, 327–344 (2021).
4. Bosi, C. et al. Pan-cancer analysis of antibody-drug conjugate targets and putative predictors of treatment response. *Eur. J. Cancer* **195**, 113379 (2023).
5. Savage, S. R. et al. Pan-cancer proteogenomics expands the landscape of therapeutic targets. *Cell* **187**, 4389–4407.e15 (2024).
6. Vogel, C. & Marcotte, E. M. Insights into the regulation of protein abundance from proteomic and transcriptomic analyses. *Nat. Rev. Genet.* **13**, 227–232 (2012).
7. Liu, Y., Beyer, A. & Aebersold, R. On the dependency of cellular protein levels on mRNA abundance. *Cell* **165**, 535–550 (2016).
8. Upadhya, S. R. & Ryan, C. J. Experimental reproducibility limits the correlation between mRNA and protein abundances in tumor proteomic profiles. *Cell Rep. Methods* **2**, 100288 (2022).
9. Upadhya, S. R. & Ryan, C. J. Antibody reliability influences observed mRNA–protein correlations in tumour samples. *Life Sci. Alliance* **6**, e202201885 (2023).
10. Liao, Y. et al. A proteogenomics data-driven knowledge base of human cancer. *Cell Syst.* **14**, 777–787.e5 (2023).
11. Nawaz, U. et al. OnCorr: a pan-cancer mRNA-protein correlation tool for precision oncology. *npj Precis. Oncol.* **10**, 128 (2026).
12. Yang, M. et al. Community assessment of the predictability of cancer protein and phosphoprotein levels from genomics and transcriptomics. *Cell Syst.* **11**, 186–195.e9 (2020).
13. Li, H., Siddiqui, O., Zhang, H. & Guan, Y. Joint learning improves protein abundance prediction in cancers. *BMC Biol.* **17**, 107 (2019).
14. Srivastava, H. et al. Protein prediction models support widespread post-transcriptional regulation of protein abundance by interacting partners. *PLoS Comput. Biol.* **18**, e1010702 (2022).
15. Li, J. et al. TCPA: a resource for cancer functional proteomics data. *Nat. Methods* **10**, 1046–1047 (2013).
16. Gonçalves, E. et al. Pan-cancer proteomic map of 949 human cell lines. *Cancer Cell* **40**, 835–849.e8 (2022).
17. Soltis, A. R. et al. Proteogenomic analysis of lung adenocarcinoma reveals tumor heterogeneity, survival determinants, and therapeutically relevant pathways. *Cell Rep. Med.* **3**, 100819 (2022).
18. Ke, G. et al. LightGBM: a highly efficient gradient boosting decision tree. Advances in Neural Information Processing Systems 30 (2017).
19. Krug, K. et al. Proteogenomic landscape of breast cancer tumorigenesis and targeted therapy. *Cell* **183**, 1436–1456.e31 (2020).
20. Clark, D. J. et al. Integrated proteogenomic characterization of clear cell renal cell carcinoma. *Cell* **179**, 964–983.e31 (2019).
21. Vasaikar, S. et al. Proteogenomic analysis of human colon cancer reveals new therapeutic opportunities. *Cell* **177**, 1035–1049.e19 (2019).
22. Wang, L. B. et al. Proteogenomic and metabolomic characterization of human glioblastoma. *Cancer Cell* **39**, 509–528.e20 (2021).
23. Huang, C. et al. Proteogenomic insights into the biology and treatment of hPV-negative head and neck squamous cell carcinoma. *Cancer Cell* **39**, 361–379.e16 (2021).
24. Satpathy, S. et al. A proteogenomic portrait of lung squamous cell carcinoma. *Cell* **184**, 4348–4371.e40 (2021).
25. Gillette, M. A. et al. Proteogenomic characterization reveals therapeutic vulnerabilities in lung adenocarcinoma. *Cell* **182**, 200–225.e35 (2020).
26. McDermott, J. E. et al. Proteogenomic characterization of ovarian HGSC implicates mitotic kinases, replication stress in observed chromosomal instability. *Cell Rep. Med.* **1**, 100004 (2020).
27. Cao, L. et al. Proteogenomic characterization of pancreatic ductal adenocarcinoma. *Cell* **184**, 5031–5052.e26 (2021).
28. Dou, Y. et al. Proteogenomic characterization of endometrial carcinoma. *Cell* **180**, 729–748.e26 (2020).
29. Broad DepMap. DepMap 24Q4 Public. figshare+ https://doi.org/10.25452/figshare.plus.27993248.v1 (2024).
30. Ghandi, M. et al. Next-generation characterization of the cancer cell line encyclopedia. *Nature* **569**, 503–508 (2019).
31. Nusinow, D. P. et al. Quantitative proteomics of the cancer cell line encyclopedia. *Cell* **180**, 387–402.e16 (2020).
32. Uhlen, M. et al. A pathology atlas of the human cancer transcriptome. *Science* **357**, eaan2507 (2017).
33. Wilcoxon, F. Individual comparisons by ranking methods. *Biometrics Bull.* **1**, 80 (1945).
34. Lei, J., G’Sell, M., Rinaldo, A., Tibshirani, R. J. & Wasserman, L. Distribution-free predictive inference for regression. *J. Am. Stat. Assoc.* **113**, 1094–1111 (2018).
35. Harris, C. R. et al. Array programming with NumPy. *Nature* **585**, 357–362 (2020).
36. McKinney, W. Data structures for statistical computing in Python. *Proceedings of the Python in Science Conference* 56–61 (2010).
37. Virtanen, P. et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nat. Methods* **17**, 261–272 (2020).
38. Pedregosa, F. et al. Scikit-learn: machine learning in Python. J. Mach. Learn. Res. 12, 2825–2830 (2011).
39. Hunter, J. D. Matplotlib: a 2D graphics environment. *Comput. Sci. Eng.* **9**, 90–95 (2007).
40. Holm, S. A simple sequentially rejective multiple test procedure. *Scand. J. Stat.* **6**, 65–70 (1979).

## Acknowledgements

The results shown here are in whole or part based upon data generated by the TCGA Research Network: https://www.cancer.gov/tcga. Data used in this publication were generated by the Clinical Proteomic Tumor Analysis Consortium (NCI/NIH). This research received no specific grant from any funding agency in the public, commercial or not-for-profit sectors.

## Author contributions

J.K. designed the study, wrote the analysis plans and the preregistrations, directed and reviewed the computational work, analysed the results and wrote the manuscript. J.C. supervised the study, reviewed the analysis plans and the manuscript, and provided critical revision. Both authors approved the submitted version.

## Competing interests

The authors declare no competing interests.

## Figure legends

**Figure 1. Development of the concordance forecast and the anchor-residual predictor on CPTAC tumors.** (a) Spearman correlation between predicted and observed concordance in each held-out cancer type in leave-one-cancer-out analysis, for the forecast (F1s-LightGBM) and the RNA-variance baseline. (b) Mean difference between the correlation of the anchor-residual prediction with protein and that of the gene's own mRNA (ρ_AR − ρ_S0) in each held-out cancer type; genes per cancer type: 6,778–7,426. Cancer types are ordered by the forecast correlation in (a). The bars in (b) are grey. Development analyses, not preregistered tests.


**Figure 2. External evaluation of the concordance forecast.** (a) Spearman correlation between predicted and observed concordance on TCGA RPPA (4,045 antibody–cancer pairs; 134 antibodies), ProCan cell lines (6,118 protein–lineage pairs; 2,746 proteins) and APOLLO lung adenocarcinoma (7,196 genes). The series are Forecast (preregistered), RNA variance (preregistered), CPTAC mean lookup (post hoc) and Cascade lookup (post hoc), with the markers shown in the panel. The preregistered series have 95% intervals. Lookup series are point estimates. (b) Preregistered differences in Spearman correlation with 95% intervals: forecast minus RNA-variance baseline (H3, H4 and P2). (c) Post hoc differences: forecast minus cascade lookup on all pairs and on pairs with a matched-cancer CPTAC correlation, and hybrid rule minus forecast. Intervals are 2.5th–97.5th percentiles of 2,000 resamples of antibodies, proteins or genes.

**Figure 3. Calibration of 80% forecast intervals with locally measured genes.** Mean coverage of nominal 80% intervals over 100 draws (five folds by 20 repeats) with K = 10, 20, 50 or 100 calibration genes, for no adjustment, a shift by the mean error δ, and the shift plus an empirical quantile, on RPPA (a) and ProCan (b), both post hoc, and APOLLO (c; pre-specified descriptive). The panel legend lists the three methods. On the APOLLO panel the unadjusted series is marked already within the band. Error bars, 2.5th–97.5th percentiles of draws; grey band, coverage 0.75–0.85; dashed line, 0.80.

**Figure 4. External tests of the anchor-residual predictor.** (a) Mean unit-level differences with 95% intervals for the preregistered tests H1 (RPPA, 123 antibodies), H2 (ProCan, 2,623 proteins), P1 (APOLLO, 7,029 genes), and Q1 and Q2 (APOLLO, 6,322 genes). (b) Grey squares: post hoc comparisons with the anchor-free baselines on RPPA (100 antibodies) and ProCan (2,399 proteins), without intervals. (c) Median difference ρ_AR − ρ_S0 per context against RNA-profile similarity to CPTAC, for 31 RPPA cancer types and 3 ProCan lineages (post hoc; circles) and APOLLO (preregistered; diamond). A horizontal line marks zero. Dashed lines mark the range of similarity across CPTAC holdout cancer types (0.883–0.979).


**Figure 5. Concordance forecasts for antibody–drug conjugate targets.** (a) Forecast concordance for 29 ADC targets in 33 TCGA cancer types (post hoc). The row label is CLDN18* and the footnote states that gene-level assays do not separate CLDN18.2. The colour scale spans the plotted forecasts. (b) Forecast against observed concordance for targets whose protein was measured: RPPA (124 target–cancer pairs; Spearman correlation 0.308), ProCan frozen lineage forecasts (6 pairs; 0.829) and APOLLO (21 targets; 0.584), with CPTAC leave-one-cancer-out predictions (228 pairs; 0.554) for reference. A separate ProCan series pairs the TCGA-disease forecast with the ProCan measurement (12 pairs; 0.793) and is not shown. Open circle, CLDN18 in APOLLO. The diagonal is marked y = x.


## Tables

**Table 1. Data sets.**

| Data set | Platform | Samples | Units analysed | Role |
|---|---|---|---|---|
| CPTAC, 10 cancer types | RNA-seq; mass spectrometry | 1,027 tumors | 7,432 genes (P set) | Training and development |
| CCLE lung and skin lineages | RNA-seq; TMT mass spectrometry | 74 and 33 cell lines | 6,911 and 5,947 genes | Development |
| TCGA RPPA, 31 cancer types | RNA-seq; reverse-phase protein array | 7,754 samples in the matrix | 134 antibodies (4,045 antibody–cancer pairs) | Preregistered tests H1, H3 |
| ProCan, 3 lineages | RNA-seq; DIA mass spectrometry | 146 cell lines in the 3 lineages | 2,746 proteins (6,118 protein–lineage pairs) | Preregistered tests H2, H4 |
| APOLLO lung adenocarcinoma | RNA-seq; TMT mass spectrometry | 87 tumors | 7,196 genes | Preregistered tests P1, P2, Q1, Q2 |
| TCGA PanCanAtlas, 33 cancer types | RNA-seq | 9,871 tumors | 29 ADC targets and 19,211 genes | ADC atlas (post hoc) |



Notes: Units analysed are those of the concordance tests (H3, H4, P2). The anchor-residual tests used 123 RPPA antibodies (H1), 2,623 ProCan proteins (H2), 7,029 APOLLO genes (P1) and 6,322 APOLLO genes (Q1, Q2).

**Table 2. Preregistered tests.**

| Test | Data set | Comparison | Units | Estimate (95% CI) | P | Holm decision |
|---|---|---|---|---|---|---|
| H1 | TCGA RPPA | ρ_AR − ρ_S0, mean | 123 antibodies | 0.0339 (0.0242–0.0432) | 1.5 × 10^−13 | Rejected |
| H2 | ProCan | ρ_AR − ρ_S0, mean | 2,623 proteins | 0.000932 (−0.00236 to 0.00423) | 0.10 | Not rejected |
| H3 | TCGA RPPA | Forecast − RNA-variance baseline, Spearman | 134 antibodies (4,045 pairs) | 0.205 (0.110–0.295) | 0.0005 | Rejected |
| H4 | ProCan | Forecast − RNA-variance baseline, Spearman | 2,746 proteins (6,118 pairs) | 0.276 (0.245–0.304) | 0.0005 | Rejected |
| P1 | APOLLO | ρ_AR − ρ_S0, mean | 7,029 genes | 0.0850 (0.0822–0.0878) | < 10^−300 | Rejected |
| P2 | APOLLO | Forecast − RNA-variance baseline, Spearman | 7,196 genes | 0.255 (0.238–0.273) | 0.0005 | Rejected |
| Q1 | APOLLO | ρ_AR − ρ_R1, mean | 6,322 genes | 0.0228 (0.0212–0.0244) | 4.5 × 10^−203 | Rejected |
| Q2 | APOLLO | ρ_AR − ρ_R2, mean | 6,322 genes | 0.0627 (0.0597–0.0659) | < 10^−300 | Rejected |



Notes: Wilcoxon signed-rank tests (one-sided) for H1, H2, P1, Q1 and Q2; resampling tests of the difference in Spearman correlation for H3, H4 and P2 (2,000 resamples; smallest attainable P = 0.0005). Intervals are 2.5th–97.5th percentiles of 2,000 unit resamples. Holm families: H1–H4; P1 and P2; Q1 and Q2. For P1 and Q2, SciPy returned a P value of 0 because the normal approximation underflowed; those entries are reported as P < 10^−300.

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

Notes: Differences in Spearman correlation with observed concordance; 95% intervals from 2,000 resamples of antibodies (RPPA), proteins (ProCan) or genes (APOLLO). All pairs: 4,045 (134 antibodies), 6,118 (2,746 proteins) and 7,196 genes. L1, pairs with a matched-cancer CPTAC correlation: 1,206 (127 antibodies), 3,356 (2,091 proteins) and 6,982 genes. On L1 pairs, the cascade lookup equals the matched lookup. Covered, pairs whose gene had stored CPTAC correlations, so that the forecast's gene-level features were gene-specific (6,382 APOLLO genes; all RPPA and ProCan L1 pairs); not covered, pairs whose gene had none (600 APOLLO genes, all outside the P set). The split was defined after the L1 results had been seen. The CPTAC mean lookup used the correlations stored at development; the cascade, matched and pooled lookups used a recomputed set (Methods).
