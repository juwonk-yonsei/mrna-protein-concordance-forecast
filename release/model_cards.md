# Model cards

## F1s-LightGBM

Training data: CPTAC paired tumor RNA and protein, 10 cancers, 72880 gene-by-cancer pairs in the development table. The stored final predictions used the model fit on all development cancers. This stage did not refit it.

Features: rho_mean, rho_sd, rho_min, rho_max, n_defined_cancers, holdout_rna_mean_percentile, holdout_rna_sd_percentile, train_rna_mean_percentile_mean, train_protein_rate_mean.

Hyperparameters, read from the methods record rather than from a saved booster file: n_estimators 300, learning_rate 0.05, num_leaves 31, min_child_samples 20, subsample 0.8, subsample_freq 1, colsample_bytree 0.8, random_state 20261001, n_jobs -1.

Preregistered performance: Spearman correlation with observed rho exceeded B_sd on RPPA (difference 0.205), ProCan (0.276), and APOLLO (0.255).

Calibration recommendation: post hoc, δ plus a quantile, K at least 50, for an 80% interval whose mean coverage sits in 0.75–0.85 on RPPA, ProCan, and APOLLO.

Known limits: the unadjusted 80% interval covers about 0.51 of RPPA pairs. Within one antibody, across cancers, the median gain over B_sd includes zero. Inside P on APOLLO, the gain over F0 includes zero. On the RPPA L1 layer, the gain over the matched CPTAC correlation includes zero.

Intended use: research only. Not a clinical test and not a substitute for a measured correlation when that measurement exists and the question is the observed value.

Export: the fitted booster object was not stored. A text-model export would require a refit. No refit was run, and no text model file was written.

## Anchor-residual

Training data and features are the stored stage-5 ridge on anchors and 50 residual genes. The frozen coefficient file is `results/stage6/frozen/r1r2_model.npz`, SHA-256 `4c1b0f6f7a75cbb91b1043f01005d154e01ec2ba74d28da0402918beccf6d4cc`, 121062358 bytes. It was not refit. A byte copy is `release/models/r1r2_model.npz` when that copy is present. The match test against the frozen file is a SHA-256 comparison, not a new prediction run, because the coefficient array is the model.

Preregistered performance: mean delta 0.0339 on RPPA and 0.0850 on APOLLO. ProCan mean delta 0.000932 was not significant.

Known limits: hematopoietic ProCan lineages. Research use only.
