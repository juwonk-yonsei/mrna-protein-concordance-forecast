# Project B release v1.3.0

Local copy of the code, analysis plans, frozen outputs and Scientific Reports files for the mRNA–protein concordance forecast. v1.0.0, v1.1.0 and v1.2.0 are unchanged and are not in this folder.

The code is under the MIT License. Derived result tables and model files are under Creative Commons Attribution 4.0 International. Raw data matrices are not redistributed.

Raw matrices are not in this folder. The manuscript Data availability section lists the public locations. This folder does not redistribute them.

Also absent: sample identifiers, the planted-signal and null checks, and any file over 100 MB. The three large files are listed in `zenodo_only/files.tsv`.

Authors: Juwon Kang and Junjeong Choi. A GitHub URL and a Zenodo DOI are not filled in.

## Forecast table v3

`results/stage12/forecast_all_genes_v3.tsv.gz` adds `model_gene_features` to the v2 table. The other columns are the v2 columns.

Of the 19,211 genes, 11,762 have model_gene_features 0: the gene-level features are the default values, so the forecast rests on the gene's RNA distribution in the query cohort.

If model_gene_features is 0 and cptac_pooled_rho is available, the pooled CPTAC correlation was more accurate than the forecast in a post hoc comparison on APOLLO and should be preferred.

## Forecast model

One model file is released, `release/models/F1s_lightgbm.txt`. A post hoc grid over 1, 16 and 72 threads and over automatic, column-wise and row-wise histograms did not reproduce the frozen RPPA or ProCan predictions within 1e-12. The LightGBM version recorded before the RPPA and ProCan fit and after the APOLLO fit is 4.7.0, so no other version was installed. The cause was not identified, and no second model file was written.

## Commands

Set `B_CPU_CAP=none` before a project script. Redraw the figures from `results/figures/source` and the stored result tables with:

```
B_CPU_CAP=none python scripts/36_stage10_figures.py
```

Scripts, in filename order:

- `scripts/01_download.py`: 
- `scripts/02_cptac_audit.py`: 
- `scripts/03_aux_audit.py`: 
- `scripts/04_models.py`: 
- `scripts/05_maps.py`: Section 2. Antibody names and sample IDs only. No RPPA or ProCan values.
- `scripts/06_develop.py`: Section 3. AR and RNA-proxy reliability. CPTAC only.
- `scripts/07_freeze.py`: Section 4. Freeze RNA-proxy forecasts. Does not read locked values.
- `scripts/08_stage3.py`: Stage 3 development checks and freezes. Does not read locked values.
- `scripts/09_confirm.py`: Confirmation metrics.
- `scripts/10_confirm_B4.py`: Stage-4 confirmation.
- `scripts/11_posthoc_B4.py`: Post-hoc tables. Not part of the preregistered confirmation.
- `scripts/12_stage5_r1r2.py`: Refit stage-1 R1 and R2 on all ten CPTAC cancers and score RPPA and ProCan.
- `scripts/12_stage5_scope.py`: Stage-5 post-hoc sensitivity, H2 decomposition, and training-context distance.
- `scripts/12_stage5_tables.py`: Stage-5 post-hoc tables that use only the stage-4 confirmation pair files.
- `scripts/13_stage5_atlas.py`: ADC target reliability atlas for TCGA primary tumors.
- `scripts/14_stage5_lock.py`: Lock the APOLLO LUAD protein matrix and record the cohort survey.
- `scripts/15_stage6_rna.py`: APOLLO LUAD RNA processing and protein-header unit bounds.
- `scripts/16_stage6_power.py`: Development-only power and CPTAC LUAD reference values for the B6 preregistration.
- `scripts/17_confirm_B6.py`: APOLLO LUAD confirmation.
- `scripts/18_check_B6.py`: Independent check of a B6 confirmation directory.
- `scripts/19_stage6_freeze.py`: Freeze R1/R2 weights and APOLLO predictions. Does not read APOLLO protein values.
- `scripts/20_stage6_posthoc.py`: Post-hoc tables after the APOLLO confirmation output was committed.
- `scripts/21_stage6_raw.py`: Post-hoc sensitivity on the public non-imputed APOLLO matrix.
- `scripts/22_stage6_figures.py`: Draft figures and tables for the APOLLO confirmation.
- `scripts/23_stage7.py`: Post-hoc F0 comparison, calibration table, and ADC atlas final table.
- `scripts/24_stage7_figures.py`: Final figures 3, 5, and 6 at journal size.
- `scripts/25_stage8_lookup.py`: Post-hoc matched-context lookup L* and the HYB rule.
- `scripts/26_check_refs.py`: Record Crossref checks. Writes manuscript/checks only after the queries return.
- `scripts/27_check_numbers.py`: Check the stage-8 headline numbers in the manuscript against the result tables.
- `scripts/28_stage8_figures.py`: Redraw figures 1, 2, and 4 at journal size, and add F1s-L* to figure 3.
- `scripts/29_check_ledger.py`: Check manuscript numbers against manuscript/checks/ledger_B.tsv.
- `scripts/30_stage9_f0_lstar.py`: Post-hoc F0 minus L* on the stored stage-8 lookup table.
- `scripts/31_stage9_export.py`: Refit F1s-LightGBM for a public model file and compare it with frozen predictions.
- `scripts/32_stage10_lpool.py`: Post-hoc pooled lookup versus matched lookup on stored stage-8 rows.
- `scripts/33_stage10_export.py`: Export the public F1s-LightGBM text model from a frozen-parameter refit.
- `scripts/34_stage10_forecast.py`: All-gene TCGA forecasts from the exported F1s model and the ADC-atlas feature code.
- `scripts/35_stage10_diag.py`: Rebuild RPPA and ProCan F1s features with the freeze-commit code.
- `scripts/36_stage10_figures.py`: Draw the Scientific Reports figures from stored result tables.
- `scripts/37_stage10_ledger.py`: Type each number in the Scientific Reports manuscript and check it.
- `scripts/38_stage10_docx.py`: Build the Scientific Reports submission files.
- `scripts/39_stage10_release.py`: Assemble the local release folder and the review pack.
- `scripts/40_stage11_l1_split.py`: Post-hoc split of L1 pairs by whether F0 is the global mean.
- `scripts/41_stage11_forecast_v2.py`: Add observed CPTAC lookup columns to a copy of the all-gene forecast table.
- `scripts/42_stage11_b4_refit.py`: Refit F1s-LightGBM on the training rows stored at the B4 freeze commit.
- `scripts/43_stage12_repro.py`: Post-hoc LightGBM thread and histogram grid against the frozen forecasts.
- `scripts/44_stage12_forecast_v3.py`: Add model_gene_features to the all-gene forecast. No model is fitted.
- `scripts/45_stage12_adc.py`: Post-hoc openFDA and Drugs@FDA search for approved antibody-drug conjugates.
- `scripts/46_stage12_release.py`: Assemble release v1.1.0 and review pack v2. v1.0.0 is left in place.
- `scripts/47_stage13_release.py`: Assemble release v1.2.0 and review pack v3. v1.0.0 and v1.1.0 stay in place.
- `scripts/48_stage14_pack.py`: Build review pack v4. No new release: the training-row comparison is OPTION X and no license or DOI was supplied.
- `scripts/49_stage15_roundtrip.py`: CSV round-trip check for the frozen F1s training table.
- `scripts/50_stage15_pack.py`: Build review pack v5. No new release: the CSV round trip is OPTION Z and no license or DOI was supplied.
- `scripts/51_stage16_release.py`: Assemble release v1.3.0 and review pack v6. Earlier releases stay in place.
- `scripts/__init__.py`: 
- `scripts/b2_lib.py`: Shared stage-2 helpers. Importing this module removes the stage-1 CPU cap.
- `scripts/b3_lib.py`: Stage-3 helpers. Importing this module removes the stage-1 CPU cap.
- `scripts/common.py`: 
- `scripts/loader.py`: Shared loader.

`scripts_freeze/08_stage3.py` is the file from commit `224779714bffd58b02059f51de48d92d1324a24f`.

## Plan commits

```
2cbaafbabcb7328e38322405959e1ede2328d7e0 2026-10-15 20:55:16 +0900 Lock the CSV round-trip check before any fit.
4fdb98c75cc30e4c3a629d955d990586779b7959 2026-10-15 15:37:30 +0900 Commit the stage-14 training-row comparison before any fit.
546f04904132a258bbb0adca014fde1998b640b5 2026-10-15 14:54:49 +0900 Commit the stage-13 version check before any refit.
32d4e7191a140668a27cd716d728e63daa73aa71 2026-10-15 14:07:53 +0900 Commit the stage-12 reproduction grid before any refit.
773bd579f9c2bb4eb65e1e0f0cf955fbbe71a6eb 2026-10-15 11:47:50 +0900 Commit the stage-11 L1 split and forecast-v2 plan before any calculation.
68a32e5974eb0426032c8310d13ec7d716fc27c7 2026-10-14 10:38:30 +0900 Commit the stage-10 Lpool and all-gene forecast plan before any calculation.
87967f098ea1434f39c3cd3109d2edeea1c1a0ef 2026-10-02 15:50:05 +0900 Preregister the APOLLO LUAD confirmation hypotheses.
f7befd8324bcef63d89cebb57de36c0232b5b18f 2026-10-02 09:57:46 +0900 Preregister stage-4 hypotheses before any external baseline or unlock.
7860a43c446c433f3a49563f993970d68afce264 2026-10-01 23:16:00 +0900 Freeze the stage-3 plan and scripts before any F1s or AR fit.
```

Result tables keep their original column headers. `excluded.tsv` lists files that were not copied.
