[AUTHOR TO FILL: submission date]

The Editors
Scientific Reports

Dear Editors,

Please consider our manuscript, "A preregistered model forecasts mRNA–protein concordance in unmeasured cancer contexts across three proteomic platforms", for publication as an Article in Scientific Reports.

Most tumor cohorts have transcriptomes but no proteomes, so mRNA is used as a stand-in for protein, for example when judging whether an antibody–drug conjugate target is expressed. Whether that is reasonable depends on how well a gene's mRNA tracks its protein across tumors, and this concordance is known only where both were measured. We trained a model that forecasts the concordance from RNA data alone, using 1,027 tumors of 10 cancer types from CPTAC, and tested it on three external platforms in preregistered analyses: reverse-phase protein arrays of TCGA tumors, DIA proteomes of cancer cell lines and a TMT lung adenocarcinoma cohort.

The forecast ranked observed concordance better than an RNA-variance baseline on all three platforms. In post hoc comparisons, it was at least as accurate as looking up observed CPTAC correlations, pooled across cancer types or taken from the matched type, for genes whose CPTAC data it had been given; for genes whose CPTAC correlations the model had not used, a lookup was better, and we report this openly. A short calibration with 50 locally measured genes made its 80% intervals hold on all three platforms. We provide forecasts for antibody–drug conjugate targets across 33 TCGA cancer types and for 19,211 genes.

We believe the work suits Scientific Reports because its value lies in a technically careful answer to a practical question that many cancer researchers face. We report every preregistered test, including one that failed, label all post hoc analyses as such, and release the code, the frozen predictions, the trained model and the forecast tables ([AUTHOR TO FILL: Zenodo DOI]).

As required by Springer Nature policy, the use of large language models in this study is documented in the Methods. Claude Opus 5.5 (Anthropic). A coding agent in the Cursor editor running Grok 4.7 (xAI). The manuscript has not been published and is not under consideration elsewhere. All authors have approved the submission and declare no competing interests.

Sincerely,

Junjeong Choi
College of Pharmacy, Yonsei Institute of Pharmaceutical Sciences, Yonsei University, Incheon 21983, Republic of Korea
junjeong@yonsei.ac.kr
