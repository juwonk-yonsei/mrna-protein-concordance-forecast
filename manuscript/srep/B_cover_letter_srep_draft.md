<!-- B_cover_letter_srep_draft.md — cover letter draft for Scientific Reports, written by Claude. [FILL]/[AUTHOR TO FILL] items are for the researcher or the agent as noted in B 지시서 10. -->

[AUTHOR TO FILL: date, e.g. 20 October 2026]

The Editors
Scientific Reports

Dear Editors,

Please consider our manuscript, "A preregistered model forecasts mRNA–protein concordance in unmeasured cancer contexts across three proteomic platforms", for publication as an Article in Scientific Reports.

Most tumor cohorts have transcriptomes but no proteomes, so mRNA is used as a stand-in for protein, for example when judging whether an antibody–drug conjugate target is expressed. Whether that is reasonable depends on how well a gene's mRNA tracks its protein across tumors, and this concordance is known only where both were measured. We trained a model that forecasts the concordance from RNA data alone, using 1,027 tumors of 10 cancer types from CPTAC, and tested it on three external platforms in preregistered analyses: reverse-phase protein arrays of TCGA tumors, DIA proteomes of cancer cell lines and a TMT lung adenocarcinoma cohort.

The forecast ranked observed concordance better than an RNA-variance baseline on all three platforms. In post hoc comparisons, it was at least as accurate as looking up observed CPTAC correlations, even from the matched cancer type, and a short calibration with 50 locally measured genes made its 80% intervals hold on all three platforms. We provide forecasts for antibody–drug conjugate targets across 33 TCGA cancer types [FILL: "and for N genes" if OPTION A].

We believe the work suits Scientific Reports because its value lies in a technically careful answer to a practical question that many cancer researchers face. We report every preregistered test, including one that failed, label all post hoc analyses as such, and release the code, the frozen predictions, the trained model and the forecast tables ([FILL: Zenodo DOI]).

As required by Springer Nature policy, the use of large language models in this study is documented in the Methods. The manuscript has not been published and is not under consideration elsewhere. All authors have approved the submission and declare no competing interests. [AUTHOR TO CONFIRM]

[AUTHOR TO FILL: optional suggested reviewers with affiliations and emails; optional excluded reviewers]

Sincerely,

[AUTHOR TO FILL: name, affiliation, email]
