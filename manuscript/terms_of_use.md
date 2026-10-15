# Terms received on 2026-10-15

Addresses are those in `logs/download_log.tsv`, the manuscript Data availability section, or `release/terms_to_check.md`. Sentences below are copied from the response body. Where the body had no license or citation sentence, that is stated and no sentence is supplied.

No received sentence forbids redistribution of the derived result tables or the model files.

## 1. CPTAC / Proteomic Data Commons

Requested URL: https://kb.linkedomics.org/
Date: 2026-10-15
HTTP status: 200
Final URL: https://kb.linkedomics.org/
Bytes stored: 20758

Verbatim citation sentence on the page:

> Please cite: A proteogenomics data-driven knowledge base of human cancer, Yuxing Liao, Sara R. Savage, Yongchao Dou, Zhiao Shi, Xinpei Yi, Wen Jiang, Jonathan T. Lei, Bing Zhang, Cell Systems, 2023.

The page links that title to https://doi.org/10.1016/j.cels.2023.07.007.

Footer on the same page:

> © 2022-2023 Zhang Lab

The download page https://kb.linkedomics.org/download returned HTTP 200 (70031 bytes). Its visible text includes "CPTAC Pan-cancer Data Download". No license sentence was in that body.

Requested URL: https://pdc.cancer.gov/pdc/
Date: 2026-10-15
HTTP status: 200
Final URL: https://pdc.cancer.gov/pdc/
Bytes stored: 56105

The body is an application shell. No license, citation, or terms sentence was in the HTML.

Manuscript contrast: reference 10 is Liao et al., Cell Systems 2023. The acknowledgements name the Clinical Proteomic Tumor Analysis Consortium (NCI/NIH). No further citation sentence was received. The manuscript was not edited.

## 2. The Cancer Proteome Atlas

Requested URL: https://tcpaportal.org/tcpa/
Date: 2026-10-15
HTTP status: 200
Final URL: https://tcpa.drbioright.org/rppa500/main.html
Bytes stored: 9324

Meta element on the response:

> Copyright content="2013, MDACC"

Body sentence:

> Javascript required to use this web app.

No redistribution sentence was in the body.

Manuscript contrast: reference 15 is Li et al., Nature Methods 2013. No further citation or redistribution sentence was received. The manuscript was not edited.

## 3. PRIDE PXD030304

Requested URL: https://www.ebi.ac.uk/pride/archive/projects/PXD030304
Date: 2026-10-15
HTTP status: 200
Final URL: https://www.ebi.ac.uk/pride/archive/projects/PXD030304
Bytes stored: 7222

The body is a PRIDE application shell. No license sentence was in the HTML.

Manuscript contrast: the ProCan paragraph names PRIDE PXD030304 and cites Gonçalves et al. as reference 16. No license sentence was received to compare. The manuscript was not edited.

## 4. Gonçalves et al. 2022 data availability

Requested URL: https://doi.org/10.1016/j.ccell.2022.06.010
Date: 2026-10-15
HTTP status: 200
Final URL: Elsevier linkinghub retrieve page for pii S1535610822002744
Bytes stored: 2661

The body title is "Redirecting". No data-availability sentence was in the response.

The redirect target written in that page is https://www.sciencedirect.com/science/article/pii/S1535610822002744. A follow-up request returned HTTP 403, curl exit 56, 0 bytes. Data-availability sentence: 열지 못함.

A page already stored from the download log (https://depmap.sanger.ac.uk/documentation/datasets/proteomics/, file `data/raw/docs/sanger_proteomics.html`, downloaded 2026-10-01, HTTP status in the log: ok) contains:

> Publication reference: Gonçalves, E., Poulos, R.C., Cai, Z., et al., Pan-cancer proteomic map of 949 human cell lines, Cancer Cell, (2022)

and names raw .wiff files at https://www.ebi.ac.uk/pride/archive/projects/PXD030304. That stored page has no license sentence. It is not the paper's data-availability section.

Manuscript contrast: reference 16 is Gonçalves et al., Cancer Cell 2022. The paper's data-availability sentence was not received, so no missing-sentence comparison was made. The manuscript was not edited.

## 5. DepMap 24Q4

Requested URL: https://depmap.org/portal/data_page/?tab=allData
Date: 2026-10-15
HTTP status: 200
Bytes stored: 4846

The body title is "DepMap — Verification". The page is a Cloudflare Turnstile challenge. No license or citation sentence was in the body.

Requested URL: https://doi.org/10.25452/figshare.plus.27993248.v1
Date: 2026-10-15
HTTP status: 202
Final URL: https://plus.figshare.com/articles/dataset/DepMap_24Q4_Public/27993248/1
Bytes stored: 0

A direct request of that figshare URL also returned HTTP 202 and 0 bytes. License sentence: 열지 못함.

The file already stored from the download log (`data/raw/depmap/README_24Q4.txt`, URL https://ndownloader.figshare.com/files/51065795, downloaded 2026-10-01) has no line containing license, CC BY, citation, terms, or copyright.

Manuscript contrast: reference 29 is Broad DepMap, DepMap 24Q4 Public, with that DOI. No citation requirement was received. The manuscript was not edited.

## 6. Gygi CCLE proteome

Requested URL: https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz
Date: 2026-10-15
HTTP status: 200
Content type: application/x-gzip
Bytes stored: 400000 (first part of the data file)

The response is the data file, not a terms page. No use sentence was in the response.

The paper DOI named in the manuscript, https://doi.org/10.1016/j.cell.2019.12.023, returned HTTP 200 and a page titled "Redirecting" to pii S0092867419313856. The ScienceDirect follow-up was not opened (the Gonçalves follow-up to the same host returned HTTP 403). Data-availability sentence: 열지 못함.

Manuscript contrast: reference 31 is Nusinow et al., Cell 2020. No use sentence was received to compare. The manuscript was not edited.
