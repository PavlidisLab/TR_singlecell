This directory contains scripts for pre-processing scRNA-seq data and generating
the aggregate coexpression matrices used for analysis. The CPM directory contains
scripts for counts per million normalization, while lognorm contains scripts
for Seurat's log normalization.

Data was obtained from https://cellxgene.cziscience.com/ or from GEO scrapes.

Data from cxg is further divided by whether the corresponding Seurat object was 
from a single download, or multiple downloads requiring merging. This information
is tracked in the master metadata sheet. The script expects the ID of the 
experiment (which is matched to the download directory) and the species. A bash 
script was used to batch run multiple experiments, which expects an input .tsv
containing the IDs and species to run:

```
Rscript R/preprocessing_scripts/CPM/preprocess_cellxgene_cpm.R Posner2022 Mouse
Rscript R/preprocessing_scripts/CPM/preprocess_multi_cellxgene_cpm.R GSE168215 Human
  
bash R/preprocessing_scripts/CPM/batch_preprocess_cellxgene_cpm.sh input.tsv
bash R/preprocessing_scripts/CPM/batch_preprocess_multi_cellxgene_cpm.sh multi_input.tsv
```

Data from GEO each had an individual processing script due to heterogeneity in
the downloaded file input (eg, Seurat RDS objects versus text files versus
sparse matrices) and metadata structure. A batch script was used to batch run
multiple experiments, and also expects an input .tsv containing the IDs to run:

```
Rscript R/preprocessing_scripts/CPM/GSE102827_cpm.R
bash R/preprocessing_scripts/CPM/batch_preprocess_geo_cpm.sh input.tsv
```

There are 6 datasets for which the acquired data from GEO already was processed
with no available raw count data: 

GSE129788
GSE132364
GSE195445
GSE211963
GSE231924
GSE195445Human
GSE195445Mouse


Originally Seurat's log normalization was used on all datasets with raw count 
matrices, but we switched over to CPM for a few reasons: consistency with 
prior works like Harris 2021, Eric's work on showing that log emphasizes 
correlations of lowly expressed genes, the CPM data had (very slightly) better
performance on the curation benchmark, and, most importantly for analysis here, 
the log norm data had a larger magnitude of Topk overlap in the null comparisons
compared to CPM.
