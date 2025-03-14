## Identifying Reproducible Transcription Regulator Coexpression Patterns with Single Cell Transcriptomics
[https://www.biorxiv.org/content/10.1101/2024.02.15.580581v3]
TODO: update when PLOS CB link is live!

This project looks to identify which gene partners are most commonly coexpressed
with each TR in human and mouse, across a large corpus of single cell RNA-seq
data. This information is compared to literature curated targets from 
low-throughput experiments, as well as aggregated ChIP-seq binding scores.

A main deliverable is the summarized/ranked information for each TR, which can
be found at: <https://borealisdata.ca/dataset.xhtml?persistentId=doi:10.5683/SP3/HJ1B24>

NOTE: Analysis relies on objects (i.e., the underlying scRNA-seq datasets and 
resulting gene x gene coexpression matrices) that live on the Pavlab servers and
are not easily shareable to a data repo given their size. Please contact me if 
you have any questions or requests!


Pavlab: The data and the associated download scripts are found:

```/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/amorin_sc_datasets```

Once scRNA-seq data was acquired, the scripts used to preprocess and build
an aggregate coexpression network for each dataset are found in this repo:

```R/preprocessing_scripts/CPM/```

The Unibind ChIP-seq data used for binding evidence were downloaded and 
summarized using:
<https://github.com/PavlidisLab/Unibind_analysis/>

For convenience, the Borealis repo above also contains the gene x experiment 
binding matrices used to generate the aggregate binding summaries.
