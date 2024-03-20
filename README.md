## Identifying Reproducible Transcription Regulator Coexpression Patterns with Single Cell Transcriptomics
[https://www.biorxiv.org/content/10.1101/2024.02.15.580581v3]

Analysis relies on objects (i.e., the underlying scRNA-seq datasets and resulting
gene x gene coexpression matrices) that live on the Pavlab servers and are not easily 
shareable to a data repo given their size. Please contact if you have any 
questions or requests!

The summarized rankings are provided in Borealis [Repo will be made public upon publication]

Pavlab: The data and the associated download scripts are found:
/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/amorin_sc_datasets

Once scRNA-seq data was acquired, the scripts used to preprocess each and build
an aggregate coexpression network are found in this repo:
R/preprocessing_scripts/CPM/

The Unibind ChIP-seq data used for binding evidence were downloaded and 
summarized using:
https://github.com/PavlidisLab/Unibind_analysis/
For convenience, the Borealis repo also contains the gene x experiment binding 
matrices used to generate the aggregate binding summaries.
