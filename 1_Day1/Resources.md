__TABLE OF CONTENT__

- [Dataset information](#dataset-information)
- [Analysis container](#analysis-container)
- [Manuals on scRNASeq data analysis](#manuals-on-scrnaseq-data-analysis)
- [Information on tools and workflows](#information-on-tools-and-workflows)
  * [Scanpy ecosystem](#scanpy-ecosystem)
  * [Normalization](#normalization)
  * [Feature selection](#feature-selection)
  * [Batch correction](#batch-correction)
  * [Pseudobulk approaches](#pseudobulk-calculation)
  * [Metacell approaches](#metacell-calculation)
- [Useful Literature](#useful-literature)


<br>

## Dataset information

* Paper Reference: [__Chromatin and gene-regulatory dynamics of the developing human cerebral cortex at single-cell resolution__](https://doi.org/10.1016/j.cell.2021.07.039)
* DOI: [https://doi.org/10.1016/j.cell.2021.07.039](https://doi.org/10.1016/j.cell.2021.07.039)
* Original Code - Resources: [Dataset Interactive Viewer](https://scbrainregulation.su.domains/)
* You can download the starting files from [GEO entry GSE162170](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162170) selecting from the supplementary files section _GSE162170_rna_counts.tsv_ and _GSE162170_rna_cell_metadata.txt_ 
* To generate the h5ad used as input for the [first notebook](1_FiltNormBatch.ipynb), you can follow the workflow illustrated [here](Compiled/0_AssembleAdata.html)

<br>

## Analysis container

* To retrieve the singularity image used for the analyses of Day 1, you can use this singularity command that will pull the image from dockerhub : `singularity pull docker://cerebrassu/downstream:cerebra-24.06.1`


<br>

## Manuals on scRNASeq data analysis

1. [Orchestrating Single-Cell Analysis with Bioconductor](https://bioconductor.org/books/release/OSCA/)

2. [Single-cell best practices](https://www.sc-best-practices.org/preamble.html)

<br>


## Information on tools and workflows

### Scanpy ecosystem

* [Scanpy documentation](https://scanpy.readthedocs.io/en/stable/)
* [Scanpy tutorials](https://scanpy.readthedocs.io/en/stable/tutorials.html)
* [Scanpy repo](https://github.com/scverse/scanpy)

### Normalization

* [Scran paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0947-7)
* [Scran R vignette](https://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html)
* [Scran tutorial in Scanpy](https://github.com/theislab/single-cell-tutorial/blob/master/latest_notebook/Case-study_Mouse-intestinal-epithelium_1906.ipynb)
* [Benchmarking of normalization methods](https://www.nature.com/articles/s41592-023-01814-1)

### Feature selection

* [Triku paper](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac017/6547682)
* [Triku documentation](https://triku.readthedocs.io/en/latest/triku-work.html)

### Batch correction

* [Harmony Quick Start](https://portals.broadinstitute.org/harmony/articles/quickstart.html)
* [Harmony implementation in scanpy](https://scanpy.readthedocs.io/en/stable/generated/scanpy.external.pp.harmony_integrate.html)
* [Benchmarking of batch correction methods](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1850-9)
* [Benchmarking atlas-level data integration](https://www.nature.com/articles/s41592-021-01336-8)

### PseudoBulk calculation

* [Decoupler pseudobulk workflow](https://decoupler-py.readthedocs.io/en/latest/notebooks/pseudobulk.html#Generation-of-pseudo-bulk-profiles)
* [Decoupler general](https://decoupler-py.readthedocs.io/en/1.3.2/index.html)
* [Pertpy](https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/differential_gene_expression.html)
* [Confronting false discoveries in single-cell differential expression](https://www.nature.com/articles/s41467-021-25960-2)
* [Benchmarking study](https://academic.oup.com/bib/article/23/5/bbac286/6649780)

### Metacell calculation

* [SEACells paper](https://www.nature.com/articles/s41587-023-01716-9)
* [SEACells GitHub](https://github.com/dpeerlab/SEACells)
* [SEACells workflow](https://github.com/dpeerlab/SEACells/blob/main/notebooks/SEACell_computation.ipynb)
* [Metacells2 paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02667-1)
* [Metacells2 website](https://www.weizmann.ac.il/math/tanay/research-activities/metacell-2)
* [Metacells2 documentation](https://metacells.readthedocs.io/en/latest/)
* [Review on metacell approaches, pros and cons](https://www.embopress.org/doi/full/10.1038/s44320-024-00045-6)


<br>

## Useful literature

### Dimensionality reduction and UMAP debate 
* [Basic guidelines for single cell analysis](https://www.embopress.org/doi/full/10.15252/msb.20188746) here you can find the first indications and guidelines (2019 is already old in the single cell field but the principles still holds) for single cell data visualization
* [UMAP documentation](https://umap-learn.readthedocs.io/en/latest/index.html) UMAP documentation is a very rich resource to properly check how it works and how it should be used
* [UMAP for single cell visualization](https://www.nature.com/articles/nbt.4314) Original article introducing UMAP to visualise single cell data
* UMAP criticisms [link1](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011288) [link2](https://www.sciencedirect.com/science/article/pii/S2405471223002090)
* [Recent benchmarck](https://openproblems.bio/results/dimensionality_reduction/) The last is the link to the specific recent benchmarck on dimensionality reduction, but this is part of a great resource for the field, across many different tasks for single cell analysis [Preprint](https://doi.org/10.21203/rs.3.rs-4181617/v1) [Repository](https://github.com/openproblems-bio/openproblems?tab=readme-ov-file)
has context menu
* [Differential abundance analysis with Milo](https://www.nature.com/articles/s41587-021-01033-z/figures/1)
