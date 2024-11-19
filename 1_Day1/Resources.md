__TABLE OF CONTENT__

- [Dataset information](#dataset-information)
- [Manuals on scRNASeq data analysis](#manuals-on-scrnaseq-data-analysis)
- [Information on tools and workflows](#information-on-tools-and-workflows)
  * [Scanpy ecosystem](#scanpy-ecosystem)
  * [Normalization](#normalization)
  * [Feature selection](#feature-selection)
  * [Batch correction](#batch-correction)
  * [Pseudobulk approaches](#pseudobulk-calculation)
  * [Metacell approaches](#metacell-calculation)



<br>

## Dataset information

* Paper Reference: [__Chromatin and gene-regulatory dynamics of the developing human cerebral cortex at single-cell resolution__](https://doi.org/10.1016/j.cell.2021.07.039)
* DOI: [https://doi.org/10.1016/j.cell.2021.07.039](https://doi.org/10.1016/j.cell.2021.07.039)
* Original Code - Resources: [Dataset Interactive Viewer](https://scbrainregulation.su.domains/)


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

### PseudoBulk calculation

* [Decoupler pseudobulk workflow](https://decoupler-py.readthedocs.io/en/latest/notebooks/pseudobulk.html#Generation-of-pseudo-bulk-profiles)
* [Decoupler general](https://decoupler-py.readthedocs.io/en/1.3.2/index.html)
* [Pertpy](https://pertpy.readthedocs.io/en/stable/tutorials/notebooks/differential_gene_expression.html)

### Metacell calculation

* [SEACells paper](https://www.nature.com/articles/s41587-023-01716-9)
* [SEACells GitHub](https://github.com/dpeerlab/SEACells)
* [SEACells workflow](https://github.com/dpeerlab/SEACells/blob/main/notebooks/SEACell_computation.ipynb)
* [Metacells2 paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02667-1)
* [Metacells2 website](https://www.weizmann.ac.il/math/tanay/research-activities/metacell-2)
* [Metacells2 documentation](https://metacells.readthedocs.io/en/latest/)
