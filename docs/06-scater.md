---
output: html_document
---

## Bioconductor and `scater`

We strongly recommend all new comers and even experienced high-throughput data analysts to use well developed and maintained [Bioconductor methods and classes](https://www.bioconductor.org/developers/how-to/commonMethodsAndClasses/).

### `SingleCellExperiment` class

(this chapter is adapted from the information on the [official class page](http://bioconductor.org/packages/SingleCellExperiment))

`SingleCellExperiment` is a S4 class for storing data from single-cell experiments. This includes specialized methods to store and retrieve spike-in information, dimensionality reduction coordinates and size factors for each cell, along with the usual metadata for genes and libraries.

In practice, an object of this class can be created using the its constructor:

```r
library(SingleCellExperiment)
counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
rownames(counts) <- paste("gene", 1:10, sep = "")
colnames(counts) <- paste("cell", 1:10, sep = "")
sce <- SingleCellExperiment(
    assays = list(counts = counts),
    rowData = data.frame(gene_names = paste("gene_name", 1:10, sep = "")),
    colData = data.frame(cell_names = paste("cell_name", 1:10, sep = ""))
)
sce
```

```
## class: SingleCellExperiment 
## dim: 10 10 
## metadata(0):
## assays(1): counts
## rownames(10): gene1 gene2 ... gene9 gene10
## rowData names(1): gene_names
## colnames(10): cell1 cell2 ... cell9 cell10
## colData names(1): cell_names
## reducedDimNames(0):
## spikeNames(0):
```

In the `SingleCellExperiment`, users can assign arbitrary names to entries of assays. To assist interoperability between packages, we provide some suggestions for what the names should be for particular types of data:

* __counts__: Raw count data, e.g., number of reads or transcripts for a particular gene.
* __normcounts__: Normalized values on the same scale as the original counts. For example, counts divided by cell-specific size factors that are centred at unity.
* __logcounts__: Log-transformed counts or count-like values. In most cases, this will be defined as log-transformed normcounts, e.g., using log base 2 and a pseudo-count of 1.
* __cpm__: Counts-per-million. This is the read count for each gene in each cell, divided by the library size of each cell in millions.
* __tpm__: Transcripts-per-million. This is the number of transcripts for each gene in each cell, divided by the total number of transcripts in that cell (in millions).

Each of these suggested names has an appropriate getter/setter method for convenient manipulation of the `SingleCellExperiment`. For example, we can take the (very specifically named) `counts` slot, normalise it and assign it to `normcounts` instead:


```r
normcounts(sce) <- log2(counts(sce) + 1)
sce
```

```
## class: SingleCellExperiment 
## dim: 10 10 
## metadata(0):
## assays(2): counts normcounts
## rownames(10): gene1 gene2 ... gene9 gene10
## rowData names(1): gene_names
## colnames(10): cell1 cell2 ... cell9 cell10
## colData names(1): cell_names
## reducedDimNames(0):
## spikeNames(0):
```

```r
dim(normcounts(sce))
```

```
## [1] 10 10
```

```r
head(normcounts(sce))
```

```
##          cell1    cell2    cell3    cell4    cell5    cell6    cell7
## gene1 3.584963 3.584963 3.584963 3.584963 3.459432 3.459432 3.169925
## gene2 3.459432 3.321928 3.459432 4.000000 3.459432 3.700440 3.169925
## gene3 3.584963 3.807355 3.584963 3.906891 3.000000 2.807355 3.584963
## gene4 3.169925 3.906891 3.584963 3.584963 3.584963 3.700440 3.584963
## gene5 2.584963 2.584963 3.700440 3.807355 3.459432 3.584963 3.459432
## gene6 3.321928 3.584963 2.807355 3.169925 3.321928 3.584963 3.906891
##          cell8    cell9   cell10
## gene1 2.584963 3.584963 3.459432
## gene2 3.169925 3.169925 3.807355
## gene3 3.321928 3.000000 3.459432
## gene4 3.807355 3.700440 3.169925
## gene5 3.321928 3.700440 3.321928
## gene6 3.321928 2.807355 3.169925
```

### `scater` package

[`scater`](http://bioconductor.org/packages/scater/) is a R package for single-cell RNA-seq analysis [@McCarthy2017-kb]. The package contains several useful methods for quality control, visualisation and pre-processing of data prior to further downstream analysis.

`scater` features the following functionality:

* Automated computation of QC metrics
* Transcript quantification from read data with pseudo-alignment
* Data format standardisation
* Rich visualizations for exploratory analysis
* Seamless integration into the Bioconductor universe
* Simple normalisation methods

We highly recommend to use `scater` for all single-cell RNA-seq analyses and `scater` is the basis of the first part of the course.

As illustrated in the figure below, `scater` will help you with quality control, filtering and normalization of your expression matrix following mapping and alignment. <span style="color:red">Keep in mind that this figure represents the original version of `scater` where an `SCESet` class was used. In the newest version this figure is still correct, except that `SCESet` can be substituted with the `SingleCellExperiment` class.</span>


![](figures/scater_qc_workflow.png)
