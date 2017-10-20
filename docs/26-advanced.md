---
output: html_document
---

# Advanced exercises

For the final part of the course we would like you to work on more open ended problems. The goal is to carry out the type of analyses that you would be doing for an actual research project. 

Participants who have their own dataset that they are interested in should feel free to work with them. 

For other participants we recommend downloading a dataset from the [conquer](http://imlspenticton.uzh.ch:3838/conquer/) resource (<b>co</b>nsistent <b>qu</b>antification of <b>e</b>xternal <b>r</b>na-seq data). [conquer](http://imlspenticton.uzh.ch:3838/conquer/) uses [Salmon](http://salmon.readthedocs.io/en/latest/salmon.html) to quantify the transcript abundances in a given sample. For a given organism, the fasta files containing cDNA and ncRNA sequences from Ensembl are complemented with ERCC spike-in sequences, and a Salmon quasi-mapping index is built for the entire catalog. Then [Salmon](http://salmon.readthedocs.io/en/latest/salmon.html) is run to estimate the abundance of each transcript. The abundances estimated by [Salmon](http://salmon.readthedocs.io/en/latest/salmon.html) are summarised and provided to the user in the form of a `MultiAssayExperiment` object. This object can be downloaded via the buttons in the _MultiAssayExperiment_ column. The provided `MultiAssayExperiment` object contains two "experiments", corresponding to the gene-level and transcript-level values.

The gene-level experiment contains four "assays":

* TPM
* count
* count_lstpm (count-scale length-scaled TPMs)
* avetxlength (the average transcript length, which can be used as offsets in count models based on the count assay, see [here](http://f1000research.com/articles/4-1521/)).

The transcript-level experiment contains three "assays":

* TPM
* count
* efflength (the effective length estimated by [Salmon](http://salmon.readthedocs.io/en/latest/salmon.html))

The `MultiAssayExperiment` also contains the phenotypic data (in the `pData` slot), as well as some metadata for the data set (the genome, the organism and the Salmon index that was used for the quantification).

Here we will show you how to create an SCESet from a `MultiAssayExperiment` object. For example, if you download `Shalek2013` dataset you will be able to create an SCESet using the following code:


```r
library(MultiAssayExperiment)
library(SummarizedExperiment)
library(scater)
d <- readRDS("~/Desktop/GSE41265.rds")
cts <- assays(experiments(d)[["gene"]])[["count_lstpm"]]
tpms <- assays(experiments(d)[["gene"]])[["TPM"]]
phn <- pData(d)
sceset <- newSCESet(
    countData = cts, 
    tpmData = tpms,
    phenoData = new("AnnotatedDataFrame", data = as.data.frame(phn))
)
```

You can also see that several different QC metrics have already been pre-calculated on the [conquer](http://imlspenticton.uzh.ch:3838/conquer/) website.

Here are some suggestions for questions that you can explore:

* There are two mESC datasets from different labs (i.e. Xue and Kumar). Can you merge them and remove the batch effects?

* Clustering and pseudotime analysis look for different patterns among cells. How might you tell which is more appropriate for your dataset?

* One of the main challenging in hard clustering is to identify the appropriate value for k. Can you use one or more of the clustering tools to explore the different hierarchies available? What are good mathematical and/or biological criteria for determining k?

* The choice of normalization strategy matters, but how do you determine which is the best method? Explore the effect of different normalizations on downstream analyses.

* scRNA-seq datasets are high-dimensional and since most dimensions (ie genes) are not informative. Consequently, dimensionality reduction and feature selection are important when analyzing and visualizing the data. Consider the effect of different feature selection methods and dimensionality reduction on clustering and/or pseudotime inference.

* One of the main challenges after clustering cells is to interpret the biological relevance of the subpopulations. One approach is to identify [gene ontology](http://geneontology.org/) terms that are enriched for the set of marker genes. Identify marker genes (e.g. using SC3 or M3Drop) and explore the ontology terms using [gProfiler](http://biit.cs.ut.ee/gprofiler/), [WebGestalt](http://www.webgestalt.org/) or [DAVID](https://david.ncifcrf.gov/).

* Similarly, when ordering cells according to pseudotime we would like to understand what underlying cellular processes are changing over time. Identify a set of changing genes from the aligned cells and use ontology terms to characterize them.
