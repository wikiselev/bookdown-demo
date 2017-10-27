---
output: html_document
---

## Data visualization

### Introduction

In this chapter we will continue to work with the filtered `Tung` dataset produced in the previous chapter. We will explore different ways of visualizing the data to allow you to asses what happened to the expression matrix after the quality control step. `scater` package provides several very useful functions to simplify visualisation. 

One important aspect of single-cell RNA-seq is to control for batch effects. Batch effects are technical artefacts that are added to the samples during handling. For example, if two sets of samples were prepared in different labs or even on different days in the same lab, then we may observe greater similarities between the samples that were handled together. In the worst case scenario, batch effects may be [mistaken](http://f1000research.com/articles/4-121/v1) for true biological variation. The `Tung` data allows us to explore these issues in a controlled manner since some of the salient aspects of how the samples were handled have been recorded. Ideally, we expect to see batches from the same individual grouping together and distinct groups corresponding to each individual. 




```r
library(SingleCellExperiment)
library(scater)
options(stringsAsFactors = FALSE)
umi <- readRDS("tung/umi.rds")
umi.qc <- umi[rowData(umi)$use, colData(umi)$use]
endog_genes <- !rowData(umi.qc)$is_feature_control
```

### PCA plot {#visual-pca}

The easiest way to overview the data is by transforming it using the principal component analysis and then visualize the first two principal components.

[Principal component analysis (PCA)](https://en.wikipedia.org/wiki/Principal_component_analysis) is a statistical procedure that uses a transformation to convert a set of observations into a set of values of linearly uncorrelated variables called principal components (PCs). The number of principal components is less than or equal to the number of original variables.

Mathematically, the PCs correspond to the eigenvectors of the covariance matrix. The eigenvectors are sorted by eigenvalue so that the first principal component accounts for as much of the variability in the data as possible, and each succeeding component in turn has the highest variance possible under the constraint that it is orthogonal to the preceding components (the figure below is taken from [here](http://www.nlpca.org/pca_principal_component_analysis.html)).

<div class="figure" style="text-align: center">
<img src="figures/pca.png" alt="Schematic representation of PCA dimensionality reduction" width="100%" />
<p class="caption">(\#fig:clust-pca)Schematic representation of PCA dimensionality reduction</p>
</div>

#### Before QC

Without log-transformation:

```r
plotPCA(
    umi[endog_genes, ],
    exprs_values = "counts",
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual"
)
```

<div class="figure" style="text-align: center">
<img src="09-exprs-overview_files/figure-html/expr-overview-pca-before-qc1-1.png" alt="PCA plot of the tung data" width="90%" />
<p class="caption">(\#fig:expr-overview-pca-before-qc1)PCA plot of the tung data</p>
</div>

With log-transformation:

```r
plotPCA(
    umi[endog_genes, ],
    exprs_values = "logcounts_raw",
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual"
)
```

<div class="figure" style="text-align: center">
<img src="09-exprs-overview_files/figure-html/expr-overview-pca-before-qc2-1.png" alt="PCA plot of the tung data" width="90%" />
<p class="caption">(\#fig:expr-overview-pca-before-qc2)PCA plot of the tung data</p>
</div>

Clearly log-transformation is benefitial for our data - it reduces the variance on the first principal component and already separates some biological effects. Moreover, it makes the distribution of the expression values more normal. In the following analysis and chapters we will be using log-transformed raw counts by default.

__However, note that just a log-transformation is not enough to account for different technical factors between the cells (e.g. sequencing depth). Therefore, please do not use `logcounts_raw` for your downstream analysis, instead as a minimum suitable data use the `logcounts` slot of the `SingleCellExperiment` object, which not just log-transformed, but also normalised by library size (e.g. CPM normalisation). In the course we use `logcounts_raw` only for demonstration purposes!__

#### After QC


```r
plotPCA(
    umi.qc[endog_genes, ],
    exprs_values = "logcounts_raw",
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual"
)
```

<div class="figure" style="text-align: center">
<img src="09-exprs-overview_files/figure-html/expr-overview-pca-after-qc-1.png" alt="PCA plot of the tung data" width="90%" />
<p class="caption">(\#fig:expr-overview-pca-after-qc)PCA plot of the tung data</p>
</div>

Comparing Figure \@ref(fig:expr-overview-pca-before-qc2) and Figure \@ref(fig:expr-overview-pca-after-qc), it is clear that after quality control the NA19098.r2 cells no longer form a group of outliers.

By default only the top 500 most variable genes are used by scater to calculate the PCA. This can be adjusted by changing the `ntop` argument. 

__Exercise 1__
How do the PCA plots change if when all 14,214 genes are used? Or when only top 50 genes are used? Why does the fraction of variance accounted for by the first PC change so dramatically?

__Hint__ Use `ntop` argument of the `plotPCA` function.

__Our answer__

<div class="figure" style="text-align: center">
<img src="09-exprs-overview_files/figure-html/expr-overview-pca-after-qc-exercise1-1-1.png" alt="PCA plot of the tung data (14214 genes)" width="90%" />
<p class="caption">(\#fig:expr-overview-pca-after-qc-exercise1-1)PCA plot of the tung data (14214 genes)</p>
</div>

<div class="figure" style="text-align: center">
<img src="09-exprs-overview_files/figure-html/expr-overview-pca-after-qc-exercise1-2-1.png" alt="PCA plot of the tung data (50 genes)" width="90%" />
<p class="caption">(\#fig:expr-overview-pca-after-qc-exercise1-2)PCA plot of the tung data (50 genes)</p>
</div>

If your answers are different please compare your code with [ours](https://github.com/hemberg-lab/scRNA.seq.course/blob/master/07-exprs-overview.Rmd) (you need to search for this exercise in the opened file).

### tSNE map {#visual-tsne}

An alternative to PCA for visualizing scRNASeq data is a tSNE plot. [tSNE](https://lvdmaaten.github.io/tsne/) (t-Distributed Stochastic Neighbor Embedding) combines dimensionality reduction (e.g. PCA) with random walks on the nearest-neighbour network to map high dimensional data (i.e. our 14,214 dimensional expression matrix) to a 2-dimensional space while preserving local distances between cells. In contrast with PCA, tSNE is a stochastic algorithm which means running the method multiple times on the same dataset will result in different plots. Due to the non-linear and stochastic nature of the algorithm, tSNE is more difficult to intuitively interpret tSNE. To ensure reproducibility, we fix the "seed" of the random-number generator in the code below so that we always get the same plot. 


#### Before QC


```r
plotTSNE(
    umi[endog_genes, ],
    exprs_values = "logcounts_raw",
    perplexity = 130,
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    rand_seed = 123456
)
```

<div class="figure" style="text-align: center">
<img src="09-exprs-overview_files/figure-html/expr-overview-tsne-before-qc-1.png" alt="tSNE map of the tung data" width="90%" />
<p class="caption">(\#fig:expr-overview-tsne-before-qc)tSNE map of the tung data</p>
</div>

#### After QC


```r
plotTSNE(
    umi.qc[endog_genes, ],
    exprs_values = "logcounts_raw",
    perplexity = 130,
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    rand_seed = 123456
)
```

<div class="figure" style="text-align: center">
<img src="09-exprs-overview_files/figure-html/expr-overview-tsne-after-qc-1.png" alt="tSNE map of the tung data" width="90%" />
<p class="caption">(\#fig:expr-overview-tsne-after-qc)tSNE map of the tung data</p>
</div>

Interpreting PCA and tSNE plots is often challenging and due to their stochastic and non-linear nature, they are less intuitive. However, in this case it is clear that they provide a similar picture of the data. Comparing Figure \@ref(fig:expr-overview-tsne-before-qc) and \@ref(fig:expr-overview-tsne-after-qc), it is again clear that the samples from NA19098.r2 are no longer outliers after the QC filtering.

Furthermore tSNE requires you to provide a value of `perplexity` which reflects the number of neighbours used to build the nearest-neighbour network; a high value creates a dense network which clumps cells together while a low value makes the network more sparse allowing groups of cells to separate from each other. `scater` uses a default perplexity of the total number of cells divided by five (rounded down).

You can read more about the pitfalls of using tSNE [here](http://distill.pub/2016/misread-tsne/).

__Exercise 2__
How do the tSNE plots change when a perplexity of 10 or 200 is used? How does the choice of perplexity affect the interpretation of the results?

__Our answer__

<div class="figure" style="text-align: center">
<img src="09-exprs-overview_files/figure-html/expr-overview-tsne-after-qc-exercise2-1-1.png" alt="tSNE map of the tung data (perplexity = 10)" width="90%" />
<p class="caption">(\#fig:expr-overview-tsne-after-qc-exercise2-1)tSNE map of the tung data (perplexity = 10)</p>
</div>

<div class="figure" style="text-align: center">
<img src="09-exprs-overview_files/figure-html/expr-overview-tsne-after-qc-exercise2-2-1.png" alt="tSNE map of the tung data (perplexity = 200)" width="90%" />
<p class="caption">(\#fig:expr-overview-tsne-after-qc-exercise2-2)tSNE map of the tung data (perplexity = 200)</p>
</div>

### Big Exercise

Perform the same analysis with read counts of the Blischak data. Use `tung/reads.rds` file to load the reads SCESet object. Once you have finished please compare your results to ours (next chapter).

### sessionInfo()


```
## R version 3.4.2 (2017-09-28)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Debian GNU/Linux 9 (stretch)
## 
## Matrix products: default
## BLAS/LAPACK: /usr/lib/libopenblasp-r0.2.19.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C             
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats4    methods   stats     graphics  grDevices utils    
## [8] datasets  base     
## 
## other attached packages:
##  [1] scater_1.5.21               ggplot2_2.2.1              
##  [3] SingleCellExperiment_0.99.4 SummarizedExperiment_1.6.5 
##  [5] DelayedArray_0.2.7          matrixStats_0.52.2         
##  [7] Biobase_2.36.2              GenomicRanges_1.28.6       
##  [9] GenomeInfoDb_1.12.3         IRanges_2.10.5             
## [11] S4Vectors_0.14.7            BiocGenerics_0.22.1        
## [13] knitr_1.17                 
## 
## loaded via a namespace (and not attached):
##  [1] viridis_0.4.0           edgeR_3.18.1           
##  [3] bit64_0.9-7             viridisLite_0.2.0      
##  [5] shiny_1.0.5             assertthat_0.2.0       
##  [7] highr_0.6               blob_1.1.0             
##  [9] GenomeInfoDbData_0.99.0 vipor_0.4.5            
## [11] yaml_2.1.14             RSQLite_2.0            
## [13] backports_1.1.1         lattice_0.20-34        
## [15] glue_1.1.1              limma_3.32.10          
## [17] digest_0.6.12           XVector_0.16.0         
## [19] colorspace_1.3-2        cowplot_0.8.0          
## [21] htmltools_0.3.6         httpuv_1.3.5           
## [23] Matrix_1.2-7.1          plyr_1.8.4             
## [25] XML_3.98-1.9            pkgconfig_2.0.1        
## [27] biomaRt_2.32.1          bookdown_0.5           
## [29] zlibbioc_1.22.0         xtable_1.8-2           
## [31] scales_0.5.0            Rtsne_0.13             
## [33] tibble_1.3.4            lazyeval_0.2.0         
## [35] magrittr_1.5            mime_0.5               
## [37] memoise_1.1.0           evaluate_0.10.1        
## [39] beeswarm_0.2.3          shinydashboard_0.6.1   
## [41] tools_3.4.2             data.table_1.10.4-3    
## [43] stringr_1.2.0           munsell_0.4.3          
## [45] locfit_1.5-9.1          AnnotationDbi_1.38.2   
## [47] bindrcpp_0.2            compiler_3.4.2         
## [49] rlang_0.1.2             rhdf5_2.20.0           
## [51] grid_3.4.2              RCurl_1.95-4.8         
## [53] tximport_1.4.0          rjson_0.2.15           
## [55] labeling_0.3            bitops_1.0-6           
## [57] rmarkdown_1.6           gtable_0.2.0           
## [59] DBI_0.7                 reshape2_1.4.2         
## [61] R6_2.2.2                gridExtra_2.3          
## [63] dplyr_0.7.4             bit_1.1-12             
## [65] bindr_0.1               rprojroot_1.2          
## [67] stringi_1.1.5           ggbeeswarm_0.6.0       
## [69] Rcpp_0.12.13
```
