---
knit: bookdown::preview_chapter
---

## Identifying confounding factors

### Introduction

There is a large number of potential confounders, artifacts and biases in sc-RNA-seq data. One of the main challenges in analyzing scRNA-seq data stems from the fact that it is difficult to carry out a true technical replicate (why?) to distinguish biological and technical variability. In the previous chapters we considered batch effects and in this chapter we will continue to explore how experimental artifacts can be identified and removed. We will continue using the `scater` package since it provides a set of methods specifically for quality control of experimental and explanatory variables. Moreover, we will continue to work with the Blischak data that was used in the previous chapter.




```r
library(scater, quietly = TRUE)
options(stringsAsFactors = FALSE)
umi <- readRDS("tung/umi.rds")
umi.qc <- umi[rowData(umi)$use, colData(umi)$use]
endog_genes <- !rowData(umi.qc)$is_feature_control
```

The `umi.qc` dataset contains filtered cells and genes. Our next step is to explore technical drivers of variability in the data to inform data normalisation before downstream analysis.

### Correlations with PCs

Let's first look again at the PCA plot of the QCed dataset:

```r
plotPCA(
    umi.qc[endog_genes, ],
    exprs_values = "logcounts_raw",
    colour_by = "batch",
    size_by = "total_features"
)
```

<div class="figure" style="text-align: center">
<img src="11-confounders_files/figure-html/confound-pca-1.png" alt="PCA plot of the tung data" width="90%" />
<p class="caption">(\#fig:confound-pca)PCA plot of the tung data</p>
</div>

`scater` allows one to identify principal components that correlate with experimental and QC variables of interest (it ranks principle components by $R^2$ from a linear model regressing PC value against the variable of interest).

Let's test whether some of the variables correlate with any of the PCs.

#### Detected genes


```r
plotQC(
    umi.qc[endog_genes, ],
    type = "find-pcs",
    exprs_values = "logcounts_raw",
    variable = "total_features"
)
```

<div class="figure" style="text-align: center">
<img src="11-confounders_files/figure-html/confound-find-pcs-total-features-1.png" alt="PC correlation with the number of detected genes" width="90%" />
<p class="caption">(\#fig:confound-find-pcs-total-features)PC correlation with the number of detected genes</p>
</div>

Indeed, we can see that `PC1` can be almost completely explained by the number of detected genes. In fact, it was also visible on the PCA plot above. This is a well-known issue in scRNA-seq and was described [here](http://biorxiv.org/content/early/2015/12/27/025528).

### Explanatory variables

`scater` can also compute the marginal $R^2$ for each variable when fitting a linear model regressing expression values for each gene against just that variable, and display a density plot of the gene-wise marginal $R^2$ values for the variables.


```r
plotQC(
    umi.qc[endog_genes, ],
    type = "expl",
    exprs_values = "logcounts_raw",
    variables = c(
        "total_features",
        "total_counts",
        "batch",
        "individual",
        "pct_counts_ERCC",
        "pct_counts_MT"
    )
)
```

<div class="figure" style="text-align: center">
<img src="11-confounders_files/figure-html/confound-find-expl-vars-1.png" alt="Explanatory variables" width="90%" />
<p class="caption">(\#fig:confound-find-expl-vars)Explanatory variables</p>
</div>

This analysis indicates that the number of detected genes (again) and also the sequencing depth (number of counts) have substantial explanatory power for many genes, so these variables are good candidates for conditioning out in a normalisation step, or including in downstream statistical models. Expression of ERCCs also appears to be an important explanatory variable and one notable feature of the above plot is that batch explains more than individual. What does that tell us about the technical and biological variability of the data?

### Other confounders

In addition to correcting for batch, there are other factors that one
may want to compensate for. As with batch correction, these
adjustments require extrinsic information. One popular method is
[scLVM](https://github.com/PMBio/scLVM) which allows you to identify
and subtract the effect from processes such as cell-cycle or
apoptosis.

In addition, protocols may differ in terms of their coverage of each transcript, 
their bias based on the average content of __A/T__ nucleotides, or their ability to capture short transcripts.
Ideally, we would like to compensate for all of these differences and biases.

### Exercise

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
## [1] stats4    parallel  methods   stats     graphics  grDevices utils    
## [8] datasets  base     
## 
## other attached packages:
##  [1] scater_1.5.21               SingleCellExperiment_0.99.4
##  [3] SummarizedExperiment_1.6.5  DelayedArray_0.2.7         
##  [5] matrixStats_0.52.2          GenomicRanges_1.28.6       
##  [7] GenomeInfoDb_1.12.3         IRanges_2.10.5             
##  [9] S4Vectors_0.14.7            ggplot2_2.2.1              
## [11] Biobase_2.36.2              BiocGenerics_0.22.1        
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
## [31] scales_0.5.0            tibble_1.3.4           
## [33] lazyeval_0.2.0          magrittr_1.5           
## [35] mime_0.5                memoise_1.1.0          
## [37] evaluate_0.10.1         beeswarm_0.2.3         
## [39] shinydashboard_0.6.1    tools_3.4.2            
## [41] data.table_1.10.4-3     stringr_1.2.0          
## [43] munsell_0.4.3           locfit_1.5-9.1         
## [45] AnnotationDbi_1.38.2    bindrcpp_0.2           
## [47] compiler_3.4.2          rlang_0.1.2            
## [49] rhdf5_2.20.0            grid_3.4.2             
## [51] RCurl_1.95-4.8          tximport_1.4.0         
## [53] rjson_0.2.15            labeling_0.3           
## [55] bitops_1.0-6            rmarkdown_1.6          
## [57] gtable_0.2.0            DBI_0.7                
## [59] reshape2_1.4.2          R6_2.2.2               
## [61] gridExtra_2.3           dplyr_0.7.4            
## [63] bit_1.1-12              bindr_0.1              
## [65] rprojroot_1.2           stringi_1.1.5          
## [67] ggbeeswarm_0.6.0        Rcpp_0.12.13
```
