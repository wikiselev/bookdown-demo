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

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{18-confounders_files/figure-latex/confound-pca-1} 

}

\caption{PCA plot of the tung data}(\#fig:confound-pca)
\end{figure}

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

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{18-confounders_files/figure-latex/confound-find-pcs-total-features-1} 

}

\caption{PC correlation with the number of detected genes}(\#fig:confound-find-pcs-total-features)
\end{figure}

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

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{18-confounders_files/figure-latex/confound-find-expl-vars-1} 

}

\caption{Explanatory variables}(\#fig:confound-find-expl-vars)
\end{figure}

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
## R version 3.4.3 (2017-11-30)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Debian GNU/Linux 9 (stretch)
## 
## Matrix products: default
## BLAS: /usr/lib/openblas-base/libblas.so.3
## LAPACK: /usr/lib/libopenblasp-r0.2.19.so
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
##  [1] scater_1.6.3               SingleCellExperiment_1.0.0
##  [3] SummarizedExperiment_1.8.1 DelayedArray_0.4.1        
##  [5] matrixStats_0.53.1         GenomicRanges_1.30.3      
##  [7] GenomeInfoDb_1.14.0        IRanges_2.12.0            
##  [9] S4Vectors_0.16.0           ggplot2_2.2.1             
## [11] Biobase_2.38.0             BiocGenerics_0.24.0       
## [13] knitr_1.20                
## 
## loaded via a namespace (and not attached):
##  [1] viridis_0.5.0          httr_1.3.1             edgeR_3.20.9          
##  [4] bit64_0.9-7            viridisLite_0.3.0      shiny_1.0.5           
##  [7] assertthat_0.2.0       blob_1.1.0             vipor_0.4.5           
## [10] GenomeInfoDbData_1.0.0 yaml_2.1.17            progress_1.1.2        
## [13] pillar_1.2.1           RSQLite_2.0            backports_1.1.2       
## [16] lattice_0.20-34        glue_1.2.0             limma_3.34.9          
## [19] digest_0.6.15          XVector_0.18.0         colorspace_1.3-2      
## [22] cowplot_0.9.2          htmltools_0.3.6        httpuv_1.3.6.1        
## [25] Matrix_1.2-7.1         plyr_1.8.4             XML_3.98-1.10         
## [28] pkgconfig_2.0.1        biomaRt_2.34.2         bookdown_0.7          
## [31] zlibbioc_1.24.0        xtable_1.8-2           scales_0.5.0          
## [34] tibble_1.4.2           lazyeval_0.2.1         magrittr_1.5          
## [37] mime_0.5               memoise_1.1.0          evaluate_0.10.1       
## [40] beeswarm_0.2.3         shinydashboard_0.6.1   tools_3.4.3           
## [43] data.table_1.10.4-3    prettyunits_1.0.2      stringr_1.3.0         
## [46] munsell_0.4.3          locfit_1.5-9.1         AnnotationDbi_1.40.0  
## [49] bindrcpp_0.2           compiler_3.4.3         rlang_0.2.0           
## [52] rhdf5_2.22.0           grid_3.4.3             RCurl_1.95-4.10       
## [55] tximport_1.6.0         rjson_0.2.15           labeling_0.3          
## [58] bitops_1.0-6           rmarkdown_1.8          gtable_0.2.0          
## [61] DBI_0.7                reshape2_1.4.3         R6_2.2.2              
## [64] gridExtra_2.3          dplyr_0.7.4            bit_1.1-12            
## [67] bindr_0.1              rprojroot_1.3-2        ggbeeswarm_0.6.0      
## [70] stringi_1.1.6          Rcpp_0.12.15           xfun_0.1
```
