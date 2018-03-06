---
output: html_document
---

## Expression QC (Reads)




```r
library(SingleCellExperiment)
library(scater)
options(stringsAsFactors = FALSE)
```


```r
reads <- read.table("tung/reads.txt", sep = "\t")
anno <- read.table("tung/annotation.txt", sep = "\t", header = TRUE)
```


```r
head(reads[ , 1:3])
```

```
##                 NA19098.r1.A01 NA19098.r1.A02 NA19098.r1.A03
## ENSG00000237683              0              0              0
## ENSG00000187634              0              0              0
## ENSG00000188976             57            140              1
## ENSG00000187961              0              0              0
## ENSG00000187583              0              0              0
## ENSG00000187642              0              0              0
```

```r
head(anno)
```

```
##   individual replicate well      batch      sample_id
## 1    NA19098        r1  A01 NA19098.r1 NA19098.r1.A01
## 2    NA19098        r1  A02 NA19098.r1 NA19098.r1.A02
## 3    NA19098        r1  A03 NA19098.r1 NA19098.r1.A03
## 4    NA19098        r1  A04 NA19098.r1 NA19098.r1.A04
## 5    NA19098        r1  A05 NA19098.r1 NA19098.r1.A05
## 6    NA19098        r1  A06 NA19098.r1 NA19098.r1.A06
```


```r
reads <- SingleCellExperiment(
    assays = list(counts = as.matrix(reads)), 
    colData = anno
)
```


```r
keep_feature <- rowSums(counts(reads) > 0) > 0
reads <- reads[keep_feature, ]
```


```r
isSpike(reads, "ERCC") <- grepl("^ERCC-", rownames(reads))
isSpike(reads, "MT") <- rownames(reads) %in% 
    c("ENSG00000198899", "ENSG00000198727", "ENSG00000198888",
    "ENSG00000198886", "ENSG00000212907", "ENSG00000198786",
    "ENSG00000198695", "ENSG00000198712", "ENSG00000198804",
    "ENSG00000198763", "ENSG00000228253", "ENSG00000198938",
    "ENSG00000198840")
```


```r
reads <- calculateQCMetrics(
    reads,
    feature_controls = list(
        ERCC = isSpike(reads, "ERCC"), 
        MT = isSpike(reads, "MT")
    )
)
```


```r
hist(
    reads$total_counts,
    breaks = 100
)
abline(v = 1.3e6, col = "red")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{15-exprs-qc-reads_files/figure-latex/total-counts-hist-reads-1} 

}

\caption{Histogram of library sizes for all cells}(\#fig:total-counts-hist-reads)
\end{figure}


```r
filter_by_total_counts <- (reads$total_counts > 1.3e6)
```


```r
table(filter_by_total_counts)
```

```
## filter_by_total_counts
## FALSE  TRUE 
##   180   684
```


```r
hist(
    reads$total_features,
    breaks = 100
)
abline(v = 7000, col = "red")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{15-exprs-qc-reads_files/figure-latex/total-features-hist-reads-1} 

}

\caption{Histogram of the number of detected genes in all cells}(\#fig:total-features-hist-reads)
\end{figure}


```r
filter_by_expr_features <- (reads$total_features > 7000)
```


```r
table(filter_by_expr_features)
```

```
## filter_by_expr_features
## FALSE  TRUE 
##   116   748
```


```r
plotPhenoData(
    reads,
    aes_string(
        x = "total_features",
        y = "pct_counts_MT",
        colour = "batch"
    )
)
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{15-exprs-qc-reads_files/figure-latex/mt-vs-counts-reads-1} 

}

\caption{Percentage of counts in MT genes}(\#fig:mt-vs-counts-reads)
\end{figure}


```r
plotPhenoData(
    reads,
    aes_string(
        x = "total_features",
        y = "pct_counts_ERCC",
        colour = "batch"
    )
)
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{15-exprs-qc-reads_files/figure-latex/ercc-vs-counts-reads-1} 

}

\caption{Percentage of counts in ERCCs}(\#fig:ercc-vs-counts-reads)
\end{figure}


```r
filter_by_ERCC <- 
    reads$batch != "NA19098.r2" & reads$pct_counts_ERCC < 25
table(filter_by_ERCC)
```

```
## filter_by_ERCC
## FALSE  TRUE 
##   103   761
```

```r
filter_by_MT <- reads$pct_counts_MT < 30
table(filter_by_MT)
```

```
## filter_by_MT
## FALSE  TRUE 
##    18   846
```


```r
reads$use <- (
    # sufficient features (genes)
    filter_by_expr_features &
    # sufficient molecules counted
    filter_by_total_counts &
    # sufficient endogenous RNA
    filter_by_ERCC &
    # remove cells with unusual number of reads in MT genes
    filter_by_MT
)
```


```r
table(reads$use)
```

```
## 
## FALSE  TRUE 
##   258   606
```


```r
reads <- plotPCA(
    reads,
    size_by = "total_features", 
    shape_by = "use",
    pca_data_input = "pdata",
    detect_outliers = TRUE,
    return_SCE = TRUE
)
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{15-exprs-qc-reads_files/figure-latex/auto-cell-filt-reads-1} 

}

\caption{PCA plot used for automatic detection of cell outliers}(\#fig:auto-cell-filt-reads)
\end{figure}


```r
table(reads$outlier)
```

```
## 
## FALSE  TRUE 
##   756   108
```


```r
library(limma)
```

```
## 
## Attaching package: 'limma'
```

```
## The following object is masked from 'package:scater':
## 
##     plotMDS
```

```
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
```

```r
auto <- colnames(reads)[reads$outlier]
man <- colnames(reads)[!reads$use]
venn.diag <- vennCounts(
    cbind(colnames(reads) %in% auto,
    colnames(reads) %in% man)
)
vennDiagram(
    venn.diag,
    names = c("Automatic", "Manual"),
    circle.col = c("blue", "green")
)
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{15-exprs-qc-reads_files/figure-latex/cell-filt-comp-reads-1} 

}

\caption{Comparison of the default, automatic and manual cell filters}(\#fig:cell-filt-comp-reads)
\end{figure}


```r
plotQC(reads, type = "highest-expression")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{15-exprs-qc-reads_files/figure-latex/top50-gene-expr-reads-1} 

}

\caption{Number of total counts consumed by the top 50 expressed genes}(\#fig:top50-gene-expr-reads)
\end{figure}


```r
filter_genes <- apply(
    counts(reads[, colData(reads)$use]), 
    1, 
    function(x) length(x[x > 1]) >= 2
)
rowData(reads)$use <- filter_genes
```


```r
table(filter_genes)
```

```
## filter_genes
## FALSE  TRUE 
##  2664 16062
```


```r
dim(reads[rowData(reads)$use, colData(reads)$use])
```

```
## [1] 16062   606
```


```r
assay(reads, "logcounts_raw") <- log2(counts(reads) + 1)
reducedDim(reads) <- NULL
```


```r
saveRDS(reads, file = "tung/reads.rds")
```

By comparing Figure \@ref(fig:cell-filt-comp) and Figure \@ref(fig:cell-filt-comp-reads), it is clear that the reads based filtering removed more cells than the UMI based analysis. If you go back and compare the results you should be able to conclude that the ERCC and MT filters are more strict for the reads-based analysis.


```r
sessionInfo()
```

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
## [1] parallel  stats4    methods   stats     graphics  grDevices utils    
## [8] datasets  base     
## 
## other attached packages:
##  [1] limma_3.34.9               scater_1.6.3              
##  [3] ggplot2_2.2.1              SingleCellExperiment_1.0.0
##  [5] SummarizedExperiment_1.8.1 DelayedArray_0.4.1        
##  [7] matrixStats_0.53.1         Biobase_2.38.0            
##  [9] GenomicRanges_1.30.3       GenomeInfoDb_1.14.0       
## [11] IRanges_2.12.0             S4Vectors_0.16.0          
## [13] BiocGenerics_0.24.0        knitr_1.20                
## 
## loaded via a namespace (and not attached):
##   [1] ggbeeswarm_0.6.0       minqa_1.2.4            colorspace_1.3-2      
##   [4] mvoutlier_2.0.9        rjson_0.2.15           modeltools_0.2-21     
##   [7] class_7.3-14           mclust_5.4             rprojroot_1.3-2       
##  [10] XVector_0.18.0         pls_2.6-0              cvTools_0.3.2         
##  [13] MatrixModels_0.4-1     flexmix_2.3-14         bit64_0.9-7           
##  [16] AnnotationDbi_1.40.0   mvtnorm_1.0-7          sROC_0.1-2            
##  [19] splines_3.4.3          tximport_1.6.0         robustbase_0.92-8     
##  [22] nloptr_1.0.4           robCompositions_2.0.6  pbkrtest_0.4-7        
##  [25] kernlab_0.9-25         cluster_2.0.6          shinydashboard_0.6.1  
##  [28] shiny_1.0.5            rrcov_1.4-3            compiler_3.4.3        
##  [31] httr_1.3.1             backports_1.1.2        assertthat_0.2.0      
##  [34] Matrix_1.2-7.1         lazyeval_0.2.1         htmltools_0.3.6       
##  [37] quantreg_5.35          prettyunits_1.0.2      tools_3.4.3           
##  [40] bindrcpp_0.2           gtable_0.2.0           glue_1.2.0            
##  [43] GenomeInfoDbData_1.0.0 reshape2_1.4.3         dplyr_0.7.4           
##  [46] Rcpp_0.12.15           trimcluster_0.1-2      sgeostat_1.0-27       
##  [49] nlme_3.1-129           fpc_2.1-11             lmtest_0.9-35         
##  [52] laeken_0.4.6           xfun_0.1               stringr_1.3.0         
##  [55] lme4_1.1-15            mime_0.5               XML_3.98-1.10         
##  [58] edgeR_3.20.9           DEoptimR_1.0-8         zoo_1.8-1             
##  [61] zlibbioc_1.24.0        MASS_7.3-45            scales_0.5.0          
##  [64] VIM_4.7.0              SparseM_1.77           rhdf5_2.22.0          
##  [67] RColorBrewer_1.1-2     yaml_2.1.17            memoise_1.1.0         
##  [70] gridExtra_2.3          biomaRt_2.34.2         reshape_0.8.7         
##  [73] stringi_1.1.6          RSQLite_2.0            pcaPP_1.9-73          
##  [76] e1071_1.6-8            boot_1.3-18            prabclus_2.2-6        
##  [79] rlang_0.2.0            pkgconfig_2.0.1        bitops_1.0-6          
##  [82] evaluate_0.10.1        lattice_0.20-34        bindr_0.1             
##  [85] labeling_0.3           cowplot_0.9.2          bit_1.1-12            
##  [88] GGally_1.3.2           plyr_1.8.4             magrittr_1.5          
##  [91] bookdown_0.7           R6_2.2.2               DBI_0.7               
##  [94] pillar_1.2.1           mgcv_1.8-23            RCurl_1.95-4.10       
##  [97] sp_1.2-7               nnet_7.3-12            tibble_1.4.2          
## [100] car_2.1-6              rmarkdown_1.8          viridis_0.5.0         
## [103] progress_1.1.2         locfit_1.5-9.1         grid_3.4.3            
## [106] data.table_1.10.4-3    blob_1.1.0             diptest_0.75-7        
## [109] vcd_1.4-4              digest_0.6.15          xtable_1.8-2          
## [112] httpuv_1.3.6.1         munsell_0.4.3          beeswarm_0.2.3        
## [115] viridisLite_0.3.0      vipor_0.4.5
```
