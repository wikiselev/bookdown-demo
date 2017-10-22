---
output: html_document
---

# Imputation



```r
library(scImpute)
library(SC3)
library(scater)
library(SingleCellExperiment)
set.seed(1234567)
```


As discussed previously, one of the main challenges when analyzing scRNA-seq data is the presence of zeros, or dropouts. The dropouts are assumed to have arisen for three possible reasons:

* The gene was not expressed in the cell and hence there are no transcripts to sequence
* The gene was expressed, but for some reason the transcripts were lost somewhere prior to sequencing
* The gene was expressed and transcripts were captured and turned into cDNA, but the sequencing depth was not sufficient to produce any reads.

Thus, dropouts could be result of experimental shortcomings, and if this is the case then we would like to provide computational corrections. One possible solution is to impute the dropouts in the expression matrix. To be able to impute gene expression values, one must have an underlying model. However, since we do not know which dropout events are technical artefacts and which correspond to the transcript being truly absent, imputation is a difficult challenges.

To the best of our knowledge, there are currently two different imputation methods available: MAGIC [@vanDijk2017] and scImpute [@Li2017]. Since [MAGIC](https://github.com/pkathail/magic) is only available for Python or Matlab, we will only cover [scImpute](https://github.com/Vivianstats/scImpute) in this chapter.

To test scImpute, we use the default parameters and we apply it to the Pollen dataset that we have worked with before. scImpute takes a .csv or .txt file as an input:


```
## [1] "reading in raw count matrix ..."
## [1] "estimating mixture models ..."
## [1] "imputing dropout values ..."
## [1] "writing imputed count matrix ..."
```

```
## [1] 0
```

Now we can compare the results with original data by considering a PCA plot

<img src="22-imputation_files/figure-html/unnamed-chunk-4-1.png" width="672" style="display: block; margin: auto;" />

Compare this result to the original data in Fig. X, Chapter Y. What are the most significant differences?

To evaluate the impact of the imputation, we use `SC3` to cluster the imputed matrix

```
## Estimating k...
```

```
## [1] 8
```

```
## Setting SC3 parameters...
```

```
## Setting a range of k...
```

```
## Calculating distances between the cells...
```

```
## Performing transformations and calculating eigenvectors...
```

```
## Performing k-means clustering...
```

```
## Calculating consensus matrix...
```

<img src="22-imputation_files/figure-html/unnamed-chunk-5-1.png" width="672" style="display: block; margin: auto;" /><img src="22-imputation_files/figure-html/unnamed-chunk-5-2.png" width="672" style="display: block; margin: auto;" />

Based on the PCA and the clustering results, do you think that imputation is a good idea for the Pollen dataset?

## sessionInfo()


```
## R version 3.4.2 (2017-09-28)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Debian GNU/Linux buster/sid
## 
## Matrix products: default
## BLAS: /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats4    methods   parallel  stats     graphics  grDevices utils    
## [8] datasets  base     
## 
## other attached packages:
##  [1] scater_1.5.20               ggplot2_2.2.1              
##  [3] SC3_1.5.5                   SingleCellExperiment_0.99.4
##  [5] SummarizedExperiment_1.6.5  DelayedArray_0.2.7         
##  [7] matrixStats_0.52.2          Biobase_2.36.2             
##  [9] GenomicRanges_1.28.6        GenomeInfoDb_1.12.3        
## [11] IRanges_2.10.5              S4Vectors_0.14.7           
## [13] BiocGenerics_0.22.1         scImpute_0.0.2             
## [15] glmnet_2.0-13               foreach_1.4.3              
## [17] Matrix_1.2-11               knitr_1.17                 
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6            bit64_0.9-7            
##  [3] doParallel_1.0.11       RColorBrewer_1.1-2     
##  [5] rprojroot_1.2           tools_3.4.2            
##  [7] backports_1.1.1         doRNG_1.6.6            
##  [9] R6_2.2.2                vipor_0.4.5            
## [11] KernSmooth_2.23-15      DBI_0.7                
## [13] lazyeval_0.2.0          colorspace_1.3-2       
## [15] gridExtra_2.3           bit_1.1-12             
## [17] compiler_3.4.2          pkgmaker_0.22          
## [19] labeling_0.3            bookdown_0.5           
## [21] caTools_1.17.1          scales_0.5.0           
## [23] DEoptimR_1.0-8          mvtnorm_1.0-6          
## [25] robustbase_0.92-7       stringr_1.2.0          
## [27] digest_0.6.12           rmarkdown_1.6          
## [29] XVector_0.16.0          pkgconfig_2.0.1        
## [31] rrcov_1.4-3             htmltools_0.3.6        
## [33] WriteXLS_4.0.0          limma_3.32.10          
## [35] rlang_0.1.2             RSQLite_2.0            
## [37] shiny_1.0.5             bindr_0.1              
## [39] gtools_3.5.0            dplyr_0.7.4            
## [41] RCurl_1.95-4.8          magrittr_1.5           
## [43] GenomeInfoDbData_0.99.0 ggbeeswarm_0.6.0       
## [45] Rcpp_0.12.13            munsell_0.4.3          
## [47] viridis_0.4.0           edgeR_3.18.1           
## [49] stringi_1.1.5           yaml_2.1.14            
## [51] zlibbioc_1.22.0         rhdf5_2.20.0           
## [53] gplots_3.0.1            plyr_1.8.4             
## [55] grid_3.4.2              blob_1.1.0             
## [57] gdata_2.18.0            shinydashboard_0.6.1   
## [59] lattice_0.20-35         cowplot_0.8.0          
## [61] locfit_1.5-9.1          rjson_0.2.15           
## [63] rngtools_1.2.4          reshape2_1.4.2         
## [65] codetools_0.2-15        biomaRt_2.32.1         
## [67] glue_1.1.1              XML_3.98-1.9           
## [69] evaluate_0.10.1         data.table_1.10.4-2    
## [71] httpuv_1.3.5            gtable_0.2.0           
## [73] assertthat_0.2.0        mime_0.5               
## [75] xtable_1.8-2            e1071_1.6-8            
## [77] class_7.3-14            pcaPP_1.9-72           
## [79] viridisLite_0.2.0       tibble_1.3.4           
## [81] pheatmap_1.0.8          iterators_1.0.8        
## [83] beeswarm_0.2.3          AnnotationDbi_1.38.2   
## [85] registry_0.3            memoise_1.1.0          
## [87] tximport_1.4.0          bindrcpp_0.2           
## [89] cluster_2.0.6           ROCR_1.0-7
```
