---
output: html_document
---

## Imputation



```r
library(scImpute)
library(SC3)
library(scater)
library(SingleCellExperiment)
library(mclust)
set.seed(1234567)
```


As discussed previously, one of the main challenges when analyzing scRNA-seq data is the presence of zeros, or dropouts. The dropouts are assumed to have arisen for three possible reasons:

* The gene was not expressed in the cell and hence there are no transcripts to sequence
* The gene was expressed, but for some reason the transcripts were lost somewhere prior to sequencing
* The gene was expressed and transcripts were captured and turned into cDNA, but the sequencing depth was not sufficient to produce any reads.

Thus, dropouts could be result of experimental shortcomings, and if this is the case then we would like to provide computational corrections. One possible solution is to impute the dropouts in the expression matrix. To be able to impute gene expression values, one must have an underlying model. However, since we do not know which dropout events are technical artefacts and which correspond to the transcript being truly absent, imputation is a difficult challenges.

To the best of our knowledge, there are currently two different imputation methods available: MAGIC [@Van_Dijk2017-bh] and scImpute [@Li2017-tz]. [MAGIC](https://github.com/pkathail/magic) is only available for Python or Matlab, but we will run it from within R.

### scImpute

To test `scImpute`, we use the default parameters and we apply it to the Deng dataset that we have worked with before. scImpute takes a .csv or .txt file as an input:


```r
deng <- readRDS("deng/deng-reads.rds")
write.csv(counts(deng), "deng.csv")
scimpute(
    count_path = "deng.csv",
    infile = "csv",
    outfile = "txt", 
    out_dir = "./",
    Kcluster = 10,
    ncores = 2
)
```

```
## [1] "reading in raw count matrix ..."
## [1] "number of genes in raw count matrix 22431"
## [1] "number of cells in raw count matrix 268"
## [1] "inferring cell similarities ..."
## [1] "cluster sizes:"
##  [1] 12  9 26  5  9 57 58 43 17 22
## [1] "estimating dropout probability for type 1 ..."
## [1] "imputing dropout values for type 1 ..."
## [1] "estimating dropout probability for type 2 ..."
## [1] "imputing dropout values for type 2 ..."
## [1] "estimating dropout probability for type 3 ..."
## [1] "imputing dropout values for type 3 ..."
## [1] "estimating dropout probability for type 4 ..."
## [1] "imputing dropout values for type 4 ..."
## [1] "estimating dropout probability for type 5 ..."
## [1] "imputing dropout values for type 5 ..."
## [1] "estimating dropout probability for type 6 ..."
## [1] "imputing dropout values for type 6 ..."
## [1] "estimating dropout probability for type 7 ..."
## [1] "imputing dropout values for type 7 ..."
## [1] "estimating dropout probability for type 8 ..."
## [1] "imputing dropout values for type 8 ..."
## [1] "estimating dropout probability for type 9 ..."
## [1] "imputing dropout values for type 9 ..."
## [1] "estimating dropout probability for type 10 ..."
## [1] "imputing dropout values for type 10 ..."
## [1] "writing imputed count matrix ..."
```

```
##  [1]  17  18  88 111 126 177 186 229 244 247
```

Now we can compare the results with original data by considering a PCA plot


```r
res <- read.table("scimpute_count.txt")
colnames(res) <- NULL
res <- SingleCellExperiment(
    assays = list(logcounts = log2(as.matrix(res) + 1)), 
    colData = colData(deng)
)
plotPCA(
    res, 
    colour_by = "cell_type2"
)
```

<img src="21-imputation_files/figure-html/unnamed-chunk-4-1.png" width="672" style="display: block; margin: auto;" />

Compare this result to the original data in Chapter \@ref(clust-methods). What are the most significant differences?

To evaluate the impact of the imputation, we use `SC3` to cluster the imputed matrix

```r
res <- sc3_estimate_k(res)
metadata(res)$sc3$k_estimation
```

```
## [1] 6
```

```r
res <- sc3(res, ks = 10, gene_filter = FALSE)
```

```r
adjustedRandIndex(colData(deng)$cell_type2, colData(res)$sc3_10_clusters)
```

```
## [1] 0.4687332
```

```r
plotPCA(
    res, 
    colour_by = "sc3_10_clusters"
)
```

<img src="21-imputation_files/figure-html/unnamed-chunk-5-1.png" width="672" style="display: block; margin: auto;" />

__Exercise:__ Based on the PCA and the clustering results, do you think that imputation using `scImpute` is a good idea for the Deng dataset?

### MAGIC


```r
system("python3 utils/MAGIC.py -d deng.csv -o MAGIC_count.csv --cell-axis columns -l 1 --pca-non-random csv")
```


```r
res <- read.csv("MAGIC_count.csv", header = TRUE)
rownames(res) <- res[,1]
res <- res[,-1]
res <- t(res)
res <- SingleCellExperiment(
    assays = list(logcounts = res), 
    colData = colData(deng)
)
plotPCA(
    res, 
    colour_by = "cell_type2"
)
```

<img src="21-imputation_files/figure-html/unnamed-chunk-7-1.png" width="672" style="display: block; margin: auto;" />

Compare this result to the original data in Chapter \@ref(clust-methods). What are the most significant differences?

To evaluate the impact of the imputation, we use `SC3` to cluster the imputed matrix

```r
res <- sc3_estimate_k(res)
metadata(res)$sc3$k_estimation
```

```
## [1] 4
```

```r
res <- sc3(res, ks = 10, gene_filter = FALSE)
```

```r
adjustedRandIndex(colData(deng)$cell_type2, colData(res)$sc3_10_clusters)
```

```
## [1] 0.3752866
```

```r
plotPCA(
    res, 
    colour_by = "sc3_10_clusters"
)
```

<img src="21-imputation_files/figure-html/unnamed-chunk-8-1.png" width="672" style="display: block; margin: auto;" />
__Exercise:__ What is the difference between `scImpute` and `MAGIC` based on the PCA and clustering analysis? Which one do you think is best to use?


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
## [1] stats4    methods   parallel  stats     graphics  grDevices utils    
## [8] datasets  base     
## 
## other attached packages:
##  [1] mclust_5.3                  scater_1.5.21              
##  [3] ggplot2_2.2.1               SC3_1.5.6                  
##  [5] SingleCellExperiment_0.99.4 SummarizedExperiment_1.6.5 
##  [7] DelayedArray_0.2.7          matrixStats_0.52.2         
##  [9] Biobase_2.36.2              GenomicRanges_1.28.6       
## [11] GenomeInfoDb_1.12.3         IRanges_2.10.5             
## [13] S4Vectors_0.14.7            BiocGenerics_0.22.1        
## [15] scImpute_0.0.4              doParallel_1.0.11          
## [17] iterators_1.0.8             foreach_1.4.3              
## [19] penalized_0.9-50            survival_2.40-1            
## [21] kernlab_0.9-25              knitr_1.17                 
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6            bit64_0.9-7            
##  [3] RColorBrewer_1.1-2      rprojroot_1.2          
##  [5] tools_3.4.2             backports_1.1.1        
##  [7] doRNG_1.6.6             R6_2.2.2               
##  [9] vipor_0.4.5             KernSmooth_2.23-15     
## [11] DBI_0.7                 lazyeval_0.2.0         
## [13] colorspace_1.3-2        gridExtra_2.3          
## [15] bit_1.1-12              compiler_3.4.2         
## [17] pkgmaker_0.22           labeling_0.3           
## [19] bookdown_0.5            caTools_1.17.1         
## [21] scales_0.5.0            DEoptimR_1.0-8         
## [23] mvtnorm_1.0-6           robustbase_0.92-7      
## [25] stringr_1.2.0           digest_0.6.12          
## [27] rmarkdown_1.6           XVector_0.16.0         
## [29] pkgconfig_2.0.1         rrcov_1.4-3            
## [31] htmltools_0.3.6         WriteXLS_4.0.0         
## [33] limma_3.32.10           rlang_0.1.2            
## [35] RSQLite_2.0             shiny_1.0.5            
## [37] bindr_0.1               gtools_3.5.0           
## [39] dplyr_0.7.4             RCurl_1.95-4.8         
## [41] magrittr_1.5            GenomeInfoDbData_0.99.0
## [43] Matrix_1.2-7.1          ggbeeswarm_0.6.0       
## [45] Rcpp_0.12.13            munsell_0.4.3          
## [47] viridis_0.4.0           edgeR_3.18.1           
## [49] stringi_1.1.5           yaml_2.1.14            
## [51] zlibbioc_1.22.0         rhdf5_2.20.0           
## [53] gplots_3.0.1            plyr_1.8.4             
## [55] grid_3.4.2              blob_1.1.0             
## [57] gdata_2.18.0            shinydashboard_0.6.1   
## [59] lattice_0.20-34         cowplot_0.8.0          
## [61] splines_3.4.2           locfit_1.5-9.1         
## [63] rjson_0.2.15            rngtools_1.2.4         
## [65] reshape2_1.4.2          codetools_0.2-15       
## [67] biomaRt_2.32.1          glue_1.1.1             
## [69] XML_3.98-1.9            evaluate_0.10.1        
## [71] data.table_1.10.4-3     httpuv_1.3.5           
## [73] gtable_0.2.0            assertthat_0.2.0       
## [75] mime_0.5                xtable_1.8-2           
## [77] e1071_1.6-8             class_7.3-14           
## [79] pcaPP_1.9-72            viridisLite_0.2.0      
## [81] tibble_1.3.4            pheatmap_1.0.8         
## [83] beeswarm_0.2.3          AnnotationDbi_1.38.2   
## [85] registry_0.3            memoise_1.1.0          
## [87] tximport_1.4.0          bindrcpp_0.2           
## [89] cluster_2.0.6           ROCR_1.0-7
```
