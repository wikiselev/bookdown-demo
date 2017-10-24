---
output: html_document
---




```r
set.seed(1234567)
```

## `Seurat`

[Seurat](http://satijalab.org/seurat/) is a popular R package that is designed for QC, analysis, and exploration of single cell RNA-seq data, i.e. many of the tasks covered in this course. Although the authors provide several [tutorials](http://satijalab.org/seurat/get_started.html), here we provide a brief overview by following an [example](http://satijalab.org/seurat/pbmc3k_tutorial.html) created by the authors of `Seurat` (2,800 Peripheral Blood Mononuclear Cells). We mostly use default values in various function calls, for more details please consult the documentation and the authors. We start by loading the Deng data that we have used before:

```r
deng <- readRDS("deng/deng-reads.rds")
```

First, we initialize the Seurat object with the raw (non-normalized data). Keep all genes expressed in >= 3 cells. Keep all cells with at least 200 detected genes:

```r
library(SingleCellExperiment)
library(Seurat)
seuset <- CreateSeuratObject(
    raw.data = counts(deng),
    min.cells = 3, 
    min.genes = 200
)
```


```r
VlnPlot(object = seuset, features.plot = c("nGene", "nUMI"), nCol = 3)
```

<img src="19-seurat_files/figure-html/unnamed-chunk-5-1.png" width="672" style="display: block; margin: auto;" />

```r
GenePlot(object = seuset, gene1 = "nUMI", gene2 = "nGene")
```

<img src="19-seurat_files/figure-html/unnamed-chunk-5-2.png" width="672" style="display: block; margin: auto;" />

```r
seuset <- NormalizeData(object = seuset, normalization.method = "LogNormalize", 
    scale.factor = 10000)

seuset <- FindVariableGenes(
    object = seuset
)
```

<img src="19-seurat_files/figure-html/unnamed-chunk-5-3.png" width="672" style="display: block; margin: auto;" />

```r
length(x = seuset@var.genes)
```

```
## [1] 1459
```

```r
seuset <- ScaleData(object = seuset, vars.to.regress = c("nUMI"))
```

```
## [1] "Regressing out nUMI"
##   |                                                                         |                                                                 |   0%  |                                                                         |=                                                                |   1%  |                                                                         |=                                                                |   2%  |                                                                         |==                                                               |   2%  |                                                                         |==                                                               |   3%  |                                                                         |==                                                               |   4%  |                                                                         |===                                                              |   4%  |                                                                         |===                                                              |   5%  |                                                                         |====                                                             |   6%  |                                                                         |=====                                                            |   7%  |                                                                         |=====                                                            |   8%  |                                                                         |======                                                           |   8%  |                                                                         |======                                                           |   9%  |                                                                         |======                                                           |  10%  |                                                                         |=======                                                          |  10%  |                                                                         |=======                                                          |  11%  |                                                                         |=======                                                          |  12%  |                                                                         |========                                                         |  12%  |                                                                         |========                                                         |  13%  |                                                                         |=========                                                        |  14%  |                                                                         |==========                                                       |  15%  |                                                                         |==========                                                       |  16%  |                                                                         |===========                                                      |  16%  |                                                                         |===========                                                      |  17%  |                                                                         |===========                                                      |  18%  |                                                                         |============                                                     |  18%  |                                                                         |============                                                     |  19%  |                                                                         |=============                                                    |  20%  |                                                                         |==============                                                   |  21%  |                                                                         |==============                                                   |  22%  |                                                                         |===============                                                  |  22%  |                                                                         |===============                                                  |  23%  |                                                                         |===============                                                  |  24%  |                                                                         |================                                                 |  24%  |                                                                         |================                                                 |  25%  |                                                                         |=================                                                |  26%  |                                                                         |==================                                               |  27%  |                                                                         |==================                                               |  28%  |                                                                         |===================                                              |  28%  |                                                                         |===================                                              |  29%  |                                                                         |===================                                              |  30%  |                                                                         |====================                                             |  30%  |                                                                         |====================                                             |  31%  |                                                                         |====================                                             |  32%  |                                                                         |=====================                                            |  32%  |                                                                         |=====================                                            |  33%  |                                                                         |======================                                           |  34%  |                                                                         |=======================                                          |  35%  |                                                                         |=======================                                          |  36%  |                                                                         |========================                                         |  36%  |                                                                         |========================                                         |  37%  |                                                                         |========================                                         |  38%  |                                                                         |=========================                                        |  38%  |                                                                         |=========================                                        |  39%  |                                                                         |==========================                                       |  40%  |                                                                         |===========================                                      |  41%  |                                                                         |===========================                                      |  42%  |                                                                         |============================                                     |  42%  |                                                                         |============================                                     |  43%  |                                                                         |============================                                     |  44%  |                                                                         |=============================                                    |  44%  |                                                                         |=============================                                    |  45%  |                                                                         |==============================                                   |  46%  |                                                                         |===============================                                  |  47%  |                                                                         |===============================                                  |  48%  |                                                                         |================================                                 |  48%  |                                                                         |================================                                 |  49%  |                                                                         |================================                                 |  50%  |                                                                         |=================================                                |  50%  |                                                                         |=================================                                |  51%  |                                                                         |=================================                                |  52%  |                                                                         |==================================                               |  52%  |                                                                         |==================================                               |  53%  |                                                                         |===================================                              |  54%  |                                                                         |====================================                             |  55%  |                                                                         |====================================                             |  56%  |                                                                         |=====================================                            |  56%  |                                                                         |=====================================                            |  57%  |                                                                         |=====================================                            |  58%  |                                                                         |======================================                           |  58%  |                                                                         |======================================                           |  59%  |                                                                         |=======================================                          |  60%  |                                                                         |========================================                         |  61%  |                                                                         |========================================                         |  62%  |                                                                         |=========================================                        |  62%  |                                                                         |=========================================                        |  63%  |                                                                         |=========================================                        |  64%  |                                                                         |==========================================                       |  64%  |                                                                         |==========================================                       |  65%  |                                                                         |===========================================                      |  66%  |                                                                         |============================================                     |  67%  |                                                                         |============================================                     |  68%  |                                                                         |=============================================                    |  68%  |                                                                         |=============================================                    |  69%  |                                                                         |=============================================                    |  70%  |                                                                         |==============================================                   |  70%  |                                                                         |==============================================                   |  71%  |                                                                         |==============================================                   |  72%  |                                                                         |===============================================                  |  72%  |                                                                         |===============================================                  |  73%  |                                                                         |================================================                 |  74%  |                                                                         |=================================================                |  75%  |                                                                         |=================================================                |  76%  |                                                                         |==================================================               |  76%  |                                                                         |==================================================               |  77%  |                                                                         |==================================================               |  78%  |                                                                         |===================================================              |  78%  |                                                                         |===================================================              |  79%  |                                                                         |====================================================             |  80%  |                                                                         |=====================================================            |  81%  |                                                                         |=====================================================            |  82%  |                                                                         |======================================================           |  82%  |                                                                         |======================================================           |  83%  |                                                                         |======================================================           |  84%  |                                                                         |=======================================================          |  84%  |                                                                         |=======================================================          |  85%  |                                                                         |========================================================         |  86%  |                                                                         |=========================================================        |  87%  |                                                                         |=========================================================        |  88%  |                                                                         |==========================================================       |  88%  |                                                                         |==========================================================       |  89%  |                                                                         |==========================================================       |  90%  |                                                                         |===========================================================      |  90%  |                                                                         |===========================================================      |  91%  |                                                                         |===========================================================      |  92%  |                                                                         |============================================================     |  92%  |                                                                         |============================================================     |  93%  |                                                                         |=============================================================    |  94%  |                                                                         |==============================================================   |  95%  |                                                                         |==============================================================   |  96%  |                                                                         |===============================================================  |  96%  |                                                                         |===============================================================  |  97%  |                                                                         |===============================================================  |  98%  |                                                                         |================================================================ |  98%  |                                                                         |================================================================ |  99%  |                                                                         |=================================================================| 100%
## [1] "Scaling data matrix"
##   |                                                                         |                                                                 |   0%  |                                                                         |=================================================================| 100%
```

```r
seuset <- RunPCA(
    object = seuset, 
    pc.genes = seuset@var.genes, 
    do.print = TRUE, 
    pcs.print = 1:5, 
    genes.print = 5
)
```

```
## [1] "PC1"
## [1] "Lrpap1" "Actg1"  "Snx2"   "Gm2a"   "Psap"  
## [1] ""
## [1] "Gm13023" "Gm10436" "Zbed3"   "Oog1"    "Trim61" 
## [1] ""
## [1] ""
## [1] "PC2"
## [1] "Lrrfip1"  "BC053393" "Id2"      "Krt8"     "Gsta4"   
## [1] ""
## [1] "Gm11517" "Obox6"   "Trim43b" "Trim43c" "Pdxk"   
## [1] ""
## [1] ""
## [1] "PC3"
## [1] "Obox3"   "G2e3"    "Sfi1"    "Gm9125"  "Gm12789"
## [1] ""
## [1] "Fam46c" "Casp8"  "Gja4"   "Pnpla3" "Bhmt"  
## [1] ""
## [1] ""
## [1] "PC4"
## [1] "Upp1"   "Col4a1" "Fn1"    "Baz2b"  "Tdgf1" 
## [1] ""
## [1] "Slc5a6"  "Chka"    "Cldn7"   "Slc4a11" "Gpd1l"  
## [1] ""
## [1] ""
## [1] "PC5"
## [1] "BC051665"      "3110003A17Rik" "Gsta1"         "Acpp"         
## [5] "4930506M07Rik"
## [1] ""
## [1] "Slc20a1" "Abca3"   "Abcd4"   "Cdh1"    "Lrp2"   
## [1] ""
## [1] ""
```

```r
VizPCA(object = seuset, pcs.use = 1:2)
```

<img src="19-seurat_files/figure-html/unnamed-chunk-5-4.png" width="672" style="display: block; margin: auto;" />

```r
PCAPlot(object = seuset, dim.1 = 1, dim.2 = 2)
```

<img src="19-seurat_files/figure-html/unnamed-chunk-5-5.png" width="672" style="display: block; margin: auto;" />

```r
seuset <- ProjectPCA(object = seuset, do.print = FALSE)

PCHeatmap(
    object = seuset, 
    pc.use = 1:6, 
    cells.use = 500, 
    do.balanced = TRUE, 
    label.columns = FALSE,
    use.full = FALSE
)
```

<img src="19-seurat_files/figure-html/unnamed-chunk-5-6.png" width="672" style="display: block; margin: auto;" />

```r
seuset <- JackStraw(
    object = seuset, 
    num.replicate = 100, 
    do.print = FALSE
)

JackStrawPlot(object = seuset, PCs = 1:6)
```

```
## Warning: Removed 6126 rows containing missing values (geom_point).
```

<img src="19-seurat_files/figure-html/unnamed-chunk-5-7.png" width="672" style="display: block; margin: auto;" />

```r
PCElbowPlot(object = seuset)
```

<img src="19-seurat_files/figure-html/unnamed-chunk-5-8.png" width="672" style="display: block; margin: auto;" />

```r
seuset <- FindClusters(
    object = seuset, 
    reduction.type = "pca", 
    dims.use = 1:5, 
    resolution = 0.6, 
    print.output = 0, 
    save.SNN = TRUE
)

table(seuset@ident)
```

```
## 
##  0  1  2  3 
## 88 74 60 46
```

```r
mclust::adjustedRandIndex(colData(deng)$cell_type2, seuset@ident)
```

```
## [1] 0.3820495
```

```r
seuset <- RunTSNE(
    object = seuset,
    dims.use = 1:5,
    do.fast = TRUE
)
```

__Exercise 12__: Compare the results between `SC3` and `Seurat`.

__Our solution__:

```r
colData(deng)$Seurat <- as.character(seuset@ident)
sc3_plot_expression(deng, k = 10, show_pdata = "Seurat")
```


Seurat can also find marker genes, e.g. marker genes for cluster 2:

```r
markers <- FindMarkers(seuset, 2)
FeaturePlot(seuset, 
            head(rownames(markers)), 
            cols.use = c("lightgrey", "blue"), 
            nCol = 3)
```

<img src="19-seurat_files/figure-html/unnamed-chunk-7-1.png" width="672" style="display: block; margin: auto;" />

__Exercise 13__: Compare marker genes provided by `Seurat` and `SC3`.

### sessionInfo()


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
## [1] parallel  stats4    methods   stats     graphics  grDevices utils    
## [8] datasets  base     
## 
## other attached packages:
##  [1] Seurat_2.1.0                Matrix_1.2-11              
##  [3] cowplot_0.8.0               ggplot2_2.2.1              
##  [5] SingleCellExperiment_0.99.4 SummarizedExperiment_1.6.5 
##  [7] DelayedArray_0.2.7          matrixStats_0.52.2         
##  [9] Biobase_2.36.2              GenomicRanges_1.28.6       
## [11] GenomeInfoDb_1.12.3         IRanges_2.10.5             
## [13] S4Vectors_0.14.7            BiocGenerics_0.22.1        
## [15] knitr_1.17                 
## 
## loaded via a namespace (and not attached):
##   [1] backports_1.1.1         Hmisc_4.0-3            
##   [3] VGAM_1.0-4              NMF_0.20.6             
##   [5] sn_1.5-0                plyr_1.8.4             
##   [7] igraph_1.1.2            lazyeval_0.2.0         
##   [9] splines_3.4.2           gridBase_0.4-7         
##  [11] digest_0.6.12           foreach_1.4.3          
##  [13] htmltools_0.3.6         lars_1.2               
##  [15] gdata_2.18.0            magrittr_1.5           
##  [17] checkmate_1.8.5         cluster_2.0.6          
##  [19] doParallel_1.0.11       mixtools_1.1.0         
##  [21] ROCR_1.0-7              sfsmisc_1.1-1          
##  [23] recipes_0.1.0           gower_0.1.2            
##  [25] dimRed_0.1.0            R.utils_2.5.0          
##  [27] colorspace_1.3-2        dplyr_0.7.4            
##  [29] RCurl_1.95-4.8          bindr_0.1              
##  [31] survival_2.41-3         iterators_1.0.8        
##  [33] ape_4.1                 glue_1.1.1             
##  [35] DRR_0.0.2               registry_0.3           
##  [37] gtable_0.2.0            ipred_0.9-6            
##  [39] zlibbioc_1.22.0         XVector_0.16.0         
##  [41] kernlab_0.9-25          ddalpha_1.3.1          
##  [43] prabclus_2.2-6          DEoptimR_1.0-8         
##  [45] scales_0.5.0            mvtnorm_1.0-6          
##  [47] rngtools_1.2.4          Rcpp_0.12.13           
##  [49] dtw_1.18-1              xtable_1.8-2           
##  [51] htmlTable_1.9           tclust_1.3-1           
##  [53] foreign_0.8-69          proxy_0.4-17           
##  [55] mclust_5.3              SDMTools_1.1-221       
##  [57] Formula_1.2-2           tsne_0.1-3             
##  [59] lava_1.5.1              prodlim_1.6.1          
##  [61] htmlwidgets_0.9         FNN_1.1                
##  [63] gplots_3.0.1            RColorBrewer_1.1-2     
##  [65] fpc_2.1-10              acepack_1.4.1          
##  [67] modeltools_0.2-21       ica_1.0-1              
##  [69] pkgconfig_2.0.1         R.methodsS3_1.7.1      
##  [71] flexmix_2.3-14          nnet_7.3-12            
##  [73] caret_6.0-77            labeling_0.3           
##  [75] rlang_0.1.2             reshape2_1.4.2         
##  [77] munsell_0.4.3           tools_3.4.2            
##  [79] ranger_0.8.0            ggridges_0.4.1         
##  [81] evaluate_0.10.1         stringr_1.2.0          
##  [83] yaml_2.1.14             ModelMetrics_1.1.0     
##  [85] robustbase_0.92-7       caTools_1.17.1         
##  [87] purrr_0.2.4             bindrcpp_0.2           
##  [89] pbapply_1.3-3           nlme_3.1-131           
##  [91] R.oo_1.21.0             RcppRoll_0.2.2         
##  [93] compiler_3.4.2          ggjoy_0.4.0            
##  [95] tibble_1.3.4            stringi_1.1.5          
##  [97] lattice_0.20-35         trimcluster_0.1-2      
##  [99] diffusionMap_1.1-0      data.table_1.10.4-2    
## [101] bitops_1.0-6            irlba_2.3.1            
## [103] R6_2.2.2                latticeExtra_0.6-28    
## [105] bookdown_0.5            KernSmooth_2.23-15     
## [107] gridExtra_2.3           codetools_0.2-15       
## [109] MASS_7.3-47             gtools_3.5.0           
## [111] assertthat_0.2.0        CVST_0.2-1             
## [113] pkgmaker_0.22           rprojroot_1.2          
## [115] withr_2.0.0             mnormt_1.5-5           
## [117] GenomeInfoDbData_0.99.0 diptest_0.75-7         
## [119] grid_3.4.2              rpart_4.1-11           
## [121] timeDate_3012.100       tidyr_0.7.2            
## [123] class_7.3-14            rmarkdown_1.6          
## [125] segmented_0.5-2.2       Rtsne_0.13             
## [127] numDeriv_2016.8-1       scatterplot3d_0.3-40   
## [129] lubridate_1.6.0         base64enc_0.1-3
```
