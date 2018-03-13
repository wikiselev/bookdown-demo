---
output: html_document
---

## Search scRNA-Seq data



```r
library(scfind)
library(SingleCellExperiment)
set.seed(1234567)
```

### About 

`scfind` is a tool that allows one to search single cell RNA-Seq collections (Atlas) using lists of genes, e.g. searching for cells and cell-types where a specific set of genes are expressed. `scfind` is a [Bioconductor package](http://bioconductor.org/packages/scfind). Cloud implementation of `scfind` with a large collection of datasets is available on our [website](http://www.hemberg-lab.cloud/scfind).

<div class="figure" style="text-align: center">
<img src="figures/scfind.png" alt="scfind can be used to search large collection of scRNA-seq data by a list of gene IDs." width="80%" />
<p class="caption">(\#fig:unnamed-chunk-3)scfind can be used to search large collection of scRNA-seq data by a list of gene IDs.</p>
</div>

### Dataset

We will run `scfind` on the same human pancreas dataset as in the previous chapter. `scfind` also operates on `SingleCellExperiment` class:

```r
muraro <- readRDS("pancreas/muraro.rds")
```

### Gene Index

Now we need to create a gene index using our dataset:

```r
cellIndex <- buildCellIndex(muraro)
```

The gene index contains for each gene indexes of the cells where it is expressed. This is similar to sparsification of the expression matrix. In addition to this the index is also compressed in a way that it can accessed very quickly. We estimated that one can achieve 5-10 compression factor with this method.

By default the `cell_type1` column of the `colData` slot of the `SingleCellExperiment` object is used to define cell types, however it can also be defined manually using the `cell_type_column` argument of the `buildCellTypeIndex` function (check `?buildCellTypeIndex`).

### Marker genes

Now let's define lists of cell type specific marker genes. We will use the marker genes identified in the original publication, namely in Figure 1:

```r
# these genes are taken from fig. 1
muraro_alpha <- c("GCG", "LOXL4", "PLCE1", "IRX2", "GC", "KLHL41", 
                  "CRYBA2", "TTR", "TM4SF4", "RGS4")
muraro_beta <- c("INS", "IAPP", "MAFA", "NPTX2", "DLK1", "ADCYAP1", 
                 "PFKFB2", "PDX1", "TGFBR3", "SYT13")
muraro_gamma <- c("PPY", "SERTM1", "CARTPT", "SLITRK6", "ETV1", 
                  "THSD7A", "AQP3", "ENTPD2", "PTGFR", "CHN2")
muraro_delta <- c("SST", "PRG4", "LEPR", "RBP4", "BCHE", "HHEX", 
                  "FRZB", "PCSK1", "RGS2", "GABRG2")
```

### Search cells by a gene list

`findCell` function returns a list of p-values corresponding to all cell types in a given dataset. It also outputs a list of cells in which genes from the given gene list are co-expressed. We will run it on all lists of marker genes defined above:

```r
res <- findCell(cellIndex, muraro_alpha)
barplot(-log10(res$p_values), ylab = "-log10(pval)", las = 2)
```

<img src="32-search_files/figure-html/unnamed-chunk-7-1.png" width="672" style="display: block; margin: auto;" />

```r
head(res$common_exprs_cells)
```

```
##   cell_id cell_type
## 1       1     alpha
## 2       3     alpha
## 3       7     alpha
## 4       9     alpha
## 5      15     alpha
## 6      20     alpha
```

__Exercise 1__

Perform a search by _beta_, _delta_ and _gamma_ gene lists and explore the results.

<img src="32-search_files/figure-html/unnamed-chunk-8-1.png" width="672" style="display: block; margin: auto;" />

```
##   cell_id cell_type
## 1      71     alpha
## 2      72      beta
## 3      92      beta
## 4      96      beta
## 5      98      beta
## 6     102      beta
```

<img src="32-search_files/figure-html/unnamed-chunk-8-2.png" width="672" style="display: block; margin: auto;" />

```
##   cell_id cell_type
## 1      40     delta
## 2     212     delta
## 3     225     delta
## 4     253     delta
## 5     330     delta
## 6     400     delta
```

<img src="32-search_files/figure-html/unnamed-chunk-8-3.png" width="672" style="display: block; margin: auto;" />

```
##   cell_id cell_type
## 1      53     alpha
## 2     102      beta
## 3     255     gamma
## 4     305     gamma
## 5     525     gamma
## 6     662     gamma
```


__Exercise 2__

Load the `segerstolpe` and search it using _alpha_, _beta_, _delta_ and _gamma_ gene lists identified in `muraro` dataset.

<img src="32-search_files/figure-html/unnamed-chunk-9-1.png" width="672" style="display: block; margin: auto;" />

```
##   cell_id cell_type
## 1      18     alpha
## 2      20     alpha
## 3      24     alpha
## 4      32     alpha
## 5      43     alpha
## 6      48     alpha
```

<img src="32-search_files/figure-html/unnamed-chunk-9-2.png" width="672" style="display: block; margin: auto;" />

```
##   cell_id     cell_type
## 1      15 co-expression
## 2      58          beta
## 3     300          beta
## 4     390 co-expression
## 5     504 co-expression
## 6     506          beta
```

<img src="32-search_files/figure-html/unnamed-chunk-9-3.png" width="672" style="display: block; margin: auto;" />

```
##   cell_id     cell_type
## 1     170         delta
## 2     715         delta
## 3    1039 co-expression
## 4    1133         delta
## 5    1719         delta
## 6    1721         delta
```

<img src="32-search_files/figure-html/unnamed-chunk-9-4.png" width="672" style="display: block; margin: auto;" />

```
##   cell_id cell_type
## 1      47     gamma
## 2     458     gamma
## 3     476     gamma
## 4     600     gamma
## 5     606     gamma
## 6     622     gamma
```

### sessionInfo()


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
##  [1] SingleCellExperiment_1.0.0 SummarizedExperiment_1.8.1
##  [3] DelayedArray_0.4.1         matrixStats_0.53.1        
##  [5] Biobase_2.38.0             GenomicRanges_1.30.3      
##  [7] GenomeInfoDb_1.14.0        IRanges_2.12.0            
##  [9] S4Vectors_0.16.0           BiocGenerics_0.24.0       
## [11] scfind_1.0.0               knitr_1.20                
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.15           highr_0.6              plyr_1.8.4            
##  [4] pillar_1.2.1           compiler_3.4.3         XVector_0.18.0        
##  [7] bindr_0.1              bitops_1.0-6           tools_3.4.3           
## [10] zlibbioc_1.24.0        digest_0.6.15          bit_1.1-12            
## [13] tibble_1.4.2           evaluate_0.10.1        lattice_0.20-34       
## [16] pkgconfig_2.0.1        rlang_0.2.0            Matrix_1.2-7.1        
## [19] yaml_2.1.17            xfun_0.1               bindrcpp_0.2          
## [22] GenomeInfoDbData_1.0.0 stringr_1.3.0          dplyr_0.7.4           
## [25] rprojroot_1.3-2        grid_3.4.3             glue_1.2.0            
## [28] R6_2.2.2               hash_2.2.6             rmarkdown_1.8         
## [31] bookdown_0.7           reshape2_1.4.3         magrittr_1.5          
## [34] backports_1.1.2        htmltools_0.3.6        assertthat_0.2.0      
## [37] stringi_1.1.6          RCurl_1.95-4.10
```

