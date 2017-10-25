---
output: html_document
---

## Normalization practice (Reads)



<div class="figure" style="text-align: center">
<img src="14-exprs-norm-reads_files/figure-html/norm-pca-raw-reads-1.png" alt="PCA plot of the tung data" width="90%" />
<p class="caption">(\#fig:norm-pca-raw-reads)PCA plot of the tung data</p>
</div>

<div class="figure" style="text-align: center">
<img src="14-exprs-norm-reads_files/figure-html/norm-pca-cpm-reads-1.png" alt="PCA plot of the tung data after CPM normalisation" width="90%" />
<p class="caption">(\#fig:norm-pca-cpm-reads)PCA plot of the tung data after CPM normalisation</p>
</div>
<div class="figure" style="text-align: center">
<img src="14-exprs-norm-reads_files/figure-html/norm-ours-rle-cpm-reads-1.png" alt="Cell-wise RLE of the tung data" width="90%" />
<p class="caption">(\#fig:norm-ours-rle-cpm-reads)Cell-wise RLE of the tung data</p>
</div>


```
## Warning in normalizeSCE(object, exprs_values = exprs_values,
## return_norm_as_exprs = return_norm_as_exprs): spike-in transcripts in
## 'ERCC' should have their own size factors
```

<div class="figure" style="text-align: center">
<img src="14-exprs-norm-reads_files/figure-html/norm-pca-tmm-reads-1.png" alt="PCA plot of the tung data after TMM normalisation" width="90%" />
<p class="caption">(\#fig:norm-pca-tmm-reads)PCA plot of the tung data after TMM normalisation</p>
</div>
<div class="figure" style="text-align: center">
<img src="14-exprs-norm-reads_files/figure-html/norm-ours-rle-tmm-reads-1.png" alt="Cell-wise RLE of the tung data" width="90%" />
<p class="caption">(\#fig:norm-ours-rle-tmm-reads)Cell-wise RLE of the tung data</p>
</div>


```
## Warning in .local(object, ...): spike-in transcripts in 'ERCC' should have
## their own size factors
```

<div class="figure" style="text-align: center">
<img src="14-exprs-norm-reads_files/figure-html/norm-pca-lsf-umi-1.png" alt="PCA plot of the tung data after LSF normalisation" width="90%" />
<p class="caption">(\#fig:norm-pca-lsf-umi)PCA plot of the tung data after LSF normalisation</p>
</div>

<div class="figure" style="text-align: center">
<img src="14-exprs-norm-reads_files/figure-html/norm-ours-rle-scran-reads-1.png" alt="Cell-wise RLE of the tung data" width="90%" />
<p class="caption">(\#fig:norm-ours-rle-scran-reads)Cell-wise RLE of the tung data</p>
</div>


```
## Warning in normalizeSCE(object, exprs_values = exprs_values,
## return_norm_as_exprs = return_norm_as_exprs): spike-in transcripts in
## 'ERCC' should have their own size factors
```

<div class="figure" style="text-align: center">
<img src="14-exprs-norm-reads_files/figure-html/norm-pca-rle-reads-1.png" alt="PCA plot of the tung data after RLE normalisation" width="90%" />
<p class="caption">(\#fig:norm-pca-rle-reads)PCA plot of the tung data after RLE normalisation</p>
</div>

<div class="figure" style="text-align: center">
<img src="14-exprs-norm-reads_files/figure-html/norm-ours-rle-rle-reads-1.png" alt="Cell-wise RLE of the tung data" width="90%" />
<p class="caption">(\#fig:norm-ours-rle-rle-reads)Cell-wise RLE of the tung data</p>
</div>


```
## Warning in normalizeSCE(object, exprs_values = exprs_values,
## return_norm_as_exprs = return_norm_as_exprs): spike-in transcripts in
## 'ERCC' should have their own size factors
```

<div class="figure" style="text-align: center">
<img src="14-exprs-norm-reads_files/figure-html/norm-pca-uq-reads-1.png" alt="PCA plot of the tung data after UQ normalisation" width="90%" />
<p class="caption">(\#fig:norm-pca-uq-reads)PCA plot of the tung data after UQ normalisation</p>
</div>
<div class="figure" style="text-align: center">
<img src="14-exprs-norm-reads_files/figure-html/norm-ours-rle-uq-reads-1.png" alt="Cell-wise RLE of the tung data" width="90%" />
<p class="caption">(\#fig:norm-ours-rle-uq-reads)Cell-wise RLE of the tung data</p>
</div>

<div class="figure" style="text-align: center">
<img src="14-exprs-norm-reads_files/figure-html/norm-pca-downsample-reads-1.png" alt="PCA plot of the tung data after downsampling" width="90%" />
<p class="caption">(\#fig:norm-pca-downsample-reads)PCA plot of the tung data after downsampling</p>
</div>
<div class="figure" style="text-align: center">
<img src="14-exprs-norm-reads_files/figure-html/norm-ours-rle-downsample-reads-1.png" alt="Cell-wise RLE of the tung data" width="90%" />
<p class="caption">(\#fig:norm-ours-rle-downsample-reads)Cell-wise RLE of the tung data</p>
</div>









<div class="figure" style="text-align: center">
<img src="14-exprs-norm-reads_files/figure-html/norm-pca-tpm-reads-1.png" alt="PCA plot of the tung data after TPM normalisation" width="90%" />
<p class="caption">(\#fig:norm-pca-tpm-reads)PCA plot of the tung data after TPM normalisation</p>
</div>



<div class="figure" style="text-align: center">
<img src="14-exprs-norm-reads_files/figure-html/norm-pca-fpkm-reads-1.png" alt="PCA plot of the tung data after FPKM normalisation" width="90%" />
<p class="caption">(\#fig:norm-pca-fpkm-reads)PCA plot of the tung data after FPKM normalisation</p>
</div>

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
## [1] stats4    parallel  methods   stats     graphics  grDevices utils    
## [8] datasets  base     
## 
## other attached packages:
##  [1] knitr_1.17                  scran_1.5.14               
##  [3] BiocParallel_1.10.1         scater_1.5.20              
##  [5] SingleCellExperiment_0.99.4 SummarizedExperiment_1.6.5 
##  [7] DelayedArray_0.2.7          matrixStats_0.52.2         
##  [9] GenomicRanges_1.28.6        GenomeInfoDb_1.12.3        
## [11] IRanges_2.10.5              S4Vectors_0.14.7           
## [13] ggplot2_2.2.1               Biobase_2.36.2             
## [15] BiocGenerics_0.22.1         scRNA.seq.funcs_0.1.0      
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6            bit64_0.9-7            
##  [3] rprojroot_1.2           dynamicTreeCut_1.63-1  
##  [5] tools_3.4.2             backports_1.1.1        
##  [7] DT_0.2                  R6_2.2.2               
##  [9] hypergeo_1.2-13         vipor_0.4.5            
## [11] DBI_0.7                 lazyeval_0.2.0         
## [13] colorspace_1.3-2        gridExtra_2.3          
## [15] moments_0.14            bit_1.1-12             
## [17] compiler_3.4.2          orthopolynom_1.0-5     
## [19] labeling_0.3            bookdown_0.5           
## [21] scales_0.5.0            stringr_1.2.0          
## [23] digest_0.6.12           rmarkdown_1.6          
## [25] XVector_0.16.0          pkgconfig_2.0.1        
## [27] htmltools_0.3.6         highr_0.6              
## [29] limma_3.32.10           htmlwidgets_0.9        
## [31] rlang_0.1.2             RSQLite_2.0            
## [33] FNN_1.1                 shiny_1.0.5            
## [35] bindr_0.1               zoo_1.8-0              
## [37] dplyr_0.7.4             RCurl_1.95-4.8         
## [39] magrittr_1.5            GenomeInfoDbData_0.99.0
## [41] Matrix_1.2-11           Rcpp_0.12.13           
## [43] ggbeeswarm_0.6.0        munsell_0.4.3          
## [45] viridis_0.4.0           stringi_1.1.5          
## [47] yaml_2.1.14             edgeR_3.18.1           
## [49] MASS_7.3-47             zlibbioc_1.22.0        
## [51] rhdf5_2.20.0            Rtsne_0.13             
## [53] plyr_1.8.4              grid_3.4.2             
## [55] blob_1.1.0              shinydashboard_0.6.1   
## [57] contfrac_1.1-11         lattice_0.20-35        
## [59] cowplot_0.8.0           splines_3.4.2          
## [61] locfit_1.5-9.1          igraph_1.1.2           
## [63] rjson_0.2.15            reshape2_1.4.2         
## [65] biomaRt_2.32.1          XML_3.98-1.9           
## [67] glue_1.1.1              evaluate_0.10.1        
## [69] data.table_1.10.4-2     deSolve_1.20           
## [71] httpuv_1.3.5            gtable_0.2.0           
## [73] assertthat_0.2.0        mime_0.5               
## [75] xtable_1.8-2            viridisLite_0.2.0      
## [77] tibble_1.3.4            elliptic_1.3-7         
## [79] AnnotationDbi_1.38.2    beeswarm_0.2.3         
## [81] memoise_1.1.0           tximport_1.4.0         
## [83] bindrcpp_0.2            statmod_1.4.30
```

