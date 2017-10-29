---
knit: bookdown::preview_chapter
---

## Identifying confounding factors (Reads)



\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{12-confounders-reads_files/figure-latex/confound-pca-reads-1} 

}

\caption{PCA plot of the tung data}(\#fig:confound-pca-reads)
\end{figure}

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{12-confounders-reads_files/figure-latex/confound-find-pcs-total-features-reads-1} 

}

\caption{PC correlation with the number of detected genes}(\#fig:confound-find-pcs-total-features-reads)
\end{figure}

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{12-confounders-reads_files/figure-latex/confound-find-expl-vars-reads-1} 

}

\caption{Explanatory variables}(\#fig:confound-find-expl-vars-reads)
\end{figure}


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
##  [1] knitr_1.17                  scater_1.5.21              
##  [3] SingleCellExperiment_0.99.4 SummarizedExperiment_1.6.5 
##  [5] DelayedArray_0.2.7          matrixStats_0.52.2         
##  [7] GenomicRanges_1.28.6        GenomeInfoDb_1.12.3        
##  [9] IRanges_2.10.5              S4Vectors_0.14.7           
## [11] ggplot2_2.2.1               Biobase_2.36.2             
## [13] BiocGenerics_0.22.1        
## 
## loaded via a namespace (and not attached):
##  [1] viridis_0.4.0           edgeR_3.18.1           
##  [3] bit64_0.9-7             viridisLite_0.2.0      
##  [5] shiny_1.0.5             assertthat_0.2.0       
##  [7] blob_1.1.0              GenomeInfoDbData_0.99.0
##  [9] vipor_0.4.5             yaml_2.1.14            
## [11] RSQLite_2.0             backports_1.1.1        
## [13] lattice_0.20-34         glue_1.1.1             
## [15] limma_3.32.10           digest_0.6.12          
## [17] XVector_0.16.0          colorspace_1.3-2       
## [19] cowplot_0.8.0           htmltools_0.3.6        
## [21] httpuv_1.3.5            Matrix_1.2-7.1         
## [23] plyr_1.8.4              XML_3.98-1.9           
## [25] pkgconfig_2.0.1         biomaRt_2.32.1         
## [27] bookdown_0.5            zlibbioc_1.22.0        
## [29] xtable_1.8-2            scales_0.5.0           
## [31] tibble_1.3.4            lazyeval_0.2.0         
## [33] magrittr_1.5            mime_0.5               
## [35] memoise_1.1.0           evaluate_0.10.1        
## [37] beeswarm_0.2.3          shinydashboard_0.6.1   
## [39] tools_3.4.2             data.table_1.10.4-3    
## [41] stringr_1.2.0           munsell_0.4.3          
## [43] locfit_1.5-9.1          AnnotationDbi_1.38.2   
## [45] bindrcpp_0.2            compiler_3.4.2         
## [47] rlang_0.1.2             rhdf5_2.20.0           
## [49] grid_3.4.2              RCurl_1.95-4.8         
## [51] tximport_1.4.0          rjson_0.2.15           
## [53] labeling_0.3            bitops_1.0-6           
## [55] rmarkdown_1.6           gtable_0.2.0           
## [57] DBI_0.7                 reshape2_1.4.2         
## [59] R6_2.2.2                gridExtra_2.3          
## [61] dplyr_0.7.4             bit_1.1-12             
## [63] bindr_0.1               rprojroot_1.2          
## [65] stringi_1.1.5           ggbeeswarm_0.6.0       
## [67] Rcpp_0.12.13
```
