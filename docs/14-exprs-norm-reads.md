---
output: html_document
---

## Normalization practice (Reads)



\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{14-exprs-norm-reads_files/figure-latex/norm-pca-raw-reads-1} 

}

\caption{PCA plot of the tung data}(\#fig:norm-pca-raw-reads)
\end{figure}

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{14-exprs-norm-reads_files/figure-latex/norm-pca-cpm-reads-1} 

}

\caption{PCA plot of the tung data after CPM normalisation}(\#fig:norm-pca-cpm-reads)
\end{figure}
\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{14-exprs-norm-reads_files/figure-latex/norm-ours-rle-cpm-reads-1} 

}

\caption{Cell-wise RLE of the tung data}(\#fig:norm-ours-rle-cpm-reads)
\end{figure}


```
## Warning in normalizeSCE(object, exprs_values = exprs_values, return_log
## = return_log, : spike-in transcripts in 'ERCC' should have their own size
## factors
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{14-exprs-norm-reads_files/figure-latex/norm-pca-rle-reads-1} 

}

\caption{PCA plot of the tung data after RLE normalisation}(\#fig:norm-pca-rle-reads)
\end{figure}

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{14-exprs-norm-reads_files/figure-latex/norm-ours-rle-rle-reads-1} 

}

\caption{Cell-wise RLE of the tung data}(\#fig:norm-ours-rle-rle-reads)
\end{figure}


```
## Warning in normalizeSCE(object, exprs_values = exprs_values, return_log
## = return_log, : spike-in transcripts in 'ERCC' should have their own size
## factors
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{14-exprs-norm-reads_files/figure-latex/norm-pca-uq-reads-1} 

}

\caption{PCA plot of the tung data after UQ normalisation}(\#fig:norm-pca-uq-reads)
\end{figure}
\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{14-exprs-norm-reads_files/figure-latex/norm-ours-rle-uq-reads-1} 

}

\caption{Cell-wise RLE of the tung data}(\#fig:norm-ours-rle-uq-reads)
\end{figure}


```
## Warning in normalizeSCE(object, exprs_values = exprs_values, return_log
## = return_log, : spike-in transcripts in 'ERCC' should have their own size
## factors
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{14-exprs-norm-reads_files/figure-latex/norm-pca-tmm-reads-1} 

}

\caption{PCA plot of the tung data after TMM normalisation}(\#fig:norm-pca-tmm-reads)
\end{figure}
\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{14-exprs-norm-reads_files/figure-latex/norm-ours-rle-tmm-reads-1} 

}

\caption{Cell-wise RLE of the tung data}(\#fig:norm-ours-rle-tmm-reads)
\end{figure}


```
## Warning in .local(object, ...): spike-in transcripts in 'ERCC' should have
## their own size factors
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{14-exprs-norm-reads_files/figure-latex/norm-pca-lsf-umi-1} 

}

\caption{PCA plot of the tung data after LSF normalisation}(\#fig:norm-pca-lsf-umi)
\end{figure}

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{14-exprs-norm-reads_files/figure-latex/norm-ours-rle-scran-reads-1} 

}

\caption{Cell-wise RLE of the tung data}(\#fig:norm-ours-rle-scran-reads)
\end{figure}

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{14-exprs-norm-reads_files/figure-latex/norm-pca-downsample-reads-1} 

}

\caption{PCA plot of the tung data after downsampling}(\#fig:norm-pca-downsample-reads)
\end{figure}
\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{14-exprs-norm-reads_files/figure-latex/norm-ours-rle-downsample-reads-1} 

}

\caption{Cell-wise RLE of the tung data}(\#fig:norm-ours-rle-downsample-reads)
\end{figure}









\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{14-exprs-norm-reads_files/figure-latex/norm-pca-tpm-reads-1} 

}

\caption{PCA plot of the tung data after TPM normalisation}(\#fig:norm-pca-tpm-reads)
\end{figure}



\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{14-exprs-norm-reads_files/figure-latex/norm-pca-fpkm-reads-1} 

}

\caption{PCA plot of the tung data after FPKM normalisation}(\#fig:norm-pca-fpkm-reads)
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
##  [1] knitr_1.17                  scran_1.5.15               
##  [3] BiocParallel_1.10.1         scater_1.5.21              
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
## [27] htmltools_0.3.6         limma_3.32.10          
## [29] htmlwidgets_0.9         rlang_0.1.2            
## [31] RSQLite_2.0             FNN_1.1                
## [33] shiny_1.0.5             bindr_0.1              
## [35] zoo_1.8-0               dplyr_0.7.4            
## [37] RCurl_1.95-4.8          magrittr_1.5           
## [39] GenomeInfoDbData_0.99.0 Matrix_1.2-7.1         
## [41] Rcpp_0.12.13            ggbeeswarm_0.6.0       
## [43] munsell_0.4.3           viridis_0.4.0          
## [45] stringi_1.1.5           yaml_2.1.14            
## [47] edgeR_3.18.1            MASS_7.3-45            
## [49] zlibbioc_1.22.0         rhdf5_2.20.0           
## [51] Rtsne_0.13              plyr_1.8.4             
## [53] grid_3.4.2              blob_1.1.0             
## [55] shinydashboard_0.6.1    contfrac_1.1-11        
## [57] lattice_0.20-34         cowplot_0.8.0          
## [59] splines_3.4.2           locfit_1.5-9.1         
## [61] igraph_1.1.2            rjson_0.2.15           
## [63] reshape2_1.4.2          biomaRt_2.32.1         
## [65] XML_3.98-1.9            glue_1.1.1             
## [67] evaluate_0.10.1         data.table_1.10.4-3    
## [69] deSolve_1.20            httpuv_1.3.5           
## [71] gtable_0.2.0            assertthat_0.2.0       
## [73] mime_0.5                xtable_1.8-2           
## [75] viridisLite_0.2.0       tibble_1.3.4           
## [77] elliptic_1.3-7          AnnotationDbi_1.38.2   
## [79] beeswarm_0.2.3          memoise_1.1.0          
## [81] tximport_1.4.0          bindrcpp_0.2           
## [83] statmod_1.4.30
```

