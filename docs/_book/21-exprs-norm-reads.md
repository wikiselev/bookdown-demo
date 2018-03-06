---
output: html_document
---

## Normalization practice (Reads)



\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{21-exprs-norm-reads_files/figure-latex/norm-pca-raw-reads-1} 

}

\caption{PCA plot of the tung data}(\#fig:norm-pca-raw-reads)
\end{figure}

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{21-exprs-norm-reads_files/figure-latex/norm-pca-cpm-reads-1} 

}

\caption{PCA plot of the tung data after CPM normalisation}(\#fig:norm-pca-cpm-reads)
\end{figure}
\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{21-exprs-norm-reads_files/figure-latex/norm-ours-rle-cpm-reads-1} 

}

\caption{Cell-wise RLE of the tung data}(\#fig:norm-ours-rle-cpm-reads)
\end{figure}


```
## Warning in normalizeSCE(object, exprs_values = exprs_values, return_log
## = return_log, : spike-in transcripts in 'ERCC' should have their own size
## factors
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{21-exprs-norm-reads_files/figure-latex/norm-pca-rle-reads-1} 

}

\caption{PCA plot of the tung data after RLE normalisation}(\#fig:norm-pca-rle-reads)
\end{figure}

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{21-exprs-norm-reads_files/figure-latex/norm-ours-rle-rle-reads-1} 

}

\caption{Cell-wise RLE of the tung data}(\#fig:norm-ours-rle-rle-reads)
\end{figure}


```
## Warning in normalizeSCE(object, exprs_values = exprs_values, return_log
## = return_log, : spike-in transcripts in 'ERCC' should have their own size
## factors
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{21-exprs-norm-reads_files/figure-latex/norm-pca-uq-reads-1} 

}

\caption{PCA plot of the tung data after UQ normalisation}(\#fig:norm-pca-uq-reads)
\end{figure}
\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{21-exprs-norm-reads_files/figure-latex/norm-ours-rle-uq-reads-1} 

}

\caption{Cell-wise RLE of the tung data}(\#fig:norm-ours-rle-uq-reads)
\end{figure}


```
## Warning in normalizeSCE(object, exprs_values = exprs_values, return_log
## = return_log, : spike-in transcripts in 'ERCC' should have their own size
## factors
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{21-exprs-norm-reads_files/figure-latex/norm-pca-tmm-reads-1} 

}

\caption{PCA plot of the tung data after TMM normalisation}(\#fig:norm-pca-tmm-reads)
\end{figure}
\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{21-exprs-norm-reads_files/figure-latex/norm-ours-rle-tmm-reads-1} 

}

\caption{Cell-wise RLE of the tung data}(\#fig:norm-ours-rle-tmm-reads)
\end{figure}


```
## Warning in .local(object, ...): spike-in transcripts in 'ERCC' should have
## their own size factors
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{21-exprs-norm-reads_files/figure-latex/norm-pca-lsf-umi-1} 

}

\caption{PCA plot of the tung data after LSF normalisation}(\#fig:norm-pca-lsf-umi)
\end{figure}

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{21-exprs-norm-reads_files/figure-latex/norm-ours-rle-scran-reads-1} 

}

\caption{Cell-wise RLE of the tung data}(\#fig:norm-ours-rle-scran-reads)
\end{figure}

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{21-exprs-norm-reads_files/figure-latex/norm-pca-downsample-reads-1} 

}

\caption{PCA plot of the tung data after downsampling}(\#fig:norm-pca-downsample-reads)
\end{figure}
\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{21-exprs-norm-reads_files/figure-latex/norm-ours-rle-downsample-reads-1} 

}

\caption{Cell-wise RLE of the tung data}(\#fig:norm-ours-rle-downsample-reads)
\end{figure}
















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
##  [1] knitr_1.20                 scran_1.6.8               
##  [3] BiocParallel_1.12.0        scater_1.6.3              
##  [5] SingleCellExperiment_1.0.0 SummarizedExperiment_1.8.1
##  [7] DelayedArray_0.4.1         matrixStats_0.53.1        
##  [9] GenomicRanges_1.30.3       GenomeInfoDb_1.14.0       
## [11] IRanges_2.12.0             S4Vectors_0.16.0          
## [13] ggplot2_2.2.1              Biobase_2.38.0            
## [15] BiocGenerics_0.24.0        scRNA.seq.funcs_0.1.0     
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6           bit64_0.9-7            progress_1.1.2        
##  [4] httr_1.3.1             rprojroot_1.3-2        dynamicTreeCut_1.63-1 
##  [7] tools_3.4.3            backports_1.1.2        DT_0.4                
## [10] R6_2.2.2               hypergeo_1.2-13        vipor_0.4.5           
## [13] DBI_0.7                lazyeval_0.2.1         colorspace_1.3-2      
## [16] gridExtra_2.3          prettyunits_1.0.2      moments_0.14          
## [19] bit_1.1-12             compiler_3.4.3         orthopolynom_1.0-5    
## [22] labeling_0.3           bookdown_0.7           scales_0.5.0          
## [25] stringr_1.3.0          digest_0.6.15          rmarkdown_1.8         
## [28] XVector_0.18.0         pkgconfig_2.0.1        htmltools_0.3.6       
## [31] limma_3.34.9           htmlwidgets_1.0        rlang_0.2.0           
## [34] RSQLite_2.0            FNN_1.1                shiny_1.0.5           
## [37] bindr_0.1              zoo_1.8-1              dplyr_0.7.4           
## [40] RCurl_1.95-4.10        magrittr_1.5           GenomeInfoDbData_1.0.0
## [43] Matrix_1.2-7.1         Rcpp_0.12.15           ggbeeswarm_0.6.0      
## [46] munsell_0.4.3          viridis_0.5.0          stringi_1.1.6         
## [49] yaml_2.1.17            edgeR_3.20.9           MASS_7.3-45           
## [52] zlibbioc_1.24.0        rhdf5_2.22.0           Rtsne_0.13            
## [55] plyr_1.8.4             grid_3.4.3             blob_1.1.0            
## [58] shinydashboard_0.6.1   contfrac_1.1-11        lattice_0.20-34       
## [61] cowplot_0.9.2          locfit_1.5-9.1         pillar_1.2.1          
## [64] igraph_1.1.2           rjson_0.2.15           reshape2_1.4.3        
## [67] biomaRt_2.34.2         XML_3.98-1.10          glue_1.2.0            
## [70] evaluate_0.10.1        data.table_1.10.4-3    deSolve_1.20          
## [73] httpuv_1.3.6.1         gtable_0.2.0           assertthat_0.2.0      
## [76] xfun_0.1               mime_0.5               xtable_1.8-2          
## [79] viridisLite_0.3.0      tibble_1.4.2           elliptic_1.3-7        
## [82] AnnotationDbi_1.40.0   beeswarm_0.2.3         memoise_1.1.0         
## [85] tximport_1.6.0         bindrcpp_0.2           statmod_1.4.30
```

