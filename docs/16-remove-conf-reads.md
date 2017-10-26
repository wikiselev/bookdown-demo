---
output: html_document
---

# Dealing with confounders (Reads)




```r
library(scRNA.seq.funcs)
library(RUVSeq)
library(scater)
library(SingleCellExperiment)
library(scran)
library(kBET)
library(sva) # Combat
library(edgeR)
set.seed(1234567)
options(stringsAsFactors = FALSE)
reads <- readRDS("tung/reads.rds")
reads.qc <- reads[rowData(reads)$use, colData(reads)$use]
endog_genes <- !rowData(reads.qc)$is_feature_control
erccs <- rowData(reads.qc)$is_feature_control

qclust <- quickCluster(reads.qc, min.size = 30)
reads.qc <- computeSumFactors(reads.qc, sizes = 15, clusters = qclust)
reads.qc <- normalize(reads.qc)
```

## Remove Unwanted Variation

### RUVg


```r
ruvg <- RUVg(counts(reads.qc), erccs, k = 1)
assay(reads.qc, "ruvg1") <- log2(
    t(t(ruvg$normalizedCounts) / colSums(ruvg$normalizedCounts) * 1e6) + 1
)
ruvg <- RUVg(counts(reads.qc), erccs, k = 10)
assay(reads.qc, "ruvg10") <- log2(
    t(t(ruvg$normalizedCounts) / colSums(ruvg$normalizedCounts) * 1e6) + 1
)
```

### RUVs


```r
scIdx <- matrix(-1, ncol = max(table(reads.qc$individual)), nrow = 3)
tmp <- which(reads.qc$individual == "NA19098")
scIdx[1, 1:length(tmp)] <- tmp
tmp <- which(reads.qc$individual == "NA19101")
scIdx[2, 1:length(tmp)] <- tmp
tmp <- which(reads.qc$individual == "NA19239")
scIdx[3, 1:length(tmp)] <- tmp
cIdx <- rownames(reads.qc)
ruvs <- RUVs(counts(reads.qc), cIdx, k = 1, scIdx = scIdx, isLog = FALSE)
assay(reads.qc, "ruvs1") <- log2(
    t(t(ruvs$normalizedCounts) / colSums(ruvs$normalizedCounts) * 1e6) + 1
)
ruvs <- RUVs(counts(reads.qc), cIdx, k = 10, scIdx = scIdx, isLog = FALSE)
assay(reads.qc, "ruvs10") <- log2(
    t(t(ruvs$normalizedCounts) / colSums(ruvs$normalizedCounts) * 1e6) + 1
)
```

### Combat


```r
combat_data <- logcounts(reads.qc)
mod_data <- as.data.frame(t(combat_data))
# Basic batch removal
mod0 = model.matrix(~ 1, data = mod_data) 
# Preserve biological variability
mod1 = model.matrix(~ reads.qc$individual, data = mod_data) 
# adjust for total genes detected
mod2 = model.matrix(~ reads.qc$total_features, data = mod_data)
assay(reads.qc, "combat") <- ComBat(
    dat = t(mod_data), 
    batch = factor(reads.qc$batch), 
    mod = mod0,
    par.prior = TRUE,
    prior.plots = FALSE
)
```

```
## Standardizing Data across genes
```

__Exercise 1__


```
## Standardizing Data across genes
```

### mnnCorrect 


```r
do_mnn <- function(data.qc) {
    batch1 <- data.qc[, data.qc$replicate == "r1"]
    batch2 <- data.qc[, data.qc$replicate == "r2"]
    batch3 <- data.qc[, data.qc$replicate == "r3"]
    
    if (ncol(batch2) > 0) {
        x = mnnCorrect(
          logcounts(batch1), logcounts(batch2), logcounts(batch3),  
          k = 20,
          sigma = 0.1,
          cos.norm = TRUE,
          svd.dim = 2
        )
        return(cbind(x$corrected[[1]], x$corrected[[2]], x$corrected[[3]]))
    } else {
        x = mnnCorrect(
          logcounts(batch1), logcounts(batch3),  
          k = 20,
          sigma = 0.1,
          cos.norm = TRUE,
          svd.dim = 2
        )
        return(cbind(x$corrected[[1]], x$corrected[[2]]))
    }
}

indi1 <- do_mnn(reads.qc[, reads.qc$individual == "NA19098"])
indi2 <- do_mnn(reads.qc[, reads.qc$individual == "NA19101"])
indi3 <- do_mnn(reads.qc[, reads.qc$individual == "NA19239"])

assay(reads.qc, "mnn") <- cbind(indi1, indi2, indi3);

# For a balanced design: 
#assay(reads.qc, "mnn") <- mnnCorrect(
#    list(B1 = logcounts(batch1), B2 = logcounts(batch2), B3 = logcounts(batch3)),  
#    k = 20,
#    sigma = 0.1,
#    cos.norm = TRUE,
#    svd.dim = 2
#)
```

### GLM


```r
glm_fun <- function(g, batch, indi) {
  model <- glm(g ~ batch + indi)
  model$coef[1] <- 0 # replace intercept with 0 to preserve reference batch.
  return(model$coef)
}
effects <- apply(
    logcounts(reads.qc), 
    1, 
    glm_fun, 
    batch = reads.qc$batch, 
    indi = reads.qc$individual
)
corrected <- logcounts(reads.qc) - t(effects[as.numeric(factor(reads.qc$batch)), ])
assay(reads.qc, "glm") <- corrected
```

__Exercise 2__



### Effectiveness 1


```r
for(n in assayNames(reads.qc)) {
    print(
        plotPCA(
            reads.qc[endog_genes, ],
            colour_by = "batch",
            size_by = "total_features",
            shape_by = "individual",
            exprs_values = n
        ) +
        ggtitle(n)
    )
}
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-10-1.png" width="672" style="display: block; margin: auto;" /><img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-10-2.png" width="672" style="display: block; margin: auto;" /><img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-10-3.png" width="672" style="display: block; margin: auto;" /><img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-10-4.png" width="672" style="display: block; margin: auto;" /><img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-10-5.png" width="672" style="display: block; margin: auto;" /><img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-10-6.png" width="672" style="display: block; margin: auto;" /><img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-10-7.png" width="672" style="display: block; margin: auto;" /><img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-10-8.png" width="672" style="display: block; margin: auto;" /><img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-10-9.png" width="672" style="display: block; margin: auto;" /><img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-10-10.png" width="672" style="display: block; margin: auto;" /><img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-10-11.png" width="672" style="display: block; margin: auto;" /><img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-10-12.png" width="672" style="display: block; margin: auto;" />

### Effectiveness 2


```r
res <- list()
for(n in assayNames(reads.qc)) {
	res[[n]] <- suppressWarnings(calc_cell_RLE(assay(reads.qc, n), erccs))
}
par(mar=c(6,4,1,1))
boxplot(res, las=2)
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-11-1.png" width="672" style="display: block; margin: auto;" />

### Effectiveness 3


```r
for(n in assayNames(reads.qc)) {
    print(
        plotQC(
            reads.qc[endog_genes, ],
            type = "expl",
            exprs_values = n,
            variables = c(
                "total_features",
                "total_counts",
                "batch",
                "individual",
                "pct_counts_ERCC",
                "pct_counts_MT"
            )
        ) +
        ggtitle(n)
    )
}
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-12-1.png" width="672" style="display: block; margin: auto;" /><img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-12-2.png" width="672" style="display: block; margin: auto;" /><img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-12-3.png" width="672" style="display: block; margin: auto;" /><img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-12-4.png" width="672" style="display: block; margin: auto;" /><img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-12-5.png" width="672" style="display: block; margin: auto;" />

```
## Warning: Removed 6 rows containing non-finite values (stat_density).
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-12-6.png" width="672" style="display: block; margin: auto;" />

```
## Warning: Removed 54 rows containing non-finite values (stat_density).
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-12-7.png" width="672" style="display: block; margin: auto;" /><img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-12-8.png" width="672" style="display: block; margin: auto;" /><img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-12-9.png" width="672" style="display: block; margin: auto;" /><img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-12-10.png" width="672" style="display: block; margin: auto;" /><img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-12-11.png" width="672" style="display: block; margin: auto;" />

```
## Warning in self$trans$transform(x): NaNs produced
```

```
## Warning: Transformation introduced infinite values in continuous x-axis
```

```
## Warning: Removed 129 rows containing non-finite values (stat_density).
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-12-12.png" width="672" style="display: block; margin: auto;" />

### Effectiveness 4


```r
compare_kBET_results <- function(sce){
    indiv <- unique(sce$individual)
    norms <- assayNames(sce) # Get all normalizations
    results <- list()
    for (i in indiv){ 
        for (j in norms){
            tmp <- kBET(
                df = t(assay(sce[,sce$individual== i], j)), 
                batch = sce$batch[sce$individual==i], 
                heuristic = TRUE, 
                verbose = FALSE, 
                addTest = FALSE, 
                plot = FALSE)
            results[[i]][[j]] <- tmp$summary$kBET.observed[1]
        }
    }
    return(as.data.frame(results))
}

eff_debatching <- compare_kBET_results(reads.qc)
```


```r
require("reshape2")
require("RColorBrewer")
# Plot results
dod <- melt(as.matrix(eff_debatching),  value.name = "kBET")
colnames(dod)[1:2] <- c("Normalisation", "Individual")

colorset <- c('gray', brewer.pal(n = 9, "RdYlBu"))

ggplot(dod, aes(Normalisation, Individual, fill=kBET)) +  
    geom_tile() +
    scale_fill_gradient2(
        na.value = "gray",
        low = colorset[2],
        mid=colorset[6],
        high = colorset[10],
        midpoint = 0.5, limit = c(0,1)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) + 
    theme(
        axis.text.x = element_text(
            angle = 45, 
            vjust = 1, 
            size = 12, 
            hjust = 1
        )
    ) + 
    ggtitle("Effect of batch regression methods per individual")
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-14-1.png" width="672" style="display: block; margin: auto;" />

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
##  [1] RColorBrewer_1.1-2          reshape2_1.4.2             
##  [3] sva_3.24.4                  genefilter_1.58.1          
##  [5] mgcv_1.8-22                 nlme_3.1-131               
##  [7] kBET_0.99.0                 scran_1.5.14               
##  [9] scater_1.5.21               SingleCellExperiment_0.99.4
## [11] ggplot2_2.2.1               RUVSeq_1.10.0              
## [13] edgeR_3.18.1                limma_3.32.10              
## [15] EDASeq_2.10.0               ShortRead_1.34.2           
## [17] GenomicAlignments_1.12.2    SummarizedExperiment_1.6.5 
## [19] DelayedArray_0.2.7          matrixStats_0.52.2         
## [21] Rsamtools_1.28.0            GenomicRanges_1.28.6       
## [23] GenomeInfoDb_1.12.3         Biostrings_2.44.2          
## [25] XVector_0.16.0              IRanges_2.10.5             
## [27] S4Vectors_0.14.7            BiocParallel_1.10.1        
## [29] Biobase_2.36.2              BiocGenerics_0.22.1        
## [31] scRNA.seq.funcs_0.1.0       knitr_1.17                 
## 
## loaded via a namespace (and not attached):
##  [1] Rtsne_0.13              ggbeeswarm_0.6.0       
##  [3] colorspace_1.3-2        rjson_0.2.15           
##  [5] hwriter_1.3.2           dynamicTreeCut_1.63-1  
##  [7] rprojroot_1.2           DT_0.2                 
##  [9] bit64_0.9-7             AnnotationDbi_1.38.2   
## [11] splines_3.4.2           R.methodsS3_1.7.1      
## [13] tximport_1.4.0          DESeq_1.28.0           
## [15] geneplotter_1.54.0      annotate_1.54.0        
## [17] cluster_2.0.6           R.oo_1.21.0            
## [19] shinydashboard_0.6.1    shiny_1.0.5            
## [21] compiler_3.4.2          backports_1.1.1        
## [23] assertthat_0.2.0        Matrix_1.2-11          
## [25] lazyeval_0.2.0          htmltools_0.3.6        
## [27] tools_3.4.2             igraph_1.1.2           
## [29] bindrcpp_0.2            gtable_0.2.0           
## [31] glue_1.1.1              GenomeInfoDbData_0.99.0
## [33] dplyr_0.7.4             Rcpp_0.12.13           
## [35] rtracklayer_1.36.6      stringr_1.2.0          
## [37] mime_0.5                hypergeo_1.2-13        
## [39] statmod_1.4.30          XML_3.98-1.9           
## [41] zlibbioc_1.22.0         MASS_7.3-47            
## [43] zoo_1.8-0               scales_0.5.0           
## [45] aroma.light_3.6.0       rhdf5_2.20.0           
## [47] yaml_2.1.14             memoise_1.1.0          
## [49] gridExtra_2.3           biomaRt_2.32.1         
## [51] latticeExtra_0.6-28     stringi_1.1.5          
## [53] RSQLite_2.0             orthopolynom_1.0-5     
## [55] GenomicFeatures_1.28.5  contfrac_1.1-11        
## [57] rlang_0.1.2             pkgconfig_2.0.1        
## [59] moments_0.14            bitops_1.0-6           
## [61] evaluate_0.10.1         lattice_0.20-35        
## [63] bindr_0.1               labeling_0.3           
## [65] htmlwidgets_0.9         cowplot_0.8.0          
## [67] bit_1.1-12              deSolve_1.20           
## [69] plyr_1.8.4              magrittr_1.5           
## [71] bookdown_0.5            R6_2.2.2               
## [73] DBI_0.7                 survival_2.41-3        
## [75] RCurl_1.95-4.8          tibble_1.3.4           
## [77] rmarkdown_1.6           viridis_0.4.0          
## [79] locfit_1.5-9.1          grid_3.4.2             
## [81] data.table_1.10.4-2     FNN_1.1                
## [83] blob_1.1.0              digest_0.6.12          
## [85] xtable_1.8-2            httpuv_1.3.5           
## [87] elliptic_1.3-7          R.utils_2.5.0          
## [89] munsell_0.4.3           beeswarm_0.2.3         
## [91] viridisLite_0.2.0       vipor_0.4.5
```
