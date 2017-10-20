---
output: html_document
---




```r
set.seed(1234567)
```

## SEURAT

Let's load the data and look at it:

```r
pollen <- readRDS("pollen/pollen.rds")
```

Here we follow an [example](http://satijalab.org/seurat/get_started.html) created by the authors of `SEURAT` (8,500 Pancreas cells). We mostly use default values in various function calls, for more details please consult the documentation and the authors:





Initialize the Seurat object with the raw (non-normalized data).  Keep all genes expressed in >= 3 cells. Keep all cells with at least 200 detected genes:

```r
library(Seurat)
seuset <- CreateSeuratObject(
    raw.data = logcounts(pollen), 
    min.cells = 3, 
    min.genes = 200
)

VlnPlot(object = seuset, features.plot = c("nGene", "nUMI"), nCol = 3)

GenePlot(object = seuset, gene1 = "nUMI", gene2 = "nGene")

seuset <- FindVariableGenes(
    object = seuset, 
    mean.function = ExpMean, 
    dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, 
    x.high.cutoff = 3, 
    y.cutoff = 0.5
)

length(x = seuset@var.genes)

seuset <- ScaleData(object = seuset, vars.to.regress = c("nUMI"))

seuset <- RunPCA(
    object = seuset, 
    pc.genes = seuset@var.genes, 
    do.print = TRUE, 
    pcs.print = 1:5, 
    genes.print = 5
)

VizPCA(object = seuset, pcs.use = 1:2)

PCAPlot(object = seuset, dim.1 = 1, dim.2 = 2)

seuset <- ProjectPCA(object = seuset, do.print = FALSE)

PCHeatmap(
    object = seuset, 
    pc.use = 1:6, 
    cells.use = 500, 
    do.balanced = TRUE, 
    label.columns = FALSE,
    use.full = FALSE
)

seuset <- JackStraw(
    object = seuset, 
    num.replicate = 100, 
    do.print = FALSE
)

JackStrawPlot(object = seuset, PCs = 1:6)

PCElbowPlot(object = seuset)

seuset <- FindClusters(
    object = seuset, 
    reduction.type = "pca", 
    dims.use = 1:5, 
    resolution = 0.6, 
    print.output = 0, 
    save.SNN = TRUE
)

table(seuset@ident)

mclust::adjustedRandIndex(colData(pollen)$cell_type1, seuset@ident)

# this did not work!
# Error in La.svd(x, nu, nv) : LAPACK routines cannot be loaded
# seuset <- RunTSNE(
#     object = seuset, 
#     dims.use = 1:5, 
#     do.fast = TRUE
# )
```

__Exercise 12__: Compare the results between `SC3` and `SEURAT`.

__Our solution__:

```r
pData(pollen)$SEURAT <- as.character(pollen_seurat@ident)
sc3_plot_expression(pollen, k = 11, show_pdata = "SEURAT")
```


Seurat can also find marker genes, e.g. marker genes for cluster 2:

```r
markers <- FindMarkers(pollen_seurat, 2)
FeaturePlot(pollen_seurat, 
            head(rownames(markers)), 
            cols.use = c("lightgrey", "blue"), 
            nCol = 3)
```

__Exercise 13__: Compare marker genes provided by `SEURAT` and `SC3`.

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
## [1] stats     graphics  grDevices utils     datasets  base     
## 
## other attached packages:
## [1] knitr_1.17
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.13               XVector_0.16.0            
##  [3] magrittr_1.5               GenomicRanges_1.28.6      
##  [5] BiocGenerics_0.22.1        zlibbioc_1.22.0           
##  [7] IRanges_2.10.5             lattice_0.20-35           
##  [9] stringr_1.2.0              GenomeInfoDb_1.12.3       
## [11] tools_3.4.2                grid_3.4.2                
## [13] SummarizedExperiment_1.6.5 parallel_3.4.2            
## [15] Biobase_2.36.2             matrixStats_0.52.2        
## [17] htmltools_0.3.6            yaml_2.1.14               
## [19] rprojroot_1.2              digest_0.6.12             
## [21] bookdown_0.5               Matrix_1.2-11             
## [23] GenomeInfoDbData_0.99.0    S4Vectors_0.14.7          
## [25] bitops_1.0-6               RCurl_1.95-4.8            
## [27] evaluate_0.10.1            rmarkdown_1.6             
## [29] DelayedArray_0.2.7         stringi_1.1.5             
## [31] compiler_3.4.2             methods_3.4.2             
## [33] backports_1.1.1            stats4_3.4.2
```
