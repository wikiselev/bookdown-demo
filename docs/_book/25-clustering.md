---
output: html_document
---

## Clustering example {#clust-methods}




```r
library(pcaMethods)
library(pcaReduce)
library(SC3)
library(scater)
library(SingleCellExperiment)
library(pheatmap)
library(mclust)
set.seed(1234567)
```

To illustrate clustering of scRNA-seq data, we consider the `Deng` dataset of cells from developing mouse embryo [@Deng2014-mx]. We have preprocessed the dataset and created a `SingleCellExperiment` object in advance. We have also annotated the cells with the cell types identified in the original publication (it is the `cell_type2` column in the `colData` slot).

### Deng dataset

Let's load the data and look at it:

```r
deng <- readRDS("deng/deng-reads.rds")
deng
```

```
## class: SingleCellExperiment 
## dim: 22431 268 
## metadata(0):
## assays(2): counts logcounts
## rownames(22431): Hvcn1 Gbp7 ... Sox5 Alg11
## rowData names(10): feature_symbol is_feature_control ...
##   total_counts log10_total_counts
## colnames(268): 16cell 16cell.1 ... zy.2 zy.3
## colData names(30): cell_type2 cell_type1 ... pct_counts_ERCC
##   is_cell_control
## reducedDimNames(0):
## spikeNames(1): ERCC
```

Let's look at the cell type annotation:

```r
table(colData(deng)$cell_type2)
```

```
## 
##     16cell      4cell      8cell early2cell earlyblast  late2cell 
##         50         14         37          8         43         10 
##  lateblast   mid2cell   midblast         zy 
##         30         12         60          4
```

A simple PCA analysis already separates some strong cell types and provides some insights in the data structure:

```r
plotPCA(deng, colour_by = "cell_type2")
```

<img src="25-clustering_files/figure-html/unnamed-chunk-5-1.png" width="672" style="display: block; margin: auto;" />
As you can see, the early cell types separate quite well, but the three blastocyst timepoints are more difficult to distinguish.

### SC3

Let's run `SC3` clustering on the Deng data. The advantage of the `SC3` is that it can directly ingest a `SingleCellExperiment` object.

Now let's image we do not know the number of clusters _k_ (cell types). `SC3` can estimate a number of clusters for you:

```r
deng <- sc3_estimate_k(deng)
```

```
## Estimating k...
```

```r
metadata(deng)$sc3$k_estimation
```

```
## [1] 6
```

Interestingly, the number of cell types predicted by `SC3` is smaller than in the original data annotation. However, early, mid and late stages of different cell types together, we will have exactly 6 cell types. We store the merged cell types in `cell_type1` column of the `colData` slot:

```r
plotPCA(deng, colour_by = "cell_type1")
```

<img src="25-clustering_files/figure-html/unnamed-chunk-7-1.png" width="672" style="display: block; margin: auto;" />

Now we are ready to run `SC3` (we also ask it to calculate biological properties of the clusters): 

```r
deng <- sc3(deng, ks = 10, biology = TRUE)
```

```
## Setting SC3 parameters...
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

```
## Calculating biology...
```

`SC3` result consists of several different outputs (please look in [@Kiselev2016-bq] and [SC3 vignette](http://bioconductor.org/packages/release/bioc/vignettes/SC3/inst/doc/my-vignette.html) for more details). Here we show some of them:

Consensus matrix:

```r
sc3_plot_consensus(deng, k = 10, show_pdata = "cell_type2")
```

<img src="25-clustering_files/figure-html/unnamed-chunk-9-1.png" width="672" style="display: block; margin: auto;" />

Silhouette plot:

```r
sc3_plot_silhouette(deng, k = 10)
```

<img src="25-clustering_files/figure-html/unnamed-chunk-10-1.png" width="672" style="display: block; margin: auto;" />

Heatmap of the expression matrix:

```r
sc3_plot_expression(deng, k = 10, show_pdata = "cell_type2")
```

<img src="25-clustering_files/figure-html/unnamed-chunk-11-1.png" width="672" style="display: block; margin: auto;" />

Identified marker genes:

```r
sc3_plot_markers(deng, k = 10, show_pdata = "cell_type2")
```

<img src="25-clustering_files/figure-html/unnamed-chunk-12-1.png" width="672" style="display: block; margin: auto;" />

PCA plot with highlighted `SC3` clusters:

```r
plotPCA(deng, colour_by = "sc3_10_clusters")
```

<img src="25-clustering_files/figure-html/unnamed-chunk-13-1.png" width="672" style="display: block; margin: auto;" />

Compare the results of `SC3` clustering with the original publication cell type labels:

```r
adjustedRandIndex(colData(deng)$cell_type2, colData(deng)$sc3_10_clusters)
```

```
## [1] 0.7705208
```

__Note__ `SC3` can also be run in an interactive `Shiny` session:

```r
sc3_interactive(deng)
```

This command will open `SC3` in a web browser.

__Note__ Due to direct calculation of distances `SC3` becomes very slow when the number of cells is $>5000$. For large datasets containing up to $10^5$ cells we recomment using `Seurat` (see chapter \@ref(seurat-chapter)).

* __Exercise 1__: Run `SC3` for $k$ from 8 to 12 and explore different clustering solutions in your web browser.

* __Exercise 2__: Which clusters are the most stable when $k$ is changed from 8 to 12? (Look at the "Stability" tab)

* __Exercise 3__: Check out differentially expressed genes and marker genes for the obtained clusterings. Please use $k=10$.

* __Exercise 4__: Change the marker genes threshold (the default is 0.85). Does __SC3__ find more marker genes?

### pcaReduce

`pcaReduce` operates directly on the expression matrix. It is recommended to use a gene filter and log transformation before running `pcaReduce`. We will use the default `SC3` gene filter (note that the `exprs` slot of a `scater` object is log-transformed by default).


```r
# use the same gene filter as in SC3
input <- logcounts(deng[rowData(deng)$sc3_gene_filter, ])
```

There are several parameters used by `pcaReduce`:
* `nbt` defines a number of `pcaReduce` runs (it is stochastic and may have different solutions after different runs)
* `q` defines number of dimensions to start clustering with. The output will contain partitions for all $k$ from 2 to q+1.
* `method` defines a method used for clustering. `S` - to perform sampling based merging, `M` - to perform merging based on largest probability.

We will run `pcaReduce` 1 time:

```r
# run pcaReduce 1 time creating hierarchies from 1 to 30 clusters
pca.red <- PCAreduce(t(input), nbt = 1, q = 30, method = 'S')[[1]]
```


```r
colData(deng)$pcaReduce <- as.character(pca.red[,32 - 10])
plotPCA(deng, colour_by = "pcaReduce")
```

<img src="25-clustering_files/figure-html/unnamed-chunk-18-1.png" width="672" style="display: block; margin: auto;" />

__Exercise 5__: Run pcaReduce for $k=2$ and plot a similar PCA plot. Does it look good?

__Hint__: When running pcaReduce for different $k$s you do not need to rerun PCAreduce function, just use already calculated `pca.red` object.

__Our solution__:
<div class="figure" style="text-align: center">
<img src="25-clustering_files/figure-html/clust-pca-reduce2-1.png" alt="Clustering solutions of pcaReduce method for $k=2$." width="672" />
<p class="caption">(\#fig:clust-pca-reduce2)Clustering solutions of pcaReduce method for $k=2$.</p>
</div>

__Exercise 6__: Compare the results between `pcaReduce` and the original publication cell types for $k=10$.

__Our solution__:

```
## [1] 0.4216031
```

### tSNE + kmeans

[tSNE](https://lvdmaaten.github.io/tsne/) plots that we saw before (\@ref(visual-tsne)) when used the __scater__ package are made by using the [Rtsne](https://cran.r-project.org/web/packages/Rtsne/index.html) and [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) packages. Here we will do the same:

```r
deng <- plotTSNE(deng, rand_seed = 1, return_SCE = TRUE)
```

<div class="figure" style="text-align: center">
<img src="25-clustering_files/figure-html/clust-tsne-1.png" alt="tSNE map of the patient data" width="672" />
<p class="caption">(\#fig:clust-tsne)tSNE map of the patient data</p>
</div>

Note that all points on the plot above are black. This is different from what we saw before, when the cells were coloured based on the annotation. Here we do not have any annotation and all cells come from the same batch, therefore all dots are black.

Now we are going to apply _k_-means clustering algorithm to the cloud of points on the tSNE map. How many groups do you see in the cloud?

We will start with $k=8$:

```r
colData(deng)$tSNE_kmeans <- as.character(kmeans(deng@reducedDims$TSNE, centers = 8)$clust)
plotTSNE(deng, rand_seed = 1, colour_by = "tSNE_kmeans")
```

<div class="figure" style="text-align: center">
<img src="25-clustering_files/figure-html/clust-tsne-kmeans2-1.png" alt="tSNE map of the patient data with 8 colored clusters, identified by the k-means clustering algorithm" width="672" />
<p class="caption">(\#fig:clust-tsne-kmeans2)tSNE map of the patient data with 8 colored clusters, identified by the k-means clustering algorithm</p>
</div>

__Exercise 7__: Make the same plot for $k=10$.

__Exercise 8__: Compare the results between `tSNE+kmeans` and the original publication cell types. Can the results be improved by changing the `perplexity` parameter?

__Our solution__:

```
## [1] 0.3701639
```

As you may have noticed, both `pcaReduce` and `tSNE+kmeans` are stochastic
and give different results every time they are run. To get a better
overview of the solutions, we need to run the methods multiple times. `SC3` is also stochastic, but thanks to the consensus step, it is more robust and less likely to produce different outcomes.

### SNN-Cliq

Here we run SNN-cliq with te default parameters provided in the author's example:


```r
distan <- "euclidean"
par.k <- 3
par.r <- 0.7
par.m <- 0.5
# construct a graph
scRNA.seq.funcs::SNN(
    data = t(input),
    outfile = "snn-cliq.txt",
    k = par.k,
    distance = distan
)
# find clusters in the graph
snn.res <- 
    system(
        paste0(
            "python utils/Cliq.py ", 
            "-i snn-cliq.txt ",
            "-o res-snn-cliq.txt ",
            "-r ", par.r,
            " -m ", par.m
        ),
        intern = TRUE
    )
cat(paste(snn.res, collapse = "\n"))
```

```
## input file snn-cliq.txt
## find 66 quasi-cliques
## merged into 29 clusters
## unique assign done
```

```r
snn.res <- read.table("res-snn-cliq.txt")
# remove files that were created during the analysis
system("rm snn-cliq.txt res-snn-cliq.txt")

colData(deng)$SNNCliq <- as.character(snn.res[,1])
plotPCA(deng, colour_by = "SNNCliq")
```

<img src="25-clustering_files/figure-html/unnamed-chunk-21-1.png" width="672" style="display: block; margin: auto;" />

__Exercise 9__: Compare the results between `SNN-Cliq` and the original publication cell types.

__Our solution__:

```
## [1] 0.2629731
```

### SINCERA

As mentioned in the previous chapter [SINCERA](https://research.cchmc.org/pbge/sincera.html) is based on hierarchical clustering. One important thing to keep in mind is that it performs a gene-level z-score transformation before doing clustering:


```r
# perform gene-by-gene per-sample z-score transformation
dat <- apply(input, 1, function(y) scRNA.seq.funcs::z.transform.helper(y))
# hierarchical clustering
dd <- as.dist((1 - cor(t(dat), method = "pearson"))/2)
hc <- hclust(dd, method = "average")
```

If the number of cluster is not known [SINCERA](https://research.cchmc.org/pbge/sincera.html) can identify __k__ as the minimum height of the hierarchical tree that generates no more than a specified number of singleton clusters (clusters containing only 1 cell)

```r
num.singleton <- 0
kk <- 1
for (i in 2:dim(dat)[2]) {
    clusters <- cutree(hc, k = i)
    clustersizes <- as.data.frame(table(clusters))
    singleton.clusters <- which(clustersizes$Freq < 2)
    if (length(singleton.clusters) <= num.singleton) {
        kk <- i
    } else {
        break;
    }
}
cat(kk)
```

```
## 6
```

Let's now visualize the SINCERA results as a heatmap:

```r
pheatmap(
    t(dat),
    cluster_cols = hc,
    cutree_cols = kk,
    kmeans_k = 100,
    show_rownames = FALSE
)
```

<div class="figure" style="text-align: center">
<img src="25-clustering_files/figure-html/clust-sincera-1.png" alt="Clustering solutions of SINCERA method using found $k$" width="672" />
<p class="caption">(\#fig:clust-sincera)Clustering solutions of SINCERA method using found $k$</p>
</div>

__Exercise 10__: Compare the results between `SINCERA` and the original publication cell types.

__Our solution__:

```
## [1] 0.3823537
```

__Exercise 11__: Is using the singleton cluster criteria for finding __k__ a good idea?

### sessionInfo()


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
##  [1] pheatmap_1.0.8             scater_1.6.3              
##  [3] SingleCellExperiment_1.0.0 SummarizedExperiment_1.8.1
##  [5] DelayedArray_0.4.1         matrixStats_0.53.1        
##  [7] GenomicRanges_1.30.3       GenomeInfoDb_1.14.0       
##  [9] IRanges_2.12.0             S4Vectors_0.16.0          
## [11] ggplot2_2.2.1              SC3_1.7.7                 
## [13] pcaReduce_1.0              mclust_5.4                
## [15] mnormt_1.5-5               pcaMethods_1.70.0         
## [17] Biobase_2.38.0             BiocGenerics_0.24.0       
## [19] knitr_1.20                
## 
## loaded via a namespace (and not attached):
##   [1] Rtsne_0.13             ggbeeswarm_0.6.0       colorspace_1.3-2      
##   [4] rjson_0.2.15           class_7.3-14           rprojroot_1.3-2       
##   [7] XVector_0.18.0         bit64_0.9-7            AnnotationDbi_1.40.0  
##  [10] mvtnorm_1.0-7          scRNA.seq.funcs_0.1.0  codetools_0.2-15      
##  [13] tximport_1.6.0         doParallel_1.0.11      robustbase_0.92-8     
##  [16] cluster_2.0.6          shinydashboard_0.6.1   shiny_1.0.5           
##  [19] rrcov_1.4-3            compiler_3.4.3         httr_1.3.1            
##  [22] backports_1.1.2        assertthat_0.2.0       Matrix_1.2-7.1        
##  [25] lazyeval_0.2.1         limma_3.34.9           htmltools_0.3.6       
##  [28] prettyunits_1.0.2      tools_3.4.3            bindrcpp_0.2          
##  [31] gtable_0.2.0           glue_1.2.0             GenomeInfoDbData_1.0.0
##  [34] reshape2_1.4.3         dplyr_0.7.4            doRNG_1.6.6           
##  [37] Rcpp_0.12.15           gdata_2.18.0           iterators_1.0.9       
##  [40] xfun_0.1               stringr_1.3.0          mime_0.5              
##  [43] hypergeo_1.2-13        rngtools_1.2.4         gtools_3.5.0          
##  [46] WriteXLS_4.0.0         statmod_1.4.30         XML_3.98-1.10         
##  [49] edgeR_3.20.9           DEoptimR_1.0-8         MASS_7.3-45           
##  [52] zlibbioc_1.24.0        scales_0.5.0           rhdf5_2.22.0          
##  [55] RColorBrewer_1.1-2     yaml_2.1.17            memoise_1.1.0         
##  [58] gridExtra_2.3          pkgmaker_0.22          biomaRt_2.34.2        
##  [61] stringi_1.1.6          RSQLite_2.0            highr_0.6             
##  [64] pcaPP_1.9-73           foreach_1.4.4          orthopolynom_1.0-5    
##  [67] e1071_1.6-8            contfrac_1.1-11        caTools_1.17.1        
##  [70] moments_0.14           rlang_0.2.0            pkgconfig_2.0.1       
##  [73] bitops_1.0-6           evaluate_0.10.1        lattice_0.20-34       
##  [76] ROCR_1.0-7             bindr_0.1              labeling_0.3          
##  [79] cowplot_0.9.2          bit_1.1-12             deSolve_1.20          
##  [82] plyr_1.8.4             magrittr_1.5           bookdown_0.7          
##  [85] R6_2.2.2               gplots_3.0.1           DBI_0.7               
##  [88] pillar_1.2.1           RCurl_1.95-4.10        tibble_1.4.2          
##  [91] KernSmooth_2.23-15     rmarkdown_1.8          viridis_0.5.0         
##  [94] progress_1.1.2         locfit_1.5-9.1         grid_3.4.3            
##  [97] data.table_1.10.4-3    blob_1.1.0             digest_0.6.15         
## [100] xtable_1.8-2           httpuv_1.3.6.1         elliptic_1.3-7        
## [103] munsell_0.4.3          registry_0.5           beeswarm_0.2.3        
## [106] viridisLite_0.3.0      vipor_0.4.5
```
