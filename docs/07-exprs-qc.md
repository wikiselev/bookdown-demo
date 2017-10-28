---
output: html_document
---

## Expression QC (UMI) {#exprs-qc}

### Introduction

Once gene expression has been quantified it is summarized as an __expression matrix__ where each row corresponds to a gene (or transcript) and each column corresponds to a single cell. This matrix should be examined to remove poor quality cells which were not detected in either read QC or mapping QC steps. Failure to remove low quality cells at this
stage may add technical noise which has the potential to obscure
the biological signals of interest in the downstream analysis. 

Since there is currently no standard method for performing scRNASeq the expected values for the various QC measures that will be presented here can vary substantially from experiment to experiment. Thus, to perform QC we will be looking for cells which are outliers with respect to the rest of the dataset rather than comparing to independent quality standards. Consequently, care should be taken when comparing quality metrics across datasets collected using different protocols.


### Tung dataset

To illustrate cell QC, we consider a
[dataset](http://jdblischak.github.io/singleCellSeq/analysis/) of
 induced pluripotent stem cells generated from three different individuals [@Tung2017-ba] in [Yoav Gilad](http://giladlab.uchicago.edu/)'s lab at the
University of Chicago. The experiments were carried out on the
Fluidigm C1 platform and to facilitate the quantification both unique
molecular identifiers (UMIs) and ERCC _spike-ins_ were used. The data files are located in the `tung` folder in your working directory. These files are the copies of the original files made on the 15/03/16. We will use these copies for reproducibility purposes.




```r
library(SingleCellExperiment)
library(scater)
options(stringsAsFactors = FALSE)
```

Load the data and annotations:

```r
molecules <- read.table("tung/molecules.txt", sep = "\t")
anno <- read.table("tung/annotation.txt", sep = "\t", header = TRUE)
```

Inspect a small portion of the expression matrix

```r
head(molecules[ , 1:3])
```

```
##                 NA19098.r1.A01 NA19098.r1.A02 NA19098.r1.A03
## ENSG00000237683              0              0              0
## ENSG00000187634              0              0              0
## ENSG00000188976              3              6              1
## ENSG00000187961              0              0              0
## ENSG00000187583              0              0              0
## ENSG00000187642              0              0              0
```

```r
head(anno)
```

```
##   individual replicate well      batch      sample_id
## 1    NA19098        r1  A01 NA19098.r1 NA19098.r1.A01
## 2    NA19098        r1  A02 NA19098.r1 NA19098.r1.A02
## 3    NA19098        r1  A03 NA19098.r1 NA19098.r1.A03
## 4    NA19098        r1  A04 NA19098.r1 NA19098.r1.A04
## 5    NA19098        r1  A05 NA19098.r1 NA19098.r1.A05
## 6    NA19098        r1  A06 NA19098.r1 NA19098.r1.A06
```

The data consists of 3 individuals and 3 replicates and therefore has 9 batches in total.

We standardize the analysis by using both `SingleCellExperiment` (SCE) and `scater` packages. First, create the SCE object:

```r
umi <- SingleCellExperiment(
    assays = list(counts = as.matrix(molecules)), 
    colData = anno
)
```

Remove genes that are not expressed in any cell:

```r
keep_feature <- rowSums(counts(umi) > 0) > 0
umi <- umi[keep_feature, ]
```

Define control features (genes) - ERCC spike-ins and mitochondrial genes ([provided](http://jdblischak.github.io/singleCellSeq/analysis/qc-filter-ipsc.html) by the authors):

```r
isSpike(umi, "ERCC") <- grepl("^ERCC-", rownames(umi))
isSpike(umi, "MT") <- rownames(umi) %in% 
    c("ENSG00000198899", "ENSG00000198727", "ENSG00000198888",
    "ENSG00000198886", "ENSG00000212907", "ENSG00000198786",
    "ENSG00000198695", "ENSG00000198712", "ENSG00000198804",
    "ENSG00000198763", "ENSG00000228253", "ENSG00000198938",
    "ENSG00000198840")
```

Calculate the quality metrics:

```r
umi <- calculateQCMetrics(
    umi,
    feature_controls = list(
        ERCC = isSpike(umi, "ERCC"), 
        MT = isSpike(umi, "MT")
    )
)
```


### Cell QC

#### Library size

Next we consider the total number of RNA molecules detected per
sample (if we were using read counts rather than UMI counts this would
be the total number of reads). Wells with few reads/molecules are likely to have been broken or failed to capture a cell, and should thus be removed.


```r
hist(
    umi$total_counts,
    breaks = 100
)
abline(v = 25000, col = "red")
```

<div class="figure" style="text-align: center">
<img src="07-exprs-qc_files/figure-html/total-counts-hist-1.png" alt="Histogram of library sizes for all cells" width="90%" />
<p class="caption">(\#fig:total-counts-hist)Histogram of library sizes for all cells</p>
</div>

__Exercise 1__

1. How many cells does our filter remove?

2. What distribution do you expect that the
total number of molecules for each cell should follow?

__Our answer__


```
## filter_by_total_counts
## FALSE  TRUE 
##    46   818
```

#### Detected genes

In addition to ensuring sufficient sequencing depth for each sample, we also want to make sure that the reads are distributed across the transcriptome. Thus, we count the total number of unique genes detected in each sample.


```r
hist(
    umi$total_features,
    breaks = 100
)
abline(v = 7000, col = "red")
```

<div class="figure" style="text-align: center">
<img src="07-exprs-qc_files/figure-html/total-features-hist-1.png" alt="Histogram of the number of detected genes in all cells" width="90%" />
<p class="caption">(\#fig:total-features-hist)Histogram of the number of detected genes in all cells</p>
</div>

From the plot we conclude that most cells have between 7,000-10,000 detected genes,
which is normal for high-depth scRNA-seq. However, this varies by
experimental protocol and sequencing depth. For example, droplet-based methods
or samples with lower sequencing-depth typically detect fewer genes per cell. The most notable feature in the above plot is the __"heavy tail"__ on the left hand side of the
distribution. If detection rates were equal across the cells then the
distribution should be approximately normal. Thus we remove those
cells in the tail of the distribution (fewer than 7,000 detected genes).

__Exercise 2__

How many cells does our filter remove?

__Our answer__


```
## filter_by_expr_features
## FALSE  TRUE 
##   116   748
```

#### ERCCs and MTs

Another measure of cell quality is the ratio between ERCC _spike-in_
RNAs and endogenous RNAs. This ratio can be used to estimate the total amount
of RNA in the captured cells. Cells with a high level of _spike-in_ RNAs
had low starting amounts of RNA, likely due to the cell being
dead or stressed which may result in the RNA being degraded.


```r
plotPhenoData(
    umi,
    aes_string(
        x = "total_features",
        y = "pct_counts_MT",
        colour = "batch"
    )
)
```

<div class="figure" style="text-align: center">
<img src="07-exprs-qc_files/figure-html/mt-vs-counts-1.png" alt="Percentage of counts in MT genes" width="90%" />
<p class="caption">(\#fig:mt-vs-counts)Percentage of counts in MT genes</p>
</div>


```r
plotPhenoData(
    umi,
    aes_string(
        x = "total_features",
        y = "pct_counts_ERCC",
        colour = "batch"
    )
)
```

<div class="figure" style="text-align: center">
<img src="07-exprs-qc_files/figure-html/ercc-vs-counts-1.png" alt="Percentage of counts in ERCCs" width="90%" />
<p class="caption">(\#fig:ercc-vs-counts)Percentage of counts in ERCCs</p>
</div>

The above analysis shows that majority of the cells from NA19098.r2 batch have a very high ERCC/Endo ratio. Indeed, it has been shown by the authors that this batch contains cells of smaller size. 

__Exercise 3__

Create filters for removing batch NA19098.r2 and cells with high expression of mitochondrial genes (>10% of total counts in a cell).

__Our answer__


```
## filter_by_ERCC
## FALSE  TRUE 
##    96   768
```

```
## filter_by_MT
## FALSE  TRUE 
##    31   833
```

__Exercise 4__

What would you expect to see in the ERCC vs counts plot if you were examining a dataset containing cells of different sizes (eg. normal & senescent cells)?

__Answer__

You would expect to see a group corresponding to the smaller cells (normal) with a higher fraction of ERCC reads than a separate group corresponding to the larger cells (senescent).

### Cell filtering

#### Manual

Now we can define a cell filter based on our previous analysis:


```r
umi$use <- (
    # sufficient features (genes)
    filter_by_expr_features &
    # sufficient molecules counted
    filter_by_total_counts &
    # sufficient endogenous RNA
    filter_by_ERCC &
    # remove cells with unusual number of reads in MT genes
    filter_by_MT
)
```


```r
table(umi$use)
```

```
## 
## FALSE  TRUE 
##   207   657
```

#### Automatic

Another option available in `scater` is to conduct PCA on a set of QC metrics and then use automatic outlier detection to identify potentially problematic cells. 

By default, the following metrics are used for PCA-based outlier detection:

* __pct_counts_top_100_features__
* __total_features__
* __pct_counts_feature_controls__
* __n_detected_feature_controls__
* __log10_counts_endogenous_features__
* __log10_counts_feature_controls__

`scater` first creates a matrix where the rows represent cells and the columns represent the different QC metrics. Here, the PCA plot provides a 2D representation of cells ordered by their quality metrics. The outliers are then detected using methods from the mvoutlier package. 


```r
umi <- plotPCA(
    umi,
    size_by = "total_features", 
    shape_by = "use",
    pca_data_input = "pdata",
    detect_outliers = TRUE,
    return_SCE = TRUE
)
```

<div class="figure" style="text-align: center">
<img src="07-exprs-qc_files/figure-html/auto-cell-filt-1.png" alt="PCA plot used for automatic detection of cell outliers" width="90%" />
<p class="caption">(\#fig:auto-cell-filt)PCA plot used for automatic detection of cell outliers</p>
</div>


```r
table(umi$outlier)
```

```
## 
## FALSE  TRUE 
##   819    45
```

### Compare filterings

__Exercise 5__

Compare the default, automatic and manual cell filters. Plot a Venn diagram of the outlier cells from these filterings.

__Hint__: Use `vennCounts` and `vennDiagram` functions from the [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) package to make a Venn diagram.

__Answer__


```
## 
## Attaching package: 'limma'
```

```
## The following object is masked from 'package:scater':
## 
##     plotMDS
```

```
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
```

<div class="figure" style="text-align: center">
<img src="07-exprs-qc_files/figure-html/cell-filt-comp-1.png" alt="Comparison of the default, automatic and manual cell filters" width="90%" />
<p class="caption">(\#fig:cell-filt-comp)Comparison of the default, automatic and manual cell filters</p>
</div>

### Gene analysis

#### Gene expression

In addition to removing cells with poor quality, it is usually a good idea to exclude genes where we suspect that technical artefacts may have skewed the results. Moreover, inspection of the gene expression profiles may provide insights about how the experimental procedures could be improved.

It is often instructive to consider the number of reads consumed by the top 50 expressed genes.


```r
plotQC(umi, type = "highest-expression")
```

<div class="figure" style="text-align: center">
<img src="07-exprs-qc_files/figure-html/top50-gene-expr-1.png" alt="Number of total counts consumed by the top 50 expressed genes" width="90%" />
<p class="caption">(\#fig:top50-gene-expr)Number of total counts consumed by the top 50 expressed genes</p>
</div>

The distributions are relatively flat indicating (but not guaranteeing!) good coverage of the full transcriptome of these cells. However, there are several spike-ins in the top 15 genes which suggests a greater dilution of the spike-ins may be preferrable if the experiment is to be repeated.


#### Gene filtering

It is typically a good idea to remove genes whose expression level is considered __"undetectable"__. We define a gene as  detectable if at least two cells contain more than 1 transcript from the gene. If we were considering read counts rather than UMI counts a reasonable threshold is to require at least five reads in at least two cells. However, in both cases the threshold strongly depends on the sequencing depth. It is important to keep in mind that genes must be filtered after cell filtering since some genes may only be detected in poor quality cells (__note__ `colData(umi)$use` filter applied to the `umi` dataset).


```r
filter_genes <- apply(
    counts(umi[ , colData(umi)$use]), 
    1, 
    function(x) length(x[x > 1]) >= 2
)
rowData(umi)$use <- filter_genes
```


```r
table(filter_genes)
```

```
## filter_genes
## FALSE  TRUE 
##  4660 14066
```

Depending on the cell-type, protocol and sequencing depth, other cut-offs may be appropriate.


### Save the data

Dimensions of the QCed dataset (do not forget about the gene filter we defined above):

```r
dim(umi[rowData(umi)$use, colData(umi)$use])
```

```
## [1] 14066   657
```

Let's create an additional slot with log-transformed counts (we will need it in the next chapters) and remove saved PCA results from the `reducedDim` slot:

```r
assay(umi, "logcounts_raw") <- log2(counts(umi) + 1)
reducedDim(umi) <- NULL
```

Save the data:

```r
saveRDS(umi, file = "tung/umi.rds")
```

### Big Exercise

Perform exactly the same QC analysis with read counts of the same Blischak data. Use `tung/reads.txt` file to load the reads. Once you have finished please compare your results to ours (next chapter).

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
## [1] parallel  stats4    methods   stats     graphics  grDevices utils    
## [8] datasets  base     
## 
## other attached packages:
##  [1] limma_3.32.10               scater_1.5.21              
##  [3] ggplot2_2.2.1               SingleCellExperiment_0.99.4
##  [5] SummarizedExperiment_1.6.5  DelayedArray_0.2.7         
##  [7] matrixStats_0.52.2          Biobase_2.36.2             
##  [9] GenomicRanges_1.28.6        GenomeInfoDb_1.12.3        
## [11] IRanges_2.10.5              S4Vectors_0.14.7           
## [13] BiocGenerics_0.22.1         knitr_1.17                 
## 
## loaded via a namespace (and not attached):
##   [1] ggbeeswarm_0.6.0        minqa_1.2.4            
##   [3] colorspace_1.3-2        mvoutlier_2.0.8        
##   [5] rjson_0.2.15            modeltools_0.2-21      
##   [7] class_7.3-14            mclust_5.3             
##   [9] rprojroot_1.2           XVector_0.16.0         
##  [11] pls_2.6-0               cvTools_0.3.2          
##  [13] MatrixModels_0.4-1      flexmix_2.3-14         
##  [15] bit64_0.9-7             AnnotationDbi_1.38.2   
##  [17] mvtnorm_1.0-6           sROC_0.1-2             
##  [19] splines_3.4.2           tximport_1.4.0         
##  [21] robustbase_0.92-7       nloptr_1.0.4           
##  [23] robCompositions_2.0.6   pbkrtest_0.4-7         
##  [25] kernlab_0.9-25          cluster_2.0.6          
##  [27] shinydashboard_0.6.1    shiny_1.0.5            
##  [29] rrcov_1.4-3             compiler_3.4.2         
##  [31] backports_1.1.1         assertthat_0.2.0       
##  [33] Matrix_1.2-7.1          lazyeval_0.2.0         
##  [35] htmltools_0.3.6         quantreg_5.34          
##  [37] tools_3.4.2             bindrcpp_0.2           
##  [39] gtable_0.2.0            glue_1.1.1             
##  [41] GenomeInfoDbData_0.99.0 reshape2_1.4.2         
##  [43] dplyr_0.7.4             Rcpp_0.12.13           
##  [45] trimcluster_0.1-2       sgeostat_1.0-27        
##  [47] nlme_3.1-129            fpc_2.1-10             
##  [49] lmtest_0.9-35           laeken_0.4.6           
##  [51] stringr_1.2.0           lme4_1.1-14            
##  [53] mime_0.5                XML_3.98-1.9           
##  [55] edgeR_3.18.1            DEoptimR_1.0-8         
##  [57] zoo_1.8-0               zlibbioc_1.22.0        
##  [59] MASS_7.3-45             scales_0.5.0           
##  [61] VIM_4.7.0               rhdf5_2.20.0           
##  [63] SparseM_1.77            RColorBrewer_1.1-2     
##  [65] yaml_2.1.14             memoise_1.1.0          
##  [67] gridExtra_2.3           biomaRt_2.32.1         
##  [69] reshape_0.8.7           stringi_1.1.5          
##  [71] RSQLite_2.0             highr_0.6              
##  [73] pcaPP_1.9-72            e1071_1.6-8            
##  [75] boot_1.3-18             prabclus_2.2-6         
##  [77] rlang_0.1.2             pkgconfig_2.0.1        
##  [79] bitops_1.0-6            evaluate_0.10.1        
##  [81] lattice_0.20-34         bindr_0.1              
##  [83] labeling_0.3            cowplot_0.8.0          
##  [85] bit_1.1-12              GGally_1.3.2           
##  [87] plyr_1.8.4              magrittr_1.5           
##  [89] bookdown_0.5            R6_2.2.2               
##  [91] DBI_0.7                 mgcv_1.8-22            
##  [93] RCurl_1.95-4.8          sp_1.2-5               
##  [95] nnet_7.3-12             tibble_1.3.4           
##  [97] car_2.1-5               rmarkdown_1.6          
##  [99] viridis_0.4.0           locfit_1.5-9.1         
## [101] grid_3.4.2              data.table_1.10.4-3    
## [103] blob_1.1.0              diptest_0.75-7         
## [105] vcd_1.4-3               digest_0.6.12          
## [107] xtable_1.8-2            httpuv_1.3.5           
## [109] munsell_0.4.3           beeswarm_0.2.3         
## [111] viridisLite_0.2.0       vipor_0.4.5
```
