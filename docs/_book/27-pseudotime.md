---
# knit: bookdown::preview_chapter
output: html_document
---

## Pseudotime analysis



```r
library(SingleCellExperiment)
library(TSCAN)
library(M3Drop)
library(monocle)
library(destiny)
library(SLICER)
library(ouija)
library(scater)
library(ggplot2)
library(ggthemes)
library(ggbeeswarm)
library(corrplot)
set.seed(1)
```

In many situations, one is studying a process where cells change
continuously. This includes, for example, many differentiation processes
taking place during development: following a stimulus, cells
will change from one cell-type to another. Ideally, we would like to
monitor the expression levels of an individual cell over
time. Unfortunately, such monitoring is not possible with scRNA-seq
since the cell is lysed (destroyed) when the RNA is extracted.

Instead, we must sample at multiple time-points and obtain snapshots
of the gene expression profiles. Since some of the cells will proceed
faster along the differentiation than others, each snapshot may
contain cells at varying points along the developmental
progression. We use statistical methods to order the cells along one
or more trajectories which represent the underlying developmental
trajectories, this ordering is referred to as "pseudotime".

In this chapter we will consider five different tools: Monocle, TSCAN,
destiny, SLICER and ouija for ordering cells according to their pseudotime
development. To illustrate the methods we will be using a dataset on
mouse embryonic development [@Deng2014-mx]. The dataset consists of
268 cells from 10 different time-points of early mouse development. In this case, there is no need for pseudotime alignment since the cell labels provide information about the development trajectory. Thus, the labels allow us to establish a ground truth so that we can evaluate and compare the different methods.

A recent review by Cannoodt et al provides a detailed summary of the
various computational methods for trajectory inference from
single-cell transcriptomics [@Cannoodt2016-uj]. They discuss several
tools, but unfortunately for our purposes many of these tools do not
have complete or well-maintained implementations, and/or are not
implemented in R.

Cannoodt et al cover:

* SCUBA - Matlab implementation
* Wanderlust - Matlab (and requires registration to even download)
* Wishbone - Python
* SLICER - R, but package only available on Github
* SCOUP - C++ command line tool
* Waterfall - R, but one R script in supplement
* Mpath - R pkg, but available as tar.gz on Github; function
documentation but no vignette/workflow
* Monocle - Bioconductor package
* TSCAN - Bioconductor package

Unfortunately only two tools discussed (Monocle and TSCAN) meet the
gold standard of open-source software hosted in a reputable repository.

The following figures from the paper summarise some of the features of
the various tools.

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{figures/cannoodt_pseudotime_properties} 

}

\caption{Descriptions of trajectory inference methods for single-cell transcriptomics data (Fig. 2 from Cannoodt et al, 2016).}(\#fig:pseudotime-methods-description)
\end{figure}

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{figures/cannoodt_pseudotime_methods} 

}

\caption{Characterization of trajectory inference methods for single-cell transcriptomics data (Fig. 3 from Cannoodt et al, 2016).}(\#fig:pseudotime-methods)
\end{figure}


### First look at Deng data

Let us take a first look at the Deng data, without yet applying sophisticated pseudotime methods. As the plot below shows, simple PCA does a very good job of displaying the structure in these data. It is only once we reach the blast cell types ("earlyblast", "midblast", "lateblast") that PCA struggles to separate the distinct cell types.


```r
deng_SCE <- readRDS("deng/deng-reads.rds")
deng_SCE$cell_type2 <- factor(
    deng_SCE$cell_type2,
    levels = c("zy", "early2cell", "mid2cell", "late2cell",
                        "4cell", "8cell", "16cell", "earlyblast",
                        "midblast", "lateblast")
)
cellLabels <- deng_SCE$cell_type2
deng <- counts(deng_SCE)
colnames(deng) <- cellLabels
deng_SCE <- plotPCA(deng_SCE, colour_by = "cell_type2", 
                    return_SCE = TRUE)
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/data-overview-1} \end{center}

PCA, here, provides a useful baseline for assessing different pseudotime methods. For a very naive pseudotime we can just take the co-ordinates of the first principal component.


```r
deng_SCE$PC1 <- reducedDim(deng_SCE, "PCA")[,1]
ggplot(as.data.frame(colData(deng_SCE)), aes(x = PC1, y = cell_type2, 
                              colour = cell_type2)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_tableau() + theme_classic() +
    xlab("First principal component") + ylab("Timepoint") +
    ggtitle("Cells ordered by first principal component")
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/pca-pseudotime-1} \end{center}

As the plot above shows, PC1 struggles to correctly order cells early and late in the developmental timecourse, but overall does a relatively good job of ordering cells by developmental time.   

Can bespoke pseudotime methods do better than naive application of PCA?
 

### TSCAN

TSCAN combines clustering with pseudotime analysis. First it clusters the cells using `mclust`, which is based on a mixture of normal distributions. Then it builds a minimum spanning tree to connect the clusters. The branch of this tree that connects the largest number of clusters is the main branch which is used to determine pseudotime.

First we will try to use all genes to order the cells.


```r
procdeng <- TSCAN::preprocess(deng)
colnames(procdeng) <- 1:ncol(deng)
dengclust <- TSCAN::exprmclust(procdeng, clusternum = 10)
TSCAN::plotmclust(dengclust)
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/tscan-all-genes-1} \end{center}

```r
dengorderTSCAN <- TSCAN::TSCANorder(dengclust, orderonly = FALSE)
pseudotime_order_tscan <- as.character(dengorderTSCAN$sample_name)
deng_SCE$pseudotime_order_tscan <- NA
deng_SCE$pseudotime_order_tscan[as.numeric(dengorderTSCAN$sample_name)] <- 
    dengorderTSCAN$Pseudotime
```

Frustratingly, TSCAN only provides pseudotime values for 221 of 268 cells, silently returning missing values for non-assigned cells.

Again, we examine which timepoints have been assigned to each state:


```r
cellLabels[dengclust$clusterid == 10]
```

```
##  [1] late2cell late2cell late2cell late2cell late2cell late2cell late2cell
##  [8] late2cell late2cell late2cell
## 10 Levels: zy early2cell mid2cell late2cell 4cell 8cell ... lateblast
```

```r
ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_order_tscan, 
           y = cell_type2, colour = cell_type2)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_tableau() + theme_classic() +
    xlab("TSCAN pseudotime") + ylab("Timepoint") +
    ggtitle("Cells ordered by TSCAN pseudotime")
```

```
## Warning: Removed 47 rows containing missing values (position_quasirandom).
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/tscan-vs-truth-1} \end{center}

TSCAN gets the development trajectory the "wrong way around", in the sense that later pseudotime values correspond to early timepoints and vice versa. This is not inherently a problem (it is easy enough to reverse the ordering to get the intuitive interpretation of pseudotime), but overall it would be a stretch to suggest that TSCAN performs better than PCA on this dataset. (As it is a PCA-based method, perhaps this is not entirely surprising.)


__Exercise 1__ Compare results for different numbers of clusters (`clusternum`).

### monocle

Monocle skips the clustering stage of TSCAN and directly builds a
minimum spanning tree on a reduced dimension representation of the
cells to connect all cells. Monocle then identifies the longest path
in this tree to determine pseudotime. If the data contains diverging
trajectories (i.e. one cell type differentiates into two different
cell-types), monocle can identify these. Each of the resulting forked paths is
defined as a separate cell state.

Unfortunately, Monocle does not work when all the genes are used, so
we must carry out feature selection. First, we use M3Drop:

```r
m3dGenes <- as.character(
    M3DropFeatureSelection(deng)$Gene
)
```

```
## Warning in bg__calc_variables(expr_mat): Warning: Removing 1134 undetected
## genes.
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/m3d-select-genes-1} \end{center}

```r
d <- deng[which(rownames(deng) %in% m3dGenes), ]
d <- d[!duplicated(rownames(d)), ]
```

Now run monocle:

```r
colnames(d) <- 1:ncol(d)
geneNames <- rownames(d)
rownames(d) <- 1:nrow(d)
pd <- data.frame(timepoint = cellLabels)
pd <- new("AnnotatedDataFrame", data=pd)
fd <- data.frame(gene_short_name = geneNames)
fd <- new("AnnotatedDataFrame", data=fd)

dCellData <- newCellDataSet(d, phenoData = pd, featureData = fd, expressionFamily = tobit())
dCellData <- setOrderingFilter(dCellData, which(geneNames %in% m3dGenes))
dCellData <- estimateSizeFactors(dCellData)
dCellDataSet <- reduceDimension(dCellData, pseudo_expr = 1)
dCellDataSet <- orderCells(dCellDataSet, reverse = FALSE)
plot_cell_trajectory(dCellDataSet)
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/monocle-all-genes-1} \end{center}

```r
# Store the ordering
pseudotime_monocle <-
    data.frame(
        Timepoint = phenoData(dCellDataSet)$timepoint,
        pseudotime = phenoData(dCellDataSet)$Pseudotime,
        State = phenoData(dCellDataSet)$State
    )
rownames(pseudotime_monocle) <- 1:ncol(d)
pseudotime_order_monocle <-
    rownames(pseudotime_monocle[order(pseudotime_monocle$pseudotime), ])
```

We can again compare the inferred pseudotime to the known sampling timepoints.

```r
deng_SCE$pseudotime_monocle <- pseudotime_monocle$pseudotime
ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_monocle, 
           y = cell_type2, colour = cell_type2)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_tableau() + theme_classic() +
    xlab("monocle pseudotime") + ylab("Timepoint") +
    ggtitle("Cells ordered by monocle pseudotime")
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/monocle-vs-truth-1} \end{center}

Monocle - at least with its default settings - performs poorly on these data. The "late2cell" group is completely separated from the "zy", "early2cell" and "mid2cell" cells (though these are correctly ordered), and there is no separation at all of "4cell", "8cell", "16cell" or any blast cell groups.


### Diffusion maps

[Diffusion maps](https://en.wikipedia.org/wiki/Diffusion_map) were introduced by [Ronald Coifman and Stephane Lafon](http://www.sciencedirect.com/science/article/pii/S1063520306000546), and the underlying idea is to assume that the data are samples from a diffusion process. The method infers the low-dimensional manifold by estimating the eigenvalues and eigenvectors for the diffusion operator related to the data.

[Haghverdi et al](http://biorxiv.org/content/biorxiv/early/2015/08/04/023309.full.pdf) have applied the diffusion maps concept to the analysis of single-cell RNA-seq data to create an R package called [destiny](http://bioconductor.org/packages/destiny).

We will take the ranko prder of cells in the first diffusion map component as "diffusion map pseudotime" here.


```r
deng <- logcounts(deng_SCE)
colnames(deng) <- cellLabels
dm <- DiffusionMap(t(deng))

tmp <- data.frame(DC1 = eigenvectors(dm)[,1],
                  DC2 = eigenvectors(dm)[,2],
                  Timepoint = deng_SCE$cell_type2)
ggplot(tmp, aes(x = DC1, y = DC2, colour = Timepoint)) +
    geom_point() + scale_color_tableau() + 
    xlab("Diffusion component 1") + 
    ylab("Diffusion component 2") +
    theme_classic()
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/destiny-deng-1} \end{center}

```r
deng_SCE$pseudotime_diffusionmap <- rank(eigenvectors(dm)[,1])
ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_diffusionmap, 
           y = cell_type2, colour = cell_type2)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_tableau() + theme_classic() +
    xlab("Diffusion map pseudotime (first diffusion map component)") +
    ylab("Timepoint") +
    ggtitle("Cells ordered by diffusion map pseudotime")
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/destiny-deng-2} \end{center}

Like the other methods, using the first diffusion map component from destiny as pseudotime does a good job at ordering the early time-points (if we take high values as "earlier" in developement), but it is unable to distinguish the later ones.

__Exercise 2__ Do you get a better resolution between the later time points by considering additional eigenvectors?

__Exercise 3__ How does the ordering change if you only use the genes identified by M3Drop?


### SLICER

The SLICER method is an algorithm for constructing trajectories that
describe gene expression changes during a sequential biological
process, just as Monocle and TSCAN are. SLICER is designed to capture
highly nonlinear gene expression changes, automatically select genes
related to the process, and detect multiple branch and loop features
in the trajectory [@Welch2016-jr]. The SLICER R package is available
from its [GitHub repository](https://github.com/jw156605/SLICER) and
can be installed from there using the `devtools` package.

We use the `select_genes` function in SLICER to automatically select
the genes to use in builing the cell trajectory. The function uses
"neighbourhood variance" to identify genes that vary smoothly, rather
than fluctuating randomly, across the set of cells. Following this, we
determine which value of "k" (number of nearest neighbours) yields an embedding that
most resembles a trajectory. Then we estimate the [locally linear
embedding](https://en.wikipedia.org/wiki/Nonlinear_dimensionality_reduction) of the cells.


```r
library("lle")
slicer_genes <- select_genes(t(deng))
k <- select_k(t(deng[slicer_genes,]), kmin = 30, kmax=60)
```

```
## finding neighbours
## calculating weights
## computing coordinates
## finding neighbours
## calculating weights
## computing coordinates
## finding neighbours
## calculating weights
## computing coordinates
## finding neighbours
## calculating weights
## computing coordinates
## finding neighbours
## calculating weights
## computing coordinates
## finding neighbours
## calculating weights
## computing coordinates
## finding neighbours
## calculating weights
## computing coordinates
```

```r
slicer_traj_lle <- lle(t(deng[slicer_genes,]), m = 2, k)$Y
```

```
## finding neighbours
## calculating weights
## computing coordinates
```

```r
reducedDim(deng_SCE, "LLE") <- slicer_traj_lle
plotReducedDim(deng_SCE, use_dimred = "LLE", colour_by = "cell_type2") +
    xlab("LLE component 1") + ylab("LLE component 2") +
    ggtitle("Locally linear embedding of cells from SLICER")
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/slicer-analyis-1} \end{center}

With the locally linear embedding computed we can construct a
k-nearest neighbour graph that is fully connected. This plot displays
a (yellow) circle for each cell, with the cell ID number overlaid in
blue. Here we show the graph computed using 10 nearest
neighbours. Here, SLICER appears to detect one major trajectory with
one branch.


```r
slicer_traj_graph <- conn_knn_graph(slicer_traj_lle, 10)
plot(slicer_traj_graph, main = "Fully connected kNN graph from SLICER")
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/slicer-build-graph-1} \end{center}

From this graph we can identify "extreme" cells that are candidates
for start/end cells in the trajectory.


```r
ends <- find_extreme_cells(slicer_traj_graph, slicer_traj_lle)
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/slicer-1} \end{center}

```r
start <- ends[1]
```

Having defined a start cell we can order the cells in the estimated pseudotime.


```r
pseudotime_order_slicer <- cell_order(slicer_traj_graph, start)
branches <- assign_branches(slicer_traj_graph, start)

pseudotime_slicer <-
    data.frame(
        Timepoint = cellLabels,
        pseudotime = NA,
        State = branches
    )
pseudotime_slicer$pseudotime[pseudotime_order_slicer] <-
    1:length(pseudotime_order_slicer)
deng_SCE$pseudotime_slicer <- pseudotime_slicer$pseudotime
```

We can again compare the inferred pseudotime to the known sampling
timepoints. SLICER does not provide a pseudotime value per se, just an
ordering of cells.


```r
ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_slicer, 
           y = cell_type2, colour = cell_type2)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_tableau() + theme_classic() +
    xlab("SLICER pseudotime (cell ordering)") +
    ylab("Timepoint") +
    theme_classic()
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/slicer-vs-truth-1} \end{center}

Like the previous method, SLICER here provides a good ordering for the
early time points. It places "16cell" cells before "8cell" cells, but provides better ordering for blast cells than many of the earlier methods.

__Exercise 4__ How do the results change for different k? (e.g. k = 5) What about changing the number of nearest neighbours in
the call to `conn_knn_graph`?

__Exercise 5__ How does the ordering change if you use a different set
of genes from those chosen by SLICER (e.g. the genes identified by M3Drop)?

### Ouija

Ouija (http://kieranrcampbell.github.io/ouija/) takes a different approach from the pseudotime estimation methods we have looked at so far. Earlier methods have all been "unsupervised", which is to say that apart from perhaps selecting informative genes we do not supply the method with any prior information about how we expect certain genes or the trajectory as a whole to behave. 

Ouija, in contrast, is a probabilistic framework that allows for interpretable learning of single-cell pseudotimes using only small panels of marker genes. This method:

* infers pseudotimes from a small number of marker genes letting you understand why the pseudotimes have been learned in terms of those genes;
* provides parameter estimates (with uncertainty) for interpretable gene regulation behaviour (such as the peak time or the upregulation time); 
* has a Bayesian hypothesis test to find genes regulated before others along the trajectory; 
* identifies metastable states, ie discrete cell types along the continuous trajectory.

We will supply the following marker genes to Ouija (with timepoints where they are expected to be highly expressed):

* Early timepoints: Dazl, Rnf17, Sycp3, Nanog, Pou5f1, Fgf8, Egfr, Bmp5, Bmp15
* Mid timepoints: Zscan4b, Foxa1, Prdm14, Sox21
* Late timepoints: Creb3, Gpx4, Krt8, Elf5, Eomes, Cdx2, Tdgf1, Gdf3

With Ouija we can model genes as either exhibiting monotonic up or down regulation (known as switch-like behaviour), or transient behaviour where the gene briefly peaks. By default, Ouija assumes all genes exhibit switch-like behaviour (the authors assure us not to worry if we get it wrong - the noise model means incorrectly specifying a transient gene as switch-like has minimal effect).

Here we can "cheat" a little and check that our selected marker genes do actually identify different timepoints of the differentiation process.


```r
ouija_markers_down <- c("Dazl", "Rnf17", "Sycp3", "Fgf8", 
                        "Egfr", "Bmp5", "Bmp15", "Pou5f1")
ouija_markers_up <- c("Creb3", "Gpx4", "Krt8", "Elf5", "Cdx2", 
                      "Tdgf1", "Gdf3", "Eomes")
ouija_markers_transient <- c("Zscan4b", "Foxa1", "Prdm14", "Sox21")
ouija_markers <- c(ouija_markers_down, ouija_markers_up, 
                   ouija_markers_transient)
plotExpression(deng_SCE, ouija_markers, x = "cell_type2", colour_by = "cell_type2") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/ouija-response-type-1} \end{center}

In order to fit the pseudotimes wesimply call `ouija`, passing in the expected response types. Note that if no response types are provided then they are all assumed to be switch-like by default, which we will do here. The input to Ouija can be a cell-by-gene matrix of non-negative expression values, or an ExpressionSet object, or, happily, by selecting the `logcounts` values from a SingleCellExperiment object.

We can apply prior information about whether genes are up- or down-regulated across the differentiation process, and also provide prior information about when the switch in expression or a peak in expression is likely to occur. 

We can fit the Ouija model using either:

* Hamiltonian Monte Carlo (HMC) - full MCMC inference where gradient information of the log-posterior is used to “guide” the random walk through the parameter space, or
* Automatic Differentiation Variational Bayes (ADVI or simply VI) - approximate inference where the KL divergence to an approximate distribution is minimised.

In general, HMC will provide more accurate inference with approximately correct posterior variance for all parameters. However, VB is orders of magnitude quicker than HMC and while it may underestimate posterior variance, the Ouija authors suggest that anecdotally it often performs as well as HMC for discovering posterior pseudotimes.

To help the Ouija model, we provide it with prior information about the strength of switches for up- and down-regulated genes. By setting switch strength to -10 for down-regulated genes and 10 for up-regulated genes with a prior strength standard deviation of 0.5 we are telling the model that we are confident about the expected behaviour of these genes across the differentiation process.



```r
options(mc.cores = parallel::detectCores())
response_type <- c(rep("switch", length(ouija_markers_down) + 
                           length(ouija_markers_up)), 
                   rep("transient", length(ouija_markers_transient)))
switch_strengths <- c(rep(-10, length(ouija_markers_down)),
                      rep(10, length(ouija_markers_up)))
switch_strength_sd <- c(rep(0.5, length(ouija_markers_down)),
                      rep(0.5, length(ouija_markers_up)))
garbage <- capture.output(
    oui_vb <- ouija(deng_SCE[ouija_markers,],
                    single_cell_experiment_assay = "logcounts", 
                    response_type = response_type,
                    switch_strengths = switch_strengths,
                    switch_strength_sd = switch_strength_sd,
                    inference_type = "vb")
)

print(oui_vb)
```

```
## A Ouija fit with 268 cells and 20 marker genes 
## Inference type:  Variational Bayes 
## (Gene behaviour) Switch/transient: 16 / 4
```

We can plot the gene expression over pseudotime along with the maximum a posteriori (MAP) estimates of the mean function (the sigmoid or Gaussian transient function) using the plot_expression function. 


```r
plot_expression(oui_vb)
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/ouija-plot-exprs-1} \end{center}

We can also visualise when in the trajectory gene regulation behaviour occurs, either in the form of the switch time or the peak time (for switch-like or transient genes) using the plot_switch_times and plot_transient_times functions:


```r
plot_switch_times(oui_vb)
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/ouija-plot-switch-times-1} \end{center}

```r
plot_peak_times(oui_vb)
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/ouija-plot-switch-times-2} \end{center}

Identify metastable states using consistency matrices.


```r
cmo <- consistency_matrix(oui_vb)
plot_consistency(oui_vb)
```

```
## Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font
## metrics unknown for Unicode character U+2192

## Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font
## metrics unknown for Unicode character U+2192
```

```
## Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): conversion
## failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted for <e2>
```

```
## Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): conversion
## failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted for <86>
```

```
## Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): conversion
## failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted for <92>
```

```
## Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font
## metrics unknown for Unicode character U+2192

## Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): font
## metrics unknown for Unicode character U+2192
```

```
## Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): conversion
## failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted for <e2>
```

```
## Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): conversion
## failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted for <86>
```

```
## Warning in grid.Call(C_stringMetric, as.graphicsAnnot(x$label)): conversion
## failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted for <92>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <e2>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <86>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <92>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <e2>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <86>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <92>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <e2>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <86>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <92>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <e2>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <86>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <92>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <e2>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <86>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <92>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <e2>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <86>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <92>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <e2>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <86>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <92>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <e2>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <86>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <92>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <e2>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <86>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <92>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <e2>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <86>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <92>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <e2>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <86>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <92>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <e2>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <86>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <92>
```

```
## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x
## $y, : conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot
## substituted for <e2>
```

```
## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x
## $y, : conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot
## substituted for <86>
```

```
## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x
## $y, : conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot
## substituted for <92>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <e2>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <86>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <92>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <e2>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <86>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot substituted
## for <92>
```

```
## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x
## $y, : conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot
## substituted for <e2>
```

```
## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x
## $y, : conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot
## substituted for <86>
```

```
## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x
## $y, : conversion failure on 'Pseudotime order →' in 'mbcsToSbcs': dot
## substituted for <92>
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/ouija-consistency-1} \end{center}

```r
cell_classifications <- cluster_consistency(cmo)
```


```r
map_pst <- map_pseudotime(oui_vb)
ouija_pseudotime <- data.frame(map_pst, cell_classifications)

ggplot(ouija_pseudotime, aes(x = map_pst, y = cell_classifications)) +
  geom_point() +
  xlab("MAP pseudotime") +
  ylab("Cell classification")
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/ouija-pseudotime-1} \end{center}

```r
deng_SCE$pseudotime_ouija <- ouija_pseudotime$map_pst
deng_SCE$ouija_cell_class <- ouija_pseudotime$cell_classifications

ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_ouija, 
           y = cell_type2, colour = cell_type2)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_tableau() + theme_classic() +
    xlab("Ouija pseudotime") +
    ylab("Timepoint") +
    theme_classic()
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/ouija-pseudotime-2} \end{center}

Ouija does quite well in the ordering of the cells here, although it can be sensitive to the choice of marker genes and prior information supplied. How do the results change if you select different marker genes or change the priors?

Ouija identifies four metastable states here, which we might annotate as  "zygote/2cell", "4/8/16 cell", "blast1" and "blast2".


```r
ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = as.factor(ouija_cell_class), 
           y = pseudotime_ouija, colour = cell_type2)) +
    geom_boxplot() + 
    coord_flip() +
    scale_color_tableau() + theme_classic() +
    xlab("Ouija cell classification") +
    ylab("Ouija pseudotime") +
    theme_classic()
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/ouija-states-1} \end{center}

A common analysis is to work out the regulation orderings of genes. For example, is gene A upregulated before gene B? Does gene C peak before the downregulation of gene D? Ouija answers these questions in terms of a Bayesian hypothesis test of whether the difference in regulation timing (either switch time or peak time) is significantly different to 0. This is collated using the gene_regulation function.


```r
gene_regs <- gene_regulation(oui_vb)
head(gene_regs)
```

```
## # A tibble: 6 x 7
## # Groups:   label, gene_A [6]
##   label        gene_A gene_B mean_difference lower_95 upper_95 significant
##   <chr>        <chr>  <chr>            <dbl>    <dbl>    <dbl> <lgl>      
## 1 Bmp15 - Cdx2 Bmp15  Cdx2           -0.0434 -0.0912   0.0110  FALSE      
## 2 Bmp15 - Cre~ Bmp15  Creb3           0.278   0.220    0.330   TRUE       
## 3 Bmp15 - Elf5 Bmp15  Elf5           -0.656  -0.688   -0.618   TRUE       
## 4 Bmp15 - Eom~ Bmp15  Eomes           0.0766  0.00433  0.153   TRUE       
## 5 Bmp15 - Fox~ Bmp15  Foxa1          -0.0266 -0.0611   0.00287 FALSE      
## 6 Bmp15 - Gdf3 Bmp15  Gdf3            0.0704  0.00963  0.134   TRUE
```

What conclusions can you draw from the gene regulation output from Ouija?

If you have time, you might try the HMC inference method and see if that changes the Ouija results in any way.


### Comparison of the methods

How do the trajectories inferred by TSCAN, Monocle, Diffusion Map, SLICER and Ouija compare?

TSCAN and Diffusion Map methods get the trajectory the "wrong way round", so we'll adjust that for these comparisons.


```r
df_pseudotime <- as.data.frame(
    colData(deng_SCE)[, grep("pseudotime", colnames(colData(deng_SCE)))]
)
colnames(df_pseudotime) <- gsub("pseudotime_", "", 
                                colnames(df_pseudotime))
df_pseudotime$PC1 <- deng_SCE$PC1
df_pseudotime$order_tscan <- -df_pseudotime$order_tscan
df_pseudotime$diffusionmap <- -df_pseudotime$diffusionmap

corrplot.mixed(cor(df_pseudotime, use = "na.or.complete"), 
               order = "hclust", tl.col = "black",
               main = "Correlation matrix for pseudotime results",
               mar = c(0, 0, 3.1, 0))
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/compare-results-1} \end{center}

We see here that Ouija, TSCAN and SLICER all give trajectories that are similar and strongly correlated with PC1. Diffusion Map is less strongly correlated with these methods, and Monocle gives very different results.


__Exercise 6__: Compare destiny and SLICER to TSCAN, Monocle and Ouija in more depth. Where and how do they differ?

### Expression of genes through time

Each package also enables the visualization of expression through pseudotime. Following individual genes is very helpful for identifying genes that play an important role in the differentiation process. We illustrate the procedure using the `Rhoa` gene.

We have added the pseudotime values computed with all methods here to
the `colData` slot of an `SCE` object. Having done that, the full
plotting capabilities of the `scater` package can be used to
investigate relationships between gene expression, cell populations
and pseudotime. This is particularly useful for the packages such as
SLICER that do not provide plotting functions.


__Principal components__

```r
plotExpression(deng_SCE, "Rhoa", x = "PC1", 
               colour_by = "cell_type2", show_violin = FALSE,
               show_smooth = TRUE)
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/Rhoa-pc1-1} \end{center}

__TSCAN__

```r
plotExpression(deng_SCE, "Rhoa", x = "pseudotime_order_tscan", 
               colour_by = "cell_type2", show_violin = FALSE,
               show_smooth = TRUE)
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/Rhoa-tscan-1} \end{center}

__Monocle__

```r
plotExpression(deng_SCE, "Rhoa", x = "pseudotime_monocle", 
               colour_by = "cell_type2", show_violin = FALSE,
               show_smooth = TRUE)
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/Rhoa-monocle-1} \end{center}

__Diffusion Map__

```r
plotExpression(deng_SCE, "Rhoa", x = "pseudotime_diffusionmap", 
               colour_by = "cell_type2", show_violin = FALSE,
               show_smooth = TRUE)
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/Rhoa-diff-map-1} \end{center}


__SLICER__

```r
plotExpression(deng_SCE, "Rhoa", x = "pseudotime_slicer", 
               colour_by = "cell_type2", show_violin = FALSE,
               show_smooth = TRUE)
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/Rhoa-slicer-1} \end{center}

__Ouija__

```r
plotExpression(deng_SCE, "Rhoa", x = "pseudotime_ouija", 
               colour_by = "cell_type2", show_violin = FALSE,
               show_smooth = TRUE)
```



\begin{center}\includegraphics{27-pseudotime_files/figure-latex/Rhoa-ouija-1} \end{center}

How many of these methods outperform the naive approach of using the first principal component to represent pseudotime for these data?

__Exercise 7__: Repeat the exercise using a subset of the genes, e.g. the set of highly variable genes that can be obtained using `Brennecke_getVariableGenes()`

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
##  [1] splines   parallel  stats4    methods   stats     graphics  grDevices
##  [8] utils     datasets  base     
## 
## other attached packages:
##  [1] bindrcpp_0.2               rstan_2.17.3              
##  [3] StanHeaders_2.17.2         lle_1.1                   
##  [5] snowfall_1.84-6.1          snow_0.4-2                
##  [7] MASS_7.3-45                scatterplot3d_0.3-40      
##  [9] corrplot_0.84              ggbeeswarm_0.6.0          
## [11] ggthemes_3.4.0             scater_1.6.3              
## [13] ouija_0.99.0               Rcpp_0.12.15              
## [15] SLICER_0.2.0               destiny_2.6.1             
## [17] monocle_2.6.3              DDRTree_0.1.5             
## [19] irlba_2.3.2                VGAM_1.0-5                
## [21] ggplot2_2.2.1              Matrix_1.2-7.1            
## [23] M3Drop_3.05.00             numDeriv_2016.8-1         
## [25] TSCAN_1.16.0               SingleCellExperiment_1.0.0
## [27] SummarizedExperiment_1.8.1 DelayedArray_0.4.1        
## [29] matrixStats_0.53.1         Biobase_2.38.0            
## [31] GenomicRanges_1.30.3       GenomeInfoDb_1.14.0       
## [33] IRanges_2.12.0             S4Vectors_0.16.0          
## [35] BiocGenerics_0.24.0        knitr_1.20                
## 
## loaded via a namespace (and not attached):
##   [1] utf8_1.1.3             shinydashboard_0.6.1   R.utils_2.6.0         
##   [4] tidyselect_0.2.4       lme4_1.1-15            RSQLite_2.0           
##   [7] AnnotationDbi_1.40.0   htmlwidgets_1.0        grid_3.4.3            
##  [10] combinat_0.0-8         Rtsne_0.13             munsell_0.4.3         
##  [13] codetools_0.2-15       statmod_1.4.30         colorspace_1.3-2      
##  [16] fastICA_1.2-1          rstudioapi_0.7         robustbase_0.92-8     
##  [19] vcd_1.4-4              tensor_1.5             VIM_4.7.0             
##  [22] TTR_0.23-3             labeling_0.3           slam_0.1-42           
##  [25] splancs_2.01-40        tximport_1.6.0         bbmle_1.0.20          
##  [28] GenomeInfoDbData_1.0.0 polyclip_1.6-1         bit64_0.9-7           
##  [31] pheatmap_1.0.8         rhdf5_2.22.0           rprojroot_1.3-2       
##  [34] coda_0.19-1            xfun_0.1               R6_2.2.2              
##  [37] RcppEigen_0.3.3.4.0    locfit_1.5-9.1         bitops_1.0-6          
##  [40] spatstat.utils_1.8-0   assertthat_0.2.0       scales_0.5.0          
##  [43] nnet_7.3-12            beeswarm_0.2.3         gtable_0.2.0          
##  [46] goftest_1.1-1          rlang_0.2.0            MatrixModels_0.4-1    
##  [49] lazyeval_0.2.1         acepack_1.4.1          checkmate_1.8.5       
##  [52] inline_0.3.14          yaml_2.1.17            reshape2_1.4.3        
##  [55] abind_1.4-5            backports_1.1.2        httpuv_1.3.6.1        
##  [58] Hmisc_4.1-1            tensorA_0.36           tools_3.4.3           
##  [61] bookdown_0.7           cubature_1.3-11        gplots_3.0.1          
##  [64] RColorBrewer_1.1-2     proxy_0.4-21           MCMCglmm_2.25         
##  [67] plyr_1.8.4             progress_1.1.2         base64enc_0.1-3       
##  [70] zlibbioc_1.24.0        purrr_0.2.4            RCurl_1.95-4.10       
##  [73] densityClust_0.3       prettyunits_1.0.2      rpart_4.1-10          
##  [76] alphahull_2.1          deldir_0.1-14          reldist_1.6-6         
##  [79] viridis_0.5.0          cowplot_0.9.2          zoo_1.8-1             
##  [82] ggrepel_0.7.0          cluster_2.0.6          magrittr_1.5          
##  [85] data.table_1.10.4-3    SparseM_1.77           lmtest_0.9-35         
##  [88] RANN_2.5.1             mime_0.5               evaluate_0.10.1       
##  [91] xtable_1.8-2           XML_3.98-1.10          smoother_1.1          
##  [94] pbkrtest_0.4-7         mclust_5.4             gridExtra_2.3         
##  [97] biomaRt_2.34.2         HSMMSingleCell_0.112.0 compiler_3.4.3        
## [100] tibble_1.4.2           crayon_1.3.4           KernSmooth_2.23-15    
## [103] minqa_1.2.4            R.oo_1.21.0            htmltools_0.3.6       
## [106] mgcv_1.8-23            corpcor_1.6.9          Formula_1.2-2         
## [109] tidyr_0.8.0            DBI_0.7                boot_1.3-18           
## [112] car_2.1-6              cli_1.0.0              sgeostat_1.0-27       
## [115] R.methodsS3_1.7.1      gdata_2.18.0           bindr_0.1             
## [118] igraph_1.1.2           pkgconfig_2.0.1        foreign_0.8-67        
## [121] laeken_0.4.6           sp_1.2-7               vipor_0.4.5           
## [124] XVector_0.18.0         stringr_1.3.0          digest_0.6.15         
## [127] spatstat.data_1.2-0    rmarkdown_1.8          htmlTable_1.11.2      
## [130] edgeR_3.20.9           curl_3.1               shiny_1.0.5           
## [133] gtools_3.5.0           quantreg_5.35          rjson_0.2.15          
## [136] nloptr_1.0.4           nlme_3.1-129           viridisLite_0.3.0     
## [139] limma_3.34.9           pillar_1.2.1           lattice_0.20-34       
## [142] httr_1.3.1             DEoptimR_1.0-8         survival_2.40-1       
## [145] glue_1.2.0             xts_0.10-1             qlcMatrix_0.9.5       
## [148] FNN_1.1                spatstat_1.55-0        bit_1.1-12            
## [151] class_7.3-14           stringi_1.1.6          blob_1.1.0            
## [154] memoise_1.1.0          latticeExtra_0.6-28    caTools_1.17.1        
## [157] dplyr_0.7.4            e1071_1.6-8            ape_5.0               
## [160] tripack_1.3-8
```
