---
output: html_document
---

# Biological Analysis

## Clustering Introduction



Once we have normalized the data and removed confounders we can carry out analyses that are relevant to the biological questions at hand. The exact nature of the analysis depends on the dataset. Nevertheless, there are a few aspects that are useful in a wide range of contexts and we will be discussing some of them in the next few chapters. We will start with the clustering of scRNA-seq data.

### Introduction

One of the most promising applications of scRNA-seq is _de novo_ discovery
and annotation of cell-types based on transcription
profiles. Computationally, this is a hard problem as it amounts to
__unsupervised clustering__. That is, we need to identify groups of
cells based on the similarities of the transcriptomes without any
prior knowledge of the labels. Moreover, in most situations we do not even know the number of clusters _a priori_. The problem is made even more challenging
due to the high level of noise (both technical and biological) and the large number of dimensions
(i.e. genes). 

### Dimensionality reductions

When working with large datasets, it can often be beneficial to apply
some sort of dimensionality reduction method. By projecting
the data onto a lower-dimensional sub-space, one is often able to
significantly reduce the amount of noise. An additional benefit is
that it is typically much easier to visualize the data in a 2 or
3-dimensional subspace. We have already discussed PCA (chapter \@ref(visual-pca)) and t-SNE (chapter \@ref(visual-pca)).

### Clustering methods

__Unsupervised clustering__ is useful in many different applications and
it has been widely studied in machine learning. Some of the most
popular approaches are __hierarchical clustering__, __k-means clustering__ and __graph-based clustering__.

#### Hierarchical clustering

In [hierarchical clustering](https://en.wikipedia.org/wiki/Hierarchical_clustering), one can use either a bottom-up or a
top-down approach. In the former case, each cell is initially assigned to
its own cluster and pairs of clusters are subsequently merged to
create a hieararchy:

<div class="figure" style="text-align: center">
<img src="figures/hierarchical_clustering1.png" alt="Raw data" width="30%" />
<p class="caption">(\#fig:clust-hierarch-raw)Raw data</p>
</div>

<div class="figure" style="text-align: center">
<img src="figures/hierarchical_clustering2.png" alt="The hierarchical clustering dendrogram" width="50%" />
<p class="caption">(\#fig:clust-hierarch-dendr)The hierarchical clustering dendrogram</p>
</div>

With a top-down strategy, one instead starts with
all observations in one cluster and then recursively split each
cluster to form a hierarchy. One of the
advantages of this strategy is that the method is deterministic.

#### k-means

In [_k_-means clustering](https://en.wikipedia.org/wiki/K-means_clustering), the goal is to partition _N_ cells into _k_
different clusters. In an iterative manner, cluster centers are
assigned and each cell is assigned to its nearest cluster:

<div class="figure" style="text-align: center">
<img src="figures/k-means.png" alt="Schematic representation of the k-means clustering" width="100%" />
<p class="caption">(\#fig:clust-k-means)Schematic representation of the k-means clustering</p>
</div>

Most methods for scRNA-seq analysis includes a _k_-means step at some point.

#### Graph-based methods

Over the last two decades there has been a lot of interest in
analyzing networks in various domains. One goal is to identify groups
or modules of nodes in a network.

<div class="figure" style="text-align: center">
<img src="figures/graph_network.jpg" alt="Schematic representation of the graph network" width="100%" />
<p class="caption">(\#fig:clust-graph)Schematic representation of the graph network</p>
</div>

Some of these methods can be applied
to scRNA-seq data by building a graph where each node represents a cell. Note that constructing the graph and assigning weights to the edges is not trivial. One advantage of graph-based methods is that some of them are very efficient and can be applied to networks containing millions of nodes.

### Challenges in clustering

* What is the number of clusters _k_?
* What is a cell type?
* __Scalability__: in the last few years the number of cells in scRNA-seq experiments has grown by several orders of magnitude from ~$10^2$ to ~$10^6$
* Tools are not user-friendly

### Tools for scRNA-seq data

#### [SINCERA](https://research.cchmc.org/pbge/sincera.html)

* SINCERA [@Guo2015-ok] is based on hierarchical clustering
* Data is converted to _z_-scores before clustering
* Identify _k_ by finding the first singleton cluster in the hierarchy

#### [pcaReduce](https://github.com/JustinaZ/pcaReduce)

pcaReduce [@Zurauskiene2016-kg] combines PCA, _k_-means and “iterative” hierarchical clustering. Starting from a large number of clusters pcaReduce iteratively merges similar clusters; after each merging event it removes the principle component explaning the least variance in the data.

#### [SC3](http://bioconductor.org/packages/SC3/)

<div class="figure" style="text-align: center">
<img src="figures/sc3.png" alt="SC3 pipeline" width="100%" />
<p class="caption">(\#fig:clust-sc3)SC3 pipeline</p>
</div>

* SC3 [@Kiselev2016-bq] is based on PCA and spectral dimensionality reductions
* Utilises _k_-means
* Additionally performs the consensus clustering

#### tSNE + k-means

* Based on __tSNE__ maps
* Utilises _k_-means

#### [SNN-Cliq](http://bioinfo.uncc.edu/SNNCliq/)

`SNN-Cliq` [@Xu2015-vf] is a graph-based method. First the method identifies the k-nearest-neighbours of each cell according to the _distance_ measure. This is used to calculate the number of Shared Nearest Neighbours (SNN) between each pair of cells. A graph is built by placing an edge between two cells If they have at least one SNN. Clusters are defined as groups of cells with many edges between them using a "clique" method. SNN-Cliq requires several parameters to be defined manually.

#### Seurat clustering

[`Seurat`](https://github.com/satijalab/seurat) clustering is based on a _community detection_ approach similar to `SNN-Cliq` and to one previously proposed for analyzing CyTOF data [@Levine2015-fk]. Since `Seurat` has become more like an all-in-one tool for scRNA-seq data analysis we dedicate a separate chapter to discuss it in more details (chapter \@ref(seurat-chapter)).

### Comparing clustering

To compare two sets of clustering labels we can use [adjusted Rand index](https://en.wikipedia.org/wiki/Rand_index). The index is a measure of the similarity between two data clusterings. Values of the adjusted Rand index lie in $[0;1]$ interval, where $1$ means that two clusterings are identical and $0$ means the level of similarity expected by chance.
