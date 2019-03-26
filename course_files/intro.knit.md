---
output: html_document
---

# Introduction to single-cell RNA-seq



## Bulk RNA-seq

* A major breakthrough (replaced microarrays) in the late 00's and has been widely used since
* Measures the __average expression level__ for each gene across a large population of input cells
* Useful for comparative transcriptomics, e.g. samples of the same tissue from different species
* Useful for quantifying expression signatures from ensembles, e.g. in disease studies
* __Insufficient__ for studying heterogeneous systems, e.g. early development studies, complex tissues (brain)
* Does __not__ provide insights into the stochastic nature of gene expression
    
## scRNA-seq

* A __new__ technology, first publication by [@Tang2009-bu]
* Did not gain widespread popularity until [~2014](https://www.ohio.edu/bioinformatics/upload/Single-Cell-RNA-seq-Method-of-the-Year-2013.pdf) when new protocols and lower sequencing costs made it more accessible
* Measures the __distribution of expression levels__ for each gene across a population of cells
* Allows to study new biological questions in which __cell-specific changes in transcriptome are important__, e.g. cell type identification, heterogeneity of cell responses, stochasticity of gene expression, inference of gene regulatory networks across the cells.
* Datasets range __from $10^2$ to $10^6$ cells__ and increase in size every year
* Currently there are several different protocols in use, e.g. SMART-seq2 [@Picelli2013-sb], CELL-seq [@Hashimshony2012-kd] and Drop-seq [@Macosko2015-ix]
* There are also commercial platforms available, including the [Fluidigm C1](https://www.fluidigm.com/products/c1-system), [Wafergen ICELL8](https://www.wafergen.com/products/icell8-single-cell-system) and the [10X Genomics Chromium](https://www.10xgenomics.com/single-cell/)
* Several computational analysis methods from bulk RNA-seq __can__ be used
* __In most cases__ computational analysis requires adaptation of the existing methods or development of new ones

## Workflow

<div class="figure" style="text-align: center">
<img src="figures/RNA-Seq_workflow-5.pdf.jpg" alt="Single cell sequencing (taken from Wikipedia)" width="90%" />
<p class="caption">(\#fig:intro-rna-seq-workflow)Single cell sequencing (taken from Wikipedia)</p>
</div>

Overall, experimental scRNA-seq protocols are similar to the methods used for bulk RNA-seq. We will be discussing some of the most common approaches in the next chapter.

## Computational Analysis

This course is concerned with the computational analysis of the data
obtained from scRNA-seq experiments. The first steps (yellow) are general for any highthroughput sequencing data. Later steps (orange) require a mix of existing RNASeq analysis methods and novel methods to address the technical difference of scRNASeq. Finally the biological interpretation (blue) __should__ be analyzed with methods specifically developed for scRNASeq.

<div class="figure" style="text-align: center">
<img src="figures/flowchart.png" alt="Flowchart of the scRNA-seq analysis" width="65%" />
<p class="caption">(\#fig:intro-flowchart)Flowchart of the scRNA-seq analysis</p>
</div>

There are several reviews of the scRNA-seq analysis available including [@Stegle2015-uv].

Today, there are also several different platforms available for carrying out one or more steps in the flowchart above. These include:

* [Falco](https://github.com/VCCRI/Falco/) a single-cell RNA-seq processing framework on the cloud.
* [SCONE](https://github.com/YosefLab/scone) (Single-Cell Overview of Normalized Expression), a package for single-cell RNA-seq data quality control and normalization.
* [Seurat](http://satijalab.org/seurat/) is an R package designed for QC, analysis, and exploration of single cell RNA-seq data.
* [ASAP](https://asap.epfl.ch/) (Automated Single-cell Analysis Pipeline) is an interactive web-based platform for single-cell analysis.
* [Bioconductor](https://master.bioconductor.org/packages/release/workflows/html/simpleSingleCell.html) is a open-source, open-development software project for the analysis of high-throughput genomics data, including packages for the analysis of single-cell data.


## Challenges

The main difference between bulk and single cell RNA-seq is that each sequencing library represents a single cell, instead of a population of cells. Therefore, significant attention has to be paid to comparison of the results from different cells (sequencing libraries). The main sources of discrepancy between the libraries are:

* __Amplification__ (up to 1 million fold)
* __Gene 'dropouts'__ in which a gene is observed at a moderate expression level in one cell but is not detected in another cell [@Kharchenko2014-ts].

In both cases the discrepancies are introduced due to low starting amounts of transcripts since the RNA comes from one cell only. Improving the transcript capture efficiency and reducing the amplification bias are currently active areas of research. However, as we shall see in this course, it is possible to alleviate some of these issues through proper normalization and corrections.

