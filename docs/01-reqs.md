---
output: html_document
---

# Technical requirements

## R-based

This course is based on the popular programming language [R](https://www.r-project.org/). However, one of the methods that we describe ([SNN-Cliq](http://bioinfo.uncc.edu/SNNCliq/)) is only partly R-based. It makes a simple _python_ call from R and requires a user to have write permissions to the working directory.

## Docker image

If you do not want to install all the packages required for the course manually, you can run a course docker image which contains all the required packages.

Make sure Docker is installed on your system. If not, please follow [these instructions](https://docs.docker.com/engine/installation/). To run the course docker image:

```
docker run -it quay.io/hemberg-group/scrna-seq-course:latest R
```

It will download the course docker image (may take some time) and start a new R session in a docker container with all packages installed and all data files available.

## Manual installation

If you are not using a docker image of the course, then to be able to run all code chunks of the course you need to clone or download the [course GitHub repository](https://github.com/hemberg-lab/scRNA.seq.course) and start an R session in the cloned folder. You will also need to install the following R packages (ordered by purposes):

### General

[devtools](https://cran.r-project.org/web/packages/devtools/index.html)

```r
install.packages("devtools")
```

[BiocInstaller](http://bioconductor.org/packages/BiocInstaller)

```r
source('https://bioconductor.org/biocLite.R')
biocLite('BiocInstaller')
```

[scRNA.seq.funcs](https://github.com/hemberg-lab/scRNA.seq.funcs) - R package containing some special functions used in this course:

```r
devtools::install_github("hemberg-lab/scRNA.seq.funcs")
```

### Plotting

[pheatmap](https://cran.r-project.org/web/packages/pheatmap/index.html)

```r
install.packages("pheatmap")
```

[limma](https://bioconductor.org/packages/limma)

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("limma")
```

### QC and normalisation

[scater](http://bioconductor.org/packages/scater)

```r
source('https://bioconductor.org/biocLite.R')
biocLite('scater')
```

[mvoutlier](https://cran.r-project.org/web/packages/mvoutlier/index.html) - for an automatic outlier detection used by the [scater](https://github.com/davismcc/scater) package.

```r
install.packages("mvoutlier")
```

[statmod](https://cran.r-project.org/web/packages/statmod/index.html) - a dependency for [mvoutlier](https://cran.r-project.org/web/packages/mvoutlier/index.html).

```r
install.packages("statmod")
```

[scran](http://bioconductor.org/packages/scran)

```r
source('https://bioconductor.org/biocLite.R')
biocLite('scran')
```

[RUVSeq](https://bioconductor.org/packages/RUVSeq)

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("RUVSeq")
```

### Clustering

[pcaReduce](https://github.com/JustinaZ/pcaReduce)

```r
devtools::install_github("JustinaZ/pcaReduce")
```

[pcaMethods](http://bioconductor.org/pcaMethods) is a [pcaReduce](https://github.com/JustinaZ/pcaReduce) dependency:

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("pcaMethods")
```

[SC3](http://bioconductor.org/packages/SC3)

```r
source("https://bioconductor.org/biocLite.R")
biocLite("SC3")
```

[SEURAT](https://github.com/satijalab/seurat)

```r
devtools::install_github('satijalab/seurat')
```

### Dropouts

[M3Drop](https://github.com/tallulandrews/M3Drop)

```r
devtools::install_github('tallulandrews/M3Drop')
```

### Pseudotime

[TSCAN](https://bioconductor.org/packages/TSCAN)

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("TSCAN")
```

[monocle](https://bioconductor.org/packages/monocle)

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("monocle")
```

[destiny](https://bioconductor.org/packages/destiny)

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("destiny")
```

[SLICER](https://github.com/jw156605/SLICER)

```r
devtools::install_github('jw156605/SLICER')
```

### Differential Expression

[ROCR](https://cran.r-project.org/web/packages/ROCR/index.html)

```r
install.packages("ROCR")
```

[edgeR](https://bioconductor.org/packages/edgeR)

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
```

[DESeq2](https://bioconductor.org/packages/DESeq2)

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
```

[MAST](https://bioconductor.org/packages/MAST)

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("MAST")
```

[scde](http://hms-dbmi.github.io/scde/) (optional)

```r
devtools::install_github("hms-dbmi/scde", build_vignettes = FALSE)
```

* Installation on Mac OS X may require [this additional gfortran library](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/):
```
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```

* See the [help page](http://hms-dbmi.github.io/scde/help.html) for additional support.

## Extra tools

[MultiAssayExperiment](https://bioconductor.org/packages/MultiAssayExperiment) for working with [conquer](http://imlspenticton.uzh.ch:3838/conquer/) datasets:

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("MultiAssayExperiment")
```

[SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment) for working with [conquer](http://imlspenticton.uzh.ch:3838/conquer/) datasets:

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("SummarizedExperiment")
```

