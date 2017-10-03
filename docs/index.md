--- 
title: "Analysis of single cell RNA-seq data"
author: "Vladimir Kiselev, Tallulah Andrews, Davis McCarthy and Martin Hemberg"
date: "2017-10-03"
knit: "bookdown::render_book"
documentclass: book
bibliography: [book.bib]
biblio-style: apalike
link-citations: yes
always_allow_html: yes
---

# About the course

Today it is possible to obtain genome-wide transcriptome data from single cells using high-throughput sequencing (scRNA-seq). The main advantage of scRNA-seq is that the cellular resolution and the genome wide scope makes it possible to address issues that are intractable using other methods, e.g. bulk RNA-seq or single-cell RT-qPCR. However, to analyze scRNA-seq data, novel methods are required and some of the underlying assumptions for the methods developed for bulk RNA-seq experiments are no longer valid.

In this course we will discuss some of the questions that can be addressed using scRNA-seq as well as the available computational and statistical methods avialable. The course is taught through the University of Cambridge <a href="http://training.csx.cam.ac.uk/bioinformatics/" target="blank">Bioinformatics training unit</a>, but the material found on these pages is meant to be used for anyone interested in learning about computational analysis of scRNA-seq data. The course is taught twice per year and the material here is updated prior to each event.

The number of computational tools is increasing rapidly and we are doing our best to keep up to date with what is available. One of the main constraints for this course is that we would like to use tools that are implemented in R and that run reasonably fast. Moreover, we will also confess to being somewhat biased towards methods that have been developed either by us or by our friends and colleagues. 

## Video

<iframe width="560" height="315" src="https://www.youtube.com/embed/i58Fk6R03PA?list=PLEyKDyF1qdOaoWBu8jNwN4o5z6wlm3aES" frameborder="0" allowfullscreen></iframe>

## Registration  

Please follow this link and register for the __"Analysis of single cell RNA-seq data"__ course:
<a href="http://training.csx.cam.ac.uk/bioinformatics/search" target="blank">http://training.csx.cam.ac.uk/bioinformatics/search</a>

## GitHub
<a href="https://github.com/hemberg-lab/scRNA.seq.course" target="blank">https://github.com/hemberg-lab/scRNA.seq.course</a>

## Docker image

The course can be reproduced without any package installation by running the course docker image which contains all the required packages.

Make sure Docker is installed on your system. If not, please follow [these instructions](https://docs.docker.com/engine/installation/). To run the course docker image:

```
docker run -it quay.io/hemberg-group/scrna-seq-course:latest R
```

It will download the course docker image (may take some time) and start a new R session in a docker container with all packages installed and all data files available.

Note: It is currently not possible to plot figures interactively from the Docker, however figures can still be written to file for inspection.

## License
<b>GPL-3</b>

## Prerequisites

The course is intended for those who have basic familiarity with Unix and the R scripting language.

We will also assume that you are familiar with mapping and analysing bulk RNA-seq data as well as with the commonly available computational tools.

We recommend attending the [Introduction to RNA-seq and ChIP-seq data analysis](http://training.csx.cam.ac.uk/bioinformatics/search) or the [Analysis of high-throughput sequencing data with Bioconductor](http://training.csx.cam.ac.uk/bioinformatics/search) before attending this course.

## Contact

If you have any __comments__, __questions__ or __suggestions__ about the material, please contact <a href="mailto:vladimir.yu.kiselev@gmail.com">Vladimir Kiselev</a>.
