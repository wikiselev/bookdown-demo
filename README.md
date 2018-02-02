# About the course

Today it is possible to obtain genome-wide transcriptome data from single cells using high-throughput sequencing (scRNA-seq). The main advantage of scRNA-seq is that the cellular resolution and the genome wide scope makes it possible to address issues that are intractable using other methods, e.g. bulk RNA-seq or single-cell RT-qPCR. However, to analyze scRNA-seq data, novel methods are required and some of the underlying assumptions for the methods developed for bulk RNA-seq experiments are no longer valid.

In this course we will discuss some of the questions that can be addressed using scRNA-seq as well as the available computational and statistical methods avialable. The course is taught through the University of Cambridge <a href="http://training.csx.cam.ac.uk/bioinformatics/" target="blank">Bioinformatics training unit</a>, but the material found on these pages is meant to be used for anyone interested in learning about computational analysis of scRNA-seq data. The course is taught twice per year and the material here is updated prior to each event.

The number of computational tools is increasing rapidly and we are doing our best to keep up to date with what is available. One of the main constraints for this course is that we would like to use tools that are implemented in R and that run reasonably fast. Moreover, we will also confess to being somewhat biased towards methods that have been developed either by us or by our friends and colleagues. 

## Web page

__[http://hemberg-lab.github.io/scRNA.seq.course](http://hemberg-lab.github.io/scRNA.seq.course)__

## Video

This video was recorded in November 2017, at that time the course contained less chapters than the current version.

__[https://www.youtube.com/embed/IrlNcJwPClQ?list=PLEyKDyF1qdObdFBc3JncwXAnMUHlcd0ap](https://www.youtube.com/embed/IrlNcJwPClQ?list=PLEyKDyF1qdObdFBc3JncwXAnMUHlcd0ap)__

## Registration  

Please follow this link and register for the __"Analysis of single cell RNA-seq data"__ course:
<a href="http://training.csx.cam.ac.uk/bioinformatics/search" target="blank">http://training.csx.cam.ac.uk/bioinformatics/search</a>

## GitHub
<a href="https://github.com/hemberg-lab/scRNA.seq.course" target="blank">https://github.com/hemberg-lab/scRNA.seq.course</a>

## Docker image

[![Docker Repository on Quay](https://quay.io/repository/hemberg-group/scrna-seq-course-rstudio/status "Docker Repository on Quay")](https://quay.io/repository/hemberg-group/scrna-seq-course-rstudio)

The course can be reproduced without any package installation by running the course docker RStudio image which contains all the required packages.

Make sure Docker is installed on your system. If not, please follow [these instructions](https://docs.docker.com/engine/installation/). To run the course RStudio docker image:

```
docker run -d -p 8787:8787 quay.io/hemberg-group/scrna-seq-course-rstudio
```

This downloads the docker image (may take some time) and start a new Rstudio session in a docker container with all packages installed and all data files available.

Then visit `localhost:8787` in your browser and log in with `username:password` as `rstudio:rstudio`. Now you are ready to go!

More details on how ot run RStudio docker with different options can be found [here](https://hub.docker.com/r/rocker/rstudio-stable/).

## Manual installation

If you are not using a docker image of the course, then to be able to run all code chunks of the course you need to clone or download the [course GitHub repository](https://github.com/hemberg-lab/scRNA.seq.course) and start an R session in the cloned folder. You will also need to install all packages listed in the course docker files: [Dockerfile1](https://github.com/hemberg-lab/scRNA.seq.course-base/blob/master/Dockerfile) and [Dockerfile2](https://github.com/hemberg-lab/scRNA.seq.course-R/blob/master/Dockerfile).

Alternatively, you can just install packages listed in a chapter of interest.

## License
All of the course material is licensed under <b>GPL-3</b>. Anyone is welcome to go through the material in order to learn about analysis of scRNA-seq data. If you plan to use the material for your own teaching, we would appreciate if you tell us about it in addition to providing a suitable citation.

## Prerequisites

The course is intended for those who have basic familiarity with Unix and the R scripting language.

We will also assume that you are familiar with mapping and analysing bulk RNA-seq data as well as with the commonly available computational tools.

We recommend attending the [Introduction to RNA-seq and ChIP-seq data analysis](http://training.csx.cam.ac.uk/bioinformatics/search) or the [Analysis of high-throughput sequencing data with Bioconductor](http://training.csx.cam.ac.uk/bioinformatics/search) before attending this course.

## Contact

If you have any __comments__, __questions__ or __suggestions__ about the material, please contact <a href="mailto:vladimir.yu.kiselev@gmail.com">Vladimir Kiselev</a>.
