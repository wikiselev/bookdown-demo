FROM rocker/r-base

RUN apt-get update \
        && apt-get install -y --no-install-recommends \
            libcurl4-gnutls-dev \
            libssl-dev \
            libxml2-dev \
            libcairo2 \
            pandoc \
            pandoc-citeproc \
            r-cran-rjava \
            python \
            texlive-full

# install R packages
RUN Rscript -e "install.packages('devtools')"

RUN Rscript -e "install.packages('bookdown')"
RUN Rscript -e "install.packages('knitr')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('BiocInstaller')"
RUN Rscript -e "devtools::install_github('hemberg-lab/scRNA.seq.funcs')"

RUN Rscript -e "install.packages('pheatmap')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('limma')"

RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('scater')"
RUN Rscript -e "install.packages('statmod')"
RUN Rscript -e "install.packages('mvoutlier')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('scran')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('RUVSeq')"

RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('SC3')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('pcaMethods')"
RUN Rscript -e "devtools::install_github('JustinaZ/pcaReduce')"
RUN Rscript -e "devtools::install_github('satijalab/seurat')"
  
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('M3Drop')"

RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('TSCAN')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('monocle')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('destiny')"
RUN Rscript -e "devtools::install_github('jw156605/SLICER')"

RUN Rscript -e "install.packages('ROCR')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('DESeq2')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('edgeR')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('MAST')"

## optional
# RUN Rscript -e "devtools::install_github('hms-dbmi/scde', build_vignettes = FALSE)"

RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('MultiAssayExperiment')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('SummarizedExperiment')"

# add our scripts
ADD . /

# run scripts
CMD bash build.sh
