FROM rocker/r-base

RUN apt-get update -y --no-install-recommends \ 
        && apt-get -y install -f \
            libssl-dev \
            libcurl4-openssl-dev \
            libxml2-dev \
            libcairo2 \
            pandoc \
            pandoc-citeproc \
            r-cran-rjava \
            python \
            python3.6 \
            python3-pip \
            git
#            texlive-full

# install MAGIC
RUN git clone git://github.com/pkathail/magic.git \
        && cd magic \
        && pip3 install numpy \
        && pip3 install argparse \
        && pip3 install .

# install R packages
RUN Rscript -e "install.packages('devtools')"

RUN Rscript -e "install.packages('bookdown')"
RUN Rscript -e "install.packages('knitr')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('BiocInstaller')"
RUN Rscript -e "devtools::install_github('hemberg-lab/scRNA.seq.funcs')"

RUN Rscript -e "install.packages('pheatmap')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('limma')"

RUN Rscript -e "devtools::install_github('drisso/SingleCellExperiment')"
RUN Rscript -e "devtools::install_github('grimbough/Rhdf5lib')"
RUN Rscript -e "devtools::install_github('LTLA/beachmat')"
# RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('scater')"
RUN Rscript -e "devtools::install_github('davismcc/scater')"
RUN Rscript -e "install.packages('statmod')"
RUN Rscript -e "install.packages('mvoutlier')"
# RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('scran')"
RUN Rscript -e "devtools::install_github('MarioniLab/scran')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('RUVSeq')"

RUN Rscript -e "install.packages('irr')"

RUN Rscript -e "install.packages('penalized')"
RUN Rscript -e "devtools::install_github('Vivianstats/scImpute')"

RUN Rscript -e "devtools::install_github('theislab/kBET')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('sva')"

# RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('SC3')"
RUN Rscript -e "devtools::install_github('hemberg-lab/SC3')"
RUN Rscript -e "source('https://bioconductor.org/biocLite.R');biocLite('pcaMethods')"
RUN Rscript -e "devtools::install_github('JustinaZ/pcaReduce')"
RUN Rscript -e "install.packages('Seurat')"
  
RUN Rscript -e "devtools::install_github('tallulandrews/M3Drop')"

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
