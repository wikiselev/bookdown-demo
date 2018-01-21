FROM rocker/rstudio

RUN apt-get update -y --no-install-recommends \ 
    && apt-get -y install -f \
        libssl-dev \
        libcurl4-openssl-dev \
        libxml2-dev \
        libcairo2 \
        pandoc \
        pandoc-citeproc \
        r-cran-rjava \
        libmariadb-client-lgpl-dev \
        python \
#       python3.6 \
        python3-pip \
        python3-tk \
        git
#       texlive-full

# Install STAR
RUN git clone https://github.com/alexdobin/STAR.git /opt/STAR && \
    ln -s /opt/STAR/bin/Linux_x86_64/STAR /usr/local/bin/STAR && \
    ln -s /opt/STAR/bin/Linux_x86_64/STARlong /usr/local/bin/STARlong

# Install SAMTools
RUN curl -fsSL https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 -o /opt/samtools-1.3.1.tar.bz2 && \
    tar xvjf /opt/samtools-1.3.1.tar.bz2 -C /opt/ && \
    cd /opt/samtools-1.3.1 && \
    make && \
    make install && \
    rm /opt/samtools-1.3.1.tar.bz2

# Install featureCounts
RUN curl -fsSL http://downloads.sourceforge.net/project/subread/subread-1.5.1/subread-1.5.1-Linux-x86_64.tar.gz -o /opt/subread-1.5.1-Linux-x86_64.tar.gz && \
    tar xvzf /opt/subread-1.5.1-Linux-x86_64.tar.gz -C /opt/ && \
    ln -s /opt/subread-1.5.1-Linux-x86_64/bin/featureCounts /usr/local/bin/featureCounts && \
    rm /opt/subread-1.5.1-Linux-x86_64.tar.gz

# install MAGIC
RUN git clone git://github.com/pkathail/magic.git \
    && cd magic \
    && pip3 install numpy \
    && pip3 install argparse \
    && pip3 install .

# install R packages
RUN echo 'install.packages(c("devtools", "bookdown", "knitr", "pheatmap", "statmod", "mvoutlier", "mclust", "dplyr", "penalized", "cluster", "Seurat", "KernSmooth", "mgcv", "ROCR", "googleVis"))' > /opt/packages.r && \
    echo 'source("https://bioconductor.org/biocLite.R")' > /opt/packages.r && \
    echo 'biocLite()' >> /opt/packages.r && \
    echo 'biocLite(c("limma", "SingleCellExperiment", "Rhdf5lib", "beachmat", "scater", "scran", "RUVSeq", "sva", "SC3", "pcaMethods", "TSCAN", "monocle", "destiny", "DESeq2", "edgeR", "MAST", "scfind", "scmap", "MultiAssayExperiment", "SummarizedExperiment"))' >> /opt/packages.r && \
    echo 'devtools::install_github(c("hemberg-lab/scRNA.seq.funcs", "Vivianstats/scImpute", "theislab/kBET", "JustinaZ/pcaReduce", "tallulandrews/M3Drop", "jw156605/SLICER"))' > /opt/packages.r && \
    Rscript /opt/packages.r && \
    mkdir -p  /usr/local/lib/R/site-library

# add our scripts
ADD . /home/rstudio/

# run scripts
CMD cd /home/rstudio && bash build.sh
