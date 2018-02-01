FROM quay.io/hemberg-group/scrna-seq-course-base:latest

# install R packages
RUN echo 'install.packages(c("devtools", "bookdown", "knitr", "pheatmap", "statmod", "mvoutlier", "mclust", "dplyr", "penalized", "cluster", "Seurat", "KernSmooth", "mgcv", "ROCR", "googleVis", "tidyverse", "ggfortify"))' > /opt/packages.r && \
    echo 'source("https://bioconductor.org/biocLite.R")' >> /opt/packages.r && \
    echo 'biocLite()' >> /opt/packages.r && \
    echo 'biocLite(c("limma", "SingleCellExperiment", "Rhdf5lib", "beachmat", "scater", "scran", "RUVSeq", "sva", "SC3", "pcaMethods", "TSCAN", "monocle", "destiny", "DESeq2", "edgeR", "MAST", "scfind", "scmap", "MultiAssayExperiment", "SummarizedExperiment"))' >> /opt/packages.r && \
    echo 'devtools::install_github(c("hemberg-lab/scRNA.seq.funcs", "Vivianstats/scImpute", "theislab/kBET", "JustinaZ/pcaReduce", "tallulandrews/M3Drop", "jw156605/SLICER"))' >> /opt/packages.r && \
    Rscript /opt/packages.r

# add our scripts
ADD . /home/rstudio/

# run scripts
CMD cd /home/rstudio && bash build.sh
