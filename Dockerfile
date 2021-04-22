FROM jupyter/base-notebook

# versions of software
ARG r_cran_version="cran40"
ARG rstudio_version="1.2.5019"
ARG fastqc_version="0.11.9"            
ARG kallisto_version="0.46.2"          
ARG samtools_version="1.11"            
ARG trimgalore_version="0.6.6"         
ARG bedtools_version="2.29.2"           
ARG featurecounts_version="2.0.1"       

#old software verions
#fastqc=0.11.5
#kallisto=0.43.1
#samtools=1.9
#trimgalore=0.4.5
#bedtools=2.27.1
#featurecounts=1.5.1

USER root

# Install OS packages
RUN apt-get update && apt-get install -yq --no-install-recommends \
    # basic packages 
    htop vim emacs unzip git wget nano rsync \
    parallel sshpass lsb-release libhdf5-dev hdf5-tools libigraph0-dev \
    # R dependecies
    libgsl0-dev libxml2-dev libboost-all-dev libssl-dev libhdf5-dev curl \
    libudunits2-dev libgdal-dev libgeos-dev libproj-dev build-essential xorg-dev \
    libreadline-dev libc6-dev libclang-8-dev zlib1g-dev libbz2-dev liblzma-dev \
    libcurl4-openssl-dev libcairo2-dev libpango1.0-dev tcl-dev tk-dev openjdk-8-jdk \
    gfortran libncurses5-dev libncursesw5-dev procps texlive libv8-dev libgit2-dev \
    # RStudio dependencies - https://github.com/rstudio/rstudio/tree/master/dependencies/linux
    gdebi-core ant clang cmake debsigs dpkg-sig expect fakeroot gnupg1 libacl1-dev \
    libattr1-dev libbz2-dev libcap-dev libclang-6.0-dev libclang-dev \
    libegl1-mesa libfuse2 libgl1-mesa-dev libgtk-3-0 libpam-dev libpq-dev \
    libsqlite3-dev libuser1-dev libxslt1-dev lsof ninja-build patchelf pkg-config \
    psmisc rrdtool software-properties-common uuid-dev wget zlib1g-dev && \
    # clean up
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# clean conda cache
RUN  conda clean --all --yes

# give jovyan sudo permissions
RUN sed -i -e "s/Defaults    requiretty.*/ #Defaults    requiretty/g" /etc/sudoers && \
    echo "jovyan ALL= (ALL) NOPASSWD: ALL" >> /etc/sudoers.d/jovyan

# Install R - https://cran.r-project.org/bin/linux/ubuntu/README.html
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    echo "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -c | awk '{print $2}')-$r_cran_version/" | sudo tee -a /etc/apt/sources.list && \
    add-apt-repository ppa:c2d4u.team/c2d4u4.0+ && \
    apt-get update && apt-get install -yq --no-install-recommends \
        r-base \
        r-base-dev \
    && apt-get clean \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
    rm -rf /var/lib/apt/lists/*

# Install RStudio
RUN RSTUDIO_PKG=rstudio-server-${rstudio_version}-amd64.deb && \
    cd /tmp && \
    wget -q https://download2.rstudio.org/server/bionic/amd64/${RSTUDIO_PKG} && \
    gdebi -n /tmp/${RSTUDIO_PKG} && \
    rm /tmp/${RSTUDIO_PKG}
ENV PATH="${PATH}:/usr/lib/rstudio-server/bin"
ENV LD_LIBRARY_PATH="/usr/lib/R/lib:/lib:/usr/lib/x86_64-linux-gnu:/usr/lib/jvm/java-8-openjdk-amd64/jre/lib/amd64/server"

# Install FastQC
RUN curl -fsSL http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${fastqc_version}.zip -o /opt/fastqc_v${fastqc_version}.zip && \
    unzip /opt/fastqc_v${fastqc_version}.zip -d /opt/ && \
    chmod 755 /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc && \
    rm /opt/fastqc_v${fastqc_version}.zip

# Install Kallisto
RUN curl -fsSL https://github.com/pachterlab/kallisto/releases/download/v${kallisto_version}/kallisto_linux-v${kallisto_version}.tar.gz -o /opt/kallisto_linux-v${kallisto_version}.tar.gz && \
    tar xvzf /opt/kallisto_linux-v${kallisto_version}.tar.gz -C /opt/ && \
    ln -s /opt/kallisto_linux-v${kallisto_version}/kallisto /usr/local/bin/kallisto && \
    rm /opt/kallisto_linux-v${kallisto_version}.tar.gz

# Install STAR
RUN git clone https://github.com/alexdobin/STAR.git /opt/STAR && \
    ln -s /opt/STAR/bin/Linux_x86_64/STAR /usr/local/bin/STAR && \
    ln -s /opt/STAR/bin/Linux_x86_64/STARlong /usr/local/bin/STARlong

# Install SAMTools
RUN curl -fsSL https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 -o /opt/samtools-${samtools_version}.tar.bz2 && \
    tar xvjf /opt/samtools-${samtools_version}.tar.bz2 -C /opt/ && \
    cd /opt/samtools-${samtools_version} && \
    make && \
    make install && \
    rm /opt/samtools-${samtools_version}.tar.bz2

# Install featureCounts
RUN curl -fsSL http://downloads.sourceforge.net/project/subread/subread-${featurecounts_version}/subread-${featurecounts_version}-Linux-x86_64.tar.gz -o /opt/subread-${featurecounts_version}-Linux-x86_64.tar.gz && \
    tar xvf /opt/subread-${featurecounts_version}-Linux-x86_64.tar.gz -C /opt/ && \
    ln -s /opt/subread-${featurecounts_version}-Linux-x86_64/bin/featureCounts /usr/local/bin/featureCounts && \
    rm /opt/subread-${featurecounts_version}-Linux-x86_64.tar.gz

# Install TrimGalore
RUN mkdir /opt/TrimGalore && \
    curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/${trimgalore_version}.zip -o /opt/TrimGalore/trim_galore_v${trimgalore_version}.zip && \
    unzip /opt/TrimGalore/trim_galore_v${trimgalore_version}.zip -d /opt/TrimGalore && \
    ln -s /opt/TrimGalore/trim_galore /usr/local/bin/trim_galore && \
    rm /opt/TrimGalore/trim_galore_v${trimgalore_version}.zip

# Install bedtools2
RUN curl -fsSL https://github.com/arq5x/bedtools2/releases/download/v${bedtools_version}/bedtools-${bedtools_version}.tar.gz -o /opt/bedtools-${bedtools_version}.tar.gz && \
    tar xvzf /opt/bedtools-${bedtools_version}.tar.gz -C /opt/ && \
    cd /opt/bedtools2 && \
    make && \
    cd - && \
    cp /opt/bedtools2/bin/* /usr/local/bin && \
    rm /opt/bedtools-${bedtools_version}.tar.gz

# Set Rprofile for binary installs
RUN echo 'options(repos = c(REPO_NAME = "https://packagemanager.rstudio.com/all/__linux__/focal/latest"))' > /home/jovyan/.Rprofile

# Install base R packages
RUN Rscript -e 'install.packages(c("devtools", "IRkernel", "Rmagic", "BiocManager"))' && \
    Rscript -e 'IRkernel::installspec(prefix="/opt/conda")' && \
    Rscript -e 'BiocManager::install()'

# Install other CRAN
RUN Rscript -e 'install.packages(c( \
    "Seurat", "rJava", "umap", "bookdown", "cluster", "KernSmooth", \
    "ROCR", "googleVis", "ggbeeswarm", "SLICER", "ggfortify", \
    "tidyverse", "pheatmap", "plyr", "dplyr", "readr", "reshape", \
    "reshape2", "reticulate", "viridis", "ggplot2", "ggthemes", "cowplot", \
    "ggforce", "ggridges", "ggrepel", "gplots", "igraph", "car", \
    "ggpubr", "httpuv", "xtable", "sourcetools", "modeltools", "R.oo", \
    "R.methodsS3", "shiny", "later", "checkmate", "bibtex", "lsei", \
    "bit", "segmented", "mclust", "flexmix", "prabclus", "diptest", "mvtnorm", \
    "robustbase", "kernlab", "trimcluster", "proxy", "R.utils", "htmlwidgets", \
    "hexbin", "crosstalk", "promises", "acepack", "zoo", "npsurv", "iterators", \
    "snow", "bit64", "permute", "mixtools", "lars", "ica", "fpc", "ape", \
    "pbapply", "irlba", "dtw", "plotly", "metap", "lmtest", "fitdistrplus", "png", \
    "foreach", "vegan", "tidyr", "withr", "magrittr", "rmpi", "knitr", \
    "statmod", "mvoutlier", "penalized", "mgcv", "corrplot", \
    "lsa", "uwot", "optparse", "DrImpute", "alluvial"))'

# Install Bioconductor packages
RUN Rscript -e 'BiocManager::install(c( \
    "graph", "RBGL", "gtools", "xtable", "pcaMethods", "limma", "SingleCellExperiment", \
    "Rhdf5lib", "scater", "scran", "RUVSeq", "sva", "SC3", "TSCAN", "monocle", "destiny", \
    "DESeq2", "edgeR", "MAST", "scmap", "biomaRt", "MultiAssayExperiment", "SummarizedExperiment", \
    "beachmat", "DropletUtils", "EnsDb.Hsapiens.v86", "batchelor"))'

# Install github packages
RUN Rscript -e 'devtools::install_github(c( \
    "immunogenomics/harmony", "tallulandrews/M3Drop", "hemberg-lab/scRNA.seq.funcs", \
    "Vivianstats/scImpute", "theislab/kBET", "kieranrcampbell/ouija", "hemberg-lab/scfind", \
    "cole-trapnell-lab/monocle3"))'

# Install python packages
RUN pip install --upgrade --no-cache \
    cutadapt magic-impute awscli==1.16.14 \
    jupyter-server-proxy \
    jupyter-rsession-proxy rpy2

# JupyterLab extension to launch registered applications in the python package
RUN jupyter labextension install @jupyterlab/server-proxy

# fix permissions
RUN fix-permissions $CONDA_DIR && \
    fix-permissions /home/$NB_USER

# download data and extra files from S3
#COPY ./poststart.sh /home/jovyan

# add course files
#COPY course_files /home/jovyan
#RUN chmod -R 777 /home/jovyan

# maker user sudo
ENV GRANT_SUDO yes
# use jupyter lab by default
ENV JUPYTER_ENABLE_LAB yes

USER $NB_UID
