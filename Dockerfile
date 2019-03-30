FROM quay.io/cellgeni/cellgeni-jupyter:v0.4.9

USER root

# pre-requisites
RUN apt-get update && apt-get install -yq --no-install-recommends \
    libncurses5-dev \
    libncursesw5-dev \
    procps

# Install FastQC
RUN curl -fsSL http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip -o /opt/fastqc_v0.11.5.zip && \
    unzip /opt/fastqc_v0.11.5.zip -d /opt/ && \
    chmod 755 /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc && \
    rm /opt/fastqc_v0.11.5.zip

# Install Kallisto
RUN curl -fsSL https://github.com/pachterlab/kallisto/releases/download/v0.43.1/kallisto_linux-v0.43.1.tar.gz -o /opt/kallisto_linux-v0.43.1.tar.gz && \
    tar xvzf /opt/kallisto_linux-v0.43.1.tar.gz -C /opt/ && \
    ln -s /opt/kallisto_linux-v0.43.1/kallisto /usr/local/bin/kallisto && \
    rm /opt/kallisto_linux-v0.43.1.tar.gz

# Install STAR
RUN git clone https://github.com/alexdobin/STAR.git /opt/STAR && \
    ln -s /opt/STAR/bin/Linux_x86_64/STAR /usr/local/bin/STAR && \
    ln -s /opt/STAR/bin/Linux_x86_64/STARlong /usr/local/bin/STARlong

# Install SAMTools
RUN curl -fsSL https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 -o /opt/samtools-1.9.tar.bz2 && \
    tar xvjf /opt/samtools-1.9.tar.bz2 -C /opt/ && \
    cd /opt/samtools-1.9 && \
    make && \
    make install && \
    rm /opt/samtools-1.9.tar.bz2

# Install featureCounts
RUN curl -fsSL http://downloads.sourceforge.net/project/subread/subread-1.5.1/subread-1.5.1-Linux-x86_64.tar.gz -o /opt/subread-1.5.1-Linux-x86_64.tar.gz && \
    tar xvf /opt/subread-1.5.1-Linux-x86_64.tar.gz -C /opt/ && \
    ln -s /opt/subread-1.5.1-Linux-x86_64/bin/featureCounts /usr/local/bin/featureCounts && \
    rm /opt/subread-1.5.1-Linux-x86_64.tar.gz

# Install cutadapt and MAGIC and awscli (to download data)
RUN pip install cutadapt magic-impute awscli

# Install TrimGalore
RUN mkdir /opt/TrimGalore && \
    curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.4.5.zip -o /opt/TrimGalore/trim_galore_v0.4.5.zip && \
    unzip /opt/TrimGalore/trim_galore_v0.4.5.zip -d /opt/TrimGalore && \
    ln -s /opt/TrimGalore/trim_galore /usr/local/bin/trim_galore && \
    rm /opt/TrimGalore/trim_galore_v0.4.5.zip

# Install bedtools2
RUN curl -fsSL https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz -o /opt/bedtools-2.27.1.tar.gz && \
    tar xvzf /opt/bedtools-2.27.1.tar.gz -C /opt/ && \
    cd /opt/bedtools2 && \
    make && \
    cd - && \
    cp /opt/bedtools2/bin/* /usr/local/bin && \
    rm /opt/bedtools-2.27.1.tar.gz

# install CRAN packages
RUN apt-get update && apt-get install -yq --no-install-recommends \
    r-cran-knitr \
    r-cran-statmod \
    r-cran-mvoutlier \
    r-cran-penalized \
    r-cran-mgcv \
    r-cran-corrplot

# Install other CRAN
RUN Rscript -e 'install.packages(c("bookdown", "cluster", "KernSmooth", "ROCR", "googleVis", "ggbeeswarm", "SLICER", "ggfortify", "mclust"))'

# install github packages
# see here for with_libpaths description:
# https://stackoverflow.com/questions/24646065/how-to-specify-lib-directory-when-installing-development-version-r-packages-from
# (do not install anything in the home directory, it will be wiped out when a volume is mounted to the docker container)
RUN Rscript -e 'withr::with_libpaths(new = "/usr/lib/R/site-library/", devtools::install_github(c("hemberg-lab/scRNA.seq.funcs", "Vivianstats/scImpute", "theislab/kBET", "JustinaZ/pcaReduce", "kieranrcampbell/ouija")))'

# Install Bioconductor packages
RUN Rscript -e 'BiocManager::install(c("MultiAssayExperiment", "SummarizedExperiment"), version = "3.8")'

USER $NB_UID

# add our scripts
ADD course_files /home/jovyan

# download data and extra files from S3
COPY ./poststart.sh /
