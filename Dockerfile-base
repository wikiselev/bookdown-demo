FROM rocker/rstudio

RUN apt-get update -y --no-install-recommends && \ 
    apt-get -y install -f \
        libssl-dev \
        libcurl4-openssl-dev \
        libxml2-dev \
        libcairo2 \
        pandoc \
        pandoc-citeproc \
        r-cran-rjava \
        libmariadb-client-lgpl-dev \
        python \
        python3-pip \
        python3-tk \
        git \
        curl \
        texlive-full

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

# Install cutadapt
RUN pip3 install cutadapt

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

# install MAGIC
RUN git clone git://github.com/pkathail/magic.git && \
    cd magic && \
    pip3 install numpy && \
    pip3 install argparse && \
    pip3 install .
