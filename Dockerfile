# Set the base image to u
FROM ubuntu:22.10

ENV DEBIAN_FRONTEND noninteractive

ARG SAMTOOLSVER=1.15
ARG FASTQC_VER="0.11.9"
ARG BBTOOLSVER=38.96
ARG GATK_VER=4.2.6.0

RUN apt-get update && apt-get upgrade  -y && apt-get install -y wget
RUN apt-get update && apt-get install -y unzip
RUN apt-get update && apt-get install -y git
RUN apt-get update && apt-get install -y cmake
RUN apt-get update && apt-get install -y build-essential
RUN apt-get update && apt-get install -y gcc-multilib

RUN apt-get update && apt-get install -y perl
RUN apt-get update && apt-get install -y curl

RUN apt-get update && \
    apt-get install --no-install-recommends -y \
    autoconf \
    automake \
    bzip2 \
    lbzip2 \
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    libncurses5-dev \
    libssl-dev \
    pbzip2 \
    pigz \
    zlib1g-dev && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get autoclean


WORKDIR /opt
RUN wget --progress=dot:giga https://github.com/samtools/samtools/releases/download/${SAMTOOLSVER}/samtools-${SAMTOOLSVER}.tar.bz2 && \
    tar -xjf samtools-${SAMTOOLSVER}.tar.bz2 && \
    rm samtools-${SAMTOOLSVER}.tar.bz2
WORKDIR /opt/samtools-${SAMTOOLSVER}
RUN ./configure && \
    make && \
    make install
WORKDIR /opt

RUN wget --progress=dot:giga https://sourceforge.net/projects/bbmap/files/BBMap_${BBTOOLSVER}.tar.gz && \
    tar -xzf BBMap_${BBTOOLSVER}.tar.gz && \
    rm BBMap_${BBTOOLSVER}.tar.gz && \
    ln -s /opt/bbmap/*.sh /usr/local/bin/


RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VER}.zip && \
    unzip fastqc_v${FASTQC_VER}.zip && \
    rm fastqc_v${FASTQC_VER}.zip && \
    chmod 755 FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc


RUN apt-get update && apt-get install -y bowtie2
RUN apt-get update && apt-get install -y freebayes
RUN apt-get update && apt-get install -y bwa


RUN wget https://github.com/AstraZeneca-NGS/VarDict/archive/refs/heads/master.zip && \
   unzip master.zip && \
   chmod 755 VarDict-master/* && \
   ln -s /opt/VarDict-master/* /usr/local/bin/

RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
	tar jxf bcftools-1.9.tar.bz2 && \
	rm bcftools-1.9.tar.bz2 && \
	cd bcftools-1.9 && \
	make install


RUN wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip && \
     unzip snpEff_latest_core.zip && rm -r snpEff_latest_core.zip && \
     ln -s /opt/snpEff/  /usr/local/bin/snpEff


RUN apt-get update -y && apt-get install -y \
    libnss-sss \
    openjdk-8-jre \
    less \
    vim

RUN mkdir -p /usr/local/bin/picard \
    && curl -SL https://github.com/broadinstitute/picard/releases/download/2.21.5/picard.jar \
    > /usr/local/bin/picard/picard.jar

RUN chmod 0644 /usr/local/bin/picard/picard.jar

RUN apt-get update && apt-get install -y --no-install-recommends \
         python3.5 \
         python3-pip \
         && \
         apt-get clean && \
         rm -rf /var/lib/apt/lists/*


RUN apt-get update && \
  	apt-get install -y openjdk-11-jre && \
  	rm -rf /var/lib/apt/lists/*

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN wget https://github.com/broadinstitute/gatk/releases/download/4.1.4.1/gatk-4.1.4.1.zip \
 	&& unzip gatk-4.1.4.1.zip \
 	&& rm gatk-4.1.4.1.zip -f \
 	&& cd gatk-4.1.4.1 \
 	&& ./gatk --list


RUN pip3 install -U pip setuptools
RUN pip3 install pandas
RUN pip3 install matplotlib
RUN pip3 install pysam
RUN pip3 install xlrd==1.2.0
RUN pip3 install pyfaidx
RUN pip3 install seaborn
RUN apt-get update && apt-get install -y VarDict

ENV PATH="/opt/gatk-4.1.4.1/:${PATH}"
