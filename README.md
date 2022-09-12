# Nextflow Next-generation Sequence-analysis Toolkit (Nf-NeST) : A standardized bioinformatics framework for analyzing SNPs in next-generation sequencing data


# Requirements
- Unix-like operating system (Linux, macOS, etc)
- Java 8
- Docker
- Nextflow

# Installation

1. Download git repository:

```
git clone https://github.com/yyr4/Nf-NeST.git

```
2. Install Nextflow (version 20.07.x or higher):

```
 curl -s https://get.nextflow.io | bash

```
3. If you don't have it already install Docker in your computer, download and install from [here](https://docs.docker.com/).  

4. Launch the pipeline execution:
```
 ./nextflow run yyr4/Nf-NeST -with-docker
```