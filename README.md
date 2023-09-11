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
cd Nf-NeST
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin

```
3. If you don't have it already install Docker in your computer, download and install from [here](https://docs.docker.com/).  



# Input

- Input raw data are fastq Paired reads separated by R1 and R2. Example test data are given in the folder test. User can run their own analysis by Specifies the location of the  reads FASTQ file (--reads option)

```

nextflow run main.nf --reads '/home/dataset/*{R1,R2}*.fastq.gz' -profile docker

```
- Referance fasta file of targeted gene amplicones.(mars_pf_ref.fasta)
- Bed file is a tab delimited text file which contains genomics coordinates and associated annotations (mars_pf.bed)
- Variant of Interest file from WHO reported known snps(voinew3.csv)
- SNPEff custom annotation database files in genebank format(6Genes_ref/genes.gbk)
- SNPEff's config file (6Genes_ref/snpEff.config)


# Output

- Output folder will be created under Nf-NeST. The results can be found under **output/Summary/** folder.
