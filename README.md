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
chmod +x nextflow

```
3. If you don't have it already install Docker in your computer, download and install from [here](https://docs.docker.com/).  

     - Docker Nf-NeST comes with a Doker image that can be run pipeline with virtual environment. To setup up Docker image, run the following command from the Nf-NeST directory.

```
cd Nf-NeST
docker pull dpatel89/nfnest_image:latest
 
```
4. Launch the pipeline execution:

```
nextflow run main.nf -with-docker dpatel89/nfnest_image

```

# Input

- Input raw data are fastq Paired reads separated by R1 and R2. Example test data are given in the folder test. User can run their own analysis by Specifies the location of the  reads FASTQ file (--reads option)
  
```
nextflow run main.nf --reads '/home/dataset/*{R1,R2}*.fastq.gz' -with-docker dpatel89/nfnest_image

```
- Referance fasta file of targeted gene amplicones.(mars_pf_ref.fasta)
- Bed file which is a tab delimited text file contains genomics coordinates and associated annotations (mars_pf.bed)
- Variant of Interest file of WHO reported known snps(voinew3.csv)


# Output

- Output folder will be created under Nf-NeST. The results can be found under **output/Snpfilter/** folder. Separate folder for each samples will be created with reportable and novel mutations csv file.  

