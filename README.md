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

- If the samples are Plasmodium falciparum, run the below command. Pfalciparum reference files are provided in Pfalciparum_MaRS folder.

```

nextflow run main.nf --reads '/home/dataset/*{R1,R2}*.fastq.gz' -profile docker

```


  - Referance fasta file of targeted gene amplicones.(mars_pf_ref.fasta)
  - gff file is a tab deliminated text file which contains annotation feature describing the gene and protein sequences (mars_pf.gff)
  - Variant of Interest file from WHO reported known snps(voinew3.csv)
  - SNPEff custom annotation database files in genebank format(pf_3D7_snpEff_db/genes.gbk)
  - SNPEff's config file (pf_3D7_snpEff_db/snpEff.config)


* To find the drug resisatnce mutation of Plasmodium Vivax specimens, use the below command. Pvivax reference files are provided in Pvivax_MaRS folder.

```
nextflow run main.nf --reads "/home/dataset/*{R1,R2}*.fastq.gz" --ref Pvivax_MaRS/pv_sal_Ref/PV-salI.fasta --gff Pvivax_MaRS/pv_sal_Ref/PV_sal1.gff --voi Pvivax_MaRS/pv_sal_Ref/voinew3.csv --dbName pv_sal_snpEff_db --out PV_out -profile docker

```

**Note**:
- The current workflow was designed for Plasmodium falciparum. To run the PVivax samples, user need to specify the snpeff config file in the nextflow.config file (params.snpeff_config = "$baseDir/Pvivax_MaRS/pv_sal_snpEff_db"). Also, Uncomment the line 5 (path "pv_sal_snpEff_db", type: 'dir', emit: buildDB) in annotation module file from path "$baseDirT/modules/annotation.nf" 


# Output

- Output folder will be created under Nf-NeST. The results can be found under **output/Summary/** folder.

* Nf-NeST produces various table reports and summarization figures which are described below.

    - All_final_snp.csv : It reports all the Snps (Known + Novel) found in the sample
 
    - DMS_EPI_report.csv : This is a report that is used for reporting the Domestic imported malaria cases.

    - Reportable_snps.csv: It reports all the WHO reported and other known mutations which are listed in voinew3.csv file, (minimum coverage is >= 5 and min Allele frequency >=5%)
 
    - Introns_final_snp.csv: It reports all  intronic variants found in all the samples in the study
 
    - Novel_snps.csv: As name suggested, it reports all the novel mutations which are not listed in voinew3.csv file and have minimum confidence 2 or more (i.e, that snp is reported by 2 or more variant caller), minimum coverage is >= 5 and min Allele frequency >=5%

    - Reads_Metrics_Samples.csv:  It reports total Reads aligned per gene per sample

    - Reportable_Per_SNP_depth.pdf: Read depth of coverage for single nucleotide polymorphisms (SNPs) associated with malaria drug resistance

    - SNPs-Reportable.pdf: Bar graph depicting the wild type, major and minor allele frequencies of associated and/or confirmed resistance SNPs . Allele frequencies are indicated on the x axis, and the variants of interest are listed along the y-axis. The color coding indicates the type of mutation found in the samples; blue is for wild type, green for minor allele mutation and red for major allele mutation.

    - SNPs-Novel-missense.pdf and SNPs-Novel-synonymous.pdf: Bar graph depicting the wild type, major and minor allele frequencies of Novel missense and non-Synonymous mutations.
