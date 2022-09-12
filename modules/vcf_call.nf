process Sam_sort {

    tag { "Sam_sort ${sample_id}"}
    publishDir "${params.out}/bam_out", mode:'copy'

    // input is reads and refenace
    input:

    tuple val(sample_id), path(sam)


    // output is samfile
    output:
    tuple val (sample_id) ,path("${sample_id}_sorted.bam"), emit: Bam_file


    //script //  bwa mem ref.fa read1.fq read2.fq > aln-pe.sam
    script:
    """

    samtools view -S -b ${sample_id}_align.sam > ${sample_id}_align.bam
    samtools sort -o "${sample_id}_sorted.bam" ${sample_id}_align.bam
    """

}


process Picard_add_read{

    tag { "Picard_add_read${sample_id}"}
    publishDir "${params.out}/picard_out", mode:'copy'


    input:
    tuple val (sample_id) ,path("${sample_id}_align.sam")

    output:
    tuple val (sample_id) ,path("${sample_id}_picard_readgroup.bam"), emit: Picard_out_bam



    script:

    """

    java -jar /usr/local/bin/picard/picard.jar AddOrReplaceReadGroups  I=${sample_id}_align.sam O=${sample_id}_picard_readgroup.bam \\
    SORT_ORDER=coordinate RGLB=ExomeSeq  RGPL=Illumina RGPU=Miseq  RGSM=${sample_id} CREATE_INDEX=True


    """

}

process VCF_call {
  tag  { "VCF_call ${sample_id }"}
  publishDir "${params.out}/vcf_files", mode:'copy'
  cpus params.cpus
  memory params.memory


  input:
    tuple path(ref), path("*"), val (sample_id), path("${sample_id}_picard_readgroup.bam")

  output:
    //tuple val (sample_id) ,path("${sample_id}.mpileup"),      emit: BCF_out
    //tuple val (sample_id) ,path("${sample_id}_bcf.vcf"),      emit: Variants_out
    tuple val (sample_id) ,path("${sample_id}_samtools.vcf"),path("${sample_id}_Freebayes.vcf"),  path("${sample_id}_gatk.vcf") , emit: variants



  script:

    """
    bcftools mpileup -O b -o ${sample_id}.mpileup -f ${ref} ${sample_id}_picard_readgroup.bam
    bcftools call --ploidy 1 -m -v -o ${sample_id}_bcf.vcf ${sample_id}.mpileup
    vcfutils.pl varFilter ${sample_id}_bcf.vcf > ${sample_id}_samtools.vcf
    freebayes -f  ${ref} -F 0.01 -E 3 --report-all-haplotype-alleles --haplotype-length 1 ${sample_id}_picard_readgroup.bam > ${sample_id}_Freebayes.vcf
    gatk HaplotypeCaller --native-pair-hmm-threads 8 -R ${ref} -I ${sample_id}_picard_readgroup.bam  --min-base-quality-score 0 -O ${sample_id}_gatk.vcf
    """

}
