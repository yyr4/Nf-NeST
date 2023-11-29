process getcoverage {

    publishDir "${params.out}/samtoolcoverage/${sample_id}", mode : "copy"

    input:
     tuple val (sample_id) ,path("${sample_id}_picard_readgroup.bam")

    output:
    tuple val(sample_id), path("${sample_id}_coverage.txt" ), emit: Samtools_coverage
    tuple val(sample_id), path("${sample_id}_finalcoverageall.txt")
    tuple val(sample_id), path("${sample_id}_depth.txt") ,      emit: samtools_depth


    script:
    """
    samtools index ${sample_id}_picard_readgroup.bam
    samtools coverage ${sample_id}_picard_readgroup.bam -o ${sample_id}_coverage.txt
    samtools depth ${sample_id}_picard_readgroup.bam -aa -o ${sample_id}_depth.txt

    cat *coverage.txt > ${sample_id}_finalcoverageall.txt

    """
}


process WT_cov {
    publishDir "${params.out}/WT_cov/${sample_id}", mode : "copy"
    tag  { "WT_cov ${sample_id }"}
    errorStrategy 'ignore'
    input:

      tuple path(ref),path(gff), path(voi), val(sample_id), path("${sample_id}_depth.txt")



    output:
      tuple val(sample_id), path("${sample_id}_coverage.csv"), emit: WT_coverage

    script:

        """
        Wt_cov.py -R ${ref} -G ${gff} -N '${sample_id}_depth.txt' -V ${voi} > ${sample_id}_coverage.csv

        """
}


process Trim_Stats {

      // publishDir "${params.out}/Reads", pattern: "*.csv", mode : "copy"

      input:
      tuple val(sample_id), path("${sample_id}.stats.txt"), path("${sample_id}_coverage.txt")


      output:

      tuple val(sample_id), path("${sample_id}_readcoverage.csv"), emit: Reads_cov

      script:
      """
      Reads_stats.py  -S "${sample_id}.stats.txt" -N "${sample_id}_coverage.txt" > "${sample_id}_readcoverage.csv"

      """




}

process Reads_merge{
    publishDir "${params.out}/Summary/", mode : "copy"


    input:

      file("*")

    output:

      file ("Reads_Metrics_Samples.csv")


     script:


       """
       awk 'FNR==1 && NR!=1{next;}{print}' *.csv > Reads_Metrics_Samples.csv
       """
}
