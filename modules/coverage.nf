process getcoverage {

    publishDir "${params.out}/samtoolcoverage/${sample_id}", mode : "copy"

    input:
     tuple val (sample_id) ,path("${sample_id}_picard_readgroup.bam")

    output:
    tuple val(sample_id), path("${sample_id}_coverage.txt" )
    tuple val(sample_id), path("${sample_id}_finalcoverageall.txt")
    tuple val(sample_id), path("${sample_id}_depth.txt") ,      emit: samtools_depth


    script:
    """
    samtools index ${sample_id}_picard_readgroup.bam
    samtools coverage ${sample_id}_picard_readgroup.bam -o ${sample_id}_coverage.txt
    samtools depth ${sample_id}_picard_readgroup.bam -a -o ${sample_id}_depth.txt

    cat *coverage.txt > ${sample_id}_finalcoverageall.txt

    """
}


process WT_cov {
    publishDir "${params.out}/WT_cov/${sample_id}", mode : "copy"
    tag  { "WT_cov ${sample_id }"}

    input:

      tuple path(ref),path(bed), path(voi), val(sample_id), path("${sample_id}_depth.txt")



    output:
      tuple val(sample_id), path("${sample_id}_coverage.csv"), emit: WT_coverage

    script:

        """
        Wt_cov.py -R ${ref} -B ${bed} -N '${sample_id}_depth.txt' -V ${voi} > ${sample_id}_coverage.csv

        """
}
