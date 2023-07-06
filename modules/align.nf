process BWA_align{

    tag { "BWA_align ${sample_id}"}
    // publishDir "${params.out}/bwa_align", pattern: "*.sam",  mode:'copy'

    // input is reads and refenace
    input:

    tuple path(ref), path("*"), val(sample_id), path(reads)


    // output is samfile
    output:
    tuple val (sample_id) ,path("${sample_id}_align.sam"), emit: Bwa_samfile


    //script //  bwa mem ref.fa read1.fq read2.fq > aln-pe.sam
    script:
    """
    bwa mem  ${ref} ${reads}  > ${sample_id}_align.sam
    """

}
