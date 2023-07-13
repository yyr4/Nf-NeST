// Adapter and quality trimming
process Trim_reads {
    // errorStrategy 'ignore'

    tag { "Trim_reads${sample_id}"}

    // publishDir "${params.out}/Trimmed_Fastq/${sample_id}", pattern: "*.fastq", mode : "copy"
    publishDir "${params.out}/Trimmed_Fastq/${sample_id}/Stats", pattern: "*.txt", mode : "copy"

    input:
    tuple path (adapter), val(sample_id), path(reads)



    output:
    tuple val(sample_id), path("*.bbduk.fastq"), emit: Trimmed_fastq
    // file "${sample_id}.stats.txt"
    tuple val(sample_id), path("${sample_id}.stats.txt"), emit: Trimmed_stats


    script:


     """

     bbduk.sh  -Xmx1g  ktrimright=t k=27 hdist=1 edist=0 ref=${adapter} \\
     qtrim=rl trimq=35 minlength=100 trimbyoverlap=t minoverlap=24 qin=33 in=${reads[0]} in2=${reads[1]} \\
     out=${sample_id}.R1.bbduk.fastq out2=${sample_id}.R2.bbduk.fastq stats=${sample_id}.stats.txt

    """


}
