process BWA_index {

    tag { "BWA_index ${ref}"}
    label 'process_low'

    publishDir "${params.out}/bwa_index", mode:'copy'

    input:
    path ref

    output:

    tuple path(ref), path("*"),    emit: Index_file

    script:
    """
    bwa index ${ref} > ref.index
    java -jar /usr/local/bin/picard/picard.jar CreateSequenceDictionary R=${ref}
    samtools faidx ${ref}

    """

}
