process PreFastqC {

    tag { "PreFastqC ${reads}"}

    publishDir "${params.out}/PrefastqC/HTML", pattern: "*.html",  mode:'copy'
    publishDir "${params.out}/PrefastqC/ZIP", pattern: "*.zip",  mode:'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*")

    script:
    """
    fastqc ${reads}

    """

}

process PostFastqC{

    tag { "PostFastqC ${ref}"}

    publishDir "${params.out}/PostfastqC/HTML", pattern: "*.html",  mode:'copy'
    publishDir "${params.out}/PostfastqC/ZIP", pattern: "*.zip",  mode:'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*")

    script:
    """
    echo ${sample_id}
    fastqc ${reads}

    """

}
