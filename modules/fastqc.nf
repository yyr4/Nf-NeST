process PreFastqC {

    tag { "PreFastqC ${reads}"}


    publishDir "${params.out}/PrefastqC/ZIP", pattern: "*.zip",  mode:'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path ("*_fastqc.zip")

    script:
    """
    fastqc ${reads}

    """

}

process pre_multiqc{


     publishDir "${params.out}/multiqc1/", pattern: "*.html",  mode:'copy'


     input:
     path ("*")

     output:
     path("multiqc_report.html")

     script:
     """

      multiqc .

     """

}

process PostFastqC{

    //tag { "PostFastqC ${reads}"}

    publishDir "${params.out}/PostfastqC/",  mode:'copy'


    input:
    tuple val(sample_id), path(reads)

    output:
    path ("*_fastqc.zip")

    script:
    """

    fastqc ${reads}

    """

}

process multiqc{


     publishDir "${params.out}/multiqc/", pattern: "*.html",  mode:'copy'


     input:
     path ("*")

     output:
     path("multiqc_report.html")

     script:
     """

      multiqc .

     """

}
