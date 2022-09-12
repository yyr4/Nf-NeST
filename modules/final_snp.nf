process Snpfilter {
    publishDir "${params.out}/Snpfilter/${sample_id}", mode : "copy"
    tag  { "Snpfilter ${sample_id }"}

    input:
      tuple  val(sample_id), path ("${sample_id}_merge.csv"), path("${sample_id}_coverage.csv")


    output:

      tuple val(sample_id), path("${sample_id}_final_snp.csv"), path("${sample_id}_Reportable.csv"), path("${sample_id}_Novel.csv"),   emit: snp_report

    script:


       """
       final_snpfilter.py  -A ${sample_id}_merge.csv -C ${sample_id}_coverage.csv -V $params.voi > ${sample_id}_final_snp.csv


       """


}
