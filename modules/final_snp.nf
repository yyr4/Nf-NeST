process Snpfilter {
    publishDir "${params.out}/Snpfilter/${sample_id}", mode : "copy"
    tag  { "Snpfilter ${sample_id }"}

    input:

      tuple path(voi), val(sample_id), path ("${sample_id}_merge.csv"), path("${sample_id}_coverage.csv")


    output:

      tuple val(sample_id), path("${sample_id}_final_snp.csv"),   emit: snp_report

    script:


       """
       final_snpfilter.py  -A ${sample_id}_merge.csv -C ${sample_id}_coverage.csv -V ${voi} > ${sample_id}_final_snp.csv


       """


}

process Summary_merge{
    publishDir "${params.out}/Summary_merge/", mode : "copy"


    input:

      file("*")

    output:

      file ("merged_final_snp.csv")


     script:


       """
       awk 'FNR==1 && NR!=1{next;}{print}' *.csv > merged_final_snp.csv
       """
}


process Summary {
    publishDir "${params.out}/Summary/", mode : "copy"


    input:

      path "merged_final_snp.csv"


    output:

     path "All_final_snp.csv", emit: Summary_Snp
     path "Reportable_snps.csv", emit: Reportable_snps
     path "Novel_snps.csv", emit: Novel_snps

    script:


       """

       summary.py -f merged_final_snp.csv

       """


}


process Dataviz_Reportable_snps {
    publishDir "${params.out}/Summary/", mode : "copy"


    input:

      path x


    output:

      file ('Reportable_Per_SNP_depth.pdf')
      file ('SNPs-Reportable.pdf')


    script:


       """
        Dataviz_Reportable_snps.py -r $x

       """

}
process DataViz_Novel_snps {
    publishDir "${params.out}/Summary/", mode : "copy"


    input:
      path y

    output:

      file ('SNPs-Novel-missense.pdf')
      file ('SNPs-Novel-synonymous.pdf')

    script:


       """
        Dataviz_Novel_snps.py -n $y
       """

}
