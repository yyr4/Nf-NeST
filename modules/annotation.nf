process annotation {
      tag  { "annotation ${sample_id }"}
      publishDir "${params.out}/snpEff_annotation", mode:'copy'


      input:

      tuple path(snpeff_config), val (sample_id), path("${sample_id}_samtools.vcf"), path("${sample_id}_Freebayes.vcf"), path("${sample_id}_gatk.vcf")


      output:
      tuple val(sample_id), path("${sample_id}_samtools_ann.vcf"), path("${sample_id}_Freebayes_ann.vcf"), path("${sample_id}_gatk_ann.vcf"),  emit: var_annotation



      script:
      """
      mkdir /usr/local/bin/snpEff/data
      cp -r ${snpeff_config} /usr/local/bin/snpEff/data

      java -jar /usr/local/bin/snpEff/snpEff.jar -c /usr/local/bin/snpEff/snpEff.config -hgvs1LetterAa -noShiftHgvs ${snpeff_config} ${sample_id}_Freebayes.vcf > ${sample_id}_Freebayes_ann.vcf
      java -jar /usr/local/bin/snpEff/snpEff.jar -c /usr/local/bin/snpEff/snpEff.config -hgvs1LetterAa -noShiftHgvs ${snpeff_config} ${sample_id}_gatk.vcf > ${sample_id}_gatk_ann.vcf
      java -jar /usr/local/bin/snpEff/snpEff.jar -c /usr/local/bin/snpEff/snpEff.config -hgvs1LetterAa -noShiftHgvs ${snpeff_config} ${sample_id}_samtools.vcf > ${sample_id}_samtools_ann.vcf



      """
}



process vartype {
    publishDir "${params.out}/vartype", mode : "copy"

    input:
        tuple val(sample_id), path("${sample_id}_samtools_ann.vcf"), path("${sample_id}_Freebayes_ann.vcf"), path("${sample_id}_gatk_ann.vcf")

    output:
        tuple val(sample_id), path("${sample_id}_samtools_vartype.vcf"), path("${sample_id}_freeBayes_vartype.vcf"), path("${sample_id}_HaplotypeCaller_vartype.vcf"), emit: vartype_annotation

    script:
        """
        java -jar /usr/local/bin/snpEff/SnpSift.jar varType ${sample_id}_samtools_ann.vcf > ${sample_id}_samtools_vartype.vcf
        java -jar /usr/local/bin/snpEff/SnpSift.jar varType ${sample_id}_Freebayes_ann.vcf > ${sample_id}_freeBayes_vartype.vcf
        java -jar /usr/local/bin/snpEff/SnpSift.jar varType ${sample_id}_gatk_ann.vcf > ${sample_id}_HaplotypeCaller_vartype.vcf
        """
}
