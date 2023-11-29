process buildsnpeff_db{
      output:

      path "pf_3D7_snpEff_db", type: 'dir', emit: buildDB
      //path "pv_sal_snpEff_db", type: 'dir', emit: buildDB

      script:
      """

      mkdir -p ${params.dbName}/${params.dbName}
      cp ${params.snpeff_config}/snpEff.config ${params.dbName}
      cp ${params.snpeff_config}/genes.gbk ${params.dbName}/${params.dbName}


      java -jar /usr/local/bin/snpEff/snpEff.jar build -c ${params.dbName}/snpEff.config -genbank -v ${params.dbName}

      """
}
process annotation {
      tag  { "annotation ${sample_id }"}
      // publishDir "${params.out}/snpEff_annotation", mode:'copy'


      input:

      tuple path(snpeff_config), val (sample_id), path("${sample_id}_samtools.vcf"), path("${sample_id}_Freebayes.vcf"), path("${sample_id}_gatk.vcf"), path("${sample_id}_vardict.vcf")
      each(dbDir)

      output:
      tuple val(sample_id), path("${sample_id}_samtools_ann.vcf"), path("${sample_id}_Freebayes_ann.vcf"), path("${sample_id}_gatk_ann.vcf"), path("${sample_id}_vardict_ann.vcf"), emit: var_annotation



      script:
      """

      java -jar /usr/local/bin/snpEff/snpEff.jar -c ${dbDir}/snpEff.config -hgvs1LetterAa -noShiftHgvs ${params.dbName} ${sample_id}_Freebayes.vcf > ${sample_id}_Freebayes_ann.vcf
      java -jar /usr/local/bin/snpEff/snpEff.jar -c ${dbDir}/snpEff.config -hgvs1LetterAa -noShiftHgvs ${params.dbName} ${sample_id}_gatk.vcf > ${sample_id}_gatk_ann.vcf
      java -jar /usr/local/bin/snpEff/snpEff.jar -c ${dbDir}/snpEff.config -hgvs1LetterAa -noShiftHgvs ${params.dbName} ${sample_id}_samtools.vcf > ${sample_id}_samtools_ann.vcf
      java -jar /usr/local/bin/snpEff/snpEff.jar -c ${dbDir}/snpEff.config -hgvs1LetterAa -noShiftHgvs ${params.dbName} ${sample_id}_vardict.vcf > ${sample_id}_vardict_ann.vcf



      """
}



process vartype {
    publishDir "${params.out}/vartype", mode : "copy"

    input:
        tuple val(sample_id), path("${sample_id}_samtools_ann.vcf"), path("${sample_id}_Freebayes_ann.vcf"), path("${sample_id}_gatk_ann.vcf"), path("${sample_id}_vardict_ann.vcf")

    output:
        tuple val(sample_id), path("${sample_id}_samtools_vartype.vcf"), path("${sample_id}_freeBayes_vartype.vcf"), path("${sample_id}_HaplotypeCaller_vartype.vcf"), path("${sample_id}_vardict_vartype.vcf"), emit: vartype_annotation

    script:
        """
        java -jar /usr/local/bin/snpEff/SnpSift.jar varType ${sample_id}_samtools_ann.vcf > ${sample_id}_samtools_vartype.vcf
        java -jar /usr/local/bin/snpEff/SnpSift.jar varType ${sample_id}_Freebayes_ann.vcf > ${sample_id}_freeBayes_vartype.vcf
        java -jar /usr/local/bin/snpEff/SnpSift.jar varType ${sample_id}_gatk_ann.vcf > ${sample_id}_HaplotypeCaller_vartype.vcf
        java -jar /usr/local/bin/snpEff/SnpSift.jar varType ${sample_id}_vardict_ann.vcf > ${sample_id}_vardict_vartype.vcf
        """
}
