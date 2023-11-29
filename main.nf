nextflow.enable.dsl=2

/*
 * Input parameters: read pairs
 * Params are stored in the params.config file
 */

 params.ref = "$baseDir/Ref/mars_pf_ref.fasta"
 params.reads = "$baseDir/Test/*{R1,R2}*.fastq.gz"
 params.adapter = "$baseDir/Ref/adapters.fa"
 params.out = "$baseDir/output"
 params.snpeff_config = "$baseDir/6Genes_ref"
 params.gff = "$baseDir/Ref/mars_pf.gff"
 //params.bed = "$baseDir/Ref/mars_pf.bed"
 params.voi = "$baseDir/Ref/voinew3.csv"
 params.memory = "8g"
 params.cpus = 1
 params.dbName = "PV_db_ref"

// import modules


include { Trim_reads; } from './modules/trim_read'
include { PreFastqC; pre_multiqc; PostFastqC; multiqc} from './modules/fastqc'
include {BWA_index} from './modules/index'
include {BWA_align} from './modules/align'
include {Sam_sort; Picard_add_read; VCF_call; Get_Bed} from './modules/vcf_call'
include {buildsnpeff_db; annotation; vartype  } from './modules/annotation'
include {vcf_to_DF; csv_merge} from './modules/csv_merge'
include {getcoverage; WT_cov ; Trim_Stats; Reads_merge} from './modules/coverage'
include {Snpfilter; Summary_merge; Summary; Dataviz_Reportable_snps; DataViz_Novel_snps; Introns_merge } from './modules/final_snp'



workflow {

    // cerate a channel for ref and reads
    ref_ch =  Channel.fromPath(params.ref, checkIfExists: true)
    read_ch = Channel.fromFilePairs(params.reads, checkIfExists: true).filter{ it.size()>0 }
    adapter_ch = Channel.fromPath(params.adapter, checkIfExists: true)
    snpeff_ch = Channel.fromPath(params.snpeff_config, checkIfExists: true)
    gff_ch = Channel.fromPath( params.gff, checkIfExists: true )
    voi_ch = Channel.fromPath( params.voi, checkIfExists: true )
  //  bed_ch = Channel.fromPath( params.bed, checkIfExists: true )

    buildsnpeff_db()
    //PreFastqC(read_ch)
    Trim_reads(adapter_ch.combine(read_ch))
    PostFastqC(Trim_reads.out.Trimmed_fastq)
    multiqc(PostFastqC.out.collect())
    BWA_index(ref_ch)
    BWA_align(BWA_index.out.collect().combine(Trim_reads.out.Trimmed_fastq))
    Sam_sort(BWA_align.out.Bwa_samfile)
    Picard_add_read(BWA_align.out.Bwa_samfile)
    Get_Bed (gff_ch)
    VCF_call(BWA_index.out.combine(Get_Bed.out).combine(Picard_add_read.out.Picard_out_bam))
    annotation(snpeff_ch.combine(VCF_call.out.variants), buildsnpeff_db.out.buildDB)
    vartype(annotation.out.var_annotation)
    getcoverage(Picard_add_read.out.Picard_out_bam)
    Trim_Stats(Trim_reads.out.Trimmed_stats.join(getcoverage.out.Samtools_coverage))
    WT_cov(ref_ch.combine(gff_ch).combine(voi_ch).combine(getcoverage.out.samtools_depth))
    Reads_merge(Trim_Stats.out.Reads_cov.collect())
    vcf_to_DF(vartype.out.vartype_annotation)
    csv_merge(vcf_to_DF.out.csv_annotate)
    Snpfilter(voi_ch.combine(csv_merge.out.CSV_merge.join(WT_cov.out.WT_coverage)))
    Summary_merge(Snpfilter.out.snp_report.collect())
    Summary(Summary_merge.out)
    Introns_merge(csv_merge.out.Introns_file.collect())
    Dataviz_Reportable_snps(Summary.out.Reportable_snps)
    DataViz_Novel_snps(Summary.out.Novel_snps)



  }
