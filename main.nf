nextflow.enable.dsl=2

/*
 * Input parameters: read pairs
 * Params are stored in the params.config file
 */



// import modules

include { PreFastqC} from './modules/fastqc'
include { Trim_reads } from './modules/trim_read'
include { PostFastqC} from './modules/fastqc'
include {BWA_index} from './modules/index'
include {BWA_align} from './modules/align'
include {Sam_sort; Picard_add_read; VCF_call} from './modules/vcf_call'
include {annotation; vartype  } from './modules/annotation'
include {vcf_to_DF; csv_merge} from './modules/csv_merge'
include {getcoverage; WT_cov  } from './modules/coverage'
include {Snpfilter} from './modules/final_snp'



workflow {

    // cerate a channel for ref and reads
    ref_ch =  Channel.fromPath(params.ref, checkIfExists: true)
    read_ch = Channel.fromFilePairs(params.reads, checkIfExists: true, size: -1)
    adapter_ch = Channel.fromPath(params.adapter, checkIfExists: true)
    snpeff_ch = Channel.fromPath(params.snpeff_config, checkIfExists: true)
    // pyscripts = Channel.fromPath(params.pyscripts, checkIfExists: true)

    PreFastqC(read_ch)
    Trim_reads(adapter_ch.combine(read_ch))
    PostFastqC(Trim_reads.out.Trimmed_fastq)
    BWA_index(ref_ch)
    BWA_align(BWA_index.out.collect().combine(Trim_reads.out.Trimmed_fastq))
    Sam_sort(BWA_align.out.Bwa_samfile)
    Picard_add_read(BWA_align.out.Bwa_samfile)
    VCF_call(BWA_index.out.combine(Picard_add_read.out.Picard_out_bam))
    annotation(snpeff_ch.combine(VCF_call.out.variants))
    vartype(annotation.out.var_annotation)
    vcf_to_DF(vartype.out.vartype_annotation)
    csv_merge(vcf_to_DF.out.csv_annotate)
    getcoverage(Picard_add_read.out.Picard_out_bam)
    WT_cov(getcoverage.out.samtools_depth)
    Snpfilter(csv_merge.out.CSV_merge.join(WT_cov.out.WT_coverage))

  }
