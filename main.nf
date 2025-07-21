nextflow.enable.dsl=2

//params.help = false
params.threads = 1
params.output_dir = './results'
params.DEG = true
params.strandness = 1
params.pairend = 1

//print usage
if (params.help) {
    log.info ''
    log.info 'Pipeline to call DAP-seq peaks with MACS3 and peak annotation and motif analyses'
    log.info '--------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run ./ -params-file params.yml'
    log.info ''
    log.info 'Options:'
    log.info '	--help	Show this message and exit.'
    log.info '  --bowtie_idx bowtie_idx    bowtie index prefix for analysis.'
    log.info '	--fq_sheet fq_sheet    A tab-delimited file storing the samples information, with three columns: sample_id,fq1,fq2,single_end,control.'
    log.info '	--fasta    Genome fasta file for the analyzing species.'
    log.info '	--gtf    Genome gtf file for the analyzing species.'
    log.info '  --output_dir OUTDIR   Name for directory for saving the results. Default: results/'
    log.info '  --fq_dir raw_data    The folder where the raw .fq files are.'
    log.info '  --DEG    Whether to perform DEG anaysis. Default: false'
    log.info '  --DEG_design    Path to the file of design for group comparison in DEG analysis.'
    log.info '  --strandness STRAND   Library strandness for featureCounts: 0 (unstranded), 1 (stranded), 2 (reverse stranded). Default: 1'
    log.info '  --pairend    1 (pair-end), 0 (single-end). Default: 1'
    exit 1
}

include {INPUT_CHECK} from "./subworkflow/input_check"
include {TRIMGALORE} from "./module/Trimgalore"
include {FASTQC} from "./module/Fastqc"
include {BAM2BW} from "./module/Bam2bw"
include {BOWTIE2MAP} from "./module/Bowtie2Mapping"
include {MARK_DUPLICATES} from './module/markduplicates'
include {FEATURECOUNTS} from './module/featurecounts'
include {REPORT} from './module/Report'
include {DESEQ} from './module/DESEQ2'


// Define a process to collect completion signals
process COMPLETION_CHECK {
    input:
        val(count)
        val(output_dir)
    
    output:
        val(true), emit: ready
    
    exec:
        println "All $count samples have completed processing"
}


// SIMPLIFIED AND CORRECTED APPROACH:

workflow {
    INPUT_CHECK(params.fq_sheet)
    
    INPUT_CHECK.out.reads.count().view { "=== TOTAL INPUT SAMPLES: $it ===" }
    
    TRIMGALORE(INPUT_CHECK.out.reads)
    BOWTIE2MAP(TRIMGALORE.out.reads)

    FASTQC(TRIMGALORE.out.reads)
    //MARK_DUPLICATES(BOWTIE2MAP.out.bam)
    
    //MARK_DUPLICATES
    BOWTIE2MAP
        .out
        .bam
        .join(BOWTIE2MAP.out.bai, by: [0])
        .set {ch_genome_bam_bai}

    ch_genome_bam_bai.count().view { "=== TOTAL BAM SAMPLES: $it ===" }

    BAM2BW(ch_genome_bam_bai)

    BAM2BW.out.count().set{ sample_count }

    COMPLETION_CHECK(sample_count, params.output_dir)

    // Collect all BAM files from BOWTIE2MAP and pass them together to FEATURECOUNTS
    BOWTIE2MAP.out.bam
        .map { meta, bam -> bam }  // Extract only the BAM file path, discard metadata
        .collect()
        .set { ch_collected_bams }


    // Debug: view the collected BAM files
    ch_collected_bams.view { "=== COLLECTED BAM FILES: $it ===" }

    FEATURECOUNTS(ch_collected_bams, params.gtf, params.strandness, params.pairend)


    REPORT(params.fq_sheet, params.output_dir, COMPLETION_CHECK.out.ready, FEATURECOUNTS.out.counts)

    REPORT.out.counts.view { "=== report out counts: $it ==="}

    if (params.DEG) {
        DESEQ(REPORT.out.counts)
    }

}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}