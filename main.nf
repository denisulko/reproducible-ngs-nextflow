#!/usr/bin/env nextflow
/*
 * ============================================================================
 * Reproducible Germline NGS Variant Calling & MultiQC Telemetry Pipeline
 * Nextflow DSL2 End-to-End Orchestrator
 * ============================================================================
 */

nextflow.enable.dsl = 2

include { FASTP_TRIM               } from './modules/local/fastp_trim.nf'
include { BWA_MEM2_ALIGN           } from './modules/local/bwa_align.nf'
include { SAMTOOLS_MARKDUP_STATS   } from './modules/local/bwa_align.nf'
include { BCFTOOLS_CALL_FILTER     } from './modules/local/variant_call.nf'
include { MULTIQC_COHORT_AGGREGATE } from './modules/local/multiqc_aggregate.nf'

workflow GERMLINE_NGS_PIPELINE {
    take:
    ch_samplesheet // channel: [ val(meta), path(reads_r1), path(reads_r2) ]
    ch_fasta       // path: GRCh38 reference FASTA
    ch_targets     // path: Exome/WGS target regions BED

    main:
    // Step 1: Adapter trimming, poly-G removal, and Phred Q30 pre/post QC
    FASTP_TRIM(ch_samplesheet)

    // Step 2: BWA-MEM2 coordinate-sorted read alignment against GRCh38
    BWA_MEM2_ALIGN(FASTP_TRIM.out.trimmed_reads, ch_fasta)

    // Step 3: Duplicate marking, coverage depth profiling, and samtools stats
    SAMTOOLS_MARKDUP_STATS(BWA_MEM2_ALIGN.out.aligned_bam, ch_targets)

    // Step 4: Germline SNV/Indel calling & hard quality filtering (QUAL>=30, DP>=15)
    BCFTOOLS_CALL_FILTER(SAMTOOLS_MARKDUP_STATS.out.dedup_bam, ch_fasta, ch_targets)

    // Step 5: Cohort-level MultiQC & SLA telemetry aggregation
    ch_all_qc = FASTP_TRIM.out.fastp_json
        .map { _meta, json -> json }
        .mix(SAMTOOLS_MARKDUP_STATS.out.flagstat_log.map { _meta, log -> log })
        .mix(BCFTOOLS_CALL_FILTER.out.vcf_stats.map { _meta, stats -> stats })
        .collect()

    MULTIQC_COHORT_AGGREGATE(ch_all_qc)

    emit:
    filtered_vcfs = BCFTOOLS_CALL_FILTER.out.filtered_vcf
    cohort_report = MULTIQC_COHORT_AGGREGATE.out.summary_csv
}

workflow {
    Channel
        .fromPath(params.input_samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                sample_id : row.sample_id,
                batch_id  : row.batch_id,
                assay     : row.assay_type,
                sex       : row.sex
            ]
            tuple(meta, file(row.fastq_r1), file(row.fastq_r2))
        }
        .set { ch_input_samples }

    GERMLINE_NGS_PIPELINE(
        ch_input_samples,
        file(params.fasta_ref),
        file(params.target_bed)
    )
}
