process MULTIQC_COHORT_AGGREGATE {
    tag "cohort_multiqc"
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path qc_files

    output:
    path "multiqc_cohort_summary.csv", emit: summary_csv
    path "multiqc_report.json", emit: report_json
    path "pipeline_execution_report.md", emit: report_md

    script:
    """
    python ${projectDir}/scripts/run_multiqc_aggregator.py \\
        --output-dir .
    """
}
