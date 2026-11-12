process FASTP_TRIM {
    tag "${meta.sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: "*.json"

    input:
    tuple val(meta), path(reads_r1), path(reads_r2)

    output:
    tuple val(meta), path("${meta.sample_id}_R1.trimmed.fq.gz"), path("${meta.sample_id}_R2.trimmed.fq.gz"), emit: trimmed_reads
    tuple val(meta), path("${meta.sample_id}.fastp.json"), emit: fastp_json
    tuple val(meta), path("${meta.sample_id}.fastp.html"), emit: fastp_html

    script:
    """
    fastp \\
        --in1 ${reads_r1} \\
        --in2 ${reads_r2} \\
        --out1 ${meta.sample_id}_R1.trimmed.fq.gz \\
        --out2 ${meta.sample_id}_R2.trimmed.fq.gz \\
        --detect_adapter_for_pe \\
        --trim_poly_g \\
        --qualified_quality_phred 20 \\
        --length_required 50 \\
        --thread ${task.cpus} \\
        --json ${meta.sample_id}.fastp.json \\
        --html ${meta.sample_id}.fastp.html
    """
}
