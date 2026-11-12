process BWA_MEM2_ALIGN {
    tag "${meta.sample_id}"
    label 'process_high'

    input:
    tuple val(meta), path(trimmed_r1), path(trimmed_r2)
    path fasta_ref

    output:
    tuple val(meta), path("${meta.sample_id}.sorted.bam"), emit: aligned_bam

    script:
    def rg = "@RG\\tID:${meta.sample_id}\\tSM:${meta.sample_id}\\tPL:ILLUMINA\\tLB:${meta.batch_id}\\tPU:${meta.assay}"
    """
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R '${rg}' \\
        ${fasta_ref} \\
        ${trimmed_r1} ${trimmed_r2} \\
        | samtools sort -@ 4 -o ${meta.sample_id}.sorted.bam -
    """
}

process SAMTOOLS_MARKDUP_STATS {
    tag "${meta.sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/alignments", mode: 'copy', pattern: "*.flagstat.txt"

    input:
    tuple val(meta), path(sorted_bam)
    path target_bed

    output:
    tuple val(meta), path("${meta.sample_id}.dedup.bam"), path("${meta.sample_id}.dedup.bam.bai"), emit: dedup_bam
    tuple val(meta), path("${meta.sample_id}.flagstat.txt"), emit: flagstat_log

    script:
    """
    samtools markdup -@ ${task.cpus} -s ${sorted_bam} ${meta.sample_id}.dedup.bam
    samtools index -@ ${task.cpus} ${meta.sample_id}.dedup.bam
    samtools flagstat -@ ${task.cpus} ${meta.sample_id}.dedup.bam > ${meta.sample_id}.flagstat.txt
    """
}
