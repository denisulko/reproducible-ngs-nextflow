process BCFTOOLS_CALL_FILTER {
    tag "${meta.sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple val(meta), path(dedup_bam), path(dedup_bai)
    path fasta_ref
    path target_bed

    output:
    tuple val(meta), path("${meta.sample_id}.filtered.vcf.gz"), emit: filtered_vcf
    tuple val(meta), path("${meta.sample_id}.vcf.stats.txt"), emit: vcf_stats

    script:
    """
    bcftools mpileup \\
        --threads ${task.cpus} \\
        -Ou -f ${fasta_ref} \\
        -R ${target_bed} \\
        -a FORMAT/DP,FORMAT/AD \\
        ${dedup_bam} \\
    | bcftools call \\
        --threads ${task.cpus} \\
        -mv -Ou \\
    | bcftools filter \\
        -i 'QUAL>=30 && INFO/DP>=15' \\
        -Oz -o ${meta.sample_id}.filtered.vcf.gz

    bcftools index -t ${meta.sample_id}.filtered.vcf.gz
    bcftools stats ${meta.sample_id}.filtered.vcf.gz > ${meta.sample_id}.vcf.stats.txt
    """
}
