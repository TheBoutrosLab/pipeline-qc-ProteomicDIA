process collect_QCMetrics {
    container params.docker_image_pyopenms
    containerOptions params.project_script_mount_options
    publishDir "${params.output_dir_base}/output", mode: 'copy'
    label 'process_low'
    input:
        path(input_metrics)
    output:
        path("${params.dataset_id}_QC_metrics.tsv"), emit: metrics
    script:
    """
    python3 ${projectDir}/src/CollectQCMetrics.py --inputs ${input_metrics.join(' ')} --output ${params.dataset_id}_QC_metrics.tsv
    """
}
