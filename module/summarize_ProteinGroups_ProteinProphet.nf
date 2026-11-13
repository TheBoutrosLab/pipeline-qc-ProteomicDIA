process summarize_ProteinGroups_ProteinProphet {
    container params.docker_image_pyopenms
    containerOptions params.project_script_mount_options

    publishDir "${params.output_dir_base}/output/ProteinGroups", mode: 'copy'

    label 'process_low'

    input:
        path(input_file)

    output:
        path("${params.dataset_id}_ProteinGroups_metrics.tsv"), emit: metrics

    script:
    """
    python3 ${projectDir}/src/SummarizeProteinProphetGroups.py \\
        --input ${input_file} \\
        --dataset-id ${params.dataset_id} \\
        --output ${params.dataset_id}_ProteinGroups_metrics.tsv
    """
}
