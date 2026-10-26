process summarize_ProteinGroups_ProtLoom {
    container params.docker_image_pyopenms
    containerOptions params.project_script_mount_options
    publishDir "${params.output_dir_base}/output/ProteinGroups", mode: 'copy'
    label 'process_low'
    input:
        path(input_files)
    output:
        path("${params.dataset_id}_ProteinGroups_metrics.tsv"), emit: metrics
    script:
    """
    python3 ${projectDir}/src/SummarizeProtLoomGroups.py --input-files ${input_files.join(' ')} --dataset-id ${params.dataset_id} --output ${params.dataset_id}_ProteinGroups_metrics.tsv
    """
}
