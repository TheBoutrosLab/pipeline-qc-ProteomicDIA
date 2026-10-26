process summarize_EasyPQP_library {
    container params.docker_image_pyopenms
    containerOptions params.project_script_mount_options
    publishDir "${params.output_dir_base}/output/EasyPQP", mode: 'copy'
    label 'process_low'
    input:
        path(input_library)
        path(input_summary)
    output:
        path("${params.dataset_id}_EasyPQP_metrics.tsv"), emit: metrics
    script:
    """
    python3 ${projectDir}/src/SummarizeEasyPQPLibrary.py --library ${input_library} --export-summary ${input_summary} --dataset-id ${params.dataset_id} --output ${params.dataset_id}_EasyPQP_metrics.tsv
    """
}
