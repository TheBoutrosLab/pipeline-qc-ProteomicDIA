process summarize_DIANN_report {
    container params.docker_image_pyopenms
    containerOptions params.project_script_mount_options
    publishDir "${params.output_dir_base}/output/DIA-NN", mode: 'copy'
    label 'process_low'
    input:
        path(input_file)
        path(input_samples)
    output:
        path("${params.dataset_id}_DIANN_metrics.tsv"), emit: metrics
    script:
    """
    python3 ${projectDir}/src/SummarizeDIANNReport.py --input ${input_file} --sample-csv ${input_samples} --dataset-id ${params.dataset_id} --output ${params.dataset_id}_DIANN_metrics.tsv
    """
}
