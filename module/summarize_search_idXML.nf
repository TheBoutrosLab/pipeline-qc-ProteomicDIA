process summarize_search_idXML {
    container params.docker_image_pyopenms
    containerOptions params.project_script_mount_options

    publishDir "${params.output_dir_base}/output/search_idXML", mode: 'copy'

    label 'process_low'

    input:
        tuple(
            val(sample_name),
            path(input_file)
        )

    output:
        path("${sample_name}_search_idXML_metrics.tsv"), emit: metrics

    script:
    """
    python3 ${projectDir}/src/SummarizeSearchIdXML.py \\
        --sample-name ${sample_name} \\
        --input ${input_file} \\
        --dataset-id ${params.dataset_id} \\
        --output ${sample_name}_search_idXML_metrics.tsv
    """
}
