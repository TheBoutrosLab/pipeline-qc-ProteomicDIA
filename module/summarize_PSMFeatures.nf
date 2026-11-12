process summarize_PSMFeatures {
    container params.docker_image_pyopenms
    containerOptions params.project_script_mount_options
    publishDir { "${params.output_dir_base}/output/PSMFeatures/${sample_name}/${fraction_id ?: 'all'}" }, mode: 'copy'
    label 'process_low'
    input:
        tuple val(sample_name), val(search_engine), val(route), val(fraction_id), path(input_file)
    output:
        path("${sample_name}_${search_engine}_${route}_${fraction_id ?: 'all'}_PSMFeatures_metrics.tsv"), emit: metrics
    script:
    """
    python3 ${projectDir}/src/SummarizePSMFeatures.py \\
        --sample-name ${sample_name} \\
        --search-engine ${search_engine} \\
        --route ${route} \\
        --fraction-id '${fraction_id}' \\
        --input ${input_file} \\
        --dataset-id ${params.dataset_id} \\
        --output ${sample_name}_${search_engine}_${route}_${fraction_id ?: 'all'}_PSMFeatures_metrics.tsv
    """
}
