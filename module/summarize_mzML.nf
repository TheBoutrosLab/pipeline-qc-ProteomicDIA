process summarize_mzML {
    container params.docker_image_pyopenms
    containerOptions params.project_script_mount_options

    publishDir { "${params.output_dir_base}/output/mzML/${sample_name}/${fraction_id}" }, mode: 'copy'

    label 'process_low'

    input:
        tuple val(sample_name), val(fraction_id), path(input_file)

    output:
        path("${sample_name}_${fraction_id}_mzML_metrics.tsv"), emit: metrics

    script:
    """
    python3 ${projectDir}/src/SummarizeMzML.py \\
        --sample-name ${sample_name} \\
        --fraction-id ${fraction_id} \\
        --input ${input_file} \\
        --dataset-id ${params.dataset_id} \\
        --output ${sample_name}_${fraction_id}_mzML_metrics.tsv
    """
}
