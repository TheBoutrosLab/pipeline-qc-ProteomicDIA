process extract_MS1_matrix {
    container params.docker_image_pyopenms
    containerOptions params.project_script_mount_options

    publishDir { "${params.output_dir_base}/output/mzML/${sample_name}/${fraction_id}" }, mode: 'copy'

    label 'process_low'

    input:
        tuple(
            val(sample_name),
            val(fraction_id),
            path(input_file)
        )

    output:
        tuple(
            val(sample_name),
            val(fraction_id),
            path("${sample_name}_${fraction_id}_MS1_matrix.tsv"),
            path("${sample_name}_${fraction_id}_MS1_matrix.json"),
            emit: matrix
        )

    script:
    """
    python3 ${projectDir}/src/ExtractMS1Matrix.py \\
        --sample-name ${sample_name} \\
        --fraction-id ${fraction_id} \\
        --input ${input_file} \\
        --matrix ${sample_name}_${fraction_id}_MS1_matrix.tsv \\
        --metadata ${sample_name}_${fraction_id}_MS1_matrix.json
    """
}
