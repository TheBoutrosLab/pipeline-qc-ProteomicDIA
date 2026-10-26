process extract_MS1_matrix {
    container params.docker_image_pyopenms
    containerOptions params.project_script_mount_options
    publishDir "${params.output_dir_base}/output/mzML", mode: 'copy'
    label 'process_low'
    input:
        tuple val(sample_name), path(input_file)
    output:
        tuple val(sample_name), path("${sample_name}_MS1_matrix.tsv"),
              path("${sample_name}_MS1_matrix.json"), emit: matrix
    script:
    """
    python3 ${projectDir}/src/ExtractMS1Matrix.py --sample-name ${sample_name} --input ${input_file} --matrix ${sample_name}_MS1_matrix.tsv --metadata ${sample_name}_MS1_matrix.json
    """
}
