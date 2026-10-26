process plot_MS1_landscape {
    container params.docker_image_plotly
    containerOptions params.project_script_mount_options
    publishDir "${params.output_dir_base}/output/mzML", mode: 'copy'
    label 'process_low'
    input:
        tuple val(sample_name), path(matrix_file), path(metadata_file)
    output:
        path("${sample_name}_MS1_landscape.html"), emit: html
    script:
    """
    python3 ${projectDir}/src/PlotMS1Landscape.py --matrix ${matrix_file} --metadata ${metadata_file} --output ${sample_name}_MS1_landscape.html
    """
}
