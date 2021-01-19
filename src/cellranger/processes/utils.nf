nextflow.enable.dsl=2

if(!params.containsKey("test")) {
    binDir = "${workflow.projectDir}/src/cellranger/bin/"
} else {
    binDir = ""
}

process SC__CELLRANGER__PREPARE_FOLDER {

    publishDir "${params.global.outdir}/data/raw/cellranger_fastq_folders", mode: 'symlink', overwrite: true
    label 'compute_resources__minimal'

    input:
        tuple val(sampleId), val(fastqs)

    output:
        tuple val(sampleId), path("${sampleId}_s*")

    script:
        def cmd = ''
        for(int i = 0; i < fastqs.size(); i++) {
            cmd += "mkdir ${sampleId}_s${i+1}; "
            cmd += "ln -s ${fastqs[i].join(' ')} ${sampleId}_s${i+1}; "
        }
        """
        $cmd
        """

}
