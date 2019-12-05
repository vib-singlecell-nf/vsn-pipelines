nextflow.preview.dsl=2

if(!params.containsKey("test")) {
    binDir = "${workflow.projectDir}/src/cellranger/bin/"
} else {
    binDir = ""
}

process SC__CELLRANGER__PREPARE_FOLDER {

    clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
    // publishDir "${params.outdir}/data", mode: 'link', overwrite: true

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
