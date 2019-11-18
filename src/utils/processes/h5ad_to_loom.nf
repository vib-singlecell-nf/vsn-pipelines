nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  binDir = "${workflow.projectDir}/src/utils/bin/"
} else {
  binDir = ""
}


process SC__H5AD_TO_LOOM {

  container params.sc.scanpy.container
  clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  publishDir "${params.outdir}/loom", mode: 'link', overwrite: true

  input:
    tuple val(id), file(f)
  output:
    tuple val(id), file("${id}.SC__H5AD_TO_LOOM.loom")
  script:
    """
    ${binDir}h5ad_to_loom.py \
         $f \
         "${id}.SC__H5AD_TO_LOOM.loom" 
    """
}

process SC__H5AD_TO_FILTERED_LOOM {

  container params.sc.scanpy.container
  clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  publishDir "${params.outdir}/loom", mode: 'link', overwrite: true

  input:
    tuple val(id), file(f)
  output:
    tuple val(id), file("${id}.filtered.loom")
  script:
    """
    ${binDir}h5ad_to_filtered_loom.py \
         $f \
         "${id}.filtered.loom"
    """
}

