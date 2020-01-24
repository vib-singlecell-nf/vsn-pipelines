nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/utils/bin/" : ""

process SC__UTILS__EXTRACT_FEATURE_METADATA {

    container params.sc.scanpy.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=2 -l walltime=1:00:00 -A ${params.global.qsubaccount}"

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), path("${sampleId}.SC__UTILS__EXTRACT_FEATURE_METADATA.tsv")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.utils.extract_feature_metadata)
		processParams = sampleParams.local
        columnNamesAsArguments = processParams.columnNames.collect({ '--column-name' + ' ' + it }).join(' ')
        """
        ${binDir}sc_h5ad_extract_metadata.py \
            --axis feature \
            ${columnNamesAsArguments} \
            $f \
            "${sampleId}.SC__UTILS__EXTRACT_FEATURE_METADATA.tsv"
        """

}
