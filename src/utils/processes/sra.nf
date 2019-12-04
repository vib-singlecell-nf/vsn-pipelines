nextflow.preview.dsl=2

import java.nio.file.Files
import java.nio.file.Paths

if(!params.containsKey("test")) {
    binDir = "${workflow.projectDir}/src/utils/bin/"
} else {
    binDir = ""
}

processParams = params.sc.utils.sra_metadata

process GET_SRA_DB {

    // container params.sc.utils.container
    publishDir "${processParams.sraDbOutDir}", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=1 -l walltime=1:00:00 -A ${params.qsubaccount}"

    output:
    file("SRAmetadb.sqlite")
    
    script:
    """
    pysradb metadb \
        --out-dir "."
    """

}

process SRA_TO_METADATA {

    // container params.sc.utils.container
    publishDir "${params.global.outdir}/metadata", mode: 'link', overwrite: true
    clusterOptions "-l nodes=1:ppn=1 -l walltime=1:00:00 -A ${params.qsubaccount}"

    input:
    tuple val(sraId), val(sampleFilters)
    file(sraDb)
    
    output:
    file "${sraId}_metadata.tsv"
    
    script:
    if(sraDb.name != 'NO_FILE') {
        sraDbAsArgument = "--sra-db ${sraDb}"
    } else {
        sraDbAsArgument = (processParams.containsKey('sraDb') && processParams.sraDb != '') ? '--sra-db ' + processParams.sraDb : ''
    }
    def sampleFiltersAsArguments = sampleFilters.collect({ '--sample-filter' + ' ' + it }).join(' ')
        """
        ${binDir}sra_to_metadata.py \
        ${sraId} \
        ${sraDbAsArgument} \
            ${sampleFiltersAsArguments} \
        --output "${sraId}_metadata.tsv"
    """

}
        """

}
