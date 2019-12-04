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

def normalizeSRAFastqs(fastqPaths, index, sampleName) {
    /*
     * Rename samples SRRXXXXXX_[1|2].fastq.gz to more comprehensive file name ${sampleName}_S1_L00${index+1}_R[1|2]_001.fastq.gz
     * Here we follow 10xGenomics file naming convention
     */
    return fastqPaths.collect {
        fastqPath ->
            (full, parentDir, srrId, readType) = (fastqPath =~ /(.+)\/(SRR[0-9]*)_([1-2]).fastq.gz/)[0]
            normalizedFastQName = "${sampleName}_S1_L00${index+1}_R${readType}_001.fastq.gz"
            return [fastqPath, normalizedFastQName]
    }
}

process NORMALIZE_SRA_FASTQS {

    // publishDir "${params.global.outdir}/data", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=1 -l walltime=1:00:00 -A ${params.qsubaccount}"

    input:
    tuple val(sampleId), val(fastqs)
    
    output:
    tuple val(sampleId), path("*.fastq.gz")
    
    script:
    def normalizedFastqs = fastqs
        .withIndex()
        .collect {
            fastq, index -> normalizeSRAFastqs(fastq, index, sampleId)
        }
        .flatten()
        .collate(2)
    def cmd = ''
    for(int i = 0; i < normalizedFastqs.size(); i++)
        cmd += "ln -s ${normalizedFastqs[i][0]} ${normalizedFastqs[i][1]}; " 
    """
    $cmd
    """

}
