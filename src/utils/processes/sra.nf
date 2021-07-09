nextflow.enable.dsl=2

import java.nio.file.Files
import java.nio.file.Paths

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/utils/bin" : Paths.get(workflow.scriptFile.getParent().getParent().toString(), "utils/bin")


processParams = params.utils.sra_metadata

process GET_SRA_DB {

    container params.utils.container
    publishDir "${processParams.sraDbOutDir}", mode: 'link', overwrite: true
    label 'compute_resources__default'

    output:
        file("SRAmetadb.sqlite")

    script:
        """
        pysradb metadb \
            --out-dir "."
        """

}

process SRA_TO_METADATA {

    container params.utils.container
    publishDir "${params.global.outdir}/metadata", mode: 'link', overwrite: true
    errorStrategy 'retry'
    maxRetries 5
    label 'compute_resources__default'

    input:
        tuple val(sraId), val(sampleFilters)
        file(sraDb)

    output:
        file "${sraId}_metadata.tsv"

    script:
        if(processParams.mode == 'db') {
            if(sraDb.name != 'NO_FILE') {
                sraDbAsArgument = "--sra-db ${sraDb}"
            } else {
                if(!processParams.containsKey('sraDb') || processParams.sraDb == '')
                    throw new Exception("The db modue requires sraDb to be specified")
                sraDbAsArgument = '--sra-db ' + processParams.sraDb
            }
        } else if(processParams.mode == 'web') {
            sraDbAsArgument = ''
        } else {
            throw new Exception("The "+ processParams.mode +" mode does not exist. Choose one of: web, db.")
        }
        def sampleFiltersAsArguments = sampleFilters.collect({ '--sample-filter' + ' "' + it + '"'}).join(' ')
        """
        ${binDir}/sra_to_metadata.py \
            ${sraId} \
            ${sraDbAsArgument} \
            ${sampleFiltersAsArguments} \
            --output "${sraId}_metadata.tsv"
        """

}

def normalizeSRAFastQ(fastQPath, sampleName, fastqReadSuffixes) {
    /*
     * Rename samples SRRXXXXXX_[1-9].fastq.gz to more comprehensive file name ${sampleName}_S1_L001_${fastqReadSuffixes[\1-1]}_001.fastq.gz
     * Here we follow 10xGenomics file naming convention
     */
    (full, srrId, readType) = (fastQPath =~ /(SRR[0-9]*)_([1-9]).fastq.gz/)[0]

    if(readType.toInteger()-1 >= fastqReadSuffixes.size()) {
        throw new Exception("Read suffix for the current FASTQ file from "+ srrId +"SRA ID with index "+ readType + " cannot be extracted from params.utils.sra_normalize_fastqs.fastq_read_suffixes.")
    }
    readSuffix = fastqReadSuffixes[readType.toInteger()-1]
    normalizedFastQName = "${sampleName}_S1_L001_${readSuffix}_001.fastq.gz"
    return [fastQPath, normalizedFastQName]
}

process NORMALIZE_SRA_FASTQS {

    // publishDir "${params.global.outdir}/data", mode: 'symlink'
    label 'compute_resources__default'

    input:
        tuple val(sampleId), file(fastqs)

    output:
        tuple val(sampleId), path("*.fastq.gz")

    script:
        def normalizedFastqs = fastqs
            .collect {
                fastq -> normalizeSRAFastQ(fastq, sampleId, params.utils.sra_normalize_fastqs.fastq_read_suffixes)
            }
        def cmd = ''
        for(int i = 0; i < normalizedFastqs.size(); i++)
            cmd += "ln -s ${normalizedFastqs[i][0]} ${normalizedFastqs[i][1]}; " 
        """
        $cmd
        """

}
