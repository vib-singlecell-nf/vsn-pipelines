nextflow.enable.dsl=2

include {
    EDIRECT__SRAID_TO_SAMPLENAME
} from '../processes/sra_fastq_urls.nf'

workflow SRA_FASTQ_URLS {

    take:
        sraProjectId
        sampleNamesToRetrieve

    main:
        Channel
            .fromSRA(sraProjectId)
            .map { it[0] }
            .set { sraIDs }
        sraIDsToSample = EDIRECT__SRAID_TO_SAMPLENAME( sraIDs )
        sraFastqUrls = sraIDsToSample
            .join(sraIDs)
            .map { it -> tuple(it[0],it[1],"ftp://ftp.sra.ebi.ac.uk/" + it[2])} 

        if(!params.containsKey('quiet')) sraFastqUrls.view()

    emit:
        sraFastqUrls

}
