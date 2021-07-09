nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the utils module:
////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    DOWNLOAD_FASTQS_FROM_SRA_ACC_ID;
} from "../processes/downloadFastQ" params(params)
include {
    FIX_AND_COMPRESS_SRA_FASTQS;
} from "../../singlecelltoolkit/processes/fix_and_compress_fastqs" params(params)


workflow SRATOOLKIT__DOWNLOAD_FASTQS {

    take:
        // Expects (sraId, sampleId)
        data

    main:
        out = data | \
            DOWNLOAD_FASTQS_FROM_SRA_ACC_ID | \
            FIX_AND_COMPRESS_SRA_FASTQ

    emit:
        // Returns (sraId, *.fastq.gz)
        out

}