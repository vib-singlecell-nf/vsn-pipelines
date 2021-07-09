nextflow.enable.dsl=2

import static groovy.json.JsonOutput.*



//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:
include {
    getBaseName;
} from '../utils/processes/files.nf'
include {
    FILTER_AND_ANNOTATE_AND_CLEAN;
} from '../utils/workflows/filterAnnotateClean.nf' params(params)
include { 
    INIT;
} from '../utils/workflows/utils' params(params)
INIT(params)
include {
    SC__FILE_CONVERTER;
    SC__FILE_CONCATENATOR;
} from '../utils/processes/utils' params(params)
include {
    getDataChannel
} from '../channels/channels' params(params)
include {
    QC_FILTER;
} from './workflows/qc_filter.nf' params(params)
include {
    SINGLE_SAMPLE
} from './workflows/single_sample.nf' params(params)

workflow single_sample {

    main:
        // run the pipeline
        getDataChannel | \
            SC__FILE_CONVERTER | \
            SINGLE_SAMPLE

}

workflow single_sample_qc {

    take:
        data

    main:
        out = SC__FILE_CONVERTER( data )
        out = FILTER_AND_ANNOTATE_AND_CLEAN( out )
        QC_FILTER( out )

}

workflow multi_sample_qc {

    take:
        data

    main:
        if(!params?.sc?.scanpy?.filter) {
            throw new Exception("VSN ERROR: Missing params.sc.scanpy.filter config.")
        }
        if(!params?.sc?.file_concatenator) {
            throw new Exception("VSN ERROR: Missing params.sc.file_concatenator config.")
        }

        out = data | \
            SC__FILE_CONVERTER | \
            FILTER_AND_ANNOTATE_AND_CLEAN

        out = QC_FILTER( out ).filtered
        out = SC__FILE_CONCATENATOR( 
            out.map {
                it -> it[1]
            }.toSortedList( 
                { a, b -> getBaseName(a, "SC") <=> getBaseName(b, "SC") }
            )
        )

    emit:
        filtered = out

}