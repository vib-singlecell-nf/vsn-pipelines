import static groovy.json.JsonOutput.*

nextflow.enable.dsl=2

include {
    INIT;
} from './src/utils/workflows/utils' params(params)

INIT(params)

include {
    SC__FILE_CONVERTER;
} from './src/utils/processes/utils' params(params)

include {
    getDataChannel;
} from './src/channels/channels' params(params)


include {
    CELLRANGER_ATAC
} from './src/cellranger-atac/main.nf' params(params)

include {
    ATAC_PREPROCESS;
} from './workflows/atac/preprocess.nf' params(params)

include {
    ATAC_PREPROCESS_RAPID;
} from './workflows/atac/preprocess_rapid.nf' params(params)

include {
    ATAC_QC_PREFILTER;
} from './workflows/atac/qc_filtering.nf' params(params)

include {
    ATAC_PREPROCESS_WITH_METADATA;
} from './workflows/atac/preprocess.nf' params(params)

include {
    get_bam;
    BAP__BARCODE_MULTIPLET_WF;
} from './src/bap/main.nf' params(params)

include {
    freemuxlet as FREEMUXLET;
} from './workflows/popscle' params(params)

/*
    ATAC-seq pipelines
*/


// runs mkfastq, then cellranger-atac count:
workflow cellranger_atac {

    CELLRANGER_ATAC(
        file(params.tools.cellranger_atac.mkfastq.csv),
        file(params.tools.cellranger_atac.mkfastq.runFolder),
        file(params.tools.cellranger_atac.count.reference)
    )

}


workflow atac_preprocess {

    // generic ATAC-seq preprocessing pipeline: adapter trimming, mapping, fragments file generation
    ATAC_PREPROCESS(file(params.data.atac_preprocess.metadata))

}


workflow atac_preprocess_bap {

    ATAC_PREPROCESS(file(params.data.atac_preprocess.metadata)) |
        get_bam |
        BAP__BARCODE_MULTIPLET_WF

}

workflow atac_preprocess_rapid {

    ATAC_PREPROCESS_RAPID(file(params.data.atac_preprocess.metadata))

}

/*
 QC
 */
workflow atac_qc_filtering {

    getDataChannel | ATAC_QC_PREFILTER

}

workflow atac_preprocess_with_qc {

    // generic ATAC-seq preprocessing pipeline: adapter trimming, mapping, fragments file generation
    pp = ATAC_PREPROCESS(file(params.data.atac_preprocess.metadata))
    ATAC_QC_PREFILTER(pp.bam.mix(pp.fragments))

}

workflow atac_preprocess_freemuxlet {

    // generic ATAC-seq preprocessing pipeline: adapter trimming, mapping, fragments file generation
    ATAC_PREPROCESS_WITH_METADATA(file(params.tools.atac.preprocess.metadata))
    FREEMUXLET(ATAC_PREPROCESS_WITH_METADATA.out.bam)
}

workflow bap {

    getDataChannel | BAP__BARCODE_MULTIPLET_WF

}