import static groovy.json.JsonOutput.*

nextflow.preview.dsl=2

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

/* 
    ATAC-seq pipelines
*/


// runs mkfastq, then cellranger-atac count:
workflow cellranger_atac {

    include {
        CELLRANGER_ATAC
    } from './src/cellranger-atac/main.nf' params(params)
    
    CELLRANGER_ATAC(
        file(params.sc.cellranger_atac.mkfastq.csv),
        file(params.sc.cellranger_atac.mkfastq.runFolder),
        file(params.sc.cellranger_atac.count.reference)
    )

}


workflow cistopic {

    include {
        cistopic as CISTOPIC
    } from './src/cistopic/main' params(params)
    
    getDataChannel | CISTOPIC

}


workflow atac_preprocess {

    // generic ATAC-seq preprocessing pipeline: adapter trimming, mapping, fragments file generation
    include {
        ATAC_PREPROCESS_WITH_METADATA;
    } from './workflows/atac/preprocess.nf' params(params)

    ATAC_PREPROCESS_WITH_METADATA(file(params.data.atac_preprocess.metadata))

}

workflow atac_qc_filtering {

    include {
        ATAC_QC_PREFILTER;
    } from './workflows/atac/qc_filtering.nf' params(params)

    getDataChannel | ATAC_QC_PREFILTER

}

workflow atac_preprocess_freemuxlet {

    // generic ATAC-seq preprocessing pipeline: adapter trimming, mapping, fragments file generation
    include {
        ATAC_PREPROCESS_WITH_METADATA;
    } from './workflows/atac/preprocess.nf' params(params)
    include {
        freemuxlet as FREEMUXLET;
    } from './workflows/popscle' params(params)

    ATAC_PREPROCESS_WITH_METADATA(file(params.sc.atac.preprocess.metadata))
    FREEMUXLET(ATAC_PREPROCESS_WITH_METADATA.out.bam)
}

