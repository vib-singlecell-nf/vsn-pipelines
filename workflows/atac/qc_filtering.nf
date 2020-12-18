nextflow.preview.dsl=2

//////////////////////////////////////////////////////
// process imports:
include { SC__ARCHR__CREATE_ARROW_UNFILTERED; } from './../../src/archr/processes/createArrow_unfiltered.nf' params(params)
include { SC__ARCHR__CELL_CALLING; } from './../../src/archr/processes/cell_calling.nf' params(params)

include { SC__PYCISTOPIC__MACS2_CALL_PEAKS; } from './../../src/pycistopic/processes/macs2_call_peaks.nf' params(params)
include { SC__PYCISTOPIC__COMPUTE_QC_STATS; } from './../../src/pycistopic/processes/compute_qc_stats.nf' params(params)

include {
    PUBLISH as PUBLISH_METADATA;
    PUBLISH as PUBLISH_PEAKS;
} from "../../src/utils/workflows/utils.nf" params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow ATAC_QC_PREFILTER {

    take:
        data

    main:

        data.branch {
            fragments: it[3] == 'fragments'
            bam: it[3] == 'bam'
        }
        .set{ data_split }

        peaks = SC__PYCISTOPIC__MACS2_CALL_PEAKS(data_split.bam)
        PUBLISH_PEAKS(peaks.map { it -> tuple(it[0], it[1]) }, 'peaks', 'narrowPeak', 'macs2', false)

        data_split.fragments.join(peaks)
                  .map { it -> tuple(it[0], it[1], it[2], it[4]) }
                  .set{ fragpeaks }

        qc_stats = SC__PYCISTOPIC__COMPUTE_QC_STATS(fragpeaks)
        PUBLISH_METADATA(qc_stats.map { it -> tuple(it[0], it[1]) }, 'metadata.tsv', 'gz', 'pycistopic', false)

}


workflow ATAC_QC_FILTERING {

    take:
        data

    main:
        
        data.branch {
            fragments: it[3] == 'fragments'
            bam: it[3] == 'bam'
        }
        .set{ data_split }

        SC__ARCHR__CREATE_ARROW_UNFILTERED(data_split.fragments) |
            SC__ARCHR__CELL_CALLING

}

