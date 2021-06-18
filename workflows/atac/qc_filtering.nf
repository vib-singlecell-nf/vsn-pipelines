nextflow.enable.dsl=2

//////////////////////////////////////////////////////
// process imports:
include { SC__ARCHR__CREATE_ARROW_UNFILTERED; } from './../../src/archr/processes/createArrow_unfiltered.nf'
include { SC__ARCHR__CELL_CALLING; } from './../../src/archr/processes/cell_calling.nf'

include { PYCISTOPIC__BIOMART_ANNOT; } from './../../src/pycistopic/processes/biomart_annot.nf'
include { PYCISTOPIC__MACS2_CALL_PEAKS; } from './../../src/pycistopic/processes/macs2_call_peaks.nf'
include { PYCISTOPIC__COMPUTE_QC_STATS; } from './../../src/pycistopic/processes/compute_qc_stats.nf'
include {
    PYCISTOPIC__QC_REPORT;
    REPORT_TO_HTML;
} from './../../src/pycistopic/processes/call_cells.nf'

include {
    SIMPLE_PUBLISH as PUBLISH_PEAKS;
    SIMPLE_PUBLISH as PUBLISH_SUMMITS;
} from '../../src/utils/processes/utils.nf'

//////////////////////////////////////////////////////
//  Define the workflow 

workflow ATAC_QC_PREFILTER {

    take:
        data

    main:

        data.branch {
            fragments: it[2] == 'fragments'
            bam:       it[2] == 'bam'
        }
        .set{ data_split }

        biomart = PYCISTOPIC__BIOMART_ANNOT()

        peaks = PYCISTOPIC__MACS2_CALL_PEAKS(data_split.bam.map { it -> tuple(it[0], it[1][0], it[1][1] ) } )
        PUBLISH_PEAKS(peaks.map { it -> tuple(it[0], it[1]) }, '.peaks.narrowPeak', 'macs2')
        PUBLISH_SUMMITS(peaks.map { it -> tuple(it[0], it[2]) }, '.summits.bed', 'macs2')

        data_split.fragments.map { it -> tuple(it[0], it[1][0], it[1][1] ) }
                            .join(peaks)
                            //.map { it -> [ tuple(*it[0..1], it[3]) ] }
                            .map { it -> ["${it[0]},${it[1]},${it[3]}"] }
                            .collect()
                            .set { fragpeaks }

        qc_stats = PYCISTOPIC__COMPUTE_QC_STATS(fragpeaks, biomart)

        PYCISTOPIC__QC_REPORT(
            file(workflow.projectDir + params.tools.pycistopic.call_cells.report_ipynb),
            data_split.fragments.map { it -> it[0] }.collect(), // all sampleIds
            qc_stats,
            params.global.project_name + "__pycisTopic_QC_report"
        ) \
        | map { it -> it[0] }
        | REPORT_TO_HTML

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

