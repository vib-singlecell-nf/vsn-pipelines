nextflow.enable.dsl=2

//////////////////////////////////////////////////////
// process imports:
include { SC__ARCHR__CREATE_ARROW_UNFILTERED; } from './../../src/archr/processes/createArrow_unfiltered.nf'
include { SC__ARCHR__CELL_CALLING; } from './../../src/archr/processes/cell_calling.nf'

include { PYCISTOPIC__BIOMART_ANNOT; } from './../../src/pycistopic/processes/biomart_annot.nf'
include { PYCISTOPIC__MACS2_CALL_PEAKS; } from './../../src/pycistopic/processes/macs2_call_peaks.nf'
include {
    rename_fragments;
    PYCISTOPIC__COMPUTE_QC_STATS;
} from './../../src/pycistopic/processes/compute_qc_stats.nf'
include {
    SCTK__SATURATION;
    SCTK__SATURATION as SCTK__SATURATION_BC_WL;
} from './../../src/singlecelltoolkit/processes/saturation.nf'
include {
    PYCISTOPIC__QC_REPORT;
    REPORT_TO_HTML;
} from './../../src/pycistopic/processes/call_cells.nf'

include {
    SIMPLE_PUBLISH as PUBLISH_PEAKS;
    SIMPLE_PUBLISH as PUBLISH_SUMMITS;
    SIMPLE_PUBLISH as PUBLISH_SATURATION_TSV;
    SIMPLE_PUBLISH as PUBLISH_SATURATION_PNG;
    SIMPLE_PUBLISH as PUBLISH_SATURATION_BC_WL_TSV;
    SIMPLE_PUBLISH as PUBLISH_SATURATION_BC_WL_PNG;
} from '../../src/utils/processes/utils.nf'


//////////////////////////////////////////////////////
//  Define the workflow 


workflow cellranger_output_to_bam_fragments {
    /*
    Cell Ranger ATAC::
       possorted_bam.bam, fragments.tsv.gz
    Cell Ranger ARC::
       atac_possorted_bam.bam, atac_fragments.tsv.gz
    */

    take:
        data // standard data channel [ sampleId, path, type, format]

    main:

        bam = data.map{ it -> tuple(it[0], [
                                    file(it[1]+"/*possorted*bam.bam")[0],
                                    file(it[1]+"/*possorted*bam.bam.bai")[0],
                                    ]) }
        fragments = data.map{ it -> tuple(it[0], [
//                                    *file(it[1]+"/*fragments.tsv.gz"),
//                                    *file(it[1]+"/*fragments.tsv.gz.tbi"),
                                    file(it[1]+"/*fragments.tsv.gz")[0],
                                    file(it[1]+"/*fragments.tsv.gz.tbi")[0],
                                    ]) }

        if(!params.containsKey('quiet')) bam.view()
        if(!params.containsKey('quiet')) fragments.view()

    emit:
        bam
        fragments

}


workflow ATAC_QC_PREFILTER {

    take:
        data

    main:

        data.branch {
            fragments:  it[2] == 'fragments'
            bam:        it[2] == 'bam'
            cellranger: it[2] == '10x_atac_cellranger_mex_outs'
        }
        .set{ data_split }

        // split the cellranger data into separate bam and fragments channels
        data_split.cellranger \
            | cellranger_output_to_bam_fragments
            | set { data_cr }
        /* 'mix' the separate bam and fragments channels with the
           cellranger bam and fragments files, and use these channels going
           forward
         */
        bam = data_split.bam.mix(data_cr.bam)
        /* for fragments, rename the files to include the sample ID
           prefix (necessary for cellranger inputs, which all have the same
           file name). This is not currently necessary for the bam files since
           they are processed in separate processes.
         */
        /*fragments = rename_fragments(
            data_split.fragments.mix(data_cr.fragments)
            )
        */
        fragments = data_split.fragments.mix(data_cr.fragments)


        biomart = PYCISTOPIC__BIOMART_ANNOT()

        peaks = PYCISTOPIC__MACS2_CALL_PEAKS(bam.map { it -> tuple(it[0], it[1][0], it[1][1] ) } )
        PUBLISH_PEAKS(peaks.map { it -> tuple(it[0], it[1]) }, '.peaks.narrowPeak', 'macs2')
        PUBLISH_SUMMITS(peaks.map { it -> tuple(it[0], it[2]) }, '.summits.bed', 'macs2')

        /* pycisTopic qc: pass every fragment/peak file into a single process
           together. These will be formatted as a string "sampleId,fragments,peak",
           which is parsed in the python script. The fragments and peaks files
           are staged separately
        */
        fragments.map { it -> tuple(it[0], it[1][0].getName(), it[1][1].getName() ) } // [sampleId, fragments, fragments.tbi]
                 .join(peaks.map{ it -> tuple(it[0], it[1].getName()) }) // combine with peaks for each sample
                 .map { it -> ["${it[0]},${it[1]},${it[3]}"] } // join as string
                 .collect() // collapse to a single channel element
                 .set { fragpeaks }

        qc_stats = PYCISTOPIC__COMPUTE_QC_STATS(fragpeaks,
                                                biomart,
                                                fragments.map { it -> it[1][0] }.collect(),
                                                peaks.map { it -> it[1]}.collect()
                                                )

        PYCISTOPIC__QC_REPORT(
            file(workflow.projectDir + params.tools.pycistopic.call_cells.report_ipynb),
            fragments.map { it -> it[0] }.collect(), // all sampleIds
            qc_stats,
            params.global.project_name + "__pycisTopic_QC_report"
        ) \
        | map { it -> it[0] }
        | REPORT_TO_HTML

        /* saturation */
        if(! params.tools.singlecelltoolkit.saturation.skip) {
            SCTK__SATURATION(fragments.map { it -> tuple(it[0], it[1][0], it[1][1] ) }, '', '')
            SCTK__SATURATION_BC_WL(fragments.map { it -> tuple(it[0], it[1][0], it[1][1] ) },
                PYCISTOPIC__QC_REPORT.out, 'RUN')

            /* publish saturation outputs */
            PUBLISH_SATURATION_TSV(SCTK__SATURATION.out.map { it -> tuple(it[0], it[1]) }, '.sampling_stats.tsv', 'singlecelltoolkit/saturation')
            PUBLISH_SATURATION_PNG(SCTK__SATURATION.out.map { it -> tuple(it[0], it[2]) }, '.saturation.png', 'singlecelltoolkit/saturation')

            PUBLISH_SATURATION_BC_WL_TSV(SCTK__SATURATION_BC_WL.out.map { it -> tuple(it[0], it[1]) }, '.sampling_stats.tsv', 'singlecelltoolkit/saturation_bc_wl')
            PUBLISH_SATURATION_BC_WL_PNG(SCTK__SATURATION_BC_WL.out.map { it -> tuple(it[0], it[2]) }, '.saturation.png', 'singlecelltoolkit/saturation_bc_wl')
        }
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

        SC__ARCHR__CREATE_ARROW_UNFILTERED(data) |
            SC__ARCHR__CELL_CALLING

}

