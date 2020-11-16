nextflow.preview.dsl=2

//////////////////////////////////////////////////////
// process imports:
/*
include { SC__ARCHR__CREATE_ARROW_UNFILTERED; } from './../../src/archr/processes/createArrow_unfiltered.nf' params(params)
include { SC__ARCHR__CELL_CALLING; } from './../../src/archr/processes/cell_calling.nf' params(params)

include { SC__PYCISTOPIC__MACS2_CALL_PEAKS; } from './../../src/pycistopic/processes/macs2_call_peaks.nf' params(params)
*/

//////////////////////////////////////////////////////
//  Define the workflow 

workflow ATAC_QC_PREFILTER {

    take:
        data

    main:

        data.view()
        //SC__PYCISTOPIC__MACS2_CALL_PEAKS(data)

        //SC__PYCISTOPIC__PREFILTER(data)

}


/*
workflow ATAC_QC_FILTERING {

    take:
        data

    main:
        
        data.view()
        //data.branch {
        //    fragments: it[3] == 'fragments'
        //    bam: it[3] == 'bam'
        //}
        //.set{ data_split }
        //data_split.fragments.view()

        SC__ARCHR__CREATE_ARROW_UNFILTERED(data) |
            SC__ARCHR__CELL_CALLING

    /*
        cell_calling
        doublet_freemuxlet
        bap_multiplets
    */
//}

