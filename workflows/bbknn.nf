/*
 * Remote run test
 * Source:
 * 
 * Steps considered: 

 */ 
import static groovy.json.JsonOutput.*

nextflow.preview.dsl=2

// print all parameters:
// println(prettyPrint(toJson( params )))


//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include QC_FILTER from '../src/scanpy/workflows/qc_filter.nf' params(params)
include SC__FILE_CONCATENATOR from '../src/utils/processes/utils.nf' params(params.sc.file_concatenator + params.global + params)
include NORMALIZE_TRANSFORM from '../src/scanpy/workflows/normalize_transform.nf' params(params + params.global)
include HVG_SELECTION from '../src/scanpy/workflows/hvg_selection.nf' params(params + params.global)
include DIM_REDUCTION from '../src/scanpy/workflows/dim_reduction.nf' params(params + params.global)
include CLUSTER_IDENTIFICATION from '../src/scanpy/workflows/cluster_identification.nf' params(params + params.global)
include SC__H5AD_TO_FILTERED_LOOM from '../src/utils/processes/h5adToLoom.nf' params(params + params.global)
include SC__H5AD_TO_LOOM from '../../utils/processes/h5adToLoom.nf' params(params)
include BEC_BBKNN from '../src/scanpy/workflows/bec_bbknn.nf' params(params)

// data channel to start from 10x data:
include getChannel as getTenXChannel from '../src/channels/tenx.nf' params(params)

// reporting:
include SC__SCANPY__MERGE_REPORTS from '../src/scanpy/processes/reports.nf' params(params + params.global)
include SC__SCANPY__REPORT_TO_HTML from '../src/scanpy/processes/reports.nf' params(params + params.global)


workflow bbknn_base {

    get:
        data

    main:
        // run the pipeline
        QC_FILTER( data ) // Remove concat 
        SC__FILE_CONCATENATOR( QC_FILTER.out.filtered.map{it -> it[1]}.collect() )
        NORMALIZE_TRANSFORM( SC__FILE_CONCATENATOR.out )
        HVG_SELECTION( NORMALIZE_TRANSFORM.out )

        //// include all pre-merge dim reductions. These will be replaced in the bbknn step
        DIM_REDUCTION( HVG_SELECTION.out.scaled )
        CLUSTER_IDENTIFICATION( DIM_REDUCTION.out.dimred )
        BEC_BBKNN(
            SC__FILE_CONCATENATOR.out.join(CLUSTER_IDENTIFICATION.out.marker_genes)
        )

        // conversion
        //// convert h5ad to X (here we choose: loom format)
        filteredloom = SC__H5AD_TO_FILTERED_LOOM( SC__FILE_CONCATENATOR.out )
        scopeloom = SC__H5AD_TO_LOOM(
            QC_FILTER.out.filteredloom.join(BEC_BBKNN.out.data)
        )

        // collect the reports:
        ipynbs = HVG_SELECTION.out.report
            .join(BEC_BBKNN.out.cluster_report)
            .join(BEC_BBKNN.out.bbknn_report)
            .map{ tuple( it[0], it.drop(1) ) }
        // reporting:
        SC__SCANPY__MERGE_REPORTS(ipynbs, "merged_report")
        SC__SCANPY__REPORT_TO_HTML(SC__SCANPY__MERGE_REPORTS.out)

    emit:
        filteredloom
        scopeloom

}

workflow bbknn_standalone {

    main:
        data = getTenXChannel( params.global.tenx_folder ).view()
        bbknn_base( data )

    emit:
        filteredloom = bbknn_base.out.filteredloom
        scopeloom = bbknn_base.out.scopeloom

}

workflow bbknn {

    get:
        data

    main:
        bbknn_base( data )

    emit:
        bbknn_base.out.filteredloom
        bbknn_base.out.scopeloom

}
