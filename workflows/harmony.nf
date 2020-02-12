nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include '../src/utils/processes/utils.nf' params(params.sc.file_concatenator + params.global + params)

include QC_FILTER from '../src/scanpy/workflows/qc_filter.nf' params(params)
include NORMALIZE_TRANSFORM from '../src/scanpy/workflows/normalize_transform.nf' params(params + params.global)
include HVG_SELECTION from '../src/scanpy/workflows/hvg_selection.nf' params(params + params.global)
include DIM_REDUCTION from '../src/scanpy/workflows/dim_reduction.nf' params(params + params.global)
// CLUSTER_IDENTIFICATION
include '../src/scanpy/processes/cluster.nf' params(params + params.global)
include '../src/scanpy/workflows/cluster_identification.nf' params(params + params.global) // Don't only import a specific process (the function needs also to be imported)
include BEC_HARMONY from '../src/harmony/workflows/bec_harmony.nf' params(params)

include SC__H5AD_TO_FILTERED_LOOM from '../src/utils/processes/h5adToLoom.nf' params(params + params.global)
include FILE_CONVERTER from '../src/utils/workflows/fileConverter.nf' params(params)

// data channel to start from 10x data:
include getChannel as getTenXChannel from '../src/channels/tenx.nf' params(params)

// reporting:
include UTILS__GENERATE_WORKFLOW_CONFIG_REPORT from '../src/utils/processes/reports.nf' params(params)
include SC__SCANPY__MERGE_REPORTS from '../src/scanpy/processes/reports.nf' params(params + params.global)
include SC__SCANPY__REPORT_TO_HTML from '../src/scanpy/processes/reports.nf' params(params + params.global)


workflow harmony_base {

    take:
        data

    main:
        // Run the pipeline
        QC_FILTER( data ) // Remove concat 
        SC__FILE_CONCATENATOR( QC_FILTER.out.filtered.map{it -> it[1]}.collect() )
        NORMALIZE_TRANSFORM( SC__FILE_CONCATENATOR.out )
        HVG_SELECTION( NORMALIZE_TRANSFORM.out )
        DIM_REDUCTION( HVG_SELECTION.out.scaled )

        // Perform the clustering step w/o batch effect correction (for comparison matter)
        clusterIdentificationPreBatchEffectCorrection = CLUSTER_IDENTIFICATION( 
            NORMALIZE_TRANSFORM.out,
            DIM_REDUCTION.out.dimred_pca_tsne_umap,
            "Pre Batch Effect Correction"
        )

        // Perform the batch effect correction
        BEC_HARMONY(
            NORMALIZE_TRANSFORM.out,
            // include only PCA since Harmony will correct this
            DIM_REDUCTION.out.dimred_pca.map { it -> tuple(it[0], it[1]) },
            clusterIdentificationPreBatchEffectCorrection.marker_genes
        )
        
        // Conversion
        // Convert h5ad to X (here we choose: loom format)
        filteredloom = SC__H5AD_TO_FILTERED_LOOM( SC__FILE_CONCATENATOR.out )
        scopeloom = FILE_CONVERTER(
            BEC_HARMONY.out.data.groupTuple(),
            'loom',
            SC__FILE_CONCATENATOR.out
        )

        project = CLUSTER_IDENTIFICATION.out.marker_genes.map { it -> it[0] }
        UTILS__GENERATE_WORKFLOW_CONFIG_REPORT(
            file(workflow.projectDir + params.utils.workflow_configuration.report_ipynb)
        )
        // collect the reports:
        ipynbs = project.combine(
            UTILS__GENERATE_WORKFLOW_CONFIG_REPORT.out
        ).join(
            HVG_SELECTION.out.report
        ).join(
            BEC_HARMONY.out.cluster_report
        ).combine(
            BEC_HARMONY.out.harmony_report,
            by: 0
        ).map { 
            tuple( it[0], it.drop(1) ) 
        }
        // reporting:
        def clusteringParams = SC__SCANPY__CLUSTERING_PARAMS( clean(params.sc.scanpy.clustering) )
        SC__SCANPY__MERGE_REPORTS(
            ipynbs,
            "merged_report",
            clusteringParams.isBenchmarkMode()
        )
        SC__SCANPY__REPORT_TO_HTML(SC__SCANPY__MERGE_REPORTS.out)

    emit:
        filteredloom
        scopeloom

}

workflow harmony_standalone {

    main:
        data = getTenXChannel( params.data.tenx.cellranger_outs_dir_path ).view()
        harmony_base( data )

    emit:
        filteredloom = harmony_base.out.filteredloom
        scopeloom = harmony_base.out.scopeloom

}

workflow harmony {

    take:
        data

    main:
        harmony_base( data )

    emit:
        filteredloom = harmony_base.out.filteredloom
        scopeloom = harmony_base.out.scopeloom

}
