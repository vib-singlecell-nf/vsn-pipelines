nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the utils module:
include {
    clean;
} from '../src/utils/processes/utils.nf' params(params)
include {
    PUBLISH;
} from '../src/utils/workflows/utils.nf'
include {
    SC__H5AD_TO_FILTERED_LOOM;
} from '../src/utils/processes/h5adToLoom.nf' params(params)
include {
    UTILS__GENERATE_WORKFLOW_CONFIG_REPORT;
} from '../src/utils/processes/reports.nf' params(params)
include {
    FILE_CONVERTER;
} from '../src/utils/workflows/fileConverter.nf' params(params)


include {
    star as STAR;
} from '../workflows/star.nf' params(params)

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from external tools:
include {
    QC_FILTER;
} from '../src/scanpy/workflows/qc_filter.nf' params(params)
include {
    NORMALIZE_TRANSFORM;
} from '../src/scanpy/workflows/normalize_transform.nf' params(params)
include {
    HVG_SELECTION;
} from '../src/scanpy/workflows/hvg_selection.nf' params(params)
include {
    SC__SCANPY__REGRESS_OUT;
} from '../src/scanpy/processes/regress_out.nf' params(params)
include {
    NEIGHBORHOOD_GRAPH;
} from '../src/scanpy/workflows/neighborhood_graph.nf' params(params)
include {
    DIM_REDUCTION_PCA;
} from '../src/scanpy/workflows/dim_reduction_pca.nf' params(params)
include {
    DIM_REDUCTION_TSNE_UMAP;
} from '../src/scanpy/workflows/dim_reduction.nf' params(params)
include {
    SC__SCANPY__CLUSTERING_PARAMS;
} from '../src/scanpy/processes/cluster.nf' params(params)
include {
    CLUSTER_IDENTIFICATION;
} from '../src/scanpy/workflows/cluster_identification.nf' params(params)
include {
    SC__DIRECTS__SELECT_DEFAULT_CLUSTERING
} from '../src/directs/processes/selectDefaultClustering.nf'

// reporting:
include {
    SC__SCANPY__MERGE_REPORTS;
} from '../src/scanpy/processes/reports.nf' params(params)
include {
    SC__SCANPY__REPORT_TO_HTML;
} from '../src/scanpy/processes/reports.nf' params(params)

workflow single_sample_star {
    
    /*******************************************
     * Data processing
     */
    data = STAR()
    samples = data.map { it -> it[0] }
    UTILS__GENERATE_WORKFLOW_CONFIG_REPORT(
        file(workflow.projectDir + params.utils.workflow_configuration.report_ipynb)
    )
    out = FILTER_AND_ANNOTATE_AND_CLEAN( data )

    if(params.getToolParams("scanpy").containsKey("filter")) {
        out = QC_FILTER( out ).filtered // Remove concat
    }    
    NORMALIZE_TRANSFORM( out )
    HVG_SELECTION( NORMALIZE_TRANSFORM.out )
    if(params.getToolParams("scanpy").containsKey("regress_out")) {
        preprocessed_data = SC__SCANPY__REGRESS_OUT( HVG_SELECTION.out.scaled )
    } else {
        preprocessed_data = HVG_SELECTION.out.scaled
    }
    DIM_REDUCTION_PCA( preprocessed_data )
    NEIGHBORHOOD_GRAPH( DIM_REDUCTION_PCA.out )
    DIM_REDUCTION_TSNE_UMAP( NEIGHBORHOOD_GRAPH.out )
    CLUSTER_IDENTIFICATION(
        NORMALIZE_TRANSFORM.out,
        DIM_REDUCTION_TSNE_UMAP.out.dimred_tsne_umap,
        "No Batch Effect Correction"
    )

    // Conversion
    // Convert h5ad to X (here we choose: loom format)
    filteredloom = SC__H5AD_TO_FILTERED_LOOM( QC_FILTER.out.filtered )
    scopeloom = FILE_CONVERTER(
        CLUSTER_IDENTIFICATION.out.marker_genes,
        'SINGLE_SAMPLE_STAR.final_output',
        'loom',
        QC_FILTER.out.filtered
    )

    // Define the parameters for clustering
    def clusteringParams = SC__SCANPY__CLUSTERING_PARAMS( clean(params.getToolParams("scanpy").clustering) )

    // Select a default clustering when in parameter exploration mode
    if(params.sc.containsKey("directs") && clusteringParams.isParameterExplorationModeOn()) {
        scopeloom = SC__DIRECTS__SELECT_DEFAULT_CLUSTERING( scopeloom )
    }

    // Publishing
    PUBLISH( 
        CLUSTER_IDENTIFICATION.out.marker_genes.map {
            it -> tuple(it[0], it[1], null)
        },
        "single_sample_star.final_output",
        "h5ad",
        null,
        clusteringParams.isParameterExplorationModeOn()
    )

    /*******************************************
     * Reporting
     */

    ipynbs = QC_FILTER.out.report.map {
        it -> tuple(it[0], it[1])
    }.mix(
        samples.combine(UTILS__GENERATE_WORKFLOW_CONFIG_REPORT.out),
        HVG_SELECTION.out.report.map {
            it -> tuple(it[0], it[1])
        },
        DIM_REDUCTION_TSNE_UMAP.out.report.map {
            it -> tuple(it[0], it[1])
        }
    ).join(
        CLUSTER_IDENTIFICATION.out.report,
        by: 0
    )

    if(!clusteringParams.isParameterExplorationModeOn()) {
        ipynbs = ipynbs.map {
            it -> tuple(it[0], it[1..it.size()-2], null)
        }
    } else {
        ipynbs = ipynbs.map {
            it -> tuple(it[0], it[1..it.size()-2], it[it.size()-1])
        }
    }

    SC__SCANPY__MERGE_REPORTS(
        ipynbs,
        "merged_report",
        clusteringParams.isParameterExplorationModeOn()
    )
    SC__SCANPY__REPORT_TO_HTML(SC__SCANPY__MERGE_REPORTS.out)

    emit:
        filteredloom
        scopeloom

}
