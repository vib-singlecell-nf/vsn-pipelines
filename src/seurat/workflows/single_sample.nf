nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the utils module:
include {
    PUBLISH;
} from '../../utils/workflows/utils.nf' params(params)
include {
    UTILS__GENERATE_WORKFLOW_CONFIG_REPORT;
} from '../../utils/processes/reports.nf' params(params)

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    FILTER;
} from './filter.nf' params(params)
include {
    NORMALIZE;
    NORMALIZE_SCALE_SCT;
} from './normalize.nf' params(params)
include {
    HVG_SELECTION;
} from './hvg_selection.nf' params(params)
include {
    DIM_REDUCTION_PCA;
} from './dim_reduction_pca.nf' params(params)
include {
    DIM_REDUCTION_TSNE_UMAP;
} from './dim_reduction.nf' params(params)
include {
    NEIGHBORHOOD_GRAPH;
} from './neighborhood_graph.nf' params(params)
include {
    CLUSTERING
} from './clustering.nf' params(params)
include {
    DIFFERENTIAL_GENE_EXPRESSION;
} from './deg.nf' params(params)

//////////////////////////////////////////////////////
// Define the workflow

workflow single_sample {
    take:
        data
    
    main:
        // out = FILTER_AND_ANNOTATE_AND_CLEAN( data )
        filtered = params.tools.seurat?.filter ? FILTER( data ).filtered : data
        
        if(params.tools.seurat.normalization.method == "SCT") {
            // Normalize, scale + HVG in 1 step using SCT
            normalized_scaled = NORMALIZE_SCALE_SCT( filtered )
        } else {
            NORMALIZE( filtered ) | HVG_SELECTION
            normalized_scaled = HVG_SELECTION.out.scaled
        }
        
        DIM_REDUCTION_PCA( normalized_scaled )
        NEIGHBORHOOD_GRAPH( DIM_REDUCTION_PCA.out )
        CLUSTERING( NEIGHBORHOOD_GRAPH.out )
        DIM_REDUCTION_TSNE_UMAP( CLUSTERING.out )
        DIFFERENTIAL_GENE_EXPRESSION( CLUSTERING.out )

        UTILS__GENERATE_WORKFLOW_CONFIG_REPORT(
            file(workflow.projectDir + params.utils.workflow_configuration.report_ipynb)
        )

        PUBLISH(
            DIM_REDUCTION_TSNE_UMAP.out.dimred_tsne_umap,
            params.global.project_name+".single_sample_seurat.test_object",
            "Rds",
            'seurat',
            false
        )

    emit:
        data
}