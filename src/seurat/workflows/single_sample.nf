nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the utils module:
include {
    PUBLISH;
} from '../../utils/workflows/utils.nf' params(params)
include {
    FILE_CONVERTER as FILE_CONVERTER_TO_SCOPE;
} from '../../utils/workflows/fileConverter.nf' params(params)
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
        filtered = params.tools.seurat?.filter ? FILTER( data ).filtered : data
        
        if (params.containsKey("sct")) {
            // Normalize, scale + HVG in 1 step using SCT
            NORMALIZE_SCALE_SCT( filtered )
            normalized_scaled = NORMALIZE_SCALE_SCT.out.scaled
        } else {
            NORMALIZE( filtered ) | HVG_SELECTION
            normalized_scaled = HVG_SELECTION.out.scaled
        }
        
        DIM_REDUCTION_PCA( normalized_scaled )
        NEIGHBORHOOD_GRAPH( DIM_REDUCTION_PCA.out )
        DIM_REDUCTION_TSNE_UMAP( NEIGHBORHOOD_GRAPH.out )
        CLUSTERING( DIM_REDUCTION_TSNE_UMAP.out.dimred_tsne_umap )
        DIFFERENTIAL_GENE_EXPRESSION( CLUSTERING.out.clustered )

        UTILS__GENERATE_WORKFLOW_CONFIG_REPORT(
            file(workflow.projectDir + params.utils.workflow_configuration.report_ipynb)
        )

        PUBLISH(
            CLUSTERING.out.clustered,
            params.global.project_name+".single_sample_seurat.final_output",
            "Rds",
            'seurat',
            false
        )

        FILE_CONVERTER_TO_SCOPE(
            CLUSTERING.out.clustered,
            'SINGLE_SAMPLE_SEURAT.final_output',
            'seuratRdsToSCopeLoom',
            null
        )

    emit:
        filtered_seurat_rds = filtered
        scope_loom = FILE_CONVERTER_TO_SCOPE.out
        seurat_rds = CLUSTERING.out.clustered
        marker_genes = DIFFERENTIAL_GENE_EXPRESSION.out.marker_genes
        marker_genes_xlsx = DIFFERENTIAL_GENE_EXPRESSION.out.marker_genex_xlsx

}