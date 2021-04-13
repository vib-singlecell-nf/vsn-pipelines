nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    SC__SEURAT__MARKER_GENES;
} from '../processes/marker_genes.nf' params(params)
include {
    SC__SEURAT__MARKER_GENES_TO_XLSX;
} from '../processes/utils.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow DIFFERENTIAL_GENE_EXPRESSION {

    take:
        data
    
    main:
        SC__SEURAT__MARKER_GENES( data ) | SC__SEURAT__MARKER_GENES_TO_XLSX

    emit:
        marker_genes = SC__SEURAT__MARKER_GENES.out
        marker_genex_xlsx = SC__SEURAT__MARKER_GENES_TO_XLSX.out
}