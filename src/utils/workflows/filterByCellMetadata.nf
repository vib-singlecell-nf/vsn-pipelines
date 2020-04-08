nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  process imports:
include './../processes/utils.nf' params(params)
include SC__PREPARE_OBS_FILTER from './../processes/h5adSubset' params(params)
include SC__APPLY_OBS_FILTER from './../processes/h5adSubset' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow FILTER_BY_CELL_METADATA {

    take:
        // Expects (sampleId, h5ad)
        data
        // Expects name of tool ([string] || null)
        tool

    main:
        def workflowParams = isParamNull(tool) ? params.sc.cell_filter : params.sc[tool].cell_filter
        Channel
            .from(workflowParams.filters)
            .set{ filters }
        SC__PREPARE_OBS_FILTER(
            data.combine(filters),
            isParamNull(tool) ? 'NULL' : tool
        )
        out = SC__APPLY_OBS_FILTER(
            SC__PREPARE_OBS_FILTER.out.groupTuple(),
            isParamNull(tool) ? 'NULL' : tool
        )

    emit:
        out

}
