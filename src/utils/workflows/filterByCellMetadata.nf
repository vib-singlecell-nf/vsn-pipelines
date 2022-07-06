nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  Process imports:
include {
    isParamNull;
    getToolParams;
} from './../processes/utils.nf' params(params)
include {
    SC__PREPARE_OBS_FILTER;
} from './../processes/h5adSubset' params(params)
include {
    SC__APPLY_OBS_FILTER;
} from './../processes/h5adSubset' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow FILTER_BY_CELL_METADATA {

    take:
        // Expects (sampleId, h5ad) : Channel
        data
        // Describes: name of tool
        // Expects tool: (string || null)
        // Values
        // - tool != null:
        //   - The given tool is performing itself a cell-based filtering
        //   - params.tools[tool] should exist
        // - tool == null:
        //   - params.utils.cell_filter should exist
        tool

    main:
        def workflowParams = isParamNull(tool) ?
            params.utils.cell_filter :
            getToolParams(params.tools, tool)["cell_filter"]
        Channel
            .from(workflowParams.filters)
            .set{ filters }
        SC__PREPARE_OBS_FILTER(
            data.map {
                it -> tuple(it[0], it[1])
            }.combine(filters),
            isParamNull(tool) ? 'NULL' : tool
        )
        out = SC__APPLY_OBS_FILTER(
            SC__PREPARE_OBS_FILTER.out.groupTuple(),
            isParamNull(tool) ? 'NULL' : tool
        )

    emit:
        out

}
