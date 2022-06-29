nextflow.enable.dsl=2

import java.nio.file.Paths
import groovy.transform.TupleConstructor
import groovyx.gpars.dataflow.DataflowBroadcast
import nextflow.util.ArrayTuple
import nextflow.script.ScriptBinding

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/seurat/bin" : Paths.get(workflow.scriptFile.getParent().toString(), "bin")

include {
	isParamNull;
} from '../../utils/processes/utils.nf'

@TupleConstructor
class SC__SEURAT__NEIGHBORHOOD_GRAPH_PARAMS {
    Script env = null;
    Map params = null;
    LinkedHashMap configParams = null;

    void setEnv(env) {
        this.env = env
    }

    void setParams(params) {
        this.params = params
    }

    void setConfigProcessParams(params) {
        this.configProcessParams = params
    }

    String getNPcsAsArgument(nPcs) {
        if (!this.env.isParamNull(nPcs) && this.configParams.containsKey('nPcs'))
            throw new Exception("SC__SEURAT__NEIGHBORHOOD_GRAPH: nPcs is both statically (" + this.configParams["nPcs"] + ") and dynamically (" + nPcs + ") set. Choose one.")
        if (!this.env.isParamNull(nPcs))
            return '--n-pcs ' + nPcs.replaceAll("\n","")

        return this.configParams.containsKey('nPcs') ? '--n-pcs ' + this.configParams.nPcs : ''
    }
}

process SC__SEURAT__NEIGHBORHOOD_GRAPH {

    container params.tools.seurat.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__default'

    input:
        tuple \
            val(sampleId), \
            path(f), \
            val(nPcs)
    
    output:
        tuple \
            val(sampleId), \
            path("${sampleId}.SC__SEURAT__NEIGHBORHOOD_GRAPH.${processParams.off}"), \
            val(nPcs)

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.seurat.neighborhood_graph)
        processParams = sampleParams.local

        def _processParams = new SC__SEURAT__NEIGHBORHOOD_GRAPH_PARAMS()
        _processParams.setEnv(this)
        _processParams.setParams(params)
        _processParams.setConfigParams(processParams)
        """
        ${binDir}/nn/sc_neighborhood_graph.R \
            --input $f \
            --output ${sampleId}.SC__SEURAT__NEIGHBORHOOD_GRAPH.${processParams.off} \
            ${_processParams.getNPcsAsArgument(nPcs)}
        """
}
