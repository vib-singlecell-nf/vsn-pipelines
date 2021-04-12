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
class SC__SEURAT__DIM_REDUCTION_PARAMS {
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

    String getNCompsAsArgument(nComps) {
        if (!this.env.isParamNull(nComps))
            nComps = nComps.replaceAll("\n", "")
        if (!this.env.isParamNull(nComps) && this.configParams.containsKey('nComps'))
            throw new Exception("SC__SEURAT__DIM_REDUCTION: nComps is both statically (" + this.configParams["nComps"] + ") and dynamically (" + nComps + ") set. Choose one.")
        if (!this.env.isParamNull(nComps))
            return '--n-comps ' + nComps

        return this.configParams.containsKey('nComps') ? '--n-comps ' + this.configParams.nComps : ''
    }

    String getNPcsAsArgument(nPcs) {
        if (!this.env.isParamNull(nPcs))
            nPcs = nPcs.replaceAll("\n","")
        if (!this.env.isParamNull(nPcs) && this.configParams.containsKey('nPcs'))
            throw new Exception("SC__SEURAT__DIM_REDUCTION: nPcs is both statically (" + this.configParams["nPcs"] + ") and dynamically (" + nPcs + ") set for method (" + this.configParams["method"] + "). Choose one.")
        if (!this.env.isParamNull(nPcs))
            return '--n-pcs ' + nPcs

        return this.configParams.containsKey('nPcs') ? '--n-pcs ' + this.configParams.nPcs : ''
    }
}

process SC__SEURAT__DIM_REDUCTION {

    container params.tools.seurat.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__default'

    input:
        tuple \
            val(sampleId), \
            path(f), \
            val(nComps)

    output:
        tuple \
            val(sampleId), \
            path("${sampleId}.SC__SEURAT__DIM_REDUCTION_${method}.${processParams.off}"), \
            val(nComps)

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.seurat.dim_reduction.get(params.method))
        processParams = sampleParams.local
        method = processParams.method.replaceAll('-','').toUpperCase()

        def _processParams = new SC__SEURAT__DIM_REDUCTION_PARAMS()
        _processParams.setEnv(this)
        _processParams.setParams(params)
        _processParams.setConfigParams(processParams)
        """
        ${binDir}/dim_reduction/sc_dim_reduction.R \
            --input $f \
            --output ${sampleId}.SC__SEURAT__DIM_REDUCTION_${method}.${processParams.off} \
            --method ${processParams.method} \
            --seed ${params.global.seed} \
            ${(processParams.containsKey('algorithm')) ? '--algorithm ' + processParams.algorithm : ''} \
            ${_processParams.getNCompsAsArgument(nComps)} \
            ${_processParams.getNPcsAsArgument(nComps)} \
            ${(processParams.containsKey('nNeighbors')) ? '--n-neighbors ' + processParams.nNeighbors : ''}
        """
}
