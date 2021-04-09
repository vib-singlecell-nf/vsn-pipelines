nextflow.enable.dsl=2

import java.nio.file.Paths
import groovy.transform.TupleConstructor
import groovyx.gpars.dataflow.DataflowBroadcast
import nextflow.util.ArrayTuple
import nextflow.script.ScriptBinding

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/seurat/bin" : Paths.get(workflow.scriptFile.getParent().toString(), "bin")

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
        if (!this.env.isParamNull(nComps) && this.configParams.containsKey('nComps'))
            throw new Exception("SC__SCANPY__DIM_REDUCTION: nComps is both statically (" + this.configParams["nComps"] + ") and dynamically (" + nComps + ") set. Choose one.")
        if (!this.env.isParamNull(nComps))
            return '--n-comps ' + nComps.replaceAll("\n","")

    return this.configParams.containsKey('nComps') ? '--n-comps ' + this.configParams.nComps : ''
    }
}

process SC__SEURAT__DIM_REDUCTION {

    container params.tools.seurat.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__default'

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), path("${sampleId}.SC__SEURAT__DIM_REDUCTION_${method}.${processParams.off}")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.seurat.dim_reduction.get(params.method))
        processParams = sampleParams.local
        method = processParams.method.replaceAll('-','').toUpperCase()

        // Cannot call constructor with parameter if nComps is not provided (aka NULL), type do not match
        // def _processParams = new SC__SEURAT__DIM_REDUCTION_PARAMS()
        // _processParams.setEnv(this)
        // _processParams.setParams(params)
        // _processParams.setConfigParams(processParams)
        """
        ${binDir}/dim_reduction/sc_dim_reduction.R \
            --input $f \
            --output ${sampleId}.SC__SEURAT__DIM_REDUCTION_${method}.${processParams.off} \
            --method ${processParams.method} \
            --seed ${params.global.seed} \
            ${(processParams.containsKey('algorithm')) ? '--algorithm ' + processParams.algorithm : ''} \
            ${(processParams.containsKey('nComps')) ? '--n-comps ' + processParams.nComps : ''} \
            ${(processParams.containsKey('nPcs')) ? '--n-pcs ' + processParams.nPcs : ''} \
            ${(processParams.containsKey('nNeighbors')) ? '--n-neighbors ' + processParams.nNeighbors : ''}
        """
}
