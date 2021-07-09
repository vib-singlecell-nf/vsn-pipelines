nextflow.enable.dsl=2

import java.nio.file.Paths

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/pcacv/bin" : Paths.get(workflow.scriptFile.getParent().getParent().toString(), "pcacv/bin")


process PCACV__FIND_OPTIMAL_NPCS {
    
    container params.tools.pcacv.container
    publishDir "${params.global.outdir}/data/pcacv", mode: 'link'
    label 'compute_resources__pcacv'

    input:
        tuple \
            val(sampleId), \
            path(f)

    output:
        tuple \
            val(sampleId), \
            stdout, \
            emit: optimalNumberPC
        tuple \
            val(sampleId), \
            path("${sampleId}.PCACV__FIND_OPTIMAL_NPCS.*"), \
            emit: files

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.pcacv.find_optimal_npcs)
        processParams = sampleParams.local
        """
        export OPENBLAS_NUM_THREADS=1
        ${binDir}/run_pca_cv.R \
            --input-file ${f} \
            --seed ${params.global.seed} \
            ${(processParams.containsKey('accessor')) ? '--accessor "' + processParams.accessor.replace('$','\\$') + '"': ''} \
            ${(processParams.containsKey('useVariableFeatures')) ? '--use-variable-features ' + processParams.useVariableFeatures: ''} \
            ${(processParams.containsKey('kFold')) ? '--k-fold ' + processParams.kFold: ''} \
            ${(processParams.containsKey('fromNPC')) ? '--from-n-pc ' + processParams.fromNPC: ''} \
            ${(processParams.containsKey('toNPC')) ? '--to-n-pc ' + processParams.toNPC: ''} \
            ${(processParams.containsKey('byNPC')) ? '--by-n-pc ' + processParams.byNPC: ''} \
            ${(processParams.containsKey('nPCFallback')) ? '--n-pc-fallback ' + processParams.nPCFallback: ''} \
            ${(processParams.containsKey('maxIters')) ? '--max-iters ' + processParams.maxIters: ''} \
            --n-cores ${task.cpus} \
            ${(processParams.containsKey('defaultSVD')) ? '--default-svd ' + processParams.defaultSVD: ''} \
            ${(processParams.containsKey('verbose')) ? '--verbose ' + processParams.verbose: ''} \
            --output-prefix "${sampleId}.PCACV__FIND_OPTIMAL_NPCS" \
            > .command.log 2>&1
        cat "${sampleId}.PCACV__FIND_OPTIMAL_NPCS.OPTIMAL_NPCS.txt"
        """

}
