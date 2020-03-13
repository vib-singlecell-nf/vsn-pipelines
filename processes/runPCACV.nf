nextflow.preview.dsl=2

if(!params.containsKey("test")) {
    binDir = "${workflow.projectDir}/src/pcacv/bin/"
} else {
    binDir = ""
}

process PCACV__FIND_OPTIMAL_NPCS {
    
    container params.pcacv.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.global.threads} -l walltime=1:00:00 -A ${params.global.qsubaccount}"

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), stdout, emit: optimalNumberPC
        tuple val(sampleId), path("${sampleId}.PCACV__FIND_OPTIMAL_NPCS.*")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.pcacv.find_optimal_npcs)
		processParams = sampleParams.local
        """
        ${binDir}run_pca_cv.R \
            --input-file ${f} \
            ${'--seed ' + params.global.seed} \
            ${(processParams.containsKey('accessor')) ? '--accessor "' + processParams.accessor.replace('$','\\$') + '"': ''} \
            ${(processParams.containsKey('useVariableFeatures')) ? '--use-variable-features ' + processParams.useVariableFeatures: ''} \
            ${(processParams.containsKey('kFold')) ? '--k-fold ' + processParams.libraries: ''} \
			${(processParams.containsKey('fromNPC')) ? '--from-n-pc ' + processParams.fromNPC: ''} \
			${(processParams.containsKey('toNPC')) ? '--to-n-pc ' + processParams.toNPC: ''} \
			${(processParams.containsKey('byNPC')) ? '--by-n-pc ' + processParams.byNPC: ''} \
            ${(processParams.containsKey('maxIters')) ? '--max-iters ' + processParams.libraries: ''} \
			${(processParams.containsKey('nCores')) ? '--n-cores ' + processParams.nCores: ''} \
			${(processParams.containsKey('defaultSVD')) ? '--default-svd ' + processParams.defaultSVD: ''} \
            ${(processParams.containsKey('verbose')) ? '--verbose ' + processParams.verbose: ''} \
            --output-prefix "${sampleId}.PCACV__FIND_OPTIMAL_NPCS" \
            > .command.log 2>&1
        cat "${sampleId}.PCACV__FIND_OPTIMAL_NPCS.OPTIMAL_NPCS.txt"
        """

}
