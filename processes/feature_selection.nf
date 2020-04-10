nextflow.preview.dsl=2

import java.nio.file.Paths

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scanpy/bin" : Paths.get(workflow.scriptFile.getParent().toString(), "bin")

process SC__SCANPY__FIND_HIGHLY_VARIABLE_GENES {

  	container params.sc.scanpy.container
  	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true

  	input:
    	tuple val(sampleId), path(f)

  	output:
    	tuple val(sampleId), path("${sampleId}.SC__SCANPY__FIND_HIGHLY_VARIABLE_GENES.${processParams.off}")

  	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scanpy.feature_selection)
		processParams = sampleParams.local
		"""
		${binDir}/feature_selection/sc_find_variable_genes.py \
			--method ${processParams.method} \
			${(processParams.containsKey('minMean')) ? '--min-mean ' + processParams.minMean : ''} \
			${(processParams.containsKey('maxMean')) ? '--max-mean ' + processParams.maxMean : ''} \
			${(processParams.containsKey('minDisp')) ? '--min-disp ' + processParams.minDisp : ''} \
			${(processParams.containsKey('maxDisp')) ? '--max-disp ' + processParams.maxDisp : ''} \
			$f \
			"${sampleId}.SC__SCANPY__FIND_HIGHLY_VARIABLE_GENES.${processParams.off}"
		"""

}

process SC__SCANPY__SUBSET_HIGHLY_VARIABLE_GENES {

  	container params.sc.scanpy.container
  	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true

  	input:
    	tuple val(sampleId), path(f)

  	output:
    	tuple val(sampleId), path("${sampleId}.SC__SCANPY__SUBSET_HIGHLY_VARIABLE_GENES.${processParams.off}")

  	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scanpy.feature_selection)
		processParams = sampleParams.local
		"""
		${binDir}/feature_selection/sc_subset_variable_genes.py \
			$f \
			"${sampleId}.SC__SCANPY__SUBSET_HIGHLY_VARIABLE_GENES.${processParams.off}"
		"""

}

