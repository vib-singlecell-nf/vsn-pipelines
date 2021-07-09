nextflow.enable.dsl=2

import java.nio.file.Paths
import groovy.transform.TupleConstructor

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scrublet/bin" : Paths.get(workflow.scriptFile.getParent().getParent().toString(), "scrublet/bin")

include {
	isParamNull;
} from '../../utils/processes/utils.nf'

@TupleConstructor()
class SC__SCRUBLET__DOUBLET_DETECTION_PARAMS {

	Script env = null;
	Map params = null;
	LinkedHashMap configParams = null;
	// Parameters definiton
	String off = null;

	void setEnv(env) {
		this.env = env
	}

	void setParams(params) {
		this.params = params
	}

	void setConfigProcessParams(params) {
		this.configProcessParams = params
	}

	String getNPrinCompsAsArgument(nPrinComps) {
		// Check if nComps is both dynamically and if statically set
		if(!this.env.isParamNull(nPrinComps) && this.configParams.containsKey('nPrinComps'))
			throw new Exception("SC__SCRUBLET__DOUBLET_DETECTION: nPrinComps is both statically (" + this.configParams["nPrinComps"] + ") and dynamically (" + nPrinComps + ") set. Choose one.")
		if(!this.env.isParamNull(nPrinComps))
			return '--n-prin-comps ' + nPrinComps.replaceAll("\n","")
		return this.configParams.containsKey('nPrinComps') ? '--n-prin-comps ' + this.configParams.nPrinComps: ''
	}

}

def SC__SCRUBLET__DOUBLET_DETECTION_PARAMS(params) {
	return (new SC__SCRUBLET__DOUBLET_DETECTION_PARAMS(params))
}

process SC__SCRUBLET__DOUBLET_DETECTION {

	container params.tools.scrublet.container
	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__mem'

	input:
		tuple \
			val(sampleId), \
			path(adataRaw), \
            path(adataWithHvgInfo), \
			val(stashedParams), \
			val(nPrinComps)

	output:
		tuple \
			val(sampleId), \
			path("${sampleId}.SC__SCRUBLET__DOUBLET_DETECTION.ScrubletDoubletTable.tsv"), \
			path("${sampleId}.SC__SCRUBLET__DOUBLET_DETECTION.ScrubletObject.pklz"), \
			val(stashedParams), \
			val(nPrinComps)

	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.tools.scrublet.doublet_detection)
		processParams = sampleParams.local
		def _processParams = new SC__SCRUBLET__DOUBLET_DETECTION_PARAMS()
		_processParams.setEnv(this)
		_processParams.setParams(params)
		_processParams.setConfigParams(processParams)
		"""
		${binDir}/sc_doublet_detection.py \
			${processParams?.useVariableFeatures ? '--use-variable-features ' + processParams.useVariableFeatures : ''} \
            ${processParams?.syntheticDoubletUmiSubsampling ? '--synthetic-doublet-umi-subsampling ' + processParams.syntheticDoubletUmiSubsampling : ''} \
            ${processParams?.minCounts ? '--min-counts ' + processParams.minCounts : ''} \
            ${processParams?.minCells ? '--min-cells ' + processParams.minCells : ''} \
            ${processParams?.minGeneVariabilityPctl ? '--min-gene-variability-pctl ' + processParams.minGeneVariabilityPctl : ''} \
            ${processParams?.logTransform ? '--log-transform ' + processParams.logTransform : ''} \
            ${processParams?.meanCenter ? '--mean-center ' + processParams.meanCenter : ''} \
            ${processParams?.normalizeVariance ? '--normalize-variance ' + processParams.normalizeVariance : ''} \
            ${_processParams.getNPrinCompsAsArgument(nPrinComps)} \
            ${processParams?.technology ? '--technology ' + processParams.technology : ''} \
			${processParams?.useVariableFeatures ? '--h5ad-with-variable-features-info ' + adataWithHvgInfo : ''} \
			${processParams?.threshold ? '--threshold ' + processParams.threshold : ''} \
			--output-prefix "${sampleId}.SC__SCRUBLET__DOUBLET_DETECTION" \
			$adataRaw
		"""

}
