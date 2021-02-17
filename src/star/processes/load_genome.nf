nextflow.enable.dsl=2

process SC__STAR__LOAD_GENOME {

  	container params.getToolParams("star").container
    label 'compute_resources__default'

	input:
		file(starIndex)

	output:
		val starIndexLoaded 

	script:
		starIndexLoaded = true
		"""
		STAR \
			--genomeLoad LoadAndExit \
			--genomeDir ${starIndex}
		"""

}
