nextflow.preview.dsl=2

toolParams = params.sc.cellranger

process SC__CELLRANGER__MKFASTQ {

	publishDir "${params.global.outdir}/fastqs", saveAs: { filename -> dirname = filename =~ /(.*)_fastqOut/; "${dirname[0][1]}" }, mode: 'link', overwrite: true
  	container toolParams.container

  	input:
		file(csv)
    	file(runFolder)

  	output:
    	file "*_fastqOut"

  	script:
		"""
		cellranger mkfastq \
			--run=${runFolder} \
			--csv=${csv} \
			${(toolParams.mkfastq.containsKey('runID')) ? '--id ' + toolParams.mkfastq.runID: ''} \
			${(toolParams.mkfastq.containsKey('samplesheet')) ? '--samplesheet ' + toolParams.mkfastq.samplesheet: ''} \
			${(toolParams.mkfastq.containsKey('ignoreDualIndex')) ? '--ignore-dual-index ' + toolParams.mkfastq.ignoreDualIndex: ''} \
			${(toolParams.mkfastq.containsKey('qc')) ? '--qc ' + toolParams.mkfastq.qc: ''} \
			${(toolParams.mkfastq.containsKey('lanes')) ? '--lanes ' + toolParams.mkfastq.lanes: ''} \
			${(toolParams.mkfastq.containsKey('useBasesMask')) ? '--use-bases-mask ' + toolParams.mkfastq.useBasesMask: ''} \
			${(toolParams.mkfastq.containsKey('deleteUndetermined')) ? '--delete-undetermined ' + toolParams.mkfastq.deleteUndetermined: ''} \
			${(toolParams.mkfastq.containsKey('outputDir')) ? '--output-dir ' + toolParams.mkfastq.outputDir: ''} \
			${(toolParams.mkfastq.containsKey('project')) ? '--project ' + toolParams.mkfastq.project: ''} \
			${(toolParams.mkfastq.containsKey('jobMode')) ? '--jobmode ' + toolParams.mkfastq.jobMode: ''} \
			${(toolParams.mkfastq.containsKey('localCores')) ? '--localcores ' + toolParams.mkfastq.localCores: ''} \
			${(toolParams.mkfastq.containsKey('localMem')) ? '--localmem ' + toolParams.mkfastq.localMem: ''}
		
		for sample in \$(tail -n+2 ${csv} | cut -f2 -d','); do
			ln -s ${(params.global.containsKey('outputDir')) ? params.global.outputDir + "*/\${sample}" : "*/outs/fastq_path/*/\${sample}"} \${sample}_fastqOut
		done
		"""

}
