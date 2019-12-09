nextflow.preview.dsl=2

process SC__CELLRANGER__MKFASTQ {

	publishDir "${params.global.outdir}/fastqs", saveAs: { filename -> dirname = filename =~ /(.*)_fastqOut/; "${dirname[0][1]}" }, mode: 'link', overwrite: true
  	container params.sc.cellranger.container

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
			${(params.sc.cellranger.mkfastq.containsKey('runID')) ? '--id ' + params.sc.cellranger.mkfastq.runID: ''} \
			${(params.sc.cellranger.mkfastq.containsKey('samplesheet')) ? '--samplesheet ' + params.sc.cellranger.mkfastq.samplesheet: ''} \
			${(params.sc.cellranger.mkfastq.containsKey('ignoreDualIndex')) ? '--ignore-dual-index ' + params.sc.cellranger.mkfastq.ignoreDualIndex: ''} \
			${(params.sc.cellranger.mkfastq.containsKey('qc')) ? '--qc ' + params.sc.cellranger.mkfastq.qc: ''} \
			${(params.sc.cellranger.mkfastq.containsKey('lanes')) ? '--lanes ' + params.sc.cellranger.mkfastq.lanes: ''} \
			${(params.sc.cellranger.mkfastq.containsKey('useBasesMask')) ? '--use-bases-mask ' + params.sc.cellranger.mkfastq.useBasesMask: ''} \
			${(params.sc.cellranger.mkfastq.containsKey('deleteUndetermined')) ? '--delete-undetermined ' + params.sc.cellranger.mkfastq.deleteUndetermined: ''} \
			${(params.sc.cellranger.mkfastq.containsKey('outputDir')) ? '--output-dir ' + params.sc.cellranger.mkfastq.outputDir: ''} \
			${(params.sc.cellranger.mkfastq.containsKey('project')) ? '--project ' + params.sc.cellranger.mkfastq.project: ''} \
			${(params.sc.cellranger.mkfastq.containsKey('jobMode')) ? '--jobmode ' + params.sc.cellranger.mkfastq.jobMode: ''} \
			${(params.sc.cellranger.mkfastq.containsKey('localCores')) ? '--localcores ' + params.sc.cellranger.mkfastq.localCores: ''} \
			${(params.sc.cellranger.mkfastq.containsKey('localMem')) ? '--localmem ' + params.sc.cellranger.mkfastq.localMem: ''}
		
		for sample in \$(tail -n+2 ${csv} | cut -f2 -d','); do
			ln -s ${(params.global.containsKey('outputDir')) ? params.global.outputDir + "*/\${sample}" : "*/outs/fastq_path/*/\${sample}"} \${sample}_fastqOut
		done
		"""

}
