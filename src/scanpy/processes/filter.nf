nextflow.preview.dsl=2

if(!params.containsKey("test")) {
  	binDir = "${workflow.projectDir}/src/scanpy/bin/"
} else {
  	binDir = ""
}

process SC__SCANPY__COMPUTE_QC_STATS {

  	container params.sc.scanpy.container
  	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"

  	input:
    tuple val(id), file(f)
  	
	output:
    tuple val(id), file("${id}.SC__SCANPY__COMPUTE_QC_STATS.${processParams.off}")
  	
	script:
	processParams = params.sc.scanpy.filter
    """
    ${binDir}filter/sc_cell_gene_filtering.py \
    	compute \
      	$f \
		${id}.SC__SCANPY__COMPUTE_QC_STATS.${processParams.off} \
		${(processParams.containsKey('cellFilterMinNCounts')) ? '--min-n-counts ' + processParams.cellFilterMinNCounts : ''} \
		${(processParams.containsKey('cellFilterMaxNCounts')) ? '--max-n-counts ' + processParams.cellFilterMaxNCounts : ''} \
		${(processParams.containsKey('cellFilterMinNGenes')) ? '--min-n-genes ' + processParams.cellFilterMinNGenes : ''} \
		${(processParams.containsKey('cellFilterMaxNGenes')) ? '--max-n-genes ' + processParams.cellFilterMaxNGenes : ''} \
		${(processParams.containsKey('cellFilterMaxPercentMito')) ? '--max-percent-mito ' + processParams.cellFilterMaxPercentMito : ''} \
		${(processParams.containsKey('geneFilterMinNCells')) ? '--min-number-cells ' + processParams.geneFilterMinNCells : ''}
    """
}


process SC__SCANPY__GENE_FILTER {

    container params.sc.scanpy.container
    clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true

    input:
    tuple val(id), file(f)
    
	output:
    tuple val(id), file("${id}.SC__SCANPY__GENE_FILTER.${processParams.off}")
    
	script:
	processParams = params.sc.scanpy.filter
    """
    ${binDir}filter/sc_cell_gene_filtering.py \
        genefilter \
        $f \
        ${id}.SC__SCANPY__GENE_FILTER.${processParams.off} \
        ${(processParams.containsKey('geneFilterMinNCells')) ? '--min-number-cells ' + processParams.geneFilterMinNCells : ''}
    """
}


process SC__SCANPY__CELL_FILTER {

    container params.sc.scanpy.container
    clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true

    input:
    tuple val(id), file(f)
    
	output:
    tuple val(id), file("${id}.SC__SCANPY__CELL_FILTER.${processParams.off}")
    
	script:
	processParams = params.sc.scanpy.filter
    """
    ${binDir}filter/sc_cell_gene_filtering.py \
        cellfilter \
        $f \
        ${id}.SC__SCANPY__CELL_FILTER.${processParams.off} \
        ${(processParams.containsKey('cellFilterMinNCounts')) ? '--min-n-counts ' + processParams.cellFilterMinNCounts : ''} \
        ${(processParams.containsKey('cellFilterMaxNCounts')) ? '--max-n-counts ' + processParams.cellFilterMaxNCounts : ''} \
        ${(processParams.containsKey('cellFilterMinNGenes')) ? '--min-n-genes ' + processParams.cellFilterMinNGenes : ''} \
        ${(processParams.containsKey('cellFilterMaxNGenes')) ? '--max-n-genes ' + processParams.cellFilterMaxNGenes : ''} \
        ${(processParams.containsKey('cellFilterMaxPercentMito')) ? '--max-percent-mito ' + processParams.cellFilterMaxPercentMito : ''}
    """
}
