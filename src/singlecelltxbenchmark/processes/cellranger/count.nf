nextflow.preview.dsl=2

process SC__CELLRANGER__COUNT {

  publishDir "${params.outdir}/counts", mode: 'symlink'
  container params.container

  input:
    file(transcriptome)
    file(fastqs)

  output:
    tuple(_sampleName, file("${_sampleName}/outs") )

  script:
    sample = fastqs.getName() =~ /(.*)_fastqOut/
    _sampleName = sample[0][1]
    """
    cellranger count \
        --id=${_sampleName} \
        --sample=${_sampleName} \
        --fastqs=${fastqs} \
        --transcriptome=${transcriptome} \
        ${(params.containsKey('libraries')) ? '--libraries ' + params.libraries: ''} \
        ${(params.containsKey('featureRef')) ? '--feature-ref ' + params.featureRef: ''} \
        ${(params.containsKey('expectCells')) ? '--expect-cells ' + params.expectCells: ''} \
        ${(params.containsKey('forceCells')) ? '--force-cells ' + params.forceCells: ''} \
        ${(params.containsKey('nosecondary')) ? '--nosecondary ' + params.nosecondary: ''} \
        ${(params.containsKey('noLibraries')) ? '--no-libraries ' + params.noLibraries: ''} \
        ${(params.containsKey('chemistry')) ? '--chemistry ' + params.chemistry: ''} \
        ${(params.containsKey('r1Length')) ? '--r1-length ' + params.r1Length: ''} \
        ${(params.containsKey('r2Length')) ? '--r2-length ' + params.r2Length: ''} \
        ${(params.containsKey('lanes')) ? '--lanes ' + params.lanes: ''} \
        ${(params.containsKey('localCores')) ? '--localcores ' + params.localCores: ''} \
        ${(params.containsKey('localMem')) ? '--localmem ' + params.localMem: ''} \
        ${(params.containsKey('indicies')) ? '--indicies ' + params.indicies: ''} 
    """
}
