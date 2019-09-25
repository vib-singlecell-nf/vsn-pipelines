nextflow.preview.dsl=2

process SC__CELLRANGER__MKFASTQ {

  publishDir "${params.outdir}/fastqs", saveAs: { filename -> dirname = filename =~ /(.*)_fastqOut/; "${dirname[0][1]}" }, mode: 'symlink'
  container params.containers.cellranger

  input:
    file(samplesheet)
    file(runFolder)

  output:
    file "*_fastqOut"

  script:
    """
    cellranger mkfastq \
        --run=${runFolder} \
        --samplesheet=${samplesheet} \
        ${(params.containsKey('runID')) ? '--id ' + params.runID: ''} \
        ${(params.containsKey('csv')) ? '--csv ' + params.csv: ''} \
        ${(params.containsKey('ignoreDualIndex')) ? '--ignore-dual-index ' + params.ignoreDualIndex: ''} \
        ${(params.containsKey('qc')) ? '--qc ' + params.qc: ''} \
        ${(params.containsKey('lanes')) ? '--lanes ' + params.lanes: ''} \
        ${(params.containsKey('useBasesMask')) ? '--use-bases-mask ' + params.useBasesMask: ''} \
        ${(params.containsKey('deleteUndetermined')) ? '--delete-undetermined ' + params.deleteUndetermined: ''} \
        ${(params.containsKey('outputDir')) ? '--output-dir ' + params.outputDir: ''} \
        ${(params.containsKey('project')) ? '--project ' + params.project: ''} \
        ${(params.containsKey('jobMode')) ? '--jobmode ' + params.jobMode: ''} \
        ${(params.containsKey('localCores')) ? '--localcores ' + params.localCores: ''} \
        ${(params.containsKey('localMem')) ? '--localmem ' + params.localMem: ''}
    
    for sample in \$(tail -n+2 ${samplesheet} | cut -f2 -d','); do
        ln -s ${(params.containsKey('outputDir')) ? params.outputDir + "*/\${sample}" : "*/outs/fastq_path/*/\${sample}"} \${sample}_fastqOut
    done
    """
}
