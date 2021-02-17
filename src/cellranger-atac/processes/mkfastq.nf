nextflow.enable.dsl=2

toolParams = params.tools.cellranger_atac

process SC__CELLRANGER_ATAC__MKFASTQ {

    publishDir "${params.global.outdir}/fastqs", saveAs: { filename -> dirname = filename =~ /(.*)_fastqOut/; "${dirname[0][1]}" }, mode: 'link', overwrite: true
    container toolParams.container
    label 'compute_resources__cellranger_mkfastq'

    input:
        file(csv)
        file(runFolder)

    output:
        file "*_fastqOut"

    script:
        """
        cellranger-atac mkfastq \
            --run=${runFolder} \
            --csv=${csv} \
            ${(toolParams.mkfastq.containsKey('samplesheet')) ? '--samplesheet ' + toolParams.mkfastq.samplesheet: ''} \
            ${(toolParams.mkfastq.containsKey('qc')) ? '--qc ' + toolParams.mkfastq.qc: ''} \
            ${(toolParams.mkfastq.containsKey('lanes')) ? '--lanes ' + toolParams.mkfastq.lanes: ''} \
            ${(toolParams.mkfastq.containsKey('useBasesMask')) ? '--use-bases-mask ' + toolParams.mkfastq.useBasesMask: ''} \
            ${(toolParams.mkfastq.containsKey('deleteUndetermined')) ? '--delete-undetermined ' + toolParams.mkfastq.deleteUndetermined: ''} \
            ${(toolParams.mkfastq.containsKey('outputDir')) ? '--output-dir ' + toolParams.mkfastq.outputDir: ''} \
            ${(toolParams.mkfastq.containsKey('project')) ? '--project ' + toolParams.mkfastq.project: ''} \
            --jobmode=local \
            --localcores=${task.cpus} \
            --localmem=${task.memory.toGiga()}
        
        for sample in \$(tail -n+2 ${csv} | cut -f2 -d','); do
            ln -s ${(params.global.containsKey('outputDir')) ? params.global.outputDir + "*/\${sample}" : "*/outs/fastq_path/*/\${sample}"} \${sample}_fastqOut
        done
        """
}

