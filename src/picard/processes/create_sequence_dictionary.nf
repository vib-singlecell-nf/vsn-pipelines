nextflow.enable.dsl=2

process PICARD__CREATE_SEQUENCE_DICTIONARY {

    container params.getToolParams("picard").container
    publishDir "${params.global.outdir}/00.refdata", mode: 'symlink'
    label 'compute_resources__default'

    input:
        file(genome)
        file(tmpDir)
    
    output:
        file "${genome.baseName}.dict"

    script:
        """
        java -Djava.io.tmpdir=$tmpDir -jar \
            /picard.jar \
                CreateSequenceDictionary \
                    R=${genome} \
                    O=${genome.baseName}.dict
        """

}
