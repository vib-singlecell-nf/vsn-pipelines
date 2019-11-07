nextflow.preview.dsl=2

process PICARD__CREATE_SEQUENCE_DICTIONARY {

    container params.picard.container
    publishDir "${params.outdir}/00.refdata", mode: 'symlink'
    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=1:00:00 -A ${params.qsubaccount}"

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