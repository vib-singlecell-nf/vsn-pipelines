nextflow.preview.dsl=2

process PICARD__CREATE_SEQUENCE_DICTIONARY {
    publishDir "${params.outdir}/00.refdata", mode: 'symlink'

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=1:00:00 -A ${params.qsubaccount}"

    input:
        file(genome)
    output:
        file "${genome.baseName}.dict"
    script:
        """
        module load Java/1.8.0_192
        java -Djava.io.tmpdir=$DWMAX/tmp -jar \
            /ddn1/vol1/staging/leuven/stg_00002/lcb/dwmax/software/genius/picard/2.20.6/bin/picard.jar \
                CreateSequenceDictionary \
                    R=${genome} \
                    O=${genome.baseName}.dict
        """
}