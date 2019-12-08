nextflow.preview.dsl=2

process SC__STAR__BUILD_INDEX {

    container params.sc.star.container
    clusterOptions "-l nodes=1:ppn=${processParams.runThreadN} -l walltime=1:00:00 -A ${params.global.qsubaccount}"

    input:
        file(annotation)
        file(genome)

    output:
        file("STAR_index")

    script:
        processParams = params.sc.star.build_genome
        """
        mkdir STAR_index
        STAR \
            --runThreadN ${processParams.runThreadN} \
            --runMode genomeGenerate \
            --genomeDir STAR_index \
            --genomeFastaFiles ${genome} \
            --sjdbGTFfile ${annotation} \
            --sjdbOverhang ${processParams.sjdbOverhang} \
            --genomeSAindexNbases ${processParams.genomeSAindexNbases} # Suggested by STAR (default: 14), otherwise keeps on hanging
        """

}
