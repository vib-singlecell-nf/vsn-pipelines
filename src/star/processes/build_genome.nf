nextflow.preview.dsl=2

process SC__STAR__BUILD_INDEX {

    container params.sc.star.container

    clusterOptions "-l nodes=1:ppn=${params.runThreadN} -l walltime=1:00:00 -A ${params.qsubaccount}"

    input:
        file(annotation)
        file(genome)
    output:
        file("STAR_index")
    script:
        """
        mkdir STAR_index
        STAR \
            --runThreadN ${params.runThreadN} \
            --runMode genomeGenerate \
            --genomeDir STAR_index \
            --genomeFastaFiles ${genome} \
            --sjdbGTFfile ${annotation} \
            --sjdbOverhang ${params.sjdbOverhang} \
            --genomeSAindexNbases ${params.genomeSAindexNbases} # Suggested by STAR (default: 14), otherwise keeps on hanging
        """
}

