nextflow.enable.dsl=2

process SC__STAR__BUILD_INDEX {

    container params.tools.star.container
    label 'compute_resources__star_build_genome'

    input:
        file(annotation)
        file(genome)

    output:
        file("STAR_index")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.star.build_genome)
		processParams = sampleParams.local
        """
        mkdir STAR_index
        STAR \
            --runThreadN ${task.cpus} \
            --runMode genomeGenerate \
            --genomeDir STAR_index \
            --genomeFastaFiles ${genome} \
            --sjdbGTFfile ${annotation} \
            --sjdbOverhang ${processParams.sjdbOverhang} \
            --genomeSAindexNbases ${processParams.genomeSAindexNbases} # Suggested by STAR (default: 14), otherwise keeps on hanging
        """

}
