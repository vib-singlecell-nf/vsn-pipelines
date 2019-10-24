nextflow.preview.dsl=2

process SC__STAR__BUILD_INDEX {

    clusterOptions "-l nodes=1:ppn=${params.threads} -l walltime=1:00:00 -A ${params.qsubaccount}"

    input:
        file(annotation)
        file(genome)
    output:
        file("STAR_index")
    script:
        """
        mkdir STAR_index
        singularity run \
            -B /ddn1/vol1/staging/leuven/stg_00002/:/ddn1/vol1/staging/leuven/stg_00002/ \
            -B /staging/leuven/stg_00002/:/staging/leuven/stg_00002/ \
            -B /ddn1/vol1/staging/leuven/res_00001/genomes/:/ddn1/vol1/staging/leuven/res_00001/genomes/ \
            -B /staging/leuven/res_00001/genomes/:/staging/leuven/res_00001/genomes/ \
                /staging/leuven/res_00001/software/STAR/2.7.1a/STAR_2.7.1a.sif \
                --runThreadN ${params.threads} \
                --runMode genomeGenerate \
                --genomeDir STAR_index \
                --genomeFastaFiles ${genome} \
                --sjdbGTFfile ${annotation} \
                --sjdbOverhang ${params.sjdbOverhang} \
                --genomeSAindexNbases ${params.genomeSAindexNbases} # Suggested by STAR (default: 14), otherwise keeps on hanging
        """
}

