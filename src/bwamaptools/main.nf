nextflow.enable.dsl=2

import java.nio.file.Paths

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include {
    BWAMAPTOOLS__BWA_MEM_PE as BWA_MEM_PE;
} from './processes/mapping.nf' params(params)
include {
    BWAMAPTOOLS__INDEX_BAM as INDEX_BAM;
} from './processes/index.nf' params(params)
include {
    BWAMAPTOOLS__MAPPING_SUMMARY as MAPPING_SUMMARY;
} from './processes/mapping_summary.nf' params(params)
include {
    PUBLISH as PUBLISH_BAM;
    PUBLISH as PUBLISH_BAM_INDEX;
    PUBLISH as PUBLISH_MAPPING_SUMMARY;
} from "../utils/workflows/utils.nf" params(params)

//////////////////////////////////////////////////////
// Define the workflow

workflow get_bwa_index {

    take:
        fasta_path

    main:

        bwa_fasta = Channel.fromPath(fasta_path)

        bwa_index_path = Paths.get(
                                   Paths.get(fasta_path).getParent().toString(),
                                   "*.{amb,ann,bwt,fai,flat,gdx,pac,sa,0123,bwt.2bit.64}"
                                   )
        bwa_index = Channel.fromPath(bwa_index_path,
                                     glob: true,
                                     type: 'file',
                                     )
                           .ifEmpty { exit 1, "ERROR: Could not find bwa indices from: ${bwa_index_path}." }
                           .collect()
                           .toList()

        data_channel = bwa_fasta.combine(bwa_index)

    emit:
        data_channel

}


workflow BWA_MAPPING_PE {

    take:
        data // a channel of [val(sampleId), path(fastq_PE1), path(fastq_PE2)]

    main:
        /* 
           1) create a channel linking bwa index files from genome.fa in params, and
           2) combine this channel with the items in the data channel
        */
        bwa_inputs = get_bwa_index(params.tools.bwamaptools.bwa_fasta).combine(data)

        BWA_MEM_PE(bwa_inputs) |
            INDEX_BAM |
            MAPPING_SUMMARY

        // publish output:
        PUBLISH_BAM(INDEX_BAM.out, 'bwa.out.possorted', 'bam', 'bam', false)
        PUBLISH_BAM_INDEX(INDEX_BAM.out.map{it -> tuple(it[0], it[2])}, 'bwa.out.possorted.bam', 'bai', 'bam', false)
        PUBLISH_MAPPING_SUMMARY(MAPPING_SUMMARY.out, 'mapping_stats', 'tsv', 'bam', false)

    emit:
        INDEX_BAM.out

}

