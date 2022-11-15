nextflow.enable.dsl=2

// process imports
include {
    SCTK__EXTRACT_HYDROP_ATAC_BARCODE as SCTK__EXTRACT_HYDROP_ATAC_BARCODE_2x384;
    SCTK__EXTRACT_HYDROP_ATAC_BARCODE as SCTK__EXTRACT_HYDROP_ATAC_BARCODE_3x96;
} from './../../src/singlecelltoolkit/processes/extract_hydrop_atac_barcode.nf'
include {
    TRIMGALORE__TRIM;
} from './../../src/trimgalore/processes/trim.nf'
include {
    FASTP__ADAPTER_TRIMMING as FASTP__TRIM;
} from './../../src/fastp/processes/adapter_trimming.nf'

include {
    SIMPLE_PUBLISH as PUBLISH_FASTQS_TRIMLOG_PE1;
    SIMPLE_PUBLISH as PUBLISH_FASTQS_TRIMLOG_PE2;
    SIMPLE_PUBLISH as PUBLISH_FASTQS_TRIMLOG_FASTP;
    SIMPLE_PUBLISH as PUBLISH_FRAGMENTS;
    SIMPLE_PUBLISH as PUBLISH_FRAGMENTS_INDEX;
} from '../../src/utils/processes/utils.nf'

// workflow imports:
include {
    BWA_MAPPING_PE;
} from './../../src/bwamaptools/main.nf'
include {
    PICARD__MERGE_SAM_FILES_AND_SORT;
} from './../../src/gatk/processes/merge_sam_files.nf'
include {
    SAMTOOLS__SORT_BAM;
} from './../../src/samtools/processes/sort_bam.nf'
include {
    BAM_TO_FRAGMENTS;
    DETECT_BARCODE_MULTIPLETS;
} from './../../src/barcard/main.nf'
//} from './../../src/sinto/main.nf'


include {
    barcode_correction as bc_correct_standard;
    barcode_correction as bc_correct_hydrop_2x384;
    barcode_correction as bc_correct_hydrop_3x96;
    biorad_bc as bc_correct_biorad;
} from './../../src/singlecelltoolkit/main.nf'


//////////////////////////////////////////////////////
//  Define the workflow

workflow ATAC_PREPROCESS_RAPID {

    take:
        metadata

    main:

        // import metadata
        data = Channel.from(metadata)
                      .splitCsv(
                          header:true,
                          sep: '\t'
                          )
                      .map {
                          row -> tuple(
                              row.sample_name + "___" + file(row.fastq_PE1_path)
                                  .getSimpleName()
                                  .replaceAll(row.sample_name,""),
                              row.technology,
                              file(row.fastq_PE1_path, checkIfExists: true),
                              row.fastq_barcode_path,
                              file(row.fastq_PE2_path, checkIfExists: true)
                              )
                      }
                      .branch {
                        biorad:       it[1] == 'biorad'
                        hydrop_3x96:  it[1] == 'hydrop_3x96'
                        hydrop_2x384: it[1] == 'hydrop_2x384'
                        standard:     true // capture all other technology types here
                      }

        /* standard data
           barcode correction */
        bc_correct_standard(data.standard)

        /* HyDrop ATAC
           extract barcode and correct */
           // HyDrop 3x96
        SCTK__EXTRACT_HYDROP_ATAC_BARCODE_3x96(data.hydrop_3x96, '3x96') \
            | bc_correct_hydrop_3x96
           // HyDrop 2x384
        SCTK__EXTRACT_HYDROP_ATAC_BARCODE_2x384(data.hydrop_2x384, '2x384') \
            | bc_correct_hydrop_2x384

        /* BioRad data
           extract barcode and correct */
        bc_correct_biorad(data.biorad)

        /* downstream steps */
        bc_correct_standard.out
            .mix(bc_correct_hydrop_3x96.out)
            .mix(bc_correct_hydrop_2x384.out)
            .mix(bc_correct_biorad.out) \
            | adapter_trimming \
            | mapping

    emit:
        // emit in a format compatible with getDataChannel output:
        bam = mapping.out.bam.map { it -> tuple(it[0], [it[1], it[2]], 'bam') }
        fragments = mapping.out.fragments.map { it -> tuple(it[0], [it[1], it[2]], 'fragments') }
}


/* sub-workflows used above */

workflow adapter_trimming {

    take:
        fastq_dex

    main:

        // run adapter trimming:
        switch(params.atac_preprocess_tools.adapter_trimming_method) {
            case 'Trim_Galore':
                fastq_dex_trim = TRIMGALORE__TRIM(fastq_dex);
                PUBLISH_FASTQS_TRIMLOG_PE1(fastq_dex_trim.map{ it -> tuple(it[0], it[3]) }, '.R1.trimming_report.txt', 'reports/trim');
                PUBLISH_FASTQS_TRIMLOG_PE2(fastq_dex_trim.map{ it -> tuple(it[0], it[4]) }, '.R2.trimming_report.txt', 'reports/trim');
                break;
            case 'fastp':
                fastq_dex_trim = FASTP__TRIM(fastq_dex);
                PUBLISH_FASTQS_TRIMLOG_FASTP(fastq_dex_trim.map{ it -> tuple(it[0], it[3]) }, '.fastp.trimming_report.html', 'reports/trim');
                break;
        }

    emit:
        fastq_dex_trim

}


workflow mapping {

    take:
        fastq_dex_trim

    main:

        // map with bwa mem:
        aligned_bam = BWA_MAPPING_PE(
            fastq_dex_trim.map { it -> tuple(it[0].split("___")[0], // [val(unique_sampleId),
                                             *it[0..2] ) // val(sampleId), path(fastq_PE1), path(fastq_PE2)]
                  })

        // split by sample size:
        aligned_bam.map{ it -> tuple(it[0].split("___")[0], it[1]) } // [ sampleId, bam ]
                   .groupTuple()
                   .branch {
                       to_merge: it[1].size() > 1
                       no_merge: it[1].size() == 1
                   }
                   .set { aligned_bam_size_split }

        // merge samples with multiple files:
        bam_merged = PICARD__MERGE_SAM_FILES_AND_SORT(aligned_bam_size_split.to_merge)

        // re-combine with single files:
        bam_merged.mix(aligned_bam_size_split.no_merge.map { it -> tuple(it[0], *it[1]) })
           .set { bam }
//           .set { aligned_bam_sample_merged }

        // sort the bam (not duplicate marked!)
        //sorted_bam = SAMTOOLS__SORT_BAM(merged_bam)

        // generate a fragments file:
        fragments = BAM_TO_FRAGMENTS(bam)

        DETECT_BARCODE_MULTIPLETS(fragments)

        // publish fragments output:
        PUBLISH_FRAGMENTS(fragments.map{ it -> tuple(it[0..1]) }, '.fragments.tsv.gz', 'fragments')
        PUBLISH_FRAGMENTS_INDEX(fragments.map{ it -> tuple(it[0],it[2]) }, '.fragments.tsv.gz.tbi', 'fragments')

    emit:
        bam
        fragments

}
