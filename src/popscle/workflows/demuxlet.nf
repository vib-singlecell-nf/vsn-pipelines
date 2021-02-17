nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include {
    DSC_PILEUP_FILTERED
} from './dsc_pileup.nf' params(params)
include {
    SC__POPSCLE__FREEMUXLET;
    SC__POPSCLE__DEMUXLET;
} from '../processes/demuxlet.nf' params(params)

//////////////////////////////////////////////////////
// Define the workflow

workflow get_bam_barcodes_from_cellranger_rna {
    take:
        data // cellranger-rna output channel [ sampleId, path ]
    main:
        out = data.map{ it -> tuple(it[0],
                                    file(it[1]+"/possorted*bam.bam"),
                                    file(it[1]+"/filtered_*_bc_matrix/barcodes.tsv*")
                                    ) }
    emit:
        out
}

workflow get_bam_barcodes_from_cellranger_atac {
    take:
        data // cellranger-atac output channel [ sampleId, path ]
    main:
        out = data.map{ it -> tuple(it[0],
                                    file(it[1]+"/possorted*bam.bam"),
                                    file(it[1]+"/filtered_peak_bc_matrix/barcodes.tsv*")
                                    ) }
    emit:
        out
}

workflow get_bam_barcodes_from_cellranger_arc {
    take:
        data // cellranger-arc output channel [ sampleId, path ] // not yet implemented
    main:
        out = data.map{ it -> tuple(it[0],
                                    file(it[1]+"/atac_possorted*bam.bam"),
                                    file(it[1]+"/filtered_feature_bc_matrix/barcodes.tsv*")
                                    ) }
    emit:
        out
}


workflow data_channel_to_bam_barcodes {

    take:
        data // standard data channel [ sampleId, path, type, format, group ]

    main:

        data.branch {
              bam:     it[2] == 'bam'
              tsv:     it[2] == 'tsv'
              cr_rna:  it[2] == '10x_cellranger_mex_outs'
              cr_atac: it[2] == '10x_atac_cellranger_mex_outs'
              cr_arc:  it[2] == '10x_arc_cellranger_mex_outs'
            }
            .set{datab}


        // cellranger rna:
        data_cr_rna = get_bam_barcodes_from_cellranger_rna(datab.cr_rna)

        // cellranger atac:
        data_cr_atac = get_bam_barcodes_from_cellranger_atac(datab.cr_atac)

        // cellranger rc:
        data_cr_arc = get_bam_barcodes_from_cellranger_arc(datab.cr_arc)

        // combine the bam+barcodes channels:
        datab.bam.join(datab.tsv)
                 .map{ it -> tuple(it[0],    // sampleId
                                   it[1][0], // bam
                                   it[5]     // tsv (barcodes)
                                   ) }
                 // combine outputs:
                 .mix(data_cr_rna)
                 .mix(data_cr_atac)
                 .mix(data_cr_arc)
                 .set{out}

    emit:
        out
}


workflow FREEMUXLET {

    take:
        data

    main:
        DSC_PILEUP_FILTERED(data) |
        SC__POPSCLE__FREEMUXLET
    
    emit:
        SC__POPSCLE__FREEMUXLET.out
}

workflow DEMUXLET {

    take:
        data

    main:
        vcf = file(params.tools.popscle.vcf)
        DSC_PILEUP_FILTERED(data)
        SC__POPSCLE__DEMUXLET(DSC_PILEUP_FILTERED.out, vcf)
    
    emit:
        SC__POPSCLE__DEMUXLET.out
}

