nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:
include SC__SCRUBLET__DOUBLET_DETECTION from "../processes/doublet_detection" params(params)

//////////////////////////////////////////////////////
//  Import from external modules:
include ANNOTATE_BY_CELL_METADATA from '../../utils/workflows/annotateByCellMetadata.nf' params(params)
include FILTER_BY_CELL_METADATA from '../../utils/workflows/filterByCellMetadata' params(params)

workflow DOUBLET_REMOVAL {

    take:
        // Expects (sampleId, adataRaw, adataWithHvgInfo, stashedParams | null, nPrinComps | null)
        data

    main:
        SC__SCRUBLET__DOUBLET_DETECTION(
            data
        )

        ANNOTATE_BY_CELL_METADATA(
            data.map {
                it -> tuple(it[0], it[1])
            },
            SC__SCRUBLET__DOUBLET_DETECTION.out.map {
                it -> tuple(it[0], it[1])
            },
            "scrublet"
        )

        FILTER_BY_CELL_METADATA(
            ANNOTATE_BY_CELL_METADATA.out,
            "scrublet"
        )
        // SC__SCRUBLET__REPORT( 
        // )

    // emit:
    //     doublet_removed = FILTER_BY_CELL_METADATA.out
    //     report = SC__SCRUBLET__REPORT.out.ipynb

}