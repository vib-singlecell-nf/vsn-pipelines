nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the utils module:
include {
    ANNOTATE_BY_CELL_METADATA as ANNOTATE_BY_CELL_METADATA;
} from '../../utils/workflows/annotateByCellMetadata.nf' params(params)
include {
    FILTER_BY_CELL_METADATA
} from '../../utils/workflows/filterByCellMetadata' params(params)
include {
    UTILS__REPORT_TO_HTML
} from '../../utils/processes/reports.nf' params(params)

include {
    PUBLISH as PUBLISH_SCRUBLET_OBJECT;
    PUBLISH as PUBLISH_H5AD_DOUBLETS_ANNOTATED;
    PUBLISH as PUBLISH_H5AD_DOUBLETS_REMOVED;
 } from '../../utils/workflows/utils.nf' params(params)

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    SC__SCRUBLET__DOUBLET_DETECTION
} from "../processes/doublet_detection" params(params)
include {
    SC__SCRUBLET__DOUBLET_DETECTION_REPORT
} from "../processes/reports" params(params)


workflow DOUBLET_REMOVAL {

    take:
        // Expects (sampleId, adataRaw, adataWithHvgInfo, stashedParams || null, nPrinComps || null)
        data
        // Expects (sampleId, adataWithEmbeddings, stashedParams)
        finalProcessedData

    main:
        SC__SCRUBLET__DOUBLET_DETECTION(
            data
        )
        // Publish *.ScrubletObject.pklz
        PUBLISH_SCRUBLET_OBJECT(
            SC__SCRUBLET__DOUBLET_DETECTION.out.map { 
                it -> tuple(it[0], it[2])
            },
            null,
            "h5ad",
            "scrublet",
            false
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

        // Publish
        PUBLISH_H5AD_DOUBLETS_ANNOTATED(
            ANNOTATE_BY_CELL_METADATA.out,
            null,
            "h5ad",
            "scrublet",
            false
        )

        FILTER_BY_CELL_METADATA(
            ANNOTATE_BY_CELL_METADATA.out,
            "scrublet"
        )

        // Publish
        PUBLISH_H5AD_DOUBLETS_REMOVED(
            FILTER_BY_CELL_METADATA.out,
            null,
            "h5ad",
            "scrublet",
            false
        )
        

        SC__SCRUBLET__DOUBLET_DETECTION_REPORT(
            file(workflow.projectDir + params.tools.scrublet.doublet_detection.report_ipynb),
                SC__SCRUBLET__DOUBLET_DETECTION.out.map {
                // Extract the Scrublet object file
                it -> tuple(it[0], it[2])
            }.join(
                // Get the h5ad with Scrublet info
                ANNOTATE_BY_CELL_METADATA.out
            ).join(
                finalProcessedData.map {
                    // Extract the Scrublet object file
                    it -> tuple(it[0], it[1])
                }
            )
        )
        UTILS__REPORT_TO_HTML(
            SC__SCRUBLET__DOUBLET_DETECTION_REPORT.out
        )

    emit:
        doublet_detection = SC__SCRUBLET__DOUBLET_DETECTION.out
        data_annotated = ANNOTATE_BY_CELL_METADATA.out
        data_doublets_removed = FILTER_BY_CELL_METADATA.out
        report = SC__SCRUBLET__DOUBLET_DETECTION_REPORT.out

}