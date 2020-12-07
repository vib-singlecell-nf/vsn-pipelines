nextflow.preview.dsl=2

import java.nio.file.Paths

//////////////////////////////////////////////////////
//  process imports:

include {
    SC__CELLRANGER_ATAC__COUNT_WITH_METADATA
} from './../processes/count' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow CELLRANGER_ATAC_COUNT_WITH_METADATA {

    take:
        metadata

    main:
        // Define the sampleId
        data = Channel.from(
            metadata
        ).splitCsv(
            header:true,
            sep: '\t'
        ).map {
            row -> tuple(
                row.processing_date + "__" + row.sample_name + "__" + row.short_uuid,
                row.fastqs_sample_prefix,
                file(Paths.get(row.fastqs_parent_dir_path, row.fastqs_dir_name)),
                // Begin CellRanger parameters
                row.expect_cells
            )
        }
        SC__CELLRANGER_ATAC__COUNT_WITH_METADATA( data )

    emit:
        SC__CELLRANGER_ATAC__COUNT_WITH_METADATA.out

}

