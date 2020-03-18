nextflow.preview.dsl=2

import java.nio.file.Paths

//////////////////////////////////////////////////////
//  process imports:

include SC__CELLRANGER__COUNT_WITH_LIBRARIES   from './../processes/count'    params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow CELLRANGER_COUNT_WITH_LIBRARIES {

    take:
        transcriptome
        featureRef
        libraries

    main:
        // Define the sampleId
        data = Channel.from(
            libraries
        ).splitCsv(
            header:true,
            limit: 1
        ).map {
            row -> tuple(
                row.sample,
                featureRef,
                libraries
            )
        }
        SC__CELLRANGER__COUNT_WITH_LIBRARIES( transcriptome, data )

    emit:
        SC__CELLRANGER__COUNT_WITH_LIBRARIES.out

}
