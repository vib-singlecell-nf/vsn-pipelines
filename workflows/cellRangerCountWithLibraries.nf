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
        if(libraries.contains(',')) {
            libraries = Arrays.asList(libraries.split(',')); 
        }
        Channel
            .fromPath(libraries, checkIfExists: true)
            .map {
                path -> file("${path}")
            } 
            .multiMap {
            csv: it.splitCsv(header:true, limit: 1)
            files: file(it)
            }
            .set { data }
        
        data = data.csv.merge( data.files )
            .map {
            tuple(
                it[0].sample,
                it[1]
            )

        } 

        SC__CELLRANGER__COUNT_WITH_LIBRARIES( transcriptome, featureRef, data )

}
