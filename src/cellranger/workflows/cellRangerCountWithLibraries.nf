nextflow.enable.dsl=2

import java.nio.file.Paths

//////////////////////////////////////////////////////
//  process imports:

include {
    SC__CELLRANGER__COUNT_WITH_LIBRARIES;
} from './../processes/count' params(params)
include {
    SC__CELLRANGER__PREFLIGHT;
} from './../processes/preflight' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow CELLRANGER_COUNT_WITH_LIBRARIES {

    take:
        transcriptome
        featureRef
        librariesFiles

    main:
        if (!(librariesFiles instanceof Map) && libraryFiles) {
            poolName = params.global.containsKey('project_name') ? params.global.project_name : ''
            if(librariesFiles.contains(',')) {
                librariesFiles = Arrays.asList(librariesFiles.split(','))
                .collectEntries { f -> 
                    [(file(f).baseName): f]
                }
            } else {
                librariesFiles = [(file(librariesFiles).baseName): librariesFiles]
            }
        }

        Channel.from(
            librariesFiles.collect { samplePoolName, libraryFiles -> 
                if(libraryFiles.contains(',')) {
                    libraryFiles = Arrays.asList(libraryFiles.split(',')).collect { f -> file(f) } 
                } else {
                    libraryFiles = [file(libraryFiles)]
                }
                libFData = libraryFiles.collectMany { libraryFile ->
                    libraryFile
                    .splitCsv(skip: 1)
                    .collect {
                        return tuple(file(it[0]), it[1], it[2])
                    }
                }
                return tuple(samplePoolName, libFData)
            }
        )
        .map { pool -> 
            transposed =  GroovyCollections.transpose(pool[1])
            return tuple(
                pool[0], 
                transposed[0],
                transposed[1],
                transposed[2]              
            )
        }
        .set { data }
        
        SC__CELLRANGER__PREFLIGHT()
        SC__CELLRANGER__COUNT_WITH_LIBRARIES( transcriptome, featureRef, data )

    emit:
        SC__CELLRANGER__COUNT_WITH_LIBRARIES.out
}


