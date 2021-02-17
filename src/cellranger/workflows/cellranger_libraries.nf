nextflow.enable.dsl=2

import java.nio.file.Paths

//////////////////////////////////////////////////////
//  process imports:

include {
    MKFASTQ
} from './mkfastq' params(params)
include {
    SC__CELLRANGER__PREFLIGHT;
} from './../processes/preflight' params(params)
include {
    SC__CELLRANGER__COUNT_WITH_LIBRARIES
} from './../processes/count' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow CELLRANGER_LIBRARIES {

    take:
        mkfastq_csv
        runFolder
        transcriptome
        featureRef

    main:
        // Sanity Checking
        libMap = params.tools.cellranger.librariesMap
        if (! (libMap instanceof Map)) {
            throw new Exception("When running the full cellranger pipeline with libraries, you must specify the librariesMap (see docs).")
        }

        librariesFiles = params.tools.cellranger.count.libraries

        if (!(librariesFiles instanceof Map) && librariesFiles) {
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

        // Get any existing libraries files and process them
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
        .set { oldRunsData }
        
        // Run MKFASTQ on current run
        SC__CELLRANGER__PREFLIGHT()

        data = MKFASTQ(mkfastq_csv, runFolder)

        // Get Library info for MKFASTQ run from params
        Set libSamples = []
        samplesToPool = [:]
        Channel.fromList(
            libMap.collectMany { samplePoolName, samplePool ->          
                samplePool.collect { 
                    libSamples.add(it.sampleName)
                    samplesToPool[it.sampleName] = samplePoolName
                    return tuple(it.sampleName, it.assay)
                }
            }
        )
        .concat(data)
        .groupTuple(size: 2)
        .map { sample, info ->
            return tuple(samplesToPool[sample], tuple(file(info[1]), sample, info[0]))
        }
        .groupTuple()
        .set { curRunData }

        // Make sure that all of the samples listed in the libraries file in the config were outputs of MKFASTQ
        curRunData.collect{ pool -> pool[1].collect { it[1] } }
        .map {
            if (! it.containsAll(libSamples)) {
                diff = libSamples.plus(it)
                diff.removeAll(it)
                throw new Exception("Could not find sample(s) ${diff} from librariesMap in output of MKFASTQ. Ensure your librariesMap and samplesheet match.")
            }
        }

        // Finally join the old and new data ready for counting
        curRunData.concat(oldRunsData)
        .groupTuple()
        .map { pool -> 
            transposed =  GroovyCollections.transpose(pool[1].collectMany { it })
            return tuple(
                pool[0], 
                transposed[0],
                transposed[1],
                transposed[2]              
            )
        }
        .set { data }

        SC__CELLRANGER__COUNT_WITH_LIBRARIES( transcriptome, featureRef, data )
    
    emit:
        SC__CELLRANGER__COUNT_WITH_LIBRARIES.out

}