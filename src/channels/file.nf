nextflow.enable.dsl=2

include {
    extractSample
} from '../utils/processes/files.nf'

workflow getChannel {

    take:
        glob
        sampleSuffixWithExtension // Suffix after the sample name in the file paths
        groups

    main:
        // Check whether multiple globs are provided
        if(glob.contains(',')) {
            glob = Arrays.asList(glob.split(',')); 
        }
        data_channel = Channel
            .fromPath(glob, checkIfExists: true)
            .map {
                path -> tuple(
                    *extractSample(
                        "${path}",
                        sampleSuffixWithExtension,
                        groups
                    ),
                    file("${path}")
                )
            }.map {
                // reorder: sample ID, file path, tag
                it -> tuple(it[0], it[2], it[1])
            }

    emit:
        data_channel

}

workflow getChannelWithIndex {

    take:
        glob
        sampleSuffixWithExtension // Suffix after the sample name in the file paths
        indexFileExtension // file extension of the paired index file (e.g. '.bai', '.tbi')
        groups

    main:
        // Check whether multiple globs are provided
        if(glob.contains(',')) {
            glob = Arrays.asList(glob.split(',')); 
        }
        data_channel = Channel
            .fromPath(glob, checkIfExists: true)
            .map {
                path -> tuple(*extractSample("${path}", sampleSuffixWithExtension, groups), file("${path}"), file("${path}${indexFileExtension}"))
            }
            .map {
                // reorder: sample ID, [file path, file index path], tag
                it -> tuple(it[0], [it[2],it[3]], it[1])
            }

    emit:
        data_channel

}

workflow getChannelFromFilePath {

    take:
        filePath
        sampleSuffixWithExtension // Suffix after the sample name in the file paths
        groups
    
    main:
        data_channel = Channel.of(
            tuple(filePath)
            )
            .map {
                it -> tuple(*extractSample("${it[0]}", sampleSuffixWithExtension, groups), file("${it[0]}"))
            }
            .map {
                // reorder: sample ID, file path, tag
                it -> tuple(it[0], it[2], it[1])
            }

    emit:
        data_channel

}

