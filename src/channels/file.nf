nextflow.preview.dsl=2

include {
    extractSample
} from '../utils/processes/files.nf'

workflow getChannel {

    take:
        glob
        sampleSuffixWithExtension // Suffix after the sample name in the file paths

    main:
        // Check whether multiple globs are provided
        if(glob.contains(',')) {
            glob = Arrays.asList(glob.split(',')); 
        }
        channel = Channel
            .fromPath(glob, checkIfExists: true)
            .map {
                path -> tuple(extractSample( "${path}", sampleSuffixWithExtension ), file("${path}"))
            }

    emit:
        channel

}

workflow getChannelWithIndex {

    take:
        glob
        sampleSuffixWithExtension // Suffix after the sample name in the file paths
        indexFileExtension // file extension of the paired index file (e.g. '.bai', '.tbi')

    main:
        // Check whether multiple globs are provided
        if(glob.contains(',')) {
            glob = Arrays.asList(glob.split(',')); 
        }
        channel = Channel
            .fromPath(glob, checkIfExists: true)
            .map {
                path -> tuple(extractSample( "${path}", sampleSuffixWithExtension ), file("${path}"), file("${path}${indexFileExtension}"))
            }

    emit:
        channel

}

workflow getChannelFromFilePath {

    take:
        filePath
        sampleSuffixWithExtension // Suffix after the sample name in the file paths
    
    main:
        channel = Channel.of(
            tuple(filePath)
        ).map {
            it -> tuple(extractSample( "${it[0]}", sampleSuffixWithExtension ), file("${it[0]}"))
        }

    emit:
        channel

}
