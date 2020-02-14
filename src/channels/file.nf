nextflow.preview.dsl=2

include '../utils/processes/files.nf'

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
