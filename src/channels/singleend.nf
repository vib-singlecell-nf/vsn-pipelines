nextflow.preview.dsl=2

def extractSample(path) {
    pattern = /(.+)\/(.+)_R[1-2](.*)\.fastq(\.gz)?/
    (full, parentDir, id, whateverSuffix, compressionExtension) = (path =~ pattern)[0]

    return id
}

workflow getChannel {
    get:
        glob
    main:
        // Check whether multiple globs are provided
        if(glob.contains(',')) {
            glob = Arrays.asList(glob.split(',')); 
        }
        channel = Channel
            .fromPath(glob, checkIfExists: true)
            .map {
                path -> tuple(extractSample( "${path}" ), path("${path}"))
            }
    emit:
        channel
}
