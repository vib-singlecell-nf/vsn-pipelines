nextflow.preview.dsl=2

def extractSample(path) {
    (full, parentDir, id, whateverSuffix, compressionExtension) = (path =~ /(.+)\/(.+)_R1(.*)\.fastq(\.gz)?/)[0]
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
                path -> tuple(extractSample( "${path}" ), file("${path}"))
            }
    emit:
        channel
}
