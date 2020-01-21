nextflow.preview.dsl=2

def extractSample(path, suffix) {
    if(!path.endsWith(".h5ad"))
        throw new Exception("Wrong channel used for data: "+ path)
    // Extract the sample name based on the given path and on the given suffix
    suffix = suffix.replace(".","\\.")
    pattern = /(.+)\/(.+)${suffix}/
    (full, parentDir, id) = (path =~ pattern)[0]
    return id
}

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
