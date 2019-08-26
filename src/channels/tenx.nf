nextflow.preview.dsl=2

def extractSample(path) {
    (full, parentDir, id) = (path =~ /(.+)\/(.+)\/outs\/filtered_feature_bc_matrix$/)[0]
    return id
}

workflow getChannel(glob) {
    emit:
        Channel
            .fromPath(glob, type: 'dir')
            .map { 
                path -> tuple(extractSample( "${path}" ), file("${path}"))
            }
}