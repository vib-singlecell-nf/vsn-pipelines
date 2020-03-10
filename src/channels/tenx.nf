nextflow.preview.dsl=2

CELLRANGER_H5_REGEX = /(.+)\/(.+)\/outs\/(.+)\.h5/
CELLRANGER_MEX_REGEX = /(.+)\/(.+)\/outs/

def extractSampleFromH5(path) {
    if (!(path ==~ CELLRANGER_H5_REGEX))
        throw new Exception("Incorrect Cell Ranger .h5 path. The parameter params.data.tenx.cellranger_h5 in the config file should point to the .h5 file.")
    // Allow to detect data generated by CellRanger prior and post to version 3.
    (full, parentDir, id, filename) = (path =~ CELLRANGER_H5_REGEX)[0]
    return id
}

workflow getH5Channel {

    take:
        glob

    main:
        // Check whether multiple globs are provided
        if(glob.contains(',')) {
            glob = Arrays.asList(glob.split(',')); 
        }
        channel = Channel
            .fromPath(glob, type: 'file', checkIfExists: true)
            .map {
                filePath -> tuple(extractSampleFromH5( "${filePath}" ), file("${filePath}"))
            }

    emit:
        channel

}

def extractSampleFromMEX(path) {
    // Allow to detect data generated by CellRanger prior and post to version 3.
    if (!(path ==~ CELLRANGER_MEX_REGEX))
        throw new Exception("Incorrect Cell Ranger MEX path. The parameter params.data.tenx.cellranger_mex in the config file should point to the outs folder.")
    (full, parentDir, id) = (path =~ CELLRANGER_MEX_REGEX)[0]
    return id
}

workflow getMEXChannel {

    take:
        glob

    main:
        // Check whether multiple globs are provided
        if(glob.contains(',')) {
            glob = Arrays.asList(glob.split(',')); 
        }
        channel = Channel
            .fromPath(glob, type: 'dir', checkIfExists: true)
            .map {
                filePath -> tuple(extractSampleFromMEX( "${filePath}" ), file("${filePath}"))
            }

    emit:
        channel

}
