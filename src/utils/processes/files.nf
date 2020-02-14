def getBaseName(file) {
    (full, filename, process, ext) = ( file.getName() =~ /(.+)\.SC(.+)\.(.+)/)[0]
    return filename
}

def extractSample(path, suffix) {
    if(!path.endsWith(".h5ad"))
        throw new Exception("Wrong channel used for data: "+ path)
    // Extract the sample name based on the given path and on the given suffix
    suffix = suffix.replace(".","\\.")
    pattern = /(.+)\/(.+)${suffix}/
    (full, parentDir, id) = (path =~ pattern)[0]
    return id
}
