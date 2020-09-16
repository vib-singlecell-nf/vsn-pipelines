def getBaseName(file) {
    (full, filename, process, ext) = ( file.getName() =~ /(.+)\.SC(.+)\.(.+)/)[0]
    return filename
}

def extractSample(path, suffix) {
    // Extract the sample name based on the given path and on the given suffix
    if(suffix instanceof String)
        suffix = [suffix]
    suffix = suffix.collect { it.replace(".","\\.") }
    for (String sufx : suffix) {
        pattern = /(.+)\/(.+)${sufx}/
        res = (path =~ pattern)
        if(res.size() == 0) continue
        if(res.size() == 1) {
            (full, parentDir, id) = res[0]
            return id
        }
    }
    throw new Exception("VSN ERROR: the suffix params couldn't match any of the file paths. Make sure the suffix param exist in the file paths.")
}
