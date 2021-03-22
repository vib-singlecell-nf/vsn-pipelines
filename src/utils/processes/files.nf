
def getBaseName(file, suffix) {
    // Default value suffix = "SC" does not work! Weird...
    res = (file.getName() =~ /(.+)\.${suffix}(.+)\.(.+)/)
    if(res.size() == 0) {
        throw new Exception("VSN ERROR: Cannot get base name.")
    }
    (full, filename, process, ext) = res[0]
    return filename
}

def extractSample(path, suffix, groups) {
    // Extract the sample name based on the given path and on the given suffix
    def _suffix = suffix instanceof String ? [suffix] : suffix
    _suffix = _suffix.collect { it.replace(".","\\.") }
    for(int i = 0; i<_suffix.size(); i++) {
        def sufx = _suffix[i]
        
        def pattern = /(.+)\/(.+)${sufx}/
        def res = (path =~ pattern)
        if(res.size() == 0) continue
        if(res.size() == 1) {
            def (full, parentDir, id) = res[0]
            if(groups != null) {
                return new Tuple(id, groups[i])
            } else {
                return new Tuple(id, 'NULL')
            }
        }
    }
    throw new Exception("VSN ERROR: the suffix params couldn't match any of the file paths. Make sure the suffix param exist in the file paths.")
}
