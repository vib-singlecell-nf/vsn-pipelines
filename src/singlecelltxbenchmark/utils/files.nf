def getBaseName(file) {
    (full, filename, process, ext) = ( file.getName() =~ /(.+)\.SC(.+)\.(.+)/)[0]
    return filename
}