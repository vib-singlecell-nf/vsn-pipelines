def _getBaseName(file) {
    (full, filename, process, ext) = (file.getName() =~ /(.+)\.SC(.+)\.(.+)/)[0]
    return filename
}

workflow getBaseName(file) {
    emit:
        _getBaseName(file)
}