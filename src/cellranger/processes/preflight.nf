nextflow.enable.dsl=2

toolParams = params.tools.cellranger

process SC__CELLRANGER__PREFLIGHT {

    exec:
        if (toolParams.containsKey('container')) {
            if (toolParams.container == '/path/to/cellranger/cellranger' || toolParams.container == '') {
                throw new Exception("You must specify a container image for Cellranger!")
            }   
        } else {
            throw new Exception("No container entry")
        }
        
}
