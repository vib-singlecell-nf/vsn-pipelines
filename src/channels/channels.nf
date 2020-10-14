nextflow.preview.dsl=2

include {
    getOutsChannel as getTenXCellRangerOutsChannel;
    getH5Channel as getTenXCellRangerH5Channel;
    getMEXChannel as getTenXCellRangerMEXChannel;
} from './tenx' params(params)
include getChannel as getFileChannel from './file' params(params)

boolean isCollectionOrArray(object) {    
    [Collection, Object[]].any { it.isAssignableFrom(object.getClass()) }
}

def isOuts = { glob ->
    // Avoid file() which will resolve the glob
    def _isOuts = { arr -> arr.collect { new File(it).getName() == "outs" }.sum(0, { it ? 1 : 0 }) == glob.size() }
    if(glob.contains(','))
        return _isOuts(Arrays.asList(glob.split(',')))
    if(isCollectionOrArray(glob))
        return _isOuts(glob)
    return new File(glob).getName() == "outs"
}

workflow getDataChannel {

    main:
        data = Channel.empty()
        // Define the output file format
        outputFileFormat = "h5ad"
        // This allows to use it dynamically (i.e.: addParams)
        if(params.containsKey("off")) {
            outputFileFormat = params.off
        } else {
            // If not dynamically set, we use h5ad by default
            outputFileFormat = "h5ad"
            if(params.sc.file_converter.containsKey("off")) {
                outputFileFormat = params.sc.file_converter.off
            }
        }

        if(params.data.containsKey("tenx") && params.data.tenx.containsKey("cellranger_mex")) {
            if(isOuts(params.data.tenx.cellranger_mex)) {
                data = data.concat(
                    getTenXCellRangerOutsChannel(
                        params.data.tenx.cellranger_mex
                    ).map {
                        it -> tuple(it[0], it[1], "10x_cellranger_mex_outs", outputFileFormat)
                    }
                ).view()
            } else {
                data = data.concat(
                    getTenXCellRangerMEXChannel(
                        params.data.tenx.cellranger_mex
                    ).map {
                        it -> tuple(it[0], it[1], "10x_cellranger_mex", outputFileFormat)
                    }
                ).view()
            }
        }
        if(params.data.containsKey("tenx_atac") && params.data.tenx_atac.containsKey("cellranger_mex")) {
            if(isOuts(params.data.tenx_atac.cellranger_mex)) {
                data = data.concat(
                    getTenXCellRangerOutsChannel(
                        params.data.tenx_atac.cellranger_mex
                    ).map {
                        it -> tuple(it[0], it[1], "10x_atac_cellranger_mex_outs", outputFileFormat)
                    }
                ).view()
            } else {
                data = data.concat(
                    getTenXCellRangerMEXChannel(
                        params.data.tenx_atac.cellranger_mex
                    ).map {
                        it -> tuple(it[0], it[1], "10x_atac_cellranger_mex", outputFileFormat)
                    }
                ).view()
            }
        }
        if(params.data.containsKey("tenx") && params.data.tenx.containsKey("cellranger_h5")) {
            if(isOuts(params.data.tenx.cellranger_h5)) {
                data = data.concat(
                    getTenXCellRangerOutsChannel(
                        params.data.tenx.cellranger_h5
                    ).map {
                        it -> tuple(it[0], it[1], "10x_cellranger_h5_outs", outputFileFormat)
                    }
                ).view()
            } else {
                data = data.concat(
                    getTenXCellRangerH5Channel( 
                        params.data.tenx.cellranger_h5
                    ).map {
                        it -> tuple(it[0], it[1], "10x_cellranger_h5", outputFileFormat)
                    }
                ).view()
            }
        }
        if(params.data.containsKey("h5ad")) {
            dataH5ad = params.data.h5ad
            def filePaths = null
            def suffix = null
            if(!dataH5ad.containsKey("file_paths") && !dataH5ad.containsKey("suffix")) {
                filePaths = dataH5ad.collect { k,v -> v["file_paths"] }.flatten()
                suffix = dataH5ad.collect { k,v -> v["suffix"] }.flatten()
            } else {
                filePaths = dataH5ad.file_paths
                suffix = dataH5ad.suffix
            }
            data = data.concat(
                getFileChannel(
                    filePaths,
                    suffix
                ).map {
                    it -> tuple(it[0], it[1], "h5ad", outputFileFormat)
                }
            ).view()
        }
        if(params.data.containsKey("loom")) {
            data = data.concat(
                getFileChannel( 
                    params.data.loom.file_paths,
                    params.data.loom.suffix
                ).map {
                    it -> tuple(it[0], it[1], "loom", outputFileFormat)
                }
            ).view()
        }
        if(params.data.containsKey("tsv")) {
            data = data.concat(
                getFileChannel( 
                    params.data.tsv.file_paths,
                    params.data.tsv.suffix
                ).map {
                    it -> tuple(it[0], it[1], "tsv", outputFileFormat)
                }
            ).view()
        }
        if(params.data.containsKey("csv")) {
            data = data.concat(
                getFileChannel( 
                    params.data.csv.file_paths,
                    params.data.csv.suffix
                ).map {
                    it -> tuple(it[0], it[1], "csv", outputFileFormat)
                }
            ).view()
        }
        if(params.data.containsKey("seurat_rds")) {
            data = data.concat(
                getFileChannel( 
                    params.data.seurat_rds.file_paths,
                    params.data.seurat_rds.suffix
                ).map {
                    it -> tuple(it[0], it[1], "seurat_rds", outputFileFormat)
                }
            ).view()
        }
        data.ifEmpty { exit 1, "Pipeline cannot run: no data provided." }

    emit:
        data

}
