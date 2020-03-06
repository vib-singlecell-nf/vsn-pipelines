nextflow.preview.dsl=2

include getMEXChannel as getTenXCellRangerMEXChannel from './tenx' params(params)
include getH5Channel as getTenXCellRangerH5Channel from './tenx' params(params)
include getChannel as getFileChannel from './file' params(params)

workflow getDataChannel {

    main:
        data = Channel.empty()
        if(params.data.containsKey("tenx") && params.data.tenx.containsKey("cellranger_mex")) {
            data = data.concat(
                getTenXCellRangerMEXChannel(
                    params.data.tenx.cellranger_mex
                ).map {
                    it -> tuple(it[0], it[1], "10x_cellranger_mex", "h5ad")
                }
            ).view()
        }
        if(params.data.containsKey("tenx_atac") && params.data.tenx_atac.containsKey("cellranger_mex")) {
            data = data.concat(
                getTenXCellRangerMEXChannel(
                    params.data.tenx_atac.cellranger_mex
                ).map {
                    it -> tuple(it[0], it[1], "10x_atac_cellranger_mex", "cistopic_rds")
                }
            ).view()
        }
        if(params.data.containsKey("tenx") && params.data.tenx.containsKey("cellranger_h5")) {
            data = data.concat(
                getTenXCellRangerH5Channel( 
                    params.data.tenx.cellranger_h5
                ).map {
                    it -> tuple(it[0], it[1], "10x_cellranger_h5", "h5ad")
                }
            ).view()
        }
        if(params.data.containsKey("h5ad")) {
            data = data.concat(
                getFileChannel( 
                    params.data.h5ad.file_paths,
                    params.data.h5ad.suffix
                ).map {
                    it -> tuple(it[0], it[1], "h5ad", "h5ad")
                }
            ).view()
        }
        if(params.data.containsKey("tsv")) {
            data = data.concat(
                getFileChannel( 
                    params.data.tsv.file_paths,
                    params.data.tsv.suffix
                ).map {
                    it -> tuple(it[0], it[1], "tsv", "h5ad")
                }
            ).view()
        }
        if(params.data.containsKey("csv")) {
            data = data.concat(
                getFileChannel( 
                    params.data.csv.file_paths,
                    params.data.csv.suffix
                ).map {
                    it -> tuple(it[0], it[1], "csv", "h5ad")
                }
            ).view()
        }
        if(params.data.containsKey("seurat_rds")) {
            data = data.concat(
                getFileChannel( 
                    params.data.seurat_rds.file_paths,
                    params.data.seurat_rds.suffix
                ).map {
                    it -> tuple(it[0], it[1], "seurat_rds", "h5ad")
                }
            ).view()
        }
        data.ifEmpty { exit 1, "Pipeline cannot run: no data provided." }

    emit:
        data

}
