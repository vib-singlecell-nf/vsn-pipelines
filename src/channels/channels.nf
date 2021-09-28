nextflow.enable.dsl=2

include {
    getOutsChannel as getTenXCellRangerOutsChannel;
    getH5Channel as getTenXCellRangerH5Channel;
    getMEXChannel as getTenXCellRangerMEXChannel;
} from './tenx' params(params)
include {
    getChannel as getFileChannel;
    getChannelWithIndex as getFileChannelWithIndex;
} from './file' params(params)
include {
    isCollectionOrArray;
} from './../utils/processes/utils'


def isOuts = { glob ->
    // Avoid file() which will resolve the glob
    def _isOuts = { arr -> arr.collect { new File(it).getName() == "outs" }.sum(0, { it ? 1 : 0 }) == arr.size() }
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
            if(params.utils.file_converter.containsKey("off")) {
                outputFileFormat = params.utils.file_converter.off
            }
        }

        if(params.data.containsKey("tenx") && params.data.tenx.containsKey("cellranger_mex")) {
            if(isOuts(params.data.tenx.cellranger_mex)) {
                data = data.concat(
                    getTenXCellRangerOutsChannel(
                        params.data.tenx.cellranger_mex
                    ).map {
                        it -> tuple(it[0], it[1], "10x_cellranger_mex_outs", outputFileFormat, 'NULL')
                    }
                )
            } else {
                data = data.concat(
                    getTenXCellRangerMEXChannel(
                        params.data.tenx.cellranger_mex
                    ).map {
                        it -> tuple(it[0], it[1], "10x_cellranger_mex", outputFileFormat, 'NULL')
                    }
                )
            }
        }
        if(params.data.containsKey("tenx_atac") && params.data.tenx_atac.containsKey("cellranger_mex")) {
            if(isOuts(params.data.tenx_atac.cellranger_mex)) {
                data = data.concat(
                    getTenXCellRangerOutsChannel(
                        params.data.tenx_atac.cellranger_mex
                    ).map {
                        it -> tuple(it[0], it[1], "10x_atac_cellranger_mex_outs", outputFileFormat, 'NULL')
                    }
                )
            } else {
                data = data.concat(
                    getTenXCellRangerMEXChannel(
                        params.data.tenx_atac.cellranger_mex
                    ).map {
                        it -> tuple(it[0], it[1], "10x_atac_cellranger_mex", outputFileFormat, 'NULL')
                    }
                )
            }
        }
        if(params.data.containsKey("tenx_arc") && params.data.tenx_arc.containsKey("cellranger_mex")) {
            if(isOuts(params.data.tenx_arc.cellranger_mex)) {
                data = data.concat(
                    getTenXCellRangerOutsChannel(
                        params.data.tenx_arc.cellranger_mex
                    ).map {
                        it -> tuple(it[0], it[1], "10x_arc_cellranger_mex_outs", outputFileFormat, 'NULL')
                    }
                )
            } else {
                data = data.concat(
                    getTenXCellRangerMEXChannel(
                        params.data.tenx_arc.cellranger_mex
                    ).map {
                        it -> tuple(it[0], it[1], "10x_arc_cellranger_mex", outputFileFormat, 'NULL')
                    }
                )
            }
        }
        if(params.data.containsKey("tenx") && params.data.tenx.containsKey("cellranger_h5")) {
            if(isOuts(params.data.tenx.cellranger_h5)) {
                data = data.concat(
                    getTenXCellRangerOutsChannel(
                        params.data.tenx.cellranger_h5
                    ).map {
                        it -> tuple(it[0], it[1], "10x_cellranger_h5_outs", outputFileFormat, 'NULL')
                    }
                )
            } else {
                data = data.concat(
                    getTenXCellRangerH5Channel( 
                        params.data.tenx.cellranger_h5
                    ).map {
                        it -> tuple(it[0], it[1], "10x_cellranger_h5", outputFileFormat, 'NULL')
                    }
                )
            }
        }
        if(params.data.containsKey("h5ad")) {
            dataH5ad = params.data.h5ad
            def filePaths = null
            def suffix = null
            def groups = null
            if(!dataH5ad.containsKey("file_paths") && !dataH5ad.containsKey("suffix")) {
                filePaths = dataH5ad.collect { k,v -> v["file_paths"] }.flatten()
                suffix = dataH5ad.collect { k,v -> v["suffix"] }.flatten()
                numGroups = dataH5ad.size()
                numUniqueSuffixes = suffix.unique { a, b -> a <=> b }.size()
                if(numGroups != numUniqueSuffixes) {
                    throw new Exception("VSN ERROR: Number of data groups ("+ numGroups +") should be equal to the number of unique suffixes ("+ numUniqueSuffixes +"). Fix the group definitions in params.data.")
                }
                // Check if any suffix overlap
                for(int i=0; i<suffix.size(); i++) {
                    for(int j=0; j<suffix.size(); j++) {
                        if(i==j) continue
                        if(suffix[i].length() > suffix[j].length()) {
                            if(suffix[i].contains(suffix[j])) {
                                throw new Exception("VSN ERROR: some suffixes are overlapping ("+ suffix[j] +" with "+ suffix[i] +")")
                            }
                        } else {
                            if(suffix[j].contains(suffix[i])) {
                                throw new Exception("VSN ERROR: some suffixes are overlapping ("+ suffix[i] +" with "+ suffix[j] +")")
                            }
                        }
                    }  
                }
                // optional
                groups = dataH5ad.collect { 
                    k,v -> "group" in v ? v["group"] : 'NULL' 
                }
            } else {
                filePaths = dataH5ad.file_paths
                suffix = dataH5ad.suffix
                groups = null
            }
            data = data.concat(
                getFileChannel(
                    filePaths,
                    suffix,
                    groups
                ).map {
                    it -> tuple(
                        it[0], // sample ID
                        it[1], // file path
                        "h5ad", // input format
                        outputFileFormat,
                        it[2] // group
                    )
                }
            )
        }
        if(params.data.containsKey("loom")) {
            data = data.concat(
                getFileChannel( 
                    params.data.loom.file_paths,
                    params.data.loom.suffix,
                    'NULL'
                ).map {
                    it -> tuple(it[0], it[1], "loom", outputFileFormat, 'NULL')
                }
            )
        }
        if(params.data.containsKey("tsv")) {
            data = data.concat(
                getFileChannel( 
                    params.data.tsv.file_paths,
                    params.data.tsv.suffix,
                    'NULL'
                ).map {
                    it -> tuple(it[0], it[1], "tsv", outputFileFormat, 'NULL')
                }
            )
        }
        if(params.data.containsKey("csv")) {
            data = data.concat(
                getFileChannel( 
                    params.data.csv.file_paths,
                    params.data.csv.suffix,
                    'NULL'
                ).map {
                    it -> tuple(it[0], it[1], "csv", outputFileFormat, 'NULL')
                }
            )
        }
        if(params.data.containsKey("seurat_rds")) {
            data = data.concat(
                getFileChannel( 
                    params.data.seurat_rds.file_paths,
                    params.data.seurat_rds.suffix,
                    'NULL'
                ).map {
                    it -> tuple(it[0], it[1], "seurat_rds", outputFileFormat, 'NULL')
                }
            )
        }
        if(params.data.containsKey("fragments")) {
            data = data.concat(
                getFileChannelWithIndex(
                    params.data.fragments.file_paths,
                    params.data.fragments.suffix,
                    params.data.fragments.index_extension,
                    'NULL'
                ).map {
                    it -> tuple(it[0], it[1], "fragments", outputFileFormat, 'NULL')
                }
            )
        }
        if(params.data.containsKey("bam")) {
            data = data.concat(
                getFileChannelWithIndex(
                    params.data.bam.file_paths,
                    params.data.bam.suffix,
                    params.data.bam.index_extension,
                    'NULL'
                ).map {
                    it -> tuple(it[0], it[1], "bam", outputFileFormat, 'NULL' )
                }
            )
        }
        data.ifEmpty { exit 1, "Pipeline cannot run: no data provided." }

        if(!params.containsKey('quiet')) data.view()

    emit:
        data

}
