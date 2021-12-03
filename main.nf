import static groovy.json.JsonOutput.*

nextflow.enable.dsl=2

include { 
    INIT;
} from './src/utils/workflows/utils' params(params)

INIT(params)

include {
    SC__FILE_CONVERTER;
} from './src/utils/processes/utils' params(params)

include {
    getDataChannel;
} from './src/channels/channels' params(params)

// run multi-sample with bbknn, output a scope loom file
workflow bbknn {

    include {
        bbknn as BBKNN;
    } from './workflows/bbknn' params(params)
    include {
        PUBLISH as PUBLISH_SCOPE;
        PUBLISH as PUBLISH_SCANPY;
    } from "./src/utils/workflows/utils" params(params)

    getDataChannel | BBKNN

    if(params.utils?.publish) {
        PUBLISH_SCOPE(
            BBKNN.out.scopeloom,
            "BBKNN",
            "loom",
            null,
            false
        )
        PUBLISH_SCANPY(
            BBKNN.out.scanpyh5ad,
            "BBKNN",
            "h5ad",
            null,
            false
        )
    }
}

// run multi-sample with mnncorrect, output a scope loom file
workflow mnncorrect {

    include {
        mnncorrect as MNNCORRECT;
    } from './workflows/mnncorrect' params(params)
    include {
        PUBLISH as PUBLISH_SCOPE;
        PUBLISH as PUBLISH_SCANPY;
    } from "./src/utils/workflows/utils" params(params)

    getDataChannel | MNNCORRECT

    if(params.utils?.publish) {
        PUBLISH_SCOPE(
            MNNCORRECT.out.scopeloom,
            "MNNCORRECT",
            "loom",
            null,
            false
        )
        PUBLISH_SCANPY(
            MNNCORRECT.out.scanpyh5ad,
            "MNNCORRECT",
            "h5ad",
            null,
            false
        )
    }

}

// run multi-sample with bbknn, output a scope loom file
workflow harmony {

    include {
        harmony as HARMONY;
    } from './workflows/harmony' params(params)
    include {
        PUBLISH as PUBLISH_SCOPE;
        PUBLISH as PUBLISH_SCANPY;
    } from "./src/utils/workflows/utils" params(params)

    batchVariables = params.tools.harmony.varsUse
    outputSuffix = params.utils?.publish?.annotateWithBatchVariableName ? "HARMONY" + "_BY_" + batchVariables.join("_").toUpperCase() : "HARMONY"

    getDataChannel | HARMONY

    if(params.utils?.publish) {
        PUBLISH_SCOPE(
            HARMONY.out.scopeloom,
            outputSuffix,
            "loom",
            null,
            false
        )
        PUBLISH_SCANPY(
            HARMONY.out.scanpyh5ad,
            outputSuffix,
            "h5ad",
            null,
            false
        )
    }

}

workflow harmony_only {

    include {
        BEC_HARMONY as HARMONY;
    } from './src/harmony/workflows/harmony_only' params(params)
    include {
        PUBLISH as PUBLISH_HARMONY;
    } from "./src/utils/workflows/utils" params(params)
    

    batchVariables = params.tools.harmony.varsUse
    outputSuffix = params.utils?.publish?.annotateWithBatchVariableName ? "HARMONY" + "_BY_" + batchVariables.join("_").toUpperCase() : "HARMONY"

    getDataChannel | HARMONY

    if(params.utils?.publish) {
        PUBLISH_HARMONY(
            HARMONY.out,
            params.utils?.publish?.annotateWithBatchVariableName ? "HARMONY" + "_BY_" +  batchVariables.join("_").toUpperCase() : "HARMONY",
            "h5ad",
            null,
            false
        )
    }

}

// run multi-sample with bbknn, then scenic from the filtered output:
workflow bbknn_scenic {

    include {
        bbknn as BBKNN;
    } from './workflows/bbknn' params(params)
    include {
        scenic_append as SCENIC_APPEND; 
    } from './src/scenic/main' params(params)
    include {
        PUBLISH as PUBLISH_BBKNN;
        PUBLISH as PUBLISH_BBKNN_SCENIC;
    } from "./src/utils/workflows/utils" params(params)

    getDataChannel | BBKNN

    if(params.utils?.publish) {
        PUBLISH_BBKNN(
            BBKNN.out.scanpyh5ad,
            "BBKNN",
            "h5ad",
            null,
            false
        )
    }

    SCENIC_APPEND(
        BBKNN.out.filteredloom,
        BBKNN.out.scopeloom
    )

    if(params.utils?.publish) {
        PUBLISH_BBKNN_SCENIC(
            SCENIC_APPEND.out,
            "BBKNN_SCENIC",
            "loom",
            null,
            false
        )
    }

}

// run multi-sample with harmony, then scenic from the filtered output:
workflow harmony_scenic {

    include {
        harmony as HARMONY;
    } from './workflows/harmony' params(params)
    include {
        scenic_append as SCENIC_APPEND;
    } from './src/scenic/main' params(params)
    include {
        PUBLISH as PUBLISH_HARMONY;
        PUBLISH as PUBLISH_HARMONY_SCENIC;
    } from "./src/utils/workflows/utils" params(params)

    batchVariables = params.tools.harmony.varsUse
    outputSuffix = params.utils?.publish?.annotateWithBatchVariableName ? "HARMONY" + "_BY_" + batchVariables.join("_").toUpperCase() : "HARMONY"

    getDataChannel | HARMONY

    if(params.utils?.publish) {
        PUBLISH_HARMONY(
            HARMONY.out.scanpyh5ad,
            outputSuffix,
            "h5ad",
            null,
            false
        )
    }

    SCENIC_APPEND( 
        HARMONY.out.filteredloom,
        HARMONY.out.scopeloom 
    )

    if(params.utils?.publish) {
        PUBLISH_HARMONY_SCENIC(
            SCENIC_APPEND.out,
            outputSuffix + "_SCENIC",
            "loom",
            null,
            false
        )
    }

}


// run single_sample, output a scope loom file
workflow single_sample {

    include {
        single_sample as SINGLE_SAMPLE;
    } from './workflows/single_sample' params(params)
    include {
        PUBLISH as PUBLISH_SINGLE_SAMPLE_SCOPE;
        PUBLISH as PUBLISH_SINGLE_SAMPLE_SCANPY;
    } from "./src/utils/workflows/utils" params(params)

    getDataChannel | SINGLE_SAMPLE

    if(params.utils?.publish) {
        PUBLISH_SINGLE_SAMPLE_SCOPE(
            SINGLE_SAMPLE.out.scopeloom,
            "SINGLE_SAMPLE",
            "loom",
            null,
            false
        )
        PUBLISH_SINGLE_SAMPLE_SCANPY(
            SINGLE_SAMPLE.out.scanpyh5ad,
            "SINGLE_SAMPLE",
            "h5ad",
            null,
            false
        )
    }  

}

// run single_sample QC
workflow single_sample_qc {

    include {
        single_sample_qc as SINGLE_SAMPLE_QC;
    } from './src/scanpy/main' params(params)
    include {
        PUBLISH;
    } from "./src/utils/workflows/utils" params(params)

    getDataChannel | SINGLE_SAMPLE_QC

    if(params.utils?.publish) {
        PUBLISH(
            SINGLE_SAMPLE_QC.out.filtered,
            "SINGLE_SAMPLE_QC",
            "h5ad",
            null,
            false
        )
    }  

}

workflow multi_sample_qc {

    include {
        multi_sample_qc as MULTI_SAMPLE_QC;
    } from './src/scanpy/main' params(params)
    include {
        PUBLISH;
    } from "./src/utils/workflows/utils" params(params)

    getDataChannel | MULTI_SAMPLE_QC

    if(params.utils?.publish) {
        PUBLISH(
            MULTI_SAMPLE_QC.out.filtered,
            "MULTI_SAMPLE_QC",
            "h5ad",
            null,
            false
        )
    }  

}

workflow multi_sample {

    include {
        multi_sample as MULTI_SAMPLE;
    } from './workflows/multi_sample' params(params)

    getDataChannel | MULTI_SAMPLE
    include {
        PUBLISH as PUBLISH_SCOPE;
        PUBLISH as PUBLISH_SCANPY;
    } from "./src/utils/workflows/utils" params(params)

    if(params.utils?.publish) {
        PUBLISH_SCOPE(
            MULTI_SAMPLE.out.scopeloom,
            "MULTI_SAMPLE",
            "loom",
            null,
            false
        )
        PUBLISH_SCANPY(
            MULTI_SAMPLE.out.scanpyh5ad,
            "MULTI_SAMPLE",
            "h5ad",
            null,
            false
        )
    }

}

// run single_sample, then scenic from the filtered output:
workflow single_sample_scenic {

    include {
        scenic_append as SCENIC_APPEND;
    } from './src/scenic/main' params(params)
    include {
        single_sample as SINGLE_SAMPLE;
    } from './workflows/single_sample' params(params)
    include {
        PUBLISH as PUBLISH_SCANPY;
        PUBLISH as PUBLISH_SCOPE;
    } from "./src/utils/workflows/utils" params(params)

    getDataChannel | SINGLE_SAMPLE

    if(params.utils?.publish) {
        PUBLISH_SCANPY(
            SINGLE_SAMPLE.out.scanpyh5ad,
            "SINGLE_SAMPLE",
            "h5ad",
            null,
            false
        )
    }

    SCENIC_APPEND(
        SINGLE_SAMPLE.out.filteredloom,
        SINGLE_SAMPLE.out.scopeloom
    )

    if(params.utils?.publish) {
        PUBLISH_SCOPE(
            SCENIC_APPEND.out,
            "SINGLE_SAMPLE_SCENIC",
            "loom",
            null,
            false
        )
    }

}

workflow pcacv {

    include {
        PCACV__FIND_OPTIMAL_NPCS;
    } from './src/pcacv/processes/runPCACV' params(params)
    getDataChannel().map {
        it -> tuple(it[0], it[1])
    }
    .set{ data }
    PCACV__FIND_OPTIMAL_NPCS(data)
}

workflow single_sample_scrublet {

    include {
        SINGLE_SAMPLE as SCANPY__SINGLE_SAMPLE;
    } from './src/scanpy/workflows/single_sample' params(params)
    include {
        DOUBLET_REMOVAL as SCRUBLET__DOUBLET_REMOVAL;
    } from "./src/scrublet/workflows/doublet_removal" params(params)
    include {
        ANNOTATE_BY_CELL_METADATA;
    } from './src/utils/workflows/annotateByCellMetadata.nf' params(params)
    include {
        PUBLISH as PUBLISH_SINGLE_SAMPLE_SCRUBLET;
    } from './src/utils/workflows/utils.nf' params(params)
    include {
        SC__H5AD_TO_LOOM
    } from './src/utils/processes/h5adToLoom.nf' params(params)

    data = getDataChannel | SC__FILE_CONVERTER
    SCANPY__SINGLE_SAMPLE( data )
    SCRUBLET__DOUBLET_REMOVAL(
        data.join( SCANPY__SINGLE_SAMPLE.out.dr_pca_data ),
        SCANPY__SINGLE_SAMPLE.out.final_processed_data
    )
    // Annotate the final processed file with doublet information inferred from Scrublet
    ANNOTATE_BY_CELL_METADATA(
        SCANPY__SINGLE_SAMPLE.out.final_processed_data.map {
            it -> tuple(it[0], it[1])
        },
        SCRUBLET__DOUBLET_REMOVAL.out.doublet_detection.map {
            it -> tuple(it[0], it[1])
        },
        "scrublet"
    )
    SC__H5AD_TO_LOOM(
        SCANPY__SINGLE_SAMPLE.out.filtered_data.map {
            it -> tuple(it[0], it[1])
        }.join(
            ANNOTATE_BY_CELL_METADATA.out
        )
    )

    if(params.utils?.publish) {
        PUBLISH_SINGLE_SAMPLE_SCRUBLET(
            SC__H5AD_TO_LOOM.out,
            "SINGLE_SAMPLE_SCRUBLET",
            "loom",
            null,
            false
        )
    }

}

workflow decontx {

    include {
        decontx as CELDA__DECONTX;
    } from "./src/celda/main" params(params)
    // Run DecontX on the data
    CELDA__DECONTX()

}

workflow single_sample_decontx {

    include {
        SINGLE_SAMPLE as SCANPY__SINGLE_SAMPLE;
    } from './src/scanpy/workflows/single_sample' params(params)
    include {
        decontx as CELDA__DECONTX;
    } from "./src/celda/main" params(params)
    include {
        ANNOTATE_BY_CELL_METADATA_BY_PAIR;
    } from './src/utils/workflows/annotateByCellMetadata.nf' params(params)
    include {
        PUBLISH;
    } from './src/utils/workflows/utils.nf' params(params)
    include {
        SC__H5AD_TO_LOOM;
    } from './src/utils/processes/h5adToLoom.nf' params(params)

    data = getDataChannel \
        | SC__FILE_CONVERTER
    // Run Single-sample pipeline on the data
    SCANPY__SINGLE_SAMPLE( data )
    // Run DecontX on the data
    CELDA__DECONTX()

    // Annotate the final processed file with doublet information inferred from Scrublet
    ANNOTATE_BY_CELL_METADATA_BY_PAIR(
        SCANPY__SINGLE_SAMPLE.out.final_processed_data,
        CELDA__DECONTX.out.outlier_table,
        "celda.decontx"
    )
    SC__H5AD_TO_LOOM(
        SCANPY__SINGLE_SAMPLE.out.filtered_data.map {
            it -> tuple(it[0], it[1])
        }.join(
            ANNOTATE_BY_CELL_METADATA_BY_PAIR.out
        )
    )
    if(params.utils?.publish) {
        PUBLISH(
            SC__H5AD_TO_LOOM.out,
            "SINGLE_SAMPLE_CELDA_DECONTX_"+ params.tools.celda.decontx.strategy.toUpperCase(),
            "loom",
            null,
            false
        )
    }

}

workflow single_sample_decontx_scrublet {
    include {
        SINGLE_SAMPLE as SCANPY__SINGLE_SAMPLE;
    } from './src/scanpy/workflows/single_sample' params(params)
    include {
        decontx as CELDA__DECONTX;
    } from "./src/celda/main" params(params)
    include {
        DOUBLET_REMOVAL as SCRUBLET__DOUBLET_REMOVAL;
    } from "./src/scrublet/workflows/doublet_removal" params(params)
    include {
        ANNOTATE_BY_CELL_METADATA_BY_PAIR;
    } from './src/utils/workflows/annotateByCellMetadata.nf' params(params)
    include {
        PUBLISH as PUBLISH_SINGLE_SAMPLE_ANNOTATED;
        PUBLISH as PUBLISH_CELDA_DECONTX_SCRUBLET
    } from './src/utils/workflows/utils.nf' params(params)
    include {
        SC__H5AD_TO_LOOM;
    } from './src/utils/processes/h5adToLoom.nf' params(params)

    data = getDataChannel \
        | SC__FILE_CONVERTER
    // Run Single-sample pipeline on the INPUT DATA
    SCANPY__SINGLE_SAMPLE( data )
    // Run DecontX on the INPUT DATA
    CELDA__DECONTX()
    // Run Scrublet on the DecontX FILTERED DATA
    SCRUBLET__DOUBLET_REMOVAL(
        CELDA__DECONTX.out.decontx_processed.join( SCANPY__SINGLE_SAMPLE.out.dr_pca_data ),
        SCANPY__SINGLE_SAMPLE.out.final_processed_data
    )

    // Annotate the final processed file with doublet information inferred from Scrublet
    ANNOTATE_BY_CELL_METADATA_BY_PAIR(
        SCANPY__SINGLE_SAMPLE.out.final_processed_data,
        SCRUBLET__DOUBLET_REMOVAL.out.doublet_detection.map {
            it -> tuple(it[0], it[1])
        }.mix(
            CELDA__DECONTX.out.outlier_table
        ).groupTuple(),
        "scrublet" // Works because decontx requires the same params in cell_annotate
    )
    SC__H5AD_TO_LOOM(
        SCANPY__SINGLE_SAMPLE.out.filtered_data.map {
            it -> tuple(it[0], it[1])
        }.join(
            ANNOTATE_BY_CELL_METADATA_BY_PAIR.out
        )
    )
    if(params.utils?.publish) {
        // Publish loom file where INPUT DATA:
        // - processed by single_sample pipeline and ->
        // - annotated by DecontX potential outliers
        // - annotated by Scrublet potential doublets
        PUBLISH_SINGLE_SAMPLE_ANNOTATED(
            SC__H5AD_TO_LOOM.out,
            "SINGLE_SAMPLE_ANNOTATED",
            "loom",
            null,
            false
        )
        // Publish h5ad file where INPUT DATA:
        // - processed by DecontX (either through filter or correct strategy) and ->
        // - potential doublets removed by Scrublet 
        PUBLISH_CELDA_DECONTX_SCRUBLET(
            SCRUBLET__DOUBLET_REMOVAL.out.data_doublets_removed,
            "CELDA_DECONTX_"+ params.tools.celda.decontx.strategy.toUpperCase() +"_SCRUBLET",
            "h5ad",
            null,
            false
        )
    }
}

workflow soupx {

    include {
        soupx as SOUPX__DECONTX;
    } from "./src/soupx/main" params(params)
    // Run DecontX on the data
    SOUPX__DECONTX()

}

workflow single_sample_soupx {

    include {
        SINGLE_SAMPLE as SCANPY__SINGLE_SAMPLE;
    } from './src/scanpy/workflows/single_sample' params(params)
    include {
        soupx as SOUPX;
    } from "./src/soupx/main" params(params)
    include {
        ANNOTATE_BY_CELL_METADATA_BY_PAIR;
    } from './src/utils/workflows/annotateByCellMetadata.nf' params(params)
    include {
        PUBLISH;
    } from './src/utils/workflows/utils.nf' params(params)
    include {
        SC__H5AD_TO_LOOM;
    } from './src/utils/processes/h5adToLoom.nf' params(params)

    data = getDataChannel \
        | SC__FILE_CONVERTER
    // Run Single-sample pipeline on the INPUT DATA:
    SCANPY__SINGLE_SAMPLE( data )
    // Run SoupX on the INPUT DATA:
    SOUPX()

    SC__H5AD_TO_LOOM(
        SCANPY__SINGLE_SAMPLE.out.filtered_data.map {
            it -> tuple(it[0], it[1])
        }
    )
    if(params.utils?.publish) {
        PUBLISH(
            SC__H5AD_TO_LOOM.out,
            "SINGLE_SAMPLE",
            "loom",
            null,
            false
        )
    }

}

workflow single_sample_soupx_scrublet {
    include {
        SINGLE_SAMPLE as SCANPY__SINGLE_SAMPLE;
    } from './src/scanpy/workflows/single_sample' params(params)
    include {
        soupx as SOUPX;
    } from "./src/soupx/main" params(params)
    include {
        DOUBLET_REMOVAL as SCRUBLET__DOUBLET_REMOVAL;
    } from "./src/scrublet/workflows/doublet_removal" params(params)
    include {
        ANNOTATE_BY_CELL_METADATA_BY_PAIR;
    } from './src/utils/workflows/annotateByCellMetadata.nf' params(params)
    include {
        PUBLISH as PUBLISH_SINGLE_SAMPLE_ANNOTATED;
        PUBLISH as PUBLISH_SOUPX_SCRUBLET
    } from './src/utils/workflows/utils.nf' params(params)
    include {
        SC__H5AD_TO_LOOM;
    } from './src/utils/processes/h5adToLoom.nf' params(params)

    data = getDataChannel \
        | SC__FILE_CONVERTER
    // Run Single-sample pipeline on the INPUT DATA
    SCANPY__SINGLE_SAMPLE( data )
    // Run Soupx on the INPUT DATA
    SOUPX()
    // Run Scrublet on the DecontX FILTERED DATA
    SCRUBLET__DOUBLET_REMOVAL(
        SOUPX.out.soupx_processed.join( SCANPY__SINGLE_SAMPLE.out.dr_pca_data ),
        SCANPY__SINGLE_SAMPLE.out.final_processed_data
    )

    // Annotate the final processed file with doublet information inferred from Scrublet
    ANNOTATE_BY_CELL_METADATA_BY_PAIR(
        SCANPY__SINGLE_SAMPLE.out.final_processed_data,
        SCRUBLET__DOUBLET_REMOVAL.out.doublet_detection.map {
            it -> tuple(it[0], it[1])
        },
        "scrublet"
    )
    SC__H5AD_TO_LOOM(
        SCANPY__SINGLE_SAMPLE.out.filtered_data.map {
            it -> tuple(it[0], it[1])
        }.join(
            ANNOTATE_BY_CELL_METADATA_BY_PAIR.out
        )
    )
    if(params.utils?.publish) {
        // Publish loom file where INPUT DATA:
        // - processed by single_sample pipeline and ->
        // - annotated by Scrublet potential doublets
        PUBLISH_SINGLE_SAMPLE_ANNOTATED(
            SC__H5AD_TO_LOOM.out,
            "SINGLE_SAMPLE_ANNOTATED",
            "loom",
            null,
            false
        )
        // Publish h5ad file where INPUT DATA:
        // - processed by SoupX and ->
        // - potential doublets removed by Scrublet 
        PUBLISH_SOUPX_SCRUBLET(
            SCRUBLET__DOUBLET_REMOVAL.out.data_doublets_removed,
            "SOUPX_SCRUBLET",
            "h5ad",
            null,
            false
        )
    }
}

// run single_sample, then scenic from the previous input (not standalone):
workflow pipe_single_sample_scenic {

    take:
        data
    main:
        include {
            scenic_append as SCENIC_APPEND;
        } from './src/scenic/main' params(params)
        include {
            single_sample as SINGLE_SAMPLE;
        } from './workflows/single_sample' params(params)
        include {
            PUBLISH as PUBLISH_P_SINGLE_SAMPLE_SCENIC
        } from "./src/utils/workflows/utils" params(params)

        data | SINGLE_SAMPLE
        SCENIC_APPEND(
            SINGLE_SAMPLE.out.filteredloom,
            SINGLE_SAMPLE.out.scopeloom
        )

        if(params.utils?.publish) {
            PUBLISH_P_SINGLE_SAMPLE_SCENIC(
                SCENIC_APPEND.out,
                "P_SINGLE_SAMPLE_SCENIC",
                "loom",
                null,
                false
            )
        }

}


// run scenic directly from an existing loom file:
workflow scenic {

    include {
        scenic as SCENIC;
    } from './src/scenic/main' params(params)
    include {
        PUBLISH as PUBLISH_SCENIC;
    } from "./src/utils/workflows/utils" params(params)

    SCENIC( 
        Channel.of( tuple(params.global.project_name, file(params.tools.scenic.filteredLoom))) 
    )

    if(params.utils?.publish) {
        PUBLISH_SCENIC(
            SCENIC.out,
            "SCENIC",
            "loom",
            null,
            false
        )
    }

}


// runs mkfastq, then CellRanger count:
workflow cellranger {

    include {
        CELLRANGER;
    } from './src/cellranger/main' params(params)

    CELLRANGER(
        file(params.tools.cellranger.mkfastq.csv),
        file(params.tools.cellranger.mkfastq.runFolder),
        file(params.tools.cellranger.count.transcriptome)
    )

    emit:
        CELLRANGER.out
}

workflow cellranger_libraries {

    include {
        CELLRANGER_LIBRARIES;
    } from './src/cellranger/workflows/cellranger_libraries' params(params)

    CELLRANGER_LIBRARIES(
        file(params.tools.cellranger.mkfastq.csv),
        file(params.tools.cellranger.mkfastq.runFolder),
        file(params.tools.cellranger.count.transcriptome),
        file(params.tools.cellranger.count.featureRef)
    )

    emit:
        CELLRANGER_LIBRARIES.out

}

workflow cellranger_count_metadata {

    include {
        CELLRANGER_COUNT_WITH_METADATA;
    } from './src/cellranger/workflows/cellRangerCountWithMetadata' params(params)

    CELLRANGER_COUNT_WITH_METADATA(
        file(params.tools.cellranger.count.transcriptome),
        file(params.tools.cellranger.count.metadata)
    )
    emit:
        CELLRANGER_COUNT_WITH_METADATA.out

}

workflow cellranger_count_metadata_single_sample_scenic {

    cellranger_count_metadata | \
        map {
            it -> tuple(it[0], it[1], "10x_cellranger_mex", "h5ad")
        } | \
        pipe_single_sample_scenic

}

workflow cellranger_count_libraries {

    include {
        CELLRANGER_COUNT_WITH_LIBRARIES;
    } from './src/cellranger/workflows/cellRangerCountWithLibraries' params(params)

    CELLRANGER_COUNT_WITH_LIBRARIES(
        file(params.tools.cellranger.count.transcriptome),
        file(params.tools.cellranger.count.featureRef),
        params.tools.cellranger.count.libraries
    )

    emit:
        CELLRANGER_COUNT_WITH_LIBRARIES.out

}

workflow cellranger_count_demuxlet {
    include {
        get_bam_barcodes_from_cellranger_rna;
        DEMUXLET;
    } from './src/popscle/workflows/demuxlet.nf' params(params)
    include {
        SC__CELLRANGER__COUNT as CELLRANGER_COUNT;
    } from './src/cellranger/processes/count'
    if (params.tools.cellranger.count.fastqs instanceof Map) {
        // Remove default key
        Channel.from(params.tools.cellranger.count.fastqs.findAll {
            it.key != 'default' 
        }.collect { k, v -> 
            // Split possible multiple file paths
            if(v.contains(',')) {
                v = Arrays.asList(v.split(',')).collect { fqs -> file(fqs) }
            } else {
                v = file(v)
            }
            tuple(k, v) 
        })
        // Group fastqs per sample
        .groupTuple()
        .map {
            it -> tuple(it[0], *it[1])
        }
        .set { fastq_data }           
    }
    data = CELLRANGER_COUNT(
        params.tools.cellranger.count.transcriptome,
        fastq_data
    )
    get_bam_barcodes_from_cellranger_rna(data) |
        DEMUXLET
}

workflow freemuxlet {
    include {
        data_channel_to_bam_barcodes;
        FREEMUXLET;
    } from './src/popscle/workflows/demuxlet.nf' params(params)
    
    getDataChannel |
        data_channel_to_bam_barcodes |
        FREEMUXLET
}

workflow demuxlet {
    include {
        data_channel_to_bam_barcodes;
        DEMUXLET;
    } from './src/popscle/workflows/demuxlet.nf' params(params)

    getDataChannel |
        data_channel_to_bam_barcodes |
        DEMUXLET
}

// runs mkfastq, CellRanger count, then single_sample:
workflow single_sample_cellranger {

    include {
        single_sample as SINGLE_SAMPLE;
    } from './workflows/single_sample' params(params)

    data = cellranger()
    SINGLE_SAMPLE(
        data.map {
            tuple(it[0], it[1], "10x_cellranger_mex", "h5ad")
            }
    )

}

workflow cellranger_multi_sample {

    include { 
        multi_sample as MULTI_SAMPLE;
    } from './workflows/multi_sample' params(params)

    data = cellranger()
    MULTI_SAMPLE(
        data.map {
            tuple(it[0], it[1], "10x_cellranger_mex", "h5ad")
        }
    )

}

workflow cellranger_multi_sample_demuxlet {

    include {
        multi_sample as MULTI_SAMPLE;
    } from './workflows/multi_sample' params(params)
    include {
        get_bam_barcodes_from_cellranger_rna;
        DEMUXLET;
    } from './src/popscle/workflows/demuxlet.nf' params(params)

    data = cellranger()
    MULTI_SAMPLE(        
        data.map {
            tuple(it[0], it[1], "10x_cellranger_mex", "h5ad")
        }
    )
    get_bam_barcodes_from_cellranger_rna(data) |
        DEMUXLET

}

workflow cellranger_libraries_multi_sample {

    include {
        multi_sample as MULTI_SAMPLE;
    } from './workflows/multi_sample' params(params)

    data = cellranger_libraries()
    MULTI_SAMPLE(        
        data.map {
            tuple(it[0], it[1], "10x_cellranger_mex", "h5ad")
        }
    )
}

workflow cellranger_libraries_freemuxlet_multi_sample {

    include {
        multi_sample as MULTI_SAMPLE;
    } from './workflows/multi_sample' params(params)
    include {
        get_bam_barcodes_from_cellranger_rna;
        FREEMUXLET;
    } from './src/popscle/workflows/demuxlet.nf' params(params)

    data = cellranger_libraries()
    MULTI_SAMPLE(
        data.map {
            tuple(it[0], it[1], "10x_cellranger_mex", "h5ad")
            }
    )
    get_bam_barcodes_from_cellranger_rna(data) |
        FREEMUXLET

}

workflow cellranger_libraries_demuxlet_multi_sample {

    include {
        multi_sample as MULTI_SAMPLE;
    } from './workflows/multi_sample' params(params)
    include {
        get_bam_barcodes_from_cellranger_rna;
        DEMUXLET;
    } from './src/popscle/workflows/demuxlet.nf' params(params)

    data = cellranger_libraries()
    MULTI_SAMPLE(
        data.map {
            tuple(it[0], it[1], "10x_cellranger_mex", "h5ad")
            }
    )
    get_bam_barcodes_from_cellranger_rna(data) |
        DEMUXLET
}

workflow star {

    include {
        star as STAR;
    } from './workflows/star' params(params)
    STAR()

}


workflow single_sample_star {

    include {
        single_sample_star as SINGLE_SAMPLE_STAR;
    } from './workflows/single_sample_star' params(params)

    SINGLE_SAMPLE_STAR()

}

workflow nemesh {

    include {
        nemesh as NEMESH;
    } from './workflows/nemesh' params(params)

    NEMESH()

}

workflow sra {

    main:
        include {
            getChannel as getSRAChannel;
        } from './src/channels/sra' params(params)
        include {
            DOWNLOAD_FROM_SRA;
        } from './src/utils/workflows/downloadFromSRA' params(params)
        include {
            PUBLISH;
        } from "./src/utils/workflows/utils" params(params)
        
        // Run 
        out = DOWNLOAD_FROM_SRA( getSRAChannel( params.data.sra ) )
        if(params.utils?.publish) {
            PUBLISH(
                out.transpose(),
                null,
                null,
                "sra",
                false
            )
        }
    emit:
        out

}

workflow sra_cellranger_bbknn {

    main:
        include {
            SC__CELLRANGER__PREPARE_FOLDER;
            SC__CELLRANGER__COUNT;
        } from './src/cellranger/processes/utils' params(params)
        include {
            bbknn as BBKNN;
        } from './workflows/bbknn' params(params)
        
        // Run 
        out = sra()
        SC__CELLRANGER__PREPARE_FOLDER( out.groupTuple() )
        SC__CELLRANGER__COUNT(
            file(params.tools.cellranger.count.transcriptome),
            SC__CELLRANGER__PREPARE_FOLDER.out
        )
        BBKNN( 
            SC__CELLRANGER__COUNT.out.map {
                it -> tuple(it[0], it[1], "10x_cellranger_mex", "h5ad")
            }
        )

    emit:
        filteredLoom = BBKNN.out.filteredloom
        scopeLoom = BBKNN.out.scopeloom

}

workflow sra_cellranger_bbknn_scenic {

    include {
        scenic_append as SCENIC_APPEND;
    } from './src/scenic/main' params(params)
    include {
        PUBLISH as PUBLISH_SRA_CELLRANGER_BBKNN_SCENIC;
    } from "./src/utils/workflows/utils" params(params)

    sra_cellranger_bbknn()
    SCENIC_APPEND(
        sra_cellranger_bbknn.out.filteredLoom,
        sra_cellranger_bbknn.out.scopeLoom
    )

    if(params.utils?.publish) {
        PUBLISH_SRA_CELLRANGER_BBKNN_SCENIC(
            SCENIC_APPEND.out,
            "SRA_CELLRANGER_BBKNN_SCENIC",
            "loom",
            null,
            false
        )
    }  

}

/**
 * Utility workflows
 */

workflow cell_annotate {

    include {
        STATIC__ANNOTATE_BY_CELL_METADATA as ANNOTATE_BY_CELL_METADATA;
    } from './src/utils/workflows/annotateByCellMetadata' params(params)
    include {
        PUBLISH;
    } from "./src/utils/workflows/utils" params(params)

    // Run
    getDataChannel | \
        SC__FILE_CONVERTER
    ANNOTATE_BY_CELL_METADATA( 
        SC__FILE_CONVERTER.out, 
        null,
    )
    PUBLISH(
        ANNOTATE_BY_CELL_METADATA.out,
        "ANNOTATE_BY_CELL_METADATA",
        "h5ad",
        "utils",
        false
    )

}

workflow cell_annotate_filter {

    main:
        _cell_annotate_filter(true)

}

workflow _cell_annotate_filter {

    take:
        // Expects publish : boolean
        publish
    main:
        include {
            STATIC__ANNOTATE_BY_CELL_METADATA as ANNOTATE_BY_CELL_METADATA;
        } from './src/utils/workflows/annotateByCellMetadata' params(params)
        include {
            FILTER_BY_CELL_METADATA
        } from './src/utils/workflows/filterByCellMetadata' params(params)
        include {
            PUBLISH as PUBLISH_H5AD_CELL_ANNOTATED;
            PUBLISH as PUBLISH_H5AD_CELL_FILTERED;
        } from "./src/utils/workflows/utils" params(params)

        // Run
        getDataChannel | \
            SC__FILE_CONVERTER

        if(!params.utils?.cell_annotate)
            throw new Exception("VSN ERROR: The cell_annotate param is missing in params.utils.")

        // Annotate & publish
        ANNOTATE_BY_CELL_METADATA( 
            SC__FILE_CONVERTER.out, 
            null,
        )
        if(params.utils.cell_annotate.containsKey("publish") && params.utils.cell_annotate.publish) {
            PUBLISH_H5AD_CELL_ANNOTATED(
                ANNOTATE_BY_CELL_METADATA.out,
                "ANNOTATE_BY_CELL_METADATA",
                "h5ad",
                "utils",
                false
            )
        }

        if(!params.utils?.cell_filter)
            throw new Exception("VSN ERROR: The cell_filter param is missing in params.utils.")

        // Filter (& clean) & publish
        FILTER_BY_CELL_METADATA(
            ANNOTATE_BY_CELL_METADATA.out,
            null
        )

        if(params.utils.cell_filter?.publish) {
            PUBLISH_H5AD_CELL_FILTERED(
                FILTER_BY_CELL_METADATA.out,
                "FILTER_BY_CELL_METADATA",
                "h5ad",
                "utils",
                false
            )
        }
        if(params.utils?.publish && publish) {
            PUBLISH_H5AD_CELL_FILTERED(
                FILTER_BY_CELL_METADATA.out,
                "CELL_ANNOTATE_FILTER",
                "h5ad",
                "utils",
                false
            )
        }
    emit:
        out = FILTER_BY_CELL_METADATA.out

}

workflow cell_annotate_filter_and_sample_annotate {

    main:
        // Import
        include {
            SC__H5AD_BEAUTIFY;
        } from './src/utils/processes/h5adUpdate.nf' params(params)
        include {
            hasMetadataFilePath;
            SC__ANNOTATE_BY_SAMPLE_METADATA
        } from './src/utils/processes/h5adAnnotate.nf' params(params)
        include {
            PUBLISH;
        } from "./src/utils/workflows/utils" params(params)


        // Run
        out = _cell_annotate_filter(false)

        // Annotate cells based on an indexed sample-based metadata table
        if(!params.utils?.sample_annotate)
            throw new Exception("VSN ERROR: The sample_annotate param is missing in params.utils.")

        if (!hasMetadataFilePath(params.utils.sample_annotate)) {
            throw new Exception("VSN ERROR: The metadataFilePath param is missing in sample_annotate.")
        }
        out = SC__ANNOTATE_BY_SAMPLE_METADATA( out )

        if(params.utils.file_cleaner) {
            out = SC__H5AD_BEAUTIFY( out )
        }

        if(params.utils?.publish) {
            PUBLISH(
                out,
                "CELL_ANNOTATE_FILTER_AND_SAMPLE_ANNOTATE",
                "h5ad",
                "utils",
                false
            )
        }

}

workflow filter_and_annotate_and_clean {

    include {
        FILTER_AND_ANNOTATE_AND_CLEAN
    } from './src/utils/workflows/filterAnnotateClean' params(params)
    include {
        PUBLISH;
    } from "./src/utils/workflows/utils" params(params)

    // Run
    getDataChannel | \
        SC__FILE_CONVERTER | \
        FILTER_AND_ANNOTATE_AND_CLEAN

    PUBLISH(
        FILTER_AND_ANNOTATE_AND_CLEAN.out,
        "FILTER_AND_ANNOTATE_AND_CLEAN",
        "h5ad",
        "utils",
        false
    )

}
