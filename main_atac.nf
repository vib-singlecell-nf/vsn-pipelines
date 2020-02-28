import static groovy.json.JsonOutput.*

nextflow.preview.dsl=2

if(!params.global.containsKey('seed')) {
    params.seed = workflow.manifest.version.replaceAll("\\.","").toInteger()

    Channel.from('').view {
            """
------------------------------------------------------------------
\u001B[32m No seed detected in the config \u001B[0m
\u001B[32m To ensure reproducibility the seed has been set to ${params.seed} \u001B[0m
------------------------------------------------------------------
            """
    }
}

def paramsCopy = params.findAll({!["parseConfig", "parse-config"].contains(it.key)})
params.manifestAsJSON = toJson(workflow.manifest)
params.paramsAsJSON = toJson(paramsCopy)

/* 
    ATAC-seq pipelines
*/


// runs mkfastq, then cellranger-atac count:
workflow cellranger_atac {

    include CELLRANGER_ATAC from './src/cellranger-atac/main.nf' params(params)
    CELLRANGER_ATAC(
        file(params.sc.cellranger_atac.mkfastq.csv),
        file(params.sc.cellranger_atac.mkfastq.runFolder),
        file(params.sc.cellranger_atac.count.reference)
    )

}

