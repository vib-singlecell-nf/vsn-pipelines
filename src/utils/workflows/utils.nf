nextflow.preview.dsl=2

import static groovy.json.JsonOutput.*

//////////////////////////////////////////////////////
//  Process imports:

include {
    includeConfig;
} from "./../processes/config.nf" params(params)
include {
    isParamNull;
    COMPRESS_HDF5;
    SC__PUBLISH;
    SC__PUBLISH as SC__PUBLISH_PROXY;
} from "./../processes/utils.nf" params(params)
include {
    SC__H5AD_CLEAN
} from "./../processes/h5adUpdate.nf" params(params)

formatsAllowed = ['h5ad', 'loom']
taggedFilesToClean = ['final_output']

//////////////////////////////////////////////////////
//  Define the workflow 

workflow PUBLISH {

    take: 
        data
        fileOutputSuffix
        fileOutputFormat
        toolName
        isParameterExplorationModeOn

    main:
        // Clean
        if(fileOutputSuffix != null && taggedFilesToClean.any { fileOutputSuffix.contains(it) } && isParameterExplorationModeOn) {
            if(!formatsAllowed.contains(fileOutputFormat))
                throw new Exception("The format " + fileOutputFormat + " is currently not allowed to be published.")
            if(fileOutputFormat == "h5ad") {
                out = SC__H5AD_CLEAN(
                    data
                )
            } else {
                out = data
            }
        } else {
            out = data
        }

        // Compress only if part of formatsAllowed
        if(fileOutputSuffix != null && formatsAllowed.any { fileOutputFormat.contains(it) }) {
            out = COMPRESS_HDF5(
                out.map {
                    // if stashedParams not there, just put null 3rd arg
                    it -> tuple(it[0], it[1], it.size() > 2 ? it[2]: null)
                },
                isParamNull(fileOutputSuffix) ? 'NULL' : fileOutputSuffix,
                isParamNull(toolName) ? 'NULL' : toolName,
            )
        }

        // Proxy to avoid file name collision
        SC__PUBLISH_PROXY(
            data.map {
                // if stashedParams not there, just put null 3rd arg
                it -> tuple(it[0], it[1], it.size() > 2 ? it[2]: null)
            },
            isParamNull(fileOutputSuffix) ? 'NULL' : fileOutputSuffix,
            isParamNull(toolName) ? 'NULL' : toolName,
            isParameterExplorationModeOn
        )

        // Publish
        SC__PUBLISH(
            out.map {
                // if stashedParams not there, just put null 3rd arg
                it -> tuple(it[0], it[1], it.size() > 2 ? it[2]: null)
            },
            isParamNull(fileOutputSuffix) ? 'NULL' : fileOutputSuffix,
            isParamNull(toolName) ? 'NULL' : toolName,
            isParameterExplorationModeOn
        )

    emit:
        SC__PUBLISH_PROXY.out

}

workflow COMBINE_BY_PARAMS {

    take:
        // Expects (sampleId, data, unstashedParams)
        A
        // Expects (sampleId, data, [stashedParams])
        B
        params

    main:
        if(params != null && params.isParameterExplorationModeOn()) {
            out = A.concat(
                B.map { it -> tuple(it[0], it[1], *it[2]) } // Unstash params
            ).map {
                it -> tuple(it[2..(it.size()-1)], it[0], it[1]) // Stash params
            }.groupTuple(
                by: [0, params.numParams()-1]
            ).map { 
                it -> tuple(it[1], *it[2], it[0]) 
            }
        } else {
            out = A.join(B)
        }

    emit:
        out

}

def setSeed(params) {
    if(!params.global.containsKey('seed')) {
        params.global.seed = workflow.manifest.version.replaceAll("\\.","").toInteger()

        Channel.from('').view {
                """
------------------------------------------------------------------
\u001B[32m No seed detected in the config \u001B[0m
\u001B[32m To ensure reproducibility the seed has been set to ${params.global.seed} \u001B[0m
------------------------------------------------------------------
                """
        }
    } else {
        Channel.from('').view {
                """
------------------------------------------------------------------
\u001B[32m Custom seed detected in the config \u001B[0m
\u001B[32m Seed is set to ${params.global.seed} \u001B[0m
------------------------------------------------------------------
                """
        }
    }
}

def INIT(params) {

    // Set the seed
    setSeed(params)
    if(!params.containsKey("misc") || !params.misc.containsKey("test")) {
        includeConfig(params, 'conf/test.config')
        params.misc.test.enabled = false
    }
    // Save manifest and params for notebook
    // Remove any closure attached to the config (this is for backward compatibility)
    def paramsCopy = params.findAll({!["parseConfig", "parse-config"].contains(it.key)})
    params.misc.manifestAsJSON = toJson(workflow.manifest)
    params.misc.paramsAsJSON = toJson(paramsCopy)
    // Include generic configs
    includeConfig(params, 'conf/generic.config')
    includeConfig(params, 'src/utils/conf/workflow_report.config')
    return params

}
