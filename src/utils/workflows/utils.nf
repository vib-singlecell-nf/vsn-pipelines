nextflow.enable.dsl=2

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
} from "./../processes/utils.nf" params(params)

formatsAllowed = ['h5ad', 'loom']
taggedFilesToSkipPublishing = ['final_output']

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
        out = data

        if(fileOutputSuffix != null && !taggedFilesToSkipPublishing.any { fileOutputSuffix.contains(it) }) {
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
        }

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
        SC__PUBLISH.out

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

        if(!params.containsKey('quiet')) {
            Channel.from('').view {
                """
------------------------------------------------------------------
\u001B[32m No seed detected in the config \u001B[0m
\u001B[32m To ensure reproducibility the seed has been set to ${params.global.seed} \u001B[0m
------------------------------------------------------------------
                """
            }
        }
    } else {
        if(!params.containsKey('quiet')) {
            Channel.from('').view {
                """
------------------------------------------------------------------
\u001B[32m Custom seed detected in the config \u001B[0m
\u001B[32m Seed is set to ${params.global.seed} \u001B[0m
------------------------------------------------------------------
                """
            }
        }
        // If seed is of type String, it should be converted to an Integer because R doesn't not allow to have seeds of type character (see set.seed)
        if (params.global.seed instanceof String) {
            params.global.seed = params.global.seed.hashCode().abs()
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
    params.logDir = params.global.outdir + '/nextflow_log'
    includeConfig(params, 'src/utils/conf/workflow_report.config')

}
