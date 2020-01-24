import java.nio.file.Paths
import groovy.transform.Memoized
import nextflow.script.ScriptBinding
import nextflow.config.ConfigParser
import static groovy.json.JsonOutput.*


def updateParams(params, resolvedParams, setter) {
    resolvedParams.each { k, v ->
        if(setter == null) {
            if(resolvedParams[k] instanceof Map)
                updateParams(params, resolvedParams[k], params."${k}")
        } else {
            setter."${k}" = resolvedParams[k] instanceof Map ? updateParams(params, resolvedParams[k], setter."${k}") : v
        }
    }
}

@Memoized
def resolveParams(Map params, boolean verbose) {
    if(!params.containsKey("strategy"))
        return params
    if(params.strategy != "min")
        return params
    def isRootDir = workflow.projectDir.getParent().getName() == "vib-singlecell-nf"
    def config = new ConfigParser().setBinding([params: params])
    def co = new ConfigObject()
    co.putAll(params)
    co.flatten().each { key, val ->
        if(key.endsWith("configVersion")) {
            // Extract the tool name based on the key
            def tool = key.split("\\.")[-2]
            // Build the path the versioned config of the current tool
            def toolBaseDir = isRootDir ? Paths.get(workflow.projectDir.toRealPath(), "src", tool) : workflow.projectDir.toRealPath()
            config = config.parse(Paths.get(toolBaseDir.toString(), "conf/min/base/${val}.config"))
        }
    }
    // Update the strategy since params has been resolved
    config.params.strategy = "max"
    updateParams(params, config.params, null)
    if(verbose)
        println(prettyPrint(toJson(params)))
    return params
}
