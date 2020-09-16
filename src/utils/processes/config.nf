import java.nio.file.Paths
import groovy.transform.Memoized
import nextflow.script.ScriptBinding
import nextflow.config.ConfigParser
import static groovy.json.JsonOutput.*


def updateParams(params, resolvedParams, setter) {
    resolvedParams.each { k, v ->
        if(setter == null) {
            if(v instanceof Map) {
                if(!params.containsKey(k))
                    params."${k}" = [:]
                updateParams(params, v, params."${k}")
            } else {
                params."${k}" = v
            }
        } else {
            if(!setter.containsKey(k))
                setter."${k}" = [:]
            setter."${k}" = v instanceof Map ? updateParams(params, v, setter."${k}") : v
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

def includeConfig(Map params, String configRelativeFilePath) {
    def repoFilePath = workflow.scriptFile.getParent()
    def isMainRepo = repoFilePath.getName() == "vsn-pipelines"
    def config = new ConfigParser().setBinding([params: params])
    def co = new ConfigObject()
    def toolBaseDir = isMainRepo ? repoFilePath.toRealPath().toString() : repoFilePath.getParent().getParent().toRealPath().toString()
    config = config.parse(Paths.get(toolBaseDir, configRelativeFilePath))
    updateParams(params, config.params, null)
    return params
}
