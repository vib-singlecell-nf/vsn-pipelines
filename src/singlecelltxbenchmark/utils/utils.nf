// Could be optimized
def groupParams(params) {
    paramsMapGroupedByProcess = [:]

    params.each {
        if(it.key.contains('___')) {
            keySplitted = it.key.split('___')
            if(paramsMapGroupedByProcess.containsKey(keySplitted[0])) {
                paramsMapGroupedByProcess.get(keySplitted[0]).put(keySplitted[1], it.value)
            } else {
                tmp = [:]
                tmp.put(keySplitted[1], it.value)
                paramsMapGroupedByProcess.put(keySplitted[0], tmp)
            }
        }
    }

    // Global parameters (i.e.: not containing ___)
    // Populate the process parameters with the global parameters
    params.each {
        key = it.key
        value = it.value
        if(!key.contains('___')) {
            paramsMapGroupedByProcess.each {
                paramsMapGroupedByProcess.get(it.key).put(key, value)
            }
        }
    }

    return paramsMapGroupedByProcess
}

import static groovy.json.JsonOutput.*

def prettyPrint(object) {
    println(prettyPrint(toJson(object)))
}