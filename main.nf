//
// Version:
// Test:
// Command: 
// 
//
/*
 * Remote run test
 * Source:
 * 
 * Steps considered: 

 */ 
import static groovy.json.JsonOutput.*

nextflow.preview.dsl=2


params.each { println "${it}" }
println "${workflow.projectDir}"

//////////////////////////////////////////////////////
//  Define the parameters for current testing proces

include SCENIC from './src/scenic/main.nf' params(params)

include CELLRANGER from './src/cellranger/main.nf' params(params)


workflow {
    // SCENIC( file( params.sc.scenic.filteredloom ) )
    CELLRANGER()
}










