nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include SC__TEMPLATE__PROCESS1 from './processes/process1.nf' params(params)


//////////////////////////////////////////////////////
// Define the workflow

workflow template {

    take:
        data

    main:
        SC__TEMPLATE__PROCESS1(data)

}
