nextflow.preview.dsl=2

//////////////////////////////////////////////////////
//  Import sub-workflows from the modules:

include SC__FILE_CONVERTER from '../utils/processes/utils.nf' params(params)

include SC__TEMPLATE__PROCESS1 from './processes/process1.nf' params(params)


//////////////////////////////////////////////////////
// Define the workflow

workflow template {

    take:
        data

    main:
        data = SC__FILE_CONVERTER(data)
        data.view()

        SC__TEMPLATE__PROCESS1(data)

}

