nextflow.enable.dsl=2

include {
    SC__DIRECTS__SELECT_DEFAULT_CLUSTERING
} from './processes/selectDefaultClustering.nf'

// Test 1: SC__DIRECTS__SELECT_DEFAULT_CLUSTERING (from processes/)
// Time: ?
// Command: 
//  cd tests
//  nextflow config .. -profile test__select_default_clustering > test__select_default_clustering.config
//  nextflow -C test__select_default_clustering.config run ../main.test.nf --test SC__DIRECTS__SELECT_DEFAULT_CLUSTERING

workflow {

    main:
        switch(params.test) {
            case "SC__DIRECTS__SELECT_DEFAULT_CLUSTERING":
                test = Channel.of(tuple('TEST', params.tools.directs.inputLoom, null))
                SC__DIRECTS__SELECT_DEFAULT_CLUSTERING( test )
            break;
            default:
                throw new Exception("The test parameters should be specified.")
            break;
        }

}
