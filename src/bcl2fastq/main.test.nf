nextflow.enable.dsl=2

include {
    getDataChannel;
} from './../channels/channels.nf' params(params)



workflow {

    main:
        switch(params.test) {
            case "BCL2FASTQ_DEMULTIPLEX":
                include {
                    BCL2FASTQ__DEMULTIPLEX;
                } from './processes/demultiplex'
            break;
            case "getdata":
                getDataChannel | view
        }
}