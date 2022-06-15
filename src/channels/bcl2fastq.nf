nextflow.enable.dsl=2


workflow getChannel {
    take:
        runFolder
        sampleSheet
    
    main:
        runFolder_channel = Channel
            .fromPath(runFolder, type: 'dir', checkIfExists: true)
        sampleSheet_channel = Channel
            .fromPath(sampleSheet, checkIfExists: true)
        data_channel = runFolder_channel
                       .combine(sampleSheet_channel)
        
    emit:
        data_channel
}