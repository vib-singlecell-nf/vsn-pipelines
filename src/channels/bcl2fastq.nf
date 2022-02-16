nextflow.enable.dsl=2


workflow getChannel {
    take:
        runFolder
        sampleSheet
    
    main:
        runFolder_channel = channel
            .fromPath(runFolder, type: 'dir', checkIfExists: true)
        sampleSheet_channel = channel
            .fromPath(sampleSheet, checkIfExists: true)
        data_channel = runFolder_channel | combine(sampleSheet_channel)
        
    emit:
        data_channel
}