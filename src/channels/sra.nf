nextflow.preview.dsl=2

workflow getChannel {

    take:
        // Expects sra Map [[id: "id1", samples: ["glob1", ...]], ...]
        sra

    main:
        channel = Channel.fromList(
            sra
        ).map {
            it -> tuple(it.id, it.samples)
        }.view()

    emit:
        channel

}