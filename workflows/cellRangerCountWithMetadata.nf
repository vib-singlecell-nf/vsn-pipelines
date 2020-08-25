nextflow.preview.dsl=2

import java.nio.file.Paths

//////////////////////////////////////////////////////
//  process imports:

include {
    SC__CELLRANGER__COUNT_WITH_METADATA;
} from './../processes/count' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow CELLRANGER_COUNT_WITH_METADATA {

    def getFastQsFilePaths = { fastqsParentDirPath, fastqsDirName ->
        if(fastqsDirName == "n/a" || fastqsDirName == "n.a." || fastqsDirName == "none" || fastqsDirName == "null")
            if(!fastqsParentDirPath.contains(','))
                return file(fastqsParentDirPath)
            return Arrays.asList(glob.split(',')).collect { file(it) }
        return file(Paths.get(fastqsParentDirPath, fastqsDirName))
    }

    take:
        transcriptome
        metadata

    main:
        // Define the sampleId
        data = Channel.from(
            metadata
        ).splitCsv(
            header:true,
            sep: '\t'
        ).map {
            row -> tuple(
                row.containsKey("short_uuid") ? row.short_uuid + "__" + row.sample_name : row.sample_name,
                row.fastqs_sample_prefix,
                row.fastqs_parent_dir_path,
                row.fastqs_dir_name,
                // Begin CellRanger parameters
                row.expect_cells
            )
        }.groupTuple(
            by: 0
        ).map {
            it -> tuple(
                it[0],
                *it[1].unique(), 
                getFastQsFilePaths(
                    it[2],
                    it[3]
                ),
                *it[4].unique())
        }
        SC__CELLRANGER__COUNT_WITH_METADATA( transcriptome, data )

    emit:
        SC__CELLRANGER__COUNT_WITH_METADATA.out

}
