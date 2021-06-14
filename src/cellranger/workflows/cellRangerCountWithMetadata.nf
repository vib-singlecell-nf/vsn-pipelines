nextflow.enable.dsl=2

import java.nio.file.Paths

EMPTY_VALUES = ["n/a","n.a.","none","null"]

//////////////////////////////////////////////////////
//  process imports:

include {
    SC__CELLRANGER__COUNT_WITH_METADATA;
} from './../processes/count' params(params)
include {
    SC__CELLRANGER__PREFLIGHT;
} from './../processes/preflight' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow CELLRANGER_COUNT_WITH_METADATA {

    def getFastQsFilePaths = { fastqsParentDirPath, fastqsDirName ->
        if(fastqsParentDirPath instanceof List)
            return fastqsParentDirPath
        if(fastqsDirName in EMPTY_VALUES)
            if(!fastqsParentDirPath.contains(','))
                return file(fastqsParentDirPath)
            return Arrays.asList(fastqsParentDirPath.split(',')).collect { file(it) }
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
                row.containsKey("expect_cells") ? row.expect_cells : 'NULL',
                row.containsKey("chemistry") ? row.chemistry : 'NULL'
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
                *it[4..(it.size()-1)].collect { 
                    arg -> 
                    argBySample = arg.unique()
                    if(argBySample.size() != 1)
                        throw new Exception("Multiple values of a Cell Ranger parameter detected for sample: "+ it[0])
                    return argBySample[0]
                }
            )
        }
        
        SC__CELLRANGER__PREFLIGHT()
        SC__CELLRANGER__COUNT_WITH_METADATA( transcriptome, data )

    emit:
        SC__CELLRANGER__COUNT_WITH_METADATA.out

}
