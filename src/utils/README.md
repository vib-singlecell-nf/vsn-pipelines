# Utils module

## Cell-based metadata annotation

The profile `utils_cell_annotate` should be added when generating the main config using `nextflow config`. This will add the following entry in the config:

```
params {
    tools {
        cell_annotate {
            iff = '10x_cellranger_mex'
            off = 'h5ad'
            cellMetaDataFilePath = ''
            indexColumnName = ''
            sampleColumnName = ''
            annotationColumnNames = ['']
        }
    }
}
```
Then, the following parameters should be updated to use the module feature:

- `cellMetaDataFilePath` is a .tsv file (with header) with at least 2 columns: a column containing all the cell IDs and an annotation column.
- `indexColumnName` is the column name from `cellMetaDataFilePath` containing the cell IDs information.
- `sampleColumnName` is the column name from `cellMetaDataFilePath` containing the sample ID/name information.
- `annotationColumnNames` is an array of columns names from `cellMetaDataFilePath` containing different annotation metadata to add.

## Sample-based metadata annotation
The profile `utils_sample_annotate` should be added when generating the main config using nextflow config. This will add the following entry in the config:

```
params {
    tools {
        sample_annotate {
            iff = '10x_cellranger_mex'
            off = 'h5ad' 
            type = 'sample' 
            metadataFilePath = 'data/10x/1k_pbmc/metadata.tsv'
        }
    }
}
```
Then, the following parameters should be updated to use the module feature:

- `metadataFilePath` is a .tsv file (with header) with at least 2 columns where the first column need to match the sample IDs. Any other columns will be added as annotation in the final loom i.e.: all the cells related to their sample will get annotated with their given annotations.

| id  | chemistry | ... |
| ------------- | ------------- | ------------- |
| 1k_pbmc_v2_chemistry  | v2  | ... |
| 1k_pbmc_v3_chemistry  | v3  | ... |

Sample-annotating the samples using this system will allow any user to query all the annotation using the SCope portal. This is especially relevant when samples needs to be compared across specific annotations (check compare tab with SCope).
