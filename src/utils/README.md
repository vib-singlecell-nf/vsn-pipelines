# Utils Module

## Cell-based Metadata Annotation

The profile `utils_cell_annotate` should be added when generating the main config using `nextflow config`. This will add the following entry in the config:

```
params {
    sc {
        cell_annotate {
            iff = '10x_mtx'
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

- `cellMetaDataFilePath` is a TSV file (with header) with at least 2 columns: a column containing all the cell IDs and an annotation column.
- `indexColumnName` is the column name from `cellMetaDataFilePath` containing the cell IDs information.
- `sampleColumnName` is the column name from `cellMetaDataFilePath` containing the sample ID/name information.
- `annotationColumnNames` is an array of columns names from `cellMetaDataFilePath` containing different annotation metadata to add.
