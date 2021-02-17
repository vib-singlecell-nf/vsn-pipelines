Input Data Formats
===================

Depending on the type of data you run the pipeline with, one or more appropriate profiles should be set when running ``nextflow config``.
These profiles are indicated in the sections below.

Specifying multiple samples
***************************

All the input data parameters are compatible with the following features:

- Glob patterns

.. code::

    "data/10x/1k_pbmc/1k_pbmc_*/outs/"

- Comma separated paths (paths can contain glob patterns)

.. code::

    "data/10x/1k_pbmc/1k_pbmc_v2_chemistry/outs/, data/10x/1k_pbmc/1k_pbmc_v3_chemistry/outs/"

- Array of paths (paths can contain glob patterns)

.. code::

    [
        "data/10x/1k_pbmc/1k_pbmc_v2_chemistry/outs/",
        "data/10x/1k_pbmc/1k_pbmc_v3_chemistry/outs/"
    ]

----

.. _using_10x_datasets:

Cell Ranger (10x Genomics)
**************************

Data from a standard Cell Ranger output directory can be easily ingested into the pipeline by using the proper input channel (``tenx_mex`` or ``tenx_h5``, depending on which file should be used).
Multiple samples can be selected by providing the path to this directory using glob patterns.

.. code::

    /home/data/
    └── cellranger
        ├── sample_A
        │   └── outs
        │       ├── filtered_feature_bc_matrix
        │       │   ├── barcodes.tsv
        │       │   ├── genes.tsv
        │       │   └── matrix.mtx
        │       └── filtered_feature_bc_matrix.h5
        └── sample_B
            └── outs
                ├── filtered_feature_bc_matrix
                │   ├── barcodes.tsv
                │   ├── genes.tsv
                │   └── matrix.mtx
                └── filtered_feature_bc_matrix.h5


MEX
___

To use the Cell Ranger Market Exchange (**MEX**) files, use the following profile when generating the config file::

    -profile tenx

This profile adds the following parameter (``params.data.tenx.cellranger_mex``) into the generated .config file::

    [...]
    data {
        tenx {
            cellranger_mex = "/home/data/cellranger/sample*/outs/"
        }
    }
    [...]


H5
__

To use the Cell Ranger ``h5`` file as input, use the following profile::

    -profile tenx_h5

This profile adds the ``params.data.tenx.cellranger_h5`` parameter into the generated .config file::

    [...]
    data {
        tenx {
            cellranger_h5 = "/home/data/cellranger/sample*/outs/"
        }
    }
    [...]


Input file detection
____________________

Setting the input directory appropriately, using a glob in the directory path in place of the sample names, will collect all the samples listed in the ``filtered_[feature|gene]_bc_matrix`` directories listed above.
For example, in ``params.data.tenx``, setting::

    cellranger_mex = "/home/data/cellranger/sample*/outs/"

or

.. code::

    cellranger_h5 = "/home/data/cellranger/sample*/outs/"

will recursively find all 10x samples in that directory.

The pipeline will use either the ``outs/filtered_feature_bc_matrix/`` or the ``outs/raw_feature_bc_matrix/`` depending on the setting of the ``params.utils.file_converter.useFilteredMatrix`` (``true`` uses filtered; ``false`` uses raw).

----

H5AD (Scanpy)
*************
Use the following profile when generating the config file::

    -profile h5ad


In the generated .config file, make sure the ``file_paths`` parameter is set with the paths to the ``.h5ad`` files::

    [...]
    data {
        h5ad {
            file_paths = "data/1k_pbmc_v*_chemistry_SUFFIX.SC__FILE_CONVERTER.h5ad"
            suffix = "_SUFFIX.SC__FILE_CONVERTER.h5ad"
        }
    }
    [...]

- The ``suffix`` parameter is used to infer the sample name from the file paths (it is removed from the input file path to derive a sample name).

In case there are multiple .h5ad files that need to be processed with different suffixes, the multi-labelled strategy should be used to define the h5ad parameter::

    [...]
    data {
        h5ad {
            GROUP1 {
                file_paths = "[path-to-group1-files]/*.SUFFIX1.h5ad"
                suffix = ".SUFFIX1.h5ad"
                group = ["technology", "10x"]
            }
            GROUP2 {
                file_paths = "[path-to-group1-files]/*.SUFFIX2.h5ad"
                suffix = ".SUFFIX2.h5ad"
                group = ["technology", "smart-seq2"]
            }
        }
    }
    [...]

Notes: 

- ``GROUP1``, ``GROUP2`` are just example names here. They can be replaced by any value as long as they are alphanumeric (underscores are allowed).
- All the different `suffix` defined should unique.
- ``file_paths`` and ``suffix`` do allow list of paths/globs in the multi-labelled strategy.
- ``group`` [optional] should be an array of 2 elements where first element define the group name and the second the group value. This will add cell-based annotation for each group of files

----

Loom
****

Use the following profile when generating the config file::

    -profile loom


In the generated .config file, make sure the ``file_paths`` parameter is set with the paths to the ``.loom`` files::

    [...]
    data {
        loom {
            file_paths = "data/1k_pbmc_v*_chemistry_SUFFIX.SC__FILE_CONVERTER.loom"
            suffix = "_SUFFIX.SC__FILE_CONVERTER.loom"
        }
    }
    [...]

- The ``suffix`` parameter is used to infer the sample name from the file paths (it is removed from the input file path to derive a sample name).

----

Seurat Rds
**********

Use the following profile when generating the config file::

    -profile seurat_rds


In the generated .config file, make sure the ``file_paths`` parameter is set with the paths to the ``.Rds`` files::

    [...]
    data {
        seurat_rds {
            file_paths = "data/1k_pbmc_v*_chemistry_SUFFIX.SC__FILE_CONVERTER.Rds"
            suffix = "_SUFFIX.SC__FILE_CONVERTER.Rds"
        }
    }
    [...]

- The pipelines expect a Seurat v3 object contained in the .Rds file. (Seurat v2 objects are currently not supported).
- The ``suffix`` parameter is used to infer the sample name from the file paths (it is removed from the input file path to derive a sample name).

----

TSV
***
Use the following profile when generating the config file::

    -profile tsv


In the generated .config file, make sure the ``file_paths`` parameter is set with the paths to the ``.tsv`` files::

    [...]
    data {
        h5ad {
            file_paths = "data/1k_pbmc_v*_chemistry_SUFFIX.SC__FILE_CONVERTER.tsv"
            suffix = "_SUFFIX.SC__FILE_CONVERTER.tsv"
        }
    }
    [...]

- The ``suffix`` parameter is used to infer the sample name from the file paths (it is removed from the input file path to derive a sample name).

----

CSV
***
Use the following profile when generating the config file::

    -profile csv


In the generated .config file, make sure the ``file_paths`` parameter is set with the paths to the ``.csv`` files::

    [...]
    data {
        h5ad {
            file_paths = "data/1k_pbmc_v*_chemistry_SUFFIX.SC__FILE_CONVERTER.csv"
            suffix = "_SUFFIX.SC__FILE_CONVERTER.csv"
        }
    }
    [...]

- The ``suffix`` parameter is used to infer the sample name from the file paths (it is removed from the input file path to derive a sample name).


