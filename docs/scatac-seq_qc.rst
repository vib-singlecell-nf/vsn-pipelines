scATAC-seq QC and Cell Calling
==============================

This workflow uses the Python implementation of `cisTopic <https://github.com/aertslab/cisTopic>`_ (pycisTopic) to perform quality control and cell calling.
The inputs here are a fragments and bam file for each sample.

This workflow is currently available in the ``develop_atac`` branch (use ``nextflow pull vib-singlecell-nf/vsn-pipelines -r develop_atac`` to sync this branch).

----


Running the workflow
********************

Technical considerations
------------------------

1. Direct the Nextflow work directory to an alternate path (e.g. a scratch drive) using the ``NXF_WORK`` environmental variable::

    nwork=/path/to/scratch/example_project
    mkdir $nwork
    export NXF_WORK=$nwork

Note that if you start a new shell, ``NXF_WORK`` must be set again, or the pipeline will not resume properly.


2. Important for pycisTopic Ray issues: the system default temp location may become full.
   A workaround is to include a volume mapping to the alternate ``/tmp`` ``-B /alternate/path/to/tmp:/tmp`` using the volume mount options in Docker or Singularity.
   For example in the container engine options:
  - Singularity run options: ``runOptions = '--cleanenv -H $PWD -B /data,/alternate/path/to/tmp:/tmp'``
  - Docker run options: ``runOptions = '-i -v /data:/data -v /alternate/path/to/tmp:/tmp'``


3. Use the ``--quiet`` flag with ``nextflow run`` to suppress the printing of each file that is detected by the pipeline.

----

Configuration
-------------

For each sample, this pipeline take a bam and a fragments file.
These can be specified separately, or from a Cell Ranger ATAC/ARC ``outs/`` path.

Input with independent bam and fragments files
______________________________________________

Use the profiles ``bam`` and ``fragments``::

    nextflow config vib-singlecell-nf/vsn-pipelines/main_atac.nf \
        -profile atac_qc_filtering,bam,fragments,vsc,pycistopic_hg38 \
        > atac_qc.config

Preset profiles are available for human (``pycistopic_hg38``), mouse (``pycistopic_mm10``), and fly (``pycistopic_dmel``).
Or, these profiles can be omitted and set manually in the config (biomart, macs2).


Input data (bam and fragments files) are specified in the data section::

    data {
       fragments {
          file_paths = '/staging/leuven/stg_00002/lcb/cflerin/analysis/asap/20210527_hydrop-atac_asabr/atac_preprocess/out_run1/data/fragments/ASA__*tsv.gz'
          suffix = '.sinto.fragments.tsv.gz'
          index_extension = '.tbi'
       }
       bam {
          file_paths = '/staging/leuven/stg_00002/lcb/cflerin/analysis/asap/20210527_hydrop-atac_asabr/atac_preprocess/out_run1/data/bam/ASA*bam'
          suffix = '.bwa.out.possorted.bam'
          index_extension = '.bai'
       }
    }


Multiple files can be specified with ``*`` in ``file_paths`` or by separating the paths with a comma.

.. warning::

    The ``suffix`` for both bam and fragments will be removed from the filename to get sample IDs.
    The sample names obtained must match between bam and fragments for the files to be paired properly in the workflow.


Input with Cell Ranger ATAC data
________________________________

Use the ``tenx_atac`` profile::

    nextflow config vib-singlecell-nf/vsn-pipelines/main_atac.nf \
        -profile atac_qc_filtering,tenx_atac,vsc,pycistopic_hg38 \
        > atac_qc.config

Input data (the Cell Ranger ``outs/`` path) are specified in the data section::

    data {
        tenx_atac {
            cellranger_mex = '/data/cellranger_atac_2.0/*/outs,/data/processed/cellranger_arc_2.0.0/*/outs'
        }
    }

Multiple files can be specified with ``*`` in ``tenx_atac`` or by separating the paths with a comma.


Input directly from the preprocessing pipeline
______________________________________________

It is also possible to run these QC steps directly after the ``atac_preprocess`` pipeline, with a single command.
In this case, all the appropriate configuration profiles must be included at the configuration start::

    nextflow config vib-singlecell-nf/vsn-pipelines/main_atac.nf \
        -profile atac_preprocess,atac_qc_filtering,pycistopic_hg38,vsc \
        > atac_preprocess_and_qc.config

Note that here, we do not include ``bam`` and ``fragments`` profiles that specify the input data locations to the QC steps since these are piped directly from the preprocessing pipeline.
One caveat to this is that it could potentially make it harder to run the qc pipeline with ``-resume`` later on, especially if the Nextflow ``work/`` directory is not saved due to disk space concerns.

To execute the preprocessing and mapping pipeline in one step, use the ``atac_preprocess_with_qc`` entry point::

    nextflow -C atac_preprocess_and_qc.config run \
        vib-singlecell-nf/vsn-pipelines/main_atac.nf \
        -entry atac_preprocess_with_qc -resume --quiet


----

Execution
---------

After configuring, the workflow can be run with:

.. code:: bash

    nextflow -C atac_qc.config run \
        vib-singlecell-nf/vsn-pipelines/main_atac.nf \
        -entry atac_qc_filtering --quiet -resume

After completing, view the report in ``out/notebooks/<project_name>__pycisTopic_QC_report.html``. To change the filtering settings, use the ``params.tools.pycistopic.call_cells`` section.

Adjusting the filter settings
-----------------------------

In the pycisTopic parameters, filter settings can be applied in this section::

    pycistopic {
        call_cells {
            report_ipynb = '/src/pycistopic/bin/pycisTopic_qc_report_template.ipynb'
            use_density_coloring_on_scatterplot = true
            use_detailed_title_on_scatterplot = true
            filter_frags_lower = '1000'
            filter_frags_upper = ''
            filter_tss_lower = '8'
            filter_tss_upper = ''
            filter_frip_lower = ''
            filter_frip_upper = ''
            filter_dup_rate_lower = ''
            filter_dup_rate_upper = ''
        }
    }

If a setting is empty (``''``), this filter will not be applied.
If set to a single value (i.e. ``filter_frags_lower=1000``), this will apply this filter value to all samples.
To use sample-specific filters, this can be written as::

    filter_frags_lower = [
      'default': 1000,
      'Sample_1': 1500,
      'Sample_2': 2000,
    ]

The ``default`` setting (optional) is applied to all samples not listed in array.
If this default setting is missing, no filter will be applied to samples not listed in the array (all barcodes kept).

After setting the filters, the pipeline can be re-run to apply the new filters (use ``-resume``).

The additional settings control the output of the scatter plots in the report:
* ``use_density_coloring_on_scatterplot``: Slower when turned on; it can be helpful to set this to ``false`` until the proper thresholds are determined.
* ``use_detailed_title_on_scatterplot``: Adds the cell count and median values after filtering to the title of each plot.

----

Output
******

An example output tree is shown here.

.. code:: bash

    out/
    ├── data
    │   ├── macs2
    │   │   ├── sample_1.peaks.narrowPeak
    │   │   ├── sample_1.summits.bed
    │   │   ├── sample_2.peaks.narrowPeak
    │   │   └── sample_2.summits.bed
    │   └── pycistopic
    │       └── qc
    │           ├── benchmark_library_downsampled__metadata.pickle
    │           ├── benchmark_library_downsampled__profile_data.pickle
    │           ├── selected_barcodes
    │           │   ├── sample_1.cell_barcodes.txt
    │           │   └── sample_2.cell_barcodes.txt
    │           └── selected_barcodes_nFrag
    │               ├── sample_1.barcodes_nFrag_thr.txt
    │               └── sample_2.barcodes_nFrag_thr.txt
    └── notebooks
        ├── example_project__pycisTopic_QC_report.html
        └── example_project__pycisTopic_QC_report.ipynb


* ``macs2``: contains the narrowPeak and bed file for each sample.
* ``pycistopic``:
  * ``qc``: contains Python objects (in pickle format) for the metadata and profile data computed by pycisTopic.
    * ``selected_barcodes``: contains a text file with selected cell barcodes (one per line) based on the thresholds set in the config file.
    * ``selected_barcodes_nFrag``: contains a text file with barcodes (one per line) that have unique fragment counts greater than the ``params.tools.pycistopic.compute_qc_stats.n_frag`` setting in the pycisTopic parameters.

