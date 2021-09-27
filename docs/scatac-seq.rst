scATAC-seq Preprocessing
========================


This pipeline takes fastq files from paired end single cell ATAC-seq, and applies preprocessing steps to align the reads to a reference genome, and produce a bam file and scATAC-seq fragments file.

This workflow is currently available in the ``develop_atac`` branch (use ``nextflow pull vib-singlecell-nf/vsn-pipelines -r develop_atac`` to sync this branch).

----

Pipeline Steps
**************

The full steps are:

- Barcode correction:

  * For 'standard' and 'multiome' samples (e.g. 10x Genomics or similar) correction is performed against a whitelist by 
    `this method <https://github.com/aertslab/single_cell_toolkit/blob/master/correct_barcode_in_fastq.sh>`_ 
    from `aertslab/single_cell_toolkit <https://github.com/aertslab/single_cell_toolkit>`_.
  * For 'biorad' samples, barcode correction is performed by
    `this script <https://github.com/aertslab/single_cell_toolkit/blob/master/extract_and_correct_biorad_barcode_in_fastq.sh>`_
    in our `aertslab/single_cell_toolkit <https://github.com/aertslab/single_cell_toolkit>`_
    (previously, this was done with `BAP <https://github.com/caleblareau/bap>`_).

- Fastq barcoding: Add the barcode sequence to the comment field of the fastq sequence identifier.
  Uses methods from `aertslab/single_cell_toolkit <https://github.com/aertslab/single_cell_toolkit>`_.
- Read/adapter trimming 
  (`Trim_Galore <https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/>`_
  or `fastp <https://github.com/OpenGene/fastp>`_).
- Mapping to a reference genome:

  * ``bwa mem`` is used with default parameters, with a choice of the original 
    `bwa mem <https://github.com/lh3/bwa>`_, or `bwa-mem2 <https://github.com/bwa-mem2/bwa-mem2>`_.
- Mark PCR and optical duplicates (`MarkDuplicates (Picard) <https://gatk.broadinstitute.org/hc/en-us/articles/360057439771-MarkDuplicates-Picard->`_ 
  or `MarkDuplicatesSpark (GATK) <https://gatk.broadinstitute.org/hc/en-us/articles/360057438771-MarkDuplicatesSpark>`_).
- Estimate library complexity with 
  `EstimateLibraryComplexity (Picard) <https://gatk.broadinstitute.org/hc/en-us/articles/360057438451-EstimateLibraryComplexity-Picard->`_.
- A fragments file is created using `Sinto <https://github.com/timoast/sinto>`_.

----

Pipeline Details
****************

Input Metadata
--------------

The input to this pipeline is a (tab-delimited) metadata table with the sample ID, sequencing technology, and locations of the fastq files.
Note that the fastq file fields must be full paths; this is not shown here for clarity:

.. list-table:: Metadata Table
    :widths: 10 10 10 10 10
    :header-rows: 1

    * - sample_name
      - technology
      - fastq_PE1_path
      - fastq_barcode_path
      - fastq_PE2_path
    * - sample_1
      - standard
      - sample_1_R1.fastq.gz
      - sample_1_R2.fastq.gz
      - sample_1_R3.fastq.gz
    * - sample_2
      - multiome
      - sample_2_R1.fastq.gz
      - sample_2_R2.fastq.gz
      - sample_2_R3.fastq.gz
    * - sample_3
      - biorad
      - sample_3_R1.fastq.gz
      -  
      - sample_3_R3.fastq.gz
    * - sample_4
      - revcomp_wl
      - sample_4_R1.fastq.gz
      - sample_2_R2.fastq.gz
      - sample_4_R3.fastq.gz
    * - sample_5
      - hydrop_3x96
      - sample_5_R1.fastq.gz
      - sample_5_R2.fastq.gz
      - sample_5_R3.fastq.gz
    * - sample_6
      - hydrop_2x384
      - sample_5_R1.fastq.gz
      - sample_5_R2.fastq.gz
      - sample_5_R3.fastq.gz

The columns represent:

- ``sample_name`` Sample name for labeling the sample in the pipeline and output files. This can be any arbitrary string.
- ``technology``: This controls the barcode correction and processing methods to use for the fastq files. Currently only the ``biorad`` option involves different processing steps. Otherwise, the value in this field (e.g. ``standard``, ``multiome``) controls which barcode whitelist is used for correction. See below for additional details.
- ``fastq_PE1_path``: The full path to the fastq file for the first read in a pair.
- ``fastq_barcode_path``: The full path to the fastq file containing the barcodes. This column can be blank/empty depending on the technology setting (e.g. ``biorad``).
- ``fastq_PE2_path``: The full path to the fastq file for the second read in a pair.


Fastq input
-----------

Fastq input for each sample can be given as a single set of files (R1, R2, R3), or it can be multiple files, in the case of samples which have been split across multiple sequencing lanes.
A combination of these cases can be processed together in one metadata file.

Within the pipeline, all of the reads from each set of files will be considered and labeled as one read group.
Read group names are taken from the first read in the fastq::

    @A01044:19:HLYKFDRXX:1:2101:4291:1000 1:N:0:ACTCAGAC

will produce a RG in the bam::

    @RG     ID:A01044:19:HLYKFDRXX:1        SM:sample_1        LB:A01044:19:HLYKFDRXX:1__sample_1 PL:ILLUMINA


Single fastq input
__________________

In this situation, there is only one set of reads per sample, the metadata file will look very similar to one of the rows from above.
There will be one read group in the final bam file.


Split fastq input
_________________

In this case, multiple fastq files (in rows) for each sample can be given.
In this example, there are two sets of fastqs for ``sample_1`` that were run on two separate lanes.
Note that the sample ID is the same for both rows:

.. list-table:: Metadata Table with split fastqs
    :widths: 10 10 10 10 10
    :header-rows: 1

    * - sample_name
      - technology
      - fastq_PE1_path
      - fastq_barcode_path
      - fastq_PE2_path
    * - sample_1
      - standard
      - sample_1_S1_L001_R1_001.fastq.gz
      - sample_1_S1_L001_R2_001.fastq.gz
      - sample_1_S1_L001_R3_001.fastq.gz
    * - sample_1
      - standard
      - sample_1_S1_L002_R1_001.fastq.gz
      - sample_1_S1_L002_R2_001.fastq.gz
      - sample_1_S1_L002_R3_001.fastq.gz

In this situation, each set of fastqs will be processed separately for the barcode correction, barcode addition to the fastq comment field, adaptor trimming, and mapping steps.
Following mapping, each mapped bam is merged and duplicates are marked using the full data.
Downstream steps are done with the merged data.


Generating the metadata file
----------------------------

Note that there is an easy way to create the metadata from the file paths for each sample by using the following bash command (expand to view).
Special thanks here to Gert Hulselmans for expanding the capabilities of this function.

.. raw:: html

   <details>
   <summary><a>metadata generator</a></summary>

.. code-block:: none

    create_atac_metadata () {
        local sample="${1}";
        local technology="${2}";
        local fastq_prefix="${3}";
        local read_labels="${4}";
        if [ "${sample}" == "header" ]; then
            printf 'sample_name\ttechnology\tfastq_PE1_path\tfastq_barcode_path\tfastq_PE2_path\n';
            return 0;
        fi
        if [ ${#@} -ne 4 ] ; then
            printf 'Usage: create_atac_metadata sample technology fastq_prefix read_labels\n\n';
            printf 'Arguments:\n';
            printf '    sample:       sample name\n';
            printf '    technology:   "standard", "hydrop_3x96", "hydrop_2x384", or "biorad"\n';
            printf '    fastq_prefix: path prefix to FASTQ files.\n';
            printf '    read_labels:  comma separated read labels for R1, R2 and R3 that select: R1,R2,R3.\n';
            return 1;
        fi
        read_labels_arr=(${read_labels//,/ });
        # Get R1, R2 and R3 FASTQ filenames for
        R1=(${fastq_prefix}*${read_labels_arr[0]}*.{fastq,fq,fastq.gz,fq.gz})
        R2=(${fastq_prefix}*${read_labels_arr[1]}*.{fastq,fq,fastq.gz,fq.gz})
        R3=(${fastq_prefix}*${read_labels_arr[2]}*.{fastq,fq,fastq.gz,fq.gz})
        for i in "${!R1[@]}" ; do
            # Check if R1 FASTQ file exist (and is not just a glob like "${sample}*R1*.fq").
            if [ -e "${R1[i]}" ] ; then
                printf "${sample}\t${technology}\t${R1[i]}\t${R2[i]}\t${R3[i]}\n";
            fi
        done
    }

To run use the options:

#. Sample ID (if this parameter is "header", it will print the metadata header and stop)
#. Technology (e.g. "standard")
#. The "file prefix" full path to your fastq files, matching the common portions of the file names (without any glob ``*`` expansions)
#. The "read labels" to indicate how the files are named and match the remainder of the file names (e.g. "R1,R2,R3", "R1,UMI,R2", etc.)

.. code-block:: none

    create_atac_metadata header > auto_metadata.tsv
    create_atac_metadata sample_1 standard /path/to/sample_1_subset_S R1,R2,R3 >> auto_metadata.tsv
    create_atac_metadata sample_2 standard /path/to/sample_2_subset_S R1,R2,R3 >> auto_metadata.tsv
    create_atac_metadata sample_5 hydrop_3x96 /path/to/sample_5_ R1,R2,R3 >> auto_metadata.tsv
    create_atac_metadata sample_6 hydrop_2x384 /path/to/sample_6_ R1,R2,R3 >> auto_metadata.tsv

.. raw:: html

   </details>

----

Technology types
----------------

The "technology" field in the metadata table controls two things:

1. **How technology-specific pipeline steps are applied.**
   Currently there are three specific settings (``biorad``, ``hydrop_3x96``, and ``hydrop_2x384``) that use alternate pipelines processes (to extract and correct the barcode sequence from the input fastqs).
   Using any other keyword is allowed, and samples will be run with the standard pipeline steps (barcode correction against a whitelist).

2. **Which whitelist is used for barcode correction.**
   The "technology" field must match a key in the ``params.tools.singlecelltoolkit.barcode_correction.whitelist`` parameter list in the config file for that sample to be associated with a particular barcode whitelist.
   The "technology" field and whitelist key name can be set to any arbitrary string (e.g. ``standard``), with the exception of the technology-specific keywords above.

The main modes are:

``standard`` 
____________

The ``standard`` setting is the main pipeline mode.
It assumes a typical 10x Genomics style format with two read pair fastqs and a barcode fastq (note that in the example here, the barcode correction has already been performed, writing the ``CB`` tag into the comment of the barcode fastq)::

    $ zcat sample_1_R1.fastq.gz | head -n 4
    @A00311:74:HMLK5DMXX:1:1101:2013:1000 1:N:0:ACTCAGAC
    NTTGTCTCAGCACCCCCCGACATGGATTCAGGCTGTCTCTTATACACATC
    +
    #FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

    $ zcat sample_1_R2.fastq.gz | head -n 4
    @A00311:74:HMLK5DMXX:1:1101:2013:1000 2:N:0:ACTCAGAC CB:Z:CTGTTCGCAAAGCATA
    CTGTTCGCAAAGCATA
    +
    F:FFFFFFFFFFFFFF

    $ zcat sample_1_R3.fastq.gz | head -n 4
    @A00311:74:HMLK5DMXX:1:1101:2013:1000 3:N:0:ACTCAGAC
    CCTGAATCCATGTCGGGGGGTGCTGAGACAAGCTGTCTCTTATACACAT
    +
    FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

The barcoding step here uses a 
`helper script <https://github.com/aertslab/single_cell_toolkit/blob/master/barcode_10x_scatac_fastqs.sh>`_
from `aertslab/single_cell_toolkit <https://github.com/aertslab/single_cell_toolkit>`_
which transforms this input into two paired fastq files with the barcode information embedded in the fastq comments field::

    $ zcat sample_1_dex_R1.fastq.gz | head -n 4
    @A00311:74:HMLK5DMXX:1:1101:2013:1000 CR:Z:CTGTTCGCAAAGCATA     CY:Z:F:FFFFFFFFFFFFFF   CB:Z:CTGTTCGCAAAGCATA
    NTTGTCTCAGCACCCCCCGACATGGATTCAGGCTGTCTCTTATACACATC
    +
    #FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

    $ zcat sample_1_dex_R2.fastq.gz | head -n 4
    @A00311:74:HMLK5DMXX:1:1101:2013:1000 CR:Z:CTGTTCGCAAAGCATA     CY:Z:F:FFFFFFFFFFFFFF   CB:Z:CTGTTCGCAAAGCATA
    CCTGAATCCATGTCGGGGGGTGCTGAGACAAGCTGTCTCTTATACACAT
    +
    FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF


``multiome``/alternate
______________________

The ``multiome`` or alternately-named settings work with the same pipeline steps as ``standard`` with the exception of the whitelist used for barcode correction.
The whitelists are supplied in the params file (``params.tools.singlecelltoolkit.barcode_correction.whitelist``).
This can be used to supply alternate whitelists for certain samples, for example if you need to supply a reverse complemented whitelist for samples run in certain sequencing machines.


``hydrop_3x96``/``hydrop_2x384``
__________

The HyDrop settings (either ``hydrop_3x96`` or ``hydrop_2x384`` depending on the library preparation used) processes data generated by the HyDrop ATAC protocol
(see `hydrop.aertslab.org <https://hydrop.aertslab.org/>`_ and `the associated preprint <https://doi.org/10.1101/2021.06.04.447104>`_).
This approach differs from the standard pipeline in only the initial step, which is to extract and process the HyDrop barcodes from the sequencing output.
Here, `this script <https://github.com/aertslab/single_cell_toolkit/blob/master/extract_hydrop_atac_barcode_from_R2_fastq.sh>`_ is used to take the R2 read from the sequencer::

    $ zcat sample_5_R2.fastq.gz | head -n 4
    @VH00445:5:AAAL5KYM5:1:1101:63923:1019 2:N:0:ACACGTGGAC
    CACTGGTGGTAGGGTACTCGGACAAGTGGAGCAGTAGCTGAAGTGTAGAAG
    +
    CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

and transform it into::

    $ zcat sample_5_hydrop_barcode_R2.fastq.gz
    @VH00445:5:AAAL5KYM5:1:1101:63923:1019 2:N:0:ACACGTGGAC
    CACTGGTGGTGACAAGTGGAAAGTGTAGAA
    +
    CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

The two HyDrop modes (``hydrop_3x96``, ``hydrop_2x384``) differ only in the way the initial barcode extraction script works.
Following this, they are processed in the same way as the standard pipeline, including whitelist-based barcode correction (note that the two HyDrop modes require different barcode whitelists to be used here).

``biorad`` 
__________

The ``biorad`` setting processes BioRad data using 
`this script <https://github.com/aertslab/single_cell_toolkit/blob/master/extract_and_correct_biorad_barcode_in_fastq.sh>`_
in our `aertslab/single_cell_toolkit <https://github.com/aertslab/single_cell_toolkit>`_
(previously, this was done with `BAP <https://github.com/caleblareau/bap>`_).
This takes input data::

    $ zcat sample_2_R1.fastq.gz | head -n 4
    @A00794:327:HTJ55DRXX:1:2101:1154:1016 1:N:0:TAAGGCGA
    GATCACCATATGCATGACATTCACGAGTCACTGAGTAACGCCTCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCTGCAATGGCTGGAGCACACCCCATACTCATTCTGGTCTCCTT
    +
    FFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFF:FFFFFFFFF:FFFFFFFFFFFFFFFFF,FF,FFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFF:F:FFFF,FFFFFF

    $ zcat sample_2_R2.fastq.gz | head -n 4
    @A00794:327:HTJ55DRXX:1:2101:1154:1016 2:N:0:TAAGGCGA
    GTGTTTGGCTGAGGAAAGTGTGTGAAGCCCCGATATGTGA
    +
    FFF,FFF:FFF:FF,FFFFF:F:FFFFFFFFFFF,,F:FF

and directly produces paired fastq files with the barcode added in the fastq comments field::

    $ zcat sample_2_dex_R1.fastq.gz | head -n 4
    @A00794:327:HTJ55DRXX:1:2101:1154:1016 CR:Z:GATCACCATTCACGTAACGCC       CY:Z:FFFFFFFFFFFFFF:FFFFFF      CB:Z:GATCACCATTCACGTAACGCC      br:Z:0,0,0_0,0,0,1
    CTGCAATGGCTGGAGCACACCCCATACTCATTCTGGTCTCCTT
    +
    F:FFFFFFFFFFFFFFFFFFFFFFFFFFF:F:FFFF,FFFFFF

    $ zcat sample_2_dex_R2.fastq.gz | head -n 4
    @A00794:327:HTJ55DRXX:1:2101:1154:1016 CR:Z:GATCACCATTCACGTAACGCC       CY:Z:FFFFFFFFFFFFFF:FFFFFF      CB:Z:GATCACCATTCACGTAACGCC      br:Z:0,0,0_0,0,0,1
    GTGTTTGGCTGAGGAAAGTGTGTGAAGCCCCGATATGTGA
    +
    FFF,FFF:FFF:FF,FFFFF:F:FFFFFFFFFFF,,F:FF

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


2. Temporary directory mapping.
   For large BAM files, the system default temp location may become full.
   A workaround is to include a volume mapping to the alternate ``/tmp`` ``-B /alternate/path/to/tmp:/tmp`` using the volume mount options in Docker or Singularity.
   For example in the container engine options:
  - Singularity run options: ``runOptions = '--cleanenv -H $PWD -B /data,/alternate/path/to/tmp:/tmp'``
  - Docker run options: ``runOptions = '-i -v /data:/data -v /alternate/path/to/tmp:/tmp'``


Configuration
-------------

To generate a config file, use the ``atac_preprocess`` profile along with ``docker`` or ``singularity``.
Note that the full path to ``vib-singlecell-nf/vsn-pipelines/main_atac.nf`` must be used:

.. code:: bash

    nextflow config \
        vib-singlecell-nf/vsn-pipelines/main_atac.nf \
        -profile atac_preprocess,singularity \
        > atac_preprocess.config


.. note::

    It is also possible to run the pycisTopic QC steps directly after this ``atac_preprocess`` pipeline, with a single command.
    Please see
    `here <scatac-seq_qc.html#input-directly-from-the-preprocessing-pipeline>`_
    for details on how to run with this configuration.


Parameters
----------

The ATAC-specific parameters are described here.
The important parameters to verify are:

- ``params.data.atac_preprocess.metadata``: the path to the metadata file.
- ``params.tools.bwamaptools.bwa_fasta``: the path to the bwa reference fasta file.
  This should be already indexed with ``bwa index``, and the index files located in the same directory as the fasta file. Note that ``bwa`` and ``bwa-mem2`` use different indexes that are not interchangeable.
- ``params.tools.singlecelltoolkit.barcode_correction.whitelist``: Whitelists for barcode correction are supplied here.
  The whitelists are matched to samples based on the parameter key here ('standard', 'multiome', 'hydrop_3x96', 'hydrop_2x384', etc.) and the technology field listed for each sample in the metadata file.
  Barcode whitelists can (optionally) be gzipped.
  There are currently no checks performed to ensure that the sample barcodes have any overlap to the whitelist (the barcode correction reports should be checked for this).


Choice of tools
_______________

Several steps have options for the choice of method to use.
These options are controlled within ``params.atac_preprocess_tools``.

- Adapter trimming (``adapter_trimming_method``): Can be either of ``Trim_Galore`` (default), or ``fastp``.
- Duplicate marking (``mark_duplicates_method``): Can be either of ``MarkDuplicates`` (Picard tools, default) or ``MarkDuplicatesSpark`` (GATK).
  We currently recommend Picard MarkDuplicates because it has the capability to perform barcode-aware marking of PCR duplicates.
  MarkDuplicatesSpark has the advantage of parallelization, however it requires a large SSD to use for temporary files.

Additionally:

- Mapping: Use parameter ``params.tools.bwamaptools.bwa_version`` to select either ``bwa`` or ``bwa-mem2``. These should give virtually identical results, however ``bwa-mem2``, while faster, has used more memory in our tests. Note that the index (``bwa_index``) is not interchangeable between the versions.


Optional parameters
___________________

- Within ``params.tools.sinto.fragments``:

  - One of (but not both) ``barcodetag`` or ``barcode_regex`` needs to be set to tell Sinto where to find the barcodes in the bam file. The default is to use ``barcodetag`` of ``CB``.
  - ``mapq``: Controls quality filtering settings for generating the fragments file. Discards reads with quality score lower than this number (default 30).


Execution
---------

After configuring, the workflow can be run with:

.. code:: bash

    nextflow -C atac_preprocess.config run \
        vib-singlecell-nf/vsn-pipelines/main_atac.nf \
        -entry atac_preprocess -resume


Output
******

An example output tree is shown here.

.. code:: bash

    out/
    ├── data
    │   ├── bam
    │   │   ├── sample_1.bwa.out.possorted.bam
    │   │   ├── sample_1.bwa.out.possorted.bam.bai
    │   │   ├── sample_2.bwa.out.possorted.bam
    │   │   └── sample_2.bwa.out.possorted.bam.bai
    │   ├── fragments
    │   │   ├── sample_1.sinto.fragments.tsv.gz
    │   │   ├── sample_1.sinto.fragments.tsv.gz.tbi
    │   │   ├── sample_2.sinto.fragments.tsv.gz
    │   │   └── sample_2.sinto.fragments.tsv.gz.tbi
    │   └── reports
    │       ├── barcode
    │       │   ├── sample_1____S7_R1_001.corrected.bc_stats.log
    │       │   └── sample_2____S8_R1_001.corrected.bc_stats.log
    │       ├── mapping_stats
    │       │   ├── sample_1.mapping_stats.tsv
    │       │   └── sample_2.mapping_stats.tsv
    │       ├── mark_duplicates
    │       │   ├── sample_1.library_complexity_metrics.txt
    │       │   ├── sample_1.mark_duplicates_metrics.txt
    │       │   ├── sample_2.library_complexity_metrics.txt
    │       │   └── sample_2.mark_duplicates_metrics.txt
    │       └── trim
    │           ├── sample_1____S7_R1_001.fastp.trimming_report.html
    │           └── sample_2____S8_R1_001.fastp.trimming_report.html
    └── nextflow_reports
        ├── execution_report.html
        ├── execution_timeline.html
        ├── execution_trace.txt
        └── pipeline_dag.dot

----

