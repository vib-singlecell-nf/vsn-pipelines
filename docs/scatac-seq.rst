scATAC-seq Preprocessing
========================


This pipeline takes fastq files from paired end single cell ATAC-seq, and applies preprocessing steps to align the reads to a reference genome, and produce a bam file and scATAC-seq fragments file.

This workflow is currently available in the ``develop_atac`` branch (use the ``-r develop_atac`` option when running ``nextflow pull`` and ``nextflow run``).

----

Pipeline Steps
**************

The full steps are:

- Barcode correction:

  * For 'standard' and 'multiome' samples (e.g. 10x Genomics or similar) correction is performed against a whitelist by 
    `this method <https://github.com/aertslab/single_cell_toolkit/blob/master/correct_barcode_in_fastq.sh>`_ 
    from `aertslab/single_cell_toolkit <https://github.com/aertslab/single_cell_toolkit>`_.
  * For 'biorad' samples, barcode correction is performed by `BAP <https://github.com/caleblareau/bap>`_.

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

Pipeline Input Metadata
***********************

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

Note that there is an easy way to create the metadata from the file paths for each sample by using the following bash command (expand to view):

.. raw:: html

   <details>
   <summary><a>metadata generator</a></summary>

.. code-block:: none

    create_atac_metadata() {
        local sample="${1}"
        local technology="${2}"
        local file_prefix="${3}"
        local read_labels="${4}"
        if [ "${sample}" == "header" ]; then
            echo -e "sample_name\ttechnology\tfastq_PE1_path\tfastq_barcode_path\tfastq_PE2_path"
            return 1
        fi
        read_labels_arr=(${read_labels//,/ })
        R1=(${file_prefix}*${read_labels_arr[0]}*)
        R2=(${file_prefix}*${read_labels_arr[1]}*)
        R3=(${file_prefix}*${read_labels_arr[2]}*)
        for i in "${!R1[@]}"; do
            echo -e "${sample}\t${technology}\t${R1[i]}\t${R2[i]}\t${R3[i]}";
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

.. raw:: html

   </details>

----

Technology
----------

The "technology" field in the metadata table controls how technology-specific pipeline steps are applied, as well as which whitelist is used for barcode correction.
Currently the only the ``biorad`` setting uses alternate pipelines processes (to extract and correct the barcode sequence from the two input fastqs).
Except for the ``biorad`` setting, the samples will be processed in the standard pipeline (barcode correction against a whitelist).

The "technology" field can be set to any string (e.g. ``standard``), but note that the entry in this field must match the barcode label given in the ``params.tools.singlecelltoolkit.barcode_correction.whitelist`` parameter.
Commonly used default settings are:

``standard`` 
____________

The ``standard`` setting assumes a typical 10x Genomics style format with two read pair fastqs and a barcode fastq (note here that the barcode correction has already been performed, writing the ``CB`` tag into the comment of the barcode fastq)::

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


``biorad`` 
__________

The ``biorad`` setting processes BioRad data using `BAP <https://github.com/caleblareau/bap/wiki/Working-with-BioRad-data>`_.
This takes input data::

    $ zcat sample_2_R1.fastq.gz | head -n 4
    @NB551608:167:HNYFJBGXC:1:11101:11281:1033 1:N:0:TAAGGCGA
    GCGTANACGTATGCATGACGGAAGTTAGTCACTGAGTCAGCAATCGTCGGCAGCGTCAGATGAGTNTAAGAGACAGGGTCAGGATGCGAGATTGACGGCTGCAATAACTAATAGGAAC
    +
    AAAAA#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEE<EEEE6EA#6E<66AAEEEEEAEEEEEEEEEEEEAEEAEEEEEEEEE<EEEEEEEEEEE/E

    $ zcat sample_2_R2.fastq.gz | head -n 4
    @NB551608:167:HNYFJBGXC:1:11101:11281:1033 2:N:0:TAAGGCGA
    NNGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    +
    ##A####################################


And produces paired fastq files with the barcode integrated into the read name (with a ``_`` delimiter)::

    $ zcat sample_2_dex_R1.fastq.gz | head -n 4
    @GCGTAGAGGAAGTTTCAGCAA_NB551608:167:HNYFJBGXC:1:11101:11281:1033 1:N:0:TAAGGCGA
    GGTCAGGATGCGAGATTGACGGCTGCAATAACTAATAGGAAC
    +
    EEAEEEEEEEEEEEEAEEAEEEEEEEEE<EEEEEEEEEEE/E

    $ zcat sample_2_dex_R2.fastq.gz | head -n 4
    @GCGTAGAGGAAGTTTCAGCAA_NB551608:167:HNYFJBGXC:1:11101:11281:1033 2:N:0:TAAGGCGA
    NNGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    +
    ##A####################################


----

Running the workflow
********************

Configuration
-------------

To generate a config file, use the ``atac_preprocess`` profile along with ``docker`` or ``singularity``.
Note that the full path to ``vib-singlecell-nf/vsn-pipelines/main_atac.nf`` must be used:

.. code:: bash

    nextflow config \
        vib-singlecell-nf/vsn-pipelines/main_atac.nf \
        -profile atac_preprocess,singularity \
        > atac_preprocess.config


Parameters
----------

The ATAC-specific parameters are described here.
The important parameters to verify are:

- ``params.data.atac_preprocess.metadata``: the path to the metadata file.
- ``params.tools.bwamaptools.bwa_fasta``: the path to the bwa reference fasta file. This should be already indexed with ``bwa index``, and the index files located in the same directory as the fasta file. Note that ``bwa`` and ``bwa-mem2`` use different indexes that are not interchangeable.
- ``params.tools.singlecelltoolkit.barcode_correction.whitelist``: Whitelists for barcode correction are supplied here. The whitelists are matched to samples based on the parameter key here ('standard', 'multiome', etc.) and the technology field listed for each sample in the metadata file.

Choice of tools
_______________

Several steps have options for the choice of method to use.
These options are controlled within ``params.atac_preprocess_tools``.

- Adapter trimming (``adapter_trimming_method``): Can be either of ``Trim_Galore`` (default), or ``fastp``.
- Duplicate marking (``mark_duplicates_method``): Can be either of ``MarkDuplicates`` (Picard tools, default) or ``MarkDuplicatesSpark`` (GATK). We currently recommend Picard MarkDuplicates because it has the capability to perform barcode-aware marking of PCR duplicates. MarkDuplicatesSpark has the advantage of parallelization, however it requires a large SSD to use for temporary files.

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

----

Other considerations
--------------------

Temporary directory mapping
___________________________

For large BAM files, the system default temp location may become full.
   A workaround is to include a volume mapping to the alternate ``/tmp`` ``-B /alternate/path/to/tmp:/tmp`` using the volume mount options in Docker or Singularity.
   For example in the container engine options:
  - Singularity run options: ``runOptions = '--cleanenv -H $PWD -B /data,/alternate/path/to/tmp:/tmp'``
  - Docker run options: ``runOptions = '-i -v /data:/data -v /alternate/path/to/tmp:/tmp'``

Alternate Nextflow work location
________________________________

Direct the Nextflow work directory to an alternate path (e.g. a scratch drive) using the ``NXF_WORK`` environmental variable::

    nwork=/path/to/scratch/example_project
    mkdir $nwork
    export NXF_WORK=$nwork

