scATAC-seq Pipelines
====================

----

scATAC-seq preprocessing
************************

This pipeline takes fastq files from paired end single cell ATAC-seq, and applies preprocessing steps to align the reads to a reference genome, and produce a bam file and scATAC-seq fragments file.
The full steps are:

- Barcode correction:

  * For 'standard' and 'multiome' samples (e.g. 10x Genomics) correction is performed against a whitelist by `this script <https://github.com/aertslab/single_cell_toolkit/blob/master/correct_barcode_in_fastq.sh>`_.
  * For 'biorad' samples, barcode correction is performed by `BAP <https://github.com/caleblareau/bap>`_.

- Debarcoding: Add the barcode sequence to the beginning of the fastq sequence identifier
- Read/adapter trimming
- Mapping to a reference genome:

  * ``bwa mem`` is used with default parameters.
  * Duplicates are marked with ``samtools markdup``.
  * Droplet barcodes are included in the BAM file with the ``CR`` tag (by default). No barcode correction is performed.

- A fragments file is created using `Sinto <https://github.com/timoast/sinto>`_.

Input
-----

The input to this pipeline is a (tab-delimited) metadata table with the sample ID, sequencing technology, and locations of the fastq files:

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

The columns represent:

- ``sample_name`` Sample name for labeling the sample in the pipeline and output files. This can be any arbitrary string.
- ``technology``: This described the barcode correction and processing methods to use for the fastq files. Current options are ``standard``, ``multiome``, or ``biorad``. See below for additional details.
- ``fastq_PE1_path``: The full path to the fastq file for the first read in a pair.
- ``fastq_barcode_path``: The full path to the fastq file containing the barcodes. This column can be blank/empty depending on the technology setting.
- ``fastq_PE2_path``: The full path to the fastq file for the second read in a pair.

Technology
----------

This controls how both barcode correction and debarcoding is applied to the input fastq files.
Available options are:

``standard`` 
____________

The ``standard`` setting assumes a typical 10x Genomics style format with two read pair fastqs and a barcode fastq:

.. code:: none

    $ zcat sample_1_R1.fastq.gz | head -n 4
    @A00311:74:HMLK5DMXX:1:1101:2013:1000 1:N:0:ACTCAGAC
    NTTGTCTCAGCACCCCCCGACATGGATTCAGGCTGTCTCTTATACACATC
    +
    #FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

    $ zcat sample_1_R2.fastq.gz | head -n 4
    @A00311:74:HMLK5DMXX:1:1101:2013:1000 2:N:0:ACTCAGAC
    CTGTTCGCAAAGCATA
    +
    F:FFFFFFFFFFFFFF

    $ zcat sample_1_R3.fastq.gz | head -n 4
    @A00311:74:HMLK5DMXX:1:1101:2013:1000 3:N:0:ACTCAGAC
    CCTGAATCCATGTCGGGGGGTGCTGAGACAAGCTGTCTCTTATACACAT
    +
    FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

The debarcoding step here uses a 
`helper script <https://github.com/aertslab/single_cell_toolkit/blob/master/debarcode_10x_scatac_fastqs.sh>`_
which transforms this input into two paired fastq files with the barcode integrated into the read name:

.. code:: none

    $ zcat sample_1_dex_R1.fastq.gz | head -n 4
    @CTGTTCGCAAAGCATA:A00311:74:HMLK5DMXX:1:1101:2013:1000 1:N:0:ACTCAGAC
    NTTGTCTCAGCACCCCCCGACATGGATTCAGGCTGTCTCTTATACACATC
    +
    #FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

    $ zcat sample_1_dex_R2.fastq.gz | head -n 4
    @CTGTTCGCAAAGCATA:A00311:74:HMLK5DMXX:1:1101:2013:1000 3:N:0:ACTCAGAC
    CCTGAATCCATGTCGGGGGGTGCTGAGACAAGCTGTCTCTTATACACAT
    +
    FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF


``multiome``
____________

The ``multiome`` setting works the same as ``standard`` with the exception of the whitelist used for barcode correction.
The whitelists are supplied in the params file (``params.tools.singlecelltoolkit.barcode_correction.whitelist``).


``biorad`` 
__________

The ``biorad`` setting processes BioRad data using `BAP <https://github.com/caleblareau/bap/wiki/Working-with-BioRad-data>`_.
This takes input data:

.. code:: none

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


And produces paired fastq files with the barcode integrated into the read name (with a ``_`` delimiter):

.. code:: none

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


Running the workflow
--------------------

To generate a config file, use the ``atac_preprocess`` profile along with ``docker`` or ``singularity``.
Note that the full path to ``vib-singlecell-nf/vsn-pipelines/main_atac.nf`` must be used:

.. code:: bash

    nextflow config \
        vib-singlecell-nf/vsn-pipelines/main_atac.nf \
        -profile atac_preprocess,singularity \
        > atac_preprocess.config

The ATAC-specific parameters are described here.
The important parameters to change are:

- ``params.data.atac_preprocess.metadata``: the path to the metadata file.
- ``params.tools.bwamaptools.bwa_fasta``: the path to the bwa reference fasta file. This should be already indexed with ``bwa index``, and the index files located in the same directory as the fasta file.
- ``params.tools.singlecelltoolkit.barcode_correction.whitelist``: Whitelists for barcode correction are supplied here. The whitelists are matched to samples based on the parameter key here ('standard', 'multiome') and the technology field listed for each sample in the metadata file.

Optional parameters to change:

- Within ``params.tools.bwamaptools.add_barcode_as_tag``:

  - ``tag``: controls the naming of the barcode tag added to the bam (``CR`` by default).
  - ``delimiter_to_split_qname``: Controls which delimiter to split the bam read name field to get the barcode. By default it uses the regex ``'[:|_]'`` to split on both ``:`` and ``|``.

- Within ``params.tools.sinto.fragments``:

  - One of (but not both) ``barcodetag`` or ``barcode_regex`` needs to be set to tell Sinto where to find the barcodes in the bam file. The default is to use ``barcodetag`` of ``CR``.
  - ``mapq``: Controls quality filtering settings for generating the fragments file. Discards reads with quality score lower than this number (default 30).
  - ``temp_dir``: Controls where temp files are stored during fragments processing. For large BAM files, the system default temp location may become full. An alternate temp path can be specified here. Be sure to also include this temp path in the global volume mounts for Docker/Singularity in the config file.


After configuring, the workflow can be run with:

.. code:: bash

    nextflow -C atac_preprocess.config run \
        vib-singlecell-nf/vsn-pipelines/main_atac.nf \
        -entry atac_preprocess -resume

----
