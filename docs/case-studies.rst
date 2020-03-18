Case Studies
=============


Kurmangaliyev et al., 2019
--------------------------

A single cell analysis of transcriptional control of neuronal connectivity in Drosophila,
based on `Kurmangaliyev et al., 2019 <https://elifesciences.org/articles/50822>`_.

.. note:: `Full tutorial here <https://vsn-pipelines-examples.readthedocs.io/en/latest/Kurmangaliyev.html>`_.

This case study illustrates the following steps:

1. **Input data** is loaded directly from the `Sequence Read Archive (SRA) <https://www.ncbi.nlm.nih.gov/sra>`_ by giving an SRA identifier to the ``sra`` input channel.
2. Cell Ranger is run to generate expression counts
3. Multiple samples are combined, and **batch effect correction** is performed with both BBKNN and Harmony (in separate pipeline runs).
4. **Gene regulatory network inference** is performed using the SCENIC pipeline. The SCENIC append mode is used to include the SCENIC results with both independent batch effect correction methods, to avoid re-running SCENIC.


Hung et al., 2019
-----------------

A single cell analysis of the adult Drosophila midgut, based on
`Hung et al., 2019 <https://vsn-pipelines-examples.readthedocs.io/en/latest/PBMC10k.html>`_.

.. note:: `Full tutorial here <https://vsn-pipelines-examples.readthedocs.io/en/latest/Hung.html>`_.

This case study illustrates the following steps:

1. **Input data** is loaded directly from the `Sequence Read Archive (SRA) <https://www.ncbi.nlm.nih.gov/sra>`_ by giving an SRA identifier to the ``sra`` input channel.
2. Cell Ranger is run to generate expression counts
3. Multiple samples are combined, and **batch effect correction** is performed with BBKNN
4. **Gene regulatory network inference** is performed using the SCENIC pipeline.


PBMC10k
-------

An analysis of a sample dataset from 10x Genomics consisting of 10,000 PBMCs from a healthy human donor.

.. note:: `Full tutorial here <https://vsn-pipelines-examples.readthedocs.io/en/latest/PBMC10k.html>`_.

This case study illustrates the following steps:

1. **Input data** is filtered Cell Ranger counts downloaded from the 10x Genomics support website.
2. The single sample is run through the standard ``single_sample`` pipeline.
3. **Gene regulatory network inference** is performed using the SCENIC pipeline and integrated with the highly variable genes analysis.

