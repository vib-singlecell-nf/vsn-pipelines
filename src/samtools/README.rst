
Samtools Docker images
======================

This directory contains Dockerfiles for base images used here and for other images in the VSN Pipelines repository.


To build the Base image
-----------------------

This base image is based on ``debian:buster-slim`` and has a compiled verison of 
`zlib-ng <https://github.com/zlib-ng/zlib-ng>`_ for faster compression and decompression.

Image tag format: simple version numbers (0.1, 0.2, ...).

.. code:: bash

    docker build -t vibsinglecellnf/samtools:base-0.3 . -f Dockerfile.samtools-base
    podman build -t vibsinglecellnf/samtools:base-0.3 . -f Dockerfile.samtools-base

This base image is used in several other images within VSN::
    
- samtools [this directory]


To build the Samtools image
---------------------------

This uses the base image above and adds Samtools and HTSlib

Image tag format: ``<base image version>-<samtools release version>``.

.. code:: bash

    docker build -t vibsinglecellnf/samtools:0.3-1.15.1 .
    podman build -t vibsinglecellnf/samtools:0.3-1.15.1 .

This samtools image is used in several other images within VSN::
    
- singlecelltoolkit
- bwamaptools
- popscle


