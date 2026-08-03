ParaChrSNP
==========

ParaChrSNP is a containerized Snakemake workflow for chromosome-wise SNP
discovery and population-genomic analysis from paired-end resequencing data.
It connects read quality control, alignment, duplicate removal, per-sample and
per-chromosome variant calling, cohort joint calling, filtering, annotation,
format conversion, and optional population-genomic analyses in one workflow.

.. image:: Parachrsnp.png
     :alt: ParaChrSNP workflow
     :width: 100%
     :align: center

Performance benchmarks
----------------------

ParaChrSNP was benchmarked with four plant genomes at 10x, 30x, 50x, and
100x sequencing depths on a server with 512 GB RAM and 64 CPU cores. The
table reports the observed running time in minutes (min).

.. list-table:: ParaChrSNP running time at different sequencing depths
   :header-rows: 1
   :widths: 24 16 15 15 15 15
   :align: center

   * - Species
     - Genome size
     - 10x (min)
     - 30x (min)
     - 50x (min)
     - 100x (min)
   * - Arabidopsis
     - 0.12 Gb
     - 2.72
     - 8.56
     - 13.98
     - 25.75
   * - Watermelon
     - 0.37 Gb
     - 8.16
     - 22.75
     - 39.71
     - 82.64
   * - Soybean
     - 0.98 Gb
     - 24.82
     - 74.26
     - 128.81
     - 251.53
   * - Tea
     - 2.96 Gb
     - 71.64
     - 201.31
     - 394.05
     - 746.34

Actual running time depends on the number of samples, read length, storage
performance, CPU and memory resources, selected joint-calling method, and
enabled optional modules.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation/index
   quick_start/index
   usage/index
   web_interface/index
   changelog/index
   faq/index
