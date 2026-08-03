.. Parachrsnp documentation master file, created by
   sphinx-quickstart on Mon Aug  3 19:06:42 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Parachrsnp
========================

ParaChrSNP is a containerized Snakemake workflow for chromosome-wise SNP discovery and population-genomic analysis from paired-end resequencing data. It streamlines the full process from read quality control and alignment to cohort-level variant calling, filtering, annotation, format conversion, optional genotype imputation, CNV detection and downstream population analyses. The central design of ParaChrSNP is per-sample by chromosome parallel variant calling, which splits large variant-calling tasks into independent chromosome-level units and then integrates them into a cohort VCF for downstream analysis.

.. image:: Parachrsnp.png
     :alt: ParaChrSNP workflow
     :width: 100%
     :align: center

.. toctree::
   :maxdepth: 2
   :caption: Contents:

