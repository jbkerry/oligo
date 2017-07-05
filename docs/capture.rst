#########
Capture-C
#########

Description
===========

.. automodule:: capture
   :platform: Unix
   
Functions: :func:`gen_oligos`

.. note::
    
    For full functionality, `capture` should be run from the command line in order to test the efficiency of the generated oligos. This involves a pipeline that incorporates functions from the :doc:`tools <tools>` module.

Usage
=====

When run from the command line, `capture.py` takes the following parameters

.. code-block:: bash
   :caption: Parameters for capture.py

   usage: capture.py [-h] -f FASTA -g GENOME -b BED [-o OLIGO] [-e ENZYME]
                  [-s STAR_INDEX] [--blat]

   optional arguments:
      -h, --help            show this help message and exit
      -f FASTA, --fasta FASTA
                            Path to reference genome fasta.
      -g GENOME, --genome GENOME
                            Genome build e.g. 'mm10' or 'hg38'.
      -b BED, --bed BED     Path to bed file with capture viewpoint coordinates
      -o OLIGO, --oligo OLIGO
                            The size (in bp) of the oligo to design, default=70
      -e ENZYME, --enzyme ENZYME
                            Name of restriction enzyme, default=DpnII
      -s STAR_INDEX, --star_index STAR_INDEX
                            Path to STAR index directory. Omit this option if
                            running with BLAT (--blat)
      --blat                Detect off-targets using BLAT instead of STAR.

Requirements
------------

| Python 3.4+, pysam, numpy, biopython
| RepeatMasker
| STAR or BLAT (:ref:`see below <star-blat>`)

Specifics
---------

**Bed file** (-b, \\--bed)
    A 4-column, tab-delimited bed file containing the coordinates and names of the viewpoints you want to capture from. This must be in the format `chr`, `start`, `stop`, `viewpoint_name`. Typically, the coordinates should each span 1bp as shown below:

.. csv-table::
   :align: right
   
   chr7, 20205, 20206, geneX
   chr8, 1310000, 1310001, geneY

.. warning::

   Names in the last column should be unique so that the oligos can be unambiguosuly linked back to a named viewpoint.

**Restriction enzyme** (-e, \\--enzyme)
    The restriction enzyme being used in the Capture-C experiment. This determines the recogition sequence used to define the fragment boundaries and hence the starts and ends of the oligos. The current version supports `DpnII` (GATC), `NlaIII` (CATG) and `HindIII` (AAGCTT).
    If this option is omitted, `DpnII` will be used by default.

.. _star-blat:

**STAR** (-s, \\--s_idx) or **BLAT** (\\--blat)
    To check for off-targets, either the sequence aligner STAR or the BLAST-like Alignment Tool (BLAT) can be used. By default, STAR is used, unless `capture.py` is run with the `\\--blat` flag. Since BLAT is more widely used to detect off-target binding events, this might be preferred
    by the user. However, BLAT can be particulary slow for large designs, especially for the human reference genomes. STAR's exceptional speed is better suited for designs with >500 viewpoints. If the `\\--blat` flag is not selected, the path to the STAR index must be supplied
    after the `-s` (or `\\--s_idx`) flag.
   
Examples
========

Below are examples using the `capture` module for different scenarios

Functions
=========

.. autofunction:: gen_oligos()




