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

Examples
--------

Below are examples using the `capture` pipeline for different scenarios

.. code-block:: bash
    :caption: 50bp oligos for NlaIII fragments in hg19 build, using STAR to check off-target binding

    python capture.py -f ~/hg19/Sequence/genome.fa -g hg19 -b viewpoints.bed -o 50 -e NlaIII -s ~/hg19/STAR/
    
.. code-block:: bash
    :caption: 70bp oligos for HindIII fragments in mm10 build, using BLAT to check off-target binding

    python capture.py -f ~/mm10/Sequence/genome.fa -g mm10 -b mouse_viewpoints.bed -e HindIII --blat

Specifics
---------

**Bed file** (-b, \\--bed)
    A 4-column, tab-delimited bed file containing the coordinates and names of the viewpoints you want to capture from. This must be in the format `chr`, `start`, `stop`, `viewpoint_name`. Typically, the coordinates each span 1bp as shown below:

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

**STAR** (-s, \\--star_index) or **BLAT** (\\--blat)
    To check for off-target binding, either the sequence aligner STAR or the BLAST-like Alignment Tool (BLAT) can be used. By default, STAR is used, unless `capture.py` is run with the `\\--blat` flag. Since BLAT is more widely used to detect off-target binding events, this might be preferred
    by the user. However, BLAT can be particulary slow for large designs, especially for the human reference genomes. STAR's exceptional speed is better suited for designs with >500 viewpoints. If the `\\--blat` flag is not selected, the path to the STAR index must be supplied
    after the `-s` (or `\\--star_index`) flag.

Requirements
------------

| Python 3.4+, pysam, numpy, biopython
| RepeatMasker
| STAR or BLAT (:ref:`see above <star-blat>`)

Functions
=========

As well as being run as a full pipeline from the command line, the `oligo` modules have been written such that the individual functions can be easily run in a python shell. The pipeline runs the functions in the following order:

#. :func:`capture.gen_oligos`
#. :func:`tools.write_oligos`
#. :func:`tools.check_off_target`
#. :func:`tools.get_density`

Below is a detailed list of functions in the `capture` module:

.. autofunction:: gen_oligos()




