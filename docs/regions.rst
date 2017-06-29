#######
Regions
#######

.. toctree::


Description
===========

.. automodule:: regions
   :platform: Unix
   
Functions: :func:`gen_oligos`

.. note::
    
    For full functionality, regions.py should be run from the command line in order to test the efficiency of the generated oligos. This involves a pipeline that incorporates methods from the tools module.
    This is the intended use but can also be run in its constituent parts and modules are listed throughout the doc pages

Usage
=====

When run from the command line, *Regions* takes the following parameters

.. code-block:: bash
   :caption: Parameters for regions.py

   usage: regions.py [-h] -f FASTA -g GENOME -b BED [-o OLIGO] [-e ENZYME]
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

-f
    Path to your reference genome fasta file. This is typically a file called **genome.fa**.
-g
    The genome build which you are generating oligos for, e.g. **mm10** or **hg38**. This must match with the genome file you have supplied for the previous option.
-b
    A tab-delimited 4-column bed file containing the coordinates of viewpoint you want to capture from. This must be in the format **chr**, **start**, **stop**, **viewpoint_name**. The coordinates are typcially 1bp e.g.


+---------+---------+---------+---------+
| chr7    | 20000   | 20001   | geneX   |
+---------+---------+---------+---------+ 


.. autofunction:: gen_oligos()

.. function:: gen_oligos()
   Explanation here
   
The :func:`gen_oligos` can be used for generate oligos

Examples
========

Below are examples using the *Regions* module for different scenarios

