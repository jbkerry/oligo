#########
Capture-C
#########

.. currentmodule:: design

:doc:`design.Capture <capture_class>`

Description
===========

The image below shows a schematic of how the `Capture` pipeline designs oligos adjacent to the first restriction site of a specified restriction enzyme (DpnII in this example), on the left- and right-hand sides. In this case the user has supplied viewpoint
coordinates at chr2:5500000-5500001 and chr5:63223000-63223001.

.. figure:: _static/capture_oligo_gen.png

    Schematic of oligo design by `Capture`
    
It is possible for the designed oligos to overlap, up to within 1bp of each other, when the restriction fragment is less than twice the length of the oligo. If the viewpoint coordinate is in a fragment with length less than the specified oligo length, no oligos will be generated for that fragment.
If the fragment length exactly equals the oligo length, only one oligo will be generated.

.. note::
    
    For full functionality, `Capture` should be run from the command line in order to test the efficiency of the generated oligos. This involves a pipeline that incorporates methods from the :doc:`Tools <tools_class>` class.

Usage
=====

When run from the command line, ``oligo capture`` takes the following parameters:

.. option:: -h, --help
    
    (flag) Show this help message and exit
    
.. option:: -f <reference fasta>, --fasta <reference fasta>
    
    (str) The path to the reference genome fasta
    
.. option:: -g <genome>, --genome <genome>

    ({mm9, mm10, hg18, hg19, hg38}) The name of the genome build
    
.. option:: -b <bed file>, --bed <bed file>

    (str) The path to the :ref:`bed file <bed-file>` containing the capture viewpoint coordinates
    
.. option:: -o <oligo length>, --oligo <oligo length>

    (int, optional) The length (bp) of the oligos to design, default=70
    
.. option:: -e <enzyme>, --enzyme <enzyme>

    ({DpnII, NlaIII, HindIII}, optional) Name of the :ref:`restriction enzyme <enzyme>` to be used
    for fragment digestion, default=DpnII
    
.. option:: -s <STAR index>, --star_index <STAR index>

    (str) The path to the STAR index directory; omit this option if running with BLAT (:option:`--blat`)
    
.. option:: --blat

    (flag) Detect off-target binding using :ref:`BLAT instead of STAR <star-blat>`

Examples
--------

Below are examples using the `Capture` pipeline for different scenarios

.. code-block:: bash
    :caption: 50bp oligos for NlaIII fragments in hg19 build, using STAR to check off-target binding

    $ python -m oligo capture -f ~/hg19/Sequence/genome.fa -g hg19 -b viewpoints.bed -o 50 -e NlaIII -s ~/hg19/STAR/
    
.. code-block:: bash
    :caption: 70bp oligos for HindIII fragments in mm10 build, using BLAT to check off-target binding

    $ python -m oligo capture -f ~/mm10/Sequence/genome.fa -g mm10 -b mouse_viewpoints.bed -e HindIII --blat

Specifics
---------

.. _bed-file:

**Bed file** (:option:`-b`, :option:`--bed`)
    A 4-column, tab-delimited bed file containing the coordinates and names of the viewpoints you want to capture from. This must be in the format `chr`, `start`, `stop`, `viewpoint_name` (the file should not have a
    header row). Typically, the coordinates each span 1bp as shown below:

.. csv-table::
   :class: center-table
   
   chr7, 20205, 20206, geneX
   chr8, 1310000, 1310001, geneY

.. caution::

   Names in the last column should be unique so that the oligos can be unambiguosuly linked back to a named viewpoint.
   
.. _enzyme:

**Restriction enzyme** (:option:`-e`, :option:`--enzyme`)
    The restriction enzyme being used in the Capture-C experiment. This determines the recogition sequence used to define the fragment boundaries and hence the starts and ends of the oligos. The current version supports `DpnII` (GATC), `NlaIII` (CATG) and `HindIII` (AAGCTT).
    If this option is omitted, `DpnII` will be used by default.

.. _star-blat:

**STAR** (:option:`-s`, :option:`--star_index`) or **BLAT** (:option:`--blat`)
    To check for off-target binding, either the sequence aligner STAR or the BLAST-like Alignment Tool (BLAT) can be used. By default, STAR is used, unless `design.py Capture` is run with the :option:`--blat` flag. Since BLAT is more widely used to detect off-target binding events, this might be preferred
    by the user. However, BLAT can be particulary slow for large designs, especially for the human reference genomes. STAR's exceptional speed is better suited for designs with >500 viewpoints (1000 oligos). If the :option:`--blat` flag is not selected, the path to the STAR index must be supplied
    after the :option:`-s` (or :option:`--star_index`) flag.

API
===

As well as being run as a full pipeline from the command line, the `design` classes have been written such that the individual methods can be easily run in a python shell. The `Capture` pipeline implements methods from :doc:`design.Capture <capture_class>`.
The following examples show the order in which the class methods are implemented:

.. code-block:: python
    :caption: Create an instance of the Capture class

    >>> from oligo.design import Capture
    >>> c = Capture(genome='mm9', fa='mm9_genome.fa')
    Loading reference fasta file...
        ...complete
    
.. code-block:: python
    :caption: Generate oligos and write to fasta file

    >>> c.gen_oligos(bed='capture_sites.bed').write_oligos()
    Generating oligos...
        ...complete.
    Oligos stored in the oligo_seqs attribute
    Wrote oligos to oligo_seqs.fa
    
.. code-block:: python
    :caption: Check for repeats in oligo sequences and align oligos to genome

    >>> c.detect_repeats().align_to_genome(s_idx='/mm9/STAR')
    Checking for repeat sequences in oligos, with RepeatMasker...
        ...complete. Output written to oligo_seqs.fa.out
    Aligning oligos to the genome, with STAR...
        ...complete. Output written to oligos_Aligned.out.sam

.. code-block:: python
    :caption: Calculate longest repeat length and number of off-target alignments for each oligo; write to oligo_info.txt

    >>> c.extract_repeats().calculate_density().write_oligo_info()
    Repeat scores calculated
    Density scores calculated
    Oligo information written to oligo_info.txt

See :doc:`design.Capture <capture_class>` for more detailed information

.. centered:: :doc:`Top of Page <capture>`

