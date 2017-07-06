.. _top:

#####
oligo
#####

* :ref:`Overview <overview>`
* :ref:`Output <output>`
    * :ref:`Missing Values <missing>`
* :ref:`Choosing Good Oligos <filtering>`
    * :ref:`Density Score <density>`
    * :ref:`Repeat Length <repeat>`
* :ref:`genindex`

.. _overview:

.. include:: ../README.rst
    :end-before: For full documentation

More detailed usage information can be found in the individual pages, via the navigation on the left. A schematic of the pipeline workflows is shown below.

.. image:: _static/oligo_flow.png

.. container:: top-link

    :ref:`top of page <top>`

.. _output:

Output
======

Although a number of files are output from the pipelines, the important one is `oligo_info.txt`. This is a tab-delimited text file that contains the following information for every oligo:

**chr**
    the chromsome that the oligo/fragment is on
**start**
    the bp start coordinate of the oligo
**stop**
    the bp stop coordinate of the oligo
**fragment_start**
    the bp start coordinate of the fragment
**fragment_stop**
    the bp stop coordinate of the fragment
**side_of_fragment**
    the side of the fragment that the oligo is situated on (left or right)
**sequence**
    the full DNA sequence of the oligo
**total_number_of_alignments**
    the total number of times the oligo was found in the BLAT/STAR alignment file
**density_score**
    the base-pair average and length-normalised number of times the oligo was found in the BLAT/STAR alignment file (see :ref:`Choosing Good Oligos <filtering>` for more details)
**repeat_length**
    the length of the longest simple sequence repeat found in the oligo
**repeat_class**
    the class of the longest simple sequence repeat found in the oligo
**GC%**
    the GC percentage of the oligo sequence
**associations**
    the position/gene name associated with this oligo; this is the viewpoint name supplied in the bed file (4th column)
    
.. _missing:

Missing Values
--------------

Due to particular differences between the three pipelines, and in order to keep a consistent output format between the three, there are instances where some values in the file will be purposefully missing.

**fragment_start, fragment_stop, side_of_fragment**
    These values will be replaced with a '.' for the Tiled Capture pipeline when run in FISH mode, and for the CRISPR Off-Target pipeline, as both of these pipelines are restriction fragment-independent
**total_number_of_alignments, density_score**
    If any of the pipelines are run using STAR to check for off-target binding there may be some oligos that have a '0' value for both of these columns. This represents a bad oligo because
    it either could not be mapped by STAR or it mapped so many times that it passed the threshold that means a read is ignored. These oligos should be filtered out (see :ref:`Choosing Good Oligos <filtering>` for more details).
    The best density_score value for an oligo is 1.0 as this means it mapped to only one (the expected) location in the genome.
**repeat_length, repeat_class**
    If oligos have '0', and 'N/A' values for these columns, respectively, it means that no simple sequence repeats were detected in that oligo
**association**
    This value will be replaced with a '.' for the Tiled Capture pipeline as these oligos are generated for adjacent sites across one large region and not for different viewpoints associated with unique names 

.. _filtering:

Choosing Good Oligos
====================

.. note::

    Cut-offs for efficient oligos:

    | 1 <= **density_score** <= 30 (if using BLAT)
    | 1 <= **density_score** <= 50 (if using STAR)
    | **repeat_length** <= oligo_length/4
    
.. _density:

Density Score
-------------

When performing a Capture experiment it is important that the oligo is not susceptible to off-target binding as this would result in the pulldown of unwanted material, generating false positives. To assess the degree of off-target binding, each oligo is aligned against the
genome using either STAR or BLAT. The number of times an oligo aligns to the genome is used to assess its degree of off-target binding and this value is represented by the density score. Density score is calculated by counting the base-pair coverage of each alignment for an oligo
and dividing that value by the length of the oligo. A schematic of this calculation is shown below, for two 70bp oligos.

.. image:: _static/density_score.png

Therefore, a lower density score, the better, with the exception of values below 1. A value between 0 and 1 means no perfect alignment was found so the pulldown will not be 100% efficient. A score of 0 means that either no alignments were found so this oligo will not capture
anything or, in the case of using STAR, can also mean that it mapped so many times that it passed the threshold for a read to be ignored. A safe cut-off for density score to ensure efficient oligo capture without problematic off-target binding is greater than or equal to 1 and
less than or equal to 30 (if using BLAT), or less than or equal to 50 (if using STAR).

.. _repeat:

Repeat Length
-------------

The presence of simple sequence repeats in an oligo can also cause a greater degree of off-target binding due to ambiguous presence of these repeats throughout the genome. Usually if an oligo has a long repeat within it, it will also have a high density score. However, it
is still best to filter on repeat length as well. We filter for oligos that have a repeat length less than or equal to a quarter of the length of the oligo, so for an 80bp oligo only those with a repeat length <= 20 would be accepted. Unlike density score, 0 is the best
value as this means that the oligo does not contain any simple sequence repeats.

.. centered:: :ref:`Top of Page <top>`

.. toctree::
    :hidden:
   
    Capture-C <capture>
    Tiled Capture <tiled>
    CRISPR Off-Target <off_target>
    tools
    