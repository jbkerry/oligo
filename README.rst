Capture Oligo Design
====================

The *oligo* python package consists of three modules for Capture-C and FISH oligo design.

* `capture.py <http://oligo.rtfd.io/en/latest/capture.html>`_: designs oligos for a standard Capture-C experiment. The user supplies a list of viewpoint coordinates, and oligos are generated adjacent to the flanking recognition sequence of a specified restriction enzyme.
* `tiled.py <http://oligo.rtfd.io/en/latest/tiled.html>`_: designs oligos for multiple adjacent restriction fragments across a specified region of a chromosome, or for the entire chromosome. If *tiled* is run in FISH mode, oligos are generated independent of restriction fragments and
  are instead generated for a user-specified step size, in an adjacent manner.
* `off_target.py <http://oligo.rtfd.io/en/latest/off_target.html>`_: designs oligos to capture DNA surrounding potential CRISPR off-target cut sites to allow for efficient sequencing to determine off-target activity.

These three modules all generate oligo sequences, based on different underlying behaviours. When run from the command line, the modules utilise functions from the `tools <http://oligo.rtfd.io/en/latest/tools.html>`_ module, in a pipeline that checks
the off-target binding and repeat content of the oligos. This information is output in a file called *oligo_info.txt*; oligo sequences are written to a FASTA file called *oligo_seqs.fa*

**Requires:**

* `Python 3 <https://docs.python.org/3/>`_ (with `pysam <http://pysam.readthedocs.io/en/latest>`_, `numpy <http://www.numpy.org/>`_ & `biopython <http://biopython.org/wiki/Biopython>`_)
* `RepeatMasker <http://www.repeatmasker.org/>`_
* `STAR <https://github.com/alexdobin/STAR>`_ or `BLAT <https://genome.ucsc.edu/FAQ/FAQblat.html>`_

.. important::

    Paths to directories containing executables for STAR, BLAT and RepeatMasker must be set in the `config.txt` file before using the pipelines

For full documentation, see http://oligo.readthedocs.io

**There is currently no stable release of the oligo package**
