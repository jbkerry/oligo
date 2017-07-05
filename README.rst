Capture Oligo Design
====================

The *oligo* python package consists of three modules for Capture-C and FISH oligo design.

* `capture.py <http://oligo.rtfd.io/en/latest/capture.html>`_: designs oligos for a standard Capture-C experiment. The user supplies a list of viewpoint coordinates and oligos are generated adjacent to the flanking recognition sequence of a specified restriction enzyme.
* `tiled.py <http://oligo.rtfd.io/en/latest/tiled.html>`_: designs oligos for multiple adjacent restriction fragments across a specified region of a chromosome, or for the entire chromosome. If *tiled* is run in FISH mode, oligos are generated independent of restriction fragments and
  are instead generated for a user-specified step size, in an end-to-end manner.
* `off_target.py <http://oligo.rtfd.io/en/latest/off_target.html>`_: designs oligos to capture DNA surrounding potential CRISPR off-target cut sites to allow for efficient sequencing to determine off-target activity.

The three modules listed above all output a FASTA file called *oligo_seqs.fa* which contains the oligo sequences. When run from the command line, the modules utilise functions from the *tools* module in a pipeline that checks the off-target binding and repeat content of the oligos.
This information is output in a file called *oligo_info.txt*.

For full documentation, see http://oligo.readthedocs.io (currently incomplete)

**There is currently no stable release of the oligo package**
