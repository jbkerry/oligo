Capture Oligo Design
====================

The *oligo* python package consists of three modules for Capture-C and FISH oligo design.

* [Capture](http://oligo.rtfd.io/en/latest/capture.html): designs oligos for a standard Capture-C experiment. The user supplies a list of viewpoint coordinates and oligos are generated adjacent to the flanking recognition sequence of a specified restriction enzyme.
* **Tiled**: designs oligos for multiple adjacent restriction fragments across a specified region of a chromosome, or for the entire chromosome. If *tiled* is run in FISH mode, oligos are generated independent of restriction fragments and
  are instead generated for a user-specified step size, in an end-to-end manner.
* **Off-target**: designs oligos to capture DNA surrounding potential CRISPR off-target cut sites to allow for efficient sequencing to determine off-target activity.

The three modules listed above all output a FASTA file called oligo_seqs.fa which contains the oligo sequences. When run from the command line, all of the modules listed above utilise functions from the *tools* module.

For full documentation, see http://oligo.readthedocs.io (currently incomplete)

**There is currently no stable release of the oligo package**
