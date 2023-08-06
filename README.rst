####################
Capture Oligo Design
####################

The *design.py* module consists of three classes for Capture-C oligo design.

* `Capture <http://oligo.rtfd.io/en/latest/capture.html>`_: designs oligos for a standard Capture-C experiment. The user supplies a list of viewpoint coordinates, and oligos are generated adjacent to the flanking recognition sequence of a specified restriction enzyme.
* `Tiled <http://oligo.rtfd.io/en/latest/tiled.html>`_: designs oligos for multiple adjacent restriction fragments across a specified region of a chromosome, or for the entire chromosome. If *Tiled* is run in contiguous mode, oligos are generated independent of restriction fragments and
  are instead generated for a user-specified step size, in an adjacent manner.
* `OffTarget <http://oligo.rtfd.io/en/latest/off_target.html>`_: designs oligos to capture DNA surrounding potential CRISPR off-target cut sites to allow for efficient sequencing to determine off-target activity.

These three classes all generate oligo sequences, based on different underlying behaviours. When run from the command line, *design.py* uses methods from the `Tools <http://oligo.rtfd.io/en/latest/tools_class.html>`_ class in the *tools.py* module, in a pipeline that checks
the off-target binding and repeat content of the oligos. This information is output in a file called *oligo_info.txt*; oligo sequences are written to a FASTA file called *oligo_seqs.fa*

Installation
============
Local
-----

To install ``oligo`` on your local machine, it is recommended to first create a new Python environment ( >=3.8 ) using your preferred method e.g. ``conda``, ``pyenv`` etc. Once the new environment is activated, install ``oligo``
via ``pip``. Note that the package is called ``oligo`` but has a published name of ``oligo-capture`` to ensure it had a unique name in the public repository. 

.. code-block:: bash

  $ pip install oligo-capture

Ensure it has installed correctly by running the following command and verifying that you see the installed version in the standard output

.. code-block:: bash

  $ python -m oligo --version
  oligo v0.2.0

``oligo`` can be run with one of three subcommands

+----------------+------------+
| ``capture``    | Header 2   |
+----------------+------------+
| ``tiled``      | column 2   |
+----------------+------------+
| ``off-target`` | ...        |
+----------------+------------+

The subcommand follows the ``oligo`` command and options for the subcommand are then specified afterwards.
Below, is an example using the ``off-target`` subcommand:

.. code-block:: bash

  $ python -m oligo off-target -f /path/to/human/genome.fa -g hg38 

Docker
------

put here

**Requires:**

* `Python <https://docs.python.org/3/>`_ 2.7 or >=3.4  (with `pysam <http://pysam.readthedocs.io/en/latest>`_, `numpy <http://www.numpy.org/>`_, `pandas <http://pandas.pydata.org/>`_ & `biopython <http://biopython.org/wiki/Biopython>`_)
* `RepeatMasker <http://www.repeatmasker.org/>`_
* `STAR <https://github.com/alexdobin/STAR>`_ or `BLAT <https://genome.ucsc.edu/FAQ/FAQblat.html>`_

For full documentation, see http://oligo.readthedocs.io
