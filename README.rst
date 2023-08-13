####################
Capture Oligo Design
####################

For full documentation, see http://oligo.readthedocs.io

``oligo`` provides functionality to automate primer design for DNA capture experiments, providing the user with details about efficiency of the primers generated.

.. contents:: Table of Contents
   :depth: 2

Installation
============
Local
-----
oligo
^^^^^
To install ``oligo`` on your local machine, it is recommended to first create a new Python environment ( >=3.8 ) using your preferred method e.g. ``conda``, ``pyenv`` etc. Once the new environment is activated, install ``oligo``
via ``pip``. Note that the package is called ``oligo`` but has a published name of ``oligo-capture`` to ensure it had a unique name in the public repository. 

.. code-block:: bash

  $ pip install oligo-capture

Ensure it has installed correctly by running the following command and verifying that you see the installed version in the standard output

.. code-block:: bash

  $ python -m oligo --version
  oligo v0.2.0

Dependencies
^^^^^^^^^^^^

Before running the full ``oligo`` pipeline you will need to install RepeatMasker and either BLAT or STAR, depending on how you intend to run ``oligo`` (see below for more details)

*RepeatMasker*
  ``oligo`` uses RepeatMasker (RM) to determine if oligos contain simple sequence repeats as these can reduce the efficiency of the oligo for targeted capture. Follow the instructions
  on `RepeatMasker home page <http://www.repeatmasker.org/>`_ for installing RM on your local system. The most recent version of RM that ``oligo`` has been tested with is v4.1.5. As detailed
  on the RM home page, its installation depends on a Sequence Search Engine (it is recommended to use **HMMER** for ``oligo``) and Tandem Repeat Finder (TRF). For the Repeat Database, RM ships
  with the curated set of Dfam 3.7 Transposable Elements which is sufficient but users are free to use the full set if required; further instructions are on the RM home page.

.. highlights::

  The RM home page mentions that it requires the Python library ``h5py``, however this is listed as a dependency of the ``oligo`` package so will already be installed in your Python environment
  from when you ran the ``pip install`` step.

*BLAT*
  ``oligo`` uses the BLAST-Like Alignment Tool (BLAT) to determine any off-target binding sites of an oligo within the genome, in addition to its intended binding site. An oligo that binds
  to multiple regions will have a reduced score since it will perform a less-specific capture. `BLAT <https://genome.ucsc.edu/FAQ/FAQblat.html>`_ executables can be found by going to
  `<http://hgdownload.soe.ucsc.edu/admin/exe/>`_ and locating the BLAT directory in the for your systems archtecture. For example, for Linux.x86 architecture, ``rsync`` should be used
  to get the BLAT executables on your system:

.. code-block:: bash

  $ rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/blat/ ./

See `<http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/>`_ for more details.

*STAR*
  As an alternative to BLAT, ``oligo`` allows users to use the Spliced Transcripts Alignment to a Reference (STAR) alignment program for increased speed when determining multiple binding
  events for oligo sequences. BLAT is more widely used to detect off-target binding events, however BLAT can be particulary slow for large designs, especially for the human
  reference genomes. STAR’s exceptional speed is better suited for designs with >1000 oligos. If you think you would prefer to use STAR instead, visit the
  `STAR GitHub page <https://github.com/alexdobin/STAR>`_ for instructions on how to install it.

Docker
------

Due to ``oligo`` requiring various third-party software, it can instead be run from a pre-made Docker image that has everything needed already installed. This should make the setup much
easier for users as well as reducing the need to install lots of software on their local machines. Running via Docker is obviously less flexible in terms of the configuration of the
third-party software but has been built with the most common use cases in mind and reducing the image size to as small as possible, without losing any of requirements ``oligo`` uses from
the third-party software.

First pull the latest oligo image onto your local machine:

.. code-block:: bash

  $ docker pull jbkerry/oligo:0.2.0-alpha

The docker entrypoint is set to run ``oligo`` with the config file already set up to point to the install executables of BLAT and RepeatMasker so users can run the image, starting with
the ``oligo`` sub-command that is required.

#TODO: mounting directories

The example command used above is shown again below but this time, using the Docker image:

.. code-block:: bash

  $ docker run -v /local/path/oligo_results:/results -v /local/human:/genome jbkerry/oligo:0.2.0-alpha off-target -f /genome/genome.fa -g hg38

Installation specifics
^^^^^^^^^^^^^^^^^^^^^^
Below is a list of the versions and alterations that have been made to the standard installs of third-party software for the ``oligo`` Docker image:
  * RepeatMasker v4.1.5
  
    * Dfam.h5 library has been replace with an HMM matrix containing only mouse- and human-specific transposable elements* in order to reduce the size of the Docker image
  * HMMER v3.3.2
  * Tandem Repeat Finder v4.09.1
  * BLAT v37.x1
  
These HMM matrices were generated with the following two commands (``famdb.py`` comes bundled with the latest versions of RepeatMasker):

.. code-block:: bash

  $ ./famdb.py -i Libraries/RepeatMaskerLib.h5 families --format hmm 'Homo sapiens' --include-class-in-name >humans.hmm
  $ ./famdb.py -i Libraries/RepeatMaskerLib.h5 families --format hmm 'Mus musculus' --include-class-in-name >mouse.hmm

The Dockerfile in the ``oligo`` GitHub repository can be referenced for details of the how the Docker imges was built. Some reference data files that get copied into the image at build
time are not present in the repository but can be provided to the user if needed.

Usage
=====

``oligo`` can be run with one of three subcommands

* `capture <http://oligo.rtfd.io/en/latest/capture.html>`_: designs oligos for a standard Capture-C experiment. The user supplies a list of viewpoint coordinates, and oligos are generated adjacent to the flanking recognition sequence of a specified restriction enzyme.
* `tiled <http://oligo.rtfd.io/en/latest/tiled.html>`_: designs oligos for multiple adjacent restriction fragments across a specified region of a chromosome, or for the entire chromosome. If ``tiled`` is run in contiguous mode, oligos are generated independent of restriction fragments and
  are instead generated for a user-specified step size, in an adjacent manner.
* `off-target <http://oligo.rtfd.io/en/latest/off_target.html>`_: designs oligos to capture DNA surrounding potential CRISPR off-target cut sites to allow for efficient sequencing to determine off-target activity.

These subcommands all generate oligo sequences, based on different underlying behaviours. Methods from the `Tools <http://oligo.rtfd.io/en/latest/tools_class.html>`_ class in the ``oligo.tools`` module are then used to check
the off-target binding and repeat content of the oligos. This information is output in a file called *oligo_info.txt*; oligo sequences are written to a FASTA file called *oligo_seqs.fa*

**Example**

The subcommand follows the ``oligo`` command and options for the subcommand are then specified afterwards.
Below, is an example using the ``off-target`` subcommand:

.. code-block:: bash

  $ python -m oligo off-target -f /path/to/human/genome.fa -g hg38 
