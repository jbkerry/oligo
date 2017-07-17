#####
Tools
#####

.. container:: subtitle

    tools.py

Description
===========

.. automodule:: tools
    :platform: Unix

If you run any of the other three modules as a pipeline, i.e. from the command line, you won't need to use the `tools` module as this is incorporated into the pipelines. However, the functions in `tools` are
listed below in case you prefer to run individual parts of the pipeline from a Python shell.

Functions
=========

.. autofunction:: write_oligos
.. autofunction:: check_off_target
.. autofunction:: get_density

Private functions
=================

.. py:function:: _get_gc
.. py:function:: _get_repeats
.. py:function:: _write_file

.. centered:: :doc:`Top of Page <tools>`