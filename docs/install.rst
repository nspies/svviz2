.. _installation:

Installation
============

Prerequisites
-------------

svviz2 requires python version 3.3 or greater. Some of the packages that are installed by svviz2 also require a C compiler to be installed. 

Additional, optional software:

- To perform tandem repeat detection, download `tandem repeats finder <http://tandem.bu.edu/trf/trf.download.html>`_, rename the binary to "trf" and move it into your ``$PATH``
- To visualize the dotplots, the `rpy2 <https://rpy2.bitbucket.io>`_ package must be installed
- To convert visualizations to pdf format, either `inkscape <https://inkscape.org/>`_, `rsvg-convert <https://github.com/GNOME/librsvg>`_ or (macOS only) `webkitToPDF <https://github.com/nspies/webkitToPDF>`_ must be installed into your ``$PATH``
- While svviz2 includes bwa internally, you may wish to install `bwa mem <https://github.com/lh3/bwa>`_ in order to create reference index files


Installing into a virtual environment
-------------------------------------

The preferred method is to use `pip` to install svviz2 into a virtualenv:

::

    python3 -m venv </path/to/svviz2_venv>
    source </path/to/svviz2_venv/bin/activate>
    pip3 install -U git+git://github.com/nspies/svviz2.git

**An important note**

svviz2 compiles against the installed version of pysam. If you update pysam, you very likely will need to reinstall svviz2.