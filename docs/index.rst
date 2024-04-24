.. MRChem documentation master file, created by
   sphinx-quickstart on Tue Jan 26 15:03:29 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==================================
Welcome to MRCPP's documentation!
==================================

The MultiResolution Computation Program Package (MRCPP) is a general purpose
numerical mathematics library based on multiresolution analysis and the
multiwavelet basis which provide low-scaling algorithms as well as rigorous
error control in numerical computations.

The code is being developed at the `Hylleraas Centre for Quantum Molecular
Sciences <http://www.ctcc.no/>`_ at
`UiT - The Arctic University of Norway <http://en.uit.no>`_.

--------------------------------------------------------------------------------

The code can be found on `GitHub <https://github.com/MRChemSoft/mrcpp>`_.


.. toctree::
   :maxdepth: 1
   :caption: Installation

   install.rst

.. toctree::
   :maxdepth: 1
   :caption: Application Program Interface

   mrcpp_api/introduction
   mrcpp_api/mwfunctions
   mrcpp_api/mwoperators
   mrcpp_api/gaussians
   mrcpp_api/parallel
   mrcpp_api/printer
   mrcpp_api/plotter
   mrcpp_api/timer

.. toctree::
   :maxdepth: 2
   :caption: Programmer's Manual

   programmers_manual/overview
   programmers_manual/clang_tidy.rst
