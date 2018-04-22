------------------
Obtaining the code
------------------

The latest version of MRCPP is available on `GitHub
<http://github.com/MRChemSoft/mrcpp>`_::

    $ git clone git@github.com:MRChemSoft/mrcpp.git


-----------------
Building the code
-----------------


Prerequisites
-------------

* g++-4.8 (or equivalent)
* `CMake <http://cmake.org>`_ version 3.3 or higher.
* `Eigen <http://eigen.tuxfamily.org>`_ version 3.1 or higher.
* BLAS (optional)


Configuration
-------------

The configuration and build process is managed through CMake, and a ``setup``
script is provided for the configuration step. MRCPP's only dependency is
Eigen3, and the path to this library must be specified::

    $ export EIGEN3_ROOT=/path/to/eigen3
    $ ./setup [options] [<builddir>]

The setup script will create a directory called *<builddir>* (default ``build``)
and run CMake. There are several options available for the setup, and the most
important are:

``--cxx=<CXX>``
  C++ compiler [default: g++]
``--omp``
  Enable OpenMP parallelization [default: False]
``--mpi``
  Enable MPI parallelization [default: False]
``--enable-tests``
  Enable tests [default: True]
``--type=<TYPE>``
  Set the CMake build type (debug, release, relwithdebinfo, minsizerel) [default: release]
``--prefix=<PATH>``
  Set the install path for make install
``-h --help``
  List all options


Compilation
-----------

After successful configuration, the code is compiled using the ``make`` command
in the *<builddir>* directory::

    $ cd <builddir>
    $ make


-------------
Running tests
-------------

A set of tests is provided with the code to verify that the code compiled
properly. To compile the test suite, add the ``--enable-tests`` option to
setup, then run the tests with ``make test``::

    $ ./setup --enable-tests build
    $ cd build
    $ make
    $ make test


----------------
Running examples
----------------

In addition to the test suite, the code comes with a number of small code
snippets that demonstrate the features and the API of the library. These are
located in the ``examples`` directory. To compile the example codes, add the
``enable-examples`` option to setup. This will add a ``examples`` directory in
your *<builddir>* where the executables will be located. E.g. to compile and run
the MW projection example::

    $ ./setup --enable-examples build-serial
    $ cd build-serial
    $ make
    $ cd examples
    $ ./projection

The shared memory parallelization (OpenMP) is controlled by the environment
variable ``OMP_NUM_THREADS`` (make sure you have compiled with the ``--omp``
option to setup). E.g. to compile and run the Poisson solver example using 10
CPU cores::

    $ ./setup --enable-examples --omp build-omp
    $ cd build-omp
    $ make
    $ cd examples
    $ OMP_NUM_THREADS=10 ./poisson

To run in MPI parallel, use the ``mpirun`` (or equivalent) command (make sure
you have compiled with the ``--mpi`` option to setup, and used MPI compatible
compilers, e.g. ``--cxx=mpicxx``). Only examples with an `mpi` prefix will be
affected by running in MPI::

    $ ./setup --cxx=mpicxx --enable-examples --mpi build-mpi
    $ cd build-mpi
    $ make
    $ cd examples
    $ mpirun -np 4 ./mpi_send_tree

To run in hybrid OpenMP/MPI parallel, simply combine the two above::

    $ ./setup --cxx=mpicxx --enable-examples --omp --mpi build-hybrid
    $ cd build-hybrid
    $ make
    $ cd examples
    $ export OMP_NUM_THREADS=5
    $ mpirun -np 4 ./mpi_send_tree

Note that the core of MRCPP is *only* OpenMP parallelized. All MPI data or work
distribution must be done manually in the application program, using the tools
provided by MRCPP (see the Parallel section of the API).

----------
Pilot code
----------

Finally, MRCPP comes with a personal sandbox where you can experiment and test
new ideas, without messing around in the git repository. In the ``pilot/``
directory you will find a skeleton code called ``mrcpp.cpp.sample``. To trigger
a build, re-name (copy) this file to ``mrcpp.cpp``::

    $ cd pilot
    $ cp mrcpp.cpp.sample mrcpp.cpp

Now a corresponding executable will be build in ``<builddir>/pilot/``. Feel
free to do whatever you like in your own pilot code, but please don't add this
file to git. Also, please don't commit any changes to the existing examples
(unless you know what you're doing).

---------------------
MRCPP as a dependency
---------------------

Building MRCPP provides CMake configuration files exporting the libraries and
headers as targets to be consumed by third-party projects also using CMake::

    $ ./setup --prefix=$HOME/Software/mrcpp
    $ cd build
    $ make
    $ ctest
    $ make install

Now libraries, headers and CMake configuration files can be found under the
given prefix::

    mrcpp/
    ├── include/
    │   └── MRCPP/
    ├── lib64/
    │   ├── libmrcpp.a
    │   ├── libmrcpp.so -> libmrcpp.so.1*
    │   └── libmrcpp.so.1*
    └── share/
        └── cmake/

As an example, the ``pilot`` sample can be built with the following ``CMakeLists.txt``:

.. literalinclude:: snippets/CMakeLists.txt

This will set up the include paths and library paths correctly.
During configuration you will have to specify *where* the CMake configuration
file for MRCPP is located::

   $ cmake -H. -Bbuild -DMRCPP_DIR=$HOME/Software/share/cmake/MRCPP
