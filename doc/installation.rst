============
Installation
============


-------------------
Build prerequisites
-------------------

* g++-4.8 (or equivalent)
* `CMake-2.8 <http://cmake.org>`_
* `Eigen-3.1 <http://eigen.tuxfamily.org>`_


-------------------------------
Obtaining and building the code
-------------------------------

The latest version of MRCPP is available on `GitHub
<http://github.com/MRChemSoft/mrcpp>`_::

    $ git clone git@github.com:MRChemSoft/mrcpp.git

A test suite is provided (requires the ``--enable-tests`` flag to setup) to make
sure that everything is compiled properly. To build the code using OpenMP
parallelization with the GNU tool chain::

    $ cd mrcpp
    $ ./setup --omp --enable-tests
    $ cd build
    $ make
    $ make test

In order to get MPI functionalities, you need to specify a compiler with MPI
options::

    $ ./setup --omp --mpi --cxx=mpicxx


-------------------------
Running the test programs
-------------------------

In addition to the test suite, the code comes with a number of simple examples
that demonstrate the features and the API of the library.

The shared memory parallelization (OpenMP) is controlled by the environment
variable ``OMP_NUM_THREADS`` (make sure you have compiled with the ``--omp``
option to setup), and to run on 10 CPUs::

    $ export OMP_NUM_THREADS 10
    $ ./mrcpp

To run in MPI parallel, use the ``mpirun`` (or equivalent) command (make sure
you have compiled with the ``--mpi`` option to setup, and used MPI compatible
compilers, e.g. ``--cxx=mpicxx``)::

    $ mpirun -np 10 ./mrcpp

To run in hybrid OpenMP/MPI parallel, simply combine the two above::

    $ export OMP_NUM_THREADS 10
    $ mpirun -np 10 ./mrcpp

Note that the core of MRCPP is *only* OpenMP parallelized. All MPI data or work
distribution must be done manually in the application program, using the tools
provided by MRCPP.

