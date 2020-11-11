![MRCPP logo](https://github.com/MRChemSoft/mrcpp/raw/master/docs/gfx/logo.png)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3606670.svg)](https://doi.org/10.5281/zenodo.3606670)
[![License](https://img.shields.io/badge/license-%20LGPLv3-blue.svg)](../master/LICENSE)
[![Documentation Status](https://readthedocs.org/projects/mrcpp/badge/?version=latest)](http://mrcpp.readthedocs.io/en/latest/?badge=latest)
![Build and test MRCPP](https://github.com/MRChemSoft/mrcpp/workflows/Build%20and%20test%20MRCPP/badge.svg)
[![CircleCI](https://circleci.com/gh/MRChemSoft/mrcpp/tree/master.svg?style=svg)](https://circleci.com/gh/MRChemSoft/mrcpp)
[![codecov](https://codecov.io/gh/MRChemSoft/mrcpp/branch/master/graph/badge.svg)](https://codecov.io/gh/MRChemSoft/mrcpp)

The MultiResolution Computation Program Package (MRCPP) is a general
purpose numerical mathematics library based on multiresolution analysis
and the multiwavelet basis which provide low-scaling algorithms as well
as rigorous error control in numerical computations.

The code is being developed at the Hylleraas Centre for Quantum Molecular
Sciences at UiT - The Arctic University of Norway.

### User support: [mrchem.slack.com](https://join.slack.com/t/mrchem/shared_invite/enQtNTI3MjMzNjM0NTk0LWNkODZjNTMwYmM4NmRmODExMjQzMDc3NThlMzNmNmIyNWQwM2YwOGY0OWY4NmNmNzE4ZmM2NzgxYzUzNDg3NDM)
### Documentation: [mrcpp.readthedocs.io](http://mrcpp.readthedocs.io)


## Installation

For optimal performance it is recommended to build from source, as the packaged
builds are quite generic without architecture specific optimizations.

### From source including code examples

To build MRCPP from source with MPI+OpenMP parallelization:

    $ git clone git@github.com:MRChemSoft/mrcpp.git
    $ cd mrcpp
    $ ./setup --prefix=<install-dir> --enable-examples --omp --mpi --cxx=<mpi-compiler> <build-dir>
    $ cd <build-dir>
    $ make
    $ make test
    $ make install

The `--enable-examples` will trigger compilation of the code snippets in the
`examples/` directory. Each example will get a separate executable under
`<build-dir>/bin`, but these are not installed.

For more information on different kinds of builds, see
[installation instructions](http://mrcpp.readthedocs.io/en/latest/install.html).


### Using Conda

[![Anaconda-Server Badge](https://anaconda.org/conda-forge/mrcpp/badges/version.svg)](https://anaconda.org/conda-forge/mrcpp)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/mrcpp/badges/latest_release_date.svg)](https://anaconda.org/conda-forge/mrcpp)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/mrcpp/badges/downloads.svg)](https://anaconda.org/conda-forge/mrcpp)

To install MRCPP in a Conda environment `myenv`:

    $ conda create -n myenv
    $ conda activate myenv
    $ conda install -c conda-forge mrcpp                # latest version (OpenMP only)
    $ conda install -c conda-forge mrcpp=1.3.6          # tagged version (OpenMP only)
    $ conda install -c conda-forge mrcpp=*=*openmpi*    # latest version (MPI+OpenMP)
    $ conda install -c conda-forge mrcpp=*=*mpich*      # latest version (MPI+OpenMP)

To list all available versions

    $ conda search -c conda-forge mrcpp

### Using Spack

To install MRCPP in a Spack environment `myenv`:

    $ spack env create myenv
    $ spack env activate myenv
    $ spack install mrcpp                               # latest version (MPI+OpenMP)
    $ spack install mrcpp @1.3.6                        # tagged version (MPI+OpenMP)
    $ spack install mrcpp -mpi                          # latest version (OpenMP only)

For information on available Spack builds:

    $ spack info mrcpp


### Using EasyBuild

To install MRCPP in an EasyBuild/Lmod environment (only MPI+OpenMP version
available):

    $ eb MRCPP-<version>-<toolchain> --fetch
    $ eb MRCPP-<version>-<toolchain> --robot
    $ module load MRCPP/<version>-<toolchain>

See
[EasyBuild](https://github.com/easybuilders/easybuild-easyconfigs/tree/develop/easybuild/easyconfigs/m/MRCPP)
for available `<versions>` and `<toolchains>`.
