# Change log

## Version 1.4.2 2023-01-20

### Fixed

- Invalid read in coef block for nodes without coefs
- Node count error for negative scales with PBC
- Missing includes resulting in compile errors for recent gcc


## Version 1.4.1 2022-01-05

### Changed

- Catch version 2.13.5

### Fixed

- Buffer overflow in NodeAllocator
- Polynomial integration bound error message


## Version 1.4.0 2021-10-13

### Added

- New constructors for BoundingBox and MRA, accepting box=[-L,L] argument
- Add FunctionTree::evalf_precise() which evaluates both scaling and wavelet parts
- Possibility for empty tree skeletons without allocated coefficients
- Possibility to manually set location of MW filters at configure time

### Changed

- Eigen version 3.4.0
- Improvements under the hood of the Tree classes (non-API):
  - remove GenNode and ProjectedNode
  - simplified NodeAllocator class
  - new improved OMP locking

### Fixed

- Get node center and (lower/upper) bounds


## Version 1.3.6 2020-10-27

### Fixed

- Put OpenMP and MPI in CMake INTERFACE
- Fix --arch-flags CMake option
- Fix kramdown security issue
- Document variuos packaging options

## Version 1.3.5 2020-10-05

### Fixed

- Remove final qualifier GaussFunc/GaussPoly (fixing pybind11 bindings)

## Version 1.3.4 2020-10-02

### Fixed

- Allow to remove arch-specific optimization flags
- Undo export of compiler flags with MRCPP CMake target
- Encapsulate OpenMP within MRCPP and export MRCPP_HAS_OMP flag
- Encapsulate MPI within MRCPP and export MRCPP_HAS_MPI flag
- Fix faulty OpenMP locks
- Minor inlining optimizations

## Version 1.3.3 2020-09-23

### Fixed

- Export compiler flags with MRCPP CMake target


## Version 1.3.2 2020-09-14

### Fixed

- Miscellaneous fixes for building on conda-forge


## Version 1.3.1 2020-09-02

### Fixed

- Uninitialized (max)norms used in space-varying precision algorithms


## Version 1.3.0 2020-08-28

### Added

- Function mapping functionality
- Possibility for faster operator application by space-varying precision


## Version 1.2.0 2020-04-13

### Added

- Projecton of vector function onto FunctionTreeVector (API)
- Construction of semi-periodic Gaussians (non-API)
- CircleCI test build config

### Changed

- New versioning scheme and branching model, documented in README
- Much improved API documentation
- Proper installation of MW filters
- FunctionTrees are now evaluated as zero outside the domain

### Fixed

- Faulty cube plot grid generation
- Pull Eigen from gitlab instead of deprecated github mirror
- Incompatible compiler flags between debug and release
- Failing periodic unit tests
- Plotting range verification


## Version 1.1.0 2019-10-30

### Added

- Option in the API for building trees using absolute MW precision.
- Functionality to estimate absolute node overlaps for screening purposes.

### Fixed

- Added missing b-spline derivative header to MWOperators include.
- Make sure that Eigen3 is found along MRCPP.
