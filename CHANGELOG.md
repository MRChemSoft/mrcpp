# Change log

## Version 1.3.3 2020-09-23

### Fixed

- Export compiler flags with MRCPP Cmake target


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
