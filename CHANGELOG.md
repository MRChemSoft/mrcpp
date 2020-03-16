# Change log

## Version 1.2.0-alpha3 2020-03-16

### Added

- Projecton of vector function onto FunctionTreeVector (API)
- Construction of semi-periodic Gaussians (non-API)

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


## Version 1.1.0 2019-10-30

### Added

- Option in the API for building trees using absolute MW precision.
- Functionality to estimate absolute node overlaps for screening purposes.

### Fixed

- Added missing b-spline derivative header to MWOperators include.
- Make sure that Eigen3 is found along MRCPP.
