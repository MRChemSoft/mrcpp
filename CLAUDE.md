# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

MRCPP is a C++ multiresolution/multiwavelet numerical library (C++ library, used by MRChem). Build system is CMake (Autocmake-generated, with CPM for fetching Catch2). Optional MPI and OpenMP.

## Build and test

Configure with the `setup` wrapper script (a CMake front-end), then build in the generated directory:

    ./setup --prefix=<install-dir> [--omp] [--mpi] [--cxx=<compiler>] [--enable-examples] <build-dir>
    cd <build-dir> && make && make test

- Tests are Catch2 v3, compiled into a single executable `mrcpp-tests` (under `<build-dir>/bin`). `ENABLE_TESTS` is ON by default.
- Run a single test case or tag directly: `<build-dir>/bin/mrcpp-tests "[mw_tree]"` or `mrcpp-tests "MWTree: Fetching nodes"`. Tests are tagged in `TEST_CASE("name", "[specific_tag], [group_tag], [dir]")`; CTest registers one test per tag (see `add_Catch_test` in `cmake/custom/tests.cmake`), so `ctest -R <tag>` also works.
- `--enable-examples` builds each file in `examples/` as its own executable in `<build-dir>/bin` (not installed).
- Several build dirs exist locally (`_build`, `build`, `build_dev`, `build-debug`, ...) and `setup_dev.sh`; these are untracked scratch — don't commit them.
- Formatting: `.clang-format` at the repo root (4-space indent). Docs use Doxygen (`docs/Doxyfile.in`) + Sphinx/RTD (`docs/`, `docs/programmers_manual/*.rst`).

## Architecture

Source in `src/` is split by layer; the public headers are re-exported through `api/` (`api/MWFunctions`, `MWOperators`, `Gaussians`, `Plotter`, `Printer`, `Timer`, `Parallel`) and `api/mrcpp_declarations.h` holds forward declarations. When adding a public class, update the corresponding `api/` header and the relevant `CMakeLists.txt` in `src/<dir>` (sources are listed explicitly).

- `src/core` — numerical foundations: scaling bases (`LegendreBasis`, `InterpolatingBasis`), multiwavelet filters (`MWFilter`, `FilterCache`), Gauss quadrature, cross-correlation coefficient tables (`CrossCorrelation*`, `SchrodingerEvolution_CrossCorrelation`). Singleton-style caches derive from `ObjectCache`; filter/cross-correlation data files come from `share/` (see `cmake/custom/mw_filters.cmake`).
- `src/trees` — the data structures. `MWTree` / `MWNode` are the base; `FunctionTree`/`FunctionNode` hold function coefficients, `OperatorTree`/`OperatorNode` (and `CornerOperatorTree`) hold operator matrices. `MultiResolutionAnalysis` bundles a `BoundingBox` (root-node layout) with a scaling basis and max scale; trees are always built against an MRA. `NodeAllocator` manages node/coefficient memory pools, `NodeBox`/`HilbertPath`/`NodeIndex` handle spatial indexing, `TreeIterator` traverses.
- `src/functions` — analytic functions (`RepresentableFunction`, `GaussFunc`, `GaussExp`, `GaussPoly`, `Polynomial`, `LegendrePoly`, `BoysFunction`, special functions) used as projection sources.
- `src/treebuilders` — algorithms that build/modify trees using a **Calculator + Adaptor + TreeBuilder** pattern: `TreeBuilder` drives refinement; a `*Calculator` computes node coefficients (`ProjectionCalculator`, `ConvolutionCalculator`, `AdditionCalculator`, `MultiplicationCalculator`, ...); a `*Adaptor` decides which nodes to split (`WaveletAdaptor`, `AnalyticAdaptor`, `BandWidthAdaptor`, ...). Free-function entry points live in `project.cpp`, `apply.cpp`, `complex_apply.cpp`, `add.cpp`, `multiply.cpp`, `grid.cpp`, `map.cpp`.
- `src/operators` — `MWOperator` subclasses (`Poisson`, `Helmholtz`, `Heat`, `TimeEvolution`, `ABGV`/`PH`/`BS` derivatives, `Identity`). Convolution operators are built by projecting a kernel (`*Kernel`) with `CrossCorrelationCalculator` into an `OperatorTree`, then applied via `apply()`.
- `src/utils` — parallel helpers (MPI/OpenMP), `Printer`, `Timer`, `Plotter`-support, math utils, tree-wide diagnostics.

Tests in `tests/` mirror the `src/` layout (`core`, `functions`, `operators`, `trees`, `treebuilders`, `scalingfactor`), each with its own `CMakeLists.txt` that must list new test files; shared helpers are in `tests/factory_functions.h`.

## Conventions / workflow

- **Coding standards for new code are in `CODING_STANDARDS.md`** (naming and brace/indent style). Follow them for all new or modified lines; do not reformat or rename legacy code unprompted. Background survey: `docs/CODING_STANDARDS_AUDIT.md`.

- Branching model (see `CONTRIBUTING.md`): `master` is the development line and carries a `-alpha` version label in `VERSION`; releases are on `release/X.Y` branches; features target `master`, PATCH fixes are applied linearly on release branches. Semantic versioning.
- Source files carry the LGPL license header block at the top.
- The current branch `coding_standards` is about code-style/documentation cleanup (Doxygen comments, `docs/programmers_manual`, clang-tidy notes in `docs/programmers_manual/clang_tidy.rst`).
