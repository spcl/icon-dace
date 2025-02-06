<!--
ICON

---------------------------------------------------------------
Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
Contact information: icon-model.org

See AUTHORS.TXT for a list of authors
See LICENSES/ for license information
SPDX-License-Identifier: CC-BY-4.0
---------------------------------------------------------------
-->

# Math-Support library
This repository is an external library of ICON collecting low-level modules of ICON supporting math operations.

[![pipeline status](https://gitlab.dkrz.de/icon-libraries/libmath-support/badges/master/pipeline.svg)](https://gitlab.dkrz.de/icon-libraries/libmath-support/pipelines/latest?ref=master)  
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://gitlab.dkrz.de/icon-libraries/libmath-support/-/blob/master/LICENSES/BSD-3-Clause.txt)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://gitlab.dkrz.de/icon-libraries/libmath-support/-/blob/master/LICENSES/CC-BY-4.0.txt)
[![License: CC0-1.0](https://img.shields.io/badge/License-CC0_1.0-lightgrey.svg)](https://gitlab.dkrz.de/icon-libraries/libmath-support/-/blob/master/LICENSES/CC0-1.0.txt)

## Requirements
The following packages/libraries are required for `libmath-support`.
- Fortran compiler
- CMake 3.18+

The following requirements are optional
- OpenACC(nvhpc) for GPU support
- [`fprettify`](https://github.com/pseewald/fprettify) for Fortran code formatting
- [REUSE](https://reuse.software) v3.0.0+ for licensing

---
# For Users

## What modules are in the `libmath-support` library?
The `libmath-support` library includes some general modules that are used in ICON to support math operations but are independent of the data types in ICON. Here is a list of the supported modules.
- `mo_gridman_constants`: constants relevant for implementational issues
- `mo_lib_grid_geometry_info`: grid-related basic structures and geometry parameters
- `mo_lib_loopindices`: determines start and end indices of do loops
- `mo_math_constants`: main mathematical constants
- `mo_math_types`: basis math types
- `mo_math_utilities`: various mathematical algorithms

## How to link modules from `libmath-support`?
- The library will be built into `libmath-support/build/src/libmath-support.so`
- The module files (Fortran) will be built under `libmath-support/build/src/mod/`

---
# For developers

## Some notes for developers
- The `libmath-support` library is only configured by CMake.
  - Tips and standards on CMake [ICON developer wiki/CMake recommendations and requirements](https://gitlab.dkrz.de/icon/wiki/-/wikis/CMake-recommendations-and-requirements)
- The `libmath-support` library uses `fprettify` for formatting Fortran codes. Run `make format` before you commit.
- The `libmath-support` library is unit tested. (work in progress) All merge request changes are preferable to have a unit test. See [icon-c/Wiki/Testing and building of ICON C/Unit test frameworks](https://gitlab.dkrz.de/icon/icon-c/-/wikis/ICON-C-Phase-0/Testing-and-building-of-ICON-C#unit-test-frameworks) for more information on unit testing.
  - __Fortran__ unit tests are written in [FortUTF](https://github.com/artemis-beta/FortUTF). [Assertions Documentation](https://github.com/artemis-beta/FortUTF/blob/main/docs/assertions.md)
- Fortran preprocessing is automatically applied for files with `.F90` extensions. See [\#4](https://gitlab.dkrz.de/icon-libraries/libfortran-support/-/issues/4) for more details.

## How to add modules in `libmath-support`?
1. Put your module file under `libmath-support/src`
2. Add your file to the library configuration list at `libmath-support/src/CMakeLists.txt`.
```
add_library(math-support
  mo_gridman_constants.f90
  mo_lib_grid_geometry_info.f90
  ...
# Add your module to the list. The list is in alphabetical order.
  mo_math_types.f90
  mo_math_utilities.F90
  ...
)
```
3. Try to compile the code.
```
mkdir build
cd build
cmake ..
! Or the following for GPU support with OpenACC
cmake -DMS_ENABLE_OPENACC=ON ..
make
```
4. Format the code by `make format`.
5. Make sure your code is tested. Check [developer note](#some-notes-for-developers).
6. Check license by `reuse lint`. Check [requirements](#requirements).
    - Code snippets should have license BSD-3-Clause
    - Documentations should have license CC-BY-4.0
    - Files unrelated to the library itself should have license CC0-1.0

## How to contribute
Please open a merge request and select one of our templates: __[feature/bugfix]__. Detailed instructions on how to proceed are provided there.

## Contact
This repository is mainly maintained by the following maintainers:
- __Pradipta Samanta__ (samanta@dkrz.de)
- __Yen-Chen Chen__ (yen-chen.chen@kit.edu)
- __Sergey Kosukhin__ (sergey.kosukhin@mpimet.mpg.de)

This repository is owned by the `icon-libraries` group, contacts about general ICON library questions:
- __Terry Cojean__ (terry.cojean@kit.edu)
- __Will Sawyer__ (william.sawyer@cscs.ch)
- __Florian Prill__ (florian.prill@dwd.de)
- __Luis Kornblueh__ (luis.kornblueh@mpimet.mpg.de)
