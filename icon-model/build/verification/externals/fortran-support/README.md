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

# Fortran-support library
This repository is an external library of ICON collecting low-level supporting modules of ICON.

[![Latest Release](https://gitlab.dkrz.de/icon-libraries/libfortran-support/-/badges/release.svg)](https://gitlab.dkrz.de/icon-libraries/libfortran-support/-/releases)
[![pipeline status](https://gitlab.dkrz.de/icon-libraries/libfortran-support/badges/master/pipeline.svg?key_text=Pipeline&key_width=55)](https://gitlab.dkrz.de/icon-libraries/libfortran-support/pipelines/latest?ref=master)
[![coverage report](https://gitlab.dkrz.de/icon-libraries/libfortran-support/badges/master/coverage.svg?key_text=Test%20coverage&key_width=90)](https://gitlab.dkrz.de/icon-libraries/libfortran-support/pipelines/latest?ref=master)  
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://gitlab.dkrz.de/icon-libraries/libfortran-support/-/blob/master/LICENSES/BSD-3-Clause.txt)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://gitlab.dkrz.de/icon-libraries/libfortran-support/-/blob/master/LICENSES/CC-BY-4.0.txt)
[![License: CC0-1.0](https://img.shields.io/badge/License-CC0_1.0-lightgrey.svg)](https://gitlab.dkrz.de/icon-libraries/libfortran-support/-/blob/master/LICENSES/CC0-1.0.txt)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://gitlab.dkrz.de/icon-libraries/libfortran-support/-/blob/master/LICENSES/MIT.txt)

## Requirements
The following packages/libraries are required for `libfortran-support`.
- Fortran compiler
- C compiler
- CMake 3.18+

The following requirements are optional
- C++ compiler (testing with [GoogleTest](https://github.com/google/googletest))
- OpenACC(nvhpc) for GPU support
- [`fprettify`](https://github.com/pseewald/fprettify) for Fortran code formatting
- `clang-format` for C/C++ code formatting
- [REUSE](https://reuse.software) v3.0.0+ for licensing
- [Ragel State Machine Compiler](http://www.colm.net/open-source/ragel/) 7.0+ for C code generation from `.rl` files, used for
  - `nml_annotate.c`
  - `util_arithmetic_expr.c`
  - `util_string_parse.c`

---
# For Users

## What modules are in the `libfortran-support` library?
The `libfortran-support` library includes some general Fortran supporting modules that are used in ICON but are independent of the data types in ICON. Here is a list of the supported modules.
- `mo_exception`: message logger for ICON
- `mo_expression`: expression parsing
- `mo_fortran_tools`: basic array allocation, copy, etc. subroutines with GPU(OpenACC) support
- `mo_hash_table`: hash table operations
- `mo_io_units`: io unit definitions
- `mo_namelist`: open/close namelist files
- `mo_octree`: octree data structure and operations
- `mo_random_number_generators`: generators for pseudo random numbers (should be moved to math-support once available)
- `mo_simple_dump`: array value dumping
- `mo_util_backtrace`: function backtrace
- `mo_util_file`: file operations
- `mo_util_libc`: standard C functions interface
- `mo_util_nml`: read/annotate namelist files
- `mo_util_rusage`: RSS information list
- `mo_util_sort`: array sorting
- `mo_util_stride`: data stride
- `mo_util_string`: string handling
- `mo_util_string_parse`: string parsing for arithmetic expressions
- `mo_util_system`: system functions (exit, abort, ...)
- `mo_util_table`: build and use table
- `mo_util_texthash`: hash text operations
- `mo_util_timer`: timer functions

## How to link modules from `fortran-support`?
- The library will be built into `fortran-support/build/src/libfortran-support.so`
- The module files (Fortran) will be built under `fortran-support/build/src/mod/`

---
# For developers

## Some notes for developers
- The `fortran-support` library is only configured by CMake.
  - Tips and standards on CMake [ICON developer wiki/CMake recommendations and requirements](https://gitlab.dkrz.de/icon/wiki/-/wikis/CMake-recommendations-and-requirements)
- The `fortran-support` library uses `fprettify` for formatting Fortran codes. Run `make format` before you commit.
- The `fortran-support` library is unit tested. All merge request changes are required to have a unit test. See [icon-c/Wiki/Testing and building of ICON C/Unit test frameworks](https://gitlab.dkrz.de/icon/icon-c/-/wikis/ICON-C-Phase-0/Testing-and-building-of-ICON-C#unit-test-frameworks) for more information on unit testing.
  - __Fortran__ unit tests are written in [FortUTF](https://github.com/artemis-beta/FortUTF). [Assertions Documentation](https://github.com/artemis-beta/FortUTF/blob/main/docs/assertions.md)
  - __C__ unit tests are written in [GoogleTest](https://github.com/google/googletest). [User's Guide](https://google.github.io/googletest/)
- Fortran preprocessing is automatically applied for files with `.F90` extensions. See [\#4](https://gitlab.dkrz.de/icon-libraries/libfortran-support/-/issues/4) for more details.

## How to add modules in `fortran-support`?
1. Put your module file under `fortran-support/src`
2. Add your file to the library configuration list at `fortran-support/src/CMakeLists.txt`.
```
add_library(fortran-support
  mo_exception.F90
  mo_expression.F90
  ...
# Add your module to the list. The list is in alphabetical order.
  mo_util_nml.F90
  mo_util_rusage.F90
  mo_util_sort.F90
  ...
)
```
3. Try to compile the code.
```
mkdir build
cd build
cmake ..
! Or the following for GPU support with OpenACC
cmake -DFS_ENABLE_OPENACC=ON ..
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
- __Yen-Chen Chen__ (yen-chen.chen@kit.edu)
- __Jonas Jucker__ (jonas.jucker@env.ethz.ch)
- __Will Sawyer__ (william.sawyer@cscs.ch)

This repository is owned by the `icon-libraries` group, contacts about general ICON library questions:
- __Terry Cojean__ (terry.cojean@kit.edu)
- __Will Sawyer__ (william.sawyer@cscs.ch)
- __Florian Prill__ (florian.prill@dwd.de)
- __Luis Kornblueh__ (luis.kornblueh@mpimet.mpg.de)
