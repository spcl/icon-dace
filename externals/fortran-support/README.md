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

## Requirements
The following packages/libraries are required for `libfortran-support`.
- Fortran compiler
- C compiler
- C++ compiler
- CMake 3.18+

The following requirements are optional

- CMake 3.20+ for Fortran unit testing
- `fprettify` for Fortran code formatting
- `clang-format` for C/C++ code formatting
- Ragel State Machine Compiler 7.0+ for C code generation from `.rl` files, used for
  - `nml_annotate.c`
  - `util_arithmetic_expr.c`
  - `util_string_parse.c`

---
# For Users

## What modules are in the `libfortran-support` library?
The `libfortran-support` library includes some general Fortran supporting modules that are used in ICON but are independent of the data types in ICON. Here is a list of the supported modules.
- `mo_exception`: message logger for ICON
- `mo_expression`: expression parsing
- `mo_hash_table`: hash table operations
- `mo_io_units`: io unit definitions
- `mo_namelist`: open/close namelist files
- `mo_octree`: octree data structure and operations
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
  - Tips and standards on CMake https://gitlab.dkrz.de/icon/wiki/-/wikis/CMake-recommendations-and-requirements
- The `fortran-support` library uses `fprettify` for formatting Fortran codes. Run `make format` before you commit.
- The `fortran-support` library is unit tested. (work in progress) All merge request changes are preferable to have a unit test.
- Fortran preprocessing is automatically applied for files with `.F90` extensions. See [\#4](https://gitlab.dkrz.de/icon-libraries/libfortran-support/-/issues/4) for more details.

## How to add modules in `fortran-support`?
1. Put your module file under `fortran-support/src`
2. Add your file to the library configuration list at `fortran-support/src/CMakeLists.txt`.
```
add_library(fortran-support
  mo_exception.F90
  mo_expression.F90
  ...
# Add your module to the list. The list is in alphebatic order.
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
make
```
4. Format the code by `make format`.
5. Make sure your code is tested. For more information on unit tests, check out this [GitLab WIKI page](https://gitlab.dkrz.de/icon/icon-c/-/wikis/ICON-C-Phase-0/Testing-and-building-of-ICON-C#unit-test-frameworks)

## How to contribute
Please open a merge request and select one of our templates for new features or bugfixes. Detailed instructions on how to proceed are provided there.

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
