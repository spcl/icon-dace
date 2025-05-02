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

# Math-Interpolation library
This repository is an external library of ICON collecting modules of ICON that does the mathematical operations of interpolation and reconsctuction of scalar and vector fields

[![pipeline status](https://gitlab.dkrz.de/icon-libraries/libmath-interpolation/badges/master/pipeline.svg)](https://gitlab.dkrz.de/icon-libraries/libmath-interpolation/pipelines/latest?ref=master)

## Requirements
The following packages/libraries are required for `libmath-interpolation`.
- Fortran compiler
- CMake 3.18+

The following requirements are optional

- `fprettify` for Fortran code formatting
- `OpenACC` to run the code on GPU nodes

---
# For Users

## What modules are in the `libmath-interpolation` library?
The `libmath-interpolation` library includes modules that interpolates and reconstructs scalar and vector fields used in ICON. Here is a list of the supported modules.
- `mo_lib_interpolation_scalar`: contains the implementation of interpolation and reconstruction for scalar fields
- `mo_lib_interpolation_vector`: contains the implementation of interpolation and reconstruction for vector fields
- `mo_lib_intp_rbf`: perform vector RBF reconstruction

## How to link modules from `libmath-interpolation`?
- The library will be built into `libmath-interpolation/build/src/libmath-interpolation.a`
- The module files (Fortran) will be built under `libmath-interpolation/build/src/mod/`

---
# For developers

## Some notes for developers
- The `libmath-interpolation` library is only configured by CMake.
  - Tips and standards on CMake https://gitlab.dkrz.de/icon/wiki/-/wikis/CMake-recommendations-and-requirements
- The `libmath-interpolation` library will use `fprettify` in future for formatting Fortran codes, although it does not use it currently.
- The `libmath-interpolation` library is unit tested(work in progress). All merge request changes are preferable to have a unit test.
- Fortran preprocessing is automatically applied for files with `.F90` extensions.

## How to add modules in `libmath-interpolation`?
1. Put your module file under `libmath-interpolation/src`
2. Add your file to the library configuration list at `libmath-interpolation/src/CMakeLists.txt`.
```
add_library(math-interpolation
  mo_lib_interpolation_scalar.F90
  ...
# Add your module to the list. The list is in alphabetical order.
  ...
  mo_lib_intp_rbf.F90
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
- __Pradipta Samanta__ (samanta@dkrz.de)
- __Yen-Chen Chen__ (yen-chen.chen@kit.edu)
- __Sergey Kosukhin__ (sergey.kosukhin@mpimet.mpg.de)

This repository is owned by the `icon-libraries` group, contacts about general ICON library questions:
- __Terry Cojean__ (terry.cojean@kit.edu)
- __Will Sawyer__ (william.sawyer@cscs.ch)
- __Florian Prill__ (florian.prill@dwd.de)
- __Luis Kornblueh__ (luis.kornblueh@mpimet.mpg.de)
