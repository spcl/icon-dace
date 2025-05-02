# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/io/shared/mo_cdi_pio_interface.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/build/verification//"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/io/shared/mo_cdi_pio_interface.f90"
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE mo_cdi_pio_interface
  IMPLICIT NONE
  PRIVATE
  INTEGER, PUBLIC :: cdi_base_namespace
  INTEGER, PUBLIC :: nml_io_cdi_pio_namespace
  INTEGER, PUBLIC :: nml_io_cdi_pio_client_comm
  INTEGER, PUBLIC :: nml_io_cdi_pio_conf_handle
END MODULE mo_cdi_pio_interface
