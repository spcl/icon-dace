!
! Basic class used in turbulent mixing package (tmx)
!
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

MODULE mo_surrogate_class

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_surrogate

  TYPE, ABSTRACT :: t_surrogate
  END TYPE t_surrogate

END MODULE mo_surrogate_class
