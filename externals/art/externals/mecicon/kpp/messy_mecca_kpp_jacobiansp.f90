! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Sparse Jacobian Data Structures File
!
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: GPL-3.0-only  
! ---------------------------------------------------------------

MODULE messy_mecca_kpp_JacobianSP

  PUBLIC
  SAVE


! Sparse Jacobian Data


  INTEGER, PARAMETER, DIMENSION(66) :: LU_IROW = (/ &
       1,  1,  2,  2,  2,  3,  3,  3,  4,  4,  4,  4, &
       4,  4,  5,  5,  5,  5,  5,  5,  6,  6,  6,  6, &
       6,  6,  7,  7,  7,  8,  8,  8,  8,  8,  8,  8, &
       9,  9,  9,  9,  9,  9,  9,  9,  9, 10, 10, 10, &
      10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 12, &
      12, 12, 12, 12, 12, 12 /)

  INTEGER, PARAMETER, DIMENSION(66) :: LU_ICOL = (/ &
       1, 12,  2,  6,  9,  3, 10, 11,  3,  4,  7, 10, &
      11, 12,  1,  5,  8,  9, 11, 12,  2,  6,  7,  9, &
      10, 11,  7,  9, 10,  6,  7,  8,  9, 10, 11, 12, &
       2,  5,  6,  7,  8,  9, 10, 11, 12,  3,  4,  7, &
       9, 10, 11, 12,  3,  5,  8,  9, 10, 11, 12,  1, &
       4,  7,  9, 10, 11, 12 /)

  INTEGER, PARAMETER, DIMENSION(13) :: LU_CROW = (/ &
       1,  3,  6,  9, 15, 21, 27, 30, 37, 46, 53, 60, &
      67 /)

  INTEGER, PARAMETER, DIMENSION(13) :: LU_DIAG = (/ &
       1,  3,  6, 10, 16, 22, 27, 32, 42, 50, 58, 66, &
      67 /)


END MODULE messy_mecca_kpp_JacobianSP

