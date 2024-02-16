!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
!>
!!  Module determines kinds for different precisions.
!!
!!  Module determines kinds for different precisions
!!  Number model from which the SELECTED_*\\_KIND are requested: <br>
!!
!! @f{tabular}{{r@{\hspace*{3em}}c@{\hspace*{3em}}c}
!!                     &4 byte REAL     &8 byte REAL        \\\
!!        CRAY:        &-               &precision =   13   \\\
!!                     &                &exponent  = 2465   \\\
!!        IEEE:        &precision = 6   &precision =   15   \\\
!!                     &exponent  = 37  &exponent  =  307
!! @f}
!! \\medskip
!!
!!  Most likely this are the only possible models.
!!
!! @par Revision History
!!  Working precision and comments are added by Luis Kornblueh (2001)
!!  Modification by Luis Kornblueh (2010-02-16):
!!  Working precision selection is parameterized by values, which can be as
!!  well used for selection of MPI data types associated with the respective
!!  Fortran datatypes.
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
!!
MODULE mo_kind

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  !--------------------------------------------------------------------
  !
  ! Floating point section
  ! ----------------------
  !
  INTEGER, PARAMETER :: ps =   6
  INTEGER, PARAMETER :: rs =  37
  !
  INTEGER, PARAMETER :: pd =  12
  INTEGER, PARAMETER :: rd = 307
  !
  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(ps,rs) !< single precision
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(pd,rd) !< double precision
  !
  INTEGER, PARAMETER :: wp = dp                        !< selected working precision
  !
#ifdef __MIXED_PRECISION
  INTEGER, PARAMETER :: vp = sp
#else
  INTEGER, PARAMETER :: vp = wp
#endif
  !
  ! Integer section
  ! ---------------
  !
  INTEGER, PARAMETER :: pi2 =  4
  INTEGER, PARAMETER :: pi4 =  9
  INTEGER, PARAMETER :: pi8 = 14  ! could be larger, but SX cannot do some operations otherwise
  !
  INTEGER, PARAMETER :: i2 = SELECTED_INT_KIND(pi2)   !< at least 2 byte integer
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(pi4)   !< at least 4 byte integer
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(pi8)   !< at least 8 byte integer
  !
  !
  ! The following variable is made available internally only. configure needs to detect
  ! the addressing mode and according to this mo_kind has to be updated by an preprocessor
  ! directive and #include <config.h>. This needs some changes.
  !
  INTEGER, PARAMETER :: wi = i4                       !< selected working precission
  !
  PUBLIC :: sp, dp, wp, vp, i4, i8, i2
  !
  !--------------------------------------------------------------------

END MODULE mo_kind

