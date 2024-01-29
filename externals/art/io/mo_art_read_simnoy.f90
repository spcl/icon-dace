!
! mo_art_read_simnoy
! This module sets initial values simplified N2O/NOy chemistry
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

MODULE mo_art_read_simnoy
     USE mo_kind,                 ONLY: wp
  
  IMPLICIT NONE
  
  PRIVATE
  
  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_read_simnoy'

! ---------------------------------
! define local table parameter
! ---------------------------------

  INTEGER, PARAMETER    :: lparm = 20        ! heights
  INTEGER, PARAMETER    :: kparm = 18        ! latitudes
  INTEGER, PARAMETER    :: mparm = 12        ! month
  INTEGER, PARAMETER    :: nparm = 5         ! tables
  REAL(wp), PARAMETER   :: lparm_min = 14.   ! <lowest height with available
                                           !  values in table, (km)
  REAL(wp), PARAMETER   :: lparm_dist = 2.   ! <constant distance between
                                              !  table height levels

  REAL(wp), DIMENSION(lparm,kparm,mparm,nparm)   :: tparm_noy
  
  
  PUBLIC :: n2onoy_read, tparm_noy
  PUBLIC :: lparm, kparm, mparm, nparm, lparm_min, lparm_dist

  CONTAINS
!!
!!-------------------------------------------------------------------------
!!
! -----------------------------------------
! --- READ LINNOY tables from ... 
! --- for simplified N2O/NOy chemistry
! -----------------------------------------

SUBROUTINE n2onoy_read()

!<
! SUBROUTINE n2onoy_read
! reads values needed for initializing 
! simplified N2O/NOy chemistry
! based on Olsen et al, 2001
! tables from JPL 98 (?)
! Part of module: mo_art_read_simnoy
! Author: Christopher Diekmann
!>

  CHARACTER(LEN=50)     :: title
  CHARACTER(LEN=120)    :: n2onoyinitfile
  INTEGER               :: n,m,k,l

  n2onoyinitfile='Simnoy2002.dat'

  IF(n2onoyinitfile == '') THEN
      tparm_noy(:,:,:,:) = 0.
  ELSE
    OPEN(61, FILE=n2onoyinitfile, STATUS='old')
    READ(61,901) title
    
    DO n=1,nparm
      READ(61,901) title
      DO m=1,mparm
        DO k=1,kparm
          READ(61,902) (tparm_noy(l,k,m,n),l=1,lparm)
        ENDDO
      ENDDO
    ENDDO

    901 FORMAT(A50)
    902 FORMAT(20x,6e10.3/(8e10.3))
    CLOSE(61)
  ENDIF

END SUBROUTINE n2onoy_read
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_read_simnoy
