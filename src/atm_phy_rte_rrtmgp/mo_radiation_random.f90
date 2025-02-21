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

MODULE mo_radiation_random
  USE mo_kind, ONLY: dp, i8
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PARAMETER:: big_endian = (transfer(1_i8, 1) == 0)

  INTEGER, PARAMETER :: seed_size = 4
  
  PUBLIC :: get_random, seed_size, get_random_rank3
  INTERFACE get_random_rank3
    MODULE PROCEDURE kissvec_mask_rank3, kissvec_all_rank3
  END INTERFACE get_random_rank3
  INTERFACE get_random
    MODULE PROCEDURE kissvec_masked, kissvec_all 
  END INTERFACE get_random
  
CONTAINS








! George Masaglia's KISS random number generator
! adapted from public domain code available at http://www.fortran.com/kiss.f90

  SUBROUTINE kissvec_all(kproma, kbdim, seed, harvest)
    INTEGER, INTENT(IN) :: kproma, kbdim
    INTEGER, INTENT(INOUT) :: seed(:,:) ! Dimension kbdim, seed_size or bigger
    REAL(DP), INTENT(INOUT) :: harvest(KBDIM) 
    
    INTEGER(i8) :: kiss
    INTEGER :: jk 
    
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO jk = 1, kproma
      seed(jk,1) = transfer(69069_i8 * seed(jk,1) + 1327217885,1)
      seed(jk,2) = (ieor (seed(jk,2), ishft (seed(jk,2),  13)))
      seed(jk,2) = (ieor (seed(jk,2), ishft (seed(jk,2),  - 17)))
      seed(jk,2) = (ieor (seed(jk,2), ishft (seed(jk,2),  5)))
      seed(jk,3) = 18000 * iand (seed(jk,3), 65535) + ishft (seed(jk,3), - 16)
      seed(jk,4) = 30903 * iand (seed(jk,4), 65535) + ishft (seed(jk,4), - 16)
      kiss       = int(seed(jk,1), i8) + seed(jk,2) + ishft (seed(jk,3), 16) + seed(jk,4)
      harvest(jk) = transfer(kiss,1)*2.328306e-10_dp + 0.5_dp
    END DO
  END SUBROUTINE kissvec_all

  SUBROUTINE kissvec_all_rank3(kproma, kbdim, klev, ksamps, seed, harvest)
    INTEGER, INTENT(IN) :: kproma, kbdim, klev, ksamps
    INTEGER, INTENT(INOUT) :: seed(:,:) ! Dimension kbdim, seed_size or bigger
    REAL(DP), INTENT(INOUT) :: harvest(KBDIM,klev,ksamps) 
    
    INTEGER(i8) :: kiss
    INTEGER :: i,j,k
    
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO k = 1, ksamps
    !$ACC LOOP SEQ
    DO j = klev, 1, -1
    !$ACC LOOP GANG VECTOR
    DO i = 1, kproma
      seed(i,1) = transfer(69069_i8 * seed(i,1) + 1327217885,1)
      seed(i,2) = (ieor (seed(i,2), ishft (seed(i,2),  13)))
      seed(i,2) = (ieor (seed(i,2), ishft (seed(i,2),  - 17)))
      seed(i,2) = (ieor (seed(i,2), ishft (seed(i,2),  5)))
      seed(i,3) = 18000 * iand (seed(i,3), 65535) + ishft (seed(i,3), - 16)
      seed(i,4) = 30903 * iand (seed(i,4), 65535) + ishft (seed(i,4), - 16)
      kiss      = int(seed(i,1), i8) + seed(i,2) + ishft (seed(i,3), 16) + seed(i,4)
      harvest(i,j,k) = transfer(kiss,1)*2.328306e-10_dp + 0.5_dp
    END DO
    END DO
    END DO
    !$ACC END PARALLEL
  END SUBROUTINE kissvec_all_rank3

  SUBROUTINE kissvec_mask_rank3(kproma, kbdim, klev, ksamps, seed, &
    mask, harvest)
    INTEGER, INTENT(IN) :: kproma, kbdim, klev, ksamps
    INTEGER, INTENT(INOUT) :: seed(:,:) ! Dimension kbdim, seed_size or bigger
    LOGICAL,  INTENT(IN) :: mask(KBDIM,klev)    
    REAL(DP), INTENT(INOUT) :: harvest(KBDIM,klev,ksamps) 
    
    INTEGER(i8) :: kiss
    INTEGER :: i,j,k
    
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO k = 1, ksamps
    !$ACC LOOP SEQ
    DO j = klev, 1, -1
    !$ACC LOOP GANG VECTOR
    DO i = 1, kproma
      IF (mask(i,j)) THEN
        seed(i,1) = transfer(69069_i8 * seed(i,1) + 1327217885,1)
        seed(i,2) = (ieor (seed(i,2), ishft (seed(i,2),  13)))
        seed(i,2) = (ieor (seed(i,2), ishft (seed(i,2),  - 17)))
        seed(i,2) = (ieor (seed(i,2), ishft (seed(i,2),  5)))
        seed(i,3) = 18000 * iand (seed(i,3), 65535) + ishft (seed(i,3), - 16)
        seed(i,4) = 30903 * iand (seed(i,4), 65535) + ishft (seed(i,4), - 16)
        kiss      = int(seed(i,1), i8) + seed(i,2) + ishft (seed(i,3), 16) + seed(i,4)
        harvest(i,j,k) = transfer(kiss,1)*2.328306e-10_dp + 0.5_dp
      ELSE
        harvest(i,j,k) = 0._dp
      ENDIF
    END DO
    END DO
    END DO
    !$ACC END PARALLEL
  END SUBROUTINE kissvec_mask_rank3


!NOTE: not invoked in runtime, see mo_radiation_cld_sampling.f90
  SUBROUTINE kissvec_masked(kproma, kbdim, seed, mask, harvest)
    INTEGER, INTENT(IN) :: kproma, kbdim
    INTEGER, INTENT(INOUT) :: seed(:,:) ! Dimension kbdim, seed_size or bigger
    LOGICAL,  INTENT(IN) :: mask(KBDIM)    
    REAL(DP), INTENT(INOUT) :: harvest(KBDIM) 
    
    INTEGER(i8) :: kiss 
    INTEGER     :: jk 
    
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO jk = 1, kproma
      IF(mask(jk)) THEN  
        seed(jk,1) = transfer(69069_i8 * seed(jk,1) + 1327217885,1)
        seed(jk,2) = (ieor (seed(jk,2), ishft (seed(jk,2),  13)))
        seed(jk,2) = (ieor (seed(jk,2), ishft (seed(jk,2),  - 17)))
        seed(jk,2) = (ieor (seed(jk,2), ishft (seed(jk,2),  5)))
        seed(jk,3) = 18000 * iand (seed(jk,3), 65535) + ishft (seed(jk,3), - 16)
        seed(jk,4) = 30903 * iand (seed(jk,4), 65535) + ishft (seed(jk,4), - 16)
        kiss       = int(seed(jk,1), i8) + seed(jk,2) + ishft (seed(jk,3), 16) + seed(jk,4)
        harvest(jk) = transfer(kiss,1)*2.328306e-10_dp + 0.5_dp
      ELSE  
        harvest(jk) = 0._dp
      END IF 
    END DO
  END SUBROUTINE kissvec_masked
  
END MODULE mo_radiation_random
