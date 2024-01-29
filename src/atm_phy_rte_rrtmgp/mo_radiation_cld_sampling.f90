! Module to provide Monte Carlo samples based on cloud fraction
!
! Remarks
!   This module draws Monte Carlo samples from profiles of the atmospheric
!   state and on overlap assumption ( 1=max-ran, 2=maximum, 3=random).
!   Users provide profiles of cloud fraction and in-cloud optical thickness.
!   A single sample is drawn; cloud fraction is either 1 or 0 and optical
!   thickness is set to 0 where cloud fraction is 0.
!
! Origin
!   Written by Robert Pincus; simplifed from code originally written for the
!   GFDL atmospheric model AM2.
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

MODULE mo_radiation_cld_sampling 

  USE mo_kind, ONLY: wp
  USE mo_exception, ONLY: finish
  USE mo_radiation_random, ONLY: get_random, get_random_rank3
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sample_cld_state

  INTEGER, PARAMETER :: mode = 1 !< 1=max-ran, 2=maximum, 3=random 
CONTAINS

!c.f. rrtmg_lw's mcica_subcol_gen_lw.f90:generate_stochastic_clouds

  SUBROUTINE sample_cld_state(kproma, kbdim, klev, ksamps, &
    rnseeds, cld_frc, is_cloudy)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: &
      kproma, kbdim, klev, ksamps
    INTEGER,  INTENT(INOUT) :: rnseeds(:,:)
    REAL(wp), INTENT(IN) :: cld_frc(kbdim,klev)
    LOGICAL,  INTENT(INOUT) :: is_cloudy(kbdim,klev,ksamps)

    REAL(wp) :: rank(kbdim,klev,ksamps), one_minus(KBDIM,klev)
    INTEGER  :: jk, js, jl

    !$ACC DATA PRESENT(rnseeds, cld_frc, is_cloudy) &
    !$ACC   CREATE(rank, one_minus)

    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
    one_minus(1:kproma,:) = 1.0_wp - cld_frc(1:kproma,:)
    !$ACC END KERNELS

    ! Here is_cloudy(:,:,1) indicates whether any cloud is present 
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
    is_cloudy(1:kproma,1:klev,1) = cld_frc(1:kproma,1:klev) > 0._wp
    !$ACC END KERNELS

    SELECT CASE(mode) 
      ! Maximum-random overlap
      CASE(1) 
            ! mask means we compute random numbers only when cloud is present 
        CALL get_random_rank3(kproma, kbdim, klev, ksamps, rnseeds, &
          is_cloudy(:,:,1), rank)
        ! There may be a better way to structure this calculation...
        !$ACC PARALLEL DEFAULT(PRESENT) FIRSTPRIVATE(ksamps, kproma, klev) ASYNC(1)
        !$ACC LOOP SEQ
        DO jk = 2, klev
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO js = 1, ksamps
            DO jl = 1, kproma
              rank(jl,jk,js) = MERGE( &
                ! Max overlap:
                rank(jl,jk-1,js), & 
                ! ... or random overlap in the clear sky portion:  
                rank(jl,jk,js) * one_minus(jl,jk-1), & 
                ! depending on whether or not you have cloud in the layer above 
                rank(jl,jk-1,js) > one_minus(jl,jk-1) )
            END DO
          END DO
        END DO
        !$ACC END PARALLEL

      ! Max overlap means every cell in a column is identical 
      CASE(2) 
        DO js = 1, ksamps
          CALL get_random(kproma, kbdim, rnseeds, rank(:, 1, js))

          ! rank(1:kproma,2:klev,js) = SPREAD(rank(1:kproma,1,js), &
          !   DIM=2, NCOPIES=(klev-1))
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) FIRSTPRIVATE(klev, kproma, js) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO jk = 2, klev
            DO jl = 1, kproma
              rank(jl,jk,js) = rank(jl,1,js)
            END DO
          END DO
          !$ACC END PARALLEL LOOP
        END DO 

      ! Random overlap means every cell is independent
      CASE(3) 
        ! DO js = 1, ksamps
        !   DO jk = klev, 1, -1
        !     ! mask means we compute random numbers only when cloud is present 
        !     CALL get_random(kproma, kbdim, rnseeds, is_cloudy(:,jk,1), &
        !       rank(:,jk,js))
        !   END DO 
        ! END DO
        !
        ! Seems this is essentially the same as the previous commented out part
        CALL get_random_rank3(kproma, kbdim, klev, ksamps, rnseeds, &
          is_cloudy(:,:,1), rank)
      CASE DEFAULT
        CALL finish('In sample_cld_state: unknown overlap assumption') 
    END SELECT
    ! Now is_cloudy indicates whether the sample (ks) is cloudy or not. 
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) FIRSTPRIVATE(ksamps, klev, kproma) GANG VECTOR COLLAPSE(3) ASYNC(1)
    DO js = 1, ksamps
      DO jk = 1, klev
        DO jl = 1, kproma
          is_cloudy(jl,jk,js) = &
            rank(jl,jk,js) > one_minus(jl,jk)
        END DO
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC WAIT
    !$ACC END DATA
  END SUBROUTINE sample_cld_state

END MODULE mo_radiation_cld_sampling
