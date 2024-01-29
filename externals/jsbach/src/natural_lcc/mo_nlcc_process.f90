!> Contains the routines for the nlcc processes
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_nlcc_process
#ifndef __NO_JSBACH__

  USE mo_jsb_control, ONLY: debug_on
  USE mo_kind,        ONLY: wp
  USE mo_exception,   ONLY: message, finish, message_text
  USE mo_jsb_physical_constants, ONLY: tmelt

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_climbuf, bioclim_limits, scale_fpc, potential_tree_fpc
  PUBLIC :: fpc_to_cover_fract_pot, cover_fract_pot_to_cover_fract, desert_fraction, fpc_daily

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_process'
  REAL(wp), PARAMETER :: fract_small = 1.e-10_wp   !! smallest fraction of PFT tiles
  REAL(wp), PARAMETER :: act_fpc_min = 0.005_wp    !! minimum actual FPC to be considered to calculate potential FPC
  REAL(wp), PARAMETER :: sum_npp_min = 1.e-12_wp   !! minimum total npp for tree cover
  REAL(wp), PARAMETER :: tree_fpc_max = 1.0_wp     !! maximum FPC of woody PFTs
  REAL(wp), PARAMETER :: npp_nonlinearity = 1.5_wp !! parameter controlling the non-linearity in the dynamic equation with
                                                   !! respect to NPP
  REAL(wp), PARAMETER :: desert_extend = 0.65_wp   !! parameter controlling the extend of desert
                                                   !! (the lower the value the more desert)
  REAL(wp), PARAMETER :: desert_margin = 2._wp     !! parameter controlling the transition from vegetated land to desert
                                                   !! (the higher the value the sharper the transition)
  REAL(wp), PARAMETER :: tau_desert = 50._wp       !! time constant [year] by which the extension of deserts is adapted

CONTAINS
  !
  ! ===============================================================================================================================
  !>
  !! Compute running averages of climate variables for nlcc
  !!

  ELEMENTAL PURE SUBROUTINE calc_climbuf (lstart, init_running_means, new_day, new_month, new_year, delta_time, t_air, &
                                          seconds_day, seconds_month, temp_sum_day, temp_sum_month, &
                                          min_mmtemp_of_yr, max_mmtemp_of_yr, min_mmtemp20, max_mmtemp20, &
                                          gdd_sum_year, gdd_prev_year)

    LOGICAL, INTENT(in) :: lstart ! first time step in a simulation
    LOGICAL, INTENT(in) :: init_running_means ! initialize the min/max temperature climatology
    LOGICAL, INTENT(in) :: new_day ! first time step in a day
    LOGICAL, INTENT(in) :: new_month ! firest time step in a month
    LOGICAL, INTENT(in) :: new_year ! first time step in a year
    REAL(wp), INTENT(in) :: delta_time ! time step length [s]
    REAL(wp), INTENT(in) :: t_air ! air temperature of lowest atmosphere level [K]
    REAL(wp), INTENT(inout) :: seconds_day ! sum of seconds till new_day 
    REAL(wp), INTENT(inout) :: seconds_month ! sum of seconds till new_month
    REAL(wp), INTENT(inout) :: temp_sum_day
    REAL(wp), INTENT(inout) :: temp_sum_month
    REAL(wp), INTENT(inout) :: min_mmtemp_of_yr ! temperature of coldest month of previous year [K]
    REAL(wp), INTENT(inout) :: max_mmtemp_of_yr ! temperature of warmest month of previous year [K]
    REAL(wp), INTENT(inout) :: min_mmtemp20 ! temperature of coldest month climatology [K]
    REAL(wp), INTENT(inout) :: max_mmtemp20 ! temperature of warmest month climatology [K]
    REAL(wp), INTENT(inout) :: gdd_sum_year ! growing degree days of current year [C]
    REAL(wp), INTENT(inout) :: gdd_prev_year ! growing degree days of previous year [C]

!   local variables
    REAL(wp) :: prev_day_mean_temp   ! air temperature of previous day [K]
    REAL(wp) :: prev_month_mean_temp ! air temperature of previous month [K]


    IF (new_year) THEN
       gdd_prev_year = gdd_sum_year
       gdd_sum_year = 0._wp
    ENDIF

    IF (.NOT. new_day .OR. lstart) THEN
       temp_sum_day = temp_sum_day + t_air * delta_time
    ELSE
       prev_day_mean_temp = temp_sum_day / seconds_day
       temp_sum_day = t_air * delta_time
       gdd_sum_year = gdd_sum_year + MAX(0._wp,prev_day_mean_temp - tmelt - 5._wp)
    ENDIF

    IF (.NOT. new_month .OR. lstart) THEN
       temp_sum_month = temp_sum_month + t_air * delta_time
    ELSE
       prev_month_mean_temp = temp_sum_month / seconds_month
       temp_sum_month = t_air * delta_time
       ! --- minimum/maximum monthly mean temperature (20 year climatology)
       ! --- calculate coldest/warmest month of the year so far
       min_mmtemp_of_yr = min(min_mmtemp_of_yr, prev_month_mean_temp)
       max_mmtemp_of_yr = max(max_mmtemp_of_yr, prev_month_mean_temp)
       IF (new_year) THEN
          ! --- build 20yr running mean of coldest and warmest monthly temperature
          IF (init_running_means) THEN
             min_mmtemp20 = min_mmtemp_of_yr
             max_mmtemp20 = max_mmtemp_of_yr
          ELSE
             min_mmtemp20 = (min_mmtemp20 * 19._wp + min_mmtemp_of_yr) / 20._wp
             max_mmtemp20 = (max_mmtemp20 * 19._wp + max_mmtemp_of_yr) / 20._wp
          ENDIF
          ! --- reset for current year
          min_mmtemp_of_yr =  1000._wp
          max_mmtemp_of_yr = -1000._wp
       END IF ! new_year
    END IF ! new_month

    ! count time in seconds that has passed till summing of variables started
    seconds_day = seconds_day + delta_time
    IF (new_day) seconds_day = delta_time
    seconds_month = seconds_month + delta_time
    IF (new_month) seconds_month = delta_time

  END SUBROUTINE calc_climbuf

! --- bioclim_limits ---------------------------------------------------------------------------------------

  ELEMENTAL PURE SUBROUTINE bioclim_limits (is_dynamic, lct_tcmin, lct_tcmax, lct_twmax, lct_min_temprange, lct_gddmin, &
                                            min_mmtemp20, max_mmtemp20, gdd_prev_year, bio_exist)

!
! calculation of bioclimatic limits for each PFT based on modified
! LPJ lookup-table
!

! input
    LOGICAL,  INTENT(in)  :: is_dynamic        ! mask indicating whether a PFT is subject to vegetation dynamics
    REAL(WP), INTENT(in)  :: lct_tcmin         ! PFT-specific minimum coldest-month temperature limit [C]
    REAL(WP), INTENT(in)  :: lct_tcmax         ! PFT-specific maximum coldest-month temperature limit [C]
    REAL(WP), INTENT(in)  :: lct_twmax         ! PFT-specific maximum warmest-month temperature limit [C]
    REAL(WP), INTENT(in)  :: lct_min_temprange ! PFT-specific minimum difference of 20-year average warme st
                                               ! minus coldest month temperature [C]
    REAL(WP), INTENT(in)  :: lct_gddmin        ! PFT-specific minimum GDD limit
    REAL(wp), INTENT(in)  :: min_mmtemp20      ! Minimum monthly mean temp. (20yr climatology) [K]
    REAL(wp), INTENT(in)  :: max_mmtemp20      ! Maximum monthly mean temp. (20yr climatology) [K]
    REAL(wp), INTENT(in)  :: gdd_prev_year     ! GDD of previous year

! output
    REAL(wp), INTENT(inout) :: bio_exist

!------------------------------------------------------------------------------

    bio_exist = 1._wp

    IF (is_dynamic .AND. (min_mmtemp20 - tmelt < lct_tcmin .OR. min_mmtemp20 - tmelt > lct_tcmax  &
       .OR. gdd_prev_year < lct_gddmin .OR. max_mmtemp20 - tmelt > lct_twmax                      &
       .OR. (max_mmtemp20 - min_mmtemp20) < lct_min_temprange)) bio_exist = 0._wp
      
  END SUBROUTINE bioclim_limits


! --- potential_tree_fpc ---------------------------------------------------------------------------

  SUBROUTINE potential_tree_fpc(nidx, ntiles, is_dynamic, is_woody, bio_exist, &
                                npp_ave, act_fpc, pot_fpc)

!
! subroutine calculates potential FPC (in absence of disturbances) based on NPP for each PFT
! called once a year!
!

! input
    INTEGER,  INTENT(in)  :: nidx            ! vector length
    INTEGER,  INTENT(in)  :: ntiles          ! number of tiles (different land cover type mosaics per gridcell)
    LOGICAL,  INTENT(in)  :: is_dynamic(:,:) ! flag indicating dynamic PFTs
    LOGICAL,  INTENT(in)  :: is_woody(:,:)   ! flag indicating woody PFTs
    REAL(wp), INTENT(in)  :: bio_exist(:,:)  ! bio_exist=1. indicates that all bio-climatic limits are fulfilled
    REAL(wp), INTENT(in)  :: npp_ave(:,:)    ! NPP averaged over 5 years
    REAL(wp), INTENT(in)  :: act_fpc(:,:)    ! actual foliage projective cover

! output
    REAL(wp), INTENT(out) :: pot_fpc(:,:)

! local variables
    INTEGER           :: i, itile
    REAL(wp)          :: sum_npp(nidx)
      
! summing up NPP
    sum_npp(:) = 0._wp

    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (is_woody(i,itile) .AND. is_dynamic(i,itile) .AND. &
           npp_ave(i,itile) > sum_npp_min .AND. bio_exist(i,itile) > 0.5_wp) &
           sum_npp(i) = sum_npp(i) + npp_ave(i,itile) ** npp_nonlinearity * MAX(act_fpc_min,act_fpc(i,itile))
       END DO
    END DO

!-- determine potential FPC
    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (is_woody(i,itile) .AND. is_dynamic(i,itile) .AND. npp_ave(i,itile) > sum_npp_min .AND. &
              bio_exist(i,itile) > 0.5_wp) THEN
             pot_fpc(i,itile) = npp_ave(i,itile) ** npp_nonlinearity / sum_npp(i) * &
                                tree_fpc_max * MAX(act_fpc_min,act_fpc(i,itile))
          ELSE
             pot_fpc(i,itile) = 0._wp
          END IF
       END DO
    END DO

  END SUBROUTINE potential_tree_fpc


! --- fpc_to_cover_fract_pot -------------------------------------------------------------------------

  SUBROUTINE fpc_to_cover_fract_pot (nidx, ntiles, is_dynamic, is_woody, act_fpc, bare_fpc, cover_fract_pot)
!
! Convert fractional plant cover (act_fpc) calculated within nlcc to jsbach
! cover fractions reflecting natural vegetation, only (cover_fract_pot).
! The only difference between act_fpc and cover_fract_pot is caused by bare_fpc,
! which is added to the grass cover fractions.
!

! input
    INTEGER,  INTENT(in)  :: nidx            ! vector length
    INTEGER,  INTENT(in)  :: ntiles          ! number of tiles (different land cover type mosaics per gridcell)
    LOGICAL,  INTENT(in)  :: is_dynamic(:,:) !
    LOGICAL,  INTENT(in)  :: is_woody(:,:)   !
    REAL(wp), INTENT(in)  :: act_fpc(:,:)    ! actual fpc in nlcc
    REAL(wp), INTENT(in)  :: bare_fpc(:)     ! nlcc bare soil fraction
! output
    REAL(wp), INTENT(out) :: cover_fract_pot(:,:) ! cover fractions if there was only natural vegetation

! local variables
    INTEGER   :: i, itile
    REAL(wp)  :: sum_grass_fpc(nidx)     ! part of the vegetated area covered by grass

! -- initialization
    cover_fract_pot(:,:) = 0._wp

! -- find out fpc of grass
    sum_grass_fpc(:) = 0._wp

    DO i = 1,nidx
      DO itile = 1,ntiles
         IF (.NOT. is_woody(i,itile) .AND. is_dynamic(i,itile)) &
            sum_grass_fpc(i) = sum_grass_fpc(i) + act_fpc(i,itile)
      END DO
    END DO

!-- calculate new cover_fractions
    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (is_dynamic(i,itile)) THEN
             IF (is_woody(i,itile)) THEN
                cover_fract_pot(i,itile) = act_fpc(i,itile)
             ELSE
                cover_fract_pot(i,itile) = act_fpc(i,itile) * (1._wp + bare_fpc(i) / MAX(sum_grass_fpc(i),EPSILON(1._wp)))
             END IF
          END IF
       END DO
    END DO

    CALL scale_cover_fract (nidx, ntiles, is_dynamic(:,:), cover_fract_pot(:,:))

    IF (ANY(SUM(cover_fract_pot(:,:),DIM=2) > 1._wp + ntiles*EPSILON(1._wp)) .OR. &
        ANY(SUM(cover_fract_pot(:,:),DIM=2) < 1._wp - ntiles*EPSILON(1._wp))) THEN
       WRITE (message_text,*) 'sum of cover_fract_pot /= 1: ', &
            MINVAL(SUM(cover_fract_pot(:,:),DIM=2)), MAXVAL(SUM(cover_fract_pot(:,:),DIM=2)), &
            MAXLOC(SUM(cover_fract_pot(:,:),DIM=2))
       CALL finish ('fpc_to_cover_fract_pot', message_text)
    END IF
    DO itile = 1,ntiles
       IF (ANY(is_dynamic(:,itile) .AND. cover_fract_pot(:,itile) < fract_small)) THEN
          WRITE (message_text,*) 'cover_fract_pot too small: ', MINVAL(cover_fract_pot(:,itile)), itile
          CALL finish ('fpc_to_cover_fract_pot', message_text)
       END IF
    END DO

  END SUBROUTINE fpc_to_cover_fract_pot


! --- cover_fract_pot_to_cover_fract ---------------------------------------------------------------------------

  SUBROUTINE cover_fract_pot_to_cover_fract (nidx, ntiles, is_dynamic, is_woody, is_pasture, is_crop, &
                                             cover_fract_pot, cover_fract_pot_previous_year, cover_fract)

!
! convert the fractional plant cover of potential natural vegetation to plant cover fractions including agricultural areas.
! Pasture is first established on grasslands, whereas crops are established at the expense of all natural cover types.
!
! In runs with land use transitions, the transitions from natural forests or grass lands to crops and pastures are predefined.
! It is not wanted, that the dynamic vegetation interferes with these transitions by e.g. re-establishing forests where they
! had just been removed. Thus the changes in forest fraction calculated by the dynamic vegetation are not instantanously 
! applied. This results in an imbalance between actual cover_fractions and the ones the dynamic vegetation calculates for
! instantanous establishement of crops and pasture (cover_fract_inst).
! The final cover_fractions are based on cover_fract_inst with the exception, that the forest fraction is based on the forest
! fraction after land use transitions (cover_fract entering the routine), and only the natural change in cover_fractions in
! the direction to reduce the inbalance is taken into account (see below).
!
! input
    INTEGER,  INTENT(in)  :: nidx                 ! vector length
    INTEGER,  INTENT(in)  :: ntiles               ! number of tiles (different land cover type mosaics per gridcell)
    LOGICAL,  INTENT(in)  :: is_dynamic(:,:)      ! PFTs taking part in competition
    LOGICAL,  INTENT(in)  :: is_woody(:,:)        ! PFTs that build wood
    LOGICAL,  INTENT(in)  :: is_pasture(:,:)      ! logical mask for pasture
    LOGICAL,  INTENT(in)  :: is_crop(:,:)         ! logical mask for crops
    REAL(wp), INTENT(in)  :: cover_fract_pot(:,:) ! cover fractions if there was only natural vegetation
    REAL(wp), INTENT(in)  :: cover_fract_pot_previous_year(:,:) ! last years cover_fract_pot

! inout
    REAL(wp), INTENT(inout) :: cover_fract(:,:)  ! cover fraction within jsbach

! local variables

    INTEGER   :: i, itile
    REAL(wp)  :: nwoody(nidx)                   ! number of woody PFTs (as REAL)
    REAL(wp)  :: nexcluded(nidx)                ! number of non-dynamic PFTs (crops, pastures)
    REAL(wp)  :: cover_fract_inst(nidx,ntiles)  ! cover fractions that would result from an instantanous implementation of
                                                ! land use change and vegetation dynamics
    REAL(wp)  :: excluded_fract(nidx)           ! fraction of the vegetated part of the grid box not considered for vegetation
                                                ! dynamics (e.g. agricultural areas)
    REAL(wp)  :: sum_woody_fract(nidx)          ! fraction of the vegetated part of the grid box covered by woody plants
    REAL(wp)  :: sum_grass_fract(nidx)          ! fraction of the vegetated part of the grid box covered by grass
    REAL(wp)  :: sum_woody_fract_inst(nidx)     ! fraction of the vegetated part of the grid box covered by woody plants
                                                ! based on cover_fract_inst
    REAL(wp)  :: sum_woody_fract_pot(nidx)      ! woody fraction of the vegetated part of the grid box considering natural
                                                ! vegetation, only
    REAL(wp)  :: sum_grass_fract_pot(nidx)      ! grass fraction of the vegetated part of the grid box considering natural
                                                ! vegetation, only
    REAL(wp)  :: sum_woody_fract_old(nidx)      ! fraction of the vegetated part of the grid box covered by woody plants
                                                ! (based on act_fpc_previous_year) that would have resulted from an
                                                ! instantaneous establishment of crop land and pasture
    REAL(wp)  :: sum_woody_fract_act(nidx)      ! fraction of the vegetated part of the grid box covered by woody plants
    REAL(wp)  :: sum_grass_fract_act(nidx)      ! fraction of the vegetated part of the grid box covered by grasses
    REAL(wp)  :: delta_woody(nidx)              ! deviation of the fraction of the vegetated part of the grid box covered 
                                                ! by woody plants from the fraction that would result from an instantaneous
                                                !  establishment of crop land and pasture

    !-- find out number of woody types and number of exculded types
    nwoody(:) = 0._wp
    nexcluded(:) = 0._wp
    DO i=1,nidx
       DO itile = 1,ntiles
          IF (is_dynamic(i,itile) .AND. is_woody(i,itile)) nwoody(i) = nwoody(i) + 1._wp
          IF (.NOT. is_dynamic(i,itile)) nexcluded(i) = nexcluded(i) + 1._wp
       ENDDO
    ENDDO

    !-- Sum up cover fractions and of woody types, grasses, pasture etc. The fractions (cover_fract) were last updated by
    !   last years land use change, and do not yet include current vegetation dynamics.
    sum_woody_fract(:) = 0._wp
    sum_grass_fract(:) = 0._wp
    excluded_fract(:) = 0._wp

    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (is_woody(i,itile) .AND. is_dynamic(i,itile)) THEN
             sum_woody_fract(i) = sum_woody_fract(i) + cover_fract(i,itile)
          ELSE IF (.NOT. is_woody(i,itile) .AND. is_dynamic(i,itile)) THEN
             sum_grass_fract(i) = sum_grass_fract(i) + cover_fract(i,itile)
          ELSE IF (.NOT. is_dynamic(i,itile)) THEN
             excluded_fract(i) = excluded_fract(i) + cover_fract(i,itile)
          END IF
       END DO
    END DO

    !-- Sum up woody and grass fractions of the potential natural vegetation based on current vegetation dynamics
    sum_woody_fract_pot(:) = 0._wp
    sum_grass_fract_pot(:) = 0._wp

    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (is_woody(i,itile) .AND. is_dynamic(i,itile)) THEN
             sum_woody_fract_pot(i) = sum_woody_fract_pot(i) + cover_fract_pot(i,itile)
          ELSE IF (.NOT. is_woody(i,itile) .AND. is_dynamic(i,itile)) THEN
             sum_grass_fract_pot(i) = sum_grass_fract_pot(i) + cover_fract_pot(i,itile)
          END IF
       END DO
    END DO

    !-- Calculate the woody fraction resulting from last years potential vegetation assuming an instantanous
    !   establishment of crops and pastures
    CALL calc_cover_fract_inst(nidx, ntiles, is_woody(:,:), is_dynamic(:,:), is_pasture(:,:), is_crop(:,:), &
                               cover_fract_pot_previous_year(:,:), cover_fract(:,:), cover_fract_inst(:,:))

    sum_woody_fract_old(:) = 0._wp
    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (is_woody(i,itile) .AND. is_dynamic(i,itile)) &
               sum_woody_fract_old(i) = sum_woody_fract_old(i) + cover_fract_inst(i,itile)
       END DO
    END DO

    !-- Calculate the woody fraction resulting from current potential vegetation assuming an instantanous
    !   establishment of crops and pastures
    CALL calc_cover_fract_inst(nidx, ntiles, is_woody(:,:), is_dynamic(:,:), is_pasture(:,:), is_crop(:,:), &
                               cover_fract_pot(:,:), cover_fract(:,:), cover_fract_inst(:,:))

    sum_woody_fract_inst(:) = 0._wp
    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (is_woody(i,itile) .AND. is_dynamic(i,itile)) THEN
             sum_woody_fract_inst(i) = sum_woody_fract_inst(i) + cover_fract_inst(i,itile)
          END IF
       END DO
    END DO

    !-- Calculation of the final woody fraction. A difference between sum_woody_fract and sum_woody_fract_old indicates
    !   an imbalance of the state of the dynamic vegetation and the cover fractions calculated from landuse transitions.
    !   On the other hand, the difference between sum_woody_fract_inst and sum_woody_fract_old show the trend of the woody
    !   fraction caused by vegetation dynamics.
    !   The idea of this routine is, to allow shifts in forest fractions only if they reduce the imbalance between nlcc
    !   and landuse change.
    DO i = 1,nidx
       IF (sum_woody_fract(i) - sum_woody_fract_old(i) > 0._wp) THEN
          ! -- woody fraction of cover_fract entering the routine (last updated by last years land use change) is greater
          !    than the the woody fraction resulting from last years potential vegetation and instantanous establishment of
          !    crops and pastures. Thus there is an inbalance, the woody fraction of cover_fract is too high.
          IF (sum_woody_fract_inst(i) - sum_woody_fract_old(i) > 0._wp) THEN
             ! -- vegetation dynamics even enlarge woody fraction
             !    ==> no change in actual woody fraction, unless the enlargement by vegetation dynamics exceeds the initial
             !        woody surplus
             delta_woody(i) = MAX(0._wp, sum_woody_fract_inst(i) - sum_woody_fract(i))
          ELSE
             ! -- vegetation dynamics shrink woody fraction
             delta_woody(i) = sum_woody_fract_inst(i) - sum_woody_fract_old(i)
          END IF
       ELSE
          ! -- The woody fraction of cover_fract is too small
          IF (sum_woody_fract_inst(i) - sum_woody_fract_old(i) > 0._wp) THEN
             ! -- vegetation dynamics enlarge woody fraction
             delta_woody(i) = sum_woody_fract_inst(i) - sum_woody_fract_old(i)
          ELSE
             ! -- vegetation dynamics shrink woody fraction
             !    ==> no change in actual woody fraction, unless the reduction by vegetation dynamics exceeds the initial
             !        woody deficit
             delta_woody(i) = MIN(0._wp, sum_woody_fract_inst(i) - sum_woody_fract(i))
          END IF
       END IF

       ! -- in some cases delta_woody can not be used to calculate the actual woody fraction. If there is no grass left
       !    after instantanous establishement of agricultural areas, the woody fractions of cover_fract_inst and
       !    cover_fract_inst_old are identical and delta_woody is zero, even if vegetation dynamics change considerably.
       !    The following formulation is needed to assure that sum_woody_fract_act is smaller than sum_woody_fract_pot and
       !    sum_grass_fract_act is smaller than sum_grass_fract_pot in all cases.
       sum_woody_fract_act(i) = MIN(sum_woody_fract_pot(i),&
                                MAX(nwoody(i)*fract_small, sum_woody_fract_pot(i)-excluded_fract(i)+nexcluded(i)*fract_small, &
                                    sum_woody_fract(i) + delta_woody(i)))
    END DO

    sum_grass_fract_act(:) = 1._wp - sum_woody_fract_act(:) - excluded_fract(:)

    !-- Use actual woody fraction to determine all actual fractions
    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (is_dynamic(i,itile)) THEN
             IF (is_woody(i,itile)) THEN
                cover_fract(i,itile) = cover_fract_pot(i,itile) * sum_woody_fract_act(i) / sum_woody_fract_pot(i)
             ELSE
                cover_fract(i,itile) = cover_fract_pot(i,itile) * sum_grass_fract_act(i) / sum_grass_fract_pot(i)
             END IF
          END IF
       END DO
    END DO

    CALL scale_cover_fract(nidx, ntiles, is_dynamic(:,:), cover_fract(:,:))

    !-- Stop simulation, if sum of cover fractions is not equal to 1 or a cover fraction is smaller than fract_small
    IF (ANY(SUM(cover_fract(:,:),DIM=2) > 1._wp + REAL(ntiles) * EPSILON(1._wp)) .OR. &
        ANY(SUM(cover_fract(:,:),DIM=2) < 1._wp - REAL(ntiles) * EPSILON(1._wp))) THEN
       WRITE (message_text,*) 'sum of cover_fract /= 1: ', &
            MINVAL(SUM(cover_fract(:,:),DIM=2)), MAXVAL(SUM(cover_fract(:,:),DIM=2)), &
            MAXLOC(SUM(cover_fract(:,:),DIM=2))
       CALL finish ('cover_fract_pot_to_cover_fract', message_text)
    END IF
    DO itile = 1,ntiles
       IF (ANY(cover_fract(:,itile) < fract_small)) THEN
          WRITE (message_text,*) 'cover_fract too small: ', MINVAL(cover_fract(:,itile))
          CALL finish ('cover_fract_pot_to_cover_fract', message_text)
       END IF
    END DO

  END SUBROUTINE cover_fract_pot_to_cover_fract


! --- desert_fraction --------------------------------------------------------

  SUBROUTINE desert_fraction(nidx, ntiles, init_running_means, accelerate_nlcc,   &
                             is_woody, is_dynamic,                                &
                             act_fpc, bare_fpc,                                   &
                             sla, max_green_bio,                                  &
                             sum_green_bio_memory, desert_fpc)
! !DESCRIPTION:
!
! Calculation of desert fraction from green biomass
!
!------------------------------------------------------------------------------

! input
    INTEGER,  INTENT(in)  :: nidx              ! vector length
    INTEGER,  INTENT(in)  :: ntiles            ! number of tiles (different land cover type mosaics per gridcell)
    LOGICAL,  INTENT(in)  :: init_running_means! initialisation of sum_green_bio_memory
    REAL(wp), INTENT(in)  :: accelerate_nlcc   ! acceleration factor for vegetation and desert dynamics
    LOGICAL,  INTENT(in)  :: is_woody(:,:)     ! flag to label a woody PFT
    LOGICAL,  INTENT(in)  :: is_dynamic(:,:)   ! flag to label PFTs that take part in the natural vegetation dynamics
    REAL(wp), INTENT(in)  :: act_fpc(:,:)      ! actual FPC
    REAL(wp), INTENT(in)  :: bare_fpc(:)       ! bare FPC
    REAL(wp), INTENT(in)  :: sla(:,:)          ! specific leaf area
    REAL(wp), INTENT(in)  :: max_green_bio(:,:)! maximum value of green biomass within a year

! input/output
    REAL(wp), INTENT(inout) :: sum_green_bio_memory(:) ! vegetated fraction calculated from green biomass

! output
    REAL(wp), INTENT(out) :: desert_fpc(:) ! desert FPC

! local variables
    INTEGER           :: i, itile
    REAL(wp)          :: sum_green_bio(nidx)
    REAL(wp)          :: sum_act_fpc(nidx)
    REAL(wp)          :: sum_grass_fpc(nidx)
    INTEGER           :: n_grass_pft(nidx)

!------------------------------------------------------------------------------

    sum_green_bio(:) = 0._wp
    sum_act_fpc(:) = 0._wp
    sum_grass_fpc(:) = 0._wp
    n_grass_pft(:) = 0

    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (.NOT. is_woody(i,itile) .AND. is_dynamic(i,itile)) THEN
             n_grass_pft(i) = n_grass_pft(i) + 1
             sum_grass_fpc(i) = sum_grass_fpc(i) + act_fpc(i,itile)
          END IF
       END DO
    END DO

    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (is_woody(i,itile) .AND. is_dynamic(i,itile)) THEN
             sum_green_bio(i) = sum_green_bio(i) + MAX(0._wp, act_fpc(i,itile) * (1._wp -          &
               exp(-desert_extend * ((max_green_bio(i,itile) * sla(i,itile)) ** desert_margin) /   &
               (3._wp ** (desert_margin - 1._wp)))))
             sum_act_fpc(i) = sum_act_fpc(i) + act_fpc(i,itile)
          ELSE IF (.NOT. is_woody(i,itile) .AND. is_dynamic(i,itile)) THEN
             IF (sum_grass_fpc(i) > EPSILON(1._wp)) THEN
                sum_green_bio(i) = sum_green_bio(i) +                                              &
                 MAX(0._wp, act_fpc(i,itile) * (1._wp + bare_fpc(i) / sum_grass_fpc(i)) * (1._wp - &
                 exp(-desert_extend * ((max_green_bio(i,itile) * sla(i,itile)) ** desert_margin) / &
                (3._wp ** (desert_margin - 1._wp)))))
                sum_act_fpc(i) = sum_act_fpc(i) + act_fpc(i,itile) * (1._wp + bare_fpc(i) / sum_grass_fpc(i))
             ELSE
                sum_green_bio(i) = sum_green_bio(i) +                                              &
                 MAX(0._wp, (bare_fpc(i) / REAL(n_grass_pft(i))) * (1._wp -                        &
                 exp(-desert_extend * ((max_green_bio(i,itile) * sla(i,itile)) ** desert_margin) / &
                 (3._wp ** (desert_margin - 1._wp)))))
                sum_act_fpc(i) = sum_act_fpc(i) + (bare_fpc(i) / REAL(n_grass_pft(i)))
             END IF
          END IF
       END DO
    END DO

    WHERE (sum_act_fpc(:) > EPSILON(1._wp))
       sum_green_bio(:) = sum_green_bio(:) / sum_act_fpc(:)
    END WHERE
    IF (init_running_means) THEN
       sum_green_bio_memory(:) = sum_green_bio(:)
    ELSE
       sum_green_bio_memory(:) = (sum_green_bio_memory(:) * (tau_desert/accelerate_nlcc - 1._wp) + sum_green_bio(:)) &
                                 / (tau_desert/accelerate_nlcc)
    END IF

    desert_fpc(:) = MIN(1._wp - 2._wp * EPSILON(1._wp),MAX(0._wp,1._wp - sum_green_bio_memory(:)))

  END SUBROUTINE desert_fraction

! --- fpc_daily ---------------------------------------------------------------------------

  SUBROUTINE fpc_daily(nidx, ntiles, accelerate_nlcc, is_woody, is_dynamic, tau_Cpool_woods,    &
                       npp_ave, bio_exist, burned_fract, damaged_fract, pot_fpc, &
                       act_fpc, bare_fpc, fuel)

! Calculation of daily FPC for trees, shrubs, and herbaceous plants

! input
    INTEGER,         INTENT(in) :: nidx            ! vector length
    INTEGER,         INTENT(in) :: ntiles          ! number of tiles (different land cover type mosaics per gridcell)
    REAL(wp),        INTENT(in) :: accelerate_nlcc ! factor to accellerate vegetation and desert dynamics
    LOGICAL,         INTENT(in) :: is_woody(:,:)
    LOGICAL,         INTENT(in) :: is_dynamic(:,:)
    REAL(wp),        INTENT(in) :: tau_Cpool_woods(:,:) ! lifetime of Cpool_woods [days]
    REAL(wp),        INTENT(in) :: npp_ave(:,:)
    REAL(wp),        INTENT(in) :: bio_exist(:,:)
    REAL(wp),        INTENT(in) :: burned_fract(:,:)
    REAL(wp),        INTENT(in) :: damaged_fract(:,:)
    REAL(wp),        INTENT(in) :: pot_fpc(:,:)

! inout
    REAL(wp), INTENT(inout) :: act_fpc(:,:)
    REAL(wp), INTENT(inout) :: bare_fpc(:)
    REAL(wp), INTENT(inout), optional :: fuel(:)

! local variables
    INTEGER     :: i, j, itile
    REAL(wp)    :: act_fpc_previous_day(nidx,ntiles)
    REAL(wp)    :: total_act_fpc(nidx)
    REAL(wp)    :: grass_sum(nidx)
    REAL(wp)    :: non_woody_fpc(nidx)
    REAL(wp)    :: woody_estab_fpc(nidx)
    REAL(wp)    :: theta, veg_inc, tau_pft
    REAL(wp), PARAMETER  ::  min_dist_woody = 0.002_wp / 365._wp !! mortality of woody plants (not caused by fire/windthrow,
                                                                 !! but by other unresolved processes)
    REAL(wp), PARAMETER  ::  grass_mortality = 0.01_wp / 365._wp !! grass mortality (not caused by fire, but by other
                                                                 !! unresolved processes)
    REAL(wp), PARAMETER  ::  tau_Cpool_woods_to_tau_pft = 1._wp

    ! initialisations

    act_fpc_previous_day(:,:) = act_fpc(:,:)
    non_woody_fpc(:) = bare_fpc(:)

    DO itile = 1,ntiles
       DO i = 1,nidx
!!$ TR          IF (is_dynamic(i,itile) .AND. is_woody(i,itile)) burned_local(i,itile) = 3.e-5_wp
!!$ TR          IF (is_dynamic(i,itile) .AND. .NOT. is_woody(i,itile)) burned_local(i,itile) = 1.e-4_wp
!!$ TR          IF (is_dynamic(i,itile) .AND. is_woody(i,itile)) damaged_local(i,itile) = 1.e-5_wp
          IF (.NOT. is_dynamic(i,itile)) act_fpc(i,itile) = 0._wp
          IF (.NOT. is_woody(i,itile)) non_woody_fpc(i) = non_woody_fpc(i) + act_fpc(i,itile)
       END DO
    END DO

    !-- adjust grass fractions to changes in woody type fpc

    !-- reset temporary arrays
    grass_sum(:)     = 0._wp

    !-- temporary sums for grass types
    DO itile = 1,ntiles
       DO i = 1,nidx
          IF (.NOT. is_woody(i,itile) .AND. is_dynamic(i,itile) .AND. bio_exist(i,itile) > 0.5_wp .AND. &
              npp_ave(i,itile) > EPSILON(1._wp)) THEN
            grass_sum(i) = grass_sum(i) + npp_ave(i,itile) * MAX(act_fpc_min,act_fpc(i,itile))
          END IF
       END DO
    END DO

!!$    !-- area wind break of woody types (trees and shrubs) and burned
!!$    IF (PRESENT(burned_frac)) THEN
!!$      CALL disturbed_frac(lctlib, nidx,kidx0,kidx1,ntiles,DIST_FIRE, &
!!$                          dynamic_pft(:) .OR. (dist_opts%lburn_pasture .AND. lctlib%pasture_pft(:)), with_yasso, &
!!$                          act_fpc,veg_fract_correction,surf,cbal,climbuf,burned_frac,burned_frac_diag,fuel)
!!$      CALL disturbed_frac(lctlib, nidx,kidx0,kidx1,ntiles,DIST_WINDBREAK,woody_pft(:) .AND. dynamic_pft(:), with_yasso, &
!!$                          act_fpc,veg_fract_correction,surf,cbal,climbuf,damaged_frac) 
!!$    END IF

    !-- dynamic equation for act_fpc, woody types

    DO i = 1,nidx
      veg_inc = 0._wp
      DO itile = 1,ntiles
        tau_pft = tau_Cpool_woods(i,itile) * tau_Cpool_woods_to_tau_pft
        IF (is_dynamic(i,itile) .AND. is_woody(i,itile)) THEN
          veg_inc = veg_inc + MAX(fract_small,act_fpc(i,itile) + accelerate_nlcc * &
          ((pot_fpc(i,itile) - act_fpc(i,itile)) / tau_pft &
          - (burned_fract(i,itile) + damaged_fract(i,itile) + min_dist_woody) * act_fpc_previous_day(i,itile))) &
          - act_fpc(i,itile)
        ELSEIF(is_dynamic(i,itile) .AND. .NOT. is_woody(i,itile)) THEN
          IF (bio_exist(i,itile) > 0.5_wp .AND. npp_ave(i,itile) > EPSILON(1._wp)) THEN
            veg_inc = veg_inc + MAX(fract_small,act_fpc(i,itile) - accelerate_nlcc * &
            (burned_fract(i,itile) * act_fpc_previous_day(i,itile) + &
            grass_mortality * act_fpc(i,itile) - &
            bare_fpc(i) / grass_sum(i) * npp_ave(i,itile) * MAX(act_fpc_min,act_fpc(i,itile)) / tau_pft)) &
            - act_fpc(i,itile)
          ELSE
            veg_inc = veg_inc + MAX(fract_small,act_fpc(i,itile) - accelerate_nlcc * &
            (burned_fract(i,itile) * act_fpc_previous_day(i,itile) + &
            act_fpc(i,itile) / tau_pft)) &
            - act_fpc(i,itile)
          ENDIF
        ENDIF
      END DO
      theta = 1._wp
      IF (bare_fpc(i) < (veg_inc + fract_small)) theta = 0._wp
      DO itile = 1,ntiles
        tau_pft = tau_Cpool_woods(i,itile) * tau_Cpool_woods_to_tau_pft
        IF (is_dynamic(i,itile) .AND. is_woody(i,itile)) THEN
          act_fpc(i,itile) = MAX(fract_small,act_fpc(i,itile) + accelerate_nlcc * &
            ((theta * pot_fpc(i,itile) - act_fpc(i,itile)) / tau_pft                &
            - (burned_fract(i,itile) + damaged_fract(i,itile) + min_dist_woody) * act_fpc_previous_day(i,itile)))
        ENDIF
      ENDDO
    ENDDO

    !-- dynamic equation for act_fpc, grass types

    DO i = 1,nidx
      DO itile = 1,ntiles
        tau_pft = tau_Cpool_woods(i,itile) * tau_Cpool_woods_to_tau_pft
        IF (is_dynamic(i,itile) .AND. .NOT. is_woody(i,itile)) THEN
          IF (bio_exist(i,itile) > 0.5_wp .AND. npp_ave(i,itile) > EPSILON(1._wp)) THEN
            act_fpc(i,itile) = MAX(fract_small,act_fpc(i,itile) - accelerate_nlcc * &
               (burned_fract(i,itile) * act_fpc_previous_day(i,itile) + &
               grass_mortality * act_fpc(i,itile) - &
               bare_fpc(i) / grass_sum(i) * npp_ave(i,itile) * MAX(act_fpc_min,act_fpc(i,itile)) / tau_pft))
          ELSE
            act_fpc(i,itile) = MAX(fract_small,act_fpc(i,itile) - accelerate_nlcc * &
               (burned_fract(i,itile) * act_fpc_previous_day(i,itile) + &
               act_fpc(i,itile) / tau_pft))
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    !-- calculate total area fraction covered by dynamic vegetation

    total_act_fpc(:) = 0._wp

    DO itile = 1,ntiles
       DO i = 1,nidx
          IF (is_dynamic(i,itile)) total_act_fpc(i) = total_act_fpc(i) + act_fpc(i,itile)
       END DO
    END DO

    !-- calculation of bare soil fraction

    bare_fpc(:) = MAX(0._wp,1._wp - total_act_fpc(:))

  END SUBROUTINE fpc_daily

! --- scale_fpc ---------------------------------------------------------------------------

  SUBROUTINE scale_fpc (routine, nidx, iblk, ntiles, is_dynamic, act_fpc, bare_fpc)
!
! Rescaling of fractional plant coverage to assure that
!  - the sum of act_fpc and bare soil is exactly one
!  - all tiles have a minimum vegetated fraction
!
! input
    CHARACTER, INTENT(in) :: routine         ! name of the calling routine
    INTEGER,  INTENT(in)  :: nidx            ! vector length
    INTEGER,  INTENT(in)  :: iblk            ! number of block
    INTEGER,  INTENT(in)  :: ntiles          ! number of tiles (different land cover type mosaics per gridcell)
    LOGICAL,  INTENT(in)  :: is_dynamic(:,:) ! PFTs taking part in competition

! inout
    REAL(wp), INTENT(inout) :: act_fpc(:,:) ! actual fpc in nlcc
    REAL(wp), INTENT(inout) :: bare_fpc(:)  ! fraction of bare ground

! local variables
    LOGICAL   :: just_return
    INTEGER   :: i, itile       ! indices
    INTEGER   :: nsparce(nidx)  ! number of PFTs with less then fract_small vegetation
    REAL(wp)  :: sum_fpc(nidx)  ! sum of the different cover fractions
    REAL(wp)  :: excess         ! extra fraction, that needs to be redistibuted
    REAL(wp)  :: scalable       ! sum of fraction of tiles that can be scaled
    REAL(wp)  :: rescale        ! factor to rescale fraction of tiles

    ! test, if scaling is necessary
    just_return = .true.
    sum_fpc(:) = bare_fpc(:)
    DO i = 1,nidx
       DO itile = 1,ntiles
          sum_fpc(i) = sum_fpc(i) + act_fpc(i,itile)
          IF (is_dynamic(i,itile) .AND. &
              act_fpc(i,itile) < fract_small) just_return = .false.
       END DO
       IF (ABS(1._wp - sum_fpc(i)) > 3._wp * EPSILON(1._wp)) just_return = .false.
    END DO

    IF (just_return) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine),'Attention: act_fpc is rescaled!')

    ! sum act_fpc and bare_fpc
    sum_fpc(:) = bare_fpc(:)
    nsparce(:) = 0

    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (is_dynamic(i,itile)) THEN
             IF (act_fpc(i,itile) > fract_small) THEN
                sum_fpc(i) = sum_fpc(i) + act_fpc(i,itile)
             ELSE
                nsparce(i) = nsparce(i) + 1
                act_fpc(i,itile) = fract_small
             END IF
          ELSE
             act_fpc(i,itile) = 0._wp
          END IF
       END DO
    END DO

    sum_fpc(:) = sum_fpc(:) + REAL(nsparce(:)) * fract_small

    ! scaling of bare_fpc
    WHERE (sum_fpc(:) <= 1._wp - EPSILON(1._wp))
       bare_fpc(:) = bare_fpc(:) + (1._wp - sum_fpc(:))
    ELSEWHERE (sum_fpc(:) >= 1._wp + EPSILON(1._wp))
       bare_fpc(:) = bare_fpc(:) / sum_fpc(:)
    END WHERE

    ! scaling of act_fpc
    DO i = 1,nidx
       excess = 0._wp
       DO itile = 1,ntiles
          IF (is_dynamic(i,itile)) THEN
             IF (sum_fpc(i) >= 1._wp + EPSILON(1._wp) .AND. act_fpc(i,itile) > fract_small * sum_fpc(i)) THEN
                act_fpc(i,itile) = act_fpc(i,itile) / sum_fpc(i)
             ELSE IF (sum_fpc(i) >= 1._wp + EPSILON(1._wp)) THEN
                excess = excess + fract_small - act_fpc(i,itile) / sum_fpc(i)
                act_fpc(i,itile) = fract_small
             END IF
          END IF
       END DO
       scalable = bare_fpc(i)
       DO itile = 1,ntiles
          IF (is_dynamic(i,itile)) THEN
             IF (act_fpc(i,itile) > 2._wp * fract_small) scalable = scalable + act_fpc(i,itile)
          END IF
       END DO
       rescale = (scalable - excess) / scalable
       bare_fpc(i) = bare_fpc(i) * rescale
       DO itile = 1,ntiles
          IF (is_dynamic(i,itile)) THEN
             IF (act_fpc(i,itile) > 2._wp * fract_small) act_fpc(i,itile) = act_fpc(i,itile) * rescale
          END IF
       END DO
    END DO

  END SUBROUTINE scale_fpc


! --- calc_cover_fract_inst ---------------------------------------------------------------------------

  SUBROUTINE calc_cover_fract_inst(nidx, ntiles, is_woody, is_dynamic, is_pasture, is_crop, &
                                   cover_fract_pot, cover_fract, cover_fract_inst)
!
! convert fractional plant cover calculated within nlcc to jsbach cover fractions
! that would result from an instantaneous establishment of crop land and pasture
! (pasture is preferrentially established on grasslands)
!
! input
    INTEGER,  INTENT(in)  :: nidx              ! vector length
    INTEGER,  INTENT(in)  :: ntiles            ! number of tiles (different land cover type mosaics per gridcell)
    LOGICAL,  INTENT(in)  :: is_woody(:,:)     ! PFTs that build wood
    LOGICAL,  INTENT(in)  :: is_dynamic(:,:)   ! PFTs taking part in competition
    LOGICAL,  INTENT(in)  :: is_pasture(:,:)   ! PFTs that are pasture
    LOGICAL,  INTENT(in)  :: is_crop(:,:)      ! PFTs that are crop
    REAL(wp), INTENT(in)  :: cover_fract_pot(:,:)  ! potential cover fractions (without land use)
    REAL(wp), INTENT(in)  :: cover_fract(:,:)  ! cover fraction within jsbach

! output 
    REAL(wp), INTENT(out) :: cover_fract_inst(:,:) ! cover fraction within jsbach resulting from an
                                                   ! instantaneous establishment of crop land and pasture
! local variables
    INTEGER   :: i, itile
    REAL(wp)  :: grass(nidx)             ! number of grass types (as real)
    REAL(wp)  :: sum_grass_fract(nidx)   ! potential fraction of grass if there was no land use
    REAL(wp)  :: sum_pasture_fract(nidx) ! fraction of the vegetated part of the grid box covered with pasture
    REAL(wp)  :: sum_crop_fract(nidx)    ! fraction of the vegetated part of the grid box covered with crops

    !-- find out number of grass types
    grass(:) = 0._wp
    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (.NOT. is_woody(i,itile) .AND. is_dynamic(i,itile)) grass(i) = grass(i) + 1._wp
       ENDDO
       grass(i) = MAX(1._wp,grass(i))
    ENDDO

    !-- sum up FPC of grasses, pastures and crops
    sum_grass_fract(:) = 0._wp
    sum_pasture_fract(:) = 0._wp
    sum_crop_fract(:) = 0._wp
    DO i = 1,nidx
      DO itile = 1,ntiles
         IF (.NOT. is_woody(i,itile) .AND. is_dynamic(i,itile)) &
            sum_grass_fract(i) = sum_grass_fract(i) + cover_fract_pot(i,itile)
         IF (is_pasture(i,itile)) sum_pasture_fract(i) = sum_pasture_fract(i) + cover_fract(i,itile)
         IF (is_crop(i,itile)) sum_crop_fract(i) = sum_crop_fract(i) + cover_fract(i,itile)
      END DO
    END DO

    !-- establishment of pastures (at the expense of grass land)
    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (is_dynamic(i,itile)) THEN
             IF (sum_pasture_fract(i) <= sum_grass_fract(i)) THEN
                ! There are enough grasslands
                IF (is_woody(i,itile)) THEN
                   cover_fract_inst(i,itile) = cover_fract_pot(i,itile)
                ELSE
                   cover_fract_inst(i,itile) = cover_fract_pot(i,itile) &
                      * (sum_grass_fract(i) - sum_pasture_fract(i) + grass(i) * fract_small) / sum_grass_fract(i)
                END IF
             ELSE
                ! There are not enough grasslands for all pastures, so that woody types are also reduced
                IF (is_woody(i,itile)) THEN
                   cover_fract_inst(i,itile) = cover_fract_pot(i,itile) &
                      * (1._wp - sum_pasture_fract(i) - 3._wp * fract_small) &
                      / (1._wp - sum_grass_fract(i) - 3._wp * fract_small)
                ELSE
                   cover_fract_inst(i,itile) = fract_small
                END IF
             END IF
          ELSE IF (is_pasture(i,itile)) THEN
             cover_fract_inst(i,itile) = cover_fract(i,itile)
          ELSE
             cover_fract_inst(i,itile) = fract_small
          END IF
       END DO
    END DO

    !-- establishment of crops (at the expense of all natural pfts)
    DO itile = 1,ntiles
       WHERE (is_crop(:,itile))
          cover_fract_inst(:,itile) = cover_fract(:,itile)
       ELSEWHERE (.NOT. is_pasture(:,itile))
          cover_fract_inst(:,itile) = cover_fract_inst(:,itile) * (1._wp - sum_pasture_fract(:) - sum_crop_fract(:)) &
             / (1._wp - sum_pasture_fract(:) - fract_small)
       END WHERE
    END DO

    CALL scale_cover_fract (nidx, ntiles, is_dynamic(:,:), cover_fract_inst(:,:))

    IF (ANY(SUM(cover_fract_inst(:,:),DIM=2) > 1._wp + REAL(ntiles) * EPSILON(1._wp)) .OR. &
        ANY(SUM(cover_fract_inst(:,:),DIM=2) < 1._wp - REAL(ntiles) * EPSILON(1._wp))) THEN
       WRITE (message_text,*) 'sum of cover_fract_inst /= 1: ', &
            MINVAL(SUM(cover_fract_inst(:,:),DIM=2)), MAXVAL(SUM(cover_fract_inst(:,:),DIM=2)), &
            MAXLOC(SUM(cover_fract_inst(:,:),DIM=2))
       CALL finish ('calc_cover_fract_inst', message_text)
    END IF
    DO itile = 1,ntiles
       IF (ANY(cover_fract_inst(:,itile) < fract_small)) THEN
          WRITE (message_text,*) 'cover_fract_inst too small: ', MINVAL(cover_fract_inst(:,itile)), &
               MINLOC(cover_fract_inst(:,itile))
          CALL finish ('calc_cover_fract_inst', message_text)
       END IF
    END DO

  END SUBROUTINE calc_cover_fract_inst


!------------------------------------------------------------------------------

SUBROUTINE scale_cover_fract (nidx, ntiles, is_dynamic, cover_fract)
!
! !DESCRIPTION:
!
! Rescaling of cover fractions to assure that
!  - the sum of cover fractions is one
!  - all tiles have at least a minimum vegetated fraction
!
!------------------------------------------------------------------------------
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(in)  :: nidx               ! vector length
    INTEGER,  INTENT(in)  :: ntiles             ! number of tiles
    LOGICAL,  INTENT(in)  :: is_dynamic(:,:)    ! mask for land surface treated by jsbach
!
! !IN- and OUTPUT PARAMETERS:
! 
    REAL(wp), INTENT(inout) :: cover_fract(:,:) ! vegetated fraction

! !LOCAL VARIABLES:
!
    INTEGER   :: i, iter
    INTEGER   :: niter             ! number of iterations needed
    INTEGER   :: nsparce(nidx)     ! number of PFTs with a vegetated fraction of less then fract_small
    REAL(wp)  :: sum_fract(nidx)   ! sum of all cover fractions
    REAL(wp)  :: excluded_fract(nidx) ! sum of all cover fractions

!------------------------------------------------------------------------------

    ! Make sure, crop and pasture have a cover fraction of at least fract_small

    WHERE (.NOT. is_dynamic(:,:))
       cover_fract(:,:) = MAX(fract_small,cover_fract(:,:))
    END WHERE

    ! The more tiles we have, the more iterations are needed. For binary identical results
    ! the number of iterations must not change.
    niter = ntiles
    DO iter = 1, niter

       sum_fract(:) = 0._wp
       excluded_fract(:) = 0._wp
       nsparce(:) = 0

       DO i = 1,ntiles
          WHERE (cover_fract(:,i) > fract_small .AND. is_dynamic(:,i))
             sum_fract(:) = sum_fract(:) + cover_fract(:,i)
          ELSEWHERE (is_dynamic(:,i))
             nsparce(:) = nsparce(:) + 1
          ELSEWHERE
             excluded_fract(:) = excluded_fract(:) + cover_fract(:,i)
          END WHERE
       END DO
       DO i = 1,ntiles
          WHERE (cover_fract(:,i) > fract_small .AND. is_dynamic(:,i))
             cover_fract(:,i) = cover_fract(:,i) * (1._wp - excluded_fract(:) - REAL(nsparce,wp)*fract_small) / sum_fract(:)
          ELSEWHERE (is_dynamic(:,i))
             cover_fract(:,i) = fract_small
          ELSEWHERE
             cover_fract(:,i) = MAX(fract_small,cover_fract(:,i))
          END WHERE
       END DO

    END DO

    IF (ANY(cover_fract(:,:) < fract_small)) THEN
       WRITE(message_text,*) 'cover_fract still smaller ', fract_small, ' after ', ntiles, ' iterations:', &
            MINVAL(MERGE(cover_fract(:,:), 1._wp, cover_fract > 0._wp))
       CALL message ('scale_cover_fract', message_text)

       ! sometimes cover_fract remains slightly smaller than fract_small for numerical reasons
       WHERE (cover_fract(:,:) < fract_small .AND. cover_fract(:,:) >= fract_small - REAL(ntiles,wp)*EPSILON(1._wp))
          cover_fract(:,:) = fract_small
       END WHERE
       IF (ANY(cover_fract(:,:) < fract_small)) THEN
          WRITE(message_text,*) 'cover_fract still smaller ', fract_small, ' after ', ntiles, ' iterations:', &
               MINVAL(MERGE(cover_fract(:,:), 1._wp, cover_fract > 0._wp))
          CALL finish ('scale_cover_fract', message_text)
       END IF
    END IF

  END SUBROUTINE scale_cover_fract

#endif
END MODULE mo_nlcc_process
