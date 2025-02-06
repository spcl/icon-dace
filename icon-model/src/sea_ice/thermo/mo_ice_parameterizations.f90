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

! Provide an implementation of various sea ice parametrizations.

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_ice_parameterizations

  USE mo_kind,                ONLY: wp
  USE mo_run_config,          ONLY: dtime
  USE mo_dynamics_config,     ONLY: nold
  USE mo_exception,           ONLY: finish
  USE mo_impl_constants,      ONLY: max_char_length
  USE mo_util_dbg_prnt,       ONLY: dbg_print

  USE mo_model_domain,        ONLY: t_patch
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range

  USE mo_physical_constants,  ONLY: rhoi, rhos, rho_ref, clw,             &
    &                               alb_sno_vis, alb_sno_nir, alb_ice_vis, alb_ice_nir
  USE mo_sea_ice_nml,         ONLY: i_ice_albedo, i_Qio_type, i_ice_therm, albi, albim, albsm, albs, Cd_io, Ch_io
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_sea_ice_types,       ONLY: t_sea_ice
  USE mo_fortran_tools,       ONLY: set_acc_host_or_device


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_ice_albedo
  PUBLIC :: oce_ice_heatflx
  PUBLIC :: ice_draft_and_flooding
  PUBLIC :: ice_cut_off

  CHARACTER(len=12)           :: str_module    = 'IceShared'  ! Output of module for 1 line debug

CONTAINS

  !-------------------------------------------------------------------------------
  !>
  !! ! set_ice_albedo: set ice albedo
  !!
  SUBROUTINE set_ice_albedo(i_startidx_c, i_endidx_c, nbdim, kice, Tsurf, hi, hs, &
      & albvisdir, albvisdif, albnirdir, albnirdif, lacc)

    INTEGER, INTENT(IN)  :: i_startidx_c, i_endidx_c, nbdim, kice
    REAL(wp),INTENT(IN)  :: Tsurf(nbdim,kice)
    REAL(wp),INTENT(IN)  :: hi   (nbdim,kice)
    REAL(wp),INTENT(IN)  :: hs   (nbdim,kice)
    REAL(wp),INTENT(OUT) :: albvisdir  (nbdim,kice)
    REAL(wp),INTENT(OUT) :: albvisdif  (nbdim,kice)
    REAL(wp),INTENT(OUT) :: albnirdir  (nbdim,kice)
    REAL(wp),INTENT(OUT) :: albnirdif  (nbdim,kice)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc


    !Local variables
    REAL(wp), PARAMETER :: albtrans   = 0.5_wp
    REAL(wp)            :: albflag, frac_snow
    INTEGER             :: jc,k
    LOGICAL :: lzacc
    !-------------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    SELECT CASE (i_ice_albedo)
    CASE (1)
      ! This is Uwe's albedo expression from the old budget function
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO k=1,kice
        !$ACC LOOP GANG VECTOR PRIVATE(albflag)
        DO jc = i_startidx_c,i_endidx_c

          albflag =  1.0_wp/ ( 1.0_wp+albtrans * (Tsurf(jc,k))**2 )

          IF ( hi(jc,k) > 0._wp ) THEN
            IF ( hs(jc,k) > 1.e-2_wp ) THEN
              albvisdir(jc,k) =  albflag * albsm + (1.0_wp-albflag) * albs
            ELSE
              albvisdir(jc,k) =  albflag * albim + (1.0_wp-albflag) * albi
            ENDIF
          ELSE
            albvisdir(jc,k) = 0._wp
          ENDIF

        ENDDO
      ENDDO
      !$ACC END PARALLEL

      ! all albedos are the same
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO k=1,kice
        !$ACC LOOP GANG VECTOR
        DO jc = 1,nbdim
          albvisdif(jc,k) = albvisdir(jc,k)
          albnirdir(jc,k) = albvisdir(jc,k)
          albnirdif(jc,k) = albvisdir(jc,k)
        END DO
      END DO
      !$ACC END PARALLEL

    CASE (2)
      ! This is the CCSM 3 albedo scheme
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
!PREVENT_INCONSISTENT_IFORT_FMA
      DO k=1,kice
        !$ACC LOOP GANG VECTOR PRIVATE(frac_snow)
        DO jc = i_startidx_c,i_endidx_c
          frac_snow = hs(jc,k)/( hs(jc,k)+0.02_wp )
          IF ( Tsurf(jc,k) > -1._wp ) THEN
            albvisdir(jc,k) = frac_snow*( alb_sno_vis - 0.100_wp*(Tsurf(jc,k)+1._wp) ) &
              &     + (1._wp-frac_snow)*( alb_ice_vis - 0.075_wp*(Tsurf(jc,k)+1._wp) )
            albnirdir(jc,k) = frac_snow*( alb_sno_nir - 0.150_wp*(Tsurf(jc,k)+1._wp) ) &
              &     + (1._wp-frac_snow)*( alb_ice_nir - 0.075_wp*(Tsurf(jc,k)+1._wp) )
          ELSE
            albvisdir(jc,k) = frac_snow*alb_sno_vis + (1._wp-frac_snow)*alb_ice_vis
            albnirdir(jc,k) = frac_snow*alb_sno_nir + (1._wp-frac_snow)*alb_ice_nir
          ENDIF
        ENDDO
      ENDDO
      !$ACC END PARALLEL

      ! diffuse and direct albedos are the same
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO k=1,kice
        !$ACC LOOP GANG VECTOR
        DO jc = 1, nbdim
          albvisdif(jc,k) = albvisdir(jc,k)
          albnirdif(jc,k) = albnirdir(jc,k)
        END DO
      END DO
      !$ACC END PARALLEL

    END SELECT

  END SUBROUTINE set_ice_albedo

  !-------------------------------------
  !
  ! oce_ice_heatflx
  !
  ! Calculates the heat flux from the uppermost water layer into the ice.
  !
  ! Currently (as in growth.f90): all energy available in upper ocean grid cell
  ! is supplied to the ice and the upper ocean temperature is held at the
  ! freezing point. This is not very physical.
  !
  ! Positive flux upwards.

SUBROUTINE oce_ice_heatflx (p_patch, p_os, ice, lacc)
    TYPE(t_patch), INTENT(IN), TARGET       :: p_patch
    TYPE(t_hydro_ocean_state), INTENT(IN)   :: p_os
    TYPE(t_sea_ice),           INTENT(INOUT):: ice
    LOGICAL, INTENT(IN), OPTIONAL           :: lacc

    ! Local
    INTEGER :: jb, jc, i_startidx_c, i_endidx_c, k
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp) :: u_star
    REAL(wp), POINTER  :: sst(:,:), u_diag(:,:), v_diag(:,:)
    LOGICAL :: lzacc

    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ice_parameterizations:oce_ice_heatflx'

    CALL set_acc_host_or_device(lzacc, lacc)

    all_cells   => p_patch%cells%all
    sst         => p_os%p_prog(nold(1))%tracer(:,1,:,1)
    u_diag      => p_os%p_diag%u(:,1,:)
    v_diag      => p_os%p_diag%v(:,1,:)

!    ! Initialization
!    ice%zheatOceI(:,:,:) = 0.0_wp ! initialized in ice_zero

    !---------DEBUG DIAGNOSTICS-------------------------------------------
     CALL dbg_print('O-I-HeatFlx: SST'       ,sst           ,str_module,4, in_subset=p_patch%cells%owned)
     CALL dbg_print('O-I-HeatFlx: zUnderIce' ,ice%zUnderIce ,str_module,4, in_subset=p_patch%cells%owned)
!    CALL dbg_print('O-I-HeatFlx: Tfw'       ,ice%Tfw       ,str_module,4, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

    IF (i_ice_therm == 3) RETURN ! if analytical atmospheric fluxes, NO heatflux from the "ocean" is applied

    SELECT CASE ( i_Qio_type )

    CASE (0)
!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, jc) SCHEDULE(dynamic)
      DO jb = 1,p_patch%nblks_c
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
        DO jc = i_startidx_c,i_endidx_c
          ice%zHeatOceI(jc,:,jb) = 0.0_wp
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDDO

    CASE (1)
!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, jc) SCHEDULE(dynamic)
      DO jb = 1,p_patch%nblks_c
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
        DO jc = i_startidx_c,i_endidx_c
          IF (ice%hi(jc,1,jb) > 0._wp) THEN
            ! energy of warm water below ice covered part of grid area only is used for melting
            ice%zHeatOceI(jc,:,jb) = ( sst(jc,jb) - ice%Tfw(jc,jb) ) * ice%zUnderIce(jc,jb) * clw*rho_ref/dtime
          ENDIF
        ENDDO
        !$ACC END PARALLEL LOOP
      END DO
!ICON_OMP_END_PARALLEL_DO

   CASE(2)
!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, jc) SCHEDULE(dynamic)
    DO jb = 1,p_patch%nblks_c
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
      DO jc = i_startidx_c,i_endidx_c
        IF (ice%hi(jc,1,jb) > 0._wp) THEN
            ! melting energy depends on velocity difference between water and ice, bulk formulation
            u_star = SQRT(Cd_io*( (u_diag(jc,jb)-ice%u(jc,jb))**2 + &
                         &        (v_diag(jc,jb)-ice%v(jc,jb))**2 ) )
            ice%zHeatOceI(jc,:,jb) = ( sst(jc,jb) - ice%Tfw(jc,jb) ) *rho_ref*clw*Ch_io*u_star
        ENDIF
      ENDDO
      !$ACC END PARALLEL LOOP
    END DO
!ICON_OMP_END_PARALLEL_DO

    CASE (3)
!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, jc, k) SCHEDULE(dynamic)
      DO jb = 1,p_patch%nblks_c
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
        DO jc = i_startidx_c,i_endidx_c
          DO k=1,ice%kice
          IF (ice%hi(jc,k,jb) > 0._wp) THEN
            ! ALL energy of warm water over the whole grid area is used for melting ice - divide by concentration and distribute
            ! SLO/EO 2014-04-11 - this is wrong, must be accompanied by correction elsewhere, since open part of water is still losing heat
            ice%zHeatOceI(jc,k,jb) = ( sst(jc,jb) - ice%Tfw(jc,jb) ) * ice%zUnderIce(jc,jb) * clw*rho_ref / &
                                   & (dtime*ice%concSum(jc,jb)) * ( ice%conc(jc,k,jb)/ice%concSum(jc,jb) )
          ELSE  ! hi<=0
            ! Correction needed if ice was melted in the last timestep
            ice%zHeatOceI(jc,k,jb) = 0.0_wp
          ENDIF
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP
      END DO
!ICON_OMP_END_PARALLEL_DO

    CASE DEFAULT
      CALL finish(TRIM(routine), 'Invalid i_Qio_type')

    END SELECT

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('O-I-HeatFlx: zHeatOceI ' ,ice%zHeatOceI , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('O-I-HeatFlx: conc      ' ,ice%conc      , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('O-I-HeatFlx: concSum   ' ,ice%concSum   , str_module, 4, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE oce_ice_heatflx


  !-------------------------------------------------------------------------------
  !>
  !! ! ice_draft_and_flooding: calculate draft and conversion of snow below waterlevel to ice
  !!
  !! Currently not quite physical: Snow is pushed together to form new ice, hence snow thickness
  !! decreases more than ice thickness by rhoi/rhos ( analogue to old growth.f90 sea-ice model )
  !! Salt content of snow ice is equal to that of normal ice, salt is removed from the ocean
  !!
  SUBROUTINE ice_draft_and_flooding (p_patch, ice, lacc)

    TYPE(t_patch),      INTENT(IN), TARGET  :: p_patch
    TYPE (t_sea_ice),   INTENT(INOUT)       :: ice
    LOGICAL, INTENT(IN), OPTIONAL           :: lacc

    ! Local
    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER                       :: jb, k, jc, i_startidx_c, i_endidx_c
    REAL(wp)                      :: hi_from_flood ! Increase in ice thickness due to flooding      [m]
    LOGICAL                       :: lzacc
    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    all_cells            => p_patch%cells%all

!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, k, jc, hi_from_flood) SCHEDULE(dynamic)
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) IF(lzacc)
      DO k=1,ice%kice
        DO jc = i_startidx_c,i_endidx_c

            ! Archimedes principle
            ice%draft       (jc,k,jb) = ( rhoi*ice%hi(jc,k,jb) + rhos*ice%hs(jc,k,jb) )/rho_ref
            ! Increase in ice thickness due to flooding
            hi_from_flood             = ice%draft(jc,k,jb) - MIN( ice%draft(jc,k,jb), ice%hi(jc,k,jb) )
            ! Thickness of snow that is converted into ice (in units of snow thickness/density)
            ice%snow_to_ice (jc,k,jb) = hi_from_flood*rhoi/rhos

            ! update hi and hs
            ice%hs          (jc,k,jb) = ice%hs(jc,k,jb) - ice%snow_to_ice(jc,k,jb)
            ice%hi          (jc,k,jb) = ice%hi(jc,k,jb) + hi_from_flood

        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!ICON_OMP_END_PARALLEL_DO

!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, jc) SCHEDULE(dynamic)
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
      DO jc = i_startidx_c,i_endidx_c
        ice%draftave(jc,jb) = sum(ice%draft(jc,:,jb) * ice%conc(jc,:,jb))
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!ICON_OMP_END_PARALLEL_DO


    ! do not update here
!    ice%zUnderIce   (:,:) = v_base%del_zlev_m(1) + p_os%p_prog(nold(1))%h(:,:) - ice%draftave(:,:)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('Flooding: hi aft flooding'  , ice%hi        ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('Flooding: hs aft flooding'  , ice%hs        ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('Flooding: snow_to_ice'      , ice%snow_to_ice, str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('Flooding: draft aft flood'  , ice%draft     , str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('Flooding: draftave aft'     , ice%draftave  , str_module, 3, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE ice_draft_and_flooding

  !-------------------------------------------------------------------------------
  !>
  !! ! ice_cut_off: Fix overshoots in concentration that appear due to ice convergence
  !! NB: currently ONLY for the one-ice-class case, kice=1
  !!
  SUBROUTINE ice_cut_off (p_patch, p_ice, lacc)

    TYPE(t_patch),TARGET,   INTENT(IN)    :: p_patch
    TYPE(t_sea_ice),        INTENT(INOUT) :: p_ice
    LOGICAL, INTENT(IN), OPTIONAL         :: lacc

    ! Local
    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER :: jb, jc, i_startidx_c, i_endidx_c
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    all_cells            => p_patch%cells%all

    ! Fix overshoots in concentration that can appear due to ice convergence in areas with conc ~ 1
!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, i_startidx_c, i_endidx_c)
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
      DO jc = i_startidx_c,i_endidx_c
        IF (( p_ice%conc(jc,1,jb) > 1._wp )) THEN
          p_ice%conc  (jc,1,jb) = 1._wp
          ! New ice and snow thickness
          p_ice%hi   (jc,1,jb) = p_ice%vol (jc,1,jb)/( p_ice%conc(jc,1,jb)*p_patch%cells%area(jc,jb) )
          p_ice%hs   (jc,1,jb) = p_ice%vols(jc,1,jb)/( p_ice%conc(jc,1,jb)*p_patch%cells%area(jc,jb) )
        END IF
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!ICON_OMP_END_PARALLEL_DO

    ! Fix undershoots - ONLY for the one-ice-class case
    ! Quick fix, should be reformulated to occur at the advection stage
!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, i_startidx_c, i_endidx_c)
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
      DO jc = i_startidx_c,i_endidx_c
        IF (( p_ice%conc(jc,1,jb) < TINY(1._wp) ) .OR. ( p_ice%hi(jc,1,jb) < TINY(1._wp) )) THEN
          p_ice%conc  (jc,1,jb) = 0._wp
          p_ice%hi    (jc,1,jb) = 0._wp
          p_ice%hs    (jc,1,jb) = 0._wp
          p_ice%vol   (jc,1,jb) = 0._wp
          p_ice%vols  (jc,1,jb) = 0._wp
        END IF
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!ICON_OMP_END_PARALLEL_DO

!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, i_startidx_c, i_endidx_c)
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
      DO jc = i_startidx_c,i_endidx_c
        p_ice%concSum(jc,jb) = SUM(p_ice%conc(jc,:,jb))
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE ice_cut_off

END MODULE mo_ice_parameterizations
