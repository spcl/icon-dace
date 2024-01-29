!
! mo_art_surface_value
! This module provides the calculation of surface values of a tracer
! Based on: Rieger et al. (2015): ICON-ART 1.0 - a new online-coupled model system
!                                 from the global to regional scale
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

MODULE mo_art_surface_value
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_parallel_config,               ONLY: nproma
  USE mo_impl_constants,                ONLY: min_rlcell_int
  USE mo_impl_constants_grf,            ONLY: grf_bdywidth_c
  USE mo_loopindices,                   ONLY: get_indices_c
  USE mo_physical_constants,            ONLY: grav
  USE mo_nonhydro_types,                ONLY: t_nh_metrics, t_nh_diag, t_nh_prog
  USE mo_nwp_phy_types,                 ONLY: t_nwp_phy_diag
  USE mo_fortran_tools,                 ONLY: assert_acc_device_only

! ART
  USE mo_art_config,                    ONLY: art_config
  USE mo_art_diag_types,                ONLY: t_art_diag, art_diag_tracer_index
  USE mo_art_modes_linked_list,         ONLY: p_mode_state, t_mode
  USE mo_art_modes,                     ONLY: t_fields_2mom, t_fields_radio
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_emiss_types,               ONLY: t_art_emiss_type_container
  USE mo_art_impl_constants,            ONLY: IART_ACC_DRYDEPO

  IMPLICIT NONE

  PRIVATE

  PUBLIC   :: art_surface_value

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_surface_value( p_patch, p_prog, p_metrics, p_diag, prm_diag,      &
  &                           art_diag, p_rho, dt, jg, jb, vdep, surface_value,  &
  &                           cfactor_opt, lacc )
!<
! SUBROUTINE art_surface_value
! This subroutine calculates an artificial surface value
! Part of Module: mo_art_surface_value
! Based on: Rieger et al. (2015): ICON-ART 1.0 - a new online-coupled model system
!                                 from the global to regional scale
! Author: Jochen Foerstner, DWD
! Initial Release: 2013-12-03
! Modifications:
! 2016-11-30: Michael Weimer, KIT
! - included defintion of surface value for emissions
! YYYY-MM-DD: <name>, <instituttion>
! - ...
!>
  TYPE(t_patch), TARGET, INTENT(in) ::  &  !< patch on which computation
    &  p_patch                             !< is performed

  TYPE(t_nh_prog),  INTENT(in)      :: p_prog                !< the prog vars

  TYPE(t_nh_metrics), INTENT(in)    :: p_metrics
  TYPE(t_nh_diag), INTENT(in)       :: p_diag                !< the diag vars
  TYPE(t_nwp_phy_diag), INTENT(in)  :: prm_diag              !< atm phys vars

  TYPE(t_art_diag), POINTER         :: art_diag              !< Pointer to ART diagnostic fields
  REAL(wp), INTENT(in)              :: p_rho(:,:,:)          !< air density
  REAL(wp), INTENT(in)              :: dt                    !< time step

  INTEGER, INTENT(in)               :: jg, jb

  REAL(wp), INTENT(in)              :: vdep(:,:,:)           !< deposition velocity

  REAL(wp), INTENT(out)             :: surface_value(:,:,:)  !< surface value according to
                                                             !    transfer coeff.

  REAL(wp), INTENT(in), OPTIONAL    :: cfactor_opt(:,:)      !< optional conversion factor

  LOGICAL, OPTIONAL, INTENT(IN)     :: lacc

  !Local variables
  ! ---------------------------------------

  TYPE(t_mode), POINTER       :: &
    &  this_mode                    !< Pointer to current mode

  TYPE(t_art_emiss_type_container) :: &
    &  emiss

  REAL(wp) :: tch_vabs(nproma), &    ! Units: m/s
    &         ref_height(nproma)     ! Units: m
                                     ! dim:(nproma)

  REAL(wp) :: cfactor(nproma)        ! Units: -
                                     ! dim:(nproma)


  REAL(wp) :: vdep_ratio, numerator, denominator,  conv_fac, &
    &         conc_atmosphere, surface_value_tmp, depo_update_fac

  INTEGER  :: nlev,                                              & !<number of full and half levels
    &         i_nchdom, i_rlstart, i_rlend, i_startblk, i_endblk

  INTEGER  :: num_radioact,                  &
    &         jsp, jc, i_startidx, i_endidx, &
    &         itracer, i, j,                 &
    &         ierror

  LOGICAL  :: lmass_2mom(art_config(jg)%nturb_tracer)
  INTEGER  :: iradioact(art_config(jg)%nturb_tracer)

  INTEGER :: idx_diag               !< index of tracer in diagnostics container

  CALL assert_acc_device_only('art_surface_value', lacc)

  !!
  !!-------------BEGIN SUBROUTINE -----------------------------------------------
  !!

  !$ACC DATA CREATE(cfactor, ref_height, tch_vabs) &
  !$ACC   PRESENT(art_diag, prm_diag, p_diag, p_metrics, p_prog, p_rho, surface_value, vdep)

  IF ( PRESENT( cfactor_opt ) ) THEN
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = 1, nproma
      cfactor(jc) = cfactor_opt(jc,jb)
    END DO
    !$ACC END PARALLEL
  ELSE
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = 1, nproma
      cfactor(jc) = 1.0_wp
    END DO
    !$ACC END PARALLEL
  END IF

  lmass_2mom(:) = .FALSE.
  iradioact(:) = -1

  num_radioact = 0

  IF (art_config(jg)%lart_aerosol) THEN

    this_mode => p_mode_state(jg)%p_mode_list%p%first_mode

    DO WHILE(ASSOCIATED(this_mode))
      ! Select type of mode
      SELECT TYPE (fields=>this_mode%fields)

        TYPE IS (t_fields_2mom)
        DO i = 1, fields%ntr-1
          DO jsp = 1, art_config(jg)%nturb_tracer
            itracer = p_prog%turb_tracer(jb,jsp)%idx_tracer
            IF (itracer == fields%itr3(i)) THEN
              lmass_2mom(jsp) = .TRUE.
            ENDIF
          ENDDO
        ENDDO

        TYPE IS (t_fields_radio)
        num_radioact = num_radioact+1
        DO jsp = 1, art_config(jg)%nturb_tracer
          itracer = p_prog%turb_tracer(jb,jsp)%idx_tracer
          IF (itracer == fields%itr) iradioact(jsp) = num_radioact
        END DO

      END SELECT
      this_mode => this_mode%next_mode
    END DO

  END IF

  nlev      = p_patch%nlev        !< Number of vertical full levels

  !Get all cell enitities, except halos
  i_nchdom   = MAX(1,p_patch%n_childdom)
  i_rlstart  = grf_bdywidth_c+1
  i_rlend    = min_rlcell_int
  i_startblk = p_patch%cells%start_blk(i_rlstart,1)
  i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

  CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,          &
    &                 i_startidx, i_endidx, i_rlstart, i_rlend )

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
  !$ACC LOOP GANG VECTOR
  DO jc = i_startidx, i_endidx
    tch_vabs(jc)   = prm_diag%tch(jc,jb) * SQRT( p_diag%u(jc,nlev,jb)**2 + p_diag%v(jc,nlev,jb)**2)
    ref_height(jc) = prm_diag%gz0(jc,jb) / grav
    ref_height(jc) = 10.0_wp * MIN( MAX( ref_height(jc), 0.01_wp ), 0.1_wp )
    ref_height(jc) = 2.0_wp * ref_height(jc) * p_metrics%inv_ddqz_z_full(jc,nlev,jb)
  ENDDO !jc
  !$ACC END PARALLEL

  ! DIAGNOSTIC: acc_drydepo (accumulated dry deposition of art-tracer)
  DO jsp = 1, art_config(jg)%nturb_tracer
    IF ( iradioact(jsp) == -1 ) THEN
      ! get tracer indices in tracer- and in diagnostics-container for diagnostic 'acc_drydepo'
      idx_diag = art_diag_tracer_index(IART_ACC_DRYDEPO, p_prog%turb_tracer(jb,jsp)%idx_tracer)

      IF ( lmass_2mom(jsp) ) THEN
        conv_fac = 1.e-9_wp
      ELSE
        conv_fac = 1.0_wp
      ENDIF

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(vdep_ratio, conc_atmosphere, numerator, denominator, surface_value_tmp, depo_update_fac)
      DO jc = i_startidx, i_endidx
        vdep_ratio  = vdep(jc,jb,jsp) / tch_vabs(jc)
        ! save concentration in lowest model layer
        conc_atmosphere = cfactor(jc) * p_prog%turb_tracer(jb,jsp)%ptr(jc,nlev)
        numerator   = conc_atmosphere * ( 1.0_wp - vdep_ratio*ref_height(jc) )
        denominator = 1.0_wp + vdep_ratio * ( 1.0_wp-ref_height(jc) )
        surface_value_tmp = numerator / denominator
        depo_update_fac = dt * p_rho(jc,nlev,jb) * tch_vabs(jc)
        surface_value(jc,jb,jsp) = MAX( surface_value_tmp, conc_atmosphere - conc_atmosphere / depo_update_fac )
        ! dry deposition shall be positive, thus substract negative surface flux for accumulation
        IF ( idx_diag > 0 ) THEN
           art_diag%acc_drydepo(jc,jb,idx_diag) =                     &
                &    art_diag%acc_drydepo(jc,jb,idx_diag)             &
                &    + depo_update_fac * conv_fac                     &
                &    * ( conc_atmosphere - surface_value(jc,jb,jsp) )
        END IF
      ENDDO !jc
      !$ACC END PARALLEL

    !NOTE: for the future: the following ELSE block may be replaced by the upper one by simply adding the
    !      radioactive tracer names to the acc_drydepo-block in the diagnostics.xml file
    ELSE
      ! radioactive tracer
      DO jc = i_startidx, i_endidx
        vdep_ratio  = vdep(jc,jb,jsp) / tch_vabs(jc)
        ! save concentration in lowest model layer
        conc_atmosphere = cfactor(jc) * p_prog%turb_tracer(jb,jsp)%ptr(jc,nlev)
        numerator   = conc_atmosphere * ( 1.0_wp - vdep_ratio*ref_height(jc) )
        denominator = 1.0_wp + vdep_ratio * ( 1.0_wp-ref_height(jc) )
        surface_value_tmp = numerator / denominator
        depo_update_fac = dt * p_rho(jc,nlev,jb) * tch_vabs(jc)
        surface_value(jc,jb,jsp) = MAX( surface_value_tmp, conc_atmosphere - conc_atmosphere / depo_update_fac )
        ! Accumulated dry deposition [Bq/m2]
        art_diag%radioact(iradioact(jsp))%drydepo(jc,jb) =       &
          &    art_diag%radioact(iradioact(jsp))%drydepo(jc,jb)  &
          &    + depo_update_fac                                 &
          &    * ( conc_atmosphere - surface_value(jc,jb,jsp) )
      ENDDO !jc
    END IF


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Treatment of emissions in turbulence scheme
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IF (art_config(jg)%lart_emiss_turbdiff) THEN
      ! in case of emissions to be added to the turbulence scheme the surface
      ! value is overwritten and treated as flux (see mo_art_turbdiff_interface)
      ! Please be aware that the above surface value (especially for
      ! aerosols!!!) cannot be combined with the below calculation. For
      ! aerosols, use lart_emiss_turbdiff = .false. in art_nml instead
      surface_value(:,:,jsp) = 0.0_wp
      itracer = p_prog%turb_tracer(jb,jsp)%idx_tracer

      CALL p_art_data(jg)%emiss%get(itracer,emiss,ierror)

      IF (ierror == 0) THEN
        IF (emiss%bioonl%idx > 0) THEN
          DO jc = i_startidx, i_endidx
            surface_value(jc,jb,jsp) = surface_value(jc,jb,jsp)      &
                          &            - emiss%bioonl%cbio(jc,jb)    &
                          &            * emiss%bioonl%scaling_factor
          END DO
        END IF


        IF ((emiss%num_types_prescribed == 0) .AND. (emiss%bioonl%idx <= 0)   &
            &  .AND. (emiss%std%mode == 1)) THEN
          DO jc = i_startidx, i_endidx
            surface_value(jc,jb,jsp) = surface_value(jc,jb,jsp) - emiss%std%val
          END DO
        END IF

        IF (emiss%num_types_prescribed > 0) THEN
          DO j = 1,emiss%num_types_prescribed
            DO jc = i_startidx, i_endidx
              surface_value(jc,jb,jsp) = surface_value(jc,jb,jsp)            &
                    &                    - emiss%types(j)%emiss_2d(jc,jb,2)  &
                    &                    * emiss%types(j)%scaling_factor
            END DO
          END DO
        END IF

      END IF
    END IF
  ENDDO !jsp

  !$ACC WAIT
  !$ACC END DATA

END SUBROUTINE art_surface_value
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_surface_value
