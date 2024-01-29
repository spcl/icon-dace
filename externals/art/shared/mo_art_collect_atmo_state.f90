!
! mo_art_collect_atmo_state
! This module collects the atmospheric state and physics dependent
! variables
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

MODULE mo_art_collect_atmo_state
  USE mo_kind,                          ONLY: wp
  USE mo_grid_config,                   ONLY: nroot
  USE mo_var,                           ONLY: t_var
  USE mo_var_list,                      ONLY: t_var_list_ptr, find_tracer_by_index
  USE mo_var_groups,                    ONLY: var_groups_dyn
  USE mo_impl_constants,                ONLY: min_rlcell_int, inwp, iaes
  USE mo_impl_constants_grf,            ONLY: grf_bdywidth_c
  USE mo_physical_constants,            ONLY: grav, amo3, amd, amch4
  USE mo_nonhydro_types,                ONLY: t_nh_state, t_nh_prog
  USE mo_exception,                     ONLY: message,finish,message_text
  USE mo_fortran_tools,                 ONLY: set_acc_host_or_device

  USE mo_model_domain,                  ONLY: p_patch
  USE mo_parallel_config,               ONLY: nproma
  USE mo_art_config,                    ONLY: art_config

  USE mo_ext_data_types,                ONLY: t_external_data
  USE mo_nwp_phy_types,                 ONLY: t_nwp_phy_diag

  USE mo_atm_phy_nwp_config,            ONLY: atm_phy_nwp_config
  USE mo_newcld_optics,                 ONLY: rrtm_reimax => reimax,  &
                                          &   rrtm_reimin => reimin,  &
                                          &   rrtm_relmin => relmin,  &
                                          &   rrtm_relmax => relmax
  USE mo_radiation_cloud_optics,        ONLY: rte_rrtmgp_reimax => reimax,  &
                                          &   rte_rrtmgp_reimin => reimin,  &
                                          &   rte_rrtmgp_relmin => relmin,  &
                                          &   rte_rrtmgp_relmax => relmax
#ifndef __NO_AES__
  USE mo_aes_phy_memory,                ONLY: prm_field
#endif
  USE mo_aes_phy_config,                ONLY: aes_phy_config
  USE mo_aes_rad_config,                ONLY: aes_rad_config
  USE mo_bc_ozone,                      ONLY: read_bc_ozone, ext_ozone
  USE mo_o3_util,                       ONLY: o3_pl2ml, o3_timeint
  USE mo_bc_greenhouse_gases,           ONLY: ghg_ch4mmr
#ifndef __NO_JSBACH__
  USE mo_jsb_interface,                 ONLY: jsbach_get_var
#endif

  USE mtime,                            ONLY: datetime

  USE mo_art_config,                    ONLY: art_config

  !ART
  USE mo_art_modes_linked_list,         ONLY: p_mode_state, t_mode
  USE mo_art_modes,                     ONLY: t_fields_2mom
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_atmo_data,                 ONLY: t_art_atmo
  USE mo_art_chem_data,                 ONLY: t_art_chem_indices
  USE mo_art_sza,                       ONLY: art_calc_sza
  USE mo_art_diagnostics,               ONLY: art_calc_tropopause
  USE mo_art_wrapper_routines,          ONLY: art_get_indices_c
  USE mo_art_read_initial,              ONLY: art_init_aero, art_init_chem
#ifdef __ART_GPL
  USE mo_art_init_full_chemistry,       ONLY: art_init_full_chemistry
#endif
  USE mo_art_init_chemtracer,           ONLY: art_init_chemtracer, art_init_chemtracer_ext, &
                                          &   art_set_init_tracer_trop_strat
  USE mo_art_init_linoz_gems,           ONLY: calc_o3_gems_linoz
  USE mo_art_chem_utils,                ONLY: art_convert_tracers_vmr_mmr

#include "add_var_acc_macro.inc"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_collect_atmo_state'

  PUBLIC :: art_collect_atmo_state_nwp, art_update_atmo_state_nwp
  PUBLIC :: art_collect_atmo_state_aes, art_update_atmo_state_aes
  PUBLIC :: td_and_t2rh
  PUBLIC :: art_deallocate_atmo_state
  PUBLIC :: art_collect_indices
  PUBLIC :: art_init_tracer_values_nwp, art_init_tracer_values_aes
  PUBLIC :: art_collect_radiation_properties

CONTAINS

SUBROUTINE art_collect_indices(jg)
!<
! SUBROUTINE art_collect_indices
! This routine collects the dimensions from the host model, and other indices
! because they have to be available already at an early stage of the model
! Part of Module: mo_art_collect_atmo_state
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-05
! Modifications:
!>
  IMPLICIT NONE

  INTEGER, INTENT(in)     ::  &
    &  jg                   !< patch id
  ! local variables
  INTEGER ::                  &
    &  i_startidx, i_endidx,  &
    &  jc, jb,                &
    &  i_rlstart, i_rlend,    & !< temporal indices
    &  i_nchdom
  TYPE(t_art_atmo),POINTER    :: &
     &  art_atmo           !< Pointer to ART atmospheric fields

  i_nchdom   = MAX(1,p_patch(jg)%n_childdom)
  i_rlstart = grf_bdywidth_c+1
  i_rlend   = min_rlcell_int

  art_atmo => p_art_data(jg)%atmo

  art_atmo%i_startblk  =  p_patch(jg)%cells%start_blk(i_rlstart,1)
  art_atmo%i_endblk    =  p_patch(jg)%cells%end_blk(i_rlend,i_nchdom)
  art_atmo%nblks       =  p_patch(jg)%nblks_c
  art_atmo%nlev        =  p_patch(jg)%nlev
  art_atmo%npromz      =  p_patch(jg)%npromz_c
  art_atmo%nproma      =  nproma
  art_atmo%nlevp1      =  p_patch(jg)%nlevp1
  art_atmo%nbisect     =  p_patch(jg)%level
  art_atmo%nroot       =  nroot

  ALLOCATE(art_atmo%lat(art_atmo%nproma,art_atmo%nblks))
  ALLOCATE(art_atmo%lon(art_atmo%nproma,art_atmo%nblks))

  !$ACC ENTER DATA CREATE(art_atmo%lat, art_atmo%lon)

  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
    DO jc = i_startidx, i_endidx
      art_atmo%lat(jc,jb) = p_patch(jg)%cells%center(jc,jb)%lat
      art_atmo%lon(jc,jb) = p_patch(jg)%cells%center(jc,jb)%lon
    END DO
  END DO

  art_atmo%cell_area   => p_patch(jg)%cells%area(:,:)
END SUBROUTINE art_collect_indices

SUBROUTINE art_collect_atmo_state(jg,p_nh_state)
!<
! SUBROUTINE art_collect_atmo_state
! This routine collects all physical variables that are not dependent on the
! physics package used in ICON. In addition, the fields that are dependent on
! the physics are allocated here. They are used as allocatable arrays because
! they have to be calculated in one of the physics packages as they are not
! available there.
! Part of Module: mo_art_collect_atmo_state
! Author: Jennifer Schroeter, KIT
! Initial Release: around 2018-10
! Modifications:
!>
  IMPLICIT NONE

  INTEGER, INTENT(in) ::   &
    &  jg                  !< patch id
  TYPE(t_nh_state), INTENT(in) ::  &
    &  p_nh_state          !< diagnostic variables

  ! local variable
  TYPE(t_art_atmo),POINTER    :: &
    &  art_atmo            !< Pointer to ART atmospheric fields


  art_atmo => p_art_data(jg)%atmo

  ALLOCATE(p_art_data(jg)%atmo%theta(art_atmo%nproma, art_atmo%nlev, art_atmo%nblks))
  ALLOCATE(p_art_data(jg)%atmo%llsm(art_atmo%nproma, art_atmo%nblks))
  ALLOCATE(p_art_data(jg)%atmo%gz0(art_atmo%nproma,art_atmo%nblks))
  ALLOCATE(p_art_data(jg)%atmo%swflx_par_sfc(art_atmo%nproma,art_atmo%nblks))
  ALLOCATE(p_art_data(jg)%atmo%rh_2m(art_atmo%nproma,art_atmo%nblks))
  ALLOCATE(p_art_data(jg)%atmo%o3_clim(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks))

  ALLOCATE(p_art_data(jg)%atmo%sza(art_atmo%nproma,art_atmo%nblks))
  ALLOCATE(p_art_data(jg)%atmo%sza_deg(art_atmo%nproma,art_atmo%nblks))
  ALLOCATE(p_art_data(jg)%atmo%ptropo(art_atmo%nproma, art_atmo%nblks))
  ALLOCATE(p_art_data(jg)%atmo%ktrpwmop1(art_atmo%nproma, art_atmo%nblks))
  ALLOCATE(p_art_data(jg)%atmo%ktrpwmo(art_atmo%nproma, art_atmo%nblks))
  ALLOCATE(p_art_data(jg)%atmo%ktrpwmop1_real(art_atmo%nproma, art_atmo%nblks))

  ALLOCATE(p_art_data(jg)%atmo%swflxsfc(art_atmo%nproma,art_atmo%nblks))
  ALLOCATE(p_art_data(jg)%atmo%t_2m(art_atmo%nproma,art_atmo%nblks))

  art_atmo%temp     => p_nh_state%diag%temp
  __acc_attach(art_atmo%temp)
  art_atmo%tempv    => p_nh_state%diag%tempv
  art_atmo%temp_ifc => p_nh_state%diag%temp_ifc
  art_atmo%pres     => p_nh_state%diag%pres
  __acc_attach(art_atmo%pres)

  art_atmo%z_mc     => p_nh_state%metrics%z_mc(:,:,:)
  art_atmo%z_ifc    => p_nh_state%metrics%z_ifc(:,:,:)
  art_atmo%dz       => p_nh_state%metrics%ddqz_z_full(:,:,:)
  __acc_attach(art_atmo%dz)
  art_atmo%o3_ext   => art_atmo%o3_clim

  art_atmo%pres_ifc => p_nh_state%diag%pres_ifc(:,:,:)
  art_atmo%pres_sfc => p_nh_state%diag%pres_sfc(:,:)
  art_atmo%dpres_mc => p_nh_state%diag%dpres_mc(:,:,:)
  art_atmo%u        => p_nh_state%diag%u
  art_atmo%v        => p_nh_state%diag%v
  art_atmo%vor      => p_nh_state%diag%vor

  !$ACC ENTER DATA CREATE(p_art_data(jg)%atmo%gz0, p_art_data(jg)%atmo%ktrpwmop1_real, p_art_data(jg)%atmo%ktrpwmop1) &
  !$ACC   CREATE(p_art_data(jg)%atmo%ktrpwmo, p_art_data(jg)%atmo%llsm, p_art_data(jg)%atmo%ptropo) &
  !$ACC   CREATE(p_art_data(jg)%atmo%rh_2m, p_art_data(jg)%atmo%swflxsfc, p_art_data(jg)%atmo%swflx_par_sfc) &
  !$ACC   CREATE(p_art_data(jg)%atmo%sza, p_art_data(jg)%atmo%sza_deg, p_art_data(jg)%atmo%theta) &
  !$ACC   CREATE(p_art_data(jg)%atmo%t_2m)

END SUBROUTINE art_collect_atmo_state


SUBROUTINE art_collect_atmo_state_nwp(jg,mtime_current,p_nh_state, ext_data, prm_diag, p_prog)
!<
! SUBROUTINE art_collect_atmo_state_nwp
! This routine collects all physical variables that are dependent on the
! physics package used in ICON (in this case NWP physics)
! Part of Module: mo_art_collect_atmo_state
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-05
! Modifications:
!>
  IMPLICIT NONE

  INTEGER, INTENT(in) ::   &
    &  jg                !< patch id
  TYPE(datetime), POINTER ::  &
    &  mtime_current     !< current date and time
  TYPE(t_nh_state), INTENT(in)  ::   &
    &  p_nh_state        !< diagnostic variables, independent of physics package
  TYPE(t_external_data), INTENT(in) ::  &
    &  ext_data          !< external data (boundary etc.) for NWP physics
  TYPE(t_nwp_phy_diag), INTENT(in)  ::   &
    &  prm_diag          !< diagnostic fields for NWP physics
  TYPE(t_nh_prog), INTENT(in) ::  &
    &  p_prog            !< prognostic variables (such as dynamics)

  ! Local variables
  TYPE (t_art_atmo), POINTER    :: &
    &  art_atmo           !< Pointer to ART atmospheric fields

  CALL art_collect_atmo_state(jg, p_nh_state)

  art_atmo => p_art_data(jg)%atmo

  art_atmo%llsm          =  ext_data%atm%llsm_atm_c
  art_atmo%fr_land       => ext_data%atm%fr_land
  art_atmo%fr_glac       => ext_data%atm%fr_glac
  art_atmo%clc           => prm_diag%clc
  art_atmo%tot_cld       => prm_diag%tot_cld
  art_atmo%acdnc         => prm_diag%acdnc
  art_atmo%o3_field_icon => ext_data%atm%o3
  art_atmo%albedo        => prm_diag%albdif
  art_atmo%u_10m         => prm_diag%u_10m
  art_atmo%v_10m         => prm_diag%v_10m
  art_atmo%rain_gsp_rate => prm_diag%rain_gsp_rate
  __acc_attach(art_atmo%rain_gsp_rate)
  art_atmo%rain_con_rate => prm_diag%rain_con_rate
  __acc_attach(art_atmo%rain_con_rate)
  art_atmo%tch           => prm_diag%tch
  art_atmo%tcm           => prm_diag%tcm
  art_atmo%lai           => ext_data%atm%lai

  CALL art_update_atmo_state_nwp(jg,mtime_current, p_prog, prm_diag)
END SUBROUTINE art_collect_atmo_state_nwp

SUBROUTINE art_collect_radiation_properties(iforcing, jg)
  INTEGER, INTENT(in) :: iforcing
  INTEGER, INTENT(in) :: jg

  TYPE (t_art_atmo), POINTER :: art_atmo

  art_atmo => p_art_data(jg)%atmo

  IF (iforcing == inwp) THEN
    SELECT CASE (atm_phy_nwp_config(jg)%inwp_radiation)
      CASE(1)
        art_atmo%relmin = rrtm_relmin
        art_atmo%relmax = rrtm_relmax
        art_atmo%reimin = rrtm_reimin
        art_atmo%reimax = rrtm_reimax
      CASE(3)
        art_atmo%relmin = rte_rrtmgp_relmin
        art_atmo%relmax = rte_rrtmgp_relmax
        art_atmo%reimin = rte_rrtmgp_reimin
        art_atmo%reimax = rte_rrtmgp_reimax
    END SELECT
  ENDIF

  IF (iforcing == iaes) THEN
#ifdef __NO_AES__
    CALL finish('mo_art_collect_atmo_state:art_collect_radiation_properties', &
              & ' For AES support remove --disable-aes and reconfigure')
#else
    art_atmo%relmin = rte_rrtmgp_relmin
    art_atmo%relmax = rte_rrtmgp_relmax
    art_atmo%reimin = rte_rrtmgp_reimin
    art_atmo%reimax = rte_rrtmgp_reimax
#endif
  ENDIF

END SUBROUTINE art_collect_radiation_properties

SUBROUTINE art_collect_atmo_state_aes(jg,mtime_current,p_nh_state, p_prog)
!<
! SUBROUTINE art_collect_atmo_state_aes
! This routine collects all physical variables that are dependent on the
! physics package used in ICON (in this case AES physics)
! Part of Module: mo_art_collect_atmo_state
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-05
! Modifications:
!>
  IMPLICIT NONE

  INTEGER, INTENT(in) :: &
    &  jg               !< patch id
  TYPE(datetime), POINTER ::   &
    &  mtime_current    !< current date and time
  TYPE(t_nh_state), INTENT(IN)  ::   &
    &  p_nh_state       !< diagnostic variables independent of physics
  TYPE(t_nh_prog), TARGET, INTENT(in) :: &
    &  p_prog           !< prognostic variables

  ! Local variables
  INTEGER :: &
    &  i_startidx, i_endidx, jb, jc
  TYPE(t_art_atmo),POINTER    :: &
     &  art_atmo           !< Pointer to ART atmospheric fields

#ifdef __NO_AES__
  CALL finish('mo_art_collect_atmo_state:art_collect_atmo_state_aes',    &
            & ' For AES support remove --disable-aes and reconfigure')
#else
  CALL art_collect_atmo_state(jg, p_nh_state)

  art_atmo => p_art_data(jg)%atmo
  art_atmo%fr_land  => prm_field(jg)%lsmask

  DO jb = art_atmo%i_startblk,art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

    DO jc = i_startidx, i_endidx
      art_atmo%llsm(jc,jb) = (art_atmo%fr_land(jc,jb) > 0.5_wp)
    END DO
  END DO

  art_atmo%fr_glac       => prm_field(jg)%glac
  art_atmo%clc           => prm_field(jg)%aclc
  art_atmo%tot_cld       => prm_field(jg)%qtrc_phy
  art_atmo%acdnc         => prm_field(jg)%acdnc
  art_atmo%o3_field_icon => prm_field(jg)%o3
  art_atmo%albedo        => prm_field(jg)%albedo
  art_atmo%u_10m         => prm_field(jg)%uas
  art_atmo%v_10m         => prm_field(jg)%vas
  art_atmo%rain_gsp_rate => prm_field(jg)%rsfl
  __acc_attach(art_atmo%rain_gsp_rate)
  art_atmo%tch           => prm_field(jg)%cfh(:,UBOUND(prm_field(jg)%cfh,2),:)
  art_atmo%tcm           => prm_field(jg)%cfm(:,UBOUND(prm_field(jg)%cfm,2),:)

! Leaf area index
#ifdef __NO_JSBACH__
  art_atmo%lai  => NULL()
#else
  IF (aes_phy_config(jg)%ljsb) THEN
    ! LAI from JSBACH, currently not available
    CALL jsbach_get_var('pheno_lai',1,ptr2d=art_atmo%lai)
  ELSE
    art_atmo%lai  => NULL()
  ENDIF
#endif

  IF ((.NOT. ASSOCIATED(art_atmo%lai))                      &
        & .AND. (ASSOCIATED(p_art_data(jg)%ext%land%pft))) THEN

    CALL finish('mo_art_collect_atmo_state:art_collect_atmo_state_aes',    &
            &   'Found emiss_onlBIO for at least one tracer'                 &
            & //' although LAI is not available.')
  END IF

  CALL art_update_atmo_state_aes(jg, mtime_current, p_prog)
#endif
END SUBROUTINE art_collect_atmo_state_aes
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

SUBROUTINE art_update_atmo_state(jg, mtime_current, p_prog, lacc)
!<
! SUBROUTINE art_update_atmo_state
! This routine updates all time-dependent physical variables that are not dependent on the
! physics package used in ICON. This is needed because the prognostic variables
! in ICON change according to the time level (nnow and nnew) that is used.
! Part of Module: mo_art_collect_atmo_state
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-05
! Modifications:
!>
  INTEGER, INTENT(in) :: &
    &  jg              !< patch id
  TYPE(datetime), POINTER :: &
    &  mtime_current   !< current date and time
  TYPE(t_nh_prog), TARGET, INTENT(in) :: &
    &  p_prog          !< prognostic variables
  LOGICAL, INTENT(in), OPTIONAL :: lacc ! If true, use openacc
  ! local
  REAL(wp), PARAMETER :: &
    &  p0 = 100000._wp,  &  !< pressure for potential temperature
    &  kappa = 0.286_wp     !< constant for potential temperature
  INTEGER ::       &
    &  jc,jk,jb,   &
    &  i_startidx, &
    &  i_endidx
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo        !< ART atmo fields
  LOGICAL :: lzacc             ! OpenACC flag
  CALL set_acc_host_or_device(lzacc, lacc)

  art_atmo => p_art_data(jg)%atmo

  ! The p_prog structure differs from time step to time step (nnew and nnow
  ! indices), so the pointers have to be updated at each time step
  art_atmo%rho     => p_prog%rho
  __acc_attach(art_atmo%rho)
  art_atmo%exner   => p_prog%exner
  __acc_attach(art_atmo%exner)
  art_atmo%tke     => p_prog%tke
  __acc_attach(art_atmo%tke)
  art_atmo%theta_v => p_prog%theta_v
  __acc_attach(art_atmo%theta_v)

  !$ACC DATA PRESENT(art_atmo) IF(lzacc)

  DO jb = art_atmo%i_startblk,art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1, art_atmo%nlev
      DO jc = i_startidx, i_endidx
        art_atmo%theta(jc,jk,jb) = art_atmo%temp(jc,jk,jb)    &
                           &       * (p0 / art_atmo%pres(jc,jk,jb)) ** kappa
      END DO
    END DO
    !$ACC END PARALLEL
  END DO

  !$ACC END DATA

  ! Update SZA
  CALL art_calc_sza(mtime_current, jg, lacc=lzacc)
END SUBROUTINE art_update_atmo_state

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

SUBROUTINE art_update_atmo_state_nwp(jg, mtime_current, p_prog, prm_diag, lacc)
!<
! SUBROUTINE art_update_atmo_state_nwp
! This routine updates all time-dependent physical variables that depend on the
! physics package used in ICON (in this case NWP). This is needed because the prognostic variables
! in ICON change according to the time level (nnow and nnew) that is used.
! Part of Module: mo_art_collect_atmo_state
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-05
! Modifications:
!>
  INTEGER, INTENT(in) :: &
    &  jg             !< patch id
  TYPE(datetime), POINTER :: &
    &  mtime_current  !< current date and time
  TYPE(t_nh_prog), TARGET, INTENT(in) :: &
    &  p_prog         !< prognostic variables
  TYPE(t_nwp_phy_diag), INTENT(in) :: &
    &  prm_diag       !< diagnostic fields in NWP physics
  LOGICAL, INTENT(in), OPTIONAL :: lacc ! If true, use openacc
  ! local
  INTEGER ::        &
    &  jb, jc,      &
    &  istart,      &
    &  iend
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo       !< ART atmo fields
  LOGICAL :: lzacc             ! OpenACC flag
  CALL set_acc_host_or_device(lzacc, lacc)

  art_atmo => p_art_data(jg)%atmo

  CALL art_update_atmo_state(jg, mtime_current, p_prog, lacc=lzacc)

  !$ACC DATA PRESENT(art_atmo, prm_diag) IF(lzacc)

  ! The acllocatable arrays that depend on time have to be updated each time
  ! step as well
  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, istart, iend)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO jc = istart, iend
      art_atmo%swflx_par_sfc(jc,jb) = prm_diag%swflx_par_sfc(jc,jb)
      art_atmo%swflxsfc(jc,jb)      = prm_diag%swflxsfc(jc,jb)
      art_atmo%t_2m(jc,jb)          = prm_diag%t_2m(jc,jb)
      art_atmo%rh_2m(jc,jb)         = prm_diag%rh_2m(jc,jb)
      art_atmo%gz0(jc,jb)           = prm_diag%gz0(jc,jb)
    END DO
    !$ACC END PARALLEL
  END DO

  !$ACC END DATA

  ! update ozone climatology
  IF (art_config(jg)%lart_chem) THEN
    IF (art_config(jg)%lart_chemtracer) THEN
      CALL calc_o3_gems_linoz(p_patch(jg),mtime_current,art_atmo%o3_clim)
    END IF
  END IF
  CALL art_calc_tropopause(jg,art_atmo%nlev-20,10,lacc=lzacc)
END SUBROUTINE art_update_atmo_state_nwp

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

SUBROUTINE art_update_atmo_state_aes(jg, mtime_current, p_prog)
!<
! SUBROUTINE art_update_atmo_state_aes
! This routine updates all time-dependent physical variables that depend on the
! physics package used in ICON (in this case AES physics). This is needed
! because the prognostic variables in ICON change according to the time level
! (nnow and nnew) that is used.
! Part of Module: mo_art_collect_atmo_state
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-05
! Modifications:
!>
  INTEGER, INTENT(in) :: &
    &  jg              !< patch id
  TYPE(datetime), POINTER :: &
    &  mtime_current   !< current date and time
  TYPE(t_nh_prog), TARGET, INTENT(in) :: &
    &  p_prog          !< prognostic variables
  ! local
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo        !< ART atmo fields
  INTEGER :: &
    &  jb, istart, iend, &
    &  jc
  REAL(wp), ALLOCATABLE :: &
    &  zo3_timint(:,:)

#ifdef __NO_AES__
  CALL finish('mo_art_collect_atmo_state:art_update_atmo_state_aes',    &
            & ' For AES support remove --disable-aes and reconfigure')
#else
  art_atmo => p_art_data(jg)%atmo

  CALL art_update_atmo_state(jg, mtime_current, p_prog)

  ! The allocatable arrays that depend on time have to be updated each time
  ! step as well

  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, istart, iend)

    DO jc = istart, iend
      art_atmo%swflx_par_sfc(jc,jb) = prm_field(jg)%rpds_dir(jc,jb) + prm_field(jg)%rpds_dif(jc,jb)
      art_atmo%swflxsfc(jc,jb)      = prm_field(jg)%rvds_dir(jc,jb) + prm_field(jg)%rvds_dif(jc,jb) +  &
        &                             prm_field(jg)%rnds_dir(jc,jb) + prm_field(jg)%rnds_dif(jc,jb)
      art_atmo%t_2m(jc,jb)          = prm_field(jg)%tas(jc,jb)
      art_atmo%rh_2m(jc,jb)         = td_and_t2rh(prm_field(jg)%dew2(jc,jb), prm_field(jg)%tas(jc,jb))
      art_atmo%gz0(jc,jb)           = grav * prm_field(jg)%z0m(jc,jb)
    END DO
  END DO

  IF (art_config(jg)%lart_chem) THEN
    IF (art_config(jg)%lart_chemtracer) THEN
      IF (ALLOCATED(ext_ozone)) THEN
        !! O3 climatology
        ALLOCATE(zo3_timint(art_atmo%nproma,ext_ozone(jg)%nplev_o3))

        DO jb = art_atmo%i_startblk, art_atmo%i_endblk
          CALL art_get_indices_c(jg, jb, istart, iend)

          CALL o3_timeint(jcs = istart, jce = iend, kbdim = art_atmo%nproma,  &
               &          nlev_pres=ext_ozone(jg)%nplev_o3,              &
               &          ext_o3=ext_ozone(jg)%o3_plev(:,:,jb,:),        &
               &          current_date=mtime_current,                    &
               &          o3_time_int=zo3_timint                 )

          CALL o3_pl2ml ( jcs = istart,jce = iend, kbdim = art_atmo%nproma, &
                 &        nlev_pres = ext_ozone(jg)%nplev_o3,    &
                 &        klev = art_atmo%nlev,                  &
                 &        pfoz = ext_ozone(jg)%plev_full_o3,     &
                 &        phoz = ext_ozone(jg)%plev_half_o3,     &
                 &        ppf  = art_atmo%pres(:,:,jb),  &
                 &        pph  = art_atmo%pres_ifc(:,:,jb),  &
                 &        o3_time_int = zo3_timint,              &
                 &        o3_clim     = art_atmo%o3_clim(:,:,jb))

        ENDDO
      END IF
    END IF
  END IF

  CALL art_calc_tropopause(jg)
#endif
END SUBROUTINE art_update_atmo_state_aes

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

SUBROUTINE art_deallocate_atmo_state(jg)
!<
! SUBROUTINE art_deallocate_atmo_state
! Deallocation of the structure
! Part of Module: mo_art_collect_atmo_state
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-05
! Modifications:
!>
  IMPLICIT NONE
  INTEGER, INTENT(in) :: &
    &  jg          !< patch id
  ! local variables
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo    !< ART atmo fields

  art_atmo => p_art_data(jg)%atmo


  IF (ALLOCATED(art_atmo%o3_clim))        DEALLOCATE(art_atmo%o3_clim)

  IF (ALLOCATED(art_atmo%gz0)) THEN
    !$ACC EXIT DATA DELETE(art_atmo%gz0)
    DEALLOCATE(art_atmo%gz0)
  END IF
  IF (ALLOCATED(art_atmo%swflx_par_sfc)) THEN
    !$ACC EXIT DATA DELETE(art_atmo%swflx_par_sfc)
    DEALLOCATE(art_atmo%swflx_par_sfc)
  END IF
  IF (ALLOCATED(art_atmo%t_2m)) THEN
    !$ACC EXIT DATA DELETE(art_atmo%t_2m)
    DEALLOCATE(art_atmo%t_2m)
  END IF
  IF (ALLOCATED(art_atmo%lat)) THEN
    !$ACC EXIT DATA DELETE(art_atmo%lat)
    DEALLOCATE(art_atmo%lat)
  END IF
  IF (ALLOCATED(art_atmo%lon)) THEN
    !$ACC EXIT DATA DELETE(art_atmo%lon)
    DEALLOCATE(art_atmo%lon)
  END IF
  IF (ALLOCATED(art_atmo%llsm)) THEN
    !$ACC EXIT DATA DELETE(art_atmo%llsm)
    DEALLOCATE(art_atmo%llsm)
  END IF
  IF (ALLOCATED(art_atmo%rh_2m)) THEN
    !$ACC EXIT DATA DELETE(art_atmo%rh_2m)
    DEALLOCATE(art_atmo%rh_2m)
  END IF
  IF (ALLOCATED(art_atmo%swflxsfc)) THEN
    !$ACC EXIT DATA DELETE(art_atmo%swflxsfc)
    DEALLOCATE(art_atmo%swflxsfc)
  END IF
  IF (ALLOCATED(art_atmo%theta)) THEN
    !$ACC EXIT DATA DELETE(art_atmo%theta)
    DEALLOCATE(art_atmo%theta)
  END IF

  !$ACC EXIT DATA DELETE(art_atmo%sza)
  DEALLOCATE(art_atmo%sza)
  NULLIFY(art_atmo%sza)
  !$ACC EXIT DATA DELETE(art_atmo%sza_deg)
  DEALLOCATE(art_atmo%sza_deg)
  NULLIFY(art_atmo%sza_deg)
  !$ACC EXIT DATA DELETE(art_atmo%ptropo)
  DEALLOCATE(art_atmo%ptropo)
  NULLIFY(art_atmo%ptropo)
  !$ACC EXIT DATA DELETE(art_atmo%ktrpwmop1_real)
  DEALLOCATE(art_atmo%ktrpwmop1_real)
  NULLIFY(art_atmo%ktrpwmop1_real)
  !$ACC EXIT DATA DELETE(art_atmo%ktrpwmop1)
  DEALLOCATE(art_atmo%ktrpwmop1)
  NULLIFY(art_atmo%ktrpwmop1)
  !$ACC EXIT DATA DELETE(art_atmo%ktrpwmo)
  DEALLOCATE(art_atmo%ktrpwmo)
  NULLIFY(art_atmo%ktrpwmo)

  NULLIFY(art_atmo%temp)
  NULLIFY(art_atmo%tempv)
  NULLIFY(art_atmo%temp_ifc)
  NULLIFY(art_atmo%pres)
  NULLIFY(art_atmo%exner)
  NULLIFY(art_atmo%vor)
  NULLIFY(art_atmo%rho)
  NULLIFY(art_atmo%theta_v)
  NULLIFY(art_atmo%z_mc)
  NULLIFY(art_atmo%z_ifc)
  NULLIFY(art_atmo%dz)
  NULLIFY(art_atmo%pres_ifc)
  NULLIFY(art_atmo%dpres_mc)
  NULLIFY(art_atmo%o3_field_icon)
  NULLIFY(art_atmo%o3_ext)
  NULLIFY(art_atmo%fr_land)
  NULLIFY(art_atmo%cell_area)
  NULLIFY(art_atmo%fr_glac)
  NULLIFY(art_atmo%clc)
  NULLIFY(art_atmo%tot_cld)
  NULLIFY(art_atmo%acdnc)
  NULLIFY(art_atmo%albedo)
  NULLIFY(art_atmo%lai)
  NULLIFY(art_atmo%u)
  NULLIFY(art_atmo%v)
  NULLIFY(art_atmo%tke)
  NULLIFY(art_atmo%pres_sfc)
  NULLIFY(art_atmo%u_10m)
  NULLIFY(art_atmo%v_10m)
  NULLIFY(art_atmo%rain_gsp_rate)
  NULLIFY(art_atmo%rain_con_rate)
  NULLIFY(art_atmo%tch)
  NULLIFY(art_atmo%tcm)
END SUBROUTINE art_deallocate_atmo_state

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

ELEMENTAL FUNCTION td_and_t2rh(td, t) RESULT (rh)
  REAL(wp) :: rh                ! result in [%]
  REAL(wp), INTENT(in) :: td    ! input expected in [K]
  REAL(wp), INTENT(in) :: t     ! input expected in [K]
  REAL(wp), PARAMETER :: b =  17.62_wp
  REAL(wp), PARAMETER :: c = -30.03_wp
  REAL(wp), PARAMETER :: t0 = 273.15_wp

  rh = 100.0_wp*(EXP((b*(td-t0))/(c+td))/EXP((b*(t-t0))/(c+t)))

END FUNCTION td_and_t2rh

!!
!!-------------------------------------------------------------------------
!!

SUBROUTINE art_init_tracer_values(jg, tracer, p_prog_list)
!<
! SUBROUTINE art_init_tracer_values
! Initialisation of tracer values. At this stage, the atmospheric state is
! already available so that tracers can be initialised temperature and pressure
! dependent.
! Part of Module: mo_art_collect_atmo_state
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-11
! Modifications:
!>
  IMPLICIT NONE

  INTEGER,INTENT(in)          :: &
    &  jg                          !< number of model domain
  REAL(wp), POINTER, INTENT(inout)  :: &
    &  tracer(:,:,:,:)             !< Tracer mixing ratios [kg kg-1]
  TYPE(t_var_list_ptr), INTENT(IN)      ::  &
    &  p_prog_list                 !< current prognostic state list

  ! local variables
  TYPE(t_art_atmo),POINTER    :: &
    &  art_atmo                    !< Pointer to ART atmo fields
  TYPE(t_mode),POINTER        :: &
    &  current_mode                !< pointer to loop through mode structure
  TYPE(t_var), POINTER        :: &
    &  var_list_elem               !< pointer to element of variable list
  INTEGER                     :: &
    &  grp_id, jsp, tr_idx, ierror

  art_atmo => p_art_data(jg)%atmo

  IF (art_config(jg)%lart_chem) THEN

    IF (art_config(jg)%lart_chemtracer) THEN
      CALL art_init_chemtracer(jg, tracer)
    END IF

    IF (art_config(jg)%lart_mecca) THEN
#ifdef __ART_GPL
      CALL art_init_full_chemistry(jg, p_prog_list, tracer)
#else
      CALL finish('mo_art_collect_atmo_state', &
          &     'You have activated MECCA which is software under GPL licence. Deactivate MECCA or activate GPL software for ART.')
#endif
    END IF

    SELECT CASE(art_config(jg)%iart_init_gas)
      CASE(0)
        ! Nothing to do here, tracers are initialized with 0.0_wp automatically
      CASE(1)
        ! At this point climatological gas profiles may be included as a standard initialization
      CASE(4)
        CALL message('', 'ART: Initializing chemical tracers using method given by .xml')

        CALL art_init_chemtracer_ext(p_patch(jg),p_prog_list, tracer)

      CASE(5)
        CALL message('','ART: Initializing chemical tracers')

        CALL art_init_chem(jg, p_prog_list, tracer)

      CASE DEFAULT
        CALL finish('mo_art_init:art_init', &
          &      'Unknown initialization action for gases.')
    END SELECT

    IF (art_config(jg)%iart_init_gas /= 5) THEN
      CALL art_convert_tracers_vmr_mmr(tracer, p_prog_list)
    END IF

  END IF
  ! ----------------------------------
  ! --- 2.0 Set tracer initial values from files
  ! ----------------------------------
  IF (art_config(jg)%lart_aerosol) THEN
    SELECT CASE(art_config(jg)%iart_init_aero)
      CASE(0)
        ! CASE(0) is the new standard case; a fix number concentration of 100 #/kg per mode is
        ! prescribed. This number concentration together with the corresponding mass mixing ratios
        ! is set in the mode meta data in set_meta_init_nmb_mass_conc, here only the tracer fields
        ! are initialized with these values

        ! id of the group of tracers initialized together with the meteorological first guess fields
        grp_id = var_groups_dyn%group_id('tracer_fg_in')

        ! Loop through modes
        current_mode => p_mode_state(jg)%p_mode_list%p%first_mode
        DO WHILE(ASSOCIATED(current_mode))

        ! Select type of mode
          SELECT TYPE (fields => current_mode%fields)
            CLASS IS (t_fields_2mom)

              ! check if the number concentration tracer was not already initialized with the first guess,
              ! i.e. is not in the respective group
              var_list_elem => find_tracer_by_index(p_prog_list, fields%itr0)
              IF (.NOT. var_list_elem%info%in_group(grp_id)) THEN
                WRITE (message_text,'(3A,E12.5,A)') 'iart_init_aero is 0: mode ',TRIM(fields%name), &
                  &                     ' is initialized with a fix number concentration of ',      &
                  &                     fields%info%init_nmb_conc, ' #/kg'
                CALL message (TRIM(routine)//':art_init', message_text)

                ! write number concentration to tracer field
                tracer(:,:,:,fields%itr0) = fields%info%init_nmb_conc
              END IF

              DO jsp=1,fields%ntr-1
                ! check if the mass mixing ratio tracer was not already initialized with the first guess,
                ! i.e. is not in the respective group
                var_list_elem => find_tracer_by_index(p_prog_list, fields%itr3(jsp))
                IF (.NOT. var_list_elem%info%in_group(grp_id)) THEN
                  ! write mass mixing ratio to tracer field
                  tracer(:,:,:,fields%itr3(jsp)) = fields%info%init_mass_conc(jsp)
                END IF
              ENDDO !jsp

            CLASS DEFAULT
              ! modes other than t_fields_2mom are initzialized with 0
!              CALL finish('mo_art_init','initializing tracer in CASE(0) failed')
          END SELECT !fields
          current_mode => current_mode%next_mode
        ENDDO !ASSOCIATED(current_mode)

      CASE(1)
        ! At this point climatological aerosol profiles may be included as standard initialization

      CASE(2)
        ! Initialize with fixed value different from 0.0_wp (using 10.0_wp)
        CALL p_art_data(jg)%dict_tracer%get('so4_sol_ait',tr_idx,ierror)
        IF( ierror /= 0 ) THEN
          CALL finish('mo_art_init','initializing so4_sol_ait failed')
        END IF
        tracer(:,:,:,tr_idx) = 0.1_wp

        ! Initialize with fixed value different from 0.0_wp (using 10.0_wp)
        CALL p_art_data(jg)%dict_tracer%get('nh4_sol_ait',tr_idx,ierror)
        IF( ierror /= 0 ) THEN
          CALL finish('mo_art_init','initializing nh4_sol_ait failed')
        END IF
         tracer(:,:,:,tr_idx) = 0.04_wp !

        ! Initialize with fixed value different from 0.0_wp (using 10.0_wp)
        CALL p_art_data(jg)%dict_tracer%get('no3_sol_ait',tr_idx,ierror)
        IF( ierror /= 0 ) THEN
          CALL finish('mo_art_init','initializing no3_sol_ait failed')
        END IF
         tracer(:,:,:,tr_idx) = 0.01_wp !

        ! Initialize with fixed value different from 0.0_wp (using 10.0_wp)
        CALL p_art_data(jg)%dict_tracer%get('nmb_sol_ait',tr_idx,ierror)
        IF( ierror /= 0 ) THEN
          CALL finish('mo_art_init','initializing nmb_sol_ait failed')
        END IF

        ! Initialize with fixed value different from 0.0_wp (using 10.0_wp)
        CALL p_art_data(jg)%dict_tracer%get('so4_sol_acc',tr_idx,ierror)
        IF( ierror /= 0 ) THEN
          CALL finish('mo_art_init','initializing so4_sol_acc failed')
        END IF
        tracer(:,:,:,tr_idx) = 0.9_wp
!       tracer(:,:,:,tr_idx) = 0.09_wp

        ! Initialize with fixed value different from 0.0_wp (using 10.0_wp)
        CALL p_art_data(jg)%dict_tracer%get('nh4_sol_acc',tr_idx,ierror)
        IF( ierror /= 0 ) THEN
          CALL finish('mo_art_init','initializing nh4_sol_acc failed')
        END IF
         tracer(:,:,:,tr_idx) = 0.33_wp !

        ! Initialize with fixed value different from 0.0_wp (using 10.0_wp)
        CALL p_art_data(jg)%dict_tracer%get('no3_sol_acc',tr_idx,ierror)
        IF( ierror /= 0 ) THEN
          CALL finish('mo_art_init','initializing no3_sol_acc failed')
        END IF
         tracer(:,:,:,tr_idx) = 0.09_wp !

        ! Initialize with fixed value different from 0.0_wp (using 10.0_wp)
        CALL p_art_data(jg)%dict_tracer%get('nmb_sol_acc',tr_idx,ierror)
        IF( ierror /= 0 ) THEN
          CALL finish('mo_art_init','initializing nmb_sol_ait failed')
        END IF
 !       tracer(:,:,:,tr_idx) = 3.06e-09_wp / 3.14_wp / ((0.07e-06_wp**3)*EXP(
 !       4.5_wp * (LOG( 1.7_wp )**2)))
 !       print *,'ANZAHL ACC', tracer(:,:,:,tr_idx)

      CASE(5)
        ! Set aerosol tracer initial values from file
        CALL message('','ART: Initializing Aerosol')

        CALL art_init_aero(jg,p_prog_list,tracer)

      CASE(6)
        ! read background climatology from echam-ham
        CALL message('','ART: Initializing aerosol with echam-ham climatology')

        CALL art_init_aero(jg, p_prog_list, tracer, use_echam_climatology=.true.)
        
      CASE DEFAULT
        CALL finish('mo_art_init:art_init', &
          &      'Unknown initialization action for aerosol.')
    END SELECT
  ENDIF !lart_aerosol

  NULLIFY(art_atmo)
  
END SUBROUTINE art_init_tracer_values

!!
!!-------------------------------------------------------------------------
!!

SUBROUTINE art_init_tracer_values_nwp(jg, tracer, current_date, p_prog_list)
!<
! SUBROUTINE art_init_tracer_values_nwp
! Initialisation of tracer values. At this stage, the atmospheric state is
! already available so that tracers can be initialised temperature and pressure
! dependent.
! Part of Module: mo_art_collect_atmo_state
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-05
! Modifications:
!>
  IMPLICIT NONE

  INTEGER,INTENT(in)          :: &
    &  jg                          !< number of model domain
  REAL(wp), POINTER, INTENT(inout)  :: &
    &  tracer(:,:,:,:)             !< Tracer mixing ratios [kg kg-1]
  TYPE(datetime),POINTER,INTENT(IN) :: &
    &  current_date                !< current date and time
  TYPE(t_var_list_ptr), INTENT(IN)      ::  &
    &  p_prog_list                 !< current prognostic state list

  ! local variables
  TYPE(t_art_atmo),POINTER    :: &
    &  art_atmo                     !< Pointer to ART atmo fields
  TYPE(t_art_chem_indices),POINTER    :: &
    &  art_indices                  !< Pointer to ART tracer indices

  art_atmo => p_art_data(jg)%atmo
  art_indices  => p_art_data(jg)%chem%indices

  IF (art_config(jg)%lart_chem) THEN
    IF (art_config(jg)%lart_chemtracer) THEN
      art_atmo%o3_clim(:,:,:) = 0._wp
      CALL calc_o3_gems_linoz(p_patch(jg),current_date,art_atmo%o3_clim)
    END IF

    IF (art_indices%iTRO3  /= 0) THEN
      tracer(:,:,:,art_indices%iTRO3) = art_atmo%o3_clim(:,:,:) * amd / amo3
    ENDIF

    ! ----------------------------------
    ! CH4 initialisation
    ! ----------------------------------
    ! this value depends on the start date of ICON run:
    ! mean value from MACC reanalysis of 01/01/2004
    ! (ATTENTION: only valid for this initialisation time!)
    ! Values are valid for 90 level version of ICON

    IF (art_indices%iTRCH4   /=0) THEN
      CALL art_set_init_tracer_trop_strat(jg,                             &
                      &                   tracer(:,:,:,art_indices%iTRCH4),  &
                      &                   1.67e-06_wp,                    &
                      &                   9.9e-07_wp,                     &
                      &                   5)
    END IF

  END IF

  CALL art_init_tracer_values(jg,tracer,p_prog_list)


END SUBROUTINE art_init_tracer_values_nwp

!!
!!-------------------------------------------------------------------------
!!

SUBROUTINE art_init_tracer_values_aes(jg, tracer, current_date, p_prog_list)
!<
! SUBROUTINE art_init_tracer_values_aes
! Initialisation of tracer values. At this stage, the atmospheric state is
! already available so that tracers can be initialised temperature and pressure
! dependent.
! Part of Module: mo_art_collect_atmo_state
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-05
! Modifications:
!>
  IMPLICIT NONE

  INTEGER,INTENT(in)          :: &
    &  jg                          !< number of model domain
  REAL(wp), POINTER, INTENT(inout)  :: &
    &  tracer(:,:,:,:)             !< Tracer mixing ratios [kg kg-1]
  TYPE(datetime),POINTER,INTENT(in) ::   &
    &  current_date                !< current date and time
  TYPE(t_var_list_ptr), INTENT(IN)    :: &
    &  p_prog_list                 !< current prognostic state list


  ! local variables
  TYPE(t_art_atmo),POINTER    :: &
    &  art_atmo                     !< Pointer to ART atmo fields
  TYPE(t_art_chem_indices), POINTER :: &
    &  art_indices                     !< Pointer to ART indicesetrised chem fields

  INTEGER ::       &
    &  jb, jk, jc, &  !< loop indices
    &  istart,iend
  REAL(wp), ALLOCATABLE :: &
    &  zo3_timint(:,:)   !< time interpolated Ozone from climatology
  REAL(wp), PARAMETER :: &
    &  xp(3) = (/1.25e-01_wp,  683.0_wp, -1.43_wp/)  !< factors for CH4
  REAL(wp)::          &
    &  zx_m, zx_d,    &  !< additional factors for CH4
    &  gas_scenario      !< volume mixing ratio of CH4

#ifdef __NO_AES__
  CALL finish('mo_art_collect_atmo_state:art_init_tracer_values_aes',    &
            & ' For AES support remove --disable-aes and reconfigure')
#else
  art_atmo   => p_art_data(jg)%atmo
  art_indices  => p_art_data(jg)%chem%indices

  IF (art_config(jg)%lart_chem) THEN
    IF (art_config(jg)%lart_chemtracer) THEN
      !! O3 initialisation
      CALL read_bc_ozone(current_date%date%year,p_patch(jg),aes_rad_config(jg)%irad_o3)

      ALLOCATE(zo3_timint(art_atmo%nblks,ext_ozone(jg)%nplev_o3))

      art_atmo%o3_clim(:,:,:) = 0.0_wp

      DO jb = art_atmo%i_startblk, art_atmo%i_endblk
        CALL art_get_indices_c(jg, jb, istart, iend)

        CALL o3_timeint(jcs = 1, jce = iend, kbdim = art_atmo%nproma,  &
             &          nlev_pres=ext_ozone(jg)%nplev_o3,              &
             &          ext_o3=ext_ozone(jg)%o3_plev(:,:,jb,:),        &
             &          current_date=current_date,                     &
             &          o3_time_int=zo3_timint                 )


        CALL o3_pl2ml ( 1,jce = iend, kbdim = art_atmo%nproma, &
               &        nlev_pres = ext_ozone(jg)%nplev_o3,    &
               &        klev = art_atmo%nlev,                  &
               &        pfoz = ext_ozone(jg)%plev_full_o3,     &
               &        phoz = ext_ozone(jg)%plev_half_o3,     &
               &        ppf  = art_atmo%pres(:,:,jb),   &
               &        pph  =  art_atmo%pres_ifc(:,:,jb), &
               &        o3_time_int = zo3_timint,              &
               &        o3_clim     = art_atmo%o3_clim(:,:,jb))

      ENDDO


      IF (art_indices%iTRO3 /= 0) THEN
        DO jb = art_atmo%i_startblk, art_atmo%i_endblk
          CALL art_get_indices_c(jg, jb, istart, iend)

          DO jk = 1,art_atmo%nlev
            DO jc = istart,iend
              tracer(jc,jk,jb,art_indices%iTRO3) = art_atmo%o3_clim(jc,jk,jb) / (amo3/amd)
            END DO
          END DO
        END DO
      END IF


      ! passive ozone initialisation
      IF (art_indices%iTRO3_pas /= 0) THEN
        DO jb = art_atmo%i_startblk, art_atmo%i_endblk
          CALL art_get_indices_c(jg, jb, istart, iend)

          DO jk = 1,art_atmo%nlev
            DO jc = istart,iend
              tracer(jc,jk,jb,art_indices%iTRO3_pas) = art_atmo%o3_clim(jc,jk,jb) / (amo3/amd)
            END DO
          END DO
        END DO
      END IF

      ! CH4 initialisation in case of AES physics
      IF (art_indices%iTRCH4 /= 0) THEN
        DO jb = art_atmo%i_startblk, art_atmo%i_endblk
          CALL art_get_indices_c(jg, jb, istart, iend)

          DO jk = 1,art_atmo%nlev
            DO jc = istart,iend

              gas_scenario = ghg_ch4mmr/(amch4/amd)

              zx_m = (gas_scenario+xp(1)*gas_scenario)*0.5_wp
              zx_d = (gas_scenario-xp(1)*gas_scenario)*0.5_wp
              tracer(jc,jk,jb, art_indices%iTRCH4) = (1-(zx_d/zx_m)    &
                &                   * TANH(LOG(art_atmo%pres(jc,jk,jb) &
                &                    /xp(2)) /xp(3))) * zx_m

            ENDDO
          ENDDO
        ENDDO
      END IF
    END IF
  END IF

  CALL art_init_tracer_values(jg,tracer,p_prog_list)
#endif 

END SUBROUTINE art_init_tracer_values_aes


END MODULE mo_art_collect_atmo_state
