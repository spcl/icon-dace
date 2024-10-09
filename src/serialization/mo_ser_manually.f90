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

! Serialize any scalar and array component in selected global derived type
! variables and that were not initialized with add_var or add_ref.
!
! HOW TO:
! 1) Add a USE-ONLY statement for each derived type variable
! 2) Add ser_component() calls to ser_manually() for each component

MODULE mo_ser_manually

#ifdef SERIALIZE
  USE mo_ser_common,  ONLY: ser_component, t_ser_options
  USE mo_kind,        ONLY: sp, dp

#ifdef _OPENACC
  USE openacc,               ONLY: acc_is_present
#endif

  ! Global variables with scalars that are going to be serialized
  USE mo_nonhydro_state,     ONLY: p_nh_state
  USE mo_nwp_phy_state,      ONLY: phy_params
  USE mo_ext_data_state,     ONLY: ext_data

  USE mo_nh_pzlev_config,    ONLY: nh_pzlev_config
    IMPLICIT NONE

  PUBLIC :: ser_manually

  PRIVATE

  ! Test dataset
  TYPE :: t_test
    INTEGER :: int
    REAL(dp) :: double
    REAL(sp) :: float
    LOGICAL :: log_a, log_b
  END TYPE t_test
  TYPE(t_test), TARGET :: p_test
  LOGICAL :: lscalars_test_was_run = .FALSE.


  CONTAINS

  SUBROUTINE ser_manually(abs_threshold, rel_threshold, ser_mode, domain)
    INTEGER, INTENT(IN) :: abs_threshold
    INTEGER, INTENT(IN) :: rel_threshold
    INTEGER, INTENT(IN) :: ser_mode
    INTEGER, INTENT(IN) :: domain
    TYPE(t_ser_options) :: o

    o%abs_threshold = abs_threshold
    o%rel_threshold = rel_threshold
    o%lopenacc = .TRUE. ! all variables tested here should be available on GPU
    o%ser_mode = ser_mode
    o%domain = domain

    ! Make sure that serialbox and openACC work with scalars
    CALL ser_test_scalars()

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! List all scalars that should be serialized here !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL ser_component(o, "p_nh_state(jg)%metrics%zd_listdim", p_nh_state(domain)%metrics%zd_listdim)
    CALL ser_component(o, "p_nh_state(jg)%metrics%pg_listdim", p_nh_state(domain)%metrics%pg_listdim)
    CALL ser_component(o, "p_nh_state(jg)%metrics%nudge_c_dim", p_nh_state(domain)%metrics%nudge_c_dim)
    CALL ser_component(o, "p_nh_state(jg)%metrics%nudge_e_dim", p_nh_state(domain)%metrics%nudge_e_dim)
    CALL ser_component(o, "p_nh_state(jg)%metrics%bdy_halo_c_dim", p_nh_state(domain)%metrics%bdy_halo_c_dim)
    CALL ser_component(o, "p_nh_state(jg)%metrics%bdy_mflx_e_dim", p_nh_state(domain)%metrics%bdy_mflx_e_dim)

    CALL ser_component(o, "phy_params(jg)%kcon1", phy_params(domain)%kcon1)
    CALL ser_component(o, "phy_params(jg)%kcon2", phy_params(domain)%kcon2)
    CALL ser_component(o, "phy_params(jg)%tau", phy_params(domain)%tau)
    CALL ser_component(o, "phy_params(jg)%mfcfl", phy_params(domain)%mfcfl)
    CALL ser_component(o, "phy_params(jg)%tau0", phy_params(domain)%tau0)
    CALL ser_component(o, "phy_params(jg)%rhebc_land", phy_params(domain)%rhebc_land)
    CALL ser_component(o, "phy_params(jg)%rhebc_ocean", phy_params(domain)%rhebc_ocean)
    CALL ser_component(o, "phy_params(jg)%rhebc_land_trop", phy_params(domain)%rhebc_land_trop)
    CALL ser_component(o, "phy_params(jg)%rhebc_ocean_trop", phy_params(domain)%rhebc_ocean_trop)
    CALL ser_component(o, "phy_params(jg)%texc", phy_params(domain)%texc)
    CALL ser_component(o, "phy_params(jg)%qexc", phy_params(domain)%qexc)
    CALL ser_component(o, "phy_params(jg)%rcucov", phy_params(domain)%rcucov)
    CALL ser_component(o, "phy_params(jg)%rcucov_trop", phy_params(domain)%rcucov_trop)
    CALL ser_component(o, "phy_params(jg)%entrorg", phy_params(domain)%entrorg)
    CALL ser_component(o, "phy_params(jg)%detrpen", phy_params(domain)%detrpen)
    CALL ser_component(o, "phy_params(jg)%entrdd", phy_params(domain)%entrdd)
    CALL ser_component(o, "phy_params(jg)%entstpc1", phy_params(domain)%entstpc1)
    CALL ser_component(o, "phy_params(jg)%entstpc2", phy_params(domain)%entstpc2)
    CALL ser_component(o, "phy_params(jg)%rprcon", phy_params(domain)%rprcon)
    CALL ser_component(o, "phy_params(jg)%rdepths", phy_params(domain)%rdepths)
    CALL ser_component(o, "phy_params(jg)%lmfscv", phy_params(domain)%lmfscv)
    CALL ser_component(o, "phy_params(jg)%lmfmid", phy_params(domain)%lmfmid)
    CALL ser_component(o, "phy_params(jg)%lmfpen", phy_params(domain)%lmfpen)
    CALL ser_component(o, "phy_params(jg)%lmfdsnow", phy_params(domain)%lmfdsnow)
    CALL ser_component(o, "phy_params(jg)%lgrayzone_deepconv", phy_params(domain)%lgrayzone_deepconv)
    CALL ser_component(o, "phy_params(jg)%klaunch", phy_params(domain)%klaunch)
    CALL ser_component(o, "phy_params(jg)%ngwdlim", phy_params(domain)%ngwdlim)
    CALL ser_component(o, "phy_params(jg)%ngwdtop", phy_params(domain)%ngwdtop)
    CALL ser_component(o, "phy_params(jg)%nktopg", phy_params(domain)%nktopg)
    CALL ser_component(o, "phy_params(jg)%gkwake", phy_params(domain)%gkwake)
    CALL ser_component(o, "phy_params(jg)%gkdrag", phy_params(domain)%gkdrag)
    CALL ser_component(o, "phy_params(jg)%gfrcrit", phy_params(domain)%gfrcrit)
    CALL ser_component(o, "phy_params(jg)%grcrit", phy_params(domain)%grcrit)
    CALL ser_component(o, "phy_params(jg)%mean_charlen", phy_params(domain)%mean_charlen)
    CALL ser_component(o, "phy_params(jg)%k060", phy_params(domain)%k060)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! List all arrays that should be serialized here  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! outdated example as prep_adv now uses the variable list mechanism.
!    IF(ALLOCATED(prep_adv)) THEN
!#ifdef _OPENACC
!      ! prep_adv can be tested only if variable is already present on device
!      IF( acc_is_present(prep_adv(domain:domain)) .AND. &
!        & acc_is_present(prep_adv(domain)%vn_traj) ) THEN
!#else
!      IF(.TRUE.) THEN
!#endif
!        CALL ser_component(o, "mass_flx_me", prep_adv(domain)%mass_flx_me)
!        CALL ser_component(o, "mass_flx_ic", prep_adv(domain)%mass_flx_ic)
!        CALL ser_component(o, "vn_traj",     prep_adv(domain)%vn_traj)
!        CALL ser_component(o, "q_int",       prep_adv(domain)%q_int)
!        CALL ser_component(o, "q_ubc",       prep_adv(domain)%q_ubc)
!      ENDIF
!    ENDIF

    CALL ser_component(o, "atm%list_seaice%ncount", ext_data(domain)%atm%list_seaice%ncount)
    CALL ser_component(o, "atm%list_seaice%idx", ext_data(domain)%atm%list_seaice%idx)
    CALL ser_component(o, "atm%list_lake%ncount", ext_data(domain)%atm%list_lake%ncount)
    CALL ser_component(o, "atm%list_lake%idx", ext_data(domain)%atm%list_lake%idx)
    CALL ser_component(o, "atm%list_land%ncount", ext_data(domain)%atm%list_land%ncount)
    CALL ser_component(o, "atm%list_land%idx", ext_data(domain)%atm%list_land%idx)
    CALL ser_component(o, "atm%list_seawtr%ncount", ext_data(domain)%atm%list_seawtr%ncount)
    CALL ser_component(o, "atm%list_seawtr%idx", ext_data(domain)%atm%list_seawtr%idx)
    CALL ser_component(o, "atm%emis_rad", ext_data(domain)%atm%emis_rad)
    CALL ser_component(o, "atm%z0_lcc", ext_data(domain)%atm%z0_lcc)
    CALL ser_component(o, "atm%z0_lcc_min", ext_data(domain)%atm%z0_lcc_min)
    CALL ser_component(o, "atm%plcovmax_lcc", ext_data(domain)%atm%plcovmax_lcc)
    CALL ser_component(o, "atm%laimax_lcc", ext_data(domain)%atm%laimax_lcc)
    CALL ser_component(o, "atm%rootdmax_lcc", ext_data(domain)%atm%rootdmax_lcc)
    CALL ser_component(o, "atm%stomresmin_lcc", ext_data(domain)%atm%stomresmin_lcc)
    CALL ser_component(o, "atm%snowalb_lcc", ext_data(domain)%atm%snowalb_lcc)
    CALL ser_component(o, "atm%snowtile_lcc", ext_data(domain)%atm%snowtile_lcc)
    CALL ser_component(o, "atm%t_cl", ext_data(domain)%atm%t_cl)

    IF(ALLOCATED(nh_pzlev_config(domain)%plevels%values)) THEN
      CALL ser_component(o, "plevels%nvalues", nh_pzlev_config(domain)%plevels%nvalues)
      CALL ser_component(o, "plevels%values", nh_pzlev_config(domain)%plevels%values)
    ENDIF
    IF(ALLOCATED(nh_pzlev_config(domain)%zlevels%values)) THEN
      CALL ser_component(o, "zlevels%nvalues", nh_pzlev_config(domain)%zlevels%nvalues)
      CALL ser_component(o, "zlevels%values", nh_pzlev_config(domain)%zlevels%values)
    ENDIF
    IF(ALLOCATED(nh_pzlev_config(domain)%ilevels%values)) THEN
      CALL ser_component(o, "ilevels%nvalues", nh_pzlev_config(domain)%ilevels%nvalues)
      CALL ser_component(o, "ilevels%values", nh_pzlev_config(domain)%ilevels%values)
    ENDIF
    IF(ASSOCIATED(nh_pzlev_config(domain)%p3d) .AND. SIZE(nh_pzlev_config(domain)%p3d) > 0) &
      CALL ser_component(o, "p3d", nh_pzlev_config(domain)%p3d)
    IF(ASSOCIATED(nh_pzlev_config(domain)%z3d) .AND. SIZE(nh_pzlev_config(domain)%z3d) > 0) &
      CALL ser_component(o, "z3d", nh_pzlev_config(domain)%z3d)
    IF(ASSOCIATED(nh_pzlev_config(domain)%i3d) .AND. SIZE(nh_pzlev_config(domain)%i3d) > 0) &
      CALL ser_component(o, "i3d", nh_pzlev_config(domain)%i3d)

  END SUBROUTINE ser_manually

! supporting subroutines:

  SUBROUTINE ser_test_scalars()
    ! test openACC update mechanism for scalars
    !
    ! The main objective of this test is to make sure, that references of the
    ! scalar components in p_test are passed to ser_scalar_* such that the
    ! ACC UPDATE inside ser_scalar_* affects the components of p_test.
    TYPE(t_ser_options) :: o

    IF(lscalars_test_was_run) return

    o%abs_threshold = 8
    o%rel_threshold = 8
    o%lopenacc = .TRUE.

    p_test%int = 23
    p_test%float = 3.14
    p_test%double = 3.14159
    p_test%log_a = .TRUE.
    p_test%log_b = .FALSE.

    !$ACC DATA COPYIN(p_test)

#ifdef SERIALIZE_READ_REFERENCE
    ! manipulate GPU
    !$ACC SERIAL DEFAULT(PRESENT)
    p_test%int = 42
    p_test%float = 2.71
    p_test%double = 2.71828
    p_test%log_a = .FALSE.
    p_test%log_b = .TRUE.
    !$ACC END SERIAL

    o%ser_mode = 3 ! compare

    ! These should report 5 differences (which will be captured by the __TEST__ mode)
    CALL ser_component(o, "__TEST__error__p_test%int",    p_test%int)
    CALL ser_component(o, "__TEST__error__p_test%float",  p_test%float)
    CALL ser_component(o, "__TEST__error__p_test%double", p_test%double)
    CALL ser_component(o, "__TEST__error__p_test%log_a",  p_test%log_a)
    CALL ser_component(o, "__TEST__error__p_test%log_b",  p_test%log_b)

    o%ser_mode = 1 ! read serialized values
#else
    o%ser_mode = 0 ! write
#endif
    CALL ser_component(o, "__TEST__error__p_test%int",    p_test%int)
    CALL ser_component(o, "__TEST__error__p_test%float",  p_test%float)
    CALL ser_component(o, "__TEST__error__p_test%double", p_test%double)
    CALL ser_component(o, "__TEST__error__p_test%log_a",  p_test%log_a)
    CALL ser_component(o, "__TEST__error__p_test%log_b",  p_test%log_b)

#ifdef SERIALIZE_READ_REFERENCE
    o%ser_mode = 3 ! compare values, that were read in, to reference
#endif

    ! These should report no differences
    CALL ser_component(o, "__TEST__noerror__p_test%int",    p_test%int)
    CALL ser_component(o, "__TEST__noerror__p_test%float",  p_test%float)
    CALL ser_component(o, "__TEST__noerror__p_test%double", p_test%double)
    CALL ser_component(o, "__TEST__noerror__p_test%log_a",  p_test%log_a)
    CALL ser_component(o, "__TEST__noerror__p_test%log_b",  p_test%log_b)

    !$ACC END DATA

#ifdef SERIALIZE_READ_REFERENCE
    ! Turn this test off after it was run in read-and-compare mode.
    ! This test is always run in write mode to make sure that the comparison
    ! find compare values for any save point.
    lscalars_test_was_run = .TRUE.
#endif

  END SUBROUTINE ser_test_scalars

#endif

END MODULE mo_ser_manually
