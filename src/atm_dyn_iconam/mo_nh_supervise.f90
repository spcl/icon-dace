!
! Computes total integrals and maxwinds in the nonhydrostatic model.
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
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nh_supervise

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_nonhydro_types,      ONLY: t_nh_state, t_nh_prog, t_nh_diag
  USE mo_model_domain,        ONLY: t_patch
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_cell
  USE mo_grid_config,         ONLY: l_limited_area, grid_sphere_radius
  USE mo_math_constants,      ONLY: pi
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_run_config,          ONLY: dtime, msg_level, output_mode,           &
    &                               ltransport, ntracer, lforcing, iforcing, &
    &                               iqm_max
  USE mo_impl_constants,      ONLY: SUCCESS, inwp, iaes, &
    &                               min_rlcell_int, min_rledge_int, iheldsuarez
  USE mo_physical_constants,  ONLY: cvd
  USE mo_mpi,                 ONLY: my_process_is_stdio, get_my_mpi_all_id, &
    &                               process_mpi_stdio_id
  USE mo_io_units,            ONLY: find_next_free_unit
  USE mo_sync,                ONLY: global_sum_array, global_max
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_fortran_tools,       ONLY: init, set_acc_host_or_device, assert_acc_device_only
  USE mo_nh_diagnose_pres_temp,ONLY: calc_qsum

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname   = 'mo_nh_supervise'

  ! Needed by supervise_total_integrals_nh to keep data between steps
  REAL(wp), ALLOCATABLE, SAVE :: z_tracer_mass_0(:)      ! tracer specific total mass at first step

  INTEGER :: n_file_ti = -1, n_file_tti = -1,  check_total_quant_fileid = -1       ! file identifiers

  ! --- Print-out of max winds to an ASCII file.
  !     This requires namelist setting 'run_nml::output = "maxwinds"'
  CHARACTER(len=*), PARAMETER :: maxwinds_filename = "maxwinds.log"   !< file name
  INTEGER                     :: maxwinds_funit                       !< file unit

  PUBLIC :: init_supervise_nh
  PUBLIC :: finalize_supervise_nh
  PUBLIC :: supervise_total_integrals_nh
  PUBLIC :: print_maxwinds
  PUBLIC :: compute_dpsdt
  PUBLIC :: compute_hcfl


CONTAINS

  !-----------------------------------------------------------------------------
  !> init_supervise_nh
  !
  !  Initialization routine for this module (eg. opening of files)
  !
  SUBROUTINE init_supervise_nh( )
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::init_supervise_nh"
    INTEGER :: istat

    ! --- Print-out of max winds to an ASCII file.
    !     This requires namelist setting 'run_nml::output = "maxwinds"'
    IF (output_mode%l_maxwinds .AND. my_process_is_stdio()) THEN
      maxwinds_funit = find_next_free_unit(100,1000)
      CALL message(routine,"Open log file "//maxwinds_filename//" for writing.")
      OPEN(UNIT=maxwinds_funit, FILE=maxwinds_filename, ACTION="write", &
        &  FORM='FORMATTED', IOSTAT=istat)
      IF (istat/=SUCCESS) &
        &  CALL finish(routine, 'could not open '//maxwinds_filename)
    END IF
  END SUBROUTINE init_supervise_nh


  !-----------------------------------------------------------------------------
  !> finalize_supervise_nh
  !
  !  Clean-up routine for this module (eg. closing of files)
  !
  SUBROUTINE finalize_supervise_nh( )
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::finalize_supervise_nh"
    INTEGER :: istat

    ! --- close max winds ASCII file.
    IF (output_mode%l_maxwinds .AND. my_process_is_stdio()) THEN
      CLOSE(maxwinds_funit, IOSTAT=istat)
      IF (istat/=SUCCESS) &
        &  CALL finish(routine,'could not close '//maxwinds_filename)
    END IF
  END SUBROUTINE finalize_supervise_nh
  

  !-----------------------------------------------------------------------------
  !! supervise_total_integrals_nh
  
  SUBROUTINE supervise_total_integrals_nh( k_step, patch, nh_state, int_state, ntimlev, ntimlev_rcf, l_last_step, lacc)

    INTEGER,                  INTENT(IN) :: k_step            ! actual time step
    TYPE(t_patch),            INTENT(IN) :: patch             ! Patch
    TYPE(t_nh_state), TARGET, INTENT(IN) :: nh_state          ! NH State
    TYPE(t_int_state),        INTENT(IN) :: int_state         ! Interpolation State
    INTEGER,                  INTENT(IN) :: ntimlev           ! time level
    INTEGER,                  INTENT(IN) :: ntimlev_rcf       ! rcf time level
    LOGICAL,                  INTENT(IN) :: l_last_step
    LOGICAL, OPTIONAL,        INTENT(IN) :: lacc              ! If true, use openacc

    REAL(wp), SAVE :: z_total_mass_0    !< total air mass including vapor and condensate at first time step [kg]
    REAL(wp), SAVE :: z_dry_mass_0      !< dry air mass at first time step [kg]
    REAL(wp), SAVE :: z_total_energy_0  !< total energy at first time step [kg]
    !
    REAL(wp) :: z_total_mass            !< total air mass, including vapor and condensate partial densities [kg]
    REAL(wp) :: z_dry_mass              !< dry air mass (i.e. excluding vapor and condensate) [kg]
    REAL(wp) :: z_kin_energy            !< kinetic energy [Nm]
    REAL(wp) :: z_pot_energy            !< potential energy [Nm]
    REAL(wp) :: z_int_energy            !< internal energy [Nm]
    REAL(wp) :: z_total_energy          !< kinetic + potential + internal energy
    REAL(wp) :: z_mean_surfp            !< mean surface pressure [Pa]
    !
    REAL(wp) :: z_total_mass_re         !< changes in total mass compared to first step
    REAL(wp) :: z_dry_mass_re           !< changes in dry mass compared to first step
    REAL(wp) :: z_total_energy_re       !< changes in total energy compared to first step
    !
    REAL(wp) :: z_kin_energy_re         !< percentage of kinetic energy out of total energy
    REAL(wp) :: z_int_energy_re         !< percentage of internal energy out of total energy
    REAL(wp) :: z_pot_energy_re         !< percentage of potential energy out of total energy
    REAL(wp) :: z_volume
#ifndef NOMPI
    REAL(wp) ::   z_total_mass_2d(nproma,patch%nblks_c), &
      &           z_dry_mass_2d(nproma,patch%nblks_c),   &
      &           z_kin_energy_2d(nproma,patch%nblks_c), &
      &           z_int_energy_2d(nproma,patch%nblks_c), &
      &           z_pot_energy_2d(nproma,patch%nblks_c), &
      &           z_surfp_2d(nproma,patch%nblks_c)

#endif

    ! tracer fields
    REAL(wp):: z_tracer_mass(ntracer)       ! tracer specific total mass
    REAL(wp):: z_tracer_mass_re(ntracer)    ! changes in tracer specific total mass compared to first step
    REAL(wp):: z_aux_tracer(nproma,patch%nblks_c)
    REAL(wp):: z_total_tracer_mass          ! total tracer mass
    REAL(wp):: z_total_tracer_mass_re       ! changes in total tracer mass compared to first step
    REAL(wp):: z_water_mass                 ! total mass of water (including vapor)
    REAL(wp):: z_water_mass_re              ! changes in total water mass compared to first step
    !
    REAL(wp), SAVE :: z_total_tracer_mass_0 ! total tracer mass at first step
    REAL(wp), SAVE :: z_water_mass_0        ! total mass of water at first step

    INTEGER :: jb, jk, jc, jt               ! loop indices
    INTEGER :: nlen, npromz_c, nblks_c
    INTEGER :: nlev                           ! number of full levels
    INTEGER :: ist                            ! status variable
    TYPE(t_nh_prog), POINTER :: prog       ! prog state
    TYPE(t_nh_prog), POINTER :: prog_rcf   ! prog_rcf state
    TYPE(t_nh_diag), POINTER :: diag       ! diag state
    REAL(wp) :: z_qsum(nproma,patch%nlev,patch%nblks_c)  ! total condensate including vapour
    REAL(wp) :: z_elapsed_time
    REAL(wp) :: z_ekin(nproma,patch%nlev,patch%nblks_c)

    REAL(wp) :: max_vn, max_w
    INTEGER  :: max_vn_level, max_vn_process, max_w_level, max_w_process
    INTEGER  :: water_tracer_list(iqm_max)    ! list of all water tracer IDs (including vapor)

    CHARACTER(*), PARAMETER :: routine = modname//"::supervise_total_integrals_nh"
    !-----------------------------------------------------------------------------

    CALL assert_acc_device_only(routine, lacc)

    ! store IDs of all water tracers in a list, including vapor
    IF ( lforcing .AND. iforcing /= iheldsuarez ) THEN
      water_tracer_list = (/(jt, jt=1,iqm_max)/)
    ENDIF

    !$ACC DATA CREATE(z_ekin, z_qsum, z_aux_tracer) COPYIN(water_tracer_list)
#ifndef NOMPI
    !$ACC DATA CREATE(z_total_mass_2d, z_dry_mass_2d, z_kin_energy_2d, z_int_energy_2d, z_pot_energy_2d, z_surfp_2d) COPYIN(water_tracer_list)
#endif

    IF (.NOT. ALLOCATED (z_tracer_mass_0)) THEN
      ALLOCATE (z_tracer_mass_0(ntracer), STAT=ist)
      IF(ist/=SUCCESS)THEN
        CALL finish (routine, 'allocation of z_tracer_mass_0 failed')
      ENDIF
      z_tracer_mass_0 = 0.0_wp
    END IF


    z_elapsed_time = dtime*REAL(k_step,wp)/3600.0_wp

    prog     => nh_state%prog(ntimlev)
    prog_rcf => nh_state%prog(ntimlev_rcf)
    diag     => nh_state%diag

    nblks_c   = patch%nblks_c
    npromz_c  = patch%npromz_c

    ! number of vertical levels
    nlev = patch%nlev

    IF (iforcing <= 1) THEN ! u and v are not diagnosed regularly if physics is turned off
      CALL rbf_vec_interpol_cell(prog%vn,patch,int_state,diag%u,diag%v, opt_acc_async=.TRUE.)
    ENDIF


!$OMP PARALLEL
    CALL init(z_qsum, opt_acc_async=.TRUE.)
!$OMP BARRIER

!$OMP DO PRIVATE(jb,jk,jc,nlen)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO jc = 1, nlen
          z_ekin(jc,jk,jb) = 0.5_wp*(diag%u(jc,jk,jb)**2 + diag%v(jc,jk,jb)**2) +  &
            0.25_wp*(prog%w(jc,jk,jb)**2 + prog%w(jc,jk+1,jb)**2)
        ENDDO
      ENDDO
      !$ACC END PARALLEL

      IF ( ltransport .AND. lforcing .AND. iforcing /= iheldsuarez ) THEN
        ! compute total condensate INCLUDING vapor
        !
        CALL calc_qsum (tracer      = prog_rcf%tracer(:,:,:,:), & !in
          &             qsum        = z_qsum(:,:,jb),           & !inout
          &             tracer_list = water_tracer_list(:),     & !in
          &             jb          = jb,                       & !in
          &             i_startidx  = 1,                        & !in
          &             i_endidx    = nlen,                     & !in
          &             slev        = 1,                        & !in
          &             slev_moist  = 1,                        & !in
          &             nlev        = nlev)                       !in
      ENDIF

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

#ifdef NOMPI

#ifdef _OPENACC
    CALL finish(routine, "The NOMPI has not been tested with OpenACC.")
#endif

    z_total_mass  = 0.0_wp
    z_dry_mass    = 0.0_wp
    z_kin_energy  = 0.0_wp
    z_pot_energy  = 0.0_wp
    z_int_energy  = 0.0_wp
    z_mean_surfp  = 0.0_wp
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(z_volume) &
      !$ACC   REDUCTION(+, z_total_mass, z_kin_energy, z_int_energy, z_pot_energy, z_dry_mass)
      DO jk = 1, nlev
        DO jc = 1, nlen
          z_volume = patch%cells%area(jc,jb)*nh_state%metrics%ddqz_z_full(jc,jk,jb) &
            &       /REAL(patch%n_patch_cells_g,wp)*nh_state%metrics%deepatmo_vol_mc(jk)
          z_total_mass = z_total_mass + &
            &              prog%rho(jc,jk,jb)*z_volume
          z_dry_mass   = z_dry_mass +   &
            &              prog%rho(jc,jk,jb)*(1._wp - z_qsum(jc,jk,jb))*z_volume
          z_kin_energy = z_kin_energy + &
            &              prog%rho(jc,jk,jb)*z_ekin(jc,jk,jb)*z_volume
          z_int_energy = z_int_energy + &
            &              cvd*prog%exner(jc,jk,jb)*prog%rho(jc,jk,jb)*prog%theta_v(jc,jk,jb)*z_volume
          z_pot_energy = z_pot_energy + &
            &              prog%rho(jc,jk,jb)*nh_state%metrics%geopot(jc,jk,jb)*z_volume
        ENDDO
      ENDDO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) REDUCTION(+, z_mean_surfp)
      DO jc = 1, nlen
        z_mean_surfp = z_mean_surfp + diag%pres_sfc(jc,jb)*patch%cells%area(jc,jb) /  &
          &  (4._wp*grid_sphere_radius**2*pi)
      ENDDO
      !$ACC END PARALLEL
      !$ACC WAIT
    ENDDO
    z_total_energy = z_int_energy+z_kin_energy+z_pot_energy

#else

!$OMP PARALLEL
    CALL init(z_total_mass_2d, lacc=.TRUE., opt_acc_async=.TRUE.)
    CALL init(z_dry_mass_2d, lacc=.TRUE., opt_acc_async=.TRUE.)
    CALL init(z_kin_energy_2d, lacc=.TRUE., opt_acc_async=.TRUE.)
    CALL init(z_int_energy_2d, lacc=.TRUE., opt_acc_async=.TRUE.)
    CALL init(z_pot_energy_2d, lacc=.TRUE., opt_acc_async=.TRUE.)
    CALL init(z_surfp_2d, lacc=.TRUE., opt_acc_async=.TRUE.)
!$OMP BARRIER

!$OMP DO PRIVATE(jb,jk,jc,nlen,z_volume)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk = 1,nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(z_volume)
        DO jc = 1, nlen
          z_volume = patch%cells%area(jc,jb)*nh_state%metrics%ddqz_z_full(jc,jk,jb) &
            & /REAL(patch%n_patch_cells_g,wp)*nh_state%metrics%deepatmo_vol_mc(jk)
          z_total_mass_2d(jc,jb) = z_total_mass_2d(jc,jb)&
            & +prog%rho(jc,jk,jb)*z_volume
          z_dry_mass_2d(jc,jb) = z_dry_mass_2d(jc,jb)&
            & +prog%rho(jc,jk,jb)*(1._wp-z_qsum(jc,jk,jb))*z_volume
          z_kin_energy_2d(jc,jb) = z_kin_energy_2d(jc,jb)&
            & +prog%rho(jc,jk,jb)*z_ekin(jc,jk,jb)*z_volume
          z_int_energy_2d(jc,jb) = z_int_energy_2d(jc,jb)&
            & +cvd*prog%exner(jc,jk,jb)*prog%rho(jc,jk,jb)*prog%theta_v(jc,jk,jb)*z_volume
          z_pot_energy_2d(jc,jb) = z_pot_energy_2d(jc,jb)&
            & +prog%rho(jc,jk,jb)*nh_state%metrics%geopot(jc,jk,jb)*z_volume
        ENDDO
      ENDDO
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = 1, nlen
        z_surfp_2d(jc,jb) = diag%pres_sfc(jc,jb)*patch%cells%area(jc,jb)
        IF(.NOT. patch%cells%decomp_info%owner_mask(jc,jb)) THEN
          z_total_mass_2d(jc,jb) = 0._wp
          z_dry_mass_2d(jc,jb)   = 0._wp
          z_kin_energy_2d(jc,jb) = 0._wp
          z_int_energy_2d(jc,jb) = 0._wp
          z_pot_energy_2d(jc,jb) = 0._wp
          z_surfp_2d(jc,jb)      = 0._wp
        ENDIF
      ENDDO
      !$ACC END PARALLEL
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    z_mean_surfp   = global_sum_array( z_surfp_2d, lacc=.TRUE. )/(4._wp*grid_sphere_radius**2*pi)
    z_total_mass   = global_sum_array( z_total_mass_2d, lacc=.TRUE. )
    z_dry_mass     = global_sum_array( z_dry_mass_2d, lacc=.TRUE. )
    z_kin_energy   = global_sum_array( z_kin_energy_2d, lacc=.TRUE. )
    z_int_energy   = global_sum_array( z_int_energy_2d, lacc=.TRUE. )
    z_pot_energy   = global_sum_array( z_pot_energy_2d, lacc=.TRUE. )
    z_total_energy = z_int_energy+z_kin_energy+z_pot_energy

#endif

    IF (k_step == 1) THEN
      z_total_mass_0   = z_total_mass
      z_dry_mass_0     = z_dry_mass
      z_total_energy_0 = z_total_energy

      ! Open the datafile and store initial values of mass and energy
      IF (my_process_is_stdio()) THEN
        CALL open_total_integral_files(z_total_mass_0, z_dry_mass_0, z_total_energy_0)
      ENDIF
    ENDIF

    ! percentage of single energies out of the total energy
    z_kin_energy_re = z_kin_energy/z_total_energy*100.0_wp
    z_int_energy_re = z_int_energy/z_total_energy*100.0_wp
    z_pot_energy_re = z_pot_energy/z_total_energy*100.0_wp
    ! changes compared to first step
    z_total_mass_re   = z_total_mass   /z_total_mass_0    -1.0_wp
    z_dry_mass_re     = z_dry_mass     /z_dry_mass_0      -1.0_wp
    z_total_energy_re = z_total_energy /z_total_energy_0  -1.0_wp

    IF (my_process_is_stdio()) THEN
      if (n_file_ti >= 0) &
        WRITE(n_file_ti,'(i8,7e20.12)') &
        &   k_step, z_total_mass_re, z_dry_mass_re, z_total_energy_re, z_kin_energy_re, &
        &   z_int_energy_re, z_pot_energy_re, z_mean_surfp
      IF (l_last_step) THEN
        if (n_file_ti >= 0) &
          CLOSE(n_file_ti)
      ENDIF
    ENDIF


    IF ( ltransport .OR. ANY((/inwp,iaes/)==iforcing) ) THEN

      z_tracer_mass(:) = 0.0_wp ! init must not used be used (this is no GPU variable)
      z_water_mass = 0.0_wp

      DO jt=1, ntracer
!$OMP PARALLEL
        CALL init(z_aux_tracer(:,:), lacc=.TRUE., opt_acc_async=.TRUE.) ! reinitialize for each jt
!$OMP BARRIER

!$OMP DO PRIVATE(jb,jk,jc,nlen,z_volume)
        DO jb = 1, nblks_c
          IF (jb /= nblks_c) THEN
            nlen = nproma
          ELSE
            nlen = npromz_c
          ENDIF

          ! compute tracer mass in each vertical column
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
          !$ACC LOOP SEQ
          DO jk = 1, nlev
            !$ACC LOOP GANG VECTOR PRIVATE(z_volume)
            DO jc = 1, nlen
              z_volume = patch%cells%area(jc,jb)             &
                &    * nh_state%metrics%ddqz_z_full(jc,jk,jb) &
                &    * prog%rho(jc,jk,jb) * nh_state%metrics%deepatmo_vol_mc(jk)

              z_aux_tracer(jc,jb) = z_aux_tracer(jc,jb)    &
                &    + prog_rcf%tracer(jc,jk,jb,jt) * z_volume
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDDO
!$OMP END DO

!$OMP DO PRIVATE(jb,jc)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jb = 1, nblks_c
          DO jc = 1, nproma
            IF(.NOT.patch%cells%decomp_info%owner_mask(jc,jb)) z_aux_tracer(jc,jb) = 0._wp
          ENDDO
        ENDDO
        !$ACC END PARALLEL
!$OMP END DO
!$OMP END PARALLEL
        z_tracer_mass(jt) = global_sum_array(z_aux_tracer(:,:), lacc=.TRUE.)
      ENDDO  ! ntracer

      z_total_tracer_mass = SUM(z_tracer_mass(1:ntracer))
      !
      IF ( ANY((/inwp,iaes/)==iforcing) ) THEN
        ! when run in physics mode, compute total water mass (including vapor)
        z_water_mass = SUM(z_tracer_mass((/water_tracer_list/)))
      ENDIF

      ! Save total tracer mass at first time step
      IF (k_step == 1) THEN
        z_tracer_mass_0(:)     = z_tracer_mass(:)
        z_total_tracer_mass_0  = z_total_tracer_mass
        z_water_mass_0         = z_water_mass
      ELSE
        ! compute mass change relative to time step 1
        !
        IF (z_total_tracer_mass_0 == 0._wp) THEN
          z_total_tracer_mass_re = 0._wp
        ELSE
          z_total_tracer_mass_re = z_total_tracer_mass /z_total_tracer_mass_0 -1._wp
        ENDIF
        IF (z_water_mass_0 == 0._wp) THEN
          z_water_mass_re = 0._wp
        ELSE
          z_water_mass_re = z_water_mass /z_water_mass_0 -1._wp
        ENDIF
        !
        DO jt=1,ntracer
          IF (z_tracer_mass_0(jt) == 0._wp) THEN
            z_tracer_mass_re(jt) = 0._wp
          ELSE
            z_tracer_mass_re(jt) = z_tracer_mass(jt) /z_tracer_mass_0(jt) - 1._wp
          ENDIF
        ENDDO
      ENDIF

      ! write to file
      IF (my_process_is_stdio()) THEN
        IF (n_file_tti > 0) THEN
           ! total tracer mass
           WRITE(n_file_tti,'(i8,f22.8,a22,2e22.12)')             &
           k_step, z_elapsed_time, 'tot', z_total_tracer_mass, z_total_tracer_mass_re
           ! water mass
           WRITE(n_file_tti,'(i8,f22.8,a22,2e22.12)')             &
           k_step, z_elapsed_time, 'water', z_water_mass, z_water_mass_re
           !
          DO jt=1,ntracer
            WRITE(n_file_tti,'(i8,f22.8,i22,2e22.12)')          &
            k_step, z_elapsed_time, jt, z_tracer_mass(jt), z_tracer_mass_re(jt)
          ENDDO
        ENDIF
      ENDIF

      IF (l_last_step) THEN
        IF (n_file_ti >= 0) CLOSE(n_file_ti)
        IF (n_file_tti > 0) CLOSE(n_file_tti)

        DEALLOCATE(z_tracer_mass_0, STAT=ist)
        IF(ist/=SUCCESS)THEN
          CALL finish (routine, 'deallocation of z_tracer_mass_0 failed')
        ENDIF
      ENDIF

    ENDIF    ! ltransport


    ! write additional check quantities (check_global_quantities)
    CALL calculate_maxwinds(patch, prog%vn, prog%w, &
        & max_vn, max_vn_level, max_vn_process,        &
        & max_w, max_w_level, max_w_process, lacc=.TRUE. )

    IF (my_process_is_stdio()) THEN
      IF (check_total_quant_fileid >= 0) THEN
        WRITE(check_total_quant_fileid,'(i8,2e20.12)') &
          &  k_step, max_vn, max_w
        IF (l_last_step) THEN
          CLOSE(check_total_quant_fileid)
          check_total_quant_fileid = -1
        ENDIF
      ENDIF
    ENDIF

    !$ACC WAIT
#ifndef NOMPI
    !$ACC END DATA ! z_total_mass_2d, z_kin_energy_2d, z_int_energy_2d, z_pot_energy_2d, z_surfp_2d, z_dry_mass_2d
#endif
    !$ACC END DATA

  END SUBROUTINE supervise_total_integrals_nh
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE open_total_integral_files(total_mass_0, dry_mass_0, total_energy_0)

    REAL(wp), INTENT(IN) :: total_mass_0     !< total mass intial value
    REAL(wp), INTENT(IN) :: dry_mass_0       !< dry mass initial value
    REAL(wp), INTENT(IN) :: total_energy_0   !< total energy initial value
    INTEGER :: istat

    n_file_ti = find_next_free_unit(100,1000)
    ! write(0,*) "n_file_ti", n_file_ti
    OPEN(UNIT=n_file_ti,FILE='total_integrals.dat',ACTION="write", &
         FORM='FORMATTED',IOSTAT=istat)
    IF (istat/=SUCCESS) THEN
      CALL finish('supervise_total_integrals_nh','could not open total_integrals.dat')
    ENDIF
    WRITE (n_file_ti,'(A8,3A20)')'TIMESTEP',&
      '      total mass m0,',&
      '     dry mass mdry0,',&
      '     total energy e0'
    WRITE(n_file_ti,'(i8,3e20.12)') &
      &   (/1/), total_mass_0, dry_mass_0, total_energy_0
    WRITE(n_file_ti,*)' '

    WRITE (n_file_ti,'(A8,7A20)')'TIMESTEP',&
      '            m/m0 -1,',&
      '      mdry/mdry0 -1,',&
      '            e/e0 -1,',&
      '             % kine,',&
      '             % inne,',&
      '             % pote,',&
      '    mean surf press.'

    ! Open the datafile for tracer diagnostic
    IF (ltransport .OR. ( iforcing == inwp ) .OR. ( iforcing == iaes )) THEN

      n_file_tti = find_next_free_unit(100,1000)
      OPEN(UNIT=n_file_tti,FILE='tracer_total_integrals.dat',ACTION="write", &
           FORM='FORMATTED',IOSTAT=istat)
      IF (istat/=SUCCESS) THEN
        CALL finish('supervise_total_integrals_nh','could not open tracer_total_integrals.dat')
      ENDIF
      WRITE (n_file_tti,'(A8,4A22)')'TIMESTEP',&
        '    ELAPSED TIME (hr),',&
        '        TRACER NR (#),',&
        '     TRACER MASS (kg),',&
        '                m/m0-1'
    ENDIF

    check_total_quant_fileid = find_next_free_unit(100,1000)
    OPEN(UNIT=check_total_quant_fileid,FILE='check_global_quantities.dat', &
      ACTION="write",FORM='FORMATTED',IOSTAT=istat)
    IF (istat/=SUCCESS) THEN
      CALL finish('supervise_total_integrals_nh','could not open check_global_quantities.dat')
    ENDIF
    WRITE (check_total_quant_fileid, '(A8,2A20)') &
      'TIMESTEP',           &
      '            Max vn,',&
      '            Max w'

  END SUBROUTINE open_total_integral_files

  !-------------------------------------------------------------------------
  !>
  !! Computation of maximum horizontal and vertical wind speed for runtime diagnostics
  !! Was included in mo_nh_stepping before
  !!
  SUBROUTINE print_maxwinds(patch, vn, w, lacc)

    TYPE(t_patch), INTENT(IN) :: patch    ! Patch
    REAL(wp),      INTENT(IN) :: vn(:,:,:), w(:,:,:) ! horizontal and vertical wind speed
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    ! local variables
    REAL(wp) :: max_vn, max_w
    INTEGER  :: max_vn_level, max_vn_process, max_w_level, max_w_process
    LOGICAL :: lzacc ! non-optional version of lacc

    !-----------------------------------------------------------------------
    CALL set_acc_host_or_device(lzacc, lacc)

    CALL calculate_maxwinds(patch, vn, w,       &
        & max_vn, max_vn_level, max_vn_process, &
        & max_w, max_w_level, max_w_process, lacc=lzacc )


    !--- Get max over all PEs
    IF (msg_level >= 7) THEN

      ! print a detailed information on global maxima:
      ! containing the process ID and the level where the
      ! maximum occurred.

      IF (msg_level >= 13) THEN
        WRITE(message_text,'(a,i3,a,2(e18.10,a,i5,a,i3,a))') 'MAXABS VN, W in domain', patch%id, ':', &
          & max_vn, " (on proc #", max_vn_process, ", level ", max_vn_level, "), ", &
          & max_w,  " (on proc #", max_w_process,  ", level ", max_w_level,  "), "
      ELSE
        WRITE(message_text,'(a,i3,a,2(e18.10,a,i3,a))') 'MAXABS VN, W in domain', patch%id, ':', &
          & max_vn, " at level ",  max_vn_level, ", ", &
          & max_w,  " at level ",  max_w_level,  ", "
      END IF

    ELSE

      ! on PE0 print a short information on global maxima:
      WRITE(message_text,'(a,2e18.10)') 'MAXABS VN, W ', max_vn, max_w

    END IF
    CALL message('',message_text)

    ! --- Print-out of max winds to an ASCII file.
    ! 
    !     This requires namelist setting 'run_nml::output = "maxwinds"'

    IF (output_mode%l_maxwinds .AND. my_process_is_stdio()) THEN
      WRITE(maxwinds_funit,'(a,i3,a,2(e18.10,a,i3,a))') 'MAXABS VN, W in domain', patch%id, ':', &
        & max_vn, " at level ",  max_vn_level, ", ", &
        & max_w,  " at level ",  max_w_level,  ", "
    END IF

  END SUBROUTINE print_maxwinds

  !-------------------------------------------------------------------------
  !>
  !! Computation of maximum horizontal and vertical wind speed for runtime diagnostics
  !! Was included in mo_nh_stepping before
  !!
  SUBROUTINE calculate_maxwinds(patch, vn, w,   &
        & max_vn, max_vn_level, max_vn_process, &
        & max_w, max_w_level, max_w_process, lacc )

    TYPE(t_patch), INTENT(IN)  :: patch    ! Patch
    REAL(wp),      INTENT(IN)  :: vn(:,:,:), w(:,:,:) ! horizontal and vertical wind speed
    REAL(wp),      INTENT(OUT) :: max_vn, max_w
    INTEGER,       INTENT(OUT) :: max_vn_level, max_vn_process, max_w_level, max_w_process
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    ! local variables
    REAL(wp) :: vn_aux(patch%edges%end_blk(min_rledge_int,MAX(1,patch%n_childdom)),patch%nlev)
    REAL(wp) :: w_aux (patch%cells%end_blk(min_rlcell_int,MAX(1,patch%n_childdom)),patch%nlevp1)
    REAL(wp) :: vn_aux_lev(patch%nlev), w_aux_lev(patch%nlevp1), vmax(2), vn_aux_tmp, w_aux_tmp

    INTEGER  :: istartblk_c, istartblk_e, iendblk_c, iendblk_e, i_startidx, i_endidx
    INTEGER  :: jb, jk, jg
#if defined( __INTEL_COMPILER ) || defined( _OPENACC ) || defined (__SX__)
    INTEGER  :: jec
#endif
    INTEGER  :: proc_id(2), keyval(2)
    LOGICAL :: lzacc ! non-optional version of lacc

    !-----------------------------------------------------------------------
    CALL set_acc_host_or_device(lzacc, lacc)

    istartblk_c = patch%cells%start_block(grf_bdywidth_c+1)
    istartblk_e = patch%edges%start_block(grf_bdywidth_e+1)
    iendblk_c   = patch%cells%end_block(min_rlcell_int)
    iendblk_e   = patch%edges%end_block(min_rledge_int)
    jg          = patch%id

    !$ACC DATA PRESENT(vn, w, patch) COPYOUT(vn_aux, w_aux) IF(lzacc)

    IF (jg > 1 .OR. l_limited_area) THEN
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      vn_aux(1:MIN(istartblk_e,iendblk_e),:) = 0._wp
      w_aux (1:MIN(istartblk_c,iendblk_c),:) = 0._wp
      !$ACC END KERNELS
    ENDIF

!$OMP PARALLEL
#if defined( __INTEL_COMPILER ) || defined (__SX__)
!$OMP DO PRIVATE(jb, jk, jec, i_startidx, i_endidx, vn_aux_tmp) ICON_OMP_DEFAULT_SCHEDULE
#else
!$OMP DO PRIVATE(jb, jk, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
#endif
    DO jb = istartblk_e, iendblk_e

      CALL get_indices_e(patch, jb, istartblk_e, iendblk_e, i_startidx, i_endidx, &
                         grf_bdywidth_e+1, min_rledge_int)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG PRIVATE(vn_aux_tmp)
!$NEC novector
      DO jk = 1, patch%nlev
#if defined( __INTEL_COMPILER ) || defined( _OPENACC ) || defined (__SX__)
        vn_aux_tmp = 0._wp
        !$ACC LOOP VECTOR REDUCTION(MAX: vn_aux_tmp)
        DO jec = i_startidx,i_endidx
          vn_aux_tmp = MAX(vn_aux_tmp, -vn(jec,jk,jb), vn(jec,jk,jb))
        ENDDO
        vn_aux(jb,jk) = vn_aux_tmp
#else
        vn_aux(jb,jk) = MAXVAL(ABS(vn(i_startidx:i_endidx,jk,jb)))
#endif
      ENDDO
      !$ACC END PARALLEL
    END DO
!$OMP END DO

#if defined( __INTEL_COMPILER ) || defined (__SX__)
!$OMP DO PRIVATE(jb, jk, jec, i_startidx, i_endidx, w_aux_tmp) ICON_OMP_DEFAULT_SCHEDULE
#else
!$OMP DO PRIVATE(jb, jk, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
#endif
    DO jb = istartblk_c, iendblk_c

      CALL get_indices_c(patch, jb, istartblk_c, iendblk_c, i_startidx, i_endidx, &
                         grf_bdywidth_c+1, min_rlcell_int)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG PRIVATE(w_aux_tmp)
!$NEC novector
      DO jk = 1, patch%nlevp1
#if defined( __INTEL_COMPILER ) || defined( _OPENACC ) || defined (__SX__)
        w_aux_tmp = 0._wp
        !$ACC LOOP VECTOR REDUCTION(MAX: w_aux_tmp)
        DO jec = i_startidx,i_endidx
          w_aux_tmp = MAX(w_aux_tmp, -w(jec,jk,jb), w(jec,jk,jb))
        ENDDO
        w_aux(jb,jk) = w_aux_tmp
#else
        w_aux(jb,jk) = MAXVAL(ABS(w(i_startidx:i_endidx,jk,jb)))
#endif
      ENDDO
      !$ACC END PARALLEL
    END DO
!$OMP END DO

    !$ACC WAIT
    !$ACC END DATA

! At this point vn_aux and w_aux reside on the host.  
! Avoid doing MAXVAL with OpenACC -- this is not well supported!
#ifndef __SX__
!$OMP DO PRIVATE(jk) ICON_OMP_DEFAULT_SCHEDULE

    DO jk = 1, patch%nlev
      vn_aux_lev(jk) = MAXVAL(vn_aux(:,jk))
      w_aux_lev (jk) = MAXVAL(w_aux(:,jk))
    END DO

!$OMP END DO
#else
    vn_aux_lev = 0._wp
    w_aux_lev  = 0._wp
!$OMP DO PRIVATE(jb,jk) REDUCTION(max:vn_aux_lev) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = istartblk_e, iendblk_e
      DO jk = 1, patch%nlev
        vn_aux_lev(jk) = MAX(vn_aux_lev(jk),vn_aux(jb,jk))
      ENDDO
    ENDDO
!$OMP END DO
!$OMP DO PRIVATE(jb,jk) REDUCTION(max:w_aux_lev) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = istartblk_c, iendblk_c
      DO jk = 1, patch%nlev
        w_aux_lev(jk) = MAX(w_aux_lev(jk),w_aux(jb,jk))
      ENDDO
    ENDDO
!$OMP END DO
#endif

!$OMP END PARALLEL

    ! Add surface level for w
    jk = patch%nlevp1
    w_aux_lev (jk) = MAXVAL(w_aux(:,jk))

    !--- Get max over all PEs
    vmax(1)   = MAXVAL(vn_aux_lev(:))
    keyval(1) = MAXLOC(vn_aux_lev(:),1)
    vmax(2)   = MAXVAL(w_aux_lev(:))
    keyval(2) = MAXLOC(w_aux_lev(:),1)

    proc_id(:) = get_my_mpi_all_id()

    vmax       = global_max(vmax, proc_id=proc_id, keyval=keyval, iroot=process_mpi_stdio_id)

    max_vn         = vmax(1)
    max_vn_level   = keyval(1)
    max_vn_process = proc_id(1)

    max_w          = vmax(2)
    max_w_level    = keyval(2)
    max_w_process  = proc_id(2)

  END SUBROUTINE calculate_maxwinds


  !-------------------------------------------------------------------------
  !>
  !! Computes the maximum horizontal CFL number for the number of model levels (counted from the model top)
  !! selected by namelist. This is used for an adaptive adjustment of the physics-dynamics time step ratio
  !!
  SUBROUTINE compute_hcfl(p_patch, vn, dtime, nlev_hcfl, max_hcfl)

    TYPE(t_patch), INTENT(IN)  :: p_patch
    REAL(wp),      INTENT(IN)  :: vn(:,:,:), dtime
    INTEGER,       INTENT(IN)  :: nlev_hcfl
    REAL(wp),      INTENT(OUT) :: max_hcfl

    INTEGER  :: istartblk, iendblk, i_startidx, i_endidx
    INTEGER  :: jb, jk, je
    REAL(wp) :: max_hcfl_blk(p_patch%nblks_e), max_hcfl_tmp ! workaround for old NVIDIA compilers

    !-----------------------------------------------------------------------

    istartblk = p_patch%edges%start_block(grf_bdywidth_e+1)
    iendblk   = p_patch%edges%end_block(min_rledge_int)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx, max_hcfl_tmp) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = istartblk, iendblk

      CALL get_indices_e(p_patch, jb, istartblk, iendblk, i_startidx, i_endidx, &
                         grf_bdywidth_e+1, min_rledge_int)

      max_hcfl_tmp = 0._wp
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) REDUCTION(MAX: max_hcfl_tmp)
      DO jk = 1, nlev_hcfl
        DO je = i_startidx, i_endidx
          max_hcfl_tmp = MAX(max_hcfl_tmp,dtime*ABS(vn(je,jk,jb))*p_patch%edges%inv_dual_edge_length(je,jb))
        ENDDO
      ENDDO
      !$ACC END PARALLEL

      max_hcfl_blk(jb) = max_hcfl_tmp

    END DO
!$OMP END DO
!$OMP END PARALLEL


    max_hcfl = MAXVAL(max_hcfl_blk(istartblk:iendblk))

  END SUBROUTINE compute_hcfl

  !>
  !! Compute surface pressure time tendency abs(dpsdt)
  !!
  !! Compute surface pressure time tendency. If desired, 
  !! a spacial average is computed for the domain given 
  !! and written to the log file. 
  !! 
  SUBROUTINE compute_dpsdt (pt_patch, dt, pt_diag, lacc)

    TYPE(t_patch),       INTENT(IN)    :: pt_patch     !< grid/patch info
    REAL(wp),            INTENT(IN)    :: dt           !< time step [s]
    TYPE(t_nh_diag),     INTENT(INOUT) :: pt_diag      !< the diagnostic variables
    LOGICAL,  OPTIONAL,  INTENT(IN)    :: lacc         !< running on GPU if lacc=.TRUE.

    ! local
    INTEGER :: jc, jb                         !< loop indices
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx

    REAL(wp) :: dps_blk(pt_patch%nblks_c), dps_blk_scal
    INTEGER  :: npoints_blk(pt_patch%nblks_c), npoints_blk_scal
    REAL(wp) :: dpsdt_avg                     !< spatial average of ABS(dpsdt)
    INTEGER  :: npoints
    LOGICAL  :: lzacc
  !-------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_block(rl_start)
    i_endblk   = pt_patch%cells%end_block(rl_end)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP PRIVATE(jc)
      DO jc = i_startidx, i_endidx
        pt_diag%ddt_pres_sfc(jc,jb) = (pt_diag%pres_sfc(jc,jb)-pt_diag%pres_sfc_old(jc,jb))/dt
        pt_diag%pres_sfc_old(jc,jb) = pt_diag%pres_sfc(jc,jb)
      ENDDO
      !$ACC END PARALLEL
    ENDDO
!$OMP END DO

    IF (msg_level >= 11 .AND. .NOT. p_test_run) THEN
      ! dpsdt diagnostic - omitted in the case of a parallelization test (p_test_run) because this
      ! is a purely diagnostic quantity, for which it does not make sense to implement an order-invariant
      ! summation
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,dps_blk_scal,npoints_blk_scal)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        dps_blk_scal = 0._wp
        npoints_blk_scal = 0
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) REDUCTION(+: dps_blk_scal, npoints_blk_scal) IF(lzacc)
        DO jc = i_startidx, i_endidx
          dps_blk_scal = dps_blk_scal + ABS(pt_diag%ddt_pres_sfc(jc,jb))
          npoints_blk_scal = npoints_blk_scal + 1
        ENDDO
        !$ACC END PARALLEL LOOP
        !$ACC WAIT(1)
        dps_blk(jb) = dps_blk_scal
        npoints_blk(jb) = npoints_blk_scal
      ENDDO
!$OMP END DO

!$OMP MASTER
      dpsdt_avg = SUM(dps_blk(i_startblk:i_endblk))
      npoints   = SUM(npoints_blk(i_startblk:i_endblk))
      dpsdt_avg = global_sum_array(dpsdt_avg, opt_iroot=process_mpi_stdio_id)
      npoints   = global_sum_array(npoints  , opt_iroot=process_mpi_stdio_id)
      IF (my_process_is_stdio()) THEN
        dpsdt_avg = dpsdt_avg/(REAL(npoints,wp))
        ! Exclude initial time step where pres_sfc_old is zero
        IF (dpsdt_avg < 10000._wp/dt) THEN
          WRITE(message_text,'(a,f12.6,a,i3)') 'average |dPS/dt| =',dpsdt_avg,' Pa/s in domain',pt_patch%id
          CALL message('nwp_nh_interface: ', message_text)
       ENDIF
      ENDIF
!$OMP END MASTER
    ENDIF  ! msg_level
!$OMP END PARALLEL
    !$ACC WAIT(1)

  END SUBROUTINE compute_dpsdt

END MODULE mo_nh_supervise


