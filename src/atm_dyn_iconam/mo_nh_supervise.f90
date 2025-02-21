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

MODULE mo_nh_supervise

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_nonhydro_types,      ONLY: t_nh_state, t_nh_prog, t_nh_diag
  USE mo_model_domain,        ONLY: t_patch
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_cell
  USE mo_grid_config,         ONLY: n_dom, l_limited_area, grid_sphere_radius
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
  USE mo_fortran_tools,       ONLY: init, set_acc_host_or_device, assert_lacc_equals_i_am_accel_node, assert_acc_device_only

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname   = 'mo_nh_supervise'

  ! Needed by supervise_total_integrals_nh to keep data between steps
  REAL(wp), ALLOCATABLE, SAVE :: z_total_tracer_old(:)
  REAL(wp), ALLOCATABLE, SAVE :: z_total_tracer_0(:)

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
    TYPE(t_patch),            INTENT(IN) :: patch(n_dom)      ! Patch
    TYPE(t_nh_state), TARGET, INTENT(IN) :: nh_state(n_dom)   ! NH State
    TYPE(t_int_state),        INTENT(IN) :: int_state(n_dom)  ! Interpolation State
    INTEGER,                  INTENT(IN) :: ntimlev(n_dom)    ! time level
    INTEGER,                  INTENT(IN) :: ntimlev_rcf(n_dom)! rcf time level
    LOGICAL,                  INTENT(IN) :: l_last_step
    LOGICAL, OPTIONAL,        INTENT(IN) :: lacc ! If true, use openacc

    REAL(wp), SAVE :: z_total_mass_0    !< total air mass including vapor and condensate at first time step [kg]
    REAL(wp), SAVE :: z_total_drymass_0 !< total dry air mass at first time step [kg]
    REAL(wp), SAVE :: z_total_energy_0  !< total energy at first time step [kg]
    !
    REAL(wp) :: z_total_mass            !< total air mass, including vapor and condensate partial densities [kg]
    REAL(wp) :: z_total_drymass         !< total dry air mass (i.e. excluding vapor and condensate) [kg]
    REAL(wp) :: z_kin_energy            !< kinetic energy [Nm]
    REAL(wp) :: z_pot_energy            !< potential energy [Nm]
    REAL(wp) :: z_int_energy            !< internal energy [Nm]
    REAL(wp) :: z_total_energy          !< kinetic + potential + internal energy
    REAL(wp) :: z_mean_surfp            !< mean surface pressure [Pa]
    !
    REAL(wp) :: z_total_mass_re         !< changes in total mass compared to first step
    REAL(wp) :: z_total_drymass_re      !< changes in total dry mass compared to first step
    REAL(wp) :: z_total_energy_re       ! changes in total energy compared to first step
    !
    REAL(wp) :: z_kin_energy_re         !< percentage of kinetic energy out of total energy
    REAL(wp) :: z_int_energy_re         !< percentage of internal energy out of total energy
    REAL(wp) :: z_pot_energy_re         !< percentage of potential energy out of total energy
    REAL(wp) :: z_volume

    REAL (wp):: z_total_tracer(ntracer)       ! total tracer mass
    REAL (wp):: z_aux_tracer(nproma,patch(1)%nblks_c)
    REAL (wp):: z_rel_err_tracer_s1(ntracer) ! relative error of total tracer
    REAL (wp):: z_rel_err_tracer(ntracer)    ! relative error of total tracer

    INTEGER :: jg, jb, jk, jc, jt             ! loop indices
    INTEGER :: nlen, npromz_c, nblks_c
    INTEGER :: nlev                           ! number of full levels
    INTEGER :: ist                            ! status variable
    TYPE(t_nh_prog), POINTER :: prog       ! prog state
    TYPE(t_nh_prog), POINTER :: prog_rcf   ! prog_rcf state
    TYPE(t_nh_diag), POINTER :: diag       ! diag state
    REAL(wp) :: z_qsum(nproma,patch(1)%nlev,patch(1)%nblks_c)  ! total condensate including vapour
    REAL(wp) :: z_total_moist, z_elapsed_time
    REAL(wp) :: z_ekin(nproma,patch(1)%nlev,patch(1)%nblks_c)

    REAL(wp) :: max_vn, max_w
    INTEGER  :: max_vn_level, max_vn_process, max_w_level, max_w_process

    !-----------------------------------------------------------------------------

    CALL assert_lacc_equals_i_am_accel_node('mo_nh_stepping:supervise_total_integrals_nh', lacc)
    CALL assert_acc_device_only('mo_nh_stepping:supervise_total_integrals_nh', lacc)

    !$ACC DATA CREATE(z_ekin, z_qsum, z_aux_tracer)


    ! Hack [ha]:
    IF (.NOT. ALLOCATED (z_total_tracer_old)) THEN
      ALLOCATE (z_total_tracer_old(ntracer), STAT=ist)
      IF(ist/=SUCCESS)THEN
        CALL finish ('mo_nh_stepping:supervise_total_integrals_nh', &
          'allocation of z_total_tracer_old failed')
      ENDIF
      ALLOCATE (z_total_tracer_0(ntracer), STAT=ist)
      IF(ist/=SUCCESS)THEN
        CALL finish ('mo_nh_stepping:supervise_total_integrals_nh', &
          'allocation of z_total_tracer_0 failed')
      ENDIF
      z_total_tracer_old = 0.0_wp
      z_total_tracer_0   = 0.0_wp
    END IF

    !  write(0,*) 'k_step=', k_step
    ! Open the datafile
    IF (k_step == 1 .AND. my_process_is_stdio()) THEN
      CALL open_total_integral_files()
    ENDIF


    z_elapsed_time = dtime*REAL(k_step,wp)/3600.0_wp

    jg = 1 ! It does not make sense to double-count nested domains!

    prog     => nh_state(jg)%prog(ntimlev(jg))
    prog_rcf => nh_state(jg)%prog(ntimlev_rcf(jg))
    diag     => nh_state(jg)%diag

    nblks_c   = patch(jg)%nblks_c
    npromz_c  = patch(jg)%npromz_c

    ! number of vertical levels
    nlev = patch(jg)%nlev

    IF (iforcing <= 1) THEN ! u and v are not diagnosed regularly if physics is turned off
      CALL rbf_vec_interpol_cell(prog%vn,patch(jg),int_state(jg),diag%u,diag%v, opt_acc_async=.TRUE.)
    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,jt,nlen)
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

      ! compute total condensate INCLUDING vapour
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO jc = 1, nlen
          z_qsum(jc,jk,jb) = 0._wp
        ENDDO
      ENDDO
      IF ( ltransport .AND. lforcing .AND. iforcing /= iheldsuarez ) THEN
        !$ACC LOOP SEQ
        DO jt=1, iqm_max
        !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
          DO jk = 1, nlev
            DO jc = 1, nlen
              z_qsum(jc,jk,jb) = z_qsum(jc,jk,jb) + prog_rcf%tracer(jc,jk,jb,jt)
            ENDDO  !jc
          ENDDO  !jk
        ENDDO  !jt
      ENDIF
      !$ACC END PARALLEL
    ENDDO
!$OMP END DO
!$OMP END PARALLEL







    z_total_mass    = 0.0_wp
    z_total_drymass = 0.0_wp
    z_kin_energy    = 0.0_wp
    z_pot_energy    = 0.0_wp
    z_int_energy    = 0.0_wp
    z_mean_surfp    = 0.0_wp
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(z_volume) &
      !$ACC   REDUCTION(+, z_total_mass, z_kin_energy, z_int_energy, z_pot_energy, z_total_drymass)
      DO jk = 1, nlev
        DO jc = 1, nlen
          z_volume = patch(jg)%cells%area(jc,jb)*nh_state(jg)%metrics%ddqz_z_full(jc,jk,jb) &
            &       /patch(jg)%n_patch_cells_g*nh_state(jg)%metrics%deepatmo_vol_mc(jk)
          z_total_mass = z_total_mass + &
            &              prog%rho(jc,jk,jb)*z_volume
          z_kin_energy = z_kin_energy + &
            &              prog%rho(jc,jk,jb)*z_ekin(jc,jk,jb)*z_volume
          z_int_energy = z_int_energy + &
            &              cvd*prog%exner(jc,jk,jb)*prog%rho(jc,jk,jb)*prog%theta_v(jc,jk,jb)*z_volume
          z_pot_energy = z_pot_energy + &
            &              prog%rho(jc,jk,jb)*nh_state(jg)%metrics%geopot(jc,jk,jb)*z_volume
          z_total_drymass = z_total_drymass + &
            &              prog%rho(jc,jk,jb)*(1._wp - z_qsum(jc,jk,jb))*z_volume
        ENDDO
      ENDDO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) REDUCTION(+, z_mean_surfp)
      DO jc = 1, nlen
        z_mean_surfp = z_mean_surfp + diag%pres_sfc(jc,jb)*patch(jg)%cells%area(jc,jb) /  &
          &  (4._wp*grid_sphere_radius**2*pi)
      ENDDO
      !$ACC END PARALLEL
      !$ACC WAIT
    ENDDO
    z_total_energy = z_int_energy+z_kin_energy+z_pot_energy


    IF (k_step == 1) THEN
      z_total_mass_0    = z_total_mass
      z_total_drymass_0 = z_total_drymass
      z_total_energy_0  = z_total_energy
    ENDIF

    ! percentage of single energies out of the total energy
    z_kin_energy_re = z_kin_energy/z_total_energy*100.0_wp
    z_int_energy_re = z_int_energy/z_total_energy*100.0_wp
    z_pot_energy_re = z_pot_energy/z_total_energy*100.0_wp
    ! changes compared to first step
    z_total_mass_re   = z_total_mass   /z_total_mass_0    -1.0_wp
    z_total_drymass_re= z_total_drymass/z_total_drymass_0 -1.0_wp
    z_total_energy_re = z_total_energy /z_total_energy_0  -1.0_wp

    IF (my_process_is_stdio()) THEN
      if (n_file_ti >= 0) &
        WRITE(n_file_ti,'(i8,7e20.12)') &
        &   k_step, z_total_mass_re, z_total_drymass_re, z_total_energy_re, z_kin_energy_re, &
        &   z_int_energy_re, z_pot_energy_re, z_mean_surfp
      IF (l_last_step) THEN
        if (n_file_ti >= 0) &
          CLOSE(n_file_ti)
      ENDIF
    ENDIF

    IF (ltransport  .OR. ( iforcing == inwp ) .OR. ( iforcing == iaes )) THEN

      z_total_tracer(:) = 0.0_wp ! init must not used be used (this is no GPU variable)
      z_total_moist = 0.0_wp

      DO jt=1, ntracer

        CALL init(z_aux_tracer(:,:), opt_acc_async=.TRUE.) ! reinitialize for each jt

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
              z_volume = patch(jg)%cells%area(jc,jb)             &
                &    * nh_state(jg)%metrics%ddqz_z_full(jc,jk,jb) &
                &    * prog%rho(jc,jk,jb) * nh_state(jg)%metrics%deepatmo_vol_mc(jk)

              z_aux_tracer(jc,jb) = z_aux_tracer(jc,jb)    &
                &    + prog_rcf%tracer(jc,jk,jb,jt) * z_volume
            ENDDO
          ENDDO
          !$ACC END PARALLEL

        ENDDO

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jb = 1, nblks_c
          DO jc = 1, nproma
            IF(.NOT.patch(jg)%cells%decomp_info%owner_mask(jc,jb)) z_aux_tracer(jc,jb) = 0._wp
          ENDDO
        ENDDO
        !$ACC END PARALLEL
        z_total_tracer(jt) = global_sum_array(z_aux_tracer(:,:), lacc=.TRUE.)

      ENDDO  ! ntracer

      IF (lforcing) THEN
        ! when run in physics mode, compute total mass of all tracers
        ! (water substances)
        z_total_moist = (SUM(z_total_tracer(:)))
      ENDIF


      ! Save total tracer mass at first time step
      IF (k_step == 1) THEN
        z_total_tracer_0(:)    = z_total_tracer(:)
        z_rel_err_tracer_s1(:) = 0._wp
        z_rel_err_tracer(:)    = 0._wp

      ELSE

        ! compute relative mass error
        ! a) relative to time step 1 (z_rel_err_tracer_s1)
        ! b) relative to the previous time step n-1 (z_rel_err_tracer)
        DO jt=1, ntracer

          IF (z_total_tracer_old(jt) == 0._wp) THEN
            z_rel_err_tracer(jt) = 0._wp
          ELSE
            z_rel_err_tracer(jt) = (z_total_tracer(jt)/z_total_tracer_old(jt)) - 1._wp
          ENDIF
          IF (z_total_tracer_0(jt) == 0._wp) THEN
            z_rel_err_tracer_s1(jt) = 0._wp
          ELSE
            IF (lforcing) THEN
              z_rel_err_tracer_s1(jt) = (z_total_moist/z_total_tracer_0(jt)) - 1._wp
            ELSE
              z_rel_err_tracer_s1(jt) = (z_total_tracer(jt)/z_total_tracer_0(jt)) - 1._wp
            ENDIF
          ENDIF

        ENDDO

      ENDIF

      ! save total tracer mass for the next step
      z_total_tracer_old(:) = z_total_tracer(:)

      DO jt=1, ntracer
        if (n_file_tti > 0) &
          WRITE(n_file_tti,'(i21,f22.8,i21,3e40.16)')             &
          k_step, z_elapsed_time, jt, z_total_tracer(jt), &
          z_rel_err_tracer(jt), z_rel_err_tracer_s1(jt)
      ENDDO

      IF (l_last_step) THEN
        if (n_file_ti >= 0) &
          CLOSE(n_file_ti)
        if (n_file_tti > 0) &
          CLOSE(n_file_tti)

        DEALLOCATE(z_total_tracer_old, z_total_tracer_0, STAT=ist)
        IF(ist/=SUCCESS)THEN
          CALL finish ('mo_nh_stepping:supervise_total_integrals_nh', &
            'deallocation of z_total_tracer_old, z_total_tracer_0 failed')
        ENDIF
      ENDIF

    ENDIF    ! ltransport

    ! write additional check quantities (check_global_quantities)
    CALL calculate_maxwinds(patch(1), prog%vn, prog%w, &
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
    !$ACC END DATA

  END SUBROUTINE supervise_total_integrals_nh
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE open_total_integral_files( )

    INTEGER :: istat

    n_file_ti = find_next_free_unit(100,1000)
    ! write(0,*) "n_file_ti", n_file_ti
    OPEN(UNIT=n_file_ti,FILE='total_integrals.dat',ACTION="write", &
         FORM='FORMATTED',IOSTAT=istat)
    IF (istat/=SUCCESS) THEN
      CALL finish('supervise_total_integrals_nh','could not open total_integrals.dat')
    ENDIF
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

      !!$       ALLOCATE(z_total_tracer_old(ntracer), STAT=ist)
      !!$       IF(ist/=SUCCESS)THEN
      !!$          CALL finish ('mo_nh_stepping:supervise_total_integrals_nh', &
      !!$               'allocation of z_total_tracer_old failed')
      !!$       ENDIF
      !!$       ALLOCATE(z_total_tracer_0(ntracer), STAT=ist)
      !!$       IF(ist/=SUCCESS)THEN
      !!$          CALL finish ('mo_nh_stepping:supervise_total_integrals_nh', &
      !!$               'allocation of z_total_tracer_0 failed')
      !!$       ENDIF


      n_file_tti = find_next_free_unit(100,1000)
      OPEN(UNIT=n_file_tti,FILE='tracer_total_integrals.dat',ACTION="write", &
           FORM='FORMATTED',IOSTAT=istat)
      IF (istat/=SUCCESS) THEN
        CALL finish('supervise_total_integrals_nh','could not open tracer_total_integrals.dat')
      ENDIF
      WRITE (n_file_tti,'(A22,A22,A22,3A40)') &
        ' TIMESTEP            ,',&
        ' ELAPSED TIME    (hr),',&
        ' TRACER NR        (#),',&
        ' TOTAL TRACER   (kg),',&
        ' RELATIVE ERROR to step N-1(TRACER)',&
        ' RELATIVE ERROR to step 1 (TRACER)'
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

    INTEGER  :: i_nchdom, istartblk_c, istartblk_e, iendblk_c, iendblk_e, i_startidx, i_endidx
    INTEGER  :: jb, jk, jg
    INTEGER  :: proc_id(2), keyval(2)
    LOGICAL :: lzacc ! non-optional version of lacc

    !-----------------------------------------------------------------------
    CALL set_acc_host_or_device(lzacc, lacc)

    i_nchdom    = MAX(1,patch%n_childdom)
    istartblk_c = patch%cells%start_blk(grf_bdywidth_c+1,1)
    istartblk_e = patch%edges%start_blk(grf_bdywidth_e+1,1)
    iendblk_c   = patch%cells%end_blk(min_rlcell_int,i_nchdom)
    iendblk_e   = patch%edges%end_blk(min_rledge_int,i_nchdom)
    jg          = patch%id

    !$ACC DATA PRESENT(vn, w, patch) COPYOUT(vn_aux, w_aux) IF(lzacc)

    IF (jg > 1 .OR. l_limited_area) THEN
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      vn_aux(1:MIN(istartblk_e,iendblk_e),:) = 0._wp
      w_aux (1:MIN(istartblk_c,iendblk_c),:) = 0._wp
      !$ACC END KERNELS
    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = istartblk_e, iendblk_e

      CALL get_indices_e(patch, jb, istartblk_e, iendblk_e, i_startidx, i_endidx, &
                         grf_bdywidth_e+1, min_rledge_int)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG PRIVATE(vn_aux_tmp)
!$NEC novector
      DO jk = 1, patch%nlev
        vn_aux(jb,jk) = MAXVAL(ABS(vn(i_startidx:i_endidx,jk,jb)))
      ENDDO
      !$ACC END PARALLEL
    END DO
!$OMP END DO

!$OMP DO PRIVATE(jb, jk, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = istartblk_c, iendblk_c

      CALL get_indices_c(patch, jb, istartblk_c, iendblk_c, i_startidx, i_endidx, &
                         grf_bdywidth_c+1, min_rlcell_int)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG PRIVATE(w_aux_tmp)
!$NEC novector
      DO jk = 1, patch%nlevp1
        w_aux(jb,jk) = MAXVAL(ABS(w(i_startidx:i_endidx,jk,jb)))
      ENDDO
      !$ACC END PARALLEL
    END DO
!$OMP END DO

    !$ACC WAIT
    !$ACC END DATA

! At this point vn_aux and w_aux reside on the host.  
! Avoid doing MAXVAL with OpenACC -- this is not well supported!
!$OMP DO PRIVATE(jk) ICON_OMP_DEFAULT_SCHEDULE

    DO jk = 1, patch%nlev
      vn_aux_lev(jk) = MAXVAL(vn_aux(:,jk))
      w_aux_lev (jk) = MAXVAL(w_aux(:,jk))
    END DO

!$OMP END DO

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


