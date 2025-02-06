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

! This module contains the interface routines for calling the RTTOV library
! for generating synthetic satellite images

MODULE mo_rttov_interface
  !
  !
  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: SUCCESS
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_model_domain,        ONLY: t_patch, t_grid_cells, p_patch_local_parent, p_patch
  USE mo_intp_data_strc,      ONLY: t_int_state, p_int_state_local_parent
  USE mo_grf_intp_data_strc,  ONLY: t_gridref_state, p_grf_state_local_parent
  USE mo_grf_nudgintp,        ONLY: interpol_scal_nudging
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: msg_level, iqv, iqc, iqi, iqs, ltimer
  USE mo_grid_config,         ONLY: l_limited_area
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
  USE mo_nonhydro_state,      ONLY: p_nh_state
  USE mo_ext_data_types,      ONLY: t_external_data, t_external_atmos
  USE mo_nwp_lnd_types,       ONLY: t_lnd_state, t_lnd_prog, t_lnd_diag
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_impl_constants,      ONLY: min_rlcell_int, RTTOV_BT_CL, RTTOV_BT_CS, &
    &                               RTTOV_RAD_CL, RTTOV_RAD_CS, vname_len
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_impl_constants_grf,  ONLY: grf_bdyintp_start_c, grf_fbk_start_c
  USE mo_communication,       ONLY: exchange_data, exchange_data_mult
  USE mo_nh_vert_interp,      ONLY: prepare_lin_intp, prepare_extrap, lin_intp
  USE mo_nh_vert_interp_ipz,  ONLY: z_at_plevels
  USE mo_satad,               ONLY: qsat_rho
  USE mo_physical_constants,  ONLY: rd, vtmpc1, earth_radius, grav
  USE mo_cufunctions,         ONLY: foealfa
  USE mo_math_constants,      ONLY: rad2deg
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_synsat
  USE mo_synsat_config,       ONLY: sat_compute, num_sensors, jpnsat, mchans,        &
    &                               total_numchans, total_channels, instruments,     &
    &                               addclouds, nlev_rttov, num_images, iwc2effdiam,  &
    &                               iceshape, zenmax10, get_synsat_name
  USE mo_name_list_output_config, ONLY: is_variable_in_output
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_mpi,                 ONLY: p_pe, p_comm_work, p_io, num_work_procs, p_barrier, &
    &                               get_my_mpi_all_id, i_am_accel_node
  USE mo_fortran_tools,       ONLY: assert_acc_device_only
#ifdef __USE_RTTOV
  USE mo_rtifc,               ONLY: rtifc_set_opts, rtifc_init, rtifc_fill_input, &
                                    rtifc_direct, rtifc_errmsg,                   &
                                    NO_ERROR, default_gas_units, gas_unit_ppmvdry
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rttov_initialize
  PUBLIC :: rttov_finalize
  PUBLIC :: copy_rttov_ubc, rttov_driver


  REAL(wp), PARAMETER :: conv_qv = 1.60771704e6_wp ! conversion factor kg/kg -> ppmv for humidity

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_rttov_interface'

  !> parameter: output verbosity (for debugging purposes)
  INTEGER, PARAMETER :: dbg_level = 0

  !> LOGICAL switches: which of the available channels and synsat images are actually computed?
  LOGICAL, ALLOCATABLE :: lsynsat_product(:), lsynsat_chan(:,:)

  !> Number of channels that are actually computed
  INTEGER  :: numchans(jpnsat)

  !> RTTOV-indices of sensors
  INTEGER  :: rttov_indx(jpnsat)

  !> RTTOV-indices of channels
  INTEGER, ALLOCATABLE :: chan_rtidx(:,:)

  !> computed channel index (necessary if only a small number of channels is selected)
  INTEGER, ALLOCATABLE :: chan_idx(:,:) ! (1:numchans, isens)

CONTAINS

  !> Initialize RTTOV interface.
  !
  SUBROUTINE rttov_initialize()
    CHARACTER(LEN=*), PARAMETER    :: routine = modname//"::rttov_initialize"
    INTEGER                    :: n_chans(num_sensors)
    INTEGER                    :: channels(mchans, num_sensors)
    CHARACTER(LEN=vname_len) :: shortname
    CHARACTER(LEN=128)         :: longname
    INTEGER                    :: isens, ierrstat, k, isynsat
#ifdef __USE_RTTOV
    INTEGER                    :: istatus
#endif

    ! --- determine, which of the satellite images have been actually
    !     requested by the user:

    ALLOCATE(lsynsat_product(4*mchans), &
      &      lsynsat_chan(num_sensors, mchans), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    lsynsat_product(:) = .FALSE.
    lsynsat_chan(:,:)  = .FALSE.
    DO isens = 1, num_sensors
      IF (sat_compute(isens)%lcloud_tem) THEN
        DO k = 1, total_numchans(isens)
          ! determine name of the synsat_idx'th satellite image:
          CALL get_synsat_name(.FALSE., .TRUE., k, shortname, longname)
          IF (dbg_level > 1)  WRITE (0,*) TRIM(shortname), " (k=", k, ")"
          isynsat = (k-1)*4+RTTOV_BT_CL
          ! check if image is in any of our output namelists
          lsynsat_product(isynsat) =  &
            &    is_variable_in_output(var_name=TRIM(shortname))
          IF (lsynsat_product(isynsat))  lsynsat_chan(isens, k) = .TRUE.
        ENDDO
      ENDIF
      IF (sat_compute(isens)%lclear_tem) THEN
        DO k = 1, total_numchans(isens)
          CALL get_synsat_name(.FALSE., .FALSE., k, shortname, longname)
          IF (dbg_level > 1)  WRITE (0,*) TRIM(shortname), " (k=", k, ")"
          isynsat = (k-1)*4+RTTOV_BT_CS
          lsynsat_product(isynsat) =  &
            &    is_variable_in_output(var_name=TRIM(shortname))
          IF (lsynsat_product(isynsat))  lsynsat_chan(isens, k) = .TRUE.
        ENDDO
      ENDIF
      IF (sat_compute(isens)%lcloud_rad) THEN
        DO k = 1, total_numchans(isens)
          CALL get_synsat_name(.TRUE., .TRUE., k, shortname, longname)
          IF (dbg_level > 1)  WRITE (0,*) TRIM(shortname), " (k=", k, ")"
          isynsat = (k-1)*4+RTTOV_RAD_CL
          lsynsat_product(isynsat) =  &
            &    is_variable_in_output(var_name=TRIM(shortname))
          IF (lsynsat_product(isynsat))  lsynsat_chan(isens, k) = .TRUE.
        ENDDO
      ENDIF
      IF (sat_compute(isens)%lclear_rad) THEN
        DO k = 1, total_numchans(isens)
          CALL get_synsat_name(.TRUE., .FALSE., k, shortname, longname)
          IF (dbg_level > 1)  WRITE (0,*) TRIM(shortname), " (k=", k, ")"
          isynsat = (k-1)*4+RTTOV_RAD_CS
          lsynsat_product(isynsat) =  &
            &    is_variable_in_output(var_name=TRIM(shortname))
          IF (lsynsat_product(isynsat))  lsynsat_chan(isens, k) = .TRUE.
        ENDDO
      ENDIF
    END DO
    !$ACC ENTER DATA COPYIN(lsynsat_product)

    IF (dbg_level > 1) THEN
      WRITE (0,*) routine, ": lsynsat_product = ", lsynsat_product
    END IF

    ! --- set up a mapping: computed channel -> list of available channels

    ALLOCATE(chan_idx(mchans,num_sensors), chan_rtidx(mchans,num_sensors), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    chan_idx(:,:) = 0    
    numchans(:)   = 0
    DO isens = 1, num_sensors
      channels(:,isens) = 0
      n_chans(isens)    = 0
      DO k=1,mchans
        IF (lsynsat_chan(isens, k)) THEN
          n_chans(isens) = n_chans(isens) + 1
          channels(n_chans(isens),isens) = total_channels(k,isens)
          chan_idx(n_chans(isens),isens) = k
        END IF
      END DO
      numchans(isens) = n_chans(isens)
    END DO


    IF (dbg_level > 1) THEN
      WRITE (0,*) routine, ": n_chans  = ", n_chans
      WRITE (0,*) routine, ": chan_idx = ", chan_idx
    END IF

    ! --- initialize RTTOV

#ifdef __USE_RTTOV
    IF (ANY(n_chans(1:num_sensors) > 0)) THEN    
      IF (dbg_level > 2) THEN
        CALL p_barrier(p_comm_work)
        WRITE (0,*) routine, ": CALL to rttov_init"
      END IF

      DO isens = 1, num_sensors
        CALL rtifc_set_opts(                  &
          iopt           = rttov_indx(isens), &
          new            = .true.,            &
          addclouds      = .true.,            &
          apply_reg_lims = .true.,            &
          addsolar       = .false.,           &
          addrefrac      = .true.,            &
          addinterp      = .false.,           &
          use_q2m        = .true.,            &
          cloud_overlap  = 2,                 & !Simplified, much faster scheme
          do_checkinput  = .false.,           & !do_checkinput is expensive on vector machine
          init           = .true.)
      END DO

      CALL rtifc_init(  &
        instruments,    &
        channels,       &
        n_chans,        &
        chan_rtidx,     &
        rttov_indx,     &
        p_pe,           &
        num_work_procs, &
        p_io,           &
        p_comm_work,    &
        istatus      )

      IF (istatus /= NO_ERROR) THEN
        WRITE(message_text,'(a)') TRIM(rtifc_errmsg(istatus))
        CALL finish(routine ,message_text)
      ENDIF
      IF (dbg_level > 2) THEN
        CALL p_barrier(p_comm_work)
        WRITE (0,*) routine, ": CALL to rttov_init done."
      END IF
    END IF
#endif
  END SUBROUTINE rttov_initialize


  !> Deallocate some data structures for RTTOV
  !
  SUBROUTINE rttov_finalize()
    CHARACTER(LEN=*), PARAMETER    :: routine = modname//"::rttov_finalize"
    INTEGER :: ierrstat
    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(lsynsat_product)

    DEALLOCATE(lsynsat_product, lsynsat_chan, chan_idx, chan_rtidx, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END SUBROUTINE rttov_finalize


!>
!! Driver routine to call the RTTOV library for computing synthetic satellite images
!! For the time being, using a reduced radiation grid is mandatory
!!
SUBROUTINE rttov_driver (prm_diag, p_lnd_state, ext_data, jg, jgp, nnow, lacc)

  TYPE(t_nwp_phy_diag),  TARGET, INTENT(IN) :: prm_diag     ! for dom jg
  TYPE(t_lnd_state),     TARGET, INTENT(IN) :: p_lnd_state  ! for dom jg
  TYPE(t_external_data), TARGET, INTENT(IN) :: ext_data     ! for dom jg
  INTEGER, INTENT(IN) :: jg, jgp ! grid ID and parent grid ID
  INTEGER, INTENT(IN) :: nnow    ! time level nnow valid for long time step
  LOGICAL, OPTIONAL,   INTENT(IN)   :: lacc ! If true, use openacc

  TYPE(t_grid_cells),     POINTER :: p_gcp
  TYPE(t_patch),          POINTER :: p_pp
  TYPE(t_nh_prog),        POINTER :: p_nh_prog
  TYPE(t_nh_diag),        POINTER :: p_nh_diag
  TYPE(t_nh_metrics),     POINTER :: p_nh_metrics
  TYPE(t_external_atmos), POINTER :: p_extdata
  TYPE(t_lnd_prog),       POINTER :: p_lnd_prog
  TYPE(t_lnd_diag),       POINTER :: p_lnd_diag

  REAL(wp), DIMENSION(nproma, nlev_rttov, p_patch_local_parent(jg)%nblks_c) :: &
    temp_rttov, qv_rttov, qc_rttov, qcc_rttov, qi_rttov, qs_rttov, clc_rttov

  REAL(wp) :: rg_synsat(nproma, num_images, p_patch_local_parent(jg)%nblks_c)

  REAL(wp), DIMENSION(nproma, p_patch_local_parent(jg)%nblks_c) :: &
    rg_cosmu0, rg_psfc, rg_hsfc, rg_tsfc, rg_t2m, rg_qv2m, rg_u10m, rg_v10m

  INTEGER, DIMENSION(nproma, p_patch_local_parent(jg)%nblks_c) :: rg_stype, rg_wtype


  ! Local variables for RTTOV calls (with RTTOV-specific memory layout)
  REAL(wp), DIMENSION(nlev_rttov,nproma) :: temp, pres, qv

  REAL(wp) ::  cld(6,nlev_rttov-1,nproma), clc(nlev_rttov-1,nproma)

  INTEGER,  DIMENSION(1:nproma*MAXVAL(numchans(:))) :: iprof, ichan
  REAL(wp), DIMENSION(MAXVAL(numchans(:)),nproma) :: emiss, T_b, T_b_clear, rad, rad_clear

  REAL(wp) :: pres_rttov(nlev_rttov), r_sat, sat_a(nproma), sat_z(nproma), alpha_e, r_atm, lon

  INTEGER :: idg(nproma), ish(nproma)

  INTEGER :: jb, jc, jk, i_startblk, i_endblk, is, ie, j, k, rlstart
  INTEGER :: nlev_rg, isens, n_profs, ncalc, &
    &        istatus, isynsat

  ! Pointer to the acutally used variant of clc:
  REAL(wp), DIMENSION(:,:,:), POINTER ::  ptr_clc => NULL()


  CALL assert_acc_device_only('rttov_driver', lacc)

  ! first, check if nothing to do:
  IF (MAXVAL(numchans(:)) == 0)  RETURN

  IF (ltimer) CALL timer_start(timer_synsat)

  ! Skip empty patches which may occur with processor splitting or vector-host offloading
  ! This needs to be done after the timer call in order to avoid timer inconsistencies
  IF (p_patch(jg)%n_patch_cells == 0) THEN
    IF (ltimer) CALL timer_stop(timer_synsat)
    RETURN
  ENDIF

  nlev_rg = p_patch(jgp)%nlev

  p_gcp        => p_patch_local_parent(jg)%cells
  p_pp         => p_patch_local_parent(jg)
  p_nh_prog    => p_nh_state(jg)%prog(nnow)
  p_nh_diag    => p_nh_state(jg)%diag
  p_nh_metrics => p_nh_state(jg)%metrics
  p_extdata    => ext_data%atm
  p_lnd_prog   => p_lnd_state%prog_lnd(nnow)
  p_lnd_diag   => p_lnd_state%diag_lnd

  ! Decide which field for cloud cover has to be used:
  IF (atm_phy_nwp_config(jg)%luse_clc_rad) THEN
    ptr_clc => prm_diag%clc_rad
  ELSE
    ptr_clc => prm_diag%clc
  END IF

  ! distance of satellite from middle of the earth
  r_sat       = 35880.e3_wp + earth_radius

#ifdef __USE_RTTOV
  ! set backward compatibility mode for input unit of QV
  default_gas_units = gas_unit_ppmvdry
#endif

  ! Define RTTOV levels
  CALL define_rttov_levels (nlev_rttov, pres_rttov)

  !$ACC DATA &
  !$ACC   CREATE(temp_rttov, qv_rttov, qc_rttov, qcc_rttov, qi_rttov, qs_rttov, clc_rttov) &
  !$ACC   CREATE(rg_stype, rg_wtype) &
  !$ACC   CREATE(rg_cosmu0, rg_psfc, rg_hsfc, rg_tsfc, rg_t2m, rg_qv2m, rg_u10m, rg_v10m) &
  !$ACC   CREATE(temp, pres, qv) &
  !$ACC   CREATE(cld, clc) &
  !$ACC   COPYIN(pres_rttov)

  CALL prepare_rttov_input(jg, jgp, nlev_rg, p_nh_metrics%z_ifc, p_nh_diag%pres,               &
    p_nh_diag%dpres_mc, p_nh_diag%temp, prm_diag%tot_cld, ptr_clc,                             &
    p_nh_prog%tracer(:,:,:,iqs), prm_diag%con_udd(:,:,:,3), prm_diag%con_udd(:,:,:,7),         &
    p_extdata%fr_land, p_extdata%fr_lake, p_lnd_diag%fr_seaice, prm_diag%cosmu0,               &
    p_nh_diag%pres_sfc, p_lnd_prog%t_g, prm_diag%t_2m, prm_diag%qv_2m,                         &
    prm_diag%u_10m, prm_diag%v_10m, prm_diag%buffer_rttov,                                     &
    temp_rttov, qv_rttov, qc_rttov, qcc_rttov, qi_rttov, qs_rttov, clc_rttov,                  &
    rg_stype, rg_wtype, rg_cosmu0, rg_psfc, rg_hsfc, rg_tsfc, rg_t2m, rg_qv2m, rg_u10m, rg_v10m)

  ! Call RTTOV library
  !
  rlstart    = grf_bdyintp_start_c

  i_startblk = p_gcp%start_block(rlstart)
  i_endblk   = p_gcp%end_block(min_rlcell_int)

#ifdef __USE_RTTOV
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, is, ie, rlstart, min_rlcell_int)


    ! Choices for RTTOV ice cloud scheme and ice crystal shape
    idg(:) = iwc2effdiam ! McFarquhar et al. (2003); seems to the only numerically stable option
    ish(:) = iceshape    ! hexagonal crystals

    ! Copy input variables into RTTOV buffer
    IF (dbg_level > 2)  WRITE (0,*) "Copy input variables into RTTOV buffer"

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO jk = 1, nlev_rttov
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = is, ie
        pres(jk,jc) = pres_rttov(jk)
        temp(jk,jc) = temp_rttov(jc,jk,jb)
        qv(jk,jc)   = qv_rttov(jc,jk,jb)
      ENDDO
    ENDDO

    DO jk = 1, nlev_rttov-1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = is, ie
        ! cld(1) = stratiform cloud water, computed as total - convective
        cld(1,jk,jc) = MAX(0._wp, qc_rttov(jc,jk,jb) - qcc_rttov(jc,jk,jb))
        cld(2,jk,jc) = 0._wp
        ! cld(3) = convective cloud water
        cld(3,jk,jc) = qcc_rttov(jc,jk,jb)
        cld(4,jk,jc) = 0._wp
        cld(5,jk,jc) = 0._wp
        ! cld(6) = cloud ice
        cld(6,jk,jc) = qi_rttov(jc,jk,jb)
        ! clc = cloud fraction

        clc(jk,jc) = MIN(1._wp-1.e-8_wp,clc_rttov(jc,jk,jb))
        IF (ANY(cld(:,jk,jc) > 0._wp)) THEN
          clc(jk,jc) = MAX(1.e-8_wp,clc(jk,jc))
        ENDIF

      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC UPDATE ASYNC(1) &
    !$ACC   HOST(pres, temp, qv, rg_t2m, rg_qv2m, rg_psfc, rg_hsfc, rg_u10m) &
    !$ACC   HOST(rg_v10m, rg_tsfc, rg_stype, rg_wtype, rg_cosmu0, cld, clc)
    !$ACC WAIT(1)

    CALL rtifc_fill_input(                                &
         istatus,                                         &
         rttov_indx(1:num_sensors),                       &
         press      = pres(:,is:ie) ,                     &
         temp       = temp(:,is:ie),                      &
         humi       = qv(:,is:ie),                        &
         t2m        = rg_t2m(is:ie,jb),                   &
         q2m        = rg_qv2m(is:ie,jb),                  &
         psurf      = rg_psfc(is:ie,jb),                  &
         hsurf      = rg_hsfc(is:ie,jb),                  &
         u10m       = rg_u10m(is:ie,jb),                  &
         v10m       = rg_v10m(is:ie,jb),                  &
         stemp      = rg_tsfc(is:ie,jb),                  &
         stype      = rg_stype(is:ie,jb),                 &
         watertype  = rg_wtype(is:ie,jb),                 &
         lat        = p_gcp%center(is:ie,jb)%lat*rad2deg, &
         lon        = p_gcp%center(is:ie,jb)%lon*rad2deg, &
         sat_zen    = (/(0.0_wp, jc=is,ie)/),             &
         sun_zen    = rg_cosmu0(is:ie,jb),                & ! actually unused for addsolar=.false.
         cloud      = cld(:,:,is:ie),                     &
         cfrac      = clc(:,is:ie),                       &
         ice_scheme = ish(is:ie),                         &
         idg        = idg(is:ie),                         &
         ivect      = 1,                                  &
         pe         = p_pe)


    IF (istatus /= NO_ERROR) THEN
      WRITE(message_text,'(a)') TRIM(rtifc_errmsg(istatus))
      CALL finish('RTTOV fill_input ERROR' ,message_text)
    ENDIF

    sensor_loop: DO isens = 1, num_sensors

      ! Set up instrument, channel and profile numbers for RTTOV

      n_profs = ie-is+1
      ncalc = n_profs * numchans(isens)
      DO j = 1, numchans(isens)
        iprof(j:ncalc:numchans(isens)) = (/ (k, k=1,n_profs) /)
        ! note: ichan contains the list of channel indices wrt. the
        ! list provided to rttov_init:
        ichan(j:ncalc:numchans(isens)) =  (/ (chan_rtidx(j,isens), k=1,n_profs) /)
      ENDDO

      ! Set/compute some sensor dependent quantities

      DO jc = is, ie 
        ! Since the emissitivy is intent(inout) in RTTOV, we have to 
        ! reinitialize it
        emiss(:, jc-is+1) = 0.
        DO k = 1,  numchans(isens)
          emiss(k, jc-is+1) = sat_compute(isens)%emissivity(chan_idx(k,isens))
        ENDDO
        lon = p_gcp%center(jc,jb)%lat - sat_compute(isens)%longitude
        ! Calculate the satellite zenith angle 
        alpha_e  = ACOS(COS(p_gcp%center(jc,jb)%lat) * COS(lon))
        r_atm    = SQRT(r_sat**2 + earth_radius**2 -2*r_sat*earth_radius*COS(alpha_e))
        sat_z(jc) = ASIN(SIN(alpha_e)*r_sat/r_atm) * rad2deg
        sat_z(jc) = MIN(ABS(sat_z(jc)), zenmax10)

        ! Calculate the satellite azimuth angle
        sat_a(jc) = rad2deg * (1. + ATAN2(TAN(lon), SIN(p_gcp%center(jc,jb)%lat)))
      ENDDO


      IF (dbg_level > 2) THEN
        WRITE (0,*) "PE ", get_my_mpi_all_id(), " :: CALL to rttov_direct_ifc: isens = ", isens, &
          &         "; iprof = ", iprof(1:ncalc),                                                &
          &         ", ichan=", ichan(1:ncalc), ", emiss=", emiss(1:numchans(isens), 1:n_profs), &
          &         ", numchans=", numchans(isens)
      END IF
      CALL rtifc_direct(                                          &
             rttov_indx(isens),                                   &
             iprof(1:ncalc),                                      &
             ichan(1:ncalc),                                      &
             emiss(1:numchans(isens), 1:n_profs),                 &
             T_b(1:numchans(isens), 1:n_profs),                   &
             istatus,                                             &
             sat_azi   = sat_a(is:ie),                            &
             sat_zen   = sat_z(is:ie),                            &
             T_b_clear = T_b_clear(1:numchans(isens), 1:n_profs), &
             rad       = rad      (1:numchans(isens), 1:n_profs), &
             radclear  = rad_clear(1:numchans(isens), 1:n_profs), &
             pe        = p_pe)
      IF (istatus /= NO_ERROR) THEN
        WRITE(message_text,'(a)') TRIM(rtifc_errmsg(istatus))
        CALL finish('RTTOV synsat calc ERROR' ,message_text)
      ENDIF

      IF (dbg_level > 2)  WRITE (0,*) "PE ", get_my_mpi_all_id(), " :: copy result into rg_synsat array"
      IF (sat_compute(isens)%lcloud_tem) THEN
        DO k = 1, numchans(isens)
          isynsat = (chan_idx(k,isens)-1)*4+RTTOV_BT_CL
          IF (lsynsat_product(isynsat)) THEN
            IF (dbg_level > 2)  WRITE (0,*) "copy synsat product ", isynsat, " into place."
            rg_synsat(is:ie, isynsat, jb) = T_b(k, 1:n_profs) 
          ELSE
            rg_synsat(is:ie, isynsat, jb) = 0._wp
          END IF
        ENDDO
      ENDIF
      IF (sat_compute(isens)%lclear_tem) THEN
        DO k = 1, numchans(isens)
          isynsat = (chan_idx(k,isens)-1)*4+RTTOV_BT_CS
          IF (lsynsat_product(isynsat)) THEN
            IF (dbg_level > 2)  WRITE (0,*) "copy synsat product ", isynsat, " into place."
            rg_synsat(is:ie, isynsat, jb) = T_b_clear(k, 1:n_profs) 
          ELSE
            rg_synsat(is:ie, isynsat, jb) = 0._wp
          END IF
        ENDDO
      ENDIF
      IF (sat_compute(isens)%lcloud_rad) THEN
        DO k = 1, numchans(isens)
          isynsat = (chan_idx(k,isens)-1)*4+RTTOV_RAD_CL
          IF (lsynsat_product(isynsat)) THEN
            IF (dbg_level > 2)  WRITE (0,*) "copy synsat product ", isynsat, " into place."
            rg_synsat(is:ie, isynsat, jb) = Rad(k, 1:n_profs) 
          ELSE
            rg_synsat(is:ie, isynsat, jb) = 0._wp
          END IF
        ENDDO
      ENDIF
      IF (sat_compute(isens)%lclear_rad) THEN
        DO k = 1, numchans(isens)
          isynsat = (chan_idx(k,isens)-1)*4+RTTOV_RAD_CS
          IF (lsynsat_product(isynsat)) THEN
            IF (dbg_level > 2)  WRITE (0,*) "copy synsat product ", isynsat, " into place."
            rg_synsat(is:ie, isynsat, jb) = Rad_clear(k, 1:n_profs) 
          ELSE
            rg_synsat(is:ie, isynsat, jb) = 0._wp
          END IF
        ENDDO
      ENDIF

    END DO sensor_loop

  ENDDO
#endif
  !$ACC END DATA

  IF (dbg_level > 2) THEN
    CALL p_barrier(p_comm_work)
    WRITE (0,*) "CALL to downscale_rttov_output"
  END IF

  !$ACC DATA COPYIN(rg_synsat)
  CALL downscale_rttov_output(jg, num_images, rg_synsat, &
  &                           prm_diag%synsat_arr, lsynsat_product)
  !$ACC WAIT
  !$ACC END DATA

  IF (dbg_level > 2) THEN
    CALL p_barrier(p_comm_work)
    WRITE (0,*) "Leave downscale_rttov_output"
  END IF

  IF (ltimer) CALL timer_stop(timer_synsat)
END SUBROUTINE rttov_driver


!>
!! This routine interpolates the input fields for RTTOV to the reduced radiation grid,
!! combined with a vertical interpolation to the predefined RTTOV model levels
!!
!!
SUBROUTINE prepare_rttov_input(jg, jgp, nlev_rg, z_ifc, pres, dpres, temp, tot_cld, clc, qs, &
  det_rate, con_upd, fr_land, fr_lake, fr_seaice, cosmu0, psfc, tsfc, t2m, qv2m, u10m, v10m, &
  buffer_rttov, temp_rttov, qv_rttov, qc_rttov, qcc_rttov, qi_rttov, qs_rttov, clc_rttov,    &
  rg_stype, rg_wtype, rg_cosmu0, rg_psfc, rg_hsfc, rg_tsfc, rg_t2m, rg_qv2m, rg_u10m, rg_v10m)


  ! Input grid parameters
  INTEGER, INTENT(IN)  :: jg, jgp  ! domain IDs of main and reduced grids
  INTEGER, INTENT(IN)  :: nlev_rg  ! number of model levels on reduced grid

  ! Input fields (on full grid)
  REAL(wp), INTENT(IN) ::                                                                  &
    z_ifc(:,:,:), pres(:,:,:), dpres(:,:,:), temp(:,:,:), tot_cld(:,:,:,:), clc(:,:,:),    &
    qs(:,:,:), det_rate(:,:,:), con_upd(:,:,:), fr_land(:,:), fr_lake(:,:), fr_seaice(:,:),&
    cosmu0(:,:), tsfc(:,:), psfc(:,:), t2m(:,:), qv2m(:,:), u10m(:,:), v10m(:,:)

  ! Buffer for auxiliary temperature and pressure levels above the vertical nest interface
  REAL(wp), POINTER, INTENT(IN) :: buffer_rttov(:,:,:) ! pointer keyword is necessary, as it can be a nullified pointer.

  ! Corresponding output fields (on reduced grid)
  REAL(wp), INTENT(OUT) ::                                                                  &
    temp_rttov(:,:,:), qv_rttov(:,:,:), qc_rttov(:,:,:), qcc_rttov(:,:,:), qi_rttov(:,:,:), &
    qs_rttov(:,:,:), clc_rttov(:,:,:), rg_cosmu0(:,:), rg_tsfc(:,:), rg_psfc(:,:),          &
    rg_hsfc(:,:), rg_t2m(:,:), rg_qv2m(:,:), rg_u10m(:,:), rg_v10m(:,:)

  INTEGER, INTENT(OUT) :: rg_stype(:,:), rg_wtype(:,:)

  ! Local auxiliary fields
  REAL(wp), DIMENSION(nproma,nlev_rg,p_patch_local_parent(jg)%nblks_c)  ::   &
    rg_pres, rg_temp, rg_qv, rg_z_mc

  REAL(wp), DIMENSION(nproma,nlev_rg+1,p_patch_local_parent(jg)%nblks_c) ::  &
    rg_qc_vi, rg_qcc_vi, rg_qi_vi, rg_qs_vi, rg_clc_vi, rg_z_ifc

  REAL(wp), DIMENSION(nproma,nlev_rg) :: zdpres, zclc, zqc, zqcc, zqi, zqs

  REAL(wp), DIMENSION(nproma,nlev_rttov,p_patch_local_parent(jg)%nblks_c) ::  &
    wfac_lin, wfac_lin_i, z3d_rttov, pres3d_rttov, qc_vi_rttov, qcc_vi_rttov, &
    qi_vi_rttov, qs_vi_rttov, clc_vi_rttov

  REAL(wp), DIMENSION(nproma,p_patch_local_parent(jg)%nblks_c) :: wfacpbl1, wfacpbl2

  INTEGER,  DIMENSION(nproma,nlev_rttov,p_patch_local_parent(jg)%nblks_c) :: &
    idx0_lin, idx0_lin_i

  INTEGER,  DIMENSION(nproma,p_patch_local_parent(jg)%nblks_c) :: &
    kpbl1, kpbl2, bot_idx_lin, bot_idx_lin_i


  ! Convenience pointers
  TYPE(t_grid_cells), POINTER     :: p_gcp
  TYPE(t_patch),      POINTER     :: p_pp

  INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk
  REAL(wp), POINTER :: p_fbkwgt(:,:,:)

  ! Indices
  INTEGER :: jb, jc, jk, jk1, &
             i_startblk, i_endblk, i_startidx, i_endidx, nblks_c, npromz_c

  INTEGER :: nlev, nlevp1      !< number of full and half levels
  INTEGER :: nshift, nlevp1_rg

  REAL(wp) :: qv_aux, pres_rttov(nlev_rttov), zfr_land, zfr_lake, zfr_si, zrho, &
              zdet_rate, zcon_upd

  !-----------------------------------------------------------------------

  IF (msg_level >= 10) THEN
    WRITE(message_text,'(a,i2,a,i2)') 'vertical interpolation and upscaling of RTTOV input',&
      jg,' =>',jgp
    CALL message('prepare_rttov_input',message_text)
  ENDIF

  ! Number of levels of the full grid
  nlev   = p_patch(jg)%nlev
  nlevp1 = p_patch(jg)%nlevp1

  ! nlev_rg carries the number of model levels of the reduced grid,
  ! which may be larger than nlev
  nlevp1_rg = nlev_rg + 1
  nshift = nlev_rg - nlev ! resulting shift parameter


  p_gcp => p_patch_local_parent(jg)%cells
  p_pp  => p_patch_local_parent(jg)

  ! Set pointers to index and coefficient fields for cell-based variables
  iidx => p_gcp%child_idx
  iblk => p_gcp%child_blk

  p_fbkwgt => p_grf_state_local_parent(jg)%fbk_wgt_bln


  ! Define RTTOV levels
  CALL define_rttov_levels (nlev_rttov, pres_rttov)

  !$ACC DATA &
  !$ACC   CREATE(rg_pres, rg_temp, rg_qv, rg_z_mc) &
  !$ACC   CREATE(rg_qc_vi, rg_qcc_vi, rg_qi_vi, rg_qs_vi, rg_clc_vi, rg_z_ifc) &
  !$ACC   CREATE(zdpres, zclc, zqc, zqcc, zqi, zqs) &
  !$ACC   CREATE(wfac_lin, wfac_lin_i, z3d_rttov, pres3d_rttov, qc_vi_rttov, qcc_vi_rttov) &
  !$ACC   CREATE(qi_vi_rttov, qs_vi_rttov, clc_vi_rttov) &
  !$ACC   CREATE(wfacpbl1, wfacpbl2) &
  !$ACC   CREATE(idx0_lin, idx0_lin_i) &
  !$ACC   CREATE(kpbl1, kpbl2, bot_idx_lin, bot_idx_lin_i) &
  !$ACC   COPYIN(pres_rttov)
!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

  ! Upscaling of RTTOV input variables from compute grid to reduced grid
  i_startblk = p_gcp%start_block(grf_bdyintp_start_c)
  i_endblk   = p_gcp%end_block(min_rlcell_int)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jk1,zdpres,zrho,qv_aux,zclc,zqc,zqcc,zqi,zqs, &
!$OMP            zfr_land,zfr_lake,zfr_si,zdet_rate,zcon_upd)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk,          &
      i_startidx, i_endidx, grf_bdyintp_start_c, min_rlcell_int)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO jk = 1, nlev_rttov
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = i_startidx, i_endidx
        pres3d_rttov(jc,jk,jb) = pres_rttov(jk)*100._wp  ! convert in Pa
      ENDDO
    ENDDO

    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jc = i_startidx, i_endidx
      rg_clc_vi(jc,nlevp1_rg,jb) = 0._wp
      rg_qc_vi(jc,nlevp1_rg,jb)  = 0._wp
      rg_qcc_vi(jc,nlevp1_rg,jb) = 0._wp
      rg_qi_vi(jc,nlevp1_rg,jb)  = 0._wp
      rg_qs_vi(jc,nlevp1_rg,jb)  = 0._wp

      rg_z_ifc(jc,nlevp1_rg,jb) =                                     &
        z_ifc(iidx(jc,jb,1),nlevp1,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        z_ifc(iidx(jc,jb,2),nlevp1,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        z_ifc(iidx(jc,jb,3),nlevp1,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        z_ifc(iidx(jc,jb,4),nlevp1,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)
    ENDDO

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
!DIR$ IVDEP
      DO jk = 1, nlev
        jk1 = jk + nshift
#else
    !$ACC LOOP SEQ
    DO jk = 1, nlev
      jk1 = jk + nshift
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(qv_aux, zrho, zdet_rate, zcon_upd)
      DO jc = i_startidx, i_endidx
#endif

        rg_z_ifc(jc,jk1,jb) =                                       &
          z_ifc(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          z_ifc(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          z_ifc(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          z_ifc(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        rg_pres(jc,jk1,jb) =                                       &
          pres(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          pres(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          pres(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          pres(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        zdpres(jc,jk1) =                                            &
          dpres(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          dpres(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          dpres(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          dpres(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        rg_temp(jc,jk1,jb) =                                       &
          temp(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          temp(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          temp(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          temp(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        qv_aux =                                                          &
          tot_cld(iidx(jc,jb,1),jk,iblk(jc,jb,1),iqv)*p_fbkwgt(jc,jb,1) + &
          tot_cld(iidx(jc,jb,2),jk,iblk(jc,jb,2),iqv)*p_fbkwgt(jc,jb,2) + &
          tot_cld(iidx(jc,jb,3),jk,iblk(jc,jb,3),iqv)*p_fbkwgt(jc,jb,3) + &
          tot_cld(iidx(jc,jb,4),jk,iblk(jc,jb,4),iqv)*p_fbkwgt(jc,jb,4)

        ! limit qv to saturation and convert into volume mixing ratio (ppmv)
        zrho = rg_pres(jc,jk1,jb)/(rd*(rg_temp(jc,jk1,jb)+vtmpc1*qv_aux))
        qv_aux = MIN(qv_aux, qsat_rho(rg_temp(jc,jk1,jb),zrho))
        rg_qv(jc,jk1,jb) = qv_aux / (1._wp - qv_aux) * conv_qv

        zclc(jc,jk1) =                                            &
          clc(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          clc(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          clc(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          clc(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        zqc(jc,jk1) =                                                     &
          tot_cld(iidx(jc,jb,1),jk,iblk(jc,jb,1),iqc)*p_fbkwgt(jc,jb,1) + &
          tot_cld(iidx(jc,jb,2),jk,iblk(jc,jb,2),iqc)*p_fbkwgt(jc,jb,2) + &
          tot_cld(iidx(jc,jb,3),jk,iblk(jc,jb,3),iqc)*p_fbkwgt(jc,jb,3) + &
          tot_cld(iidx(jc,jb,4),jk,iblk(jc,jb,4),iqc)*p_fbkwgt(jc,jb,4)

        ! convective cloud water (computation taken over from cloud cover scheme)
        !
        ! averaged convective updraft detrainment rate
        zdet_rate =                                                    &
          det_rate(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          det_rate(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          det_rate(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          det_rate(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        ! averaged updraft condensate
        zcon_upd =                                                    &
          con_upd(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          con_upd(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          con_upd(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          con_upd(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        ! diagnosed averaged convective cloud water, must not exceed total QC
        zqcc(jc,jk1) = MIN((zdet_rate/zrho) / (zdet_rate/zrho + 1.0_wp/1500._wp) * &
                       zcon_upd * foealfa(rg_temp(jc,jk1,jb)), zqc(jc,jk1) )

        zqi(jc,jk1) =                                                     &
          tot_cld(iidx(jc,jb,1),jk,iblk(jc,jb,1),iqi)*p_fbkwgt(jc,jb,1) + &
          tot_cld(iidx(jc,jb,2),jk,iblk(jc,jb,2),iqi)*p_fbkwgt(jc,jb,2) + &
          tot_cld(iidx(jc,jb,3),jk,iblk(jc,jb,3),iqi)*p_fbkwgt(jc,jb,3) + &
          tot_cld(iidx(jc,jb,4),jk,iblk(jc,jb,4),iqi)*p_fbkwgt(jc,jb,4)

        zqs(jc,jk1) =                                            &
          qs(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          qs(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          qs(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          qs(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      ENDDO
    ENDDO


    ! Fill levels above model top
    !$ACC LOOP SEQ
    DO jk = 1, nshift
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(qv_aux)
      DO jc = i_startidx, i_endidx
        rg_z_ifc(jc,jk,jb) = buffer_rttov(jc,jk,jb)
        rg_pres (jc,jk,jb) = buffer_rttov(jc,nshift+jk,jb)
        zdpres(jc,jk)      = buffer_rttov(jc,2*nshift+jk,jb)
        rg_temp (jc,jk,jb) = buffer_rttov(jc,3*nshift+jk,jb)
        qv_aux             = buffer_rttov(jc,4*nshift+jk,jb)
        !
        ! convert qv into volume mixing ratio (ppmv)
        rg_qv(jc,jk,jb)    = qv_aux / (1._wp - qv_aux) * conv_qv
        !
        ! for the cloud variables, we assume that the top of the nested model domain is not
        ! lower than the upper boundary of moisture physics calculation
        zclc(jc,jk) = 0._wp
        zqc(jc,jk)  = 0._wp
        zqcc(jc,jk) = 0._wp
        zqi(jc,jk)  = 0._wp
        zqs(jc,jk)  = 0._wp
      ENDDO
    ENDDO

    ! Compute rg_z_mc
    !$ACC LOOP SEQ
    DO jk = 1, nlev_rg
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = i_startidx, i_endidx
        rg_z_mc(jc,jk,jb) = 0.5_wp*(rg_z_ifc(jc,jk,jb)+rg_z_ifc(jc,jk+1,jb))
      ENDDO
    ENDDO

    ! Vertically integrate the cloud quantities
    !$ACC LOOP SEQ
    DO jk = nlev_rg, 1, -1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = i_startidx, i_endidx
        rg_clc_vi(jc,jk,jb) = rg_clc_vi(jc,jk+1,jb) + zdpres(jc,jk)*zclc(jc,jk)
        rg_qc_vi(jc,jk,jb)  = rg_qc_vi(jc,jk+1,jb)  + zdpres(jc,jk)*zqc(jc,jk)
        rg_qcc_vi(jc,jk,jb) = rg_qcc_vi(jc,jk+1,jb) + zdpres(jc,jk)*zqcc(jc,jk)
        rg_qi_vi(jc,jk,jb)  = rg_qi_vi(jc,jk+1,jb)  + zdpres(jc,jk)*zqi(jc,jk)
        rg_qs_vi(jc,jk,jb)  = rg_qs_vi(jc,jk+1,jb)  + zdpres(jc,jk)*zqs(jc,jk)
      ENDDO
    ENDDO

    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zfr_land, zfr_lake, zfr_si, qv_aux)
    DO jc = i_startidx, i_endidx

      rg_hsfc(jc,jb) = rg_z_ifc(jc,nlevp1_rg,jb) * 0.001_wp ! convert in km

      zfr_land =                                                 &
        fr_land(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        fr_land(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        fr_land(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        fr_land(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      zfr_lake =                                                 &
        fr_lake(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        fr_lake(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        fr_lake(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        fr_lake(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      zfr_si =                                                     &
        fr_seaice(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        fr_seaice(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        fr_seaice(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        fr_seaice(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      IF (zfr_land > 0.5_wp) THEN
        rg_stype(jc,jb) = 0
        rg_wtype(jc,jb) = 0
      ELSE IF (zfr_si > 0.5_wp) THEN
        rg_stype(jc,jb) = 2
        rg_wtype(jc,jb) = 1
      ELSE 
        rg_stype(jc,jb) = 1
        IF (zfr_lake > 0.5_wp) THEN
          rg_wtype(jc,jb) = 0
        ELSE
          rg_wtype(jc,jb) = 1
        ENDIF
      ENDIF
      
      rg_cosmu0(jc,jb) =                                        &
        cosmu0(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        cosmu0(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        cosmu0(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        cosmu0(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      rg_tsfc(jc,jb) =                                        &
        tsfc(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        tsfc(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        tsfc(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        tsfc(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      rg_psfc(jc,jb) =  MIN(1099.9_wp, 0.01_wp* (             & ! convert in hPa
        psfc(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        psfc(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        psfc(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        psfc(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4) ) )

      rg_t2m(jc,jb) =                                        &
        t2m(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        t2m(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        t2m(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        t2m(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      qv_aux =                                                &
        qv2m(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        qv2m(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        qv2m(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        qv2m(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      ! convert qv into volume mixing ratio (ppmv)
      rg_qv2m(jc,jb) = qv_aux / (1._wp - qv_aux) * conv_qv

      rg_u10m(jc,jb) =                                        &
        u10m(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        u10m(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        u10m(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        u10m(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      rg_v10m(jc,jb) =                                        &
        v10m(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        v10m(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        v10m(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        v10m(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

    ENDDO
    !$ACC END PARALLEL

  ENDDO
!$OMP END DO

  ! Fill surface values in the boundary interpolation zone
  IF (l_limited_area .AND. jg == 1) THEN

    i_startblk = p_gcp%start_block(grf_bdyintp_start_c)
    i_endblk   = p_gcp%end_block(grf_fbk_start_c+1)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_pp, jb, i_startblk, i_endblk,          &
        i_startidx, i_endidx, grf_bdyintp_start_c, grf_fbk_start_c+1)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        rg_tsfc(jc,jb) = rg_temp(jc,nlev_rg,jb)
        rg_t2m(jc,jb)  = rg_temp(jc,nlev_rg,jb)
        rg_qv2m(jc,jb) = rg_qv  (jc,nlev_rg,jb)
        rg_u10m(jc,jb) = 0._wp
        rg_v10m(jc,jb) = 0._wp
      ENDDO
      !$ACC END PARALLEL

    ENDDO
!$OMP END DO
  ENDIF
!$OMP END PARALLEL


  ! Compute interpolation coefficients

  nblks_c  = p_gcp%end_block(min_rlcell_int)
  npromz_c = p_gcp%end_index(min_rlcell_int)

  CALL prepare_extrap(rg_z_mc, nblks_c, npromz_c, nlev_rg,    &
                      kpbl1, wfacpbl1, kpbl2, wfacpbl2, lacc=.TRUE.)

  CALL z_at_plevels(rg_pres, rg_temp, rg_z_mc, pres3d_rttov, z3d_rttov, &
                    nblks_c, npromz_c, nlev_rg, nlev_rttov,             &
                    kpbl1, wfacpbl1, kpbl2, wfacpbl2, lacc=.TRUE.      )

  CALL prepare_lin_intp(rg_z_mc, z3d_rttov, nblks_c, npromz_c, nlev_rg,             &
                        nlev_rttov, wfac_lin, idx0_lin, bot_idx_lin, lacc=.TRUE., lextrap=.FALSE.)

  CALL prepare_lin_intp(rg_z_ifc, z3d_rttov, nblks_c, npromz_c, nlevp1_rg,                &
                        nlev_rttov, wfac_lin_i, idx0_lin_i, bot_idx_lin_i, lacc=.TRUE., lextrap=.FALSE.)

  ! Execute interpolation

  CALL lin_intp(rg_temp, temp_rttov,  nblks_c, npromz_c, nlev_rg,   &
                nlev_rttov, wfac_lin, idx0_lin, bot_idx_lin,        &
                wfacpbl1, kpbl1, wfacpbl2, kpbl2, l_loglin=.FALSE., &
                l_pd_limit=.TRUE., l_extrapol=.FALSE., lacc=.TRUE. )

  CALL lin_intp(rg_qv, qv_rttov,  nblks_c, npromz_c, nlev_rg,       &
                nlev_rttov, wfac_lin, idx0_lin, bot_idx_lin,        &
                wfacpbl1, kpbl1, wfacpbl2, kpbl2, l_loglin=.TRUE.,  &
                l_pd_limit=.TRUE., l_extrapol=.FALSE., lacc=.TRUE. )

  CALL lin_intp(rg_qc_vi, qc_vi_rttov, nblks_c, npromz_c, nlevp1_rg, &
                nlev_rttov, wfac_lin, idx0_lin, bot_idx_lin,         &
                wfacpbl1, kpbl1, wfacpbl2, kpbl2, l_loglin=.TRUE.,   &
                l_pd_limit=.TRUE., l_extrapol=.FALSE., lacc=.TRUE.  )

  CALL lin_intp(rg_qcc_vi, qcc_vi_rttov, nblks_c, npromz_c, nlevp1_rg, &
                nlev_rttov, wfac_lin, idx0_lin, bot_idx_lin,           &
                wfacpbl1, kpbl1, wfacpbl2, kpbl2, l_loglin=.TRUE.,     &
                l_pd_limit=.TRUE., l_extrapol=.FALSE., lacc=.TRUE.    )

  CALL lin_intp(rg_qi_vi, qi_vi_rttov, nblks_c, npromz_c, nlevp1_rg, &
                nlev_rttov, wfac_lin, idx0_lin, bot_idx_lin,         &
                wfacpbl1, kpbl1, wfacpbl2, kpbl2, l_loglin=.TRUE.,   &
                l_pd_limit=.TRUE., l_extrapol=.FALSE., lacc=.TRUE.  )

  CALL lin_intp(rg_qs_vi, qs_vi_rttov, nblks_c, npromz_c, nlevp1_rg, &
                nlev_rttov, wfac_lin, idx0_lin, bot_idx_lin,         &
                wfacpbl1, kpbl1, wfacpbl2, kpbl2, l_loglin=.TRUE.,   &
                l_pd_limit=.TRUE., l_extrapol=.FALSE., lacc=.TRUE.  )

  CALL lin_intp(rg_clc_vi, clc_vi_rttov, nblks_c, npromz_c, nlevp1_rg, &
                nlev_rttov, wfac_lin, idx0_lin, bot_idx_lin,           &
                wfacpbl1, kpbl1, wfacpbl2, kpbl2, l_loglin=.TRUE.,     &
                l_pd_limit=.TRUE., l_extrapol=.FALSE., lacc=.TRUE.    )


  ! Calculate layer averages of cloud variables from vertically integrated fields
  ! The unit required by RTTOV is g/m**3 (!!)

  i_startblk = p_gcp%start_block(grf_bdyintp_start_c)
  i_endblk   = p_gcp%end_block(min_rlcell_int)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk,          &
      i_startidx, i_endidx, grf_bdyintp_start_c, min_rlcell_int)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO jk = 1, nlev_rttov-1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = i_startidx, i_endidx
        qc_rttov(jc,jk,jb)  = 1.e3_wp*(qc_vi_rttov(jc,jk,jb)-qc_vi_rttov(jc,jk+1,jb)) / &
                              (grav*(z3d_rttov(jc,jk,jb)-z3d_rttov(jc,jk+1,jb)))
        IF (ABS(qc_rttov(jc,jk,jb)) < 1.e-8_wp) qc_rttov(jc,jk,jb) = 0._wp
        qcc_rttov(jc,jk,jb) = 1.e3_wp*(qcc_vi_rttov(jc,jk,jb)-qcc_vi_rttov(jc,jk+1,jb)) / &
                              (grav*(z3d_rttov(jc,jk,jb)-z3d_rttov(jc,jk+1,jb)))
        IF (ABS(qcc_rttov(jc,jk,jb)) < 1.e-8_wp) qcc_rttov(jc,jk,jb) = 0._wp
        qi_rttov(jc,jk,jb)  = 1.e3_wp*(qi_vi_rttov(jc,jk,jb)-qi_vi_rttov(jc,jk+1,jb)) / &
                              (grav*(z3d_rttov(jc,jk,jb)-z3d_rttov(jc,jk+1,jb)))
        IF (ABS(qi_rttov(jc,jk,jb)) < 1.e-8_wp) qi_rttov(jc,jk,jb) = 0._wp
        qs_rttov(jc,jk,jb)  = 1.e3_wp*(qs_vi_rttov(jc,jk,jb)-qs_vi_rttov(jc,jk+1,jb)) / &
                              (grav*(z3d_rttov(jc,jk,jb)-z3d_rttov(jc,jk+1,jb)))
        IF (ABS(qs_rttov(jc,jk,jb)) < 1.e-8_wp) qs_rttov(jc,jk,jb) = 0._wp
        clc_rttov(jc,jk,jb) = MIN(1._wp,(clc_vi_rttov(jc,jk,jb)-clc_vi_rttov(jc,jk+1,jb)) / &
                              (100._wp*(pres_rttov(jk+1)-pres_rttov(jk))) )
        IF (ABS(clc_rttov(jc,jk,jb)) < 1.e-8_wp) clc_rttov(jc,jk,jb) = 0._wp
        ! Fix clouds at levels below ground
        IF (rg_psfc(jc,jb) < pres_rttov(jk)) THEN
          qc_rttov(jc,jk,jb)  = 0._wp
          qcc_rttov(jc,jk,jb) = 0._wp
          qi_rttov(jc,jk,jb)  = 0._wp
          qs_rttov(jc,jk,jb)  = 0._wp
          clc_rttov(jc,jk,jb) = 0._wp  
        ENDIF
      ENDDO
    ENDDO
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jc = i_startidx, i_endidx
      qc_rttov(jc,nlev_rttov,jb)  = 0._wp
      qcc_rttov(jc,nlev_rttov,jb) = 0._wp
      qi_rttov(jc,nlev_rttov,jb)  = 0._wp
      qs_rttov(jc,nlev_rttov,jb)  = 0._wp
      clc_rttov(jc,nlev_rttov,jb) = 0._wp
    ENDDO
    !$ACC END PARALLEL
  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  !$ACC WAIT
  !$ACC END DATA

END SUBROUTINE prepare_rttov_input



!>
!! Back-interpolation of synthetic satellite images from reduced grid to compute grid
!!
SUBROUTINE downscale_rttov_output(jg, nimg, rg_satimg, satimg, l_enabled)

  ! Input grid parameters
  INTEGER, INTENT(IN)  :: jg     ! domain ID of main grid
  INTEGER, INTENT(IN)  :: nimg   ! number of sat images (= number of "levels" in satimg field)

  ! Array with sat images on reduced grid
  REAL(wp), INTENT(INOUT) :: rg_satimg(:,:,:)

  ! Array with sat images on full grid
  REAL(wp), INTENT(OUT) :: satimg(:,:,:)

  ! LOGICAL field: skip image if .FALSE.
  LOGICAL, INTENT(IN) :: l_enabled(:)

  ! Indices
  INTEGER :: i_chidx, jb, jk, jc, i_startblk, i_endblk, i_startidx, i_endidx
  INTEGER :: jc1, jc2, jc3, jc4, jb1, jb2, jb3, jb4
  INTEGER :: rl_start, rl_end

  ! Convenience pointers
  TYPE(t_patch),         POINTER :: p_pp
  TYPE(t_int_state),     POINTER :: p_int
  TYPE(t_gridref_state), POINTER :: p_grf

  INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk

  p_pp  => p_patch_local_parent(jg)
  p_int => p_int_state_local_parent(jg)
  p_grf => p_grf_state_local_parent(jg)

  ! pointers to child index/block
  iidx => p_pp%cells%child_idx
  iblk => p_pp%cells%child_blk

  i_chidx  = p_patch(jg)%parent_child_index

  ! Synchronize sat image array before interpolation
  IF (dbg_level > 2) THEN
    CALL p_barrier(p_comm_work)
    WRITE (0,*) "Synchronize sat image array before interpolation"
  END IF
  CALL exchange_data(p_pat=p_pp%comm_pat_c, lacc=.TRUE., recv=rg_satimg)

  ! Execute interpolation from reduced grid to full grid
  IF (dbg_level > 2) THEN
    CALL p_barrier(p_comm_work)
    WRITE (0,*) "Execute interpolation from reduced grid to full grid"
  END IF
  CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), 0, 1, 1, &
    &                         rg_satimg, satimg, overshoot_fac=1.0_wp,opt_l_enabled=l_enabled)
  IF (dbg_level > 2) THEN
    CALL p_barrier(p_comm_work)
    WRITE (0,*) "CALL to interpol_scal_nudging done."
  END IF

  ! Fill nest boundary points by copying the values from the reduced grid to the full grid
  rl_start = -1
  rl_end   = -2
  i_startblk = p_pp%cells%start_block(rl_start)
  i_endblk   = p_pp%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jc1,jc2,jc3,jc4,jb1,jb2,jb3,jb4)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk,       &
                       i_startidx, i_endidx, rl_start, rl_end)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR PRIVATE(jc1, jc2, jc3, jc4, jb1, jb2, jb3, jb4)
    DO jc = i_startidx, i_endidx
      jc1 = iidx(jc,jb,1)
      jc2 = iidx(jc,jb,2)
      jc3 = iidx(jc,jb,3)
      jc4 = iidx(jc,jb,4)
      jb1 = iblk(jc,jb,1)
      jb2 = iblk(jc,jb,2)
      jb3 = iblk(jc,jb,3)
      jb4 = iblk(jc,jb,4)

      DO jk = 1,nimg
        IF (.NOT. l_enabled(jk))  CYCLE
        satimg(jc1,jk,jb1) = rg_satimg(jc,jk,jb)
        satimg(jc2,jk,jb2) = rg_satimg(jc,jk,jb)
        satimg(jc3,jk,jb3) = rg_satimg(jc,jk,jb)
        satimg(jc4,jk,jb4) = rg_satimg(jc,jk,jb)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

END SUBROUTINE downscale_rttov_output


!>
!! Define the pressure levels (interface/half levels, to be precise) used by RTTOV
!! Unit is hPa!
!!
SUBROUTINE define_rttov_levels (nlev_rt, pres_rttov)

  INTEGER, INTENT(IN) :: nlev_rt
  REAL(wp), INTENT(OUT) :: pres_rttov(:)


      SELECT CASE (nlev_rt)
      CASE (43) ! RTTOV9 levels and RTTOV10 levels (in RTTOV-9 compat mode,
                ! using RTTOV9 coefficients)
        pres_rttov =                                                        &
         (/           0.10,   0.29,   0.69,   1.42,  2.611,  4.407,   6.95, &
             10.37,  14.81,  20.40,  27.26,  35.51,  45.29,  56.73,  69.97, &
             85.18, 102.05, 122.04, 143.84, 167.95, 194.36, 222.94, 253.71, &
            286.60, 321.50, 358.28, 396.81, 436.95, 478.54, 521.46, 565.54, &
            610.60, 656.43, 702.73, 749.12, 795.09, 839.95, 882.80, 922.46, &
            957.44, 985.88,1005.43,1013.25  /)
      CASE (44) ! RTTOV10 levels (using RTTOV9 coefficients)
        pres_rttov =                                                        &
         (/  0.005,   0.10,   0.29,   0.69,   1.42,  2.611,  4.407,   6.95, &
             10.37,  14.81,  20.40,  27.26,  35.51,  45.29,  56.73,  69.97, &
             85.18, 102.05, 122.04, 143.84, 167.95, 194.36, 222.94, 253.71, &
            286.60, 321.50, 358.28, 396.81, 436.95, 478.54, 521.46, 565.54, &
            610.60, 656.43, 702.73, 749.12, 795.09, 839.95, 882.80, 922.46, &
            957.44, 985.88,1005.43,1013.25  /)
      CASE (50) ! RTTOV-10 levels (in RTTOV-9 compat mode)
        pres_rttov =                                                         &
         (/             0.01,    0.10,    0.20,     0.50,    0.80,    1.20,  &
               1.60,    2.20,    2.70,    3.50,     4.20,    5.00,    6.95,  &
              10.37,   14.81,   20.40,   27.26,    35.51,   45.29,   56.73,  &
              69.97,   85.18,  102.05,  122.04,   143.84,  167.95,  194.36,  &
             222.94,  253.71,  286.60,  321.50,   358.38,  396.81,  436.95,  &
             478.54,  521.46,  565.54,  610.60,   656.43,  702.73,  749.12,  &
             795.09,  839.95,  882.80,  922.46,   957.44,  985.88, 1005.43,  &
            1025.00, 1050.00 /)
      CASE (51) ! RTTOV-10 levels
        pres_rttov =                                                         &
         (/    0.005,   0.01,    0.10,    0.20,     0.50,    0.80,    1.20,  &
               1.60,    2.20,    2.70,    3.50,     4.20,    5.00,    6.95,  &
              10.37,   14.81,   20.40,   27.26,    35.51,   45.29,   56.73,  &
              69.97,   85.18,  102.05,  122.04,   143.84,  167.95,  194.36,  &
             222.94,  253.71,  286.60,  321.50,   358.38,  396.81,  436.95,  &
             478.54,  521.46,  565.54,  610.60,   656.43,  702.73,  749.12,  &
             795.09,  839.95,  882.80,  922.46,   957.44,  985.88, 1005.43,  &
            1025.00, 1050.00 /)
      CASE (53) ! RTTOV-10 levels (in RTTOV-9 compat mode)
        pres_rttov =                                                         &
         (/             0.0131,    0.0304,    0.0644,    0.1263,    0.2324,  &
             0.4052,    0.6749,    1.0801,    1.6691,    2.5011,    3.6462,  &
             5.1864,    7.2150,    9.8368,   13.1672,   17.3308,   22.4601,  &
            28.6937,   36.1735,   45.0430,   55.4433,   67.5109,   81.3744,  &
            97.1505,  114.9420,  134.8320,  156.8850,  181.1390,  207.6090,  &
           236.2780,  267.1010,  300.0000,  334.8650,  371.5530,  409.8890,  &
           449.6680,  490.6520,  532.5770,  575.1540,  618.0710,  660.9960,  &
           703.5860,  745.4840,  786.3280,  825.7550,  863.4050,  898.9280,  &
           931.9850,  962.2590,  989.4510, 1013.2900, 1033.5400, 1050.0000  /)
      CASE (54) ! RTTOV-10 levels
        pres_rttov =                                                         &
         (/  0.0050,    0.0131,    0.0304,    0.0644,    0.1263,    0.2324,  &
             0.4052,    0.6749,    1.0801,    1.6691,    2.5011,    3.6462,  &
             5.1864,    7.2150,    9.8368,   13.1672,   17.3308,   22.4601,  &
            28.6937,   36.1735,   45.0430,   55.4433,   67.5109,   81.3744,  &
            97.1505,  114.9420,  134.8320,  156.8850,  181.1390,  207.6090,  &
           236.2780,  267.1010,  300.0000,  334.8650,  371.5530,  409.8890,  &
           449.6680,  490.6520,  532.5770,  575.1540,  618.0710,  660.9960,  &
           703.5860,  745.4840,  786.3280,  825.7550,  863.4050,  898.9280,  &
           931.9850,  962.2590,  989.4510, 1013.2900, 1033.5400, 1050.0000  /)
      ! case (100) hires coefficients. only available for iasi, airs, cris
      ! case (101)
      CASE default
        call finish('define_rttov_levels', 'failed to determine pressure levels')
      END SELECT

END SUBROUTINE define_rttov_levels


!>
!! This routine copies additional model levels to the local parent grid if the RTTOV library
!! is called on a vertically nested grid
!!
SUBROUTINE copy_rttov_ubc (jg, jgc, prm_diag, lacc)
  ! Input grid parameters
  INTEGER,              INTENT(in) :: jg, jgc
  TYPE(t_nwp_phy_diag), INTENT(in) :: prm_diag(:)
  LOGICAL, OPTIONAL,    INTENT(in) :: lacc ! If true, use openacc

  ! Local fields

  INTEGER :: nshift

  CALL assert_acc_device_only("copy_rttov_ubc", lacc=lacc)

  nshift = p_patch(jgc)%nshift
  CALL exchange_data_mult(p_pat=p_patch_local_parent(jgc)%comm_pat_glb_to_loc_c,                                   &
       lacc=.TRUE.,                                                                                                &
       nfields=5, ndim2tot=5*nshift,                                                                               &
       RECV1=prm_diag(jgc)%buffer_rttov(:,         1:  nshift,:), SEND1=p_nh_state(jg)%metrics%z_ifc(:,1:nshift,:),&
       RECV2=prm_diag(jgc)%buffer_rttov(:,  nshift+1:2*nshift,:), SEND2=p_nh_state(jg)%diag%pres(:,1:nshift,:),    &
       RECV3=prm_diag(jgc)%buffer_rttov(:,2*nshift+1:3*nshift,:), SEND3=p_nh_state(jg)%diag%dpres_mc(:,1:nshift,:),&
       RECV4=prm_diag(jgc)%buffer_rttov(:,3*nshift+1:4*nshift,:), SEND4=p_nh_state(jg)%diag%temp(:,1:nshift,:),    &
       RECV5=prm_diag(jgc)%buffer_rttov(:,4*nshift+1:5*nshift,:), SEND5=prm_diag(jg)%tot_cld(:,1:nshift,:,iqv)     )
END SUBROUTINE copy_rttov_ubc


END MODULE mo_rttov_interface



