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

! Module to provide interface to radiation routines.
!
! Remarks
!   This module contains routines that provide the interface between ECHAM
!   and the radiation code.  Mostly it organizes and calculates the
!   information necessary to call the radiative transfer solvers.
!
!         Hauke Schmidt, MPI-M, Hamburg (2009-12-18): Few modifications to
!              allow specific solar irradiance for AMIP-type and preindustrial
!              simulations.
!         Luis Kornblueh, MPI-M, Hamburg (2010-04-06): Never ever use write
!              directly
!         Martin Schultz, FZJ, Juelich (2010-04-13):
!              Extracted public parameters into new module mo_radiation_parameters
!              to avoid circular dependencies in submodels
!                                      (2010-06-03):
!              Added submodel calls, decl_sun_cur
!         Roland Wirth, DWD, Offenbach (2021-09-16):
!              Add output fluxes for near-IR, visible, and diffuse radiation.
!
! Origin
!   Major segments of this code combines and rewrites (for the ICON standard)
!   code previously contained in the ECHAM5 routines rad_int.f90,
!   radiation.f90 and prerad.f90.  Modifications were also made to provide
!   a cleaner interface to the aerosol and cloud properties. 

MODULE mo_radiation

  USE mo_bc_greenhouse_gases,  ONLY: ghg_co2mmr, ghg_ch4mmr, ghg_n2ommr, ghg_cfcmmr
  USE mo_aerosol_util,         ONLY: zaea_rrtm,zaes_rrtm,zaeg_rrtm
  USE mo_kind,                 ONLY: wp, i8
  USE mo_exception,            ONLY: finish

  USE mo_model_domain,         ONLY: t_patch
  USE mo_nonhydro_state,       ONLY: p_nh_state

  USE mo_math_constants,       ONLY: pi, rpi, rad2deg, dbl_eps
  USE mo_math_types,           ONLY: t_geographical_coordinates
  USE mo_physical_constants,   ONLY: grav,  rd,    avo,   amd,  amw,  &
    &                                amco2, amch4, amn2o, amo3, amo2, &
    &                                amc11, amc12, stbo
  USE mo_time_config,          ONLY: time_config
  USE mo_radiation_config,     ONLY: tsi_radt,   ssi_radt,            &
    &                                irad_co2,   mmr_co2,             &
    &                                irad_ch4,   mmr_ch4,   vpp_ch4,  &
    &                                irad_n2o,   mmr_n2o,   vpp_n2o,  &
    &                                irad_o2,    mmr_o2,              &
    &                                irad_cfc11, mmr_cfc11,           &
    &                                irad_cfc12, mmr_cfc12,           &
    &                                isolrad, izenith,                &
    &                                cos_zenith_fixed, islope_rad,    &
    &                                irad_aero, iRadAeroNone,         &
    &                                iRadAeroConst, iRadAeroTegen,    &
    &                                iRadAeroART, iRadAeroConstKinne, &
    &                                iRadAeroKinne, iRadAeroVolc,     &
    &                                iRadAeroKinneVolc
  USE mo_lnd_nwp_config,       ONLY: isub_seaice, isub_lake, isub_water
  USE mo_extpar_config,        ONLY: ext_atm_attr
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_newcld_optics,        ONLY: newcld_optics
  USE mo_bc_aeropt_kinne,      ONLY: set_bc_aeropt_kinne
  USE mo_bc_aeropt_cmip6_volc, ONLY: add_bc_aeropt_cmip6_volc
  USE mo_lrtm_par,             ONLY: jpband => nbndlw, jpxsec => maxxsec
  USE mo_lrtm,                 ONLY: lrtm
  USE mo_srtm_config,          ONLY: jpsw, jpinpx
  USE mo_srtm,                 ONLY: srtm_srtm_224gp
  USE mo_timer,                ONLY: timers_level, timer_start, timer_stop,  &
    &                                timer_radiation, timer_rrtm_prep,       &
    &                                timer_lrtm, timer_srtm, timer_rrtm_post
  USE mo_nh_testcases_nml,     ONLY: zenithang
#ifdef __ICON_ART
  USE mo_art_radiation_interface, ONLY: art_rad_aero_interface
#endif
  USE mtime,                   ONLY: datetime, newDatetime, timedelta, newTimedelta, &
       &                             getPTStringFromMS, OPERATOR(+),                 &
       &                             NO_OF_MS_IN_A_MINUTE, NO_OF_MS_IN_A_HOUR,       &
       &                             getDayOfYearFromDatetime, MAX_TIMEDELTA_STR_LEN,&
       &                             deallocateTimedelta, deallocateDatetime,        &
       &                             NO_OF_MS_IN_A_SECOND
  
  USE mo_grid_config,          ONLY: l_scm_mode
  USE mo_scm_nml,              ONLY: lon_scm, lat_scm 
  USE mo_fortran_tools,        ONLY: set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pre_radiation_nwp
  PUBLIC :: radiation_nwp
  PUBLIC :: radheat
  PUBLIC :: pre_radiation_nwp_steps


  ! --- radiative transfer parameters
  !
  REAL(wp), PARAMETER :: diff   = 1.66_wp   !< LW Diffusivity Factor

CONTAINS

  SUBROUTINE pre_radiation_nwp_steps( &
    & kbdim,cosmu0_dark,p_inc_rad,p_inc_radheat,p_sim_time,pt_patch,zsmu0,zsct, lacc)

    INTEGER, INTENT(IN)   :: &
      & kbdim

    REAL(wp), INTENT(IN)   :: &
      & cosmu0_dark, &
      & p_inc_rad, &                                                   !radiation time step in seconds
      & p_inc_radheat, &
      & p_sim_time

    TYPE(t_patch), INTENT(IN), TARGET :: pt_patch                      ! patch

    REAL(wp), INTENT(OUT), OPTIONAL   :: zsct                          ! solar constant (at time of year)
    REAL(wp), INTENT(OUT)             :: zsmu0(kbdim,pt_patch%nblks_c) ! cosine of zenith angle
    LOGICAL                           :: lacc ! accelerator flag
    LOGICAL                           :: l_tsi_recalculated

    REAL(wp) ::                    &
      & p_sim_time_rad,            &
      & zstunde,                   &
      & ztwo, ztho  ,              &
      & zdek,                      &
      & zsocof, zeit0,             &
      & zsct_h

    REAL(wp), SAVE ::              &
      & zsct_save, zdtzgl,         &
      & zdeksin,zdekcos


    REAL(wp) ::                          &
      & zsinphi(kbdim,pt_patch%nblks_c) ,&
      & zcosphi(kbdim,pt_patch%nblks_c) ,&
      & zeitrad(kbdim,pt_patch%nblks_c) ,&
      & z_cosmu0(kbdim,pt_patch%nblks_c)

    INTEGER :: jj, itaja

    INTEGER :: ie,jb,jc,jmu0,n_zsct,nsteps

    INTEGER :: n_cosmu0pos(kbdim,pt_patch%nblks_c)

    INTEGER , SAVE :: itaja_zsct_previous = 0

    TYPE(datetime), POINTER :: current => NULL()
    TYPE(timedelta), POINTER :: td => NULL()
    CHARACTER(len=MAX_TIMEDELTA_STR_LEN) :: td_string 
        
    TYPE(t_geographical_coordinates), TARGET, ALLOCATABLE :: scm_center(:,:)
    TYPE(t_geographical_coordinates), POINTER             :: ptr_center(:,:)


#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: zsinphi,zcosphi,zeitrad,z_cosmu0,n_cosmu0pos
#endif

    ! In case of CMIP irradiation (isolrad==2) is used, tsi_rad changes during the day
    ! (ssi_time_interpolation), in which case it makes sense to redo the scaling
    ! by orbit elements (recalculate zsct) at each radiation time step and not
    ! just on the first radiation step of the day.
    l_tsi_recalculated = .FALSE.
    IF (isolrad==2) l_tsi_recalculated = .TRUE.

    ! SCM: read lat/lon for horizontally uniform zenith angle
    IF ( l_scm_mode ) THEN
#ifdef _OPENACC
      if (lacc) CALL finish("Single column mode is untested with OpenACC")
#endif
      ALLOCATE(scm_center(SIZE(pt_patch%cells%center,1),SIZE(pt_patch%cells%center,2)))
      !$ACC ENTER DATA CREATE(scm_center) IF(lacc)
      DO jb = 1,pt_patch%nblks_c
        ie = MERGE(kbdim, pt_patch%npromz_c, jb /= pt_patch%nblks_c)
        !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lacc)
        scm_center(1:ie,jb)%lat = lat_scm * pi/180.
        scm_center(1:ie,jb)%lon = lon_scm * pi/180.
        !$ACC END KERNELS
      ENDDO
      ptr_center => scm_center
    ELSE
      ptr_center => pt_patch%cells%center
    ENDIF

    IF (izenith == 0) THEN
    ! local insolation = constant = global mean insolation (ca. 340 W/m2)
    ! zenith angle = 0,
#ifdef _OPENACC
      IF (lacc) CALL finish('pre_radiation_nwp','Only ported on gpu for izenith == 3 and 4')
#endif
      DO jb = 1, pt_patch%nblks_c
        ie = MERGE(kbdim, pt_patch%npromz_c, jb /= pt_patch%nblks_c)
        zsmu0(1:ie,jb) = 1._wp ! sun in zenith everywhere
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi_radt/4._wp ! scale ztsi to get the correct global mean insolation
    ELSEIF(izenith == 1) THEN
    ! circular non-seasonal orbit,
    ! perpetual equinox,
    ! no diurnal cycle,
    ! local time always 12:00
    ! --> sin(time of day)=1 ) and zenith angle depends on latitude only
#ifdef _OPENACC
      IF (lacc) CALL finish('pre_radiation_nwp','Only ported on gpu for izenith == 3 and 4')
#endif
      DO jb = 1, pt_patch%nblks_c
        ie = MERGE(kbdim, pt_patch%npromz_c, jb /= pt_patch%nblks_c)
        zsmu0(1:ie,jb) = COS( ptr_center(1:ie,jb)%lat )
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi_radt * rpi ! because sun is always in local noon, the TSI needs to be
      ! scaled by 1/pi to get the correct global mean insolation
    ELSEIF (izenith == 2) THEN
    ! circular non-seasonal orbit,
    ! perpetual equinox,
    ! no diurnal cycle,
    ! local time always  07:14:15 or 16:45:45
    ! --> sin(time of day)=1/pi and zenith angle depends on latitude only
#ifdef _OPENACC
      IF (lacc) CALL finish('pre_radiation_nwp','Only ported on gpu for izenith == 3 and 4')
#endif
      DO jb = 1, pt_patch%nblks_c
        ie = MERGE(kbdim, pt_patch%npromz_c, jb /= pt_patch%nblks_c)
        zsmu0(1:ie,jb) = COS( ptr_center(1:ie,jb)%lat ) * rpi
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi_radt
    ELSEIF (izenith == 3) THEN  !Second: case izenith==3 (time (but no date) needed)
    ! circular non-seasonal orbit,
    ! perpetual equinox,
    ! with diurnal cycle,

      !$ACC DATA CREATE(n_cosmu0pos, z_cosmu0) COPYOUT(zsmu0) IF(lacc)

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lacc)
      zsmu0(:,:)=0.0_wp
      n_cosmu0pos(:,:) = 0
      !$ACC END KERNELS
      nsteps = NINT(p_inc_rad/p_inc_radheat)

      DO jmu0=1,nsteps

        p_sim_time_rad = p_sim_time + (REAL(jmu0,wp)-0.5_wp)*p_inc_radheat

        current => newDatetime(time_config%tc_exp_startdate)
        CALL getPTStringFromMS(INT(1000.0_wp*p_sim_time_rad,i8), td_string)
        td => newTimedelta(td_string)
        current = time_config%tc_exp_startdate + td
        jj = INT(current%date%year)
        itaja = getDayOfYearFromDateTime(current)
        zstunde = current%time%hour+( &
             &    REAL(current%time%minute*NO_OF_MS_IN_A_MINUTE &
             &        +current%time%second*NO_OF_MS_IN_A_SECOND &
             &        +current%time%ms,wp)/REAL(NO_OF_MS_IN_A_HOUR,wp))
        CALL deallocateDatetime(current)
        CALL deallocateTimedelta(td)

        DO jb = 1, pt_patch%nblks_c
          ie = MERGE(kbdim, pt_patch%npromz_c, jb /= pt_patch%nblks_c)

          !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lacc)
          z_cosmu0(1:ie,jb) = -COS( ptr_center(1:ie,jb)%lat ) &
            & *COS( ptr_center(1:ie,jb)%lon                &
            &      +zstunde/24._wp* 2._wp*pi )
          !$ACC END KERNELS

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
          !$ACC LOOP GANG VECTOR
          DO jc = 1,ie
            IF ( z_cosmu0(jc,jb) > -1.e-5_wp ) THEN
              zsmu0(jc,jb) = zsmu0(jc,jb) + MAX(1.e-3_wp,z_cosmu0(jc,jb))**2
              n_cosmu0pos(jc,jb) = n_cosmu0pos(jc,jb) + 1
            ENDIF
          ENDDO
          !$ACC END PARALLEL

        ENDDO !jb

      ENDDO!jmu0

      DO jb = 1, pt_patch%nblks_c

        ie = MERGE(kbdim, pt_patch%npromz_c, jb /= pt_patch%nblks_c)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
        !$ACC LOOP GANG VECTOR
!DIR$ SIMD
        DO jc = 1,ie
          IF (n_cosmu0pos(jc,jb) > 0) THEN
            zsmu0(jc,jb) = SQRT(zsmu0(jc,jb)/REAL(n_cosmu0pos(jc,jb),wp))
          ELSE
            zsmu0(jc,jb) = cosmu0_dark
          ENDIF
        ENDDO
        !$ACC END PARALLEL

      ENDDO !jb

      IF (PRESENT(zsct)) zsct = tsi_radt
      !$ACC WAIT(1)
      !$ACC END DATA

    ELSEIF (izenith == 4) THEN
    ! elliptical seasonal orbit,
    !  with diurnal cycle

      !$ACC DATA CREATE(n_cosmu0pos, z_cosmu0, zsinphi, zcosphi, zeitrad) COPYOUT(zsmu0) PRESENT(ptr_center) IF(lacc)

      zsct_h = 0.0_wp
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lacc)
      zsmu0(:,:)=0.0_wp
      n_cosmu0pos(:,:) = 0
      !$ACC END KERNELS
      n_zsct = 0

      nsteps = NINT(p_inc_rad/p_inc_radheat)

      DO jmu0=1,nsteps

        p_sim_time_rad = p_sim_time + (REAL(jmu0,wp)-0.5_wp)*p_inc_radheat

        current => newDatetime(time_config%tc_exp_startdate)
        CALL getPTStringFromMS(INT(1000.0_wp*p_sim_time_rad,i8), td_string)
        td => newTimedelta(td_string)
        current = time_config%tc_exp_startdate + td
        jj = INT(current%date%year)
        itaja = getDayOfYearFromDateTime(current)
        zstunde = current%time%hour+( &
             &    REAL(current%time%minute*NO_OF_MS_IN_A_MINUTE &
             &        +current%time%second*NO_OF_MS_IN_A_SECOND &
             &        +current%time%ms,wp)/REAL(NO_OF_MS_IN_A_HOUR,wp))
        CALL deallocateDatetime(current)
        CALL deallocateTimedelta(td)

        IF ( itaja /= itaja_zsct_previous .OR. l_tsi_recalculated ) THEN
          itaja_zsct_previous = itaja

          ztwo    = 0.681_wp + 0.2422_wp*REAL(jj-1949,wp)-REAL((jj-1949)/4,wp)
          ztho    = 2._wp*pi*( REAL(itaja, wp) -1.0_wp + ztwo )/365.2422_wp
          zdtzgl  = 0.000075_wp + 0.001868_wp*COS(      ztho) - 0.032077_wp*SIN(      ztho) &
            - 0.014615_wp*COS(2._wp*ztho) - 0.040849_wp*SIN(2._wp*ztho)
          zdek    = 0.006918_wp - 0.399912_wp*COS(      ztho) + 0.070257_wp*SIN(      ztho) &
            - 0.006758_wp*COS(2._wp*ztho) + 0.000907_wp*SIN(2._wp*ztho) &
            - 0.002697_wp*COS(3._wp*ztho) + 0.001480_wp*SIN(3._wp*ztho)

          zdeksin = SIN (zdek)
          zdekcos = COS (zdek)

          IF ( PRESENT(zsct) ) THEN

            zsocof  = 1.000110_wp + 0.034221_wp*COS(   ztho) + 0.001280_wp*SIN(   ztho) &
              + 0.000719_wp*COS(2._wp*ztho) + 0.000077_wp*SIN(2._wp*ztho)
            zsct_save = zsocof*tsi_radt
            zsct_h = zsct_h + zsct_save
            n_zsct = n_zsct + 1

          ENDIF

        ENDIF

        zeit0   = pi*(zstunde-12._wp)/12._wp + zdtzgl

        DO jb = 1, pt_patch%nblks_c
          ie = MERGE(kbdim, pt_patch%npromz_c, jb /= pt_patch%nblks_c)

          !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lacc)
          zsinphi(1:ie,jb)      = SIN (ptr_center(1:ie,jb)%lat)
          zcosphi(1:ie,jb)      = SQRT(1.0_wp - zsinphi(1:ie,jb)**2)
          zeitrad(1:ie,jb)      = zeit0 + ptr_center(1:ie,jb)%lon
          z_cosmu0(1:ie,jb)     = zdeksin * zsinphi(1:ie,jb) + zdekcos * zcosphi(1:ie,jb) * &
            COS(zeitrad(1:ie,jb))
          !$ACC END KERNELS

          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lacc)
          DO jc = 1,ie
            IF ( z_cosmu0(jc,jb) > -1.e-5_wp ) THEN
              zsmu0(jc,jb) = zsmu0(jc,jb) + MAX(1.e-3_wp,z_cosmu0(jc,jb))**2
              n_cosmu0pos(jc,jb) = n_cosmu0pos(jc,jb) + 1
            ENDIF
          ENDDO
          !$ACC END PARALLEL LOOP

        ENDDO

      ENDDO !jmu0
      !$ACC WAIT(1)

      DO jb = 1, pt_patch%nblks_c

        ie = MERGE(kbdim, pt_patch%npromz_c, jb /= pt_patch%nblks_c)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
        !$ACC LOOP GANG VECTOR
        DO jc = 1,ie
          IF ( n_cosmu0pos(jc,jb) > 0 ) THEN
            ! The averaged cosine of zenith angle is limited to 0.05 in order to avoid
            ! numerical trouble along the day-night boundary near the model top
            zsmu0(jc,jb) = MAX(0.05_wp,SQRT(zsmu0(jc,jb)/REAL(n_cosmu0pos(jc,jb),wp)))
          ELSE
            zsmu0(jc,jb) = cosmu0_dark
          ENDIF
        ENDDO
        !$ACC END PARALLEL

      ENDDO !jb

      IF ( PRESENT(zsct) ) THEN
        IF ( n_zsct > 0 ) THEN
          zsct = zsct_h/REAL(n_zsct,wp)
        ELSE
          zsct = zsct_save
        ENDIF
      ENDIF
      !$ACC END DATA

    ELSEIF (izenith == 5) THEN
     ! Radiative convective equilibrium
     ! circular non-seasonal orbit,
     ! perpetual equinox,
     ! no diurnal cycle,
     ! the product tsi*cos(zenith angle) should equal 340 W/m2
     ! see Popke et al. 2013 and Cronin 2013
#ifdef _OPENACC
      IF (lacc) CALL finish('pre_radiation_nwp','Only ported on gpu for izenith == 3 and 4')
#endif
      DO jb = 1, pt_patch%nblks_c
        ie = MERGE(kbdim, pt_patch%npromz_c, jb /= pt_patch%nblks_c)
        zsmu0(1:ie,jb) = COS(zenithang*pi/180._wp)
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi_radt ! no rescale tsi was adjstd in atm_phy_nwp w ssi_rce

    ELSEIF (izenith == 6) THEN
     ! Prescribed cos(solar zenith angle), for single column model (SCM)
#ifdef _OPENACC
      IF (lacc) CALL finish('pre_radiation_nwp','Only ported on gpu for izenith == 3 and 4')
#endif
      DO jb = 1, pt_patch%nblks_c
        ie = MERGE(kbdim, pt_patch%npromz_c, jb /= pt_patch%nblks_c)
        zsmu0(1:ie,jb) = cos_zenith_fixed
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi_radt ! no rescale tsi was adjstd in atm_phy_nwp w ssi_rce

    ENDIF


    IF (l_scm_mode) THEN
      DEALLOCATE(scm_center)
      !$ACC WAIT(1)
      !$ACC EXIT DATA DELETE(scm_center)
    ENDIF

  END SUBROUTINE pre_radiation_nwp_steps

  SUBROUTINE pre_radiation_nwp(kbdim,p_inc_rad,p_sim_time,pt_patch,zsmu0,zsct,&
    &                          slope_ang,slope_azi,horizon,cosmu0_slp,shading_mask,lacc)

    INTEGER, INTENT(IN)   :: &
      & kbdim

    REAL(wp), INTENT(IN)   :: &
      & p_inc_rad, & !radiation time step in seconds
      & p_sim_time

    TYPE(t_patch), INTENT(IN), TARGET :: pt_patch    ! Patch

    REAL(wp), INTENT(OUT), OPTIONAL   :: zsct                  ! solar constant (at time of year)
    REAL(wp), INTENT(OUT)             :: zsmu0(kbdim,pt_patch%nblks_c)   ! Cosine of zenith angle
    ! Optional fields for slope-dependent surface radiation: slope angle, slope azimuth, and slope-dependent cosine of zenith angle
    REAL(wp), INTENT(IN), DIMENSION(kbdim,pt_patch%nblks_c) :: slope_ang,slope_azi
    REAL(wp), INTENT(IN), DIMENSION(kbdim,pt_patch%nblks_c,ext_atm_attr(pt_patch%id)%nhori) :: horizon
    REAL(wp), INTENT(OUT),DIMENSION(kbdim,pt_patch%nblks_c) :: cosmu0_slp, shading_mask
    LOGICAL, OPTIONAL,           INTENT(in)   :: lacc            !< GPU flag

! Local arrays
! --------------------------------------------------------------------
    REAL(wp), DIMENSION(kbdim) :: zsinphi, zcosphi, zeitrad, czra, szra, csang, ssang, csazi, ssazi, &
                                  zha_sun,    & ! interpolated horizon angle below the sun [degree]
                                  zphi_sun,   & ! sun azimuth angle [rad]
                                  ztheta_sun, & ! sun elevation angle [rad]
                                  ztheta        ! slope angle [rad]
! Local scalars
! --------------------------------------------------------------------
    REAL(wp) ::                    &
      & p_sim_time_rad, x1, x2,    &
      & zihor, zstunde,            &
      & ztwo  , ztho  ,            &
      & zdtzgl, zdek  ,            &
      & zsocof, zeit0 ,            &
      & zdeksin,zdekcos

    INTEGER :: &
      & k, ii, jc, jj, itaja, jb, ie, jg

    LOGICAL :: &
      & lslope_aspect, lzacc      !switches
    LOGICAL :: l_tsi_recalculated

!------------------------------------------------------------------------------

    INTEGER , SAVE :: itaja_zsct_previous = 0
    REAL(wp), SAVE :: zsct_save

    TYPE(datetime), POINTER :: current => NULL()
    TYPE(timedelta), POINTER :: td => NULL()
    CHARACTER(len=MAX_TIMEDELTA_STR_LEN) :: td_string

    TYPE(t_geographical_coordinates), TARGET, ALLOCATABLE :: scm_center(:,:)
    TYPE(t_geographical_coordinates), POINTER             :: ptr_center(:,:)

    INTEGER :: nhori

    jg = pt_patch%id

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: zsinphi,zcosphi,zeitrad,czra,szra,csang,ssang,csazi,ssazi,zha_sun,zphi_sun,ztheta_sun,ztheta
#endif

    nhori = ext_atm_attr(jg)%nhori   ! number of sectors for horizon

    ! In case of CMIP irradiation (isolrad==2) is used, tsi_rad changes during the day
    ! (ssi_time_interpolation), in which case it makes sense to redo the scaling
    ! by orbit elements (recalculate zsct) at each radiation time step and not
    ! just on the first radiation step of the day.
    l_tsi_recalculated = .FALSE.
    IF (isolrad==2) l_tsi_recalculated = .TRUE.

    ! SCM: read lat/lon for horizontally uniform zenith angle
    IF (l_scm_mode) THEN
      ALLOCATE(scm_center(SIZE(pt_patch%cells%center,1),SIZE(pt_patch%cells%center,2)))
      DO jb = 1,pt_patch%nblks_c
        ie = MERGE(kbdim, pt_patch%npromz_c, jb /= pt_patch%nblks_c)
        scm_center(1:ie,jb)%lat = lat_scm * pi/180.
        scm_center(1:ie,jb)%lon = lon_scm * pi/180.
      ENDDO
      ptr_center => scm_center
    ELSE
      ptr_center => pt_patch%cells%center
    ENDIF

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA CREATE(zsinphi, zcosphi, zeitrad, czra, szra, csang, ssang, csazi, ssazi, zha_sun, zphi_sun, ztheta_sun, ztheta) IF(lzacc)

    !First case: izenith==0 to izenith==2 and izenith=6 (no date and time needed)
    IF (izenith == 0) THEN
     ! for testing: provisional setting of cos(zenith angle) and TSI
     ! The global mean insolation is TSI/4 (ca. 340 W/m2)
      DO jb = 1, pt_patch%nblks_c
        ie = MERGE(kbdim, pt_patch%npromz_c, jb /= pt_patch%nblks_c)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO jc= 1, ie
          zsmu0(jc,jb) = 1._wp ! sun in zenith everywhere
        ENDDO
        !$ACC END PARALLEL
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi_radt/4._wp ! scale ztsi to get the correct global mean insolation
    ELSEIF(izenith == 1) THEN
      ! circular non-seasonal orbit, zenith angle dependent on latitude only,
      ! no diurnal cycle (always at 12:00 local time --> sin(time of day)=1 )
      DO jb = 1, pt_patch%nblks_c
        ie = MERGE(kbdim, pt_patch%npromz_c, jb /= pt_patch%nblks_c)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO jc= 1, ie
          zsmu0(jc,jb) = COS( ptr_center(jc,jb) %lat )
        ENDDO
        !$ACC END PARALLEL
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi_radt/pi ! because sun is always in local noon, the TSI needs to be
                                       ! scaled by 1/pi to get the correct global mean insolation
    ELSEIF (izenith == 2) THEN
      ! circular non-seasonal orbit, no diurnal cycle
      ! at 07:14:15 or 16:45:45 local time (--> sin(time of day)=1/pi )
      DO jb = 1, pt_patch%nblks_c
        ie = MERGE(kbdim, pt_patch%npromz_c, jb /= pt_patch%nblks_c)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO jc= 1, ie
          zsmu0(jc,jb) = COS( ptr_center(jc,jb)%lat ) * rpi
        ENDDO
        !$ACC END PARALLEL
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi_radt

    ELSEIF (izenith == 6) THEN
      ! Prescribed cos(solar zenith angle), for single column model
#ifdef _OPENACC
      IF (lacc) CALL finish('pre_radiation_nwp','Not ported on gpu for izenith == 6')
#endif
      DO jb = 1, pt_patch%nblks_c
        ie = MERGE(kbdim, pt_patch%npromz_c, jb /= pt_patch%nblks_c)
        DO jc= 1, ie
          zsmu0(jc,jb) = cos_zenith_fixed
        ENDDO
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi_radt ! no rescale

    ELSE

      p_sim_time_rad = p_sim_time + 0.5_wp*p_inc_rad

      current => newDatetime(time_config%tc_exp_startdate)
      CALL getPTStringFromMS(INT(1000.0_wp*p_sim_time_rad,i8), td_string)
      td => newTimedelta(td_string)
      current = time_config%tc_exp_startdate + td
      jj = INT(current%date%year)
      itaja = getDayOfYearFromDateTime(current)
      zstunde = current%time%hour+( &
          &    REAL(current%time%minute*NO_OF_MS_IN_A_MINUTE &
          &        +current%time%second*NO_OF_MS_IN_A_SECOND &
          &        +current%time%ms,wp)/REAL(NO_OF_MS_IN_A_HOUR,wp))
      CALL deallocateDatetime(current)
      CALL deallocateTimedelta(td)

      !Second case: izenith==3 (time (but no date) needed)
      IF (izenith == 3) THEN

        DO jb = 1, pt_patch%nblks_c
          ie = MERGE(kbdim, pt_patch%npromz_c, jb /= pt_patch%nblks_c)
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG VECTOR
          DO jc= 1, ie
            zsmu0(jc,jb) = -COS( ptr_center(jc,jb)%lat ) &
              & *COS( ptr_center(jc,jb)%lon                &
              &      +zstunde * (1._wp/24._wp) * 2._wp * pi )
          ENDDO
          !$ACC END PARALLEL
        ENDDO
        IF (PRESENT(zsct)) zsct = tsi_radt

      !Third case: izenith=4 (time and date needed)
      ELSEIF (izenith == 4) THEN

        ztwo    = 0.681_wp + 0.2422_wp*REAL(jj-1949,wp)-REAL((jj-1949)/4,wp)
        ztho    = 2._wp*pi*( REAL(itaja, wp) -1.0_wp + ztwo )/365.2422_wp
        zdtzgl  = 0.000075_wp + 0.001868_wp*COS(      ztho) - 0.032077_wp*SIN(      ztho) &
          - 0.014615_wp*COS(2._wp*ztho) - 0.040849_wp*SIN(2._wp*ztho)
        zdek    = 0.006918_wp - 0.399912_wp*COS(      ztho) + 0.070257_wp*SIN(      ztho) &
          - 0.006758_wp*COS(2._wp*ztho) + 0.000907_wp*SIN(2._wp*ztho) &
          - 0.002697_wp*COS(3._wp*ztho) + 0.001480_wp*SIN(3._wp*ztho)
        zeit0   = pi*(zstunde-12._wp)/12._wp + zdtzgl
        zdeksin = SIN (zdek)
        zdekcos = COS (zdek)

        IF ( PRESENT(zsct) ) THEN
          !decide whether new zsct calculation is necessary
          IF ( itaja /= itaja_zsct_previous .OR. l_tsi_recalculated ) THEN
            itaja_zsct_previous = itaja
            zsocof  = 1.000110_wp + 0.034221_wp*COS(   ztho) + 0.001280_wp*SIN(   ztho) &
              + 0.000719_wp*COS(2._wp*ztho) + 0.000077_wp*SIN(2._wp*ztho)
            zsct_save = zsocof*tsi_radt
          ENDIF
          zsct = zsct_save
        ENDIF

        lslope_aspect = .TRUE.


        IF (islope_rad(jg) >= 2) zihor = REAL(INT(360.0_wp/nhori),wp)

        DO jb = 1, pt_patch%nblks_c
          ie = MERGE(kbdim, pt_patch%npromz_c, jb /= pt_patch%nblks_c)
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG VECTOR
          DO jc = 1, ie
            zsinphi(jc)      = SIN (ptr_center(jc,jb)%lat)
            zcosphi(jc)      = SQRT(1.0_wp - zsinphi(jc)**2)
            zeitrad(jc)      = zeit0 + ptr_center(jc,jb)%lon
            czra   (jc)      = COS(zeitrad(jc))
            zsmu0  (jc,jb)   = zdeksin * zsinphi(jc) + zdekcos * zcosphi(jc) * czra(jc)
          ENDDO
          !$ACC END PARALLEL

          IF (islope_rad(jg) == 1 .OR. islope_rad(jg) == 3) THEN
            !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
            !$ACC LOOP GANG VECTOR
            DO jc = 1, ie
              szra(jc)  = SIN(zeitrad(jc))
              csang(jc) = COS(slope_ang(jc,jb))
              ssang(jc) = SIN(slope_ang(jc,jb))
              csazi(jc) = COS(slope_azi(jc,jb))
              ssazi(jc) = SIN(slope_azi(jc,jb))
              cosmu0_slp(jc,jb) = zdeksin * ( zsinphi(jc)*csang(jc) - zcosphi(jc)*csazi(jc)*ssang(jc) ) + &
                zdekcos * ( zcosphi(jc)*czra(jc)*csang(jc) + zsinphi(jc)*czra(jc)*csazi(jc)*ssang(jc) + &
                szra(jc)*ssazi(jc)*ssang(jc) )
              IF (cosmu0_slp(jc,jb) < 1.e-3_wp) cosmu0_slp(jc,jb) = 0._wp
            ENDDO
            !$ACC END PARALLEL
          ENDIF

          IF (islope_rad(jg) == 3) THEN ! add orographic shading
            !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
            !$ACC LOOP GANG VECTOR PRIVATE(x1, ii, k)
            DO jc = 1, ie
              ztheta_sun(jc) = ASIN(zsmu0(jc,jb)) ! solar elevation angle
              x1 = MIN(1._wp, MAX(-1._wp, (zsmu0(jc,jb)*zsinphi(jc)-zdeksin)/(COS(ztheta_sun(jc))*zcosphi(jc)) ))
              zphi_sun(jc) = SIGN(acos(x1),SIN(zeitrad(jc)))+pi  ! solar azimuth angle

              ! elevation angle of orography in the direction of the sun, obtained
              ! by linear interpolation between the precomputed horizon angles
              ii = MIN(nhori-1, INT(rad2deg*zphi_sun(jc)/zihor))
              k = MOD(ii+1,nhori)
              zha_sun(jc) = (horizon(jc,jb,k+1) *(rad2deg*zphi_sun(jc)-zihor*ii) +     &
                             horizon(jc,jb,ii+1)*(zihor*(ii+1)-rad2deg*zphi_sun(jc))) / zihor

              ! compute shading mask
              shading_mask(jc,jb) = MERGE(0._wp, 1._wp, zha_sun(jc) > rad2deg*ztheta_sun(jc))
            ENDDO
            !$ACC END PARALLEL
          ENDIF

          IF (islope_rad(jg) == 2) THEN
            !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
            !$ACC LOOP GANG VECTOR PRIVATE(x1, x2, ii, k)
            DO jc = 1, ie
              szra(jc)  = SIN(zeitrad(jc))

              ! The following code is from COSMO V5.10 in subroutine calc_rad_corrections of
              ! radiation_utilities.f90 where :
              ! thetain  (:) is replaced by the slope_ang(jc,jb)
              ! phi      (:) is replaced by the slope_azi(jc,jb)

              ! sun elevation angle
              ztheta_sun(jc) = ASIN(zsmu0(jc,jb))

              ! sun azimuth angle
              x1 = zdekcos * szra(jc) / COS(ztheta_sun(jc))
              x2 = ( zsinphi(jc) * zdekcos * czra(jc) -   &
                   & zcosphi(jc) * zdeksin ) / COS(ztheta_sun(jc))

              IF (x2 < -1.0_wp) x2 = -1.0_wp
              IF (x2 >  1.0_wp) x2 =  1.0_wp
              zphi_sun(jc) = ACOS(x2)
              IF (x1 < 0) zphi_sun(jc) = - zphi_sun(jc)
              zphi_sun(jc) = zphi_sun(jc) + pi

              ! sun elevation angle corrected by refraction (empiric formula, A.H.Siemer(1988))
              ztheta_sun(jc) = ztheta_sun(jc) + (1.569000_wp - ztheta_sun(jc)) / &
                                 &     (185.5_wp + 3620.0_wp * ztheta_sun(jc))

              ! night or day
              IF (ztheta_sun(jc) < 0.0_wp) THEN
                ztheta_sun(jc) = 0.0_wp
              ENDIF

              ! the horizon has a spatial resolution of 360/nhori degrees.
              ! a distance weighted linear interpolation is 
              ! done between the two neighboring points.
              ii = INT(rad2deg*zphi_sun(jc)/zihor)
              IF (ii >= nhori) THEN
                ii = nhori - 1
              ENDIF

              k = MOD(ii+1,nhori)

              zha_sun(jc) = (horizon(jc,jb,k+1) *(rad2deg*zphi_sun(jc)-zihor*ii) +     &
                &            horizon(jc,jb,ii+1)*(zihor*(ii+1)-rad2deg*zphi_sun(jc)))/ &
                &            zihor

              ! compute shadowmask
              shading_mask(jc,jb) = MERGE(0._wp, 1._wp, zha_sun(jc) > rad2deg*ztheta_sun(jc))
 
              ! compute correction factor and multiply by sin(ztheta_sun) to get cosmu0_slp
              ! slope angle and aspect switched off
              IF (.NOT. lslope_aspect) THEN
                ztheta(jc) = 0.0_wp
              ELSE
                ztheta(jc) = slope_ang(jc,jb) ! thetain(ip)
              ENDIF

              IF (ztheta_sun(jc) > 0.01_wp) THEN
              ! Mueller and Scherer (2005) formula (MWR)
              ! New formula (lower correction, theoretically correct derived)
                cosmu0_slp(jc,jb) = SIN(ztheta_sun(jc)) *                    &
                  ( COS(ztheta(jc)) + (SIN(ztheta(jc))/TAN(ztheta_sun(jc)))* &
                                        COS(zphi_sun(jc) - slope_azi(jc,jb)) )
              ELSE
                cosmu0_slp(jc,jb) = SIN(ztheta_sun(jc))
              ENDIF

              ! Consistency check to avoid negative corrections:
              ! active in situations with low sun elevation (slope angles > sun elevation)
              ! when the slope aspect is greater than the daily maxima or smaller than the
              ! daily minima of the sun azimuth angle (during the sunshine time, a kind
              ! of self shading effect).
              IF (cosmu0_slp(jc,jb) < 0.0_wp) THEN
                cosmu0_slp(jc,jb) = 0.0_wp
              ENDIF

            ENDDO
            !$ACC END PARALLEL

          ENDIF
        ENDDO

      ELSEIF (izenith == 5) THEN
      ! Radiative convective equilibrium
      ! circular non-seasonal orbit,
      ! perpetual equinox,
      ! no diurnal cycle,
      ! the product tsi*cos(zenith angle) should equal 340 W/m2
      ! see Popke et al. 2013 and Cronin 2013
        DO jb = 1, pt_patch%nblks_c
          ie = MERGE(kbdim, pt_patch%npromz_c, jb /= pt_patch%nblks_c)
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG VECTOR
          DO jc = 1, ie
            zsmu0(jc,jb) = COS(zenithang*pi/180._wp)
          ENDDO
          !$ACC END PARALLEL
        ENDDO
        IF (PRESENT(zsct)) zsct = tsi_radt ! no rescale tsi was adjstd in atm_phy_nwp w ssi_rce

      ENDIF
  
    ENDIF

    ! Avoid uninitialized output fields
    DO jb = 1, pt_patch%nblks_c
      ie = MERGE(kbdim, pt_patch%npromz_c, jb /= pt_patch%nblks_c)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = 1, ie
        IF (islope_rad(jg) <= 1) shading_mask(jc,jb) = 1._wp
        IF (islope_rad(jg) == 0) cosmu0_slp(jc,jb) = 0._wp ! not used
      ENDDO
      !$ACC END PARALLEL
    ENDDO

    IF (l_scm_mode) THEN
      DEALLOCATE(scm_center)
    ENDIF
    !$ACC WAIT
    !$ACC END DATA ! CREATE(zsinphi,zcosphi,zeitrad,czra,szra,csang,ssang,csazi,ssazi,zha_sun,zphi_sun,ztheta_sun,ztheta)

  END SUBROUTINE pre_radiation_nwp



  !-----------------------------------------------------------------------------
  !>
  !! Organizes the calls to the radiation solver
  !!
  !! Remarks This routine organises the input/output for the radiation
  !! computation.  The state of radiatively active constituents is set as the
  !! input. Output are flux transmissivities and emissivities at all the half
  !! levels of the grid (respectively ratio solar flux/solar input and ratio
  !! thermal flux/local black-body flux). This output will be used in radheat
  !! at all time steps until the next full radiation time step.
  !
  SUBROUTINE radiation_nwp(                                                &
    ! input
    &  current_date                                                        &
    & ,jg, jb, icpl_reff                                                   &
    & ,jcs               ,jce               ,kbdim                         &
    & ,klev              ,klevp1                                           &
    & ,ktype             ,zland           ,zglac            ,cos_mu0       &
    & ,alb_vis_dir       ,alb_nir_dir     ,alb_vis_dif      ,alb_nir_dif   &
    & ,emis_rad                                                            &
    & ,tk_sfc            ,pp_hl            ,pp_fl                          &
    & ,tk_fl             ,qm_vap          ,qm_liq           ,qm_ice        &
    & ,qm_o3                                                               &
    & ,cdnc              ,reff_liq        ,reff_frz         ,cld_frc       &
    & ,zaeq1, zaeq2, zaeq3, zaeq4, zaeq5, dust_tunefac                     &
    ! output
    & ,cld_cvr, flx_lw_net, flx_uplw_sfc, trsol_net, trsol_up_toa          &
    & ,trsol_up_sfc, trsol_dn_sfc_diffus, trsol_clr_sfc, trsol_nir_sfc     &
    & ,trsol_vis_sfc, trsol_par_sfc, fr_nir_sfc_diffus                     &
    & ,fr_vis_sfc_diffus, fr_par_sfc_diffus                                &
    & ,lwflx_clr_sfc                                                       &
    ! optional output: 3D flux output
    & ,flx_lw_dn      ,flx_sw_dn     ,flx_lw_up     ,flx_sw_up             &
    & ,flx_lw_dn_clr  ,flx_sw_dn_clr ,flx_lw_up_clr ,flx_sw_up_clr         )

    ! input
    ! -----
    !
    TYPE(datetime), POINTER, INTENT(in) :: current_date !< current date
    INTEGER, INTENT(in)   :: &
      &  jg,                 & !< domain index
      &  jb,                 & !< block index
      &  icpl_reff,          & !< option for couplig radiation with reff
      &  jcs,                & !< start index for loop over block
      &  jce,                & !< end   index for loop over block
      &  kbdim,              & !< dimension of block over cells
      &  klev,               & !< number of full levels = number of layers
      &  klevp1                !< number of half levels = number of layer interfaces
    INTEGER, INTENT(in)   :: &
      &  ktype(kbdim)          !< type of convection

    REAL(wp), INTENT(in)  :: &
      &  cos_mu0(kbdim),     & !< cos of zenith angle
      &  alb_vis_dir(kbdim), & !< surface albedo for visible range and direct light
      &  alb_nir_dir(kbdim), & !< surface albedo for NIR range and direct light
      &  alb_vis_dif(kbdim), & !< surface albedo for visible range and diffuse light
      &  alb_nir_dif(kbdim), & !< surface albedo for NIR range and diffuse light
      &  emis_rad(kbdim),    & !< longwave surface emissivity
      &  tk_sfc(kbdim),      & !< Surface temperature
      &  pp_hl(kbdim,klevp1),& !< pressure at half levels [Pa]
      &  pp_fl(kbdim,klev),  & !< Pressure at full levels [Pa]
      &  tk_fl(kbdim,klev),  & !< Temperature on full levels [K]
      &  qm_vap(kbdim,klev), & !< Water vapor mixing ratio
      &  qm_liq(kbdim,klev), & !< Liquid water mixing ratio
      &  qm_ice(kbdim,klev), & !< Ice water mixing ratio
      &  cld_frc(kbdim,klev),& !< Cloud fraction
      &  zaeq1(kbdim,klev) , & !< aerosol continental
      &  zaeq2(kbdim,klev) , & !< aerosol maritime
      &  zaeq3(kbdim,klev) , & !< aerosol mineral dust
      &  zaeq4(kbdim,klev) , & !< aerosol urban
      &  zaeq5(kbdim,klev) , & !< aerosol stratospheric background
      &  dust_tunefac(kbdim,jpband) !< LW tuning factor for dust aerosol


    REAL(wp), INTENT(in), POINTER, OPTIONAL  :: &
      &  reff_liq(:,:),      & !< Effective radius liquid phase [m]
      &  reff_frz(:,:)         !< Effective radius frozen phase [m]

    REAL(wp), INTENT(in),POINTER  :: &
      &  cdnc(:,:),          & !< Cloud drop number concentration
      &  zland(:),           & !< land-sea mask. (1. = land, 0. = sea/lakes)
      &  zglac(:)              !< fraction of land covered by glaciers

    ! output
    ! ------
    !
    REAL(wp), INTENT(inout) ::      &
      &  cld_cvr(kbdim),            & !< Cloud cover in a column
      &  flx_lw_net(kbdim,klevp1),  & !< Net longwave radiative flux [W/m**2]
      &  flx_uplw_sfc(kbdim),       & !< Srfc upward lw flux  [Wm2]
      &  trsol_net (kbdim,klevp1),  & !< Net solar transmissivity (=net radiative flux normalized by irradiance)
      &  trsol_up_toa(kbdim),       & !< TOA upward shortwave transmissivity (normalized upward flux)
      &  trsol_up_sfc(kbdim),       & !< Surface upward shortwave transmissivity (normalized upward flux)
      &  trsol_dn_sfc_diffus(kbdim),& !< Surface downward diffuse shortwave transmissivity (normalized diffuse downward flux)
      &  trsol_clr_sfc(kbdim),      & !< Surface net clear-sky solar transmissivity at surface
      &  lwflx_clr_sfc(kbdim),      & !< Longwave net clear-sky surface flux
      &  trsol_nir_sfc(kbdim),      & !< Surface transmissivity for near-IR part of solar radiation
      &  trsol_vis_sfc(kbdim),      & !< Surface transmissivity for visible part of solar radiation
      &  trsol_par_sfc(kbdim),      & !< Surface transmissivity for photosynthetically active part of solar radiation
      &  fr_nir_sfc_diffus(kbdim),  & !< Diffuse fraction of near-IR part of solar radiation
      &  fr_vis_sfc_diffus(kbdim),  & !< Diffuse fraction of visible part of solar radiation
      &  fr_par_sfc_diffus(kbdim)     !< Diffuse fraction of photosynthetically active part of solar radiation

    INTEGER  :: jk, jl

    REAL(wp) ::                     &
      &  pp_sfc(kbdim),             & !< surface pressure [Pa}
      &  tk_hl(kbdim,klevp1),       & !< Tempeature at half levels [Pa]
      &  xq_vap(kbdim,klev),        & !< Water vapor mixing ratio
      &  xq_liq(kbdim,klev),        & !< Liquid water mixing ratio
      &  xq_ice(kbdim,klev),        & !< Ice mixing ratio
      &  cld_frc_sec(kbdim,klev),   & !< secure cloud fraction [m2/m2]
      &  xm_co2(kbdim,klev),        & !< CO2 mixing ratio
      &  qm_o3(kbdim,klev),         & !< Ozone mixing ratio
      &  xm_o2(kbdim,klev),         & !< O2 mixing ratio
      &  xm_ch4(kbdim,klev),        & !< Methane mixing ratio
      &  xm_n2o(kbdim,klev),        & !< Nitrous Oxide mixing ratio
      &  xm_cfc11(kbdim,klev),      & !< CFC 11 mixing ratio
      &  xm_cfc12(kbdim,klev),      & !< CFC 12 mixing ratio
      &  flx_uplw_sfc_clr(kbdim),   & !< Srfc upward lw flux (clear sky) [Wm2]
      &  flx_upsw_sfc_clr(kbdim),   & !< Srfc upward sw flux (clear sky) [Wm2]
      &  flx_upsw_sfc(kbdim),       & !< Srfc upward sw flux (total) [Wm2]
      &  flx_upsw_toa(kbdim),       & !< TOA  upward sw flux (total) [Wm2]
      &  flx_dnsw_diff_sfc(kbdim),  & !< Srfc diffuse downward sw flux (total) [Wm2]
      &  flx_par_sfc(kbdim),        & !< Srfc downward PAR flux [Wm2]
      &  fr_vis_sfc(kbdim),         & !< Srfc downward visible flux fraction [Wm2]
      &  flx_sw_net(kbdim,klevp1),  & !< Net dwnwrd SW flux [Wm2]
      &  flx_sw_dn_sfc(kbdim),      & !< Downward SW flux [Wm2]
      &  flx_lw_net_clr(kbdim,klevp1),& !< Net dn LW flux (clear sky) [Wm2]
      &  flx_sw_net_clr(kbdim,klevp1)   !< Net dn SW flux (clear sky) [Wm2]
    !optional output: 3D flux output
    REAL(wp),  INTENT(inout), OPTIONAL :: &    
      &  flx_lw_dn(:,:),            & !< Downward LW flux (all-sky) [Wm2]
      &  flx_sw_dn(:,:),            & !< Downward SW flux (all-sky) [Wm2]
      &  flx_lw_up(:,:),            & !< Upward LW flux   (all-sky) [Wm2]
      &  flx_sw_up(:,:),            & !< Upward SW flux   (all-sky) [Wm2]
      &  flx_lw_dn_clr(:,:),        & !< Downward LW flux (clear sky) [Wm2]
      &  flx_sw_dn_clr(:,:),        & !< Downward SW flux (clear sky) [Wm2]
      &  flx_lw_up_clr(:,:),        & !< Upward LW flux   (clear sky) [Wm2]
      &  flx_sw_up_clr(:,:)           !< Upward SW flux   (clear sky) [Wm2]
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: pp_sfc,tk_hl,xq_vap,xq_liq,xq_ice,cld_frc_sec
!DIR$ ATTRIBUTES ALIGN : 64 :: xm_co2,xm_o2,xm_ch4,xm_cfc11,xm_cfc12
!DIR$ ATTRIBUTES ALIGN : 64 :: xm_n2o
!DIR$ ATTRIBUTES ALIGN : 64 :: flx_uplw_sfc_clr,flx_upsw_sfc_clr,flx_upsw_sfc
!DIR$ ATTRIBUTES ALIGN : 64 :: flx_upsw_toa,flx_dnsw_diff_sfc,flx_par_sfc
!DIR$ ATTRIBUTES ALIGN : 64 :: fr_vis_sfc
!DIR$ ATTRIBUTES ALIGN : 64 :: flx_sw_net,flx_lw_net_clr,flx_sw_net_clr
#endif

    IF (timers_level > 5) CALL timer_start(timer_radiation)

    !
    ! 1.1 p, T, q(vap,liq,ice) and clouds
    ! -----------------------------------
    !
    ! --- Pressure (surface and distance between half levels)
    !
    pp_sfc(1:jce)        = pp_hl(1:jce,klevp1)
    !
    ! --- temperature at half levels
    !
    DO jk=2,klev
      DO jl = 1, jce
        tk_hl(jl,jk) = (tk_fl(jl,jk-1)*pp_fl(jl,jk-1)*( pp_fl(jl,jk)          &
          &    - pp_hl(jl,jk) ) + tk_fl(jl,jk)*pp_fl(jl,jk)*( pp_hl(jl,jk)    &
          &    - pp_fl(jl,jk-1))) /(pp_hl(jl,jk)*(pp_fl(jl,jk) -pp_fl(jl,jk-1)))
      END DO
    END DO
    DO jl = 1, jce
 !     tk_hl(jl,klevp1) = tk_sfc(jl)
      tk_hl(jl,klevp1) = tk_fl(jl,klev) + (pp_hl(jl,klevp1) - pp_fl(jl,klev)) * &
         (tk_fl(jl,klev-1) - tk_fl(jl,klev))/(pp_fl(jl,klev-1) - pp_fl(jl,klev))
      tk_hl(jl,1)      = tk_fl(jl,1)+(pp_hl(jl,1)-pp_fl(jl,1))*(tk_fl(jl,1) - tk_hl(jl,2))  &
        &                / (pp_fl(jl,1)-pp_hl(jl,2))
    END DO
    !
    ! --- phases of water substance
    !
    xq_vap(1:jce,:) = MAX(qm_vap(1:jce,:),0.0_wp)
    xq_liq(1:jce,:) = MAX(qm_liq(1:jce,:),0.0_wp)       ! cloud liquid
    xq_ice(1:jce,:) = MAX(qm_ice(1:jce,:),0.0_wp)       ! cloud ice
    !
    ! --- cloud cover
    !
    DO jk = 1, klev
      DO jl = 1, jce

        ! no cloud water -> no cloud fraction
        IF (xq_liq(jl,jk) > 0.0_wp .OR. xq_ice(jl,jk) > 0.0_wp) THEN
          cld_frc_sec(jl,jk) = cld_frc(jl,jk)
        ELSE
          cld_frc_sec(jl,jk) = 0.0_wp
        END IF

        ! cloud fraction <= 100% !
        cld_frc_sec(jl,jk) = MIN( cld_frc_sec(jl,jk), 1.0_wp )

      END DO
    END DO

    cld_cvr(1:jce) = 1.0_wp - cld_frc_sec(1:jce,1)
    DO jk = 2, klev
      cld_cvr(1:jce) = cld_cvr(1:jce)                                                &
        &              *(1.0_wp-MAX(cld_frc_sec(1:jce,jk),cld_frc_sec(1:jce,jk-1)))  &
        &              /(1.0_wp-MIN(cld_frc_sec(1:jce,jk-1),1.0_wp-EPSILON(1.0_wp)))
    END DO
    cld_cvr(1:jce) = 1.0_wp-cld_cvr(1:jce)
    !
    ! 1.2 Non-water tracers
    ! ---------------------
    !
    ! --- gas profiles in [ppm]
    !
    xm_co2   (1:jce,:) = gas_profile(jce, klev, irad_co2,        &
      &                              mmr_gas = mmr_co2,          &
      &                              gas_scenario = ghg_co2mmr   )
    xm_ch4   (1:jce,:) = gas_profile(jce, klev, irad_ch4,        &
      &                              mmr_gas = mmr_ch4,          &
      &                              gas_scenario = ghg_ch4mmr,  &
      &                              pressure = pp_fl,           &
      &                              xp = vpp_ch4                )
    xm_n2o   (1:jce,:) = gas_profile(jce, klev, irad_n2o,        &
      &                              mmr_gas = mmr_n2o,          &
      &                              gas_scenario = ghg_n2ommr,  &
      &                              pressure = pp_fl,           &
      &                              xp = vpp_n2o                )
    xm_o2    (1:jce,:) = gas_profile(jce, klev, irad_o2,         &
      &                              mmr_gas = mmr_o2            )
#ifdef __SX__
    xm_cfc11 (1:jce,:) = gas_profile(jce, klev, irad_cfc11,        &
      &                              mmr_gas = REAL(mmr_cfc11,wp), &
      &                              gas_scenario = ghg_cfcmmr(1)  )
    xm_cfc12 (1:jce,:) = gas_profile(jce, klev, irad_cfc12,        &
      &                              mmr_gas = REAL(mmr_cfc12,wp), &
      &                              gas_scenario = ghg_cfcmmr(2)  )
#else
    xm_cfc11 (1:jce,:) = gas_profile(jce, klev, irad_cfc11,        &
      &                              mmr_gas = mmr_cfc11 ,         &
      &                              gas_scenario = ghg_cfcmmr(1)  ) 
    xm_cfc12 (1:jce,:) = gas_profile(jce, klev, irad_cfc12,        &
      &                              mmr_gas = mmr_cfc12,          &
      &                              gas_scenario = ghg_cfcmmr(2)  )
#endif
    !

    ! 2.0 Call interface to radiation solver
    ! --------------------------------------
    !
    CALL rrtm_interface(                                                    &
      ! input
      & current_date                                                       ,&
      & jg              ,jb              ,icpl_reff                        ,&
      & jcs             ,jce             ,kbdim           ,klev            ,&
      & ktype           ,zland           ,zglac                            ,&
      & cos_mu0                                                            ,&
      & alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif     ,&
      & emis_rad                                                           ,&
      & pp_fl           ,pp_hl           ,pp_sfc          ,tk_fl           ,&
      & tk_hl           ,tk_sfc          ,xq_vap                           ,&
      & xq_liq          ,xq_ice                                            ,&
      & cdnc            ,reff_liq        ,reff_frz                         ,&
      & cld_frc_sec                                                        ,&
      & qm_o3           ,xm_co2          ,xm_ch4                           ,&
      & xm_n2o          ,xm_cfc11        ,xm_cfc12        ,xm_o2           ,&
      & zaeq1,zaeq2,zaeq3,zaeq4,zaeq5                                      ,&
      ! output
      & flx_lw_net      ,flx_sw_net      ,flx_lw_net_clr  ,flx_sw_net_clr  ,&
      & flx_uplw_sfc    ,flx_upsw_sfc    ,flx_uplw_sfc_clr,flx_upsw_sfc_clr,&
      ! optional output: 3D flux output
      & flx_lw_dn       ,flx_sw_dn       ,flx_lw_up       ,flx_sw_up       ,&
      & flx_lw_dn_clr   ,flx_sw_dn_clr   ,flx_lw_up_clr   ,flx_sw_up_clr   ,&
      ! optional arguments
      & dust_tunefac=dust_tunefac ,flx_dnsw_diff_sfc=flx_dnsw_diff_sfc     ,&
      & flx_upsw_toa=flx_upsw_toa ,flx_dnpar_sfc=flx_par_sfc               ,&
      & vis_frc_sfc=fr_vis_sfc, nir_dff_frc_sfc=fr_nir_sfc_diffus          ,&
      & vis_dff_frc_sfc=fr_vis_sfc_diffus, par_dff_frc_sfc=fr_par_sfc_diffus)


    !
    ! 3.0 Diagnostics
    ! ---------------
    !
    ! --- Diagnose transmissivities from solar fluxes

    trsol_net (1:jce,1:klevp1) = flx_sw_net (1:jce,1:klevp1)/SPREAD(cos_mu0(1:jce)*tsi_radt,2,klevp1)
    trsol_clr_sfc(1:jce)       = flx_sw_net_clr(1:jce,klevp1)/(cos_mu0(1:jce)*tsi_radt)
    trsol_up_sfc(1:jce)        = flx_upsw_sfc  (1:jce)/(cos_mu0(1:jce)*tsi_radt)
    trsol_up_toa(1:jce)        = flx_upsw_toa  (1:jce)/(cos_mu0(1:jce)*tsi_radt)
    trsol_dn_sfc_diffus(1:jce) = flx_dnsw_diff_sfc(1:jce)/(cos_mu0(1:jce)*tsi_radt)

    flx_sw_dn_sfc(1:jce)       = flx_sw_net(1:jce,klevp1) + flx_upsw_sfc(1:jce)

    trsol_nir_sfc(1:jce)       = (1._wp - fr_vis_sfc(1:jce)) * flx_sw_dn_sfc(1:jce)/(cos_mu0(1:jce)*tsi_radt)
    trsol_vis_sfc(1:jce)       = fr_vis_sfc(1:jce) * flx_sw_dn_sfc(1:jce)/(cos_mu0(1:jce)*tsi_radt)
    trsol_par_sfc(1:jce)       = flx_par_sfc(1:jce)/(cos_mu0(1:jce)*tsi_radt)

    lwflx_clr_sfc(1:jce)       = flx_lw_net_clr(1:jce,klevp1)

    IF (timers_level > 5) CALL timer_stop(timer_radiation)

  END SUBROUTINE radiation_nwp

  !---------------------------------------------------------------------------
  !>
  !! GAS_PROFILE:  Determines Gas distributions based on case specification
  !!
  !! Description: This routine calculates the gas distributions for one of
  !! five cases:  (0) no gas present; (1) prognostic gas; (2) specified
  !! mixing ratio; (3) mixing ratio decaying with height given profile;
  !! (4) scenario run with different mixing ratio.
  !
  FUNCTION gas_profile (jce, klev, igas, mmr_gas, gas_scenario, mmr_gas_v, &
    &                   gas_scenario_v, gas_val, xp, pressure)

    INTEGER,  INTENT (in)           :: jce, klev, igas
    REAL(wp), INTENT (in), OPTIONAL :: mmr_gas, gas_scenario
    REAL(wp), INTENT (in), OPTIONAL :: pressure(:,:), xp(3)
    REAL(wp), INTENT (in), OPTIONAL :: mmr_gas_v(:,:)
    REAL(wp), INTENT (in), OPTIONAL :: gas_scenario_v(:,:)
    REAL(wp), INTENT (in), OPTIONAL :: gas_val(:,:)

    REAL(wp) :: gas_profile(jce,klev)
    REAL(wp) :: zx_d, zx_m
    LOGICAL  :: gas_initialized

    gas_initialized = .FALSE.
    SELECT CASE (igas)
    CASE (0)
      gas_profile(1:jce,:) = 0.0_wp
      gas_initialized = .TRUE.
    CASE (1)
      IF (PRESENT(gas_val)) THEN
        gas_profile(1:jce,:) = MAX(gas_val(1:jce,:), EPSILON(1.0_wp))
        gas_initialized = .TRUE.
      END IF
    CASE (2)
      IF (PRESENT(mmr_gas)) THEN
        gas_profile(1:jce,:) = mmr_gas
        gas_initialized = .TRUE.
      ELSE IF (PRESENT(mmr_gas_v)) THEN
        gas_profile(1:jce,:) = mmr_gas_v(1:jce,:)
        gas_initialized = .TRUE.
      END IF
    CASE (3)
      IF (PRESENT(mmr_gas) .AND. PRESENT(xp) .AND. PRESENT(pressure)) THEN
        zx_m = (mmr_gas+xp(1)*mmr_gas)*0.5_wp
        zx_d = (mmr_gas-xp(1)*mmr_gas)*0.5_wp
        gas_profile(1:jce,:)=(1._wp-(zx_d/zx_m)*TANH(LOG(pressure(1:jce,:)   &
          &                     /xp(2)) /xp(3))) * zx_m
        gas_initialized = .TRUE.
      END IF
    CASE (4)
      IF (PRESENT(gas_scenario)) THEN
        ! TODO: (Hauke Schmidt)
        ! If the respective parameters are present, a vertical
        ! profile is calculated as in option (3). This allows a seamless
        ! continuation of preindustrial control with scenarios. The treatment here is
        ! inconsistent with having two different options for the constant
        ! concentration cases (2 without and 3 with profile). However, instead
        ! of adding a fifth option, it seems more advisable to clean up the
        ! complete handling of radiation switches, later.
        IF (PRESENT(xp) .AND. PRESENT(pressure)) THEN
          zx_m = (gas_scenario+xp(1)*gas_scenario)*0.5_wp
          zx_d = (gas_scenario-xp(1)*gas_scenario)*0.5_wp
          gas_profile(1:jce,:)=(1._wp-(zx_d/zx_m)*TANH(LOG(pressure(1:jce,:)   &
            &                     /xp(2)) /xp(3))) * zx_m
        ELSE
          gas_profile(1:jce,:) = gas_scenario
        END IF
        gas_initialized = .TRUE.
      ELSE IF (PRESENT(gas_scenario_v)) THEN
        gas_profile(1:jce,:) = gas_scenario_v(1:jce,:)
        gas_initialized = .TRUE.
      END IF
    END SELECT
    IF (.NOT. gas_initialized) &
      CALL finish('radiation','gas_profile options not supported')

  END FUNCTION gas_profile
  !-----------------------------------------------------------------------------
  !>
  !! Arranges input and calls rrtm sw and lw routines
  !!
  !! Remarks
  !!   Because the RRTM indexes vertical levels differently than ECHAM a chief
  !!   function of this routine is to reorder the input in the vertical.  In
  !!   addition some cloud physical properties are prescribed, which are
  !!   required to derive cloud optical properties
  !!
  !! The gases are passed into RRTM via two multi-constituent arrays:
  !!   zwkl and wx_r. zwkl has JPINPX species and  wx_r has JPXSEC species
  !!   The species are identifed as follows.
  !!     ZWKL [#/cm2]          WX_R [#/cm2]
  !!    index = 1 => H20     index = 1 => n/a
  !!    index = 2 => CO2     index = 2 => CFC11
  !!    index = 3 =>  O3     index = 3 => CFC12
  !!    index = 4 => N2O     index = 4 => n/a
  !!    index = 5 => n/a
  !!    index = 6 => CH4
  !!    index = 7 => O2
  !


  SUBROUTINE rrtm_interface(                                              &
    ! input
    & current_date                                                       ,&
    & jg              ,jb              ,icpl_reff                        ,&
    & jcs             ,jce             ,kbdim           ,klev            ,&
    & ktype           ,zland           ,zglac                            ,&
    & pmu0                                                               ,&
    & alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif     ,&
    & emis_rad                                                           ,&
    & pp_fl           ,pp_hl           ,pp_sfc          ,tk_fl           ,&
    & tk_hl           ,tk_sfc          ,xm_vap                           ,&
    & xm_liq          ,xm_ice                                            ,&
    & cdnc            ,reff_liq        ,reff_frz                         ,&
    & cld_frc                                                            ,&
    & xm_o3           ,xm_co2          ,xm_ch4                           ,&
    & xm_n2o          ,xm_cfc11        ,xm_cfc12        ,xm_o2           ,&
    & zaeq1, zaeq2, zaeq3, zaeq4, zaeq5                                  ,&
    ! output
    & flx_lw_net      ,flx_sw_net      ,flx_lw_net_clr  ,flx_sw_net_clr  ,&
    & flx_uplw_sfc    ,flx_upsw_sfc    ,flx_uplw_sfc_clr,flx_upsw_sfc_clr,&
    ! optional output: 3D flux output
    & flx_lw_dn       ,flx_sw_dn       ,flx_lw_up       ,flx_sw_up       ,&
    & flx_lw_dn_clr   ,flx_sw_dn_clr   ,flx_lw_up_clr   ,flx_sw_up_clr   ,&
    ! optional input
    & dust_tunefac                                                       ,&
    ! optional output
    & flx_dnsw_diff_sfc,flx_upsw_toa   ,flx_dnpar_sfc                    ,&
    & vis_frc_sfc      ,nir_dff_frc_sfc,vis_dff_frc_sfc ,par_dff_frc_sfc  )

    TYPE(datetime), POINTER, INTENT(in) :: current_date
    
    INTEGER,INTENT(in)  ::                &
      &  jg,                              & !< domain index
      &  jb,                              & !< block index
      &  icpl_reff,                       & !< in option for radiation reff coupling
      &  jcs,                             & !< number of skipped columns
      &  jce,                             & !< number of columns
      &  kbdim,                           & !< first dimension of 2-d arrays
      &  klev                               !< number of levels

    INTEGER,INTENT(in)  ::                &
      &  ktype(kbdim)                       !< type of convection

    REAL(wp),INTENT(in) ::                &
      &  pmu0(kbdim),                     & !< mu0 for solar zenith angle
      &  alb_vis_dir(kbdim),              & !< surface albedo for vis range and dir light
      &  alb_nir_dir(kbdim),              & !< surface albedo for NIR range and dir light
      &  alb_vis_dif(kbdim),              & !< surface albedo for vis range and dif light
      &  alb_nir_dif(kbdim),              & !< surface albedo for NIR range and dif light
      &  emis_rad(kbdim),                 & !< longwave surface emissivity
      &  pp_fl(kbdim,klev),               & !< full level pressure in Pa
      &  pp_hl(kbdim,klev+1),             & !< half level pressure in Pa
      &  pp_sfc(kbdim),                   & !< surface pressure in Pa
      &  tk_fl(kbdim,klev),               & !< full level temperature in K
      &  tk_hl(kbdim,klev+1),             & !< half level temperature in K
      &  tk_sfc(kbdim),                   & !< surface temperature in K
      &  xm_vap(kbdim,klev),              & !< specific humidity in g/g
      &  xm_liq(kbdim,klev),              & !< specific liquid water content
      &  xm_ice(kbdim,klev),              & !< specific ice content in g/g
      &  cld_frc(kbdim,klev),             & !< fractional cloud cover
      &  xm_o3(kbdim,klev),               & !< o3 mass mixing ratio
      &  xm_co2(kbdim,klev),              & !< co2 mass mixing ratio
      &  xm_ch4(kbdim,klev),              & !< ch4 mass mixing ratio
      &  xm_n2o(kbdim,klev),              & !< n2o mass mixing ratio
      &  xm_cfc11(kbdim,klev),            & !< cfc 11 mass mixing ratio
      &  xm_cfc12(kbdim,klev),            & !< cfc 12 mass mixing ratio
      &  xm_o2(kbdim,klev),               & !< o2  mass mixing ratio
      &  zaeq1(kbdim,klev),               & !< aerosol continental
      &  zaeq2(kbdim,klev),               & !< aerosol maritime
      &  zaeq3(kbdim,klev),               & !< aerosol urban
      &  zaeq4(kbdim,klev),               & !< aerosol volcano ashes
      &  zaeq5(kbdim,klev)                  !< aerosol stratospheric background

    REAL(wp),INTENT(in),POINTER  ::                &
      &  cdnc(:,:),                       & !< cloud nuclei concentration
      &  zland(:),                        & !< land-sea mask. (1. = land, 0. = sea/lakes)
      &  zglac(:)                           !< fraction of land covered by glaciers


    REAL(wp),INTENT(in), POINTER, OPTIONAL :: &
      &  reff_liq(:,:),              & !< effective radius liquid phase in m
      &  reff_frz(:,:)                 !< effective radius frozen phase in m

    REAL(wp), INTENT(out) ::              &
      &  flx_lw_net(kbdim,klev+1),        & !< net downward LW flux profile,
      &  flx_sw_net(kbdim,klev+1),        & !< net downward SW flux profile,
      &  flx_lw_net_clr(kbdim,klev+1),    & !< clrsky downward LW flux profile,
      &  flx_sw_net_clr(kbdim,klev+1),    & !< clrsky downward SW flux profile,
      &  flx_uplw_sfc(kbdim),             & !< sfc LW upward flux,
      &  flx_upsw_sfc(kbdim),             & !< sfc SW upward flux,
      &  flx_uplw_sfc_clr(kbdim),         & !< clrsky sfc LW upward flux,
      &  flx_upsw_sfc_clr(kbdim)            !< clrsky sfc SW upward flux,

    REAL(wp), INTENT(in),  OPTIONAL ::    dust_tunefac(kbdim,jpband) ! LW absorption tuning factor for dust

    REAL(wp), INTENT(out), OPTIONAL ::    &
      &  flx_dnsw_diff_sfc(kbdim),        & !< sfc SW diffuse downward flux,
      &  flx_upsw_toa(kbdim),             & !< TOA SW upward flux,
      &  flx_dnpar_sfc(kbdim),            & !< PAR downward sfc flux
      &  vis_frc_sfc(kbdim),              & !< Visible fraction of net surface SW radiation
      &  nir_dff_frc_sfc(kbdim),          & !< Diffuse fraction of downward surface near-infrared radiation at surface
      &  vis_dff_frc_sfc(kbdim),          & !< Diffuse fraction of downward surface visible radiation at surface
      &  par_dff_frc_sfc(kbdim),          & !< Diffuse fraction of downward surface PAR
      ! 3D flux output
      &  flx_lw_dn(:,:),                  & !< Downward LW flux profile, all-sky
      &  flx_sw_dn(:,:),                  & !< Downward SW flux profile, all-sky
      &  flx_lw_up(:,:),                  & !< Upward LW flux profile, all-sky
      &  flx_sw_up(:,:),                  & !< Upward SW flux profile, all-sky
      &  flx_lw_dn_clr(:,:),              & !< Downward LW flux profile, clear sky
      &  flx_sw_dn_clr(:,:),              & !< Downward SW flux profile, clear sky
      &  flx_lw_up_clr(:,:),              & !< Upward LW flux profile, clear sky
      &  flx_sw_up_clr(:,:)                 !< Upward SW flux profile, clear sky

    INTEGER  :: jk, jl, jp, jkb, jspec,   & !< loop indicies
      &  icldlyr(kbdim,klev)                !< index for clear or cloudy

    REAL(wp) ::                           &
      &  zsemiss(kbdim,jpband),           & !< LW surface emissivity by band
      &  pm_sfc(kbdim),                   & !< surface pressure in hPa
      &  amm,                             & !< molecular weight of moist air
      &  delta,                           & !< pressure thickness
      &  zscratch                           !< scratch array

    REAL(wp) :: z_sum_aea, z_sum_aes !help variables for aerosol

    !
    ! --- vertically reversed _vr variables
    !
    REAL(wp) ::                           &
      &  col_dry_vr(kbdim,klev),          & !< number of molecules/cm2 of
      &  pm_fl_vr(kbdim,klev),            & !< full level pressure [hPa]
      &  tk_fl_vr(kbdim,klev),            & !< full level temperature [K]
      &  tk_hl_vr(kbdim,klev+1),          & !< half level temperature [K]
      &  cld_frc_vr(kbdim,klev),          & !< secure cloud fraction
      &  ziwgkg_vr(kbdim,klev),           & !< specific ice water content
      &  ziwc_vr(kbdim,klev),             & !< ice water content per volume
      &  ziwp_vr(kbdim,klev),             & !< ice water path in g/m2
      &  zlwgkg_vr(kbdim,klev),           & !< specific liquid water content
      &  zlwp_vr(kbdim,klev),             & !< liquid water path in g/m2
      &  zlwc_vr(kbdim,klev),             & !< liquid water content per
      &  wx_vr(kbdim,jpxsec,klev),        & !< number of molecules/cm2 of
      &  cld_tau_lw_vr(kbdim,klev,jpband),& !< LW optical thickness of clouds
      &  cld_tau_sw_vr(kbdim,jpsw,klev),  & !< extincion
      &  cld_cg_sw_vr(kbdim,jpsw,klev),   & !< asymmetry factor
      &  cld_piz_sw_vr(kbdim,jpsw,klev),  & !< single scattering albedo
      &  aer_tau_lw_vr(kbdim,klev,jpband),& !< LW optical thickness of aerosols
      &  aer_tau_sw_vr(kbdim,klev,jpsw),  & !< aerosol optical thickness
      &  aer_cg_sw_vr(kbdim,klev,jpsw),   & !< aerosol asymmetry factor
      &  aer_piz_sw_vr(kbdim,klev,jpsw),  & !< aerosol single scattering albedo
      &  flx_uplw_vr(kbdim,klev+1),       & !< upward flux, total sky
      &  flx_uplw_clr_vr(kbdim,klev+1),   & !< upward flux, clear sky
      &  flx_dnlw_vr(kbdim,klev+1),       & !< downward flux, total sky
      &  flx_dnlw_clr_vr(kbdim,klev+1),   & !< downward flux, clear sky
      &  flx_upsw(kbdim,klev+1),          & !< upward flux total sky
      &  flx_upsw_clr(kbdim,klev+1),      & !< upward flux clear sky
      &  flx_dnsw(kbdim,klev+1),          & !< downward flux total sky
      &  flx_dnsw_clr(kbdim,klev+1),      & !< downward flux clear sky
      ! Alberto: Using pointers does not work here because the routine is called inside an OMP
      ! region and produce a race condition when memory is allocated.
      &  cdnc_vr(kbdim,klev)     ,        & !< cloud nuclei concentration
      &  reff_liq_vr(kbdim,klev) ,        & !< effective radius liquid phase [m]
      &  reff_frz_vr(kbdim,klev)            !< effective radius frozen phase [m]



    REAL(wp), TARGET ::                   &
      &  wkl_vr(kbdim,jpinpx,klev)        !< number of molecules/cm2 of
    REAL(wp) :: tune_dust(kbdim,jpband)  ! local variable for LW absorption tuning of dust

    CHARACTER(LEN=3)     :: c_irad_aero


    LOGICAL:: l_coupled_reff          ! Use the effective radius from microphysics
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: icldlyr,zsemiss,pm_sfc
!DIR$ ATTRIBUTES ALIGN : 64 :: col_dry_vr,pm_fl_vr,tk_fl_vr,tk_hl_vr
!DIR$ ATTRIBUTES ALIGN : 64 :: cdnc_vr,cld_frc_vr,ziwgkg_vr,ziwc_vr,ziwp_vr
!DIR$ ATTRIBUTES ALIGN : 64 :: reff_liq_vr,reff_frz_vr
!DIR$ ATTRIBUTES ALIGN : 64 :: zlwgkg_vr,zlwp_vr,zlwc_vr,wkl_vr,wx_vr
!DIR$ ATTRIBUTES ALIGN : 64 :: cld_tau_lw_vr,cld_tau_sw_vr,cld_cg_sw_vr,cld_piz_sw_vr
!DIR$ ATTRIBUTES ALIGN : 64 :: aer_tau_lw_vr,aer_tau_sw_vr,aer_cg_sw_vr,aer_piz_sw_vr
!DIR$ ATTRIBUTES ALIGN : 64 :: flx_uplw_vr,flx_uplw_clr_vr,flx_dnlw_vr
!DIR$ ATTRIBUTES ALIGN : 64 :: flx_dnlw_clr_vr,flx_upsw,flx_upsw_clr,flx_dnsw
!DIR$ ATTRIBUTES ALIGN : 64 :: flx_dnsw_clr
!DIR$ ATTRIBUTES ALIGN : 64 :: tune_dust
#endif


    ! Initialize output variables
    flx_lw_net(:,:)     = 0._wp
    flx_lw_net_clr(:,:) = 0._wp
    flx_sw_net(:,:)     = 0._wp
    flx_sw_net_clr(:,:) = 0._wp
    flx_uplw_sfc(:)     = 0._wp
    flx_uplw_sfc_clr(:) = 0._wp
    flx_upsw_sfc(:)     = 0._wp
    flx_upsw_sfc_clr(:) = 0._wp
    
    l_coupled_reff = icpl_reff > 0

    IF (atm_phy_nwp_config(jg)%l_3d_rad_fluxes) THEN
      IF (PRESENT(flx_lw_dn))     flx_lw_dn(:,:)      = 0._wp
      IF (PRESENT(flx_sw_dn))     flx_sw_dn(:,:)      = 0._wp
      IF (PRESENT(flx_lw_up))     flx_lw_up(:,:)      = 0._wp
      IF (PRESENT(flx_sw_up))     flx_sw_up(:,:)      = 0._wp
      IF (PRESENT(flx_lw_dn_clr)) flx_lw_dn_clr(:,:)  = 0._wp
      IF (PRESENT(flx_sw_dn_clr)) flx_sw_dn_clr(:,:)  = 0._wp
      IF (PRESENT(flx_lw_up_clr)) flx_lw_up_clr(:,:)  = 0._wp
      IF (PRESENT(flx_sw_up_clr)) flx_sw_up_clr(:,:)  = 0._wp
    END IF

    IF (PRESENT(flx_dnsw_diff_sfc)) flx_dnsw_diff_sfc(:) = 0._wp
    IF (PRESENT(flx_upsw_toa))      flx_upsw_toa(:)      = 0._wp
    IF (PRESENT(flx_dnpar_sfc))     flx_dnpar_sfc(:)     = 0._wp
    IF (PRESENT(vis_frc_sfc))       vis_frc_sfc(:)       = 0._wp
    IF (PRESENT(nir_dff_frc_sfc))   nir_dff_frc_sfc(:)   = 0._wp
    IF (PRESENT(vis_dff_frc_sfc))   vis_dff_frc_sfc(:)   = 0._wp
    IF (PRESENT(par_dff_frc_sfc))   par_dff_frc_sfc(:)   = 0._wp
    
    !
    ! 1.0 Constituent properties
    !--------------------------------

    IF (timers_level > 7) CALL timer_start(timer_rrtm_prep)

    !
    ! --- control for infintesimal cloud fractions
    !
    DO jk = 1, klev
      !
      jkb = klev+1-jk
      cld_frc_vr(1:jce,jk)  = cld_frc(1:jce,jkb)

      DO jl=1,jce
        IF (cld_frc_vr(jl,jk) > 2.0_wp*EPSILON(1.0_wp)) THEN
          ! only clouds > 2 epsilon are made visible to radiation
          icldlyr  (jl,jk) = 1
          ziwgkg_vr(jl,jk) = xm_ice(jl,jkb)*1000.0_wp/cld_frc_vr(jl,jk)
          zlwgkg_vr(jl,jk) = xm_liq(jl,jkb)*1000.0_wp/cld_frc_vr(jl,jk)
        ELSE
          ! clouds <= 2 epsilon are ade invisble to radiation
          icldlyr  (jl,jk) = 0
          ziwgkg_vr(jl,jk) = 0.0_wp
          zlwgkg_vr(jl,jk) = 0.0_wp
        ENDIF
      END DO
    END DO
    !
    ! --- main constituent reordering
    !
    DO jl = 1, jce
      tk_hl_vr(jl,klev+1) = tk_hl(jl,1)
      pm_sfc(jl)          = 0.01_wp*pp_sfc(jl)
    END DO

    DO jk = 1, klev
      jkb = klev+1-jk
      ! initialization
      wkl_vr(:,:,jk) = 0.0_wp
      wx_vr (:,:,jk) = 0.0_wp
      DO jl = 1, jce
        delta = pp_hl(jl,jkb+1)-pp_hl(jl,jkb)
        !
        ! --- thermodynamic arrays
        !
        pm_fl_vr(jl,jk) = 0.01_wp*pp_fl(jl,jkb)
        tk_hl_vr(jl,jk) = tk_hl(jl,jkb+1)
        tk_fl_vr(jl,jk) = tk_fl(jl,jkb)
        !
        ! --- cloud properties
        !
        zscratch       = pp_fl(jl,jkb)/tk_fl(jl,jkb)
        ziwc_vr(jl,jk) = ziwgkg_vr(jl,jk)*zscratch/rd
        ziwp_vr(jl,jk) = ziwgkg_vr(jl,jk)*delta/grav
        zlwc_vr(jl,jk) = zlwgkg_vr(jl,jk)*zscratch/rd
        zlwp_vr(jl,jk) = zlwgkg_vr(jl,jk)*delta/grav
        !
        ! --- radiatively active gases
        !
        wkl_vr(jl,1,jk)   = xm_vap(jl,jkb)*amd/amw
        wkl_vr(jl,2,jk)   = xm_co2(jl,jkb)*amd/amco2
        wkl_vr(jl,3,jk)   = xm_o3(jl,jkb) *amd/amo3
        wkl_vr(jl,4,jk)   = xm_n2o(jl,jkb)*amd/amn2o
        wkl_vr(jl,6,jk)   = xm_ch4(jl,jkb)*amd/amch4
        wkl_vr(jl,7,jk)   = xm_o2 (jl,jkb)*amd/amo2
        amm               = (1.0_wp-wkl_vr(jl,1,jk))*amd + wkl_vr(jl,1,jk)*amw
        col_dry_vr(jl,jk) = (0.01_wp*delta)*10.0_wp*avo/grav/amm / (1.0_wp+wkl_vr(jl,1,jk))
        !
        ! --- alternate treatment for cfcs
        !
        wx_vr(jl,2,jk) = col_dry_vr(jl,jk)*xm_cfc11(jl,jkb)*amd/amc11*1.e-20_wp
        wx_vr(jl,3,jk) = col_dry_vr(jl,jk)*xm_cfc12(jl,jkb)*amd/amc12*1.e-20_wp
      END DO
        !
        ! --- extra cloud properties for rad coupling
        ! 
      IF ( l_coupled_reff ) THEN
        DO jl = 1, jce
          reff_liq_vr(jl,jk) = reff_liq(jl,jkb)
          reff_frz_vr(jl,jk) = reff_frz(jl,jkb)      
        END DO
      ELSE
        DO jl = 1, jce
          cdnc_vr(jl,jk) = cdnc(jl,jkb)*1.e-6_wp
        END DO
      END IF
    END DO
    DO jp = 1, 7
      wkl_vr(1:jce,jp,:)=col_dry_vr(1:jce,:)*wkl_vr(1:jce,jp,:)
    END DO
    !
    ! 2.0 Surface Properties
    ! --------------------------------
    zsemiss(1:jce,:) = SPREAD(emis_rad(1:jce),2,jpband)
    !
    ! 3.0 Particulate Optical Properties
    ! --------------------------------

    SELECT CASE (irad_aero)
    CASE (iRadAeroNone,iRadAeroConst)
      aer_tau_lw_vr(:,:,:) = 0.0_wp
      aer_tau_sw_vr(:,:,:) = 0.0_wp
      aer_piz_sw_vr(:,:,:) = 1.0_wp
      aer_cg_sw_vr(:,:,:)  = 0.0_wp
    CASE (iRadAeroTegen)
      IF (PRESENT(dust_tunefac)) THEN
        tune_dust(1:jce,1:jpband) = dust_tunefac(1:jce,1:jpband)
      ELSE
        tune_dust(1:jce,1:jpband) = 1._wp
      ENDIF
      DO jspec=1,jpband
        DO jk=1,klev
          jkb = klev+1-jk
          DO jl = 1,jce
            ! LW opt thickness of aerosols
            aer_tau_lw_vr(jl,jk,jspec) =  zaeq1(jl,jkb) * zaea_rrtm(jspec,1) &
              &                         + zaeq2(jl,jkb) * zaea_rrtm(jspec,2) &
              &   + tune_dust(jl,jspec) * zaeq3(jl,jkb) * zaea_rrtm(jspec,3) &
              &                         + zaeq4(jl,jkb) * zaea_rrtm(jspec,4) &
              &                         + zaeq5(jl,jkb) * zaea_rrtm(jspec,5)
          ENDDO
        ENDDO
      ENDDO
      DO jspec=1+jpband,jpband+jpsw
        DO jk=1,klev
          jkb = klev+1-jk
          DO jl = 1,jce

            z_sum_aea = zaeq1(jl,jkb) * zaea_rrtm(jspec,1) &
              &       + zaeq2(jl,jkb) * zaea_rrtm(jspec,2) &
              &       + zaeq3(jl,jkb) * zaea_rrtm(jspec,3) &
              &       + zaeq4(jl,jkb) * zaea_rrtm(jspec,4) &
              &       + zaeq5(jl,jkb) * zaea_rrtm(jspec,5)

            z_sum_aes = zaeq1(jl,jkb) * zaes_rrtm(jspec,1) &
              &       + zaeq2(jl,jkb) * zaes_rrtm(jspec,2) &
              &       + zaeq3(jl,jkb) * zaes_rrtm(jspec,3) &
              &       + zaeq4(jl,jkb) * zaes_rrtm(jspec,4) &
              &       + zaeq5(jl,jkb) * zaes_rrtm(jspec,5)

            ! sw aerosol optical thickness
            aer_tau_sw_vr(jl,jk,jspec-jpband) = z_sum_aea + z_sum_aes

            ! sw aerosol single scattering albedo
            aer_piz_sw_vr(jl,jk,jspec-jpband) = z_sum_aes / MAX(dbl_eps, z_sum_aea + z_sum_aes)

            ! sw aerosol asymmetry factor
            aer_cg_sw_vr(jl,jk,jspec-jpband) =                                  &
              & (   zaeq1(jl,jkb) * zaes_rrtm(jspec,1) * zaeg_rrtm(jspec,1)   &
              &   + zaeq2(jl,jkb) * zaes_rrtm(jspec,2) * zaeg_rrtm(jspec,2)   &
              &   + zaeq3(jl,jkb) * zaes_rrtm(jspec,3) * zaeg_rrtm(jspec,3)   &
              &   + zaeq4(jl,jkb) * zaes_rrtm(jspec,4) * zaeg_rrtm(jspec,4)   &
              &   + zaeq5(jl,jkb) * zaes_rrtm(jspec,5) * zaeg_rrtm(jspec,5) ) / MAX(dbl_eps, z_sum_aes)
          ENDDO
        ENDDO
      ENDDO
#ifdef __ICON_ART
    CASE (iRadAeroART)
      CALL art_rad_aero_interface(zaeq1,zaeq2,zaeq3,zaeq4,zaeq5, &
        &                         zaea_rrtm,zaes_rrtm,zaeg_rrtm, &
        &                         jg,jb,1,klev,1,jce,jpband,jpsw,&
        &                         aer_tau_lw_vr,                 &
        &                         aer_tau_sw_vr,                 &
        &                         aer_piz_sw_vr,                 &
        &                         aer_cg_sw_vr)
#endif
    CASE (iRadAeroConstKinne,iRadAeroKinne)
       ! this is for rrtm radiation in the NWP part, we do not introduce
       ! the simple plumes here
       CALL set_bc_aeropt_kinne( current_date                        ,&
        & jg                                                        ,&
        & 1, jce           ,kbdim                 ,klev             ,&
        & jb               ,jpsw                  ,jpband           ,&
        & p_nh_state(jg)% metrics% z_mc(:,:,jb)                     ,&
        & p_nh_state(jg)% metrics% ddqz_z_full(:,:,jb)              ,&
        & aer_tau_sw_vr    ,aer_piz_sw_vr         , aer_cg_sw_vr    ,&
        & aer_tau_lw_vr                                              )
    CASE (iRadAeroVolc)
      ! this is for rrtm radiation in the NWP part, we do not introduce
      ! the simple plumes here
      ! set zero aerosol before adding CMIP6 volcanic aerosols
      aer_tau_lw_vr(:,:,:) = 0.0_wp
      aer_tau_sw_vr(:,:,:) = 0.0_wp
      aer_piz_sw_vr(:,:,:) = 1.0_wp
      aer_cg_sw_vr(:,:,:)  = 0.0_wp
      CALL add_bc_aeropt_cmip6_volc( current_date ,jg,            1 ,&
        & jce              ,kbdim                 ,klev             ,&
        & jb               ,jpsw                  ,jpband           ,&
        & p_nh_state(jg)% metrics% z_mc(:,:,jb)                     ,&
        & p_nh_state(jg)% metrics% ddqz_z_full(:,:,jb)              ,&
        & aer_tau_sw_vr    ,aer_piz_sw_vr         ,aer_cg_sw_vr     ,&
        & aer_tau_lw_vr                                              )
      
   CASE (iRadAeroKinneVolc)
      CALL set_bc_aeropt_kinne( current_date                        ,&
        & jg                                                        ,&
        & 1,   jce         ,kbdim                 ,klev             ,&
        & jb               ,jpsw                  ,jpband           ,&
        & p_nh_state(jg)% metrics% z_mc(:,:,jb)                     ,&
        & p_nh_state(jg)% metrics% ddqz_z_full(:,:,jb)              ,&
        & aer_tau_sw_vr    ,aer_piz_sw_vr         ,aer_cg_sw_vr     ,&
        & aer_tau_lw_vr                                              )
      CALL add_bc_aeropt_cmip6_volc( current_date ,jg,            1 ,&      
        & jce              ,kbdim                 ,klev             ,&
        & jb               ,jpsw                  ,jpband           ,&
        & p_nh_state(jg)% metrics% z_mc(:,:,jb)                     ,&
        & p_nh_state(jg)% metrics% ddqz_z_full(:,:,jb)              ,&
        & aer_tau_sw_vr    ,aer_piz_sw_vr         ,aer_cg_sw_vr     ,&
        & aer_tau_lw_vr                                              )

   CASE DEFAULT
      WRITE (c_irad_aero,'(i3)') irad_aero
      CALL finish ('rrtm_interface of mo_radition','irad_aero= '// &
                   TRIM(ADJUSTL(c_irad_aero))//' does not exist')
    END SELECT

    IF (jce >= jcs) THEN
      CALL newcld_optics(                                                       &
        & jce          ,kbdim        ,klev         ,jpband       ,jpsw         ,&
        & zglac        ,zland        ,ktype        ,icldlyr      ,tk_fl_vr     ,&
        & zlwp_vr      ,ziwp_vr      ,zlwc_vr      ,ziwc_vr      ,cdnc_vr      ,&
        & reff_liq_vr  ,reff_frz_vr  ,icpl_reff                                ,& 
        & cld_tau_lw_vr,cld_tau_sw_vr,cld_piz_sw_vr,cld_cg_sw_vr                )
    ENDIF

    IF (timers_level > 7) CALL timer_stop(timer_rrtm_prep)

    !
    ! 4.0 Radiative Transfer Routines
    ! --------------------------------
    IF (timers_level > 7) CALL timer_start(timer_lrtm)
    IF (jce >= jcs) THEN
      CALL lrtm(                                                                &
        !    input
        &    jce             ,klev                                             ,&
        &    pm_fl_vr        ,pm_sfc          ,tk_fl_vr        ,tk_hl_vr       ,&
        &    tk_sfc          ,wkl_vr          ,wx_vr           ,col_dry_vr     ,&
        &    zsemiss         ,cld_frc_vr      ,cld_tau_lw_vr   ,aer_tau_lw_vr  ,&
        !    output
        &    flx_uplw_vr     ,flx_dnlw_vr     ,flx_uplw_clr_vr,flx_dnlw_clr_vr )
    ENDIF

    IF (timers_level > 7) CALL timer_stop(timer_lrtm)

    IF (timers_level > 7) CALL timer_start(timer_srtm)
    IF (jce >= jcs) THEN
      CALL srtm_srtm_224gp(                                                     &
        !    input
        &    jce             ,kbdim           ,klev            ,jpsw           ,&
        &    alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif    ,&
        &    pm_fl_vr        ,tk_fl_vr        ,pmu0                            ,&
        &    col_dry_vr      ,wkl_vr                                           ,&
        &    cld_frc_vr      ,cld_tau_sw_vr   ,cld_cg_sw_vr    ,cld_piz_sw_vr  ,&
        &    aer_tau_sw_vr   ,aer_cg_sw_vr    ,aer_piz_sw_vr                   ,&
        &    ssi_radt                                                          ,&
        !    output
        &    flx_dnsw        ,flx_upsw        ,flx_dnsw_clr    ,flx_upsw_clr,   &
        !    optional output
        &    flxd_dff_sfc=flx_dnsw_diff_sfc ,                                   &
        &    flxd_par_sfc    = flx_dnpar_sfc,                                   &
        &    vis_frc_sfc     = vis_frc_sfc,                                     &
        &    nir_dff_frc_sfc = nir_dff_frc_sfc,                                 &
        &    vis_dff_frc_sfc = vis_dff_frc_sfc,                                 &
        &    par_dff_frc_sfc = par_dff_frc_sfc                                  )
    ENDIF
    IF (timers_level > 7) CALL timer_stop(timer_srtm)


    ! 5.0 Post Processing
    ! --------------------------------
    IF (timers_level > 7) CALL timer_start(timer_rrtm_post)

    DO jk = 1, klev+1
      jkb = klev+2-jk
      DO jl = jcs, jce
        flx_lw_net(jl,jk)     = flx_dnlw_vr(jl,jkb) - flx_uplw_vr(jl,jkb)
        flx_lw_net_clr(jl,jk) = flx_dnlw_clr_vr(jl,jkb) - flx_uplw_clr_vr(jl,jkb)
        flx_sw_net(jl,jk)     = flx_dnsw(jl,jk) - flx_upsw(jl,jk)
        flx_sw_net_clr(jl,jk) = flx_dnsw_clr(jl,jk) - flx_upsw_clr(jl,jk)
      END DO
    END DO

    IF (atm_phy_nwp_config(jg)%l_3d_rad_fluxes) THEN
      DO jk = 1, klev+1
        jkb = klev+2-jk
        DO jl = jcs, jce
          flx_lw_dn(jl,jk)      = flx_dnlw_vr(jl,jkb) 
          flx_sw_dn(jl,jk)      = flx_dnsw(jl,jk)
          flx_lw_up(jl,jk)      = flx_uplw_vr(jl,jkb)    
          flx_sw_up(jl,jk)      = flx_upsw(jl,jk)        
          flx_lw_dn_clr(jl,jk)  = flx_dnlw_clr_vr(jl,jkb)
          flx_sw_dn_clr(jl,jk)  = flx_dnsw_clr(jl,jk)    
          flx_lw_up_clr(jl,jk)  = flx_uplw_clr_vr(jl,jkb)
          flx_sw_up_clr(jl,jk)  = flx_upsw_clr(jl,jk)    
        END DO
      END DO
    END IF

    DO jl = jcs, jce
      flx_uplw_sfc(jl)     = flx_uplw_vr(jl,1)
      flx_uplw_sfc_clr(jl) = flx_uplw_clr_vr(jl,1)
      flx_upsw_sfc(jl)     = flx_upsw(jl,klev+1)
      flx_upsw_sfc_clr(jl) = flx_upsw_clr(jl,klev+1)
    END DO
    IF (PRESENT(flx_upsw_toa)) THEN
      DO jl = jcs, jce
        flx_upsw_toa(jl) = flx_upsw(jl,1)
      END DO
    END IF
    !


    IF (timers_level > 7) CALL timer_stop(timer_rrtm_post)

  END SUBROUTINE rrtm_interface


  !-----------------------------------------------------------------------------
  !>
  !! Compute shortwave and longwave heating rates
  !!
  !! The radheat subroutine computes the radiative heating rates resulting from
  !! the divergence of the vertical profiles of longwave and shortwave net fluxes.
  !!
  !! - Shortwave net flux profiles are computed from:
  !!   - the vertical profiles of net transmissivity
  !!   - the solar incoming flux at TOA
  !! - Longwave net flux profiles are given as input
  !! - Specific heat depends on the moisture in the air
  !!

  SUBROUTINE radheat (jcs, jce, jg, kbdim, &
    &                 klev  , klevp1,  &
    &                 ntiles        ,  &
    &                 ntiles_wtr    ,  &
    &                 pmair         ,  &
    &                 pqv           ,  &
    &                 pcd           ,  &
    &                 pcv           ,  &
    &                 pi0           ,  &
    &                 pemiss        ,  &
    &                 ptsfc         ,  &
    &                 ptsfc_t       ,  & ! optional: tile-specific ground temperature
    &                 ptsfctrad     ,  &
    &                 pqc           ,  & ! optional ! must be present if opt_nh_corr=.true.
    &                 pqi           ,  & ! optional ! must be present if opt_nh_corr=.true.
    &                 ppres_ifc     ,  & ! optional ! must be present if opt_nh_corr=.true.
    &                 albedo, albedo_t,& ! optional: albedo fields
    &                 lwflx_up_sfc_rs, & ! optional: longwave upward flux at surface
    &                 trsol_up_toa,    & ! optional: normalized shortwave upward flux at the top of the atmosphere
    &                 trsol_up_sfc,    & ! optional: normalized shortwave upward flux at the surface
    &                 trsol_nir_sfc,   & ! optional: normalized near-infrared downward flux at the surface
    &                 trsol_vis_sfc,   & ! optional: normalized visible downward flux at the surface
    &                 trsol_par_sfc,   & ! optional: normalized photosynthetically active downward flux at the surface
    &                 trsol_dn_sfc_diff,&! optional: normalized shortwave diffuse downward radiative flux at the surface
    &                 trsol_clr_sfc,   & ! optional: normalized shortwave clear-sky net radiative flux at the surface
    &                 list_land_count, & ! optional: number of land points
    &                 list_land_idx,   & ! optional: index list of land points
    &                 list_seawtr_count,&! optional: number of water points
    &                 list_seawtr_idx, & ! optional: index list of water points
    &                 list_seaice_count,&! optional: number of seaice points
    &                 list_seaice_idx, & ! optional: index list of seaice points
    &                 list_lake_count, & ! optional: number of lake points
    &                 list_lake_idx,   & ! optional: index list of lake points
    &                 gp_count_t,      & ! optional: number of land points per tile 
    &                 idx_lst_t,       & ! optional: index list of land points per tile
    &                 cosmu0,          & ! optional: cosine of zenith angle
    &                 cosmu0_slp,      & ! optional: slope-dependent cosine of zenith angle
    &                 shading_mask,    & ! optional: orographic shading mask
    &                 skyview,         & ! optional: skyview factor for islope_rad=2
    &                 opt_nh_corr   ,  & ! optional: switch for applying corrections for NH model
    &                 use_trsolclr_sfc,& ! optional: use clear-sky surface transmissivity passed on input
    &                 ptrmsw        ,  &
    &                 pflxlw        ,  &
    &                 pdtdtradsw    ,  &
    &                 pdtdtradlw    ,  &
    &                 pflxsfcsw     ,  &
    &                 pflxsfcsw_os   ,  &
    &                 pflxsfcsw_tan_os   ,  &
    &                 pflxsfclw     ,  &
    &                 pflxsfcsw_t   ,  &
    &                 pflxsfclw_t   ,  &
    &                 pflxtoasw     ,  &
    &                 pflxtoalw     ,  &
    &                 lwflx_up_sfc  ,  &
    &                 swflx_up_toa  ,  &
    &                 swflx_up_sfc  ,  &
    &                 swflx_up_sfc_os,  &
    &                 swflx_up_sfc_tan_os,  &
    &                 swflx_nir_sfc ,  &
    &                 swflx_vis_sfc ,  &
    &                 swflx_par_sfc ,  &
    &                 swflx_par_sfc_tan_os,  &
    &                 swflx_clr_sfc ,  &
    &                 swflx_dn_sfc_diff, &
    &                 lacc             )

    INTEGER,  INTENT(in)  ::    &
      &     jcs, jce, jg, kbdim,    &
      &     klev, klevp1, ntiles, ntiles_wtr

    REAL(wp), INTENT(in)  ::           &
      &     pmair      (kbdim,klev),   & ! mass of air in layer                     [kg/m2]
      &     pqv        (kbdim,klev),   & ! specific humidity at t-dt                [kg/kg]
      &     pcd,                       & ! specific heat of dry air                 [J/kg/K]
      &     pcv,                       & ! specific heat of vapor                   [J/kg/K]
      &     pi0        (kbdim),        & ! local solar incoming flux at TOA         [W/m2]
      &     pemiss     (kbdim),        & ! lw sfc emissivity
      &     ptsfc      (kbdim),        & ! surface temperature at t                 [K]
      &     ptsfctrad  (kbdim),        & ! surface temperature at trad              [K]
      &     ptrmsw     (kbdim,klevp1), & ! shortwave transmissivity at trad         []
      &     pflxlw     (kbdim,klevp1)    ! longwave net flux at trad                [W/m2]

    REAL(wp), INTENT(in), OPTIONAL  ::         &
      &     pqc       (kbdim,klev),  & ! specific cloud water               [kg/kg]
      &     pqi       (kbdim,klev),  & ! specific cloud ice                 [kg/kg]
      &     ppres_ifc (kbdim,klevp1),& ! pressure at interfaces             [Pa]
      &     ptsfc_t   (:,:),         & ! tile-specific surface temperature at t  [K]
                                       ! dim: (kbdim,,ntiles+ntiles_wtr)
      &     cosmu0    (kbdim),       & ! cosine of solar zenith angle (w.r.t. plain surface)
      &     cosmu0_slp(kbdim),       & ! slope-dependent cosine of solar zenith angle
      &     shading_mask(kbdim),     & ! shading mask
      &     skyview   (:),           & ! skyview factor for islope_rad=2
      &     albedo    (kbdim),       & ! grid-box average albedo
      &     albedo_t  (:,:),         & ! tile-specific albedo
                                       ! dim: (kbdim,ntiles+ntiles_wtr)
      &     lwflx_up_sfc_rs(kbdim),& ! longwave upward flux at surface calculated at radiation time steps
      &     trsol_up_toa(kbdim),   & ! normalized shortwave upward flux at the top of the atmosphere
      &     trsol_up_sfc(kbdim),   & ! normalized shortwave upward flux at the surface
      &     trsol_nir_sfc(kbdim),  & ! normalized near-infrared downward flux at the surface
      &     trsol_vis_sfc(kbdim),  & ! normalized visible downward flux at the surface
      &     trsol_par_sfc(kbdim),  & ! normalized photosynthetically active downward flux at the surface
      &     trsol_clr_sfc(kbdim),  & ! normalized shortwave clear-sky net radiative flux at the surface
      &     trsol_dn_sfc_diff(kbdim) ! normalized shortwave diffuse downward radiative flux at the surface

    INTEGER, INTENT(in), OPTIONAL  ::   &
      &     list_land_count,            &  ! number of land points
      &     list_land_idx(kbdim),       &  ! index list of land points
      &     list_seawtr_count,          &  ! number of water points
      &     list_seawtr_idx(kbdim),     &  ! index list of water points
      &     list_lake_count,            &  ! number of lake points
      &     list_lake_idx(kbdim),       &  ! index list of lake points
      &     list_seaice_count,          &  ! number of seaice points
      &     list_seaice_idx(kbdim),     &  ! index list of seaice points
      &     gp_count_t(:),              &  ! number of land points per tile
                                           ! dim: (ntiles)
      &     idx_lst_t(:,:)                 ! index list of land points per tile
                                           ! dim: (kbdim,ntiles)

    LOGICAL, INTENT(in), OPTIONAL   ::  &
      &     opt_nh_corr, use_trsolclr_sfc

    REAL(wp), INTENT(inout) ::       &
      &     pdtdtradsw (kbdim,klev), & ! shortwave temperature tendency           [K/s]
      &     pdtdtradlw (kbdim,klev)    ! longwave temperature tendency            [K/s]

    REAL(wp), INTENT(inout), OPTIONAL :: &
      &     pflxsfcsw (kbdim),       & ! shortwave surface net flux [W/m2]
      &     pflxsfcsw_os (kbdim),     & ! shortwave surface net flux uncorr. [W/m2]
      &     pflxsfcsw_tan_os (kbdim),     & ! shortwave surface net flux uncorr. [W/m2]
      &     pflxsfclw (kbdim),       & ! longwave  surface net flux [W/m2]
      &     pflxsfcsw_t(:,:),        & ! tile-specific shortwave surface net flux [W/m2]
                                       ! dim: (kbdim,ntiles+ntiles_wtr)
      &     pflxsfclw_t(:,:),        & ! tile-specific longwave surface net flux [W/m2]
                                       ! dim: (kbdim,ntiles+ntiles_wtr)
      &     pflxtoasw (kbdim), &       ! shortwave toa net flux [W/m2]
      &     pflxtoalw (kbdim), &       ! longwave  toa net flux [W/m2]
      &     lwflx_up_sfc(kbdim), &     ! longwave upward flux at surface [W/m2]
      &     swflx_up_toa(kbdim), &     ! shortwave upward flux at the top of the atmosphere [W/m2]
      &     swflx_up_sfc(kbdim), &     ! shortwave upward flux at the surface [W/m2]
      &     swflx_up_sfc_os(kbdim), &  ! shortwave upward flux at the surface [W/m2] incl. orographic shading
      &     swflx_up_sfc_tan_os(kbdim), & ! shortwave upward flux at the surface [W/m2] incl. slope-dependent and orographic shading
      &     swflx_nir_sfc(kbdim), &    ! near-infrared downward flux at the surface [W/m2]
      &     swflx_vis_sfc(kbdim), &    ! visible downward flux at the surface [W/m2]
      &     swflx_par_sfc(kbdim), &    ! photosynthetically active downward flux at the surface [W/m2]
      &     swflx_par_sfc_tan_os(kbdim), &    ! photosynthetically active downward flux at the surface [W/m2]
      &     swflx_clr_sfc(kbdim), &    ! clear-sky net shortwave flux at the surface [W/m2]
      &     swflx_dn_sfc_diff(kbdim)   ! shortwave diffuse downward radiative flux at the surface [W/m2]

    LOGICAL, OPTIONAL, INTENT(in) :: lacc ! GPU flag

    LOGICAL :: lzacc

    ! Local arrays
    REAL(wp) ::                    &
      &     zflxsw (kbdim,klevp1), &
      &     zflxlw (kbdim,klevp1), &
      &     zconv  (kbdim,klev)  , &
      &     tqv    (kbdim)       , &
      &     dlwem_o_dtg(kbdim)   , &
      &     lwfac1 (kbdim)       , &
      &     lwfac2 (kbdim)       , &
      &     intclw (kbdim,klevp1), &
      &     intcli (kbdim,klevp1), &
      &     dlwflxall_o_dtg(kbdim,klevp1)

    REAL(wp) :: dflxsw_o_dalb(kbdim), trsolclr(kbdim), logtqv, shading_corr(kbdim), slopeshade_corr(kbdim)

    ! local scalars
    REAL(wp) :: dpresg, pfaclw, intqctot, dlwflxclr_o_dtg, solrad, angle_ratio, &
      &         swfac1, swfac2

    REAL(wp), PARAMETER  :: pscal = 1._wp/4000._wp ! pressure scale for longwave correction

    INTEGER :: jc,jk,jt,ic

    LOGICAL  :: l_nh_corr, lcalc_trsolclr, lcalc_clrflx

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: zflxsw,zflxlw,zconv,tqv
!DIR$ ATTRIBUTES ALIGN : 64 :: dlwem_o_dtg,lwfac1,lwfac2,intclw,intcli
!DIR$ ATTRIBUTES ALIGN : 64 :: dlwflxall_o_dtg,dflxsw_o_dalb
!DIR$ ATTRIBUTES ALIGN : 64 :: trsolclr,shading_corr,slopeshade_corr
#endif
    IF ( PRESENT(opt_nh_corr) ) THEN
      l_nh_corr = opt_nh_corr
      IF (l_nh_corr .AND. .NOT.(PRESENT(pqc).AND.PRESENT(pqi).AND.PRESENT(ppres_ifc))) THEN
        CALL finish('radheat','error: if opt_nh_corr, pqc, pqi, and ppres_ifc must be present.')
      ENDIF
    ELSE
      l_nh_corr = .FALSE.
    ENDIF

    IF (islope_rad(jg) >= 2 .AND. .NOT. (PRESENT(skyview) )) THEN
      ! we need skyview
      CALL finish('radheat', 'I/O field skyview is missing')
    ENDIF

    IF ( (islope_rad(jg)> 0).AND. .NOT. &
        (PRESENT(pflxsfcsw_os)    .AND. PRESENT(pflxsfcsw_tan_os) .AND. PRESENT(swflx_par_sfc_tan_os) .AND. &
         PRESENT(swflx_up_sfc_os) .AND. PRESENT(swflx_up_sfc_tan_os)) ) THEN
      CALL finish('radheat', 'output fields for slope-dependent shortwave fluxes are missing')
    ENDIF

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA CREATE(zflxsw, zflxlw, zconv, tqv, dlwem_o_dtg) &
    !$ACC   CREATE(lwfac1, lwfac2, intclw, intcli, dlwflxall_o_dtg) &
    !$ACC   CREATE(dflxsw_o_dalb, trsolclr) &
    !$ACC   CREATE(shading_corr, slopeshade_corr) IF(lzacc)

    lcalc_trsolclr = .TRUE.
    IF (PRESENT(use_trsolclr_sfc) .AND. PRESENT(trsol_clr_sfc)) THEN
      IF (use_trsolclr_sfc) lcalc_trsolclr = .FALSE.
    ENDIF
    IF (PRESENT(swflx_clr_sfc)) THEN
      lcalc_clrflx = .TRUE.
    ELSE
      lcalc_clrflx = .FALSE.
    ENDIF

 
    ! Conversion factor for heating rates
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1, klev
      DO jc = jcs, jce
        zconv(jc,jk) = 1._wp/(pmair(jc,jk)*(pcd+(pcv-pcd)*pqv(jc,jk)))
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    ! preset of slope correction factor for solar radiation
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs, jce
      shading_corr(jc) = 1._wp
      slopeshade_corr(jc) = 1._wp
    ENDDO
    !$ACC END PARALLEL

    ! lev == 1        => TOA
    ! lev in [2,klev] => Atmosphere
    ! lev == klevp1   => Surface
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1, klevp1
      DO jc = jcs, jce 
        zflxsw(jc,jk)      = ptrmsw(jc,jk) * pi0(jc)
      END DO
    END DO
    !$ACC END PARALLEL

    ! Longwave fluxes
    ! - TOA
!    zflxlw(jcs:jce,1)      = pflxlw(jcs:jce,1)
!    ! - Atmosphere
!    zflxlw(jcs:jce,2:klev) = pflxlw(jcs:jce,2:klev)

    IF (l_nh_corr) THEN !

      IF (islope_rad(jg) >= 1) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR PRIVATE(solrad, angle_ratio)
        DO jc = jcs, jce
          solrad = MAX(0.1_wp,ptrmsw(jc,klevp1)/(1._wp-albedo(jc)),trsol_dn_sfc_diff(jc))
          angle_ratio = MIN(10._wp,cosmu0_slp(jc)/MAX(1.e-5_wp,cosmu0(jc)))

          shading_corr(jc)    = (trsol_dn_sfc_diff(jc) + (solrad-trsol_dn_sfc_diff(jc))*shading_mask(jc))/solrad
          slopeshade_corr(jc) = (trsol_dn_sfc_diff(jc) + (solrad-trsol_dn_sfc_diff(jc))*angle_ratio*shading_mask(jc))/solrad
        ENDDO
        !$ACC END PARALLEL
      ENDIF

      ! Additional shortwave fluxes for NWP requirements
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = jcs, jce
        swflx_up_toa(jc)      = pi0(jc)*trsol_up_toa(jc)
        swflx_up_sfc(jc)      = pi0(jc)*trsol_up_sfc(jc)
        swflx_dn_sfc_diff(jc) = pi0(jc)*trsol_dn_sfc_diff(jc)
        swflx_nir_sfc(jc)     = pi0(jc)*trsol_nir_sfc(jc)
        swflx_vis_sfc(jc)     = pi0(jc)*trsol_vis_sfc(jc)
        swflx_par_sfc(jc)     = pi0(jc)*trsol_par_sfc(jc)
      ENDDO
      !$ACC END PARALLEL


      IF (islope_rad(jg)> 0) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO jc = jcs, jce
          swflx_up_sfc_os(jc)      = pi0(jc)*trsol_up_sfc(jc)*shading_corr(jc)
          swflx_up_sfc_tan_os(jc)  = pi0(jc)*trsol_up_sfc(jc)*slopeshade_corr(jc)
          swflx_par_sfc_tan_os(jc) = pi0(jc)*trsol_par_sfc(jc)*slopeshade_corr(jc)
        ENDDO
        !$ACC END PARALLEL
      ENDIF


      IF (lcalc_clrflx) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO jc = jcs, jce
          swflx_clr_sfc(jc) = pi0(jc)*trsol_clr_sfc(jc)
        ENDDO
        !$ACC END PARALLEL
      ENDIF

      ! Correction of longwave fluxes for changes in ground temperature
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = 1, kbdim
        tqv(jc)            = 0._wp
        intclw(jc,klevp1)  = 0._wp
        intcli(jc,klevp1)  = 0._wp
      ENDDO
      !$ACC END PARALLEL
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO jk = klev,1,-1
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(dpresg)
        DO jc = jcs, jce
          dpresg        = (ppres_ifc(jc,jk+1) - ppres_ifc(jc,jk))/grav
          tqv(jc)       = tqv(jc)+pqv(jc,jk)*dpresg
          intclw(jc,jk) = intclw(jc,jk+1)+pqc(jc,jk)*dpresg
          intcli(jc,jk) = intcli(jc,jk+1)+pqi(jc,jk)*dpresg
        ENDDO
      ENDDO

      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(logtqv)
      DO jc = jcs, jce
        logtqv = LOG(MAX(1._wp,tqv(jc)))
        dlwem_o_dtg(jc) = pemiss(jc)*4._wp*stbo*ptsfc(jc)**3
        lwflx_up_sfc(jc) = lwflx_up_sfc_rs(jc) + dlwem_o_dtg(jc)*(ptsfc(jc) - ptsfctrad(jc))
        lwfac2(jc) = 0.92_wp*EXP(-0.07_wp*logtqv)
        lwfac1(jc) = MERGE(1.677_wp, 0.4388_wp, tqv(jc) > 15._wp) &
             * EXP(MERGE(-0.72_wp, -0.225_wp, tqv(jc) > 15._wp) *logtqv)
      ENDDO

      !$ACC LOOP SEQ
      DO jk = 1,klevp1
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(pfaclw, intqctot, dlwflxclr_o_dtg)
        DO jc = jcs, jce
          pfaclw = lwfac1(jc)+(lwfac2(jc)-lwfac1(jc))*EXP(-SQRT((ppres_ifc(jc,klevp1)- &
            ppres_ifc(jc,jk))*pscal))
          intqctot = MIN(0.30119_wp,MAX(1.008e-3_wp,intclw(jc,jk)+0.2_wp*intcli(jc,jk)))

          dlwflxclr_o_dtg = -dlwem_o_dtg(jc)*pfaclw

          ! derivative of LW flux w.r.t. ground temperature
          dlwflxall_o_dtg(jc,jk) = dlwflxclr_o_dtg*(1._wp-(6.9_wp+LOG(intqctot))/5.7_wp)
          ! Now apply the correction
          zflxlw(jc,jk) =  pflxlw(jc,jk) + dlwflxall_o_dtg(jc,jk) * (ptsfc(jc) - ptsfctrad(jc))
        ENDDO
      ENDDO
      !$ACC END PARALLEL

      ! Disaggregation of longwave and shortwave fluxes for tile approach
      IF (ntiles > 1) THEN
        IF (lcalc_trsolclr) THEN
          ! parameterization of clear-air solar transmissivity in order to use the same
          ! formulation as in mo_phys_nest_utilities:downscale_rad_output
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG VECTOR
          DO jc = jcs, jce
            trsolclr(jc) = MAX(0.02_wp,0.8_wp*cosmu0(jc)/(0.25_wp*tqv(jc))**0.15)**0.333_wp*&
             (1._wp-albedo(jc))**(1._wp-0.2_wp*cosmu0(jc)+0.1_wp*MIN(10._wp,tqv(jc))**0.33_wp)
          ENDDO
          !$ACC END PARALLEL
        ELSE
          ! use clear-air solar transmissivity passed as argument
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG VECTOR
          DO jc = jcs, jce
            trsolclr(jc) = trsol_clr_sfc(jc)
          ENDDO
          !$ACC END PARALLEL
        ENDIF

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR PRIVATE(swfac1, swfac2)
        DO jc = jcs, jce
          swfac1 = EXP(0.36_wp*LOG( MAX(1.e-3_wp,ptrmsw(jc,klevp1))/MAX(1.e-3_wp,trsolclr(jc)) ))
          swfac2 = EXP(0.1_wp*LOG( MAX(0.25_wp,3._wp*cosmu0(jc)) ))

          ! derivative of SW surface flux w.r.t. albedo
          dflxsw_o_dalb(jc) = - zflxsw(jc,klevp1)*swfac1/((1._wp-albedo(jc))*swfac2)
        ENDDO
        !$ACC END PARALLEL
 
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP SEQ
        DO jt = 1,ntiles
          !$ACC LOOP GANG VECTOR PRIVATE(jc)
          DO ic = 1, gp_count_t(jt)
            jc = idx_lst_t(ic,jt)
            pflxsfcsw_t(jc,jt) = slopeshade_corr(jc) * MAX(0.1_wp*zflxsw(jc,klevp1), zflxsw(jc,klevp1) + &
                                 dflxsw_o_dalb(jc)*(albedo_t(jc,jt)-albedo(jc)))
            pflxsfclw_t(jc,jt) = zflxlw(jc,klevp1) + dlwflxall_o_dtg(jc,klevp1)* &
                                 (ptsfc_t(jc,jt)-ptsfc(jc))
          ENDDO
        ENDDO
        !$ACC END PARALLEL

        ! seaice points
        !
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR PRIVATE(jc)
        DO ic = 1, list_seaice_count
          jc = list_seaice_idx(ic)
          pflxsfcsw_t(jc,isub_seaice) = MAX(0.1_wp*zflxsw(jc,klevp1), zflxsw(jc,klevp1) &
            &                  + dflxsw_o_dalb(jc)*(albedo_t(jc,isub_seaice)-albedo(jc)))
          pflxsfclw_t(jc,isub_seaice) = zflxlw(jc,klevp1) + dlwflxall_o_dtg(jc,klevp1) &
            &                  * (ptsfc_t(jc,isub_seaice)-ptsfc(jc))
        ENDDO
        !$ACC END PARALLEL

        ! lake points
        !
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR PRIVATE(jc)
        DO ic = 1, list_lake_count
          jc = list_lake_idx(ic)
          pflxsfcsw_t(jc,isub_lake) = MAX(0.1_wp*zflxsw(jc,klevp1), zflxsw(jc,klevp1) &
            &                  + dflxsw_o_dalb(jc)*(albedo_t(jc,isub_lake)-albedo(jc)))
          pflxsfclw_t(jc,isub_lake) = zflxlw(jc,klevp1) + dlwflxall_o_dtg(jc,klevp1) &
            &                  * (ptsfc_t(jc,isub_lake)-ptsfc(jc))
        ENDDO
        !$ACC END PARALLEL

        ! (open) water points (needed for A-O coupling)
        !
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR PRIVATE(jc)
        DO ic = 1, list_seawtr_count
          jc = list_seawtr_idx(ic)
          pflxsfcsw_t(jc,isub_water) = MAX(0.1_wp*zflxsw(jc,klevp1), zflxsw(jc,klevp1) &
            &                  + dflxsw_o_dalb(jc)*(albedo_t(jc,isub_water)-albedo(jc)))
          pflxsfclw_t(jc,isub_water) = zflxlw(jc,klevp1) + dlwflxall_o_dtg(jc,klevp1) &
            &                  * (ptsfc_t(jc,isub_water)-ptsfc(jc))
        ENDDO
        !$ACC END PARALLEL

      ELSE IF (PRESENT(pflxsfcsw_t) .AND. PRESENT(pflxsfclw_t)) THEN

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR PRIVATE(jc)
        DO ic = 1, list_land_count
          jc = list_land_idx(ic)
          pflxsfcsw_t(jc,1) = slopeshade_corr(jc) * zflxsw(jc,klevp1)
          pflxsfclw_t(jc,1) = zflxlw(jc,klevp1)
        ENDDO
        !$ACC END PARALLEL

        ! seaice points
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR PRIVATE(jc)
        DO ic = 1, list_seaice_count
          jc = list_seaice_idx(ic)
          pflxsfcsw_t(jc,1) = zflxsw(jc,klevp1)
          pflxsfclw_t(jc,1) = zflxlw(jc,klevp1)
        ENDDO
        !$ACC END PARALLEL

        ! lake points
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR PRIVATE(jc)
        DO ic = 1, list_lake_count
          jc = list_lake_idx(ic)
          pflxsfcsw_t(jc,1) = zflxsw(jc,klevp1)
          pflxsfclw_t(jc,1) = zflxlw(jc,klevp1)
        ENDDO
        !$ACC END PARALLEL

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR PRIVATE(jc)
        DO ic = 1, list_seawtr_count
          jc = list_seawtr_idx(ic)
          pflxsfcsw_t(jc,1) = zflxsw(jc,klevp1)
          pflxsfclw_t(jc,1) = zflxlw(jc,klevp1)
        ENDDO
        !$ACC END PARALLEL

      ENDIF ! ntiles


    ENDIF

    !
    !
    !     4.2  Fluxes and heating rates except for lowest layer
    !
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1, klev
      DO jc = jcs, jce
        pdtdtradsw(jc,jk) = (zflxsw(jc,jk)-zflxsw(jc,jk+1)) * &
          & zconv(jc,jk)
        pdtdtradlw(jc,jk) = (zflxlw(jc,jk)-zflxlw(jc,jk+1)) * &
          & zconv(jc,jk)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !
    !     4.3 net fluxes at surface
    !
    IF ( PRESENT(pflxsfcsw) ) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = jcs, jce
        pflxsfcsw(jc) = zflxsw(jc,klevp1)
      ENDDO
      !$ACC END PARALLEL
    ENDIF

    IF (islope_rad(jg)> 0) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = jcs, jce
        pflxsfcsw_os(jc)     = shading_corr(jc) * zflxsw(jc,klevp1)
        pflxsfcsw_tan_os(jc) = slopeshade_corr(jc) * zflxsw(jc,klevp1)
      ENDDO
      !$ACC END PARALLEL
    ENDIF

    IF ( PRESENT(pflxsfclw) ) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = jcs, jce
        pflxsfclw(jc) = zflxlw(jc,klevp1)
      ENDDO
      !$ACC END PARALLEL
    ENDIF

    !
    !     4.4 net sw flux at toa
    !
    IF ( PRESENT(pflxtoasw) ) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = jcs, jce
        pflxtoasw(jc) = zflxsw(jc,1)
      ENDDO
      !$ACC END PARALLEL
    ENDIF
    IF ( PRESENT(pflxtoalw) ) THEN 
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = jcs, jce
        pflxtoalw(jc) = zflxlw(jc,1)
      ENDDO
      !$ACC END PARALLEL
    ENDIF

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE radheat

END MODULE mo_radiation
