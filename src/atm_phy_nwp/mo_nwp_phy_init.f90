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

! Initialize the physical schemes at start time

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nwp_phy_init

  USE mo_kind,                ONLY: wp
  USE mo_math_constants,      ONLY: rad2deg
  USE mo_physical_constants,  ONLY: grav, rd_o_cpd, cpd, p0ref, rd, p0sl_bg,         &
    &                               dtdz_standardatm, lh_v=>alv, o3mr2gg, tf_salt
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag,t_nwp_phy_tend
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_ext_data_init,       ONLY: diagnose_ext_aggr, vege_clim
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_vertical_coord_table,ONLY: vct_a
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: min_rlcell, min_rlcell_int, io3_ape,            &
    &                               MODE_COMBINED, MODE_IFSANA, icosmo, ismag,      &
    &                               iprog, igme, SUCCESS,                           &
    &                               MODE_COSMO, MODE_ICONVREMAP, iss, iorg, ibc,    &
    &                               iso4, idu, LSS_JSBACH, LSS_TERRA, ivdiff
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_parallel_config,     ONLY: nproma
  USE mo_fortran_tools,       ONLY: copy
  USE mo_run_config,          ONLY: ltestcase, iqv, iqc, inccn, ininpot, msg_level 
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config, lrtm_filename,               &
    &                               cldopt_filename, icpl_aero_conv, iprog_aero
  USE mo_extpar_config,       ONLY: ext_o3_attr, itype_vegetation_cycle

  !radiation
  USE mo_newcld_optics,       ONLY: setup_newcld_optics
  USE mo_lrtm_setup,          ONLY: lrtm_setup
  USE mo_radiation_config,    ONLY: irad_aero, iRadAeroTegen, iRadAeroART,            &
    &                               iRadAeroCAMSclim, iRadAeroCAMStd,                 &
    &                               iRadAeroConstKinne, iRadAeroKinne, iRadAeroVolc,  &
    &                               iRadAeroKinneVolc,  iRadAeroKinneVolcSP,          &
    &                               iRadAeroKinneSP,                                  &
    &                               ssi_radt, tsi_radt,irad_o3, rad_csalbw,           &
    &                               ghg_filename, irad_co2, irad_cfc11, irad_cfc12,   &
    &                               irad_n2o, irad_ch4, isolrad
  USE mo_nwp_aerosol,         ONLY: nwp_aerosol_init
  USE mo_srtm_config,         ONLY: setup_srtm, ssi_amip, ssi_coddington
  USE mo_aerosol_util,        ONLY: init_aerosol_props_tegen_rrtm,                  &
    &                               zaea_rrtm, zaes_rrtm, zaeg_rrtm
  USE mo_o3_util,             ONLY: o3_pl2ml!, o3_zl2ml
#ifdef __ECRAD
  USE mo_nwp_ecrad_init,      ONLY: setup_ecrad
  USE mo_ecrad,               ONLY: ecrad_conf, IGasModelIFSRRTMG,                  &
    &                               ecrad_ssi_default, ecrad_ssi_coddington
  USE mo_aerosol_util,        ONLY: init_aerosol_props_tegen_ecrad
#endif

  USE mo_2mom_mcrph_driver,   ONLY: two_moment_mcrph_init
  USE microphysics_1mom_schemes, ONLY: microphysics_1mom_init
  USE mo_sbm_util,            ONLY: sbm_init 

#ifdef __ICON_ART
  USE mo_art_clouds_interface,ONLY: art_clouds_interface_2mom_init
#endif
  USE mo_cpl_aerosol_microphys, ONLY: lookupcreate_segalkhain, specccn_segalkhain_simple, &
                                      ncn_from_tau_aerosol_speccnconst

  ! convection
  USE mo_cuparameters,        ONLY: sucst,  sucumf,    &
    &                               su_yoethf,         &
    &                               sucldp, suphli,    &
    &                               suvdf , suvdfs

  ! turbulence
  USE mo_turbdiff_config,     ONLY: turbdiff_config
  USE turb_data,              ONLY: get_turbdiff_param, lsflcnd, &
                                    impl_s, impl_t,              &
                                    impl_weight,                 &
                                    z0_ice, & 
                                    imode_pat_len, pat_len, ndim
  USE turb_transfer,          ONLY: turbtran
  USE turb_diffusion,         ONLY: turbdiff
  USE mo_nwp_vdiff_interface, ONLY: nwp_vdiff_init, nwp_vdiff_update_seaice_list
  USE mo_turb_vdiff_config,   ONLY: vdiff_config

  USE mo_nwp_sfc_utils,       ONLY: nwp_surface_init, init_snowtile_lists, init_sea_lists, &
    &                               aggregate_tg_qvs, copy_lnd_prog_now2new
  USE mo_lnd_nwp_config,      ONLY: ntiles_total, lsnowtile, ntiles_water, llake, &
    &                               lseaice, zml_soil, nlev_soil, dzsoil_icon => dzsoil
  USE sfc_flake_data,         ONLY: h_Ice_min_flk, tpl_T_f
  USE sfc_terra_data,         ONLY: csalbw, cpwp, cfcap
  USE mo_satad,               ONLY: sat_pres_water, &  !! saturation vapor pressure w.r.t. water
    &                               sat_pres_ice, &    !! saturation vapor pressure w.r.t. ice
    &                               spec_humi          !! Specific humidity

  USE data_gwd,               ONLY: sugwwms

  USE mo_master_config,       ONLY: isRestart
  USE mo_nwp_parameters,      ONLY: t_phy_params

  USE mo_initicon_config,     ONLY: init_mode, lread_tke, itype_sma
  USE mo_apt_routines,        ONLY: init_apt_fields, apply_landalb_tuning, apply_sma
  USE mo_nwp_tuning_config,   ONLY: tune_zceff_min, tune_v0snow, tune_zvz0i, tune_icesedi_exp, tune_box_liq_sfc_fac
  USE mo_cuparameters,        ONLY: sugwd
  USE mtime,                  ONLY: datetime, MAX_DATETIME_STR_LEN, &
    &                               datetimeToString, newDatetime, deallocateDatetime
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights,         &
    &                                  calculate_time_interpolation_weights
  USE mo_timer,               ONLY: timers_level, timer_start, timer_stop,   &
    &                               timer_init_nwp_phy, timer_phys_reff, timer_upatmo
  USE mo_bc_greenhouse_gases, ONLY: read_bc_greenhouse_gases, bc_greenhouse_gases_time_interpolation
  USE mo_nwp_reff_interface,  ONLY: init_reff
  USE mo_upatmo_config,       ONLY: upatmo_config
  USE mo_upatmo_impl_const,   ONLY: iUpatmoPrcStat, iUpatmoStat
#ifndef __NO_ICON_UPATMO__
  USE mo_upatmo_phy_setup,    ONLY: init_upatmo_phy_nwp
#endif

  USE mo_ape_params,          ONLY: ape_sst
  USE mo_nh_testcases_nml,    ONLY: nh_test_name, ape_sst_case, th_cbl, sol_const
  USE mo_grid_config,         ONLY: l_scm_mode
  USE mo_scm_nml,             ONLY: i_scm_netcdf, lscm_read_tke, lscm_read_z0, &
                                    scm_sfc_temp, scm_sfc_qv
  USE mo_nh_torus_exp,        ONLY: read_soil_profile_nc

  USE mo_cover_koe,           ONLY: cover_koe_config
  USE mo_bc_aeropt_kinne,     ONLY: read_bc_aeropt_kinne
  USE mo_bc_aeropt_cmip6_volc,ONLY: read_bc_aeropt_cmip6_volc
  USE mo_bc_aeropt_splumes,   ONLY: setup_bc_aeropt_splumes
  USE mo_bc_ozone,            ONLY: read_bc_ozone
  USE mo_bc_solar_irradiance, ONLY: read_bc_solar_irradiance

  USE mo_sppt_state,          ONLY: sppt
  USE mo_sppt_config,         ONLY: sppt_config
  USE mo_sppt_util,           ONLY: init_rn

  IMPLICIT NONE

  PRIVATE


  PUBLIC  :: init_nwp_phy, init_cloud_aero_cpl, clim_cdnc

  CHARACTER(len=*), PARAMETER :: modname = 'mo_nwp_phy_init'

CONTAINS

SUBROUTINE init_nwp_phy ( p_patch, p_metrics,             &
                       &  p_prog_now,  p_diag,            &
                       &  prm_diag,prm_nwp_tend,          &
                       &  p_prog_lnd_now, p_prog_lnd_new, &
                       &  p_prog_wtr_now, p_prog_wtr_new, &
                       &  p_diag_lnd,                     &
                       &  ext_data, phy_params, ini_date, &
                       &  lnest_start, lreset             )

  TYPE(t_patch),        TARGET,INTENT(in)    :: p_patch
  TYPE(t_nh_metrics),          INTENT(in)    :: p_metrics
  TYPE(t_nh_prog),      TARGET,INTENT(inout) :: p_prog_now !!the prognostic variables
  TYPE(t_nh_diag),      TARGET,INTENT(inout) :: p_diag  !!the diagostic variables
  TYPE(t_external_data),       INTENT(inout) :: ext_data
  TYPE(t_nwp_phy_diag),        INTENT(inout) :: prm_diag
  TYPE(t_nwp_phy_tend), TARGET,INTENT(inout) :: prm_nwp_tend
  TYPE(t_lnd_prog),            INTENT(inout) :: p_prog_lnd_now, p_prog_lnd_new
  TYPE(t_wtr_prog),            INTENT(inout) :: p_prog_wtr_now, p_prog_wtr_new
  TYPE(t_lnd_diag),            INTENT(inout) :: p_diag_lnd
  TYPE(t_phy_params),          INTENT(inout) :: phy_params
  TYPE(datetime),              POINTER       :: ini_date     ! current datetime (mtime)
  LOGICAL, INTENT(IN), OPTIONAL              :: lnest_start, lreset

  INTEGER             :: jk, jk1
  REAL(wp)            :: rsltn   ! horizontal resolution
  REAL(wp)            :: pref(p_patch%nlev)
  REAL(wp)            :: zlat, zprat, zn1, zn2, zcdnc
  REAL(wp)            :: zpres, zpres0
  REAL(wp)            :: l_hori(nproma)
  REAL(wp)            :: scale_fac ! scale factor used only for RCE cases
  REAL(wp) :: zvari   (nproma,p_patch%nlevp1,0:ndim)  !< storage for internal variables and gradients in turbdif and turbtran
  REAL(wp) :: zrhon   (nproma,p_patch%nlevp1)         !< storage for half-level density in turbdiff and turbtran
  REAL(wp)            :: w_so_profile(nlev_soil)    !soil moisture for SCM initialization
  REAL(wp)            :: t_so_profile(nlev_soil+1)  !soil temperature for SCM initialization
  REAL(wp)            :: t_g_in                     !T,g for SCM initialization

  INTEGER             :: icur_date    ! current date converted to integer

  ! Reference atmosphere parameters
  REAL(wp), PARAMETER :: htropo = 11000._wp       ! [m]    tropopause height
  REAL(wp), PARAMETER :: t00    = 288.15_wp       ! [m]    temperature at sea level

  REAL(wp), PARAMETER :: grav_o_rd = grav / rd
  REAL(wp), PARAMETER :: cpd_o_rd  = cpd  / rd

  REAL(wp), PARAMETER :: pr800  = 800._wp / 1013.25_wp
  REAL(wp), PARAMETER :: pr400  = 400._wp / 1013.25_wp
  REAL(wp), PARAMETER :: pr700  = 700._wp / 1013.25_wp

  REAL(wp) :: ttropo, ptropo, temp, zfull

  REAL(wp) :: dz1, dz2, dz3, fact_z0rough
  REAL(wp), ALLOCATABLE :: zrefpres(:,:,:)   ! ref press computed from ref exner
  REAL(wp), ALLOCATABLE :: zreftemp(:,:,:)   ! ref temp computed from ref exner
  REAL(wp), ALLOCATABLE :: zpres_sfc(:,:)    ! ref sfc press
  REAL(wp), ALLOCATABLE :: zpres_ifc(:,:,:)  ! ref press at interfaces

  LOGICAL :: lland, lglac, lshallow, ldetrain_prec, lgrayzone_dc, lrestune_off, lmflimiter_off
  LOGICAL :: lstoch_expl, lstoch_sde,lstoch_deep,lvvcouple,lvv_shallow_deep
  LOGICAL :: ltkeinp_loc, lgz0inp_loc  !< turbtran switches
  LOGICAL :: linit_mode, lturb_init, lreset_mode
  LOGICAL :: lupatmo_phy
  LOGICAL :: l_filename_year
  LOGICAL, DIMENSION(nproma) :: l_lake, l_sice

  INTEGER :: jb,ic,jc,jt,jg,ist,nzprv
  INTEGER :: nlev, nlevp1, nlevcm    !< number of full, half and canopy levels
  INTEGER :: nshift                  !< shift with respect to global grid
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk    !> blocks
  INTEGER :: i_startidx, i_endidx    !! slices
  INTEGER :: i_nchdom                !! domain index
  INTEGER :: lc_class,i_lc_si

  INTEGER :: ierrstat=0
  CHARACTER (LEN=25) :: eroutine
  CHARACTER (LEN=80) :: errormsg

  INTEGER :: nblks_c

  INTEGER :: k1500m                  ! index of first half level above 1500m
  INTEGER :: istatus=0

  REAL(wp) :: hag                    ! height above ground
  REAL(wp) :: h650_standard, h850_standard, h950_standard  ! height of 850hPa and 950hPa level in m 

  REAL(wp) :: N_cn0,z0_nccn,z1e_nccn,N_in0,z0_nin,z1e_nin     ! for CCN and IN in case of gscp=5

  CHARACTER(len=*), PARAMETER ::  &
     routine = modname//':init_nwp_phy'

  CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: datetime_string, yyyymmdd

  IF (timers_level > 3) CALL timer_start(timer_init_nwp_phy)

  ! Local control variable for extended turbulence initializations
  IF (ANY((/MODE_IFSANA,MODE_COMBINED,MODE_COSMO/) == init_mode) .OR. &
      init_mode == MODE_ICONVREMAP .AND. .NOT. lread_tke) THEN
    lturb_init = .TRUE.
  ELSE
    lturb_init = .FALSE.
  ENDIF

  ! This is needed for correct flow control when a nested domain is initialized after restarting
  IF (PRESENT(lnest_start)) THEN
    linit_mode = lnest_start
    lturb_init = .TRUE.
  ELSE
    linit_mode = .NOT. isRestart()
  ENDIF

  IF (PRESENT(lreset)) THEN
    lreset_mode = lreset
  ELSE
    lreset_mode = .FALSE.
  ENDIF

  i_nchdom  = MAX(1,p_patch%n_childdom)

  ierrstat=0
  eroutine=''
  errormsg=''

  ! number of vertical levels
  nlev   = p_patch%nlev
  nlevp1 = p_patch%nlevp1
  jg     = p_patch%id

  nshift = p_patch%nshift_total

  i_lc_si= ext_data%atm%i_lc_snow_ice

  nblks_c = p_patch%nblks_c

  dz1 = 0.0_wp
  dz2 = 0.0_wp
  dz3 = 0.0_wp

  ! Initialization of upper-atmosphere physics 
  ! only in case of no reset and if the upatmo physics are switched on
  ! (upper-atmosphere physics are not integrated into the IAU iterations)
  lupatmo_phy = (.NOT. lreset_mode) .AND. &
    & upatmo_config(jg)%nwp_phy%l_phy_stat( iUpatmoPrcStat%enabled ) 

  IF ( nh_test_name == 'RCE' .OR. nh_test_name == 'RCE_Tconst' .OR. nh_test_name == 'RCE_Tprescr' .OR. nh_test_name == 'RCE_bubble' ) THEN
    ! allocate storage var for press to be used in o3_pl2ml
    ALLOCATE (zrefpres(nproma,nlev,nblks_c),STAT=istatus)
    IF(istatus/=SUCCESS)THEN
      CALL finish(routine, &
                 'allocation of zrefpres failed')
    END IF
    ALLOCATE (zreftemp(nproma,nlev,nblks_c),STAT=istatus)
    IF(istatus/=SUCCESS)THEN
      CALL finish(routine, &
                 'allocation of zreftemp failed')
    END IF
    ALLOCATE (zpres_sfc(nproma,nblks_c),STAT=istatus)
    IF(istatus/=SUCCESS)THEN
      CALL finish(routine, &
                 'allocation of zpres_sfc failed')
    END IF
    ALLOCATE (zpres_ifc(nproma,nlevp1,nblks_c),STAT=istatus)
    IF(istatus/=SUCCESS)THEN
      CALL finish(routine, &
                 'allocation of zpres_ifc failed')
    END IF
    zrefpres = 0.0_wp
    zreftemp = 0.0_wp
    zpres_sfc = 0.0_wp
    zpres_ifc = 0.0_wp
  END IF

  ! for both restart and non-restart runs. Could not be included into
  ! mo_ext_data_state/init_index_lists due to its dependence on fr_seaice
  IF (msg_level >= 12) CALL message('mo_nwp_phy_init:', 'initialize sea-ice lists, call init_sea_lists')
  CALL init_sea_lists(p_patch, lseaice, p_diag_lnd%fr_seaice(:,:), ext_data)
  IF (atm_phy_nwp_config(p_patch%id)%inwp_turb == ivdiff) THEN
    ! TODO: move this to a more appropriate place.
    !$ACC UPDATE DEVICE(p_diag_lnd%fr_seaice) ASYNC(1)
    CALL nwp_vdiff_update_seaice_list ( &
        & p_patch, p_diag_lnd%fr_seaice(:,:), ext_data%atm%list_sea, ext_data%atm%list_seaice, &
        & lacc=.TRUE. &
      )
  END IF

  IF (.NOT. lreset_mode) THEN
    IF (itype_vegetation_cycle >= 2) CALL vege_clim (p_patch, ext_data, p_diag)
    CALL init_apt_fields(p_patch, p_diag, prm_diag, ext_data, p_diag_lnd, p_prog_wtr_now)
    CALL apply_landalb_tuning (p_patch, prm_diag, ext_data)
  ENDIF

  IF (itype_sma >= 2 .AND. .NOT. isRestart()) &
    CALL apply_sma (p_patch, p_diag, ext_data, p_diag_lnd, p_prog_lnd_now)

  ! Diagnose aggregated external parameter fields
  ! (mainly for output purposes)
  ! aggregated sai needed below for organize_turbdiff
  ! This routine is called after init_sea_lists,  
  ! in order to have all tile-related index lists available.
  !
  CALL diagnose_ext_aggr (p_patch, ext_data)

  IF (linit_mode) THEN ! initialize field for time-dependent LW emissivity
!$OMP PARALLEL
    CALL copy(ext_data%atm%emis_rad, prm_diag%lw_emiss, lacc=.FALSE.)
!$OMP END PARALLEL
  ENDIF

  rl_start = 1 ! Initialization should be done for all points
  rl_end   = min_rlcell

  i_startblk = p_patch%cells%start_blk(rl_start,1)
  i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

    ! mask fields to distinguish between tropics and extratropics (needed for some tuning measures)
    !
    DO jc = i_startidx,i_endidx
      zlat = ABS(p_patch%cells%center(jc,jb)%lat*rad2deg)
      IF (zlat < 25._wp) THEN
        prm_diag%tropics_mask(jc,jb) = 1._wp
      ELSE IF (zlat > 30._wp) THEN
        prm_diag%tropics_mask(jc,jb) = 0._wp
      ELSE
        prm_diag%tropics_mask(jc,jb) = (30._wp-zlat)/5._wp
      ENDIF
      IF (zlat < 12.5_wp) THEN
        prm_diag%innertropics_mask(jc,jb) = 1._wp
      ELSE IF (zlat > 17.5_wp) THEN
        prm_diag%innertropics_mask(jc,jb) = 0._wp
      ELSE
        prm_diag%innertropics_mask(jc,jb) = (17.5_wp-zlat)/5._wp
      ENDIF
      IF (zlat > 60._wp) THEN
        prm_diag%sso_lat_mask(jc,jb) = 0._wp
      ELSE IF (zlat < 30._wp) THEN
        prm_diag%sso_lat_mask(jc,jb) = 1._wp
      ELSE
        prm_diag%sso_lat_mask(jc,jb) = COS(3._wp*(zlat-30._wp)/rad2deg)
      ENDIF
    ENDDO

    ! pat_len field for circulation term in turbdiff
    DO jc = i_startidx,i_endidx
      prm_diag%pat_len(jc,jb) = 300._wp*EXP(1.5_wp*LOG(MAX(1.e-2_wp,(ext_data%atm%sso_stdh_raw(jc,jb)-150._wp)/300._wp)))
    ENDDO

  ENDDO

  IF (linit_mode) THEN

    IF (ltestcase) THEN

      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
             &  i_startidx, i_endidx, rl_start, rl_end)

        IF ( (nh_test_name == 'APE_nwp' .OR. nh_test_name == 'dcmip_tc_52' .OR.        &
            & nh_test_name == 'RCEMIP_analytical'  .OR.                                &
            & nh_test_name == 'RCE_Tconst' .OR. nh_test_name == 'RCE_Tprescr' &
            & .OR. nh_test_name == 'CBL_flxconst') ) THEN

          ! t_g = ape_sst1

          DO jc = i_startidx, i_endidx
            zlat = p_patch%cells%center(jc,jb)%lat
            p_prog_lnd_now%t_g  (jc,jb)   = ape_sst(ape_sst_case,zlat) ! set SST
            p_prog_lnd_new%t_g  (jc,jb)   = ape_sst(ape_sst_case,zlat)
            p_prog_lnd_now%t_g_t(jc,jb,1) = ape_sst(ape_sst_case,zlat)
            p_prog_lnd_new%t_g_t(jc,jb,1) = ape_sst(ape_sst_case,zlat)
            ! Humidity at water surface = humidity at saturation
            p_diag_lnd%qv_s(jc,jb)     = &
              &  spec_humi(sat_pres_water(p_prog_lnd_now%t_g(jc,jb)),p_diag%pres_sfc(jc,jb))
            p_diag_lnd%qv_s_t(jc,jb,1) = p_diag_lnd%qv_s(jc,jb)
          END DO


        ELSE IF ( nh_test_name == 'wk82' ) THEN

          DO jc = i_startidx, i_endidx
            p_prog_lnd_now%t_g (jc,jb) = p_diag%temp  (jc,nlev,jb)*  &
                      ((p_diag%pres_sfc(jc,jb))/p_diag%pres(jc,nlev,jb))**rd_o_cpd
            p_prog_lnd_new%t_g (jc,jb) = p_prog_lnd_now%t_g (jc,jb)
            p_prog_lnd_now%t_g_t(jc,jb,1) = p_prog_lnd_now%t_g (jc,jb)
            p_prog_lnd_new%t_g_t(jc,jb,1) = p_prog_lnd_now%t_g (jc,jb)

            p_diag_lnd%qv_s     (jc,jb) = &
              &  spec_humi(sat_pres_water(p_prog_lnd_now%t_g (jc,jb)),p_diag%pres_sfc(jc,jb))
            p_diag_lnd%qv_s    (jc,jb) = MIN (p_diag_lnd%qv_s(jc,jb) ,   &
              &   p_prog_now%tracer(jc,nlev,jb,iqv))
            p_diag_lnd%qv_s_t(jc,jb,1) = p_diag_lnd%qv_s(jc,jb)
          END DO

        ELSE IF ( ( nh_test_name == 'RCE' .OR. nh_test_name == 'RCE_Tconst' &
               & .OR. nh_test_name == 'RCE_Tprescr' )                       &
      &            .AND. atm_phy_nwp_config(jg)%inwp_turb/=ismag) THEN !

          DO jc = i_startidx, i_endidx
            p_prog_lnd_now%t_g  (jc,jb)   = th_cbl(1)
            p_prog_lnd_new%t_g  (jc,jb)   = p_prog_lnd_now%t_g (jc,jb)
            p_prog_lnd_now%t_g_t(jc,jb,1) = p_prog_lnd_now%t_g (jc,jb)
            p_prog_lnd_new%t_g_t(jc,jb,1) = p_prog_lnd_now%t_g (jc,jb)
            p_diag_lnd%qv_s(jc,jb)     = &
              &  spec_humi(sat_pres_water(p_prog_lnd_now%t_g (jc,jb)),p_diag%pres_sfc(jc,jb))
            p_diag_lnd%qv_s_t(jc,jb,1) = p_diag_lnd%qv_s(jc,jb)
          END DO

        !initialize soil moisture and temperature and Tskin for SCM case
        ELSEIF (l_scm_mode .AND. (atm_phy_nwp_config(jg)%inwp_surface == 1) .AND. &
      &          (i_scm_netcdf == 1)) THEN 

          !IF (i_scm_netcdf==2) THEN !unified format
          ! CALL read_soil_profile_nc_uf(w_so_profile,t_so_profile,t_g_in) 
          !ELSE
           CALL read_soil_profile_nc(w_so_profile,t_so_profile,t_g_in) 
          !ENDIF

          DO jt = 1, ntiles_total
            DO jc = i_startidx, i_endidx
              DO jk = 1,nlev_soil+1 
                p_prog_lnd_now%t_so_t(jc,jk,jb,jt) = t_so_profile(jk)
                p_prog_lnd_new%t_so_t(jc,jk,jb,jt) = t_so_profile(jk)
              END DO
              DO jk = 1,nlev_soil
                ! soil moisture from SCM input needs to be converted from kg/m2 to m
                p_prog_lnd_now%w_so_t(jc,jk,jb,jt) = w_so_profile(jk)/1000._wp 
                p_prog_lnd_new%w_so_t(jc,jk,jb,jt) = w_so_profile(jk)/1000._wp 
              END DO
            END DO
          END DO

          ! set T_G
          DO jc = i_startidx, i_endidx
            p_prog_lnd_now%t_g(jc,jb) = t_g_in
          END DO

          ! cold start when T_SO=0
          IF ( MAXVAL(p_prog_lnd_now%t_so_t) == 0.0_wp ) THEN
            WRITE(*,*) 'init_nwp_phy: perform cold start of T_SO and W_SO because they are 0.0 at input'
            ! set soil temperature to T_G
            DO jt = 1, ntiles_total
              DO jk = 1,nlev_soil+1
                DO jc = i_startidx, i_endidx
                  p_prog_lnd_now%t_so_t(jc,jk,jb,jt) = p_prog_lnd_now%t_g(jc,jb)
                  p_prog_lnd_new%t_so_t(jc,jk,jb,jt) = p_prog_lnd_now%t_g(jc,jb)
                END DO
              END DO
            ! set soil water content to 0.5*(fcap+pwp), use soil type loam (5)
              DO jk = 1,nlev_soil
                DO jc = i_startidx, i_endidx
                  p_prog_lnd_now%w_so_t(jc,jk,jb,jt) = 0.5_wp * (cfcap(5)+cpwp(5)) * dzsoil_icon(jk)
                  p_prog_lnd_new%w_so_t(jc,jk,jb,jt) = 0.5_wp * (cfcap(5)+cpwp(5)) * dzsoil_icon(jk)
                END DO
              END DO
            END DO
            WRITE(*,*) "init_nwp_phy after cold start: w_so ", p_prog_lnd_now%w_so_t(i_startidx,:,jb,:)
            WRITE(*,*) "init_nwp_phy after cold start: t_so ", p_prog_lnd_now%t_so_t(i_startidx,:,jb,:)
          END IF

          ! initialize T skin as top soil layer (as done in mo_initicon_utils.f90)
          DO jt = 1, ntiles_total
            DO jc = i_startidx, i_endidx
              p_prog_lnd_now%t_sk_t(jc,jb,jt) =  p_prog_lnd_now%t_so_t(jc,1,jb,jt)
              p_prog_lnd_new%t_sk_t(jc,jb,jt) =  p_prog_lnd_new%t_so_t(jc,1,jb,jt)
            END DO
          END DO

        ELSEIF ( l_scm_mode .AND. (atm_phy_nwp_config(jg)%inwp_surface == 0) .AND. &
      &          (scm_sfc_temp==1) .AND. (scm_sfc_qv==3) .AND. (i_scm_netcdf > 0) ) THEN 

          CALL read_soil_profile_nc(t_g_in=t_g_in) 
      
      	  ! set T_G
          DO jc = i_startidx, i_endidx
            p_prog_lnd_now%t_g(jc,jb) = t_g_in
            p_prog_lnd_new%t_g(jc,jb) = t_g_in
            p_prog_lnd_now%t_g_t(jc,jb,1) = p_prog_lnd_now%t_g  (jc,jb)
            p_prog_lnd_new%t_g_t(jc,jb,1) = p_prog_lnd_now%t_g  (jc,jb)

            p_diag_lnd%qv_s(jc,jb) = spec_humi( sat_pres_water(t_g_in), p_diag%pres_sfc(jc,jb)) 
            p_diag_lnd%qv_s_t(jc,jb,1) = p_diag_lnd%qv_s(jc,jb)
          END DO

        ELSE  ! any other testcase

          ! t_g  =  t(nlev)
          ! qv_ s= qv(nlev)
          ! KF increase the surface values to obtain fluxes

          DO jc = i_startidx, i_endidx
            p_prog_lnd_now%t_g  (jc,jb)   = p_diag%temp (jc,nlev,jb)!+0.2_wp
            p_prog_lnd_new%t_g  (jc,jb)   = p_diag%temp (jc,nlev,jb)!+0.2_wp
            p_prog_lnd_now%t_g_t(jc,jb,1) = p_prog_lnd_now%t_g  (jc,jb)
            p_prog_lnd_new%t_g_t(jc,jb,1) = p_prog_lnd_now%t_g  (jc,jb)
            ! KF NOTE: as long as we have only water as lower boundary
            ! this is the same setting as for APE
            p_diag_lnd%qv_s    (jc,jb) = &
              & spec_humi(sat_pres_water(p_prog_lnd_now%t_g (jc,jb)),p_diag%pres_sfc(jc,jb))
            p_diag_lnd%qv_s_t(jc,jb,1) = p_diag_lnd%qv_s(jc,jb)
          END DO
        ENDIF ! nh_test_name

        ! Copy t_g to t_seasfc for idealized cases with surface scheme (would be undefined otherwise)
        IF ( atm_phy_nwp_config(jg)%inwp_surface == 1 ) THEN
          DO jc = i_startidx, i_endidx
            p_diag_lnd%t_seasfc(jc,jb) = p_prog_lnd_now%t_g(jc,jb)
          ENDDO
        ENDIF

      ENDDO  ! jb

    ELSE ! For real-case simulations, initialize also qv_s and the tile-based fields

      IF (init_mode /= MODE_IFSANA) THEN
        ! t_g_t and qv_s_t are initialized in read_dwdfg_sfc, calculate the aggregated values
        ! needed for example for initializing the turbulence fields
        CALL aggregate_tg_qvs( p_patch = p_patch,                      & ! in
          &                    frac_t  = ext_data%atm%frac_t(:,:,:),   & ! in
          &                    t_g_t   = p_prog_lnd_now%t_g_t(:,:,:),  & ! in
          &                    qv_s_t  = p_diag_lnd%qv_s_t(:,:,:),     & ! in
          &                    t_g     = p_prog_lnd_now%t_g(:,:),      & ! inout
          &                    qv_s    = p_diag_lnd%qv_s(:,:)          ) ! inout
      ENDIF

      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
             &  i_startidx, i_endidx, rl_start, rl_end)

        IF (init_mode /= MODE_IFSANA) THEN

          DO jc = i_startidx, i_endidx
            p_prog_lnd_new%t_g(jc,jb) =  p_prog_lnd_now%t_g(jc,jb)
          ENDDO

          DO jt = 1, ntiles_total+ntiles_water
            DO jc = i_startidx, i_endidx
              p_prog_lnd_new%t_g_t(jc,jb,jt) = p_prog_lnd_now%t_g_t(jc,jb,jt)
            END DO
          END DO


        ELSE  ! init_mode == MODE_IFSANA

          ! MODE_IFSANA
          ! t_g:
          ! Note, that in copy_prepicon2prog the entire t_g field is initialized with
          ! t_skin. Lake points are re-initialized with MIN(306.15_wp,tskin).
          !
          ! Here, t_g is re-initialized over sea water points with t_seasfc.
          ! Thus:
          ! t_g = tskin (from IFS), for land, lake and seaice points
          ! t_g = t_seasfc for open water
          !
          ! If l_sst_in==FALSE, then t_seasfc=t_skin (with a limiter), so nothing important happens
          !
          ! qv_s:
          ! Over the sea and over the ice, qv_s is set to the saturated value
          ! Over the land we take the minimum of the saturated value and the value
          ! at the first main level above ground
          !

          ! t_g_t, qv_s and qv_s_t are not initialized in case of MODE_IFSANA
!$NEC ivdep
          DO ic=1, ext_data%atm%list_seawtr%ncount(jb)
            jc = ext_data%atm%list_seawtr%idx(ic,jb)
            IF (lseaice) THEN
              ! all points are open water points
              p_prog_lnd_now%t_g(jc,jb) = p_diag_lnd%t_seasfc(jc,jb)
            ELSE
              ! only points with fr_seaice(jc,jb) <= 0.5_wp are open water points and thus
              ! re-initialized with t_seasfc
              IF (p_diag_lnd%fr_seaice(jc,jb) <= 0.5_wp) THEN   ! water point
                p_prog_lnd_now%t_g(jc,jb) = p_diag_lnd%t_seasfc(jc,jb)
              ENDIF
            ENDIF
            p_diag_lnd%qv_s    (jc,jb)    = &
              & spec_humi(sat_pres_water(p_prog_lnd_now%t_g(jc,jb)),p_diag%pres_sfc(jc,jb))
          END DO

          DO ic=1, ext_data%atm%list_seaice%ncount(jb)
            jc = ext_data%atm%list_seaice%idx(ic,jb)
            p_diag_lnd%qv_s    (jc,jb)    = &
              & spec_humi(sat_pres_ice(p_prog_lnd_now%t_g(jc,jb)),p_diag%pres_sfc(jc,jb))
          END DO

          DO ic=1, ext_data%atm%list_lake%ncount(jb)
            jc = ext_data%atm%list_lake%idx(ic,jb)
            ! lake points already initialized in mo_initicon_utils:copy_initicon2prog_sfc
            p_diag_lnd%qv_s    (jc,jb)    = &
              & spec_humi(sat_pres_water(p_prog_lnd_now%t_g(jc,jb)),p_diag%pres_sfc(jc,jb))
          END DO

          DO ic=1, ext_data%atm%list_land%ncount(jb)
            jc = ext_data%atm%list_land%idx(ic,jb)
            p_diag_lnd%qv_s(jc,jb) = &
              &  spec_humi(sat_pres_water(p_prog_lnd_now%t_g (jc,jb)),p_diag%pres_sfc(jc,jb))
            p_diag_lnd%qv_s(jc,jb) = MIN (p_diag_lnd%qv_s(jc,jb), &
              &                    p_prog_now%tracer(jc,nlev,jb,iqv))
          END DO

          DO jc = i_startidx, i_endidx
            p_prog_lnd_new%t_g(jc,jb) = p_prog_lnd_now%t_g(jc,jb)
          ENDDO


          DO jt = 1, ntiles_total+ntiles_water

            DO jc = i_startidx, i_endidx
              p_prog_lnd_now%t_g_t(jc,jb,jt) = p_prog_lnd_now%t_g(jc,jb)
              p_prog_lnd_new%t_g_t(jc,jb,jt) = p_prog_lnd_now%t_g(jc,jb)
              p_diag_lnd%qv_s_t(jc,jb,jt) = p_diag_lnd%qv_s(jc,jb)
            ENDDO
          ENDDO
        END IF  ! init_mode == MODE_IFSANA

      ENDDO  ! jb

    ENDIF  ! ltestcase

    CALL message(modname, 'initialized surface temp and humidity')


  ELSE  ! in case of restart
    !
    ! necessary, because only t_g(nnow_rcf) is written to the restart file
    ! with the following copy statement the ocean points of t_g(nnew_rcf) are
    ! filled with the correct values.

    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        &  i_startidx, i_endidx, rl_start, rl_end)

      DO jc = i_startidx, i_endidx
        p_prog_lnd_new%t_g (jc,jb) = p_prog_lnd_now%t_g (jc,jb)
      ENDDO
      DO jt = 1, ntiles_total+ntiles_water
        DO jc = i_startidx, i_endidx
          p_prog_lnd_new%t_g_t(jc,jb,jt) = p_prog_lnd_now%t_g_t(jc,jb,jt)
        ENDDO
      ENDDO
    ENDDO

    IF (ltestcase .AND. ( nh_test_name == 'RCE' .OR. nh_test_name == 'RCE_Tconst' .OR. nh_test_name == 'RCE_Tprescr' ) .AND. &
        ANY( (/ismag,iprog/) /= atm_phy_nwp_config(jg)%inwp_turb)) THEN !
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &  i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          p_prog_lnd_now%t_g (jc,jb) = th_cbl(1)
          p_prog_lnd_new%t_g (jc,jb) = p_prog_lnd_now%t_g (jc,jb)
          p_diag_lnd%qv_s    (jc,jb) = &
          & spec_humi(sat_pres_water(p_prog_lnd_now%t_g (jc,jb)),p_diag%pres_sfc(jc,jb))
        ENDDO
      ENDDO
    ENDIF

  END IF



  ! index of first half level with height >= 1500 m (above boundary layer)
  k1500m = 1

  !--------------------------------------------------------------
  !>reference pressure according to U.S. standard atmosphere
  ! (with the caveat that the stratosphere is assumed isothermal, which does not hurt
  !  because pref is used for determining model level indices referring to pressures
  !  >= 60 hPa)
  !--------------------------------------------------------------
  ttropo = t00 + dtdz_standardatm*htropo
  ptropo = p0sl_bg*(ttropo/t00)**(-grav/(rd*dtdz_standardatm))
  DO jk = nlev, 1, -1
    jk1 = jk + nshift
    IF ( vct_a(jk1) >= 1500.0_wp .AND. vct_a(jk1+1) < 1500.0_wp ) THEN
      k1500m = jk  ! not jk1 (may result in out-of-bounds)!
    END IF
    zfull = 0.5_wp*(vct_a(jk1) + vct_a(jk1+1))
    IF (zfull < htropo) THEN
      temp = t00 + dtdz_standardatm*zfull
      pref(jk) = p0sl_bg*(temp/t00)**(-grav/(rd*dtdz_standardatm))
    ELSE
      pref(jk) = ptropo*EXP(-grav*(zfull-htropo)/(rd*ttropo))
    ENDIF
  ENDDO

  ! Reference pressure for vertical distribution of aerosol optical depths
  ! Uses coarse-scale orography if available, otherwise use model orography
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

    IF (itype_vegetation_cycle >= 2) THEN
      DO jc = i_startidx,i_endidx
        temp = t00 + dtdz_standardatm*ext_data%atm%topo_t2mclim(jc,jb)
        prm_diag%pref_aerdis(jc,jb) = p0sl_bg*(temp/t00)**(-grav/(rd*dtdz_standardatm))
      ENDDO
    ELSE
      DO jc = i_startidx,i_endidx
        temp = t00 + dtdz_standardatm*ext_data%atm%topography_c(jc,jb)
        prm_diag%pref_aerdis(jc,jb) = p0sl_bg*(temp/t00)**(-grav/(rd*dtdz_standardatm))
      ENDDO
    ENDIF

  ENDDO

  ! start filling phy_params
  !
  !--------------------------------------------------------------
  !< characteristic gridlength needed by convection and turbulence
  !--------------------------------------------------------------
!   CALL sphere_cell_mean_char_length (p_patch%n_patch_cells_g, phy_params%mean_charlen)
  ! read it directly from the patch%geometry_info
  phy_params%mean_charlen = p_patch%geometry_info%mean_characteristic_length
!   write(0,*) "=============================================="
!   write(0,*) "mean_charlen=", phy_params%mean_charlen, &
!     & p_patch%geometry_info%mean_characteristic_length
!   write(0,*) "=============================================="

  ! compute level index corresponding to the HAG of the 60hPa level 
  ! (currently only needed by mo_nwp_diagnosis:cal_cape_cin) 
  phy_params%k060=1
  DO jk=nlev,1,-1
    IF(pref(jk) >  60.e2_wp) phy_params%k060=jk
  ENDDO


  !------------------------------------------
  !< call for cloud microphysics
  !------------------------------------------

  SELECT CASE ( atm_phy_nwp_config(jg)%inwp_gscp )

  CASE (1,2,3)  ! cloud microphysics from COSMO (V 5.0)
    IF (msg_level >= 12)  CALL message(modname, 'init microphysics')

      CALL microphysics_1mom_init( &
        igscp    = atm_phy_nwp_config(jg)%inwp_gscp, &
        tune_zceff_min   = tune_zceff_min,               &
        tune_v0snow      = tune_v0snow,                  &
        tune_zvz0i       = tune_zvz0i,                   &
        tune_icesedi_exp = tune_icesedi_exp,             &
        tune_mu_rain        = atm_phy_nwp_config(1)%mu_rain,&
        tune_rain_n0_factor = atm_phy_nwp_config(1)%rain_n0_factor)
  

  CASE (4,7) !two moment microphysics
    IF (msg_level >= 12)  CALL message(modname, 'init microphysics: two-moment')

    
    ! Provide turbulent length to microphysics to calculate the dissipation factor.
    ! The Prandtl constant (0.4) arises from different definitions for tur_len: microphysics follows Mellor-Yamada.
    atm_phy_nwp_config(jg)%cfg_2mom%turb_len = turbdiff_config(jg)%tur_len * 0.4_wp

    IF (jg == 1) CALL two_moment_mcrph_init(igscp=atm_phy_nwp_config(jg)%inwp_gscp, msg_level=msg_level, &
         &                                  cfg_2mom=atm_phy_nwp_config(jg)%cfg_2mom)

    ! Init of number concentrations moved to mo_initicon_io.f90 !!!

  CASE (8) !sbm micrphysics
    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init microphysics: sbm')

    IF (jg == 1) CALL sbm_init(p_patch, p_prog_now, ext_data%atm%fr_land, p_metrics%ddqz_z_full)

  CASE (5) !two moment microphysics
    IF (msg_level >= 12)  CALL message(modname, 'init microphysics: two-moment')

    IF (jg == 1) CALL two_moment_mcrph_init(atm_phy_nwp_config(jg)%inwp_gscp,&
         &                                  N_cn0,z0_nccn,z1e_nccn,N_in0,z0_nin,z1e_nin,msg_level, &
         &                                  cfg_2mom=atm_phy_nwp_config(jg)%cfg_2mom)

    ! Init of number concentrations moved to mo_initicon_io.f90 !!!

    ! This needs to be coupled to the aerosols (e.g., Tegen climatology)
    IF (linit_mode) THEN ! Initial condition for CCN and IN fields
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jk1,jc,i_startidx,i_endidx,zfull) ICON_OMP_GUIDED_SCHEDULE
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, rl_start, rl_end)
        DO jk=1,nlev
          DO jc=i_startidx,i_endidx
            jk1 = jk + nshift
            zfull = 0.5_wp*(vct_a(jk1)+vct_a(jk1+1))
            IF(zfull > z0_nccn) THEN
              p_prog_now%tracer(jc,jk,jb,inccn) = N_cn0*EXP((z0_nccn-zfull)/z1e_nccn)
            ELSE
              p_prog_now%tracer(jc,jk,jb,inccn) = N_cn0
            END IF
            IF(zfull > z0_nin) THEN
              p_prog_now%tracer(jc,jk,jb,ininpot)  = N_in0*EXP((z0_nin -zfull)/z1e_nin)
            ELSE
              p_prog_now%tracer(jc,jk,jb,ininpot)  = N_in0
            END IF
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    END IF
#ifdef __ICON_ART
  CASE (6) ! two-moment scheme with prognostic cloud droplet number
           ! and chemical composition taken from the ART extension
    IF (msg_level >= 12)  CALL message(modname, 'init microphysics: ART two-moment')
    
    IF (jg == 1) CALL art_clouds_interface_2mom_init(msg_level,cfg_2mom=atm_phy_nwp_config(jg)%cfg_2mom)

    ! Init of number concentrations moved to mo_initicon_io.f90 !!!
#endif
  END SELECT

  ! Fill parameters for cover_koe
  ! Set physics options in cloud cover derived type
  ! This should be moved to init routine in an interface module
  cover_koe_config(jg)%icldscheme  = atm_phy_nwp_config(jg)%inwp_cldcover
  cover_koe_config(jg)%inwp_turb   = atm_phy_nwp_config(jg)%inwp_turb
  cover_koe_config(jg)%inwp_gscp   = atm_phy_nwp_config(jg)%inwp_gscp
  cover_koe_config(jg)%inwp_cpl_re = atm_phy_nwp_config(jg)%icpl_rad_reff
  cover_koe_config(jg)%inwp_reff   = atm_phy_nwp_config(jg)%icalc_reff
  cover_koe_config(jg)%lsgs_cond   = atm_phy_nwp_config(jg)%lsgs_cond
  cover_koe_config(jg)%tune_box_liq_sfc_fac = tune_box_liq_sfc_fac(jg)

#ifdef _OPENACC
  SELECT CASE( cover_koe_config(jg)%icldscheme )
    CASE( 3 )
      CALL finish(routine, "inwp_cldcover=3 is not possible with OpenACC.")
    CASE( 4 )
      CALL finish(routine, "inwp_cldcover=4 is not possible with OpenACC.")
  END SELECT
#endif

  !$ACC ENTER DATA CREATE(cover_koe_config(jg:jg))
  !$ACC UPDATE DEVICE(cover_koe_config(jg:jg)) ASYNC(1) ! This updates all components of cover_koe_config as they are statically allocated

  ! Initiate parameters for reff calculations
  IF (atm_phy_nwp_config(jg)%icalc_reff > 0) THEN
    IF (timers_level > 10) CALL timer_start(timer_phys_reff)
    CALL init_reff ( prm_diag, p_patch, p_prog_now) 
    IF (timers_level > 10) CALL timer_stop(timer_phys_reff)
  END IF

    

  ! Compute lookup tables for aerosol-microphysics coupling
  IF (jg == 1 .AND. (atm_phy_nwp_config(jg)%icpl_aero_gscp > 0 .OR. icpl_aero_conv > 0)) &
    CALL lookupcreate_segalkhain()

  !------------------------------------------
  !< radiation
  !------------------------------------------
  SELECT CASE ( atm_phy_nwp_config(jg)%inwp_radiation )
  CASE (1, 4)

!    prm_diag%lfglac (:,:) = ext_data%atm%soiltyp(:,:) == 1  !soiltyp=ice

    SELECT CASE(atm_phy_nwp_config(jg)%inwp_radiation)
      CASE(1) ! RRTM init
        !
        IF (msg_level >= 12)  CALL message(modname, 'init RRTM')
        !
        CALL setup_srtm
        CALL lrtm_setup(lrtm_filename)
        CALL setup_newcld_optics(cldopt_filename)

        IF ( irad_aero == iRadAeroTegen .OR. irad_aero == iRadAeroART) THEN
          CALL init_aerosol_props_tegen_rrtm
        ELSE
          zaea_rrtm(:,:) = 0.0_wp
          zaes_rrtm(:,:) = 0.0_wp
          zaeg_rrtm(:,:) = 0.0_wp
        ENDIF

        !------------------------------------------------------------
        ! Initialize solar flux in SW bands and solar constant (W/m2)
        !------------------------------------------------------------
        SELECT CASE(isolrad)
          CASE(0)       ! Use ssi values for AMIP-type CMIP5 simulation
            ssi_radt(:) = ssi_amip(:)
          CASE(1)       ! Use ssi values from Coddington et al (2016)
            ssi_radt(:) = ssi_coddington(:)
        END SELECT
        tsi_radt    = SUM(ssi_radt(:))

        ! In case of Aqua planet or RCE experiment:
        IF ( nh_test_name == 'APE_nwp' .OR. nh_test_name == 'dcmip_tc_52' ) THEN
          ssi_radt(:) = ssi_radt(:)*1365._wp/tsi_radt
          tsi_radt    = 1365._wp
        ENDIF  ! APE

        IF ( nh_test_name == 'RCE'         .OR. nh_test_name == 'RCE_Tconst'         .OR. &
           & nh_test_name == 'RCE_Tprescr' .OR. nh_test_name == 'RCEMIP_analytical') THEN
          scale_fac   = sol_const/1361.371_wp   ! computed relative to amip (1361)
          ssi_radt(:) = scale_fac*ssi_amip(:)
          tsi_radt    = SUM(ssi_radt(:))
        ENDIF
        !
      CASE(4) ! ecRad init
#ifdef __ECRAD
        IF (msg_level >= 12)  CALL message(modname, 'init ECRAD')
        !
        ! Do ecrad initialization only once
        IF (.NOT.lreset_mode .AND. jg==1) THEN
          CALL setup_ecrad(ecrad_conf)
          !
          ! Setup Tegen aerosol needs to be done only once for all domains
          IF (irad_aero == iRadAeroTegen .OR. irad_aero == iRadAeroART) THEN
            IF (ecrad_conf%i_gas_model_sw == IGasModelIFSRRTMG .AND. ecrad_conf%i_gas_model_lw == IGasModelIFSRRTMG) THEN
              CALL init_aerosol_props_tegen_ecrad(ecrad_conf, .TRUE.)
            ELSE IF (ecrad_conf%i_gas_model_sw .NE. ecrad_conf%i_gas_model_lw ) THEN
              CALL finish(routine, "Differing gas models for LW and SW are currently unsupported. ")
            ELSE
              CALL init_aerosol_props_tegen_ecrad(ecrad_conf, .FALSE.)
            ENDIF !ecrad_conf%i_gas_model
          ENDIF !irad_aero==iRadAeroTegen .OR. iRadAeroART
        ENDIF ! .NOT.lreset_mode .AND. jg==1
        !
        ! Domain-specific aerosol setups
        IF (ANY( irad_aero == (/iRadAeroConstKinne, iRadAeroKinneVolcSP, iRadAeroKinneSP/) )) THEN
          ! Only the background aerosol (pre-industry) is read in:
          l_filename_year = .FALSE.
          CALL read_bc_aeropt_kinne(ini_date, p_patch, l_filename_year, ecrad_conf%n_bands_lw, ecrad_conf%n_bands_sw)
        ENDIF

        CALL nwp_aerosol_init(ini_date, p_patch)

        IF (ANY( irad_aero == (/iRadAeroKinne,iRadAeroKinneVolc/) )) THEN
          ! Transient Kinne aerosol:
          l_filename_year = .TRUE.
          CALL read_bc_aeropt_kinne(ini_date, p_patch, l_filename_year, ecrad_conf%n_bands_lw, ecrad_conf%n_bands_sw)
        ENDIF
        IF (ANY( irad_aero == (/iRadAeroVolc,iRadAeroKinneVolc,iRadAeroKinneVolcSP/) )) THEN
          ! Volcanic aerosol from CMIP6
          CALL read_bc_aeropt_cmip6_volc(ini_date, ecrad_conf%n_bands_lw, ecrad_conf%n_bands_sw)
        ENDIF
        IF (ANY( irad_aero == (/iRadAeroKinneVolcSP,iRadAeroKinneSP/) )) THEN
          ! Simple plume anthropogenic aerosol
          CALL setup_bc_aeropt_splumes
        ENDIF
        !
        ! Read ozone transient data
        IF (irad_o3 == 5) CALL read_bc_ozone(ini_date%date%year,p_patch,irad_o3,vmr2mmr_opt=o3mr2gg)

        ! cloud_num_fac is used in clim_cdnc, but is only available after the 1st call of init_slowphys
        ! however, clim_cdnc has to be called once before the 1st call of init_slowphys
        IF (atm_phy_nwp_config(jg)%lscale_cdnc .AND. linit_mode) THEN
          prm_diag%cloud_num_fac(:,:) = 1._wp
        ENDIF

        !------------------------------------------------------------
        ! Initialize solar flux in SW bands and solar constant (W/m2)
        !------------------------------------------------------------
        SELECT CASE (isolrad)
          CASE(0)       ! Use default ssi values from ecRad
            ssi_radt(:) = ecrad_ssi_default(:)
          CASE(1)       ! 1: Use ssi values from Coddington et al (2016)
            ssi_radt(:) = ecrad_ssi_coddington(:)
          CASE(2)       ! 2: Use ssi values from external file
            CALL read_bc_solar_irradiance(ini_date%date%year,.TRUE.)
            ssi_radt(:) = 0._wp
        END SELECT
        tsi_radt    = SUM(ssi_radt(:))

        ! In case of Aqua planet or RCE experiment:
        IF ( nh_test_name == 'APE_nwp' .OR. nh_test_name == 'dcmip_tc_52' ) THEN
          ssi_radt(:) = ssi_radt(:)*1365._wp/tsi_radt
          tsi_radt    = 1365._wp
        ENDIF  ! APE
        IF ( nh_test_name == 'RCE'         .OR. nh_test_name == 'RCE_Tconst'         .OR. &
           & nh_test_name == 'RCE_Tprescr' .OR. nh_test_name == 'RCEMIP_analytical') THEN
          scale_fac   = sol_const/1361.371_wp   ! computed relative to amip (1361)
          ssi_radt(:) = scale_fac*ssi_amip(:)
          tsi_radt    = SUM(ssi_radt(:))
        ENDIF
        !
#else
        CALL finish(routine,  &
          &      'atm_phy_nwp_config(jg)%inwp_radiation = 4 needs -D__ECRAD.')
#endif
    END SELECT

    rl_start = 1  ! Initialization should be done for all points
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,jk,zprat,zpres,lland,lglac,zn1,&
!$OMP zn2,zcdnc,dz1,dz2,dz3) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
           &  i_startidx, i_endidx, rl_start, rl_end)

      ! Initialize cloud droplet number concentration (acdnc)
      ! like in mo_echam_phy_init.f90
      DO jk = 1,nlev
        DO jc = i_startidx, i_endidx
          zpres = p0ref * (p_metrics%exner_ref_mc(jc,jk,jb))**(cpd/rd)
          zprat=(MIN(8._wp,80000._wp/zpres))**2

          lland = ext_data%atm%llsm_atm_c(jc,jb)
          lglac = ext_data%atm%soiltyp(jc,jb) == 1
          IF (lland.AND.(.NOT.lglac)) THEN
            zn1= 50._wp
            zn2=220._wp
          ELSE
            zn1= 50._wp
            zn2= 80._wp
          ENDIF
          IF (zpres < 80000._wp) THEN
            zcdnc=1.e6_wp*(zn1+(zn2-zn1)*(EXP(1._wp-zprat)))
          ELSE
            zcdnc=zn2*1.e6_wp
          ENDIF
          prm_diag%acdnc(jc,jk,jb) = zcdnc
        END DO !jc
      END DO   !jk
      IF ( nh_test_name == 'RCE' .OR. nh_test_name == 'RCE_Tconst'   &
         & .OR. nh_test_name == 'RCE_Tprescr'                        ) THEN
        DO jk = 1,nlev
          DO jc = i_startidx, i_endidx
            !--- computation of reference pressure field from the reference exner field
            zrefpres(jc,jk,jb) = p0ref * (p_metrics%exner_ref_mc(jc,jk,jb))**(cpd/rd)
            ! here we choose to use temp to compute sfc pres instead of tempv
            zreftemp(jc,jk,jb) = p_metrics%theta_ref_mc(jc,jk,jb)*p_metrics%exner_ref_mc(jc,jk,jb)
          END DO !jc
        END DO   !jk
      END IF
      IF ( nh_test_name == 'RCE' .OR. nh_test_name == 'RCE_Tconst'   &
         & .OR. nh_test_name == 'RCE_Tprescr'                        ) THEN
        ! a ref press field needs to be computed for testcases with a
        ! constant ozone.  the reference field allows the ozone to be
        ! interpolated at a restart without changing due to a changing p field.
        DO jc = i_startidx,i_endidx
          ! we also need pres at interface levels; first we need sfc press...
          ! Height differences between surface and third-lowest main level
          dz1 = p_metrics%z_ifc(jc,nlev,jb)   - p_metrics%z_ifc(jc,nlevp1,jb)
          dz2 = p_metrics%z_ifc(jc,nlev-1,jb) - p_metrics%z_ifc(jc,nlev,jb)
          dz3 = p_metrics%z_mc (jc,nlev-2,jb) - p_metrics%z_ifc(jc,nlev-1,jb)
          ! Compute surface pressure starting from three lowest levels
          zpres_sfc(jc,jb) = p0ref * EXP( cpd_o_rd*LOG(p_metrics%exner_ref_mc(jc,nlev-2,jb))  + &
                             grav_o_rd*(dz1/zreftemp(jc,nlev,jb) + dz2/zreftemp(jc,nlev-1,jb) + &
                             dz3/zreftemp(jc,nlev-2,jb)) )

          zpres_ifc(jc,nlevp1,jb) = zpres_sfc(jc,jb)
        END DO !jc

        ! compute interface from nlev-1 to TOA
        DO jk = nlev,2,-1
          DO jc = 1, i_endidx
            ! pressure at interface levels
            zpres_ifc(jc,jk,jb) = SQRT(zrefpres(jc,jk,jb)*zrefpres(jc,jk-1,jb) )
          END DO
          DO jc = i_startidx, i_endidx !pres at top ifc = pres at top model lev
            zpres_ifc(jc,1,jb) = zrefpres(jc,1,jb)
          END DO
        END DO
      END IF


    !------------------------------------------
    ! APE ozone profile, vertical setting needed only once for NH
    !------------------------------------------
      !IF (irad_o3 == io3_ape .AND. linit_mode) THEN
      IF (irad_o3 == io3_ape ) THEN

!        CALL o3_zl2ml(p_patch%nblks_c,p_patch%npromz_c,        & !
!          &           nlev_o3,      nlev,                      & ! vertical levels in/out
!          &           zf_aux,   p_metrics%z_mc,                & ! vertical in/out
!          &           ext_data%atm_td%o3(:,:,:,nmonths),p_prog%tracer(:,:,:,io3))! o3Field in/out

        IF ( nh_test_name == 'RCE' .OR. nh_test_name == 'RCE_Tconst' &
           & .OR. nh_test_name == 'RCE_Tprescr'                      ) THEN
          CALL o3_pl2ml (jcs=i_startidx, jce=i_endidx,     &
            & kbdim=nproma,                                &
            & nlev_pres = ext_o3_attr(jg)%nlev_o3,         &
            & klev= nlev ,                                 &
            & pfoz = ext_data%atm_td%pfoz(:),              &
            & phoz = ext_data%atm_td%phoz(:),              &! in o3-levs
            & ppf = zrefpres (:,:,jb),                  &! in  pres
            & pph = zpres_ifc(:,:,jb),               &! in  pres_halfl
            & o3_time_int = ext_data%atm_td%o3(:,:,jb,ext_o3_attr(jg)%nmonths),     &! in
            & o3_clim     = ext_data%atm%o3(:,:,jb) )         ! OUT
        ELSE ! default behaviour
          CALL o3_pl2ml (jcs=i_startidx, jce=i_endidx,     &
            & kbdim=nproma,                                &
            & nlev_pres = ext_o3_attr(jg)%nlev_o3,         &
            & klev= nlev ,                                 &
            & pfoz = ext_data%atm_td%pfoz(:),              &
            & phoz = ext_data%atm_td%phoz(:),              &! in o3-levs
            & ppf = p_diag%pres (:,:,jb),                  &! in  pres
            & pph = p_diag%pres_ifc(:,:,jb),               &! in  pres_halfl
            & o3_time_int = ext_data%atm_td%o3(:,:,jb,ext_o3_attr(jg)%nmonths),     &! in
            & o3_clim     = ext_data%atm%o3(:,:,jb) )         ! OUT
        ENDIF
      ENDIF

    ENDDO      !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    DO ist = 1, UBOUND(csalbw,1)
      rad_csalbw(ist) = csalbw(ist) / (2.0_wp * zml_soil(1))
    ENDDO

  END SELECT !inwp_radiation

  IF ( nh_test_name == 'RCE' .OR. nh_test_name == 'RCE_Tconst' .OR. &
     & nh_test_name == 'RCE_Tprescr'                                ) THEN
    DEALLOCATE (zrefpres)
    DEALLOCATE (zreftemp)
    DEALLOCATE (zpres_sfc)
    DEALLOCATE (zpres_ifc)
  END IF

  !----------------------------------------------------------------
  !< initializations needed both for convection and inwp_cldcover=1
  !----------------------------------------------------------------

  IF ( atm_phy_nwp_config(jg)%inwp_convection == 1 .OR. &
    &  atm_phy_nwp_config(jg)%inwp_cldcover == 1   .OR. &
    &  atm_phy_nwp_config(jg)%inwp_surface == 1 )   THEN

    !This has to be done here because not only convection, but also inwp_cldcover == 1
    !uses mo_cufunctions's foealfa. Therefore, the parameters of the function foealfa
    !have to be initialized by calls of sucst and su_yoethf.

    ! get current date in iso-format "yyyymmddThhmmssZ" (String)
    CALL datetimeToString(ini_date, datetime_string)
    ! convert first 8 characters to integer (yyyy-mm-dd)
    WRITE (yyyymmdd, '(a,a,a)')  datetime_string(1:4), datetime_string(6:7), datetime_string(9:10)
    READ  (yyyymmdd,'(i8)') icur_date

    CALL sucst(54,icur_date,0,0)
    CALL su_yoethf

  ENDIF


  !------------------------------------------
  !< call for convection
  !------------------------------------------

  IF ( atm_phy_nwp_config(jg)%inwp_convection == 1 ) THEN

    IF (msg_level >= 12)  CALL message(modname, 'init convection')

    ! Please take care for scale-dependent initializations!
    ! rsltn = Average mesh size of ICON grid
    ! needed for RTAU - CAPE calculation
    ! adapted for more general gemoetries
    rsltn = p_patch%geometry_info%mean_characteristic_length


!    WRITE(message_text,'(i3,i10,f20.10)') jg, nsmax, phy_params%mean_charlen
!    CALL message('nwp_phy_init, nsmax=', message_text)

    lshallow   = atm_phy_nwp_config(jg)%lshallowconv_only
    lgrayzone_dc = atm_phy_nwp_config(jg)%lgrayzone_deepconv
    ldetrain_prec = atm_phy_nwp_config(jg)%ldetrain_conv_prec
    lrestune_off = atm_phy_nwp_config(jg)%lrestune_off
    lmflimiter_off = atm_phy_nwp_config(jg)%lmflimiter_off
    lstoch_expl = atm_phy_nwp_config(jg)%lstoch_expl
    lstoch_sde = atm_phy_nwp_config(jg)%lstoch_sde
    lstoch_deep = atm_phy_nwp_config(jg)%lstoch_deep
    lvvcouple = atm_phy_nwp_config(jg)%lvvcouple
    lvv_shallow_deep = atm_phy_nwp_config(jg)%lvv_shallow_deep
    
    CALL sucumf(rsltn,nlev,phy_params,lshallow,lgrayzone_dc,ldetrain_prec,lrestune_off, &
         & lmflimiter_off,lstoch_expl,lstoch_sde,lstoch_deep,lvvcouple,lvv_shallow_deep, &
         & pref)
    CALL suphli
    CALL suvdf
    CALL suvdfs
    CALL sucldp

    CALL message(modname, 'convection initialized')
  ELSE
    ! initialize parameters that are accessed outside the convection scheme
    phy_params%rcucov           = 0._wp
    phy_params%rcucov_trop      = 0._wp
    phy_params%rhebc_land       = 0._wp
    phy_params%rhebc_ocean      = 0._wp
    phy_params%rhebc_land_trop  = 0._wp
    phy_params%rhebc_ocean_trop = 0._wp
    phy_params%entrorg          = 0._wp
    phy_params%texc             = 0._wp
    phy_params%qexc             = 0._wp
  ENDIF

  ! Initialize fields k850 and k950, which are required for computing the
  ! convective contribution to wind gusts
  !
  ! k800, k400 will be used for inwp_convection==0 as well. 
  ! k700 is used for LHN data assimilation
  ! Thus we need to make sure that they are initialized.
  prm_diag%k650(:,:) = nlev
  prm_diag%k850(:,:) = nlev
  prm_diag%k950(:,:) = nlev
  prm_diag%k800(:,:) = nlev
  prm_diag%k400(:,:) = nlev
  prm_diag%k700(:,:) = nlev

  rl_start = 1  ! Initialization should be done for all points
  rl_end   = min_rlcell

  i_startblk = p_patch%cells%start_blk(rl_start,1)
  i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

  ! height of 850 and 950hPa surface for US standard atmosphere in m
  ! For derivation, see documentation of US standard atmosphere
  h650_standard = 3590.69_wp
  h850_standard = 1457.235199_wp
  h950_standard = 540.3130233_wp

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx,hag,zpres,zpres0) ICON_OMP_DEFAULT_SCHEDULE

  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
      &  i_startidx, i_endidx, rl_start, rl_end)

    DO jc=i_startidx, i_endidx

      DO jk=nlev, 1, -1
        ! height above ground
        hag = p_metrics%z_mc(jc,jk,jb)-ext_data%atm%topography_c(jc,jb)

        IF (hag < h650_standard) THEN
          prm_diag%k650(jc,jb) = jk
        ENDIF
        IF (hag < h950_standard) THEN
          prm_diag%k950(jc,jb) = jk
        ENDIF
        IF (hag < h850_standard) THEN
          prm_diag%k850(jc,jb) = jk
        ELSE
          EXIT
        ENDIF
      ENDDO
      ! security measure
      prm_diag%k650(jc,jb) = MAX(prm_diag%k650(jc,jb),2)
      prm_diag%k950(jc,jb) = MAX(prm_diag%k950(jc,jb),2)
      prm_diag%k850(jc,jb) = MAX(prm_diag%k850(jc,jb),2)

      ! analogous initialization of k800 and k400, based on reference pressure
      ! because this is more meaningful for k400 in the presence of very high orography
      zpres0 = p0ref * (p_metrics%exner_ref_mc(jc,nlev,jb))**(cpd/rd)
      DO jk=nlev-1, 2, -1
        zpres = p0ref * (p_metrics%exner_ref_mc(jc,jk,jb))**(cpd/rd)
        IF (zpres/zpres0 >= pr800) prm_diag%k800(jc,jb) = jk
        IF (zpres/zpres0 >= pr700) prm_diag%k700(jc,jb) = jk
        IF (zpres/zpres0 >= pr400*SQRT(p0ref/zpres0)) THEN
          prm_diag%k400(jc,jb) = jk
        ELSE
          EXIT
        ENDIF
      ENDDO
   ENDDO  ! jc

  ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  !------------------------------------------
  !< surface initialization (including seaice)
  !------------------------------------------

  SELECT CASE (atm_phy_nwp_config(jg)%inwp_surface)
  CASE (LSS_TERRA)
    IF (linit_mode) THEN
      CALL nwp_surface_init(p_patch, ext_data, p_prog_lnd_now, p_prog_lnd_new, &
        &                   p_prog_wtr_now, p_prog_wtr_new, p_diag_lnd, p_diag, prm_diag)
    ELSE
      IF ( lsnowtile ) THEN ! restart mode with snowtiles
        CALL init_snowtile_lists(p_patch, ext_data, p_diag_lnd)
      ENDIF
    ENDIF

    ! Copy timelevel now to timelevel new for land state. This has no impact on the prognostic
    ! results but ensures that the output over non-prognostic grid points (water) is the same
    ! for even and odd multiples of the advection time step
    CALL copy_lnd_prog_now2new(p_patch, p_prog_lnd_now, p_prog_lnd_new)
    IF (msg_level >= 12)  CALL message(modname, 'init TERRA')
  CASE (LSS_JSBACH)
    CALL nwp_vdiff_init( &
        & p_patch, vdiff_config(jg), prm_diag, ext_data, p_prog_lnd_now, p_diag_lnd, &
        & prm_diag%nwp_vdiff_state &
      )
  END SELECT


  !------------------------------------------
  !< setup for turbulence
  !------------------------------------------

  ! initialize gz0 (roughness length * g)
  !
  IF ( ANY( (/icosmo,igme,ismag,iprog/)==atm_phy_nwp_config(jg)%inwp_turb ) .AND. &
       linit_mode ) THEN !seclected turbulence scheme requires '%gz0'-initialization

    ! gz0 is initialized if we do not start from an own first guess
    IF (lturb_init) THEN !no first guess present

      IF (msg_level >= 12)  CALL message(modname, 'init roughness length')

      IF (turbdiff_config(jg)%lconst_z0) THEN
        ! constant z0 for idealized tests
        prm_diag%gz0(:,:) = grav * turbdiff_config(jg)%const_z0
      ELSE IF (atm_phy_nwp_config(jg)%itype_z0 == 1) THEN
        ! default
        prm_diag%gz0(:,:) = grav * ext_data%atm%z0(:,:)

      ELSE IF (atm_phy_nwp_config(jg)%itype_z0 >= 2) THEN !use tile-specific information

        rl_start = grf_bdywidth_c + 1 ! land-cover classes are not set for nest-boundary points
        rl_end   = min_rlcell_int

        i_startblk = p_patch%cells%start_blk(rl_start,1)
        i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)


        ! Scaling factor for SSO contribution to roughness length ("Erdmann Heise formula")
        fact_z0rough = 1.e-5_wp*ATAN(phy_params%mean_charlen/2250._wp)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,ic,jt,i_startidx,i_endidx,lc_class) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
            &                i_startidx, i_endidx, rl_start, rl_end)

          ! specify land-cover-related roughness length over land points:

!$NEC ivdep
          DO ic = 1, ext_data%atm%list_land%ncount(jb)
            jc = ext_data%atm%list_land%idx(ic,jb)
            prm_diag%gz0(jc,jb) = 0._wp
          END DO

          !Notes:
          !In '%gz0', snow-cover dependent tile-contributions of provided tile-specific z0-values
          ! are added together for each grid-point being considered as a land-point
          ! (defined by "%fr_land>frlnd_thrhld", s. SUB 'init_index_lists').
          !The resulting grid-point mean of '%gz0' is going to be used below for the initialization calls of
          ! 'turbtran' and 'turbdiff'.
          !For non-land points (water- or ice-surfaces), '%gz0' is initialized by SUB 'turbtran'.

          DO jt = 1, ntiles_total
!$NEC ivdep
            DO ic = 1, ext_data%atm%list_land%ncount(jb)
              jc = ext_data%atm%list_land%idx(ic,jb)
              lc_class = MAX(1,ext_data%atm%lc_class_t(jc,jb,jt)) ! to avoid segfaults
              prm_diag%gz0(jc,jb) = prm_diag%gz0(jc,jb) + ext_data%atm%frac_t(jc,jb,jt)*grav*( &
                  (1._wp-p_diag_lnd%snowfrac_t(jc,jb,jt))*ext_data%atm%z0_lcc(lc_class) +      &
                         p_diag_lnd%snowfrac_t(jc,jb,jt) *0.5_wp*ext_data%atm%z0_lcc(i_lc_si) ) !i_lc_si: snow/ice class
            ENDDO
          ENDDO
          IF (atm_phy_nwp_config(jg)%itype_z0 == 3) THEN
!$NEC ivdep
            DO ic = 1, ext_data%atm%list_land%ncount(jb)
              jc = ext_data%atm%list_land%idx(ic,jb)
              prm_diag%gz0(jc,jb) = prm_diag%gz0(jc,jb) + &
                                    grav*MIN(fact_z0rough*ext_data%atm%sso_stdh_raw(jc,jb)**2,7.5_wp)
            ENDDO
          ENDIF
          DO jt = ntiles_total+1, ntiles_total+ntiles_water ! required if there are mixed land-water points
!$NEC ivdep
            DO ic = 1, ext_data%atm%list_land%ncount(jb)
              jc = ext_data%atm%list_land%idx(ic,jb)
              lc_class = MAX(1,ext_data%atm%lc_class_t(jc,jb,jt)) ! to avoid segfaults
              prm_diag%gz0(jc,jb) = prm_diag%gz0(jc,jb) + &
                                    ext_data%atm%frac_t(jc,jb,jt) * grav*ext_data%atm%z0_lcc(lc_class)
            ENDDO
          ENDDO
        ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

      ENDIF !use tile-specific information

    END IF !no first guess present
  ENDIF !seclected turbulence scheme requires '%gz0'-initialization

  IF ( ANY( (/icosmo/)==atm_phy_nwp_config(jg)%inwp_turb ) ) THEN

    ! allocate and init implicit weights for tridiagonal solver
    ALLOCATE( turbdiff_config(jg)%impl_weight(nlevp1), STAT=istatus )
    ! note that impl_weight => turbdiff_config(jg)%impl_weight
    !$ACC ENTER DATA CREATE(turbdiff_config(jg:jg))
    !$ACC ENTER DATA CREATE(turbdiff_config(jg)%impl_weight)
    IF(istatus/=SUCCESS)THEN
      CALL finish(routine, &
                 'allocation of impl_weight failed')
    ENDIF

    CALL get_turbdiff_param(jg)

    ! using an over implicit value (impl_s) near surface,
    ! reduced to in general slightly off-centered value (impl_t)
    ! in about 1500 m height
    DO jk = 1, k1500m
      impl_weight(jk) = impl_t
    END DO
    DO jk = k1500m+1, nlev
      impl_weight(jk) = impl_t &
                      + (impl_s-impl_t) * (jk-k1500m) / REAL(nlev-k1500m, wp)
    END DO
    impl_weight(nlevp1) = impl_s
    ! impl_weight is never changed
    !$ACC UPDATE DEVICE(turbdiff_config(jg)%impl_weight) ASYNC(1)

! computing l_pat: cannot be done in mo_ext_data_init, because it seems we do not have
!                  phy_params(jg)%mean_charlen available then?

    rl_start = 1 ! Initialization is done also for nest boundary points
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)

      ! effective length scale of circulation patterns
      DO jc = i_startidx, i_endidx
        IF (ext_data%atm%fr_land(jc,jb)<=0.5_wp) THEN
          ext_data%atm%l_pat (jc,jb) = 0.0_wp
        ELSE !a point being dominated by a land surface
          IF (imode_pat_len == 2) THEN !using the standard deviat. of SGS orography as a lower limit 
                                       ! of 'l_pat'
            ext_data%atm%l_pat (jc,jb) = MIN( pat_len, prm_diag%pat_len(jc,jb)          )
          ELSE !using a constant 'l_pat'
            ext_data%atm%l_pat (jc,jb) = pat_len
          ENDIF
          ext_data%atm%l_pat (jc,jb) = (phy_params%mean_charlen * ext_data%atm%l_pat (jc,jb)) / &
             &                         (phy_params%mean_charlen + ext_data%atm%l_pat (jc,jb))
        ENDIF
      ENDDO
    ENDDO

  ENDIF

  ! Initialize turbulence models
  !
  IF ( ( ANY( (/icosmo/)==atm_phy_nwp_config(jg)%inwp_turb ) ) .AND. linit_mode ) THEN

    IF (msg_level >= 12)  CALL message(modname, 'init COSMO turbulence')

    rl_start = 1 ! Initialization is done also for nest boundary points
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)


#ifndef __PGI
!FIXME: PGI + OpenMP produce deadlock in this loop. Compiler bug suspected
!$OMP PARALLEL DO PRIVATE(jb,jk,i_startidx,i_endidx,ic,jc,jt, &
!$OMP            ltkeinp_loc,lgz0inp_loc,nlevcm,l_hori,nzprv,zvari,zrhon, &
!$OMP            l_lake,l_sice, &
!$OMP            ierrstat, errormsg, eroutine) ICON_OMP_DEFAULT_SCHEDULE
#endif

    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)

      IF (.not. ltestcase) THEN
       IF (lturb_init) THEN

        ltkeinp_loc = .FALSE.  ! initialize TKE field
        lgz0inp_loc = .FALSE.  ! initialize gz0 field (water points only)

       ELSE
        !
        ! TKE and gz0 are not re-initialized, but re-used from the first guess
        !
        ltkeinp_loc = .TRUE.   ! do NOT re-initialize TKE field (read from FG)
        lgz0inp_loc = .TRUE.   ! do NOT re-initialize gz0 field (read from FG)

        ! Note that TKE in turbtran/turbdiff is defined as the turbulence velocity scale
        ! TVS=SQRT(2*TKE). The TKE is limited to 5.e-5 here because it may be zero on lateral
        ! boundary points for the limited-area mode, which would cause a crash in the initialization
        ! performed here but hs no impact on the results otherwise.
        !
        DO jk =1,nlevp1
          DO jc = i_startidx, i_endidx
            p_prog_now%tke(jc,jk,jb)= SQRT(2.0_wp*MAX(5.e-5_wp,p_prog_now%tke(jc,jk,jb)))
          ENDDO
        ENDDO
       ENDIF
      ELSE !ltestcase

       IF (lscm_read_tke) THEN
        ltkeinp_loc = .TRUE.   ! do NOT re-initialize TKE field (read from FG)
        ! Note that TKE in turbtran/turbdiff is defined as the turbulence velocity scale
        ! TVS=SQRT(2*TKE). The TKE is limited to 5.e-5 here because it may be zero on lateral
        ! boundary points for the limited-area mode, which would cause a crash in the initialization
        ! performed here but hs no impact on the results otherwise.
        !
        DO jk =1,nlevp1
          DO jc = i_startidx, i_endidx
            p_prog_now%tke(jc,jk,jb)= SQRT(2.0_wp*MAX(5.e-5_wp,p_prog_now%tke(jc,jk,jb)))
          ENDDO
        ENDDO
       ELSE
        ltkeinp_loc = .FALSE.  ! initialize TKE field
       ENDIF
       IF (lscm_read_z0) THEN
        lgz0inp_loc = .TRUE.   ! do NOT re-initialize gz0 field (read from FG)
       ELSE
        lgz0inp_loc = .FALSE.  ! initialize gz0 field (water points only)
       ENDIF

      ENDIF

      l_hori(i_startidx:i_endidx)=phy_params%mean_charlen

      DO jc = i_startidx, i_endidx
        l_lake(jc)=(ext_data%atm%depth_lk(jc,jb)>0._wp)
         l_sice(jc) = ext_data%atm%fr_land(jc,jb)<=0.5_wp .AND. & ! a non-land point, and:
                      MERGE( MERGE( (p_prog_wtr_now%h_ice(jc,jb)>=h_Ice_min_flk), & !either more than marginal lake-ice depth from lake-scheme
                                    (p_prog_lnd_now%t_g  (jc,jb)< tpl_T_f      ), & ! or frozen lake due to surface-temperature,
                                    llake   ), & !dependent on whether lake scheme is active or not
                             MERGE( (p_prog_wtr_now%h_ice(jc,jb)>0._wp  ), & !either sea-ice present from sea-ice scheme
                                    (p_prog_lnd_now%t_g  (jc,jb)<tf_salt), & !or frozen salty sea due to surface temperature
                                    lseaice ), & !dependent on whether sea-ice scheme is active or not
                             l_lake(jc) ) !seperate treatment of lake points and non-lake points
      END DO

      !Note(MR):
      !Since there should be no '%h_ice' at non-land points, 'l_sice' indicates all grid-points coverd by sea- or lake-ice,
      ! while frozen land points e.g. covered by snow or glaciers) are excluded.
      !This is additionally ensured by only considering points with "%fr_land(jc,jb)<=0.5p".

      nlevcm = 3 !so far no vertically resolved roughness layer

      ! set nzprv to 1 and use it as nprv (which is INOUT in turbtran and turbdiff
      !  but this feature is not used here during initialization)
      nzprv  = 1

      ! turbtran: only surface-layer turbulence
      CALL turbtran (                                                          & !only surface-layer turbulence
!
        &  iini=1, ltkeinp=ltkeinp_loc, lgz0inp=lgz0inp_loc, lstfnct=.TRUE. ,  &
        &          lsrflux=.TRUE., lnsfdia=.TRUE., lrunscm=.FALSE.,            & !incl. near-surf. diagn. and surf.-flux calcul.
        &          ladsshr=.FALSE.,                                            & !no additional NTC-shear
!
        &  dt_tke=atm_phy_nwp_config(jg)%dt_fastphy  ,                         &
        &  nprv=nzprv, ntur=1, ntim=1,                                         &
        &  nvec=nproma, ke=2, ke1=3, kcm=nlevcm, iblock=jb,                    &
        &  ivstart=i_startidx, ivend=i_endidx,                                 &
!
        &  l_hori=l_hori, hhl=p_metrics%z_ifc(:,nlev-1:nlevp1,jb),             &
        &  l_pat=ext_data%atm%l_pat(:,jb),                                     &
        &  fr_land=ext_data%atm%fr_land(:,jb),                                 &
        &  l_lake=l_lake(:),                                                   & !in (lake surfaces)
        &  l_sice=l_sice(:),                                                   & !in (ice  surfaces)
!
        &  rlamh_fac=prm_diag%rlamh_fac_t(:,jb,1),                             &
        &  sai=ext_data%atm%sai(:,jb),                                         &
        &  urb_isa=ext_data%atm%urb_isa_t(:,jb,1),                             &
        &  gz0=prm_diag%gz0(:,jb),                                             &
!
        &  t_g=p_prog_lnd_now%t_g(:,jb),                                       &
        &  qv_s=p_diag_lnd%qv_s(:,jb),                                         &
        &  ps=p_diag%pres_sfc(:,jb),                                           &
!
        &  u=p_diag%u(:,nlev-1:nlev,jb),                                       &
        &  v=p_diag%v(:,nlev-1:nlev,jb),                                       &
        &  t=p_diag%temp(:,nlev-1:nlev,jb),                                    &
        &  qv=p_prog_now%tracer(:,nlev-1:nlev,jb,iqv),                         &
        &  qc=p_prog_now%tracer(:,nlev-1:nlev,jb,iqc),                         &
        &  epr=p_prog_now%exner(:,nlev-1:nlev,jb),                             &
!
        &  tcm=prm_diag%tcm(:,jb),                                             &
        &  tch=prm_diag%tch(:,jb),                                             &
        &  tvm=prm_diag%tvm(:,jb),                                             &
        &  tvh=prm_diag%tvh(:,jb),                                             &
        &  tfm=prm_diag%tfm(:,jb),                                             &
        &  tfh=prm_diag%tfh(:,jb),                                             &
        &  tfv=prm_diag%tfv(:,jb),                                             &
        &  tkr=prm_diag%tkr_t(:,jb,1),                                         &
        &  tke=p_prog_now%tke(:,nlev-1:nlevp1,jb),                             &
        &  tkvm=prm_diag%tkvm(:,nlev-1:nlevp1,jb),                             &
        &  tkvh=prm_diag%tkvh(:,nlev-1:nlevp1,jb),                             &
        &  rcld=prm_diag%rcld(:,nlev-1:nlevp1,jb),                             &
!
        &  t_2m=prm_diag%t_2m(:,jb),                                           &
        &  qv_2m=prm_diag%qv_2m(:,jb),                                         &
        &  td_2m=prm_diag%td_2m(:,jb),                                         &
        &  rh_2m=prm_diag%rh_2m(:,jb),                                         &
        &  u_10m=prm_diag%u_10m(:,jb),                                         &
        &  v_10m=prm_diag%v_10m(:,jb),                                         &
!
        &  shfl_s=prm_diag%shfl_s(:,jb),                                       &
        &  qvfl_s=prm_diag%qhfl_s(:,jb)                                        &
        &                                                                      ) !end of 'turbtran' call

      prm_diag%lhfl_s(i_startidx:i_endidx,jb) = prm_diag%qhfl_s(i_startidx:i_endidx,jb) * lh_v

      nlevcm = nlevp1 !so far no vertically resolved roughness layer

      CALL turbdiff ( &
!
        &  iini=1,                                                  & ! separate initialization before the time loop
        &  ltkeinp=ltkeinp_loc,                                     & ! 
        &  lstfnct=.TRUE. ,                                         & !
        &  l3dturb=.FALSE.,                                         & ! not yet arranged for ICON
        &  lrunsso=(atm_phy_nwp_config(jg)%inwp_sso > 0),           & ! running COSMO SSO scheme
        &  lruncnv=(atm_phy_nwp_config(jg)%inwp_convection > 0),    & ! running convection
        &  lrunscm=.FALSE.,                                         & ! no single column model
!
        &  dt_var=atm_phy_nwp_config(jg)%dt_fastphy,                &
        &  dt_tke=atm_phy_nwp_config(jg)%dt_fastphy,                &
        &  nprv=nzprv, ntur=1, ntim=1,                              &
        &  nvec=nproma, ke=nlev, ke1=nlevp1, kcm=nlevcm, iblock=jb, &
        &  ivstart=i_startidx, ivend=i_endidx,                      &
!
        &  l_hori=l_hori,                                           &
        &  hhl=p_metrics%z_ifc(:,:,jb),                             &
        &  dp0=p_diag%dpres_mc(:,:,jb),                             &
        &  gz0=prm_diag%gz0(:,jb),                                  &
        &  l_pat = ext_data%atm%l_pat(:,jb),                        &
!
        &  t_g=p_prog_lnd_now%t_g(:,jb),                            &
        &  qv_s=p_diag_lnd%qv_s(:,jb),                              &
        &  ps=p_diag%pres_sfc(:,jb),                                &
!
        &  u=p_diag%u(:,:,jb),                                      &
        &  v=p_diag%v(:,:,jb),                                      &
        &  w=p_prog_now%w(:,:,jb),                                  &
        &  t=p_diag%temp(:,:,jb),                                   &
        &  qv=p_prog_now%tracer(:,:,jb,iqv),                        &
        &  qc=p_prog_now%tracer(:,:,jb,iqc),                        &
        &  prs=p_diag%pres(:,:,jb),                                 &
        &  rhoh=p_prog_now%rho(:,:,jb),                             &
        &  rhon=zrhon(:,:),                                         &
        &  epr=p_prog_now%exner(:,:,jb),                            &
!
        &  impl_weight=turbdiff_config(jg)%impl_weight,             &
        &  tvm=prm_diag%tvm(:,jb),                                  &
        &  tvh=prm_diag%tvh(:,jb),                                  &
        &  tfm=prm_diag%tfm(:,jb),                                  &
        &  tfh=prm_diag%tfh(:,jb),                                  &
        &  tfv=prm_diag%tfv(:,jb),                                  &
        &  tke=p_prog_now%tke(:,:,jb),                              &
        &  tkvm=prm_diag%tkvm(:,:,jb),                              &
        &  tkvh=prm_diag%tkvh(:,:,jb),                              &
        &  tprn=prm_diag%tprn(:,:,jb),                              &
        &  rcld=prm_diag%rcld(:,:,jb),                              &
!
        &  trop_mask=prm_diag%tropics_mask(:,jb),                   & !in
        &  innertrop_mask=prm_diag%innertropics_mask(:,jb),         & !in
!
        &  tketens=prm_nwp_tend%ddt_tke(:,:,jb),                    &
        &  ut_sso=REAL(prm_nwp_tend%ddt_u_sso(:,:,jb),wp),          &
        &  vt_sso=REAL(prm_nwp_tend%ddt_v_sso(:,:,jb),wp),          &
        &  u_tens=prm_nwp_tend%ddt_u_turb(:,:,jb),                  &
        &  v_tens=prm_nwp_tend%ddt_v_turb(:,:,jb),                  &
!
        &  zvari=zvari                                              & !out
        &                                                           ) !end of 'turbdiff' call

      ! preparation for concentration boundary condition. Usually inactive for standard ICON runs.
      IF ( .NOT. lsflcnd ) THEN
        prm_diag%lhfl_s(i_startidx:i_endidx,jb) = prm_diag%qhfl_s(i_startidx:i_endidx,jb) * lh_v
      END IF

      ! tile-specific quantities needed by turbtran
      !
      DO jt = 1, ntiles_total+ntiles_water

        IF (jt<=ntiles_total) THEN !for a not-tiled model run or for a land (sub)tile
           prm_diag%gz0_t(:,jb,jt) = prm_diag%gz0(:,jb) !use '%gz0' for the tile-specific '%gz0_t' at each point
        ELSE !for non-land tiles
           DO jc = i_startidx, i_endidx !for all points
             IF (ext_data%atm%fr_land(jc,jb) <= 0.5_wp ) THEN !'%gz0' has been calculated by 'turbtran'
               prm_diag%gz0_t(jc,jb,jt) = prm_diag%gz0(jc,jb) !use this '%gz0' for the tile-specific '%gz0_t'

               !Note:
               !In this case, '%gz0' already refers to a non-land surface and is a good estimate for '%gz0_t'.
             ELSE !'%gz0' has not been calculated by 'turbtran'
               prm_diag%gz0_t(jc,jb,jt) = grav*z0_ice !use the value for sea-ice as applied in 'turbtran'
             END IF
           END DO
        END IF

        prm_diag%tvs_s_t (:,jb,jt) = p_prog_now%tke(:,nlevp1,jb)  !here: SQRT(2*TKE)
        prm_diag%tkvm_s_t(:,jb,jt) = prm_diag%tkvm(:,nlevp1,jb)
        prm_diag%tkvh_s_t(:,jb,jt) = prm_diag%tkvh(:,nlevp1,jb)
        prm_diag%rcld_s_t(:,jb,jt) = prm_diag%rcld(:,nlevp1,jb)
        prm_diag%tvm_t   (:,jb,jt) = prm_diag%tvm(:,jb)
        prm_diag%tvh_t   (:,jb,jt) = prm_diag%tvh(:,jb)
        IF (jt>1) prm_diag%tkr_t   (:,jb,jt) = prm_diag%tkr_t(:,jb,1) !necessary if "imode_trancnf>=4"
      ENDDO


      ! Note that TKE in turbtran/turbdiff is defined as the turbulence velocity scale
      ! TVS=SQRT(2*TKE)
      !
      DO jk =1,nlevp1
        p_prog_now%tke(i_startidx:i_endidx,jk,jb)= 0.5_wp                        &
          &                                * (p_prog_now%tke(i_startidx:i_endidx,jk,jb))**2
      ENDDO

    ENDDO  ! jb

    IF (msg_level >= 12)  CALL message(modname, 'Cosmo turbulence initialized')

  ELSE IF (  atm_phy_nwp_config(jg)%inwp_turb == igme .AND. linit_mode) THEN

    IF (msg_level >= 12)  CALL message(modname, 'init GME turbulence')

    rl_start = grf_bdywidth_c + 1 ! land-cover classes are not set for nest-boundary points
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jt,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)

      ! paranoia: Make sure that rcld is initialized  (needed by cloud cover scheme)
      prm_diag%rcld(:,:,jb) = 0._wp

      DO jt = 1, ntiles_total+ntiles_water
        prm_diag%gz0_t(i_startidx:i_endidx,jb,jt) = prm_diag%gz0(i_startidx:i_endidx,jb)
      ENDDO
      !Note:
      !Although only '%gz0' (as a grid-surface mean) is required for the turbulence part,
      ! the tile-specific roughnes-length '%gz0_t' is still used in the soil/surfce part!
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  ELSE IF ( ANY( (/ismag,iprog/) == atm_phy_nwp_config(jg)%inwp_turb) .AND. linit_mode ) THEN

    CALL message(modname, 'init LES turbulence')

    IF(atm_phy_nwp_config(jg)%inwp_surface == 0)THEN
      IF (turbdiff_config(jg)%lconst_z0) THEN
        prm_diag%gz0_t(:,:,1) = grav * turbdiff_config(jg)%const_z0
      ELSE
        CALL finish(routine, 'Only constant roughness length allowed idealized LES cases!')
      END IF
    ELSE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, rl_start, rl_end)

        DO jt = 1, ntiles_total+ntiles_water
          prm_diag%gz0_t(i_startidx:i_endidx,jb,jt) = prm_diag%gz0(i_startidx:i_endidx,jb)
        ENDDO
        !Note:
        !Although only '%gz0' (as a grid-surface mean) is required for the turbulence part,
        ! the tile-specific roughnes-length '%gz0_t' is still used in the soil/surfce part!
      ENDDO
    END IF

  END IF
  ! No handling of VDIFF here. VDIFF does its initialization together with the slow physics.
  ! Configured in event setup in src/configure_model/mo_atm_phy_nwp_config.f90.

  ! Gravity wave drag scheme
  !
  IF ( atm_phy_nwp_config(jg)%inwp_gwd == 1 ) THEN  ! IFS gwd scheme

    CALL sugwwms(nflevg=nlev, ppref=pref, klaunch=phy_params%klaunch)
    CALL message(modname, 'non-orog GWs initialized')

  END IF

  ! SSO scheme
  !
  CALL sugwd(nlev, pref, phy_params, jg)
  IF (linit_mode) prm_diag%ktop_envel(:,:) = nlev

  ! read time-dependent boundary conditions from file
  ! well mixed greenhouse gases, horizontally constant
  IF(ANY((/irad_co2,irad_cfc11,irad_cfc12,irad_n2o,irad_ch4/) == 4)) THEN
    ! read annual means
    CALL read_bc_greenhouse_gases(ghg_filename)
    CALL bc_greenhouse_gases_time_interpolation(ini_date, print_report=(msg_level >= 5))
    ! interpolation to the current date and time takes place in the radiation interface
  ENDIF

#ifndef __NO_ICON_UPATMO__
  ! Upper-atmosphere physics
  !
  IF (lupatmo_phy) THEN
    IF (upatmo_config(jg)%l_status( iUpatmoStat%timer )) CALL timer_start(timer_upatmo)
    CALL init_upatmo_phy_nwp( mtime_datetime    = ini_date,          & !in
      &                       p_patch           = p_patch,           & !in
      &                       p_metrics         = p_metrics,         & !in
      &                       p_prog            = p_prog_now,        & !in
      &                       p_diag            = p_diag,            & !in
      &                       prm_nwp_diag      = prm_diag,          & !in
      &                       prm_nwp_tend      = prm_nwp_tend,      & !in
      &                       nproma            = nproma             ) !in
    IF (upatmo_config(jg)%l_status( iUpatmoStat%timer )) CALL timer_stop(timer_upatmo)
  ENDIF
#endif

  ! SPPT
  IF (linit_mode) THEN

    IF (sppt_config(jg)%lsppt) THEN

      ! Initate, i.e. generate random patterns during initiation
      CALL init_rn(p_patch, ini_date, sppt_config(jg), &
        &          sppt(jg)%rn_2d_now, sppt(jg)%rn_2d_new)

      CALL message(modname, 'Initialisation of SPPT completed')

    ENDIF

  ENDIF


  IF (timers_level > 3) CALL timer_stop(timer_init_nwp_phy)

END SUBROUTINE init_nwp_phy


  SUBROUTINE init_cloud_aero_cpl ( mtime_date, p_patch, p_metrics, ext_data, prm_diag)

    TYPE(datetime),   POINTER               :: mtime_date
    TYPE(t_patch),               INTENT(in) :: p_patch
    TYPE(t_nh_metrics),          INTENT(in) :: p_metrics
    TYPE(t_external_data),       INTENT(in) :: ext_data

    TYPE(t_nwp_phy_diag),        INTENT(inout) :: prm_diag

    INTEGER          :: imo1, imo2
    INTEGER          :: rl_start, rl_end, i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER          :: jb, jc, jg, nlev

    REAL(wp) :: wgt, zncn(nproma, p_patch%nlev)
    
    TYPE(t_time_interpolation_weights) :: current_time_interpolation_weights

    TYPE(datetime), POINTER :: mtime_hour
    
    jg = p_patch%id
    nlev = p_patch%nlev

    IF (ALL (irad_aero /= (/iRadAeroTegen, iRadAeroART, iRadAeroCAMSclim, iRadAeroCAMStd/))) RETURN
    IF (atm_phy_nwp_config(jg)%icpl_aero_gscp /= 1 .AND. icpl_aero_conv /= 1) RETURN

    
    mtime_hour => newDatetime(mtime_date)
    mtime_hour%time%minute = 0
    mtime_hour%time%second = 0
    mtime_hour%time%ms     = 0          
    current_time_interpolation_weights = calculate_time_interpolation_weights(mtime_hour)
    call deallocateDatetime(mtime_hour)
    imo1 = current_time_interpolation_weights%month1
    imo2 = current_time_interpolation_weights%month2
    wgt = current_time_interpolation_weights%weight2

    rl_start = 1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,zncn)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
 
      IF (iprog_aero == 0) THEN
        DO jc = i_startidx, i_endidx

          prm_diag%aerosol(jc,iss,jb) = ext_data%atm_td%aer_ss(jc,jb,imo1) + &
            ( ext_data%atm_td%aer_ss(jc,jb,imo2)   - ext_data%atm_td%aer_ss(jc,jb,imo1)   ) * wgt
          prm_diag%aerosol(jc,iorg,jb) = ext_data%atm_td%aer_org(jc,jb,imo1) + &
            ( ext_data%atm_td%aer_org(jc,jb,imo2)  - ext_data%atm_td%aer_org(jc,jb,imo1)  ) * wgt
          prm_diag%aerosol(jc,ibc,jb) = ext_data%atm_td%aer_bc(jc,jb,imo1) + &
            ( ext_data%atm_td%aer_bc(jc,jb,imo2)   - ext_data%atm_td%aer_bc(jc,jb,imo1)   ) * wgt
          prm_diag%aerosol(jc,iso4,jb) = ext_data%atm_td%aer_so4(jc,jb,imo1) + &
            ( ext_data%atm_td%aer_so4(jc,jb,imo2)  - ext_data%atm_td%aer_so4(jc,jb,imo1)  ) * wgt
          prm_diag%aerosol(jc,idu,jb) = ext_data%atm_td%aer_dust(jc,jb,imo1) + &
            ( ext_data%atm_td%aer_dust(jc,jb,imo2) - ext_data%atm_td%aer_dust(jc,jb,imo1) ) * wgt

        ENDDO
      ENDIF

      CALL ncn_from_tau_aerosol_speccnconst (nproma, nlev, i_startidx, i_endidx, nlev, nlev, &
        p_metrics%z_ifc(:,:,jb), prm_diag%aerosol(:,iss,jb), prm_diag%aerosol(:,iso4,jb),    &
        prm_diag%aerosol(:,iorg,jb), prm_diag%aerosol(:,idu,jb), zncn)

      CALL specccn_segalkhain_simple (nproma, i_startidx, i_endidx, zncn(:,nlev), prm_diag%cloud_num(:,jb))

      ! Impose lower limit on cloud_num over land
      DO jc = i_startidx, i_endidx
        IF (ext_data%atm%llsm_atm_c(jc,jb) .OR. ext_data%atm%llake_c(jc,jb)) &
          prm_diag%cloud_num(jc,jb) = MAX(175.e6_wp,prm_diag%cloud_num(jc,jb))
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE init_cloud_aero_cpl

  !------------------------------------------------
  ! Use climatological data of cloud droplet number 
  ! Satellite based data are provided in EXTPAR 
  !------------------------------------------------

  SUBROUTINE clim_cdnc(mtime_date, p_patch, ext_data, prm_diag)

    TYPE(datetime)       , INTENT(in)    :: mtime_date
    TYPE(t_patch)        , INTENT(in)    :: p_patch
    TYPE(t_external_data), INTENT(in)    :: ext_data

    TYPE(t_nwp_phy_diag) , INTENT(inout) :: prm_diag

    INTEGER  :: imo1, imo2
    INTEGER  :: rl_start, rl_end, i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: jb, jc

    REAL(wp) :: wgt

    TYPE(t_time_interpolation_weights) :: current_time_interpolation_weights

    TYPE(datetime), POINTER            :: mtime_hour

    CALL message('mo_nwp_phy_init:', 'Use climatological cdnc')

    mtime_hour => newDatetime(mtime_date)
    mtime_hour%time%minute = 0
    mtime_hour%time%second = 0
    mtime_hour%time%ms     = 0
    current_time_interpolation_weights = calculate_time_interpolation_weights(mtime_hour)
    call deallocateDatetime(mtime_hour)
    imo1 = current_time_interpolation_weights%month1
    imo2 = current_time_interpolation_weights%month2
    wgt  = current_time_interpolation_weights%weight2
    rl_start = 1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          ! Calculate the weighted average of monthly cloud droplet number
          prm_diag%cloud_num(jc,jb) = ( ext_data%atm_td%cdnc(jc,jb,imo1) + &
                   ( ext_data%atm_td%cdnc(jc,jb,imo2) - ext_data%atm_td%cdnc(jc,jb,imo1) ) * wgt )

          ! scaling of external cdnc with a scaling factor derived from the simple plumes
          IF ( atm_phy_nwp_config(p_patch%id)%lscale_cdnc ) THEN
              prm_diag%cloud_num(jc,jb) = prm_diag%cloud_num_fac(jc,jb) * prm_diag%cloud_num(jc,jb)
          ENDIF
        ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE clim_cdnc

END MODULE mo_nwp_phy_init

