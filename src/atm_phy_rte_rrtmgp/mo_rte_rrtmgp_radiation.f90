! Module to provide interface to radiation routines.
!
! Remarks
!   This module contains routines that provide the interface between ECHAM
!   and the radiation code.  Mostly it organizes and calculates the
!   information necessary to call the radiative transfer solvers of the
!   the rte-rrtmgp version (R. Pinus).
!
!         Sebastian Rast, MPI-M, Hamburg (2019-08-06): Renamings for
!              separation of this routine from former rte_rrtmgp code.
!
! Origin
!     This code is based on a former rte_rrtmgp implementation.
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

MODULE mo_rte_rrtmgp_radiation

  USE mo_kind,                ONLY: wp, i8
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_model_domain,        ONLY: t_patch

  USE mo_physical_constants,  ONLY: rae
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_aes_phy_config,      ONLY: aes_phy_tc
  USE mo_aes_rad_config,      ONLY: aes_rad_config
  USE mo_coupling_config,     ONLY: is_coupled_to_aero

  USE mo_orbit,               ONLY: orbit_kepler, orbit_vsop87, get_orbit_times
  USE mo_radiation_solar_data, ONLY: ssi_default, ssi_amip,                    &
                                     ssi_cmip5_picontrol, ssi_cmip6_picontrol, &
                                     ssi_RCEdiurnOn, ssi_RCEdiurnOff,          &
                                     ssi_radt, tsi_radt, ssi_RCEmip_analytical,&
                                     psctm,                                    &  
                                     ssi_factor
  USE mo_solar_parameters,    ONLY: solar_parameters
  USE mo_cloud_gas_profiles,  ONLY: gas_profiles, cloud_profiles, snow_profiles
  USE mo_radiation_general,   ONLY: nbndsw

  USE mo_rte_rrtmgp_interface,ONLY : rte_rrtmgp_interface

  USE mtime, ONLY: datetime, getTotalSecondsTimeDelta !, datetimeToString

  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: pre_rte_rrtmgp_radiation, rte_rrtmgp_radiation

  CONTAINS

  SUBROUTINE pre_rte_rrtmgp_radiation( p_patch,          datetime_radiation, &
                                     & current_datetime, ltrig_rad,          &
                                     & amu0_x,           rdayl_x,            &
                                     & amu0m_x,          rdaylm_x            )
  !-----------------------------------------------------------------------------
  !>
  ! Prepares information for radiation call

    TYPE(t_patch),           INTENT(in) :: p_patch
    TYPE(datetime), POINTER, INTENT(in) :: datetime_radiation, & !< date and time of radiative transfer calculation
         &                                 current_datetime       !< current time step
    LOGICAL,                 INTENT(in) :: ltrig_rad !< .true. if SW radiative transfer calculation has to be done at current time step
    REAL(wp),                INTENT(out) :: amu0_x(:,:), rdayl_x(:,:), &
         &                                  amu0m_x(:,:), rdaylm_x(:,:)

    LOGICAL  :: l_write_solar
    INTEGER(i8) :: icurrentyear
    INTEGER  :: icurrentmonth, i
    INTEGER, SAVE :: iprevmonth=-9999
    REAL(wp) :: rasc_sun, decl_sun, dist_sun
    REAL(wp) :: time_of_day, orbit_date
    REAL(wp) :: time_of_day_rt, orbit_date_rt
    REAL(wp) :: dt_ext
    REAL(wp) :: tsi

    ! Shortcuts to components of aes_rad_config
    !
    LOGICAL , POINTER :: l_orbvsop87, ldiur, l_sph_symm_irr, lyr_perp
    INTEGER , POINTER :: isolrad, icosmu0, yr_perp
    REAL(wp), POINTER :: fsolrad, cecc, cobld, clonp
    INTEGER  :: jg
!     CHARACTER(LEN=32)               :: datestring, datestring2
!     CHARACTER(LEN=*), PARAMETER     :: method_name="pre_rte_rrtmgp_radiation"

    !
    jg = p_patch%id
    isolrad        => aes_rad_config(jg)% isolrad
    fsolrad        => aes_rad_config(jg)% fsolrad
    l_orbvsop87    => aes_rad_config(jg)% l_orbvsop87
    cecc           => aes_rad_config(jg)% cecc
    cobld          => aes_rad_config(jg)% cobld
    clonp          => aes_rad_config(jg)% clonp
    lyr_perp       => aes_rad_config(jg)% lyr_perp
    yr_perp        => aes_rad_config(jg)% yr_perp
    ldiur          => aes_rad_config(jg)% ldiur
    l_sph_symm_irr => aes_rad_config(jg)% l_sph_symm_irr
    icosmu0        => aes_rad_config(jg)% icosmu0

    !
    ! 1.0 Compute orbital parameters for current time step
    ! --------------------------------
    !
    ! "time_of_day"  defines the local noon longitude for the time of this time step.
    ! "orbit_date" is not used. Instead always "orbit_date_rt" is used.
    CALL get_orbit_times(l_orbvsop87, current_datetime, lyr_perp, yr_perp, time_of_day, orbit_date)
    !
    ! "time_of_day_rt" defines the local noon longitude for the time of the
    ! radiative transfer calculation.
    ! "orbit_date_rt" defines the orbit position at the radiation time.
    ! This orbit position is kept constant through the radiation interval.
    CALL get_orbit_times(l_orbvsop87, datetime_radiation, lyr_perp, yr_perp, time_of_day_rt, orbit_date_rt)


    ! Compute the orbital parameters of Earth for "orbit_date_rt".
    IF (l_orbvsop87) THEN 
      CALL orbit_vsop87 (             orbit_date_rt, rasc_sun, decl_sun, dist_sun)
    ELSE
      CALL orbit_kepler (cecc, cobld, clonp, orbit_date_rt, rasc_sun, decl_sun, dist_sun)
    END IF

!     CALL datetimeToString(current_datetime,   datestring)
!     CALL datetimeToString(datetime_radiation, datestring2)
!     WRITE(message_text,'(4a)') ' current datetime: ', datestring, '; radiation datetime: ', datestring2
!     CALL message (method_name, message_text)


    dt_ext = 0.0_wp
    ! Compute cos(zenith angle) "amu0_x" for the current time "time_of_day", if the
    ! SW fluxes should be adjusted to the current sun (icosmu0=1:4), or the radiation
    ! time "time_of_day_rt", if the heating shall be computed for the sun position
    ! used for the radiative transfer (icosmu0=0).
    ! In both cases use the orbit parameters "decl_sun" and "dist_sun" valid for
    ! the "orbit_date_rt".
    !
    !DA TODO: move to the GPU?
    SELECT CASE (icosmu0)
    CASE (0)
       CALL solar_parameters(decl_sun,        time_of_day_rt,  &
            &                icosmu0,         dt_ext,          &
            &                ldiur,           l_sph_symm_irr,  &
            &                p_patch,                          &
            &                amu0_x,          rdayl_x          )
    CASE (1:4)
       CALL solar_parameters(decl_sun,        time_of_day,     &
            &                icosmu0,         dt_ext,          &
            &                ldiur,           l_sph_symm_irr,  &
            &                p_patch,                          &
            &                amu0_x,          rdayl_x          )
    CASE DEFAULT
       CALL finish('mo_rte_rrtmgp_radiation/pre_rte_rrtmgp_radiation','invalid icosmu0, must be in 0:4')
    END SELECT
    !
    !
    ! 2.0 Prepare time dependent quantities for rad (on radiation timestep)
    ! --------------------------------
    IF (ltrig_rad) THEN

      ! Compute cos(zenith angle) "amu0m_x" for the radiation time "time_of_day_rt"
      ! and the orbit parameters "decl_sun" and "dist_sun" valid for "orbit_date_rt".
      !
      ! "amu0m_x" is needed for the incoming SW flux field at the top of the atmosphere
      ! and for the optical paths.
      !
      SELECT CASE (icosmu0)
      CASE (0)
         dt_ext = 0.0_wp
      CASE (1:4)
         ! Extend the sunlit area for the radiative transfer calculations over an extended area
         ! including a rim of width dt_rad/2/86400*2pi (in radian) around the sunlit hemisphere.
         dt_ext = getTotalSecondsTimeDelta(aes_phy_tc(p_patch%id)%dt_rad,current_datetime)
      END SELECT
      !
      CALL solar_parameters(decl_sun,        time_of_day_rt,   &
           &                icosmu0,         dt_ext,           &
           &                ldiur,           l_sph_symm_irr,   &
           &                p_patch,                           &
           &                amu0m_x,         rdaylm_x          )
      !
      ! Consider curvature of the atmosphere for high zenith angles:
      ! The atmospheric path for a zenith angle mu0 through a spherical shell of
      ! thickness H and an inner radius R (and ratio rae=H/R) is shorter than
      ! a path at equal zenith angle through a plane parallel layer of thickness H.
      !
      WHERE (rdaylm_x == 1.0_wp) 
         amu0m_x(:,:)  = rae/(SQRT(amu0m_x(:,:)**2+rae*(rae+2.0_wp))-amu0m_x(:,:))
      END WHERE

      !++jsr&hs
      ! 3.0 Prepare possibly time dependent total solar and spectral irradiance
      ! --------------------------------
      ! ATTENTION: 
      ! This part requires some further work. Currently, a solar constant of
      ! 1361.371 is used as default. This is the TSI averaged over the
      ! years 1979 to 1988, and should be used for AMIP type runs. If lcouple is
      ! true, a solar constant of 1360.875 is used, the average for the years 1844
      ! to 1856. This should be used for a preindustrial control run.
      ! The spectral distribution of this TSI is currently also prescribed for
      ! these two cases depending on the lcouple switch.
      ! For transient CMIP5 simulations, the available time
      ! varying TSI and SSI has to be read in and used here.

      SELECT CASE (isolrad)
      CASE (0)
        tsi = SUM(ssi_default)
        ssi_factor = ssi_default
      CASE (1)
        tsi=tsi_radt
        ssi_factor=ssi_radt
      ! solar irradiance was read in aes_phy_bcs_global
      CASE (2)
        tsi = SUM(ssi_cmip5_picontrol)
        ssi_factor = ssi_cmip5_picontrol
      CASE (3)
        tsi = SUM(ssi_amip)
        ssi_factor = ssi_amip
      CASE (4)
        tsi = SUM(ssi_RCEdiurnOn)
        ssi_factor = ssi_RCEdiurnOn
      CASE (5)
        tsi = SUM(ssi_RCEdiurnOff)
        ssi_factor = ssi_RCEdiurnOff
      CASE (6)
        tsi = SUM(ssi_cmip6_picontrol)
        ssi_factor = ssi_cmip6_picontrol
      CASE (7)
        tsi = SUM(ssi_RCEmip_analytical)
        ssi_factor = ssi_RCEmip_analytical
      CASE default
        WRITE (message_text, '(a,i2,a)') &
             'isolrad = ', isolrad, ' in radctl namelist is not supported'
        CALL message('pre_radiation', message_text)
      END SELECT
      psctm(jg) = tsi/dist_sun**2 * fsolrad
      ssi_factor(:) = ssi_factor(:)/tsi

      ! output of solar constant every month

      icurrentmonth = datetime_radiation%date%month
      icurrentyear = datetime_radiation%date%year
      l_write_solar = icurrentmonth/=iprevmonth
      iprevmonth = icurrentmonth
      IF (l_write_solar) THEN
        CALL message('','')
        WRITE (message_text,'(a,i0,a,i2.2,a,f6.1)') &
             'Total solar constant [W/m^2] for ',      &
             icurrentyear, '-',                        &
             icurrentmonth, ' = ', tsi
        CALL message('',message_text)
        CALL message('','')
        DO i = 1, nbndsw
          WRITE (message_text,'(a,i2,a,f7.5)') &
               '   solar constant fraction: band ', i, &
               ' = ', ssi_factor(i)
          CALL message('',message_text)
        END DO
      END IF
    END IF ! ltrig_rad

  END SUBROUTINE pre_rte_rrtmgp_radiation
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  SUBROUTINE rte_rrtmgp_radiation( &
    jg, jb, jcs, jce, nproma, klev, ntracer, &
    & loland         ,&!< in  land-sea mask. (1. = land, 0. = sea/lakes)
    & loglac         ,&!< in  fraction of land covered by glaciers
    & this_datetime  ,&!< in  actual time step
    & pcos_mu0       ,&!< in  cosine of solar zenith angle
    & daylght_frc    ,&!< in  daylight fraction; with diurnal cycle 0 or 1, with zonal mean in [0,1]
    & alb_vis_dir    ,&!< in  surface albedo for visible range, direct
    & alb_nir_dir    ,&!< in  surface albedo for near IR range, direct
    & alb_vis_dif    ,&!< in  surface albedo for visible range, diffuse
    & alb_nir_dif    ,&!< in  surface albedo for near IR range, diffuse
    & emissivity     ,&!< in surface longwave emissivity
    & tk_sfc         ,&!< in  grid box mean surface temperature
    & zf             ,&!< in  geometric height at full level      [m]
    & zh             ,&!< in  geometric height at half level      [m]
    & dz             ,&!< in  geometric height thickness of layer [m]
    & pp_hl          ,&!< in  pressure at half levels at t-dt [Pa]
    & pp_fl          ,&!< in  pressure at full levels at t-dt [Pa]
    & rad_2d         ,&!< inout arbitrary 2d field in radiation for output
    & tk_fl          ,&!< in  tk_fl  = temperature at full level at t-dt
    & xm_air         ,&!< in  air mass in layer [kg/m2]
    & xq_trc         ,&!< in  tracer  mass fraction [kg/kg]
    & xv_ozn         ,&!< out ozone volume mixing ratio [mol/mol]
    !
    & reff_ice       ,&!< inout  effective radius of cloud ice (Fu needles) [m]
    & tau_ice        ,&!< inout  cloud ice optical depth, integrated over all bands
    & reff_snow      ,&!< inout  effective radius of snow (Fu needles) [m]
    & tau_snow       ,&!< inout  snow optical depth, integrated over all bands 
    & cdnc           ,&!< in  cloud droplet number concentration
    !
    & lw_dnw_clr     ,&!< out clear-sky downward longwave  at all levels
    & lw_upw_clr     ,&!< out clear-sky upward   longwave  at all levels
    & sw_dnw_clr     ,&!< out clear-sky downward shortwave at all levels
    & sw_upw_clr     ,&!< out clear-sky upward   shortwave at all levels
    & lw_dnw         ,&!< out all-sky   downward longwave  at all levels
    & lw_upw         ,&!< out all-sky   upward   longwave  at all levels
    & sw_dnw         ,&!< out all-sky   downward shortwave at all levels
    & sw_upw         ,&!< out all-sky   upward   shortwave at all levels
    !
    & vis_dn_dir_sfc ,&!< out all-sky downward direct visible radiation at surface
    & par_dn_dir_sfc ,&!< out all-sky downward direct PAR     radiation at surface
    & nir_dn_dir_sfc ,&!< out all-sky downward direct near-IR radiation at surface
    & vis_dn_dff_sfc ,&!< out all-sky downward diffuse visible radiation at surface
    & par_dn_dff_sfc ,&!< out all-sky downward diffuse PAR     radiation at surface
    & nir_dn_dff_sfc ,&!< out all-sky downward diffuse near-IR radiation at surface
    & vis_up_sfc     ,&!< out all-sky upward visible radiation at surface
    & par_up_sfc     ,&!< out all-sky upward PAR     radiation at surfac
    & nir_up_sfc     ,&!< out all-sky upward near-IR radiation at surface
    & aer_aod_533    ,&!< out  aerosol optical density at 533 nm
    & aer_ssa_533    ,&!< out  single scattering albedo at 533 nm
    & aer_asy_533    ,&!< out  asymmetrie factor at 533 nm
    & aer_aod_2325   ,&!< out  aerosol optical density at 2325 nm
    & aer_ssa_2325   ,&!< out  single scattering albedo at 2325 nm
    & aer_asy_2325   ,&!< out  asymmetrie factor at 2325 nm
    & aer_aod_9731    &!< out  aerosol optical density at 9731 nm
    &                 )

    INTEGER, INTENT(in)     :: &
    jg, jb, jcs, jce, nproma, klev, ntracer

    LOGICAL, INTENT(IN)     :: &
    & loland(:),           & !< land mask
    & loglac(:)              !< glacier mask

    TYPE(datetime), POINTER :: this_datetime !< actual time step

    REAL(wp), INTENT(IN)    :: &
    & pcos_mu0(:),         & !< cosine of solar zenith angle
    & daylght_frc(:),      & !< daylight fraction; with diurnal cycle 0 or 1, with zonal mean in [0,1]
    & alb_vis_dir(:),      & !< surface albedo for visible range and direct light
    & alb_nir_dir(:),      & !< surface albedo for NIR range and direct light
    & alb_vis_dif(:),      & !< surface albedo for visible range and diffuse light
    & alb_nir_dif(:),      & !< surface albedo for NIR range and diffuse light
    & emissivity(:),       & !< surface longwave emissivity
    & tk_sfc(:),           & !< Surface temperature
    & zf(:,:),          & !< geometric height at full level      [m]
    & zh(:,:),          & !< geometric height at half level      [m]
    & dz(:,:),          & !< geometric height thickness of layer [m]
    & pp_hl(:,:),       & !< pressure at half levels [Pa]
    & pp_fl(:,:),       & !< Pressure at full levels [Pa]
    & tk_fl(:,:),       & !< Temperature on full levels [K]
    & xm_air(:,:),      & !< air mass in layer [kg/m2]
    & xq_trc(:,:,:),    & !< tracer mass fraction [kg/kg]
    & cdnc(:,:)           !< Cloud drop number concentration
    REAL(wp), INTENT(INOUT) :: &
    & rad_2d(:),        & !< arbitrary 2d field in radiation for output
    & reff_ice(:,:),    & !< effective radius of cloud ice [m]
    & tau_ice(:,:),     & !< optical depth of cloud ice, integraded over all bands     
    & reff_snow(:,:),   & !< effective radius of snow [m]
    & tau_snow(:,:),    & !< optical depth snow, integrated over all bands     
    & xv_ozn(:,:)         !< ozone volume mixing ratio  [mol/mol]

    REAL(wp), TARGET, INTENT(INOUT)   :: &
    & lw_dnw_clr(:,:),& !< Clear-sky downward longwave  at all levels
    & lw_upw_clr(:,:),& !< Clear-sky upward   longwave  at all levels
    & sw_dnw_clr(:,:),& !< Clear-sky downward shortwave at all levels
    & sw_upw_clr(:,:),& !< Clear-sky upward   shortwave at all levels
    & lw_dnw(:,:),    & !< All-sky   downward longwave  at all levels
    & lw_upw(:,:),    & !< All-sky   upward   longwave  at all levels
    & sw_dnw(:,:),    & !< All-sky   downward shortwave at all levels
    & sw_upw(:,:)       !< All-sky   upward   shortwave at all levels

    REAL (wp), INTENT (INOUT) :: &
    & vis_dn_dir_sfc(:)    , & !< Diffuse downward flux surface visible radiation 
    & par_dn_dir_sfc(:)    , & !< Diffuse downward flux surface PAR
    & nir_dn_dir_sfc(:)    , & !< Diffuse downward flux surface near-infrared radiation
    & vis_dn_dff_sfc(:)    , & !< Direct  downward flux surface visible radiation 
    & par_dn_dff_sfc(:)    , & !< Direct  downward flux surface PAR
    & nir_dn_dff_sfc(:)    , & !< Direct  downward flux surface near-infrared radiation
    & vis_up_sfc    (:)    , & !< Upward  flux surface visible radiation 
    & par_up_sfc    (:)    , & !< Upward  flux surface PAR
    & nir_up_sfc    (:)    , & !< Upward  flux surface near-infrared radiation
    & aer_aod_533   (:,:)  , & !< aerosol optical density at 533 nm
    & aer_ssa_533   (:,:)  , & !< single scattering albedo at 533 nm
    & aer_asy_533   (:,:)  , & !< asymmetrie factor at 533 nm
    & aer_aod_2325  (:,:)  , & !< aerosol optical density at 2325 nm
    & aer_ssa_2325  (:,:)  , & !< single scattering albedo at 2325 nm
    & aer_asy_2325  (:,:)  , & !< asymmetrie factor at 2325 nm
    & aer_aod_9731  (:,:)      !< aerosol optical density at 9731 nm


    REAL (wp) ::      &
         !< general remark: volume mixing ratio always means
         !< (moles tracer)/(moles dry air)
    & xvmr_vap(nproma,klev),           & !< water vapor volume mixing ratio
    & xm_liq(nproma,klev),             & !< cloud water mass in layer [kg/m2]
    & xm_ice(nproma,klev),             & !< cloud ice   mass in layer [kg/m2]
    & xm_snw(nproma,klev),             & !< snow        mass in layer [kg/m2]
    & xvmr_co2(nproma,klev),           & !< CO2 volume mixing ratio
    & xvmr_o3(nproma,klev),            & !< O3  volume mixing ratio
    & xvmr_o2(nproma,klev),            & !< O2  volume mixing ratio
    & xvmr_ch4(nproma,klev),           & !< CH4 volume mixing ratio
    & xvmr_n2o(nproma,klev),           & !< N2O volume mixing ratio
    & xvmr_cfc(nproma,klev,2),         & !< CFC volume mixing ratio
    & xc_frc(nproma,klev)                !< cloud fraction

    REAL (wp) :: tk_hl(nproma,klev+1)

    REAL (wp) :: pp_sfc(nproma)

    INTEGER   :: jl, jk

    !$ACC DATA PRESENT(xv_ozn) &
    !$ACC   CREATE(pp_sfc, tk_hl, xm_liq, xm_ice, xc_frc, xm_snw) &
    !$ACC   CREATE(xvmr_vap, xvmr_co2, xvmr_o3, xvmr_o2, xvmr_ch4) &
    !$ACC   CREATE(xvmr_n2o, xvmr_cfc)

    CALL calculate_temperature_pressure(jcs, jce, nproma, klev, &
      pp_hl(:,:), pp_fl(:,:), tk_fl(:,:), tk_sfc(:), pp_sfc, tk_hl)

    CALL gas_profiles ( jg,           jb,                    jcs,             &
                      & jce,          nproma,                klev,            &
                      & ntracer,      this_datetime,         pp_hl,           &
                      & pp_fl,        xq_trc,                                 &
                      & xvmr_vap,     xvmr_co2,              xvmr_o3,         &
                      & xvmr_o2,      xvmr_ch4,              xvmr_n2o,        &
                      & xvmr_cfc                                              )
    !
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO jk=1,klev
      DO jl=jcs,jce
        xv_ozn(jl,jk) = xvmr_o3(jl,jk)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    CALL cloud_profiles ( jg,           jcs,            jce,              &
         &                klev,         xq_trc, xm_air,                   &
         &                xm_liq,       xm_ice,         xc_frc            )

    CALL snow_profiles  ( jg,           jcs,            jce,              &
         &                klev,         xq_trc, xm_air, xm_snw            )

    CALL rte_rrtmgp_interface(jg, jb, jcs, jce, nproma, klev             ,&
      aes_rad_config(jg)%irad_aero, is_coupled_to_aero()                 ,&
      psctm(jg), ssi_factor, loland, loglac, this_datetime               ,&
      pcos_mu0        ,daylght_frc                                       ,&
      alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif     ,&
      emissivity                                                         ,&
      zf              ,zh              ,dz                               ,&
      pp_sfc          ,pp_fl           ,pp_hl                            ,&
      tk_sfc          ,tk_fl           ,tk_hl                            ,&
      rad_2d                                                             ,&
      xvmr_vap        ,xm_liq          ,xm_ice                           ,&
      reff_ice        ,tau_ice         ,reff_snow       ,tau_snow        ,&
      cdnc            ,xc_frc          ,xm_snw                           ,&
      xvmr_co2        ,xvmr_ch4        ,xvmr_n2o        ,xvmr_cfc        ,&
      xvmr_o3         ,xvmr_o2                                           ,&
      lw_upw          ,lw_upw_clr      ,lw_dnw          ,lw_dnw_clr      ,&
      sw_upw          ,sw_upw_clr      ,sw_dnw          ,sw_dnw_clr      ,&
      vis_dn_dir_sfc  ,par_dn_dir_sfc  ,nir_dn_dir_sfc                   ,&
      vis_dn_dff_sfc  ,par_dn_dff_sfc  ,nir_dn_dff_sfc                   ,&
      vis_up_sfc      ,par_up_sfc      ,nir_up_sfc                       ,&
      aer_aod_533     ,aer_ssa_533     ,aer_asy_533                      ,&
      aer_aod_2325    ,aer_ssa_2325    ,aer_asy_2325                     ,&
      aer_aod_9731                                                       ) 
    !$ACC WAIT
    !$ACC END DATA
    !-------------------------------------------------------------------

  END SUBROUTINE rte_rrtmgp_radiation
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE calculate_temperature_pressure (               &
      & jcs           ,&!< in  start index for loop over block
      & jce           ,&!< in  end   index for loop over block
      & kbdim         ,&!< in  dimension of block over cells
      & klev          ,&!< in  number of full levels = number of layers
      & pp_hl         ,&! in
      & pp_fl         ,&
      & tk_fl         ,&
      & tk_sfc        ,&
      & pp_sfc        ,& ! out
      & tk_hl   )

  
    INTEGER, INTENT(in) :: &
      & jcs               ,&
      & jce               ,&
      & kbdim             ,&
      & klev              
    ! in 
    REAL(wp), INTENT(in) :: &
      & pp_hl(kbdim,klev+1),&
      & pp_fl(kbdim,klev)  ,&
      & tk_fl(kbdim,klev)  ,&
      & tk_sfc(kbdim)

    !out      
    REAL(wp), INTENT(inout) :: &
      & pp_sfc (kbdim)        ,&
      & tk_hl  (kbdim,klev+1)  
    
    ! local counters
    INTEGER              :: jk, jl
    !
    !$ACC DATA PRESENT(pp_hl, tk_fl, pp_fl, tk_sfc, pp_sfc, tk_hl)
    !
    ! 1.0 calculate variable input parameters (location and state variables)
    ! --------------------------------
    ! 
    ! --- Pressure (surface and distance between half levels)
    !
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
    pp_sfc(jcs:jce)   = pp_hl(jcs:jce,klev+1)
    !$ACC END KERNELS
    !
    ! --- temperature at half levels
    !
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO jk=2,klev
      DO jl = jcs, jce
        tk_hl(jl,jk) = (tk_fl(jl,jk-1)*pp_fl(jl,jk-1)*( pp_fl(jl,jk)          &
             & - pp_hl(jl,jk) ) + tk_fl(jl,jk)*pp_fl(jl,jk)*( pp_hl(jl,jk)    &
             & - pp_fl(jl,jk-1))) /(pp_hl(jl,jk)*(pp_fl(jl,jk) -pp_fl(jl,jk-1)))
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
    tk_hl(jcs:jce,klev+1) = tk_sfc(jcs:jce)
    !$ACC END KERNELS

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO jl = jcs, jce
      tk_hl(jl,1)      = tk_fl(jl,1)-pp_fl(jl,1)*(tk_fl(jl,1) - tk_hl(jl,2))  &
            &             / (pp_fl(jl,1)-pp_hl(jl,2))
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC END DATA
  END SUBROUTINE calculate_temperature_pressure 
  !-------------------------------------------------------------------

  END MODULE mo_rte_rrtmgp_radiation
