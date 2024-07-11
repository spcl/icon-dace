!
! mo_art_emission_biomBurn
! This module provides routines for preparing biomass burning paramtererization by using landuse
! data and creating diurnal cycle.
! All changes made after initial creation are marked in comments. Main parts equals to the
! original version of S. Freitas.
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

MODULE mo_art_emission_biomBurn
! ICON
  USE mo_kind,                          ONLY: wp
  USE mtime,                            ONLY: datetime
  USE mo_impl_constants,                ONLY: MAX_CHAR_LENGTH
  USE mo_math_constants,                ONLY: pi
! ART
  USE mo_art_external_types,            ONLY: t_art_biomBurn_properties
  USE mo_art_emission_plumerise,        ONLY: smk_pr_driver
  USE mo_art_emission_BBPlume,          ONLY: BB_plume, a_priori_BB_plume
  USE mo_art_impl_constants,            ONLY: UNDEF_REAL_ART


  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=MAX_CHAR_LENGTH) :: thisroutine = 'mo_art_emission_biomBurn'

  PUBLIC :: art_emission_biomBurn_prepare, art_emission_biomBurn, art_emission_biomBurn_bbplume

CONTAINS

SUBROUTINE art_emission_biomBurn_prepare(lu_frac_crop_irrig,lu_frac_crop_rain,lu_frac_crop_mos,  &
  &                                      lu_frac_veg_mos,lu_frac_forest_b_eg,lu_frac_forest_b_d, &
  &                                      lu_frac_woodland,lu_forest_n_eg,lu_frac_forest_n_d,     &
  &                                      lu_frac_forest_bn,lu_frac_shrub_mos,lu_frac_shrub_eg,   &
  &                                      lu_frac_shrub,lu_frac_grass,lu_frac_sparse,             &
  &                                      lu_frac_forest_rf,lu_frac_forest_pf,lu_frac_grass_rf,   &
  &                                      dc_hflux_min_res,dc_hflux_max_res,dc_burnt_area_res,    &
  &                                      dc_emis_res,                                            &
  &                                      jb, istart, iend)
!>
!! SUBROUTINE art_emission_biomBurn_prepare
!! - This routine calculates a diurnal cycle of possible heatfluxes, emissions and burned areas
!!   due to biomass burning activities based on landuse data. It's used as a prepatation for
!!   each grid point to make it capable for biomass burning.
!! -
!! Based on Walter et al. - The importance of plume rise on the concentrations and
!!                          atmospheric impacts of biomass burning aerosol
!! Part of Module: mo_art_emssion_biomBurn
!! Author: Jonas Straub, KIT
!! Initial Release: 2018-02-19
!!
!! Modifications:
!! YYYY-MM-DD: <name>, <institution>
!! - ...
!<

  INTEGER,INTENT(IN)        :: &
    &  jb,                     & !< Block index
    &  istart, iend              !< Start and end index of nproma loop

  REAL(wp),INTENT(IN)       :: &
    &  lu_frac_crop_irrig(:),  & !< Land-cover fraction per classification,
    &  lu_frac_crop_rain(:),   & !  exact meaning of each classification is defined in
    &  lu_frac_crop_mos(:),    & !  mo_ext_data_types, respectively database: GLOBCOVER2009
    &  lu_frac_veg_mos(:),     & !
    &  lu_frac_forest_b_eg(:), & !  range: 0 <= lu_frac_ <= 1
    &  lu_frac_forest_b_d(:),  & !    0: no plant of this classification is present
    &  lu_frac_woodland(:),    & !    1: whole grid is covered by plants of this classification
    &  lu_forest_n_eg(:),      & !
    &  lu_frac_forest_n_d(:),  & !
    &  lu_frac_forest_bn(:),   & !
    &  lu_frac_shrub_mos(:),   & !
    &  lu_frac_shrub_eg(:),    & !
    &  lu_frac_shrub(:),       & !
    &  lu_frac_grass(:),       & !
    &  lu_frac_sparse(:),      & !
    &  lu_frac_forest_rf(:),   & !
    &  lu_frac_forest_pf(:),   & !
    &  lu_frac_grass_rf(:)       !

  REAL(wp),INTENT(INOUT)  :: &
    &  dc_hflux_min_res(:,:,:),     & !< diurnal cycle for heat flux min W m^-2
    &  dc_hflux_max_res(:,:,:),     & !< diurnal cycle for heat flux max W m^-2
    &  dc_burnt_area_res(:,:,:),    & !< diurnal cycle for burnt area m^2
    &  dc_emis_res(:,:,:)             !< (basic) diurnal cycle for emissions
                                    !  (multiplied later with GFAS FRP emissions)

  ! local
  INTEGER                  :: &
    &  hour,                  & !< loop variable
    &  jc                       !< loop variable for nproma

  INTEGER,PARAMETER        :: &
    &  bio_types=3              !< amount of different biome types (=3 according Freitas, 2006)

  REAL(wp)                 :: &
    &  hflux_min_res,         & !< minimum heatflux according to vegetation type(s) kW m^-2
    &  hflux_max_res,         & !< maximum heatflux according to vegetation type(s) kW m^-2
    &  dc_wgh_hflux,          & !< diurnal weighted heatflux
    &  dc_wgh_emis,           & !< diurnal weighted emission
    &  sum_dc_t_emis,         & !< factor to normalize diurnal cycle tropical forest
    &  sum_dc_s_emis,         & !< factor to normalize diurnal cycle savannas
    &  sum_dc_g_emis            !< factor to normalize diurnal cycl grassland+cropland


  REAL(wp),POINTER         :: &
    &  dc_basic(:),           & !< diurnal hourly basic profile
    &  dc_t_emis(:),          & !< diurnal hourly profile of burning tropical forest
    &  dc_s_emis(:),          & !< diurnal hourly profile of burning savannas
    &  dc_g_emis(:),          & !< diurnal hourly profile of burning grassland+cropland
    &  dc_t_hflux(:),         & !< heatflux and burnt area, tropical forest
    &  dc_s_hflux(:),         & !< heatflux and burnt area, savanna
    &  dc_g_hflux(:),         & !< heatflux and burnt area, grassland/cropland
    &  frp_veg(:),            & !< biome type cluster 1-3
    &  normveg(:)               !< normalize merged vegetation class fractions, all 3 biome types

  REAL(wp),PARAMETER       :: &
    &  hmin_class_1 = 30._wp, & !< heat flux limits of Freitas et al. (2006), table 1
    &  hmax_class_1 = 80._wp, & !  in kW /m2
    &  hmin_class_2 = 4.4_wp, & !
    &  hmax_class_2 = 23._wp, & !
    &  hmin_class_3 = 3.3_wp, & !
    &  hmax_class_3 = 3.3_wp    !
  ! ----------------------------------------------------------------

  ! -------------------------------------------
  ! --- 0. Allocate & initialize local storage
  ! -------------------------------------------
  ALLOCATE( dc_basic  (24) )
  ALLOCATE( dc_t_emis (24) )
  ALLOCATE( dc_s_emis (24) )
  ALLOCATE( dc_g_emis (24) )
  ALLOCATE( dc_t_hflux(24) )
  ALLOCATE( dc_s_hflux(24) )
  ALLOCATE( dc_g_hflux(24) )
  ALLOCATE( frp_veg(bio_types) )
  ALLOCATE( normveg(bio_types) )

  dc_basic          = 0.0_wp
  dc_t_emis         = 0.0_wp
  dc_s_emis         = 0.0_wp
  dc_g_emis         = 0.0_wp
  dc_t_hflux        = 0.0_wp
  dc_s_hflux        = 0.0_wp
  dc_g_hflux        = 0.0_wp
  frp_veg           = 0.0_wp
  normveg           = 0.0_wp


  ! ----------------------------------
  ! --- 1. Create basic profile
  ! ----------------------------------
  ! diurnal cycle basic profile (gaussian distribution with sigma=2.5 and mu=12.5) - hour=1-24
  ! means 00-23 local solar time (LST)
  DO hour=1,24
    dc_basic(hour) = (1._wp/(2.5_wp*SQRT(2._wp*pi)))     &
      &              * EXP( (1._wp-1.5_wp)                &
      &                     * ( ((REAL(hour,wp)-1._wp)-12.5_wp)**2 &
      &                         / 2.5_wp**2 ) )

  ! Transform basic profile according to Zhang and Kondragunta (2008):
  ! - peak fire activities (PFA): 10-15 LST
  ! - (A) forest  : 52.1% of total daily burnt areas in PFA
  ! - (B) savannas: 61.8% of total daily burnt areas in PFA
  ! - (C) (grassland+cropland)/2: 73.7% of total daily burnt areas in PFA
  ! The three biome types are selected according to the available heat flux
  ! values of Freitas et al. (2006), Table 1.

  ! DIURNAL CYCLE BASIC PROFILE
  ! Tropical forest: add 0.0378  to fullfill (A)
    dc_t_emis(hour) = (0.0387_wp  + dc_basic(hour))
  ! Savannas       : add 0.01754 to fullfill (B)
    dc_s_emis(hour) = (0.01754_wp + dc_basic(hour))
  ! Grassland      : add 0.00308 to fullfill (C)
    dc_g_emis(hour) = (0.00308_wp + dc_basic(hour))
  ENDDO

  ! Transform basic profile to basic profile for emission (_emis)
  ! and heat flux and burnt area (_hflux)
  sum_dc_t_emis = sum(dc_t_emis)
  sum_dc_s_emis = sum(dc_s_emis)
  sum_dc_g_emis = sum(dc_g_emis)
  DO hour=1,24
  ! DIURNAL CYCLE BASIC PROFILE FOR EMISSIONS(+)
  ! Tropical forest: normalize to 1
    dc_t_emis(hour) = dc_t_emis(hour)/sum_dc_t_emis
  ! Savannas       : normalize to 1
    dc_s_emis(hour) = dc_s_emis(hour)/sum_dc_s_emis
  ! Grassland      : normalize to 1
    dc_g_emis(hour) = dc_g_emis(hour)/sum_dc_g_emis

  ! DIURNAL CYCLE FOR BASIC PROFILE HEAT FLUX AND BURNT AREA(*)
  ! Tropical forest: scale with 9.89  to set max to 1
    dc_t_hflux(hour) = 9.89_wp  * dc_t_emis(hour)
  ! Savannas       : scale with 8.17  to set max to 1
    dc_s_hflux(hour) = 8.17_wp  * dc_s_emis(hour)
  ! Grassland      : scale with 6.735 to set max to 1
    dc_g_hflux(hour) = 6.735_wp * dc_g_emis(hour)

  !(+) Later in art_emission_biomBurn: this basic profile is multiplied with the daily sum of
  !                                    the emissions.
  !(*) Later here: this basic profile is multiplied with the burnt area and the resulting
  !                heatflux min and max according to the weighting sum of the three biome types.
  !                The max is set to 1 to realize that the diurnal maximum is equal to hflux_max
  !                or hflux_min.
  ENDDO !hour

  ! ----------------------------------
  ! --- 2. Merge landuse classes to the three biome types of Freitas et al. (2006)
  ! ----------------------------------
  !  Landuse class attribution is made for GLOBCOVER2009 datset!

  ! ----------------------------------
  ! --- Loop over all horizontal gridpoints
  ! ----------------------------------
  DO jc = istart, iend

    ! Tropical forest (all tree covers taken as tropical forest)
    frp_veg(1) = lu_frac_forest_b_eg(jc) + lu_frac_forest_b_d(jc) + lu_frac_woodland(jc)   &
      &        + lu_forest_n_eg(jc)      + lu_frac_forest_n_d(jc) + lu_frac_forest_bn(jc)  &
      &        + lu_frac_forest_rf(jc)   + lu_frac_forest_pf(jc)
    ! Woody savanna - cerrado
    frp_veg(2) = lu_frac_shrub(jc)
    ! Grassland - pasture - cropland
    frp_veg(3) = lu_frac_crop_irrig(jc)  + lu_frac_crop_rain(jc)  + lu_frac_crop_mos(jc)   &
      &        + lu_frac_veg_mos(jc)     + lu_frac_shrub_mos(jc)  + lu_frac_shrub_eg(jc)   &
      &        + lu_frac_grass(jc)       + lu_frac_sparse(jc)     + lu_frac_grass_rf(jc)

    ! go on if nonzero vegetation fraction available (otherwise the GFAS fire pixel is not
    !                                                 accepted)
    IF(frp_veg(1)>0._wp .or. frp_veg(2)>0._wp .or. frp_veg(3)>0._wp) THEN
      ! normalize merged vegetation class fractions to 1 and ...
      normveg(1)  = frp_veg(1) / (frp_veg(1) + frp_veg(2) + frp_veg(3))
      normveg(2)  = frp_veg(2) / (frp_veg(1) + frp_veg(2) + frp_veg(3))
      normveg(3)  = frp_veg(3) / (frp_veg(1) + frp_veg(2) + frp_veg(3))

      ! ...calculate resulting heat flux limits (weighting sum) kW m^-2
      hflux_min_res = hmin_class_1*normveg(1) &  !forest
        &           + hmin_class_2*normveg(2) &  !savanna
        &           + hmin_class_3*normveg(3)    !grassland
      hflux_max_res = hmax_class_1*normveg(1) &  !forest
        &           + hmax_class_2*normveg(2) &  !savanna
        &           + hmax_class_3*normveg(3)    !grassland

      ! calculate resulting diurnal cycle
      DO hour=1,24
        dc_wgh_hflux = dc_t_hflux(hour)*normveg(1) &  !forest
          &          + dc_s_hflux(hour)*normveg(2) &  !savanna
          &          + dc_g_hflux(hour)*normveg(3)    !grassland
        dc_wgh_emis  = dc_t_emis(hour) *normveg(1) &  !forest
          &          + dc_s_emis(hour) *normveg(2) &  !savanna
          &          + dc_g_emis(hour) *normveg(3)    !grassland

        ! the diurnal cycles for hflux min/max and burnt areas can be used directly in
        ! art_emission_biomBurn
        ! change from kW m^-2 to SI-Unit W m^-2 for the heatfluxes
        dc_hflux_min_res(jc,jb,hour)  = dc_wgh_hflux * hflux_min_res * 1000._wp
        dc_hflux_max_res(jc,jb,hour)  = dc_wgh_hflux * hflux_max_res * 1000._wp 
        dc_burnt_area_res(jc,jb,hour) = dc_wgh_hflux * 500000._wp ! use 50 ha (500000 m2)
        ! the diurnal cycles for emis have to be multiplied with the daily sum
        ! of emissions (kg m-2 d-1)
        dc_emis_res(jc,jb,hour)       = dc_wgh_emis
      ENDDO ! hour


    ELSE
      ! vegetation fraction is zero - no emission possible
      dc_hflux_min_res(jc,jb,:)  = 0._wp
      dc_hflux_max_res(jc,jb,:)  = 0._wp
      dc_burnt_area_res(jc,jb,:) = 0._wp
      dc_emis_res(jc,jb,:)       = 0._wp
    ENDIF ! nonzero vegetation fraction

  END DO !jc

  DEALLOCATE( dc_basic ,dc_t_emis,dc_s_emis,dc_g_emis,  &
    &         dc_t_hflux,dc_s_hflux,dc_g_hflux,frp_veg,normveg )


END SUBROUTINE art_emission_biomBurn_prepare
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_emission_biomBurn(t,p,u,v,qv,z_mc,z_ifc,lon,area,dz,current_date,  &
  &                              dc_hflux_min_res,dc_hflux_max_res,dc_burnt_area_res, &
  &                              dc_emis_res,flux_bc,emiss_rate,nlev,nlevp1,jb,istart,iend )

!>
!! SUBROUTINE art_emission_biomBurn
!! - Maintain the smoke plume height calculation and initialize the tracer release
!!   into the atmosphere.
!! -
!! Based on Walter et al. - The importance of plume rise on the concentrations and
!!                          atmospheric impacts of biomass burning aerosol
!! Part of Module: mo_art_emission_biomBurn
!! Author: Jonas Straub, KIT
!! Initial Release: 2018-02-19
!!
!! Modifications:
!! YYYY-MM-DD: <name>, <institution>
!! - ...
!<

  TYPE(datetime),INTENT(IN) :: &
  &  current_date                !< Date and time information

  INTEGER,INTENT(IN)        :: &
    &  jb,                     & !< Block index
    &  istart, iend,           & !< Start and end index of nproma loop
    &  nlev,                   & !< Number of full levels (equals index of lowest full level)
    &  nlevp1                    !< Number of half levels (equals index of lowest half level)

  REAL(wp),INTENT(INOUT)    :: &
    &  emiss_rate(:,:)           !< Emission rates [mug m-3 s-1], dimensions: nproma, nlev

  REAL(wp),INTENT(IN)       :: &
    ! dimensions: (nproma)
    &  lon(:),                 & !< longitude (rad)
    &  area(:),                & !< area of triangle (m2)
    ! dimensions: (nproma, nlev)
    &  t(:,:),                 & !< air temperature  (K)
    &  p(:,:),                 & !< pressure (Pa)
    &  u(:,:),                 & !< zonal windspeed (m/s)
    &  v(:,:),                 & !< meridional windspeed (m/s)
    &  z_mc(:,:),              & !< geometric height at full levels (m)
    &  z_ifc(:,:),             & !< geometric height at half levels (m)
    &  qv(:,:),                & !< specific humidity (on model levels) (kg/kg)
    &  dz(:,:),                & !< Height of model layer (m)
    ! dimensions: (nproma, hour)
    &  dc_hflux_min_res(:,:,:),  & !< diurnal cycle for heat flux min (W m^-2)
    &  dc_hflux_max_res(:,:,:),  & !< diurnal cycle for heat flux max (W m^-2)
    &  dc_burnt_area_res(:,:,:), & !< diurnal cycle for burnt area (m^2)
    &  dc_emis_res(:,:,:),       & !< (basic) diurnal cycle for emissions
                                 !  (multiplied later with GFAS FRP emissions)
    &  flux_bc(:)                !< wildfires flux of black carbon, proportional to FRP (kg/(m2 s)

  !local
  INTEGER                   :: &
    &  kfrp,                   & !< loop variable: vertical layer index for plumerise variables
    &  jk,                     & !< loop variables: vertical layer index for ICON variables,
    &  jc                        !< loop variables: nproma index

  REAL(wp)                      :: &
    ! unit is not ICON-conform, but is expected in plumerise model!
    &  ztop_pr,                & !< top height for plume, return-value of plumerise model (m)
    &  zbot_pr,                & !< bottom height for plume, return-value of plumerise model(m)
    &  topo,                   & !< lowest halflevel (m)
    &  heat_flux_top,          & !< local heatflux for smoke plume top height (W m^-2)
    &  heat_flux_bot,          & !< local heatflux for smoke plume bottom height (W m^-2)
    &  burnt_area                !< local burnt area at current hour, based on 50ha (m^2)

  REAL(wp)                  :: &
    &  ztop,                   & !< postprocessing corrected top height (m)
    &  zbot,                   & !< postprocessing correcteed bottom height (m)
    &  timediff_UTC,           & !< local timeshift to UTC
    &  ls_time                   !< local solar time

  REAL(wp),POINTER              :: &
    ! unit is not ICON-conform, but is expected in plumerise model!
    &  tenv(:),                & !< environmental atmospheric state for plumerise model.
    &  penv(:),                & !  layer description is vice versa than in ICON
    &  uenv(:),                & !  naming: <variable>env
    &  venv(:),                & !
    &  hmlev(:),               & !  full levels height; corresponds to z_mc
    &  qvenv(:)                  !
  ! ----------------------------------------------------------------

  ! ----------------------------------
  ! --- Loop over all horizontal gridpoints
  ! ----------------------------------
  DO jc = istart, iend

    ! Check if there is a detected fire at this gridpoint
    IF ( flux_bc(jc) > 0.005e-9_wp .and. dc_emis_res(jc,jb,13) > 0._wp ) THEN

      ! -------------------------------------------
      ! --- 0. Allocate & initialize local storage
      ! -------------------------------------------
      ALLOCATE( tenv (nlev) )
      ALLOCATE( penv (nlev) )
      ALLOCATE( uenv (nlev) )
      ALLOCATE( venv (nlev) )
      ALLOCATE( hmlev(nlev) )
      ALLOCATE( qvenv(nlev) )

      tenv  = UNDEF_REAL_ART
      penv  = UNDEF_REAL_ART
      uenv  = UNDEF_REAL_ART
      venv  = UNDEF_REAL_ART
      hmlev = UNDEF_REAL_ART
      qvenv = UNDEF_REAL_ART

      ! -------------------------------------------
      ! --- 1. Generate vertical inverted environment variable field for plumerise model
      ! -------------------------------------------
      DO kfrp=1,nlev
        ! transform plumerise-index to ICON-index
        jk = nlev+1-kfrp
        tenv (kfrp) = t(jc,jk)
        penv (kfrp) = p(jc,jk)
        uenv (kfrp) = u(jc,jk)
        venv (kfrp) = v(jc,jk)
        hmlev(kfrp) = z_mc(jc,jk)
        qvenv(kfrp) = qv(jc,jk)
      ENDDO !jk

      ! height of lowest halflevel
      topo = z_ifc(jc,nlevp1)
      ! local timeshift to UTC, in hours
      timediff_UTC = lon(jc) / (2._wp*pi) * 24._wp
      ! local solar time, sun is at the zenith at ls_time=12.0
      ls_time = MODULO(current_date%time%hour + timediff_UTC + 24._wp, 24._wp)

      ! derive local biomass burning properties for the plumerise model
      heat_flux_top = dc_hflux_max_res(jc,jb,CEILING(ls_time))
      burnt_area    = dc_burnt_area_res(jc,jb,CEILING(ls_time))
      heat_flux_bot = dc_hflux_min_res(jc,jb,CEILING(ls_time))


      ! -------------------------------------------
      ! --- 2. calculate top and bottom height of smoke plume
      ! -------------------------------------------
      CALL smk_pr_driver(nlev,heat_flux_top,heat_flux_bot,burnt_area,tenv,penv,uenv,venv, &
        &                hmlev,qvenv,topo,ztop_pr,zbot_pr)

      ! make sure vertical plume expansion is > 0m, default: 100m
      ztop = 1.0_wp*ztop_pr ! transform accuracy to (wp)
      zbot = 1.0_wp*zbot_pr

      IF (ztop == 0._wp) ztop = ztop + 100._wp
      IF (zbot == ztop)  zbot = ztop - 100._wp

      ! skip heights on half levels or something like that?!
      ! maybe for later calculatations in calc_wgt()
      ztop = ztop+topo
      zbot = zbot+topo

      ! -------------------------------------------
      ! --- 3. Calculate and set vertical distribution of emitted fire aerosol
      ! -------------------------------------------

      CALL art_calc_wgt(ztop,zbot,nlev,z_ifc(jc,:),emiss_rate(jc,:),                     &
        &               dc_emis_res(jc,jb,CEILING(ls_time)),flux_bc(jc),area(jc),dz(jc,:) )

      DEALLOCATE(tenv,penv,uenv,venv,hmlev,qvenv)

    ENDIF ! flux_bc


  ENDDO !jc


END SUBROUTINE art_emission_biomBurn
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_emission_biomBurn_bbplume(t,p,u,v,qv,z_mc,z_ifc,lon,area,dz,current_date,  &
  &                              dc_hflux_min_res,dc_hflux_max_res,dc_burnt_area_res, &
  &                              dc_emis_res,flux_bc,emiss_rate,nlev,nlevp1,jb,istart,iend )

!>
!! SUBROUTINE art_emission_biomBurn_bbplume
!! - Maintain the smoke plume height calculation and initialize the tracer release
!!   into the atmosphere.
!! - features updated implementation of the plume rise model
!! Based on Walter et al. - The importance of plume rise on the concentrations and
!!                          atmospheric impacts of biomass burning aerosol
!! Part of Module: mo_art_emission_biomBurn
!<

  TYPE(datetime),INTENT(IN) :: &
  &  current_date                !< Date and time information

  INTEGER,INTENT(IN)        :: &
    &  jb,                     & !< Block index
    &  istart, iend,           & !< Start and end index of nproma loop
    &  nlev,                   & !< Number of full levels (equals index of lowest full level)
    &  nlevp1                    !< Number of half levels (equals index of lowest half level)

  REAL(wp),INTENT(INOUT)    :: &
    &  emiss_rate(:,:)           !< Emission rates [mug m-3 s-1], dimensions: nproma, nlev

  REAL(wp),INTENT(IN)       :: &
    ! dimensions: (nproma)
    &  lon(:),                 & !< longitude (rad)
    &  area(:),                & !< area of triangle (m2)
    ! dimensions: (nproma, nlev)
    &  t(:,:),                 & !< air temperature  (K)
    &  p(:,:),                 & !< pressure (Pa)
    &  u(:,:),                 & !< zonal windspeed (m/s)
    &  v(:,:),                 & !< meridional windspeed (m/s)
    &  z_mc(:,:),              & !< geometric height at full levels (m)
    &  z_ifc(:,:),             & !< geometric height at half levels (m)
    &  qv(:,:),                & !< specific humidity (on model levels) (kg/kg)
    &  dz(:,:),                & !< Height of model layer (m)
    ! dimensions: (nproma, hour)
    &  dc_hflux_min_res(:,:,:),  & !< diurnal cycle for heat flux min (W m^-2)
    &  dc_hflux_max_res(:,:,:),  & !< diurnal cycle for heat flux max (W m^-2)
    &  dc_burnt_area_res(:,:,:), & !< diurnal cycle for burnt area (m^2)
    &  dc_emis_res(:,:,:),       & !< (basic) diurnal cycle for emissions
                                 !  (multiplied later with GFAS FRP emissions)
    &  flux_bc(:)                !< wildfires flux of black carbon, proportional to FRP (kg/(m2 s)

  !local
  INTEGER                   :: &
    &  step = 5,               & !<  uneven integer - thinning icon data passed to plume rise model
    &  toplvl,                 & !< index of highst level in thinned data (from surface to space)
    &  k,                      & !< loop index
    &  kfrp,                   & !< loop variable: vertical layer index for plumerise variables
    &  jk,                     & !< loop variables: vertical layer index for ICON variables,
    &  jc                        !< loop variables: nproma index

  INTEGER, allocatable      :: &
    &  mask1(:),               & ! mask to thin icon data passed toi plume rise model
    &  mask2(:)                  ! mask to thin on half levels

  REAL(wp)                      :: &
    ! unit is not ICON-conform, but is expected in plumerise model!
    &  ztop_pr,                & !< top height for plume, return-value of plumerise model (m)
    &  zbot_pr,                & !< bottom height for plume, return-value of plumerise model (m)
    &  topo,                   & !< lowest halflevel (m)
    &  heat_flux_top,          & !< local heatflux for smoke plume top height (W m^-2)
    &  heat_flux_bot,          & !< local heatflux for smoke plume bottom height (W m^-2)
    &  burnt_area                !< local burnt area at current hour, based on 50ha (m^2)

  REAL(wp)                  :: &
    &  ztop,                   & !< postprocessing corrected top height (m)
    &  zbot,                   & !< postprocessing correcteed bottom height (m)
    &  timediff_UTC,           & !< local timeshift to UTC
    &  ls_time                   !< local solar time

  LOGICAL                   :: & ! logical to check with a priori condition:
    & l_cpr_top,               & ! - whether a plume top computation is necessary
    & l_cpr_bot                  ! - whether a plume bot computation is necessary
           !
  ! ----------------------------------------------------------------

  ! ----------------------------------
  ! --- Loop over all horizontal gridpoints
  ! ----------------------------------
  DO jc = istart, iend

    ! Check if there is a detected fire at this gridpoint
    IF ( flux_bc(jc) > 5.0e-14_wp .and. dc_emis_res(jc,jb,13) > 0._wp ) THEN

      ! height of lowest halflevel
      topo = z_ifc(jc,nlevp1)
      ! local timeshift to UTC, in hours
      timediff_UTC = lon(jc) / (2._wp*pi) * 24._wp
      ! local solar time, sun is at the zenith at ls_time=12.0 - make sure that ls_time never equals zero
      ls_time = MODULO(current_date%time%hour + timediff_UTC + 24._wp, 24._wp) + 1.e-12_wp

      ! derive local biomass burning properties for the plumerise model
      heat_flux_top = dc_hflux_max_res(jc,jb,CEILING(ls_time)) !(W m^-2)
      burnt_area    = dc_burnt_area_res(jc,jb,CEILING(ls_time))
      heat_flux_bot = dc_hflux_min_res(jc,jb,CEILING(ls_time)) !(W m^-2)

      ! thinning icon vertical data to cope with large vertical velocities in BB_plume
      ! also invert order from 'space to surface' to 'surface to space'
      toplvl=(nlevp1-1)/step ! divide available levels to sparce data
      allocate(mask1(toplvl),mask2(toplvl+1)) ! masks to thin data
      mask1 = (/(k,k=CEILING(1.5_wp*step)+(toplvl-2)*step,CEILING(1.5_wp*step),-step),1/)
      mask2 = (/(k,k=1+toplvl*step,1,-step)/) ! keep lowest level and skip with step through vertical

      l_cpr_top = .FALSE.
      l_cpr_bot = .FALSE.

      CALL a_priori_BB_plume(heat_f = heat_flux_top, & ! (W m^-2)
           &                   area =    burnt_area, & ! (m^2)
           &                     pe =    p(jc,nlev), & ! surface pressure
           &                  l_cpr =     l_cpr_top)
      CALL a_priori_BB_plume(heat_f = heat_flux_bot, & ! (W m^-2)
           &                   area =    burnt_area, & ! (m^2)
           &                     pe =    p(jc,nlev), & ! surface pressure
           &                  l_cpr =     l_cpr_bot)
      ! BB_plume expects fields to be ordered from surface to space
      IF (l_cpr_top) THEN
         CALL BB_plume(nlev =                      toplvl, &
              &      heat_f =               heat_flux_top, & ! (W m^-2)
              &        area =                  burnt_area, & ! (m^2)
              &       Te_in =       t(jc,mask1(1:toplvl)), & ! (K)
              &       pe_in =       p(jc,mask1(1:toplvl)), & ! (Pa)
              &          ue =       u(jc,mask1(1:toplvl)), & ! (m s^-1)
              &          ve =       v(jc,mask1(1:toplvl)), & ! (m s^-1)
              &      qve_in =      qv(jc,mask1(1:toplvl)), & ! (kg kg^-1)
              &     z_mc_in =    z_mc(jc,mask1(1:toplvl)), & ! (m)
              &    z_ifc_in = z_ifc(jc,mask2(1:toplvl+1)), & ! (m)
              &        ztop =                    ztop_pr)   ! (m)
      ELSE
         ztop_pr = 0._wp
      END IF

      IF (l_cpr_bot) THEN
         CALL BB_plume(nlev =                      toplvl, &
              &      heat_f =               heat_flux_bot, & ! (W m^-2)
              &        area =                  burnt_area, & ! (m^2)
              &       Te_in =       t(jc,mask1(1:toplvl)), & ! (K)
              &       pe_in =       p(jc,mask1(1:toplvl)), & ! (Pa)
              &          ue =       u(jc,mask1(1:toplvl)), & ! (m s^-1)
              &          ve =       v(jc,mask1(1:toplvl)), & ! (m s^-1)
              &      qve_in =      qv(jc,mask1(1:toplvl)), & ! (kg kg^-1)
              &     z_mc_in =    z_mc(jc,mask1(1:toplvl)), & ! (m)
              &    z_ifc_in = z_ifc(jc,mask2(1:toplvl+1)), & ! (m)
              &        ztop =                    zbot_pr)
      ELSE
         zbot_pr = 0._wp
      END IF

      ztop = ztop_pr
      zbot = zbot_pr

      ! keep the limiters from original plumerise
      IF (ztop <= 100._wp) ztop = 100._wp
      IF (zbot >= ztop)  zbot = ztop - 100._wp

      ! add topography
      ztop = ztop+topo
      zbot = zbot+topo


      CALL art_calc_wgt(ztop,zbot,nlev,z_ifc(jc,:),emiss_rate(jc,:),                     &
        &               dc_emis_res(jc,jb,CEILING(ls_time)),flux_bc(jc),area(jc),dz(jc,:) )

      deallocate(mask1,mask2)

    ENDIF ! flux_bc


  ENDDO !jc


END SUBROUTINE art_emission_biomBurn_bbplume
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_calc_wgt( ztop,zbot,nlev,z_ifc,emiss_rate,dc_emis_res,flux_bc,area,dz )

!>
!! SUBROUTINE art_calc_wgt
!! - computes weighting function for vertical distribution of smoke plume particles
!! -
!! Based on Walter et al. - The importance of plume rise on the concentrations and
!!                          atmospheric impacts of biomass burning aerosol
!! Part of Module: mo_art_emission_biomBurn
!! Author: Jonas Straub, KIT
!! Initial Release: 2018-02-21
!!
!! Modifications:
!! YYYY-MM-DD: <name>, <institution>
!! - ...
!<

  INTEGER,INTENT(IN)     :: &
    &  nlev                   !< Number of levels (equals index of lowest full level)

  REAL(wp),INTENT(INOUT)    :: &
    &  emiss_rate(:)           !< Emission rates (mug m-3 s-1), dimension: nlev

  REAL(wp),INTENT(IN)    :: &
    &  ztop,                & !< postprocessing corrected top height  (m)
    &  zbot,                & !< postprocessing correcteed bottom height  (m)
    &  z_ifc(:),            & !< geometric height at half levels (m)
    &  area,                & !< area of triangle (m2)
    &  dz(:),               & !< Height of model layer (m)
    &  dc_emis_res,         & !< (basic) diurnal cycle for emissions,
                              !  at current gridpoint and local time
    &  flux_bc                !< wildfires flux of black carbon, ~ FRP (kg/(m2 s)
  ! local
  INTEGER                :: &
    &  half,                & !< index for halflevel
    &  full,                & !< index for fulllevel
    &  jk                     !< loop index for vertical layer

  REAL(wp)               :: &
    &  ztopoflayer,         & !< loop variable, loop iterates from bottom of plume
                              !  t^3 top. ztopoflayer displays the currently uppermost layer
    &  zstartop,            & !< dimensionless height of upper grid cell level
    &  zstarbot               !< dimensionless height of lower grid cell level

  REAL(wp),POINTER       :: &
    &  wgt(:)                 !< vertical weighting function for emission distribution
  ! ----------------------------------------------------------------

  ! Check if there is a detected fire at this gridpoint
  IF ( dc_emis_res > 0._wp ) THEN
    ! -------------------------------------------
    ! --- 0. Allocate & initialize local storage
    ! -------------------------------------------
    ALLOCATE( wgt(nlev) )

    half        = 1
    ztopoflayer = 0._wp
    zstarbot    = 0._wp
    wgt(:)      = 0._wp

    ! -------------------------------------------
    ! --- 1. Find levels of plume
    ! -------------------------------------------
    ! find halflevel which is immediately below zbot
    DO WHILE (z_ifc(half) > zbot)
      half = half + 1
      IF (half > nlev) THEN
        EXIT
      ENDIF
    ENDDO

    ! half is now halflevel above zbot
    IF (half > 1) THEN
      half= half-1
    ENDIF

    ! zbot is within fulllevel full
    full = half
    ztopoflayer = z_ifc(half)

    ! -------------------------------------------
    ! --- 2. Computes vertical weighting function for particle concentration in plume
    ! -------------------------------------------
    ! from bottom to top of layer
    DO WHILE (ztopoflayer < ztop)
      zstartop = (ztopoflayer - zbot) / (ztop - zbot)
      IF (full > 0 .and. full <= nlev .and. size(wgt)>0) THEN
        wgt(full) = 3._wp*zstartop**2 - 2._wp*zstartop**3  &
          &       - (3._wp*zstarbot**2 - 2._wp*zstarbot**3)
      ENDIF
      zstarbot = zstartop
      half = half - 1
      full = full - 1
      ztopoflayer = z_ifc(half)
    ENDDO

    ! fullfill boundary condition
    zstartop  = 1._wp
    IF (full > 0 .and. full <= nlev .and. size(wgt)>0) THEN
      wgt(full) = 3._wp*zstartop**2 - 2._wp*zstartop**3  &
        &       - (3._wp*zstarbot**2 - 2._wp*zstarbot**3)
    ENDIF

    ! -------------------------------------------
    ! --- 3. outputrate
    ! -------------------------------------------
    DO jk=1,nlev
      ! = bc flux * vertical weighting * diurnal cycle * burning_area / cell_volume * 10e9 * daily emitted mass / hourly emitted mass
      emiss_rate(jk) = flux_bc * wgt(jk) * dc_emis_res / dz(jk) * 1.e9_wp * 24._wp
      ! in (mug m-3 s-1)
    ENDDO

    DEALLOCATE( wgt )

  ENDIF !dc_emiss_ret

END SUBROUTINE art_calc_wgt


END MODULE mo_art_emission_biomBurn
