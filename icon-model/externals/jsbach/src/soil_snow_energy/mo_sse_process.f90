!> Contains the routines for the soil and snow energy processes
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_sse_process
#ifndef __NO_JSBACH__

  USE mo_jsb_physical_constants, ONLY: rhoh2o
  USE mo_kind,      ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: relative_humidity_soil, init_soil_temperature, &
    & calc_vol_heat_capacity, calc_thermal_conductivity, calc_soil_temperature, calc_snow_temperature, &
    & calc_snow_abcoeff, calc_soiltemp_old, Get_liquid_max

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_process'

CONTAINS

  ! --------------------------------------------------------------------------------------------------------- !
  !>
  !! Calculates the volumetric heat capacity of soil (de Vries, 1963)
  !!
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  FUNCTION calc_vol_heat_capacity( &
    & vol_heat_cap_soil_dry_mineral,              &
    & vol_heat_cap_soil_dry_organic,              &
    & porosity,                                   &
    & fract_organic,                              &
    & fract_water,                                &
    & fract_ice,                                  &
    & layer_depth,                                &
    & dz                                          &
    & ) RESULT(vol_heat_cap)

    !$ACC ROUTINE SEQ

    USE mo_sse_constants, ONLY: vol_hcap_water, vol_hcap_ice, vol_hcap_air, vol_hcap_bedrock

    REAL(wp), INTENT(in) :: porosity, vol_heat_cap_soil_dry_mineral, vol_heat_cap_soil_dry_organic, fract_organic, &
                         &  fract_water, fract_ice, layer_depth, dz
    REAL(wp)             :: vol_heat_cap

    IF (layer_depth > 0._wp) THEN  ! above bedrock
      vol_heat_cap = (1._wp - porosity)                                            & ! vol. fraction of solid (mineral + organic)
        &              * ( (1._wp - fract_organic) * vol_heat_cap_soil_dry_mineral & ! heat cap of mineral soil
        &                  +        fract_organic  * vol_heat_cap_soil_dry_organic & ! heat cap of organic soil
        &                )                                                         &
        &            + fract_water * vol_hcap_water                                & ! heat cap of liquid water
        &            + fract_ice   * vol_hcap_ice                                  & ! heat cap of frozen water
        ! @todo : porosity - fract_water - fract_ice should not become lower than zero.
        &            + MAX((porosity - fract_water - fract_ice),0._wp) * vol_hcap_air  ! heat cap of remaining free pores
      vol_heat_cap = ( vol_heat_cap * layer_depth + vol_hcap_bedrock * (dz - layer_depth) ) / dz ! soil and bedrock mixture
    ELSE
      vol_heat_cap = vol_hcap_bedrock
    END IF

  END FUNCTION calc_vol_heat_capacity

  ! --------------------------------------------------------------------------------------------------------- !
  !>
  !! Calculates the soil thermal conductivity after Johansen 1977.
  !!
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  FUNCTION calc_thermal_conductivity( &
    & heat_cond_soil_mineral,                        &
    & water, ice,                                    &
    & porosity,                                      &
    & fract_organic,                                 &
    & vol_porosity_org, hcond_org, hcond_dry_org,    &
    & layer_depth,                                   &
    & dz                                             &
    & ) RESULT(heat_cond)

    USE mo_jsb_physical_constants, ONLY: dens_soil
    USE mo_sse_constants,          ONLY: hcond_water, hcond_ice, hcond_bedrock


    REAL(wp), INTENT(in) :: heat_cond_soil_mineral, water, ice, porosity, &
                          & fract_organic, vol_porosity_org, hcond_org, hcond_dry_org, layer_depth, dz
    REAL(wp)             :: heat_cond

    REAL(wp) ::           &
      & dens_bulk,        & !< Bulk density of dry soil (including pores)
      & porosity_mineral, & !< Volumentric porosity of mineral soil part
      & hcond_dry_min,    & !< Dry thermal conductivity soil (including mineral soil and pores)
      & hcond_dry,        & !< Dry thermal conductivity of soil (including mineral and organic soil and pores)
      & sat_fact,         & !< Scaling factor to account for empty pore space
      & saturation,       & !< Degree of saturation
      & kersten,          & !< Kersten number
      & hcond_soil,       & !< Thermal conductivity of soil solids (mineral and organic)
      & hcond_sat,        & !< Thermal conductivity for saturated soil
      & fract_water,      & !< Effective fractional water volume for saturated soil computation
      & fract_ice           !< Effective fractional ice volume for saturated soil computation

    !$ACC ROUTINE SEQ

    IF (layer_depth > 0._wp) THEN  ! above bedrock
      ! Thermal conductivity of dry soil
      porosity_mineral = (porosity - fract_organic * vol_porosity_org) / (1._wp - fract_organic)
      dens_bulk = dens_soil * (1._wp - porosity_mineral)
      hcond_dry_min = (0.135_wp * dens_bulk + 64.7_wp) / (dens_soil - 0.947_wp * dens_bulk)
      hcond_dry = (1._wp - fract_organic) * hcond_dry_min + fract_organic * hcond_dry_org

      IF (water + ice <= EPSILON(1._wp)) THEN
        heat_cond = hcond_dry
      ELSE


        ! Thermal conductivity of soil solids
        hcond_soil = (1._wp - fract_organic) * heat_cond_soil_mineral + fract_organic * hcond_org

        ! Thermal conductivity of saturated soil
        saturation = (water + ice) / layer_depth / porosity  ! fractional volume of water/ice rel. to pore volume

        ! ! --- Alternate formulation ---
        ! ! @todo: in jsbach3, frozen soil is checked by t_soil_sl1 <= tmelt, however, the computations in calc_soil_temperature
        ! ! are such that some ice can remain in a layer even if t_soil_sl1 > tmelt.
        ! IF (ice > 0._wp) THEN  ! Frozen case
        !   kersten = saturation
        !   fract_water = water / layer_depth
        !   fract_ice   = porosity - fract_water
        !   hcond_sat =   hcond_soil      ** (1._wp-porosity) &
        !     &         * hcond_water     ** fract_water      &
        !     &         * hcond_ice       ** fract_ice
        ! ELSE                   ! Unfrozen case
        !   ! ! This is the case for fine unfrozen soil
        !   ! IF (saturation > 0.1_wp) THEN
        !   !   kersten = LOG10(saturation) + 1._wp
        !   ! ELSE
        !   !   kersten = 0._wp
        !   ! END IF
        !   ! This is the case for coarse unfrozen soil
        !   IF (saturation > 0.05_wp) THEN
        !     kersten = 0.7_wp * LOG10(saturation) + 1._wp
        !   ELSE
        !     kersten = 0._wp
        !   END IF
        !   fract_water = porosity
        !   hcond_sat =   hcond_soil      ** (1._wp-porosity) &
        !     &         * hcond_water     ** porosity
        ! END IF

        ! @todo: in jsbach3, frozen soil is checked by t_soil_sl1 <= tmelt, however, the computations in calc_soil_temperature
        ! are such that some ice can remain in a layer even if t_soil_sl1 > tmelt.
        IF (ice > 0._wp) THEN  ! Frozen case
          kersten = saturation
        ELSE                   ! Unfrozen case
          ! This is the case for fine unfrozen soil
          IF (saturation > 0.1_wp) THEN
            kersten = LOG10(saturation) + 1._wp
          ELSE
            kersten = 0._wp
          END IF
          ! ! This is the case for coarse unfrozen soil
          ! IF (saturation > 0.05_wp) THEN
          !   kersten = 0.7_wp * LOG10(saturation) + 1._wp
          ! ELSE
          !   kersten = 0._wp
          ! END IF
        END IF
        ! The scaling is changed compared to earlier versions. The current implementation interpretes the fractions
        ! of the individual components such that fract_water + fract_ice equals porosity, which is also supported by
        ! Nikoosokhan et al (2016), equation 4
        sat_fact    = MAX(0._wp, porosity / ((water + ice) / layer_depth)) ! Scaling factor to upscale water/ice to saturation
        fract_water = (water / layer_depth) * sat_fact ! Water fraction rel. to completely saturated soil
        fract_ice   = (  ice / layer_depth) * sat_fact ! Ice fraction rel. to completely saturated soil
        hcond_sat   = hcond_soil      ** (1._wp-fract_water-fract_ice) &
          &         * hcond_water     ** fract_water                   &
          &         * hcond_ice       ** fract_ice

        heat_cond = (hcond_sat - hcond_dry) * kersten + hcond_dry

        ! Finally, consider soil and bedrock mixture
        heat_cond = ( heat_cond * layer_depth + hcond_bedrock * (dz - layer_depth) ) / dz

        heat_cond = MAX(heat_cond, hcond_dry)
      END IF
    ELSE
      heat_cond = hcond_bedrock
    END IF

  END FUNCTION calc_thermal_conductivity

  ! --------------------------------------------------------------------------------------------------------- !
  !>
  !! Calculates the relative soil humidity wrt field capacity
  !!
  REAL(wp) FUNCTION relative_humidity_soil(ws, fc)

    !$ACC ROUTINE SEQ

  USE mo_jsb_math_constants, ONLY: pi

    REAL(wp), INTENT(in) :: &
      & ws,                 & !< moisture content [m], e.g. of uppermost soil layer
      fc                      !< field capacity (maximum moisture content) [m], e.g. of uppermost soil layer

    REAL(wp) :: moisture

    moisture = MIN(ws, fc)
    IF (moisture > 0._wp) THEN
      relative_humidity_soil = 0.5_wp * (1._wp - COS(moisture * pi / fc)) ! relative_humidity increases with higher moisture and
                                                                          !  lower fc. Note, fc has to be  >= than moisture!
    ELSE
      relative_humidity_soil = 0._wp
    END IF

  END FUNCTION relative_humidity_soil

  ! --------------------------------------------------------------------------------------------------------- !
  !>
  !! Initialize soil and surface temperature
  !!
  SUBROUTINE init_soil_temperature(cmid, tclim, tsoil)
    ! Method:
    !
    ! Starting from the tslclim field temperatures are set in relation
    ! to the depth of the soil layer and position of the initial
    ! day in the annual cycle.
    !
    ! tsl is at 0.07 m
    ! thickness of layers 0.065, 0.254, 0.913, 2.902, 5.700 m (s. soiltemp)
    !
    USE mo_jsb_time,       ONLY: t_datetime, get_time_start, get_year, deallocateDatetime, &
                                 & get_year_length, get_month_length, get_day_length, get_year_day
    USE mo_jsb_math_constants, ONLY: pi

    REAL(wp), INTENT(in)  :: cmid(:)     !< Depths of soil layer mid points
    REAL(wp), INTENT(in)  :: tclim(:,:)  !< Climatological monthly surface temperature
    REAL(wp), INTENT(out) :: tsoil(:,:)  !< Soil temperatures

    INTEGER                   :: nsoil, ndim
    REAL(wp)                  :: zkap, zsqrt
    REAL(wp), ALLOCATABLE     :: zrange(:), zdmax(:)
    INTEGER,  ALLOCATABLE     :: jmax(:)
    TYPE(t_datetime), POINTER :: inidate
    INTEGER                   :: iniyear
    INTEGER                   :: im, i, ic
    REAL(wp)                  :: year_len, day_in_year
    REAL(wp)                  :: month_len(12), month_mid(12)

    REAL(wp), PARAMETER :: cmid_offset = 0.07_wp

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_soil_temperature'

    ndim  = SIZE(tsoil,1)
    nsoil = SIZE(tsoil,2)

    ALLOCATE(zrange(ndim), jmax(ndim), zdmax(ndim))

    zkap = 7.5e-7_wp

    inidate => get_time_start()
    iniyear = get_year(inidate)

    year_len    = REAL(get_year_length(iniyear), wp)
    day_in_year = INT(get_year_day(inidate))
    CALL deallocateDatetime(inidate)

    DO im=1,12
      month_len(im) = REAL(get_month_length(iniyear, im), wp)
      month_mid(im) = month_len(im) * 0.5_wp
    END DO

    zsqrt = SQRT(zkap * year_len * get_day_length() / pi)

    ! Calendar month of maximum surface temperature
    jmax(:) = MAXLOC(tclim(:,:), DIM=2)

    ! Difference between months of maximum and minimum surface temperature
    zrange(:) = MAXVAL(tclim(:,:), DIM=2) - MINVAL(tclim(:,:), DIM=2)

    ! Calendar day of maximum surface temperature
    zdmax(:) = 0._wp
    DO ic=1,ndim
      DO im=1,jmax(ic)
        zdmax(ic) = zdmax(ic) + month_len(im)
      END DO
      zdmax(ic) = zdmax(ic) - month_mid(jmax(ic))
    END DO

    ! Initialize soil temperatures
    DO i=1,nsoil
      tsoil(:,i) = SUM(tclim(:,:), DIM=2) / 12._wp                            &
                 & + 0.5_wp * zrange(:) * EXP(-(cmid(i)-cmid_offset) / zsqrt) &
                 &   * COS(2._wp * pi * (day_in_year - zdmax(:)) / year_len - (cmid(i)-cmid_offset) / zsqrt)
    END DO

    DEALLOCATE(zrange, jmax, zdmax)

  END SUBROUTINE init_soil_temperature

  ! --------------------------------------------------------------------------------------------------------- !
  !>
  !! Calculates the soil layer temperatures, a/b coefficients, ground heat flux and ground heat capacity.
  !!
  SUBROUTINE calc_soil_temperature(       &
    & nc, nsoil, dz, mids,                &
    & delta_time, lstart,                 &
    & l_freeze, l_supercool,              &
    & t_soil_top,                         &
    & vol_heat_cap, heat_cond,            &
    & t_soil_sl,                          &
    & t_soil_acoef, t_soil_bcoef,         &
    & hcap_grnd, grnd_hflx,               &
    & thaw_depth,                         &
    ! optional
    & ws_max, matric_pot, bclapp, ws,     &
    & ice, frozen, melt, liquid_max       &
    & )

    USE mo_jsb_physical_constants, ONLY: tmelt, alf, rhoh2o, rhoi

    IMPLICIT NONE

    !-----------------------------------------------------------------------
    !  arguments

    INTEGER,  INTENT(in)    ::  nc                       !! length of the vector
    INTEGER,  INTENT(in)    ::  nsoil                    !! number of soil layers
    REAL(wp), INTENT(in)    ::  delta_time               !! Time step
    LOGICAL,  INTENT(in)    ::  lstart                   !! if T, start of experiment
    LOGICAL,  INTENT(in)    ::  l_freeze, l_supercool
    REAL(wp), INTENT(in)    ::  dz(:)                    !! soil layer thickness [m]
    REAL(wp), INTENT(in)    ::  mids(:)                  !! depth of mids of soil layers [m]
    REAL(wp), INTENT(in)    ::  t_soil_top(:)            !! surface temperature at top of soil [K]
    REAL(wp), INTENT(in)    ::  vol_heat_cap(:,:)        !! soil heat capacity [J/m^3K]
    REAL(wp), INTENT(in)    ::  heat_cond(:,:)           !! soil thermal conductivity [J/m/s/K]

    REAL(wp), INTENT(inout) ::  t_soil_sl(:,:)           !! soil temperature [K]
    REAL(wp), INTENT(inout) ::  t_soil_acoef(:,:)        !! soil A coefficient of Richtmyer and Morton scheme
    REAL(wp), INTENT(inout) ::  t_soil_bcoef(:,:)        !! soil B coefficient of Richtmyer and Morton scheme

    REAL(wp), INTENT(out)   ::  hcap_grnd(:)             !! heat capacity of the ground [J m-2 K-1]
    REAL(wp), INTENT(out)   ::  grnd_hflx(:)             !! ground heat flux [J m-2 s-1]
    REAL(wp), INTENT(out)   ::  thaw_depth(:)            !! Thawing depth [m]
    REAL(wp), INTENT(in), OPTIONAL    ::  ws_max(:,:)       !! Maximum water storage [m]
    REAL(wp), INTENT(in), OPTIONAL    ::  matric_pot(:,:)   !! Soil matric potential
    REAL(wp), INTENT(in), OPTIONAL    ::  bclapp(:,:)       !! Clapp and Hornberger parameter
    REAL(wp), INTENT(inout), OPTIONAL ::  ws(:,:)           !! Water storage [m]
    REAL(wp), INTENT(inout), OPTIONAL ::  ice(:,:)          !! Ice storage [m]
    REAL(wp), INTENT(out), OPTIONAL   ::  frozen(:,:)       !! Water flux from freezing soil water [kg m-2 s-1]
    REAL(wp), INTENT(out), OPTIONAL   ::  melt(:,:)         !! Water flux from melting soil ice    [kg m-2 s-1]
    REAL(wp), INTENT(out), OPTIONAL   ::  liquid_max(:,:)   !! Water that remains liquid at sub-zero temps [m]

    !------------------------------------------------------------------
    !  local Variables

    INTEGER  :: ic, is
    REAL(wp) :: z1(nc,nsoil)
    REAL(wp) :: zd1(nsoil-1)
    REAL(wp) :: zdz1(nc,nsoil-1), zdz2(nc,nsoil)
    REAL(wp) :: heat_cap(nc,nsoil)

    !-----------------------------------------------------------------------------------------------
    !  Computation of useful constants

    !$ACC DATA CREATE(z1, zd1, zdz1, zdz2, heat_cap)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO is = 1,nsoil-1
       zd1(is) = 1._wp / (mids(is+1) - mids(is))
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO ic = 1, nc
      DO is = 1, nsoil
        heat_cap(ic,is) = dz(is) * vol_heat_cap(ic,is)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !-----------------------------------------------------------------------------------------------
    !  Computation of soil temperaturs
    !    from surface temperatur and A and B coefficients computed at the previous time step
    !-----------------------------------------------------------------------------------------------

    IF (PRESENT(ws)) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO is = 1, nsoil
        DO ic = 1, nc
          frozen    (ic,is) = 0._wp
          melt      (ic,is) = 0._wp
          liquid_max(ic,is) = 0._wp
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      thaw_depth(ic) = 0._wp
    END DO
    !$ACC END PARALLEL LOOP

    IF (.NOT. lstart) THEN

      ! uppermost layer
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic=1,nc
        t_soil_sl(ic,1) = t_soil_top(ic)
      END DO
      !$ACC END PARALLEL LOOP

      ! deeper layers
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO is = 1, nsoil-1
        !$ACC LOOP GANG VECTOR
        DO ic=1,nc
          t_soil_sl(ic,is+1) = t_soil_acoef(ic,is) + t_soil_bcoef(ic,is) * t_soil_sl(ic,is)
        END DO
        !$ACC END LOOP
      END DO
      !$ACC END LOOP
      !$ACC END PARALLEL

      ! Phase change calculations
      IF (l_freeze .AND. PRESENT(ws)) THEN

        IF (l_supercool) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO is=1,nsoil
            DO ic=1,nc
              liquid_max(ic,is) = Get_liquid_max( &
                & t_soil_sl (ic,is), &
                & ws_max    (ic,is), &
                & matric_pot(ic,is), &
                & bclapp    (ic,is)  &
                & )
            END DO
          END DO
          !$ACC END PARALLEL LOOP
        ELSE
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO is=1,nsoil
            DO ic=1,nc
              liquid_max(ic,is) = 0._wp
            END DO
          END DO
          !$ACC END PARALLEL LOOP
        END IF

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO ic=1,nc
          DO is=1,nsoil
            IF (t_soil_sl(ic,is) < tmelt .AND. ws(ic,is) > liquid_max(ic,is)) THEN
              frozen(ic,is) = MIN(ws(ic,is) - liquid_max(ic,is), &
                &               heat_cap(ic,is) * (tmelt - t_soil_sl(ic,is)) / (alf * rhoh2o) &
                &              )
            ELSE IF (t_soil_sl(ic,is) > tmelt .AND. ice(ic,is) > 0._wp) THEN
              melt(ic,is) = MIN(ice(ic,is), heat_cap(ic,is) * (t_soil_sl(ic,is) - tmelt) / (alf * rhoi))
            END IF
          END DO
        END DO
        !$ACC END PARALLEL LOOP

        ! Note: only either frozen or melt is non-zero at any given point
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO is=1,nsoil
          DO ic=1,nc
            ws(ic,is) = ws(ic,is) - frozen(ic,is) + melt(ic,is)
            ice(ic,is) = ice(ic,is) + frozen(ic,is) - melt(ic,is)
            t_soil_sl(ic,is) = t_soil_sl(ic,is) + frozen(ic,is) * alf * rhoh2o / heat_cap(ic,is) - &
              & melt(ic,is) * alf * rhoi / heat_cap(ic,is)

            frozen(ic,is) = frozen(ic,is) * rhoh2o / delta_time
            melt  (ic,is) = melt  (ic,is) * rhoi   / delta_time
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      END IF

      ! Thawing depth
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO ic=1,nc
        thaw_depth(ic) = mids(nsoil)   ! Initialize with maximum thawing depth
        IF (t_soil_sl(ic,1) <= tmelt) THEN ! Top soil layer is frozen
          thaw_depth(ic) = 0._wp
        END IF
        !$ACC LOOP SEQ
        DO is=2,nsoil
          IF  (       thaw_depth(ic)  == mids(nsoil) & ! still the inital value
            &   .AND. t_soil_sl(ic,is-1) >  tmelt &
            &   .AND. t_soil_sl(ic,is  ) <= tmelt &
            & ) THEN
            thaw_depth(ic) = mids(is-1) + (mids(is)-mids(is-1)) * (t_soil_sl(ic,is-1) - &
              & tmelt) / (t_soil_sl(ic,is-1) - t_soil_sl(ic,is))
          END IF
        END DO
        !$ACC END LOOP
      END DO
      !$ACC END LOOP
      !$ACC END PARALLEL

    ELSE ! lstart
      IF (l_freeze .AND. PRESENT(ws)) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO is = 1, nsoil
          DO ic=1,nc
            liquid_max(ic,is) = 0._wp
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END IF

    !------------------------------------------------------------------------------------------------
    !  Computation of the Richtmyer and Morton A and B coefficients for the next time step
    !------------------------------------------------------------------------------------------------
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO is = 1, nsoil
      DO ic=1,nc
        zdz2(ic,is) = heat_cap(ic,is) / delta_time
        z1(ic,is)   = 0._wp
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO is = 1, nsoil-1
      DO ic=1,nc
        zdz1(ic,is) = zd1(is) * heat_cond(ic,is)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    ! lowest soil layer
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic=1,nc
      z1(ic,nsoil-1) = zdz2(ic,nsoil) + zdz1(ic,nsoil-1)
      t_soil_acoef(ic,nsoil-1) = zdz2(ic,nsoil) * t_soil_sl(ic,nsoil) / z1(ic,nsoil-1)
      t_soil_bcoef(ic,nsoil-1) = zdz1(ic,nsoil-1) / z1(ic,nsoil-1)
      END DO
    !$ACC END PARALLEL LOOP

    ! soil layers above
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO ic=1,nc
      !$ACC LOOP SEQ
      DO is = nsoil-1, 2, -1
        z1(ic,is-1) = 1._wp / (zdz2(ic,is) + zdz1(ic,is-1) + zdz1(ic,is) * (1._wp - t_soil_bcoef(ic,is)))
        t_soil_acoef(ic,is-1) = (t_soil_sl(ic,is) * zdz2(ic,is) + zdz1(ic,is) * t_soil_acoef(ic,is)) * z1(ic,is-1)
        t_soil_bcoef(ic,is-1) = zdz1(ic,is-1) * z1(ic,is-1)
      END DO
      !$ACC END LOOP
    END DO
    !$ACC END LOOP
    !$ACC END PARALLEL

    !------------------------------------------------------------------------------------------------
    ! Computation of surface diffusive heat flux from the ground and heat capacity of the ground
    !------------------------------------------------------------------------------------------------
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic=1,nc
      grnd_hflx(ic) = zdz1(ic,1) * (t_soil_acoef(ic,1) + (t_soil_bcoef(ic,1) - 1._wp) * t_soil_sl(ic,1))
      hcap_grnd(ic) = (zdz2(ic,1) * delta_time + delta_time * (1._wp - t_soil_bcoef(ic,1)) * zdz1(ic,1))
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE calc_soil_temperature

  ! --------------------------------------------------------------------------------------------------------- !
  !>
  !! Calculates the snow layer temperatures.
  !!
  SUBROUTINE calc_snow_temperature(         &
    & nc, nsnow,                            &
    & lstart,                               &
    & itop_old, itop,                       &
    & t_unfilt, t_soil_acoef, t_soil_bcoef, &
    & t_snow, t_snow_acoef, t_snow_bcoef,   &
    & t_soil_top)

    IMPLICIT NONE

    !-----------------------------------------------------------------------
    !  arguments

    INTEGER,  INTENT(in)    :: nc                         !! length of the vector
    INTEGER,  INTENT(in)    :: nsnow                      !! number of (fixed) snow layers
    LOGICAL,  INTENT(in)    :: lstart                     !! if T, start of experiment
    INTEGER,  INTENT(in)    :: itop(:), itop_old(:)       !! Top snow layer for previous and new time step
    REAL(wp), INTENT(in)    :: t_unfilt(:)                !! Surface temperature at top of snow [K]
    REAL(wp), INTENT(in)    :: t_soil_acoef(:)            !! Top layer soil A coefficient of Richtmyer and Morton scheme
    REAL(wp), INTENT(in)    :: t_soil_bcoef(:)            !! Top layer soil B coefficient of Richtmyer and Morton scheme

    REAL(wp), INTENT(inout) :: t_snow(:,:)                !! Snow temperature [K]
    REAL(wp), INTENT(inout) :: t_snow_acoef(:,:)          !! Snow A coefficient of Richtmyer and Morton scheme
    REAL(wp), INTENT(inout) :: t_snow_bcoef(:,:)          !! Snow B coefficient of Richtmyer and Morton scheme

    REAL(wp), INTENT(out)   :: t_soil_top(:)              !! New temperature of top soil layer

    !------------------------------------------------------------------
    !  local Variables

    INTEGER  :: ic, is

    ! Initialize snow layer temperatures with surface temperature
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO ic=1,nc
      DO is=1,nsnow
        t_snow(ic,is) = t_unfilt(ic)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !-----------------------------------------------------------------------------------------------
    !  Computation of new snow temperatures
    !-----------------------------------------------------------------------------------------------
    !
    IF (.NOT. lstart) THEN
      ! Initialize a/b coefficients for snow layers which did not exists in the last time step
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO ic=1,nc
        DO is=1,nsnow
          IF (itop(ic) <= is .AND. itop_old(ic) > nsnow) THEN         ! New snow layers on formerly snowfree soil
            t_snow_acoef(ic,is) = t_soil_acoef(ic)                    ! --> use old top soil coefficients
            t_snow_bcoef(ic,is) = t_soil_bcoef(ic)
          ELSE IF (itop(ic) <= is .AND. is <= itop_old(ic)) THEN      ! zero+ new snow layers: from new top layer to top prev layer
            t_snow_acoef(ic,is) = t_snow_acoef(ic,itop_old(ic))       ! --> use old top snow coefficients
            t_snow_bcoef(ic,is) = t_snow_bcoef(ic,itop_old(ic))
          END IF
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      ! Compute snow layer temperatures from a/b coefficients, except for the top snow layer which
      ! retains the surface temperature
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO ic=1,nc
        !$ACC LOOP SEQ
        DO is=2,nsnow
          IF (is > itop(ic)) THEN                                     ! Snow layers below new top layer
            t_snow(ic,is) = t_snow_acoef(ic,is-1) + t_snow_bcoef(ic,is-1) * t_snow(ic,is-1)
          END IF
        END DO
        !$ACC END LOOP
      END DO
      !$ACC END LOOP
      !$ACC END PARALLEL

      ! Compute new temperature for top soil layer ( or use the surface temperature where there was or is no snow)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic=1,nc
        IF (MAX(itop(ic), itop_old(ic)) <= nsnow) THEN
          t_soil_top(ic) = t_snow_acoef(ic,nsnow) + t_snow_bcoef(ic,nsnow) * t_snow(ic,nsnow)
        ELSE
          t_soil_top(ic) = t_unfilt(ic)
        END IF
      END DO
      !$ACC END PARALLEL LOOP

    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic=1,nc
        t_soil_top(ic) = t_unfilt(ic)
      END DO
      !$ACC END PARALLEL LOOP

    END IF

  END SUBROUTINE calc_snow_temperature

  ! --------------------------------------------------------------------------------------------------------- !
  !>
  !! Calculates the snow layer a/b coefficients, its ground heat flux, and heat capacity.
  !!
  SUBROUTINE calc_snow_abcoeff(           &
    & nc, nsnow,                          &
    & dz_soil,                            &
    & delta_time,                         &
    & dz, itop,                           &
    & vol_heat_cap, heat_cond,            &
    & vol_heat_cap_soil, heat_cond_soil,  &
    & t_soil_sl,                          &
    & t_soil_acoef, t_soil_bcoef,         &
    & t_snow, t_snow_acoef, t_snow_bcoef, &
    & hcap_grnd, grnd_hflx)

    IMPLICIT NONE

    !-----------------------------------------------------------------------
    !  arguments

    INTEGER,  INTENT(in)    :: nc                         !! length of the vector
    INTEGER,  INTENT(in)    :: nsnow                      !! number of (fixed) snow layers
    REAL(wp), INTENT(in)    :: delta_time                 !! Time step
    REAL(wp), INTENT(in)    :: dz(:,:)                    !! (dynamic) snow layer thickness [m]
    REAL(wp), INTENT(in)    :: dz_soil(2)                 !! Thickness of two uppermost soil layers [m]
    INTEGER,  INTENT(in)    :: itop(:)                    !! Top snow layer for previous and new time step
    REAL(wp), INTENT(in)    :: vol_heat_cap(:)            !! snow volumetric heat capacity [J/m^3K]
    REAL(wp), INTENT(in)    :: heat_cond(:)               !! snow thermal conductivity [J/m/s/K]
    REAL(wp), INTENT(in)    :: vol_heat_cap_soil(:)       !! top soil layer volumetric heat capacity
    REAL(wp), INTENT(in)    :: heat_cond_soil(:)          !! top soil layer heat conductivity
    REAL(wp), INTENT(in)    :: t_soil_sl(:,:)             !! Temperature of soil layers [K]
    REAL(wp), INTENT(in)    :: t_soil_acoef(:,:)          !! top soil layer A coefficient of Richtmyer and Morton scheme
    REAL(wp), INTENT(in)    :: t_soil_bcoef(:,:)          !! top soil layer B coefficient of Richtmyer and Morton scheme
    REAL(wp), INTENT(in)    :: t_snow(:,:)                !! snow temperature [K]

    REAL(wp), INTENT(out)   :: t_snow_acoef(:,:)          !! snow A coefficient of Richtmyer and Morton scheme
    REAL(wp), INTENT(out)   :: t_snow_bcoef(:,:)          !! snow B coefficient of Richtmyer and Morton scheme
    REAL(wp), INTENT(out)   :: hcap_grnd(:)               !! heat capacity of the top snow layer [J m-2 K-1]
    REAL(wp), INTENT(out)   :: grnd_hflx(:)               !! ground heat flux [J m-2 s-1]

    !------------------------------------------------------------------
    !  local Variables

    INTEGER  :: ic, is, itop_new !, itop_prev
    REAL(wp) :: zmid(nc,nsnow+2)
    REAL(wp) :: z1(nc)
    REAL(wp) :: zd1(nc,nsnow+1)
    REAL(wp) :: zdz1(nc,nsnow+1), zdz2(nc,nsnow+1)

    !$ACC DATA CREATE(zmid, z1, zd1, zdz1, zdz2)

    !$ACC PARALLEL LOOP GANG VECTOR ASYNC(1)
    DO ic=1,nc
      grnd_hflx(ic) = 0._wp
      hcap_grnd(ic) = 0._wp

      !-----------------------------------------------------------------------------------------------
      !  Computation of useful constants
      !
      ! dz and zmid are zero for unused snow layers
      !
      ! mid of snow layers
      zmid(ic,1) = dz(ic,1) / 2._wp
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO ic = 1, nc
      !$ACC LOOP SEQ
      DO is=2,nsnow
        zmid(ic,is) = zmid(ic,is-1) + 0.5_wp * (dz(ic,is-1) + dz(ic,is))
        IF (zmid(ic,is) > zmid(ic,is-1)) THEN
          zd1(ic,is-1) = 1._wp / (zmid(ic,is) - zmid(ic,is-1))
        END IF
      END DO
      !$ACC END LOOP
    END DO
    !$ACC END LOOP
    !$ACC END PARALLEL

    ! mid of two uppermost soil layers
    !$ACC PARALLEL LOOP GANG VECTOR ASYNC(1)
    DO ic=1,nc
      zmid(ic,nsnow+1) = zmid(ic,nsnow) + 0.5_wp * (dz(ic,nsnow) + dz_soil(1))
      zmid(ic,nsnow+2) = zmid(ic,nsnow+1) + 0.5_wp * (dz_soil(1) + dz_soil(2))
      zd1(ic,nsnow) = 1._wp / (zmid(ic,nsnow+1) - zmid(ic,nsnow))
      zd1(ic,nsnow+1) = 1._wp / (zmid(ic,nsnow+2) - zmid(ic,nsnow+1))
    END DO
    !$ACC END PARALLEL LOOP

    !
    !------------------------------------------------------------------------------------------------
    !  Computation of the Richtmyer and Morton A and B coefficients for the next time step
    !------------------------------------------------------------------------------------------------

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR PRIVATE(itop_new)
    DO ic=1,nc
      itop_new = itop(ic)
      !$ACC LOOP SEQ
      DO is = 1, itop_new - 1
        t_snow_acoef(ic,is) = 0._wp
        t_snow_bcoef(ic,is) = 0._wp
      END DO
      !$ACC END LOOP
    END DO
    !$ACC END LOOP
    !$ACC END PARALLEL

    !$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) VECTOR COLLAPSE(2) ASYNC(1)
    DO ic=1,nc
      DO is = 1,nsnow
        IF (is >= itop(ic)) THEN
          zdz2(ic,is) = vol_heat_cap(ic) * dz(ic,is) / delta_time
          zdz1(ic,is) = zd1(ic,is) * heat_cond(ic)
        END IF
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic=1,nc
      zdz2(ic,nsnow+1) = vol_heat_cap_soil(ic) * dz_soil(1) / delta_time
      zdz1(ic,nsnow+1) = zd1(ic,nsnow+1) * heat_cond_soil(ic)

      ! lowest snow layer
      IF (itop(ic) <= nsnow) THEN  ! Snow present
        z1(ic) = 1._wp / (zdz2(ic,nsnow+1) + zdz1(ic,nsnow) + zdz1(ic,nsnow+1) * (1._wp - t_soil_bcoef(ic,1)))
        t_snow_acoef(ic,nsnow) = (t_soil_sl(ic,1) * zdz2(ic,nsnow+1) + zdz1(ic,nsnow+1) * t_soil_acoef(ic,1)) * z1(ic)
        t_snow_bcoef(ic,nsnow) = zdz1(ic,nsnow) * z1(ic)
      ELSE
        t_snow_acoef(ic,nsnow) = 0._wp
        t_snow_bcoef(ic,nsnow) = 0._wp
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    ! snow layers above up to top snow layer
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO ic=1,nc
      !$ACC LOOP SEQ
      DO is = nsnow, 2, -1
        IF (is > itop(ic)) THEN
          z1(ic) = 1._wp / (zdz2(ic,is) + zdz1(ic,is-1) + zdz1(ic,is) * (1._wp - t_snow_bcoef(ic,is)))
          t_snow_acoef(ic,is-1) = (t_snow(ic,is) * zdz2(ic,is) + zdz1(ic,is) * t_snow_acoef(ic,is)) * z1(ic)
          t_snow_bcoef(ic,is-1) = zdz1(ic,is-1) * z1(ic)
        END IF
      END DO
      !$ACC END LOOP
    END DO
    !$ACC END LOOP
    !$ACC END PARALLEL

    !------------------------------------------------------------------------------------------------
    ! Computation of surface diffusive heat flux from the ground and heat capacity of the ground
    !------------------------------------------------------------------------------------------------
    !
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO ic=1,nc
      DO is=1,nsnow
        IF (itop(ic) == is) THEN ! Snow present
          grnd_hflx(ic) = zdz1(ic,is) * (t_snow_acoef(ic,is) + (t_snow_bcoef(ic,is) - 1._wp) * t_snow(ic,is))
          hcap_grnd(ic) = (zdz2(ic,is) * delta_time + delta_time * (1._wp - t_snow_bcoef(ic,is)) * zdz1(ic,is))
        END IF
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE calc_snow_abcoeff

  ! --------------------------------------------------------------------------------------------------------- !
  !>
  !! Legacy routine to calculate soil temperatures, a/b coefficient, ground heat flux, and ground heat capacity
  !!
  SUBROUTINE calc_soiltemp_old( &
    & nc,nsoil,cdel,cmid,     &
    & delta_time,               &
    & lstart,                   &
    & pts,                      &
    & psn,                      &
    & psocond, prgcgn,          &
    & pgrndc, pgrndd,           &
    & ptsoil,                   &
    & pgrndcapc, pgrndhflx,     &
    & ldglac                    &
    & )
  !
  !   AUTHOR:  FREDERIC HOURDIN     30/01/92
  !
  !            ADAPTED TO THE LMD-GCM BY JAN POLCHER  26/02/92
  !            ADAPTED TO THE ECHAM-GCM BY JAN-PETER SCHULZ, MPI  03/02/96
  !
  !            J.-P. SCHULZ   MPI - OCTOBER 1997 :
  !               ROUTINE USED FOR IMPLEMENTATION OF AN IMPLICIT
  !               COUPLING BETWEEN LAND SURFACE AND ATMOSPHERE IN THE
  !               ECHAM4 GCM.
  !            U.SCHLESE DKRZ - NOVEMBER 1999  MODIFIED FOR ECHAM5
  !            U.Schlese DKRZ - February 2000  new soil temperatures
  !            L Kornblueh, MPI, January 2003, removed MERGE
  !
  !            Adapted to JSBACH by Thomas Raddatz, Mai 2004
  !            Adapted to ICON  by Thomas Raddatz, Sep 2011
  !   OBJECTIVE:  COMPUTATION OF:
  !               THE GROUND TEMPERATURE EVOLUTION
  !               THE GROUND SPECIFIC HEAT "CAPCAL"
  !               THE SURFACE DIFFUSIVE FLUX FROM GROUND "F0"
  !
  !
  !   METHOD:  IMPLICIT TIME INTEGRATION
  !
  !   CONSECUTIVES GROUND TEMPERATURES ARE RELATED BY:
  !           T(K+1) = C(K) + D(K)*T(K)  (1)
  !   THE COEFFICIENTS C (=pgrndc) AND D (=pgrndd) ARE COMPUTED AT THE
  !   T-DT TIME-STEP.
  !   ROUTINE STRUCTURE:
  !   1)NEW TEMPERATURES ARE COMPUTED  USING (1)
  !   2)C AND D COEFFICIENTS ARE COMPUTED FROM THE NEW TEMPERATURE
  !     PROFILE FOR THE T+DT TIME-STEP
  !   3)THE COEFFICIENTS A AND B ARE COMPUTED WHERE THE DIFFUSIVE
  !     FLUXES AT THE T+DT TIME-STEP IS GIVEN BY
  !            FDIFF = A + B TS(T+DT)
  !     OR     FDIFF = F0 + CAPCAL (TS(T+DT)-TS(T))/DT
  !            WITH F0 = A + B (TS(T))
  !                 CAPCAL = B*DT
  !
  !     ------------------------------------------------------------------
  !
  !   DECLARATIONS:
  !
  !!$ TR USE mo_jsbach_constants   , ONLY: RhoH2O
  !!$ TR USE mo_time_control       , ONLY: lstart
  !
  !-----------------------------------------------------------------------
  !  ARGUMENTS
  !
    INTEGER, Intent(in)  ::  nc                   !! length of the vector
    INTEGER, Intent(in)  ::  nsoil                !! number of soil layers (fixed to 5)
    REAL(wp), Intent(in) ::  cdel(nsoil)
    REAL(wp), Intent(in) ::  cmid(nsoil)
    REAL(wp), Intent(in) ::  delta_time           !! time step
    LOGICAL, INTENT(in)  ::  lstart
    REAL(wp), Intent(in)     ::  pts(nc)            !! surface temperature at top of soil [K]
    REAL(wp), Intent(in)     ::  psn(nc)            !! equivalent snow depth [m water]
    REAL(wp), Intent(in)     ::  psocond(nc)        !! soil heat conductivity
    REAL(wp), Intent(in)     ::  prgcgn(nc)         !! soil heat capacity [J/m^3K]
    REAL(wp), Intent(inout)  ::  pgrndc(nc,nsoil)   !!
    REAL(wp), Intent(inout)  ::  pgrndd(nc,nsoil)   !!
    REAL(wp), Intent(inout)  ::  ptsoil(nc,nsoil)   !! soil temperature [K]
    REAL(wp), Intent(out)    ::  pgrndcapc(nc)      !!
    REAL(wp), Intent(out)    ::  pgrndhflx(nc)      !! ground heat flux
    LOGICAL, Intent(in)  ::  ldglac(nc)         !! glacier mask
  !
  !     ------------------------------------------------------------------
  !
  !  local Variables
  !
    INTEGER :: is, i
    REAL(wp) :: zso_cond(nc), zso_capa(nc)
    REAL(wp) :: z1(nc)
    REAL(wp) :: zd1(nsoil-1)
    REAL(wp) :: zdz1(nc,nsoil-1),   zdz2(nc,nsoil)
    REAL(wp) :: zkappa(nc,nsoil), zcapa(nc,nsoil)
    REAL(wp) :: zsnow_h(nc), zx1(nc), zx2(nc)
    REAL(wp) :: zrici, zdifiz, zsn_cond, zsn_dens, zsn_capa
  !
  !     ------------------------------------------------------------------
  !
  !*    1.  SPECIFYING THE DEPTHS OF THE TEMPERATURE LEVELS.
  !
  !*    1.1 SOME CONSTANTS.
  !
    zrici = 2.09e+06_wp                                 !! volumetric heat capacity of ice [j/m**3/k]
    zdifiz = 12.e-07_wp                                 !! temperature diffusivity of ice  [m**2/s]
    zsn_cond = 0.31_wp                                  !! snow thermal conductivity [j/s/m/k]
    zsn_dens = 330.0_wp                                 !! snow density              [kg/m**3]
    zsn_capa = 634500.0_wp                              !! snow  heat capacity   [j/m**3/k]
  !
  !*    1.2 COMPUTING SOME USEFUL CONSTANTS.
  !
    !$ACC DATA CREATE(zso_cond, zso_capa, z1, zd1, zdz1, zdz2, zkappa, zcapa, zsnow_h, zx1, zx2)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO is = 1,nsoil-1
       zd1(is) = 1._wp / (cmid(is+1) - cmid(is))
    END DO
    !$ACC END PARALLEL LOOP
  !
  !*    1.3 COMPUTE OF THE SOIL THERMAL CONDUCTIVITY [J/S/M/K] FROM
  !*        THE SOIL TEMPERATURE DIFFUSIVITY [M**2/S].
  !
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO i = 1, nc
      IF (ldglac(i)) THEN
        zso_capa(i) = zrici
        zso_cond(i) = zso_capa(i) * zdifiz
      ELSE
        zso_capa(i) = prgcgn(i)
        zso_cond(i) = psocond(i)
      END IF
    END DO
    !$ACC END PARALLEL LOOP
  !
  !*    1.4 PRE-SET THERMAL CONDUCTIVITY AT ALL LEVELS.
  !
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO i = 1, nc
      DO is = 1,nsoil
        zkappa(i,is) = zso_cond(i)
        zcapa(i,is)  = zso_capa(i)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

  !
  !   --------------------------------------------------------------
  !   COMPUTATION OF THE GROUND TEMPERATURES USING THE CGRD AND DGRD
  !   COEFFICIENTS COMPUTED AT THE PREVIOUS TIME-STEP
  !   --------------------------------------------------------------
  !
  !   Upper layer
  !
    IF (.NOT. lstart) THEN

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO i = 1, nc
        ptsoil(i,1) = pts(i)
      END DO
      !$ACC END PARALLEL LOOP
  !
  !   Deeper layers
  !

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO i = 1, nc
        DO is = 1,nsoil-1
          ptsoil(i,is+1) = pgrndc(i,is) + pgrndd(i,is) * ptsoil(i,is)
        END DO
      END DO
    !$ACC END PARALLEL LOOP

    END IF
  !
  !   ---------------------------------------------------------------
  !   COMPUTATION OF THE CGRD AND DGRD COEFFICIENTS FOR THE NEXT STEP
  !   ---------------------------------------------------------------
  !
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO i = 1, nc
      zsnow_h(i) = psn(i) * RhoH2O / zsn_dens
  !
  !*       Special treatment for first layer
  !
      IF ( zsnow_h(i) > cmid(2) ) THEN
        zcapa(i,1) = zsn_capa
        zkappa(i,1) = zsn_cond
      ELSE IF( zsnow_h(i) > 0.0_wp .AND. zsnow_h(i) <= cmid(2) ) THEN
        zx1(i) = zsnow_h(i) / cmid(2)
        zx2(i) = ( cmid(2) - zsnow_h(i)) / cmid(2)
        zcapa(i,1) = zx1(i) * zsn_capa + zx2(i) * zso_capa(i)
        zkappa(i,1) = 1._wp / ( zx1(i) / zsn_cond + zx2(i) / zso_cond(i) )
      ELSE
        zcapa(i,1) = zso_capa(i)
        zkappa(i,1) = zso_cond(i)
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO i = 1, nc
      !$ACC LOOP SEQ
      DO is = 2, nsoil - 2
        IF ( zsnow_h(i) > cmid(is+1) ) THEN
          zcapa(i,is) = zsn_capa
          zkappa(i,is) = zsn_cond
        ELSE IF ( zsnow_h(i) > cmid(is) .AND. zsnow_h(i) <= cmid(is+1) ) THEN
          zx1(i) = (zsnow_h(i) - cmid(is)) * zd1(is)
          zx2(i) = ( cmid(is+1) - zsnow_h(i)) * zd1(is)
          zcapa(i,is) = zx1(i) * zsn_capa + zx2(i) * zso_capa(i)
          zkappa(i,is) = 1._wp / ( zx1(i) / zsn_cond + zx2(i) / zso_cond(i) )
        ELSE
          zcapa(i,is) = zso_capa(i)
          zkappa(i,is) = zso_cond(i)
        END IF
      END DO
    END DO
    !$ACC END PARALLEL
  !
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO i = 1, nc
      DO is=1,nsoil
        zdz2(i,is) = zcapa(i,is) * cdel(is) / delta_time
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO i = 1, nc
      DO is=1,nsoil-1
        zdz1(i,is) = zd1(is) * zkappa(i,is)
      END DO
    END DO
    !$ACC END PARALLEL LOOP
    !
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO i = 1, nc
      z1(i) = zdz2(i,nsoil) + zdz1(i,nsoil-1)
      pgrndc(i,nsoil-1) = zdz2(i,nsoil) * ptsoil(i,nsoil) / z1(i)
      pgrndd(i,nsoil-1) = zdz1(i,nsoil-1) / z1(i)
    END DO
    !$ACC END PARALLEL LOOP
    !

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO i = 1, nc
      !$ACC LOOP SEQ
      DO is=nsoil-1,2,-1
        z1(i) = 1._wp / (zdz2(i,is) + zdz1(i,is-1) + zdz1(i,is) * (1._wp - pgrndd(i,is)))
        pgrndc(i,is-1) = (ptsoil(i,is) * zdz2(i,is) + zdz1(i,is) * pgrndc(i,is)) * z1(i)
        pgrndd(i,is-1) = zdz1(i,is-1) * z1(i)
      END DO
      !$ACC END LOOP
    END DO
    !$ACC END LOOP
    !$ACC END PARALLEL
  !
  !   ---------------------------------------------------------
  !   COMPUTATION OF THE SURFACE DIFFUSIVE FLUX FROM GROUND AND
  !   CALORIFIC CAPACITY OF THE GROUND:
  !   ---------------------------------------------------------
  !
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO i = 1, nc
      pgrndhflx(i) = zdz1(i,1) * (pgrndc(i,1) + (pgrndd(i,1) - 1._wp) * ptsoil(i,1))
      pgrndcapc(i) = (zdz2(i,1) * delta_time + delta_time * (1._wp - pgrndd(i,1)) * zdz1(i,1))
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE calc_soiltemp_old

  ! --------------------------------------------------------------------------------------------------------- !
  !>
  !! Calulates potentially supercooled water within soil layers
  !!
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  FUNCTION Get_liquid_max( &
    & temp, w_max, matric_pot, bclapp     &
    & ) RESULT(liquid_max)

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY: tmelt, alf, grav

    REAL(wp), INTENT(in) :: &
      & temp, w_max, matric_pot, bclapp
    REAL(wp) :: liquid_max

    IF (temp < tmelt) THEN
      ! supercooled water equation (Niu&Yang,2006)
      liquid_max =                                                               &
        & w_max                                                                  & ! Maximum water storage
        & * ( alf * (tmelt - temp) / MAX(1.e-10_wp,  grav * temp * (-matric_pot)) ) &
        &   ** ( -1._wp / MAX(1._wp, bclapp))
    ELSE
      liquid_max = 0._wp
    END IF

  END FUNCTION Get_liquid_max

#endif
END MODULE mo_sse_process
