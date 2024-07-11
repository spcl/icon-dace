!> QUINCY soil-turbulence calculation
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### calculate the task update_soil_turbulence, i.e., surface condition & turbulence, drag coefficient
!>
MODULE mo_q_spq_update_turbulence
#ifndef __NO_QUINCY__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_soil_turbulence

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_spq_update_turbulence'

CONTAINS

  ! ======================================================================================================= !
  !>update soil moisture and theta for the plant only model - NOTE update/improve docu
  !>
  ! two cases considered:
  !   a) dynamic water uptake calculation
  !   b) prescribed values from namelist
  SUBROUTINE update_soil_turbulence(tile, options)

    USE mo_jsb_class,              ONLY: Get_model
    USE mo_jsb_process_class,      ONLY: A2L_, SPQ_, VEG_
    USE mo_jsb_tile_class,         ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,         ONLY: t_jsb_task_options
    USE mo_jsb_model_class,        ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,       ONLY: t_lctlib_element
    USE mo_jsb_math_constants,     ONLY: dtime
    USE mo_atmland_constants,      ONLY: min_wind, max_wind
    USE mo_atmland_util,           ONLY: calc_spec_humidity_sat
    USE mo_jsb_physical_constants, ONLY: grav, r_gas_vapour, r_gas_dryair, cpd, cpvd1, rd_o_cpd, rvd1
    USE mo_spq_constants,          ONLY: height_wind, w_snow_min
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(SPQ_)
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(SPQ_)
    dsl4jsb_Use_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER       :: model                 !< the model
    TYPE(t_lctlib_element), POINTER       :: lctlib                !< land-cover-type library - parameter across pft's
    REAL(wp), DIMENSION(options%nc)       :: ztvd, ztvir, q_sat_surf, ztvs, zg, zgh, zdu2, zril, zcdnl, zchln, &
                                             zscfl, zucfl, zucfhl, zcons, zcfncl, zcfnchl, zcfml, zcfhl, &
                                             csat, cair, z0m, z0h, zchnl, zlev, wind_air
    REAL(wp), DIMENSION(options%nc)       :: fract_fpc, blended_height
    REAL(wp)                              :: zcons8, zcons9, zcons11, zcons12
    INTEGER                               :: iblk, ics, ice, nc,ic !< dimensions
    REAL(wp)                              :: steplen               !< ...
    REAL(wp)                              :: alpha                 !< implicitness factor
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_soil_turbulence'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(SPQ_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(SPQ_)
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_
    dsl4jsb_Real2D_onChunk                :: t_air
    dsl4jsb_Real2D_onChunk                :: q_air
    dsl4jsb_Real2D_onChunk                :: wind_10m
    dsl4jsb_Real2D_onChunk                :: press_srf
    dsl4jsb_Real2D_onChunk                :: drag_srf
    dsl4jsb_Real2D_onChunk                :: ga
    dsl4jsb_Real2D_onChunk                :: pch
    dsl4jsb_Real2D_onChunk                :: t_acoef
    dsl4jsb_Real2D_onChunk                :: t_bcoef
    dsl4jsb_Real2D_onChunk                :: q_acoef
    dsl4jsb_Real2D_onChunk                :: q_bcoef
    ! SPQ_ 2D
    dsl4jsb_Real2D_onChunk                :: fact_q_air
    dsl4jsb_Real2D_onChunk                :: fact_qsat_srf
    ! SPQ_ 3D
    dsl4jsb_Real3D_onChunk                :: t_soil_sl
    dsl4jsb_Real3D_onChunk                :: w_snow_snl
    dsl4jsb_Real3D_onChunk                :: t_snow_snl
    ! VEG_
    dsl4jsb_Real2D_onChunk                :: height
    dsl4jsb_Real2D_onChunk                :: lai
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    steplen = options%steplen
    alpha   = options%alpha
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(SPQ_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts(1)%lib_id)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(SPQ_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(SPQ_)
    dsl4jsb_Get_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)                ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, q_air)                ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, wind_10m)             ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, press_srf)            ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, drag_srf)             ! out
    dsl4jsb_Get_var2D_onChunk(A2L_, ga)                   ! out
    dsl4jsb_Get_var2D_onChunk(A2L_, pch)                  ! out
    dsl4jsb_Get_var2D_onChunk(A2L_, t_acoef)              ! out
    dsl4jsb_Get_var2D_onChunk(A2L_, t_bcoef)              ! out
    dsl4jsb_Get_var2D_onChunk(A2L_, q_acoef)              ! out
    dsl4jsb_Get_var2D_onChunk(A2L_, q_bcoef)              ! out
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SPQ_, fact_q_air)           ! in
    dsl4jsb_Get_var2D_onChunk(SPQ_, fact_qsat_srf)        ! in
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SPQ_, t_soil_sl)            ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, t_snow_snl)            ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, w_snow_snl)            ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, height)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)                  ! in
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.8 correct height for bare soil effect on roughness
    !>
    ! note: the '* 1.0_wp' for calc of blended_height(:) is the minimum/default height in meter
    !       using the VEG_ min_height parameter (= 0.1_wp) instead would change simulation results
    fract_fpc(:)      = 1.0_wp - EXP(- 0.5_wp * lai(:))
    blended_height(:) = height(:) * fract_fpc(:) + (1._wp - fract_fpc(:)) * 1.0_wp

    !>0.9 init local variables and constants
    !>
    ! Some constants from echam/mo_surface_land.f90
    csat(:)     = fact_qsat_srf(:)
    cair(:)     = fact_q_air(:)
    wind_air(:) = MIN(MAX(wind_10m(:),min_wind),max_wind)
    zlev(:)     = height_wind + 3._wp/10._wp * height(:)
    z0m(:)      = MAX(blended_height(:) / 10.0_wp, 0.1_wp)
    z0h(:)      = z0m(:) / 10.0_wp

    ! geopotential of the surface layer (see echam's auxhybc.f90 & geopot.f90)
    ! If height_wind is set, then the measurement height + an offset from the average vegetation height
    ! is used. Otherwise, the code defaults to the half-level of ECHAM's lowest atmospheric layer
    zg(:)       = zlev(:) * grav
    zgh(:)      = zg(:)                          ! height_humidity and height_temperature have to be equal
    zdu2(:)     = wind_air(:)**2._wp

    !>1.0 ...
    !>
    ! virtual potential air temperature (see mo_surface_boundary.f90)
    ! according to Saucier, WJ Principles of Meteoroligical Analyses
    ! tv = t * (1 + 0.61 * q * t)    ! virtual temperature
    ! td = t * ( 1000 / p_mb ) ^ R/cdp  ! potential temperature
    ! tvd = tair * (100000/p_pa)^rd_o_cpd * 1 + rvd1 * q_air) ! virtual potential temperature
    ztvd(:)       = t_air(:) * ( 100000._wp/press_srf(:))**rd_o_cpd * ( 1._wp + rvd1 * q_air(:) )
    ztvir(:)      = t_air(:) * ( 1._wp + rvd1 * q_air(:) )
    ! virtual potential surface temperature
    WHERE (w_snow_snl(:,1) > w_snow_min)
      q_sat_surf(:) = calc_spec_humidity_sat(t_snow_snl(:,1),press_srf(:))
      ztvs(:)       = t_snow_snl(:,1) * ( 100000._wp/press_srf(:))**rd_o_cpd * &
                    ( 1._wp + rvd1 * ( csat(:) * q_sat_surf(:) + ( 1._wp - cair(:) ) * q_air(:)))
    ELSEWHERE
      q_sat_surf(:) = calc_spec_humidity_sat(t_soil_sl(:,1),press_srf(:))
      ztvs(:)       = t_soil_sl(:,1) * ( 100000._wp/press_srf(:))**rd_o_cpd * &
                    ( 1._wp + rvd1 * ( csat(:) * q_sat_surf(:) + ( 1._wp - cair(:) ) * q_air(:)))
    ENDWHERE

    !>2.0 ...
    !>
    ! Richardson number (dry, Brinkop & Roeckner 1995, Tellus)
    ! ztvd, ztvs are now virtual potential temperatures, changed by Thomas Raddatz 07.2014
    ! SZ: limited zril to < 1.5 to avoid complete decoupling
    ! effectively this limits drag to values larger than about 10% of the neutral drag
    zril(:)       = MIN(zg(:) * ( ztvd(:) - ztvs(:) ) / ( zdu2(:) * (ztvd(:) + ztvs(:)) / 2._wp ), 0.5_wp)

    DO ic = 1, nc
      CALL calc_drag_coefficient(steplen, &
                                alpha, &
                                t_air(ic),q_air(ic),wind_air(ic),press_srf(ic), &
                                zg(ic),zgh(ic),z0m(ic),z0h(ic), &
                                zril(ic), &                 ! inout
                                drag_srf(ic),pch(ic),ga(ic))  ! out
    END DO

    ! Richtmeyer-Morton coefficient for explicit coupling
    t_acoef(:) = 0.0_wp
    t_bcoef(:) = t_air(:) * cpd * (1._wp + cpvd1 * q_air(:)) + zgh
    q_acoef(:) = 0.0_wp
    q_bcoef(:) = q_air(:)
  END SUBROUTINE update_soil_turbulence

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_soil_turbulence
  !
  !-----------------------------------------------------------------------------------------------------
  !> Routine to calculate drag coefficient for a given surface condition and turbulence
  !!
  !! Calculates drag coefficient and aerodynamic conductance
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_drag_coefficient( &
                         steplen, &
                         alpha, &
                         t_air,q_air,wind_air,press_srf, &
                         zg,zgh,z0m,z0h,zril, &
                         drag_srf,pch,ga)

    USE mo_phy_schemes,             ONLY: heat_transfer_coef
    USE mo_jsb_physical_constants,  ONLY: grav, r_gas_vapour, r_gas_dryair, von_karman, rvd1

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp),                        INTENT(in)    :: steplen                !< step length for this process
    REAL(wp),                        INTENT(in)    :: alpha                  !< implicitness factor
    REAL(wp),                        INTENT(in)    :: t_air, &
                                                      q_air, &
                                                      wind_air, &
                                                      press_srf, &
                                                      zg, &
                                                      zgh, &
                                                      z0m, &
                                                      z0h, &
                                                      zril
    REAL(wp),                        INTENT(out)   :: drag_srf, &
                                                      pch, &
                                                      ga
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)                              :: zcdnl, zchln, &
                                             zscfl, zucfl, zucfhl, zcons, zcfncl, zcfnchl, zcfml, zcfhl, &
                                             zchnl,zdu2
    REAL(wp)                              :: zcons8,zcons9,zcons11, zcons12
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_drag_coefficient'


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    ! 0.9 init local variables and constants
    zcons8     = 2._wp * 5.0_wp
    zcons9     = 3._wp * 5.0_wp
    zcons11    = 3._wp * 5.0_wp * 5.0_wp
    zcons12    = alpha * steplen * grav / r_gas_dryair
    zdu2       = wind_air**2._wp

    ! Neutral drag coefficient for momentum and heat
    zcdnl = (von_karman / LOG(1._wp + zg / (grav *z0m )))**2._wp
    zchnl = von_karman**2._wp / (LOG(1._wp + zg / (grav * z0m )) * &
                           LOG( ( grav * z0m + zgh ) / (grav * z0h )) )

    ! account for stable/unstable case: helper variables
    zscfl  = SQRT ( 1._wp + 5._wp * ABS(zril))
    zucfl  = 1._wp / (1._wp + zcons11 * zcdnl * SQRT(ABS(zril) * (1._wp  + zg / (grav * z0m))))
    zucfhl = 1._wp / (1._wp + zcons11 * zchnl * SQRT(ABS(zril) * (1._wp  + zgh / (grav * z0h))))

    ! ignoring cloud water correction (see mo_surface_land.f90)
    zcons    = zcons12 * press_srf / ( t_air * (1._wp + rvd1 * q_air))
    zcfncl   = zcons * SQRT(zdu2) * zcdnl
    zcfnchl  = zcons * SQRT(zdu2) * zchnl

    ! Stable / Unstable case
    IF ( zril > 0._wp ) THEN
       zcfml = zcfncl  / (1._wp + zcons8 * zril / zscfl)
       zcfhl = zcfnchl / (1._wp + zcons9 * zril * zscfl)
    ELSE
       zcfml = zcfncl  * (1._wp - zcons8 * zril * zucfl)
       zcfhl = zcfnchl * (1._wp - zcons9 * zril * zucfhl)
    END IF

    ! total surface drag (ECHAM)
    drag_srf = zcfhl
    pch      = zcfhl / zcfnchl * zchnl
    ! aerodynamic conductance (m/s)
    ga       = heat_transfer_coef(drag_srf, steplen, alpha) * wind_air
  END SUBROUTINE calc_drag_coefficient

#endif
END MODULE mo_q_spq_update_turbulence
