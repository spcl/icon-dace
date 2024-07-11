!> QUINCY routines for shortwave surface radiation budget
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
!>#### routines for routines for the shortwave surface radiation budget, e.g., albedo and short-wave radiation
!>
MODULE mo_q_rad_process
#ifndef __NO_QUINCY__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_math_constants,    ONLY: eps8, eps1, seconds_per_day
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: q_radiation

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_rad_process'

CONTAINS

  ! ======================================================================================================= !
  !>calculates surface radiation reflection and absorption
  !>
  ! @TODO currently no soil and snow albedo calculation in the q_radiation routine
  SUBROUTINE q_radiation(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: A2L_, Q_RAD_, SPQ_, VEG_
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid
    USE mo_q_rad_parameters,      ONLY: rad2ppfd, def_alb_vis_soil, def_alb_nir_soil, rfr_ratio_toc, k_r2fr_chl, &
      &                                 def_alb_vis_snow, def_alb_nir_snow
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(Q_RAD_)
    dsl4jsb_Use_memory(SPQ_)
    dsl4jsb_Use_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER   :: model                  !< the model
    TYPE(t_lctlib_element), POINTER   :: lctlib                 !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER   :: vgrid_canopy_q_assimi  !< Vertical grid
    TYPE(t_jsb_vgrid),      POINTER   :: vgrid_snow_spq         !< Vertical grid
    INTEGER                           :: ncanopy                !< number of canopy layers
    INTEGER                           :: icanopy                !< loop 1,ncanopy
    INTEGER                           :: nsnow                  !< number of snow layers
    INTEGER                           :: iblk, ics, ice, nc     !< dimensions
    INTEGER                           :: ic                     !< looping
    REAL(wp)                          :: dtime                  !< timestep
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':q_radiation'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(Q_RAD_)
    dsl4jsb_Def_memory(SPQ_)
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Q_RAD_ 2D
    dsl4jsb_Real2D_onChunk :: sw_srf_net
    dsl4jsb_Real2D_onChunk :: swvis_srf_net
    dsl4jsb_Real2D_onChunk :: swnir_srf_net
    dsl4jsb_Real2D_onChunk :: rad_srf_net
    dsl4jsb_Real2D_onChunk :: alb_vis
    dsl4jsb_Real2D_onChunk :: alb_nir
    dsl4jsb_Real2D_onChunk :: alb_vis_soil
    dsl4jsb_Real2D_onChunk :: alb_nir_soil
    dsl4jsb_Real2D_onChunk :: alb_vis_snow
    dsl4jsb_Real2D_onChunk :: alb_nir_snow
    dsl4jsb_Real2D_onChunk :: alb_vis_can
    dsl4jsb_Real2D_onChunk :: alb_nir_can
    dsl4jsb_Real2D_onChunk :: arad_vis_soil
    dsl4jsb_Real2D_onChunk :: arad_nir_soil
    dsl4jsb_Real2D_onChunk :: arad_vis_can
    dsl4jsb_Real2D_onChunk :: arad_nir_can
    dsl4jsb_Real2D_onChunk :: appfd
    dsl4jsb_Real2D_onChunk :: rfr_ratio_boc
    dsl4jsb_Real2D_onChunk :: albedo
    dsl4jsb_Real2D_onChunk :: albedo_noon
    dsl4jsb_Real2D_onChunk :: net_radiation
    ! Q_RAD_ 3D
    dsl4jsb_Real3D_onChunk :: ppfd_sunlit_cl
    dsl4jsb_Real3D_onChunk :: ppfd_shaded_cl
    dsl4jsb_Real3D_onChunk :: arad_sunlit_vis_cl
    dsl4jsb_Real3D_onChunk :: arad_shaded_vis_cl
    dsl4jsb_Real3D_onChunk :: fleaf_sunlit_vis_cl
    dsl4jsb_Real3D_onChunk :: arad_sunlit_nir_cl
    dsl4jsb_Real3D_onChunk :: arad_shaded_nir_cl
    dsl4jsb_Real3D_onChunk :: fleaf_sunlit_nir_cl
    ! A2L_
    dsl4jsb_Real2D_onChunk :: swpar_srf_down
    dsl4jsb_Real2D_onChunk :: fract_par_diffuse
    dsl4jsb_Real2D_onChunk :: cos_zenith_angle
    dsl4jsb_Real2D_onChunk :: swnir_srf_down
    dsl4jsb_Real2D_onChunk :: swvis_srf_down
    dsl4jsb_Real2D_onChunk :: local_time_day_seconds
    ! SPQ_ 3D
    dsl4jsb_Real3D_onChunk :: w_snow_snl
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk :: lai
    ! VEG_ 3D
    dsl4jsb_Real3D_onChunk :: leaf_nitrogen_cl
    dsl4jsb_Real3D_onChunk :: fn_chl_cl
    dsl4jsb_Real3D_onChunk :: fleaf_sunlit_cl
    dsl4jsb_Real3D_onChunk :: lai_cl
    dsl4jsb_Real3D_onChunk :: cumm_lai_cl
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(Q_RAD_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! run only at PFT tiles (i.e. leafs of vegetation tile) because this routine is tailored towards vegetation only
    !  in the "usecase_pfts" Q_RAD_ is defined to run at leafs, which includes also other tiles than PFT
    IF (tile%lcts(1)%lib_id == 0) THEN
      IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Run only at PFT tiles. '//TRIM(tile%name)//' is not a PFT. ')
      RETURN
    ENDIF
    ! ----------------------------------------------------------------------------------------------------- !
    model                 => Get_model(tile%owner_model_id)
    lctlib                => model%lctlib(tile%lcts(1)%lib_id)
    vgrid_canopy_q_assimi => Get_vgrid('q_canopy_layer')
    ncanopy               =  vgrid_canopy_q_assimi%n_levels
    vgrid_snow_spq        => Get_vgrid('snow_layer_spq')
    nsnow                 =  vgrid_snow_spq%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(Q_RAD_)
    dsl4jsb_Get_memory(SPQ_)
    dsl4jsb_Get_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var2D_onChunk(A2L_, swpar_srf_down)           ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, fract_par_diffuse)        ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, cos_zenith_angle)         ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, swnir_srf_down)           ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, swvis_srf_down)           ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, local_time_day_seconds)   ! in
    ! ---------------------------
    ! Q_RAD_ 2D
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, sw_srf_net)               ! out
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, swvis_srf_net)            ! out
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, swnir_srf_net)            ! out
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, rad_srf_net)              ! out
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, alb_vis)                  ! out
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, alb_nir)                  ! out
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, alb_vis_soil)             ! out
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, alb_nir_soil)             ! out
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, alb_vis_snow)             ! out
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, alb_nir_snow)             ! out
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, alb_vis_can)              ! out
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, alb_nir_can)              ! out
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, arad_vis_soil)            ! out
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, arad_nir_soil)            ! out
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, arad_vis_can)             ! out
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, arad_nir_can)             ! out
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, appfd)                    ! out
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, rfr_ratio_boc)            ! out
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, albedo)                   ! out
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, albedo_noon)              ! out
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, net_radiation)            ! out
    ! Q_RAD_ 3D
    dsl4jsb_Get_var3D_onChunk(Q_RAD_, ppfd_sunlit_cl)           ! out
    dsl4jsb_Get_var3D_onChunk(Q_RAD_, ppfd_shaded_cl)           ! out
    dsl4jsb_Get_var3D_onChunk(Q_RAD_, arad_sunlit_vis_cl)       ! out
    dsl4jsb_Get_var3D_onChunk(Q_RAD_, arad_shaded_vis_cl)       ! out
    dsl4jsb_Get_var3D_onChunk(Q_RAD_, fleaf_sunlit_vis_cl)      ! out
    dsl4jsb_Get_var3D_onChunk(Q_RAD_, arad_sunlit_nir_cl)       ! out
    dsl4jsb_Get_var3D_onChunk(Q_RAD_, arad_shaded_nir_cl)       ! out
    dsl4jsb_Get_var3D_onChunk(Q_RAD_, fleaf_sunlit_nir_cl)      ! out
    ! ---------------------------
    ! SPQ_ 3D
    dsl4jsb_Get_var3D_onChunk(SPQ_, w_snow_snl)               ! in
    ! ---------------------------
    ! VEG_ 2D
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)                      ! in
    ! VEG_ 3D
    dsl4jsb_Get_var3D_onChunk(VEG_, leaf_nitrogen_cl)         ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_chl_cl)                ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_cl)          ! out
    dsl4jsb_Get_var3D_onChunk(VEG_, lai_cl)                   ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, cumm_lai_cl)              ! in
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 Snow and soil albedo
    !>
    ! for now, prescribed, should be based on JSBACH 4 routines, and include the effect of snow on the
    ! soil's albedo
    alb_vis_soil(:) = def_alb_vis_soil
    alb_nir_soil(:) = def_alb_nir_soil
    alb_vis_snow(:) = def_alb_vis_snow
    alb_nir_snow(:) = def_alb_nir_snow

    !> 2.0 calculate VIS radiation budget (absorbed radiation on canopy layers and soil + albedo)
    ! Note: shortwave diffuse fraction is assumed to be identical to PAR diffuse radiation
    CALL calc_sw_srf_net(nc, ncanopy, nsnow, &
                     lctlib%sigma_vis, lctlib%sigma_nir, lctlib%omega_clumping, lctlib%crown_shape_factor, &
                     "vis", &
                     swvis_srf_down(:), fract_par_diffuse(:), cos_zenith_angle(:), &
                     lai(:), lai_cl(:,:), cumm_lai_cl(:,:),&
                     alb_vis_soil(:), alb_vis_snow(:), w_snow_snl(:,:), &
                     arad_sunlit_vis_cl(:,:),arad_shaded_vis_cl(:,:),fleaf_sunlit_vis_cl(:,:), &
                     arad_vis_can(:),arad_vis_soil(:), &
                     alb_vis_can(:),alb_vis(:),swvis_srf_net(:))

    !> 3.0 calculate NIR radiation budget (absorbed radiation on canopy layers and soil + albedo)
    !!
    !! Following Weiss & Norman, Agri. For. Met. 1985, the fraction of VIS and NIR
    !! diffuse radiation are very similar
    ! Note: shortwave diffuse fraction is assumed to be identical to PAR diffuse radiation
    CALL calc_sw_srf_net(nc, ncanopy, nsnow, &
                     lctlib%sigma_vis, lctlib%sigma_nir, lctlib%omega_clumping, lctlib%crown_shape_factor, &
                     "nir", &
                     swnir_srf_down(:), fract_par_diffuse(:), cos_zenith_angle(:), &
                     lai(:), lai_cl(:,:), cumm_lai_cl(:,:), &
                     alb_nir_soil(:), alb_nir_snow, w_snow_snl(:,:), &
                     arad_sunlit_nir_cl(:,:),arad_shaded_nir_cl(:,:),fleaf_sunlit_nir_cl(:,:), &
                     arad_nir_can(:),arad_nir_soil(:), &
                     alb_nir_can(:),alb_nir(:),swnir_srf_net(:))

    !> 4.0 calculate albedo and net radiation for total shortwave radiation
    !!
    DO ic = 1, nc
      IF ((swvis_srf_down(ic) + swnir_srf_down(ic)) > eps8) THEN
        albedo(ic) = ( alb_vis(ic) * swvis_srf_down(ic) + alb_nir(ic) * swnir_srf_down(ic) ) &
          &           / (swvis_srf_down(ic) + swnir_srf_down(ic))
      ELSE
        ! albedo is not defined when radiation <= 0 and set to zero
        albedo(ic) = 0.0_wp
      ENDIF

      ! at noon: set albedo_noon to current albedo (noon equals half a day plus one timestep)
      IF (ABS(local_time_day_seconds(ic) - ((seconds_per_day / 2._wp) + dtime)) < eps1) THEN
        albedo_noon(ic) = albedo(ic)
      ENDIF
    END DO

    sw_srf_net(:) = swvis_srf_net(:) + swnir_srf_net(:)

    !> 5.0 diagnose net radiation (temporarily here for diagnostic purposes)
    !!
    !! Only for now longwave balance is zero
    rad_srf_net(:)   = sw_srf_net(:)
    net_radiation(:) = rad_srf_net(:)

    !> 6.0 convert radiation to Photosynthetically Active Photon Flux Density (ppfd)
    !!
    fleaf_sunlit_cl(:,:) = fleaf_sunlit_vis_cl(:,:)
    WHERE(swvis_srf_down(:) > eps8)
      appfd(:) = arad_vis_can(:) * swpar_srf_down(:)/swvis_srf_down(:) * rad2ppfd
    ELSEWHERE
      appfd(:) = 0.0_wp
    ENDWHERE

    DO icanopy=1,ncanopy
      DO ic = 1, nc
        IF (swvis_srf_down(ic) > eps8) THEN
          ppfd_sunlit_cl(ic,icanopy) = arad_sunlit_vis_cl(ic,icanopy) * swpar_srf_down(ic)/swvis_srf_down(ic) * rad2ppfd
          ppfd_shaded_cl(ic,icanopy) = arad_shaded_vis_cl(ic,icanopy) * swpar_srf_down(ic)/swvis_srf_down(ic) * rad2ppfd
        ELSE
          ppfd_sunlit_cl(ic,icanopy) = 0.0_wp
          ppfd_shaded_cl(ic,icanopy) = 0.0_wp
        END IF
      END DO
    END DO

    !> 7.0 calculate forest floor red to far-red ratio
    rfr_ratio_boc(:) = rfr_ratio_toc * &
                       exp( k_r2fr_chl * SUM(leaf_nitrogen_cl(:,:) * fn_chl_cl(:,:) * lai_cl(:,:), DIM=2))

  END SUBROUTINE q_radiation


  !-----------------------------------------------------------------------------------------------------
  !> calculate shortwave radiation budget of the tile
  !!
  !! @par info on routine
  !!  This routine has been adapted from ADFs implementation of the radiation calculations
  !!  by Spitters, C J T. 1986. \n Separating the Diffuse and Direct Component of Global Radiation
  !!  and Its Implications for Modeling Canopy Photosynthesis .2. Calculation of Canopy Photosynthesis.
  !!  Agricultural and Forest Meteorology 38 (1-3): 231-42. \n
  !!  Equations numbers refer to the Spitters paper \n
  !!  - The routine has been extended to diagnose albedo in the visible and NIR range. The implicit assumption
  !!  is made that albedo in the PAR and VIS range are identical \n
  !!  - The routine has been adjusted to a two-stream approximation scheme by adding the upward flux
  !!  from soil reflection \n
  !!  - The routine has been adjusted to account for canopy clumping (Campbell & Norman 1998) \n
  !! @todo: calc_sw_srf_net does not account non-leaf surfaces
  !!
  !! @par Spitters Eq. 12 is wrong
  !!  accordingly to Wang2003, "radf4" was replaced by "-kbl" when calculating "idrdra"
  !!
  !! @par canopy clumping factor
  !!  to disable clumping set lctlib(:)%omega_clumping = 1.0_wp
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_sw_srf_net( nc                  , &   ! in
                              ncanopy             , &   ! in
                              nsnow               , &   ! in
                              lctlib_sigma_vis    , &   ! in
                              lctlib_sigma_nir    , &   ! in
                              lctlib_omega_clumping     , &   ! in
                              lctlib_crown_shape_factor , &   ! in
                              vis_nir_id          , &   ! in
                              sw_srf_down         , &   ! in
                              fract_par_diffuse   , &   ! in
                              cos_zenith_angle    , &   ! in
                              lai                 , &   ! in
                              lai_cl              , &   ! in
                              cumm_lai_cl         , &   ! in
                              alb_soil            , &   ! in
                              alb_snow            , &   ! in
                              w_snow_snl          , &   ! in
                              arad_sunlit_cl      , &   ! out
                              arad_shaded_cl      , &   ! out
                              fleaf_sunlit_cl     , &   ! out
                              arad_can            , &   ! out
                              arad_soil           , &   ! out
                              alb_can             , &   ! out
                              alb                 , &   ! out
                              sw_srf_net)               ! out

    USE mo_q_rad_parameters,        ONLY: min_cos_zenith_angle_rad, kbl0_vis, kdf0_vis, kbl0_nir, kdf0_nir, rho2sbeta, &
                                          k_csf
    USE mo_spq_constants,           ONLY: w_snow_min
    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER                                , INTENT(in)  :: nc                 !< dimensions
    INTEGER                                , INTENT(in)  :: ncanopy            !< number of canopy layers
    INTEGER                                , INTENT(in)  :: nsnow
    REAL(wp)                               , INTENT(in)  :: lctlib_sigma_vis   !< lctlib parameter
    REAL(wp)                               , INTENT(in)  :: lctlib_sigma_nir   !< lctlib parameter
    REAL(wp)                               , INTENT(in)  :: lctlib_omega_clumping     !< lctlib parameter
    REAL(wp)                               , INTENT(in)  :: lctlib_crown_shape_factor !< lctlib parameter
    CHARACTER(len=3)                       , INTENT(in)  :: vis_nir_id         !< vis nir identifier for this function only
    REAL(wp), DIMENSION(nc)                , INTENT(in)  :: sw_srf_down        !< shortwave downward radiation flux (W m-2)
    REAL(wp), DIMENSION(nc)                , INTENT(in)  :: fract_par_diffuse
      !< fraction of diffuse radiation (--). Note: shortwave diffuse fraction is assumed to be identical to PAR diffuse radiation
    REAL(wp), DIMENSION(nc)                , INTENT(in)  :: cos_zenith_angle   !< solar angle (radians)
    REAL(wp), DIMENSION(nc)                , INTENT(in)  :: lai                !< leaf area index (m-2 m-2)
    REAL(wp), DIMENSION(nc,ncanopy)        , INTENT(in)  :: lai_cl             !< leaf area index of canopy layer (m-2 m-2)
    REAL(wp), DIMENSION(nc,ncanopy)        , INTENT(in)  :: cumm_lai_cl        !< cummulative leaf area index above centre of canopy layer
    REAL(wp), DIMENSION(nc)                , INTENT(in)  :: alb_soil           !< soil albedo
    REAL(wp), DIMENSION(nc)                , INTENT(in)  :: alb_snow           !< snow albedo
    REAL(wp), DIMENSION(nc,nsnow)          , INTENT(in)  :: w_snow_snl         !< snow
    REAL(wp), DIMENSION(nc,ncanopy)        , INTENT(out) :: arad_sunlit_cl     !< absorbed radiation on sunlit leaves (W m-2)
    REAL(wp), DIMENSION(nc,ncanopy)        , INTENT(out) :: arad_shaded_cl     !< absorbed radiation of shaded leaves (W m-2)
    REAL(wp), DIMENSION(nc,ncanopy)        , INTENT(out) :: fleaf_sunlit_cl    !< fraction of leaves that is sunlit (--)
    REAL(wp), DIMENSION(nc)                , INTENT(out) :: arad_can           !< absorbed radiation by the canopy (W m-2)
    REAL(wp), DIMENSION(nc)                , INTENT(out) :: arad_soil          !< absorbed radiation by the soil (W m-2)
    REAL(wp), DIMENSION(nc)                , INTENT(out) :: alb_can            !< albedo of the canopy
    REAL(wp), DIMENSION(nc)                , INTENT(out) :: alb                !< albedo of the tile
    REAL(wp), DIMENSION(nc)                , INTENT(out) :: sw_srf_net         !< net shortwave radiation of the tile
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc)                 :: tmp_cos_zenith_angle            ! local variable of solar angle (radians)
    REAL(wp), DIMENSION(nc)                 :: omega                           ! clumping factor
    REAL(wp), DIMENSION(nc)                 :: rho                             ! fraction reflected by green canopy (--)
    REAL(wp), DIMENSION(nc)                 :: kdf, kbl                        ! extinction coefficients for radiation over black canopy
    REAL(wp), DIMENSION(nc)                 :: irad_dir_down, &                ! downward direct and diffuse radiation
                                               irad_dif_down
    REAL(wp), DIMENSION(nc)                 :: irad_dir_down_bc, &             ! light intensity at the bottom of the canopy
                                               irad_dif_down_bc
    REAL(wp), DIMENSION(nc)                 :: irad_dif_up                     ! upward diffuse shortwave radiation
    REAL(wp), DIMENSION(nc)                 :: idfa                            ! diffuse radiation penetrating foliage at lc
    REAL(wp), DIMENSION(nc)                 :: idra                            ! direct radiation penetrating sunlit foliage at lc
    REAL(wp), DIMENSION(nc)                 :: idrdra                          ! direct radiation remaining direct at lc
    REAL(wp), DIMENSION(nc)                 :: arad_can_up                     ! absorbed radiation from the upward stream
    REAL(wp), DIMENSION(nc,ncanopy)         :: arad_shaded_up_cl               ! idem for canopy layers
    REAL(wp), DIMENSION(nc)                 :: radf1,radf2,radf3,radf4,radf5   ! dummy variables
    REAL(wp)                                :: temp1,temp2,radf3b              ! dummy variables
    REAL(wp)                                :: lctlib_sigma, &                 ! hlp for local lctlib variable, differentiating vis/nir
                                               kbl0, kdf0                      ! hlp for local rad_constants variables, differentiating vis/nir
    INTEGER                                 :: icanopy                         ! iteration over canopy layers
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_sw_srf_net'


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 0.9 set local variables according to vis/nir function call and initialise output
    !!
    SELECT CASE(vis_nir_id)
      CASE("vis")
        lctlib_sigma = lctlib_sigma_vis
        kbl0         = kbl0_vis                 ! extinction coefficient over black leaves
        kdf0         = kdf0_vis                 ! extinction coefficient for diffuse radiation
      CASE("nir")
        lctlib_sigma = lctlib_sigma_nir
        kbl0         = kbl0_nir                 ! extinction coefficient over black leaves
        kdf0         = kdf0_nir                 ! extinction coefficient for diffuse radiation
    END SELECT
    alb(:)                  = 0.0_wp
    alb_can(:)              = 0.0_wp
    arad_can(:)             = 0.0_wp
    arad_can_up(:)          = 0.0_wp    ! local var
    arad_soil(:)            = 0.0_wp
    arad_sunlit_cl(:,:)     = 0.0_wp
    arad_shaded_cl(:,:)     = 0.0_wp
    arad_shaded_up_cl(:,:)  = 0.0_wp    ! local var
    fleaf_sunlit_cl(:,:)    = 0.0_wp
    sw_srf_net(:)           = 0.0_wp


    !> 1.0 calculations
    !!

    !> enforce tmp_cos_zenith_angle to be larger than 'min_cos_zenith_angle_rad'
    tmp_cos_zenith_angle(:) = cos_zenith_angle(:)
    WHERE(cos_zenith_angle(:) < min_cos_zenith_angle_rad) tmp_cos_zenith_angle(:) = min_cos_zenith_angle_rad

    !> Eq.1 constants for albedo calculations
    radf3b = (1.0_wp - lctlib_sigma)
    temp1  = SQRT(radf3b)
    temp2  = 2.0_wp * (1.0_wp - temp1)/(1.0_wp + temp1)

    !> @par IF only do calculations if daylight is available, either direct or diffuse (in case of negative cos_zenith_angle)
    WHERE(sw_srf_down(:) > eps8)

       !> Downwelling direct radiation (Note: shortwave diffuse fraction is assumed to be identical to PAR diffuse radiation)
       irad_dir_down(:) = (1.0_wp - fract_par_diffuse(:)) * sw_srf_down(:)
       !> Downwelling diffuse radiation
       irad_dif_down(:) = fract_par_diffuse(:) * sw_srf_down(:)

       !> @par IF lai > eps8
       WHERE(lai(:) > eps8)

          !> Eq.1  Fraction of radiation reflected by green closed canopy (fraction)
          rho(:) = temp2/(1.0_wp + rho2sbeta * tmp_cos_zenith_angle(:))

          !> clumping factor (eq 15.35) Campbell & Norman (1998) as a function of sun-angle
          omega(:) = lctlib_omega_clumping / (lctlib_omega_clumping + (1.0_wp - lctlib_omega_clumping) * &
                  exp(-k_csf * ACOS(tmp_cos_zenith_angle(:))**lctlib_crown_shape_factor))

          !> Extinction coefficient for radiation over black leaves in canopy
          kdf(:)   = kdf0 * temp1 * omega(:) ! temp1 == SQRT(1.0_wp - lctlib_sigma)
          kbl(:)   = kbl0 * kdf(:)/(0.8_wp * temp1 * tmp_cos_zenith_angle(:))

          !> Eq.10 & 11 Some factors for radiation calculations that change with time
          radf1(:) = (1.0_wp - rho(:)) * irad_dif_down(:) * kdf(:)            ! eq 10
          radf2(:) = (1.0_wp - rho(:)) * temp1 * irad_dir_down(:) * kbl(:)    ! eq 11
          !> Eq.12 absorbed direct radiation
          radf3(:) = radf3b * irad_dir_down(:) * kbl(:)                    ! eq 12
          radf4(:) = -temp1 * kbl(:)

          !> light intensity at bottom of canopy
          irad_dif_down_bc(:) = (1.0_wp - rho(:)) * irad_dif_down(:) * exp(-kdf(:)  * lai(:))
          irad_dir_down_bc(:) = (1.0_wp - rho(:)) * irad_dir_down(:) * exp(radf4(:) * lai(:))

          !> diffuse upward shortwave flux from soil reflection and absorbed radiation of the soil
          !> first check if there is snow
          WHERE (w_snow_snl(:,1) < w_snow_min)
            irad_dif_up(:) = alb_soil(:) * (irad_dif_down_bc(:) + irad_dir_down_bc(:))
            arad_soil(:)   = (1.0_wp - alb_soil(:)) * (irad_dif_down_bc(:) + irad_dir_down_bc(:))
          ELSEWHERE
            irad_dif_up(:) = alb_snow(:) * (irad_dif_down_bc(:) + irad_dir_down_bc(:))
            arad_soil(:)   = (1.0_wp - alb_snow(:)) * (irad_dif_down_bc(:) + irad_dir_down_bc(:))
          ENDWHERE
          !> factor for the absorption calculation of the upward stream
          radf5(:) = (1.0_wp - rho(:)) * irad_dif_up(:) * kdf(:)            ! eq 10

          !> upward shortwave diffuse flux reaching the top of the canopy
          irad_dif_up(:) = (1.0_wp - rho(:)) * irad_dif_up(:) * exp(-kdf(:)*lai(:))

       ELSEWHERE
          !> case of no LAI: absorped radiation and upward flux from soil only
          WHERE (w_snow_snl(:,1) < w_snow_min)
            arad_soil(:)   = (1.0_wp - alb_soil(:) ) * (irad_dif_down(:) + irad_dir_down(:))
            irad_dif_up(:) = alb_soil(:)             * (irad_dif_down(:) + irad_dir_down(:))
          ELSEWHERE
            arad_soil(:)   = (1.0_wp - alb_snow(:) ) * (irad_dif_down(:) + irad_dir_down(:))
            irad_dif_up(:) = alb_snow(:)             * (irad_dif_down(:) + irad_dir_down(:))
          ENDWHERE
       ENDWHERE

    ENDWHERE

    !> @par loop over canopy layers. All fluxes are evaluated at the centre of the canopy layer
    DO icanopy = 1,ncanopy
       !> @par IF leaf area index of canopy layer > eps8: Integrate absorbed radiation over actually existing canopy layers
       WHERE(sw_srf_down(:) > eps8 .AND. lai_cl(:,icanopy) > eps8)

          !> @par IF radf2 <= eps8
          WHERE(radf2(:) <= eps8) ! case of only diffuse radiation
             !> Eq.10 Mean radiation penetrating foliage at lc
             arad_shaded_cl(:,icanopy)    = radf1(:) * exp(-kdf(:) * cumm_lai_cl(:,icanopy)) + &    ! eq 10
                                            radf5(:) * exp(-kdf(:) * (MAX(lai(:)-cumm_lai_cl(:,icanopy),0.0_wp)))
             !> radiation penetrating foliage from below for albedo cacumm_lai_cl(icanopy)ulation
             arad_shaded_up_cl(:,icanopy) = radf5(:) * exp(-kdf(:) * (MAX(lai(:)-cumm_lai_cl(:,icanopy),0.0_wp)))
          !> @par ELSE Case of direct and diffuse radiation
          ELSEWHERE
             !> Sunlit foliage at cumm_lai_cl(icanopy) (fraction)
             !! kbl: Extinction coefficient for radiation over black leaves in canopy
             fleaf_sunlit_cl(:,icanopy) = exp(-kbl(:) * cumm_lai_cl(:,icanopy))
             !> Eq.10 Diffuse radiation penetrating foliage at cumm_lai_cl(icanopy)
             idfa(:) = radf1(:) * exp(-kdf(:) * cumm_lai_cl(:,icanopy)) + &
                       radf5(:) * exp(-kdf(:) * (MAX(lai(:)-cumm_lai_cl(:,icanopy),0.0_wp)))
             !> radiation penetrating foliage from below for albedo accumulation ?
             arad_shaded_up_cl(:,icanopy) = radf5(:) * exp(-kdf(:) * (MAX(lai(:)-cumm_lai_cl(:,icanopy),0.0_wp)))
             !> Direct radiation penetrating sunlit foliage at cumm_lai_cl(icanopy)
             idra(:) = radf2(:) * exp(radf4(:) * cumm_lai_cl(:,icanopy))
             !> Eq.12 Mean direct radiation remaining direct at cumm_lai_cl(icanopy)
             !! radf4 replaced by -kbl according to Wang2003 (error in Spitters 1986 Eqn.12)
             idrdra(:) = radf3(:) * exp(-kbl(:) * cumm_lai_cl(:,icanopy))
             !> Radiation penetrating shaded foliage at cumm_lai_cl(icanopy) (umol/m2/s)
             !! Eq.13 Sum of diffuse radiation and diffuse component of direct flux
             arad_shaded_cl(:,icanopy) = idfa(:) + MAX((idra(:) - idrdra(:)),0.0_wp)
             !> Radiation penetrating sunlit foliage at lc (umol/m2/s)
             !! sum of diffuse and direct radiation
             arad_sunlit_cl(:,icanopy) = arad_shaded_cl(:,icanopy) + radf3(:)
          ENDWHERE

          !> total absorbed radiation
          arad_can(:) = arad_can(:) + ( arad_shaded_cl(:,icanopy) * (1.0_wp-fleaf_sunlit_cl(:,icanopy)) + &
                                        arad_sunlit_cl(:,icanopy) * fleaf_sunlit_cl(:,icanopy) )       * &
                                        lai_cl(:,icanopy)
          !> for albedo calculation: radiation absorbed from the upward stream
          arad_can_up(:) = arad_can_up(:) + arad_shaded_up_cl(:,icanopy) * lai_cl(:,icanopy)
       ENDWHERE
    ENDDO

    !> @par IF only do calculations if daylight is available, either direct or diffuse (in case of negative cos_zenith_angle)
    WHERE(sw_srf_down(:) > eps8)
       !> @par IF lai > eps8
       WHERE(lai(:) > eps8)
          !> albedo of the canopy, as seen from the top of the canopy
          !! inferred from the difference between top and bottom of the
          !! canopy light intensities minus the amount of absorbed downward radiation
          alb_can(:) = 1.0_wp - (irad_dif_down_bc(:) + irad_dir_down_bc(:) + (arad_can(:) - arad_can_up(:))) / &
                                (irad_dir_down(:) + irad_dif_down(:))

          !> total upward shortwave radiation flux at the top of the canopy:
          !! soil-based flux + reflectance of the canopy
          irad_dif_up(:) = irad_dif_up(:) + &
                           alb_can(:) * (irad_dir_down(:) + irad_dif_down(:) - (irad_dif_down_bc(:) + irad_dir_down_bc(:)))
       ENDWHERE

       !> albedo and net shortwave radiation of the tile
       alb(:)        = irad_dif_up(:) / (irad_dir_down(:) + irad_dif_down(:))
       sw_srf_net(:) = (1.0_wp - alb(:)) * sw_srf_down(:)

    ENDWHERE

  END SUBROUTINE calc_sw_srf_net

#endif
END MODULE mo_q_rad_process
