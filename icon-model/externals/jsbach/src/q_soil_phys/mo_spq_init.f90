!> QUINCY soil-physics variables init
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
!>#### initialization of soil-physics-quincy memory variables using, e.g., ic & bc input files
!>
MODULE mo_spq_init
#ifndef __NO_QUINCY__

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: message_text, finish, message
  USE mo_jsb_control,             ONLY: debug_on
  USE mo_jsb_math_constants,      ONLY: eps8
  USE mo_jsb_process_class,       ONLY: SPQ_

#ifdef __QUINCY_STANDALONE__
  USE mo_qs_process_init_util,    ONLY: init_spq_soil_properties_spp1685_sites
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: spq_init

  TYPE t_spq_init_vars
    REAL(wp), POINTER ::                &
      & soil_sand               (:,:) => NULL(), &
      & soil_silt               (:,:) => NULL(), &
      & soil_clay               (:,:) => NULL(), &
      & soil_sat_water_content  (:,:) => NULL(), &
      & soil_bulk_density       (:,:) => NULL()
  END TYPE t_spq_init_vars

  TYPE(t_spq_init_vars) :: spq_init_vars

  CHARACTER(len=*), PARAMETER :: modname = 'mo_spq_init'

CONTAINS

  ! ======================================================================================================= !
  !> Intialize SPQ_ process
  !>
  SUBROUTINE spq_init(tile)

    USE mo_jsb_tile_class,         ONLY: t_jsb_tile_abstract

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':spq_init'

#ifdef __QUINCY_STANDALONE__
    CALL spq_qs_read_init_vars(tile)
    CALL spq_init_ic_bc(tile)
    CALL spq_finalize_init_vars()
#else
    ! call only at root tile, i.e., the only tile without associated parent tile (to avoid unnessary i/o)
    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      CALL spq_read_init_vars(tile)
    END IF

    CALL spq_init_ic_bc(tile)

    IF (tile%Is_last_process_tile(SPQ_)) THEN
      CALL spq_finalize_init_vars()
    END IF
#endif
  END SUBROUTINE spq_init

  ! ======================================================================================================= !
  !> Intialize SPQ_ process from ic and bc input files
  !>
  SUBROUTINE spq_init_ic_bc(tile)
    USE mo_jsb_class,              ONLY: Get_model
    USE mo_jsb_tile_class,         ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,        ONLY: t_jsb_model
    USE mo_quincy_model_config,    ONLY: QLAND, QPLANT, QSOIL, QCANOPY, Q_TEST_CANOPY, Q_TEST_RADIATION
    USE mo_jsb_grid_class,         ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,               ONLY: Get_grid, Get_vgrid
    USE mo_util,                   ONLY: soil_depth_to_layers_2d
    USE mo_jsb_physical_constants, ONLY: Tzero
    USE mo_spq_constants,          ONLY: soil_heat_cap, soil_therm_cond, kdiff_sat_max, &
      &                                  k_pwp_s, k_pwp_c, k_pwp_sc, k_pwp_a, k_pwp_at, k_pwp_bt, &
      &                                  k_fc_s, k_fc_c, k_fc_sc, k_fc_a, k_fc_at, k_fc_bt, k_fc_ct, &
      &                                  k_sat_s, k_sat_c, k_sat_sc, k_sat_a, k_sat_at, k_sat_bt, k_sat_ct, k_sat_dt, &
      &                                  wpot_pwp, wpot_fc, wpot_fc
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(SPQ_)
    dsl4jsb_Use_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER         :: model              !< the model
    TYPE(t_jsb_grid),       POINTER         :: hgrid              !< Horizontal grid
    TYPE(t_jsb_vgrid),      POINTER         :: vgrid_soil_sb      !< Vertical grid
    TYPE(t_jsb_vgrid),      POINTER         :: vgrid_snow_spq     !< Vertical grid snow
    INTEGER                                 :: nsoil_sb           !< number of soil layers as used/defined by the SB_ process
    INTEGER                                 :: nsnow              !< number of snow layers
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)   :: soil_awc
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)   :: hlp2
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: hlp1
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: theta_pwp_sl
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: theta_fc_sl
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: theta_sat_sl
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: sand_spq_init_sl
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: silt_spq_init_sl
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: clay_spq_init_sl
    CHARACTER(len=3)                        :: site_ID_spp1685
    INTEGER                                 :: is, ic, iblk
    INTEGER                                 :: nproma, nblks
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':spq_init_ic_bc'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(SPQ_)
    dsl4jsb_Def_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! SPQ_ 2D
    dsl4jsb_Real2D_onDomain    :: soil_depth
    dsl4jsb_Real2D_onDomain    :: root_depth
    dsl4jsb_Real2D_onDomain    :: t_srf_new
    dsl4jsb_Real2D_onDomain    :: t_srf_old
    dsl4jsb_Real2D_onDomain    :: temp_srf_eff_4
    dsl4jsb_Real2D_onDomain    :: qsat_star
    dsl4jsb_Real2D_onDomain    :: s_star
    dsl4jsb_Real2D_onDomain    :: fact_q_air
    dsl4jsb_Real2D_onDomain    :: fact_qsat_srf
    dsl4jsb_Real2D_onDomain    :: z0h
    dsl4jsb_Real2D_onDomain    :: z0m
    dsl4jsb_Real2D_onDomain    :: w_skin
    dsl4jsb_Real2D_onDomain    :: w_soil_root
    dsl4jsb_Real2D_onDomain    :: w_soil_root_pwp
    dsl4jsb_Real2D_onDomain    :: w_soil_root_fc
    dsl4jsb_Real2D_onDomain    :: w_soil_root_theta
    dsl4jsb_Real2D_onDomain    :: w_soil_root_pot
    dsl4jsb_Real2D_onDomain    :: num_sl_above_bedrock
    ! SPQ_ 3D
    dsl4jsb_Real3D_onDomain    :: soil_depth_sl          ! (soil_depth_sl == soil_lay_width_sl)
    dsl4jsb_Real3D_onDomain    :: soil_lay_width_sl
    dsl4jsb_Real3D_onDomain    :: soil_lay_depth_center_sl
    dsl4jsb_Real3D_onDomain    :: soil_lay_depth_ubound_sl
    dsl4jsb_Real3D_onDomain    :: soil_lay_depth_lbound_sl
    dsl4jsb_Real3D_onDomain    :: bulk_dens_sl
    dsl4jsb_Real3D_onDomain    :: heat_capa_sl
    dsl4jsb_Real3D_onDomain    :: therm_cond_sl
    dsl4jsb_Real3D_onDomain    :: w_soil_pwp_sl
    dsl4jsb_Real3D_onDomain    :: w_soil_fc_sl
    dsl4jsb_Real3D_onDomain    :: w_soil_sat_sl
    dsl4jsb_Real3D_onDomain    :: w_soil_sl
    dsl4jsb_Real3D_onDomain    :: w_soil_pot_sl
    dsl4jsb_Real3D_onDomain    :: w_ice_sl
    dsl4jsb_Real3D_onDomain    :: saxtonA
    dsl4jsb_Real3D_onDomain    :: saxtonB
    dsl4jsb_Real3D_onDomain    :: saxtonC
    dsl4jsb_Real3D_onDomain    :: kdiff_sat_sl
    dsl4jsb_Real3D_onDomain    :: sand_sl
    dsl4jsb_Real3D_onDomain    :: silt_sl
    dsl4jsb_Real3D_onDomain    :: clay_sl
    dsl4jsb_Real3D_onDomain    :: volume_min_sl
    dsl4jsb_Real3D_onDomain    :: t_snow_snl
    dsl4jsb_Real3D_onDomain    :: t_soil_sl
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(SPQ_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model          => Get_model(tile%owner_model_id)
    hgrid          => Get_grid(model%grid_id)
    vgrid_soil_sb  => Get_vgrid('soil_layer_sb')
    vgrid_snow_spq => Get_vgrid('snow_layer_spq')
    nblks          =  hgrid%nblks
    nproma         =  hgrid%nproma
    nsoil_sb       =  vgrid_soil_sb%n_levels
    nsnow          =  vgrid_snow_spq%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(SPQ_)
    dsl4jsb_Get_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    ALLOCATE(soil_awc(nproma, nblks))
    ALLOCATE(hlp2(nproma, nblks))
    ALLOCATE(hlp1(nproma, nsoil_sb, nblks))
    ALLOCATE(theta_pwp_sl(nproma, nsoil_sb, nblks))
    ALLOCATE(theta_fc_sl(nproma, nsoil_sb, nblks))
    ALLOCATE(theta_sat_sl(nproma, nsoil_sb, nblks))
    ALLOCATE(sand_spq_init_sl(nproma, nsoil_sb, nblks))
    ALLOCATE(clay_spq_init_sl(nproma, nsoil_sb, nblks))
    ALLOCATE(silt_spq_init_sl(nproma, nsoil_sb, nblks))
    ! ----------------------------------------------------------------------------------------------------- !
    ! SPQ_ 2D
    dsl4jsb_Get_var2D_onDomain(SPQ_, soil_depth)
    dsl4jsb_Get_var2D_onDomain(SPQ_, root_depth)
    dsl4jsb_Get_var2D_onDomain(SPQ_, t_srf_new)
    dsl4jsb_Get_var2D_onDomain(SPQ_, t_srf_old)
    dsl4jsb_Get_var2D_onDomain(SPQ_, temp_srf_eff_4)
    dsl4jsb_Get_var2D_onDomain(SPQ_, qsat_star)
    dsl4jsb_Get_var2D_onDomain(SPQ_, s_star)
    dsl4jsb_Get_var2D_onDomain(SPQ_, fact_q_air)
    dsl4jsb_Get_var2D_onDomain(SPQ_, fact_qsat_srf)
    dsl4jsb_Get_var2D_onDomain(SPQ_, z0h)
    dsl4jsb_Get_var2D_onDomain(SPQ_, z0m)
    dsl4jsb_Get_var2D_onDomain(SPQ_, w_skin)
    dsl4jsb_Get_var2D_onDomain(SPQ_, w_soil_root)
    dsl4jsb_Get_var2D_onDomain(SPQ_, w_soil_root_pwp)
    dsl4jsb_Get_var2D_onDomain(SPQ_, w_soil_root_fc)
    dsl4jsb_Get_var2D_onDomain(SPQ_, w_soil_root_theta)
    dsl4jsb_Get_var2D_onDomain(SPQ_, w_soil_root_pot)
    dsl4jsb_Get_var2D_onDomain(SPQ_, num_sl_above_bedrock)
    ! SPQ_ 3D
    dsl4jsb_Get_var3D_onDomain(SPQ_, soil_depth_sl)              ! out   (soil_depth_sl == soil_lay_width_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, soil_lay_width_sl)          ! out
    dsl4jsb_Get_var3D_onDomain(SPQ_, soil_lay_depth_center_sl)   ! out
    dsl4jsb_Get_var3D_onDomain(SPQ_, soil_lay_depth_ubound_sl)   ! out
    dsl4jsb_Get_var3D_onDomain(SPQ_, soil_lay_depth_lbound_sl)   ! out
    dsl4jsb_Get_var3D_onDomain(SPQ_, bulk_dens_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, heat_capa_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, therm_cond_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, w_soil_pwp_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, w_soil_fc_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, w_soil_sat_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, w_soil_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, w_soil_pot_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, w_ice_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, saxtonA)
    dsl4jsb_Get_var3D_onDomain(SPQ_, saxtonB)
    dsl4jsb_Get_var3D_onDomain(SPQ_, saxtonC)
    dsl4jsb_Get_var3D_onDomain(SPQ_, kdiff_sat_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, sand_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, silt_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, clay_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, volume_min_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, t_snow_snl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, t_soil_sl)
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.9 snow init
    !>
    t_snow_snl(:,:,:) = Tzero

    !>1.0 Total soil depth until bedrock
    !>
    ! QUINCY soil_depth for one site (read from namelist)
    soil_depth(:,:) = dsl4jsb_Config(SPQ_)%soil_depth
    ! for all points of the domain soil_depth(:,:) is at least 0.1
    soil_depth(:,:) = MAX(0.1_wp, MERGE(soil_depth(:,:), 0._wp, soil_depth(:,:) >= 0._wp))

    !>2.0 calc actual soil layer depths based on the fixed layer thicknesses from namelist and the bedrock depth
    !>
    ! NOTE: soil_depth_to_layers_2d() returns a 3D variable
    soil_depth_sl(:,:,:) = soil_depth_to_layers_2d(soil_depth(:,:), &     ! Total soil depth until bedrock (from textures), 2D variable
      &                                            vgrid_soil_sb%dz(:))   ! Soil layer thicknesses from namelist
    ! pass these values to a SPQ_ variable with a correct name
    soil_lay_width_sl(:,:,:) = soil_depth_sl(:,:,:)

    !>  2.1 calc more metrics of the soil layers
    !>
    ! i)  lower & upper boundary
    ! ii) depth at the center of each layer
    DO iblk = 1,nblks
      DO is = 1,nsoil_sb
        ! lower & upper boundary
        IF (is == 1) THEN
          soil_lay_depth_lbound_sl(:, is, iblk) = 0.0_wp
          soil_lay_depth_ubound_sl(:, is, iblk) = soil_lay_width_sl(:, is, iblk)
        ELSE
          soil_lay_depth_lbound_sl(:, is, iblk) = soil_lay_depth_lbound_sl(:, is-1, iblk) &
            &                                        + soil_lay_width_sl(:, is-1, iblk)
          soil_lay_depth_ubound_sl(:, is, iblk) = soil_lay_depth_ubound_sl(:, is-1, iblk) &
            &                                        + soil_lay_width_sl(:, is, iblk)
        ENDIF
        ! soil-layer center
        soil_lay_depth_center_sl(:, is, iblk) = (soil_lay_depth_lbound_sl(:, is, iblk) &
          &                                        + soil_lay_depth_ubound_sl(:, is, iblk)) &
          &                                        * 0.5_wp
      END DO
    END DO

    !>   2.2 get number of soil layers above bedrock for each gridcell
    !>
    !>    use this variable for looping over soil layers,
    !>    excluding soil layers with a width smaller/equal eps8
    !>
    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1,nsoil_sb
          IF (soil_lay_width_sl(ic, is, iblk) > eps8) THEN
            num_sl_above_bedrock(ic, iblk) = REAL(is, wp)
          ELSE
            EXIT  ! stop looping over further soil layers with thickness < eps8
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    !> 3.0 init soil-property variables and land2atmosphere variables
    !>

    !>   3.1 init with one site-specific value for all soil layers
    !>
    bulk_dens_sl(:,:,:)   = SPREAD(spq_init_vars%soil_bulk_density(:,:), DIM = 2, ncopies = nsoil_sb)
    heat_capa_sl(:,:,:)   = soil_heat_cap
    therm_cond_sl(:,:,:)  = soil_therm_cond
    sand_sl(:,:,:)        = SPREAD(spq_init_vars%soil_sand(:,:), DIM = 2, ncopies = nsoil_sb)
    silt_sl(:,:,:)        = SPREAD(spq_init_vars%soil_silt(:,:), DIM = 2, ncopies = nsoil_sb)
    clay_sl(:,:,:)        = SPREAD(spq_init_vars%soil_clay(:,:), DIM = 2, ncopies = nsoil_sb)
    volume_min_sl(:,:,:)  = 1._wp
    ! land2atmosphere
    !   these variables are initialized with the mem%Add_var() routine in addition (same values) !
    t_soil_sl(:,:,:)      = 273.15_wp
    t_srf_new(:,:)        = 273.15_wp
    t_srf_old(:,:)        = 273.15_wp
    temp_srf_eff_4(:,:)   = 273.15_wp ** 4.0_wp
    qsat_star(:,:)        = 0.0075_wp   ! value from JSBACH
    s_star(:,:)           = 2.9E5_wp    ! value from JSBACH
    fact_q_air(:,:)       = 0.5_wp      ! value from JSBACH
    fact_qsat_srf(:,:)    = 0.5_wp      ! value from JSBACH
    z0h(:,:)              = 1.0_wp      ! value from JSBACH
    z0m(:,:)              = 1.0_wp      ! value from JSBACH

#ifdef __QUINCY_STANDALONE__
    !>     3.1.1 soil-layer specific init for sites of the "SPP 1685" site-set
    !>
    !>     works only with 15 soil layer
    !>
    IF (model%config%flag_spp1685) THEN
      DO iblk = 1,nblks
        CALL init_spq_soil_properties_spp1685_sites( &
          & nproma, &                           ! in
          & nsoil_sb, &
          & model%config%lon, &                 ! in ('site_longitude')
          & bulk_dens_sl(:, :, iblk), &         ! out
          & sand_sl(:, :, iblk), &
          & silt_sl(:, :, iblk), &
          & clay_sl(:, :, iblk) )               ! out
      END DO
    END IF
#endif

    !>   3.2 ensure a certain minimum proportion for sand, clay and silt
    !>    needed for the Saxton equations
    !>
    ! this is to ensure applying the empirical functions by Saxton within the range of the texture they were designed for \n
    ! avoiding ranges where calculation would return negative theta_pwp_sl or relatively large kdiff_sat_sl values
    !
    ! Saxton equations with too small soil-texture fractions may lead to unrealistic percolation rates and soil water potentials
    !
    ! for more details see ticket #560
    !
    ! using the modified *_spq_init_sl soil texture for "saxton functions here in spq_init_ic_bc" only
    sand_spq_init_sl(:,:,:) = sand_sl(:,:,:)
    clay_spq_init_sl(:,:,:) = clay_sl(:,:,:)
    silt_spq_init_sl(:,:,:) = silt_sl(:,:,:)
    DO iblk = 1,nblks
      CALL calc_soil_texture_saxton_compatible( &
        & nproma, &                         ! in
        & nsoil_sb, &                       ! in
        & sand_spq_init_sl(:, :, iblk), &   ! inout
        & silt_spq_init_sl(:, :, iblk), &   ! inout
        & clay_spq_init_sl(:, :, iblk) )    ! inout
    END DO

    !> 4.0 soil hydraulic properties, based on Saxton and Rawls 2006, SSSA
    !>

    !>   4.1 fractional water holding capacity at permanent wilting point
    !>   NOTE that this assumes that the soil hydraulic properties are identical between the layers
    !>
    hlp1(:,:,:)         = k_pwp_s * sand_spq_init_sl(:,:,:) + &
                          k_pwp_c * clay_spq_init_sl(:,:,:) + &
                          k_pwp_sc * sand_spq_init_sl(:,:,:) * clay_spq_init_sl(:,:,:) + &
                          k_pwp_a
    theta_pwp_sl(:,:,:) = k_pwp_at + (1.0_wp + k_pwp_bt) * hlp1(:,:,:)

    !>   4.2 fractional water holding capacity at field capacity
    !>
    hlp1(:,:,:)         = k_fc_s * sand_spq_init_sl(:,:,:) + &
                          k_fc_c * clay_spq_init_sl(:,:,:) + &
                          k_fc_sc * sand_spq_init_sl(:,:,:) * clay_spq_init_sl(:,:,:) + &
                          k_fc_a
    theta_fc_sl(:,:,:)  = k_fc_at + (1.0_wp + k_fc_bt) * hlp1(:,:,:) + k_fc_ct * hlp1(:,:,:) ** 2._wp

    !>   4.3 fractional water holding capacity at saturation
    !>
    hlp1(:,:,:)         = k_sat_s * sand_spq_init_sl(:,:,:) + &
                          k_sat_c * clay_spq_init_sl(:,:,:) + &
                          k_sat_sc * sand_spq_init_sl(:,:,:) * clay_spq_init_sl(:,:,:) + &
                          k_sat_a
    theta_sat_sl(:,:,:) = theta_fc_sl(:,:,:) + k_sat_at + (1.0_wp + k_sat_bt) * hlp1(:,:,:) - &
                          k_sat_ct * sand_spq_init_sl(:,:,:) + k_sat_dt

    !>   4.4 parameters describing the moisture-tension curve
    !>   note that for convenience, B has a negative sign and A is in MPa
    saxtonB(:,:,:)      = - (LOG(wpot_pwp) - LOG(wpot_fc)) / (LOG(theta_fc_sl(:,:,:)) - LOG(theta_pwp_sl(:,:,:)))
    saxtonA(:,:,:)      = - (EXP( LOG(wpot_fc) - saxtonB(:,:,:) * LOG(theta_fc_sl(:,:,:)))) / 1000.0_wp

    !>   4.5 parameters describing the moisture-diffusivity curve
    !>
    kdiff_sat_sl(:,:,:) = kdiff_sat_max * (theta_sat_sl(:,:,:) - theta_fc_sl(:,:,:)) ** (3.0_wp + 1._wp / saxtonB(:,:,:))
    saxtonC(:,:,:)      = 3.0_wp - 2._wp * saxtonB(:,:,:)

    !> 5.0 calc depth of roots
    !>   given by prescribed AWC (soil water in the rooting zone [m])
    !>
    !>   the init of 'root_fraction_sl(:,:,:)' (in sec 6.3) depends on this init of the root_depth(:)
    !>
    soil_awc(:,:)         = 0._wp
    w_soil_root_pwp(:,:)  = 0._wp
    w_soil_root_fc(:,:)   = 0._wp
    root_depth(:,:)       = 0._wp
    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1,INT(num_sl_above_bedrock(ic, iblk))
          hlp1(ic, is, iblk) = soil_lay_width_sl(ic, is, iblk) &
            &                     * (theta_fc_sl(ic, is, iblk) - theta_pwp_sl(ic, is, iblk))
          soil_awc(ic, iblk)    = soil_awc(ic, iblk) + hlp1(ic, is, iblk)
          IF (soil_awc(ic, iblk) < (spq_init_vars%soil_sat_water_content(ic, iblk) / 1000._wp)) THEN
            root_depth(ic, iblk)      = root_depth(ic, iblk) + soil_lay_width_sl(ic, is, iblk)
            w_soil_root_pwp(ic, iblk) = w_soil_root_pwp(ic, iblk) + theta_pwp_sl(ic, is, iblk) &
              &                        * soil_lay_width_sl(ic, is, iblk)
            w_soil_root_fc(ic, iblk)  = w_soil_root_fc(ic, iblk) + theta_fc_sl(ic, is, iblk) &
              &                        * soil_lay_width_sl(ic, is, iblk)
          ELSE
            hlp2(ic, iblk)            = MAX(0._wp, spq_init_vars%soil_sat_water_content(ic, iblk) &
              &                        / 1000._wp - soil_awc(ic, iblk) + hlp1(ic, is, iblk))
            root_depth(ic, iblk)      = root_depth(ic, iblk) + hlp2(ic, iblk) &
              &                        / (theta_fc_sl(ic, is, iblk) - theta_pwp_sl(ic, is, iblk))
            w_soil_root_pwp(ic, iblk) = w_soil_root_pwp(ic, iblk) + theta_pwp_sl(ic, is, iblk) * hlp2(ic, iblk) &
              &                        / (theta_fc_sl(ic, is, iblk) - theta_pwp_sl(ic, is, iblk))
            w_soil_root_fc(ic, iblk)  = w_soil_root_fc(ic, iblk) + theta_fc_sl(ic, is, iblk) * hlp2(ic, iblk) &
              &                        / (theta_fc_sl(ic, is, iblk) - theta_pwp_sl(ic, is, iblk))
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    !>   5.1 water in the surface reservoir [m]
    !>
    w_skin(:,:) = 0.0_wp

    !>   5.2 adjust initial water amount to prescribed theta
    !>
    w_soil_root_theta(:,:)  = dsl4jsb_Config(SPQ_)%soil_theta_prescribe
    w_soil_root(:,:)        = w_soil_root_pwp(:,:) &
      &                     + dsl4jsb_Config(SPQ_)%soil_theta_prescribe * (w_soil_root_fc(:,:) - w_soil_root_pwp(:,:))
    w_soil_root_pot(:,:)    = saxtonA(:, 1, :) * (MIN(w_soil_root_fc,w_soil_root(:,:)) / root_depth(:,:)) ** saxtonB(:, 1, :)

    !>   5.3 root fraction
    !>    root_fraction_sl(:,:,:) is a variable of veg_mem and depends on the init of root_depth(:,:)
    !> ... moved to veg_init!
    !>

    !> 6.0 soil water and ice pwp, fc and sat per layer depth
    !>
    w_soil_pwp_sl(:,:,:)  = 0.0_wp
    w_soil_fc_sl(:,:,:)   = 0.0_wp
    w_soil_sat_sl(:,:,:)  = 0.0_wp
    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1,INT(num_sl_above_bedrock(ic, iblk))
          IF (soil_depth_sl(ic, is, iblk) > 0.0_wp) THEN
            w_soil_pwp_sl(ic, is, iblk) = theta_pwp_sl(ic, is, iblk) * soil_depth_sl(ic, is, iblk)
            w_soil_fc_sl(ic, is, iblk)  = theta_fc_sl(ic, is, iblk) * soil_depth_sl(ic, is, iblk)
            w_soil_sat_sl(ic, is, iblk) = theta_sat_sl(ic, is, iblk) * soil_depth_sl(ic, is, iblk)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    !>   6.1 adjust initial water amount to prescribed theta
    !>   avoid water stress (i.e., dry soil) of plants at the first timestep
    !>
    w_soil_sl(:,:,:) = w_soil_pwp_sl(:,:,:) &
      &                + (dsl4jsb_Config(SPQ_)%soil_theta_prescribe * (w_soil_fc_sl(:,:,:) - w_soil_pwp_sl(:,:,:)))

    !>   6.2 soil water potential per layer
    !>
    WHERE(soil_depth_sl(:,:,:) > eps8)
      hlp1(:,:,:)          = MIN(w_soil_fc_sl(:,:,:), w_soil_sl(:,:,:)) / soil_depth_sl(:,:,:)
      w_soil_pot_sl(:,:,:) = saxtonA(:,:,:) * hlp1(:,:,:) ** saxtonB(:,:,:)
    ELSEWHERE
      w_soil_pot_sl(:,:,:) = 0.0_wp
    ENDWHERE

    ! Deallocate local allocatable variables
    DEALLOCATE(soil_awc)
    DEALLOCATE(hlp2)
    DEALLOCATE(hlp1)
    DEALLOCATE(theta_pwp_sl)
    DEALLOCATE(theta_fc_sl)
    DEALLOCATE(theta_sat_sl)
    DEALLOCATE(sand_spq_init_sl)
    DEALLOCATE(clay_spq_init_sl)
    DEALLOCATE(silt_spq_init_sl)

  END SUBROUTINE spq_init_ic_bc

  ! ======================================================================================================= !
  !> de-allocate SPQ_ init vars
  !>
  SUBROUTINE spq_finalize_init_vars
    DEALLOCATE( &
      & spq_init_vars%soil_sand, &
      & spq_init_vars%soil_silt, &
      & spq_init_vars%soil_clay, &
      & spq_init_vars%soil_sat_water_content, &
      & spq_init_vars%soil_bulk_density)
  END SUBROUTINE spq_finalize_init_vars

#ifdef __QUINCY_STANDALONE__
  SUBROUTINE spq_qs_read_init_vars(tile)
    USE mo_jsb_model_class,    ONLY: t_jsb_model
    USE mo_jsb_class,          ONLY: get_model
    USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
    USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
    USE mo_jsb_io,             ONLY: missval

    dsl4jsb_Use_config(SPQ_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(SPQ_)

    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_soil_sb
    INTEGER                     :: nproma
    INTEGER                     :: nblks
    INTEGER                     :: nsoil_sb
    CHARACTER(len=*), PARAMETER :: routine = modname//':spq_qs_read_init_vars'

    IF (debug_on()) CALL message(routine, 'Reading/setting SPQ_ init vars')

    model  => Get_model(tile%owner_model_id)
    hgrid  => Get_grid(model%grid_id)
    nproma = hgrid%Get_nproma()
    nblks  = hgrid%Get_nblks()
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    nsoil_sb      =  vgrid_soil_sb%n_levels

    ALLOCATE( &
      & spq_init_vars%soil_sand               (nproma, nblks),       &
      & spq_init_vars%soil_silt               (nproma, nblks),       &
      & spq_init_vars%soil_clay               (nproma, nblks),       &
      & spq_init_vars%soil_sat_water_content  (nproma, nblks),       &
      & spq_init_vars%soil_bulk_density       (nproma, nblks)        &
      & )

    spq_init_vars%soil_sand              (:,:) = missval
    spq_init_vars%soil_silt              (:,:) = missval
    spq_init_vars%soil_clay              (:,:) = missval
    spq_init_vars%soil_sat_water_content (:,:) = missval
    spq_init_vars%soil_bulk_density      (:,:) = missval

    dsl4jsb_Get_config(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !

    !> sand
    !>
    spq_init_vars%soil_sand(:,:) = dsl4jsb_Config(SPQ_)%soil_sand

    !> silt
    !>
    spq_init_vars%soil_silt(:,:) = dsl4jsb_Config(SPQ_)%soil_silt

    !> clay
    !>
    spq_init_vars%soil_clay(:,:) = dsl4jsb_Config(SPQ_)%soil_clay

    !> bulk density
    !>
    spq_init_vars%soil_bulk_density(:,:) = dsl4jsb_Config(SPQ_)%bulk_density

    !> saturated water content (awc)
    !>
    spq_init_vars%soil_sat_water_content(:,:) = dsl4jsb_Config(SPQ_)%soil_awc_prescribe
  END SUBROUTINE spq_qs_read_init_vars
#else
  ! ======================================================================================================= !
  !> read SPQ_ init vars from input file
  !>
  SUBROUTINE spq_read_init_vars(tile)
    USE mo_jsb_model_class,    ONLY: t_jsb_model
    USE mo_jsb_class,          ONLY: get_model
    USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
    USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
    USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_jsb_io,             ONLY: missval

    dsl4jsb_Use_config(SPQ_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(SPQ_)

    REAL(wp), POINTER           :: ptr_2D(:, :)
    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_soil_sb
    TYPE(t_input_file)          :: input_file
    INTEGER                     :: nproma
    INTEGER                     :: nblks
    INTEGER                     :: nsoil_sb
    CHARACTER(len=*), PARAMETER :: routine = modname//':spq_read_init_vars'

    IF (debug_on()) CALL message(routine, 'Reading/setting SPQ_ init vars')

    model  => Get_model(tile%owner_model_id)
    hgrid  => Get_grid(model%grid_id)
    nproma = hgrid%Get_nproma()
    nblks  = hgrid%Get_nblks()
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    nsoil_sb      =  vgrid_soil_sb%n_levels

    ALLOCATE( &
      & spq_init_vars%soil_sand               (nproma, nblks),       &
      & spq_init_vars%soil_silt               (nproma, nblks),       &
      & spq_init_vars%soil_clay               (nproma, nblks),       &
      & spq_init_vars%soil_sat_water_content  (nproma, nblks),       &
      & spq_init_vars%soil_bulk_density       (nproma, nblks)        &
      & )

    spq_init_vars%soil_sand              (:,:) = missval
    spq_init_vars%soil_silt              (:,:) = missval
    spq_init_vars%soil_clay              (:,:) = missval
    spq_init_vars%soil_sat_water_content (:,:) = missval
    spq_init_vars%soil_bulk_density      (:,:) = missval

    dsl4jsb_Get_config(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SPQ_)%bc_quincy_soil_filename), model%grid_id)


    !> sand
    !>
    ptr_2D => input_file%Read_2d(  &
      & variable_name='SNDPPT',    &
      & fill_array = spq_init_vars%soil_sand)
    ! convert from % to [0,1]
    spq_init_vars%soil_sand(:,:) = spq_init_vars%soil_sand(:,:) / 100.0_wp
    spq_init_vars%soil_sand(:,:) = MIN(1.0_wp, MERGE(ptr_2D, 0.3_wp, ptr_2D > 0.0_wp))

    !> silt
    !>
    ptr_2D => input_file%Read_2d(  &
      & variable_name='SLTPPT',    &
      & fill_array = spq_init_vars%soil_silt)
    ! convert from % to [0,1]
    spq_init_vars%soil_silt(:,:) = spq_init_vars%soil_silt(:,:) / 100.0_wp
    spq_init_vars%soil_silt(:,:) = MIN(1.0_wp, MERGE(ptr_2D, 0.4_wp, ptr_2D > 0.0_wp))

    !> clay
    !>
    spq_init_vars%soil_clay(:,:) = MIN(1.0_wp, 1.0_wp - spq_init_vars%soil_sand(:,:) + spq_init_vars%soil_silt(:,:))

    !> bulk density
    !>
    ptr_2D => input_file%Read_2d(  &
      & variable_name='BLDFIE',    &
      & fill_array = spq_init_vars%soil_bulk_density)
    spq_init_vars%soil_bulk_density(:,:) = MERGE(ptr_2D, 1000.0_wp, ptr_2D > 0.0_wp)

    !> saturated water content (awc)
    !>
    ! values in bc file are negative, but need to be positive
    ptr_2D => input_file%Read_2d(  &
      & variable_name='AWCtS',     &
      & fill_array = spq_init_vars%soil_sat_water_content)
    ! multiply with -1 to get positive values
    spq_init_vars%soil_sat_water_content(:,:) = -1.0_wp * spq_init_vars%soil_sat_water_content(:,:)
    ! replace values larger than 2000 with 300
    ! that would cover missing values (e.g. grid cells over ocean with value 9999, after multiplication with -1)
    spq_init_vars%soil_sat_water_content(:,:) = MERGE(ptr_2D, 300.0_wp, ptr_2D < 2000.0_wp)
    ! replace values larger than 800 with 800
    spq_init_vars%soil_sat_water_content(:,:) = MERGE(ptr_2D, 800.0_wp, ptr_2D < 800.0_wp)
    ! replace values less than 10 with 10
    spq_init_vars%soil_sat_water_content(:,:) = MERGE(ptr_2D, 10.0_wp, ptr_2D > 10.0_wp)

    CALL input_file%Close()
  END SUBROUTINE spq_read_init_vars
#endif

  ! ======================================================================================================= !
  !>ensure a certain minimum proportion for sand, clay and silt
  !>
  !> Saxton equations with too small soil-texture fractions may lead to unrealistic percolation rates and soil water potentials
  !>
  !> using these modified soil texture for "saxton functions in spq_init_ic_bc" only
  !>
  SUBROUTINE calc_soil_texture_saxton_compatible( &
    & nproma, &
    & nsoil_sb, &
    & sand_sl, &
    & silt_sl, &
    & clay_sl )

    INTEGER,                                INTENT(in)    :: nproma        !< dimensions
    INTEGER,                                INTENT(in)    :: nsoil_sb      !< number of soil layers
    REAL(wp), DIMENSION(nproma, nsoil_sb),  INTENT(inout) :: sand_sl       !< sand proportion per layer
    REAL(wp), DIMENSION(nproma, nsoil_sb),  INTENT(inout) :: silt_sl       !< silt proportion per layer
    REAL(wp), DIMENSION(nproma, nsoil_sb),  INTENT(inout) :: clay_sl       !< clay proportion per layer
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), DIMENSION(nproma, nsoil_sb) :: hlp_sand_sl, &
                                             hlp_silt_sl, &
                                             hlp_clay_sl
    REAL(wp), DIMENSION(nproma, nsoil_sb) :: sand_original_sl, &
                                             silt_original_sl, &
                                             clay_original_sl
    REAL(wp)                              :: min_prop_soil_texture, &
                                             max_prop_soil_texture
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_soil_texture_saxton_compatible'
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.9 init local variables
    !>
    min_prop_soil_texture = 0.2_wp
    max_prop_soil_texture = 1.0_wp - (2.0_wp * min_prop_soil_texture)
    hlp_sand_sl(:,:)      = 0.0_wp
    hlp_clay_sl(:,:)      = 0.0_wp
    hlp_silt_sl(:,:)      = 0.0_wp
    sand_original_sl(:,:) = sand_sl(:,:)
    clay_original_sl(:,:) = clay_sl(:,:)
    silt_original_sl(:,:) = silt_sl(:,:)

    !>1.0 ensure a minimum proportion 'min_prop_soil_texture' for sand, clay and silt
    !>  this is to ensure applying the empirical functions by Saxton within the range of the texture they were designed for
    !>
    WHERE (sand_sl(:,:) < min_prop_soil_texture .AND. clay_sl(:,:) < min_prop_soil_texture)
      sand_sl(:,:)     = min_prop_soil_texture
      clay_sl(:,:)     = min_prop_soil_texture
      silt_sl(:,:)     = max_prop_soil_texture
    ELSEWHERE (sand_sl(:,:) < min_prop_soil_texture .AND. silt_sl(:,:) < min_prop_soil_texture)
      sand_sl(:,:)     = min_prop_soil_texture
      silt_sl(:,:)     = min_prop_soil_texture
      clay_sl(:,:)     = max_prop_soil_texture
    ELSEWHERE (clay_sl(:,:) < min_prop_soil_texture .AND. silt_sl(:,:) < min_prop_soil_texture)
      clay_sl(:,:)     = min_prop_soil_texture
      silt_sl(:,:)     = min_prop_soil_texture
      sand_sl(:,:)     = max_prop_soil_texture
    ELSEWHERE (sand_sl(:,:) < min_prop_soil_texture)
      hlp_sand_sl(:,:) = min_prop_soil_texture - sand_sl(:,:)
      ! proportion of hlp_sand_sl(:,:) to be drawn from clay and silt
      hlp_clay_sl(:,:) = clay_sl(:,:) + (1.0_wp - (clay_sl(:,:) + silt_sl(:,:))) / 2.0_wp
      hlp_silt_sl(:,:) = silt_sl(:,:) + (1.0_wp - (clay_sl(:,:) + silt_sl(:,:))) / 2.0_wp
      ! where new clay proportion in soil would be less than "min_prop_soil_texture + 0.01" use only from silt
      WHERE     ( (clay_sl(:,:) - hlp_sand_sl(:,:) * hlp_clay_sl(:,:)) < (min_prop_soil_texture + 0.01) )
        silt_sl(:,:) = silt_sl(:,:) - hlp_sand_sl(:,:)
      ! where new silt proportion in soil would be less than "min_prop_soil_texture + 0.01" use only from clay
      ELSEWHERE ( (silt_sl(:,:) - hlp_sand_sl(:,:) * hlp_silt_sl(:,:)) < (min_prop_soil_texture + 0.01) )
        clay_sl(:,:) = clay_sl(:,:) - hlp_sand_sl(:,:)
      ELSEWHERE
        clay_sl(:,:) = clay_sl(:,:) - hlp_sand_sl(:,:) * hlp_clay_sl(:,:)
        silt_sl(:,:) = silt_sl(:,:) - hlp_sand_sl(:,:) * hlp_silt_sl(:,:)
      ENDWHERE
      sand_sl(:,:) = min_prop_soil_texture
    ELSEWHERE (clay_sl(:,:) < min_prop_soil_texture)
      hlp_clay_sl(:,:) = min_prop_soil_texture - clay_sl(:,:)
      ! proportion of hlp_clay_sl(:,:) to be drawn from sand and silt
      hlp_sand_sl(:,:) = sand_sl(:,:) + (1.0_wp - (sand_sl(:,:) + silt_sl(:,:))) / 2.0_wp
      hlp_silt_sl(:,:) = silt_sl(:,:) + (1.0_wp - (sand_sl(:,:) + silt_sl(:,:))) / 2.0_wp
      ! where new sand proportion in soil would be less than "min_prop_soil_texture + 0.01" use only from silt
      WHERE     ( (sand_sl(:,:) - hlp_clay_sl(:,:) * hlp_sand_sl(:,:)) < (min_prop_soil_texture + 0.01) )
        silt_sl(:,:) = silt_sl(:,:) - hlp_clay_sl(:,:)
      ! where new silt proportion in soil would be less than "min_prop_soil_texture + 0.01" use only from sand
      ELSEWHERE ( (silt_sl(:,:) - hlp_clay_sl(:,:) * hlp_silt_sl(:,:)) < (min_prop_soil_texture + 0.01) )
        sand_sl(:,:) = sand_sl(:,:) - hlp_clay_sl(:,:)
      ELSEWHERE
        sand_sl(:,:) = sand_sl(:,:) - hlp_clay_sl(:,:) * hlp_sand_sl(:,:)
        silt_sl(:,:) = silt_sl(:,:) - hlp_clay_sl(:,:) * hlp_silt_sl(:,:)
      ENDWHERE
      clay_sl(:,:) = min_prop_soil_texture
    ELSEWHERE (silt_sl(:,:) < min_prop_soil_texture)
      hlp_silt_sl(:,:) = min_prop_soil_texture - silt_sl(:,:)
      ! proportion of hlp_silt_sl(:,:) to be drawn from sand and clay
      hlp_sand_sl(:,:) = sand_sl(:,:) + (1.0_wp - (sand_sl(:,:) + clay_sl(:,:))) / 2.0_wp
      hlp_clay_sl(:,:) = clay_sl(:,:) + (1.0_wp - (sand_sl(:,:) + clay_sl(:,:))) / 2.0_wp
      ! where new sand proportion in soil would be less than "min_prop_soil_texture + 0.01" use only from clay
      WHERE     ( (sand_sl(:,:) - hlp_silt_sl(:,:) * hlp_sand_sl(:,:)) < (min_prop_soil_texture + 0.01) )
        clay_sl(:,:) = clay_sl(:,:) - hlp_silt_sl(:,:)
      ! where new clay proportion in soil would be less than "min_prop_soil_texture + 0.01" use only from sand
      ELSEWHERE ( (clay_sl(:,:) - hlp_silt_sl(:,:) * hlp_clay_sl(:,:)) < (min_prop_soil_texture + 0.01) )
        sand_sl(:,:) = sand_sl(:,:) - hlp_silt_sl(:,:)
      ELSEWHERE
        sand_sl(:,:) = sand_sl(:,:) - hlp_silt_sl(:,:) * hlp_sand_sl(:,:)
        clay_sl(:,:) = clay_sl(:,:) - hlp_silt_sl(:,:) * hlp_clay_sl(:,:)
      ENDWHERE
      silt_sl(:,:) = min_prop_soil_texture
    ENDWHERE

    !>  1.1 ensure all adds up to exactly 1.0_wp
    !>
    clay_sl(:,:) = 1.0_wp
    clay_sl(:,:) = clay_sl(:,:) - (sand_sl(:,:) + silt_sl(:,:))

    !>  1.2 output any changes in soil texture to stdout or log
    !>      only for QS QUINCY standalone
    !>
#ifdef __QUINCY_STANDALONE__
    IF (ABS(SUM(sand_sl(:,:) - sand_original_sl(:,:))) > eps8) THEN
      CALL message(TRIM(routine), ' Modified the soil texture sand input for saxton calculation.')
      print*,TRIM(routine), " Modified the soil texture sand input for saxton calculation: ", sand_sl(:,:) - sand_original_sl(:,:)
    ENDIF
    IF (ABS(SUM(clay_sl(:,:) - clay_original_sl(:,:))) > eps8) THEN
      CALL message(TRIM(routine), ' Modified the soil texture clay input for saxton calculation.')
      print*,TRIM(routine), " Modified the soil texture clay input for saxton calculation: ", clay_sl(:,:) - clay_original_sl(:,:)
    ENDIF
    IF (ABS(SUM(silt_sl(:,:) - silt_original_sl(:,:))) > eps8) THEN
      CALL message(TRIM(routine), ' Modified the soil texture silt input for saxton calculation.')
      print*,TRIM(routine), " Modified the soil texture silt input for saxton calculation: ", silt_sl(:,:) - silt_original_sl(:,:)
    ENDIF
#endif
  END SUBROUTINE calc_soil_texture_saxton_compatible

#endif
END MODULE mo_spq_init
