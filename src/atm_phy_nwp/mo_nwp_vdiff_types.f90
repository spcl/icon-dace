!
! Types for the interface to the VDIFF turbulence scheme and JSBACH land-surface scheme.
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

MODULE mo_nwp_vdiff_types

  USE mo_cdi, ONLY: &
      & DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, GRID_UNSTRUCTURED, TSTEP_CONSTANT, &
      & TSTEP_INSTANT
  USE mo_cdi_constants, ONLY: GRID_CELL, GRID_UNSTRUCTURED_CELL
  USE mo_cf_convention, ONLY: t_cf_var
  USE mo_coupling_config, ONLY: is_coupled_to_ocean
  USE mo_grib2, ONLY: t_grib2_var, grib2_var
  USE mo_io_config, ONLY: lnetcdf_flt64_output
  USE mo_kind, ONLY: wp
  USE mo_var_groups, ONLY: groups
  USE mo_var_list, ONLY: t_var_list_ptr, add_var, add_ref
  USE mo_zaxis_type, ONLY: ZA_REFERENCE, ZA_SURFACE

  USE mtime, ONLY: datetime

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_nwp_vdiff_albedos
  PUBLIC :: t_nwp_vdiff_sea_state
  PUBLIC :: t_nwp_vdiff_state

  PUBLIC :: SFC_LAND
  PUBLIC :: SFC_WATER
  PUBLIC :: SFC_ICE
  PUBLIC :: SFC_NUM
  PUBLIC :: SFC_NAMES

  PUBLIC :: SFT_LAND
  PUBLIC :: SFT_LWTR
  PUBLIC :: SFT_LICE
  PUBLIC :: SFT_SWTR
  PUBLIC :: SFT_SICE
  PUBLIC :: SFT_NUM
  PUBLIC :: SFT_NAMES
  PUBLIC :: SFT_CLASS

  ! VDIFF surface classes.
  INTEGER, PARAMETER :: SFC_LAND = 1 !< VDIFF Index of land surface.
  INTEGER, PARAMETER :: SFC_WATER = 2 !< VDIFF Index of open water.
  INTEGER, PARAMETER :: SFC_ICE = 3 !< VDIFF Index of ice over water.
  INTEGER, PARAMETER :: SFC_NUM = 3 !< VDIFF Number of surface types.

  !> Surface class names.
  CHARACTER(len=*), PARAMETER :: SFC_NAMES(SFC_NUM) = [ 'land ', 'water', 'ice  ' ]

  ! Internal surface types.
  INTEGER, PARAMETER :: SFT_LAND = 1 !< Land surface type.
  INTEGER, PARAMETER :: SFT_LWTR = 2 !< Lake water surface type.
  INTEGER, PARAMETER :: SFT_LICE = 3 !< Lake ice surface type.
  INTEGER, PARAMETER :: SFT_SWTR = 4 !< Sea water surface type.
  INTEGER, PARAMETER :: SFT_SICE = 5 !< Sea ice surface type.
  INTEGER, PARAMETER :: SFT_NUM = 5 !< Number of surface types.

  !> Surface type names.
  CHARACTER(len=*), PARAMETER :: SFT_NAMES(SFT_NUM) = [ 'land', 'lwtr', 'lice', 'swtr', 'sice' ]
  !> Surface type classes.
  INTEGER, PARAMETER :: SFT_CLASS(SFT_NUM) = [ SFC_LAND, SFC_WATER, SFC_ICE, SFC_WATER, SFC_ICE ]

  !> Albedos of different surfaces per surface and radiation type.
  TYPE t_nwp_vdiff_albedos
    !> Albedo for near-IR diffuse radiation [1] (nproma,nblks_c,SFT_NUM).
    REAL(wp), CONTIGUOUS, POINTER :: alb_nir_dif(:,:,:) => NULL()
    !> Albedo for near-IR direct radiation [1] (nproma,nblks_c,SFT_NUM).
    REAL(wp), CONTIGUOUS, POINTER :: alb_nir_dir(:,:,:) => NULL()
    !> Albedo for visible diffuse radiation [1] (nproma,nblks_c,SFT_NUM).
    REAL(wp), CONTIGUOUS, POINTER :: alb_vis_dif(:,:,:) => NULL()
    !> Albedo for visible direct radiation [1] (nproma,nblks_c,SFT_NUM).
    REAL(wp), CONTIGUOUS, POINTER :: alb_vis_dir(:,:,:) => NULL()

    !> Long-wave emissivity [1] (nproma,nblks_c,SFT_NUM).
    REAL(wp), CONTIGUOUS, POINTER :: lw_emissivity(:,:,:) => NULL()

    !> Deallocate fields when structure is destroyed.
    LOGICAL :: deallocate = .FALSE.
  CONTAINS
    PROCEDURE :: init => nwp_vdiff_albedos_init

    PROCEDURE :: d2h => nwp_vdiff_albedos_d2h
    PROCEDURE :: h2d => nwp_vdiff_albedos_h2d

    PROCEDURE, PRIVATE :: nwp_vdiff_albedos_copy
    GENERIC :: ASSIGNMENT(=) => nwp_vdiff_albedos_copy

    FINAL :: nwp_vdiff_albedos_destroy
  END TYPE t_nwp_vdiff_albedos


  TYPE t_nwp_vdiff_sea_state

    !> Surface albedos.
    !> Albedo for near-IR diffuse radiation [1] (nproma,nblks_c,SFT_SWTR:SFT_SICE).
    REAL(wp), CONTIGUOUS, POINTER :: alb_nir_dif(:,:,:) => NULL()
    !> Albedo for near-IR direct radiation [1] (nproma,nblks_c,SFT_SWTR:SFT_SICE).
    REAL(wp), CONTIGUOUS, POINTER :: alb_nir_dir(:,:,:) => NULL()
    !> Albedo for visible diffuse radiation [1] (nproma,nblks_c,SFT_SWTR:SFT_SICE).
    REAL(wp), CONTIGUOUS, POINTER :: alb_vis_dif(:,:,:) => NULL()
    !> Albedo for visible direct radiation [1] (nproma,nblks_c,SFT_SWTR:SFT_SICE).
    REAL(wp), CONTIGUOUS, POINTER :: alb_vis_dir(:,:,:) => NULL()

    !> Long-wave emissivity [1] (nproma,nblks_c,SFT_SWTR:SFT_SICE).
    REAL(wp), CONTIGUOUS, POINTER :: lw_emissivity(:,:,:) => NULL()

    !> Zonal ocean surface velocity [m/s] (nproma,nblks_c).
    REAL(wp), CONTIGUOUS, POINTER :: ocean_u(:,:) => NULL()
    !> Meridional ocean surface velocity [m/s] (nproma,nblks_c).
    REAL(wp), CONTIGUOUS, POINTER :: ocean_v(:,:) => NULL()

    !> Natural CO2 flux over sea [kg/m**2(ocean)/s] (nproma,nblks_c).
    REAL(wp), CONTIGUOUS, POINTER :: flx_co2_natural_sea(:,:) => NULL()

    !> Difference between initial sea surface temperature and the climatological field at
    !! reference time.
    REAL(wp), CONTIGUOUS, POINTER :: t_seasfc_offset(:,:) => NULL()

    !> Reference time for initial sea surface temperature.
    TYPE(datetime) :: time_ref_t_seasfc

    !> Time for last sea surface temperature update.
    TYPE(datetime) :: time_last_update_t_seasfc

  CONTAINS

    PROCEDURE :: init => nwp_vdiff_sea_state_init

    PROCEDURE :: d2h => nwp_vdiff_sea_state_d2h
    PROCEDURE :: h2d => nwp_vdiff_sea_state_h2d

    FINAL :: nwp_vdiff_sea_state_destroy

  END TYPE t_nwp_vdiff_sea_state


  !>
  !! Memory structure for #nwp_vdiff. Contains fields that need to be saved between time steps.
  !! Fields get propagated in time by nwp_vdiff.
  TYPE t_nwp_vdiff_state
    !> Exchange coefficient for heat on half-levels above the surface [m**2/s]
    !! (nproma,nlev,nblks_c).
    REAL(wp), CONTIGUOUS, POINTER :: exchange_coeff_h(:,:,:) => NULL()
    !> Exchange coefficient for momentum on half-levels above the surface [m**2/s]
    !! (nproma,nlev,nblks_c).
    REAL(wp), CONTIGUOUS, POINTER :: exchange_coeff_m(:,:,:) => NULL()
    !> Surface exchange coefficient for heat over surface classes [m/s]
    !! (nproma,nblks_c,SFC_NUM).
    REAL(wp), CONTIGUOUS, POINTER :: exchange_coeff_h_sfc(:,:,:) => NULL()
    !> Surface exchange coefficient for momentum over surface classes [m/s]
    !! (nproma,nblks_c,SFC_NUM).
    REAL(wp), CONTIGUOUS, POINTER :: exchange_coeff_m_sfc(:,:,:) => NULL()

    !> Lake-ice fraction from jsbach [1] (nproma,nblks_c).
    REAL(wp), CONTIGUOUS, POINTER :: fr_ice_on_lake(:,:) => NULL()

    !> Temperature of each surface type [K] (nproma,nblks_c,SFT_NUM).
    REAL(wp), CONTIGUOUS, POINTER :: temp_sft(:,:,:) => NULL()

    !> Roughness length for momentum over surface classes [m] (nproma,nblks_c,SFC_NUM).
    REAL(wp), CONTIGUOUS, POINTER :: z0m_sfc(:,:,:) => NULL()
    !> Roughness length for heat over land [m] (nproma,nblks_c).
    REAL(wp), CONTIGUOUS, POINTER :: z0h_land(:,:) => NULL()

    !> Natural CO2 flux over land [kg/m**2/s] (nproma,nblks_c).
    REAL(wp), CONTIGUOUS, POINTER :: flx_co2_natural_land(:,:) => NULL()

    !> Total turbulent energy [J/kg?] (nproma,nlev,nblks_c)
    REAL(wp), CONTIGUOUS, POINTER :: total_turbulence_energy(:,:,:) => NULL()

    !> Air coefficient for moisture flux (C_air in jsbach (2.56)) [1] (nproma,nblks_c).
    REAL(wp), CONTIGUOUS, POINTER :: fact_q_air(:,:) => NULL()
    !> Surface coefficient for moisture flux (C_sat in jsbach (2.56)) [1] (nproma,nblks_c).
    REAL(wp), CONTIGUOUS, POINTER :: fact_qsat_srf(:,:) => NULL()

    !> Grid-box mean friction velocity at previous time step [m/s] (nproma,nblks_c).
    REAL(wp), CONTIGUOUS, POINTER :: ustar(:,:) => NULL()

    !> Convective velocity scale per surface class [m/s] (nproma,nblks_c,SFC_NUM).
    REAL(wp), CONTIGUOUS, POINTER :: wstar_sfc(:,:,:) => NULL()

    !> Variance of virtual potential temperature [K**2?] (nproma,nlev,nblks_c).
    REAL(wp), CONTIGUOUS, POINTER :: theta_v_var(:,:,:) => NULL()
    !> Total water variance [(kg/kg)**2] (nproma,nlev,nblks_c).
    REAL(wp), CONTIGUOUS, POINTER :: total_water_var(:,:,:) => NULL()

    !> Latent heat flux over surface type [W/m**2] (nproma,nblks_c,SFT_NUM).
    REAL(wp), CONTIGUOUS, POINTER :: flx_heat_latent_sft(:,:,:) => NULL()
    !> Sensible heat flux over surface type [W/m**2] (nproma,nblks_c,SFT_NUM).
    REAL(wp), CONTIGUOUS, POINTER :: flx_heat_sensible_sft(:,:,:) => NULL()

    !> Sea state.
    TYPE(t_nwp_vdiff_sea_state) :: sea_state

  CONTAINS
    PROCEDURE :: init => nwp_vdiff_state_init

    PROCEDURE :: d2h => nwp_vdiff_state_d2h
    PROCEDURE :: h2d => nwp_vdiff_state_h2d

    FINAL :: nwp_vdiff_state_destroy
  END TYPE t_nwp_vdiff_state

CONTAINS

  !> Initialize a `t_nwp_vdiff_state` structure, adding the variables to the given varlist.
  SUBROUTINE nwp_vdiff_state_init (self, nproma, nlev, nblks_c, varlist)

    CLASS(t_nwp_vdiff_state), INTENT(OUT) :: self !< Object to initialize.
    INTEGER, INTENT(IN) :: nproma !< Block size.
    INTEGER, INTENT(IN) :: nlev !< Number of levels.
    INTEGER, INTENT(IN) :: nblks_c !< Number of cell blocks.
    TYPE(t_var_list_ptr), INTENT(INOUT) :: varlist !< Varlist to append variables.

    TYPE(t_cf_var) :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: datatype_flt
    INTEGER :: grib2_bits

    INTEGER :: shape2d(2) !< Shape of 2D fields.
    INTEGER :: shape3d(3) !< Shape of 3D fields.
    INTEGER :: shape2d_sft(3) !< Shape of 2D fields with surface-type index.
    INTEGER :: shape2d_sfc(3) !< Shape of 2D fields with surface-class index.

    INTEGER :: i

    REAL(wp), POINTER :: p2d(:,:)

    shape2d = [ nproma, nblks_c ]
    shape3d = [ nproma, nlev, nblks_c ]
    shape2d_sft = [ nproma, nblks_c, SFT_NUM ]
    shape2d_sfc = [ nproma, nblks_c, SFC_NUM ]

    IF (lnetcdf_flt64_output) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    END IF

    grib2_bits = DATATYPE_PACK16

    !$ACC ENTER DATA ASYNC(1) COPYIN(self)

    ! self%exchange_coeff_h(nproma,nlev,nblks_c)
    cf_desc = t_cf_var('exchange_coeff_h', 'm2 s-1', 'Exchange coefficient for heat', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'exchange_coeff_h', self%exchange_coeff_h, &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
        & ldims=shape3d, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.FALSE., in_group=groups('vdiff') &
      )
    !$ACC ENTER DATA ASYNC(1) ATTACH(self%exchange_coeff_h)

    ! self%exchange_coeff_m(nproma,nlev,nblks_c)
    cf_desc = t_cf_var( &
        & 'exchange_coeff_m', 'm2 s-1', 'Exchange coefficient for momentum', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'exchange_coeff_m', self%exchange_coeff_m, &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
        & ldims=shape3d, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.FALSE., in_group=groups('vdiff') &
      )
    !$ACC ENTER DATA ASYNC(1) ATTACH(self%exchange_coeff_m)

    ! self%exchange_coeff_h_sfc(nproma,nblks_c,SFC_NUM)
    cf_desc = t_cf_var('exchange_coeff_h_sfc', 'm s-1', &
        & 'Surface exchange coefficient for heat over surface classes', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'exchange_coeff_h_sfc', self%exchange_coeff_h_sfc, &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
        & ldims=shape2d_sfc, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.FALSE., loutput=.FALSE., lcontainer=.TRUE. &
      )
    !$ACC ENTER DATA ASYNC(1) ATTACH(self%exchange_coeff_h_sfc)

    DO i = 1, SFC_NUM
      cf_desc = t_cf_var('exchange_coeff_h_' // TRIM(SFC_NAMES(i)), 'm s-1', &
          & 'Surface exchange coefficient for heat over ' // TRIM(SFC_NAMES(i)), datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_ref(varlist, 'exchange_coeff_h_sfc', 'exchange_coeff_h_' // TRIM(SFC_NAMES(i)), &
          & p2d, GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ref_idx=i, &
          & ldims=shape2d, loutput=.TRUE., lrestart=.FALSE., in_group=groups('vdiff-sft') &
        )
    END DO

    ! self%exchange_coeff_m_sfc(nproma,nblks_c,SFC_NUM)
    cf_desc = t_cf_var('exchange_coeff_m_sfc', 'm s-1', &
        & 'Surface exchange coefficient for momentum over surface classes', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'exchange_coeff_m_sfc', self%exchange_coeff_m_sfc, &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
        & ldims=shape2d_sfc, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.FALSE., loutput=.FALSE., lcontainer=.TRUE. &
      )
    !$ACC ENTER DATA ASYNC(1) ATTACH(self%exchange_coeff_m_sfc)

    DO i = 1, SFC_NUM
      cf_desc = t_cf_var('exchange_coeff_m_' // TRIM(SFC_NAMES(i)), 'm s-1', &
          & 'Surface exchange coefficient for momentum over ' // TRIM(SFC_NAMES(i)), datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_ref(varlist, 'exchange_coeff_m_sfc', 'exchange_coeff_m_' // TRIM(SFC_NAMES(i)), &
          & p2d, GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ref_idx=i, &
          & ldims=shape2d, loutput=.TRUE., lrestart=.FALSE., in_group=groups('vdiff-sft') &
        )
    END DO

    ! self%fact_q_air(nproma,nblks_c)
    cf_desc = t_cf_var('fact_q_air', '-', 'Atmosphere factor for humidity flux', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'fact_q_air', self%fact_q_air, &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
        & ldims=shape2d, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.TRUE., in_group=groups('vdiff') &
      )
    !$ACC ENTER DATA ASYNC(1) ATTACH(self%fact_q_air)

    ! self%fact_qsat_srf(nproma,nblks_c)
    cf_desc = t_cf_var('fact_qsat_srf', '-', 'Surface factor for humidity flux', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'fact_qsat_srf', self%fact_qsat_srf, &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
        & ldims=shape2d, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.TRUE., in_group=groups('vdiff') &
      )
    !$ACC ENTER DATA ASYNC(1) ATTACH(self%fact_qsat_srf)

    ! self%flx_co2_natural_land(nproma,nblks_c)
    cf_desc = t_cf_var('flx_co2_natural_land', 'kg m-2 s-1', 'Natural land surface CO2 flux', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'flx_co2_natural_land', self%flx_co2_natural_land, &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
        & ldims=shape2d, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.TRUE., in_group=groups('vdiff') &
      )
    !$ACC ENTER DATA ASYNC(1) ATTACH(self%flx_co2_natural_land)

    ! self%flx_heat_latent_sft(nproma,nblks_c,SFT_NUM)
    cf_desc = t_cf_var('flx_heat_latent_sft', 'W m-2', &
        & 'Latent heat flux over surface types', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'flx_heat_latent_sft', self%flx_heat_latent_sft, &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
        & ldims=shape2d_sft, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.FALSE., loutput=.FALSE., lcontainer=.TRUE. &
      )
    !$ACC ENTER DATA ASYNC(1) ATTACH(self%flx_heat_latent_sft)

    DO i = 1, SFT_NUM
      cf_desc = t_cf_var('flx_heat_latent_' // SFT_NAMES(i), 'W m-2', &
          & 'Latent heat flux over ' // SFT_NAMES(i), datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_ref(varlist, 'flx_heat_latent_sft', 'flx_heat_latent_' // SFT_NAMES(i), p2d, &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ref_idx=i, ldims=shape2d, &
          & loutput=.TRUE., lrestart=.TRUE., in_group=groups('vdiff-sft') &
        )
    END DO

    ! self%flx_heat_sensible_sft(nproma,nblks_c,SFT_NUM)
    cf_desc = t_cf_var('flx_heat_sensible_sft', 'W m-2', &
        & 'Sensible heat flux over surface types', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'flx_heat_sensible_sft', self%flx_heat_sensible_sft, &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
        & ldims=shape2d_sft, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.FALSE., loutput=.FALSE., lcontainer=.TRUE. &
      )
    !$ACC ENTER DATA ASYNC(1) ATTACH(self%flx_heat_sensible_sft)

    DO i = 1, SFT_NUM
      cf_desc = t_cf_var('flx_heat_sensible_' // SFT_NAMES(i), 'W m-2', &
          & 'Sensible heat flux over ' // SFT_NAMES(i), datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_ref(varlist, 'flx_heat_sensible_sft', 'flx_heat_sensible_' // SFT_NAMES(i), p2d, &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ref_idx=i, ldims=shape2d, &
          & loutput=.TRUE., lrestart=.TRUE., in_group=groups('vdiff-sft') &
        )
    END DO

    ! self%fr_ice_on_lake(nproma,nblks_c), &
    cf_desc = t_cf_var('fr_ice_on_lake', 'm2(ice) m-2(lake)', 'Fraction of frozen lake area', &
        & datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'fr_ice_on_lake', self%fr_ice_on_lake, &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
        & ldims=shape2d, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.TRUE., in_group=groups('vdiff') &
      )
    !$ACC ENTER DATA ASYNC(1) ATTACH(self%fr_ice_on_lake)

    ! self%temp_sft(nproma,nblks_c,SFT_NUM), &
    cf_desc = t_cf_var('temp_sft', 'K', 'Surface temperature over surface types', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'temp_sft', self%temp_sft, &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
        & ldims=shape2d_sft, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.FALSE., loutput=.FALSE., lcontainer=.TRUE. &
      )
    !$ACC ENTER DATA ASYNC(1) ATTACH(self%temp_sft)

    DO i = 1, SFT_NUM
      cf_desc = t_cf_var('temp_s_' // SFT_NAMES(i), 'K', &
          & 'Surface temperature over ' // SFT_NAMES(i), datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_ref(varlist, 'temp_sft', 'temp_s_' // SFT_NAMES(i), p2d, &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ref_idx=i, ldims=shape2d, &
          & loutput=.TRUE., lrestart=.TRUE., in_group=groups('vdiff-sft') &
        )
    END DO

    ! self%theta_v_var(nproma,nlev,nblks_c), &
    cf_desc = t_cf_var('theta_v_var', 'K2', 'Virtual potential temperature variance', &
        & datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'theta_v_var', self%theta_v_var, &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
        & ldims=shape3d, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.FALSE. &
        & & ! In principle lrestart=.TRUE., in_group=groups('vdiff'), but currently all-zero.
      )
    !$ACC ENTER DATA ASYNC(1) ATTACH(self%theta_v_var)

    ! self%total_turbulence_energy(nproma,nlev,nblks_c), &
    cf_desc = t_cf_var('total_turbulence_energy', 'J kg-1', 'Total turbulence energy', &
        & datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'total_turbulence_energy', self%total_turbulence_energy, &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
        & ldims=shape3d, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.TRUE., in_group=groups('vdiff') &
      )
    !$ACC ENTER DATA ASYNC(1) ATTACH(self%total_turbulence_energy)

    ! self%total_water_var(nproma,nlev,nblks_c), &
    cf_desc = t_cf_var('total_water_var', 'kg2 kg-2', 'Total water variance', &
        & datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'total_water_var', self%total_water_var, &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
        & ldims=shape3d, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.FALSE. &
        & & ! In principle lrestart=.TRUE., in_group=groups('vdiff'), but currently all-zero.
      )
    !$ACC ENTER DATA ASYNC(1) ATTACH(self%total_water_var)

    ! self%ustar(nproma,nblks_c), &
    cf_desc = t_cf_var('ustar', 'm s-1', 'Friction velocity', &
        & datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'ustar', self%ustar, &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
        & ldims=shape2d, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.TRUE., in_group=groups('vdiff') &
      )
    !$ACC ENTER DATA ASYNC(1) ATTACH(self%ustar)

    ! self%wstar_sfc(nproma,nblks_c,SFC_NUM), &
    cf_desc = t_cf_var('wstar_sfc', 'm s-1', 'Convective velocity scale per surface class', &
        & datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'wstar_sfc', self%wstar_sfc, &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
        & ldims=shape2d_sfc, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.FALSE., loutput=.FALSE., lcontainer=.TRUE. &
      )
    !$ACC ENTER DATA ASYNC(1) ATTACH(self%wstar_sfc)

    DO i = 1, SFC_NUM
      cf_desc = t_cf_var('wstar_' // TRIM(SFC_NAMES(i)), 'm s-1', &
          & 'Convective velocity scale over ' // TRIM(SFC_NAMES(i)), datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_ref(varlist, 'wstar_sfc', 'wstar_' // TRIM(SFC_NAMES(i)), p2d, &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ref_idx=i, ldims=shape2d, &
          & loutput=.TRUE., lrestart=.TRUE., in_group=groups('vdiff-sft') &
        )
    END DO

    ! self%z0h_land(nproma,nblks_c), &
    cf_desc = t_cf_var('z0h_land', 'm', 'Roughness length for heat over land', &
        & datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'z0h_land', self%z0h_land, &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
        & ldims=shape2d, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.TRUE., in_group=groups('vdiff') &
      )
    !$ACC ENTER DATA ASYNC(1) ATTACH(self%z0h_land)

    ! self%z0m_sfc(nproma,nblks_c,SFC_NUM), &
    cf_desc = t_cf_var('z0m_sfc', 'm', 'Momentum roughness length over surface classes', &
        & datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'z0m_sfc', self%z0m_sfc, &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
        & ldims=shape2d_sfc, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.FALSE., loutput=.FALSE., lcontainer=.TRUE. &
      )
    !$ACC ENTER DATA ASYNC(1) ATTACH(self%z0m_sfc)

    DO i = 1, SFC_NUM
      cf_desc = t_cf_var('z0m_' // SFC_NAMES(i), 'm', &
          & 'Momentum roughness length over ' // SFC_NAMES(i), datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_ref(varlist, 'z0m_sfc', 'z0m_' // SFC_NAMES(i), p2d, &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ref_idx=i, ldims=shape2d, &
          & loutput=.TRUE., lrestart=.TRUE., in_group=groups('vdiff-sft') &
        )
    END DO

    CALL self%sea_state%init(nproma, nblks_c, varlist)

  END SUBROUTINE nwp_vdiff_state_init

  !> Perform a device-to-host update of the state structure.
  SUBROUTINE nwp_vdiff_state_d2h (self)

    CLASS(t_nwp_vdiff_state), INTENT(INOUT) :: self

    CALL self%sea_state%d2h()

    !$ACC UPDATE ASYNC(1) &
    !$ACC   HOST(self%exchange_coeff_h) &
    !$ACC   HOST(self%exchange_coeff_m) &
    !$ACC   HOST(self%exchange_coeff_h_sfc) &
    !$ACC   HOST(self%exchange_coeff_m_sfc) &
    !$ACC   HOST(self%fact_q_air) &
    !$ACC   HOST(self%fact_qsat_srf) &
    !$ACC   HOST(self%flx_co2_natural_land) &
    !$ACC   HOST(self%flx_heat_latent_sft) &
    !$ACC   HOST(self%flx_heat_sensible_sft) &
    !$ACC   HOST(self%fr_ice_on_lake) &
    !$ACC   HOST(self%temp_sft) &
    !$ACC   HOST(self%theta_v_var) &
    !$ACC   HOST(self%total_turbulence_energy) &
    !$ACC   HOST(self%total_water_var) &
    !$ACC   HOST(self%ustar) &
    !$ACC   HOST(self%wstar_sfc) &
    !$ACC   HOST(self%z0m_sfc) &
    !$ACC   HOST(self%z0h_land)
    !$ACC WAIT(1)

  END SUBROUTINE nwp_vdiff_state_d2h

  !> Perform a host-to-device update of the state structure.
  SUBROUTINE nwp_vdiff_state_h2d (self)

    CLASS(t_nwp_vdiff_state), INTENT(INOUT) :: self

    CALL self%sea_state%h2d()

    !$ACC UPDATE ASYNC(1) &
    !$ACC   DEVICE(self%exchange_coeff_h) &
    !$ACC   DEVICE(self%exchange_coeff_m) &
    !$ACC   DEVICE(self%exchange_coeff_h_sfc) &
    !$ACC   DEVICE(self%exchange_coeff_m_sfc) &
    !$ACC   DEVICE(self%fact_q_air) &
    !$ACC   DEVICE(self%fact_qsat_srf) &
    !$ACC   DEVICE(self%flx_co2_natural_land) &
    !$ACC   DEVICE(self%flx_heat_latent_sft) &
    !$ACC   DEVICE(self%flx_heat_sensible_sft) &
    !$ACC   DEVICE(self%fr_ice_on_lake) &
    !$ACC   DEVICE(self%temp_sft) &
    !$ACC   DEVICE(self%theta_v_var) &
    !$ACC   DEVICE(self%total_turbulence_energy) &
    !$ACC   DEVICE(self%total_water_var) &
    !$ACC   DEVICE(self%ustar) &
    !$ACC   DEVICE(self%wstar_sfc) &
    !$ACC   DEVICE(self%z0m_sfc) &
    !$ACC   DEVICE(self%z0h_land)

  END SUBROUTINE nwp_vdiff_state_h2d


  !> Destroy a `t_nwp_vdiff_state` structure.
  SUBROUTINE nwp_vdiff_state_destroy (self)

    TYPE(t_nwp_vdiff_state), INTENT(INOUT) :: self !< Object to destroy.

    ! This routine gets called when entering init. Thanks, Fortran!
    IF (.NOT. ASSOCIATED(self%fact_q_air)) RETURN

    !$ACC WAIT

    !$ACC EXIT DATA FINALIZE DELETE(self)

  END SUBROUTINE nwp_vdiff_state_destroy


  !> Initialize a `t_nwp_vdiff_albedos` structure. This routine is OpenMP orphaned.
  SUBROUTINE nwp_vdiff_albedos_init (self, kproma, nblks_c)

    CLASS(t_nwp_vdiff_albedos), INTENT(OUT) :: self !< Object to initialize.
    INTEGER, INTENT(IN) :: kproma !< Block size.
    INTEGER, INTENT(IN) :: nblks_c !< Number of cell blocks.

    !$OMP SINGLE
      ALLOCATE( &
          & self%alb_nir_dif(kproma,nblks_c,SFT_NUM), &
          & self%alb_nir_dir(kproma,nblks_c,SFT_NUM), &
          & self%alb_vis_dif(kproma,nblks_c,SFT_NUM), &
          & self%alb_vis_dir(kproma,nblks_c,SFT_NUM), &
          & self%lw_emissivity(kproma,nblks_c,SFT_NUM) &
        )
      self%deallocate = .TRUE.
    !$OMP END SINGLE

    !$OMP WORKSHARE
      self%alb_nir_dif(:,:,:) = 0._wp
      self%alb_nir_dir(:,:,:) = 0._wp
      self%alb_vis_dif(:,:,:) = 0._wp
      self%alb_vis_dir(:,:,:) = 0._wp
      self%lw_emissivity(:,:,:) = 0._wp
    !$OMP END WORKSHARE

    !$ACC ENTER DATA ASYNC(1) COPYIN(self, self%alb_nir_dif, self%alb_nir_dir) &
    !$ACC   COPYIN(self%alb_vis_dif, self%alb_vis_dir, self%lw_emissivity)

  END SUBROUTINE nwp_vdiff_albedos_init


  !> Perform a device-to-host update of the albedo structure.
  SUBROUTINE nwp_vdiff_albedos_d2h (self)

    CLASS(t_nwp_vdiff_albedos), INTENT(INOUT) :: self !< Object.

    ! Suppress unused warnings.
    IF (ASSOCIATED(self%lw_emissivity)) THEN ; END IF

    !$ACC UPDATE ASYNC(1) &
    !$ACC   HOST(self%alb_nir_dif) &
    !$ACC   HOST(self%alb_nir_dir) &
    !$ACC   HOST(self%alb_vis_dif) &
    !$ACC   HOST(self%alb_vis_dir) &
    !$ACC   HOST(self%lw_emissivity)
    !$ACC WAIT(1)

  END SUBROUTINE nwp_vdiff_albedos_d2h

  !> Perform a host-to-device update of the albedo structure.
  SUBROUTINE nwp_vdiff_albedos_h2d (self)

    CLASS(t_nwp_vdiff_albedos), INTENT(INOUT) :: self !< Object.

    ! Suppress unused warnings.
    IF (ASSOCIATED(self%lw_emissivity)) THEN ; END IF

    !$ACC UPDATE ASYNC(1) &
    !$ACC   DEVICE(self%alb_nir_dif) &
    !$ACC   DEVICE(self%alb_nir_dir) &
    !$ACC   DEVICE(self%alb_vis_dif) &
    !$ACC   DEVICE(self%alb_vis_dir) &
    !$ACC   DEVICE(self%lw_emissivity)

  END SUBROUTINE nwp_vdiff_albedos_h2d

  !> Copy assignment for `t_nwp_vdiff_albedos`. This is done on the device if OpenACC is enabled.
  SUBROUTINE nwp_vdiff_albedos_copy (self, src)

    CLASS(t_nwp_vdiff_albedos), INTENT(INOUT) :: self !< Target.
    CLASS(t_nwp_vdiff_albedos), INTENT(IN) :: src !< Source.

    INTEGER :: i1, i2, i3
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(3)
      DO i3 = LBOUND(self%alb_nir_dif, 3), UBOUND(self%alb_nir_dif, 3)
        DO i2 = LBOUND(self%alb_nir_dif, 2), UBOUND(self%alb_nir_dif, 2)
          DO i1 = LBOUND(self%alb_nir_dif, 1), UBOUND(self%alb_nir_dif, 1)
            self%alb_nir_dif(i1,i2,i3) = src%alb_nir_dif(i1,i2,i3)
            self%alb_nir_dir(i1,i2,i3) = src%alb_nir_dir(i1,i2,i3)
            self%alb_vis_dif(i1,i2,i3) = src%alb_vis_dif(i1,i2,i3)
            self%alb_vis_dir(i1,i2,i3) = src%alb_vis_dir(i1,i2,i3)
            self%lw_emissivity(i1,i2,i3) = src%lw_emissivity(i1,i2,i3)
          END DO
        END DO
      END DO
    !$ACC END PARALLEL

  END SUBROUTINE nwp_vdiff_albedos_copy

  !> Destroy a `t_nwp_vdiff_albedos` structure.
  SUBROUTINE nwp_vdiff_albedos_destroy (self)

    TYPE(t_nwp_vdiff_albedos), INTENT(INOUT) :: self !< Object to destroy.

    ! This routine gets called when entering init. Thanks, Fortran!
    IF (.NOT. ASSOCIATED(self%alb_nir_dif)) RETURN

    !$ACC WAIT
    IF (.NOT. self%deallocate) THEN
      !$ACC EXIT DATA FINALIZE DELETE(self)
    ELSE
      !$ACC EXIT DATA FINALIZE &
      !$ACC   DELETE(self%alb_nir_dif) &
      !$ACC   DELETE(self%alb_nir_dir) &
      !$ACC   DELETE(self%alb_vis_dif) &
      !$ACC   DELETE(self%alb_vis_dir) &
      !$ACC   DELETE(self%lw_emissivity) &
      !$ACC   DELETE(self)

      DEALLOCATE( &
          & self%alb_nir_dif, &
          & self%alb_nir_dir, &
          & self%alb_vis_dif, &
          & self%alb_vis_dir, &
          & self%lw_emissivity &
        )
    END IF

  END SUBROUTINE nwp_vdiff_albedos_destroy


  SUBROUTINE nwp_vdiff_sea_state_init (self, kproma, nblks_c, varlist)

    USE mo_time_config, ONLY: time_config
    USE mo_util_mtime, ONLY: assumePrevMidnight

    CLASS(t_nwp_vdiff_sea_state), INTENT(OUT) :: self !< Object to initialize.
    INTEGER, INTENT(IN) :: kproma !< Block size.
    INTEGER, INTENT(IN) :: nblks_c !< Number of cell blocks.
    TYPE(t_var_list_ptr), INTENT(INOUT) :: varlist !< Varlist to append variables.

    TYPE(t_cf_var) :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: datatype_flt
    INTEGER :: grib2_bits

    INTEGER :: shape2d(2) !< Shape of 2D fields.
    INTEGER :: shape2d_sft(3) !< Shape of 2D fields with surface-type index.

    INTEGER :: i

    REAL(wp), POINTER :: p2d(:,:)

    shape2d = [ kproma, nblks_c ]
    shape2d_sft = [ kproma, nblks_c, 2 ]

    IF (lnetcdf_flt64_output) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    END IF

    grib2_bits = DATATYPE_PACK16

    self%time_ref_t_seasfc = assumePrevMidnight(time_config%tc_exp_startdate)
    self%time_last_update_t_seasfc = self%time_ref_t_seasfc
    !$ACC ENTER DATA ASYNC(1) COPYIN(self)

    cf_desc = t_cf_var('t_seasfc_offset', 'K', 'Initial sea surface temperature', &
        & datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 't_seasfc_offset', self%t_seasfc_offset, &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
        & ldims=shape2d_sft, isteptype=TSTEP_CONSTANT, lopenacc=.TRUE., &
        & lrestart=.TRUE., loutput=.TRUE. &
      )
    !$ACC ENTER DATA ASYNC(1) ATTACH(self%t_seasfc_offset)

    cf_desc = t_cf_var('alb_nir_dif', '-', 'Surface albedo for NIR diffuse radiation', &
        & datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'alb_nir_dif', self%alb_nir_dif, &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
        & ldims=shape2d_sft, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.FALSE., loutput=.FALSE., lcontainer=.TRUE. &
      )

    cf_desc = t_cf_var('alb_nir_dir', '-', 'Surface albedo for NIR direct radiation', &
        & datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'alb_nir_dir', self%alb_nir_dir, &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
        & ldims=shape2d_sft, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.FALSE., loutput=.FALSE., lcontainer=.TRUE. &
      )

    cf_desc = t_cf_var('alb_vis_dif', '-', 'Surface albedo for VIS diffuse radiation', &
        & datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'alb_vis_dif', self%alb_vis_dif, &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
        & ldims=shape2d_sft, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.FALSE., loutput=.FALSE., lcontainer=.TRUE. &
      )

    cf_desc = t_cf_var('alb_vis_dir', '-', 'Surface albedo for VIS direct radiation', &
        & datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'alb_vis_dir', self%alb_vis_dir, &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
        & ldims=shape2d_sft, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.FALSE., loutput=.FALSE., lcontainer=.TRUE. &
      )

    cf_desc = t_cf_var('lw_emissivity', '-', 'Surface emissivity for long-wave radiation', &
        & datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(varlist, 'lw_emissivity', self%lw_emissivity, &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
        & ldims=shape2d_sft, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
        & lrestart=.FALSE., loutput=.FALSE., lcontainer=.TRUE. &
      )

    self%alb_nir_dif(1:,1:,SFT_SWTR:) => self%alb_nir_dif(:,:,:)
    self%alb_nir_dir(1:,1:,SFT_SWTR:) => self%alb_nir_dir(:,:,:)
    self%alb_vis_dif(1:,1:,SFT_SWTR:) => self%alb_vis_dif(:,:,:)
    self%alb_vis_dir(1:,1:,SFT_SWTR:) => self%alb_vis_dir(:,:,:)
    self%lw_emissivity(1:,1:,SFT_SWTR:) => self%lw_emissivity(:,:,:)

    !$ACC ENTER DATA ASYNC(1) &
    !$ACC   ATTACH(self%alb_nir_dif) &
    !$ACC   ATTACH(self%alb_nir_dir) &
    !$ACC   ATTACH(self%alb_vis_dif) &
    !$ACC   ATTACH(self%alb_vis_dir) &
    !$ACC   ATTACH(self%lw_emissivity)

    DO i = SFT_SWTR, SFT_SICE
      cf_desc = t_cf_var('alb_nir_dif_' // SFT_NAMES(i), '-', &
          & SFT_NAMES(i) // ' surface albedo for NIR diffuse radiation', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_ref(varlist, 'alb_nir_dif', 'alb_nir_dif_' // SFT_NAMES(i), p2d, &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ref_idx=i-SFT_SWTR+1, &
          & ldims=shape2d, loutput=.TRUE., lrestart=.TRUE., in_group=groups('vdiff-sft') &
        )

      cf_desc = t_cf_var('alb_nir_dir_' // SFT_NAMES(i), '-', &
          & SFT_NAMES(i) // ' surface albedo for NIR direct radiation', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_ref(varlist, 'alb_nir_dir', 'alb_nir_dir_' // SFT_NAMES(i), p2d, &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ref_idx=i-SFT_SWTR+1, &
          & ldims=shape2d, loutput=.TRUE., lrestart=.TRUE., in_group=groups('vdiff-sft') &
        )

      cf_desc = t_cf_var('alb_vis_dif_' // SFT_NAMES(i), '-', &
          & SFT_NAMES(i) // ' surface albedo for VIS diffuse radiation', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_ref(varlist, 'alb_vis_dif', 'alb_vis_dif_' // SFT_NAMES(i), p2d, &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ref_idx=i-SFT_SWTR+1, &
          & ldims=shape2d, loutput=.TRUE., lrestart=.TRUE., in_group=groups('vdiff-sft') &
        )

      cf_desc = t_cf_var('alb_vis_dir_' // SFT_NAMES(i), '-', &
          & SFT_NAMES(i) // ' surface albedo for VIS direct radiation', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_ref(varlist, 'alb_vis_dir', 'alb_vis_dir_' // SFT_NAMES(i), p2d, &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ref_idx=i-SFT_SWTR+1, &
          & ldims=shape2d, loutput=.TRUE., lrestart=.TRUE., in_group=groups('vdiff-sft') &
        )

      cf_desc = t_cf_var('lw_emissivity_' // SFT_NAMES(i), '-', &
          & SFT_NAMES(i) // ' surface emissivity for long-wave radiation', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_ref(varlist, 'lw_emissivity', 'lw_emissivity_' // SFT_NAMES(i), p2d, &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ref_idx=i-SFT_SWTR+1, &
          & ldims=shape2d, loutput=.TRUE., lrestart=.TRUE., in_group=groups('vdiff-sft') &
        )
    END DO

    IF (is_coupled_to_ocean()) THEN
      cf_desc = t_cf_var('flx_co2_natural_sea', 'kg m-2 s-1', 'Natural sea surface CO2 flux', &
          & datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var(varlist, 'flx_co2_natural_sea', self%flx_co2_natural_sea, &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
          & ldims=shape2d, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
          & lrestart=.TRUE., in_group=groups('vdiff') &
        )
      !$ACC ENTER DATA ASYNC(1) ATTACH(self%flx_co2_natural_sea)

      cf_desc = t_cf_var('ocean_u', 'm s-1', 'Zonal ocean surface velocity', &
          & datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var(varlist, 'ocean_u', self%ocean_u, &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
          & ldims=shape2d, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
          & lrestart=.TRUE., in_group=groups('vdiff') &
        )
      !$ACC ENTER DATA ASYNC(1) ATTACH(self%ocean_u)

      cf_desc = t_cf_var('ocean_v', 'm s-1', 'Meridional ocean surface velocity', &
          & datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, grib2_bits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var(varlist, 'ocean_v', self%ocean_v, &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
          & ldims=shape2d, isteptype=TSTEP_INSTANT, lopenacc=.TRUE., &
          & lrestart=.TRUE., in_group=groups('vdiff') &
        )
      !$ACC ENTER DATA ASYNC(1) ATTACH(self%ocean_v)
    END IF

  END SUBROUTINE nwp_vdiff_sea_state_init


  SUBROUTINE nwp_vdiff_sea_state_d2h (self)

    CLASS(t_nwp_vdiff_sea_state), INTENT(INOUT) :: self

    !$ACC UPDATE ASYNC(1) &
    !$ACC   HOST(self%alb_nir_dif) &
    !$ACC   HOST(self%alb_nir_dir) &
    !$ACC   HOST(self%alb_vis_dif) &
    !$ACC   HOST(self%alb_vis_dir) &
    !$ACC   HOST(self%lw_emissivity) &
    !$ACC   HOST(self%t_seasfc_offset) &
    !$ACC   HOST(self%time_ref_t_seasfc) &
    !$ACC   HOST(self%time_last_update_t_seasfc)

    !$ACC UPDATE ASYNC(1) IF(ASSOCIATED(self%ocean_u)) &
    !$ACC   HOST(self%flx_co2_natural_sea) &
    !$ACC   HOST(self%ocean_u) &
    !$ACC   HOST(self%ocean_v)

    !$ACC WAIT(1)

  END SUBROUTINE nwp_vdiff_sea_state_d2h

  !> Perform a host-to-device update of the state structure.
  SUBROUTINE nwp_vdiff_sea_state_h2d (self)

    CLASS(t_nwp_vdiff_sea_state), INTENT(INOUT) :: self

    !$ACC UPDATE ASYNC(1) &
    !$ACC   DEVICE(self%alb_nir_dif) &
    !$ACC   DEVICE(self%alb_nir_dir) &
    !$ACC   DEVICE(self%alb_vis_dif) &
    !$ACC   DEVICE(self%alb_vis_dir) &
    !$ACC   DEVICE(self%lw_emissivity) &
    !$ACC   DEVICE(self%t_seasfc_offset) &
    !$ACC   DEVICE(self%time_ref_t_seasfc) &
    !$ACC   DEVICE(self%time_last_update_t_seasfc)

    !$ACC UPDATE ASYNC(1) IF(ASSOCIATED(self%ocean_u)) &
    !$ACC   DEVICE(self%flx_co2_natural_sea) &
    !$ACC   DEVICE(self%ocean_u) &
    !$ACC   DEVICE(self%ocean_v)

  END SUBROUTINE nwp_vdiff_sea_state_h2d


  !> Destroy a `t_nwp_vdiff_sea_state` structure.
  SUBROUTINE nwp_vdiff_sea_state_destroy (self)

    TYPE(t_nwp_vdiff_sea_state), INTENT(INOUT) :: self !< Object to destroy.

    ! This routine gets called when entering init. Thanks, Fortran!
    IF (.NOT. ASSOCIATED(self%t_seasfc_offset)) RETURN

    !$ACC WAIT

    !$ACC EXIT DATA FINALIZE DELETE(self)

  END SUBROUTINE nwp_vdiff_sea_state_destroy

END MODULE mo_nwp_vdiff_types
