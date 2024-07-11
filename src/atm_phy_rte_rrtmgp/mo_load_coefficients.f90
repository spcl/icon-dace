! The original file has been published under the BSD-3-Clause license with the copyright
!
! Copyright (c) 2015-2018, Atmospheric and Environmental Research and Regents of the University of Colorado.
! All rights reserved.
!
! The file has been modified and released with the following license information:
!
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
! The gas optics class used by RRMTGP needs to be initialized with data stored in a netCDF file.
!    RRTMGP itself doesn't include methods for reading the data so we don't conflict with users'
!    local environment. This module provides a straight-forward implementation of reading the data
!    and calling gas_optics%load().
!
! -------------------------------------------------------------------------------------------------!

module mo_load_coefficients
  !
  ! Modules for working with rte and rrtmgp
  !
  use mo_rte_kind,           only: wp, wl
  use mo_exception,          only: finish
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_cloud_optics,       only: ty_cloud_optics
  use mo_netcdf,             only: nf90_nowrite, nf90_noerr
  use mo_netcdf_parallel,    only: &
    p_nf90_open, p_nf90_close, p_nf90_inq_varid, &
    p_nf90_inq_dimid, p_nf90_inquire_dimension, &
    p_nf90_get_var

  ! --------------------------------------------------
  implicit none
  private
  integer :: fid
  interface load_and_init
    module procedure load_and_init_gas_optics_rrtmgp, load_and_init_cloud_optics_rrtmgp
  end interface

  public :: load_and_init, stop_on_err

contains
  subroutine stop_on_err(msg)
    character(len=*), intent(in) :: msg

    IF (msg /= '') CALL finish('ERROR in RRTMGP', msg)
  end subroutine
  !--------------------------------------------------------------------------------------------------------------------
  ! read optical coefficients from NetCDF file
  subroutine load_and_init_gas_optics_rrtmgp(kdist, filename, available_gases)
    class(ty_gas_optics_rrtmgp), intent(inout) :: kdist
    character(len=*),     intent(in   ) :: filename
    class(ty_gas_concs),  intent(in   ) :: available_gases ! Which gases does the host model have available?
    ! --------------------------------------------------
    !
    ! Variables that will be passed to gas_optics%load()
    !
    character(len=32), dimension(:), allocatable :: gas_names
    integer,  dimension(:,:,:),      allocatable :: key_species
    integer,  dimension(:,:  ),      allocatable :: band2gpt
    real(wp), dimension(:,:  ),      allocatable :: band_lims
    real(wp)                                     :: press_ref_trop, temp_ref_p, temp_ref_t
    real(wp), dimension(:      ),    allocatable :: press_ref
    real(wp), dimension(:      ),    allocatable :: temp_ref
    real(wp), dimension(:,:,:  ),    allocatable :: vmr_ref
    real(wp), dimension(:,:,:,:),    allocatable :: kmajor

    character(len=32), dimension(:),  allocatable :: gas_minor, identifier_minor
    character(len=32), dimension(:),  allocatable :: minor_gases_lower,               minor_gases_upper
    integer, dimension(:,:),          allocatable :: minor_limits_gpt_lower,          minor_limits_gpt_upper
    logical(wl), dimension(:),        allocatable :: minor_scales_with_density_lower, minor_scales_with_density_upper
    character(len=32), dimension(:),  allocatable :: scaling_gas_lower,               scaling_gas_upper
    logical(wl), dimension(:),        allocatable :: scale_by_complement_lower,       scale_by_complement_upper
    integer, dimension(:),            allocatable :: kminor_start_lower,              kminor_start_upper
    real(wp), dimension(:,:,:),       allocatable :: kminor_lower,                    kminor_upper

    real(wp), dimension(:,:,:  ), allocatable :: rayl_lower, rayl_upper
    real(wp), dimension(:      ), allocatable :: solar_quiet, solar_facular, solar_sunspot
    real(wp)                                  :: tsi_default, mg_default, sb_default
    real(wp), dimension(:,:    ), allocatable :: totplnk
    real(wp), dimension(:,:,:,:), allocatable :: planck_frac
    real(wp), dimension(:,:    ), allocatable :: optimal_angle_fit
    ! -----------------
    !
    ! Book-keeping variables
    !
    integer :: nf_status
    integer :: ntemps,          &
               npress,          &
               nabsorbers,      &
               nextabsorbers,   &
               nminorabsorbers, &
               nmixingfracs,    &
               nlayers,         &
               nbnds,           &
               ngpts,           &
               npairs,          &
               nminor_absorber_intervals_lower, &
               nminor_absorber_intervals_upper, &
               ncontributors_lower, &
               ncontributors_upper, &
               ninternalSourcetemps, &
               string_len , &
               nfit_coeffs
    ! --------------------------------------------------
    !
    ! How big are the various arrays?
    !
    call open_file(trim(fileName), nf_status)
    if(nf_status /= nf90_noerr) then
      call stop_on_err("load_and_init(): can't open file " // trim(fileName))
    endif
    ntemps            = get_dim_size('temperature')
    npress            = get_dim_size('pressure')
    nabsorbers        = get_dim_size('absorber')
    nminorabsorbers   = get_dim_size('minor_absorber')
    nextabsorbers     = get_dim_size('absorber_ext')
    nmixingfracs      = get_dim_size('mixing_fraction')
    nlayers           = get_dim_size('atmos_layer')
    nbnds             = get_dim_size('bnd')
    ngpts             = get_dim_size('gpt')
    npairs            = get_dim_size('pair')
    nminor_absorber_intervals_lower &
                      = get_dim_size('minor_absorber_intervals_lower')
    nminor_absorber_intervals_upper  &
                      = get_dim_size('minor_absorber_intervals_upper')
    ninternalSourcetemps &
                      = get_dim_size('temperature_Planck')
    ncontributors_lower = get_dim_size('contributors_lower')
    ncontributors_upper = get_dim_size('contributors_upper')
    string_len          = get_dim_size('string_len')
    IF (string_len == 0) THEN ! SW file uses a different variable name (???)
      string_len        = get_dim_size('string32')
    END IF
    nfit_coeffs         = get_dim_size('fit_coeffs') ! Will be 0 for SW
    ! -----------------
    !
    ! Read the many arrays
    !
    call read_char('gas_names', string_len, nabsorbers, gas_names, nf_status)
    call read_int_3('key_species', 2, nlayers, nbnds, key_species, nf_status)
    call read_double_2('bnd_limits_wavenumber', 2, nbnds, band_lims, nf_status)
    call read_int_2('bnd_limits_gpt', 2, nbnds, band2gpt, nf_status)
    call read_double_1('press_ref', npress, press_ref, nf_status)
    call read_double_1('temp_ref',  ntemps, temp_ref, nf_status)
    call read_double_0('absorption_coefficient_ref_P', temp_ref_p, nf_status)
    call read_double_0('absorption_coefficient_ref_T', temp_ref_t, nf_status)
    call read_double_0('press_ref_trop', press_ref_trop, nf_status)
    call read_double_3('kminor_lower', ncontributors_lower, nmixingfracs, ntemps, &
      kminor_lower, nf_status)
    call read_double_3('kminor_upper', ncontributors_upper, nmixingfracs, ntemps, &
      kminor_upper, nf_status)
    call read_char('gas_minor', string_len, nminorabsorbers, gas_minor, nf_status)
    call read_char('identifier_minor', string_len, nminorabsorbers, identifier_minor, nf_status)
    call read_char('minor_gases_lower', string_len, nminor_absorber_intervals_lower, &
      minor_gases_lower, nf_status)
    call read_char('minor_gases_upper', string_len, nminor_absorber_intervals_upper, &
      minor_gases_upper, nf_status)
    call read_int_2('minor_limits_gpt_lower', npairs, &
      nminor_absorber_intervals_lower, minor_limits_gpt_lower, nf_status)
    call read_int_2('minor_limits_gpt_upper', npairs, &
      nminor_absorber_intervals_upper, minor_limits_gpt_upper, nf_status)
    call read_logical_1('minor_scales_with_density_lower', &
      nminor_absorber_intervals_lower, minor_scales_with_density_lower, nf_status)
    call read_logical_1('minor_scales_with_density_upper', &
      nminor_absorber_intervals_upper, minor_scales_with_density_upper, nf_status)
    call read_logical_1('scale_by_complement_lower', &
      nminor_absorber_intervals_lower, scale_by_complement_lower, nf_status)
    call read_logical_1('scale_by_complement_upper', &
      nminor_absorber_intervals_upper, scale_by_complement_upper, nf_status)
    call read_char('scaling_gas_lower', string_len, nminor_absorber_intervals_lower, &
      scaling_gas_lower, nf_status)
    call read_char('scaling_gas_upper', string_len, nminor_absorber_intervals_upper, &
      scaling_gas_upper, nf_status)
    call read_int_1('kminor_start_lower', nminor_absorber_intervals_lower, &
      kminor_start_lower, nf_status)
    call read_int_1('kminor_start_upper', nminor_absorber_intervals_upper, &
      kminor_start_upper, nf_status)
    call read_double_3('vmr_ref', nlayers, nextabsorbers, ntemps, vmr_ref, nf_status)

    call read_double_4('kmajor',  ngpts, nmixingfracs,  npress+1, ntemps, kmajor, &
      nf_status)
    if(var_exists('rayl_lower')) then
      call read_double_3('rayl_lower', ngpts, nmixingfracs, ntemps, rayl_lower, &
        nf_status)
      call read_double_3('rayl_upper', ngpts, nmixingfracs, ntemps, rayl_upper, &
        nf_status)
    end if
    ! --------------------------------------------------
    !
    ! Initialize the gas optics class with data. The calls look slightly different depending
    !   on whether the radiation sources are internal to the atmosphere (longwave) or external (shortwave)
    ! gas_optics%load() returns a string; a non-empty string indicates an error.
    !
    if(var_exists('totplnk')) then
      !
      ! If there's a totplnk variable in the file it's a longwave (internal sources) type
      !
      call read_double_2('totplnk', ninternalSourcetemps, nbnds, totplnk, nf_status)
      call read_double_4('plank_fraction', ngpts, nmixingfracs, npress+1, ntemps, &
        planck_frac, nf_status)
      call read_double_2('optimal_angle_fit', nfit_coeffs, nbnds, optimal_angle_fit, nf_status)

      call stop_on_err(kdist%load(available_gases, &
                                  gas_names,   &
                                  key_species, &
                                  band2gpt,    &
                                  band_lims,   &
                                  press_ref,   &
                                  press_ref_trop, &
                                  temp_ref,    &
                                  temp_ref_p, temp_ref_t,     &
                                  vmr_ref, kmajor,            &
                                  kminor_lower, kminor_upper, &
                                  gas_minor,identifier_minor, &
                                  minor_gases_lower, minor_gases_upper, &
                                  minor_limits_gpt_lower, &
                                  minor_limits_gpt_upper, &
                                  minor_scales_with_density_lower, &
                                  minor_scales_with_density_upper, &
                                  scaling_gas_lower, scaling_gas_upper, &
                                  scale_by_complement_lower, &
                                  scale_by_complement_upper, &
                                  kminor_start_lower, &
                                  kminor_start_upper, &
                                  totplnk, planck_frac,       &
                                  rayl_lower, rayl_upper, &
                                  optimal_angle_fit))
    else
      !
      ! Solar source doesn't have an dependencies yet
      !
      call read_double_1('solar_source_quiet',   ngpts, solar_quiet,   nf_status)
      call read_double_1('solar_source_facular', ngpts, solar_facular, nf_status)
      call read_double_1('solar_source_sunspot', ngpts, solar_sunspot, nf_status)
      call read_double_0('tsi_default', tsi_default, nf_status)
      call read_double_0('mg_default',  mg_default,  nf_status)
      call read_double_0('sb_default',  sb_default,  nf_status)

      call stop_on_err(kdist%load(available_gases, &
                                  gas_names,   &
                                  key_species, &
                                  band2gpt,    &
                                  band_lims,   &
                                  press_ref,   &
                                  press_ref_trop, &
                                  temp_ref,    &
                                  temp_ref_p, temp_ref_t,     &
                                  vmr_ref, kmajor,            &
                                  kminor_lower, kminor_upper, &
                                  gas_minor,identifier_minor,&
                                  minor_gases_lower, minor_gases_upper, &
                                  minor_limits_gpt_lower, &
                                  minor_limits_gpt_upper, &
                                  minor_scales_with_density_lower, &
                                  minor_scales_with_density_upper, &
                                  scaling_gas_lower, scaling_gas_upper, &
                                  scale_by_complement_lower, &
                                  scale_by_complement_upper, &
                                  kminor_start_lower, &
                                  kminor_start_upper, &
                                  solar_quiet, solar_facular, solar_sunspot, &
                                  tsi_default, mg_default, sb_default, &
                                  rayl_lower, rayl_upper))
    end if
    ! --------------------------------------------------
    call close_file(nf_status)
  end subroutine load_and_init_gas_optics_rrtmgp
!  ------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------
!
! read cloud optical property LUT coefficients from NetCDF file
!
subroutine load_and_init_cloud_optics_rrtmgp(cloud_spec, cld_coeff_file)
  class(ty_cloud_optics),     intent(inout) :: cloud_spec
  character(len=*),           intent(in   ) :: cld_coeff_file
  ! -----------------
  ! Local variables
  integer :: nf_status
  integer :: ncid, nband, nrghice, nsize_liq, nsize_ice

  real(wp), dimension(:,:), allocatable                :: band_lims_wvn
  ! Lookup table interpolation constants
  real(wp) :: radliq_lwr          ! liquid particle size lower bound for interpolation
  real(wp) :: radliq_upr          ! liquid particle size upper bound for interpolation
  real(wp) :: radliq_fac          ! constant for calculating interpolation indices for liquid
  real(wp) :: radice_lwr          ! ice particle size lower bound for interpolation
  real(wp) :: radice_upr          ! ice particle size upper bound for interpolation
  real(wp) :: radice_fac          ! constant for calculating interpolation indices for ice
  ! LUT coefficients
  real(wp), dimension(:,:),   allocatable :: lut_extliq   ! extinction: liquid
  real(wp), dimension(:,:),   allocatable :: lut_ssaliq   ! single scattering albedo: liquid
  real(wp), dimension(:,:),   allocatable :: lut_asyliq   ! asymmetry parameter: liquid
  real(wp), dimension(:,:,:), allocatable :: lut_extice   ! extinction: ice
  real(wp), dimension(:,:,:), allocatable :: lut_ssaice   ! single scattering albedo: ice
  real(wp), dimension(:,:,:), allocatable :: lut_asyice   ! asymmetry parameter: ice
  ! -----------------
  ! Open cloud optical property coefficient file
  call open_file(trim(cld_coeff_file), nf_status)
  if(nf_status /= nf90_noerr) then
    call stop_on_err("load_and_init(): can't open file " // trim(cld_coeff_file))
  endif

  ! Read LUT coefficient dimensions
  nband     = get_dim_size('nband')
  nrghice   = get_dim_size('nrghice')
  nsize_liq = get_dim_size('nsize_liq')
  nsize_ice = get_dim_size('nsize_ice')

  allocate(band_lims_wvn(2, nband))
  call read_double_2('bnd_limits_wavenumber', 2, nband, band_lims_wvn, nf_status)

  ! Read LUT constants
  call read_double_0('radliq_lwr', radliq_lwr, nf_status)
  call read_double_0('radliq_upr', radliq_upr, nf_status)
  call read_double_0('radliq_fac', radliq_fac, nf_status)
  call read_double_0('radice_lwr', radice_lwr, nf_status)
  call read_double_0('radice_upr', radice_upr, nf_status)
  call read_double_0('radice_fac', radice_fac, nf_status)

  ! Allocate cloud property lookup table input arrays
  allocate(lut_extliq(nsize_liq, nband), &
           lut_ssaliq(nsize_liq, nband), &
           lut_asyliq(nsize_liq, nband), &
           lut_extice(nsize_ice, nband, nrghice), &
           lut_ssaice(nsize_ice, nband, nrghice), &
           lut_asyice(nsize_ice, nband, nrghice))

  ! Read LUT coefficients
  call read_double_2('lut_extliq',  nsize_liq, nband, lut_extliq, nf_status)
  call read_double_2('lut_ssaliq',  nsize_liq, nband, lut_ssaliq, nf_status)
  call read_double_2('lut_asyliq',  nsize_liq, nband, lut_asyliq, nf_status)
  call read_double_3('lut_extice',  nsize_ice, nband, nrghice, lut_extice, nf_status)
  call read_double_3('lut_ssaice',  nsize_ice, nband, nrghice, lut_ssaice, nf_status)
  call read_double_3('lut_asyice',  nsize_ice, nband, nrghice, lut_asyice, nf_status)

  call close_file(nf_status)

  call stop_on_err(cloud_spec%load(band_lims_wvn,                      &
                                   radliq_lwr, radliq_upr, radliq_fac, &
                                   radice_lwr, radice_upr, radice_fac, &
                                   lut_extliq, lut_ssaliq, lut_asyliq, &
                                   lut_extice, lut_ssaice, lut_asyice))
end subroutine load_and_init_cloud_optics_rrtmgp

!  ------------------------------------------------------------------------------
  subroutine open_file(filename, status)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: status

    status = p_nf90_open(trim(filename), nf90_nowrite, fid)
  end subroutine open_file
  ! --------------------------------------------------
  subroutine close_file(status)
    integer, intent(inout) :: status
    status = p_nf90_close(fid)
  end subroutine close_file
  ! --------------------------------------------------
  function var_exists(name)
    character(len=*), intent(in) :: name
    logical :: var_exists

    integer :: dummy
    var_exists = (p_nf90_inq_varid(fid, trim(name), dummy) == nf90_noerr)
  end function var_exists
  ! --------------------------------------------------
  function get_dim_size(name) result (ret)
    character(len=*), intent(in) :: name
    integer :: ret

    integer :: dimid, status

    status = p_nf90_inq_dimid(fid, trim(name), dimid)
    if (status /= nf90_noerr) then
      ret = 0
    else
      status = p_nf90_inquire_dimension(fid, dimid, len = ret)
      if (status /= nf90_noerr) then
        ret = 0
      endif
    endif
  end function get_dim_size
  ! --------------------------------------------------
  subroutine read_char(name, string_len, nx, ret, status)
    character(len=*), intent(in) :: name
    integer, intent(in) :: string_len, nx
    character(len=*), allocatable, intent(out) :: ret(:)
    integer, intent(out) :: status

    integer :: varid, i, l
    character(len=string_len), allocatable :: tmp(:)
    character(len=string_len) :: hold

    status = p_nf90_inq_varid(fid, trim(name), varid)
    if (status /= nf90_noerr) return
    allocate(ret(nx))
    allocate(tmp(nx))
    status = p_nf90_get_var(fid, varid, tmp, [1,1], [string_len,nx])

    do i = 1, nx
      l = len(trim(tmp(i)))
      hold = ' '
      hold(1:l) = tmp(i)
      ret(i) = hold
    enddo

  end subroutine read_char
  ! --------------------------------------------------
  subroutine read_double_0(name, ret, status)
    character(len=*), intent(in) :: name
    real(wp), intent(out) :: ret
    integer, intent(out) :: status

    integer :: varid
    real(wp) :: tmp(1)

    status = p_nf90_inq_varid(fid, trim(name), varid)
    if (status /= nf90_noerr) return
    status = p_nf90_get_var(fid, varid, tmp, [1], [1])
    ret = tmp(1)

  end subroutine read_double_0
  ! --------------------------------------------------
  subroutine read_double_1(name, nx, ret, status)
    character(len=*), intent(in) :: name
    integer, intent(in) :: nx
    real(wp), allocatable, intent(out) :: ret(:)
    integer, intent(out) :: status

    integer :: varid

    status = p_nf90_inq_varid(fid, trim(name), varid)
    if (status /= nf90_noerr) return
    allocate(ret(nx))
    status = p_nf90_get_var(fid, varid, ret, [1], [nx])

  end subroutine read_double_1
  ! --------------------------------------------------
  subroutine read_double_2(name, nx, ny, ret, status)
    character(len=*), intent(in) :: name
    integer, intent(in) :: nx, ny
    real(wp), allocatable, intent(out) :: ret(:,:)
    integer, intent(out) :: status

    integer :: varid

    status = p_nf90_inq_varid(fid, trim(name), varid)
    if (status /= nf90_noerr) return
    allocate(ret(nx,ny))
    status = p_nf90_get_var(fid, varid, ret, [1,1], [nx,ny])

  end subroutine read_double_2
  ! --------------------------------------------------
  subroutine read_double_3(name, nx, ny, nz, ret, status)
    character(len=*), intent(in) :: name
    integer, intent(in) :: nx, ny, nz
    real(wp), allocatable, intent(out) :: ret(:,:,:)
    integer, intent(out) :: status

    integer :: varid

    status = p_nf90_inq_varid(fid, trim(name), varid)
    if (status /= nf90_noerr) return
    allocate(ret(nx,ny,nz))
    status = p_nf90_get_var(fid, varid, ret, [1,1,1], [nx,ny,nz])

  end subroutine read_double_3
  ! --------------------------------------------------
  subroutine read_double_4(name, nx, ny, nz, nt, ret, status)
    character(len=*), intent(in) :: name
    integer, intent(in) :: nx, ny, nz, nt
    real(wp), allocatable, intent(out) :: ret(:,:,:,:)
    integer, intent(out) :: status

    integer :: varid

    status = p_nf90_inq_varid(fid, trim(name), varid)
    if (status /= nf90_noerr) return
    allocate(ret(nx,ny,nz,nt))
    status = p_nf90_get_var(fid, varid, ret, [1,1,1,1], [nx,ny,nz,nt])

  end subroutine read_double_4
  ! --------------------------------------------------
  subroutine read_int_1(name, nx, ret, status)
    character(len=*), intent(in) :: name
    integer, intent(in) :: nx
    integer, allocatable, intent(out) :: ret(:)
    integer, intent(out) :: status

    integer :: varid

    status = p_nf90_inq_varid(fid, trim(name), varid)
    if (status /= nf90_noerr) return
    allocate(ret(nx))
    status = p_nf90_get_var(fid, varid, ret, [1], [nx])

  end subroutine read_int_1
  ! --------------------------------------------------
  subroutine read_int_2(name, nx, ny, ret, status)
    character(len=*), intent(in) :: name
    integer, intent(in) :: nx, ny
    integer, allocatable, intent(out) :: ret(:,:)
    integer, intent(out) :: status

    integer :: varid

    status = p_nf90_inq_varid(fid, trim(name), varid)
    if (status /= nf90_noerr) return
    allocate(ret(nx,ny))
    status = p_nf90_get_var(fid, varid, ret, [1,1], [nx,ny])

  end subroutine read_int_2
  ! --------------------------------------------------
  subroutine read_int_3(name, nx, ny, nz, ret, status)
    character(len=*), intent(in) :: name
    integer, intent(in) :: nx, ny, nz
    integer, allocatable, intent(out) :: ret(:,:,:)
    integer, intent(out) :: status

    integer :: varid

    status = p_nf90_inq_varid(fid, trim(name), varid)
    if (status /= nf90_noerr) return
    allocate(ret(nx,ny,nz))
    status = p_nf90_get_var(fid, varid, ret, [1,1,1], [nx,ny,nz])

  end subroutine read_int_3
  ! --------------------------------------------------
  subroutine read_logical_1(name, nx, ret, status)
    character(len=*), intent(in) :: name
    integer, intent(in) :: nx
    logical(wl), allocatable, intent(out) :: ret(:)
    integer, intent(out) :: status
    integer, dimension(nx) :: tmp

    integer :: varid

    status = p_nf90_inq_varid(fid, trim(name), varid)
    if (status /= nf90_noerr) return
    status = p_nf90_get_var(fid, varid, tmp, [1], [nx])
    if (status /= nf90_noerr) return
    allocate(ret(nx))
    ret = .false.
    where (tmp /= 0) ret = .true.

  end subroutine read_logical_1
  ! --------------------------------------------------
end module
