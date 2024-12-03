! # 1 "radiation/radiation_ecckd_interface.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_ecckd_interface.f90"
! radiation_ecckd_interface.f90 - interface to ecckd gas optics model
!
! (c) copyright 2020- ecmwf.
!
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.
!
! in applying this licence, ecmwf does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! author:  robin hogan
! email:   r.j.hogan@ecmwf.int
! license: see the copying file for details
!

module radiation_ecckd_interface

  implicit none

  public  :: setup_gas_optics, set_gas_units, gas_optics !, planck_function

contains

  !---------------------------------------------------------------------
  ! setup the ecckd generalized gas optics model
  subroutine setup_gas_optics(config)

    use parkind1, only : jprb
    use radiation_config
    use ecradhook,  only : lhook, dr_hook, jphook

    type(config_type), intent(inout), target :: config

    integer :: jj
    
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_ecckd_interface:setup_gas_optics',0,hook_handle)

    if (config%do_sw .and. config%i_gas_model_sw == igasmodelecckd) then

      ! read shortwave ecckd gas optics netcdf file
      call config%gas_optics_sw%read(trim(config%gas_optics_sw_file_name), &
           &                         config%iverbosesetup)

      ! copy over relevant properties
      config%n_g_sw     = config%gas_optics_sw%ng

      if (config%do_cloud_aerosol_per_sw_g_point) then
        ! bands and g points are the same
        config%n_bands_sw = config%n_g_sw
      else
        ! bands are groups of g points and span a continuous region of
        ! wavenumber space
        config%n_bands_sw = config%gas_optics_sw%spectral_def%nband
      end if

      allocate(config%i_band_from_g_sw          (config%n_g_sw))
      allocate(config%i_band_from_reordered_g_sw(config%n_g_sw))
      allocate(config%i_g_from_reordered_g_sw   (config%n_g_sw))
        
      if (config%do_cloud_aerosol_per_sw_g_point) then
        config%i_band_from_g_sw           = [ (jj, jj = 1,config%n_g_sw) ]
        config%i_band_from_reordered_g_sw = [ (jj, jj = 1,config%n_g_sw) ]
      else
        config%i_band_from_g_sw &
             &  = config%gas_optics_sw%spectral_def%i_band_number
        config%i_band_from_reordered_g_sw &
             &  = config%gas_optics_sw%spectral_def%i_band_number
      end if
      config%i_g_from_reordered_g_sw      = [ (jj, jj = 1,config%n_g_sw) ]

      if (config%iverbosesetup >= 2) then
        call config%gas_optics_sw%print()
      end if

      ! override solar spectral irradiance, if filename provided
      if (config%use_spectral_solar_cycle) then
        call config%gas_optics_sw%read_spectral_solar_cycle(config%ssi_file_name, &
             &           config%iverbosesetup, config%use_updated_solar_spectrum)
      end if

    end if

    if (config%do_lw .and. config%i_gas_model_lw == igasmodelecckd) then

      ! read longwave ecckd gas optics netcdf file
      call config%gas_optics_lw%read(trim(config%gas_optics_lw_file_name), &
           &                         config%iverbosesetup)

      ! copy over relevant properties
      config%n_g_lw     = config%gas_optics_lw%ng

      if (config%do_cloud_aerosol_per_lw_g_point) then
        ! bands and g points are the same
        config%n_bands_lw = config%n_g_lw
      else
        ! bands are groups of g points and span a continuous region of
        ! wavenumber space
        config%n_bands_lw = config%gas_optics_lw%spectral_def%nband
      end if

      allocate(config%i_band_from_g_lw          (config%n_g_lw))
      allocate(config%i_band_from_reordered_g_lw(config%n_g_lw))
      allocate(config%i_g_from_reordered_g_lw   (config%n_g_lw))

      if (config%do_cloud_aerosol_per_lw_g_point) then
        config%i_band_from_g_lw           = [ (jj, jj = 1,config%n_g_lw) ]
        config%i_band_from_reordered_g_lw = [ (jj, jj = 1,config%n_g_lw) ]
      else
        config%i_band_from_g_lw &
             &  = config%gas_optics_lw%spectral_def%i_band_number
        config%i_band_from_reordered_g_lw &
             &  = config%gas_optics_lw%spectral_def%i_band_number
      end if
      config%i_g_from_reordered_g_lw      = [ (jj, jj = 1,config%n_g_lw) ]

      if (config%iverbosesetup >= 2) then
        call config%gas_optics_lw%print()
      end if

    end if

    ! the i_spec_* variables are used solely for storing spectral
    ! data, and this can either be by band or by g-point
    if (config%do_save_spectral_flux) then
      if (config%do_save_gpoint_flux) then
        config%n_spec_sw = config%n_g_sw
        config%n_spec_lw = config%n_g_lw
        config%i_spec_from_reordered_g_sw => config%i_g_from_reordered_g_sw
        config%i_spec_from_reordered_g_lw => config%i_g_from_reordered_g_lw
      else
        config%n_spec_sw = config%n_bands_sw
        config%n_spec_lw = config%n_bands_lw
        config%i_spec_from_reordered_g_sw => config%i_band_from_reordered_g_sw
        config%i_spec_from_reordered_g_lw => config%i_band_from_reordered_g_lw
      end if
    else
      config%n_spec_sw = 0
      config%n_spec_lw = 0
      nullify(config%i_spec_from_reordered_g_sw)
      nullify(config%i_spec_from_reordered_g_lw)
    end if

    if (lhook) call dr_hook('radiation_ecckd_interface:setup_gas_optics',1,hook_handle)
    
  end subroutine setup_gas_optics


  !---------------------------------------------------------------------
  ! scale gas mixing ratios according to required units
  subroutine set_gas_units(gas)

    use radiation_gas, only : gas_type, ivolumemixingratio
    use ecradhook,       only : lhook, dr_hook, jphook
    
    type(gas_type),    intent(inout) :: gas

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_ecckd_interface:set_gas_units',0,hook_handle)

    call gas%set_units(ivolumemixingratio)

    if (lhook) call dr_hook('radiation_ecckd_interface:set_gas_units',1,hook_handle)

  end subroutine set_gas_units


  !---------------------------------------------------------------------
  ! compute gas optical depths, shortwave scattering, planck function
  ! and incoming shortwave radiation at top-of-atmosphere
  subroutine gas_optics(ncol,nlev,istartcol,iendcol, &
       &  config, single_level, thermodynamics, gas, & 
       &  od_lw, od_sw, ssa_sw, lw_albedo, planck_hl, lw_emission, &
       &  incoming_sw)

    use parkind1, only : jprb
    use ecradhook,  only : lhook, dr_hook, jphook

    use radiation_config,         only : config_type, igasmodelecckd
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_single_level,   only : single_level_type
    use radiation_gas_constants,  only : nmaxgases
    use radiation_gas

    integer,                  intent(in) :: ncol               ! number of columns
    integer,                  intent(in) :: nlev               ! number of levels
    integer,                  intent(in) :: istartcol, iendcol ! range of cols to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(gas_type),           intent(in) :: gas

    ! longwave albedo of the surface
    real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
         &  intent(in), optional :: lw_albedo

    ! gaseous layer optical depth in longwave and shortwave, and
    ! shortwave single scattering albedo (i.e. fraction of extinction
    ! due to rayleigh scattering) at each g-point
    real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_lw
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_sw, ssa_sw

    ! the planck function (emitted flux from a black body) at half
    ! levels at each longwave g-point
    real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol), &
         &   intent(out), optional :: planck_hl
    ! planck function for the surface (w m-2)
    real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
         &   intent(out), optional :: lw_emission

    ! the incoming shortwave flux into a plane perpendicular to the
    ! incoming radiation at top-of-atmosphere in each of the shortwave
    ! g-points
    real(jprb), dimension(config%n_g_sw,istartcol:iendcol), &
         &   intent(out), optional :: incoming_sw

    ! temperature at full levels (k)
    real(jprb) :: temperature_fl(istartcol:iendcol,nlev)

    real(jprb) :: concentration_scaling(nmaxgases)
    
    logical :: is_volume_mixing_ratio
    
    integer :: jcol, jlev, jg

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_ecckd_interface:gas_optics',0,hook_handle)

    !temperature_fl(istartcol:iendcol,:) &
    !     &  = 0.5_jprb * (thermodynamics%temperature_hl(istartcol:iendcol,1:nlev) &
    !     &               +thermodynamics%temperature_hl(istartcol:iendcol,2:nlev+1))
 
    temperature_fl(istartcol:iendcol,:) &
         &  = (thermodynamics%temperature_hl(istartcol:iendcol,1:nlev) &
         &     *thermodynamics%pressure_hl(istartcol:iendcol,1:nlev) &
         &    +thermodynamics%temperature_hl(istartcol:iendcol,2:nlev+1) &
         &     *thermodynamics%pressure_hl(istartcol:iendcol,2:nlev+1)) &
         &  / (thermodynamics%pressure_hl(istartcol:iendcol,1:nlev) &
         &    +thermodynamics%pressure_hl(istartcol:iendcol,2:nlev+1))

    ! check that the gas concentrations are stored in volume mixing
    ! ratio with no scaling; if not, return a vector of scalings
    call gas%assert_units(ivolumemixingratio, scale_factor=1.0_jprb, &
         &                istatus=is_volume_mixing_ratio)
    if (.not. is_volume_mixing_ratio) then
      call gas%get_scaling(ivolumemixingratio, concentration_scaling)
    else
      concentration_scaling = 1.0_jprb
    end if
    
    if (config%do_sw .and. config%i_gas_model_sw == igasmodelecckd) then

      if (is_volume_mixing_ratio) then
        call config%gas_optics_sw%calc_optical_depth(ncol,nlev,istartcol,iendcol, &
             &  nmaxgases, thermodynamics%pressure_hl, &
             &  temperature_fl, gas%mixing_ratio, &
             &  od_sw, rayleigh_od_fl=ssa_sw)
      else
        call config%gas_optics_sw%calc_optical_depth(ncol,nlev,istartcol,iendcol, &
             &  nmaxgases, thermodynamics%pressure_hl, &
             &  temperature_fl, gas%mixing_ratio, &
             &  od_sw, rayleigh_od_fl=ssa_sw, opt_concentration_scaling=concentration_scaling)
      end if
      
      ! at this point od_sw = absorption optical depth and ssa_sw =
      ! rayleigh optical depth: convert to total optical depth and
      ! single-scattering albedo
      do jcol = istartcol,iendcol
        do jlev = 1, nlev
          do jg = 1, config%n_g_sw
            od_sw(jg,jlev,jcol)  = od_sw(jg,jlev,jcol) + ssa_sw(jg,jlev,jcol)
            ssa_sw(jg,jlev,jcol) = ssa_sw(jg,jlev,jcol) / od_sw(jg,jlev,jcol)
          end do
        end do
      end do

      if (present(incoming_sw)) then
        if (single_level%spectral_solar_cycle_multiplier == 0.0_jprb) then
          call config%gas_optics_sw%calc_incoming_sw(single_level%solar_irradiance, &
               &        incoming_sw)
        else
          call config%gas_optics_sw%calc_incoming_sw(single_level%solar_irradiance, &
               &        incoming_sw, single_level%spectral_solar_cycle_multiplier)
        end if
      end if

    end if

    if (config%do_lw .and. config%i_gas_model_lw == igasmodelecckd) then

      if (is_volume_mixing_ratio) then
        call config%gas_optics_lw%calc_optical_depth(ncol,nlev,istartcol,iendcol, &
             &  nmaxgases, thermodynamics%pressure_hl, &
             &  temperature_fl, gas%mixing_ratio, &
             &  od_lw)
      else
        call config%gas_optics_lw%calc_optical_depth(ncol,nlev,istartcol,iendcol, &
             &  nmaxgases, thermodynamics%pressure_hl, &
             &  temperature_fl, gas%mixing_ratio, &
             &  od_lw, opt_concentration_scaling=concentration_scaling)
      end if

      ! calculate the planck function for each g point
      do jcol = istartcol,iendcol
        call config%gas_optics_lw%calc_planck_function(nlev+1, &
             &  thermodynamics%temperature_hl(jcol,:), planck_hl(:,:,jcol))
      end do
      call config%gas_optics_lw%calc_planck_function(iendcol+1-istartcol, &
           &  single_level%skin_temperature(istartcol:iendcol), &
           &  lw_emission(:,:))
!nec$ forced_collapse
      lw_emission = lw_emission * (1.0_jprb - lw_albedo)

    end if

    if (lhook) call dr_hook('radiation_ecckd_interface:gas_optics',1,hook_handle)
    
  end subroutine gas_optics

  ! !---------------------------------------------------------------------
  ! ! externally facing function for computing the planck function
  ! ! without reference to any gas profile; typically this would be used
  ! ! for computing the emission by a surface.
  ! subroutine planck_function(config, temperature, planck_surf)

  !   use parkind1,                 only : jprb
  !   use radiation_config,         only : config_type

  !   type(config_type), intent(in) :: config
  !   real(jprb),        intent(in) :: temperature

  !   ! planck function of the surface (w m-2)
  !   real(jprb), dimension(config%n_g_lw), intent(out) :: planck_surf

  ! end subroutine planck_function

end module radiation_ecckd_interface
! #define __atomic_acquire 2
! #define __char_bit__ 8
! #define __float_word_order__ __order_little_endian__
! #define __order_little_endian__ 1234
! #define __order_pdp_endian__ 3412
! #define __gfc_real_10__ 1
! #define __finite_math_only__ 0
! #define __gnuc_patchlevel__ 0
! #define __gfc_int_2__ 1
! #define __sizeof_int__ 4
! #define __sizeof_pointer__ 8
! #define __gfortran__ 1
! #define __gfc_real_16__ 1
! #define __stdc_hosted__ 0
! #define __no_math_errno__ 1
! #define __sizeof_float__ 4
! #define __pic__ 2
! #define _language_fortran 1
! #define __sizeof_long__ 8
! #define __gfc_int_8__ 1
! #define __dynamic__ 1
! #define __sizeof_short__ 2
! #define __gnuc__ 13
! #define __sizeof_long_double__ 16
! #define __biggest_alignment__ 16
! #define __atomic_relaxed 0
! #define _lp64 1
! #define __ecrad_little_endian 1
! #define __gfc_int_1__ 1
! #define __order_big_endian__ 4321
! #define __byte_order__ __order_little_endian__
! #define __sizeof_size_t__ 8
! #define __pic__ 2
! #define __sizeof_double__ 8
! #define __atomic_consume 1
! #define __gnuc_minor__ 3
! #define __gfc_int_16__ 1
! #define __lp64__ 1
! #define __atomic_seq_cst 5
! #define __sizeof_long_long__ 8
! #define __atomic_acq_rel 4
! #define __atomic_release 3
! #define __version__ "13.3.0"

