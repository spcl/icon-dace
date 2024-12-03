! # 1 "radiation/radiation_monochromatic.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_monochromatic.f90"
! radiation_interface.f90 - monochromatic gas/cloud optics for testing
!
! (c) copyright 2014- ecmwf.
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
!
! modifications
!   2017-04-11  r. hogan  receive "surface" dummy argument
!   2017-09-13  r. hogan  revert
!   2018-08-29  r. hogan  particulate single-scattering albedo / asymmetry from namelist

module radiation_monochromatic

  implicit none

  public  :: setup_gas_optics, gas_optics, set_gas_units, &
       &     setup_cloud_optics, cloud_optics,            &
       &     setup_aerosol_optics, add_aerosol_optics

contains

  ! provides elemental function "delta_eddington"

! # 1 "radiation/radiation_delta_eddington.h" 1
! radiation_delta_eddington.h - delta-eddington scaling -*- f90 -*-
!
! (c) copyright 2015- ecmwf.
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
!
! this file is intended to be included inside a module to ensure that
! these simple functions may be inlined

!---------------------------------------------------------------------
! perform in-place delta-eddington scaling of the phase function
elemental subroutine delta_eddington(od, ssa, g)

  use parkind1, only : jprb
  
  ! total optical depth, single scattering albedo and asymmetry
  ! factor
  real(jprb), intent(inout) :: od, ssa, g
  
  ! fraction of the phase function deemed to be in the forward lobe
  ! and therefore treated as if it is not scattered at all
  real(jprb) :: f
  
  f   = g*g
  od  = od * (1.0_jprb - ssa*f)
  ssa = ssa * (1.0_jprb - f) / (1.0_jprb - ssa*f)
  g   = g / (1.0_jprb + g)
  
end subroutine delta_eddington


!---------------------------------------------------------------------
! perform in-place delta-eddington scaling of the phase function, but
! using extensive variables (i.e. the scattering optical depth,
! scat_od, rather than the single-scattering albedo, and the
! scattering-optical-depth-multiplied-by-asymmetry-factor, scat_od_g,
! rather than the asymmetry factor.
elemental subroutine delta_eddington_extensive(od, scat_od, scat_od_g)

  !$acc routine seq

  use parkind1, only : jprb

  ! total optical depth, scattering optical depth and asymmetry factor
  ! multiplied by the scattering optical depth
  real(jprb), intent(inout) :: od, scat_od, scat_od_g

  ! fraction of the phase function deemed to be in the forward lobe
  ! and therefore treated as if it is not scattered at all
  real(jprb) :: f, g

  if (scat_od > 0.0_jprb) then
    g = scat_od_g / scat_od
  else
    g = 0.0
  end if

  f         = g*g
  od        = od - scat_od * f
  scat_od   = scat_od * (1.0_jprb - f)
  scat_od_g = scat_od * g / (1.0_jprb + g)
  
end subroutine delta_eddington_extensive


!---------------------------------------------------------------------
! array version of delta_eddington_extensive, more likely to vectorize
 subroutine delta_eddington_extensive_vec(ng, od, scat_od, scat_od_g)

  use parkind1, only : jprb

  ! total optical depth, scattering optical depth and asymmetry factor
  ! multiplied by the scattering optical depth
  integer,                   intent(in)    :: ng
  real(jprb), dimension(ng), intent(inout) :: od, scat_od, scat_od_g

  ! fraction of the phase function deemed to be in the forward lobe
  ! and therefore treated as if it is not scattered at all
  real(jprb) :: f, g
  integer :: j

  do j = 1,ng
    g            = scat_od_g(j) / max(scat_od(j), 1.0e-24)
    f            = g*g
    od(j)        = od(j) - scat_od(j) * f
    scat_od(j)   = scat_od(j) * (1.0_jprb - f)
    scat_od_g(j) = scat_od(j) * g / (1.0_jprb + g)
  end do
  
end subroutine delta_eddington_extensive_vec


!---------------------------------------------------------------------
! perform in-place delta-eddington scaling of the phase function,
! using the scattering optical depth rather than the single scattering
! albedo
elemental subroutine delta_eddington_scat_od(od, scat_od, g)

  use parkind1, only : jprb
  
  ! total optical depth, scattering optical depth and asymmetry factor
  real(jprb), intent(inout) :: od, scat_od, g

  ! fraction of the phase function deemed to be in the forward lobe
  ! and therefore treated as if it is not scattered at all
  real(jprb) :: f

  !$acc routine seq

  f       = g*g
  od      = od - scat_od * f
  scat_od = scat_od * (1.0_jprb - f)
  g       = g / (1.0_jprb + g)

end subroutine delta_eddington_scat_od


!---------------------------------------------------------------------
! revert delta-eddington-scaled quantities in-place, back to their
! original state
elemental subroutine revert_delta_eddington(od, ssa, g)

  use parkind1, only : jprb
  
  ! total optical depth, single scattering albedo and asymmetry
  ! factor
  real(jprb), intent(inout) :: od, ssa, g
  
  ! fraction of the phase function deemed to be in the forward lobe
  ! and therefore treated as if it is not scattered at all
  real(jprb) :: f
  
  g   = g / (1.0_jprb - g)
  f   = g*g
  ssa = ssa / (1.0_jprb - f + f*ssa);
  od  = od / (1.0_jprb - ssa*f)
  
end subroutine revert_delta_eddington
! # 32 "radiation/radiation_monochromatic.f90" 2

  !---------------------------------------------------------------------
  ! setup the arrays in the config object corresponding to the
  ! monochromatic gas optics model.  the directory argument is not
  ! used, since no look-up tables need to be loaded.
  subroutine setup_gas_optics(config, directory)

    use radiation_config, only : config_type
    
    type(config_type), intent(inout) :: config
    character(len=*),  intent(in)    :: directory

    ! in the monochromatic model we have simply one band and g-point
    ! in both the longwave and shortwave parts of the spectrum
    config%n_g_sw     = 1
    config%n_g_lw     = 1
    config%n_bands_sw = 1
    config%n_bands_lw = 1

    ! allocate arrays
    allocate(config%i_band_from_g_sw          (config%n_g_sw))
    allocate(config%i_band_from_g_lw          (config%n_g_lw))
    allocate(config%i_band_from_reordered_g_sw(config%n_g_sw))
    allocate(config%i_band_from_reordered_g_lw(config%n_g_lw))
    allocate(config%i_g_from_reordered_g_sw(config%n_g_sw))
    allocate(config%i_g_from_reordered_g_lw(config%n_g_lw))

    ! indices are trivial...
    config%i_band_from_g_sw           = 1
    config%i_band_from_g_lw           = 1
    config%i_g_from_reordered_g_sw    = 1
    config%i_g_from_reordered_g_lw    = 1
    config%i_band_from_reordered_g_sw = 1
    config%i_band_from_reordered_g_lw = 1

  end subroutine setup_gas_optics


  !---------------------------------------------------------------------
  ! dummy routine for scaling gas mixing ratios
  subroutine set_gas_units(gas)

    use radiation_gas,           only : gas_type
    type(gas_type),    intent(inout) :: gas

  end subroutine set_gas_units


  !---------------------------------------------------------------------
  ! dummy setup routine for cloud optics: in fact, no setup is
  ! required for monochromatic case
  subroutine setup_cloud_optics(config)

    use radiation_config, only : config_type
    type(config_type), intent(inout) :: config

  end subroutine setup_cloud_optics


  !---------------------------------------------------------------------
  ! dummy subroutine since no aerosols are represented in
  ! monochromatic case
  subroutine setup_aerosol_optics(config)

    use radiation_config,              only : config_type
    type(config_type), intent(inout) :: config

  end subroutine setup_aerosol_optics


  !---------------------------------------------------------------------
  ! compute gas optical depths, shortwave scattering, planck function
  ! and incoming shortwave radiation at top-of-atmosphere
  subroutine gas_optics(ncol,nlev,istartcol,iendcol, &
       config, single_level, thermodynamics, gas, lw_albedo, & 
       od_lw, od_sw, ssa_sw, planck_hl, lw_emission, &
       incoming_sw)

    use parkind1,                 only : jprb
    use radiation_config,         only : config_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_single_level,   only : single_level_type
    use radiation_gas,            only : gas_type
    use radiation_constants,      only : pi, stefanboltzmann

    ! inputs
    integer, intent(in) :: ncol               ! number of columns
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(gas_type),           intent(in) :: gas

    ! longwave albedo of the surface
    real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
         &  intent(in) :: lw_albedo

    ! outputs

    ! gaseous layer optical depth in longwave and shortwave, and
    ! shortwave single scattering albedo (i.e. fraction of extinction
    ! due to rayleigh scattering) at each g-point
    real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_lw
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_sw, ssa_sw

    ! the planck function (emitted flux from a black body) at half
    ! levels and at the surface at each longwave g-point
    real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol), intent(out) :: &
         &   planck_hl
    real(jprb), dimension(config%n_g_lw,istartcol:iendcol), intent(out) :: &
         &   lw_emission

    ! the incoming shortwave flux into a plane perpendicular to the
    ! incoming radiation at top-of-atmosphere in each of the shortwave
    ! g-points
    real(jprb), dimension(config%n_g_sw,istartcol:iendcol), intent(out) :: &
         &   incoming_sw
    
    ! ratio of the optical depth of the entire atmospheric column that
    ! is in the current layer
    real(jprb), dimension(istartcol:iendcol) :: extinction_fraction

    ! in the monochromatic model, the absorption by the atmosphere is
    ! assumed proportional to the mass in each layer, so is defined in
    ! terms of a total zenith optical depth and then distributed with
    ! height according to the pressure.
    !real(jprb), parameter :: total_od_sw = 0.10536_jprb ! transmittance 0.9
    !real(jprb), parameter :: total_od_lw = 1.6094_jprb  ! transmittance 0.2

    integer :: jlev

    do jlev = 1,nlev
      ! the fraction of the total optical depth in the current layer
      ! is proportional to the fraction of the mass of the atmosphere
      ! in the current layer, computed from pressure assuming
      ! hydrostatic balance
      extinction_fraction = &
           &   (thermodynamics%pressure_hl(istartcol:iendcol,jlev+1) &
           &   -thermodynamics%pressure_hl(istartcol:iendcol,jlev)) &
           &   /thermodynamics%pressure_hl(istartcol:iendcol,nlev)
      
      ! assign longwave and shortwave optical depths
      od_lw(1,jlev,:) = config%mono_lw_total_od*extinction_fraction
      od_sw(1,jlev,:) = config%mono_sw_total_od*extinction_fraction
    end do

    ! shortwave single-scattering albedo is almost entirely rayleigh
    ! scattering
    ssa_sw = 0.999999_jprb

    ! entire shortwave spectrum represented in one band
    incoming_sw(1,:) = single_level%solar_irradiance

    if (single_level%is_simple_surface) then
      if (config%mono_lw_wavelength <= 0.0_jprb) then
        ! entire longwave spectrum represented in one band
        lw_emission(1,istartcol:iendcol) &
             &  = stefanboltzmann * single_level%skin_temperature(istartcol:iendcol)**4 &
             &  * single_level%lw_emissivity(istartcol:iendcol,1)
        do jlev = 1,nlev+1
          planck_hl(1,jlev,istartcol:iendcol) = stefanboltzmann * thermodynamics%temperature_hl(istartcol:iendcol,jlev)**4
        end do
      else
        ! single wavelength: multiply by pi to convert w sr-1 m-3 to w m-3
        lw_emission(1,istartcol:iendcol) = pi*planck_function(config%mono_lw_wavelength, &
             &             single_level%skin_temperature(istartcol:iendcol)) &
             &  * single_level%lw_emissivity(istartcol:iendcol,1)
        do jlev = 1,nlev+1
          planck_hl(1,jlev,istartcol:iendcol) = pi*planck_function(config%mono_lw_wavelength, &
               &             thermodynamics%temperature_hl(istartcol:iendcol,jlev))
        end do
      end if
    else
      lw_emission = transpose(single_level%lw_emission)
    end if

  end subroutine gas_optics


  !---------------------------------------------------------------------
  ! compute cloud optical depth, single-scattering albedo and
  ! g factor in the longwave and shortwave
  subroutine cloud_optics(nlev,istartcol,iendcol, &
       &   config, thermodynamics, cloud, & 
       &   od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
       &   od_sw_cloud, ssa_sw_cloud, g_sw_cloud)

    use parkind1,                 only : jprb
    use radiation_config,         only : config_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_cloud,          only : cloud_type
    use radiation_constants,      only : accelduetogravity, &
         &   densityliquidwater, densitysolidice

    ! inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in) :: config
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(cloud_type),   intent(in) :: cloud

    ! outputs

    ! layer optical depth, single scattering albedo and g factor of
    ! clouds in each longwave band, where the latter two
    ! variables are only defined if cloud longwave scattering is
    ! enabled (otherwise both are treated as zero).
    real(jprb), dimension(config%n_bands_lw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_lw_cloud
    real(jprb), dimension(config%n_bands_lw_if_scattering,nlev,istartcol:iendcol), &
         &   intent(out) :: ssa_lw_cloud, g_lw_cloud

    ! layer optical depth, single scattering albedo and g factor of
    ! clouds in each shortwave band
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_sw_cloud, ssa_sw_cloud, g_sw_cloud

    ! in-cloud liquid and ice water path in a layer, in kg m-2
    real(jprb), dimension(nlev,istartcol:iendcol) :: lwp_kg_m2, iwp_kg_m2

    integer  :: jlev, jcol

    ! factor to convert from gridbox-mean mass mixing ratio to
    ! in-cloud water path
    real(jprb) :: factor

    ! convert cloud mixing ratio into liquid and ice water path
    ! in each layer
    do jlev = 1, nlev
      do jcol = istartcol, iendcol
        ! factor to convert from gridbox-mean mass mixing ratio to
        ! in-cloud water path involves the pressure difference in
        ! pa, acceleration due to gravity and cloud fraction
        ! adjusted to avoid division by zero.
        factor = ( thermodynamics%pressure_hl(jcol,jlev+1)    &
             &    -thermodynamics%pressure_hl(jcol,jlev  )  ) &
             &   / (accelduetogravity &
             &   * max(epsilon(1.0_jprb), cloud%fraction(jcol,jlev)))
        lwp_kg_m2(jlev,jcol) = factor * cloud%q_liq(jcol,jlev)
        iwp_kg_m2(jlev,jcol) = factor * cloud%q_ice(jcol,jlev)
      end do
    end do

    ! geometric optics approximation: particles treated as much larger
    ! than the wavelength in both shortwave and longwave
    od_sw_cloud(1,:,:) &
         &   = (3.0_jprb/(2.0_jprb*densityliquidwater)) &
         &   * lwp_kg_m2 / transpose(cloud%re_liq(istartcol:iendcol,:)) &
         &   + (3.0_jprb / (2.0_jprb * densitysolidice)) &
         &   * iwp_kg_m2 / transpose(cloud%re_ice(istartcol:iendcol,:))
    od_lw_cloud(1,:,:) = lwp_kg_m2 * 137.22_jprb &
         &   + (3.0_jprb / (2.0_jprb * densitysolidice)) &
         &   * iwp_kg_m2 / transpose(cloud%re_ice(istartcol:iendcol,:))

    if (config%iverbose >= 4) then
      do jcol = istartcol,iendcol
        write(*,'(a,i0,a,f7.3,a,f7.3)') 'profile ', jcol, ': shortwave optical depth = ', &
             &  sum(od_sw_cloud(1,:,jcol)*cloud%fraction(jcol,:)), &
             &  ', longwave optical depth = ', &
             &  sum(od_lw_cloud(1,:,jcol)*cloud%fraction(jcol,:))
        !    print *, 'lwp = ', sum(lwp_kg_m2(:,istartcol)*cloud%fraction(istartcol,:))
      end do
    end if

    ssa_sw_cloud = config%mono_sw_single_scattering_albedo
    g_sw_cloud   = config%mono_sw_asymmetry_factor

    ! in-place delta-eddington scaling
    call delta_eddington(od_sw_cloud, ssa_sw_cloud, g_sw_cloud)

    if (config%do_lw_cloud_scattering) then
      ssa_lw_cloud = config%mono_lw_single_scattering_albedo
      g_lw_cloud   = config%mono_lw_asymmetry_factor
      ! in-place delta-eddington scaling
      call delta_eddington(od_lw_cloud, ssa_lw_cloud, g_lw_cloud)
    end if

  end subroutine cloud_optics


  !---------------------------------------------------------------------
  ! dummy subroutine since no aerosols are represented in
  ! monochromatic case
  subroutine add_aerosol_optics(nlev,istartcol,iendcol, &
       &  config, thermodynamics, gas, aerosol, & 
       &  od_lw, ssa_lw, g_lw, od_sw, ssa_sw, g_sw)

    use parkind1,                      only : jprb

    use radiation_config,              only : config_type
    use radiation_thermodynamics,      only : thermodynamics_type
    use radiation_gas,                 only : gas_type
    use radiation_aerosol,             only : aerosol_type

    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in), target :: config
    type(thermodynamics_type),intent(in)  :: thermodynamics
    type(gas_type),           intent(in)  :: gas
    type(aerosol_type),       intent(in)  :: aerosol
    ! optical depth, single scattering albedo and asymmetry factor of
    ! the atmosphere (gases on input, gases and aerosols on output)
    ! for each g point. note that longwave ssa and asymmetry and
    ! shortwave asymmetry are all zero for gases, so are not yet
    ! defined on input and are therefore intent(out).
    real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), intent(inout) :: od_lw
    real(jprb), dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol), &
         &  intent(out) :: ssa_lw, g_lw
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(inout) &
         &  :: od_sw, ssa_sw
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(out) :: g_sw

    g_sw(:,:,istartcol:iendcol) = 0.0_jprb

    if (config%do_lw_aerosol_scattering) then
      ssa_lw(:,:,istartcol:iendcol) = 0.0_jprb
      g_lw(:,:,istartcol:iendcol)   = 0.0_jprb
    end if

  end subroutine add_aerosol_optics

  !---------------------------------------------------------------------
  ! planck function in terms of wavelength
  elemental function planck_function(wavelength, temperature)

    use parkind1,            only : jprb

    use radiation_constants, only : boltzmannconstant, planckconstant, &
         &                          speedoflight

    real(jprb), intent(in) :: wavelength  ! metres
    real(jprb), intent(in) :: temperature ! kelvin

    ! output in w sr-1 m-3
    real(jprb)             :: planck_function

    if (temperature > 0.0_jprb) then
      planck_function = 2.0_jprb * planckconstant * speedoflight**2 &
           &   / (wavelength**5 &
           &   * (exp(planckconstant*speedoflight &
           &         /(wavelength*boltzmannconstant*temperature)) - 1.0_jprb))
    else
      planck_function = 0.0_jprb
    end if

  end function planck_function

end module radiation_monochromatic
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

