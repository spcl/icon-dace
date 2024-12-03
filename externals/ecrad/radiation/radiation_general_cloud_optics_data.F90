! # 1 "radiation/radiation_general_cloud_optics_data.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_general_cloud_optics_data.f90"
! radiation_general_cloud_optics_data.f90 - type to store generalized cloud optical properties
!
! (c) copyright 2019- ecmwf.
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


! # 1 "radiation/ecrad_config.h" 1
! ecrad_config.h - preprocessor definitions to configure compilation ecrad -*- f90 -*-
!
! (c) copyright 2023- ecmwf.
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
! this file should be included in fortran source files that require
! different optimizations or settings for different architectures and
! platforms.  feel free to maintain a site-specific version of it.

! the following settings turn on optimizations specific to the
! long-vector nec sx (the short-vector x86-64 architecture is assumed
! otherwise). 




  



  




! in the ifs, an mpi version of easy_netcdf capability is used so that
! only one mpi task reads the data files and shares with the other
! tasks. the mpi version is not used for writing files.

!#define easy_netcdf_read_mpi 1
! # 18 "radiation/radiation_general_cloud_optics_data.f90" 2

module radiation_general_cloud_optics_data

  use parkind1, only : jprb

  implicit none

  public

  !---------------------------------------------------------------------
  ! this type holds the configuration information to compute optical
  ! properties for a particular type of cloud or hydrometeor in one of
  ! the shortwave or longwave
  type general_cloud_optics_type
    ! band-specific (or g-point-specific) values as a look-up table
    ! versus effective radius dimensioned (nband,n_effective_radius)

    ! extinction coefficient per unit mass (m2 kg-1)
    real(jprb), allocatable, dimension(:,:) :: &
         &  mass_ext

    ! single-scattering albedo and asymmetry factor (dimensionless)
    real(jprb), allocatable, dimension(:,:) :: &
         &  ssa, asymmetry

    ! number of effective radius coefficients, start value and
    ! interval in look-up table
    integer    :: n_effective_radius = 0
    real(jprb) :: effective_radius_0, d_effective_radius

    ! name of cloud/precip type and scattering model
    ! (e.g. "mie_droplet", "fu-muskatel_ice"). these are used to
    ! generate the name of the data file from which the coefficients
    ! are read.
    character(len=511) :: type_name

    ! do we use bands or g-points?
    logical :: use_bands = .false.

   contains
     procedure :: setup => setup_general_cloud_optics
     procedure :: add_optical_properties
     procedure :: save => save_general_cloud_optics_data

  end type general_cloud_optics_type

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
! # 68 "radiation/radiation_general_cloud_optics_data.f90" 2

  !---------------------------------------------------------------------
  ! setup cloud optics coefficients by reading them from a file
  subroutine setup_general_cloud_optics(this, file_name, specdef, &
       &                                use_bands, use_thick_averaging, &
       &                                weighting_temperature, &
       &                                iverbose)

    use ecradhook,                       only : lhook, dr_hook, jphook



    use easy_netcdf,          only : netcdf_file

    use radiation_spectral_definition, only : spectral_definition_type
    use radiation_io,                  only : nulout, nulerr, radiation_abort

    class(general_cloud_optics_type), intent(inout)    :: this
    character(len=*), intent(in)               :: file_name
    type(spectral_definition_type), intent(in) :: specdef
    logical, intent(in), optional              :: use_bands, use_thick_averaging
    real(jprb), intent(in), optional           :: weighting_temperature ! k
    integer, intent(in), optional              :: iverbose

    ! spectral properties read from file, dimensioned (wavenumber,
    ! n_effective_radius)
    real(jprb), dimension(:,:), allocatable :: mass_ext, & ! m2 kg-1
         &                                     ssa, asymmetry

    ! reflectance of an infinitely thick cloud, needed for thick
    ! averaging
    real(jprb), dimension(:,:), allocatable :: ref_inf

    ! coordinate variables from file
    real(jprb), dimension(:), allocatable :: wavenumber       ! cm-1
    real(jprb), dimension(:), allocatable :: effective_radius ! m

    ! matrix mapping optical properties in the file to values per
    ! g-point or band, such that in the thin-averaging case,
    ! this%mass_ext=matmul(mapping,file%mass_ext), so mapping is
    ! dimensioned (ngpoint,nwav)
    real(jprb), dimension(:,:), allocatable :: mapping

    ! the netcdf file containing the coefficients
    type(netcdf_file)  :: file

    real(jprb) :: diff_spread
    integer    :: iverb
    integer    :: nre  ! number of effective radii
    integer    :: nwav ! number of wavenumbers describing cloud

    logical    :: use_bands_local, use_thick_averaging_local

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_general_cloud_optics_data:setup',0,hook_handle)

    ! set local values of optional inputs
    if (present(iverbose)) then
      iverb = iverbose
    else
      iverb = 2
    end if

    if (present(use_bands)) then
      use_bands_local = use_bands
    else
      use_bands_local = .false.
    end if

    if (present(use_thick_averaging)) then
      use_thick_averaging_local = use_thick_averaging
    else
      use_thick_averaging_local = .false.
    end if

    ! open the scattering file and configure the way it is read
    call file%open(trim(file_name), iverbose=iverb)
    !call file%transpose_matrices()

    ! read coordinate variables
    call file%get('wavenumber', wavenumber)
    call file%get('effective_radius', effective_radius)

    ! read the band-specific coefficients
    call file%get('mass_extinction_coefficient', mass_ext)
    call file%get('single_scattering_albedo', ssa)
    call file%get('asymmetry_factor', asymmetry)

    ! close scattering file
    call file%close()

    ! check effective radius is evenly spaced
    nre = size(effective_radius)
    ! fractional range of differences, should be near zero for evenly
    ! spaced data
    diff_spread = (maxval(effective_radius(2:nre)-effective_radius(1:nre-1))  &
         &        -minval(effective_radius(2:nre)-effective_radius(1:nre-1))) &
         &      /  minval(abs(effective_radius(2:nre)-effective_radius(1:nre-1)))
    if (diff_spread > 0.01_jprb) then
      write(nulerr, '(a,a,a)') '*** error: effective_radius in ', &
           &  trim(file_name), ', is not evenly spaced to 1%'
      call radiation_abort('radiation configuration error')
    end if

    ! set up effective radius coordinate variable
    this%n_effective_radius = nre
    this%effective_radius_0 = effective_radius(1)
    this%d_effective_radius = effective_radius(2) - effective_radius(1)

    nwav = size(wavenumber)

    ! define the mapping matrix
    call specdef%calc_mapping(wavenumber, mapping, &
         weighting_temperature=weighting_temperature, use_bands=use_bands)

    ! thick averaging should be performed on delta-eddington scaled
    ! quantities (it makes no difference to thin averaging)
    call delta_eddington(mass_ext, ssa, asymmetry)

    ! thin averaging
    this%mass_ext  = matmul(mapping, mass_ext)
    this%ssa       = matmul(mapping, mass_ext*ssa) / this%mass_ext
    this%asymmetry = matmul(mapping, mass_ext*ssa*asymmetry) / (this%mass_ext*this%ssa)

    if (use_thick_averaging_local) then
      ! thick averaging as described by edwards and slingo (1996),
      ! modifying only the single-scattering albedo
      allocate(ref_inf(nwav, nre))

      ! eqs. 18 and 17 of edwards & slingo (1996)
      ref_inf = sqrt((1.0_jprb - ssa) / (1.0_jprb - ssa*asymmetry))
      ref_inf = (1.0_jprb - ref_inf) / (1.0_jprb + ref_inf)
      ! here the left-hand side is actually the averaged ref_inf
      this%ssa = matmul(mapping, ref_inf)
      ! eq. 19 of edwards and slingo (1996)
      this%ssa = 4.0_jprb * this%ssa / ((1.0_jprb + this%ssa)**2 &
           &  - this%asymmetry * (1.0_jprb - this%ssa)**2)

      deallocate(ref_inf)
    end if

    deallocate(mapping)

    ! revert back to unscaled quantities
    call revert_delta_eddington(this%mass_ext, this%ssa, this%asymmetry)

    if (iverb >= 2) then
      write(nulout,'(a,a)') '  file: ', trim(file_name)
      if (present(weighting_temperature)) then
        write(nulout,'(a,f7.1,a)') '  weighting temperature: ', weighting_temperature, ' k'
      else
        write(nulout,'(a,f7.1,a)') '  weighting temperature: ', specdef%reference_temperature, ' k'
      end if
      if (use_thick_averaging_local) then
        write(nulout,'(a)') '  ssa averaging: optically thick limit'
      else
        write(nulout,'(a)') '  ssa averaging: optically thin limit'
      end if
      if (use_bands_local) then
        write(nulout,'(a,i0,a)') '  spectral discretization: ', specdef%nband, ' bands'
      else
        write(nulout,'(a,i0,a)') '  spectral discretization: ', specdef%ng, ' g-points'
      end if
      write(nulout,'(a,i0,a,f6.1,a,f6.1,a)') '  effective radius look-up: ', nre, ' points in range ', &
           &  effective_radius(1)*1.0e6_jprb, '-', effective_radius(nre)*1.0e6_jprb, ' um'
      write(nulout,'(a,i0,a,i0,a)') '  wavenumber range: ', int(specdef%min_wavenumber()), '-', &
           &  int(specdef%max_wavenumber()), ' cm-1'
    end if

    if (lhook) call dr_hook('radiation_general_cloud_optics_data:setup',1,hook_handle)

  end subroutine setup_general_cloud_optics


  !---------------------------------------------------------------------
  ! add the optical properties of a particular cloud type to the
  ! accumulated optical properties of all cloud types
  subroutine add_optical_properties(this, ng, nlev, ncol, &
       &                            cloud_fraction, &
       &                            water_path, effective_radius, &
       &                            od, scat_od, scat_asymmetry)

    use ecradhook, only : lhook, dr_hook, jphook

    class(general_cloud_optics_type), intent(in) :: this

    ! number of g points, levels and columns
    integer, intent(in) :: ng, nlev, ncol

    ! properties of present cloud type, dimensioned (ncol,nlev)
    real(jprb), intent(in) :: cloud_fraction(:,:)
    real(jprb), intent(in) :: water_path(:,:)       ! kg m-2
    real(jprb), intent(in) :: effective_radius(:,:) ! m

    ! optical properties which are additive per cloud type,
    ! dimensioned (ng,nlev,ncol)
    real(jprb), intent(inout), dimension(ng,nlev,ncol) &
         &  :: od             ! optical depth of layer
    real(jprb), intent(inout), dimension(ng,nlev,ncol), optional &
         &  :: scat_od, &     ! scattering optical depth of layer
         &     scat_asymmetry ! scattering optical depth x asymmetry factor

    real(jprb) :: od_local

    real(jprb) :: re_index, weight1, weight2
    integer :: ire

    integer :: jcol, jlev, jg

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_general_cloud_optics_data:add_optical_properties',0,hook_handle)

    if (present(scat_od)) then
      do jcol = 1,ncol
        do jlev = 1,nlev
          if (cloud_fraction(jcol, jlev) > 0.0_jprb) then
            re_index = max(1.0_jprb, min(1.0_jprb + (effective_radius(jcol,jlev)-this%effective_radius_0) &
                 &              / this%d_effective_radius, this%n_effective_radius-0.0001_jprb))
            ire = int(re_index)
            weight2 = re_index - ire
            weight1 = 1.0_jprb - weight2
            do jg = 1, ng
              od_local = water_path(jcol, jlev) * (weight1*this%mass_ext(jg,ire) &
                 &                                +weight2*this%mass_ext(jg,ire+1))
              od(jg,jlev,jcol) = od(jg,jlev,jcol) + od_local
              od_local = od_local * (weight1*this%ssa(jg,ire) &
                 &                  +weight2*this%ssa(jg,ire+1))
              scat_od(jg,jlev,jcol) = scat_od(jg,jlev,jcol) + od_local
              scat_asymmetry(jg,jlev,jcol) = scat_asymmetry(jg,jlev,jcol) &
                 & + od_local * (weight1*this%asymmetry(jg,ire) &
                 &              +weight2*this%asymmetry(jg,ire+1))
            end do

          end if
        end do
      end do
    else
      ! no scattering: return the absorption optical depth
      do jcol = 1,ncol
        do jlev = 1,nlev
          if (water_path(jcol, jlev) > 0.0_jprb) then
            re_index = max(1.0_jprb, min(1.0_jprb + (effective_radius(jcol,jlev)-this%effective_radius_0) &
                 &              / this%d_effective_radius, this%n_effective_radius-0.0001_jprb))
            ire = int(re_index)
            weight2 = re_index - ire
            weight1 = 1.0_jprb - weight2
            od(:,jlev,jcol) = od(:,jlev,jcol) &
                 &  + water_path(jcol, jlev) * (weight1*this%mass_ext(:,ire) &
                 &                             +weight2*this%mass_ext(:,ire+1)) &
                 &  * (1.0_jprb - (weight1*this%ssa(:,ire)+weight2*this%ssa(:,ire+1)))
          end if
        end do
      end do
    end if

    if (lhook) call dr_hook('radiation_general_cloud_optics_data:add_optical_properties',1,hook_handle)

  end subroutine add_optical_properties


  !---------------------------------------------------------------------
  ! return the planck function (in w m-2 (cm-1)-1) for a given
  ! wavenumber (cm-1) and temperature (k), ensuring double precision
  ! for internal calculation
  elemental function calc_planck_function_wavenumber(wavenumber, temperature)

    use parkind1,            only : jprb, jprd
    use radiation_constants, only : speedoflight, boltzmannconstant, planckconstant

    real(jprb), intent(in) :: wavenumber  ! cm-1
    real(jprb), intent(in) :: temperature ! k
    real(jprb) :: calc_planck_function_wavenumber

    real(jprd) :: freq ! hz
    real(jprd) :: planck_fn_freq ! w m-2 hz-1

    freq = 100.0_jprd * real(speedoflight,jprd) * real(wavenumber,jprd)
    planck_fn_freq = 2.0_jprd * real(planckconstant,jprd) * freq**3 &
         &  / (real(speedoflight,jprd)**2 * (exp(real(planckconstant,jprd)*freq &
         &     /(real(boltzmannconstant,jprd)*real(temperature,jprd))) - 1.0_jprd))
    calc_planck_function_wavenumber = real(planck_fn_freq * 100.0_jprd * real(speedoflight,jprd), jprb)

  end function calc_planck_function_wavenumber

  !---------------------------------------------------------------------
  ! save cloud optical properties in the named file
  subroutine save_general_cloud_optics_data(this, file_name, iverbose)

    use ecradhook,     only : lhook, dr_hook, jphook
    use easy_netcdf, only : netcdf_file

    class(general_cloud_optics_type), intent(in) :: this
    character(len=*),                 intent(in) :: file_name
    integer,                optional, intent(in) :: iverbose

    ! object for output netcdf file
    type(netcdf_file) :: out_file

    real(jprb) :: effective_radius(this%n_effective_radius)
    integer :: ire
    
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_general_cloud_optics_data:save',0,hook_handle)

    ! create the file
    call out_file%create(trim(file_name), iverbose=iverbose)

    ! define dimensions
    call out_file%define_dimension("band", size(this%mass_ext,1))
    call out_file%define_dimension("effective_radius", this%n_effective_radius)

    ! put global attributes
    call out_file%put_global_attributes( &
         &   title_str="optical properties of "//trim(this%type_name) &
         &   //" hydrometeors using the spectral intervals of ecrad", &
         &   source_str="ecrad offline radiation model")

    ! define variables
    call out_file%define_variable("effective_radius", units_str="m", &
         &  long_name="effective radius", dim1_name="effective_radius")
    call out_file%define_variable("mass_extinction_coefficient", units_str="m2 kg-1", &
         &  long_name="mass-extinction coefficient", &
         &  dim2_name="effective_radius", dim1_name="band")
    call out_file%define_variable("single_scattering_albedo", units_str="1", &
         &  long_name="single scattering albedo", &
         &  dim2_name="effective_radius", dim1_name="band")
    call out_file%define_variable("asymmetry_factor", units_str="1", &
         &  long_name="asymmetry factor", &
         &  dim2_name="effective_radius", dim1_name="band")

    ! define effective radius
    do ire = 1,this%n_effective_radius
      effective_radius(ire) = this%effective_radius_0 + this%d_effective_radius*(ire-1)
    end do
    
    ! write variables
    call out_file%put("effective_radius", effective_radius)
    call out_file%put("mass_extinction_coefficient", this%mass_ext)
    call out_file%put("single_scattering_albedo", this%ssa)
    call out_file%put("asymmetry_factor", this%asymmetry)
    
    call out_file%close()

    if (lhook) call dr_hook('radiation_general_cloud_optics_data:save',1,hook_handle)

  end subroutine save_general_cloud_optics_data
  
end module radiation_general_cloud_optics_data
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

