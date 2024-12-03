! # 1 "radiation/radiation_aerosol_optics_data.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_aerosol_optics_data.f90"
! radiation_aerosol_optics_data.f90 - type to store aerosol optical properties
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
! modifications
!   2017-10-23  r. hogan  renamed single-character variables
!   2018-04-20  a. bozzo  read optical properties at selected wavelengths


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
! # 20 "radiation/radiation_aerosol_optics_data.f90" 2

module radiation_aerosol_optics_data

  use parkind1,      only : jprb
  use radiation_io,  only : nulerr, radiation_abort

  implicit none
  public

  private :: get_line

  ! the following provide possible values for
  ! aerosol_optics_type%iclass, which is used to map the user's type
  ! index to the hydrophobic or hydrophilic types loaded from the
  ! aerosol optics file. initially iclass is equal to
  ! aerosolclassundefined, which will throw an error if ever the user
  ! tries to use this aerosol type. the user may specify that an
  ! aerosol type is to be ignored in the radiation calculation, in
  ! which case iclass will be set equal to aerosolclassignored.
  enum, bind(c)
     enumerator iaerosolclassundefined,   iaerosolclassignored, &
          &     iaerosolclasshydrophobic, iaerosolclasshydrophilic
  end enum

  integer, parameter :: nmaxstringlength = 2000
  integer, parameter :: nmaxlinelength   = 200

  !---------------------------------------------------------------------
  ! this type holds the configuration information to compute
  ! aerosol optical properties
  type aerosol_optics_type
     ! a vector of length ntype, iclass maps user-defined types on to
     ! the hydrophilic or hydrophobic aerosol classes using the
     ! enumerators above
     integer, allocatable, dimension(:) :: iclass

     ! also a vector of length ntype, itype maps user-defined types on
     ! to specific hydrophilic or hydrophobic aerosol types
     integer, allocatable, dimension(:) :: itype

     ! wavenumber (cm-1) upper and lower bounds of each spectral
     ! interval, which if used in the rrtmg gas optics scheme should
     ! match its band bounds
     real(jprb), allocatable, dimension(:) :: wavenumber1_sw, wavenumber2_sw
     real(jprb), allocatable, dimension(:) :: wavenumber1_lw, wavenumber2_lw

     ! scattering properties are provided separately in the shortwave
     ! and longwave for hydrophobic and hydrophilic aerosols.
     ! hydrophobic aerosols are dimensioned (nband,n_type_phobic):
     real(jprb), allocatable, dimension(:,:) :: &
          &  mass_ext_sw_phobic, & ! mass-extinction coefficient (m2 kg-1)
          &  ssa_sw_phobic,      & ! single scattering albedo
          &  g_sw_phobic,        & ! asymmetry factor
!          &  ssa_g_sw_phobic,    & ! ssa*g
          &  mass_ext_lw_phobic, & ! mass-extinction coefficient (m2 kg-1)
!          &  mass_abs_lw_phobic, & ! mass-absorption coefficient (m2 kg-1)
          &  ssa_lw_phobic,      & ! single scattering albedo
          &  g_lw_phobic           ! asymmetry factor

     ! hydrophilic aerosols are dimensioned (nband,nrh,n_type_philic):
     real(jprb), allocatable, dimension(:,:,:) :: &
          &  mass_ext_sw_philic, & ! mass-extinction coefficient (m2 kg-1)
          &  ssa_sw_philic,      & ! single scattering albedo
          &  g_sw_philic,        & ! asymmetry factor
 !         &  ssa_g_sw_philic,    & ! ssa*g
          &  mass_ext_lw_philic, & ! mass-extinction coefficient (m2 kg-1)
 !         &  mass_abs_lw_philic, & ! mass-absorption coefficient (m2 kg-1)
          &  ssa_lw_philic,      & ! single scattering albedo
          &  g_lw_philic           ! asymmetry factor

     ! wavelengths at which monochromatic properties are stored,
     ! dimensioned (n_mono_wl), units metres
     real(jprb), allocatable :: wavelength_mono(:)

     ! scattering properties at selected monochromatic wavelengths
     ! (n_mono_wl,n_type_phobic)
     real(jprb), allocatable, dimension(:,:) :: &
          &  mass_ext_mono_phobic, & ! mass-extinction coefficient (m2 kg-1)
          &  ssa_mono_phobic,      & ! single scattering albedo
          &  g_mono_phobic,        & ! asymmetry factor
          &  lidar_ratio_mono_phobic ! lidar ratio
     ! ...hydrophilic aerosols dimensioned (n_mono_wl,nrh,n_type_philic):
     real(jprb), allocatable, dimension(:,:,:) :: &
          &  mass_ext_mono_philic, & ! mass-extinction coefficient (m2 kg-1)
          &  ssa_mono_philic,      & ! single scattering albedo
          &  g_mono_philic,        & ! asymmetry factor
          &  lidar_ratio_mono_philic ! lidar ratio

     ! for hydrophilic aerosols, the lower bounds of the relative
     ! humidity bins is a vector of length nrh:
     real(jprb), allocatable, dimension(:) :: &
          &  rh_lower    ! dimensionless (1.0 = 100% humidity)

     ! strings describing the aerosol types
     character(len=nmaxstringlength) :: description_phobic_str = ' '
     character(len=nmaxstringlength) :: description_philic_str = ' '

     ! the number of user-defined aerosol types
     integer :: ntype

     ! the number of hydrophobic and hydrophilic types read from the
     ! aerosol optics file
     integer :: n_type_phobic = 0
     integer :: n_type_philic = 0

     ! number of relative humidity bins
     integer :: nrh = 0

     ! number of longwave and shortwave bands of the data in the file,
     ! and monochromatic wavelengths
     integer :: n_bands_lw = 0, n_bands_sw = 0, n_mono_wl = 0

     ! do we have any hydrophilic types?
     logical :: use_hydrophilic = .true.

     ! do we have monochromatic optical properties
     logical :: use_monochromatic = .false.

   contains
     procedure :: setup => setup_aerosol_optics
     procedure :: save  => save_aerosol_optics
     procedure :: allocate
     procedure :: initialize_types
     procedure :: set_hydrophobic_type
     procedure :: set_hydrophilic_type
     procedure :: set_empty_type
     procedure :: set_types
     procedure :: calc_rh_index
     procedure :: print_description

  end type aerosol_optics_type

contains


  !---------------------------------------------------------------------
  ! setup aerosol optics coefficients by reading them from a file
  subroutine setup_aerosol_optics(this, file_name, iverbose)

    use ecradhook,              only : lhook, dr_hook, jphook



    use easy_netcdf,          only : netcdf_file

    use radiation_io,         only : nulerr, radiation_abort

    class(aerosol_optics_type), intent(inout) :: this
    character(len=*), intent(in)              :: file_name
    integer, intent(in), optional             :: iverbose

    ! the netcdf file containing the aerosol optics data
    type(netcdf_file)  :: file

    real(jprb), allocatable :: wavelength_tmp(:)
    integer            :: iverb

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics_data:setup',0,hook_handle)

    if (present(iverbose)) then
       iverb = iverbose
    else
       iverb = 2
    end if

    ! open the aerosol scattering file and configure the way it is
    ! read
    call file%open(trim(file_name), iverbose=iverb)

    if (file%exists('mass_ext_sw_hydrophilic')) then
      this%use_hydrophilic = .true.
    else
      this%use_hydrophilic = .false.
    end if

    ! read the wavenumber bounds
    call file%get('wavenumber1_sw', this%wavenumber1_sw)
    call file%get('wavenumber2_sw', this%wavenumber2_sw)
    call file%get('wavenumber1_lw', this%wavenumber1_lw)
    call file%get('wavenumber2_lw', this%wavenumber2_lw)

    ! read the raw scattering data
    call file%get('mass_ext_sw_hydrophobic', this%mass_ext_sw_phobic)
    call file%get('ssa_sw_hydrophobic',      this%ssa_sw_phobic)
    call file%get('asymmetry_sw_hydrophobic',this%g_sw_phobic)
    call file%get('mass_ext_lw_hydrophobic', this%mass_ext_lw_phobic)
    call file%get('ssa_lw_hydrophobic',      this%ssa_lw_phobic)
    call file%get('asymmetry_lw_hydrophobic',this%g_lw_phobic)

    call file%get_global_attribute('description_hydrophobic', &
         &                         this%description_phobic_str)

    ! precompute ssa*g for the shortwave and mass-absorption for the
    ! longwave - tbd fix
    !allocate(this%ssa_g_sw_phobic(...

    if (this%use_hydrophilic) then
      call file%get('mass_ext_sw_hydrophilic', this%mass_ext_sw_philic)
      call file%get('ssa_sw_hydrophilic',      this%ssa_sw_philic)
      call file%get('asymmetry_sw_hydrophilic',this%g_sw_philic)
      call file%get('mass_ext_lw_hydrophilic', this%mass_ext_lw_philic)
      call file%get('ssa_lw_hydrophilic',      this%ssa_lw_philic)
      call file%get('asymmetry_lw_hydrophilic',this%g_lw_philic)

      call file%get('relative_humidity1',      this%rh_lower)

      call file%get_global_attribute('description_hydrophilic', &
           &                         this%description_philic_str)
    end if

    ! read the monochromatic scattering data at selected wavelengths
    ! if available in the input file
    if (file%exists('mass_ext_mono_hydrophobic')) then
      this%use_monochromatic = .true.

      if (allocated(this%wavelength_mono)) then
        ! user has provided required monochromatic wavelengths, which
        ! must match those in the file (in the more recent "general"
        ! aerosol optics, interpolation provides optical properties at
        ! the requested wavelengths)
        call file%get('wavelength_mono', wavelength_tmp)
        if (size(wavelength_tmp) /= size(this%wavelength_mono)) then
          write(nulerr,'(a,i0,a,i0,a)') '*** error: ', size(this%wavelength_mono), &
               &  ' monochromatic wavelengths requested but ', &
               &  size(wavelength_tmp), ' in file'
          call radiation_abort('radiation configuration error')
        end if
        if (any(abs(this%wavelength_mono-wavelength_tmp) &
               &  / this%wavelength_mono > 0.01_jprb)) then
          write(nulerr,'(a,a)') '*** error: requested monochromatic wavelengths', &
               &  'must all be within 1% of values in file'
          call radiation_abort('radiation configuration error')
        end if
      else
        ! user has not provided required wavelengths, so we save the
        ! monochromatic wavelengths in the file
        call file%get('wavelength_mono', this%wavelength_mono)
      end if

      call file%get('mass_ext_mono_hydrophobic', this%mass_ext_mono_phobic)
      call file%get('ssa_mono_hydrophobic',      this%ssa_mono_phobic)
      call file%get('asymmetry_mono_hydrophobic',this%g_mono_phobic)
      call file%get('lidar_ratio_mono_hydrophobic',this%lidar_ratio_mono_phobic)
      if (this%use_hydrophilic) then
        call file%get('mass_ext_mono_hydrophilic', this%mass_ext_mono_philic)
        call file%get('ssa_mono_hydrophilic',      this%ssa_mono_philic)
        call file%get('asymmetry_mono_hydrophilic',this%g_mono_philic)
        call file%get('lidar_ratio_mono_hydrophilic',this%lidar_ratio_mono_philic)
      end if
    else
      this%use_monochromatic = .false.
    end if

    ! close aerosol scattering file
    call file%close()

    ! get array sizes
    this%n_bands_lw    = size(this%mass_ext_lw_phobic, 1)
    this%n_bands_sw    = size(this%mass_ext_sw_phobic, 1)
    if (this%use_monochromatic) then
      this%n_mono_wl   = size(this%mass_ext_mono_phobic, 1)
    else
      this%n_mono_wl   = 0
    end if
    this%n_type_phobic = size(this%mass_ext_lw_phobic, 2)

    if (this%use_hydrophilic) then
      this%n_type_philic = size(this%mass_ext_lw_philic, 3)
      this%nrh           = size(this%mass_ext_lw_philic, 2)

      ! check agreement of dimensions
      if (size(this%mass_ext_lw_philic,1) /= this%n_bands_lw) then
        write(nulerr,'(a,a)') '*** error: mass extinction for hydrophilic and hydrophobic ', &
             &                'aerosol have different numbers of longwave bands'
        call radiation_abort('radiation configuration error')
      end if
      if (size(this%mass_ext_sw_philic,1) /= this%n_bands_sw) then
        write(nulerr,'(a,a)') '*** error: mass extinction for hydrophilic and hydrophobic ', &
             &                'aerosol have different numbers of shortwave bands'
        call radiation_abort('radiation configuration error')
      end if
      if (size(this%rh_lower) /= this%nrh) then
        write(nulerr,'(a)') '*** error: size(relative_humidity1) /= size(mass_ext_sw_hydrophilic,2)'
        call radiation_abort('radiation configuration error')
      end if

    else
      this%n_type_philic = 0
      this%nrh           = 0
    end if

    if (lhook) call dr_hook('radiation_aerosol_optics_data:setup',1,hook_handle)

  end subroutine setup_aerosol_optics


  !---------------------------------------------------------------------
  ! initialize the arrays describing the user's aerosol types
  subroutine initialize_types(this, ntype)

    class(aerosol_optics_type), intent(inout) :: this
    integer,                    intent(in)    :: ntype

    ! allocate memory for mapping arrays
    this%ntype = ntype
    allocate(this%iclass(ntype))
    allocate(this%itype(ntype))

    this%iclass = iaerosolclassundefined
    this%itype  = 0

  end subroutine initialize_types

  !---------------------------------------------------------------------
  ! allocate arrays for aerosol optics data type
  subroutine allocate(this, n_type_phobic, n_type_philic, nrh, &
       &              n_bands_lw, n_bands_sw, n_mono_wl)

    use ecradhook,     only : lhook, dr_hook, jphook

    class(aerosol_optics_type), intent(inout) :: this
    integer, intent(in) :: n_type_phobic, n_type_philic, nrh
    integer, intent(in) :: n_bands_lw, n_bands_sw, n_mono_wl

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics_data:allocate',0,hook_handle)

    this%n_type_phobic = n_type_phobic
    this%n_type_philic = n_type_philic
    this%nrh           = nrh
    this%n_bands_lw    = n_bands_lw
    this%n_bands_sw    = n_bands_sw
    this%n_mono_wl     = n_mono_wl

    if (n_type_philic > 0) then
      this%use_hydrophilic = .true.
    else
      this%use_hydrophilic = .false.
    end if

    if (n_bands_sw > 0) then
      allocate(this%mass_ext_sw_phobic(n_bands_sw, n_type_phobic))
      allocate(this%ssa_sw_phobic(n_bands_sw, n_type_phobic))
      allocate(this%g_sw_phobic(n_bands_sw, n_type_phobic))
    end if
    if (n_bands_lw > 0) then
      allocate(this%mass_ext_lw_phobic(n_bands_lw, n_type_phobic))
      allocate(this%ssa_lw_phobic(n_bands_lw, n_type_phobic))
      allocate(this%g_lw_phobic(n_bands_lw, n_type_phobic))
    end if
    if (n_mono_wl > 0) then
      allocate(this%mass_ext_mono_phobic(n_mono_wl, n_type_phobic))
      allocate(this%ssa_mono_phobic(n_mono_wl, n_type_phobic))
      allocate(this%g_mono_phobic(n_mono_wl, n_type_phobic))
      allocate(this%lidar_ratio_mono_phobic(n_mono_wl, n_type_phobic))
    end if

    if (n_type_philic > 0 .and. nrh > 0) then
      if (n_bands_sw > 0) then
        allocate(this%mass_ext_sw_philic(n_bands_sw, nrh, n_type_philic))
        allocate(this%ssa_sw_philic(n_bands_sw, nrh, n_type_philic))
        allocate(this%g_sw_philic(n_bands_sw, nrh, n_type_philic))
      end if
      if (n_bands_lw > 0) then
        allocate(this%mass_ext_lw_philic(n_bands_lw, nrh, n_type_philic))
        allocate(this%ssa_lw_philic(n_bands_lw, nrh, n_type_philic))
        allocate(this%g_lw_philic(n_bands_lw, nrh, n_type_philic))
      end if
      if (n_mono_wl > 0) then
        allocate(this%mass_ext_mono_philic(n_mono_wl, nrh, n_type_philic))
        allocate(this%ssa_mono_philic(n_mono_wl, nrh, n_type_philic))
        allocate(this%g_mono_philic(n_mono_wl, nrh, n_type_philic))
        allocate(this%lidar_ratio_mono_philic(n_mono_wl, nrh, n_type_philic))
      end if
    end if

    if (lhook) call dr_hook('radiation_aerosol_optics_data:allocate',1,hook_handle)

  end subroutine allocate


  !---------------------------------------------------------------------
  ! save aerosol optical properties in the named file
  subroutine save_aerosol_optics(this, file_name, iverbose)

    use ecradhook,              only : lhook, dr_hook, jphook
    use easy_netcdf,          only : netcdf_file

    class(aerosol_optics_type), intent(inout) :: this
    character(len=*),           intent(in)    :: file_name
    integer,          optional, intent(in)    :: iverbose

    ! object for output netcdf file
    type(netcdf_file) :: out_file

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics_data:save',0,hook_handle)

    ! create the file
    call out_file%create(trim(file_name), iverbose=iverbose)

    ! define dimensions
    call out_file%define_dimension("band_lw", this%n_bands_lw)
    call out_file%define_dimension("band_sw", this%n_bands_sw)
    call out_file%define_dimension("hydrophilic", this%n_type_philic)
    call out_file%define_dimension("hydrophobic", this%n_type_phobic)
    call out_file%define_dimension("relative_humidity", this%nrh)
    !if (this%use_monochromatic) then
    !  call out_file%define_dimension("wavelength_mono", this%n_mono_wl)
    !end if

    ! put global attributes
    call out_file%put_global_attributes( &
         &   title_str="aerosol optical properties in the spectral intervals of the gas-optics scheme for ecrad", &
         &   source_str="ecrad offline radiation model")
    call out_file%put_global_attribute( &
         &  "description_hydrophobic", this%description_phobic_str)
    call out_file%put_global_attribute( &
         &  "description_hydrophilic", this%description_philic_str)

    ! define variables
    call out_file%define_variable("mass_ext_sw_hydrophobic", units_str="m2 kg-1", &
         &  long_name="shortwave mass-extinction coefficient of hydrophobic aerosols", &
         &  dim2_name="hydrophobic", dim1_name="band_sw")
    call out_file%define_variable("ssa_sw_hydrophobic", units_str="1", &
         &  long_name="shortwave single scattering albedo of hydrophobic aerosols", &
         &  dim2_name="hydrophobic", dim1_name="band_sw")
    call out_file%define_variable("asymmetry_sw_hydrophobic", units_str="1", &
         &  long_name="shortwave asymmetry factor of hydrophobic aerosols", &
         &  dim2_name="hydrophobic", dim1_name="band_sw")

    call out_file%define_variable("mass_ext_lw_hydrophobic", units_str="m2 kg-1", &
         &  long_name="longwave mass-extinction coefficient of hydrophobic aerosols", &
         &  dim2_name="hydrophobic", dim1_name="band_lw")
    call out_file%define_variable("ssa_lw_hydrophobic", units_str="1", &
         &  long_name="longwave single scattering albedo of hydrophobic aerosols", &
         &  dim2_name="hydrophobic", dim1_name="band_lw")
    call out_file%define_variable("asymmetry_lw_hydrophobic", units_str="1", &
         &  long_name="longwave asymmetry factor of hydrophobic aerosols", &
         &  dim2_name="hydrophobic", dim1_name="band_lw")

    call out_file%define_variable("mass_ext_sw_hydrophilic", units_str="m2 kg-1", &
         &  long_name="shortwave mass-extinction coefficient of hydrophilic aerosols", &
         &  dim3_name="hydrophilic", dim2_name="relative_humidity", dim1_name="band_sw")
    call out_file%define_variable("ssa_sw_hydrophilic", units_str="1", &
         &  long_name="shortwave single scattering albedo of hydrophilic aerosols", &
         &  dim3_name="hydrophilic", dim2_name="relative_humidity", dim1_name="band_sw")
    call out_file%define_variable("asymmetry_sw_hydrophilic", units_str="1", &
         &  long_name="shortwave asymmetry factor of hydrophilic aerosols", &
         &  dim3_name="hydrophilic", dim2_name="relative_humidity", dim1_name="band_sw")

    call out_file%define_variable("mass_ext_lw_hydrophilic", units_str="m2 kg-1", &
         &  long_name="longwave mass-extinction coefficient of hydrophilic aerosols", &
         &  dim3_name="hydrophilic", dim2_name="relative_humidity", dim1_name="band_lw")
    call out_file%define_variable("ssa_lw_hydrophilic", units_str="1", &
         &  long_name="longwave single scattering albedo of hydrophilic aerosols", &
         &  dim3_name="hydrophilic", dim2_name="relative_humidity", dim1_name="band_lw")
    call out_file%define_variable("asymmetry_lw_hydrophilic", units_str="1", &
         &  long_name="longwave asymmetry factor of hydrophilic aerosols", &
         &  dim3_name="hydrophilic", dim2_name="relative_humidity", dim1_name="band_lw")

    ! write variables
    call out_file%put("mass_ext_sw_hydrophobic", this%mass_ext_sw_phobic)
    call out_file%put("ssa_sw_hydrophobic", this%ssa_sw_phobic)
    call out_file%put("asymmetry_sw_hydrophobic", this%g_sw_phobic)
    call out_file%put("mass_ext_lw_hydrophobic", this%mass_ext_lw_phobic)
    call out_file%put("ssa_lw_hydrophobic", this%ssa_lw_phobic)
    call out_file%put("asymmetry_lw_hydrophobic", this%g_lw_phobic)
    call out_file%put("mass_ext_sw_hydrophilic", this%mass_ext_sw_philic)
    call out_file%put("ssa_sw_hydrophilic", this%ssa_sw_philic)
    call out_file%put("asymmetry_sw_hydrophilic", this%g_sw_philic)
    call out_file%put("mass_ext_lw_hydrophilic", this%mass_ext_lw_philic)
    call out_file%put("ssa_lw_hydrophilic", this%ssa_lw_philic)
    call out_file%put("asymmetry_lw_hydrophilic", this%g_lw_philic)

    call out_file%close()

    if (lhook) call dr_hook('radiation_aerosol_optics_data:save',1,hook_handle)

  end subroutine save_aerosol_optics


  !---------------------------------------------------------------------
  ! map user type "itype" onto stored hydrophobic type "i_type_phobic"
  subroutine set_hydrophobic_type(this, itype, i_type_phobic)

    use ecradhook,     only : lhook, dr_hook, jphook

    class(aerosol_optics_type), intent(inout) :: this
    integer, intent(in)                       :: itype, i_type_phobic
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics_data:set_hydrophobic_type',0,hook_handle)

    if (itype < 1 .or. itype > this%ntype) then
      write(nulerr,'(a,i0)') '*** error: aerosol type must be in the range 1 to ', &
           &                  this%ntype
       call radiation_abort('error setting up aerosols')
    end if
    if (i_type_phobic < 1 .or. i_type_phobic > this%n_type_phobic) then
      write(nulerr,'(a,i0)') '*** error: hydrophobic type must be in the range 1 to ', &
           &                  this%n_type_phobic
      call radiation_abort('error setting up aerosols')
    end if

    this%iclass(itype) = iaerosolclasshydrophobic
    this%itype (itype) = i_type_phobic

    if (lhook) call dr_hook('radiation_aerosol_optics_data:set_hydrophobic_type',1,hook_handle)

  end subroutine set_hydrophobic_type


  !---------------------------------------------------------------------
  ! map user type "itype" onto stored hydrophilic type "i_type_philic"
  subroutine set_hydrophilic_type(this, itype, i_type_philic)

    use ecradhook,     only : lhook, dr_hook, jphook

    class(aerosol_optics_type), intent(inout) :: this
    integer, intent(in)                       :: itype, i_type_philic
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics_data:set_hydrophilic_type',0,hook_handle)

    if (.not. this%use_hydrophilic) then
      write(nulerr,'(a)') '*** error: attempt to set hydrophilic aerosol type when no such types present'
      call radiation_abort('error setting up aerosols')
    end if

    if (itype < 1 .or. itype > this%ntype) then
      write(nulerr,'(a,i0)') '*** error: aerosol type must be in the range 1 to ', &
            &          this%ntype
      call radiation_abort('error setting up aerosols')
    end if
    if (i_type_philic < 1 .or. i_type_philic > this%n_type_philic) then
      write(nulerr,'(a,i0)') '*** error: hydrophilic type must be in the range 1 to ', &
           &                 this%n_type_philic
      call radiation_abort('error setting up aerosols')
    end if

    this%iclass(itype) = iaerosolclasshydrophilic
    this%itype (itype) = i_type_philic

    if (lhook) call dr_hook('radiation_aerosol_optics_data:set_hydrophilic_type',1,hook_handle)

  end subroutine set_hydrophilic_type


  !---------------------------------------------------------------------
  ! set a user type "itype" to be ignored in the radiation scheme
  subroutine set_empty_type(this, itype)

    use ecradhook,     only : lhook, dr_hook, jphook

    class(aerosol_optics_type), intent(inout) :: this
    integer, intent(in)                       :: itype
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics_data:set_empty_type',0,hook_handle)

    if (itype < 1 .or. itype > this%ntype) then
      write(nulerr,'(a,i0)') '*** error: aerosol type must be in the range 1 to ', &
           &                 this%ntype
      call radiation_abort('error setting up aerosols')
    end if

    this%iclass(itype) = iaerosolclassignored

    if (lhook) call dr_hook('radiation_aerosol_optics_data:set_empty_type',1,hook_handle)

  end subroutine set_empty_type


  !---------------------------------------------------------------------
  ! set user types "itypes" to map onto the stored hydrophobic and
  ! hydrophilic types according to its sign and value, with a value of
  ! 0 indicating that this type is to be ignored.  thus if itypes=(/
  ! 3, 4, -6, 0 /) then user types 1 and 2 map on to hydrophobic types
  ! 3 and 4, user type 3 maps on to hydrophilic type 6 and user type 4
  ! is ignored.
  subroutine set_types(this, itypes)

    use ecradhook,     only : lhook, dr_hook, jphook

    class(aerosol_optics_type), intent(inout) :: this
    integer, dimension(:), intent(in)         :: itypes

    integer :: jtype
    integer :: istart, iend
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics_data:set_types',0,hook_handle)

    istart = lbound(itypes,1)
    iend   = ubound(itypes,1)

    do jtype = istart, iend
      if (itypes(jtype) == 0) then
        call this%set_empty_type(jtype)
      else if (itypes(jtype) > 0) then
        call this%set_hydrophobic_type(jtype, itypes(jtype))
      else
        call this%set_hydrophilic_type(jtype, -itypes(jtype))
      end if
    end do

    if (lhook) call dr_hook('radiation_aerosol_optics_data:set_types',1,hook_handle)

  end subroutine set_types


  !---------------------------------------------------------------------
  ! return an index to the relative-humdity array, or zero if no
  ! hydrophilic types are present. this function does so little that
  ! it is best to remove the dr hook call.
  function calc_rh_index(this, rh)

    !use ecradhook,     only : lhook, dr_hook, jphook

    class(aerosol_optics_type), intent(in)    :: this
    real(jprb),                 intent(in)    :: rh
    integer                                   :: calc_rh_index
    !real(jphook) :: hook_handle

    !if (lhook) call dr_hook('radiation_aerosol_optics_data:calc_rh_index',0,hook_handle)

    if (.not. this%use_hydrophilic) then
      calc_rh_index = 0
    else if (rh > this%rh_lower(this%nrh)) then
      calc_rh_index = this%nrh
    else
      calc_rh_index = 1
      do while (rh > this%rh_lower(calc_rh_index + 1))
        calc_rh_index = calc_rh_index + 1
      end do
    end if

    !if (lhook) call dr_hook('radiation_aerosol_optics_data:calc_rh_index',1,hook_handle)

  end function calc_rh_index


  !---------------------------------------------------------------------
  ! print a description of the aerosol types to nulout
  subroutine print_description(this, i_type_map)

    use radiation_io, only : nulout

    class(aerosol_optics_type), intent(in) :: this
    integer,                    intent(in) :: i_type_map(:)

    integer :: jtype

    if (size(i_type_map) > 0) then
      write(nulout,'(a)') 'aerosol mapping:'
    else
      write(nulout,'(a)') 'no aerosol types in radiation scheme'
    end if

    do jtype = 1,size(i_type_map)
      if (i_type_map(jtype) > 0) then
        write(nulout,'(i4,a,a)') jtype, ' -> hydrophobic type ', &
             &  trim(get_line(this%description_phobic_str, i_type_map(jtype)))
      else if (i_type_map(jtype) < 0) then
        write(nulout,'(i4,a,a)') jtype, ' -> hydrophilic type ', &
             &  trim(get_line(this%description_philic_str, -i_type_map(jtype)))
      else
        write(nulout,'(i4,a)') jtype, ' is unused'
      end if
    end do

  end subroutine print_description


  !---------------------------------------------------------------------
  ! private helper function for print_description
  pure function get_line(str,iline) result(line_str)
    character(len=*), intent(in)  :: str
    integer,          intent(in)  :: iline
    character(len=nmaxlinelength) :: line_str

    integer :: istart, iend, i_start_new, ioffset, ilength, i_line_current
    logical :: is_fail

    i_line_current = 1
    istart = 1
    iend = len(str)
    is_fail = .false.
    line_str = ' '

    ! find index of first character
    do while (i_line_current < iline)
      i_start_new = scan(str(istart:iend), new_line(' '))
      if (i_start_new == 0) then
        is_fail = .true.
        cycle
      else
        istart = istart + i_start_new
      end if
      i_line_current = i_line_current + 1
    end do

    if (.not. is_fail) then
      ! find index of last character
      ioffset = scan(str(istart:iend), new_line(' '))
      if (ioffset == 0) then
        ilength = len(trim(str(istart:iend)))
      else
        ilength = ioffset - 1
      end if

      if (ilength > nmaxlinelength) then
        ilength = nmaxlinelength
      end if
      iend = istart + ilength - 1

      line_str = str(istart:iend)
    else
      write(line_str,'(i0,a)') iline, ': <unknown>'
    end if

  end function get_line

end module radiation_aerosol_optics_data
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

