! # 1 "radiation/radiation_spectral_definition.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_spectral_definition.f90"
! radiation_spectral_definition.f90 - derived type to describe a spectral definition
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
! # 18 "radiation/radiation_spectral_definition.f90" 2

module radiation_spectral_definition

  use parkind1,    only : jprb

  implicit none

  public

  real(jprb), parameter :: solarreferencetemperature       = 5777.0_jprb ! k
  real(jprb), parameter :: terrestrialreferencetemperature = 273.15_jprb ! k

  !---------------------------------------------------------------------
  ! a derived type describing the contribution of the g points of a
  ! correlated k-distribution gas-optics model from each part of the
  ! spectrum. this is used primarily to map the cloud and aerosol
  ! optical properties on to the gas g points.
  type spectral_definition_type
    
    ! spectral mapping of g points

    ! number of wavenumber intervals
    integer :: nwav = 0
    ! number of k terms / g points
    integer :: ng   = 0
    ! start and end wavenumber (cm-1), dimensioned (nwav)
    real(jprb), allocatable :: wavenumber1(:)
    real(jprb), allocatable :: wavenumber2(:)
    ! fraction of each g point in each wavenumber interval,
    ! dimensioned (nwav, ng)
    real(jprb), allocatable :: gpoint_fraction(:,:)

    ! spectral weighting information for generating mappings to/from
    ! different spectral grids: this can be in terms of a reference
    ! temperature (k) to generate a planck function, or the
    ! solar_spectral_irradiance (w m-2) if available in the gas-optics
    ! file.
    real(jprb) :: reference_temperature = -1.0_jprb
    real(jprb), allocatable :: solar_spectral_irradiance(:)
    
    ! band information

    ! number of bands
    integer :: nband = 0
    ! lower and upper bounds of wavenumber bands (cm-1), dimensioned
    ! (nband)
    real(jprb), allocatable :: wavenumber1_band(:)
    real(jprb), allocatable :: wavenumber2_band(:)
    ! band (one based) to which each g point belongs
    integer,    allocatable :: i_band_number(:)

  contains
    procedure :: read => read_spectral_definition
    procedure :: allocate_bands_only
    procedure :: deallocate
    procedure :: find => find_wavenumber
    procedure :: calc_mapping
    procedure :: calc_mapping_from_bands
    procedure :: calc_mapping_from_wavenumber_bands
    procedure :: print_mapping_from_bands
    procedure :: min_wavenumber
    procedure :: max_wavenumber

  end type spectral_definition_type

contains

  !---------------------------------------------------------------------
  ! read the description of a spectral definition from a netcdf
  ! file of the type used to describe an ecckd model
  subroutine read_spectral_definition(this, file)




    use easy_netcdf,          only : netcdf_file

    use ecradhook,     only : lhook, dr_hook, jphook

    class(spectral_definition_type), intent(inout) :: this
    type(netcdf_file),               intent(inout) :: file

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_spectral_definition:read',0,hook_handle)

    ! read spectral mapping of g points
    call file%get('wavenumber1', this%wavenumber1)
    call file%get('wavenumber2', this%wavenumber2)
    call file%get('gpoint_fraction', this%gpoint_fraction)

    ! read band information
    call file%get('wavenumber1_band', this%wavenumber1_band)
    call file%get('wavenumber2_band', this%wavenumber2_band)
    call file%get('band_number', this%i_band_number)

    ! read spectral weighting information
    if (file%exists('solar_spectral_irradiance')) then
      ! this is on the same grid as wavenumber1,2
      call file%get('solar_spectral_irradiance', &
           &        this%solar_spectral_irradiance)
    end if
    if (file%exists('solar_irradiance')) then
      ! shortwave default temperature
      this%reference_temperature = solarreferencetemperature
    else
      ! longwave reference temperature
      this%reference_temperature = terrestrialreferencetemperature
    end if
    
    ! band number is 0-based: add 1
    this%i_band_number = this%i_band_number + 1

    this%nwav  = size(this%wavenumber1)
    this%ng    = size(this%gpoint_fraction, 2);
    this%nband = size(this%wavenumber1_band)

    if (lhook) call dr_hook('radiation_spectral_definition:read',1,hook_handle)

  end subroutine read_spectral_definition


  !---------------------------------------------------------------------
  ! store a simple band description by copying over the reference
  ! temperature and the lower and upper wavenumbers of each band
  subroutine allocate_bands_only(this, reference_temperature, wavenumber1, wavenumber2)

    use ecradhook,     only : lhook, dr_hook, jphook

    class(spectral_definition_type), intent(inout) :: this
    real(jprb),                      intent(in)    :: reference_temperature    ! k
    real(jprb),        dimension(:), intent(in)    :: wavenumber1, wavenumber2 ! cm-1

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_spectral_definition:allocate_bands_only',0,hook_handle)

    call this%deallocate()

    this%nband = size(wavenumber1)
    allocate(this%wavenumber1_band(this%nband))
    allocate(this%wavenumber2_band(this%nband))
    this%wavenumber1_band = wavenumber1
    this%wavenumber2_band = wavenumber2
    this%reference_temperature = reference_temperature
    
    if (lhook) call dr_hook('radiation_spectral_definition:allocate_bands_only',1,hook_handle)

  end subroutine allocate_bands_only


  !---------------------------------------------------------------------
  ! deallocate memory inside a spectral definition object
  subroutine deallocate(this)

    class(spectral_definition_type), intent(inout) :: this
    
    this%nwav  = 0
    this%ng    = 0
    this%nband = 0
    this%reference_temperature = -1.0_jprb

    if (allocated(this%wavenumber1))      deallocate(this%wavenumber1)
    if (allocated(this%wavenumber2))      deallocate(this%wavenumber2)
    if (allocated(this%wavenumber1_band)) deallocate(this%wavenumber1_band)
    if (allocated(this%wavenumber2_band)) deallocate(this%wavenumber2_band)
    if (allocated(this%gpoint_fraction))  deallocate(this%gpoint_fraction)
    if (allocated(this%i_band_number))    deallocate(this%i_band_number)

  end subroutine deallocate


  !---------------------------------------------------------------------
  ! find the index to the highest wavenumber in the spectral
  ! definition that is lower than or equal to "wavenumber", used for
  ! implementing look-up tables
  pure function find_wavenumber(this, wavenumber)
    class(spectral_definition_type), intent(in) :: this
    real(jprb),                      intent(in) :: wavenumber ! cm-1
    integer                                     :: find_wavenumber

    if (wavenumber < this%wavenumber1(1) .or. wavenumber > this%wavenumber2(this%nwav)) then
      ! wavenumber not present
      find_wavenumber = 0
    else
      find_wavenumber = 1
      do while (wavenumber > this%wavenumber2(find_wavenumber) &
           &    .and. find_wavenumber < this%nwav)
        find_wavenumber = find_wavenumber + 1
      end do
    end if
  end function find_wavenumber


  !---------------------------------------------------------------------
  ! compute a mapping matrix "mapping" that can be used in an
  ! expression y=matmul(mapping,x) where x is a variable containing
  ! optical properties at each input "wavenumber", and y is this
  ! variable mapped on to the spectral intervals in the spectral
  ! definition "this". 
  subroutine calc_mapping(this, wavenumber, mapping, weighting_temperature, use_bands)

    use ecradhook,      only : lhook, dr_hook, jphook
    use radiation_io, only : nulerr, radiation_abort

    class(spectral_definition_type), intent(in)    :: this
    real(jprb),                      intent(in)    :: wavenumber(:) ! cm-1
    real(jprb), allocatable,         intent(inout) :: mapping(:,:)
    real(jprb), optional,            intent(in)    :: weighting_temperature ! k
    logical,    optional,            intent(in)    :: use_bands

    ! spectral weights to apply, same length as wavenumber above
    real(jprb), dimension(:), allocatable :: weight, planck_weight

    ! wavenumbers (cm-1) marking triangle of influence of a cloud
    ! spectral point
    real(jprb) :: wavenum0, wavenum1, wavenum2

    integer    :: nwav ! number of wavenumbers describing cloud

    ! indices to wavenumber intervals in spectral definition structure
    integer    :: isd, isd0, isd1, isd2

    ! wavenumber index
    integer    :: iwav
    
    ! loop indices
    integer    :: jg, jwav, jband

    logical    :: use_bands_local

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_spectral_definition:calc_mapping',0,hook_handle)

    if (present(use_bands)) then
      use_bands_local = use_bands
    else
      use_bands_local = .false.
    end if

    nwav = size(wavenumber)

    if (allocated(mapping)) then
      deallocate(mapping)
    end if
    
    ! define the mapping matrix
    if (use_bands_local) then
      ! cloud properties per band

      allocate(mapping(this%nband, nwav))
      allocate(weight(nwav))

      ! planck weight uses the wavenumbers of the cloud points
      allocate(planck_weight(nwav))
      if (present(weighting_temperature)) then
        if (weighting_temperature > 0.0_jprb) then
          planck_weight = calc_planck_function_wavenumber(wavenumber, &
               &                          weighting_temperature)
        else
          ! legacy mode: unweighted average
          planck_weight = 1.0_jprb
        end if
      else
        planck_weight = calc_planck_function_wavenumber(wavenumber, &
             &                          this%reference_temperature)
      end if

      do jband = 1,this%nband
        weight = 0.0_jprb
        do jwav = 1,nwav
          ! work out wavenumber range for which this cloud wavenumber
          ! will be applicable
          if (wavenumber(jwav) >= this%wavenumber1_band(jband) &
               & .and. wavenumber(jwav) <= this%wavenumber2_band(jband)) then
            if (jwav > 1) then
              wavenum1 = max(this%wavenumber1_band(jband), &
                   &  0.5_jprb*(wavenumber(jwav-1)+wavenumber(jwav)))
            else
              wavenum1 = this%wavenumber1_band(jband)
            end if
            if (jwav < nwav) then
              wavenum2 = min(this%wavenumber2_band(jband), &
                   &  0.5_jprb*(wavenumber(jwav)+wavenumber(jwav+1)))
            else
              wavenum2 = this%wavenumber2_band(jband)
            end if
            ! this cloud wavenumber is weighted by the wavenumber
            ! range of its applicability multiplied by the planck
            ! function at an appropriate temperature
            weight(jwav) = (wavenum2-wavenum1) * planck_weight(jwav)
          end if
        end do
        if (sum(weight) <= 0.0_jprb) then
          ! no cloud wavenumbers lie in the band; interpolate to
          ! central wavenumber of band instead
          if (wavenumber(1) >= this%wavenumber2_band(jband)) then
            ! band is entirely below first cloudy wavenumber
            weight(1) = 1.0_jprb
          else if (wavenumber(nwav) <= this%wavenumber1_band(jband)) then
            ! band is entirely above last cloudy wavenumber
            weight(nwav) = 1.0_jprb
          else
            ! find interpolating points
            iwav = 2
            do while (wavenumber(iwav) < this%wavenumber2_band(jband))
              iwav = iwav+1
            end do
            weight(iwav-1) = planck_weight(iwav-1) * (wavenumber(iwav) &
                 &  - 0.5_jprb*(this%wavenumber2_band(jband)+this%wavenumber1_band(jband)))
            weight(iwav) = planck_weight(iwav) * (-wavenumber(iwav-1) &
                 &  + 0.5_jprb*(this%wavenumber2_band(jband)+this%wavenumber1_band(jband)))
          end if
        end if
        mapping(jband,:) = weight / sum(weight)
      end do

      deallocate(weight)
      deallocate(planck_weight)

    else
      ! cloud properties per g-point

      if (this%ng == 0) then
        write(nulerr,'(a)') '*** error: requested cloud/aerosol mapping per g-point but only available per band'
        call radiation_abort('radiation configuration error')
      end if

      allocate(mapping(this%ng, nwav))
      allocate(weight(this%nwav))
      allocate(planck_weight(this%nwav))

      if (allocated(this%solar_spectral_irradiance)) then
        planck_weight = this%solar_spectral_irradiance
      else
        planck_weight = calc_planck_function_wavenumber(0.5_jprb &
             &  * (this%wavenumber1 + this%wavenumber2), &
             &  this%reference_temperature)
      end if

      mapping = 0.0_jprb
      ! loop over wavenumbers representing cloud
      do jwav = 1,nwav
        ! clear the weights. the weight says for one wavenumber in the
        ! cloud file, what is its fractional contribution to each of
        ! the spectral-definition intervals
        weight = 0.0_jprb

        ! cloud properties are linearly interpolated between each of
        ! the nwav cloud points; therefore, the influence of a
        ! particular cloud point extends as a triangle between
        ! wavenum0 and wavenum2, peaking at wavenum1
        wavenum1 = wavenumber(jwav)
        isd1 = this%find(wavenum1)
        if (isd1 < 1) then
          cycle
        end if
        if (jwav > 1) then
          wavenum0 = wavenumber(jwav-1)

          ! map triangle under (wavenum0,0) to (wavenum1,1) to the
          ! wavenumbers in this%gpoint_fraction
          isd0 = this%find(wavenum0)
          if (isd0 == isd1) then
            ! triangle completely within the range
            ! this%wavenumber1(isd0)-this%wavenumber2(isd0)
            weight(isd0) = 0.5_jprb*(wavenum1-wavenum0) &
                 &       / (this%wavenumber2(isd0)-this%wavenumber1(isd0))
          else
            if (isd0 >= 1) then
              ! left part of triangle
              weight(isd0) = 0.5_jprb * (this%wavenumber2(isd0)-wavenum0)**2 &
                   &       / ((this%wavenumber2(isd0)-this%wavenumber1(isd0)) &
                   &         *(wavenum1-wavenum0))
            end if
            ! right part of triangle (trapezium)
!            weight(isd1) = 0.5_jprb * (wavenum1-this%wavenumber1(isd1)) &
!                 &       * (wavenum1 + this%wavenumber1(isd1) - 2.0_jprb*wavenum0) &
!                 &       / (wavenum1-wavenum0)
            weight(isd1) = 0.5_jprb * (1.0_jprb &
                 &  + (this%wavenumber1(isd1)-wavenum1)/(wavenum1-wavenum0)) &
                 &  * (wavenum1-this%wavenumber1(isd1)) &
                 &  / (this%wavenumber2(isd1)-this%wavenumber1(isd1))
            if (isd1-isd0 > 1) then
              do isd = isd0+1,isd1-1
                ! intermediate trapezia
                weight(isd) = 0.5_jprb * (this%wavenumber1(isd)+this%wavenumber2(isd) &
                     &                    - 2.0_jprb*wavenum0) &
                     &      / (wavenum1-wavenum0)
              end do
            end if
          end if

        else
          ! first cloud wavenumber: all wavenumbers in the spectral
          ! definition below this will use the first one
          if (isd1 >= 1) then
            weight(1:isd1-1) = 1.0_jprb
            weight(isd1) = (wavenum1-this%wavenumber1(isd1)) &
                 &       / (this%wavenumber2(isd1)-this%wavenumber1(isd1))
          end if
        end if

        if (jwav < nwav) then
          wavenum2 = wavenumber(jwav+1)

          ! map triangle under (wavenum1,1) to (wavenum2,0) to the
          ! wavenumbers in this%gpoint_fraction
          isd2 = this%find(wavenum2)

          if (isd1 == isd2) then
            ! triangle completely within the range
            ! this%wavenumber1(isd1)-this%wavenumber2(isd1)
            weight(isd1) = weight(isd1) + 0.5_jprb*(wavenum2-wavenum1) &
                 &       / (this%wavenumber2(isd1)-this%wavenumber1(isd1))
          else
            if (isd2 >= 1 .and. isd2 <= this%nwav) then
              ! right part of triangle
              weight(isd2) = weight(isd2) + 0.5_jprb * (wavenum2-this%wavenumber1(isd2))**2 &
                   &       / ((this%wavenumber2(isd2)-this%wavenumber1(isd2)) &
                   &         *(wavenum2-wavenum1))
            end if
            ! left part of triangle (trapezium)
!            weight(isd1) = weight(isd1) + 0.5_jprb * (this%wavenumber2(isd1)-wavenum1) &
!                 &       * (wavenum1 + this%wavenumber2(isd1) - 2.0_jprb*wavenum2) &
!                 &       / (wavenum2-wavenum1)
            weight(isd1) = weight(isd1) + 0.5_jprb * (1.0_jprb &
                 &  + (wavenum2-this%wavenumber2(isd1)) / (wavenum2-wavenum1)) &
                 &  * (this%wavenumber2(isd1)-wavenum1) &
                 &  / (this%wavenumber2(isd1)-this%wavenumber1(isd1))
            if (isd2-isd1 > 1) then
              do isd = isd1+1,isd2-1
                ! intermediate trapezia
                weight(isd) = weight(isd) + 0.5_jprb * (2.0_jprb*wavenum2 &
                     & - this%wavenumber1(isd) - this%wavenumber2(isd)) &
                     &      / (wavenum2-wavenum1)
              end do
            end if
          end if

        else
          ! last cloud wavenumber: all wavenumbers in the spectral
          ! definition above this will use the last one
          if (isd1 <= this%nwav) then
            weight(isd1+1:this%nwav) = 1.0_jprb
            weight(isd1) = (this%wavenumber2(isd1)-wavenum1) &
                 &       / (this%wavenumber2(isd1)-this%wavenumber1(isd1))
          end if
        end if

        weight = weight * planck_weight

        do jg = 1,this%ng
          mapping(jg, jwav) = sum(weight * this%gpoint_fraction(:,jg))
        end do

      end do

      deallocate(weight)
      deallocate(planck_weight)

      ! normalize mapping matrix
      do jg = 1,this%ng
        mapping(jg,:) = mapping(jg,:) * (1.0_jprb/sum(mapping(jg,:)))
      end do

    end if

    if (lhook) call dr_hook('radiation_spectral_definition:calc_mapping',1,hook_handle)

  end subroutine calc_mapping


  !---------------------------------------------------------------------
  ! under normal operation (if use_fluxes is .false. or not present),
  ! compute a mapping matrix "mapping" that can be used in an
  ! expression y=matmul(mapping^t,x) where x is a variable containing
  ! optical properties in input bands (e.g. albedo in shortwave albedo
  ! bands), and y is this variable mapped on to the spectral intervals
  ! in the spectral definition "this". note that "mapping" is here
  ! transposed from the convention in the calc_mapping routine.  under
  ! the alternative operation (if use_fluxes is present and .true.),
  ! the mapping works in the reverse sense: if y contains fluxes in
  ! each ecrad band or g-point, then x=matmul(mapping,y) would return
  ! fluxes in x averaged to user-supplied "input" bands. in this
  ! version, the bands are described by their wavelength bounds
  ! (wavelength_bound, which must be increasing and exclude the end
  ! points) and the index of the mapping matrix that each band
  ! corresponds to (i_intervals, which has one more element than
  ! wavelength_bound and can have duplicated values if an
  ! albedo/emissivity value is to be associated with more than one
  ! discontinuous ranges of the spectrum).
  subroutine calc_mapping_from_bands(this, &
       &  wavelength_bound, i_intervals, mapping, use_bands, use_fluxes)

    use ecradhook,      only : lhook, dr_hook, jphook
    use radiation_io, only : nulerr, radiation_abort

    class(spectral_definition_type), intent(in)    :: this
    ! monotonically increasing wavelength bounds (m) between
    ! intervals, not including the outer bounds (which are assumed to
    ! be zero and infinity)
    real(jprb),                      intent(in)    :: wavelength_bound(:)
    ! the albedo band indices corresponding to each interval
    integer,                         intent(in)    :: i_intervals(:)
    real(jprb), allocatable,         intent(inout) :: mapping(:,:)
    logical,    optional,            intent(in)    :: use_bands
    logical,    optional,            intent(in)    :: use_fluxes

    ! planck function and central wavenumber of each wavenumber
    ! interval of the spectral definition
    real(jprb) :: planck(this%nwav)         ! w m-2 (cm-1)-1
    real(jprb) :: wavenumber_mid(this%nwav) ! cm-1

    real(jprb), allocatable :: mapping_denom(:,:)

    real(jprb) :: wavenumber1_bound, wavenumber2_bound

    ! to work out weights we sample the planck function at five points
    ! in the interception between an input interval and a band, and
    ! use the trapezium rule
    integer, parameter :: nsample = 5
    integer :: isamp
    real(jprb), dimension(nsample) :: wavenumber_sample, planck_sample
    real(jprb), parameter :: weight_sample(nsample) &
         &        = [0.5_jprb, 1.0_jprb, 1.0_jprb, 1.0_jprb, 0.5_jprb]

    ! index of input value corresponding to each wavenumber interval
    integer :: i_input(this%nwav)

    ! number of albedo/emissivity values that will be provided, some
    ! of which may span discontinuous intervals in wavelength space
    integer :: ninput

    ! number of albedo/emissivity intervals represented, where some
    ! may be grouped to have the same value of albedo/emissivity (an
    ! example is in the thermal infrared where classically the ifs has
    ! ninput=2 and ninterval=3, since only two emissivities are
    ! provided representing (1) the infrared window, and (2) the
    ! intervals to each side of the infrared window.
    integer :: ninterval

    logical    :: use_bands_local, use_fluxes_local

    ! loop indices
    integer    :: jg, jband, jin, jint, jwav

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_spectral_definition:calc_mapping_from_bands',0,hook_handle)

    if (present(use_bands)) then
      use_bands_local = use_bands
    else
      use_bands_local = .false.
    end if

    if (present(use_fluxes)) then
      use_fluxes_local = use_fluxes
    else
      use_fluxes_local = .false.
    end if

    ! count the number of input intervals 
    ninterval = size(i_intervals)
    ninput    = maxval(i_intervals)
    
    if (allocated(mapping)) then
      deallocate(mapping)
    end if
    
    ! check wavelength is monotonically increasing
    if (ninterval > 2) then
      do jint = 2,ninterval-1
        if (wavelength_bound(jint) <= wavelength_bound(jint-1)) then
          write(nulerr, '(a)') '*** error: wavelength bounds must be monotonically increasing'
          call radiation_abort()
        end if
      end do
    end if

    ! define the mapping matrix
    if (use_bands_local) then
      ! require properties per band

      allocate(mapping(ninput, this%nband))
      mapping = 0.0_jprb

      if (use_fluxes_local) then
        allocate(mapping_denom(ninput, this%nband))
        mapping_denom = 0.0_jprb
      end if

      do jband = 1,this%nband
        do jint = 1,ninterval
          if (jint == 1) then
            ! first input interval in wavelength space: lower
            ! wavelength bound is 0 m, so infinity cm-1
            wavenumber2_bound = this%wavenumber2_band(jband)
          else
            wavenumber2_bound = min(this%wavenumber2_band(jband), &
                 &                  0.01_jprb/wavelength_bound(jint-1))
          end if

          if (jint == ninterval) then
            ! final input interval in wavelength space: upper
            ! wavelength bound is infinity m, so 0 cm-1
            wavenumber1_bound = this%wavenumber1_band(jband)
          else
            wavenumber1_bound = max(this%wavenumber1_band(jband), &
                 &                  0.01_jprb/wavelength_bound(jint))

          end if

          if (wavenumber2_bound > wavenumber1_bound) then
            ! current input interval contributes to current band;
            ! compute the weight of the contribution in proportion to
            ! an approximate calculation of the integral of the planck
            ! function over the relevant part of the spectrum
            wavenumber_sample = wavenumber1_bound + [(isamp,isamp=0,nsample-1)] &
                 &  * (wavenumber2_bound-wavenumber1_bound) / real(nsample-1,jprb)
            planck_sample = calc_planck_function_wavenumber(wavenumber_sample, &
                 &                                 this%reference_temperature)
            mapping(i_intervals(jint),jband) = mapping(i_intervals(jint),jband) &
                 &  + sum(planck_sample*weight_sample) * (wavenumber2_bound-wavenumber1_bound)
            if (use_fluxes_local) then
              ! compute an equivalent sample containing the entire ecrad band
              wavenumber_sample = this%wavenumber1_band(jband) + [(isamp,isamp=0,nsample-1)] &
                   &  * (this%wavenumber2_band(jband)-this%wavenumber1_band(jband)) &
                   &  / real(nsample-1,jprb)
              planck_sample = calc_planck_function_wavenumber(wavenumber_sample, &
                   &                                 this%reference_temperature)
              mapping_denom(i_intervals(jint),jband) = mapping_denom(i_intervals(jint),jband) &
                 &  + sum(planck_sample*weight_sample) * (this%wavenumber2_band(jband)-this%wavenumber1_band(jband))
            end if
          end if

        end do
      end do

      if (use_fluxes_local) then
        mapping = mapping / max(1.0e-12_jprb, mapping_denom)
        deallocate(mapping_denom)
      end if

    else
      ! require properties per g-point

      if (this%ng == 0) then
        write(nulerr,'(a)') '*** error: requested surface mapping per g-point but only available per band'
        call radiation_abort('radiation configuration error')
      end if

      allocate(mapping(ninput,this%ng))
      mapping = 0.0_jprb

      wavenumber_mid = 0.5_jprb * (this%wavenumber1 + this%wavenumber2)
      if (allocated(this%solar_spectral_irradiance)) then
        planck = this%solar_spectral_irradiance
      else
        planck = calc_planck_function_wavenumber(wavenumber_mid, &
             &                       this%reference_temperature)
      end if

! # 723 "radiation/radiation_spectral_definition.f90"

      ! loop through all intervals
      do jint = 1,ninterval
        ! loop through the wavenumbers for gpoint_fraction
        do jwav = 1,this%nwav
          if (jint == 1) then
            ! first input interval in wavelength space: lower
            ! wavelength bound is 0 m, so infinity cm-1
            wavenumber2_bound = this%wavenumber2(jwav)
          else
            wavenumber2_bound = min(this%wavenumber2(jwav), &
                 &                  0.01_jprb/wavelength_bound(jint-1))
          end if

          if (jint == ninterval) then
            ! final input interval in wavelength space: upper
            ! wavelength bound is infinity m, so 0 cm-1
            wavenumber1_bound = this%wavenumber1(jwav)
          else
            wavenumber1_bound = max(this%wavenumber1(jwav), &
                 &                  0.01_jprb/wavelength_bound(jint))

          end if

          if (wavenumber2_bound > wavenumber1_bound) then
            ! overlap between input interval and gpoint_fraction
            ! interval: compute the weight of the contribution in
            ! proportion to an approximate calculation of the integral
            ! of the planck function over the relevant part of the
            ! spectrum
            mapping(i_intervals(jint),:) = mapping(i_intervals(jint),:) + this%gpoint_fraction(jwav,:) &
                 &  * (planck(jwav) * (wavenumber2_bound - wavenumber1_bound) &
                 &                  / (this%wavenumber2(jwav)-this%wavenumber1(jwav)))
          end if
        end do
      end do
      if (use_fluxes_local) then
        do jg = 1,this%ng
          mapping(:,jg) = mapping(:,jg) / sum(this%gpoint_fraction(:,jg) * planck)
        end do
      end if


      
    end if

    if (.not. use_fluxes_local) then
      ! normalize mapping matrix
      do jg = 1,size(mapping,dim=2)
        mapping(:,jg) = mapping(:,jg) * (1.0_jprb/sum(mapping(:,jg)))
      end do
    end if

    if (lhook) call dr_hook('radiation_spectral_definition:calc_mapping_from_bands',1,hook_handle)

  end subroutine calc_mapping_from_bands


  !---------------------------------------------------------------------
  ! as calc_mapping_from_bands but in terms of wavenumber bounds from
  ! wavenumber1 to wavenumber2
  subroutine calc_mapping_from_wavenumber_bands(this, &
       &  wavenumber1, wavenumber2, mapping, use_bands, use_fluxes)

    use ecradhook,      only : lhook, dr_hook, jphook

    class(spectral_definition_type), intent(in)    :: this
    real(jprb), intent(in)    :: wavenumber1(:), wavenumber2(:)
    real(jprb), allocatable,         intent(inout) :: mapping(:,:)
    logical,    optional,            intent(in)    :: use_bands
    logical,    optional,            intent(in)    :: use_fluxes

    ! monotonically increasing wavelength bounds (m) between
    ! intervals, not including the outer bounds (which are assumed to
    ! be zero and infinity)
    real(jprb) :: wavelength_bound(size(wavenumber1)-1)
    ! the albedo band indices corresponding to each interval
    integer    :: i_intervals(size(wavenumber1))

    ! lower wavelength bound (m) of each band
    real(jprb) :: wavelength1(size(wavenumber1))

    logical    :: is_band_unassigned(size(wavenumber1))

    ! number of albedo/emissivity intervals represented, where some
    ! may be grouped to have the same value of albedo/emissivity (an
    ! example is in the thermal infrared where classically the ifs has
    ! ninput=2 and ninterval=3, since only two emissivities are
    ! provided representing (1) the infrared window, and (2) the
    ! intervals to each side of the infrared window.
    integer :: ninterval

    ! index to next band in order of increasing wavelength
    integer :: inext

    ! loop indices
    integer :: jint

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_spectral_definition:calc_mapping_from_wavenumber_bands',0,hook_handle)

    wavelength1 = 0.01_jprb / wavenumber2
    ninterval = size(wavelength1)
    
    is_band_unassigned = .true.

    do jint = 1,ninterval
      inext = minloc(wavelength1, dim=1, mask=is_band_unassigned)
      if (jint > 1) then
        wavelength_bound(jint-1) = wavelength1(inext)
      end if
      is_band_unassigned(inext) = .false.
      i_intervals(jint) = inext
    end do

    call calc_mapping_from_bands(this, wavelength_bound, i_intervals, mapping, use_bands, use_fluxes)

    if (lhook) call dr_hook('radiation_spectral_definition:calc_mapping_from_wavenumber_bands',1,hook_handle)

  end subroutine calc_mapping_from_wavenumber_bands


  !---------------------------------------------------------------------
  ! print out the mapping computed by calc_mapping_from_bands
  subroutine print_mapping_from_bands(this, mapping, use_bands)

    use radiation_io, only : nulout

    class(spectral_definition_type), intent(in) :: this
    real(jprb), allocatable,         intent(in) :: mapping(:,:) ! (ninput,nband/ng)
    logical,    optional,            intent(in) :: use_bands

    logical :: use_bands_local

    integer :: nin, nout
    integer :: jin, jout

    if (present(use_bands)) then
      use_bands_local = use_bands
    else
      use_bands_local = .false.
    end if

    nin = size(mapping,1)
    nout = size(mapping,2)

    if (nin <= 1) then
      write(nulout, '(a)') '  all spectral intervals will use the same albedo/emissivity'
    else if (use_bands_local) then
      write(nulout, '(a,i0,a,i0,a)') '  mapping from ', nin, ' values to ', nout, ' bands (wavenumber ranges in cm-1)'
      if (nout <= 40) then
        do jout = 1,nout
          write(nulout,'(i6,a,i6,a)',advance='no') nint(this%wavenumber1_band(jout)), ' to', &
               &                        nint(this%wavenumber2_band(jout)), ':'
          do jin = 1,nin
            write(nulout,'(f5.2)',advance='no') mapping(jin,jout)
          end do
          write(nulout, '()')
        end do
      else
        do jout = 1,30
          write(nulout,'(i6,a,i6,a)',advance='no') nint(this%wavenumber1_band(jout)), ' to', &
               &                        nint(this%wavenumber2_band(jout)), ':'
          do jin = 1,nin
            write(nulout,'(f5.2)',advance='no') mapping(jin,jout)
          end do
          write(nulout, '()')
        end do
        write(nulout,'(a)') '  ...'
        write(nulout,'(i6,a,i6,a)',advance='no') nint(this%wavenumber1_band(nout)), ' to', &
             &                        nint(this%wavenumber2_band(nout)), ':'
        do jin = 1,nin
          write(nulout,'(f5.2)',advance='no') mapping(jin,nout)
        end do
        write(nulout, '()')
      end if
    else
      write(nulout, '(a,i0,a,i0,a)') '  mapping from ', nin, ' values to ', nout, ' g-points'
      if (nout <= 40) then
        do jout = 1,nout
          write(nulout,'(i3,a)',advance='no') jout, ':'
          do jin = 1,nin
            write(nulout,'(f5.2)',advance='no') mapping(jin,jout)
          end do
          write(nulout, '()')
        end do
      else
        do jout = 1,30
          write(nulout,'(i3,a)',advance='no') jout, ':'
          do jin = 1,nin
            write(nulout,'(f5.2)',advance='no') mapping(jin,jout)
          end do
          write(nulout, '()')
        end do
        write(nulout,'(a)') '  ...'
        write(nulout,'(i3,a)',advance='no') nout, ':'
        do jin = 1,nin
          write(nulout,'(f5.2)',advance='no') mapping(jin,nout)
        end do
        write(nulout, '()')
      end if
    end if

  end subroutine print_mapping_from_bands


  !---------------------------------------------------------------------
  ! return the minimum wavenumber of this object in cm-1
  pure function min_wavenumber(this)

    class(spectral_definition_type), intent(in)    :: this
    real(jprb) :: min_wavenumber

    if (this%nwav > 0) then
      min_wavenumber = this%wavenumber1(1)
    else
      min_wavenumber = minval(this%wavenumber1_band)
    end if

  end function min_wavenumber


  !---------------------------------------------------------------------
  ! return the maximum wavenumber of this object in cm-1
  pure function max_wavenumber(this)

    class(spectral_definition_type), intent(in)    :: this
    real(jprb) :: max_wavenumber

    if (this%nwav > 0) then
      max_wavenumber = this%wavenumber1(this%nwav)
    else
      max_wavenumber = maxval(this%wavenumber2_band)
    end if

  end function max_wavenumber


  !---------------------------------------------------------------------
  ! return the planck function (in w m-2 (cm-1)-1) for a given
  ! wavenumber (cm-1) and temperature (k), ensuring double precision
  ! for internal calculation.  if temperature is 0 or less then unity
  ! is returned; since this function is primarily used to weight an
  ! integral by the planck function, a temperature of 0 or less means
  ! no weighting is to be applied.
  elemental function calc_planck_function_wavenumber(wavenumber, temperature)

    use parkind1,            only : jprb, jprd
    use radiation_constants, only : speedoflight, boltzmannconstant, planckconstant

    real(jprb), intent(in) :: wavenumber  ! cm-1
    real(jprb), intent(in) :: temperature ! k
    real(jprb) :: calc_planck_function_wavenumber

    real(jprd) :: freq ! hz
    real(jprd) :: planck_fn_freq ! w m-2 hz-1

    if (temperature > 0.0_jprd) then
      freq = 100.0_jprd * real(speedoflight,jprd) * real(wavenumber,jprd)
      planck_fn_freq = 2.0_jprd * real(planckconstant,jprd) * freq**3 &
           &  / (real(speedoflight,jprd)**2 * (exp(real(planckconstant,jprd)*freq &
           &     /(real(boltzmannconstant,jprd)*real(temperature,jprd))) - 1.0_jprd))
      calc_planck_function_wavenumber = real(planck_fn_freq * 100.0_jprd * real(speedoflight,jprd), jprb)
    else
      calc_planck_function_wavenumber = 1.0_jprb
    end if

  end function calc_planck_function_wavenumber

end module radiation_spectral_definition
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

