! # 1 "radiation/radiation_ecckd.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_ecckd.f90"
! radiation_ecckd.f90 - ecckd generalized gas optics model
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
! # 18 "radiation/radiation_ecckd.f90" 2

module radiation_ecckd

  use parkind1, only : jprb
  use radiation_gas_constants
  use radiation_ecckd_gas
  use radiation_spectral_definition, only : spectral_definition_type

  implicit none

  public

  !---------------------------------------------------------------------
  ! this derived type contains all the data needed to describe a
  ! correlated k-distribution gas optics model created using the ecckd
  ! tool
  type ckd_model_type

    ! gas information

    ! number of gases
    integer :: ngas = 0
    ! array of individual-gas data objects
    type(ckd_gas_type), allocatable :: single_gas(:)
    ! mapping from the "gas codes" in the radiation_gas_constants
    ! module to an index to the single_gas array, where zero means gas
    ! not present (or part of a composite gas)
    integer :: i_gas_mapping(0:nmaxgases)

    ! coordinates of main look-up table for absorption coeffts

    ! number of pressure and temperature points
    integer :: npress = 0
    integer :: ntemp  = 0
    ! natural logarithm of first (lowest) pressure (pa) and increment
    real(jprb) :: log_pressure1, d_log_pressure
    ! first temperature profile (k), dimensioned (npress)
    real(jprb), allocatable :: temperature1(:)
    ! temperature increment (k)
    real(jprb) :: d_temperature

    ! look-up table for planck function

    ! number of entries
    integer :: nplanck = 0
    ! temperature of first element of look-up table and increment (k)
    real(jprb), allocatable :: temperature1_planck
    real(jprb), allocatable :: d_temperature_planck
    ! planck function (black body flux into a horizontal plane) in w
    ! m-2, dimensioned (ng,nplanck)
    real(jprb), allocatable :: planck_function(:,:)

    ! normalized solar irradiance in each g point, dimensioned (ng)
    real(jprb), allocatable :: norm_solar_irradiance(:)

    ! normalized amplitude of variations in the solar irradiance
    ! through the solar cycle in each g point, dimensioned (ng).
    ! since the user always provides the solar irradiance si
    ! integrated across the spectrum, this variable must sum to zero:
    ! this ensures that the solar irradiance in each g-point is
    ! ssi=si*(norm_solar_irradiance +
    ! a*norm_amplitude_solar_irradiance) for any a, where a denotes
    ! the amplitude of deviations from the mean solar spectrum,
    ! typically between -1.0 and 1.0 and provided by
    ! single_level%solar_spectral_multiplier.
    real(jprb), allocatable :: norm_amplitude_solar_irradiance(:)
    
    ! rayleigh molar scattering coefficient in m2 mol-1 in each g
    ! point
    real(jprb), allocatable :: rayleigh_molar_scat(:)

    ! ! spectral mapping of g points

    ! ! number of wavenumber intervals
    ! integer :: nwav = 0

    ! number of k terms / g points
    integer :: ng   = 0

    ! spectral definition describing bands and g points
    type(spectral_definition_type) :: spectral_def

    ! shortwave: true, longwave: false
    logical :: is_sw

  contains

    procedure :: read => read_ckd_model
    procedure :: read_spectral_solar_cycle
! vectorized version of the optical depth look-up performs better on
! nec, but slower on x86



    procedure :: calc_optical_depth => calc_optical_depth_ckd_model

    procedure :: print => print_ckd_model
    procedure :: calc_planck_function
    procedure :: calc_incoming_sw
!    procedure :: deallocate => deallocate_ckd_model

  end type ckd_model_type


contains

  !---------------------------------------------------------------------
  ! read a complete ecckd gas optics model from a netcdf file
  ! "filename"
  subroutine read_ckd_model(this, filename, iverbose)




    use easy_netcdf,          only : netcdf_file

    !use radiation_io, only : nulerr, radiation_abort
    use ecradhook,              only : lhook, dr_hook, jphook

    class(ckd_model_type), intent(inout) :: this
    character(len=*),      intent(in)    :: filename
    integer, optional,     intent(in)    :: iverbose

    type(netcdf_file) :: file

    real(jprb), allocatable :: pressure_lut(:)
    real(jprb), allocatable :: temperature_full(:,:)
    real(jprb), allocatable :: temperature_planck(:)

    character(len=512) :: constituent_id

    integer :: iverbose_local

    ! loop counters
    integer :: jgas, jjgas

    integer :: istart, inext, nchar, i_gas_code

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_ecckd:read_ckd_model',0,hook_handle)

    if (present(iverbose)) then
      iverbose_local = iverbose
    else
      iverbose_local = 3
    end if

    call file%open(trim(filename), iverbose=iverbose_local)

    ! read temperature and pressure coordinate variables
    call file%get('pressure', pressure_lut)
    this%log_pressure1 = log(pressure_lut(1))
    this%npress = size(pressure_lut)
    this%d_log_pressure = log(pressure_lut(2)) - this%log_pressure1
    call file%get('temperature', temperature_full)
    allocate(this%temperature1(this%npress));
    this%temperature1 = temperature_full(:,1)
    this%d_temperature = temperature_full(1,2)-temperature_full(1,1)
    this%ntemp = size(temperature_full,2)
    deallocate(temperature_full)
    
    ! read planck function, or solar irradiance and rayleigh
    ! scattering coefficient
    if (file%exists('solar_irradiance')) then
      this%is_sw = .true.
      call file%get('solar_irradiance', this%norm_solar_irradiance)
      this%norm_solar_irradiance = this%norm_solar_irradiance &
           &  / sum(this%norm_solar_irradiance)
      call file%get('rayleigh_molar_scattering_coeff', &
           &  this%rayleigh_molar_scat)
    else
      this%is_sw = .false.
      call file%get('temperature_planck', temperature_planck)
      this%nplanck = size(temperature_planck)
      this%temperature1_planck = temperature_planck(1)
      this%d_temperature_planck = temperature_planck(2) - temperature_planck(1)
      deallocate(temperature_planck)
      call file%get('planck_function', this%planck_function)
    end if

    ! read the spectral definition information into a separate
    ! derived type
    call this%spectral_def%read(file);
    this%ng = this%spectral_def%ng

    ! read gases
    call file%get('n_gases', this%ngas)
    allocate(this%single_gas(this%ngas))
    call file%get_global_attribute('constituent_id', constituent_id)
    nchar = len(trim(constituent_id))
    istart = 1
    this%i_gas_mapping = 0
    do jgas = 1, this%ngas
      if (jgas < this%ngas) then
        inext = istart + scan(constituent_id(istart:nchar), ' ')
      else
        inext = nchar+2
      end if
      ! find gas code
      i_gas_code = 0
      do jjgas = 1, nmaxgases
        if (constituent_id(istart:inext-2) == trim(gaslowercasename(jjgas))) then
          i_gas_code = jjgas
          exit
        end if
      end do
      ! if (i_gas_code == 0) then
      !   write(nulerr,'(a,a,a)') '*** error: gas "', constituent_id(istart:inext-2), &
      !        & '" not understood'
      !   call radiation_abort('radiation configuration error')
      ! end if
      this%i_gas_mapping(i_gas_code) = jgas;
      call this%single_gas(jgas)%read(file, constituent_id(istart:inext-2), i_gas_code)
      istart = inext
    end do

    call file%close()

    if (lhook) call dr_hook('radiation_ecckd:read_ckd_model',1,hook_handle)

  end subroutine read_ckd_model

  !---------------------------------------------------------------------
  ! print a description of the correlated k-distribution model to the
  ! "nulout" unit
  subroutine print_ckd_model(this)

    use radiation_io, only : nulout
    use radiation_gas_constants

    class(ckd_model_type), intent(in)  :: this

    integer :: jgas
    
    if (this%is_sw) then
      write(nulout,'(a)',advance='no') 'ecckd shortwave gas optics model: '
    else
      write(nulout,'(a)',advance='no') 'ecckd longwave gas optics model: '
    end if

    write(nulout,'(i0,a,i0,a,i0,a,i0,a)') &
         &  nint(this%spectral_def%wavenumber1(1)), '-', &
         &  nint(this%spectral_def%wavenumber2(size(this%spectral_def%wavenumber2))), &
         &  ' cm-1, ', this%ng, ' g-points in ', this%spectral_def%nband, ' bands'
    write(nulout,'(a,i0,a,i0,a,i0,a)') '  look-up table sizes: ', this%npress, ' pressures, ', &
         &  this%ntemp, ' temperatures, ', this%nplanck, ' planck-function entries'
    write(nulout, '(a)') '  gases:'
    do jgas = 1,this%ngas
      if (this%single_gas(jgas)%i_gas_code > 0) then
        write(nulout, '(a,a,a)', advance='no') '    ', &
             &  trim(gasname(this%single_gas(jgas)%i_gas_code)), ': '
      else
        write(nulout, '(a)', advance='no') '    composite of well-mixed background gases: '
      end if
      select case (this%single_gas(jgas)%i_conc_dependence)
        case (iconcdependencenone)
          write(nulout, '(a)') 'no concentration dependence'
        case (iconcdependencelinear)
          write(nulout, '(a)') 'linear concentration dependence'
        case (iconcdependencerelativelinear)
          write(nulout, '(a,e14.6)') 'linear concentration dependence relative to a mole fraction of ', &
               &  this%single_gas(jgas)%reference_mole_frac
        case (iconcdependencelut)
          write(nulout, '(a,i0,a,e14.6,a,e13.6)') 'look-up table with ', this%single_gas(jgas)%n_mole_frac, &
               &  ' log-spaced mole fractions in range ', exp(this%single_gas(jgas)%log_mole_frac1), &
               &  ' to ', exp(this%single_gas(jgas)%log_mole_frac1 &
               &           + this%single_gas(jgas)%n_mole_frac*this%single_gas(jgas)%d_log_mole_frac)
      end select
    end do

  end subroutine print_ckd_model


  !---------------------------------------------------------------------
  ! read the amplitude of the spectral variations associated with the
  ! solar cycle and map to g-points
  subroutine read_spectral_solar_cycle(this, filename, iverbose, use_updated_solar_spectrum)




    use easy_netcdf,          only : netcdf_file

    use radiation_io,         only : nulout, nulerr, radiation_abort
    use ecradhook,              only : lhook, dr_hook, jphook

    ! reference total solar irradiance (w m-2)
    real(jprb), parameter :: referencetsi = 1361.0_jprb
   
    class(ckd_model_type), intent(inout) :: this
    character(len=*),      intent(in)    :: filename
    integer, optional,     intent(in)    :: iverbose
    ! do we update the mean solar spectral irradiance for each g-point
    ! based on the contents of the file?
    logical, optional,     intent(in)    :: use_updated_solar_spectrum

    type(netcdf_file) :: file

    ! solar spectral irradiance, its amplitude and wavenumber
    ! coordinate variable, read from netcdf file
    real(jprb), allocatable :: wavenumber(:) ! cm-1
    real(jprb), allocatable :: ssi(:) ! w m-2 cm
    real(jprb), allocatable :: ssi_amplitude(:) ! w m-2 cm

    ! as above but on the wavenumber grid delimited by
    ! this%wavenumber1 and this%wavenumber2
    real(jprb), allocatable :: ssi_grid(:)
    real(jprb), allocatable :: ssi_amplitude_grid(:)
    real(jprb), allocatable :: wavenumber_grid(:)
    
    ! old normalized solar irradiance in case it gets changed and we
    ! need to report the amplitude of the change
    real(jprb), allocatable :: old_norm_solar_irradiance(:)
    
    real(jprb) :: dwav_grid
    
    ! number of input wavenumbers, and number on ecckd model's grid
    integer :: nwav, nwav_grid
    ! corresponding loop indices
    integer :: jwav, jwav_grid, jg

    integer :: iband
    
    integer :: iverbose_local

    real(jphook) :: hook_handle
    
    if (lhook) call dr_hook('radiation_ecckd:read_spectral_solar_cycle',0,hook_handle)

    if (present(iverbose)) then
      iverbose_local = iverbose
    else
      iverbose_local = 3
    end if

    call file%open(trim(filename), iverbose=iverbose_local)

    call file%get('wavenumber', wavenumber)
    call file%get('mean_solar_spectral_irradiance', ssi)
    call file%get('ssi_solar_cycle_amplitude', ssi_amplitude)

    call file%close()

    nwav = size(wavenumber)
    
    nwav_grid = size(this%spectral_def%wavenumber1)
    allocate(ssi_grid(nwav_grid))
    allocate(ssi_amplitude_grid(nwav_grid))
    allocate(wavenumber_grid(nwav_grid))
    wavenumber_grid = 0.5_jprb * (this%spectral_def%wavenumber1+this%spectral_def%wavenumber2)
    dwav_grid = this%spectral_def%wavenumber2(1)-this%spectral_def%wavenumber1(1)
    
    ssi_grid = 0.0_jprb
    ssi_amplitude_grid = 0.0_jprb

    ! interpolate input ssi to regular wavenumber grid
    do jwav_grid = 1,nwav_grid
      do jwav = 1,nwav-1
        if (wavenumber(jwav) < wavenumber_grid(jwav_grid) &
             &  .and. wavenumber(jwav+1) >= wavenumber_grid(jwav_grid)) then
          ! linear interpolation - this is not perfect
          ssi_grid(jwav_grid) = (ssi(jwav)*(wavenumber(jwav+1)-wavenumber_grid(jwav_grid)) &
               &                +ssi(jwav+1)*(wavenumber_grid(jwav_grid)-wavenumber(jwav))) &
               &                * dwav_grid / (wavenumber(jwav+1)-wavenumber(jwav))
          ssi_amplitude_grid(jwav_grid) = (ssi_amplitude(jwav)*(wavenumber(jwav+1)-wavenumber_grid(jwav_grid)) &
               &                +ssi_amplitude(jwav+1)*(wavenumber_grid(jwav_grid)-wavenumber(jwav))) &
               &                * dwav_grid / (wavenumber(jwav+1)-wavenumber(jwav))
          exit
        end if
      end do
    end do

    ! optionally update the solar irradiances in each g-point, and the
    ! spectral solar irradiance on the wavenumber grid corresponding
    ! to gpoint_fraction
    allocate(old_norm_solar_irradiance(nwav_grid))
    old_norm_solar_irradiance = this%norm_solar_irradiance
    if (present(use_updated_solar_spectrum)) then
      if (use_updated_solar_spectrum) then
        if (.not. allocated(this%spectral_def%solar_spectral_irradiance)) then
          write(nulerr,'(a)') 'cannot use_updated_solar_spectrum unless gas optics model is from ecckd >= 1.4'
          call radiation_abort()
        end if
        this%norm_solar_irradiance = old_norm_solar_irradiance &
             &  * matmul(ssi_grid,this%spectral_def%gpoint_fraction) &
             &  / matmul(this%spectral_def%solar_spectral_irradiance,this%spectral_def%gpoint_fraction)
        this%norm_solar_irradiance = this%norm_solar_irradiance / sum(this%norm_solar_irradiance)
        this%spectral_def%solar_spectral_irradiance = ssi_grid
      end if
    end if
    
    ! map on to g-points
    this%norm_amplitude_solar_irradiance &
         &  = this%norm_solar_irradiance &
         &  * matmul(ssi_amplitude_grid, this%spectral_def%gpoint_fraction) &
         &  / matmul(ssi_grid,this%spectral_def%gpoint_fraction)
    
    ! remove the mean from the solar-cycle fluctuations, since the
    ! user will scale with total solar irradiance
    this%norm_amplitude_solar_irradiance &
         &  = (this%norm_solar_irradiance+this%norm_amplitude_solar_irradiance) &
         &  / sum(this%norm_solar_irradiance+this%norm_amplitude_solar_irradiance) &
         &  - this%norm_solar_irradiance

    ! print the spectral solar irradiance per g point, and solar cycle amplitude 
    if (iverbose_local >= 2) then
      write(nulout,'(a,f6.1,a)') 'g-point, solar irradiance for nominal tsi = ', &
           &  referencetsi, ' w m-2, solar cycle amplitude (at solar maximum), update to original solar irradiance'
      iband = 0
      do jg = 1,this%ng
        if (this%spectral_def%i_band_number(jg) > iband) then
          iband = this%spectral_def%i_band_number(jg)
          write(nulout, '(i2,f10.4,f7.3,a,f8.4,a,i2,a,f7.1,a,f7.1,a)') &
               &  jg, referencetsi*this%norm_solar_irradiance(jg), &
               &  100.0_jprb * this%norm_amplitude_solar_irradiance(jg) &
               &  / this%norm_solar_irradiance(jg), '% ', &
               &  100.0_jprb * (this%norm_solar_irradiance(jg) &
               &  / old_norm_solar_irradiance(jg) - 1.0_jprb), '% band ', iband, ': ', &
               &  this%spectral_def%wavenumber1_band(iband), '-', &
               &  this%spectral_def%wavenumber2_band(iband), ' cm-1'
        else
          write(nulout, '(i2,f10.4,f7.3,a,f8.4,a)') jg, referencetsi*this%norm_solar_irradiance(jg), &
               &  100.0_jprb * this%norm_amplitude_solar_irradiance(jg) &
               &  / this%norm_solar_irradiance(jg), '% ', &
               &  100.0_jprb * (this%norm_solar_irradiance(jg) &
               &  / old_norm_solar_irradiance(jg) - 1.0_jprb), '%'
        end if
      end do
    end if
    
    if (lhook) call dr_hook('radiation_ecckd:read_spectral_solar_cycle',1,hook_handle)
    
  end subroutine read_spectral_solar_cycle
  

  !---------------------------------------------------------------------
  ! compute layerwise optical depth for each g point for ncol columns
  ! at nlev layers
  subroutine calc_optical_depth_ckd_model(this, ncol, nlev, istartcol, iendcol, nmaxgas, &
       &  pressure_hl, temperature_fl, mole_fraction_fl, &
       &  optical_depth_fl, rayleigh_od_fl, opt_concentration_scaling)

    use ecradhook,             only : lhook, dr_hook, jphook
    use radiation_constants, only : accelduetogravity

    ! input variables

    class(ckd_model_type), intent(in), target  :: this
    ! number of columns, levels and input gases
    integer,               intent(in)  :: ncol, nlev, nmaxgas, istartcol, iendcol
    ! pressure at half levels (pa), dimensioned (ncol,nlev+1)
    real(jprb),            intent(in)  :: pressure_hl(ncol,nlev+1)
    ! temperature at full levels (k), dimensioned (ncol,nlev)
    real(jprb),            intent(in)  :: temperature_fl(istartcol:iendcol,nlev)
    ! gas mole fractions at full levels (mol mol-1), dimensioned (ncol,nlev,nmaxgas)
    real(jprb),            intent(in)  :: mole_fraction_fl(ncol,nlev,nmaxgas)
    ! optional concentration scaling of each gas
    real(jprb), optional,  intent(in)  :: opt_concentration_scaling(nmaxgas)
    
    ! output variables

    ! layer absorption optical depth for each g point
    real(jprb),            intent(out) :: optical_depth_fl(this%ng,nlev,istartcol:iendcol)
    ! in the shortwave only, the rayleigh scattering optical depth
    real(jprb),  optional, intent(out) :: rayleigh_od_fl(this%ng,nlev,istartcol:iendcol)

    ! local variables

    real(jprb), pointer :: molar_abs(:,:,:), molar_abs_conc(:,:,:,:)

    ! natural logarithm of pressure at full levels
    real(jprb) :: log_pressure_fl(nlev)

    ! optical depth of single gas at one point in space versus
    ! spectral interval
    !real(jprb) :: od_single_gas(this%ng)

    real(jprb) :: multiplier(nlev), simple_multiplier(nlev), global_multiplier, temperature1
    real(jprb) :: scaling
    real(jprb) :: concentration_scaling(nmaxgas)

    ! indices and weights in temperature, pressure and concentration interpolation
    real(jprb) :: pindex1, tindex1, cindex1
    real(jprb) :: pw1(nlev), pw2(nlev), tw1(nlev), tw2(nlev), cw1(nlev), cw2(nlev)
    integer    :: ip1(nlev), it1(nlev), ic1(nlev)

    ! natural logarithm of mole fraction at one point
    real(jprb) :: log_conc

    ! minimum mole fraction in look-up-table
    real(jprb) :: mole_frac1

    integer :: jcol, jlev, jgas, igascode

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_ecckd:calc_optical_depth',0,hook_handle)

    global_multiplier = 1.0_jprb / (accelduetogravity * 0.001_jprb * airmolarmass)
    concentration_scaling = 1.0_jprb
    if (present(opt_concentration_scaling)) concentration_scaling = opt_concentration_scaling

    do jcol = istartcol,iendcol

      log_pressure_fl = log(0.5_jprb * (pressure_hl(jcol,1:nlev)+pressure_hl(jcol,2:nlev+1)))

      do jlev = 1,nlev
        ! find interpolation points in pressure
        pindex1 = (log_pressure_fl(jlev)-this%log_pressure1) &
             &    / this%d_log_pressure
        pindex1 = 1.0_jprb + max(0.0_jprb, min(pindex1, this%npress-1.0001_jprb))
        ip1(jlev) = int(pindex1)
        pw2(jlev) = pindex1 - ip1(jlev)
        pw1(jlev) = 1.0_jprb - pw2(jlev)

        ! find interpolation points in temperature
        temperature1 = pw1(jlev)*this%temperature1(ip1(jlev)) &
             &       + pw2(jlev)*this%temperature1(ip1(jlev)+1)
        tindex1 = (temperature_fl(jcol,jlev) - temperature1) &
             &    / this%d_temperature
        tindex1 = 1.0_jprb + max(0.0_jprb, min(tindex1, this%ntemp-1.0001_jprb))
        it1(jlev) = int(tindex1)
        tw2(jlev) = tindex1 - it1(jlev)
        tw1(jlev) = 1.0_jprb - tw2(jlev)

        ! concentration multiplier
        simple_multiplier(jlev) = global_multiplier &
             &  * (pressure_hl(jcol,jlev+1) - pressure_hl(jcol,jlev))
      end do

      optical_depth_fl(:,:,jcol) = 0.0_jprb
      
      do jgas = 1,this%ngas

        associate (single_gas => this%single_gas(jgas))
          igascode = this%single_gas(jgas)%i_gas_code
          
          select case (single_gas%i_conc_dependence)
            
          case (iconcdependencelinear)
            molar_abs => this%single_gas(jgas)%molar_abs
            multiplier = simple_multiplier * mole_fraction_fl(jcol,:,igascode)

            multiplier = multiplier * concentration_scaling(igascode)
            
            do jlev = 1,nlev
              optical_depth_fl(:,jlev,jcol) = optical_depth_fl(:,jlev,jcol) &
                   &        + (multiplier(jlev)*tw1(jlev)) * (pw1(jlev) * molar_abs(:,ip1(jlev),it1(jlev)) &
                   &                +pw2(jlev) * molar_abs(:,ip1(jlev)+1,it1(jlev))) &
                   &        + (multiplier(jlev)*tw2(jlev)) * (pw1(jlev) * molar_abs(:,ip1(jlev),it1(jlev)+1) &
                   &                +pw2(jlev) * molar_abs(:,ip1(jlev)+1,it1(jlev)+1))
            end do

          case (iconcdependencerelativelinear)
            molar_abs => this%single_gas(jgas)%molar_abs

            multiplier = simple_multiplier &
                   &  * (mole_fraction_fl(jcol,:,igascode)*concentration_scaling(igascode) &
                   &     - single_gas%reference_mole_frac)
            
            do jlev = 1,nlev
              optical_depth_fl(:,jlev,jcol) = optical_depth_fl(:,jlev,jcol) &
                   &        + (multiplier(jlev)*tw1(jlev)) * (pw1(jlev) * molar_abs(:,ip1(jlev),it1(jlev)) &
                   &                +pw2(jlev) * molar_abs(:,ip1(jlev)+1,it1(jlev))) &
                   &        + (multiplier(jlev)*tw2(jlev)) * (pw1(jlev) * molar_abs(:,ip1(jlev),it1(jlev)+1) &
                   &                +pw2(jlev) * molar_abs(:,ip1(jlev)+1,it1(jlev)+1))
            end do

          case (iconcdependencenone)
            ! composite gases
            molar_abs => this%single_gas(jgas)%molar_abs
            do jlev = 1,nlev
              optical_depth_fl(:,jlev,jcol) = optical_depth_fl(:,jlev,jcol) &
                   &  + (simple_multiplier(jlev)*tw1(jlev)) * (pw1(jlev) * molar_abs(:,ip1(jlev),it1(jlev)) &
                   &                              +pw2(jlev) * molar_abs(:,ip1(jlev)+1,it1(jlev))) &
                   &  + (simple_multiplier(jlev)*tw2(jlev)) * (pw1(jlev) * molar_abs(:,ip1(jlev),it1(jlev)+1) &
                   &                              +pw2(jlev) * molar_abs(:,ip1(jlev)+1,it1(jlev)+1))
            end do

          case (iconcdependencelut)

            scaling = concentration_scaling(igascode)
            
            ! logarithmic interpolation in concentration space
            molar_abs_conc => this%single_gas(jgas)%molar_abs_conc
            mole_frac1 = exp(single_gas%log_mole_frac1)
            do jlev = 1,nlev
              ! take care of mole_fraction == 0
              log_conc = log(max(mole_fraction_fl(jcol,jlev,igascode)*scaling, mole_frac1))
              cindex1  = (log_conc - single_gas%log_mole_frac1) / single_gas%d_log_mole_frac
              cindex1  = 1.0_jprb + max(0.0_jprb, min(cindex1, single_gas%n_mole_frac-1.0001_jprb))
              ic1(jlev) = int(cindex1)
              cw2(jlev) = cindex1 - ic1(jlev)
              cw1(jlev) = 1.0_jprb - cw2(jlev)
            end do
              ! od_single_gas = cw1 * (tw1 * (pw1 * molar_abs_conc(:,ip1,it1,ic1) &
              !      &                       +pw2 * molar_abs_conc(:,ip1+1,it1,ic1)) &
              !      &                +tw2 * (pw1 * molar_abs_conc(:,ip1,it1+1,ic1) &
              !      &                       +pw2 * molar_abs_conc(:,ip1+1,it1+1,ic1))) &
              !      &         +cw2 * (tw1 * (pw1 * molar_abs_conc(:,ip1,it1,ic1+1) &
              !      &                       +pw2 * molar_abs_conc(:,ip1+1,it1,ic1+1)) &
              !      &                +tw2 * (pw1 * molar_abs_conc(:,ip1,it1+1,ic1+1) &
              !      &                       +pw2 * molar_abs_conc(:,ip1+1,it1+1,ic1+1)))
            do jlev = 1,nlev
              optical_depth_fl(:,jlev,jcol) = optical_depth_fl(:,jlev,jcol) &
                   &  + (simple_multiplier(jlev) * mole_fraction_fl(jcol,jlev,igascode) * scaling) * ( &
                   &      (cw1(jlev) * tw1(jlev) * pw1(jlev)) * molar_abs_conc(:,ip1(jlev),it1(jlev),ic1(jlev)) &
                   &     +(cw1(jlev) * tw1(jlev) * pw2(jlev)) * molar_abs_conc(:,ip1(jlev)+1,it1(jlev),ic1(jlev)) &
                   &     +(cw1(jlev) * tw2(jlev) * pw1(jlev)) * molar_abs_conc(:,ip1(jlev),it1(jlev)+1,ic1(jlev)) &
                   &     +(cw1(jlev) * tw2(jlev) * pw2(jlev)) * molar_abs_conc(:,ip1(jlev)+1,it1(jlev)+1,ic1(jlev)) &
                   &     +(cw2(jlev) * tw1(jlev) * pw1(jlev)) * molar_abs_conc(:,ip1(jlev),it1(jlev),ic1(jlev)+1) &
                   &     +(cw2(jlev) * tw1(jlev) * pw2(jlev)) * molar_abs_conc(:,ip1(jlev)+1,it1(jlev),ic1(jlev)+1) &
                   &     +(cw2(jlev) * tw2(jlev) * pw1(jlev)) * molar_abs_conc(:,ip1(jlev),it1(jlev)+1,ic1(jlev)+1) &
                   &     +(cw2(jlev) * tw2(jlev) * pw2(jlev)) * molar_abs_conc(:,ip1(jlev)+1,it1(jlev)+1,ic1(jlev)+1))
            end do
          end select
            
        end associate

      end do

      ! ensure the optical depth is not negative
      optical_depth_fl(:,:,jcol) = max(0.0_jprb, optical_depth_fl(:,:,jcol))

      ! rayleigh scattering
      if (this%is_sw .and. present(rayleigh_od_fl)) then
        do jlev = 1,nlev
          rayleigh_od_fl(:,jlev,jcol) = global_multiplier &
               &  * (pressure_hl(jcol,jlev+1) - pressure_hl(jcol,jlev)) * this%rayleigh_molar_scat
        end do
      end if

    end do

    if (lhook) call dr_hook('radiation_ecckd:calc_optical_depth',1,hook_handle)

  end subroutine calc_optical_depth_ckd_model


  !---------------------------------------------------------------------
  ! vectorized variant of above routine
  subroutine calc_optical_depth_ckd_model_vec(this, ncol, nlev, istartcol, iendcol, nmaxgas, &
       &  pressure_hl, temperature_fl, mole_fraction_fl, &
       &  optical_depth_fl, rayleigh_od_fl, opt_concentration_scaling)

    use ecradhook,             only : lhook, dr_hook, jphook
    use radiation_constants, only : accelduetogravity

    ! input variables

    class(ckd_model_type), intent(in), target  :: this
    ! number of columns, levels and input gases
    integer,               intent(in)  :: ncol, nlev, nmaxgas, istartcol, iendcol
    ! pressure at half levels (pa), dimensioned (ncol,nlev+1)
    real(jprb),            intent(in)  :: pressure_hl(ncol,nlev+1)
    ! temperature at full levels (k), dimensioned (ncol,nlev)
    real(jprb),            intent(in)  :: temperature_fl(istartcol:iendcol,nlev)
    ! gas mole fractions at full levels (mol mol-1), dimensioned (ncol,nlev,nmaxgas)
    real(jprb),            intent(in)  :: mole_fraction_fl(ncol,nlev,nmaxgas)
    ! optional concentration scaling of each gas
    real(jprb), optional,  intent(in)  :: opt_concentration_scaling(nmaxgas)
    
    ! output variables

    ! layer absorption optical depth for each g point
    real(jprb),            intent(out) :: optical_depth_fl(this%ng,nlev,istartcol:iendcol)
    ! in the shortwave only, the rayleigh scattering optical depth
    real(jprb),  optional, intent(out) :: rayleigh_od_fl(this%ng,nlev,istartcol:iendcol)

    ! local variables

    real(jprb), pointer :: molar_abs(:,:,:), molar_abs_conc(:,:,:,:)

    ! natural logarithm of pressure at full levels
    real(jprb) :: log_pressure_fl

    ! optical depth of single gas at one point in space versus
    ! spectral interval
    !real(jprb) :: od_single_gas(this%ng)

    real(jprb) :: multiplier, simple_multiplier(ncol,nlev), global_multiplier, temperature1
    real(jprb) :: scaling
    real(jprb) :: concentration_scaling(nmaxgas)

    ! indices and weights in temperature, pressure and concentration interpolation
    real(jprb) :: pindex1, tindex1, cindex1
    real(jprb) :: pw1(ncol,nlev), pw2(ncol,nlev), tw1(ncol,nlev), tw2(ncol,nlev), cw1(ncol,nlev), cw2(ncol,nlev)
    integer    :: ip1(ncol,nlev), it1(ncol,nlev), ic1(ncol,nlev)

    ! natural logarithm of mole fraction at one point
    real(jprb) :: log_conc

    ! minimum mole fraction in look-up-table
    real(jprb) :: mole_frac1

    ! layer absorption optical depth for each g point (memory layout adjusted to vectorization)
    real(jprb) :: od_fl(ncol,this%ng,nlev)

    integer :: jcol, jlev, jgas, igascode, jg

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_ecckd:calc_optical_depth_vec',0,hook_handle)

    global_multiplier = 1.0_jprb / (accelduetogravity * 0.001_jprb * airmolarmass)
    concentration_scaling = 1._jprb
    if (present(opt_concentration_scaling)) concentration_scaling = opt_concentration_scaling

    od_fl(:,:,:) = 0.0_jprb

    do jlev = 1,nlev
      do jcol = istartcol,iendcol

        log_pressure_fl = log(0.5_jprb * (pressure_hl(jcol,jlev)+pressure_hl(jcol,jlev+1)))

        ! find interpolation points in pressure
        pindex1 = (log_pressure_fl-this%log_pressure1) &
             &    / this%d_log_pressure
        pindex1 = 1.0_jprb + max(0.0_jprb, min(pindex1, this%npress-1.0001_jprb))
        ip1(jcol,jlev) = int(pindex1)
        pw2(jcol,jlev) = pindex1 - ip1(jcol,jlev)
        pw1(jcol,jlev) = 1.0_jprb - pw2(jcol,jlev)

        ! find interpolation points in temperature
        temperature1 = pw1(jcol,jlev)*this%temperature1(ip1(jcol,jlev)) &
             &       + pw2(jcol,jlev)*this%temperature1(ip1(jcol,jlev)+1)
        tindex1 = (temperature_fl(jcol,jlev) - temperature1) &
             &    / this%d_temperature
        tindex1 = 1.0_jprb + max(0.0_jprb, min(tindex1, this%ntemp-1.0001_jprb))
        it1(jcol,jlev) = int(tindex1)
        tw2(jcol,jlev) = tindex1 - it1(jcol,jlev)
        tw1(jcol,jlev) = 1.0_jprb - tw2(jcol,jlev)

        ! concentration multiplier
        simple_multiplier(jcol,jlev) = global_multiplier &
             &  * (pressure_hl(jcol,jlev+1) - pressure_hl(jcol,jlev))
      end do
    end do

    do jgas = 1,this%ngas

      associate (single_gas => this%single_gas(jgas))
        igascode = this%single_gas(jgas)%i_gas_code
          
        select case (single_gas%i_conc_dependence)
            
        case (iconcdependencelinear)
          molar_abs => this%single_gas(jgas)%molar_abs

          do jlev = 1,nlev
            do jg = 1, this%ng
              do jcol = istartcol,iendcol
                multiplier = simple_multiplier(jcol,jlev) * mole_fraction_fl(jcol,jlev,igascode) &
                   &                * concentration_scaling(igascode)

                od_fl(jcol,jg,jlev) = od_fl(jcol,jg,jlev) &
                   &        + (multiplier*tw1(jcol,jlev)) * (pw1(jcol,jlev) * molar_abs(jg,ip1(jcol,jlev),it1(jcol,jlev)) &
                   &                +pw2(jcol,jlev) * molar_abs(jg,ip1(jcol,jlev)+1,it1(jcol,jlev))) &
                   &        + (multiplier*tw2(jcol,jlev)) * (pw1(jcol,jlev) * molar_abs(jg,ip1(jcol,jlev),it1(jcol,jlev)+1) &
                   &                +pw2(jcol,jlev) * molar_abs(jg,ip1(jcol,jlev)+1,it1(jcol,jlev)+1))
              end do
            end do
          end do

        case (iconcdependencerelativelinear)
          molar_abs => this%single_gas(jgas)%molar_abs

          do jlev = 1,nlev
            do jg = 1, this%ng
              do jcol = istartcol,iendcol
                 multiplier = simple_multiplier(jcol,jlev)  * (mole_fraction_fl(jcol,jlev,igascode) &
                   & * concentration_scaling(igascode) &
                   & - single_gas%reference_mole_frac)

                od_fl(jcol,jg,jlev) = od_fl(jcol,jg,jlev) &
                   &        + (multiplier*tw1(jcol,jlev)) * (pw1(jcol,jlev) * molar_abs(jg,ip1(jcol,jlev),it1(jcol,jlev)) &
                   &                +pw2(jcol,jlev) * molar_abs(jg,ip1(jcol,jlev)+1,it1(jcol,jlev))) &
                   &        + (multiplier*tw2(jcol,jlev)) * (pw1(jcol,jlev) * molar_abs(jg,ip1(jcol,jlev),it1(jcol,jlev)+1) &
                   &                +pw2(jcol,jlev) * molar_abs(jg,ip1(jcol,jlev)+1,it1(jcol,jlev)+1))
              end do
            end do
          end do

        case (iconcdependencenone)
          ! composite gases
          molar_abs => this%single_gas(jgas)%molar_abs

          do jlev = 1,nlev
            do jg = 1, this%ng
              do jcol = istartcol,iendcol
                od_fl(jcol,jg,jlev) = od_fl(jcol,jg,jlev) &
                   &        + (simple_multiplier(jcol,jlev)*tw1(jcol,jlev)) *   &
                   &                (pw1(jcol,jlev) * molar_abs(jg,ip1(jcol,jlev),it1(jcol,jlev)) &
                   &                +pw2(jcol,jlev) * molar_abs(jg,ip1(jcol,jlev)+1,it1(jcol,jlev))) &
                   &        + (simple_multiplier(jcol,jlev)*tw2(jcol,jlev)) *   &
                   &                (pw1(jcol,jlev) * molar_abs(jg,ip1(jcol,jlev),it1(jcol,jlev)+1) &
                   &                +pw2(jcol,jlev) * molar_abs(jg,ip1(jcol,jlev)+1,it1(jcol,jlev)+1))
              end do
            end do
          end do

          case (iconcdependencelut)
            scaling = concentration_scaling(igascode)
            ! logarithmic interpolation in concentration space
            molar_abs_conc => this%single_gas(jgas)%molar_abs_conc
            mole_frac1 = exp(single_gas%log_mole_frac1)

            do jlev = 1,nlev
              do jcol = istartcol,iendcol
                ! take care of mole_fraction == 0
                log_conc = log(max(mole_fraction_fl(jcol,jlev,igascode)*scaling, mole_frac1))
                cindex1  = (log_conc - single_gas%log_mole_frac1) / single_gas%d_log_mole_frac
                cindex1  = 1.0_jprb + max(0.0_jprb, min(cindex1, single_gas%n_mole_frac-1.0001_jprb))
                ic1(jcol,jlev) = int(cindex1)
                cw2(jcol,jlev) = cindex1 - ic1(jcol,jlev)
                cw1(jcol,jlev) = 1.0_jprb - cw2(jcol,jlev)
              end do
            end do

          do jlev = 1,nlev
            do jg = 1, this%ng
!nec$ select_vector
              do jcol = istartcol,iendcol

                od_fl(jcol,jg,jlev) = od_fl(jcol,jg,jlev) &
                   &  + (simple_multiplier(jcol,jlev) * mole_fraction_fl(jcol,jlev,igascode) *scaling) * ( &
                   &      (cw1(jcol,jlev) * tw1(jcol,jlev) * pw1(jcol,jlev)) * &
                   &      molar_abs_conc(jg,ip1(jcol,jlev),it1(jcol,jlev),ic1(jcol,jlev)) &
                   &     +(cw1(jcol,jlev) * tw1(jcol,jlev) * pw2(jcol,jlev)) * &
                   &      molar_abs_conc(jg,ip1(jcol,jlev)+1,it1(jcol,jlev),ic1(jcol,jlev)) &
                   &     +(cw1(jcol,jlev) * tw2(jcol,jlev) * pw1(jcol,jlev)) * &
                   &      molar_abs_conc(jg,ip1(jcol,jlev),it1(jcol,jlev)+1,ic1(jcol,jlev)) &
                   &     +(cw1(jcol,jlev) * tw2(jcol,jlev) * pw2(jcol,jlev)) * &
                   &      molar_abs_conc(jg,ip1(jcol,jlev)+1,it1(jcol,jlev)+1,ic1(jcol,jlev)) &
                   &     +(cw2(jcol,jlev) * tw1(jcol,jlev) * pw1(jcol,jlev)) * &
                   &      molar_abs_conc(jg,ip1(jcol,jlev),it1(jcol,jlev),ic1(jcol,jlev)+1) &
                   &     +(cw2(jcol,jlev) * tw1(jcol,jlev) * pw2(jcol,jlev)) * &
                   &      molar_abs_conc(jg,ip1(jcol,jlev)+1,it1(jcol,jlev),ic1(jcol,jlev)+1) &
                   &     +(cw2(jcol,jlev) * tw2(jcol,jlev) * pw1(jcol,jlev)) * &
                   &      molar_abs_conc(jg,ip1(jcol,jlev),it1(jcol,jlev)+1,ic1(jcol,jlev)+1) &
                   &     +(cw2(jcol,jlev) * tw2(jcol,jlev) * pw2(jcol,jlev)) * &
                   &      molar_abs_conc(jg,ip1(jcol,jlev)+1,it1(jcol,jlev)+1,ic1(jcol,jlev)+1))
              end do
            end do
          end do
        end select
            
      end associate

      ! ensure the optical depth is not negative
      do jcol = istartcol,iendcol
        do jlev = 1,nlev
          do jg = 1, this%ng
            optical_depth_fl(jg,jlev,jcol) = max(0.0_jprb, od_fl(jcol,jg,jlev))
          end do
        end do
      end do

      ! rayleigh scattering
      if (this%is_sw .and. present(rayleigh_od_fl)) then
        do jcol = istartcol,iendcol
          do jlev = 1,nlev
            do jg = 1, this%ng
              rayleigh_od_fl(jg,jlev,jcol) = global_multiplier &
               &  * (pressure_hl(jcol,jlev+1) - pressure_hl(jcol,jlev)) * this%rayleigh_molar_scat(jg)
            end do
          end do
        end do
      end if

    end do

    if (lhook) call dr_hook('radiation_ecckd:calc_optical_depth_vec',1,hook_handle)

  end subroutine calc_optical_depth_ckd_model_vec

  
  !---------------------------------------------------------------------
  ! calculate the planck function integrated across each of the g
  ! points of this correlated k-distribution model, for a given
  ! temperature, where planck function is defined as the flux emitted
  ! by a black body (rather than radiance)
  subroutine calc_planck_function(this, nt, temperature, planck)

    class(ckd_model_type), intent(in)  :: this
    integer,    intent(in)  :: nt
    real(jprb), intent(in)  :: temperature(:) ! k
    real(jprb), intent(out) :: planck(this%ng,nt) ! w m-2

    real(jprb) :: tindex1, tw1, tw2
    integer    :: it1, jt

    do jt = 1,nt
      tindex1 = (temperature(jt) - this%temperature1_planck) &
           &   * (1.0_jprb / this%d_temperature_planck)
      if (tindex1 >= 0) then
        ! normal interpolation, and extrapolation for high temperatures
        tindex1 = 1.0_jprb + tindex1
        it1 = min(int(tindex1), this%nplanck-1)
        tw2 = tindex1 - it1
        tw1 = 1.0_jprb - tw2
        planck(:,jt) = tw1 * this%planck_function(:,it1) &
             &       + tw2 * this%planck_function(:,it1+1)
      else
        ! interpolate linearly to zero
        planck(:,jt) = this%planck_function(:,1) &
             &       * (temperature(jt)/this%temperature1_planck)
      end if
    end do

  end subroutine calc_planck_function
  

  !---------------------------------------------------------------------
  ! return the spectral solar irradiance integrated over each g point
  ! of a solar correlated k-distribution model, given the
  ! total_solar_irradiance
  subroutine calc_incoming_sw(this, total_solar_irradiance, &
       &                      spectral_solar_irradiance, &
       &                      solar_spectral_multiplier)

    use radiation_io, only : nulerr, radiation_abort

    class(ckd_model_type), intent(in)    :: this
    real(jprb),            intent(in)    :: total_solar_irradiance ! w m-2
    real(jprb),            intent(inout) :: spectral_solar_irradiance(:,:) ! w m-2
    real(jprb), optional,  intent(in)    :: solar_spectral_multiplier

    if (.not. present(solar_spectral_multiplier)) then
      spectral_solar_irradiance &
           &  = spread(total_solar_irradiance * this%norm_solar_irradiance, &
           &           2, size(spectral_solar_irradiance,2))
    else if (allocated(this%norm_amplitude_solar_irradiance)) then
      spectral_solar_irradiance &
           &  = spread(total_solar_irradiance * (this%norm_solar_irradiance &
           &   + solar_spectral_multiplier*this%norm_amplitude_solar_irradiance), &
           &           2, size(spectral_solar_irradiance,2))
    else if (solar_spectral_multiplier == 0.0_jprb) then
      spectral_solar_irradiance &
           &  = spread(total_solar_irradiance * this%norm_solar_irradiance, &
           &           2, size(spectral_solar_irradiance,2))
    else
      write(nulerr, '(a)') '*** error in calc_incoming_sw: no information present on solar cycle'
      call radiation_abort()
    end if

  end subroutine calc_incoming_sw

end module radiation_ecckd
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

