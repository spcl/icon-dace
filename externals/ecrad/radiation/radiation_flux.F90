! # 1 "radiation/radiation_flux.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_flux.f90"
! this file has been modified for the use in icon

! radiation_flux.f90 - derived type to store the output fluxes
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
!   2017-09-08  r. hogan  store g-point fluxes
!   2017-10-23  r. hogan  renamed single-character variables
!   2019-01-08  r. hogan  added "indexed_sum_profile"
!   2019-01-14  r. hogan  out_of_physical_bounds calls routine in radiation_config
!   2021-01-20  r. hogan  added heating_rate_out_of_physical_bounds function
!   2022-12-07  r. hogan  added top-of-atmosphere spectral output


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
! # 26 "radiation/radiation_flux.f90" 2

module radiation_flux

  use parkind1, only : jprb

  implicit none
  public

  !---------------------------------------------------------------------
  ! this derived type contains the output from the radiation
  ! calculation.  currently this is solely flux profiles, but in
  ! future surface fluxes in each band may be stored in order that the
  ! calling program can compute surface-radiation such as
  ! photosynthetically active radiation and uv index.
  type flux_type
     ! all the following are broad-band fluxes in w m-2 with
     ! dimensions (ncol,nlev+1).  note that only those fluxes that are
     ! requested will be used, so clear-sky and direct-beam arrays may
     ! not be allocated
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_up, lw_dn, &   ! upwelling and downwelling longwave
          &  sw_up, sw_dn, &   ! upwelling and downwelling shortwave
          &  sw_dn_direct, &   ! direct-beam shortwave into a horizontal plane
          &  lw_up_clear, lw_dn_clear, & ! clear-sky quantities...
          &  sw_up_clear, sw_dn_clear, &
          &  sw_dn_direct_clear
     ! as above but fluxes in each spectral band in w m-2 with
     ! dimensions (nband,ncol,nlev+1).  these are only allocated if
     ! config%do_save_spectral_flux==.true.
     real(jprb), allocatable, dimension(:,:,:) :: &
          &  lw_up_band, lw_dn_band, &   ! upwelling and downwelling longwave
          &  sw_up_band, sw_dn_band, &   ! upwelling and downwelling shortwave
          &  sw_dn_direct_band, &        ! direct-beam shortwave
          &  lw_up_clear_band, lw_dn_clear_band, & ! clear-sky quantities...
          &  sw_up_clear_band, sw_dn_clear_band, &
          &  sw_dn_direct_clear_band
     ! surface downwelling quantities at each g point, dimensioned
     ! (ng,ncol), that are always saved by the solver, except for the
     ! clear-sky ones that are only produced if
     ! config%do_clear==.true.
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_dn_surf_g, lw_dn_surf_clear_g, &
          &  sw_dn_diffuse_surf_g, sw_dn_direct_surf_g, &
          &  sw_dn_diffuse_surf_clear_g, sw_dn_direct_surf_clear_g
     ! top-of-atmosphere quantities at each g point, dimensioned
     ! (ng,ncol), that are always saved by the solver, except for the
     ! clear-sky ones that are only produced if
     ! config%do_clear==.true.
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_up_toa_g, lw_up_toa_clear_g, &
          &  sw_dn_toa_g, sw_up_toa_g, sw_up_toa_clear_g
     ! shortwave downwelling spectral fluxes in w m-2 at the surface,
     ! from which quantities such as photosynthetically active and uv
     ! radiation can be computed. only allocated if
     ! config%do_surface_sw_spectral_flux==.true.  note that the
     ! clear-sky quantities are only computed if
     ! config%do_clear==.true., but direct fluxes are computed whether
     ! or not do_direct==.true.. the dimensions are (nband,ncol).
     real(jprb), allocatable, dimension(:,:) :: &
          &  sw_dn_surf_band, sw_dn_direct_surf_band, &
          &  sw_dn_surf_clear_band, sw_dn_direct_surf_clear_band
     ! top-of-atmosphere spectral fluxes in w m-2. only allocated if
     ! config%do_toa_spectral_flux=.true.. note that the clear-sky
     ! quantities are only computed if config%do_clear==.true.. the
     ! dimensions are (nband,ncol).
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_up_toa_band, lw_up_toa_clear_band, &
          &  sw_dn_toa_band, sw_up_toa_band, sw_up_toa_clear_band
     ! surface downwelling fluxes in w m-2 at the spectral resolution
     ! needed by any subsequent canopy radiative transfer.  if
     ! config%use_canopy_full_spectrum_[sw|lw] then these will be at
     ! g-point resolution; otherwise they will be at
     ! config%n_albedo_bands and config%n_emiss_bands resolution.
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_dn_surf_canopy, &
          &  sw_dn_diffuse_surf_canopy, sw_dn_direct_surf_canopy

     ! diagnosed cloud cover from the short- and long-wave solvers
     real(jprb), allocatable, dimension(:) :: &
          &  cloud_cover_lw, cloud_cover_sw
     ! longwave derivatives needed by hogan and bozzo (2015) method
     ! for approximate longwave updates in between the full radiation
     ! calls: rate of change of upwelling broad-band flux with respect
     ! to surface value, dimensioned (ncol,nlev+1)
     real(jprb), allocatable, dimension(:,:) :: &
          &  lw_derivatives

   contains
     procedure :: allocate   => allocate_flux_type
     procedure :: deallocate => deallocate_flux_type
     procedure :: calc_surface_spectral
     procedure :: calc_toa_spectral
     procedure :: out_of_physical_bounds
     procedure :: heating_rate_out_of_physical_bounds




  end type flux_type

! added for dwd (2020)



      logical, parameter :: use_indexed_sum_vec = .false.


contains

  !---------------------------------------------------------------------
  ! allocate arrays for flux profiles, using config to define which
  ! fluxes are needed.  the arrays are dimensioned for columns between
  ! istartcol, iendcol and levels from 1 to nlev+1
  subroutine allocate_flux_type(this, config, istartcol, iendcol, nlev)

    use ecradhook,          only : lhook, dr_hook, jphook
    use radiation_io,     only : nulerr, radiation_abort
    use radiation_config, only : config_type

    integer, intent(in)             :: istartcol, iendcol, nlev
    class(flux_type), intent(inout) :: this
    type(config_type), intent(in)   :: config

    integer                         :: jcol
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_flux:allocate',0,hook_handle)

    ! allocate longwave arrays
    if (config%do_lw) then
      allocate(this%lw_up(istartcol:iendcol,nlev+1))
      allocate(this%lw_dn(istartcol:iendcol,nlev+1))
      !$acc enter data create(this%lw_up, this%lw_dn) async(1)
      if (config%do_clear) then
        allocate(this%lw_up_clear(istartcol:iendcol,nlev+1))
        allocate(this%lw_dn_clear(istartcol:iendcol,nlev+1))
        !$acc enter data create(this%lw_up_clear, this%lw_dn_clear) async(1)
      end if
      
      if (config%do_save_spectral_flux) then
        if (config%n_spec_lw == 0) then
          write(nulerr,'(a)') '*** error: number of lw spectral points to save not yet defined ' &
               & // 'so cannot allocate spectral flux arrays'
          call radiation_abort()
        end if
        
        allocate(this%lw_up_band(config%n_spec_lw,istartcol:iendcol,nlev+1))
        allocate(this%lw_dn_band(config%n_spec_lw,istartcol:iendcol,nlev+1))
        !$acc enter data create(this%lw_up_band, this%lw_dn_band) async(1)
        if (config%do_clear) then
          allocate(this%lw_up_clear_band(config%n_spec_lw, &
               &                         istartcol:iendcol,nlev+1))
          allocate(this%lw_dn_clear_band(config%n_spec_lw, &
               &                         istartcol:iendcol,nlev+1))
          !$acc enter data create(this%lw_up_clear_band, this%lw_dn_clear_band) &
          !$acc   async(1)
        end if
      end if
      
      if (config%do_lw_derivatives) then
        allocate(this%lw_derivatives(istartcol:iendcol,nlev+1))
        !$acc enter data create(this%lw_derivatives) async(1)
      end if

      if (config%do_toa_spectral_flux) then
        if (config%n_bands_lw == 0) then
          write(nulerr,'(a)') '*** error: number of lw bands not yet defined ' &
               & // 'so cannot allocate toa spectral flux arrays'
          call radiation_abort()
        end if
        allocate(this%lw_up_toa_band(config%n_bands_lw, istartcol:iendcol))
        !$acc enter data create(this%lw_up_toa_band) async(1)
        if (config%do_clear) then
          allocate(this%lw_up_toa_clear_band(config%n_bands_lw, istartcol:iendcol))
          !$acc enter data create(this%lw_up_toa_clear_band) async(1)
        end if
      end if
 
      ! allocate g-point downwelling fluxes at surface, and toa fluxes
      allocate(this%lw_dn_surf_g(config%n_g_lw,istartcol:iendcol))
      allocate(this%lw_up_toa_g (config%n_g_lw,istartcol:iendcol))
      !$acc enter data create(this%lw_dn_surf_g, this%lw_up_toa_g) async(1)
      if (config%do_clear) then
        allocate(this%lw_dn_surf_clear_g(config%n_g_lw,istartcol:iendcol))
        allocate(this%lw_up_toa_clear_g (config%n_g_lw,istartcol:iendcol))
        !$acc enter data create(this%lw_dn_surf_clear_g, this%lw_up_toa_clear_g) async(1)
      end if

      if (config%do_canopy_fluxes_lw) then
        ! downward fluxes at top of canopy at the spectral resolution
        ! used in the canopy radiative transfer scheme
        allocate(this%lw_dn_surf_canopy(config%n_canopy_bands_lw,istartcol:iendcol))
        !$acc enter data create(this%lw_dn_surf_canopy) async(1)
      end if
    end if
    
    ! allocate shortwave arrays
    if (config%do_sw) then
      allocate(this%sw_up(istartcol:iendcol,nlev+1))
      allocate(this%sw_dn(istartcol:iendcol,nlev+1))
      !$acc enter data create(this%sw_up, this%sw_dn) async(1)
      if (config%do_sw_direct) then
        allocate(this%sw_dn_direct(istartcol:iendcol,nlev+1))
        !$acc enter data create(this%sw_dn_direct) async(1)
      end if
      if (config%do_clear) then
        allocate(this%sw_up_clear(istartcol:iendcol,nlev+1))
        allocate(this%sw_dn_clear(istartcol:iendcol,nlev+1))
        !$acc enter data create(this%sw_up_clear, this%sw_dn_clear) async(1)
        if (config%do_sw_direct) then
          allocate(this%sw_dn_direct_clear(istartcol:iendcol,nlev+1))
          !$acc enter data create(this%sw_dn_direct_clear) async(1)
        end if
      end if
      
      if (config%do_save_spectral_flux) then
        if (config%n_spec_sw == 0) then
          write(nulerr,'(a)') '*** error: number of sw spectral points to save not yet defined ' &
               & // 'so cannot allocate spectral flux arrays'
          call radiation_abort()
        end if
        
        allocate(this%sw_up_band(config%n_spec_sw,istartcol:iendcol,nlev+1))
        allocate(this%sw_dn_band(config%n_spec_sw,istartcol:iendcol,nlev+1))
        !$acc enter data create(this%sw_up_band, this%sw_dn_band) async(1)
        
        if (config%do_sw_direct) then
          allocate(this%sw_dn_direct_band(config%n_spec_sw, &
               &                          istartcol:iendcol,nlev+1))
          !$acc enter data create(this%sw_dn_direct_band) async(1)
        end if
        if (config%do_clear) then
          allocate(this%sw_up_clear_band(config%n_spec_sw, &
               &                         istartcol:iendcol,nlev+1))
          allocate(this%sw_dn_clear_band(config%n_spec_sw, &
               &                         istartcol:iendcol,nlev+1))
          !$acc enter data create(this%sw_up_clear_band, this%sw_dn_clear_band) &
          !$acc   async(1)
          if (config%do_sw_direct) then
            allocate(this%sw_dn_direct_clear_band(config%n_spec_sw, &
                 &                                istartcol:iendcol, nlev+1))
            !$acc enter data create(this%sw_dn_direct_clear_band) async(1)
          end if
        end if
      end if
      
      if (config%do_surface_sw_spectral_flux) then
        if (config%n_bands_sw == 0) then
          write(nulerr,'(a)') '*** error: number of sw bands not yet defined ' &
               & // 'so cannot allocate toa spectral flux arrays'
          call radiation_abort()
        end if
        allocate(this%sw_dn_surf_band(config%n_bands_sw,istartcol:iendcol))
        allocate(this%sw_dn_direct_surf_band(config%n_bands_sw,istartcol:iendcol))
        !$acc enter data create(this%sw_dn_surf_band, this%sw_dn_direct_surf_band) &
        !$acc   async(1)
        if (config%do_clear) then
          allocate(this%sw_dn_surf_clear_band(config%n_bands_sw, &
               &                              istartcol:iendcol))
          allocate(this%sw_dn_direct_surf_clear_band(config%n_bands_sw, &
               &                                     istartcol:iendcol))
          !$acc enter data create(this%sw_dn_surf_clear_band, &
          !$acc   this%sw_dn_direct_surf_clear_band) async(1)
        end if
      end if

      if (config%do_toa_spectral_flux) then
        if (config%n_bands_sw == 0) then
          write(nulerr,'(a)') '*** error: number of sw bands not yet defined ' &
               & // 'so cannot allocate surface spectral flux arrays'
          call radiation_abort()
        end if
        allocate(this%sw_dn_toa_band(config%n_bands_sw, istartcol:iendcol))
        allocate(this%sw_up_toa_band(config%n_bands_sw, istartcol:iendcol))
        !$acc enter data create(this%sw_dn_toa_band, this%sw_up_toa_band) async(1)
        if (config%do_clear) then
          allocate(this%sw_up_toa_clear_band(config%n_bands_sw, istartcol:iendcol))
          !$acc enter data create(this%sw_up_toa_clear_band) async(1)
        end if
      end if
      
      ! allocate g-point downwelling fluxes at surface, and toa fluxes
      allocate(this%sw_dn_diffuse_surf_g(config%n_g_sw,istartcol:iendcol))
      allocate(this%sw_dn_direct_surf_g (config%n_g_sw,istartcol:iendcol))
      allocate(this%sw_dn_toa_g         (config%n_g_sw,istartcol:iendcol))
      allocate(this%sw_up_toa_g         (config%n_g_sw,istartcol:iendcol))
      !$acc enter data &
      !$acc   create(this%sw_dn_diffuse_surf_g,this%sw_dn_direct_surf_g) &
      !$acc   create(this%sw_dn_toa_g,this%sw_up_toa_g) &
      !$acc   async(1)
      if (config%do_clear) then
        allocate(this%sw_dn_diffuse_surf_clear_g(config%n_g_sw,istartcol:iendcol))
        allocate(this%sw_dn_direct_surf_clear_g (config%n_g_sw,istartcol:iendcol))
        allocate(this%sw_up_toa_clear_g         (config%n_g_sw,istartcol:iendcol))
        !$acc enter data create(this%sw_dn_diffuse_surf_clear_g, &
        !$acc   this%sw_dn_direct_surf_clear_g, this%sw_up_toa_clear_g) async(1)
      end if

      if (config%do_canopy_fluxes_sw) then
        ! downward fluxes at top of canopy at the spectral resolution
        ! used in the canopy radiative transfer scheme
        allocate(this%sw_dn_diffuse_surf_canopy(config%n_canopy_bands_sw,istartcol:iendcol))
        allocate(this%sw_dn_direct_surf_canopy (config%n_canopy_bands_sw,istartcol:iendcol))
        !$acc enter data create(this%sw_dn_diffuse_surf_canopy, &
        !$acc   this%sw_dn_direct_surf_canopy) async(1)
      end if
    end if
    
    ! allocate cloud cover arrays
    allocate(this%cloud_cover_lw(istartcol:iendcol))
    allocate(this%cloud_cover_sw(istartcol:iendcol))
    !$acc enter data create(this%cloud_cover_lw, this%cloud_cover_sw) async(1)

    !$acc wait ! accwa (nvhpc 22.7) crashes otherwise

    ! some solvers may not write to cloud cover, so we initialize to
    ! an unphysical value
    !$acc parallel default(none) present(this) async(1)
    !$acc loop gang vector
    do jcol = istartcol,iendcol
      this%cloud_cover_lw(jcol) = -1.0_jprb
      this%cloud_cover_sw(jcol) = -1.0_jprb
    end do
    !$acc end parallel

    if (lhook) call dr_hook('radiation_flux:allocate',1,hook_handle)
    
  end subroutine allocate_flux_type


  !---------------------------------------------------------------------
  ! deallocate flux arrays
  subroutine deallocate_flux_type(this)

    use ecradhook,          only : lhook, dr_hook, jphook

    class(flux_type), intent(inout) :: this
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_flux:deallocate',0,hook_handle)

    if (allocated(this%lw_up)) then
      !$acc exit data delete(this%lw_up) async(1)
      !$acc exit data delete(this%lw_dn) async(1) if(allocated(this%lw_dn))
      !$acc exit data delete(this%lw_up_clear) async(1) if(allocated(this%lw_up_clear))
      !$acc exit data delete(this%lw_dn_clear) async(1) if(allocated(this%lw_dn_clear))
      !$acc wait
      deallocate(this%lw_up)
      if (allocated(this%lw_dn))       deallocate(this%lw_dn)
      if (allocated(this%lw_up_clear)) deallocate(this%lw_up_clear)
      if (allocated(this%lw_dn_clear)) deallocate(this%lw_dn_clear)
    end if

    if (allocated(this%sw_up)) then
      !$acc exit data delete(this%sw_up) async(1)
      !$acc exit data delete(this%sw_dn) async(1) if(allocated(this%sw_dn))
      !$acc exit data delete(this%sw_up_clear) async(1) if(allocated(this%sw_up_clear))
      !$acc exit data delete(this%sw_dn_clear) async(1) if(allocated(this%sw_dn_clear))
      !$acc exit data delete(this%sw_dn_direct) async(1) if(allocated(this%sw_dn_direct))
      !$acc exit data delete(this%sw_dn_direct_clear) async(1) if(allocated(this%sw_dn_direct_clear))
      !$acc wait
      deallocate(this%sw_up)
      if (allocated(this%sw_dn))        deallocate(this%sw_dn)
      if (allocated(this%sw_up_clear))  deallocate(this%sw_up_clear)
      if (allocated(this%sw_dn_clear))  deallocate(this%sw_dn_clear)
      if (allocated(this%sw_dn_direct)) deallocate(this%sw_dn_direct)
      if (allocated(this%sw_dn_direct_clear)) &
           &   deallocate(this%sw_dn_direct_clear)
    end if

    if (allocated(this%lw_up_band)) then
      !$acc exit data delete(this%lw_up_band) async(1)
      !$acc exit data delete(this%lw_dn_band) async(1) if(allocated(this%lw_dn_band))
      !$acc exit data delete(this%lw_up_clear_band) async(1) if(allocated(this%lw_up_clear_band))
      !$acc exit data delete(this%lw_dn_clear_band) async(1) if(allocated(this%lw_dn_clear_band))
      !$acc wait
      deallocate(this%lw_up_band)
      if (allocated(this%lw_dn_band))       deallocate(this%lw_dn_band)
      if (allocated(this%lw_up_clear_band)) deallocate(this%lw_up_clear_band)
      if (allocated(this%lw_dn_clear_band)) deallocate(this%lw_dn_clear_band)
    end if
    
    if (allocated(this%sw_up_band)) then
      !$acc exit data delete(this%sw_up_band) async(1)
      !$acc exit data delete(this%sw_dn_band) async(1) if(allocated(this%sw_dn_band))
      !$acc exit data delete(this%sw_up_clear_band) async(1) if(allocated(this%sw_up_clear_band))
      !$acc exit data delete(this%sw_dn_clear_band) async(1) if(allocated(this%sw_dn_clear_band))
      !$acc exit data delete(this%sw_dn_direct_band) async(1) if(allocated(this%sw_dn_direct_band))
      !$acc exit data delete(this%sw_dn_direct_clear_band) async(1) &
      !$acc   if(allocated(this%sw_dn_direct_clear_band))
      !$acc wait
      deallocate(this%sw_up_band)
      if (allocated(this%sw_dn_band))        deallocate(this%sw_dn_band)
      if (allocated(this%sw_up_clear_band))  deallocate(this%sw_up_clear_band)
      if (allocated(this%sw_dn_clear_band))  deallocate(this%sw_dn_clear_band)
      if (allocated(this%sw_dn_direct_band)) deallocate(this%sw_dn_direct_band)
      if (allocated(this%sw_dn_direct_clear_band)) &
           &   deallocate(this%sw_dn_direct_clear_band)      
    end if

    if (allocated(this%sw_dn_surf_band)) then
      !$acc exit data delete(this%sw_dn_surf_band) async(1)
      !$acc exit data delete(this%sw_dn_direct_surf_band) async(1) if(allocated(this%sw_dn_direct_surf_band))
      !$acc wait
      deallocate(this%sw_dn_surf_band)
      deallocate(this%sw_dn_direct_surf_band)
    end if
    if (allocated(this%sw_dn_surf_clear_band)) then
      !$acc exit data delete(this%sw_dn_surf_clear_band) async(1)
      !$acc exit data delete(this%sw_dn_direct_surf_clear_band) async(1) &
      !$acc   if(allocated(this%sw_dn_direct_surf_clear_band))
      !$acc wait
      deallocate(this%sw_dn_surf_clear_band)
      deallocate(this%sw_dn_direct_surf_clear_band)
    end if

    !$acc exit data delete(this%lw_dn_surf_canopy) async(1) if(allocated(this%lw_dn_surf_canopy))
    !$acc exit data delete(this%sw_dn_diffuse_surf_canopy) async(1) &
    !$acc   if(allocated(this%sw_dn_diffuse_surf_canopy))
    !$acc exit data delete(this%sw_dn_direct_surf_canopy) async(1) &
    !$acc   if(allocated(this%sw_dn_direct_surf_canopy))
    !$acc wait
    if (allocated(this%lw_dn_surf_canopy)) deallocate(this%lw_dn_surf_canopy)
    if (allocated(this%sw_dn_diffuse_surf_canopy)) deallocate(this%sw_dn_diffuse_surf_canopy)
    if (allocated(this%sw_dn_direct_surf_canopy)) deallocate(this%sw_dn_direct_surf_canopy)

    if (allocated(this%cloud_cover_sw)) then
      !$acc exit data delete(this%cloud_cover_sw) wait(1)
      deallocate(this%cloud_cover_sw)
    end if
    if (allocated(this%cloud_cover_lw)) then
      !$acc exit data delete(this%cloud_cover_lw) wait(1)
      deallocate(this%cloud_cover_lw)
    end if

    if (allocated(this%lw_derivatives)) then
      !$acc exit data delete(this%lw_derivatives) wait(1)
      deallocate(this%lw_derivatives)
    end if

    !$acc exit data delete(this%lw_dn_surf_g) async(1) if(allocated(this%lw_dn_surf_g))
    !$acc exit data delete(this%lw_dn_surf_clear_g) async(1) if(allocated(this%lw_dn_surf_clear_g))
    !$acc exit data delete(this%sw_dn_diffuse_surf_g) async(1) if(allocated(this%sw_dn_diffuse_surf_g))
    !$acc exit data delete(this%sw_dn_direct_surf_g) async(1) if(allocated(this%sw_dn_direct_surf_g))
    !$acc exit data delete(this%sw_dn_diffuse_surf_clear_g) async(1) &
    !$acc   if(allocated(this%sw_dn_diffuse_surf_clear_g))
    !$acc exit data delete(this%sw_dn_direct_surf_clear_g) async(1) &
    !$acc   if(allocated(this%sw_dn_direct_surf_clear_g))
    !$acc wait
    if (allocated(this%lw_dn_surf_g))               deallocate(this%lw_dn_surf_g)
    if (allocated(this%lw_dn_surf_clear_g))         deallocate(this%lw_dn_surf_clear_g)
    if (allocated(this%sw_dn_diffuse_surf_g))       deallocate(this%sw_dn_diffuse_surf_g)
    if (allocated(this%sw_dn_direct_surf_g))        deallocate(this%sw_dn_direct_surf_g)
    if (allocated(this%sw_dn_diffuse_surf_clear_g)) deallocate(this%sw_dn_diffuse_surf_clear_g)
    if (allocated(this%sw_dn_direct_surf_clear_g))  deallocate(this%sw_dn_direct_surf_clear_g)

    !$acc exit data delete(this%lw_up_toa_g) async(1) if(allocated(this%lw_up_toa_g))
    !$acc exit data delete(this%sw_up_toa_g) async(1) if(allocated(this%sw_up_toa_g))
    !$acc exit data delete(this%sw_dn_toa_g) async(1) if(allocated(this%sw_dn_toa_g))
    !$acc exit data delete(this%lw_up_toa_clear_g) async(1) if(allocated(this%lw_up_toa_clear_g))
    !$acc exit data delete(this%sw_up_toa_clear_g) async(1) if(allocated(this%sw_up_toa_clear_g))
    !$acc wait
    if (allocated(this%lw_up_toa_g))                deallocate(this%lw_up_toa_g)
    if (allocated(this%sw_up_toa_g))                deallocate(this%sw_up_toa_g)
    if (allocated(this%sw_dn_toa_g))                deallocate(this%sw_dn_toa_g)
    if (allocated(this%lw_up_toa_clear_g))          deallocate(this%lw_up_toa_clear_g)
    if (allocated(this%sw_up_toa_clear_g))          deallocate(this%sw_up_toa_clear_g)

    if (lhook) call dr_hook('radiation_flux:deallocate',1,hook_handle)

  end subroutine deallocate_flux_type
  

  !---------------------------------------------------------------------
  ! calculate surface downwelling fluxes in each band using the
  ! downwelling surface fluxes at each g point
  subroutine calc_surface_spectral(this, config, istartcol, iendcol)

    use ecradhook,          only : lhook, dr_hook, jphook



    use radiation_config, only : config_type

    class(flux_type),  intent(inout) :: this
    type(config_type), intent(in)    :: config
    integer,           intent(in)    :: istartcol, iendcol

    integer :: jcol, jband, jalbedoband, nalbedoband

    ! longwave surface downwelling in each band needed to compute
    ! canopy fluxes
    real(jprb) :: lw_dn_surf_band(config%n_bands_lw,istartcol:iendcol)

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_flux:calc_surface_spectral',0,hook_handle)








    !$acc data present(config, this)

    if (config%do_sw .and. config%do_surface_sw_spectral_flux) then

      if (use_indexed_sum_vec) then
        call indexed_sum_vec(this%sw_dn_direct_surf_g, &
             &               config%i_band_from_reordered_g_sw, &
             &               this%sw_dn_direct_surf_band, istartcol, iendcol)
        call indexed_sum_vec(this%sw_dn_diffuse_surf_g, &
             &               config%i_band_from_reordered_g_sw, &
             &               this%sw_dn_surf_band, istartcol, iendcol)
        do jcol = istartcol,iendcol
          this%sw_dn_surf_band(:,jcol) &
               &  = this%sw_dn_surf_band(:,jcol) &
               &  + this%sw_dn_direct_surf_band(:,jcol)
        end do
      else
        !$acc parallel default(none) num_gangs(iendcol-istartcol+1) num_workers(1) &
        !$acc   vector_length(32*((config%n_g_sw-1)/32+1)) async(1)
        !$acc loop gang
        do jcol = istartcol,iendcol
          call indexed_sum(this%sw_dn_direct_surf_g(:,jcol), &
               &           config%i_band_from_reordered_g_sw, &
               &           this%sw_dn_direct_surf_band(:,jcol))
          call indexed_sum(this%sw_dn_diffuse_surf_g(:,jcol), &
               &           config%i_band_from_reordered_g_sw, &
               &           this%sw_dn_surf_band(:,jcol))
          this%sw_dn_surf_band(:,jcol) &
               &  = this%sw_dn_surf_band(:,jcol) &
               &  + this%sw_dn_direct_surf_band(:,jcol)
        end do
        !$acc end parallel
      end if

      if (config%do_clear) then
        if (use_indexed_sum_vec) then
          call indexed_sum_vec(this%sw_dn_direct_surf_clear_g, &
               &               config%i_band_from_reordered_g_sw, &
               &               this%sw_dn_direct_surf_clear_band, istartcol, iendcol)
          call indexed_sum_vec(this%sw_dn_diffuse_surf_clear_g, &
               &               config%i_band_from_reordered_g_sw, &
               &               this%sw_dn_surf_clear_band, istartcol, iendcol)
          do jcol = istartcol,iendcol
            this%sw_dn_surf_clear_band(:,jcol) &
                 &  = this%sw_dn_surf_clear_band(:,jcol) &
                 &  + this%sw_dn_direct_surf_clear_band(:,jcol)
          end do
        else
          !$acc parallel default(none) num_gangs(iendcol-istartcol+1) num_workers(1) &
          !$acc   vector_length(32*(config%n_g_sw-1)/32+1) async(1)
          !$acc loop gang
          do jcol = istartcol,iendcol
            call indexed_sum(this%sw_dn_direct_surf_clear_g(:,jcol), &
                 &           config%i_band_from_reordered_g_sw, &
                 &           this%sw_dn_direct_surf_clear_band(:,jcol))
            call indexed_sum(this%sw_dn_diffuse_surf_clear_g(:,jcol), &
                 &           config%i_band_from_reordered_g_sw, &
                 &           this%sw_dn_surf_clear_band(:,jcol))
            this%sw_dn_surf_clear_band(:,jcol) &
                 &  = this%sw_dn_surf_clear_band(:,jcol) &
                 &  + this%sw_dn_direct_surf_clear_band(:,jcol)
          end do
          !$acc end parallel
        end if
      end if

    end if ! do_surface_sw_spectral_flux

    ! fluxes in bands required for canopy radiative transfer
    if (config%do_sw .and. config%do_canopy_fluxes_sw) then
      if (config%use_canopy_full_spectrum_sw) then
        this%sw_dn_diffuse_surf_canopy(:,istartcol:iendcol) = this%sw_dn_diffuse_surf_g(:,istartcol:iendcol)
        this%sw_dn_direct_surf_canopy (:,istartcol:iendcol) = this%sw_dn_direct_surf_g (:,istartcol:iendcol)
      else if (config%do_nearest_spectral_sw_albedo) then
        if (use_indexed_sum_vec) then
          call indexed_sum_vec(this%sw_dn_direct_surf_g, &
               &               config%i_albedo_from_band_sw(config%i_band_from_reordered_g_sw), &
               &               this%sw_dn_direct_surf_canopy, istartcol, iendcol)
          call indexed_sum_vec(this%sw_dn_diffuse_surf_g, &
               &               config%i_albedo_from_band_sw(config%i_band_from_reordered_g_sw), &
               &               this%sw_dn_diffuse_surf_canopy, istartcol, iendcol)
        else
          do jcol = istartcol,iendcol
            call indexed_sum(this%sw_dn_direct_surf_g(:,jcol), &
                 &           config%i_albedo_from_band_sw(config%i_band_from_reordered_g_sw), &
                 &           this%sw_dn_direct_surf_canopy(:,jcol))
            call indexed_sum(this%sw_dn_diffuse_surf_g(:,jcol), &
                 &           config%i_albedo_from_band_sw(config%i_band_from_reordered_g_sw), &
                 &           this%sw_dn_diffuse_surf_canopy(:,jcol))
          end do
        end if
      else
        ! more accurate calculations using weights, but requires
        ! this%sw_dn_[direct_]surf_band to be defined, i.e.
        ! config%do_surface_sw_spectral_flux == .true.
        nalbedoband = size(config%sw_albedo_weights,1)
        this%sw_dn_diffuse_surf_canopy(:,istartcol:iendcol) = 0.0_jprb
        this%sw_dn_direct_surf_canopy (:,istartcol:iendcol) = 0.0_jprb
        do jband = 1,config%n_bands_sw
          do jalbedoband = 1,nalbedoband
            if (config%sw_albedo_weights(jalbedoband,jband) /= 0.0_jprb) then
              ! initially, "diffuse" is actually "total"
              this%sw_dn_diffuse_surf_canopy(jalbedoband,istartcol:iendcol) &
                   &  = this%sw_dn_diffuse_surf_canopy(jalbedoband,istartcol:iendcol) &
                   &  + config%sw_albedo_weights(jalbedoband,jband) &
                   &    * this%sw_dn_surf_band(jband,istartcol:iendcol)
              this%sw_dn_direct_surf_canopy(jalbedoband,istartcol:iendcol) &
                   &  = this%sw_dn_direct_surf_canopy(jalbedoband,istartcol:iendcol) &
                   &  + config%sw_albedo_weights(jalbedoband,jband) &
                   &    * this%sw_dn_direct_surf_band(jband,istartcol:iendcol)
            end if
          end do
        end do
        ! subtract the direct from total to get diffuse
        this%sw_dn_diffuse_surf_canopy(:,istartcol:iendcol) &
             &  = this%sw_dn_diffuse_surf_canopy(:,istartcol:iendcol) &
             &  - this%sw_dn_direct_surf_canopy(:,istartcol:iendcol)
      end if

    end if ! do_canopy_fluxes_sw

    if (config%do_lw .and. config%do_canopy_fluxes_lw) then
      if (config%use_canopy_full_spectrum_lw) then
        this%lw_dn_surf_canopy(:,istartcol:iendcol) = this%lw_dn_surf_g(:,istartcol:iendcol)
      else if (config%do_nearest_spectral_lw_emiss) then
        if (use_indexed_sum_vec) then
          call indexed_sum_vec(this%lw_dn_surf_g, &
               &               config%i_emiss_from_band_lw(config%i_band_from_reordered_g_lw), &
               &               this%lw_dn_surf_canopy, istartcol, iendcol)
        else
          do jcol = istartcol,iendcol
            call indexed_sum(this%lw_dn_surf_g(:,jcol), &
                 &           config%i_emiss_from_band_lw(config%i_band_from_reordered_g_lw), &
                 &           this%lw_dn_surf_canopy(:,jcol))
          end do
        end if
      else
        ! compute fluxes in each longwave emissivity interval using
        ! weights; first sum over g points to get the values in bands
        if (use_indexed_sum_vec) then
          call indexed_sum_vec(this%lw_dn_surf_g, &
               &               config%i_band_from_reordered_g_lw, &
               &               lw_dn_surf_band, istartcol, iendcol)
        else
          do jcol = istartcol,iendcol
            call indexed_sum(this%lw_dn_surf_g(:,jcol), &
                 &           config%i_band_from_reordered_g_lw, &
                 &           lw_dn_surf_band(:,jcol))
          end do
        end if
        nalbedoband = size(config%lw_emiss_weights,1)
        this%lw_dn_surf_canopy(:,istartcol:iendcol) = 0.0_jprb
        do jband = 1,config%n_bands_lw
          do jalbedoband = 1,nalbedoband
            if (config%lw_emiss_weights(jalbedoband,jband) /= 0.0_jprb) then
              this%lw_dn_surf_canopy(jalbedoband,istartcol:iendcol) &
                   &  = this%lw_dn_surf_canopy(jalbedoband,istartcol:iendcol) &
                   &  + config%lw_emiss_weights(jalbedoband,jband) &
                   &    * lw_dn_surf_band(jband,istartcol:iendcol)
            end if
          end do
        end do
      end if
    end if

    !$acc end data

    if (lhook) call dr_hook('radiation_flux:calc_surface_spectral',1,hook_handle)

  end subroutine calc_surface_spectral


  !---------------------------------------------------------------------
  ! calculate top-of-atmosphere fluxes in each band using the fluxes
  ! at each g point
  subroutine calc_toa_spectral(this, config, istartcol, iendcol)

    use ecradhook,          only : lhook, dr_hook, jphook
    use radiation_config, only : config_type





    class(flux_type),  intent(inout) :: this
    type(config_type), intent(in)    :: config
    integer,           intent(in)    :: istartcol, iendcol

    integer :: jcol, jband

    real(jphook) :: hook_handle
    
    if (lhook) call dr_hook('radiation_flux:calc_toa_spectral',0,hook_handle)








    if (config%do_sw .and. config%do_toa_spectral_flux) then

      if (use_indexed_sum_vec) then
        call indexed_sum_vec(this%sw_dn_toa_g, &
             &               config%i_band_from_reordered_g_sw, &
             &               this%sw_dn_toa_band, istartcol, iendcol)
        call indexed_sum_vec(this%sw_up_toa_g, &
             &               config%i_band_from_reordered_g_sw, &
             &               this%sw_up_toa_band, istartcol, iendcol)
      else
        do jcol = istartcol,iendcol
          call indexed_sum(this%sw_dn_toa_g(:,jcol), &
               &           config%i_band_from_reordered_g_sw, &
               &           this%sw_dn_toa_band(:,jcol))
          call indexed_sum(this%sw_up_toa_g(:,jcol), &
               &           config%i_band_from_reordered_g_sw, &
               &           this%sw_up_toa_band(:,jcol))
        end do
      end if
      
      if (config%do_clear) then
        if (use_indexed_sum_vec) then
          call indexed_sum_vec(this%sw_up_toa_clear_g, &
               &               config%i_band_from_reordered_g_sw, &
               &               this%sw_up_toa_clear_band, istartcol, iendcol)
        else
          do jcol = istartcol,iendcol
            call indexed_sum(this%sw_up_toa_clear_g(:,jcol), &
                 &               config%i_band_from_reordered_g_sw, &
                 &               this%sw_up_toa_clear_band(:,jcol))
          end do
        end if
      end if
    end if

    if (config%do_lw .and. config%do_toa_spectral_flux) then

      if (use_indexed_sum_vec) then
        call indexed_sum_vec(this%lw_up_toa_g, &
             &               config%i_band_from_reordered_g_lw, &
             &               this%lw_up_toa_band, istartcol, iendcol)
      else
        do jcol = istartcol,iendcol
          call indexed_sum(this%lw_up_toa_g(:,jcol), &
               &           config%i_band_from_reordered_g_lw, &
               &           this%lw_up_toa_band(:,jcol))
        end do
      end if
      
      if (config%do_clear) then
        if (use_indexed_sum_vec) then
          call indexed_sum_vec(this%lw_up_toa_clear_g, &
               &               config%i_band_from_reordered_g_lw, &
               &               this%lw_up_toa_clear_band, istartcol, iendcol)
        else
          do jcol = istartcol,iendcol
            call indexed_sum(this%lw_up_toa_clear_g(:,jcol), &
                 &               config%i_band_from_reordered_g_lw, &
                 &               this%lw_up_toa_clear_band(:,jcol))
          end do
        end if
      end if
    end if
    
    if (lhook) call dr_hook('radiation_flux:calc_toa_spectral',1,hook_handle)

  end subroutine calc_toa_spectral
  
    
  !---------------------------------------------------------------------
  ! return .true. if the most important flux variables are out of a
  ! physically sensible range, optionally only considering columns
  ! between istartcol and iendcol
  function out_of_physical_bounds(this, istartcol, iendcol) result(is_bad)

    use ecradhook,          only : lhook, dr_hook, jphook
    use radiation_check,  only : out_of_bounds_2d

    class(flux_type), intent(inout) :: this
    integer, optional,intent(in) :: istartcol, iendcol
    logical                      :: is_bad

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_flux:out_of_physical_bounds',0,hook_handle)

    is_bad =    out_of_bounds_2d(this%lw_up, 'lw_up', 10.0_jprb, 900.0_jprb, .false., i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%lw_dn, 'lw_dn', 0.0_jprb,  800.0_jprb, .false., i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%sw_up, 'sw_up', 0.0_jprb, 1500.0_jprb, .false., i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%sw_dn, 'sw_dn', 0.0_jprb, 1500.0_jprb, .false., i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%sw_dn_direct, 'sw_dn_direct', 0.0_jprb, 1500.0_jprb, .false., i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%lw_derivatives, 'lw_derivatives', 0.0_jprb, 1.0_jprb, .false., i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%sw_dn_surf_band, 'sw_dn_surf_band', 0.0_jprb, 1500.0_jprb, &
         &                       .false., j1=istartcol, j2=iendcol) &
         & .or. out_of_bounds_2d(this%sw_dn_surf_clear_band, 'sw_dn_surf_clear_band', 0.0_jprb, 1500.0_jprb, &
         &                       .false., j1=istartcol, j2=iendcol)

    if (lhook) call dr_hook('radiation_flux:out_of_physical_bounds',1,hook_handle)

  end function out_of_physical_bounds
  
  !---------------------------------------------------------------------
  ! return .true. if the heating rates are out of a physically
  ! sensible range, optionally only considering columns between
  ! istartcol and iendcol. this function allocates and deallocates
  ! memory due to the requirements for inputs of out_of_bounds_2d.
  function heating_rate_out_of_physical_bounds(this, nlev, istartcol, iendcol, pressure_hl) result(is_bad)
    
    use radiation_check, only : out_of_bounds_2d
    use radiation_constants, only : accelduetogravity

    ! "cp" (j kg-1 k-1)
    real(jprb), parameter :: specificheatdryair = 1004.0

    class(flux_type), intent(inout) :: this
    integer, intent(in) :: istartcol, iendcol, nlev
    logical                      :: is_bad
    
    real(jprb), intent(in) :: pressure_hl(:,:)

    real(jprb), allocatable :: hr_k_day(:,:)

    real(jprb) :: scaling(istartcol:iendcol,nlev)
    
    allocate(hr_k_day(istartcol:iendcol,nlev))

    scaling = -(24.0_jprb * 3600.0_jprb * accelduetogravity / specificheatdryair) &
         &  / (pressure_hl(istartcol:iendcol,2:nlev+1) - pressure_hl(istartcol:iendcol,1:nlev))
    ! shortwave
    hr_k_day = scaling * (this%sw_dn(istartcol:iendcol,2:nlev+1) - this%sw_up(istartcol:iendcol,2:nlev+1) &
         &               -this%sw_dn(istartcol:iendcol,1:nlev)   + this%sw_up(istartcol:iendcol,1:nlev))
    is_bad = out_of_bounds_2d(hr_k_day, 'sw_heating_rate_k_day', 0.0_jprb, 200.0_jprb, &
         &                    .false., i1=istartcol, i2=iendcol)

    ! longwave
    hr_k_day = scaling * (this%lw_dn(istartcol:iendcol,2:nlev+1) - this%lw_up(istartcol:iendcol,2:nlev+1) &
         &               -this%lw_dn(istartcol:iendcol,1:nlev)   + this%lw_up(istartcol:iendcol,1:nlev))
    is_bad = is_bad .or. out_of_bounds_2d(hr_k_day, 'lw_heating_rate_k_day', -250.0_jprb, 150.0_jprb, &
         &                                .false., i1=istartcol, i2=iendcol)

    deallocate(hr_k_day)

  end function heating_rate_out_of_physical_bounds


  !---------------------------------------------------------------------
  ! sum elements of "source" into "dest" according to index "ind".
  ! "source" and "ind" should have the same size and bounds, and no
  ! element of "ind" should refer outside the bounds of "dest".  this
  ! version increments existing contents of "dest".
  pure subroutine add_indexed_sum(source, ind, dest)

    real(jprb), intent(in)    :: source(:)
    integer,    intent(in)    :: ind(:)
    real(jprb), intent(inout) :: dest(:)

    integer :: ig, jg, istart, iend

    istart = lbound(source,1)
    iend   = ubound(source,1)

    do jg = istart, iend
      ig = ind(jg)
      dest(ig) = dest(ig) + source(jg)
    end do

  end subroutine add_indexed_sum


  !---------------------------------------------------------------------
  ! as "add_indexed_sum" but this version overwrites existing contents
  ! of "dest"
  pure subroutine indexed_sum(source, ind, dest)

    real(jprb), intent(in)  :: source(:)
    integer,    intent(in)  :: ind(:)
    real(jprb), intent(out) :: dest(:)

    integer :: ig, jg, istart, iend

    !$acc routine vector

    dest = 0.0

    istart = lbound(source,1)
    iend   = ubound(source,1)

    !$acc loop vector
    do jg = istart, iend
      ig = ind(jg)
      !$acc atomic update
      dest(ig) = dest(ig) + source(jg)
      !$acc end atomic
    end do

  end subroutine indexed_sum

  !---------------------------------------------------------------------
  ! vectorized version of "add_indexed_sum"
  subroutine indexed_sum_vec(source, ind, dest, ist, iend)

    real(jprb), intent(in)  :: source(:,:)
    integer,    intent(in)  :: ind(:)
    real(jprb), intent(out) :: dest(:,:)
    integer,    intent(in)  :: ist, iend

    integer :: ig, jg, jc

    dest = 0.0

    do jg = lbound(source,1), ubound(source,1)
      ig = ind(jg)
      do jc = ist, iend
        dest(ig,jc) = dest(ig,jc) + source(jg,jc)
      end do
    end do

  end subroutine indexed_sum_vec

  !---------------------------------------------------------------------
  ! as "add_indexed_sum" but a whole vertical profiles
  pure subroutine add_indexed_sum_profile(source, ind, dest)

    real(jprb), intent(in)  :: source(:,:)
    integer,    intent(in)  :: ind(:)
    real(jprb), intent(out) :: dest(:,:)

    integer :: ig, jg, istart, iend, jlev, nlev

    istart = lbound(source,1)
    iend   = ubound(source,1)
    nlev   = size(source,2)

    do jlev = 1,nlev
      do jg = istart, iend
        ig = ind(jg)
        dest(ig,jlev) = dest(ig,jlev) + source(jg,jlev)
      end do
    end do

  end subroutine add_indexed_sum_profile


  !---------------------------------------------------------------------
  ! as "indexed_sum" but a whole vertical profiles
  pure subroutine indexed_sum_profile(source, ind, dest)

    real(jprb), intent(in)  :: source(:,:)
    integer,    intent(in)  :: ind(:)
    real(jprb), intent(out) :: dest(:,:)

    integer :: ig, jg, istart, iend, jlev, nlev

    dest = 0.0

    istart = lbound(source,1)
    iend   = ubound(source,1)
    nlev   = size(source,2)

    do jlev = 1,nlev
      do jg = istart, iend
        ig = ind(jg)
        dest(ig,jlev) = dest(ig,jlev) + source(jg,jlev)
      end do
    end do

  end subroutine indexed_sum_profile
  
! # 1091 "radiation/radiation_flux.f90"

end module radiation_flux
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

