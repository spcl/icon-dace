! # 1 "radiation/radiation_aerosol.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_aerosol.f90"
! this file has been modified for the use in icon

! radiation_aerosol.f90 - derived type describing aerosol
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
!   2018-04-15  r. hogan  add "direct" option
!   2019-01-14  r. hogan  added out_of_physical_bounds routine

module radiation_aerosol

  use parkind1, only : jprb
  use radiation_io, only : nulerr, radiation_abort

  implicit none
  public

  !---------------------------------------------------------------------
  ! type describing the aerosol content in the atmosphere
  type aerosol_type
     ! the mass mixing ratio of config%n_aerosol_types different
     ! aerosol types dimensioned
     ! (ncol,istartlev:iendlev,config%n_aerosol_types), where ncol is
     ! the number of columns, istartlev:iendlev is the range of model
     ! levels where aerosols are present
     real(jprb), allocatable, dimension(:,:,:) :: &
          &  mixing_ratio  ! mass mixing ratio (kg/kg)

     ! alternatively, if is_direct=true, the optical properties are
     ! provided directly and are dimensioned
     ! (nband,istartlev:iendlev,ncol)
     real(jprb), allocatable, dimension(:,:,:) :: &
          &  od_sw, ssa_sw, g_sw, & ! shortwave optical properties
          &  od_lw, ssa_lw, g_lw    ! longwave optical properties

     ! range of levels in which the aerosol properties are provided
     integer :: istartlev, iendlev

     ! are the optical properties going to be provided directly by the
     ! user?
     logical :: is_direct = .false.

   contains
     procedure :: allocate        => allocate_aerosol_arrays
     procedure :: allocate_direct => allocate_aerosol_arrays_direct
     procedure :: deallocate      => deallocate_aerosol_arrays
     procedure :: out_of_physical_bounds




  end type aerosol_type

contains

  !---------------------------------------------------------------------
  ! allocate array for describing aerosols, although in the offline
  ! code these are allocated when they are read from the netcdf file
  subroutine allocate_aerosol_arrays(this, ncol, istartlev, iendlev, ntype)

    use ecradhook,     only : lhook, dr_hook, jphook

    class(aerosol_type), intent(inout) :: this
    integer, intent(in)                :: ncol  ! number of columns
    integer, intent(in)                :: istartlev, iendlev ! level range
    integer, intent(in)                :: ntype ! number of aerosol types
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol:allocate',0,hook_handle)

    allocate(this%mixing_ratio(ncol,istartlev:iendlev,ntype))
    !$acc enter data create(this%mixing_ratio)
    this%is_direct = .false.
    this%istartlev = istartlev
    this%iendlev   = iendlev






    if (lhook) call dr_hook('radiation_aerosol:allocate',1,hook_handle)

  end subroutine allocate_aerosol_arrays


  !---------------------------------------------------------------------
  ! allocate arrays for describing aerosol optical properties
  subroutine allocate_aerosol_arrays_direct(this, config, &
       &                                    ncol, istartlev, iendlev)

    use ecradhook,          only : lhook, dr_hook, jphook
    use radiation_config, only : config_type

    class(aerosol_type), intent(inout) :: this
    type(config_type),   intent(in)    :: config
    integer, intent(in)                :: ncol  ! number of columns
    integer, intent(in)                :: istartlev, iendlev ! level range
    integer                            :: jband, jlev, jcol

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol:allocate_direct',0,hook_handle)

    this%is_direct = .true.
    this%istartlev = istartlev
    this%iendlev   = iendlev

    if (config%do_sw) then
      allocate(this%od_sw (config%n_bands_sw,istartlev:iendlev,ncol))
      allocate(this%ssa_sw(config%n_bands_sw,istartlev:iendlev,ncol))
      allocate(this%g_sw  (config%n_bands_sw,istartlev:iendlev,ncol))
      !$acc enter data create(this%od_sw, this%ssa_sw, this%g_sw) async(1)
    end if

    if (config%do_lw) then
      allocate(this%od_lw (config%n_bands_lw,istartlev:iendlev,ncol))
      allocate(this%ssa_lw(config%n_bands_lw,istartlev:iendlev,ncol))
      allocate(this%g_lw  (config%n_bands_lw,istartlev:iendlev,ncol))
      !$acc enter data create(this%od_lw, this%ssa_lw, this%g_lw) async(1)
      ! if longwave scattering by aerosol is not to be represented,
      ! then the user may wish to just provide absorption optical
      ! depth in od_lw, in which case we must set the following two
      ! variables to zero

      !$acc wait ! accwa (nvhpc 22.7) crashes otherwise

      !$acc parallel default(none) present(this, config) async(1)
      !$acc loop gang vector collapse(3)
      do jcol = 1,ncol
        do jlev = istartlev,iendlev
          do jband = 1,config%n_bands_lw
            this%ssa_lw(jband,jlev,jcol) = 0.0_jprb
            this%g_lw(jband,jlev,jcol) = 0.0_jprb
          end do
        end do
      end do
      !$acc end parallel
    end if

    if (lhook) call dr_hook('radiation_aerosol:allocate_direct',1,hook_handle)

  end subroutine allocate_aerosol_arrays_direct


  !---------------------------------------------------------------------
  ! deallocate array
  subroutine deallocate_aerosol_arrays(this)

    use ecradhook,     only : lhook, dr_hook, jphook

    class(aerosol_type), intent(inout) :: this

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol:deallocate',0,hook_handle)

    !$acc exit data delete(this%mixing_ratio) async(1) if(allocated(this%mixing_ratio))
    !$acc exit data delete(this%od_sw) async(1) if(allocated(this%od_sw))
    !$acc exit data delete(this%ssa_sw) async(1) if(allocated(this%ssa_sw))
    !$acc exit data delete(this%g_sw) async(1) if(allocated(this%g_sw))
    !$acc exit data delete(this%od_lw) async(1) if(allocated(this%od_lw))
    !$acc exit data delete(this%ssa_lw) async(1) if(allocated(this%ssa_lw))
    !$acc exit data delete(this%g_lw) async(1) if(allocated(this%g_lw))
    !$acc wait
    if (allocated(this%mixing_ratio)) deallocate(this%mixing_ratio)
    if (allocated(this%od_sw))        deallocate(this%od_sw)
    if (allocated(this%ssa_sw))       deallocate(this%ssa_sw)
    if (allocated(this%g_sw))         deallocate(this%g_sw)
    if (allocated(this%od_lw))        deallocate(this%od_lw)
    if (allocated(this%ssa_lw))       deallocate(this%ssa_lw)
    if (allocated(this%g_lw))         deallocate(this%g_lw)
 
    if (lhook) call dr_hook('radiation_aerosol:deallocate',1,hook_handle)

  end subroutine deallocate_aerosol_arrays


  !---------------------------------------------------------------------
  ! return .true. if variables are out of a physically sensible range,
  ! optionally only considering columns between istartcol and iendcol
  function out_of_physical_bounds(this, istartcol, iendcol, do_fix) result(is_bad)

    use ecradhook,          only : lhook, dr_hook, jphook
    use radiation_check,  only : out_of_bounds_3d

    class(aerosol_type),   intent(inout) :: this
    integer,      optional,intent(in) :: istartcol, iendcol
    logical,      optional,intent(in) :: do_fix
    logical                           :: is_bad

    logical    :: do_fix_local

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol:out_of_physical_bounds',0,hook_handle)

    if (present(do_fix)) then
      do_fix_local = do_fix
    else
      do_fix_local = .false.
    end if

    is_bad =    out_of_bounds_3d(this%mixing_ratio, 'aerosol%mixing_ratio', &
         &                       0.0_jprb, 1.0_jprb, do_fix_local, i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_3d(this%od_sw, 'aerosol%od_sw', &
         &                       0.0_jprb, 100.0_jprb, do_fix_local, k1=istartcol, k2=iendcol) &
         & .or. out_of_bounds_3d(this%od_lw, 'aerosol%od_lw', &
         &                       0.0_jprb, 100.0_jprb, do_fix_local, k1=istartcol, k2=iendcol) &
         & .or. out_of_bounds_3d(this%ssa_sw, 'aerosol%ssa_sw', &
         &                       0.0_jprb, 1.0_jprb, do_fix_local, k1=istartcol, k2=iendcol) &
         & .or. out_of_bounds_3d(this%ssa_lw, 'aerosol%ssa_lw', &
         &                       0.0_jprb, 1.0_jprb, do_fix_local, k1=istartcol, k2=iendcol) &
         & .or. out_of_bounds_3d(this%g_sw, 'aerosol%g_sw', &
         &                       0.0_jprb, 1.0_jprb, do_fix_local, k1=istartcol, k2=iendcol) &
         & .or. out_of_bounds_3d(this%g_lw, 'aerosol%g_lw', &
         &                       0.0_jprb, 1.0_jprb, do_fix_local, k1=istartcol, k2=iendcol)

    if (lhook) call dr_hook('radiation_aerosol:out_of_physical_bounds',1,hook_handle)

  end function out_of_physical_bounds

! # 267 "radiation/radiation_aerosol.f90"
  
end module radiation_aerosol
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

