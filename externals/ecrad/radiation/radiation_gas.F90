! # 1 "radiation/radiation_gas.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_gas.f90"
! this file has been modified for the use in icon

! radiation_gas.f90 - derived type to store the gas mixing ratios
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
!   2019-01-14  r. hogan  added out_of_physical_bounds routine

module radiation_gas

  use parkind1, only : jprb
  use radiation_gas_constants

  implicit none
  public

  ! available units
  enum, bind(c)
    enumerator imassmixingratio, ivolumemixingratio
  end enum

  !---------------------------------------------------------------------
  ! this derived type describes the gaseous composition of the
  ! atmosphere; gases may be stored as part of a 3d array (if their
  ! variation with height/column is to be represented) or one 1d array
  ! (if they are to be assumed globally well-mixed).
  type gas_type
    ! units of each stored gas (or 0 if not present)
    integer :: iunits(nmaxgases) = 0

    ! scaling factor that should be applied to each stored gas to get
    ! a dimensionless result, e.g. if iunits=ivolumemixingratio then
    ! 1.0e-6 is used to indicate the units are actually ppmv: need to
    ! multiply by 1e-6 to get mol/mol.
    real(jprb) :: scale_factor(nmaxgases) = 1.0_jprb

    ! mixing ratios of variable gases, dimensioned (ncol, nlev,
    ! nmaxgases)
    real(jprb), allocatable, dimension(:,:,:) :: mixing_ratio

    ! flag to indicate whether a gas is present
    logical :: is_present(nmaxgases) = .false.

    ! flag to indicate whether a gas is well mixed
    logical :: is_well_mixed(nmaxgases) = .false.

    integer :: ntype          = 0 ! number of gas types described

    integer :: ncol           = 0 ! number of columns in mixing_ratio
    integer :: nlev           = 0 ! number of levels  in mixing_ratio

    ! a list of length ntype of gases whose volume mixing ratios have
    ! been provided
    integer :: icode(nmaxgases) = 0

   contains
     procedure :: allocate   => allocate_gas
     procedure :: deallocate => deallocate_gas
     procedure :: put        => put_gas
     procedure :: put_well_mixed => put_well_mixed_gas
     procedure :: scale      => scale_gas
     procedure :: set_units  => set_units_gas
     procedure :: assert_units => assert_units_gas
     procedure :: get        => get_gas
     procedure :: get_scaling
     procedure :: reverse    => reverse_gas
     procedure :: out_of_physical_bounds





  end type gas_type

contains


  !---------------------------------------------------------------------
  ! allocate a derived type for holding gas mixing ratios given the
  ! number of columns and levels
  subroutine allocate_gas(this, ncol, nlev)

    use ecradhook, only : lhook, dr_hook, jphook

    class(gas_type), intent(inout) :: this
    integer,         intent(in)    :: ncol, nlev

    integer      :: jcol, jlev, jgas
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_gas:allocate',0,hook_handle)

    call this%deallocate()

    allocate(this%mixing_ratio(ncol, nlev, nmaxgases))
    !$acc enter data create(this%mixing_ratio) async(1)

    !$acc parallel default(none) present(this) async(1)
    !$acc loop gang vector collapse(3)
    do jgas = 1,nmaxgases
      do jlev = 1, nlev
        do jcol = 1,ncol
          this%mixing_ratio(jcol,jlev,jgas) = 0.0_jprb
        end do
      end do
    end do
    !$acc end parallel

    this%ncol = ncol
    this%nlev = nlev
    !$acc update device(this%ncol, this%nlev) async(1)

    if (lhook) call dr_hook('radiation_gas:allocate',1,hook_handle)

  end subroutine allocate_gas


  !---------------------------------------------------------------------
  ! deallocate memory and reset arrays
  subroutine deallocate_gas(this)

    use ecradhook, only : lhook, dr_hook, jphook

    class(gas_type), intent(inout) :: this

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_gas:deallocate',0,hook_handle)

    if (allocated(this%mixing_ratio)) then
       !$acc exit data delete(this%mixing_ratio) wait(1)
       deallocate(this%mixing_ratio)
    end if

    this%iunits = 0
    this%scale_factor = 0.0_jprb
    this%is_present = .false.
    this%is_well_mixed = .false.
    this%ntype = 0
    this%ncol = 0
    this%nlev = 0
    this%icode = 0

    if (lhook) call dr_hook('radiation_gas:deallocate',1,hook_handle)

  end subroutine deallocate_gas


  !---------------------------------------------------------------------
  ! put gas mixing ratio corresponding to gas id "igas" with units
  ! "iunits"
  subroutine put_gas(this, igas, iunits, mixing_ratio, scale_factor, &
       istartcol)

    use ecradhook,        only : lhook, dr_hook, jphook
    use radiation_io,   only : nulerr, radiation_abort

    class(gas_type),      intent(inout) :: this
    integer,              intent(in)    :: igas
    integer,              intent(in)    :: iunits
    real(jprb),           intent(in)    :: mixing_ratio(:,:)
    real(jprb), optional, intent(in)    :: scale_factor
    integer,    optional, intent(in)    :: istartcol

    integer :: i1, i2, jc, jk


    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_gas:put',0,hook_handle)

    ! check inputs
    if (igas <= igasnotpresent .or. iunits > nmaxgases) then
      write(nulerr,'(a,i0,a,i0,a,i0)') '*** error: provided gas id (', &
           &   igas, ') must be in the range ', igasnotpresent+1, ' to ', &
           &   nmaxgases
      call radiation_abort()
    end if
    if (iunits < imassmixingratio .or. iunits > ivolumemixingratio) then
      write(nulerr,'(a,i0,a,i0,a,i0)') '*** error: provided gas units (', &
           &   iunits, ') must be in the range ', imassmixingratio, ' to ', &
           &   ivolumemixingratio
      call radiation_abort()
    end if

    if (.not. allocated(this%mixing_ratio)) then
      write(nulerr,'(a,i0,a,i0,a,i0)') '*** error: attempt to put data to unallocated radiation_gas object'
      call radiation_abort()
    end if

    if (present(istartcol)) then
      i1 = istartcol
    else
      i1 = 1
    end if

    i2 = i1 + size(mixing_ratio,1) - 1

    if (i1 < 1 .or. i2 < 1 .or. i1 > this%ncol .or. i2 > this%ncol) then
      write(nulerr,'(a,i0,a,i0,a,i0)') '*** error: attempt to put columns indexed ', &
           &   i1, ' to ', i2, ' to array indexed 1 to ', this%ncol
      call radiation_abort()
    end if

    if (size(mixing_ratio,2) /= this%nlev) then
      write(nulerr,'(a,i0,a)') &
           &  '*** error: gas mixing ratio expected to have ', this%nlev, &
           &  ' levels'
      call radiation_abort()
    end if

    if (.not. this%is_present(igas)) then
      ! gas not present until now
      this%ntype = this%ntype + 1
      this%icode(this%ntype) = igas
      !$acc update device(this%icode(this%ntype:this%ntype)) async(1)
    end if
    this%is_present(igas) = .true.
    this%iunits(igas) = iunits
    this%is_well_mixed(igas) = .false.

    !$acc parallel default(none) present(this, mixing_ratio) async(1)
    !$acc loop gang vector collapse(2)
    do jk = 1,this%nlev
      do jc = i1,i2
        this%mixing_ratio(jc,jk,igas) = mixing_ratio(jc-i1+1,jk)
      end do
    end do
    !$acc end parallel
    if (present(scale_factor)) then
      this%scale_factor(igas) = scale_factor
    else
      this%scale_factor(igas) = 1.0_jprb
    end if

    if (lhook) call dr_hook('radiation_gas:put',1,hook_handle)

  end subroutine put_gas


  !---------------------------------------------------------------------
  ! put well-mixed gas mixing ratio corresponding to gas id "igas"
  ! with units "iunits"
  subroutine put_well_mixed_gas(this, igas, iunits, mixing_ratio, &
       scale_factor, istartcol, iendcol)

    use ecradhook,        only : lhook, dr_hook, jphook
    use radiation_io,   only : nulerr, radiation_abort

    class(gas_type),      intent(inout) :: this
    integer,              intent(in)    :: igas
    integer,              intent(in)    :: iunits
    real(jprb),           intent(in)    :: mixing_ratio
    real(jprb), optional, intent(in)    :: scale_factor
    integer,    optional, intent(in)    :: istartcol, iendcol

    real(jphook) :: hook_handle

    integer :: i1, i2, jc, jk

    if (lhook) call dr_hook('radiation_gas:put_well_mixed',0,hook_handle)

    ! check inputs
    if (igas <= igasnotpresent .or. igas > nmaxgases) then
      write(nulerr,'(a,i0,a,i0,a,i0)') '*** error: provided gas id (', &
           &   igas, ') must be in the range ', igasnotpresent+1, ' to ', &
           &   nmaxgases
      call radiation_abort()
    end if
    if (iunits < imassmixingratio .or. iunits > ivolumemixingratio) then
      write(nulerr,'(a,i0,a,i0,a,i0)') '*** error: provided gas units (', &
           &   iunits, ') must be in the range ', imassmixingratio, ' to ', &
           &   ivolumemixingratio
      call radiation_abort()
    end if

    if (.not. allocated(this%mixing_ratio)) then
      write(nulerr,'(a)') '*** error: attempt to put well-mixed gas data to unallocated radiation_gas object'
      call radiation_abort()
    end if

    if (present(istartcol)) then
      i1 = istartcol
    else
      i1 = 1
    end if

    if (present(iendcol)) then
      i2 = iendcol
    else
      i2 = this%ncol
    end if

    if (i1 < 1 .or. i2 < 1 .or. i1 > this%ncol .or. i2 > this%ncol) then
      write(nulerr,'(a,i0,a,i0,a,i0)') '*** error: attempt to put columns indexed ', &
           &   i1, ' to ', i2, ' to array indexed 1 to ', this%ncol
      call radiation_abort()
    end if

    if (.not. this%is_present(igas)) then
      ! gas not present until now
      this%ntype = this%ntype + 1
      this%icode(this%ntype) = igas
      !$acc update device(this%icode(this%ntype:this%ntype)) async(1)

    end if
    ! map uses a negative value to indicate a well-mixed value
    this%is_present(igas)              = .true.
    this%iunits(igas)                  = iunits
    this%is_well_mixed(igas)           = .true.

    !$acc parallel default(none) present(this) async(1)
    !$acc loop gang vector collapse(2) 
    do jk = 1,this%nlev
      do jc = i1,i2
        this%mixing_ratio(jc,jk,igas) = mixing_ratio
      end do
    end do
    !$acc end parallel
    if (present(scale_factor)) then
      this%scale_factor(igas) = scale_factor
    else
      this%scale_factor(igas) = 1.0_jprb
    end if

    if (lhook) call dr_hook('radiation_gas:put_well_mixed',1,hook_handle)

  end subroutine put_well_mixed_gas


  !---------------------------------------------------------------------
  ! scale gas concentrations, e.g. igas=ico2 and set scale_factor=2 to
  ! double co2.  note that this does not perform the scaling
  ! immediately, but changes the scale factor for the specified gas,
  ! ready to be used in set_units_gas.
  subroutine scale_gas(this, igas, scale_factor, lverbose)

    use radiation_io, only : nulout

    class(gas_type),      intent(inout) :: this
    integer,              intent(in)    :: igas
    real(jprb),           intent(in)    :: scale_factor
    logical,    optional, intent(in)    :: lverbose

    if (scale_factor /= 1.0_jprb) then
      this%scale_factor(igas) = this%scale_factor(igas) * scale_factor
      if (present(lverbose)) then
        if (lverbose) then
          write(nulout,'(a,a,a,f0.3)') '  scaling ', trim(gasname(igas)), &
               &  ' concentration by ', scale_factor
        end if
      end if
    end if

  end subroutine scale_gas


  !---------------------------------------------------------------------
  ! scale the gas concentrations so that they have the units "iunits"
  ! and are therefore ready to be used by the gas optics model within
  ! ecrad with no further scaling.  the existing scale_factor for each
  ! gas is applied.  if "igas" is present then apply only to gas with
  ! id "igas", otherwise to all gases. optional argument scale_factor
  ! specifies scaling that any subsequent access would need to apply
  ! to get a dimensionless result (consistent with definition of
  ! gas_type). so say that your gas optics model requires gas
  ! concentrations in ppmv, specify iunits=ivolumemixingratio and
  ! scale_factor=1.0e-6. if the gas concentrations were currently
  ! dimensionless volume mixing ratios, then the values would be
  ! internally divided by 1.0e-6.
  recursive subroutine set_units_gas(this, iunits, igas, scale_factor)
    class(gas_type),      intent(inout) :: this
    integer,              intent(in)    :: iunits
    integer,    optional, intent(in)    :: igas
    real(jprb), optional, intent(in)    :: scale_factor

    integer :: jg, jcol, jlev

    ! scaling factor to convert from old to new
    real(jprb) :: sf

    ! new scaling factor to store inside the gas object
    real(jprb) :: new_sf

    if (present(scale_factor)) then
      ! "sf" is the scaling to be applied now to the numbers (and may
      ! be modified below), "new_sf" is the value to be stored along
      ! with the numbers, informing subsequent routines how much you
      ! would need to multiply the numbers by to get a dimensionless
      ! result.
      sf     = 1.0_jprb / scale_factor
      new_sf = scale_factor
    else
      sf     = 1.0_jprb
      new_sf = 1.0_jprb
    end if

    if (present(igas)) then
      if (this%is_present(igas)) then
        if (iunits == imassmixingratio &
             &   .and. this%iunits(igas) == ivolumemixingratio) then
          sf = sf * gasmolarmass(igas) / airmolarmass
        else if (iunits == ivolumemixingratio &
             &   .and. this%iunits(igas) == imassmixingratio) then
          sf = sf * airmolarmass / gasmolarmass(igas)
        end if
        sf = sf * this%scale_factor(igas)

        if (sf /= 1.0_jprb) then
          !$acc parallel default(none) present(this, this%mixing_ratio, igas) async(1)
          !$acc loop gang vector collapse(2)
          do jlev = 1,this%nlev
            do jcol = 1,this%ncol
              this%mixing_ratio(jcol,jlev,igas) = this%mixing_ratio(jcol,jlev,igas) * sf
            enddo
          enddo
          !$acc end parallel
        end if
        ! store the new units and scale factor for this gas inside the
        ! gas object
        this%iunits(igas) = iunits
        this%scale_factor(igas) = new_sf
      end if
    else
      do jg = 1,this%ntype
        call this%set_units(iunits, igas=this%icode(jg), scale_factor=new_sf)
      end do
    end if

  end subroutine set_units_gas

  
  !---------------------------------------------------------------------
  ! return a vector indicating the scaling that one would need to
  ! apply to each gas in order to obtain the dimension units in
  ! "iunits" (which can be ivolumemixingratio or imassmixingratio)
  subroutine get_scaling(this, iunits, scaling)
    class(gas_type), intent(in)  :: this
    integer,         intent(in)  :: iunits
    real(jprb),      intent(out) :: scaling(nmaxgases)
    integer :: jg
    
    scaling = this%scale_factor
    do jg = 1,nmaxgases
      if (iunits == imassmixingratio .and. this%iunits(jg) == ivolumemixingratio) then
        scaling(jg) = scaling(jg) * gasmolarmass(jg) / airmolarmass
      else if (iunits == ivolumemixingratio .and. this%iunits(jg) == imassmixingratio) then
        scaling(jg) = scaling(jg) * airmolarmass / gasmolarmass(jg)
      end if
    end do
    
  end subroutine get_scaling

  
  !---------------------------------------------------------------------
  ! assert that gas mixing ratio units are "iunits", applying to gas
  ! with id "igas" if present, otherwise to all gases. otherwise the
  ! program will exit, except if the optional argument "istatus" is
  ! provided in which case it will return true if the units are
  ! correct and false if they are not. optional argument scale factor
  ! specifies any subsequent multiplication to apply; for ppmv one
  ! would use iunits=ivolumemixingratio and scale_factor=1.0e6.
  recursive subroutine assert_units_gas(this, iunits, igas, scale_factor, istatus)

    use radiation_io,   only : nulerr, radiation_abort

    class(gas_type),      intent(in)  :: this
    integer,              intent(in)  :: iunits
    integer,    optional, intent(in)  :: igas
    real(jprb), optional, intent(in)  :: scale_factor
    logical,    optional, intent(out) :: istatus

    integer :: jg

    real(jprb) :: sf

    if (present(scale_factor)) then
      sf = scale_factor
    else
      sf = 1.0_jprb
    end if

    if (present(istatus)) then
      istatus = .true.
    end if
    
    if (present(igas)) then
      if (this%is_present(igas)) then
        if (iunits /= this%iunits(igas)) then
          if (present(istatus)) then
            istatus = .false.
          else
            write(nulerr,'(a,a,a)') '*** error: ', trim(gasname(igas)), &
                 &  ' is not in the required units'
            call radiation_abort()
          end if
        else if (sf /= this%scale_factor(igas)) then
          if (present(istatus)) then
            istatus = .false.
          else
            write(nulerr,'(a,a,a,e12.4,a,e12.4)') '*** error: ', gasname(igas), &
                 &  ' scaling of ', this%scale_factor(igas), &
                 &  ' does not match required ', sf
            call radiation_abort()
          end if
        end if
      end if
    else
      do jg = 1,this%ntype
        call this%assert_units(iunits, igas=this%icode(jg), scale_factor=sf, istatus=istatus)
      end do
    end if

  end subroutine assert_units_gas


  !---------------------------------------------------------------------
  ! get gas mixing ratio corresponding to gas id "igas" with units
  ! "iunits" and return as a 2d array of dimensions (ncol,nlev).  the
  ! array will contain zeros if the gas is not stored.
  subroutine get_gas(this, igas, iunits, mixing_ratio, scale_factor, &
       &   istartcol)

    use ecradhook,        only : lhook, dr_hook, jphook
    use radiation_io,   only : nulerr, radiation_abort

    class(gas_type),      intent(in)  :: this
    integer,              intent(in)  :: igas
    integer,              intent(in)  :: iunits
    real(jprb),           intent(out) :: mixing_ratio(:,:)
    real(jprb), optional, intent(in)  :: scale_factor
    integer,    optional, intent(in)  :: istartcol

    real(jprb)                        :: sf
    integer                           :: i1, i2

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_gas:get',0,hook_handle)






    if (present(scale_factor)) then
      sf = scale_factor
    else
      sf = 1.0_jprb
    end if

    if (present(istartcol)) then
      i1 = istartcol
    else
      i1 = 1
    end if

    i2 = i1 + size(mixing_ratio,1) - 1

    if (i1 < 1 .or. i2 < 1 .or. i1 > this%ncol .or. i2 > this%ncol) then
      write(nulerr,'(a,i0,a,i0,a,i0)') '*** error: attempt to get columns indexed ', &
           &   i1, ' to ', i2, ' from array indexed 1 to ', this%ncol
      call radiation_abort()
    end if

    if (size(mixing_ratio,2) /= this%nlev) then
      write(nulerr,'(a,i0,a)') &
           &  '*** error: gas destination array expected to have ', this%nlev, &
           &  ' levels'
      call radiation_abort()
    end if

    if (.not. this%is_present(igas)) then
      mixing_ratio = 0.0_jprb
    else
      if (iunits == imassmixingratio &
           &   .and. this%iunits(igas) == ivolumemixingratio) then
        sf = sf * gasmolarmass(igas) / airmolarmass
      else if (iunits == ivolumemixingratio &
           &   .and. this%iunits(igas) == imassmixingratio) then
        sf = sf * airmolarmass / gasmolarmass(igas)
      end if
      sf = sf * this%scale_factor(igas)

      if (sf /= 1.0_jprb) then
        mixing_ratio = this%mixing_ratio(i1:i2,:,igas) * sf
      else
        mixing_ratio = this%mixing_ratio(i1:i2,:,igas)
      end if
    end if

    if (lhook) call dr_hook('radiation_gas:get',1,hook_handle)

  end subroutine get_gas


  !---------------------------------------------------------------------
  ! copy data to "gas_rev", reversing the height ordering of the gas
  ! data
  subroutine reverse_gas(this, istartcol, iendcol, gas_rev)

    class(gas_type), intent(in) :: this
    integer,        intent(in)  :: istartcol, iendcol
    type(gas_type), intent(out) :: gas_rev

    gas_rev%iunits = this%iunits
    gas_rev%scale_factor = this%scale_factor
    gas_rev%is_present = this%is_present
    gas_rev%is_well_mixed = this%is_well_mixed
    gas_rev%ntype = this%ntype
    gas_rev%ncol = this%ncol
    gas_rev%nlev = this%nlev
    gas_rev%icode = this%icode

    if (allocated(gas_rev%mixing_ratio)) deallocate(gas_rev%mixing_ratio)

    if (allocated(this%mixing_ratio)) then
      allocate(gas_rev%mixing_ratio(istartcol:iendcol,this%nlev,nmaxgases))
      gas_rev%mixing_ratio(istartcol:iendcol,:,:) &
           &  = this%mixing_ratio(istartcol:iendcol,this%nlev:1:-1,:)
    end if

  end subroutine reverse_gas

  !---------------------------------------------------------------------
  ! return .true. if variables are out of a physically sensible range,
  ! optionally only considering columns between istartcol and iendcol
  function out_of_physical_bounds(this, istartcol, iendcol, do_fix) result(is_bad)

    use ecradhook,          only : lhook, dr_hook, jphook
    use radiation_check,  only : out_of_bounds_3d

    class(gas_type),   intent(inout) :: this
    integer,  optional,intent(in) :: istartcol, iendcol
    logical,  optional,intent(in) :: do_fix
    logical                       :: is_bad

    logical    :: do_fix_local

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_gas:out_of_physical_bounds',0,hook_handle)

    if (present(do_fix)) then
      do_fix_local = do_fix
    else
      do_fix_local = .false.
    end if

    is_bad = out_of_bounds_3d(this%mixing_ratio, 'gas%mixing_ratio', &
         &                    0.0_jprb, 1.0_jprb, do_fix_local, i1=istartcol, i2=iendcol)

    if (lhook) call dr_hook('radiation_gas:out_of_physical_bounds',1,hook_handle)

  end function out_of_physical_bounds

! # 691 "radiation/radiation_gas.f90"

end module radiation_gas
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

