! # 1 "radiation/radiation_thermodynamics.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_thermodynamics.f90"
! this file has been modified for the use in icon

! radiation_thermodynamics.f90 - derived type for pressure & temperature
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
!   2017-05-11  r. hogan  fix startcol/endcol for get_layer_mass
!   2019-01-14  r. hogan  added out_of_physical_bounds routine
!   2019-01-14  r. hogan  capped h2o_sat_liq at 1

module radiation_thermodynamics

  use parkind1, only : jprb

  implicit none
  public

  !---------------------------------------------------------------------
  ! derived type for storing pressure and temperature at half and full levels
  type thermodynamics_type
     real(jprb), allocatable, dimension(:,:) :: &
          &  pressure_hl,    & ! (ncol,nlev+1) pressure (pa)
          &  temperature_hl, & ! (ncol,nlev+1) temperature (k)
          &  pressure_fl,    & ! (ncol,nlev) pressure (pa)
          &  temperature_fl    ! (ncol,nlev) temperature (k)

     ! the following is a function of pressure and temperature: you
     ! can calculate it according to your favourite formula, or the
     ! calc_saturation_wrt_liquid subroutine can be used to do this
     ! approximately
     real(jprb), allocatable, dimension(:,:) :: &
          &  h2o_sat_liq ! (ncol,nlev) specific humidity at liquid
                         ! saturation (kg/kg)

     ! using the interpolation method for temperature and pressure from half levels
     ! to full levels that is used in the subroutine gas_optics in radiation_ifs_rrtm
     ! can result in values for ind1 that exceed the bounds of absa within 
     ! srtm_taumol16 for ecrad inside icon. this can be avoided by directly
     ! passing pressure_fl and temperature_fl from icon to ecrad. with 
     ! rrtm_pass_temppres_fl = .true., the fields pressure_fl and temperature_fl
     ! are allocated and used within gas_optics in radiation_ifs_rrtm.
     logical :: &
          &  rrtm_pass_temppres_fl 
   contains
     procedure :: allocate   => allocate_thermodynamics_arrays
     procedure :: deallocate => deallocate_thermodynamics_arrays
     procedure :: calc_saturation_wrt_liquid
     procedure :: get_layer_mass
     procedure :: get_layer_mass_column
     procedure :: out_of_physical_bounds





  end type thermodynamics_type

contains


  !---------------------------------------------------------------------
  ! allocate variables with specified dimensions
  subroutine allocate_thermodynamics_arrays(this, ncol, nlev, &
       &                                    use_h2o_sat, rrtm_pass_temppres_fl)

    use ecradhook,  only : lhook, dr_hook, jphook

    class(thermodynamics_type), intent(inout) :: this
    integer, intent(in)           :: ncol  ! number of columns
    integer, intent(in)           :: nlev  ! number of levels
    logical, intent(in), optional :: use_h2o_sat ! allocate h2o_sat_liq?
    logical, intent(in), optional :: rrtm_pass_temppres_fl ! directly pass temperature
                                                           ! and pressure on full levels

    logical :: use_h2o_sat_local

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_thermodynamics:allocate',0,hook_handle)

    allocate(this%pressure_hl(ncol,nlev+1))
    allocate(this%temperature_hl(ncol,nlev+1))
    !$acc enter data create(this%pressure_hl, this%temperature_hl) async(1)

    use_h2o_sat_local = .false.
    if (present(use_h2o_sat)) then
      use_h2o_sat_local = use_h2o_sat
    end if

    this%rrtm_pass_temppres_fl = .false.
    if (present(rrtm_pass_temppres_fl)) then
      this%rrtm_pass_temppres_fl = rrtm_pass_temppres_fl
    end if

    if (this%rrtm_pass_temppres_fl) then
      allocate(this%pressure_fl(ncol,nlev))
      allocate(this%temperature_fl(ncol,nlev))
      !$acc enter data create(this%pressure_fl, this%temperature_fl) async(1)
    end if
    
    if (use_h2o_sat_local) then
      allocate(this%h2o_sat_liq(ncol,nlev))
      !$acc enter data create(this%h2o_sat_liq) async(1)
    end if    

    if (lhook) call dr_hook('radiation_thermodynamics:allocate',1,hook_handle)

  end subroutine allocate_thermodynamics_arrays


  !---------------------------------------------------------------------
  ! deallocate variables
  subroutine deallocate_thermodynamics_arrays(this)

    use ecradhook,  only : lhook, dr_hook, jphook

    class(thermodynamics_type), intent(inout) :: this

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_thermodynamics:deallocate',0,hook_handle)

    if (allocated(this%pressure_hl)) then
      !$acc exit data delete(this%pressure_hl) wait(1)
      deallocate(this%pressure_hl)
    end if
    if (allocated(this%temperature_hl)) then
      !$acc exit data delete(this%temperature_hl) wait(1)
      deallocate(this%temperature_hl)
    end if
    if (allocated(this%pressure_fl)) then
      !$acc exit data delete(this%pressure_fl) wait(1)
      deallocate(this%pressure_fl)
    end if
    if (allocated(this%temperature_fl)) then
      !$acc exit data delete(this%temperature_fl) wait(1)
      deallocate(this%temperature_fl)
    end if
    if (allocated(this%h2o_sat_liq)) then
      !$acc exit data delete(this%h2o_sat_liq) wait(1)
      deallocate(this%h2o_sat_liq)
    end if

    if (lhook) call dr_hook('radiation_thermodynamics:deallocate',1,hook_handle)
  
  end subroutine deallocate_thermodynamics_arrays


  !---------------------------------------------------------------------
  ! calculate approximate saturation with respect to liquid
  subroutine calc_saturation_wrt_liquid(this,istartcol,iendcol)

    use ecradhook,  only : lhook, dr_hook, jphook

    class(thermodynamics_type), intent(inout) :: this
    integer, intent(in)                       :: istartcol, iendcol

    ! pressure and temperature at full levels
    real(jprb) :: pressure, temperature

    ! vapour pressure (pa)
    real(jprb) :: e_sat

    integer :: ncol, nlev ! dimension sizes
    integer :: jcol, jlev ! loop indices for column and level

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_thermodynamics:calc_saturation_wrt_liquid',0,hook_handle)

    ncol = size(this%pressure_hl,1)
    nlev = size(this%pressure_hl,2) - 1

    if (.not. allocated(this%h2o_sat_liq)) then
      allocate(this%h2o_sat_liq(ncol,nlev))
      !$acc enter data create(this%h2o_sat_liq) async(1)
    end if

    !$acc parallel default(none) present(this) async(1)
    !$acc loop gang vector collapse(2) private(pressure, temperature, e_sat)
    do jlev = 1,nlev
       do jcol = istartcol,iendcol
          pressure = 0.5 * (this%pressure_hl(jcol,jlev)+this%pressure_hl(jcol,jlev+1))
          temperature = 0.5 * (this%temperature_hl(jcol,jlev)+this%temperature_hl(jcol,jlev+1))
          e_sat = 6.11e2_jprb * exp( 17.269_jprb * (temperature-273.16_jprb) / (temperature-35.86_jprb) )
          ! this formula can go above 1 at low pressure so needs to be
          ! capped
          this%h2o_sat_liq(jcol,jlev) = min(1.0_jprb, 0.622_jprb * e_sat / pressure)
       end do
    end do
    !$acc end parallel

    if (lhook) call dr_hook('radiation_thermodynamics:calc_saturation_wrt_liquid',1,hook_handle)

  end subroutine calc_saturation_wrt_liquid


  !---------------------------------------------------------------------
  ! calculate the dry mass of each layer, neglecting humidity effects.
  ! the first version is for all columns.
  subroutine get_layer_mass(this,istartcol,iendcol,layer_mass)

    use ecradhook,              only : lhook, dr_hook, jphook
    use radiation_constants,  only : accelduetogravity

    class(thermodynamics_type), intent(in)  :: this
    integer,                    intent(in)  :: istartcol, iendcol
    real(jprb),                 intent(out) :: layer_mass(:,:)

    integer    :: nlev
    real(jprb) :: inv_g

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_thermodynamics:get_layer_mass',0,hook_handle)

    nlev  = ubound(this%pressure_hl,2) - 1
    inv_g = 1.0_jprb / accelduetogravity

    layer_mass(istartcol:iendcol,1:nlev) &
         &  = ( this%pressure_hl(istartcol:iendcol,2:nlev+1) &
         &     -this%pressure_hl(istartcol:iendcol,1:nlev  )  ) &
         &  * inv_g 
    
    if (lhook) call dr_hook('radiation_thermodynamics:get_layer_mass',1,hook_handle)

  end subroutine get_layer_mass

  !---------------------------------------------------------------------
  ! calculate the dry mass of each layer, neglecting humidity effects.
  ! the second version is for one column, the one numbered "icol".
  subroutine get_layer_mass_column(this, icol, layer_mass)

    use ecradhook,              only : lhook, dr_hook, jphook
    use radiation_constants,  only : accelduetogravity

    class(thermodynamics_type), intent(in)  :: this
    integer,                    intent(in)  :: icol
    real(jprb),                 intent(out) :: layer_mass(:)

    integer    :: nlev
    real(jprb) :: inv_g

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_thermodynamics:get_layer_mass_column',0,hook_handle)

    nlev  = ubound(this%pressure_hl,2) - 1
    inv_g = 1.0_jprb / accelduetogravity

    layer_mass = ( this%pressure_hl(icol,2:nlev+1) &
             &    -this%pressure_hl(icol,1:nlev  )  ) &
             &   * inv_g
    
    if (lhook) call dr_hook('radiation_thermodynamics:get_layer_mass_column',1,hook_handle)

  end subroutine get_layer_mass_column


  !---------------------------------------------------------------------
  ! estimate the separation between the mid-points of model layers
  ! given the half-level pressure and temperature.  this is not in
  ! terms of the "thermodynamics" type as it is useful for computing
  ! overlap decorrelation lengths and hence cloud cover outside the
  ! radiation scheme.
  subroutine get_layer_separation(pressure_hl, temperature_hl, layer_separation)

    use ecradhook,              only : lhook, dr_hook, jphook
    use radiation_constants,  only : gasconstantdryair, accelduetogravity

    ! pressure (pa) and temperature (k) at half-levels, dimensioned
    ! (ncol,nlev+1) where ncol is the number of columns and nlev is
    ! the number of model levels
    real(jprb), dimension(:,:), intent(in)  :: pressure_hl, temperature_hl

    ! layer separation in metres, dimensioned (ncol,nlev-1)
    real(jprb), dimension(:,:), intent(out) :: layer_separation

    ! ratio of gas constant for dry air to acceleration due to gravity
    real(jprb), parameter :: r_over_g = gasconstantdryair / accelduetogravity

    ! loop indices and array bounds
    integer    :: jlev
    integer    :: i1, i2, nlev

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_thermodynamics:get_layer_separation',0,hook_handle)

    i1   = lbound(pressure_hl,1)
    i2   = ubound(pressure_hl,1)
    nlev =   size(pressure_hl,2)-1

    if (pressure_hl(i1,2) > pressure_hl(i1,1)) then
      ! pressure is increasing with index (order of layers is
      ! top-of-atmosphere to surface). in case pressure_hl(:,1)=0, we
      ! don't take the logarithm of the first pressure in each column.
      layer_separation(i1:i2,1) = r_over_g * temperature_hl(i1:i2,2) &
           &                    * log(pressure_hl(i1:i2,3)/pressure_hl(i1:i2,2))
      
      ! for other layers we take the separation between midpoints to
      ! be half the separation between the half-levels at the edge of
      ! the two adjacent layers
      do jlev = 2,nlev-1
        layer_separation(i1:i2,jlev) = (0.5_jprb * r_over_g) * temperature_hl(i1:i2,jlev+1) &
             &                    * log(pressure_hl(i1:i2,jlev+2)/pressure_hl(i1:i2,jlev))

      end do

    else
      ! pressure is decreasing with index (order of layers is surface
      ! to top-of-atmosphere).  in case pressure_hl(:,nlev+1)=0, we
      ! don't take the logarithm of the last pressure in each column.

      do jlev = 1,nlev-2
        layer_separation(i1:i2,jlev) = (0.5_jprb * r_over_g) * temperature_hl(i1:i2,jlev+1) &
             &                    * log(pressure_hl(i1:i2,jlev)/pressure_hl(i1:i2,jlev+2))

      end do
      layer_separation(i1:i2,nlev-1) = r_over_g * temperature_hl(i1:i2,nlev) &
           &                    * log(pressure_hl(i1:i2,nlev-1)/pressure_hl(i1:i2,nlev))

    end if

    if (lhook) call dr_hook('radiation_thermodynamics:get_layer_separation',1,hook_handle)    

  end subroutine get_layer_separation


  !---------------------------------------------------------------------
  ! return .true. if variables are out of a physically sensible range,
  ! optionally only considering columns between istartcol and iendcol
  function out_of_physical_bounds(this, istartcol, iendcol, do_fix) result(is_bad)

    use ecradhook,          only : lhook, dr_hook, jphook
    use radiation_check,  only : out_of_bounds_2d

    class(thermodynamics_type), intent(inout) :: this
    integer,           optional,intent(in) :: istartcol, iendcol
    logical,           optional,intent(in) :: do_fix
    logical                                :: is_bad

    logical    :: do_fix_local

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_thermodynamics:out_of_physical_bounds',0,hook_handle)

    if (present(do_fix)) then
      do_fix_local = do_fix
    else
      do_fix_local = .false.
    end if

    ! dangerous to cap pressure_hl as then the pressure difference across a layer could be zero
    is_bad =    out_of_bounds_2d(this%pressure_hl, 'pressure_hl', 0.0_jprb, 110000.0_jprb, &
         &                       .false., i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%temperature_hl, 'temperature_hl', 100.0_jprb,  400.0_jprb, &
         &                       do_fix_local, i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%h2o_sat_liq, 'h2o_sat_liq', 0.0_jprb, 1.0_jprb, &
         &                       do_fix_local, i1=istartcol, i2=iendcol)

    if (lhook) call dr_hook('radiation_thermodynamics:out_of_physical_bounds',1,hook_handle)

  end function out_of_physical_bounds

! # 425 "radiation/radiation_thermodynamics.f90"
  
end module radiation_thermodynamics
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

