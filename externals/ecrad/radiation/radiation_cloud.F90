! # 1 "radiation/radiation_cloud.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_cloud.f90"
! this file has been modified for the use in icon

! radiation_cloud.f90 - derived type to store cloud/precip properties
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
!   2019-01-14  r. hogan  added inv_inhom_effective_size variable
!   2019-01-14  r. hogan  added out_of_physical_bounds routine
!   2019-06-14  r. hogan  added capability to store any number of cloud/precip types

module radiation_cloud

  use parkind1, only : jprb

  implicit none
  public

  !---------------------------------------------------------------------
  ! the intention is that all variables describing clouds and
  ! radiatively-active precipitation are contained in this derived
  ! type, and if cloud variables are to be added in future, they can
  ! be added to this type without requiring extra variables to be
  ! passed between subroutines elsewhere in the program.
  type cloud_type
    ! for maximum flexibility, an arbitrary number "ntype" of
    ! hydrometeor types can be stored, dimensioned (ncol,nlev,ntype)
    integer                                   :: ntype = 0
    real(jprb), allocatable, dimension(:,:,:) :: &
         &  mixing_ratio, &  ! mass mixing ratio (kg/kg)
         &  effective_radius ! (m)

    ! for backwards compatibility, we also allow for the two
    ! traditional cloud types, liquid cloud droplets and ice cloud
    ! particles, dimensioned (ncol,nlev)
    real(jprb), pointer, dimension(:,:) :: &
         &  q_liq,  q_ice,  & ! mass mixing ratio (kg/kg)
         &  re_liq, re_ice    ! effective radius (m)

    ! for the moment, the different types of hydrometeor are assumed
    ! to be mixed with each other, so there is just one cloud fraction
    ! variable varying from 0 to 1
    real(jprb), allocatable, dimension(:,:) :: fraction

    ! the fractional standard deviation of cloud optical depth in the
    ! cloudy part of the gridbox.  in the tripleclouds representation
    ! of cloud inhomogeneity, this is implemented by splitting the
    ! cloudy part of the gridbox into two equal-area regions, one
    ! with the cloud optical depth scaled by 1+fractional_std and the
    ! other scaled by 1-fractional_std. this variable is dimensioned
    ! (ncol,nlev)
    real(jprb), allocatable, dimension(:,:) :: fractional_std

    ! the inverse of the effective horizontal size of the clouds in
    ! the gridbox, used to compute the cloud edge length per unit
    ! gridbox area for use in representing 3d effects. this variable
    ! is dimensioned (ncol,nlev).
    real(jprb), allocatable, dimension(:,:) :: inv_cloud_effective_size ! m-1

    ! similarly for the in-cloud heterogeneities, used to compute the
    ! edge length between the optically thin and thick cloudy regions
    ! of the gridbox.
    real(jprb), allocatable, dimension(:,:) :: inv_inhom_effective_size ! m-1

    ! the following variable describes the overlap of cloud boundaries
    ! in adjacent layers, with dimensions (ncol,nlev-1): 1 corresponds
    ! to maximum overlap and 0 to random overlap. depending on the
    ! ecrad configuration, it may be the "alpha" overlap parameter of
    ! hogan and illingworth (2000) or the "beta" overlap parameter of
    ! shonk et al. (2010).
    real(jprb), allocatable, dimension(:,:) :: overlap_param

  contains
    procedure :: allocate   => allocate_cloud_arrays
    procedure :: deallocate => deallocate_cloud_arrays
    procedure :: set_overlap_param_fix
    procedure :: set_overlap_param_var
    generic   :: set_overlap_param => set_overlap_param_fix, set_overlap_param_var
    procedure :: set_overlap_param_approx
    procedure :: create_fractional_std
    procedure :: create_inv_cloud_effective_size
    procedure :: create_inv_cloud_effective_size_eta
    procedure :: param_cloud_effective_separation_eta
    procedure :: crop_cloud_fraction
    procedure :: out_of_physical_bounds





  end type cloud_type

contains

  !---------------------------------------------------------------------
  ! allocate arrays for describing clouds and precipitation, although
  ! in the offline code these are allocated when they are read from
  ! the netcdf file
  subroutine allocate_cloud_arrays(this, ncol, nlev, ntype, use_inhom_effective_size)

    use ecradhook,     only : lhook, dr_hook, jphook





    class(cloud_type), intent(inout), target :: this
    integer, intent(in)              :: ncol   ! number of columns
    integer, intent(in)              :: nlev   ! number of levels
    ! number of cloud/precip particle types.  if not present then the
    ! older cloud behaviour is assumed: two types are present, (1)
    ! liquid and (2) ice, and they can be accessed via q_liq, q_ice,
    ! re_liq and re_ice.
    integer, intent(in), optional    :: ntype
    logical, intent(in), optional    :: use_inhom_effective_size

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloud:allocate',0,hook_handle)

    if (present(ntype)) then
      this%ntype = ntype
    else
      this%ntype = 2
    end if
    !$acc update device(this%ntype)
    allocate(this%mixing_ratio(ncol,nlev,this%ntype))
    allocate(this%effective_radius(ncol,nlev,this%ntype))
    !$acc enter data create(this%mixing_ratio) async(1)
    !$acc enter data create(this%effective_radius) async(1)
    nullify(this%q_liq)
    nullify(this%q_ice)
    nullify(this%re_liq)
    nullify(this%re_ice)
    if (.not. present(ntype)) then
      ! older interface in which only liquid and ice are supported
      this%q_liq  => this%mixing_ratio(:,:,1)
      this%q_ice  => this%mixing_ratio(:,:,2)
      this%re_liq => this%effective_radius(:,:,1)
      this%re_ice => this%effective_radius(:,:,2)






    end if

    allocate(this%fraction(ncol,nlev))
    allocate(this%overlap_param(ncol,nlev-1))
    allocate(this%fractional_std(ncol,nlev))
    allocate(this%inv_cloud_effective_size(ncol,nlev))
    !$acc enter data create(this%fraction) async(1)
    !$acc enter data create(this%overlap_param) async(1)
    !$acc enter data create(this%fractional_std) async(1)
    !$acc enter data create(this%inv_cloud_effective_size) async(1)

    if (present(use_inhom_effective_size)) then
      if (use_inhom_effective_size) then
        allocate(this%inv_inhom_effective_size(ncol,nlev))
        !$acc enter data create(this%inv_inhom_effective_size) async(1)
      end if
    end if

    if (lhook) call dr_hook('radiation_cloud:allocate',1,hook_handle)

  end subroutine allocate_cloud_arrays


  !---------------------------------------------------------------------
  ! deallocate arrays
  subroutine deallocate_cloud_arrays(this)

    use ecradhook,     only : lhook, dr_hook, jphook

    class(cloud_type), intent(inout) :: this

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloud:deallocate',0,hook_handle)

    nullify(this%q_liq)
    nullify(this%q_ice)
    nullify(this%re_liq)
    nullify(this%re_ice)

    !$acc exit data delete(this%mixing_ratio) async(1) if(allocated(this%mixing_ratio))
    !$acc exit data delete(this%effective_radius) async(1) if(allocated(this%effective_radius))
    !$acc exit data delete(this%fraction) async(1) if(allocated(this%fraction))
    !$acc exit data delete(this%overlap_param) async(1) if(allocated(this%overlap_param))
    !$acc exit data delete(this%fractional_std) async(1) if(allocated(this%fractional_std))
    !$acc exit data delete(this%inv_cloud_effective_size) async(1) if(allocated(this%inv_cloud_effective_size))
    !$acc exit data delete(this%inv_inhom_effective_size) async(1) if(allocated(this%inv_inhom_effective_size))
    !$acc wait
    if (allocated(this%mixing_ratio))     deallocate(this%mixing_ratio)
    if (allocated(this%effective_radius)) deallocate(this%effective_radius)
    if (allocated(this%fraction))         deallocate(this%fraction)
    if (allocated(this%overlap_param))    deallocate(this%overlap_param)
    if (allocated(this%fractional_std))   deallocate(this%fractional_std)
    if (allocated(this%inv_cloud_effective_size)) &
         &  deallocate(this%inv_cloud_effective_size)
    if (allocated(this%inv_inhom_effective_size)) &
         &  deallocate(this%inv_inhom_effective_size)

    if (lhook) call dr_hook('radiation_cloud:deallocate',1,hook_handle)

  end subroutine deallocate_cloud_arrays


  !---------------------------------------------------------------------
  ! compute and store the overlap parameter from the provided overlap
  ! decorrelation length (in metres).  if istartcol and/or iendcol are
  ! provided then only columns in this range are computed.  if the
  ! overlap_param array has not been allocated then it will be
  ! allocated to be of the correct size relative to the pressure
  ! field. this version assumes a fixed decorrelation_length for all
  ! columns.
  subroutine set_overlap_param_fix(this, thermodynamics, decorrelation_length, &
       &  istartcol, iendcol)

    use ecradhook,                  only : lhook, dr_hook, jphook
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_constants,      only : gasconstantdryair, accelduetogravity

    class(cloud_type),         intent(inout) :: this
    type(thermodynamics_type), intent(in)    :: thermodynamics
    real(jprb),                intent(in)    :: decorrelation_length ! m
    integer,         optional, intent(in)    :: istartcol, iendcol

    ! ratio of gas constant for dry air to acceleration due to gravity
    real(jprb), parameter :: r_over_g = gasconstantdryair / accelduetogravity

    ! process only columns i1 to i2, which will be istartcol to
    ! iendcol if they were provided
    integer :: i1, i2

    integer :: ncol, nlev

    integer :: jcol, jlev

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloud:set_overlap_param_fix',0,hook_handle)

    ! pressure at half-levels, pressure_hl, is defined at nlev+1
    ! points
    ncol = size(thermodynamics%pressure_hl,dim=1)
    nlev = size(thermodynamics%pressure_hl,dim=2)-1

    if (present(istartcol)) then
      i1 = istartcol
    else
      i1 = 1
    end if

    if (present(iendcol)) then
      i2 = iendcol
    else
      i2 = ncol
    end if

    if (.not. allocated(this%overlap_param)) then
      ! if pressure is of size (ncol,nlev+1) then overlap_param is of
      ! size (ncol,nlev-1), since overlap parameter is only defined here
      ! for interfaces between model layers, not for the interface to
      ! space or the surface
      allocate(this%overlap_param(ncol, nlev-1))
      !$acc enter data create(this%overlap_param) async(1)
    end if

    !$acc data present(this, thermodynamics)

    !$acc update host(thermodynamics%pressure_hl(i1,1:2)) wait(1)
    if (thermodynamics%pressure_hl(i1,2) > thermodynamics%pressure_hl(i1,1)) then
      ! pressure is increasing with index (order of layers is
      ! top-of-atmosphere to surface). in case pressure_hl(:,1)=0, we
      ! don't take the logarithm of the first pressure in each column.
      !$acc parallel default(none) async(1)
      !$acc loop gang(static:1) vector
      do jcol = i1,i2
        this%overlap_param(jcol,1) = exp(-(r_over_g/decorrelation_length) &
             &                            * thermodynamics%temperature_hl(jcol,2) &
             &                            *log(thermodynamics%pressure_hl(jcol,3) &
             &                                /thermodynamics%pressure_hl(jcol,2)))
      end do

      !$acc loop seq
      do jlev = 2,nlev-1
        !$acc loop gang(static:1) vector
        do jcol = i1,i2
          this%overlap_param(jcol,jlev) = exp(-(0.5_jprb*r_over_g/decorrelation_length) &
              &                            * thermodynamics%temperature_hl(jcol,jlev+1) &
              &                            *log(thermodynamics%pressure_hl(jcol,jlev+2) &
              &                                /thermodynamics%pressure_hl(jcol,jlev)))
        end do
      end do
      !$acc end parallel

    else
       ! pressure is decreasing with index (order of layers is surface
       ! to top-of-atmosphere).  in case pressure_hl(:,nlev+1)=0, we
       ! don't take the logarithm of the last pressure in each column.
      !$acc parallel default(none) async(1)
      !$acc loop seq
      do jlev = 1,nlev-2
        !$acc loop gang(static:1) vector
        do jcol = i1,i2
          this%overlap_param(jcol,jlev) = exp(-(0.5_jprb*r_over_g/decorrelation_length) &
              &                            * thermodynamics%temperature_hl(jcol,jlev+1) &
              &                            *log(thermodynamics%pressure_hl(jcol,jlev) &
              &                                /thermodynamics%pressure_hl(jcol,jlev+2)))
        end do
      end do

      !$acc loop gang(static:1) vector
      do jcol = i1,i2
        this%overlap_param(jcol,nlev-1) = exp(-(r_over_g/decorrelation_length) &
            &                            * thermodynamics%temperature_hl(jcol,nlev) &
            &                            *log(thermodynamics%pressure_hl(jcol,nlev-1) &
            &                                /thermodynamics%pressure_hl(jcol,nlev)))
      end do
      !$acc end parallel
    end if

    !$acc end data

    if (lhook) call dr_hook('radiation_cloud:set_overlap_param_fix',1,hook_handle)

  end subroutine set_overlap_param_fix


  !---------------------------------------------------------------------
  ! compute and store the overlap parameter from the provided overlap
  ! decorrelation length (in metres), which may vary with column. only
  ! columns from istartcol to iendcol are computed.  if the
  ! overlap_param array has not been allocated then it will be
  ! allocated to be of the correct size relative to the pressure
  ! field.
  subroutine set_overlap_param_var(this, thermodynamics, decorrelation_length, &
       &                           istartcol, iendcol)

    use ecradhook,                  only : lhook, dr_hook, jphook
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_constants,      only : gasconstantdryair, accelduetogravity




    class(cloud_type),         intent(inout) :: this
    type(thermodynamics_type), intent(in)    :: thermodynamics
    integer,                   intent(in)    :: istartcol, iendcol
    real(jprb),                intent(in)    :: decorrelation_length(istartcol:iendcol) ! m

    ! ratio of gas constant for dry air to acceleration due to gravity
    real(jprb), parameter :: r_over_g = gasconstantdryair / accelduetogravity

    integer :: ncol, nlev

    integer :: jcol, jlev

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloud:set_overlap_param_var',0,hook_handle)

    ! pressure at half-levels, pressure_hl, is defined at nlev+1
    ! points
    ncol = size(thermodynamics%pressure_hl,dim=1)
    nlev = size(thermodynamics%pressure_hl,dim=2)-1

    if (.not. allocated(this%overlap_param)) then
      ! if pressure is of size (ncol,nlev+1) then overlap_param is of
      ! size (ncol,nlev-1), since overlap parameter is only defined here
      ! for interfaces between model layers, not for the interface to
      ! space or the surface
      allocate(this%overlap_param(ncol, nlev-1))
      !$acc enter data create(this%overlap_param) async(1)
    end if

    !$acc data present(this, thermodynamics, decorrelation_length)
    !$acc update host(thermodynamics%pressure_hl(istartcol,1:2)) wait(1)
    if (thermodynamics%pressure_hl(istartcol,2) > thermodynamics%pressure_hl(istartcol,1)) then
      ! pressure is increasing with index (order of layers is
      ! top-of-atmosphere to surface). in case pressure_hl(:,1)=0, we
      ! don't take the logarithm of the first pressure in each column.
      !$acc parallel default(none) async(1)
      !$acc loop gang(static:1) vector
      do jcol = istartcol,iendcol
        this%overlap_param(jcol,1) = exp(-(r_over_g/decorrelation_length(jcol)) &
             &                            * thermodynamics%temperature_hl(jcol,2) &
             &                            *log(thermodynamics%pressure_hl(jcol,3) &
             &                                /thermodynamics%pressure_hl(jcol,2)))
      end do

      !$acc loop seq
      do jlev = 2,nlev-1
        !$acc loop gang(static:1) vector
        do jcol = istartcol,iendcol
          this%overlap_param(jcol,jlev) = exp(-(0.5_jprb*r_over_g/decorrelation_length(jcol)) &
              &                            * thermodynamics%temperature_hl(jcol,jlev+1) &
              &                            *log(thermodynamics%pressure_hl(jcol,jlev+2) &
              &                                /thermodynamics%pressure_hl(jcol,jlev)))
        end do
      end do
      !$acc end parallel

    else
       ! pressure is decreasing with index (order of layers is surface
       ! to top-of-atmosphere).  in case pressure_hl(:,nlev+1)=0, we
       ! don't take the logarithm of the last pressure in each column.
      !$acc parallel default(none) async(1)
      !$acc loop seq
      do jlev = 1,nlev-2
        !$acc loop gang(static:1) vector
        do jcol = istartcol,iendcol
          this%overlap_param(jcol,jlev) = exp(-(0.5_jprb*r_over_g/decorrelation_length(jcol)) &
              &                            * thermodynamics%temperature_hl(jcol,jlev+1) &
              &                            *log(thermodynamics%pressure_hl(jcol,jlev) &
              &                                /thermodynamics%pressure_hl(jcol,jlev+2)))
        end do
      end do

      !$acc loop gang(static:1) vector
      do jcol = istartcol,iendcol
        this%overlap_param(jcol,nlev-1) = exp(-(r_over_g/decorrelation_length(jcol)) &
            &                            * thermodynamics%temperature_hl(jcol,nlev) &
            &                            *log(thermodynamics%pressure_hl(jcol,nlev-1) &
            &                                /thermodynamics%pressure_hl(jcol,nlev)))
      end do
      !$acc end parallel
    end if

    !$acc end data

    if (lhook) call dr_hook('radiation_cloud:set_overlap_param_var',1,hook_handle)

  end subroutine set_overlap_param_var


  !---------------------------------------------------------------------
  ! compute and store the overlap parameter from the provided overlap
  ! decorrelation length (in metres).  if istartcol and/or iendcol are
  ! provided then only columns in this range are computed.  if the
  ! overlap_param array has not been allocated then it will be
  ! allocated to be of the correct size relative to the pressure
  ! field. this is the approximate method as it assumes a fixed
  ! atmospheric scale height, which leads to differences particularly
  ! in low cloud.
  subroutine set_overlap_param_approx(this, thermodynamics, decorrelation_length, &
       &  istartcol, iendcol)

    use ecradhook,                  only : lhook, dr_hook, jphook
    use radiation_thermodynamics, only : thermodynamics_type

    class(cloud_type),         intent(inout) :: this
    type(thermodynamics_type), intent(in)    :: thermodynamics
    real(jprb),                intent(in)    :: decorrelation_length ! m
    integer,         optional, intent(in)    :: istartcol, iendcol

    ! to convert decorrelation length (m) to overlap parameter between
    ! layers, we need an estimate for the thickness of the layer. this
    ! is found using the pressure difference between the edges of the
    ! layer, along with the approximate scale height of the atmosphere
    ! (m) given here:
    real(jprb), parameter :: scale_height = 8000.0_jprb

    ! process only columns i1 to i2, which will be istartcol to
    ! iendcol if they were provided
    integer :: i1, i2

    integer :: ncol, nlev

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloud:set_overlap_param_approx',0,hook_handle)

    ! pressure at half-levels, pressure_hl, is defined at nlev+1
    ! points
    ncol = size(thermodynamics%pressure_hl,dim=1)
    nlev = size(thermodynamics%pressure_hl,dim=2)-1

    if (present(istartcol)) then
      i1 = istartcol
    else
      i1 = 1
    end if

    if (present(iendcol)) then
      i2 = iendcol
    else
      i2 = ncol
    end if

    if (.not. allocated(this%overlap_param)) then
      ! if pressure is of size (ncol,nlev+1) then overlap_param is of
      ! size (ncol,nlev-1), since overlap parameter is only defined here
      ! for interfaces between model layers, not for the interface to
      ! space or the surface
      allocate(this%overlap_param(ncol, nlev-1))
    end if

    if (thermodynamics%pressure_hl(i1,2) > thermodynamics%pressure_hl(i1,1)) then
       ! pressure is increasing with index (order of layers is
       ! top-of-atmosphere to surface). in case pressure_hl(:,1)=0, we
       ! don't take the logarithm of the first pressure in each
       ! column.
       this%overlap_param(i1:i2,:) = exp(-(scale_height/decorrelation_length) &
            &  * ( log(thermodynamics%pressure_hl(i1:i2,3:nlev+1) &
            &         /thermodynamics%pressure_hl(i1:i2,2:nlev  )) ) )
    else
       ! pressure is decreasing with index (order of layers is surface
       ! to top-of-atmosphere).  in case pressure_hl(:,nlev+1)=0, we
       ! don't take the logarithm of the last pressure in each column.
       this%overlap_param(i1:i2,:) = exp(-(scale_height/decorrelation_length) &
            &  * ( log(thermodynamics%pressure_hl(i1:i2,1:nlev-1) &
            &         /thermodynamics%pressure_hl(i1:i2,2:nlev  )) ) )
    end if

    if (lhook) call dr_hook('radiation_cloud:set_overlap_param_approx',1,hook_handle)

  end subroutine set_overlap_param_approx


  !---------------------------------------------------------------------
  ! create a matrix of constant fractional standard deviations
  ! (dimensionless)
  subroutine create_fractional_std(this, ncol, nlev, frac_std)

    use ecradhook,                  only : lhook, dr_hook, jphook

    class(cloud_type), intent(inout) :: this
    integer,           intent(in)    :: ncol, nlev
    real(jprb),        intent(in)    :: frac_std

    integer :: jcol, jlev

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloud:create_fractional_std',0,hook_handle)

    if (allocated(this%fractional_std)) then
      !$acc exit data delete(this%fractional_std) wait(1)
      !$acc wait ! accwa (nvhpc 22.7, nvhpc 22.5) crashes otherwise
       deallocate(this%fractional_std)
    end if
    
    allocate(this%fractional_std(ncol, nlev))
    !$acc enter data create(this%fractional_std) async(1)

    !$acc parallel default(none) present(this) async(1)
    !$acc loop gang vector collapse(2)
    do jlev = 1, nlev
      do jcol = 1, ncol
      this%fractional_std(jcol, jlev) = frac_std
      end do
    end do
    !$acc end parallel

    if (lhook) call dr_hook('radiation_cloud:create_fractional_std',1,hook_handle)

  end subroutine create_fractional_std


  !---------------------------------------------------------------------
  ! create a matrix of constant inverse cloud effective size (m-1)
  subroutine create_inv_cloud_effective_size(this, ncol, nlev, inv_eff_size)

    use ecradhook,                  only : lhook, dr_hook, jphook

    class(cloud_type), intent(inout) :: this
    integer,           intent(in)    :: ncol, nlev
    real(jprb),        intent(in)    :: inv_eff_size

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloud:create_inv_cloud_effective_size',0,hook_handle)

    if (allocated(this%inv_cloud_effective_size)) then
       deallocate(this%inv_cloud_effective_size)
    end if
    
    allocate(this%inv_cloud_effective_size(ncol, nlev))

    this%inv_cloud_effective_size = inv_eff_size

    if (lhook) call dr_hook('radiation_cloud:create_inv_cloud_effective_size',1,hook_handle)

  end subroutine create_inv_cloud_effective_size


  !---------------------------------------------------------------------
  ! create a matrix of inverse cloud effective size (m-1) according to
  ! the value of eta (=pressure divided by surface pressure)
  subroutine create_inv_cloud_effective_size_eta(this, ncol, nlev, &
       &  pressure_hl, inv_eff_size_low, inv_eff_size_mid, inv_eff_size_high, &
       &  eta_low_mid, eta_mid_high, istartcol, iendcol)

    use ecradhook,                  only : lhook, dr_hook, jphook

    class(cloud_type), intent(inout) :: this
    integer,           intent(in)    :: ncol, nlev
    ! pressure on half levels (pa)
    real(jprb),        intent(in)    :: pressure_hl(:,:)
    ! inverse effective size for low, mid and high cloud (m-1)
    real(jprb),        intent(in)    :: inv_eff_size_low
    real(jprb),        intent(in)    :: inv_eff_size_mid
    real(jprb),        intent(in)    :: inv_eff_size_high
    ! eta values at low-mid and mid-high interfaces
    real(jprb),        intent(in)    :: eta_low_mid, eta_mid_high
    integer, optional, intent(in)    :: istartcol, iendcol

    ! ratio of layer midpoint pressure to surface pressure
    real(jprb) :: eta(nlev)

    ! indices of column, level and surface half-level
    integer :: jcol, isurf

    ! local values of istartcol, iendcol
    integer :: i1, i2

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloud:create_inv_cloud_effective_size_eta',0,hook_handle)

    if (allocated(this%inv_cloud_effective_size)) then
      deallocate(this%inv_cloud_effective_size)
    end if
    
    allocate(this%inv_cloud_effective_size(ncol, nlev))

    if (present(istartcol)) then
      i1 = istartcol
    else
      i1 = 1
    end if

    if (present(iendcol)) then
      i2 = iendcol
    else
      i2 = ncol
    end if

    ! locate the surface half-level
    if (pressure_hl(1,1) > pressure_hl(1,2)) then
      isurf = 1
    else
      isurf = nlev+1
    end if

    do jcol = i1,i2
      eta = (pressure_hl(jcol,1:nlev)+pressure_hl(jcol,2:nlev+1)) &
           &  * (0.5_jprb / pressure_hl(jcol,isurf))
      where (eta > eta_low_mid)
        this%inv_cloud_effective_size(jcol,:) = inv_eff_size_low
      elsewhere (eta > eta_mid_high)
        this%inv_cloud_effective_size(jcol,:) = inv_eff_size_mid
      elsewhere
        this%inv_cloud_effective_size(jcol,:) = inv_eff_size_high
      end where
    end do

    if (lhook) call dr_hook('radiation_cloud:create_inv_cloud_effective_size_eta',1,hook_handle)

  end subroutine create_inv_cloud_effective_size_eta


  !---------------------------------------------------------------------
  ! create a matrix of inverse cloud and inhomogeneity effective size
  ! (m-1) parameterized according to the value of eta (=pressure
  ! divided by surface pressure): effective_separation =
  ! coeff_a + coeff_b*exp(-(eta**power)).  
  subroutine param_cloud_effective_separation_eta(this, ncol, nlev, &
       &  pressure_hl, separation_surf, separation_toa, power, &
       &  inhom_separation_factor, istartcol, iendcol)

    use ecradhook,                  only : lhook, dr_hook, jphook

    class(cloud_type), intent(inout) :: this
    integer,           intent(in)    :: ncol, nlev
    ! pressure on half levels (pa)
    real(jprb),        intent(in)    :: pressure_hl(:,:)
    ! separation distances at surface and top-of-atmosphere, and power
    ! on eta
    real(jprb),           intent(in) :: separation_surf ! m
    real(jprb),           intent(in) :: separation_toa ! m
    real(jprb),           intent(in) :: power
    real(jprb), optional, intent(in) :: inhom_separation_factor
    integer,    optional, intent(in) :: istartcol, iendcol

    ! ratio of layer midpoint pressure to surface pressure
    real(jprb) :: eta(nlev)

    ! effective cloud separation (m)
    real(jprb) :: eff_separation(nlev)

    ! coefficients used to compute effective separation distance
    real(jprb) :: coeff_e, coeff_a, coeff_b, inhom_sep_factor

    ! indices of column, level and surface half-level
    integer :: jcol, isurf

    ! local values of istartcol, iendcol
    integer :: i1, i2

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloud:param_cloud_effective_separation_eta',0,hook_handle)

    if (present(inhom_separation_factor)) then
      inhom_sep_factor = inhom_separation_factor
    else
      inhom_sep_factor = 1.0_jprb
    end if

    coeff_e = 1.0_jprb - exp(-1.0_jprb)
    coeff_b = (separation_toa - separation_surf) / coeff_e
    coeff_a = separation_toa - coeff_b

    if (allocated(this%inv_cloud_effective_size)) then
      deallocate(this%inv_cloud_effective_size)
    end if
     if (allocated(this%inv_inhom_effective_size)) then
      deallocate(this%inv_inhom_effective_size)
    end if
   
    allocate(this%inv_cloud_effective_size(ncol, nlev))
    allocate(this%inv_inhom_effective_size(ncol, nlev))

    if (present(istartcol)) then
      i1 = istartcol
    else
      i1 = 1
    end if

    if (present(iendcol)) then
      i2 = iendcol
    else
      i2 = ncol
    end if

    ! locate the surface half-level
    if (pressure_hl(1,1) > pressure_hl(1,2)) then
      isurf = 1
    else
      isurf = nlev+1
    end if

    do jcol = i1,i2
      eta = (pressure_hl(jcol,1:nlev)+pressure_hl(jcol,2:nlev+1)) &
           &  * (0.5_jprb / pressure_hl(jcol,isurf))
      eff_separation = coeff_a + coeff_b * exp(-eta**power)
      this%inv_cloud_effective_size(jcol,:) = 1.0_jprb / (eff_separation &
           &  * sqrt(max(1.0e-5_jprb,this%fraction(jcol,:)*(1.0_jprb-this%fraction(jcol,:)))))
      this%inv_inhom_effective_size(jcol,:) = 1.0_jprb / (eff_separation * inhom_sep_factor &
           &  * sqrt(max(1.0e-5_jprb,0.5_jprb*this%fraction(jcol,:)*(1.0_jprb-0.5_jprb*this%fraction(jcol,:)))))
    end do

    if (lhook) call dr_hook('radiation_cloud:param_cloud_effective_separation_eta',1,hook_handle)

  end subroutine param_cloud_effective_separation_eta


  !---------------------------------------------------------------------
  ! remove "ghost" clouds: those with a cloud fraction that is too
  ! small to treat sensibly (e.g. because it implies that the
  ! "in-cloud" water content is too high), or with a cloud water
  ! content that is too small.  we do this in one place to ensure that
  ! all subsequent subroutines can assume that if cloud_fraction > 0.0
  ! then cloud is really present and should be treated.
  subroutine crop_cloud_fraction(this, istartcol, iendcol, &
       &    cloud_fraction_threshold, cloud_mixing_ratio_threshold)
    
    use ecradhook, only : lhook, dr_hook, jphook

    class(cloud_type), intent(inout) :: this
    integer,           intent(in)    :: istartcol, iendcol

    integer :: nlev
    integer :: jcol, jlev, jh

    real(jprb) :: cloud_fraction_threshold, cloud_mixing_ratio_threshold
    real(jprb) :: sum_mixing_ratio(istartcol:iendcol)

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloud:crop_cloud_fraction',0,hook_handle)

    nlev = size(this%fraction,2)

    !$acc parallel default(none) present(this) create(sum_mixing_ratio) async(1)
    !$acc loop seq
    do jlev = 1,nlev
      !$acc loop gang(static:1) vector
      do jcol = istartcol,iendcol
        sum_mixing_ratio(jcol) = 0.0_jprb
      end do
      !$acc loop seq
      do jh = 1, this%ntype
        !$acc loop gang(static:1) vector
        do jcol = istartcol,iendcol
          sum_mixing_ratio(jcol) = sum_mixing_ratio(jcol) + this%mixing_ratio(jcol,jlev,jh)
        end do
      end do
      !$acc loop gang(static:1) vector
      do jcol = istartcol,iendcol
        if (this%fraction(jcol,jlev)        < cloud_fraction_threshold &
             &  .or. sum_mixing_ratio(jcol) < cloud_mixing_ratio_threshold) then
          this%fraction(jcol,jlev) = 0.0_jprb
        end if
      end do
    end do
    !$acc end parallel

    if (lhook) call dr_hook('radiation_cloud:crop_cloud_fraction',1,hook_handle)

  end subroutine crop_cloud_fraction


  !---------------------------------------------------------------------
  ! return .true. if variables are out of a physically sensible range,
  ! optionally only considering columns between istartcol and iendcol
  function out_of_physical_bounds(this, istartcol, iendcol, do_fix) result(is_bad)

    use ecradhook,          only : lhook, dr_hook, jphook
    use radiation_check, only : out_of_bounds_2d, out_of_bounds_3d

    class(cloud_type), intent(inout) :: this
    integer,  optional,intent(in) :: istartcol, iendcol
    logical,  optional,intent(in) :: do_fix
    logical                       :: is_bad

    logical    :: do_fix_local

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloud:out_of_physical_bounds',0,hook_handle)

    if (present(do_fix)) then
      do_fix_local = do_fix
    else
      do_fix_local = .false.
    end if

    is_bad =    out_of_bounds_3d(this%mixing_ratio, 'cloud%mixing_ratio', 0.0_jprb, 1.0_jprb, &
         &                       do_fix_local, i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_3d(this%effective_radius, 'cloud%effective_radius', 0.0_jprb, 0.1_jprb, &
         &                       do_fix_local, i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%fraction, 'cloud%fraction', 0.0_jprb, 1.0_jprb, &
         &                       do_fix_local, i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%fractional_std, 'fractional_std', 0.0_jprb, 10.0_jprb, &
         &                       do_fix_local, i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%inv_cloud_effective_size, 'inv_cloud_effective_size', &
         &                       0.0_jprb, 1.0_jprb, do_fix_local, i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%inv_inhom_effective_size, 'inv_inhom_effective_size', &
         &                       0.0_jprb, 1.0_jprb, do_fix_local, i1=istartcol, i2=iendcol) &
         & .or. out_of_bounds_2d(this%overlap_param, 'overlap_param', -0.5_jprb, 1.0_jprb, &
         &                       do_fix_local, i1=istartcol, i2=iendcol)

    if (lhook) call dr_hook('radiation_cloud:out_of_physical_bounds',1,hook_handle)

  end function out_of_physical_bounds

! # 909 "radiation/radiation_cloud.f90"

end module radiation_cloud
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

