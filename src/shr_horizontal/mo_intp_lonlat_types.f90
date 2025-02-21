! Contains the the interpolation data structures for the regular grid.
!
! A hierarchy of data structures is needed since one can use
!
! * several lon-lat grids, where each of them is applied to
!   * several ICON domains, where each of them requires
!     * several interpolation methods for the different variables.
!
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE mo_intp_lonlat_types

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlvert, max_dom, SUCCESS,     &
    &                               HINTP_TYPE_NONE, HINTP_TYPE_LONLAT_RBF,                   &
    &                               HINTP_TYPE_LONLAT_NNB, HINTP_TYPE_LONLAT_BCTR
  USE mo_math_types,          ONLY: t_cartesian_coordinates, t_geographical_coordinates
  USE mo_math_utilities,      ONLY: gc2cc
  USE mo_lonlat_grid,         ONLY: t_lon_lat_grid, OPERATOR(==)
  USE mo_communication,       ONLY: t_comm_gather_pattern
  USE mo_interpol_config,     ONLY: rbf_vec_dim_c, rbf_dim_c2l, l_mono_c2l
  USE mo_model_domain,        ONLY: t_patch
  USE mo_communication,       ONLY: idx_1d




  IMPLICIT NONE


  PUBLIC

  !> module name  
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_intp_lonlat_types'

  !> level of output verbosity (for debugging purposes)
  INTEGER, PARAMETER  :: dbg_level = 0
  

  !> General set of interpolation weights and stencil indices.
  !
  TYPE :: t_intp_coeff
  
    INTEGER,  ALLOCATABLE  :: stencil(:,:)      ! array defining number of entries in stencil
  
    INTEGER,  ALLOCATABLE  :: idx(:,:,:)        ! index array defining the
                                                ! stencil (stencilsize,nproma,nblks_lonlat)
  
    INTEGER,  ALLOCATABLE  :: blk(:,:,:)        ! ... dito for the blocks
  
  CONTAINS
    PROCEDURE            t_intp_coeff_init                             !< constructor
    PROCEDURE, PUBLIC :: finalize       => t_intp_coeff_finalize       !< destructor
    PROCEDURE, PUBLIC :: maxstencilsize => t_intp_coeff_maxstencilsize
  
  END TYPE t_intp_coeff
  
  
  !> Scalar interpolation weights.
  !
  TYPE, EXTENDS(t_intp_coeff) :: t_intp_scalar_coeff
    REAL(wp), ALLOCATABLE  :: coeff(:,:,:)      ! array containing interpolation 
                                                ! weights (stencilsize,nproma,nblks_lonlat)

    ! monotonicity can be enforced by demanding that the interpolated
    ! value is not higher or lower than the stencil point values:
    LOGICAL                :: l_cutoff

    ! array containing additional data sites (only in debugging mode)
    REAL(wp),ALLOCATABLE :: v(:,:,:,:)

  CONTAINS
    PROCEDURE, PUBLIC :: init        => t_intp_scalar_coeff_init
    PROCEDURE, PUBLIC :: finalize    => t_intp_scalar_coeff_finalize   !< destructor
    PROCEDURE         :: t_intp_scalar_interpolate_i
    PROCEDURE         :: t_intp_scalar_interpolate_r
    GENERIC,   PUBLIC :: interpolate => t_intp_scalar_interpolate_r, t_intp_scalar_interpolate_i

    PROCEDURE, PUBLIC :: visualize_stencil => t_intp_scalar_visualize
  
  END TYPE t_intp_scalar_coeff

  
  !> Vector interpolation weights.
  !
  TYPE, EXTENDS(t_intp_coeff) :: t_intp_vec_coeff
    REAL(wp), ALLOCATABLE  :: coeff(:,:,:,:)      ! array containing interpolation 
                                                  ! weights (stencilsize,2,nproma,nblks_lonlat)
  
  CONTAINS
    PROCEDURE, PUBLIC :: init     => t_intp_vec_coeff_init
    PROCEDURE, PUBLIC :: finalize => t_intp_vec_coeff_finalize   !< destructor
    PROCEDURE, PUBLIC :: interpolate => t_intp_vec_interpolate_r
  
  END TYPE t_intp_vec_coeff
  
  
  !> Data structure containing coefficients for (optional)
  !  interpolation onto a particular lon-lat grid.
  !
  !  This data structure is needed since there are multiple
  !  interpolation methods available.
  !
  TYPE t_lon_lat_intp
  
    ! --- Radial Basis Function (RBF) interpolation
    TYPE(t_intp_vec_coeff) :: rbf_vec
  
    ! --- nearest-neighbor interpolation
    !
    !     (we need interpolation weights also here to blank out points)
  
    TYPE(t_intp_scalar_coeff) :: nnb
  
    ! --- barycentric interpolation
  
    TYPE(t_intp_scalar_coeff) :: baryctr ! stencil and weights for barycentric interpolation
  
    ! --- direct RBF interpolation from cell centers to lon-lat points:
  
    TYPE(t_intp_scalar_coeff) :: rbf_c2l
  
    ! --- other data fields
  
    ! coordinates of the lon-lat points (nproma,nb nblks_lonlat)
    TYPE(t_geographical_coordinates), ALLOCATABLE :: ll_coord(:,:)
  
    INTEGER               :: nthis_local_pts            ! number of points local to this PE
    REAL(wp)              :: rbf_scale                  ! RBF shape parameter
    LOGICAL               :: l_initialized = .FALSE.
  
    ! data field for distributed computations (available on all PEs)
    INTEGER, ALLOCATABLE  :: global_idx(:)    ! for each lon-lat point on this PE: global idx
  
  CONTAINS
    PROCEDURE, PUBLIC :: init          => t_lon_lat_intp_init
    PROCEDURE, PUBLIC :: finalize      => t_lon_lat_intp_finalize
    PROCEDURE, PUBLIC :: nblks_lonlat  => t_lon_lat_intp_nblks
    PROCEDURE, PUBLIC :: npromz_lonlat => t_lon_lat_intp_npromz
    PROCEDURE, PUBLIC :: contract      => t_lon_lat_intp_contract
    PROCEDURE            t_lon_lat_intp_interpolate_rvec
    PROCEDURE            t_lon_lat_intp_interpolate_r
    PROCEDURE            t_lon_lat_intp_interpolate_i
    GENERIC,   PUBLIC :: interpolate   => t_lon_lat_intp_interpolate_rvec, &
      &                                   t_lon_lat_intp_interpolate_r,    &
      &                                   t_lon_lat_intp_interpolate_i
  
  END TYPE t_lon_lat_intp
  
  
  !> Collects all interpolation coefficients together with the
  !  corresponding lon-lat grid and a communication pattern for
  !  gathering.
  !
  !  This data type is required since a single lon-lat grid needs
  !  interpolation data for each ICON domain.
  !
  TYPE t_lon_lat_data
    TYPE(t_lon_lat_grid)                 :: grid
    LOGICAL                              :: l_dom (max_dom)
    TYPE(t_lon_lat_intp)                 :: intp  (max_dom)
    TYPE(t_comm_gather_pattern)          :: p_pat (max_dom)

  CONTAINS
    PROCEDURE, PUBLIC :: init     => t_lon_lat_data_init      !< constructor
    PROCEDURE, PUBLIC :: finalize => t_lon_lat_data_finalize  !< destructor

  END TYPE t_lon_lat_data

  
  !> Global list of lon-lat grids and interpolation coefficients.
  !
  !  All lon-lat grids needed for output are stored in this
  !  registry. Several subroutines collaborate by this data structure:
  !  coefficient computation, interpolation algorithm, I/O.
  !
  !  @note On the compute PEs information is *local*, i.e. only for
  !        each lon-lat grid only those parts "owned" by the process
  !        are stored.
  !
  TYPE t_lon_lat_list  
 
    !> Global list of lon-lat grids and interpolation coefficients.
    TYPE (t_lon_lat_data), ALLOCATABLE :: list(:)
  
    !> Actual no. of lon-lat grids currently used in this model
    INTEGER               :: ngrids

  CONTAINS
    PROCEDURE, PUBLIC :: init          => t_lon_lat_list_init           !< constructor
    PROCEDURE, PUBLIC :: finalize      => t_lon_lat_list_finalize       !< destructor
    PROCEDURE, PUBLIC :: get_ID        => t_lon_lat_list_get_ID
    PROCEDURE, PUBLIC :: add_new_grid  => t_lon_lat_list_add_new_grid
  
  END TYPE t_lon_lat_list


  ! MODULE VARIABLES --------------------------------------------------------------

  TYPE (t_lon_lat_list), TARGET, SAVE  :: lonlat_grids
  
CONTAINS

  !--------------------------------------------------------------------------------
  !> Constructor for general interpolation coefficient data structure.
  !
  SUBROUTINE t_intp_coeff_init(this, stencilsize, nproma, nblks)
    CLASS(t_intp_coeff), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: stencilsize, nproma, nblks
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::t_intp_coeff_init"
    INTEGER :: ierr

    ALLOCATE (this%stencil(             nproma, nblks), &
        &     this%idx(    stencilsize, nproma, nblks), &
        &     this%blk(    stencilsize, nproma, nblks), &
        &     STAT=ierr )
    IF (ierr /= SUCCESS)  CALL finish (routine, 'Allocation for coeffs failed!')
    !$ACC ENTER DATA CREATE(this%stencil, this%idx, this%blk)

    this%stencil(:,:) = stencilsize
    this%idx(:,:,:)   = -1
    this%blk(:,:,:)   = -1

    !$ACC UPDATE DEVICE(this%stencil) ASYNC(1)
    !$ACC UPDATE DEVICE(this%idx) ASYNC(1)
    !$ACC UPDATE DEVICE(this%blk) ASYNC(1)

  END SUBROUTINE t_intp_coeff_init


  !--------------------------------------------------------------------------------
  !> Destructor for general interpolation coefficient data structure.
  !
  SUBROUTINE t_intp_coeff_finalize(this)
    CLASS(t_intp_coeff), INTENT(INOUT) :: this
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::t_intp_coeff_finalize"
    INTEGER :: ist

    !$ACC WAIT(1)
    IF (ALLOCATED(this%idx)) THEN
      DEALLOCATE (this%idx, STAT=ist)
      IF (ist /= SUCCESS)  CALL finish (routine, 'deallocation for barycentric stencil indices failed')
      !$ACC EXIT DATA DELETE(this%idx)
    END IF
    IF (ALLOCATED(this%blk)) THEN
      DEALLOCATE (this%blk, STAT=ist)
      IF (ist /= SUCCESS)  CALL finish (routine, 'deallocation for barycentric stencil indices failed')
      !$ACC EXIT DATA DELETE(this%blk)
    END IF
    IF (ALLOCATED(this%stencil)) THEN
      DEALLOCATE (this%stencil, STAT=ist)
      IF (ist /= SUCCESS)  CALL finish (routine, 'deallocation for barycentric stencil indices failed')
      !$ACC EXIT DATA DELETE(this%stencil)
    END IF
  END SUBROUTINE t_intp_coeff_finalize


  !--------------------------------------------------------------------------------
  !> @return Stencil size of this interpolation.
  !
  FUNCTION t_intp_coeff_maxstencilsize(this)
    CLASS(t_intp_coeff), INTENT(IN) :: this
    INTEGER :: t_intp_coeff_maxstencilsize
    IF (.NOT. ALLOCATED(this%idx)) THEN
      t_intp_coeff_maxstencilsize = 0
    ELSE
      t_intp_coeff_maxstencilsize = SIZE(this%idx,1)
    END IF
  END FUNCTION t_intp_coeff_maxstencilsize


  !--------------------------------------------------------------------------------
  !> Constructor for general interpolation coefficient data structure.
  !
  SUBROUTINE t_intp_scalar_coeff_init(this, stencilsize, nproma, nblks, l_cutoff)
    CLASS(t_intp_scalar_coeff), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: stencilsize, nproma, nblks
    LOGICAL, INTENT(IN) :: l_cutoff
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::t_intp_scalar_coeff_init"
    INTEGER :: ierr

    CALL t_intp_coeff_init(this, stencilsize, nproma, nblks)

    ALLOCATE (this%coeff(stencilsize, nproma, nblks), STAT=ierr )
    IF (ierr /= SUCCESS)  CALL finish (routine, 'Allocation for coeffs failed!')
    !$ACC ENTER DATA CREATE(this%coeff)
    this%coeff(:,:,:) = 0._wp
    !$ACC UPDATE DEVICE(this%coeff) ASYNC(1)

    this%l_cutoff = l_cutoff

    IF (dbg_level > 5) THEN
      ALLOCATE (this%v(3,stencilsize, nproma, nblks), STAT=ierr )
      IF (ierr /= SUCCESS)  CALL finish (routine, 'Allocation for coeffs failed!')
      !$ACC ENTER DATA CREATE(this%v)
      this%v(:,:,:,:) = 0._wp
      !$ACC UPDATE DEVICE(this%v) ASYNC(1)
    END IF
  END SUBROUTINE t_intp_scalar_coeff_init


  !--------------------------------------------------------------------------------
  !> Destructor for general interpolation coefficient data structure.
  !
  SUBROUTINE t_intp_scalar_coeff_finalize(this)
    CLASS(t_intp_scalar_coeff), INTENT(INOUT) :: this
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::t_intp_scalar_coeff_finalize"
    INTEGER :: ist

    CALL t_intp_coeff_finalize(this)
    IF (ALLOCATED(this%coeff)) THEN
      DEALLOCATE (this%coeff, STAT=ist)
      IF (ist /= SUCCESS) CALL finish (routine, 'deallocation for barycentric lon-lat coefficients failed')
      !$ACC EXIT DATA DELETE(this%coeff)
    END IF
    IF (ALLOCATED(this%v)) THEN
      DEALLOCATE (this%v, STAT=ist)
      IF (ist /= SUCCESS) CALL finish (routine, 'deallocation for barycentric lon-lat coefficients failed')
      !$ACC EXIT DATA DELETE(this%v)
    END IF
  END SUBROUTINE t_intp_scalar_coeff_finalize


  !--------------------------------------------------------------------------------
  !> Constructor for general interpolation coefficient data structure.
  !
  SUBROUTINE t_intp_vec_coeff_init(this, stencilsize, nproma, nblks)
    CLASS(t_intp_vec_coeff), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: stencilsize, nproma, nblks
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::t_intp_vec_coeff_init"
    INTEGER :: ierr

    CALL t_intp_coeff_init(this, stencilsize, nproma, nblks)

    ALLOCATE (this%coeff(stencilsize, 2, nproma, nblks), STAT=ierr )
    IF (ierr /= SUCCESS)  CALL finish (routine, 'Allocation for coeffs failed!')
    !$ACC ENTER DATA CREATE(this%coeff)
    this%coeff(:,:,:,:) = 0._wp
    !$ACC UPDATE DEVICE(this%coeff) ASYNC(1)
  END SUBROUTINE t_intp_vec_coeff_init


  !--------------------------------------------------------------------------------
  !> Destructor for general interpolation coefficient data structure.
  !
  SUBROUTINE t_intp_vec_coeff_finalize(this)
    CLASS(t_intp_vec_coeff), INTENT(INOUT) :: this
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::t_intp_vec_coeff_finalize"
    INTEGER :: ist

    CALL t_intp_coeff_finalize(this)
    IF (ALLOCATED(this%coeff)) THEN
      DEALLOCATE (this%coeff, STAT=ist)
      IF (ist /= SUCCESS) CALL finish (routine, 'deallocation for barycentric lon-lat coefficients failed')
      !$ACC EXIT DATA DELETE(this%coeff)
    END IF
  END SUBROUTINE t_intp_vec_coeff_finalize


  !--------------------------------------------------------------------------------
  !> Constructor for lon-lat interpolation data structure.
  !
  SUBROUTINE t_lon_lat_intp_init(this, nproma)
    CLASS(t_lon_lat_intp), INTENT(INOUT) :: this    
    INTEGER,               INTENT(IN)    :: nproma
    INTEGER :: nblks_lonlat

    nblks_lonlat = this%nblks_lonlat(nproma)

    ! --- Radial Basis Function (RBF) interpolation
    CALL this%rbf_vec%init(rbf_vec_dim_c, nproma, nblks_lonlat)
    
    ! --- direct RBF interpolation from cell centers to lon-lat points:
    CALL this%rbf_c2l%init(rbf_dim_c2l, nproma, nblks_lonlat, l_mono_c2l)
    
    ! --- nearest-neighbor interpolation
    !
    !     (we need interpolation weights also here to blank out points)
    CALL this%nnb%init(1, nproma, nblks_lonlat, .FALSE.)
    
    ! --- barycentric interpolation
    CALL this%baryctr%init(3, nproma, nblks_lonlat, .FALSE.)

    ! global_idx is needed in device
    !$ACC ENTER DATA CREATE(this%global_idx)
  END SUBROUTINE t_lon_lat_intp_init


  !--------------------------------------------------------------------------------
  !> Destructor for lon-lat interpolation data structure.
  !
  SUBROUTINE t_lon_lat_intp_finalize(this)
    CLASS(t_lon_lat_intp), INTENT(INOUT) :: this    
    
    CHARACTER(*), PARAMETER :: routine = modname//"::t_lon_lat_intp_finalize"
    INTEGER :: ist

    IF (.NOT. this%l_initialized) RETURN

    DEALLOCATE (this%global_idx, STAT=ist )
    IF (ist /= SUCCESS)  CALL finish (routine, 'deallocation for lon-lat coefficients failed')
    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(this%global_idx)
    ! -- deallocate rbf_vec data structure
    CALL this%rbf_vec%finalize()
    
    ! -- deallocate rbf_c2l data structure
    CALL this%rbf_c2l%finalize()
    
    ! -- deallocate nnb data structure
    CALL this%nnb%finalize()
    
    ! -- deallocate baryctr data structure
    CALL this%baryctr%finalize()
    
    ! -- array with lon-lat coordinates
    DEALLOCATE(this%ll_coord, stat=ist)
    IF (ist /= SUCCESS)  CALL finish (routine, 'Deallocation of array with lon-lat coordinates!')
    !$ACC EXIT DATA DELETE(this%ll_coord)

    this%l_initialized = .FALSE.
  END SUBROUTINE t_lon_lat_intp_finalize


  !--------------------------------------------------------------------------------
  !> Destructor for general interpolation coefficient data structure.
  !
  FUNCTION t_lon_lat_intp_nblks(this, nproma)
    INTEGER :: t_lon_lat_intp_nblks
    CLASS(t_lon_lat_intp), INTENT(IN) :: this
    INTEGER,               INTENT(IN) :: nproma

    t_lon_lat_intp_nblks  = (this%nthis_local_pts - 1)/nproma + 1
  END FUNCTION t_lon_lat_intp_nblks


  !--------------------------------------------------------------------------------
  !> Destructor for general interpolation coefficient data structure.
  !
  FUNCTION t_lon_lat_intp_npromz(this, nproma)
    INTEGER :: t_lon_lat_intp_npromz
    CLASS(t_lon_lat_intp), INTENT(IN) :: this
    INTEGER,               INTENT(IN) :: nproma
    t_lon_lat_intp_npromz = this%nthis_local_pts - (this%nblks_lonlat(nproma)-1)*nproma
  END FUNCTION t_lon_lat_intp_npromz


  !--------------------------------------------------------------------------------
  !> Resize arrays with lon-lat interpolation data from global to
  !  local size. This can be performed as soon as the setup phase is
  !  finished.
  SUBROUTINE t_lon_lat_intp_contract(this)
    CLASS(t_lon_lat_intp), INTENT(INOUT) :: this
    
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::t_lon_lat_intp_contract"
    INTEGER, ALLOCATABLE  :: tmp_global_idx(:)
    INTEGER               :: errstat, nlocal_pts
    
    nlocal_pts  = this%nthis_local_pts
    ! first allocate temporary storage and copy fields:
    ALLOCATE(tmp_global_idx(nlocal_pts), STAT=errstat )
    IF (errstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed')
    !$ACC ENTER DATA CREATE(tmp_global_idx)
    tmp_global_idx(1:nlocal_pts) = this%global_idx(1:nlocal_pts)
    !$ACC UPDATE DEVICE(tmp_global_idx) ASYNC(1)
! Note: IF_PRESENT not allowed in EXIT DATA, but global_idx should be on DEVICE
    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(this%global_idx)
    CALL MOVE_ALLOC(tmp_global_idx, this%global_idx)
  END SUBROUTINE t_lon_lat_intp_contract


  !--------------------------------------------------------------------------------
  !> Constructor: Allocate lon-lat interpolation data structure
  !
  SUBROUTINE t_lon_lat_data_init(this, grid, rbf_scale)
    CLASS(t_lon_lat_data), INTENT(INOUT) :: this
    TYPE (t_lon_lat_grid), INTENT(IN)    :: grid
    REAL(wp),              INTENT(IN)    :: rbf_scale

    this%l_dom(:)         = .FALSE.
    this%grid             = grid
    this%intp%rbf_scale   = rbf_scale
  END SUBROUTINE t_lon_lat_data_init


  !--------------------------------------------------------------------------------
  !> Destructor: Clean-up of interpolation data structure
  !
  SUBROUTINE t_lon_lat_data_finalize(this)
    CLASS(t_lon_lat_data), INTENT(INOUT) :: this
    ! local variables
    INTEGER :: jg
    
    DO jg=1,SIZE(this%intp)
      CALL this%intp(jg)%finalize()
    END DO
  END SUBROUTINE t_lon_lat_data_finalize
  

  !--------------------------------------------------------------------------------
  !> Constructor: Setup of lon-lat registry
  !
  SUBROUTINE t_lon_lat_list_init(this)
    CLASS(t_lon_lat_list), INTENT(INOUT) :: this
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::t_lon_lat_list_init"
    INTEGER,      PARAMETER :: INITIAL_SIZE = 20
    INTEGER :: ist
    
    ! not much to do yet...
    this%ngrids = 0
    ALLOCATE(this%list(INITIAL_SIZE), STAT=ist)  ! Lists are not on the device
    IF (ist /= SUCCESS)  CALL finish (routine, 'Allocate failed!')
  END SUBROUTINE t_lon_lat_list_init


  !--------------------------------------------------------------------------------
  !> Destructor: Frees all data allocated by lon-lat grid list
  !
  SUBROUTINE t_lon_lat_list_finalize(this)
    CLASS(t_lon_lat_list), INTENT(INOUT) :: this
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::t_lon_lat_list_finalize"
    INTEGER :: i, ist
    
    DO i=1, this%ngrids
      CALL this%list(i)%finalize()
    END DO
    DEALLOCATE(this%list, STAT=ist)              ! Lists are not on the device
    IF (ist /= SUCCESS)  CALL finish (routine, 'Deallocate failed!')
    this%ngrids = 0
  END SUBROUTINE t_lon_lat_list_finalize


  !--------------------------------------------------------------------------------
  ! @return index in lon-lat grid list for "grid", -1 if unknown.
  !
  FUNCTION t_lon_lat_list_get_ID(this, grid)
    INTEGER :: t_lon_lat_list_get_ID
    CLASS(t_lon_lat_list), INTENT(IN)    :: this
    TYPE(t_lon_lat_grid),  INTENT(IN)    :: grid
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::t_lon_lat_list_get_ID"
    INTEGER :: i

    t_lon_lat_list_get_ID = -1 ! default value: "not found"
    IF (.NOT. ALLOCATED(this%list)) RETURN

    DO i=1,this%ngrids
      IF (this%list(i)%grid == grid) THEN
        t_lon_lat_list_get_ID = i
        EXIT
      END IF
    END DO
  END FUNCTION t_lon_lat_list_get_ID


  !--------------------------------------------------------------------------------
  !> @return index with "free" lon-lat grid data
  !
  FUNCTION t_lon_lat_list_add_new_grid(this)
    INTEGER :: t_lon_lat_list_add_new_grid
    CLASS(t_lon_lat_list), INTENT(INOUT) :: this
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::t_lon_lat_list_add_new_grid"
    TYPE(t_lon_lat_data), ALLOCATABLE :: temp(:)
    INTEGER :: new_size

    IF (this%ngrids == SIZE(this%list)) THEN
      ! Pre-allocated number of lon-lat grids exceeded. We need to
      ! resize the grid data structure:
      new_size = 2*SIZE(this%list)
      ALLOCATE(temp(new_size))
      temp(1:SIZE(this%list)) = this%list(:)
      CALL move_alloc(temp, this%list)
    END IF
    this%ngrids = this%ngrids + 1
    t_lon_lat_list_add_new_grid = this%ngrids
  END FUNCTION t_lon_lat_list_add_new_grid


  !--------------------------------------------------------------------------------
  !> Performs nearest neighbor interpolation, INTEGER implementation
  !
  SUBROUTINE t_intp_scalar_interpolate_i( this, p_cell_in, nproma, nblks_lonlat, npromz_lonlat, &
    &                                     p_out, opt_slev, opt_elev)

    ! Indices of source points and interpolation coefficients
    CLASS(t_intp_scalar_coeff), TARGET, INTENT(IN)           :: this
    ! input cell-based variable for which gradient at cell center is computed
    INTEGER,                            INTENT(IN)           :: p_cell_in(:,:,:) ! dim: (nproma,nlev,nblks)
    INTEGER,                            INTENT(IN)           :: nproma, nblks_lonlat, npromz_lonlat
    ! reconstructed scalar value at lon-lat point, dim: (nproma,nlev,nblks_lonlat)
    INTEGER,                            INTENT(INOUT)        :: p_out(:,:,:)
    ! optional vertical start/end level
    INTEGER,                            INTENT(IN), OPTIONAL :: opt_slev, opt_elev

    ! LOCAL VARIABLES
    CHARACTER(*), PARAMETER :: routine = modname//"::t_intp_scalar_interpolate_i"
    INTEGER :: slev, elev,                 & ! vertical start and end level
      &        i_startidx, i_endidx,       & ! start/end index
      &        jc, jb, jk                    ! integer over lon-lat points, levels
    INTEGER,  DIMENSION(:,:,:), POINTER :: iidx, iblk
    REAL(wp), DIMENSION(:,:,:), POINTER :: ptr_coeff

    INTEGER :: stencilsize   ! lon-lat grid blocking info
    LOGICAL :: l_cutoff

    slev = 1
    elev = UBOUND(p_cell_in,2)
    ! check optional arguments
    IF ( PRESENT(opt_slev) ) slev = opt_slev
    IF ( PRESENT(opt_elev) ) elev = opt_elev

    iidx        => this%idx
    iblk        => this%blk
    ptr_coeff   => this%coeff
    stencilsize =  this%maxstencilsize()
    l_cutoff    =  this%l_cutoff

    ! consistency check
    IF (stencilsize /= 1) THEN
      CALL finish(routine, "Not implemented!")
    END IF
    ! consistency check
    IF (l_cutoff) THEN
      CALL finish(routine, "Not implemented!")
    END IF

    !$ACC DATA PRESENT(p_cell_in, ptr_coeff, p_out, iidx, iblk) IF(i_am_accel_node)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc), SCHEDULE(runtime)
    DO jb = 1,nblks_lonlat
      i_startidx = 1
      i_endidx   = nproma
      IF (jb == nblks_lonlat) i_endidx = npromz_lonlat
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)




      DO jk = slev, elev
        DO jc = i_startidx, i_endidx

          p_out(jc,jk,jb) = INT(ptr_coeff(1,jc,jb)) *  &
            & p_cell_in(iidx(1,jc,jb), jk, iblk(1,jc,jb))
        ENDDO
      ENDDO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    !$ACC END DATA

  END SUBROUTINE t_intp_scalar_interpolate_i
   

  !--------------------------------------------------------------------------------
  !> Performs vector RBF reconstruction at lon-lat grid points.
  !  This routine takes a cell-based variable as input and reconstructs
  !  the gradient (of the FD approximation) at the lon-lat grid points.
  !
  !  More precisely, the routine "rbf_interpol_c2grad_lonlat" combines
  !  the computation of gradients from scalar variables,
  !  CALL grad_fd_norm( p_cell_in(:,:,:), &
  !    &                ptr_patch, grad_norm_psi_e )
  !
  !  and the application of interpolation coefficients,
  !  CALL rbf_vec_interpol_lonlat( grad_norm_psi_e, ptr_patch, ptr_int, grad_x,  &
  !    &                           grad_y, ptr_int%tri_idx(:,:,:),        &
  !    &                           nblks_lonlat, npromz_lonlat )
  !
  SUBROUTINE t_intp_scalar_interpolate_r( this, p_cell_in, nproma, nblks_lonlat, npromz_lonlat, &
    &                                     p_out, opt_slev, opt_elev)

    ! Indices of source points and interpolation coefficients
    CLASS(t_intp_scalar_coeff), TARGET, INTENT(IN)           :: this
    ! input triangle-cell-based variable for which gradient at cell center is computed
    REAL(wp),                           INTENT(IN)           :: p_cell_in(:,:,:)
    INTEGER,                            INTENT(IN)           :: nproma, nblks_lonlat, npromz_lonlat
    ! reconstructed scalar value at lon-lat point, dim: (nproma,nlev,nblks_lonlat)
    REAL(wp),                           INTENT(INOUT)        :: p_out(:,:,:)     
    ! optional vertical start/end level
    INTEGER,                            INTENT(IN), OPTIONAL :: opt_slev, opt_elev

    ! Local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::t_intp_scalar_interpolate_r"
    INTEGER                             :: slev, elev,               &  ! vertical start and end level
      &                                    jc, jb, jk,               &  ! integer over lon-lat points, levels
      &                                    i_startidx, i_endidx         ! start/end index
    REAL(wp)                            :: vmin, vmax
    INTEGER,  DIMENSION(:,:,:), POINTER :: iidx, iblk
    REAL(wp), DIMENSION(:,:,:), POINTER :: ptr_coeff
    INTEGER                             :: stencilsize   ! lon-lat grid blocking info
    LOGICAL                             :: l_cutoff
    INTEGER                             :: ji  ! For the min-max loops
    REAL(wp)                            :: tmp ! For the min-max loops

    slev = 1
    elev = UBOUND(p_cell_in,2)
    ! check optional arguments
    IF ( PRESENT(opt_slev) ) slev = opt_slev
    IF ( PRESENT(opt_elev) ) elev = opt_elev

    iidx        => this%idx
    iblk        => this%blk
    ptr_coeff   => this%coeff
    stencilsize =  this%maxstencilsize()
    l_cutoff    =  this%l_cutoff

    ! consistency check
    IF (ANY(stencilsize == (/ 1,3 /)) .AND. l_cutoff) THEN
      CALL finish(routine, "Not implemented!")
    END IF

    !$ACC DATA PRESENT(p_cell_in, ptr_coeff, p_out, iidx, iblk) IF(i_am_accel_node)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc,vmin,vmax), SCHEDULE(runtime)

    DO jb = 1,nblks_lonlat

      i_startidx = 1
      i_endidx   = nproma
      IF (jb == nblks_lonlat) i_endidx = npromz_lonlat

      ! we have to duplicate code here for different stencil sizes,
      ! otherwise we would break vectorization...
      SELECT CASE(stencilsize)

      CASE(1)

        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)




        DO jk = slev, elev
          DO jc = i_startidx, i_endidx


            p_out(jc,jk,jb) = ptr_coeff(1,jc,jb) *  &
              & p_cell_in(iidx(1,jc,jb), jk, iblk(1,jc,jb))

          END DO
        END DO

      CASE(3)

        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)




        DO jk = slev, elev
          DO jc = i_startidx, i_endidx


            p_out(jc,jk,jb) = &
          &    ptr_coeff(1,jc,jb)*                                            &
          &    p_cell_in(iidx(1,jc,jb), jk, iblk(1,jc,jb))   &
          &  + ptr_coeff(2,jc,jb)*                                            &
          &    p_cell_in(iidx(2,jc,jb), jk, iblk(2,jc,jb))   &
          &  + ptr_coeff(3,jc,jb)*                                            &
          &    p_cell_in(iidx(3,jc,jb), jk, iblk(3,jc,jb))

          END DO
        END DO

      CASE(4)

        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)




        DO jk = slev, elev
          DO jc = i_startidx, i_endidx


            p_out(jc,jk,jb) =                                                   &
              ptr_coeff(1 ,jc,jb)*p_cell_in(iidx(1 ,jc,jb),jk,iblk(1 ,jc,jb)) + &
              ptr_coeff(2 ,jc,jb)*p_cell_in(iidx(2 ,jc,jb),jk,iblk(2 ,jc,jb)) + &
              ptr_coeff(3 ,jc,jb)*p_cell_in(iidx(3 ,jc,jb),jk,iblk(3 ,jc,jb)) + &
              ptr_coeff(4 ,jc,jb)*p_cell_in(iidx(4 ,jc,jb),jk,iblk(4 ,jc,jb))

            ! monotonicity can be enforced by demanding that the interpolated
            ! value is not higher or lower than the stencil point values.

            ! Cf. the "lmono" implementation in the GME:
            ! D. Majewski, "Documentation of the new global model (GME)
            !               of the DWD" (1996)
            IF (l_cutoff) THEN
              vmin = p_cell_in(iidx(1 ,jc,jb),jk,iblk(1 ,jc,jb))
              vmax = vmin
!NEC$ unroll_complete
              DO ji = 2, 4
                tmp = p_cell_in(iidx(ji,jc,jb),jk,iblk(ji,jc,jb))
                vmin = MIN(vmin, tmp)
                vmax = MAX(vmax, tmp)
              END DO

              p_out(jc,jk,jb) = MAX( MIN(p_out(jc,jk,jb), vmax), vmin )
            END IF

          ENDDO
        ENDDO

      CASE(10)

        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)




        DO jk = slev, elev
          DO jc = i_startidx, i_endidx


            p_out(jc,jk,jb) =                                                   &
              ptr_coeff(1 ,jc,jb)*p_cell_in(iidx(1 ,jc,jb),jk,iblk(1 ,jc,jb)) + &
              ptr_coeff(2 ,jc,jb)*p_cell_in(iidx(2 ,jc,jb),jk,iblk(2 ,jc,jb)) + &
              ptr_coeff(3 ,jc,jb)*p_cell_in(iidx(3 ,jc,jb),jk,iblk(3 ,jc,jb)) + &
              ptr_coeff(4 ,jc,jb)*p_cell_in(iidx(4 ,jc,jb),jk,iblk(4 ,jc,jb)) + &
              ptr_coeff(5 ,jc,jb)*p_cell_in(iidx(5 ,jc,jb),jk,iblk(5 ,jc,jb)) + &
              ptr_coeff(6 ,jc,jb)*p_cell_in(iidx(6 ,jc,jb),jk,iblk(6 ,jc,jb)) + &
              ptr_coeff(7 ,jc,jb)*p_cell_in(iidx(7 ,jc,jb),jk,iblk(7 ,jc,jb)) + &
              ptr_coeff(8 ,jc,jb)*p_cell_in(iidx(8 ,jc,jb),jk,iblk(8 ,jc,jb)) + &
              ptr_coeff(9 ,jc,jb)*p_cell_in(iidx(9 ,jc,jb),jk,iblk(9 ,jc,jb)) + &
              ptr_coeff(10,jc,jb)*p_cell_in(iidx(10,jc,jb),jk,iblk(10,jc,jb))

            ! monotonicity can be enforced by demanding that the interpolated
            ! value is not higher or lower than the stencil point values.

            ! Cf. the "lmono" implementation in the GME:
            ! D. Majewski, "Documentation of the new global model (GME)
            !               of the DWD" (1996)
            IF (l_cutoff) THEN
              vmin = p_cell_in(iidx(1 ,jc,jb),jk,iblk(1 ,jc,jb))
              vmax = vmin
!NEC$ unroll_complete
              DO ji = 2, 10
                tmp = p_cell_in(iidx(ji,jc,jb),jk,iblk(ji,jc,jb))
                vmin = MIN(vmin, tmp)
                vmax = MAX(vmax, tmp)
              END DO

              p_out(jc,jk,jb) = MAX( MIN(p_out(jc,jk,jb), vmax), vmin )
            END IF

          ENDDO
        ENDDO

      CASE(13)

        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)




        DO jk = slev, elev
          DO jc = i_startidx, i_endidx


            p_out(jc,jk,jb) =                                                   &
              ptr_coeff(1 ,jc,jb)*p_cell_in(iidx(1 ,jc,jb),jk,iblk(1 ,jc,jb)) + &
              ptr_coeff(2 ,jc,jb)*p_cell_in(iidx(2 ,jc,jb),jk,iblk(2 ,jc,jb)) + &
              ptr_coeff(3 ,jc,jb)*p_cell_in(iidx(3 ,jc,jb),jk,iblk(3 ,jc,jb)) + &
              ptr_coeff(4 ,jc,jb)*p_cell_in(iidx(4 ,jc,jb),jk,iblk(4 ,jc,jb)) + &
              ptr_coeff(5 ,jc,jb)*p_cell_in(iidx(5 ,jc,jb),jk,iblk(5 ,jc,jb)) + &
              ptr_coeff(6 ,jc,jb)*p_cell_in(iidx(6 ,jc,jb),jk,iblk(6 ,jc,jb)) + &
              ptr_coeff(7 ,jc,jb)*p_cell_in(iidx(7 ,jc,jb),jk,iblk(7 ,jc,jb)) + &
              ptr_coeff(8 ,jc,jb)*p_cell_in(iidx(8 ,jc,jb),jk,iblk(8 ,jc,jb)) + &
              ptr_coeff(9 ,jc,jb)*p_cell_in(iidx(9 ,jc,jb),jk,iblk(9 ,jc,jb)) + &
              ptr_coeff(10,jc,jb)*p_cell_in(iidx(10,jc,jb),jk,iblk(10,jc,jb)) + &
              ptr_coeff(11,jc,jb)*p_cell_in(iidx(11,jc,jb),jk,iblk(11,jc,jb)) + &
              ptr_coeff(12,jc,jb)*p_cell_in(iidx(12,jc,jb),jk,iblk(12,jc,jb)) + &
              ptr_coeff(13,jc,jb)*p_cell_in(iidx(13,jc,jb),jk,iblk(13,jc,jb))

            ! monotonicity can be enforced by demanding that the interpolated
            ! value is not higher or lower than the stencil point values.

            ! Cf. the "lmono" implementation in the GME:
            ! D. Majewski, "Documentation of the new global model (GME)
            !               of the DWD" (1996)
            IF (l_cutoff) THEN
              vmin = p_cell_in(iidx(1 ,jc,jb),jk,iblk(1 ,jc,jb))
              vmax = vmin
!NEC$ unroll_complete
              DO ji = 2, 13
                tmp = p_cell_in(iidx(ji,jc,jb),jk,iblk(ji,jc,jb))
                vmin = MIN(vmin, tmp)
                vmax = MAX(vmax, tmp)
              END DO

              p_out(jc,jk,jb) = MAX( MIN(p_out(jc,jk,jb), vmax), vmin )
            END IF





          ENDDO
        ENDDO


      END SELECT

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    !$ACC END DATA

  END SUBROUTINE t_intp_scalar_interpolate_r


  !--------------------------------------------------------------------------------
  !> Writes the interpolation stencil (and weights) for a particular
  !  destination grid point to an ASCII file in the *.asy format
  !  ("Asymptote" is a widely known vector graphics language, see
  !  http://asymptote.sourceforge.net).
  !
  !  Afterwards, run the command-line
  !   asy -f pdf <filename>
  !  to produce the PDF result. 
  !
  SUBROUTINE t_intp_scalar_visualize(this, point_gc, jc,jb, ptr_patch, filename)
    ! Indices of source points and interpolation coefficients
    CLASS(t_intp_scalar_coeff), TARGET,  INTENT(IN)  :: this
    TYPE(t_geographical_coordinates), INTENT(IN)  :: point_gc
    ! index of lon-lat grid point: point_gc = ll_coord(jc,jb):
    INTEGER,                          INTENT(IN)  :: jc,jb
    TYPE(t_patch),                    INTENT(IN)  :: ptr_patch
    CHARACTER(LEN=*),                 INTENT(IN)  :: filename
    ! local variables
    LOGICAL,  PARAMETER :: pdf_output    = .TRUE.
    REAL(wp), PARAMETER :: scale_factor  = 10000._wp
    INTEGER                          :: out_unit, i, j, cidx,cblk, vidx,vblk, global_idx
    TYPE(t_cartesian_coordinates)    :: camera, ll_point, vv

    out_unit=20
    OPEN (unit=out_unit,file=TRIM(filename),action="write",status="replace")

    ! --- print asy preamble
    WRITE (out_unit,*) "import three;import graph3;" 
    IF (pdf_output) THEN
      WRITE (out_unit,*) 'settings.prc = false; settings.tex = "pdflatex"; settings.render = 0;'
    ELSE
      WRITE (out_unit,*) 'settings.outformat=""; settings.render = -1;'
    END IF
    WRITE (out_unit,*) "size(10cm);size3(50cm,50cm,50cm);scale(true);"
    WRITE (out_unit,*) "defaultpen(0.5 + fontsize(10.0) + Helvetica());"

    ! --- set camera position
    ll_point = gc2cc(point_gc)
    camera   = ll_point
    camera%x(:) = scale_factor * camera%x(:)
    WRITE (out_unit,*) "currentprojection=orthographic(",camera%x(1),",",camera%x(2),",",camera%x(3),");"

    ! --- print out additional vertices (if available):
    IF (ALLOCATED(this%v)) THEN
      WRITE (out_unit,*) "triple[] vpoints;"
      DO i=1,3
        vv%x(:) = this%v(:,i,jc,jb)
        CALL print_point("vpoints.push", vv)
      END DO
      WRITE (out_unit,*) "  draw(vpoints[0]--vpoints[1], red+0.5bp+dashed);"
      WRITE (out_unit,*) "  draw(vpoints[1]--vpoints[2], red+0.5bp+dashed);"
      WRITE (out_unit,*) "  draw(vpoints[2]--vpoints[0], red+0.5bp+dashed);"
    END IF

    ! --- print out vertex list and cell center points
    WRITE (out_unit,*) "triple[] points; triple[] cpoints; string[] wgt; string[] cell_name;"
    CALL print_point("points.push", ll_point)
    DO i=1,this%stencil(jc,jb)
      cidx = this%idx(i,jc,jb)
      cblk = this%blk(i,jc,jb)
      global_idx = idx_1d(cidx,cblk)
      IF (global_idx < 1) THEN
        WRITE (0,*) "cidx,cblk = ", cidx, cblk
      END IF
      WRITE (out_unit,*) 'cell_name.push("', ptr_patch%cells%decomp_info%glb_index(global_idx), '");'
      CALL print_point("cpoints.push", gc2cc(ptr_patch%cells%center(cidx,cblk)))
      WRITE (out_unit,'(a,F10.4,a)') 'wgt.push("',this%coeff(i,jc,jb),'");'
      DO j=1,3
        vidx = ptr_patch%cells%vertex_idx(cidx,cblk,j)
        vblk = ptr_patch%cells%vertex_blk(cidx,cblk,j)
        CALL print_point("points.push", gc2cc(ptr_patch%verts%vertex(vidx,vblk)))
      END DO
    END DO

    ! --- print out triangulation
    WRITE (out_unit,*) "int[][] trn = {"
    DO i=1,this%stencil(jc,jb)
      WRITE (out_unit,*) "{"
      DO j=1,3
        WRITE (out_unit,*) 3*(i-1)+j, ", "
      END DO
      IF (i<this%stencil(jc,jb)) THEN
        WRITE (out_unit,*) "},"
      ELSE
        WRITE (out_unit,*) "}"
      END IF
    END DO
    WRITE (out_unit,*) "};"

    ! --- draw vertices
    WRITE (out_unit,*)  "for(int i=1; i < points.length; ++i) {"
    WRITE (out_unit,*)  "  dot(points[i], black);"
    WRITE (out_unit,*)  "}"
    WRITE (out_unit,*)  "for(int i=0; i < cpoints.length; ++i) {"
    WRITE (out_unit,*)  "  dot(cpoints[i], blue);"
    WRITE (out_unit,*)  "  draw(cpoints[i]--points[0], blue+1.0bp);"
    WRITE (out_unit,*)  '  label(scale(0.5)*minipage("$\textnormal{wgt:} {"+wgt[i]+"}$"),', &
      &                 "        0.5*(cpoints[i]+points[0]),2NE, blue);"
    WRITE (out_unit,*)  '  label(scale(0.5)*minipage("${"+cell_name[i]+"}$"),', &
      &                 "        cpoints[i],2NE, black);"
    WRITE (out_unit,*)  "}"

    ! --- draw triangles
    WRITE (out_unit,*) "for(int i=0; i < trn.length; ++i) {"
    WRITE (out_unit,*) "  draw(points[trn[i][0]]--points[trn[i][1]], black+1.0bp);"
    WRITE (out_unit,*) "  draw(points[trn[i][1]]--points[trn[i][2]], black+1.0bp);"
    WRITE (out_unit,*) "  draw(points[trn[i][2]]--points[trn[i][0]], black+1.0bp);"
    WRITE (out_unit,*) "}"

    ! --- draw destination point
    WRITE (out_unit,*)  "dot(points[0], red);"

    CLOSE(out_unit)

  CONTAINS

    SUBROUTINE print_point(prefix, p)
      CHARACTER(LEN=*), INTENT(IN) :: prefix
      TYPE(t_cartesian_coordinates), INTENT(IN) :: p
      TYPE(t_cartesian_coordinates) :: vv
      vv = p; vv%x(:) = scale_factor * vv%x(:)
      WRITE (out_unit,'(a,3(F10.4,a))') TRIM(prefix)//"((",vv%x(1),",",vv%x(2),",",vv%x(3),"));"
    END SUBROUTINE print_point
  
  END SUBROUTINE t_intp_scalar_visualize


  !--------------------------------------------------------------------------------
  !> Performs vector RBF reconstruction at lon-lat grid points.
  !
  ! This routine is based on mo_intp_rbf::rbf_vec_interpol_cell()
  !
  SUBROUTINE t_intp_vec_interpolate_r( this, p_vn_in, nproma, nblks_lonlat, npromz_lonlat, &
    &                                  grad_x, grad_y, opt_slev, opt_elev)

    ! Indices of source points and interpolation coefficients
    CLASS(t_intp_vec_coeff), TARGET, INTENT(IN)           :: this
    ! input normal components of (velocity) vectors at edge midpoints
    REAL(wp),                        INTENT(IN)           :: p_vn_in(:,:,:)  ! dim: (nproma,nlev,nblks_e)
    INTEGER,                         INTENT(IN)           :: nproma, nblks_lonlat, npromz_lonlat
    ! reconstructed x/y-components of velocity vector, dim: (nproma,nlev,nblks_lonlat)
    REAL(wp),                        INTENT(INOUT)        :: grad_x(:,:,:), grad_y(:,:,:)
    ! optional vertical start/end level:
    INTEGER,                         INTENT(in), OPTIONAL :: opt_slev, opt_elev

    ! LOCAL VARIABLES
    INTEGER :: slev, elev,                 & ! vertical start and end level
      &        i_startidx, i_endidx,       & ! start/end index
      &        jc, jb, jk                    ! integer over lon-lat points, levels

    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk
    REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff

    slev = 1
    elev = UBOUND(p_vn_in,2)
    ! check optional arguments
    IF ( PRESENT(opt_slev) ) slev = opt_slev
    IF ( PRESENT(opt_elev) ) elev = opt_elev

    iidx      => this%idx
    iblk      => this%blk
    ptr_coeff => this%coeff

    !$ACC DATA PRESENT(p_vn_in, ptr_coeff, grad_x, grad_y, iidx, iblk) IF(i_am_accel_node)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc), SCHEDULE(runtime)
    DO jb = 1,nblks_lonlat

      i_startidx = 1
      i_endidx   = nproma
      IF (jb == nblks_lonlat) i_endidx = npromz_lonlat

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) COLLAPSE(2) ASYNC(1) IF(i_am_accel_node)




      DO jk = slev, elev
        DO jc = i_startidx, i_endidx


          grad_x(jc,jk,jb) =                                               &
            ptr_coeff(1,1,jc,jb)*p_vn_in(iidx(1,jc,jb),jk,iblk(1,jc,jb)) + &
            ptr_coeff(2,1,jc,jb)*p_vn_in(iidx(2,jc,jb),jk,iblk(2,jc,jb)) + &
            ptr_coeff(3,1,jc,jb)*p_vn_in(iidx(3,jc,jb),jk,iblk(3,jc,jb)) + &
            ptr_coeff(4,1,jc,jb)*p_vn_in(iidx(4,jc,jb),jk,iblk(4,jc,jb)) + &
            ptr_coeff(5,1,jc,jb)*p_vn_in(iidx(5,jc,jb),jk,iblk(5,jc,jb)) + &
            ptr_coeff(6,1,jc,jb)*p_vn_in(iidx(6,jc,jb),jk,iblk(6,jc,jb)) + &
            ptr_coeff(7,1,jc,jb)*p_vn_in(iidx(7,jc,jb),jk,iblk(7,jc,jb)) + &
            ptr_coeff(8,1,jc,jb)*p_vn_in(iidx(8,jc,jb),jk,iblk(8,jc,jb)) + &
            ptr_coeff(9,1,jc,jb)*p_vn_in(iidx(9,jc,jb),jk,iblk(9,jc,jb))

          grad_y(jc,jk,jb) =                                               &
            ptr_coeff(1,2,jc,jb)*p_vn_in(iidx(1,jc,jb),jk,iblk(1,jc,jb)) + &
            ptr_coeff(2,2,jc,jb)*p_vn_in(iidx(2,jc,jb),jk,iblk(2,jc,jb)) + &
            ptr_coeff(3,2,jc,jb)*p_vn_in(iidx(3,jc,jb),jk,iblk(3,jc,jb)) + &
            ptr_coeff(4,2,jc,jb)*p_vn_in(iidx(4,jc,jb),jk,iblk(4,jc,jb)) + &
            ptr_coeff(5,2,jc,jb)*p_vn_in(iidx(5,jc,jb),jk,iblk(5,jc,jb)) + &
            ptr_coeff(6,2,jc,jb)*p_vn_in(iidx(6,jc,jb),jk,iblk(6,jc,jb)) + &
            ptr_coeff(7,2,jc,jb)*p_vn_in(iidx(7,jc,jb),jk,iblk(7,jc,jb)) + &
            ptr_coeff(8,2,jc,jb)*p_vn_in(iidx(8,jc,jb),jk,iblk(8,jc,jb)) + &
            ptr_coeff(9,2,jc,jb)*p_vn_in(iidx(9,jc,jb),jk,iblk(9,jc,jb))

        ENDDO
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    !$ACC END DATA

  END SUBROUTINE t_intp_vec_interpolate_r


  !--------------------------------------------------------------------------------
  !> REAL fields: Driver routine for the interpolation of
  !  cell-based variables at lon-lat grid points.
  !
  SUBROUTINE t_lon_lat_intp_interpolate_r( this, name, p_cell_in, nproma, p_lonlat_out, hintp_type)

    ! Indices of source points and interpolation coefficients
    CLASS(t_lon_lat_intp), INTENT(IN)           :: this
    CHARACTER(LEN=*),      INTENT(IN)           :: name                !< variable name
    ! input cell-based variable for which gradient at cell center is computed
    REAL(wp),              INTENT(IN)           :: p_cell_in(:,:,:)    !< dim: (nproma,nlev,nblks_c)
    INTEGER,               INTENT(IN)           :: nproma
    ! output lon-lat-based variable
    REAL(wp),              INTENT(INOUT)        :: p_lonlat_out(:,:,:) !< dim: (nproma,nlev,nblks_lonlat)
    ! horizontal interpolation type
    INTEGER,               INTENT(IN)           :: hintp_type
    
    ! Local Parameters:
    CHARACTER(*), PARAMETER :: routine = modname//"::t_lon_lat_intp_interpolate_r"
    INTEGER  :: nblks_lonlat, npromz_lonlat

    nblks_lonlat  = this%nblks_lonlat(nproma)
    npromz_lonlat = this%npromz_lonlat(nproma)
    
    SELECT CASE(hintp_type)
      
    CASE (HINTP_TYPE_LONLAT_RBF)
      ! RBF interpolation
      CALL this%rbf_c2l%interpolate( p_cell_in(:,:,:), nproma, nblks_lonlat, npromz_lonlat, &
        &                            p_lonlat_out(:,:,:))
      
    CASE (HINTP_TYPE_LONLAT_NNB)
      ! Nearest-neighbor interpolation
      CALL this%nnb%interpolate( p_cell_in(:,:,:), nproma, nblks_lonlat, npromz_lonlat, &
        &                        p_lonlat_out(:,:,:))
      
    CASE (HINTP_TYPE_LONLAT_BCTR)
      ! Barycentric interpolation
      CALL this%baryctr%interpolate( p_cell_in(:,:,:), nproma, nblks_lonlat, npromz_lonlat, &
        &                            p_lonlat_out(:,:,:))
      
    CASE DEFAULT
      CALL finish(routine, "Internal error with variable "//TRIM(name))
      
    END SELECT
    
  END SUBROUTINE t_lon_lat_intp_interpolate_r


  !--------------------------------------------------------------------------------
  !> INTEGER fields: Driver routine for the interpolation of
  !  cell-based variables at lon-lat grid points.
  !
  SUBROUTINE t_lon_lat_intp_interpolate_i( this, name, p_cell_in, nproma, p_lonlat_out, hintp_type)

    ! Indices of source points and interpolation coefficients
    CLASS(t_lon_lat_intp), TARGET, INTENT(IN)           :: this
    CHARACTER(LEN=*),              INTENT(IN)           :: name            !< variable name
    ! input cell-based variable for which gradient at cell center is computed
    INTEGER,                       INTENT(IN)           :: p_cell_in(:,:,:)    ! dim: (nproma,nlev,nblks_c)
    INTEGER,                       INTENT(IN)           :: nproma
    ! output lon-lat-based variable
    INTEGER,                       INTENT(INOUT)        :: p_lonlat_out(:,:,:) ! dim: (nproma,nlev,nblks_lonlat)
    ! horizontal interpolation type
    INTEGER,                       INTENT(IN)           :: hintp_type
    
    ! Local Parameters:
    CHARACTER(*), PARAMETER :: routine = modname//"::interpol_lonlat_int"
    INTEGER  :: nblks_lonlat, npromz_lonlat

    nblks_lonlat  = this%nblks_lonlat(nproma)
    npromz_lonlat = this%npromz_lonlat(nproma)
    
    SELECT CASE(hintp_type)
      
    CASE (HINTP_TYPE_LONLAT_NNB)
      ! Nearest-neighbor interpolation:
      CALL this%nnb%interpolate( p_cell_in(:,:,:), nproma, nblks_lonlat, npromz_lonlat, &
        &                        p_lonlat_out(:,:,:))
      
    CASE DEFAULT
      CALL finish(routine, "Internal error with variable "//TRIM(name))
      
    END SELECT
    
  END SUBROUTINE t_lon_lat_intp_interpolate_i


  !--------------------------------------------------------------------------------
  !> Performs vector RBF reconstruction at lon-lat grid points.
  !
  ! This routine is based on mo_intp_rbf::rbf_vec_interpol_cell()
  !
  SUBROUTINE t_lon_lat_intp_interpolate_rvec( this, p_vn_in, nproma, &
    &                                         grad_x, grad_y, hintp_type, opt_slev, opt_elev)

    ! Indices of source points and interpolation coefficients
    CLASS(t_lon_lat_intp), TARGET,   INTENT(IN)           :: this
    ! input normal components of (velocity) vectors at edge midpoints
    REAL(wp),                        INTENT(IN)           :: p_vn_in(:,:,:)  !< dim: (nproma,nlev,nblks_e)
    INTEGER,                         INTENT(IN)           :: nproma
    ! reconstructed x/y-components of velocity vector, dim: (nproma,nlev,nblks_lonlat)
    REAL(wp),                        INTENT(INOUT)        :: grad_x(:,:,:), grad_y(:,:,:)
    ! horizontal interpolation type
    INTEGER,                         INTENT(IN)           :: hintp_type
    ! optional vertical start/end level:
    INTEGER,                         INTENT(in), OPTIONAL :: opt_slev, opt_elev
    ! Local Parameters:
    CHARACTER(*), PARAMETER :: routine = modname//"::t_lon_lat_intp_interpolate_rvec"
    INTEGER  :: nblks_lonlat, npromz_lonlat

    nblks_lonlat  = this%nblks_lonlat(nproma)
    npromz_lonlat = this%npromz_lonlat(nproma)

    SELECT CASE(hintp_type)
      
    CASE (HINTP_TYPE_LONLAT_RBF)
      ! Nearest-neighbor interpolation:
      CALL this%rbf_vec%interpolate(p_vn_in, nproma, nblks_lonlat, npromz_lonlat, &
        &                           grad_x, grad_y, opt_slev, opt_elev)
     
    CASE DEFAULT
      CALL finish(routine, "Internal error.")
      
    END SELECT
    
  END SUBROUTINE t_lon_lat_intp_interpolate_rvec


END MODULE mo_intp_lonlat_types
