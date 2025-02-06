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

! Meta-data type definitions for ICON variables.

MODULE mo_var_metadata_types
  USE ISO_C_BINDING,            ONLY: C_SIZE_T, C_LOC, c_ptr
  USE mo_kind,                  ONLY: dp, wp, sp
  USE mo_impl_constants,        ONLY: TLEV_NNOW, vname_len, &
    & VINTP_METHOD_LIN, HINTP_TYPE_LONLAT_RBF
  USE mo_grib2,                 ONLY: t_grib2_var
  USE mo_action_types,          ONLY: t_var_action
  USE mo_cf_convention,         ONLY: t_cf_var
  USE mo_tracer_metadata_types, ONLY: t_tracer_meta
  USE mo_model_domain,          ONLY: t_subset_range
  USE mo_var_groups,            ONLY: MAX_GROUPS
  USE mo_cdi,                   ONLY: TSTEP_INSTANT, CDI_UNDEFID
  USE mo_exception,             ONLY: finish
  USE mo_util_stride,           ONLY: util_get_ptrdiff
  USE mo_util_libc,             ONLY: memcpy_f

  IMPLICIT NONE
  PRIVATE

  CHARACTER(*), PARAMETER :: modname = 'mo_var_metadata_types'

  ! ---------------------------------------------------------------
  ! CONSTANTS

  ! list of vertical interpolation types
  ! 
  ! A variable can have any combination of this which means that it
  ! can be interpolated vertically in these different ways.
  CHARACTER(len=1), PARAMETER :: VINTP_TYPE_LIST(3) = &
    (/ "Z",  "P", "I" /)


  ! list of available post-op's (small arithmetic operations on
  ! fields). The implementation is placed in "mo_post_op.f90".o
  INTEGER, PARAMETER, PUBLIC   :: POST_OP_NONE      = -1  !< trivial post-op ("do nothing")
  INTEGER, PARAMETER, PUBLIC   :: POST_OP_SCALE     =  1  !< multiply by scalar factor "arg1"
  INTEGER, PARAMETER, PUBLIC   :: POST_OP_RHO       =  2  !< multiply by rho to get densities instead
  INTEGER, PARAMETER, PUBLIC   :: POST_OP_LUC       =  3  !< convert landuse classes from internal values 
                                                          !< to GRIB2 values (table 4.243) and vice versa. 
  INTEGER, PARAMETER, PUBLIC   :: POST_OP_LIN2DBZ   =  4  !< convert linear values to dbz: dbzval = 10*log10(val)
  INTEGER, PARAMETER, PUBLIC   :: POST_OP_OFFSET    =  5  !< add offset value arg1

  ! list of available variable classes
  INTEGER, PARAMETER, PUBLIC :: CLASS_DEFAULT       = 0
  INTEGER, PARAMETER, PUBLIC :: CLASS_TILE          = 1   !< variable contains tile-specific information
  INTEGER, PARAMETER, PUBLIC :: CLASS_TILE_LAND     = 2   !< variable contains tile-specific information
                                                          !< but is restricted to land-tiles only
  INTEGER, PARAMETER, PUBLIC :: CLASS_SYNSAT        = 3
  INTEGER, PARAMETER, PUBLIC :: CLASS_CHEM          = 4   !< atmospheric chemical constituent (PDT 40)
  INTEGER, PARAMETER, PUBLIC :: CLASS_CHEM_STAT     = 5   !< atmospheric chemical constituent (PDT 42)
                                                          !< statistical process
  INTEGER, PARAMETER, PUBLIC :: CLASS_CHEM_OPTP     = 6   !< atmospheric chemical constituent (PDT 48)
                                                          !< optical properties
  INTEGER, PARAMETER, PUBLIC :: CLASS_DISTR         = 7   !< variable based on a distribuition function (PDT 57)
  INTEGER, PARAMETER, PUBLIC :: CLASS_DISTR_STAT    = 8   !< variable based on a distribuition function (PDT 67)
                                                          !< statistical process

  ! ---------------------------------------------------------------
  ! META-DATA TYPE DEFINITIONS

  TYPE :: t_union_vals
    REAL(dp) :: rval = 0._dp
    REAL(sp) :: sval = 0._sp
    INTEGER  :: ival = 0
    LOGICAL  :: lval = .FALSE.
  END type t_union_vals

  !> data specific for pz-level interpolation.
  TYPE :: t_vert_interp_meta
    ! meta data containing the groups to which a variable belongs
    LOGICAL  :: vert_intp_type(SIZE(VINTP_TYPE_LIST)) = .FALSE.
    INTEGER  :: vert_intp_method                      = VINTP_METHOD_LIN
    LOGICAL  :: l_hires_intp                          = .FALSE., &
         &      l_restore_fricred                     = .FALSE., &
         &      l_loglin                              = .FALSE., &
         &      l_extrapol                            = .TRUE.,  &
         &      l_satlimit                            = .FALSE., &
         &      l_restore_pbldev                      = .FALSE., &
         &      l_pd_limit                            = .FALSE., &
         &      l_restore_sfcinv, l_hires_corr
    REAL(wp) :: lower_limit = 0._wp, extrapol_dist
  END TYPE t_vert_interp_meta

  !> data specific for horizontal interpolation.
  TYPE :: t_hor_interp_meta
    INTEGER :: hor_intp_type = HINTP_TYPE_LONLAT_RBF ! NONE/RBF/Nearest-Neighbor/...
    INTEGER :: fallback_type = HINTP_TYPE_LONLAT_RBF ! replaces "hor_intp_type" if this is not feasible
    INTEGER :: lonlat_id     = 0 ! lon-lat grid (ID in global list)
  END TYPE t_hor_interp_meta


  !> This type defines small arithmetic operations ("post-ops") as
  !  post-processing tasks.
  !
  !  These post-processing tasks are restricted to point-wise
  !  operations (no halo synchronization) of a single field, like
  !  value scaling.
  !
  !  @note The "post-ops" are performed at output time and DO NOT
  !        MODIFY THE FIELD ITSELF.
  !
  TYPE t_post_op_meta
    INTEGER                    :: ipost_op_type = POST_OP_NONE !< type of post-processing operation
    LOGICAL                    :: lnew_cf       = .FALSE.
    TYPE(t_cf_var)             :: new_cf        = t_cf_var('', '', '', -1) !< CF information of modified field
    LOGICAL                    :: lnew_grib2    = .FALSE.
    TYPE(t_grib2_var)          :: new_grib2 !< GRIB2 information of modified field
    TYPE(t_union_vals)         :: arg1 !< post-op argument (e.g. scaling factor)
  END TYPE t_post_op_meta


  TYPE :: t_var_metadata
    INTEGER(C_SIZE_T)          :: c_interop = 1234 ! just a C-interoperable dummy member
    CHARACTER(len=vname_len)   :: name        = ''             ! variable name
    INTEGER                    :: var_class   = CLASS_DEFAULT  ! variable type
    INTEGER                    :: data_type   = -1             ! variable data type: REAL_T, SINGLE_T, INT_T, BOOL_T
    TYPE(t_cf_var)             :: cf          = t_cf_var('', '', '', -1)  ! CF convention information 
    TYPE(t_grib2_var)          :: grib2  ! GRIB2 related information
    LOGICAL                    :: allocated   = .FALSE.        ! allocation status
    INTEGER                    :: ndims       = 0              ! number of dimensions used
    INTEGER                    :: used_dimensions(5) = 0       ! final dimensions of variable
    LOGICAL                    :: lrestart    = .FALSE.        ! write field to restart
    LOGICAL                    :: loutput     = .TRUE.         ! write field to output
    INTEGER                    :: isteptype   = TSTEP_INSTANT  ! Type of statistical processing
    TYPE(t_union_vals)         :: resetval                     ! reset value for accumulated fields
    LOGICAL                    :: lrestart_cont = .FALSE.      ! continue if not in restart file     
    LOGICAL                    :: lrestart_read = .FALSE.      ! field has been set from restart file
    TYPE(t_union_vals)         :: initval                      ! value if not in restart file
    LOGICAL                    :: lcontainer   = .FALSE.       ! true, if this is a container
    LOGICAL                    :: lcontained   = .FALSE.       ! true, if this is in a container
    INTEGER                    :: ncontained   = 0             ! index in container
    INTEGER                    :: maxcontained = 0             ! container size   
    INTEGER                    :: var_ref_pos  = -1            ! for containers: dimension index for references
    INTEGER                    :: hgrid        = -1            ! CDI horizontal grid type
    INTEGER                    :: vgrid        = -1            ! CDI vertical grid type
    TYPE(t_subset_range)       :: subset                       ! subset for latter field access
    INTEGER                    :: dom          = -1            ! domain (used to be pointer to varlist%patch_id 
    INTEGER                    :: tlev_source  = TLEV_NNOW     ! Information where to find the actual
    !                                                     timelevel for timelevel dependent variables:        
    !                                                      = 0 : nnow
    !                                                      = 1 : nnow_rcf
    !                                                      ... more may follow
    INTEGER                    :: cdiVarID     = CDI_UNDEFID
    INTEGER                    :: cdiGridID    = CDI_UNDEFID
    TYPE(t_post_op_meta)       :: post_op               !<  "post-op" (small arithmetic operations) for this variable
    TYPE(t_var_action)         :: action_list
    ! Metadata for vertical/horizontal interpolation
    !
    ! Note that setting these parameters to non-default values does
    ! not mean that interpolation is actually performed for this
    ! variables (this is controlled by namelist settings) but only
    ! that this is possible!
    TYPE(t_vert_interp_meta)   :: vert_interp 
    TYPE(t_hor_interp_meta)    :: hor_interp 
    ! meta data containing the groups to which a variable belongs
    LOGICAL :: in_group(MAX_GROUPS)

    ! Flag: defines, if this field is updated by the internal
    ! post-processing scheduler
    INTEGER :: l_pp_scheduler_task             = 0

    ! Metadata for missing value masking

    LOGICAL                    :: lmiss = .FALSE.           ! flag: true, if variable should be initialized with missval
    TYPE(t_union_vals)         :: missval                   ! missing value
    LOGICAL                    :: lmask_boundary = .FALSE.  ! flag: true, if interpolation zone should be masked *in output*
    ! Index of tracer in tracer and in diagnostics container
    INTEGER                    :: idx_tracer   = -1         !< index of tracer in tracer container
    INTEGER                    :: idx_diag     = -1         !< index of tracer in diagnostics container
    LOGICAL                    :: lopenacc     = .FALSE.    ! Variable exists on GPU
  END TYPE t_var_metadata

  TYPE t_var_metadata_ptr
    TYPE(t_var_metadata), POINTER :: p => NULL()
  END TYPE t_var_metadata_ptr
  ! The type t_var_metadata_dynamic is (in contrast to t_var_metadata) not transfered to the output PE.
  ! This allows for dynamical objects inside t_var_metadata_dynamic like pointers or allocatables.
  TYPE t_var_metadata_dynamic
    CLASS(t_tracer_meta), ALLOCATABLE :: tracer      ! Tracer-specific metadata
  END TYPE t_var_metadata_dynamic

  PUBLIC :: VINTP_TYPE_LIST
  PUBLIC :: t_union_vals
  PUBLIC :: t_var_metadata, t_var_metadata_ptr, var_metadata_get_size
  PUBLIC :: var_metadata_toBinary, var_metadata_fromBinary
  PUBLIC :: t_var_metadata_dynamic
  PUBLIC :: t_vert_interp_meta
  PUBLIC :: t_hor_interp_meta
  PUBLIC :: t_post_op_meta

CONTAINS

  INTEGER FUNCTION var_metadata_get_size() RESULT(isize)
    INTEGER, SAVE :: isize_saved = 0
    LOGICAL, SAVE :: init = .FALSE.

    IF (init) THEN
      isize = isize_saved
    ELSE
      isize = (size_byte() + 3) / 4
      isize_saved = isize
      init = .TRUE.
    END IF
  END FUNCTION var_metadata_get_size

  INTEGER(C_SIZE_T) FUNCTION size_byte() RESULT(bsize)
    TYPE(t_var_metadata), TARGET :: dummy(2)
    INTEGER(C_SIZE_T), SAVE :: bsize_saved = 0
    LOGICAL, SAVE :: init = .FALSE.

    IF (init) THEN
      bsize = bsize_saved
    ELSE
      bsize = util_get_ptrdiff(dummy(1)%c_interop, dummy(2)%c_interop)
      bsize_saved = bsize
      init = .TRUE.
    END IF
  END FUNCTION size_byte

  FUNCTION var_metadata_toBinary(info, isize) RESULT(bin)
    TYPE(t_var_metadata), INTENT(IN), TARGET :: info
    INTEGER, INTENT(IN) :: isize
    INTEGER, TARGET :: bin(isize)
    TYPE(c_ptr) :: dummy_cptr

    IF (var_metadata_get_size() .GT. isize) &
      & CALL finish("var_metadata_toBinary", "size mismatch")
    dummy_cptr = memcpy_f(C_LOC(bin), C_LOC(info), size_byte())
  END FUNCTION var_metadata_toBinary

  FUNCTION var_metadata_fromBinary(bin, isize) RESULT(info)
    INTEGER, INTENT(IN) :: isize
    INTEGER, INTENT(IN), TARGET :: bin(isize)
    TYPE(t_var_metadata), TARGET :: info
    TYPE(c_ptr) :: dummy_cptr

    IF (var_metadata_get_size() .GT. isize) &
      & CALL finish("var_metadata_fromBinary", "size mismatch")
    dummy_cptr = memcpy_f(C_LOC(info), C_LOC(bin), size_byte())
  END FUNCTION var_metadata_fromBinary

END MODULE mo_var_metadata_types
