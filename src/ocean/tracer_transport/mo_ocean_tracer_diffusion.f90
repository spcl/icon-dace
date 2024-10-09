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

! Contains the implementation of velocity and tracer diffusion for the ICON ocean model.

!----------------------------
#include "omp_definitions.inc"
#include "icon_definitions.inc"
#include "iconfor_dsl_definitions.inc"
!----------------------------
MODULE mo_ocean_tracer_diffusion
  
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_ocean_nml,           ONLY: n_zlev
  USE mo_run_config,          ONLY: dtime
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_ocean_tracer_transport_types,  ONLY: t_ocean_tracer
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sync,                ONLY: sync_c, sync_e, sync_v, sync_patch_array, sync_patch_array_mult
  USE mo_exception,           ONLY: finish, warning !, message_text, message
  USE mtime,                  ONLY: MAX_DATETIME_STR_LEN
  USE mo_ocean_time_events,   ONLY: getCurrentDate_to_String
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_timer,               ONLY: timer_start, timer_stop, timers_level, timer_dif_vert
  USE mo_fortran_tools,       ONLY: set_acc_host_or_device
 
  IMPLICIT NONE
  
  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_ocean_tracer_diffusion'

  ! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'
  CHARACTER(LEN=12)           :: str_module    = 'oceDiffusion'  ! Output of module for 1 line debug
  INTEGER :: idt_src       = 1               ! Level of detail for 1 line debug
  
  LOGICAL :: eliminate_upper_diag = .false.
  !
  ! PUBLIC INTERFACE
  !
  INTEGER, PARAMETER :: top=1
  PUBLIC :: tracer_diffusion_horz
  PUBLIC :: tracer_diffusion_vertical_implicit
  
  CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: oldDateString = ""

CONTAINS
  
  !-------------------------------------------------------------------------
  !Subroutine computes the horizontal diffusive flux of an arbitrary tracer.
   SUBROUTINE tracer_diffusion_horz(patch_3D, trac_in, diff_flx, k_t, subset_range, lacc)
    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
    REAL(wp), INTENT(in)              :: trac_in(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp), INTENT(inout)           :: diff_flx(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp), OPTIONAL                :: k_t(:,:,:) !mixing coefficient for tracer    
    TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    !
    !Local variables
    INTEGER :: level, blockNo, edge_index
    INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: start_edge_index, end_edge_index, start_block, end_block
    LOGICAL :: is_in_domain
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    LOGICAL :: lzacc

    ! Pointers needed for GPU/OpenACC
    INTEGER,  DIMENSION(:,:),   POINTER :: dolic_e
    INTEGER,  DIMENSION(:,:,:), POINTER :: cell_idx
    INTEGER,  DIMENSION(:,:,:), POINTER :: cell_blk
    REAL(wp), DIMENSION(:,:,:), POINTER :: prism_thick_e
    REAL(wp), DIMENSION(:,:),   POINTER :: inv_dual_edge_length
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_ocediffusion:tracer_diffusion_horz')
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2d(1)
    edges_in_domain => patch_2D%edges%in_domain
    !-------------------------------------------------------------------------------
    dolic_e => patch_3D%p_patch_1d(1)%dolic_e
    cell_idx => patch_2D%edges%cell_idx
    cell_blk => patch_2D%edges%cell_blk
    prism_thick_e => patch_3D%p_patch_1d(1)%prism_thick_e
    inv_dual_edge_length => patch_2D%edges%inv_dual_edge_length
    !-------------------------------------------------------------------------------
    start_block = edges_in_domain%start_block
    end_block = edges_in_domain%end_block

    IF (PRESENT(subset_range)) THEN
      is_in_domain = subset_range%is_in_domain
    END IF

    CALL set_acc_host_or_device(lzacc, lacc)

  IF(PRESENT(k_t))THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, edge_index, level, &
!ICON_OMP il_c1, ib_c1, il_c2, ib_c2) ICON_OMP_DEFAULT_SCHEDULE
    !$ACC DATA COPYIN(k_t) IF(lzacc)
    DO blockNo = start_block, end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      diff_flx(:,:,blockNo) = 0.0_wp
      !$ACC END KERNELS

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO edge_index = start_edge_index, end_edge_index
        !Get indices of two adjacent triangles
        il_c1 = cell_idx(edge_index,blockNo,1)
        ib_c1 = cell_blk(edge_index,blockNo,1)
        il_c2 = cell_idx(edge_index,blockNo,2)
        ib_c2 = cell_blk(edge_index,blockNo,2)

        DO level=1,  dolic_e(edge_index,blockNo)
          
          diff_flx(edge_index,level,blockNo) = &
            &   k_t(edge_index,level,blockNo) &
            & * prism_thick_e(edge_index,level,blockNo)  &
            & * (trac_in(il_c2,level,ib_c2) - trac_in(il_c1,level,ib_c1))       &
            & * inv_dual_edge_length(edge_index,blockNo)
            
        ENDDO
        
      ENDDO
      !$ACC END PARALLEL
    ENDDO
    !$ACC WAIT(1)
    !$ACC END DATA
!ICON_OMP_END_PARALLEL_DO

    IF (PRESENT(subset_range)) THEN
      IF (.NOT. is_in_domain) &
        & CALL sync_patch_array(sync_e, patch_2D, diff_flx)
    ENDIF
  ELSEIF(.NOT.PRESENT(k_t))THEN  
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, edge_index, level, &
!ICON_OMP il_c1, ib_c1, il_c2, ib_c2) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = start_block, end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      diff_flx(:,:,blockNo) = 0.0_wp
      !$ACC END KERNELS

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO edge_index = start_edge_index, end_edge_index
        !Get indices of two adjacent triangles
        il_c1 = cell_idx(edge_index,blockNo,1)
        ib_c1 = cell_blk(edge_index,blockNo,1)
        il_c2 = cell_idx(edge_index,blockNo,2)
        ib_c2 = cell_blk(edge_index,blockNo,2)

        DO level=1, dolic_e(edge_index,blockNo)
          
          diff_flx(edge_index,level,blockNo) = &
            & prism_thick_e(edge_index,level,blockNo)  &
            & * (trac_in(il_c2,level,ib_c2) - trac_in(il_c1,level,ib_c1))       &
            & * inv_dual_edge_length(edge_index,blockNo)
            
        ENDDO
        
      ENDDO
      !$ACC END PARALLEL
    ENDDO
    !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO

    IF (PRESENT(subset_range)) THEN
      IF (.NOT. is_in_domain) &
        & CALL sync_patch_array(sync_e, patch_2D, diff_flx)
    ENDIF

  ENDIF  
    
  END SUBROUTINE tracer_diffusion_horz
  !-------------------------------------------------------------------------

  
  
  !-------------------------------------------------------------------------
  !Subroutine computes the horizontal diffusive flux of an arbitrary tracer.
  !  Note: this was an attempt to create a vectorized version, which proved to be slower
  !  To be further examined.
  !  We disable it for now
!    SUBROUTINE tracer_diffusion_horz(patch_3D, trac_in, diff_flx, k_t)
!     TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
!     REAL(wp), INTENT(in)              :: trac_in(nproma,n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
!     REAL(wp), INTENT(inout)           :: diff_flx(nproma,n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
!     REAL(wp), OPTIONAL                :: k_t(:,:,:) !mixing coefficient for tracer    
! !     TYPE(t_subset_range), TARGET, INTENT(in), OPTIONAL :: subset_range
!     !
!     !Local variables
!     INTEGER :: level, blockNo, edge_index
!     INTEGER :: il_c1(nproma), ib_c1(nproma), il_c2(nproma), ib_c2(nproma)
!     INTEGER :: start_edge_index, end_edge_index
!     TYPE(t_subset_range), POINTER :: edges_in_domain
!     TYPE(t_patch), POINTER :: patch_2D
!     integer :: imx
!     ! CHARACTER(len=max_char_length), PARAMETER :: &
!     !        & routine = ('mo_ocediffusion:tracer_diffusion_horz')
!     !-------------------------------------------------------------------------------
!     patch_2D        => patch_3D%p_patch_2d(1)
!     edges_in_domain => patch_2D%edges%in_domain
!     !-------------------------------------------------------------------------------
!     
!   IF(PRESENT(k_t))THEN
! !ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, edge_index, level, &
! !ICON_OMP il_c1, ib_c1, il_c2, ib_c2) ICON_OMP_DEFAULT_SCHEDULE
!     DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!       CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
!       diff_flx(:,:,blockNo) = 0.0_wp
!       
!       DO edge_index = start_edge_index, end_edge_index          
!         !Get indices of two adjacent triangles
!         il_c1(edge_index) = patch_2D%edges%cell_idx(edge_index,blockNo,1)
!         ib_c1(edge_index) = patch_2D%edges%cell_blk(edge_index,blockNo,1)
!         il_c2(edge_index) = patch_2D%edges%cell_idx(edge_index,blockNo,2)
!         ib_c2(edge_index) = patch_2D%edges%cell_blk(edge_index,blockNo,2)
!       ENDDO
! 
! #ifndef __LVECTOR__
!       imx =  MAXVAL(patch_3d%p_patch_1d(1)%dolic_e(start_edge_index:end_edge_index,blockNo))     
! !     DO level=1,  MAXVAL(patch_3d%p_patch_1d(1)%dolic_e(start_edge_index:end_edge_index,blockNo))     
! !NEC$ outerloop_unroll(4)
!       DO level=1,  imx
!         DO edge_index = start_edge_index, end_edge_index
!           IF (patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo) < level) CYCLE          
! #else     
!       DO edge_index = start_edge_index, end_edge_index
!         DO level=1,  patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)
! #endif
! 
!           diff_flx(edge_index,level,blockNo) = &
!             &   k_t(edge_index,level,blockNo) &
!             & * patch_3D%p_patch_1d(1)%prism_thick_e(edge_index,level,blockNo)  &
!             & * (trac_in(il_c2(edge_index),level,ib_c2(edge_index))     -       &
!             &   trac_in(il_c1(edge_index),level,ib_c1(edge_index)))       &
!             & * patch_2D%edges%inv_dual_edge_length(edge_index,blockNo)
!             
!         ENDDO
!         
!       ENDDO
!     ENDDO
! !ICON_OMP_END_PARALLEL_DO
! 
! !     IF (PRESENT(subset_range)) THEN
! !       IF (.NOT. subset_range%is_in_domain) &
! !         & CALL sync_patch_array(sync_e, patch_2D, diff_flx)
! !     ENDIF    
!     
!   ELSEIF(.NOT.PRESENT(k_t))THEN  
! !ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, edge_index, level, &
! !ICON_OMP il_c1, ib_c1, il_c2, ib_c2) ICON_OMP_DEFAULT_SCHEDULE
!     DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
!       CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
!       diff_flx(:,:,blockNo) = 0.0_wp
!       DO edge_index = start_edge_index, end_edge_index
!         !Get indices of two adjacent triangles
!         il_c1(edge_index) = patch_2D%edges%cell_idx(edge_index,blockNo,1)
!         ib_c1(edge_index) = patch_2D%edges%cell_blk(edge_index,blockNo,1)
!         il_c2(edge_index) = patch_2D%edges%cell_idx(edge_index,blockNo,2)
!         ib_c2(edge_index) = patch_2D%edges%cell_blk(edge_index,blockNo,2)
!       ENDDO
!       
! #ifndef __LVECTOR__
! !     DO level=1,  MAXVAL(patch_3d%p_patch_1d(1)%dolic_e(start_edge_index:end_edge_index,blockNo))     
!       imx = MAXVAL(patch_3d%p_patch_1d(1)%dolic_e(start_edge_index:end_edge_index,blockNo))
! !NEC$ outerloop_unroll(4)
!       DO level = 1, imx
!         DO edge_index = start_edge_index, end_edge_index
!           IF (patch_3d%p_patch_1d(1)%dolic_e(edge_index,blockNo) < level) CYCLE          
! #else     
!       DO edge_index = start_edge_index, end_edge_index
!         DO level=1,  patch_3D%p_patch_1d(1)%dolic_e(edge_index,blockNo)
! #endif          
!           diff_flx(edge_index,level,blockNo) = &
!             & patch_3D%p_patch_1d(1)%prism_thick_e(edge_index,level,blockNo)  &
!             & * (trac_in(il_c2(edge_index),level,ib_c2(edge_index)) -         &
!             &    trac_in(il_c1(edge_index),level,ib_c1(edge_index)))       &
!             & * patch_2D%edges%inv_dual_edge_length(edge_index,blockNo)
!             
!         ENDDO
!         
!       ENDDO
!     ENDDO
! !ICON_OMP_END_PARALLEL_DO
! 
! !     IF (PRESENT(subset_range)) THEN
! !       IF (.NOT. subset_range%is_in_domain) &
! !         & CALL sync_patch_array(sync_e, patch_2D, diff_flx)
! !     ENDIF 
! 
!   ENDIF  
!     
!   END SUBROUTINE tracer_diffusion_horz
  !-------------------------------------------------------------------------


  !------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE tracer_diffusion_vertical_implicit( &
    & patch_3d,                  &
    & ocean_tracer,              &
    & a_v,                       &
    & h, lacc) !,  &
    ! & diff_column)

    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_ocean_tracer), TARGET :: ocean_tracer
    REAL(wp), INTENT(inout)              :: a_v(:,:,:)
    REAL(wp), INTENT(in)                 :: h(:,:)
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    !
    INTEGER :: cell_block, start_index, end_index
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    LOGICAL :: lzacc

    start_timer(timer_dif_vert,4)

    !-----------------------------------------------------------------------
    cells_in_domain       =>  patch_3d%p_patch_2d(1)%cells%in_domain
    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    IF (oldDateString /= getCurrentDate_to_String()) THEN 
      eliminate_upper_diag = .not. eliminate_upper_diag
      oldDateString = getCurrentDate_to_String()
!       IF (my_process_is_stdio()) THEN
!         write (0,*) "eliminate_upper_diag changed to ", eliminate_upper_diag, " at ", TRIM(oldDateString)
!       ENDIF
    ENDIF
    !-----------------------------------------------------------------------

#ifdef __LVECTOR__
! vectorized for NEC
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index) ICON_OMP_DEFAULT_SCHEDULE
      DO cell_block = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, cell_block, start_index, end_index)

        CALL tracer_diffusion_vertical_implicit_onBlock_lvector( &
          & patch_3d,                       &
          & ocean_tracer,                   &
          & a_v(:,:,cell_block),            &
          & h(:,cell_block),                &
          & cell_block, start_index, end_index, lacc=lzacc)

      END DO
!ICON_OMP_END_PARALLEL_DO
#else
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index) ICON_OMP_DEFAULT_SCHEDULE
      DO cell_block = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, cell_block, start_index, end_index)

        CALL tracer_diffusion_vertical_implicit_onBlock_cpuvector( &
          & patch_3d,                       &
          & ocean_tracer,                   &
          & a_v(:,:,cell_block),            &
          & h(:,cell_block),                &
          & cell_block, start_index, end_index, lacc=lzacc)

      END DO
!ICON_OMP_END_PARALLEL_DO
#endif


    stop_timer(timer_dif_vert,4)
    
  END SUBROUTINE tracer_diffusion_vertical_implicit
  !------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !!Subroutine implements implicit vertical diffusion for scalar fields.
  !>
  !! sbr identical to sbr above but now with homogeneous boundary conditions
  !!
  !!
  !! The result ocean_tracer%concetration is calculated on domain_cells
  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE tracer_diffusion_vertical_implicit_onBlock( &
    & patch_3d,                &
    & ocean_tracer,            &
    & a_v,                     &
    & h,                       &
    & blockNo, start_index, end_index) !,  &
    ! & diff_column)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_ocean_tracer), TARGET :: ocean_tracer
    REAL(wp), INTENT(inout)              :: a_v(:,:)
    REAL(wp), INTENT(in)                 :: h(:)
    INTEGER, INTENT(in) :: blockNo, start_index, end_index
    !
    !
    REAL(wp) :: inv_prism_thickness(1:n_zlev), inv_prisms_center_distance(1:n_zlev)
    REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)! , nb(1:n_zlev)
    REAL(wp) :: fact(1:n_zlev)
    REAL(wp) :: column_tracer(1:n_zlev)
    REAL(wp), POINTER :: field_column(:,:,:)
    INTEGER :: bottom_level
    INTEGER :: cell_index, level
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    REAL(wp) :: top_cell_thickness

!     REAL(wp) :: tmp, tmp_add
    !-----------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2d%cells%in_domain
    field_column    => ocean_tracer%concentration
    !-----------------------------------------------------------------------
!     dt_inv = 1.0_wp/dtime
    
    DO cell_index = start_index, end_index
      bottom_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)

      IF (bottom_level < 2 ) CYCLE ! nothing to diffuse

      top_cell_thickness = &
              & patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(cell_index,1,blockNo) + h(cell_index)

      inv_prism_thickness(1) = 1.0_wp /  top_cell_thickness
      
      inv_prisms_center_distance(2) = 1.0_wp / ( 0.5_wp * &
              & (top_cell_thickness + patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(cell_index,2,blockNo)))
      
      DO level=2,bottom_level
        inv_prism_thickness(level)        = patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_index,level,blockNo)
       ENDDO
      DO level=3,bottom_level
        inv_prisms_center_distance(level) = patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(cell_index,level,blockNo)
      ENDDO
      
      DO level=1,bottom_level
       column_tracer(level) = field_column(cell_index,level,blockNo)
      ENDDO

      !------------------------------------
      ! Fill triangular matrix
      ! b is diagonal, a is the upper diagonal, c is the lower
      !   top level
      a(1) = 0.0_wp
      c(1) = -a_v(cell_index,2) * inv_prism_thickness(1) * inv_prisms_center_distance(2)*dtime
      b(1) = 1.0_wp - c(1)
      DO level = 2, bottom_level-1
        a(level) = - a_v(cell_index,level)   * inv_prism_thickness(level) * inv_prisms_center_distance(level)*dtime
        c(level) = - a_v(cell_index,level+1) * inv_prism_thickness(level) * inv_prisms_center_distance(level+1)*dtime
        b(level) = 1.0_wp - a(level) - c(level)
      END DO
      ! bottom
      a(bottom_level) = -a_v(cell_index,bottom_level) * &
        & inv_prism_thickness(bottom_level) * inv_prisms_center_distance(bottom_level)*dtime
      b(bottom_level) = 1.0_wp - a(bottom_level)
      c(bottom_level) = 0.0_wp

      IF (eliminate_upper_diag) THEN
        ! solve the tridiagonal matrix by eliminating c (the upper diagonal) 
        DO level=bottom_level-1,1,-1
          fact(level)=c(level)/b(level+1)
          b(level)=b(level)-a(level+1)*fact(level)
          c(level) = 0.0_wp
          column_tracer(level) = column_tracer(level) - fact(level)*column_tracer(level+1)
        ENDDO
        
        ocean_tracer%concentration(cell_index,1,blockNo) = column_tracer(1)/b(1)
        DO level=2,bottom_level
         field_column(cell_index,level,blockNo) = (column_tracer(level) - &
            a(level)* field_column(cell_index,level-1,blockNo)) / b(level)    
        ENDDO
        
      ELSE
        ! solve the tridiagonal matrix by eliminating a (the lower diagonal) 
        DO level=2, bottom_level
          fact(level)=a(level)/b(level-1)
          b(level)=b(level)-c(level-1)*fact(level)
          a(level) = 0.0_wp
          column_tracer(level) = column_tracer(level) - fact(level)*column_tracer(level-1)
        ENDDO
        ocean_tracer%concentration(cell_index,bottom_level,blockNo) = column_tracer(bottom_level)/b(bottom_level)
        DO level=bottom_level-1,1,-1
         field_column(cell_index,level,blockNo) = (column_tracer(level) - &
            c(level)* field_column(cell_index,level+1,blockNo)) / b(level)    
        ENDDO                 
      
      ENDIF
    
    ENDDO ! cell_index
    
  END SUBROUTINE tracer_diffusion_vertical_implicit_onBlock
  !------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE tracer_diffusion_vertical_implicit_onBlock_lvector( &
    & patch_3d,                &
    & ocean_tracer,            &
    & a_v,                     &
    & h,                       &
    & blockNo, start_index, end_index, lacc) !,  &
    ! & diff_column)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_ocean_tracer), TARGET :: ocean_tracer
    REAL(wp), INTENT(inout)              :: a_v(:,:)
    REAL(wp), INTENT(in)                 :: h(:)
    INTEGER, INTENT(in) :: blockNo, start_index, end_index
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    !
    !
    REAL(wp) :: inv_prism_thickness(nproma,1:n_zlev), inv_prisms_center_distance(nproma,1:n_zlev)
    REAL(wp) :: a(nproma,1:n_zlev), b(nproma,1:n_zlev), c(nproma,1:n_zlev)! , nb(1:n_zlev)
    REAL(wp) :: fact(nproma,1:n_zlev)
    REAL(wp) :: column_tracer(nproma,1:n_zlev)
    REAL(wp), POINTER :: field_column(:,:,:)
    INTEGER :: bottom_level(nproma), max_bottom_level
    INTEGER :: cell_index, level
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    REAL(wp) :: top_cell_thickness(nproma)
    LOGICAL :: lzacc
    CHARACTER(len=*), PARAMETER :: routine = modname//':tracer_diffusion_vertical_implicit_onBlock_lvector'

    !-----------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2d%cells%in_domain
    field_column    => ocean_tracer%concentration
    !-----------------------------------------------------------------------
!     dt_inv = 1.0_wp/dtime

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA COPYIN(a_v, h, patch_3d%p_patch_1d(1)%dolic_c, patch_3d%p_patch_1d(1)%inv_prism_thick_c) &
    !$ACC   COPYIN(patch_3d%p_patch_1d(1)%inv_prism_center_dist_c, patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_c) &
    !$ACC   CREATE(a, b, bottom_level, c, column_tracer, fact, inv_prism_thickness) &
    !$ACC   CREATE(inv_prisms_center_distance, top_cell_thickness) &
    !$ACC   COPY(field_column) IF(lzacc)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO cell_index = start_index, end_index
      bottom_level(cell_index) = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
      IF (bottom_level(cell_index) < 2 ) CYCLE ! nothing to diffuse

      top_cell_thickness(cell_index) = &
              & patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(cell_index,1,blockNo) + h(cell_index)

      inv_prism_thickness(cell_index,1) = 1.0_wp /  top_cell_thickness(cell_index)
      
      inv_prisms_center_distance(cell_index,2) = 1.0_wp / ( 0.5_wp * &
          & (top_cell_thickness(cell_index) + patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c(cell_index,2,blockNo))) 
    ENDDO
    !$ACC END PARALLEL
    !$ACC WAIT(1)

    max_bottom_level = -1
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) REDUCTION(MAX: max_bottom_level) IF(lzacc)
    DO level = start_index, end_index
      max_bottom_level = MAX(max_bottom_level, bottom_level(level))
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO level=1,max_bottom_level
      DO cell_index = start_index, end_index
        IF (bottom_level(cell_index) < 2 .OR. level > bottom_level(cell_index)) CYCLE ! nothing to diffuse
          
        column_tracer(cell_index,level) = field_column(cell_index,level,blockNo)
        IF (level >= 2) &  ! it should work without this
            inv_prism_thickness(cell_index, level)   = patch_3d%p_patch_1d(1)%inv_prism_thick_c(cell_index,level,blockNo)
        IF (level >= 3) &  ! it should work without this
            inv_prisms_center_distance(cell_index, level) = patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(cell_index,level,blockNo)          
         
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !------------------------------------
    ! Fill triangular matrix
    ! b is diagonal, a is the upper diagonal, c is the lower
    
    !  top level
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO cell_index = start_index, end_index
      IF (bottom_level(cell_index) < 2) CYCLE ! nothing to diffuse
      a(cell_index,1) = 0.0_wp
      c(cell_index,1) = -a_v(cell_index,2) * inv_prism_thickness(cell_index,1) * inv_prisms_center_distance(cell_index,2)*dtime
      b(cell_index,1) = 1.0_wp - c(cell_index,1)
    ENDDO

    !$ACC LOOP SEQ
    DO level=2,max_bottom_level-1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO cell_index = start_index, end_index
      IF (level >= bottom_level(cell_index)) CYCLE ! nothing to diffuse        
        a(cell_index,level) = - a_v(cell_index,level)   * inv_prism_thickness(cell_index,level) * inv_prisms_center_distance(cell_index,level)*dtime
        c(cell_index,level) = - a_v(cell_index,level+1) * inv_prism_thickness(cell_index,level) * inv_prisms_center_distance(cell_index,level+1)*dtime
        b(cell_index,level) = 1.0_wp - a(cell_index,level) - c(cell_index,level)
      END DO
    ENDDO
    
    ! bottom
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO cell_index = start_index, end_index
      IF (bottom_level(cell_index) < 2) CYCLE ! nothing to diffuse        
      a(cell_index,bottom_level(cell_index)) = -a_v(cell_index,bottom_level(cell_index)) * &
        & inv_prism_thickness(cell_index,bottom_level(cell_index)) * inv_prisms_center_distance(cell_index,bottom_level(cell_index))*dtime
      b(cell_index,bottom_level(cell_index)) = 1.0_wp - a(cell_index,bottom_level(cell_index))
      c(cell_index,bottom_level(cell_index)) = 0.0_wp
    ENDDO
            
    IF (eliminate_upper_diag) THEN
      ! solve the tridiagonal matrix by eliminating c (the upper diagonal) 
      !$ACC LOOP SEQ
      DO level=max_bottom_level-1,1,-1
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO cell_index = start_index, end_index
          IF (level >= bottom_level(cell_index) .OR. bottom_level(cell_index) < 2) CYCLE ! nothing to diffuse        
          
          fact(cell_index,level)=c(cell_index,level)/b(cell_index,level+1)
          b(cell_index,level)=b(cell_index,level)-a(cell_index,level+1)*fact(cell_index,level)
          c(cell_index,level) = 0.0_wp
          column_tracer(cell_index,level) = column_tracer(cell_index,level) - fact(cell_index,level)*column_tracer(cell_index,level+1)
          
        ENDDO
      ENDDO

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO cell_index = start_index, end_index
        IF (bottom_level(cell_index) < 2) CYCLE ! nothing to diffuse
        field_column(cell_index,1,blockNo) = column_tracer(cell_index,1)/b(cell_index,1)
      ENDDO

      !$ACC LOOP SEQ
      DO level=2,max_bottom_level
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO cell_index=start_index,end_index
          IF (level > bottom_level(cell_index) .OR. bottom_level(cell_index) < 2) CYCLE ! nothing to diffuse        
          
          field_column(cell_index,level,blockNo) = (column_tracer(cell_index,level) - &
              a(cell_index,level)* field_column(cell_index,level-1,blockNo)) / b(cell_index,level)    
        ENDDO
      ENDDO
      
    ELSE
    
      ! solve the tridiagonal matrix by eliminating a (the lower diagonal)        
      !$ACC LOOP SEQ
      DO level=2, max_bottom_level
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO cell_index=start_index,end_index
          IF (level > bottom_level(cell_index) .OR. bottom_level(cell_index) < 2) CYCLE ! nothing to diffuse        
        
          fact(cell_index,level)=a(cell_index,level)/b(cell_index,level-1)
          b(cell_index,level)=b(cell_index,level)-c(cell_index,level-1)*fact(cell_index,level)
          a(cell_index,level) = 0.0_wp
          column_tracer(cell_index,level) = column_tracer(cell_index,level) - fact(cell_index,level)*column_tracer(cell_index,level-1)
        ENDDO
      ENDDO

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO cell_index=start_index,end_index
        IF (bottom_level(cell_index) < 2) CYCLE ! nothing to diffuse        
        field_column(cell_index,bottom_level(cell_index),blockNo) = &
            column_tracer(cell_index,bottom_level(cell_index))/b(cell_index,bottom_level(cell_index))
      ENDDO

      !$ACC LOOP SEQ
      DO level=max_bottom_level-1,1,-1
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO cell_index=start_index,end_index
          IF (level >= bottom_level(cell_index) .OR. bottom_level(cell_index) < 2) CYCLE ! nothing to diffuse        
        
          field_column(cell_index,level,blockNo) = (column_tracer(cell_index,level) - &
            c(cell_index,level) * field_column(cell_index,level+1,blockNo)) / b(cell_index,level)    
        
        ENDDO
      ENDDO
      
    ENDIF  ! eliminate_upper_diag
    !$ACC END PARALLEL
    !$ACC WAIT(1)
    !$ACC END DATA
  END SUBROUTINE tracer_diffusion_vertical_implicit_onBlock_lvector
  !------------------------------------------------------------------------
 
  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE tracer_diffusion_vertical_implicit_onBlock_cpuvector( &
    & patch_3d,                &
    & ocean_tracer,            &
    & a_v,                     &
    & h,                       &
    & blockNo, start_index, end_index, lacc) !,  &
    ! & diff_column)
    
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(t_ocean_tracer), TARGET :: ocean_tracer
    REAL(wp), INTENT(inout)              :: a_v(:,:)
    REAL(wp), INTENT(in)                 :: h(:)
    INTEGER, INTENT(in) :: blockNo, start_index, end_index
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    !
    !
    REAL(wp) :: inv_prism_thickness(nproma,1:n_zlev), inv_prisms_center_distance(nproma,1:n_zlev)
    REAL(wp) :: a(nproma,1:n_zlev), b(nproma,1:n_zlev), c(nproma,1:n_zlev)! , nb(1:n_zlev)
    REAL(wp) :: fact(nproma,1:n_zlev)
    REAL(wp) :: column_tracer(nproma,1:n_zlev)
    REAL(wp), POINTER :: field_column(:,:,:)
    INTEGER :: bottom_level(nproma), max_bottom_level
    INTEGER :: cell_index, level
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2d
    REAL(wp) :: top_cell_thickness(nproma)
    LOGICAL :: lzacc

!     REAL(wp) :: tmp, tmp_add

    ! Pointers needed for GPU/OpenACC
    INTEGER, POINTER :: dolic_c(:,:)
    REAL(wp), POINTER :: inv_prism_thick_c(:,:,:), inv_prism_center_dist_c(:,:,:)
    REAL(wp), POINTER :: prism_thick_flat_sfc_c(:,:,:)

    !-----------------------------------------------------------------------
    patch_2d        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2d%cells%in_domain
    field_column    => ocean_tracer%concentration
    !-----------------------------------------------------------------------
    dolic_c => patch_3d%p_patch_1d(1)%dolic_c
    inv_prism_thick_c => patch_3d%p_patch_1d(1)%inv_prism_thick_c
    inv_prism_center_dist_c => patch_3d%p_patch_1d(1)%inv_prism_center_dist_c
    prism_thick_flat_sfc_c => patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c
    !-----------------------------------------------------------------------
!     dt_inv = 1.0_wp/dtime

    CALL set_acc_host_or_device(lzacc, lacc)

#ifdef NAGFOR
    inv_prism_thickness(:,:) = 0.0_wp
    inv_prisms_center_distance(:,:) = 0.0_wp
#endif

    !$ACC DATA CREATE(a, b, bottom_level, c, column_tracer, fact, inv_prism_thickness) &
    !$ACC   CREATE(inv_prisms_center_distance, top_cell_thickness) IF(lzacc)

    max_bottom_level = -1
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) REDUCTION(MAX: max_bottom_level) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO cell_index = start_index, end_index
      bottom_level(cell_index) = dolic_c(cell_index,blockNo)
      max_bottom_level = MAX(max_bottom_level, bottom_level(cell_index))
      IF (bottom_level(cell_index) < 2 ) CYCLE ! nothing to diffuse

      top_cell_thickness(cell_index) = &
              & prism_thick_flat_sfc_c(cell_index,1,blockNo) + h(cell_index)

      inv_prism_thickness(cell_index,1) = 1.0_wp /  top_cell_thickness(cell_index)
      
      inv_prisms_center_distance(cell_index,2) = 1.0_wp / ( 0.5_wp * &
          & (top_cell_thickness(cell_index) + prism_thick_flat_sfc_c(cell_index,2,blockNo))) 
          
      DO level=2,bottom_level(cell_index)
        inv_prism_thickness(cell_index, level)        = inv_prism_thick_c(cell_index,level,blockNo)
       ENDDO
      DO level=3,bottom_level(cell_index)
        inv_prisms_center_distance(cell_index,level) = inv_prism_center_dist_c(cell_index,level,blockNo)
      ENDDO
      
      DO level=1,bottom_level(cell_index)
       column_tracer(cell_index, level) = field_column(cell_index,level,blockNo)
      ENDDO
          
    ENDDO
    !$ACC END PARALLEL
    !$ACC WAIT(1)

    !------------------------------------
    ! Fill triangular matrix
    ! b is diagonal, a is the upper diagonal, c is the lower

    !  top level
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO cell_index = start_index, end_index
!       IF (bottom_level(cell_index) < 2) CYCLE ! nothing to diffuse
      a(cell_index,1) = 0.0_wp
      c(cell_index,1) = -a_v(cell_index,2) * inv_prism_thickness(cell_index,1) * inv_prisms_center_distance(cell_index,2)*dtime
    ENDDO

    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO cell_index = start_index, end_index
      b(cell_index,1) = 1.0_wp - c(cell_index,1)
    ENDDO

    !$ACC LOOP SEQ
    DO level=2, max_bottom_level-1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO cell_index = start_index, end_index
        a(cell_index,level) = - a_v(cell_index,level)   * inv_prism_thickness(cell_index,level) * inv_prisms_center_distance(cell_index,level)*dtime
        c(cell_index,level) = - a_v(cell_index,level+1) * inv_prism_thickness(cell_index,level) * inv_prisms_center_distance(cell_index,level+1)*dtime
      END DO
    ENDDO

    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO cell_index = start_index, end_index
      DO level=2,bottom_level(cell_index)-1
        b(cell_index,level) = 1.0_wp - a(cell_index,level) - c(cell_index,level)
      END DO
    ENDDO

    ! bottom
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO cell_index = start_index, end_index
      IF (bottom_level(cell_index) < 2) CYCLE ! nothing to diffuse        
      a(cell_index,bottom_level(cell_index)) = -a_v(cell_index,bottom_level(cell_index)) * &
        & inv_prism_thickness(cell_index,bottom_level(cell_index)) * inv_prisms_center_distance(cell_index,bottom_level(cell_index))*dtime
      c(cell_index,bottom_level(cell_index)) = 0.0_wp
      b(cell_index,bottom_level(cell_index)) = 1.0_wp - a(cell_index,bottom_level(cell_index))
    ENDDO

    IF (eliminate_upper_diag) THEN
      ! solve the tridiagonal matrix by eliminating c (the upper diagonal)
      !$ACC LOOP SEQ
      DO level=max_bottom_level-1,1,-1
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO cell_index = start_index, end_index
          IF (level >= bottom_level(cell_index) .OR. bottom_level(cell_index) < 2) CYCLE ! nothing to diffuse        
          
          fact(cell_index,level)=c(cell_index,level)/b(cell_index,level+1)
          b(cell_index,level)=b(cell_index,level)-a(cell_index,level+1)*fact(cell_index,level)
          c(cell_index,level) = 0.0_wp
          column_tracer(cell_index,level) = column_tracer(cell_index,level) - fact(cell_index,level)*column_tracer(cell_index,level+1)
          
        ENDDO
      ENDDO

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO cell_index = start_index, end_index
        IF (bottom_level(cell_index) < 2) CYCLE ! nothing to diffuse
        field_column(cell_index,1,blockNo) = column_tracer(cell_index,1)/b(cell_index,1)
      ENDDO

      !$ACC LOOP SEQ
      DO level=2,max_bottom_level
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO cell_index=start_index,end_index
          IF (level > bottom_level(cell_index) .OR. bottom_level(cell_index) < 2) CYCLE ! nothing to diffuse        
          
          field_column(cell_index,level,blockNo) = (column_tracer(cell_index,level) - &
              a(cell_index,level)* field_column(cell_index,level-1,blockNo)) / b(cell_index,level)    
        ENDDO
      ENDDO
      
    ELSE

      ! solve the tridiagonal matrix by eliminating a (the lower diagonal)
      !$ACC LOOP SEQ
      DO level=2, max_bottom_level
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO cell_index=start_index,end_index
          IF (level > bottom_level(cell_index) .OR. bottom_level(cell_index) < 2) CYCLE ! nothing to diffuse        
        
          fact(cell_index,level)=a(cell_index,level)/b(cell_index,level-1)
          b(cell_index,level)=b(cell_index,level)-c(cell_index,level-1)*fact(cell_index,level)
          a(cell_index,level) = 0.0_wp
          column_tracer(cell_index,level) = column_tracer(cell_index,level) - fact(cell_index,level)*column_tracer(cell_index,level-1)
        ENDDO
      ENDDO

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO cell_index=start_index,end_index
        IF (bottom_level(cell_index) < 2) CYCLE ! nothing to diffuse        
        field_column(cell_index,bottom_level(cell_index),blockNo) = &
            column_tracer(cell_index,bottom_level(cell_index))/b(cell_index,bottom_level(cell_index))
      ENDDO

      !$ACC LOOP SEQ
      DO level=max_bottom_level-1,1,-1
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO cell_index=start_index,end_index
          IF (level >= bottom_level(cell_index) .OR. bottom_level(cell_index) < 2) CYCLE ! nothing to diffuse        
        
          field_column(cell_index,level,blockNo) = (column_tracer(cell_index,level) - &
            c(cell_index,level) * field_column(cell_index,level+1,blockNo)) / b(cell_index,level)    
        
        ENDDO
      ENDDO
      
    ENDIF  ! eliminate_upper_diag
    !$ACC END PARALLEL
    !$ACC WAIT(1)
    !$ACC END DATA
  END SUBROUTINE tracer_diffusion_vertical_implicit_onBlock_cpuvector
  !------------------------------------------------------------------------

END MODULE mo_ocean_tracer_diffusion
