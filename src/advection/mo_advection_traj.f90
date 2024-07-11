! Some utilities which are specific to the transport algorithm.
! Routines are mostly dealing with the computation of backward
! trajectories and/or departure regions.
!
! Module contains some functions and procedures which are specifically related
! to the transport schemes. These subroutines or functions are needed at
! various places within the transport scheme. Therefore outsourcing these
! routines protects from possible circular dependencies.
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_advection_traj

  USE mo_kind,                ONLY: wp, vp, i4
  USE mo_exception,           ONLY: finish
  USE mo_model_domain,        ONLY: t_patch
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_parallel_config,     ONLY: nproma
  USE mo_loopindices,         ONLY: get_indices_e
  USE mo_impl_constants,      ONLY: min_rledge_int, SUCCESS
  USE mo_timer,               ONLY: timer_start, timer_stop, timers_level, timer_back_traj
  USE mo_advection_utils,     ONLY: t_list2D

! Some compiler NVHPC versions have issues with the fast atomic path
#ifdef _OPENACC
#define _USE_FAST_ATOMIC
#endif

#if defined (_USE_FAST_ATOMIC) && (__NVCOMPILER_MAJOR__ == 21 && __NVCOMPILER_MINOR__ == 3)
#undef _USE_FAST_ATOMIC
#endif

#if defined (_USE_FAST_ATOMIC) && (__NVCOMPILER_MAJOR__ == 99 && __NVCOMPILER_PATCHLEVEL__ <= 224325)
#undef _USE_FAST_ATOMIC
#endif


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: btraj_dreg
  PUBLIC :: t_back_traj
  PUBLIC :: btraj_compute_o1

  TYPE t_back_traj
    ! line indices of cell centers in which the calculated barycenters are located
    ! dim: (nproma,nlev,p_patch%nblks_e)
    INTEGER, CONTIGUOUS, POINTER :: cell_idx(:,:,:) => NULL()
    !
    ! block indices of cell centers in which the calculated barycenters are located
    ! dim: (nproma,nlev,p_patch%nblks_e)
    INTEGER, CONTIGUOUS, POINTER :: cell_blk(:,:,:) => NULL()
    !
    ! distance vectors cell center --> barycenter of advected area (geographical coordinates)
    ! dim: (nproma,nlev,p_patch%nblks_e,2)
    REAL(vp), CONTIGUOUS, POINTER :: distv_bary(:,:,:,:) => NULL()

  CONTAINS
    !
    PROCEDURE :: construct
    PROCEDURE :: destruct
    
  END TYPE t_back_traj

  CHARACTER(len=*), PARAMETER :: modname = 'mo_advection_traj'

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Allocate object components
  !!
  !! Allocates all components of the object of class t_back_traj
  !!
  SUBROUTINE construct(obj, nproma, nlev, nblks, ncoord)
    CLASS(t_back_traj) :: obj
    INTEGER, INTENT(IN) :: nproma
    INTEGER, INTENT(IN) :: nlev
    INTEGER, INTENT(IN) :: nblks
    INTEGER, INTENT(IN) :: ncoord
    !
    ! local
    INTEGER :: ist

    CHARACTER(len=*), PARAMETER :: routine = modname//':construct'

    ALLOCATE(obj%cell_idx(nproma,nlev,nblks), &
      &      obj%cell_blk(nproma,nlev,nblks), &
      &      obj%distv_bary(nproma,nlev,nblks,ncoord), STAT=ist)
    IF (ist /= SUCCESS) CALL finish(routine, 'allocation failed')

    !$ACC ENTER DATA CREATE(obj)
    !$ACC ENTER DATA CREATE(obj%cell_idx, obj%cell_blk, obj%distv_bary)

  END SUBROUTINE construct



  !-------------------------------------------------------------------------
  !>
  !! Deallocate object components
  !!
  !! Deallocates all components of the object of class t_back_traj
  !!
  SUBROUTINE destruct(obj)
    CLASS(t_back_traj) :: obj
    !
    ! local
    INTEGER :: ist

    CHARACTER(len=*), PARAMETER :: routine = modname//':destruct'

    IF (ASSOCIATED(obj%cell_idx)) THEN

      !$ACC WAIT(1)
      !$ACC EXIT DATA DELETE(obj%cell_idx, obj%cell_blk, obj%distv_bary)
      !$ACC EXIT DATA DELETE(obj)

      DEALLOCATE(obj%cell_idx, obj%cell_blk, obj%distv_bary, STAT=ist)
      IF (ist /= SUCCESS) CALL finish(routine, 'deallocation failed')

    ENDIF

  END SUBROUTINE destruct



  !-------------------------------------------------------------------------
  !>
  !! Computation of first order backward trajectories for FFSL transport scheme
  !!
  !! Computes backward trajectories in order to determine the barycenter of the
  !! departure region. Here, a simple first order method is used. Computations are
  !! performed on a plane tangent to the edge midpoint. Coordinate axes point into
  !! the local normal and tangential direction.
  !! Once the barycenter of the departure region is known, the distance vector
  !! between the circumcenter of the upstream cell and the barycenter is computed.
  !! In a final step, this vector is transformed into a rotated coordinate system
  !! which has its origin at the circumcenter. The coordinate axes point to the local
  !! east and local north. Note that this subroutine has specifically been designed
  !! for the MIURA scheme with second order (linear) reconstruction of the subgrid 
  !! distribution.
  !!
  SUBROUTINE btraj_compute_o1( btraj, ptr_p, ptr_int, p_vn, p_vt, p_dthalf, &
    &                          opt_rlstart, opt_rlend, opt_slev, opt_elev,  &
    &                          opt_acc_async )

    TYPE(t_back_traj), INTENT(INOUT) :: btraj

    TYPE(t_patch), TARGET, INTENT(in) ::      &  !< patch on which computation is performed
         &  ptr_p

    TYPE(t_int_state), TARGET, INTENT(in) ::  &  !< pointer to data structure for interpolation
         &  ptr_int

    REAL(wp), INTENT(IN)    ::  &  !< normal component of velocity vector at edge midpoints
         &  p_vn(:,:,:)

    REAL(wp), INTENT(IN)    ::  &  !< tangential component of velocity vector at
         &  p_vt(:,:,:)               !< edge midpoints

    REAL(wp), INTENT(IN)    ::  &  !< $0.5 \Delta t$
         &  p_dthalf

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
         &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
         &  opt_rlend                     !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
         &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
         &  opt_elev

    LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async   !< optional async OpenACC

    REAL(wp) :: z_ntdistv_bary_1, z_ntdistv_bary_2      !< cell center --> barycenter in 'normal' and
                                                        !< 'tangential' coordinates.

    INTEGER :: je, jk, jb        !< index of edge, vert level, block
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rlstart, i_rlend
    INTEGER :: slev, elev        !< vertical start and end level
! These convenience pointers are needed to avoid PGI trying to copy derived type instance btraj back from device to host
    INTEGER, CONTIGUOUS, POINTER  :: p_cell_idx(:,:,:), p_cell_blk(:,:,:)
    REAL(vp), CONTIGUOUS, POINTER :: p_distv_bary(:,:,:,:)
    LOGICAL :: lvn_pos

    !-------------------------------------------------------------------------

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF

    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = ptr_p%nlev
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 4
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 2
    ENDIF


    i_startblk = ptr_p%edges%start_block(i_rlstart)
    i_endblk   = ptr_p%edges%end_block(i_rlend)

    !-------------------------------------------------------------------------
    IF (timers_level > 5) THEN
      CALL timer_start(timer_back_traj)
    ENDIF

! These convenience pointers are needed to avoid PGI trying to copy derived type instance btraj back from device to host
    p_cell_idx   => btraj%cell_idx
    p_cell_blk   => btraj%cell_blk
    p_distv_bary => btraj%distv_bary

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,z_ntdistv_bary_1,z_ntdistv_bary_2,lvn_pos) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_p, jb, i_startblk, i_endblk,        &
           i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR TILE(64, 2)
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
          !
          ! Calculate backward trajectories
          !

          ! position of barycenter in normal direction
          ! pos_barycenter_1 = - p_vn(je,jk,jb) * p_dthalf

          ! position of barycenter in tangential direction
          ! pos_barycenter_2 = - p_vt(je,jk,jb) * p_dthalf

          ! logical auxiliary for MERGE operations: .TRUE. for vn >= 0
          lvn_pos = p_vn(je,jk,jb) >= 0._wp

          ! If vn > 0 (vn < 0), the upwind cell is cell 1 (cell 2)

          ! line and block indices of neighbor cell with barycenter
          p_cell_idx(je,jk,jb) = &
             &   MERGE(ptr_p%edges%cell_idx(je,jb,1),ptr_p%edges%cell_idx(je,jb,2),lvn_pos)

          p_cell_blk(je,jk,jb) = &
             &   MERGE(ptr_p%edges%cell_blk(je,jb,1),ptr_p%edges%cell_blk(je,jb,2),lvn_pos)


          ! Calculate the distance cell center --> barycenter for the cell,
          ! in which the barycenter is located. The distance vector points
          ! from the cell center to the barycenter.
          z_ntdistv_bary_1 =  - ( p_vn(je,jk,jb) * p_dthalf     &
               & + MERGE(ptr_int%pos_on_tplane_e(je,1,1,jb),        &
               &         ptr_int%pos_on_tplane_e(je,2,1,jb),lvn_pos))

          z_ntdistv_bary_2 =  - ( p_vt(je,jk,jb) * p_dthalf     &
               & + MERGE(ptr_int%pos_on_tplane_e(je,1,2,jb),        &
               &         ptr_int%pos_on_tplane_e(je,2,2,jb),lvn_pos))

          ! In a last step, transform this distance vector into a rotated
          ! geographical coordinate system with its origin at the circumcenter
          ! of the upstream cell. Coordinate axes point to local East and local
          ! North.

          ! component in longitudinal direction
          p_distv_bary(je,jk,jb,1) =                                                       &
               &   z_ntdistv_bary_1*MERGE(ptr_p%edges%primal_normal_cell(je,jb,1)%v1,         &
               &                           ptr_p%edges%primal_normal_cell(je,jb,2)%v1,lvn_pos) &
               & + z_ntdistv_bary_2*MERGE(ptr_p%edges%dual_normal_cell(je,jb,1)%v1,           &
               &                           ptr_p%edges%dual_normal_cell(je,jb,2)%v1,lvn_pos)

          ! component in latitudinal direction
          p_distv_bary(je,jk,jb,2) =                                                       &
               &   z_ntdistv_bary_1*MERGE(ptr_p%edges%primal_normal_cell(je,jb,1)%v2,         &
               &                           ptr_p%edges%primal_normal_cell(je,jb,2)%v2,lvn_pos) &
               & + z_ntdistv_bary_2*MERGE(ptr_p%edges%dual_normal_cell(je,jb,1)%v2,           &
               &                           ptr_p%edges%dual_normal_cell(je,jb,2)%v2,lvn_pos)

        ENDDO ! loop over edges
      ENDDO   ! loop over vertical levels
      !$ACC END PARALLEL

    END DO    ! loop over blocks

!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (timers_level > 5) CALL timer_stop(timer_back_traj)

    IF ( PRESENT(opt_acc_async) ) THEN
      IF ( opt_acc_async ) THEN
        RETURN
      END IF
    END IF
    !$ACC WAIT(1)

  END SUBROUTINE btraj_compute_o1
 
 

  !-------------------------------------------------------------------------
  !>
  !! Computation of first order backward trajectories for FFSL transport scheme
  !!
  !! Computes backward trajectories in order to determine an approximation to the
  !! departure region. Here, the departure region is approximated by a rhomboid 
  !! with the help of first order accurate backward trajectories which start at 
  !! edge vertices. Computations are performed on a plane tangent to the edge 
  !! midpoint. Base vectors of this coordinate system (S1) point into the 
  !! local tangential and normal direction. Once the departure region vertices 
  !! are known w.r.t. S1, they are transformed into second coordinate frame (S2) 
  !! which follows from S1 by translation and rotation. The origin of S2 is 
  !! located at the cell circumcenter of the upstream cell, with the base 
  !! vectors pointing to local east and local north.
  !! So far, we take care that the departure region vertices are stored in
  !! counterclockwise order. This ensures that the following gaussian 
  !! quadrature is positive definite.
  !!
  !! This subroutine may be combined with any reconstruction method 
  !! for the subgrid distribution.
  !!
  !! NOTE_1: Since we are only interested in the departure region average rather than 
  !!       the departure region integral, counterclockwise numbering is not strictly 
  !!       necessary. Maybe we should remove the computational overhead of counterclockwise 
  !!       numbering at some time. However, the vertices must not be numbered in 
  !!       random order. Care must be taken that the points are numbered either 
  !!       clockwise or counterclockwise. 
  !!       
  !! Note_2: The coordinates for 2 of the 4 vertices are time independent. However, 
  !!       tests indicated that re-computing these coordinates is faster than fetching 
  !!       precomputed ones from memory. 
  !!
  SUBROUTINE btraj_dreg( ptr_p, ptr_int, p_vn, p_vt, p_dt, lcounterclock, &
       &                   p_cell_idx, p_cell_blk, p_coords_dreg_v,       &
       &                   opt_rlstart, opt_rlend, opt_slev, opt_elev,    &
       &                   opt_falist )

    TYPE(t_patch), TARGET, INTENT(IN) ::     &  !< patch on which computation is performed
         &  ptr_p

    TYPE(t_int_state), TARGET, INTENT(IN) :: &  !< pointer to data structure for interpolation
         &  ptr_int

    REAL(wp), INTENT(IN)    ::  &  !< normal component of velocity vector at edge midpoints
         &  p_vn(:,:,:)

    REAL(wp), INTENT(IN)    ::  &  !< tangential component of velocity vector at
         &  p_vt(:,:,:)               !< edge midpoints

    REAL(wp), INTENT(IN)    ::  &  !< time step $\Delta t$
         &  p_dt

    LOGICAL, INTENT(IN)     ::  &  !< if TRUE, flux area vertices are ordered 
         &  lcounterclock            !< counterclockwise. If FALSE, some are ordered
                                     !< counterclockwise, some clockwise

    REAL(vp), INTENT(OUT) ::    &  !< coordinates of departure region vertices. The origin
         &  p_coords_dreg_v(:,:,:,:,:)!< of the coordinate system is at the circumcenter of
                                      !< the upwind cell. Base vectors point to local East
                                      !< and North. (geographical coordinates)
                                      !< dim: (nproma,4,2,nlev,ptr_p%nblks_e)

    INTEGER, INTENT(OUT)  ::    &  !< line and block indices of upwind cell
         &  p_cell_idx(:,:,:)   !< dim: (nproma,nlev,ptr_p%nblks_e)
    INTEGER, INTENT(OUT)  ::    &  !< line and block indices of upwind cell
         &  p_cell_blk(:,:,:)   !< dim: (nproma,nlev,ptr_p%nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
         &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
         &  opt_rlend                     !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
         &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
         &  opt_elev

    TYPE(t_list2D), INTENT(INOUT), OPTIONAL :: & !< list with points for which a local
         &  opt_falist                      !< polynomial approximation is insufficient 
                                            !< and a piecewise approximation is needed, 
                                            !< instead

    REAL(wp) ::            &       !< coordinates of departure points 
         &  depart_pts(2,2)        !< in edge-based coordinate system

    REAL(wp) ::            &       !< coordinates of departure region vertices
         &  pos_dreg_vert_c(4,2)   !< as seen from translated coordinate system.
                                   !< origin at circumcenter of upwind cell

    REAL(wp) ::            &       !< position on tangential plane depending
         &  pos_on_tplane_e(2)     !< on the sign of vn

    REAL(wp) ::            &       !< primal and dual normals of cell lying
         &  pn_cell_1, pn_cell_2, dn_cell_1, dn_cell_2    !< in the direction of vn

    REAL(wp) ::            &       !< edge vertices
         &  edge_verts(nproma,2,2)

    INTEGER :: je, jk, jb          !< index of edge, vert level, block
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rlstart, i_rlend
    INTEGER :: slev, elev          !< vertical start and end level
    LOGICAL :: lvn_pos             !< vn > 0: TRUE/FALSE
    LOGICAL :: lvn_sys_pos(nproma,ptr_p%nlev)   !< vn*system_orient > 0

    ! for index list generation
    LOGICAL :: llist_gen           !< if TRUE, generate index list
    INTEGER :: ie, ie_capture      !< counter, and its captured value

    INTEGER, PARAMETER :: &
         &  gang_size = 128, pts_per_thread = 20
    INTEGER :: num_gangs, jg, jl
    INTEGER :: gang_ie(1)
    INTEGER :: gang_captured_ie(1)                    !< OpenACC captures
    INTEGER(i4) :: gang_eidx(pts_per_thread*gang_size)!< OpenACC gang-local lists
    INTEGER(i4) :: gang_elev(pts_per_thread*gang_size)

    REAL(wp):: traj_length         !< backward trajectory length [m]
    REAL(wp):: e2c_length          !< edge-upwind cell circumcenter length [m]
    !-------------------------------------------------------------------------

    IF (timers_level > 5) THEN
      CALL timer_start(timer_back_traj)
    ENDIF

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF

    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = ptr_p%nlev
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 4
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 2
    ENDIF

    IF (PRESENT(opt_falist)) THEN
      llist_gen = .TRUE.
    ELSE
      llist_gen = .FALSE.
    ENDIF


    i_startblk = ptr_p%edges%start_block(i_rlstart)
    i_endblk   = ptr_p%edges%end_block(i_rlend)

    !$ACC DATA PRESENT(ptr_p, ptr_int, p_vn, p_vt, p_coords_dreg_v, p_cell_idx, p_cell_blk) &
    !$ACC   NO_CREATE(opt_falist) CREATE(edge_verts, lvn_sys_pos)

    IF (llist_gen) THEN
      !$ACC KERNELS ASYNC(1)
      opt_falist%len(:) = 0
      !$ACC END KERNELS
    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,ie,i_startidx,i_endidx,traj_length,e2c_length, &
!$OMP depart_pts,pos_dreg_vert_c,pos_on_tplane_e,pn_cell_1,pn_cell_2,dn_cell_1,dn_cell_2,lvn_pos,&
!$OMP lvn_sys_pos,edge_verts,ie_capture) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_p, jb, i_startblk, i_endblk, &
           i_startidx, i_endidx, i_rlstart, i_rlend)


      ! get local copy of edge vertices
      !$ACC PARALLEL ASYNC(1)
      !$ACC LOOP GANG VECTOR
!NEC$ ivdep
      DO je = i_startidx, i_endidx
        edge_verts(je,1:2,1:2) = ptr_int%pos_on_tplane_e(je,3:4,1:2,jb)
      ENDDO
      !$ACC END PARALLEL

      ! logical switch for merge options regarding the counterclockwise numbering
      IF (lcounterclock) THEN

        !$ACC PARALLEL ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = slev, elev
          DO je = i_startidx, i_endidx
            lvn_sys_pos(je,jk) = p_vn(je,jk,jb)*ptr_p%edges%tangent_orientation(je,jb) >= 0._wp
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ELSE
        !$ACC PARALLEL ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = slev, elev
          DO je = i_startidx, i_endidx
            lvn_sys_pos(je,jk) = .FALSE.
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDIF

      ! generate list of points that require special treatment
      !
      IF (llist_gen) THEN
        ie = 0
        ! OpenACC/GPU - using a two-stage atomics below: each gang has one atomic counter in its shared memory.
        ! There is then a global counter that is only updated by a single thread from each gang

        ! Nvidia compilers only support 65536 gangs
        ! So if nproma*nlev > 65536 * gang_size * pts_per_thread the FAST_ATOMIC code won't work
#ifdef __NVCOMPILER_MAJOR__
        IF ( (elev-slev+1)*(i_endidx-i_startidx+1) > 65535 * gang_size * pts_per_thread ) THEN
          CALL finish(modname//":btraj_dreg", & 
            "Too many grid points for the fast algorithm, please undefine _USE_FAST_ATOMIC and recompile")
        END IF
#endif

        num_gangs = ( (elev-slev+1)*(i_endidx-i_startidx+1) + gang_size-1) / gang_size
        !$ACC PARALLEL ASYNC(1) &
        !$ACC   PRIVATE(gang_ie, gang_captured_ie, gang_elev, gang_eidx) &
#ifndef _USE_FAST_ATOMIC
        !$ACC   COPY(ie) &
#endif
        !$ACC   NUM_GANGS(num_gangs) VECTOR_LENGTH(gang_size)
        gang_ie(1) = 0

        !$ACC LOOP GANG VECTOR PRIVATE(lvn_pos, traj_length, e2c_length, ie_capture) COLLAPSE(2)
        DO jk = slev, elev
          DO je = i_startidx, i_endidx
            ! logical switch for MERGE operations: .TRUE. for p_vn >= 0
            lvn_pos     = p_vn(je,jk,jb) >= 0._wp

            ! compute length of backward trajectory
            traj_length = SQRT(p_vn(je,jk,jb)**2 + p_vt(je,jk,jb)**2) * p_dt

            ! distance from edge midpoint to upwind cell circumcenter [m]
            e2c_length  = MERGE(ptr_p%edges%edge_cell_length(je,jb,1),       &
              &                 ptr_p%edges%edge_cell_length(je,jb,2),lvn_pos)

            IF (traj_length > 1.25_wp*e2c_length) THEN   ! add point to index list
#ifndef _USE_FAST_ATOMIC
              ! Default code path 
              !$ACC ATOMIC CAPTURE
              ie = ie + 1
              ie_capture = ie
              !$ACC END ATOMIC
              opt_falist%eidx(ie_capture,jb) = je
              opt_falist%elev(ie_capture,jb) = jk
#else
              ! Optimized OpenACC path
              !$ACC ATOMIC CAPTURE
              gang_ie(1) = gang_ie(1) + 1
              ie_capture = gang_ie(1)
              !$ACC END ATOMIC
              gang_eidx(ie_capture) = je
              gang_elev(ie_capture) = jk
#endif
            ENDIF
          ENDDO ! loop over edges
        ENDDO   ! loop over vertical levels

#ifndef _USE_FAST_ATOMIC
        !$ACC END PARALLEL

        ! Default code path
        ! store list dimension
        !$ACC KERNELS ASYNC(1) COPYIN(ie)
        opt_falist%len(jb) = ie
        !$ACC END KERNELS
#else
        ! Optimized OpenACC path
        ! Copy gang-local lists into the global array
        !$ACC ATOMIC CAPTURE
        gang_captured_ie(1) = opt_falist%len(jb)
        opt_falist%len(jb)  = opt_falist%len(jb) + gang_ie(1)
        !$ACC END ATOMIC

        !$ACC LOOP GANG
        DO jg = 1, num_gangs
          !$ACC LOOP VECTOR
          DO jl = 1, gang_size*pts_per_thread
            IF (jl <= gang_ie(1)) THEN
              opt_falist%eidx(gang_captured_ie(1)+jl, jb) = gang_eidx(jl)
              opt_falist%elev(gang_captured_ie(1)+jl, jb) = gang_elev(jl)
            END IF
          END DO
        END DO
        !$ACC END PARALLEL
#endif
      ENDIF

      !$ACC PARALLEL ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(depart_pts, pos_dreg_vert_c, pos_on_tplane_e) &
      !$ACC   PRIVATE(lvn_pos, pn_cell_1, pn_cell_2, dn_cell_1, dn_cell_2)
      DO jk = slev, elev
!DIR$ IVDEP, PREFERVECTOR
!$NEC ivdep
        DO je = i_startidx, i_endidx



          ! departure region and correct counterclockwise numbering of vertices
          !--------------------------------------------------------------------
          !
          ! Quadrilaterals show the position of the departure region, depending 
          ! on the sign of vn.
          !
          !        -1                          +1            : system orientation
          !
          !  3\--------------\2       4\--------------\3
          !    \              \         \              \         [vn < 0]
          !     \      N       \         \      N       \
          !      \      \       \         \      \       \      <- edge normal
          !   4,2 \------\-------\1      1 \------\-------\2,4  <- triangle edge
          !        \              \         \              \
          !         \              \         \              \
          !          \              \         \              \   [vn > 0]
          !         3 \--------------\4      2 \--------------\3


          ! Determine the upwind cell
          ! Cell indices are chosen such that the direction from cell 1 to cell 2
          ! is the positive direction of the normal vector N.
          ! Vertex indices are chosen such that the direction from vertex 1 to vertex 2
          ! is the positive direction of the tangential vector T.
          !
          ! If (T,N,Z) form a right-hand system, the system orientation is 1.
          ! If (T,N,Z) form a left-hand system, the system orientation is -1.

          !
          ! Calculate backward trajectories, starting at the two edge vertices
          ! (arrival points). It is assumed that the velocity vector is constant 
          ! along the edge.
          !

          ! logical switch for MERGE operations: .TRUE. for p_vn >= 0
          lvn_pos     = p_vn(je,jk,jb) >= 0._wp


          ! get line and block indices of upwind cell
          p_cell_idx(je,jk,jb) = MERGE(ptr_p%edges%cell_idx(je,jb,1),       &
            &                          ptr_p%edges%cell_idx(je,jb,2),lvn_pos)
          p_cell_blk(je,jk,jb) = MERGE(ptr_p%edges%cell_blk(je,jb,1),       &
            &                          ptr_p%edges%cell_blk(je,jb,2),lvn_pos)


          ! departure points of the departure cell. Point 1 belongs to edge-vertex 1, 
          ! point 2 belongs to edge_vertex 2.
          !
          ! position of vertex 4 (vn > 0) / vertex 2(vn < 0) in normal direction
          depart_pts(1,1)      = edge_verts(je,1,1) - p_vn(je,jk,jb) * p_dt

          ! position of vertex 4 (vn > 0) / vertex 2(vn < 0) in tangential direction
          depart_pts(1,2)      = edge_verts(je,1,2) - p_vt(je,jk,jb) * p_dt

          ! position of vertex 3 in normal direction
          depart_pts(2,1)      = edge_verts(je,2,1) - p_vn(je,jk,jb) * p_dt

          ! position of vertex 3 in tangential direction
          depart_pts(2,2)      = edge_verts(je,2,2) - p_vt(je,jk,jb) * p_dt



          ! determine correct position on tangential plane
          pos_on_tplane_e(1) = MERGE(ptr_int%pos_on_tplane_e(je,1,1,jb), &
               &                     ptr_int%pos_on_tplane_e(je,2,1,jb),lvn_pos)

          pos_on_tplane_e(2) = MERGE(ptr_int%pos_on_tplane_e(je,1,2,jb), &
               &                     ptr_int%pos_on_tplane_e(je,2,2,jb),lvn_pos)

          ! Calculate position of departure region vertices in a translated
          ! coordinate system. The origin is located at the circumcenter
          ! of the upwind cell. The distance vectors point from the cell center
          ! to the vertices.
          !
          ! Take care of correct counterclockwise numbering below
          !
          pos_dreg_vert_c(1,1:2) = edge_verts(je,1,1:2) - pos_on_tplane_e(1:2)

          pos_dreg_vert_c(2,1)   = MERGE(depart_pts(1,1),edge_verts(je,2,1),lvn_sys_pos(je,jk)) &
               &                 - pos_on_tplane_e(1)
          pos_dreg_vert_c(2,2)   = MERGE(depart_pts(1,2),edge_verts(je,2,2),lvn_sys_pos(je,jk)) &
               &                 - pos_on_tplane_e(2)

          pos_dreg_vert_c(3,1:2) = depart_pts(2,1:2) - pos_on_tplane_e(1:2)

          pos_dreg_vert_c(4,1)   = MERGE(edge_verts(je,2,1),depart_pts(1,1),lvn_sys_pos(je,jk)) &
               &                 - pos_on_tplane_e(1)
          pos_dreg_vert_c(4,2)   = MERGE(edge_verts(je,2,2),depart_pts(1,2),lvn_sys_pos(je,jk)) &
               &                 - pos_on_tplane_e(2)

          ! In a last step, these distance vectors are transformed into a rotated
          ! geographical coordinate system, which still has its origin at the circumcenter
          ! of the upwind cell. Now the coordinate axes point to local East and local
          ! North.
          !
          ! Determine primal and dual normals of the cell lying in the direction of vn
          pn_cell_1 = MERGE(ptr_p%edges%primal_normal_cell(je,jb,1)%v1,       &
               &             ptr_p%edges%primal_normal_cell(je,jb,2)%v1,lvn_pos)

          pn_cell_2 = MERGE(ptr_p%edges%primal_normal_cell(je,jb,1)%v2,       &
               &             ptr_p%edges%primal_normal_cell(je,jb,2)%v2,lvn_pos)

          dn_cell_1 = MERGE(ptr_p%edges%dual_normal_cell(je,jb,1)%v1,       &
               &             ptr_p%edges%dual_normal_cell(je,jb,2)%v1,lvn_pos)

          dn_cell_2 = MERGE(ptr_p%edges%dual_normal_cell(je,jb,1)%v2,       &
               &             ptr_p%edges%dual_normal_cell(je,jb,2)%v2,lvn_pos)

          ! components in longitudinal direction
          p_coords_dreg_v(je,1:4,1,jk,jb) =                                         &
               & pos_dreg_vert_c(1:4,1) * pn_cell_1 + pos_dreg_vert_c(1:4,2) * dn_cell_1

          ! components in latitudinal direction
          p_coords_dreg_v(je,1:4,2,jk,jb) =                                         &
               & pos_dreg_vert_c(1:4,1) * pn_cell_2 + pos_dreg_vert_c(1:4,2) * dn_cell_2

        ENDDO ! loop over edges
      ENDDO   ! loop over vertical levels
      !$ACC END PARALLEL
    END DO    ! loop over blocks

    !$ACC UPDATE HOST(opt_falist%len) ASYNC(1) IF(llist_gen)
    !$ACC WAIT(1)
    !$ACC END DATA

!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (timers_level > 5) CALL timer_stop(timer_back_traj)

  END SUBROUTINE btraj_dreg

END MODULE mo_advection_traj

