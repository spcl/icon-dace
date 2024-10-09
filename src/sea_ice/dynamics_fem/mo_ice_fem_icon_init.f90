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

! Contains allocation, initialization and setup of the FEM mesh
! in relation to the ICON grid.

!----------------------------
#include "omp_definitions.inc"
!----------------------------
!-----------------------------------------------------------------------
!
!  ! 1) Various initialization routines for the FEM sea-ice model
!  ! 2) Sync routines used in the EVP solver
!
!-----------------------------------------------------------------------
MODULE mo_ice_fem_icon_init
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_exception,           ONLY: message
  USE mo_impl_constants,      ONLY: max_char_length
  USE mo_impl_constants,      ONLY: success

  USE mo_parallel_config,     ONLY: nproma
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_math_constants,      ONLY: rad2deg, deg2rad

  USE mo_sync,                ONLY: sync_v, sync_patch_array
  USE mo_fortran_tools,       ONLY: set_acc_host_or_device

  IMPLICIT NONE

  PUBLIC  :: init_fem_wgts
  PUBLIC  :: init_fem_wgts_extra
  PUBLIC  :: destruct_fem_wgts
  PUBLIC  :: ice_fem_grid_init
  PUBLIC  :: ice_fem_grid_post

  PUBLIC  :: exchange_nod2D

  INTERFACE exchange_nod2d
    MODULE PROCEDURE exchange_nod2d
    MODULE PROCEDURE exchange_nod2Dx2
    MODULE PROCEDURE exchange_nod2Dx3
  END INTERFACE exchange_nod2d

  PRIVATE :: basisfunctions_nod
  PRIVATE :: calc_f_rot

  TYPE(t_patch), POINTER, PRIVATE :: fem_patch

  REAL(wp), ALLOCATABLE, PUBLIC   :: c2v_wgt(:,:,:), v2c_wgt(:,:,:)
  REAL(wp), ALLOCATABLE, PUBLIC  :: rot_mat(:,:,:,:), rot_mat_3D(:,:)
  ! Longitude and latitude of the north pole of the rotated grid
  ! These values put the north pole on the Indonesian/Malasian island Kalimantan and the south pole
  ! in north-west Brazil, near the border to Venezuela and Colombia.
  REAL(wp), PARAMETER, PUBLIC :: pollon = 114._wp*deg2rad, pollat = 0._wp
  !  the North pole is centered on Greenland (40W, 75N) and the South pole on Antarctica
  !  was used to test consistency of the results with respect to the change in the rotated grid
  !  REAL(wp), PARAMETER, PRIVATE :: pollon = -40._wp*deg2rad, pollat = 75._wp*deg2rad

CONTAINS

  !-------------------------------------------------------------------------
  !
  !> Constructor of weights for FEM sea-ice model
  !! We calculate only c2v_wgt (PRIVATE), which is used by cells2verts_scalar
  !!
  SUBROUTINE init_fem_wgts(p_patch_3D)
    TYPE(t_patch_3D), TARGET, INTENT(IN)    :: p_patch_3D

    !Local variables
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ice_fem_interface:init_fem_wgts'

    ! Patch
    TYPE(t_patch), POINTER :: p_patch
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

    ! Indexing
    INTEGER :: ist
    INTEGER :: jb, jv, ji
    INTEGER :: i_startidx_v, i_endidx_v
    INTEGER :: cell_index, cell_block

    ! The denominator
    REAL(wp) :: rdeno, rsum

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    p_patch => p_patch_3D%p_patch_2D(1)
    iidx => p_patch%verts%cell_idx
    iblk => p_patch%verts%cell_blk

    ! Weights for the cells2verts_scalar call
    ALLOCATE(c2v_wgt(nproma,6,p_patch%nblks_v),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating c2v_wgt failed')
    ENDIF
    c2v_wgt=0.0_wp

    ! loop through all patch verts
    !
!ICON_OMP_PARALLEL_DO PRIVATE(jb,i_startidx_v,i_endidx_v,jv,ji,rsum,     &
!ICON_OMP        rdeno,cell_index,cell_block) ICON_OMP_DEFAULT_SCHEDULE
  DO jb =  p_patch%verts%in_domain%start_block, p_patch%verts%in_domain%end_block
    CALL get_index_range(p_patch%verts%in_domain, jb, i_startidx_v, i_endidx_v)

    DO jv =  i_startidx_v, i_endidx_v

        rsum = 0.0_wp

        DO ji = 1, 6
    !    in case of pentagons and boundary nodes:
    !    cell_index = 0 when use_duplicated_connectivity = .false. (always the case for the ocean)
    !    cell_index = dummy values (duplicates) when use_duplicated_connectivity = .true.
          cell_index = iidx(jv,jb,ji)
          cell_block = iblk(jv,jb,ji)
          IF (cell_index > 0)                                      &
            & rsum = rsum +                                        &
            &      p_patch_3D%wet_c  (cell_index,1,cell_block) *   &
            &      p_patch%cells%area(cell_index,  cell_block)

        ENDDO

        rdeno =  1._wp/MAX( TINY(rdeno), rsum )

        DO ji = 1, 6
          cell_index = iidx(jv,jb,ji)
          cell_block = iblk(jv,jb,ji)
          IF (cell_index > 0)                                 &
            & c2v_wgt(jv,ji,jb) =                             &
            &   p_patch_3D%wet_c  (cell_index,1,cell_block) * &
            &   p_patch%cells%area(cell_index,  cell_block) * rdeno

        ENDDO

      ENDDO
 
    ENDDO
!ICON_OMP_END_PARALLEL_DO

    CALL message (TRIM(routine), 'end')        

  END SUBROUTINE init_fem_wgts

  !-------------------------------------------------------------------------
  !
  !> Construct interpolation weights for FEM sea-ice model:
  !! v2c_wgt: coeffs for interpolation from verts to cells by area weighting
  !! Identical to verts_aw_cells (see mo_intp_data_strc for description).
  !! Code copied from mo_intp_coeffs_lsq_bln.
  !!
  SUBROUTINE init_fem_wgts_extra(ptr_patch)

    TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

    !Local variables
    INTEGER :: i_startidx_c, i_endidx_c, nblks_c, ist
    INTEGER :: jc, je, jb, jv  ! integer over edges, blocks and levels
    INTEGER :: ile, ibe, ilv1,ilv2,ibv1,ibv2,ilv,ibv, idx_ce

    !-------------------------------------------------------------------------
    nblks_c = ptr_patch%nblks_c

    ! For the triangular ICON grid: geometry_info%cell_type = 3 (each cell has 3 vertices)
    ALLOCATE (v2c_wgt(nproma, ptr_patch%geometry_info%cell_type, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
        CALL finish ('mo_ice_fem_interface:init_fem_wgts_extra',                     &
        &             'allocation for v2c_wgt failed')
    ENDIF

    v2c_wgt = 0.0_wp

!ICON_OMP_PARALLEL_DO PRIVATE(jb,jc,je,jv,i_startidx_c,i_endidx_c,ile,ibe,idx_ce,&
!ICON_OMP                     ilv1,ilv2,ibv1,ibv2,ilv,ibv) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = ptr_patch%cells%all%start_block, ptr_patch%cells%all%end_block

    CALL get_index_range(ptr_patch%cells%all, jb, i_startidx_c, i_endidx_c)
    DO jc = i_startidx_c,i_endidx_c

       DO je = 1, ptr_patch%cells%num_edges(jc,jb)

          ile = ptr_patch%cells%edge_idx(jc,jb,je)
          ibe = ptr_patch%cells%edge_blk(jc,jb,je)
          IF ( ptr_patch%edges%cell_idx(ile,ibe,1) == jc .AND. &
               ptr_patch%edges%cell_blk(ile,ibe,1) == jb ) THEN
               idx_ce = 1
          ELSE
               idx_ce = 2
          ENDIF

          ilv1 = ptr_patch%edges%vertex_idx(ile,ibe,1)
          ibv1 = ptr_patch%edges%vertex_blk(ile,ibe,1)
          ilv2 = ptr_patch%edges%vertex_idx(ile,ibe,2)
          ibv2 = ptr_patch%edges%vertex_blk(ile,ibe,2)

          DO jv = 1, ptr_patch%cells%num_edges(jc,jb)
            ilv = ptr_patch%cells%vertex_idx(jc,jb,jv)
            ibv = ptr_patch%cells%vertex_blk(jc,jb,jv)

            IF (ilv == ilv1 .AND. ibv == ibv1) THEN
              v2c_wgt(jc,jv,jb) =   &
                v2c_wgt(jc,jv,jb) + &
                0.5_wp/ptr_patch%cells%area(jc,jb) *             &
                ptr_patch%edges%edge_cell_length(ile,ibe,idx_ce)*&
                ptr_patch%edges%edge_vert_length(ile,ibe,1)
            ENDIF
            IF (ilv == ilv2 .AND. ibv == ibv2) THEN
              v2c_wgt(jc,jv,jb)  =  &
                v2c_wgt(jc,jv,jb) + &
                0.5_wp/ptr_patch%cells%area(jc,jb) *             &
                ptr_patch%edges%edge_cell_length(ile,ibe,idx_ce)*&
                ptr_patch%edges%edge_vert_length(ile,ibe,2)
            ENDIF

          ENDDO

       ENDDO

    ENDDO !loop over all cells

  ENDDO   !loop over all blocks
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE init_fem_wgts_extra

  !-------------------------------------------------------------------------
  !
  !> Destructor of patch for FEM sea-ice model, deallocates c2v_wgt
  !!
  SUBROUTINE destruct_fem_wgts

    !Local variables
    INTEGER :: ist
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ice_fem_interface:destruct_fem_wgts'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    DEALLOCATE(c2v_wgt,STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'deallocating c2v_wgt failed')
    ENDIF

    IF (allocated(v2c_wgt)) DEALLOCATE(v2c_wgt,STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'deallocating v2c_wgt failed')
    ENDIF

    DEALLOCATE(rot_mat,STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'deallocating rot_mat failed')
    ENDIF

    DEALLOCATE(rot_mat_3D,STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'deallocating rot_mat_3D failed')
    ENDIF

    CALL message (TRIM(routine), 'end')        

  END SUBROUTINE destruct_fem_wgts
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !> Initialise the FEM grid
  !! This replaces the routine Sergey used to read the grid from file. We give values to
  !! coord_nod2D, nod2D, index_nod2D, elem2D_nodes and elem2D from the ICON grid.
  !!
  SUBROUTINE ice_fem_grid_init(p_patch_3D)

  USE mo_math_utilities,      ONLY: rotate_latlon, rotate_latlon_vec!, disp_new_vect

    ! For the AWI FEM
    USE mo_ice_fem_mesh,           ONLY: coord_nod2D, nod2D, index_nod2D, elem2D_nodes, elem2D, nod2D_elems
    USE mo_physical_constants, ONLY: earth_angular_velocity
    USE mo_ice_fem_mesh,           ONLY: coriolis_nod2D

    TYPE(t_patch_3D), TARGET, INTENT(in) :: p_patch_3D

    ! local variables
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ice_fem_interface:ice_fem_grid_init'

    ! Patch
    TYPE(t_patch), POINTER :: p_patch

    ! Indexing
    INTEGER :: k, jb, jc, jv
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: i_startidx_v, i_endidx_v
    INTEGER :: cell_index, cell_block

    ! Temporary variables/buffers and flags
    INTEGER :: verts(nproma,p_patch_3D%p_patch_2D(1)%nblks_v)
    INTEGER :: elems(nproma,p_patch_3D%p_patch_2D(1)%nblks_e)
    INTEGER :: buffy(nproma*p_patch_3D%p_patch_2D(1)%nblks_v)
    INTEGER :: ist
    REAL(wp) :: lat, lon, cos_d, sin_d
    REAL(wp) :: cos_lonp, sin_lonp, cos_latp, sin_latp

    ! Masking for the halo-points
    INTEGER :: halo_mask(nproma, p_patch_3D%p_patch_2D(1)%nblks_v)

!--------------------------------------------------------------------------------------------------

    p_patch => p_patch_3D%p_patch_2D(1)
    fem_patch => p_patch

    ! nod2D is the number of nodes = number of vertices
    nod2D = p_patch%n_patch_verts
    ALLOCATE(coord_nod2D(2,nod2D), index_nod2D(nod2D),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating coord_nod2D failed')
    ENDIF
    ALLOCATE(coriolis_nod2D(nod2D),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating coriolis_nod2D failed')
    ENDIF
    ALLOCATE(rot_mat(nproma,p_patch%nblks_v,2,2),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating rot_mat failed')
    ENDIF
    ALLOCATE(rot_mat_3D(3,3),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating rot_mat_3D failed')
    ENDIF


    ! Go through all the vertices and assign coordinates and land mask to coord_nod2D and
    ! index_nod2D
    k=0
    DO jb = p_patch%verts%all%start_block, p_patch%verts%all%end_block
      CALL get_index_range(p_patch%verts%all, jb, i_startidx_v, i_endidx_v) 
      DO jv = i_startidx_v,i_endidx_v
        k=k+1
        verts(jv,jb) = k
        lat = p_patch%verts%vertex(jv,jb)%lat
        lon = p_patch%verts%vertex(jv,jb)%lon

        ! Calculate the rotation matrix once
        ! see p.27 at http://www.cosmo-model.org/content/model/documentation/core/cosmoDyncsNumcs.pdf
        CALL rotate_latlon_vec( lon, lat, pollon, pollat, sin_d, cos_d )
        rot_mat(jv,jb,1,:) = (/ cos_d, sin_d /)
        rot_mat(jv,jb,2,:) = (/-sin_d, cos_d /)

        ! Use a rotated grid with the north pole at (pollon,pollat)
        CALL rotate_latlon(lat, lon, pollat, pollon)

!==========================================================================
        coriolis_nod2D(k) = calc_f_rot(lat,lon,pollat,earth_angular_velocity)

        ! x-coords in degrees
        coord_nod2D(1,k) = lon*rad2deg
        ! y-coords in degrees
        coord_nod2D(2,k) = lat*rad2deg
        ! border (and land) points (0 or 1)
        index_nod2D(k)   =  MIN(1, MAX(0, 1+p_patch_3D%surface_vertex_sea_land_mask(jv,jb)))
      ENDDO
    ENDDO

    ! the rotation matrix is a modified version of r2*r1 from disp_new_vect
    cos_lonp=COS(pollon)
    sin_lonp=SIN(pollon)
    cos_latp=COS(pollat)
    sin_latp=SIN(pollat)

    rot_mat_3D(1,:) =  (/  cos_lonp*sin_latp,  sin_lonp*sin_latp, -cos_latp /)
    rot_mat_3D(2,:) =  (/ -sin_lonp ,      cos_lonp,       0.0_wp /)
    rot_mat_3D(3,:) =  (/  cos_lonp*cos_latp,  sin_lonp*cos_latp,  sin_latp /)

    ! elem2D is the number of elements = number of cells
    elem2D = p_patch%n_patch_cells
    ALLOCATE(elem2D_nodes(3,elem2D),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating coord_nod2D failed')
    ENDIF

    ! Establish the connectivity between elements and nodes. The array elem2D_nodes lists the
    ! nodes each element is made up of.
    k=0
    DO jb = p_patch%cells%all%start_block, p_patch%cells%all%end_block
      CALL get_index_range(p_patch%cells%all, jb, i_startidx_c, i_endidx_c) 
      DO jc = i_startidx_c,i_endidx_c
        k=k+1
        elem2D_nodes(1,k) =   &
          &   verts(p_patch%cells%vertex_idx(jc,jb,1),p_patch%cells%vertex_blk(jc,jb,1))
        elem2D_nodes(2,k) =   &
          &   verts(p_patch%cells%vertex_idx(jc,jb,2),p_patch%cells%vertex_blk(jc,jb,2))
        elem2D_nodes(3,k) =   &
          &   verts(p_patch%cells%vertex_idx(jc,jb,3),p_patch%cells%vertex_blk(jc,jb,3))
      ENDDO
    ENDDO

    ! allocate verts->cells connectivity index for the restructured solver loop
    ALLOCATE(nod2D_elems(6,nod2D),STAT=ist) ! 6 neighbouring cells by default, can be less for pentagons and boundary nodes
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating nod2D_elems failed')
    ENDIF
    ! initialize with zeros
    nod2D_elems = 0

    k=0
    DO jb = p_patch%cells%all%start_block, p_patch%cells%all%end_block
      CALL get_index_range(p_patch%cells%all, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c,i_endidx_c
        k=k+1
        elems(jc,jb) = k
      ENDDO
    ENDDO
    ! Establish the connectivity between elements and nodes.
    ! Array nod2D_elems lists the elements that share a given node.
    k=0
    DO jb = p_patch%verts%all%start_block, p_patch%verts%all%end_block
      CALL get_index_range(p_patch%verts%all, jb, i_startidx_v, i_endidx_v)
      DO jv = i_startidx_v, i_endidx_v
        k=k+1
        DO jc = 1, 6
          cell_index = p_patch%verts%cell_idx(jv,jb,jc)
          cell_block = p_patch%verts%cell_blk(jv,jb,jc)

          IF (cell_index > 0)                                      &
             nod2D_elems(jc,k) = elems(cell_index,cell_block)

        END DO
      END DO
    END DO

    ! Mask out halopoints
    ! This is neccesary so that Sergey's EVP routine doesn't try to calculate velocities on the
    ! halo-points. In his code ponts with index_nod2D=0 are used for calculation and the velocity
    ! at points with index_nod2D=1 (land borders) is set to zero (no-slip). The velocity at points
    ! with index_nod2D>1 will therefore not be calculated and not set to zero.
    halo_mask = 0
    DO jb = p_patch%verts%not_owned%start_block, p_patch%verts%not_owned%end_block
      CALL get_index_range(p_patch%verts%not_owned, jb, i_startidx_v, i_endidx_v)
        DO jv = i_startidx_v, i_endidx_v
          halo_mask(jv,jb) = 2
        ENDDO
    ENDDO

    buffy = RESHAPE(halo_mask(:,:), SHAPE(buffy))
    index_nod2D = index_nod2D+buffy(1:SIZE(index_nod2D))

  END SUBROUTINE ice_fem_grid_init

  !-------------------------------------------------------------------------
  !
  !> Copy basis functions arrays (bafux, bafuy) to (bafux_nod, bafuy_nod)
  !! which are indexed for each node, rather than for each element
  !! Those arrays are necessary for the OMP-parallelized version of stress2rhs_omp
  !! Where we loop over the nodes to update the rhs value in the momentum equatoins.
  !!
  SUBROUTINE basisfunctions_nod

    USE mo_ice_fem_mesh,           ONLY: nod2D, elem2D_nodes, nod2D_elems,    &
                                   & bafux, bafuy, bafux_nod, bafuy_nod

    ! Indexing
    INTEGER :: je, jn, k
    INTEGER :: ist
    ! Temporary var
    INTEGER ::  elem, elnodes(3)
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ice_fem_interface:basisfunctions_nod'
    LOGICAL :: found, kill

    ! allocate bafux_nod, bafuy_nod
    ! 6 neighbouring cells by default, can be less for pentagons and boundary nodes
    ALLOCATE(bafux_nod(6,nod2D),bafuy_nod(6,nod2D),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating bafux_nod/bafuy_nod failed')
    ENDIF
    ! initialize with zeros
    bafux_nod = 0.0_wp
    bafuy_nod = 0.0_wp
    kill = .false.
!ICON_OMP_PARALLEL_DO PRIVATE(jn,je,elem,elnodes,k,found) ICON_OMP_DEFAULT_SCHEDULE
    DO jn=1, nod2D
         IF (kill) CYCLE
         DO je=1,6
            elem=nod2D_elems(je,jn)

            IF (elem > 0) THEN
              elnodes=elem2D_nodes(:,elem)

              ! find the node for the giveN element elem, that matches jn
              ! just brute force search
              found = .false.
              DO k=1,3
                IF (jn==elnodes(k)) THEN
                    found = .true.
                    exit
                ENDIF
              END DO

              IF (found) THEN
                  bafux_nod(je,jn) = bafux(k,elem)
                  bafuy_nod(je,jn) = bafuy(k,elem)
              ELSE
                kill = .true.
!                CALL finish(routine, "FEM element-vertex connectivity inconsistency")
              ENDIF

            END IF
         END DO
    END DO
!ICON_OMP_END_PARALLEL_DO
    if (kill) &
       CALL finish(routine, "FEM element-vertex connectivity inconsistency")
  END SUBROUTINE basisfunctions_nod

  !-------------------------------------------------------------------------
  !> Calculates the Coriolis parameters at a point on the rotated grid with
  !> the new pole at (pollon, pollat).
  !-------------------------------------------------------------------------
  !!
  function calc_f_rot(lat, lon, pollat, omega)

  implicit none

  real(wp), intent(in)  :: lat, lon, pollat, omega
  real(wp)              :: calc_f_rot

  calc_f_rot = 2._wp*omega*(-cos(lon)*cos(lat)*cos(pollat)+sin(lat)*sin(pollat))

  end function calc_f_rot

  !-------------------------------------------------------------------------
  !
  !> Synchronisation routine for the FEM
  !! This replaces the routine Sergey used to synchronize arrays across CPUs. It is just a wrapper
  !! around sync_patch_array using the PRIVATE variable fem_patch as patch and reshaping the
  !! variable before and after syncing.
  !!
  SUBROUTINE exchange_nod2D(u_ice)

    REAL(wp), INTENT(INOUT) :: u_ice(fem_patch%n_patch_verts)

    ! Temporary buffer
    REAL(wp) :: u_(nproma, fem_patch%nblks_v)

    CALL copy_fem2icon(u_ice, u_, 1, 1)

    CALL sync_patch_array(SYNC_V, fem_patch, u_(:,:))

    CALL copy_icon2fem(u_, u_ice, 1, 1)
  END SUBROUTINE exchange_nod2D

  !-------------------------------------------------------------------------
  ! copy ice velocities to ICON variables
  SUBROUTINE copy_fem2icon(u_ice, u_, nlev, lev_idx, lacc)
    INTEGER, INTENT(in) :: nlev, lev_idx
    REAL(wp), INTENT(IN) :: u_ice(fem_patch%n_patch_verts)
    REAL(wp), INTENT(INOUT) :: u_(nproma, nlev, fem_patch%nblks_v)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    INTEGER :: jv, jb, npad, nlast
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) IF(lzacc)
    DO jb = 1, fem_patch%nblks_v-1
      DO jv = 1, nproma
        u_(jv, lev_idx, jb) = u_ice((jb-1)*nproma + jv)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    npad = nproma*fem_patch%nblks_v - fem_patch%n_patch_verts
    nlast = nproma - npad
    jb = fem_patch%nblks_v
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
    DO jv = 1, nlast
      u_(jv, lev_idx, jb) = u_ice((jb-1)*nproma + jv)
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
    DO jv = nlast+1,nproma
      u_(jv, lev_idx, jb) = -9999._wp
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE copy_fem2icon

  ! Reshape and copy ice velocity to FEM model variables
  SUBROUTINE copy_icon2fem(u_, u_ice, nlev, lev_idx, lacc)
    INTEGER, INTENT(in) :: nlev, lev_idx
    REAL(wp), INTENT(IN) :: u_(nproma, nlev, fem_patch%nblks_v)
    REAL(wp), INTENT(OUT) :: u_ice(fem_patch%n_patch_verts)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    LOGICAL :: lzacc
    INTEGER :: jv, jb, npad, nlast

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) IF(lzacc)
    DO jb = 1, fem_patch%nblks_v-1
      DO jv = 1, nproma
        u_ice((jb-1)*nproma + jv) = u_(jv, lev_idx, jb)
      END DO
    END DO
    !$ACC END PARALLEL LOOP
    npad = nproma*fem_patch%nblks_v - fem_patch%n_patch_verts
    nlast = nproma - npad
    jb = fem_patch%nblks_v
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
    DO jv = 1, nlast
      u_ice((jb-1)*nproma + jv) = u_(jv, lev_idx, jb)
    END DO
    !$ACC END PARALLEL LOOP
  END SUBROUTINE copy_icon2fem

  SUBROUTINE exchange_nod2Dx2(u1_ice, u2_ice, lacc)
    REAL(wp), INTENT(INOUT) :: u1_ice(fem_patch%n_patch_verts), &
         u2_ice(fem_patch%n_patch_verts)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc
    
    LOGICAL :: lzacc
    ! Temporary buffer
    REAL(wp) :: u_(nproma, 2, fem_patch%nblks_v)

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA CREATE(u_) IF(lzacc)

    CALL copy_fem2icon(u1_ice, u_, 2, 1, lacc=lzacc)
    CALL copy_fem2icon(u2_ice, u_, 2, 2, lacc=lzacc)

    CALL sync_patch_array(SYNC_V, fem_patch, u_)

    CALL copy_icon2fem(u_, u1_ice, 2, 1, lacc=lzacc)
    CALL copy_icon2fem(u_, u2_ice, 2, 2, lacc=lzacc)

    !$ACC END DATA

  END SUBROUTINE exchange_nod2Dx2

  SUBROUTINE exchange_nod2Dx3(u1_ice, u2_ice, u3_ice)
    REAL(wp), INTENT(INOUT) :: u1_ice(fem_patch%n_patch_verts), &
         u2_ice(fem_patch%n_patch_verts), u3_ice(fem_patch%n_patch_verts)

    ! Temporary buffer
    REAL(wp) :: u_(nproma, 3, fem_patch%nblks_v)

    CALL copy_fem2icon(u1_ice, u_, 3, 1)
    CALL copy_fem2icon(u2_ice, u_, 3, 2)
    CALL copy_fem2icon(u3_ice, u_, 3, 3)

    CALL sync_patch_array(SYNC_V, fem_patch, u_)

    CALL copy_icon2fem(u_, u1_ice, 3, 1)
    CALL copy_icon2fem(u_, u2_ice, 3, 2)
    CALL copy_icon2fem(u_, u3_ice, 3, 3)

  END SUBROUTINE exchange_nod2Dx3

  !-------------------------------------------------------------------------
  !
  !> Post-initialisation work for the FEM grid
  !! In order for the FEM model to pass the p-test we need to sort the lists myList_elem2D and
  !! myList_nod2D according to the global index. This is neccesary so that summation in stress2rhs
  !! (in ice_dyn_fem/ice_evp.f90) is always done in the same order on all CPUs.
  !! In addition to this we also need to set lmass_matrix to the dual area size, since lmass_matrix
  !! does not synchronise across CPUs the way it is calculated here. We also set voltriangle to the
  !! cell area, since the cell area is more accurately calculated by ICON than the FEM code.
  !!
  SUBROUTINE ice_fem_grid_post(p_patch)
    USE mo_ice_fem_types,       ONLY: lmass_matrix
    USE mo_ice_fem_mesh,        ONLY: myList_nod2D, myList_elem2D, coriolis_nod2D, voltriangle
   !USE mo_ice_fem_mesh,        ONLY: nod2D, elem2D, cos_elem2D, sin_elem2D

    USE mo_util_sort,           ONLY: quicksort

    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch

    ! Local variables
    ! Patch and ranges
    TYPE(t_subset_range), POINTER :: all_verts
    ! Indexing
    INTEGER :: i_startidx_v, i_endidx_v
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: k, jb, jv, jc

    ! Global lists
    INTEGER :: globList_nod2D(p_patch%n_patch_verts)
    INTEGER :: globList_elem2D(p_patch%n_patch_cells)

    ! Temporary variables/buffers
    REAL(wp) :: buffy_v(nproma*p_patch%nblks_v)
    REAL(wp) :: buffy_c(nproma*p_patch%alloc_cell_blocks)

  !-------------------------------------------------------------------------

    ! Copy basis functions arrays (bafux, bafuy) to (bafux_nod, bafuy_nod)
    ! Used by the OMP-parallelized loop in stress2rhs_omp
    CALL basisfunctions_nod

    all_verts => p_patch%verts%all
    
    ! Sort the list myList_nod2D using the global index
    ! Note: myList_nod2D==(1:nod2D) and myList_elem2D==(1:elem2D)
    k=0
    DO jb = 1,p_patch%nblks_v
      CALL get_index_range(all_verts, jb, i_startidx_v, i_endidx_v)
      DO jv = i_startidx_v,i_endidx_v
        k=k+1
        globList_nod2D(k) = p_patch%verts%decomp_info%glb_index(k)
      ENDDO
    ENDDO

    CALL quicksort(globList_nod2D,myList_nod2D)

    ! Sort the list myList_elem2D using the global index
    k=0
    DO jb = p_patch%cells%all%start_block, p_patch%cells%all%end_block
      CALL get_index_range(p_patch%cells%all, jb, i_startidx_c, i_endidx_c) 
      DO jc = i_startidx_c,i_endidx_c
        k=k+1
        globList_elem2D(k) = p_patch%cells%decomp_info%glb_index(k)
      ENDDO
    ENDDO

    CALL quicksort(globList_elem2D,myList_elem2D)

    ! Re-create lmass_matrix as one third of the sum of areas of triangles that share vertex n
    ! Or in ICON-terms as the dual area

    buffy_v = RESHAPE(p_patch%verts%dual_area(:,:), SHAPE(buffy_v))
    lmass_matrix = buffy_v(1:SIZE(lmass_matrix))

    ! voltriangle should also be the cell area
    buffy_c = RESHAPE(p_patch%cells%area(:,:), SHAPE(buffy_c))
    voltriangle = buffy_c(1:SIZE(voltriangle))

    ! IS NOT NECESSARY, already handled by the previously called intitialization routines
    ! But for later restructuring -- keep this way, instead of recalculating
    buffy_v = RESHAPE(p_patch%verts%f_v(:,:), SHAPE(buffy_v))
    coriolis_nod2D = buffy_v(1:SIZE(coriolis_nod2D))

  END SUBROUTINE ice_fem_grid_post

END MODULE mo_ice_fem_icon_init
