! Contains the implementation of the top and bottom ocean boundary conditions
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

MODULE mo_ocean_check_tools
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_impl_constants,     ONLY: max_char_length, sea_boundary, boundary, sea, land, land_boundary, min_dolic
  USE mo_model_domain,       ONLY: t_patch, t_patch_3D
  USE mo_ocean_nml,          ONLY: n_zlev, no_tracer
  USE mo_dynamics_config,    ONLY: nold,nnew
  USE mo_run_config,         ONLY: dtime
  USE mo_exception,          ONLY: finish, message, message_text
  USE mo_util_dbg_prnt,      ONLY: dbg_print, c_i, c_b, near_proc_id
  USE mo_grid_config,        ONLY: n_dom
  USE mo_grid_subset,        ONLY: t_subset_range, get_index_range
  USE mo_ocean_types,          ONLY: t_hydro_ocean_state
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_run_config,         ONLY: dtime, nsteps
  USE mo_math_constants,     ONLY: pi
  USE mo_mpi,                ONLY: get_my_global_mpi_id, p_pe
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: init_oce_index
  PUBLIC :: ocean_check_level_sea_land_mask
  PUBLIC :: check_ocean_subsets

CONTAINS
  !-------------------------------------------------------------------------
  !>
  !! Initialization of indices for some output on ocean variables
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE init_oce_index (patch_2d, patch_3d, pstate_oce, p_ext_data)

    TYPE(t_patch),             TARGET, INTENT(in)     :: patch_2d(n_dom)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)       :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout)  :: pstate_oce(n_dom)
    TYPE(t_external_data),     TARGET, INTENT(in)     :: p_ext_data(n_dom)

    INTEGER :: jg, jt, i, islmval, idolic
    INTEGER :: c_k, ne_b(3), ne_i(3), nc_b(3), nc_i(3), nv_b(3), nv_i(3)
    REAL(wp) :: zlon, zlat, bathy
    CHARACTER(LEN=90) :: form4ar

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_ocean_check_tools:init_oce_index'

    !CALL message(TRIM(routine), 'Start' )

    IF (n_dom > 1 ) THEN
      CALL finish(TRIM(routine), ' N_DOM > 1 is not allowed')
    END IF
    jg = n_dom

    ! set time level to old level (nold(1)=3)
    jt = nold(jg)

    c_k = 1

    WRITE(message_text,'(a)') 'Output diverse variables:'
    CALL message (TRIM(routine), message_text)

    WRITE(message_text,'(3(a,i8))')          &
      & '  nzlev   =', n_zlev,             &
      & ', n_dom   =', n_dom,              &
      & ', nproma  =', nproma
    CALL message ('', message_text)
    WRITE(message_text,'(3(a,i8))')          &
      & '  nblks_c =', patch_2d(jg)%nblks_c, &
      & ', nblks_e =', patch_2d(jg)%nblks_e, &
      & ', nblks_v =', patch_2d(jg)%nblks_v
    CALL message ('', message_text)
    WRITE(message_text,'(3(a,i8))')          &
      & '  nold(jg)=', nold(jg),           &
      & ', nnew(jg)=', nnew(jg),           &
      & ', Nsteps  =', nsteps
    CALL message ('', message_text)
    WRITE(message_text,'(5(a,i8))')          &
      & '  SEA=',     sea,                 &
      & ', SEA_BOUNDARY=',sea_boundary,    &
      & ', BOUNDARY=',boundary,            &
      & ', LAND_BOUNDARY=',land_boundary,  &
      & ', LAND=,',   land!,                &
    CALL message ('', message_text)
    WRITE(message_text,'(4(a,g18.6))')       &
      & '  Time Step=',dtime
    CALL message ('', message_text)

    !------------------------------------------------------------------
    ! Test output of indices
    !------------------------------------------------------------------

    ! ! values for the blocking: cells
    ! rl_start = 1
    ! rl_end = min_rlcell
    ! i_startblk = patch_2D(jg)%cells%start_blk(rl_start,1)
    ! i_endblk   = patch_2D(jg)%cells%end_blk(rl_end,1)
    ! ! first block
    ! CALL get_indices_c(patch_2D(jg), i_startblk, i_startblk, i_endblk, i_startidxf, i_endidxf, &
    !   &                            rl_start, rl_end)
    ! ! last block
    ! CALL get_indices_c(patch_2D(jg), i_endblk, i_startblk, i_endblk, i_startidxl, i_endidxl, &
    !   &                            rl_start, rl_end)
    ! WRITE(message_text,'(8(a,i5))')                         &
    !   &   ' Cells: rl_start=',rl_start,' rl_end=',rl_end,   &
    !   &   ' stblk=',i_startblk,' endblk=',i_endblk,         &
    !   &   ' fst blk: sidx=',i_startidxf,' eidx=',i_endidxf, &
    !   &   ' lst blk: sidx=',i_startidxl,' eidx=',i_endidxl
    ! CALL message ('', message_text)

    ! ! values for the blocking: edges
    ! rl_end = min_rledge
    ! i_startblk = patch_2D(jg)%edges%start_blk(rl_start,1)
    ! i_endblk   = patch_2D(jg)%edges%end_blk(rl_end,1)
    ! ! first block
    ! CALL get_indices_e(patch_2D(jg), i_startblk, i_startblk, i_endblk, i_startidxf, i_endidxf, &
    !   &                            rl_start, rl_end)
    ! ! last block
    ! CALL get_indices_e(patch_2D(jg), i_endblk, i_startblk, i_endblk, i_startidxl, i_endidxl, &
    !   &                            rl_start, rl_end)
    ! WRITE(message_text,'(8(a,i5))')                         &
    !   &   ' Edges: rl_start=',rl_start,' rl_end=',rl_end,   &
    !   &   ' stblk=',i_startblk,' endblk=',i_endblk,         &
    !   &   ' fst blk: sidx=',i_startidxf,' eidx=',i_endidxf, &
    !   &   ' lst blk: sidx=',i_startidxl,' eidx=',i_endidxl
    ! CALL message ('', message_text)

    ! ! values for the blocking: verts
    ! rl_end = min_rlvert
    ! i_startblk = patch_2D(jg)%verts%start_blk(rl_start,1)
    ! i_endblk   = patch_2D(jg)%verts%end_blk(rl_end,1)
    ! ! first block
    ! CALL get_indices_v(patch_2D(jg), i_startblk, i_startblk, i_endblk, i_startidxf, i_endidxf, &
    !   &                            rl_start, rl_end)
    ! ! last block
    ! CALL get_indices_v(patch_2D(jg), i_endblk, i_startblk, i_endblk, i_startidxl, i_endidxl, &
    !   &                            rl_start, rl_end)
    ! WRITE(message_text,'(8(a,i5))')                         &
    !   &   ' Verts: rl_start=',rl_start,' rl_end=',rl_end,   &
    !   &   ' stblk=',i_startblk,' endblk=',i_endblk,         &
    !   &   ' fst blk: sidx=',i_startidxf,' eidx=',i_endidxf, &
    !   &   ' lst blk: sidx=',i_startidxl,' eidx=',i_endidxl
    ! CALL message ('', message_text)

    !------------------------------------------------------------------
    ! Check parameters
    !------------------------------------------------------------------

    IF (p_pe /= near_proc_id) RETURN

    ! slm and coordinates of this point:
    islmval = patch_3d%lsm_c(c_i,c_k,c_b)
    idolic  = patch_3d%p_patch_1d(1)%dolic_c  (c_i,    c_b)
    bathy   = p_ext_data(jg)%oce%bathymetry_c (c_i,c_b)
    zlat = patch_2d(jg)%cells%center(c_i,c_b)%lat * 180.0_wp / pi
    zlon = patch_2d(jg)%cells%center(c_i,c_b)%lon * 180.0_wp / pi

    ! output format
    99 FORMAT(a,i4,a,i4,a,i3,a,i3,3(a,f9.2))
    97 FORMAT(a,i1,a,i4,a,i4,a,i3,a,i3,3(a,f9.2))
    form4ar = '(4(a,g20.9))'

    CALL message (TRIM(routine), 'Conditions at test cell (C), including bathymetry:')
    WRITE(message_text,99) ' Cell C: block=',c_b,'  index=',c_i,             &
      & '  lsm_c=', islmval,'  dolic_c=',idolic,'  bathy_c=', bathy, &
      & '  lat=',zlat,'  lon=',zlon
    CALL message (' ', message_text)
    IF(no_tracer>=1)THEN
      WRITE(message_text,form4ar)  &
        & ' Elev. h at Cell C    =', pstate_oce(jg)%p_prog(jt)%h(c_i,c_b),            &
        & '  Tracer 1 =', pstate_oce(jg)%p_prog(jt)%tracer(c_i,c_k,c_b,1), &
        & '  Test level at cell C = ',REAL(c_k,wp)
      CALL message (' ', message_text)
    ELSEIF(no_tracer==0)THEN
      WRITE(message_text,form4ar)  &
        & ' Elev. h at Cell C    =', pstate_oce(jg)%p_prog(jt)%h(c_i,c_b),            &
      !&                  '  Tracer 1 =', pstate_oce(jg)%p_prog(jt)%tracer(c_i,c_k,c_b,1), &
        & '  Test level (jk) at cell C = ',REAL(c_k,wp)
      CALL message (' ', message_text)
    ENDIF

    !------------------------------------------------------------------
    ! find and print corresponding edges/verts of test cell
    !------------------------------------------------------------------

    DO i = 1, 3 ! 3 edges of cell C at (ne_i,ne_b)
      ! slm and coordinates of edges
      ne_b(i)=patch_2d(jg)%cells%edge_blk(c_i,c_b,i)
      ne_i(i)=patch_2d(jg)%cells%edge_idx(c_i,c_b,i)
      islmval = patch_3d%lsm_e  (ne_i(i),c_k,ne_b(i))
      idolic  = patch_3d%p_patch_1d(1)%dolic_e    (ne_i(i),ne_b(i))
      bathy   = p_ext_data(jg)%oce%bathymetry_e (ne_i(i),ne_b(i))
      zlat    = patch_2d(jg)%edges%center         (ne_i(i),ne_b(i))%lat * 180.0_wp / pi
      zlon    = patch_2d(jg)%edges%center         (ne_i(i),ne_b(i))%lon * 180.0_wp / pi
      ! output
      WRITE(message_text,97) ' Edge E',i,' block=',ne_b(i),'  index=',ne_i(i),              &
        & '  lsm_e=', islmval,'  dolic_e=',idolic,'  bathy_e=', bathy, &
        & '  lat=',zlat,'  lon=',zlon
      CALL message (' ', message_text)
    END DO

    DO i = 1, 3 ! 3 vertices of cell C at (nv_i,nv_b)
      ! slm and coordinates of vertices
      nv_b(i)=patch_2d(jg)%cells%vertex_blk(c_i,c_b,i)
      nv_i(i)=patch_2d(jg)%cells%vertex_idx(c_i,c_b,i)
      islmval = patch_3d%lsm_c(c_i,c_k,c_b)
      idolic  = patch_3d%p_patch_1d(1)%dolic_c  (c_i,    c_b)
      bathy   = p_ext_data(jg)%oce%bathymetry_c (c_i,c_b)
      zlat    = patch_2d(jg)%edges%center         (nv_i(i),nv_b(i))%lat * 180.0_wp / pi
      zlon    = patch_2d(jg)%edges%center         (nv_i(i),nv_b(i))%lon * 180.0_wp / pi
      ! output
      WRITE(message_text,97) ' Vert V',i,' block=',nv_b(i),'  index=',nv_i(i),              &
        & '  lsm_c=', islmval,'  dolic_c=',idolic,'  bathy_c=', bathy, &
        & '  lat=',zlat,'  lon=',zlon
      CALL message (' ', message_text)
    END DO

    DO i = 1, 3 ! 3 neighbours of cell C at (nc_i,nc_b)
      ! slm and coordinates of neighbouring cells
      ! #slo# - careful at index-boundaries: in ocean grid exist cells without neighbours !!
      nc_b(i)=patch_2d(jg)%cells%neighbor_blk(c_i,c_b,i)
      nc_i(i)=patch_2d(jg)%cells%neighbor_idx(c_i,c_b,i)
      IF ( nc_i(i) == 0 .OR. nc_b(i) == 0) THEN
        nc_i(i) = c_i
        nc_b(i) = c_b
        WRITE(message_text,'(a)') ' Neighbor Cell is on LAND - NOT DEFINED'
      ELSE
        islmval = patch_3d%lsm_c  (nc_i(i),c_k,nc_b(i))
        idolic  = patch_3d%p_patch_1d(1)%dolic_c    (nc_i(i),    nc_b(i))
        bathy   = p_ext_data(jg)%oce%bathymetry_c (nc_i(i),nc_b(i))
        zlat    = patch_2d(jg)%cells%center         (nc_i(i),nc_b(i))%lat * 180.0_wp / pi
        zlon    = patch_2d(jg)%cells%center         (nc_i(i),nc_b(i))%lon * 180.0_wp / pi
        WRITE(message_text,97) ' Neighbor  C',i,' =',nc_b(i),'  index=',nc_i(i),            &
          & '  lsm_c=', islmval,'  dolic_c=',idolic,'  bathy_c=', bathy, &
          & '  lat=',zlat,'  lon=',zlon
      END IF
      ! output
      CALL message (' ', message_text)
    END DO

  END SUBROUTINE init_oce_index
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
!<Optimize:inUse>
  SUBROUTINE ocean_check_level_sea_land_mask( patch_3d )
    !
    TYPE(t_patch_3D ),TARGET, INTENT(IN):: patch_3d

    INTEGER :: block, idx, start_idx, end_idx, level
    INTEGER :: cell1_idx, cell1_blk, cell2_idx, cell2_blk
    INTEGER :: edge_levels, cell1_levels, cell2_levels
    TYPE(t_subset_range), POINTER :: all_cells, edges_in_domain
    TYPE(t_patch), POINTER        :: patch_2d
    CHARACTER(len=*), PARAMETER :: method_name='mo_ocean_check_tools:ocean_check_level_sea_land_mask'
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2D(1)
    all_cells => patch_2d%cells%all
    edges_in_domain => patch_2d%edges%in_domain
    !-----------------------------------------------------------------------
    ! check cells
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, start_idx, end_idx)
      DO idx = start_idx, end_idx

        DO level=1, patch_3d%p_patch_1d(1)%dolic_c(idx, block)
          IF(patch_3d%lsm_c(idx, level, block) > sea_boundary) THEN
            write(0,*) " dolic_c=", patch_3d%p_patch_1d(1)%dolic_c(idx, block), " at level=", level, &
              & "lsm_c= ", patch_3d%lsm_c(idx, level, block)
            CALL finish(method_name,"Inconsistent cell levels vs sea_land_mask")
          ENDIF
        ENDDO

        DO level=patch_3d%p_patch_1d(1)%dolic_c(idx, block)+1, n_zlev
          IF(patch_3d%lsm_c(idx, level, block) < boundary)THEN
            write(0,*) " dolic_c=", patch_3d%p_patch_1d(1)%dolic_c(idx, block), " at level=", level, &
              & "lsm_c= ", patch_3d%lsm_c(idx, level, block)
            CALL finish(method_name,"Inconsistent cell levels vs sea_land_mask")
          ENDIF
        ENDDO

      ENDDO
    ENDDO

    !-----------------------------------------------------------------------
    ! check edges
    DO block = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, block, start_idx, end_idx)
      DO idx = start_idx, end_idx

        cell1_idx = patch_2d%edges%cell_idx(idx, block, 1)
        cell1_blk = patch_2d%edges%cell_blk(idx, block, 1)
        cell2_idx = patch_2d%edges%cell_idx(idx, block, 2)
        cell2_blk = patch_2d%edges%cell_blk(idx, block, 2)

        ! check sea
        DO level=1, patch_3d%p_patch_1d(1)%dolic_e(idx, block)
          IF(patch_3d%lsm_e(idx, level, block) > sea_boundary) THEN
            write(0,*) " dolic_e=", patch_3d%p_patch_1d(1)%dolic_e(idx, block), " at level=", level, &
              & "lsm_e= ", patch_3d%lsm_e(idx, level, block)
            CALL finish(method_name,"Inconsistent edge levels vs sea_land_mask")
          ENDIF

          ! check neigboring cells
          IF (cell1_idx < 1 .or. cell2_idx < 1) &
            & CALL finish(method_name,"Sea edge with one cell missing")
          IF (patch_3d%lsm_c(cell1_idx, level, cell1_blk) >= boundary .or. &
              & patch_3d%lsm_c(cell2_idx, level, cell2_blk) >= boundary) THEN
            write(0,*) " at level=", level
            CALL finish(method_name,"Sea edge with one cell land")
          ENDIF

        ENDDO ! levels

        ! check depth consistancy
        IF (patch_3d%p_patch_1d(1)%dolic_e(idx, block) > 0) THEN
          edge_levels = patch_3d%p_patch_1d(1)%dolic_e(idx, block)
          cell1_levels = patch_3d%p_patch_1d(1)%dolic_c(cell1_idx, cell1_blk)
          cell2_levels = patch_3d%p_patch_1d(1)%dolic_c(cell2_idx, cell2_blk)
          IF (edge_levels /= MIN(cell1_levels, cell2_levels)) &
            CALL finish(method_name,"edge levels /= min cell levels")

          IF (patch_3D%p_patch_1D(1)%prism_thick_e(idx, edge_levels, block) /= &
            MIN(patch_3D%p_patch_1D(1)%prism_thick_c(cell1_idx, edge_levels, cell1_blk), &
                patch_3D%p_patch_1D(1)%prism_thick_c(cell2_idx, edge_levels, cell2_blk))) &
            CALL finish(method_name,"edge depth /= min cell depths")
        ENDIF
        

        ! check land
        DO level=patch_3d%p_patch_1d(1)%dolic_e(idx, block)+1, n_zlev
          IF(patch_3d%lsm_e(idx, level, block) < boundary)THEN
            write(0,*) " dolic_e=", patch_3d%p_patch_1d(1)%dolic_e(idx, block), " at level=", level, &
              & "lsm_e= ", patch_3d%lsm_e(idx, level, block)
            CALL finish(method_name,"Inconsistent edge levels vs sea_land_mask")
          ENDIF

          ! check neigboring cells
          IF (cell1_idx > 0 .and. cell2_idx > 0) THEN
            IF (patch_3d%lsm_c(cell1_idx, level, cell1_blk) < boundary .and. &
                & patch_3d%lsm_c(cell2_idx, level, cell2_blk) < boundary) THEN
              write(0,*) " at level=", level
              CALL finish(method_name,"Land edge with two cells sea")
            ENDIF
          ENDIF

        ENDDO

      ENDDO
    ENDDO

  END SUBROUTINE ocean_check_level_sea_land_mask
  !-------------------------------------------------------------------------

  !------------------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE check_ocean_subsets(patch_3d)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d

    TYPE(t_patch),     POINTER :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells, all_edges, all_verts

    INTEGER :: BLOCK, startidx, endidx, idx
    !-----------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells   => patch_2d%cells%ALL
    all_edges   => patch_2d%edges%ALL
    all_verts   => patch_2d%verts%ALL

    !    IF (patch_2D%edges%in_domain%no_of_holes > 0) THEN
    !      DO block = patch_2D%edges%in_domain%start_block, patch_2D%edges%in_domain%end_block
    !        CALL get_index_range(patch_2D%edges%in_domain, block, startidx, endidx)
    !        DO idx = startidx, endidx
    !          IF (patch_2D%edges%decomp_info%halo_level(idx, block) /= 0) THEN
    !            write(0,*) get_my_global_mpi_id(), ":", idx, block, " edge in hole ", patch_3D%surface_edge_sea_land_mask(idx, block), &
    !              & patch_2D%edges%decomp_info%halo_level(idx, block)
    !          ELSE
    !            write(0,*) get_my_global_mpi_id(), ":", idx, block, " edge in set ", patch_3D%surface_edge_sea_land_mask(idx, block), &
    !              & patch_2D%edges%decomp_info%halo_level(idx, block)
    !          ENDIF
    !        ENDDO
    !      ENDDO
    !
    !      CALL finish("patch%edges%in_domain", "no_of_holes > 0, gmres for the ocean requires no_of_holes=0")
    !    ENDIF

    !    IF (patch_2D%cells%in_domain%no_of_holes > 0 .and. .false.) THEN
    IF (patch_2d%cells%in_domain%no_of_holes > 0 ) THEN
      DO BLOCK = patch_2d%cells%in_domain%start_block, patch_2d%cells%in_domain%end_block
        CALL get_index_range(patch_2d%cells%in_domain, BLOCK, startidx, endidx)
        DO idx = startidx, endidx
          IF (patch_2d%cells%decomp_info%halo_level(idx, BLOCK) /= 0) THEN
            WRITE(0,*) get_my_global_mpi_id(), ":", idx, BLOCK, " cell in hole ", patch_3d%surface_cell_sea_land_mask(idx, BLOCK), &
              & patch_2d%cells%decomp_info%halo_level(idx, BLOCK)
          ELSE
            WRITE(0,*) get_my_global_mpi_id(), ":", idx, BLOCK, " cell in set ", patch_3d%surface_cell_sea_land_mask(idx, BLOCK), &
              & patch_2d%cells%decomp_info%halo_level(idx, BLOCK)
          ENDIF
        ENDDO
      ENDDO

      CALL finish("patch%cells%in_domain", "no_of_holes > 0, gmres for the ocean requires no_of_holes=0")
    ENDIF

  END SUBROUTINE check_ocean_subsets
  !------------------------------------------------------------------------------------


END MODULE mo_ocean_check_tools
