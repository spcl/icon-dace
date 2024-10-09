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

! The module <i>mo_util_dbg_prnt</i> prints out max and min as well as single
! cell and neighbouring values of 2- and 3-dim arrays of the icon core for
! debug purposes

!----------------------------
#include "icon_definitions.inc"
!----------------------------
MODULE mo_util_dbg_prnt
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                   ONLY: wp
  USE mo_mpi,                    ONLY: my_process_is_stdio, p_pe, get_my_mpi_all_id !, p_comm_work, p_bcast
  USE mo_io_units,               ONLY: nerr
  USE mo_parallel_config,        ONLY: nproma, p_test_run
  USE mo_impl_constants,         ONLY: max_char_length
  USE mo_timer,                  ONLY: timer_start, timer_stop, timer_dbg_prnt, timers_level
  USE mo_sync,                   ONLY: sync_c, sync_patch_array, global_max
  USE mo_grid_subset,            ONLY: t_subset_range, get_index_range
  USE mo_dbg_nml,                ONLY: str_mod_tst, dim_mod_tst, dbg_lon_in, dbg_lat_in, &
    & idbg_mxmn, idbg_val, idbg_slev, idbg_elev, &
    & idbg_idx, idbg_blk
  USE mo_math_constants,         ONLY: pi
  USE mo_exception,              ONLY: message, message_text
  USE mo_model_domain,           ONLY: t_patch
  USE mo_grid_subset,            ONLY: t_subset_range
  USE mo_statistics,             ONLY: global_minmaxmean
  USE mo_icon_comm_interface,    ONLY: icon_comm_barrier
  USE mo_model_domain,           ONLY: t_patch_3d
  
  IMPLICIT NONE
  
  PRIVATE

  ! Public subroutines:
  PUBLIC :: init_dbg_index
  PUBLIC :: dbg_print
  PUBLIC :: debug_verts_on_edges
  PUBLIC :: debug_print_MaxMinMean, debug_printValue

  ! Public variables: should be removed!
  PUBLIC :: c_i, c_b, nc_i, nc_b, near_proc_id
  !PUBLIC :: v_subdom_cell !, v_suball_cell, v_subset_edge  !  part of subset to store
  
  ! indices of cells and neighbours for debug output at single cell
  INTEGER :: c_b, c_i, ne_b(3), ne_i(3), nc_b(3), nc_i(3), nv_b(3), nv_i(3), near_proc_id
  INTEGER :: loc_nblks_c, loc_nblks_e, loc_nblks_v
  LOGICAL :: p_test_run_bac
  !TYPE(t_subset_range) :: v_subdom_cell !, v_suball_cell, v_subset_edge  !  part of subset to store

  INTERFACE dbg_print
    MODULE PROCEDURE dbg_print_2d
    MODULE PROCEDURE dbg_print_3d
  END INTERFACE
  
  
CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  SUBROUTINE debug_verts_on_edges( vertex_variable, patch_3d, variable_name )
    !
    REAL(wp) :: vertex_variable(:,:,:)
    TYPE(t_patch_3d ),TARGET, INTENT(in):: patch_3d
    CHARACTER(LEN=*) :: variable_name
    
    INTEGER :: start_edge_index, end_edge_index, jb, jk, edge_index
    TYPE(t_subset_range), POINTER :: edges_owned
    TYPE(t_patch), POINTER :: patch_2d
    REAL(wp), POINTER :: edge_added_values(:,:,:)
    CHARACTER(LEN=*), PARAMETER :: method_name='debug_verts_on_edges'
    !-----------------------------------------------------------------------
    patch_2d   => patch_3d%p_patch_2d(1)
    edges_owned => patch_2d%edges%owned
    ALLOCATE(edge_added_values(nproma,edges_owned%max_vertical_levels,patch_3d%p_patch_2d(1)%nblks_e))
    !-------------------------------------------------------------
    DO jb = edges_owned%start_block, edges_owned%end_block
      CALL get_index_range(edges_owned, jb, start_edge_index, end_edge_index)
      DO edge_index = start_edge_index, end_edge_index
        DO jk = 1, edges_owned%max_vertical_levels
          edge_added_values(edge_index,jk,jb) = &
            & vertex_variable(patch_2d%edges%vertex_idx(edge_index, jb, 1), jk, &
            & patch_2d%edges%vertex_blk(edge_index, jb, 1)) + &
            & vertex_variable(patch_2d%edges%vertex_idx(edge_index, jb, 2), jk, &
            & patch_2d%edges%vertex_blk(edge_index, jb, 2))
        END DO
      END DO
    END DO
    
    CALL dbg_print(variable_name, edge_added_values, method_name, 4, &
      & in_subset = patch_2d%edges%owned)
    
    DEALLOCATE(edge_added_values)
    
  END SUBROUTINE debug_verts_on_edges
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !! Initialization of indices for debug output
  !
  ! TODO: parallelize
  !
  SUBROUTINE init_dbg_index (ppatch)

    TYPE(t_patch),             TARGET, INTENT(inout) :: ppatch

    INTEGER :: i
    REAL(wp) :: zlon, zlat, zarea, zlength
    
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_util_dbg_prnt:init_dbg_index'
    
    CALL message(TRIM(routine), 'Start' )
    
    ! test submit via eclipse workshop / 2013-02-11
    
    ! fill subset for use in dbg_print without passing the patch in every call
    !    v_subdom_cell = ppatch%cells%in_domain
    
    ! local module variables for check of cells/edges/verts
    loc_nblks_c =ppatch%nblks_c
    loc_nblks_e =ppatch%nblks_e
    loc_nblks_v =ppatch%nblks_v
    
    !  For a correct dicision of cells/edges/verts the number of points on domain would be better
    !    but is not available as local dimension
    !  In order to keep a difference in number of blocks please use a low number for nproma
    !    i.e nproma=1 for debugging
    !  loc_patch_c =ppatch%n_patch_cells
    !  loc_patch_e =ppatch%n_patch_edges
    !  loc_patch_v =ppatch%n_patch_verts

    ! module index/block for one cell output
    IF ((idbg_idx /= 0 ) .OR. (idbg_blk /= 0 )) THEN
      c_i = idbg_idx
      c_b = idbg_blk
      near_proc_id = p_pe
    ELSE
      ! search for block/index of debug output cell at lat/lon
      ! given by namelist dbg_index_nml - not yet parallelized
      CALL find_latlonindex (ppatch, dbg_lat_in, dbg_lon_in, c_i, c_b, near_proc_id)
    END IF

    IF (p_pe /= near_proc_id) RETURN

    zlat = ppatch%cells%center(c_i,c_b)%lat * 180.0_wp / pi
    zlon = ppatch%cells%center(c_i,c_b)%lon * 180.0_wp / pi

    !------------------------------------------------------------------
    ! print test cell
    !------------------------------------------------------------------

    ! output format
    99 FORMAT(     2(a,i4),2(a,f9.2),a,f13.2)
    97 FORMAT(a,i1,2(a,i4),2(a,f9.2),a,f13.2)

    IF (idbg_val > 0) THEN
      zarea = ppatch%cells%area(c_i,c_b)*1.0e-6_wp ! in km2
      CALL message (TRIM(routine), 'Conditions at test cell (C), and edges/verts/neighbors:', all_print=.TRUE.)
      WRITE(message_text,99) ' Cell C: block=',c_b,'  index=',c_i,               &
        & '  lat=',zlat,'  lon=',zlon,                        &
        & '  cell-area  =', zarea
      CALL message (' ', message_text, all_print=.TRUE.)
    END IF

    !------------------------------------------------------------------
    ! find and print corresponding edges/verts/neighbors of test cell
    !------------------------------------------------------------------

    DO i = 1, 3 ! 3 edges of cell C at (ne_i,ne_b)
      ne_b(i) = ppatch%cells%edge_blk(c_i,c_b,i)
      ne_i(i) = ppatch%cells%edge_idx(c_i,c_b,i)

      IF (idbg_val > 0) THEN
        zlat    = ppatch%edges%center(ne_i(i),ne_b(i))%lat * 180.0_wp / pi
        zlon    = ppatch%edges%center(ne_i(i),ne_b(i))%lon * 180.0_wp / pi
        zlength = ppatch%edges%primal_edge_length(ne_i(i),ne_b(i))*0.001_wp  ! in km
        ! output
        WRITE(message_text,97) ' Edge E',i,' block=',ne_b(i),'  index=',ne_i(i), &
          & '  lat=',zlat,'  lon=',zlon,                      &
          & '  edge-length=',zlength
        CALL message (' ', message_text, all_print=.TRUE.)
      END IF
    END DO

    DO i = 1, 3 ! 3 vertices of cell C at (nv_i,nv_b)
      nv_b(i) = ppatch%cells%vertex_blk(c_i,c_b,i)
      nv_i(i) = ppatch%cells%vertex_idx(c_i,c_b,i)
      IF (idbg_val > 0) THEN
        zlat    = ppatch%edges%center(nv_i(i),nv_b(i))%lat * 180.0_wp / pi
        zlon    = ppatch%edges%center(nv_i(i),nv_b(i))%lon * 180.0_wp / pi
        ! output
        WRITE(message_text,97) ' Vert V',i,' block=',nv_b(i),'  index=',nv_i(i), &
          & '  lat=',zlat,'  lon=',zlon
        CALL message (' ', message_text, all_print=.TRUE.)
      END IF
    END DO

    DO i = 1, 3 ! 3 neighbours of cell C at (nc_i,nc_b)
      nc_b(i)=ppatch%cells%neighbor_blk(c_i,c_b,i)
      nc_i(i)=ppatch%cells%neighbor_idx(c_i,c_b,i)

      IF ( nc_i(i) == 0 .OR. nc_b(i) == 0) THEN
        nc_i(i) = c_i
        nc_b(i) = c_b
        WRITE(message_text,'(a)') ' Neighbor Cell is NOT DEFINED'
      ELSE
        zlat = ppatch%cells%center(nc_i(i),nc_b(i))%lat * 180.0_wp / pi
        zlon = ppatch%cells%center(nc_i(i),nc_b(i))%lon * 180.0_wp / pi
        zarea= ppatch%cells%area(nc_i(i),nc_b(i))*1.0e-6_wp  ! in km2
        WRITE(message_text,97) ' Neighbor  C',i,' =',nc_b(i),'  index=',nc_i(i),  &
          & '  lat=',zlat,'  lon=',zlon,                       &
          & '  cell-area  =', zarea
      END IF

      ! output
      IF (idbg_val > 0) THEN
        CALL message (' ', message_text, all_print=.TRUE.)
      END IF
    END DO

  END SUBROUTINE init_dbg_index
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Search for a cell center at given longitude and latitude
  !! provided in namelist dbg_index_nml
  !!
  !! TODO: parallelize
  !
  SUBROUTINE find_latlonindex (ppatch, plat_in, plon_in, iidx, iblk, proc_id)

    TYPE(t_patch), TARGET, INTENT(inout)     :: ppatch
    REAL(wp),              INTENT(in)     :: plat_in       ! cell latitude to search for
    REAL(wp),              INTENT(in)     :: plon_in       ! cell longitude to search for
    INTEGER,               INTENT(out)    :: iidx          ! index of nearest cell
    INTEGER,               INTENT(out)    :: iblk          ! block of nearest cell
    INTEGER,               INTENT(out)    :: proc_id       ! process where nearest cell is found
    
    INTEGER :: jb, jc, i_startidx, i_endidx, mproc_id   !, mpi_comm
    REAL(wp) :: zlon, zlat, zdist, zdist_cmp, xctr
    REAL(wp) :: zdst_c(nproma,ppatch%nblks_c)
    TYPE(t_subset_range), POINTER :: owned_cells !,cells_in_domain, all_cells
    
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_util_dbg_prnt:find_latlonindex'
    
    CALL message(TRIM(routine), 'Start' )
    
    !all_cells       => ppatch%cells%all
    !cells_in_domain => ppatch%cells%in_domain
    owned_cells      => ppatch%cells%owned
    
    ! initial distance to compare
    zdist_cmp   = 100000.0_wp
    zdst_c(:,:) = 100000.0_wp
    
    !  loop over owned cells
    DO jb = owned_cells%start_block, owned_cells%end_block
      CALL get_index_range(owned_cells, jb, i_startidx, i_endidx)
      
      DO jc = i_startidx, i_endidx
        
        zlat    = ppatch%cells%center(jc,jb)%lat * 180.0_wp / pi
        zlon    = ppatch%cells%center(jc,jb)%lon * 180.0_wp / pi
        
        zdist       = SQRT((zlat-plat_in)*(zlat-plat_in) + (zlon-plon_in)*(zlon-plon_in))
        zdst_c(jc,jb) = zdist
        IF (zdist < zdist_cmp) THEN
          iblk = jb
          iidx = jc
          zdist_cmp = zdist
        END IF
        
      END DO
    END DO
    
    CALL sync_patch_array(sync_c, ppatch, zdst_c(:,:))
    
    ! find PE with minimum distance
    ! disable p_test_run since global_max will be different
    
    p_test_run_bac = p_test_run
    p_test_run = .FALSE.
    
    ! comparing zdist_cmp over all domains
    ! local variable mproc_id must be pre-set to p_pe in all domains
    mproc_id = p_pe
    !write(20+p_pe,*) ' 0max: mproc_id=',mproc_id,'  p_pe=',p_pe,'  zdst=',zdst_c(iidx,iblk),'  idx/blk=',iidx,iblk,' cmp=',zdist_cmp
    xctr = global_max(-zdist_cmp, mproc_id)
    !write(20+p_pe,*) ' 1max: mproc_id=',mproc_id,'  p_pe=',p_pe,'  zdst=',zdst_c(iidx,iblk),'  idx/blk=',iidx,iblk,' cmp=',zdist_cmp
    
    ! broadcast indices and proc_id to global, mainly for output - not the best, better write info directly
    !mpi_comm = p_comm_work
    !CALL p_bcast(iidx     , mproc_id, mpi_comm)
    !CALL p_bcast(iblk     , mproc_id, mpi_comm)
    !CALL p_bcast(zdist_cmp, mproc_id, mpi_comm)
    !CALL p_bcast(zlat_min , mproc_id, mpi_comm)  !  must be set accordingly above
    !CALL p_bcast(zlon_min , mproc_id, mpi_comm)  !  must be set accordingly above
    
    proc_id = mproc_id
    p_test_run = p_test_run_bac
    
    99 FORMAT(3a,i4,a,i4,3(a,f9.3))
    98 FORMAT(2a,3(a,f9.3))
    
    ! write info directly by mproc_id: needs barrier to avoid merging messages from pe_io and mproc_id
    ! write(0,*) get_my_mpi_all_id(), "enter icon_comm_barrier:"
    CALL icon_comm_barrier(for_patch=ppatch)
    CALL icon_comm_barrier(for_patch=ppatch) ! twice since the pipelining of output is still sometimes out of order
    IF (p_pe .EQ. mproc_id) THEN
      !IF (my_process_is_stdio()) THEN
      zlat = ppatch%cells%center(iidx,iblk)%lat * 180.0_wp / pi
      zlon = ppatch%cells%center(iidx,iblk)%lon * 180.0_wp / pi
      WRITE(0,98) ' ',TRIM(routine),' Found  cell nearest to         latitude=', plat_in,'  longitude=',plon_in
      WRITE(0,99) ' ',TRIM(routine),' Found  block=',iblk,'  index=',iidx,'  latitude=',zlat,'  longitude=',zlon
      WRITE(0,'(3a,i3,a,f9.3)') ' ',TRIM(routine),' FOUND: proc_id for nearest cell is=',mproc_id, &
        & '; distance in degrees =', zdist_cmp
      !  WRITE(125,98) ' ',TRIM(routine),' Found  cell nearest to         latitude=', plat_in,'  longitude=',plon_in
      !  WRITE(125,99) ' ',TRIM(routine),' Found  block=',iblk,'  index=',iidx,'  latitude=',zlat,'  longitude=',zlon
      !  WRITE(125,'(3a,i3,a,f9.3)') ' ',TRIM(routine),' FOUND: proc_id for nearest cell is=',mproc_id, &
      !    &                  '; distance in degrees =', zdist_cmp
      ! WRITE(125,'(3a,i3,a,i3,a,f9.3)') ' ',TRIM(routine),' FOUND: at block=',iblk,' index=',iidx, &
      !   &                  ' distance in degrees =', zdist_cmp
      ! WRITE(125,'(3a,2i3,a,f9.3)') ' ',TRIM(routine),' FOUND: using MINLOC: at idx/blk=', &
      !   &                  MINLOC(zdst_c(:,:)),' distance in degrees =',MINVAL(zdst_c(:,:))
    END IF
    CALL icon_comm_barrier(for_patch=ppatch)
    CALL icon_comm_barrier(for_patch=ppatch) ! twice since the pipelining of output is still sometimes out of order
    ! write(0,*) get_my_mpi_all_id(), "leave icon_comm_barrier:"
    
  END SUBROUTINE find_latlonindex
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !! Print out min and max or a specific cell value and neighbors of a 3-dim array.
  !!
  !! Reduce writing effort for a simple print.
  !! The amount of prints is controlled by comparison of a fixed level of detail
  !! for output (inDetail_level) with variables idbg_mxmn/idbg_val  that are
  !! given via namelist dbg_index_nml
  !
  SUBROUTINE dbg_print_3d( description, p_array, place, inDetail_level, in_subset )
    
    CHARACTER(LEN=*),      INTENT(in) :: description    ! description of array
    REAL(wp),              INTENT(in) :: p_array(:,:,:) ! 3-dim array for debugging
    CHARACTER(LEN=*),      INTENT(in) :: place    ! defined string for source of current array
    INTEGER,               INTENT(in) :: inDetail_level    ! source level from module for print output
    TYPE(t_subset_range),  TARGET, OPTIONAL :: in_subset
    
    ! local variables
    CHARACTER(LEN=27) ::  strout
    CHARACTER(LEN=12) ::  strmod
    INTEGER ::  slev, elev, elev_val, elev_mxmn
    INTEGER ::  iout, icheck_str_mod, jstr, i, jk, nlev, ndimblk
    REAL(wp)          :: minmaxmean(3)

    IF ((idbg_val < inDetail_level) .AND. (idbg_mxmn < inDetail_level)) RETURN

    start_detail_timer(timer_dbg_prnt,10)    
    

    ! ! valid g-format without offset of decimal point
    ! 981 FORMAT(a,a12,':',a27,' C:',i3,  g26.18,3(a,i0,a,  g20.12))
    ! 982 FORMAT(a,a12,':',a27,'  :',i3,   26x,  3(a,i0,a,  g20.12))
    ! 991 FORMAT(a,a12,':',a27,'  :',i3, 2g26.18)
    ! ! valid e-format with first digit > zero
    ! 981 FORMAT(a,a12,':',a27,' C:',i3, 1pe26.18,3(a,i0,a,1pe20.12))
    ! 982 FORMAT(a,a12,':',a27,'  :',i3,    26x,  3(a,i0,a,1pe20.12))
    ! 991 FORMAT(a,a12,':',a27,'  :',i3,1p2e26.18)
    981 FORMAT(a,a12,':',a27,' C:',i3, 1pg26.18,3(a,i0,a,1pg20.12))
    982 FORMAT(a,a12,':',a27,'  :',i3,    26x,  3(a,i0,a,1pg20.12))
    992 FORMAT(a,a12,':',a27,'  :',i3, 1pg26.18, 1pg26.18, 1pg26.18)
    
    ! check print output level inDetail_level (1-5) with namelist given value (idbg_val)
    ! for output at given index
    
    !
    ! All calculations are done inside this IF only
    !
    
    IF (idbg_val >= inDetail_level) THEN
      !
      ! dimensions
      !                           !  index 1:nproma
      nlev    = SIZE(p_array,2)   !  vertical dimension (levels)
      ndimblk = SIZE(p_array,3)   !  blocks 1:nblks for cells/edges/verts
      
      ! output channel: stderr
      iout = nerr
      
      ! compare defined source string with namelist-given output string
      icheck_str_mod = 0
      DO jstr = 1, dim_mod_tst
        IF (place == str_mod_tst(jstr) .OR. str_mod_tst(jstr) == 'all') &
          & icheck_str_mod = 1
      END DO
      
      ! if place not found in str_mod_tst - no output
      IF (icheck_str_mod == 0) THEN
        stop_detail_timer(timer_dbg_prnt,10)
        RETURN
      ENDIF

      strout=TRIM(description)
      strmod=TRIM(place)
      
      ! check start and end index for output of vertical levels via namelist
      slev = 1
      IF (idbg_slev > 1)    slev = idbg_slev
      IF (slev      > nlev) slev = nlev
      elev = nlev
      IF (idbg_elev < nlev) elev = idbg_elev
      
      ! idbg_val<4: one level output only (slev), apart from permanent output (init)
      elev_val = elev
      IF (idbg_val < 4 .AND. inDetail_level > 0) elev_val = slev
      
      DO jk = slev, elev_val
        
        ! LL ERROR Note: this should be rewritten.
        ! it's not safe to use the number of blocks to identify the what is the grid entity
        ! write value at index
        IF (ndimblk == loc_nblks_c) THEN
          IF (p_pe == near_proc_id) &
            & WRITE(iout,981) '        VALUE ', strmod, strout, jk, p_array(c_i,jk,c_b), &
            & (' C',i,':',p_array(nc_i(i),jk,nc_b(i)),i=1,3)
        ELSE IF (ndimblk == loc_nblks_e) THEN
          IF (p_pe == near_proc_id) &
            & WRITE(iout,982) '        VALUE ', strmod, strout, jk, &
            & (' E',i,':',p_array(ne_i(i),jk,ne_b(i)),i=1,3)
        ELSE IF (ndimblk == loc_nblks_v) THEN
          IF (p_pe == near_proc_id) &
            & WRITE(iout,982) '        VALUE ', strmod, strout, jk, &
            & (' V',i,':',p_array(nv_i(i),jk,nv_b(i)),i=1,3)
        END IF
        
      END DO
      
    END IF
    
    ! check print output level inDetail_level (1-5) with namelist given value (idbg_mxmn)
    ! for MIN/MAX output:
    
    !
    ! All calculations are done inside this IF only
    !
    IF (idbg_mxmn >= inDetail_level ) THEN
      !
      ! dimensions
      !                           !  index 1:nproma
      nlev    = SIZE(p_array,2)   !  vertical dimension (levels)
      ndimblk = SIZE(p_array,3)   !  blocks 1:nblks for cells/edges/verts
      
      ! output channel: stderr
      iout = nerr
      
      ! compare defined source string with namelist-given output string
      icheck_str_mod = 0
      DO jstr = 1, dim_mod_tst
        IF (place == str_mod_tst(jstr) .OR. str_mod_tst(jstr) == 'all') &
          & icheck_str_mod = 1
      END DO
      
      ! if place not found in str_mod_tst - no output
      IF (icheck_str_mod == 0) THEN
        stop_detail_timer(timer_dbg_prnt,10)
        RETURN
      ENDIF
      
      strout=TRIM(description)
      strmod=TRIM(place)
      
      ! check start and end index for output of vertical levels via namelist
      slev = 1
      IF (idbg_slev > 1)    slev = idbg_slev
      IF (slev      > nlev) slev = nlev
      elev = nlev
      IF (idbg_elev < nlev) elev = idbg_elev
      
      ! idbg_mxmn<4: one level output only (slev), independent of elev_val
      elev_mxmn = elev
    !  IF (idbg_mxmn < 4 .AND. inDetail_level > 0) elev_mxmn = slev
      
      ! print out maximum and minimum value
      ! ctrn=minval(p_array(:, slev:elev_mxmn, :))
      DO jk = slev, elev_mxmn
        
        IF (PRESENT(in_subset)) THEN
          minmaxmean(:) = global_minmaxmean(values=p_array(:,:,:), in_subset=in_subset, start_level=jk, end_level=jk)
        ELSE
          minmaxmean(:) = global_minmaxmean(values=p_array(:,:,:), start_level=jk, end_level=jk)
        ENDIF
        
        IF (my_process_is_stdio()) &
          & WRITE(iout,992) ' MAX/MIN/MEAN ', strmod, strout, jk, minmaxmean(2), minmaxmean(1), minmaxmean(3)
        
        
        ! location of max/min - parallelize!
        ! WRITE(iout,983) ' LOC ',strout,jk, &
        !   &              MAXLOC(p_array(1:nproma,jk,1:ndimblk)),     &
        !   &              MINLOC(p_array(1:nproma,jk,1:ndimblk))
        ! 983 FORMAT(a,a12,':',a27,'  :',i3, 4i4)
        
      END DO
      
    END IF
    
    stop_detail_timer(timer_dbg_prnt,10)
    
  END SUBROUTINE dbg_print_3d
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! Print out min and max or a specific cell value and neighbors of a 2-dim array.
  !!
  SUBROUTINE dbg_print_2d( description, p_array, place, inDetail_level, in_subset )
    
    CHARACTER(LEN=*),      INTENT(in) :: description    ! description of array
    REAL(wp),              INTENT(in) :: p_array(:,:)   ! 2-dim array for debugging
    CHARACTER(LEN=*),      INTENT(in) :: place    ! defined string for source of current array
    INTEGER,               INTENT(in) :: inDetail_level    ! source level from module for print output
    TYPE(t_subset_range),  TARGET, OPTIONAL :: in_subset
    
    ! local variables
    CHARACTER(LEN=27) ::  strout
    CHARACTER(LEN=12) ::  strmod
    INTEGER ::  iout, icheck_str_mod, jstr, i, jk, ndimblk
    REAL(wp)          ::  minmaxmean(3)

    IF ((idbg_val < inDetail_level) .AND. (idbg_mxmn < inDetail_level)) RETURN
    
    start_detail_timer(timer_dbg_prnt,10)
    
    ! dimensions - first dimension is nproma
    ndimblk = SIZE(p_array,2)
    
    ! output channel: stderr
    iout = nerr
    
    ! compare defined source string with namelist-given output string
    icheck_str_mod = 0
    DO jstr = 1, dim_mod_tst
      IF (place == str_mod_tst(jstr) .OR. str_mod_tst(jstr) == 'all') &
        & icheck_str_mod = 1
    END DO
    
    ! if place not found in str_mod_tst - no output
    IF (icheck_str_mod == 0) THEN
      stop_detail_timer(timer_dbg_prnt,10)
      RETURN
    ENDIF
    

    ! ! valid e-format with first digit > zero
    ! 981 FORMAT(a,a12,':',a27,' C:',i3, 1pe26.18,3(a,i0,a,1pe20.12))
    ! 982 FORMAT(a,a12,':',a27,'  :',i3,    26x,  3(a,i0,a,1pe20.12))
    ! 991 FORMAT(a,a12,':',a27,'  :',i3,1p2e26.18)
    
    981 FORMAT(a,a12,':',a27,' C:',i3, 1pg26.18,3(a,i0,a,1pg20.12))
    982 FORMAT(a,a12,':',a27,'  :',i3,    26x,  3(a,i0,a,1pg20.12))
    992 FORMAT(a,a12,':',a27,'  :',i3, 1pg26.18, 1pg26.18, 1pg26.18)
    
    strout=TRIM(description)
    strmod=TRIM(place)
    
    ! surface level output only
    jk = 0
    
    
    ! check print output level inDetail_level (1-5) with namelist given value (idbg_val)
    ! for output at given index
    
    IF (idbg_val >= inDetail_level) THEN
      
      !write(iout,*) ' ndimblk and loc_nblks = ',ndimblk, loc_nblks_c, loc_nblks_e, loc_nblks_v
      
      ! write value at index
      IF (ndimblk == loc_nblks_c) THEN
        IF (p_pe == near_proc_id) &
          & WRITE(iout,981) '        VALUE ', strmod, strout, jk, p_array(c_i,c_b), &
          & (' C',i,':',p_array(nc_i(i),nc_b(i)),i=1,3)
      ELSE IF (ndimblk == loc_nblks_e) THEN
        IF (p_pe == near_proc_id) &
          & WRITE(iout,982) '        VALUE ', strmod, strout, jk, &
          & (' E',i,':',p_array(ne_i(i),ne_b(i)),i=1,3)
      ELSE IF (ndimblk == loc_nblks_v) THEN
        IF (p_pe == near_proc_id) &
          & WRITE(iout,982) '        VALUE ', strmod, strout, jk, &
          & (' V',i,':',p_array(nv_i(i),nv_b(i)),i=1,3)
      END IF
      
    END IF
    
    ! check print output level inDetail_level (1-5) with namelist given value (idbg_mxmn)
    ! for MIN/MAX output:
    
    IF (idbg_mxmn >= inDetail_level ) THEN

      IF (PRESENT(in_subset)) THEN
        minmaxmean(:) = global_minmaxmean(values=p_array(:,:), in_subset=in_subset)
      ELSE
        minmaxmean(:) = global_minmaxmean(values=p_array(:,:))
      ENDIF

      IF (my_process_is_stdio()) &
        & WRITE(iout,992) ' MAX/MIN/MEAN ', strmod, strout, jk, minmaxmean(2), minmaxmean(1), minmaxmean(3)
      
      ! WRITE(0,*) ' MAX/MIN/MEAN ', minmaxmean(2), minmaxmean(1), minmaxmean(3)
      
    END IF
    
    stop_detail_timer(timer_dbg_prnt,10)
    
  END SUBROUTINE dbg_print_2d
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Print out given  min, mean and max
  SUBROUTINE debug_print_MaxMinMean( description, minmaxmean, place, inDetail_level )

    CHARACTER(LEN=*),      INTENT(in) :: description    ! description of array
    REAL(wp),              INTENT(in) :: minmaxmean(3)   ! 2-dim array for debugging
    CHARACTER(LEN=*),      INTENT(in) :: place    ! defined string for source of current array
    INTEGER,               INTENT(in) :: inDetail_level    ! source level from module for print output

    992 FORMAT(a,a12,':',a27,'  :', 1pg26.18, 1pg26.18, 1pg26.18)

    IF (idbg_mxmn >= inDetail_level) THEN
      IF (my_process_is_stdio()) &
        & WRITE(nerr,992) ' MAX/MIN/MEAN ',TRIM(place), TRIM(description), minmaxmean(2), minmaxmean(1), minmaxmean(3)
    ENDIF

  END SUBROUTINE debug_print_MaxMinMean
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Print out given  min, mean and max
  SUBROUTINE debug_printValue( description, val, value1, value2, detail_level )

    CHARACTER(LEN=*),      INTENT(in) :: description    
    REAL(wp),              INTENT(in) :: val
    REAL(wp), OPTIONAL,    INTENT(in) :: value1, value2    
    INTEGER,               INTENT(in) :: detail_level    


!     992 FORMAT(a,a27,'  :', 1pg26.18)
!     993 FORMAT(a,a27,'  :', 1pg26.18, ", ", 1pg26.18, ", ", 1pg26.18)
!     994 FORMAT(a,a27,'  :', 1pg26.18, ", ", 1pg26.18, ", ", 1pg26.18)

    IF (idbg_mxmn < detail_level) RETURN
    IF (.NOT. my_process_is_stdio()) RETURN

    IF (PRESENT(value2)) THEN
       WRITE(nerr,*) TRIM(description), ": ", val, value1, value2
    ELSEIF (PRESENT(value1)) THEN
       WRITE(nerr,*) TRIM(description), ": ", val, value1
    ELSE
       WRITE(nerr,*) TRIM(description), ": ", val
    ENDIF

  END SUBROUTINE debug_printValue


END MODULE mo_util_dbg_prnt

