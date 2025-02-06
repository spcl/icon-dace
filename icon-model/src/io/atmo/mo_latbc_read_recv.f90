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

! This module provides basic methods to send
! and receive file data using asynchronous
! communication.

MODULE mo_latbc_read_recv

#ifndef NOMPI
  USE mpi
#endif

  USE mo_kind,               ONLY: sp, i8
  USE mo_exception,          ONLY: finish, message, message_text
  USE mo_impl_constants,     ONLY: SUCCESS
  USE mo_cdi_constants,      ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE
  USE mo_mpi,                ONLY: p_pe_work,  &
    &                              num_work_procs, p_real_sp               
  USE mo_util_cdi,           ONLY: get_cdi_varID
  USE mo_async_latbc_types,  ONLY: t_patch_data, t_latbc_data
  USE mo_reorder_info,       only: t_reorder_info
  USE mo_cdi,                ONLY: streamInqVlist, vlistInqVarZaxis, vlistInqVarGrid, gridInqSize, zaxisInqSize, &
                                 & streamReadVarSliceF
  USE mo_limarea_config,     ONLY: latbc_config
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: prefetch_cdi_2d 
  PUBLIC  :: prefetch_cdi_3d
  PUBLIC  :: compute_data_receive

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_latbc_read_recv::'

CONTAINS

  !-------------------------------------------------------------------------
  !> Read 3D dataset from file.
  ! 
  !  Note: This implementation uses a 2D buffer.
  ! 
  SUBROUTINE prefetch_cdi_3d(streamID, varname, latbc_data, nlevs, hgrid, ioff)
#ifndef NOMPI
    INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(INOUT) :: ioff(0:)
#else
    INTEGER, INTENT(INOUT)                        :: ioff(0:)
#endif
    INTEGER,                 INTENT(IN)  :: streamID     !< ID of CDI file stream
    CHARACTER(len=*),        INTENT(IN)  :: varname      !< Var name of field to be read
    !> patch data containing information for prefetch
    TYPE(t_latbc_data), TARGET, INTENT(IN) :: latbc_data
    INTEGER,                 INTENT(OUT) :: nlevs        !< return value: no. of vertical levels
    INTEGER,                 INTENT(IN)  :: hgrid        !< stored variable location indication
    
    ! local constants:
    CHARACTER(len=*), PARAMETER :: routine = modname//'prefetch_cdi_3d'
    ! local variables:
    INTEGER                         :: vlistID, varID, zaxisID, gridID,   &
      &                                jk, ierrstat, dimlen(2), nmiss
    REAL(sp), ALLOCATABLE :: read_buf(:) ! temporary local array for reading

    INTEGER                         :: nread
    TYPE(t_reorder_info), POINTER :: p_ri

#ifndef NOMPI
    ! allocate a buffer for one vertical level
    IF (hgrid == GRID_UNSTRUCTURED_CELL) THEN

      IF (latbc_config%lsparse_latbc) THEN
        nread = SIZE(latbc_data%global_index%cells)
      ELSE
        nread = latbc_data%patch_data%n_patch_cells_g
      END IF
      p_ri => latbc_data%patch_data%cells

    ELSE IF (hgrid == GRID_UNSTRUCTURED_EDGE) THEN

      IF (latbc_config%lsparse_latbc) THEN
        nread = SIZE(latbc_data%global_index%edges)
      ELSE
        nread = latbc_data%patch_data%n_patch_edges_g
      END IF
      p_ri => latbc_data%patch_data%edges

    ELSE
      CALL finish(routine, "invalid grid type")
    ENDIF

    ! allocate read buffer (which might be smaller than domain, if
    ! only boudary rows are read from file.
    ALLOCATE(read_buf(nread), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! get var ID
    vlistID   = streamInqVlist(streamID)
    varID     = get_cdi_varID(streamID, name=TRIM(varname)) 
    zaxisID   = vlistInqVarZaxis(vlistID, varID)
    gridID    = vlistInqVarGrid(vlistID, varID)
    dimlen(1) = gridInqSize(gridID)
    dimlen(2) = zaxisInqSize(zaxisID)
    nlevs     = dimlen(2) ! vertical levels of netcdf file

    ! Check variable dimensions:
    IF (dimlen(1) /= SIZE(read_buf)) THEN
       WRITE(message_text,'(a,2i4,a,i0)') "Horizontal cells: ", dimlen(1), SIZE(read_buf), &
         &                           "nlev: ", nlevs
       CALL message(routine, message_text)
       CALL finish(routine, "Incompatible dimensions!")
    END IF

    ! fixme: send more than 1 level at a time
    DO jk=1, nlevs

      ! read record as 1D field
      CALL streamReadVarSliceF(streamID, varID, jk-1, read_buf, nmiss)

      ! send 2d buffer using MPI_PUT
      CALL prefetch_proc_send(latbc_data%patch_data%mem_win%mpi_win, read_buf, 1, p_ri, ioff)
    ENDDO ! jk=1,nlevs 
  
    DEALLOCATE(read_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")
#endif
  
  END SUBROUTINE prefetch_cdi_3d

  !-------------------------------------------------------------------------
  !> Read 2D dataset from file, implementation for REAL fields
  !
  SUBROUTINE prefetch_cdi_2d (streamID, varname, latbc_data, hgrid, ioff)
#ifndef NOMPI
    INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(INOUT) :: ioff(0:)
#else
    INTEGER, INTENT(INOUT)                        :: ioff(0:)
#endif
    CHARACTER(len=*), PARAMETER :: routine = modname//'prefetch_cdi_2d'
    INTEGER,                 INTENT(IN) :: streamID     !< ID of CDI file stream
    CHARACTER(len=*),        INTENT(IN) :: varname      !< Varname of field to be read
    TYPE(t_latbc_data),      INTENT(IN) :: latbc_data   !< patch data containing information for prefetch 
    INTEGER,                 INTENT(IN) :: hgrid        !< stored variable location indication
    INTEGER :: nlevs_read

    CALL prefetch_cdi_3d(streamID, varname, latbc_data, nlevs_read, hgrid, ioff)
    IF (nlevs_read /= 1) CALL finish(routine, "Invalid number of vertical levels for "//TRIM(varname)//"!")
  END SUBROUTINE prefetch_cdi_2d

  !-------------------------------------------------------------------------
  !> compute processor copy 2D or 3D dataset from memory window buffer.
  ! 
  SUBROUTINE compute_data_receive (hgrid, nlevs, var_out, eoff, patch_data)
 
    INTEGER,             INTENT(IN)    :: hgrid          !< stored variable location indication
    INTEGER,             INTENT(IN)    :: nlevs          !< vertical levels of netcdf file
    REAL(sp),            INTENT(INOUT) :: var_out(:,:,:) !< output field
    INTEGER(i8),         INTENT(INOUT) :: eoff
    TYPE(t_patch_data),  TARGET, INTENT(IN)   :: patch_data
    TYPE(t_reorder_info), POINTER             :: p_ri
    ! local constants:
    CHARACTER(len=*), PARAMETER :: &
         routine = modname//'compute_data_receive'
    ! local variables:
    INTEGER     :: j, jl, jb, jk, mpi_error      

#ifndef NOMPI
    ! Get pointer to appropriate reorder_info
    SELECT CASE (hgrid)
    CASE(GRID_UNSTRUCTURED_CELL)
       p_ri => patch_data%cells
    CASE(GRID_UNSTRUCTURED_EDGE)
       p_ri => patch_data%edges
    CASE default
       CALL finish(routine,'unknown grid type')
    END SELECT

    CALL MPI_Win_lock(MPI_LOCK_SHARED, p_pe_work, MPI_MODE_NOCHECK, patch_data%mem_win%mpi_win, mpi_error)

    ! initialize output field:
    var_out(:,:,:) = 0._sp

    DO jk=1, nlevs
      DO j = 1,  p_ri%n_own 
        jb = p_ri%own_blk(j) ! Block index in distributed patch
        jl = p_ri%own_idx(j) ! Line  index in distributed patch
        var_out(jl,jk,jb) = REAL(patch_data%mem_win%mem_ptr_sp(eoff+INT(j,i8)),sp)
      ENDDO
      eoff = eoff + INT(p_ri%n_own,i8) 
    END DO ! jk=1,nlevs

    CALL MPI_Win_unlock(p_pe_work, patch_data%mem_win%mpi_win, mpi_error)
#endif

  END SUBROUTINE compute_data_receive
 
  !-------------------------------------------------------------------------
  !> routine on the prefetch PE to write variable values to memory window
  !
  !  @note This subroutine is called by prefetch PE only.
  !  Initial revision by M. Pondkule, DWD (2014-05-27) 
  !
  ! fixme: allow actually using more than nlevs = 1
  SUBROUTINE prefetch_proc_send(win, var1_sp, nlevs, p_ri, ioff)
#ifndef NOMPI
    INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(INOUT) :: ioff(0:)
#else
    INTEGER, INTENT(INOUT)                        :: ioff(0:)
#endif
    ! local variables  
    INTEGER, INTENT(IN) :: win
    REAL(sp),                   INTENT(IN) :: var1_sp(:)
    INTEGER,                    INTENT(IN) :: nlevs
    TYPE(t_reorder_info), INTENT(in) :: p_ri

    CHARACTER(len=*), PARAMETER :: routine = modname//'prefetch_proc_send'
    INTEGER                        :: voff(0:num_work_procs-1) 
    INTEGER                        :: nv_off_np(0:num_work_procs)
    REAL(sp), ALLOCATABLE          :: var3_sp(:)
    INTEGER                        :: np, nval, nv_off, mpi_error, ierrstat
    INTEGER                        :: dst_start, dst_end, src_start, src_end
#ifndef NOMPI
#ifdef DO_NOT_COMBINE_PUT_AND_NOCHECK
    INTEGER, PARAMETER :: lock_assert = 0
#else
    INTEGER, PARAMETER :: lock_assert = MPI_MODE_NOCHECK
#endif
#endif


    IF (.NOT. ALLOCATED(p_ri%pe_own)) THEN
      CALL finish(routine, "Internal error: data structure p_ri%pe_own unallocated!")
    END IF

    ! compute the total offset for each PE
    nv_off       = 0
    nv_off_np(0) = 0
    DO np = 0, num_work_procs-1
       voff(np)        = nv_off
       nval            = p_ri%pe_own(np) * nlevs
       nv_off          = nv_off + nval
       nv_off_np(np+1) = nv_off_np(np) + p_ri%pe_own(np)
    END DO

    ALLOCATE(var3_sp(p_ri%n_glb), STAT=ierrstat) ! Must be allocated to exact size
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      
!$OMP PARALLEL
!$OMP DO PRIVATE(dst_start, dst_end, src_start, src_end)
          DO np = 0, num_work_procs-1
            dst_start = nv_off_np(np)+1
            dst_end   = nv_off_np(np+1)
            src_start = voff(np)+1
            src_end   = voff(np)+p_ri%pe_own(np)
            IF ((src_end-src_start) /= (dst_end-dst_start)) THEN
              WRITE (0,*) "(src_end-src_start+1) = ", (src_end-src_start+1)
              WRITE (0,*) "(dst_end-dst_start+1) = ", (dst_end-dst_start+1)
              CALL finish(routine, "internal error!")
            END IF
            var3_sp(dst_start:dst_end) = var1_sp(p_ri%reorder_index(src_start:src_end))
          ENDDO
!$OMP END DO
!$OMP END PARALLEL

    nv_off  = 0_i8
    nval = 0

#ifndef NOMPI
    DO np = 0, num_work_procs-1

       IF (p_ri%pe_own(np) == 0) CYCLE

       nval = p_ri%pe_own(np)*nlevs

       ! fixme: use mpi_win_lock_all here
       CALL MPI_Win_lock(MPI_LOCK_EXCLUSIVE, np, lock_assert, win, mpi_error)
         
       ! consistency check:
       IF (SIZE(var3_sp) < nv_off+nval) THEN
         WRITE (0,*) "SIZE(var3_sp) = ", SIZE(var3_sp), " < nv_off+nval = ", nv_off+nval
         CALL finish(routine, "Internal error!")
       END IF

       CALL MPI_PUT(var3_sp(nv_off+1), nval, p_real_sp, np, ioff(np), nval, p_real_sp, &
         & win, mpi_error) !MPI_WIN_NULL) 
         
       CALL MPI_Win_unlock(np, win, mpi_error)

       ! Update the offset in var1
       nv_off = nv_off + nval

       ! Update the offset in the memory window on compute PEs
       ioff(np) = ioff(np) + INT(nval,i8) ! replace with mpi_kind

    ENDDO
#endif

    DEALLOCATE(var3_sp, STAT=ierrstat) ! Must be allocated to exact size
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

  END SUBROUTINE prefetch_proc_send

  !------------------------------------------------------------------------------------------------

END MODULE mo_latbc_read_recv
