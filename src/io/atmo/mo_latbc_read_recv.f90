! This module provides basic methods to send
! and receive file data using asynchronous
! communication.
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

MODULE mo_latbc_read_recv





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



    INTEGER, INTENT(INOUT)                        :: ioff(0:)

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

  
  END SUBROUTINE prefetch_cdi_3d

  !-------------------------------------------------------------------------
  !> Read 2D dataset from file, implementation for REAL fields
  !
  SUBROUTINE prefetch_cdi_2d (streamID, varname, latbc_data, hgrid, ioff)
    INTEGER, INTENT(INOUT)                        :: ioff(0:)
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


  END SUBROUTINE compute_data_receive
 
  !-------------------------------------------------------------------------
  !> routine on the prefetch PE to write variable values to memory window
  !
  !  @note This subroutine is called by prefetch PE only.
  !  Initial revision by M. Pondkule, DWD (2014-05-27) 
  !
  ! fixme: allow actually using more than nlevs = 1
  SUBROUTINE prefetch_proc_send(win, var1_sp, nlevs, p_ri, ioff)
    INTEGER, INTENT(INOUT)                        :: ioff(0:)
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


    DEALLOCATE(var3_sp, STAT=ierrstat) ! Must be allocated to exact size
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

  END SUBROUTINE prefetch_proc_send

  !------------------------------------------------------------------------------------------------

END MODULE mo_latbc_read_recv
