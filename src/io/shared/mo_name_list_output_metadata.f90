! Module handling the transfer of variable meta info to asynchronous I/O PEs.
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

MODULE mo_name_list_output_metadata

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_f_pointer
  USE mo_exception,                         ONLY: finish
  USE mo_var_metadata_types,                ONLY: t_var_metadata, var_metadata_get_size, &
    & var_metadata_toBinary, var_metadata_fromBinary
  USE mo_name_list_output_types,            ONLY: t_mem_win
  USE mo_mpi,                               ONLY: p_int, p_comm_work_io,      &
    &                                             p_comm_rank
  USE mo_dynamics_config,                   ONLY: nnow, nnow_rcf, nnew, nnew_rcf
  USE mo_impl_constants,                    ONLY: TLEV_NNOW, TLEV_NNOW_RCF, TLEV_NNEW, TLEV_NNEW_RCF

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: metainfo_allocate_memory_window
  PUBLIC :: metainfo_write_to_memwin
  PUBLIC :: metainfo_get_from_buffer
  PUBLIC :: metainfo_get_timelevel

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_name_list_output_metadata'

CONTAINS

  !-------------------------------------------------------------------------------------------------
  !> Allocates an MPI memory window for the meta info of the variables fields.
  !   - allocation for asynchronous I/O mode only.
  !   - allocation on I/O PEs and PE #0 only.
  !
  SUBROUTINE metainfo_allocate_memory_window(memwin, nvars)

    TYPE(t_mem_win),      INTENT(INOUT) :: memwin ! MPI memory window
    INTEGER,              INTENT(IN)    :: nvars  ! total no. of variables
  END SUBROUTINE metainfo_allocate_memory_window


  !-------------------------------------------------------------------------------------------------
  !> Store a variable's meta-info to a memory window.
  !  This subroutine does nothing on PEs except compute PE #0 or
  !  if we are not running in asynchronous I/O mode.
  !
  SUBROUTINE metainfo_write_to_memwin(memwin, ivar, info)
    TYPE(t_mem_win),      INTENT(INOUT) :: memwin ! MPI memory window
    INTEGER,              INTENT(IN)    :: ivar   ! index of variable (corresponds to data memwin)
    TYPE(t_var_metadata), INTENT(IN)    :: info   ! meta data for variable
    CHARACTER(*), PARAMETER :: routine = modname//"::metainfo_write_to_memwin"

    IF (.NOT. ASSOCIATED(memwin%mem_ptr_metainfo_pe0)) &
      & CALL finish(routine, "Internal error!")

    ! copy the info object into the memory window
    memwin%mem_ptr_metainfo_pe0(:, ivar) = var_metadata_toBinary(info, &
      & SIZE(memwin%mem_ptr_metainfo_pe0, 1))
  END SUBROUTINE metainfo_write_to_memwin


  !-------------------------------------------------------------------------------------------------
  !> Retrieve a variable's meta-info from memory window.
  !  This subroutine does nothing on compute PEs or
  !  if we are not running in asynchronous I/O mode.
  !
  SUBROUTINE metainfo_get_from_buffer(buf, info)
    INTEGER, INTENT(IN) :: buf(:)
    TYPE(t_var_metadata), INTENT(OUT) :: info   ! meta data for variable

    ! copy the info object from the memory window
    info = var_metadata_fromBinary(buf, SIZE(buf))
  END SUBROUTINE metainfo_get_from_buffer

  !-------------------------------------------------------------------------------------------------
  !> Return the output timelevel for given variables info object
  !
  INTEGER FUNCTION metainfo_get_timelevel(info,domain) RESULT(timelevel)
    TYPE(t_var_metadata), INTENT(IN) :: info
    INTEGER, INTENT(IN)              :: domain
    CHARACTER(*), PARAMETER       :: routine = modname//"::metainfo_get_timelevel"

    SELECT CASE (info%tlev_source)
    CASE(TLEV_NNOW);     timelevel = nnow(domain)
    CASE(TLEV_NNOW_RCF); timelevel = nnow_rcf(domain)
    CASE(TLEV_NNEW);     timelevel = nnew(domain)
    CASE(TLEV_NNEW_RCF); timelevel = nnew_rcf(domain)
    CASE DEFAULT
      CALL finish(routine,'Unsupported tlev_source')
    END SELECT
  END FUNCTION metainfo_get_timelevel
END MODULE mo_name_list_output_metadata
