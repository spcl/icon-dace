! This module provides basic methods for reading
! a NetCDF file in parallel or sequential in a transparent way.
!
! Contains routines for reading data from netcdf-Files of various shape.
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
#define define_fill_target REAL(wp), TARGET, OPTIONAL
!#define define_fill_target REAL(wp), TARGET, ALLOCATABLE, OPTIONAL
!#define define_fill_target REAL(wp), POINTER, OPTIONAL
#define define_fill_target_int INTEGER, TARGET, OPTIONAL
!#define define_fill_target_int INTEGER, TARGET, ALLOCATABLE, OPTIONAL
!#define define_fill_target_int INTEGER, POINTER, OPTIONAL

MODULE mo_read_netcdf_broadcast_2

  USE mo_kind,               ONLY: sp, dp, wp
  USE mo_exception,          ONLY: message, warning, finish
  USE mo_impl_constants,     ONLY: success
  USE mo_parallel_config,    ONLY: nproma
  USE mo_io_units,           ONLY: filename_max
  USE mo_mpi,                ONLY: my_process_is_mpi_workroot, p_comm_work, &
    &                              process_mpi_root_id, p_bcast
  USE mo_communication,      ONLY: t_scatterPattern
  USE mo_fortran_tools,      ONLY: t_ptr_2d, t_ptr_2d_int, t_ptr_3d, t_ptr_3d_int
  USE mo_netcdf_errhandler,  ONLY: nf
  USE mo_netcdf
  !-------------------------------------------------------------------------

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: netcdf_open_input, netcdf_close

  PUBLIC :: netcdf_read_att_int
  PUBLIC :: netcdf_read_inq_varexists
  PUBLIC :: netcdf_read_0D_real
  PUBLIC :: netcdf_read_0D_int
  PUBLIC :: netcdf_read_1D
  PUBLIC :: netcdf_read_1D_extdim_time
  PUBLIC :: netcdf_read_1D_extdim_extdim_time
  PUBLIC :: netcdf_read_extdim_slice_extdim_extdim_extdim
  PUBLIC :: netcdf_read_2D_int
  PUBLIC :: netcdf_read_2D
  PUBLIC :: netcdf_read_REAL_2D_all
  PUBLIC :: netcdf_read_2D_time
  PUBLIC :: netcdf_read_REAL_3D_all
  PUBLIC :: netcdf_read_3D
  PUBLIC :: netcdf_read_3D_time
  PUBLIC :: netcdf_read_2D_extdim
  PUBLIC :: netcdf_read_2D_extdim_int
  PUBLIC :: netcdf_read_3D_extdim
  PUBLIC :: netcdf_get_missValue
  PUBLIC :: t_p_scatterPattern

  TYPE t_p_scatterPattern
    CLASS(t_scatterPattern), POINTER :: p
  END TYPE t_p_scatterPattern

  INTERFACE netcdf_read_att_int
    MODULE PROCEDURE netcdf_read_ATT_INT
  END INTERFACE netcdf_read_att_int

  INTERFACE netcdf_read_0D_real
    MODULE PROCEDURE netcdf_read_REAL_0D
  END INTERFACE netcdf_read_0D_real

  INTERFACE netcdf_read_0D_int
    MODULE PROCEDURE netcdf_read_INT_0D
  END INTERFACE netcdf_read_0D_int

  INTERFACE netcdf_read_1D
    MODULE PROCEDURE netcdf_read_REAL_1D
  END INTERFACE netcdf_read_1D

  INTERFACE netcdf_read_1D_extdim_time
    MODULE PROCEDURE netcdf_read_REAL_1D_extdim_time
  END INTERFACE netcdf_read_1D_extdim_time

  INTERFACE netcdf_read_1D_extdim_extdim_time
    MODULE PROCEDURE netcdf_read_REAL_1D_extdim_extdim_time
  END INTERFACE netcdf_read_1D_extdim_extdim_time

  INTERFACE netcdf_read_extdim_slice_extdim_extdim_extdim
    MODULE PROCEDURE netcdf_read_REAL_extdim_slice_extdim_extdim_extdim
  END INTERFACE netcdf_read_extdim_slice_extdim_extdim_extdim

  INTERFACE netcdf_read_2D_int
    MODULE PROCEDURE netcdf_read_INT_2D
    MODULE PROCEDURE netcdf_read_INT_2D_multivar
  END INTERFACE netcdf_read_2D_int

  INTERFACE netcdf_read_2D
    MODULE PROCEDURE netcdf_read_REAL_2D
    MODULE PROCEDURE netcdf_read_REAL_2D_multivar
  END INTERFACE netcdf_read_2D

  INTERFACE netcdf_read_2D_time
    MODULE PROCEDURE netcdf_read_REAL_2D_time
  END INTERFACE netcdf_read_2D_time

  INTERFACE netcdf_read_2D_extdim
    MODULE PROCEDURE netcdf_read_REAL_2D_extdim
    MODULE PROCEDURE netcdf_read_REAL_2D_extdim_multivar
  END INTERFACE netcdf_read_2D_extdim

  INTERFACE netcdf_read_2D_extdim_int
    MODULE PROCEDURE netcdf_read_INT_2D_extdim
    MODULE PROCEDURE netcdf_read_INT_2D_extdim_multivar
  END INTERFACE netcdf_read_2D_extdim_int

  INTERFACE netcdf_read_3D
    MODULE PROCEDURE netcdf_read_REAL_3D
  END INTERFACE netcdf_read_3D

  INTERFACE netcdf_read_3D_time
    MODULE PROCEDURE netcdf_read_REAL_3D_time
  END INTERFACE netcdf_read_3D_time

  INTERFACE netcdf_read_3D_extdim
    MODULE PROCEDURE netcdf_read_REAL_3D_extdim
  END INTERFACE netcdf_read_3D_extdim

  INTEGER, PARAMETER :: MAX_VAR_DIMS = 16 ! NF90_MAX_VAR_DIMS

CONTAINS

  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_ATT_INT(file_id, variable_name, attribute_name) result(res)

    INTEGER                      :: res

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    CHARACTER(LEN=*), INTENT(IN) :: attribute_name

    INTEGER :: varid, var_type, var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)
    INTEGER :: zlocal(1)

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_ATT_INT'


    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, &
        &                 var_size, var_dim_name)

      CALL nf(nf90_get_att(file_id, varid, attribute_name, zlocal(:)), &
              attribute_name)
    ENDIF

    ! broadcast...
    CALL p_bcast(zlocal, process_mpi_root_id, p_comm_work)

    res=zlocal(1)

  END FUNCTION netcdf_read_ATT_INT
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_inq_varexists(file_id, variable_name) result(ret)
    LOGICAL                      :: ret
    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name

    INTEGER                      :: err, varid
    CHARACTER(LEN=*), PARAMETER  :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_inq_varexists'

    IF( my_process_is_mpi_workroot()  ) THEN
      err = nf90_inq_varid(file_id, variable_name, varid)
    ENDIF

    CALL p_bcast(err, process_mpi_root_id, p_comm_work)

    ret = (err == nf90_noerr)
  END FUNCTION netcdf_read_inq_varexists
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_0D(file_id, variable_name) result(res)

    REAL(wp)            :: res

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name

    INTEGER :: varid, var_type, var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)
    REAL(wp)     :: zlocal(1)

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_REAL_0D'


    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, &
        &                 var_size, var_dim_name)

      CALL nf(nf90_get_var(file_id, varid, zlocal(:)), variable_name)
    ENDIF

    ! broadcast...
    CALL p_bcast(zlocal, process_mpi_root_id, p_comm_work)

    res=zlocal(1)

  END FUNCTION netcdf_read_REAL_0D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_INT_0D(file_id, variable_name) result(res)

    INTEGER                      :: res

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name

    INTEGER :: varid, var_type, var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)
    INTEGER :: zlocal(1)

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_INT_0D'


    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, &
        &                 var_size, var_dim_name)

      CALL nf(nf90_get_var(file_id, varid, zlocal(:)), variable_name)
    ENDIF

    ! broadcast...
    CALL p_bcast(zlocal, process_mpi_root_id, p_comm_work)

    res=zlocal(1)

  END FUNCTION netcdf_read_INT_0D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_1D(file_id, variable_name, fill_array) &
    result(res)

    REAL(wp), POINTER            :: res(:)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:)

    INTEGER :: varid, var_type, var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)
    INTEGER :: return_status

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_REAL_1D'

    ! trivial return value.
    NULLIFY(res)

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, &
        &                 var_size, var_dim_name)

      ! check if the dims look ok
      IF (var_dims /= 1 ) THEN
        write(0,*) "var_dims = ", var_dims, " var_size=", var_size(1)
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

    ENDIF

    ! we need to sync the var_size...
    CALL p_bcast(var_size(1:2), process_mpi_root_id, p_comm_work)

    IF (PRESENT(fill_array)) THEN
      res => fill_array
    ELSE
      ALLOCATE( res(var_size(1)), stat=return_status )
      IF (return_status /= success) THEN
        CALL finish (method_name, 'ALLOCATE( netcdf_read_REAL_1D )')
      ENDIF
    ENDIF

    ! check if the size is correct
    IF (SIZE(res,1) < var_size(1)) &
      CALL finish(method_name, "allocated size < var_size")
    IF (SIZE(res,1) > var_size(1)) &
      CALL warning(method_name, "allocated size > var_size")

    IF( my_process_is_mpi_workroot()) THEN
      CALL nf(nf90_get_var(file_id, varid, res(:)), variable_name)
    ENDIF

    ! broadcast...
    CALL p_bcast(res, process_mpi_root_id, p_comm_work)

  END FUNCTION netcdf_read_REAL_1D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_1D_extdim_time(file_id, variable_name, &
    &                                      fill_array, dim_names, &
    &                                      start_timestep, end_timestep) &
    result(res)

    REAL(wp), POINTER            :: res(:,:,:)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_names(:)
    INTEGER, INTENT(IN), OPTIONAL:: start_timestep, end_timestep

    INTEGER :: varid, var_type, var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)
    INTEGER :: file_time_steps, time_steps, start_time, end_time
    INTEGER :: start_read_index(3), count_read_index(3)
    INTEGER :: idim
    INTEGER :: return_status

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_REAL_1D_extdim_time'

    ! trivial return value.
    NULLIFY(res)

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, &
        &                 var_size, var_dim_name)

      ! check if the dims look ok
      IF (var_dims /= 3 ) THEN
        WRITE(0,*) "var_dims = ", var_dims
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

      IF (PRESENT(dim_names)) THEN
        DO idim = 1, 2
          IF (dim_names(idim) /= var_dim_name(idim)) THEN
            WRITE(0,*) 'dim_name(',idim,')=',TRIM(ADJUSTL(var_dim_name(idim))),&
                       ' but dimension name', &
                       TRIM(ADJUSTL(dim_names(idim))),' was expected'
            CALL finish(method_name, 'dimension name mismatch')
          END IF
        END DO
      END IF

      IF (var_dim_name(3) /= 'time' ) THEN
        WRITE(0,*) 'no time dimension found, the last dimension must be time'
        CALL finish(method_name, 'no time dimension found')
      END IF
    ENDIF

    ! we need to sync the var_size...
    CALL p_bcast(var_size(1:3), process_mpi_root_id, p_comm_work)
    file_time_steps=var_size(3)

    ! calculate time range
    IF (PRESENT(start_timestep)) THEN
      start_time = start_timestep
    ELSE
      start_time = 1
    ENDIF
    IF (PRESENT(end_timestep)) THEN
      end_time = end_timestep
    ELSE
      end_time = file_time_steps
    ENDIF
!!$    use_time_range = (start_time /= 1) .OR. (end_time /= file_time_steps)
    time_steps = end_time - start_time + 1
!    write(0,*) "start,end time=", start_time, end_time
    IF (time_steps < 1) &
      & CALL finish(method_name, "number of time steps < 1")

    IF (PRESENT(fill_array)) THEN
      res => fill_array
    ELSE
      ALLOCATE( res(var_size(1),var_size(2),time_steps), stat=return_status )
      IF (return_status /= success) THEN
        CALL finish (method_name, 'ALLOCATE( res )')
      ENDIF
      res(:,:,:)=0.0_wp
    ENDIF

!!$    ! check if the size is correct
!!$    IF (SIZE(netcdf_read_REAL_1D,1) < var_size(1)) &
!!$      CALL finish(method_name, "allocated size < var_size")
!!$    IF (SIZE(netcdf_read_REAL_1D,1) > var_size(1)) &
!!$      CALL warning(method_name, "allocated size > var_size")

    IF( my_process_is_mpi_workroot()) THEN
      start_read_index = (/1,1,start_time/)
      count_read_index = (/var_size(1),var_size(2),time_steps/)
      CALL nf(nf90_get_var(file_id, varid, res(:,:,:), &
        &                  start_read_index, count_read_index), &
        &     variable_name)
    ENDIF

    ! broadcast...
    CALL p_bcast(res, process_mpi_root_id, p_comm_work)

  END FUNCTION netcdf_read_REAL_1D_extdim_time
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_1D_extdim_extdim_time( &
    file_id, variable_name, fill_array, dim_names, start_timestep, &
    end_timestep) result(res)

    REAL(wp), POINTER            :: res(:,:,:,:)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_names(:)
    INTEGER, INTENT(IN), OPTIONAL:: start_timestep, end_timestep

    INTEGER :: varid, var_type, var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)
    INTEGER :: file_time_steps, time_steps, start_time, end_time
    INTEGER :: start_read_index(4), count_read_index(4)
    INTEGER :: idim
    INTEGER :: return_status

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_REAL_1D_extdim_extdim_time'

    ! trivial return value.
    NULLIFY(res)

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, &
        &                 var_size, var_dim_name)

      ! check if the dims look ok
      IF (var_dims /= 4 ) THEN
        WRITE(0,*) "var_dims = ", var_dims
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

      IF (PRESENT(dim_names)) THEN
        DO idim = 1, 3
          IF (dim_names(idim) /= var_dim_name(idim)) THEN
            WRITE(0,*) 'dim_name(',idim,')=',TRIM(ADJUSTL(var_dim_name(idim))),&
                       ' but dimension name ', &
                       TRIM(ADJUSTL(dim_names(idim))),' was expected'
            CALL finish(method_name, 'dimension name mismatch')
          END IF
        END DO
      END IF

      IF (var_dim_name(4) /= 'time' ) THEN
        WRITE(0,*) 'no time dimension found, the last dimension must be time'
        CALL finish(method_name, 'no time dimension found')
      END IF
    ENDIF

    ! we need to sync the var_size...
    CALL p_bcast(var_size(1:4), process_mpi_root_id, p_comm_work)
    file_time_steps=var_size(4)

    ! calculate time range
    IF (PRESENT(start_timestep)) THEN
      start_time = start_timestep
    ELSE
      start_time = 1
    ENDIF
    IF (PRESENT(end_timestep)) THEN
      end_time = end_timestep
    ELSE
      end_time = file_time_steps
    ENDIF
!!$    use_time_range = (start_time /= 1) .OR. (end_time /= file_time_steps)
    time_steps = end_time - start_time + 1
!    write(0,*) "start,end time=", start_time, end_time
    IF (time_steps < 1) &
      & CALL finish(method_name, "number of time steps < 1")

    IF (PRESENT(fill_array)) THEN
      res => fill_array
    ELSE
      ALLOCATE( res(var_size(1),var_size(2), &
                var_size(3),time_steps), stat=return_status )
      IF (return_status /= success) &
        CALL finish (method_name, &
          &          'ALLOCATE( netcdf_read_REAL_1D_extdim_extdim_time )')
      res(:,:,:,:)=0.0_wp
    ENDIF

!!$    ! check if the size is correct
!!$    IF (SIZE(netcdf_read_REAL_1D,1) < var_size(1)) &
!!$      CALL finish(method_name, "allocated size < var_size")
!!$    IF (SIZE(netcdf_read_REAL_1D,1) > var_size(1)) &
!!$      CALL warning(method_name, "allocated size > var_size")

    IF( my_process_is_mpi_workroot()) THEN
      start_read_index = (/1,1,1,start_time/)
      count_read_index = (/var_size(1),var_size(2),var_size(3),time_steps/)
      CALL nf(nf90_get_var(file_id, varid, res(:,:,:,:), &
        &                  start_read_index, count_read_index), &
        &     variable_name)
    ENDIF

    ! broadcast...
    CALL p_bcast(res, process_mpi_root_id, p_comm_work)

  END FUNCTION netcdf_read_REAL_1D_extdim_extdim_time
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_extdim_slice_extdim_extdim_extdim( &
    file_id, variable_name, fill_array, dim_names, start_extdim1, &
    end_extdim1) result(res)

    REAL(wp), POINTER            :: res(:,:,:,:)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_names(:)
    INTEGER, INTENT(IN), OPTIONAL:: start_extdim1, end_extdim1

    INTEGER :: varid, var_type, var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)
    INTEGER :: file_extdim1_steps, dim1_steps, start_dim1, end_dim1
    INTEGER :: start_read_index(4), count_read_index(4)
    INTEGER :: idim
    INTEGER :: return_status

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_REAL_extdim_slice_extdim_extdim_extdim'

    ! trivial return value.
    NULLIFY(res)

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, &
        &                 var_size, var_dim_name)

      ! check if the dims look ok
      IF (var_dims /= 4 ) THEN
        WRITE(0,*) "var_dims = ", var_dims
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

      IF (PRESENT(dim_names)) THEN
        DO idim = 1, 4
          IF (TRIM(dim_names(idim)) /= TRIM(var_dim_name(idim))) THEN
            WRITE(0,*) 'dim_name(',idim,')=',TRIM(ADJUSTL(var_dim_name(idim))),&
                       ' but dimension name ', &
                       TRIM(ADJUSTL(dim_names(idim))),' was expected'
            CALL finish(method_name, 'dimension name mismatch')
          END IF
        END DO
      END IF
    END IF

    ! we need to sync the var_size...
    CALL p_bcast(var_size(1:4), process_mpi_root_id, p_comm_work)
    file_extdim1_steps=var_size(1)

    ! calculate slice of 4th dimension
    IF (PRESENT(start_extdim1)) THEN
      start_dim1 = start_extdim1
    ELSE
      start_dim1 = 1
    ENDIF
    IF (PRESENT(end_extdim1)) THEN
      end_dim1 = end_extdim1
    ELSE
      end_dim1 = file_extdim1_steps
    ENDIF
!!$    use_time_range = (start_time /= 1) .OR. (end_time /= file_time_steps)
    dim1_steps = end_dim1 - start_dim1 + 1
!    write(0,*) "start,end time=", start_time, end_time
    IF (dim1_steps < 1) &
      & CALL finish(method_name, "number of slices of dimension 1 < 1")

    IF (PRESENT(fill_array)) THEN
      res => fill_array
    ELSE
      ALLOCATE( res(dim1_steps,var_size(2),var_size(3), &
                var_size(4)), stat=return_status )
      IF (return_status /= success) &
        CALL finish (method_name, &
          &          'ALLOCATE( netcdf_read_REAL_extdim_slice_extdim_extdim_extdim )')
      res(:,:,:,:)=0.0_wp
    ENDIF

!!$    ! check if the size is correct
!!$    IF (SIZE(netcdf_read_REAL_1D,1) < var_size(1)) &
!!$      CALL finish(method_name, "allocated size < var_size")
!!$    IF (SIZE(netcdf_read_REAL_1D,1) > var_size(1)) &
!!$      CALL warning(method_name, "allocated size > var_size")

    IF( my_process_is_mpi_workroot()) THEN
      start_read_index = (/start_dim1,1,1,1/)
      count_read_index = (/dim1_steps,var_size(2),var_size(3),var_size(4)/)
      CALL nf(nf90_get_var(file_id, varid, res(:,:,:,:), &
        &                  start_read_index, count_read_index), &
        &     variable_name)
    ENDIF

    ! broadcast...
    CALL p_bcast(res, process_mpi_root_id, p_comm_work)

  END FUNCTION netcdf_read_REAL_extdim_slice_extdim_extdim_extdim
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_INT_2D(file_id, variable_name, fill_array, &
    &                         n_g, scatter_pattern) result(res)
    INTEGER, POINTER                 :: res(:,:)

    INTEGER, INTENT(IN)              :: file_id
    CHARACTER(LEN=*), INTENT(IN)     :: variable_name
    define_fill_target_int           :: fill_array(:,:)
    INTEGER, INTENT(IN)              :: n_g
    CLASS(t_scatterPattern), POINTER :: scatter_pattern

    TYPE(t_ptr_2d_int) :: fill_arrays(1)
    TYPE(t_p_scatterPattern) :: scatter_pattern_(1)
    TYPE(t_ptr_2d_int) :: results(1)

    scatter_pattern_(1)%p => scatter_pattern

    IF (PRESENT(fill_array)) THEN
      fill_arrays(1)%p => fill_array
      results = netcdf_read_INT_2D_multivar(file_id=file_id, &
        &                                   variable_name=variable_name,&
        &                                   n_vars=1, &
        &                                   fill_arrays=fill_arrays, &
        &                                   n_g=n_g, &
        &                                   scatter_patterns=scatter_pattern_)
    ELSE
      results = netcdf_read_INT_2D_multivar(file_id=file_id, &
        &                                   variable_name=variable_name,&
        &                                   n_vars=1, n_g=n_g, &
        &                                   scatter_patterns=scatter_pattern_)
    END IF

    res => results(1)%p

  END FUNCTION netcdf_read_INT_2D

  FUNCTION netcdf_read_INT_2D_multivar(file_id, variable_name, n_vars, &
    &                                  fill_arrays, n_g, scatter_patterns) &
    result(res)

    INTEGER, INTENT(IN)                    :: n_vars
    INTEGER, INTENT(IN)                    :: file_id
    CHARACTER(LEN=*), INTENT(IN)           :: variable_name
    TYPE(t_ptr_2d_int), OPTIONAL        :: fill_arrays(n_vars)
    INTEGER, INTENT(IN)                    :: n_g
    TYPE(t_p_scatterPattern),INTENT(INOUT) :: scatter_patterns(n_vars)

    TYPE(t_ptr_2d_int)                  :: res(n_vars)

    INTEGER :: varid, var_type, var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)
    INTEGER :: return_status, i
    INTEGER, ALLOCATABLE :: tmp_array(:)

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_INT_2D_multivar'

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, &
        &                 var_size, var_dim_name)

      ! check if the dims look ok
      IF (var_dims /= 1 .OR. var_size(1) /= n_g) THEN
        write(0,*) "var_dims = ", var_dims, " var_size=", var_size, " n_g=", n_g
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

      IF (var_type /= NF90_INT) CALL finish(method_name, "invalid var_type")

    ENDIF

    ALLOCATE(tmp_array(MERGE(n_g, 0, my_process_is_mpi_workroot())), &
      &      stat=return_status )
    IF (return_status /= success) CALL finish(method_name,'ALLOCATE(tmp_array)')

    IF( my_process_is_mpi_workroot()) &
      CALL nf(nf90_get_var(file_id, varid, tmp_array(:)), variable_name)

    DO i = 1, n_vars
      IF (PRESENT(fill_arrays)) THEN
        res(i)%p => fill_arrays(i)%p
      ELSE
        ALLOCATE( res(i)%p(nproma, &
          &                   (scatter_patterns(i)%p%myPointCount - 1)/nproma + 1), &
          &       stat=return_status )
        IF (return_status /= success) THEN
          CALL finish (method_name, 'ALLOCATE( res )')
        ENDIF
        res(i)%p(:,:) = 0
      ENDIF

      CALL scatter_patterns(i)%p%distribute(tmp_array, res(i)%p, .FALSE.)
    END DO

  END FUNCTION netcdf_read_INT_2D_multivar
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_2D_all(file_id, variable_name, fill_array) &
    result(res)

    REAL(wp), POINTER            :: res(:,:)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:)

    INTEGER :: varid, var_type, var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)
    INTEGER :: return_status

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_REAL_2D_all'

    ! trivial return value.
    NULLIFY(res)

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, &
        &                 var_size, var_dim_name)

      ! check if the dims look ok
      IF (var_dims /= 2 ) THEN
        write(0,*) "var_dims = ", var_dims, " var_size=(", var_size(1), var_size(2), ")"
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

    ENDIF

    ! we need to sync the var_size...
    CALL p_bcast(var_size(1:2), process_mpi_root_id, p_comm_work)

    IF (PRESENT(fill_array)) THEN
      res => fill_array
    ELSE
      ALLOCATE( res(var_size(1),var_size(2)), stat=return_status )
      IF (return_status /= success) THEN
        CALL finish (method_name, 'ALLOCATE( netcdf_read_REAL_2D_all )')
      ENDIF
    ENDIF

    ! check if the size is correct
    IF (SIZE(res,1) < var_size(1)) &
      CALL finish(method_name, "allocated size < var_size")
    IF (SIZE(res,1) > var_size(1)) &
      CALL warning(method_name, "allocated size > var_size")
    IF (SIZE(res,2) < var_size(2)) &
      CALL finish(method_name, "allocated size < var_size")
    IF (SIZE(res,2) > var_size(2)) &
      CALL warning(method_name, "allocated size > var_size")

    IF( my_process_is_mpi_workroot()) THEN
      CALL nf(nf90_get_var(file_id, varid, res(:,:)), variable_name)
    ENDIF

    ! broadcast...
    CALL p_bcast(res, process_mpi_root_id, p_comm_work)

  END FUNCTION netcdf_read_REAL_2D_all
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_2D(file_id, variable_name, fill_array, &
    &                          n_g, scatter_pattern) result(res)

    REAL(wp), POINTER                :: res(:,:)

    INTEGER, INTENT(IN)              :: file_id
    CHARACTER(LEN=*), INTENT(IN)     :: variable_name
    define_fill_target               :: fill_array(:,:)
    INTEGER, INTENT(IN)              :: n_g
    CLASS(t_scatterPattern), POINTER :: scatter_pattern

    TYPE(t_ptr_2d) :: fill_arrays(1)
    TYPE(t_p_scatterPattern) :: scatter_pattern_(1)
    TYPE(t_ptr_2d) :: results(1)

    scatter_pattern_(1)%p => scatter_pattern

    IF (PRESENT(fill_array)) THEN
      fill_arrays(1)%p => fill_array
      results = netcdf_read_REAL_2D_multivar(file_id=file_id, &
        &                                    variable_name=variable_name,&
        &                                    n_vars=1, &
        &                                    fill_arrays=fill_arrays, &
        &                                    n_g=n_g, &
        &                                    scatter_patterns=scatter_pattern_)
    ELSE
      results = netcdf_read_REAL_2D_multivar(file_id=file_id, &
        &                                    variable_name=variable_name,&
        &                                    n_vars=1, n_g=n_g, &
        &                                    scatter_patterns=scatter_pattern_)
    END IF

    res => results(1)%p

  END FUNCTION netcdf_read_REAL_2D

  FUNCTION netcdf_read_REAL_2D_multivar(file_id, variable_name, n_vars, &
    &                                   fill_arrays, n_g, scatter_patterns) &
    result(res)

    INTEGER, INTENT(IN)                     :: n_vars
    INTEGER, INTENT(IN)                     :: file_id
    CHARACTER(LEN=*), INTENT(IN)            :: variable_name
    TYPE(t_ptr_2d), OPTIONAL          :: fill_arrays(n_vars)
    INTEGER, INTENT(IN)                     :: n_g
    TYPE(t_p_scatterPattern), INTENT(INOUT) :: scatter_patterns(n_vars)

    TYPE(t_ptr_2d)                    :: res(n_vars)

    INTEGER :: varid, var_type(1), var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)
    INTEGER :: return_status, i
    REAL(sp), ALLOCATABLE :: tmp_array_sp(:)
    REAL(dp), ALLOCATABLE :: tmp_array_dp(:)

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_REAL_2D_multivar'

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type(1), var_dims, &
        &                 var_size, var_dim_name)

      ! check if the dims look ok
      IF (var_dims /= 1 .OR. var_size(1) /= n_g) THEN
        write(0,*) "var_dims = ", var_dims, " var_size=", var_size, " n_g=", n_g
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF
    ENDIF

    CALL p_bcast(var_type, process_mpi_root_id, p_comm_work)

    IF( my_process_is_mpi_workroot()) THEN
      SELECT CASE(var_type(1))
        CASE(NF90_DOUBLE)
          ALLOCATE(tmp_array_dp(n_g), stat=return_status )
          CALL nf(nf90_get_var(file_id, varid, tmp_array_dp), variable_name)
        CASE(NF90_FLOAT, NF90_INT)
          ALLOCATE(tmp_array_sp(n_g), stat=return_status )
          CALL nf(nf90_get_var(file_id, varid, tmp_array_sp), variable_name)
        CASE default
          CALL finish(method_name, "incompatible var_type")
      END SELECT
    ELSE
      ALLOCATE(tmp_array_dp(0), tmp_array_sp(0), stat=return_status)
    ENDIF
    IF (return_status /= success) CALL finish(method_name,'ALLOCATE(tmp_array)')

    DO i = 1, n_vars
      IF (PRESENT(fill_arrays)) THEN
        res(i)%p => fill_arrays(i)%p
      ELSE
        ALLOCATE( res(i)%p(nproma, &
          &                   (scatter_patterns(i)%p%myPointCount - 1)/nproma + 1), &
          &       stat=return_status )
        IF (return_status /= success) THEN
          CALL finish (method_name, 'ALLOCATE( res )')
        ENDIF
        res(i)%p(:,:) = 0.0_wp
      ENDIF

      IF (var_type(1) == NF90_DOUBLE) THEN
        CALL scatter_patterns(i)%p%distribute(tmp_array_dp, res(i)%p, .FALSE.)
      ELSE
        CALL scatter_patterns(i)%p%distribute(tmp_array_sp, res(i)%p, .FALSE.)
      END IF
    END DO

  END FUNCTION netcdf_read_REAL_2D_multivar
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, n_g) fortran-style: O3(n_g, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, blocks, time)
  FUNCTION netcdf_read_REAL_2D_time(file_id, variable_name, fill_array, &
    &                               n_g, scatter_pattern, start_timestep, &
    &                               end_timestep) result(res)
    REAL(wp), POINTER                :: res(:,:,:)

    INTEGER, INTENT(IN)              :: file_id
    CHARACTER(LEN=*), INTENT(IN)     :: variable_name
    define_fill_target               :: fill_array(:,:,:)
    INTEGER, INTENT(IN)              :: n_g
    CLASS(t_scatterPattern), POINTER :: scatter_pattern
    INTEGER, INTENT(in), OPTIONAL    :: start_timestep, end_timestep

    res => netcdf_read_REAL_2D_extdim( &
      & file_id=file_id, variable_name=variable_name, fill_array=fill_array, &
      & n_g=n_g, scatter_pattern=scatter_pattern, start_extdim=start_timestep, &
      & end_extdim=end_timestep, extdim_name="time" )

  END FUNCTION netcdf_read_REAL_2D_time
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, ncells) fortran-style: O3(ncells, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, blocks, time)
  FUNCTION netcdf_read_REAL_2D_extdim(file_id, variable_name, &
    &                                 fill_array, n_g, scatter_pattern, &
    &                                 start_extdim, end_extdim, &
    &                                 extdim_name ) result(res)

    REAL(wp), POINTER                      :: res(:,:,:)

    INTEGER, INTENT(IN)                    :: file_id
    CHARACTER(LEN=*), INTENT(IN)           :: variable_name
    define_fill_target                     :: fill_array(:,:,:)
    INTEGER, INTENT(IN)                    :: n_g
    CLASS(t_scatterPattern), POINTER       :: scatter_pattern
    INTEGER, INTENT(in), OPTIONAL          :: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name

    TYPE(t_ptr_3d) :: fill_arrays(1)
    TYPE(t_p_scatterPattern) :: scatter_pattern_(1)
    TYPE(t_ptr_3d) :: results(1)

    scatter_pattern_(1)%p => scatter_pattern

    IF (PRESENT(fill_array)) THEN
      fill_arrays(1)%p => fill_array
      results = netcdf_read_REAL_2D_extdim_multivar( &
        file_id=file_id, variable_name=variable_name, n_vars=1, &
        fill_arrays=fill_arrays,  n_g=n_g, scatter_patterns=scatter_pattern_, &
        start_extdim=start_extdim, end_extdim=end_extdim, &
        extdim_name=extdim_name)
    ELSE
      results = netcdf_read_REAL_2D_extdim_multivar( &
        file_id=file_id, variable_name=variable_name, n_vars=1, n_g=n_g, &
        scatter_patterns=scatter_pattern_, start_extdim=start_extdim, &
        end_extdim=end_extdim, extdim_name=extdim_name)
    END IF

    res => results(1)%p

  END FUNCTION netcdf_read_REAL_2D_extdim

  FUNCTION netcdf_read_REAL_2D_extdim_multivar(file_id, variable_name,  &
    &                                          n_vars, fill_arrays, n_g,&
    &                                          scatter_patterns, start_extdim, &
    &                                          end_extdim, extdim_name) &
    result(res)

    INTEGER, INTENT(IN)                     :: n_vars
    INTEGER, INTENT(IN)                     :: file_id
    CHARACTER(LEN=*), INTENT(IN)            :: variable_name
    TYPE(t_ptr_3d), OPTIONAL          :: fill_arrays(n_vars)
    INTEGER, INTENT(IN)                     :: n_g
    TYPE(t_p_scatterPattern), INTENT(INOUT) :: scatter_patterns(n_vars)
    INTEGER, INTENT(in), OPTIONAL           :: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL  :: extdim_name

    TYPE(t_ptr_3d)          :: res(n_vars)

    INTEGER :: varid, var_type(1), var_dims
    INTEGER, TARGET :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)

    INTEGER :: file_time_steps, time_steps, start_time, end_time
    INTEGER :: start_read_index(2), count_read_index(2)

    INTEGER :: return_status, i, t
    REAL(sp), ALLOCATABLE :: tmp_array_sp(:)
    REAL(dp), ALLOCATABLE :: tmp_array_dp(:)
    REAL(wp), POINTER :: tmp_res(:,:)
    CLASS(t_scatterPattern), POINTER :: scatter_pattern_

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_REAL_2D_extdim_multivar'

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type(1), var_dims,&
        &                 var_size, var_dim_name)

      ! check if we have indeed 3 dimensions
      IF (var_dims /= 2 ) THEN
        WRITE(0,*) variable_name, ": var_dims = ", var_dims
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

      ! check if the input has the right shape/size
      IF ( var_size(1) /= n_g) THEN
        WRITE(0,*) variable_name, ": var_dims = ", var_dims, " var_size=", &
          &        var_size, " n_g=", n_g
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

    ENDIF

    CALL p_bcast(var_type, process_mpi_root_id, p_comm_work)

    ! we need to sync the var_size...
    CALL p_bcast(var_size(1:2), process_mpi_root_id, p_comm_work)
    file_time_steps      = var_size(2)

    ! calculate time range
    IF (PRESENT(start_extdim)) THEN
      start_time = start_extdim
    ELSE
      start_time = 1
    ENDIF
    IF (PRESENT(end_extdim)) THEN
      end_time = end_extdim
    ELSE
      end_time = file_time_steps
    ENDIF
    time_steps = end_time - start_time + 1
!    write(0,*) "start,end time=", start_time, end_time
    IF (time_steps < 1) &
      & CALL finish(method_name, "ext dim size < 1")

    DO i = 1, n_vars
      !-----------------------
      IF (PRESENT(fill_arrays)) THEN
        res(i)%p => fill_arrays(i)%p
        IF (SIZE(res(i)%p,3) < time_steps) &
          CALL finish(method_name, "allocated size < time_steps")
      ELSE
        ALLOCATE(res(i)%p(nproma, &
          &                  (scatter_patterns(i)%p%myPointCount - 1) / nproma + 1, &
          &                  time_steps), stat=return_status)
        IF (return_status /= success) THEN
          CALL finish (method_name, 'ALLOCATE( res )')
        ENDIF
        res(i)%p(:,:,:) = 0.0_wp
      ENDIF
    END DO

    IF( my_process_is_mpi_workroot()) THEN
      SELECT CASE(var_type(1))
        CASE(NF90_DOUBLE)
          ALLOCATE(tmp_array_dp(n_g), stat=return_status )
        CASE(NF90_FLOAT, NF90_INT)
          ALLOCATE(tmp_array_sp(n_g), stat=return_status )
        CASE default
          CALL finish(method_name, "incompatible var_type")
      END SELECT
    ELSE
      ALLOCATE(tmp_array_dp(0), tmp_array_sp(0), stat=return_status)
    ENDIF
    IF (return_status /= success) CALL finish(method_name,'ALLOCATE(tmp_array)')

    DO t = 1, time_steps
      IF( my_process_is_mpi_workroot()) THEN

        start_read_index = (/ 1, start_time+t-1 /)
        count_read_index = (/ n_g, 1/)
        IF (var_type(1) == NF90_DOUBLE) THEN
          CALL nf(nf90_get_var(file_id, varid, tmp_array_dp(:), &
            &                  start_read_index, count_read_index), &
            &     variable_name)
        ELSE
          CALL nf(nf90_get_var(file_id, varid, tmp_array_sp(:), &
            &                  start_read_index, count_read_index), &
            &     variable_name)
        END IF
      ENDIF


      DO i = 1, n_vars
        tmp_res => res(i)%p(:,:,LBOUND(res(i)%p, 3)+t-1)
        ! this is needed by PGI, it causes an internal compiler error if
        ! scatter_patterns(i)%p is used directly
        scatter_pattern_ => scatter_patterns(i)%p
        IF (var_type(1) == NF90_DOUBLE) THEN
          CALL scatter_pattern_%distribute(tmp_array_dp, tmp_res, .FALSE.)
        ELSE
          CALL scatter_pattern_%distribute(tmp_array_sp, tmp_res, .FALSE.)
        END IF
      END DO
    END DO

  END FUNCTION netcdf_read_REAL_2D_extdim_multivar
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, ncells) fortran-style: O3(ncells, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, blocks, time)
  FUNCTION netcdf_read_INT_2D_extdim(file_id, variable_name, &
    &                                fill_array, n_g, scatter_pattern, &
    &                                start_extdim, end_extdim, &
    &                                extdim_name ) result(res)

    INTEGER, POINTER                       :: res(:,:,:)

    INTEGER, INTENT(IN)                    :: file_id
    CHARACTER(LEN=*), INTENT(IN)           :: variable_name
    define_fill_target_int                 :: fill_array(:,:,:)
    INTEGER, INTENT(IN)                    :: n_g
    CLASS(t_scatterPattern), POINTER       :: scatter_pattern
    INTEGER, INTENT(in), OPTIONAL          :: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name

    TYPE(t_ptr_3d_int) :: fill_arrays(1)
    TYPE(t_p_scatterPattern) :: scatter_pattern_(1)
    TYPE(t_ptr_3d_int) :: results(1)

    scatter_pattern_(1)%p => scatter_pattern

    IF (PRESENT(fill_array)) THEN
      fill_arrays(1)%p => fill_array
      results = netcdf_read_INT_2D_extdim_multivar( &
        file_id=file_id, variable_name=variable_name, n_vars=1, &
        fill_arrays=fill_arrays,  n_g=n_g, scatter_patterns=scatter_pattern_, &
        start_extdim=start_extdim, end_extdim=end_extdim, &
        extdim_name=extdim_name)
    ELSE
      results = netcdf_read_INT_2D_extdim_multivar( &
        file_id=file_id, variable_name=variable_name, n_vars=1, n_g=n_g, &
        scatter_patterns=scatter_pattern_, start_extdim=start_extdim, &
        end_extdim=end_extdim, extdim_name=extdim_name)
    END IF

    res => results(1)%p

  END FUNCTION netcdf_read_INT_2D_extdim

  FUNCTION netcdf_read_INT_2D_extdim_multivar(file_id, variable_name, &
    &                                         n_vars, fill_arrays, n_g, &
    &                                         scatter_patterns, start_extdim, &
    &                                         end_extdim, extdim_name ) &
    result(res)

    INTEGER, INTENT(IN)                     :: n_vars
    INTEGER, INTENT(IN)                     :: file_id
    CHARACTER(LEN=*), INTENT(IN)            :: variable_name
    TYPE(t_ptr_3d_int), OPTIONAL         :: fill_arrays(n_vars)
    INTEGER, INTENT(IN)                     :: n_g
    TYPE(t_p_scatterPattern), INTENT(INOUT) :: scatter_patterns(n_vars)
    INTEGER, INTENT(in), OPTIONAL           :: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL  :: extdim_name

    TYPE(t_ptr_3d_int)          :: res(n_vars)

    INTEGER :: varid, var_type, var_dims
    INTEGER, TARGET :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)

    INTEGER :: file_time_steps, time_steps, start_time, end_time
    INTEGER :: start_read_index(2), count_read_index(2)

    INTEGER :: return_status, i, t
    INTEGER, ALLOCATABLE :: tmp_array(:)
    INTEGER, POINTER :: tmp_res(:,:)

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_INT_2D_extdim_multivar'

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, &
        &                 var_size, var_dim_name)

      ! check if we have indeed 3 dimensions
      IF (var_dims /= 2 ) THEN
        WRITE(0,*) variable_name, ": var_dims = ", var_dims
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

      ! check if the input has the right shape/size
      IF ( var_size(1) /= n_g) THEN
        WRITE(0,*) variable_name, ": var_dims = ", var_dims, " var_size=", &
          &        var_size, " n_g=", n_g
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

    ENDIF

    ! we need to sync the var_size...
    CALL p_bcast(var_size(1:2), process_mpi_root_id, p_comm_work)
    file_time_steps      = var_size(2)

    ! calculate time range
    IF (PRESENT(start_extdim)) THEN
      start_time = start_extdim
    ELSE
      start_time = 1
    ENDIF
    IF (PRESENT(end_extdim)) THEN
      end_time = end_extdim
    ELSE
      end_time = file_time_steps
    ENDIF
    time_steps = end_time - start_time + 1
!    write(0,*) "start,end time=", start_time, end_time
    IF (time_steps < 1) &
      & CALL finish(method_name, "ext dim size < 1")

    ALLOCATE(tmp_array(MERGE(n_g, 0, my_process_is_mpi_workroot())), &
      &      stat=return_status )
    IF (return_status /= success) THEN
      CALL finish (method_name, 'ALLOCATE( tmp_array )')
    ENDIF

    DO i = 1, n_vars
      IF (PRESENT(fill_arrays)) THEN
        res(i)%p => fill_arrays(i)%p
        IF (SIZE(res(i)%p,3) < time_steps) &
          CALL finish(method_name, "allocated size < time_steps")
      ELSE
        ALLOCATE(res(i)%p(nproma, &
          &                  (scatter_patterns(i)%p%myPointCount - 1) / nproma + 1, &
          &                  time_steps), stat=return_status)
        IF (return_status /= success) CALL finish (method_name, 'ALLOCATE(res)')
        res(i)%p(:,:,:) = 0
      END IF
    END DO

    DO t = 1, time_steps

      IF( my_process_is_mpi_workroot()) THEN

        start_read_index = (/ 1, start_time + t - 1 /)
        count_read_index = (/ n_g, 1 /)
        CALL nf(nf90_get_var(file_id, varid, tmp_array(:), &
          &                  start_read_index, count_read_index), &
          &     variable_name)
      ENDIF

      DO i = 1, n_vars
        tmp_res => res(i)%p(:,:,LBOUND(res(i)%p, 3)+t-1)
        CALL scatter_patterns(i)%p%distribute(tmp_array, tmp_res, .FALSE.)
      END DO
    END DO

  END FUNCTION netcdf_read_INT_2D_extdim_multivar
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  FUNCTION netcdf_read_REAL_3D_all(file_id, variable_name, fill_array) &
    result(res)

    REAL(wp), POINTER            :: res(:,:,:)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:)

    INTEGER :: varid, var_type, var_dims
    INTEGER :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)
    INTEGER :: return_status

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_REAL_3D_all'

    ! trivial return value.
    NULLIFY(res)

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type, var_dims, &
        &                 var_size, var_dim_name)

      ! check if the dims look ok
      IF (var_dims /= 3 ) THEN
        write(0,*) "var_dims = ", var_dims, " var_size=(", var_size(1), var_size(2), var_size(3), ")"
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

    ENDIF

    ! we need to sync the var_size...
    CALL p_bcast(var_size(1:3), process_mpi_root_id, p_comm_work)

    IF (PRESENT(fill_array)) THEN
      res => fill_array
    ELSE
      ALLOCATE( res(var_size(1),var_size(2),var_size(3)), stat=return_status )
      IF (return_status /= success) THEN
        CALL finish (method_name, 'ALLOCATE( netcdf_read_REAL_3D_all )')
      ENDIF
    ENDIF

    ! check if the size is correct
    IF (SIZE(res,1) < var_size(1)) &
      CALL finish(method_name, "allocated size < var_size")
    IF (SIZE(res,1) > var_size(1)) &
      CALL warning(method_name, "allocated size > var_size")
    IF (SIZE(res,2) < var_size(2)) &
      CALL finish(method_name, "allocated size < var_size")
    IF (SIZE(res,2) > var_size(2)) &
      CALL warning(method_name, "allocated size > var_size")
    IF (SIZE(res,3) < var_size(3)) &
      CALL finish(method_name, "allocated size < var_size")
    IF (SIZE(res,3) > var_size(3)) &
      CALL warning(method_name, "allocated size > var_size")

    IF( my_process_is_mpi_workroot()) THEN
      CALL nf(nf90_get_var(file_id, varid, res(:,:,:)), variable_name)
    ENDIF

    ! broadcast...
    CALL p_bcast(res, process_mpi_root_id, p_comm_work)

  END FUNCTION netcdf_read_REAL_3D_all
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O2(levels, n_g) fortran-style: O2(n_g, levels)
  ! The fill_array  has the structure:
  !       fill_array(nproma, levels, blocks)
  FUNCTION netcdf_read_REAL_3D(file_id, variable_name, fill_array, n_g, &
    &                          scatter_pattern, levelsdim_name) result(res)

    REAL(wp), POINTER  :: res(:,:,:)

    INTEGER, INTENT(IN)              :: file_id
    CHARACTER(LEN=*), INTENT(IN)     :: variable_name
    define_fill_target               :: fill_array(:,:,:)
    INTEGER, INTENT(IN)              :: n_g
    CLASS(t_scatterPattern), POINTER :: scatter_pattern
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: levelsdim_name

    INTEGER :: varid, var_type(1), var_dims
    INTEGER, TARGET :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)

    INTEGER :: file_vertical_levels

    INTEGER :: return_status, i
    REAL(sp), ALLOCATABLE :: tmp_array_sp(:)
    REAL(dp), ALLOCATABLE :: tmp_array_dp(:)
    REAL(wp), POINTER :: res_level(:,:)

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_REAL_3D'

    NULLIFY(res)

    IF( my_process_is_mpi_workroot() ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type(1), var_dims, &
        &                 var_size, var_dim_name)

      ! check if we have indeed 2 dimensions
      IF (var_dims /= 2 ) THEN
        WRITE(0,*) variable_name, ": var_dims = ", var_dims
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

      ! check if the input has the right shape/size
      IF ( var_size(1) /= n_g) THEN
        WRITE(0,*) variable_name, ": var_dims = ", var_dims, " var_size=", &
          &        var_size, " n_g=", n_g
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

    ENDIF

    CALL p_bcast(var_type, process_mpi_root_id, p_comm_work)

    ! we need to sync the var_size...
    CALL p_bcast(var_size(1:2), process_mpi_root_id, p_comm_work)
    file_vertical_levels = var_size(2)

    !-----------------------
    IF (PRESENT(fill_array)) THEN
      res => fill_array
    ELSE
      ALLOCATE( res (nproma, file_vertical_levels, &
        &            (scatter_pattern%myPointCount - 1)/nproma + 1), &
        &       stat=return_status )
      IF (return_status /= success) THEN
        CALL finish (method_name, 'ALLOCATE( res )')
      ENDIF
      res(:,:,:) = 0.0_wp
    ENDIF

    IF( my_process_is_mpi_workroot()) THEN
      SELECT CASE(var_type(1))
        CASE(NF90_DOUBLE)
          ALLOCATE(tmp_array_dp(n_g), stat=return_status )
        CASE(NF90_FLOAT, NF90_INT)
          ALLOCATE(tmp_array_sp(n_g), stat=return_status )
        CASE default
          CALL finish(method_name, "incompatible var_type")
      END SELECT
    ELSE
      ALLOCATE(tmp_array_dp(0), tmp_array_sp(0), stat=return_status)
    ENDIF
    IF (return_status /= success) CALL finish(method_name,'ALLOCATE(tmp_array)')

    IF (file_vertical_levels /= SIZE(res,2)) &
      CALL finish(method_name, 'file_vertical_levels /= SIZE(fill_array,2)')
    !-----------------------

    DO i = 1, file_vertical_levels
      IF( my_process_is_mpi_workroot()) THEN
        IF (var_type(1) == NF90_DOUBLE) THEN
          CALL nf(nf90_get_var(file_id, varid, tmp_array_dp(:), &
            &                  (/1,i/), (/n_g,1/)), &
            &     variable_name)
        ELSE
          CALL nf(nf90_get_var(file_id, varid, tmp_array_sp(:), &
            &                  (/1,i/), (/n_g,1/)), &
            &     variable_name)
        END IF
      ENDIF

      res_level => res(:,i,:)
      IF (var_type(1) == NF90_DOUBLE) THEN
        CALL scatter_pattern%distribute(tmp_array_dp, res_level, .FALSE.)
      ELSE
        CALL scatter_pattern%distribute(tmp_array_sp, res_level, .FALSE.)
      END IF
    END DO

  END FUNCTION netcdf_read_REAL_3D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, levels, ncells) fortran-style: O3(ncells, levels, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, levels, blocks, time)
  FUNCTION netcdf_read_REAL_3D_time(file_id, variable_name, fill_array, &
    &                               n_g, scatter_pattern, start_timestep, &
    &                               end_timestep, levelsdim_name) &
    & result(res)

    REAL(wp), POINTER  :: res(:,:,:,:)

    INTEGER, INTENT(IN)              :: file_id
    CHARACTER(LEN=*), INTENT(IN)     :: variable_name
    define_fill_target               :: fill_array(:,:,:,:)
    INTEGER, INTENT(IN)              :: n_g
    CLASS(t_scatterPattern), POINTER :: scatter_pattern
    INTEGER, INTENT(in), OPTIONAL    :: start_timestep, end_timestep
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: levelsdim_name

    res => netcdf_read_REAL_3D_extdim(   &
      & file_id=file_id,                 &
      & variable_name=variable_name,     &
      & fill_array=fill_array, n_g=n_g,  &
      & scatter_pattern=scatter_pattern, &
      & start_extdim=start_timestep,     &
      & end_extdim=end_timestep,         &
      & levelsdim_name=levelsdim_name,   &
      & extdim_name="time")

  END FUNCTION netcdf_read_REAL_3D_time
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, levels, n_g) fortran-style: O3(n_g, levels, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, levels, blocks, time)
  FUNCTION netcdf_read_REAL_3D_extdim(file_id, variable_name, &
    &                                 fill_array, n_g, scatter_pattern, &
    &                                 start_extdim, end_extdim, &
    &                                 levelsdim_name, extdim_name ) &
    &  result(res)

    REAL(wp), POINTER  :: res(:,:,:,:)

    INTEGER, INTENT(IN)              :: file_id
    CHARACTER(LEN=*), INTENT(IN)     :: variable_name
    define_fill_target               :: fill_array(:,:,:,:)
    INTEGER, INTENT(IN)              :: n_g
    CLASS(t_scatterPattern), POINTER :: scatter_pattern
    INTEGER, INTENT(in), OPTIONAL    :: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name, levelsdim_name

    INTEGER :: varid, var_type(1), var_dims
    INTEGER, TARGET :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: var_dim_name(MAX_VAR_DIMS)

    INTEGER :: file_vertical_levels, file_time_steps, time_steps, start_time, &
      &        end_time
    INTEGER :: start_read_index(3), count_read_index(3)

    INTEGER :: return_status, i, tt
    REAL(sp), ALLOCATABLE :: tmp_array_sp(:)
    REAL(dp), ALLOCATABLE :: tmp_array_dp(:)
    REAL(wp), POINTER  :: res_level(:,:)

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_netcdf_broadcast_2:netcdf_read_REAL_3D_extdim'

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL netcdf_inq_var(file_id, variable_name, varid, var_type(1), var_dims, &
        &                 var_size, var_dim_name)

      ! check if we have indeed 3 dimensions
      IF (var_dims /= 3 ) THEN
        WRITE(0,*) variable_name, ": var_dims = ", var_dims
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

      ! check if the input has the right shape/size
      IF ( var_size(1) /= n_g) THEN
        WRITE(0,*) variable_name, ": var_dims = ", var_dims, " var_size=", &
          &        var_size, " n_g=", n_g
        CALL finish(method_name, "Dimensions mismatch")
      ENDIF

    ENDIF

    CALL p_bcast(var_type, process_mpi_root_id, p_comm_work)

    ! we need to sync the var_size...
    CALL p_bcast(var_size(1:3), process_mpi_root_id, p_comm_work)
    file_vertical_levels = var_size(2)
    file_time_steps      = var_size(3)

    ! calculate time range
    IF (PRESENT(start_extdim)) THEN
      start_time = start_extdim
    ELSE
      start_time = 1
    ENDIF
    IF (PRESENT(end_extdim)) THEN
      end_time = end_extdim
    ELSE
      end_time = file_time_steps
    ENDIF

    time_steps = end_time - start_time + 1
    IF (time_steps < 1) &
      & CALL finish(method_name, "extdim size < 1")

    IF( my_process_is_mpi_workroot()) THEN
      SELECT CASE(var_type(1))
        CASE(NF90_DOUBLE)
          ALLOCATE(tmp_array_dp(n_g), stat=return_status )
        CASE(NF90_FLOAT, NF90_INT)
          ALLOCATE(tmp_array_sp(n_g), stat=return_status )
        CASE default
          CALL finish(method_name, "incompatible var_type")
      END SELECT
    ELSE
      ALLOCATE(tmp_array_dp(0), tmp_array_sp(0), stat=return_status)
    ENDIF
    IF (return_status /= success) CALL finish(method_name,'ALLOCATE(tmp_array)')

    IF (PRESENT(fill_array)) THEN
      res => fill_array(:,:,:,1:time_steps)
    ELSE
      ALLOCATE(res(nproma, file_vertical_levels, &
        &          (scatter_pattern%myPointCount - 1)/nproma + 1, time_steps), &
        &      stat=return_status )
      IF (return_status /= success) THEN
        CALL finish (method_name, 'ALLOCATE( res )')
      ENDIF
      res(:,:,:,:) = 0.0_wp
    ENDIF

    IF (SIZE(res,4) < time_steps) &
      CALL finish(method_name, "allocated size < time_steps")

    IF (file_vertical_levels /= SIZE(res,2)) &
      CALL finish(method_name, 'file_vertical_levels /= SIZE(fill_array,2)')
    !-----------------------

    DO tt=1, time_steps
      DO i=1, file_vertical_levels
        IF( my_process_is_mpi_workroot()) THEN
          start_read_index = (/ 1, i, tt + start_time - 1 /)
          count_read_index      = (/ n_g, 1, 1 /)
          IF (var_type(1) == NF90_DOUBLE) THEN
            CALL nf(nf90_get_var(file_id, varid, tmp_array_dp(:), &
              &                  start_read_index, count_read_index), &
              &     variable_name)
          ELSE
            CALL nf(nf90_get_var(file_id, varid, tmp_array_sp(:), &
              &                  start_read_index, count_read_index), &
              &     variable_name)
          END IF
        ENDIF
        
        res_level => res(:,i,:,LBOUND(res, 4)+tt-1)
        IF (var_type(1) == NF90_DOUBLE) THEN
          CALL scatter_pattern%distribute(tmp_array_dp, res_level, .FALSE.)
        ELSE
          CALL scatter_pattern%distribute(tmp_array_sp, res_level, .FALSE.)
        END IF
      END DO
    END DO

  END FUNCTION netcdf_read_REAL_3D_extdim
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION netcdf_open_input(filename)
    CHARACTER(LEN=*), INTENT(IN) :: filename

    INTEGER :: file_id

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL nf(nf90_open(TRIM(filename), nf90_nowrite, file_id), TRIM(filename))
    ELSE
      file_id = -1 ! set it to an invalid value
    ENDIF

    netcdf_open_input = file_id

  END FUNCTION netcdf_open_input
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION netcdf_close(file_id)
    INTEGER, INTENT(IN) :: file_id

    netcdf_close = -1
    IF( my_process_is_mpi_workroot()  ) THEN
        netcdf_close = nf90_close(file_id)
    ENDIF

  END FUNCTION netcdf_close
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE netcdf_inq_var(file_id, name, varid, var_type, var_dims, var_size, &
    &                       var_dim_name)
    INTEGER, INTENT(IN) :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: name

    INTEGER, INTENT(OUT) :: varid, var_type, var_dims
    INTEGER, INTENT(OUT) :: var_size(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max), INTENT(OUT) :: var_dim_name(MAX_VAR_DIMS)

    INTEGER  :: number_of_attributes
    INTEGER :: var_dims_reference(MAX_VAR_DIMS)
    CHARACTER(LEN=filename_max) :: check_var_name
    INTEGER :: i

    IF ( .NOT. my_process_is_mpi_workroot() ) RETURN

    CALL nf(nf90_inq_varid(file_id, name, varid), name)
    CALL nf(nf90_inquire_variable(file_id, varid, check_var_name, var_type, var_dims, &
      &                           var_dims_reference, number_of_attributes), &
      &     check_var_name)
    DO i=1, var_dims
      CALL nf(nf90_inquire_dimension (file_id, var_dims_reference(i), len = var_size(i)), &
        &     check_var_name)
      CALL nf(nf90_inquire_dimension(file_id, var_dims_reference(i), name = var_dim_name(i)), &
        &     check_var_name)
    ENDDO
!     write(0,*) " Read var_dims, var_size:",  var_dims, var_size
!     write(0,*) " check_var_name:",  check_var_name
!     write(0,*) " name:", name

  END SUBROUTINE netcdf_inq_var
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE netcdf_get_missValue(file_id, variable_name, has_missValue, missValue)
    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    LOGICAL                      :: has_missValue
    REAL(wp)                     :: missValue

    REAL(dp) :: readMissValue
    REAL(wp) :: broadcastValue(2)    
    INTEGER :: varid, return_status

    IF( my_process_is_mpi_workroot()  ) THEN
 
      ! write(0,*) "netcdf_get_missValue...", TRIM(variable_name) 
      CALL nf(nf90_inq_varid(file_id, variable_name, varid), variable_name)
      ! write(0,*) TRIM(variable_name), " id=", varid
  
      return_status = nf90_get_att(file_id, varid, "missing_value", readMissValue)
      IF (return_status == nf90_noerr) THEN
        has_missValue = .true.
        broadcastValue(1) = 1.0_wp
      ELSE
        has_missValue = .false.
        readMissValue = 0.0_wp
        broadcastValue(1) = 0.0_wp
      ENDIF
      broadcastValue(2) = readMissValue
      write(0,*)  TRIM(variable_name), "read miss=", has_missValue, readMissValue

    ENDIF

    CALL p_bcast(broadcastValue, process_mpi_root_id, p_comm_work)
    
    IF (broadcastValue(1) == 0.0_wp) THEN
      has_missValue = .false.
    ELSE
      has_missValue = .true.
    ENDIF
    missValue = broadcastValue(2)

    ! write(0,*)  TRIM(variable_name), " miss=", has_missValue, missValue
     
  END SUBROUTINE netcdf_get_missValue

END MODULE mo_read_netcdf_broadcast_2
