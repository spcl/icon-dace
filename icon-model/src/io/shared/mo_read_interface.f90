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

! This module provides basic methods for reading
! a NetCDF file in parallel or sequential in a transparent way.

#define define_fill_target REAL(wp), TARGET, OPTIONAL
#define define_return_pointer REAL(wp), POINTER, OPTIONAL
#define define_fill_target_int INTEGER, TARGET, OPTIONAL
#define define_return_pointer_int INTEGER, POINTER, OPTIONAL

MODULE mo_read_interface

  USE mo_kind
  USE mo_exception,          ONLY: finish, message_text, message
  USE mo_io_config,         ONLY:  read_netcdf_broadcast_method, &
    & read_netcdf_distribute_method,  default_read_method
  USE mo_read_netcdf_broadcast_2, ONLY: netcdf_open_input, netcdf_close, &
    &                                   netcdf_read_2D, netcdf_read_2D_int, &
    &                                   netcdf_read_REAL_2D_all, &
    &                                   netcdf_read_2D_extdim, &
    &                                   netcdf_read_2D_extdim_int, &
    &                                   netcdf_read_REAL_3D_all, &
    &                                   netcdf_read_3D_extdim, &
    &                                   read_0D_real => netcdf_read_0D_real, &
    &                                   netcdf_read_1D, netcdf_read_3D, &
    &                                   netcdf_read_1D_extdim_time, &
    &                                   netcdf_read_1D_extdim_extdim_time, &
    &                                   netcdf_read_extdim_slice_extdim_extdim_extdim, &
    &                                   t_p_scatterPattern, &
    &                                   netcdf_get_missValue, &
    &                                   netcdf_read_inq_varexists
  USE mo_read_netcdf_distributed, ONLY: t_distrib_read_data, distrib_nf_open, &
    &                                   distrib_read, distrib_nf_close, &
    &                                   distrib_inq_var_dims, idx_lvl_blk, &
    &                                   idx_blk_time, distrib_nf_inq_varexists
  USE mo_fortran_tools, ONLY: t_ptr_2d, t_ptr_2d_int, t_ptr_3d, t_ptr_3d_int, &
    & t_ptr_4d
  USE mo_model_domain, ONLY: t_patch
  USE mo_parallel_config, ONLY: nproma, p_test_run
  USE mo_model_domain, ONLY: t_patch
  USE mo_communication, ONLY: t_scatterPattern
  USE mo_mpi, ONLY: p_comm_work_test, p_comm_work, p_io, &
       &            my_process_is_mpi_workroot, p_bcast
  USE mo_impl_constants, ONLY: on_cells, on_vertices, on_edges
  USE mo_netcdf_errhandler, ONLY: nf
  USE mo_netcdf
  !-------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_netcdf_broadcast_method, read_netcdf_distribute_method
  PUBLIC :: t_stream_id, p_t_patch

  PUBLIC :: openInputFile, closeFile

  PUBLIC :: read_0D_real
  PUBLIC :: read_1D
  PUBLIC :: read_1D_extdim_time
  PUBLIC :: read_1D_extdim_extdim_time
  PUBLIC :: read_extdim_slice_extdim_extdim_extdim
  PUBLIC :: read_2D
  PUBLIC :: read_bcast_REAL_2D
  PUBLIC :: read_2D_int
  PUBLIC :: read_2D_time
  PUBLIC :: read_2D_1time
  PUBLIC :: read_2D_1lev_1time
  PUBLIC :: read_3D
  PUBLIC :: read_bcast_REAL_3D
  PUBLIC :: read_3D_1time
  PUBLIC :: read_3D_time
  PUBLIC :: read_2D_extdim
  PUBLIC :: read_2D_extdim_int
  PUBLIC :: read_3D_extdim
  PUBLIC :: read_inq_varexists

  PUBLIC :: nf  ! temporary hack, wich allows to USE nf via this module. Currently required for ART.

  PUBLIC :: on_cells, on_vertices, on_edges

  !--------------------------------------------------------

  TYPE t_read_info

    ! data required for read_netcdf_distribute_method
    TYPE(t_distrib_read_data), POINTER :: dist_read_info

    ! data required for read_netcdf_broadcast_method
    INTEGER :: n_l !< number of local points
    INTEGER :: n_g !< number of global points
    CLASS(t_scatterPattern), POINTER :: scatter_pattern
  END TYPE

  TYPE p_t_patch
    TYPE(t_patch), POINTER :: p
  END TYPE p_t_patch

  TYPE t_stream_id
    INTEGER  :: file_id             ! netcdf file id, or similar
    INTEGER  :: return_status       ! the latest operation return status
    INTEGER  :: input_method        ! read_netcdf_broadcast_method,
                                    ! read_netcdf_distribute_method, etc

    TYPE(t_read_info), ALLOCATABLE :: read_info(:,:)
    
  END TYPE t_stream_id
  !--------------------------------------------------------

  INTERFACE openInputFile
    MODULE PROCEDURE openInputFile_dist
    MODULE PROCEDURE openInputFile_dist_multivar
    MODULE PROCEDURE openInputFile_bcast
  END INTERFACE openInputFile

  INTERFACE closeFile
    MODULE PROCEDURE closeFile_dist
    MODULE PROCEDURE closeFile_bcast
  END INTERFACE closeFile

  !--------------------------------------------------------
  ! naming convention:
  !--------------------------------------------------------
  ! - data distribution type
  !   - bcast -> all processes get the same data
  !   - dist  -> based on decomposition information each process gets only a
  !              part of the data
  ! - data type
  !   - REAL -> single or double precision real
  !   - INT -> integer
  ! - dimensions
  !   - 0D                    -> single value (only bcast)
  !   - 1D                    -> single column (dimensions: (nlev))
  !   - 1D_extdim_time        -> single column
  !                              (dimensions: (nlev,nextdim,ntime))
  !   - 1D_extdim_extime_time -> single column
  !                              (dimensions: (nlev,nextdim_a,nextdim_b,ntime))
  !   - 2D                    -> single slice (dimensions: (nproma,nblk))
  !   - 2D_1time              -> single slice; file contains time dimension with
  !                              size 1 (dimensions: (nproma,nblk))
  !   - 2D_1lev_1time         -> single slice; file contains level and time
  !                              dimension with size 1
  !                              (dimensions: (nproma,nblk))
  !   - 2D_time               -> single slice with multiple time steps
  !                              (dimensions: (nproma,nblk,ntime))
  !   - 2D_extdim             -> single slice with additional dimension
  !                              (dimensions: (nproma,nblk,nextdim))
  !   - 3D                    -> 3D data field (dimensions: (nproma,nlev,nblk))
  !   - 3D_1time              -> 3D data field; file contains time dimension
  !                              with size 1 (dimensions: (nproma,nlev,nblk))
  !   - 3D_time               -> 3D data field with multiple time steps
  !                              (dimensions: (nproma,nlev,nblk))
  !   - 3D_extdim             -> 3D data field with additional dimension
  !                              (dimensions: (nproma,nlev,nblk))
  ! - multivar
  !   - reads single array from file and distributes it to multiple array with
  !     independent decompositions
  !--------------------------------------------------------

  INTERFACE read_1D
    MODULE PROCEDURE read_bcast_REAL_1D
  END INTERFACE read_1D

  INTERFACE read_1D_extdim_time
    MODULE PROCEDURE read_bcast_REAL_1D_extdim_time
  END INTERFACE read_1D_extdim_time

  INTERFACE read_1D_extdim_extdim_time
    MODULE PROCEDURE read_bcast_REAL_1D_extdim_extdim_time
  END INTERFACE read_1D_extdim_extdim_time

  INTERFACE read_extdim_slice_extdim_extdim_extdim
    MODULE PROCEDURE read_bcast_REAL_extdim_slice_extdim_extdim_extdim
  END INTERFACE read_extdim_slice_extdim_extdim_extdim
  
  INTERFACE read_2D_int
    MODULE PROCEDURE read_dist_INT_2D
    MODULE PROCEDURE read_dist_INT_2D_multivar
  END INTERFACE read_2D_int

  INTERFACE read_2D
    MODULE PROCEDURE read_dist_REAL_2D
    MODULE PROCEDURE read_dist_REAL_2D_multivar
  END INTERFACE read_2D

  INTERFACE read_2D_1time
    MODULE PROCEDURE read_dist_REAL_2D_1time
  END INTERFACE read_2D_1time

  INTERFACE read_2D_1lev_1time
    MODULE PROCEDURE read_dist_REAL_2D_1lev_1time
  END INTERFACE read_2D_1lev_1time

  INTERFACE read_2D_time
    MODULE PROCEDURE read_dist_REAL_2D_time
  END INTERFACE read_2D_time

  INTERFACE read_2D_extdim
    MODULE PROCEDURE read_dist_REAL_2D_extdim
    MODULE PROCEDURE read_dist_REAL_2D_extdim_multivar
  END INTERFACE read_2D_extdim

  INTERFACE read_2D_extdim_int
    MODULE PROCEDURE read_dist_INT_2D_extdim
    MODULE PROCEDURE read_dist_INT_2D_extdim_multivar
  END INTERFACE read_2D_extdim_int

  INTERFACE read_3D
    MODULE PROCEDURE read_dist_REAL_3D
  END INTERFACE read_3D

  INTERFACE read_3D_1time
    MODULE PROCEDURE read_dist_REAL_3D_1time
  END INTERFACE read_3D_1time

  INTERFACE read_3D_time
    MODULE PROCEDURE read_dist_REAL_3D_time
  END INTERFACE read_3D_time

  INTERFACE read_3D_extdim
    MODULE PROCEDURE read_dist_REAL_3D_extdim
  END INTERFACE read_3D_extdim

  INTERFACE read_inq_varexists
    MODULE PROCEDURE inq_varexists_bcast
    MODULE PROCEDURE inq_varexists_dist
  END INTERFACE read_inq_varexists

CONTAINS

  FUNCTION inq_varexists_bcast(file_id, variable_name) result(ret)
    INTEGER,           INTENT(INOUT) :: file_id
    CHARACTER(LEN=*),  INTENT(IN   ) :: variable_name
    LOGICAL                          :: ret
    ret = netcdf_read_inq_varexists(file_id, variable_name)
  END FUNCTION inq_varexists_bcast

  FUNCTION inq_varexists_dist(stream_id, variable_name) result(ret)
    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    CHARACTER(LEN=*),  INTENT(IN   ) :: variable_name
    LOGICAL                          :: ret
    CHARACTER(LEN=*), PARAMETER  :: method_name = &
      'mo_read_interface:inq_varexists'

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      ret = netcdf_read_inq_varexists(stream_id%file_id, variable_name)
    CASE (read_netcdf_distribute_method)
      ret = distrib_nf_inq_varexists(stream_id%file_id, variable_name)
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT
  END FUNCTION inq_varexists_dist

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_bcast_REAL_1D(file_id, variable_name, fill_array, &
    &                           return_pointer)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:)
    define_return_pointer        :: return_pointer(:)

    REAL(wp), POINTER            :: tmp_pointer(:)
    CHARACTER(LEN=NF90_MAX_NAME)   :: variable_name_
    CHARACTER(LEN=*), PARAMETER  :: method_name = &
      'mo_read_interface:read_bcast_REAL_1D'

    ! make variable name available on all processes
    CALL bcast_varname(variable_name, variable_name_)

    ! check whether fill_array and/or return_pointer was provided
    IF (.NOT. (PRESENT(fill_array) .OR. PRESENT(return_pointer))) &
      CALL finish(method_name, "invalid arguments")

    ! there is only one implementation for this read routine type
    tmp_pointer => netcdf_read_1D(file_id, variable_name_, fill_array)
    IF (PRESENT(return_pointer)) return_pointer => tmp_pointer

  END SUBROUTINE read_bcast_REAL_1D
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_bcast_REAL_2D(file_id, variable_name, fill_array, &
    &                           return_pointer)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:)
    define_return_pointer        :: return_pointer(:,:)

    REAL(wp), POINTER            :: tmp_pointer(:,:)
    CHARACTER(LEN=NF90_MAX_NAME)   :: variable_name_
    CHARACTER(LEN=*), PARAMETER  :: method_name = &
      'mo_read_interface:read_bcast_REAL_2D'

    ! make variable name available on all processes
    CALL bcast_varname(variable_name, variable_name_)

    ! check whether fill_array and/or return_pointer was provided
    IF (.NOT. (PRESENT(fill_array) .OR. PRESENT(return_pointer))) &
      CALL finish(method_name, "invalid arguments")

    ! there is only one implementation for this read routine type
    tmp_pointer => netcdf_read_REAL_2D_all(file_id, variable_name_, fill_array)
    IF (PRESENT(return_pointer)) return_pointer => tmp_pointer

  END SUBROUTINE read_bcast_REAL_2D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_bcast_REAL_1D_extdim_time(file_id, variable_name, &
    &                                       fill_array, return_pointer, &
    &                                       dim_names, start_timestep, &
    &                                       end_timestep)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:)
    define_return_pointer        :: return_pointer(:,:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_names(:)
    INTEGER, INTENT(IN), OPTIONAL:: start_timestep, end_timestep

    REAL(wp), POINTER            :: tmp_pointer(:,:,:)
    CHARACTER(LEN=NF90_MAX_NAME)   :: variable_name_
    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:read_bcast_REAL_1D_extdim_time'

    ! make variable name available on all processes
    CALL bcast_varname(variable_name, variable_name_)

    ! check whether fill_array and/or return_pointer was provided
    IF (.NOT. (PRESENT(fill_array) .OR. PRESENT(return_pointer))) &
      CALL finish(method_name, "invalid arguments")

    ! there is only one implementation for this read routine type
    tmp_pointer => netcdf_read_1D_extdim_time(file_id, variable_name_, &
      &                                       fill_array, dim_names, &
      &                                       start_timestep, end_timestep)
    IF (PRESENT(return_pointer)) return_pointer => tmp_pointer

  END SUBROUTINE read_bcast_REAL_1D_extdim_time
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_bcast_REAL_1D_extdim_extdim_time( &
    file_id, variable_name, fill_array, return_pointer, dim_names, &
    start_timestep, end_timestep)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:,:)
    define_return_pointer        :: return_pointer(:,:,:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_names(:)
    INTEGER, INTENT(IN), OPTIONAL:: start_timestep, end_timestep

    REAL(wp), POINTER            :: tmp_pointer(:,:,:,:)
    CHARACTER(LEN=NF90_MAX_NAME)   :: variable_name_
    CHARACTER(LEN=*), PARAMETER  :: method_name = &
      'mo_read_interface:read_bcast_REAL_1D_extdim_extdim_time'

    ! make variable name available on all processes
    CALL bcast_varname(variable_name, variable_name_)

    ! check whether fill_array and/or return_pointer was provided
    IF (.NOT. (PRESENT(fill_array) .OR. PRESENT(return_pointer))) &
      CALL finish(method_name, "invalid arguments")

    ! there is only one implementation for this read routine type
    tmp_pointer => netcdf_read_1D_extdim_extdim_time(file_id, variable_name_, &
      &                                              fill_array, dim_names, &
      &                                              start_timestep, end_timestep)
    IF (PRESENT(return_pointer)) return_pointer => tmp_pointer

  END SUBROUTINE read_bcast_REAL_1D_extdim_extdim_time
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_bcast_REAL_extdim_slice_extdim_extdim_extdim ( &
    file_id, variable_name, fill_array, return_pointer, dim_names, &
    start_extdim1, end_extdim1)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:,:)
    define_return_pointer        :: return_pointer(:,:,:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_names(:)
    INTEGER, INTENT(IN), OPTIONAL:: start_extdim1, end_extdim1

    REAL(wp), POINTER            :: tmp_pointer(:,:,:,:)
    CHARACTER(LEN=NF90_MAX_NAME)   :: variable_name_
    CHARACTER(LEN=*), PARAMETER  :: method_name = &
      'mo_read_interface:read_bcast_REAL_extdim_extdim_extdim_extdim_slice'

    ! make variable name available on all processes
    CALL bcast_varname(variable_name, variable_name_)

    ! check whether fill_array and/or return_pointer was provided
    IF (.NOT. (PRESENT(fill_array) .OR. PRESENT(return_pointer))) &
      CALL finish(method_name, "invalid arguments")

    ! there is only one implementation for this read routine type
    tmp_pointer => netcdf_read_extdim_slice_extdim_extdim_extdim( &
                          file_id, variable_name_,                &
      &                   fill_array, dim_names,                  &
      &                   start_extdim1, end_extdim1              )
    IF (PRESENT(return_pointer)) return_pointer => tmp_pointer

  END SUBROUTINE read_bcast_REAL_extdim_slice_extdim_extdim_extdim
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_dist_INT_2D_multivar(stream_id, location, variable_name, &
    &                                  n_var, fill_array, return_pointer)

    TYPE(t_stream_id), INTENT(INOUT)       :: stream_id
    INTEGER, INTENT(IN)                    :: location
    CHARACTER(LEN=*), INTENT(IN)           :: variable_name
    INTEGER, INTENT(IN)                    :: n_var
    TYPE(t_ptr_2d_int), OPTIONAL        :: fill_array(:)
    TYPE(t_ptr_2d_int), OPTIONAL        :: return_pointer(:)

    TYPE(t_ptr_2d_int)                  :: tmp_return(n_var)
    TYPE(t_p_scatterPattern)               :: scatter_patterns(n_var)
    INTEGER                                :: n_g, i
    TYPE(t_ptr_2d_int), ALLOCATABLE     :: var_data_2d(:)
    CHARACTER(LEN=NF90_MAX_NAME)             :: variable_name_
    CHARACTER(LEN=*), PARAMETER            :: method_name = &
      'mo_read_interface:read_dist_INT_2D_multivar'
    TYPE(t_distrib_read_data) :: io_data(n_var)
    ! make variable name available on all processes
    CALL bcast_varname(variable_name, variable_name_)

    ! check whether fill_array and/or return_pointer was provided
    IF (.NOT. (PRESENT(fill_array) .OR. PRESENT(return_pointer))) &
      CALL finish(method_name, "invalid arguments")

    IF (SIZE(stream_id%read_info, 2) /= n_var) &
      CALL finish(method_name, &
        "number of distributions in stream_id does not match number of variables")

    CALL check_dimensions(stream_id%file_id, variable_name_, 1, &
      &                   (/stream_id%read_info(location, 1)%n_g/), location)

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      DO i = 1, n_var
        scatter_patterns(i)%p => &
          stream_id%read_info(location, i)%scatter_pattern
      END DO
      n_g = stream_id%read_info(location, 1)%n_g
      tmp_return = netcdf_read_2D_int(stream_id%file_id, variable_name_, n_var, &
          &                           fill_array, n_g, scatter_patterns)
      IF (PRESENT(return_pointer)) return_pointer(1:n_var) = tmp_return
    CASE (read_netcdf_distribute_method)

      ALLOCATE(var_data_2d(n_var))

      ! gather pointers of all output fields
      IF (PRESENT(fill_array)) THEN
        DO i = 1, n_var
          var_data_2d(i)%p => fill_array(i)%p
        END DO
      ELSE
        DO i = 1, n_var
          ALLOCATE(var_data_2d(i)%p(nproma, &
            (stream_id%read_info(location, i)%n_l - 1)/nproma + 1))
          var_data_2d(i)%p(:,:) = 0
        END DO
      ENDIF
      IF (PRESENT(return_pointer)) THEN
        DO i = 1, n_var
          return_pointer(i)%p => var_data_2d(i)%p
        END DO
      END IF
      DO i = 1, n_var
        io_data(i)%basic_data_index &
             = stream_id%read_info(location, i)%dist_read_info%basic_data_index
        io_data(i)%pat => stream_id%read_info(location, i)%dist_read_info%pat
      END DO
      CALL distrib_read(stream_id%file_id, variable_name_, var_data_2d, io_data)
      DEALLOCATE(var_data_2d)
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END SUBROUTINE read_dist_INT_2D_multivar

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_dist_INT_2D(stream_id, location, variable_name, fill_array, &
    &                         return_pointer)

    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    INTEGER, INTENT(IN)              :: location
    CHARACTER(LEN=*), INTENT(IN)     :: variable_name
    define_fill_target_int           :: fill_array(:,:)
    define_return_pointer_int        :: return_pointer(:,:)
    TYPE(t_ptr_2d_int) :: tmp_ptr(1)
    INTEGER, POINTER                 :: tmp_pointer(:,:)
    CHARACTER(LEN=NF90_MAX_NAME)       :: variable_name_
    CHARACTER(LEN=*), PARAMETER      :: method_name = &
      'mo_read_interface:read_dist_INT_2D'

    ! make variable name available on all processes
    CALL bcast_varname(variable_name, variable_name_)

    ! check whether fill_array and/or return_pointer was provided
    IF (.NOT. (PRESENT(fill_array) .OR. PRESENT(return_pointer))) &
      CALL finish(method_name, "invalid arguments")

    CALL check_dimensions(stream_id%file_id, variable_name_, 1, &
      &                   (/stream_id%read_info(location, 1)%n_g/), location)

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      tmp_pointer => &
        netcdf_read_2D_int(stream_id%file_id, variable_name_, fill_array, &
        stream_id%read_info(location, 1)%n_g, &
        stream_id%read_info(location, 1)%scatter_pattern)
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
    CASE (read_netcdf_distribute_method)
      IF (PRESENT(fill_array)) THEN
        tmp_pointer => fill_array
      ELSE
        ALLOCATE(tmp_pointer(nproma, &
          (stream_id%read_info(location, 1)%n_l - 1)/nproma + 1))
        tmp_pointer(:,:) = 0
      ENDIF
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
      tmp_ptr(1)%p => tmp_pointer
      CALL distrib_read(stream_id%file_id, variable_name_, tmp_ptr, &
        &               (/stream_id%read_info(location, 1)%dist_read_info/))
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END SUBROUTINE read_dist_INT_2D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_dist_REAL_2D_multivar(stream_id, location, variable_name, &
    &                                   n_var, fill_array, return_pointer)

    TYPE(t_stream_id), INTENT(INOUT)  :: stream_id
    INTEGER, INTENT(IN)               :: location
    CHARACTER(LEN=*), INTENT(IN)      :: variable_name
    INTEGER, INTENT(IN)               :: n_var
    TYPE(t_ptr_2d), OPTIONAL    :: fill_array(:)
    TYPE(t_ptr_2d), OPTIONAL    :: return_pointer(:)

    TYPE(t_ptr_2d)              :: tmp_return(n_var)
    TYPE(t_p_scatterPattern)          :: scatter_patterns(n_var)
    INTEGER                           :: n_g, i
    TYPE(t_ptr_2d), ALLOCATABLE :: var_data_2d(:)
    CHARACTER(LEN=NF90_MAX_NAME)        :: variable_name_
    CHARACTER(LEN=*), PARAMETER       :: method_name = &
      'mo_read_interface:read_dist_REAL_2D_multivar'
    TYPE(t_distrib_read_data) :: io_data(n_var)

    ! make variable name available on all processes
    CALL bcast_varname(variable_name, variable_name_)

    ! check whether fill_array and/or return_pointer was provided
    IF (.NOT. (PRESENT(fill_array) .OR. PRESENT(return_pointer))) &
      CALL finish(method_name, "invalid arguments")

    IF (SIZE(stream_id%read_info, 2) /= n_var) &
      CALL finish(method_name, &
        "number of stream ids does not match number of variables")

    CALL check_dimensions(stream_id%file_id, variable_name_, 1, &
      &                   (/stream_id%read_info(location, 1)%n_g/), location)

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      DO i = 1, n_var
        scatter_patterns(i)%p => &
          stream_id%read_info(location, i)%scatter_pattern
      END DO
      n_g = stream_id%read_info(location, 1)%n_g
      tmp_return = netcdf_read_2D(stream_id%file_id, variable_name_, n_var, &
            &                     fill_array, n_g, scatter_patterns)
      IF (PRESENT(return_pointer)) return_pointer(1:n_var) = tmp_return
    CASE (read_netcdf_distribute_method)

      ALLOCATE(var_data_2d(n_var))

      ! gather pointers of all output fields
      IF (PRESENT(fill_array)) THEN
        DO i = 1, n_var
          var_data_2d(i)%p => fill_array(i)%p
        END DO
      ELSE
        DO i = 1, n_var
          ALLOCATE(var_data_2d(i)%p(nproma, &
            (stream_id%read_info(location, i)%n_l - 1)/nproma + 1))
          var_data_2d(i)%p(:,:) = 0.0_wp
        END DO
      ENDIF
      IF (PRESENT(return_pointer)) THEN
        DO i = 1, n_var
          return_pointer(i)%p => var_data_2d(i)%p
        END DO
      END IF
      DO i = 1, n_var
        io_data(i)%basic_data_index &
             = stream_id%read_info(location, i)%dist_read_info%basic_data_index
        io_data(i)%pat => stream_id%read_info(location, i)%dist_read_info%pat
      END DO
      CALL distrib_read(stream_id%file_id, variable_name_, var_data_2d, io_data)
      DEALLOCATE(var_data_2d)
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END SUBROUTINE read_dist_REAL_2D_multivar

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_dist_REAL_2D(stream_id, location, variable_name, fill_array, &
    &                          return_pointer)

    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    INTEGER, INTENT(IN)          :: location
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:)
    define_return_pointer        :: return_pointer(:,:)
    TYPE(t_ptr_2d)               :: tmp_ptr(1)
    REAL(wp), POINTER            :: tmp_pointer(:,:)
    CHARACTER(LEN=NF90_MAX_NAME)   :: variable_name_
    CHARACTER(LEN=*), PARAMETER  :: method_name = &
      'mo_read_interface:read_dist_REAL_2D'

    ! make variable name available on all processes
    CALL bcast_varname(variable_name, variable_name_)

    ! check whether fill_array and/or return_pointer was provided
    IF (.NOT. (PRESENT(fill_array) .OR. PRESENT(return_pointer))) &
      CALL finish(method_name, "invalid arguments")

    CALL check_dimensions(stream_id%file_id, variable_name_, 1, &
      &                   (/stream_id%read_info(location, 1)%n_g/), location)

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      tmp_pointer => &
        netcdf_read_2D(stream_id%file_id, variable_name_, fill_array, &
        stream_id%read_info(location, 1)%n_g, &
        stream_id%read_info(location, 1)%scatter_pattern)
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
    CASE (read_netcdf_distribute_method)
      IF (PRESENT(fill_array)) THEN
        tmp_pointer => fill_array
      ELSE
        ALLOCATE(tmp_pointer(nproma, &
          (stream_id%read_info(location, 1)%n_l - 1)/nproma + 1))
        tmp_pointer(:,:) = 0.0_wp
      ENDIF
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
      tmp_ptr(1)%p => tmp_pointer
      CALL distrib_read(stream_id%file_id, variable_name_, tmp_ptr, &
        &               (/stream_id%read_info(location, 1)%dist_read_info/))
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END SUBROUTINE read_dist_REAL_2D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, n) fortran-style: O3(n, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, blocks)
  ! Since the array in the file has a time dimension we can map this case to
  ! read_dist_REAL_2D_extdim. However fill_array lacks the time dimension. To
  ! add this dimension we use an assumed-size array in
  ! read_dist_REAL_2D_1time_. In order to use assumed-size in this case
  ! we need the shape of the original fill_array. This is determined by
  ! read_dist_REAL_2D_1time.
  SUBROUTINE read_dist_REAL_2D_1time(stream_id, location, variable_name, &
    & fill_array, return_pointer,                                        &
    & has_missValue, missValue)
    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    INTEGER, INTENT(IN)              :: location
    CHARACTER(LEN=*), INTENT(IN)     :: variable_name
    define_fill_target               :: fill_array(:,:)
    define_return_pointer            :: return_pointer(:,:)  
    LOGICAL, OPTIONAL                :: has_missValue
    REAL(wp), OPTIONAL               :: missValue    
    
    CHARACTER(LEN=*), PARAMETER      :: method_name = &
      'mo_read_interface:read_dist_REAL_2D_1time'

    IF (PRESENT(fill_array)) THEN
      CALL read_dist_REAL_2D_1time_(stream_id, location, variable_name, &
        & SHAPE(fill_array), fill_array,                                &
        & return_pointer,                                               &
        & has_missValue=has_missValue,                                  &
        & missValue=missValue)
    ELSE
      CALL read_dist_REAL_2D_1time_(stream_id, location, variable_name, &
        &  (/0,0/), return_pointer=return_pointer,                      &
        & has_missValue=has_missValue,                                  &
        & missValue=missValue)
    END IF

  END SUBROUTINE read_dist_REAL_2D_1time


  SUBROUTINE read_dist_REAL_2D_1time_(stream_id, location, variable_name, &
    & array_shape, fill_array, return_pointer,                            &
    & has_missValue, missValue)
    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    INTEGER, INTENT(IN)              :: location
    CHARACTER(LEN=*), INTENT(IN)     :: variable_name
    INTEGER, INTENT(IN)              :: array_shape(2)
    define_fill_target               :: fill_array(array_shape(1), &
      &                                            array_shape(2), 1)
    define_return_pointer            :: return_pointer(:,:)
    LOGICAL, OPTIONAL                :: has_missValue
    REAL(wp), OPTIONAL               :: missValue    
    
    CHARACTER(LEN=*), PARAMETER      :: method_name = &
      'mo_read_interface:read_dist_REAL_2D_1time_'

    REAL(wp), POINTER :: return_pointer_(:,:,:)

    ! Since fill_array now has a time dimension we can call
    ! read_dist_REAL_2D_extdim
    IF (PRESENT(return_pointer)) THEN
      CALL read_dist_REAL_2D_extdim( &
        & stream_id=stream_id, location=location, variable_name=variable_name, &
        & fill_array=fill_array, return_pointer=return_pointer_, &
        & start_extdim=1, end_extdim=1, extdim_name="time" ,     &
        & has_missValue=has_missValue,                           &
        & missValue=missValue)

      ALLOCATE(return_pointer(SIZE(return_pointer_,1),SIZE(return_pointer_,2)))
      return_pointer(:,:) = return_pointer_(:,:,1)
      DEALLOCATE(return_pointer_)
    ELSE
      CALL read_dist_REAL_2D_extdim( &
        & stream_id=stream_id, location=location, variable_name=variable_name, &
        & fill_array=fill_array, start_extdim=1, end_extdim=1, &
        & extdim_name="time",                                  &
        & has_missValue=has_missValue,                           &
        & missValue=missValue)
    END IF

  END SUBROUTINE read_dist_REAL_2D_1time_
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, levels, n) fortran-style: O3(n, levels, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, blocks)
  ! Since the array in the file has a level and time dimension we can map this
  ! case to read_dist_REAL_3D_extdim. However fill_array lacks the level and
  ! time dimension. To add these dimensions we use an assumed-size array in
  ! read_dist_REAL_2D_1lev_1time_. In order to use assumed-size in this case we
  ! need the shape of the original fill_array. This is determined by
  ! read_dist_REAL_2D_1lev_1time.
  SUBROUTINE read_dist_REAL_2D_1lev_1time(stream_id, location, variable_name, &
    &                                     fill_array, return_pointer)
    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    INTEGER, INTENT(IN)              :: location
    CHARACTER(LEN=*), INTENT(IN)     :: variable_name
    define_fill_target               :: fill_array(:,:)
    define_return_pointer            :: return_pointer(:,:)
    CHARACTER(LEN=*), PARAMETER      :: method_name = &
      'mo_read_interface:read_dist_REAL_2D_1lev_1time'

    IF (PRESENT(fill_array)) THEN
      CALL read_dist_REAL_2D_1lev_1time_(stream_id, location, variable_name, &
        &                                SHAPE(fill_array), fill_array, &
        &                                return_pointer)
    ELSE
      CALL read_dist_REAL_2D_1lev_1time_(stream_id, location, variable_name, &
        &                                (/0,0/), return_pointer=return_pointer)
    END IF

  END SUBROUTINE read_dist_REAL_2D_1lev_1time

  SUBROUTINE read_dist_REAL_2D_1lev_1time_(stream_id, location, variable_name, &
    &                                      array_shape, fill_array, &
    &                                      return_pointer)
    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    INTEGER, INTENT(IN)              :: location
    CHARACTER(LEN=*), INTENT(IN)     :: variable_name
    INTEGER, INTENT(IN)              :: array_shape(2)
    define_fill_target               :: fill_array(array_shape(1), 1, &
      &                                            array_shape(2), 1)
    define_return_pointer            :: return_pointer(:,:)
    
    CHARACTER(LEN=*), PARAMETER      :: method_name = &
      'mo_read_interface:read_dist_REAL_2D_1lev_1time_'

    REAL(wp), POINTER :: return_pointer_(:,:,:,:)

    ! Since fill_array now has a level and time dimension we can call
    ! read_dist_REAL_3D_extdim
    IF (PRESENT(return_pointer)) THEN
      CALL read_dist_REAL_3D_extdim( &
        & stream_id=stream_id, location=location, variable_name=variable_name, &
        & fill_array=fill_array, return_pointer=return_pointer_, &
        & start_extdim=1, end_extdim=1, extdim_name="time" )

      ALLOCATE(return_pointer(SIZE(return_pointer_,1),SIZE(return_pointer_,3)))
      return_pointer(:,:) = return_pointer_(:,1,:,1)
      DEALLOCATE(return_pointer_)
    ELSE
      CALL read_dist_REAL_3D_extdim( &
        & stream_id=stream_id, location=location, variable_name=variable_name, &
        & fill_array=fill_array, start_extdim=1, end_extdim=1, &
        & extdim_name="time" )
    END IF

  END SUBROUTINE read_dist_REAL_2D_1lev_1time_
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, n) fortran-style: O3(n, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, blocks, time)
  ! We can map this case to read_dist_REAL_2D_extdim.
  SUBROUTINE read_dist_REAL_2D_time(stream_id, location, variable_name, &
    &  fill_array, return_pointer, start_timestep,  &
    &  end_timestep,                                &
    &  has_missValue, missValue)
    
    TYPE(t_stream_id), INTENT(INOUT) :: stream_id
    INTEGER, INTENT(IN)              :: location
    CHARACTER(LEN=*), INTENT(IN)     :: variable_name
    define_fill_target               :: fill_array(:,:,:)
    define_return_pointer            :: return_pointer(:,:,:)
    INTEGER, INTENT(in), OPTIONAL    :: start_timestep, end_timestep
    LOGICAL, OPTIONAL                :: has_missValue
    REAL(wp), OPTIONAL               :: missValue
    
    
    CHARACTER(LEN=*), PARAMETER      :: method_name = &
      'mo_read_interface:read_dist_REAL_2D_time'

    CALL read_dist_REAL_2D_extdim(&
      & stream_id=stream_id, location=location, variable_name=variable_name, &
      & fill_array=fill_array, return_pointer=return_pointer, &
      & start_extdim=start_timestep, end_extdim=end_timestep, &
      & extdim_name="time",                                   &
      & has_missValue=has_missValue,                          &
      & missValue=missValue)

  END SUBROUTINE read_dist_REAL_2D_time
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, n) fortran-style: O3(n, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, blocks, time)
  SUBROUTINE read_dist_REAL_2D_extdim(stream_id, location, variable_name, &
    &  fill_array, return_pointer, start_extdim,                          &
    &  end_extdim, extdim_name,                                           &
    &  has_missValue, missValue)

    TYPE(t_stream_id), INTENT(INOUT)       :: stream_id
    INTEGER, INTENT(IN)                    :: location
    CHARACTER(LEN=*), INTENT(IN)           :: variable_name
    define_fill_target                     :: fill_array(:,:,:)
    define_return_pointer                  :: return_pointer(:,:,:)
    INTEGER, INTENT(in), OPTIONAL          :: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name
    LOGICAL, OPTIONAL                      :: has_missValue
    REAL(wp), OPTIONAL                     :: missValue
    REAL(wp), POINTER                      :: tmp_pointer(:,:,:)
    TYPE(t_ptr_3d)                         :: tmp_ptr(1)
    INTEGER                                :: var_dimlen(2), var_start(2), &
      &                                       var_end(2), var_ndims
    CHARACTER(LEN=NF90_MAX_NAME)             :: variable_name_
    CHARACTER(LEN=*), PARAMETER            :: method_name = &
      'mo_read_interface:read_dist_REAL_2D_extdim'

    ! make variable name available on all processes
    CALL bcast_varname(variable_name, variable_name_)


    var_dimlen(:) = (/stream_id%read_info(location, 1)%n_g, -1/)
    IF (PRESENT(fill_array)) THEN
      var_dimlen(2) = SIZE(fill_array, 3)
    END IF

    ! check whether fill_array and/or return_pointer was provided
    IF (.NOT. (PRESENT(fill_array) .OR. PRESENT(return_pointer))) &
      CALL finish(method_name, "invalid arguments")

    IF (PRESENT(start_extdim) .NEQV. PRESENT(end_extdim)) &
      CALL finish(method_name, "invalid arguments")

    var_start(:) = (/1, 1/)
    var_end(:) = var_dimlen(:)

    IF (PRESENT(start_extdim)) THEN
      var_start(2) = start_extdim
      var_end(2) = end_extdim
    END IF

    IF (PRESENT(extdim_name)) THEN
      CALL check_dimensions(stream_id%file_id, variable_name_, 2, var_dimlen, &
        &                   location, (/extdim_name/), &
        &                   ref_var_dim_start=var_start, &
        &                   ref_var_dim_end=var_end)
    ELSE
      CALL check_dimensions(stream_id%file_id, variable_name_, 2, var_dimlen, &
        &                   location, ref_var_dim_start=var_start, &
        &                   ref_var_dim_end=var_end)
    END IF

    IF (PRESENT(has_missValue) .AND. PRESENT(missValue)) THEN
      CALL netcdf_get_missValue(stream_id%file_id, variable_name_, has_missValue, missValue)
    ENDIF
    
    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      tmp_pointer => &
         & netcdf_read_2D_extdim(stream_id%file_id, variable_name_, fill_array,&
         & stream_id%read_info(location, 1)%n_g, &
         & stream_id%read_info(location, 1)%scatter_pattern, start_extdim, &
         & end_extdim, extdim_name )
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
    CASE (read_netcdf_distribute_method)
      IF (PRESENT(fill_array)) THEN
        tmp_pointer => fill_array
      ELSE
        IF (PRESENT(start_extdim)) THEN
          var_dimlen(2) = end_extdim - start_extdim + 1
        ELSE
          CALL distrib_inq_var_dims(stream_id%file_id, variable_name_, &
            &                       var_ndims, var_dimlen)
          var_end(2) = var_dimlen(2)
        END IF
        ALLOCATE(tmp_pointer(nproma, &
          (stream_id%read_info(location, 1)%n_l - 1)/nproma + 1, var_dimlen(2)))
        tmp_pointer(:,:,:) = 0.0_wp
      ENDIF
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
      tmp_ptr(1)%p => tmp_pointer
      CALL distrib_read(stream_id%file_id, variable_name_, tmp_ptr, &
        &               (/stream_id%read_info(location, 1)%dist_read_info/), &
        &               edim=var_dimlen(2:2), dimo=idx_blk_time, &
        &               start_ext_dim=var_start(2:2), end_ext_dim=var_end(2:2))
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END SUBROUTINE read_dist_REAL_2D_extdim
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, n) fortran-style: O3(n, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, blocks, time)
  SUBROUTINE read_dist_REAL_2D_extdim_multivar( &
    stream_id, location, variable_name, n_var, fill_array, return_pointer, &
    start_extdim, end_extdim, extdim_name )

    TYPE(t_stream_id), INTENT(INOUT)       :: stream_id
    INTEGER, INTENT(IN)                    :: location
    CHARACTER(LEN=*), INTENT(IN)           :: variable_name
    INTEGER, INTENT(IN)                    :: n_var
    TYPE(t_ptr_3d), OPTIONAL         :: fill_array(:)
    TYPE(t_ptr_3d), OPTIONAL         :: return_pointer(:)
    INTEGER, INTENT(in), OPTIONAL          :: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name

    TYPE(t_ptr_3d)                   :: tmp_return(n_var)
    TYPE(t_p_scatterPattern)               :: scatter_patterns(n_var)
    INTEGER                                :: n_g, i
    TYPE(t_ptr_3d), ALLOCATABLE      :: var_data_3d(:)
    INTEGER                                :: var_dimlen(2), var_ndims, &
      &                                       var_start(2), var_end(2)
    CHARACTER(LEN=NF90_MAX_NAME)             :: variable_name_
    CHARACTER(LEN=*), PARAMETER            :: method_name = &
      'mo_read_interface:read_dist_REAL_2D_extdim_multivar'
    TYPE(t_distrib_read_data) :: io_data(n_var)

    ! make variable name available on all processes
    CALL bcast_varname(variable_name, variable_name_)

    ! check whether fill_array and/or return_pointer was provided
    IF (.NOT. (PRESENT(fill_array) .OR. PRESENT(return_pointer))) &
      CALL finish(method_name, "invalid arguments")

    IF (PRESENT(start_extdim) .NEQV. PRESENT(end_extdim)) &
      CALL finish(method_name, "invalid arguments")

    IF (SIZE(stream_id%read_info, 2) /= n_var) &
      CALL finish(method_name, &
        "number of stream ids does not match number of variables")

    var_dimlen(:) = (/stream_id%read_info(location, 1)%n_g, -1/)
    IF (PRESENT(fill_array)) THEN
      var_dimlen(2) = SIZE(fill_array(1)%p, 3)
    END IF

    var_start(:) = (/1, 1/)
    var_end(:) = var_dimlen(:)

    IF (PRESENT(start_extdim)) THEN
      var_start(2) = start_extdim
      var_end(2) = end_extdim
    END IF

    IF (PRESENT(extdim_name)) THEN
      CALL check_dimensions(stream_id%file_id, variable_name_, 2, var_dimlen, &
        &                   location, (/extdim_name/), &
        &                   ref_var_dim_start=var_start, &
        &                   ref_var_dim_end=var_end)
    ELSE
      CALL check_dimensions(stream_id%file_id, variable_name_, 2, var_dimlen, &
        &                   location, ref_var_dim_start=var_start, &
        &                   ref_var_dim_end=var_end)
    END IF
    
    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      DO i = 1, n_var
        scatter_patterns(i)%p => &
          stream_id%read_info(location, i)%scatter_pattern
      END DO
      n_g = stream_id%read_info(location, 1)%n_g
      tmp_return = netcdf_read_2D_extdim( &
        stream_id%file_id, variable_name_, n_var, fill_array, n_g, &
        scatter_patterns, start_extdim, end_extdim, extdim_name )
      IF (PRESENT(return_pointer)) return_pointer(1:n_var) = tmp_return
    CASE (read_netcdf_distribute_method)

      ALLOCATE(var_data_3d(n_var))

      IF (PRESENT(start_extdim)) THEN
        var_dimlen(2) = end_extdim - start_extdim + 1
      ELSE
        IF (.NOT. PRESENT(fill_array)) &
          CALL distrib_inq_var_dims(stream_id%file_id, variable_name_, &
            &                       var_ndims, var_dimlen)
         var_end(2) = var_dimlen(2)
      END IF

      ! gather pointers of all output fields
      IF (PRESENT(fill_array)) THEN
        DO i = 1, n_var
          var_data_3d(i)%p => fill_array(i)%p
        END DO
      ELSE
        DO i = 1, n_var
          ALLOCATE(var_data_3d(i)%p(nproma, &
            (stream_id%read_info(location, 1)%n_l - 1)/nproma + 1, var_dimlen(2)))
          var_data_3d(i)%p(:,:,:) = 0.0_wp
        END DO
      ENDIF
      IF (PRESENT(return_pointer)) THEN
        DO i = 1, n_var
          return_pointer(i)%p => var_data_3d(i)%p
        END DO
      END IF
      DO i = 1, n_var
        io_data(i)%basic_data_index &
             = stream_id%read_info(location, i)%dist_read_info%basic_data_index
        io_data(i)%pat => stream_id%read_info(location, i)%dist_read_info%pat
      END DO

      CALL distrib_read(stream_id%file_id, variable_name_, var_data_3d, &
        &               io_data, edim=var_dimlen(2:2), dimo=idx_blk_time, &
        &               start_ext_dim=var_start(2:2), end_ext_dim=var_end(2:2))
      DEALLOCATE(var_data_3d)
    CASE DEFAULT
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END SUBROUTINE read_dist_REAL_2D_extdim_multivar
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, n) fortran-style: O3(n, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, blocks, time)
  SUBROUTINE read_dist_INT_2D_extdim(stream_id, location, variable_name, &
    &                                fill_array, return_pointer, start_extdim, &
    &                                end_extdim, extdim_name )

    TYPE(t_stream_id), INTENT(INOUT)       :: stream_id
    INTEGER, INTENT(IN)                    :: location
    CHARACTER(LEN=*), INTENT(IN)           :: variable_name
    define_fill_target_int                 :: fill_array(:,:,:)
    define_return_pointer_int              :: return_pointer(:,:,:)
    INTEGER, INTENT(in), OPTIONAL          :: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name
    TYPE(t_ptr_3d_int)                     :: tmp_ptr(1)  
    INTEGER, POINTER                       :: tmp_pointer(:,:,:)
    INTEGER                                :: var_dimlen(2), var_start(2), &
      &                                       var_end(2), var_ndims
    CHARACTER(LEN=NF90_MAX_NAME)             :: variable_name_
    CHARACTER(LEN=*), PARAMETER            :: method_name = &
      'mo_read_interface:read_dist_INT_2D_extdim'

    ! make variable name available on all processes
    CALL bcast_varname(variable_name, variable_name_)

    var_dimlen(:) = (/stream_id%read_info(location, 1)%n_g, -1/)
    IF (PRESENT(fill_array)) THEN
      var_dimlen(2) = SIZE(fill_array, 3)
    END IF

    ! check whether fill_array and/or return_pointer was provided
    IF (.NOT. (PRESENT(fill_array) .OR. PRESENT(return_pointer))) &
      CALL finish(method_name, "invalid arguments")

    IF (PRESENT(start_extdim) .NEQV. PRESENT(end_extdim)) &
      CALL finish(method_name, "invalid arguments")

    var_start(:) = (/1, 1/)
    var_end(:) = var_dimlen(:)

    IF (PRESENT(start_extdim)) THEN
      var_start(2) = start_extdim
      var_end(2) = end_extdim
    END IF

    IF (PRESENT(extdim_name)) THEN
      CALL check_dimensions(stream_id%file_id, variable_name_, 2, var_dimlen, &
        &                   location, (/extdim_name/), &
        &                   ref_var_dim_start=var_start, &
        &                   ref_var_dim_end=var_end)
    ELSE
      CALL check_dimensions(stream_id%file_id, variable_name_, 2, var_dimlen, &
        &                   location, ref_var_dim_start=var_start, &
        &                   ref_var_dim_end=var_end)
    END IF

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      tmp_pointer => &
         & netcdf_read_2D_extdim_int(stream_id%file_id, variable_name_, &
         & fill_array, stream_id%read_info(location, 1)%n_g, &
         & stream_id%read_info(location, 1)%scatter_pattern, start_extdim, &
         & end_extdim, extdim_name )
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
    CASE (read_netcdf_distribute_method)
      IF (PRESENT(fill_array)) THEN
        tmp_pointer => fill_array
      ELSE
        IF (PRESENT(start_extdim)) THEN
          var_dimlen(2) = end_extdim - start_extdim + 1
        ELSE
          CALL distrib_inq_var_dims(stream_id%file_id, variable_name_, &
            &                       var_ndims, var_dimlen)
          var_end(2) = var_dimlen(2)
        END IF
        ALLOCATE(tmp_pointer(nproma, &
          (stream_id%read_info(location, 1)%n_l - 1)/nproma + 1, var_dimlen(2)))
        tmp_pointer(:,:,:) = 0
      ENDIF
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
      tmp_ptr(1)%p => tmp_pointer
      CALL distrib_read(stream_id%file_id, variable_name_, tmp_ptr, &
        &               (/stream_id%read_info(location, 1)%dist_read_info/), &
        &               edim=var_dimlen(2:2), dimo=idx_blk_time, &
        &               start_ext_dim=var_start(2:2), end_ext_dim=var_end(2:2))
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END SUBROUTINE read_dist_INT_2D_extdim
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, n) fortran-style: O3(n, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, blocks, time)
  SUBROUTINE read_dist_INT_2D_extdim_multivar( &
    stream_id, location, variable_name, n_var, fill_array, return_pointer, &
    start_extdim, end_extdim, extdim_name )

    TYPE(t_stream_id), INTENT(INOUT)       :: stream_id
    INTEGER, INTENT(IN)                    :: location
    CHARACTER(LEN=*), INTENT(IN)           :: variable_name
    INTEGER, INTENT(IN)                    :: n_var
    TYPE(t_ptr_3d_int), OPTIONAL        :: fill_array(:)
    TYPE(t_ptr_3d_int), OPTIONAL        :: return_pointer(:)
    INTEGER, INTENT(in), OPTIONAL          :: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name

    TYPE(t_ptr_3d_int)                  :: tmp_return(n_var)
    TYPE(t_p_scatterPattern)               :: scatter_patterns(n_var)
    INTEGER                                :: n_g, i
    TYPE(t_ptr_3d_int), ALLOCATABLE     :: var_data_3d(:)
    INTEGER                                :: var_dimlen(2), var_ndims, &
      &                                       var_start(2), var_end(2)
    CHARACTER(LEN=NF90_MAX_NAME)             :: variable_name_
    CHARACTER(LEN=*), PARAMETER            :: method_name = &
      'mo_read_interface:read_dist_INT_2D_extdim_multivar'
    TYPE(t_distrib_read_data) :: io_data(n_var)

    ! make variable name available on all processes
    CALL bcast_varname(variable_name, variable_name_)

    ! check whether fill_array and/or return_pointer was provided
    IF (.NOT. (PRESENT(fill_array) .OR. PRESENT(return_pointer))) &
      CALL finish(method_name, "invalid arguments")

    IF (PRESENT(start_extdim) .NEQV. PRESENT(end_extdim)) &
      CALL finish(method_name, "invalid arguments")

    IF (SIZE(stream_id%read_info, 2) /= n_var) &
      CALL finish(method_name, &
        "number of stream ids does not match number of variables")

    var_dimlen(:) = (/stream_id%read_info(location, 1)%n_g, -1/)
    IF (PRESENT(fill_array)) THEN
      var_dimlen(2) = SIZE(fill_array(1)%p, 3)
    END IF

    var_start(:) = (/1, 1/)
    var_end(:) = var_dimlen(:)

    IF (PRESENT(start_extdim)) THEN
      var_start(2) = start_extdim
      var_end(2) = end_extdim
    END IF

    IF (PRESENT(extdim_name)) THEN
      CALL check_dimensions(stream_id%file_id, variable_name_, 2, var_dimlen, &
        &                   location, (/extdim_name/), &
        &                   ref_var_dim_start=var_start, &
        &                   ref_var_dim_end=var_end)
    ELSE
      CALL check_dimensions(stream_id%file_id, variable_name_, 2, var_dimlen, &
        &                   location, ref_var_dim_start=var_start, &
        &                   ref_var_dim_end=var_end)
    END IF

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      DO i = 1, n_var
        scatter_patterns(i)%p => &
          stream_id%read_info(location, i)%scatter_pattern
      END DO
      n_g = stream_id%read_info(location, 1)%n_g
      tmp_return = netcdf_read_2D_extdim_int( &
        stream_id%file_id, variable_name_, n_var, fill_array, n_g, &
        scatter_patterns, start_extdim, end_extdim, extdim_name )
      IF (PRESENT(return_pointer)) return_pointer(1:n_var) = tmp_return
    CASE (read_netcdf_distribute_method)

      ALLOCATE(var_data_3d(n_var))

      IF (PRESENT(start_extdim)) THEN
        var_dimlen(2) = end_extdim - start_extdim + 1
      ELSE
        IF (.NOT. PRESENT(fill_array)) &
          CALL distrib_inq_var_dims(stream_id%file_id, variable_name_, &
            &                       var_ndims, var_dimlen)
        var_end(2) = var_dimlen(2)
      END IF

      ! gather pointers of all output fields
      IF (PRESENT(fill_array)) THEN
        DO i = 1, n_var
          var_data_3d(i)%p => fill_array(i)%p
        END DO
      ELSE
        DO i = 1, n_var
          ALLOCATE(var_data_3d(i)%p(nproma, &
            (stream_id%read_info(location, 1)%n_l - 1)/nproma + 1, var_dimlen(2)))
          var_data_3d(i)%p(:,:,:) = 0
        END DO
      ENDIF
      IF (PRESENT(return_pointer)) THEN
        DO i = 1, n_var
          return_pointer(i)%p => var_data_3d(i)%p
        END DO
      END IF
      DO i = 1, n_var
        io_data(i)%basic_data_index &
             = stream_id%read_info(location, i)%dist_read_info%basic_data_index
        io_data(i)%pat => stream_id%read_info(location, i)%dist_read_info%pat
      END DO
      CALL distrib_read(stream_id%file_id, variable_name_, var_data_3d, &
        &               io_data, edim=var_dimlen(2:2), dimo=idx_blk_time, &
        &               start_ext_dim=var_start(2:2), end_ext_dim=var_end(2:2))
      DEALLOCATE(var_data_3d)
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END SUBROUTINE read_dist_INT_2D_extdim_multivar
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_bcast_REAL_3D(file_id, variable_name, fill_array, &
    &                           return_pointer)

    INTEGER, INTENT(IN)          :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    define_fill_target           :: fill_array(:,:,:)
    define_return_pointer        :: return_pointer(:,:,:)

    REAL(wp), POINTER            :: tmp_pointer(:,:,:)
    CHARACTER(LEN=NF90_MAX_NAME)   :: variable_name_
    CHARACTER(LEN=*), PARAMETER  :: method_name = &
      'mo_read_interface:read_bcast_REAL_3D'

    ! make variable name available on all processes
    CALL bcast_varname(variable_name, variable_name_)

    ! check whether fill_array and/or return_pointer was provided
    IF (.NOT. (PRESENT(fill_array) .OR. PRESENT(return_pointer))) &
      CALL finish(method_name, "invalid arguments")

    ! there is only one implementation for this read routine type
    tmp_pointer => netcdf_read_REAL_3D_all(file_id, variable_name_, fill_array)
    IF (PRESENT(return_pointer)) return_pointer => tmp_pointer

  END SUBROUTINE read_bcast_REAL_3D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O2(levels, n) fortran-style: O2(n, levels)
  ! The fill_array  has the structure:
  !       fill_array(nproma, levels, blocks)
  SUBROUTINE read_dist_REAL_3D(stream_id, location, variable_name, fill_array, &
    &                          return_pointer, levelsDimName)

    TYPE(t_stream_id), INTENT(INOUT)       :: stream_id
    INTEGER, INTENT(IN)                    :: location
    CHARACTER(LEN=*), INTENT(IN)           :: variable_name
    define_fill_target                     :: fill_array(:,:,:)
    define_return_pointer                  :: return_pointer(:,:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: levelsDimName
    TYPE(t_ptr_3d)                         :: tmp_ptr(1)
    INTEGER                                :: var_ndims, var_dimlen(2)
    REAL(wp), POINTER                      :: tmp_pointer(:,:,:)
    CHARACTER(LEN=NF90_MAX_NAME)             :: variable_name_
    CHARACTER(LEN=*), PARAMETER            :: method_name = &
      'mo_read_interface:read_dist_REAL_3D'

    ! make variable name available on all processes
    CALL bcast_varname(variable_name, variable_name_)

    var_dimlen(:) = (/stream_id%read_info(location, 1)%n_g, -1/)
    IF (PRESENT(fill_array)) THEN
      var_dimlen(2) = SIZE(fill_array, 2)
    END IF

    ! check whether fill_array and/or return_pointer was provided
    IF (.NOT. (PRESENT(fill_array) .OR. PRESENT(return_pointer))) &
      CALL finish(method_name, "invalid arguments")

    IF (PRESENT(levelsDimName)) THEN
      CALL check_dimensions(stream_id%file_id, variable_name_, 2, var_dimlen, &
        &                   location, (/levelsDimName/))
    ELSE
      CALL check_dimensions(stream_id%file_id, variable_name_, 2, var_dimlen, &
        &                   location)
    END IF

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      tmp_pointer => &
         & netcdf_read_3D(stream_id%file_id, variable_name_, &
         & fill_array, stream_id%read_info(location, 1)%n_g, &
         & stream_id%read_info(location, 1)%scatter_pattern, levelsDimName)
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
    CASE (read_netcdf_distribute_method)
      IF (PRESENT(fill_array)) THEN
        tmp_pointer => fill_array
      ELSE
        CALL distrib_inq_var_dims(stream_id%file_id, variable_name_, &
          &                       var_ndims, var_dimlen)
        ALLOCATE(tmp_pointer(nproma, var_dimlen(2), &
          (stream_id%read_info(location, 1)%n_l - 1)/nproma + 1))
        tmp_pointer(:,:,:) = 0.0_wp
      ENDIF
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer
      tmp_ptr(1)%p => tmp_pointer
      CALL distrib_read(stream_id%file_id, variable_name_, tmp_ptr, &
        &               (/stream_id%read_info(location, 1)%dist_read_info/), &
        &               edim=var_dimlen(2:2), dimo=idx_lvl_blk)
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END SUBROUTINE read_dist_REAL_3D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, levels, n) fortran-style: O3(n, levels, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, levels, blocks)
  ! Since the array in the file has a time dimension we can map this case to
  ! read_dist_REAL_3D_extdim. However fill_array lacks the time dimension. To
  ! add this dimension we use an assumed-size array in read_dist_REAL_3D_1time_.
  ! In order to use assumed-size in this case we need the shape of the original
  ! fill_array. This is determined by read_dist_REAL_3D_1time.
  SUBROUTINE read_dist_REAL_3D_1time(stream_id, location, variable_name, &
    & fill_array, return_pointer, levelsDimName,                         &
    & has_missValue, missValue)

    TYPE(t_stream_id), INTENT(INOUT)       :: stream_id
    INTEGER, INTENT(IN)                    :: location
    CHARACTER(LEN=*), INTENT(IN)           :: variable_name
    define_fill_target                     :: fill_array(:,:,:)
    define_return_pointer                  :: return_pointer(:,:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: levelsDimName
    LOGICAL, OPTIONAL                      :: has_missValue
    REAL(wp), OPTIONAL                     :: missValue
    
    IF (PRESENT(fill_array)) THEN
      CALL read_dist_REAL_3D_1time_(stream_id, location, variable_name, &
        & SHAPE(fill_array), fill_array, &
        & return_pointer, levelsDimName, &
        & has_missValue, missValue)
    ELSE
      CALL read_dist_REAL_3D_1time_(stream_id, location, variable_name, &
        & (/0,0,0/), return_pointer=return_pointer, &
        & levelsDimName=levelsDimName,              &
        & has_missValue=has_missValue,              &
        & missValue=missValue)
    END IF

  END SUBROUTINE read_dist_REAL_3D_1time

  SUBROUTINE read_dist_REAL_3D_1time_(stream_id, location, variable_name, &
    & array_shape, fill_array, return_pointer, &
    & levelsDimName,                           &
    & has_missValue, missValue)

    TYPE(t_stream_id), INTENT(INOUT)       :: stream_id
    INTEGER, INTENT(IN)                    :: location
    CHARACTER(LEN=*), INTENT(IN)           :: variable_name
    INTEGER, INTENT(IN)                    :: array_shape(3)
    define_fill_target                     :: fill_array(array_shape(1), &
      &                                                  array_shape(2), &
      &                                                  array_shape(3), 1)
    define_return_pointer                  :: return_pointer(:,:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: levelsDimName
    LOGICAL, OPTIONAL                      :: has_missValue
    REAL(wp), OPTIONAL                     :: missValue
    
    CHARACTER(LEN=*), PARAMETER            :: method_name = &
      'mo_read_interface:read_dist_REAL_3D_1time_'

    REAL(wp), POINTER :: return_pointer_(:,:,:,:)

    ! Since fill_array now has a time dimension we can call
    ! read_dist_REAL_3D_extdim
    IF (PRESENT(return_pointer)) THEN

      CALL read_dist_REAL_3D_extdim(  &
        & stream_id=stream_id,                 &
        & location=location,                   &
        & variable_name=variable_name,         &
        & fill_array=fill_array,               &
        & return_pointer=return_pointer_,      &
        & start_extdim=1,                      &
        & end_extdim=1,                        &
        & levelsDimName=levelsDimName,         &
        & extdim_name="time",                  &
        & has_missValue=has_missValue,         &
        & missValue=missValue)

      return_pointer => return_pointer_(:,:,:,1)
      ALLOCATE(return_pointer(SIZE(return_pointer_,1), &
        &                     SIZE(return_pointer_,2), &
        &                     SIZE(return_pointer_,3)))
      return_pointer(:,:,:) = return_pointer_(:,:,:,1)
      DEALLOCATE(return_pointer_)
    ELSE

      CALL read_dist_REAL_3D_extdim(  &
        & stream_id=stream_id,                 &
        & location=location,                   &
        & variable_name=variable_name,         &
        & fill_array=fill_array,               &
        & start_extdim=1,                      &
        & end_extdim=1,                        &
        & levelsDimName=levelsDimName,         &
        & extdim_name="time",                  &
        & has_missValue=has_missValue,         &
        & missValue=missValue)
    END IF

  END SUBROUTINE read_dist_REAL_3D_1time_
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, levels, n) fortran-style: O3(n, levels, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, levels, blocks, time)
  ! We can map this case to read_dist_REAL_3D_extdim.
  SUBROUTINE read_dist_REAL_3D_time(stream_id, location, variable_name, &
    & fill_array, return_pointer, start_timestep,  &
    & end_timestep, levelsDimName,                 &
    & has_missValue, missValue)

    TYPE(t_stream_id), INTENT(INOUT)       :: stream_id
    INTEGER, INTENT(IN)                    :: location
    CHARACTER(LEN=*), INTENT(IN)           :: variable_name
    define_fill_target                     :: fill_array(:,:,:,:)
    define_return_pointer                  :: return_pointer(:,:,:,:)
    INTEGER, INTENT(in), OPTIONAL          :: start_timestep, end_timestep
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: levelsDimName
    LOGICAL, OPTIONAL                      :: has_missValue
    REAL(wp), OPTIONAL                     :: missValue
    
    CHARACTER(LEN=*), PARAMETER            :: method_name = &
      'mo_read_interface:read_dist_REAL_3D_time'

    CALL read_dist_REAL_3D_extdim( &
      & stream_id=stream_id,                 &
      & location=location,                   &
      & variable_name=variable_name,         &
      & fill_array=fill_array,               &
      & return_pointer=return_pointer,       &
      & start_extdim=start_timestep,         &
      & end_extdim=end_timestep,             &
      & levelsDimName=levelsDimName,         &
      & extdim_name="time",                  &
      & has_missValue=has_missValue,         &
      & missValue=missValue)

  END SUBROUTINE read_dist_REAL_3D_time
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  ! By default the netcdf input has the structure :
  !      c-style(ncdump): O3(time, levels, n) fortran-style: O3(n, levels, time)
  ! The fill_array  has the structure:
  !       fill_array(nproma, levels, blocks, time)
  SUBROUTINE read_dist_REAL_3D_extdim(stream_id, location, variable_name, &
    & fill_array, return_pointer, start_extdim,&
    & end_extdim, levelsDimName, extdim_name,  &
    & has_missValue, missValue)

    TYPE(t_stream_id), INTENT(INOUT)       :: stream_id
    INTEGER, INTENT(IN)                    :: location
    CHARACTER(LEN=*), INTENT(IN)           :: variable_name
    define_fill_target                     :: fill_array(:,:,:,:)
    define_return_pointer                  :: return_pointer(:,:,:,:)
    INTEGER, INTENT(in), OPTIONAL          :: start_extdim, end_extdim
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name, levelsDimName
    LOGICAL, OPTIONAL                      :: has_missValue
    REAL(wp), OPTIONAL                     :: missValue

    INTEGER                                :: var_ndims, var_dimlen(3), &
      &                                       var_start(3), var_end(3)
    TYPE(t_ptr_4d)                         :: tmp_pointer(1)
    CHARACTER(LEN=128)                     :: temp_string_array(2)
    CHARACTER(LEN=NF90_MAX_NAME)             :: variable_name_
    CHARACTER(LEN=*), PARAMETER            :: method_name = &
      'mo_read_interface:read_dist_REAL_3D_extdim'

    ! make variable name available on all processes
    CALL bcast_varname(variable_name, variable_name_)

    var_dimlen(:) = (/stream_id%read_info(location, 1)%n_g, -1, -1/)
    IF (PRESENT(fill_array)) THEN
      var_dimlen(2) = SIZE(fill_array, 2)
      var_dimlen(3) = SIZE(fill_array, 4)
    END IF

    ! check whether fill_array and/or return_pointer was provided
    IF (.NOT. (PRESENT(fill_array) .OR. PRESENT(return_pointer))) &
      CALL finish(method_name, "invalid arguments")

    IF (PRESENT(start_extdim) .NEQV. PRESENT(end_extdim)) &
      CALL finish(method_name, "invalid arguments")

    var_start(:) = (/1, 1, 1/)
    var_end(:) = var_dimlen(:)

    IF (PRESENT(start_extdim)) THEN
      var_start(3) = start_extdim
      var_end(3) = end_extdim
    END IF

    IF (PRESENT(levelsDimName) .AND. PRESENT(extdim_name)) THEN
      temp_string_array(1) = levelsDimName
      temp_string_array(2) = extdim_name
      CALL check_dimensions(stream_id%file_id, variable_name_, 3, var_dimlen, &
        &                   location, temp_string_array, &
        &                   ref_var_dim_start=var_start, &
        &                   ref_var_dim_end=var_end)
    ELSE
      CALL check_dimensions(stream_id%file_id, variable_name_, 3, var_dimlen, &
        &                   location, ref_var_dim_start=var_start, &
        &                   ref_var_dim_end=var_end)
    END IF

    IF (PRESENT(has_missValue) .AND. PRESENT(missValue)) THEN
      CALL netcdf_get_missValue(stream_id%file_id, variable_name_, has_missValue, missValue)
    ENDIF

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      tmp_pointer(1)%p => &
         & netcdf_read_3D_extdim(stream_id%file_id, variable_name_, &
         & fill_array, stream_id%read_info(location, 1)%n_g, &
         & stream_id%read_info(location, 1)%scatter_pattern, &
         & start_extdim, end_extdim, levelsDimName, extdim_name )
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer(1)%p
    CASE (read_netcdf_distribute_method)
      IF (PRESENT(fill_array)) THEN
        tmp_pointer(1)%p => fill_array
      ELSE
        CALL distrib_inq_var_dims(stream_id%file_id, variable_name_, &
          &                       var_ndims, var_dimlen)
        IF (PRESENT(start_extdim)) var_dimlen(3) = end_extdim - start_extdim + 1
        ALLOCATE(tmp_pointer(1)%p(nproma, var_dimlen(2), &
          (stream_id%read_info(location, 1)%n_l - 1)/nproma + 1, var_dimlen(3)))
        tmp_pointer(1)%p(:,:,:,:) = 0.0_wp
      ENDIF
      var_end(:) = var_start(:) + var_dimlen(:) - 1
      IF (PRESENT(return_pointer)) return_pointer => tmp_pointer(1)%p

      CALL distrib_read(stream_id%file_id, variable_name_, tmp_pointer, &
        &               [stream_id%read_info(location, 1)%dist_read_info], &
        &               edim=var_dimlen(2:3), start_ext_dim=var_start(2:3), &
        &               end_ext_dim=var_end(2:3))
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT


  END SUBROUTINE read_dist_REAL_3D_extdim
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE openInputFile_dist_multivar(stream_id, filename, &
       patches, input_method)
    TYPE(t_stream_id), INTENT(out) :: stream_id
    CHARACTER(LEN=*), INTENT(IN) :: filename
    TYPE(p_t_patch), TARGET, INTENT(IN) :: patches(:)
    INTEGER, OPTIONAL, INTENT(IN) :: input_method

    INTEGER :: n_var, i

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:openInputFile_dist_multivar'

    n_var = SIZE(patches)

    IF (n_var < 1) CALL finish(method_name, "invalid number of patches")

    IF (PRESENT(input_method)) THEN
      stream_id%input_method = input_method
    ELSE
      stream_id%input_method = default_read_method
    END IF

    ALLOCATE(stream_id%read_info(3, n_var))

    DO i = 1, n_var

      IF ((patches(1)%p%n_patch_cells_g /= patches(i)%p%n_patch_cells_g) .OR. &
          (patches(1)%p%n_patch_edges_g /= patches(i)%p%n_patch_edges_g) .OR. &
          (patches(1)%p%n_patch_verts_g /= patches(i)%p%n_patch_verts_g))     &
        CALL finish(method_name, "patches do not match")

      stream_id%read_info(on_cells, i)%n_g = &
        patches(i)%p%n_patch_cells_g
      stream_id%read_info(on_edges, i)%n_g = &
        patches(i)%p%n_patch_edges_g
      stream_id%read_info(on_vertices, i)%n_g = &
        patches(i)%p%n_patch_verts_g

      stream_id%read_info(on_cells, i)%n_l = &
        patches(i)%p%n_patch_cells
      stream_id%read_info(on_edges, i)%n_l = &
        patches(i)%p%n_patch_edges
      stream_id%read_info(on_vertices, i)%n_l = &
        patches(i)%p%n_patch_verts
    END DO

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)

      stream_id%file_id = netcdf_open_input(filename)

      DO i = 1, n_var
        stream_id%read_info(on_cells, i)%scatter_pattern => &
          patches(i)%p%comm_pat_scatter_c
        NULLIFY(stream_id%read_info(on_cells, i)%dist_read_info)
        stream_id%read_info(on_edges, i)%scatter_pattern => &
          patches(i)%p%comm_pat_scatter_e
        NULLIFY(stream_id%read_info(on_edges, i)%dist_read_info)
        stream_id%read_info(on_vertices, i)%scatter_pattern => &
          patches(i)%p%comm_pat_scatter_v
        NULLIFY(stream_id%read_info(on_vertices, i)%dist_read_info)
      END DO

    CASE (read_netcdf_distribute_method)

      stream_id%file_id = distrib_nf_open(TRIM(filename))

      DO i = 1, n_var
        stream_id%read_info(on_cells, i)%dist_read_info => &
          patches(i)%p%cells%dist_io_data
        NULLIFY(stream_id%read_info(on_cells, i)%scatter_pattern)
        stream_id%read_info(on_vertices, i)%dist_read_info => &
          patches(i)%p%verts%dist_io_data
        NULLIFY(stream_id%read_info(on_vertices, i)%scatter_pattern)
        stream_id%read_info(on_edges, i)%dist_read_info => &
          patches(i)%p%edges%dist_io_data
        NULLIFY(stream_id%read_info(on_edges, i)%scatter_pattern)
      END DO

    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

  END SUBROUTINE openInputFile_dist_multivar

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE openInputFile_dist(stream_id, filename, patch, input_method)
    TYPE(t_stream_id), INTENT(out) :: stream_id
    CHARACTER(LEN=*), INTENT(IN) :: filename
    TYPE(t_patch), TARGET, INTENT(IN) :: patch
    INTEGER, OPTIONAL, INTENT(IN) :: input_method

    TYPE(p_t_patch) :: patch_(1)

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:openInputFile_dist'

    patch_(1)%p => patch

    CALL openInputFile_dist_multivar(stream_id, filename, patch_, input_method)

  END SUBROUTINE openInputFile_dist
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !>
  SUBROUTINE openInputFile_bcast(stream_id, filename)

    INTEGER, INTENT(out) :: stream_id
    CHARACTER(LEN=*), INTENT(IN) :: filename

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:openInputFile_bcast'

    stream_id = netcdf_open_input(filename)

  END SUBROUTINE openInputFile_bcast
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE closeFile_dist(stream_id)
    TYPE(t_stream_id), INTENT(INOUT) :: stream_id

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:closeFile_dist'

    SELECT CASE(stream_id%input_method)
    CASE (read_netcdf_broadcast_method)
      stream_id%return_status = netcdf_close(stream_id%file_id)
    CASE (read_netcdf_distribute_method)
      CALL distrib_nf_close(stream_id%file_id)
    CASE default
      CALL finish(method_name, "unknown input_method")
    END SELECT

    DEALLOCATE(stream_id%read_info)

  END SUBROUTINE closeFile_dist
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE closeFile_bcast(file_id, return_status)
    INTEGER, INTENT(IN) :: file_id
    INTEGER, OPTIONAL, INTENT(OUT) :: return_status

    INTEGER :: ret
    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:closeFile_bcast'

    ret = netcdf_close(file_id)
    IF (PRESENT(return_status)) return_status = ret

  END SUBROUTINE closeFile_bcast

  SUBROUTINE check_dimensions(file_id, variable_name, ref_var_ndims, &
                              ref_var_dimlen, location, extdim_name, &
                              ref_var_dim_start, ref_var_dim_end)

    INTEGER, INTENT(IN)                    :: file_id
    CHARACTER(LEN=NF90_MAX_NAME), INTENT(IN) :: variable_name
    INTEGER, INTENT(IN)                    :: ref_var_ndims, &
      &                                       ref_var_dimlen(ref_var_ndims)
    INTEGER, OPTIONAL, INTENT(IN)          :: ref_var_dim_start(ref_var_ndims)
    INTEGER, OPTIONAL, INTENT(IN)          :: ref_var_dim_end(ref_var_ndims)
    INTEGER, INTENT(IN)                    :: location
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extdim_name(2:ref_var_ndims)

    INTEGER :: varid, var_ndims, var_dimlen(NF90_MAX_VAR_DIMS), &
      &        var_dimids(NF90_MAX_VAR_DIMS)
    CHARACTER(LEN=NF90_MAX_NAME) :: var_dim_name(NF90_MAX_VAR_DIMS)
    INTEGER :: i, tlen

    CHARACTER(LEN=*), PARAMETER :: method_name = &
      'mo_read_interface:check_dimensions'

    IF (file_id == -1) RETURN

    IF (PRESENT(ref_var_dim_start) .NEQV. PRESENT(ref_var_dim_end)) &
      CALL finish(method_name, "invalid arguments")

    tlen = LEN_TRIM(variable_name)
    CALL nf(nf90_inq_varid(file_id, variable_name(1:tlen), varid), &
      &     method_name // "("//variable_name(1:tlen)//")")
    CALL nf(nf90_inquire_variable(file_id, varid, ndims = var_ndims), &
      &     method_name // "("//variable_name(1:tlen)//")")
    CALL nf(nf90_inquire_variable(file_id, varid, dimids = var_dimids), method_name)
    DO i = 1, var_ndims
      CALL nf(nf90_inquire_dimension(file_id, var_dimids(i), len = var_dimlen(i)), method_name)
      CALL nf(nf90_inquire_dimension(file_id, var_dimids(i), name = var_dim_name(i)), method_name)
    END DO

    IF (var_ndims /= ref_var_ndims ) THEN
      WRITE(0,*) variable_name(1:tlen), ": var_ndims = ", var_ndims
      CALL finish(method_name, "Dimensions mismatch")
    ENDIF

    IF (PRESENT(ref_var_dim_start)) THEN
      DO i = 1, var_ndims
        IF ((ref_var_dim_start(i) /= -1) .AND. &
          & (ref_var_dim_end(i) /= -1)) THEN
          IF (ref_var_dim_end(i) < ref_var_dim_start(i)) THEN
            WRITE(0,*) variable_name(1:tlen), ": ref_var_dim_start(:) = ", &
              &        ref_var_dim_start(:), "; ref_var_dim_end(:) = ", &
              &        ref_var_dim_end(:)
            CALL finish(method_name, "invalid start end")
          END IF
          IF ((ref_var_dim_start(i) == 0) .OR. &
            & (ref_var_dim_start(i) > var_dimlen(i))) THEN
            WRITE(0,*) variable_name(1:tlen), ": ref_var_dim_start(:) = ", &
              &        ref_var_dim_start(:), "; var_dimlen(:) = ", &
              &        var_dimlen(1:ref_var_ndims)
            CALL finish(method_name, "invalid start")
          END IF
          IF ((ref_var_dim_end(i) > var_dimlen(i))) THEN
            WRITE(0,*) variable_name(1:tlen), ": ref_var_dim_end(:) = ", &
              &        ref_var_dim_end(:), "; var_dimlen(:) = ", &
              &        var_dimlen(1:ref_var_ndims)
            CALL finish(method_name, "invalid end")
          END IF
        END IF
      END DO
    ELSE
      DO i = 1, var_ndims
        IF ((ref_var_dimlen(i) /= -1) .AND. &
          & (ref_var_dimlen(i) /= var_dimlen(i))) THEN
          WRITE(0,*) variable_name(1:tlen), ": ref_var_dimlen(:) = ", &
            &        ref_var_dimlen(:), "; var_dimlen(:) = ", &
            &        var_dimlen(1:ref_var_ndims)
          CALL finish(method_name, "Dimensions mismatch")
        END IF
      END DO
    END IF

    ! check if the dim have reasonable names

    SELECT CASE(location)
      CASE (on_cells)
        IF (.NOT. ((TRIM(var_dim_name(1)) == 'cell') .OR. &
          &        (TRIM(var_dim_name(1)) == 'ncells'))) THEN
          write(0,*) TRIM(var_dim_name(1))
          WRITE(message_text,*) variable_name(1:tlen), " ", &
            &                   TRIM(var_dim_name(1)), " /= std_cells_dim_name"
          CALL finish(method_name, message_text)
        ENDIF
      CASE (on_vertices)
        IF (.NOT. ((TRIM(var_dim_name(1)) == 'vertex') .OR. &
          &        (TRIM(var_dim_name(1)) == 'nverts') .OR. &
          &        (TRIM(var_dim_name(1)) == 'ncells_3'))) THEN
          write(0,*) TRIM(var_dim_name(1))
          WRITE(message_text,*) variable_name(1:tlen), " ", TRIM(var_dim_name(1)), &
            &                   " /= std_verts_dim_name"
          CALL finish(method_name, message_text)
        ENDIF
      CASE (on_edges)
        IF (.NOT. ((TRIM(var_dim_name(1)) == 'edge') .OR. &
          &        (TRIM(var_dim_name(1)) == 'nedges') .OR. &
          &        (TRIM(var_dim_name(1)) == 'ncells_2'))) THEN
          write(0,*) TRIM(var_dim_name(1))
          WRITE(message_text,*) variable_name(1:tlen), " ", TRIM(var_dim_name(1)), &
            &                   " /= std_edge_dim_name"
          CALL finish(method_name, message_text)
        ENDIF
    END SELECT

    IF (PRESENT(extdim_name)) THEN
      DO i = 2, ref_var_ndims
        IF (TRIM(extdim_name(i)) /= TRIM(var_dim_name(i))) THEN
          WRITE(message_text,*) variable_name(1:tlen), ":", &
            &                   TRIM(extdim_name(i)), "/=",  &
            &                   TRIM(var_dim_name(i))
          CALL finish(method_name, TRIM(message_text))
        ENDIF
      END DO
    END IF

  END SUBROUTINE check_dimensions

  !-------------------------------------------------------------------------

  SUBROUTINE bcast_varname(string_in, string_out)

    CHARACTER(LEN=*), INTENT(IN) :: string_in
    CHARACTER(LEN=NF90_MAX_NAME), INTENT(OUT) :: string_out

    IF (my_process_is_mpi_workroot()) THEN
      IF (LEN(string_in) > NF90_MAX_NAME) &
        CALL finish("bcast_varname", "invalid string length")

      string_out = string_in
    END IF

    CALL p_bcast(string_out, p_io, &
      &          MERGE(p_comm_work_test, p_comm_work, p_test_run))
  END SUBROUTINE bcast_varname

END MODULE mo_read_interface
