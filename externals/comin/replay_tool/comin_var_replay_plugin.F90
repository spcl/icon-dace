!> @file comin_var_replay_plugin.F90
!! @brief A plugin that replays the values of variables from a netcdf
!! file (recorded with `comin_var_recorder_plugin`).
!
!  @authors 01/2024 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
MODULE comin_var_replay_plugin

  USE comin_plugin_interface, ONLY: t_comin_var_ptr, t_comin_plugin_info,                          &
       &                            t_comin_descrdata_global, t_comin_var_descriptor,              &
       &                            t_comin_setup_version_info, comin_current_get_plugin_info,     &
       &                            comin_plugin_finish, comin_parallel_get_host_mpi_rank,         &
       &                            comin_setup_get_version, comin_descrdata_get_global,           &
       &                            comin_var_request_add, comin_callback_register, comin_var_get, &
       &                            comin_metadata_set,                                            &
       &                            EP_ATM_TIMELOOP_START, COMIN_FLAG_WRITE,                       &
       &                            EP_SECONDARY_CONSTRUCTOR,                                      &
       &                            comin_parallel_get_host_mpi_comm, comin_error_check,           &
       &                            COMIN_METADATA_TYPEID_INTEGER, COMIN_METADATA_TYPEID_REAL,     &
       &                            COMIN_METADATA_TYPEID_CHARACTER, COMIN_METADATA_TYPEID_LOGICAL
  USE netcdf,                 ONLY: nf90_open, nf90_inq_dimid, nf90_inq_ncid,                      &
       &                            nf90_inquire, nf90_inquire_variable, nf90_inq_varids,          &
       &                            nf90_get_att, nf90_get_var, NF90_GLOBAL, NF90_NOWRITE,         &
       &                            NF90_NOERR, nf90_inquire_attribute, NF90_DOUBLE, NF90_INT,     &
       &                            NF90_BYTE, NF90_CHAR, nf90_inq_attname
  USE netcdf_utils,           ONLY: nf90
  USE iso_c_binding,          ONLY: C_DOUBLE
  USE utils,                  ONLY: int2string

  IMPLICIT NONE

  TYPE :: var_ptr_ptr
    TYPE(t_comin_var_ptr), POINTER :: ptr
  END TYPE var_ptr_ptr

  ! entry point where the variable values are injected
  INTEGER, PARAMETER :: ep = EP_ATM_TIMELOOP_START

  INTEGER :: ncid, host_rank, host_comm, host_comm_size
  INTEGER, ALLOCATABLE :: domain_ncids(:)
  INTEGER :: event_dimid, event_counter = 1

  TYPE(t_comin_plugin_info)                 :: plugin_info
  TYPE(t_comin_descrdata_global), POINTER   :: global

  INTEGER,                      ALLOCATABLE :: nvars(:)
  TYPE(t_comin_var_descriptor), ALLOCATABLE :: var_descrs(:,:)
  TYPE(var_ptr_ptr),            ALLOCATABLE :: var_ptrs(:,:)
  INTEGER,                      ALLOCATABLE :: var_varids(:,:)

CONTAINS

  SUBROUTINE comin_main() BIND(C)

    INTEGER                          :: ierr, file_comin_version(3), status, i, j, jg
    INTEGER                          :: file_host_rank, file_host_size, tracer
    INTEGER                          :: nAtts, atttype, attlen, attvalue_int
    REAL(C_DOUBLE)                   :: attvalue_real
    TYPE(t_comin_setup_version_info) :: comin_version
    CHARACTER(LEN=64)                :: name, attname
    CHARACTER(LEN=:), ALLOCATABLE    :: filename_prefix, attvalue_char
    TYPE(t_comin_var_descriptor)     :: tmp_descr

    CALL comin_current_get_plugin_info(plugin_info)

    host_rank = comin_parallel_get_host_mpi_rank()
    host_comm = comin_parallel_get_host_mpi_comm()
    CALL MPI_COMM_SIZE(host_comm, host_comm_size, ierr)
    IF (LEN_TRIM(plugin_info%options) == 0) THEN
      filename_prefix = "vars_"
    ELSE
      filename_prefix = plugin_info%options
    END IF
    WRITE(0,*) "Reading variables from " // filename_prefix // int2string(host_rank) // ".nc"
    CALL nf90(nf90_open(filename_prefix // int2string(host_rank) // ".nc", &
         & NF90_NOWRITE, ncid))

    comin_version = comin_setup_get_version()
    CALL nf90(nf90_get_att(ncid, NF90_GLOBAL, "comin_version", file_comin_version))
    IF (file_comin_version(1) /= comin_version%version_no_major .OR. &
      & file_comin_version(2) /= comin_version%version_no_minor) &
      &   CALL comin_plugin_finish("comin_var_replay_plugin", "Incompatible comin versions")

    CALL nf90(nf90_get_att(ncid, NF90_GLOBAL, "host_comm_rank", file_host_rank))
    IF(file_host_rank /= host_rank) &
      CALL comin_plugin_finish("comin_var_replay_plugin", "host rank number does not match in file")
    CALL nf90(nf90_get_att(ncid, NF90_GLOBAL, "host_comm_size", file_host_size))
    IF(file_host_size /= host_comm_size) &
      CALL comin_plugin_finish("comin_var_replay_plugin", "host comm size does not match in file")

    CALL nf90(nf90_inq_dimid(ncid, "event", event_dimid))

    global => comin_descrdata_get_global()
    ALLOCATE(domain_ncids(global%n_dom))
    ALLOCATE(nvars(global%n_dom))
    DO jg=1,global%n_dom
      status = nf90_inq_ncid(ncid, "domain_"// int2string(jg), domain_ncids(jg))
      IF(status /= NF90_NOERR) CYCLE
      CALL nf90(nf90_inquire(domain_ncids(jg), nVariables=nvars(jg)))
    END DO

    ALLOCATE(var_descrs(global%n_dom,MAXVAL(nvars)))
    ALLOCATE(var_ptrs(global%n_dom,  MAXVAL(nvars)))
    ALLOCATE(var_varids(global%n_dom,MAXVAL(nvars)))

    DO jg=1,global%n_dom
      CALL nf90(nf90_inq_varids(domain_ncids(jg), nvars(jg), var_varids(jg, 1:nvars(jg))))
      DO i=1,nvars(jg)
        CALL nf90(nf90_inquire_variable(domain_ncids(jg), var_varids(jg, i), name, nAtts=nAtts))
        var_descrs(jg, i) = t_comin_var_descriptor(name=name, id=jg)
        CALL nf90(nf90_get_att(domain_ncids(jg), var_varids(jg, i), "metadata::tracer", tracer))
        IF (tracer /= 0 .AND. jg > 1) CYCLE
        tmp_descr = var_descrs(jg, i)
        IF (tracer /= 0) tmp_descr%id = -1
        CALL comin_var_request_add(tmp_descr, .FALSE.)
        DO j=1,nAtts
          CALL nf90(nf90_inq_attname(domain_ncids(jg), var_varids(jg, i), j, attname))
          IF(INDEX(attname, "metadata::") /= 1) CYCLE
          CALL nf90(nf90_inquire_attribute(domain_ncids(jg), var_varids(jg, i), attname, atttype, attlen))
          SELECT CASE (atttype)
          CASE(NF90_DOUBLE)
            CALL nf90(nf90_get_att(domain_ncids(jg), var_varids(jg, i), attname, attvalue_real))
            CALL comin_metadata_set(tmp_descr, attname(11:), attvalue_real)
          CASE(NF90_INT)
            CALL nf90(nf90_get_att(domain_ncids(jg), var_varids(jg, i), attname, attvalue_int))
            CALL comin_metadata_set(tmp_descr, attname(11:), attvalue_int)
          CASE(NF90_BYTE)
            CALL nf90(nf90_get_att(domain_ncids(jg), var_varids(jg, i), attname, attvalue_int))
            CALL comin_metadata_set(tmp_descr, attname(11:), attvalue_int /= 0)
          CASE(NF90_CHAR)
            ALLOCATE(CHARACTER(LEN=attlen) :: attvalue_char)
            CALL nf90(nf90_get_att(domain_ncids(jg), var_varids(jg, i), attname, attvalue_char))
            CALL comin_metadata_set(tmp_descr, attname(11:), attvalue_char)
            DEALLOCATE(attvalue_char)
          END SELECT
        END DO
      END DO
    END DO

    CALL comin_callback_register(EP_SECONDARY_CONSTRUCTOR, comin_var_replay_secondary_ctor)
    CALL comin_callback_register(ep, comin_var_replay)
  END SUBROUTINE comin_main

  SUBROUTINE comin_var_replay_secondary_ctor() BIND(C)
    INTEGER :: i, jg
    INTEGER :: pos_jc, pos_jk, pos_jb
    DO jg=1,global%n_dom
      DO i=1,nvars(jg)
        CALL comin_var_get([ep], var_descrs(jg, i), COMIN_FLAG_WRITE, var_ptrs(jg, i)%ptr)
        IF( .NOT. ASSOCIATED(var_ptrs(jg, i)%ptr)) THEN
          CALL comin_plugin_finish("comin_var_replay_plugin", "Cannot get_vat " // var_descrs(jg, i)%name)
        ENDIF
        ! Check if the index positions are the same
        CALL nf90(nf90_get_att(domain_ncids(jg), var_varids(jg, i), "pos_jc", pos_jc))
        IF (pos_jc /= var_ptrs(jg, i)%ptr%pos_jc) &
          CALL comin_plugin_finish("comin_var_replay_plugin", &
                                   "pos_jc does not match for var " // TRIM(var_descrs(jg, i)%name))
        CALL nf90(nf90_get_att(domain_ncids(jg), var_varids(jg, i), "pos_jk", pos_jk))
        IF (pos_jk /= var_ptrs(jg, i)%ptr%pos_jk) &
          CALL comin_plugin_finish("comin_var_replay_plugin", &
                                   "pos_jk does not match for var " // TRIM(var_descrs(jg, i)%name))
        CALL nf90(nf90_get_att(domain_ncids(jg), var_varids(jg, i), "pos_jb", pos_jb))
        IF (pos_jb /= var_ptrs(jg, i)%ptr%pos_jb) &
          CALL comin_plugin_finish("comin_var_replay_plugin", &
                                   "pos_jb does not match for var " // TRIM(var_descrs(jg, i)%name))
      END DO
    END DO
  END SUBROUTINE comin_var_replay_secondary_ctor

  SUBROUTINE comin_var_replay() BIND(C)
    INTEGER :: i, jg
    DO jg=1,global%n_dom
      DO i=1,nvars(jg)
        CALL nf90(nf90_get_var(domain_ncids(jg), var_varids(jg, i), var_ptrs(jg, i)%ptr%ptr, &
             &                 [1, 1, 1, 1, 1, event_counter]))
      END DO
    END DO
    event_counter = event_counter+1
  END SUBROUTINE comin_var_replay

END MODULE comin_var_replay_plugin
