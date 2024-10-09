!> @file comin_var_recorder_plugin.F90
!! @brief A plugin that records the values of variables to a netcdf file
!
!  @authors 01/2024 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
MODULE comin_var_recorder_plugin

  USE comin_plugin_interface, ONLY: t_comin_var_ptr, t_comin_plugin_info,                      &
       &                            t_comin_descrdata_global,                                  &
       &                            t_comin_var_descriptor, t_comin_setup_version_info,        &
       &                            comin_current_get_plugin_info,                             &
       &                            comin_parallel_get_host_mpi_rank, comin_setup_get_version, &
       &                            comin_callback_get_ep_name, comin_plugin_finish,           &
       &                            comin_callback_register, comin_descrdata_get_global,       &
       &                            comin_var_get, t_comin_var_descriptor,                     &
       &                            comin_metadata_get_typeid, comin_metadata_get,             &
       &                            EP_ATM_TIMELOOP_END, COMIN_FLAG_READ, EP_DESTRUCTOR,       &
       &                            EP_SECONDARY_CONSTRUCTOR,                                  &
       &                            comin_parallel_get_host_mpi_comm, comin_error_check,       &
       &                            COMIN_METADATA_TYPEID_INTEGER, COMIN_METADATA_TYPEID_REAL, &
       &                            COMIN_METADATA_TYPEID_CHARACTER,                           &
       &                            COMIN_METADATA_TYPEID_LOGICAL, COMIN_METADATA_TYPEID_UNDEFINED, &
       &                            t_comin_var_metadata_iterator, comin_metadata_get_iterator
  USE netcdf,                 ONLY: nf90_create, nf90_put_att, nf90_def_dim, nf90_def_grp,     &
       &                            nf90_put_var, nf90_close, NF90_GLOBAL, NF90_NETCDF4,       &
       &                            NF90_UNLIMITED, NF90_DOUBLE
  USE netcdf_utils,           ONLY: nf90, nf90_utils_def_var
  USE iso_c_binding,          ONLY: C_DOUBLE
  USE utils,                  ONLY: int2string

  IMPLICIT NONE

  TYPE :: var_ptr_ptr
    TYPE(t_comin_var_ptr), POINTER :: ptr
  END TYPE var_ptr_ptr

  ! entry point where variable values are captured.
  INTEGER, PARAMETER :: ep = EP_ATM_TIMELOOP_END

  INTEGER :: ncid
  INTEGER, ALLOCATABLE :: domain_ncids(:)
  INTEGER :: event_dimid, event_counter = 0

  TYPE(t_comin_plugin_info)                 :: plugin_info
  TYPE(t_comin_descrdata_global), POINTER   :: global

  TYPE(t_comin_var_descriptor), ALLOCATABLE :: var_descrs(:,:)
  TYPE(var_ptr_ptr),            ALLOCATABLE :: var_ptrs(:,:)
  INTEGER,                      ALLOCATABLE :: var_varids(:,:)

CONTAINS

  SUBROUTINE comin_main() BIND(C)

    INTEGER                          :: ierr
    INTEGER                          :: host_rank, host_comm, host_comm_size
    TYPE(t_comin_setup_version_info) :: comin_version
    CHARACTER(LEN=:), ALLOCATABLE    :: ep_name

    CALL comin_current_get_plugin_info(plugin_info)
    host_rank = comin_parallel_get_host_mpi_rank()
    host_comm = comin_parallel_get_host_mpi_comm()
    CALL MPI_COMM_SIZE(host_comm, host_comm_size, ierr)
    CALL nf90(nf90_create("vars_" // int2string(host_rank) // ".nc", &
                          NF90_NETCDF4, ncid))

    comin_version = comin_setup_get_version()
    CALL nf90(nf90_put_att(ncid, NF90_GLOBAL, "comin_version", &
         & [comin_version%version_no_major, &
         &  comin_version%version_no_minor, &
         &  comin_version%version_no_patch]))
    CALL nf90(nf90_put_att(ncid, NF90_GLOBAL, "ep", ep))
    CALL comin_callback_get_ep_name(ep, ep_name)
    CALL nf90(nf90_put_att(ncid, NF90_GLOBAL, "ep_name", &
         & ep_name // "(" // int2string(ep) // ")"))
    CALL nf90(nf90_put_att(ncid, NF90_GLOBAL, "host_comm_rank", host_rank))
    CALL nf90(nf90_put_att(ncid, NF90_GLOBAL, "host_comm_size", host_comm_size))

    CALL nf90(nf90_def_dim(ncid, "event", NF90_UNLIMITED, event_dimid))

    CALL comin_callback_register(EP_SECONDARY_CONSTRUCTOR, var_recorder_secondary_ctor)
    CALL comin_callback_register(ep, comin_var_recorder_callback)
    CALL comin_callback_register(EP_DESTRUCTOR, comin_var_recorder_destructor)

  END SUBROUTINE comin_main

  SUBROUTINE var_recorder_secondary_ctor() BIND(C)
    INTEGER :: jg, i
    CHARACTER(LEN=:), ALLOCATABLE :: var_names(:)
    TYPE(t_comin_var_metadata_iterator) :: metadata_it

    CALL split_string(plugin_info%options, var_names)

    global => comin_descrdata_get_global()
    ALLOCATE(var_descrs(global%n_dom,SIZE(var_names)))
    ALLOCATE(var_ptrs(global%n_dom,SIZE(var_names)))
    ALLOCATE(var_varids(global%n_dom,SIZE(var_names)))
    ALLOCATE(domain_ncids(global%n_dom))
    DO jg=1, global%n_dom
      CALL nf90(nf90_def_grp(ncid, "domain_"//int2string(jg), domain_ncids(jg)))
      DO i=1,SIZE(var_names)
        var_descrs(jg, i) = t_comin_var_descriptor(name=TRIM(var_names(i)), id=jg)
        CALL comin_var_get([ep], var_descrs(jg, i), COMIN_FLAG_READ, var_ptrs(jg, i)%ptr)

        IF(.NOT. ASSOCIATED(var_ptrs(jg, i)%ptr)) THEN
          CALL comin_plugin_finish("comin_var_recorder_plugin", &
          &                        "Pointer for variable " // var_names(i) // " not associated")
        ENDIF
        var_varids(jg, i) = nf90_utils_def_var(domain_ncids(jg), TRIM(var_names(i)), &
             & NF90_DOUBLE, SHAPE(var_ptrs(jg, i)%ptr%ptr), [event_dimid])
        CALL nf90(nf90_put_att(domain_ncids(jg), var_varids(jg, i), &
                               "pos_jc", var_ptrs(jg, i)%ptr%pos_jc))
        CALL nf90(nf90_put_att(domain_ncids(jg), var_varids(jg, i), &
                               "pos_jk", var_ptrs(jg, i)%ptr%pos_jk))
        CALL nf90(nf90_put_att(domain_ncids(jg), var_varids(jg, i), &
                               "pos_jb", var_ptrs(jg, i)%ptr%pos_jb))
        CALL nf90(nf90_put_att(domain_ncids(jg), var_varids(jg, i), &
                               "pos_jn", var_ptrs(jg, i)%ptr%pos_jn))
        CALL nf90(nf90_put_att(domain_ncids(jg), var_varids(jg, i), &
                               "ncontained", var_ptrs(jg, i)%ptr%ncontained))
        CALL nf90(nf90_put_att(domain_ncids(jg), var_varids(jg, i), &
                               "lcontained", MERGE(1,0,var_ptrs(jg, i)%ptr%lcontainer)))

        CALL comin_metadata_get_iterator(var_descrs(jg, i), metadata_it)
        DO WHILE(.NOT. metadata_it%is_end())
          CALL save_metadata(domain_ncids(jg), var_varids(jg, i), var_descrs(jg, i), metadata_it%key())
          CALL metadata_it%next()
        END DO
      END DO
    END DO
  END SUBROUTINE var_recorder_secondary_ctor

  SUBROUTINE comin_var_recorder_callback() BIND(C)
    INTEGER :: jg, i

    event_counter = event_counter + 1

    DO jg=1, global%n_dom
      DO i=1,SIZE(var_ptrs, 2)
        CALL nf90(nf90_put_var(domain_ncids(jg), var_varids(jg, i), &
             & var_ptrs(jg, i)%ptr%ptr, [1,1,1,1,1,event_counter]))
      END DO
    END DO

  END SUBROUTINE comin_var_recorder_callback

  SUBROUTINE comin_var_recorder_destructor() BIND(C)
    CALL nf90(nf90_close(ncid))
  END SUBROUTINE comin_var_recorder_destructor

  SUBROUTINE save_metadata(ncid, varid, var_descr, key)
    INTEGER, INTENT(IN) :: ncid, varid
    TYPE(t_comin_var_descriptor), INTENT(IN) :: var_descr
    CHARACTER(LEN=*), INTENT(IN) :: key

    LOGICAL :: logical_value
    INTEGER(kind=1) :: l2i_value
    INTEGER :: int_value
    CHARACTER(LEN=:), ALLOCATABLE :: char_value
    REAL(C_DOUBLE) :: real_value

    SELECT CASE (comin_metadata_get_typeid(var_descr, key))
    CASE (COMIN_METADATA_TYPEID_LOGICAL)
      CALL comin_metadata_get(var_descr, key, logical_value)
      l2i_value = MERGE(1_1, 0_1, logical_value)
      CALL nf90(nf90_put_att(ncid, varid, "metadata::" // key, l2i_value))
    CASE(COMIN_METADATA_TYPEID_INTEGER)
      CALL comin_metadata_get(var_descr, key, int_value)
      CALL nf90(nf90_put_att(ncid, varid, "metadata::" // key, int_value))
    CASE(COMIN_METADATA_TYPEID_REAL)
      CALL comin_metadata_get(var_descr, key, real_value)
      CALL nf90(nf90_put_att(ncid, varid, "metadata::" // key, real_value))
    CASE(COMIN_METADATA_TYPEID_CHARACTER)
      CALL comin_metadata_get(var_descr, key, char_value)
      CALL nf90(nf90_put_att(ncid, varid, "metadata::" // key, char_value))
    CASE(COMIN_METADATA_TYPEID_UNDEFINED)
      ! No problem if metadata was not set at all
    CASE DEFAULT
      CALL comin_plugin_finish("comin_var_recorder_plugin", "Unknown metadata: " // key)
    END SELECT
  END SUBROUTINE save_metadata

  SUBROUTINE split_string(in, out)
    CHARACTER(LEN=*),              INTENT(IN)  :: in
    CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: out(:)

    INTEGER :: l, idx = 1, no_str = 0, max_len = 0

    DO
      l = INDEX(in(idx:), ",")-1
      IF (l < 0) THEN
        l = LEN(in)-idx+1
      ENDIF
      max_len = MAX(max_len, l)
      idx = idx+l+1
      no_str = no_str +1
      IF (idx > LEN(in)) EXIT
    END DO

    ALLOCATE(CHARACTER(LEN=max_len) :: out(no_str))
    no_str = 1
    idx = 1
    DO
      l = INDEX(in(idx:), ",")-1
      IF(l < 0) THEN
        out(no_str) = in(idx:)
        EXIT
      ELSE
        out(no_str) = in(idx:idx+l-1)
      ENDIF
      idx = idx+l+1
      no_str = no_str +1
    END DO
  END SUBROUTINE split_string

END MODULE comin_var_recorder_plugin
