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

MODULE mo_ser_all

#ifdef SERIALIZE
  USE m_serialize, ONLY: fs_add_savepoint_metainfo, fs_create_savepoint
  USE utils_ppser, ONLY: ppser_savepoint
  USE mo_impl_constants,     ONLY: REAL_T, SINGLE_T, INT_T, BOOL_T, VNAME_LEN
  USE mo_exception,          ONLY: warning, message
  USE mo_ser_common,         ONLY: init, t_ser_options, ser_component, open_compare_file, close_compare_file

  USE mtime,                 ONLY: datetime, datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_time_config,        ONLY: time_config
  USE mo_dynamics_config,    ONLY: nnow, nnew, nnow_rcf, nnew_rcf
  USE mo_var_metadata_types, ONLY: t_var_metadata
  USE mo_model_domain,       ONLY: t_patch
  USE mo_var_list,           ONLY: t_var_list_ptr
  USE mo_var,                ONLY: t_var
  USE mo_var_list_register,  ONLY: vlr_get
  USE mo_run_config,         ONLY: iforcing, ldass_lhn
  USE mo_impl_constants,     ONLY: inwp
  USE mo_ser_nml,            ONLY: ser_initialization, ser_output_diag, ser_output_opt, ser_output_diag_dyn, &
                                   ser_latbc_data, ser_nesting_save_progvars, &
                                   ser_dynamics, ser_diffusion, ser_nesting_compute_tendencies, &
                                   ser_nesting_boundary_interpolation, ser_nesting_relax_feedback, &
                                   ser_step_advection, ser_turbtrans, ser_turbdiff, &
                                   ser_physics, ser_lhn, ser_nudging, ser_all_debug, ser_surface, &
                                   ser_microphysics, ser_convection, ser_cover, ser_radiation, &
                                   ser_radheat, ser_gwdrag, ser_time_loop_end, ser_physics_init, &
                                   ser_reset_to_initial_state, ser_set_reff
  USE mo_ser_manually,       ONLY: ser_manually
  USE mo_mpi,                ONLY: get_my_mpi_work_id
  USE mo_sppt_config,        ONLY: sppt_config
#endif

  IMPLICIT NONE

  PUBLIC :: serialize_all ! This is the only component that has to be available without SERIALIZE

  CONTAINS

#ifdef SERIALIZE
  PRIVATE

  SUBROUTINE char_to_hash(c, a)
      CHARACTER(len=*), INTENT(in) :: c
      INTEGER, INTENT(out) :: a
      INTEGER :: i

      INTEGER :: p = 31
      INTEGER :: m = 1e5 + 5
      INTEGER :: p_pow

      a = 0
      p_pow = 1
      DO i=1,len_trim(c)
          a = MOD(a + (ichar(c(i:i)) + 1) * p_pow, m)
          p_pow = MOD(p_pow * p, m)
      END DO

  END SUBROUTINE


  SUBROUTINE ser_var_list(name, abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain, substr, timelev)
    CHARACTER(len=*), INTENT(IN)            :: name    ! name of output var_list
    INTEGER, INTENT(IN)                     :: abs_threshold
    INTEGER, INTENT(IN)                     :: rel_threshold
    LOGICAL, INTENT(IN)                     :: lupdate_cpu
    INTEGER, INTENT(IN)                     :: ser_mode
    INTEGER, INTENT(IN), OPTIONAL           :: domain  ! domain index to append
    CHARACTER(len=*), INTENT(IN), OPTIONAL  :: substr  ! String after domain, before timelev
    INTEGER, INTENT(IN), OPTIONAL           :: timelev ! timelev index to append

    TYPE(t_var_list_ptr)          :: list
    TYPE(t_var_metadata), POINTER :: info
    TYPE(t_var),          POINTER :: element
    TYPE(t_ser_options)           :: o
    CHARACTER(len=VNAME_LEN)      :: listname

    INTEGER :: ii

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: ser_name
    INTEGER :: listhash

    listname = TRIM(name)

    o%abs_threshold = abs_threshold
    o%rel_threshold = rel_threshold
    o%lupdate_cpu = lupdate_cpu
    o%ser_mode = ser_mode

    ! Append domain index
    IF( PRESENT(domain) ) THEN
      WRITE(listname, '(a,i2.2)') TRIM(listname), domain
      o%domain = domain
    END IF

    ! Append substr
    IF( PRESENT(substr) ) THEN
      listname = TRIM(listname)//TRIM(substr)
    END IF

    ! Append timelev index
    IF( PRESENT(timelev) ) THEN
      WRITE(listname, '(a,i2.2)') TRIM(listname), timelev
    END IF


    CALL vlr_get(list, listname)


    IF ( ASSOCIATED(list%p) ) THEN
      for_all_list_elements: DO ii = 1, list%p%nvars
        element => list%p%vl(ii)%p
        info    => element%info
        call char_to_hash(listname, listhash)
        write(ser_name, '(a,i5.5)') TRIM(info%name)//'_', listhash
        o%lopenacc  = info%lopenacc

        ! Contained variables are (/should be) already serialized through their container.
        IF ( info%lcontained ) THEN 
          ! just check meta data
          o%lopenacc = .FALSE.
          o%lupdate_cpu = .FALSE.
          call ser_component(o, TRIM(ser_name)//"_contained", info%used_dimensions)
          CYCLE for_all_list_elements
        END IF

        ! only compare shape of empty variables
        IF ( ANY(info%used_dimensions==0) ) THEN
          o%lopenacc = .FALSE.
          o%lupdate_cpu = .FALSE.
          call ser_component(o, TRIM(ser_name)//"_shape", info%used_dimensions)
          CYCLE for_all_list_elements
        END IF

        SELECT CASE(info%data_type)
          CASE (REAL_T)
            SELECT CASE ( info%ndims )
              CASE(1)
                call ser_component(o, TRIM(ser_name), element%r_ptr(:,1,1,1,1))
              CASE(2)
                call ser_component(o, TRIM(ser_name), element%r_ptr(:,:,1,1,1))
              CASE(3)
                call ser_component(o, TRIM(ser_name), element%r_ptr(:,:,:,1,1))
              CASE(4)
                call ser_component(o, TRIM(ser_name), element%r_ptr(:,:,:,:,1))
              END SELECT
          CASE (SINGLE_T)
            SELECT CASE ( info%ndims )
              CASE(1)
                call ser_component(o, TRIM(ser_name), element%s_ptr(:,1,1,1,1))
              CASE(2)
                call ser_component(o, TRIM(ser_name), element%s_ptr(:,:,1,1,1))
              CASE(3)
                call ser_component(o, TRIM(ser_name), element%s_ptr(:,:,:,1,1))
              CASE(4)
                call ser_component(o, TRIM(ser_name), element%s_ptr(:,:,:,:,1))
            END SELECT
          CASE (INT_T)
            SELECT CASE ( info%ndims )
              CASE(1)
                call ser_component(o, TRIM(ser_name), element%i_ptr(:,1,1,1,1))
              CASE(2)
                call ser_component(o, TRIM(ser_name), element%i_ptr(:,:,1,1,1))
              CASE(3)
                call ser_component(o, TRIM(ser_name), element%i_ptr(:,:,:,1,1))
              CASE(4)
                call ser_component(o, TRIM(ser_name), element%i_ptr(:,:,:,:,1))
            END SELECT
          CASE (BOOL_T)
            SELECT CASE ( info%ndims )
              CASE(1)
                call ser_component(o, TRIM(ser_name), element%l_ptr(:,1,1,1,1))
              CASE(2)
                call ser_component(o, TRIM(ser_name), element%l_ptr(:,:,1,1,1))
              CASE(3)
                call ser_component(o, TRIM(ser_name), element%l_ptr(:,:,:,1,1))
              CASE(4)
                call ser_component(o, TRIM(ser_name), element%l_ptr(:,:,:,:,1))
            END SELECT
        END SELECT

      ENDDO for_all_list_elements
    END IF

  END SUBROUTINE ser_var_list
#endif

  SUBROUTINE serialize_all(nproma, jg, savepoint_base, is_input, opt_lupdate_cpu, opt_id, opt_dt)
    USE mtime, ONLY: datetime

    INTEGER, INTENT(IN) :: nproma, jg
    CHARACTER(LEN=*), INTENT(IN) :: savepoint_base
    LOGICAL, INTENT(IN) :: is_input

    LOGICAL, INTENT(IN), OPTIONAL :: opt_lupdate_cpu
    INTEGER, INTENT(IN), OPTIONAL :: opt_id
    ! use this to pass a datetime that describes the exact current sub-timestep of a nested domain.
    TYPE(datetime), INTENT(IN), OPTIONAL, POINTER :: opt_dt 

#ifdef SERIALIZE
    LOGICAL :: lupdate_cpu
    INTEGER :: id, ser_mode
    INTEGER, POINTER :: ser_setting(:)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: date
    CHARACTER(LEN=LEN(savepoint_base)+4) :: savepoint_name
    CHARACTER(LEN=100) :: compare_file_name
    LOGICAL :: do_serialization
    INTEGER :: rel_threshold, abs_threshold

    do_serialization = .TRUE.

    ! make sure CPU and GPU are synchronous
    !$ACC WAIT

    SELECT CASE(savepoint_base)
      CASE("initialization")
        ser_setting => ser_initialization
        ! this is serialized just once
        IF(is_input) &
          CALL warning('mo_ser:all:serialize_all:initialization', 'initialization should not be used as input')
      CASE("output_diag_dyn")
        ser_setting => ser_output_diag_dyn
      CASE("output_diag")
        ser_setting => ser_output_diag
      CASE("output_opt")
        ser_setting => ser_output_opt
      CASE("latbc_data")
        ser_setting => ser_latbc_data
      CASE("nesting_save_progvars")
        ser_setting => ser_nesting_save_progvars
      CASE("dynamics")
        ser_setting => ser_dynamics
      CASE("diffusion")
        ser_setting => ser_diffusion
      CASE("step_advection")
        ser_setting => ser_step_advection
      CASE("nesting_compute_tendencies")
        ser_setting => ser_nesting_compute_tendencies
      CASE("nesting_boundary_interpolation")
        ser_setting => ser_nesting_boundary_interpolation
      CASE("nesting_relax_feedback")
        ser_setting => ser_nesting_relax_feedback
      CASE("physics")
        ser_setting => ser_physics
      CASE("physics_init")
        ser_setting => ser_physics_init
      CASE("lhn")
        ser_setting => ser_lhn
      CASE("nudging")
        ser_setting => ser_nudging
      CASE("surface")
        ser_setting => ser_surface
      CASE("microphysics")
        ser_setting => ser_microphysics
      CASE("turbdiff")
        ser_setting => ser_turbdiff
      CASE("turbtrans")
        ser_setting => ser_turbtrans
      CASE("convection")
        ser_setting => ser_convection
      CASE("cover")
        ser_setting => ser_cover
      CASE("set_reff")
        ser_setting => ser_set_reff
      CASE("radiation")
        ser_setting => ser_radiation
      CASE("radheat")
        ser_setting => ser_radheat
      CASE("gwdrag")
        ser_setting => ser_gwdrag
      CASE("time_loop_end")
        ser_setting => ser_time_loop_end
      CASE("reset_to_initial_state")
        ser_setting => ser_reset_to_initial_state
      CASE DEFAULT
        CALL warning('SER','Use default ser_all_debug settings for savepoint_base = '//savepoint_base)
        ser_setting => ser_all_debug
   END SELECT

    rel_threshold = ser_setting(2)
    abs_threshold = ser_setting(3)

    IF(ser_setting(1)/=0) THEN

       IF(is_input) THEN
#ifdef SERIALIZE_CREATE_REFERENCE
          ser_mode = 0 ! write
#elif SERIALIZE_READ_REFERENCE
          ser_mode = 1 ! read
#elif SERIALIZE_PERTURB_REFERENCE
          ser_mode = 2 ! read perturb
#endif
          savepoint_name = savepoint_base//"-in"
       ELSE
          savepoint_name = savepoint_base//"-out"
          ser_setting(1) = ser_setting(1) - 1 ! Reduce # of times to serialize
#ifdef SERIALIZE_CREATE_REFERENCE
          ser_mode = 0 ! write
#elif SERIALIZE_READ_REFERENCE
          ser_mode = 3 ! compare
#elif SERIALIZE_PERTURB_REFERENCE
          ser_mode = 3 ! compare
#endif
       ENDIF

       ! only update CPU data if requested and for write mode
       IF(PRESENT(opt_lupdate_cpu) .AND. ser_mode==0) THEN
           lupdate_cpu = opt_lupdate_cpu
       ELSE
           lupdate_cpu = .FALSE.
       ENDIF

       IF(PRESENT(opt_id)) THEN
           id = opt_id
       ELSE
           id = 0
       ENDIF

       CALL warning('SER:'//TRIM(savepoint_name),'Serialization is active!')

       IF(PRESENT(opt_dt)) THEN
         CALL datetimeToString(opt_dt, date)
       ELSE
         ! time_config%tc_current_date stores the current time of DOM01 only
         CALL datetimeToString(time_config%tc_current_date, date)
       END IF

       WRITE(compare_file_name, "(A,A1,A,A1,I2.2,A4,I2.2,A5,I2.2)") savepoint_name, "_", TRIM(date), &
         "_", id, "_dom", jg, "_rank", get_my_mpi_work_id()
       CALL message('SER working on',compare_file_name)

#ifdef _OPENACC
       IF(lupdate_cpu) THEN
         CALL warning('GPU:'//TRIM(savepoint_name),'GPU HOST synchronization forced by serialization!')
       ENDIF
#endif

       IF(ser_mode == 3) THEN
         CALL open_compare_file(TRIM(compare_file_name))
       ELSE

       ENDIF

       CALL init('icon')
       call fs_create_savepoint(TRIM(savepoint_name), ppser_savepoint)
       call fs_add_savepoint_metainfo(ppser_savepoint, 'nproma', nproma)
       call fs_add_savepoint_metainfo(ppser_savepoint, 'date', TRIM(date))
       call fs_add_savepoint_metainfo(ppser_savepoint, 'id', id)
       call fs_add_savepoint_metainfo(ppser_savepoint, 'domain', jg)

       IF(iforcing == inwp) THEN
           ! Serialize NWP fields
           CALL ser_var_list('prm_diag_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg)
           CALL ser_var_list('prm_tend_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg)
           CALL ser_var_list('lnd_prog_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg, substr='_and_timelev_', timelev=nnow(jg))
           CALL ser_var_list('lnd_prog_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg, substr='_and_timelev_', timelev=nnew(jg))
           CALL ser_var_list('lnd_diag_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg)
           CALL ser_var_list('wtr_prog_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg, substr='_and_timelev_', timelev=nnow(jg))
           CALL ser_var_list('wtr_prog_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg, substr='_and_timelev_', timelev=nnew(jg))
           CALL ser_var_list('ext_data_atm_td_D', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg)
           CALL ser_var_list('ext_data_atm_D', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg)

           ! Serialize sppt
           IF (sppt_config(jg)%lsppt) THEN
            CALL ser_var_list('sppt_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg)
           END IF
       ENDIF

       IF(ldass_lhn) THEN
           ! Serialize radar data and LHN fields
           CALL ser_var_list('radar_data_ct_dom_',    abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg)
           CALL ser_var_list('radar_data_td_dom_',    abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg)
           CALL ser_var_list('lhn_fields_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg)
       ENDIF

       ! Serialize dynamics fields
       CALL ser_var_list('nh_state_metrics_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg)
       CALL ser_var_list('nh_state_diag_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg)
       CALL ser_var_list('nh_state_prog_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg, substr='_and_timelev_', timelev=nnew(jg)) !p_prog
       CALL ser_var_list('nh_state_prog_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg, substr='_and_timelev_', timelev=nnow(jg)) !p_prog
       IF(nnow_rcf(jg) /= nnew(jg) .AND. nnow_rcf(jg) /= nnow(jg)) THEN
         CALL ser_var_list('nh_state_prog_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg, substr='_and_timelev_', timelev=nnow_rcf(jg)) !p_prog_now_rcf
       ENDIF
       IF(nnew_rcf(jg) /= nnew(jg) .AND. nnew_rcf(jg) /= nnow_rcf(jg) .AND. nnew_rcf(jg) /= nnow(jg)) THEN
         CALL ser_var_list('nh_state_prog_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg, substr='_and_timelev_', timelev=nnew_rcf(jg)) !p_prog_rcf
       ENDIF

       ! Serialize fields that prepare advection
       CALL ser_var_list('prepadv_of_domain_', abs_threshold, rel_threshold, lupdate_cpu, ser_mode, domain=jg)

       ! Serialize selected scalars and arrays that are not part of any list
       CALL ser_manually(abs_threshold, rel_threshold, lupdate_cpu, ser_mode, jg)

       IF(ser_mode == 3) THEN
         CALL close_compare_file()
       ENDIF


    ENDIF !do_serialization
#endif

  END SUBROUTINE serialize_all

END MODULE mo_ser_all
