!
! @brief Main program for the ICON atmospheric model
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

MODULE mo_test_communication

  USE mo_kind,                ONLY: wp, sp, dp, vp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_mpi,                 ONLY: work_mpi_barrier, p_n_work, p_pe_work, &
    &                               p_comm_work, p_comm_rank, p_comm_size, &
    &                               p_barrier
  USE mpi
  USE mo_timer,               ONLY: ltimer, new_timer, timer_start, timer_stop, &
    & print_timer, activate_sync_timers, timers_level, timer_barrier, timer_radiaton_recv
  USE mo_parallel_config,     ONLY: nproma, icon_comm_method
  USE mo_communication,       ONLY: t_comm_gather_pattern, t_comm_pattern, &
    &                               delete_comm_pattern, &
    &                               setup_comm_gather_pattern, &
    &                               delete_comm_gather_pattern, exchange_data, &
    &                               exchange_data_4de1, &
    &                               exchange_data_mult, &
    &                               exchange_data_mult_mixprec, &
    &                               t_comm_allgather_pattern, &
    &                               setup_comm_allgather_pattern, &
    &                               delete_comm_allgather_pattern, &
    &                               exchange_data_grf, t_p_comm_pattern, &
    &                               t_comm_pattern_collection, &
    &                               delete_comm_pattern_collection
  USE mo_communication_factory, ONLY: setup_comm_pattern, &
    &                               setup_comm_pattern_collection
  USE mo_decomposition_tools, ONLY: t_glb2loc_index_lookup, &
    &                               init_glb2loc_index_lookup, &
    &                               set_inner_glb_index, &
    &                               deallocate_glb2loc_index_lookup

  USE mo_master_control,      ONLY: get_my_process_name

  USE mo_model_domain,        ONLY: p_patch
  USE mo_atmo_model,          ONLY: construct_atmo_model, destruct_atmo_model
  USE mo_math_gradients,      ONLY: grad_fd_norm

  !nh utils
  USE mo_nonhydro_state,      ONLY: p_nh_state
  USE mo_atmo_nonhydrostatic, ONLY: construct_atmo_nonhydrostatic, destruct_atmo_nonhydrostatic
  USE mo_async_latbc_types,   ONLY: t_latbc_data

  USE mo_parallel_config,    ONLY: iorder_sendrecv
  USE mo_sync,               ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array, &
    &                              sync_patch_array_mult, sync_patch_array_4de1
  USE mo_icon_comm_lib

  USE mo_icon_testbed_config, ONLY: testbed_model, test_halo_communication, &
    & testbed_iterations, calculate_iterations, test_gather_communication, &
    & test_exchange_communication, test_bench_exchange_data_mult
  USE mo_grid_config, ONLY: n_dom_start




!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

PUBLIC :: test_communication
PUBLIC :: halo_communication_3D_testbed, gather_communication_testbed
PUBLIC :: exchange_communication_testbed, exchange_communication_grf_testbed
PUBLIC :: bench_exchange_data_mult

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE test_communication(namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    TYPE(t_latbc_data) :: latbc !< data structure for async latbc prefetching

    CHARACTER(*), PARAMETER :: method_name = "mo_test_communication:test_communication"


    !---------------------------------------------------------------------

    write(0,*) TRIM(get_my_process_name()), ': Start of ', method_name

    ltimer = .false.
    timers_level = 0
    activate_sync_timers = .false.
    CALL construct_atmo_model(namelist_filename,shr_namelist_filename)
    CALL construct_atmo_nonhydrostatic(latbc)
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! DO the tests
!     ltimer = .true.
!     activate_sync_timers = .true.
    ltimer = .false.
    timers_level = 0
    activate_sync_timers = .false.

    SELECT CASE(testbed_model)

    CASE(test_halo_communication)
      CALL halo_communication_3D_testbed()

    CASE(test_gather_communication)
      CALL gather_communication_testbed()

    CASE(test_exchange_communication)
      CALL exchange_communication_testbed()
      CALL exchange_communication_grf_testbed()





    CASE(test_bench_exchange_data_mult)
      CALL bench_exchange_data_mult()

    CASE default
      CALL finish(method_name, "Unrecognized communication testbed_model")

    END SELECT
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! print the timers
!    IF (my_process_is_stdio()) THEN
    CALL message("===================", "=======================")
    WRITE(message_text,*) "Communication Iterations=", testbed_iterations
    CALL message(method_name, TRIM(message_text))
    CALL print_timer()
!    ENDIF
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! Carry out the shared clean-up processes
    !---------------------------------------------------------------------
    CALL destruct_atmo_nonhydrostatic(latbc)
    CALL destruct_atmo_model()
!     CALL destruct_icon_communication()
    CALL message(TRIM(method_name),'clean-up finished')



  END SUBROUTINE test_communication
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE halo_communication_3D_testbed()

    !---------------------------------------------------------------------
    ! test the 3D sync
    !---------------------------------------------------------------------
    iorder_sendrecv = 1
    CALL test_oldsync_cells_3D("sync_1_1")
    CALL test_oldsync_edges_3D("sync_1_1")

    iorder_sendrecv = 3
    CALL test_oldsync_cells_3D("sync_1_3")
    CALL test_oldsync_edges_3D("sync_1_3")

    icon_comm_method = 1
    CALL test_iconcom_cells_3D("comm_1")
    CALL test_iconcom_edges_3D("comm_1")

    icon_comm_method = 2
    CALL test_iconcom_cells_3D("com_2")
    CALL test_iconcom_edges_3D("com_2")

    icon_comm_method = 3
    CALL test_iconcom_cells_3D("com_3")
    CALL test_iconcom_edges_3D("com_3")

    icon_comm_method = 102
    CALL test_iconcom_cells_3D("com_102")
    CALL test_iconcom_edges_3D("com_102")

    icon_comm_method = 103
    CALL test_iconcom_cells_3D("com_103")
    CALL test_iconcom_edges_3D("com_103")

    RETURN
  END SUBROUTINE halo_communication_3D_testbed
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE test_iconcom_cells_3D(timer_descr)
    CHARACTER(len=*), INTENT(in) :: timer_descr

    ! 3D variables
    REAL(wp), POINTER :: pnt_3D_cells_1(:,:,:), pnt_3D_cells_2(:,:,:), pnt_3D_cells_3(:,:,:), &
      & pnt_3D_cells_4(:,:,:)

    INTEGER :: timer_3D_cells_1, timer_3D_cells_2, timer_3D_cells_3, timer_3D_cells_4

    pnt_3D_cells_1 => p_nh_state(n_dom_start)%prog(1)%w(:,:,:)
    pnt_3D_cells_2 => p_nh_state(n_dom_start)%prog(1)%rho(:,:,:)
    pnt_3D_cells_3 => p_nh_state(n_dom_start)%prog(1)%exner(:,:,:)
    pnt_3D_cells_4 => p_nh_state(n_dom_start)%prog(1)%theta_v(:,:,:)

    pnt_3D_cells_1(:,:,:) = 0.0_wp
    pnt_3D_cells_2(:,:,:) = 0.0_wp
    pnt_3D_cells_3(:,:,:) = 0.0_wp
    pnt_3D_cells_4(:,:,:) = 0.0_wp

    timer_3D_cells_1  = new_timer  (timer_descr//"_3dcells_1")
    CALL test_iconcom_3D(p_patch(n_dom_start)%sync_cells_not_in_domain, &
      & var1=pnt_3D_cells_1, &
      & timer_id=timer_3D_cells_1)

    timer_3D_cells_2  = new_timer  (timer_descr//"_3dcells_2")
    CALL test_iconcom_3D(p_patch(n_dom_start)%sync_cells_not_in_domain, &
      & var1=pnt_3D_cells_1,var2=pnt_3D_cells_2, &
      & timer_id=timer_3D_cells_2)

    timer_3D_cells_3  = new_timer  (timer_descr//"_3dcells_3")
    CALL test_iconcom_3D(p_patch(n_dom_start)%sync_cells_not_in_domain, &
      & var1=pnt_3D_cells_1,var2=pnt_3D_cells_2, &
      & var3=pnt_3D_cells_3, &
      & timer_id=timer_3D_cells_3)

    timer_3D_cells_4  = new_timer  (timer_descr//"_3dcells_4")
    CALL test_iconcom_3D(p_patch(n_dom_start)%sync_cells_not_in_domain, &
      & var1=pnt_3D_cells_1,var2=pnt_3D_cells_2, &
      & var3=pnt_3D_cells_3,var4=pnt_3D_cells_4, &
      & timer_id=timer_3D_cells_4)

  END SUBROUTINE test_iconcom_cells_3D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE test_iconcom_edges_3D(timer_descr)
    CHARACTER(len=*), INTENT(in) :: timer_descr

    ! 3D variables
    REAL(vp), POINTER :: pnt_3D_edges_1(:,:,:), pnt_3D_edges_3(:,:,:)
    REAL(wp), POINTER :: pnt_3D_edges_2(:,:,:)

    INTEGER :: timer_3D_edges_1, timer_3D_edges_2, timer_3D_edges_3

    pnt_3D_edges_1 => p_nh_state(n_dom_start)%diag%ddt_vn_phy(:,:,:)
    pnt_3D_edges_2 => p_nh_state(n_dom_start)%diag%mass_fl_e(:,:,:)
    pnt_3D_edges_3 => p_nh_state(n_dom_start)%diag%vt(:,:,:)
!     pnt_3D_edges_4 => p_nh_state(n_dom_start)%diag%hfl_tracer(:,:,:,1)
    pnt_3D_edges_1(:,:,:) = 0.0_vp
    pnt_3D_edges_2(:,:,:) = 0.0_wp
    pnt_3D_edges_3(:,:,:) = 0.0_vp
!     pnt_3D_edges_4(:,:,:) = 0.0_wp

    timer_3D_edges_1  = new_timer  (timer_descr//"_3dedges_1")
    CALL test_iconcom_3D(p_patch(n_dom_start)%sync_cells_not_in_domain, &
      & var1=pnt_3D_edges_2, timer_id=timer_3D_edges_1)

!    timer_3D_edges_1  = new_timer  (timer_descr//"_3dedges_1")
!    CALL test_iconcom_3D(p_patch(n_dom_start)%sync_cells_not_in_domain, &
!      & var1=pnt_3D_edges_1, timer_id=timer_3D_edges_1)

!    timer_3D_edges_2  = new_timer  (timer_descr//"_3dedges_2")
!    CALL test_iconcom_3D(p_patch(n_dom_start)%sync_cells_not_in_domain, &
!      & var1=pnt_3D_edges_1,var2=pnt_3D_edges_2, &
!      & timer_id=timer_3D_edges_2)

!    timer_3D_edges_3  = new_timer  (timer_descr//"_3dedges_3")
!    CALL test_iconcom_3D(p_patch(n_dom_start)%sync_cells_not_in_domain, &
!      & var1=pnt_3D_edges_1,var2=pnt_3D_edges_2, &
!      & var3=pnt_3D_edges_3, &
!      & timer_id=timer_3D_edges_3)

!     timer_3D_edges_4  = new_timer  (timer_descr//"_3dedges_4")
!     CALL test_iconcom_3D(p_patch(n_dom_start)%sync_cells_not_in_domain, &
!       & var1=pnt_3D_edges_1,var2=pnt_3D_edges_2, &
!       & var3=pnt_3D_edges_3,var4=pnt_3D_edges_4, &
!       & timer_id=timer_3D_edges_4)

  END SUBROUTINE test_iconcom_edges_3D
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE test_oldsync_cells_3D(timer_descr)
    CHARACTER(len=*), INTENT(in) :: timer_descr

    ! 3D variables
    REAL(wp), POINTER :: pnt_3D_cells_1(:,:,:), pnt_3D_cells_2(:,:,:), pnt_3D_cells_3(:,:,:), &
      & pnt_3D_cells_4(:,:,:)
    REAL(wp), ALLOCATABLE :: pnt_4D_cells(:,:,:,:)

    INTEGER :: timer_3D_cells_1, timer_3D_cells_2, timer_3D_cells_3, &
      &        timer_3D_cells_4, timer_4DE1_cells

    INTEGER :: timer_sync_1_2_3D_cells

    INTEGER :: i, j, k, l, m, n
    CHARACTER(len=128) :: str_i

    pnt_3D_cells_1 => p_nh_state(n_dom_start)%prog(1)%w(:,:,:)
    pnt_3D_cells_2 => p_nh_state(n_dom_start)%prog(1)%rho(:,:,:)
    pnt_3D_cells_3 => p_nh_state(n_dom_start)%prog(1)%exner(:,:,:)
    pnt_3D_cells_4 => p_nh_state(n_dom_start)%prog(1)%theta_v(:,:,:)
    pnt_3D_cells_1(:,:,:) = 0.0_wp
    pnt_3D_cells_2(:,:,:) = 0.0_wp
    pnt_3D_cells_3(:,:,:) = 0.0_wp
    pnt_3D_cells_4(:,:,:) = 0.0_wp

    timer_3D_cells_1  = new_timer  (timer_descr//"_3dcells_1")
    CALL test_oldsync_3D(SYNC_C, var1=pnt_3D_cells_1, &
      & timer_id=timer_3D_cells_1)

    timer_3D_cells_2  = new_timer  (timer_descr//"_3dcells_2")
    CALL test_oldsync_3D(SYNC_C, var1=pnt_3D_cells_1,var2=pnt_3D_cells_2, &
      & timer_id=timer_3D_cells_2)

    timer_3D_cells_3  = new_timer  (timer_descr//"_3dcells_3")
    CALL test_oldsync_3D(SYNC_C, var1=pnt_3D_cells_1,var2=pnt_3D_cells_2, &
      & var3=pnt_3D_cells_3, &
      & timer_id=timer_3D_cells_3)

    timer_3D_cells_4  = new_timer  (timer_descr//"_3dcells_4")
    CALL test_oldsync_3D(SYNC_C, var1=pnt_3D_cells_1,var2=pnt_3D_cells_2, &
      & var3=pnt_3D_cells_3,var4=pnt_3D_cells_4, &
      & timer_id=timer_3D_cells_4)

    DO i = 1, 16
      ALLOCATE(pnt_4D_cells(i, SIZE(pnt_3D_cells_1, 1), &
        &                      SIZE(pnt_3D_cells_1, 2), &
        &                      SIZE(pnt_3D_cells_1, 3)))

      n = 1
      DO j = 1, SIZE(pnt_4D_cells, 1)
        DO k = 1, SIZE(pnt_4D_cells, 2)
          DO l = 1, SIZE(pnt_4D_cells, 3)
            DO m = 1, SIZE(pnt_4D_cells, 4)
              pnt_4D_cells(j, k, l, m) = n
              n = n + 1
            END DO
          END DO
        END DO
      END DO
      write (str_i, '(I4)') i
      timer_4DE1_cells = new_timer(timer_descr//"_4de1cells_"//ADJUSTL(TRIM(str_i)))
      CALL test_oldsync_4DE1(SYNC_C, var=pnt_4D_cells, nfields=i, &
        &                    timer_id=timer_4DE1_cells)
      DEALLOCATE(pnt_4D_cells)
    END DO

  END SUBROUTINE test_oldsync_cells_3D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE test_oldsync_edges_3D(timer_descr)
    CHARACTER(len=*), INTENT(in) :: timer_descr

    ! 3D variables
    REAL(vp), POINTER :: pnt_3D_edges_1(:,:,:), pnt_3D_edges_3(:,:,:)
    REAL(wp), POINTER :: pnt_3D_edges_2(:,:,:)

    INTEGER :: timer_3D_edges_1, timer_3D_edges_2, timer_3D_edges_3

    pnt_3D_edges_1 => p_nh_state(n_dom_start)%diag%ddt_vn_phy(:,:,:)
    pnt_3D_edges_2 => p_nh_state(n_dom_start)%diag%mass_fl_e(:,:,:)
    pnt_3D_edges_3 => p_nh_state(n_dom_start)%diag%vt(:,:,:)
!     pnt_3D_edges_4 => p_nh_state(n_dom_start)%diag%hfl_tracer(:,:,:,1)
    pnt_3D_edges_1(:,:,:) = 0.0_wp
    pnt_3D_edges_2(:,:,:) = 0.0_wp
    pnt_3D_edges_3(:,:,:) = 0.0_wp
!     pnt_3D_edges_4(:,:,:) = 0.0_wp

    timer_3D_edges_1  = new_timer  (timer_descr//"_3dedges_1")
    CALL test_oldsync_3D(SYNC_E, var1=pnt_3D_edges_2, &
      & timer_id=timer_3D_edges_1)

!    timer_3D_edges_1  = new_timer  (timer_descr//"_3dedges_1")
!    CALL test_oldsync_3D(SYNC_E, var1=pnt_3D_edges_1, &
!      & timer_id=timer_3D_edges_1)

!    timer_3D_edges_2  = new_timer  (timer_descr//"_3dedges_2")
!    CALL test_oldsync_3D(SYNC_E, var1=pnt_3D_edges_1,var2=pnt_3D_edges_2, &
!      & timer_id=timer_3D_edges_2)

!    timer_3D_edges_3  = new_timer  (timer_descr//"_3dedges_3")
!    CALL test_oldsync_3D(SYNC_E, var1=pnt_3D_edges_1,var2=pnt_3D_edges_2, &
!      & var3=pnt_3D_edges_3, &
!      & timer_id=timer_3D_edges_3)

!     timer_3D_edges_4  = new_timer  (timer_descr//"_3dedges_4")
!     CALL test_oldsync_3D(SYNC_E, var1=pnt_3D_edges_1,var2=pnt_3D_edges_2, &
!       & var3=pnt_3D_edges_3,var4=pnt_3D_edges_4, &
!       & timer_id=timer_3D_edges_4)

  END SUBROUTINE test_oldsync_edges_3D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE test_communication_2D

    ! 2D variables
    REAL(wp), POINTER :: pnt_2D_cells(:,:), pnt_2D_edges(:,:), pnt_2D_verts(:,:)

    INTEGER :: timer_sync_1_1_2D_cells, timer_sync_1_1_2D_edges, timer_sync_1_1_2D_verts, &
      & timer_sync_1_1_2D_all
    INTEGER :: timer_sync_1_2_2D_cells, timer_sync_1_2_2D_edges, timer_sync_1_2_2D_verts, &
      & timer_sync_1_2_2D_all
!     INTEGER :: timer_sync_2_1_2D_cells, timer_sync_2_1_2D_edges, timer_sync_2_1_2D_verts, &
!       & timer_sync_2_1_2D_all
!     INTEGER :: timer_sync_2_2_2D_cells, timer_sync_2_2_2D_edges, timer_sync_2_2_2D_verts, &
!       & timer_sync_2_2_2D_all
!     INTEGER :: timer_sync_3_1_2D_cells, timer_sync_3_1_2D_edges, timer_sync_3_1_2D_verts, &
!       & timer_sync_3_1_2D_all
!     INTEGER :: timer_sync_3_2_2D_cells, timer_sync_3_2_2D_edges, timer_sync_3_2_2D_verts, &
!       & timer_sync_3_2_2D_all

    INTEGER :: timer_iconcom_2D_cells, timer_iconcom_2D_edges, timer_iconcom_2D_verts, &
      & timer_iconcom_2D_all, timer_iconcom_2D_comb, timer_iconcom_2D_keep


    ! other
    INTEGER :: comm_1, comm_2, comm_3

    INTEGER :: i


    CHARACTER(*), PARAMETER :: method_name = "mo_test_communication:test_communication"

    !---------------------------------------------------------------------

    pnt_2D_cells => p_patch(n_dom_start)%cells%area(:,:)
    pnt_2D_edges => p_patch(n_dom_start)%edges%primal_edge_length(:,:)
    pnt_2D_verts => p_patch(n_dom_start)%verts%dual_area(:,:)
    !---------------------------------------------------------------------
    ! test the 2D sync
    !---------------------------------------------------------------------
    iorder_sendrecv = 1
    timer_sync_1_1_2D_cells  = new_timer  ("sync_1_1_2D_cells")
    timer_sync_1_1_2D_edges  = new_timer  ("sync_1_1_2D_edges")
    timer_sync_1_1_2D_verts  = new_timer  ("sync_1_1_2D_verts")
    timer_sync_1_1_2D_all    = new_timer  ("sync_1_1_2D_all")

    CALL work_mpi_barrier()
    CALL test_sync_2D( SYNC_C, pnt_2D_cells, timer_sync_1_1_2D_cells)
    CALL work_mpi_barrier()
    CALL test_sync_2D( SYNC_E, pnt_2D_edges, timer_sync_1_1_2D_edges)
    CALL work_mpi_barrier()
    CALL test_sync_2D( SYNC_V, pnt_2D_verts, timer_sync_1_1_2D_verts)
    CALL work_mpi_barrier()
    CALL test_sync_2D_all(pnt_2D_cells, pnt_2D_edges, pnt_2D_verts, timer_sync_1_1_2D_all)

    !---------------------------------------------------------------------
    iorder_sendrecv = 2
    timer_sync_1_2_2D_cells  = new_timer  ("sync_1_2_2D_cells")
    timer_sync_1_2_2D_edges  = new_timer  ("sync_1_2_2D_edges")
    timer_sync_1_2_2D_verts  = new_timer  ("sync_1_2_2D_verts")
    timer_sync_1_2_2D_all    = new_timer  ("sync_1_2_2D_all")

    CALL work_mpi_barrier()
    CALL test_sync_2D( SYNC_C, pnt_2D_cells, timer_sync_1_2_2D_cells)
    CALL work_mpi_barrier()
    CALL test_sync_2D( SYNC_E, pnt_2D_edges, timer_sync_1_2_2D_edges)
    CALL work_mpi_barrier()
    CALL test_sync_2D( SYNC_V, pnt_2D_verts, timer_sync_1_2_2D_verts)
    CALL work_mpi_barrier()
    CALL test_sync_2D_all(pnt_2D_cells, pnt_2D_edges, pnt_2D_verts, timer_sync_1_2_2D_all)

    !---------------------------------------------------------------------
!     itype_comm = 2
!     iorder_sendrecv = 1
!     timer_sync_2_1_2D_cells  = new_timer  ("sync_2_1_2D_cells")
!     timer_sync_2_1_2D_edges  = new_timer  ("sync_2_1_2D_edges")
!     timer_sync_2_1_2D_verts  = new_timer  ("sync_2_1_2D_verts")
!     timer_sync_2_1_2D_all    = new_timer  ("sync_2_1_2D_all")
!
!     CALL test_sync_2D( SYNC_C, pnt_2D_cells, timer_sync_2_1_2D_cells)
!     CALL test_sync_2D( SYNC_E, pnt_2D_edges, timer_sync_2_1_2D_edges)
!     CALL test_sync_2D( SYNC_V, pnt_2D_verts, timer_sync_2_1_2D_verts)
!     CALL test_sync_2D_all(pnt_2D_cells, pnt_2D_edges, pnt_2D_verts, timer_sync_2_1_2D_all)
!
!     !---------------------------------------------------------------------
!     itype_comm = 2
!     iorder_sendrecv = 2
!     timer_sync_2_2_2D_cells  = new_timer  ("sync_2_2_2D_cells")
!     timer_sync_2_2_2D_edges  = new_timer  ("sync_2_2_2D_edges")
!     timer_sync_2_2_2D_verts  = new_timer  ("sync_2_2_2D_verts")
!     timer_sync_2_2_2D_all    = new_timer  ("sync_2_2_2D_all")
!
!     CALL test_sync_2D( SYNC_C, pnt_2D_cells, timer_sync_2_2_2D_cells)
!     CALL test_sync_2D( SYNC_E, pnt_2D_edges, timer_sync_2_2_2D_edges)
!     CALL test_sync_2D( SYNC_V, pnt_2D_verts, timer_sync_2_2_2D_verts)
!     CALL test_sync_2D_all(pnt_2D_cells, pnt_2D_edges, pnt_2D_verts, timer_sync_2_2_2D_all)
!
!     !---------------------------------------------------------------------
!     itype_comm = 3
!     iorder_sendrecv = 1
!     timer_sync_3_1_2D_cells  = new_timer  ("sync_3_1_2D_cells")
!     timer_sync_3_1_2D_edges  = new_timer  ("sync_3_1_2D_edges")
!     timer_sync_3_1_2D_verts  = new_timer  ("sync_3_1_2D_verts")
!     timer_sync_3_1_2D_all    = new_timer  ("sync_3_1_2D_all")
!
!     CALL test_sync_2D( SYNC_C, pnt_2D_cells, timer_sync_3_1_2D_cells)
!     CALL test_sync_2D( SYNC_E, pnt_2D_edges, timer_sync_3_1_2D_edges)
!     CALL test_sync_2D( SYNC_V, pnt_2D_verts, timer_sync_3_1_2D_verts)
!     CALL test_sync_2D_all(pnt_2D_cells, pnt_2D_edges, pnt_2D_verts, timer_sync_3_1_2D_all)
!
!     !---------------------------------------------------------------------
!     itype_comm = 3
!     iorder_sendrecv = 2
!     timer_sync_3_2_2D_cells  = new_timer  ("sync_3_2_2D_cells")
!     timer_sync_3_2_2D_edges  = new_timer  ("sync_3_2_2D_edges")
!     timer_sync_3_2_2D_verts  = new_timer  ("sync_3_2_2D_verts")
!     timer_sync_3_2_2D_all    = new_timer  ("sync_3_2_2D_all")
!
!     CALL test_sync_2D( SYNC_C, pnt_2D_cells, timer_sync_3_2_2D_cells)
!     CALL test_sync_2D( SYNC_E, pnt_2D_edges, timer_sync_3_2_2D_edges)
!     CALL test_sync_2D( SYNC_V, pnt_2D_verts, timer_sync_3_2_2D_verts)
!     CALL test_sync_2D_all(pnt_2D_cells, pnt_2D_edges, pnt_2D_verts, timer_sync_3_2_2D_all)


    !---------------------------------------------------------------------
    ! test the 2D icon_comm_lib
    !---------------------------------------------------------------------
    timer_iconcom_2D_cells  = new_timer  ("iconcom_2D_cells")
    timer_iconcom_2D_edges  = new_timer  ("iconcom_2D_edges")
    timer_iconcom_2D_verts  = new_timer  ("iconcom_2D_verts")
    timer_iconcom_2D_all    = new_timer  ("iconcom_2D_all")
    timer_iconcom_2D_comb   = new_timer  ("iconcom_2D_comb")
    timer_iconcom_2D_keep   = new_timer  ("iconcom_2D_keep")

    ! test the 2D iconcom on cells

    CALL work_mpi_barrier()
    CALL timer_start(timer_iconcom_2D_cells)
    DO i=1,testbed_iterations
       CALL icon_comm_sync(pnt_2D_cells, p_patch(n_dom_start)%sync_cells_not_in_domain)
    ENDDO
    CALL timer_stop(timer_iconcom_2D_cells)

    ! test the 2D iconcom on edges
    CALL work_mpi_barrier()
    CALL timer_start(timer_iconcom_2D_edges)
    DO i=1,testbed_iterations
       CALL icon_comm_sync(pnt_2D_edges, p_patch(n_dom_start)%sync_edges_not_owned)
    ENDDO
    CALL timer_stop(timer_iconcom_2D_edges)

    ! test the 2D iconcom on verts
    CALL work_mpi_barrier()
    CALL timer_start(timer_iconcom_2D_verts)
    DO i=1,testbed_iterations
       CALL icon_comm_sync(pnt_2D_verts, p_patch(n_dom_start)%sync_verts_not_owned)
    ENDDO
    CALL timer_stop(timer_iconcom_2D_verts)

    ! test the 2D iconcom on all
    CALL work_mpi_barrier()
    CALL timer_start(timer_iconcom_2D_all)
    DO i=1,testbed_iterations
       CALL icon_comm_sync(pnt_2D_cells, p_patch(n_dom_start)%sync_cells_not_in_domain)
       CALL icon_comm_sync(pnt_2D_edges, p_patch(n_dom_start)%sync_edges_not_owned)
       CALL icon_comm_sync(pnt_2D_verts, p_patch(n_dom_start)%sync_verts_not_owned)
    ENDDO
    CALL timer_stop(timer_iconcom_2D_all)

    ! test the 2D iconcom on combined
    CALL work_mpi_barrier()
    CALL timer_start(timer_iconcom_2D_comb)
    DO i=1,testbed_iterations

       comm_1 = new_icon_comm_variable(pnt_2D_cells, &
         &  p_patch(n_dom_start)%sync_cells_not_in_domain, &
         &  status=is_ready, scope=until_sync)

       comm_2 = new_icon_comm_variable(pnt_2D_edges, &
         & p_patch(n_dom_start)%sync_edges_not_owned, status=is_ready, &
         & scope=until_sync)

       comm_3 = new_icon_comm_variable(pnt_2D_verts, &
         & p_patch(n_dom_start)%sync_verts_not_owned, status=is_ready, &
         & scope=until_sync)

       CALL icon_comm_sync_all()

    ENDDO
    CALL timer_stop(timer_iconcom_2D_comb)

    ! test the 2D iconcom on keeping the communicators
    comm_1 = new_icon_comm_variable(pnt_2D_cells, p_patch(n_dom_start)%sync_cells_not_in_domain)
    comm_2 = new_icon_comm_variable(pnt_2D_edges, &
      & p_patch(n_dom_start)%sync_edges_not_owned)
    comm_3 = new_icon_comm_variable(pnt_2D_verts, &
      & p_patch(n_dom_start)%sync_verts_not_owned)

    CALL work_mpi_barrier()
    CALL timer_start(timer_iconcom_2D_keep)
    DO i=1,testbed_iterations

       CALL icon_comm_var_is_ready(comm_1)
       CALL icon_comm_var_is_ready(comm_2)
       CALL icon_comm_var_is_ready(comm_3)
       CALL icon_comm_sync_all()

    ENDDO
    CALL timer_stop(timer_iconcom_2D_keep)
    CALL delete_icon_comm_variable(comm_1)
    CALL delete_icon_comm_variable(comm_2)
    CALL delete_icon_comm_variable(comm_3)
    !---------------------------------------------------------------------
  END SUBROUTINE test_communication_2D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE test_iconcom_3D(comm_pattern, var1, var2, var3, var4, timer_id)
    INTEGER :: comm_pattern, timer_id
    REAL(wp) , POINTER:: var1(:,:,:)
    REAL(wp) , POINTER, OPTIONAL :: var2(:,:,:), var3(:,:,:), var4(:,:,:)

    INTEGER :: i

    CALL work_mpi_barrier()

    IF (.NOT. PRESENT(var2)) THEN
      DO i=1,testbed_iterations
        CALL timer_start(timer_id)
        CALL icon_comm_sync(var1, comm_pattern)
        CALL timer_stop(timer_id)
      ENDDO
      RETURN
    ENDIF

    IF (PRESENT(var4)) THEN
      DO i=1,testbed_iterations
        CALL timer_start(timer_id)
        CALL icon_comm_sync(var1, var2, var3, var4, comm_pattern)
        CALL timer_stop(timer_id)
      ENDDO

   ELSEIF (PRESENT(var3)) THEN
      DO i=1,testbed_iterations
        CALL timer_start(timer_id)
        CALL icon_comm_sync(var1, var2, var3,  comm_pattern)
        CALL timer_stop(timer_id)
      ENDDO
   ELSE
      DO i=1,testbed_iterations
        CALL timer_start(timer_id)
        CALL icon_comm_sync(var1, var2,  comm_pattern)
        CALL timer_stop(timer_id)
      ENDDO
   ENDIF

   CALL work_mpi_barrier()

  END SUBROUTINE test_iconcom_3D
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE test_oldsync_3D(comm_pattern, var1, var2, var3, var4, timer_id)
    INTEGER :: comm_pattern, timer_id
    REAL(wp) , POINTER:: var1(:,:,:)
    REAL(wp) , POINTER, OPTIONAL :: var2(:,:,:), var3(:,:,:), var4(:,:,:)

    INTEGER :: i

    CALL work_mpi_barrier()

    IF (.NOT. PRESENT(var2)) THEN
      DO i=1,testbed_iterations
        CALL timer_start(timer_id)
        CALL sync_patch_array( comm_pattern, p_patch(n_dom_start), var1 )
        CALL timer_stop(timer_id)
      ENDDO
      RETURN
    ENDIF

    IF (PRESENT(var4)) THEN
      DO i=1,testbed_iterations
        CALL timer_start(timer_id)
        CALL sync_patch_array_mult(comm_pattern, p_patch(n_dom_start), 4, var1, var2, var3, var4)
        CALL timer_stop(timer_id)
      ENDDO

   ELSEIF (PRESENT(var3)) THEN
      DO i=1,testbed_iterations
        CALL timer_start(timer_id)
        CALL sync_patch_array_mult(comm_pattern, p_patch(n_dom_start), 3, var1, var2, var3)
        CALL timer_stop(timer_id)
      ENDDO
   ELSE
      DO i=1,testbed_iterations
        CALL timer_start(timer_id)
        CALL sync_patch_array_mult(comm_pattern, p_patch(n_dom_start), 2, var1, var2)
        CALL timer_stop(timer_id)
      ENDDO
   ENDIF

   CALL work_mpi_barrier()

  END SUBROUTINE test_oldsync_3D
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE test_oldsync_4DE1(comm_pattern, var, nfields, timer_id)
    INTEGER :: comm_pattern, nfields, timer_id
    REAL(wp) :: var(:,:,:,:)

    INTEGER :: i

    CALL work_mpi_barrier()
    DO i=1,testbed_iterations
      CALL timer_start(timer_id)
      CALL sync_patch_array_4de1(comm_pattern, p_patch(n_dom_start), nfields, var)
      CALL timer_stop(timer_id)
    ENDDO

   CALL work_mpi_barrier()

  END SUBROUTINE test_oldsync_4DE1
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE test_sync_2D(comm_pattern, var, timer)
    INTEGER :: comm_pattern, timer
    REAL(wp) , POINTER:: var(:,:)

    INTEGER :: i

    CALL timer_start(timer)
    DO i=1,testbed_iterations
      CALL sync_patch_array( comm_pattern, p_patch(n_dom_start), var )
      CALL do_calculations()
    ENDDO
    CALL timer_stop(timer)
  END SUBROUTINE test_sync_2D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE test_sync_2D_all(cell_var, edge_var, vert_var, timer)
    INTEGER :: comm_pattern, timer
    REAL(wp) , POINTER:: cell_var(:,:), edge_var(:,:), vert_var(:,:)

    INTEGER :: i

    CALL timer_start(timer)
    DO i=1,testbed_iterations
      CALL sync_patch_array(SYNC_C , p_patch(n_dom_start), cell_var )
      CALL sync_patch_array(SYNC_E , p_patch(n_dom_start), edge_var )
      CALL sync_patch_array(SYNC_V , p_patch(n_dom_start), vert_var )
    ENDDO
    CALL timer_stop(timer)
  END SUBROUTINE test_sync_2D_all
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE test_sync_3D(comm_pattern, var, timer)
    INTEGER :: comm_pattern, timer
    REAL(wp) , POINTER:: var(:,:,:)

    INTEGER :: i

    CALL timer_start(timer)
    DO i=1,testbed_iterations
      CALL sync_patch_array( comm_pattern, p_patch(n_dom_start), var )
      CALL do_calculations()
    ENDDO
    CALL timer_stop(timer)
  END SUBROUTINE test_sync_3D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE test_sync_3D_all(cell_var, edge_var, vert_var, timer)

    REAL(wp) , POINTER:: cell_var(:,:,:), edge_var(:,:,:), vert_var(:,:,:)

    INTEGER :: i, timer

    CALL timer_start(timer)
    DO i=1,testbed_iterations
      CALL sync_patch_array(SYNC_C , p_patch(n_dom_start), cell_var )
      CALL do_calculations()
      CALL sync_patch_array(SYNC_E , p_patch(n_dom_start), edge_var )
      CALL do_calculations()
      CALL sync_patch_array(SYNC_V , p_patch(n_dom_start), vert_var )
      CALL do_calculations()
    ENDDO
    CALL timer_stop(timer)
  END SUBROUTINE test_sync_3D_all
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE do_calculations()

    INTEGER :: i

    DO i=1,calculate_iterations
      CALL grad_fd_norm (p_nh_state(n_dom_start)%prog(1)%theta_v(:,:,:), &
        & p_patch(n_dom_start), p_nh_state(n_dom_start)%prog(1)%vn(:,:,:))
    ENDDO

  END SUBROUTINE do_calculations
  !-------------------------------------------------------------------------
!     pnt_3D_verts_1 => p_hydro_state(n_dom_start)%diag%rel_vort(:,:,:)
!
!
!     !---------------------------------------------------------------------
!     ! Call cmmunication methods
!   ! INTEGER :: itype_comm = 1
!   ! 1 = synchronous communication with local memory for exchange buffers
!   ! 2 = synchronous communication with global memory for exchange buffers
!   ! 3 = asynchronous communication within dynamical core with global memory
!   !     for exchange buffers (not yet implemented)
!
!   ! Order of send/receive sequence in exchange routines
!  ! INTEGER :: iorder_sendrecv = 1
!   ! 1 = irecv, send
!   ! 2 = isend, recv
!
!
!
!     timer_sync_1_1_3D_edges  = new_timer  ("sync_1_1_3D_edges")
!     timer_sync_1_1_3D_verts  = new_timer  ("sync_1_1_3D_verts")
!     timer_sync_1_1_3D_all    = new_timer  ("sync_1_1_3D_all")
!
!     CALL work_mpi_barrier()
!     CALL test_sync_3D( SYNC_C, pnt_3D_cells_1, timer_sync_1_1_3D_cells)
!     CALL work_mpi_barrier()
!     CALL test_sync_3D( SYNC_E, pnt_3D_edges_1, timer_sync_1_1_3D_edges)
!     CALL work_mpi_barrier()
!     CALL test_sync_3D( SYNC_V, pnt_3D_verts_1, timer_sync_1_1_3D_verts)
!     CALL work_mpi_barrier()
!     CALL test_sync_3D_all(pnt_3D_cells_1, pnt_3D_edges_1, pnt_3D_verts_1, timer_sync_1_1_3D_all)
!
!     !---------------------------------------------------------------------
!     itype_comm = 1
!     iorder_sendrecv = 2
!     timer_sync_1_2_3D_cells  = new_timer  ("sync_1_2_3D_cells")
!     timer_sync_1_2_3D_edges  = new_timer  ("sync_1_2_3D_edges")
!     timer_sync_1_2_3D_verts  = new_timer  ("sync_1_2_3D_verts")
!     timer_sync_1_2_3D_all    = new_timer  ("sync_1_2_3D_all")
!
!     CALL work_mpi_barrier()
!     CALL test_sync_3D( SYNC_C, pnt_3D_cells_1, timer_sync_1_2_3D_cells)
!     CALL work_mpi_barrier()
!     CALL test_sync_3D( SYNC_E, pnt_3D_edges_1, timer_sync_1_2_3D_edges)
!     CALL work_mpi_barrier()
!     CALL test_sync_3D( SYNC_V, pnt_3D_verts_1, timer_sync_1_2_3D_verts)
!     CALL work_mpi_barrier()
!     CALL test_sync_3D_all(pnt_3D_cells_1, pnt_3D_edges_1, pnt_3D_verts_1, timer_sync_1_2_3D_all)
!
!     !---------------------------------------------------------------------
! !     itype_comm = 2
! !     iorder_sendrecv = 1
! !     timer_sync_2_1_3D_cells  = new_timer  ("sync_2_1_3D_cells")
! !     timer_sync_2_1_3D_edges  = new_timer  ("sync_2_1_3D_edges")
! !     timer_sync_2_1_3D_verts  = new_timer  ("sync_2_1_3D_verts")
! !     timer_sync_2_1_3D_all    = new_timer  ("sync_2_1_3D_all")
! !
! !     CALL test_sync_3D( SYNC_C, pnt_3D_cells_1, timer_sync_2_1_3D_cells)
! !     CALL test_sync_3D( SYNC_E, pnt_3D_edges_1, timer_sync_2_1_3D_edges)
! !     CALL test_sync_3D( SYNC_V, pnt_3D_verts_1, timer_sync_2_1_3D_verts)
! !     CALL test_sync_3D_all(pnt_3D_cells_1, pnt_3D_edges_1, pnt_3D_verts_1, timer_sync_2_1_3D_all)
! !
! !     !---------------------------------------------------------------------
! !     itype_comm = 2
! !     iorder_sendrecv = 2
! !     timer_sync_2_2_3D_cells  = new_timer  ("sync_2_2_3D_cells")
! !     timer_sync_2_2_3D_edges  = new_timer  ("sync_2_2_3D_edges")
! !     timer_sync_2_2_3D_verts  = new_timer  ("sync_2_2_3D_verts")
! !     timer_sync_2_2_3D_all    = new_timer  ("sync_2_2_3D_all")
! !
! !     CALL test_sync_3D( SYNC_C, pnt_3D_cells_1, timer_sync_2_2_3D_cells)
! !     CALL test_sync_3D( SYNC_E, pnt_3D_edges_1, timer_sync_2_2_3D_edges)
! !     CALL test_sync_3D( SYNC_V, pnt_3D_verts_1, timer_sync_2_2_3D_verts)
! !     CALL test_sync_3D_all(pnt_3D_cells_1, pnt_3D_edges_1, pnt_3D_verts_1, timer_sync_2_2_3D_all)
! !
! !     !---------------------------------------------------------------------
! !     itype_comm = 3
! !     iorder_sendrecv = 1
! !     timer_sync_3_1_3D_cells  = new_timer  ("sync_3_1_3D_cells")
! !     timer_sync_3_1_3D_edges  = new_timer  ("sync_3_1_3D_edges")
! !     timer_sync_3_1_3D_verts  = new_timer  ("sync_3_1_3D_verts")
! !     timer_sync_3_1_3D_all    = new_timer  ("sync_3_1_3D_all")
! !
! !     CALL test_sync_3D( SYNC_C, pnt_3D_cells_1, timer_sync_3_1_3D_cells)
! !     CALL test_sync_3D( SYNC_E, pnt_3D_edges_1, timer_sync_3_1_3D_edges)
! !     CALL test_sync_3D( SYNC_V, pnt_3D_verts_1, timer_sync_3_1_3D_verts)
! !     CALL test_sync_3D_all(pnt_3D_cells_1, pnt_3D_edges_1, pnt_3D_verts_1, timer_sync_3_1_3D_all)
! !
! !     !---------------------------------------------------------------------
! !     itype_comm = 3
! !     iorder_sendrecv = 2
! !     timer_sync_3_2_3D_cells  = new_timer  ("sync_3_2_3D_cells")
! !     timer_sync_3_2_3D_edges  = new_timer  ("sync_3_2_3D_edges")
! !     timer_sync_3_2_3D_verts  = new_timer  ("sync_3_2_3D_verts")
! !     timer_sync_3_2_3D_all    = new_timer  ("sync_3_2_3D_all")
! !
! !     CALL test_sync_3D( SYNC_C, pnt_3D_cells_1, timer_sync_3_2_3D_cells)
! !     CALL test_sync_3D( SYNC_E, pnt_3D_edges_1, timer_sync_3_2_3D_edges)
! !     CALL test_sync_3D( SYNC_V, pnt_3D_verts_1, timer_sync_3_2_3D_verts)
! !     CALL test_sync_3D_all(pnt_3D_cells_1, pnt_3D_edges_1, pnt_3D_verts_1, timer_sync_3_2_3D_all)
!
!
!
!     !---------------------------------------------------------------------
!     ! test the 3D icon_comm_lib
!     !---------------------------------------------------------------------
!     timer_iconcom_3D_cells  = new_timer  ("iconcom_3D_cells")
!     timer_iconcom_3D_edges  = new_timer  ("iconcom_3D_edges")
!     timer_iconcom_3D_verts  = new_timer  ("iconcom_3D_verts")
!     timer_iconcom_3D_all    = new_timer  ("iconcom_3D_all")
!     timer_iconcom_3D_comb   = new_timer  ("iconcom_3D_comb")
!     timer_iconcom_3D_keep   = new_timer  ("iconcom_3D_keep")
!
!     ! test the 3D iconcom on cells
!     CALL work_mpi_barrier()
!     CALL timer_start(timer_iconcom_3D_cells)
!     DO i=1,testbed_iterations
!       CALL icon_comm_sync(pnt_3D_cells_1, cells_not_in_domain)
!       CALL do_calculations()
!     ENDDO
!     CALL timer_stop(timer_iconcom_3D_cells)
!
!     CALL work_mpi_barrier()
!
!     ! test the 3D iconcom on edges
!     CALL work_mpi_barrier()
!     CALL timer_start(timer_iconcom_3D_edges)
!     DO i=1,testbed_iterations
!       CALL icon_comm_sync(pnt_3D_edges_1, edges_not_owned)
!       CALL do_calculations()
!     ENDDO
!     CALL timer_stop(timer_iconcom_3D_edges)
!
!     ! test the 3D iconcom on verts
!     CALL timer_start(timer_iconcom_3D_verts)
!     DO i=1,testbed_iterations
!       CALL icon_comm_sync(pnt_3D_verts_1, verts_not_owned)
!       CALL do_calculations()
!     ENDDO
!     CALL timer_stop(timer_iconcom_3D_verts)
!
!     ! test the 3D iconcom on all
!     CALL work_mpi_barrier()
!     CALL timer_start(timer_iconcom_3D_all)
!     DO i=1,testbed_iterations
!        CALL icon_comm_sync(pnt_3D_cells_1, cells_not_in_domain)
!        CALL icon_comm_sync(pnt_3D_edges_1, edges_not_owned)
!        CALL icon_comm_sync(pnt_3D_verts_1, verts_not_owned)
!     ENDDO
!     CALL timer_stop(timer_iconcom_3D_all)
!
!     ! test the 3D iconcom on combined
!     CALL work_mpi_barrier()
!     CALL timer_start(timer_iconcom_3D_comb)
!     DO i=1,testbed_iterations
!
!       comm_1 = new_icon_comm_variable(pnt_3D_cells_1, cells_not_in_domain, &
!         & p_patch(n_dom_start), status=is_ready, scope=until_sync)
!
!       comm_2 = new_icon_comm_variable(pnt_3D_edges_1, &
!         & edges_not_owned, p_patch(n_dom_start), status=is_ready, scope=until_sync)
!
!       comm_3 = new_icon_comm_variable(pnt_3D_verts_1, &
!         & verts_not_owned, p_patch(n_dom_start), status=is_ready, scope=until_sync)
!
!       CALL icon_comm_sync_all()
!
!       CALL do_calculations()
!       CALL do_calculations()
!       CALL do_calculations()
!
!     ENDDO
!     CALL timer_stop(timer_iconcom_3D_comb)
!
!     ! test the 3D iconcom on keeping the communicators
!     comm_1 = new_icon_comm_variable(pnt_3D_cells_1, cells_not_in_domain, &
!       & p_patch(n_dom_start))
!     comm_2 = new_icon_comm_variable(pnt_3D_edges_1, &
!       & edges_not_owned)
!     comm_3 = new_icon_comm_variable(pnt_3D_verts_1, &
!       & verts_not_owned)
!
!     CALL work_mpi_barrier()
!     CALL timer_start(timer_iconcom_3D_keep)
!     DO i=1,testbed_iterations
!
!       CALL icon_comm_var_is_ready(comm_1)
!       CALL icon_comm_var_is_ready(comm_2)
!       CALL icon_comm_var_is_ready(comm_3)
!       CALL icon_comm_sync_all()
!
!       CALL do_calculations()
!       CALL do_calculations()
!       CALL do_calculations()
!
!     ENDDO
!     CALL timer_stop(timer_iconcom_3D_keep)
!     CALL delete_icon_comm_variable(comm_3)
!     CALL delete_icon_comm_variable(comm_2)
!     CALL delete_icon_comm_variable(comm_1)
    !---------------------------------------------------------------------

  SUBROUTINE exchange_communication_testbed(test_gpu)

    LOGICAL, OPTIONAL, INTENT(IN) :: test_gpu

    INTEGER, ALLOCATABLE :: owner_local_src(:), owner_local_dst(:), &
      &                     glb_index_src(:), glb_index_dst(:)
    INTEGER :: i, j, n, local_size_src, local_size_dst, global_size
    TYPE(t_glb2loc_index_lookup) :: send_glb2loc_index
    CLASS(t_comm_pattern), POINTER :: comm_pattern

    INTEGER :: nlev
    REAL(wp), ALLOCATABLE :: in_array_r_2d(:,:), in_array_r_3d(:,:,:)
    INTEGER, ALLOCATABLE ::  in_array_i_2d(:,:), in_array_i_3d(:,:,:)
    LOGICAL, ALLOCATABLE ::  in_array_l_2d(:,:), in_array_l_3d(:,:,:)
    REAL(wp), ALLOCATABLE :: out_array_r_2d(:,:), out_array_r_3d(:,:,:)
    INTEGER, ALLOCATABLE ::  out_array_i_2d(:,:), out_array_i_3d(:,:,:)
    LOGICAL, ALLOCATABLE ::  out_array_l_2d(:,:), out_array_l_3d(:,:,:)
    REAL(wp), ALLOCATABLE :: add_array_r_2d(:,:), add_array_r_3d(:,:,:)
    INTEGER, ALLOCATABLE ::  add_array_i_2d(:,:), add_array_i_3d(:,:,:)
    REAL(wp), ALLOCATABLE :: ref_out_array_r_2d(:,:), ref_out_array_r_3d(:,:,:)
    INTEGER, ALLOCATABLE ::  ref_out_array_i_2d(:,:), ref_out_array_i_3d(:,:,:)
    LOGICAL, ALLOCATABLE ::  ref_out_array_l_2d(:,:), ref_out_array_l_3d(:,:,:)
    REAL(dp), ALLOCATABLE :: in_array_dp_4d(:,:,:,:), &
      &                      out_array_dp_4d(:,:,:,:), &
      &                      add_array_dp_4d(:,:,:,:), &
      &                      ref_out_array_dp_4d(:,:,:,:)
    REAL(sp), ALLOCATABLE :: in_array_sp_4d(:,:,:,:), &
      &                      out_array_sp_4d(:,:,:,:), &
      &                      add_array_sp_4d(:,:,:,:), &
      &                      ref_out_array_sp_4d(:,:,:,:)

    CHARACTER(*), PARAMETER :: method_name = &
      "mo_test_communication:exchange_communication_testbed"

    ! generate communication pattern
    local_size_src = 10 * nproma
    local_size_dst = 16 * nproma
    global_size = local_size_src * p_n_work
    ALLOCATE(owner_local_src(local_size_src), glb_index_src(local_size_src), &
      &      owner_local_dst(local_size_dst), glb_index_dst(local_size_dst))

    owner_local_src(1:local_size_src:2) = -1
    owner_local_src(2:local_size_src:2) = p_pe_work

    owner_local_dst(1:local_size_dst:2) = -1
    owner_local_dst(2:3*nproma:2) = MOD(p_pe_work + p_n_work - 1, p_n_work)
    owner_local_dst(3*nproma+2:13*nproma:2) = p_pe_work
    owner_local_dst(13*nproma+2:16*nproma:2) = MOD(p_pe_work + 1, p_n_work)
    DO i = 1, local_size_src
      glb_index_src(i) = local_size_src * p_pe_work + i
    END DO
    DO i = 1, local_size_dst
      glb_index_dst(i) = MOD(local_size_src * p_pe_work - 3 * nproma - 1 + i + &
      &                      global_size, global_size) + 1
    END DO
    CALL init_glb2loc_index_lookup(send_glb2loc_index, global_size)

    CALL set_inner_glb_index(send_glb2loc_index, glb_index_src, &
      &                      (/(i, i=1, local_size_src)/))

    CALL setup_comm_pattern(local_size_dst, owner_local_dst, glb_index_dst, &
      &                     send_glb2loc_index, local_size_src, &
      &                     owner_local_src, glb_index_src, comm_pattern)

    nlev = 7

    ALLOCATE(in_array_r_2d(nproma, 10), in_array_r_3d(nproma, nlev, 10), &
      &      in_array_i_2d(nproma, 10), in_array_i_3d(nproma, nlev, 10), &
      &      in_array_l_2d(nproma, 10), in_array_l_3d(nproma, nlev, 10), &
      &      out_array_r_2d(nproma, 16), out_array_r_3d(nproma, nlev, 16), &
      &      out_array_i_2d(nproma, 16), out_array_i_3d(nproma, nlev, 16), &
      &      out_array_l_2d(nproma, 16), out_array_l_3d(nproma, nlev, 16), &
      &      add_array_r_2d(nproma, 16), add_array_r_3d(nproma, nlev, 16), &
      &      add_array_i_2d(nproma, 16), add_array_i_3d(nproma, nlev, 16), &
      &      ref_out_array_r_2d(nproma, 16), &
      &      ref_out_array_r_3d(nproma, nlev, 16), &
      &      ref_out_array_i_2d(nproma, 16), &
      &      ref_out_array_i_3d(nproma, nlev, 16), &
      &      ref_out_array_l_2d(nproma, 16), &
      &      ref_out_array_l_3d(nproma, nlev, 16))

    in_array_r_2d = RESHAPE(glb_index_src, (/nproma, 10/))
    in_array_i_2d = RESHAPE(glb_index_src, (/nproma, 10/))
    in_array_l_2d = .TRUE.
    DO i = 1, nlev
      in_array_r_3d(:,i,:) = in_array_r_2d + (i - 1) * global_size
      in_array_i_3d(:,i,:) = in_array_i_2d + (i - 1) * global_size
    END DO
    in_array_l_3d = .TRUE.

    out_array_r_2d = -1
    out_array_i_2d = -1
    out_array_l_2d = .FALSE.
    out_array_r_3d = -1
    out_array_i_3d = -1
    out_array_l_3d = .FALSE.

    ref_out_array_r_2d = RESHAPE(MERGE(-1, glb_index_dst, &
      &                                owner_local_dst == -1), (/nproma, 16/))
    ref_out_array_i_2d = RESHAPE(MERGE(-1, glb_index_dst, &
      &                                owner_local_dst == -1), (/nproma, 16/))
    ref_out_array_l_2d = RESHAPE(MERGE(.FALSE., .TRUE., owner_local_dst == -1), &
      &                          (/nproma, 16/))
    DO i = 1, nlev
      ref_out_array_r_3d(:,i,:) = MERGE(-1._wp, ref_out_array_r_2d + &
        &                               (i - 1) * global_size, &
        &                               -1 == RESHAPE(owner_local_dst, &
        &                                             (/nproma, 16/)))
      ref_out_array_i_3d(:,i,:) = MERGE(-1, ref_out_array_i_2d + &
        &                               (i - 1) * global_size, &
        &                               -1 == RESHAPE(owner_local_dst, &
        &                                             (/nproma, 16/)))
      ref_out_array_l_3d(:,i,:) = ref_out_array_l_2d
    END DO

    CALL check_exchange(in_array_r_2d = in_array_r_2d, &
      &                 in_array_r_3d = in_array_r_3d, &
      &                 in_array_i_2d = in_array_i_2d, &
      &                 in_array_i_3d = in_array_i_3d, &
      &                 in_array_l_2d = in_array_l_2d, &
      &                 in_array_l_3d = in_array_l_3d, &
      &                 out_array_r_2d = out_array_r_2d, &
      &                 out_array_r_3d = out_array_r_3d, &
      &                 out_array_i_2d = out_array_i_2d, &
      &                 out_array_i_3d = out_array_i_3d, &
      &                 out_array_l_2d = out_array_l_2d, &
      &                 out_array_l_3d = out_array_l_3d, &
      &                 ref_out_array_r_2d = ref_out_array_r_2d, &
      &                 ref_out_array_r_3d = ref_out_array_r_3d, &
      &                 ref_out_array_i_2d = ref_out_array_i_2d, &
      &                 ref_out_array_i_3d = ref_out_array_i_3d, &
      &                 ref_out_array_l_2d = ref_out_array_l_2d, &
      &                 ref_out_array_l_3d = ref_out_array_l_3d, &
      &                 comm_pattern = comm_pattern)

    add_array_r_2d = MERGE(-1._wp, ref_out_array_r_2d, &
      &                    -1 == RESHAPE(owner_local_dst, (/nproma, 16/)))
    add_array_i_2d = MERGE(-1, ref_out_array_i_2d, &
      &                    -1 == RESHAPE(owner_local_dst, (/nproma, 16/)))
    DO i = 1, nlev
      add_array_r_3d(:,i,:) = MERGE(-1._wp, ref_out_array_r_3d(:,i,:), &
        &                           -1 == RESHAPE(owner_local_dst, (/nproma, 16/)))
      add_array_i_3d(:,i,:) = MERGE(-1, ref_out_array_i_3d(:,i,:), &
        &                           -1 == RESHAPE(owner_local_dst, (/nproma, 16/)))
    END DO

    ref_out_array_r_2d = MERGE(-1._wp, 2._wp * ref_out_array_r_2d, &
      &                        -1 == RESHAPE(owner_local_dst, (/nproma, 16/)))
    ref_out_array_i_2d = MERGE(-1, 2 * ref_out_array_i_2d, &
      &                    -1 == RESHAPE(owner_local_dst, (/nproma, 16/)))
    DO i = 1, nlev
      ref_out_array_r_3d(:,i,:) = MERGE(-1._wp, 2._wp * ref_out_array_r_3d(:,i,:), &
        &                               -1 == RESHAPE(owner_local_dst, &
        &                               (/nproma, 16/)))
      ref_out_array_i_3d(:,i,:) = MERGE(-1, 2 * ref_out_array_i_3d(:,i,:), &
        &                               -1 == RESHAPE(owner_local_dst, &
        &                               (/nproma, 16/)))
    END DO

    CALL check_exchange(in_array_r_2d = in_array_r_2d, &
      &                 in_array_r_3d = in_array_r_3d, &
      &                 in_array_i_2d = in_array_i_2d, &
      &                 in_array_i_3d = in_array_i_3d, &
      &                 in_array_l_2d = in_array_l_2d, &
      &                 in_array_l_3d = in_array_l_3d, &
      &                 out_array_r_2d = out_array_r_2d, &
      &                 out_array_r_3d = out_array_r_3d, &
      &                 out_array_i_2d = out_array_i_2d, &
      &                 out_array_i_3d = out_array_i_3d, &
      &                 out_array_l_2d = out_array_l_2d, &
      &                 out_array_l_3d = out_array_l_3d, &
      &                 add_array_r_2d = add_array_r_2d, &
      &                 add_array_r_3d = add_array_r_3d, &
      &                 add_array_i_2d = add_array_i_2d, &
      &                 add_array_i_3d = add_array_i_3d, &
      &                 ref_out_array_r_2d = ref_out_array_r_2d, &
      &                 ref_out_array_r_3d = ref_out_array_r_3d, &
      &                 ref_out_array_i_2d = ref_out_array_i_2d, &
      &                 ref_out_array_i_3d = ref_out_array_i_3d, &
      &                 ref_out_array_l_2d = ref_out_array_l_2d, &
      &                 ref_out_array_l_3d = ref_out_array_l_3d, &
      &                 comm_pattern = comm_pattern)

    CALL delete_comm_pattern(comm_pattern)

    ! generate communication pattern
    local_size_dst = 10 * nproma
    DEALLOCATE(owner_local_dst, glb_index_dst)
    ALLOCATE(owner_local_dst(local_size_dst), glb_index_dst(local_size_dst))

    owner_local_dst(1:local_size_dst:2) = -1
    owner_local_dst(2:7*nproma:2) = p_pe_work
    owner_local_dst(7*nproma+2:local_size_dst:2) = &
      MOD(p_pe_work + p_n_work + 1, p_n_work)
    DO i = 1, local_size_dst
      glb_index_dst(i) = MOD(local_size_src * p_pe_work + 3 * nproma - 1 + i + &
      &                      global_size, global_size) + 1
    END DO

    CALL setup_comm_pattern(local_size_dst, owner_local_dst, glb_index_dst, &
      &                     send_glb2loc_index, local_size_src, &
      &                     owner_local_src, glb_index_src, comm_pattern)

    DEALLOCATE(out_array_r_2d, out_array_r_3d, out_array_i_2d, out_array_i_3d, &
      &        out_array_l_2d, out_array_l_3d, add_array_r_2d, add_array_r_3d, &
      &        add_array_i_2d, add_array_i_3d, ref_out_array_r_2d, &
      &        ref_out_array_r_3d, ref_out_array_i_2d, ref_out_array_i_3d, &
      &        ref_out_array_l_2d, ref_out_array_l_3d)
    ALLOCATE(out_array_r_2d(nproma, 10), out_array_r_3d(nproma, nlev, 10), &
      &      out_array_i_2d(nproma, 10), out_array_i_3d(nproma, nlev, 10), &
      &      out_array_l_2d(nproma, 10), out_array_l_3d(nproma, nlev, 10), &
      &      add_array_r_2d(nproma, 10), add_array_r_3d(nproma, nlev, 10), &
      &      add_array_i_2d(nproma, 10), add_array_i_3d(nproma, nlev, 10), &
      &      ref_out_array_r_2d(nproma, 10), &
      &      ref_out_array_r_3d(nproma, nlev, 10), &
      &      ref_out_array_i_2d(nproma, 10), &
      &      ref_out_array_i_3d(nproma, nlev, 10), &
      &      ref_out_array_l_2d(nproma, 10), &
      &      ref_out_array_l_3d(nproma, nlev, 10))

    out_array_r_2d = -1
    out_array_i_2d = -1
    out_array_l_2d = .FALSE.
    out_array_r_3d = -1
    out_array_i_3d = -1
    out_array_l_3d = .FALSE.

    ref_out_array_r_2d = RESHAPE(MERGE(-1, glb_index_dst, &
      &                                owner_local_dst == -1), (/nproma, 10/))
    ref_out_array_i_2d = RESHAPE(MERGE(-1, glb_index_dst, &
      &                                owner_local_dst == -1), (/nproma, 10/))
    ref_out_array_l_2d = RESHAPE(MERGE(.FALSE., .TRUE., owner_local_dst == -1), &
      &                          (/nproma, 10/))
    DO i = 1, nlev
      ref_out_array_r_3d(:,i,:) = MERGE(-1._wp, ref_out_array_r_2d + &
        &                               (i - 1) * global_size, &
        &                               -1 == RESHAPE(owner_local_dst, &
        &                                             (/nproma, 10/)))
      ref_out_array_i_3d(:,i,:) = MERGE(-1, ref_out_array_i_2d + &
        &                               (i - 1) * global_size, &
        &                               -1 == RESHAPE(owner_local_dst, &
        &                                             (/nproma, 10/)))
      ref_out_array_l_3d(:,i,:) = ref_out_array_l_2d
    END DO

    CALL check_exchange(in_array_r_2d = in_array_r_2d, &
      &                 in_array_r_3d = in_array_r_3d, &
      &                 in_array_i_2d = in_array_i_2d, &
      &                 in_array_i_3d = in_array_i_3d, &
      &                 in_array_l_2d = in_array_l_2d, &
      &                 in_array_l_3d = in_array_l_3d, &
      &                 out_array_r_2d = out_array_r_2d, &
      &                 out_array_r_3d = out_array_r_3d, &
      &                 out_array_i_2d = out_array_i_2d, &
      &                 out_array_i_3d = out_array_i_3d, &
      &                 out_array_l_2d = out_array_l_2d, &
      &                 out_array_l_3d = out_array_l_3d, &
      &                 ref_out_array_r_2d = ref_out_array_r_2d, &
      &                 ref_out_array_r_3d = ref_out_array_r_3d, &
      &                 ref_out_array_i_2d = ref_out_array_i_2d, &
      &                 ref_out_array_i_3d = ref_out_array_i_3d, &
      &                 ref_out_array_l_2d = ref_out_array_l_2d, &
      &                 ref_out_array_l_3d = ref_out_array_l_3d, &
      &                 comm_pattern = comm_pattern)

    add_array_r_2d = MERGE(-1._wp, ref_out_array_r_2d, &
      &                    -1 == RESHAPE(owner_local_dst, (/nproma, 10/)))
    add_array_i_2d = MERGE(-1, ref_out_array_i_2d, &
      &                    -1 == RESHAPE(owner_local_dst, (/nproma, 10/)))
    DO i = 1, nlev
      add_array_r_3d(:,i,:) = MERGE(-1._wp, ref_out_array_r_3d(:,i,:), &
        &                           -1 == RESHAPE(owner_local_dst, (/nproma, 10/)))
      add_array_i_3d(:,i,:) = MERGE(-1, ref_out_array_i_3d(:,i,:), &
        &                           -1 == RESHAPE(owner_local_dst, (/nproma, 10/)))
    END DO

    ref_out_array_r_2d = MERGE(-1._wp, 2._wp * ref_out_array_r_2d, &
      &                        -1 == RESHAPE(owner_local_dst, (/nproma, 10/)))
    ref_out_array_i_2d = MERGE(-1, 2 * ref_out_array_i_2d, &
      &                    -1 == RESHAPE(owner_local_dst, (/nproma, 10/)))
    DO i = 1, nlev
      ref_out_array_r_3d(:,i,:) = MERGE(-1._wp, 2._wp * ref_out_array_r_3d(:,i,:), &
        &                               -1 == RESHAPE(owner_local_dst, &
        &                               (/nproma, 10/)))
      ref_out_array_i_3d(:,i,:) = MERGE(-1, 2 * ref_out_array_i_3d(:,i,:), &
        &                               -1 == RESHAPE(owner_local_dst, &
        &                               (/nproma, 10/)))
    END DO

    CALL check_exchange(in_array_r_2d = in_array_r_2d, &
      &                 in_array_r_3d = in_array_r_3d, &
      &                 in_array_i_2d = in_array_i_2d, &
      &                 in_array_i_3d = in_array_i_3d, &
      &                 in_array_l_2d = in_array_l_2d, &
      &                 in_array_l_3d = in_array_l_3d, &
      &                 out_array_r_2d = out_array_r_2d, &
      &                 out_array_r_3d = out_array_r_3d, &
      &                 out_array_i_2d = out_array_i_2d, &
      &                 out_array_i_3d = out_array_i_3d, &
      &                 out_array_l_2d = out_array_l_2d, &
      &                 out_array_l_3d = out_array_l_3d, &
      &                 add_array_r_2d = add_array_r_2d, &
      &                 add_array_r_3d = add_array_r_3d, &
      &                 add_array_i_2d = add_array_i_2d, &
      &                 add_array_i_3d = add_array_i_3d, &
      &                 ref_out_array_r_2d = ref_out_array_r_2d, &
      &                 ref_out_array_r_3d = ref_out_array_r_3d, &
      &                 ref_out_array_i_2d = ref_out_array_i_2d, &
      &                 ref_out_array_i_3d = ref_out_array_i_3d, &
      &                 ref_out_array_l_2d = ref_out_array_l_2d, &
      &                 ref_out_array_l_3d = ref_out_array_l_3d, &
      &                 comm_pattern = comm_pattern)

    out_array_r_2d = in_array_r_2d
    out_array_i_2d = in_array_i_2d
    out_array_l_2d = in_array_l_2d
    out_array_r_3d = in_array_r_3d
    out_array_i_3d = in_array_i_3d
    out_array_l_3d = in_array_l_3d

    ref_out_array_r_2d = RESHAPE(MERGE(glb_index_src, glb_index_dst, &
      &                                owner_local_dst == -1), (/nproma, 10/))
    ref_out_array_i_2d = RESHAPE(MERGE(glb_index_src, glb_index_dst, &
      &                                owner_local_dst == -1), (/nproma, 10/))
    ref_out_array_l_2d = .TRUE.
    DO i = 1, nlev
      ref_out_array_r_3d(:,i,:) = ref_out_array_r_2d + (i - 1) * global_size
      ref_out_array_i_3d(:,i,:) = ref_out_array_i_2d + (i - 1) * global_size
    END DO
    ref_out_array_l_3d = .TRUE.

    CALL check_exchange(out_array_r_2d = out_array_r_2d, &
      &                 out_array_r_3d = out_array_r_3d, &
      &                 out_array_i_2d = out_array_i_2d, &
      &                 out_array_i_3d = out_array_i_3d, &
      &                 out_array_l_2d = out_array_l_2d, &
      &                 out_array_l_3d = out_array_l_3d, &
      &                 ref_out_array_r_2d = ref_out_array_r_2d, &
      &                 ref_out_array_r_3d = ref_out_array_r_3d, &
      &                 ref_out_array_i_2d = ref_out_array_i_2d, &
      &                 ref_out_array_i_3d = ref_out_array_i_3d, &
      &                 ref_out_array_l_2d = ref_out_array_l_2d, &
      &                 ref_out_array_l_3d = ref_out_array_l_3d, &
      &                 comm_pattern = comm_pattern)

    out_array_r_2d = in_array_r_2d
    out_array_i_2d = in_array_i_2d
    out_array_l_2d = in_array_l_2d
    out_array_r_3d = in_array_r_3d
    out_array_i_3d = in_array_i_3d
    out_array_l_3d = in_array_l_3d

    add_array_r_2d = MERGE(-1._wp, ref_out_array_r_2d, &
      &                    -1 == RESHAPE(owner_local_dst, (/nproma, 10/)))
    add_array_i_2d = MERGE(-1, ref_out_array_i_2d, &
      &                    -1 == RESHAPE(owner_local_dst, (/nproma, 10/)))
    DO i = 1, nlev
      add_array_r_3d(:,i,:) = MERGE(-1._wp, ref_out_array_r_3d(:,i,:), &
        &                           -1 == RESHAPE(owner_local_dst, (/nproma, 10/)))
      add_array_i_3d(:,i,:) = MERGE(-1, ref_out_array_i_3d(:,i,:), &
        &                           -1 == RESHAPE(owner_local_dst, (/nproma, 10/)))
    END DO

    ref_out_array_r_2d = RESHAPE(MERGE(glb_index_src, 2 * glb_index_dst, &
      &                                owner_local_dst == -1), (/nproma, 10/))
    ref_out_array_i_2d = RESHAPE(MERGE(glb_index_src, 2 * glb_index_dst, &
      &                                owner_local_dst == -1), (/nproma, 10/))
    DO i = 1, nlev
      ref_out_array_r_3d(:,i,:) = MERGE(ref_out_array_r_3d(:,i,:), &
        &                               2._wp * ref_out_array_r_3d(:,i,:), &
        &                               -1 == RESHAPE(owner_local_dst, &
        &                               (/nproma, 10/)))
      ref_out_array_i_3d(:,i,:) = MERGE(ref_out_array_i_3d(:,i,:), &
        &                               2 * ref_out_array_i_3d(:,i,:), &
        &                               -1 == RESHAPE(owner_local_dst, &
        &                               (/nproma, 10/)))
    END DO

    CALL check_exchange(out_array_r_2d = out_array_r_2d, &
      &                 out_array_r_3d = out_array_r_3d, &
      &                 out_array_i_2d = out_array_i_2d, &
      &                 out_array_i_3d = out_array_i_3d, &
      &                 out_array_l_2d = out_array_l_2d, &
      &                 out_array_l_3d = out_array_l_3d, &
      &                 add_array_r_2d = add_array_r_2d, &
      &                 add_array_r_3d = add_array_r_3d, &
      &                 add_array_i_2d = add_array_i_2d, &
      &                 add_array_i_3d = add_array_i_3d, &
      &                 ref_out_array_r_2d = ref_out_array_r_2d, &
      &                 ref_out_array_r_3d = ref_out_array_r_3d, &
      &                 ref_out_array_i_2d = ref_out_array_i_2d, &
      &                 ref_out_array_i_3d = ref_out_array_i_3d, &
      &                 ref_out_array_l_2d = ref_out_array_l_2d, &
      &                 ref_out_array_l_3d = ref_out_array_l_3d, &
      &                 comm_pattern = comm_pattern)

    DO n = 1, 16

      ALLOCATE(in_array_dp_4d(n, nproma, nlev, 10), &
        &      out_array_dp_4d(n, nproma, nlev, 10), &
        &      ref_out_array_dp_4d(n, nproma, nlev, 10))

      DO i = 1, nlev
        DO j = 1, n
          in_array_dp_4d(j,:,i,:) = RESHAPE(glb_index_src, (/nproma, 10/)) + &
            &                               (i - 1) * global_size + 0.1_wp * j
        END DO
      END DO
      out_array_dp_4d = -1
      DO i = 1, nlev
        DO j = 1, n
          ref_out_array_dp_4d(j,:,i,:) = &
            RESHAPE(MERGE(-1._wp, glb_index_dst + (i - 1) * global_size + &
            &             0.1_wp * j, owner_local_dst == -1), (/nproma, 10/))
        END DO
      END DO

      CALL check_exchange_4de1(in_array=in_array_dp_4d, &
        &                      out_array=out_array_dp_4d, &
        &                      ref_out_array=ref_out_array_dp_4d, &
        &                      comm_pattern=comm_pattern)

      out_array_dp_4d = in_array_dp_4d
      DO i = 1, nlev
        DO j = 1, n
          ref_out_array_dp_4d(j,:,i,:) = &
            RESHAPE(MERGE(glb_index_src, glb_index_dst, &
            &             owner_local_dst == -1), (/nproma, 10/)) + &
            &             (i - 1) * global_size + 0.1_wp * j
        END DO
      END DO

      CALL check_exchange_4de1(out_array=out_array_dp_4d, &
        &                      ref_out_array=ref_out_array_dp_4d, &
        &                      comm_pattern=comm_pattern)

      DEALLOCATE(in_array_dp_4d, out_array_dp_4d, ref_out_array_dp_4d)
    END DO

    ALLOCATE(in_array_dp_4d(nproma,nlev,10,20), &
      &      add_array_dp_4d(nproma,nlev,10,20), &
      &      out_array_dp_4d(nproma,nlev,10,20), &
      &      ref_out_array_dp_4d(nproma,nlev,10,20), &
      &      in_array_sp_4d(nproma,nlev,10,20), &
      &      add_array_sp_4d(nproma,nlev,10,20), &
      &      out_array_sp_4d(nproma,nlev,10,20), &
      &      ref_out_array_sp_4d(nproma,nlev,10,20))

    DO n = 1, 20
      DO i = 1, nlev
        in_array_dp_4d(:,i,:,n) = &
          RESHAPE(MERGE(-1._wp, glb_index_src + (i - 1) * global_size + &
          &             n * 0.1_wp, -1 == owner_local_src), (/nproma, 10/))
        ref_out_array_dp_4d(:,i,:,n) = &
          RESHAPE(MERGE(-1._wp, glb_index_dst + (i - 1) * global_size + &
          &             n * 0.1_wp, -1 == owner_local_dst), (/nproma, 10/))
      END DO
    END DO
    in_array_sp_4d = in_array_dp_4d
    ref_out_array_sp_4d = ref_out_array_dp_4d

    out_array_dp_4d = -1._dp
    CALL check_exchange_mult_4d( &
      out_array1 = out_array_dp_4d(:,:,:,1), &
      in_array1 = in_array_dp_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_dp_4d(:,:,:,1), &
      out_array2 = out_array_dp_4d(:,:,:,2), &
      in_array2 = in_array_dp_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_dp_4d(:,:,:,2), &
      out_array3 = out_array_dp_4d(:,:,:,3), &
      in_array3 = in_array_dp_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_dp_4d(:,:,:,3), &
      out_array4 = out_array_dp_4d(:,:,:,4), &
      in_array4 = in_array_dp_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_dp_4d(:,:,:,4), &
      out_array5 = out_array_dp_4d(:,:,:,5), &
      in_array5 = in_array_dp_4d(:,:,:,5), &
      ref_out_array5 = ref_out_array_dp_4d(:,:,:,5), &
      out_array6 = out_array_dp_4d(:,:,:,6), &
      in_array6 = in_array_dp_4d(:,:,:,6), &
      ref_out_array6 = ref_out_array_dp_4d(:,:,:,6), &
      out_array7 = out_array_dp_4d(:,:,:,7), &
      in_array7 = in_array_dp_4d(:,:,:,7), &
      ref_out_array7 = ref_out_array_dp_4d(:,:,:,7), &
      out_array4d = out_array_dp_4d(:,:,:,8:), &
      in_array4d = in_array_dp_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_dp_4d(:,:,:,8:), &
      comm_pattern = comm_pattern)
    out_array_dp_4d = -1._dp
    out_array_sp_4d = -1._sp
    CALL check_exchange_mixprec_4d_dp( &
      out_array1_dp = out_array_dp_4d(:,:,:,1), &
      in_array1_dp = in_array_dp_4d(:,:,:,1), &
      ref_out_array1_dp = ref_out_array_dp_4d(:,:,:,1), &
      out_array2_dp = out_array_dp_4d(:,:,:,2), &
      in_array2_dp = in_array_dp_4d(:,:,:,2), &
      ref_out_array2_dp = ref_out_array_dp_4d(:,:,:,2), &
      out_array3_dp = out_array_dp_4d(:,:,:,3), &
      in_array3_dp = in_array_dp_4d(:,:,:,3), &
      ref_out_array3_dp = ref_out_array_dp_4d(:,:,:,3), &
      out_array4_dp = out_array_dp_4d(:,:,:,4), &
      in_array4_dp = in_array_dp_4d(:,:,:,4), &
      ref_out_array4_dp = ref_out_array_dp_4d(:,:,:,4), &
      out_array5_dp = out_array_dp_4d(:,:,:,5), &
      in_array5_dp = in_array_dp_4d(:,:,:,5), &
      ref_out_array5_dp = ref_out_array_dp_4d(:,:,:,5), &
      out_array4d_dp = out_array_dp_4d(:,:,:,8:), &
      in_array4d_dp = in_array_dp_4d(:,:,:,8:), &
      ref_out_array4d_dp = ref_out_array_dp_4d(:,:,:,8:), &
      out_array1_sp = out_array_sp_4d(:,:,:,1), &
      in_array1_sp = in_array_sp_4d(:,:,:,1), &
      ref_out_array1_sp = ref_out_array_sp_4d(:,:,:,1), &
      out_array2_sp = out_array_sp_4d(:,:,:,2), &
      in_array2_sp = in_array_sp_4d(:,:,:,2), &
      ref_out_array2_sp = ref_out_array_sp_4d(:,:,:,2), &
      out_array3_sp = out_array_sp_4d(:,:,:,3), &
      in_array3_sp = in_array_sp_4d(:,:,:,3), &
      ref_out_array3_sp = ref_out_array_sp_4d(:,:,:,3), &
      out_array4_sp = out_array_sp_4d(:,:,:,4), &
      in_array4_sp = in_array_sp_4d(:,:,:,4), &
      ref_out_array4_sp = ref_out_array_sp_4d(:,:,:,4), &
      out_array5_sp = out_array_sp_4d(:,:,:,5), &
      in_array5_sp = in_array_sp_4d(:,:,:,5), &
      ref_out_array5_sp = ref_out_array_sp_4d(:,:,:,5), &
      out_array4d_sp = out_array_sp_4d(:,:,:,8:), &
      in_array4d_sp = in_array_sp_4d(:,:,:,8:), &
      ref_out_array4d_sp = ref_out_array_sp_4d(:,:,:,8:), &
      comm_pattern = comm_pattern)
    out_array_dp_4d = in_array_dp_4d
    CALL check_exchange_mult_4d( &
      out_array1 = out_array_dp_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_dp_4d(:,:,:,1), &
      out_array2 = out_array_dp_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_dp_4d(:,:,:,2), &
      out_array3 = out_array_dp_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_dp_4d(:,:,:,3), &
      out_array4 = out_array_dp_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_dp_4d(:,:,:,4), &
      out_array5 = out_array_dp_4d(:,:,:,5), &
      ref_out_array5 = ref_out_array_dp_4d(:,:,:,5), &
      out_array6 = out_array_dp_4d(:,:,:,6), &
      ref_out_array6 = ref_out_array_dp_4d(:,:,:,6), &
      out_array7 = out_array_dp_4d(:,:,:,7), &
      ref_out_array7 = ref_out_array_dp_4d(:,:,:,7), &
      out_array4d = out_array_dp_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_dp_4d(:,:,:,8:), &
      comm_pattern = comm_pattern)
    out_array_dp_4d = in_array_dp_4d
    out_array_sp_4d = in_array_sp_4d
    CALL check_exchange_mixprec_4d_dp( &
      out_array1_dp = out_array_dp_4d(:,:,:,1), &
      ref_out_array1_dp = ref_out_array_dp_4d(:,:,:,1), &
      out_array2_dp = out_array_dp_4d(:,:,:,2), &
      ref_out_array2_dp = ref_out_array_dp_4d(:,:,:,2), &
      out_array3_dp = out_array_dp_4d(:,:,:,3), &
      ref_out_array3_dp = ref_out_array_dp_4d(:,:,:,3), &
      out_array4_dp = out_array_dp_4d(:,:,:,4), &
      ref_out_array4_dp = ref_out_array_dp_4d(:,:,:,4), &
      out_array5_dp = out_array_dp_4d(:,:,:,5), &
      ref_out_array5_dp = ref_out_array_dp_4d(:,:,:,5), &
      out_array4d_dp = out_array_dp_4d(:,:,:,8:), &
      ref_out_array4d_dp = ref_out_array_dp_4d(:,:,:,8:), &
      out_array1_sp = out_array_sp_4d(:,:,:,1), &
      ref_out_array1_sp = ref_out_array_sp_4d(:,:,:,1), &
      out_array2_sp = out_array_sp_4d(:,:,:,2), &
      ref_out_array2_sp = ref_out_array_sp_4d(:,:,:,2), &
      out_array3_sp = out_array_sp_4d(:,:,:,3), &
      ref_out_array3_sp = ref_out_array_sp_4d(:,:,:,3), &
      out_array4_sp = out_array_sp_4d(:,:,:,4), &
      ref_out_array4_sp = ref_out_array_sp_4d(:,:,:,4), &
      out_array5_sp = out_array_sp_4d(:,:,:,5), &
      ref_out_array5_sp = ref_out_array_sp_4d(:,:,:,5), &
      out_array4d_sp = out_array_sp_4d(:,:,:,8:), &
      ref_out_array4d_sp = ref_out_array_sp_4d(:,:,:,8:), &
      comm_pattern = comm_pattern)

    out_array_dp_4d = -1._dp
    CALL check_exchange_mult_4d( &
      out_array1 = out_array_dp_4d(:,:,:,1), &
      in_array1 = in_array_dp_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_dp_4d(:,:,:,1), &
      out_array2 = out_array_dp_4d(:,:,:,2), &
      in_array2 = in_array_dp_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_dp_4d(:,:,:,2), &
      out_array3 = out_array_dp_4d(:,:,:,3), &
      in_array3 = in_array_dp_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_dp_4d(:,:,:,3), &
      out_array4 = out_array_dp_4d(:,:,:,4), &
      in_array4 = in_array_dp_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_dp_4d(:,:,:,4), &
      out_array5 = out_array_dp_4d(:,:,:,5), &
      in_array5 = in_array_dp_4d(:,:,:,5), &
      ref_out_array5 = ref_out_array_dp_4d(:,:,:,5), &
      out_array6 = out_array_dp_4d(:,:,:,6), &
      in_array6 = in_array_dp_4d(:,:,:,6), &
      ref_out_array6 = ref_out_array_dp_4d(:,:,:,6), &
      out_array7 = out_array_dp_4d(:,:,:,7), &
      in_array7 = in_array_dp_4d(:,:,:,7), &
      ref_out_array7 = ref_out_array_dp_4d(:,:,:,7), &
      out_array4d = out_array_dp_4d(:,:,:,8:), &
      in_array4d = in_array_dp_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_dp_4d(:,:,:,8:), &
      nshift = 0, comm_pattern = comm_pattern)
    out_array_dp_4d = -1._dp
    out_array_sp_4d = -1._sp
    CALL check_exchange_mixprec_4d_dp( &
      out_array1_dp = out_array_dp_4d(:,:,:,1), &
      in_array1_dp = in_array_dp_4d(:,:,:,1), &
      ref_out_array1_dp = ref_out_array_dp_4d(:,:,:,1), &
      out_array2_dp = out_array_dp_4d(:,:,:,2), &
      in_array2_dp = in_array_dp_4d(:,:,:,2), &
      ref_out_array2_dp = ref_out_array_dp_4d(:,:,:,2), &
      out_array3_dp = out_array_dp_4d(:,:,:,3), &
      in_array3_dp = in_array_dp_4d(:,:,:,3), &
      ref_out_array3_dp = ref_out_array_dp_4d(:,:,:,3), &
      out_array4_dp = out_array_dp_4d(:,:,:,4), &
      in_array4_dp = in_array_dp_4d(:,:,:,4), &
      ref_out_array4_dp = ref_out_array_dp_4d(:,:,:,4), &
      out_array5_dp = out_array_dp_4d(:,:,:,5), &
      in_array5_dp = in_array_dp_4d(:,:,:,5), &
      ref_out_array5_dp = ref_out_array_dp_4d(:,:,:,5), &
      out_array4d_dp = out_array_dp_4d(:,:,:,8:), &
      in_array4d_dp = in_array_dp_4d(:,:,:,8:), &
      ref_out_array4d_dp = ref_out_array_dp_4d(:,:,:,8:), &
      out_array1_sp = out_array_sp_4d(:,:,:,1), &
      in_array1_sp = in_array_sp_4d(:,:,:,1), &
      ref_out_array1_sp = ref_out_array_sp_4d(:,:,:,1), &
      out_array2_sp = out_array_sp_4d(:,:,:,2), &
      in_array2_sp = in_array_sp_4d(:,:,:,2), &
      ref_out_array2_sp = ref_out_array_sp_4d(:,:,:,2), &
      out_array3_sp = out_array_sp_4d(:,:,:,3), &
      in_array3_sp = in_array_sp_4d(:,:,:,3), &
      ref_out_array3_sp = ref_out_array_sp_4d(:,:,:,3), &
      out_array4_sp = out_array_sp_4d(:,:,:,4), &
      in_array4_sp = in_array_sp_4d(:,:,:,4), &
      ref_out_array4_sp = ref_out_array_sp_4d(:,:,:,4), &
      out_array5_sp = out_array_sp_4d(:,:,:,5), &
      in_array5_sp = in_array_sp_4d(:,:,:,5), &
      ref_out_array5_sp = ref_out_array_sp_4d(:,:,:,5), &
      out_array4d_sp = out_array_sp_4d(:,:,:,8:), &
      in_array4d_sp = in_array_sp_4d(:,:,:,8:), &
      ref_out_array4d_sp = ref_out_array_sp_4d(:,:,:,8:), &
      nshift = 0, comm_pattern = comm_pattern)
    out_array_dp_4d = in_array_dp_4d
    CALL check_exchange_mult_4d( &
      out_array1 = out_array_dp_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_dp_4d(:,:,:,1), &
      out_array2 = out_array_dp_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_dp_4d(:,:,:,2), &
      out_array3 = out_array_dp_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_dp_4d(:,:,:,3), &
      out_array4 = out_array_dp_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_dp_4d(:,:,:,4), &
      out_array5 = out_array_dp_4d(:,:,:,5), &
      ref_out_array5 = ref_out_array_dp_4d(:,:,:,5), &
      out_array6 = out_array_dp_4d(:,:,:,6), &
      ref_out_array6 = ref_out_array_dp_4d(:,:,:,6), &
      out_array7 = out_array_dp_4d(:,:,:,7), &
      ref_out_array7 = ref_out_array_dp_4d(:,:,:,7), &
      out_array4d = out_array_dp_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_dp_4d(:,:,:,8:), &
      nshift = 0, comm_pattern = comm_pattern)
    out_array_dp_4d = in_array_dp_4d
    out_array_sp_4d = in_array_sp_4d
    CALL check_exchange_mixprec_4d_dp( &
      out_array1_dp = out_array_dp_4d(:,:,:,1), &
      ref_out_array1_dp = ref_out_array_dp_4d(:,:,:,1), &
      out_array2_dp = out_array_dp_4d(:,:,:,2), &
      ref_out_array2_dp = ref_out_array_dp_4d(:,:,:,2), &
      out_array3_dp = out_array_dp_4d(:,:,:,3), &
      ref_out_array3_dp = ref_out_array_dp_4d(:,:,:,3), &
      out_array4_dp = out_array_dp_4d(:,:,:,4), &
      ref_out_array4_dp = ref_out_array_dp_4d(:,:,:,4), &
      out_array5_dp = out_array_dp_4d(:,:,:,5), &
      ref_out_array5_dp = ref_out_array_dp_4d(:,:,:,5), &
      out_array4d_dp = out_array_dp_4d(:,:,:,8:), &
      ref_out_array4d_dp = ref_out_array_dp_4d(:,:,:,8:), &
      out_array1_sp = out_array_sp_4d(:,:,:,1), &
      ref_out_array1_sp = ref_out_array_sp_4d(:,:,:,1), &
      out_array2_sp = out_array_sp_4d(:,:,:,2), &
      ref_out_array2_sp = ref_out_array_sp_4d(:,:,:,2), &
      out_array3_sp = out_array_sp_4d(:,:,:,3), &
      ref_out_array3_sp = ref_out_array_sp_4d(:,:,:,3), &
      out_array4_sp = out_array_sp_4d(:,:,:,4), &
      ref_out_array4_sp = ref_out_array_sp_4d(:,:,:,4), &
      out_array5_sp = out_array_sp_4d(:,:,:,5), &
      ref_out_array5_sp = ref_out_array_sp_4d(:,:,:,5), &
      out_array4d_sp = out_array_sp_4d(:,:,:,8:), &
      ref_out_array4d_sp = ref_out_array_sp_4d(:,:,:,8:), &
      nshift = 0, comm_pattern = comm_pattern)

    ref_out_array_dp_4d(:,:5,:,:) = -1._dp
    in_array_dp_4d(:,:5,:,:) = -1._dp
    out_array_dp_4d = -1._dp
    CALL check_exchange_mult_4d( &
      out_array1 = out_array_dp_4d(:,:,:,1), &
      in_array1 = in_array_dp_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_dp_4d(:,:,:,1), &
      out_array2 = out_array_dp_4d(:,:,:,2), &
      in_array2 = in_array_dp_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_dp_4d(:,:,:,2), &
      out_array3 = out_array_dp_4d(:,:,:,3), &
      in_array3 = in_array_dp_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_dp_4d(:,:,:,3), &
      out_array4 = out_array_dp_4d(:,:,:,4), &
      in_array4 = in_array_dp_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_dp_4d(:,:,:,4), &
      out_array5 = out_array_dp_4d(:,:,:,5), &
      in_array5 = in_array_dp_4d(:,:,:,5), &
      ref_out_array5 = ref_out_array_dp_4d(:,:,:,5), &
      out_array6 = out_array_dp_4d(:,:,:,6), &
      in_array6 = in_array_dp_4d(:,:,:,6), &
      ref_out_array6 = ref_out_array_dp_4d(:,:,:,6), &
      out_array7 = out_array_dp_4d(:,:,:,7), &
      in_array7 = in_array_dp_4d(:,:,:,7), &
      ref_out_array7 = ref_out_array_dp_4d(:,:,:,7), &
      out_array4d = out_array_dp_4d(:,:,:,8:), &
      in_array4d = in_array_dp_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_dp_4d(:,:,:,8:), &
      nshift = 5, comm_pattern = comm_pattern)
    ref_out_array_dp_4d(:,:5,:,:) = -1._dp
    ref_out_array_sp_4d(:,:5,:,:) = -1._sp
    in_array_dp_4d(:,:5,:,:) = -1._dp
    in_array_sp_4d(:,:5,:,:) = -1._sp
    out_array_dp_4d = -1._dp
    out_array_sp_4d = -1._sp
    CALL check_exchange_mixprec_4d_dp( &
      out_array1_dp = out_array_dp_4d(:,:,:,1), &
      in_array1_dp = in_array_dp_4d(:,:,:,1), &
      ref_out_array1_dp = ref_out_array_dp_4d(:,:,:,1), &
      out_array2_dp = out_array_dp_4d(:,:,:,2), &
      in_array2_dp = in_array_dp_4d(:,:,:,2), &
      ref_out_array2_dp = ref_out_array_dp_4d(:,:,:,2), &
      out_array3_dp = out_array_dp_4d(:,:,:,3), &
      in_array3_dp = in_array_dp_4d(:,:,:,3), &
      ref_out_array3_dp = ref_out_array_dp_4d(:,:,:,3), &
      out_array4_dp = out_array_dp_4d(:,:,:,4), &
      in_array4_dp = in_array_dp_4d(:,:,:,4), &
      ref_out_array4_dp = ref_out_array_dp_4d(:,:,:,4), &
      out_array5_dp = out_array_dp_4d(:,:,:,5), &
      in_array5_dp = in_array_dp_4d(:,:,:,5), &
      ref_out_array5_dp = ref_out_array_dp_4d(:,:,:,5), &
      out_array4d_dp = out_array_dp_4d(:,:,:,8:), &
      in_array4d_dp = in_array_dp_4d(:,:,:,8:), &
      ref_out_array4d_dp = ref_out_array_dp_4d(:,:,:,8:), &
      out_array1_sp = out_array_sp_4d(:,:,:,1), &
      in_array1_sp = in_array_sp_4d(:,:,:,1), &
      ref_out_array1_sp = ref_out_array_sp_4d(:,:,:,1), &
      out_array2_sp = out_array_sp_4d(:,:,:,2), &
      in_array2_sp = in_array_sp_4d(:,:,:,2), &
      ref_out_array2_sp = ref_out_array_sp_4d(:,:,:,2), &
      out_array3_sp = out_array_sp_4d(:,:,:,3), &
      in_array3_sp = in_array_sp_4d(:,:,:,3), &
      ref_out_array3_sp = ref_out_array_sp_4d(:,:,:,3), &
      out_array4_sp = out_array_sp_4d(:,:,:,4), &
      in_array4_sp = in_array_sp_4d(:,:,:,4), &
      ref_out_array4_sp = ref_out_array_sp_4d(:,:,:,4), &
      out_array5_sp = out_array_sp_4d(:,:,:,5), &
      in_array5_sp = in_array_sp_4d(:,:,:,5), &
      ref_out_array5_sp = ref_out_array_sp_4d(:,:,:,5), &
      out_array4d_sp = out_array_sp_4d(:,:,:,8:), &
      in_array4d_sp = in_array_sp_4d(:,:,:,8:), &
      ref_out_array4d_sp = ref_out_array_sp_4d(:,:,:,8:), &
      nshift = 5, comm_pattern = comm_pattern)
    out_array_dp_4d = in_array_dp_4d
    CALL check_exchange_mult_4d( &
      out_array1 = out_array_dp_4d(:,:,:,1), &
      ref_out_array1 = ref_out_array_dp_4d(:,:,:,1), &
      out_array2 = out_array_dp_4d(:,:,:,2), &
      ref_out_array2 = ref_out_array_dp_4d(:,:,:,2), &
      out_array3 = out_array_dp_4d(:,:,:,3), &
      ref_out_array3 = ref_out_array_dp_4d(:,:,:,3), &
      out_array4 = out_array_dp_4d(:,:,:,4), &
      ref_out_array4 = ref_out_array_dp_4d(:,:,:,4), &
      out_array5 = out_array_dp_4d(:,:,:,5), &
      ref_out_array5 = ref_out_array_dp_4d(:,:,:,5), &
      out_array6 = out_array_dp_4d(:,:,:,6), &
      ref_out_array6 = ref_out_array_dp_4d(:,:,:,6), &
      out_array7 = out_array_dp_4d(:,:,:,7), &
      ref_out_array7 = ref_out_array_dp_4d(:,:,:,7), &
      out_array4d = out_array_dp_4d(:,:,:,8:), &
      ref_out_array4d = ref_out_array_dp_4d(:,:,:,8:), &
      nshift = 5, comm_pattern = comm_pattern)
    out_array_dp_4d = in_array_dp_4d
    out_array_sp_4d = in_array_sp_4d
    CALL check_exchange_mixprec_4d_dp( &
      out_array1_dp = out_array_dp_4d(:,:,:,1), &
      ref_out_array1_dp = ref_out_array_dp_4d(:,:,:,1), &
      out_array2_dp = out_array_dp_4d(:,:,:,2), &
      ref_out_array2_dp = ref_out_array_dp_4d(:,:,:,2), &
      out_array3_dp = out_array_dp_4d(:,:,:,3), &
      ref_out_array3_dp = ref_out_array_dp_4d(:,:,:,3), &
      out_array4_dp = out_array_dp_4d(:,:,:,4), &
      ref_out_array4_dp = ref_out_array_dp_4d(:,:,:,4), &
      out_array5_dp = out_array_dp_4d(:,:,:,5), &
      ref_out_array5_dp = ref_out_array_dp_4d(:,:,:,5), &
      out_array4d_dp = out_array_dp_4d(:,:,:,8:), &
      ref_out_array4d_dp = ref_out_array_dp_4d(:,:,:,8:), &
      out_array1_sp = out_array_sp_4d(:,:,:,1), &
      ref_out_array1_sp = ref_out_array_sp_4d(:,:,:,1), &
      out_array2_sp = out_array_sp_4d(:,:,:,2), &
      ref_out_array2_sp = ref_out_array_sp_4d(:,:,:,2), &
      out_array3_sp = out_array_sp_4d(:,:,:,3), &
      ref_out_array3_sp = ref_out_array_sp_4d(:,:,:,3), &
      out_array4_sp = out_array_sp_4d(:,:,:,4), &
      ref_out_array4_sp = ref_out_array_sp_4d(:,:,:,4), &
      out_array5_sp = out_array_sp_4d(:,:,:,5), &
      ref_out_array5_sp = ref_out_array_sp_4d(:,:,:,5), &
      out_array4d_sp = out_array_sp_4d(:,:,:,8:), &
      ref_out_array4d_sp = ref_out_array_sp_4d(:,:,:,8:), &
      nshift = 5, comm_pattern = comm_pattern)

    DEALLOCATE(in_array_dp_4d, add_array_dp_4d, out_array_dp_4d, &
      &        ref_out_array_dp_4d, in_array_sp_4d, add_array_sp_4d, &
      &        out_array_sp_4d, ref_out_array_sp_4d)

    CALL delete_comm_pattern(comm_pattern)
    CALL deallocate_glb2loc_index_lookup(send_glb2loc_index)

  CONTAINS

    SUBROUTINE check_exchange(in_array_r_2d, in_array_r_3d, &
      &                       in_array_i_2d, in_array_i_3d, &
      &                       in_array_l_2d, in_array_l_3d, &
      &                       out_array_r_2d, out_array_r_3d, &
      &                       out_array_i_2d, out_array_i_3d, &
      &                       out_array_l_2d, out_array_l_3d, &
      &                       add_array_r_2d, add_array_r_3d, &
      &                       add_array_i_2d, add_array_i_3d, &
      &                       ref_out_array_r_2d, ref_out_array_r_3d, &
      &                       ref_out_array_i_2d, ref_out_array_i_3d, &
      &                       ref_out_array_l_2d, ref_out_array_l_3d, &
      &                       comm_pattern)

      REAL(wp), OPTIONAL, INTENT(IN) :: in_array_r_2d(:,:), in_array_r_3d(:,:,:)
      INTEGER, OPTIONAL, INTENT(IN) ::  in_array_i_2d(:,:), in_array_i_3d(:,:,:)
      LOGICAL, OPTIONAL, INTENT(IN) ::  in_array_l_2d(:,:), in_array_l_3d(:,:,:)
      REAL(wp), INTENT(INOUT) :: out_array_r_2d(:,:), out_array_r_3d(:,:,:)
      INTEGER, INTENT(INOUT) ::  out_array_i_2d(:,:), out_array_i_3d(:,:,:)
      LOGICAL, INTENT(INOUT) ::  out_array_l_2d(:,:), out_array_l_3d(:,:,:)
      REAL(wp), OPTIONAL, INTENT(IN) :: add_array_r_2d(:,:), add_array_r_3d(:,:,:)
      INTEGER, OPTIONAL, INTENT(IN) ::  add_array_i_2d(:,:), add_array_i_3d(:,:,:)
      REAL(wp), INTENT(IN) :: ref_out_array_r_2d(:,:), ref_out_array_r_3d(:,:,:)
      INTEGER, INTENT(IN) ::  ref_out_array_i_2d(:,:), ref_out_array_i_3d(:,:,:)
      LOGICAL, INTENT(IN) ::  ref_out_array_l_2d(:,:), ref_out_array_l_3d(:,:,:)
      CLASS(t_comm_pattern), POINTER, INTENT(INOUT) :: comm_pattern









      !$ACC DATA COPYIN(add_array_r_2d) &
      !$ACC   IF(present(add_array_r_2d) .AND. i_am_accel_node)
      !$ACC DATA COPYIN(in_array_r_2d) &
      !$ACC   IF(present(in_array_r_2d) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array_r_2d) IF(i_am_accel_node)
      CALL exchange_data(p_pat=comm_pattern, recv=out_array_r_2d, &
        &                send=in_array_r_2d, add=add_array_r_2d, &
        &                l_recv_exists=.TRUE.)
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      IF (ANY(out_array_r_2d /= ref_out_array_r_2d)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result r_2d"
        CALL finish(method_name, message_text)
      END IF

      !$ACC DATA COPYIN(add_array_r_3d) &
      !$ACC   IF(present(add_array_r_3d) .AND. i_am_accel_node)
      !$ACC DATA COPYIN(in_array_r_3d) &
      !$ACC   IF(present(in_array_r_3d) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array_r_3d) IF(i_am_accel_node)
      CALL exchange_data(p_pat=comm_pattern, recv=out_array_r_3d, &
        &                send=in_array_r_3d, add=add_array_r_3d)
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      IF (ANY(out_array_r_3d /= ref_out_array_r_3d)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result r_3d"
        CALL finish(method_name, message_text)
      END IF

      !$ACC DATA COPYIN(add_array_i_2d) &
      !$ACC   IF(present(add_array_i_2d) .AND. i_am_accel_node)
      !$ACC DATA COPYIN(in_array_i_2d) &
      !$ACC   IF(present(in_array_i_2d) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array_i_2d) IF(i_am_accel_node)
      CALL exchange_data(p_pat=comm_pattern, recv=out_array_i_2d, &
        &                send=in_array_i_2d, add=add_array_i_2d, &
        &                l_recv_exists=.TRUE.)
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      IF (ANY(out_array_i_2d /= ref_out_array_i_2d)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result i_2d"
        CALL finish(method_name, message_text)
      END IF

      !$ACC DATA COPYIN(add_array_i_3d) &
      !$ACC   IF(present(add_array_i_3d) .AND. i_am_accel_node)
      !$ACC DATA COPYIN(in_array_i_3d) &
      !$ACC   IF(present(in_array_i_3d) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array_i_3d) IF(i_am_accel_node)
      CALL exchange_data(p_pat=comm_pattern, recv=out_array_i_3d, &
        &                send=in_array_i_3d, add=add_array_i_3d)
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      IF (ANY(out_array_i_3d /= ref_out_array_i_3d)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result i_3d"
        CALL finish(method_name, message_text)
      END IF

      !$ACC DATA COPYIN(in_array_l_2d) &
      !$ACC   IF(present(in_array_l_2d) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array_l_2d) IF(i_am_accel_node)
      CALL exchange_data(p_pat=comm_pattern, recv=out_array_l_2d, &
        &                send=in_array_l_2d, l_recv_exists=.TRUE.)
      !$ACC END DATA
      !$ACC END DATA
      IF (ANY(out_array_l_2d .NEQV. ref_out_array_l_2d)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result l_2d"
        CALL finish(method_name, message_text)
      END IF

      !$ACC DATA COPYIN(in_array_l_3d) &
      !$ACC   IF(present(in_array_l_3d) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array_l_3d) IF(i_am_accel_node)
      CALL exchange_data(p_pat=comm_pattern, recv=out_array_l_3d, &
        &                send=in_array_l_3d)
      !$ACC END DATA
      !$ACC END DATA
      IF (ANY(out_array_l_3d .NEQV. ref_out_array_l_3d)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result l_3d"
        CALL finish(method_name, message_text)
      END IF





    END SUBROUTINE check_exchange

    SUBROUTINE check_exchange_4de1(in_array, out_array, ref_out_array, &
      &                            comm_pattern)

      REAL(wp), OPTIONAL, INTENT(IN) :: in_array(:,:,:,:)
      REAL(wp), INTENT(INOUT) :: out_array(:,:,:,:)
      REAL(wp), INTENT(IN) :: ref_out_array(:,:,:,:)
      CLASS(t_comm_pattern), POINTER, INTENT(INOUT) :: comm_pattern

      INTEGER :: nfields, ndim2tot

      nfields = SIZE(out_array, 1)
      ndim2tot = nfields * SIZE(out_array, 3)









      !$ACC DATA COPY(out_array) COPYIN(in_array) IF(i_am_accel_node)
      CALL exchange_data_4de1(comm_pattern, nfields, ndim2tot, out_array, &
        &                     in_array)
      !$ACC END DATA

      IF (ANY(out_array /= ref_out_array)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result 4de1"
        CALL finish(method_name, message_text)
      END IF





    END SUBROUTINE check_exchange_4de1

    SUBROUTINE check_exchange_mult_4d( &
      & out_array1, in_array1, ref_out_array1, &
      & out_array2, in_array2, ref_out_array2, &
      & out_array3, in_array3, ref_out_array3, &
      & out_array4, in_array4, ref_out_array4, &
      & out_array5, in_array5, ref_out_array5, &
      & out_array6, in_array6, ref_out_array6, &
      & out_array7, in_array7, ref_out_array7, &
      & out_array4d, in_array4d, ref_out_array4d, &
      & nshift, comm_pattern)

      REAL(dp), INTENT(INOUT), OPTIONAL ::  &
        out_array1(:,:,:), out_array2(:,:,:), out_array3(:,:,:), &
        out_array4(:,:,:), out_array5(:,:,:), out_array6(:,:,:), &
        out_array7(:,:,:)
      REAL(dp), INTENT(IN), OPTIONAL ::  &
        in_array1(:,:,:), in_array2(:,:,:), in_array3(:,:,:), &
        in_array4(:,:,:), in_array5(:,:,:), in_array6(:,:,:), &
        in_array7(:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL ::  &
        ref_out_array1(:,:,:), ref_out_array2(:,:,:), &
        ref_out_array3(:,:,:), ref_out_array4(:,:,:), &
        ref_out_array5(:,:,:), ref_out_array6(:,:,:), &
        ref_out_array7(:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL :: out_array4d(:,:,:,:)
      REAL(dp), INTENT(IN   ), OPTIONAL :: in_array4d(:,:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL :: ref_out_array4d(:,:,:,:)

      INTEGER, OPTIONAL, INTENT(IN) :: nshift
      CLASS(t_comm_pattern), POINTER, INTENT(INOUT) :: comm_pattern

      CALL check_exchange_mult_dp( &
        & out_array1=out_array1, in_array1=in_array1, ref_out_array1=ref_out_array1, &
        & out_array2=out_array2, in_array2=in_array2, ref_out_array2=ref_out_array2, &
        & out_array3=out_array3, in_array3=in_array3, ref_out_array3=ref_out_array3, &
        & out_array4=out_array4, in_array4=in_array4, ref_out_array4=ref_out_array4, &
        & out_array5=out_array5, in_array5=in_array5, ref_out_array5=ref_out_array5, &
        & out_array6=out_array6, in_array6=in_array6, ref_out_array6=ref_out_array6, &
        & out_array7=out_array7, in_array7=in_array7, ref_out_array7=ref_out_array7, &
        & out_array4d=out_array4d, in_array4d=in_array4d, ref_out_array4d=ref_out_array4d, &
        & nshift=nshift, comm_pattern=comm_pattern)

      CALL check_exchange_mult_dp( &
        & out_array1=out_array1, in_array1=in_array1, ref_out_array1=ref_out_array1, &
        & out_array2=out_array2, in_array2=in_array2, ref_out_array2=ref_out_array2, &
        & out_array3=out_array3, in_array3=in_array3, ref_out_array3=ref_out_array3, &
        & out_array4=out_array4, in_array4=in_array4, ref_out_array4=ref_out_array4, &
        & out_array5=out_array5, in_array5=in_array5, ref_out_array5=ref_out_array5, &
        & out_array6=out_array6, in_array6=in_array6, ref_out_array6=ref_out_array6, &
        & out_array7=out_array7, in_array7=in_array7, ref_out_array7=ref_out_array7, &
        ! & out_array4d=out_array4d, in_array4d=in_array4d, ref_out_array4d=ref_out_array4d, &
        & nshift=nshift, comm_pattern=comm_pattern)

    END SUBROUTINE check_exchange_mult_4d

    SUBROUTINE check_exchange_mult_dp( &
      & out_array1, in_array1, ref_out_array1, &
      & out_array2, in_array2, ref_out_array2, &
      & out_array3, in_array3, ref_out_array3, &
      & out_array4, in_array4, ref_out_array4, &
      & out_array5, in_array5, ref_out_array5, &
      & out_array6, in_array6, ref_out_array6, &
      & out_array7, in_array7, ref_out_array7, &
      & out_array4d, in_array4d, ref_out_array4d, &
      & nshift, comm_pattern)

      REAL(dp), INTENT(INOUT), OPTIONAL ::  &
        out_array1(:,:,:), out_array2(:,:,:), out_array3(:,:,:), &
        out_array4(:,:,:), out_array5(:,:,:), out_array6(:,:,:), &
        out_array7(:,:,:)
      REAL(dp), INTENT(IN), OPTIONAL ::  &
        in_array1(:,:,:), in_array2(:,:,:), in_array3(:,:,:), &
        in_array4(:,:,:), in_array5(:,:,:), in_array6(:,:,:), &
        in_array7(:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL ::  &
        ref_out_array1(:,:,:), ref_out_array2(:,:,:), &
        ref_out_array3(:,:,:), ref_out_array4(:,:,:), &
        ref_out_array5(:,:,:), ref_out_array6(:,:,:), &
        ref_out_array7(:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL :: out_array4d(:,:,:,:)
      REAL(dp), INTENT(IN   ), OPTIONAL :: in_array4d(:,:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL :: ref_out_array4d(:,:,:,:)

      INTEGER, OPTIONAL, INTENT(IN) :: nshift
      CLASS(t_comm_pattern), POINTER, INTENT(INOUT) :: comm_pattern

      CALL check_exchange_mult( &
        & out_array1=out_array1, in_array1=in_array1, ref_out_array1=ref_out_array1, &
        & out_array2=out_array2, in_array2=in_array2, ref_out_array2=ref_out_array2, &
        & out_array3=out_array3, in_array3=in_array3, ref_out_array3=ref_out_array3, &
        & out_array4=out_array4, in_array4=in_array4, ref_out_array4=ref_out_array4, &
        & out_array5=out_array5, in_array5=in_array5, ref_out_array5=ref_out_array5, &
        & out_array6=out_array6, in_array6=in_array6, ref_out_array6=ref_out_array6, &
        & out_array7=out_array7, in_array7=in_array7, ref_out_array7=ref_out_array7, &
        & out_array4d=out_array4d, in_array4d=in_array4d, ref_out_array4d=ref_out_array4d, &
        & nshift=nshift, comm_pattern=comm_pattern)
      CALL check_exchange_mult( &
        & out_array1=out_array1, in_array1=in_array1, ref_out_array1=ref_out_array1, &
        & out_array2=out_array2, in_array2=in_array2, ref_out_array2=ref_out_array2, &
        & out_array3=out_array3, in_array3=in_array3, ref_out_array3=ref_out_array3, &
        & out_array4=out_array4, in_array4=in_array4, ref_out_array4=ref_out_array4, &
        & out_array5=out_array5, in_array5=in_array5, ref_out_array5=ref_out_array5, &
        & out_array6=out_array6, in_array6=in_array6, ref_out_array6=ref_out_array6, &
        ! & out_array7=out_array7, in_array7=in_array7, ref_out_array7=ref_out_array7, &
        & out_array4d=out_array4d, in_array4d=in_array4d, ref_out_array4d=ref_out_array4d, &
        & nshift=nshift, comm_pattern=comm_pattern)
      CALL check_exchange_mult( &
        & out_array1=out_array1, in_array1=in_array1, ref_out_array1=ref_out_array1, &
        & out_array2=out_array2, in_array2=in_array2, ref_out_array2=ref_out_array2, &
        & out_array3=out_array3, in_array3=in_array3, ref_out_array3=ref_out_array3, &
        & out_array4=out_array4, in_array4=in_array4, ref_out_array4=ref_out_array4, &
        & out_array5=out_array5, in_array5=in_array5, ref_out_array5=ref_out_array5, &
        ! & out_array6=out_array6, in_array6=in_array6, ref_out_array6=ref_out_array6, &
        ! & out_array7=out_array7, in_array7=in_array7, ref_out_array7=ref_out_array7, &
        & out_array4d=out_array4d, in_array4d=in_array4d, ref_out_array4d=ref_out_array4d, &
        & nshift=nshift, comm_pattern=comm_pattern)
      CALL check_exchange_mult( &
        & out_array1=out_array1, in_array1=in_array1, ref_out_array1=ref_out_array1, &
        & out_array2=out_array2, in_array2=in_array2, ref_out_array2=ref_out_array2, &
        & out_array3=out_array3, in_array3=in_array3, ref_out_array3=ref_out_array3, &
        & out_array4=out_array4, in_array4=in_array4, ref_out_array4=ref_out_array4, &
        ! & out_array5=out_array5, in_array5=in_array5, ref_out_array5=ref_out_array5, &
        ! & out_array6=out_array6, in_array6=in_array6, ref_out_array6=ref_out_array6, &
        ! & out_array7=out_array7, in_array7=in_array7, ref_out_array7=ref_out_array7, &
        & out_array4d=out_array4d, in_array4d=in_array4d, ref_out_array4d=ref_out_array4d, &
        & nshift=nshift, comm_pattern=comm_pattern)
      CALL check_exchange_mult( &
        & out_array1=out_array1, in_array1=in_array1, ref_out_array1=ref_out_array1, &
        & out_array2=out_array2, in_array2=in_array2, ref_out_array2=ref_out_array2, &
        & out_array3=out_array3, in_array3=in_array3, ref_out_array3=ref_out_array3, &
        ! & out_array4=out_array4, in_array4=in_array4, ref_out_array4=ref_out_array4, &
        ! & out_array5=out_array5, in_array5=in_array5, ref_out_array5=ref_out_array5, &
        ! & out_array6=out_array6, in_array6=in_array6, ref_out_array6=ref_out_array6, &
        ! & out_array7=out_array7, in_array7=in_array7, ref_out_array7=ref_out_array7, &
        & out_array4d=out_array4d, in_array4d=in_array4d, ref_out_array4d=ref_out_array4d, &
        & nshift=nshift, comm_pattern=comm_pattern)
      CALL check_exchange_mult( &
        & out_array1=out_array1, in_array1=in_array1, ref_out_array1=ref_out_array1, &
        & out_array2=out_array2, in_array2=in_array2, ref_out_array2=ref_out_array2, &
        ! & out_array3=out_array3, in_array3=in_array3, ref_out_array3=ref_out_array3, &
        ! & out_array4=out_array4, in_array4=in_array4, ref_out_array4=ref_out_array4, &
        ! & out_array5=out_array5, in_array5=in_array5, ref_out_array5=ref_out_array5, &
        ! & out_array6=out_array6, in_array6=in_array6, ref_out_array6=ref_out_array6, &
        ! & out_array7=out_array7, in_array7=in_array7, ref_out_array7=ref_out_array7, &
        & out_array4d=out_array4d, in_array4d=in_array4d, ref_out_array4d=ref_out_array4d, &
        & nshift=nshift, comm_pattern=comm_pattern)
      CALL check_exchange_mult( &
        & out_array1=out_array1, in_array1=in_array1, ref_out_array1=ref_out_array1, &
        ! & out_array2=out_array2, in_array2=in_array2, ref_out_array2=ref_out_array2, &
        ! & out_array3=out_array3, in_array3=in_array3, ref_out_array3=ref_out_array3, &
        ! & out_array4=out_array4, in_array4=in_array4, ref_out_array4=ref_out_array4, &
        ! & out_array5=out_array5, in_array5=in_array5, ref_out_array5=ref_out_array5, &
        ! & out_array6=out_array6, in_array6=in_array6, ref_out_array6=ref_out_array6, &
        ! & out_array7=out_array7, in_array7=in_array7, ref_out_array7=ref_out_array7, &
        & out_array4d=out_array4d, in_array4d=in_array4d, ref_out_array4d=ref_out_array4d, &
        & nshift=nshift, comm_pattern=comm_pattern)
      CALL check_exchange_mult( &
        ! & out_array1=out_array1, in_array1=in_array1, ref_out_array1=ref_out_array1, &
        ! & out_array2=out_array2, in_array2=in_array2, ref_out_array2=ref_out_array2, &
        ! & out_array3=out_array3, in_array3=in_array3, ref_out_array3=ref_out_array3, &
        ! & out_array4=out_array4, in_array4=in_array4, ref_out_array4=ref_out_array4, &
        ! & out_array5=out_array5, in_array5=in_array5, ref_out_array5=ref_out_array5, &
        ! & out_array6=out_array6, in_array6=in_array6, ref_out_array6=ref_out_array6, &
        ! & out_array7=out_array7, in_array7=in_array7, ref_out_array7=ref_out_array7, &
        & out_array4d=out_array4d, in_array4d=in_array4d, ref_out_array4d=ref_out_array4d, &
        & nshift=nshift, comm_pattern=comm_pattern)

    END SUBROUTINE check_exchange_mult_dp

    SUBROUTINE check_exchange_mult( &
      & out_array1, in_array1, ref_out_array1, &
      & out_array2, in_array2, ref_out_array2, &
      & out_array3, in_array3, ref_out_array3, &
      & out_array4, in_array4, ref_out_array4, &
      & out_array5, in_array5, ref_out_array5, &
      & out_array6, in_array6, ref_out_array6, &
      & out_array7, in_array7, ref_out_array7, &
      & out_array4d, in_array4d, ref_out_array4d, &
      & nshift, comm_pattern)

      REAL(dp), INTENT(INOUT), OPTIONAL ::  &
        out_array1(:,:,:), out_array2(:,:,:), out_array3(:,:,:), &
        out_array4(:,:,:), out_array5(:,:,:), out_array6(:,:,:), &
        out_array7(:,:,:)
      REAL(dp), INTENT(IN), OPTIONAL ::  &
        in_array1(:,:,:), in_array2(:,:,:), in_array3(:,:,:), &
        in_array4(:,:,:), in_array5(:,:,:), in_array6(:,:,:), &
        in_array7(:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL ::  &
        ref_out_array1(:,:,:), ref_out_array2(:,:,:), &
        ref_out_array3(:,:,:), ref_out_array4(:,:,:), &
        ref_out_array5(:,:,:), ref_out_array6(:,:,:), &
        ref_out_array7(:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL :: out_array4d(:,:,:,:)
      REAL(dp), INTENT(IN   ), OPTIONAL :: in_array4d(:,:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL :: ref_out_array4d(:,:,:,:)

      INTEGER, OPTIONAL, INTENT(IN) :: nshift
      CLASS(t_comm_pattern), POINTER, INTENT(INOUT) :: comm_pattern


      REAL(dp), ALLOCATABLE ::  &
        tmp_out_array1(:,:,:), tmp_out_array2(:,:,:), tmp_out_array3(:,:,:), &
        tmp_out_array4(:,:,:), tmp_out_array5(:,:,:), tmp_out_array6(:,:,:), &
        tmp_out_array7(:,:,:), tmp_out_array4d(:,:,:,:)

      INTEGER :: nfields, ndim2tot, ndim2, kshift

      kshift = 0
      IF (PRESENT(nshift)) kshift = nshift

      nfields = 0
      ndim2tot = 0
      IF (PRESENT(out_array4d)) THEN
        ALLOCATE(tmp_out_array4d(SIZE(out_array4d,1),SIZE(out_array4d,2),&
                                    SIZE(out_array4d,3),SIZE(out_array4d,4)))
        tmp_out_array4d = out_array4d
        nfields = nfields + SIZE(out_array4d, 4)
        ndim2 = SIZE(out_array4d, 2)
        ndim2tot = ndim2tot + &
          &        MERGE(1, ndim2 - kshift, ndim2 == 1) * SIZE(out_array4d, 4)
      END IF
      IF (PRESENT(out_array1)) THEN
        ALLOCATE(tmp_out_array1(SIZE(out_array1,1),SIZE(out_array1,2),&
                                   SIZE(out_array1,3)))
        tmp_out_array1 = out_array1
        nfields = nfields + 1
        ndim2 = SIZE(out_array1, 2)
        ndim2tot = ndim2tot + MERGE(1, ndim2 - kshift, ndim2 == 1)
        IF (PRESENT(out_array2)) THEN
          ALLOCATE(tmp_out_array2(SIZE(out_array2,1),SIZE(out_array2,2),&
                                     SIZE(out_array2,3)))
          tmp_out_array2 = out_array2
          nfields = nfields + 1
          ndim2 = SIZE(out_array2, 2)
          ndim2tot = ndim2tot + MERGE(1, ndim2 - kshift, ndim2 == 1)
          IF (PRESENT(out_array3)) THEN
            ALLOCATE(tmp_out_array3(SIZE(out_array3,1),SIZE(out_array3,2),&
                                       SIZE(out_array3,3)))
            tmp_out_array3 = out_array3
            nfields = nfields + 1
            ndim2 = SIZE(out_array3, 2)
            ndim2tot = ndim2tot + MERGE(1, ndim2 - kshift, ndim2 == 1)
            IF (PRESENT(out_array4)) THEN
              ALLOCATE(tmp_out_array4(SIZE(out_array4,1),SIZE(out_array4,2),&
                                         SIZE(out_array4,3)))
              tmp_out_array4 = out_array4
              nfields = nfields + 1
              ndim2 = SIZE(out_array4, 2)
              ndim2tot = ndim2tot + MERGE(1, ndim2 - kshift, ndim2 == 1)
              IF (PRESENT(out_array5)) THEN
                ALLOCATE(tmp_out_array5(SIZE(out_array5,1),SIZE(out_array5,2),&
                                           SIZE(out_array5,3)))
                tmp_out_array5 = out_array5
                nfields = nfields + 1
                ndim2 = SIZE(out_array5, 2)
                ndim2tot = ndim2tot + MERGE(1, ndim2 - kshift, ndim2 == 1)
                IF (PRESENT(out_array6)) THEN
                  ALLOCATE(tmp_out_array6(SIZE(out_array6,1),SIZE(out_array6,2),&
                                             SIZE(out_array6,3)))
                  tmp_out_array6 = out_array6
                  nfields = nfields + 1
                  ndim2 = SIZE(out_array6, 2)
                  ndim2tot = ndim2tot + MERGE(1, ndim2 - kshift, ndim2 == 1)
                  IF (PRESENT(out_array7)) THEN
                    ALLOCATE(tmp_out_array7(SIZE(out_array7,1),SIZE(out_array7,2),&
                                               SIZE(out_array7,3)))
                    tmp_out_array7 = out_array7
                    nfields = nfields + 1
                    ndim2 = SIZE(out_array7, 2)
                    ndim2tot = ndim2tot + MERGE(1, ndim2 - kshift, ndim2 == 1)
                  END IF
                END IF
              END IF
            END IF
          END IF
        END IF
      END IF









      !$ACC DATA COPYIN(in_array1) IF(present(in_array1) .AND. i_am_accel_node)
      !$ACC DATA COPYIN(in_array2) IF(present(in_array2) .AND. i_am_accel_node)
      !$ACC DATA COPYIN(in_array3) IF(present(in_array3) .AND. i_am_accel_node)
      !$ACC DATA COPYIN(in_array4) IF(present(in_array4) .AND. i_am_accel_node)
      !$ACC DATA COPYIN(in_array5) IF(present(in_array5) .AND. i_am_accel_node)
      !$ACC DATA COPYIN(in_array6) IF(present(in_array6) .AND. i_am_accel_node)
      !$ACC DATA COPYIN(in_array7) IF(present(in_array7) .AND. i_am_accel_node)
      !$ACC DATA COPYIN(in_array4d) IF(present(in_array4d) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array1) IF(present(out_array1) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array2) IF(present(out_array2) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array3) IF(present(out_array3) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array4) IF(present(out_array4) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array5) IF(present(out_array5) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array6) IF(present(out_array6) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array7) IF(present(out_array7) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array4d) IF(present(out_array4d) .AND. i_am_accel_node)
      IF (nfields > 0) THEN
        CALL exchange_data_mult( &
          p_pat=comm_pattern, nfields=nfields, ndim2tot=ndim2tot, &
          recv1=out_array1, send1=in_array1, &
          recv2=out_array2, send2=in_array2, &
          recv3=out_array3, send3=in_array3, &
          recv4=out_array4, send4=in_array4, &
          recv5=out_array5, send5=in_array5, &
          recv6=out_array6, send6=in_array6, &
          recv7=out_array7, send7=in_array7, &
          recv4d=out_array4d, send4d=in_array4d, &
          nshift=kshift)
      END IF
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA





      IF (PRESENT(out_array1)) THEN
        IF (ANY(out_array1 /= ref_out_array1)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result mult dp_3d_1"
          CALL finish(method_name, message_text)
        END IF
        out_array1 = tmp_out_array1
      END IF
      IF (PRESENT(out_array2)) THEN
        IF (ANY(out_array2 /= ref_out_array2)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result mult dp_3d_2"
          CALL finish(method_name, message_text)
        END IF
        out_array2 = tmp_out_array2
      END IF
      IF (PRESENT(out_array3)) THEN
        IF (ANY(out_array3 /= ref_out_array3)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result mult dp_3d_3"
          CALL finish(method_name, message_text)
        END IF
        out_array3 = tmp_out_array3
      END IF
      IF (PRESENT(out_array4)) THEN
        IF (ANY(out_array4 /= ref_out_array4)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result mult dp_3d_4"
          CALL finish(method_name, message_text)
        END IF
        out_array4 = tmp_out_array4
      END IF
      IF (PRESENT(out_array5)) THEN
        IF (ANY(out_array5 /= ref_out_array5)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result mult dp_3d_5"
          CALL finish(method_name, message_text)
        END IF
        out_array5 = tmp_out_array5
      END IF
      IF (PRESENT(out_array6)) THEN
        IF (ANY(out_array6 /= ref_out_array6)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result mult dp_3d_6"
          CALL finish(method_name, message_text)
        END IF
        out_array6 = tmp_out_array6
      END IF
      IF (PRESENT(out_array7)) THEN
        IF (ANY(out_array7 /= ref_out_array7)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result mult dp_3d_7"
          CALL finish(method_name, message_text)
        END IF
        out_array7 = tmp_out_array7
      END IF
      IF (PRESENT(out_array4d)) THEN
        IF (ANY(out_array4d /= ref_out_array4d)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result mult dp_4d"
          CALL finish(method_name, message_text)
        END IF
        out_array4d = tmp_out_array4d
      END IF
    END SUBROUTINE check_exchange_mult

    SUBROUTINE check_exchange_mixprec_4d_dp( &
      & out_array1_dp, in_array1_dp, ref_out_array1_dp, &
      & out_array2_dp, in_array2_dp, ref_out_array2_dp, &
      & out_array3_dp, in_array3_dp, ref_out_array3_dp, &
      & out_array4_dp, in_array4_dp, ref_out_array4_dp, &
      & out_array5_dp, in_array5_dp, ref_out_array5_dp, &
      & out_array4d_dp, in_array4d_dp, ref_out_array4d_dp, &
      & out_array1_sp, in_array1_sp, ref_out_array1_sp, &
      & out_array2_sp, in_array2_sp, ref_out_array2_sp, &
      & out_array3_sp, in_array3_sp, ref_out_array3_sp, &
      & out_array4_sp, in_array4_sp, ref_out_array4_sp, &
      & out_array5_sp, in_array5_sp, ref_out_array5_sp, &
      & out_array4d_sp, in_array4d_sp, ref_out_array4d_sp, &
      & nshift, comm_pattern)

      REAL(dp), INTENT(INOUT), OPTIONAL ::  &
        out_array1_dp(:,:,:), out_array2_dp(:,:,:), out_array3_dp(:,:,:), &
        out_array4_dp(:,:,:), out_array5_dp(:,:,:)
      REAL(dp), INTENT(IN), OPTIONAL ::  &
        in_array1_dp(:,:,:), in_array2_dp(:,:,:), in_array3_dp(:,:,:), &
        in_array4_dp(:,:,:), in_array5_dp(:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL ::  &
        ref_out_array1_dp(:,:,:), ref_out_array2_dp(:,:,:), &
        ref_out_array3_dp(:,:,:), ref_out_array4_dp(:,:,:), &
        ref_out_array5_dp(:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL :: out_array4d_dp(:,:,:,:)
      REAL(dp), INTENT(IN   ), OPTIONAL :: in_array4d_dp(:,:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL :: ref_out_array4d_dp(:,:,:,:)
      REAL(sp), INTENT(INOUT), OPTIONAL ::  &
        out_array1_sp(:,:,:), out_array2_sp(:,:,:), out_array3_sp(:,:,:), &
        out_array4_sp(:,:,:), out_array5_sp(:,:,:)
      REAL(sp), INTENT(IN), OPTIONAL ::  &
        in_array1_sp(:,:,:), in_array2_sp(:,:,:), in_array3_sp(:,:,:), &
        in_array4_sp(:,:,:), in_array5_sp(:,:,:)
      REAL(sp), INTENT(INOUT), OPTIONAL ::  &
        ref_out_array1_sp(:,:,:), ref_out_array2_sp(:,:,:), &
        ref_out_array3_sp(:,:,:), ref_out_array4_sp(:,:,:), &
        ref_out_array5_sp(:,:,:)
      REAL(sp), INTENT(INOUT), OPTIONAL :: out_array4d_sp(:,:,:,:)
      REAL(sp), INTENT(IN   ), OPTIONAL :: in_array4d_sp(:,:,:,:)
      REAL(sp), INTENT(INOUT), OPTIONAL :: ref_out_array4d_sp(:,:,:,:)

      INTEGER, OPTIONAL, INTENT(IN) :: nshift
      CLASS(t_comm_pattern), POINTER, INTENT(INOUT) :: comm_pattern

      CALL check_exchange_mixprec_4d_sp( &
        & out_array1_dp=out_array1_dp, in_array1_dp=in_array1_dp, ref_out_array1_dp=ref_out_array1_dp, &
        & out_array2_dp=out_array2_dp, in_array2_dp=in_array2_dp, ref_out_array2_dp=ref_out_array2_dp, &
        & out_array3_dp=out_array3_dp, in_array3_dp=in_array3_dp, ref_out_array3_dp=ref_out_array3_dp, &
        & out_array4_dp=out_array4_dp, in_array4_dp=in_array4_dp, ref_out_array4_dp=ref_out_array4_dp, &
        & out_array5_dp=out_array5_dp, in_array5_dp=in_array5_dp, ref_out_array5_dp=ref_out_array5_dp, &
        & out_array4d_dp=out_array4d_dp, in_array4d_dp=in_array4d_dp, ref_out_array4d_dp=ref_out_array4d_dp, &
        & out_array1_sp=out_array1_sp, in_array1_sp=in_array1_sp, ref_out_array1_sp=ref_out_array1_sp, &
        & out_array2_sp=out_array2_sp, in_array2_sp=in_array2_sp, ref_out_array2_sp=ref_out_array2_sp, &
        & out_array3_sp=out_array3_sp, in_array3_sp=in_array3_sp, ref_out_array3_sp=ref_out_array3_sp, &
        & out_array4_sp=out_array4_sp, in_array4_sp=in_array4_sp, ref_out_array4_sp=ref_out_array4_sp, &
        & out_array5_sp=out_array5_sp, in_array5_sp=in_array5_sp, ref_out_array5_sp=ref_out_array5_sp, &
        & out_array4d_sp=out_array4d_sp, in_array4d_sp=in_array4d_sp, ref_out_array4d_sp=ref_out_array4d_sp, &
        & nshift=nshift, comm_pattern=comm_pattern)
      CALL check_exchange_mixprec_4d_sp( &
        & out_array1_dp=out_array1_dp, in_array1_dp=in_array1_dp, ref_out_array1_dp=ref_out_array1_dp, &
        & out_array2_dp=out_array2_dp, in_array2_dp=in_array2_dp, ref_out_array2_dp=ref_out_array2_dp, &
        & out_array3_dp=out_array3_dp, in_array3_dp=in_array3_dp, ref_out_array3_dp=ref_out_array3_dp, &
        & out_array4_dp=out_array4_dp, in_array4_dp=in_array4_dp, ref_out_array4_dp=ref_out_array4_dp, &
        & out_array5_dp=out_array5_dp, in_array5_dp=in_array5_dp, ref_out_array5_dp=ref_out_array5_dp, &
        ! & out_array4d_dp=out_array4d_dp, in_array4d_dp=in_array4d_dp, ref_out_array4d_dp=ref_out_array4d_dp, &
        & out_array1_sp=out_array1_sp, in_array1_sp=in_array1_sp, ref_out_array1_sp=ref_out_array1_sp, &
        & out_array2_sp=out_array2_sp, in_array2_sp=in_array2_sp, ref_out_array2_sp=ref_out_array2_sp, &
        & out_array3_sp=out_array3_sp, in_array3_sp=in_array3_sp, ref_out_array3_sp=ref_out_array3_sp, &
        & out_array4_sp=out_array4_sp, in_array4_sp=in_array4_sp, ref_out_array4_sp=ref_out_array4_sp, &
        & out_array5_sp=out_array5_sp, in_array5_sp=in_array5_sp, ref_out_array5_sp=ref_out_array5_sp, &
        & out_array4d_sp=out_array4d_sp, in_array4d_sp=in_array4d_sp, ref_out_array4d_sp=ref_out_array4d_sp, &
        & nshift=nshift, comm_pattern=comm_pattern)

    END SUBROUTINE check_exchange_mixprec_4d_dp

    SUBROUTINE check_exchange_mixprec_4d_sp( &
      & out_array1_dp, in_array1_dp, ref_out_array1_dp, &
      & out_array2_dp, in_array2_dp, ref_out_array2_dp, &
      & out_array3_dp, in_array3_dp, ref_out_array3_dp, &
      & out_array4_dp, in_array4_dp, ref_out_array4_dp, &
      & out_array5_dp, in_array5_dp, ref_out_array5_dp, &
      & out_array4d_dp, in_array4d_dp, ref_out_array4d_dp, &
      & out_array1_sp, in_array1_sp, ref_out_array1_sp, &
      & out_array2_sp, in_array2_sp, ref_out_array2_sp, &
      & out_array3_sp, in_array3_sp, ref_out_array3_sp, &
      & out_array4_sp, in_array4_sp, ref_out_array4_sp, &
      & out_array5_sp, in_array5_sp, ref_out_array5_sp, &
      & out_array4d_sp, in_array4d_sp, ref_out_array4d_sp, &
      & nshift, comm_pattern)

      REAL(dp), INTENT(INOUT), OPTIONAL ::  &
        out_array1_dp(:,:,:), out_array2_dp(:,:,:), out_array3_dp(:,:,:), &
        out_array4_dp(:,:,:), out_array5_dp(:,:,:)
      REAL(dp), INTENT(IN), OPTIONAL ::  &
        in_array1_dp(:,:,:), in_array2_dp(:,:,:), in_array3_dp(:,:,:), &
        in_array4_dp(:,:,:), in_array5_dp(:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL ::  &
        ref_out_array1_dp(:,:,:), ref_out_array2_dp(:,:,:), &
        ref_out_array3_dp(:,:,:), ref_out_array4_dp(:,:,:), &
        ref_out_array5_dp(:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL :: out_array4d_dp(:,:,:,:)
      REAL(dp), INTENT(IN   ), OPTIONAL :: in_array4d_dp(:,:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL :: ref_out_array4d_dp(:,:,:,:)
      REAL(sp), INTENT(INOUT), OPTIONAL ::  &
        out_array1_sp(:,:,:), out_array2_sp(:,:,:), out_array3_sp(:,:,:), &
        out_array4_sp(:,:,:), out_array5_sp(:,:,:)
      REAL(sp), INTENT(IN), OPTIONAL ::  &
        in_array1_sp(:,:,:), in_array2_sp(:,:,:), in_array3_sp(:,:,:), &
        in_array4_sp(:,:,:), in_array5_sp(:,:,:)
      REAL(sp), INTENT(INOUT), OPTIONAL ::  &
        ref_out_array1_sp(:,:,:), ref_out_array2_sp(:,:,:), &
        ref_out_array3_sp(:,:,:), ref_out_array4_sp(:,:,:), &
        ref_out_array5_sp(:,:,:)
      REAL(sp), INTENT(INOUT), OPTIONAL :: out_array4d_sp(:,:,:,:)
      REAL(sp), INTENT(IN   ), OPTIONAL :: in_array4d_sp(:,:,:,:)
      REAL(sp), INTENT(INOUT), OPTIONAL :: ref_out_array4d_sp(:,:,:,:)

      INTEGER, OPTIONAL, INTENT(IN) :: nshift
      CLASS(t_comm_pattern), POINTER, INTENT(INOUT) :: comm_pattern

      CALL check_exchange_mixprec_dp( &
        & out_array1_dp=out_array1_dp, in_array1_dp=in_array1_dp, ref_out_array1_dp=ref_out_array1_dp, &
        & out_array2_dp=out_array2_dp, in_array2_dp=in_array2_dp, ref_out_array2_dp=ref_out_array2_dp, &
        & out_array3_dp=out_array3_dp, in_array3_dp=in_array3_dp, ref_out_array3_dp=ref_out_array3_dp, &
        & out_array4_dp=out_array4_dp, in_array4_dp=in_array4_dp, ref_out_array4_dp=ref_out_array4_dp, &
        & out_array5_dp=out_array5_dp, in_array5_dp=in_array5_dp, ref_out_array5_dp=ref_out_array5_dp, &
        & out_array4d_dp=out_array4d_dp, in_array4d_dp=in_array4d_dp, ref_out_array4d_dp=ref_out_array4d_dp, &
        & out_array1_sp=out_array1_sp, in_array1_sp=in_array1_sp, ref_out_array1_sp=ref_out_array1_sp, &
        & out_array2_sp=out_array2_sp, in_array2_sp=in_array2_sp, ref_out_array2_sp=ref_out_array2_sp, &
        & out_array3_sp=out_array3_sp, in_array3_sp=in_array3_sp, ref_out_array3_sp=ref_out_array3_sp, &
        & out_array4_sp=out_array4_sp, in_array4_sp=in_array4_sp, ref_out_array4_sp=ref_out_array4_sp, &
        & out_array5_sp=out_array5_sp, in_array5_sp=in_array5_sp, ref_out_array5_sp=ref_out_array5_sp, &
        & out_array4d_sp=out_array4d_sp, in_array4d_sp=in_array4d_sp, ref_out_array4d_sp=ref_out_array4d_sp, &
        & nshift=nshift, comm_pattern=comm_pattern)
      CALL check_exchange_mixprec_dp( &
        & out_array1_dp=out_array1_dp, in_array1_dp=in_array1_dp, ref_out_array1_dp=ref_out_array1_dp, &
        & out_array2_dp=out_array2_dp, in_array2_dp=in_array2_dp, ref_out_array2_dp=ref_out_array2_dp, &
        & out_array3_dp=out_array3_dp, in_array3_dp=in_array3_dp, ref_out_array3_dp=ref_out_array3_dp, &
        & out_array4_dp=out_array4_dp, in_array4_dp=in_array4_dp, ref_out_array4_dp=ref_out_array4_dp, &
        & out_array5_dp=out_array5_dp, in_array5_dp=in_array5_dp, ref_out_array5_dp=ref_out_array5_dp, &
        & out_array4d_dp=out_array4d_dp, in_array4d_dp=in_array4d_dp, ref_out_array4d_dp=ref_out_array4d_dp, &
        & out_array1_sp=out_array1_sp, in_array1_sp=in_array1_sp, ref_out_array1_sp=ref_out_array1_sp, &
        & out_array2_sp=out_array2_sp, in_array2_sp=in_array2_sp, ref_out_array2_sp=ref_out_array2_sp, &
        & out_array3_sp=out_array3_sp, in_array3_sp=in_array3_sp, ref_out_array3_sp=ref_out_array3_sp, &
        & out_array4_sp=out_array4_sp, in_array4_sp=in_array4_sp, ref_out_array4_sp=ref_out_array4_sp, &
        & out_array5_sp=out_array5_sp, in_array5_sp=in_array5_sp, ref_out_array5_sp=ref_out_array5_sp, &
        ! & out_array4d_sp=out_array4d_sp, in_array4d_sp=in_array4d_sp, ref_out_array4d_sp=ref_out_array4d_sp, &
        & nshift=nshift, comm_pattern=comm_pattern)

    END SUBROUTINE check_exchange_mixprec_4d_sp

    SUBROUTINE check_exchange_mixprec_dp( &
      & out_array1_dp, in_array1_dp, ref_out_array1_dp, &
      & out_array2_dp, in_array2_dp, ref_out_array2_dp, &
      & out_array3_dp, in_array3_dp, ref_out_array3_dp, &
      & out_array4_dp, in_array4_dp, ref_out_array4_dp, &
      & out_array5_dp, in_array5_dp, ref_out_array5_dp, &
      & out_array4d_dp, in_array4d_dp, ref_out_array4d_dp, &
      & out_array1_sp, in_array1_sp, ref_out_array1_sp, &
      & out_array2_sp, in_array2_sp, ref_out_array2_sp, &
      & out_array3_sp, in_array3_sp, ref_out_array3_sp, &
      & out_array4_sp, in_array4_sp, ref_out_array4_sp, &
      & out_array5_sp, in_array5_sp, ref_out_array5_sp, &
      & out_array4d_sp, in_array4d_sp, ref_out_array4d_sp, &
      & nshift, comm_pattern)

      REAL(dp), INTENT(INOUT), OPTIONAL ::  &
        out_array1_dp(:,:,:), out_array2_dp(:,:,:), out_array3_dp(:,:,:), &
        out_array4_dp(:,:,:), out_array5_dp(:,:,:)
      REAL(dp), INTENT(IN), OPTIONAL ::  &
        in_array1_dp(:,:,:), in_array2_dp(:,:,:), in_array3_dp(:,:,:), &
        in_array4_dp(:,:,:), in_array5_dp(:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL ::  &
        ref_out_array1_dp(:,:,:), ref_out_array2_dp(:,:,:), &
        ref_out_array3_dp(:,:,:), ref_out_array4_dp(:,:,:), &
        ref_out_array5_dp(:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL :: out_array4d_dp(:,:,:,:)
      REAL(dp), INTENT(IN   ), OPTIONAL :: in_array4d_dp(:,:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL :: ref_out_array4d_dp(:,:,:,:)
      REAL(sp), INTENT(INOUT), OPTIONAL ::  &
        out_array1_sp(:,:,:), out_array2_sp(:,:,:), out_array3_sp(:,:,:), &
        out_array4_sp(:,:,:), out_array5_sp(:,:,:)
      REAL(sp), INTENT(IN), OPTIONAL ::  &
        in_array1_sp(:,:,:), in_array2_sp(:,:,:), in_array3_sp(:,:,:), &
        in_array4_sp(:,:,:), in_array5_sp(:,:,:)
      REAL(sp), INTENT(INOUT), OPTIONAL ::  &
        ref_out_array1_sp(:,:,:), ref_out_array2_sp(:,:,:), &
        ref_out_array3_sp(:,:,:), ref_out_array4_sp(:,:,:), &
        ref_out_array5_sp(:,:,:)
      REAL(sp), INTENT(INOUT), OPTIONAL :: out_array4d_sp(:,:,:,:)
      REAL(sp), INTENT(IN   ), OPTIONAL :: in_array4d_sp(:,:,:,:)
      REAL(sp), INTENT(INOUT), OPTIONAL :: ref_out_array4d_sp(:,:,:,:)

      INTEGER, OPTIONAL, INTENT(IN) :: nshift
      CLASS(t_comm_pattern), POINTER, INTENT(INOUT) :: comm_pattern

      CALL check_exchange_mixprec_sp( &
        & out_array1_dp=out_array1_dp, in_array1_dp=in_array1_dp, ref_out_array1_dp=ref_out_array1_dp, &
        & out_array2_dp=out_array2_dp, in_array2_dp=in_array2_dp, ref_out_array2_dp=ref_out_array2_dp, &
        & out_array3_dp=out_array3_dp, in_array3_dp=in_array3_dp, ref_out_array3_dp=ref_out_array3_dp, &
        & out_array4_dp=out_array4_dp, in_array4_dp=in_array4_dp, ref_out_array4_dp=ref_out_array4_dp, &
        & out_array5_dp=out_array5_dp, in_array5_dp=in_array5_dp, ref_out_array5_dp=ref_out_array5_dp, &
        & out_array4d_dp=out_array4d_dp, in_array4d_dp=in_array4d_dp, ref_out_array4d_dp=ref_out_array4d_dp, &
        & out_array1_sp=out_array1_sp, in_array1_sp=in_array1_sp, ref_out_array1_sp=ref_out_array1_sp, &
        & out_array2_sp=out_array2_sp, in_array2_sp=in_array2_sp, ref_out_array2_sp=ref_out_array2_sp, &
        & out_array3_sp=out_array3_sp, in_array3_sp=in_array3_sp, ref_out_array3_sp=ref_out_array3_sp, &
        & out_array4_sp=out_array4_sp, in_array4_sp=in_array4_sp, ref_out_array4_sp=ref_out_array4_sp, &
        & out_array5_sp=out_array5_sp, in_array5_sp=in_array5_sp, ref_out_array5_sp=ref_out_array5_sp, &
        & out_array4d_sp=out_array4d_sp, in_array4d_sp=in_array4d_sp, ref_out_array4d_sp=ref_out_array4d_sp, &
        & nshift=nshift, comm_pattern=comm_pattern)
      CALL check_exchange_mixprec_sp( &
        & out_array1_dp=out_array1_dp, in_array1_dp=in_array1_dp, ref_out_array1_dp=ref_out_array1_dp, &
        & out_array2_dp=out_array2_dp, in_array2_dp=in_array2_dp, ref_out_array2_dp=ref_out_array2_dp, &
        & out_array3_dp=out_array3_dp, in_array3_dp=in_array3_dp, ref_out_array3_dp=ref_out_array3_dp, &
        & out_array4_dp=out_array4_dp, in_array4_dp=in_array4_dp, ref_out_array4_dp=ref_out_array4_dp, &
        ! & out_array5_dp=out_array5_dp, in_array5_dp=in_array5_dp, ref_out_array5_dp=ref_out_array5_dp, &
        & out_array4d_dp=out_array4d_dp, in_array4d_dp=in_array4d_dp, ref_out_array4d_dp=ref_out_array4d_dp, &
        & out_array1_sp=out_array1_sp, in_array1_sp=in_array1_sp, ref_out_array1_sp=ref_out_array1_sp, &
        & out_array2_sp=out_array2_sp, in_array2_sp=in_array2_sp, ref_out_array2_sp=ref_out_array2_sp, &
        & out_array3_sp=out_array3_sp, in_array3_sp=in_array3_sp, ref_out_array3_sp=ref_out_array3_sp, &
        & out_array4_sp=out_array4_sp, in_array4_sp=in_array4_sp, ref_out_array4_sp=ref_out_array4_sp, &
        & out_array5_sp=out_array5_sp, in_array5_sp=in_array5_sp, ref_out_array5_sp=ref_out_array5_sp, &
        & out_array4d_sp=out_array4d_sp, in_array4d_sp=in_array4d_sp, ref_out_array4d_sp=ref_out_array4d_sp, &
        & nshift=nshift, comm_pattern=comm_pattern)
      CALL check_exchange_mixprec_sp( &
        & out_array1_dp=out_array1_dp, in_array1_dp=in_array1_dp, ref_out_array1_dp=ref_out_array1_dp, &
        & out_array2_dp=out_array2_dp, in_array2_dp=in_array2_dp, ref_out_array2_dp=ref_out_array2_dp, &
        & out_array3_dp=out_array3_dp, in_array3_dp=in_array3_dp, ref_out_array3_dp=ref_out_array3_dp, &
        ! & out_array4_dp=out_array4_dp, in_array4_dp=in_array4_dp, ref_out_array4_dp=ref_out_array4_dp, &
        ! & out_array5_dp=out_array5_dp, in_array5_dp=in_array5_dp, ref_out_array5_dp=ref_out_array5_dp, &
        & out_array4d_dp=out_array4d_dp, in_array4d_dp=in_array4d_dp, ref_out_array4d_dp=ref_out_array4d_dp, &
        & out_array1_sp=out_array1_sp, in_array1_sp=in_array1_sp, ref_out_array1_sp=ref_out_array1_sp, &
        & out_array2_sp=out_array2_sp, in_array2_sp=in_array2_sp, ref_out_array2_sp=ref_out_array2_sp, &
        & out_array3_sp=out_array3_sp, in_array3_sp=in_array3_sp, ref_out_array3_sp=ref_out_array3_sp, &
        & out_array4_sp=out_array4_sp, in_array4_sp=in_array4_sp, ref_out_array4_sp=ref_out_array4_sp, &
        & out_array5_sp=out_array5_sp, in_array5_sp=in_array5_sp, ref_out_array5_sp=ref_out_array5_sp, &
        & out_array4d_sp=out_array4d_sp, in_array4d_sp=in_array4d_sp, ref_out_array4d_sp=ref_out_array4d_sp, &
        & nshift=nshift, comm_pattern=comm_pattern)
      CALL check_exchange_mixprec_sp( &
        & out_array1_dp=out_array1_dp, in_array1_dp=in_array1_dp, ref_out_array1_dp=ref_out_array1_dp, &
        & out_array2_dp=out_array2_dp, in_array2_dp=in_array2_dp, ref_out_array2_dp=ref_out_array2_dp, &
        ! & out_array3_dp=out_array3_dp, in_array3_dp=in_array3_dp, ref_out_array3_dp=ref_out_array3_dp, &
        ! & out_array4_dp=out_array4_dp, in_array4_dp=in_array4_dp, ref_out_array4_dp=ref_out_array4_dp, &
        ! & out_array5_dp=out_array5_dp, in_array5_dp=in_array5_dp, ref_out_array5_dp=ref_out_array5_dp, &
        & out_array4d_dp=out_array4d_dp, in_array4d_dp=in_array4d_dp, ref_out_array4d_dp=ref_out_array4d_dp, &
        & out_array1_sp=out_array1_sp, in_array1_sp=in_array1_sp, ref_out_array1_sp=ref_out_array1_sp, &
        & out_array2_sp=out_array2_sp, in_array2_sp=in_array2_sp, ref_out_array2_sp=ref_out_array2_sp, &
        & out_array3_sp=out_array3_sp, in_array3_sp=in_array3_sp, ref_out_array3_sp=ref_out_array3_sp, &
        & out_array4_sp=out_array4_sp, in_array4_sp=in_array4_sp, ref_out_array4_sp=ref_out_array4_sp, &
        & out_array5_sp=out_array5_sp, in_array5_sp=in_array5_sp, ref_out_array5_sp=ref_out_array5_sp, &
        & out_array4d_sp=out_array4d_sp, in_array4d_sp=in_array4d_sp, ref_out_array4d_sp=ref_out_array4d_sp, &
        & nshift=nshift, comm_pattern=comm_pattern)
      CALL check_exchange_mixprec_sp( &
        & out_array1_dp=out_array1_dp, in_array1_dp=in_array1_dp, ref_out_array1_dp=ref_out_array1_dp, &
        ! & out_array2_dp=out_array2_dp, in_array2_dp=in_array2_dp, ref_out_array2_dp=ref_out_array2_dp, &
        ! & out_array3_dp=out_array3_dp, in_array3_dp=in_array3_dp, ref_out_array3_dp=ref_out_array3_dp, &
        ! & out_array4_dp=out_array4_dp, in_array4_dp=in_array4_dp, ref_out_array4_dp=ref_out_array4_dp, &
        ! & out_array5_dp=out_array5_dp, in_array5_dp=in_array5_dp, ref_out_array5_dp=ref_out_array5_dp, &
        & out_array4d_dp=out_array4d_dp, in_array4d_dp=in_array4d_dp, ref_out_array4d_dp=ref_out_array4d_dp, &
        & out_array1_sp=out_array1_sp, in_array1_sp=in_array1_sp, ref_out_array1_sp=ref_out_array1_sp, &
        & out_array2_sp=out_array2_sp, in_array2_sp=in_array2_sp, ref_out_array2_sp=ref_out_array2_sp, &
        & out_array3_sp=out_array3_sp, in_array3_sp=in_array3_sp, ref_out_array3_sp=ref_out_array3_sp, &
        & out_array4_sp=out_array4_sp, in_array4_sp=in_array4_sp, ref_out_array4_sp=ref_out_array4_sp, &
        & out_array5_sp=out_array5_sp, in_array5_sp=in_array5_sp, ref_out_array5_sp=ref_out_array5_sp, &
        & out_array4d_sp=out_array4d_sp, in_array4d_sp=in_array4d_sp, ref_out_array4d_sp=ref_out_array4d_sp, &
        & nshift=nshift, comm_pattern=comm_pattern)
      CALL check_exchange_mixprec_sp( &
        ! & out_array1_dp=out_array1_dp, in_array1_dp=in_array1_dp, ref_out_array1_dp=ref_out_array1_dp, &
        ! & out_array2_dp=out_array2_dp, in_array2_dp=in_array2_dp, ref_out_array2_dp=ref_out_array2_dp, &
        ! & out_array3_dp=out_array3_dp, in_array3_dp=in_array3_dp, ref_out_array3_dp=ref_out_array3_dp, &
        ! & out_array4_dp=out_array4_dp, in_array4_dp=in_array4_dp, ref_out_array4_dp=ref_out_array4_dp, &
        ! & out_array5_dp=out_array5_dp, in_array5_dp=in_array5_dp, ref_out_array5_dp=ref_out_array5_dp, &
        & out_array4d_dp=out_array4d_dp, in_array4d_dp=in_array4d_dp, ref_out_array4d_dp=ref_out_array4d_dp, &
        & out_array1_sp=out_array1_sp, in_array1_sp=in_array1_sp, ref_out_array1_sp=ref_out_array1_sp, &
        & out_array2_sp=out_array2_sp, in_array2_sp=in_array2_sp, ref_out_array2_sp=ref_out_array2_sp, &
        & out_array3_sp=out_array3_sp, in_array3_sp=in_array3_sp, ref_out_array3_sp=ref_out_array3_sp, &
        & out_array4_sp=out_array4_sp, in_array4_sp=in_array4_sp, ref_out_array4_sp=ref_out_array4_sp, &
        & out_array5_sp=out_array5_sp, in_array5_sp=in_array5_sp, ref_out_array5_sp=ref_out_array5_sp, &
        & out_array4d_sp=out_array4d_sp, in_array4d_sp=in_array4d_sp, ref_out_array4d_sp=ref_out_array4d_sp, &
        & nshift=nshift, comm_pattern=comm_pattern)

    END SUBROUTINE check_exchange_mixprec_dp

    SUBROUTINE check_exchange_mixprec_sp( &
      & out_array1_dp, in_array1_dp, ref_out_array1_dp, &
      & out_array2_dp, in_array2_dp, ref_out_array2_dp, &
      & out_array3_dp, in_array3_dp, ref_out_array3_dp, &
      & out_array4_dp, in_array4_dp, ref_out_array4_dp, &
      & out_array5_dp, in_array5_dp, ref_out_array5_dp, &
      & out_array4d_dp, in_array4d_dp, ref_out_array4d_dp, &
      & out_array1_sp, in_array1_sp, ref_out_array1_sp, &
      & out_array2_sp, in_array2_sp, ref_out_array2_sp, &
      & out_array3_sp, in_array3_sp, ref_out_array3_sp, &
      & out_array4_sp, in_array4_sp, ref_out_array4_sp, &
      & out_array5_sp, in_array5_sp, ref_out_array5_sp, &
      & out_array4d_sp, in_array4d_sp, ref_out_array4d_sp, &
      & nshift, comm_pattern)

      REAL(dp), INTENT(INOUT), OPTIONAL ::  &
        out_array1_dp(:,:,:), out_array2_dp(:,:,:), out_array3_dp(:,:,:), &
        out_array4_dp(:,:,:), out_array5_dp(:,:,:)
      REAL(dp), INTENT(IN), OPTIONAL ::  &
        in_array1_dp(:,:,:), in_array2_dp(:,:,:), in_array3_dp(:,:,:), &
        in_array4_dp(:,:,:), in_array5_dp(:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL ::  &
        ref_out_array1_dp(:,:,:), ref_out_array2_dp(:,:,:), &
        ref_out_array3_dp(:,:,:), ref_out_array4_dp(:,:,:), &
        ref_out_array5_dp(:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL :: out_array4d_dp(:,:,:,:)
      REAL(dp), INTENT(IN   ), OPTIONAL :: in_array4d_dp(:,:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL :: ref_out_array4d_dp(:,:,:,:)
      REAL(sp), INTENT(INOUT), OPTIONAL ::  &
        out_array1_sp(:,:,:), out_array2_sp(:,:,:), out_array3_sp(:,:,:), &
        out_array4_sp(:,:,:), out_array5_sp(:,:,:)
      REAL(sp), INTENT(IN), OPTIONAL ::  &
        in_array1_sp(:,:,:), in_array2_sp(:,:,:), in_array3_sp(:,:,:), &
        in_array4_sp(:,:,:), in_array5_sp(:,:,:)
      REAL(sp), INTENT(INOUT), OPTIONAL ::  &
        ref_out_array1_sp(:,:,:), ref_out_array2_sp(:,:,:), &
        ref_out_array3_sp(:,:,:), ref_out_array4_sp(:,:,:), &
        ref_out_array5_sp(:,:,:)
      REAL(sp), INTENT(INOUT), OPTIONAL :: out_array4d_sp(:,:,:,:)
      REAL(sp), INTENT(IN   ), OPTIONAL :: in_array4d_sp(:,:,:,:)
      REAL(sp), INTENT(INOUT), OPTIONAL :: ref_out_array4d_sp(:,:,:,:)

      INTEGER, OPTIONAL, INTENT(IN) :: nshift
      CLASS(t_comm_pattern), POINTER, INTENT(INOUT) :: comm_pattern

      CALL check_exchange_mixprec( &
        & out_array1_dp=out_array1_dp, in_array1_dp=in_array1_dp, ref_out_array1_dp=ref_out_array1_dp, &
        & out_array2_dp=out_array2_dp, in_array2_dp=in_array2_dp, ref_out_array2_dp=ref_out_array2_dp, &
        & out_array3_dp=out_array3_dp, in_array3_dp=in_array3_dp, ref_out_array3_dp=ref_out_array3_dp, &
        & out_array4_dp=out_array4_dp, in_array4_dp=in_array4_dp, ref_out_array4_dp=ref_out_array4_dp, &
        & out_array5_dp=out_array5_dp, in_array5_dp=in_array5_dp, ref_out_array5_dp=ref_out_array5_dp, &
        & out_array4d_dp=out_array4d_dp, in_array4d_dp=in_array4d_dp, ref_out_array4d_dp=ref_out_array4d_dp, &
        & out_array1_sp=out_array1_sp, in_array1_sp=in_array1_sp, ref_out_array1_sp=ref_out_array1_sp, &
        & out_array2_sp=out_array2_sp, in_array2_sp=in_array2_sp, ref_out_array2_sp=ref_out_array2_sp, &
        & out_array3_sp=out_array3_sp, in_array3_sp=in_array3_sp, ref_out_array3_sp=ref_out_array3_sp, &
        & out_array4_sp=out_array4_sp, in_array4_sp=in_array4_sp, ref_out_array4_sp=ref_out_array4_sp, &
        & out_array5_sp=out_array5_sp, in_array5_sp=in_array5_sp, ref_out_array5_sp=ref_out_array5_sp, &
        & out_array4d_sp=out_array4d_sp, in_array4d_sp=in_array4d_sp, ref_out_array4d_sp=ref_out_array4d_sp, &
        & nshift=nshift, comm_pattern=comm_pattern)
      CALL check_exchange_mixprec( &
        & out_array1_dp=out_array1_dp, in_array1_dp=in_array1_dp, ref_out_array1_dp=ref_out_array1_dp, &
        & out_array2_dp=out_array2_dp, in_array2_dp=in_array2_dp, ref_out_array2_dp=ref_out_array2_dp, &
        & out_array3_dp=out_array3_dp, in_array3_dp=in_array3_dp, ref_out_array3_dp=ref_out_array3_dp, &
        & out_array4_dp=out_array4_dp, in_array4_dp=in_array4_dp, ref_out_array4_dp=ref_out_array4_dp, &
        & out_array5_dp=out_array5_dp, in_array5_dp=in_array5_dp, ref_out_array5_dp=ref_out_array5_dp, &
        & out_array4d_dp=out_array4d_dp, in_array4d_dp=in_array4d_dp, ref_out_array4d_dp=ref_out_array4d_dp, &
        & out_array1_sp=out_array1_sp, in_array1_sp=in_array1_sp, ref_out_array1_sp=ref_out_array1_sp, &
        & out_array2_sp=out_array2_sp, in_array2_sp=in_array2_sp, ref_out_array2_sp=ref_out_array2_sp, &
        & out_array3_sp=out_array3_sp, in_array3_sp=in_array3_sp, ref_out_array3_sp=ref_out_array3_sp, &
        & out_array4_sp=out_array4_sp, in_array4_sp=in_array4_sp, ref_out_array4_sp=ref_out_array4_sp, &
        ! & out_array5_sp=out_array5_sp, in_array5_sp=in_array5_sp, ref_out_array5_sp=ref_out_array5_sp, &
        & out_array4d_sp=out_array4d_sp, in_array4d_sp=in_array4d_sp, ref_out_array4d_sp=ref_out_array4d_sp, &
        & nshift=nshift, comm_pattern=comm_pattern)
      CALL check_exchange_mixprec( &
        & out_array1_dp=out_array1_dp, in_array1_dp=in_array1_dp, ref_out_array1_dp=ref_out_array1_dp, &
        & out_array2_dp=out_array2_dp, in_array2_dp=in_array2_dp, ref_out_array2_dp=ref_out_array2_dp, &
        & out_array3_dp=out_array3_dp, in_array3_dp=in_array3_dp, ref_out_array3_dp=ref_out_array3_dp, &
        & out_array4_dp=out_array4_dp, in_array4_dp=in_array4_dp, ref_out_array4_dp=ref_out_array4_dp, &
        & out_array5_dp=out_array5_dp, in_array5_dp=in_array5_dp, ref_out_array5_dp=ref_out_array5_dp, &
        & out_array4d_dp=out_array4d_dp, in_array4d_dp=in_array4d_dp, ref_out_array4d_dp=ref_out_array4d_dp, &
        & out_array1_sp=out_array1_sp, in_array1_sp=in_array1_sp, ref_out_array1_sp=ref_out_array1_sp, &
        & out_array2_sp=out_array2_sp, in_array2_sp=in_array2_sp, ref_out_array2_sp=ref_out_array2_sp, &
        & out_array3_sp=out_array3_sp, in_array3_sp=in_array3_sp, ref_out_array3_sp=ref_out_array3_sp, &
        ! & out_array4_sp=out_array4_sp, in_array4_sp=in_array4_sp, ref_out_array4_sp=ref_out_array4_sp, &
        ! & out_array5_sp=out_array5_sp, in_array5_sp=in_array5_sp, ref_out_array5_sp=ref_out_array5_sp, &
        & out_array4d_sp=out_array4d_sp, in_array4d_sp=in_array4d_sp, ref_out_array4d_sp=ref_out_array4d_sp, &
        & nshift=nshift, comm_pattern=comm_pattern)
      CALL check_exchange_mixprec( &
        & out_array1_dp=out_array1_dp, in_array1_dp=in_array1_dp, ref_out_array1_dp=ref_out_array1_dp, &
        & out_array2_dp=out_array2_dp, in_array2_dp=in_array2_dp, ref_out_array2_dp=ref_out_array2_dp, &
        & out_array3_dp=out_array3_dp, in_array3_dp=in_array3_dp, ref_out_array3_dp=ref_out_array3_dp, &
        & out_array4_dp=out_array4_dp, in_array4_dp=in_array4_dp, ref_out_array4_dp=ref_out_array4_dp, &
        & out_array5_dp=out_array5_dp, in_array5_dp=in_array5_dp, ref_out_array5_dp=ref_out_array5_dp, &
        & out_array4d_dp=out_array4d_dp, in_array4d_dp=in_array4d_dp, ref_out_array4d_dp=ref_out_array4d_dp, &
        & out_array1_sp=out_array1_sp, in_array1_sp=in_array1_sp, ref_out_array1_sp=ref_out_array1_sp, &
        & out_array2_sp=out_array2_sp, in_array2_sp=in_array2_sp, ref_out_array2_sp=ref_out_array2_sp, &
        ! & out_array3_sp=out_array3_sp, in_array3_sp=in_array3_sp, ref_out_array3_sp=ref_out_array3_sp, &
        ! & out_array4_sp=out_array4_sp, in_array4_sp=in_array4_sp, ref_out_array4_sp=ref_out_array4_sp, &
        ! & out_array5_sp=out_array5_sp, in_array5_sp=in_array5_sp, ref_out_array5_sp=ref_out_array5_sp, &
        & out_array4d_sp=out_array4d_sp, in_array4d_sp=in_array4d_sp, ref_out_array4d_sp=ref_out_array4d_sp, &
        & nshift=nshift, comm_pattern=comm_pattern)
      CALL check_exchange_mixprec( &
        & out_array1_dp=out_array1_dp, in_array1_dp=in_array1_dp, ref_out_array1_dp=ref_out_array1_dp, &
        & out_array2_dp=out_array2_dp, in_array2_dp=in_array2_dp, ref_out_array2_dp=ref_out_array2_dp, &
        & out_array3_dp=out_array3_dp, in_array3_dp=in_array3_dp, ref_out_array3_dp=ref_out_array3_dp, &
        & out_array4_dp=out_array4_dp, in_array4_dp=in_array4_dp, ref_out_array4_dp=ref_out_array4_dp, &
        & out_array5_dp=out_array5_dp, in_array5_dp=in_array5_dp, ref_out_array5_dp=ref_out_array5_dp, &
        & out_array4d_dp=out_array4d_dp, in_array4d_dp=in_array4d_dp, ref_out_array4d_dp=ref_out_array4d_dp, &
        & out_array1_sp=out_array1_sp, in_array1_sp=in_array1_sp, ref_out_array1_sp=ref_out_array1_sp, &
        ! & out_array2_sp=out_array2_sp, in_array2_sp=in_array2_sp, ref_out_array2_sp=ref_out_array2_sp, &
        ! & out_array3_sp=out_array3_sp, in_array3_sp=in_array3_sp, ref_out_array3_sp=ref_out_array3_sp, &
        ! & out_array4_sp=out_array4_sp, in_array4_sp=in_array4_sp, ref_out_array4_sp=ref_out_array4_sp, &
        ! & out_array5_sp=out_array5_sp, in_array5_sp=in_array5_sp, ref_out_array5_sp=ref_out_array5_sp, &
        & out_array4d_sp=out_array4d_sp, in_array4d_sp=in_array4d_sp, ref_out_array4d_sp=ref_out_array4d_sp, &
        & nshift=nshift, comm_pattern=comm_pattern)
      CALL check_exchange_mixprec( &
        & out_array1_dp=out_array1_dp, in_array1_dp=in_array1_dp, ref_out_array1_dp=ref_out_array1_dp, &
        & out_array2_dp=out_array2_dp, in_array2_dp=in_array2_dp, ref_out_array2_dp=ref_out_array2_dp, &
        & out_array3_dp=out_array3_dp, in_array3_dp=in_array3_dp, ref_out_array3_dp=ref_out_array3_dp, &
        & out_array4_dp=out_array4_dp, in_array4_dp=in_array4_dp, ref_out_array4_dp=ref_out_array4_dp, &
        & out_array5_dp=out_array5_dp, in_array5_dp=in_array5_dp, ref_out_array5_dp=ref_out_array5_dp, &
        & out_array4d_dp=out_array4d_dp, in_array4d_dp=in_array4d_dp, ref_out_array4d_dp=ref_out_array4d_dp, &
        ! & out_array1_sp=out_array1_sp, in_array1_sp=in_array1_sp, ref_out_array1_sp=ref_out_array1_sp, &
        ! & out_array2_sp=out_array2_sp, in_array2_sp=in_array2_sp, ref_out_array2_sp=ref_out_array2_sp, &
        ! & out_array3_sp=out_array3_sp, in_array3_sp=in_array3_sp, ref_out_array3_sp=ref_out_array3_sp, &
        ! & out_array4_sp=out_array4_sp, in_array4_sp=in_array4_sp, ref_out_array4_sp=ref_out_array4_sp, &
        ! & out_array5_sp=out_array5_sp, in_array5_sp=in_array5_sp, ref_out_array5_sp=ref_out_array5_sp, &
        & out_array4d_sp=out_array4d_sp, in_array4d_sp=in_array4d_sp, ref_out_array4d_sp=ref_out_array4d_sp, &
        & nshift=nshift, comm_pattern=comm_pattern)

    END SUBROUTINE check_exchange_mixprec_sp

    SUBROUTINE check_exchange_mixprec( &
      & out_array1_dp, in_array1_dp, ref_out_array1_dp, &
      & out_array2_dp, in_array2_dp, ref_out_array2_dp, &
      & out_array3_dp, in_array3_dp, ref_out_array3_dp, &
      & out_array4_dp, in_array4_dp, ref_out_array4_dp, &
      & out_array5_dp, in_array5_dp, ref_out_array5_dp, &
      & out_array4d_dp, in_array4d_dp, ref_out_array4d_dp, &
      & out_array1_sp, in_array1_sp, ref_out_array1_sp, &
      & out_array2_sp, in_array2_sp, ref_out_array2_sp, &
      & out_array3_sp, in_array3_sp, ref_out_array3_sp, &
      & out_array4_sp, in_array4_sp, ref_out_array4_sp, &
      & out_array5_sp, in_array5_sp, ref_out_array5_sp, &
      & out_array4d_sp, in_array4d_sp, ref_out_array4d_sp, &
      & nshift, comm_pattern)

      REAL(dp), INTENT(INOUT), OPTIONAL ::  &
        out_array1_dp(:,:,:), out_array2_dp(:,:,:), out_array3_dp(:,:,:), &
        out_array4_dp(:,:,:), out_array5_dp(:,:,:)
      REAL(dp), INTENT(IN), OPTIONAL ::  &
        in_array1_dp(:,:,:), in_array2_dp(:,:,:), in_array3_dp(:,:,:), &
        in_array4_dp(:,:,:), in_array5_dp(:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL ::  &
        ref_out_array1_dp(:,:,:), ref_out_array2_dp(:,:,:), &
        ref_out_array3_dp(:,:,:), ref_out_array4_dp(:,:,:), &
        ref_out_array5_dp(:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL :: out_array4d_dp(:,:,:,:)
      REAL(dp), INTENT(IN   ), OPTIONAL :: in_array4d_dp(:,:,:,:)
      REAL(dp), INTENT(INOUT), OPTIONAL :: ref_out_array4d_dp(:,:,:,:)
      REAL(sp), INTENT(INOUT), OPTIONAL ::  &
        out_array1_sp(:,:,:), out_array2_sp(:,:,:), out_array3_sp(:,:,:), &
        out_array4_sp(:,:,:), out_array5_sp(:,:,:)
      REAL(sp), INTENT(IN), OPTIONAL ::  &
        in_array1_sp(:,:,:), in_array2_sp(:,:,:), in_array3_sp(:,:,:), &
        in_array4_sp(:,:,:), in_array5_sp(:,:,:)
      REAL(sp), INTENT(INOUT), OPTIONAL ::  &
        ref_out_array1_sp(:,:,:), ref_out_array2_sp(:,:,:), &
        ref_out_array3_sp(:,:,:), ref_out_array4_sp(:,:,:), &
        ref_out_array5_sp(:,:,:)
      REAL(sp), INTENT(INOUT), OPTIONAL :: out_array4d_sp(:,:,:,:)
      REAL(sp), INTENT(IN   ), OPTIONAL :: in_array4d_sp(:,:,:,:)
      REAL(sp), INTENT(INOUT), OPTIONAL :: ref_out_array4d_sp(:,:,:,:)

      INTEGER, OPTIONAL, INTENT(IN) :: nshift
      CLASS(t_comm_pattern), POINTER, INTENT(INOUT) :: comm_pattern


      REAL(dp), ALLOCATABLE ::  &
        tmp_out_array1_dp(:,:,:), tmp_out_array2_dp(:,:,:), tmp_out_array3_dp(:,:,:), &
        tmp_out_array4_dp(:,:,:), tmp_out_array5_dp(:,:,:), tmp_out_array4d_dp(:,:,:,:)
      REAL(sp), ALLOCATABLE ::  &
        tmp_out_array1_sp(:,:,:), tmp_out_array2_sp(:,:,:), tmp_out_array3_sp(:,:,:), &
        tmp_out_array4_sp(:,:,:), tmp_out_array5_sp(:,:,:), tmp_out_array4d_sp(:,:,:,:)

      INTEGER :: nfields_dp, ndim2tot_dp, nfields_sp, ndim2tot_sp, ndim2, kshift

      kshift = 0
      IF (PRESENT(nshift)) kshift = nshift

      nfields_dp = 0
      ndim2tot_dp = 0
      IF (PRESENT(out_array4d_dp)) THEN
        ALLOCATE(tmp_out_array4d_dp(SIZE(out_array4d_dp,1),SIZE(out_array4d_dp,2),&
                                    SIZE(out_array4d_dp,3),SIZE(out_array4d_dp,4)))
        tmp_out_array4d_dp = out_array4d_dp
        nfields_dp = nfields_dp + SIZE(out_array4d_dp, 4)
        ndim2 = SIZE(out_array4d_dp, 2)
        ndim2tot_dp = ndim2tot_dp + &
          &        MERGE(1, ndim2 - kshift, ndim2 == 1) * SIZE(out_array4d_dp, 4)
      END IF
      IF (PRESENT(out_array1_dp)) THEN
        ALLOCATE(tmp_out_array1_dp(SIZE(out_array1_dp,1),SIZE(out_array1_dp,2),&
                                   SIZE(out_array1_dp,3)))
        tmp_out_array1_dp = out_array1_dp
        nfields_dp = nfields_dp + 1
        ndim2 = SIZE(out_array1_dp, 2)
        ndim2tot_dp = ndim2tot_dp + MERGE(1, ndim2 - kshift, ndim2 == 1)
        IF (PRESENT(out_array2_dp)) THEN
          ALLOCATE(tmp_out_array2_dp(SIZE(out_array2_dp,1),SIZE(out_array2_dp,2),&
                                     SIZE(out_array2_dp,3)))
          tmp_out_array2_dp = out_array2_dp
          nfields_dp = nfields_dp + 1
          ndim2 = SIZE(out_array2_dp, 2)
          ndim2tot_dp = ndim2tot_dp + MERGE(1, ndim2 - kshift, ndim2 == 1)
          IF (PRESENT(out_array3_dp)) THEN
            ALLOCATE(tmp_out_array3_dp(SIZE(out_array3_dp,1),SIZE(out_array3_dp,2),&
                                       SIZE(out_array3_dp,3)))
            tmp_out_array3_dp = out_array3_dp
            nfields_dp = nfields_dp + 1
            ndim2 = SIZE(out_array3_dp, 2)
            ndim2tot_dp = ndim2tot_dp + MERGE(1, ndim2 - kshift, ndim2 == 1)
            IF (PRESENT(out_array4_dp)) THEN
              ALLOCATE(tmp_out_array4_dp(SIZE(out_array4_dp,1),SIZE(out_array4_dp,2),&
                                         SIZE(out_array4_dp,3)))
              tmp_out_array4_dp = out_array4_dp
              nfields_dp = nfields_dp + 1
              ndim2 = SIZE(out_array4_dp, 2)
              ndim2tot_dp = ndim2tot_dp + MERGE(1, ndim2 - kshift, ndim2 == 1)
              IF (PRESENT(out_array5_dp)) THEN
                ALLOCATE(tmp_out_array5_dp(SIZE(out_array5_dp,1),SIZE(out_array5_dp,2),&
                                           SIZE(out_array5_dp,3)))
                tmp_out_array5_dp = out_array5_dp
                nfields_dp = nfields_dp + 1
                ndim2 = SIZE(out_array5_dp, 2)
                ndim2tot_dp = ndim2tot_dp + MERGE(1, ndim2 - kshift, ndim2 == 1)
              END IF
            END IF
          END IF
        END IF
      END IF
      nfields_sp = 0
      ndim2tot_sp = 0
      IF (PRESENT(out_array4d_sp)) THEN
        ALLOCATE(tmp_out_array4d_sp(SIZE(out_array4d_sp,1),SIZE(out_array4d_sp,2),&
                                    SIZE(out_array4d_sp,3),SIZE(out_array4d_sp,4)))
        tmp_out_array4d_sp = out_array4d_sp
        nfields_sp = nfields_sp + SIZE(out_array4d_sp, 4)
        ndim2 = SIZE(out_array4d_sp, 2)
        ndim2tot_sp = ndim2tot_sp + &
          &        MERGE(1, ndim2 - kshift, ndim2 == 1) * SIZE(out_array4d_sp, 4)
      END IF
      IF (PRESENT(out_array1_sp)) THEN
        ALLOCATE(tmp_out_array1_sp(SIZE(out_array1_sp,1),SIZE(out_array1_sp,2),&
                                   SIZE(out_array1_sp,3)))
        tmp_out_array1_sp = out_array1_sp
        nfields_sp = nfields_sp + 1
        ndim2 = SIZE(out_array1_sp, 2)
        ndim2tot_sp = ndim2tot_sp + MERGE(1, ndim2 - kshift, ndim2 == 1)
        IF (PRESENT(out_array2_sp)) THEN
          ALLOCATE(tmp_out_array2_sp(SIZE(out_array2_sp,1),SIZE(out_array2_sp,2),&
                                     SIZE(out_array2_sp,3)))
          tmp_out_array2_sp = out_array2_sp
          nfields_sp = nfields_sp + 1
          ndim2 = SIZE(out_array2_sp, 2)
          ndim2tot_sp = ndim2tot_sp + MERGE(1, ndim2 - kshift, ndim2 == 1)
          IF (PRESENT(out_array3_sp)) THEN
            ALLOCATE(tmp_out_array3_sp(SIZE(out_array3_sp,1),SIZE(out_array3_sp,2),&
                                       SIZE(out_array3_sp,3)))
            tmp_out_array3_sp = out_array3_sp
            nfields_sp = nfields_sp + 1
            ndim2 = SIZE(out_array3_sp, 2)
            ndim2tot_sp = ndim2tot_sp + MERGE(1, ndim2 - kshift, ndim2 == 1)
            IF (PRESENT(out_array4_sp)) THEN
              ALLOCATE(tmp_out_array4_sp(SIZE(out_array4_sp,1),SIZE(out_array4_sp,2),&
                                         SIZE(out_array4_sp,3)))
              tmp_out_array4_sp = out_array4_sp
              nfields_sp = nfields_sp + 1
              ndim2 = SIZE(out_array4_sp, 2)
              ndim2tot_sp = ndim2tot_sp + MERGE(1, ndim2 - kshift, ndim2 == 1)
              IF (PRESENT(out_array5_sp)) THEN
                ALLOCATE(tmp_out_array5_sp(SIZE(out_array5_sp,1),SIZE(out_array5_sp,2),&
                                           SIZE(out_array5_sp,3)))
                tmp_out_array5_sp = out_array5_sp
                nfields_sp = nfields_sp + 1
                ndim2 = SIZE(out_array5_sp, 2)
                ndim2tot_sp = ndim2tot_sp + MERGE(1, ndim2 - kshift, ndim2 == 1)
              END IF
            END IF
          END IF
        END IF
      END IF









      !$ACC DATA COPYIN(in_array1_dp) IF(present(in_array1_dp) .AND. i_am_accel_node)
      !$ACC DATA COPYIN(in_array2_dp) IF(present(in_array2_dp) .AND. i_am_accel_node)
      !$ACC DATA COPYIN(in_array3_dp) IF(present(in_array3_dp) .AND. i_am_accel_node)
      !$ACC DATA COPYIN(in_array4_dp) IF(present(in_array4_dp) .AND. i_am_accel_node)
      !$ACC DATA COPYIN(in_array5_dp) IF(present(in_array5_dp) .AND. i_am_accel_node)
      !$ACC DATA COPYIN(in_array4d_dp) IF(present(in_array4d_dp) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array1_dp) IF(present(out_array1_dp) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array2_dp) IF(present(out_array2_dp) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array3_dp) IF(present(out_array3_dp) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array4_dp) IF(present(out_array4_dp) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array5_dp) IF(present(out_array5_dp) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array4d_dp) IF(present(out_array4d_dp) .AND. i_am_accel_node)
      !$ACC DATA COPYIN(in_array1_sp) IF(present(in_array1_sp) .AND. i_am_accel_node)
      !$ACC DATA COPYIN(in_array2_sp) IF(present(in_array2_sp) .AND. i_am_accel_node)
      !$ACC DATA COPYIN(in_array3_sp) IF(present(in_array3_sp) .AND. i_am_accel_node)
      !$ACC DATA COPYIN(in_array4_sp) IF(present(in_array4_sp) .AND. i_am_accel_node)
      !$ACC DATA COPYIN(in_array5_sp) IF(present(in_array5_sp) .AND. i_am_accel_node)
      !$ACC DATA COPYIN(in_array4d_sp) IF(present(in_array4d_sp) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array1_sp) IF(present(out_array1_sp) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array2_sp) IF(present(out_array2_sp) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array3_sp) IF(present(out_array3_sp) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array4_sp) IF(present(out_array4_sp) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array5_sp) IF(present(out_array5_sp) .AND. i_am_accel_node)
      !$ACC DATA COPY(out_array4d_sp) IF(present(out_array4d_sp) .AND. i_am_accel_node)
      IF ((nfields_dp > 0) .OR. (nfields_sp > 0)) THEN
        CALL exchange_data_mult_mixprec( &
          p_pat=comm_pattern, nfields_dp=nfields_dp, ndim2tot_dp=ndim2tot_dp, &
          recv1_dp=out_array1_dp, send1_dp=in_array1_dp, &
          recv2_dp=out_array2_dp, send2_dp=in_array2_dp, &
          recv3_dp=out_array3_dp, send3_dp=in_array3_dp, &
          recv4_dp=out_array4_dp, send4_dp=in_array4_dp, &
          recv5_dp=out_array5_dp, send5_dp=in_array5_dp, &
          recv4d_dp=out_array4d_dp, send4d_dp=in_array4d_dp, &
          nfields_sp=nfields_sp, ndim2tot_sp=ndim2tot_sp, &
          recv1_sp=out_array1_sp, send1_sp=in_array1_sp, &
          recv2_sp=out_array2_sp, send2_sp=in_array2_sp, &
          recv3_sp=out_array3_sp, send3_sp=in_array3_sp, &
          recv4_sp=out_array4_sp, send4_sp=in_array4_sp, &
          recv5_sp=out_array5_sp, send5_sp=in_array5_sp, &
          recv4d_sp=out_array4d_sp, send4d_sp=in_array4d_sp, &
          nshift=kshift)
      END IF
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA
      !$ACC END DATA





      IF (PRESENT(out_array1_dp)) THEN
        IF (ANY(out_array1_dp /= ref_out_array1_dp)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result mixprec dp_3d_1"
          CALL finish(method_name, message_text)
        END IF
        out_array1_dp = tmp_out_array1_dp
      END IF
      IF (PRESENT(out_array2_dp)) THEN
        IF (ANY(out_array2_dp /= ref_out_array2_dp)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result mixprec dp_3d_2"
          CALL finish(method_name, message_text)
        END IF
        out_array2_dp = tmp_out_array2_dp
      END IF
      IF (PRESENT(out_array3_dp)) THEN
        IF (ANY(out_array3_dp /= ref_out_array3_dp)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result mixprec dp_3d_3"
          CALL finish(method_name, message_text)
        END IF
        out_array3_dp = tmp_out_array3_dp
      END IF
      IF (PRESENT(out_array4_dp)) THEN
        IF (ANY(out_array4_dp /= ref_out_array4_dp)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result mixprec dp_3d_4"
          CALL finish(method_name, message_text)
        END IF
        out_array4_dp = tmp_out_array4_dp
      END IF
      IF (PRESENT(out_array5_dp)) THEN
        IF (ANY(out_array5_dp /= ref_out_array5_dp)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result mixprec dp_3d_5"
          CALL finish(method_name, message_text)
        END IF
        out_array5_dp = tmp_out_array5_dp
      END IF
      IF (PRESENT(out_array4d_dp)) THEN
        IF (ANY(out_array4d_dp /= ref_out_array4d_dp)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result mixprec dp_4d"
          CALL finish(method_name, message_text)
        END IF
        out_array4d_dp = tmp_out_array4d_dp
      END IF
      IF (PRESENT(out_array1_sp)) THEN
        IF (ANY(out_array1_sp /= ref_out_array1_sp)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result mixprec sp_3d_1"
          CALL finish(method_name, message_text)
        END IF
        out_array1_sp = tmp_out_array1_sp
      END IF
      IF (PRESENT(out_array2_sp)) THEN
        IF (ANY(out_array2_sp /= ref_out_array2_sp)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result mixprec sp_3d_2"
          CALL finish(method_name, message_text)
        END IF
        out_array2_sp = tmp_out_array2_sp
      END IF
      IF (PRESENT(out_array3_sp)) THEN
        IF (ANY(out_array3_sp /= ref_out_array3_sp)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result mixprec sp_3d_3"
          CALL finish(method_name, message_text)
        END IF
        out_array3_sp = tmp_out_array3_sp
      END IF
      IF (PRESENT(out_array4_sp)) THEN
        IF (ANY(out_array4_sp /= ref_out_array4_sp)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result mixprec sp_3d_4"
          CALL finish(method_name, message_text)
        END IF
        out_array4_sp = tmp_out_array4_sp
      END IF
      IF (PRESENT(out_array5_sp)) THEN
        IF (ANY(out_array5_sp /= ref_out_array5_sp)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result mixprec sp_3d_5"
          CALL finish(method_name, message_text)
        END IF
        out_array5_sp = tmp_out_array5_sp
      END IF
      IF (PRESENT(out_array4d_sp)) THEN
        IF (ANY(out_array4d_sp /= ref_out_array4d_sp)) THEN
          WRITE(message_text,'(a,i0)') "Wrong exchange result mixprec sp_4d"
          CALL finish(method_name, message_text)
        END IF
        out_array4d_sp = tmp_out_array4d_sp
      END IF
    END SUBROUTINE check_exchange_mixprec

  END SUBROUTINE exchange_communication_testbed

  SUBROUTINE exchange_communication_grf_testbed(test_gpu)

    LOGICAL, OPTIONAL, INTENT(IN) :: test_gpu

    CHARACTER(*), PARAMETER :: method_name = &
      "mo_test_communication:exchange_communication_grf_testbed"

    INTEGER, ALLOCATABLE :: owner_local_src(:), owner_local_dst(:,:), &
      &                     glb_index_src(:), glb_index_dst(:)
    REAL(wp), ALLOCATABLE :: ref_recv_add(:,:)
    INTEGER :: i, j, n, local_size_src, local_size_dst, global_size, count
    TYPE(t_glb2loc_index_lookup) :: send_glb2loc_index
    TYPE(t_p_comm_pattern) :: comm_pattern(4)
    CLASS(t_comm_pattern_collection), POINTER :: comm_pattern_collection

    REAL(wp) :: recv1(nproma,2,18), recv2(nproma,4,18), &
      &         recv3(nproma,6,18), recv4(nproma,8,18), &
      &         recv5(nproma,10,18), recv6(nproma,12,18), &
      &         recv4d1(nproma,4,18,6), recv4d2(nproma,8,18,6)
    REAL(wp) :: send1(2,12*nproma,4), send2(4,12*nproma,4), &
      &         send3(6,12*nproma,4), send4(8,12*nproma,4), &
      &         send5(10,12*nproma,4), send6(12,12*nproma,4), &
      &         send4d1(4,12*nproma,4,6), send4d2(8,12*nproma,4,6)
    REAL(wp) :: ref_recv1(nproma,2,18), ref_recv2(nproma,4,18), &
      &         ref_recv3(nproma,6,18), ref_recv4(nproma,8,18), &
      &         ref_recv5(nproma,10,18), ref_recv6(nproma,12,18), &
      &         ref_recv4d1(nproma,4,18,6), ref_recv4d2(nproma,8,18,6)

    !generate communication patterns
    local_size_src = 12 * nproma
    local_size_dst = 18 * nproma
    global_size = local_size_src * p_n_work
    ALLOCATE(owner_local_src(local_size_src), glb_index_src(local_size_src), &
      &      owner_local_dst(local_size_dst,4), glb_index_dst(local_size_dst), &
      &      ref_recv_add(nproma,18))

    owner_local_src = p_pe_work
    owner_local_dst = -1
    DO i = 1, 4
      owner_local_dst(i:3*nproma:4,i) = MOD(p_pe_work + p_n_work - 1, p_n_work)
      owner_local_dst(3*nproma+i:15*nproma:4,i) = p_pe_work
      owner_local_dst(15*nproma+i:18*nproma:4,i) = MOD(p_pe_work + 1, p_n_work)
    END DO
    DO i = 1, local_size_src
      glb_index_src(i) = local_size_src * p_pe_work + i
    END DO
    DO i = 1, local_size_dst
      glb_index_dst(i) = MOD(local_size_src * p_pe_work - 3 * nproma - 1 + i + &
      &                      global_size, global_size) + 1
    END DO
    CALL init_glb2loc_index_lookup(send_glb2loc_index, global_size)
    CALL set_inner_glb_index(send_glb2loc_index, glb_index_src, &
      &                      (/(i, i=1, local_size_src)/))

    DO i = 1, 4
      CALL setup_comm_pattern(local_size_dst, owner_local_dst(:,i), &
        &                     glb_index_dst, send_glb2loc_index, &
        &                     local_size_src, owner_local_src, glb_index_src, &
        &                     comm_pattern(i)%p)
    END DO

    CALL setup_comm_pattern_collection(comm_pattern, comm_pattern_collection)

    recv1 = -1
    recv2 = -1
    recv3 = -1
    recv4 = -1
    recv5 = -1
    recv6 = -1
    recv4d1 = -1
    recv4d2 = -1

    count = 0
    DO n = 1, SIZE(send1, 1)
      DO i = 1, local_size_src
        send1(n,i,:) = (/0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp/) + count + i + &
          &            p_pe_work * local_size_src
      END DO
      count = count + local_size_src
    END DO
    DO n = 1, SIZE(send2, 1)
      DO i = 1, local_size_src
        send2(n,i,:) = (/0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp/) + count + i + &
          &            p_pe_work * local_size_src
      END DO
      count = count + local_size_src
    END DO
    DO n = 1, SIZE(send3, 1)
      DO i = 1, local_size_src
        send3(n,i,:) = (/0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp/) + count + i + &
          &            p_pe_work * local_size_src
      END DO
      count = count + local_size_src
    END DO
    DO n = 1, SIZE(send4, 1)
      DO i = 1, local_size_src
        send4(n,i,:) = (/0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp/) + count + i + &
          &            p_pe_work * local_size_src
      END DO
      count = count + local_size_src
    END DO
    DO n = 1, SIZE(send5, 1)
      DO i = 1, local_size_src
        send5(n,i,:) = (/0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp/) + count + i + &
          &            p_pe_work * local_size_src
      END DO
      count = count + local_size_src
    END DO
    DO n = 1, SIZE(send6, 1)
      DO i = 1, local_size_src
        send6(n,i,:) = (/0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp/) + count + i + &
          &            p_pe_work * local_size_src
      END DO
      count = count + local_size_src
    END DO

    count = 0
    DO j = 1, SIZE(send4d1, 4)
      DO n = 1, SIZE(send4d1, 1)
        DO i = 1, local_size_src
          send4d1(n,i,:, j) = (/0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp/) + count + i + &
            &                 p_pe_work * local_size_src
        END DO
        count = count + local_size_src
      END DO
    END DO
    DO j = 1, SIZE(send4d2, 4)
      DO n = 1, SIZE(send4d2, 1)
        DO i = 1, local_size_src
          send4d2(n,i,:, j) = (/0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp/) + count + i + &
          &                   p_pe_work * local_size_src
        END DO
        count = count + local_size_src
      END DO
    END DO

    ref_recv_add(1::4,:) = 0.1_wp
    ref_recv_add(2::4,:) = 0.2_wp
    ref_recv_add(3::4,:) = 0.3_wp
    ref_recv_add(4::4,:) = 0.4_wp

    count = 0
    DO n = 1, SIZE(ref_recv1, 2)
      ref_recv1(:,n,:) = RESHAPE(glb_index_dst + count,(/nproma,18/)) + &
        &                ref_recv_add
      count = count + local_size_src
    END DO
    DO n = 1, SIZE(ref_recv2, 2)
      ref_recv2(:,n,:) = RESHAPE(glb_index_dst + count,(/nproma,18/)) + &
        &                ref_recv_add
      count = count + local_size_src
    END DO
    DO n = 1, SIZE(ref_recv3, 2)
      ref_recv3(:,n,:) = RESHAPE(glb_index_dst + count,(/nproma,18/)) + &
        &                ref_recv_add
      count = count + local_size_src
    END DO
    DO n = 1, SIZE(ref_recv4, 2)
      ref_recv4(:,n,:) = RESHAPE(glb_index_dst + count,(/nproma,18/)) + &
        &                ref_recv_add
      count = count + local_size_src
    END DO
    DO n = 1, SIZE(ref_recv5, 2)
      ref_recv5(:,n,:) = RESHAPE(glb_index_dst + count,(/nproma,18/)) + &
        &                ref_recv_add
      count = count + local_size_src
    END DO
    DO n = 1, SIZE(ref_recv6, 2)
      ref_recv6(:,n,:) = RESHAPE(glb_index_dst + count,(/nproma,18/)) + &
        &                ref_recv_add
      count = count + local_size_src
    END DO

    count = 0
    DO i = 1, SIZE(ref_recv4d1, 4)
      DO n = 1, SIZE(ref_recv4d1, 2)
        ref_recv4d1(:,n,:,i) = RESHAPE(glb_index_dst + count,(/nproma,18/)) + &
        &                ref_recv_add
        count = count + local_size_src
      END DO
    END DO
    DO i = 1, SIZE(ref_recv4d2, 4)
      DO n = 1, SIZE(ref_recv4d2, 2)
        ref_recv4d2(:,n,:,i) = RESHAPE(glb_index_dst + count,(/nproma,18/)) + &
        &                ref_recv_add
        count = count + local_size_src
      END DO
    END DO

    CALL check_exchange_grf(comm_pattern_collection, recv1, send1, ref_recv1, &
      &                     recv2, send2, ref_recv2, recv3, send3, ref_recv3, &
      &                     recv4, send4, ref_recv4, recv5, send5, ref_recv5, &
      &                     recv6, send6, ref_recv6, recv4d1, send4d1, &
      &                     ref_recv4d1, recv4d2, send4d2, ref_recv4d2)

    CALL delete_comm_pattern_collection(comm_pattern_collection)

    CONTAINS

    SUBROUTINE check_exchange_grf( &
      p_pat_coll, recv1, send1, ref_recv1, recv2, send2, ref_recv2, recv3, &
      send3, ref_recv3, recv4, send4, ref_recv4, recv5, send5, ref_recv5, &
      recv6, send6, ref_recv6, recv4d1, send4d1, ref_recv4d1, recv4d2, &
      send4d2, ref_recv4d2)

      CLASS(t_comm_pattern_collection), POINTER, INTENT(INOUT) :: p_pat_coll

      ! recv3d (nproma,nlev,blk)
      ! recv4d (nproma,nlev,blk,nfield)
      REAL(wp), INTENT(INOUT) :: recv1(:,:,:), recv2(:,:,:), recv3(:,:,:), &
        &                        recv4(:,:,:), recv5(:,:,:), recv6(:,:,:), &
        &                        recv4d1(:,:,:,:), recv4d2(:,:,:,:)
      ! send3d (nlev,i,npat)
      ! send4d (nlev,i,npat,nfield)
      REAL(wp), INTENT(IN) ::  send1(:,:,:), send2(:,:,:), send3(:,:,:), &
        &                      send4(:,:,:), send5(:,:,:), send6(:,:,:), &
        &                      send4d1(:,:,:,:), send4d2(:,:,:,:)

      REAL(wp), INTENT(IN) :: ref_recv1(:,:,:), ref_recv2(:,:,:), &
        &                     ref_recv3(:,:,:), ref_recv4(:,:,:), &
        &                     ref_recv5(:,:,:), ref_recv6(:,:,:), &
        &                     ref_recv4d1(:,:,:,:), ref_recv4d2(:,:,:,:)

      REAL(wp), ALLOCATABLE :: tmp_recv1(:,:,:), tmp_recv2(:,:,:), &
        &                      tmp_recv3(:,:,:), tmp_recv4(:,:,:), &
        &                      tmp_recv5(:,:,:), tmp_recv6(:,:,:), &
        &                      tmp_recv4d1(:,:,:,:), tmp_recv4d2(:,:,:,:)

      INTEGER :: nfields, ndim2tot, i









      ALLOCATE(tmp_recv1(SIZE(recv1,1),SIZE(recv1,2),SIZE(recv1,3)), &
        &      tmp_recv2(SIZE(recv2,1),SIZE(recv2,2),SIZE(recv2,3)), &
        &      tmp_recv3(SIZE(recv3,1),SIZE(recv3,2),SIZE(recv3,3)), &
        &      tmp_recv4(SIZE(recv4,1),SIZE(recv4,2),SIZE(recv4,3)), &
        &      tmp_recv5(SIZE(recv5,1),SIZE(recv5,2),SIZE(recv5,3)), &
        &      tmp_recv6(SIZE(recv6,1),SIZE(recv6,2),SIZE(recv6,3)), &
        &      tmp_recv4d1(SIZE(recv4d1,1),SIZE(recv4d1,2),&
        &                  SIZE(recv4d1,3),SIZE(recv4d1,4)), &
        &      tmp_recv4d2(SIZE(recv4d2,1),SIZE(recv4d2,2),&
        &                  SIZE(recv4d2,3),SIZE(recv4d2,4)))

      tmp_recv1 = recv1
      tmp_recv2 = recv2
      tmp_recv3 = recv3
      tmp_recv4 = recv4
      tmp_recv5 = recv5
      tmp_recv6 = recv6
      tmp_recv4d1 = recv4d1
      tmp_recv4d2 = recv4d2

      i = 1
      nfields = 1
      ndim2tot = SIZE(recv1, 2)
      recv1 = tmp_recv1
      !$ACC DATA COPY(recv1) COPYIN(send1) IF(i_am_accel_node)
      CALL exchange_data_grf(p_pat_coll=p_pat_coll, nfields=nfields, &
        &                    ndim2tot=ndim2tot, recv1=recv1, send1=send1)
      !$ACC END DATA
      IF (ANY(recv1 /= ref_recv1)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result grf", i
        CALL finish(method_name, message_text)
      END IF

      i = i + 1
      nfields = 2
      ndim2tot = SIZE(recv1, 2) + SIZE(recv2, 2)
      recv1 = tmp_recv1
      recv2 = tmp_recv2
      !$ACC DATA COPY(recv1, recv2) COPYIN(send1, send2) IF(i_am_accel_node)
      CALL exchange_data_grf(p_pat_coll=p_pat_coll, nfields=nfields, &
        &                    ndim2tot=ndim2tot, recv1=recv1, send1=send1, &
        &                    recv2=recv2, send2=send2)
      !$ACC END DATA
      IF (ANY(recv1 /= ref_recv1) .OR. ANY(recv2 /= ref_recv2)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result grf", i
        CALL finish(method_name, message_text)
      END IF

      i = i + 1
      nfields = 3
      ndim2tot = SIZE(recv1, 2) + SIZE(recv2, 2) + SIZE(recv3, 2)
      recv1 = tmp_recv1
      recv2 = tmp_recv2
      recv3 = tmp_recv3
      !$ACC DATA COPY(recv1, recv2, recv3) &
      !$ACC   COPYIN(send1, send2, send3) IF(i_am_accel_node)
      CALL exchange_data_grf(p_pat_coll=p_pat_coll, nfields=nfields, &
        &                    ndim2tot=ndim2tot, recv1=recv1, send1=send1, &
        &                    recv2=recv2, send2=send2, recv3=recv3, &
        &                    send3=send3)
      !$ACC END DATA
      IF (ANY(recv1 /= ref_recv1) .OR. ANY(recv2 /= ref_recv2) .OR. &
        & ANY(recv3 /= ref_recv3)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result grf", i
        CALL finish(method_name, message_text)
      END IF

      i = i + 1
      nfields = 4
      ndim2tot = SIZE(recv1, 2) + SIZE(recv2, 2) + SIZE(recv3, 2) + &
        &        SIZE(recv4, 2)
      recv1 = tmp_recv1
      recv2 = tmp_recv2
      recv3 = tmp_recv3
      recv4 = tmp_recv4
      !$ACC DATA COPY(recv1, recv2, recv3, recv4) &
      !$ACC   COPYIN(send1, send2, send3, send4) IF(i_am_accel_node)
      CALL exchange_data_grf(p_pat_coll=p_pat_coll, nfields=nfields, &
        &                    ndim2tot=ndim2tot, recv1=recv1, send1=send1, &
        &                    recv2=recv2, send2=send2, recv3=recv3, &
        &                    send3=send3, recv4=recv4, send4=send4)
      !$ACC END DATA
      IF (ANY(recv1 /= ref_recv1) .OR. ANY(recv2 /= ref_recv2) .OR. &
        & ANY(recv3 /= ref_recv3) .OR. ANY(recv4 /= ref_recv4)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result grf", i
        CALL finish(method_name, message_text)
      END IF

      i = i + 1
      nfields = 5
      ndim2tot = SIZE(recv1, 2) + SIZE(recv2, 2) + SIZE(recv3, 2) + &
        &        SIZE(recv4, 2) + SIZE(recv5, 2)
      recv1 = tmp_recv1
      recv2 = tmp_recv2
      recv3 = tmp_recv3
      recv4 = tmp_recv4
      recv5 = tmp_recv5
      !$ACC DATA COPY(recv1, recv2, recv3, recv4, recv5) &
      !$ACC   COPYIN(send1, send2, send3, send4, send5) IF(i_am_accel_node)
      CALL exchange_data_grf(p_pat_coll=p_pat_coll, nfields=nfields, &
        &                    ndim2tot=ndim2tot, recv1=recv1, send1=send1, &
        &                    recv2=recv2, send2=send2, recv3=recv3, &
        &                    send3=send3, recv4=recv4, send4=send4, &
        &                    recv5=recv5, send5=send5)
      !$ACC END DATA
      IF (ANY(recv1 /= ref_recv1) .OR. ANY(recv2 /= ref_recv2) .OR. &
        & ANY(recv3 /= ref_recv3) .OR. ANY(recv4 /= ref_recv4) .OR. &
        & ANY(recv5 /= ref_recv5)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result grf", i
        CALL finish(method_name, message_text)
      END IF

      i = i + 1
      nfields = 6
      ndim2tot = SIZE(recv1, 2) + SIZE(recv2, 2) + SIZE(recv3, 2) + &
        &        SIZE(recv4, 2) + SIZE(recv5, 2) + SIZE(recv6, 2)
      recv1 = tmp_recv1
      recv2 = tmp_recv2
      recv3 = tmp_recv3
      recv4 = tmp_recv4
      recv5 = tmp_recv5
      recv6 = tmp_recv6
      !$ACC DATA COPY(recv1, recv2, recv3, recv4, recv5, recv6) &
      !$ACC   COPYIN(send1, send2, send3, send4, send5, send6) IF(i_am_accel_node)
      CALL exchange_data_grf(p_pat_coll=p_pat_coll, nfields=nfields, &
        &                    ndim2tot=ndim2tot, recv1=recv1, send1=send1, &
        &                    recv2=recv2, send2=send2, recv3=recv3, &
        &                    send3=send3, recv4=recv4, send4=send4, &
        &                    recv5=recv5, send5=send5, recv6=recv6, &
        &                    send6=send6)
      !$ACC END DATA
      IF (ANY(recv1 /= ref_recv1) .OR. ANY(recv2 /= ref_recv2) .OR. &
        & ANY(recv3 /= ref_recv3) .OR. ANY(recv4 /= ref_recv4) .OR. &
        & ANY(recv5 /= ref_recv5) .OR. ANY(recv6 /= ref_recv6)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result grf", i
        CALL finish(method_name, message_text)
      END IF

      i = i + 1
      nfields = SIZE(recv4d1, 4)
      ndim2tot = SIZE(recv4d1, 4) * SIZE(recv4d1, 2)
      recv4d1 = tmp_recv4d1
      !$ACC DATA COPY(recv4d1) COPYIN(send4d1) IF(i_am_accel_node)
      CALL exchange_data_grf(p_pat_coll=p_pat_coll, nfields=nfields, &
        &                    ndim2tot=ndim2tot, recv4d1=recv4d1, &
        &                    send4d1=send4d1)
      !$ACC END DATA
      IF (ANY(recv4d1 /= ref_recv4d1)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result grf", i
        CALL finish(method_name, message_text)
      END IF

      i = i + 1
      nfields = 2 * SIZE(recv4d1, 4)
      ndim2tot = SIZE(recv4d1, 4) * SIZE(recv4d1, 2) + &
        &        SIZE(recv4d2, 4) * SIZE(recv4d2, 2)
      recv4d1 = tmp_recv4d1
      recv4d2 = tmp_recv4d2
      !$ACC DATA COPY(recv4d1, recv4d2) COPYIN(send4d1, send4d2) IF(i_am_accel_node)
      CALL exchange_data_grf(p_pat_coll=p_pat_coll, nfields=nfields, &
        &                    ndim2tot=ndim2tot, recv4d1=recv4d1, &
        &                    send4d1=send4d1, recv4d2=recv4d2, send4d2=send4d2)
      !$ACC END DATA
      IF (ANY(recv4d1 /= ref_recv4d1) .OR. ANY(recv4d2 /= ref_recv4d2)) THEN
        WRITE(message_text,'(a,i0)') "Wrong exchange result grf", i
        CALL finish(method_name, message_text)
      END IF
    END SUBROUTINE check_exchange_grf

  END SUBROUTINE exchange_communication_grf_testbed

  SUBROUTINE gather_communication_testbed()

    INTEGER :: i, j

    INTEGER :: local_size, p_n_intercomm_remote, global_size
    INTEGER, ALLOCATABLE :: owner_local(:), glb_index(:)
    TYPE(t_comm_gather_pattern), TARGET :: gather_pattern
    TYPE(t_comm_allgather_pattern) :: allgather_pattern
    LOGICAL :: disable_consistency_check

    INTEGER :: nlev
    REAL(wp), ALLOCATABLE :: in_array_r_1d(:,:), in_array_r_2d(:,:,:)
    INTEGER, ALLOCATABLE :: in_array_i_1d(:,:), in_array_i_2d(:,:,:)
    REAL(wp), ALLOCATABLE :: out_array_r_1d(:), out_array_r_2d(:,:)
    REAL(wp), ALLOCATABLE :: ref_out_array_r_1d(:), ref_out_array_r_2d(:,:)
    INTEGER, ALLOCATABLE :: out_array_i_1d(:), out_array_i_2d(:,:)
    INTEGER, ALLOCATABLE :: ref_out_array_i_1d(:), ref_out_array_i_2d(:,:)
    REAL(wp) :: fill_value
    INTEGER :: p_comm_work_backup, p_pe_work_backup, p_n_work_backup
    INTEGER :: intercomm, intercomm_key, p_comm_work_new, ierr
    LOGICAL :: is_active

    CHARACTER(*), PARAMETER :: method_name = &
      "mo_test_communication:gather_communication_testbed"

    !---------------------------------------------------------------------------
    ! simple test in which each process has its own local contiguous part of the
    ! global array
    !---------------------------------------------------------------------------
    ! generate gather pattern
    local_size = 10 * nproma
    global_size = p_n_work * local_size
    ALLOCATE(owner_local(local_size), glb_index(local_size))
    DO i = 1, local_size
      owner_local(i) = p_pe_work
      glb_index(i) = local_size * p_pe_work + i
    END DO
    disable_consistency_check = .FALSE.
    CALL setup_comm_gather_pattern(global_size, owner_local, glb_index, &
      &                            gather_pattern)

    ! initialise in- and reference out data
    nlev = 5
    fill_value = -1
    ALLOCATE(in_array_r_1d(nproma, local_size / nproma), &
      &      in_array_r_2d(nproma, nlev, local_size / nproma), &
      &      in_array_i_1d(nproma, local_size / nproma), &
      &      in_array_i_2d(nproma, nlev, local_size / nproma))
    DO i = 0, local_size-1
      in_array_r_1d(MOD(i,nproma)+1, i/nproma+1) = p_pe_work * local_size + i
      in_array_i_1d(MOD(i,nproma)+1, i/nproma+1) = p_pe_work * local_size + i
      DO j = 1, nlev
        in_array_r_2d(MOD(i,nproma)+1, j, i/nproma+1) = (j - 1) * global_size + &
          &                                             p_pe_work * local_size + i
        in_array_i_2d(MOD(i,nproma)+1, j, i/nproma+1) = (j - 1) * global_size + &
          &                                             p_pe_work * local_size + i
      END DO
    END DO
    IF (p_pe_work == 0) THEN
      ALLOCATE(out_array_r_1d(global_size), &
        &      out_array_r_2d(global_size, nlev), &
        &      out_array_i_1d(global_size), &
        &      out_array_i_2d(global_size, nlev), &
        &      ref_out_array_r_1d(global_size), &
        &      ref_out_array_r_2d(global_size, nlev), &
        &      ref_out_array_i_1d(global_size), &
        &      ref_out_array_i_2d(global_size, nlev))
      DO i = 0, global_size - 1
        ref_out_array_r_1d(i+1) = i
        ref_out_array_i_1d(i+1) = i
        DO j = 1, nlev
          ref_out_array_r_2d(i+1, j) = (j - 1) * global_size + i
          ref_out_array_i_2d(i+1, j) = (j - 1) * global_size + i
        END DO
      END DO
    ELSE
      ALLOCATE(out_array_r_1d(0), &
        &      out_array_r_2d(0, 0), &
        &      out_array_i_1d(0), &
        &      out_array_i_2d(0, 0), &
        &      ref_out_array_r_1d(0), &
        &      ref_out_array_r_2d(0, 0), &
        &      ref_out_array_i_1d(0), &
        &      ref_out_array_i_2d(0, 0))
    END IF

    ! check gather pattern
    CALL check_exchange_gather(in_array_r_1d, in_array_r_2d, &
      &                        in_array_i_1d, in_array_i_2d, &
      &                        out_array_r_1d, out_array_r_2d, &
      &                        ref_out_array_r_1d, ref_out_array_r_2d, &
      &                        out_array_i_1d, out_array_i_2d, &
      &                        ref_out_array_i_1d, ref_out_array_i_2d, &
      &                        gather_pattern, 3649)
    CALL check_exchange_gather(in_array_r_1d, in_array_r_2d, &
      &                        in_array_i_1d, in_array_i_2d, &
      &                        out_array_r_1d, out_array_r_2d, &
      &                        ref_out_array_r_1d, ref_out_array_r_2d, &
      &                        out_array_i_1d, out_array_i_2d, &
      &                        ref_out_array_i_1d, ref_out_array_i_2d, &
      &                        gather_pattern, 3656, fill_value, INT(fill_value))

    ! delete gather pattern and other arrays
    DEALLOCATE(in_array_r_1d, in_array_r_2d, in_array_i_1d, in_array_i_2d)
    DEALLOCATE(out_array_r_1d, out_array_r_2d, out_array_i_1d, out_array_i_2d, &
      &        ref_out_array_r_1d, ref_out_array_r_2d, &
      &        ref_out_array_i_1d, ref_out_array_i_2d)
    CALL delete_comm_gather_pattern(gather_pattern)
    DEALLOCATE(owner_local, glb_index)

    !---------------------------------------------------------------------------
    ! simple test in which each process has its own local contiguous part of the
    ! global array plus some overlap with neighbouring processes
    !---------------------------------------------------------------------------
    IF (p_n_work > 1) THEN
      ! generate gather pattern
      local_size = 12 * nproma
      global_size = p_n_work * 10 * nproma
      ALLOCATE(owner_local(local_size), glb_index(local_size))
      DO i = 1, local_size
        owner_local(i) = p_pe_work
        glb_index(i) =  &
          MOD(global_size + (p_pe_work * 10 - 1) * nproma + i - 1, global_size) + 1
      END DO
      DO i = 1, nproma
        owner_local(i) = MOD(p_n_work + p_pe_work - 1, p_n_work)
        owner_local(11 * nproma + i) = MOD(p_pe_work + 1, p_n_work)
      END DO
      disable_consistency_check = .FALSE.

      CALL setup_comm_gather_pattern(global_size, owner_local, glb_index, &
        &                            gather_pattern)
      ! initialise in- and reference out data
      nlev = 5
      fill_value = -1
      ALLOCATE(in_array_r_1d(nproma, local_size / nproma), &
        &      in_array_r_2d(nproma, nlev, local_size / nproma), &
        &      in_array_i_1d(nproma, local_size / nproma), &
        &      in_array_i_2d(nproma, nlev, local_size / nproma))

      DO i = 0, local_size-1
        in_array_r_1d(MOD(i,nproma)+1, i/nproma+1) = &
          MOD(global_size + (p_pe_work * 10 - 1) * nproma + i, global_size)
        in_array_i_1d(MOD(i,nproma)+1, i/nproma+1) = &
          MOD(global_size + (p_pe_work * 10 - 1) * nproma + i, global_size)
        DO j = 1, nlev
          in_array_r_2d(MOD(i,nproma)+1, j, i/nproma+1) = &
            in_array_r_1d(MOD(i,nproma)+1, i/nproma+1) + (j - 1) * global_size
          in_array_i_2d(MOD(i,nproma)+1, j, i/nproma+1) = &
            in_array_i_1d(MOD(i,nproma)+1, i/nproma+1) + (j - 1) * global_size
        END DO
      END DO
      IF (p_pe_work == 0) THEN
        ALLOCATE(out_array_r_1d(global_size), &
          &      out_array_r_2d(global_size, nlev), &
          &      out_array_i_1d(global_size), &
          &      out_array_i_2d(global_size, nlev), &
          &      ref_out_array_r_1d(global_size), &
          &      ref_out_array_r_2d(global_size, nlev), &
          &      ref_out_array_i_1d(global_size), &
          &      ref_out_array_i_2d(global_size, nlev))
        DO i = 0, global_size - 1
          ref_out_array_r_1d(i+1) = i
          ref_out_array_i_1d(i+1) = i
          DO j = 1, nlev
            ref_out_array_r_2d(i+1, j) = (j - 1) * global_size + i
            ref_out_array_i_2d(i+1, j) = (j - 1) * global_size + i
          END DO
        END DO
      ELSE
        ALLOCATE(out_array_r_1d(0), &
          &      out_array_r_2d(0, 0), &
          &      out_array_i_1d(0), &
          &      out_array_i_2d(0, 0), &
          &      ref_out_array_r_1d(0), &
          &      ref_out_array_r_2d(0, 0), &
          &      ref_out_array_i_1d(0), &
          &      ref_out_array_i_2d(0, 0))
      END IF

      ! check gather pattern
      CALL check_exchange_gather(in_array_r_1d, in_array_r_2d, &
        &                        in_array_i_1d, in_array_i_2d, &
        &                        out_array_r_1d, out_array_r_2d, &
        &                        ref_out_array_r_1d, ref_out_array_r_2d, &
        &                        out_array_i_1d, out_array_i_2d, &
        &                        ref_out_array_i_1d, ref_out_array_i_2d, &
        &                        gather_pattern, 3743)
      CALL check_exchange_gather(in_array_r_1d, in_array_r_2d, &
        &                        in_array_i_1d, in_array_i_2d, &
        &                        out_array_r_1d, out_array_r_2d, &
        &                        ref_out_array_r_1d, ref_out_array_r_2d, &
        &                        out_array_i_1d, out_array_i_2d, &
        &                        ref_out_array_i_1d, ref_out_array_i_2d, &
        &                        gather_pattern, 3750, fill_value, INT(fill_value))

      ! delete gather pattern and other arrays
      DEALLOCATE(in_array_r_1d, in_array_r_2d, in_array_i_1d, in_array_i_2d)
      DEALLOCATE(out_array_r_1d, out_array_r_2d, out_array_i_1d, out_array_i_2d, &
        &        ref_out_array_r_1d, ref_out_array_r_2d, &
        &        ref_out_array_i_1d, ref_out_array_i_2d)
      CALL delete_comm_gather_pattern(gather_pattern)
      DEALLOCATE(owner_local, glb_index)
    END IF ! (p_n_work > 1)

    !---------------------------------------------------------------------------
    ! test in which only odd numbered global indices are owned by processes
    !---------------------------------------------------------------------------
    ! generate gather pattern
    local_size = 10 * nproma
    global_size = p_n_work * local_size
    ALLOCATE(owner_local(local_size), glb_index(local_size))
    DO i = 1, local_size, 2
      owner_local(i) = p_pe_work
    END DO
    DO i = 2, local_size, 2
      owner_local(i) = -1
    END DO
    DO i = 1, local_size
      glb_index(i) = local_size * p_pe_work + i
    END DO
    disable_consistency_check = .TRUE.

    CALL setup_comm_gather_pattern(global_size, owner_local, glb_index, &
      &                            gather_pattern, disable_consistency_check)

    ! initialise in- and reference out data
    nlev = 5
    fill_value = -1
    ALLOCATE(in_array_r_1d(nproma, local_size / nproma), &
      &      in_array_r_2d(nproma, nlev, local_size / nproma), &
      &      in_array_i_1d(nproma, local_size / nproma), &
      &      in_array_i_2d(nproma, nlev, local_size / nproma))
    DO i = 0, local_size-1
      in_array_r_1d(MOD(i,nproma)+1, i/nproma+1) = p_pe_work * local_size + i
      in_array_i_1d(MOD(i,nproma)+1, i/nproma+1) = p_pe_work * local_size + i
      DO j = 1, nlev
        in_array_r_2d(MOD(i,nproma)+1, j, i/nproma+1) = (j - 1) * global_size + &
          &                                             p_pe_work * local_size + i
        in_array_i_2d(MOD(i,nproma)+1, j, i/nproma+1) = (j - 1) * global_size + &
          &                                             p_pe_work * local_size + i
      END DO
    END DO
    IF (p_pe_work == 0) THEN
      ALLOCATE(out_array_r_1d(global_size), &
        &      out_array_r_2d(global_size, nlev), &
        &      out_array_i_1d(global_size), &
        &      out_array_i_2d(global_size, nlev), &
        &      ref_out_array_r_1d(global_size), &
        &      ref_out_array_r_2d(global_size, nlev), &
        &      ref_out_array_i_1d(global_size), &
        &      ref_out_array_i_2d(global_size, nlev))
      DO i = 0, global_size - 1, 2
        ref_out_array_r_1d(i+1) = i
        ref_out_array_i_1d(i+1) = i
        DO j = 1, nlev
          ref_out_array_r_2d(i+1, j) = (j - 1) * global_size + i
          ref_out_array_i_2d(i+1, j) = (j - 1) * global_size + i
        END DO
      END DO
      DO i = 1, global_size - 1, 2
        ref_out_array_r_1d(i+1) = fill_value
        ref_out_array_i_1d(i+1) = fill_value
        DO j = 1, nlev
          ref_out_array_r_2d(i+1, j) = fill_value
          ref_out_array_i_2d(i+1, j) = fill_value
        END DO
      END DO
    ELSE
      ALLOCATE(out_array_r_1d(0), &
        &      out_array_r_2d(0, 0), &
        &      out_array_i_1d(0), &
        &      out_array_i_2d(0, 0), &
        &      ref_out_array_r_1d(0), &
        &      ref_out_array_r_2d(0, 0), &
        &      ref_out_array_i_1d(0), &
        &      ref_out_array_i_2d(0, 0))
    END IF

    ! check gather pattern
    CALL check_exchange_gather(in_array_r_1d, in_array_r_2d, &
      &                        in_array_i_1d, in_array_i_2d, &
      &                        out_array_r_1d, out_array_r_2d, &
      &                        ref_out_array_r_1d, ref_out_array_r_2d, &
      &                        out_array_i_1d, out_array_i_2d, &
      &                        ref_out_array_i_1d, ref_out_array_i_2d, &
      &                        gather_pattern, 3842, fill_value, INT(fill_value))

    ! delete gather pattern and other arrays
    DEALLOCATE(in_array_r_1d, in_array_r_2d, in_array_i_1d, in_array_i_2d)
    DEALLOCATE(out_array_r_1d, out_array_r_2d, out_array_i_1d, out_array_i_2d, &
      &        ref_out_array_r_1d, ref_out_array_r_2d, &
      &        ref_out_array_i_1d, ref_out_array_i_2d)
    CALL delete_comm_gather_pattern(gather_pattern)
    DEALLOCATE(owner_local, glb_index)
    !---------------------------------------------------------------------------
    ! test in which only even numbered global indices are owned by processes
    !---------------------------------------------------------------------------
    ! generate gather pattern
    local_size = 10 * nproma
    global_size = p_n_work * local_size
    ALLOCATE(owner_local(local_size), glb_index(local_size))
    DO i = 2, local_size, 2
      owner_local(i) = p_pe_work
    END DO
    DO i = 1, local_size, 2
      owner_local(i) = -1
    END DO
    DO i = 1, local_size
      glb_index(i) = local_size * p_pe_work + i
    END DO
    disable_consistency_check = .TRUE.

    CALL setup_comm_gather_pattern(global_size, owner_local, glb_index, &
      &                            gather_pattern, disable_consistency_check)

    ! initialise in- and reference out data
    nlev = 5
    fill_value = -1
    ALLOCATE(in_array_r_1d(nproma, local_size / nproma), &
      &      in_array_r_2d(nproma, nlev, local_size / nproma), &
      &      in_array_i_1d(nproma, local_size / nproma), &
      &      in_array_i_2d(nproma, nlev, local_size / nproma))
    DO i = 0, local_size-1
      in_array_r_1d(MOD(i,nproma)+1, i/nproma+1) = p_pe_work * local_size + i
      in_array_i_1d(MOD(i,nproma)+1, i/nproma+1) = p_pe_work * local_size + i
      DO j = 1, nlev
        in_array_r_2d(MOD(i,nproma)+1, j, i/nproma+1) = (j - 1) * global_size + &
          &                                             p_pe_work * local_size + i
        in_array_i_2d(MOD(i,nproma)+1, j, i/nproma+1) = (j - 1) * global_size + &
          &                                             p_pe_work * local_size + i
      END DO
    END DO
    IF (p_pe_work == 0) THEN
      ALLOCATE(out_array_r_1d(global_size), &
        &      out_array_r_2d(global_size, nlev), &
        &      out_array_i_1d(global_size), &
        &      out_array_i_2d(global_size, nlev), &
        &      ref_out_array_r_1d(global_size), &
        &      ref_out_array_r_2d(global_size, nlev), &
        &      ref_out_array_i_1d(global_size), &
        &      ref_out_array_i_2d(global_size, nlev))
      DO i = 1, global_size - 1, 2
        ref_out_array_r_1d(i+1) = i
        ref_out_array_i_1d(i+1) = i
        DO j = 1, nlev
          ref_out_array_r_2d(i+1, j) = (j - 1) * global_size + i
          ref_out_array_i_2d(i+1, j) = (j - 1) * global_size + i
        END DO
      END DO
      DO i = 0, global_size - 1, 2
        ref_out_array_r_1d(i+1) = fill_value
        ref_out_array_i_1d(i+1) = fill_value
        DO j = 1, nlev
          ref_out_array_r_2d(i+1, j) = fill_value
          ref_out_array_i_2d(i+1, j) = fill_value
        END DO
      END DO
    ELSE
      ALLOCATE(out_array_r_1d(0), &
        &      out_array_r_2d(0, 0), &
        &      out_array_i_1d(0), &
        &      out_array_i_2d(0, 0), &
        &      ref_out_array_r_1d(0), &
        &      ref_out_array_r_2d(0, 0), &
        &      ref_out_array_i_1d(0), &
        &      ref_out_array_i_2d(0, 0))
    END IF

    ! check gather pattern
    CALL check_exchange_gather(in_array_r_1d, in_array_r_2d, &
      &                        in_array_i_1d, in_array_i_2d, &
      &                        out_array_r_1d, out_array_r_2d, &
      &                        ref_out_array_r_1d, ref_out_array_r_2d, &
      &                        out_array_i_1d, out_array_i_2d, &
      &                        ref_out_array_i_1d, ref_out_array_i_2d, &
      &                        gather_pattern, 3932, fill_value, INT(fill_value))

    ! delete gather pattern and other arrays
    DEALLOCATE(in_array_r_1d, in_array_r_2d, in_array_i_1d, in_array_i_2d)
    DEALLOCATE(out_array_r_1d, out_array_r_2d, out_array_i_1d, out_array_i_2d, &
      &        ref_out_array_r_1d, ref_out_array_r_2d, &
      &        ref_out_array_i_1d, ref_out_array_i_2d)
    CALL delete_comm_gather_pattern(gather_pattern)
    DEALLOCATE(owner_local, glb_index)
    !---------------------------------------------------------------------------
    ! simple test in which each process has its own local contiguous part of the
    ! global array plus some overlap with neighbouring processes
    ! only odd numbered global indices are owned by processes
    !---------------------------------------------------------------------------
    IF (p_n_work > 1) THEN
      ! generate gather pattern
      local_size = 12 * nproma
      global_size = p_n_work * 10 * nproma
      ALLOCATE(owner_local(local_size), glb_index(local_size))
      DO i = 1, local_size
        owner_local(i) = p_pe_work
        glb_index(i) =  &
          MOD(global_size + (p_pe_work * 10 - 1) * nproma + i - 1, global_size) + 1
      END DO
      DO i = 1, nproma
        owner_local(i) = MOD(p_n_work + p_pe_work - 1, p_n_work)
        owner_local(11 * nproma + i) = MOD(p_pe_work + 1, p_n_work)
      END DO
      DO i = 2, local_size, 2
        owner_local(i) = -1
      END DO
      disable_consistency_check = .TRUE.
      CALL setup_comm_gather_pattern(global_size, owner_local, glb_index, &
        &                            gather_pattern, disable_consistency_check)

      ! initialise in- and reference out data
      nlev = 5
      fill_value = -1
      ALLOCATE(in_array_r_1d(nproma, local_size / nproma), &
        &      in_array_r_2d(nproma, nlev, local_size / nproma), &
        &      in_array_i_1d(nproma, local_size / nproma), &
        &      in_array_i_2d(nproma, nlev, local_size / nproma))

      DO i = 0, local_size-1
        in_array_r_1d(MOD(i,nproma)+1, i/nproma+1) = &
          MOD(global_size + (p_pe_work * 10 - 1) * nproma + i, global_size)
        in_array_i_1d(MOD(i,nproma)+1, i/nproma+1) = &
          MOD(global_size + (p_pe_work * 10 - 1) * nproma + i, global_size)
        DO j = 1, nlev
          in_array_r_2d(MOD(i,nproma)+1, j, i/nproma+1) = &
            in_array_r_1d(MOD(i,nproma)+1, i/nproma+1) + (j - 1) * global_size
          in_array_i_2d(MOD(i,nproma)+1, j, i/nproma+1) = &
            in_array_i_1d(MOD(i,nproma)+1, i/nproma+1) + (j - 1) * global_size
        END DO
      END DO
      IF (p_pe_work == 0) THEN
        ALLOCATE(out_array_r_1d(global_size), &
          &      out_array_r_2d(global_size, nlev), &
          &      out_array_i_1d(global_size), &
          &      out_array_i_2d(global_size, nlev), &
          &      ref_out_array_r_1d(global_size), &
          &      ref_out_array_r_2d(global_size, nlev), &
          &      ref_out_array_i_1d(global_size), &
          &      ref_out_array_i_2d(global_size, nlev))
        DO i = 0, global_size - 1, 2
          ref_out_array_r_1d(i+1) = i
          ref_out_array_i_1d(i+1) = i
          DO j = 1, nlev
            ref_out_array_r_2d(i+1, j) = (j - 1) * global_size + i
            ref_out_array_i_2d(i+1, j) = (j - 1) * global_size + i
          END DO
        END DO
        DO i = 1, global_size - 1, 2
          ref_out_array_r_1d(i+1) = fill_value
          ref_out_array_i_1d(i+1) = fill_value
          DO j = 1, nlev
            ref_out_array_r_2d(i+1, j) = fill_value
            ref_out_array_i_2d(i+1, j) = fill_value
          END DO
        END DO
      ELSE
        ALLOCATE(out_array_r_1d(0), &
          &      out_array_r_2d(0, 0), &
          &      out_array_i_1d(0), &
          &      out_array_i_2d(0, 0), &
          &      ref_out_array_r_1d(0), &
          &      ref_out_array_r_2d(0, 0), &
          &      ref_out_array_i_1d(0), &
          &      ref_out_array_i_2d(0, 0))
      END IF

      ! check gather pattern
      CALL check_exchange_gather(in_array_r_1d, in_array_r_2d, &
        &                        in_array_i_1d, in_array_i_2d, &
        &                        out_array_r_1d, out_array_r_2d, &
        &                        ref_out_array_r_1d, ref_out_array_r_2d, &
        &                        out_array_i_1d, out_array_i_2d, &
        &                        ref_out_array_i_1d, ref_out_array_i_2d, &
        &                        gather_pattern, 4030, fill_value, INT(fill_value))

      ! initialise reference out data (for no fill_value case)
      DEALLOCATE(out_array_r_1d, out_array_r_2d, out_array_i_1d, out_array_i_2d, &
        &        ref_out_array_r_1d, ref_out_array_r_2d, &
        &        ref_out_array_i_1d, ref_out_array_i_2d)
      IF (p_pe_work == 0) THEN
        ALLOCATE(out_array_r_1d(global_size/2), &
          &      out_array_r_2d(global_size/2, nlev), &
          &      out_array_i_1d(global_size/2), &
          &      out_array_i_2d(global_size/2, nlev), &
          &      ref_out_array_r_1d(global_size/2), &
          &      ref_out_array_r_2d(global_size/2, nlev), &
          &      ref_out_array_i_1d(global_size/2), &
          &      ref_out_array_i_2d(global_size/2, nlev))
        DO i = 1, global_size/2
          ref_out_array_r_1d(i) = 2 * (i - 1)
          ref_out_array_i_1d(i) = 2 * (i - 1)
          DO j = 1, nlev
            ref_out_array_r_2d(i, j) = (j - 1) * global_size + 2 * (i - 1)
            ref_out_array_i_2d(i, j) = (j - 1) * global_size + 2 * (i - 1)
          END DO
        END DO
      ELSE
        ALLOCATE(out_array_r_1d(0), &
          &      out_array_r_2d(0, 0), &
          &      out_array_i_1d(0), &
          &      out_array_i_2d(0, 0), &
          &      ref_out_array_r_1d(0), &
          &      ref_out_array_r_2d(0, 0), &
          &      ref_out_array_i_1d(0), &
          &      ref_out_array_i_2d(0, 0))
      END IF

      ! check gather pattern
      CALL check_exchange_gather(in_array_r_1d, in_array_r_2d, &
        &                        in_array_i_1d, in_array_i_2d, &
        &                        out_array_r_1d, out_array_r_2d, &
        &                        ref_out_array_r_1d, ref_out_array_r_2d, &
        &                        out_array_i_1d, out_array_i_2d, &
        &                        ref_out_array_i_1d, ref_out_array_i_2d, &
        &                        gather_pattern, 4071)

      ! delete gather pattern and other arrays
      DEALLOCATE(in_array_r_1d, in_array_r_2d, in_array_i_1d, in_array_i_2d)
      DEALLOCATE(out_array_r_1d, out_array_r_2d, out_array_i_1d, out_array_i_2d, &
        &        ref_out_array_r_1d, ref_out_array_r_2d, &
        &        ref_out_array_i_1d, ref_out_array_i_2d)
      CALL delete_comm_gather_pattern(gather_pattern)
      DEALLOCATE(owner_local, glb_index)
    END IF ! (p_n_work > 1)


  CONTAINS

    SUBROUTINE check_exchange_gather(in_array_r_1d, in_array_r_2d, &
      &                              in_array_i_1d, in_array_i_2d, &
      &                              out_array_r_1d, out_array_r_2d, &
      &                              ref_out_array_r_1d, ref_out_array_r_2d, &
      &                              out_array_i_1d, out_array_i_2d, &
      &                              ref_out_array_i_1d, ref_out_array_i_2d, &
      &                              gather_pattern, line, fill_value_w, fill_value_i)

      REAL(wp), INTENT(IN) :: in_array_r_1d(:,:), in_array_r_2d(:,:,:)
      INTEGER, INTENT(IN) :: in_array_i_1d(:,:), in_array_i_2d(:,:,:)
      REAL(wp), INTENT(OUT) :: out_array_r_1d(:), out_array_r_2d(:,:)
      REAL(wp), INTENT(IN) :: ref_out_array_r_1d(:), ref_out_array_r_2d(:,:)
      INTEGER, INTENT(OUT) :: out_array_i_1d(:), out_array_i_2d(:,:)
      INTEGER, INTENT(IN) :: ref_out_array_i_1d(:), ref_out_array_i_2d(:,:)
      TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern
      INTEGER, INTENT(in) :: line
      REAL(wp), OPTIONAL, INTENT(IN) :: fill_value_w
      INTEGER, OPTIONAL, INTENT(IN) :: fill_value_i
      INTEGER :: i

      CALL exchange_data(in_array=in_array_r_1d, out_array=out_array_r_1d, &
        &                gather_pattern=gather_pattern, fill_value=fill_value_w)
      IF (p_pe_work == 0) THEN
        IF (ANY(out_array_r_1d /= ref_out_array_r_1d)) THEN
          WRITE(message_text,'(a,i0)') "Wrong gather result r_1d, line=", line
          CALL finish(method_name, message_text)
        END IF
      END IF
      CALL exchange_data(in_array=in_array_r_2d, out_array=out_array_r_2d, &
        &                gather_pattern=gather_pattern, fill_value=fill_value_w)
      IF (p_pe_work == 0) THEN
        IF (ANY(out_array_r_2d /= ref_out_array_r_2d)) THEN
          WRITE(message_text,'(a,i0)') "Wrong gather result r_2d, line=", line
          CALL finish(method_name, message_text)
        END IF
      END IF
      CALL exchange_data(in_array=in_array_i_1d, out_array=out_array_i_1d, &
        &                gather_pattern=gather_pattern, fill_value=fill_value_i)
      IF (p_pe_work == 0) THEN
        IF (ANY(out_array_i_1d /= ref_out_array_i_1d)) THEN
          WRITE(message_text,'(a,i0)') "Wrong gather result i_1d, line=", line
          CALL finish(method_name, message_text)
        END IF
      END IF
      CALL exchange_data(in_array=in_array_i_2d, out_array=out_array_i_2d, &
        &                gather_pattern=gather_pattern, fill_value=fill_value_i)
      IF (p_pe_work == 0) THEN
        IF (ANY(out_array_i_2d /= ref_out_array_i_2d)) THEN
          WRITE(message_text,'(a,i0)') "Wrong gather result i_2d, line=", line
          CALL finish(method_name, message_text)
        END IF
      END IF
    END SUBROUTINE

    SUBROUTINE check_exchange_allgather(in_array_r_1d, &
      &                                 in_array_i_1d, &
      &                                 out_array_r_1d, &
      &                                 ref_out_array_r_1d, &
      &                                 out_array_i_1d, &
      &                                 ref_out_array_i_1d, &
      &                                 allgather_pattern, line, &
      &                                 fill_value_w, fill_value_i)

      REAL(wp), INTENT(IN) :: in_array_r_1d(:,:)
      INTEGER, INTENT(IN) :: in_array_i_1d(:,:)
      REAL(wp), INTENT(OUT) :: out_array_r_1d(:)
      REAL(wp), INTENT(IN) :: ref_out_array_r_1d(:)
      INTEGER, INTENT(OUT) :: out_array_i_1d(:)
      INTEGER, INTENT(IN) :: ref_out_array_i_1d(:)
      TYPE(t_comm_allgather_pattern), INTENT(IN) :: allgather_pattern
      INTEGER, INTENT(in) :: line
      REAL(wp), OPTIONAL, INTENT(IN) :: fill_value_w
      INTEGER, OPTIONAL, INTENT(IN) :: fill_value_i
      INTEGER :: i

      CALL exchange_data(in_array=in_array_r_1d, out_array=out_array_r_1d, &
        &                allgather_pattern=allgather_pattern, &
        &                fill_value=fill_value_w)
      IF (p_pe_work == 0) THEN
        IF (ANY(out_array_r_1d /= ref_out_array_r_1d)) THEN
          WRITE(message_text,'(a,i0)') "Wrong gather result r_1d, line=", line
          CALL finish(method_name, message_text)
        END IF
      END IF
      CALL exchange_data(in_array=in_array_i_1d, out_array=out_array_i_1d, &
        &                allgather_pattern=allgather_pattern, &
        &                fill_value=fill_value_i)
      IF (p_pe_work == 0) THEN
        IF (ANY(out_array_i_1d /= ref_out_array_i_1d)) THEN
          WRITE(message_text,'(a,i0)') "Wrong gather result i_1d, line=", line
          CALL finish(method_name, message_text)
        END IF
      END IF
    END SUBROUTINE

  END SUBROUTINE

  SUBROUTINE bench_exchange_data_mult()

    REAL(wp), POINTER :: var1(:,:,:), var2(:,:,:), var3(:,:,:), &
                             var4(:,:,:), var5(:,:,:), var6(:,:,:), &
                             var7(:,:,:)

    INTEGER :: nlev, nblk, nfields, i
    LOGICAL, SAVE :: first_call = .TRUE.
    INTEGER, SAVE :: timer_1var, timer_2var, timer_3var, timer_4var, &
               timer_5var, timer_6var, timer_7var
    INTEGER, SAVE :: timer_1var_b, timer_2var_b, timer_3var_b, &
                     timer_4var_b, timer_5var_b, timer_6var_b, timer_7var_b
    CLASS(t_comm_pattern), POINTER :: comm_pattern

    IF (first_call) THEN
      first_call = .FALSE.
      timer_1var =   new_timer("bench 1 var no   barrier")
      timer_2var =   new_timer("bench 2 var no   barrier")
      timer_3var =   new_timer("bench 3 var no   barrier")
      timer_4var =   new_timer("bench 4 var no   barrier")
      timer_5var =   new_timer("bench 5 var no   barrier")
      timer_6var =   new_timer("bench 6 var no   barrier")
      timer_7var =   new_timer("bench 7 var no   barrier")

      timer_1var_b = new_timer("bench 1 var with barrier")
      timer_2var_b = new_timer("bench 2 var with barrier")
      timer_3var_b = new_timer("bench 3 var with barrier")
      timer_4var_b = new_timer("bench 4 var with barrier")
      timer_5var_b = new_timer("bench 5 var with barrier")
      timer_6var_b = new_timer("bench 6 var with barrier")
      timer_7var_b = new_timer("bench 7 var with barrier")
    END IF

    nlev = p_patch(n_dom_start)%nlev
    nblk = p_patch(n_dom_start)%nblks_c
    comm_pattern => p_patch(n_dom_start)%comm_pat_c

    ALLOCATE(var1(nproma, nlev, nblk), var2(nproma, nlev, nblk), &
             var3(nproma, nlev, nblk), var4(nproma, nlev, nblk), &
             var5(nproma, nlev, nblk), var6(nproma, nlev, nblk), &
             var7(nproma, nlev, nblk))

    var1 = 0
    var2 = 0
    var3 = 0
    var4 = 0
    var5 = 0
    var6 = 0
    var7 = 0

    ! warm-up
    nfields = 1
    DO i = 1, 16
      CALL exchange_data(comm_pattern, var1)
    END DO
    CALL work_mpi_barrier()
    CALL timer_start(timer_1var)
    ! benchmarking
    DO i = 1, testbed_iterations
      CALL exchange_data(comm_pattern, var1)
    END DO
    CALL work_mpi_barrier()
    CALL timer_stop(timer_1var)

    ! warm-up
    nfields = 2
    DO i = 1, 16
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2)
    END DO
    CALL work_mpi_barrier()
    CALL timer_start(timer_2var)
    ! benchmarking
    DO i = 1, testbed_iterations
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2)
    END DO
    CALL work_mpi_barrier()
    CALL timer_stop(timer_2var)

    ! warm-up
    nfields = 3
    DO i = 1, 16
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2, recv3=var3)
    END DO
    CALL work_mpi_barrier()
    CALL timer_start(timer_3var)
    ! benchmarking
    DO i = 1, testbed_iterations
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2, recv3=var3)
    END DO
    CALL work_mpi_barrier()
    CALL timer_stop(timer_3var)

    ! warm-up
    nfields = 4
    DO i = 1, 16
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2, recv3=var3, recv4=var4)
    END DO
    CALL work_mpi_barrier()
    CALL timer_start(timer_4var)
    ! benchmarking
    DO i = 1, testbed_iterations
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2, recv3=var3, recv4=var4)
    END DO
    CALL work_mpi_barrier()
    CALL timer_stop(timer_4var)

    ! warm-up
    nfields = 5
    DO i = 1, 16
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2, recv3=var3, recv4=var4, recv5=var5)
    END DO
    CALL work_mpi_barrier()
    CALL timer_start(timer_5var)
    ! benchmarking
    DO i = 1, testbed_iterations
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2, recv3=var3, recv4=var4, recv5=var5)
    END DO
    CALL work_mpi_barrier()
    CALL timer_stop(timer_5var)

    ! warm-up
    nfields = 6
    DO i = 1, 16
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2, recv3=var3, recv4=var4, recv5=var5, recv6=var6)
    END DO
    CALL work_mpi_barrier()
    CALL timer_start(timer_6var)
    ! benchmarking
    DO i = 1, testbed_iterations
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2, recv3=var3, recv4=var4, recv5=var5, recv6=var6)
    END DO
    CALL work_mpi_barrier()
    CALL timer_stop(timer_6var)

    ! warm-up
    nfields = 7
    DO i = 1, 16
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2, recv3=var3, recv4=var4, recv5=var5, recv6=var6, recv7=var7)
    END DO
    CALL work_mpi_barrier()
    CALL timer_start(timer_7var)
    ! benchmarking
    DO i = 1, testbed_iterations
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2, recv3=var3, recv4=var4, recv5=var5, recv6=var6, recv7=var7)
    END DO
    CALL work_mpi_barrier()
    CALL timer_stop(timer_7var)

    !---------------------------------------------------------------------------

    ! warm-up
    nfields = 1
    DO i = 1, 16
      CALL exchange_data(comm_pattern, var1)
    END DO
    CALL work_mpi_barrier()
    CALL timer_start(timer_1var_b)
    ! benchmarking
    DO i = 1, testbed_iterations
      CALL work_mpi_barrier()
      CALL exchange_data(comm_pattern, var1)
    END DO
    CALL work_mpi_barrier()
    CALL timer_stop(timer_1var_b)

    ! warm-up
    nfields = 2
    DO i = 1, 16
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2)
    END DO
    CALL work_mpi_barrier()
    CALL timer_start(timer_2var_b)
    ! benchmarking
    DO i = 1, testbed_iterations
      CALL work_mpi_barrier()
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2)
    END DO
    CALL work_mpi_barrier()
    CALL timer_stop(timer_2var_b)

    ! warm-up
    nfields = 3
    DO i = 1, 16
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2, recv3=var3)
    END DO
    CALL work_mpi_barrier()
    CALL timer_start(timer_3var_b)
    ! benchmarking
    DO i = 1, testbed_iterations
      CALL work_mpi_barrier()
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2, recv3=var3)
    END DO
    CALL work_mpi_barrier()
    CALL timer_stop(timer_3var_b)

    ! warm-up
    nfields = 4
    DO i = 1, 16
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2, recv3=var3, recv4=var4)
    END DO
    CALL work_mpi_barrier()
    CALL timer_start(timer_4var_b)
    ! benchmarking
    DO i = 1, testbed_iterations
      CALL work_mpi_barrier()
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2, recv3=var3, recv4=var4)
    END DO
    CALL work_mpi_barrier()
    CALL timer_stop(timer_4var_b)

    ! warm-up
    nfields = 5
    DO i = 1, 16
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2, recv3=var3, recv4=var4, recv5=var5)
    END DO
    CALL work_mpi_barrier()
    CALL timer_start(timer_5var_b)
    ! benchmarking
    DO i = 1, testbed_iterations
      CALL work_mpi_barrier()
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2, recv3=var3, recv4=var4, recv5=var5)
    END DO
    CALL work_mpi_barrier()
    CALL timer_stop(timer_5var_b)

    ! warm-up
    nfields = 6
    DO i = 1, 16
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2, recv3=var3, recv4=var4, recv5=var5, recv6=var6)
    END DO
    CALL work_mpi_barrier()
    CALL timer_start(timer_6var_b)
    ! benchmarking
    DO i = 1, testbed_iterations
      CALL work_mpi_barrier()
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2, recv3=var3, recv4=var4, recv5=var5, recv6=var6)
    END DO
    CALL work_mpi_barrier()
    CALL timer_stop(timer_6var_b)

    ! warm-up
    nfields = 7
    DO i = 1, 16
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2, recv3=var3, recv4=var4, recv5=var5, recv6=var6, recv7=var7)
    END DO
    CALL work_mpi_barrier()
    CALL timer_start(timer_7var_b)
    ! benchmarking
    DO i = 1, testbed_iterations
      CALL work_mpi_barrier()
      CALL exchange_data_mult( &
        p_pat=comm_pattern, nfields=nfields, ndim2tot=nfields*nlev, recv1=var1, &
        recv2=var2, recv3=var3, recv4=var4, recv5=var5, recv6=var6, recv7=var7)
    END DO
    CALL work_mpi_barrier()
    CALL timer_stop(timer_7var_b)

  END SUBROUTINE

END MODULE mo_test_communication

