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

! Contains the main stepping method_name the 3-dim hydrostatic ocean model.

MODULE mo_ocean_testbed_vertical_diffusion
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp, sp
  USE mo_impl_constants,         ONLY: max_char_length, beta_plane_coriolis,full_coriolis, &
    &                               SEA_BOUNDARY, BOUNDARY, SEA, min_dolic
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d, t_subset_range, t_patch_vert
  USE mo_grid_config,            ONLY: n_dom
  USE mo_grid_subset,            ONLY: get_index_range
  USE mo_sync,                   ONLY: sync_patch_array, sync_e, sync_c !, sync_v
  USE mo_ocean_nml,              ONLY: iswm_oce, n_zlev, no_tracer
  USE mo_dynamics_config,        ONLY: nold, nnew
  USE mo_io_config,              ONLY: n_checkpoints
  USE mo_run_config,             ONLY: nsteps, dtime, ltimer, output_mode, test_mode
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_ext_data_types,         ONLY: t_external_data
  !USE mo_io_units,               ONLY: filename_max
  USE mo_ocean_ab_timestepping,  ONLY: solve_free_surface_eq_ab, &
    & calc_normal_velocity_ab,  &
    & calc_vert_velocity,       &
    & update_time_indices
  USE mo_ocean_types,            ONLY: t_hydro_ocean_state, t_hydro_ocean_diag, &
    & t_hydro_ocean_prog
  USE mo_ocean_tracer_transport_types,     ONLY: t_ocean_tracer
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff! , update_diffusion_matrices
  USE mo_ocean_forcing,          ONLY: destruct_ocean_forcing
  USE mo_ocean_physics_types,    ONLY: t_ho_params
  USE mo_name_list_output,       ONLY: write_name_list_output, istime4name_list_output
  USE mo_mpi,                    ONLY: my_process_is_stdio
  USE mo_statistics
  USE mo_sea_ice_nml,            ONLY: i_ice_dyn
  USE mo_util_dbg_prnt,          ONLY: dbg_print
  USE mo_ocean_output
  USE mo_ocean_tracer_diffusion, ONLY: tracer_diffusion_vertical_implicit
  USE mo_ocean_velocity_diffusion,ONLY: velocity_diffusion_vertical_implicit
  USE mo_parallel_config,        ONLY: nproma
  USE mo_math_utility_solvers,   ONLY: apply_triangular_matrix
  USE mo_statistics
  USE mo_timer

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: test_tracer_diffusion_vertical_implicit
  PUBLIC :: test_velocity_diffusion_vert_implicit
  
  CHARACTER(len=12)           :: debug_string = 'testbed     '  ! Output of module for 1 line debug
  
  !-------------------------------------------------------------------------
  INTEGER :: vertical_diffusion_resIterations = 0
CONTAINS

  !-------------------------------------------------------------------------
  SUBROUTINE test_tracer_diffusion_vertical_implicit( patch_3d, ocean_state, physics_parameters,       &
    & operators_coefficients)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE (t_ho_params)                               :: physics_parameters
    TYPE(t_operator_coeff),   INTENT(in)          :: operators_coefficients

    INTEGER :: inner_iter, outer_iter
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_ocean_tracer), POINTER :: ocean_tracer
!     REAL(wp) :: residual(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: start_mean_tracer, mean_tracer, mean_diff
    REAL(wp) :: random_diff(1:nproma,1:n_zlev+1,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    
    ocean_tracer => ocean_state(1)%p_prog(nold(1))%tracer_collection%tracer(1)
    cells_in_domain => patch_3d%p_patch_2D(1)%cells%in_domain

    !---------------------------------------------------------------------
    ! ocean_tracer%concentration(:,:,:) = 0.0_wp
    ! ocean_tracer%concentration(:,17,:) = 1.0_wp
 
    CALL dbg_print('tracer', ocean_tracer%concentration,  debug_string, 1, in_subset=cells_in_domain)

!    CALL fill_tracer_x_height(patch_3d, ocean_state(1))
!    write(0,*) ocean_tracer%concentration_x_height(:, :, 1)
!    sum_trace = SUM(ocean_tracer%concentration_x_height(1, :, 1))
    CALL RANDOM_NUMBER(random_diff)
    random_diff = physics_parameters%A_tracer_v(:,:,:, 1) + 0.0001 * (random_diff - 0.5_wp)
    CALL dbg_print('tracer diffu', random_diff,  debug_string, 1, in_subset=cells_in_domain)
    mean_diff = total_mean(values=random_diff, weights=patch_3d%p_patch_1d(1)%prism_volume, &
      & in_subset=cells_in_domain)
    WRITE(message_text,'(f18.10)') mean_diff
    CALL message("mean_diff=", message_text)

    start_mean_tracer = total_mean(values=ocean_tracer%concentration, weights=patch_3d%p_patch_1d(1)%prism_volume, &
      & in_subset=cells_in_domain)
    WRITE(message_text,'(f18.10)') start_mean_tracer
    CALL message("start_mean_tracer=", message_text)

    DO outer_iter=1,1000
      DO inner_iter=1,1000
        CALL RANDOM_NUMBER(random_diff)
        random_diff = physics_parameters%A_tracer_v(:,:,:, 1) + 0.0001 * (random_diff - 0.5_wp)
        CALL tracer_diffusion_vertical_implicit( patch_3D, ocean_tracer, &
          & random_diff,  ocean_state(1)%p_prog(nnew(1))%h)  
      ENDDO

      WRITE(message_text,'(i6,a)') outer_iter, 'x1000 iter, tracer'
      CALL dbg_print(message_text, ocean_tracer%concentration,  debug_string, 1, in_subset=cells_in_domain)
      CALL dbg_print("tracer diffu", random_diff,  debug_string, 1, in_subset=cells_in_domain)
      
      mean_tracer = total_mean(values=ocean_tracer%concentration, weights=patch_3d%p_patch_1d(1)%prism_volume, &
        & in_subset=cells_in_domain)
      WRITE(message_text,'(f18.10,", ", 2(E22.16,", "))') mean_tracer, mean_tracer-start_mean_tracer, &
        & (mean_tracer-start_mean_tracer)/start_mean_tracer
      CALL message("mean_trace=", message_text)

    ENDDO

  END SUBROUTINE test_tracer_diffusion_vertical_implicit
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE test_velocity_diffusion_vert_implicit( patch_3d, ocean_state, physics_parameters,       &
    & operators_coefficients)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE (t_ho_params)                               :: physics_parameters
    TYPE(t_operator_coeff),   INTENT(in)          :: operators_coefficients

    INTEGER :: inner_iter, outer_iter
    TYPE(t_subset_range), POINTER :: edges_owned, cells_owned
    REAL(wp), POINTER :: vn_inout(:,:,:)
    REAL(wp) :: sum_vn
    REAL(wp) :: ref, power_2
    REAL(wp) :: binary_significant_digits = 5

    INTEGER :: timer_vdiff_old, timer_vdiff

    edges_owned => patch_3d%p_patch_2D(1)%edges%owned
    cells_owned => patch_3d%p_patch_2D(1)%cells%owned
    !---------------------------------------------------------------------
!    CALL dbg_print('edge invthi', patch_3D%p_patch_1D(1)%inv_prism_thick_e,  debug_string, 1, in_subset=edges_owned)
!    CALL dbg_print('edge incent', patch_3D%p_patch_1D(1)%inv_prism_center_dist_e,  debug_string, 1, in_subset=edges_owned)
!    CALL dbg_print('cell invthi', patch_3D%p_patch_1D(1)%inv_prism_thick_c,  debug_string, 1, in_subset=cells_owned)
!    CALL dbg_print('cell invcent', patch_3D%p_patch_1D(1)%inv_prism_center_dist_c,  debug_string, 1, in_subset=cells_owned)
!    STOP
    !---------------------------------------------------------------------
    vn_inout  => ocean_state(1)%p_diag%vn_pred

!    ref = 1._wp / velocity_VerticalDiffusion_background
!    power_2 = (LOG(ref) / LOG(2._wp))
!    write(0,*) ref, LOG(ref), power_2, CEILING(power_2 + binary_significant_digits)
!    power_2 = REAL(2**CEILING(power_2 + binary_significant_digits),wp)
!    write(0,*) "mult=", power_2
!    velocity_VerticalDiffusion_background = REAL(INT(velocity_VerticalDiffusion_background * power_2), wp) / power_2
!    physics_parameters%a_veloc_v(:,:,:) = velocity_VerticalDiffusion_background


    ! physics_parameters%a_veloc_v(:,:,:) = 1.0_wp
    CALL dbg_print('vn diffu', physics_parameters%a_veloc_v,  debug_string, 1, in_subset=edges_owned)

    timer_vdiff_old     = new_timer("vdiff old")
    timer_vdiff         = new_timer("vdiff")

    vn_inout(:,:,:) = 0.0
    vn_inout(:,1,:) = 5.0
    CALL dbg_print('vn', vn_inout,  debug_string, 1, in_subset=edges_owned)
    sum_vn = SUM(vn_inout(1, :, 1) * patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_e(1, :, 1))
    WRITE(message_text,'(f18.10)') sum_vn
    CALL message("sum=", message_text)

    ! physics_parameters%a_veloc_v(:,:,:) = physics_parameters%a_veloc_v(:,:,:) * 1.01_wp
    CALL timer_start(timer_vdiff)
    DO outer_iter=1,1000
      DO inner_iter=1,1000

        CALL velocity_diffusion_vertical_implicit( patch_3d,  &
          & vn_inout,                                  &
          & physics_parameters%a_veloc_v,              &
          & operators_coefficients)

      ENDDO

      WRITE(message_text,'(i6,a)') outer_iter, 'x1000 iter, vn'
      CALL dbg_print(message_text, vn_inout,  debug_string, 1, in_subset=edges_owned)
      sum_vn = SUM(vn_inout(1, :, 1) * patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_e(1, :, 1))
      WRITE(message_text,'(f18.10)') sum_vn
      CALL message("sum=", message_text)

    ENDDO
    CALL timer_stop(timer_vdiff)

    CALL print_timer
    CALL delete_timer(timer_vdiff_old)
    CALL delete_timer(timer_vdiff)

  END SUBROUTINE test_velocity_diffusion_vert_implicit
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE tracer_diffusion_vertical_implicit_r1( patch_3D,               &
                                           & ocean_tracer, h_c, A_v,   &
                                           & operators_coefficients, &
                                           & residual_3D  ) !,  &
                                          ! & diff_column)

    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
    TYPE(t_ocean_tracer), TARGET         :: ocean_tracer
    REAL(wp), INTENT(inout)              :: h_c(:,:)  !surface height, relevant for thickness of first cell, in
    REAL(wp), INTENT(inout)              :: A_v(:,:,:)
    TYPE(t_operator_coeff), INTENT(in),TARGET     :: operators_coefficients
    REAL(wp), INTENT(inout) :: residual_3D(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    ! REAL(wp), INTENT(inout)           :: diff_column(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    !
    !Local variables
    INTEGER :: top_level
    INTEGER :: jc, jk, jb
    INTEGER :: start_index, end_index
    REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)
    REAL(wp) :: fact(1:n_zlev)
    ! REAL(wp) :: inv_prisms_center_distance(1:n_zlev)
    ! REAL(wp) :: inv_prism_thickness(1:n_zlev)
    REAL(wp) :: dt_inv, dt
    REAL(wp), POINTER   :: field_column(:,:,:)
    REAL(wp) :: column_tracer(1:n_zlev), residual(1:n_zlev), old_tracer(1:n_zlev)
    REAL(wp) :: residual_fraction(1:n_zlev)
    REAL(wp) :: prism_thickness(1:n_zlev), inv_prism_thickness(1:n_zlev), inv_prisms_center_distance(1:n_zlev)
    REAL(wp) :: prisms_center_distance(1:n_zlev), unit(1:n_zlev)
    INTEGER  :: bottom_level
    INTEGER  :: resIter
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER         :: patch_2D
    ! CHARACTER(len=max_char_length), PARAMETER :: &
    !        & routine = ('mo_ocean_diffusion:tracer_diffusion_impl')
    !-----------------------------------------------------------------------
    patch_2D         => patch_3D%p_patch_2D(1)
    cells_in_domain => patch_2D%cells%in_domain
    field_column    => ocean_tracer%concentration
    !-----------------------------------------------------------------------
    top_level   = 1
    dt_inv = 1.0_wp/dtime
    dt = dtime
    unit(1:n_zlev) = dt

    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_index, end_index)
      DO jc = start_index, end_index
        bottom_level = patch_3D%p_patch_1D(1)%dolic_c(jc,jb)

        IF (bottom_level > 0 ) THEN

          ! recalculate coefficients
          prism_thickness(1)     = patch_3D%p_patch_1D(1)%del_zlev_m(1) + h_c(jc,jb)
          inv_prism_thickness(1) = 1.0_wp / prism_thickness(1)
          DO jk=2,bottom_level
            prism_thickness(jk)            = patch_3D%p_patch_1D(1)%del_zlev_m(jk)
            ! inv_prism_thickness(jk)        = patch_3D%p_patch_1D(1)%inv_prism_thick_c(jc,jk,jb)
            inv_prism_thickness(jk)        = 1.0_wp / patch_3D%p_patch_1D(1)%del_zlev_m(jk)
            inv_prisms_center_distance(jk) = 2.0_wp / (prism_thickness(jk-1) + prism_thickness(jk))
            prisms_center_distance(jk)     = (prism_thickness(jk-1) + prism_thickness(jk)) * 0.5_wp
            !inv_prisms_center_distance(jk) = 1.0_wp / patch_3D%p_patch_1D(1)%del_zlev_i(jk)
          ENDDO

          !------------------------------------
             ! top level
            a(1) = 0.0_wp
            c(1) = -A_v(jc,2,jb) * inv_prism_thickness(1) * inv_prisms_center_distance(2)
            b(1) = dt_inv - c(1)
            !Fill triangular matrix
            !b is diagonal a is the lower diagonal, c is the upper
            DO jk = 2, bottom_level-1
              a(jk) = - A_v(jc,jk,jb)   * inv_prism_thickness(jk) * inv_prisms_center_distance(jk)
              c(jk) = - A_v(jc,jk+1,jb) * inv_prism_thickness(jk) * inv_prisms_center_distance(jk+1)
              b(jk) = dt_inv - a(jk) - c(jk)
            END DO
            !bottom
            a(bottom_level) = -A_v(jc,bottom_level,jb) * inv_prism_thickness(bottom_level) * &
              & inv_prisms_center_distance(bottom_level)
            c(bottom_level) = 0.0_wp
            b(bottom_level) = dt_inv - a(bottom_level)

            ! without dt_iv
             ! top level
!            a(1) = 0.0_wp
!            c(1) = - A_v(jc,2,jb) * inv_prism_thickness(1) * inv_prisms_center_distance(2)
!            b(1) = 1.0_wp - c(1)
!            !Fill triangular matrix
!            !b is diagonal a is the lower diagonal, c is the upper
!            DO jk = 2, bottom_level-1
!              a(jk) = - A_v(jc,jk,jb) * inv_prism_thickness(jk) * inv_prisms_center_distance(jk)
!              c(jk) = - A_v(jc,jk+1,jb) * inv_prism_thickness(jk) * inv_prisms_center_distance(jk+1)
!              b(jk) = 1.0_wp - a(jk) - c(jk)
!            END DO
!            !bottom
!            a(bottom_level) = -dt * A_v(jc,bottom_level,jb) * inv_prism_thickness(bottom_level) * inv_prisms_center_distance(bottom_level)
!            c(bottom_level) = 0.0_wp
!            b(bottom_level) = 1.0_wp - a(bottom_level)
!            column_tracer(1:bottom_level) = field_column(jc,1:bottom_level,jb)

!
            ! get locally the column tracer / dt
            CALL apply_triangular_matrix(bottom_level,c,b,a,unit(1:bottom_level),residual_fraction)
!            write(0,*) residual
!            stop
 !           column_tracer(1:bottom_level) = column_tracer(1:bottom_level) * residual(1:bottom_level)
!            column_tracer(1:bottom_level) = field_column(jc,1:bottom_level,jb) * dt_inv / residual(1:bottom_level)
            column_tracer(1:bottom_level) = field_column(jc,1:bottom_level,jb) * dt_inv

            column_tracer(1:bottom_level) = field_column(jc,1:bottom_level,jb) * dt_inv

            ! multiply with the diagonal (1, b1, b2, ...)
            DO jk = 2, bottom_level
              column_tracer(jk) = column_tracer(jk) * b(jk-1)
              a(jk) = a(jk) *  b(jk-1)
              b(jk) = b(jk) *  b(jk-1)
              ! c(jk) = c(jk) *  b(jk-1)
              c(jk) = dt_inv * b(jk-1) - a(jk) - b(jk)
            ENDDO
            c(bottom_level) = 0.0_wp
!

          !------------------------------------

          ! solver from lapack
          !
          ! eliminate lower diagonal
          DO jk=top_level, bottom_level-1
            fact(jk+1) = a( jk+1 ) / b( jk )
            b( jk+1 ) = b( jk+1 ) - fact(jk+1) * c( jk )
            column_tracer( jk+1 ) = column_tracer( jk+1 ) - fact(jk+1) * column_tracer( jk )
 !           b( jk+1 ) = b( jk+1 ) - c( jk ) * a( jk+1 ) / b( jk )
 !           column_tracer( jk+1 ) = column_tracer( jk+1 ) - column_tracer( jk ) * a( jk+1 ) / b( jk )
          ENDDO
  !        DO jk=top_level+1, bottom_level
  !          a(jk) = 0.0_wp
  !        ENDDO
          old_tracer(:) = column_tracer(:)

          !     Back solve with the matrix U from the factorization.
          column_tracer( bottom_level ) = column_tracer( bottom_level ) / b( bottom_level )
          DO jk =  bottom_level-1, 1, -1
            column_tracer( jk ) = ( column_tracer( jk ) - c( jk ) * column_tracer( jk+1 ) ) / b( jk )
          ENDDO

          ! check residual
          a(:) = 0.0_wp
          DO resIter=1, vertical_diffusion_resIterations
            CALL apply_triangular_matrix(bottom_level,c,b,a,column_tracer,residual)
            residual(1:bottom_level) = (old_tracer(1:bottom_level) - residual(1:bottom_level)) * residual_fraction(1:bottom_level)
            !     Back solve with the matrix U from the factorization.
            column_tracer( bottom_level ) =  column_tracer( bottom_level ) + residual( bottom_level ) / b( bottom_level )
            DO jk =  bottom_level-1, 1, -1
              column_tracer( jk ) = column_tracer( jk ) + ( residual( jk ) - c( jk ) * residual( jk+1 ) ) / b( jk )
            ENDDO
          ENDDO

          DO jk = 1, bottom_level
            ocean_tracer%concentration(jc,jk,jb) = column_tracer(jk)
          ENDDO

          ! check new residual
          CALL apply_triangular_matrix(bottom_level,c,b,a,column_tracer,residual)
          DO jk = 1, bottom_level
            residual_3D(jc,jk,jb) = old_tracer(jk) - residual(jk)
          ENDDO



        ENDIF ! bottom_level > 0

      END DO ! jc = start_index, end_index
    END DO ! jb = cells_in_domain%start_block, cells_in_domain%end_block

    ! CALL sync_patch_array(SYNC_C, patch_2D, diff_column)

  END SUBROUTINE tracer_diffusion_vertical_implicit_r1
  !------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !!Subroutine implements implicit vertical diffusion for scalar fields.
  !>
  !! sbr identical to sbr above but now with homogeneous boundary conditions
  !!
  !!
  !! The result ocean_tracer%concetration is calculated on domain_cells
  !-------------------------------------------------------------------------
  SUBROUTINE tracer_diffusion_vertical_implicit_r2( patch_3D,               &
                                           & ocean_tracer, A_v,   &
                                           & operators_coefficients ) !,  &
                                          ! & diff_column)

    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
    TYPE(t_ocean_tracer), TARGET         :: ocean_tracer
    REAL(wp), INTENT(inout)              :: A_v(:,:,:)
    TYPE(t_operator_coeff), INTENT(in),TARGET        :: operators_coefficients
    !
    !
    REAL(wp) :: inv_prism_thickness(1:n_zlev), inv_prisms_center_distance(1:n_zlev)
    REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)! , nb(1:n_zlev)
    REAL(wp) :: fact(1:n_zlev)
    REAL(wp) :: column_tracer(1:n_zlev)
    REAL(wp) :: dt_inv, diagonal_product
    REAL(wp), POINTER   :: field_column(:,:,:)
    INTEGER  :: bottom_level
    INTEGER :: jc, jk, jb
    INTEGER :: start_index, end_index
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER         :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    cells_in_domain => patch_2D%cells%in_domain
    field_column    => ocean_tracer%concentration
    !-----------------------------------------------------------------------
    dt_inv = 1.0_wp/dtime

    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_index, end_index)
      DO jc = start_index, end_index
        bottom_level = patch_3D%p_patch_1D(1)%dolic_c(jc,jb)

        IF (bottom_level < 1 ) CYCLE

        DO jk=1,bottom_level
          inv_prism_thickness(jk)        = patch_3D%p_patch_1D(1)%inv_prism_thick_c(jc,jk,jb)
          inv_prisms_center_distance(jk) = patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,jk,jb)
        ENDDO

        !------------------------------------
        ! Fill triangular matrix
        ! b is diagonal, a is the lower diagonal, c is the upper
        !   top level
        a(1) = 0.0_wp
        c(1) = -A_v(jc,2,jb) * inv_prism_thickness(1) * inv_prisms_center_distance(2)
        b(1) = dt_inv - c(1)
        DO jk = 2, bottom_level-1
          a(jk) = - A_v(jc,jk,jb)   * inv_prism_thickness(jk) * inv_prisms_center_distance(jk)
          c(jk) = - A_v(jc,jk+1,jb) * inv_prism_thickness(jk) * inv_prisms_center_distance(jk+1)
          b(jk) = dt_inv - a(jk) - c(jk)
        END DO
        ! bottom
        a(bottom_level) = -A_v(jc,bottom_level,jb) * inv_prism_thickness(bottom_level) * inv_prisms_center_distance(bottom_level)
        b(bottom_level) = dt_inv - a(bottom_level)

        ! precondition: set diagonal equal to diagonal_product
        diagonal_product = PRODUCT(b(1:bottom_level))

        DO jk = 1, bottom_level
           fact(jk) = diagonal_product / b(jk)
           a(jk)  = a(jk)  * fact(jk)
           b(jk)  = diagonal_product
           c(jk)  = dt_inv * fact(jk) - a(jk) - b(jk)
           column_tracer(jk) = field_column(jc,jk,jb) * dt_inv * fact(jk)
        ENDDO
        c(bottom_level) = 0.0_wp

        !------------------------------------
        ! solver from lapack
        !
        ! eliminate upper diagonal
        DO jk=bottom_level-1, 1, -1
          fact(jk+1)  = c( jk ) / b( jk+1 )
          b( jk ) = b( jk ) - fact(jk+1) * a( jk +1 )
          column_tracer( jk ) = column_tracer( jk ) - fact(jk+1) * column_tracer( jk+1 )
        ENDDO

        !     Back solve with the matrix U from the factorization.
        column_tracer( 1 ) = column_tracer( 1 ) / b( 1 )
        DO jk =  2, bottom_level
          column_tracer( jk ) = ( column_tracer( jk ) - a( jk ) * column_tracer( jk-1 ) ) / b( jk )
        ENDDO

        DO jk = 1, bottom_level
          ocean_tracer%concentration(jc,jk,jb) = column_tracer(jk)
        ENDDO

      END DO ! jc = start_index, end_index
    END DO ! jb = cells_in_domain%start_block, cells_in_domain%end_block

  END SUBROUTINE tracer_diffusion_vertical_implicit_r2
  !------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !!Subroutine implements implicit vertical diffusion for scalar fields.
  !>
  !! sbr identical to sbr above but now with homogeneous boundary conditions
  !!
  !!
  !! The result ocean_tracer%concetration is calculated on domain_cells
  SUBROUTINE tracer_diffusion_vertical_implicit_r0( patch_3D,               &
                                           & ocean_tracer, h_c, A_v,   &
                                           & operators_coefficients )

    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
    TYPE(t_ocean_tracer), TARGET         :: ocean_tracer
    REAL(wp), INTENT(inout)              :: h_c(:,:)  !surface height, relevant for thickness of first cell, in
    REAL(wp), INTENT(inout)              :: A_v(:,:,:)
    TYPE(t_operator_coeff), INTENT(in),TARGET     :: operators_coefficients
    !
    !
    REAL(wp) :: inv_prism_thickness(1:n_zlev), inv_prisms_center_distance(1:n_zlev)
    REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)
    REAL(wp) :: fact(1:n_zlev)
    REAL(wp) :: column_tracer(1:n_zlev)
    REAL(wp) :: dt_inv
    REAL(wp), POINTER   :: field_column(:,:,:)

    INTEGER :: top_level, bottom_level
    INTEGER :: jc, jk, jb
    INTEGER :: start_index, end_index
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER         :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    cells_in_domain => patch_2D%cells%in_domain
    field_column    => ocean_tracer%concentration
    !-----------------------------------------------------------------------
    top_level   = 1
    dt_inv = 1.0_wp/dtime

    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_index, end_index)
      DO jc = start_index, end_index
        bottom_level = patch_3D%p_patch_1D(1)%dolic_c(jc,jb)

        IF (bottom_level <= 0 ) CYCLE

        DO jk=1,bottom_level
          inv_prism_thickness(jk)        = patch_3D%p_patch_1D(1)%inv_prism_thick_c(jc,jk,jb)
          inv_prisms_center_distance(jk) = patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,jk,jb)
        ENDDO

        !------------------------------------
        ! Fill triangular matrix
        ! b is diagonal, a is the lower diagonal, c is the upper
        !   top level
        a(1) = 0.0_wp
        c(1) = -A_v(jc,2,jb) * inv_prism_thickness(1) * inv_prisms_center_distance(2)
        b(1) = dt_inv - c(1)
        DO jk = 2, bottom_level-1
          a(jk) = - A_v(jc,jk,jb)   * inv_prism_thickness(jk) * inv_prisms_center_distance(jk)
          c(jk) = - A_v(jc,jk+1,jb) * inv_prism_thickness(jk) * inv_prisms_center_distance(jk+1)
          b(jk) = dt_inv - a(jk) - c(jk)
        END DO
        ! bottom
        a(bottom_level) = -A_v(jc,bottom_level,jb) * inv_prism_thickness(bottom_level) * inv_prisms_center_distance(bottom_level)
        b(bottom_level) = dt_inv - a(bottom_level)

        ! precondition: multiply with the diagonal (1, b1, b2, ...)
        column_tracer(1) = field_column(jc,1,jb) * dt_inv
        DO jk = 2, bottom_level
          a(jk) = a(jk) *  b(jk-1)
          b(jk) = b(jk) *  b(jk-1)
          c(jk) = dt_inv * b(jk-1) - a(jk) - b(jk)

          column_tracer(jk) = field_column(jc,jk,jb) * dt_inv * b(jk-1)

        ENDDO
        c(bottom_level) = 0.0_wp

        !------------------------------------
        ! solver from lapack
        !
        ! eliminate lower diagonal
        DO jk=top_level, bottom_level-1
          fact(jk+1) = a( jk+1 ) / b( jk )
          b( jk+1 ) = b( jk+1 ) - fact(jk+1) * c( jk )
          column_tracer( jk+1 ) = column_tracer( jk+1 ) - fact(jk+1) * column_tracer( jk )
        ENDDO

        !     Back solve with the matrix U from the factorization.
        column_tracer( bottom_level ) = column_tracer( bottom_level ) / b( bottom_level )
        DO jk =  bottom_level-1, 1, -1
          column_tracer( jk ) = ( column_tracer( jk ) - c( jk ) * column_tracer( jk+1 ) ) / b( jk )
        ENDDO

        DO jk = 1, bottom_level
          ocean_tracer%concentration(jc,jk,jb) = column_tracer(jk)
        ENDDO

      END DO ! jc = start_index, end_index
    END DO ! jb = cells_in_domain%start_block, cells_in_domain%end_block

  END SUBROUTINE tracer_diffusion_vertical_implicit_r0
  !------------------------------------------------------------------------

  
  !-------------------------------------------------------------------------
  !!Subroutine tests for vertical diffusion accuracy 
  !>
  !!
  !! The result ocean_tracer%concetration is calculated on domain_cells
  !-------------------------------------------------------------------------
  SUBROUTINE tracer_diffusion_vertical_implicit_r3( patch_3D,               &
                                           & ocean_tracer, A_v)!,   &
                                          ! & diff_column)

    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
    TYPE(t_ocean_tracer), TARGET         :: ocean_tracer
    REAL(wp), INTENT(inout)              :: A_v(:,:,:)
    !
    !
    REAL(wp) :: inv_prism_thickness(1:n_zlev), inv_prisms_center_distance(1:n_zlev)
    REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)! , nb(1:n_zlev)
    REAL(wp) :: fact(1:n_zlev)
    REAL(wp) :: column_tracer(1:n_zlev)
    REAL(wp) :: dt_inv, diagonal_product
    REAL(wp) :: dt_mod, mu ! Testing dt and viscosity
    REAL(wp) :: H, eta                                                   ! Stretching parameters 
    REAL(wp) :: inv_stretch_z_i(1:n_zlev), inv_stretch_z_m(1:n_zlev)     ! Stretching parameters 
    REAL(wp) :: inv_cart_thick(1:n_zlev), inv_cart_center_dist(1:n_zlev) ! Stretching parameters 
    REAL(wp), POINTER   :: field_column(:,:,:)
    INTEGER  :: bottom_level
    INTEGER :: jc, jk, jb
    INTEGER :: start_index, end_index
    INTEGER :: test_diff 
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER         :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    cells_in_domain => patch_2D%cells%in_domain
    field_column    => ocean_tracer%concentration
    !-----------------------------------------------------------------------
    dt_inv = 1.0_wp/dtime

    test_diff =  4 

    field_column(:, :, :) = 0.; 
    field_column(:,11, :) = 1.; 
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_index, end_index)
      DO jc = start_index, end_index
!    DO jb = 1, 1 
!      DO jc = 1, 1 
        bottom_level = patch_3D%p_patch_1D(1)%dolic_c(jc,jb)

        ! FIXME: Are the following assumptions true
        H   = patch_3D%column_thick_e(jc, jb)
        eta = patch_3D%p_patch_1D(1)%depth_CellInterface(jc, bottom_level + 1, jb) - H

        IF (bottom_level < 1 ) CYCLE

        DO jk=1,bottom_level
            inv_stretch_z_i(jk) = H/( H + eta ) 
            inv_stretch_z_m(jk) = H/( H + eta )   ! Identical for z* co-ords
         
            ! FIXME: Are the following assumptions true
            inv_cart_thick(jk)       = 1.0_wp/patch_3d%p_patch_1d(1)%del_zlev_m(jk)
            inv_cart_center_dist(jk) = 1.0_wp/patch_3d%p_patch_1d(1)%del_zlev_i(jk)

        END DO

        DO jk=1,bottom_level
          inv_prism_thickness(jk)        = patch_3D%p_patch_1D(1)%inv_prism_thick_c(jc,jk,jb)
          inv_prisms_center_distance(jk) = patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,jk,jb)
        ENDDO

        SELECT CASE (test_diff)  !  

        CASE (1)
        !------------------------------------
        ! Test d2q/dz2 = f
        ! If we just remove the dt_inv terms the resulting matrix is singular
        ! because of the boundary terms
        ! Therefore the boundary terms are modified to retain coercivity 
        !------------------------------------
            !------------------------------------
            ! Fill triangular matrix
            ! b is diagonal, a is the lower diagonal, c is the upper
            !   top level
            a(1) = 0.0_wp
            c(1) = -A_v(jc,2,jb) * inv_prism_thickness(1) * inv_prisms_center_distance(2)
            b(1) = - 2.*c(1)
            DO jk = 2, bottom_level-1
              a(jk) = - A_v(jc,jk,jb)   * inv_prism_thickness(jk) * inv_prisms_center_distance(jk)
              c(jk) = - A_v(jc,jk+1,jb) * inv_prism_thickness(jk) * inv_prisms_center_distance(jk+1)
              b(jk) = - a(jk) - c(jk)
            END DO
            ! bottom
            a(bottom_level) = -A_v(jc,bottom_level,jb) * &
                & inv_prism_thickness(bottom_level) * inv_prisms_center_distance(bottom_level)
            b(bottom_level) = - 2.*a(bottom_level)
            c(bottom_level) = 0.0_wp
    
            DO jk = 1, bottom_level
               column_tracer(jk) = field_column(jc,jk,jb)* A_v(jc,bottom_level,jb) 
            ENDDO

        CASE (2)
        !------------------------------------
        ! Test d2q/dz2 = f/dt
        !------------------------------------
            dt_mod = 1.0_wp*dt_inv 
            mu     = 0.01

            !------------------------------------
            ! Fill triangular matrix
            ! b is diagonal, a is the lower diagonal, c is the upper
            !   top level
            a(1) = 0.0_wp
            c(1) = -A_v(jc,2,jb) * inv_prism_thickness(1) * inv_prisms_center_distance(2)
            b(1) = dt_mod - c(1)
            DO jk = 2, bottom_level-1
              a(jk) = - A_v(jc,jk,jb)   * inv_prism_thickness(jk) * inv_prisms_center_distance(jk)
              c(jk) = - A_v(jc,jk+1,jb) * inv_prism_thickness(jk) * inv_prisms_center_distance(jk+1) 

              b(jk) = dt_mod - a(jk) - c(jk)
            END DO
            ! bottom
            a(bottom_level) = -A_v(jc,bottom_level,jb) * &
                & inv_prism_thickness(bottom_level) * inv_prisms_center_distance(bottom_level)
            b(bottom_level) = dt_mod - a(bottom_level)
            c(bottom_level) = 0.0_wp
    
            DO jk = 1, bottom_level
               column_tracer(jk) = field_column(jc,jk,jb) * dt_mod 
            ENDDO

        CASE (3)
        !------------------------------------
        ! Test d2q/dz2 = f
        ! With stretched grid 
        !------------------------------------
            !------------------------------------
            ! Fill triangular matrix
            ! b is diagonal, a is the lower diagonal, c is the upper
            !   top level
            a(1) = 0.0_wp
            c(1) = -A_v(jc,2,jb) * inv_stretch_z_i(1) * inv_cart_thick(1) &
                & * inv_stretch_z_m(2) * inv_cart_center_dist(2)
            b(1) = - 2.*c(1)
            DO jk = 2, bottom_level-1
              a(jk) = - A_v(jc,jk,jb)    * inv_stretch_z_i(jk) * inv_cart_thick(jk) &
                 & * inv_stretch_z_m(jk) * inv_cart_center_dist(jk)
              c(jk) = - A_v(jc,jk+1,jb)  * inv_stretch_z_i(jk) * inv_cart_thick(jk) &
                 & * inv_stretch_z_m(jk+ 1) * inv_cart_center_dist(jk+1)
              b(jk) = - a(jk) - c(jk)
            END DO
            ! bottom
            a(bottom_level) = -A_v(jc,bottom_level,jb) *inv_stretch_z_i(bottom_level) &
                & * inv_cart_thick(bottom_level) * inv_stretch_z_m(bottom_level ) &
                & * inv_cart_center_dist(bottom_level)
            b(bottom_level) = - 2.*a(bottom_level)
            c(bottom_level) = 0.0_wp
    
            DO jk = 1, bottom_level
               column_tracer(jk) = field_column(jc,jk,jb)* A_v(jc,bottom_level,jb)  
            ENDDO

        CASE (4)
        !------------------------------------
        ! Test d2q/dz2 = f/dt
        ! With stretched grid 
        !------------------------------------
            dt_mod = 1.0_wp*dt_inv 
            mu     = 0.001

            !------------------------------------
            ! Fill triangular matrix
            ! b is diagonal, a is the lower diagonal, c is the upper
            !   top level
            a(1) = 0.0_wp
            c(1) = -A_v(jc,2,jb) * inv_stretch_z_i(1) * inv_cart_thick(1) &
                & * inv_stretch_z_m(2) * inv_cart_center_dist(2)
            b(1) = dt_mod - c(1)
            DO jk = 2, bottom_level-1
              a(jk) = - A_v(jc,jk,jb)    * inv_stretch_z_i(jk) * inv_cart_thick(jk) &
                 & * inv_stretch_z_m(jk) * inv_cart_center_dist(jk)
              c(jk) = - A_v(jc,jk+1,jb)  * inv_stretch_z_i(jk) * inv_cart_thick(jk) &
                 & * inv_stretch_z_m(jk+ 1) * inv_cart_center_dist(jk+1)

              b(jk) = dt_mod - a(jk) - c(jk)
            END DO
            ! bottom

            a(bottom_level) = -A_v(jc,bottom_level,jb) *inv_stretch_z_i(bottom_level) &
                & * inv_cart_thick(bottom_level) * inv_stretch_z_m(bottom_level ) &
                & * inv_cart_center_dist(bottom_level)

            b(bottom_level) = dt_mod - a(bottom_level)
            c(bottom_level) = 0.0_wp
    
            DO jk = 1, bottom_level
               column_tracer(jk) = field_column(jc,jk,jb) * dt_mod 
            ENDDO



        END SELECT 

        !------------------------------------
        ! The below algorithm is generic and not format specific
        ! 

        !------------------------------------
        ! solver from lapack
        !
        ! eliminate upper diagonal
        DO jk=bottom_level-1, 1, -1
          fact(jk+1)  = c( jk ) / b( jk+1 )
          b( jk ) = b( jk ) - fact(jk+1) * a( jk +1 )
          column_tracer( jk ) = column_tracer( jk ) - fact(jk+1) * column_tracer( jk+1 )
        ENDDO

        !     Back solve with the matrix U from the factorization.
        column_tracer( 1 ) = column_tracer( 1 ) / b( 1 )
        DO jk =  2, bottom_level
          column_tracer( jk ) = ( column_tracer( jk ) - a( jk ) * column_tracer( jk-1 ) ) / b( jk )
        ENDDO

        DO jk = 1, bottom_level
          ocean_tracer%concentration(jc,jk,jb) = column_tracer(jk)
        ENDDO

      END DO ! jc = start_index, end_index
    END DO ! jb = cells_in_domain%start_block, cells_in_domain%end_block

  END SUBROUTINE tracer_diffusion_vertical_implicit_r3
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  SUBROUTINE velocity_diffusion_vertical_implicit_r1( patch_3D,               &
                                           & velocity, A_v,   &
                                           & operators_coefficients ) !,  &
                                          ! & diff_column)

    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
    REAL(wp), INTENT(inout)              :: velocity(:,:,:)   ! on edges
    REAL(wp), INTENT(inout)              :: A_v(:,:,:)
    TYPE(t_operator_coeff), INTENT(in),TARGET        :: operators_coefficients
    !
    REAL(wp) :: dt_inv
    REAL(wp) :: inv_prism_thickness(1:n_zlev), inv_prisms_center_distance(1:n_zlev)
    REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev), diagonal_product
    REAL(wp) :: column_velocity(1:n_zlev)
    REAL(wp) :: fact(1:n_zlev)

    INTEGER :: bottom_level
    INTEGER :: edge_index, jk, edge_block
    INTEGER :: start_index, end_index
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_patch), POINTER         :: patch_2D

    !-----------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    all_edges       => patch_2D%edges%all
    dt_inv = 1.0_wp/dtime
    !-----------------------------------------------------------------------

    DO edge_block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, edge_block, start_index, end_index)
      DO edge_index = start_index, end_index
        bottom_level = patch_3D%p_patch_1D(1)%dolic_e(edge_index,edge_block)

        IF (bottom_level < 1 ) CYCLE

        ! Note : the inv_prism_thick_e, inv_prism_center_dist_e should be updated in calculate_thickness
        DO jk=1, bottom_level
          inv_prism_thickness(jk)        = patch_3D%p_patch_1D(1)%inv_prism_thick_e(edge_index,jk,edge_block)
          inv_prisms_center_distance(jk) = patch_3d%p_patch_1D(1)%inv_prism_center_dist_e(edge_index,jk,edge_block)
        ENDDO

        !------------------------------------
        ! Fill triangular matrix
        ! b is diagonal, a is the lower diagonal, c is the upper
        !   top level
        a(1) = 0.0_wp
        c(1) = -A_v(edge_index,2,edge_block) * inv_prism_thickness(1) * inv_prisms_center_distance(2)
        b(1) = dt_inv - c(1)
        DO jk = 2, bottom_level-1
          a(jk) = - A_v(edge_index,jk,edge_block)   * inv_prism_thickness(jk) * inv_prisms_center_distance(jk)
          c(jk) = - A_v(edge_index,jk+1,edge_block) * inv_prism_thickness(jk) * inv_prisms_center_distance(jk+1)
          b(jk) = dt_inv - a(jk) - c(jk)
        END DO
        ! bottom
        a(bottom_level) = -A_v(edge_index,bottom_level,edge_block) * inv_prism_thickness(bottom_level) * &
          & inv_prisms_center_distance(bottom_level)
        b(bottom_level) = dt_inv - a(bottom_level)

        ! precondition: set diagonal equal to diagonal_product
        diagonal_product = PRODUCT(b(1:bottom_level))

        DO jk = 1, bottom_level
           fact(jk) = diagonal_product / b(jk)
           a(jk)  = a(jk)  * fact(jk)
           b(jk)  = diagonal_product
           c(jk)  = dt_inv * fact(jk) - a(jk) - b(jk)
           column_velocity(jk) = velocity(edge_index,jk,edge_block) * dt_inv * fact(jk)
        ENDDO
        c(bottom_level) = 0.0_wp

        !------------------------------------
        ! solver from lapack
        !
        ! eliminate upper diagonal
        DO jk = bottom_level-1, 1, -1
          fact(jk+1)  = c( jk ) / b( jk+1 )
          b( jk ) = b( jk ) - fact(jk+1) * a( jk +1 )
          column_velocity( jk ) = column_velocity( jk ) - fact(jk+1) * column_velocity( jk+1 )
        ENDDO

        !     Back solve with the matrix U from the factorization.
        column_velocity( 1 ) = column_velocity( 1 ) / b( 1 )
        DO jk = 2, bottom_level
          column_velocity( jk ) = ( column_velocity( jk ) - a( jk ) * column_velocity( jk-1 ) ) / b( jk )
        ENDDO

        DO jk = 1, bottom_level
          velocity(edge_index,jk,edge_block) = column_velocity(jk)
        ENDDO

      END DO ! edge_index = start_index, end_index
    END DO ! edge_block = cells_in_domain%start_block, cells_in_domain%end_block

  END SUBROUTINE velocity_diffusion_vertical_implicit_r1
  !------------------------------------------------------------------------

   !------------------------------------------------------------------------
  SUBROUTINE velocity_diffusion_vertical_implicit_r2( patch_3D,               &
                                           & velocity, A_v,   &
                                           & operators_coefficients ) !,  &
                                          ! & diff_column)

    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
    REAL(wp), INTENT(inout)              :: velocity(:,:,:)   ! on edges
    REAL(wp), INTENT(inout)              :: A_v(:,:,:)
    TYPE(t_operator_coeff), INTENT(in),TARGET        :: operators_coefficients
    !
    REAL(wp) :: dt_inv
    REAL(wp) :: inv_prism_thickness(1:n_zlev), inv_prisms_center_distance(1:n_zlev)
    REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev), diagonal_product
    REAL(wp) :: column_velocity(1:n_zlev)
    REAL(wp) :: fact(1:n_zlev)

    INTEGER :: bottom_level
    INTEGER :: edge_index, jk, edge_block
    INTEGER :: start_index, end_index
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_patch), POINTER         :: patch_2D

    !-----------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    all_edges       => patch_2D%edges%all
    dt_inv = 1.0_wp/dtime
    !-----------------------------------------------------------------------

    DO edge_block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, edge_block, start_index, end_index)
      DO edge_index = start_index, end_index
        bottom_level = patch_3D%p_patch_1D(1)%dolic_e(edge_index,edge_block)

        IF (bottom_level < 1 ) CYCLE

        ! Note : the inv_prism_thick_e, inv_prism_center_dist_e should be updated in calculate_thickness
        DO jk=1, bottom_level
          inv_prism_thickness(jk)        = patch_3D%p_patch_1D(1)%inv_prism_thick_e(edge_index,jk,edge_block)
          inv_prisms_center_distance(jk) = patch_3d%p_patch_1D(1)%inv_prism_center_dist_e(edge_index,jk,edge_block)
        ENDDO

        !------------------------------------
        ! Fill triangular matrix
        ! b is diagonal, a is the lower diagonal, c is the upper
        !   top level
        a(1) = 0.0_wp
        c(1) = -A_v(edge_index,2,edge_block) * inv_prism_thickness(1) * inv_prisms_center_distance(2)
        b(1) = dt_inv - c(1)
        DO jk = 2, bottom_level-1
          a(jk) = - A_v(edge_index,jk,edge_block)   * inv_prism_thickness(jk) * inv_prisms_center_distance(jk)
          c(jk) = - A_v(edge_index,jk+1,edge_block) * inv_prism_thickness(jk) * inv_prisms_center_distance(jk+1)
          b(jk) = dt_inv - a(jk) - c(jk)
        END DO
        ! bottom
        a(bottom_level) = -A_v(edge_index,bottom_level,edge_block) * inv_prism_thickness(bottom_level) * &
          & inv_prisms_center_distance(bottom_level)
        b(bottom_level) = dt_inv - a(bottom_level)

        ! precondition: set diagonal equal to diagonal_product
        diagonal_product = PRODUCT(b(1:bottom_level))

        DO jk = 1, bottom_level
           fact(jk) = diagonal_product / b(jk)
           a(jk)  = a(jk)  * fact(jk)
           b(jk)  = diagonal_product
           c(jk)  = dt_inv * fact(jk) - a(jk) - b(jk)
           column_velocity(jk) = velocity(edge_index,jk,edge_block) * dt_inv * fact(jk)
        ENDDO
        c(bottom_level) = 0.0_wp

        !------------------------------------
        ! solver from lapack
        !
        ! eliminate upper diagonal
        DO jk = bottom_level-1, 1, -1
          fact(jk+1)  = c( jk ) / b( jk+1 )
          b( jk ) = b( jk ) - fact(jk+1) * a( jk +1 )
          column_velocity( jk ) = column_velocity( jk ) - fact(jk+1) * column_velocity( jk+1 )
        ENDDO

        !     Back solve with the matrix U from the factorization.
        column_velocity( 1 ) = column_velocity( 1 ) / b( 1 )
        DO jk = 2, bottom_level
          column_velocity( jk ) = ( column_velocity( jk ) - a( jk ) * column_velocity( jk-1 ) ) / b( jk )
        ENDDO

        DO jk = 1, bottom_level
          velocity(edge_index,jk,edge_block) = column_velocity(jk)
        ENDDO
        DO jk = bottom_level, n_zlev
          velocity(edge_index,jk,edge_block) = 0.0_wp
        ENDDO


      END DO ! edge_index = start_index, end_index
    END DO ! edge_block = cells_in_domain%start_block, cells_in_domain%end_block

  END SUBROUTINE velocity_diffusion_vertical_implicit_r2
  !------------------------------------------------------------------------


  !------------------------------------------------------------------------
  ! OLD routine
!  SUBROUTINE velocity_diffusion_vertical_implicit_old( p_patch_3D,    &
!                                          & field_column,  &
!                                          & A_v,           &
!                                          & p_op_coeff) !,    &
!    !                                      & diff_column)
!    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
!    REAL(wp), INTENT(inout)           :: field_column(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%nblks_e)
!    !surface height at edges, relevant for thickness of first cell
!    REAL(wp), INTENT(inout)           :: A_v(:,:,:)
!    TYPE(t_operator_coeff), INTENT(inout), TARGET    :: p_op_coeff
!    ! REAL(wp), INTENT(inout)             :: diff_column(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%nblks_e)
!    !
!    !Local variables
!    INTEGER :: slev
!    INTEGER :: je, jk, jb
!    INTEGER :: i_startidx, i_endidx
!    ! REAL(wp) :: a(1:n_zlev), b(1:n_zlev), c(1:n_zlev)
!    !REAL(wp) :: gam(1:n_zlev), bet(1:n_zlev)
!    REAL(wp),POINTER :: a(:), b(:), c(:)
!    REAL(wp) :: z_tmp
!    REAL(wp) :: dt_inv
!    INTEGER  :: z_dolic
!    REAL(wp) :: inv_prisms_center_distance(1:n_zlev)
!    REAL(wp) :: inv_prism_thickness(1:n_zlev)
!    TYPE(t_subset_range), POINTER :: all_edges
!    TYPE(t_patch), POINTER        :: p_patch
!    ! CHARACTER(len=max_char_length), PARAMETER :: &
!    !        & routine = ('mo_ocean_diffusion:tracer_diffusion_impl')
!    !-----------------------------------------------------------------------
!    p_patch   => p_patch_3D%p_patch_2D(1)
!    all_edges => p_patch%edges%all
!    !-----------------------------------------------------------------------
!
!    slev = 1
!    dt_inv=1.0_wp/dtime
!
!    ! diff_column(1:nproma,1:n_zlev,1:p_patch%nblks_e)= 0.0_wp
!
!    field_column(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%nblks_e)&
!    &=field_column(1:nproma,1:n_zlev,1:p_patch_3D%p_patch_2D(1)%nblks_e)*dt_inv
!
!
!    !---------DEBUG DIAGNOSTICS-------------------------------------------
!    !idt_src=5  ! output print level (1-5, fix)
!    !CALL dbg_print('VelDifImplHomIn:field_col' ,field_column, str_module,idt_src, in_subset=p_patch%edges%owned)
!    !---------------------------------------------------------------------
!
!    DO jb = all_edges%start_block, all_edges%end_block
!      CALL get_index_range(all_edges, jb, i_startidx, i_endidx)
!      DO je = i_startidx, i_endidx
!        z_dolic = p_patch_3D%p_patch_1D(1)%dolic_e(je,jb)!!v_base%dolic_e(je,jb)
!
!        ! IF (p_patch_3D%lsm_e(je,1,jb) <= sea_boundary) THEN
!          IF ( z_dolic >= 2 ) THEN
!
!
!             a => p_op_coeff%matrix_vert_diff_e(je,1:z_dolic,jb,1)
!             b => p_op_coeff%matrix_vert_diff_e(je,1:z_dolic,jb,2)
!             c => p_op_coeff%matrix_vert_diff_e(je,1:z_dolic,jb,3)
!
!
!            !Apply the matrix
!            DO jk=slev, z_dolic-1
!              IF(b(jk)/=0.0_wp)THEN
!                a(jk) = a(jk)/b(jk)
!                c(jk) = c(jk)/b(jk)
!                field_column(je,jk,jb)=field_column(je,jk,jb)/b(jk)
!                b(jk)=1.0_wp
!              ENDIF
!            END DO
!
!            DO jk=slev+1, z_dolic-1
!              b(jk)                  = b(jk)-a(jk)*c(jk-1)
!              field_column(je,jk,jb) = field_column(je,jk,jb)&
!                            &-a(jk)*field_column(je,jk-1,jb)
!              c(jk)                  = c(jk)/b(jk)
!              field_column(je,jk,jb) = field_column(je,jk,jb)/b(jk)
!              b(jk)                  = 1.0_wp
!            END DO
!
!            z_tmp = b(z_dolic)-a(z_dolic)*c(z_dolic-1)
!            z_tmp = (field_column(je,z_dolic,jb)-a(z_dolic)*field_column(je,z_dolic-1,jb))/z_tmp
!
!            field_column(je,z_dolic,jb) = z_tmp
!            DO jk = z_dolic-1,1,-1
!              field_column(je,jk,jb) = field_column(je,jk,jb)-c(jk)*field_column(je,jk+1,jb)
!            END DO
!          !  DO jk = 1,z_dolic!-1
!          !    diff_column(je,jk,jb) = field_column(je,jk,jb)
!          !  END DO
!          !ELSEIF ( z_dolic < MIN_DOLIC ) THEN
!          !  diff_column(je,1:z_dolic,jb) = 0.0_wp
!          !  field_column(je,1:z_dolic,jb)= 0.0_wp
!          ENDIF
!        !ELSEIF( v_base%lsm_e(je,1,jb) > sea_boundary ) THEN
!        !  diff_column(je,1:z_dolic,jb) = field_column(je,1:z_dolic,jb)
!         !diff_column(je,:,jb) = 0.0_wp
!         !field_column(je,:,jb)= 0.0_wp
!        ! ENDIF
!      END DO
!    END DO
!  !
!  !  !---------DEBUG DIAGNOSTICS-------------------------------------------
!  !  idt_src=5  ! output print level (1-5, fix)
!  !  CALL dbg_print('VelDifImplHom: field_col'  ,field_column             ,str_module,idt_src, in_subset=p_patch%edges%owned)
!  !  CALL dbg_print('VelDifImplHom: diff_col'   ,diff_column              ,str_module,idt_src, in_subset=p_patch%edges%owned)
!  !  !---------------------------------------------------------------------
!  !
!  END SUBROUTINE velocity_diffusion_vertical_implicit_old
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !> Initialize 3D expansion coefficients.
  !!
  !! NOT USED
  !!
!  SUBROUTINE update_diffusion_matrices(   patch_3D,          &
!                                         & p_phys_param,       &
!                                         & matrix_vert_diff_e, &
!                                         & matrix_vert_diff_c)
!
!     TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
!     TYPE (t_ho_params),       INTENT(IN)    :: p_phys_param
!     REAL(wp), INTENT(INOUT) :: matrix_vert_diff_e(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%nblks_e,1:3)
!     REAL(wp), INTENT(INOUT) :: matrix_vert_diff_c(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%alloc_cell_blocks,1:3)
!    !
!    !
!    ! matrix_vert_diff_c(:,:,:,1) : lower diagonal term
!    ! matrix_vert_diff_c(:,:,:,2) : diagonal term
!    ! matrix_vert_diff_c(:,:,:,3) : upper diagonal term
!    !
!    !Local variables
!    !
!    TYPE(t_subset_range), POINTER :: all_edges
!    TYPE(t_subset_range), POINTER :: all_cells
!    REAL(wp), POINTER :: A_v(:,:,:)
!    REAL(wp) :: inv_zinv_i(1:n_zlev)
!    REAL(wp) :: inv_zinv_m(1:n_zlev)
!    REAL(wp) :: dt_inv,dz_m(1:n_zlev),dz_i(1:n_zlev)
!    INTEGER  :: je,jc,jb,jk,i_no_t
!    INTEGER  :: slev,z_dolic
!    INTEGER  :: i_startidx_e, i_endidx_e,  i_startidx_c, i_endidx_c
!    TYPE(t_patch), POINTER :: patch_2D
!    !---------------------------------------------------------
!    patch_2D     => patch_3D%p_patch_2D(1)
!    all_cells => patch_2D%cells%all
!    all_edges => patch_2D%edges%all
!    !---------------------------------------------------------
!    slev   = 1
!    dt_inv = 1.0_wp/dtime
!
!    !The vertical diffusion matrices for the tracers
!    DO i_no_t = 1,no_tracer
!
!      A_v => p_phys_param%A_tracer_v(:,:,:, i_no_t)
!
!      DO jb = all_cells%start_block, all_cells%end_block
!        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
!        DO jc = i_startidx_c, i_endidx_c
!          z_dolic = patch_3D%p_patch_1D(1)%dolic_c(jc,jb)!v_base%dolic_c(jc,jb)
!
!          IF ( patch_3D%lsm_c(jc,1,jb) <= SEA_BOUNDARY ) THEN
!            IF ( z_dolic >= 2 ) THEN
!
!              inv_zinv_i(:) = patch_3D%p_patch_1D(1)%inv_prism_center_dist_c(jc,:,jb)
!              inv_zinv_m(:) = patch_3D%p_patch_1D(1)%inv_prism_thick_c(jc,:,jb)
!                    dz_i(:) = patch_3D%p_patch_1D(1)%prism_center_dist_c(jc,:,jb)
!                    dz_m(:) = patch_3D%p_patch_1D(1)%prism_thick_c(jc,:,jb)
!
!              !first level
!              matrix_vert_diff_c(jc,slev,jb,1) = 0.0_wp
!              matrix_vert_diff_c(jc,slev,jb,3) = -A_v(jc,slev+1,jb)*inv_zinv_m(slev)*inv_zinv_i(slev+1)
!              matrix_vert_diff_c(jc,slev,jb,2) = dt_inv- matrix_vert_diff_c(jc,slev,jb,3)
!
!              !Fill triangular matrix
!              !b is diagonal a and c are upper and lower band
!              !This corresponds to the 4th indices: "2", "1" and "3"
!              DO jk = slev+1, z_dolic-1
!                matrix_vert_diff_c(jc,jk,jb,1) = -A_v(jc,jk,jb)  *inv_zinv_m(jk) *inv_zinv_i(jk)
!                matrix_vert_diff_c(jc,jk,jb,3) = -A_v(jc,jk+1,jb)*inv_zinv_m(jk) *inv_zinv_i(jk+1)
!                matrix_vert_diff_c(jc,jk,jb,2)  = dt_inv&
!                                               & -matrix_vert_diff_c(jc,jk,jb,1)&
!                                               & -matrix_vert_diff_c(jc,jk,jb,3)
!              END DO
!
!              !bottom
!              matrix_vert_diff_c(jc,z_dolic,jb,1) = -A_v(jc,z_dolic,jb)*inv_zinv_m(z_dolic)*inv_zinv_i(z_dolic)
!              matrix_vert_diff_c(jc,z_dolic,jb,3) = 0.0_wp
!              matrix_vert_diff_c(jc,z_dolic,jb,2) = dt_inv - matrix_vert_diff_c(jc,z_dolic,jb,1)
!            ENDIF
!          ENDIF
!        END DO
!      END DO
!    END DO
!   !---------------------------------------------------------
!
!
!  !Now the velocity matrix
!  A_v => p_phys_param%A_veloc_v
!
!  DO jb = all_edges%start_block, all_edges%end_block
!    CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
!    DO je = i_startidx_e, i_endidx_e
!      z_dolic = patch_3D%p_patch_1D(1)%dolic_e(je,jb)
!
!      IF (  patch_3D%lsm_e(je,1,jb) <= SEA_BOUNDARY ) THEN
!        IF ( z_dolic >= 2 ) THEN
!
!          !inv_zinv_i(:)=1.0_wp/v_base%del_zlev_i(:)
!          !inv_zinv_m(:)=1.0_wp/v_base%del_zlev_m(:)
!          !inv_zinv_i(:) = p_os%p_diag%inv_prism_center_dist_e(je,:,jb)
!          !inv_zinv_m(:) = p_os%p_diag%inv_prism_thick_e(je,:,jb)
!          inv_zinv_i(:) = patch_3D%p_patch_1D(1)%inv_prism_center_dist_e(je,:,jb)
!          inv_zinv_m(:) = patch_3D%p_patch_1D(1)%inv_prism_thick_e(je,:,jb)
!
!
!          !Fill triangular matrix
!          !b is diagonal a and c are upper and lower band
!          DO jk = slev+1, z_dolic-1
!            !a(jk) = -A_v(jc,jk,jb)  *inv_zinv_m(jk) *inv_zinv_i(jk)
!            !c(jk) = -A_v(jc,jk+1,jb)*inv_zinv_m(jk) *inv_zinv_i(jk+1)
!            !b(jk) = dt_inv-a(jk)-c(jk)
!            matrix_vert_diff_e(je,jk,jb,1) = -A_v(je,jk,jb)  *inv_zinv_m(jk) *inv_zinv_i(jk)
!            matrix_vert_diff_e(je,jk,jb,3) = -A_v(je,jk+1,jb)*inv_zinv_m(jk) *inv_zinv_i(jk+1)
!            matrix_vert_diff_e(je,jk,jb,2)  = dt_inv&
!                                            & -matrix_vert_diff_e(je,jk,jb,1)&
!                                            & -matrix_vert_diff_e(je,jk,jb,3)
!
!          END DO
!
!          ! The first row
!          !c(slev) = -A_v(jc,slev+1,jb)*inv_zinv_m(slev)*inv_zinv_i(slev+1)!*dtime
!          !a(slev) = 0.0_wp
!          !b(slev) = dt_inv- c(slev) !- a(slev)
!          matrix_vert_diff_e(je,slev,jb,3) = -A_v(je,slev+1,jb)*inv_zinv_m(slev)*inv_zinv_i(slev+1)
!          matrix_vert_diff_e(je,slev,jb,1) = 0.0_wp
!          matrix_vert_diff_e(je,slev,jb,2) = dt_inv- matrix_vert_diff_e(je,slev,jb,3)
!
!          ! The last row
!          !a(z_dolic) = -A_v(jc,z_dolic,jb)*inv_zinv_m(z_dolic)*inv_zinv_i(z_dolic)!*dtime
!          !c(z_dolic) = 0.0_wp
!          !b(z_dolic) = dt_inv - a(z_dolic)! - c(z_dolic)
!          matrix_vert_diff_e(je,z_dolic,jb,1) = -A_v(je,z_dolic,jb)*inv_zinv_m(z_dolic)*inv_zinv_i(z_dolic)
!          matrix_vert_diff_e(je,z_dolic,jb,3) = 0.0_wp
!          matrix_vert_diff_e(je,z_dolic,jb,2) = dt_inv - matrix_vert_diff_e(je,z_dolic,jb,1)
!        ENDIF
!      ENDIF
!    END DO
!  END DO
!   END SUBROUTINE update_diffusion_matrices
  !-------------------------------------------------------------------------

END MODULE mo_ocean_testbed_vertical_diffusion
