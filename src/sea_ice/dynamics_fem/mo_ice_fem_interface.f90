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

! Contains the interface needed to call AWI FEM sea ice model
! as well as advection and interpolation routines.

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_ice_fem_interface
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_run_config,          ONLY: dtime, ltimer
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_ice_momentum, timer_ice_interp, timer_ice_advection
  USE mo_exception,           ONLY: finish

! USE mo_grid_config,         ONLY: l_limited_area, n_dom   ! for now sea-ice works on global domain-only
  USE mo_parallel_config,     ONLY: nproma
  USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_impl_constants,      ONLY: sea_boundary

  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff
  USE mo_dynamics_config,     ONLY: nold
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_ocean_nml,           ONLY: atm_pressure_included_in_icedyn, ssh_in_icedyn_type , vert_cor_type
  USE mo_ocean_surface_types, ONLY: t_atmos_for_ocean, t_ocean_surface
  USE mo_physical_constants,  ONLY: grav, rho_ref, sfc_press_pascal
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_atmos_fluxes
  USE mo_sea_ice_nml,         ONLY: i_ice_advec
  USE mo_ice_fem_advection,   ONLY: fct_ice_solve, ice_TG_rhs
  USE mo_math_types,          ONLY: t_cartesian_coordinates
  USE mo_math_utilities,      ONLY: gvec2cvec, cvec2gvec
  USE mo_scalar_product,      ONLY: map_cell2edges_3D, map_edges2cell_3D
  USE mo_icon_interpolation_scalar, ONLY: verts2cells_scalar
  USE mo_ice_fem_interpolation, ONLY: map_edges2verts, map_verts2edges,                 &
                                      gvec2cvec_c_2d, cvec2gvec_c_2d,                   &
                                      rotate_cvec_v, gvec2cvec_v_fem, cvec2gvec_v_fem,  &
                                      cells2verts_scalar_seaice
  USE mo_fortran_tools,       ONLY: set_acc_host_or_device
#ifdef _OPENACC
  USE openacc, ONLY: acc_is_present 
#endif

  IMPLICIT NONE

  PUBLIC  :: ice_fem_interface

  PUBLIC  :: ice_fem_init_vel_restart
  PUBLIC  :: ice_fem_update_vel_restart

  PRIVATE :: map_icon2fem_vec
  PRIVATE :: map_fem2icon_vec
  PRIVATE :: map_icon2fem_scalar
  PRIVATE :: map_fem2icon_scalar

  CHARACTER(len=12)           :: str_module    = 'IceFemUtils'  ! Output of module for 1 line debug
  INTEGER                     :: idt_src       = 1         ! Level of detail for 1 line debug

CONTAINS

!------------------------------------------------------------------------
!
!>  Wrapper for the call to the AWI FEM ice model
!!  We first remap the neccesary inputs, then call the momentum solver (EVPdynamics)
!!  and map the resulting velocity onto edges and cell centres.
!!
  SUBROUTINE ice_fem_interface( p_patch_3D, p_ice, p_os, p_as, &
                                atmos_fluxes, p_op_coeff, p_oce_sfc, lacc )

    USE mo_ice_fem_types,     ONLY: sigma11, sigma12, sigma22
    USE mo_ice_fem_evp,       ONLY: EVPdynamics
    USE mo_ice_fem_icon_init, ONLY: c2v_wgt
    USE mo_ice_fem_types,     ONLY: m_ice, m_snow, a_ice, elevation
    USE mo_ice_fem_icon_init, ONLY: rot_mat_3D
    USE mo_ice_fem_types,     ONLY: u_w, v_w, stress_atmice_x, stress_atmice_y
    USE mo_ice_fem_mesh,      ONLY: coord_nod2D
    USE mo_ice_fem_types,     ONLY: u_ice, v_ice
!    USE mo_ice_fem_evp_old,  ONLY: EVPdynamics_old ! non-optimized, original version of the solver

    TYPE(t_patch_3D), TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_sea_ice),          INTENT(INOUT)  :: p_ice
    TYPE(t_hydro_ocean_state),INTENT(IN)     :: p_os
    TYPE (t_atmos_fluxes),    INTENT(IN)     :: atmos_fluxes
    TYPE(t_operator_coeff),   INTENT(IN)     :: p_op_coeff
    TYPE(t_atmos_for_ocean),  INTENT(INOUT)  :: p_as
    TYPE(t_subset_range),     POINTER        :: all_cells
    TYPE(t_ocean_surface),    INTENT(INOUT)  :: p_oce_sfc
    LOGICAL, INTENT(IN), OPTIONAL            :: lacc

    ! Local variables
    TYPE(t_patch), POINTER :: p_patch
    REAL(wp), ALLOCATABLE  :: ssh(:,:) ! sea surface height (input only)         [m]
    REAL(wp), ALLOCATABLE  :: ssh_reduced(:,:) ! reduced sea surface height to take slp coupling into account     [m]
    INTEGER                :: jc, jb, start_index, end_index
    LOGICAL                :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

!--------------------------------------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    all_cells => p_patch_3d%p_patch_2d(1)%cells%all

    !$ACC UPDATE HOST(p_ice%draftave) ASYNC(1) IF(lzacc .AND. acc_is_present(p_ice%draftave))
    !$ACC WAIT(1)
    ALLOCATE(ssh(SIZE(p_ice%draftave(:,:),1),SIZE(p_ice%draftave(:,:),2)))

    !$ACC DATA COPYIN(coord_nod2D, rot_mat_3D, c2v_wgt) &
    !$ACC   COPY(ssh) IF(lzacc)

    !$ACC DATA COPYIN(u_ice, v_ice) &
    !$ACC   COPY(m_ice, m_snow, a_ice, elevation, u_w, v_w, stress_atmice_x, stress_atmice_y) IF(lzacc)
 
    IF (ssh_in_icedyn_type == 1) THEN  ! Fully including ssh
      IF (vert_cor_type == 1) THEN
        !$ACC KERNELS DEFAULT(PRESENT) IF(lzacc)
        ssh(:,:) = p_os%p_prog(nold(1))%eta_c(:,:) + p_ice%draftave(:,:)
        !$ACC END KERNELS
      ELSEIF (vert_cor_type == 0) THEN
        !$ACC KERNELS DEFAULT(PRESENT) IF(lzacc)
        ssh(:,:) = p_os%p_prog(nold(1))%h(:,:)
        !$ACC END KERNELS
      ENDIF
    ELSEIF (ssh_in_icedyn_type == 0)THEN  ! Not including ssh at all
      !$ACC KERNELS DEFAULT(PRESENT) IF(lzacc)
      ssh(:,:) = 0.0_wp
      !$ACC END KERNELS
    ELSEIF (ssh_in_icedyn_type > 1)THEN
      CALL finish('Ice dynamics: ', 'FEM dynamics do not include ssh approximation yet (ssh_in_icedyn=2)!')
    ENDIF

! this is the formulation of MPIOM to let the sea dynamics feels the atm pressure
    IF ( atm_pressure_included_in_icedyn ) THEN
      ALLOCATE(ssh_reduced(SIZE(ssh,1),SIZE(ssh,2)))
      !$ACC ENTER DATA CREATE(ssh_reduced) IF(lzacc)
!ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index, jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_index, end_index)
          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
          DO jc = start_index, end_index
            ssh_reduced(jc,jb) = ssh (jc,jb) &
                 + (p_as%pao(jc,jb)-sfc_press_pascal)/(rho_ref*grav)
          ENDDO
          !$ACC END PARALLEL LOOP
      ENDDO
!ICON_OMP_END_PARALLEL_DO

    ENDIF

    IF (ltimer) CALL timer_start(timer_ice_momentum)

! Initialization of  u_ice, v_ice in case of a restart was moved to mo_hydro_ocean_run

!--------------------------------------------------------------------------------------------------
! Interpolate and/or copy ICON variables to FEM variables
!--------------------------------------------------------------------------------------------------
    IF (ltimer) CALL timer_start(timer_ice_interp)

    ! Map scalars to vertices. Obtain: m_ice, m_snow, a_ice, elevation
    IF (atm_pressure_included_in_icedyn) THEN 
      CALL dbg_print('debug femIWrap: ssh_reduced' , ssh_reduced, str_module, 4, in_subset=p_patch%cells%owned)
      CALL map_icon2fem_scalar(p_patch, p_ice, ssh_reduced, lacc=lzacc)
    ELSE 
      CALL map_icon2fem_scalar(p_patch, p_ice, ssh, lacc=lzacc)
    ENDIF

    ! Map vectors to vertices on the rotated-pole grid. Obtain: stress_atmice_x, stress_atmice_y, u_w, v_w
    CALL map_icon2fem_vec(p_patch_3D, p_os, atmos_fluxes, p_op_coeff, lacc=lzacc)

    IF (ltimer) CALL timer_stop(timer_ice_interp)

!--------------------------------------------------------------------------------------------------
! Call FEM EVP solver
!--------------------------------------------------------------------------------------------------

    !$ACC END DATA

    CALL EVPdynamics(lacc=lzacc)

    !$ACC DATA COPYIN(u_ice, v_ice) &
    !$ACC   COPY(u_w, v_w, stress_atmice_x, stress_atmice_y, elevation) IF(lzacc)

!--------------------------------------------------------------------------------------------------
! FCT advection on FEM grid. Advection on ICON grid is done in ice_slow_interface
!--------------------------------------------------------------------------------------------------

    IF (i_ice_advec == 1) THEN
        IF (ltimer) CALL timer_start(timer_ice_advection)

        call ice_TG_rhs
        call fct_ice_solve

        IF (ltimer) CALL timer_stop(timer_ice_advection)
    ENDIF

!--------------------------------------------------------------------------------------------------
! Interpolate and/or copy FEM variables back to ICON variables
!--------------------------------------------------------------------------------------------------
    IF (ltimer) CALL timer_start(timer_ice_interp)

    ! Rotate and interpolate ice velocities back to ICON grid
    CALL map_fem2icon_vec( p_patch_3D, p_ice, p_op_coeff, lacc=lzacc )
    ! If advection is on FEM grid, interp ice scalars back to ICON grid

    IF (i_ice_advec == 1) THEN
        CALL map_fem2icon_scalar( p_patch, p_ice )
    ENDIF

    IF (ltimer) CALL timer_stop(timer_ice_interp)

! Check to make sure that u_ice, v_ice are written to the restart file was moved to mo_hydro_ocean_run

    IF (ltimer) CALL timer_stop(timer_ice_momentum)

    !$ACC END DATA
    !$ACC END DATA
    DEALLOCATE(ssh)
    IF (atm_pressure_included_in_icedyn) THEN
      !$ACC EXIT DATA DELETE(ssh_reduced) IF(lzacc)
      DEALLOCATE(ssh_reduced)
    END IF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('femIWrap: ice_u' , p_ice%u, str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('femIWrap: ice_v' , p_ice%v, str_module, 4, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE ice_fem_interface

  !-------------------------------------------------------------------------
  !
  !> Initialize u_ice, v_ice in case of a restart
  !!
  SUBROUTINE ice_fem_init_vel_restart(p_patch, p_ice)
    USE mo_ice_fem_types,       ONLY: u_ice, v_ice

    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch
    TYPE(t_sea_ice),       INTENT(IN)    :: p_ice

    ! Local variables
    TYPE(t_subset_range), POINTER :: all_verts
    INTEGER :: i_startidx_v, i_endidx_v
    INTEGER :: jk, jb, jv

  !-------------------------------------------------------------------------

    all_verts => p_patch%verts%all

    jk=0
    DO jb = all_verts%start_block, all_verts%end_block
      CALL get_index_range(all_verts, jb, i_startidx_v, i_endidx_v)
      DO jv = i_startidx_v, i_endidx_v
        jk=jk+1
!            ! Strictly speaking p_ice%u_prog and p_ice%v_prog are only used by the restart files
!            ! now, so this does not need to be done every timestep, only after restart file is
!            ! read.
         u_ice(jk) = p_ice%u_prog(jv,jb)
         v_ice(jk) = p_ice%v_prog(jv,jb)
      END DO
    END DO

  END SUBROUTINE ice_fem_init_vel_restart

  !-------------------------------------------------------------------------
  !
  !> Update p_ice%u_prog with last u_ice, v_ice values before writing a restart file
  !!
  SUBROUTINE ice_fem_update_vel_restart(p_patch, p_ice, lacc)
    USE mo_ice_fem_types,       ONLY: u_ice, v_ice

    TYPE(t_patch), TARGET, INTENT(in)       :: p_patch
    TYPE(t_sea_ice),       INTENT(inout)    :: p_ice
    LOGICAL, INTENT(IN), OPTIONAL           :: lacc

    ! Local variables
    ! Patch and ranges
    TYPE(t_subset_range), POINTER :: all_verts
    ! Indexing
    INTEGER :: i_startidx_v, i_endidx_v
    INTEGER :: i_startidx_v_1, i_endidx_v_1
    INTEGER :: jk, jb, jv
    LOGICAL :: lzacc

    !-------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    all_verts => p_patch%verts%all
    i_startidx_v_1 = 0
    i_endidx_v_1   = 0

    !$ACC DATA COPYIN(all_verts, all_verts%start_block) &
    !$ACC   COPY(u_ice, v_ice) IF(lzacc)

    DO jb = all_verts%start_block, all_verts%end_block
      CALL get_index_range(all_verts, jb, i_startidx_v, i_endidx_v)
      IF (jb > all_verts%start_block) &
          CALL get_index_range(all_verts, jb-1, i_startidx_v_1, i_endidx_v_1)

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
      DO jv = i_startidx_v, i_endidx_v
        jk = jv-i_startidx_v+1 + &
            (jb-all_verts%start_block) * (i_endidx_v_1-i_startidx_v_1+1)
        ! Strictly speaking p_ice%u_prog and p_ice%v_prog are only used by the restart files now,
        ! so this does not need to be done every timestep, only before restart file is written.
        p_ice%u_prog(jv,jb) = u_ice(jk)
        p_ice%v_prog(jv,jb) = v_ice(jk)

      END DO
      !$ACC END PARALLEL LOOP
    END DO

    !$ACC END DATA

  END SUBROUTINE ice_fem_update_vel_restart

  !-------------------------------------------------------------------------
  ! #vla, development note for:
  ! ---------- map_icon2fem_vec and map_fem2icon_vec ----------
  ! Old versoin contained a polar singularity because a cartesian vector was
  ! converted to geographical coords form form and then rotated.
  ! FIX: first rotate the cartesian vector to the rotated pole coordinates
  ! and only then convert to geographic form.
  ! #vla old code was removed during clean-up, 05-2017
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !> 1) Interpolate wind stress from cell centers to vertices
  !> 2) Rotate ocean velocities (available on the dual grid, i.e. vertices)
  !-------------------------------------------------------------------------
  !!
  SUBROUTINE map_icon2fem_vec( p_patch_3D, p_os, atmos_fluxes, p_op_coeff, lacc )

    USE mo_ice_fem_icon_init, ONLY: rot_mat_3D
    USE mo_ice_fem_types,     ONLY: u_w, v_w, stress_atmice_x, stress_atmice_y!, u_ice, v_ice !, a_ice
    USE mo_ice_fem_mesh,           ONLY: coord_nod2D

    TYPE(t_patch_3D), TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_hydro_ocean_state),INTENT(IN)     :: p_os
    TYPE (t_atmos_fluxes),    INTENT(IN)     :: atmos_fluxes
    TYPE(t_operator_coeff),   INTENT(IN)     :: p_op_coeff
    LOGICAL, INTENT(IN), OPTIONAL            :: lacc

    ! Local variables
    TYPE(t_patch), POINTER :: p_patch
    LOGICAL :: lzacc
    INTEGER :: i, j

    ! Temporary variables/buffers
    TYPE(t_cartesian_coordinates) :: p_tau_n_c(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)                      :: tau_n(nproma,p_patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_cartesian_coordinates) :: p_tau_n_dual(nproma,p_patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_cartesian_coordinates) :: p_tau_n_dual_fem(nproma,p_patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_cartesian_coordinates) :: p_vn_dual_fem(nproma,p_patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_cartesian_coordinates) :: p_vn_dual_2D (nproma,p_patch_3D%p_patch_2D(1)%nblks_v)

    CALL set_acc_host_or_device(lzacc, lacc)

!--------------------------------------------------------------------------------------------------
    p_patch => p_patch_3D%p_patch_2D(1)

    !$ACC DATA CREATE(p_tau_n_c, tau_n, p_tau_n_dual, p_tau_n_dual_fem, p_vn_dual_fem, p_vn_dual_2D) IF(lzacc)

    !**************************************************************
    ! (1) Convert lat-lon wind stress to cartesian coordinates
    !**************************************************************
    CALL gvec2cvec_c_2d(p_patch_3D, atmos_fluxes%stress_x, atmos_fluxes%stress_y, p_tau_n_c, lacc=lzacc)

    !**************************************************************
    ! (2) Interpolate 3D wind stress from cell centers to edges
    !**************************************************************
#ifdef NAGFOR
    tau_n = 0.0_wp
#endif
    CALL map_cell2edges_3D(p_patch_3D, p_tau_n_c, tau_n ,p_op_coeff, 1, lacc=lzacc)

    CALL sync_patch_array(SYNC_E, p_patch, tau_n)

    !**************************************************************
    ! (3) Interpolate 3D wind stress from edges to vertices
    !**************************************************************
#ifdef NAGFOR
    p_tau_n_dual(:,:)%x(1) = 0.0_wp
    p_tau_n_dual(:,:)%x(2) = 0.0_wp
    p_tau_n_dual(:,:)%x(3) = 0.0_wp
#endif
    CALL map_edges2verts(p_patch, tau_n, p_op_coeff%edge2vert_coeff_cc, p_tau_n_dual, lacc=lzacc)

    CALL sync_patch_array(SYNC_V, p_patch, p_tau_n_dual%x(1))
    CALL sync_patch_array(SYNC_V, p_patch, p_tau_n_dual%x(2))
    CALL sync_patch_array(SYNC_V, p_patch, p_tau_n_dual%x(3))

    !**************************************************************
    ! (4) Rotate the vectors onto the rotated grid
    !     + convert back to geographic coordinates
    !**************************************************************
    ! atmospheric stress
    CALL rotate_cvec_v(p_patch, p_tau_n_dual, rot_mat_3D, p_tau_n_dual_fem, lacc=lzacc)
    CALL cvec2gvec_v_fem(p_patch, p_tau_n_dual_fem, stress_atmice_x, stress_atmice_y, lacc=lzacc)
    ! ocean velocities
    DO j=1,p_patch_3D%p_patch_2D(1)%nblks_v
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
      DO i=1,nproma
        p_vn_dual_2D(i,j) = p_os%p_diag%p_vn_dual(i,1,j)
      END DO
      !$ACC END PARALLEL LOOP
    END DO
    CALL rotate_cvec_v(p_patch, p_vn_dual_2D, rot_mat_3D, p_vn_dual_fem, lacc=lzacc)
    CALL cvec2gvec_v_fem(p_patch, p_vn_dual_fem, u_w, v_w, lacc=lzacc)

    !$ACC END DATA

  END SUBROUTINE map_icon2fem_vec

  !-------------------------------------------------------------------------
  !> 1) Rotate back ice velocities (available on the FEM grid, i.e. vertices)
  !> 2) Interpolate to cell centers and convert back to geographic coordinates
  !-------------------------------------------------------------------------
  !!
  SUBROUTINE map_fem2icon_vec( p_patch_3D, p_ice, p_op_coeff, lacc )

    USE mo_ice_fem_icon_init, ONLY: rot_mat_3D!, pollon, pollat
    USE mo_ice_fem_types,          ONLY: u_ice, v_ice
    USE mo_ice_fem_mesh,           ONLY: coord_nod2D

    TYPE(t_patch_3D), TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_sea_ice),          INTENT(INOUT)  :: p_ice
    TYPE(t_operator_coeff),   INTENT(IN)     :: p_op_coeff
    LOGICAL, INTENT(IN), OPTIONAL            :: lacc

    ! Local variables
    TYPE(t_patch), POINTER :: p_patch
    LOGICAL :: lzacc
    INTEGER :: i, j

    ! Temporary variables/buffers
    REAL(wp)                      :: vn_e_tmp(nproma, 1, p_patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_cartesian_coordinates) :: p_vn_dual_fem(nproma,p_patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_cartesian_coordinates) :: p_vn_dual(nproma,p_patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_cartesian_coordinates) :: p_vn_c_3D(nproma,1,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates) :: p_vn_c_2D(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)                      :: rot_mat_3D_trans(SIZE(rot_mat_3D,2), SIZE(rot_mat_3D,1))

    CALL set_acc_host_or_device(lzacc, lacc)

!--------------------------------------------------------------------------------------------------
    p_patch => p_patch_3D%p_patch_2D(1)
    rot_mat_3D_trans(:,:) = TRANSPOSE(rot_mat_3D(:,:))

#ifdef NAGFOR
    p_vn_c_3D(:,:,:)%x(1) = 0.0_wp
    p_vn_c_3D(:,:,:)%x(2) = 0.0_wp
    p_vn_c_3D(:,:,:)%x(3) = 0.0_wp
    p_vn_dual_fem(:,:)%x(1) = 0.0_wp
    p_vn_dual_fem(:,:)%x(2) = 0.0_wp
    p_vn_dual_fem(:,:)%x(3) = 0.0_wp
#endif

    !$ACC DATA CREATE(vn_e_tmp, p_vn_dual_fem, p_vn_dual, p_vn_c_3D, p_vn_c_2D) &
    !$ACC   COPYIN(rot_mat_3D_trans) &
    !$ACC   IF(lzacc)

    !**************************************************************
    ! (1) Rotate ice vels to ICON variables + convert to cc
    !**************************************************************
    ! Convert the lat-lon vectors to 3d cartesian
    CALL gvec2cvec_v_fem(p_patch, u_ice, v_ice, p_vn_dual_fem, lacc=lzacc)

    ! Rotate the vectors back onto the ICON grid
    CALL rotate_cvec_v(p_patch, p_vn_dual_fem, rot_mat_3D_trans, p_vn_dual, lacc=lzacc)

    !**************************************************************
    ! (2) Interpolate ice velocities to edges for advection
    !**************************************************************
    CALL map_verts2edges(p_patch, p_vn_dual, p_op_coeff%edge2vert_coeff_cc_t, p_ice%vn_e, lacc=lzacc)

    CALL sync_patch_array(SYNC_E, p_patch, p_ice%vn_e)

    !**************************************************************
    ! (3) ... and cells for drag calculation and output
    !**************************************************************
    DO j = 1,p_patch%nblks_e
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
      DO i = 1,nproma
        vn_e_tmp(i,1,j) = p_ice%vn_e(i,j)
      END DO
      !$ACC END PARALLEL LOOP
    END DO

    CALL map_edges2cell_3D( p_patch_3D, vn_e_tmp, &
      &   p_op_coeff, p_vn_c_3D, 1, 1, lacc=lzacc)

    !**************************************************************
    ! (4) Convert back to geographic coordinates
    !**************************************************************
    DO j = 1,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
      DO i = 1,nproma
        p_vn_c_2D(i,j) = p_vn_c_3D(i,1,j)
      END DO
      !$ACC END PARALLEL LOOP
    END DO
    CALL cvec2gvec_c_2d(p_patch_3D, p_vn_c_2D, p_ice%u, p_ice%v, lacc=lzacc)

    CALL sync_patch_array(SYNC_C, p_patch, p_ice%u)
    CALL sync_patch_array(SYNC_C, p_patch, p_ice%v)

    !$ACC END DATA

  END SUBROUTINE map_fem2icon_vec

  !-------------------------------------------------------------------------
  !> 1) Interpolate ice scalars from vertices to cell centers
  !> 2) Resahpe result to get vars on FEM grid
  !-------------------------------------------------------------------------
  !!
  !pgi$r opt 1
  SUBROUTINE map_icon2fem_scalar(p_patch, p_ice, ssh, lacc)

    USE mo_ice_fem_icon_init, ONLY: c2v_wgt
    USE mo_ice_fem_types,     ONLY: m_ice, m_snow, a_ice, elevation

    TYPE(t_patch), TARGET, INTENT(INOUT) :: p_patch
    TYPE(t_sea_ice),        INTENT(IN)  :: p_ice
    REAL(wp),DIMENSION(nproma,p_patch%alloc_cell_blocks), INTENT(IN) :: ssh
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    ! Temporary variables/buffers
    REAL(wp), DIMENSION (nproma, p_ice%kice, p_patch%nblks_v)   :: buffy_array
    REAL(wp), DIMENSION (nproma*p_patch%nblks_v)                :: buffy
    REAL(wp), DIMENSION (nproma,1,p_patch%alloc_cell_blocks)    :: tmp

    INTEGER :: i, j, k, l
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA CREATE(buffy_array, buffy, tmp) IF(lzacc)

    !$ACC KERNELS DEFAULT(PRESENT) IF(lzacc)
    buffy_array(:,:,:)   = 0.0_wp
    !$ACC END KERNELS


    ! Interpolate tracers to vertices
    !$ACC KERNELS DEFAULT(PRESENT) IF(lzacc)
    tmp(:,:,:) = p_ice%hi(:,:,:) * MAX(TINY(1._wp),p_ice%conc(:,:,:))
    !$ACC END KERNELS

    CALL cells2verts_scalar_seaice( tmp, p_patch, c2v_wgt, buffy_array, lacc=lzacc )

    CALL sync_patch_array(SYNC_V, p_patch, buffy_array )

    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(lzacc)
    DO k=1,p_patch%nblks_v
      DO j=1,p_ice%kice
        DO i=1,nproma
          l = i + (j-1)*nproma + (k-1)*p_ice%kice*nproma
          IF (l<=SIZE(m_ice)) m_ice(l) = buffy_array(i,j,k)
        END DO
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    CALL cells2verts_scalar_seaice( p_ice%conc, p_patch, c2v_wgt, buffy_array, lacc=lzacc )

    CALL sync_patch_array(SYNC_V, p_patch, buffy_array )

    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(lzacc)
    DO k=1,p_patch%nblks_v
      DO j=1,p_ice%kice
        DO i=1,nproma
          l = i + (j-1)*nproma + (k-1)*p_ice%kice*nproma
          IF (l<=SIZE(a_ice)) a_ice(l) = buffy_array(i,j,k)
        END DO
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC KERNELS DEFAULT(PRESENT) IF(lzacc)
    tmp(:,:,:) = p_ice%hs(:,:,:) * MAX(TINY(1._wp),p_ice%conc(:,:,:))
    !$ACC END KERNELS

    CALL cells2verts_scalar_seaice( tmp, p_patch, c2v_wgt, buffy_array, lacc=lzacc )

    CALL sync_patch_array(SYNC_V, p_patch, buffy_array )

    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(lzacc)
    DO k=1,p_patch%nblks_v
      DO j=1,p_ice%kice
        DO i=1,nproma
          l = i + (j-1)*nproma + (k-1)*p_ice%kice*nproma
          IF (l<=SIZE(m_snow)) m_snow(l) = buffy_array(i,j,k)
        END DO
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    ! Interpolate SSH to vertices
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) IF(lzacc)
    DO j=1,SIZE(ssh,2)
      DO i=1,SIZE(ssh,1)
        tmp(i,1,j) = ssh(i,j)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    CALL cells2verts_scalar_seaice( tmp, p_patch, c2v_wgt, buffy_array, lacc=lzacc )
    CALL sync_patch_array(SYNC_V, p_patch, buffy_array )

    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) IF(lzacc)
    DO k=1,p_patch%nblks_v
      DO j=1,p_ice%kice
        DO i=1,nproma
          l = i + (j-1)*nproma + (k-1)*p_ice%kice*nproma
          IF (l<=SIZE(elevation)) elevation(l) = buffy_array(i,j,k)
        END DO
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC END DATA

  END SUBROUTINE map_icon2fem_scalar


  !-------------------------------------------------------------------------
  !> 1) Transform ice vars on FEM grid back to ICON structures
  !> 2) Interpolate from vertices to cell centers
  !-------------------------------------------------------------------------
  !!
  SUBROUTINE map_fem2icon_scalar(p_patch, p_ice)

    USE mo_ice_fem_icon_init, ONLY: v2c_wgt
    USE mo_ice_fem_types,          ONLY: m_ice, m_snow, a_ice

    TYPE(t_patch), TARGET,  INTENT(INOUT)   :: p_patch
    TYPE(t_sea_ice),        INTENT(INOUT)   :: p_ice

    ! Local variables
    TYPE(t_subset_range), POINTER :: all_verts
    INTEGER  :: i_startidx_v, i_endidx_v, jb, jk, jv

    ! Temporary variables/buffers
    REAL(wp), DIMENSION (nproma, p_ice%kice, p_patch%nblks_v) :: m_ice_buff, m_snow_buff, a_ice_buff

    all_verts => p_patch%verts%all

    m_ice_buff(:,:,:)   = 0.0_wp
    m_snow_buff(:,:,:)  = 0.0_wp
    a_ice_buff(:,:,:)   = 0.0_wp

    !**************************************************************
    ! (1) Transform ice vars on FEM grid back to ICON structures
    !**************************************************************
    jk=0 ! position index for vars on FEM grid
    DO jb = all_verts%start_block, all_verts%end_block
      CALL get_index_range(all_verts, jb, i_startidx_v, i_endidx_v)
      DO jv = i_startidx_v,i_endidx_v

        jk = jk+1
        m_ice_buff(jv,1,jb)   = m_ice(jk)
        m_snow_buff(jv,1,jb)  = m_snow(jk)
        a_ice_buff(jv,1,jb)   = a_ice(jk)

      END DO
    END DO
    ! does the same as above
!    m_ice_buff  = RESHAPE(m_ice, SHAPE(m_ice_buff))
!    m_snow_buff = RESHAPE(m_snow, SHAPE(m_snow_buff))
!    a_ice_buff  = RESHAPE(a_ice, SHAPE(a_ice_buff))

    !**************************************************************
    ! (2) Interpolate FEM ice variables to cells
    !**************************************************************
    CALL verts2cells_scalar ( m_ice_buff, p_patch, v2c_wgt, p_ice%vol ) ! multiplied by cell-area below
    CALL verts2cells_scalar ( m_snow_buff, p_patch, v2c_wgt,p_ice%vols) ! multiplied by cell-area below
    CALL verts2cells_scalar ( a_ice_buff, p_patch, v2c_wgt, p_ice%conc )

    CALL sync_patch_array   ( SYNC_C, p_patch, p_ice%vol )
    CALL sync_patch_array   ( SYNC_C, p_patch, p_ice%vols )
    CALL sync_patch_array   ( SYNC_C, p_patch, p_ice%conc )

    !**************************************************************
    ! (3) Calculate ICON ice-variables
    !**************************************************************
!ICON_OMP_WORKSHARE
        WHERE ( p_ice%conc(:,1,:) > 0._wp )
          ! New ice and snow thickness
          p_ice%hi  (:,1,:) = p_ice%vol (:,1,:) / p_ice%conc(:,1,:)
          p_ice%hs  (:,1,:) = p_ice%vols(:,1,:) / p_ice%conc(:,1,:)

          ! multiply vol and vols by cell-area, as by definition
          p_ice%vol (:,1,:) = p_ice%vol (:,1,:) * p_patch%cells%area(:,:)
          p_ice%vols(:,1,:) = p_ice%vols(:,1,:) * p_patch%cells%area(:,:)
        ELSEWHERE
          p_ice%hi  (:,1,:) = 0._wp
          p_ice%hs  (:,1,:) = 0._wp
          p_ice%vol (:,1,:) = 0._wp
          p_ice%vols(:,1,:) = 0._wp
          p_ice%conc(:,1,:) = 0._wp
        ENDWHERE
!ICON_OMP_END_WORKSHARE

  END SUBROUTINE map_fem2icon_scalar

END MODULE mo_ice_fem_interface
