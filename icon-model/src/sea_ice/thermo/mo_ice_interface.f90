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

! Provide interfaces to call sea ice ice model from the ocean and atmosphere.

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_ice_interface
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2007
  !
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: dtime, ltimer
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_exception,           ONLY: finish, message
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_util_dbg_prnt,       ONLY: dbg_print

  USE mo_timer,               ONLY: timer_start, timer_stop, timer_ice_fast,   &
    &                               timers_level, timer_extra40, timer_ice_advection
  USE mtime,                  ONLY: datetime, getDayOfYearFromDateTime

  USE mo_physical_constants,  ONLY: rhoi, ki, Tf, ci
  USE mo_sea_ice_nml,         ONLY: i_ice_therm, i_ice_dyn, i_ice_advec, hci_layer
  USE mo_ocean_nml,           ONLY: atmos_flux_analytical_type, atmos_SWnet_const, atmos_sens_const
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_ocean_surface_types, ONLY: t_ocean_surface, t_atmos_for_ocean
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_atmos_fluxes
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff

  USE mo_ice_winton,          ONLY: set_ice_temp_winton
  USE mo_ice_zerolayer,       ONLY: set_ice_temp_zerolayer, set_ice_temp_zerolayer_analytical
  USE mo_ice_slow_thermo,     ONLY: ice_slow_thermo
  USE mo_ice_parameterizations, ONLY: ice_cut_off, set_ice_albedo
  USE mo_ice_fem_interface,   ONLY: ice_fem_interface
  USE mo_ice_fem_interpolation, ONLY: cvec2gvec_c_2d
  USE mo_math_types,          ONLY: t_cartesian_coordinates
  USE mo_ice_advection,       ONLY: ice_advection_upwind

  USE mo_ice_new_dynamics,   ONLY: ice_new_dynamics, &
    &                              interface_boundary_edge_marker, interface_boundary_cell_marker
  USE mo_ocean_nml,          ONLY: n_zlev

  USE mo_impl_constants,     ONLY: sea_boundary,sea, land, boundary
  USE mo_fortran_tools,      ONLY: set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ice_thermodynamics
  PUBLIC :: ice_dynamics
  PUBLIC :: ice_fast_interface
  PUBLIC :: ice_fast

  CHARACTER(len=12)           :: str_module    = 'IceInterface'  ! Output of module for 1 line debug

CONTAINS

  !-------------------------------------------------------------------------
  !
  !>
  !! ice_dynamics
  !!
  !! This function changes:
  !! p_ice      dynamics fields of sea ice
  !!
  SUBROUTINE ice_dynamics(p_patch_3D, p_ice, p_oce_sfc, atmos_fluxes, p_os, p_as, p_op_coeff, lacc)

    TYPE(t_patch_3D ),TARGET,   INTENT(IN)      :: p_patch_3D
    TYPE(t_sea_ice),            INTENT(INOUT)   :: p_ice
    TYPE(t_ocean_surface),      INTENT(INOUT)   :: p_oce_sfc
    TYPE(t_atmos_fluxes),       INTENT(IN)      :: atmos_fluxes
    TYPE(t_atmos_for_ocean),    INTENT(INOUT)   :: p_as
    TYPE(t_hydro_ocean_state),  INTENT(IN)      :: p_os
    TYPE(t_operator_coeff),     INTENT(IN)      :: p_op_coeff
    LOGICAL, INTENT(IN), OPTIONAL               :: lacc

    ! Local variables
    TYPE(t_patch),  POINTER :: p_patch
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: owned_cells
    TYPE(t_subset_range), POINTER :: all_cells
    LOGICAL :: lzacc

    TYPE(t_cartesian_coordinates) :: cvec_ice_velocity(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)


    REAL(wp) ::  boundary_cell_marker(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) ::  boundary_edge_marker(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) ::  ice_x(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) ::  ice_y(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)
    REAL(wp) ::  ice_z(nproma,p_patch_3d%p_patch_2d(1)%nblks_e)

    REAL(wp) :: nix,niy,niz,tix,tiy,tiz

    INTEGER  :: edge_block_1,edge_block_2,edge_block_3,edge_index_1,edge_index_2,edge_index_3
    INTEGER  :: edge_block_i,start_index,end_index,edge_index_i,cell_block,cell_index
    INTEGER  :: jc, jb, i_startidx_c, i_endidx_c

    CALL set_acc_host_or_device(lzacc, lacc)

    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    all_edges       => p_patch%edges%all
    owned_cells     => p_patch%cells%owned
    all_cells       => p_patch%cells%all

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('bef.icedyn: hi  ',p_ice%hi    ,str_module,2, in_subset=p_patch%cells%owned)
    CALL dbg_print('bef.icedyn: hs  ',p_ice%hs    ,str_module,2, in_subset=p_patch%cells%owned)
    CALL dbg_print('bef.icedyn: conc',p_ice%conc  ,str_module,2, in_subset=p_patch%cells%owned)

    !---------------------------------------------------------------------
    ! (1) -------------------- Sea ice dynamics --------------------------
    !---------------------------------------------------------------------
    IF ( i_ice_dyn >= 1 ) THEN

      IF (timers_level > 1) CALL timer_start(timer_extra40)

      IF ( i_ice_dyn == 1 ) THEN
        ! solve for ice velocities (AWI FEM model wrapper)
        CALL ice_fem_interface ( p_patch_3D, p_ice, p_os, p_as, atmos_fluxes, p_op_coeff, p_oce_sfc, lacc=lzacc )
      ENDIF

      IF ( i_ice_dyn == 2 ) THEN
        CALL ice_new_dynamics( p_patch_3D, p_ice, p_os, p_as, atmos_fluxes, p_op_coeff, p_oce_sfc)
      ENDIF

      CALL dbg_print('bef.iceadv: conc',p_ice%conc  ,str_module,2, in_subset=p_patch%cells%owned)

      !---------------------------- advection---FESOM
      IF (i_ice_advec == 0) THEN
        IF (ltimer) CALL timer_start(timer_ice_advection)

        CALL ice_advection_upwind( p_patch_3D, p_op_coeff, p_ice, lacc=lzacc )

        IF (ltimer) CALL timer_stop(timer_ice_advection)
      ELSEIF (i_ice_advec == 1) THEN
         CALL message('ice_dynamics', 'FESOM ICE' )  ! nothing to do here
        ! advection on FEM grid is done inside ice_fem_interface
      ENDIF

 !    ! fix possible overshoots/undershoots after advection (previously, ice_clean_up_dyn)
      CALL ice_cut_off( p_patch, p_ice, lacc=lzacc )

      IF (i_ice_dyn == 2) THEN
        !$ACC DATA CREATE(cvec_ice_velocity, boundary_cell_marker, boundary_edge_marker, ice_x, ice_y, ice_z) IF(lzacc)

        ! kartesischer Vektor auf Kantenmitte
        !$ACC KERNELS DEFAULT(PRESENT) IF(lzacc)
        cvec_ice_velocity(:,:)%x(1)=0.0_wp
        cvec_ice_velocity(:,:)%x(2)=0.0_wp
        cvec_ice_velocity(:,:)%x(3)=0.0_wp

        boundary_cell_marker(:,:,:)=0.0_wp
        boundary_edge_marker(:,:)=0.0_wp
        ice_x(:,:)=0.0_wp
        ice_y(:,:)=0.0_wp
        ice_z(:,:)=0.0_wp
        !$ACC END KERNELS

        CALL interface_boundary_cell_marker(boundary_cell_marker, p_patch_3D, p_ice)
        CALL interface_boundary_edge_marker(boundary_edge_marker,boundary_cell_marker, p_patch_3D, p_ice)

        DO edge_block_i = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, edge_block_i, start_index, end_index)
          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(nix, tix, niy, tiy, niz, tiz) IF(lzacc)
          DO edge_index_i =  start_index, end_index
            nix=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(1)
            tix=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(1)

            niy=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(2)
            tiy=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(2)

            niz=p_patch%edges%primal_cart_normal(edge_index_i,edge_block_i)%x(3)
            tiz=p_patch%edges%dual_cart_normal(edge_index_i,edge_block_i)%x(3)

            ice_x(edge_index_i,edge_block_i)= boundary_edge_marker(edge_index_i,edge_block_i)*&
                      &(nix*p_ice%vn_e(edge_index_i,edge_block_i)+tix*p_ice%vt_e(edge_index_i,edge_block_i))
            ice_y(edge_index_i,edge_block_i)= boundary_edge_marker(edge_index_i,edge_block_i)*&
                      &(niy*p_ice%vn_e(edge_index_i,edge_block_i)+tiy*p_ice%vt_e(edge_index_i,edge_block_i))
            ice_z(edge_index_i,edge_block_i)= boundary_edge_marker(edge_index_i,edge_block_i)*&
                      &(niz*p_ice%vn_e(edge_index_i,edge_block_i)+tiz*p_ice%vt_e(edge_index_i,edge_block_i))

          ENDDO
          !$ACC END PARALLEL LOOP
        ENDDO

! Mittlung auf Zellmitte

        DO cell_block = owned_cells%start_block, owned_cells%end_block
          CALL get_index_range(owned_cells, cell_block, start_index, end_index)
          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(edge_index_1, edge_block_1) &
          !$ACC   PRIVATE(edge_index_2, edge_block_2, edge_index_3, edge_block_3) IF(lzacc)
          DO cell_index = start_index, end_index

            edge_index_1 = p_patch%cells%edge_idx(cell_index, cell_block, 1)
            edge_block_1 = p_patch%cells%edge_blk(cell_index, cell_block, 1)

            edge_index_2 = p_patch%cells%edge_idx(cell_index, cell_block, 2)
            edge_block_2 = p_patch%cells%edge_blk(cell_index, cell_block, 2)

            edge_index_3 = p_patch%cells%edge_idx(cell_index, cell_block, 3)
            edge_block_3 = p_patch%cells%edge_blk(cell_index, cell_block, 3)

            cvec_ice_velocity(cell_index,cell_block)%x(1)=1.0_wp/3.0_wp*(ice_x(edge_index_1,edge_block_1)+&
                          &ice_x(edge_index_2,edge_block_2)+ice_x(edge_index_3,edge_block_3))
            cvec_ice_velocity(cell_index,cell_block)%x(2)=1.0_wp/3.0_wp*(ice_y(edge_index_1,edge_block_1)+&
                          &ice_y(edge_index_2,edge_block_2)+ice_y(edge_index_3,edge_block_3))
            cvec_ice_velocity(cell_index,cell_block)%x(3)=1.0_wp/3.0_wp*(ice_z(edge_index_1,edge_block_1)+&
                          &ice_z(edge_index_2,edge_block_2)+ice_z(edge_index_3,edge_block_3))
          ENDDO
          !$ACC END PARALLEL LOOP
        ENDDO

        CALL cvec2gvec_c_2d(p_patch_3D, cvec_ice_velocity(:,:), p_ice%u, p_ice%v)

        !$ACC END DATA
      ENDIF

      IF (timers_level > 1) CALL timer_stop(timer_extra40)

    ELSE

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
        DO jc = i_startidx_c, i_endidx_c
          p_ice%u(jc,jb) = 0._wp
          p_ice%v(jc,jb) = 0._wp
        END DO
        !$ACC END PARALLEL LOOP
      END DO

    ENDIF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('aft.icedyn: hi  ',p_ice%hi    ,str_module,2, in_subset=p_patch%cells%owned)
    CALL dbg_print('aft.icedyn: hs  ',p_ice%hs    ,str_module,2, in_subset=p_patch%cells%owned)
    CALL dbg_print('aft.icedyn: conc',p_ice%conc  ,str_module,2, in_subset=p_patch%cells%owned)

  END SUBROUTINE ice_dynamics


  !-------------------------------------------------------------------------
  !
  !>
  !! ice_thermodynamics
  !!
  !! This function changes:
  !! p_ice      slow-thermodynamics fields of sea ice
  !! p_oce_sfc  heat and fresh-water fluxes, passed to the ocean
  !!
  SUBROUTINE ice_thermodynamics(p_patch_3D, p_ice, p_oce_sfc, atmos_fluxes, p_os, p_as, p_op_coeff, lacc)

    TYPE(t_patch_3D ),TARGET,   INTENT(IN)      :: p_patch_3D
    TYPE(t_sea_ice),            INTENT(INOUT)   :: p_ice
    TYPE(t_ocean_surface),      INTENT(INOUT)   :: p_oce_sfc
    TYPE(t_atmos_fluxes),       INTENT(IN)      :: atmos_fluxes
    TYPE(t_atmos_for_ocean),    INTENT(INOUT)   :: p_as
    TYPE(t_hydro_ocean_state),  INTENT(IN)      :: p_os
    TYPE(t_operator_coeff),     INTENT(IN)      :: p_op_coeff
    LOGICAL, INTENT(IN), OPTIONAL               :: lacc

    ! Local variables
    TYPE(t_patch),  POINTER :: p_patch
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)

#ifndef _OPENMP
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('bef.icethm: hi  ',p_ice%hi    ,str_module,2, in_subset=p_patch%cells%owned)
    CALL dbg_print('bef.icethm: hs  ',p_ice%hs    ,str_module,2, in_subset=p_patch%cells%owned)
    CALL dbg_print('bef.icethm: conc',p_ice%conc  ,str_module,2, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------
#endif

    !---------------------------------------------------------------------
    ! (2) --------------- Slow sea ice thermodynamics --------------------
    !---------------------------------------------------------------------
    CALL ice_slow_thermo(p_patch_3D, p_os, atmos_fluxes, p_ice, p_oce_sfc, lacc=lzacc)

#ifndef _OPENMP
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('aft.icethm: hi  ',p_ice%hi    ,str_module,2, in_subset=p_patch%cells%owned)
    CALL dbg_print('aft.icethm: hs  ',p_ice%hs    ,str_module,2, in_subset=p_patch%cells%owned)
    CALL dbg_print('aft.icethm: conc',p_ice%conc  ,str_module,2, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------
#endif

  END SUBROUTINE ice_thermodynamics

  !-------------------------------------------------------------------------
  !
  !>
  !! ice_fast_interface: calls fast sea ice thermodynamics in ocean_surface
  !! Note: in coupled runs atmosphere calls ice_fast directly.
  !!
  !! This function changes:
  !! p_ice        - fast-thermodynamics fields (Qtop, Qbot, Tsurf)
  !! atmos_fluxes - ice albedos (albvisdir, albvisdif, albnirdir, albnirdif)
  !!
  SUBROUTINE ice_fast_interface(p_patch, p_ice, atmos_fluxes, this_datetime, lacc)

    TYPE(t_patch), TARGET,      INTENT(IN)      :: p_patch
    TYPE(t_sea_ice),            INTENT(INOUT)   :: p_ice
    TYPE(t_atmos_fluxes),       INTENT(INOUT)   :: atmos_fluxes
    TYPE(datetime), POINTER,    INTENT(IN)      :: this_datetime
    LOGICAL, INTENT(IN), OPTIONAL               :: lacc
    !
    ! local variables
    INTEGER               :: jb, i_startidx_c, i_endidx_c, i, j
    LOGICAL               :: lzacc
    REAL(wp), DIMENSION(SIZE(atmos_fluxes%lat,1), SIZE(atmos_fluxes%lat,2)) :: nonsolar
    REAL(wp), DIMENSION(SIZE(atmos_fluxes%dlatdT,1), SIZE(atmos_fluxes%dlatdT,2)) :: nonsolardT

    TYPE(t_subset_range), POINTER :: all_cells

    CALL set_acc_host_or_device(lzacc, lacc)

    !-----------------------------------------------------------------------
    all_cells       => p_patch%cells%all
    !---------------------------------------------------------------------

    !$ACC DATA CREATE(nonsolar, nonsolardT) IF(lzacc)

!ICON_OMP_PARALLEL_DO PRIVATE(jb, i_startidx_c, i_endidx_c, nonsolar, nonsolardT) SCHEDULE(dynamic)
        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)

          DO j =1,p_ice%kice
            !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
            DO i=1,nproma
              nonsolar(i,j)   = atmos_fluxes%lat(i,j,jb) + atmos_fluxes%sens(i,j,jb) + atmos_fluxes%LWnet(i,j,jb)
              nonsolardT(i,j) = atmos_fluxes%dlatdT(i,j,jb) + atmos_fluxes%dsensdT(i,j,jb) + atmos_fluxes%dLWdT(i,j,jb)
            END DO
            !$ACC END PARALLEL LOOP
          END DO

          CALL ice_fast(i_startidx_c, i_endidx_c, nproma, p_ice%kice, dtime, &
            &   p_ice% Tsurf(:,:,jb),   &          !  intent(inout)
            &   p_ice% T1   (:,:,jb),   &          !  intent(out)   dummy for zerolayer model
            &   p_ice% T2   (:,:,jb),   &          !  intent(out)   dummy for zerolayer model
            &   p_ice% hi   (:,:,jb),   &          !  intent(in)
            &   p_ice% hs   (:,:,jb),   &          !  intent(in)
            &   p_ice% Qtop (:,:,jb),   &          !  intent(out)
            &   p_ice% Qbot (:,:,jb),   &          !  intent(out)
            &   atmos_fluxes%SWnet  (:,:,jb),   &  !  following: intent(in)
            &   nonsolar, &
            &   nonsolardT, &
            &   p_ice% Tfw  (:,  jb),   &
            &   atmos_fluxes%albvisdir(:,:,jb), &  !  intent(out)
            &   atmos_fluxes%albvisdif(:,:,jb), &  !  intent(out)
            &   atmos_fluxes%albnirdir(:,:,jb), &  !  intent(out)
            &   atmos_fluxes%albnirdif(:,:,jb), &  !  intent(out)
            &   doy=getDayOfYearFromDateTime(this_datetime), &
            &   lacc = lzacc)
        END DO
!ICON_OMP_END_PARALLEL_DO

        ! provide constant heat fluxes for special analytical cases
        ! to-do: move to more appropriate place, should not be done here
        SELECT CASE (atmos_flux_analytical_type)
        CASE (102)
          !$ACC KERNELS DEFAULT(PRESENT) IF(lzacc)
          p_ice%Qtop(:,1,:) = atmos_SWnet_const
          !$ACC END KERNELS

          !$ACC KERNELS DEFAULT(PRESENT) IF(lzacc)
          p_ice%Qbot(:,1,:) = 0.0_wp
          !$ACC END KERNELS
        CASE (103)
          !$ACC KERNELS DEFAULT(PRESENT) IF(lzacc)
          p_ice%Qtop(:,1,:) = 0.0_wp
          !$ACC END KERNELS

          !$ACC KERNELS DEFAULT(PRESENT) IF(lzacc)
          p_ice%Qbot(:,1,:) = atmos_sens_const
          !$ACC END KERNELS
        END SELECT

    !$ACC END DATA

#ifndef _OPENMP
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('ice_fast: hi     ',p_ice%hi       ,str_module,3, in_subset=p_patch%cells%owned)
    CALL dbg_print('ice_fast: hs     ',p_ice%hs       ,str_module,3, in_subset=p_patch%cells%owned)

    CALL dbg_print('aft.fast: Tsurf  ',p_ice%tsurf    ,str_module,3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aft.fast: Qtop   ',p_ice%Qtop     ,str_module,3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aft.fast: Qbot   ',p_ice%Qbot     ,str_module,3, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------
#endif

  END SUBROUTINE ice_fast_interface

  !-------------------------------------------------------------------------
  !>
  !! ! ice_fast: ice thermodynamics at atmospheric time-step.
  !!   Calculates ice/snow surface temp, air-ice fluxes and sets albedos.
  !!
  !! This function changes:
  !! Tsurf and T1, T2 (winton)  - temperature of snow/ice layer(s)
  !! Qtop, Qbot                 - heat flux available for surface/bottom melting
  !! alb{vis/nir}{dir/dif}      - albedos
  !!
  SUBROUTINE ice_fast(i_startidx_c, i_endidx_c, nbdim, kice, pdtime, &
            &   Tsurf,          & ! Surface temperature [degC]
            &   T1,             & ! Temperature of upper layer [degC]
            &   T2,             & ! Temperature of lower layer [degC]
            &   hi,             & ! Ice thickness
            &   hs,             & ! Snow thickness
            &   Qtop,           & ! Energy flux available for surface melting [W/m2]
            &   Qbot,           & ! Energy flux available for bottom melting [W/m2]
            &   SWnet,          & ! Net shortwave flux [W/m^2]
            &   nonsolar,       & ! Latent and sensible heat flux and longwave radiation [W/m^2]
            &   dnonsolardT,    & ! Derivative of non-solar fluxes w.r.t. temperature [W/m^2/K]
            &   Tfw,            & ! Freezing temperature of the ocean
            &   albvisdir,      & ! Albedo VIS, direct/parallel
            &   albvisdif,      & ! Albedo VIS, diffuse
            &   albnirdir,      & ! Albedo NIR, direct/parallel
            &   albnirdif,      & ! Albedo NIR, diffuse
            &   doy,            & ! Day of the year
            &   lacc)

    INTEGER, INTENT(IN)    :: i_startidx_c, i_endidx_c, nbdim, kice
    REAL(wp),INTENT(IN)    :: pdtime
    REAL(wp),INTENT(INOUT) :: Tsurf      (nbdim,kice)
    REAL(wp),INTENT(INOUT) :: T1         (nbdim,kice)
    REAL(wp),INTENT(INOUT) :: T2         (nbdim,kice)
    REAL(wp),INTENT(IN)    :: hi         (nbdim,kice)
    REAL(wp),INTENT(IN)    :: hs         (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: Qtop       (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: Qbot       (nbdim,kice)
    REAL(wp),INTENT(IN)    :: SWnet      (nbdim,kice)
    REAL(wp),INTENT(IN)    :: nonsolar   (nbdim,kice)
    REAL(wp),INTENT(IN)    :: dnonsolardT(nbdim,kice)
    REAL(wp),INTENT(IN)    :: Tfw        (nbdim)
    REAL(wp),INTENT(OUT)   :: albvisdir  (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: albvisdif  (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: albnirdir  (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: albnirdif  (nbdim,kice)

    INTEGER, OPTIONAL,INTENT(IN)  :: doy
    LOGICAL, OPTIONAL,INTENT(IN)  :: lacc

    INTEGER :: jk, ji
    LOGICAL :: lzacc

    !-------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    IF (ltimer) CALL timer_start(timer_ice_fast)

    SELECT CASE (i_ice_therm)

    CASE (1)
      CALL set_ice_temp_zerolayer(i_startidx_c, i_endidx_c, nbdim, kice, pdtime, &
                            &   Tsurf, hi, hs, Qtop, Qbot, SWnet, nonsolar, dnonsolardT, Tfw, lacc=lzacc)

    CASE (2)
      CALL set_ice_temp_winton(i_startidx_c, i_endidx_c, nbdim, kice, pdtime, &
                    &   Tsurf, T1, T2, hi, hs, Qtop, Qbot, SWnet, nonsolar, dnonsolardT, Tfw, lacc=lzacc)

    CASE (3)
      IF ( .NOT. PRESENT(doy) ) THEN
        CALL finish(TRIM('mo_ice_interface:ice_fast'),'i_ice_therm = 3 not allowed in this context')
      ENDIF
      CALL set_ice_temp_zerolayer_analytical(i_startidx_c, i_endidx_c, nbdim, kice, &
            &   Tsurf, hi, hs, Qtop, Qbot, Tfw, doy, lacc=lzacc)

    CASE (4)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO ji = 1, kice
        !$ACC LOOP GANG VECTOR
        DO jk = 1, nbdim
          IF ( hi(jk,ji) > 0._wp ) THEN
            Tsurf(jk,ji) = min(0._wp, Tsurf(jk,ji) + (SWnet(jk,ji)+nonsolar(jk,ji) + ki/hi(jk,ji)*(Tf-Tsurf(jk,ji))) &
        &               / (ci*rhoi*hci_layer/pdtime-dnonsolardT(jk,ji)+ki/hi(jk,ji)))
          ELSE
            Tsurf(jk,ji) = Tf
          END IF
        END DO
      END DO
      !$ACC END PARALLEL

    END SELECT

    ! New albedo based on the new surface temperature
    CALL set_ice_albedo(i_startidx_c, i_endidx_c, nbdim, kice, Tsurf, hi, hs, &
      & albvisdir, albvisdif, albnirdir, albnirdif, lacc=lzacc)

    IF (ltimer) CALL timer_stop(timer_ice_fast)

   END SUBROUTINE ice_fast

END MODULE mo_ice_interface
