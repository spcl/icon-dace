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

#include "icon_definitions.inc"
   MODULE mo_hamocc_ocean_physics

    USE mo_kind,                         ONLY: wp
    USE mo_exception,                    ONLY: finish
    USE mo_param1_bgc,                   ONLY: n_bgctra, isco212
    USE mo_parallel_config,              ONLY: nproma
    USE mo_model_domain,                 ONLY: t_patch_3D, t_patch
    USE mo_grid_subset,                  ONLY: get_index_range, t_subset_range
    USE mo_timer,                        ONLY: timer_start, timer_stop, &
    &                                          timer_bgc_tracer_ab, timer_bgc_tot,&
    &                                          timers_level
    USE mo_ocean_tracer,                 ONLY: advect_ocean_tracers
    USE mo_ocean_tracer_zstar,           ONLY: advect_ocean_tracers_zstar
!    USE mo_ocean_tracer_GMRedi,          ONLY: advect_ocean_tracers_GMRedi
    USE mo_ocean_types,                  ONLY: t_hydro_ocean_state, &
    &                                          t_operator_coeff
    USE mo_ocean_tracer_transport_types, ONLY: t_tracer_collection, t_ocean_transport_state
    USE mo_ocean_surface_types,          ONLY: t_ocean_surface, t_atmos_for_ocean
    USE mo_sea_ice_types,                ONLY: t_sea_ice
    USE mo_run_config,                   ONLY: ltimer
    USE mo_ocean_nml,                    ONLY: Cartesian_Mixing, GMRedi_configuration, &
    &                                          lsediment_only, vert_cor_type, lfb_bgc_oce
    USE mo_hamocc_types,                 ONLY: t_hamocc_prog, t_hamocc_state
    USE mo_bgc_icon_comm,                ONLY: hamocc_state
    USE mo_dynamics_config,              ONLY: nold, nnew 
    USE mo_ocean_hamocc_couple_state, ONLY: t_ocean_to_hamocc_state, t_hamocc_to_ocean_state, &
      & t_hamocc_ocean_state
    USE mo_hamocc_diagnostics,     ONLY: get_monitoring 
    USE mo_bgc_bcond,              ONLY: ext_data_bgc, update_bgc_bcond
    USE mtime,                     ONLY: datetime
    USE mo_util_dbg_prnt,          ONLY: dbg_print
    USE mo_master_control,         ONLY: my_process_is_hamocc
    USE mo_control_bgc,            ONLY: bgc_zlevs, bgc_nproma

    USE mo_hamocc_nml,             ONLY: l_bgc_check,io_stdo_bgc
    USE mo_exception, ONLY: message_to_own_unit
    USE mo_hamocc_diagnostics,  ONLY: get_inventories
    USE mo_bgc_icon, ONLY: bgc_icon
    
    ! only temporary solution
    USE mo_ocean_tracer_dev,       ONLY: advect_ocean_tracers_dev, advect_ocean_tracers_GMRedi_zstar
    USE mo_ocean_physics_types,    ONLY: v_params
    USE mo_ocean_state,            ONLY: ocean_state

    ! OpenACC data movement - to be deleted later on
    USE mo_bgc_memory_types, ONLY: t_bgc_memory, t_sediment_memory, t_aggregates_memory, &
                                 & bgc_local_memory, sediment_local_memory, aggregates_memory

    USE mo_fortran_tools, ONLY: set_acc_host_or_device


    IMPLICIT NONE
    PRIVATE

    PUBLIC:: tracer_biochemistry_transport

    CONTAINS

  SUBROUTINE tracer_biochemistry_transport(hamocc_ocean_state, operators_coefficients, current_time)

    TYPE(t_hamocc_ocean_state), TARGET               :: hamocc_ocean_state
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    TYPE(datetime), POINTER, INTENT(in)              :: current_time
 
    
    TYPE(t_tracer_collection) , POINTER              :: old_tracer_collection, new_tracer_collection
    TYPE(t_ocean_to_hamocc_state), POINTER           :: ocean_to_hamocc_state
    TYPE(t_hamocc_to_ocean_state), POINTER           :: hamocc_to_ocean_state
    TYPE(t_patch_3d ),POINTER                        :: patch_3d
    TYPE(t_ocean_transport_state), POINTER           :: transport_state
    TYPE(t_hamocc_prog), POINTER                     :: hamocc_state_prog

    ! Local variables
    REAL(wp) :: pddpo(bgc_nproma, bgc_zlevs, hamocc_ocean_state%ocean_transport_state%patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: pddpo_new(bgc_nproma, bgc_zlevs, hamocc_ocean_state%ocean_transport_state%patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: ptiestu(bgc_nproma, bgc_zlevs, hamocc_ocean_state%ocean_transport_state%patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: ssh(bgc_nproma, hamocc_ocean_state%ocean_transport_state%patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: ssh_new(bgc_nproma, hamocc_ocean_state%ocean_transport_state%patch_3d%p_patch_2d(1)%alloc_cell_blocks)
 
    REAL(wp) :: stretch_e(nproma, hamocc_ocean_state%ocean_transport_state%patch_3d%p_patch_2d(1)%nblks_e)

    ! OpenACC data movement - to be deleted later on
    TYPE(t_bgc_memory), POINTER :: local_bgc_memory
    TYPE(t_sediment_memory), POINTER :: local_sediment_memory
    TYPE(t_aggregates_memory), POINTER :: local_aggregate_memory
    INTEGER :: local_memory_idx

    INTEGER :: i, jk, jb, start_index, end_index
    LOGICAL :: lzacc

    INTEGER, DIMENSION(:,:,:), POINTER :: cell_idx, cell_blk
    REAL(wp), DIMENSION(:,:), POINTER :: stretch_c
    TYPE(t_subset_range), POINTER :: edges_in_domain, all_edges 
    
    transport_state => hamocc_ocean_state%ocean_transport_state
    patch_3d => transport_state%patch_3d
    ocean_to_hamocc_state => hamocc_ocean_state%ocean_to_hamocc_state
    hamocc_to_ocean_state => hamocc_ocean_state%hamocc_to_ocean_state
    hamocc_state_prog => hamocc_state%p_prog(nold(1))
    edges_in_domain => patch_3d%p_patch_2D(1)%edges%in_domain
    all_edges => patch_3d%p_patch_2D(1)%edges%all

    ! 2023-01 dzo-DKRZ: Use OpenACC directives in called functions for GPU runs
    lzacc = .TRUE.

    !----------------------------------------------------------------------

    IF (vert_cor_type == 1) THEN ! z* coordinate
    
      stretch_c => ocean_to_hamocc_state%stretch_c(:,:)
      ! Adapt levels to changed stretching factors
      do jk = 1,bgc_zlevs
        pddpo(:,jk,:) = patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_c(:,jk,:) * &
              &           stretch_c(:,:) * patch_3d%wet_c(:,1,:)
        pddpo_new(:,jk,:) = patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_c(:,jk,:) * &
              &           ocean_to_hamocc_state%stretch_c_new(:,:) * patch_3d%wet_c(:,1,:)

        ptiestu(:,jk,:) = patch_3d%p_patch_1d(1)%depth_cellMiddle(:,jk,:) * &
              &           stretch_c(:,:) + ocean_to_hamocc_state%draftave(:,:)
      enddo

      ! compute stretch_e as the avrege of the two cells
      cell_idx  => patch_3D%p_patch_2D(1)%edges%cell_idx
      cell_blk  => patch_3D%p_patch_2D(1)%edges%cell_blk
      
!ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index, jb, i) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, start_index, end_index)
        DO i = start_index, end_index
          stretch_e(i, jb) = 0.5 * stretch_c(cell_idx(i, jb, 1), cell_blk(i, jb, 1)) + 0.5 * stretch_c(cell_idx(i, jb, 2), cell_blk(i, jb, 2))        
        END DO
      END DO 
!ICON_OMP_END_PARALLEL_DO    
      
      ! ssh is included in the adapted level thickness and depth
      ssh(:,:) = 0.0_wp
      ssh_new(:,:) = 0.0_wp

    ELSE
      pddpo(:,:,:) = patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:)
      pddpo_new(:,:,:) = patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_c(:,:,:)

      ptiestu(:,:,:) = patch_3d%p_patch_1d(1)%depth_cellMiddle(:,:,:)
      ssh(:,:) = ocean_to_hamocc_state%h_old(:,:)
      ssh_new(:,:) = ocean_to_hamocc_state%h_new(:,:)
    ENDIF


    IF (lsediment_only) THEN
      CALL offline_sediment(patch_3d, hamocc_ocean_state, current_time, ssh, pddpo, ptiestu)
      RETURN
    ENDIF

    !$ACC DATA COPY(pddpo, pddpo_new, ptiestu, ssh, ssh_new)

    ! OpenACC data movement - to be deleted later on
    local_memory_idx = 0
    local_bgc_memory => bgc_local_memory(local_memory_idx)
    local_sediment_memory => sediment_local_memory(local_memory_idx)
    local_aggregate_memory => aggregates_memory(local_memory_idx)


    IF (vert_cor_type == 1) THEN ! z* coordinate
      CALL dilute_hamocc_tracers_zstar(patch_3d, ocean_to_hamocc_state%top_dilution_coeff, &
                  &              ocean_to_hamocc_state%stretch_c, hamocc_state%p_prog(nold(1)))
    ELSE
      CALL dilute_hamocc_tracers(patch_3d, ocean_to_hamocc_state%top_dilution_coeff, hamocc_state%p_prog(nold(1)), lacc=lzacc)
    ENDIF
    !------------------------------------------------------------------------
    
    CALL update_bgc_bcond(patch_3d, ext_data_bgc, current_time, lacc=lzacc)   

    !------------------------------------------------------------------------
    ! call HAMOCC
    if(ltimer) call timer_start(timer_bgc_tot)
    CALL bgc_icon(patch_3d, hamocc_ocean_state, ssh, pddpo, ptiestu, lacc=lzacc)
    if(ltimer) call timer_stop(timer_bgc_tot)

    IF (l_bgc_check) THEN
      CALL message_to_own_unit('3. after bgc + fluxes and weathering', 'inventories', io_stdo_bgc)
      CALL get_inventories(hamocc_state, ssh, pddpo, hamocc_state%p_prog(nold(1))%tracer, patch_3d, 0._wp, 0._wp)
    ENDIF

    !------------------------------------------------------------------------
    ! transport tracers and diffuse them
    ! fill diffusion coefficients
    old_tracer_collection => hamocc_state%p_prog(nold(1))%tracer_collection
    new_tracer_collection => hamocc_state%p_prog(nnew(1))%tracer_collection

    DO i = 1, old_tracer_collection%no_of_tracers
      old_tracer_collection%tracer(i)%hor_diffusion_coeff => ocean_to_hamocc_state%hor_diffusion_coeff
      old_tracer_collection%tracer(i)%ver_diffusion_coeff => ocean_to_hamocc_state%ver_diffusion_coeff
    END DO

    start_timer(timer_bgc_tracer_ab,1)
    hamocc_state_prog => hamocc_state%p_prog(nnew(1))

!     IF ((my_process_is_hamocc()) .AND. lfb_bgc_oce) THEN
!       CALL finish("concurrent HAMOCC", "OBGC (phytoplankton and cyano) feedback is not possible with concurent.")
!     ENDIF

    IF (vert_cor_type == 1) THEN ! zstar transport routines
        IF (GMRedi_configuration == Cartesian_Mixing) THEN
          !! Note that zstar has no horizontal diffusion
          CALL advect_ocean_tracers_zstar(old_tracer_collection, new_tracer_collection, &
                 &  transport_state, operators_coefficients, stretch_e, ocean_to_hamocc_state%stretch_c, &
                 &  ocean_to_hamocc_state%stretch_c_new)
        ELSE
          IF (my_process_is_hamocc() ) THEN
            CALL finish("concurrent HAMOCC", "GMRedi is not possible at present")
          ELSE
            CALL  advect_ocean_tracers_GMRedi_zstar(old_tracer_collection, new_tracer_collection, &
                 &  ocean_state(1), transport_state, v_params, operators_coefficients, &
                 &  ocean_to_hamocc_state%stretch_c, stretch_e, ocean_to_hamocc_state%stretch_c_new)
          ENDIF
        ENDIF
    ELSE
      IF (GMRedi_configuration == Cartesian_Mixing) THEN
        CALL advect_ocean_tracers(old_tracer_collection, new_tracer_collection, transport_state, &
             &  operators_coefficients, lacc=lzacc)
      ELSE
        IF (my_process_is_hamocc() ) THEN
          CALL finish("concurrent HAMOCC", "GMRedi is not possible at present")
        ELSE
          CALL  advect_ocean_tracers_dev(old_tracer_collection, new_tracer_collection, &
            &  ocean_state(1), transport_state, v_params, operators_coefficients)
        ENDIF
      ENDIF
    ENDIF
    stop_timer(timer_bgc_tracer_ab,1)



    IF (l_bgc_check) THEN
      CALL message_to_own_unit('4. after transport', 'inventories', io_stdo_bgc)
      CALL get_inventories(hamocc_state, ssh_new, pddpo_new, hamocc_state%p_prog(nnew(1))%tracer, patch_3d, 0._wp, 0._wp)
    ENDIF
    !------------------------------------------------------------------------

    ! Update on the host the ocean variables modified inside HAMOCC
    !$ACC UPDATE HOST(hamocc_to_ocean_state%swr_fraction, hamocc_to_ocean_state%co2_flux) ASYNC(1)

     CALL get_monitoring( hamocc_state, hamocc_state%p_prog(nnew(1))%tracer, ssh_new, pddpo_new, patch_3d)
    !------------------------------------------------------------------------
    !$ACC WAIT(1)
    !$ACC END DATA

    END SUBROUTINE tracer_biochemistry_transport


  SUBROUTINE offline_sediment(patch_3d, hamocc_ocean_state, current_time, ssh, pddpo, ptiestu)

    TYPE(t_hamocc_ocean_state), TARGET               :: hamocc_ocean_state
    TYPE(datetime), POINTER, INTENT(in)              :: current_time
    TYPE(t_patch_3D ),TARGET, INTENT(IN)             :: patch_3D

    REAL(wp), INTENT(IN) :: pddpo(bgc_nproma, bgc_zlevs, patch_3d%p_patch_2d(1)%nblks_c)
    REAL(wp), INTENT(IN) :: ptiestu(bgc_nproma, bgc_zlevs, patch_3d%p_patch_2d(1)%nblks_c)
    REAL(wp), INTENT(IN) :: ssh(bgc_nproma, patch_3d%p_patch_2d(1)%nblks_c)

    CALL update_bgc_bcond( patch_3d, ext_data_bgc,  current_time)
 
    !------------------------------------------------------------------------
    ! call HAMOCC
    if(ltimer) call timer_start(timer_bgc_tot)
    CALL bgc_icon(patch_3d, hamocc_ocean_state, ssh, pddpo, ptiestu)
    if(ltimer) call timer_stop(timer_bgc_tot)

    !------------------------------------------------------------------------

    END SUBROUTINE offline_sediment


    SUBROUTINE DILUTE_HAMOCC_TRACERS(p_patch_3D, top_dilution_coeff, hamocc_state_prog, lacc)
    ! HAMOCC tracer dilution using old and new h

    TYPE(t_patch_3D ),TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_hamocc_prog), INTENT(inout)       :: hamocc_state_prog
    REAL(wp), INTENT(in)                     :: top_dilution_coeff(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    LOGICAL, INTENT(IN), OPTIONAL            :: lacc

    ! Local variables
    INTEGER:: jb, jc, i_bgc_tra, i_startidx_c, i_endidx_c
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: p_patch
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    p_patch => p_patch_3D%p_patch_2D(1)
    all_cells => p_patch%cells%all

    DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx_c, i_endidx_c
            IF (p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb) > 0) THEN
                !$ACC LOOP SEQ
                DO i_bgc_tra = 1, n_bgctra
                     hamocc_state_prog%tracer(jc,1,jb,i_bgc_tra)  = hamocc_state_prog%tracer(jc,1,jb,i_bgc_tra) &
    &                                                               * top_dilution_coeff(jc,jb)
                ENDDO
            endif
        ENDDO
        !$ACC END PARALLEL
    ENDDO

   END SUBROUTINE DILUTE_HAMOCC_TRACERS

    SUBROUTINE dilute_hamocc_tracers_zstar(p_patch_3D, top_dilution_coeff, stretch_c, hamocc_state_prog)
    ! HAMOCC tracer dilution with interpolation of subsurface levels due to changed level thicknesses

    TYPE(t_patch_3D ),TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_hamocc_prog), INTENT(inout)       :: hamocc_state_prog
    REAL(wp), INTENT(in)                     :: top_dilution_coeff(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp), INTENT(in)                     :: stretch_c(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    ! Local variables
    INTEGER:: jb, jc, jk, i_bgc_tra, i_startidx_c, i_endidx_c
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER        :: p_patch
    REAL(wp)                      :: h_old(bgc_zlevs), h_new(bgc_zlevs), h_change(bgc_zlevs)
    INTEGER                       :: nlevs

    p_patch => p_patch_3D%p_patch_2D(1)
    all_cells => p_patch%cells%all

    DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
            IF (p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb) > 0) THEN

                nlevs = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)              

                ! old thickness of cells = prism_thick * old_stretch
                ! tdc = old_stretch / stretch_c

                DO jk = 1,nlevs
                  h_old(jk) = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,jk,jb) &
                            & * top_dilution_coeff(jc,jb) * stretch_c(jc,jb)
                  h_new(jk) = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,jk,jb) &
                            & * stretch_c(jc,jb)
                ENDDO


                ! The amount by which the layers from jk and below have grown
                ! multiplied with -1, so if the lowest layer has grown by 1m,
                ! then the value of h_change(nlevs) = -1

                h_change(nlevs) = h_new(nlevs)-h_old(nlevs)
                DO jk = nlevs-1,1,-1
                  h_change(jk) = h_change(jk+1) + h_new(jk) - h_old(jk)
                ENDDO
                
                DO i_bgc_tra = 1, n_bgctra

                  IF (top_dilution_coeff(jc,jb) > 1.0_wp) THEN

                    ! First everything is put into the first layer

                    hamocc_state_prog%tracer(jc,1,jb,i_bgc_tra) = h_old(1) / (h_old(1) + h_change(1)) &
                        &                               * hamocc_state_prog%tracer(jc,1,jb,i_bgc_tra)


                    ! Now follows an interpolation to the changed level
                    ! distribution.
                    ! Due to the smaller cell thickness only the layer
                    ! below level jk influences the new concentration in jk.
                    ! Cell thickness has become smaller, so bottom cell
                    ! concentration is unaffected.

                    DO jk = 1,nlevs-1
                      hamocc_state_prog%tracer(jc,jk,jb,i_bgc_tra) = ((h_new(jk) + h_change(jk+1)) &
                        &   * hamocc_state_prog%tracer(jc,jk,jb,i_bgc_tra) - h_change(jk+1) &
                        &   * hamocc_state_prog%tracer(jc,jk+1,jb,i_bgc_tra)) &
                        &   / h_new(jk)
                    ENDDO          

                  ELSEIF (top_dilution_coeff(jc,jb) < 1.0_wp) THEN


                    ! Now first do the interpolation, as new first level
                    ! concentration is independent of subsurface cell
                    ! concentrations, but old value is needed for the
                    ! interpolation of the second level.
                    ! Cells have become thicker, hence only the level above jk
                    ! influences the new concentration in 

                    DO jk = nlevs,2,-1
                      hamocc_state_prog%tracer(jc,jk,jb,i_bgc_tra) = &
                        & ( h_change(jk) * hamocc_state_prog%tracer(jc,jk-1,jb,i_bgc_tra) &
                        &   + (h_new(jk) - h_change(jk))&
                        &   * hamocc_state_prog%tracer(jc,jk,jb,i_bgc_tra) ) &
                        &   / h_new(jk)
                    ENDDO   


                    ! Now adapt the first level. Basically with the same formula
                    ! as above, where tracer(jk-1) = 0.

                     hamocc_state_prog%tracer(jc,1,jb,i_bgc_tra) = (h_new(1) - h_change(1)) &
                       &   * hamocc_state_prog%tracer(jc,1,jb,i_bgc_tra) / h_new(1)

                  ENDIF

                ENDDO
            endif
        ENDDO
    ENDDO

   END SUBROUTINE dilute_hamocc_tracers_zstar


   END MODULE
