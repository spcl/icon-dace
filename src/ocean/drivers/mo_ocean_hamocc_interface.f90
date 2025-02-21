! Contains the main stepping routine the 3-dim hydrostatic ocean model.
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

!----------------------------
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


!--------------------------------------------------
! timers definition
!needs:
!   USE mo_timer, ONLY: timer_start, timer_stop, timers_level, <timers_names>...
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

! DSL definitions 










!---------------------
! block definitions

!---------------------
! mappings


















!---------------------
! connectivity









!---------------------
! generic types




!---------------------
! shortcuts







































!---------------------
! Upper-lower case
!----------------------------
MODULE mo_ocean_hamocc_interface
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message, finish
  USE mo_impl_constants,         ONLY: max_char_length, success
  USE mo_parallel_config,        ONLY: nproma
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d
  USE mo_ocean_nml,              ONLY: n_zlev, no_tracer, lhamocc, &
    &  Cartesian_Mixing, GMRedi_configuration
  USE mo_dynamics_config,        ONLY: nold, nnew
  USE mo_ocean_types,              ONLY: t_hydro_ocean_state, &
    & t_operator_coeff
  USE mo_ocean_tracer_transport_types, ONLY: t_tracer_collection, t_ocean_transport_state
  USE mo_ocean_surface_types,    ONLY: t_ocean_surface, t_atmos_for_ocean
  USE mo_sea_ice_types,          ONLY: t_atmos_fluxes, t_sea_ice
  USE mo_ocean_physics_types,    ONLY: t_ho_params  
  USE mo_master_config,          ONLY: isRestart
  USE mo_hamocc_ocean_physics,   ONLY: tracer_biochemistry_transport
  USE mo_ocean_hamocc_couple_state, ONLY: t_ocean_to_hamocc_state, t_hamocc_to_ocean_state, &
    & t_ocean_transport_state, t_hamocc_ocean_state, hamocc_ocean_state
  USE mtime,                     ONLY: datetime    
  USE mo_construct_icon_hamocc,  ONLY: construct_icon_hamocc, destruct_icon_hamocc, init_icon_hamocc
  USE mo_ext_data_types,         ONLY: t_external_data
 
  USE mo_timer,                  ONLY: timer_start, timer_stop, timer_total, timer_exchange_ocean_hamocc, &
    timers_level, ltimer

  USE mo_util_dbg_prnt,          ONLY: dbg_print
  
  USE mo_master_control,         ONLY: process_exists, hamocc_process, my_process_is_hamocc, &
    & ocean_process, my_process_is_ocean

  USE iso_c_binding,             ONLY: c_loc
  
  USE  mo_ocean_hamocc_communication, ONLY:   &
    setup_ocean_2_hamocc_communication, &
    setup_hamocc_2_ocean_communication, &
    exchange_data_ocean_2_hamocc, &
    exchange_data_hamocc_2_ocean, &
    free_ocean_hamocc_communication

  USE mo_sync,                   ONLY: sync_c, sync_e, sync_patch_array, sync_patch_array_mult

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: ocean_to_hamocc_construct, ocean_to_hamocc_init, ocean_to_hamocc_end, ocean_to_hamocc_interface
  PUBLIC  :: hamocc_to_ocean_init, hamocc_to_ocean_end, hamocc_to_ocean_interface
  !-------------------------------------------------------------------------
  
CONTAINS
!-------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE ocean_to_hamocc_construct(patch_3d, ext_data)

    TYPE(t_patch_3d ), INTENT(in)     :: patch_3d
    TYPE(t_external_data), TARGET, INTENT(inout) :: ext_data
   
    IF(.not. lhamocc) return
   
    IF (process_exists(hamocc_process)) RETURN ! this will be construced in the hamocc side
    
    CALL construct_icon_hamocc(patch_3d, ext_data)

  END SUBROUTINE ocean_to_hamocc_construct
  !--------------------------------------------datetime-----------------------------

  !-------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE hamocc_to_ocean_init()

   !----------------------------------------------
    IF (process_exists(ocean_process)) THEN
      ! we run in a coupled ocean-hamocc setup 
      CALL setup_hamocc_2_ocean_communication(hamocc_ocean_state%patch_3D%p_patch_2d(1), n_zlev)      
 !     CALL exchange_ocean_to_hamocc_state()
      ! sync the input from the ocean, as this is sent only for owned cells/edges
 !     CALL sync_hamocc_input()
    ENDIF
    
    CALL init_icon_hamocc (hamocc_ocean_state)  
    
  END SUBROUTINE hamocc_to_ocean_init
  !--------------------------------------------datetime-----------------------------

  !-------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE ocean_to_hamocc_init(patch_3d, ocean_state, &
    & p_as, sea_ice, p_oce_sfc, p_phys_param)

    TYPE(t_patch_3d ), INTENT(in) , TARGET :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET     :: ocean_state
    TYPE(t_atmos_for_ocean)               :: p_as
    TYPE (t_sea_ice)                      :: sea_ice
    TYPE(t_ocean_surface)                 :: p_oce_sfc
     TYPE (t_ho_params)                   :: p_phys_param
 
    TYPE(t_ocean_transport_state), TARGET   :: my_transport_state
    TYPE(t_patch), POINTER :: patch_2d
    INTEGER :: alloc_cell_blocks, nblks_e

  
    IF(.not. lhamocc) return
   
    CALL message("ocean_to_hamocc_init", "...")
    patch_2d => patch_3d%p_patch_2d(1)
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
   !----------------------------------------------
    ! dummy allocation for using the same communication patterns
    ALLOCATE(my_transport_state%h_new(nproma,alloc_cell_blocks), &
             my_transport_state%h_old(nproma,alloc_cell_blocks), &
             my_transport_state%mass_flux_e(nproma,n_zlev,nblks_e),      &
             my_transport_state%vn(nproma,n_zlev,nblks_e),      &
             my_transport_state%w(nproma,n_zlev+1,alloc_cell_blocks))
    my_transport_state%h_new = 0.0_wp
    my_transport_state%h_old = 0.0_wp
    my_transport_state%mass_flux_e = 0.0_wp
    my_transport_state%vn    = 0.0_wp
    my_transport_state%w     = 0.0_wp
    my_transport_state%patch_3d => patch_3d
    !----------------------------------------------
    
    CALL fill_ocean_to_hamocc_interface(ocean_state, my_transport_state, p_oce_sfc, p_as, &
      & sea_ice, p_phys_param)
    !----------------------------------------------
    IF (process_exists(hamocc_process)) THEN
      ! we run in a coupled ocean-hamocc setup 
!       CALL message("setup_ocean_2_hamocc_communication", "...")
      CALL setup_ocean_2_hamocc_communication(patch_2d, n_zlev)      
!      CALL exchange_ocean_to_hamocc_state()
      
    ELSE
      CALL init_icon_hamocc (hamocc_ocean_state)  
    ENDIF
    
    DEALLOCATE(my_transport_state%h_new,       &
             my_transport_state%h_old,         &
             my_transport_state%mass_flux_e,   &
             my_transport_state%vn,            &
             my_transport_state%w)

  END SUBROUTINE ocean_to_hamocc_init
  !--------------------------------------------datetime-----------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE ocean_to_hamocc_end()

    IF(.not. lhamocc) return

    IF (process_exists(hamocc_process)) THEN
      CALL free_ocean_hamocc_communication()
    ELSE
      CALL destruct_icon_hamocc()
    ENDIF
   
  END SUBROUTINE ocean_to_hamocc_end
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE hamocc_to_ocean_end()

    IF (process_exists(ocean_process)) THEN
      CALL free_ocean_hamocc_communication()
    ENDIF
      
  END SUBROUTINE hamocc_to_ocean_end
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE ocean_to_hamocc_interface(ocean_state, transport_state, p_oce_sfc, p_as, &
      & sea_ice, p_phys_param, operators_coefficients, current_time)
    TYPE(t_hydro_ocean_state), TARGET, INTENT(in)    :: ocean_state
    TYPE(t_ocean_transport_state), TARGET            :: transport_state
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_sea_ice),          INTENT(inout)          :: sea_ice
    TYPE(t_ocean_surface)                            :: p_oce_sfc
    TYPE(t_ho_params)                                :: p_phys_param
    TYPE(t_operator_coeff),   INTENT(inout)          :: operators_coefficients
    TYPE(datetime), POINTER, INTENT(in)              :: current_time

   
    IF(.not. lhamocc) return
    
    CALL fill_ocean_to_hamocc_interface(ocean_state, transport_state, p_oce_sfc, p_as, &
      & sea_ice, p_phys_param)

!     CALL sync_patch_array(sync_e,  transport_state%patch_3d%p_patch_2d(1),  &
!       & transport_state%mass_flux_e)
 
!     CALL dbg_print('mass_flux_e All'       , transport_state%mass_flux_e, "from ocean", 1,  &
!       & transport_state%patch_3d%p_patch_2d(1)%edges%all)
!     CALL dbg_print('mass_flux_e own'       , transport_state%mass_flux_e, "from ocean", 1,  &
!       & transport_state%patch_3d%p_patch_2d(1)%edges%owned)
!     CALL dbg_print('mass_flux_e dom'       , transport_state%mass_flux_e, "from ocean", 1,  &
!       & transport_state%patch_3d%p_patch_2d(1)%edges%in_domain)
     
    !------------------------------------------------------------------------
    IF (process_exists(hamocc_process)) THEN
      ! concurrent case
      CALL exchange_ocean_to_hamocc_state()
      CALL exchange_hamocc_to_ocean_state()
 
      CALL sync_ocean_input()
     
    ELSE
      ! sequential
      CALL tracer_biochemistry_transport(hamocc_ocean_state, operators_coefficients, current_time)
    ENDIF
    
  END SUBROUTINE ocean_to_hamocc_interface
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE hamocc_to_ocean_interface()
        
    !------------------------------------------------------------------------
    IF (process_exists(ocean_process)) THEN
      ! concurrent case
      CALL exchange_ocean_to_hamocc_state()
      
      CALL exchange_hamocc_to_ocean_state()
      ! sync the input from the ocean, as this is sent only for owned cells/edges
      CALL sync_hamocc_input()

    ENDIF
    
  END SUBROUTINE hamocc_to_ocean_interface
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE fill_ocean_to_hamocc_interface(ocean_state, transport_state, p_oce_sfc, p_as, &
      & sea_ice, p_phys_param)
    TYPE(t_hydro_ocean_state), TARGET, INTENT(in)    :: ocean_state
    TYPE(t_ocean_transport_state), TARGET            :: transport_state
    TYPE(t_atmos_for_ocean),  INTENT(inout)          :: p_as
    TYPE(t_sea_ice),          INTENT(inout)          :: sea_ice
    TYPE(t_ocean_surface)                            :: p_oce_sfc
    TYPE(t_ho_params)                                :: p_phys_param
   
    TYPE(t_ocean_to_hamocc_state), POINTER           :: ocean_to_hamocc_state
    TYPE(t_hamocc_to_ocean_state), POINTER           :: hamocc_to_ocean_state

    !$ACC UPDATE DEVICE(ocean_state%p_diag%swr_frac) ASYNC(1)
    !$ACC WAIT(1) ! can be removed when all ACC compute regions are ASYNC(1)

    ocean_to_hamocc_state => hamocc_ocean_state%ocean_to_hamocc_state
    hamocc_to_ocean_state => hamocc_ocean_state%hamocc_to_ocean_state
    hamocc_ocean_state%ocean_transport_state => transport_state 
    hamocc_ocean_state%patch_3d => transport_state%patch_3d
    
    ocean_to_hamocc_state%top_dilution_coeff => p_oce_sfc%top_dilution_coeff
    ocean_to_hamocc_state%h_old              => transport_state%h_old
    ocean_to_hamocc_state%h_new              => transport_state%h_new
    ocean_to_hamocc_state%h_old_withIce      =>  ocean_state%p_prog(nold(1))%h
    ocean_to_hamocc_state%ice_concentration_sum => sea_ice%concSum    
    ocean_to_hamocc_state%temperature        => ocean_state%p_prog(nold(1))%tracer(:,:,:,1)
    ocean_to_hamocc_state%salinity           => ocean_state%p_prog(nold(1))%tracer(:,:,:,2)
    ocean_to_hamocc_state%press_hyd          => ocean_state%p_diag%press_hyd   ! (agg)
    ocean_to_hamocc_state%hor_diffusion_coeff => p_phys_param%TracerDiffusion_coeff(:,:,:,2)
    ocean_to_hamocc_state%ver_diffusion_coeff => p_phys_param%a_tracer_v(:,:,:,2)
    ocean_to_hamocc_state%short_wave_flux    => p_as%fswr  ! p_oce_sfc%HeatFlux_ShortWave
    ocean_to_hamocc_state%wind10m            => p_as%fu10
    ocean_to_hamocc_state%co2_mixing_ratio   => p_as%co2 
   
    hamocc_to_ocean_state%co2_flux           => p_as%co2flx
    hamocc_to_ocean_state%swr_fraction       => ocean_state%p_diag%swr_frac

    ! Variables for zstar calculations
!     ocean_to_hamocc_state%eta_c              => ocean_state%p_prog(nold(1))%eta_c
    ocean_to_hamocc_state%stretch_c          => ocean_state%p_prog(nold(1))%stretch_c
    ocean_to_hamocc_state%stretch_c_new      => ocean_state%p_prog(nnew(1))%stretch_c
    ocean_to_hamocc_state%draftave           => sea_ice%draftave    

    !$ACC ENTER DATA COPYIN(hamocc_ocean_state, hamocc_ocean_state%hamocc_to_ocean_state, hamocc_ocean_state%hamocc_to_ocean_state%swr_fraction) &
    !$ACC   COPYIN(hamocc_ocean_state%hamocc_to_ocean_state%co2_flux)

  END SUBROUTINE fill_ocean_to_hamocc_interface
  !------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE exchange_ocean_to_hamocc_state()
 
    TYPE(t_ocean_to_hamocc_state), POINTER :: ocean_to_hamocc_state
    TYPE(t_ocean_transport_state), POINTER :: ocean_transport_state

    REAL(wp), POINTER   :: &
      &  top_dilution_coeff(:,:),              &
      &  h_old(:,:),                           &
      &  h_new(:,:),                           &
      &  h_old_withIce(:,:),                   &
      &  ice_concentration_sum(:,:),           &
      &  temperature(:,:,:),                   &
      &  salinity(:,:,:),                      &
      &  ver_diffusion_coeff(:,:,:),           &
      &  short_wave_flux(:,:),                 &
      &  wind10m(:,:),                         &
      &  co2_mixing_ratio(:,:),                &
      &  mass_flux_e(:,:,:),                   &
      &  vn(:,:,:),                            &
      &  w(:,:,:) ,                            &
      &  press_hyd(:,:,:),                     & 
      &  stretch_c(:,:),                       &   
      &  stretch_c_new(:,:),                   &
      &  draftave(:,:)  

    IF (timers_level >= 1) CALL timer_start(timer_exchange_ocean_hamocc)

    ocean_to_hamocc_state => hamocc_ocean_state%ocean_to_hamocc_state
    ocean_transport_state => hamocc_ocean_state%ocean_transport_state
    
    
    top_dilution_coeff        => ocean_to_hamocc_state%top_dilution_coeff    
    h_old                     => ocean_to_hamocc_state%h_old                
    h_new                     => ocean_to_hamocc_state%h_new                 
    h_old_withIce             => ocean_to_hamocc_state%h_old_withIce              
    ice_concentration_sum     => ocean_to_hamocc_state%ice_concentration_sum 
    temperature               => ocean_to_hamocc_state%temperature         
    salinity                  => ocean_to_hamocc_state%salinity            
    ver_diffusion_coeff       => ocean_to_hamocc_state%ver_diffusion_coeff 
    short_wave_flux           => ocean_to_hamocc_state%short_wave_flux       
    wind10m                   => ocean_to_hamocc_state%wind10m               
    co2_mixing_ratio          => ocean_to_hamocc_state%co2_mixing_ratio      
    mass_flux_e               => ocean_transport_state%mass_flux_e         
    vn                        => ocean_transport_state%vn                  
    w                         => ocean_transport_state%w
    press_hyd                 => ocean_to_hamocc_state%press_hyd
    stretch_c                 => ocean_to_hamocc_state%stretch_c    
    stretch_c_new             => ocean_to_hamocc_state%stretch_c_new
    draftave                  => ocean_to_hamocc_state%draftave     
    
    
!     CALL exchange_data_ocean_2_hamocc(                                      &
!       &  c_loc(ocean_to_hamocc_state%top_dilution_coeff(1,1)),              &
!       &  c_loc(ocean_to_hamocc_state%h_old,(1,1)),                          &
!       &  c_loc(ocean_to_hamocc_state%h_new(1,1)),                           &
!       &  c_loc(ocean_to_hamocc_state%ice_concentration_sum(1,1)),           &
!       &  c_loc(ocean_to_hamocc_state%temperature(1,1,1)),                   &
!       &  c_loc(ocean_to_hamocc_state%salinity(1,1,1)),                      &
!       &  c_loc(ocean_to_hamocc_state%ver_diffusion_coeff(1,1,1)),           &
!       &  c_loc(ocean_to_hamocc_state%short_wave_flux(1,1)),                 &
!       &  c_loc(ocean_to_hamocc_state%wind10m(1,1)),                         &
!       &  c_loc(ocean_to_hamocc_state%co2_mixing_ratio(1,1)),                &
!       &  c_loc(ocean_transport_state%mass_flux_e(1,1,1)),                   &
!       &  c_loc(ocean_transport_state%vn(1,1,1)),                            &
!       &  c_loc(ocean_transport_state%w(1,1,1)))
      
    CALL exchange_data_ocean_2_hamocc(              &
      &  c_loc(top_dilution_coeff(1,1)),            &
      &  c_loc(h_old(1,1)),                         &
      &  c_loc(h_new(1,1)),                         &
      &  c_loc(h_old_withIce(1,1)),                 &
      &  c_loc(ice_concentration_sum(1,1)),         &
      &  c_loc(temperature(1,1,1)),                 &
      &  c_loc(salinity(1,1,1)),                    &
      &  c_loc(ver_diffusion_coeff(1,1,1)),         &
      &  c_loc(short_wave_flux(1,1)),               &
      &  c_loc(wind10m(1,1)),                       &
      &  c_loc(co2_mixing_ratio(1,1)),              &
      &  c_loc(mass_flux_e(1,1,1)),                 &
      &  c_loc(vn(1,1,1)),                          &
      &  c_loc(w(1,1,1)),                           &
      &  c_loc(press_hyd(1,1,1)),                   &
      &  c_loc(stretch_c(1,1)),                     &
      &  c_loc(stretch_c_new(1,1)),                 &
      &  c_loc(draftave(1,1))  &
      &  )
      
  
    IF (timers_level >= 1) CALL timer_stop(timer_exchange_ocean_hamocc)

  END SUBROUTINE exchange_ocean_to_hamocc_state
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE exchange_hamocc_to_ocean_state()
 
    TYPE(t_hamocc_to_ocean_state), POINTER           :: hamocc_to_ocean_state

    REAL(wp), POINTER   ::     &
      &  swr_frac(:,:,:),      &
      &  co2_flux(:,:)

    
    IF (timers_level >= 1) CALL timer_start(timer_exchange_ocean_hamocc)

    hamocc_to_ocean_state => hamocc_ocean_state%hamocc_to_ocean_state

    co2_flux => hamocc_to_ocean_state%co2_flux    
    swr_frac => hamocc_to_ocean_state%swr_fraction
    
    CALL exchange_data_hamocc_2_ocean(                &
      &  c_loc(co2_flux(1,1)),                        &
      &  c_loc(swr_frac(1,1,1)))
     
    IF (timers_level >= 1) CALL timer_stop(timer_exchange_ocean_hamocc)

    CALL dbg_print('swr_fraction', hamocc_to_ocean_state%swr_fraction , "from hamocc", 2,  &
       & hamocc_ocean_state%patch_3d%p_patch_2d(1)%cells%owned)


  END SUBROUTINE exchange_hamocc_to_ocean_state
  !-------------------------------------------------------------------------
      
      
  !-------------------------------------------------------------------------
  ! sync the veriables received from hamocc, since YAXT does not communicate halos
  SUBROUTINE sync_ocean_input()

    TYPE(t_patch), POINTER :: patch_2d
     
    patch_2d => hamocc_ocean_state%patch_3D%p_patch_2d(1)

    CALL sync_patch_array(sync_c, patch_2d,  hamocc_ocean_state%hamocc_to_ocean_state%swr_fraction)
    
  END SUBROUTINE sync_ocean_input
  !-------------------------------------------------------------------------
      
  !-------------------------------------------------------------------------
  ! sync the veriables received from the ocean, since YAXT does not communicate halos
  SUBROUTINE sync_hamocc_input()

    TYPE(t_patch), POINTER :: patch_2d
    INTEGER :: alloc_cell_blocks, nblks_e
    REAL(wp), POINTER :: gather_cells_2d(:,:,:)
    TYPE(t_ocean_transport_state), POINTER :: transport_state
     
    patch_2d => hamocc_ocean_state%patch_3D%p_patch_2d(1)
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    transport_state => hamocc_ocean_state%ocean_transport_state

    !----------------------------------------------
    ! use a copy to a 3d strucure to do the communication of the 2d
    ! Not nice, but probably better...
    ALLOCATE(gather_cells_2d(nproma,11,alloc_cell_blocks))
    gather_cells_2d(:,1,:) = hamocc_ocean_state%ocean_to_hamocc_state%top_dilution_coeff(:,:)
    gather_cells_2d(:,2,:) = hamocc_ocean_state%ocean_to_hamocc_state%h_old(:,:)                 
    gather_cells_2d(:,3,:) = hamocc_ocean_state%ocean_to_hamocc_state%h_new(:,:)                 
    gather_cells_2d(:,4,:) = hamocc_ocean_state%ocean_to_hamocc_state%h_old_withIce(:,:)                 
    gather_cells_2d(:,5,:) = hamocc_ocean_state%ocean_to_hamocc_state%ice_concentration_sum(:,:) 
    gather_cells_2d(:,6,:) = hamocc_ocean_state%ocean_to_hamocc_state%short_wave_flux(:,:)      
    gather_cells_2d(:,7,:) = hamocc_ocean_state%ocean_to_hamocc_state%wind10m(:,:)               
    gather_cells_2d(:,8,:) = hamocc_ocean_state%ocean_to_hamocc_state%co2_mixing_ratio(:,:)
    gather_cells_2d(:,9,:) = hamocc_ocean_state%ocean_to_hamocc_state%stretch_c(:,:)    
    gather_cells_2d(:,10,:) = hamocc_ocean_state%ocean_to_hamocc_state%stretch_c_new(:,:)
    gather_cells_2d(:,11,:) = hamocc_ocean_state%ocean_to_hamocc_state%draftave (:,:)    
        
    CALL sync_patch_array(sync_c, patch_2d, gather_cells_2d)

    hamocc_ocean_state%ocean_to_hamocc_state%top_dilution_coeff(:,:)    = gather_cells_2d(:,1,:)
    hamocc_ocean_state%ocean_to_hamocc_state%h_old(:,:)                 = gather_cells_2d(:,2,:)
    hamocc_ocean_state%ocean_to_hamocc_state%h_new(:,:)                 = gather_cells_2d(:,3,:)
    hamocc_ocean_state%ocean_to_hamocc_state%h_old_withIce(:,:)         = gather_cells_2d(:,4,:)
    hamocc_ocean_state%ocean_to_hamocc_state%ice_concentration_sum(:,:) = gather_cells_2d(:,5,:)
    hamocc_ocean_state%ocean_to_hamocc_state%short_wave_flux(:,:)       = gather_cells_2d(:,6,:)
    hamocc_ocean_state%ocean_to_hamocc_state%wind10m(:,:)               = gather_cells_2d(:,7,:)
    hamocc_ocean_state%ocean_to_hamocc_state%co2_mixing_ratio(:,:)      = gather_cells_2d(:,8,:)
    hamocc_ocean_state%ocean_to_hamocc_state%stretch_c(:,:)             = gather_cells_2d(:,9,:)
    hamocc_ocean_state%ocean_to_hamocc_state%stretch_c_new(:,:)         = gather_cells_2d(:,10,:)
    hamocc_ocean_state%ocean_to_hamocc_state%draftave (:,:)             = gather_cells_2d(:,11,:)
    
    DEALLOCATE(gather_cells_2d)
    
    ! sync the 3D fields
    CALL sync_patch_array_mult(sync_c, patch_2d, 3,  &
      & hamocc_ocean_state%ocean_to_hamocc_state%temperature,           &
      & hamocc_ocean_state%ocean_to_hamocc_state%salinity,              &
      & hamocc_ocean_state%ocean_to_hamocc_state%press_hyd)
      
    CALL sync_patch_array_mult(sync_c, patch_2d, 2,  &
      & hamocc_ocean_state%ocean_to_hamocc_state%ver_diffusion_coeff,   &
      & hamocc_ocean_state%ocean_transport_state%w)

!     CALL dbg_print('mass_flux_e All'       , transport_state%mass_flux_e, "before sync", 1,  &
!       & transport_state%patch_3d%p_patch_2d(1)%edges%all)
!     CALL dbg_print('mass_flux_e own'       , transport_state%mass_flux_e, "before sync", 1,  &
!       & transport_state%patch_3d%p_patch_2d(1)%edges%owned)
!     CALL dbg_print('mass_flux_e dom'       , transport_state%mass_flux_e, "before sync", 1,  &
!       & transport_state%patch_3d%p_patch_2d(1)%edges%in_domain)

    CALL sync_patch_array_mult(sync_e, patch_2d, 2,  &
      & hamocc_ocean_state%ocean_transport_state%vn,                    &
      & hamocc_ocean_state%ocean_transport_state%mass_flux_e)
  
!     CALL dbg_print('mass_flux_e All'       , transport_state%mass_flux_e, "after sync", 1,  &
!       & transport_state%patch_3d%p_patch_2d(1)%edges%all)
!     CALL dbg_print('mass_flux_e own'       , transport_state%mass_flux_e, "after sync", 1,  &
!       & transport_state%patch_3d%p_patch_2d(1)%edges%owned)
!     CALL dbg_print('mass_flux_e dom'       , transport_state%mass_flux_e, "after sync", 1,  &
!       & transport_state%patch_3d%p_patch_2d(1)%edges%in_domain)
        
  END SUBROUTINE sync_hamocc_input
  !-------------------------------------------------------------------------

END MODULE mo_ocean_hamocc_interface
