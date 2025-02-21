! Contains the implementation of the tracer transport routines for the ICON ocean model.
! This comprises advection and diffusion in horizontal and vertical direction.
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
!=============================================================================================
MODULE mo_ocean_GM_Redi
  !-------------------------------------------------------------------------
  USE mo_kind,                      ONLY: wp
  USE mo_math_types,                ONLY: t_cartesian_coordinates
  USE mo_impl_constants,            ONLY: sea_boundary, sea, min_dolic
  USE mo_math_constants,            ONLY: pi, dbl_eps
  USE mo_physical_constants,        ONLY: grav, sal_ref, rho_inv, a_t, b_s, &
    & sitodbar, sfc_press_bar
  USE mo_ocean_nml,                 ONLY: n_zlev, no_tracer,                    &
    & GMRedi_configuration,&
    & GMRedi_combined, GM_only, Redi_only,Cartesian_Mixing,&
    & k_tracer_dianeutral_parameter, k_tracer_isoneutral_parameter, &
    & k_tracer_GM_kappa_parameter,EOS_TYPE,&
    & GMRedi_configuration,GMRedi_combined, GM_only,Redi_only,Cartesian_Mixing, &
    & tapering_scheme,tapering_DanaMcWilliams,tapering_Large,tapering_Griffies, &
    & vert_mix_type, vmix_pp, vmix_tke, vmix_idemix_tke,OceanReferenceDensity,ReferencePressureIndbars, & ! by_Oliver
    & S_max, S_d, S_critical, c_speed, GMRedi_usesRelativeMaxSlopes,vert_cor_type,            &
    & RossbyRadius_max, RossbyRadius_min,switch_off_diagonal_vert_expl,lvertical_GM,Nmin,     & ! by_Oliver: added lvertical_GM, Nmin
    & GMREDI_COMBINED_DIAGNOSTIC,GM_INDIVIDUAL_DIAGNOSTIC,REDI_INDIVIDUAL_DIAGNOSTIC,&
    &TEST_MODE_REDI_ONLY,TEST_MODE_GM_ONLY,LinearThermoExpansionCoefficient, LinearHalineContractionCoefficient,&
    & SWITCH_OFF_TAPERING,SWITCH_ON_REDI_BALANCE_DIAGONSTIC, SWITCH_ON_TAPERING_HORIZONTAL_DIFFUSION,&
    &SLOPE_CALC_VIA_TEMPERTURE_SALINITY,BOLUS_VELOCITY_DIAGNOSTIC,REVERT_VERTICAL_RECON_AND_TRANSPOSED,&
    & INCLUDE_SLOPE_SQUARED_IMPLICIT

  USE mo_util_dbg_prnt,             ONLY: dbg_print
  USE mo_parallel_config,           ONLY: nproma
  USE mo_dynamics_config,           ONLY: nold, nnew
  USE mo_run_config,                ONLY: dtime, ltimer
  USE mo_ocean_types,               ONLY: t_hydro_ocean_state !, v_base
  USE mo_ocean_tracer_transport_types,        ONLY: t_ocean_tracer
  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
  USE mo_exception,                 ONLY: finish, message !, message_text, message
  USE mo_ocean_boundcond,           ONLY: top_bound_cond_tracer
  USE mo_ocean_physics_types,       ONLY: t_ho_params 
  USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff
  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
  USE mo_sync,                      ONLY: sync_patch_array_mult, sync_c, sync_e, sync_patch_array
  USE mo_timer,                     ONLY: timer_start, timer_stop
  USE mo_statistics,                ONLY: global_minmaxmean
  USE mo_mpi,                       ONLY: my_process_is_stdio !global_mpi_barrier
 
  USE mo_ocean_math_operators,      ONLY: grad_fd_norm_oce_3d_onBlock, verticalDeriv_scalar_onHalfLevels_on_block,&
    & verticalDeriv_vec_midlevel_on_block,div_oce_3d_ontriangles_onblock,div_oce_3d,verticaldiv_scalar_onfulllevels
  USE mo_scalar_product,            ONLY: map_cell2edges_3d,map_edges2cell_3d, &
    & map_scalar_center2prismtop, map_scalar_prismtop2center,map_edges2cell_with_height_3d,&
    & map_vec_prismtop2center_on_block,&
    & map_scalar_center2prismtop_GM, map_scalar_prismtop2center_GM,&
    & map_vec_prismtop2center_on_block_GM
  USE mo_ocean_thermodyn,           ONLY : calc_neutralslope_coeff_func_onColumn,&
                                         & calc_neutralslope_coeff_func_onColumn_UNESCO  
  USE mo_ocean_tracer_diffusion,    ONLY: tracer_diffusion_vertical_implicit    
  USE mo_ocean_thermodyn,           ONLY: calculate_density_onColumn ! by_Oliver

  IMPLICIT NONE
  
  PRIVATE
  CHARACTER(LEN=*), PARAMETER :: this_mod_name = 'GM_Redi'
  CHARACTER(LEN=16)           :: str_module = 'GM_Redi'  ! Output of module for 1 line debug
  INTEGER                     :: idt_src    = 1          ! Level of detail for 1 line debug

  PUBLIC  :: init_GMRedi
  PUBLIC  :: prepare_GMRedi
  PUBLIC  :: calc_ocean_physics
  
  PRIVATE :: calc_combined_GentMcWilliamsRedi_flux
  PRIVATE :: calc_neutral_slopes
  PRIVATE :: calc_tapering_function
  PRIVATE :: calc_entries_mixing_tensor
  PUBLIC  :: diagnose_Redi_flux_balance
  
  PUBLIC  :: prepare_ocean_physics_zstar
  PUBLIC  :: calc_ocean_physics_zstar
  PRIVATE :: calc_combined_GentMcWilliamsRedi_flux_zstar

  
!  PRIVATE :: calc_bolus_velocity
  PRIVATE :: vertical_GM

  TYPE(t_cartesian_coordinates), POINTER   :: taper_off_diagonal_vert(:,:,:)  
  TYPE(t_cartesian_coordinates), POINTER   :: taper_off_diagonal_horz(:,:,:)
  REAL(wp), POINTER                        :: taper_diagonal_horz(:,:,:)
  REAL(wp), POINTER                        :: taper_diagonal_vert_expl(:,:,:)    
  REAL(wp), POINTER                        :: taper_diagonal_vert_impl(:,:,:)            
 
CONTAINS

  !-------------------------------------------------------------------------
  SUBROUTINE init_GMRedi(patch_3d)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d

    ALLOCATE(taper_off_diagonal_vert(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks))   
    ALLOCATE(taper_off_diagonal_horz(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks))
    ALLOCATE(taper_diagonal_horz(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks))
    ALLOCATE(taper_diagonal_vert_expl(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks))    
    ALLOCATE(taper_diagonal_vert_impl(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks))           
    
  END SUBROUTINE init_GMRedi
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates the fluxes of the isoycnical diffusion following Redi.
  !!
  !!
!<Optimize:inUse:done>
  SUBROUTINE prepare_GMRedi(patch_3d, ocean_state, param, op_coeff)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_ho_params),                 INTENT(inout) :: param
    TYPE(t_operator_coeff),            INTENT(in) :: op_coeff
   !-------------------------------------------------------------------------------

    CALL calc_neutral_slopes(patch_3d, ocean_state, param, op_coeff)

    CALL calc_tapering_function(patch_3d, param, ocean_state)

    IF(TEST_MODE_REDI_ONLY.OR.TEST_MODE_GM_ONLY)THEN
      CALL finish(TRIM('mo_ocean_GM_Redi'), 'This GMRedi_configuration is not supported')
    ELSE
      CALL calc_entries_mixing_tensor(patch_3d, ocean_state, param, op_coeff)
    ENDIF     
     
  END SUBROUTINE prepare_GMRedi
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE calculates the fluxes of the isoycnical diffusion following Redi.
  !!
  !!
  !<Optimize:inUse:done>
  SUBROUTINE calc_ocean_physics(patch_3d, ocean_state, param, op_coeff, GMRedi_flux_horz, GMRedi_flux_vert, &
    & tracer, tracer_index, typeOfTracers)
    TYPE(t_patch_3d ),TARGET, INTENT(in)     :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET        :: ocean_state
    TYPE(t_ho_params), INTENT(inout)         :: param
    TYPE(t_operator_coeff),INTENT(in)        :: op_coeff
    REAL(wp), INTENT(inout)                  :: GMRedi_flux_horz(:,:,:)
    REAL(wp), INTENT(inout)                  :: GMRedi_flux_vert(:,:,:)
    TYPE(t_ocean_tracer), TARGET             :: tracer
    INTEGER, INTENT(IN)                      :: tracer_index
    CHARACTER(LEN=*), INTENT(in)             :: typeOfTracers
    
    !Local variables
    INTEGER :: start_cell_index, end_cell_index, cell_index, blockNo,level
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    
    !-------------------------------------------------------------------------------
    patch_2D         => patch_3d%p_patch_2d(1)
    cells_in_domain  => patch_2D%cells%in_domain   
!     GMRedi_flux_horz => ocean_state%p_diag%GMRedi_flux_horz(:,:,:,tracer_index)
!     GMRedi_flux_vert => ocean_state%p_diag%GMRedi_flux_vert(:,:,:,tracer_index)
        
    SELECT CASE(GMRedi_configuration)!GMRedi_configuration==Cartesian_Mixing)RETURN
 
    CASE(GMRedi_combined)
       
        CALL calc_combined_GentMcWilliamsRedi_flux( patch_3d,        &
            & ocean_state,      &
            & param,            &
            & op_coeff,         &
            & GMRedi_flux_horz, &
            & GMRedi_flux_vert, &
            & tracer,           &
            & tracer_index,     &
            & typeOfTracers)

    CASE DEFAULT
      CALL finish(TRIM('mo_ocean_GM_Redi'), 'This GMRedi_configuration is not supported')
    
    END SELECT
    
    IF(SWITCH_ON_REDI_BALANCE_DIAGONSTIC.AND.SLOPE_CALC_VIA_TEMPERTURE_SALINITY)THEN

     CALL diagnose_Redi_flux_balance(patch_3d, ocean_state, param, op_coeff,&
     & tracer_index, typeOfTracers)
    ENDIF           
    
  END SUBROUTINE calc_ocean_physics
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE calculates the fluxes of the isoycnical diffusion following Redi and the eddy flux of GM.
  !!
  !!
  !<Optimize:inUse>
  SUBROUTINE calc_combined_GentMcWilliamsRedi_flux(patch_3d, ocean_state, param, op_coeff,&
    & GMRedi_flux_horz, GMRedi_flux_vert, tracer, tracer_index, typeOfTracers)
    TYPE(t_patch_3d ),TARGET, INTENT(in)     :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET        :: ocean_state
    TYPE(t_ho_params),      INTENT(inout)    :: param
    TYPE(t_operator_coeff), INTENT(in)       :: op_coeff
    REAL(wp), INTENT(inout)                  :: GMRedi_flux_horz(:,:,:)
    REAL(wp), INTENT(inout)                  :: GMRedi_flux_vert(:,:,:)
    TYPE(t_ocean_tracer), TARGET             :: tracer
    INTEGER, INTENT(IN)                      :: tracer_index 
    CHARACTER(LEN=*), INTENT(in)             :: typeOfTracers
    
    !Local variables
    INTEGER :: start_cell_index, end_cell_index, cell_index,level,start_level,end_level,blockNo
    INTEGER :: start_edge_index, end_edge_index, edge_index   
    INTEGER :: il_c1, ib_c1, il_c2, ib_c2  
    TYPE(t_subset_range), POINTER :: all_cells, cells_in_domain, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp) :: flux_vert_center(nproma, n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates) :: flux_vec_horz_center(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates),POINTER :: tracer_gradient_horz_vec_center(:,:,:), slopes(:,:,:)
    REAL(wp), POINTER :: tracer_gradient_vert_center(:,:,:)

    REAL(wp) :: mapped_vertical_diagonal_impl(nproma,n_zlev+1,patch_3D%p_patch_2D(1)%alloc_cell_blocks)      
    
    REAL(wp) :: nabla_T_horz(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp) :: nabla_S_horz(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_e)   
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain 
    all_cells       => patch_2D%cells%all
    edges_in_domain => patch_2D%edges%in_domain 
    slopes          => ocean_state%p_aux%slopes 

    start_level=1
      
    IF (typeOfTracers == "ocean" )THEN
    !Set pointers for tracer gradients, according to actual tracer index
      IF(tracer_index==1)THEN
	!write(0,*) "DerivTemperature_vec_center"
	tracer_gradient_horz_vec_center => ocean_state%p_aux%PgradTemperature_horz_center
	tracer_gradient_vert_center     => ocean_state%p_aux%DerivTemperature_vert_center
      ELSEIF(tracer_index==2)THEN
	!write(0,*) "DerivSalinity_vec_center"
	tracer_gradient_horz_vec_center => ocean_state%p_aux%PgradSalinity_horz_center
	tracer_gradient_vert_center     => ocean_state%p_aux%DerivSalinity_vert_center
      ENDIF
    ELSE
!write(0,*)'-------------------------------------TRACER==3-----------------------------------------'    
!ocean_state%p_prog(nold(1))%tracer_collection%tracer(tracer_index)%concentration=ocean_state%p_diag%rho_GM
!&=ocean_state%p_prog(nold(1))%tracer_collection%tracer(tracer_index-1)%concentration 

      !Here we have to provide a sbr that calculates derivatives below    
      CALL calc_tracer_derivatives( patch_3d,&
                                  & tracer%concentration,&
                                  & ocean_state, &
                                  & op_coeff, &
                                  & tracer_index, &
                                  & typeOfTracers)

      tracer_gradient_horz_vec_center => ocean_state%p_aux%PgradTracer_horz_center
      tracer_gradient_vert_center     => ocean_state%p_aux%DerivTracer_vert_center
       
    ENDIF

!ICON_OMP_DO_PARALLEL  PRIVATE(start_cell_index,end_cell_index, cell_index, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_cells%start_block, all_cells%end_block        
      CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)
      DO cell_index = start_cell_index, end_cell_index
        DO level = 1, n_zlev
          flux_vec_horz_center(cell_index,level,blockNo)%x = 0.0_wp
        ENDDO
      ENDDO
    ENDDO
!ICON_OMP_END_DO_PARALLEL
      
    !The code below has to be executed for temperature, salinity and HAMMOC tracers, i.e. for all tracers.
    !
!ICON_OMP_DO_PARALLEL PRIVATE(start_cell_index,end_cell_index, cell_index, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block        
      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)      
      DO cell_index = start_cell_index, end_cell_index
 
        !GMRedi Flux at top layer: with tapering this is the just horizontal diffusion
        DO level = start_level, MIN(patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo),start_level)
          flux_vec_horz_center(cell_index,start_level,blockNo)%x &
          &=taper_diagonal_horz(cell_index,start_level,blockNo)  &
          &*tracer_gradient_horz_vec_center(cell_index,start_level,blockNo)%x

          !flux_vert_center(cell_index,start_level,blockNo) &
          !&=taper_diagonal_vert_expl(cell_index,start_level,blockNo)&
          !&*tracer_gradient_vert_center(cell_index,start_level,blockNo)

          ! the top level flux_vert_center will be filled from the second level, if it exists
          flux_vert_center(cell_index,start_level,blockNo) = 0.0_wp
          
        ENDDO

        DO level = start_level+1, patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)
          
          !horizontal GM-Redi Flux
          IF(ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)==0.0_wp)THEN
          tracer_gradient_horz_vec_center(cell_index,level,blockNo)%x(:)=0.0_wp
          ENDIF 
          flux_vec_horz_center(cell_index,level,blockNo)%x &
          &=taper_diagonal_horz(cell_index,level,blockNo)  &
          &*tracer_gradient_horz_vec_center(cell_index,level,blockNo)%x&
          !the second term vanishes if GM-Kappa=isoneutral diffusion
          &+taper_off_diagonal_horz(cell_index,level,blockNo)%x&
          &*tracer_gradient_vert_center(cell_index,level,blockNo)

!          flux_vec_horz_center(cell_index,level,blockNo)%x &
!          &=(taper_diagonal_horz(cell_index,level,blockNo)  &
!          &*tracer_gradient_horz_vec_center(cell_index,level,blockNo)%x&
!          &+taper_off_diagonal_horz(cell_index,level,blockNo)%x&
!          &*tracer_gradient_vert_center(cell_index,level,blockNo))

              
          !vertical GM-Redi Flux: this is the part that is explicit in time
          !If namelist option "switch_off_diagonal_vert_expl=TRUE" the
          !first (diagonal) term is set to zero. Default option "switch_off_diagonal_vert_expl=FALSE"
          !Second term contains isoneutral and GM components
          flux_vert_center(cell_index,level,blockNo)= &
          &taper_diagonal_vert_expl(cell_index,level,blockNo)&
          &*tracer_gradient_vert_center(cell_index,level,blockNo)&
          &+&
          &Dot_Product(tracer_gradient_horz_vec_center(cell_index,level,blockNo)%x,&
          &            taper_off_diagonal_vert(cell_index,level,blockNo)%x)

!          flux_vert_center(cell_index,level,blockNo)&
!          &= flux_vert_center(cell_index,level,blockNo)&
!          &*patch_3d%p_patch_1D(1)%prism_thick_flat_sfc_c(cell_index,level,blockNo)

        END DO

      END DO                
    END DO
!ICON_OMP_END_DO_PARALLEL

    !Map the explicit horizontal tracer flux from cell centers to edges (where the horizontal divergence is calculated)
    !
    ! use a vector communicator
    CALL sync_patch_array_mult(sync_c, patch_2D, 3, &
        & flux_vec_horz_center(:,:,:)%x(1), flux_vec_horz_center(:,:,:)%x(2), flux_vec_horz_center(:,:,:)%x(3))
    !    
    CALL map_cell2edges_3D( patch_3D,flux_vec_horz_center, GMRedi_flux_horz, op_coeff)

!ICON_OMP_DO_PARALLEL PRIVATE(start_edge_index,end_edge_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE
        DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block     
          CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)      
          DO edge_index = start_edge_index, end_edge_index
            DO level = start_level, patch_3D%p_patch_1D(1)%dolic_e(edge_index,blockNo)
            
              GMRedi_flux_horz(edge_index,level,blockNo)&
              &=GMRedi_flux_horz(edge_index,level,blockNo)*patch_3D%p_patch_1d(1)%prism_thick_e(edge_index,level,blockNo) 
              
            END DO                  
          END DO                
        END DO
!ICON_OMP_END_DO_PARALLEL
 
    
    IF(INCLUDE_SLOPE_SQUARED_IMPLICIT)THEN
    
    !Map the (explicit) vertical tracer flux to the prism top (where the vertical divergence is calculated later)
      IF(.NOT.REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN
	CALL map_scalar_center2prismtop(patch_3d, flux_vert_center, op_coeff,GMRedi_flux_vert)
      ELSEIF(REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN
	CALL map_scalar_center2prismtop_GM(patch_3d, flux_vert_center, op_coeff,GMRedi_flux_vert)    
      ENDIF

!       ! Now we treat the vertical isoneutral coefficient that is discretized implicitely in time.  
!       ! This is only neccessary once for temperature and salinity, the HAMOCC tracers use these value  
!       IF(tracer_index<=2)THEN
!         ! 
!         !1.) Interpolate the tapered coefficient for the vertical tracer flux from prism center to prism top: 
!         !this is the diagonal part that is handled implicitely in time
!       
!         IF(.NOT.REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN        
!           CALL map_scalar_center2prismtop( patch_3d, &
!           &                              taper_diagonal_vert_impl,&
!           &                              op_coeff,                &
!           &                              ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit)        
! 
!         ELSEIF(REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN
!           CALL map_scalar_center2prismtop_GM( patch_3d, &
!           &                              taper_diagonal_vert_impl,&
!           &                              op_coeff,                &
!           &                              ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit)        
!     
!         ENDIF
      IF(typeOfTracers == "ocean" ) THEN
      
	IF( tracer_index==1)THEN
	  CALL dbg_print('Old vert coeff: A_v', param%a_tracer_v(:,:,:, tracer_index),&
	  & this_mod_name, 4, patch_2D%cells%in_domain)
	  !CALL dbg_print('Old vert coeff: A_v', ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit,&
	  !& this_mod_name, 4, patch_2D%cells%in_domain)
	ENDIF
	  !
	  !2.) Here we combine the vertical GMRedicoefficient that is treated implicitely (mapped_vertical_diagonal_impl, this
	  !term involves the slopes-squared times the isopycnal mixing coefficient and is potentially large, therefore
	  !it is discretized implicitely) with the vertical mixing coefficient from the PP-scheme. 
	  !We follow the approach in POP, where these two contributions are added
	  !(see Reference manual POP, sect 5.1.3, in particular p. 41, after eq (150)).
        
!ICON_OMP_DO_PARALLEL PRIVATE(start_cell_index,end_cell_index, cell_index, level) ICON_OMP_DEFAULT_SCHEDULE
	DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block     
	  CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)      
	  DO cell_index = start_cell_index, end_cell_index
	    DO level = start_level, patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)
	      param%a_tracer_v(cell_index,level,blockNo, tracer_index) =  &
	      & param%a_tracer_v(cell_index,level,blockNo, tracer_index) + &
	      & ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit(cell_index,level,blockNo)
	    END DO                  
	  END DO                
	END DO
!ICON_OMP_END_DO_PARALLEL

	CALL sync_patch_array(sync_c, patch_2D, param%a_tracer_v(:,:,:,tracer_index))
           
        IF(tracer_index==1)THEN
          CALL dbg_print('New vert coeff: A_v', param%a_tracer_v(:,:,:, tracer_index),&
          & this_mod_name, 4, patch_2D%cells%in_domain)
        ENDIF
        
      ENDIF!IF(typeOfTracers == "ocean" ) 
      
    ELSEIF(.NOT.INCLUDE_SLOPE_SQUARED_IMPLICIT)THEN
        
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, level) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
    
	CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
    
	DO cell_index = start_cell_index, end_cell_index
	  DO level = start_level+1, patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)
	    !vertical GM-Redi Flux
	      GMRedi_flux_vert(cell_index,level,blockNo)   &
	      &=flux_vert_center(cell_index,level,blockNo) &
	      & + tracer_gradient_vert_center(cell_index,level,blockNo)&
	      & *taper_diagonal_vert_impl(cell_index,level,blockNo)

	  END DO                  
	END DO                
      END DO
!ICON_OMP_END_DO
            
    ENDIF  
    !---------DEBUG DIAGNOSTICS-------------------------------------------
!    idt_src=1  ! output print level (1-5, fix)
!    CALL dbg_print('InGMRedi: vert center',flux_vert_center,&
!    & str_module, idt_src, in_subset=cells_in_domain)
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('InGMRedi: GMRedi_vert',GMRedi_flux_vert,&
    & str_module, idt_src, in_subset=cells_in_domain)
    CALL dbg_print('InGMRedi: GMRedi_horz',GMRedi_flux_horz,&
    & str_module, idt_src, in_subset=edges_in_domain)
    !---------------------------------------------------------------------
    

  END SUBROUTINE calc_combined_GentMcWilliamsRedi_flux
  !-------------------------------------------------------------------------
 

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE calculates the slopes required for isopycnial diffusion and eddy parametrization.
  !!
  !!    !Note that in case of 1-component fluid we use the density for the slope calculation,
  !!    !all relevant imformation is stored in the tracer%temperature structure
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE calc_neutral_slopes(patch_3d, ocean_state, param, op_coeff)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_ho_params),                 INTENT(inout) :: param
    TYPE(t_operator_coeff),            INTENT(in) :: op_coeff
    
    !Local variables
    REAL(wp) :: grad_T_horz(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp) :: grad_S_horz(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp) :: grad_T_vert(nproma, n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: grad_S_vert(nproma, n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: grad_rho_GM_horz(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp) :: grad_rho_GM_vert(nproma, n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
   
    TYPE(t_cartesian_coordinates),POINTER :: grad_T_vec(:,:,:)    
    TYPE(t_cartesian_coordinates),POINTER :: grad_S_vec(:,:,:)        
    TYPE(t_cartesian_coordinates),POINTER :: grad_rho_GM_vec(:,:,:)            
    REAL(wp),POINTER :: grad_T_vert_center(:,:,:)
    REAL(wp),POINTER :: grad_S_vert_center(:,:,:)
    REAL(wp),POINTER :: grad_rho_GM_vert_center(:,:,:)
       
    INTEGER :: level, blockNo, je, start_level, end_level
    INTEGER :: start_cell_index, end_cell_index, cell_index
    INTEGER :: start_edge_index, end_edge_index
    REAL(wp) :: alpha, beta
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain, all_cells
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp), POINTER :: pot_temp(:,:,:), salinity(:,:,:),rho_GM(:,:,:)
!     REAL(wp):: size_grad_T_horz_vec(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
!     REAL(wp):: size_grad_S_horz_vec(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)

    REAL(wp), POINTER :: depth_cellinterface(:,:,:)
    REAL(wp):: neutral_coeff(1:n_zlev, 2), salinityColumn(1:n_zlev),neutral_coeff1(1:n_zlev, 2)
    !-----------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    all_cells       => patch_2D%cells%all
    cells_in_domain => patch_2D%cells%in_domain
    edges_in_domain => patch_2D%edges%in_domain
    depth_cellinterface => patch_3D%p_patch_1d(1)%depth_cellinterface
    !-------------------------------------------------------------------------

    IF(no_tracer>=2)THEN
      IF(SLOPE_CALC_VIA_TEMPERTURE_SALINITY)THEN
        pot_temp          => ocean_state%p_prog(nold(1))%tracer_collection%tracer(1)%concentration
        grad_T_vec        => ocean_state%p_aux%PgradTemperature_horz_center
        grad_T_vert_center=> ocean_state%p_aux%DerivTemperature_vert_center
    
        salinity          => ocean_state%p_prog(nold(1))%tracer_collection%tracer(2)%concentration
        grad_S_vec        => ocean_state%p_aux%PgradSalinity_horz_center
        grad_S_vert_center=> ocean_state%p_aux%DerivSalinity_vert_center
      
      ELSE
        rho_GM                 => ocean_state%p_diag%rho_GM   
        grad_rho_GM_vec        => ocean_state%p_aux%PgradDensity_horz_center
        grad_rho_GM_vert_center=> ocean_state%p_aux%DerivDensity_vert_center

        pot_temp          => ocean_state%p_prog(nold(1))%tracer_collection%tracer(1)%concentration
        grad_T_vec        => ocean_state%p_aux%PgradTemperature_horz_center
        grad_T_vert_center=> ocean_state%p_aux%DerivTemperature_vert_center
        
        salinity          => ocean_state%p_prog(nold(1))%tracer_collection%tracer(2)%concentration
        grad_S_vec        => ocean_state%p_aux%PgradSalinity_horz_center
        grad_S_vert_center=> ocean_state%p_aux%DerivSalinity_vert_center
        
      ENDIF
    !Note that in case of 1-component fluid we use the density for the slope calculation,
    !all relevant imformation is stored in the tracer%temperature structure
    ELSEIF(no_tracer==1)THEN      
      !pot_temp          => ocean_state%p_diag%rho	
      pot_temp          => ocean_state%p_prog(nold(1))%tracer_collection%tracer(1)%concentration      
      grad_T_vec        => ocean_state%p_aux%PgradTemperature_horz_center
      grad_T_vert_center=> ocean_state%p_aux%DerivTemperature_vert_center

      rho_GM                 => ocean_state%p_diag%rho_GM   
      grad_rho_GM_vec        => ocean_state%p_aux%PgradDensity_horz_center
      grad_rho_GM_vert_center=> ocean_state%p_aux%DerivDensity_vert_center
    
    ENDIF

!     CALL sync_patch_array(sync_c, patch_2D, rho_GM)
    
    !grad_T_vec(:,:,:)%x(1)=0.0_wp
    !grad_T_vec(:,:,:)%x(2)=0.0_wp
    !grad_T_vec(:,:,:)%x(3)=0.0_wp
    !grad_T_vert_center(:,:,:)=0.0_wp

    start_level = 1

    !-------------------------------------------------------------------------------
    !1) calculation of horizontal and vertical gradient for potential temperature and salinity
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_edge_index,end_edge_index) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      grad_T_horz(:,:,blockNo) = 0.0_wp
      !1a) calculate horizontal gradient of temperature
      CALL grad_fd_norm_oce_3d_onBlock ( &
        & pot_temp, &
        & patch_3D,                           &
        & op_coeff%grad_coeff(:,:,blockNo), &
        & grad_T_horz(:,:,blockNo),           &
        & start_edge_index, end_edge_index, blockNo)

     !1b)calculate horizontal  gradient of salinity
      IF(no_tracer>=2)THEN
        grad_S_horz(:,:,blockNo) = 0.0_wp
        CALL grad_fd_norm_oce_3d_onBlock ( &
          & salinity, &
          & patch_3D,                    &
          & op_coeff%grad_coeff(:,:,blockNo), &
          & grad_S_horz(:,:,blockNo),           &
          & start_edge_index, end_edge_index, blockNo)
     ENDIF
     IF(.NOT.SLOPE_CALC_VIA_TEMPERTURE_SALINITY)THEN
     
       grad_rho_GM_horz(:,:,blockNo) = 0.0_wp
       CALL grad_fd_norm_oce_3d_onBlock ( &
        & rho_GM, &
        & patch_3D,                    &
        & op_coeff%grad_coeff(:,:,blockNo), &
        & grad_rho_GM_horz(:,:,blockNo),           &
        & start_edge_index, end_edge_index, blockNo)
     ENDIF
    END DO ! blocks
!ICON_OMP_END_DO


    !---------------------------------------------------------------------
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index) ICON_OMP_DEFAULT_SCHEDULE
 
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

      grad_T_vert(:,:,blockNo) = 0.0_wp ! this is only for the top level
      !1c) calculation of vertical derivative for temperature and salinity
      CALL verticalDeriv_scalar_onHalfLevels_on_block( patch_3d,                &
                                                  & pot_temp(:,:,blockNo),   &
                                                  & grad_T_vert(:,:,blockNo),&
                                                  & start_level+1,             &
                                                  & blockNo,                 &
                                                  & start_cell_index,        &
                                                  & end_cell_index)
      IF((no_tracer>=2))THEN
        grad_S_vert(:,:,blockNo) = 0.0_wp! this is only for the top level
        CALL verticalDeriv_scalar_onHalfLevels_on_block( patch_3d,                &
                                                    & salinity(:,:,blockNo),   &
                                                    & grad_S_vert(:,:,blockNo),&
                                                    & start_level+1,             &
                                                    & blockNo,                 &
                                                    & start_cell_index,        &
                                                    & end_cell_index)
      ENDIF
      IF(.NOT.SLOPE_CALC_VIA_TEMPERTURE_SALINITY)THEN
        !IF((no_tracer>=2))THEN 
        !  grad_rho_GM_vert(:,:,blockNo) = ocean_state%p_diag%grad_rho_PP_vert(:,:,blockNo) 
        !ELSEIF((no_tracer==1))THEN                       
          grad_rho_GM_vert(:,:,blockNo) = 0.0_wp! this is only for the top level
          CALL verticalDeriv_scalar_onHalfLevels_on_block( patch_3d,                &
                                                    & rho_GM(:,:,blockNo),   &
                                                    & grad_rho_GM_vert(:,:,blockNo),&
                                                    & start_level+1,             &
                                                    & blockNo,                 &
                                                    & start_cell_index,        &
                                                    & end_cell_index)
        grad_rho_GM_vert(:,:,blockNo)=-grad_rho_GM_vert(:,:,blockNo)                                            
        !ENDIF
      ENDIF
    END DO ! blocks
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL

!     CALL sync_patch_array(sync_e, patch_2D, grad_rho_GM_horz)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    IF(SLOPE_CALC_VIA_TEMPERTURE_SALINITY)THEN    
      CALL dbg_print('calc_slopes: grad_T_horz',grad_T_horz,&
        & str_module,idt_src, in_subset=edges_in_domain)
      IF(no_tracer==2)THEN  
      CALL dbg_print('calc_slopes: grad_S_horz',grad_S_horz,&
        & str_module,idt_src, in_subset=edges_in_domain)
      ENDIF
    ELSEIF(.NOT.SLOPE_CALC_VIA_TEMPERTURE_SALINITY)THEN
    CALL dbg_print('calc_slopes: grad_rho_GM_horz',grad_rho_GM_horz,&
      & str_module,idt_src, in_subset=edges_in_domain)
    CALL dbg_print('calc_slopes: grad_vert',grad_rho_GM_vert,&
      & str_module,idt_src, in_subset=cells_in_domain)
      
        
    ENDIF    
!    CALL sync_patch_array(sync_c, patch_2D, grad_T_vert)       
!    IF(no_tracer>=2)   CALL sync_patch_array(sync_c, patch_2D, grad_S_vert)   

   !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=4  ! output print level (1-5, fix)
    IF(SLOPE_CALC_VIA_TEMPERTURE_SALINITY)THEN      
      CALL dbg_print('neutral_slopes: grad_T_vert',grad_T_vert,&
        & str_module,idt_src, in_subset=cells_in_domain)
      IF(no_tracer==2)THEN        
      CALL dbg_print('neutral_slopes: grad_S_vert',grad_S_vert,&
        & str_module,idt_src, in_subset=cells_in_domain) 
      ENDIF
    ELSEIF(.NOT.SLOPE_CALC_VIA_TEMPERTURE_SALINITY)THEN  
    CALL dbg_print('neutral_slopes: grad_rho_GM_vert',grad_rho_GM_vert,&
      & str_module,idt_src, in_subset=cells_in_domain) 
    
    ENDIF   
  !---------------------------------------------------------------------   
   
    !2) map horizontal and vertial derivative to cell centered vector
    CALL map_edges2cell_3d(patch_3D,  &
        & grad_T_horz,                &
        & op_coeff,                   &
        & grad_T_vec,                 &
        & subset_range=cells_in_domain)

!     CALL sync_patch_array(sync_c, patch_2D, grad_T_vert)
    IF(.NOT.REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN
      CALL map_scalar_prismtop2center(patch_3d,&
        & grad_T_vert,                &
        & op_coeff,                 &
        & grad_T_vert_center)
    ELSEIF(REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN
      CALL map_scalar_prismtop2center_GM(patch_3d,&
        & grad_T_vert,                &
        & op_coeff,                 &
        & grad_T_vert_center)    
    
    ENDIF    
        
       
    IF((no_tracer>=2))THEN        
      
      CALL map_edges2cell_3d(patch_3D,  &      
          & grad_S_horz,                &
          & op_coeff,                   &
          & grad_S_vec,                 &
          & subset_range=cells_in_domain)
      IF(.NOT.REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN          
        CALL map_scalar_prismtop2center(patch_3d,&
          & grad_S_vert,                       &
          & op_coeff,                          &
          & grad_S_vert_center)
      ELSEIF(REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN      
        CALL map_scalar_prismtop2center_GM(patch_3d,&
          & grad_S_vert,                       &
          & op_coeff,                          &
          & grad_S_vert_center)  
      ENDIF              
    ENDIF
    IF(.NOT.SLOPE_CALC_VIA_TEMPERTURE_SALINITY)THEN

      CALL map_edges2cell_3d(patch_3D,  &      
          & grad_rho_GM_horz,           &
          & op_coeff,                   &
          & grad_rho_GM_vec,            &
          & subset_range=cells_in_domain)
    IF(.NOT.REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN                    
      CALL map_scalar_prismtop2center(patch_3d,&
          & grad_rho_GM_vert,                  &
          & op_coeff,                          &
          & grad_rho_GM_vert_center)
     ELSEIF(REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN  
      CALL map_scalar_prismtop2center_GM(patch_3d,&
          & grad_rho_GM_vert,                  &
          & op_coeff,                          &
          & grad_rho_GM_vert_center)     
     
     ENDIF         
    
    ENDIF    
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------    
    IF(no_tracer>=2)THEN

        IF(SLOPE_CALC_VIA_TEMPERTURE_SALINITY)THEN
!ICON_OMP_PARALLEL PRIVATE(salinityColumn)                             
    salinityColumn(1:n_zlev) = sal_ref  ! in case of absent salinty tracer!
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,neutral_coeff, &
!ICON_OMP  level) ICON_OMP_DEFAULT_SCHEDULE
          DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
            CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

            DO cell_index = start_cell_index, end_cell_index
              end_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
              IF(end_level <= min_dolic) CYCLE
                salinityColumn(1:end_level) = salinity(cell_index,1:end_level,blockNo)

                !4.1) calculate slope coefficients as thermal expansion and saline contraction coefficients
                !
                !Nonlinear EOS, slope coefficients are calculated via the McDougall-method
                IF(EOS_TYPE/=1)THEN
          
                  neutral_coeff = calc_neutralslope_coeff_func_onColumn(         &
                  & pot_temp(cell_index,1:end_level,blockNo), salinityColumn(1:end_level),  &
                  & depth_cellinterface(cell_index,2:end_level+1,blockNo), end_level)
                  
                  !Thisis just for testing and can be deleted at some point                 
                  !neutral_coeff =calc_neutralslope_coeff_func_onColumn_UNESCO(         &
                  !& pot_temp(cell_index,1:end_level,blockNo), salinityColumn(1:end_level),  &
                  !& depth_cellinterface(cell_index,2:end_level+1,blockNo), end_level)                  
                  !neutral_coeff(:,1) = LinearThermoExpansionCoefficient 
                  ! neutral_coeff(:,2) = LinearHalineContractionCoefficient
                   
                !Linear EOS: slope coefficients are equal to EOS-coefficients
                ELSEIF(EOS_TYPE==1)THEN          
                  neutral_coeff(:,1) = LinearThermoExpansionCoefficient 
                  neutral_coeff(:,2) = LinearHalineContractionCoefficient
                ENDIF
        
                DO level = start_level+1, end_level-1

                  ocean_state%p_aux%slopes(cell_index,level,blockNo)%x &
                  & = -neutral_coeff(level,1) * grad_T_vec(cell_index,level,blockNo)%x + &
                  &      neutral_coeff(level,2) * grad_S_vec(cell_index,level,blockNo)%x

                  ocean_state%p_aux%slopes_drdx(cell_index,level,blockNo)=&
                  & DOT_PRODUCT(ocean_state%p_aux%slopes(cell_index,level,blockNo)%x,&
                           &ocean_state%p_aux%slopes(cell_index,level,blockNo)%x)
                  ocean_state%p_aux%slopes_drdx(cell_index,level,blockNo)=&
                  & sqrt(ocean_state%p_aux%slopes_drdx(cell_index,level,blockNo))

                  ocean_state%p_aux%slopes_drdz(cell_index,level,blockNo)=&
                  & -(-neutral_coeff(level,1) * grad_T_vert_center(cell_index,level,blockNo)+ &
                  &      neutral_coeff(level,2) * grad_S_vert_center(cell_index,level,blockNo))
              
                  ocean_state%p_aux%slopes(cell_index,level,blockNo)%x &
                    & = -(-neutral_coeff(level,1) * grad_T_vec(cell_index,level,blockNo)%x + &
                    &      neutral_coeff(level,2) * grad_S_vec(cell_index,level,blockNo)%x)  &           
                    &   /(-neutral_coeff(level,1) * grad_T_vert_center(cell_index,level,blockNo)+ &
                    &      neutral_coeff(level,2) * grad_S_vert_center(cell_index,level,blockNo)-dbl_eps)

                  ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)=&
                    & DOT_PRODUCT(ocean_state%p_aux%slopes(cell_index,level,blockNo)%x,&
                                 &ocean_state%p_aux%slopes(cell_index,level,blockNo)%x)
                END DO
         !Perform nearest neighbor interpolation at level where slopes are not well-defined
!           ocean_state%p_aux%slopes(cell_index,start_level,blockNo)%x &
!           &= ocean_state%p_aux%slopes(cell_index,start_level+1,blockNo)%x   
!           ocean_state%p_aux%slopes(cell_index,end_level,blockNo)%x &
!           &= ocean_state%p_aux%slopes(cell_index,end_level-1,blockNo)%x        
            END DO ! cell_index = start_cell_index, end_cell_index
          END DO  ! blockNo = all_cells%start_block, all_cells%end_block
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL
    
        ELSEIF(.NOT.SLOPE_CALC_VIA_TEMPERTURE_SALINITY)THEN
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level, &
!ICON_OMP  level) ICON_OMP_DEFAULT_SCHEDULE
          DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
            CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

            DO cell_index = start_cell_index, end_cell_index
              end_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
              IF(end_level <= min_dolic) CYCLE
              !salinityColumn(1:end_level) = salinity(cell_index,1:end_level,blockNo)
                
              DO level = start_level+1, end_level-1

                ocean_state%p_aux%slopes(cell_index,level,blockNo)%x                      &
                & =   grad_rho_GM_vec(cell_index,level,blockNo)%x  &
                &    /(grad_rho_GM_vert_center(cell_index,level,blockNo)-dbl_eps)

                ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)=&
                & DOT_PRODUCT(ocean_state%p_aux%slopes(cell_index,level,blockNo)%x,&
                             &ocean_state%p_aux%slopes(cell_index,level,blockNo)%x)

                ocean_state%p_aux%slopes_drdx(cell_index,level,blockNo)=&
                & sqrt(DOT_PRODUCT(grad_rho_GM_vec(cell_index,level,blockNo)%x,&
                                 & grad_rho_GM_vec(cell_index,level,blockNo)%x))

                ocean_state%p_aux%slopes_drdz(cell_index,level,blockNo)=&
                &         grad_rho_GM_vert_center(cell_index,level,blockNo)            

              END DO                         
!Perform nearest neighbor interpolation at level where slopes are not well-defined
!           ocean_state%p_aux%slopes(cell_index,start_level,blockNo)%x &
!           &= ocean_state%p_aux%slopes(cell_index,start_level+1,blockNo)%x   
!           ocean_state%p_aux%slopes(cell_index,end_level,blockNo)%x &
!           &= ocean_state%p_aux%slopes(cell_index,end_level-1,blockNo)%x                    
            END DO ! cell_index = start_cell_index, end_cell_index
          END DO  ! blockNo = all_cells%start_block, all_cells%end_block
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL
    
        ENDIF
        
    ELSEIF(no_tracer==1)THEN!Case of single tracer
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level, &
!ICON_OMP  level) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

        DO cell_index = start_cell_index, end_cell_index
          end_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
          IF(end_level <= min_dolic) CYCLE

          DO level = start_level+1, end_level-1
 
                ocean_state%p_aux%slopes(cell_index,level,blockNo)%x                      &
                & =   grad_rho_GM_vec(cell_index,level,blockNo)%x  &
                &    /(grad_rho_GM_vert_center(cell_index,level,blockNo)-dbl_eps)

                ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)=&
                & DOT_PRODUCT(ocean_state%p_aux%slopes(cell_index,level,blockNo)%x,&
                             &ocean_state%p_aux%slopes(cell_index,level,blockNo)%x)

                ocean_state%p_aux%slopes_drdx(cell_index,level,blockNo)=&
                & sqrt(DOT_PRODUCT(grad_rho_GM_vec(cell_index,level,blockNo)%x,&
                                 & grad_rho_GM_vec(cell_index,level,blockNo)%x))

                ocean_state%p_aux%slopes_drdz(cell_index,level,blockNo)=&
                &         grad_rho_GM_vert_center(cell_index,level,blockNo)            
          END DO
          !Perform nearest neighbor interpolation at level where slopes are not well-defined
!           ocean_state%p_aux%slopes(cell_index,start_level,blockNo)%x &
!           &= ocean_state%p_aux%slopes(cell_index,start_level+1,blockNo)%x   
!           ocean_state%p_aux%slopes(cell_index,end_level,blockNo)%x &
!           &= ocean_state%p_aux%slopes(cell_index,end_level-1,blockNo)%x    
        END DO ! cell_index = start_cell_index, end_cell_index
      END DO  ! blockNo = all_cells%start_block, all_cells%end_block
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL  
    ENDIF


 !---------DEBUG DIAGNOSTICS-------------------------------------------
!     idt_src=3  ! output print level (1-5, fix)
!     CALL dbg_print('calc_slopes: sizegradT1',grad_T_vec(:,:,:)%x(1), &
!       &  str_module,idt_src, in_subset=cells_in_domain)      
!     CALL dbg_print('calc_slopes: sizegradT2',grad_T_vec(:,:,:)%x(2), &
!       &  str_module,idt_src, in_subset=cells_in_domain)      
!     CALL dbg_print('calc_slopes: sizegradT3',grad_T_vec(:,:,:)%x(3), &
!       &  str_module,idt_src, in_subset=cells_in_domain)      
!      
!      
!    IF(no_tracer>=2)  &
!      & CALL dbg_print('calc:slopes: sizegradS',size_grad_S_horz_vec,&
!      &                str_module,idt_src, in_subset=cells_in_domain)
  !---------------------------------------------------------------------   

  !---------------------------------------------------------------------
  !---------DEBUG DIAGNOSTICS-------------------------------------------
  idt_src=3  ! output print level (1-5, fix)
  CALL dbg_print('calc_slopes: squared',ocean_state%p_aux%slopes_squared,&
    & str_module,idt_src, in_subset=cells_in_domain)
  CALL dbg_print('calc_slopes: slope abs',sqrt(ocean_state%p_aux%slopes_squared(:,:,:)),&
    & str_module,idt_src, in_subset=cells_in_domain)
  CALL dbg_print('calc_slopes: slopes dz',ocean_state%p_aux%slopes_drdz,&
    & str_module,idt_src, in_subset=cells_in_domain)
  CALL dbg_print('calc_slopes: slopes dx',ocean_state%p_aux%slopes_drdx,&
    & str_module,idt_src, in_subset=cells_in_domain)
  
   !---------------------------------------------------------------------
 !  DO level= 1, n_zlev  
 !  write(0,*)'max-min vert deriv',level,maxval(grad_T_vert_center(:,level,:)),&
 !  & minval(grad_T_vert_center(:,level,:))!,maxval(grad_S_vert_center(:,level,:)),&
 !  !& minval(grad_S_vert_center(:,level,:))
 !  END DO
  !--------------------------------------------------------------------- 

  IF(BOLUS_VELOCITY_DIAGNOSTIC)THEN
    write(0,*)'calc_bolus_velocity is enabled in the namelist but deactivated in the code'
!    CALL calc_bolus_velocity(patch_3d, ocean_state, param, op_coeff)
  ENDIF
  !---------------------------------------------------------------------  
  
 
  END SUBROUTINE calc_neutral_slopes
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE calculates horizontal and vertical derivatives of tracer and the corresponding
  !!    !reconstructions on cell centers.
  !!
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE calc_tracer_derivatives(patch_3d, tracer, ocean_state, op_coeff, tracer_index, typeOfTracers)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    REAL(wp)                                         :: tracer(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_operator_coeff),            INTENT(in)    :: op_coeff
    INTEGER, INTENT(in)                              :: tracer_index
    CHARACTER(LEN=*), INTENT(in)                     :: typeOfTracers
     
    !Local variables
    TYPE(t_cartesian_coordinates),POINTER :: grad_horz_center_vec(:,:,:)    
    REAL(wp),POINTER :: grad_vert_center(:,:,:)
    REAL(wp),POINTER :: grad_vert(:,:,:)    
    REAL(wp),POINTER :: grad_horz(:,:,:)
               
    INTEGER :: level, blockNo, je, start_level, end_level
    INTEGER :: start_cell_index, end_cell_index, cell_index
    INTEGER :: start_edge_index, end_edge_index
    TYPE(t_subset_range), POINTER :: cells_in_domain, edges_in_domain, all_cells
    TYPE(t_patch), POINTER :: patch_2D
    !-----------------------------------------------------------------------
    patch_2D        => patch_3D%p_patch_2D(1)
    all_cells       => patch_2D%cells%all
    cells_in_domain => patch_2D%cells%in_domain
    edges_in_domain => patch_2D%edges%in_domain
    start_level=1
    !-------------------------------------------------------------------------

    IF (typeOfTracers == "ocean" )THEN

      IF(tracer_index==1)THEN
      
	grad_horz            => ocean_state%p_aux%temperature_grad_horz
	grad_vert            => ocean_state%p_aux%temperature_deriv_vert        

	grad_horz_center_vec => ocean_state%p_aux%PgradTemperature_horz_center
	grad_vert_center     => ocean_state%p_aux%DerivTemperature_vert_center
	    
      ELSEIF(tracer_index==2)THEN      
      
	grad_horz             => ocean_state%p_aux%salinity_grad_horz
	grad_vert             => ocean_state%p_aux%salinity_deriv_vert  
	    
	grad_horz_center_vec  => ocean_state%p_aux%PgradSalinity_horz_center
	grad_vert_center      => ocean_state%p_aux%DerivSalinity_vert_center
      ENDIF
      
    ELSE          
      grad_horz            => ocean_state%p_aux%tracer_grad_horz
      grad_vert            => ocean_state%p_aux%tracer_deriv_vert
      
      grad_horz_center_vec => ocean_state%p_aux%PgradTracer_horz_center
      grad_vert_center     => ocean_state%p_aux%DerivTracer_vert_center    
    ENDIF

    !-------------------------------------------------------------------------------
    !1) calculation of horizontal and vertical gradient for potential temperature and salinity
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_edge_index,end_edge_index) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)

      grad_horz(:,:,blockNo) = 0.0_wp
      !1a) calculate horizontal gradient
      CALL grad_fd_norm_oce_3d_onBlock ( &
        & tracer, &
        & patch_3D,                           &
        & op_coeff%grad_coeff(:,:,blockNo), &
        & grad_horz(:,:,blockNo),           &
        & start_edge_index, end_edge_index, blockNo)

    END DO ! blocks
!ICON_OMP_END_DO


    !---------------------------------------------------------------------
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index) ICON_OMP_DEFAULT_SCHEDULE
 
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

      grad_vert(:,:,blockNo) = 0.0_wp ! this is only for the top level

      !1c) calculation of vertical derivative
      CALL verticalDeriv_scalar_onHalfLevels_on_block( patch_3d,                &
                                                  & tracer(:,:,blockNo),   &
                                                  & grad_vert(:,:,blockNo),&
                                                  & start_level+1,             &
                                                  & blockNo,                 &
                                                  & start_cell_index,        &
                                                  & end_cell_index)
    END DO ! blocks
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL


    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('calc_derivatives: grad_horz',grad_horz,&
      & str_module,idt_src, in_subset=edges_in_domain)

    CALL dbg_print('calc_derivatives: grad_vert',grad_vert,&
      & str_module,idt_src, in_subset=cells_in_domain)
  !---------------------------------------------------------------------   
   
    !2) map horizontal and vertial derivative to cell centered vector
    CALL map_edges2cell_3d(patch_3D,  &
        & grad_horz,                  &
        & op_coeff,                   &
        & grad_horz_center_vec,       &
        & subset_range=cells_in_domain)


    IF(.NOT.REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN
    CALL map_scalar_prismtop2center(patch_3d,&
        & grad_vert,                &
        & op_coeff,                 &
        & grad_vert_center)
    ELSEIF(REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN
    CALL map_scalar_prismtop2center_GM(patch_3d,&
        & grad_vert,                &
        & op_coeff,                 &
        & grad_vert_center)    
    
    ENDIF    
    !------------------------------------------------------------------------------
    CALL dbg_print('calc_derivatives: grad_vert center',grad_vert_center,&
      & str_module,idt_src, in_subset=cells_in_domain)

  END SUBROUTINE calc_tracer_derivatives
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates the tapering functions used in Gent-McWilliams-Redi parametrization.
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE calc_tapering_function(patch_3d, param, ocean_state)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_ho_params),                 INTENT(inout) :: param
    !TYPE(t_operator_coeff),            INTENT(in)    :: op_coeff
    
    !Local variables
    INTEGER :: start_cell_index, end_cell_index, cell_index,level,start_level,end_level, blockNo
    INTEGER :: start_edge_index, end_edge_index, je     
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp) :: lambda
    REAL(wp) :: depth_scale, depth
    REAL(wp) :: Coriolis_abs
    REAL(wp) :: inv_S_d, slope_abs, inv_cell_characteristic_length, cell_max_slope, cell_critical_slope
    REAL(wp) :: cell_characteristic_length
    
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain 
    start_level=1
    inv_S_d = 1.0_wp / S_d
       
    !-------------------------------------------------------------------------------
!ICON_OMP_PARALLEL   
    IF (GMRedi_usesRelativeMaxSlopes) THEN
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level, slope_abs, inv_cell_characteristic_length, &
!ICON_OMP cell_max_slope, cell_critical_slope) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

        DO cell_index = start_cell_index, end_cell_index
          end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)
          inv_cell_characteristic_length = 1.0_wp / SQRT(patch_2D%cells%area(cell_index,blockNo))
          cell_characteristic_length     = SQRT(patch_2D%cells%area(cell_index,blockNo))
          
          IF(end_level >= min_dolic) THEN

            DO level = start_level, end_level

              !cell_critical_slope      = S_critical  &
              !  & * patch_3d%p_patch_1d(1)%prism_thick_c(cell_index,level,blockNo) &
              !  & * inv_cell_characteristic_length

              !cell_max_slope = S_max &
              !  & * patch_3d%p_patch_1d(1)%prism_thick_c(cell_index,level,blockNo) &
              !  & * cell_characteristic_length&
              !  &/(2.0_wp*param%k_tracer_isoneutral(cell_index,level,blockNo)*dtime)

               ! FIXME: new from Peter
               cell_critical_slope      = S_critical  &
                 & * patch_3d%p_patch_1d(1)%prism_thick_c(cell_index,level,blockNo)**2 &
                 & /(param%k_tracer_isoneutral(cell_index,level,blockNo)*dtime)


              cell_max_slope      = S_max  &
                & * patch_3d%p_patch_1d(1)%prism_thick_c(cell_index,level,blockNo) &
                & * inv_cell_characteristic_length


              slope_abs    = sqrt(ocean_state%p_aux%slopes_squared(cell_index,level,blockNo))           

              !IF(slope_abs <= cell_max_slope)THEN

                ocean_state%p_aux%taper_function_1(cell_index,level,blockNo) &
                  &= 0.5_wp*(1.0_wp + tanh((cell_max_slope - slope_abs)*inv_S_d))
                
                !ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)=tanh((cell_critical_slope - slope_abs)*inv_S_d)!&
                !&=  patch_3d%p_patch_1d(1)%prism_thick_c(cell_index,level,blockNo) &
                !& * inv_cell_characteristic_length  

              !ELSE
              !  ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)=0.0_wp
              !ENDIF
            END DO
          ENDIF
        END DO
      END DO

!ICON_OMP_END_DO
    ELSE! not GMRedi_usesRelativeMaxSlopes
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level, slope_abs, inv_cell_characteristic_length, &
!ICON_OMP cell_max_slope, cell_critical_slope) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        cell_max_slope      = S_max
        cell_critical_slope = S_critical
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

        DO cell_index = start_cell_index, end_cell_index
          end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)
          inv_cell_characteristic_length = 1.0_wp / SQRT(patch_2D%cells%area(cell_index,blockNo))

          IF(end_level >= min_dolic) THEN

            DO level = start_level, end_level
            
             slope_abs    = sqrt(ocean_state%p_aux%slopes_squared(cell_index,level,blockNo))           
              
              !IF(slope_abs <= cell_max_slope)THEN
                ocean_state%p_aux%taper_function_1(cell_index,level,blockNo) &
                  &= 0.5_wp*(1.0_wp + tanh((cell_max_slope - slope_abs)*inv_S_d))
              !ELSE
              !  ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)=0.0_wp
              !ENDIF

            END DO
          ENDIF
        END DO
      END DO
!ICON_OMP_END_DO
    ENDIF !GMRedi_usesRelativeMaxSlopes

    !tapering schemes other than Danabasoglu-McWilliams require a second
    !tapering function
    IF(tapering_scheme/=tapering_DanaMcWilliams)THEN

!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level, &
!ICON_OMP Coriolis_abs,lambda,depth_scale,depth) ICON_OMP_DEFAULT_SCHEDULE

      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

        DO cell_index = start_cell_index, end_cell_index
          end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

          IF(end_level >= min_dolic) THEN

            DO level = start_level, end_level-1
              Coriolis_abs=abs(patch_2d%cells%f_c(cell_index,blockNo))
              lambda      = min(max(15000._WP, c_speed/Coriolis_abs),100000.0_WP)
              depth_scale = lambda*sqrt(ocean_state%p_aux%slopes_squared(cell_index,level,blockNo))
              depth       = patch_3d%p_patch_1d(1)%depth_CellMiddle(cell_index,level,blockNo)

              IF(depth<=depth_scale)THEN
                ocean_state%p_aux%taper_function_2(cell_index,level,blockNo) &
                  & = 0.5_wp*(1.0_wp+sin(pi*depth/depth_scale-pi/2.0_wp))
              ELSE
                ocean_state%p_aux%taper_function_2(cell_index,level,blockNo) =1.0_wp
              ENDIF

            END DO
          ENDIF
        END DO
      END DO
!ICON_OMP_END_DO
!       CALL sync_patch_array(sync_c, patch_2D,ocean_state%p_aux%taper_function_2)
    ENDIF
!ICON_OMP_END_PARALLEL

!  Do level=1,n_zlev
!  write(0,*)'max-min taper1/2',level,&
!   &maxval( ocean_state%p_aux%taper_function_1(:,level,:)),&
!   &minval( ocean_state%p_aux%taper_function_1(:,level,:)),& 
!  &maxval( ocean_state%p_aux%taper_function_2(:,level,:)),&
!  &minval( ocean_state%p_aux%taper_function_2(:,level,:)),& 
!  &maxval( ocean_state%p_aux%slopes_squared(:,level,:)),&
! &minval( ocean_state%p_aux%slopes_squared(:,level,:)) 
!  End do        
!  stop
     CALL dbg_print('calc_tapering_fct1:', ocean_state%p_aux%taper_function_1 ,&
       & this_mod_name, 3, patch_2D%cells%in_domain)
     IF(tapering_scheme/=tapering_DanaMcWilliams) THEN
       CALL dbg_print('calc_tapering_fct2:', ocean_state%p_aux%taper_function_2 ,&
       & this_mod_name, 3, patch_2D%cells%in_domain)
     ENDIF
   
  END SUBROUTINE calc_tapering_function
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates the entries of he Gent-McWilliams mixing tensor, 
  !! !             this includes the tapering for which different options are available
  !! !
  !!         !A) Danabasoglu,G. and J. C.McWilliams, 1995:
  !!         ! Sensitivity of the global ocean circulation to 
  !!         ! parameterizations of mesoscale tracer transports
  !!         ! Journal of Climate, 8, 2967-2987
  !!         ! For steep slope regions:
  !!         !    Exponential taper applied to both the neutral and GM operator
  !!         !    parts.
  !!
  !!         !B) Large, W.G. et al. 1997
  !!         ! Sensitivity to surface forcing and boundary layer mixing in a global 
  !!         !  ocean model: annual-mean climatology
  !!         ! JPO, 27, 2418-2447
  !!         ! For steep slope regions:
  !!         !    Exponential taper applied to both neutral operator and GM operator
  !!         ! For near surface part:
  !!         !    Sine taper also applied to both neutral operator and GM operator. 
  !!         !C) Grffies, Fundamentals of Ocean Climte Models, 2004
  !!         ! For steep slope region:
  !!         ! a) no taper applied to diagonal piece of horizontal neutral operator
  !!         ! b) hyperbolic tangent(exponential) taper applied to off-diagonal piece of
  !!         !    horizontal operator and to diagonal and off-diagonal piece of vertical
  !!         !    neutral diffusion operator. a)+b) means we transfer the tracer diffusion
  !!         !    to a horizontal-vertical manner in regions of steep neutral slopes.
  !!         ! c) Exponential taper applied to GM operator.
  !!         ! For surface layer with small slope:
  !!         ! a) sine taper applied to both neutral operator and GM operator, except the
  !!         !    diagonal piece of the horizontal diffusion.
  
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE calc_entries_mixing_tensor(patch_3d, ocean_state, param, op_coeff)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_ho_params),                 INTENT(inout) :: param
    TYPE(t_operator_coeff),            INTENT(in)    :: op_coeff
    
    !Local variables
    INTEGER :: start_cell_index, end_cell_index, cell_index,level,start_level,end_level, blockNo
    INTEGER :: start_edge_index, end_edge_index, je     
    TYPE(t_subset_range), POINTER :: cells_in_domain!,edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D    
    REAL(wp), POINTER :: K_I(:,:,:), K_D(:,:,:), kappa(:,:,:)
    REAL(wp)           :: geometric_scale, geometric_scale_factor_GMR
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain 
    start_level=1
    !-------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    !
    !The dianeutral diffusivity is the number determined by the PP-scheme
    K_I           => param%k_tracer_isoneutral
    K_D           => param%k_tracer_dianeutral !param%a_tracer_v(:,:,:,tracer_index) !  
    kappa         => param%k_tracer_GM_kappa
  
    ! by_Oliver 
    ! call Danabasoglu & Marshall (2007) parameterisation    
    IF ( lvertical_GM ) THEN
      CALL vertical_GM(patch_3d, ocean_state, param, op_coeff, kappa)
      !FIXME: use Redi coeff. = kappa
      K_I = kappa
    ENDIF
 
    IF(SWITCH_OFF_TAPERING)THEN
     ocean_state%p_aux%taper_function_1(:,:,:)=1.0_wp
     ocean_state%p_aux%taper_function_2(:,:,:)=1.0_wp
    ENDIF
    taper_diagonal_horz(:,:,:)     =0.0_wp
    taper_diagonal_vert_expl(:,:,:)=0.0_wp        
    taper_diagonal_vert_impl(:,:,:)=0.0_wp    
    taper_off_diagonal_horz(:,:,:)%x(1)=0.0_wp
    taper_off_diagonal_horz(:,:,:)%x(2)=0.0_wp
    taper_off_diagonal_horz(:,:,:)%x(3)=0.0_wp    
    taper_off_diagonal_vert(:,:,:)%x(1)=0.0_wp
    taper_off_diagonal_vert(:,:,:)%x(2)=0.0_wp
    taper_off_diagonal_vert(:,:,:)%x(3)=0.0_wp
    !-------------------------------------------------------------------------------
    
    SELECT CASE(tapering_scheme)

    CASE(tapering_DanaMcWilliams)
      IF(SWITCH_ON_TAPERING_HORIZONTAL_DIFFUSION)THEN

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

        DO cell_index = start_cell_index, end_cell_index
          end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)
          
            DO level = start_level, end_level
!IF( ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)>S_max**2)THEN             
! ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)=S_max**2/ocean_state%p_aux%slopes_squared(cell_index,level,blockNo) 
! ELSE
!ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)=1.0_wp
!ENDIF                         
              !Following recommendations in Griffies Ocean Climate Model book (sect. 15.3.4.2)
              !Danabasoglou-McWilliams tapering is for horizontal flux only applied to
              !off-diagonal terms
                            
              !coefficients for horizontal fluxes
              taper_diagonal_horz(cell_index,level,blockNo)  &
              & =K_I(cell_index,level,blockNo)&
              &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)

              taper_off_diagonal_horz(cell_index,level,blockNo)%x&
              & = (K_I(cell_index,level,blockNo)-kappa(cell_index,level,blockNo))&
              &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&
              &*ocean_state%p_aux%slopes(cell_index,level,blockNo)%x
 
              !coefficients for vertical fluxes
              taper_off_diagonal_vert(cell_index,level,blockNo)%x&
              &= (K_I(cell_index,level,blockNo)+kappa(cell_index,level,blockNo))&
              &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&                
              &*ocean_state%p_aux%slopes(cell_index,level,blockNo)%x
              
              taper_diagonal_vert_impl(cell_index,level,blockNo)  &              
              &=K_I(cell_index,level,blockNo)&
              &*ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)&
              &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)

            END DO
        END DO
      END DO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL      
      
      ELSEIF(.NOT.SWITCH_ON_TAPERING_HORIZONTAL_DIFFUSION)THEN

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

        DO cell_index = start_cell_index, end_cell_index
          end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

            DO level = start_level, end_level
            
              !Following recommendations in Griffies Ocean Climate Model book (sect. 15.3.4.2)
              !Danabasoglou-McWilliams tapering is for horizontal flux only applied to
              !off-diagonal terms
                           
              !coefficients for horizontal fluxes
              taper_diagonal_horz(cell_index,level,blockNo)  &
              & = K_I(cell_index,level,blockNo)

              taper_off_diagonal_horz(cell_index,level,blockNo)%x&
              & = (K_I(cell_index,level,blockNo)-kappa(cell_index,level,blockNo))&
              &*ocean_state%p_aux%slopes(cell_index,level,blockNo)%x&
              &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)

              !coefficients for vertical fluxes
              taper_off_diagonal_vert(cell_index,level,blockNo)%x&
              &= (K_I(cell_index,level,blockNo)+kappa(cell_index,level,blockNo))&
              &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&                
              &*ocean_state%p_aux%slopes(cell_index,level,blockNo)%x

              
              taper_diagonal_vert_impl(cell_index,level,blockNo)  &              
              &=K_I(cell_index,level,blockNo)&
              &*ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)&
              &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)

            END DO

        END DO
      END DO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL      
!       CALL sync_patch_array(sync_c, patch_2D, taper_diagonal_vert_impl)
      ENDIF
      IF(switch_off_diagonal_vert_expl)THEN
!ICON_OMP_PARALLEL      
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level) ICON_OMP_DEFAULT_SCHEDULE
        DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

          DO cell_index = start_cell_index, end_cell_index
            end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

!             IF(end_level >= min_dolic) THEN
          
              DO level = start_level, end_level
                taper_diagonal_vert_expl(cell_index,level,blockNo)=0.0_wp
              END DO
!             ENDIF
          END DO
        END DO      
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL
      ELSEIF(.NOT.switch_off_diagonal_vert_expl)THEN
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level) ICON_OMP_DEFAULT_SCHEDULE
       DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

          DO cell_index = start_cell_index, end_cell_index
            end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

!             IF(end_level >= min_dolic) THEN
          
              DO level = start_level, end_level
      
                taper_diagonal_vert_expl(cell_index,level,blockNo)  &
                &=K_D(cell_index,level,blockNo)&
                &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)
               END DO
!               ENDIF
            END DO
          END DO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL
      ENDIF            

    CASE(tapering_Large)
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level, geometric_scale) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

        DO cell_index = start_cell_index, end_cell_index
          end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

          IF(end_level >= min_dolic) THEN
             !geometric_scale=1.0_wp!sqrt(patch_2D%cells%area(cell_index,blockNo))*geometric_scale_factor_GMR
            DO level = start_level, end_level
            
              !coefficients for horizontal fluxes           
              taper_diagonal_horz(cell_index,level,blockNo)  &
              &=K_I(cell_index,level,blockNo)

              taper_off_diagonal_horz(cell_index,level,blockNo)%x&
              & = ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&
              & * ocean_state%p_aux%taper_function_2(cell_index,level,blockNo)&
              & *(K_I(cell_index,level,blockNo)-kappa(cell_index,level,blockNo))&
              &*ocean_state%p_aux%slopes(cell_index,level,blockNo)%x

              !coefficients for vertical fluxes
              taper_off_diagonal_vert(cell_index,level,blockNo)%x&
              &= (K_I(cell_index,level,blockNo)+kappa(cell_index,level,blockNo))&
              &+ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&
              &*ocean_state%p_aux%taper_function_2(cell_index,level,blockNo)&
              &*ocean_state%p_aux%slopes(cell_index,level,blockNo)%x

              
              taper_diagonal_vert_impl(cell_index,level,blockNo)  &              
              &=K_I(cell_index,level,blockNo)&
              &*ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)&
              &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&
              &*ocean_state%p_aux%taper_function_2(cell_index,level,blockNo)
 
            END DO
          ENDIF
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO
! !ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level) ICON_OMP_DEFAULT_SCHEDULE     
      IF(switch_off_diagonal_vert_expl)THEN
        DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

          DO cell_index = start_cell_index, end_cell_index
            end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

            IF(end_level >= min_dolic) THEN
          
              DO level = start_level, end_level
                taper_diagonal_vert_expl(cell_index,level,blockNo)=0.0_wp
              END DO
            ENDIF
          END DO
        END DO      
      ELSEIF(.NOT.switch_off_diagonal_vert_expl)THEN
       DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

          DO cell_index = start_cell_index, end_cell_index
            end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

            IF(end_level >= min_dolic) THEN
          
              DO level = start_level, end_level
      
                taper_diagonal_vert_expl(cell_index,level,blockNo)  &
                &=K_D(cell_index,level,blockNo)&
                &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)
               END DO
              ENDIF
            END DO
          END DO
        
      ENDIF      
! !ICON_OMP_END_PARALLEL_DO
               

    CASE(tapering_Griffies)
    CALL finish(TRIM('mo_ocean_GM_Redi'), 'The tapering option = tapering_Griffies is not supported yet')    
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,level, geometric_scale) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

        DO cell_index = start_cell_index, end_cell_index
          end_level = patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

          IF(end_level >= min_dolic) THEN
            geometric_scale=1.0_wp!sqrt(patch_2D%cells%area(cell_index,blockNo))*geometric_scale_factor_GMR
            DO level = start_level, end_level
            
              taper_diagonal_horz(cell_index,level,blockNo)  &
                &=K_I(cell_index,level,blockNo)

              taper_diagonal_vert_expl(cell_index,level,blockNo)  &
                &=K_D(cell_index,level,blockNo)              

               taper_diagonal_vert_impl(cell_index,level,blockNo)  &              
                 &=(K_I(cell_index,level,blockNo)*ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)&
                 &+K_D(cell_index,level,blockNo))*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&
                 &*ocean_state%p_aux%taper_function_2(cell_index,level,blockNo)              
              
              taper_off_diagonal_horz(cell_index,level,blockNo)%x&
                = (K_I(cell_index,level,blockNo)-kappa(cell_index,level,blockNo))&
                &*ocean_state%p_aux%slopes(cell_index,level,blockNo)%x&
                &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&
                &*ocean_state%p_aux%taper_function_2(cell_index,level,blockNo)

              taper_off_diagonal_vert(cell_index,level,blockNo)%x&
                = (K_I(cell_index,level,blockNo)+kappa(cell_index,level,blockNo))&
                &*ocean_state%p_aux%slopes(cell_index,level,blockNo)%x&
                &*ocean_state%p_aux%taper_function_1(cell_index,level,blockNo)&
                &*ocean_state%p_aux%taper_function_2(cell_index,level,blockNo)
            
            END DO
          ENDIF
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO

    END SELECT

!           IF(tracer_index<=2)THEN
        ! 
        !1.) Interpolate the tapered coefficient for the vertical tracer flux from prism center to prism top: 
        !this is the diagonal part that is handled implicitely in time
      
    IF(.NOT.REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN        
      CALL map_scalar_center2prismtop( patch_3d, &
      &                              taper_diagonal_vert_impl,&
      &                              op_coeff,                &
      &                              ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit)        

    ELSEIF(REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN
      CALL map_scalar_center2prismtop_GM( patch_3d, &
      &                              taper_diagonal_vert_impl,&
      &                              op_coeff,                &
      &                              ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit)        

    ENDIF

    CALL dbg_print('mix_tensor: horz diag', taper_diagonal_horz,&
    & this_mod_name, 3, patch_2D%cells%in_domain)
    CALL dbg_print('mix_tensor: vert diag expl', taper_diagonal_vert_expl,&
    & this_mod_name, 3, patch_2D%cells%in_domain)
    CALL dbg_print('mix_tensor: vert diag impl', taper_diagonal_vert_impl,&
    & this_mod_name, 3, patch_2D%cells%in_domain)

  END SUBROUTINE calc_entries_mixing_tensor
  !-------------------------------------------------------------------------
  
  

 
  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates the fluxes of the isoycnical diffusion following Redi and the eddy flux of GM.
  !!
   !!
  !<Optimize:inUse>
  SUBROUTINE diagnose_Redi_flux_balance(patch_3d, ocean_state, param, op_coeff,&
    &tracer_index, typeOfTracers)
    TYPE(t_patch_3d ),TARGET, INTENT(in)     :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET        :: ocean_state
    TYPE(t_ho_params),      INTENT(inout)    :: param
    TYPE(t_operator_coeff), INTENT(in)       :: op_coeff
    INTEGER, INTENT(IN)                      :: tracer_index 
    CHARACTER(LEN=*), INTENT(in)             :: typeOfTracers
    
    !Local variables
    INTEGER :: start_cell_index, end_cell_index, cell_index,level,start_level,end_level,blockNo
    INTEGER :: start_edge_index, end_edge_index, je     
    TYPE(t_subset_range), POINTER :: all_cells, cells_in_domain, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp) :: flux_vert_center(nproma, n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates) :: flux_vec_horz_center(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates),POINTER :: tracer_gradient_horz_vec_center(:,:,:), slopes(:,:,:)
    REAL(wp), POINTER :: tracer_gradient_vert_center(:,:,:)

    TYPE(t_cartesian_coordinates)            :: balance1(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)    
    TYPE(t_cartesian_coordinates)            :: balance2(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)        
    TYPE(t_cartesian_coordinates)            :: taper_off_diagonal_vert(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)    
    TYPE(t_cartesian_coordinates)            :: taper_off_diagonal_horz(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)                                 :: taper_diagonal_horz(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)                                 :: taper_diagonal_vert_expl(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)    
    REAL(wp)                                 :: taper_diagonal_vert_impl(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)            
    REAL(wp) :: mapped_vertical_diagonal_impl(nproma,n_zlev+1,patch_3D%p_patch_2D(1)%alloc_cell_blocks)     
    TYPE(t_cartesian_coordinates) :: flux_sum_horz(nproma,n_zlev, patch_3d%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)                      :: flux_sum_vert(nproma,n_zlev, patch_3d%p_patch_2D(1)%alloc_cell_blocks)                      
    REAL(wp) :: flux_vert_temp(nproma,n_zlev, patch_3d%p_patch_2D(1)%alloc_cell_blocks)           
    REAL(wp) :: flux_vert_sal(nproma,n_zlev, patch_3d%p_patch_2D(1)%alloc_cell_blocks)           
    REAL(wp):: neutral_coeff(1:n_zlev, 2), salinityColumn(1:n_zlev)
    REAL(wp), POINTER :: pot_temp(:,:,:), salinity(:,:,:), tracer(:,:,:)
    REAL(wp), POINTER :: depth_cellinterface(:,:,:)
    
    REAL(wp) :: temp_array(nproma,n_zlev, patch_3d%p_patch_2D(1)%alloc_cell_blocks)  
    
    REAL(wp)        :: GMRedi_flux_horz(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp)        :: GMRedi_flux_vert(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_e)
    
    TYPE(t_ocean_tracer) :: temp_tracer_before!(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)       
    TYPE(t_ocean_tracer) :: temp_tracer_after!(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)           
    REAL(wp) :: div_diff_flux_horz(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: div_diff_flux_vert(nproma,n_zlev, patch_3d%p_patch_2d(1)%alloc_cell_blocks)                    
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain 
    all_cells       => patch_2D%cells%all
    edges_in_domain => patch_2D%edges%in_domain 
    slopes          => ocean_state%p_aux%slopes 
    depth_cellinterface => patch_3D%p_patch_1d(1)%depth_cellinterface
    
    ALLOCATE(temp_tracer_before%concentration(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks))
    ALLOCATE(temp_tracer_after%concentration(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks))      
            
    GMRedi_flux_horz(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%nblks_e)=0.0_wp
    GMRedi_flux_vert(1:nproma,1:n_zlev,1:patch_3D%p_patch_2D(1)%nblks_e)=0.0_wp
    
    balance1(:,:,:)%x(1)=0.0_wp
    balance1(:,:,:)%x(2)=0.0_wp
    balance1(:,:,:)%x(3)=0.0_wp
    balance2(:,:,:)%x(1)=0.0_wp
    balance2(:,:,:)%x(2)=0.0_wp
    balance2(:,:,:)%x(3)=0.0_wp
            
    temp_array(1:nproma,1:n_zlev,1:patch_3d%p_patch_2D(1)%alloc_cell_blocks)=0.0_wp
    
   IF(no_tracer>=2)THEN
      pot_temp          => ocean_state%p_prog(nold(1))%tracer_collection%tracer(1)%concentration
      salinity          => ocean_state%p_prog(nold(1))%tracer_collection%tracer(2)%concentration
    ELSEIF(no_tracer==1)THEN      
      pot_temp          => ocean_state%p_diag%rho
    ENDIF
    
    start_level=1

    CALL calc_entries_mixing_tensor(patch_3d, ocean_state, param, op_coeff)
!    ENDIF   
             
    IF(no_tracer<=3)THEN

      flux_sum_horz(1:nproma,1:n_zlev, 1:patch_3d%p_patch_2D(1)%alloc_cell_blocks)%x(1)=0.0_wp 
      flux_sum_horz(1:nproma,1:n_zlev, 1:patch_3d%p_patch_2D(1)%alloc_cell_blocks)%x(2)=0.0_wp         
      flux_sum_horz(1:nproma,1:n_zlev, 1:patch_3d%p_patch_2D(1)%alloc_cell_blocks)%x(3)=0.0_wp
      
      flux_sum_vert(:,:,:) =0.0_wp                
      !flux_vert_temp(:,:,:)=0.0_wp  


      IF(tracer_index==1)THEN
        
        tracer                          => pot_temp
        tracer_gradient_horz_vec_center => ocean_state%p_aux%PgradTemperature_horz_center
        tracer_gradient_vert_center     => ocean_state%p_aux%DerivTemperature_vert_center
      ELSEIF(tracer_index==2)THEN
        
        flux_vert_sal(1:nproma,1:n_zlev, 1:patch_3d%p_patch_2D(1)%alloc_cell_blocks)=0.0_wp        

        tracer                          => salinity
        tracer_gradient_horz_vec_center => ocean_state%p_aux%PgradSalinity_horz_center
        tracer_gradient_vert_center     => ocean_state%p_aux%DerivSalinity_vert_center


    ELSEIF(tracer_index>2)THEN
! write(0,*)'-----------------------------BALANCE-DIAG-TRACER==3-----------------------------------------'    
ocean_state%p_prog(nold(1))%tracer_collection%tracer(tracer_index)%concentration=ocean_state%p_diag%rho_GM
!&=ocean_state%p_prog(nold(1))%tracer_collection%tracer(tracer_index-1)%concentration 
        tracer                          => ocean_state%p_prog(nold(1))%tracer_collection%tracer(tracer_index)%concentration
      !Here we have to provide a sbr that calculates derivatives below    
      CALL calc_tracer_derivatives( patch_3d,&
                                  & ocean_state%p_prog(nold(1))%tracer_collection%tracer(tracer_index)%concentration,&
                                  & ocean_state, &
                                  & op_coeff, &
                                  & tracer_index, &
                                  & typeOfTracers)

      tracer_gradient_horz_vec_center => ocean_state%p_aux%PgradTracer_horz_center
      tracer_gradient_vert_center     => ocean_state%p_aux%DerivTracer_vert_center

      ENDIF

      DO blockNo = all_cells%start_block, all_cells%end_block        
        CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)      
        DO cell_index = start_cell_index, end_cell_index
          DO level = 1, n_zlev
            flux_vec_horz_center(cell_index,level,blockNo)%x = 0.0_wp
          ENDDO
        ENDDO
      ENDDO
!ICON_OMP_DO_PARALLEL PRIVATE(start_cell_index,end_cell_index, cell_index, level) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block        
        CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)      
        DO cell_index = start_cell_index, end_cell_index
 
          !horizontal GMRedi Flux at top layer: with the tapering this is the just horizontal diffusion
          DO level = start_level, MIN(patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo),start_level)
            flux_vec_horz_center(cell_index,start_level,blockNo)%x &
              &=taper_diagonal_horz(cell_index,start_level,blockNo) &
              &*tracer_gradient_horz_vec_center(cell_index,start_level,blockNo)%x

!             flux_vert_center(cell_index,start_level,blockNo) &
!               &=taper_diagonal_vert_expl(cell_index,start_level,blockNo)&
!               &*tracer_gradient_vert_center(cell_index,start_level,blockNo)
!            ! the top level flux_vert_center will be filled from the second level, if it exists
             flux_vert_center(cell_index,start_level,blockNo) = 0.0_wp
          ENDDO

          DO level = start_level+1, patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)
          
            !horizontal GM-Redi Flux
            flux_vec_horz_center(cell_index,level,blockNo)%x &
              &=taper_diagonal_horz(cell_index,level,blockNo)*  &
              &tracer_gradient_horz_vec_center(cell_index,level,blockNo)%x&
              !the second term vanishes if GM-Kappa=isoneutral diffusion
              &+taper_off_diagonal_horz(cell_index,level,blockNo)%x&
              &*tracer_gradient_vert_center(cell_index,level,blockNo)
              
            !vertical GM-Redi Flux: this is the part that is explicit in time
            !If namelist option "switch_off_diagonal_vert_expl=TRUE" the
            !first (diagonal) term is set to zero. Default option "switch_off_diagonal_vert_expl=FALSE"
            !Second term contains isoneutral and GM components
            !
            !For diagnostic purposes everything is handled exlicitely.
            flux_vert_center(cell_index,level,blockNo)= &
              &Dot_Product(tracer_gradient_horz_vec_center(cell_index,level,blockNo)%x,&
              &            taper_off_diagonal_vert(cell_index,level,blockNo)%x)
                
          END DO
          ! fill the top level flux_vert_center from the second level, if it exists
!           DO level = start_level+1, MIN(patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo),start_level+1)
!             flux_vert_center(cell_index,start_level,blockNo) = flux_vert_center(cell_index,start_level+1,blockNo)
!           END DO
        END DO                
      END DO
!ICON_OMP_END_DO_PARALLEL
        
    ! use a vector communicator
    CALL sync_patch_array_mult(sync_c, patch_2D, 3, &
        & flux_vec_horz_center(:,:,:)%x(1), flux_vec_horz_center(:,:,:)%x(2), flux_vec_horz_center(:,:,:)%x(3))


    IF(INCLUDE_SLOPE_SQUARED_IMPLICIT)THEN

      ! Now we treat the vertical isoneutral coefficient that is discretized implicitely in time.  
      ! This is only neccessary once for temperature and salinity, the HAMOCC tracers use these value  
      IF(tracer_index<=2)THEN
        ! 
        !1.) Interpolate the tapered coefficient for the vertical tracer flux from prism center to prism top: 
        !this is the diagonal part that is handled implicitely in time
      
        IF(.NOT.REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN        
          CALL map_scalar_center2prismtop( patch_3d, &
          &                              taper_diagonal_vert_impl,&
          &                              op_coeff,                &
          &                              ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit)        

        ELSEIF(REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN
          CALL map_scalar_center2prismtop_GM( patch_3d, &
          &                              taper_diagonal_vert_impl,&
          &                              op_coeff,                &
          &                              ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit)        
    
        ENDIF      
        !

      ENDIF!IF(tracer_index<=1)THEN 

      !Assess implicit contribution for analysis
      CALL map_cell2edges_3D( patch_3D,flux_vec_horz_center, GMRedi_flux_horz, op_coeff)
      div_diff_flux_horz(:,:,:)=0.0_wp
      CALL div_oce_3d( GMRedi_flux_horz(:,:,:),&
                   &   patch_3d, &
                   &   op_coeff%div_coeff, &
                   &   div_diff_flux_horz )

      !vertical div of explicit part of vertical GMRedi-flux
      div_diff_flux_vert(:,:,:) = 0.0_wp
      CALL verticalDiv_scalar_onFullLevels( patch_3d, &
        & GMRedi_flux_vert(:,:,:), &
        & div_diff_flux_vert)

      temp_tracer_before%concentration = tracer &
      &- (dtime /  patch_3d%p_patch_1D(1)%prism_thick_c)    &
      & * (- div_diff_flux_horz-div_diff_flux_vert)

      temp_tracer_after%concentration=temp_tracer_before%concentration

      CALL tracer_diffusion_vertical_implicit( &
          & patch_3d,                        &
          & temp_tracer_before,              &
          & ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit, &
          & ocean_state%p_prog(nnew(1))%h)

!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, level) ICON_OMP_DEFAULT_SCHEDULE
        DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      
          CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      
          DO cell_index = start_cell_index, end_cell_index
            DO level = start_level+1, patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)
              !vertical GM-Redi Flux
               GMRedi_flux_vert(cell_index,level,blockNo)   &
                &=flux_vert_center(cell_index,level,blockNo) &
                & + (temp_tracer_after%concentration(cell_index,level,blockNo)&
                &   -temp_tracer_before%concentration(cell_index,level,blockNo))/dtime 
            END DO                  
          END DO                
        END DO
!ICON_OMP_END_DO            
      
    ELSEIF(.NOT.INCLUDE_SLOPE_SQUARED_IMPLICIT)THEN
        
!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, level) ICON_OMP_DEFAULT_SCHEDULE
        DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      
          CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      
          DO cell_index = start_cell_index, end_cell_index
            DO level = start_level+1, patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)
              !vertical GM-Redi Flux
               GMRedi_flux_vert(cell_index,level,blockNo)   &
                &=flux_vert_center(cell_index,level,blockNo) &
                & + tracer_gradient_vert_center(cell_index,level,blockNo)&
                & *taper_diagonal_vert_impl(cell_index,level,blockNo)
            END DO                  
          END DO                
        END DO
!ICON_OMP_END_DO            
    ENDIF  



! !ICON_OMP_DO_PARALLEL PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,neutral_coeff, &
! !ICON_OMP  level) ICON_OMP_DEFAULT_SCHEDULE

!ICON_OMP_DO_PARALLEL PRIVATE(start_cell_index,end_cell_index, cell_index, end_level,neutral_coeff,level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)

      DO cell_index = start_cell_index, end_cell_index
        end_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
        IF(end_level <= min_dolic) CYCLE

        
         IF (no_tracer>=2)THEN
         
           salinityColumn(1:end_level) = salinity(cell_index,1:end_level,blockNo)
           
           !Nonlinear EOS 
           IF(EOS_TYPE/=1)THEN
             !neutral_coeff = calc_neutralslope_coeff_func_onColumn(         &
             !& pot_temp(cell_index,1:end_level,blockNo), salinityColumn(1:end_level),  &
             !& depth_cellinterface(cell_index,2:end_level+1,blockNo), end_level)

             neutral_coeff = calc_neutralslope_coeff_func_onColumn(         &
             & pot_temp(cell_index,1:end_level,blockNo), salinityColumn(1:end_level),  &
             & depth_cellinterface(cell_index,2:end_level+1,blockNo), end_level)

             
           !Linear EOS: slope coefficients are equal to EOS-coefficients
           ELSEIF(EOS_TYPE==1)THEN
              neutral_coeff(:,1) = LinearThermoExpansionCoefficient 
              neutral_coeff(:,2) = LinearHalineContractionCoefficient       
            ENDIF
          
         ELSEIF(no_tracer==1)THEN         
           neutral_coeff(:,1)=1.0         
         ENDIF 


         IF(tracer_index==1)THEN
           DO level = start_level+1, end_level-1
                       
            ocean_state%p_aux%diagnose_Redi_flux_temp(cell_index,level,blockNo)%x&
            &=-neutral_coeff(level,1)*flux_vec_horz_center(cell_index,level,blockNo)%x
            
            ocean_state%p_aux%diagnose_Redi_flux_vert(cell_index,level,blockNo)&
            &=-neutral_coeff(level,1)*GMRedi_flux_vert(cell_index,level,blockNo)
 
            
           END DO
         ELSEIF(tracer_index==2)THEN
          DO level = start_level+1, end_level-1         
           
           ocean_state%p_aux%diagnose_Redi_flux_sal(cell_index,level,blockNo)%x&
           &=neutral_coeff(level,2)*flux_vec_horz_center(cell_index,level,blockNo)%x  

            !horizontal balance relation alphatimes Redi flux (temperature)=beta times Redi flux(salt)           
            flux_sum_horz(cell_index,level,blockNo)%x&
            &=ocean_state%p_aux%diagnose_Redi_flux_temp(cell_index,level,blockNo)%x&
            &+ ocean_state%p_aux%diagnose_Redi_flux_sal(cell_index,level,blockNo)%x
            !length of horizontal balance vector 
            temp_array(cell_index,level,blockNo)&
            &=sqrt(dot_product(flux_sum_horz(cell_index,level,blockNo)%x,flux_sum_horz(cell_index,level,blockNo)%x))
 
            !vertical balance relation alphatimes Redi flux (temperature)=beta times Redi flux(salt)           
            flux_sum_vert(cell_index,level,blockNo)  = ocean_state%p_aux%diagnose_Redi_flux_vert(cell_index,level,blockNo)&
            &+neutral_coeff(level,2)*GMRedi_flux_vert(cell_index,level,blockNo)
            
            
            !This is about Flux horzontal tines slopes = flux vertical        
            balance1(cell_index,level,blockNo)%x=ocean_state%p_aux%slopes(cell_index,level,blockNo)%x&
            &*(-neutral_coeff(level,1)*ocean_state%p_aux%DerivTemperature_vert_center(cell_index,level,blockNo)&
            & + neutral_coeff(level,2)*ocean_state%p_aux%DerivSalinity_vert_center(cell_index,level,blockNo))
                       
            balance2(cell_index,level,blockNo)%x=&
            &-neutral_coeff(level,1)*ocean_state%p_aux%PgradTemperature_horz_center(cell_index,level,blockNo)%x&
            &+neutral_coeff(level,2)*ocean_state%p_aux%PgradSalinity_horz_center(cell_index,level,blockNo)%x
            
            !balance1(cell_index,level,blockNo)%x&
            !&=param%k_tracer_isoneutral(cell_index,level,blockNo)*balance1(cell_index,level,blockNo)%x
            !balance2(cell_index,level,blockNo)%x&
            !&=param%k_tracer_isoneutral(cell_index,level,blockNo)*balance2(cell_index,level,blockNo)%x  
            IF(ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)>0.0_wp)THEN        
            write(1234,*)'Balance: Horizontal: Vertical: Slope Squared',level,&
            temp_array(cell_index,level,blockNo),&
            flux_sum_vert(cell_index,level,blockNo),&
             !dot_product(balance1(cell_index,level,blockNo)%x,balance1(cell_index,level,blockNo)%x),&
             !dot_product(balance2(cell_index,level,blockNo)%x,balance2(cell_index,level,blockNo)%x),&             
!            &balance1(cell_index,level,blockNo)%x(1),&
!            &balance1(cell_index,level,blockNo)%x(2),&
!            &balance1(cell_index,level,blockNo)%x(3),&
!            &balance2(cell_index,level,blockNo)%x(1),&
!            &balance2(cell_index,level,blockNo)%x(2),&
!            &balance2(cell_index,level,blockNo)%x(3),&
            &sqrt(ocean_state%p_aux%slopes_squared(cell_index,level,blockNo))
            !&ocean_state%p_aux%slopes(cell_index,level,blockNo)%x(1),&
            !&ocean_state%p_aux%slopes(cell_index,level,blockNo)%x(2),&
            !&ocean_state%p_aux%slopes(cell_index,level,blockNo)%x(3)  
            
            !write(1234,*)'det B',level,&        
            !&param%k_tracer_isoneutral(cell_index,level,blockNo),ocean_state%p_aux%taper_function_1(cell_index,level,blockNo),&
            !&sqrt(ocean_state%p_aux%slopes_squared(cell_index,level,blockNo))
            ENDIF
          END DO                        
         ENDIF
          
      END DO ! cell_index = start_cell_index, end_cell_index
    END DO  ! blockNo = all_cells%start_block, all_cells%end_block
!ICON_OMP_END_DO_PARALLEL



   !Map the explicit horizontal tracer flux from cell centers to edges (where the horizontal divergence is calculated)
    !
    ! use a vector communicator
    CALL sync_patch_array_mult(sync_c, patch_2D, 3, &
        & flux_sum_horz(:,:,:)%x(1), flux_sum_horz(:,:,:)%x(2), flux_sum_horz(:,:,:)%x(3))
    !    


    !Map the (explicit) vertical tracer flux to the prism top (where the vertical divergence is calculated later)
    IF(.NOT.REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN
      CALL map_scalar_center2prismtop(patch_3d, flux_sum_vert, op_coeff,GMRedi_flux_vert)
    ELSEIF(REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN
      CALL map_scalar_center2prismtop_GM(patch_3d, flux_sum_vert, op_coeff,GMRedi_flux_vert)    
    ENDIF
    
    
     CALL map_cell2edges_3D( patch_3D,flux_sum_horz, GMRedi_flux_horz, op_coeff)
      div_diff_flux_horz(:,:,:)=0.0_wp
      CALL div_oce_3d( GMRedi_flux_horz(:,:,:),&
                   &   patch_3d, &
                   &   op_coeff%div_coeff, &
                   &   div_diff_flux_horz )

      !vertical div of explicit part of vertical GMRedi-flux
      div_diff_flux_vert(:,:,:) = 0.0_wp
      CALL verticalDiv_scalar_onFullLevels( patch_3d, &
        & flux_sum_vert(:,:,:), &
        & div_diff_flux_vert)
    
    
    
    IF(tracer_index==no_tracer)THEN

      Do level=start_level+1,n_zlev-1
        idt_src=1  ! output print level (1-5, fix)
        !write(0,*)'level',level
        ! CALL dbg_print('diag_Redi:horz',GMRedi_flux_horz(:,level,:),&
        !  & str_module, idt_src, in_subset=edges_in_domain)
        CALL dbg_print('diag_Redi:horz',temp_array(:,level,:),&
        & str_module, idt_src, in_subset=cells_in_domain)

      END DO
    
      Do level=start_level+1,n_zlev-1
        idt_src=1  ! output print level (1-5, fix)
        !write(0,*)'level',level,tracer_index
        CALL dbg_print('diag_Redi:vert',GMRedi_flux_vert(:,level,:),&
        & str_module, idt_src, in_subset=cells_in_domain)
      END DO
    
      Do level=start_level+1,n_zlev-1
        idt_src=1  ! output print level (1-5, fix)
        CALL dbg_print('div_Redi:horz',div_diff_flux_horz(:,level,:),&
        & str_module, idt_src, in_subset=cells_in_domain)

      END DO
    
      Do level=start_level+1,n_zlev-1
        idt_src=1  ! output print level (1-5, fix)
        !write(0,*)'level',level,tracer_index
        CALL dbg_print('div_Redi:vert',div_diff_flux_vert(:,level,:),&
        & str_module, idt_src, in_subset=cells_in_domain)
      END DO    
    ENDIF
    !---------------------------------------------------------------------
    
!      CALL sync_patch_array(sync_e, patch_2D, GMRedi_flux_horz(:,:,:))
!      CALL sync_patch_array(sync_c, patch_2D, GMRedi_flux_vert(:,:,:))
     
  ELSEIF( no_tracer>2)THEN
    CALL finish(TRIM('calc_GMRediflux'),&
    & 'diagnose_Redi_flux_balance beyond temperature and salinity is not impemented yet')
  ENDIF


   DEALLOCATE(temp_tracer_before%concentration)
   DEALLOCATE(temp_tracer_after%concentration) 



  END SUBROUTINE diagnose_Redi_flux_balance
  !-------------------------------------------------------------------------

!!$!-------------------------------------------------------------------------
!!$  !>
!!$  !! !  SUBROUTINE calculates the bolus velocity from the isopycnal slopes.
!!$  !! !  This is currently used for diagnostic purposes, the GM-Redi paramerization
!!$  !! !  uses the skew-flux approach.
!!$  !!
!!$  !!
!!$!<Optimize:inUse>
!!$  SUBROUTINE calc_bolus_velocity(patch_3d, ocean_state, param, op_coeff)
!!$    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
!!$    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
!!$    TYPE(t_ho_params),                 INTENT(inout) :: param
!!$    TYPE(t_operator_coeff),            INTENT(inout) :: op_coeff
!!$    
!!$    !Local variables
!!$    REAL(wp) :: div_bolus_vn(nproma, n_zlev+1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
!!$    REAL(wp) :: cell_max_slope, inv_cell_characteristic_length,slope_abs
!!$    
!!$    TYPE(t_cartesian_coordinates) :: grad_slope_vec(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)    
!!$    TYPE(t_cartesian_coordinates) :: kappa_times_slope(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
!!$    TYPE(t_cartesian_coordinates) :: slope_deriv_atcenter(nproma, n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
!!$           
!!$    INTEGER :: level, blockNo, start_level, end_level
!!$    INTEGER :: start_cell_index, end_cell_index, cell_index
!!$    TYPE(t_subset_range), POINTER :: cells_in_domain, all_cells!, edges_in_domain
!!$    TYPE(t_patch), POINTER :: patch_2D
!!$    !-----------------------------------------------------------------------
!!$    patch_2D        => patch_3D%p_patch_2D(1)
!!$    all_cells       => patch_2D%cells%all
!!$    cells_in_domain => patch_2D%cells%in_domain
!!$    !edges_in_domain => patch_2D%edges%in_domain
!!$    !-------------------------------------------------------------------------
!!$
!!$    start_level = 1
!!$    
!!$    
!!$    !Step 1: Calculate kappa times slopes
!!$    !-------------------------------------------------------------------------------
!!$!ICON_OMP_PARALLEL
!!$!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, end_level, &
!!$!ICON_OMP  level) ICON_OMP_DEFAULT_SCHEDULE
!!$          DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
!!$            CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
!!$
!!$            DO cell_index = start_cell_index, end_cell_index
!!$            
!!$              inv_cell_characteristic_length = 1.0_wp / SQRT(patch_2D%cells%area(cell_index,blockNo))
!!$              end_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo)
!!$              IF(end_level <= min_dolic) CYCLE
!!$                              
!!$              DO level = start_level+1, end_level-1
!!$
!!$                cell_max_slope      = S_max  &
!!$                & * patch_3d%p_patch_1d(1)%prism_thick_c(cell_index,level,blockNo) &
!!$                & * inv_cell_characteristic_length
!!$              
!!$                slope_abs = sqrt(ocean_state%p_aux%slopes_squared(cell_index,level,blockNo))
!!$
!!$                IF(slope_abs <= cell_max_slope)THEN
!!$                  kappa_times_slope(cell_index,level,blockNo)%x=        &
!!$                  &-param%k_tracer_GM_kappa(cell_index,level,blockNo)   &
!!$                  &+ocean_state%p_aux%slopes(cell_index,level,blockNo)%x
!!$                ELSE
!!$                  kappa_times_slope(cell_index,level,blockNo)%x=0.0_wp                      
!!$                ENDIF
!!$              END DO                         
!!$            END DO ! cell_index = start_cell_index, end_cell_index
!!$          END DO  ! blockNo = all_cells%start_block, all_cells%end_block
!!$!ICON_OMP_END_DO_NOWAIT
!!$!ICON_OMP_END_PARALLEL
!!$    !-------------------------------------------------------------------------------
!!$
!!$    !Step 2: Calculate vertical derivative of kappa times slopes and map it back to prism center(midlevel)
!!$    !-------------------------------------------------------------------------------           
!!$!ICON_OMP_PARALLEL   
!!$!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index) ICON_OMP_DEFAULT_SCHEDULE
!!$ 
!!$    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
!!$      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
!!$        
!!$      CALL verticalDeriv_vec_midlevel_on_block( patch_3d,                      &
!!$                                              & kappa_times_slope(:,:,blockNo),&
!!$                                              & grad_slope_vec(:,:,blockNo),   &
!!$                                              & start_level+1,                 &
!!$                                              & blockNo,                       &
!!$                                              & start_cell_index,              &
!!$                                              & end_cell_index)
!!$        IF(.NOT.REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN                                       
!!$          CALL map_vec_prismtop2center_on_block( patch_3d,                   &
!!$            & grad_slope_vec(:,:,blockNo),       &
!!$            & slope_deriv_atcenter(:,:,blockNo), & 
!!$            & blockNo, start_cell_index, end_cell_index)
!!$        ELSE
!!$          CALL map_vec_prismtop2center_on_block_GM( patch_3d,                &
!!$            & grad_slope_vec(:,:,blockNo),              &
!!$            & slope_deriv_atcenter(:,:,blockNo),        & 
!!$            & blockNo, start_cell_index, end_cell_index)
!!$               
!!$       
!!$       ENDIF                                                                           
!!$                                              
!!$    END DO ! blocks
!!$!ICON_OMP_END_DO_NOWAIT
!!$!ICON_OMP_END_PARALLEL        
!!$    !------------------------------------------------------------------------------
!!$
!!$    !Step 3: Map result back to edges to obtain the horizontal bolus velocity 
!!$    !------------------------------------------------------------------------------
!!$    CALL  map_cell2edges_3d( patch_3d,                   &
!!$                           & slope_deriv_atcenter,       &
!!$                           & ocean_state%p_diag%vn_bolus,&
!!$                           & op_coeff,                   &
!!$                           & start_level+1)
!!$    !------------------------------------------------------------------------------
!!$
!!$
!!$    !Step 4: Calculate vertical bolus as horizontal divergence of horizontal bolus velocity
!!$    !------------------------------------------------------------------------------
!!$!ICON_OMP_PARALLEL    
!!$!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index,end_cell_index) ICON_OMP_DEFAULT_SCHEDULE
!!$          DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
!!$            CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
!!$        
!!$            CALL div_oce_3D_onTriangles_onBlock( &
!!$              & ocean_state%p_diag%vn_bolus,     &
!!$              & patch_3D,op_coeff%div_coeff,     &
!!$              & div_bolus_vn(:,:,blockNo),       &
!!$              & blockNo,start_cell_index, end_cell_index,start_level+1,n_zlev)!,      &
!!$              !& start_level=1, end_level=n_zlev)
!!$        
!!$            DO cell_index = start_cell_index, end_cell_index          
!!$              !use bottom boundary condition for vertical velocity at bottom of prism
!!$              DO level = patch_3d%p_patch_1d(1)%dolic_c(cell_index,blockNo), 1, -1
!!$                ocean_state%p_diag%w_bolus(cell_index,level,blockNo) &
!!$                  &= ocean_state%p_diag%w_bolus(cell_index,level+1,blockNo)&
!!$                  & - div_bolus_vn(cell_index,level,blockNo)
!!$              END DO
!!$            END DO
!!$        
!!$          END DO ! blockNo
!!$!ICON_OMP_END_PARALLEL_DO 
!!$!ICON_OMP_END_PARALLEL
!!$    !------------------------------------------------------------------------------
!!$!    CALL map_edges2cell_3d(patch_3d, ocean_state%p_diag%vn_bolus, op_coeff,ocean_state%p_diag%p_vn)
!!$!    CALL Get3DVectorTo2DLocal_array3D(vector=ocean_state%p_diag%p_vn, &
!!$!      & position_local=patch_2d%cells%center, &
!!$!      & levels=patch_3d%p_patch_1d(1)%dolic_c, &
!!$!      & subset=all_cells,                      &
!!$!      & geometry_info=patch_2d%geometry_info, &
!!$!      & x=ocean_state%p_diag%u, y=ocean_state%p_diag%v)
!!$
!!$  !---------DEBUG DIAGNOSTICS-------------------------------------------
!!$  idt_src=1  ! output print level (1-5, fix)
!!$  CALL dbg_print('calc_bolus:',(ocean_state%p_diag%w_bolus(:,:,:)),&
!!$    & str_module,idt_src, in_subset=cells_in_domain)
!!$   DO level=1,n_zlev
!!$     CALL dbg_print('calc_bolus:',ocean_state%p_diag%w_bolus(:,level,:),&
!!$       & str_module,idt_src, in_subset=cells_in_domain)
!!$  END DO
!!$!stop
!!$  !---------------------------------------------------------------------
!!$
!!$  END SUBROUTINE calc_bolus_velocity
!!$  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE Parameterisation for a vertical variation of the GM coefficient (kappa)
  !! following Danabasoglu and Marshall (2007, doi:10.1016/j.ocemod.2007.03.006).
  !! Written  by Oliver Gutjahr, MPI-M (2020).
  !!
  !! The parameterisation defines a vertical variable kappa:
  !!  K = (N2 / N2_ref) * K_ref
  !! with N2 the local buoyancy frequency computed at vertical interface levels,
  !! N2_ref the reference buoyancy frequency just below the diabatic layer
  !! depth (DLD) provided N2>0 there, otherwise N2_ref is first N2 >0 below.
  !!  K_ref is the reference kappa defined in the namelist.
  !!
  !! For all shallower depth than DLD, the ratio N2/N2_ref = 1, implying no
  !! vertical variation of K; in particular in the surface diabatic layer (or
  !! mixed layer as approximation).
  !!
  !! Between the depth at which N2 = N2_ref and the ocean bottom, it is ensured
  !! that
  !!   Nmin <= (N2 / N2_ref) <= 1.0
  !! where Nmin (>0) is a specified lower limit.
  !!
  !! In unstable regions (N2 < 0) N2/N2_ref = Nmin (0.1). This ensures the
  !! vertical smoothness of the N2/N2_ref profile and relies on the enhancement
  !! of vertical mixing coefficients to alleviate any local static instability.
  !!
  !! K and N2/N2_ref are calculated at every model time step. However, it could
  !! be relaxed to only once per simulation day.
  !!
  !! There is no vertical averaging of the N2/N2_ref profile; the value at the
  !! top of the grid cell is used for the hole grid cell.
  !!
  !! The resulting K is subject to usual diffusive numerical stability
  !! criterion.
  !! Danabasoglu and Marshall (2007) allow a maximum slope of 0.3. They further
  !! set K=0 in the bottom halves of the grid cell, which further reduces K in
  !! the bottom cell by one half.
  !!
  !! We do assume here that the DLD is sufficiently approximated with the mixed
  !! layer depth (mld; see Danabasoglu et al., 2008).

  !! As input we need N2, mld,k_tracer_GM_kappa_parameter and the absolute
  !! depths of the interface levels.

SUBROUTINE vertical_GM(patch_3d, ocean_state, param, op_coeff, kappa)
    TYPE(t_patch_3d ),TARGET, INTENT(in)     :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET        :: ocean_state
    TYPE(t_ho_params),      INTENT(inout)    :: param
    TYPE(t_operator_coeff), INTENT(in)       :: op_coeff
    REAL(wp),DIMENSION(:,:,:),INTENT(INOUT)  :: kappa

    ! Local variables
    INTEGER :: start_cell_index, end_cell_index, cell_index,level,start_level,end_level,blockNo
    INTEGER :: start_edge_index, end_edge_index, je, k, jc, jk
    TYPE(t_subset_range), POINTER :: all_cells, cells_in_domain, edges_in_domain
    INTEGER                                  :: k_mld                    ! interface level of mixed layer depth
    REAL(wp),DIMENSION(:,:), POINTER         :: mld                      ! diagnosed mixed layer depth
    INTEGER, POINTER                         :: kbot(:,:)                ! bottom level
    REAL(wp),POINTER                         :: N2(:,:,:)                ! buoyancy frequency
    REAL(wp),POINTER                         :: interface_depth(:,:,:)   ! depth of interface levels
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp)   :: N2_ref(nproma,patch_3D%p_patch_2D(1)%nblks_e)          ! reference N2 value just below 
                                                                         !  the surface diabatic layer (mixed layer)
    REAL(wp)   :: Nratio(nproma,n_zlev+1,patch_3D%p_patch_2D(1)%nblks_e) ! ratio of N2/N2_ref
    !-------------------------------------------------------------------------------

    ! These are necessary to calculate N2, can be removed once N2 becomes a
    ! global variable
    REAL(wp), POINTER :: temp(:,:,:)                                        ! potential temperature
    REAL(wp), POINTER :: salt(:,:,:)                                        ! salinity 
    REAL(wp), POINTER :: dzi(:,:,:)                                         ! inverse distance of cell centres

    REAL(wp),TARGET :: Nsqr(nproma,n_zlev+1,patch_3D%p_patch_2D(1)%nblks_e) ! temporary N2
    REAL(wp) :: rho_up(n_zlev), rho_down(n_zlev)                            ! temporary density
    REAL(wp) :: pressure(n_zlev), salinity(n_zlev)                          ! temporary pressure and salinity

    !------------------------------------------------------------------------------- 

    ! Parameters for zstar
    REAL(wp) :: s_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    !stretching factor
    REAL(wp) :: e_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)    !surface height



    patch_2D        => patch_3d%p_patch_2d(1)
    all_cells       => patch_2D%cells%all
 
    cells_in_domain => patch_2D%cells%in_domain
    start_level=1

    ! set zstar related parameters
    if (vert_cor_type == 1) then
      s_c(:,:) = ocean_state%p_prog(nold(1))%stretch_c(:,:)
      e_c(:,:) = ocean_state%p_prog(nold(1))%eta_c(:,:)
    else
      s_c = 1.0_wp
      e_c = 0.0_wp
    endif


    ! Nmin is used from namelist

    ! N2 has to be calculated (pp,tke) since it is not a global variable yet. 


         ! N2 has to be calculated (for all cases)
         temp => ocean_state%p_prog(nold(1))%tracer(:,:,:,1)
         salt => ocean_state%p_prog(nold(1))%tracer(:,:,:,2)
         dzi => patch_3d%p_patch_1d(1)%inv_prism_center_dist_c
         interface_depth => patch_3d%p_patch_1d(1)%depth_CellInterface(:,:,:)
         Nsqr(:,:,:) = 0.0_wp
         DO blockNo = all_cells%start_block, all_cells%end_block
           CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)
           DO jc = start_cell_index, end_cell_index
             pressure(1:n_zlev) = (interface_depth(jc,1:n_zlev,blockNo) * s_c(jc,blockNo) &
             &                  - e_c(jc,blockNo)) * ReferencePressureIndbars 
             rho_up(1:n_zlev-1)  = calculate_density_onColumn(&
               & temp(jc,1:n_zlev-1,blockNo), &
               & salt(jc,1:n_zlev-1,blockNo), &
               & pressure(2:n_zlev), n_zlev-1)
             rho_down(2:n_zlev)  = calculate_density_onColumn(&
               & temp(jc,2:n_zlev,blockNo), &
               & salt(jc,2:n_zlev,blockNo), &
               & pressure(2:n_zlev), n_zlev-1)
             DO jk = 2, n_zlev
               Nsqr(jc,jk,blockNo) = grav/OceanReferenceDensity * (rho_down(jk) - rho_up(jk-1)) * dzi(jc,jk,blockNo)/s_c(jc,blockNo)
             END DO
           END DO
         END DO
         N2 => Nsqr(:,:,:)  
 

    ! mixed layer depth
    mld => ocean_state%p_diag%mld 

    ! bottom level
    kbot => patch_3d%p_patch_1d(1)%dolic_c(:,:)

    ! nullify
    N2_ref(:,:) = 0.0_wp
    Nratio(:,:,:) = 1.0_wp

    ! begin calculation
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, blockNo, start_cell_index,end_cell_index)

          DO cell_index = start_cell_index, end_cell_index

            !---------------------------------------------
            ! 1) find interface level of mixed layer depth
            !---------------------------------------------
            k_mld = 1

            DO level = start_level,kbot(cell_index,blockNo)

               k_mld = level

               DO k = 1,k_mld
                 IF (k+1 > kbot(cell_index,blockNo)) THEN
                   k_mld = kbot(cell_index,blockNo)
                   EXIT
                 ELSE IF ( interface_depth(cell_index,k+1,blockNo) >= mld(cell_index,blockNo) ) THEN
                   k_mld = k
                   EXIT
                 END IF   
               END DO
            
            END DO

            !---------------------------------------------
            ! 2) find N2_ref := first N2>0 below mld
            !---------------------------------------------
              DO level = k_mld+1, patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo) 

                IF ( N2(cell_index,level,blockNo) > 0 ) THEN

                  N2_ref(cell_index,blockNo) = N2(cell_index,level,blockNo)
                  EXIT
                  ! FIXME: it seems that a safety is needed when no N2 is found

                END IF
                ! FIXME: what to do if there is no N2>0 in the water column?
  

              END DO

              ! diagnostic for N2_ref
              ocean_state%p_diag%N2_ref(cell_index,blockNo) = N2_ref(cell_index,blockNo)

 

            !---------------------------------------------
            ! 3) calculate N2/N2_ref ratio (with safeties)
            !    and the vertical variable GM kappa
            !---------------------------------------------
            ! set all kappa=kappa_ref in the mixed layer (assuming N2/N2_ref = 1.0)
            kappa(cell_index,1:k_mld,blockNo)  = k_tracer_GM_kappa_parameter            

            ! modify kappa below the mixed layer with the N2 ratio
            DO level = k_mld+1,patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

              ! calculate the N2 ratio, constrained to a minimum value of Nmin
              ! and a maximum value of 1.0.
              IF ( N2_ref(cell_index,blockNo) .EQ. 0 ) THEN

                  Nratio(cell_index,level,blockNo) = Nmin
              ELSE 
                  Nratio(cell_index,level,blockNo) = MIN(MAX(N2(cell_index,level,blockNo) / N2_ref(cell_index,blockNo), Nmin),1.0_wp)

              END IF

              ! calculate vertical GM kappa 
              kappa(cell_index,level,blockNo) = Nratio(cell_index,level,blockNo) * k_tracer_GM_kappa_parameter

              !--- diagnostics 
              ! Nratio
              ocean_state%p_diag%Nratio(cell_index,level,blockNo) = Nratio(cell_index,level,blockNo)

              ! kappa
              ocean_state%p_diag%kappa_GM(cell_index,level,blockNo) = kappa(cell_index,level,blockNo)

            END DO

          END DO
    END DO


END SUBROUTINE vertical_GM

!---------------------------------------------------------------------------




  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates the fluxes of the isoycnical diffusion following Redi.
  !!
  !!  Adapted for zstar
  !!
!<Optimize:inUse:done>
  SUBROUTINE prepare_ocean_physics_zstar(patch_3d, ocean_state, param, op_coeff, stretch_c)
    TYPE(t_patch_3d ),TARGET, INTENT(in)             :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET                :: ocean_state
    TYPE(t_ho_params),                 INTENT(inout) :: param
    TYPE(t_operator_coeff),            INTENT(in) :: op_coeff
    REAL(wp),                          INTENT(IN) :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    
    INTEGER  :: jc, jb, level, start_level, end_level 
    TYPE(t_patch), POINTER :: patch_2D
    TYPE(t_subset_range), POINTER :: cells_in_domain
    INTEGER :: start_cell_index, end_cell_index
   !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain
   !-------------------------------------------------------------------------------

    CALL calc_neutral_slopes(patch_3d, ocean_state, param, op_coeff)

    !! Convert vertical gradients from d/ds to d/dz
    start_level = 1
!ICON_OMP_DO_PARALLEL PRIVATE(start_cell_index,end_cell_index, jc, level, end_level) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index
        end_level = patch_3d%p_patch_1d(1)%dolic_c(jc, jb)
        DO level = start_level, end_level
           ocean_state%p_aux%slopes(jc, level, jb)%x &
             & = stretch_c(jc, jb) * ocean_state%p_aux%slopes(jc, level, jb)%x 

           ocean_state%p_aux%slopes_drdz(jc,level,jb) &
             & = (1.0_wp / stretch_c(jc, jb)) * ocean_state%p_aux%slopes_drdz(jc,level,jb)

           ocean_state%p_aux%slopes_squared(jc, level, jb) &
             & = stretch_c(jc, jb) * stretch_c(jc, jb) * ocean_state%p_aux%slopes_squared(jc, level, jb)
        END DO

      END DO
    END DO
!ICON_OMP_END_DO_PARALLEL

    CALL calc_tapering_function(patch_3d, param, ocean_state)
    
    IF(TEST_MODE_REDI_ONLY.OR.TEST_MODE_GM_ONLY)THEN
      CALL finish(TRIM('mo_ocean_GM_Redi'), 'This GMRedi_configuration is not supported')
    ELSE
      CALL calc_entries_mixing_tensor(patch_3d, ocean_state, param, &
                      & op_coeff)
    ENDIF
     
  END SUBROUTINE prepare_ocean_physics_zstar
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! !  SUBROUTINE calculates the fluxes of the isoycnical diffusion following Redi.
  !!
  !!  Adapted for zstar
  !!
  !<Optimize:inUse:done>
  SUBROUTINE calc_ocean_physics_zstar(patch_3d, ocean_state, param, op_coeff, GMRedi_flux_horz, &
      & GMRedi_flux_vert, tracer, tracer_index, typeOfTracers, stretch_c, stretch_e)
    TYPE(t_patch_3d ),TARGET, INTENT(in)     :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET        :: ocean_state
    TYPE(t_ho_params), INTENT(inout)         :: param
    TYPE(t_operator_coeff),INTENT(in)        :: op_coeff
    REAL(wp), INTENT(inout)                  :: GMRedi_flux_horz(:,:,:)
    REAL(wp), INTENT(inout)                  :: GMRedi_flux_vert(:,:,:)
    TYPE(t_ocean_tracer), TARGET             :: tracer
    INTEGER, INTENT(IN)                      :: tracer_index
    CHARACTER(LEN=*), INTENT(in)             :: typeOfTracers
    REAL(wp),                   INTENT(IN)   :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks) 
    REAL(wp),                   INTENT(IN)   :: stretch_e(nproma, patch_3D%p_patch_2d(1)%nblks_e)
    
    !Local variables
    INTEGER :: start_cell_index, end_cell_index, cell_index, blockNo,level
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    
!    REAL(wp), POINTER, DIMENSION(:,:,:)                :: GMredi_flux_horz
!    REAL(wp), POINTER, DIMENSION(:,:,:)     :: GMredi_flux_vert
    !-------------------------------------------------------------------------------
    patch_2D         => patch_3d%p_patch_2d(1)
    cells_in_domain  => patch_2D%cells%in_domain   
!    GMredi_flux_horz => ocean_state%p_diag%GMRedi_flux_horz(:,:,:,tracer_index)
!    GMredi_flux_vert => ocean_state%p_diag%GMRedi_flux_vert(:,:,:,tracer_index)

        
    SELECT CASE(GMRedi_configuration)!GMRedi_configuration==Cartesian_Mixing)RETURN
 
    CASE(GMRedi_combined)
    
    
        CALL calc_combined_GentMcWilliamsRedi_flux_zstar( patch_3d,        &
            & ocean_state,      &
            & param,            &
            & op_coeff,         &
            & GMRedi_flux_horz, &
            & GMRedi_flux_vert, &
            & tracer,           &
            & tracer_index, typeOfTracers, &
            & stretch_c, stretch_e)
    CASE DEFAULT
      CALL finish(TRIM('mo_ocean_GM_Redi'), 'This GMRedi_configuration is not supported')
    
    END SELECT
    
    
  END SUBROUTINE calc_ocean_physics_zstar
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE calculates the fluxes of the isoycnical diffusion following Redi and the eddy flux of GM.
  !!
  !!  Adapted for zstar
  !!
  SUBROUTINE calc_combined_GentMcWilliamsRedi_flux_zstar(patch_3d, ocean_state, param, op_coeff,&
      &GMredi_flux_horz, GMredi_flux_vert, tracer, tracer_index, typeOfTracers, stretch_c, stretch_e)
    TYPE(t_patch_3d ),TARGET, INTENT(in)     :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET        :: ocean_state
    TYPE(t_ho_params),      INTENT(inout)    :: param
    TYPE(t_operator_coeff), INTENT(in)       :: op_coeff
    REAL(wp), INTENT(inout)                  :: GMredi_flux_horz(:,:,:)
    REAL(wp), INTENT(inout)                  :: GMredi_flux_vert(:,:,:)
    TYPE(t_ocean_tracer), TARGET             :: tracer
    INTEGER, INTENT(IN)                      :: tracer_index
    CHARACTER(LEN=*), INTENT(in)             :: typeOfTracers
    REAL(wp),                   INTENT(IN)   :: stretch_c(nproma, patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp),                   INTENT(IN)   :: stretch_e(nproma, patch_3D%p_patch_2d(1)%nblks_e)

    !Local variables
    INTEGER :: start_cell_index, end_cell_index, cell_index,level,start_level,end_level,blockNo,max_level
    INTEGER :: start_edge_index, end_edge_index, edge_index
    INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    TYPE(t_subset_range), POINTER :: all_cells, cells_in_domain, edges_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    REAL(wp) :: flux_vert_center(nproma, n_zlev,patch_3D%p_patch_2d(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates) :: flux_vec_horz_center(nproma,n_zlev,patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates),POINTER :: tracer_gradient_horz_vec_center(:,:,:), slopes(:,:,:)
    REAL(wp), POINTER :: tracer_gradient_vert_center(:,:,:)

    REAL(wp) :: mapped_vertical_diagonal_impl(nproma,n_zlev+1,patch_3D%p_patch_2D(1)%alloc_cell_blocks)      
    
    REAL(wp) :: nabla_T_horz(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_e)
    REAL(wp) :: nabla_S_horz(nproma, n_zlev,patch_3D%p_patch_2d(1)%nblks_e)   
    !-------------------------------------------------------------------------------
    patch_2D        => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2D%cells%in_domain 
    all_cells       => patch_2D%cells%all
    edges_in_domain => patch_2D%edges%in_domain 
    slopes          => ocean_state%p_aux%slopes 

    start_level=1

    IF (typeOfTracers == "ocean" )THEN
     !Set pointers for tracer gradients, according to actual tracer index
     IF(tracer_index==1)THEN
      !write(0,*) "DerivTemperature_vec_center"
      tracer_gradient_horz_vec_center => ocean_state%p_aux%PgradTemperature_horz_center
      tracer_gradient_vert_center     => ocean_state%p_aux%DerivTemperature_vert_center
     ELSEIF(tracer_index==2)THEN
      !write(0,*) "DerivSalinity_vec_center"
      tracer_gradient_horz_vec_center => ocean_state%p_aux%PgradSalinity_horz_center
      tracer_gradient_vert_center     => ocean_state%p_aux%DerivSalinity_vert_center
     ENDIF
    ELSE

      !Here we have to provide a sbr that calculates derivatives below    
      CALL calc_tracer_derivatives( patch_3d,&
                                  & tracer%concentration,&
                                  & ocean_state, &
                                  & op_coeff, &
                                  & tracer_index, &
                                  & typeOfTracers)


      tracer_gradient_horz_vec_center => ocean_state%p_aux%PgradTracer_horz_center
      tracer_gradient_vert_center     => ocean_state%p_aux%DerivTracer_vert_center

      !! Convert vertical gradients from d/ds to d/dz
!ICON_OMP_PARALLEL_DO PRIVATE(start_cell_index, end_cell_index, cell_index, level, end_level, max_level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      DO cell_index = start_cell_index, end_cell_index
        end_level = patch_3d%p_patch_1d(1)%dolic_c(cell_index, blockNo)
        DO level = start_level, end_level
          tracer_gradient_vert_center(cell_index, level, blockNo) &
            & = (1.0_wp/stretch_c(cell_index, blockNo)) * tracer_gradient_vert_center(cell_index, level, blockNo)
         END DO
      END DO
    END DO
!ICON_OMP_END_DO_PARALLEL

    ENDIF

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)
      DO cell_index = start_cell_index, end_cell_index
        DO level = 1, n_zlev
          flux_vec_horz_center(cell_index,level,blockNo)%x = 0.0_wp
        ENDDO
      ENDDO
    ENDDO

    !The code below has to be executed for temperature, salinity and HAMMOC tracers, i.e. for all tracers.
    !
!ICON_OMP_DO_PARALLEL PRIVATE(start_cell_index,end_cell_index, cell_index, level, max_level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
      DO cell_index = start_cell_index, end_cell_index
 
        !GMRedi Flux at top layer: with tapering this is the just horizontal diffusion
        IF(patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)>=start_level)THEN
           flux_vec_horz_center(cell_index,start_level,blockNo)%x &
                   &=taper_diagonal_horz(cell_index,start_level,blockNo)  &
                   &*tracer_gradient_horz_vec_center(cell_index,start_level,blockNo)%x

          ! the top level flux_vert_center will be filled from the second level, if it exists
          flux_vert_center(cell_index,start_level,blockNo) = 0.0_wp
        ENDIF
       DO level = start_level+1, patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)

          !horizontal GM-Redi Flux
          IF(ocean_state%p_aux%slopes_squared(cell_index,level,blockNo)==0.0_wp)THEN
          tracer_gradient_horz_vec_center(cell_index,level,blockNo)%x(:)=0.0_wp
          ENDIF
          flux_vec_horz_center(cell_index,level,blockNo)%x &
          &=taper_diagonal_horz(cell_index,level,blockNo)  &
          &*tracer_gradient_horz_vec_center(cell_index,level,blockNo)%x&
          !the second term vanishes if GM-Kappa=isoneutral diffusion
          &+taper_off_diagonal_horz(cell_index,level,blockNo)%x&
          &*tracer_gradient_vert_center(cell_index,level,blockNo)

              
          !vertical GM-Redi Flux: this is the part that is explicit in time
          !If namelist option "switch_off_diagonal_vert_expl=TRUE" the
          !first (diagonal) term is set to zero. Default option "switch_off_diagonal_vert_expl=FALSE"
          !Second term contains isoneutral and GM components
          flux_vert_center(cell_index,level,blockNo)= &
          &taper_diagonal_vert_expl(cell_index,level,blockNo)&
          &*tracer_gradient_vert_center(cell_index,level,blockNo)&
          &+&
          &Dot_Product(tracer_gradient_horz_vec_center(cell_index,level,blockNo)%x,&
          &            taper_off_diagonal_vert(cell_index,level,blockNo)%x)

        END DO

      END DO                
    END DO
!ICON_OMP_END_DO_PARALLEL

    !Map the explicit horizontal tracer flux from cell centers to edges (where the horizontal divergence is calculated)
    !
    ! use a vector communicator
    CALL sync_patch_array_mult(sync_c, patch_2D, 3, &
        & flux_vec_horz_center(:,:,:)%x(1), flux_vec_horz_center(:,:,:)%x(2), flux_vec_horz_center(:,:,:)%x(3))
    !
    CALL map_cell2edges_3D( patch_3D,flux_vec_horz_center, GMredi_flux_horz, op_coeff)

!ICON_OMP_DO_PARALLEL PRIVATE(start_edge_index,end_edge_index, edge_index, level, max_level) ICON_OMP_DEFAULT_SCHEDULE
        DO blockNo = edges_in_domain%start_block, edges_in_domain%end_block
          CALL get_index_range(edges_in_domain, blockNo, start_edge_index, end_edge_index)
          DO edge_index = start_edge_index, end_edge_index
            DO level = start_level, patch_3D%p_patch_1D(1)%dolic_e(edge_index,blockNo)
              GMredi_flux_horz(edge_index,level,blockNo)   &
              & = GMredi_flux_horz(edge_index,level,blockNo)  &
              & * patch_3D%p_patch_1d(1)%prism_thick_e(edge_index,level,blockNo)  &
              & * stretch_e(edge_index, blockNo)
            END DO
          END DO
        END DO
!ICON_OMP_END_DO_PARALLEL


    IF(INCLUDE_SLOPE_SQUARED_IMPLICIT)THEN

    !Map the (explicit) vertical tracer flux to the prism top (where the vertical divergence is calculated later)
    IF(.NOT.REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN
      CALL map_scalar_center2prismtop(patch_3d, flux_vert_center, op_coeff,GMredi_flux_vert)
    ELSEIF(REVERT_VERTICAL_RECON_AND_TRANSPOSED)THEN
      CALL map_scalar_center2prismtop_GM(patch_3d, flux_vert_center, op_coeff,GMredi_flux_vert)    
    ENDIF


      ! Now we treat the vertical isoneutral coefficient that is discretized implicitely in time.  
      ! This is only neccessary once for temperature and salinity, the HAMOCC tracers use these value  
      IF(typeOfTracers == "ocean" )THEN
        IF(tracer_index==1)THEN
          CALL dbg_print('Old vert coeff: A_v', param%a_tracer_v(:,:,:, tracer_index),&
          & this_mod_name, 4, patch_2D%cells%in_domain)
          !CALL dbg_print('Old vert coeff: A_v', ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit,&
          !& this_mod_name, 4, patch_2D%cells%in_domain)
        ENDIF
        !
        !2.) Here we combine the vertical GMRedicoefficient that is treated implicitely (mapped_vertical_diagonal_impl, this
        !term involves the slopes-squared times the isopycnal mixing coefficient and is potentially large, therefore
        !it is discretized implicitely) with the vertical mixing coefficient from the PP-scheme.
        !We follow the approach in POP, where these two contributions are added
        !(see Reference manual POP, sect 5.1.3, in particular p. 41, after eq (150)).
        !
!ICON_OMP_DO_PARALLEL PRIVATE(start_cell_index,end_cell_index, cell_index, level, max_level) ICON_OMP_DEFAULT_SCHEDULE
        DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
          CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
          DO cell_index = start_cell_index, end_cell_index
            DO level = start_level, patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)
              param%a_tracer_v(cell_index,level,blockNo, tracer_index) =  &
              & param%a_tracer_v(cell_index,level,blockNo, tracer_index) + &
              & ocean_state%p_diag%vertical_mixing_coeff_GMRedi_implicit(cell_index,level,blockNo)
            END DO
          END DO
        END DO
!ICON_OMP_END_DO_PARALLEL

    CALL sync_patch_array(sync_c, patch_2D, param%a_tracer_v(:,:,:,tracer_index))


        IF(tracer_index==1)THEN
          CALL dbg_print('New vert coeff: A_v', param%a_tracer_v(:,:,:, tracer_index),&
          & this_mod_name, 4, patch_2D%cells%in_domain)
        ENDIF
      ENDIF!IF(typeOfTracers == "ocean" )THEN

      ELSEIF(.NOT.INCLUDE_SLOPE_SQUARED_IMPLICIT)THEN

!ICON_OMP_DO PRIVATE(start_cell_index,end_cell_index, cell_index, level, max_level) ICON_OMP_DEFAULT_SCHEDULE
        DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block

          CALL get_index_range(cells_in_domain, blockNo, start_cell_index, end_cell_index)
          DO cell_index = start_cell_index, end_cell_index
            DO level = start_level+1, patch_3D%p_patch_1D(1)%dolic_c(cell_index,blockNo)
              !vertical GM-Redi Flux
               GMredi_flux_vert(cell_index,level,blockNo)   &
                &=flux_vert_center(cell_index,level,blockNo) &
                & + tracer_gradient_vert_center(cell_index,level,blockNo)&
                & *taper_diagonal_vert_impl(cell_index,level,blockNo)

            END DO
          END DO
        END DO
!ICON_OMP_END_DO

    ENDIF
    !---------DEBUG DIAGNOSTICS-------------------------------------------
    !Do level=1,n_zlev
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('InGMRedi_zstar: GMRedi_vert',GMredi_flux_vert,&
    & str_module, idt_src, in_subset=cells_in_domain)
    !END DO
    !Do level=1,n_zlev
    CALL dbg_print('InGMRedi_zstar: GMRedi_horz',GMredi_flux_horz,&
    & str_module, idt_src, in_subset=edges_in_domain)
    !END DO
    !---------------------------------------------------------------------
    

  END SUBROUTINE calc_combined_GentMcWilliamsRedi_flux_zstar
  !-------------------------------------------------------------------------
 

END MODULE mo_ocean_GM_Redi

