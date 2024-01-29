! Contains the variables to set up the ocean model.
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
#include "iconfor_dsl_definitions.inc"
!=============================================================================================
MODULE mo_ocean_types

  USE mo_kind,                ONLY: wp, sp
  USE mo_impl_constants,      ONLY: land, land_boundary, boundary, sea_boundary, sea,  &
    & success, max_char_length, min_dolic,               &
    & full_coriolis, beta_plane_coriolis,                &
    & f_plane_coriolis, zero_coriolis, halo_levels_ceiling
  USE mo_math_types,          ONLY: t_cartesian_coordinates,      &
    & t_geographical_coordinates
  USE mo_ocean_diagnostics_types, ONLY: t_ocean_monitor
  USE mo_model_domain,        ONLY: t_patch_3d
  USE mo_ocean_tracer_transport_types

  PUBLIC :: t_hydro_ocean_base
  PUBLIC :: t_hydro_ocean_state
  PUBLIC :: t_hydro_ocean_prog
  PUBLIC :: t_hydro_ocean_diag
  PUBLIC :: t_hydro_ocean_aux
  PUBLIC :: t_onCells_Pointer_3d_wp, t_onCells_HalfLevels_Pointer_wp, t_onEdges_Pointer_3d_wp
  PUBLIC :: t_oce_config
  PUBLIC :: t_verticalAdvection_ppm_coefficients
  

  PUBLIC :: t_operator_coeff
  PUBLIC :: t_solverCoeff_singlePrecision

  TYPE t_onCells_Pointer_3d_wp
    onCells :: p  ! pointer to 3D array
  END TYPE t_onCells_Pointer_3d_wp
  TYPE t_onCells_HalfLevels_Pointer_wp
    onCells_HalfLevels :: p  ! pointer to 3D array
  END TYPE t_onCells_HalfLevels_Pointer_wp
  TYPE t_onEdges_Pointer_3d_wp
    onCells :: p  ! pointer to 3D array
  END TYPE t_onEdges_Pointer_3d_wp
  TYPE t_onEdges_HalfLevels_Pointer_wp
    onEdges_HalfLevels :: p  ! pointer to 3D array
  END TYPE t_onEdges_HalfLevels_Pointer_wp
  
  !-------------------------------------------------------------------------------
  TYPE t_verticalAdvection_ppm_coefficients
    !  coefficients for the upwind_vflux_ppm vertical advection
    !  these are allocated in a block mode (ie each block allocates its own coefficients)
    !  all dimensions are (nproma, levels),
    !  although not all the levels are actually used
    onCellsBlock ::  cellHeightRatio_This_toBelow
    onCellsBlock ::  cellHeightRatio_This_toThisBelow
    onCellsBlock ::  cellHeight_2xBelow_x_RatioThis_toThisBelow
    onCellsBlock ::  cellHeightRatio_This_toThisAboveBelow
    onCellsBlock ::  cellHeightRatio_2xAboveplusThis_toThisBelow
    onCellsBlock ::  cellHeightRatio_2xBelowplusThis_toThisAbove
    onCellsBlock ::  cellHeightRatio_ThisAbove_to2xThisplusBelow
    onCellsBlock ::  cellHeightRatio_ThisBelow_to2xThisplusAbove
    onCellsBlock ::  cellHeight_inv_ThisAboveBelow2Below

  END TYPE t_verticalAdvection_ppm_coefficients
  !-------------------------------------------------------------------------------

!   TYPE t_pointer_2d_wp
!     REAL(wp),POINTER :: p(:,:)   ! pointer to 2D array
!   END TYPE t_pointer_2d_wp
  !
  !! basis types for constructing 3-dim ocean state
  !
  TYPE t_hydro_ocean_base
    
    !! The ocean uses z-coordinates in meters in the vertical.
    !! The following data are required:
    !!
    !! n_zlev: number of z-coordinate surfaces
    !! n_zlvp: number of intermediate levels (+1)
    !! n_zlvm: number of z-coordinate distances (-1)
    INTEGER :: n_zlev, n_zlvp, n_zlvm
    
    !! del_zlev_m: thickness (height) of elemental prism, defined as the
    !!             distance between top and bottom of elemental prism,
    !!             i.e. the distance between two intermediate z-coordinate
    !!             surfaces. These data are provided by the user, all other
    !!             vertical information is calculated from this array of
    !!             thicknesses.
    !!             Dimension: n_zlev
    onGrid_1D :: del_zlev_m
    
    !! zlev_m    : position of the vertical cell centers, i.e. below zero surface;
    !!             Numbering starts from surface and increases downwards to bottom.
    !!             Dimension: n_zlev
    !!             At these surfaces the horizontal velocities, vorticity, divergence
    !!             and scalar variables are evaluated.
    onGrid_1D :: zlev_m
    
    
    !! zlev_i    : vertical position of the UPPER BORDER of the vertical cell
    !!             i.e. the position of the top of elemental prisms.
    !!             Position of first surface is 0.
    !!             Dimension: n_zlvp = n_zlev + 1
    !!             The vertical velocities are evaluated at such surfaces.
    onGrid_HalfLevels1D :: zlev_i
    
    !! del_zlev_i: distance between two z-coordinate surfaces. The first is
    !!             the distance from the ocean surface = zlev_m(1)
    !!             Dimension: n_zlev
    onGrid_1D :: del_zlev_i
    
    
    ! land-sea-mask for ocean has 3 dimensions (the 2nd is the number of
    ! vertical levels)
    ! sea=-2, sea_boundary=-1, boundary (edges only)=0, land_boundary=1, land=2
    !
    ! land-sea-mask for cell centers
    ! index1=1,nproma, index2=1,n_zlev, index3=1,alloc_cell_blocks
    onCells_3D_Int :: lsm_c
    ! land-sea-mask for cell edges
    ! index1=1,nproma, index2=1,n_zlev, index3=1,nblks_e
    onEdges_3D_Int :: lsm_e
    ! land-sea-mask for cell vertices
    ! index1=1,nproma, index2=1,n_zlev, index3=1,nblks_v
    ! INTEGER, ALLOCATABLE :: lsm_v(:,:,:)
    
    
    ! To simplify the acess to the required information within these loops
    ! we store an cell and edge based version of the deepest ocean layer
    ! in column. dolic_e(edge1) and dolic_c(cell1) are identical if 'edge1'
    ! is one of the edges of 'cell1'.
    ! If the ocean bottom is flat dolic_c and dolic_e are identical and equal
    ! to the number of z-coodinate surfaces.
    
    ! index1=1,nproma, index2=1,alloc_cell_blocks
    onCells_2D_Int :: dolic_c
    ! index1=1,nproma, index2=1,nblks_e
    onEdges_2D_Int :: dolic_e
    
    ! For diagnosis like stream functions and area calculations we add surface arrays
    ! index1=1,nproma, index2=1,alloc_cell_blocks
    onCells_2D_Int :: basin_c  ! basin information Atlantic/Indian/Pacific
    onCells_2D_Int :: regio_c  ! area information like tropical Atlantic etc.
    
    ! To simply set land points to zero we store additional 3-dim wet points
    ! dimensions as in lsm_oce:
    onCells :: wet_c  ! cell centers
    onEdges :: wet_e  ! cell edges
    
  END TYPE t_hydro_ocean_base
  
  !----------------------------------------------
  ! prognostic variables
  TYPE t_hydro_ocean_prog

    onCells_2D :: h
    onCells_2D :: eta_c !! For z* computations 
    onCells_2D :: stretch_c !! For z* computations 
    
    onEdges :: vn
    onCells_tracers :: tracer
     ! Ordering of tracers:
     !   1) pot_temp:= potential temperature, Unit: [deg C]
     !   2) salinity:= salinity, Unit [psu]
    
    TYPE(t_tracer_collection) :: tracer_collection
    
    TYPE(t_onCells_Pointer_3d_wp),ALLOCATABLE :: tracer_ptr(:)  !< pointer array: one pointer for each tracer
  END TYPE t_hydro_ocean_prog
  

  TYPE t_hydro_ocean_diag

    onCells ::          &
      & rho            ,& ! density. Unit: [kg/m^3]
      & rhopot         ,& ! potential density. Unit: [kg/m^3]
      & rho_GM         ,& ! potential density. Unit: [kg/m^3]      
      & grad_rho_PP_vert ,& ! vertical insitu density gradient. Unit: [kg/m^4]      
      & div_mass_flx_c ,& !
      & u              ,& ! reconstructed zonal velocity component. Unit [m/s]
      & v              ,& ! reconstructed meridional velocity component. Unit [m/s]
      & w_prismcenter  ,&         
!       & potential_vort_c ,& ! potential vorticity averaged to triangle cells. Unit [1/s]
      & kin            ,& ! kinetic energy. Unit [m/s].
!       & div            ,& ! divergence. Unit [m/s]
      & press_hyd      ,& ! hydrostatic pressure. Unit [m]
!       & temp_insitu    ,&
      & t,s          ,& ! dummy pointer for output variabless
      & Buoyancy_Freq  ,&
      & Richardson_Number,            &
      & osaltGMRedi,           &
      & opottempGMRedi,           &
      & div_of_GMRedi_flux,           &
      & div_of_GMRedi_flux_horizontal,&
      & div_of_GMRedi_flux_vertical,  &
      & div_of_GM_flux,               &
      & div_of_Redi_flux,             &
      & kappa_GM,                     & !by_Oliver
      & Nratio,                       & !by_Oliver
      & vertical_mixing_coeff_GMRedi_implicit,&
      & w_bolus,                      &
      & opottemptend,                 &
      & osalttend,                    &
      & delta_thetao, &
      & delta_so, &
      & uT, &  !< product of temperature and u-velocity
      & uS, &  !< product of salinity and u-velocity
      & uR, &  !< product of density and u-velocity
      & uu, &  !< square of u-velocity
      & vT, &  !< product of temperature and v-velocity
      & vS, &  !< product of salinity and v-velocity
      & vR, &  !< product of density and v-velocity
      & vv, &  !< square of  v-velocity
      & wT, &  !< product of temperature and w-velocity
      & wS, &  !< product of salinity and w-velocity
      & wR, &  !< product of density and w-velocity
      & ww, &  !< square of w-velocity
      & RR, &  !< square of density
      & SS, &  !< square of salinity
      & TT, &  !< square of temperature
      & uw, &  !< product of u-velocity and w-velocity
      & vw, &  !< product of v-velocity and w-velocity
      & uv, &  !< product of u-velocity and v-velocity
      & swrab, &  !< relative swr absorption factor
      & rsdoabsorb, &  !< Net Rate of Absorption of Shortwave Energy in Ocean Layer
!      & dummy3d,&
      & swr_frac, &  !< relative swr absorption factor from hamocc (LFB_BGC_OCE)
      & sigma0, &  !<  density anomaly (desitity - 1000)
      & heat_content_liquid_water,    &
      & odensitytend


    onCells_2D :: &
      & thick_c          ,& ! individual fluid column thickness at cells. Unit [m].
      & u_vint           ,& ! barotropic zonal velocity. Unit [m*m/s]
      & v_vint           ,& ! barotropic meridional velocity. Unit [m*m/s]
      & mld              ,& ! mixed layer depth [m].
      & mlotst           ,& ! mixed layer depth [m]. (CMIP6)
      & mlotstsq        ,& ! squared mixed layer depth [m]. (CMIP6)
      & mlotst10         ,& ! 10m mixed layer depth [m]. (EERIE)
      & mlotst10sq       ,& ! 10m squared mixed layer depth [m]. (EERIE)
      & condep           ,&! convection depth index
      & swsum            ,& ! swr fraction absorbed in the surface layer
      & heatabs          ,& ! ( total subsurface ) heating due to absorption [J m-2]
!      & dummy2d          ,&
      & heatabs_t        ,&
      & heatabb          ,&
      & hflR, &  !< product of net heatflux and density
      & fwR, &  !< product of fw flux and density
      & tauxU, &  !< product of x-windstress and u-velocity
      & tauyV, &  !< product of y-windstress and v-velocity
      & heatflux_rainevaprunoff, & !< heatflux implied by fw flux (fw enters with SST)
      & heat_content_snow ,&
      & heat_content_seaice ,&
      & delta_ice, & 
      & delta_snow, &
      & heat_content_total ,&
      & zos_square     ,&
      & ssh     ,&   ! diagnostic including sea ice
      & Rossby_Radius    ,&
      & Wavespeed_baroclinic ,&
      & N2_ref           ,& ! diag: first N2>0 below mld !by_Oliver 
      & global_moc       ,& ! MOC global
      & atlantic_moc     ,& ! MOC atlantic
      & pacific_moc      ,& ! MOC pacific
      & global_hfl       ,& ! implied ocean heat transport global
      & atlantic_hfl     ,& ! implied ocean heat transport atlantic
      & pacific_hfl      ,& ! implied ocean heat transport pacific
      & global_wfl       ,& ! implied ocean freshwater transport global
      & atlantic_wfl     ,& ! implied ocean freshwater transport atlantic
      & pacific_wfl      ,& ! implied ocean freshwater transport pacific
      & global_hfbasin   ,& ! northward ocean heat transport global
      & atlantic_hfbasin ,& ! northward ocean heat transport atlantic
      & pacific_hfbasin  ,& ! northward ocean heat transport pacific
      & global_sltbasin   ,& ! northward ocean salt transport global
      & atlantic_sltbasin ,& ! northward ocean salt transport atlantic
      & pacific_sltbasin     ! northward ocean salt transport pacific

   onCells_2D :: &
      & northernHemisphere ,&
      & southernHemisphere
      
    onCells_2D :: &
      & u_50m,    &
      & v_50m,    &
      & T_50m,    &
      & vort_f_cells_50m

    onCells_Type(t_cartesian_coordinates) :: &
      & p_vn              ! reconstructed velocity at cell center in cartesian coordinates
      
    onCells_2D_Type(t_cartesian_coordinates) :: &
      & p_mass_flux_sfc_cc  ! mass flux at surface in cartesian coordinates
      
    onCells_HalfLevels_tracers ::        &
      & GMRedi_flux_vert

    onCells_HalfLevels ::        &
      & zgrad_rho      ,& ! vertical density gradient. Unit: [kg/m^2] this is allocated on n_zlev
      & w              ,& ! vertical velocity. Unit [m/s].
      & w_old          ,& ! vertical velocity from previous timestep. Unit [m/s].
!       & w_prev         ,& ! vertical velocity at cells, from previous timestep. Unit [m/s]
!       & w_time_weighted,& ! predicted normal velocity vector at cells
      & cfl_vert       ,&   ! vertical cfl values
      & w_deriv

    onEdges_tracers :: &
      & GMRedi_flux_horz 
    
    onEdges :: &
      & mass_flx_e     ,& ! individual fluid column thickness at cells. Unit [m].
      & ptp_vn         ,& ! normal velocity after mapping P^T P, not used currently
      & vn_pred        ,& ! predicted normal velocity vector at edges.
      & vn_pred_ptp    ,& ! predicted normal velocity vector at edges.
      & vn_time_weighted,&  ! predicted normal velocity vector at edges.
      & vort_e          ,& ! vorticity interpolated to triangle edges. Unit [1/s]
!       & potential_vort_e,& ! potential vorticity at triangle edges. Unit [1/s]
      & veloc_adv_horz ,& ! horizontal velocity advection
      & veloc_adv_vert ,& ! vertical velocity advection
      & laplacian_horz ,& ! horizontal diffusion of horizontal velocity
      & laplacian_vert ,& ! vertical diffusion of horizontal velocity
      & grad           ,& ! gradient of kinetic energy. Unit [m/s]
      & press_grad     ,& ! hydrostatic pressure gradient term. Unit [m/s]
      & cfl_horz       ,& ! horizontal cfl values
!       & zlim           ,& ! zalesak limiter factor
      & vn_bolus  
      
!     onEdges_HalfLevels :: &
!       & w_e            ! vertical velocity at edges. Unit [m/s]

    onVertices :: &
      & vort            ! vorticity at triangle vertices. Unit [1/s]
      
    onVertices_2D :: &
      & vort_50m,    &     ! vorticity at 20 and 70 meters
      & vort_70m

    onVertices_Type(t_cartesian_coordinates) :: &
      & p_vn_dual
    
    onEdges_2D :: &
      & h_e              ,& ! surface height at cell edges. Unit [m].
      & thick_e          ! individual fluid column thickness at edges. Unit [m].
    
    onEdges_2D :: verticallyTotal_mass_flux_e
    
    TYPE(t_ocean_monitor) :: monitor

    ! for the testbed
    !    testbed_div (test_mode: 103-113)
    onCells :: div_model, div_diff, divPtP, divPtP_diff
    
  END TYPE t_hydro_ocean_diag
  !-------------------------------------------------------------------------
    
  !-------------------------------------------------------------------------
  !
  !! auxiliary data
  !
  TYPE t_hydro_ocean_aux

    onEdges :: &
      & g_n           ,& ! explicit velocity term in Adams-Bashford time marching routines,
                                ! at timelevel n
      & g_nm1         ,& ! explicit velocity term in Adams-Bashford time marching routines,
                                ! at timelevel n-1
      & g_nimd        ,& ! explicit velocity term in Adams-Bashford time marching routines,
                                ! located at intermediate timelevel
      & tracer_grad_horz, &      !horizontal tracer gradient                          
      & temperature_grad_horz,&  !horizontal temperature gradient                          
      & salinity_grad_horz       !horizontal salinity gradient                                

    onCells_HalfLevels ::   &
      & tracer_deriv_vert,&      ! vertical tracer gradient
      & temperature_deriv_vert,& ! vertical temperature gradient
      & salinity_deriv_vert      ! vertical salinity gradient            


    onEdges_2D :: &
      & bc_bot_vn       ,& ! normal velocity boundary condition at bottom
      & bc_top_vn       ,& ! normal velocity boundary condition at top
      & bc_top_WindStress ! normal velocity boundary condition at surface

    onCells_2D :: &
      & bc_top_u        ,& ! zonal velocity boundary condition at surface
      & bc_top_v        ,& ! meridional velocity boundary condition at surface
      & bc_top_w        ,& ! vertical velocity boundary condition at surface
      & bc_bot_w        ,&   ! vertical velocity boundary condition at bottom
      & bc_tides_potential, &
      & bc_tides_load,   &
      & bc_total_top_potential, &
      & bc_SAL_potential
            
    onCells_2D_tracers :: &
      & bc_top_tracer,    &
      & bc_bot_tracer 
      
    onCells_2D ::       &
      & p_rhs_sfc_eq     ! right hand side of surface equation
      
    onCells_2D_Type(t_cartesian_coordinates) :: bc_top_veloc_cc
    
!     TYPE(t_pointer_2d_wp), ALLOCATABLE :: bc_top_tracer(:) !< pointer array: one pointer for each tracer boundary condition

    onCells_Type(t_cartesian_coordinates) :: &
      & slopes              ! neutral slopes at cell center in cartesian coordinates

    onCells :: &
      & slopes_squared,   &
      & slopes_drdz,   &
      & slopes_drdx,   & 
      & taper_function_1, &
      & taper_function_2, &
      & diagnose_Redi_flux_vert
      
    onCells_Type(t_cartesian_coordinates) :: &
      & PgradTemperature_horz_center,        & ! reconstructed temperature gradient at cell center in cartesian coordinates
      & PgradSalinity_horz_center,           & ! reconstructed salinity gradient at cell center in cartesian coordinates
      & PgradDensity_horz_center,            & ! reconstructed density gradient at cell center in cartesian coordinates      
      & diagnose_Redi_flux_temp,             & ! diagnostics used for balance of density fluxes in GM,   
      & diagnose_Redi_flux_sal,              & ! can be removed once diagnostic is no longer necessary
      & PgradTracer_horz_center                ! reconstructed tracer gradient at cell center in cartesian coordinates, for HAMMOC in GMredi
      
   onCells ::         &
      & DerivTemperature_vert_center, &  
      & DerivSalinity_vert_center,    &
      & DerivDensity_vert_center,     & 
      & DerivTracer_vert_center
      
  END TYPE t_hydro_ocean_aux

  !-------------------------------
  INTEGER,PARAMETER :: max_tracers = 50
  TYPE t_oce_config
    CHARACTER(LEN=max_char_length) :: tracer_longnames(max_tracers)
    CHARACTER(LEN=max_char_length) :: tracer_stdnames(max_tracers)
    CHARACTER(LEN=max_char_length) :: tracer_shortnames(max_tracers)
    CHARACTER(LEN=max_char_length) :: tracer_units(max_tracers)
    INTEGER                        :: tracer_codes(max_tracers)
  END TYPE t_oce_config
  


  !-------------------------------------------------------------------------------
  TYPE t_operator_coeff

    ! 1) precomputed 3D-factors for mathematical operators (for efficiency).
    !------------------------------------------------------------------------------
    mapEdgesToCells    :: div_coeff
    mapEdgesToVertices :: rot_coeff  
    onEdges            :: grad_coeff ! this should be revised
    mapCellsToEdges_2D :: averageCellsToEdges
    
    !2) Required for description of boundary around a vertex
    !------------------------------------------------------------------------------
    onVertices_3D_Int :: bnd_edges_per_vertex
    onVertices_3D_Connectivity :: vertex_bnd_edge_idx  !(nproma,nlev,nblks_v,1:NO_DUAL_EDGES-2)
    onVertices_3D_Connectivity :: vertex_bnd_edge_blk  !(nproma,nlev,nblks_v,1:NO_DUAL_EDGES-2)
    onVertices_3D_Connectivity :: boundaryEdge_Coefficient_Index   ! this is an index to the rot_coeff  (nproma,nlev,nblks_v,1:NO_DUAL_EDGES-2)
!     REAL(wp),POINTER :: orientation(:,:,:,:)!(nproma,nlev,nblks_v,1:NO_DUAL_EDGES-2)

    ! this is a edge-to-one-cell pointer, needs to be rethinked
!     onEdges_3D_Int :: upwind_cell_idx
!     onEdges_3D_Int :: upwind_cell_blk

    !3) Scalarproduct: The following arrays are required for the reconstruction process.
    !------------------------------------------------------------------------------
    !
    ! Vector pointing from cell circumcenter to edge midpoint. In the associated
    ! cell2edge_weight-array the cell2edge_vec is multiplied by some other geometric
    ! quantities (edge-length, cell area). The weight is used in the reconstruction
    ! the vector is used in the transposed reconstruction.
    ! index=1,nproma, index2=1,nblks_c, index3=1,3
    ! other choice would be index2=1,nblks_e, index3=1,2
    ! Eventually switch to other second indexing if this is more appropriate

    ! Vector pointing from vertex (dual center) to midpoint of dual edge
    ! (/= midpoint of primal edge).
    ! In the associated vertex2dualedge_mid_weight-array the vertex2dualedge_mid_vec
    ! is multiplied by some other geometric quantities (dual edge-length, dual cell
    ! area). The weight is used in the reconstruction the vector is used in the
    ! transposed reconstruction.
    ! index=1,nproma, index2=1,nblks_v, index3=1,6
    ! other choice index2=1,nblks_e, index3=1,2
    ! Eventually switch to other second indexing if this is more appropriate
    ! new constructs for mimetic core:
    mapEdgesToEdges_3D                         :: edge2edge_viacell_coeff

! Note: the following have the connectivity at the first index
    mapEdgesToEdges_2D                         :: edge2edge_viacell_coeff_top       ! the same as the top edge2edge_viacell_coeff
    mapEdgesToEdges_2D                         :: edge2edge_viacell_coeff_integrated! the other levels integrated
    mapEdgesToEdges_2D                         :: edge2edge_viacell_coeff_all       ! all the levels integrated

    mapCellsToCells_2D                         :: lhs_all                ! the left hand side operator coefficients of the height solver
    onCells_2D_Connectivity                    :: lhs_CellToCell_index   ! connectivity of the above
    onCells_2D_Connectivity                    :: lhs_CellToCell_block   ! connectivity of the above
! End Note
 
   mapEdgesToEdges                            :: edge2edge_viavert_coeff

    !coefficient for surface layer, changes in time, in contrast to other coefficients
!     TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2cell_coeff_cc_dyn(:,:,:,:)
    !TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2vert_coeff_cc_dyn(:,:,:,:)

    mapEdgesToCells_3D_Type(t_cartesian_coordinates)    :: edge2cell_coeff_cc
    mapEdgesToCells_3D_Type(t_cartesian_coordinates)    :: edge2cell_coeff_cc_t
    mapEdgesToVertices_3D_Type(t_cartesian_coordinates) :: edge2vert_coeff_cc
    mapEdgesToVertices_3D_Type(t_cartesian_coordinates) :: edge2vert_coeff_cc_t
    mapEdgesToVertices_3D_Type(t_cartesian_coordinates) :: edge2vert_vector_cc

    onCells :: fixed_vol_norm

    blockList_Type(t_verticalAdvection_ppm_coefficients) :: verticalAdvectionPPMcoeffs

!    REAL(wp), POINTER         :: matrix_vert_diff_c(:,:,:,:)
!    REAL(wp), POINTER         :: matrix_vert_diff_e(:,:,:,:)
!    TYPE(t_onCells_Pointer_3d_wp),ALLOCATABLE :: matrix_vert_diff_c_ptr(:)
!    TYPE(t_onCells_Pointer_3d_wp),ALLOCATABLE :: matrix_vert_diff_e_ptr(:)
    onEdges_3D_Int :: edges_SeaBoundaryLevel ! boundary level based on cells:
                                             ! 1=land, 0=boundary between land and sea, -1=between a boundary sea cell and and sea-cell,...,-99999
    onCells_3D_Int :: cells_SeaBoundaryLevel ! as above

  END TYPE t_operator_coeff
    
  
  TYPE t_solverCoeff_singlePrecision
    ! the same as in t_operator_coeff in single precision for using in the solver
    onEdges_2D_RealPrecision(sp) :: grad_coeff                  ! as in t_operator_coeff for the 1st level
    mapEdgesToCells_2D_RealPrecision(sp) :: div_coeff                   ! as in t_operator_coeff for the 1st level

    mapEdgesToEdges_2D_RealPrecision(sp) :: edge2edge_viacell_coeff_all  ! as in t_operator_coeff

    onEdges_2D_RealPrecision(sp) :: edge_thickness                ! as t_hydro_ocean_diag thick_e
    onCells_2D_RealPrecision(sp) :: cell_thickness                ! as t_hydro_ocean_diag thick_c

  END TYPE t_solverCoeff_singlePrecision

  !-----------------------------------------------------------
  ! array of states
  TYPE t_hydro_ocean_state

    TYPE(t_patch_3d), POINTER :: patch_3D
    TYPE(t_hydro_ocean_prog), POINTER :: p_prog(:)    ! time array of prognostic states at different time levels
    TYPE(t_hydro_ocean_diag) :: p_diag
    TYPE(t_hydro_ocean_aux)  :: p_aux
    TYPE(t_operator_coeff), POINTER :: operator_coeff
    TYPE(t_ocean_transport_state)  :: transport_state

  END TYPE t_hydro_ocean_state
  
END MODULE mo_ocean_types

