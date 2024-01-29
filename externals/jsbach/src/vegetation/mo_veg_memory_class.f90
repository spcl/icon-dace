!> vegetation process memory (QUINCY)
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### definition and init of (memory) variables for the vegetation process
!>
MODULE mo_veg_memory_class

  USE mo_kind,                   ONLY: wp
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_impl_constants,     ONLY: false
  USE mo_exception,              ONLY: message
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: VEG_TYPE, LAND_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_veg_memory, max_no_of_vars
 
  INTEGER, PARAMETER :: max_no_of_vars      = 250  !< TYPE t_jsb_var


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for veg memory
  !!
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_veg_memory

      ! basic plant 2D
      TYPE(t_jsb_var_real2d)          :: &   
                                    height            , &            !< plant height [m]
                                    ! diameter          , &            !< plant diameter [m]
                                    ! dens_ind          , &            !< density of individuals [# m-2]
                                    ! lai               , &            !< actual leaf area index [m2 m-2]
                                    ! target_cn_leaf    , &            !< target foliar C:N ratio [mol C mol-1 N]
                                    ! target_np_leaf    , &            !< target foliar N:P ratio [mol N mol-1 P]
                                    ! mean_leaf_age     , &            !< average foliage age in days [days]
                                    ! cohort_age        , &            !< age of the cohort in years, used with veg_dynamics_scheme cohort [years]
                                    ! f_n_demand        , &            !< relative plant demand for N
                                    ! f_p_demand        , &            !< relative plant demand for P
                                    ! beta_sinklim_ps   , &            !< scaling factor to account for a direct sink limitation contraint on photosynthesis [unitless]
                                    ! dphi              , &            !< plant-soil water potential gradient [MPa]
                                    growth_req_n      , &            !< moles N required for a unit of C growth under current allocation fractions [mol N m-2 -1]
                                    growth_req_p      !, &            !< moles P required for a unit of N growth under current allocation fractions [mol P m-2 -1]
                                    ! leaf_cn_direction , &            !< indicates direction of change in leaf CN ratio to maximise NPP [unitless]
                                    ! target_lai        , &            !< LAI implied from sap wood area [m2 m-2]
                                    ! kstar_labile                     !< labile pool turnover rate assuming only temperature and moisture constraints, no stoichiometry [day-1]
      ! 1.1 memory of conditions for the labile pool dynamics
      TYPE(t_jsb_var_real2d)          :: &   
                                    gpp_tlabile_mavg               , &  !< gpp  [µmol m-2 s-1]
                                    maint_respiration_tlabile_mavg , &  !< maintenance respiration [µmol m-2 s-1]
                                    growth_req_n_tlabile_mavg      , &  !< moles N required for a unit of C growth under current allocation fractions time-averaged [mol N m-2 s-1]
                                    growth_req_p_tlabile_mavg           !< moles P required for a unit of N growth under current allocation fractions time-averaged [mol P m-2 s-1]
      TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_tphen_mavg               , &  !< air temperature time-averaged [K]
                                    t_air_week_mavg                , &  !< air temperature time-averaged [K]
                                    t_air_month_mavg               , &  !< air temperature time-averaged [K]
                                    mean_leaf_age                       !< average foliage age in days [days]
      ! various var needed in update_canopy_fluxes()
      TYPE(t_jsb_var_real2d)          :: &
                                    beta_sinklim_ps, &
                                    t_jmax_opt_mavg, &
                                    t_air_tacclim_mavg
      ! fluxes
      TYPE(t_jsb_var_real2d)          :: &   
                                    maint_respiration_pot         !, & !< potential maintenance respiration [µmol CO2 m-2 s-1]
                                    ! maint_respiration             , & !< maintenance respiration [µmol CO2 m-2 s-1]
                                    ! maint_respiration_c13         , & !< maintenance respiration [µmol 13CO2 m-2 s-1]
                                    ! maint_respiration_c14         , & !< maintenance respiration [µmol 14CO2 m-2 s-1]
                                    ! growth_respiration            , & !< growth respiration [µmol CO2 m-2 s-1]
                                    ! growth_respiration_c13        , & !< growth respiration [µmol 13CO2 m-2 s-1]
                                    ! growth_respiration_c14        , & !< growth respiration [µmol 14CO2 m-2 s-1]
                                    ! n_transform_respiration       , & !< respiration associated with N transformation [µmol CO2 m-2 s-1]
                                    ! n_fixation_respiration        , & !< respiration associated with N fixation [µmol CO2 m-2 s-1]
                                    ! n_processing_respiration      , & !< respiration associated with N processing [µmol CO2 m-2 s-1]
                                    ! n_processing_respiration_c13  , & !< respiration associated with N processing [µmol 13CO2 m-2 s-1]
                                    ! n_processing_respiration_c14  , & !< respiration associated with N processing [µmol 14CO2 m-2 s-1]
                                    ! npp                           , & !< net primary production [µmol CO2 m-2 s-1]
                                    ! npp_c13                       , & !< net primary production [µmol 13CO2 m-2 s-1]
                                    ! npp_c14                       , & !< net primary production [µmol 14CO2 m-2 s-1]
                                    ! uptake_nh4                    , & !< NH4 uptake from soil [µmol N m-2 s-1]
                                    ! uptake_nh4_n15                , & !< 15NH4 uptake from soil [µmol N m-2 s-1]
                                    ! uptake_no3                    , & !< NO3 uptake from soil [µmol N m-2 s-1]
                                    ! uptake_no3_n15                , & !< 15NO3 uptake from soil [µmol N m-2 s-1]
                                    ! uptake_po4                    , & !< PO4 uptake from soil [µmol P m-2 s-1]
                                    ! n_fixation                    , & !< symbiontic N fixation [µmol N m-2 s-1]
                                    ! n_fixation_n15                , & !< symbiontic N fixation [µmol 15N m-2 s-1]
                                    ! recycling_leaf_n              , & !< net flux of N from leaf to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    ! recycling_leaf_n15            , & !< net flux of 15N from leaf to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    ! recycling_leaf_p              , & !< net flux of P from leaf to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    ! recycling_fine_root_n         , & !< net flux of N from fine roots to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    ! recycling_fine_root_n15       , & !< net flux of 15N from fine roots to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    ! recycling_fine_root_p         , & !< net flux of P from fine roots to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    ! recycling_heart_wood_n        , & !< net flux of N from sapwood to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    ! recycling_heart_wood_n15      , & !< net flux of 15N from sapwood to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    ! recycling_heart_wood_p        , & !< net flux of P from sapwood to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    ! net_growth                        !< npp minus turnover [µmol CO2 m-2 s-1]
      TYPE(t_jsb_var_real3d)          :: &
                                    fleaf_sunlit_tfrac_mavg_cl
      ! canopy layers
      TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    leaf_nitrogen_cl  , &           !< N in canopy layer [mmol m-2 LAI]
                                    fn_chl_cl         , &           !< fraction of N in chlorophyll, pigments only [unitless]
                                    fn_et_cl          , &           !< fraction of N in electron transport [unitless]
                                    fn_rub_cl         , &           !< fraction of N in Rubisco [unitless]
                                    fn_pepc_cl        , &           !< fraction of N in PEP carboylase, C4 plants only [unitless]
                                    fn_oth_cl         , &           !< fraction of N not associated with photosynthesis [unitless]
                                    fleaf_sunlit_cl   , &           !< fraction of layer which is sunlit [unitless]
                                    lai_cl            , &           !< leaf area in the canopy layer [m2 m-2]
                                    cumm_lai_cl                     !< cummulative leaf area above the mid-point of the canopy layer [m2 m-2]
      ! 2D var used for vegetation Pool/Flux replacement
      TYPE(t_jsb_var_real2d)          :: &   
                                    veg_pool_leaf_carbon          , &
                                    veg_pool_leaf_nitrogen        , &
                                    veg_pool_leaf_phosphorus      , &
                                    veg_growth_leaf_carbon        , &
                                    veg_litterfall_leaf_carbon

  CONTAINS
    PROCEDURE :: Init => Init_veg_memory
  END TYPE t_veg_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_memory_class'


CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> initialize memory for the VEG_ process
  !!
  !!  
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_veg_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1.InOut
    CLASS(t_veg_memory),  INTENT(inout), TARGET  :: mem             !< veg memory
    CHARACTER(len=*),     INTENT(in)             :: prefix          !< process name
    CHARACTER(len=*),     INTENT(in)             :: suffix          !< tile name
    INTEGER,              INTENT(in)             :: lct_ids(:)      !< vector of land cover type IDs
    INTEGER,              INTENT(in)             :: model_id        !< model ID model\%id
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_surface                ! Vertical grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_canopy                 ! Vertical grid
    INTEGER                    :: table
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_veg_memory'
    ! ---------------------------
    ! 0.5 Get Memory
    model         => Get_model(model_id)
    table         =  tables(1)
    hgrid         => Get_grid(model%grid_id)
    vgrid_surface => Get_vgrid('surface')
    vgrid_canopy  => Get_vgrid('canopy_layer_q_')



    !> add 2D & 3D variables
    !!
    IF ( One_of(LAND_TYPE, lct_ids(:)) > 0 .OR. &
       & One_of(VEG_TYPE,  lct_ids(:)) > 0) THEN

      CALL mem%Add_var('height', mem%height, &
            hgrid, vgrid_surface, &
            t_cf('height', 'm', 'plant height'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
            
      CALL mem%Add_var('growth_req_n', mem%growth_req_n, &
            hgrid, vgrid_surface, &
            t_cf('growth_req_n', 'mol N m-2 -1', 'moles N required for a unit of C growth under current allocation fractions '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('growth_req_p', mem%growth_req_p, &
            hgrid, vgrid_surface, &
            t_cf('growth_req_p', 'mol P m-2 -1', 'moles P required for a unit of N growth under current allocation fractions '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
            
      CALL mem%Add_var('gpp_tlabile_mavg', mem%gpp_tlabile_mavg, &
            hgrid, vgrid_surface, &
            t_cf('gpp_tlabile_mavg', 'µmol m-2 s-1', 'gpp '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('maint_respiration_tlabile_mavg', mem%maint_respiration_tlabile_mavg, &
            hgrid, vgrid_surface, &
            t_cf('maint_respiration_tlabile_mavg', 'µmol m-2 s-1', 'maintenance respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('growth_req_n_tlabile_mavg', mem%growth_req_n_tlabile_mavg, &
            hgrid, vgrid_surface, &
            t_cf('growth_req_n_tlabile_mavg', '', 'moles N required for a unit of C growth under current allocation fractions'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('growth_req_p_tlabile_mavg', mem%growth_req_p_tlabile_mavg, &
            hgrid, vgrid_surface, &
            t_cf('growth_req_p_tlabile_mavg', '', 'moles P required for a unit of N growth under current allocation fractions'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('t_air_tphen_mavg', mem%t_air_tphen_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_air_tphen_mavg', '', 'air temp mavg for phenology calc'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('t_air_week_mavg', mem%t_air_week_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_air_week_mavg', '', 'avrg air temperature across 7 days'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('t_air_month_mavg', mem%t_air_month_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_air_month_mavg', '', 'avrg air temperature across 30 days'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('mean_leaf_age', mem%mean_leaf_age, &
            hgrid, vgrid_surface, &
            t_cf('mean_leaf_age', '', 'average foliage age in days '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
            

      CALL mem%Add_var('beta_sinklim_ps', mem%beta_sinklim_ps, &
            hgrid, vgrid_surface, &
            t_cf('beta_sinklim_ps', '', 'scaling factor to account for a direct sink limitation contraint on photosynthesis'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('t_jmax_opt_mavg', mem%t_jmax_opt_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_jmax_opt_mavg', '', 'optimim temperature for electron transport of photosynthesis'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
            
      CALL mem%Add_var('t_air_tacclim_mavg', mem%t_air_tacclim_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_air_tacclim_mavg', '', 'average air temperature for respiration acclimation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('maint_respiration_pot', mem%maint_respiration_pot, &
            hgrid, vgrid_surface, &
            t_cf('maint_respiration_pot', 'µmol CO2 m-2 s-1', 'potential maintenance respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
            
      CALL mem%Add_var('fleaf_sunlit_tfrac_mavg_cl', mem%fleaf_sunlit_tfrac_mavg_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fleaf_sunlit_tfrac_mavg_cl', '', 'fraction of sunlit leaves'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
            
      CALL mem%Add_var('leaf_nitrogen_cl', mem%leaf_nitrogen_cl, &
            hgrid, vgrid_canopy, &
            t_cf('leaf_nitrogen_cl', 'mmol m-2 LAI', 'N in canopy layer'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('fn_chl_cl', mem%fn_chl_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fn_chl_cl', '', 'fraction of N in chlorophyll, pigments only'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('fn_et_cl', mem%fn_et_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fn_et_cl', '', 'fraction of N in electron transport '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('fn_rub_cl', mem%fn_rub_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fn_rub_cl', '', 'fraction of N in Rubisco '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('fn_pepc_cl', mem%fn_pepc_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fn_pepc_cl', '', 'fraction of N in PEP carboylase, C4 plants only'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('fn_oth_cl', mem%fn_oth_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fn_oth_cl', '', 'fraction of N not associated with photosynthesis '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('fleaf_sunlit_cl', mem%fleaf_sunlit_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fleaf_sunlit_cl', '', 'fraction of layer which is sunlit'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('lai_cl', mem%lai_cl, &
            hgrid, vgrid_canopy, &
            t_cf('lai_cl', 'm2 m-2', 'leaf area in the canopy layer '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('cumm_lai_cl', mem%cumm_lai_cl, &
            hgrid, vgrid_canopy, &
            t_cf('cumm_lai_cl', 'm2 m-2', 'cummulative leaf area above the mid-point of the canopy layer '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('veg_pool_leaf_carbon', mem%veg_pool_leaf_carbon, &
            hgrid, vgrid_surface, &
            t_cf('veg_pool_leaf_carbon', '[]', 'leaf carbon veg pool'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('veg_pool_leaf_nitrogen', mem%veg_pool_leaf_nitrogen, &
            hgrid, vgrid_surface, &
            t_cf('veg_pool_leaf_nitrogen', '[]', 'leaf nitrogen veg pool'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
            
      CALL mem%Add_var('veg_pool_leaf_phosphorus', mem%veg_pool_leaf_phosphorus, &
            hgrid, vgrid_surface, &
            t_cf('veg_pool_leaf_phosphorus', '[]', 'leaf phosphorus veg pool'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
            
      CALL mem%Add_var('veg_growth_leaf_carbon', mem%veg_growth_leaf_carbon, &
            hgrid, vgrid_surface, &
            t_cf('veg_growth_leaf_carbon', '[]', 'leaf carbon growth flux'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
            
      CALL mem%Add_var('veg_litterfall_leaf_carbon', mem%veg_litterfall_leaf_carbon, &
            hgrid, vgrid_surface, &
            t_cf('veg_litterfall_leaf_carbon', '[]', 'leaf crabon litterfall flux'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

    ENDIF  ! IF One_of(LAND_TYPE  OR VEG_TYPE, lct_ids)

  END SUBROUTINE Init_veg_memory

END MODULE mo_veg_memory_class
