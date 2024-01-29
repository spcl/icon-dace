! Data type definition for land surface scheme (TERRA)
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

MODULE mo_nwp_lnd_types

  USE mo_kind,                 ONLY: wp
  USE mo_fortran_tools,        ONLY: t_ptr_2d3d
  USE mo_var_list,             ONLY: t_var_list_ptr

  IMPLICIT NONE
  PRIVATE

  !
  !variables
  PUBLIC :: t_lnd_state  !> state vector for land scheme
  PUBLIC :: t_lnd_prog   !!  for prognostic variables
  PUBLIC :: t_wtr_prog   !!  for prognostic variables related to lake and sea ice models
  PUBLIC :: t_lnd_diag   !!  for diagnostic variables

  !
  ! prognostic variables state vector (land)
  !
  TYPE t_lnd_prog

    REAL(wp), POINTER, CONTIGUOUS :: &
    &  t_s_t          (:,:,:)   , & ! temperature of the ground surface             (  K  )
    &  t_sk_t         (:,:,:)   , & ! skin temperature                              (  K  )
    &  t_g            (:,:)     , & ! weighted surface temperature                  (  K  )
    &  t_g_t          (:,:,:)   , & ! weighted surface temperature on tiles         (  K  )
    &  w_i_t          (:,:,:)   , & ! water content of interception water           (m H2O)
    &  w_p_t          (:,:,:)   , & ! water content of pond interception water      (m H2O)
    &  w_s_t          (:,:,:)   , & ! water content of interception snow            (m H2O)
    &  t_so_t         (:,:,:,:) , & ! soil temperature (main level)                 (  K  )
    &  w_so_t         (:,:,:,:) , & ! total water content (ice + liquid water)      (m H20)
    &  w_so_ice_t     (:,:,:,:) , & ! ice content                                   (m H20)
    &  t_snow_t       (:,:,:)   , & ! temperature of the snow-surface               (  K  )
    &  w_snow_t       (:,:,:)   , & ! water content of snow                         (m H2O)
    &  rho_snow_t     (:,:,:)   , & ! snow density                                  (kg/m**3)
    &  t_snow_mult_t  (:,:,:,:) , & ! temperature of snow                           (  K  )
    &  wtot_snow_t    (:,:,:,:) , & ! total (liquid + solid) water content of snow  (m H2O)
    &  wliq_snow_t    (:,:,:,:) , & ! liquid water content in the snow              (m H2O)
    &  rho_snow_mult_t(:,:,:,:) , & ! snow density                                  (kg/m**3)
    &  dzh_snow_t     (:,:,:,:)     ! layer thickness between half levels in snow   (  m  )

    TYPE(t_ptr_2d3d), ALLOCATABLE :: t_snow_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: t_snow_mult_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: t_s_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: t_sk_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: t_gt_ptr(:) 
    TYPE(t_ptr_2d3d), ALLOCATABLE :: w_snow_ptr(:) 
    TYPE(t_ptr_2d3d), ALLOCATABLE :: rho_snow_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: rho_snow_mult_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: w_i_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: w_p_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: w_s_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: t_so_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: w_so_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: w_so_ice_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: wliq_snow_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: wtot_snow_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: dzh_snow_ptr(:)

  END TYPE t_lnd_prog


  !
  ! prognostic variables state vector (sea)
  !
  TYPE t_wtr_prog

    REAL(wp), POINTER, CONTIGUOUS :: &
    &  t_ice        (:,:)   , & ! temperature of the sea ice             (  K  )
    &  h_ice        (:,:)   , & ! depth of the sea ice                   (  m  )
    &  t_snow_si    (:,:)   , & ! temperature of the snow on sea ice     (  K  )
    &  h_snow_si    (:,:)   , & ! depth of the snow on sea ice           (  m  )
    &  alb_si       (:,:)   , & ! prognostic (diffuse) sea-ice albedo    (  -  )
    &  t_snow_lk    (:,:)   , & ! temperature of snow on lake ice        (  K  )
    &  h_snow_lk    (:,:)   , & ! depth of snow on lake ice              (  K  )
    &  t_mnw_lk     (:,:)   , & ! mean temperature of the water column   (  K  )
    &  t_wml_lk     (:,:)   , & ! temperature of the mixed-layer         (  K  )
    &  h_ml_lk      (:,:)   , & ! thickness of the mixed-layer           (  m  )
    &  t_bot_lk     (:,:)   , & ! temperature at the water-bottom sediment interface  (  K  ) 
    &  c_t_lk       (:,:)   , & ! shape factor with respect to the       (  -  )
                                ! temperature profile in lake thermocline
    &  t_b1_lk      (:,:)   , & ! temperature at the bottom of the       (  K  ) 
                                ! upper layer of the sediments 
    &  h_b1_lk      (:,:)       ! thickness of the upper layer of the    (  K  )
                                ! sediments
  END TYPE t_wtr_prog

  !
  ! diagnostic variables state vector
  !
  TYPE t_lnd_diag

    REAL(wp), POINTER, CONTIGUOUS :: &
    &  qv_s         (:,:)   , & ! specific humidity at the surface              (kg/kg)
    &  t_s          (:,:)   , & ! temperature of the ground surface             (  K  )
    &  t_sk         (:,:)   , & ! skin temperature                              (  K  )
    &  t_seasfc     (:,:)   , & ! temperature of the sea surface                (  K  )
    &  w_i          (:,:)   , & ! water content of interception water           (m H2O)
    &  w_p          (:,:)   , & ! water content of pond interception water      (m H2O)
    &  w_s          (:,:)   , & ! water content of interception snow            (m H2O)
    &  t_so         (:,:,:) , & ! soil temperature (main level)                 (  K  )
    &  w_so         (:,:,:) , & ! total water content (ice + liquid water)      (m H20)
    &  w_so_ice     (:,:,:) , & ! ice content                                   (m H20)
    &  runoff_s     (:,:)   , & ! surface water runoff; sum over forecast       (kg/m2)
    &  runoff_g     (:,:)   , & ! soil water runoff; sum over forecast          (kg/m2)
    &  resid_wso    (:,:)   , & ! residuum of soil water budget; sum over fcst  (kg/m2)
    &  fr_seaice    (:,:)   , & ! fraction of sea ice                           ( )
                                ! as partition of total area of the
                                ! grid element, but set to 0 or 1
                                ! index1=1,nproma, index2=1,nblks_c
    &  condhf_ice   (:,:)   , & ! conductive heat flux at sea-ice bottom        ( W/m2)
    &  qv_s_t       (:,:,:) , & ! specific humidity at the surface              (kg/kg)
    &  runoff_s_t   (:,:,:) , & ! surface water runoff; sum over forecast       (kg/m2)
    &  runoff_g_t   (:,:,:) , & ! soil water runoff; sum over forecast          (kg/m2)
    &  resid_wso_t  (:,:,:) , & ! residuum of soil water budget; sum over fcst  (kg/m2)
    &  runoff_s_inst_t(:,:,:),& ! surface water runoff; instantaneous value     (kg/m2)
    &  runoff_g_inst_t(:,:,:),& ! soil water runoff; instantaneous value        (kg/m2)
    &  resid_wso_inst_t(:,:,:),&! residuum of soil water; instantaneous value   (kg/m2)
    &  rstom        (:,:)   , & ! stomatal resistance                           ( s/m )
    &  rstom_t      (:,:,:) , & ! tile based stomatal resistance                ( s/m )
    &  plantevap    (:,:)   , & ! integral function of plant evaporation        (kg/m2)
    &  plantevap_t  (:,:,:) , & ! tile-based integral function of plant evaporation (kg/m2)
    &  t_snow       (:,:)   , & ! temperature of the snow-surface               (  K  )
    &  rho_snow     (:,:)   , & ! snow density                                  (kg/m**3)
    &  w_snow       (:,:)   , & ! snow water equivalent                         (m H2O)
    &  h_snow       (:,:)   , & ! snow height                                   (  m  ) 
    &  h_snow_t     (:,:,:) , & ! snow height                                   (  m  ) 
    &  freshsnow    (:,:)   , & ! indicator for age of snow in top of snow layer(  -  )
    &  freshsnow_t  (:,:,:) , & ! indicator for age of snow in top of snow layer(  -  )
    &  snow_age     (:,:)   , & ! duration of snow cover period                 (  d  )
    &  hsnow_max    (:,:)   , & ! maximum snow depth during current snow-cover period (  m  )
    &  snowfrac     (:,:)   , & ! snow-cover fraction                           (  -  )
    &  snowfrac_lc  (:,:)   , & ! snow-cover fraction                           (  -  )
    &  snowfrac_t   (:,:,:) , & ! snow-cover fraction                           (  -  )
    &  snowfrac_lc_t(:,:,:) , & ! snow-cover fraction per land-cover class      (  -  )
    &  snowfrac_lcu_t(:,:,:), & ! snow-cover fraction per land-cover class      (  -  )
                                ! (see mo_nwp_lnd_state for further explanation of the snow-cover variables)
    &  t_snow_mult  (:,:,:) , & ! temperature of snow                           (  K  )
    &  rho_snow_mult(:,:,:) , & ! snow density                                  (kg/m**3)
    &  wliq_snow    (:,:,:) , & ! liquid water content in the snow              (m H2O)
    &  wtot_snow    (:,:,:) , & ! total (liquid + solid) water content of snow  (m H2O)
    &  dzh_snow     (:,:,:)     ! layer thickness between half levels in snow   (  m  )

    !> (Optional:) Additional diagnostic fields:
    REAL(wp), POINTER ::  &
    &  smi             (:,:,:) , & ! soil moisture index                           (  -  )
    &  snow_melt       (:,:)   , & ! snow melt                                     (kg/m2)
    &  snow_melt_flux_t(:,:,:)     ! snow melt flux tiled                          (kg/m2/s)

    TYPE(t_ptr_2d3d), ALLOCATABLE :: qv_st_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: h_snow_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: plantevap_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: freshsnow_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: snowfrac_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: snowfrac_lc_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: snowfrac_lcu_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: runoff_s_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: runoff_g_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: resid_wso_ptr(:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: snow_melt_ptr(:)

  END TYPE t_lnd_diag


! complete state vector type
!
  TYPE t_lnd_state
    TYPE(t_lnd_prog), ALLOCATABLE  :: prog_lnd(:)          ! number of time levels
    TYPE(t_wtr_prog), ALLOCATABLE  :: prog_wtr(:)          ! number of time levels
    TYPE(t_var_list_ptr), ALLOCATABLE  :: lnd_prog_nwp_list(:) ! number of time levels
    TYPE(t_var_list_ptr), ALLOCATABLE  :: wtr_prog_nwp_list(:) ! number of time levels
    TYPE(t_lnd_diag)               :: diag_lnd
    TYPE(t_var_list_ptr)               :: lnd_diag_nwp_list
  END TYPE t_lnd_state
 

END MODULE mo_nwp_lnd_types
