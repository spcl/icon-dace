!> vegetation variables init (QUINCY)
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
!>#### initialization of vegetation memory variables using, e.g., ic & bc input files
!>
MODULE mo_veg_init
#ifndef __NO_JSBACH__


  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: message, finish, message_text
  USE mo_jsb_control,             ONLY: debug_on

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: veg_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_init'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Intialize vegetation process (after memory has been set up)
  !!
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE veg_init(tile)

    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_class,             ONLY: get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_process_class,     ONLY: VEG_, PHENO_
    USE mo_veg_constants,         ONLY: itree

    ! Use of process memories
    dsl4jsb_Use_memory(PHENO_)
    dsl4jsb_Use_memory(VEG_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER   :: model
    TYPE(t_lctlib_element), POINTER   :: lctlib     !< land-cover-type library - parameter across pft's
    !TYPE(t_input_file)                :: input_file
    CHARACTER(len=*), PARAMETER :: routine = modname//':veg_init'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(VEG_)
    ! Declare pointers to variables in memory
    ! PHENO_
    dsl4jsb_Real2D_onDomain    :: lai
    ! VEG_ 2D
    dsl4jsb_Real2D_onDomain    :: beta_sinklim_ps
    dsl4jsb_Real2D_onDomain    :: t_jmax_opt_mavg
    dsl4jsb_Real2D_onDomain    :: growth_req_n_tlabile_mavg
    dsl4jsb_Real2D_onDomain    :: growth_req_p_tlabile_mavg
    dsl4jsb_Real2D_onDomain    :: t_air_tacclim_mavg
    dsl4jsb_Real2D_onDomain    :: t_air_tphen_mavg
    dsl4jsb_Real2D_onDomain    :: t_air_week_mavg
    dsl4jsb_Real2D_onDomain    :: t_air_month_mavg
    dsl4jsb_Real2D_onDomain    :: veg_pool_leaf_carbon
    dsl4jsb_Real2D_onDomain    :: veg_growth_leaf_carbon
    dsl4jsb_Real2D_onDomain    :: veg_litterfall_leaf_carbon
    dsl4jsb_Real2D_onDomain    :: veg_pool_leaf_nitrogen
    dsl4jsb_Real2D_onDomain    :: veg_pool_leaf_phosphorus
    dsl4jsb_Real2D_onDomain    :: height
    ! VEG_ 3D
    dsl4jsb_Real3D_onDomain    :: leaf_nitrogen_cl
    dsl4jsb_Real3D_onDomain    :: lai_cl
    dsl4jsb_Real3D_onDomain    :: fleaf_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onDomain    :: fn_oth_cl
    ! ---------------------------
    ! 0.4 Debug Option
    IF (debug_on()) CALL message(TRIM(routine), 'Setting initial conditions of vegetation memory for tile '//TRIM(tile%name))
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    IF (tile%lcts(1)%lib_id /= 0) THEN            ! work with lctlib only if the present tile is a pft
      lctlib => model%lctlib(tile%lcts(1)%lib_id)
    ENDIF
    ! Get process config
    !
    ! Get process memories
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(VEG_)
    ! Get process variables (Set pointers to variables in memory)
    ! PHENO_
    dsl4jsb_Get_var2D_onDomain(PHENO_, lai)                           ! in
    ! VEG_ 2D
    dsl4jsb_Get_var2D_onDomain(VEG_, beta_sinklim_ps)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_jmax_opt_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, growth_req_n_tlabile_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, growth_req_p_tlabile_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_air_tacclim_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_air_tphen_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_air_week_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, t_air_month_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_, veg_pool_leaf_carbon)
    dsl4jsb_Get_var2D_onDomain(VEG_, veg_growth_leaf_carbon)
    dsl4jsb_Get_var2D_onDomain(VEG_, veg_litterfall_leaf_carbon)
    dsl4jsb_Get_var2D_onDomain(VEG_, veg_pool_leaf_nitrogen)
    dsl4jsb_Get_var2D_onDomain(VEG_, veg_pool_leaf_phosphorus)
    dsl4jsb_Get_var2D_onDomain(VEG_, height)
    ! VEG_ 3D
    dsl4jsb_Get_var3D_onDomain(VEG_, leaf_nitrogen_cl)
    dsl4jsb_Get_var3D_onDomain(VEG_, lai_cl)
    dsl4jsb_Get_var3D_onDomain(VEG_, fleaf_sunlit_tfrac_mavg_cl)
    dsl4jsb_Get_var3D_onDomain(VEG_, fn_oth_cl)
   
    
    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    beta_sinklim_ps(:,:)                = 1.0_wp
    leaf_nitrogen_cl(:,:,:)             = 0.0_wp        ! this is intentionally wrong; i.e., it will not work at the 1st timestep (which is probably night anyway)
    lai_cl(:,:,:)                       = 0.0_wp        ! this is intentionally wrong; i.e., it will not work at the 1st timestep (which is probably night anyway)
    fn_oth_cl(:,:,:)                    = 1.0_wp        ! it is necessary to start with fn_oth = 1.0 leaving the rest of the fraction 0.0
                                                        ! fn_chl_cl, fn_rub_cl, fn_et_cl, fn_pepc_cl are set to 0.0_wp by default with the add_var()    
    ! moving averages
    t_jmax_opt_mavg(:,:)                = 17.0_wp
    growth_req_n_tlabile_mavg(:,:)      = 1._wp / 25._wp
    growth_req_p_tlabile_mavg(:,:)      = 1._wp / 14._wp
    fleaf_sunlit_tfrac_mavg_cl(:,:,:)   = 0.5_wp
    t_air_tacclim_mavg(:,:)             = 283.15_wp     ! this init value may equal t_acclim_zero (283.15) to make work: mo_veg_respiration:temperature_response_respiration    
    t_air_tphen_mavg(:,:)               = 283.15_wp
    t_air_week_mavg(:,:)                = 283.15_wp
    t_air_month_mavg(:,:)               = 283.15_wp

    ! 2D var used for vegetation Pool/Flux replacement
    veg_pool_leaf_carbon(:,:)           = 0.0_wp
    veg_growth_leaf_carbon(:,:)         = 0.0_wp
    veg_litterfall_leaf_carbon(:,:)     = 0.0_wp
        
    IF (tile%lcts(1)%lib_id /= 0) THEN      ! work with lctlib only if the present tile is a pft
      veg_pool_leaf_carbon(:,:)     = lai(:,:) / lctlib%sla
      veg_pool_leaf_nitrogen(:,:)   = lai(:,:) / lctlib%sla / lctlib%cn_leaf
      veg_pool_leaf_phosphorus(:,:) = veg_pool_leaf_nitrogen(:,:) / lctlib%np_leaf

      IF(lctlib%growthform == itree)THEN
         height(:,:)                = 10.0_wp
      ELSE 
         height(:,:)                = 1.0_wp
      ENDIF
    ENDIF

  END SUBROUTINE veg_init

#endif
END MODULE mo_veg_init
