!
! This module is the interface between nwp_nh_interface to the radiation scheme RRTM.
!
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

!----------------------------

MODULE mo_nwp_rrtm_interface

  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_nwp_tuning_config,    ONLY: tune_dust_abs
  USE mo_grid_config,          ONLY: l_limited_area, nexlevs_rrg_vnest
  USE mo_exception,            ONLY: message, message_text
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_parallel_config,      ONLY: nproma, p_test_run
  USE mo_run_config,           ONLY: msg_level, iqv, iqc, iqi
  USE mo_impl_constants,       ONLY: min_rlcell_int
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c, grf_ovlparea_start_c, grf_fbk_start_c
  USE mo_kind,                 ONLY: wp
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_nwp_lnd_types,        ONLY: t_lnd_prog
  USE mo_model_domain,         ONLY: t_patch, p_patch_local_parent
  USE mo_phys_nest_utilities,  ONLY: t_upscale_fields,upscale_rad_input, downscale_rad_output
  USE mo_nonhydro_types,       ONLY: t_nh_diag
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag
  USE mo_radiation,            ONLY: radiation_nwp
  USE mo_aerosol_util,         ONLY: tune_dust
  USE mo_lrtm_par,             ONLY: nbndlw
  USE mo_sync,                 ONLY: global_max, global_min
  USE mtime,                   ONLY: datetime
  USE mo_fortran_tools,        ONLY: assert_acc_host_only

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nwp_rrtm_interface'


  PUBLIC :: nwp_rrtm_radiation
  PUBLIC :: nwp_rrtm_radiation_reduced


CONTAINS

  !---------------------------------------------------------------------------------------
  SUBROUTINE nwp_rrtm_radiation ( current_date, pt_patch, ext_data,                      &
    &  zaeq1, zaeq2, zaeq3, zaeq4, zaeq5, pt_diag, prm_diag, lnd_prog, lacc )

    CHARACTER(len=*), PARAMETER::  &
      &  routine = modname//'::nwp_rrtm_radiation'

    TYPE(datetime), POINTER, INTENT(in) :: current_date
    
    TYPE(t_patch),        TARGET,INTENT(in) :: pt_patch     !<grid/patch info.
    TYPE(t_external_data),TARGET,INTENT(in) :: ext_data

    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    REAL(wp), INTENT(in) :: &
      & zaeq1(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq2(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq3(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq4(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq5(nproma,pt_patch%nlev,pt_patch%nblks_c)

    TYPE(t_nh_diag), TARGET, INTENT(in)  :: pt_diag     !<the diagnostic variables
    TYPE(t_nwp_phy_diag),       INTENT(inout):: prm_diag
    TYPE(t_lnd_prog),           INTENT(inout):: lnd_prog

    REAL(wp):: aclcov(nproma,pt_patch%nblks_c), dust_tunefac(nproma,nbndlw)

    REAL(wp), DIMENSION(:,:), POINTER :: ptr_clc => NULL ()
    REAL(wp), DIMENSION(:,:), POINTER :: ptr_acdnc => NULL (),    &
            & ptr_reff_qc => NULL(), ptr_reff_qi => NULL()
    REAL(wp), DIMENSION(:),    POINTER :: &
            & ptr_fr_glac => NULL(), ptr_fr_land => NULL()


    ! Local scalars:
    INTEGER:: jb
    INTEGER:: jg                !domain id
    INTEGER:: nlev, nlevp1      !< number of full and half levels

    INTEGER:: rl_start, rl_end
    INTEGER:: i_startblk, i_endblk    !> blocks
    INTEGER:: i_startidx, i_endidx    !< slices
    INTEGER:: i_nchdom                !< domain index


    CALL assert_acc_host_only(routine, lacc)

    i_nchdom  = MAX(1,pt_patch%n_childdom)
    jg        = pt_patch%id

    ! number of vertical levels
    nlev   = pt_patch%nlev
    nlevp1 = pt_patch%nlevp1



    !-------------------------------------------------------------------------
    !> Radiation
    !-------------------------------------------------------------------------

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

    IF (msg_level >= 12) &
      &           CALL message(routine, 'RRTM radiation on full grid')

    IF (p_test_run) THEN
      CALL get_indices_c(pt_patch, i_startblk, i_startblk, i_endblk, &
           &                         i_startidx, i_endidx, rl_start, rl_end)
      aclcov(1:i_startidx-1,i_startblk) = 0
      prm_diag%lwflxall(1:i_startidx-1,:,i_startblk) = 0
      prm_diag%trsolall(1:i_startidx-1,:,i_startblk) = 0
      prm_diag%lwflx_up_sfc_rs(1:i_startidx-1,i_startblk) = 0
      prm_diag%trsol_up_toa(1:i_startidx-1,i_startblk) = 0
      prm_diag%trsol_up_sfc(1:i_startidx-1,i_startblk) = 0
      prm_diag%trsol_par_sfc(1:i_startidx-1,i_startblk) = 0
      prm_diag%trsol_dn_sfc_diff(1:i_startidx-1,i_startblk) = 0
      prm_diag%trsolclr_sfc(1:i_startidx-1,i_startblk) = 0

      IF (atm_phy_nwp_config(jg)%l_3d_rad_fluxes) THEN
        prm_diag%lwflx_up    (1:i_startidx-1,:,i_startblk) = 0
        prm_diag%lwflx_dn    (1:i_startidx-1,:,i_startblk) = 0
        prm_diag%swflx_up    (1:i_startidx-1,:,i_startblk) = 0
        prm_diag%swflx_dn    (1:i_startidx-1,:,i_startblk) = 0
        prm_diag%lwflx_up_clr(1:i_startidx-1,:,i_startblk) = 0
        prm_diag%lwflx_dn_clr(1:i_startidx-1,:,i_startblk) = 0
        prm_diag%swflx_up_clr(1:i_startidx-1,:,i_startblk) = 0
        prm_diag%swflx_dn_clr(1:i_startidx-1,:,i_startblk) = 0
      END IF
    END IF

!$OMP PARALLEL PRIVATE(jb,i_startidx,i_endidx,dust_tunefac,                   &
!$OMP                  ptr_clc,ptr_acdnc,ptr_fr_land,ptr_fr_glac,ptr_reff_qc,ptr_reff_qi) 
!$OMP DO ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        &                         i_startidx, i_endidx, rl_start, rl_end)



      prm_diag%tsfctrad(i_startidx:i_endidx,jb) = lnd_prog%t_g(i_startidx:i_endidx,jb)

      IF (tune_dust_abs > 0._wp) THEN
!DIR$ NOINLINE
        CALL tune_dust(pt_patch%cells%center(:,jb)%lat,pt_patch%cells%center(:,jb)%lon,i_endidx,dust_tunefac)
      ELSE
        dust_tunefac(:,:) = 1._wp
      ENDIF

      NULLIFY(ptr_clc,ptr_acdnc,ptr_fr_land,ptr_fr_glac,ptr_reff_qc,ptr_reff_qi)

      IF (atm_phy_nwp_config(jg)%luse_clc_rad) THEN  ! clc_rad has to be used instead of clc
        ptr_clc => prm_diag%clc_rad(:,:,jb)
      ELSE
        ptr_clc => prm_diag%clc(:,:,jb)
      END IF
      
      IF (atm_phy_nwp_config(jg)%icpl_rad_reff == 0) THEN ! Internal parameterization of reff
        ptr_acdnc  =>  prm_diag%acdnc(:,:,jb)
        ptr_fr_land=>  ext_data%atm%fr_land(:,jb)  !< in     land fraction
        ptr_fr_glac=>  ext_data%atm%fr_glac(:,jb)   !< in     land glacier fraction
      ELSE
        ptr_reff_qc => prm_diag%reff_qc(:,:,jb)
        ptr_reff_qi => prm_diag%reff_qi(:,:,jb)
      ENDIF

      CALL radiation_nwp(               &
                              !
                              ! input
                              ! -----
                              !
        & current_date                     ,&!< in current date
                              ! indices and dimensions
        & jg         =jg                   ,&!< in domain index
        & jb         =jb                   ,&!< in block index
        & icpl_reff  =atm_phy_nwp_config(jg)%icpl_rad_reff,&  !< in option for radiation reff coupling
        & jcs        =i_startidx           ,&!< in  start index for loop over block
        & jce        =i_endidx             ,&!< in  end   index for loop over block
        & kbdim      =nproma               ,&!< in  dimension of block over cells
        & klev       =nlev                 ,&!< in  number of full levels = number of layers
        & klevp1     =nlevp1               ,&!< in  number of half levels = number of layer ifcs
                              !
        & ktype      =prm_diag%ktype(:,jb) ,&!< in     type of convection
                              !
                              ! surface: albedo + temperature
        & zland      =ptr_fr_land   ,&!< in     land fraction
        & zglac      =ptr_fr_glac   ,&!< in     land glacier fraction
                              !
        & cos_mu0    =prm_diag%cosmu0  (:,jb) ,&!< in  cos of zenith angle mu0
        & alb_vis_dir=prm_diag%albvisdir(:,jb) ,&!< in surface albedo for visible range, direct
        & alb_nir_dir=prm_diag%albnirdir(:,jb) ,&!< in surface albedo for near IR range, direct
        & alb_vis_dif=prm_diag%albvisdif(:,jb),&!< in surface albedo for visible range, diffuse
        & alb_nir_dif=prm_diag%albnirdif(:,jb),&!< in surface albedo for near IR range, diffuse
        & emis_rad   =prm_diag%lw_emiss(:,jb),&!< in longwave surface emissivity
        & tk_sfc     =prm_diag%tsfctrad(:,jb) ,&!< in surface temperature
                              !
                              ! atmosphere: pressure, tracer mixing ratios and temperature
        & pp_hl      =pt_diag%pres_ifc  (:,:,jb)     ,&!< in  pres at half levels at t-dt [Pa]
        & pp_fl      =pt_diag%pres      (:,:,jb)     ,&!< in  pres at full levels at t-dt [Pa]
        & tk_fl      =pt_diag%temp      (:,:,jb)     ,&!< in  temperature at full level at t-dt
        & qm_vap     =prm_diag%tot_cld  (:,:,jb,iqv) ,&!< in  water vapor mass mix ratio at t-dt
        & qm_liq     =prm_diag%tot_cld  (:,:,jb,iqc) ,&!< in cloud water mass mix ratio at t-dt
        & qm_ice     =prm_diag%tot_cld  (:,:,jb,iqi) ,&!< in cloud ice mass mixing ratio at t-dt
        & qm_o3      =ext_data%atm%o3   (:,:,jb)     ,&!< in o3 mass mixing ratio at t-dt
        & cdnc       =ptr_acdnc                      ,&!< in  cloud droplet numb conc. [1/m**3]
        & reff_liq   =ptr_reff_qc                    ,&!< in effective radius liquid phase 
        & reff_frz   =ptr_reff_qi                    ,&!< in effective radius frozen phase 
        & cld_frc    =ptr_clc                        ,&!< in  cloud fraction [m2/m2]
        & zaeq1      = zaeq1(:,:,jb)                 ,&!< in aerosol continental
        & zaeq2      = zaeq2(:,:,jb)                 ,&!< in aerosol maritime
        & zaeq3      = zaeq3(:,:,jb)                 ,&!< in aerosol urban
        & zaeq4      = zaeq4(:,:,jb)                 ,&!< in aerosol volcano ashes
        & zaeq5      = zaeq5(:,:,jb)                 ,&!< in aerosol stratospheric background
        & dust_tunefac = dust_tunefac (:,:)          ,&!< in LW tuning factor for dust aerosol
                              ! output
                              ! ------
                              !
        & cld_cvr    =  aclcov             (:,jb),  &     !< out cloud cover in a column [m2/m2]
        & flx_lw_net =  prm_diag%lwflxall(:,:,jb), &      !< out terrestrial flux, all sky, net down
        & trsol_net  =  prm_diag%trsolall(:,:,jb), &      !< out solar transmissivity, all sky, net down
        & flx_uplw_sfc = prm_diag%lwflx_up_sfc_rs(:,jb),& !< out longwave upward flux at surface
        & trsol_up_toa = prm_diag%trsol_up_toa(:,jb), &   !< out upward solar transmissivity at TOA
        & trsol_up_sfc = prm_diag%trsol_up_sfc(:,jb), &   !< out upward solar transmissivity at surface
        & trsol_nir_sfc = prm_diag%trsol_nir_sfc(:,jb), & !< out downward transmissivity for near-infrared rad. at surface
        & trsol_vis_sfc = prm_diag%trsol_vis_sfc(:,jb), & !< out downward transmissivity for visible rad. at surface
        & trsol_par_sfc = prm_diag%trsol_par_sfc(:,jb), & !< out downward transmissivity for photosynthetically active rad. at surface
        & fr_nir_sfc_diffus = prm_diag%fr_nir_sfc_diff(:,jb), & !< out diffuse fraction of downward near-infrared rad. at surface
        & fr_vis_sfc_diffus = prm_diag%fr_vis_sfc_diff(:,jb), & !< out diffuse fraction of downward visible rad. at surface
        & fr_par_sfc_diffus = prm_diag%fr_par_sfc_diff(:,jb), & !< out diffuse fraction of downward photosynthetically active rad. at surface
        & trsol_dn_sfc_diffus = prm_diag%trsol_dn_sfc_diff(:,jb), &  !< out downward diffuse solar transmissivity at surface
        & trsol_clr_sfc = prm_diag%trsolclr_sfc(:,jb), &   !< out clear-sky net transmissvity at surface
        & lwflx_clr_sfc = prm_diag%lwflxclr_sfc(:,jb), &  !< out clear-sky net LW flux at surface
        !optional output: 3D flux output
        &  flx_lw_dn     = prm_diag%lwflx_dn(:,:,jb),     & !< Downward LW flux (all-sky)   [Wm2]
        &  flx_sw_dn     = prm_diag%swflx_dn(:,:,jb),     & !< Downward SW flux (all-sky)   [Wm2]
        &  flx_lw_up     = prm_diag%lwflx_up(:,:,jb),     & !< Upward LW flux   (all-sky)   [Wm2]
        &  flx_sw_up     = prm_diag%swflx_up(:,:,jb),     & !< Upward SW flux   (all-sky)   [Wm2]
        &  flx_lw_dn_clr = prm_diag%lwflx_dn_clr(:,:,jb), & !< Downward LW flux (clear sky) [Wm2]
        &  flx_sw_dn_clr = prm_diag%swflx_dn_clr(:,:,jb), & !< Downward SW flux (clear sky) [Wm2]
        &  flx_lw_up_clr = prm_diag%lwflx_up_clr(:,:,jb), & !< Upward LW flux   (clear sky) [Wm2]
        &  flx_sw_up_clr = prm_diag%swflx_up_clr(:,:,jb) )  !< Upward SW flux   (clear sky) [Wm2]


      ENDDO ! blocks

!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE nwp_rrtm_radiation
  !---------------------------------------------------------------------------------------


  !---------------------------------------------------------------------------------------
  SUBROUTINE nwp_rrtm_radiation_reduced ( current_date, pt_patch, pt_par_patch, ext_data, &
    &                                     zaeq1,zaeq2,zaeq3,zaeq4,zaeq5,    &
    &                                     pt_diag,prm_diag,lnd_prog, lacc )

    CHARACTER(len=*), PARAMETER::  &
      &  routine = modname//'::nwp_rrtm_radiation_reduced'

    TYPE(datetime), POINTER, INTENT(in) :: current_date

    TYPE(t_patch),        TARGET,INTENT(in) :: pt_patch     !<grid/patch info.
    TYPE(t_patch),        TARGET,INTENT(in) :: pt_par_patch !<grid/patch info (parent grid)
    TYPE(t_external_data),INTENT(in):: ext_data
    REAL(wp),             INTENT(in) ::               &
      & zaeq1(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq2(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq3(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq4(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq5(nproma,pt_patch%nlev,pt_patch%nblks_c)

    TYPE(t_nh_diag), TARGET,    INTENT(inout):: pt_diag     !<the diagnostic variables
    TYPE(t_nwp_phy_diag),       INTENT(inout):: prm_diag
    TYPE(t_lnd_prog),           INTENT(inout):: lnd_prog

    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    REAL(wp):: aclcov(nproma,pt_patch%nblks_c), dust_tunefac(nproma,nbndlw)
    ! For radiation on reduced grid
    ! These fields need to be allocatable because they have different dimensions for
    ! the global grid and nested grids, and for runs with/without MPI parallelization
    ! Input fields
    REAL(wp), ALLOCATABLE, TARGET:: zrg_emis_rad (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_cosmu0   (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_albvisdir(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_albnirdir(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_albvisdif(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_albnirdif(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_albdif   (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_tsfc     (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_rtype    (:,:) ! type of convection (integer)
    INTEGER,  ALLOCATABLE, TARGET:: zrg_ktype    (:,:) ! type of convection (real)

    REAL(wp), ALLOCATABLE, TARGET:: zrg_pres_ifc (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zlp_pres_ifc (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_pres     (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_temp     (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_o3       (:,:,:)
    REAL(wp), POINTER            :: zrg_reff_liq (:,:,:) => NULL()
    REAL(wp), POINTER            :: zrg_reff_frz (:,:,:) => NULL()
    REAL(wp), POINTER            :: zrg_extra_flds(:,:,:,:) => NULL()
    REAL(wp), POINTER            :: zrg_extra_2D(:,:,:) => NULL()
    REAL(wp), ALLOCATABLE, TARGET:: zrg_tot_cld  (:,:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zlp_tot_cld  (:,:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_clc      (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aeq1(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aeq2(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aeq3(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aeq4(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aeq5(:,:,:)
    ! Output fields
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aclcov   (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_lwflxall (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_trsolall (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_lwflx_up_sfc   (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_trsol_up_toa   (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_trsol_up_sfc   (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_trsol_nir_sfc  (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_trsol_vis_sfc  (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_trsol_par_sfc  (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_fr_nir_sfc_diff(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_fr_vis_sfc_diff(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_fr_par_sfc_diff(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_trsol_dn_sfc_diff(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_trsol_clr_sfc   (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_lwflx_clr_sfc   (:,:)

    REAL(wp), ALLOCATABLE, TARGET:: zrg_lwflx_up    (:,:,:)    !< longwave  3D upward   flux          
    REAL(wp), ALLOCATABLE, TARGET:: zrg_lwflx_dn    (:,:,:)    !< longwave  3D downward flux           
    REAL(wp), ALLOCATABLE, TARGET:: zrg_swflx_up    (:,:,:)    !< shortwave 3D upward   flux          
    REAL(wp), ALLOCATABLE, TARGET:: zrg_swflx_dn    (:,:,:)    !< shortwave 3D downward flux          
    REAL(wp), ALLOCATABLE, TARGET:: zrg_lwflx_up_clr(:,:,:)    !< longwave  3D upward   flux clear-sky
    REAL(wp), ALLOCATABLE, TARGET:: zrg_lwflx_dn_clr(:,:,:)    !< longwave  3D downward flux clear-sky
    REAL(wp), ALLOCATABLE, TARGET:: zrg_swflx_up_clr(:,:,:)    !< shortwave 3D upward   flux clear-sky
    REAL(wp), ALLOCATABLE, TARGET:: zrg_swflx_dn_clr(:,:,:)    !< shortwave 3D downward flux clear-sky

    ! Pointer to parent patch or local parent patch for reduced grid
    TYPE(t_patch), POINTER       :: ptr_pp

    REAL(wp), DIMENSION(:,:), POINTER :: ptr_reff_qc => NULL(), ptr_reff_qi => NULL()

    TYPE(t_upscale_fields)   :: input_extra_flds, input_extra_2D   !< pointer array for input in upscale routine

    INTEGER   ::   irg_acdnc, irg_fr_land, irg_fr_glac      ! indices of extra fields

    ! Pointer to the acutally used variant of clc:
    REAL(wp), DIMENSION(:,:,:), POINTER ::  ptr_clc => NULL()

    ! Pointers to extra fields
    REAL(wp), DIMENSION(:,:), POINTER ::  ptr_acdnc => NULL()
    REAL(wp), DIMENSION(:),   POINTER ::  ptr_fr_glac => NULL() , ptr_fr_land => NULL()

    ! Variables for debug output
    REAL(wp) :: max_albvisdir, min_albvisdir, max_albvisdif, min_albvisdif, &
                max_albdif, min_albdif, max_tsfc, min_tsfc, max_psfc, min_psfc

    REAL(wp), DIMENSION(:), ALLOCATABLE :: max_pres_ifc, max_pres, max_temp, max_acdnc, &
        max_qv, max_qc, max_qi, max_cc, min_pres_ifc, min_pres, min_temp, min_acdnc, &
        min_qv, min_qc, min_qi, min_cc, max_reff_liq, max_reff_frz, min_reff_liq, min_reff_frz 

    REAL(wp), DIMENSION(pt_patch%nlevp1) :: max_lwflx, min_lwflx, max_swtrans, min_swtrans

    ! Local scalars:
    INTEGER:: jk,jb,jf
    INTEGER:: jg                      !domain id
    INTEGER:: nlev, nlevp1, nlev_rg   !< number of full and half levels
    INTEGER:: nblks_par_c, nblks_lp_c !nblks for reduced grid
    INTEGER:: np, nl                  !< dimension variables for allocation (3d fluxes)
    INTEGER:: rl_start, rl_end
    INTEGER:: i_startblk, i_endblk    !> blocks
    INTEGER:: i_startidx, i_endidx    !< slices
    INTEGER:: i_nchdom                !< domain index
    INTEGER:: i_chidx
    LOGICAL:: l_coupled_reff          ! Use the effective radius from microphysics

    CALL assert_acc_host_only(routine, lacc)

    i_nchdom  = MAX(1,pt_patch%n_childdom)
    jg        = pt_patch%id

    ! number of vertical levels
    nlev   = pt_patch%nlev
    nlevp1 = pt_patch%nlevp1

    ! Flag for using microph. effective radius
    l_coupled_reff = atm_phy_nwp_config(jg)%icpl_rad_reff > 0


    ! Decide which field for cloud cover has to be used:
    IF (atm_phy_nwp_config(jg)%luse_clc_rad) THEN
      ptr_clc => prm_diag%clc_rad
    ELSE
      ptr_clc => prm_diag%clc
    END IF

    !-------------------------------------------------------------------------
    !> Radiation
    !-------------------------------------------------------------------------


      ! section for computing radiation on reduced grid


      IF (msg_level >= 12) &
        &       CALL message(routine, 'RRTM radiation on reduced grid')

      i_chidx     =  pt_patch%parent_child_index

      IF (jg == 1 .AND. .NOT. l_limited_area) THEN
        ptr_pp => pt_par_patch
        nblks_par_c = pt_par_patch%nblks_c
        nblks_lp_c  =  p_patch_local_parent(jg)%nblks_c
      ELSE ! Nested domain with MPI parallelization
        ptr_pp      => p_patch_local_parent(jg)
        nblks_par_c =  ptr_pp%nblks_c
        nblks_lp_c  =  ptr_pp%nblks_c
      ENDIF

      ! Add extra layer for atmosphere above model top if requested
      IF (atm_phy_nwp_config(jg)%latm_above_top) THEN
        IF (jg == 1 .OR. pt_patch%nshift == 0) THEN
          nlev_rg = nlev + 1
        ELSE ! add a specified number levels up to the top of the parent domain in case of vertical nesting
          nlev_rg = MIN(nlev+nexlevs_rrg_vnest, pt_par_patch%nlev)
        ENDIF
      ELSE
        nlev_rg = nlev
      ENDIF

      ! Set dimensions for 3D radiative flux variables
      IF (atm_phy_nwp_config(jg)%l_3d_rad_fluxes) THEN
         np = nproma
         nl = nlev_rg+1
      ELSE
         np = 1
         nl = 1
      END IF


      ALLOCATE (zrg_cosmu0   (nproma,nblks_par_c),     &
        zrg_emis_rad (nproma,nblks_par_c),             &
        zrg_albvisdir(nproma,nblks_par_c),             &
        zrg_albnirdir(nproma,nblks_par_c),             &
        zrg_albvisdif(nproma,nblks_par_c),             &
        zrg_albnirdif(nproma,nblks_par_c),             &
        zrg_albdif   (nproma,nblks_par_c),             &
        zrg_tsfc     (nproma,nblks_par_c),             &
        zrg_rtype    (nproma,nblks_par_c),             &
        zrg_ktype    (nproma,nblks_par_c),             &
        zrg_pres_ifc (nproma,nlev_rg+1,nblks_par_c),   &
        zlp_pres_ifc (nproma,nlev_rg+1,nblks_lp_c ),   &
        zrg_pres     (nproma,nlev_rg  ,nblks_par_c),   &
        zrg_temp     (nproma,nlev_rg  ,nblks_par_c),   &
        zrg_o3       (nproma,nlev_rg  ,nblks_par_c),   &
        zrg_aeq1     (nproma,nlev_rg  ,nblks_par_c),   &
        zrg_aeq2     (nproma,nlev_rg  ,nblks_par_c),   &
        zrg_aeq3     (nproma,nlev_rg  ,nblks_par_c),   &
        zrg_aeq4     (nproma,nlev_rg  ,nblks_par_c),   &
        zrg_aeq5     (nproma,nlev_rg  ,nblks_par_c),   &
        zrg_tot_cld  (nproma,nlev_rg  ,nblks_par_c,3), &
        zlp_tot_cld  (nproma,nlev_rg  ,nblks_lp_c,3),  &
        zrg_clc      (nproma,nlev_rg  ,nblks_par_c),   &
        zrg_aclcov   (nproma,          nblks_par_c),   &
        zrg_lwflx_up_sfc   (nproma,    nblks_par_c),   &
        zrg_trsol_up_toa   (nproma,    nblks_par_c),   &
        zrg_trsol_up_sfc   (nproma,    nblks_par_c),   &
        zrg_trsol_nir_sfc  (nproma,    nblks_par_c),   &
        zrg_trsol_vis_sfc  (nproma,    nblks_par_c),   &
        zrg_trsol_par_sfc  (nproma,    nblks_par_c),   &
        zrg_fr_nir_sfc_diff  (nproma,  nblks_par_c),   &
        zrg_fr_vis_sfc_diff  (nproma,  nblks_par_c),   &
        zrg_fr_par_sfc_diff  (nproma,  nblks_par_c),   &
        zrg_trsol_dn_sfc_diff(nproma,  nblks_par_c),   &
        zrg_trsol_clr_sfc    (nproma,  nblks_par_c),   &
        zrg_lwflx_clr_sfc    (nproma,  nblks_par_c),   &
        zrg_lwflxall    (nproma,nlev_rg+1,nblks_par_c),&
        zrg_trsolall    (nproma,nlev_rg+1,nblks_par_c),&
        zrg_lwflx_up    (np, nl, nblks_par_c),         &
        zrg_lwflx_dn    (np, nl, nblks_par_c),         &   
        zrg_swflx_up    (np, nl, nblks_par_c),         &
        zrg_swflx_dn    (np, nl, nblks_par_c),         &
        zrg_lwflx_up_clr(np, nl, nblks_par_c),         &
        zrg_lwflx_dn_clr(np, nl, nblks_par_c),         &
        zrg_swflx_up_clr(np, nl, nblks_par_c),         &
        zrg_swflx_dn_clr(np, nl, nblks_par_c) )


    ! Set indices for extra fields in the upscaling routine
      irg_acdnc   = 0
      irg_fr_land = 0
      irg_fr_glac = 0
    
      CALL input_extra_flds%construct(nlev_rg)  ! Extra fields in upscaling routine. 3D fields with nlev_rg
      CALL input_extra_2D%construct(1)          ! Extra fields in upscaling routine: 2D fields
    
      IF (l_coupled_reff) THEN
        ALLOCATE(zrg_reff_liq (nproma,nlev_rg,nblks_par_c),   &
             zrg_reff_frz (nproma,nlev_rg,nblks_par_c))
      ELSE
        CALL input_extra_flds%assign(prm_diag%acdnc(:,:,:), irg_acdnc)
        CALL input_extra_2D%assign(ext_data%atm%fr_land , irg_fr_land)
        CALL input_extra_2D%assign(ext_data%atm%fr_glac , irg_fr_glac)
      ENDIF

      ! Allocate output extra arrays
      IF ( input_extra_flds%ntot > 0 )  THEN
        ALLOCATE( zrg_extra_flds(nproma,input_extra_flds%nlev_rg,nblks_par_c,input_extra_flds%ntot) )
      END IF
      IF ( input_extra_2D%ntot > 0 )  THEN
        ALLOCATE( zrg_extra_2D(nproma,nblks_par_c,input_extra_2D%ntot) )
      END IF

      rl_start = 1 ! SR radiation is not set up to handle boundaries of nested domains
      rl_end   = min_rlcell_int

      i_startblk = pt_patch%cells%start_blk(rl_start,1)
      i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

      IF (p_test_run) THEN
        CALL get_indices_c(pt_patch, i_startblk, i_startblk, i_endblk, &
             &                         i_startidx, i_endidx, rl_start, rl_end)
        zrg_lwflxall(:,:,:) = 0._wp
        zrg_trsolall(:,:,:) = 0._wp
        zrg_lwflx_up_sfc(1:i_startidx-1,i_startblk) = 0
        zrg_trsol_up_toa(1:i_startidx-1,i_startblk) = 0
        zrg_trsol_up_sfc(1:i_startidx-1,i_startblk) = 0
        zrg_trsol_nir_sfc(1:i_startidx-1,i_startblk) = 0
        zrg_trsol_vis_sfc(1:i_startidx-1,i_startblk) = 0
        zrg_trsol_par_sfc(1:i_startidx-1,i_startblk) = 0
        zrg_fr_nir_sfc_diff(1:i_startidx-1,i_startblk) = 0
        zrg_fr_vis_sfc_diff(1:i_startidx-1,i_startblk) = 0
        zrg_fr_par_sfc_diff(1:i_startidx-1,i_startblk) = 0
        zrg_trsol_dn_sfc_diff(1:i_startidx-1,i_startblk) = 0
        zrg_trsol_clr_sfc(1:i_startidx-1,i_startblk) = 0
        zrg_lwflx_clr_sfc(1:i_startidx-1,i_startblk) = 0
        IF (atm_phy_nwp_config(jg)%l_3d_rad_fluxes) THEN
          zrg_lwflx_up    (:,:,:) = 0._wp
          zrg_lwflx_dn    (:,:,:) = 0._wp
          zrg_swflx_up    (:,:,:) = 0._wp
          zrg_swflx_dn    (:,:,:) = 0._wp
          zrg_lwflx_up_clr(:,:,:) = 0._wp
          zrg_lwflx_dn_clr(:,:,:) = 0._wp
          zrg_swflx_up_clr(:,:,:) = 0._wp
          zrg_swflx_dn_clr(:,:,:) = 0._wp
        ENDIF
     ENDIF

      DO jb = i_startblk, i_endblk

        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          &                       i_startidx, i_endidx, rl_start, rl_end)

        ! Loop starts with 1 instead of i_startidx because the start index is missing in RRTM
        prm_diag%tsfctrad(1:i_endidx,jb) = lnd_prog%t_g(1:i_endidx,jb)

      ENDDO ! blocks

      CALL upscale_rad_input(pt_patch%id, pt_par_patch%id,              &
        & nlev_rg,                                                      &
        & prm_diag%lw_emiss, prm_diag%cosmu0,                           &
        & prm_diag%albvisdir, prm_diag%albnirdir, prm_diag%albvisdif,   &
        & prm_diag%albnirdif, prm_diag%albdif, prm_diag%tsfctrad,       &
        & prm_diag%ktype, pt_diag%pres_ifc, pt_diag%pres,               &
        & pt_diag%temp,                                                 &
        & prm_diag%tot_cld, ptr_clc,                                    &
        & ext_data%atm%o3, zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,           &
        & zrg_emis_rad,                                                 &
        & zrg_cosmu0, zrg_albvisdir, zrg_albnirdir, zrg_albvisdif,      &
        & zrg_albnirdif, zrg_albdif, zrg_tsfc, zrg_rtype, zrg_pres_ifc, &
        & zrg_pres, zrg_temp,                                           &
        & zrg_tot_cld, zrg_clc, zrg_o3,                                 &
        & zrg_aeq1, zrg_aeq2, zrg_aeq3, zrg_aeq4, zrg_aeq5,             &
        & zlp_pres_ifc, zlp_tot_cld, prm_diag%buffer_rrg,               &
        & atm_phy_nwp_config(jg)%icpl_rad_reff,                         &
        & prm_diag%reff_qc, prm_diag%reff_qi,                           &
        & zrg_reff_liq, zrg_reff_frz,input_extra_flds, zrg_extra_flds,  &
        & input_extra_2D, zrg_extra_2D)
    
      IF (jg == 1 .AND. l_limited_area) THEN
        rl_start = grf_fbk_start_c
      ELSE
        rl_start = grf_ovlparea_start_c
      ENDIF
      rl_end   = min_rlcell_int

      i_startblk = ptr_pp%cells%start_blk(rl_start,i_chidx)
      i_endblk   = ptr_pp%cells%end_blk(rl_end,i_chidx)

      ! Debug output of radiation input fields
      IF (msg_level >= 16) THEN

        ALLOCATE(max_pres_ifc(nlev_rg), max_pres(nlev_rg), max_temp(nlev_rg), max_acdnc(nlev_rg), &
                 max_qv(nlev_rg), max_qc(nlev_rg), max_qi(nlev_rg), max_cc(nlev_rg),              &
                 min_pres_ifc(nlev_rg), min_pres(nlev_rg), min_temp(nlev_rg), min_acdnc(nlev_rg), &
                 min_qv(nlev_rg), min_qc(nlev_rg), min_qi(nlev_rg), min_cc(nlev_rg)  )
                
        max_albvisdir = 0._wp
        min_albvisdir = 1.e10_wp
        max_albvisdif = 0._wp
        min_albvisdif = 1.e10_wp
        max_albdif = 0._wp
        min_albdif = 1.e10_wp
        max_tsfc = 0._wp
        min_tsfc = 1.e10_wp
        max_psfc = 0._wp
        min_psfc = 1.e10_wp
        max_pres_ifc = 0._wp
        max_pres    = 0._wp
        max_temp     = 0._wp
        max_acdnc    = 0._wp
        max_qv = 0._wp
        max_qc = 0._wp
        max_qi = 0._wp
        max_cc  = 0._wp
        min_pres_ifc = 1.e10_wp
        min_pres    = 1.e10_wp
        min_temp     = 1.e10_wp
        min_acdnc    = 1.e10_wp
        min_qv = 1.e10_wp
        min_qc = 1.e10_wp
        min_qi = 1.e10_wp
        min_cc  = 1.e10_wp
        IF (l_coupled_reff) THEN
          ALLOCATE(  min_reff_liq(nlev_rg), max_reff_liq(nlev_rg), &
                     min_reff_frz(nlev_rg), max_reff_frz(nlev_rg) )
          max_reff_liq = 0._wp
          max_reff_frz = 0._wp
          min_reff_liq = 1.e10_wp
          min_reff_frz = 1.e10_wp
        END IF


        DO jb = i_startblk, i_endblk

         CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)

         max_albvisdir = MAX(max_albvisdir,MAXVAL(zrg_albvisdir(i_startidx:i_endidx,jb)))
         min_albvisdir = MIN(min_albvisdir,MINVAL(zrg_albvisdir(i_startidx:i_endidx,jb)))
         max_albvisdif = MAX(max_albvisdif,MAXVAL(zrg_albvisdif(i_startidx:i_endidx,jb)))
         min_albvisdif = MIN(min_albvisdif,MINVAL(zrg_albvisdif(i_startidx:i_endidx,jb)))
         max_albdif    = MAX(max_albdif,   MAXVAL(zrg_albdif(i_startidx:i_endidx,jb)))
         min_albdif    = MIN(min_albdif,   MINVAL(zrg_albdif(i_startidx:i_endidx,jb)))
         max_tsfc = MAX(max_tsfc,MAXVAL(zrg_tsfc(i_startidx:i_endidx,jb)))
         min_tsfc = MIN(min_tsfc,MINVAL(zrg_tsfc(i_startidx:i_endidx,jb)))
         max_psfc = MAX(max_psfc,MAXVAL(zrg_pres_ifc(i_startidx:i_endidx,nlev_rg+1,jb)))
         min_psfc = MIN(min_psfc,MINVAL(zrg_pres_ifc(i_startidx:i_endidx,nlev_rg+1,jb)))
         DO jk = 1, nlev_rg
          max_pres_ifc(jk) = MAX(max_pres_ifc(jk),MAXVAL(zrg_pres_ifc(i_startidx:i_endidx,jk,jb)))
          max_pres(jk)    = MAX(max_pres(jk),MAXVAL(zrg_pres     (i_startidx:i_endidx,jk,jb)))
          max_temp(jk)     = MAX(max_temp(jk),MAXVAL(zrg_temp     (i_startidx:i_endidx,jk,jb)))
          max_qv(jk) = MAX(max_qv(jk),MAXVAL(zrg_tot_cld(i_startidx:i_endidx,jk,jb,iqv)))
          max_qc(jk) = MAX(max_qc(jk),MAXVAL(zrg_tot_cld(i_startidx:i_endidx,jk,jb,iqc)))
          max_qi(jk) = MAX(max_qi(jk),MAXVAL(zrg_tot_cld(i_startidx:i_endidx,jk,jb,iqi)))
          max_cc(jk)  = MAX(max_cc(jk),MAXVAL(zrg_clc(i_startidx:i_endidx,jk,jb)))
          min_pres_ifc(jk) = MIN(min_pres_ifc(jk),MINVAL(zrg_pres_ifc(i_startidx:i_endidx,jk,jb)))
          min_pres(jk)    = MIN(min_pres(jk),MINVAL(zrg_pres     (i_startidx:i_endidx,jk,jb)))
          min_temp(jk)     = MIN(min_temp(jk),MINVAL(zrg_temp     (i_startidx:i_endidx,jk,jb)))
          min_qv(jk) = MIN(min_qv(jk),MINVAL(zrg_tot_cld(i_startidx:i_endidx,jk,jb,iqv)))
          min_qc(jk) = MIN(min_qc(jk),MINVAL(zrg_tot_cld(i_startidx:i_endidx,jk,jb,iqc)))
          min_qi(jk) = MIN(min_qi(jk),MINVAL(zrg_tot_cld(i_startidx:i_endidx,jk,jb,iqi)))
          min_cc(jk)  = MIN(min_cc(jk),MINVAL(zrg_clc(i_startidx:i_endidx,jk,jb)))
         ENDDO
         IF (irg_acdnc > 0) THEN
           DO jk = 1, nlev_rg
             max_acdnc(jk)    = MAX(max_acdnc(jk),MAXVAL(zrg_extra_flds (i_startidx:i_endidx,jk,jb,irg_acdnc)))
             min_acdnc(jk)    = MIN(min_acdnc(jk),MINVAL(zrg_extra_flds (i_startidx:i_endidx,jk,jb,irg_acdnc)))
           ENDDO
         END IF
         IF (l_coupled_reff) THEN
           DO jk = 1, nlev_rg
             max_reff_liq(jk) = MAX(max_reff_liq(jk),MAXVAL(zrg_reff_liq    (i_startidx:i_endidx,jk,jb)))
             max_reff_frz(jk) = MAX(max_reff_frz(jk),MAXVAL(zrg_reff_frz    (i_startidx:i_endidx,jk,jb)))
             min_reff_liq(jk) = MIN(min_reff_liq(jk),MINVAL(zrg_reff_liq    (i_startidx:i_endidx,jk,jb)))
             min_reff_frz(jk) = MIN(min_reff_frz(jk),MINVAL(zrg_reff_frz    (i_startidx:i_endidx,jk,jb)))
           END DO
         END IF

        ENDDO ! blocks

        max_albvisdir = global_max(max_albvisdir)
        min_albvisdir = global_min(min_albvisdir)
        max_albvisdif = global_max(max_albvisdif)
        min_albvisdif = global_min(min_albvisdif)
        max_tsfc = global_max(max_tsfc)
        min_tsfc = global_min(min_tsfc)
        max_psfc = global_max(max_psfc)
        min_psfc = global_min(min_psfc)
        max_pres_ifc = global_max(max_pres_ifc)
        max_pres    = global_max(max_pres)
        max_temp     = global_max(max_temp)
        max_qv = global_max(max_qv)
        max_qc = global_max(max_qc)
        max_qi = global_max(max_qi)
        max_cc  = global_max(max_cc)
        min_pres_ifc = global_min(min_pres_ifc)
        min_pres    = global_min(min_pres)
        min_temp     = global_min(min_temp)
        min_qv = global_min(min_qv)
        min_qc = global_min(min_qc)
        min_qi = global_min(min_qi)
        min_cc  = global_min(min_cc)
        IF (irg_acdnc > 0) THEN
          max_acdnc    = global_max(max_acdnc)
          min_acdnc    = global_min(min_acdnc)
        END IF
        IF (l_coupled_reff) THEN
          max_reff_liq = global_max(max_reff_liq)
          max_reff_frz = global_max(max_reff_frz)        
          min_reff_liq = global_min(min_reff_liq)
          min_reff_frz = global_min(min_reff_frz)
        END IF

        WRITE(message_text,'(a,4f12.8)') 'max/min alb = ', max_albvisdir, min_albvisdir, &
          max_albvisdif, min_albvisdif
        CALL message(routine, message_text)

        WRITE(message_text,'(a,2f10.3,2f10.2)') 'max/min sfc temp/pres = ', max_tsfc, min_tsfc, &
          max_psfc, min_psfc
        CALL message(routine, message_text)

        WRITE(message_text,'(a)') 'max/min pres_ifc, pres, temp, acdnc'
        CALL message(routine, message_text)

        DO jk = 1, nlev_rg
          WRITE(message_text,'(i4,4f10.2,2f10.3,2f12.1)') jk,max_pres_ifc(jk), min_pres_ifc(jk), &
            max_pres(jk), min_pres(jk), max_temp(jk), min_temp(jk), max_acdnc(jk), min_acdnc(jk)
          CALL message(routine, message_text)
        ENDDO
        IF (l_coupled_reff) THEN
          WRITE(message_text,'(a)') 'max/min reff_liq, reff_frz'
          CALL message(routine, TRIM(message_text))

          DO jk = 1, nlev_rg
            WRITE(message_text,'(i4,4e13.5)') jk,max_reff_liq(jk), min_reff_liq(jk),                &
                 max_reff_frz(jk), min_reff_frz(jk)
            CALL message(routine, TRIM(message_text))
          ENDDO
        END IF
      
        WRITE(message_text,'(a)') 'max/min QV, QC, QI, CC'
        CALL message(routine, message_text)

        DO jk = 1, nlev_rg
          WRITE(message_text,'(i4,8e13.5)') jk,max_qv(jk), min_qv(jk), max_qc(jk), min_qc(jk), &
             max_qi(jk), min_qi(jk), max_cc(jk), min_cc(jk)
          CALL message(routine, message_text)
        ENDDO

        DEALLOCATE(max_pres_ifc, max_pres, max_temp, max_acdnc, max_qv, max_qc, max_qi, max_cc, &
                   min_pres_ifc, min_pres, min_temp, min_acdnc, min_qv, min_qc, min_qi, min_cc )
        IF (l_coupled_reff) DEALLOCATE(max_reff_liq, min_reff_liq, max_reff_frz, min_reff_frz)

      ENDIF ! msg_level >= 16


!FIXME: PGI + OpenMP produce deadlock in this loop. Compiler bug suspected
!ICON_OMP PARALLEL DO PRIVATE(jb,jk,jf,i_startidx,i_endidx) ICON_OMP_GUIDED_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, &
          &                         i_startidx, i_endidx, rl_start, rl_end)


        ! Unfortunately, the coding of SR radiation is not compatible with the presence
        ! of nested domains. Therefore, the normally unused elements of the first block
        ! need to be filled with dummy values
        IF ( (jg > 1 .OR. l_limited_area) .AND. jb == i_startblk .AND. i_endidx >= i_startidx) THEN
          zrg_emis_rad  (1:i_startidx-1,jb) = zrg_emis_rad  (i_startidx,jb)
          zrg_cosmu0    (1:i_startidx-1,jb) = zrg_cosmu0    (i_startidx,jb)
          zrg_albvisdir (1:i_startidx-1,jb) = zrg_albvisdir (i_startidx,jb)
          zrg_albnirdir (1:i_startidx-1,jb) = zrg_albnirdir (i_startidx,jb)
          zrg_albvisdif (1:i_startidx-1,jb) = zrg_albvisdif (i_startidx,jb)
          zrg_albnirdif (1:i_startidx-1,jb) = zrg_albnirdif (i_startidx,jb)
          zrg_albdif    (1:i_startidx-1,jb) = zrg_albdif    (i_startidx,jb)
          zrg_tsfc      (1:i_startidx-1,jb) = zrg_tsfc      (i_startidx,jb)
          zrg_rtype     (1:i_startidx-1,jb) = zrg_rtype     (i_startidx,jb)
          zrg_pres_ifc (1:i_startidx-1,nlev_rg+1,jb) = zrg_pres_ifc (i_startidx,nlev_rg+1,jb)
          DO jf=1,input_extra_2D%ntot
            zrg_extra_2D   (1:i_startidx-1,jb,jf) = zrg_extra_2D   (i_startidx,jb,jf)
          END DO
          DO jk = 1, nlev_rg
            zrg_pres_ifc (1:i_startidx-1,jk,jb) = zrg_pres_ifc (i_startidx,jk,jb)
            zrg_pres     (1:i_startidx-1,jk,jb) = zrg_pres     (i_startidx,jk,jb)
            zrg_temp     (1:i_startidx-1,jk,jb) = zrg_temp     (i_startidx,jk,jb)
            zrg_o3       (1:i_startidx-1,jk,jb) = zrg_o3       (i_startidx,jk,jb)
            zrg_aeq1     (1:i_startidx-1,jk,jb) = zrg_aeq1     (i_startidx,jk,jb)
            zrg_aeq2     (1:i_startidx-1,jk,jb) = zrg_aeq2     (i_startidx,jk,jb)
            zrg_aeq3     (1:i_startidx-1,jk,jb) = zrg_aeq3     (i_startidx,jk,jb)
            zrg_aeq4     (1:i_startidx-1,jk,jb) = zrg_aeq4     (i_startidx,jk,jb)
            zrg_aeq5     (1:i_startidx-1,jk,jb) = zrg_aeq5     (i_startidx,jk,jb)
            zrg_tot_cld  (1:i_startidx-1,jk,jb,iqv) = zrg_tot_cld(i_startidx,jk,jb,iqv)
            zrg_tot_cld  (1:i_startidx-1,jk,jb,iqc) = zrg_tot_cld(i_startidx,jk,jb,iqc)
            zrg_tot_cld  (1:i_startidx-1,jk,jb,iqi) = zrg_tot_cld(i_startidx,jk,jb,iqi)
            zrg_clc      (1:i_startidx-1,jk,jb) = zrg_clc(i_startidx,jk,jb)
          ENDDO
          DO jf=1,input_extra_flds%ntot
            DO jk = 1,nlev_rg
              zrg_extra_flds(1:i_startidx-1,jk,jb,jf) = zrg_extra_flds(i_startidx,jk,jb,jf)
            END DO
          END DO
          IF (l_coupled_reff) THEN
            DO jk = 1, nlev_rg
              zrg_reff_liq (1:i_startidx-1,jk,jb) = zrg_reff_liq (i_startidx,jk,jb)
              zrg_reff_frz (1:i_startidx-1,jk,jb) = zrg_reff_frz (i_startidx,jk,jb)
            END DO
          END IF
        ENDIF

      END DO



!FIXME: PGI + OpenMP produce deadlock in this loop. Compiler bug suspected
!ICON_OMP PARALLEL DO PRIVATE(jb,jk,i_startidx,i_endidx,dust_tunefac,              &    
!ICON_OMP              ptr_acdnc, ptr_fr_land,ptr_fr_glac,ptr_reff_qc,ptr_reff_qi) &
!ICON_OMP ICON_OMP_GUIDED_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, &
          &                         i_startidx, i_endidx, rl_start, rl_end)

        IF (tune_dust_abs > 0._wp) THEN
!DIR$ NOINLINE
          CALL tune_dust(ptr_pp%cells%center(:,jb)%lat,ptr_pp%cells%center(:,jb)%lon,i_endidx,dust_tunefac)
        ELSE
          dust_tunefac(:,:) = 1._wp
        ENDIF

        ! Type of convection is required as INTEGER field
        zrg_ktype(1:i_endidx,jb) = NINT(zrg_rtype(1:i_endidx,jb))

      NULLIFY(ptr_acdnc,ptr_fr_land,ptr_fr_glac,ptr_reff_qc,ptr_reff_qi)

      IF (l_coupled_reff) THEN
        ptr_reff_qc => zrg_reff_liq(:,:,jb)
        ptr_reff_qi => zrg_reff_frz(:,:,jb)
      ENDIF

      IF ( irg_acdnc   > 0 ) ptr_acdnc   => zrg_extra_flds(:,:,jb,irg_acdnc)
      IF ( irg_fr_land > 0 ) ptr_fr_land => zrg_extra_2D(:,jb,irg_fr_land)
      IF ( irg_fr_glac > 0 ) ptr_fr_glac => zrg_extra_2D(:,jb,irg_fr_glac)

        CALL radiation_nwp(               &
                                !
                                ! input
                                ! -----
                                !
          & current_date                     ,&!< in current date
                                ! indices and dimensions
          & jg          =jg                  ,&!< in domain index
          & jb          =jb                  ,&!< in block index
          & icpl_reff  =atm_phy_nwp_config(jg)%icpl_rad_reff,&  !< in option for radiation reff coupling
          & jcs         =i_startidx          ,&!< in  start index for loop over block
          & jce         =i_endidx            ,&!< in  end   index for loop over block
          & kbdim       =nproma              ,&!< in  dimension of block over cells
          & klev        =nlev_rg             ,&!< in  number of full levels = number of layers
          & klevp1      =nlev_rg+1           ,&!< in  number of half levels = number of layer ifcs
                                !
          & ktype       =zrg_ktype(:,jb)     ,&!< in type of convection
                                !
          & zland       =ptr_fr_land         ,&!< in land mask,     1. over land
          & zglac       =ptr_fr_glac         ,&!< in glacier mask,  1. over land ice
                                !
          & cos_mu0     =zrg_cosmu0  (:,jb)  ,&!< in    cos of zenith angle mu0
          & alb_vis_dir=zrg_albvisdir(:,jb)  ,&!< in    surface albedo for visible range, direct
          & alb_nir_dir=zrg_albnirdir(:,jb)  ,&!< in    surface albedo for near IR range, direct
          & alb_vis_dif=zrg_albvisdif(:,jb)  ,&!< in    surface albedo for visible range, diffuse
          & alb_nir_dif=zrg_albnirdif(:,jb)  ,&!< in    surface albedo for near IR range, diffuse
          & emis_rad   =zrg_emis_rad(:,jb)   ,&!< in longwave surface emissivity
          & tk_sfc     =zrg_tsfc     (:,jb)       ,&!< in    surface temperature
                                !
                                ! atmosphere: pressure, tracer mixing ratios and temperature
          & pp_hl      =zrg_pres_ifc(:,:,jb)    ,&!< in    pressure at half levels at t-dt [Pa]
          & pp_fl      =zrg_pres    (:,:,jb)    ,&!< in    pressure at full levels at t-dt [Pa]
          & tk_fl      =zrg_temp    (:,:,jb)    ,&!< in    temperature at full level at t-dt
          & qm_vap     =zrg_tot_cld (:,:,jb,iqv),&!< in    water vapor mass mixing ratio at t-dt
          & qm_liq     =zrg_tot_cld (:,:,jb,iqc),&!< in    cloud water mass mixing ratio at t-dt
          & qm_ice     =zrg_tot_cld (:,:,jb,iqi),&!< in    cloud ice mass mixing ratio at t-dt
          & qm_o3      = zrg_o3     (:,:,jb)    ,&!< in    O3
          & cdnc       =ptr_acdnc               ,&!< in    cloud droplet numb. conc. [1/m**3]
          & reff_liq   =ptr_reff_qc             ,&!< in    effective radius liquid phase. [m]
          & reff_frz   =ptr_reff_qi             ,&!< in    effective radius frozen phase. [m]
          & cld_frc    =zrg_clc    (:,:,jb)     ,&!< in    cld_frac = cloud fraction [m2/m2]
          & zaeq1      = zrg_aeq1(:,:,jb)       ,&!< in aerosol continental
          & zaeq2      = zrg_aeq2(:,:,jb)       ,&!< in aerosol maritime
          & zaeq3      = zrg_aeq3(:,:,jb)       ,&!< in aerosol urban
          & zaeq4      = zrg_aeq4(:,:,jb)       ,&!< in aerosol volcano ashes
          & zaeq5      = zrg_aeq5(:,:,jb)       ,&!< in aerosol stratospheric background
          & dust_tunefac = dust_tunefac (:,:)   ,&!< in LW tuning factor for dust aerosol
                                !
                                ! output
                                ! ------
                                !
          & cld_cvr    =  zrg_aclcov    (:,jb), &      !< out   cloud cover in a column [m2/m2]
          & flx_lw_net =  zrg_lwflxall(:,:,jb), &      !< out terrestrial flux, all sky, net down
          & trsol_net  =  zrg_trsolall(:,:,jb), &      !< out solar transmissivity, all sky, net down
          & flx_uplw_sfc = zrg_lwflx_up_sfc(:,jb), &   !< out longwave upward flux at surface
          & trsol_up_toa = zrg_trsol_up_toa(:,jb), &   !< out upward solar transmissivity at TOA
          & trsol_up_sfc = zrg_trsol_up_sfc(:,jb), &   !< out upward solar transmissivity at surface
          & trsol_nir_sfc = zrg_trsol_nir_sfc(:,jb), & !< downward transmissivity for near-infrared rad. at surface
          & trsol_vis_sfc = zrg_trsol_vis_sfc(:,jb), & !< downward transmissivity for visible rad. at surface
          & trsol_par_sfc = zrg_trsol_par_sfc(:,jb), & !< downward transmissivity for photosynthetically active rad. at surface
          & fr_nir_sfc_diffus = zrg_fr_nir_sfc_diff(:,jb), & !< diffuse fraction of downward near-infrared rad. at surface
          & fr_vis_sfc_diffus = zrg_fr_vis_sfc_diff(:,jb), & !< diffuse fraction of downward visible rad. at surface
          & fr_par_sfc_diffus = zrg_fr_par_sfc_diff(:,jb), & !< diffuse fraction of downward photosynthetically active rad. at surface
          & trsol_dn_sfc_diffus = zrg_trsol_dn_sfc_diff(:,jb), &  !< out downward diffuse solar transmissivity at surface
          & trsol_clr_sfc = zrg_trsol_clr_sfc(:,jb), & !< out clear-sky net transmissvity at surface (used with reduced grid only)
          & lwflx_clr_sfc = zrg_lwflx_clr_sfc(:,jb), & !< out clear-sky net LW flux at surface
          !optional output: 3D flux output
          &  flx_lw_dn     = zrg_lwflx_dn(:,:,jb),     & !< Downward LW flux (all-sky)   [Wm2]
          &  flx_sw_dn     = zrg_swflx_dn(:,:,jb),     & !< Downward SW flux (all-sky)   [Wm2]
          &  flx_lw_up     = zrg_lwflx_up(:,:,jb),     & !< Upward LW flux   (all-sky)   [Wm2]
          &  flx_sw_up     = zrg_swflx_up(:,:,jb),     & !< Upward SW flux   (all-sky)   [Wm2]
          &  flx_lw_dn_clr = zrg_lwflx_dn_clr(:,:,jb), & !< Downward LW flux (clear sky) [Wm2]
          &  flx_sw_dn_clr = zrg_swflx_dn_clr(:,:,jb), & !< Downward SW flux (clear sky) [Wm2]
          &  flx_lw_up_clr = zrg_lwflx_up_clr(:,:,jb), & !< Upward LW flux   (clear sky) [Wm2]
          &  flx_sw_up_clr = zrg_swflx_up_clr(:,:,jb) )  !< Upward SW flux   (clear sky) [Wm2]

      ENDDO ! blocks

      CALL downscale_rad_output(pt_patch%id, pt_par_patch%id, nlev_rg, zrg_aclcov,                     &
        &  zrg_lwflxall, zrg_trsolall, zrg_trsol_clr_sfc, zrg_lwflx_clr_sfc, zrg_lwflx_up_sfc,         &
        &  zrg_trsol_up_toa, zrg_trsol_up_sfc, zrg_trsol_nir_sfc, zrg_trsol_vis_sfc, zrg_trsol_par_sfc,&
        &  zrg_fr_nir_sfc_diff, zrg_fr_vis_sfc_diff, zrg_fr_par_sfc_diff, zrg_trsol_dn_sfc_diff,       &
        &  zrg_tsfc, zrg_albdif, zrg_emis_rad, zrg_cosmu0, zrg_tot_cld, zlp_tot_cld, zrg_pres_ifc,     &
        &  zlp_pres_ifc, prm_diag%tsfctrad, prm_diag%albdif, aclcov, prm_diag%lwflxall,                &
        &  prm_diag%trsolall, prm_diag%lwflx_up_sfc_rs, prm_diag%trsol_up_toa,                         &
        &  prm_diag%trsol_up_sfc, prm_diag%trsol_nir_sfc, prm_diag%trsol_vis_sfc, prm_diag%trsol_par_sfc,&
        &  prm_diag%fr_nir_sfc_diff, prm_diag%fr_vis_sfc_diff, prm_diag%fr_par_sfc_diff,               &
        &  prm_diag%trsol_dn_sfc_diff, prm_diag%trsolclr_sfc, prm_diag%lwflxclr_sfc,                   &
        &  zrg_lwflx_up         , zrg_lwflx_dn         , zrg_swflx_up         , zrg_swflx_dn,          &
        &  zrg_lwflx_up_clr     , zrg_lwflx_dn_clr     , zrg_swflx_up_clr     , zrg_swflx_dn_clr,      &
        &  prm_diag%lwflx_up    , prm_diag%lwflx_dn    , prm_diag%swflx_up    , prm_diag%swflx_dn,     &
        &  prm_diag%lwflx_up_clr, prm_diag%lwflx_dn_clr, prm_diag%swflx_up_clr, prm_diag%swflx_dn_clr  )

      ! Debug output of radiation output fields
      IF (msg_level >= 16) THEN
        max_lwflx = 0._wp
        min_lwflx = 1.e10_wp
        max_swtrans = 0._wp
        min_swtrans = 1.e10_wp

        rl_start = grf_bdywidth_c + 1
        rl_end   = min_rlcell_int

        i_startblk = pt_patch%cells%start_blk(rl_start,1)
        i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

        DO jb = i_startblk, i_endblk
         CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)

         DO jk = 1, nlevp1
          max_lwflx(jk)   = MAX(max_lwflx(jk),  MAXVAL(prm_diag%lwflxall(i_startidx:i_endidx,jk,jb)))
          max_swtrans(jk) = MAX(max_swtrans(jk),MAXVAL(prm_diag%trsolall(i_startidx:i_endidx,jk,jb)))
          min_lwflx(jk)   = MIN(min_lwflx(jk),  MINVAL(prm_diag%lwflxall(i_startidx:i_endidx,jk,jb)))
          min_swtrans(jk) = MIN(min_swtrans(jk),MINVAL(prm_diag%trsolall(i_startidx:i_endidx,jk,jb)))
         ENDDO
        ENDDO ! blocks

        max_lwflx = global_max(max_lwflx)
        min_lwflx = global_min(min_lwflx)
        max_swtrans = global_max(max_swtrans)
        min_swtrans = global_min(min_swtrans)


        WRITE(message_text,'(a)') 'max/min LW flux, SW transmissivity'
        CALL message(routine, message_text)

        DO jk = 1, nlevp1
          WRITE(message_text,'(i4,2f10.3,2f10.7)') jk,max_lwflx(jk), min_lwflx(jk), &
            max_swtrans(jk), min_swtrans(jk)
          CALL message(routine, message_text)
        ENDDO

      ENDIF ! msg_level >= 16

      DEALLOCATE (zrg_cosmu0, zrg_albvisdir, zrg_albnirdir, zrg_albvisdif, zrg_albnirdif, &
        zrg_albdif, zrg_tsfc, zrg_pres_ifc, zrg_pres, zrg_temp, zrg_o3, zrg_ktype,        &
        zrg_aeq1,zrg_aeq2,zrg_aeq3,zrg_aeq4,zrg_aeq5, zrg_tot_cld, zrg_clc,               &
        zrg_aclcov, zrg_lwflxall, zrg_trsolall, zrg_lwflx_up_sfc, zrg_trsol_up_toa,       &
        zrg_trsol_up_sfc, zrg_trsol_nir_sfc, zrg_trsol_vis_sfc, zrg_trsol_par_sfc,        &
        zrg_fr_nir_sfc_diff, zrg_fr_vis_sfc_diff, zrg_fr_par_sfc_diff,                    &
        zrg_trsol_dn_sfc_diff, zrg_trsol_clr_sfc,                                         &
        zrg_lwflx_clr_sfc, zrg_emis_rad, zlp_pres_ifc, zlp_tot_cld,                       &
        zrg_lwflx_up    , zrg_lwflx_dn    , zrg_swflx_up    , zrg_swflx_dn,               &
        zrg_lwflx_up_clr, zrg_lwflx_dn_clr, zrg_swflx_up_clr, zrg_swflx_dn_clr            )
      IF (l_coupled_reff) DEALLOCATE(zrg_reff_liq,zrg_reff_frz)     
      IF (input_extra_flds%ntot > 0 ) DEALLOCATE(zrg_extra_flds)
      IF (input_extra_2D%ntot > 0   ) DEALLOCATE(zrg_extra_2D  )
    
      CALL input_extra_flds%destruct()
      CALL input_extra_2D%destruct()

  END SUBROUTINE nwp_rrtm_radiation_reduced
  !---------------------------------------------------------------------------------------

END MODULE mo_nwp_rrtm_interface


