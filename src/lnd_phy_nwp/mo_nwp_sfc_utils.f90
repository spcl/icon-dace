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

! Utility routines related to the TERRA surface model

!----------------------------
#include "omp_definitions.inc"
!----------------------------
#if defined __xlC__
@PROCESS SPILL(564)
#endif
MODULE mo_nwp_sfc_utils

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_model_domain,        ONLY: t_patch
  USE mo_physical_constants,  ONLY: tmelt, tf_salt, grav, salinity_fac, rhoh2o
  USE mo_math_constants,      ONLY: dbl_eps, rad2deg
  USE mo_impl_constants,      ONLY: min_rlcell_int, min_rlcell, LSS_JSBACH, &
    &                               MODE_IAU, ALB_SI_MISSVAL, MAX_CHAR_LENGTH
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE sfc_flake_data,         ONLY: tpl_T_r, C_T_min, rflk_depth_bs_ref, h_Ice_min_flk
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_ext_data_init,       ONLY: diagnose_ext_aggr, interpol_monthly_mean, vege_clim
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag, t_lnd_state
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_parallel_config,     ONLY: nproma
  USE mo_grid_config,         ONLY: l_limited_area, n_dom
  USe mo_extpar_config,       ONLY: itopo, itype_vegetation_cycle
  USE mo_lnd_nwp_config,      ONLY: nlev_soil, nlev_snow, ntiles_total, ntiles_water, &
    &                               lseaice, llake, lmulti_snow, idiag_snowfrac, ntiles_lnd, &
    &                               lsnowtile, isub_water, isub_seaice, isub_lake,    &
    &                               itype_interception, lterra_urb, l2lay_rho_snow, lprog_albsi, itype_trvg, &
                                    itype_snowevap, zml_soil, dzsoil, frsi_min, hice_min
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_nwp_tuning_config,   ONLY: tune_minsnowfrac
  USE mo_initicon_config,     ONLY: init_mode_soil, ltile_coldstart, init_mode, lanaread_tseasfc, use_lakeiceana, &
                                    icpl_da_snowalb
  USE mo_run_config,          ONLY: msg_level
  USE sfc_terra_init,         ONLY: terra_init
  USE sfc_flake,              ONLY: flake_init
  USE sfc_seaice,             ONLY: seaice_init_nwp, seaice_coldinit_albsi_nwp
  USE sfc_terra_data,         ONLY: cadp, cf_snow, crhosmin_ml, crhosmax_ml
  USE turb_data,              ONLY: c_lnd, c_sea, c_stm
  USE mo_satad,               ONLY: sat_pres_water, sat_pres_ice, spec_humi
  USE mo_sync,                ONLY: global_max, global_min
  USE mo_nonhydro_types,      ONLY: t_nh_diag
  USE mo_dynamics_config,     ONLY: nnow_rcf, nnew_rcf
  USE mtime,                  ONLY: datetime, MAX_DATETIME_STR_LEN, datetimeToString
  USE mo_idx_list,            ONLY: t_idx_list_blocked, compare_sets
  USE mo_util_table,          ONLY: t_table, initialize_table, add_table_column, &
    &                                set_table_entry, print_table, finalize_table
  USE mo_util_string,         ONLY: int2string
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_index_list,          ONLY: generate_index_list
  USE mo_fortran_tools,       ONLY: set_acc_host_or_device, assert_acc_device_only
  USE mo_timer,               ONLY: ltimer, timer_nh_diagnostics, timer_start, timer_stop

  USE mo_lnd_nwp_config,      ONLY: lcuda_graph_lnd

#ifdef __NVCOMPILER
  USE mo_coupling_config,     ONLY: is_coupled_to_ocean
#endif

  IMPLICIT NONE

  PRIVATE

#ifdef __SX__
! parameter for loop unrolling
INTEGER, PARAMETER :: nlsoil= 8
#endif

  REAL(KIND=wp), PARAMETER ::            &
   & csmall_hice = 0.5_wp*h_Ice_min_flk    !< small value to handle lake-ice fraction 
                                           !< (1/2 of minimum lake-ice thickness)

  PUBLIC :: nwp_surface_init
  PUBLIC :: diag_snowfrac_tg
  PUBLIC :: aggregate_landvars, aggr_landvars
  PUBLIC :: aggregate_tg_qvs
  PUBLIC :: update_idx_lists_lnd
  PUBLIC :: update_idx_lists_sea
  PUBLIC :: process_sst_and_seaice
  PUBLIC :: sst_add_climatological_incr
  PUBLIC :: update_ndvi_dependent_fields
  PUBLIC :: init_snowtile_lists
  PUBLIC :: init_sea_lists
  PUBLIC :: copy_lnd_prog_now2new
  PUBLIC :: seaice_albedo_coldstart
  

CONTAINS

  !>
  !! Initialize soil temperature in boundary zone of LAM domain
  !!
  !! Soil tempertures in the boundary zone of the LAM domain are 
  !! filled with meaningful values. This has no immediate impact on the 
  !! prognostic results. It is, however, necessary in order to minimize 
  !! GRIB truncation errors.
  !!
  SUBROUTINE init_lamlatbc_phys (p_patch, p_prog_lnd_now, p_prog_lnd_new, p_lnd_diag)

    TYPE(t_patch)     , INTENT(IN)    :: p_patch       ! patch info
    TYPE(t_lnd_prog)  , INTENT(INOUT) :: p_prog_lnd_now, p_prog_lnd_new
    TYPE(t_lnd_diag)  , INTENT(INOUT) :: p_lnd_diag

    ! local variables
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jc

  !-------------------------------------------------------------------------

    ! only boundary cells
    rl_start = 1
    rl_end   = grf_bdywidth_c

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO jc= i_startidx, i_endidx
        ! tiled
        p_prog_lnd_now%t_g_t(jc,jb,:)    = tmelt
        p_prog_lnd_now%t_s_t(jc,jb,:)    = tmelt
        p_prog_lnd_now%t_sk_t(jc,jb,:)   = tmelt
        p_prog_lnd_now%t_so_t(jc,:,jb,:) = tmelt
        !
        p_prog_lnd_new%t_g_t(jc,jb,:)    = tmelt
        p_prog_lnd_new%t_s_t(jc,jb,:)    = tmelt
        p_prog_lnd_new%t_sk_t(jc,jb,:)   = tmelt
        p_prog_lnd_new%t_so_t(jc,:,jb,:) = tmelt
        !
        ! agg
        p_prog_lnd_now%t_g(jc,jb)    = tmelt
        p_prog_lnd_new%t_g(jc,jb)    = tmelt
        p_lnd_diag%t_s(jc,jb)        = tmelt
        p_lnd_diag%t_sk(jc,jb)       = tmelt
        p_lnd_diag%t_so(jc,:,jb)     = tmelt
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE init_lamlatbc_phys

  !-------------------------------------------------------------------------
  !! Init surface model TERRA, lake model Flake and seaice model
  SUBROUTINE nwp_surface_init( p_patch, ext_data, p_prog_lnd_now,           &
    &                          p_prog_lnd_new, p_prog_wtr_now,              &
    &                          p_prog_wtr_new, p_lnd_diag, p_diag, prm_diag )



    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch       !<grid/patch info.
    TYPE(t_external_data), INTENT(INOUT) :: ext_data
    TYPE(t_lnd_prog)     , INTENT(INOUT) :: p_prog_lnd_now, p_prog_lnd_new
    TYPE(t_wtr_prog)     , INTENT(INOUT) :: p_prog_wtr_now, p_prog_wtr_new
    TYPE(t_lnd_diag)     , INTENT(INOUT) :: p_lnd_diag
    TYPE(t_nh_diag), TARGET,INTENT(inout):: p_diag        !< diag vars
    TYPE(t_nwp_phy_diag),  INTENT(in)    :: prm_diag      !< atm phys vars

    ! Local array bounds:

    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< number of child domains

    ! Local scalars:

    INTEGER :: jc,jb,isubs,jk
    LOGICAL :: lsnowtile_warmstart


    REAL(wp) :: t_snow_now_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_snow_mult_now_t(nproma, 1:nlev_snow+1, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_rhosnowini_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_s_now_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_g_t    (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_s_new_t(nproma, p_patch%nblks_c, ntiles_total)

    REAL(wp) :: t_sk_now_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_sk_new_t(nproma, p_patch%nblks_c, ntiles_total)

    REAL(wp) :: w_snow_now_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: h_snow_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: rho_snow_now_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: rho_snow_mult_now_t(nproma, 1:nlev_snow, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_so_now_t(nproma, 1:nlev_soil+1, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: t_so_new_t(nproma, 1:nlev_soil+1, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: w_so_now_t(nproma, 1:nlev_soil  , p_patch%nblks_c, ntiles_total)
    REAL(wp) :: w_so_new_t(nproma, 1:nlev_soil  , p_patch%nblks_c, ntiles_total)
    REAL(wp) :: w_so_ice_now_t(nproma, 1:nlev_soil, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: w_so_ice_new_t(nproma, 1:nlev_soil, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: wliq_snow_now_t(nproma, 1:nlev_snow, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: wtot_snow_now_t(nproma, 1:nlev_snow, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: dzh_snow_now_t(nproma, 1:nlev_snow, p_patch%nblks_c, ntiles_total)

    INTEGER  :: soiltyp_t (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: rootdp_t  (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: plcov_t  (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: z0_t     (nproma, p_patch%nblks_c, ntiles_total)

    REAL(wp) :: freshsnow_t(nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: snowfrac_t (nproma, p_patch%nblks_c, ntiles_total)
    REAL(wp) :: sso_sigma_t(nproma,  p_patch%nblks_c, ntiles_total)
    INTEGER  :: lc_class_t (nproma,  p_patch%nblks_c, ntiles_total)


    ! local fields for lake model
    !
    REAL(wp) :: fr_lake       (nproma) ! lake fraction
    REAL(wp) :: depth_lk      (nproma) ! lake depth
    REAL(wp) :: fetch_lk      (nproma) ! wind fetch over lake
    REAL(wp) :: dp_bs_lk      (nproma) ! depth of thermally active layer of bot. sediments.
    REAL(wp) :: t_bs_lk       (nproma) ! clim. temp. at bottom of thermally active layer
                                       ! of sediments
    REAL(wp) :: gamso_lk      (nproma) ! attenuation coefficient of lake water with respect
                                       ! to sol. rad.
    REAL(wp) :: t_snow_lk_now (nproma) ! temperature of snow on lake ice
    REAL(wp) :: h_snow_lk_now (nproma) ! depth of snow on lake ice
    REAL(wp) :: t_ice_now     (nproma) ! lake ice temperature
    REAL(wp) :: h_ice_now     (nproma) ! lake ice depth
    REAL(wp) :: t_mnw_lk_now  (nproma) ! mean temperature of the water column
    REAL(wp) :: t_wml_lk_now  (nproma) ! mixed-layer temperature
    REAL(wp) :: t_bot_lk_now  (nproma) ! temperature at the water-bottom sediment interface
    REAL(wp) :: c_t_lk_now    (nproma) ! shape factor (temp. profile in lake thermocline)
    REAL(wp) :: h_ml_lk_now   (nproma) ! mixed-layer thickness
    REAL(wp) :: t_b1_lk_now   (nproma) ! temperature at the bottom of the upper layer
                                       ! of the sediments
    REAL(wp) :: h_b1_lk_now   (nproma) ! thickness of the upper layer of the sediments
    REAL(wp) :: t_scf_lk_now  (nproma) ! lake surface temperature

    LOGICAL  :: lake_mask(nproma)      ! auxiliary field for re-initialization of non-lake points
    LOGICAL  :: iceana_mask(nproma)    ! mask indicating if sea ice fraction analysis should be used in flake_init

    ! local fields for sea ice model
    !
    REAL(wp) :: frsi     (nproma)   ! sea ice fraction
    REAL(wp) :: tice_now (nproma)   ! temperature of ice upper surface at previous time
    REAL(wp) :: hice_now (nproma)   ! ice thickness at previous time level
    REAL(wp) :: tsnow_now(nproma)   ! temperature of snow upper surface at previous time
    REAL(wp) :: hsnow_now(nproma)   ! snow thickness at previous time level
    REAL(wp) :: albsi_now(nproma)   ! sea-ice albedo at previous time level
    REAL(wp) :: tice_new (nproma)   ! temperature of ice upper surface at new time
    REAL(wp) :: hice_new (nproma)   ! ice thickness at new time level
    REAL(wp) :: tsnow_new(nproma)   ! temperature of snow upper surface at new time
    REAL(wp) :: hsnow_new(nproma)   ! snow thickness at new time level
    REAL(wp) :: albsi_new(nproma)   ! sea-ice albedo at new time level

    ! local field for APT factors
    REAL(wp) :: snowfrac_fac(nproma, ntiles_total)

    INTEGER  :: icount_flk          ! total number of lake points per block
    !
    INTEGER  :: icount_ice          ! total number of sea-ice points per block

    INTEGER  :: i_count, ic, i_count_snow, isubs_snow, jg
    REAL(wp) :: temp, deglat, deglon
    REAL(wp) :: zfrice_thrhld       ! fraction threshold for creating a sea-ice grid point

  !-------------------------------------------------------------------------

    i_nchdom  = MAX(1,p_patch%n_childdom)
    jg        = p_patch%id

    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    IF (lsnowtile .AND. init_mode == MODE_IAU .AND. .NOT. ltile_coldstart) THEN
      lsnowtile_warmstart = .TRUE.
    ELSE
      lsnowtile_warmstart = .FALSE.
    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,isubs,i_count,i_count_snow,icount_ice, &
!$OMP            icount_flk,temp,ic,isubs_snow,frsi,tice_now,hice_now,               &
!$OMP            tsnow_now,hsnow_now,tice_new,hice_new,tsnow_new,hsnow_new,fr_lake,  &
!$OMP            depth_lk,fetch_lk,dp_bs_lk,t_bs_lk,gamso_lk,t_snow_lk_now,          &
!$OMP            h_snow_lk_now,t_ice_now,h_ice_now,t_mnw_lk_now,t_wml_lk_now,        &
!$OMP            t_bot_lk_now,c_t_lk_now,h_ml_lk_now,t_b1_lk_now,h_b1_lk_now,        &
!$OMP            t_scf_lk_now,zfrice_thrhld,lake_mask,albsi_now,albsi_new,           &
!$OMP            iceana_mask,deglat,deglon,snowfrac_fac), SCHEDULE(guided)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      IF (itopo == 1) THEN
        DO isubs = 1, ntiles_total
          DO jc = i_startidx, i_endidx

            ! initialize climatological layer (deepest layer of t_so)
            p_prog_lnd_now%t_so_t(jc,nlev_soil+1,jb,isubs) = ext_data%atm%t_cl(jc,jb)
            p_prog_lnd_new%t_so_t(jc,nlev_soil+1,jb,isubs) = ext_data%atm%t_cl(jc,jb)

          END DO
        END DO
      ENDIF

      ! t_s_t: initialization for open water and sea-ice tiles
      ! proper values are needed to perform surface analysis
      ! open water points: set it to SST
      ! sea-ice points   : set it to tf_salt (salt-water freezing point)
      !
      ! Note that after aggregation, t_s is copied to t_so(1)
      !
      DO ic = 1, ext_data%atm%list_seawtr%ncount(jb)
        jc = ext_data%atm%list_seawtr%idx(ic,jb)
        p_prog_lnd_now%t_s_t(jc,jb,isub_water) = p_lnd_diag%t_seasfc(jc,jb)
        p_prog_lnd_new%t_s_t(jc,jb,isub_water) = p_lnd_diag%t_seasfc(jc,jb)
        p_prog_lnd_now%t_sk_t(jc,jb,isub_water) = p_lnd_diag%t_seasfc(jc,jb)
        p_prog_lnd_new%t_sk_t(jc,jb,isub_water) = p_lnd_diag%t_seasfc(jc,jb)
      ENDDO

      DO ic = 1, ext_data%atm%list_seaice%ncount(jb)
        jc = ext_data%atm%list_seaice%idx(ic,jb)
        p_prog_lnd_now%t_s_t(jc,jb,isub_seaice) = tf_salt
        p_prog_lnd_new%t_s_t(jc,jb,isub_seaice) = tf_salt
        p_prog_lnd_now%t_sk_t(jc,jb,isub_seaice) = tf_salt
        p_prog_lnd_new%t_sk_t(jc,jb,isub_seaice) = tf_salt
      ENDDO


      ! Init t_g_t for sea water points
      !
      DO ic = 1, ext_data%atm%list_seawtr%ncount(jb)
        jc = ext_data%atm%list_seawtr%idx(ic,jb)
        temp =  p_lnd_diag%t_seasfc(jc,jb)
        p_prog_lnd_now%t_g_t(jc,jb,isub_water) = temp
        p_prog_lnd_new%t_g_t(jc,jb,isub_water) = temp
        ! includes reduction of saturation pressure due to salt content
        p_lnd_diag%qv_s_t(jc,jb,isub_water)    = salinity_fac * &  
          &       spec_humi(sat_pres_water(temp ),p_diag%pres_sfc(jc,jb) )
      END DO


!---------- Copy input fields for each tile

!----------------------------------
      DO isubs = 1,ntiles_total
!----------------------------------

        i_count = ext_data%atm%lp_count_t(jb,isubs)

        IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty

        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
          t_snow_now_t(ic,jb,isubs)          =  p_prog_lnd_now%t_snow_t(jc,jb,isubs)
          t_s_now_t(ic,jb,isubs)             =  p_prog_lnd_now%t_s_t(jc,jb,isubs)
          t_s_new_t(ic,jb,isubs)             =  p_prog_lnd_new%t_s_t(jc,jb,isubs)
          t_sk_now_t(ic,jb,isubs)            =  p_prog_lnd_now%t_sk_t(jc,jb,isubs)
          t_sk_new_t(ic,jb,isubs)            =  p_prog_lnd_new%t_sk_t(jc,jb,isubs)

          w_snow_now_t(ic,jb,isubs)          =  p_prog_lnd_now%w_snow_t(jc,jb,isubs)
          rho_snow_now_t(ic,jb,isubs)        =  p_prog_lnd_now%rho_snow_t(jc,jb,isubs)

          sso_sigma_t(ic,jb,isubs)           = ext_data%atm%sso_stdh(jc,jb)
          lc_class_t(ic,jb,isubs)            = ext_data%atm%lc_class_t(jc,jb,isubs)
          freshsnow_t(ic,jb,isubs)           = p_lnd_diag%freshsnow_t(jc,jb,isubs)
          h_snow_t(ic,jb,isubs)              = p_lnd_diag%h_snow_t(jc,jb,isubs)

          soiltyp_t(ic,jb,isubs)             =  ext_data%atm%soiltyp_t(jc,jb,isubs)
          rootdp_t(ic,jb,isubs)              =  ext_data%atm%rootdp_t(jc,jb,isubs)
          plcov_t(ic,jb,isubs)               =  ext_data%atm%plcov_t(jc,jb,isubs)

          z0_t(ic,jb,isubs)                  =  prm_diag%gz0_t(jc,jb,isubs)/grav

          IF (itype_vegetation_cycle >= 2) THEN ! use climatological temperature to specify snow density on glaciers
            t_rhosnowini_t(ic,jb,isubs)      = ext_data%atm%t2m_clim_hc(jc,jb)
          ELSE
            t_rhosnowini_t(ic,jb,isubs)      = t_snow_now_t(ic,jb,isubs)
          ENDIF

          IF (icpl_da_snowalb >= 3) THEN
            snowfrac_fac(ic,isubs)  = prm_diag%snowfrac_fac(jc,jb)
          ELSE
            snowfrac_fac(ic,isubs)  = 1._wp
          ENDIF
        ENDDO

        IF(l2lay_rho_snow .OR. lmulti_snow) THEN
          DO jk=1,nlev_snow
            DO ic = 1, i_count
              jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
              rho_snow_mult_now_t(ic,jk,jb,isubs) =  p_prog_lnd_now%rho_snow_mult_t(jc,jk,jb,isubs)
            ENDDO
          ENDDO

        ENDIF

        IMSNOWI: IF(lmulti_snow) THEN

          DO jk=1,nlev_snow+1
            DO ic = 1, i_count
              jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
              t_snow_mult_now_t(ic,jk,jb,isubs)   =  p_prog_lnd_now%t_snow_mult_t(jc,jk,jb,isubs)
            ENDDO
          ENDDO

          DO jk=1,nlev_snow
            DO ic = 1, i_count
              jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
              wliq_snow_now_t(ic,jk,jb,isubs)     =  p_prog_lnd_now%wliq_snow_t    (jc,jk,jb,isubs)
              wtot_snow_now_t(ic,jk,jb,isubs)     =  p_prog_lnd_now%wtot_snow_t    (jc,jk,jb,isubs)
              dzh_snow_now_t(ic,jk,jb,isubs)      =  p_prog_lnd_now%dzh_snow_t     (jc,jk,jb,isubs)
            ENDDO
          ENDDO

        END IF  IMSNOWI

!$NEC outerloop_unroll(nlsoil+1)
        DO jk=1,nlev_soil+1
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
            t_so_now_t(ic,jk,jb,isubs)          =  p_prog_lnd_now%t_so_t(jc,jk,jb,isubs)
            t_so_new_t(ic,jk,jb,isubs)          =  p_prog_lnd_new%t_so_t(jc,jk,jb,isubs)
          ENDDO
        ENDDO

!$NEC outerloop_unroll(nlsoil)
        DO jk=1,nlev_soil
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
            w_so_now_t(ic,jk,jb,isubs)          =  p_prog_lnd_now%w_so_t(jc,jk,jb,isubs)
            w_so_new_t(ic,jk,jb,isubs)          =  p_prog_lnd_new%w_so_t(jc,jk,jb,isubs)
            w_so_ice_now_t(ic,jk,jb,isubs)      =  p_prog_lnd_now%w_so_ice_t(jc,jk,jb,isubs)
            w_so_ice_new_t(ic,jk,jb,isubs)      =  p_prog_lnd_new%w_so_ice_t(jc,jk,jb,isubs)
          ENDDO
        ENDDO

        ! Preliminary diagnosis of snow-cover fraction for initialization of split tile index list

        IF (.NOT. lsnowtile_warmstart) THEN

          CALL diag_snowfrac_tg(                           &
            &  istart = 1, iend = i_count                , & ! start/end indices
            &  lc_class  = lc_class_t        (:,jb,isubs), & ! land-cover class
            &  i_lc_urban = ext_data%atm%i_lc_urban      , & ! land-cover class index for urban areas
            &  t_snow    = t_snow_now_t      (:,jb,isubs), & ! snow temp
            &  t_soiltop = t_s_now_t         (:,jb,isubs), & ! soil top temp
            &  w_snow    = w_snow_now_t      (:,jb,isubs), & ! snow WE
            &  rho_snow  = rho_snow_now_t    (:,jb,isubs), & ! snow density
            &  freshsnow = freshsnow_t       (:,jb,isubs), & ! fresh snow fraction
            &  sso_sigma = sso_sigma_t       (:,jb,isubs), & ! sso stdev
            &  z0        = z0_t              (:,jb,isubs), & ! vegetation roughness length
            &  snowfrac_fac = snowfrac_fac   (:,   isubs), & ! APT tuning factor for snow-cover fraction
            &  snowfrac  = snowfrac_t        (:,jb,isubs), & ! OUT: snow cover fraction
            &  t_g       = t_g_t             (:,jb,isubs)  ) ! OUT: averaged ground temp

          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
            p_lnd_diag%snowfrac_lc_t(jc,jb,isubs)  = snowfrac_t(ic,jb,isubs)
            p_prog_lnd_now%t_g_t(jc,jb,isubs)      = t_g_t(ic,jb,isubs)
          ENDDO

        ENDIF

      ENDDO

      ! create index lists for snow tiles (first call)
      IF(lsnowtile .AND. .NOT. lsnowtile_warmstart) THEN ! snow is considered as separate tiles
        DO isubs = 1, ntiles_lnd
          isubs_snow = isubs + ntiles_lnd

          CALL update_idx_lists_lnd (idx_lst_lp       = ext_data%atm%idx_lst_lp_t(:,jb,isubs),         &
                                   lp_count           = ext_data%atm%lp_count_t(jb,isubs),             &
                                   idx_lst            = ext_data%atm%idx_lst_t(:,jb,isubs),            &
                                   gp_count           = ext_data%atm%gp_count_t(jb,isubs),             &
                                   idx_lst_snow       = ext_data%atm%idx_lst_t(:,jb,isubs_snow),       &
                                   gp_count_snow      = ext_data%atm%gp_count_t(jb,isubs_snow),        &
                                   lc_frac            = ext_data%atm%lc_frac_t(:,jb,isubs),            &
                                   partial_frac       = ext_data%atm%frac_t(:,jb,isubs),               &
                                   partial_frac_snow  = ext_data%atm%frac_t(:,jb,isubs_snow),          &
                                   snowtile_flag      = ext_data%atm%snowtile_flag_t(:,jb,isubs),      &
                                   snowtile_flag_snow = ext_data%atm%snowtile_flag_t(:,jb,isubs_snow), &
                                   snowfrac           = p_lnd_diag%snowfrac_lc_t(:,jb,isubs)           )

        END DO
      END IF


      DO isubs = 1,ntiles_total

        i_count = ext_data%atm%lp_count_t(jb,isubs)
        CALL terra_init(                                          &
        &  init_mode         = init_mode_soil                   , & ! coldstart/warmstart/warmstart with snow increments
        &  nvec              = nproma                           , & ! array dimensions
        &  ivstart=1,   ivend= i_count                          , & ! optional start/end indicies
        &  iblock            = jb                               , & ! actual block
        &  ke_soil=nlev_soil-1, ke_snow=nlev_snow               , & ! without lowermost (climat.) soil layer
        &   zmls             = zml_soil                         , & ! processing soil level structure
        &  soiltyp_subs      = soiltyp_t(:,jb,isubs)            , & ! type of the soil (keys 0-9)  --
        &  rootdp            = rootdp_t(:,jb,isubs)             , & ! depth of the roots                ( m  )
        &  plcov             = plcov_t(:,jb,isubs)              , & ! surface fraction covered by plants ( -  )
        &  t_snow_now        = t_snow_now_t(:,jb,isubs)         , & ! temperature of the snow-surface   (  K  )
        &  t_snow_mult_now   = t_snow_mult_now_t(:,:,jb,isubs)  , & ! temperature of the snow-surface   (  K  )
        &  t_rhosnowini      = t_rhosnowini_t(:,jb,isubs)       , & ! temperature used for snow density initialization (  K  )
        &  t_s_now           = t_s_now_t(:,jb,isubs)            , & ! temperature of the ground surface (  K  )
        &  t_s_new           = t_s_new_t(:,jb,isubs)            , & ! temperature of the ground surface (  K  )
        &  w_snow_now        = w_snow_now_t(:,jb,isubs)         , & ! water content of snow             (m H2O)
        &  h_snow            = h_snow_t(:,jb,isubs)             , & ! snow depth                        (m H2O)
        &  rho_snow_now      = rho_snow_now_t(:,jb,isubs)       , & ! snow density                      (kg/m**3)
        &  rho_snow_mult_now = rho_snow_mult_now_t(:,:,jb,isubs), & ! snow density                      (kg/m**3)
        &  t_so_now          = t_so_now_t(:,:,jb,isubs)         , & ! soil temperature (main level)     (  K  )
        &  t_so_new          = t_so_new_t(:,:,jb,isubs)         , & ! soil temperature (main level)     (  K  )
        &  w_so_now          = w_so_now_t(:,:,jb,isubs)         , & ! total water content (ice + liquid water)     (m H20)
        &  w_so_new          = w_so_new_t(:,:,jb,isubs)         , & ! total water content (ice + liquid water)     (m H20)
        &  w_so_ice_now      = w_so_ice_now_t(:,:,jb,isubs)     , & ! ice content                       (m H20)
        &  w_so_ice_new      = w_so_ice_new_t(:,:,jb,isubs)     , & ! ice content                       (m H20)
        &  wliq_snow_now     = wliq_snow_now_t(:,:,jb,isubs)    , & ! liquid water content in the snow  (m H2O)
        &  wtot_snow_now     = wtot_snow_now_t(:,:,jb,isubs)    , & ! total (liquid + solid) water content of snow (m H2O)
        &  dzh_snow_now      = dzh_snow_now_t(:,:,jb,isubs)       ) ! layer thickness between half levels in snow  (  m  )

        IF (.NOT. lsnowtile_warmstart .OR. isubs > ntiles_lnd) THEN

          CALL diag_snowfrac_tg(                                  &
            &  istart     = 1, iend = i_count                   , & ! start/end indices
            &  lc_class   = lc_class_t              (:,jb,isubs), & ! land-cover class
            &  i_lc_urban = ext_data%atm%i_lc_urban             , & ! land-cover class index for urban areas
            &  t_snow     = t_snow_now_t            (:,jb,isubs), & ! snow temperature
            &  t_soiltop  = t_sk_now_t              (:,jb,isubs), & ! soil top temperature or skin temperature
            &  w_snow     = w_snow_now_t            (:,jb,isubs), & ! snow WE
            &  rho_snow   = rho_snow_now_t          (:,jb,isubs), & ! snow depth
            &  freshsnow  = freshsnow_t             (:,jb,isubs), & ! fresh snow fraction
            &  sso_sigma  = sso_sigma_t             (:,jb,isubs), & ! sso stdev
            &  z0         = z0_t                    (:,jb,isubs), & ! vegetation roughness length
            &  snowfrac_fac = snowfrac_fac          (:,   isubs), & ! APT tuning factor for snow-cover fraction
            &  snowfrac   = snowfrac_t              (:,jb,isubs), & ! OUT: snow cover fraction
            &  t_g        = t_g_t                   (:,jb,isubs)  ) ! OUT: averaged ground temperature

        ENDIF

!  Recover fields from index list
!
        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
          p_prog_lnd_now%t_snow_t(jc,jb,isubs)   = t_snow_now_t(ic,jb,isubs)
          p_prog_lnd_now%t_s_t(jc,jb,isubs)      = t_s_now_t(ic,jb,isubs)
          p_prog_lnd_new%t_s_t(jc,jb,isubs)      = t_s_new_t(ic,jb,isubs)
          p_prog_lnd_now%t_sk_t(jc,jb,isubs)     = t_sk_now_t(ic,jb,isubs)
          p_prog_lnd_new%t_sk_t(jc,jb,isubs)     = t_sk_new_t(ic,jb,isubs)
          p_prog_lnd_now%w_snow_t(jc,jb,isubs)   = w_snow_now_t(ic,jb,isubs)
          p_lnd_diag%h_snow_t(jc,jb,isubs)       = h_snow_t(ic,jb,isubs)
          p_prog_lnd_now%rho_snow_t(jc,jb,isubs) = rho_snow_now_t(ic,jb,isubs)
          IF (.NOT. lsnowtile_warmstart) THEN
            p_lnd_diag%snowfrac_lc_t(jc,jb,isubs)  = snowfrac_t(ic,jb,isubs)
            p_lnd_diag%snowfrac_t(jc,jb,isubs)     = snowfrac_t(ic,jb,isubs)
            p_prog_lnd_now%t_g_t(jc,jb,isubs)      = t_g_t(ic,jb,isubs)
            p_prog_lnd_new%t_g_t(jc,jb,isubs)      = t_g_t(ic,jb,isubs)
          ENDIF
          IF (isubs > ntiles_lnd) THEN
            p_lnd_diag%snowfrac_lcu_t(jc,jb,isubs) = MAX(snowfrac_t(ic,jb,isubs),p_lnd_diag%snowfrac_lc_t(jc,jb,isubs))
            p_lnd_diag%snowfrac_lcu_t(jc,jb,isubs-ntiles_lnd) = p_lnd_diag%snowfrac_lcu_t(jc,jb,isubs)
          ENDIF
        ENDDO

        IF (lsnowtile .AND. isubs > ntiles_lnd .AND. .NOT. lsnowtile_warmstart) THEN
          ! copy snowfrac_t to snow-free tile (needed for index list computation)
!$NEC ivdep
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
            p_lnd_diag%snowfrac_lc_t(jc,jb,isubs-ntiles_lnd) = p_lnd_diag%snowfrac_lc_t(jc,jb,isubs)
          ENDDO
!$NEC ivdep
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
            p_prog_lnd_now%w_snow_t(jc,jb,isubs) = p_prog_lnd_now%w_snow_t(jc,jb,isubs)/            &
                                                   MAX(0.01_wp,p_lnd_diag%snowfrac_lc_t(jc,jb,isubs))
            p_lnd_diag%h_snow_t(jc,jb,isubs)     = p_lnd_diag%h_snow_t(jc,jb,isubs)/                &
                                                   MAX(0.01_wp,p_lnd_diag%snowfrac_lc_t(jc,jb,isubs))
          ENDDO

        ENDIF

        IF (itype_snowevap == 3) THEN ! set snow age to upper limit of 365 days on grid points with glaciers,
!$NEC ivdep                           ! and reset hsnow_max to 40 m in agreement with what is done in initicon
         DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
            IF (ext_data%atm%lc_class_t(jc,jb,isubs) == ext_data%atm%i_lc_snow_ice) &
              p_lnd_diag%snow_age(jc,jb)  = 365._wp
              p_lnd_diag%hsnow_max(jc,jb) = MIN(40._wp,p_lnd_diag%hsnow_max(jc,jb))
          ENDDO
        ENDIF

        IF(l2lay_rho_snow .OR. lmulti_snow) THEN
          DO jk=1,nlev_snow
            DO ic = 1, i_count
              jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
              p_prog_lnd_now%rho_snow_mult_t(jc,jk,jb,isubs) = rho_snow_mult_now_t(ic,jk,jb,isubs)
            ENDDO
          ENDDO
        ENDIF

        IMSNOWO: IF(lmulti_snow) THEN

          DO jk=1,nlev_snow+1
            DO ic = 1, i_count
              jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
              p_prog_lnd_now%t_snow_mult_t(jc,jk,jb,isubs) =  t_snow_mult_now_t(ic,jk,jb,isubs)
            ENDDO
          ENDDO

          DO jk=1,nlev_snow
            IF (lsnowtile .AND. isubs > ntiles_lnd .AND. .NOT. lsnowtile_warmstart) THEN
              DO ic = 1, i_count
                jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
                p_prog_lnd_now%wliq_snow_t(jc,jk,jb,isubs) = wliq_snow_now_t(ic,jk,jb,isubs)/ &
                                             MAX(0.01_wp,p_lnd_diag%snowfrac_lc_t(jc,jb,isubs))
                p_prog_lnd_now%wtot_snow_t(jc,jk,jb,isubs) = wtot_snow_now_t(ic,jk,jb,isubs)/ &
                                             MAX(0.01_wp,p_lnd_diag%snowfrac_lc_t(jc,jb,isubs))
                p_prog_lnd_now%dzh_snow_t(jc,jk,jb,isubs)  = dzh_snow_now_t(ic,jk,jb,isubs)/ &
                                            MAX(0.01_wp,p_lnd_diag%snowfrac_lc_t(jc,jb,isubs))
              ENDDO
            ELSE
              DO ic = 1, i_count
                jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
                p_prog_lnd_now%wliq_snow_t(jc,jk,jb,isubs) = wliq_snow_now_t(ic,jk,jb,isubs)
                p_prog_lnd_now%wtot_snow_t(jc,jk,jb,isubs) = wtot_snow_now_t(ic,jk,jb,isubs)
                p_prog_lnd_now%dzh_snow_t(jc,jk,jb,isubs)  = dzh_snow_now_t(ic,jk,jb,isubs)
              ENDDO
            ENDIF
          ENDDO

        END IF  IMSNOWO

!$NEC outerloop_unroll(nlsoil+1)
        DO jk=1,nlev_soil+1
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
            p_prog_lnd_now%t_so_t(jc,jk,jb,isubs) = t_so_now_t(ic,jk,jb,isubs)
            p_prog_lnd_new%t_so_t(jc,jk,jb,isubs) = t_so_new_t(ic,jk,jb,isubs)
          ENDDO
        ENDDO

!$NEC outerloop_unroll(nlsoil)
        DO jk=1,nlev_soil
          DO ic = 1, i_count
            jc = ext_data%atm%idx_lst_lp_t(ic,jb,isubs)
            p_prog_lnd_now%w_so_t(jc,jk,jb,isubs) = w_so_now_t(ic,jk,jb,isubs)
            p_prog_lnd_new%w_so_t(jc,jk,jb,isubs) = w_so_new_t(ic,jk,jb,isubs)
            p_prog_lnd_now%w_so_ice_t(jc,jk,jb,isubs) = w_so_ice_now_t(ic,jk,jb,isubs)
            p_prog_lnd_new%w_so_ice_t(jc,jk,jb,isubs) = w_so_ice_new_t(ic,jk,jb,isubs)
          ENDDO
        ENDDO
      END DO ! isubs


      !===================================================!
      !                                                   !
      !  Warm-start initialization for sea-ice and lake   !
      !                                                   !
      !===================================================!

      !
      ! Warm-start initialization for fresh water lake parameterization
      ! This initialization is performed irrespective of a cold-start initialization.
      !
      IF (llake) THEN

        icount_flk = ext_data%atm%list_lake%ncount(jb) ! number of lake points in block jb

        ! Collect data for lake points in 1D-arrays
        DO ic = 1, icount_flk

          jc = ext_data%atm%list_lake%idx(ic,jb)

          fr_lake      (ic) = ext_data%atm%frac_t(jc,jb,isub_lake)
          depth_lk     (ic) = ext_data%atm%depth_lk   (jc,jb)
          frsi         (ic) = p_lnd_diag%fr_seaice(jc,jb)
          fetch_lk     (ic) = ext_data%atm%fetch_lk   (jc,jb)
          dp_bs_lk     (ic) = ext_data%atm%dp_bs_lk   (jc,jb)
          t_bs_lk      (ic) = ext_data%atm%t_bs_lk    (jc,jb)
          gamso_lk     (ic) = ext_data%atm%gamso_lk   (jc,jb)
          t_snow_lk_now(ic) = p_prog_wtr_now%t_snow_lk(jc,jb)
          h_snow_lk_now(ic) = p_prog_wtr_now%h_snow_lk(jc,jb)
          t_ice_now    (ic) = p_prog_wtr_now%t_ice    (jc,jb)
          h_ice_now    (ic) = p_prog_wtr_now%h_ice    (jc,jb)
          t_mnw_lk_now (ic) = p_prog_wtr_now%t_mnw_lk (jc,jb)
          t_wml_lk_now (ic) = p_prog_wtr_now%t_wml_lk (jc,jb)
          t_bot_lk_now (ic) = p_prog_wtr_now%t_bot_lk (jc,jb)
          c_t_lk_now   (ic) = p_prog_wtr_now%c_t_lk   (jc,jb)
          h_ml_lk_now  (ic) = p_prog_wtr_now%h_ml_lk  (jc,jb)
          t_b1_lk_now  (ic) = p_prog_wtr_now%t_b1_lk  (jc,jb)
          h_b1_lk_now  (ic) = p_prog_wtr_now%h_b1_lk  (jc,jb)
          t_scf_lk_now (ic) = p_prog_lnd_now%t_g_t    (jc,jb,isub_lake)

          ! test implementation: compute approximate mask for Great Lakes
          deglat = p_patch%cells%center(jc,jb)%lat * rad2deg
          deglon = p_patch%cells%center(jc,jb)%lon * rad2deg
          IF (deglat >= 41._wp .AND. deglat <= 49._wp .AND. deglon >= -92._wp .AND. deglon <= -76._wp) THEN
            iceana_mask(ic) = .TRUE.
            IF (deglon > -85._wp .AND. deglat > 46._wp) iceana_mask(ic) = .FALSE.
            IF (deglon < -88._wp .AND. deglat < 46._wp) iceana_mask(ic) = .FALSE.
          ELSE
            iceana_mask(ic) = .FALSE.
          ENDIF
          IF (.NOT. (use_lakeiceana .AND. lanaread_tseasfc(jg)) ) iceana_mask(ic) = .FALSE.
        ENDDO


        CALL flake_init (nflkgb     = icount_flk,             & !in
          &         use_iceanalysis = iceana_mask  (:),       & !in
          &              fr_lake    = fr_lake      (:),       & !in
          &              depth_lk   = depth_lk     (:),       & !in
          &              fr_ice     = frsi         (:),       & !in
          &              fetch_lk   = fetch_lk     (:),       & !inout
          &              dp_bs_lk   = dp_bs_lk     (:),       & !inout
          &              t_bs_lk    = t_bs_lk      (:),       & !inout
          &              gamso_lk   = gamso_lk     (:),       & !inout
          &              t_snow_p   = t_snow_lk_now(:),       & !inout
          &              h_snow_p   = h_snow_lk_now(:),       & !inout
          &              t_ice_p    = t_ice_now    (:),       & !inout
          &              h_ice_p    = h_ice_now    (:),       & !inout
          &              t_mnw_lk_p = t_mnw_lk_now (:),       & !inout
          &              t_wml_lk_p = t_wml_lk_now (:),       & !inout
          &              t_bot_lk_p = t_bot_lk_now (:),       & !inout
          &              c_t_lk_p   = c_t_lk_now   (:),       & !inout
          &              h_ml_lk_p  = h_ml_lk_now  (:),       & !inout
          &              t_b1_lk_p  = t_b1_lk_now  (:),       & !inout
          &              h_b1_lk_p  = h_b1_lk_now  (:),       & !inout
          &              t_scf_lk_p = t_scf_lk_now (:)        )


        !  Recover fields from index list
        !
!$NEC ivdep
        DO ic = 1, icount_flk

          jc = ext_data%atm%list_lake%idx(ic,jb)

          ext_data%atm%fetch_lk(jc,jb) = fetch_lk(ic)
          ext_data%atm%dp_bs_lk(jc,jb) = dp_bs_lk(ic)
          ext_data%atm%t_bs_lk (jc,jb) = t_bs_lk (ic)
          ext_data%atm%gamso_lk(jc,jb) = gamso_lk(ic)

          p_prog_wtr_now%t_snow_lk(jc,jb) = t_snow_lk_now(ic)
          p_prog_wtr_now%h_snow_lk(jc,jb) = h_snow_lk_now(ic)
          p_prog_wtr_now%t_ice    (jc,jb) = t_ice_now    (ic)
          p_prog_wtr_now%h_ice    (jc,jb) = h_ice_now    (ic)
          p_prog_wtr_now%t_mnw_lk (jc,jb) = t_mnw_lk_now (ic)
          p_prog_wtr_now%t_wml_lk (jc,jb) = t_wml_lk_now (ic)
          p_prog_wtr_now%t_bot_lk (jc,jb) = t_bot_lk_now (ic)
          p_prog_wtr_now%c_t_lk   (jc,jb) = c_t_lk_now   (ic)
          p_prog_wtr_now%h_ml_lk  (jc,jb) = h_ml_lk_now  (ic)
          p_prog_wtr_now%t_b1_lk  (jc,jb) = t_b1_lk_now  (ic)
          p_prog_wtr_now%h_b1_lk  (jc,jb) = h_b1_lk_now  (ic)

          p_prog_lnd_now%t_g_t(jc,jb,isub_lake) = t_scf_lk_now(ic)

          ! for consistency, set 
          ! t_so(0) = t_wml_lk       mixed-layer temperature (273.15K if the lake is frozen)
          p_prog_lnd_now%t_s_t(jc,jb,isub_lake) = p_prog_wtr_now%t_wml_lk (jc,jb)
          p_prog_lnd_new%t_s_t(jc,jb,isub_lake) = p_prog_lnd_now%t_s_t(jc,jb,isub_lake)
          p_prog_lnd_now%t_sk_t(jc,jb,isub_lake) = p_prog_wtr_now%t_wml_lk (jc,jb)
          p_prog_lnd_new%t_sk_t(jc,jb,isub_lake) = p_prog_lnd_now%t_s_t(jc,jb,isub_lake)


          ! In addition, initialize prognostic Flake fields at time step 'new'
          p_prog_wtr_new%t_snow_lk(jc,jb) = t_snow_lk_now(ic)
          p_prog_wtr_new%h_snow_lk(jc,jb) = h_snow_lk_now(ic)
          p_prog_wtr_new%t_ice    (jc,jb) = t_ice_now    (ic)
          p_prog_wtr_new%h_ice    (jc,jb) = h_ice_now    (ic)
          p_prog_wtr_new%t_mnw_lk (jc,jb) = t_mnw_lk_now (ic)
          p_prog_wtr_new%t_wml_lk (jc,jb) = t_wml_lk_now (ic)
          p_prog_wtr_new%t_bot_lk (jc,jb) = t_bot_lk_now (ic)
          p_prog_wtr_new%c_t_lk   (jc,jb) = c_t_lk_now   (ic)
          p_prog_wtr_new%h_ml_lk  (jc,jb) = h_ml_lk_now  (ic)
          p_prog_wtr_new%t_b1_lk  (jc,jb) = t_b1_lk_now  (ic)
          p_prog_wtr_new%h_b1_lk  (jc,jb) = h_b1_lk_now  (ic)

          ! keep fr_seaice synchronized with h_ice
          ! i.e. set fr_seaice=1 for frozen lakes
          p_lnd_diag%fr_seaice(jc,jb) = MERGE(1.0_wp, 0.0_wp, &
            &                           p_prog_wtr_now%h_ice(jc,jb) > csmall_hice)
        ENDDO  ! ic

        ! Re-Initialize lake-specific fields
        ! t_wml_lk
        ! t_mnw_lk
        ! t_bot_lk
        ! c_t_lk
        ! h_ml_lk
        ! h_snow_lk
        ! t_snow_lk
        ! t_bl_lk
        ! h_bl_lk
        ! for non-lake points, because values might be inconsistent due to GRIB-packing
        !
        ! Create lake-mask, in order to re-initialize only non-lake points
        DO jc = i_startidx, i_endidx
          lake_mask(jc) = .FALSE.
        END DO
        ! set lake-mask to .TRUE. for lake points
        DO ic = 1, icount_flk
          jc = ext_data%atm%list_lake%idx(ic,jb)
          lake_mask(jc) = .TRUE.
        ENDDO

        DO jc = i_startidx, i_endidx
          IF (.NOT. lake_mask(jc)) THEN
            ! now
            p_prog_wtr_now%t_wml_lk (jc,jb) = tpl_T_r        ! temperature of maximum density
                                                             ! of fresh water
            p_prog_wtr_now%t_mnw_lk (jc,jb) = tpl_T_r
            p_prog_wtr_now%t_bot_lk (jc,jb) = tpl_T_r
            p_prog_wtr_now%c_t_lk   (jc,jb) = C_T_min        ! minimum value
            p_prog_wtr_now%h_ml_lk  (jc,jb) = 0._wp
            p_prog_wtr_now%h_snow_lk(jc,jb) = 0._wp
            p_prog_wtr_now%t_snow_lk(jc,jb) = tmelt          ! fresh water freezing point
            p_prog_wtr_now%t_b1_lk  (jc,jb) = tpl_T_r        ! reference value, bottom-sediment
                                                             ! module is switched off
            p_prog_wtr_now%h_b1_lk  (jc,jb) = rflk_depth_bs_ref ! reference value, bottom-sediment
                                                                ! is switched off
            !
            ! new
            p_prog_wtr_new%t_wml_lk (jc,jb) = tpl_T_r        ! temperature of maximum density
                                                             ! of fresh water
            p_prog_wtr_new%t_mnw_lk (jc,jb) = tpl_T_r
            p_prog_wtr_new%t_bot_lk (jc,jb) = tpl_T_r
            p_prog_wtr_new%c_t_lk   (jc,jb) = C_T_min        ! minimum value
            p_prog_wtr_new%h_ml_lk  (jc,jb) = 0._wp
            p_prog_wtr_new%h_snow_lk(jc,jb) = 0._wp
            p_prog_wtr_new%t_snow_lk(jc,jb) = tmelt          ! fresh water freezing point
            p_prog_wtr_new%t_b1_lk  (jc,jb) = tpl_T_r        ! reference value, bottom-sediment
                                                             ! module is switched off
            p_prog_wtr_new%h_b1_lk  (jc,jb) = rflk_depth_bs_ref ! reference value, bottom-sediment
                                                                ! is switched off
          ENDIF
        ENDDO

      ENDIF  ! llake


      !
      ! Warm-start for sea-ice parameterization scheme
      ! This initialization is performed irrespectively of a coldstart initialization.
      ! At this point it is assumed that both the ice thickness and the ice surface
      ! temperature are available from the previous ICON run or any other educated guess.
      !
      IF ( lseaice ) THEN

        icount_ice = ext_data%atm%list_seaice%ncount(jb) ! number of sea-ice points in block jb

        DO ic = 1, icount_ice
          jc = ext_data%atm%list_seaice%idx(ic,jb)

          frsi     (ic) = p_lnd_diag%fr_seaice(jc,jb)
          tice_now (ic) = p_prog_wtr_now%t_ice(jc,jb)
          hice_now (ic) = p_prog_wtr_now%h_ice(jc,jb)
          tsnow_now(ic) = p_prog_wtr_now%t_snow_si(jc,jb)
          hsnow_now(ic) = p_prog_wtr_now%h_snow_si(jc,jb)
          albsi_now(ic) = p_prog_wtr_now%alb_si(jc,jb)
        ENDDO  ! jc

        ! The first argument of `seaice_init_nwp` is set .TRUE. 
        ! to make sure that the ice thickness and ice surface temperature 
        ! at "new" ice points are initialized for all run types 
        ! (for non-coupled as well as coupled runs).
        !
        CALL seaice_init_nwp ( .TRUE., icount_ice, frsi,                            & ! in
          &                    tice_now, hice_now, tsnow_now, hsnow_now, albsi_now, & ! inout
          &                    tice_new, hice_new, tsnow_new, hsnow_new, albsi_new  ) ! inout


        !  Recover fields from index list
        !
        DO ic = 1, icount_ice
          jc = ext_data%atm%list_seaice%idx(ic,jb)

          ! fields at time level now may have changed, potentially!
          p_prog_wtr_now%t_ice(jc,jb)    = tice_now(ic)
          p_prog_wtr_now%h_ice(jc,jb)    = hice_now(ic)
          p_prog_wtr_now%t_snow_si(jc,jb)= tsnow_now(ic)
          p_prog_wtr_now%h_snow_si(jc,jb)= hsnow_now(ic)
          p_prog_wtr_now%alb_si(jc,jb)   = albsi_now(ic)

          p_prog_wtr_new%t_ice(jc,jb)    = tice_new(ic)
          p_prog_wtr_new%h_ice(jc,jb)    = hice_new(ic)
          p_prog_wtr_new%t_snow_si(jc,jb)= tsnow_new(ic)
          p_prog_wtr_new%h_snow_si(jc,jb)= hsnow_new(ic)
          IF (lprog_albsi) THEN
            p_prog_wtr_new%alb_si(jc,jb)   = albsi_new(ic)
          ENDIF

          p_prog_lnd_now%t_g_t(jc,jb,isub_seaice) =  tice_now(ic)
          p_prog_lnd_new%t_g_t(jc,jb,isub_seaice) =  tice_new(ic) ! == tice_now(ic)
          p_lnd_diag%qv_s_t(jc,jb,isub_seaice)    = spec_humi(sat_pres_ice(tice_now(ic)),&
          &                                   p_diag%pres_sfc(jc,jb) )
        ENDDO  ! ic


        ! Re-initialize sea-ice specific fields h_snow_si, t_snow_si for non sea-ice points,
        ! because values might be inconsistent due to GRIB-packing
        !
        IF ( ntiles_total == 1 ) THEN  ! no tile approach
          zfrice_thrhld = 0.5_wp
        ELSE
          zfrice_thrhld = frsi_min
        ENDIF
        DO jc = i_startidx, i_endidx
          IF (p_lnd_diag%fr_seaice(jc,jb) < zfrice_thrhld) THEN   ! non-seaice point
            ! now
            p_prog_wtr_now%h_snow_si(jc,jb) = 0._wp
            p_prog_wtr_now%t_snow_si(jc,jb) = tmelt   ! snow over ice is not treated explicitly
            !
            ! new
            p_prog_wtr_new%h_snow_si(jc,jb) = 0._wp
            p_prog_wtr_new%t_snow_si(jc,jb) = tmelt   ! snow over ice is not treated explicitly
          ENDIF
        ENDDO  ! jc

      ENDIF  ! lseaice


      ! Re-initialize h_ice and t_ice which are used by both the sea-ice and lake scheme.
      ! 4 cases need to be distinguished in order to make sure, that none of the prognostic
      ! points is overwritten.
      IF      (.NOT. lseaice .AND. .NOT. llake) THEN

        ! Re-initialize all points
        !
        DO jc = i_startidx, i_endidx
          ! now
          p_prog_wtr_now%h_ice  (jc,jb) = 0._wp
          p_prog_wtr_now%t_ice  (jc,jb) = tmelt   ! fresh water freezing point
          !
          ! new
          p_prog_wtr_new%h_ice    (jc,jb) = 0._wp
          p_prog_wtr_new%t_ice    (jc,jb) = tmelt   ! fresh water freezing point
        ENDDO

      ELSE IF (      lseaice .AND. .NOT. llake) THEN

        ! Re-initialize all non-sea-ice points
        !
        IF ( ntiles_total == 1 ) THEN  ! no tile approach
          zfrice_thrhld = 0.5_wp
        ELSE
          zfrice_thrhld = frsi_min
        ENDIF
        DO jc = i_startidx, i_endidx
          IF (p_lnd_diag%fr_seaice(jc,jb) < zfrice_thrhld) THEN   ! non-seaice point
            ! now
            p_prog_wtr_now%h_ice    (jc,jb) = 0._wp
            p_prog_wtr_now%t_ice    (jc,jb) = tmelt   ! fresh water freezing point
            !
            ! new
            p_prog_wtr_new%h_ice    (jc,jb) = 0._wp
            p_prog_wtr_new%t_ice    (jc,jb) = tmelt   ! fresh water freezing point
            !
            IF (lprog_albsi) THEN
              p_prog_wtr_new%alb_si (jc,jb) = ALB_SI_MISSVAL
              p_prog_wtr_now%alb_si (jc,jb) = ALB_SI_MISSVAL
            ENDIF
          ENDIF
        ENDDO

      ELSE IF (.NOT. lseaice .AND.       llake) THEN

        ! Re-initialize all non-lake points
        !
        DO jc = i_startidx, i_endidx
          IF (.NOT. lake_mask(jc)) THEN   ! non-lake point
            ! now
            p_prog_wtr_now%h_ice    (jc,jb) = 0._wp
            p_prog_wtr_now%t_ice    (jc,jb) = tmelt   ! fresh water freezing point
            !
            ! new
            p_prog_wtr_new%h_ice    (jc,jb) = 0._wp
            p_prog_wtr_new%t_ice    (jc,jb) = tmelt   ! fresh water freezing point
          ENDIF
        ENDDO

      ELSE IF (      lseaice .AND.       llake) THEN

        ! Re-initialize all points which are neither lake nor sea-ice points
        !
        IF ( ntiles_total == 1 ) THEN  ! no tile approach
          zfrice_thrhld = 0.5_wp
        ELSE
          zfrice_thrhld = frsi_min
        ENDIF
        DO jc = i_startidx, i_endidx
          IF ( .NOT. lake_mask(jc) .AND. (p_lnd_diag%fr_seaice(jc,jb) < zfrice_thrhld) ) THEN
            ! now
            p_prog_wtr_now%h_ice    (jc,jb) = 0._wp
            p_prog_wtr_now%t_ice    (jc,jb) = tmelt   ! fresh water freezing point
            !
            ! new
            p_prog_wtr_new%h_ice    (jc,jb) = 0._wp
            p_prog_wtr_new%t_ice    (jc,jb) = tmelt   ! fresh water freezing point
            !
            IF (lprog_albsi) THEN
              p_prog_wtr_new%alb_si (jc,jb) = ALB_SI_MISSVAL
              p_prog_wtr_now%alb_si (jc,jb) = ALB_SI_MISSVAL
            ENDIF
          ENDIF
        ENDDO
      ENDIF



      IF(lsnowtile) THEN ! snow is considered as separate tiles
        DO isubs = 1, ntiles_lnd

          isubs_snow = isubs + ntiles_lnd

          CALL update_idx_lists_lnd (idx_lst_lp       = ext_data%atm%idx_lst_lp_t(:,jb,isubs),         &
                                   lp_count           = ext_data%atm%lp_count_t(jb,isubs),             &
                                   idx_lst            = ext_data%atm%idx_lst_t(:,jb,isubs),            &
                                   gp_count           = ext_data%atm%gp_count_t(jb,isubs),             &
                                   idx_lst_snow       = ext_data%atm%idx_lst_t(:,jb,isubs_snow),       &
                                   gp_count_snow      = ext_data%atm%gp_count_t(jb,isubs_snow),        &
                                   lc_frac            = ext_data%atm%lc_frac_t(:,jb,isubs),            &
                                   partial_frac       = ext_data%atm%frac_t(:,jb,isubs),               &
                                   partial_frac_snow  = ext_data%atm%frac_t(:,jb,isubs_snow),          &
                                   snowtile_flag      = ext_data%atm%snowtile_flag_t(:,jb,isubs),      &
                                   snowtile_flag_snow = ext_data%atm%snowtile_flag_t(:,jb,isubs_snow), &
                                   snowfrac           = p_lnd_diag%snowfrac_lc_t(:,jb,isubs)           )

          i_count = ext_data%atm%gp_count_t(jb,isubs)
          i_count_snow = ext_data%atm%gp_count_t(jb,isubs_snow)
!$NEC ivdep
          DO ic = 1, i_count_snow
            jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)

            ! reset field for actual snow-cover for grid points / land-cover classes for which there
            ! are seperate snow-free and snow-covered tiles
            p_lnd_diag%snowfrac_t(jc,jb,isubs)   = 0._wp
            p_prog_lnd_now%w_snow_t(jc,jb,isubs) = 0._wp
            p_lnd_diag%h_snow_t(jc,jb,isubs)     = 0._wp
            p_prog_lnd_now%t_snow_t(jc,jb,isubs) = p_prog_lnd_now%t_s_t(jc,jb,isubs)
            p_prog_lnd_now%t_g_t(jc,jb,isubs)    = p_prog_lnd_now%t_s_t(jc,jb,isubs)

            ! copy rho_snow in order to get the right tile average of snow density
            p_prog_lnd_now%rho_snow_t(jc,jb,isubs) = p_prog_lnd_now%rho_snow_t(jc,jb,isubs_snow)

            ! to prevent numerical stability problems, we require at least 1 cm of snow in order to
            ! have a snow-cover fraction of 1 on snow tiles (not critical for the single-layer
            ! snow scheme, but the multi-layer snow model becomes numerically unstable within a few
            ! time steps when associating traces of snow with a snow-cover fraction of 1)
            p_lnd_diag%snowfrac_t(jc,jb,isubs_snow) = MIN(1._wp,p_lnd_diag%h_snow_t(jc,jb,isubs_snow)*100._wp)

            ! Rediagnose t_g according to the modified snow-cover fraction
            p_prog_lnd_now%t_g_t(jc,jb,isubs_snow) =  &
              p_lnd_diag%snowfrac_t(jc,jb,isubs_snow) * p_prog_lnd_now%t_snow_t(jc,jb,isubs_snow) + &
              (1._wp-p_lnd_diag%snowfrac_t(jc,jb,isubs_snow))*p_prog_lnd_now%t_s_t(jc,jb,isubs_snow)
            ! Rediagnose qv_s_t because of the rediagnosed t_g_t
            p_lnd_diag%qv_s_t(jc,jb,isubs_snow)  =                                                  &
              &                    spec_humi(sat_pres_ice  (p_prog_lnd_now%t_g_t(jc,jb,isubs_snow)),&
              &                    p_diag%pres_sfc(jc,jb) )
             IF (lmulti_snow) THEN
               p_prog_lnd_now%t_snow_mult_t(jc,nlev_snow+1,jb,isubs) = p_prog_lnd_now%t_s_t(jc,jb,isubs)
             ENDIF

          END DO

          IF (lmulti_snow) THEN
            DO jk=1,nlev_snow
              DO ic = 1, i_count_snow
                jc = ext_data%atm%idx_lst_t(ic,jb,isubs_snow)
                p_prog_lnd_now%t_snow_mult_t(jc,jk,jb,isubs) = p_prog_lnd_now%t_s_t(jc,jb,isubs)
                p_prog_lnd_now%wliq_snow_t(jc,jk,jb,isubs) = 0._wp
                p_prog_lnd_now%wtot_snow_t(jc,jk,jb,isubs) = 0._wp
                p_prog_lnd_now%dzh_snow_t (jc,jk,jb,isubs) = 0._wp
              ENDDO
            ENDDO
          ENDIF

        END DO
      END IF


      ! Remove snow and w_i on non-existing grid points. This has no impact on the prognostic
      ! results but is needed in order to have meaningful data on the tile-based fields 
      DO isubs = 1, ntiles_total
        DO jc = i_startidx, i_endidx
          IF (ext_data%atm%frac_t(jc,jb,isubs) < 1.e-10_wp) THEN
            p_lnd_diag%h_snow_t(jc,jb,isubs)     = 0._wp
            p_prog_lnd_now%w_snow_t(jc,jb,isubs) = 0._wp
            p_lnd_diag%snowfrac_t(jc,jb,isubs)   = 0._wp
            p_prog_lnd_now%w_i_t(jc,jb,isubs)    = 0._wp
          ENDIF
          IF (ext_data%atm%lc_frac_t(jc,jb,isubs) < 1.e-10_wp) THEN
            p_lnd_diag%snowfrac_lc_t(jc,jb,isubs)   = 0._wp
          ENDIF
        ENDDO
      END DO


    ENDDO  ! jb loop
!$OMP END DO
!$OMP END PARALLEL

    ! Aggregate t_g and qv_s
    ! Loop over all points (land AND water points)
    ! Aggregation has been moved to the end of the subroutine (PR)
    !
    CALL aggregate_tg_qvs( p_patch = p_patch,                      & ! in
      &                    frac_t  = ext_data%atm%frac_t(:,:,:),   & ! in
      &                    t_g_t   = p_prog_lnd_now%t_g_t(:,:,:),  & ! in
      &                    qv_s_t  = p_lnd_diag%qv_s_t(:,:,:),     & ! in
      &                    t_g     = p_prog_lnd_now%t_g(:,:),      & ! inout
      &                    qv_s    = p_lnd_diag%qv_s(:,:)          ) ! inout

    p_prog_lnd_new%t_g(:,:)  = p_prog_lnd_now%t_g(:,:)


    ! limited area mode, only
    ! Soil and surface tempertures in the boundary zone of the LAM 
    ! domain are filled with meaningful values
    IF (l_limited_area .AND. jg==1) THEN
      CALL init_lamlatbc_phys(p_patch, p_prog_lnd_now, p_prog_lnd_new, p_lnd_diag)
    ENDIF

  END SUBROUTINE nwp_surface_init


!-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Agregate t_g and qv_s_t
  !!
  !! Segregated from nwp_surface_init, to avoid recoding again
  !-------------------------------------------------------------------------

  SUBROUTINE aggregate_tg_qvs( p_patch, frac_t, t_g_t, qv_s_t, t_g, qv_s, lacc )

    TYPE(t_patch), INTENT(IN)    :: p_patch
    REAL(wp)     , INTENT(IN)    :: frac_t(:,:,:)   !< tile-specific area fraction
    REAL(wp)     , INTENT(IN)    :: t_g_t(:,:,:)    !< surface temperature (tiled)
    REAL(wp)     , INTENT(IN)    :: qv_s_t(:,:,:)   !< specific humidity at the surface (tiled)
    REAL(wp)     , INTENT(INOUT) :: t_g(:,:)        !< surface temperature (aggregated)
    REAL(wp)     , INTENT(INOUT) :: qv_s(:,:)       !< specific humidity at the surface (aggregated)

    LOGICAL, OPTIONAL, INTENT(IN) :: lacc

    ! Local array bounds:

    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices

    ! Local :

    INTEGER :: jc,jb,isubs
    REAL(wp):: t_g_s(nproma), qv_s_s(nproma)

    LOGICAL :: lzacc
!-------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    !$ACC DATA CREATE(qv_s_s, t_g_s) IF(lzacc)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,isubs,i_startidx,i_endidx,t_g_s,qv_s_s)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

       IF (ntiles_total == 1) THEN
         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
         !$ACC LOOP GANG VECTOR
         DO jc = i_startidx, i_endidx
           t_g (jc,jb) = t_g_t (jc,jb,1)
           qv_s(jc,jb) = qv_s_t(jc,jb,1)
         ENDDO
         !$ACC END PARALLEL
       ELSE ! aggregate fields over tiles
         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
         !$ACC LOOP GANG VECTOR
         DO jc = 1, nproma
           t_g_s(jc) = 0._wp
           qv_s_s(jc) = 0._wp
         END DO
         !$ACC END PARALLEL

         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
         !$ACC LOOP SEQ
         DO isubs = 1,ntiles_total+ntiles_water
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = i_startidx, i_endidx
             t_g_s (jc) = t_g_s (jc) + frac_t(jc,jb,isubs) * t_g_t (jc,jb,isubs)**4
             qv_s_s(jc) = qv_s_s(jc) + frac_t(jc,jb,isubs) * qv_s_t(jc,jb,isubs)
           ENDDO
         ENDDO

         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO jc = i_startidx, i_endidx
           t_g (jc,jb) = SQRT(SQRT(t_g_s(jc)))
           qv_s(jc,jb) = qv_s_s(jc)
         ENDDO
         !$ACC END PARALLEL

       ENDIF  ! with or without tiles

    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE aggregate_tg_qvs


!-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !! Wrapper for computation of aggregated land variables
  !!
  SUBROUTINE aggr_landvars(p_patch, ext_data, p_lnd_state, lacc)

    TYPE(t_patch),             INTENT(IN)    :: p_patch(:)     !< grid/patch info
    TYPE(t_external_data),     INTENT(IN)    :: ext_data(:)    !< external data
    TYPE(t_lnd_state),         INTENT(INOUT) :: p_lnd_state(:) !< prog and diag lnd state
    LOGICAL,                   INTENT(IN)    :: lacc

    ! Local variables
    INTEGER :: jg ! loop indices

    CALL assert_acc_device_only("aggr_landvars", lacc)

    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%domain_is_owned) CYCLE
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      IF (  atm_phy_nwp_config(jg)%inwp_surface == 1 ) THEN
        CALL aggregate_landvars( p_patch(jg), ext_data(jg),                 &
             p_lnd_state(jg)%prog_lnd(nnow_rcf(jg)), p_lnd_state(jg)%diag_lnd, &
             lacc=lacc)
      ENDIF

    ENDDO ! jg-loop

  END SUBROUTINE aggr_landvars

!-------------------------------------------------------------------------

  !! Aggregation of tile-specific, instantaneous soil and surface fields
  !!
  !! Please note that this routine is only called at output time steps 
  !! and meteogramm timesteps. Therefore, it may only contain the aggregation 
  !! of purely diagnostic instantaneous fields. I.e. the computation of 
  !! accumulated fields, or the aggregation of fields that enter the prognostic 
  !! computations is not allowed, as this will lead to erroneous results. 
  !!
  SUBROUTINE aggregate_landvars( p_patch, ext_data, lnd_prog, lnd_diag, lacc )

    TYPE(t_patch),               INTENT(in)   :: p_patch       !< grid/patch info
    TYPE(t_external_data),       INTENT(in)   :: ext_data      !< external data
    TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog      !< prog vars for sfc
    TYPE(t_lnd_diag),            INTENT(inout):: lnd_diag      !< diag vars for sfc
    LOGICAL,                     INTENT(IN)   :: lacc

    ! Local scalars:
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< number of child domains
    INTEGER :: jc, jb, jk, isubs

    REAL(wp) :: tilefrac ! fractional area covered by tile
    REAL(wp) :: rho_snow_lim !< Snow density limited to [crhosmin_ml, crhosmax_ml]

    LOGICAL :: lmask(nproma)  ! mask array (TRUE for landpoint)
    INTEGER :: icount         ! index list length per block
    INTEGER :: ic
    INTEGER :: styp           ! soil type index

    IF (ltimer) CALL timer_start(timer_nh_diagnostics)

    !$ACC DATA CREATE(lmask) PRESENT(ext_data, lnd_prog, lnd_diag, dzsoil) IF(lacc)

    i_nchdom  = MAX(1,p_patch%n_childdom)

    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,isubs,jk,tilefrac,lmask,icount,styp) ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)


      IF (ntiles_total == 1) THEN  ! just copy prognostic variables from tile 1
                                   ! to diagnostic aggregated variable

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          lnd_diag%t_snow   (jc,jb) = lnd_prog%t_snow_t   (jc,jb,1)
          lnd_diag%t_s      (jc,jb) = lnd_prog%t_s_t      (jc,jb,1)
          lnd_diag%t_sk     (jc,jb) = lnd_prog%t_sk_t     (jc,jb,1)
          lnd_diag%w_snow   (jc,jb) = lnd_prog%w_snow_t   (jc,jb,1)
          lnd_diag%rho_snow (jc,jb) = lnd_prog%rho_snow_t (jc,jb,1)
          lnd_diag%w_i      (jc,jb) = lnd_prog%w_i_t      (jc,jb,1)
          lnd_diag%h_snow   (jc,jb) = lnd_diag%h_snow_t   (jc,jb,1)
          lnd_diag%freshsnow(jc,jb) = lnd_diag%freshsnow_t(jc,jb,1)
          lnd_diag%snowfrac (jc,jb) = lnd_diag%snowfrac_t (jc,jb,1)
          lnd_diag%rstom    (jc,jb) = lnd_diag%rstom_t    (jc,jb,1)
        ENDDO

        IF (itype_interception == 2) THEN
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx, i_endidx
            lnd_diag%w_p    (jc,jb) = lnd_prog%w_p_t      (jc,jb,1)
            lnd_diag%w_s    (jc,jb) = lnd_prog%w_s_t      (jc,jb,1)
          ENDDO
        ENDIF

        IF (itype_trvg == 3) THEN
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx, i_endidx
            lnd_diag%plantevap(jc,jb) = lnd_diag%plantevap_t(jc,jb,1)
          ENDDO
        ENDIF

        DO jk=1,nlev_soil
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx, i_endidx
            lnd_diag%t_so    (jc,jk+1,jb) = lnd_prog%t_so_t    (jc,jk+1,jb,1)
            lnd_diag%w_so    (jc,jk,  jb) = lnd_prog%w_so_t    (jc,jk,  jb,1)
            lnd_diag%w_so_ice(jc,jk,  jb) = lnd_prog%w_so_ice_t(jc,jk,  jb,1)
          ENDDO
        ENDDO

        IF (l2lay_rho_snow .OR. lmulti_snow) THEN
          !$ACC LOOP SEQ
          DO jk=1,nlev_snow
            !$ACC LOOP GANG VECTOR
            DO jc = i_startidx, i_endidx
              lnd_diag%rho_snow_mult(jc,jk,jb) = lnd_prog%rho_snow_mult_t(jc,jk,jb,1)
            ENDDO
          ENDDO
        ENDIF
        !$ACC END PARALLEL

        IF (lmulti_snow) THEN
#ifdef _OPENACC
          CALL finish('aggregate_landvars', 'lmulti_snow is not ported to openACC.')
#endif
          DO jk=1,nlev_snow
            DO jc = i_startidx, i_endidx
              lnd_diag%t_snow_mult  (jc,jk,jb) = lnd_prog%t_snow_mult_t(jc,jk,jb,1)
              lnd_diag%wliq_snow    (jc,jk,jb) = lnd_prog%wliq_snow_t(jc,jk,jb,1)
              lnd_diag%wtot_snow    (jc,jk,jb) = lnd_prog%wtot_snow_t(jc,jk,jb,1)
              lnd_diag%dzh_snow     (jc,jk,jb) = lnd_prog%dzh_snow_t(jc,jk,jb,1)
            ENDDO
          ENDDO
          DO jc = i_startidx, i_endidx
            lnd_diag%t_snow_mult(jc,nlev_snow+1,jb) = lnd_prog%t_snow_mult_t(jc,nlev_snow+1,jb,1)
          ENDDO
        ENDIF

      ELSE ! aggregate fields over tiles


        ! First initialize fields to zero in order to prepare
        ! subsequent summation over the tiles
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          lnd_diag%t_snow   (jc,jb)  = 0._wp
          lnd_diag%t_s      (jc,jb)  = 0._wp
          lnd_diag%t_sk     (jc,jb)  = 0._wp
          lnd_diag%w_snow   (jc,jb)  = 0._wp
          lnd_diag%w_i      (jc,jb)  = 0._wp
          lnd_diag%h_snow   (jc,jb)  = 0._wp
          lnd_diag%freshsnow(jc,jb)  = 0._wp
          lnd_diag%snowfrac (jc,jb)  = 0._wp
          lnd_diag%rstom    (jc,jb)  = 0._wp
        ENDDO
        

        !$ACC LOOP SEQ
        DO jk = 1, nlev_soil+1
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx, i_endidx
            lnd_diag%t_so    (jc,jk,jb) = 0._wp
          ENDDO
        ENDDO
        !$ACC LOOP SEQ
        DO jk = 1, nlev_soil
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx, i_endidx
            lnd_diag%w_so    (jc,jk,jb) = 0._wp
            lnd_diag%w_so_ice(jc,jk,jb) = 0._wp
          ENDDO
        ENDDO

        IF (itype_interception == 2) THEN
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx, i_endidx
            lnd_diag%w_p    (jc,jb)  = 0._wp
            lnd_diag%w_s    (jc,jb)  = 0._wp
          ENDDO
        ENDIF

        IF (itype_trvg == 3) THEN
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx, i_endidx
            lnd_diag%plantevap(jc,jb) = 0._wp
          ENDDO
        ENDIF

        IF (l2lay_rho_snow .OR. lmulti_snow) THEN
          !$ACC LOOP SEQ
          DO jk = 1, nlev_snow
            !$ACC LOOP GANG VECTOR
            DO jc = i_startidx, i_endidx
              lnd_diag%rho_snow_mult(jc,jk,jb) = 0._wp
            ENDDO
          ENDDO
        ENDIF
        !$ACC END PARALLEL

        IF (lmulti_snow) THEN
#ifdef _OPENACC
          CALL finish('aggregate_landvars', 'lmulti_snow is not ported to openACC.')
#endif        
          DO jk = 1, nlev_snow+1
            DO jc = i_startidx, i_endidx
              lnd_diag%t_snow_mult  (jc,jk,jb) = 0._wp
            ENDDO
          ENDDO
          DO jk = 1, nlev_snow
            DO jc = i_startidx, i_endidx
              lnd_diag%wliq_snow    (jc,jk,jb) = 0._wp
              lnd_diag%wtot_snow    (jc,jk,jb) = 0._wp
              lnd_diag%dzh_snow     (jc,jk,jb) = 0._wp
            ENDDO
          ENDDO
        ENDIF


        ! Aggregate fields without water tiles
        !
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
        !$ACC LOOP SEQ
        DO isubs = 1,ntiles_total
          !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(tilefrac)
          DO jc = i_startidx, i_endidx
            !
            ! note that frac_t must be re-scaled such that SUM(frac_t(1:ntiles_lnd)) = 1
            ! therefore we multiply by inv_frland_from_tiles
            tilefrac = ext_data%atm%frac_t(jc,jb,isubs)        &
              &      * ext_data%atm%inv_frland_from_tiles(jc,jb)

            lnd_diag%t_snow(jc,jb)    = lnd_diag%t_snow(jc,jb) + tilefrac    &
              &                         * lnd_prog%t_snow_t(jc,jb,isubs)
            lnd_diag%w_snow(jc,jb)    = lnd_diag%w_snow(jc,jb) + tilefrac    &
              &                         * lnd_prog%w_snow_t(jc,jb,isubs)
            lnd_diag%w_i(jc,jb)       = lnd_diag%w_i(jc,jb) + tilefrac       &
              &                         * lnd_prog%w_i_t(jc,jb,isubs)
            lnd_diag%h_snow(jc,jb)    = lnd_diag%h_snow(jc,jb) + tilefrac    &
              &                         * lnd_diag%h_snow_t(jc,jb,isubs)
            lnd_diag%freshsnow(jc,jb) = lnd_diag%freshsnow(jc,jb) + tilefrac &
              &                         * lnd_diag%freshsnow_t(jc,jb,isubs)
            lnd_diag%snowfrac(jc,jb)  = lnd_diag%snowfrac(jc,jb) + tilefrac  &
              &                         * lnd_diag%snowfrac_t(jc,jb,isubs)
            lnd_diag%rstom(jc,jb)     = lnd_diag%rstom(jc,jb) + tilefrac     &
              &                         * lnd_diag%rstom_t(jc,jb,isubs)
          ENDDO

          IF (itype_interception == 2) THEN
            !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(tilefrac)
            DO jc = i_startidx, i_endidx
              !
              ! note that frac_t must be re-scaled such that SUM(frac_t(1:ntiles_lnd)) = 1
              ! therefore we multiply by inv_frland_from_tiles
              tilefrac = ext_data%atm%frac_t(jc,jb,isubs)        &
                &      * ext_data%atm%inv_frland_from_tiles(jc,jb)

              lnd_diag%w_p    (jc,jb) = lnd_diag%w_p(jc,jb) + tilefrac       &
                &                       * lnd_prog%w_p_t(jc,jb,isubs)
              lnd_diag%w_s    (jc,jb) = lnd_diag%w_s(jc,jb) + tilefrac       &
                &                       * lnd_prog%w_s_t(jc,jb,isubs)
            ENDDO
          ENDIF

          IF (itype_trvg == 3) THEN
            !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(tilefrac)
            DO jc = i_startidx, i_endidx
              !
              ! note that frac_t must be re-scaled such that SUM(frac_t(1:ntiles_lnd)) = 1
              ! therefore we multiply by inv_frland_from_tiles
              tilefrac = ext_data%atm%frac_t(jc,jb,isubs)        &
                &      * ext_data%atm%inv_frland_from_tiles(jc,jb)

              lnd_diag%plantevap(jc,jb) = lnd_diag%plantevap(jc,jb) + tilefrac       &
                &                       * lnd_diag%plantevap_t(jc,jb,isubs)
            ENDDO
          ENDIF

          !$ACC LOOP SEQ
          DO jk=1,nlev_soil
            !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(tilefrac)
            DO jc = i_startidx, i_endidx
              !
              ! note that frac_t must be re-scaled such that SUM(frac_t(1:ntiles_lnd)) = 1
              ! therefore we multiply by inv_frland_from_tiles
              tilefrac = ext_data%atm%frac_t(jc,jb,isubs)         &
                &      * ext_data%atm%inv_frland_from_tiles(jc,jb)

              lnd_diag%t_so(jc,jk+1,jb)    = lnd_diag%t_so(jc,jk+1,jb) + tilefrac &
                                             * lnd_prog%t_so_t(jc,jk+1,jb,isubs)
              lnd_diag%w_so(jc,jk,jb)      = lnd_diag%w_so(jc,jk,jb) + tilefrac  &
                                             * lnd_prog%w_so_t(jc,jk,jb,isubs)
              lnd_diag%w_so_ice(jc,jk,jb)  = lnd_diag%w_so_ice(jc,jk,jb) + tilefrac  &
                                             * lnd_prog%w_so_ice_t(jc,jk,jb,isubs)
            ENDDO
          ENDDO

          IF (l2lay_rho_snow .OR. lmulti_snow) THEN
            !$ACC LOOP SEQ
            DO jk=1,nlev_snow
              !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(tilefrac)
              DO jc = i_startidx, i_endidx
                tilefrac = ext_data%atm%frac_t(jc,jb,isubs)        &
                  &      * ext_data%atm%inv_frland_from_tiles(jc,jb)

                lnd_diag%rho_snow_mult(jc,jk,jb) = lnd_diag%rho_snow_mult(jc,jk,jb) + tilefrac  &
                                                   * lnd_prog%rho_snow_mult_t(jc,jk,jb,isubs)
              ENDDO
            ENDDO
          ENDIF

        ENDDO  ! isubs
        !$ACC END PARALLEL


        IF (lmulti_snow) THEN
          DO isubs = 1,ntiles_total
            DO jk=1,nlev_snow
              DO jc = i_startidx, i_endidx
                !
                ! note that frac_t must be re-scaled such that SUM(frac_t(1:ntiles_lnd)) = 1
                ! therefore we multiply by inv_frland_from_tiles
                tilefrac = ext_data%atm%frac_t(jc,jb,isubs)        &
                  &      * ext_data%atm%inv_frland_from_tiles(jc,jb)

                lnd_diag%t_snow_mult(jc,jk,jb)   = lnd_diag%t_snow_mult(jc,jk,jb) + tilefrac  &
                                                   * lnd_prog%t_snow_mult_t(jc,jk,jb,isubs)
                lnd_diag%wliq_snow(jc,jk,jb)     = lnd_diag%wliq_snow(jc,jk,jb) + tilefrac  &
                                                   * lnd_prog%wliq_snow_t(jc,jk,jb,isubs)
                lnd_diag%wtot_snow(jc,jk,jb)     = lnd_diag%wtot_snow(jc,jk,jb) + tilefrac  &
                                                   * lnd_prog%wtot_snow_t(jc,jk,jb,isubs)
                lnd_diag%dzh_snow(jc,jk,jb)      = lnd_diag%dzh_snow(jc,jk,jb) + tilefrac   &
                                                   * lnd_prog%dzh_snow_t(jc,jk,jb,isubs)
              ENDDO
            ENDDO
            !
            DO jc = i_startidx, i_endidx
              tilefrac = ext_data%atm%frac_t(jc,jb,isubs)        &
                &      * ext_data%atm%inv_frland_from_tiles(jc,jb)

              lnd_diag%t_snow_mult(jc,nlev_snow+1,jb) = lnd_diag%t_snow_mult(jc,nlev_snow+1,jb)+ &
                                        tilefrac * lnd_prog%t_snow_mult_t(jc,nlev_snow+1,jb,isubs)
            ENDDO
          ENDDO  ! isubs
        ENDIF  ! lmulti_snow


        ! Aggregate fields with water tiles
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
        !$ACC LOOP SEQ
        DO isubs = 1,ntiles_total + ntiles_water
          !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(tilefrac)
          DO jc = i_startidx, i_endidx
            tilefrac = ext_data%atm%frac_t(jc,jb,isubs)

            lnd_diag%t_s(jc,jb)       = lnd_diag%t_s(jc,jb) + tilefrac       &
              &                         * lnd_prog%t_s_t(jc,jb,isubs)
            lnd_diag%t_sk(jc,jb)      = lnd_diag%t_sk(jc,jb) + tilefrac       &
              &                         * lnd_prog%t_sk_t(jc,jb,isubs)
          ENDDO
        ENDDO

        ! diagnose rho_snow from aggregated values of w_snow and h_snow; 
        ! by convention, snow density is zero in the absence of snow
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(rho_snow_lim)
        DO jc = i_startidx, i_endidx
          rho_snow_lim             = rhoh2o*(lnd_diag%w_snow(jc,jb)/MAX(dbl_eps,lnd_diag%h_snow(jc,jb)))
          rho_snow_lim             = MAX(MIN(rho_snow_lim,crhosmax_ml),crhosmin_ml)
          lnd_diag%rho_snow(jc,jb) = MERGE(rho_snow_lim,0._wp,lnd_diag%w_snow(jc,jb)>0._wp)
        ENDDO
        !$ACC END PARALLEL

      ENDIF  ! ntiles_total == 1

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        ! copy aggregated snowfrac field to snowfrac_lc
        lnd_diag%snowfrac_lc(jc,jb) = lnd_diag%snowfrac(jc,jb)

        ! fill t_so(1) with t_s
        lnd_diag%t_so(jc,1,jb) = lnd_diag%t_s(jc,jb)
      ENDDO

      ! Fill t_so(2:nlev_soil+1) with SST for non-land points (fr_land <= frlnd_thrhld) 
      ! Note that all points with fr_land > frlnd_thrhld are stored in the 
      ! land point index list list_land
      !
      ! create mask array
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        lmask(jc) = .FALSE.
      ENDDO
      !$ACC END PARALLEL
      !
      icount = ext_data%atm%list_land%ncount(jb)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(jc) ASYNC(1) &
      !$ACC   IF(lacc)
      DO ic = 1, icount
        jc = ext_data%atm%list_land%idx(ic,jb)
        lmask(jc) = .TRUE.
      ENDDO  ! ic
      !$ACC END PARALLEL LOOP
      !
      ! fill non-land points with SST
      ! It is assured that all mixed land/water points only contain pure land temperatures
      ! Note that the climatological layer is explicitly omitted.
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
      !$ACC LOOP SEQ
      DO jk = 2, nlev_soil
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          lnd_diag%t_so(jc,jk,jb) = MERGE(lnd_diag%t_so(jc,jk,jb),lnd_diag%t_s(jc,jb),lmask(jc))
        ENDDO  ! jc
      ENDDO  ! jk
      !
      ! climatological layer: Fill with T_CL over non-land points. 
      ! Land points are already filled with T_CL
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        lnd_diag%t_so(jc,nlev_soil+1,jb) = MERGE(lnd_diag%t_so(jc,nlev_soil+1,jb), &
          &                                      ext_data%atm%t_cl(jc,jb),         &
          &                                      lmask(jc))
      ENDDO  ! jc


      ! Make sure that aggregated w_so is always larger than air dryness point 
      ! at points where the soiltype allows infiltration of water.
      !$ACC LOOP SEQ
      DO jk=1,nlev_soil
        !$ACC LOOP GANG VECTOR PRIVATE(styp)
        DO jc = i_startidx, i_endidx
          styp = ext_data%atm%soiltyp(jc,jb)
          IF ( (styp>=3) .AND. (styp<=8)) THEN   ! 3:sand; 8:peat
            lnd_diag%w_so(jc,jk,jb) = MAX(lnd_diag%w_so(jc,jk,jb),dzsoil(jk)*cadp(styp))
          ELSE
            ! make sure that aggregated w_so=0, w_so_ice=0, where no infiltration is allowed
            lnd_diag%w_so(jc,jk,jb)     = 0._wp
            lnd_diag%w_so_ice(jc,jk,jb) = 0._wp
          ENDIF
        ENDDO
      ENDDO
      !$ACC END PARALLEL

    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

    IF (ltimer) CALL timer_stop(timer_nh_diagnostics)

  END SUBROUTINE aggregate_landvars



  !! Driver routine to (re-)initialize the snowtile index lists in the case of a restart
  !!
  SUBROUTINE init_snowtile_lists(p_patch, ext_data, p_lnd_diag)

    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch       !<grid/patch info.
    TYPE(t_external_data), INTENT(INOUT) :: ext_data
    TYPE(t_lnd_diag)     , INTENT(INOUT) :: p_lnd_diag

    ! Local array bounds:

    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< number of child domains

    ! Local scalars:

    INTEGER :: jb,isubs,i_count,isubs_snow


    i_nchdom  = MAX(1,p_patch%n_childdom)

    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,isubs,i_count,isubs_snow), SCHEDULE(guided)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      DO isubs = 1, ntiles_lnd
        isubs_snow = isubs + ntiles_lnd

        i_count = ext_data%atm%lp_count_t(jb,isubs)

        IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty

        CALL update_idx_lists_lnd (idx_lst_lp       = ext_data%atm%idx_lst_lp_t(:,jb,isubs),         &
                                 lp_count           = ext_data%atm%lp_count_t(jb,isubs),             &
                                 idx_lst            = ext_data%atm%idx_lst_t(:,jb,isubs),            &
                                 gp_count           = ext_data%atm%gp_count_t(jb,isubs),             &
                                 idx_lst_snow       = ext_data%atm%idx_lst_t(:,jb,isubs_snow),       &
                                 gp_count_snow      = ext_data%atm%gp_count_t(jb,isubs_snow),        &
                                 lc_frac            = ext_data%atm%lc_frac_t(:,jb,isubs),            &
                                 partial_frac       = ext_data%atm%frac_t(:,jb,isubs),               &
                                 partial_frac_snow  = ext_data%atm%frac_t(:,jb,isubs_snow),          &
                                 snowtile_flag      = ext_data%atm%snowtile_flag_t(:,jb,isubs),      &
                                 snowtile_flag_snow = ext_data%atm%snowtile_flag_t(:,jb,isubs_snow), &
                                 snowfrac           = p_lnd_diag%snowfrac_lc_t(:,jb,isubs)           )

      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE init_snowtile_lists



  !-------------------------------------------------------------------------
  !! Initialize sea ice and open water index lists for restart and
  !! non-restart runs. Sea ice and open water points are distinguished
  !! based on fr_seaice, which is provided by external sources (e.g. analysis, ocean model, etc.).
  !!
  !! Note that fr_seaice is potentially modified.
  !! For fr_seaice in ]0,frsi_min[, it is set to 0
  !! For fr_seaice in ]1-frsi_min,1[, it is set to 1.  
  !!
  SUBROUTINE init_sea_lists(p_patch, lseaice, fr_seaice, ext_data, opt_lverbose, lacc)

    TYPE(t_patch)        , INTENT(IN)              :: p_patch        !< grid/patch info.
    LOGICAL              , INTENT(IN)              :: lseaice        !< seaice model on/off
    REAL(wp)             , INTENT(INOUT)           :: fr_seaice(:,:) !< seaice fraction
    TYPE(t_external_data), INTENT(INOUT)           :: ext_data
    LOGICAL              , INTENT(IN), OPTIONAL    :: opt_lverbose   !< trigger message() output
    LOGICAL              , INTENT(IN), OPTIONAL    :: lacc

    ! Local array bounds:

    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk

    ! Local scalars:
    !
    INTEGER :: jb, ic, jc
    INTEGER :: jg
    INTEGER :: i_count_sea, i_count_ice, i_count_water, ncount_ice, ncount_water
    INTEGER :: idx_ice(nproma), idx_water(nproma), cond_ice(nproma), cond_water(nproma) ! used for generating index list
    INTEGER :: npoints_ice, npoints_wtr, npoints_sea
    REAL(wp):: frac_sea, lc_frac_t
    REAL(wp):: san_frac_sea, san_lc_frac_t, diff_frac, max_diff ! for sanity check
    LOGICAL :: lverbose
    LOGICAL :: lzacc

    CHARACTER(len=*), PARAMETER :: routine = 'mo_nwp_sfc_utils:init_sea_lists'
!-------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA CREATE(idx_ice, idx_water, cond_ice, cond_water) IF(lzacc)

    IF ( PRESENT(opt_lverbose) ) THEN
      lverbose = opt_lverbose
    ELSE
      lverbose = .TRUE.
    ENDIF

    ! patch ID
    jg = p_patch%id

    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)


    IF (lseaice) THEN

    ! generate sea-ice and open-water index list
    !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,ic,jc,i_count_sea,i_count_ice,i_count_water,frac_sea,ncount_ice,ncount_water,idx_ice,idx_water, &
!$OMP            cond_ice,cond_water,max_diff,lc_frac_t,diff_frac,san_frac_sea,san_lc_frac_t), SCHEDULE(guided)
    DO jb = i_startblk, i_endblk


      ! This (re)initialization is needed for correct functionality of the IAU iteration
      ext_data%atm%list_seawtr%ncount(jb) = 0
      ext_data%atm%list_seaice%ncount(jb) = 0

      ! Init sub-index lists for sea points. We distinguish between open-water
      ! (i.e. ice free) points and sea-ice points. diag_lnd%fr_seaice is used
      ! to distinguish between open-water and sea-ice.
      ! These index lists are needed by the sea-ice model
      !
      i_count_sea   = ext_data%atm%list_sea%ncount(jb)
      i_count_ice   = 0
      i_count_water = 0

      IF (i_count_sea == 0) CYCLE ! skip loop if the index list for the given block is empty

      ! For fr_seaice in ]0,frsi_min[, set fr_seaice to 0
      ! For fr_seaice in ]1-frsi_min,1[, set fr_seaice to 1. 
      ! This will ensure that sea-ice and water fractions sum up exactly 
      ! to the total sea fraction.
!$NEC ivdep
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jc)
      DO ic = 1, i_count_sea
        jc = ext_data%atm%list_sea%idx(ic,jb)
        IF (fr_seaice(jc,jb) < frsi_min ) THEN
           fr_seaice(jc,jb) = 0._wp
        ENDIF
        IF (fr_seaice(jc,jb) >= (1._wp-frsi_min) ) THEN
           fr_seaice(jc,jb) = 1._wp
        ENDIF
      ENDDO  ! ic
      !$ACC END PARALLEL


      IF ( ntiles_total == 1 ) THEN  ! no tile approach
#ifdef _OPENACC
        IF (lzacc) THEN
          CALL finish('init_sea_lists', "The code path without tiling is not ported to GPU")
        ENDIF
#endif _OPENACC

        !
        ! mixed water/ice points are not allowed. A sea point can be either
        ! ice-free or completely ice-covered
        !

!$NEC ivdep
        DO ic = 1, i_count_sea

          jc = ext_data%atm%list_sea%idx(ic,jb)


          ! ext_data%atm%lc_frac_t(jc,jb,1) already set to 1 in
          ! init_index_lists for sea points


          IF ( fr_seaice(jc,jb) >= 0.5_wp ) THEN
            !
            ! sea-ice point
            !
            i_count_ice = i_count_ice + 1
            ext_data%atm%list_seaice%idx(i_count_ice,jb) = jc
            ext_data%atm%list_seaice%ncount(jb)          = i_count_ice
            ! set surface area index (needed by turbtran)
            ext_data%atm%sai_t(jc,jb,isub_seaice) = c_sea
          ELSE
            !
            ! water point: all sea points with fr_seaice < 0.5
            !
            i_count_water = i_count_water + 1
            ext_data%atm%list_seawtr%idx(i_count_water,jb) = jc
            ext_data%atm%list_seawtr%ncount(jb)            = i_count_water
          ENDIF


          ! Initialize frac_t for sea-ice and water points
          ext_data%atm%frac_t(jc,jb,1) = ext_data%atm%lc_frac_t(jc,jb,1)

        ENDDO  ! ic

      ELSE   ! tile approach

!$NEC ivdep
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jc)
        DO ic = 1, i_count_sea
          jc = ext_data%atm%list_sea%idx(ic,jb)

          ! set sea-ice area fraction (static)
          ! simply copy from water tile. time-dependent fraction will be set lateron
          ext_data%atm%lc_frac_t(jc,jb,isub_seaice)= ext_data%atm%lc_frac_t(jc,jb,isub_water)

          !
          ! seaice point
          !
          IF ( fr_seaice(jc,jb) >= frsi_min ) THEN
            cond_ice(ic) = 1
            ! Initialize frac_t for seaice
            ext_data%atm%frac_t(jc,jb,isub_seaice) = ext_data%atm%lc_frac_t(jc,jb,isub_seaice) &
              &                                    * fr_seaice(jc,jb)

!DR Note that sai at seaice points is initialized with c/=c_sea, a corresponding update
!DR of sai_t needs to be added to the procedure which updates the seaice index list.
            ! set surface area index (needed by turbtran)
            ext_data%atm%sai_t(jc,jb,isub_seaice)  = c_sea
          ELSE
            cond_ice(ic) = 0
            ext_data%atm%frac_t(jc,jb,isub_seaice) = 0._wp
          ENDIF

          !
          ! water point: all sea points with fr_seaice < (1-frsi_min)
          !
          IF ( fr_seaice(jc,jb) <= (1._wp-frsi_min) ) THEN
            cond_water(ic) = 1
            ! Update frac_t for water tile
            ext_data%atm%frac_t(jc,jb,isub_water)  = ext_data%atm%lc_frac_t(jc,jb,isub_water)  &
              &                                    * (1._wp - fr_seaice(jc,jb))
          ELSE
            cond_water(ic) = 0
            ! necessary, since frac_t(jc,jb,isub_water) has already been initialized
            ! with nonzero values in init_index_lists
            ext_data%atm%frac_t(jc,jb,isub_water)  = 0._wp
          ENDIF
        ENDDO
        !$ACC END PARALLEL

        CALL generate_index_list(cond_ice, idx_ice, 1, i_count_sea, ext_data%atm%list_seaice%ncount(jb), &
          1, opt_acc_copy_to_host=.FALSE., opt_use_acc=lzacc)
        CALL generate_index_list(cond_water, idx_water, 1, i_count_sea, ext_data%atm%list_seawtr%ncount(jb), &
          1, opt_acc_copy_to_host=.FALSE., opt_use_acc=lzacc)
        !$ACC WAIT
        !$ACC UPDATE HOST(ext_data%atm%list_seaice%ncount(jb:jb), ext_data%atm%list_seawtr%ncount(jb:jb)) IF(lzacc)
        ncount_ice = ext_data%atm%list_seaice%ncount(jb)
        ncount_water = ext_data%atm%list_seawtr%ncount(jb)

        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO ic = 1, ncount_ice
          ext_data%atm%list_seaice%idx(ic,jb) = ext_data%atm%list_sea%idx(idx_ice(ic),jb)
        ENDDO
        !$ACC END PARALLEL

        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO ic = 1, ncount_water
          ext_data%atm%list_seawtr%idx(ic,jb) = ext_data%atm%list_sea%idx(idx_water(ic),jb)
        ENDDO
        !$ACC END PARALLEL


#ifndef __SX__
        ! Sanity check
        ! Check whether fractions of seaice and non-seaice covered tiles sum up to total sea fraction. 
        max_diff = 0.0_wp
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) REDUCTION(MAX: max_diff) IF(lzacc)
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jc, frac_sea, lc_frac_t, diff_frac) REDUCTION(MAX: max_diff)
        DO ic = 1,i_count_sea
          jc = ext_data%atm%list_sea%idx(ic,jb)
          frac_sea = ext_data%atm%frac_t(jc,jb,isub_water) + ext_data%atm%frac_t(jc,jb,isub_seaice)
          lc_frac_t = ext_data%atm%lc_frac_t(jc,jb,isub_water)
          diff_frac = ABS(frac_sea - lc_frac_t)
          max_diff = MAX(max_diff, diff_frac)
        ENDDO
        !$ACC END PARALLEL
        !$ACC WAIT

        IF (max_diff > dbl_eps) THEN ! Sanity check failed
          !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
          !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jc, frac_sea, lc_frac_t, diff_frac, san_frac_sea, san_lc_frac_t)
          DO ic = 1,i_count_sea
            jc = ext_data%atm%list_sea%idx(ic,jb)
            frac_sea = ext_data%atm%frac_t(jc,jb,isub_water) + ext_data%atm%frac_t(jc,jb,isub_seaice)
            lc_frac_t = ext_data%atm%lc_frac_t(jc,jb,isub_water)
            diff_frac =  ABS(frac_sea - lc_frac_t)
            IF (diff_frac == max_diff) THEN
              san_frac_sea = frac_sea
              san_lc_frac_t = lc_frac_t
              EXIT
            ENDIF
          ENDDO
          !$ACC END PARALLEL
          WRITE(message_text,'(a,f12.9)') 'frac_seaice + frac_water: ', san_frac_sea
          CALL message('', TRIM(message_text))
          WRITE(message_text,'(a,f12.9)') 'tot frac_sea: ', san_lc_frac_t
          CALL message('', TRIM(message_text))
          CALL finish(routine, 'sea-ice + water fractions do not sum up to total sea fraction')
        ENDIF
#endif
      ENDIF   ! IF (ntiles_total == 1)

    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL


      IF ( lverbose .OR. msg_level >= 13 ) THEN
        ! Some diagnostics: compute total number of sea-ice and open water points
        !$ACC WAIT
        !$ACC UPDATE HOST(ext_data%atm%list_seaice%idx, ext_data%atm%list_seaice%ncount) &
        !$ACC   HOST(ext_data%atm%list_seawtr%idx, ext_data%atm%list_seawtr%ncount) &
        !$ACC   HOST(ext_data%atm%list_sea%idx, ext_data%atm%list_sea%ncount) IF(lzacc)
        npoints_ice = ext_data%atm%list_seaice%get_sum_global(i_startblk,i_endblk)
        npoints_wtr = ext_data%atm%list_seawtr%get_sum_global(i_startblk,i_endblk)
        npoints_sea = ext_data%atm%list_sea   %get_sum_global(i_startblk,i_endblk)
        !
        WRITE(message_text,'(a,i3,a,i10)') 'Number of (partly) seaice covered points in domain',jg, &
          &  ':',npoints_ice
        CALL message('', TRIM(message_text))
        WRITE(message_text,'(a,i3,a,i10)') 'Number of (partly) ice-free sea points in domain',jg, &
          &  ':',npoints_wtr
        CALL message('', TRIM(message_text))
        WRITE(message_text,'(a,i3,a,i10)') 'Number of sea points in domain',jg, &
             &  ':',npoints_sea
        CALL message('', TRIM(message_text))
      ENDIF

    ELSE   ! seaice model switched off

      ! copy sea points index list to open-water index list
      !
      ext_data%atm%list_seawtr%ncount(i_startblk:i_endblk)     =     &
        &                     ext_data%atm%list_sea%ncount(i_startblk:i_endblk)
      ext_data%atm%list_seawtr%idx(:,i_startblk:i_endblk) =     &
        &                     ext_data%atm%list_sea%idx(:,i_startblk:i_endblk)

      ext_data%atm%list_seaice%ncount(i_startblk:i_endblk) = 0

    ENDIF  ! lseaice

    !$ACC END DATA ! idx_ice, idx_water, cond_ice, cond_water

  END SUBROUTINE init_sea_lists


  !-------------------------------------------------------------------------

  SUBROUTINE diag_snowfrac_tg(istart, iend, lc_class, i_lc_urban, t_snow, t_soiltop, w_snow,  &
    & rho_snow, freshsnow, sso_sigma, z0, snowfrac_fac, snowfrac, t_g, meltrate,              &
    & snowfrac_u, lacc, opt_acc_async_queue)

    INTEGER, INTENT (IN) :: istart, iend ! start and end-indices of the computation

    INTEGER, INTENT (IN) :: lc_class(:)  ! list of land-cover classes
    INTEGER, INTENT (IN) :: i_lc_urban   ! land-cover class index for urban / artificial surface
    REAL(wp), DIMENSION(:), INTENT(IN) :: t_snow, t_soiltop, w_snow, rho_snow, &
      freshsnow, sso_sigma, z0, snowfrac_fac
    REAL(wp), DIMENSION(:), INTENT(IN), OPTIONAL :: meltrate ! snow melting rate in kg/(m**2*s)

    REAL(wp), DIMENSION(:), INTENT(INOUT) :: snowfrac, t_g
    REAL(wp), DIMENSION(:), INTENT(INOUT), OPTIONAL :: snowfrac_u ! unmodified snow-cover fraction

    INTEGER  :: ic
    REAL(wp) :: h_snow, snowdepth_fac, sso_fac, lc_fac, lc_limit

    LOGICAL, OPTIONAL, INTENT(in)   :: lacc          !< GPU flag
    INTEGER, OPTIONAL, INTENT(in)   :: opt_acc_async_queue

    LOGICAL :: lzacc
    INTEGER :: acc_async_queue

    CALL set_acc_host_or_device(lzacc, lacc)

    IF(PRESENT(opt_acc_async_queue)) THEN
        acc_async_queue = opt_acc_async_queue
    ELSE
        acc_async_queue = 1
    ENDIF

    ! iend has to be in a separate data section as of NV HPC 23.1
    !$ACC DATA COPYIN(iend) ASYNC(acc_async_queue) IF(lzacc)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
    SELECT CASE (idiag_snowfrac)
    CASE (1) ! old parameterization depending on SWE only
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = istart, iend
        snowfrac(ic) = MIN(1.0_wp, w_snow(ic)/cf_snow)
        t_g(ic) = t_snow(ic) + (1.0_wp - snowfrac(ic))*(t_soiltop(ic) - t_snow(ic))
      ENDDO
    CASE (2, 20) ! more advanced parameterization depending on snow depth, accounts also for vegetation and SSO
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(h_snow, sso_fac, snowdepth_fac, lc_fac, lc_limit)
      DO ic = istart, iend
        IF (w_snow(ic) <= 1.e-6_wp) THEN
          snowfrac(ic) = 0._wp
        ELSE
          h_snow = 1000._wp*w_snow(ic)/rho_snow(ic)  ! snow depth in m
          sso_fac = SQRT(0.025_wp*MAX(25._wp,sso_sigma(ic)*(1._wp-freshsnow(ic))))
          snowdepth_fac = h_snow*(17.5_wp*freshsnow(ic)+5._wp+5._wp/sso_fac*(1._wp-freshsnow(ic)))
          lc_fac   = MAX(1._wp,SQRT(10._wp*z0(ic)))
          IF (lc_class(ic) == i_lc_urban) THEN
            lc_limit = 0.875_wp ! this accounts for the effect of human activities on snow cover
          ELSE
            lc_limit = 1._wp
          ENDIF
          snowfrac(ic) = MAX(0.05_wp,MIN(lc_limit,snowfrac_fac(ic)*snowdepth_fac/lc_fac))
        ENDIF
        t_g(ic) = t_snow(ic) + (1.0_wp - snowfrac(ic))*(t_soiltop(ic) - t_snow(ic))
      ENDDO
    END SELECT

    ! For the single-layer scheme, t_soiltop represents the weighted average between the snow-covered
    ! part and the snow-free part in the case of melting snow concurrent with above-freezing soil temperatures
    IF (.NOT. lmulti_snow) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = istart, iend
        IF (t_soiltop(ic) > tmelt .AND. t_snow(ic) >= tmelt - dbl_eps) t_g(ic) = t_soiltop(ic)
      ENDDO
    ENDIF
 
    IF (PRESENT(snowfrac_u)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic = istart, iend
       snowfrac_u(ic) = snowfrac(ic) ! save unmodified snow-cover fraction
      ENDDO
    ENDIF

    SELECT CASE (idiag_snowfrac)
    CASE (20)
      ! Artificially reduce snow-cover fraction in case of melting snow in order to reduce the ubiquitous
      ! cold bias in such situations
      IF (PRESENT(meltrate)) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO ic = istart, iend
          snowfrac(ic) = MIN(snowfrac(ic),MAX(tune_minsnowfrac,1._wp-20000._wp*MAX(0._wp,z0(ic)-0.02_wp)*meltrate(ic)))
        ENDDO
      ENDIF
    END SELECT
    !$ACC END PARALLEL
    IF (.NOT. lcuda_graph_lnd) THEN
      !$ACC WAIT(acc_async_queue)
    END IF
    !$ACC END DATA

  END SUBROUTINE diag_snowfrac_tg




  !-------------------------------------------------------------------------
  !
  !! Updating of dynamic index lists
  !!
  !! Routine updates the following dynamic index lists (if required):
  !!
  !! - dynamic
  !!
  SUBROUTINE update_idx_lists_lnd (idx_lst_lp, lp_count, idx_lst, gp_count, idx_lst_snow, &
    &             gp_count_snow, lc_frac, partial_frac, partial_frac_snow, snowtile_flag, &
    &             snowtile_flag_snow, snowfrac, lacc, opt_acc_async_queue)


    INTEGER ,    INTENT(   IN) ::  &   !< static list of all land points of a tile index
      &  idx_lst_lp(:), lp_count       !< and corresponding grid point counts


    INTEGER ,    INTENT(INOUT) ::  &   !< dynamic list of all snow-free or mixed land points
      &  idx_lst(:), gp_count          !< of a tile index
                                       !< and corresponding grid point counts


    INTEGER ,    INTENT(INOUT) ::  &   !< dynamic list of all snow-covered points for land-cover
      &  idx_lst_snow(:), gp_count_snow !< classes eligible for separate treatment
                                       !< and corresponding grid point counts

    INTEGER ,    INTENT(INOUT) ::  &   !< snowtile flag field for snow-free or mixed points
      &  snowtile_flag(:)              !< 0: inactive
                                       !< 1: active
                                       !< 2: newly activated; initialization from corresponding tile required

    INTEGER ,    INTENT(INOUT) ::  &   !< snowtile flag field for snow-covered points
      &  snowtile_flag_snow(:)         !< 0: inactive
                                       !< 1: active
                                       !< 2: newly activated; initialization from corresponding tile required

    REAL(wp),    INTENT(   IN) ::  &   !< area fraction of land-cover class
      &  lc_frac(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< snow-cover fraction
      &  snowfrac(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< snow-free and snow-covered sub-tile area fraction
      &  partial_frac(:), partial_frac_snow(:)


    ! Local variables
    INTEGER  :: ic, jc
    REAL(wp) :: eps = 1.e-6_wp

    LOGICAL, OPTIONAL, INTENT(in)   :: lacc          !< GPU flag
    INTEGER, OPTIONAL, INTENT(in)   :: opt_acc_async_queue

    LOGICAL :: lzacc
    INTEGER :: acc_async_queue

    INTEGER :: cond1(lp_count), cond2(lp_count)

    !-------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)
    
    IF(PRESENT(opt_acc_async_queue)) THEN
        acc_async_queue = opt_acc_async_queue
    ELSE
        acc_async_queue = 1
    ENDIF

    !$ACC DATA &
    !$ACC   PRESENT(gp_count, gp_count_snow) &
    !$ACC   PRESENT(idx_lst_lp, idx_lst, idx_lst_snow, snowtile_flag) &
    !$ACC   PRESENT(snowtile_flag_snow, lc_frac, snowfrac, partial_frac) &
    !$ACC   PRESENT(partial_frac_snow) &
    !$ACC   CREATE(cond1, cond2) ASYNC(acc_async_queue) IF(lzacc)

!$NEC ivdep
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
    !$ACC LOOP GANG VECTOR PRIVATE(jc)
    DO ic = 1, lp_count
      jc = idx_lst_lp(ic)

      ! Reset snowfrac to 0/1 in case of very small deviations (just to be safe)
      IF (snowfrac(jc) < eps) snowfrac(jc) = 0._wp
      IF (1._wp - snowfrac(jc) < eps) snowfrac(jc) = 1._wp

      IF (snowfrac(jc) > 0._wp) THEN ! snow tile is active
        cond1(ic) = 1
        partial_frac_snow(jc) = lc_frac(jc)*snowfrac(jc)
        IF (snowtile_flag_snow(jc) == 0) THEN
          snowtile_flag_snow(jc) = 2 ! newly activated, initialization needed
        ELSE
          snowtile_flag_snow(jc) = 1
        ENDIF
      ELSE
        cond1(ic) = 0
        snowtile_flag_snow(jc) = 0
        partial_frac_snow(jc)  = 0._wp
      ENDIF

      IF (snowfrac(jc) < 1._wp) THEN ! snow-free tile is active
        cond2(ic) = 1
        partial_frac(jc) = lc_frac(jc)*(1._wp-snowfrac(jc))
        IF (snowtile_flag(jc) == 0) THEN
          snowtile_flag(jc) = 2 ! newly activated, initialization needed
        ELSE
          snowtile_flag(jc) = 1
        ENDIF
      ELSE
        cond2(ic) = 0
        snowtile_flag(jc) = 0
        partial_frac(jc)  = 0._wp
      ENDIF
    ENDDO
    !$ACC END PARALLEL

    CALL generate_index_list(cond1, idx_lst_snow, 1, lp_count, gp_count_snow, &
      opt_acc_async_queue=acc_async_queue, opt_acc_copy_to_host=.FALSE., opt_use_acc=lzacc)

    CALL generate_index_list(cond2, idx_lst,      1, lp_count, gp_count, &
      opt_acc_async_queue=acc_async_queue, opt_acc_copy_to_host=.FALSE., opt_use_acc=lzacc)

    !$ACC PARALLEL LOOP GANG VECTOR ASYNC(acc_async_queue) IF(lzacc)
    DO ic = 1,gp_count_snow
      idx_lst_snow(ic) = idx_lst_lp(idx_lst_snow(ic))
    ENDDO

    !$ACC PARALLEL LOOP GANG VECTOR ASYNC(acc_async_queue) IF(lzacc)
    DO ic = 1,gp_count
      idx_lst(ic) = idx_lst_lp(idx_lst(ic))
    ENDDO

    IF (.NOT. lcuda_graph_lnd) THEN
      !$ACC WAIT(1)
    END IF
    !$ACC END DATA

  END SUBROUTINE update_idx_lists_lnd



  !-------------------------------------------------------------------------
  !
  !! Updating of dynamic index lists for sea points
  !!
  !! Routine updates the following dynamic index lists (if required):
  !!
  !! - dynamic open-water and sea-ice point index list, which are sub-index lists
  !!   of the static sea point index list list_sea
  !!
  SUBROUTINE update_idx_lists_sea (hice_n, pres_sfc, list_seawtr_idx, list_seawtr_count,&
    &                              list_seaice_idx, list_seaice_count, frac_t_ice,      &
    &                              frac_t_water, lc_frac_t_water, fr_seaice,            &
    &                              hice_old, tice_old, albsi_now, albsi_new,            &
    &                              t_g_t_now, t_g_t_new, t_s_t_now, t_s_t_new,          &
    &                              t_sk_t_now, t_sk_t_new, qv_s_t, t_seasfc, condhf,    &
    &                              meltpot                                              )


    REAL(wp),    INTENT(IN)    ::  &   !< sea-ice depth at new time level  [m]
      &  hice_n(:)                     !< dim: (nproma)

    REAL(wp),    INTENT(IN)    ::  &   !< surface pressure                 [Pa]
      &  pres_sfc(:)

    INTEGER ,    INTENT(INOUT) ::  &   !< dynamic open water point index list
      &  list_seawtr_idx(:),       &   !< and corresponding grid point counts
      &  list_seawtr_count

    INTEGER ,    INTENT(INOUT) ::  &   !< dynamic sea-ice point index list
      &  list_seaice_idx(:),       &   !< and corresponding grid point counts
      &  list_seaice_count

    REAL(wp),    INTENT(INOUT) ::  &   !< ice-covered and ice-free ocean area fraction
      &  frac_t_ice(:), frac_t_water(:)

    REAL(wp), INTENT(IN) ::        &   !< total ocean area fraction
      &  lc_frac_t_water(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< sea-ice fraction
      &  fr_seaice(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< sea-ice depth at old time level  [m]
      &  hice_old(:)                   !< dim: (nproma)

    REAL(wp),    INTENT(INOUT) ::  &   !< sea-ice temperature at old time level  [K]
      &  tice_old(:)                   !< dim: (nproma)

    REAL(wp),    INTENT(INOUT) ::  &   !< prognostic sea ice albedo at old time level  [1]
      &  albsi_now(:)                  !< dim: (nproma)

    REAL(wp),    INTENT(INOUT) ::  &   !< prognostic sea ice albedo at new time level  [1]
      &  albsi_new(:)                  !< dim: (nproma)

    REAL(wp),    INTENT(INOUT) ::  &   !< temperature of water tile (now)  [K]
      &  t_g_t_now(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< temperature of water tile (new)  [K]
      &  t_g_t_new(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< surface temperature of water tile (now)  [K]
      &  t_s_t_now(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< surface temperature of water tile (new)  [K]
      &  t_s_t_new(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< skin temperature of water tile (now)  [K]
      &  t_sk_t_now(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< skin temperature of water tile (new)  [K]
      &  t_sk_t_new(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< surface specific humidity        [kg/kg]
      &  qv_s_t(:)

    REAL(wp),    INTENT(INOUT) ::  &   !< sea surface temperature          [kg/kg]
      &  t_seasfc(:)

    REAL(wp), OPTIONAL, INTENT(INOUT) :: &
      &  condhf(:)                     !< conductive heat flux at bottom of sea-ice [W/m^2]

    REAL(wp), OPTIONAL, INTENT(INOUT) :: &
      &  meltpot(:)                    !< melt potential at top of sea-ice [W/m^2]

    ! Local variables
    INTEGER  :: list_seaice_count_old      !< old seaice index list and count
    INTEGER  :: list_seaice_idx_old(SIZE(list_seaice_idx,1))
    INTEGER  :: ic, jc                     !< loop indices
    INTEGER  :: i_capture !< to capture thread-local value in ACC ATOMIC
    LOGICAL  :: l_update_required
    LOGICAL  :: lhave_meltpot
    LOGICAL  :: lhave_condhf
    !-------------------------------------------------------------------------


    ! Test if a list update is required at all.
    ! The melting of a seaice tile is a relatively rare event.
    l_update_required = .FALSE.
!$NEC ivdep
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) REDUCTION(.OR.: l_update_required) PRIVATE(jc)
    DO ic = 1, list_seaice_count
      jc = list_seaice_idx(ic)
      IF ( hice_n(jc) < hice_min ) l_update_required = .TRUE.
    ENDDO
    !$ACC END PARALLEL LOOP
    IF (.NOT. lcuda_graph_lnd) THEN
      !$ACC WAIT(1)
    END IF
    IF (.NOT. l_update_required) RETURN

    IF (msg_level >= 13) CALL message('update_idx_lists_sea', &
      'One or more seaice cells melted -> List update required.')

#ifdef __NVCOMPILER
    ! nvfortran does not understand passing a NULL pointer to an optional (Fortran 2008) :(
    lhave_meltpot = is_coupled_to_ocean()
    lhave_condhf = is_coupled_to_ocean()
#else
    lhave_meltpot = PRESENT(meltpot)
    lhave_condhf = PRESENT(condhf)
#endif

    !$ACC DATA PRESENT(condhf) IF(lhave_condhf)
    !$ACC DATA PRESENT(meltpot) IF(lhave_meltpot)
    !$ACC DATA CREATE(list_seaice_idx_old) &
    !$ACC   PRESENT(hice_n, pres_sfc, list_seawtr_idx) &
    !$ACC   PRESENT(list_seaice_idx, frac_t_ice) &
    !$ACC   PRESENT(frac_t_water, lc_frac_t_water, fr_seaice) &
    !$ACC   PRESENT(hice_old, tice_old, albsi_now, albsi_new) &
    !$ACC   PRESENT(t_g_t_now, t_g_t_new, t_s_t_now, t_s_t_new) &
    !$ACC   PRESENT(t_sk_t_now, t_sk_t_new, qv_s_t, t_seasfc) NO_CREATE(condhf, meltpot)

    !$ACC PARALLEL LOOP GANG VECTOR ASYNC(1) DEFAULT(PRESENT)
    DO ic = 1, list_seaice_count
      ! save old sea-ice index list and grid point count
      list_seaice_idx_old(ic) = list_seaice_idx(ic)
      ! re-initialize sea-ice index list and grid point count
      list_seaice_idx(ic) = 0
    END DO
    !$ACC END PARALLEL LOOP
    list_seaice_count_old  = list_seaice_count

    !
    ! update index list for sea-ice and open water
    !
    ! The current sea-ice model does not allow for new sea-ice points to be
    ! created during model integration. However, the number of sea-ice points 
    ! is allowed to decrease with time, i.e. sea-ice points may turn  
    ! into water points, but not vice versa.
    !
    ! Loop over old sea-ice points, only

    IF ( ntiles_total == 1 ) THEN  ! no tile approach
#ifdef _OPENACC
      CALL finish('update_idx_lists_sea', "The code path without tiling is not tested on GPU")
#endif
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) &
      !$ACC   PRESENT(list_seawtr_count, list_seaice_count) ! these are entries of vectors that are present on device
      list_seaice_count = 0 ! do the reset on accelerator if using OpenACC
!$NEC ivdep
      !$ACC LOOP GANG VECTOR PRIVATE(jc, i_capture)
      DO ic = 1, list_seaice_count_old
        jc = list_seaice_idx_old(ic)

        IF ( hice_n(jc) >= hice_min )  THEN ! still sea-ice point
          !$ACC ATOMIC CAPTURE
          list_seaice_count = list_seaice_count + 1
          i_capture = list_seaice_count
          !$ACC END ATOMIC
          list_seaice_idx(i_capture) = jc
          ! sea-ice fraction remains unchanged
        ELSE                        ! sea-ice point has turned into water point
          ! initialize new water tile
          !$ACC ATOMIC CAPTURE
          list_seawtr_count = list_seawtr_count + 1
          i_capture = list_seawtr_count
          !$ACC END ATOMIC
          list_seawtr_idx(i_capture) = jc
          ! Initialize temperature of water tile with salt water freezing point
          t_g_t_new(jc) = tf_salt ! if the SST analysis contains a meaningful water
                                  ! temperature for this point, one may also take
                                  ! the latter
          t_g_t_now(jc) = tf_salt

          t_seasfc(jc)  = tf_salt

          ! Initialize surface saturation specific humidity for new water tile
          ! includes reduction of saturation pressure due to salt content
          qv_s_t(jc) = salinity_fac * spec_humi(sat_pres_water(t_g_t_new(jc)),pres_sfc(jc))

          ! since sea-ice melted away, the sea-ice fraction is re-set to 0
          fr_seaice(jc)  = 0._wp
          !
          ! resetting of frac_t_water and frac_t_ice is not possible without
          ! tile approach, since both point to the same array location (frac_t:,:,1).
          ! frac_t=1 is required without tile approach

          ! reset sea-ice temperature and depth at old time level in order to prevent
          ! other schemes from using them incorrectly
          tice_old(jc) = tmelt
          hice_old(jc) = 0._wp

            ! also reset ice heat fluxes
          IF (lhave_condhf) condhf(jc) = 0._wp
          IF (lhave_meltpot) meltpot(jc) = 0._wp

          ! Reset prognostic sea ice albedo for consistency
          IF (lprog_albsi) THEN
            albsi_now(jc) = ALB_SI_MISSVAL
            albsi_new(jc) = ALB_SI_MISSVAL
          ENDIF
        ENDIF
      ENDDO  ! ic

      !$ACC END PARALLEL
    ELSE
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) &
      !$ACC   PRESENT(list_seawtr_count, list_seaice_count) ! these are entries of vectors that are present on device
      list_seaice_count = 0 ! do the reset on accelerator if using OpenACC
!$NEC ivdep
      !$ACC LOOP GANG VECTOR PRIVATE(jc, i_capture)
      DO ic = 1, list_seaice_count_old
        jc = list_seaice_idx_old(ic)

        IF ( hice_n(jc) >= hice_min )  THEN ! still sea-ice point
          !$ACC ATOMIC CAPTURE
          list_seaice_count = list_seaice_count + 1
          i_capture = list_seaice_count
          !$ACC END ATOMIC
          list_seaice_idx(i_capture) = jc
          ! sea-ice fraction remains unchanged
        ELSE                        ! sea-ice point has turned into water point
          ! Check whether we need to initialize a new water tile, or whether a water tile
          ! already exists for the given point:
          IF ( fr_seaice(jc) > (1._wp-frsi_min) ) THEN
            ! water tile does not exist for given point
            ! add new water tile to water-points index list and initialize
            !$ACC ATOMIC CAPTURE
            list_seawtr_count = list_seawtr_count + 1
            i_capture = list_seawtr_count
            !$ACC END ATOMIC
            list_seawtr_idx(i_capture) = jc

            ! Initialize new water tile
            !
            ! Initialize temperature with salt water freezing point
            t_g_t_new(jc) = tf_salt ! if the SST analysis contains a meaningful water
                                    ! temperature for this point, one may also take
                                    ! the latter
            t_g_t_now(jc) = tf_salt

            t_s_t_new(jc) = tf_salt ! otherwise aggregated t_so and t_s will be
                                    ! 0 at these points
            t_s_t_now(jc) = tf_salt

            t_sk_t_new(jc) = tf_salt
            t_sk_t_now(jc) = tf_salt

            t_seasfc(jc)  = tf_salt

            !
            ! Initialize surface saturation specific humidity
            ! includes reduction of saturation pressure due to salt content
            qv_s_t(jc) = salinity_fac * spec_humi(sat_pres_water(t_g_t_new(jc)),pres_sfc(jc))
          ENDIF

          ! re-set dynamic fractions of water and sea-ice
          !
          ! new sea area fraction is the sum of the current water and sea-ice area fractions
          ! to ensure bit-reproducibility with restart, we use lc_frac_t instead of frac_t_water+frac_t_ice
          frac_t_water(jc) = lc_frac_t_water(jc)
          ! since sea-ice melted away, the sea-ice fraction is re-set to 0
          fr_seaice(jc)  = 0._wp
          frac_t_ice(jc) = 0._wp
          !
          ! reset sea-ice temperature and depth at old time level in order to prevent
          ! other schemes from using them incorrectly
          tice_old(jc) = tmelt
          hice_old(jc) = 0._wp

            ! also reset ice heat fluxes
          IF (lhave_condhf) condhf(jc) = 0._wp
          IF (lhave_meltpot) meltpot(jc) = 0._wp

          ! Reset prognostic sea ice albedo for consistency
          IF (lprog_albsi) THEN
            albsi_now(jc) = ALB_SI_MISSVAL
            albsi_new(jc) = ALB_SI_MISSVAL
          ENDIF
        ENDIF
      ENDDO  ! ic
      !$ACC END PARALLEL

    ENDIF  ! IF ( ntiles_total == 1 )
    !$ACC UPDATE ASYNC(1) HOST(list_seawtr_count, list_seaice_count) ! also update index lists?
    IF (.NOT. lcuda_graph_lnd) THEN
      !$ACC WAIT(1)
    END IF
    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA

  END SUBROUTINE update_idx_lists_sea


  !> sstice_mode = MODE SSTICE_ANA_CLINC
  !!
  !! Climatological SST increments are computed, based on the climatological 
  !! SST fields read from the external parameter file. The increment is defined 
  !! as the climatological SST difference between the actual date (target_datetime) 
  !! and the experiment start date (ref_datetime). 
  !! The increment is used to update the surface temperature of the sea water tile 
  !! (t_g_t, t_s_t, t_sk_t). 
  !! Please note that the original SST field t_seasfc remains unchanged, i.e. 
  !! the SST increment is NOT added to t_seasfc.
  !!
  SUBROUTINE sst_add_climatological_incr (p_patch, ext_data, prog_lnd, diag_lnd, pres_sfc, &
    &                       ref_datetime, target_datetime)

    TYPE(t_patch),           INTENT(IN)    :: p_patch
    TYPE(t_external_data),   INTENT(INOUT) :: ext_data
    TYPE(t_lnd_prog),        INTENT(INOUT) :: prog_lnd(:)      !< prog vars for sfc
    TYPE(t_lnd_diag),        INTENT(INOUT) :: diag_lnd         !< diag vars for sfc
    REAL(wp),                INTENT(IN)    :: pres_sfc(:,:)    !< surface pressure
    TYPE(datetime),          INTENT(IN)    :: ref_datetime     !< experiment start datetime
    TYPE(datetime),          INTENT(IN)    :: target_datetime  !< actual datetime

    ! Local scalars:
    !
    INTEGER :: jb, ic, jc
    INTEGER :: jg
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk
    INTEGER :: n_now, n_new
    INTEGER :: ierr
    REAL(wp):: sst_cl_inc
    REAL(wp):: new_sst            ! updated SST value
    REAL(wp):: max_inc, min_inc   ! max/min SST increment on given PE

    REAL(wp)::  &                 ! climatological sst field for the model initialization day
      &  sst_cl_ini_day(nproma,p_patch%nblks_c)

    REAL(wp)::  &                 ! climatological sst field for the current day
      &  sst_cl_cur_day(nproma,p_patch%nblks_c)

    REAL(wp)::  &                 ! climatological SST increment
      &  sst_inc(nproma,p_patch%nblks_c)

    CHARACTER(len=MAX_DATETIME_STR_LEN) :: target_datetime_PTString
    CHARACTER(len=MAX_DATETIME_STR_LEN) :: ref_datetime_PTString
    CHARACTER(len=*), PARAMETER :: routine = 'sst_add_climatological_incr'
!_______________

    jg = p_patch%id

    n_now = nnow_rcf(jg)
    n_new = nnew_rcf(jg)

    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)


    ! MODE SSTICE_ANA_CLINC:
    ! - SST and sea ice fraction have been read from the analysis. 
    ! - The SST (t_g_t, t_s_t, t_sk_t) is updated by climatological SST increments on a daily basis. 
    ! - The sea ice fraction is not updated by any external information. 
    !   It may only change due to the melting of seaice points (see NWP seaice model).
    !
    IF (msg_level >= 13) THEN
      WRITE(message_text,'(a)') 'Update SST with climatological increments'
      CALL message(routine, TRIM(message_text))
      !
      CALL datetimeToString(target_datetime, target_datetime_PTString, ierr)
      CALL datetimeToString(ref_datetime, ref_datetime_PTString, ierr)
      WRITE(message_text,'(a,i2,a,a)') 'Target Date for DOM ',jg,': ',TRIM(target_datetime_PTString) 
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i2,a,a)') 'Reference Date for DOM ',jg,': ',TRIM(ref_datetime_PTString) 
      CALL message('', TRIM(message_text))
    ENDIF


    ! get climatological sst for initial and current day
    !
    CALL interpol_monthly_mean(p_patch, ref_datetime,         &! in
      &                        ext_data%atm_td%sst_m,         &! in
      &                        sst_cl_ini_day                 )! out

    CALL interpol_monthly_mean(p_patch, target_datetime,      &! in
      &                        ext_data%atm_td%sst_m,         &! in
      &                        sst_cl_cur_day                 )! out


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,ic,jc,sst_cl_inc,new_sst)
    DO jb=i_startblk, i_endblk

      ! loop over all open water points and add climatological increments
!NEC$ ivdep
      DO ic = 1, ext_data%atm%list_seawtr%ncount(jb)
        jc = ext_data%atm%list_seawtr%idx(ic,jb)

        sst_cl_inc = sst_cl_cur_day(jc,jb) - sst_cl_ini_day(jc,jb)
        ! make sure, that the updated SST does not drop below 
        ! the salt-water freezing point
        new_sst = MAX(tf_salt,diag_lnd%t_seasfc(jc,jb) + sst_cl_inc)
        !
        prog_lnd(n_now)%t_g_t (jc,jb,isub_water) = new_sst
        prog_lnd(n_now)%t_s_t (jc,jb,isub_water) = new_sst
        prog_lnd(n_now)%t_sk_t(jc,jb,isub_water) = new_sst
        prog_lnd(n_new)%t_g_t (jc,jb,isub_water) = new_sst
        prog_lnd(n_new)%t_s_t (jc,jb,isub_water) = new_sst
        prog_lnd(n_new)%t_sk_t(jc,jb,isub_water) = new_sst

        ! includes reduction of saturation pressure due to salt content
        diag_lnd%qv_s_t(jc,jb,isub_water) = salinity_fac *                       & 
          &   spec_humi( sat_pres_water(prog_lnd(n_now)%t_g_t(jc,jb,isub_water)),&
          &                                  pres_sfc(jc,jb) )

      ENDDO  ! ic
    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL

    ! aggregate updated t_g_t and qv_s_t
    CALL aggregate_tg_qvs( p_patch = p_patch,                       & ! in
      &                    frac_t  = ext_data%atm%frac_t(:,:,:),    & ! in
      &                    t_g_t   = prog_lnd(n_now)%t_g_t(:,:,:),  & ! in
      &                    qv_s_t  = diag_lnd%qv_s_t(:,:,:),        & ! in
      &                    t_g     = prog_lnd(n_now)%t_g(:,:),      & ! inout
      &                    qv_s    = diag_lnd%qv_s(:,:)             ) ! inout

    ! debug output
    IF (msg_level >= 13) THEN
      sst_inc(:,:) = 0._wp
      DO jb=i_startblk, i_endblk
        ! loop over all open water points and add climatological increments
        DO ic = 1, ext_data%atm%list_seawtr%ncount(jb)
          jc = ext_data%atm%list_seawtr%idx(ic,jb)
          sst_inc(jc,jb) =  sst_cl_cur_day(jc,jb) - sst_cl_ini_day(jc,jb)
        ENDDO
      ENDDO
      max_inc = MAXVAL(sst_inc(:,:))
      min_inc = MINVAL(sst_inc(:,:))
      !
      WRITE(message_text,'(2(a,i2,a,e12.5))') 'max increment DOM', jg, ': ', global_max(max_inc), &
        &                                     ', min increment DOM', jg, ': ', global_min(min_inc)
      CALL message('', TRIM(message_text))
    ENDIF

  END SUBROUTINE sst_add_climatological_incr



  !> sstice_mode = SSTICE_INST, SSTICE_CLIM, SSTICE_AVG_MONTHLY, SSTICE_AVG_DAILY, 
  !!
  !! Processes the fields fr_seaice, t_seasfc and optionally h_ice, which are 
  !! assumed to be provided from external sources (e.g. read from file or received 
  !! from coupled ocean model at specific intervals).
  !!
  !! The dynamic index lists for seaice and open water are re-generated on the basis 
  !! of fr_seaice. Based on the updated index lists, several water tile temperature 
  !! fields are updated with t_seasfc. Optionally, the h_ice field of the NWP seaice 
  !! scheme is overwritten by the optionally provided h_ice field.
  !! 
  !! Updated fields: 
  !! general        : t_g_t, t_s_t, t_sk_t, qv_s_t
  !! seaice-specific: h_ice, t_ice, t_snow_si, h_snow_si, alb_si
  !!
  SUBROUTINE process_sst_and_seaice (p_patch, fr_seaice, t_seasfc, pres_sfc, ext_data,&
    &                               prog_lnd_now, prog_lnd_new, prog_wtr_now, prog_wtr_new, &
    &                               diag_lnd, optin_h_ice, lacc)

    TYPE(t_patch),           INTENT(IN)    :: p_patch
    REAL(wp),                INTENT(INOUT) :: fr_seaice(:,:)    !< sea ice fraction from 
                                                                !  external sources 
    REAL(wp),                INTENT(IN)    :: t_seasfc(:,:)     !< sea surface temperature from
                                                                !  external sources
    REAL(wp),                INTENT(IN)    :: pres_sfc(:,:)     !< surface pressure 
    TYPE(t_external_data),   INTENT(INOUT) :: ext_data
    TYPE(t_lnd_prog),        INTENT(INOUT) :: prog_lnd_now      !< prog vars for sfc
    TYPE(t_lnd_prog),        INTENT(INOUT) :: prog_lnd_new
    TYPE(t_wtr_prog),        INTENT(INOUT) :: prog_wtr_now      !< prog vars for sfc
    TYPE(t_wtr_prog),        INTENT(INOUT) :: prog_wtr_new
    TYPE(t_lnd_diag),        INTENT(INOUT) :: diag_lnd          !< diag vars for sfc
    REAL(wp), OPTIONAL,      INTENT(IN)    :: optin_h_ice(:,:)  !< ice thickness from 
                                                                !  external sources
    LOGICAL,  OPTIONAL,      INTENT(IN)    :: lacc

    ! Local scalars:
    !
    INTEGER :: jb, ic, jc, ncount
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    REAL(wp):: t_water
    LOGICAL :: lpresent_h_ice

    REAL(wp) :: frsi     (nproma)   ! sea ice fraction
    REAL(wp) :: tice_now (nproma)   ! temperature of ice upper surface at previous time
    REAL(wp) :: hice_now (nproma)   ! ice thickness at previous time level
    REAL(wp) :: tsnow_now(nproma)   ! temperature of snow upper surface at previous time
    REAL(wp) :: hsnow_now(nproma)   ! snow thickness at previous time level
    REAL(wp) :: albsi_now(nproma)   ! sea-ice albedo at previous time level
    REAL(wp) :: tice_new (nproma)   ! temperature of ice upper surface at new time
    REAL(wp) :: hice_new (nproma)   ! ice thickness at new time level
    REAL(wp) :: tsnow_new(nproma)   ! temperature of snow upper surface at new time
    REAL(wp) :: hsnow_new(nproma)   ! snow thickness at new time level
    REAL(wp) :: albsi_new(nproma)   ! sea-ice albedo at new time level

    TYPE(t_idx_list_blocked) :: list_seaice_new, list_seaice_old
    TYPE(t_idx_list_blocked) :: list_water_new, list_water_old

    TYPE(t_idx_list_blocked) :: list_seaice_created
    TYPE(t_idx_list_blocked) :: list_seaice_destroyed
    TYPE(t_idx_list_blocked) :: list_seaice_retained
    TYPE(t_idx_list_blocked) :: list_water_created
    TYPE(t_idx_list_blocked) :: list_water_destroyed
    TYPE(t_idx_list_blocked) :: list_water_retained

    CHARACTER(len=*), PARAMETER :: routine = 'mo_nwp_sfc_interface:process_sst_and_seaice'

    LOGICAL :: lzacc

!----------------------------------------------------------------
    CALL set_acc_host_or_device(lzacc, lacc)

    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    lpresent_h_ice = PRESENT(optin_h_ice)

    IF (lseaice) THEN

      !$ACC DATA CREATE(frsi, tice_now, hice_now, tsnow_now, hsnow_now, albsi_now) &
      !$ACC   CREATE(tice_new, hice_new, tsnow_new, hsnow_new, albsi_new) IF(lzacc)

      ! allocate index lists for seaice and open water points
      !
      ! The following are not created on GPU as only needed for check_water_idx_lists,
      ! which would require an unnecessary back copy to CPU and the porting of compare_sets
      CALL list_seaice_new%construct(nproma, p_patch%nblks_c)
      CALL list_water_new%construct (nproma, p_patch%nblks_c)
      !
      CALL list_seaice_old%construct(nproma, p_patch%nblks_c)
      CALL list_water_old%construct (nproma, p_patch%nblks_c)

      ! store old seaice and open water lists
      !
      ! seaice list
      !$ACC WAIT
      !$ACC UPDATE HOST(ext_data%atm%list_seaice%ncount, ext_data%atm%list_seaice%idx) IF(lzacc)
      list_seaice_old = ext_data%atm%list_seaice
      !
      ! water list
      !$ACC UPDATE HOST(ext_data%atm%list_seawtr%ncount, ext_data%atm%list_seawtr%idx) IF(lzacc)
      list_water_old = ext_data%atm%list_seawtr


      !
      ! generate new sea-ice and open water lists, 
      ! and set frac_t and sai_t
      !
      CALL init_sea_lists(p_patch      = p_patch,           & ! in
        &                 lseaice      = lseaice,           & ! in
        &                 fr_seaice    = fr_seaice(:,:),    & ! in(out)
        &                 ext_data     = ext_data,          & ! inout
        &                 opt_lverbose = .FALSE.,           & ! in
        &                 lacc         = lzacc)               ! in


      ! store updated index lists
      ! this is not strictly necessary. We may equally well work directly with 
      ! ext_data%atm%list_seaice
      ! ext_data%atm%list_seawtr
      !
      ! seaice list
      !$ACC WAIT
      !$ACC UPDATE HOST(ext_data%atm%list_seaice%ncount, ext_data%atm%list_seaice%idx) IF(lzacc)
      list_seaice_new = ext_data%atm%list_seaice
      !
      ! water list
      !$ACC UPDATE HOST(ext_data%atm%list_seawtr%ncount, ext_data%atm%list_seawtr%idx) IF(lzacc)
      list_water_new = ext_data%atm%list_seawtr


      ! compare old and new index lists by grouping the elements 
      ! into one of the following three sublists
      ! list_XY_retained  : element exists both in old and new list 
      ! list_XY_destroyed : element exists only in old list
      ! list_XY_created   : element exists only in new list
      !
      ! The following are created on GPU after CALL compare_sets (except list_seaice_retained)
      ! to avoid unnecessary back and forth copies. All input arrays to compare_sets need to be on CPU
      ! for later call to check_water_idx_lists but not all are needed on GPU
      CALL list_water_retained%construct(nproma, p_patch%nblks_c)
      CALL list_water_destroyed%construct(nproma, p_patch%nblks_c)
      CALL list_water_created%construct(nproma, p_patch%nblks_c)
      !
      CALL list_seaice_retained%construct(nproma, p_patch%nblks_c)
      CALL list_seaice_destroyed%construct(nproma, p_patch%nblks_c)
      CALL list_seaice_created%construct(nproma, p_patch%nblks_c)

      CALL compare_sets(p_patch        = p_patch,              &
        &               list1          = list_water_old,       &
        &               list2          = list_water_new,       &
        &               list_intersect = list_water_retained,  &
        &               list1_only     = list_water_destroyed, &
        &               list2_only     = list_water_created    )

      CALL compare_sets(p_patch        = p_patch,              &
        &               list1          = list_seaice_old,      &
        &               list2          = list_seaice_new,      &
        &               list_intersect = list_seaice_retained, &
        &               list1_only     = list_seaice_destroyed,&
        &               list2_only     = list_seaice_created   )
      !$ACC ENTER DATA CREATE(list_water_retained, list_water_destroyed) &
      !$ACC   CREATE(list_water_created, list_seaice_destroyed, list_seaice_created) IF(lzacc)
      !$ACC ENTER DATA COPYIN(list_water_retained%idx, list_water_retained%ncount) &
      !$ACC   COPYIN(list_water_destroyed%idx, list_water_destroyed%ncount) &
      !$ACC   COPYIN(list_water_created%idx, list_water_created%ncount) &
      !$ACC   COPYIN(list_seaice_destroyed%idx, list_seaice_destroyed%ncount) &
      !$ACC   COPYIN(list_seaice_created%idx, list_seaice_created%ncount) IF(lzacc)


      ! update various water and seaice related fields depending on whether 
      ! a seaice/water point is 
      ! * retained
      ! * newly generated
      ! * destroyed
      !
      ! Note that for ntiles==1 the following lists are equivalent
      ! list_seaice_created == list_water_destroyed
      ! list_water_created  == list_seaice_destroyed
      ! Furthermore, isub_seaice==isub_water==1 holds. 
      ! As a consequence, for ntiles==1 the (cosmetic) updates of t_g_t, t_s_t, qv_s_t are 
      ! skipped for all list_XYZ_destroyed lists, as list_XYZ_destroyed and list_XYZ_created 
      ! operate on exactly the same fields/memory. Otherwise there is the risk that we overwrite 
      ! meaningful results.
      !
      ! For ntiles>1 isub_seaice/=isub_water. I.e. no such risk exists in case of activated tile approach.
      !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,ic,ncount,i_startidx,i_endidx,frsi,tice_now,hice_now,tsnow_now,hsnow_now,albsi_now,  &
!$OMP            tice_new,hice_new,tsnow_new,hsnow_new,albsi_new,t_water), SCHEDULE(guided)
      DO jb = i_startblk, i_endblk


        IF (lpresent_h_ice) THEN
#ifdef _OPENACC
          CALL finish('process_sst_and_seaice', "lseaice together with lpresent_h_ice is not supported on GPU")
#endif

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
            &                i_startidx, i_endidx, rl_start, rl_end)

          !
          ! Overwrite h_ice with externally provided values
          !
          ! We assume that h_ice is consistent with fr_seaice, and that this consistency 
          ! has been checked elsewhere. 
          ! 
          DO jc = i_startidx, i_endidx
            prog_wtr_now%h_ice(jc,jb) = optin_h_ice(jc,jb)
            prog_wtr_new%h_ice(jc,jb) = optin_h_ice(jc,jb)
          ENDDO
        ENDIF


        !*************************!
        !    seaice points        !
        !*************************!

        ! I) retained seaice points
        !

        ! no further action required
        ! handled by the seaice scheme


        ! II) newly generated seaice points
        !
        ! perform seaice warmstart (compare with nwp_surface_init)
        !
        ncount = list_seaice_created%ncount(jb) ! CPU and GPU values in sync
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) FIRSTPRIVATE(ncount) IF(lzacc)
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jc)
        DO ic = 1, ncount

          jc = list_seaice_created%idx(ic,jb)

          frsi     (ic) = fr_seaice(jc,jb)
          tice_now (ic) = prog_wtr_now%t_ice(jc,jb)
          hice_now (ic) = prog_wtr_now%h_ice(jc,jb)
          tsnow_now(ic) = prog_wtr_now%t_snow_si(jc,jb)
          hsnow_now(ic) = prog_wtr_now%h_snow_si(jc,jb)
          albsi_now(ic) = prog_wtr_now%alb_si(jc,jb)
        ENDDO  ! jc
        !$ACC END PARALLEL

        ! If h_ice is provided from external sources (lpresent_h_ice=.TRUE.), 
        ! it is necessary to skip the initialization of h_ice for newly 
        ! generated seaice points. 'Newly generated' as seen from an 
        ! atmospheric model perspective does not necessarily mean that 
        ! the seaice was created due to freezing. It might have been transported into the 
        ! cell by advective processes.
        !
        CALL seaice_init_nwp ( (.NOT. lpresent_h_ice),                              & ! in 
          &                    list_seaice_created%ncount(jb), frsi,                & ! in
          &                    tice_now, hice_now, tsnow_now, hsnow_now, albsi_now, & ! inout
          &                    tice_new, hice_new, tsnow_new, hsnow_new, albsi_new, & ! inout
          &                    lacc=lzacc)


        !  Recover fields from index list
        !
        ncount = list_seaice_created%ncount(jb) ! CPU and GPU values in sync
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) FIRSTPRIVATE(ncount) IF(lzacc)
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jc)
        DO ic = 1, ncount

          jc = list_seaice_created%idx(ic,jb)

          prog_wtr_now%t_ice(jc,jb)    = tice_now(ic)
          prog_wtr_now%h_ice(jc,jb)    = hice_now(ic)
          prog_wtr_now%t_snow_si(jc,jb)= tsnow_now(ic)
          prog_wtr_now%h_snow_si(jc,jb)= hsnow_now(ic)

          prog_wtr_new%t_ice(jc,jb)    = tice_new(ic)
          prog_wtr_new%h_ice(jc,jb)    = hice_new(ic)
          prog_wtr_new%t_snow_si(jc,jb)= tsnow_new(ic)
          prog_wtr_new%h_snow_si(jc,jb)= hsnow_new(ic)

          IF (lprog_albsi) THEN
            prog_wtr_now%alb_si(jc,jb) = albsi_now(ic)
            prog_wtr_new%alb_si(jc,jb) = albsi_new(ic)
          ENDIF

          prog_lnd_now%t_g_t (jc,jb,isub_seaice) = tice_now(ic)
          prog_lnd_new%t_g_t (jc,jb,isub_seaice) = tice_new(ic) ! == tice_now(ic)
          prog_lnd_now%t_s_t (jc,jb,isub_seaice) = tf_salt
          prog_lnd_new%t_s_t (jc,jb,isub_seaice) = tf_salt
          prog_lnd_now%t_sk_t(jc,jb,isub_seaice) = tf_salt
          prog_lnd_new%t_sk_t(jc,jb,isub_seaice) = tf_salt
          diag_lnd%qv_s_t(jc,jb,isub_seaice)     = spec_humi(sat_pres_ice(tice_now(ic)),&
            &                                                pres_sfc(jc,jb) )
        ENDDO  ! ic
        !$ACC END PARALLEL


        ! III) destroyed seaice points
        !      now pure water point, no seaice tile
        !
        ! re-initialize h_ice, t_ice, h_snow_si, t_snow_si, alb_si
        ncount = list_seaice_destroyed%ncount(jb) ! CPU and GPU values in sync
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) FIRSTPRIVATE(ncount) IF(lzacc)
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jc)
        DO ic = 1, ncount

          jc = list_seaice_destroyed%idx(ic,jb)

          prog_wtr_now%h_ice(jc,jb) = 0._wp
          prog_wtr_new%h_ice(jc,jb) = 0._wp
          !
          prog_wtr_now%t_ice(jc,jb) = tmelt
          prog_wtr_new%t_ice(jc,jb) = tmelt
          !
          prog_wtr_now%h_snow_si(jc,jb)= 0._wp
          prog_wtr_new%h_snow_si(jc,jb)= 0._wp
          !
          prog_wtr_now%t_snow_si(jc,jb)= tmelt
          prog_wtr_new%t_snow_si(jc,jb)= tmelt

          IF (lprog_albsi) THEN
            ! re-initialize prognostic seaice albedo with ALB_SI_MISSVAL
            prog_wtr_now%alb_si(jc,jb) = ALB_SI_MISSVAL
            prog_wtr_new%alb_si(jc,jb) = ALB_SI_MISSVAL
          ENDIF

          IF (ntiles_total > 1) THEN
            prog_lnd_now%t_g_t (jc,jb,isub_seaice) = tmelt
            prog_lnd_new%t_g_t (jc,jb,isub_seaice) = tmelt
            prog_lnd_now%t_s_t (jc,jb,isub_seaice) = tmelt
            prog_lnd_new%t_s_t (jc,jb,isub_seaice) = tmelt
            prog_lnd_now%t_sk_t(jc,jb,isub_seaice) = tmelt
            prog_lnd_new%t_sk_t(jc,jb,isub_seaice) = tmelt
            !diag_lnd%qv_s_t(jc,jb,isub_seaice)     =
          ENDIF
        ENDDO  ! ic
        !$ACC END PARALLEL


        !*************************!
        !    water points         !
        !*************************!

        ! I) retained water points
        !    update SST
        !
        ncount = list_water_retained%ncount(jb) ! CPU and GPU values in sync
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) FIRSTPRIVATE(ncount) IF(lzacc)
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jc)
        DO ic = 1, ncount

          jc = list_water_retained%idx(ic,jb)

          t_water = t_seasfc(jc,jb)
          prog_lnd_now%t_g_t (jc,jb,isub_water)= t_water
          prog_lnd_now%t_s_t (jc,jb,isub_water)= t_water
          prog_lnd_now%t_sk_t(jc,jb,isub_water)= t_water
          prog_lnd_new%t_g_t (jc,jb,isub_water)= t_water
          prog_lnd_new%t_s_t (jc,jb,isub_water)= t_water
          prog_lnd_new%t_sk_t(jc,jb,isub_water)= t_water
          !
          ! includes reduction of saturation pressure due to salt content
          diag_lnd%qv_s_t(jc,jb,isub_water)    =  salinity_fac *      &
            &                             spec_humi( sat_pres_water(t_water ),    &
            &                             pres_sfc(jc,jb) )

        ENDDO  ! ic
        !$ACC END PARALLEL


        ! II) newly generated water points
        !     update SST
        !
        ncount = list_water_created%ncount(jb) ! CPU and GPU values in sync
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) FIRSTPRIVATE(ncount) IF(lzacc)
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jc, t_water)
        DO ic = 1, ncount

          jc = list_water_created%idx(ic,jb)

          t_water = MAX(tf_salt,t_seasfc(jc,jb))
          prog_lnd_now%t_g_t (jc,jb,isub_water)= t_water
          prog_lnd_now%t_s_t (jc,jb,isub_water)= t_water
          prog_lnd_now%t_sk_t(jc,jb,isub_water)= t_water
          prog_lnd_new%t_g_t (jc,jb,isub_water)= t_water
          prog_lnd_new%t_s_t (jc,jb,isub_water)= t_water
          prog_lnd_new%t_sk_t(jc,jb,isub_water)= t_water

          !
          ! includes reduction of saturation pressure due to salt content
          diag_lnd%qv_s_t(jc,jb,isub_water)    =  salinity_fac *      &
            &                             spec_humi( sat_pres_water(t_water ),    &
            &                             pres_sfc(jc,jb) )

        ENDDO  ! ic
        !$ACC END PARALLEL


        ! The following part is skipped for the case ntiles == 1
        ! Since 
        ! * isub_water == isub_seaice == 1, and 
        ! * list_water_destroyed = list_seaice_created 
        ! we would otherwise overwrite t_g_t, t_s_t, qv_s_t, which has already been 
        ! initialized for newly generated seaice points (list_seaice_created).
        IF (ntiles_total > 1) THEN
          !
          ! III) destroyed water points
          !      now pure seaice points
          !
          ! re-initialize water temperature which is in contact with overlying seaice
          ncount = list_water_destroyed%ncount(jb) ! CPU in sync with GPU value
          !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) FIRSTPRIVATE(ncount) IF(lzacc)
          !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jc)
          DO ic = 1, ncount

            jc = list_water_destroyed%idx(ic,jb)

            prog_lnd_now%t_g_t (jc,jb,isub_water)= tf_salt
            prog_lnd_now%t_s_t (jc,jb,isub_water)= tf_salt
            prog_lnd_now%t_sk_t(jc,jb,isub_water)= tf_salt
            prog_lnd_new%t_g_t (jc,jb,isub_water)= tf_salt
            prog_lnd_new%t_s_t (jc,jb,isub_water)= tf_salt
            prog_lnd_new%t_sk_t(jc,jb,isub_water)= tf_salt
            !
            ! includes reduction of saturation pressure due to salt content
            diag_lnd%qv_s_t(jc,jb,isub_water)    =  salinity_fac *      &
              &                             spec_humi( sat_pres_water(tf_salt ),    &
              &                             pres_sfc(jc,jb) )

          ENDDO  ! ic
          !$ACC END PARALLEL
        END IF  ! ntiles_total > 1

      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      IF (msg_level >= 12) THEN
        ! sanity check and log output
        !$ACC WAIT
        !$ACC UPDATE HOST(list_water_retained%idx, list_water_retained%ncount) &
        !$ACC   HOST(list_water_destroyed%idx, list_water_destroyed%ncount) &
        !$ACC   HOST(list_water_created%idx, list_water_created%ncount) &
        !$ACC   HOST(list_seaice_destroyed%idx, list_seaice_destroyed%ncount) &
        !$ACC   HOST(list_seaice_created%idx, list_seaice_created%ncount) IF(lzacc)
        CALL check_water_idx_lists(list_seaice_old       = list_seaice_old,       &
          &                        list_seaice_new       = list_seaice_new,       &
          &                        list_water_old        = list_water_old,        &
          &                        list_water_new        = list_water_new,        &
          &                        list_water_created    = list_water_created,    &
          &                        list_water_destroyed  = list_water_destroyed,  &
          &                        list_water_retained   = list_water_retained,   &
          &                        list_seaice_created   = list_seaice_created,   &
          &                        list_seaice_destroyed = list_seaice_destroyed, &
          &                        list_seaice_retained  = list_seaice_retained,  &
          &                        i_startblk            = i_startblk,            &
          &                        i_endblk              = i_endblk,              &
          &                        ltile_approach        = (ntiles_total>1)       )
      ENDIF


      ! cleanup
      !$ACC WAIT
      !$ACC EXIT DATA DELETE(list_water_retained%idx, list_water_retained%ncount) &
      !$ACC   DELETE(list_water_destroyed%idx, list_water_destroyed%ncount) &
      !$ACC   DELETE(list_water_created%idx, list_water_created%ncount) &
      !$ACC   DELETE(list_seaice_destroyed%idx, list_seaice_destroyed%ncount) &
      !$ACC   DELETE(list_seaice_created%idx, list_seaice_created%ncount) IF(lzacc)
      !$ACC EXIT DATA DELETE(list_water_retained, list_water_destroyed) &
      !$ACC   DELETE(list_water_created, list_seaice_destroyed, list_seaice_created) IF(lzacc)
      CALL list_seaice_new%finalize()
      CALL list_seaice_old%finalize()
      CALL list_seaice_created%finalize()
      CALL list_seaice_destroyed%finalize()
      CALL list_seaice_retained%finalize()
      !
      CALL list_water_new%finalize()
      CALL list_water_old%finalize()
      CALL list_water_created%finalize()
      CALL list_water_destroyed%finalize()
      CALL list_water_retained%finalize()

      !$ACC WAIT
      !$ACC END DATA ! local IF lpresent_h_ice

    ELSE   ! seaice model switched off

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO jb = i_startblk, i_endblk

!$NEC ivdep
        !$ACC LOOP GANG VECTOR PRIVATE(jc)
        DO ic = 1, ext_data%atm%list_sea%ncount(jb)

          jc = ext_data%atm%list_sea%idx(ic,jb)
          ! only if the dominant tile is water t_g_t is set to t_seasfc
          IF ( fr_seaice(jc,jb) < 0.5_wp ) THEN
            prog_lnd_now%t_g_t (jc,jb,isub_water) = t_seasfc(jc,jb)
            prog_lnd_now%t_s_t (jc,jb,isub_water) = t_seasfc(jc,jb)
            prog_lnd_now%t_sk_t(jc,jb,isub_water) = t_seasfc(jc,jb)
            prog_lnd_new%t_g_t (jc,jb,isub_water) = t_seasfc(jc,jb)
            prog_lnd_new%t_s_t (jc,jb,isub_water) = t_seasfc(jc,jb)
            prog_lnd_new%t_sk_t(jc,jb,isub_water) = t_seasfc(jc,jb)
            ! includes reduction of saturation pressure due to salt content
            diag_lnd%qv_s_t(jc,jb,isub_water)    =  salinity_fac *        &
              &    spec_humi( sat_pres_water(t_seasfc(jc,jb) ), pres_sfc(jc,jb) )
          END IF
        ENDDO  ! ic

      ENDDO  ! jb
      !$ACC END PARALLEL

    ENDIF  ! lseaice

    ! aggregate updated t_g_t and qv_s_t
    !
    CALL aggregate_tg_qvs( p_patch = p_patch,                    & ! in
      &                    frac_t  = ext_data%atm%frac_t(:,:,:), & ! in
      &                    t_g_t   = prog_lnd_now%t_g_t(:,:,:),  & ! in
      &                    qv_s_t  = diag_lnd%qv_s_t(:,:,:),     & ! in
      &                    t_g     = prog_lnd_now%t_g(:,:),      & ! inout
      &                    qv_s    = diag_lnd%qv_s(:,:),         & ! inout
      &                    lacc    = lzacc                       )

    CALL aggregate_tg_qvs( p_patch = p_patch,                    & ! in
      &                    frac_t  = ext_data%atm%frac_t(:,:,:), & ! in
      &                    t_g_t   = prog_lnd_new%t_g_t(:,:,:),  & ! in
      &                    qv_s_t  = diag_lnd%qv_s_t(:,:,:),     & ! in
      &                    t_g     = prog_lnd_new%t_g(:,:),      & ! inout
      &                    qv_s    = diag_lnd%qv_s(:,:),         & ! inout
      &                    lacc    = lzacc                       )

  END SUBROUTINE process_sst_and_seaice


  !! After updating the SST and sea ice fraction (from external files),
  !! the dynamic index lists for seaice and open water are re-generated.
  !! Based on these index lists several seaice and water-related fields are 
  !! updated.
  !! Updated fields: 
  !! general        : t_g_t, t_s_t, qv_s_t
  !! seaice-specific: h_ice, t_ice, t_snow_si, h_snow_si, alb_si
  !!
  SUBROUTINE check_water_idx_lists(list_seaice_old, list_seaice_new, &
    &                              list_water_old, list_water_new, list_water_created, &
    &                              list_water_destroyed, list_water_retained,  &
    &                              list_seaice_created, list_seaice_destroyed, list_seaice_retained, &
    &                              i_startblk, i_endblk, ltile_approach)

    TYPE(t_idx_list_blocked), INTENT(IN) :: list_seaice_old
    TYPE(t_idx_list_blocked), INTENT(IN) :: list_seaice_new
    TYPE(t_idx_list_blocked), INTENT(IN) :: list_water_old
    TYPE(t_idx_list_blocked), INTENT(IN) :: list_water_new
    TYPE(t_idx_list_blocked), INTENT(IN) :: list_water_created
    TYPE(t_idx_list_blocked), INTENT(IN) :: list_water_destroyed
    TYPE(t_idx_list_blocked), INTENT(IN) :: list_water_retained
    TYPE(t_idx_list_blocked), INTENT(IN) :: list_seaice_created
    TYPE(t_idx_list_blocked), INTENT(IN) :: list_seaice_destroyed
    TYPE(t_idx_list_blocked), INTENT(IN) :: list_seaice_retained
    INTEGER,                  INTENT(IN) :: i_startblk, i_endblk
    LOGICAL,                  INTENT(IN) :: ltile_approach

    ! local
    INTEGER :: npoints_ice_old, npoints_ice_new
    INTEGER :: npoints_wtr_old, npoints_wtr_new
    INTEGER :: npoints_wtr_created, npoints_wtr_destroyed, npoints_wtr_retained
    INTEGER :: npoints_ice_created, npoints_ice_destroyed, npoints_ice_retained
    INTEGER :: irow
    INTEGER :: residuum_ice, residuum_wtr
    INTEGER :: residuum_wtr2ice, residuum_ice2wtr

    ! table
    TYPE(t_table)   :: table
    CHARACTER(LEN = *), PARAMETER :: pts_type      = "gridpoint type",   & 
      &                              pts_total     = "total", &
      &                              pts_diff      = "diff (new - old)",  &
      &                              pts_retained  = "retained", &
      &                              pts_created   = "created", &
      &                              pts_destroyed = "destroyed"
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: type_name

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nwp_sfc_utils: check_water_idx_lists'

    !--------------------------------------------------------------------

    npoints_ice_old = list_seaice_old%get_sum_global(i_startblk,i_endblk)
    npoints_ice_new = list_seaice_new%get_sum_global(i_startblk,i_endblk)

    npoints_wtr_old = list_water_old%get_sum_global(i_startblk,i_endblk)
    npoints_wtr_new = list_water_new%get_sum_global(i_startblk,i_endblk)

    npoints_wtr_created   = list_water_created   %get_sum_global(i_startblk,i_endblk)
    npoints_wtr_destroyed = list_water_destroyed %get_sum_global(i_startblk,i_endblk)
    npoints_wtr_retained  = list_water_retained  %get_sum_global(i_startblk,i_endblk)

    npoints_ice_created   = list_seaice_created   %get_sum_global(i_startblk,i_endblk)
    npoints_ice_destroyed = list_seaice_destroyed %get_sum_global(i_startblk,i_endblk)
    npoints_ice_retained  = list_seaice_retained  %get_sum_global(i_startblk,i_endblk)


    ! the following will only be executed by stdio process
    IF(.NOT. my_process_is_stdio()) RETURN


    ! sanity checks
    !
    ! For each gridpoint type (water, seaice) the total change of 
    ! gridpoints (new-old) must match the sum of destroyed and created grid points.
    residuum_ice = (npoints_ice_new-npoints_ice_old) - (npoints_ice_created-npoints_ice_destroyed)
    residuum_wtr = (npoints_wtr_new-npoints_wtr_old) - (npoints_wtr_created-npoints_wtr_destroyed)

    IF ( (residuum_ice/=0) .OR. (residuum_wtr/=0) ) THEN
      WRITE(message_text,'(a)') &
        & 'Change in number of total points does not match the sum of created and destroyed points)'
      CALL finish (TRIM(routine), TRIM(message_text))
    ENDIF

    ! Without tile aproach, the total number of created (destroyed) seaice points 
    ! must equal the total number of destroyed (created) water points.
    IF (.NOT. ltile_approach) THEN
      residuum_wtr2ice = npoints_ice_created - npoints_wtr_destroyed
      residuum_ice2wtr = npoints_wtr_created - npoints_ice_destroyed
      IF ( (residuum_wtr2ice /= 0) .OR. (residuum_ice2wtr /= 0) ) THEN
        WRITE(message_text,'(a)') &
          & 'Number of created (destroyed) seaice points does not match number of destroyed (created) water points'
        CALL finish (TRIM(routine), TRIM(message_text))
      ENDIF  
    ENDIF


    ! printout

    IF (msg_level >= 13) THEN
      WRITE (0,*) " " ! newline
      WRITE(message_text,'(a)') 'Sanity check for updated water/seaice index lists'
      CALL message('', TRIM(message_text))

      ! table-based output
      CALL initialize_table(table)
      ! the latter is no longer mandatory
      CALL add_table_column(table, pts_type)
      CALL add_table_column(table, pts_total)
      CALL add_table_column(table, pts_diff)
      CALL add_table_column(table, pts_retained)
      CALL add_table_column(table, pts_created)
      CALL add_table_column(table, pts_destroyed)

      ! column 1
      irow = 1
      type_name = "seaice"
      CALL set_table_entry(table,irow,pts_type     , ADJUSTL(TRIM(type_name)))
      CALL set_table_entry(table,irow,pts_total    , ADJUSTL(TRIM(int2string(npoints_ice_new)) ))
      CALL set_table_entry(table,irow,pts_diff     , ADJUSTL(TRIM(int2string(npoints_ice_new-npoints_ice_old)) ))
      CALL set_table_entry(table,irow,pts_retained , ADJUSTL(TRIM(int2string(npoints_ice_retained)) ))
      CALL set_table_entry(table,irow,pts_created  , ADJUSTL(TRIM(int2string(npoints_ice_created)) ))
      CALL set_table_entry(table,irow,pts_destroyed, ADJUSTL(TRIM(int2string(npoints_ice_destroyed)) ))

      irow = 2
      type_name = "water"
      CALL set_table_entry(table,irow,pts_type     , ADJUSTL(TRIM(type_name)))
      CALL set_table_entry(table,irow,pts_total    , ADJUSTL(TRIM(int2string(npoints_wtr_new)) ))
      CALL set_table_entry(table,irow,pts_diff     , ADJUSTL(TRIM(int2string(npoints_wtr_new-npoints_wtr_old)) ))
      CALL set_table_entry(table,irow,pts_retained , ADJUSTL(TRIM(int2string(npoints_wtr_retained)) ))
      CALL set_table_entry(table,irow,pts_created  , ADJUSTL(TRIM(int2string(npoints_wtr_created)) ))
      CALL set_table_entry(table,irow,pts_destroyed, ADJUSTL(TRIM(int2string(npoints_wtr_destroyed)) ))


      CALL print_table(table, opt_delimiter=' | ')
      CALL finalize_table(table)

      WRITE (0,*) " " ! newline
    ENDIF  ! msg_level

  END SUBROUTINE check_water_idx_lists


  !-------------------------------------------------------------------------
  !!
  !! Update fields that depend on the actual ndvi ratio
  !!
  SUBROUTINE update_ndvi_dependent_fields( p_patch, ext_data, nh_diag )

    TYPE(t_patch)        , INTENT(IN)    :: p_patch
    TYPE(t_external_data), INTENT(INOUT) :: ext_data
    TYPE(t_nh_diag),      INTENT(IN)    :: nh_diag

    ! local
    INTEGER  :: jb,jt,ic,jc, jt_in
    INTEGER  :: rl_start, rl_end
    INTEGER  :: i_startblk, i_endblk    !> blocks
    INTEGER  :: i_count
    INTEGER  :: lu_subs

    !-------------------------------------------------------------------------


     ! exclude the boundary interpolation zone of nested domains
     rl_start = grf_bdywidth_c+1
     rl_end   = min_rlcell_int

     i_startblk = p_patch%cells%start_block(rl_start)
     i_endblk   = p_patch%cells%end_block(rl_end)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jt,ic,i_count,jc,lu_subs,jt_in)
     DO jb = i_startblk, i_endblk
       IF (ext_data%atm%list_land%ncount(jb) == 0) CYCLE ! skip loop if there is no land point
       IF (ntiles_lnd == 1) THEN
         i_count = ext_data%atm%lp_count_t(jb,1)
!$NEC ivdep
         DO ic = 1, i_count
           jc = ext_data%atm%idx_lst_lp_t(ic,jb,1)
           ! plant cover
           ext_data%atm%plcov_t  (jc,jb,1)  = ext_data%atm%ndviratio(jc,jb)                                      &
             &     * MIN(ext_data%atm%ndvi_max(jc,jb),ext_data%atm%plcov_mx(jc,jb))
           ! transpiration area index
           ext_data%atm%tai_t    (jc,jb,1)  = ext_data%atm%plcov_t  (jc,jb,1)                                    &
             &                                  * ext_data%atm%lai_mx(jc,jb)
           ! surface area index
           IF (lterra_urb) THEN
             ext_data%atm%sai_t  (jc,jb,1)  = c_lnd * (1.0_wp - ext_data%atm%urb_isa_t(jc,jb,1))                 &
                                            + ext_data%atm%urb_ai_t(jc,jb,1) * ext_data%atm%urb_isa_t(jc,jb,1)
           ELSE
             ext_data%atm%sai_t  (jc,jb,1)  = c_lnd
           END IF
           ext_data%atm%sai_t(jc,jb,1) = ext_data%atm%sai_t(jc,jb,1) + ext_data%atm%tai_t  (jc,jb,1)             &
                                                                 + c_stm*ext_data%atm%plcov_t(jc,jb,1)

         END DO
       ELSE ! ntiles_lnd > 1
         DO jt=1,ntiles_lnd
           i_count = ext_data%atm%lp_count_t(jb,jt)
!$NEC ivdep
           DO ic = 1, i_count
             jc = ext_data%atm%idx_lst_lp_t(ic,jb,jt)

             ! fix for non-dominant land points
             !!! only active for 'old' extpar datasets (20131009 and earlier) !!!
             IF (ext_data%atm%fr_land(jc,jb) < 0.5_wp) THEN
               IF (ext_data%atm%ndviratio(jc,jb) <= 0.0_wp) THEN  ! here: 0=extpar_missval
                 ! ... reset ndviratio to 0.5
                 ext_data%atm%ndviratio(jc,jb) = 0.5_wp
               ENDIF
             ENDIF

             lu_subs = ext_data%atm%lc_class_t(jc,jb,jt)
             IF (lu_subs < 0) CYCLE

             ! plant cover
             ext_data%atm%plcov_t(jc,jb,jt) = ext_data%atm%ndviratio(jc,jb)                                      &
               & * MIN(ext_data%atm%ndvi_max(jc,jb),ext_data%atm%plcovmax_lcc(lu_subs))
             ! transpiration area index
             ext_data%atm%tai_t  (jc,jb,jt) = ext_data%atm%plcov_t(jc,jb,jt)                                     &
               & * ext_data%atm%laimax_lcc(lu_subs)
             ! surface area index
             IF (lterra_urb) THEN
               ext_data%atm%sai_t(jc,jb,jt) = c_lnd * (1.0_wp - ext_data%atm%urb_isa_t(jc,jb,jt))                &
                                            + ext_data%atm%urb_ai_t(jc,jb,jt) * ext_data%atm%urb_isa_t(jc,jb,jt)
             ELSE
               ext_data%atm%sai_t(jc,jb,jt) = c_lnd
             END IF
             ext_data%atm%sai_t  (jc,jb,jt) = ext_data%atm%sai_t(jc,jb,jt) + ext_data%atm%tai_t  (jc,jb,jt)      &
                                                                     + c_stm*ext_data%atm%plcov_t(jc,jb,jt)

           END DO !ic
         END DO !jt

       END IF !ntiles
!

       IF (lsnowtile) THEN ! copy static external data fields to snow tile grid points
         DO jt = ntiles_lnd+1, ntiles_total

           jt_in = jt - ntiles_lnd
!$NEC ivdep
           DO ic = 1, ext_data%atm%lp_count_t(jb,jt)
             jc = ext_data%atm%idx_lst_lp_t(ic,jb,jt)
             ext_data%atm%plcov_t(jc,jb,jt)    = ext_data%atm%plcov_t(jc,jb,jt_in)
             ext_data%atm%tai_t(jc,jb,jt)      = ext_data%atm%tai_t(jc,jb,jt_in)
             ext_data%atm%sai_t(jc,jb,jt)      = ext_data%atm%sai_t(jc,jb,jt_in)
           ENDDO !ic

         ENDDO !jt
       ENDIF !lsnowtile

     ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL


    IF (itype_vegetation_cycle >= 2) THEN
      CALL vege_clim (p_patch, ext_data, nh_diag)
    ENDIF

    CALL diagnose_ext_aggr (p_patch, ext_data)


  END SUBROUTINE update_ndvi_dependent_fields

!-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !! Copies the tile-based prognostic land-state variables from time level now to 
  !! time level new. This has no relevance for the forecast results but
  !! avoids missing values when writing output at an odd integer multiple of the global
  !! physics time step
  !!
  !-------------------------------------------------------------------------
 
  SUBROUTINE copy_lnd_prog_now2new(p_patch, p_prog_lnd_now, p_prog_lnd_new)


    TYPE(t_patch),         INTENT(IN)    :: p_patch       !<grid/patch info.
    TYPE(t_lnd_prog)     , INTENT(INOUT) :: p_prog_lnd_now
    TYPE(t_lnd_prog)     , INTENT(INOUT) :: p_prog_lnd_new
    
    ! Local array bounds:
    
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: is, ie    !< slices

    ! Local :

    INTEGER :: jb,jt
!-------------------------------------------------------------------------


    ! include nest boundary and halo points
    rl_start = 1
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jt,is,ie)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & is, ie, rl_start, rl_end)

      DO jt = 1, ntiles_total+ntiles_water
        p_prog_lnd_new%t_s_t(is:ie,jb,jt) = p_prog_lnd_now%t_s_t(is:ie,jb,jt)
        p_prog_lnd_new%t_sk_t(is:ie,jb,jt) = p_prog_lnd_now%t_sk_t(is:ie,jb,jt)
      ENDDO

      DO jt = 1, ntiles_total
        p_prog_lnd_new%w_i_t(is:ie,jb,jt)        = p_prog_lnd_now%w_i_t(is:ie,jb,jt)
        p_prog_lnd_new%t_so_t(is:ie,:,jb,jt)     = p_prog_lnd_now%t_so_t(is:ie,:,jb,jt)
        p_prog_lnd_new%w_so_t(is:ie,:,jb,jt)     = p_prog_lnd_now%w_so_t(is:ie,:,jb,jt)
        p_prog_lnd_new%w_so_ice_t(is:ie,:,jb,jt) = p_prog_lnd_now%w_so_ice_t(is:ie,:,jb,jt)
        p_prog_lnd_new%t_snow_t(is:ie,jb,jt)     = p_prog_lnd_now%t_snow_t(is:ie,jb,jt)
        p_prog_lnd_new%w_snow_t(is:ie,jb,jt)     = p_prog_lnd_now%w_snow_t(is:ie,jb,jt)
        p_prog_lnd_new%rho_snow_t(is:ie,jb,jt)   = p_prog_lnd_now%rho_snow_t(is:ie,jb,jt)
      ENDDO

      IF (l2lay_rho_snow .OR. lmulti_snow) THEN
        DO jt = 1, ntiles_total
          p_prog_lnd_new%rho_snow_mult_t(is:ie,:,jb,jt) = p_prog_lnd_now%rho_snow_mult_t(is:ie,:,jb,jt)
        ENDDO
      ENDIF

      IF (lmulti_snow) THEN
        DO jt = 1, ntiles_total
          p_prog_lnd_new%t_snow_mult_t(is:ie,:,jb,jt)   = p_prog_lnd_now%t_snow_mult_t(is:ie,:,jb,jt)
          p_prog_lnd_new%wliq_snow_t(is:ie,:,jb,jt)     = p_prog_lnd_now%wliq_snow_t(is:ie,:,jb,jt)
          p_prog_lnd_new%wtot_snow_t(is:ie,:,jb,jt)     = p_prog_lnd_now%wtot_snow_t(is:ie,:,jb,jt)
          p_prog_lnd_new%dzh_snow_t(is:ie,:,jb,jt)      = p_prog_lnd_now%dzh_snow_t(is:ie,:,jb,jt)
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE copy_lnd_prog_now2new


  !! Coldstart initialization of prognostic seaice albedo alb_si
  !!
  SUBROUTINE seaice_albedo_coldstart(p_patch, p_lnd_state, ext_data)

    TYPE(t_patch),             INTENT(IN)    :: p_patch
    TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state
    TYPE(t_external_data),     INTENT(INOUT) :: ext_data

    ! local
    INTEGER :: jb, jc, ic          ! loop counter
    INTEGER :: jg
    INTEGER :: i_rlstart, i_rlend
    INTEGER :: i_startblk, i_endblk
    !
    REAL(wp):: zfrice_thrhld
    REAL(wp):: frsi(nproma)
    REAL(wp):: t_ice_now(nproma)
    REAL(wp):: alb_si_now(nproma), alb_si_new(nproma)
    !
    TYPE(t_lnd_diag),  POINTER :: lnd_diag
    TYPE(t_wtr_prog),  POINTER :: wtr_prog_now, wtr_prog_new
  !------------------------------------------------------------

    jg = p_patch%id

    WRITE(message_text,'(a,i3)') 'Coldstart initialization of prognostic seaice albedo in domain ',jg
    CALL message('', TRIM(message_text))

    ! set frice_thrhld depending on tile usage
    IF ( ntiles_total == 1 .AND. atm_phy_nwp_config(jg)%inwp_surface /= LSS_JSBACH ) THEN  ! no tile approach
      zfrice_thrhld = 0.5_wp
    ELSE
      zfrice_thrhld = frsi_min
    ENDIF

    lnd_diag     => p_lnd_state%diag_lnd
    wtr_prog_now => p_lnd_state%prog_wtr(nnow_rcf(jg))
    wtr_prog_new => p_lnd_state%prog_wtr(nnew_rcf(jg))

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,ic,frsi,t_ice_now,alb_si_now,alb_si_new)
    DO jb = i_startblk, i_endblk
 
      DO ic = 1, ext_data%atm%list_sea%ncount(jb)
        jc = ext_data%atm%list_sea%idx(ic,jb)
        !
        frsi(ic)      = lnd_diag%fr_seaice(jc,jb)
        t_ice_now(ic) = wtr_prog_now%t_ice(jc,jb)
        ! the following 2 lines are required, because 
        ! only a subset of points in alb_si_now(1:list_sea%ncount) 
        ! are filled by the following initialization routine 
        ! (i.e. sea-ice points only).
        alb_si_now(ic)= wtr_prog_now%alb_si(jc,jb)
        alb_si_new(ic)= wtr_prog_new%alb_si(jc,jb)
      ENDDO

      CALL  seaice_coldinit_albsi_nwp (                                       &
        &                    nswgb        = ext_data%atm%list_sea%ncount(jb), & !in
        &                    frice_thrhld = zfrice_thrhld,                    & !in
        &                    frsi         = frsi(:),                          & !in
        &                    tice_p       = t_ice_now(:),                     & !in
        &                    albsi_p      = alb_si_now(:),                    & !inout
        &                    albsi_n      = alb_si_new(:)                     & !inout
        &  )

      DO ic = 1, ext_data%atm%list_sea%ncount(jb)
        jc = ext_data%atm%list_sea%idx(ic,jb)
        wtr_prog_now%alb_si(jc,jb) = alb_si_now(ic)
        wtr_prog_new%alb_si(jc,jb) = alb_si_new(ic)
      ENDDO
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE seaice_albedo_coldstart

END MODULE mo_nwp_sfc_utils
