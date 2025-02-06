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

! This module is the interface between nwp_nh_interface to the
! turbulence parameterisations:
! inwp_turb == 1 == turbulence scheme by M. Raschendorfer run in COSMO
! inwp_turb == 2 == turbulence scheme imported from the GME
! This module handles the computation of surface transfer coefficients, only.

!OPTION! -cont -msg o
! this command should fix the problem of copying arrays in a subroutine call

!----------------------------
#include "omp_definitions.inc"
!----------------------------

#if defined __xlC__
@PROCESS SPILL(988)
#endif
MODULE mo_nwp_turbtrans_interface

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: message, finish, message_text
  USE mo_model_domain,         ONLY: t_patch
  USE mo_impl_constants,       ONLY: min_rlcell_int, icosmo, igme, ismag, iprog, max_dom
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_physical_constants,   ONLY: rd_o_cpd, grav, lh_v=>alv, lh_s=>als, rd, cpd, tf_salt
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_nwp_tuning_config,    ONLY: itune_gust_diag, tune_gustlim_agl
  USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_tend
  USE mo_nwp_phy_state,        ONLY: phy_params
  USE mo_nwp_lnd_types,        ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_parallel_config,      ONLY: nproma
  USE mo_run_config,           ONLY: msg_level, iqv, iqc, iqtke
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config, lcuda_graph_turb_tran
  USE mo_nonhydrostatic_config,ONLY: kstart_moist
  USE mo_advection_config,     ONLY: advection_config
  USE turb_data,               ONLY: get_turbdiff_param, &
                                     rsur_sher, imode_suradap, rat_can, c_lnd, imode_snowsmot
  USE mo_initicon_config,      ONLY: icpl_da_sfcfric
  USE sfc_flake_data,          ONLY: h_Ice_min_flk, tpl_T_f
  USE turb_transfer,           ONLY: turbtran
  USE mo_satad,                ONLY: sat_pres_water, spec_humi
  USE mo_gme_turbdiff,         ONLY: parturs, nearsfc
  USE mo_util_phys,            ONLY: nwp_dyn_gust
  USE mo_run_config,           ONLY: ltestcase
  USE mo_lnd_nwp_config,       ONLY: ntiles_total, ntiles_lnd, ntiles_water, llake, frlnd_thrhld, &
    &                                isub_seaice, isub_lake, isub_water, lseaice, lsnowtile
  USE mo_nh_testcases_nml,     ONLY: nh_test_name
  USE mo_grid_config,          ONLY: l_scm_mode
  USE mo_scm_nml,              ONLY: scm_sfc_mom, scm_sfc_temp ,scm_sfc_qv
  USE mo_nh_torus_exp,         ONLY: set_scm_bnd
  USE mo_timer
  USE mo_run_config,           ONLY: timers_level
  USE mo_fortran_tools,        ONLY: set_acc_host_or_device
  USE mo_coupling_config,      ONLY: is_coupled_to_waves

#ifdef ICON_USE_CUDA_GRAPH
  USE mo_acc_device_management,ONLY: accGraph, accBeginCapture, accEndCapture, accGraphLaunch
  USE, INTRINSIC :: iso_c_binding
#endif

  IMPLICIT NONE

  PRIVATE


  PUBLIC  ::  nwp_turbtrans


#ifdef ICON_USE_CUDA_GRAPH
  TYPE(accGraph) :: graphs(max_dom*2)
  TYPE(c_ptr) :: lnd_prog_new_cache(max_dom*2) = C_NULL_PTR
  LOGICAL :: graph_captured
  INTEGER :: cur_graph_id, ig
#endif
  LOGICAL :: multi_queue_processing
  INTEGER :: acc_async_queue = 1

CONTAINS
  !!
  !!-------------------------------------------------------------------------
  !!
SUBROUTINE nwp_turbtrans  ( tcall_turb_jg,                     & !>in
                          & p_patch, p_metrics,                & !>in
                          & ext_data,                          & !>in
                          & p_prog,                            & !>in
                          & p_prog_rcf,                        & !>inout
                          & p_diag ,                           & !>inout
                          & prm_diag,                          & !>inout
                          & prm_nwp_tend,                      & !>inout 
                          & wtr_prog_new,                      & !>in
                          & lnd_prog_new,                      & !>inout
                          & lnd_diag,                          & !>inout
                          & lacc                               ) !>in


  TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch        !!<grid/patch info.
  TYPE(t_external_data),TARGET,INTENT(in)   :: ext_data        !< external data
  TYPE(t_nh_metrics)          ,INTENT(in)   :: p_metrics
  TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog          !<the prog vars
  TYPE(t_nh_prog)             ,INTENT(inout):: p_prog_rcf      !< current time levels
  TYPE(t_nh_diag),      TARGET,INTENT(inout):: p_diag          !<the diag vars
  TYPE(t_nwp_phy_diag),        INTENT(inout):: prm_diag        !< atm phys vars
  TYPE(t_wtr_prog),            INTENT(in)   :: wtr_prog_new    !< prog vars for wtr
  TYPE(t_lnd_prog),     TARGET,INTENT(inout):: lnd_prog_new    !< prog vars for sfc
  TYPE(t_lnd_diag),            INTENT(inout):: lnd_diag        !< diag vars for sfc
  TYPE(t_nwp_phy_tend), TARGET,INTENT(inout):: prm_nwp_tend    !< atm tend vars 
  REAL(wp),                    INTENT(in)   :: tcall_turb_jg   !< time interval for
                                                               !< turbulence
  LOGICAL, OPTIONAL,           INTENT(in)   :: lacc            !< GPU flag
  LOGICAL :: lzacc ! non-optional version of lacc

  CHARACTER(len=*),PARAMETER :: routine = 'mo_nwp_turbtrans_interface:nwp_turbtrans'

  ! Local array bounds

  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk    !< blocks
  INTEGER :: i_startidx, i_endidx    !< slices
  INTEGER :: nzprv                   !< argument to turbtran

  ! Local Indices and switches:

  INTEGER :: jc,jb,jk,jt,jg,ic,it,i_count           !< loop indices
  INTEGER :: jk_gust(nproma)

  INTEGER  :: nlev, nlevp1, nlevcm                  !< number of full, half and canopy levels
  INTEGER  :: lc_class                              !< land-cover class

  LOGICAL  :: l_lake(nproma), l_sice(nproma), &     !< lake-, ice-surface points
   l_land, l_water
  LOGICAL  :: lgz0inp_loc !< FALSE: turbtran updates gz0 at water points
                          !< TRUE : gz0 is provided externally (e.g. by the wave model) and is not updated by turbtran
  LOGICAL  :: ladsshr     !<treatment of additional shear by NTCs or LLDCs activ

  ! Local variables related to surface roughness:

  REAL(wp), DIMENSION(nproma,0:ntiles_total+ntiles_water) :: gz0_eff_t, sai_eff_t
  REAL(wp) :: z0_mod, z0_min, sai_min, fsn_flt

  REAL(wp) :: fr_land_t(nproma), h_ice_t(nproma), &
   area_frac, fact_z0rough

  ! Local fields needed to reorder turbtran input/output fields for tile approach

  ! 1D fields
  REAL(wp), DIMENSION(nproma)   :: pres_sfc_t, l_hori, rlamh_fac,   &
   urb_isa_t, t_g_t, qv_s_t   

  ! 2D half-level fields
  REAL(wp), DIMENSION(nproma,3) :: z_ifc_t

  ! 2D full level fields (no tiles)
  REAL(wp), DIMENSION(nproma,2) :: u_t, v_t, temp_t, qv_t, qc_t, epr_t

  ! 3D full-level fields (tiles)
  REAL(wp), DIMENSION(nproma,3) :: z_tvs, tvs_t !SQRT(2*TKE) turbulence velocity scale [m/s]
  REAL(wp), DIMENSION(nproma,3) :: tkvm_t, tkvh_t, rcld_t

  ! 2D fields (tiles)
  REAL(wp), DIMENSION(nproma,ntiles_total+ntiles_water) :: & !auxilary arrays with tile-dimension for 
   tfm_t, tfh_t, t_2m_t, qv_2m_t, td_2m_t, rh_2m_t           !surface-varibales with (so far) no tile-specific global arrays

  REAL(wp), DIMENSION(nproma) :: &
   tcm_t, tch_t, tfv_t, tvm_t, tvh_t, tkr_t, &
   u_10m_t, v_10m_t, shfl_s_t, lhfl_s_t, qhfl_s_t, umfl_s_t, vmfl_s_t

  REAL(wp) :: rho_s

  INTEGER :: gp_num_t(ntiles_total+ntiles_water)        ! tile-vector with associated number of grid-points

  TYPE tile_info
       INTEGER, POINTER, CONTIGUOUS :: gp_idx(:)        ! index-list of grid-points belonging to any tile
  END TYPE tile_info

  TYPE (tile_info) :: list_t(ntiles_total+ntiles_water) ! tile-vector with pointers to index-lists of associated grid-points

  INTEGER,  POINTER :: ilist(:)                         ! pointer to the index-list of grid-points belonging to any tile

!--------------------------------------------------------------
#ifdef ICON_USE_CUDA_GRAPH
    multi_queue_processing = lcuda_graph_turb_tran
#else
    multi_queue_processing = .FALSE.
#endif

  CALL set_acc_host_or_device(lzacc, lacc)

  IF (msg_level >= 15) CALL message('mo_nwp_turbtrans_interface:', 'turbulence')
  IF (timers_level > 9) CALL timer_start(timer_nwp_turbtrans)

  ! number of vertical levels
  nlev   = p_patch%nlev
  nlevp1 = p_patch%nlevp1

  ! domain
  jg = p_patch%id

#ifdef ICON_USE_CUDA_GRAPH
  IF (lzacc .AND. lcuda_graph_turb_tran) THEN
    cur_graph_id = -1
    DO ig=1,max_dom*2
      IF (C_LOC(lnd_prog_new) == lnd_prog_new_cache(ig)) THEN
        cur_graph_id = ig
        graph_captured = .TRUE.
        EXIT
      END IF
    END DO

    IF (cur_graph_id < 0) THEN
      DO ig=1,max_dom*2
        IF (lnd_prog_new_cache(ig) == C_NULL_PTR) THEN
          cur_graph_id = ig
          lnd_prog_new_cache(ig) = C_LOC(lnd_prog_new)
          graph_captured = .FALSE.
          EXIT
        END IF
      END DO
    END IF

    IF (cur_graph_id < 0) THEN
      CALL finish('mo_nwp_turbtrans_interface: ', 'error trying to allocate CUDA graph')
    END IF

    IF (graph_captured) THEN
      WRITE(message_text,'(a,i2)') 'executing CUDA graph id ', cur_graph_id
      IF (msg_level >= 14) CALL message('mo_nwp_turbtrans_interface: ', message_text)
      CALL accGraphLaunch(graphs(cur_graph_id), 1)
      !$ACC WAIT(1)
      IF (timers_level > 9) CALL timer_stop(timer_nwp_turbtrans)
      RETURN
    ELSE
      WRITE(message_text,'(a,i2)') 'starting to capture CUDA graph, id ', cur_graph_id
      IF (msg_level >= 13) CALL message('mo_nwp_turbtrans_interface: ', message_text)
      CALL accBeginCapture(1)
    END IF
  END IF
#endif

  !$ACC DATA PRESENT(p_patch, p_metrics, ext_data, p_prog, p_prog_rcf, p_diag) &
  !$ACC   PRESENT(prm_diag, prm_nwp_tend, wtr_prog_new, lnd_prog_new, lnd_diag) &
  !$ACC   PRESENT(phy_params, advection_config) &
  !$ACC   CREATE(jk_gust, tfm_t, tfh_t) &
  !$ACC   CREATE(t_2m_t, qv_2m_t, td_2m_t, rh_2m_t) &
  !$ACC   CREATE(z_tvs) &
  !$ACC   CREATE(gz0_eff_t, sai_eff_t) &
  !$ACC   CREATE(gp_num_t, list_t) &
  !$ACC   ASYNC(1) IF(lzacc)

  ! exclude boundary interpolation zone of nested domains
  rl_start = grf_bdywidth_c+1
  rl_end   = min_rlcell_int

  i_startblk = p_patch%cells%start_block(rl_start)
  i_endblk   = p_patch%cells%end_block(rl_end)


  IF ( ANY( (/icosmo/)==atm_phy_nwp_config(jg)%inwp_turb ) ) THEN
     CALL get_turbdiff_param(jg)
  ENDIF

  ! Scaling factor for SSO contribution to roughness length ("Erdmann Heise formula")
  fact_z0rough = 1.e-5_wp*ATAN(phy_params(jg)%mean_charlen/2250._wp)

  ladsshr = (rsur_sher>0._wp) !treatment of additional shear by NTCs or LLDCs active

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jt,jc,jk,ic,it,ilist,i_startidx,i_endidx,i_count, &
!$OMP fr_land_t,urb_isa_t,l_land,l_water,l_lake,l_sice,h_ice_t, &
!$OMP shfl_s_t,lhfl_s_t,qhfl_s_t,umfl_s_t,vmfl_s_t, &
!$OMP nzprv,lc_class,z_tvs,tcm_t,tch_t,tfm_t,tfh_t,tfv_t,tvm_t,tvh_t,tkr_t,l_hori, &
!$OMP z0_mod,z0_min,gz0_eff_t,sai_min,sai_eff_t,fsn_flt, &
!$OMP t_g_t,qv_s_t,t_2m_t,qv_2m_t,td_2m_t,rh_2m_t,u_10m_t,v_10m_t,pres_sfc_t, &
!$OMP u_t,v_t,temp_t,qv_t,qc_t,epr_t,tkvm_t,tkvh_t,rcld_t,tvs_t,z_ifc_t, &
!$OMP area_frac,nlevcm,jk_gust, &
!$OMP gp_num_t,list_t, &
!$OMP rho_s,rlamh_fac,lgz0inp_loc) ICON_OMP_GUIDED_SCHEDULE

  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
      & i_startidx, i_endidx, rl_start, rl_end)

   !-------------------------------------------------------------------------
   !<  turbulent transfer
   !-------------------------------------------------------------------------

   !<  NOTE: since  turbulence is a fast process it is
   !!        allowed to do a sequential updating except for wind speed
   !!        because back-and-forth interpolation would cause too large errors
   !!  (GZ, 2011-08-29): Nevertheless, tendency fields are now passed to turbdiff
   !!        to have them available for extended diagnostic output


    IF (atm_phy_nwp_config(jg)%inwp_surface == 0) THEN

      ! check dry case
      IF( atm_phy_nwp_config(jg)%inwp_satad == 0) THEN
        lnd_diag%qv_s (:,jb) = 0._wp
      ELSE IF ( ANY( (/icosmo,igme/)==atm_phy_nwp_config(jg)%inwp_turb ) ) THEN
        IF ( ltestcase .AND. nh_test_name == 'wk82') THEN

!DR Note that this must be re-checked, once turbtran is called at the very end
!DR of the fast physics part.
!DIR$ IVDEP
         !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
         !$ACC LOOP GANG VECTOR
         DO jc = i_startidx, i_endidx
          lnd_prog_new%t_g(jc,jb) = p_diag%temp(jc,nlev,jb)*  &
                      ((p_diag%pres_sfc(jc,jb))/p_diag%pres(jc,nlev,jb))**rd_o_cpd
          lnd_diag%qv_s (jc,jb) = &
             &         spec_humi(sat_pres_water(lnd_prog_new%t_g(jc,jb)),&
             &                                   p_diag%pres_sfc(jc,jb) )
          lnd_diag%qv_s(jc,jb) = MIN (lnd_diag%qv_s(jc,jb) ,p_prog_rcf%tracer(jc,nlev,jb,iqv))
         END DO
         !$ACC END PARALLEL
        ELSE
         !
         !> adjust humidity at water surface because of changed surface pressure
         !
!DIR$ IVDEP
         !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
         !$ACC LOOP GANG VECTOR
         DO jc = i_startidx, i_endidx
           lnd_diag%qv_s (jc,jb) = &
             &         spec_humi(sat_pres_water(lnd_prog_new%t_g(jc,jb)),&
             &                                   p_diag%pres_sfc(jc,jb) )
         ENDDO
         !$ACC END PARALLEL
        END IF
      ENDIF

    END IF

    ! Specify land-cover-related roughness length over land points:
    ! NOTE:  open water, lake and sea-ice points are set in turbtran

    IF (atm_phy_nwp_config(jg)%itype_z0 >= 2) THEN

      DO jt = 1, ntiles_total !loop over all land (sub-)tiles or for not-tiled grid-points
         it = MERGE( 0, jt, (ntiles_total==1) )

        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR PRIVATE(jc, lc_class, z0_mod, z0_min, sai_min, fsn_flt, area_frac)
        DO ic = 1, ext_data%atm%gp_count_t(jb,jt) !loop over all grid-points of current land (sub-)tile or over all not-tiled land-points
           jc = ext_data%atm%idx_lst_t(ic,jb,jt)

          sai_eff_t(ic,it) = ext_data%atm%sai_t(jc,jb,jt) !effective SAI equals current land-use SAI

          ! Reduction of land-cover related roughness length due to vegetation (according G.Z): 
          lc_class = MAX(1,ext_data%atm%lc_class_t(jc,jb,jt)) !land-class index (MAX-function as to avoid segfaults) 
          z0_min = ext_data%atm%z0_lcc_min(lc_class) !minimal roughness length (meets fully smoothed snow-roughness)
          z0_mod = ext_data%atm%z0_lcc    (lc_class) !roughness length according to land class
          z0_mod = MAX( z0_mod*SQRT( MAX( 0.5_wp, MIN( 1._wp, 1.3333_wp*ext_data%atm%plcov_t(jc,jb,jt) ) ) ), & !reduced value
                                 !limited to about 70% (SQRT(0.5) of the default RL and starts at at 'plcov' below 75%
                        z0_min ) !ensure that z0 does not fall below the minimum allowed value

          ! Save so far updated g*z0 (without snow-, SSO- or possible AHP-corrections)
          prm_diag%gz0_t(jc,jb,jt) = grav*z0_mod

          ! Reduction of R-length (and SAI) due to the presence of snow:
          IF (imode_snowsmot>0) THEN !any surface smoothing by snow required
            IF (imode_snowsmot>=2) THEN !extentension according to M.R.
              IF (imode_snowsmot==2) THEN !full smoothing for z0 and SAI:
                sai_min = c_lnd !only standard SAI of a non-vegetated surface (full smoothing of SAI)
              ELSE !dynamic surface smoothing by snow of both z0 and SAI dependent on snow- and roughness height
                !smoothing factor of the surface due to the snow-cover:
                fsn_flt = lnd_diag%h_snow_t(jc,jb,jt)/( rat_can*sai_eff_t(ic,it)*z0_mod + lnd_diag%h_snow_t(jc,jb,jt) )
                !estimated z0-value for a fully snow-cov. surface approaching z0_mini for large snow-depth:
                z0_min = MAX( 0._wp, z0_mod - z0_min) !z0_luse - z0_mini
                z0_min = z0_mod - z0_min*fsn_flt      !effective z0-value for a fully snow-covered surface

                !estimated sai-value for a fully snow-cov. surface approaching sai_mini for large snow-depth:
                sai_min = MAX( 0._wp, sai_eff_t(ic,it) - c_lnd ) !sai_luse - sai_mini
                sai_min = sai_eff_t(ic,it) - sai_min*fsn_flt     !effective SAI for a fully snow-cov. surface
              END IF
              sai_eff_t(ic,it) = (1._wp-lnd_diag%snowfrac_t(jc,jb,jt))*sai_eff_t(ic,it) + & !SAI including
                                        lnd_diag%snowfrac_t(jc,jb,jt) *sai_min              ! the snow-cover contribution
            END IF
            area_frac=lnd_diag%snowfrac_t(jc,jb,jt)**2 !quadratic weighting-factor (according to G.Z.)
            z0_mod = (1._wp-area_frac)*z0_mod + & !effective z0 including the
                            area_frac *z0_min     ! snow-cover contribution
          END IF

          ! Further modification of tile-specific R-length as a contribution by Adaptive Parameter Tuning (APT):
          IF (icpl_da_sfcfric >= 1) THEN !apply tuning factor for surface friction derived from DA
            z0_mod = MIN( 1.5_wp, prm_diag%sfcfric_fac(jc,jb)*z0_mod )

            ! thus it is only applied to 'gz0_eff_t' and not to '%gz0_t'!
          END IF

          ! Further modification of tile-specific R-length as a contribution by SSO:
          IF (atm_phy_nwp_config(jg)%itype_z0 == 3) THEN !add SSO contribution to effective tile-specific R-length
            z0_mod = z0_mod + MIN( fact_z0rough*ext_data%atm%sso_stdh_raw(jc,jb)**2, 7.5_wp )

            !Note:
            !Obviously, this SSO-contribution of R-lenght is only related to its aerodynamic application in SUB 'turbtran'
            ! and not to any geometric application in 'terra'.
            !Hence, this contribution is only applied to 'gz0_eff_t' and not to '%gz0_t'.
          END IF

          gz0_eff_t(ic,it) = grav*z0_mod !effective aero-dynamic R-length used in 'turbtran' only

        END DO !loop over all grid-points belonging to a land (sub-)tile
        !$ACC END PARALLEL

      END DO !loop over all land-tiles (jt) or for not-tiled grid-points

    ELSE ! (atm_phy_nwp_config(jg)%itype_z0==1): uniform tile-averaged roughness length (without any modification)

      DO jt = 1, ntiles_total !loop over all land (sub-)tiles or for not-tiled grid-points
         it = MERGE( 0, jt, (ntiles_total==1) )
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR PRIVATE(jc)
!$NEC ivdep
        DO ic = 1, ext_data%atm%gp_count_t(jb,jt) !loop over all grid-points belonging to a land (sub-)tile
           jc = ext_data%atm%idx_lst_t(ic,jb,jt)
          sai_eff_t(ic,it) = ext_data%atm%sai(jc,jb) !grid-point average of SAI
          gz0_eff_t(ic,it) =     prm_diag%gz0(jc,jb) !                  and g*z0 

          !Note(MR): 
          !In case of a not-tiled grid-point with a water surface, '%gz0(jc,jb)=%gz0_t(jc,jb,1)' is just the latest update 
          ! that has been calculated in SUB 'turbtran'.
          !Attention(MR):
          !Particularly at costal grid-points, the grid-point average may by very different from the tile-specific value,
          ! which may be in contradiction with other tile-specific parameters!
        END DO
        !$ACC END PARALLEL
      END DO

    ENDIF !blocks according to '%itype_z0'

    IF (ntiles_total==1) THEN !a model run without tiles
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR PRIVATE(jc)
!$NEC ivdep
      DO ic = 1, ext_data%atm%gp_count_t(jb,1) !loop over all grid-points belonging to a land point
         jc = ext_data%atm%idx_lst_t(ic,jb,1)

         ! Shift 'sai_eff_t'- and 'gz0_eff_t'-values, which are so far stored contiguously at tile-index "0" 
         !  for all land points, back to their grid-point indices at tile-index "1":
         sai_eff_t(jc,1) = sai_eff_t(ic,0)
         gz0_eff_t(jc,1) = gz0_eff_t(ic,0)
      END DO
      !$ACC END PARALLEL

      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx

        ! Fill gaps in the 'sai_eff_t'- and 'gz0_eff_t'-vectors, belonging to non-land points, by already present values:
        IF (ext_data%atm%fr_land(jc,jb)<=0.5_wp) THEN !it's an effective non-land surface (water or ice)
          sai_eff_t(jc,1) = ext_data%atm%sai_t(jc,jb,1)
          gz0_eff_t(jc,1) =     prm_diag%gz0_t(jc,jb,1)
        END IF
      END DO
      !$ACC END PARALLEL

      !Note(MR): 
      !This measures are necessary, since SUB 'turbtran' expects grid-point vectors in this case!
    END IF

    !Note(MR): 
    !For non-land (sub-)tiles, the global tile-specific values '%gz0_t(jc,jb,jt)' (being iteratively updates in SUB 'turbtran')
    ! and '%sai_t(jc,jb,jt)' have not been affected; and they are still being applied in the following!

    IF (atm_phy_nwp_config(jg)%inwp_sso > 0) THEN
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        IF (prm_diag%ktop_envel(jc,jb) < nlev) THEN
          jk_gust(jc) = MERGE(prm_diag%ktop_envel(jc,jb), prm_diag%ktop_envel(jc,jb)-1, itune_gust_diag == 2)
        ELSE
          jk_gust(jc) = nlev
        ENDIF
      ENDDO
      !$ACC END PARALLEL
    ELSE
      !$ACC KERNELS ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      jk_gust(:) = nlev
      !$ACC END KERNELS
    ENDIF

    IF ( ltestcase .AND. l_scm_mode .AND. lzacc .AND. &   !lzacc false in init  step
      &  ((scm_sfc_mom >= 1) .OR. (scm_sfc_temp >= 1) .OR. (scm_sfc_qv >= 1)) ) THEN
#ifdef _OPENACC
      CALL finish( TRIM(routine),'set_scm_bnd is not supported with OpenACC.')
#endif
      CALL set_scm_bnd( nvec=nproma, ivstart=i_startidx, ivend=i_endidx,   &
          & u_s          = p_diag%u(:,nlev,jb),                            & !in
          & v_s          = p_diag%v(:,nlev,jb),                            & !in
          & th_b         = p_diag%temp(:,nlev,jb)/p_prog%exner(:,nlev,jb), & !in
          & qv_b         = p_prog_rcf%tracer(:,nlev,jb,iqv),               & !in
          & pres_sfc     = p_diag%pres_sfc(:,jb),                          & !in
          & dz_bs=p_metrics%z_mc(:,nlev,jb)-p_metrics%z_ifc(:,nlevp1,jb),  & !in
          & z0m=prm_diag%gz0(:,jb)/grav,                                   & !in
          !for now z0m is assumed to be equal to z0h - GABLS1
          & z0h=prm_diag%gz0(:,jb)/grav,                                   & !in
          & prm_nwp_tend = prm_nwp_tend,                                   & !in 
          & tvm          = prm_diag%tvm(:,jb),                             & !inout
          & tvh          = prm_diag%tvh(:,jb),                             & !inout
          & shfl_s       = prm_diag%shfl_s(:,jb),                          & !out
          & qhfl_s       = prm_diag%qhfl_s(:,jb),                          & !out
          & lhfl_s       = prm_diag%lhfl_s(:,jb),                          & !out
          & umfl_s       = prm_diag%umfl_s(:,jb),                          & !out
          & vmfl_s       = prm_diag%vmfl_s(:,jb),                          & !out
          & qv_s         = lnd_diag%qv_s(:,jb),                            & !out
          & t_g          = lnd_prog_new%t_g(:,jb) )                          !out
    ENDIF


    SELECT CASE (atm_phy_nwp_config(jg)%inwp_turb)

    CASE(icosmo)

!-------------------------------------------------------------------------
!< COSMO turbulence scheme by M. Raschendorfer
!-------------------------------------------------------------------------
 
      ! note that TKE must be converted to the turbulence velocity scale SQRT(2*TKE)
      ! for turbdiff
      ! INPUT to turbtran is timestep new
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, 3
        DO jc = i_startidx, i_endidx
          z_tvs(jc,jk) = SQRT(2._wp * p_prog_rcf%tke(jc,nlev-2+jk,jb))
        ENDDO
      ENDDO
      !$ACC END PARALLEL

      ! First call of turbtran for all grid points (water points with > 50% water
      ! fraction and tile 1 of the land points)

      !----------------------------
      IF (ntiles_total == 1) THEN ! tile approach not used; use tile-averaged fields from extpar
      !----------------------------

        !$ACC DATA CREATE(l_hori, l_lake, l_sice) ASYNC(1) IF(lzacc)

        !$ACC KERNELS ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        l_hori(i_startidx:i_endidx)=phy_params(jg)%mean_charlen !horizontal grid-scale (should be dependent on location in future!)
        !$ACC END KERNELS

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          l_lake(jc) = (ext_data%atm%depth_lk(jc,jb)>0._wp) !a lake point
          l_sice(jc) = ext_data%atm%fr_land(jc,jb)<=0.5_wp .AND. & ! a non-land point, and:
                       MERGE( MERGE( (wtr_prog_new%h_ice(jc,jb)>=h_Ice_min_flk), & !either more than marginal lake-ice depth from lake-scheme
                                     (lnd_prog_new%t_g  (jc,jb)< tpl_T_f      ), & ! or frozen lake due to surface-temperature,
                                     llake   ), & !dependent on whether lake scheme is active or not
                              MERGE( (wtr_prog_new%h_ice(jc,jb)>0._wp  ), & !either sea-ice present from sea-ice scheme
                                     (lnd_prog_new%t_g  (jc,jb)<tf_salt), & !or frozen salty sea due to surface temperature
                                     lseaice ), & !dependent on whether sea-ice scheme is active or not
                              l_lake(jc) ) !seperate treatment of lake points and non-lake points
        END DO
        !$ACC END PARALLEL

        !Note(MR):
        !Since there should be no '%h_ice' at non-land points, 'l_sice' indicates all grid-points coverd by sea- or lake-ice,
        ! while frozen land points e.g. covered by snow or glaciers) are excluded.
        !This is additionally ensured by only considering points with "%fr_land(jc,jb)<=0.5p".

        IF (ladsshr) THEN !treatment of additional shear by NTCs or LLDCs active
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx, i_endidx
            prm_diag%tfv_t(jc,jb,1) = prm_diag%tfv(jc,jb)
          ENDDO
          !$ACC END PARALLEL
        END IF

        nzprv  = 1
        nlevcm = 3 !so far no vertically resolved roughness layer (kcm=ke1)

        ! turbtran
        CALL turbtran (               & ! only surface-layer turbulence
!
          &  iini=0,                  & !
          &  ltkeinp=.FALSE.,         & !
          &  lgz0inp=.FALSE.,         &
          &  lstfnct=.TRUE. ,         & ! with stability function
          &  lsrflux=.TRUE. ,         & !
          &  lnsfdia=.TRUE. ,         & ! including near-surface diagnostics
          &  lrunscm=.FALSE.,         & ! no single column run
          &  ladsshr=ladsshr,         & !treatment of additional shear by NTCs or LLDCs for surface layer active
!
          &  dt_tke=tcall_turb_jg,                                                     & !in
          &  nprv=nzprv, ntur=1, ntim=1,                                               & !in ('tke' without time-dimension)
          &  nvec=nproma, ke=2,    ke1=3,      kcm=nlevcm, iblock=jb,                  & !in
          &  ivstart=i_startidx, ivend=i_endidx,                                       & !in
!
          &  l_pat=ext_data%atm%l_pat(:,jb),                                           & !in
          &  l_hori=l_hori,                                                            & !in
          &  hhl=p_metrics%z_ifc(:,nlev-1:nlevp1,jb),                                  & !in
          &  fr_land=ext_data%atm%fr_land(:,jb),                                       & !in
          &  l_lake=l_lake(:),                                                         & !in (lake surfaces)
          &  l_sice=l_sice(:),                                                         & !in (ice  surfaces)
!
          &  rlamh_fac=prm_diag%rlamh_fac_t(:,jb,1),                                   & !in
          &  gz0=gz0_eff_t(:,1),                                                       & !inout (incl. all above modificat.)
          &  sai=sai_eff_t(:,1),                                                       & !in    (incl. all above modificat.)
          &  urb_isa=ext_data%atm%urb_isa_t(:,jb,1),                                   & !in
!
          &  t_g=lnd_prog_new%t_g(:,jb),                                               & !in
          &  qv_s=lnd_diag%qv_s(:,jb),                                                 & !in
          &  ps=p_diag%pres_sfc(:,jb),                                                 & !in
!
          &  u=p_diag%u(:,nlev-1:nlev,jb),                                             & !in
          &  v=p_diag%v(:,nlev-1:nlev,jb),                                             & !in
          &  t=p_diag%temp(:,nlev-1:nlev,jb),                                          & !in
          &  qv=p_prog_rcf%tracer(:,nlev-1:nlev,jb,iqv),                               & !in
          &  qc=p_prog_rcf%tracer(:,nlev-1:nlev,jb,iqc),                               & !in
          &  epr=p_prog%exner(:,nlev-1:nlev,jb),                                       & !in
!
          &  tcm=prm_diag%tcm_t(:,jb,1),                                               & !out
          &  tch=prm_diag%tch_t(:,jb,1),                                               & !out
          &  tvm=prm_diag%tvm_t(:,jb,1),                                               & !inout
          &  tvh=prm_diag%tvh_t(:,jb,1),                                               & !inout
          &  tfm=prm_diag%tfm(:,jb),                                                   & !inout
          &  tfh=prm_diag%tfh(:,jb),                                                   & !inout
          &  tfv=prm_diag%tfv_t(:,jb,1),                                               & !inout
          &  tkr=prm_diag%tkr_t(:,jb,1),                                               & !inout
          &  tke=z_tvs(:,:),                                                           & !inout
          &  tkvm=prm_diag%tkvm(:,nlev-1:nlevp1,jb),                                   & !inout
          &  tkvh=prm_diag%tkvh(:,nlev-1:nlevp1,jb),                                   & !inout
          &  rcld=prm_diag%rcld(:,nlev-1:nlevp1,jb),                                   & !inout

          &  t_2m=prm_diag%t_2m(:,jb),                                                 & !inout
          &  qv_2m=prm_diag%qv_2m(:,jb),                                               & !out
          &  td_2m=prm_diag%td_2m(:,jb),                                               & !out
          &  rh_2m=prm_diag%rh_2m(:,jb),                                               & !out
          &  u_10m=prm_diag%u_10m_t(:,jb,1),                                           & !out
          &  v_10m=prm_diag%v_10m_t(:,jb,1),                                           & !out

          &  shfl_s=prm_diag%shfl_s_t(:,jb,1),                                         & !out
          &  qvfl_s=prm_diag%qhfl_s_t(:,jb,1),                                         & !out
          &  umfl_s=prm_diag%umfl_s_t(:,jb,1),                                         & !out
          &  vmfl_s=prm_diag%vmfl_s_t(:,jb,1),                                         & !out
!
          &  lacc=lzacc,                                                               & !in
          &  opt_acc_async_queue=1                                                     & !in
          &                                                                            ) !end of 'turbtran' call

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          prm_diag%lhfl_s_t(jc,jb,1) = prm_diag%qhfl_s_t(jc,jb,1) &
                                      *MERGE( lh_s, & !latent heat-flux with respect to ice surface
                                              lh_v, & !                                 water
                                              l_sice(jc) ) !a frozen surface of a lake or sea-ice point
          !Notes:
          !So far, for a frozen land-surface, always 'lh_v' is used. However, '%qhfl_s_t' and '%lhfl_s_t',
          ! as well as '%shfl_s_t', are going to be overwritten for land-points by 'terra';
          ! hence, this calculation only matters for non-land points (lakes, sea-water or sea-ice).
          !Further, for atmospheric vertical diffusion, only the grid-point variables of the fluxes
          ! '%shfl_s, %qhfl_s, umfl_s, vmfl_s' are used, which are loaded by the values for tile "1" in 'mo_nwp_sfc_interface'.
          !While, at land points, '%lhfl_s' (as well as '%lhfl_s_t') is only used for model-output,
          ! '%lhfl_s_t' (as well as '%shfl_t') is used as input for the lake- and seaice-schemes.
          IF (ext_data%atm%fr_land(jc,jb)<=0.5_wp) THEN !for ice- or water-surfaces only
            prm_diag%gz0_t(jc,jb,1) = gz0_eff_t(jc,1) !updated effective value for next call of 'turbran'
          END IF

          ! copy
          prm_diag%gz0(jc,jb)   = gz0_eff_t     (jc,   1) !for TURBDIFF

          prm_diag%tcm(jc,jb)   = prm_diag%tcm_t(jc,jb,1) !only for output
          prm_diag%tch(jc,jb)   = prm_diag%tch_t(jc,jb,1) !only for output
          prm_diag%tvm(jc,jb)   = prm_diag%tvm_t(jc,jb,1) !for TURBDIFF
          prm_diag%tvh(jc,jb)   = prm_diag%tvh_t(jc,jb,1) !for TURBDIFF
          prm_diag%u_10m(jc,jb) = prm_diag%u_10m_t(jc,jb,1) !only for output
          prm_diag%v_10m(jc,jb) = prm_diag%v_10m_t(jc,jb,1) !only for output

          ! instantaneous max/min 2m temperature over tiles (trivial operation for 1 tile)
          prm_diag%t_tilemax_inst_2m(jc,jb) = prm_diag%t_2m(jc,jb)
          prm_diag%t_tilemin_inst_2m(jc,jb) = prm_diag%t_2m(jc,jb)
          prm_diag%tmax_2m(jc,jb) = MAX(prm_diag%t_2m(jc,jb), &
            &                                        prm_diag%tmax_2m(jc,jb) )
          prm_diag%tmin_2m(jc,jb) = MIN(prm_diag%t_2m(jc,jb), &
            &                                        prm_diag%tmin_2m(jc,jb) )
        ENDDO
        !$ACC END PARALLEL

        !$ACC END DATA

      !----------------------------
      ELSE ! tile approach used
      !----------------------------

        ! Vector of pointers to the index lists of grid-points belonging to each tile:
        DO jt = 1, ntiles_total ! for all land tiles
          list_t(jt)%gp_idx => ext_data%atm%idx_lst_t(:,jb,jt)
        END DO
        list_t(isub_water)%gp_idx  => ext_data%atm%list_seawtr%idx(:,jb) ! for the open-sea tile
        list_t(isub_lake)%gp_idx   => ext_data%atm%list_lake%idx  (:,jb) ! for the lake tile
        list_t(isub_seaice)%gp_idx => ext_data%atm%list_seaice%idx(:,jb) ! for the sea-ice tile

        ! Vector of grid-point counters belonging to each tile:
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO jt = 1, ntiles_total+ntiles_water
          IF (jt <= ntiles_total) THEN    ! for land tiles
            gp_num_t(jt) = ext_data%atm%gp_count_t(jb,jt)
          ELSEIF (jt == isub_water) THEN  ! for the open-sea tile
            gp_num_t(jt) = ext_data%atm%list_seawtr%ncount(jb)
          ELSEIF (jt == isub_lake) THEN   ! for the lake tile
            gp_num_t(jt) = ext_data%atm%list_lake%ncount(jb)
          ELSE !IF (jt == isub_seaice) THEN ! for the sea-ice tile
            gp_num_t(jt) = ext_data%atm%list_seaice%ncount(jb)
            !Note: So far, there is no further tile defined!
          END IF
        END DO
        !$ACC END PARALLEL

        IF (multi_queue_processing) THEN
          ! Fork into concurrent GPU-streams for each tile:
          DO jt = 1, ntiles_total + ntiles_water
            !$ACC WAIT(1) ASYNC(jt)
          END DO
        END IF

#ifndef _OPENACC
        l_land  = .FALSE. !no land  tile called yet
        l_water = .FALSE. !no water tile called yet
#endif

        ! Working loop over land tile points, sea, lake points and seaice points
        !  with a separate index list for each tile:
        !-----------------------------
        DO jt = 1, ntiles_total + ntiles_water
        !-----------------------------

          IF (is_coupled_to_waves().AND.(jt==isub_water)) THEN
            lgz0inp_loc = .TRUE.  ! gz0 at non ice-covered sea water points is provided
                                  ! from external sources (e.g. ICON-waves). No update by turbtran.
          ELSE
            lgz0inp_loc = .FALSE. ! gz0 at water points is updated by turbtran
          END IF

          IF (multi_queue_processing) acc_async_queue = jt
          !$ACC DATA CREATE(u_t, v_t, temp_t, qv_t, qc_t, epr_t, z_ifc_t, pres_sfc_t) &
          !$ACC   CREATE(l_hori, fr_land_t, l_lake, l_sice, h_ice_t, rlamh_fac) &
          !$ACC   CREATE(urb_isa_t, t_g_t, qv_s_t, i_count) &
          !$ACC   CREATE(tcm_t, tch_t, tfv_t, tvm_t, tvh_t, tkr_t, tkvm_t, tkvh_t, rcld_t, tvs_t) &
          !$ACC   CREATE(u_10m_t, v_10m_t, shfl_s_t, lhfl_s_t, qhfl_s_t, umfl_s_t, vmfl_s_t) &
          !$ACC   ASYNC(acc_async_queue) IF(lzacc)

          !$ACC KERNELS DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
          i_count = gp_num_t(jt)
          !$ACC END KERNELS
          ilist => list_t(jt)%gp_idx
          !Note(MR): 
          !Since the scalar 'i_count' is not a constant, it needs special treatment for 'cuda_graph' capturing.

#ifndef _OPENACC
          IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty
#endif

          ! Copy input fields to the local re-indexed variables:
          ! It remains to be determined which of the model levels are actually needed for non-init calls.

          !$ACC PARALLEL ASYNC(acc_async_queue) DEFAULT(PRESENT) IF(lzacc)
          !$ACC LOOP GANG VECTOR PRIVATE(jc)
!$NEC ivdep
          DO ic = 1, i_count
            jc = ilist(ic)

            z_ifc_t (ic,1:3)    = p_metrics%z_ifc   (jc,nlev-1:nlevp1,jb)
            u_t    (ic,1:2)     = p_diag%u          (jc,nlev-1:nlev  ,jb)
            v_t    (ic,1:2)     = p_diag%v          (jc,nlev-1:nlev  ,jb)
            temp_t (ic,1:2)     = p_diag%temp       (jc,nlev-1:nlev  ,jb)
            epr_t  (ic,1:2)     = p_prog%exner      (jc,nlev-1:nlev  ,jb)
            qv_t   (ic,1:2)     = p_prog_rcf%tracer (jc,nlev-1:nlev  ,jb,iqv)
            qc_t   (ic,1:2)     = p_prog_rcf%tracer (jc,nlev-1:nlev  ,jb,iqc)
            pres_sfc_t(ic)      = p_diag%pres_sfc   (jc,jb)
            IF (jt>ntiles_total) THEN !only for non-land (sub-)tiles
              gz0_eff_t(ic,jt)  = prm_diag%gz0_t    (jc,jb,jt)     ! effective value equals previous global value
              sai_eff_t(ic,jt)  = ext_data%atm%sai_t(jc,jb,jt)     ! effective value equals previous global value
            END IF

            t_g_t  (ic)         = lnd_prog_new%t_g_t(jc,jb,jt)
            qv_s_t (ic)         = lnd_diag%qv_s_t   (jc,jb,jt)
            rcld_t (ic,1:2)     = prm_diag%rcld     (jc,nlev-1:nlev,jb)
            rcld_t (ic,3)       = prm_diag%rcld_s_t (jc,jb,jt)     ! tile-specific for lowest level (to be activated)
            tvs_t  (ic,1:2)     = z_tvs             (jc,1:2)
            tvs_t  (ic,3)       = prm_diag%tvs_s_t  (jc,jb,jt)     ! tile-specific for lowest level
            urb_isa_t(ic)       = ext_data%atm%urb_isa_t(jc,jb,jt)
            tkvm_t (ic,1:2)     = prm_diag%tkvm     (jc,nlev-1:nlev,jb)
            tkvm_t (ic,3)       = prm_diag%tkvm_s_t (jc,jb,jt)     ! tile-specific for lowest level
            tkvh_t (ic,1:2)     = prm_diag%tkvh     (jc,nlev-1:nlev,jb)
            tkvh_t (ic,3)       = prm_diag%tkvh_s_t (jc,jb,jt)     ! tile-specific for lowest level
            tkr_t  (ic)         = prm_diag%tkr_t    (jc,jb,jt)     ! used for time-step iteration (if "imode_trancnf>=4")
            IF (imode_suradap>0) THEN !artific. amplification of 'tkvm(:,ke)' required for surface layer
              tfm_t(ic,jt)      = prm_diag%tfm      (jc,jb)        ! reduct.-fact. for 'tkvm(:,ke)' due to LLDCs
              tfh_t(ic,jt)      = prm_diag%tfh      (jc,jb)        ! reduct.-fact. for 'tkvh(:,ke)' due to LLDCs
            END IF
            IF (ladsshr) THEN !treatment of additional shear by NTCs or LLDCs for surface layer active
              tfv_t(ic)         = prm_diag%tfv      (jc,jb)        ! wind-shear amplificat.-factor  due to NTCs
            END IF

            tvm_t  (ic)         = prm_diag%tvm_t    (jc,jb,jt)
            tvh_t  (ic)         = prm_diag%tvh_t    (jc,jb,jt)
            rlamh_fac(ic)       = prm_diag%rlamh_fac_t(jc,jb,jt)

            IF ((jt == ntiles_total+2) .AND. llake) THEN  !only for a lake tile with active lake scheme
              h_ice_t(ic)       = wtr_prog_new%h_ice(jc,jb)
            END IF

            l_hori(ic)=phy_params(jg)%mean_charlen !horizontal grid-scale (should be dependent on location in future!)

          ENDDO
          !$ACC END PARALLEL

          ! Tile-specific assignment of 'fr_land_t, l_lake, l_sice':

#ifndef _OPENACC
          IF (.NOT.l_land .AND. (jt <= ntiles_total)) THEN !the first present land tile
            l_land = .TRUE.
#else
          IF (jt <= ntiles_total) THEN ! land tile points
#endif
            !$ACC KERNELS DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
            fr_land_t(:)  = 1._wp
            l_lake(:) = .FALSE.
            l_sice(:) = .FALSE.
            !$ACC END KERNELS
#ifndef _OPENACC
          ELSEIF (.NOT.l_water .AND. (jt > ntiles_total)) THEN !the first present non-land/water tile
            l_water = .TRUE.
#else
          ELSE !non-land points
#endif
            !$ACC KERNELS ASYNC(acc_async_queue) DEFAULT(PRESENT) IF(lzacc)
            fr_land_t(:) = 0._wp
            !$ACC END KERNELS
          END IF

          !$ACC PARALLEL ASYNC(acc_async_queue) DEFAULT(PRESENT) IF(lzacc)
          IF (jt == ntiles_total + 1) THEN ! sea points (open water)
            !$ACC LOOP GANG VECTOR
            DO ic= 1, i_count
              l_lake(ic) = .FALSE.
              l_sice(ic) = MERGE( .FALSE., (t_g_t(ic)<tpl_T_f), lseaice ) !temperature dependent sea-ice, if "lseaice=F"
            END DO
          ELSE IF (jt == ntiles_total + 2) THEN !lake points
            !$ACC LOOP GANG VECTOR
            DO ic= 1, i_count
              l_lake(ic) = .TRUE.
              l_sice(ic) = MERGE( (h_ice_t(ic)>=h_Ice_min_flk), (t_g_t(ic)< tpl_T_f), llake ) !dependent on 'llake':
                                  !either more than marginal lake-ice depth from lake-scheme, or frozen lake due to surface-temperature
            END DO
          ELSE IF (jt == ntiles_total + 3) THEN !seaice points (only present, if "lseaice=T")
            !$ACC LOOP GANG VECTOR
            DO ic= 1, i_count
              l_lake(ic) = .FALSE.
              l_sice(ic) = .TRUE.
            END DO
          ENDIF
          !$ACC END PARALLEL
          !Notes (MR):
          !'l_sice' indicates all grid-points coverd by sea- or lake-ice, while frozen land points
          ! (e.g. covered by snow or glaciers) are excluded.

          nlevcm = 3 !so far no vertically resolved roughness layer (kcm=ke1)
          nzprv  = 1

          ! turbtran
          CALL turbtran (               & ! only surface-layer turbulence
!
            &  iini=0,                  & !
            &  ltkeinp=.FALSE.,         & !
            &  lgz0inp=lgz0inp_loc,     & !
            &  lstfnct=.TRUE. ,         & ! with stability function
            &  lsrflux=.TRUE. ,         & !
            &  lnsfdia=.TRUE. ,         & ! including near-surface diagnostics
            &  lrunscm=.FALSE.,         & ! no single column run
            &  ladsshr=ladsshr,         & !treatment of additional shear by NTCs or LLDCs for surface layer active
!
            &  dt_tke=tcall_turb_jg,                                        & !in
            &  nprv=nzprv, ntur=1, ntim=1,                                  & !in
            &  nvec=nproma, ke=2, ke1=3, kcm=nlevcm, iblock=jb,             & !in
            &  ivstart=1, ivend=i_count,                                    & !in
!
            &  l_pat=ext_data%atm%l_pat(:,jb),                              & !in
            &  l_hori=l_hori,                                               & !in
            &  hhl=z_ifc_t(:,:),                                            & !in
            &  fr_land=fr_land_t(:),                                        & !in
            &  l_lake=l_lake(:),                                            & !in (lake surfaces)
            &  l_sice=l_sice(:),                                            & !in (ice  surfaces)
            &  rlamh_fac=rlamh_fac(:),                                      & !in
            &  gz0=gz0_eff_t(:,jt),                                         & !inout effective value including all modifictions
            &  sai=sai_eff_t(:,jt),                                         & !in    effective value including all modifictions
            &  urb_isa=urb_isa_t(:),                                        & !in
!
            &  t_g=t_g_t(:),                                                & !in
            &  qv_s=qv_s_t(:),                                              & !in
            &  ps=pres_sfc_t(:),                                            & !in
!
            &  u=u_t(:,:),                                                  & !in
            &  v=v_t(:,:),                                                  & !in
            &  t=temp_t(:,:),                                               & !in
            &  qv=qv_t(:,:),                                                & !in
            &  qc=qc_t(:,:),                                                & !in
            &  epr=epr_t(:,:),                                              & !in
!
            &  tcm=tcm_t(:),                                                & !out
            &  tch=tch_t(:),                                                & !out
            &  tvm=tvm_t(:),                                                & !inout
            &  tvh=tvh_t(:),                                                & !inout
            &  tfm=tfm_t(:,jt),                                             & !inout
            &  tfh=tfh_t(:,jt),                                             & !inout
            &  tfv=tfv_t(:),                                                & !inout
            &  tkr=tkr_t(:),                                                & !inout
            &  tke=tvs_t(:,:),                                              & !inout
            &  tkvm=tkvm_t(:,:),                                            & !inout
            &  tkvh=tkvh_t(:,:),                                            & !inout
            &  rcld=rcld_t(:,:),                                            & !inout
!
            &   t_2m= t_2m_t(:,jt),                                         & !out
            &  qv_2m=qv_2m_t(:,jt),                                         & !out
            &  td_2m=td_2m_t(:,jt),                                         & !out
            &  rh_2m=rh_2m_t(:,jt),                                         & !out
            &  u_10m=u_10m_t(:),                                            & !out
            &  v_10m=v_10m_t(:),                                            & !out
!
            &  shfl_s=shfl_s_t(:),                                          & !out
            &  qvfl_s=qhfl_s_t(:),                                          & !out
            &  umfl_s=umfl_s_t(:),                                          & !out
            &  vmfl_s=vmfl_s_t(:),                                          & !out
!
            &  lacc=lzacc,                                                  & !in
            &  opt_acc_async_queue=acc_async_queue                          &
            &                                                               ) !end of 'turbtran' call

          !$ACC PARALLEL ASYNC(acc_async_queue) DEFAULT(PRESENT) IF(lzacc)
          !$ACC LOOP GANG VECTOR
          DO ic= 1, i_count
            lhfl_s_t(ic) = qhfl_s_t(ic)*MERGE( lh_s, & !latent heat-flux with respect to ice surface
                                               lh_v, & !                                 water
                                               l_sice(ic) ) !frozen surface of lake or any sea tile
            !Note(MR):
            !So far, for a frozen land-surface, always 'lh_v' is used. However, '%lhfl_s_t' is going to be overwritten
            ! for land-tiles by 'terra'; hence, this calculation only matters for non-land-tiles (lakes, sea-water or sea-ice).
            !While, for land points, '%lhfl_s' (as well as '%lhfl_s_t') is only used for model-output,
            ! '%lhfl_s_t' (as well as '%shfl_t') is used as input for the lake  and sea-ice-schemes.
          END DO
          !$ACC END PARALLEL

          !$ACC PARALLEL ASYNC(acc_async_queue) DEFAULT(PRESENT) IF(lzacc)
          !$ACC LOOP GANG VECTOR
!$NEC ivdep
          DO ic= 1, i_count
            jc = ilist(ic)

            ! Store
            prm_diag%shfl_s_t(jc,jb,jt) = shfl_s_t(ic) ! not neces., updated in TERRA (aggreg. in mo_nwp_sfc_interface)
            prm_diag%lhfl_s_t(jc,jb,jt) = lhfl_s_t(ic) ! not neces., updated in TERRA (aggreg. in mo_nwp_sfc_interface)
            prm_diag%qhfl_s_t(jc,jb,jt) = qhfl_s_t(ic) ! not neces., updated in TERRA (aggreg. in mo_nwp_sfc_interface)
            prm_diag%umfl_s_t(jc,jb,jt) = umfl_s_t(ic) ! being updated in mo_nwp_sfc_interface
            prm_diag%vmfl_s_t(jc,jb,jt) = vmfl_s_t(ic) ! being updated in mo_nwp_sfc_interface
            prm_diag%u_10m_t (jc,jb,jt) = u_10m_t(ic)  ! needed by TERRA and 'turbtran' (for some empirical parameteriz.)
            prm_diag%v_10m_t (jc,jb,jt) = v_10m_t(ic)  ! needed by TERRA and 'turbtran' (for some empirical parameteriz.)
            prm_diag%tch_t   (jc,jb,jt) = tch_t(ic)    ! needed by TERRA (ToDO: to be substituted by 'tvh_t')
            prm_diag%tcm_t   (jc,jb,jt) = tcm_t(ic)    ! needed by TERRA (ToDO: to be substituted by 'tvm_t')

            prm_diag%tfv_t   (jc,jb,jt) = tfv_t(ic)    ! needed by TERRA

            prm_diag%tvm_t   (jc,jb,jt) = tvm_t(ic)    ! to be used by TERRA instead of 'tcm'
            prm_diag%tvh_t   (jc,jb,jt) = tvh_t(ic)    ! to be used by TERRA instead of 'tch'

            prm_diag%tkr_t   (jc,jb,jt) = tkr_t(ic)    ! needed for next call of 'turbtran' if "imode_trancnf>=4"

            IF (jt>ntiles_total) THEN !for non-land (sub-)tiles (ice- or water-surfaces) only
              prm_diag%gz0_t(jc,jb,jt) = gz0_eff_t(ic,jt) ! modified value needed for next call of 'turbtran'
            END IF

            prm_diag%rcld_s_t(jc,jb,jt) = rcld_t(ic,3) ! needed for next call of 'turbtran'
            prm_diag%tvs_s_t (jc,jb,jt) =  tvs_t(ic,3) ! needed for next call of 'turbtran'
            prm_diag%tkvm_s_t(jc,jb,jt) = tkvm_t(ic,3) ! needed for next call of 'turbtran'
            prm_diag%tkvh_s_t(jc,jb,jt) = tkvh_t(ic,3) ! needed for next call of 'turbtran'

          ENDDO
          !$ACC END PARALLEL

          !$ACC END DATA

        !-----------------------------
        ENDDO ! loop over tiles
        !-----------------------------

        IF (multi_queue_processing) THEN
          ! Join into sequential GPU-streams for each tile:
          DO jt = 1, ntiles_total + ntiles_water
            !$ACC WAIT(jt) ASYNC(1)
          END DO
          !Note: For tile-aggregation, sequential GPU-streams "ASYNC(1)" are mandatory!
        END IF

        ! Aggregate tile-based output fields of turbtran over tiles
        ! i) initialize fields to zero before starting the summation

        !$ACC KERNELS ASYNC(1) DEFAULT(PRESENT) IF(lzacc)

        prm_diag%gz0   (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%tcm   (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%tch   (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%tfm   (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%tfh   (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%tvm   (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%tvh   (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%t_2m  (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%qv_2m (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%td_2m (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%rh_2m (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%u_10m (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%v_10m (i_startidx:i_endidx,jb) = 0._wp

        prm_diag%t_2m_land (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%td_2m_land(i_startidx:i_endidx,jb) = 0._wp
        prm_diag%rh_2m_land(i_startidx:i_endidx,jb) = 0._wp

        z_tvs        (i_startidx:i_endidx,3)         = 0._wp
        prm_diag%tkvm(i_startidx:i_endidx,nlevp1,jb) = 0._wp
        prm_diag%tkvh(i_startidx:i_endidx,nlevp1,jb) = 0._wp
        prm_diag%rcld(i_startidx:i_endidx,nlevp1,jb) = 0._wp

        ! re-initialize dyn_gust for proper maximum computation over tiles
        ! note that dyn_gust is an instantaneous field
        prm_diag%dyn_gust(i_startidx:i_endidx,jb) = 0._wp

        prm_diag%t_tilemax_inst_2m(i_startidx:i_endidx,jb) = -999._wp
        prm_diag%t_tilemin_inst_2m(i_startidx:i_endidx,jb) = 999._wp

        !$ACC END KERNELS

        ! ii) loop over index lists

        !$ACC DATA CREATE(i_count) ASYNC(1) IF(lzacc)
        !-----------------------------
        DO jt = 1, ntiles_total + ntiles_water
        !-----------------------------

          !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          i_count = gp_num_t(jt)
          !$ACC END KERNELS
          ilist => list_t(jt)%gp_idx

#ifndef _OPENACC
          IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty
#endif

          !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
          !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jc, area_frac)
!$NEC ivdep
          DO ic = 1, i_count
            jc = ilist(ic)

            ! Aggregate
            area_frac = ext_data%atm%frac_t(jc,jb,jt)

            prm_diag%gz0(jc,jb) = prm_diag%gz0(jc,jb) + gz0_eff_t     (ic,   jt)*area_frac !for TURBDIFF

            prm_diag%tcm(jc,jb) = prm_diag%tcm(jc,jb) + prm_diag%tcm_t(jc,jb,jt)*area_frac !only for output
            prm_diag%tch(jc,jb) = prm_diag%tch(jc,jb) + prm_diag%tch_t(jc,jb,jt)*area_frac !only for output
            prm_diag%tfm(jc,jb) = prm_diag%tfm(jc,jb) +          tfm_t(ic,jt)   *area_frac !for TURBDIFF
            prm_diag%tfh(jc,jb) = prm_diag%tfh(jc,jb) +          tfh_t(ic,jt)   *area_frac !for TURBDIFF
            prm_diag%tvm(jc,jb) = prm_diag%tvm(jc,jb) + prm_diag%tvm_t(jc,jb,jt)*area_frac !for TURBDIFF
            prm_diag%tvh(jc,jb) = prm_diag%tvh(jc,jb) + prm_diag%tvh_t(jc,jb,jt)*area_frac !for TURBDIFF

            z_tvs(jc,3)                 = z_tvs        (jc,3)         + prm_diag%tvs_s_t (jc,jb,jt)*area_frac
            prm_diag%tkvm(jc,nlevp1,jb) = prm_diag%tkvm(jc,nlevp1,jb) + prm_diag%tkvm_s_t(jc,jb,jt)*area_frac !for TURBDIFF
            prm_diag%tkvh(jc,nlevp1,jb) = prm_diag%tkvh(jc,nlevp1,jb) + prm_diag%tkvh_s_t(jc,jb,jt)*area_frac !for TURBDIFF
            prm_diag%rcld(jc,nlevp1,jb) = prm_diag%rcld(jc,nlevp1,jb) + prm_diag%rcld_s_t(jc,jb,jt)*area_frac !for TURBDIFF

            prm_diag%u_10m(jc,jb) = prm_diag%u_10m(jc,jb) + prm_diag%u_10m_t(jc,jb,jt)*area_frac !only for output
            prm_diag%v_10m(jc,jb) = prm_diag%v_10m(jc,jb) + prm_diag%v_10m_t(jc,jb,jt)*area_frac !only for output
            prm_diag%t_2m (jc,jb) = prm_diag%t_2m (jc,jb) + t_2m_t (ic,jt)*area_frac !only for output
            prm_diag%qv_2m(jc,jb) = prm_diag%qv_2m(jc,jb) + qv_2m_t(ic,jt)*area_frac !only for output
            prm_diag%td_2m(jc,jb) = prm_diag%td_2m(jc,jb) + td_2m_t(ic,jt)*area_frac !only for output
            prm_diag%rh_2m(jc,jb) = prm_diag%rh_2m(jc,jb) + rh_2m_t(ic,jt)*area_frac !only for output

            prm_diag%t_tilemax_inst_2m(jc,jb) = MAX(t_2m_t(ic,jt),prm_diag%t_tilemax_inst_2m(jc,jb))
            prm_diag%t_tilemin_inst_2m(jc,jb) = MIN(t_2m_t(ic,jt),prm_diag%t_tilemin_inst_2m(jc,jb))

          ENDDO

          ! averages over land fraction of mixed land-water points
          IF (jt <= ntiles_total) THEN
            !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jc, area_frac)
!$NEC ivdep
            DO ic = 1, i_count
              jc = ilist(ic)
              area_frac = ext_data%atm%frac_t(jc,jb,jt)/ext_data%atm%fr_land(jc,jb)
              prm_diag%t_2m_land (jc,jb) = prm_diag%t_2m_land (jc,jb) +  t_2m_t(ic,jt)*area_frac
              prm_diag%td_2m_land(jc,jb) = prm_diag%td_2m_land(jc,jb) + td_2m_t(ic,jt)*area_frac
              prm_diag%rh_2m_land(jc,jb) = prm_diag%rh_2m_land(jc,jb) + rh_2m_t(ic,jt)*area_frac
            ENDDO
          ELSE
            !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jc, area_frac)
!$NEC ivdep
            DO ic = 1, i_count
              jc = ilist(ic)
              IF (ext_data%atm%fr_land(jc,jb) == 0._wp) THEN
                area_frac = ext_data%atm%frac_t(jc,jb,jt)
                prm_diag%t_2m_land (jc,jb) = prm_diag%t_2m_land (jc,jb) +  t_2m_t(ic,jt)*area_frac
                prm_diag%td_2m_land(jc,jb) = prm_diag%td_2m_land(jc,jb) + td_2m_t(ic,jt)*area_frac
                prm_diag%rh_2m_land(jc,jb) = prm_diag%rh_2m_land(jc,jb) + rh_2m_t(ic,jt)*area_frac
              ENDIF
            ENDDO
          ENDIF
          !$ACC END PARALLEL

        !-----------------------------
        ENDDO  ! loop over tiles
        !-----------------------------
        !$ACC END DATA

      ENDIF ! tiles / no tiles

      IF (itune_gust_diag < 4) THEN
        ! Dynamic gusts are diagnosed from averaged values in order to avoid artifacts along coastlines
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = i_startidx, i_endidx
          prm_diag%dyn_gust(jc,jb) =  nwp_dyn_gust (prm_diag%u_10m(jc,jb),      &
            &                                       prm_diag%v_10m(jc,jb),      &
            &                                       prm_diag%tcm  (jc,jb),      &
            &                                       p_diag%u      (jc,nlev,jb), &
            &                                       p_diag%v      (jc,nlev,jb), &
            &                                       p_diag%u(jc,jk_gust(jc),jb),&
            &                                       p_diag%v(jc,jk_gust(jc),jb),&
            &                          ext_data%atm%lc_frac_t(jc,jb,isub_water),&
            &                                  p_metrics%mask_mtnpoints_g(jc,jb))
        ENDDO
        !$ACC END PARALLEL
      ELSE ! compute only the gust limiter; the gust calculation itself is executed at the end of each averaging interval
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG(STATIC: 1)
        DO jk = nlev, kstart_moist(jg), -1
          !$ACC LOOP VECTOR
          DO jc = i_startidx, i_endidx
            IF (p_metrics%geopot_agl(jc,jk,jb) < MAX(tune_gustlim_agl(jg)*grav, &
                p_metrics%geopot_agl(jc,jk_gust(jc),jb) + 500._wp*grav)) THEN
              prm_diag%gust_lim(jc,jb) = MAX(prm_diag%gust_lim(jc,jb), SQRT(p_diag%u(jc,jk,jb)**2 + p_diag%v(jc,jk,jb)**2))
            ENDIF
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDIF

      ! Transform updated (and aggregated) turbulent velocity scale (TVS) at the surface back to TKE:
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = i_startidx, i_endidx
        p_prog_rcf%tke(jc,nlevp1,jb)= 0.5_wp*(z_tvs(jc,3))**2
      ENDDO

      ! Interpolate updated turbulent velocity scale TVS at lowest main level:
      IF (advection_config(jg)%iadv_tke > 0) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc=i_startidx, i_endidx
          p_prog_rcf%tracer(jc,nlev,jb,iqtke) = 0.5_wp* ( z_tvs(jc,2) + z_tvs(jc,3) )
        ENDDO
      ENDIF

      !$ACC END PARALLEL

    CASE(igme,ismag,iprog)

!-------------------------------------------------------------------------
!> GME turbulence scheme
!-------------------------------------------------------------------------

      ! turbulent diffusion coefficients at the surface
      CALL parturs( zsurf=p_metrics%z_ifc(:,nlevp1,jb), z1=p_metrics%z_mc(:,nlev,jb),   & !in
        &           u1=p_diag%u(:,nlev,jb), v1=p_diag%v(:,nlev,jb),                     & !in
        &           t1=p_diag%temp(:,nlev,jb), qv1=p_prog_rcf%tracer(:,nlev,jb,iqv),    & !in
        &           t_g=lnd_prog_new%t_g(:,jb), qv_s=lnd_diag%qv_s(:,jb),               & !in
        &           ps=p_diag%pres_ifc(:,nlevp1,jb),                                    & !in
        &           fr_land=ext_data%atm%fr_land(:,jb), h_ice=wtr_prog_new%h_ice(:,jb), & !in
        &           ie=nproma, i_startidx=i_startidx, i_endidx=i_endidx,                & !in
        &           tcm=prm_diag%tcm(:,jb), tch=prm_diag%tch(:,jb),                     & !out
        &           gz0=prm_diag%gz0(:,jb),       shfl_s=prm_diag%shfl_s(:,jb),         & !inout, out: for bug-emulation only
        &           lhfl_s=prm_diag%lhfl_s(:,jb), qhfl_s=prm_diag%qhfl_s(:,jb),         & !out, out
        &           umfl_s=prm_diag%umfl_s(:,jb), vmfl_s=prm_diag%vmfl_s(:,jb),         & !out, out
        &           lacc=lzacc                                                          ) !in


      !DR inside "nearsfc", lhfl_s is converted to qhfl_s via
      !DR qhfl_s = lhfl_s/lh_v. This is incorrect over snow and ice.
      !DR Shouldn't we simply pass qhfl_s ?
      !
      ! diagnose 2 m temperature, humidity, 10 m wind
      CALL nearsfc( t=p_diag%temp(:,:,jb), qv=p_prog_rcf%tracer(:,:,jb,iqv),            & !in
        &           u=p_diag%u(:,:,jb),    v=p_diag%v(:,:,jb),                          & !in
        &           zf=p_metrics%z_mc(:,:,jb), ps=p_diag%pres_ifc(:,nlevp1,jb),         & !in
        &           t_g=lnd_prog_new%t_g(:,jb),                                         & !in
        &           tcm=prm_diag%tcm(:,jb), tch=prm_diag%tch(:,jb),                     & !in
        &           gz0=prm_diag%gz0(:,jb),                                             & !in
        &           shfl_s=prm_diag%shfl_s(:,jb), lhfl_s=prm_diag%lhfl_s(:,jb),         & !in
        &           umfl_s=prm_diag%umfl_s(:,jb), vmfl_s=prm_diag%vmfl_s(:,jb),         & !in
        &           zsurf=p_metrics%z_ifc(:,nlevp1,jb),                                 & !in
        &           fr_land=ext_data%atm%fr_land(:,jb), pf1=p_diag%pres(:,nlev,jb),     & !in
        &           qv_s=lnd_diag%qv_s(:,jb), ie=nproma, ke=nlev,                       & !in
        &           i_startidx=i_startidx, i_endidx=i_endidx,                           & !in
        &           t_2m=prm_diag%t_2m(:,jb), qv_2m=prm_diag%qv_2m(:,jb),               & !out
        &           td_2m=prm_diag%td_2m(:,jb), rh_2m=prm_diag%rh_2m(:,jb),             & !out
        &           u_10m=prm_diag%u_10m(:,jb), v_10m=prm_diag%v_10m(:,jb),             & !out
        &           lacc=lzacc                                                          ) !in

      ! dynamic gusts
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        prm_diag%dyn_gust(jc,jb) = nwp_dyn_gust(prm_diag%u_10m(jc,jb), prm_diag%v_10m(jc,jb), &
          &  prm_diag%tcm(jc,jb), p_diag%u(jc,nlev,jb), p_diag%v(jc,nlev,jb),                 &
          &  p_diag%u(jc,jk_gust(jc),jb), p_diag%v(jc,jk_gust(jc),jb),                        &
          &  ext_data%atm%lc_frac_t(jc,jb,isub_water),p_metrics%mask_mtnpoints_g(jc,jb) )
      ENDDO
      !$ACC END PARALLEL

      ! instantaneous max/min 2m temperature over tiles (trivial operation for 1 tile)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        prm_diag%t_tilemax_inst_2m(jc,jb) = prm_diag%t_2m(jc,jb)
        prm_diag%t_tilemin_inst_2m(jc,jb) = prm_diag%t_2m(jc,jb)
      ENDDO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jt = 1, ntiles_total+ntiles_water
        DO jc = i_startidx, i_endidx

          ! Copy transfer coefficients to tile-based variables, which are used in TERRA
          prm_diag%tcm_t(jc,jb,jt) = prm_diag%tcm(jc,jb)
          prm_diag%tch_t(jc,jb,jt) = prm_diag%tch(jc,jb)
          ! the GME turbulence scheme does not have tfv. Set tfv=1
          prm_diag%tfv_t(jc,jb,jt) = 1._wp    !   prm_diag%tfv(jc,jb)

          ! Copy transfer u_10m/v_10m to tile-based variables, which are used in TERRA
          prm_diag%u_10m_t(jc,jb,jt) = prm_diag%u_10m(jc,jb)
          prm_diag%v_10m_t(jc,jb,jt) = prm_diag%v_10m(jc,jb)

          ! Copy sensible and latent heat fluxes to tile-based variables
          ! (needed by Flake, sea-ice model)
          prm_diag%shfl_s_t(jc,jb,jt) = prm_diag%shfl_s(jc,jb)
          prm_diag%lhfl_s_t(jc,jb,jt) = prm_diag%lhfl_s(jc,jb)
          prm_diag%qhfl_s_t(jc,jb,jt) = prm_diag%qhfl_s(jc,jb)

          prm_diag%umfl_s_t(jc,jb,jt) = prm_diag%umfl_s(jc,jb)
          prm_diag%vmfl_s_t(jc,jb,jt) = prm_diag%vmfl_s(jc,jb)
        ENDDO
      ENDDO
      !$ACC END PARALLEL


    END SELECT !inwp_turb


    ! Compute wind speed in 10m
    ! used by mo_albedo (albedo_whitecap=1)
    ! 
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO jc = i_startidx, i_endidx
      prm_diag%sp_10m(jc,jb) = SQRT(prm_diag%u_10m(jc,jb)**2 + prm_diag%v_10m(jc,jb)**2 )
    ENDDO
    !$ACC END PARALLEL

  ENDDO ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  !$ACC END DATA

#ifdef ICON_USE_CUDA_GRAPH
    IF (lzacc .AND. lcuda_graph_turb_tran) THEN
      CALL accEndCapture(1, graphs(cur_graph_id))
      WRITE(message_text,'(a,i2,a)') 'finished to capture CUDA graph, id ', cur_graph_id, ', now executing it'
      IF (msg_level >= 13) CALL message('mo_nwp_turbtrans_interface: ', message_text)
      CALL accGraphLaunch(graphs(cur_graph_id), 1)
    END IF
#endif

  !$ACC WAIT(1)
  IF (timers_level > 9) CALL timer_stop(timer_nwp_turbtrans)

END SUBROUTINE nwp_turbtrans

END MODULE mo_nwp_turbtrans_interface
