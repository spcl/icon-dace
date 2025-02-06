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

! Setup routines for adaptive parameter tuning and related computations

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_apt_routines

  USE mo_kind,                ONLY: wp
  USE mo_math_constants,      ONLY: rad2deg, pi2
  USE mo_physical_constants,  ONLY: tmelt
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
  USE mo_nwp_lnd_types,       ONLY: t_wtr_prog, t_lnd_diag, t_lnd_prog
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_state
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: min_rlcell_int, min_rlcell
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_parallel_config,     ONLY: nproma
  USE mo_dynamics_config,     ONLY: nnow
  USE mo_time_config,         ONLY: time_config
  USE mo_grid_config,         ONLY: n_dom
  USE mo_lnd_nwp_config,      ONLY: ntiles_total, ntiles_water, ntiles_lnd, &
    &                               itype_canopy, itype_lndtbl, c_soil, c_soil_urb,  &
                                    lterra_urb, itype_eisa, cr_bsmin, nlev_soil, dzsoil, depth_hl, zml_soil
  USE sfc_terra_data,         ONLY: cporv, cadp, cpwp, cfcap
  USE mo_satad,               ONLY: sat_pres_water, &  !! saturation vapor pressure w.r.t. water
    &                               spec_humi          !! Specific humidity
  USE mo_input_instructions,  ONLY: t_readInstructionListPtr, kInputSourceAna, kInputSourceAnaI

  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_cell

  USE mo_master_config,       ONLY: isRestart
  USE mo_initicon_types,      ONLY: t_initicon_state
  USE mo_initicon_config,     ONLY: icpl_da_sfcevap, dt_ana, icpl_da_snowalb, icpl_da_landalb, icpl_da_skinc, &
                                    icpl_da_sfcfric, icpl_da_tkhmin, icpl_da_seaice, scalfac_da_sfcfric

  USE mo_nwp_tuning_config,   ONLY: itune_slopecorr



  IMPLICIT NONE

  PRIVATE


  PUBLIC  :: compute_filtincs, init_apt_fields, apply_landalb_tuning, apply_sma


  CONTAINS



  !-------------------------------------------------------------------------
  !>
  !! SUBROUTINE compute_filtincs
  !!
  !! Computes filtered assimilation increments for adaptive parameter tuning
  !!
  !-------------------------------------------------------------------------
  SUBROUTINE compute_filtincs (p_patch, p_nh_state, p_int_state, initicon, inputInstructions)

    TYPE(t_patch)                 ,INTENT(IN)    :: p_patch(:)
    TYPE(t_nh_state) , TARGET     ,INTENT(IN)    :: p_nh_state(:)
    TYPE(t_int_state),             INTENT(IN)    :: p_int_state(:)
    TYPE(t_initicon_state)        ,INTENT(INOUT) :: initicon(:)
    TYPE(t_readInstructionListPtr),INTENT(IN)    :: inputInstructions(:)

    TYPE(t_nh_diag) , POINTER :: p_diag
    TYPE(t_nh_prog),  POINTER :: p_prog_now

    INTEGER :: jg, jb, jc, nlev
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startidx, i_endidx, i_startblk, i_endblk


    REAL(wp) :: rh_inc(nproma), localtime_fac, dtfac

  !-------------------------------------------------------------------------

    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      p_diag      => p_nh_state(jg)%diag
      p_prog_now  => p_nh_state(jg)%prog(nnow(jg))

      nlev      = p_patch(jg)%nlev

      rl_start  = grf_bdywidth_c+1
      rl_end    = min_rlcell_int

      i_startblk = p_patch(jg)%cells%start_block(rl_start)
      i_endblk   = p_patch(jg)%cells%end_block(rl_end)

      ! weighting factor for standard filtering time scale of 2.5 days
      dtfac      = dt_ana/216000._wp

      ! interpolate wind and its increments to mass points
      IF (icpl_da_sfcfric >= 1) THEN
        CALL rbf_vec_interpol_cell(p_prog_now%vn, p_patch(jg), p_int_state(jg), p_diag%u, p_diag%v)
        CALL rbf_vec_interpol_cell(initicon(jg)%atm_inc%vn, p_patch(jg), p_int_state(jg), &
          initicon(jg)%atm_inc%u, initicon(jg)%atm_inc%v)
      ENDIF


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,rh_inc,localtime_fac)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)


        IF (icpl_da_sfcevap == 1 .OR. icpl_da_sfcevap == 2) THEN
          ! Time-filtering of analyzed T2M bias; only if t_2m is read from analysis
          IF (ANY((/kInputSourceAna,kInputSourceAnaI/) == inputInstructions(jg)%ptr%sourceOfVar('t_2m'))) THEN
            DO jc = i_startidx, i_endidx
              p_diag%t2m_bias(jc,jb) = p_diag%t2m_bias(jc,jb) + &
                0.4_wp*(initicon(jg)%sfc_inc%t_2m(jc,jb)-p_diag%t2m_bias(jc,jb))
            ENDDO
          ENDIF
        ELSE IF (icpl_da_sfcevap >= 3) THEN
          ! Calculate time-filtered T assimilation increment at lowest model level (time scale 2.5 days);
          ! this serves as a proxy for the averaged T2M bias
          DO jc = i_startidx, i_endidx
            p_diag%t_avginc(jc,jb) = p_diag%t_avginc(jc,jb) + &
              dtfac*(initicon(jg)%atm_inc%temp(jc,nlev,jb)-p_diag%t_avginc(jc,jb))
          ENDDO
        ENDIF

        IF (icpl_da_sfcevap >= 2) THEN
          ! Calculate time-filtered RH assimilation increment at lowest model level (time scale 2.5 days);
          ! this serves as a proxy for the averaged RH2M bias
          DO jc = i_startidx, i_endidx
            rh_inc(jc) = initicon(jg)%atm_inc%qv(jc,nlev,jb)/ &
              spec_humi(sat_pres_water(p_diag%temp(jc,nlev,jb)),p_diag%pres_sfc(jc,jb))
            p_diag%rh_avginc(jc,jb) = p_diag%rh_avginc(jc,jb) + dtfac*(rh_inc(jc)-p_diag%rh_avginc(jc,jb))
          ENDDO
        ENDIF

        IF (icpl_da_skinc >= 1) THEN
          ! weighted T assimilation increment; this serves as a proxy for the bias in diurnal temperature amplitude
          DO jc = i_startidx, i_endidx
            localtime_fac = COS(p_patch(jg)%cells%center(jc,jb)%lon + pi2/86400._wp *                           &
              (time_config%tc_exp_startdate%time%hour*3600._wp+time_config%tc_exp_startdate%time%minute*60._wp) )
            p_diag%t_wgt_avginc(jc,jb) = p_diag%t_wgt_avginc(jc,jb) + &
              dtfac*(initicon(jg)%atm_inc%temp(jc,nlev,jb)*localtime_fac-p_diag%t_wgt_avginc(jc,jb))
          ENDDO
        ENDIF

        IF (icpl_da_sfcfric >= 1) THEN
          ! weighted wind speed increment for adaptive surface friction
          DO jc = i_startidx, i_endidx
            p_diag%vabs_avginc(jc,jb) = p_diag%vabs_avginc(jc,jb) + dtfac * (                     &
              SQRT( (p_diag%u(jc,nlev,jb)+initicon(jg)%atm_inc%u(jc,nlev,jb))**2 +                &
                    (p_diag%v(jc,nlev,jb)+initicon(jg)%atm_inc%v(jc,nlev,jb))**2 ) -              &
              SQRT(p_diag%u(jc,nlev,jb)**2 + p_diag%v(jc,nlev,jb)**2) - p_diag%vabs_avginc(jc,jb) )
          ENDDO
        ENDIF

        IF (icpl_da_sfcevap >= 5) THEN
          ! Daytime-weighted T and RH increments peaking at 16:30 LT, approximating the time at which the
          ! maximum temperature is reached in summer
          DO jc = i_startidx, i_endidx
            localtime_fac = MAX(0._wp, -SIN(p_patch(jg)%cells%center(jc,jb)%lon + pi2/16._wp + pi2/86400._wp *   &
              (time_config%tc_exp_startdate%time%hour*3600._wp+time_config%tc_exp_startdate%time%minute*60._wp) ))
            p_diag%t_daywgt_avginc(jc,jb) = p_diag%t_daywgt_avginc(jc,jb) +           &
              dtfac*localtime_fac*(initicon(jg)%atm_inc%temp(jc,nlev,jb)-p_diag%t_daywgt_avginc(jc,jb))
            p_diag%rh_daywgt_avginc(jc,jb) = p_diag%rh_daywgt_avginc(jc,jb) +         &
              dtfac*localtime_fac*(rh_inc(jc)-p_diag%rh_daywgt_avginc(jc,jb))
          ENDDO
        ENDIF

      ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ENDDO  ! jg

  END SUBROUTINE compute_filtincs


  !-------------------------------------------------------------------------
  !>
  !! SUBROUTINE init_apt_fields
  !!
  !! Initialization of auxiliary fields for adaptive parameter tuning
  !!
  !-------------------------------------------------------------------------
  SUBROUTINE init_apt_fields (p_patch, p_diag, prm_diag, ext_data, p_diag_lnd, p_prog_wtr_now)

    TYPE(t_patch),           INTENT(in)    :: p_patch
    TYPE(t_nh_diag),         INTENT(in)    :: p_diag
    TYPE(t_nwp_phy_diag),    INTENT(inout) :: prm_diag
    TYPE(t_external_data),   INTENT(inout) :: ext_data
    TYPE(t_lnd_diag),        INTENT(in)    :: p_diag_lnd
    TYPE(t_wtr_prog),        INTENT(in)    :: p_prog_wtr_now


    INTEGER :: jb, ic, jc, jt
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !! slices
    INTEGER :: i_count, ilu


    REAL(wp) :: dtfac, dtfac_heatc, tbias_wgt, zlat, zlon, slope(nproma), skinc_fac, trh_bias, scal, rlamh_offset, rlamh_offset_glac, &
                zrlamh

    ! adaptation factor to analysis interval, standard value and value for adaptive heat conductivity/capacity and skin conductivity
    dtfac       =  10800._wp/dt_ana
    dtfac_heatc = (10800._wp/dt_ana)**(2._wp/3._wp)

    ! offset for adaptation of rlam_heat
    rlamh_offset      = MERGE(0.4_wp, 0.2_wp, icpl_da_sfcevap <= 4)
    rlamh_offset_glac =       0.4_wp ! always use 0.4 over full glacier points

    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jt,i_startidx,i_endidx,i_count,ilu,tbias_wgt,zlat,zlon,slope,skinc_fac,trh_bias,scal,zrlamh)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

 
      ! tuning factor for rlam_heat depending on skin conductivity and analyzed T2M/RH2M bias
      IF (itype_canopy == 2 .AND. icpl_da_sfcevap >= 3) THEN
        DO jt = 1, ntiles_total + ntiles_water
          DO jc = i_startidx,i_endidx
            IF (jt <= ntiles_lnd) THEN ! snow-free land points
              prm_diag%rlamh_fac_t(jc,jb,jt) = 1._wp - 0.9_wp*MAX(0._wp, MIN(1._wp,                         &
                2.5_wp*(dtfac*(100._wp*p_diag%rh_avginc(jc,jb)-4._wp*p_diag%t_avginc(jc,jb))-rlamh_offset) ))
            ELSE IF (jt <= ntiles_total) THEN ! snow-covered land points
              zrlamh = MERGE(rlamh_offset, rlamh_offset_glac, ext_data%atm%fr_glac(jc,jb) <= 0.99_wp)
              prm_diag%rlamh_fac_t(jc,jb,jt) = 1._wp - 0.9_wp*MAX(0._wp, MIN(1._wp, &
                2.5_wp*(dtfac*(MAX(0._wp,100._wp*p_diag%rh_avginc(jc,jb))-4._wp*p_diag%t_avginc(jc,jb))-zrlamh) ))
            ELSE IF (jt == ntiles_total + ntiles_water) THEN ! seaice points
              prm_diag%rlamh_fac_t(jc,jb,jt) = 0.25_wp
            ELSE
              prm_diag%rlamh_fac_t(jc,jb,jt) = 1._wp
            ENDIF
          ENDDO
        ENDDO
      ELSE IF (itype_canopy == 2 .AND. icpl_da_sfcevap == 2) THEN
        DO jt = 1, ntiles_total + ntiles_water
          DO jc = i_startidx,i_endidx
            IF (jt <= ntiles_total) THEN
              prm_diag%rlamh_fac_t(jc,jb,jt) =                                                                         &
                1._wp - 0.9_wp*MAX(0._wp,MIN(1._wp,(60._wp-ext_data%atm%skinc_t(jc,jb,jt))/30._wp)) *                  &
                MAX(0._wp,MIN(1._wp,2.5_wp*(p_diag%t2m_bias(jc,jb)+100._wp*dtfac*p_diag%rh_avginc(jc,jb)-rlamh_offset)))
            ELSE IF (jt == ntiles_total + ntiles_water) THEN ! seaice points
              prm_diag%rlamh_fac_t(jc,jb,jt) = 0.25_wp
            ELSE
              prm_diag%rlamh_fac_t(jc,jb,jt) = 1._wp
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      IF (icpl_da_landalb >= 1 .AND. .NOT. isRestart()) THEN
        ! Tuning increment for land albedo, applied to snow-free dominant land points
        DO jc = i_startidx,i_endidx
          IF (MAXVAL(p_diag_lnd%h_snow_t(jc,jb,1:ntiles_total)) == 0._wp .AND. ext_data%atm%fr_land(jc,jb) >= 0.5_wp) THEN
            ! albedo tuning is applied when a cold (warm) bias concides with a dry (moist) bias, implying that the bias
            ! is not related to an incorrect Bowen ratio
            trh_bias = dtfac*(p_diag%t_daywgt_avginc(jc,jb) + 30._wp*(p_diag%rh_daywgt_avginc(jc,jb)+0.25_wp*p_diag%rh_avginc(jc,jb)))
            scal = MERGE(0.2_wp, 0.3_wp, trh_bias > 0._wp)
            prm_diag%landalb_inc(jc,jb) = SIGN(scal*trh_bias**2,-trh_bias)
          ENDIF
        ENDDO
      ENDIF

      IF (icpl_da_snowalb >= 1 .AND. .NOT. isRestart()) THEN
        ! Tuning factor for snow albedo, including coastal sea ice for icpl_da_snowalb>=2
        DO jc = i_startidx,i_endidx
          IF (ANY(p_diag_lnd%h_snow_t(jc,jb,1:ntiles_total) > 0._wp) .OR. p_prog_wtr_now%h_ice(jc,jb) > 0._wp) THEN
            IF (p_diag%t_avginc(jc,jb) > 0._wp) THEN
              prm_diag%snowalb_fac(jc,jb) = MAX(0.75_wp,1._wp/(1._wp+dtfac*0.8_wp*p_diag%t_avginc(jc,jb)))
            ELSE
              prm_diag%snowalb_fac(jc,jb) = MIN(4._wp/3._wp,1._wp-dtfac*0.8_wp*p_diag%t_avginc(jc,jb))
            ENDIF
          ENDIF
          IF (icpl_da_snowalb >= 2) THEN ! albedo factor is also applied to sea ice and needs to be restricted to the vicinity of land
            scal = MIN(1._wp,100._wp*ext_data%atm%fr_land_smt(jc,jb))
            prm_diag%snowalb_fac(jc,jb) = scal*prm_diag%snowalb_fac(jc,jb) + (1._wp-scal)
          ENDIF
          IF (icpl_da_snowalb >= 3) THEN ! additional adjustment of snow-cover fraction diagnosis, relevant for small snow depths
            tbias_wgt = dtfac*(p_diag%t_avginc(jc,jb)+0.5_wp*p_diag%t_wgt_avginc(jc,jb))
            IF (tbias_wgt < 0._wp) THEN
              prm_diag%snowfrac_fac(jc,jb) = 1._wp/MAX(0.25_wp, 1._wp+1.5_wp*tbias_wgt)
            ELSE
              prm_diag%snowfrac_fac(jc,jb) = MAX(0.25_wp, 1._wp-1.5_wp*tbias_wgt)
            ENDIF
          ENDIF
        ENDDO
      ENDIF

      IF (icpl_da_seaice >= 2) THEN
        ! Tuning factor for sea ice bottom heat flux
        DO jc = i_startidx,i_endidx
          prm_diag%hflux_si_fac(jc,jb) = MIN(1._wp,MAX(0._wp,-5._wp*p_diag%t_avginc(jc,jb))) * &
            MIN(1._wp,100._wp*ext_data%atm%fr_land_smt(jc,jb))
        ENDDO
      ENDIF

      IF (icpl_da_skinc >= 2) THEN
        ! Tuning factors for soil heat capacity and conductivity
        DO jc = i_startidx,i_endidx
          IF (p_diag%t_wgt_avginc(jc,jb) < 0._wp) THEN
            prm_diag%heatcond_fac(jc,jb) = MAX(0.1_wp,  1._wp+dtfac_heatc*2.5_wp*p_diag%t_wgt_avginc(jc,jb))
            prm_diag%heatcap_fac(jc,jb)  = MAX(0.25_wp, 1._wp+dtfac_heatc*2.0_wp*p_diag%t_wgt_avginc(jc,jb))
          ELSE
            prm_diag%heatcond_fac(jc,jb) = 1._wp/MAX(0.1_wp,  1._wp-dtfac_heatc*2.5_wp*p_diag%t_wgt_avginc(jc,jb)) 
            prm_diag%heatcap_fac(jc,jb)  = 1._wp/MAX(0.25_wp, 1._wp-dtfac_heatc*2.0_wp*p_diag%t_wgt_avginc(jc,jb))
          ENDIF
        ENDDO
      ENDIF

      IF (icpl_da_tkhmin >= 1) THEN
        ! Adaptive tuning of near-surface minimum vertical diffusion for heat
        DO jc = i_startidx,i_endidx
          tbias_wgt = dtfac*(p_diag%t_avginc(jc,jb)+0.5_wp*p_diag%t_wgt_avginc(jc,jb))
          IF (tbias_wgt < 0._wp) THEN
            prm_diag%tkred_sfc_h(jc,jb) = MAX(0.25_wp, 1._wp+2._wp*tbias_wgt)
          ELSE
            prm_diag%tkred_sfc_h(jc,jb) = 1._wp/SQRT(MAX(0.25_wp, 1._wp-2._wp*tbias_wgt))
          ENDIF
        ENDDO
      ENDIF

      IF (itune_slopecorr >= 1) THEN
        ! Tuning for steep slopes (relevant for ICON-D05), involving a postprocessing of tkred_sfc_h and rlamh_fac_t
        DO jc = i_startidx,i_endidx
          slope(jc) = SQRT(ext_data%atm%grad_topo(1,jc,jb)**2 + ext_data%atm%grad_topo(2,jc,jb)**2)
          prm_diag%tkred_sfc_h(jc,jb) = prm_diag%tkred_sfc_h(jc,jb)/MIN(7.5_wp,1._wp+10._wp*SQRT(MAX(0._wp,slope(jc)-0.05_wp)))
        ENDDO
        DO jt = 1, ntiles_total + ntiles_water
          DO jc = i_startidx,i_endidx
            prm_diag%rlamh_fac_t(jc,jb,jt) = prm_diag%rlamh_fac_t(jc,jb,jt)/ &
              MIN(10._wp,1._wp+15._wp*SQRT(MAX(0._wp,slope(jc)-0.05_wp)))
           ENDDO
        ENDDO
      ENDIF

      IF (icpl_da_sfcfric >= 1) THEN
        ! Tuning factor for surface friction (roughness length and SSO blocking)
        DO jc = i_startidx,i_endidx
          IF (p_diag%vabs_avginc(jc,jb) > 0._wp) THEN
            prm_diag%sfcfric_fac(jc,jb) = MAX(0.25_wp, 1._wp-scalfac_da_sfcfric*dtfac*p_diag%vabs_avginc(jc,jb))
          ELSE
            prm_diag%sfcfric_fac(jc,jb) = 1._wp/MAX(0.25_wp, 1._wp+scalfac_da_sfcfric*dtfac*p_diag%vabs_avginc(jc,jb))
          ENDIF

          zlat = p_patch%cells%center(jc,jb)%lat*rad2deg
          zlon = p_patch%cells%center(jc,jb)%lon*rad2deg

          ! exclude Antarctic glaciers
          IF (ext_data%atm%fr_glac(jc,jb) > 0.99_wp .AND. zlat < -60._wp) prm_diag%sfcfric_fac(jc,jb) = 1._wp

          ! prevent reduction of surface friction in regions where 10m wind data are blacklisted
          ! use icpl_da_sfcfric = 2 in combination without blacklisting
          IF (icpl_da_sfcfric == 1 .AND.                                                          &
             (zlon >= 30._wp .AND. zlon <= 50._wp .AND. zlat >= 40._wp .AND. zlat <= 70._wp .OR.  &
              zlon >= 50._wp .AND. zlon <= 90._wp .AND. zlat >= 55._wp .AND. zlat <= 70._wp .OR.  &
              zlon >= 90._wp .AND. zlon <= 140._wp .AND. zlat >= 50._wp .AND. zlat <= 70._wp)) THEN 
            prm_diag%sfcfric_fac(jc,jb) = MAX(1._wp, prm_diag%sfcfric_fac(jc,jb))
          ENDIF

        ENDDO
      ENDIF


      DO jt = 1, ntiles_total
        i_count = ext_data%atm%lp_count_t(jb,jt)
        IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty

        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_lp_t(ic,jb,jt)

          ilu = ext_data%atm%lc_class_t(jc,jb,jt)

          ! compute combined temperature-humidity bias as optimization target for evaporation tuning
          !
          IF (icpl_da_sfcevap >= 3 .AND. icpl_da_skinc >= 1) THEN
            trh_bias = dtfac*(125._wp*p_diag%rh_avginc(jc,jb) - &
              4._wp*(p_diag%t_avginc(jc,jb)-0.5_wp*p_diag%t_wgt_avginc(jc,jb)))
          ELSE IF (icpl_da_sfcevap >= 3) THEN
            trh_bias = dtfac*(100._wp*p_diag%rh_avginc(jc,jb)-4._wp*p_diag%t_avginc(jc,jb))
          ELSE IF (icpl_da_sfcevap >= 2) THEN
            trh_bias = p_diag%t2m_bias(jc,jb) + 100._wp*dtfac*p_diag%rh_avginc(jc,jb)
          ELSE IF (icpl_da_sfcevap == 1) THEN
            trh_bias = p_diag%t2m_bias(jc,jb)
          ELSE
            trh_bias = 0._wp
          ENDIF


          IF (icpl_da_sfcevap >= 4) THEN
            ! minimum stomata resistance and bare soil evap resistance
            scal = MERGE(1.0_wp, 1.5_wp, icpl_da_sfcevap == 4) ! option 5 uses compensating asymmetric scaling for r_bsmin and hydiffu
            IF (trh_bias < 0._wp) THEN
              ext_data%atm%rsmin2d_t(jc,jb,jt) = ext_data%atm%stomresmin_lcc(ilu)*(1._wp-0.75_wp*trh_bias)
              ext_data%atm%r_bsmin(jc,jb)      = cr_bsmin*(1._wp-scal*trh_bias)
            ELSE
              ext_data%atm%rsmin2d_t(jc,jb,jt) = ext_data%atm%stomresmin_lcc(ilu)/(1._wp+0.75_wp*trh_bias)
              ext_data%atm%r_bsmin(jc,jb)      = cr_bsmin/(1._wp+trh_bias)
            ENDIF
          ELSE IF (icpl_da_sfcevap >= 1) THEN
            IF (trh_bias < 0._wp) THEN
              ext_data%atm%rsmin2d_t(jc,jb,jt) = ext_data%atm%stomresmin_lcc(ilu)*(1._wp-0.5_wp*trh_bias)
              ext_data%atm%eai_t(jc,jb,jt)     = MERGE(c_soil_urb,c_soil,ilu == ext_data%atm%i_lc_urban) /     &
                                                 (1._wp-0.25_wp*trh_bias)
            ELSE
              ext_data%atm%rsmin2d_t(jc,jb,jt) = ext_data%atm%stomresmin_lcc(ilu)/(1._wp+0.5_wp*trh_bias)
              ext_data%atm%eai_t(jc,jb,jt)     = MIN(MERGE(c_soil_urb,c_soil,ilu == ext_data%atm%i_lc_urban) * &
                                                 (1._wp+0.25_wp*trh_bias), 2._wp)
            ENDIF

            IF (lterra_urb .AND. ((itype_eisa == 2) .OR. (itype_eisa == 3))) THEN
              ext_data%atm%eai_t(jc,jb,jt)     = ext_data%atm%eai_t(jc,jb,jt)                                  &
                                               * (1.0_wp - ext_data%atm%urb_isa_t(jc,jb,jt))
            END IF
          ENDIF

          ! Tuning factor for skin conductivity
          IF (icpl_da_skinc >= 1) THEN
            scal = MERGE(4._wp, 2.5_wp, icpl_da_skinc == 1)
            IF (p_diag%t_wgt_avginc(jc,jb) < 0._wp) THEN
              skinc_fac = MAX(0.1_wp,1._wp+dtfac_heatc*scal*p_diag%t_wgt_avginc(jc,jb))
            ELSE
              skinc_fac = 1._wp/MAX(0.1_wp,1._wp-dtfac_heatc*scal*p_diag%t_wgt_avginc(jc,jb))
            ENDIF

            zlat = p_patch%cells%center(jc,jb)%lat*rad2deg
            IF (itype_lndtbl == 4 .AND. zlat > -10._wp .AND. zlat < 42.5_wp) THEN
              ext_data%atm%skinc_t(jc,jb,jt) = skinc_fac*MIN(200._wp,ext_data%atm%skinc_lcc(ilu)*          &
                                               (1._wp+MIN(1._wp,0.4_wp*(42.5_wp-zlat),0.4_wp*(zlat+10._wp))) )
            ELSE
              ext_data%atm%skinc_t(jc,jb,jt) = skinc_fac*ext_data%atm%skinc_lcc(ilu)
            ENDIF
          ENDIF

        ENDDO
      ENDDO

      ! Tuning factor for hydraulic diffusivity
      IF (icpl_da_sfcevap >= 5 .AND. icpl_da_skinc >= 1) THEN
        DO jc = i_startidx,i_endidx

          trh_bias = dtfac*(150._wp*(p_diag%rh_daywgt_avginc(jc,jb)+0.25_wp*p_diag%rh_avginc(jc,jb)) - &
                     3._wp*p_diag%t_daywgt_avginc(jc,jb))
          ! Asymmetric scaling (smaller for reduced than for increased diffusivity) because reducing the hydraulic
          ! diffusivity too strongly may lead to a dry bias increasing with forecast lead time. As a compensation,
          ! an opposing asymetric scaling is applied for r_bsmin (see above)
          IF (trh_bias < 0._wp) THEN
            prm_diag%hydiffu_fac(jc,jb)  = MAX(0.4_wp, 1._wp+0.25_wp*trh_bias)
          ELSE
            prm_diag%hydiffu_fac(jc,jb)  = 1._wp/MAX(0.2_wp, 1._wp-0.5_wp*trh_bias)
          ENDIF

        ENDDO
      ENDIF

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE init_apt_fields

  !-------------------------------------------------------------------------
  !>
  !! SUBROUTINE apply_landalb_tuning
  !!
  !! Applies adaptive parameter tuning to the time-interpolated MODIS albedo fields
  !!
  !-------------------------------------------------------------------------
  SUBROUTINE apply_landalb_tuning (p_patch, prm_diag, ext_data)

    TYPE(t_patch),           INTENT(in)    :: p_patch
    TYPE(t_nwp_phy_diag),    INTENT(in)    :: prm_diag
    TYPE(t_external_data),   INTENT(inout) :: ext_data
 
    INTEGER :: jb, jc
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !! slices


    REAL(wp) :: albuv, albni, albdif_fac

    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    IF (icpl_da_landalb >= 1) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,albuv,albni,albdif_fac)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        DO jc = i_startidx,i_endidx

          albuv = ext_data%atm%albuv_dif(jc,jb)+prm_diag%landalb_inc(jc,jb)
          albni = ext_data%atm%albni_dif(jc,jb)+prm_diag%landalb_inc(jc,jb)
          IF (prm_diag%landalb_inc(jc,jb) > 0._wp) THEN
            albuv = MIN(albuv, ext_data%atm%albuv_dif(jc,jb)+0.1_wp,  0.5_wp)
            albni = MIN(albni, ext_data%atm%albni_dif(jc,jb)+0.15_wp, 0.7_wp)
          ELSE
            albuv = MAX(albuv, ext_data%atm%albuv_dif(jc,jb)-0.1_wp,  0.02_wp)
            albni = MAX(albni, ext_data%atm%albni_dif(jc,jb)-0.15_wp, 0.08_wp)
          ENDIF
          albdif_fac = (albuv+albni)/(ext_data%atm%albuv_dif(jc,jb)+ext_data%atm%albni_dif(jc,jb))

          ext_data%atm%albuv_dif(jc,jb) = albuv
          ext_data%atm%albni_dif(jc,jb) = albni
          ext_data%atm%alb_dif(jc,jb)   = MAX(0.05_wp,ext_data%atm%alb_dif(jc,jb)*albdif_fac)

        ENDDO

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ENDIF

  END SUBROUTINE apply_landalb_tuning

  !-------------------------------------------------------------------------
  !>
  !! SUBROUTINE apply_sma
  !!
  !! Applies ICON-internal soil moisture adjustment building upon the APT predictors
  !!
  !-------------------------------------------------------------------------
  SUBROUTINE apply_sma (p_patch, p_diag, ext_data, lnd_diag, lnd_prog)

    TYPE(t_patch)             ,INTENT(IN)    :: p_patch
    TYPE(t_nh_diag)           ,INTENT(IN)    :: p_diag
    TYPE(t_external_data)     ,INTENT(IN)    :: ext_data
    TYPE(t_lnd_diag)          ,INTENT(IN)    :: lnd_diag
    TYPE(t_lnd_prog)          ,INTENT(INOUT) :: lnd_prog



    INTEGER :: jb, jt, jk, jc, ic              ! loop indices
    INTEGER :: nblks_c                         ! number of blocks
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startidx, i_endidx
    INTEGER :: ist
    INTEGER, PARAMETER :: nlev_soil_sma=6 ! equals the number of hydrologically active layers in TERRA


    REAL(wp) :: zsoil_hl(0:nlev_soil),smi_lim,rdp,wso_inc(6),scalfac_sma
    REAL(wp), DIMENSION(nproma) :: trh_avginc,plevap_pot,bsevap_pot,delta_h2o,smi_incint,smi_int
    REAL(wp), DIMENSION(nproma,nlev_soil) :: smi,smi_inc
  !-------------------------------------------------------------------------


    zsoil_hl(0) = 0
    zsoil_hl(1:nlev_soil) = depth_hl(1:nlev_soil)

    nblks_c   = p_patch%nblks_c
    rl_start  = 1
    rl_end    = min_rlcell

    ! Scaling factor for the integrated amount of soil water change, given a combined T-RH assimilation increment of 0.01
    ! Could be changed into a namelist parameter
    scalfac_sma = 6.e-3_wp ! = 6 mm H2O


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jt,jk,ic,jc,i_startidx,i_endidx,ist,smi,smi_lim,rdp,trh_avginc,plevap_pot,bsevap_pot,&
!$OMP            delta_h2o,smi_incint,wso_inc,smi_inc,smi_int)
    DO jb = 1, nblks_c


      CALL get_indices_c(p_patch, jb, 1, nblks_c, &
                         i_startidx, i_endidx, rl_start, rl_end)

      ! soil water on snow tiles is not touched
      DO jt = 1, ntiles_lnd

        smi_int(:) = 0._wp
        smi_incint(:) = 0._wp
        smi_inc(:,:) = 0._wp

!NEC$ ivdep
        DO ic = 1, ext_data%atm%lp_count_t(jb,jt)
          jc  = ext_data%atm%idx_lst_lp_t(ic,jb,jt)
          ist = ext_data%atm%soiltyp_t(jc,jb,jt)
          SELECT CASE(ist)
            CASE (3,4,5,6,7,8) ! soil types with non-zero water content

            rdp = MAX(1.e-3_wp,ext_data%atm%rootdp_t(jc,jb,jt))

            ! Calculate vertically integrated SMI, weighted by exponential root density profile assumed in TERRA
            ! Only the interval between the wilting point (SMI=0) and the field capacity (SMI=1) is taken into account
            DO jk= 1, nlev_soil_sma
              smi(jc,jk) = (lnd_prog%w_so_t(jc,jk,jb,jt)/dzsoil(jk)-cpwp(ist))/(cfcap(ist)-cpwp(ist))
              smi_lim = MIN(1._wp,MAX(0._wp,smi(jc,jk)))
              IF (rdp > zsoil_hl(jk-1)) THEN
                smi_int(jc) = smi_int(jc) + smi_lim*EXP(-3._wp/rdp*zml_soil(jk))*dzsoil(jk)
              ENDIF
            ENDDO

            ! potential for plant evaporation
            plevap_pot(jc) = smi_int(jc)*ext_data%atm%tai_t(jc,jb,jt)

            ! potential for base-soil evaporation
            bsevap_pot(jc) = (SUM(lnd_prog%w_so_t(jc,1:3,jb,jt))/SUM(dzsoil(1:3))-cadp(ist))* &
              ext_data%atm%eai_t(jc,jb,jt)/ext_data%atm%sai_t(jc,jb,jt)

            ! combined filtered assimilation increment for T and RH
            trh_avginc(jc) = &
              MERGE(p_diag%rh_avginc(jc,jb),MIN(0._wp,p_diag%rh_avginc(jc,jb)+0.01_wp),p_diag%rh_avginc(jc,jb)>0._wp) - &
              (0.025_wp + 0.05_wp*smi_int(jc))*(p_diag%t_avginc(jc,jb)-0.333_wp*p_diag%t_wgt_avginc(jc,jb))

            ! integrated soil water increment
            delta_h2o(jc) = trh_avginc(jc)*100._wp*scalfac_sma*dt_ana/86400._wp

            ! compute vertical distribution of soil water increments
            DO jk= 1, nlev_soil_sma

              ! part relevant for bare-soil evap: consider only the upper layers ...
              IF (jk <= 4) smi_inc(jc,jk) = 4.e-3_wp*bsevap_pot(jc)/MAX(0.02_wp,zml_soil(jk))**1.5_wp
              ! ... and suppress moistening the upper two soil layers in case of dry soil in order avoid a 
              ! short-lived impact that quickly disappears during the forecast
              IF (lnd_diag%snowfrac_lc_t(jc,jb,jt) > 0.25_wp .OR. &
                  trh_avginc(jc) > 0._wp .AND. smi(jc,3) <= 0.25_wp .AND. jk <= 2) smi_inc(jc,jk) = 0._wp

              ! consider plant evaporation in the root zone
              IF (lnd_diag%snowfrac_lc_t(jc,jb,jt) <= 0.25_wp .AND. lnd_prog%t_so_t(jc,jk+1,jb,jt) > tmelt &
                  .AND. ext_data%atm%rootdp_t(jc,jb,jt) > zsoil_hl(jk-1)) THEN
                smi_inc(jc,jk) = smi_inc(jc,jk)+plevap_pot(jc)*EXP(-2.5_wp/rdp*zml_soil(jk))
              ENDIF
              smi_incint(jc) = smi_incint(jc)+smi_inc(jc,jk)*dzsoil(jk)
            ENDDO

            ! apply soil water increments; drying below the wilting point and moistening above the 
            ! field capacity are suppressed
            DO jk= 1, nlev_soil_sma
              smi_inc(jc,jk) = smi_inc(jc,jk)/MAX(1.e-10_wp,smi_incint(jc))
              IF (smi(jc,jk) < 0._wp .AND. delta_h2o(jc) < 0._wp) smi_inc(jc,jk) = 0._wp
              IF (smi(jc,jk) > 1._wp .AND. delta_h2o(jc) > 0._wp) smi_inc(jc,jk) = 0._wp
              wso_inc(jk) = smi_inc(jc,jk)*dzsoil(jk)*delta_h2o(jc)
              lnd_prog%w_so_t(jc,jk,jb,jt) = MIN(dzsoil(jk)*cporv(ist),          &
                 MAX(lnd_prog%w_so_t(jc,jk,jb,jt)+wso_inc(jk), dzsoil(jk)*cadp(ist)) )
            ENDDO
          END SELECT

        ENDDO
      ENDDO

    ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE apply_sma

END MODULE mo_apt_routines

