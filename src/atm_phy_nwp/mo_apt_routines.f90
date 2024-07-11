!
!  Setup routines for adaptive parameter tuning and related computations
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
#include "omp_definitions.inc"
!----------------------------

MODULE mo_apt_routines

  USE mo_kind,                ONLY: wp
  USE mo_math_constants,      ONLY: rad2deg, pi2
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
  USE mo_nwp_lnd_types,       ONLY: t_wtr_prog, t_lnd_diag
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_state
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: min_rlcell_int
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_parallel_config,     ONLY: nproma
  USE mo_dynamics_config,     ONLY: nnow
  USE mo_time_config,         ONLY: time_config
  USE mo_grid_config,         ONLY: n_dom
  USE mo_lnd_nwp_config,      ONLY: ntiles_total, ntiles_water, ntiles_lnd, &
    &                               itype_canopy, itype_lndtbl, c_soil, c_soil_urb,  &
                                    lterra_urb, itype_eisa, cr_bsmin
  USE mo_satad,               ONLY: sat_pres_water, &  !! saturation vapor pressure w.r.t. water
    &                               spec_humi          !! Specific humidity
  USE mo_input_instructions,  ONLY: t_readInstructionListPtr, kInputSourceAna, kInputSourceAnaI

  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_cell

  USE mo_master_config,       ONLY: isRestart
  USE mo_initicon_types,      ONLY: t_initicon_state
  USE mo_initicon_config,     ONLY: icpl_da_sfcevap, dt_ana, icpl_da_snowalb, icpl_da_skinc,          &
                                    icpl_da_sfcfric, icpl_da_tkhmin, icpl_da_seaice, scalfac_da_sfcfric

  USE mo_nwp_tuning_config,   ONLY: itune_slopecorr



  IMPLICIT NONE

  PRIVATE


  PUBLIC  :: compute_filtincs, init_apt_fields


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
    TYPE(t_readInstructionListPtr),INTENT(INOUT) :: inputInstructions(:)

    TYPE(t_nh_diag) , POINTER :: p_diag
    TYPE(t_nh_prog),  POINTER :: p_prog_now

    INTEGER :: jg, jb, jc, nlev
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startidx, i_endidx, i_startblk, i_endblk


    REAL(wp) :: rh_inc, localtime_fac

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
              dt_ana/216000._wp*(initicon(jg)%atm_inc%temp(jc,nlev,jb)-p_diag%t_avginc(jc,jb))
          ENDDO
        ENDIF

        IF (icpl_da_sfcevap >= 2) THEN
          ! Calculate time-filtered RH assimilation increment at lowest model level (time scale 2.5 days);
          ! this serves as a proxy for the averaged RH2M bias
          DO jc = i_startidx, i_endidx
            rh_inc = initicon(jg)%atm_inc%qv(jc,nlev,jb)/ &
              spec_humi(sat_pres_water(p_diag%temp(jc,nlev,jb)),p_diag%pres_sfc(jc,jb))
            p_diag%rh_avginc(jc,jb) = p_diag%rh_avginc(jc,jb) + dt_ana/216000._wp*(rh_inc-p_diag%rh_avginc(jc,jb))
          ENDDO
        ENDIF

        IF (icpl_da_skinc >= 1) THEN
          ! weighted T assimilation increment; this serves as a proxy for the bias in diurnal temperature amplitude
          DO jc = i_startidx, i_endidx
            localtime_fac = COS(p_patch(jg)%cells%center(jc,jb)%lon + pi2/86400._wp *                           &
              (time_config%tc_exp_startdate%time%hour*3600._wp+time_config%tc_exp_startdate%time%minute*60._wp) )
            p_diag%t_wgt_avginc(jc,jb) = p_diag%t_wgt_avginc(jc,jb) + &
              dt_ana/216000._wp*(initicon(jg)%atm_inc%temp(jc,nlev,jb)*localtime_fac-p_diag%t_wgt_avginc(jc,jb))
          ENDDO
        ENDIF

        IF (icpl_da_sfcfric >= 1) THEN
          ! weighted wind speed increment for adaptive surface friction
          DO jc = i_startidx, i_endidx
            p_diag%vabs_avginc(jc,jb) = p_diag%vabs_avginc(jc,jb) + dt_ana/216000._wp * (         &
              SQRT( (p_diag%u(jc,nlev,jb)+initicon(jg)%atm_inc%u(jc,nlev,jb))**2 +                &
                    (p_diag%v(jc,nlev,jb)+initicon(jg)%atm_inc%v(jc,nlev,jb))**2 ) -              &
              SQRT(p_diag%u(jc,nlev,jb)**2 + p_diag%v(jc,nlev,jb)**2) - p_diag%vabs_avginc(jc,jb) )
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


    REAL(wp) :: dtfac_heatc, tbias_wgt, zlat, zlon, slope(nproma), skinc_fac, trh_bias, scal

    ! adaptation factor to analysis interval for adaptive heat conductivity/capacity and skin conductivity
    dtfac_heatc = (10800._wp/dt_ana)**(2._wp/3._wp)

    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jt,i_startidx,i_endidx,i_count,ilu,tbias_wgt,zlat,zlon,slope,skinc_fac,trh_bias,scal)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

 
      ! tuning factor for rlam_heat depending on skin conductivity and analyzed T2M/RH2M bias
      IF (itype_canopy == 2 .AND. icpl_da_sfcevap >= 3) THEN
        DO jt = 1, ntiles_total + ntiles_water
          DO jc = i_startidx,i_endidx
            IF (jt <= ntiles_lnd) THEN ! snow-free land points
              prm_diag%rlamh_fac_t(jc,jb,jt) = 1._wp - 0.9_wp*MAX(0._wp, MIN(1._wp,                              &
                2.5_wp*(10800._wp/dt_ana*(100._wp*p_diag%rh_avginc(jc,jb)-4._wp*p_diag%t_avginc(jc,jb))-0.4_wp) ))
            ELSE IF (jt <= ntiles_total) THEN ! snow-covered land points
              prm_diag%rlamh_fac_t(jc,jb,jt) = 1._wp - 0.9_wp*MAX(0._wp, MIN(1._wp, &
                2.5_wp*(10800._wp/dt_ana*(MAX(0._wp,100._wp*p_diag%rh_avginc(jc,jb))-4._wp*p_diag%t_avginc(jc,jb))-0.4_wp) ))
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
              prm_diag%rlamh_fac_t(jc,jb,jt) =                                                                              &
                1._wp - 0.9_wp*MAX(0._wp,MIN(1._wp,(60._wp-ext_data%atm%skinc_t(jc,jb,jt))/30._wp)) *                       &
                MAX(0._wp,MIN(1._wp,2.5_wp*(p_diag%t2m_bias(jc,jb)+100._wp*10800._wp/dt_ana*p_diag%rh_avginc(jc,jb)-0.4_wp)))
            ELSE IF (jt == ntiles_total + ntiles_water) THEN ! seaice points
              prm_diag%rlamh_fac_t(jc,jb,jt) = 0.25_wp
            ELSE
              prm_diag%rlamh_fac_t(jc,jb,jt) = 1._wp
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      IF (icpl_da_snowalb >= 1 .AND. .NOT. isRestart()) THEN
        ! Tuning factor for snow albedo
        DO jc = i_startidx,i_endidx
          IF (ANY(p_diag_lnd%h_snow_t(jc,jb,1:ntiles_total) > 0._wp) .OR. p_prog_wtr_now%h_ice(jc,jb) > 0._wp) THEN
            IF (p_diag%t_avginc(jc,jb) > 0._wp) THEN
              prm_diag%snowalb_fac(jc,jb) = MAX(0.75_wp,1._wp/(1._wp+10800._wp/dt_ana*0.8_wp*p_diag%t_avginc(jc,jb)))
            ELSE
              prm_diag%snowalb_fac(jc,jb) = MIN(4._wp/3._wp,1._wp-10800._wp/dt_ana*0.8_wp*p_diag%t_avginc(jc,jb))
            ENDIF
          ENDIF
          IF (icpl_da_snowalb >= 2) THEN ! albedo factor is also applied to sea ice and needs to be restricted to the vicinity of land
            tbias_wgt = MIN(1._wp,100._wp*ext_data%atm%fr_land_smt(jc,jb))
            prm_diag%snowalb_fac(jc,jb) = tbias_wgt*prm_diag%snowalb_fac(jc,jb) + (1._wp-tbias_wgt)
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
          tbias_wgt = 10800._wp/dt_ana*(p_diag%t_avginc(jc,jb)+0.5_wp*p_diag%t_wgt_avginc(jc,jb))
          IF (tbias_wgt < 0._wp) THEN
            prm_diag%tkred_sfc_h(jc,jb) = MAX(0.25_wp, 1._wp+2._wp*tbias_wgt)
          ELSE
            prm_diag%tkred_sfc_h(jc,jb) = 1._wp/SQRT(MAX(0.25_wp, 1._wp-2._wp*tbias_wgt))
          ENDIF
        ENDDO
      ENDIF
      IF (icpl_da_sfcfric >= 1) THEN
        ! Tuning factor for surface friction (roughness length and SSO blocking)
        DO jc = i_startidx,i_endidx
          IF (p_diag%vabs_avginc(jc,jb) > 0._wp) THEN
            prm_diag%sfcfric_fac(jc,jb) = MAX(0.25_wp, 1._wp-scalfac_da_sfcfric*10800._wp/dt_ana*p_diag%vabs_avginc(jc,jb))
          ELSE
            prm_diag%sfcfric_fac(jc,jb) = 1._wp/MAX(0.25_wp, 1._wp+scalfac_da_sfcfric*10800._wp/dt_ana*p_diag%vabs_avginc(jc,jb))
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
      IF (itune_slopecorr >= 1) THEN
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

      DO jt = 1, ntiles_total
        i_count = ext_data%atm%lp_count_t(jb,jt)
        IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty

        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_lp_t(ic,jb,jt)

          ilu = ext_data%atm%lc_class_t(jc,jb,jt)


          IF (icpl_da_sfcevap >= 3 .AND. icpl_da_skinc >= 1) THEN
            trh_bias = 10800._wp/dt_ana*(125._wp*p_diag%rh_avginc(jc,jb) - &
              4._wp*(p_diag%t_avginc(jc,jb)-0.5_wp*p_diag%t_wgt_avginc(jc,jb)))
          ELSE IF (icpl_da_sfcevap >= 3) THEN
            trh_bias = 10800._wp/dt_ana*(100._wp*p_diag%rh_avginc(jc,jb)-4._wp*p_diag%t_avginc(jc,jb))
          ELSE IF (icpl_da_sfcevap >= 2) THEN
            trh_bias = p_diag%t2m_bias(jc,jb) + 100._wp*10800._wp/dt_ana*p_diag%rh_avginc(jc,jb)
          ELSE IF (icpl_da_sfcevap == 1) THEN
            trh_bias = p_diag%t2m_bias(jc,jb)
          ELSE
            trh_bias = 0._wp
          ENDIF

          IF (icpl_da_sfcevap >= 4) THEN
            IF (trh_bias < 0._wp) THEN
              ext_data%atm%rsmin2d_t(jc,jb,jt) = ext_data%atm%stomresmin_lcc(ilu)*(1._wp-0.75_wp*trh_bias)
              ext_data%atm%r_bsmin(jc,jb)      = cr_bsmin*(1._wp-trh_bias)
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

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE init_apt_fields


END MODULE mo_apt_routines

