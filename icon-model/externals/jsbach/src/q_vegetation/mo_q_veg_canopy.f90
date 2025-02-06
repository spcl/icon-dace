!> QUINCY calculate canopy layer properties
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
!>#### calculate canopy layer properties, e.g., C & N content, LAI, fractional allocation
!>
MODULE mo_q_veg_canopy
#ifndef __NO_QUINCY__

  USE mo_kind,                        ONLY: wp
  USE mo_jsb_impl_constants,          ONLY: test_false_true
  USE mo_exception,                   ONLY: message_text, finish
  USE mo_q_assimi_process,            ONLY: calc_photosynthesis

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_canopy_layers, calc_marginal_canopy_flux_increment, calc_dir_optimal_cn_leaf

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_veg_canopy'

CONTAINS

  ! ======================================================================================================= !
  !>calc canopy layers - update the canopy layer content
  !>
  !>  LAI, cummulative LAI above layer-mid-point
  !>  foliar N, fractional allocation to photosynthetic functioning
  !>
  SUBROUTINE calc_canopy_layers( &
    & nc, &
    & ncanopy, &
    & dtime, &
    & dz, &
    & lbounds, &
    & ubounds, &
    & lctlib_ps_pathway, &
    & lctlib_k0_fn_struc, &
    & lctlib_fn_oth_min, &
    & lctlib_sla, &
    & lctlib_np_leaf, &
    & lctlib_gmin, &
    & lctlib_g0, &
    & lctlib_g1, &
    & lctlib_t_jmax_omega, &
    & flag_optimal_Nfraction, &
    & canopy_cond_scheme, & ! Q_ASSIMI_ config
    & leaf_nitrogen, &
    & lai, &
    & ppfd_sunlit_cl, &
    & ppfd_shaded_cl, &
    & fleaf_sunlit_cl, &
    & fn_rub_cl, &
    & fn_et_cl, &
    & fn_pepc_cl, &
    & fn_chl_cl, &
    & fn_oth_cl, &
    & lai_cl, &
    & cumm_lai_cl, &
    & leaf_nitrogen_cl, &
    & t_air, &
    & t_acclim, &
    & press_srf, &
    & co2_mixing_ratio, &
    & aerodyn_cond, &
    & beta_air, &
    & beta_soa, &
    & beta_soil_ps, &
    & beta_sinklim_ps, &
    & beta_soil_gs, &
    & t_jmax_opt, &
    & inquire_n_fractions)

    USE mo_veg_constants,           ONLY: carbon_per_dryweight_leaf, &
      &                                   jmax2vcmax_C4, k0_fn_chl_C4, k1_fn_chl_C4, k0_fn_pepc_C4, &
      &                                   jmax2vcmax_C3, k0_fn_chl_C3, k1_fn_chl_C3, k0_fn_pepc_C3, kfn_chl
    USE mo_q_assimi_parameters,     ONLY: kn, chl2n, jmax2n, k1_fn_struc, vcmax2n
    USE mo_q_assimi_constants,      ONLY: ic4phot
    USE mo_jsb_physical_constants,  ONLY: molar_mass_C, molar_mass_N, molar_mass_P
    USE mo_jsb_math_constants,      ONLY: eps8, eps12
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                                INTENT(in)              :: nc                    !< dimensions
    INTEGER,                                INTENT(in)              :: ncanopy               !< number of canopy layers
    REAL(wp),                               INTENT(in)              :: dtime                 !< timestep length
    REAL(wp), DIMENSION(ncanopy),           INTENT(in)              :: dz                    !< canopy layer thickness
    REAL(wp), DIMENSION(ncanopy),           INTENT(in)              :: lbounds               !< lower bound of canopy layer (lbounds < ubounds; lbounds is at top of canopy layer)
    REAL(wp), DIMENSION(ncanopy),           INTENT(in)              :: ubounds               !< upper bound of canopy layer (lbounds < ubounds; ubounds is at bottom of canopy layer)
    INTEGER,                                INTENT(in)              :: lctlib_ps_pathway     !< dimensions
    REAL(wp),                               INTENT(in)              :: lctlib_k0_fn_struc    !< lctlib parameter
    REAL(wp),                               INTENT(in)              :: lctlib_fn_oth_min     !< lctlib parameter
    REAL(wp),                               INTENT(in)              :: lctlib_sla            !< lctlib parameter
    REAL(wp),                               INTENT(in)              :: lctlib_np_leaf        !< lctlib parameter
    REAL(wp),                               INTENT(in)              :: lctlib_gmin           !< lctlib parameter
    REAL(wp),                               INTENT(in)              :: lctlib_g0             !< lctlib parameter
    REAL(wp),                               INTENT(in)              :: lctlib_g1             !< lctlib parameter
    REAL(wp),                               INTENT(in)              :: lctlib_t_jmax_omega   !< lctlib parameter
    LOGICAL,                                INTENT(in)              :: flag_optimal_Nfraction!< on/off optimise leaf internal N allocation
    CHARACTER(len=*),                       INTENT(in)              :: canopy_cond_scheme    !< canopy_conductance_scheme: medlyn / ballberry
    REAL(wp), DIMENSION(nc),                INTENT(in)              :: leaf_nitrogen         !< total canopy N (mol/m2)
    REAL(wp), DIMENSION(nc),                INTENT(in)              :: lai                   !< total leaf area index
    REAL(wp), DIMENSION(nc, ncanopy),       INTENT(inout)           :: ppfd_sunlit_cl        !< absorbed Photosynthetically Active Photon Flux Density on sunlit leafes (average over time)
    REAL(wp), DIMENSION(nc, ncanopy),       INTENT(inout)           :: ppfd_shaded_cl        !< absorbed Photosynthetically Active Photon Flux Density on shaded leafes (average over time)
    REAL(wp), DIMENSION(nc, ncanopy),       INTENT(inout)           :: fleaf_sunlit_cl       !< fraction of sunlit leaves (average over time)
    REAL(wp), DIMENSION(nc, ncanopy),       INTENT(inout)           :: fn_rub_cl             !< fraction of N in Rubisco
    REAL(wp), DIMENSION(nc, ncanopy),       INTENT(inout)           :: fn_et_cl              !< fraction of N in Electron Transport
    REAL(wp), DIMENSION(nc, ncanopy),       INTENT(inout)           :: fn_pepc_cl            !< fraction of N in PepC
    REAL(wp), DIMENSION(nc, ncanopy),       INTENT(inout)           :: fn_chl_cl             !< fraction of N in Chloroplast
    REAL(wp), DIMENSION(nc, ncanopy),       INTENT(inout)           :: fn_oth_cl             !< fraction of N not associated with PS
    REAL(wp), DIMENSION(nc, ncanopy),       INTENT(out)             :: lai_cl                !< LAI of canopy layer
    REAL(wp), DIMENSION(nc, ncanopy),       INTENT(out)             :: cumm_lai_cl           !< cummlative LAI above canopy layer
    REAL(wp), DIMENSION(nc, ncanopy),       INTENT(out)             :: leaf_nitrogen_cl      !< total N in canopy layer (mmol N/m2)
    REAL(wp), DIMENSION(nc),                INTENT(in),   OPTIONAL  :: t_air                 !< air temperature (avg)
    REAL(wp), DIMENSION(nc),                INTENT(in),   OPTIONAL  :: t_acclim              !< acclimation temperature (avg)
    REAL(wp), DIMENSION(nc),                INTENT(in),   OPTIONAL  :: press_srf             !< air pressure (avg)
    REAL(wp), DIMENSION(nc),                INTENT(in),   OPTIONAL  :: co2_mixing_ratio      !< CO2 concentration (avg)
    REAL(wp), DIMENSION(nc),                INTENT(in),   OPTIONAL  :: aerodyn_cond          !< aerodynamic conductance (avg)
    REAL(wp), DIMENSION(nc),                INTENT(in),   OPTIONAL  :: beta_air              !< beta air (avg)
    REAL(wp), DIMENSION(nc),                INTENT(in),   OPTIONAL  :: beta_soa              !< beta soa (avg)
    REAL(wp), DIMENSION(nc),                INTENT(in),   OPTIONAL  :: beta_soil_ps          !< beta soil ps (avg)
    REAL(wp), DIMENSION(nc),                INTENT(in),   OPTIONAL  :: beta_sinklim_ps       !< beta sinklim ps (avg)
    REAL(wp), DIMENSION(nc),                INTENT(in),   OPTIONAL  :: beta_soil_gs          !< beta soil gs (avg)
    REAL(wp), DIMENSION(nc),                INTENT(in),   OPTIONAL  :: t_jmax_opt            !< Jmax optimum temperature (avg)
    LOGICAL,                                INTENT(in),   OPTIONAL  :: inquire_n_fractions   !< whether to return fractions or also run optimisation
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), DIMENSION(nc)         :: leaf_nitrogen_tc                              ! top-of-the-canopy N content
    REAL(wp)                        :: k0_fn_chl, k1_fn_chl                          ! factors in chlorophyll profile calculation
    REAL(wp)                        :: jmax2vcmax                                    ! jmax to vcmax ratio
    REAL(wp)                        :: f_pepc                                        ! fraction of N in PEP Carboxylase
    LOGICAL                         :: optimise_n_fractions
    INTEGER                         :: icanopy                                       ! looping over ncanopy and nc
    REAL(wp), DIMENSION(nc,ncanopy) :: mid_point_lai_cl
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_canopy_layers'
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.9 initialise output (to ensure that 0.0 is returned in case of the below statements is TRUE)
    !>
    lai_cl(:,:)           = 0._wp
    cumm_lai_cl(:,:)      = 0._wp
    leaf_nitrogen_cl(:,:) = 0._wp
    leaf_nitrogen_tc(:)   = 0._wp
    mid_point_lai_cl(:,:) = 0._wp

    !>  0.9.1  set pathway specific constants
    !>
    IF(lctlib_ps_pathway == ic4phot) THEN
        jmax2vcmax = jmax2vcmax_C4
        k0_fn_chl  = k0_fn_chl_C4
        k1_fn_chl  = k1_fn_chl_C4
        f_pepc     = k0_fn_pepc_C4
    ELSE
        jmax2vcmax = jmax2vcmax_C3
        k0_fn_chl  = k0_fn_chl_C3
        k1_fn_chl  = k1_fn_chl_C3
        f_pepc     = k0_fn_pepc_C3
    ENDIF

    !>  0.9.2 configure options for fractions optimisation
    !>
    optimise_n_fractions = flag_optimal_Nfraction
    IF (PRESENT(inquire_n_fractions)) optimise_n_fractions = .NOT.inquire_n_fractions

    !>  0.9.3  check presence of all arguments that are OPTIONAL but needed for 'CALL optimise_canopy_n_fractions()' \n
    !>         only in case optimise_n_fractions = TRUE \n
    !>         note, all variables are of DIMENSION(nc)         but compiler gives an error if '(:)' is added here
    !>
    IF(optimise_n_fractions)THEN
       IF(.NOT.PRESENT(t_air)            .OR..NOT.PRESENT(t_acclim)         .OR..NOT.PRESENT(press_srf)    .OR. &
          .NOT.PRESENT(co2_mixing_ratio) .OR..NOT.PRESENT(aerodyn_cond)     .OR. &
          .NOT.PRESENT(beta_air)         .OR..NOT.PRESENT(beta_soa)         .OR. &
          .NOT.PRESENT(beta_soil_ps)     .OR..NOT.PRESENT(beta_sinklim_ps)  .OR..NOT.PRESENT(beta_soil_gs) .OR. &
          .NOT.PRESENT(t_jmax_opt)) THEN

          WRITE(message_text,'(a)') 'called to optimise N fractions, but without [climate] variables'
          CALL finish("mo_q_veg_canopy:calc_canopy_layer",message_text)
       ENDIF
    ENDIF


    !>1.0 calculate top of the canopy N
    !>
    ! dz(1) == ubounds(1) -> thickness of uppermost canopy layer
    WHERE(lai(:) > dz(1))
       leaf_nitrogen_tc(:) = leaf_nitrogen(:) * 1000._wp * kn / (1._wp - exp(-kn * lai(:)))
    ELSEWHERE(lai(:) > eps8)  ! if LAI does not fill the first canopy layer, top canopy N = leaf N
       leaf_nitrogen_tc(:) = leaf_nitrogen(:) * 1000._wp /lai(:)
    ENDWHERE


    !>2.0 calculate LAI and N content of canopy layers
    !>

    !>  2.1 LAI of canopy layer
    !>
    DO icanopy = 1,ncanopy
      WHERE(lai(:) > eps8) ! avoid any calculation if there are no leaves
        ! LAI calculation
        WHERE(ubounds(icanopy) < lai(:) .OR. ABS(ubounds(icanopy) - lai(:)) < eps12 )
          ! ubounds contains cummulative LAI at the bottom of this layer
          ! if canopy layer is fully developed
          lai_cl(:,icanopy) = dz(icanopy)
        ELSEWHERE(lbounds(icanopy) < lai(:) .AND. ubounds(icanopy) > lai(:))
          ! if canopy layer is only partially filled assign difference between full LAI and lc to canopy layer
          lai_cl(:,icanopy) = lai(:) - lbounds(icanopy)
        ENDWHERE

        ! make sure all LAI is represented by letting last layer go to the end of the canopy
        WHERE(icanopy==ncanopy .AND. ubounds(ncanopy) < lai(:) )
          lai_cl(:,ncanopy) = lai(:) - lbounds(ncanopy)
        ENDWHERE

        ! LAI mid-point of canopy layer
        WHERE(lai_cl(:,icanopy) > 0.0_wp)
          mid_point_lai_cl(:,icanopy) = (SUM(lai_cl(:,1:icanopy), DIM=2) - lbounds(icanopy)) / 2.0_wp
        ENDWHERE
      ENDWHERE ! lai(:) > eps8
    END DO   ! loop ncanopy

    !>  2.2 cummulative LAI above mid-point of canopy layer
    !>
    ! 1st layer
    WHERE(lai(:) > eps8) ! avoid any calculation if there are no leaves
      cumm_lai_cl(:,1) = mid_point_lai_cl(:,1)
    ENDWHERE
    ! 2-x layers
    DO icanopy = 2,ncanopy
      WHERE(lai(:) > eps8) ! avoid any calculation if there are no leaves
        WHERE(lai_cl(:,icanopy) > eps8)
          cumm_lai_cl(:,icanopy) = lbounds(icanopy) + mid_point_lai_cl(:,icanopy)
        ENDWHERE
      ENDWHERE
    ENDDO

    !>  2.3 N content of the canopy layer
    !>
    DO icanopy = 1,ncanopy
      WHERE(lai(:) > eps8)                        ! avoid any calculation if there are no leaves
        WHERE(lai(:) > ubounds(1))                ! LAI is present and LAI fills more than first layer
          WHERE(cumm_lai_cl(:,icanopy) > eps12)   ! this layer has some leafs
            leaf_nitrogen_cl(:,icanopy) = leaf_nitrogen_tc(:) * exp(-kn * cumm_lai_cl(:,icanopy))
          ELSEWHERE
            leaf_nitrogen_cl(:,icanopy) = 0.0_wp  ! set N to zero if layer has no leafs
          ENDWHERE
        ELSEWHERE
          leaf_nitrogen_cl(:,icanopy)   = 0.0_wp  ! set N to zero if layer has no leafs
        ENDWHERE
      ENDWHERE
    ENDDO
    WHERE(lai(:) > eps8 .AND..NOT. lai(:) > ubounds(1))   ! LAI is present but 1st layer is not filled
      leaf_nitrogen_cl(:,1) = leaf_nitrogen_tc(:)         ! set 1st layer N to leaf_nitrogen_tc, but leave all other layers with zero N
    ENDWHERE


    !>3.0 fractionation of N within canopy laer
    !>

    !>  3.1 fraction not associated with photosynthesis
    !>
    WHERE(leaf_nitrogen_cl(:,:) > 0.0_wp)
       fn_oth_cl(:,:) = lctlib_k0_fn_struc - k1_fn_struc * leaf_nitrogen_cl(:,:)
    ELSEWHERE
       fn_oth_cl(:,:) = 1.0_wp
    END WHERE
    ! structural limitation to the indefinite increase in phot. N described by above equation
    WHERE (fn_oth_cl(:,:) < lctlib_fn_oth_min) fn_oth_cl(:,:) = lctlib_fn_oth_min


    !>  3.2 photosynthetic fractions
    !>
    ! prescribed, or dynamically calculated if flag_optimal_Nfraction = .TRUE.
    IF (.NOT.flag_optimal_Nfraction) THEN
        !> 3.2.1 case of presribed fractions
        WHERE(lai_cl(:,:) > eps8) ! canopy layers with leaves
           ! fraction in chlorophyll
           fn_chl_cl(:,:) = MIN(0.99_wp - fn_oth_cl(:,:), &
                                (k0_fn_chl - k1_fn_chl * exp(-kfn_chl * cumm_lai_cl(:,:))) / chl2n)
           ! fraction in PepC
           fn_pepc_cl(:,:) = f_pepc
           ! fraction of electron-transport is dependent on N in electron transport versus N in light harvesting
           fn_et_cl(:,:) = (1._wp - (fn_oth_cl(:,:) + fn_chl_cl(:,:) + fn_pepc_cl(:,:) )) / &
                           (1._wp + jmax2n / (jmax2vcmax * vcmax2n))
           ! fraction of Rubisco is determined as the remainder
           fn_rub_cl(:,:) = 1._wp - fn_oth_cl(:,:) - fn_chl_cl(:,:) - fn_et_cl(:,:) - fn_pepc_cl(:,:)

           ! ! convert mmol N / m2 to mg N / g DW
           ! hlp3(:,:) = leaf_nitrogen_cl(:,:) * molar_mass_N / &
           !               ((molar_mass_C / lctlib_sla)/carbon_per_dryweight_leaf)
           ! hlp4(:,:) = leaf_nitrogen_cl(:,:) / lctlib_np_leaf * molar_mass_P / &
           !               ((molar_mass_C / lctlib_sla)/carbon_per_dryweight_leaf)
           ! !DE stuff Vmax = exp(k1 + k2 * log(P) + k3 * log(N))
           ! ! implied Vmax and Jmax in nmol / g DW / s
           ! !hlp1(:,:) = exp(4.44904_wp + 0.3472_wp*log10(hlp4(:,:)) + 0.49078*log10(hlp3(:,:)))
           ! !hlp2(:,:) = exp(5.49435_wp + 0.37345_wp*log10(hlp4(:,:)) + 0.41435_wp*log10(hlp3(:,:)))
           ! ! convert to f_rub_cl: Vmax = vmax2n * frub * leaf_n | micro-mol / m2 / s
           ! !                      Vmax1 = vmax2n * frub * leaf_n / DW * 1000 | nmol / gDW / s
           ! !                      exp(X) = vmax2n * frub * leaf_n
           ! !                      frub = exp(X) / leaf_n / vmax2n
           ! !fn_rub_cl(:,:) = exp(4.44904_wp + 0.3472_wp*log(hlp4(:,:)) + 0.49078*log(hlp3(:,:))) / &
           ! !                 (leaf_nitrogen_cl(:,:) * vcmax2n * 1000._wp / (molar_mass_C / lctlib_sla / 0.48_wp))
           ! !fn_et_cl(:,:) = exp(5.49435_wp + 0.37345_wp*log(hlp4(:,:)) + 0.41435_wp*log(hlp3(:,:))) / 4.0_wp / &
           ! !                 (leaf_nitrogen_cl(:,:) * jmax2n * 1000._wp / (molar_mass_C / lctlib_sla / 0.48_wp))
           ! !fn_oth_cl(:,:) = 1.0_wp - fn_chl_cl(:,:) - fn_et_cl(:,:) - fn_pepc_cl(:,:) - fn_rub_cl(:,:)

        ELSEWHERE ! canopy layers without leaves
           leaf_nitrogen_cl(:,:)  = 0.0_wp
           fn_oth_cl(:,:)         = 1.0_wp
           fn_chl_cl(:,:)         = 0.0_wp
           fn_et_cl(:,:)          = 0.0_wp
           fn_rub_cl(:,:)         = 0.0_wp
           fn_pepc_cl(:,:)        = 0.0_wp
        ENDWHERE

        ! for canopy layer > 1
        DO icanopy = 2,ncanopy
           WHERE(lai_cl(:,icanopy) > eps8)
              WHERE(ppfd_sunlit_cl(:,icanopy) < eps8) ! no radiation, take from above
                 ppfd_sunlit_cl(:,icanopy)   = ppfd_sunlit_cl(:,icanopy-1)
                 ppfd_shaded_cl(:,icanopy)   = ppfd_shaded_cl(:,icanopy-1)
                 fleaf_sunlit_cl(:,icanopy)  = fleaf_sunlit_cl(:,icanopy-1)
              ENDWHERE
           ENDWHERE
        ENDDO

    ELSE ! IF (.NOT.flag_optimal_Nfraction)
        !>    3.2.1 optimal N fraction case
        !>
        ! for canopy layer = 1
        icanopy = 1
        WHERE(lai_cl(:,icanopy) > eps8)
          WHERE(fn_chl_cl(:,icanopy) < eps8) ! new layer with new leaves
            ! first layers: use static allocation as first guess
            fn_chl_cl(:,icanopy)  = (k0_fn_chl - k1_fn_chl * exp(-kfn_chl * cumm_lai_cl(:,icanopy))) / chl2n
            ! fraction in PepC
            fn_pepc_cl(:,icanopy) = f_pepc
            ! fraction of electron-transport is dependent on N in electron transport versus N in light harvesting
            fn_et_cl(:,icanopy)   = (1._wp - (fn_oth_cl(:,icanopy) + fn_chl_cl(:,icanopy) + fn_pepc_cl(:,icanopy) )) / &
                                    (1._wp + jmax2n/(jmax2vcmax * vcmax2n))
            ! fraction of Rubisco is determined as the remainder
            fn_rub_cl(:,icanopy)  = 1._wp - fn_oth_cl(:,icanopy) - fn_chl_cl(:,icanopy) - fn_et_cl(:,icanopy) - &
                                      fn_pepc_cl(:,icanopy)
          ENDWHERE
        ELSEWHERE
              leaf_nitrogen_cl(:,icanopy) = 0.0_wp
              fn_oth_cl(:,icanopy)        = 1.0_wp
              fn_chl_cl(:,icanopy)        = 0.0_wp
              fn_et_cl(:,icanopy)         = 0.0_wp
              fn_rub_cl(:,icanopy)        = 0.0_wp
              fn_pepc_cl(:,icanopy)       = 0.0_wp
        ENDWHERE

        ! for canopy layer > 1
        DO icanopy = 2,ncanopy
           WHERE(lai_cl(:,icanopy) > eps8)
              WHERE(fn_chl_cl(:,icanopy) < eps8) ! other layers: copy from above layer
                 fn_chl_cl(:,icanopy)        = fn_chl_cl(:,icanopy-1)
                 fn_et_cl(:,icanopy)         = fn_et_cl(:,icanopy-1)
                 fn_rub_cl(:,icanopy)        = fn_rub_cl(:,icanopy-1)
                 fn_pepc_cl(:,icanopy)       = fn_pepc_cl(:,icanopy-1)
                 ppfd_sunlit_cl(:,icanopy)   = ppfd_sunlit_cl(:,icanopy-1)
                 ppfd_shaded_cl(:,icanopy)   = ppfd_shaded_cl(:,icanopy-1)
                 fleaf_sunlit_cl(:,icanopy)  = fleaf_sunlit_cl(:,icanopy-1)
              ENDWHERE
           ELSEWHERE
              leaf_nitrogen_cl(:,icanopy) = 0.0_wp
              fn_oth_cl(:,icanopy)        = 1.0_wp
              fn_chl_cl(:,icanopy)        = 0.0_wp
              fn_et_cl(:,icanopy)         = 0.0_wp
              fn_rub_cl(:,icanopy)        = 0.0_wp
              fn_pepc_cl(:,icanopy)       = 0.0_wp
           ENDWHERE
        ENDDO ! icanopy

        IF(optimise_n_fractions) THEN
          DO icanopy = 1,ncanopy
            CALL optimise_canopy_n_fractions(  dtime                      , &
                                               lctlib_gmin                , & ! lctlib
                                               lctlib_g0                  , &
                                               lctlib_g1                  , &
                                               lctlib_t_jmax_omega        , &
                                               lctlib_ps_pathway          , & ! lctlib
                                               TRIM(canopy_cond_scheme)   , & ! Q_ASSIMI_ config
                                               lai_cl(:,icanopy)          , & ! in
                                               leaf_nitrogen_cl(:,icanopy), &
                                               ppfd_sunlit_cl(:,icanopy)  , &
                                               ppfd_shaded_cl(:,icanopy)  , &
                                               fleaf_sunlit_cl(:,icanopy) , &
                                               t_air(:)                   , &
                                               t_acclim(:)                , &
                                               press_srf(:)               , &
                                               co2_mixing_ratio(:)        , &
                                               aerodyn_cond(:)            , &
                                               beta_air(:)                , &
                                               beta_soa(:)                , &
                                               beta_soil_ps(:)            , &
                                               beta_sinklim_ps(:)         , &
                                               beta_soil_gs(:)            , &
                                               t_jmax_opt(:)              , &
                                               fn_oth_cl(:,icanopy)       , &
                                               fn_chl_cl(:,icanopy)       , & ! inout
                                               fn_rub_cl(:,icanopy)       , &
                                               fn_et_cl(:,icanopy)        , &
                                               fn_pepc_cl(:,icanopy)      )
          ENDDO
        ENDIF ! optimise_n_fractions
    ENDIF     ! end ELSE of IF (.NOT.flag_optimal_Nfraction)

  END SUBROUTINE calc_canopy_layers


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to calc_canopy_layers
  !
  ! --------------------------------------------------------------------------------------------------
  !> Subroutine to opimise fraction of leaf N allocated to each photosynthetic component
  !! Calculates carbon assimilation (an_cl) by varying fractions in either directions and determines
  !! fractions that give maximum assimilation
  !-------------------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE optimise_canopy_n_fractions( &
                                    dtime                      , &
                                    gmin                      , &   ! lctlib
                                    g0                        , &
                                    g1                        , &
                                    t_jmax_omega              , &
                                    ps_pathway                , &   ! lctlib
                                    canopy_cond_scheme        , &   ! Q_ASSIMI_ config
                                    lai_cl                    , &
                                    leaf_nitrogen_cl          , &
                                    ppfd_sunlit_cl            , &
                                    ppfd_shaded_cl            , &
                                    fleaf_sunlit_cl           , &
                                    t_air                     , &
                                    t_acclim                  , &
                                    press_srf                 , &
                                    co2_mixing_ratio          , &
                                    aerodyn_cond              , &
                                    beta_air                  , &
                                    beta_soa                  , &
                                    beta_soil_ps              , &
                                    beta_sinklim_ps           , &
                                    beta_soil_gs              , &
                                    t_jmax_opt                , &
                                    fn_oth_cl                 , &
                                    fn_chl_cl                 , &
                                    fn_rub_cl                 , &
                                    fn_et_cl                  , &
                                    fn_pepc_cl)

    USE mo_jsb_math_constants,  ONLY: one_day, eps8
    USE mo_veg_constants,       ONLY: delta_n_fraction


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp),                 INTENT(in)    :: dtime                  !< timestep length
    REAL(wp),                 INTENT(in)    :: gmin            , &    !< land-cover-type library parameter
                                               g0              , &    !< land-cover-type library parameter
                                               g1              , &    !< land-cover-type library parameter
                                               t_jmax_omega           !< land-cover-type library parameter
    INTEGER,                  INTENT(in)    :: ps_pathway             !< land-cover-type library parameter
    CHARACTER(len=*),         INTENT(IN)    :: canopy_cond_scheme     !< canopy_conductance_scheme: medlyn / ballberry
    REAL(wp),                 INTENT(in)    :: lai_cl
    REAL(wp),                 INTENT(in)    :: leaf_nitrogen_cl, &    !< total canopy N content
                                               ppfd_sunlit_cl,   &    !< Photosynthetically Active Photon Flux Density
                                                                      !!  on sunlit leaves (day-time running mean, muemol/ms/s)
                                               ppfd_shaded_cl,   &    !< Photosynthetically Active Photon Flux Density
                                                                      !!  on shaded leaves (day-time running mean, muemol/ms/s)
                                               fleaf_sunlit_cl,  &    !< fraction of sunlit leaves (day-time running mean, muemol/ms/s)
                                               t_air,            &    !< temperature (day-time running mean, K)
                                               t_acclim,         &    !< ..
                                               press_srf,        &    !< air pressure (day-time running mean, Pa)
                                               co2_mixing_ratio, &    !< CO2 mixing ratio (day-time running mean, ppm)
                                               aerodyn_cond,     &    !< ga -- aerodynamic conductance (day-time running mean, m/s)
                                               beta_air,         &    !< beta air (day-time running mean, --)
                                               beta_soa,         &    !< beta soa (day-time running mean, --)
                                               beta_soil_ps,     &    !< beta soil for PS (day-time running mean, --)
                                               beta_sinklim_ps,  &    !< beta sinklim for PS (day-time running mean, --)
                                               beta_soil_gs,     &    !< beta soil for GC (day-time running mean, --)
                                               t_jmax_opt             !< optimum temperature for Jmax (day-time running mean, degreeC)
    REAL(wp),                 INTENT(in)    :: fn_oth_cl              !< fraction of leaf N not associated with PS
    REAL(wp),                 INTENT(inout) :: fn_chl_cl,        &    !< fraction of leaf N in chlorophyll
                                               fn_rub_cl,        &    !< fraction of leaf N in rubisco
                                               fn_et_cl,         &    !< fraction of leaf N in electron transport
                                               fn_pepc_cl             !< fraction of leaf N in electron transport
    ! ---------------------------
    ! 0.2 Local
    REAL(wp) :: hlp1, hlp2, hlp3, hlp4
    REAl(wp) :: ag0_cl, an0_cl, r0_cl, ci0_cl,gs0_cl           ! photosynthesis outputs of initial guess
    REAL(wp) :: fn_chl0_cl, fn_rub0_cl, fn_et0_cl, fn_pepc0_cl ! N fractions of initial guess
    REAL(wp) :: ag_cl, an_cl, r_cl, ci_cl,  gs_cl              ! photosynthesis outputs of optimisation runs
    REAL(wp) :: fn_rub_opt_cl, fn_et_opt_cl, fn_pepc_opt_cl    ! optimal N fractions
    REAL(wp) :: m_rub, m_et, m_pepc                            ! specific PS rates
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':optimise_canopy_n_fractions'



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    ! run this routine only if there are leaves at this canopy layer
    ! this routine is called in a loop over all canopy layers (DO icanopy = 1,ncanopy)
    IF(lai_cl > eps8) THEN


      !> 1.0 Set N fractions to previous values and calculate photosynthesis for avererage conditions
      ! set fractions to previous values. In current code chl + et+ rub=1
      ! all other routines chl +et +rub +oth = 1
      fn_chl0_cl   = fn_chl_cl /(1._wp-fn_oth_cl)
      fn_rub0_cl   = fn_rub_cl /(1._wp-fn_oth_cl)
      fn_et0_cl    = fn_et_cl  /(1._wp-fn_oth_cl)
      fn_pepc0_cl  = fn_pepc_cl/(1._wp-fn_oth_cl)

      fn_chl_cl    = fn_chl0_cl
      fn_et_cl     = fn_et0_cl
      fn_rub_cl    = fn_rub0_cl
      fn_pepc_cl   = fn_pepc0_cl

      ! calculate photosynthesis with previous fractions
      CALL calc_photosynthesis(  gmin                          , &  ! lctlib
                                 g0                            , &
                                 g1                            , &
                                 t_jmax_omega                  , &
                                 ps_pathway                    , &  ! lctlib
                                 canopy_cond_scheme            , &  ! Q_ASSIMI_ config (medlyn/ballberry)
                                 t_air                         , &  ! in
                                 press_srf                     , &
                                 co2_mixing_ratio              , &
                                 aerodyn_cond                  , &
                                 t_acclim                      , &
                                 ppfd_sunlit_cl                , &
                                 ppfd_shaded_cl                , &
                                 fleaf_sunlit_cl               , &
                                 beta_air                      , &
                                 beta_soa                      , &
                                 beta_soil_ps                  , &
                                 beta_sinklim_ps               , &
                                 beta_soil_gs                  , &
                                 fn_chl0_cl*(1._wp-fn_oth_cl)  , &
                                 fn_et0_cl*(1._wp-fn_oth_cl)   , &
                                 fn_rub0_cl*(1._wp-fn_oth_cl)  , &
                                 fn_pepc0_cl*(1._wp-fn_oth_cl) , &
                                 lai_cl                        , &    ! used for: IF statement whether this routine may run or not
                                 leaf_nitrogen_cl              , &
                                 t_jmax_opt                    , &
                                 ag0_cl                        , &  ! out
                                 an0_cl                        , &
                                 r0_cl                         , &
                                 gs0_cl                        , &
                                 ci0_cl                        , &
                                 m_rub                         , &
                                 m_et                          , &
                                 m_pepc                        , &
                                 hlp1  )

      !> 2.0 Attempt to increase PS by increasing the N fraction in chlorophyll
      !!
      fn_chl_cl = MIN(fn_chl0_cl * (1.0_wp + delta_n_fraction*dtime/one_day), 1.0_wp)

      ! optimal N fraction assuming Aj=Av=Apepc
      CALL calc_coordination_n_fractions(ps_pathway, &
                                         fn_chl_cl, m_rub, m_et, m_pepc, &
                                         fn_rub_opt_cl, fn_et_opt_cl, fn_pepc_opt_cl)

      ! check if changing fraction to optimal values is possible in one timestep
      CALL calc_update_n_fractions(ps_pathway, &
                                   delta_n_fraction*dtime/one_day, &
                                   fn_chl_cl, fn_rub_opt_cl, fn_et_opt_cl, fn_pepc_opt_cl, &
                                   fn_rub_cl, fn_et_cl, fn_pepc_cl)

      ! re-calculate photosynthesis
      CALL calc_photosynthesis( &
               gmin, g0, g1, t_jmax_omega, ps_pathway, &
               canopy_cond_scheme, &  ! Q_ASSIMI_ config (medlyn/ballberry)
               t_air, press_srf, co2_mixing_ratio, aerodyn_cond, t_acclim, &
               ppfd_sunlit_cl, ppfd_shaded_cl, fleaf_sunlit_cl, &
               beta_air, beta_soa,beta_soil_ps, beta_sinklim_ps, beta_soil_gs, &
               fn_chl_cl*(1._wp-fn_oth_cl), fn_et_cl*(1._wp-fn_oth_cl), fn_rub_cl*(1._wp-fn_oth_cl), &
               fn_pepc_cl*(1._wp-fn_oth_cl),&
               lai_cl                        , &    ! used for: IF statement whether this routine may run or not
               leaf_nitrogen_cl, &
               t_jmax_opt, &
               ag_cl, an_cl, r_cl, gs_cl, ci_cl, hlp1, hlp2, hlp3, &
               hlp4  )


      !> 3.0 check if increased chlorophyll has increased photosynthesis. If not, try decreasing it.
      !!
      IF (an_cl < an0_cl) THEN

        ! decrease N fraction in chlorophyll
        fn_chl_cl  = MAX(fn_chl0_cl * (1._wp - delta_n_fraction*dtime/one_day), eps8)
        fn_et_cl   = fn_et0_cl
        fn_rub_cl  = fn_rub0_cl
        fn_pepc_cl = fn_pepc0_cl

        ! optimal N fraction assuming Aj=Av=Apepc
        CALL calc_coordination_n_fractions(ps_pathway, &
                                           fn_chl_cl,m_rub,m_et,m_pepc, &
                                           fn_rub_opt_cl,fn_et_opt_cl,fn_pepc_opt_cl)
        ! check if changing fraction to optimal values is possible in one timestep
        CALL calc_update_n_fractions(ps_pathway, &
                                     delta_n_fraction*dtime/one_day, &
                                     fn_chl_cl,fn_rub_opt_cl,fn_et_opt_cl,fn_pepc_opt_cl, &
                                     fn_rub_cl,fn_et_cl,fn_pepc_cl)

        ! re-calculate photosynthesis
        CALL calc_photosynthesis( &
                 gmin, g0, g1, t_jmax_omega, ps_pathway, &
                 canopy_cond_scheme, &  ! Q_ASSIMI_ config (medlyn/ballberry)
                 t_air, press_srf, co2_mixing_ratio, aerodyn_cond, t_acclim, &
                 ppfd_sunlit_cl, ppfd_shaded_cl, fleaf_sunlit_cl,&
                 beta_air, beta_soa,beta_soil_ps, beta_sinklim_ps,beta_soil_gs, &
                 fn_chl_cl*(1._wp-fn_oth_cl), fn_et_cl*(1._wp-fn_oth_cl), fn_rub_cl*(1._wp-fn_oth_cl), &
                 fn_pepc_cl*(1._wp-fn_oth_cl),&
                 lai_cl                        , &    ! used for: IF statement whether this routine may run or not
                 leaf_nitrogen_cl, &
                 t_jmax_opt, &
                 ag_cl, an_cl, r_cl, gs_cl, ci_cl, hlp1, hlp2, hlp3, &
                 hlp4 )

        ! if photosyntheis has still not increased, do not change fn_chl_cl
        IF (an_cl < an0_cl) THEN
          fn_chl_cl  = fn_chl0_cl
          fn_et_cl   = fn_et0_cl
          fn_rub_cl  = fn_rub0_cl
          fn_pepc_cl = fn_pepc0_cl

          ! optimal N fraction assuming Aj=Av=Apepc
          CALL calc_coordination_n_fractions(ps_pathway, &
                                            fn_chl_cl,m_rub,m_et,m_pepc, &
                                            fn_rub_opt_cl,fn_et_opt_cl,fn_pepc_opt_cl)
          ! check if changing fraction to optimal values is possible in one timestep
          CALL calc_update_n_fractions(ps_pathway, &
                                      delta_n_fraction*dtime/one_day, &
                                      fn_chl_cl,fn_rub_opt_cl,fn_et_opt_cl,fn_pepc_opt_cl, &
                                      fn_rub_cl,fn_et_cl,fn_pepc_cl)

        ENDIF ! IF (an_cl < an0_cl) after re-calculating photosynthesis
      ENDIF   ! IF (an_cl < an0_cl)


      !> 4.0 Transform optimal fractions to fractions of total canopy N
      !!
      fn_chl_cl   = fn_chl_cl *(1._wp-fn_oth_cl)
      fn_rub_cl   = fn_rub_cl *(1._wp-fn_oth_cl)
      fn_et_cl    = fn_et_cl  *(1._wp-fn_oth_cl)
      fn_pepc_cl  = fn_pepc_cl*(1._wp-fn_oth_cl)

    ENDIF ! IF(lai_cl > eps8)

  END SUBROUTINE optimise_canopy_n_fractions


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to optimise_canopy_n_fractions
  !
  ! --------------------------------------------------------------------------------------------------
  !> Subroutine to calculate fractions of rubisco, electron transport and PepC assuming them to
  !! be coordinated (Aj=Av=Apepc)
  !!
  !! @todo: veg_canopy:calc_coordination_n_fractions: unsafe division by m_et/m_pepc
  !-------------------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE calc_coordination_n_fractions(ps_pathway     , &  ! lctlib
                                                     fn_chl_cl      , &
                                                     m_rub          , &
                                                     m_et           , &
                                                     m_pepc         , &
                                                     fn_rub_opt_cl  , &
                                                     fn_et_opt_cl   , &
                                                     fn_pepc_opt_cl)

    USE mo_jsb_math_constants,          ONLY: eps8
    USE mo_q_assimi_constants,          ONLY: ic4phot

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,  INTENT(IN)              :: ps_pathway          !< land-cover-type library parameter
    REAL(wp), INTENT(in)              :: fn_chl_cl           !< fraction of N in chlorophyll
    REAL(wp), INTENT(in)              :: m_rub               !< N-specific rate constant for Rubisco
    REAL(wp), INTENT(in)              :: m_et                !< N-specific rate constant for ET
    REAL(wp), INTENT(in)              :: m_pepc              !< N-speficic rate constant for PepC (C4-only)
    REAL(wp), INTENT(out)             :: fn_rub_opt_cl  , &  !< optimal fraction of N in Rubisco
                                         fn_et_opt_cl   , &  !< optimal fraction of N in ET
                                         fn_pepc_opt_cl      !< optimal fraction of N in PepC (C4-only)
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)                    :: hlp1, hlp2
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_coordination_n_fractions'


    !> 1.0
    !! optimal N fraction assuming Aj=Av=Apepc
    !! ensure m_pepc is > zero (m_pepc is INTENT(in) so hlp1 is needed)
    IF (m_pepc < eps8) THEN
      hlp1 = eps8
    ELSE
      hlp1 = m_pepc
    END IF
    ! ensure m_et is > zero
    IF (m_et < eps8) THEN
      hlp2 = eps8
    ELSE
      hlp2 = m_et
    END IF
    ! C4 / C3 plant
    IF (ps_pathway == ic4phot) THEN
       ! fraction in Rubisco
       fn_rub_opt_cl  = (1._wp - fn_chl_cl)/(1._wp + m_rub/hlp2 + m_rub/hlp1)
       ! fraction in electron transport
       fn_et_opt_cl   = fn_rub_opt_cl * m_rub/hlp2
       ! fraction in PEPc
       fn_pepc_opt_cl = fn_rub_opt_cl * m_rub/hlp1
       IF (fn_pepc_opt_cl < eps8) fn_pepc_opt_cl = eps8
    ELSE
       ! fraction in Rubisco
       fn_rub_opt_cl  = (1._wp - fn_chl_cl)/(1._wp + (m_rub)/(hlp2))
       ! fraction in electron transport
       fn_et_opt_cl   = fn_rub_opt_cl * (m_rub)/(hlp2)
       ! fraction in PEPc
       fn_pepc_opt_cl = 0._wp
    ENDIF
    ! ensure fn_rub_opt_cl & fn_et_opt_cl are > zero
    IF (fn_rub_opt_cl < eps8) fn_rub_opt_cl = eps8
    IF (fn_et_opt_cl < eps8)  fn_et_opt_cl  = eps8

  END SUBROUTINE calc_coordination_n_fractions


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to optimise_canopy_n_fractions
  !
  ! --------------------------------------------------------------------------------------------------
  !> Subroutine to update N fractions towards their optimal values in the limit of the
  !! imposed maximum rate of change
  !-------------------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE calc_update_n_fractions(ps_pathway         , &  ! lctlib
                                               delta_n_fraction   , &  ! in
                                               fn_chl_cl          , &
                                               fn_rub_opt_cl      , &
                                               fn_et_opt_cl       , &
                                               fn_pepc_opt_cl     , &
                                               fn_rub_cl          , &  ! inout
                                               fn_et_cl           , &
                                               fn_pepc_cl)

    USE mo_jsb_math_constants,          ONLY: eps8
    USE mo_q_assimi_constants,          ONLY: ic4phot

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,    INTENT(in)          :: ps_pathway             !< land-cover-type library parameter
    REAL(wp),   INTENT(in)          :: delta_n_fraction , &   !< maximum fractional change in this time-step
                                       fn_chl_cl        , &   !< fraction of N in chlorophyll
                                       fn_rub_opt_cl    , &   !< optimal fraction of N in Rubisco
                                       fn_et_opt_cl     , &   !< optimal fraction of N in ET
                                       fn_pepc_opt_cl         !< optimal fraction of N in PepC (C4-only)
    REAL(wp),   INTENT(inout)       :: fn_rub_cl        , &   !< updated fraction of N in Rubisco
                                       fn_et_cl         , &   !< updated fraction of N in ET
                                       fn_pepc_cl             !< updated fraction of N in PepC (C4-only)
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)                        :: hlp1, hlp2, hlp3
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_update_n_fractions'



    ! 0.9 init local var
    hlp1 = ABS(fn_rub_opt_cl-fn_rub_cl) / fn_rub_cl
    hlp2 = ABS(fn_et_opt_cl -fn_et_cl)  / fn_et_cl


    ! only C4 plant have pepc, avoid dividing by zero for C3
    IF(ps_pathway == ic4phot) THEN
       hlp3 = ABS(fn_pepc_opt_cl-fn_pepc_cl)/ fn_pepc_cl
    ELSE
       hlp3 = 0.0_wp
    ENDIF


    ! check if changing fraction to optimal values is possible in one timestep
    IF (hlp1 <= delta_n_fraction .AND. &
        hlp2 <= delta_n_fraction .AND. &
        hlp3 <= delta_n_fraction) THEN
       fn_rub_cl  = fn_rub_opt_cl
       fn_et_cl   = fn_et_opt_cl
       fn_pepc_cl = fn_pepc_opt_cl
    ! if not, change by maximum amount into right direction
    ELSEIF (fn_rub_cl < fn_rub_opt_cl) THEN
       fn_rub_cl = MIN(fn_rub_cl * (1.0_wp+delta_n_fraction), fn_rub_opt_cl)
       IF (fn_pepc_cl < fn_pepc_opt_cl) THEN
          fn_pepc_cl  = MAX(fn_pepc_cl * (1.0_wp+delta_n_fraction), fn_pepc_opt_cl)
          fn_et_cl    = MAX(1.0_wp -fn_rub_cl - fn_pepc_cl - fn_chl_cl, eps8)
       ELSE
          fn_pepc_cl  = MIN(fn_pepc_cl * (1.0_wp-delta_n_fraction), fn_pepc_opt_cl)
          fn_et_cl    = MAX(1.0_wp -fn_rub_cl - fn_pepc_cl - fn_chl_cl, eps8)
       ENDIF
    ! else ...
    ELSE
       fn_rub_cl = MIN(fn_rub_cl * (1.0_wp-delta_n_fraction), fn_rub_opt_cl)
       IF (fn_pepc_cl < fn_pepc_opt_cl) THEN
          fn_pepc_cl  = MAX(fn_pepc_cl * (1.0_wp+delta_n_fraction), fn_pepc_opt_cl)
          fn_et_cl    = MAX(1.0_wp -fn_rub_cl - fn_pepc_cl - fn_chl_cl, eps8)
       ELSE
          fn_pepc_cl  = MIN(fn_pepc_cl * (1.0_wp-delta_n_fraction), fn_pepc_opt_cl)
          fn_et_cl    = MAX(1.0_wp -fn_rub_cl - fn_pepc_cl - fn_chl_cl, eps8)
       ENDIF
    ENDIF

  END SUBROUTINE calc_update_n_fractions


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to calc_growth (currently in mo_q_veg_canopy to get rid of circular dependencies)
  !
  !------------------------------------------------------------------------------------------------------
  !> Subroutine to calculate change in net photosynthesis and transpiration with small change in LAI
  !! Needed for calculation of optimal leaf and root allocation fraction. See mo_q_veg_growth
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_marginal_canopy_flux_increment( nc                      , &
                                                  ncanopy                 , &
                                                  dtime                   , &
                                                  dz                      , &
                                                  lbounds                 , &
                                                  ubounds                 , &
                                                  lctlib_ps_pathway       , &
                                                  lctlib_sla              , &
                                                  lctlib_k0_fn_struc      , &
                                                  lctlib_fn_oth_min       , &
                                                  lctlib_np_leaf          , &
                                                  lctlib_cn_leaf          , &
                                                  lctlib_gmin             , &
                                                  lctlib_g0               , &
                                                  lctlib_g1               , &
                                                  lctlib_t_jmax_omega     , &
                                                  lctlib_sigma_vis        , &
                                                  flag_optimal_Nfraction  , &
                                                  canopy_cond_scheme      , & ! Q_ASSIMI_ config
                                                  t_air                   , &
                                                  press_srf               , &
                                                  q_air                   , &
                                                  co2_mixing_ratio        , &
                                                  aerodyn_cond            , &
                                                  net_assimilation        , &
                                                  canopy_cond             , &
                                                  beta_air                , &
                                                  beta_soa                , &
                                                  beta_soil_ps            , &
                                                  beta_soil_gs            , &
                                                  veg_pool_leaf_carbon    , &
                                                  veg_pool_leaf_nitrogen  , &
                                                  lai                     , &
                                                  beta_sinklim_ps         , &
                                                  t_air_tacclim_mavg      , &
                                                  t_jmax_opt_mavg         , &
                                                  fn_chl_cl               , &
                                                  fn_et_cl                , &
                                                  fn_rub_cl               , &
                                                  fn_pepc_cl              , &
                                                  fn_oth_cl               , &
                                                  fleaf_sunlit_cl         , &
                                                  ppfd_sunlit_cl          , &
                                                  ppfd_shaded_cl          , &
                                                  swpar_srf_down          , &
                                                  fract_par_diffuse       , &
                                                  unit_npp                , &
                                                  unit_transpiration )

    USE mo_jsb_math_constants,     ONLY: eps8
    USE mo_jsb_physical_constants, ONLY: r_gas_dryair
    USE mo_q_rad_parameters,       ONLY: rad2ppfd, kbl0_vis
    USE mo_q_assimi_constants,     ONLY: ic3phot, ic4phot
    USE mo_atmland_util,           ONLY: calc_spec_humidity_sat

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                INTENT(in) :: nc                            !< dimensions
    INTEGER,                                INTENT(in) :: ncanopy                       !< number of canopy layers
    REAL(wp),                               INTENT(in) :: dtime                         !< timestep length
    REAL(wp),DIMENSION(ncanopy),            INTENT(in) :: dz                      , &   !< canopy layer thickness
                                                          lbounds                 , &   !< lower bound of canopy layer (lbounds < ubounds; lbounds is at top of canopy layer)
                                                          ubounds                       !< upper bound of canopy layer (lbounds < ubounds; ubounds is at bottom of canopy layer)
    INTEGER,                                INTENT(in) :: lctlib_ps_pathway             !< lctlib parameter
    REAL(wp),                               INTENT(in) :: lctlib_sla                    !< lctlib parameter
    REAL(wp),                               INTENT(in) :: lctlib_k0_fn_struc            !< lctlib parameter
    REAL(wp),                               INTENT(in) :: lctlib_fn_oth_min             !< lctlib parameter
    REAL(wp),                               INTENT(in) :: lctlib_np_leaf                !< lctlib parameter
    REAL(wp),                               INTENT(in) :: lctlib_cn_leaf                !< lctlib parameter
    REAL(wp),                               INTENT(in) :: lctlib_gmin                   !< lctlib parameter
    REAL(wp),                               INTENT(in) :: lctlib_g0                     !< lctlib parameter
    REAL(wp),                               INTENT(in) :: lctlib_g1                     !< lctlib parameter
    REAL(wp),                               INTENT(in) :: lctlib_t_jmax_omega           !< lctlib parameter
    REAL(wp),                               INTENT(in) :: lctlib_sigma_vis              !< lctlib parameter
    LOGICAL,                                INTENT(in) :: flag_optimal_Nfraction        !< on/off optimise leaf internal N allocation
    CHARACTER(len=*),                       INTENT(IN) :: canopy_cond_scheme            !< canopy_conductance_scheme: medlyn / ballberry
    REAL(wp),DIMENSION(nc),                 INTENT(in) :: t_air                   , &   !< air temperature (K)
                                                          press_srf               , &   !< air pressure (Pa)
                                                          q_air                   , &   !< air specific humidity (--)
                                                          co2_mixing_ratio        , &   !< co2 mixing ratio (ppm)
                                                          aerodyn_cond            , &   !< ga -- aerodynamic conductance (day-time running mean, m/s)
                                                          net_assimilation        , &   !< net assimilation rate (micro-mol CO2 / m2 / s)
                                                          canopy_cond             , &   !< canopy conductance (m/s)
                                                          beta_air                , &   !< air humidity limiting factor for photosynthesis
                                                          beta_soa                , &   !< spring air temp lim factor for Evergr. needlel. photosynthesis
                                                          beta_soil_ps            , &   !< soil moisture limiting factor for photosynthesis
                                                          beta_soil_gs                  !< soil moisture limiting factor for canopy conductance
    REAL(wp),DIMENSION(nc),                 INTENT(in) :: veg_pool_leaf_carbon    , &   !< ..
                                                          veg_pool_leaf_nitrogen  , &   !< ..
                                                          lai                     , &   !< lai
                                                          beta_sinklim_ps         , &   !< ..
                                                          t_air_tacclim_mavg      , &   !< ..
                                                          t_jmax_opt_mavg               !< ..
    REAL(wp),DIMENSION(nc,ncanopy),         INTENT(in) :: fn_chl_cl               , &   !< ..
                                                          fn_et_cl                , &   !< ..
                                                          fn_rub_cl               , &   !< ..
                                                          fn_pepc_cl              , &   !< ..
                                                          fn_oth_cl               , &   !< ..
                                                          fleaf_sunlit_cl         , &   !< ..
                                                          ppfd_sunlit_cl          , &   !< ..
                                                          ppfd_shaded_cl                !< ..
    REAl(wp),DIMENSION(nc),                 INTENT(in) :: swpar_srf_down          , &   !< ..
                                                          fract_par_diffuse             !< ..
    REAL(wp),DIMENSION(nc),                 INTENT(out):: unit_npp                , &   !< unit change of NPP per ??
                                                          unit_transpiration            !< unit change of NPP per ??
    ! ---------------------------
    ! 0.2 Local
    REAL(wp),DIMENSION(nc)                             :: leaf_nitrogen_new
    REAL(wp),DIMENSION(nc)                             :: lai_new
    REAL(wp),DIMENSION(nc)                             :: net_assimilation_new
    REAL(wp),DIMENSION(nc)                             :: canopy_cond_new
    REAL(wp),DIMENSION(nc)                             :: beta_sinklim_local
    REAL(wp),DIMENSION(nc,ncanopy)                     :: lai_new_cl                , &
                                                          cumm_lai_new_cl           , &
                                                          loc_ppfd_sunlit_cl        , &
                                                          loc_ppfd_shaded_cl        , &
                                                          loc_fleaf_sunlit_cl       , &
                                                          fn_oth_new_cl             , &
                                                          leaf_nitrogen_new_cl      , &
                                                          fn_chl_new_cl             , &
                                                          fn_et_new_cl              , &
                                                          fn_rub_new_cl             , &
                                                          fn_pepc_new_cl
    REAL(wp),DIMENSION(nc)                             :: ag                        , &
                                                          an                        , &
                                                          resp                      , &
                                                          gs                        , &
                                                          q_sat                     , &
                                                          hlp0, hlp1, hlp2          , &   ! calc_photosynthesis dummy output
                                                          hlp3, hlp4                      ! calc_photosynthesis dummy output

    INTEGER                                            :: icanopy !< looping
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_marginal_canopy_flux_increment'



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 0.9 init some local variables
    !!
    fn_oth_new_cl(:,:)        = fn_oth_cl(:,:)
    fn_chl_new_cl(:,:)        = fn_chl_cl(:,:)
    fn_et_new_cl(:,:)         = fn_et_cl(:,:)
    fn_rub_new_cl(:,:)        = fn_rub_cl(:,:)
    fn_pepc_new_cl(:,:)       = fn_pepc_cl(:,:)
    beta_sinklim_local        = 1._wp

    !> 1.0 copy current plant state into dummy variables
    !!
    loc_ppfd_sunlit_cl(:,:)   = ppfd_sunlit_cl(:,:)
    loc_ppfd_shaded_cl(:,:)   = ppfd_shaded_cl(:,:)
    loc_fleaf_sunlit_cl(:,:)  = fleaf_sunlit_cl(:,:)

 !> @todo: mo_q_veg_canopy calc_npp_increment: I think this is wrong currently, as increaseing leaf area require an
 !! increase in roots and sapwood mass and thus respiration... Also, should be based on current mass
 !! ratios, not a fixed constant. ie
 !! f_increment = f_alloc_leaf/f_alloc_root+sapwood
 !! maint_resp = maint_resp * (1 + f_alloc_root+sapwood*1/totalmassofpools)
 !! probably only relevant for BNF and uptake calculation?

    !> 2.0 Estimate effect of a small increase in leaf C on LAI and leaf N
    !!
    ! increase LAI by 1 mol of C
    lai_new(:) = lai(:) + lctlib_sla
    ! equivalent leaf N assuming any new leaf growth will have target CN
    WHERE(veg_pool_leaf_carbon(:) > eps8)
       leaf_nitrogen_new(:) = veg_pool_leaf_nitrogen(:) + veg_pool_leaf_nitrogen(:) / veg_pool_leaf_carbon(:)
    ELSEWHERE
       leaf_nitrogen_new(:) = 1._wp / lctlib_cn_leaf
    END WHERE

    ! recompute canopy layer LAI and foliar N with new numbers
    CALL calc_canopy_layers( &
                            nc                      = nc                          , & ! in
                            ncanopy                 = ncanopy                     , &
                            dtime                   = dtime                       , &
                            dz                      = dz(:)                       , &
                            lbounds                 = lbounds(:)                  , &
                            ubounds                 = ubounds(:)                  , &
                            lctlib_ps_pathway       = lctlib_ps_pathway           , &
                            lctlib_k0_fn_struc      = lctlib_k0_fn_struc          , &
                            lctlib_fn_oth_min       = lctlib_fn_oth_min           , &
                            lctlib_sla              = lctlib_sla                  , &
                            lctlib_np_leaf          = lctlib_np_leaf              , &
                            lctlib_gmin             = lctlib_gmin                 , &
                            lctlib_g0               = lctlib_g0                   , &
                            lctlib_g1               = lctlib_g1                   , &
                            lctlib_t_jmax_omega     = lctlib_t_jmax_omega         , &
                            flag_optimal_Nfraction  = flag_optimal_Nfraction      , &
                            canopy_cond_scheme      = canopy_cond_scheme          , & !   Q_ASSIMI_ config
                            leaf_nitrogen           = leaf_nitrogen_new(:)        , &
                            lai                     = lai_new(:)                  , & ! in
                            ppfd_sunlit_cl          = loc_ppfd_sunlit_cl(:,:)     , & ! inout
                            ppfd_shaded_cl          = loc_ppfd_shaded_cl(:,:)     , &
                            fleaf_sunlit_cl         = loc_fleaf_sunlit_cl(:,:)    , &
                            fn_rub_cl               = fn_rub_new_cl(:,:)          , &
                            fn_et_cl                = fn_et_new_cl(:,:)           , &
                            fn_pepc_cl              = fn_pepc_new_cl(:,:)         , &
                            fn_chl_cl               = fn_chl_new_cl(:,:)          , &
                            fn_oth_cl               = fn_oth_new_cl(:,:)          , & ! inout
                            lai_cl                  = lai_new_cl(:,:)             , & ! out
                            cumm_lai_cl             = cumm_lai_new_cl(:,:)        , & ! out
                            leaf_nitrogen_cl        = leaf_nitrogen_new_cl(:,:)   , & ! out
                            inquire_n_fractions     = .TRUE.                      )   ! optional in

    ! In case of no leaves set top layer PAR to SW value
    WHERE (lai(:) < eps8)
       loc_ppfd_sunlit_cl(:,1)  = (1._wp - lctlib_sigma_vis) * swpar_srf_down(:) * rad2ppfd
       loc_ppfd_shaded_cl(:,1)  = (1._wp - lctlib_sigma_vis) * fract_par_diffuse(:) * swpar_srf_down(:) * rad2ppfd
       loc_fleaf_sunlit_cl(:,1) = exp(-kbl0_vis * cumm_lai_new_cl(:,1))
    END WHERE

    ! calculate photosynthesis by integration over canopy layers containing leaves
    net_assimilation_new(:)  = 0.0_wp
    canopy_cond_new(:)       = 0.0_wp
    DO icanopy = 1,ncanopy
      CALL calc_photosynthesis(     lctlib_gmin                     , & ! lctlib
                                    lctlib_g0                       , &
                                    lctlib_g1                       , &
                                    lctlib_t_jmax_omega             , &
                                    lctlib_ps_pathway               , &
                                    canopy_cond_scheme              , & ! Q_ASSIMI_ config (medlyn/ballberry)
                                    t_air(:)                        , & ! in
                                    press_srf(:)                    , &
                                    co2_mixing_ratio(:)             , &
                                    aerodyn_cond(:)                 , &
                                    t_air_tacclim_mavg(:)           , &
                                    loc_ppfd_sunlit_cl(:,icanopy)   , &
                                    loc_ppfd_shaded_cl(:,icanopy)   , &
                                    loc_fleaf_sunlit_cl(:,icanopy)  , &
                                    beta_air(:)                     , &
                                    beta_soa(:)                     , &
                                    beta_soil_ps(:)                 , &
                                    beta_sinklim_local(:)           , &
                                    beta_soil_gs(:)                 , &
                                    fn_chl_new_cl(:,icanopy)        , &
                                    fn_et_new_cl(:,icanopy)         , &
                                    fn_rub_new_cl(:,icanopy)        , &
                                    fn_pepc_new_cl(:,icanopy)       , &
                                    lai_new_cl(:,icanopy)           , &   ! used for: IF statement whether this routine may run or not
                                    leaf_nitrogen_new_cl(:,icanopy) , &
                                    t_jmax_opt_mavg(:)              , &
                                    ag(:)                           , & ! out (values per canopy layer)
                                    an(:)                           , &
                                    resp(:)                         , &
                                    gs(:)                           , &
                                    hlp0(:)                         , &
                                    hlp1(:)                         , &
                                    hlp2(:)                         , &
                                    hlp3(:)                         , &
                                    hlp4(:)                          )

      net_assimilation_new(:) = net_assimilation_new(:) + an(:) * lai_new_cl(:,icanopy)
      canopy_cond_new(:)      = canopy_cond_new(:)      + gs(:) * lai_new_cl(:,icanopy)
    ENDDO


    !> 3.0 Calculate change in canopy fluxes from altered LAI
    !!
    !> 3.1 net carbon gain from increased LAI
    !!
    unit_npp(:) = net_assimilation_new(:) - net_assimilation(:)

    !> 3.2 Calculate increment in net water loss from increased LAI
    !!
    q_sat(:) = calc_spec_humidity_sat(t_air(:), press_srf(:))
    WHERE(canopy_cond(:) > eps8 .AND. canopy_cond_new(:) > eps8)
       unit_transpiration(:) = press_srf(:) / ( r_gas_dryair * t_air(:) ) * &
                                (q_sat(:) - q_air(:)) / (1._wp / canopy_cond_new(:) + 1._wp / aerodyn_cond(:) ) - &
                                press_srf(:) / ( r_gas_dryair * t_air(:) ) * &
                                (q_sat(:) - q_air(:)) / (1._wp / canopy_cond(:) + 1._wp / aerodyn_cond(:) )
    ELSEWHERE (canopy_cond_new(:) > eps8) ! this might cover only the points in the vector that are not true in the above where statement
       unit_transpiration(:) = press_srf(:) / ( r_gas_dryair * t_air(:) ) * &
                                (q_sat(:) - q_air(:)) / (1._wp / canopy_cond_new(:) + 1._wp / aerodyn_cond(:) )
    ELSEWHERE
      unit_transpiration(:) = 0.0_wp
    END WHERE

   END SUBROUTINE calc_marginal_canopy_flux_increment


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to calc_growth (currently in mo_q_veg_canopy to get rid of circular dependencies)
  !
  !------------------------------------------------------------------------------------------------------
  !> Subroutine to determine direction in which total leaf N should change in order to increase NPP
  !!
  !! Parameter leaf_cn_direction
  !----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_dir_optimal_cn_leaf(  nc                            , & ! in
                                        ncanopy                       , &
                                        dtime                         , &
                                        dz                            , &
                                        lbounds                       , &
                                        ubounds                       , &
                                        lctlib_ps_pathway             , &
                                        lctlib_k0_fn_struc            , &
                                        lctlib_fn_oth_min             , &
                                        lctlib_sla                    , &
                                        lctlib_np_leaf                , &
                                        lctlib_gmin                   , &
                                        lctlib_g0                     , &
                                        lctlib_g1                     , &
                                        lctlib_t_jmax_omega           , &
                                        lctlib_tau_leaf               , &
                                        flag_optimal_Nfraction        , &
                                        canopy_cond_scheme            , &  ! Q_ASSIMI_ config (medlyn/ballberry)
                                        growing_season                , &
                                        veg_pool_leaf_carbon          , &
                                        veg_pool_leaf_nitrogen        , &
                                        veg_pool_root_carbon          , &
                                        lai                           , &
                                        t_jmax_opt_mavg               , &
                                        target_cn_leaf                , &
                                        target_cn_fine_root           , &
                                        t_air_tacclim_mavg            , &
                                        t_air_tcnl_mavg               , &
                                        ga_tcnl_mavg                  , &
                                        press_srf_tcnl_mavg           , &
                                        co2_mixing_ratio_tcnl_mavg    , &
                                        beta_air_tcnl_mavg            , &
                                        beta_soa_tphen_mavg           , &
                                        beta_soil_ps_tcnl_mavg        , &
                                        beta_soil_gs_tcnl_mavg        , &
                                        fn_chl_cl                     , &
                                        fn_et_cl                      , &
                                        fn_rub_cl                     , &
                                        fn_pepc_cl                    , &
                                        fn_oth_cl                     , &
                                        fleaf_sunlit_tcnl_mavg_cl     , &
                                        ppfd_sunlit_tcnl_mavg_cl      , &
                                        ppfd_shaded_tcnl_mavg_cl      , &
                                        uptake_n                      , &
                                        growth_cn                     , &
                                        npp                           , &
                                        fmaint_rate                   , &
                                        labile_carbon                 , &
                                        labile_nitrogen               , & ! in
                                        leaf_cn_direction             )   ! inout

    USE mo_jsb_math_constants,            ONLY: one_day, one_year, eps4, eps8
    USE mo_veg_constants,                 ONLY: delta_n_leaf, resorp_fract_leaf
    USE mo_lnd_time_averages,             ONLY: mavg_period_tgrowth
    USE mo_q_pheno_constants,             ONLY: ievergreen, iraingreen, isummergreen, iperennial

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                INTENT(in)    :: nc                             !< dimensions
    INTEGER,                                INTENT(in)    :: ncanopy                        !< number of canopy layers
    REAL(wp),                               INTENT(in)    :: dtime                          !< timestep length
    REAL(wp),DIMENSION(ncanopy),            INTENT(in)    :: dz                         , & !< canopy layer thickness
                                                             lbounds                    , & !< lower bound of canopy layer (lbounds < ubounds; lbounds is at top of canopy layer)
                                                             ubounds                        !< upper bound of canopy layer (lbounds < ubounds; ubounds is at bottom of canopy layer)
    INTEGER,                                INTENT(in)    :: lctlib_ps_pathway              !< lctlib parameter
    REAL(wp),                               INTENT(in)    :: lctlib_k0_fn_struc             !< lctlib parameter
    REAL(wp),                               INTENT(in)    :: lctlib_fn_oth_min              !< lctlib parameter
    REAL(wp),                               INTENT(in)    :: lctlib_sla                     !< lctlib parameter
    REAL(wp),                               INTENT(in)    :: lctlib_np_leaf                 !< lctlib parameter
    REAL(wp),                               INTENT(in)    :: lctlib_gmin                    !< lctlib parameter
    REAL(wp),                               INTENT(in)    :: lctlib_g0                      !< lctlib parameter
    REAL(wp),                               INTENT(in)    :: lctlib_g1                      !< lctlib parameter
    REAL(wp),                               INTENT(in)    :: lctlib_t_jmax_omega            !< lctlib parameter
    REAL(wp),                               INTENT(in)    :: lctlib_tau_leaf                !< lctlib parameter
    LOGICAL,                                INTENT(in)    :: flag_optimal_Nfraction         !< on/off optimise leaf internal N allocation
    CHARACTER(len=*),                       INTENT(IN)    :: canopy_cond_scheme             !< canopy_conductance_scheme: medlyn / ballberry
    REAL(wp),DIMENSION(nc),                 INTENT(in)    :: growing_season                 !< growing season
    REAL(wp),DIMENSION(nc),                 INTENT(in)    :: veg_pool_leaf_carbon       , &
                                                             veg_pool_leaf_nitrogen     , &
                                                             veg_pool_root_carbon
    REAL(wp),DIMENSION(nc),                 INTENT(in)    :: lai                        , &
                                                             t_jmax_opt_mavg            , &
                                                             target_cn_leaf             , &
                                                             target_cn_fine_root        , &
                                                             t_air_tacclim_mavg         , &
                                                             t_air_tcnl_mavg            , &
                                                             ga_tcnl_mavg               , &
                                                             press_srf_tcnl_mavg        , &
                                                             co2_mixing_ratio_tcnl_mavg , &
                                                             beta_air_tcnl_mavg         , &
                                                             beta_soa_tphen_mavg        , &
                                                             beta_soil_ps_tcnl_mavg     , &
                                                             beta_soil_gs_tcnl_mavg
    REAL(wp),DIMENSION(nc,ncanopy),         INTENT(in)    :: fn_chl_cl                  , &
                                                             fn_et_cl                   , &
                                                             fn_rub_cl                  , &
                                                             fn_pepc_cl                 , &
                                                             fn_oth_cl                  , &
                                                             fleaf_sunlit_tcnl_mavg_cl  , &
                                                             ppfd_sunlit_tcnl_mavg_cl   , &
                                                             ppfd_shaded_tcnl_mavg_cl
    REAL(wp),DIMENSION(nc),                 INTENT(in)    :: uptake_n                   , &
                                                             growth_cn                  , &
                                                             npp                        , &
                                                             fmaint_rate                , &
                                                             labile_carbon              , &
                                                             labile_nitrogen
    REAL(wp),DIMENSION(nc),                 INTENT(inout) :: leaf_cn_direction
    ! ---------------------------
    ! 0.2 Local
    REAL(wp),DIMENSION(nc)                                :: leaf_nitrogen_new, root_nitrogen_new
    REAL(wp),DIMENSION(nc,ncanopy)                        :: leaf_nitrogen_new_cl, &
                                                             fn_oth_new_cl, &
                                                             fn_chl_new_cl, &
                                                             fn_et_new_cl, &
                                                             fn_rub_new_cl, &
                                                             fn_pepc_new_cl, &
                                                             lai_new_cl, &
                                                             cumm_lai_new_cl
    REAL(wp),DIMENSION(nc,ncanopy)                        :: loc_ppfd_sunlit_cl   , &
                                                             loc_ppfd_shaded_cl   , &
                                                             loc_fleaf_sunlit_cl
    REAL(wp),DIMENSION(3)                                 :: npp_direction
    REAL(wp)                                              :: anl, agl, rl
    REAL(wp)                                              :: hlp00, hlp0, hlp1        ! calc_photosynthesis dummy output
    REAL(wp)                                              :: hlp2, hlp3, hlp4         ! calc_photosynthesis dummy output
    REAL(wp)                                              :: direction
    REAL(wp), DIMENSION(nc)                               :: beta_sinklim_ps_local
    INTEGER                                               :: icanopy, dd, ic      ! looping
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_dir_optimal_cn_leaf'



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 0.9 init local var
    !!
    loc_ppfd_sunlit_cl(:,:)   = ppfd_sunlit_tcnl_mavg_cl(:,:)
    loc_ppfd_shaded_cl(:,:)   = ppfd_shaded_tcnl_mavg_cl(:,:)
    loc_fleaf_sunlit_cl(:,:)  = fleaf_sunlit_tcnl_mavg_cl(:,:)
    beta_sinklim_ps_local(:)  = 1._wp

    fn_chl_new_cl(:,:)   = fn_chl_cl(:,:)
    fn_et_new_cl(:,:)    = fn_et_cl(:,:)
    fn_rub_new_cl(:,:)   = fn_rub_cl(:,:)


    !> 1.0 ...
    !!
    DO ic=1,nc
       ! Only change leaf N during the growing season
       IF(growing_season(ic) > test_false_true)THEN
          ! Only change leaf N if leaves present
          IF(veg_pool_leaf_nitrogen(ic) > eps4)THEN
             !> Check the N status of the plant. If there is not enough N for growth, decrease leaf N
             IF ((uptake_n(ic) *  mavg_period_tgrowth * one_day * 1e-6_wp + &
                   resorp_fract_leaf * veg_pool_leaf_nitrogen(ic) * 1._wp / (lctlib_tau_leaf * one_year * one_day) * &
                   mavg_period_tgrowth * one_day * 1e-6_wp + labile_nitrogen(ic)) * growth_cn(ic) < &
                 (npp(ic)*  mavg_period_tgrowth * 1e-6_wp * one_day + labile_carbon(ic))) THEN
                leaf_cn_direction(ic) = -1._wp

             ELSE
                !> Check for C optimality i.e. if increasing leaf N increases NPP
                !!
                npp_direction(:) = 0.0_wp
                DO dd = 1,3
                   direction                 = REAL(dd,KIND=wp) - 2.0_wp
                   leaf_nitrogen_new(ic) = veg_pool_leaf_carbon(ic) / target_cn_leaf(ic) * &
                        (1.0_wp + direction * delta_n_leaf * dtime/one_day)
                   root_nitrogen_new(ic) = veg_pool_root_carbon(ic) / target_cn_fine_root(ic) * &
                        (1.0_wp + direction * delta_n_leaf * dtime/one_day)

                ! update N distribution in the canopy
                CALL calc_canopy_layers( &
                          nc                      = nc                                  , & ! in
                          ncanopy                 = ncanopy                             , &
                          dtime                   = dtime                               , &
                          dz                      = dz(:)                               , &
                          lbounds                 = lbounds(:)                          , &
                          ubounds                 = ubounds(:)                          , &
                          lctlib_ps_pathway       = lctlib_ps_pathway                   , &
                          lctlib_k0_fn_struc      = lctlib_k0_fn_struc                  , &
                          lctlib_fn_oth_min       = lctlib_fn_oth_min                   , &
                          lctlib_sla              = lctlib_sla                          , &
                          lctlib_np_leaf          = lctlib_np_leaf                      , &
                          lctlib_gmin             = lctlib_gmin                         , &
                          lctlib_g0               = lctlib_g0                           , &
                          lctlib_g1               = lctlib_g1                           , &
                          lctlib_t_jmax_omega     = lctlib_t_jmax_omega                 , &
                          flag_optimal_Nfraction  = flag_optimal_Nfraction              , &
                          canopy_cond_scheme      = canopy_cond_scheme                  , & !   Q_ASSIMI_ config
                          leaf_nitrogen           = leaf_nitrogen_new(ic)               , &
                          lai                     = lai(ic)                             , & ! in
                          ppfd_sunlit_cl          = loc_ppfd_sunlit_cl(ic,:)            , & ! inout
                          ppfd_shaded_cl          = loc_ppfd_shaded_cl(ic,:)            , &
                          fleaf_sunlit_cl         = loc_fleaf_sunlit_cl(ic,:)           , &
                          fn_rub_cl               = fn_rub_new_cl(ic,:)                 , &
                          fn_et_cl                = fn_et_new_cl(ic,:)                  , &
                          fn_pepc_cl              = fn_pepc_new_cl(ic,:)                , &
                          fn_chl_cl               = fn_chl_new_cl(ic,:)                 , &
                          fn_oth_cl               = fn_oth_new_cl(ic,:)                 , & ! inout
                          lai_cl                  = lai_new_cl(ic,:)                    , & ! out
                          cumm_lai_cl             = cumm_lai_new_cl(ic,:)               , &
                          leaf_nitrogen_cl        = leaf_nitrogen_new_cl(ic,:)          , & ! out
                          inquire_n_fractions     = .TRUE.                                ) ! optional in

                ! adjust N fractions for any changes in canopy N distribtuion
                ! for optimal fractions this is only a aproximation
                IF (lai(ic) > eps8) THEN
                   ! avoid division by zero for plants that have no leaf nitrogen
                   WHERE (leaf_nitrogen_new_cl(ic,:) > eps8 .AND. fn_oth_cl(ic,:) < 1._wp)
                      fn_rub_new_cl(ic,:)  = fn_rub_cl(ic,:)  * &
                                                   (1._wp - fn_oth_new_cl(ic,:)) / (1._wp - fn_oth_cl(ic,:))
                      fn_et_new_cl(ic,:)   = fn_et_cl(ic,:)   * &
                                                   (1._wp - fn_oth_new_cl(ic,:)) / (1._wp - fn_oth_cl(ic,:))
                      fn_pepc_new_cl(ic,:) = fn_pepc_cl(ic,:) * &
                                                   (1._wp - fn_oth_new_cl(ic,:)) / (1._wp - fn_oth_cl(ic,:))
                      fn_chl_new_cl(ic,:)  = fn_chl_cl(ic,:)  * &
                                                   (1._wp - fn_oth_new_cl(ic,:)) / (1._wp - fn_oth_cl(ic,:))
                   END WHERE
                END IF

                ! calculate photosynthesis for each layer
                DO icanopy =1, ncanopy
                   CALL calc_photosynthesis( &
                       lctlib_gmin                               , & ! lctlib
                       lctlib_g0                                 , &
                       lctlib_g1                                 , &
                       lctlib_t_jmax_omega                       , &
                       lctlib_ps_pathway                         , &
                       canopy_cond_scheme                        , &  ! Q_ASSIMI_ config (medlyn/ballberry)
                       t_air_tcnl_mavg(ic)                       , &  ! in
                       press_srf_tcnl_mavg(ic)                   , &
                       co2_mixing_ratio_tcnl_mavg(ic)            , &
                       ga_tcnl_mavg(ic)                          , &
                       t_air_tacclim_mavg(ic)                    , &
                       ppfd_sunlit_tcnl_mavg_cl(ic,icanopy)      , &
                       ppfd_shaded_tcnl_mavg_cl(ic,icanopy)      , &
                       fleaf_sunlit_tcnl_mavg_cl(ic,icanopy)     , &
                       beta_air_tcnl_mavg(ic)                    , &
                       beta_soa_tphen_mavg(ic)                   , &
                       beta_soil_ps_tcnl_mavg(ic)                , &
                       beta_sinklim_ps_local(ic)                 , &
                       beta_soil_gs_tcnl_mavg(ic)                , &
                       fn_chl_new_cl(ic,icanopy)                 , &
                       fn_et_new_cl(ic,icanopy)                  , &
                       fn_rub_new_cl(ic,icanopy)                 , &
                       fn_pepc_new_cl(ic,icanopy)                , &
                       lai_new_cl(ic,icanopy)                    , &   ! used for: IF statement whether this routine may run or not
                       leaf_nitrogen_new_cl(ic,icanopy)          , &
                       t_jmax_opt_mavg(ic)                       , &
                       agl                                       , &
                       anl                                       , &
                       rl                                        , &
                       hlp00, hlp0, hlp1, hlp2, hlp3             , &
                       hlp4    )

                npp_direction(dd) = npp_direction(dd) + anl * lai_new_cl(ic,icanopy)
                ENDDO  ! icanopy =1, ncanopy

                ! adjust for change in fine root respiration
                npp_direction(dd) = npp_direction(dd) - root_nitrogen_new(ic) * fmaint_rate(ic)
                ENDDO  ! dd = 1,3

                ! Direction of change in leaf N. Actual leaf N is modified in interfacte routine
                !! @todo:veg_canopy:calc_dir_optimal_cn_leaf:potentially dangerous use of MAXLOC - ascertain functioning in 4D
                leaf_cn_direction(ic) = (REAL(MAXLOC(npp_direction, 1),KIND=wp) - 2.0_wp)
             ENDIF  ! if enough N available

          ELSE    ! IF(veg_pool_leaf_nitrogen(ic) > eps4)
             leaf_cn_direction(ic) = 0.0_wp
          END IF  ! IF(veg_pool_leaf_nitrogen(ic) > eps4)
       END IF  ! IF(growing_season(ic) > test_false_true)
    END DO  ! ic=1,nc

  END SUBROUTINE calc_dir_optimal_cn_leaf

#endif
END MODULE mo_q_veg_canopy
