! Initialization routine for multi-layer soil model
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

MODULE sfc_terra_init

!------------------------------------------------------------------------------

  USE mo_kind,               ONLY: wp
  USE mo_physical_constants, ONLY: t0_melt => tmelt,& ! absolute zero for temperature
    &                              rho_w => rhoh2o, & ! density of liquid water (kg/m^3)
    &                              lh_f  => alf   , & ! latent heat of fusion
    &                              g     => grav      ! acceleration due to gravity

  USE mo_lnd_nwp_config,     ONLY: lmulti_snow, lana_rho_snow, &
    &                              lmelt, lmelt_var,           &
    &                              max_toplaydepth, l2lay_rho_snow

  USE sfc_terra_data,        ONLY:                                &
      eps_soil        , & ! Multi-purpose epsilon in soil model (former zepsi)
      eps_temp        , & ! Small value to check if temperatures have exceeded a fixed threshold
      crhosmax_tmin, crhosmax_ml, crhosmin_ml, zdzhs, csnow_tmin, &
      cporv, cadp, csandf, cclayf, T_ref_ice, T_star_ice, b_clay, &
      b_silt, b_sand, b_org

!==============================================================================

  IMPLICIT NONE

  PRIVATE

!==============================================================================

  PUBLIC  :: terra_init
  PUBLIC  :: get_wsnow

CONTAINS

!==============================================================================

!>
!! Performs initialization of soil model TERRA
!!
!! Performs initialization of soil model TERRA. Depending on the initialization mode 
!! chosen, a coldstart- , warmstart-  or warmstart with snow increments (IAU)-initialization
!! is performed
!!
SUBROUTINE terra_init (            & 
                init_mode,         & ! 1 = coldstart, 2 = warmstart, 3 = warmstart with snow increments (IAU)
                nvec,              & ! array dimensions
                ivstart,           & ! start index for computations in the parallel program
                ivend,             & ! end index for computations in the parallel program
                iblock           , & ! number of block
! 4. variables for the time discretization and related variables
!                ke,           & ! nsubs0=1 for single tile, nsubs0=2 for multi-tile
!                nsubs0ubs1,& ! nsubs1=1 for single tile, nsubs1=#tiles+1 for multi-tile
                ke_soil, ke_snow , &
                zmls             , & ! processing soil level structure 
                soiltyp_subs     , & ! type of the soil (keys 0-9)                     --
                rootdp           , & ! depth of the roots                            ( m  )
                plcov            , & ! fraction of surface covered by plants         ( -  )
                t_snow_now       , & ! temperature of the snow-surface               (  K  )
                t_snow_mult_now  , & ! temperature of the snow-surface               (  K  )
                t_rhosnowini     , & ! temperature used for snow density initialization on glaciers (  K  )
                t_s_now          , & ! temperature of the ground surface             (  K  )
                t_s_new          , & ! temperature of the ground surface             (  K  )
                w_snow_now       , & ! water content of snow                         (m H2O)
                h_snow           , & ! snow depth                                   (m H2O)
                rho_snow_now     , & ! snow density                                  (kg/m**3)
                rho_snow_mult_now, & ! snow density                                  (kg/m**3)
                t_so_now         , & ! soil temperature (main level)                 (  K  )
                t_so_new         , & ! soil temperature (main level)                 (  K  )
                w_so_now         , & ! total water conent (ice + liquid water)       (m H20)
                w_so_new         , & ! total water conent (ice + liquid water)       (m H20)
                w_so_ice_now     , & ! ice content                                   (m H20)
                w_so_ice_new     , & ! ice content                                   (m H20)
                wliq_snow_now    , & ! liquid water content in the snow              (m H2O)
                wtot_snow_now    , & ! total (liquid + solid) water content of snow  (m H2O)
                dzh_snow_now       & ! layer thickness between half levels in snow   (  m  )
                                  )

!-------------------------------------------------------------------------------
! Declarations
!-------------------------------------------------------------------------------

  INTEGER, INTENT(IN)  ::  &
                  init_mode    ! 1 = coldstart, 2 = warmstart, 3 = warmstart with snow increments (IAU)
  INTEGER, INTENT(IN)  ::  &
                  nvec,         & ! array dimensions
                  ivstart,      & ! start index for computations in the parallel program
                  ivend,        & ! end index for computations in the parallel program
                  iblock,       & ! number of block
                  ke_soil, ke_snow      
  REAL    (KIND = wp)    , DIMENSION(ke_soil+1), INTENT(IN) :: &
                  zmls            ! processing soil level structure 
  INTEGER, DIMENSION(nvec), INTENT(IN) :: & 
                  soiltyp_subs    ! type of the soil (keys 0-9)                     --
  REAL    (KIND = wp)    , DIMENSION(nvec), INTENT(IN) :: & 
                  rootdp          ! depth of the roots                            ( m  )
  REAL    (KIND = wp)    , DIMENSION(nvec), INTENT(IN) :: &
                  plcov           ! plant coverage                                ( -  )
  REAL    (KIND = wp)    , DIMENSION(nvec), INTENT(INOUT) :: &
                  t_snow_now              ! temperature of the snow-surface (K)
  REAL    (KIND = wp)    , DIMENSION(nvec,0:ke_snow), INTENT(INOUT) :: &
                  t_snow_mult_now      ! temperature of the snow-surface               (  K  )
  REAL    (KIND = wp)    , DIMENSION(nvec), INTENT(IN) :: &
                  t_rhosnowini         ! temperature used for snow density initialization on glaciers (K)
  REAL    (KIND = wp)    , DIMENSION(nvec), INTENT(INOUT) :: &
                  t_s_now              ! temperature of the ground surface             (  K  )
  REAL    (KIND = wp)    , DIMENSION(nvec), INTENT(OUT) :: &
                  t_s_new              ! temperature of the ground surface             (  K  )
  REAL    (KIND = wp)    , DIMENSION(nvec), INTENT(INOUT) :: &
                  w_snow_now       , & ! water content of snow                         (m H2O)
                  rho_snow_now         ! snow density                                  (kg/m**3)
  REAL    (KIND = wp)    , DIMENSION(nvec), INTENT(INOUT) :: &
                  h_snow               ! snow depth                                   (m H2O)
  REAL    (KIND = wp)    , DIMENSION(nvec,ke_snow), INTENT(INOUT) :: &
                  rho_snow_mult_now    ! snow density                                  (kg/m**3)
  REAL    (KIND = wp)    , DIMENSION(nvec,0:ke_soil+1), INTENT(INOUT) :: &
                  t_so_now             ! soil temperature (main level)                 (  K  )
  REAL    (KIND = wp)    , DIMENSION(nvec,0:ke_soil+1), INTENT(OUT) :: &
                  t_so_new             ! soil temperature (main level)                 (  K  )
  REAL    (KIND = wp)    , DIMENSION(nvec,ke_soil+1), INTENT(INOUT) :: &
                  w_so_now         , & ! total water conent (ice + liquid water)       (m H20)
                  w_so_ice_now         ! ice content                                   (m H20)
  REAL    (KIND = wp)    , DIMENSION(nvec,ke_soil+1), INTENT(OUT) :: &
                  w_so_new         , & ! total water conent (ice + liquid water)       (m H20)
                  w_so_ice_new         ! ice content                                   (m H20)
  REAL    (KIND = wp)    , DIMENSION(nvec,ke_snow), INTENT(INOUT) :: &
                  wliq_snow_now    , & ! liquid water content in the snow              (m H2O)
                  wtot_snow_now        ! total (liquid + solid) water content of snow  (m H2O)
  REAL    (KIND = wp)    , DIMENSION(nvec,ke_snow), INTENT(INOUT) :: &
                  dzh_snow_now         ! layer thickness between half levels in snow   (  m  )

!--------------------------------------------------------------------------------
! TERRA Declarations

! Local parameters:
! ----------------
  ! are all used from sfc_terra_data

! Local scalars:
! -------------

  INTEGER ::  &
    kso            , & ! loop index for soil moisture layers           
    ksn,k          , & ! loop index for snow layers
    i,ic           , & ! loop index in x-direction              
    mstyp              ! soil type index

  REAL    (KIND=wp) ::  &
    zpsi0=0.01_wp      ! air entry potential at water saturation (m)

  REAL    (KIND=wp) ::  &
    zaa, zh_snow,       &       ! utility variable
    t_zw_low, t_zw_up, zd, zd1, zd2, zd3, zd4, &
    zw_m_low, zw_m_soil, zw_m_up, zzz,         &
    weight, wso_ice_tolerance, wso_ice_equil, rwso_ice_tolerance

  REAL    (KIND=wp) ::  &
    zroota   (nvec)      , & ! root density profile parameter (1/m)
! External parameters
    zbwt     (nvec)      , & ! root depth (with artificial minimum value)
    zsandf   (nvec)      , & ! mean fraction of sand (weight percent)
    zclayf   (nvec)      , & ! mean fraction of clay (weight percent)
    zsiltf   (nvec)      , & ! mean fraction of silt (weight percent)
    zb_por   (nvec)      , & ! pore size distribution index
    zpsis    (nvec)      , & ! air entry potential (m)
    zw_m     (nvec)          ! maximum of liquid water content  (m)

  INTEGER ::  &
    m_styp   (nvec)      , & ! soil type
    zicount1 (nvec)      , &
    zicount2 (nvec)

  REAL    (KIND=wp) ::  & ! Soil and plant parameters
!
    zadp     (nvec,ke_soil+1), & ! air dryness point
    zporv    (nvec,ke_soil+1), & ! pore volume (fraction of volume)
    zedb     (nvec)              ! utility variable

  REAL    (KIND=wp) ::             &
    zw_snow_old   (nvec)         , & !
    zrho_snow_old (nvec)         , &
    h_snow_fg     (nvec)         , &
    h_snow_incr   (nvec)         , &
    sum_weight    (nvec)         , &
    zhh_snow      (nvec,ke_snow) , & ! depth of the half level snow layers
    zhm_snow      (nvec,ke_snow) , & ! depth of snow main levels
    t_new         (nvec,ke_snow) , &
    rho_new       (nvec,ke_snow) , &
    wl_new        (nvec,ke_snow) , &
    z_old         (nvec,ke_snow) , &
    dz_old        (nvec,ke_snow)

  LOGICAL :: l_redist(nvec)

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Define GPU data region
!------------------------------------------------------------------------------

  ! Subroutine parameters IN
  !$noacc data                                                             &
  !$noacc present(zmls, soiltyp_subs, rootdp, plcov                      ) &

  ! Subroutine parameters INOUT
  !$noacc present(t_snow_now, t_snow_mult_now, t_s_now                   ) &
  !$noacc present(t_rhosnowini                                           ) &
  !$noacc present(w_snow_now, h_snow, rho_snow_now, rho_snow_mult_now    ) &
  !$noacc present(t_so_now, w_so_now, w_so_ice_now, wliq_snow_now        ) &
  !$noacc present(wtot_snow_now, dzh_snow_now                            ) &

  ! Subroutine parameters OUT
  !$noacc present(t_s_new, t_so_new, w_so_new, w_so_ice_new              ) &

  ! Local arrays
  !$noacc present(m_styp, zicount1, zicount2                             ) &
  !$noacc present(zroota, zdzhs, zbwt, zsandf, zsiltf                    ) &
  !$noacc present(zclayf, zb_por, zpsis, zw_m, zadp, zporv, zedb         ) &
  !$noacc present(zw_snow_old, zrho_snow_old, h_snow_fg, h_snow_incr     ) &
  !$noacc present(sum_weight, zhh_snow, zhm_snow, t_new, rho_new         ) &
  !$noacc present(wl_new, z_old, dz_old                                  ) &

  ! Terra data module fields
  !$noacc present(cporv, cadp, csandf, cclayf)

!------------------------------------------------------------------------------
! Begin Subroutine terra_init
!------------------------------------------------------------------------------
!==============================================================================
!  Computation of the diagnostic part I of the soil parameterization scheme
!  In this part, evaporation from the surface and transpiration is calculated.
!  A multi-layer soil water distribution is calculated by simple bulk
!  parameterisation or a Penman-Monteith-type version of evaporation
!  and transpiration which can be used alternatively.
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section I.1: Initializations
!------------------------------------------------------------------------------

! Prepare basic surface properties (for land-points only)

  !$noacc parallel
  !$noacc loop gang vector private(mstyp)
  DO i = ivstart, ivend
    mstyp       = soiltyp_subs(i)           ! soil type
    m_styp(i)   = mstyp                     ! array for soil type
    zporv (i,:) = cporv(mstyp)              ! pore volume
    zadp  (i,:) = cadp (mstyp)              ! air dryness point
    zbwt  (i)   = MAX(0.001_wp,rootdp(i))   ! Artificial minimum value
                                            ! for root depth
    zroota(i)   = 3.0_wp/zbwt(i)            ! root density profile parameter (1/m)
                                            ! zroota=0. creates the original TERRA_LM
                                            ! version with constant root density
                                            ! for root depth

    ! Arrays for soil water freezing/melting
    zsandf(i)   = csandf(mstyp)
    zclayf(i)   = cclayf(mstyp)
    zsiltf(i)   = 100.0_wp -csandf(mstyp)-cclayf(mstyp) ! Residuum of sand and clay
    zpsis (i)   = -zpsi0 * 10.0_wp**(1.88_wp-0.013_wp*zsandf(i))
    zb_por(i)   = 2.91_wp+ 0.159_wp*zclayf(i)
    zedb  (i)   = 1.0_wp/zb_por(i)
  ENDDO
  !$noacc end parallel
 
  ! Further parameters for soil water freezing/melting
  t_zw_up  = 270.15_wp     ! temp -3 degC
  t_zw_low = 233.15_wp     ! temp -40 degC
  zd  = LOG((T_ref_ice-(t_zw_low-t0_melt))/T_star_ice)
  zd1 = EXP(b_sand*zd)
  zd2 = EXP(b_clay*zd)
  zd3 = EXP(b_silt*zd)
  zd4 = EXP(b_org*zd)

  ! Tolerances for the allowed deviation of w_so_ice from its equilibrium value
  wso_ice_tolerance  = 1.05_wp! 5% deviation allowed
  rwso_ice_tolerance = 1.0_wp/wso_ice_tolerance


! For ntstep=0 : Some preparations
! ================================

  !$noacc kernels
  w_so_new(:,:) = w_so_now(:,:)
  !$noacc end kernels

! Provide for a soil moisture 1 % above air dryness point, reset soil
! moisture to zero in case of ice and rock
  !$noacc parallel
  DO kso   = 1,ke_soil+1
    !$noacc loop gang vector private(mstyp)
    DO i = ivstart, ivend
      IF (m_styp(i) >= 3) THEN
        ! soil list
        w_so_now (i,kso) = MAX(w_so_now(i,kso), 1.01_wp*zadp(i,kso)*zdzhs(kso) )
        w_so_new (i,kso) = MAX(w_so_new(i,kso), 1.01_wp*zadp(i,kso)*zdzhs(kso) )
      ELSE
        ! rock and ice list
        w_so_now (i,kso) = 0.0_wp
        w_so_new (i,kso) = 0.0_wp
      ENDIF
    END DO
  END DO
  !$noacc end parallel


! adjust temperature profile in lower soil layers, if temperature of first soil
! layer was reduced due to the detection of snow (e.g. in the analysis)
! loop over grid points
  !$noacc parallel
  !$noacc loop gang vector
  DO i = ivstart, ivend
    IF (w_snow_now(i) <= eps_soil) THEN
      ! spurious snow is removed
      t_snow_now(i)=t_so_now(i,0)
      w_snow_now(i)= 0.0_wp
    ELSE
      t_snow_now(i) = MIN( t0_melt, t_snow_now(i))
    ENDIF

    t_s_now(i)  = t_so_now(i,0)
    t_s_new(i)  = t_so_now(i,0)

    ! Set level 1 to level 0 for t_so for every landpoint
    t_so_now(i,1) = t_so_now(i,0)
    t_so_new(i,1) = t_so_now(i,0)
    t_so_new(i,0) = t_so_now(i,0)
  END DO
  !$noacc end parallel

  IF(lmulti_snow) THEN
    !$noacc parallel
    DO ksn = 0, ke_snow
      !$noacc loop gang vector
      DO i = ivstart, ivend
        IF (w_snow_now(i) <= eps_soil) THEN
          t_snow_mult_now   (i,ksn) = t_so_now(i,0)
          IF (ksn > 0) THEN
            wliq_snow_now     (i,ksn) = 0.0_wp
            wtot_snow_now     (i,ksn) = 0.0_wp
            rho_snow_mult_now (i,ksn) = 0.0_wp
            dzh_snow_now      (i,ksn) = 0.0_wp
          ENDIF
        ELSE
          t_snow_mult_now(i,ksn) = MIN( t0_melt, t_snow_mult_now(i,ksn))
        ENDIF
      END DO
    END DO
    !$noacc end parallel
  ENDIF

  ! Note that a few parts need to be shifted to the .NOT. is_coldstart branch.
  ! I.e. those which are dealing with the contributions from the analysis.
  ! Currently, the multi layer soil model starts from the FG.
  IF ( init_mode == 1 ) THEN

    ! Initialization of 
    ! - snow density (if necessary)
    ! - multi layer snow variables (if necessary)
    !   * wtot_snow, dzh_snow, wliq_snow
    ! - soil ice water content w_so_ice
    ! - snow height h_snow

!   Initialization of snow density, if necessary
!   --------------------------------------------

    IF (.NOT. lana_rho_snow) THEN

       IF(lmulti_snow) THEN

         ksn = 1

         !$noacc parallel
         !$noacc loop gang vector
         DO i = ivstart, ivend
           rho_snow_mult_now(i,ksn) = 250.0_wp    ! average initial density
           t_snow_mult_now  (i,ksn) = t_snow_now(i)
           wliq_snow_now    (i,ksn) = 0.0_wp
           ! Limit top layer to max_toplaydepth
           dzh_snow_now     (i,ksn) = MIN(max_toplaydepth, w_snow_now(i)/REAL(ke_snow,wp) &
                  &                            *rho_w/rho_snow_mult_now(i,ksn))
           wtot_snow_now    (i,ksn) = dzh_snow_now(i,ksn)*rho_snow_mult_now(i,ksn)/rho_w
           t_snow_mult_now  (i,0  ) = t_snow_now(i)
         ENDDO
         !$noacc end parallel

         k = MIN(2, ke_snow-1)
         !$noacc parallel
         DO ksn = 2, ke_snow
           !$noacc loop gang vector
           DO i = ivstart, ivend
             rho_snow_mult_now(i,ksn) = 250.0_wp    ! average initial density
             t_snow_mult_now  (i,ksn) = t_snow_now(i)
             wliq_snow_now    (i,ksn) = 0.0_wp
             ! distribute the remaining snow equally among the layers
             IF (k == ksn) THEN
               dzh_snow_now   (i,ksn) = MIN(8._wp*max_toplaydepth, (w_snow_now(i)-wtot_snow_now(i,1)) &
                    &                            /REAL(ke_snow-1,wp)*rho_w/rho_snow_mult_now(i,ksn) )
             ELSE
               dzh_snow_now   (i,ksn) = (w_snow_now(i)-wtot_snow_now(i,k))/REAL(ke_snow-k,wp) &
                    &                            *rho_w/rho_snow_mult_now(i,ksn)
             ENDIF
             wtot_snow_now    (i,ksn) = dzh_snow_now(i,ksn)*rho_snow_mult_now(i,ksn)/rho_w
           ENDDO
         ENDDO
         !$noacc end parallel
       ELSE
         !$noacc parallel
         !$noacc loop gang vector
         DO i = ivstart, ivend
           rho_snow_now(i) = 250.0_wp    ! average initial density
           IF (l2lay_rho_snow) THEN
             !US is second dim really  2??
             rho_snow_mult_now(i,1) = 250.0_wp
             rho_snow_mult_now(i,2) = 250.0_wp
           ENDIF
         ENDDO
         !$noacc end parallel
       ENDIF

!US so far for COSMO, the rest is ICON

    ELSE

      IF (lmulti_snow) THEN

        ! The sum of wtot_snow, i.e. the "old" total water equivalent depth
        !$noacc parallel
        !$noacc loop gang vector
        DO i = ivstart, ivend
          zw_snow_old  (i) = 0.0_wp
          zrho_snow_old(i) = 0.0_wp
          zicount1(i)      = 0
          zicount2(i)      = 0
          sum_weight(i)    = 0.0_wp
        ENDDO
        !$noacc end parallel

        !$noacc parallel
        DO ksn = 1, ke_snow
          !$noacc loop gang vector
          DO i = ivstart, ivend
            zw_snow_old(i) = zw_snow_old(i) + wtot_snow_now(i,ksn)
            zrho_snow_old(i) = zrho_snow_old(i) + dzh_snow_now(i,ksn)
            IF (wtot_snow_now(i,ksn) > eps_soil) zicount1(i) = zicount1(i)+1
            IF (dzh_snow_now(i,ksn) > eps_soil)  zicount2(i) = zicount2(i)+1
          END DO
        END DO
        !$noacc end parallel

        !The "old" snow density
        !$noacc parallel
        !$noacc loop gang vector
        DO i = ivstart, ivend
          zrho_snow_old(i) = zw_snow_old(i) / MAX(zrho_snow_old(i),1.0E-09_wp) *rho_w
        END DO
        !$noacc end parallel

        k = MIN(2, ke_snow-1)
        !$noacc parallel
        DO ksn = 1, ke_snow
          !$noacc loop gang vector
          DO i = ivstart, ivend
            IF(zicount1(i).EQ.ke_snow .AND. zicount2(i).EQ.ke_snow) THEN
              rho_snow_mult_now(i,ksn) = wtot_snow_now(i,ksn)/dzh_snow_now(i,ksn)*rho_w
              wtot_snow_now(i,ksn) = wtot_snow_now(i,ksn)*w_snow_now(i)/zw_snow_old(i)
              dzh_snow_now (i,ksn) = dzh_snow_now (i,ksn)*w_snow_now(i)/zw_snow_old(i)
              wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn)*w_snow_now(i)/zw_snow_old(i)
            ELSE
              IF(rho_snow_now(i) .EQ. 0._wp) rho_snow_now(i) = 250._wp
              IF (ksn == 1) THEN ! Limit top layer to max_toplaydepth
                IF (t_snow_now(i) > 268.15_wp) THEN
                  rho_snow_mult_now(i,ksn) = MIN(400._wp,rho_snow_now(i))
                ELSE ! at very low temperatures, the near-surface snow is usually not subject to strong densification
                  rho_snow_mult_now(i,ksn) = MIN(rho_snow_now(i),400._wp+6._wp* &
                    (MAX(228.15_wp,t_snow_now(i))-268.15_wp)                    )
                ENDIF
                dzh_snow_now(i,ksn) = MIN(max_toplaydepth, w_snow_now(i)/REAL(ke_snow,wp) &
                  &                  *rho_w/rho_snow_mult_now(i,ksn))
              ELSE ! distribute the remaining snow equally among the layers
                rho_snow_mult_now(i,ksn) = rho_snow_now(i)
                IF (ksn == k) THEN
                  dzh_snow_now(i,ksn) = MIN(8._wp*max_toplaydepth, (w_snow_now(i)-wtot_snow_now(i,1)) &
                      &                 /REAL(ke_snow-1,wp)*rho_w/rho_snow_mult_now(i,ksn) )
                ELSE
                  dzh_snow_now(i,ksn) = (w_snow_now(i)-sum_weight(i))/REAL(ke_snow-k,wp) &
                      &                 *rho_w/rho_snow_mult_now(i,ksn)
                ENDIF
              ENDIF
              wtot_snow_now    (i,ksn) = dzh_snow_now(i,ksn)*rho_snow_mult_now(i,ksn)/rho_w
              sum_weight(i) = sum_weight(i) + wtot_snow_now(i,ksn)
              wliq_snow_now    (i,ksn) = 0.0_wp
            END IF
          END DO
        END DO
        !$noacc end parallel

        !$noacc parallel
        !$noacc loop gang vector
        DO i = ivstart, ivend
          t_snow_mult_now  (i,0  ) = t_snow_now(i)
        END DO
        !$noacc end parallel
      ELSE
        !$noacc parallel
        !$noacc loop gang vector
        DO i = ivstart, ivend
          IF(rho_snow_now(i) .EQ. 0._wp) rho_snow_now(i) = 250._wp
          IF (l2lay_rho_snow) THEN
            IF(rho_snow_mult_now(i,1) == 0._wp) rho_snow_mult_now(i,1) = rho_snow_now(i)
            IF(rho_snow_mult_now(i,2) == 0._wp) rho_snow_mult_now(i,2) = rho_snow_now(i)
          ENDIF
        END DO
        !$noacc end parallel
      END IF
    ENDIF


!   Initialization of soil ice content
!   ----------------------------------
!   Generally the combination lmelt=.true., lmelt_var=.true. has to be used.
!   lmelt_var=.false. (soil water freezing at 273.15K) should be used only
!   for special investigations because at model start a partial frozen layer
!   is only considered if the soil ice water content is read in from
!   a pre run.
!   ----------------------------------
    IF (lmelt .AND. .NOT. lmelt_var) THEN
      !$noacc parallel
      DO kso   = 1,ke_soil+1
        !$noacc loop gang vector
        DO i = ivstart, ivend
          w_so_ice_now(i,kso) = 0.0_wp
          w_so_ice_new(i,kso) = 0.0_wp

          IF (t_so_now(i,kso) < t0_melt) THEN
            w_so_ice_now(i,kso) = w_so_now(i,kso)
            w_so_ice_new(i,kso) = w_so_now(i,kso)
          END IF ! t_so(kso) < t0_melt
        END DO
      END DO
      !$noacc end parallel
    END IF           ! lmelt .AND. .NOT. lmelt_var
    IF(lmelt .AND. lmelt_var) THEN
      !$noacc parallel
      DO kso   = 1,ke_soil+1
        !$noacc loop gang vector
        DO i = ivstart, ivend
          IF (t_so_now(i,kso) < (t0_melt-eps_temp)) THEN 

            zaa    = g*zpsis(i)/lh_f

            ! J. Helmert: Soil ice parameterization according to K. Schaefer and Jafarov, E.,2016,
            !                                                    doi:10.5194/bg-13-1991-2016

            ! Liq. water content at -3 degC
            zw_m_up =  zporv(i,kso)*zdzhs(kso)*EXP(-zedb(i)*LOG((t_zw_up - t0_melt)/(t_zw_up*zaa)) )

            ! Determine liq. water content at -40 degC
            zw_m_soil = 0.01_wp*(zsandf(i)*zd1 + zclayf(i)*zd2 + zsiltf(i)*zd3)

            ! Scale soil ice content with organic soil horizon.
            ! should decrease the root zone liquid water content of frozen soil for low temperatures significantly!
            IF(zmls(kso) < rootdp(i)) THEN
              zzz = plcov(i)*(rootdp(i)-zmls(kso))/rootdp(i)
              zw_m_low = zporv(i,kso)*zdzhs(kso)*(zzz*zd4 + (1.0_wp-zzz)*zw_m_soil)
            ELSE
              zw_m_low = zporv(i,kso)*zdzhs(kso)*zw_m_soil
            END IF

            IF (t_so_now(i,kso) < t_zw_low) THEN
              zw_m(i) = zw_m_low
            ELSE IF (t_so_now(i,kso) < t_zw_up) THEN ! Logarithmic Interpolation between -3 degC and -40 degC 
              zw_m(i) = zw_m_low*EXP((t_so_now(i,kso) - t_zw_low)*(LOG(zw_m_up) - LOG(zw_m_low))/(t_zw_up-t_zw_low))
            ELSE
              zw_m(i) = zporv(i,kso)*zdzhs(kso)*EXP(-zedb(i)*LOG((t_so_now(i,kso)-t0_melt)/(t_so_now(i,kso)*zaa)) )
            END IF

            w_so_ice_now(i,kso) = MAX (0.0_wp,w_so_now(i,kso) - zw_m(i))
            w_so_ice_new(i,kso) = MAX (0.0_wp,w_so_now(i,kso) - zw_m(i))

          ELSE ! ensure that w_so_ice is zero

            w_so_ice_now(i,kso) = 0.0_wp
            w_so_ice_new(i,kso) = 0.0_wp

          END IF
        END DO
      END DO
      !$noacc end parallel
    END IF           ! lmelt .AND. lmelt_var


!   Initialization of the local array of the grid in snow
!   -----------------------------------------------------

    IF(lmulti_snow) THEN
      !$noacc parallel
      !$noacc loop gang vector
      DO i = ivstart, ivend
        sum_weight(i) = 0.0_wp
        h_snow    (i) = 0.0_wp
      END DO
      !$noacc end parallel

      !$noacc parallel
      DO ksn = 1, ke_snow
        !$noacc loop gang vector
        DO i = ivstart, ivend
          h_snow(i) = h_snow(i) + dzh_snow_now(i,ksn)
        END DO
      END DO
      !$noacc end parallel

      !$noacc parallel
      DO ksn = 1,ke_snow
        !$noacc loop gang vector
        DO i = ivstart, ivend
          IF(h_snow(i) > 0.0_wp) THEN
            t_snow_mult_now(i,ksn) = t_snow_now(i) + &
              (t_s_now(i)-t_snow_now(i))*(sum_weight(i)+dzh_snow_now(i,ksn)*0.5)/h_snow(i)
            sum_weight(i) = sum_weight(i) + dzh_snow_now(i,ksn)
          END IF
        END DO
      END DO
      !$noacc end parallel
    ELSE
      !$noacc parallel
      !$noacc loop gang vector
      DO i = ivstart, ivend
        h_snow(i) = w_snow_now(i)/MAX(crhosmin_ml,rho_snow_now(i))*rho_w
      END DO
      !$noacc end parallel
    END IF

  ELSE  ! init_mode > 1, i.e. warmstart (i.e. assimilation cycle and forecast)
        ! in this case, w_snow and rho_snow are read from the first guess, and h_snow is updated by the snow analysis;
        ! in the IAU mode, the snow depth increment (h_snow_incr) is directly provided by the snow analysis


    IF ( init_mode == 2 ) THEN  ! snow depth increment needs to be computed
      DO i = ivstart, ivend
        h_snow_fg(i)   = w_snow_now(i)/MAX(crhosmin_ml,rho_snow_now(i))*rho_w
        ! Self-consistency fix needed because of possible GRIB truncation errors
        IF (lmulti_snow .AND. h_snow_fg(i) > eps_soil) THEN
          dzh_snow_now(i,2:ke_snow) = MAX(dzh_snow_now(i,1),dzh_snow_now(i,2:ke_snow))
          zh_snow       = SUM(dzh_snow_now(i,1:ke_snow))
          w_snow_now(i) = w_snow_now(i) * zh_snow/h_snow_fg(i)
          h_snow_fg(i)  = zh_snow
        ELSE IF (lmulti_snow) THEN
          dzh_snow_now(i,:) = 0._wp
        ENDIF
        h_snow_incr(i) = h_snow(i) - h_snow_fg(i)
      ENDDO
    ELSE IF ( init_mode == 3 ) THEN  ! snow depth increment is provided from snow analysis
      IF (lmulti_snow) THEN
        DO i = ivstart, ivend
          ! Self-consistency fix needed because of possible GRIB truncation errors
          IF (dzh_snow_now(i,1) > eps_soil) THEN
            dzh_snow_now(i,2:ke_snow) = MAX(dzh_snow_now(i,1),dzh_snow_now(i,2:ke_snow))
            h_snow_fg(i) = SUM(dzh_snow_now(i,1:ke_snow))
          ELSE
            dzh_snow_now(i,:) = 0._wp
            h_snow_fg(i)      = 0._wp
          ENDIF
          h_snow_incr(i) = h_snow(i) - h_snow_fg(i)
        ENDDO
      ELSE
        DO i = ivstart, ivend
          h_snow_fg(i) = h_snow(i)
          h_snow_incr(i) = 0._wp
        ENDDO
      ENDIF
    ENDIF

    ! Provisional fix to get rid of isolated extremely cold points in parallel routine (2014-10-16)
    ! This is applied only to the lowest prognostic soil layer, i.e. ke_soil
    DO i = ivstart, ivend
      IF (MIN(t_so_now(i,ke_soil+1),t_so_now(i,ke_soil-1)) - t_so_now(i,ke_soil) > 5._wp) THEN
        t_so_now(i,ke_soil) = MIN(t_so_now(i,ke_soil+1),t_so_now(i,ke_soil-1)) - 5._wp
      ENDIF
      IF (t_so_now(i,ke_soil) - MAX(t_so_now(i,ke_soil+1),t_so_now(i,ke_soil-1)) > 5._wp) THEN
        t_so_now(i,ke_soil) = MAX(t_so_now(i,ke_soil+1),t_so_now(i,ke_soil-1)) + 5._wp
      ENDIF
    ENDDO

    ! Ensure that w_so_ice stays within 5% of its equilibrium value
    ! In addition, it must not exceed w_so
    IF(lmelt .AND. lmelt_var) THEN
      DO kso   = 1,ke_soil+1
        DO i = ivstart, ivend
          IF (t_so_now(i,kso) < (t0_melt-eps_temp)) THEN 

            zaa     = g*zpsis(i)/lh_f

            ! J. Helmert: Soil ice parameterization according to K. Schaefer and Jafarov, E.,2016,
            !                                                    doi:10.5194/bg-13-1991-2016

            ! Liq. water content at -3 degC
            zw_m_up =  zporv(i,kso)*zdzhs(kso)*EXP(-zedb(i)*LOG((t_zw_up - t0_melt)/(t_zw_up*zaa)) )

            ! Determine liq. water content at -40 degC
            zw_m_soil = 0.01_wp*(zsandf(i)*zd1 + zclayf(i)*zd2 + zsiltf(i)*zd3)

            ! Scale soil ice content with organic soil horizon.
            ! should decrease the root zone liquid water content of frozen soil for low temperatures significantly!
            IF(zmls(kso) < rootdp(i)) THEN
              zzz = plcov(i)*(rootdp(i)-zmls(kso))/rootdp(i)
              zw_m_low = zporv(i,kso)*zdzhs(kso)*(zzz*zd4 + (1.0_wp-zzz)*zw_m_soil)
            ELSE
              zw_m_low = zporv(i,kso)*zdzhs(kso)*zw_m_soil
            END IF

            IF (t_so_now(i,kso) < t_zw_low) THEN
              zw_m(i) = zw_m_low
            ELSE IF (t_so_now(i,kso) < t_zw_up) THEN ! Logarithmic Interpolation between -3 degC and -40 degC 
              zw_m(i) = zw_m_low*EXP((t_so_now(i,kso) - t_zw_low)*(LOG(zw_m_up) - LOG(zw_m_low))/(t_zw_up-t_zw_low))
            ELSE
              zw_m(i) = zporv(i,kso)*zdzhs(kso)*EXP(-zedb(i)*LOG((t_so_now(i,kso)-t0_melt)/(t_so_now(i,kso)*zaa)) )
            END IF

            wso_ice_equil = MAX (0.0_wp,w_so_now(i,kso) - zw_m(i))
            w_so_ice_now(i,kso) = MIN(w_so_ice_now(i,kso),  wso_ice_tolerance*wso_ice_equil, w_so_now(i,kso))
            w_so_ice_now(i,kso) = MAX(w_so_ice_now(i,kso), rwso_ice_tolerance*wso_ice_equil)
          ELSE ! ensure that w_so_ice is zero
            w_so_ice_now(i,kso) = 0.0_wp
          END IF
        END DO
      END DO
    END IF

    ! For the single-layer snow case, the prognostic variables now need to be updated. This requires appropriate
    ! assumptions about the new snow density

    IF(lmulti_snow) THEN ! update of prognostic variables of multi-layer snow model
      l_redist(:) = .FALSE.
      DO i = ivstart, ivend
        IF (h_snow(i) <= eps_soil) THEN  ! snow has disappeared or was already absent
          rho_snow_now(i) = 250._wp
          w_snow_now(i)   = 0._wp

          dzh_snow_now (i,1:ke_snow) = 0._wp
          wtot_snow_now(i,1:ke_snow) = 0._wp
          rho_snow_mult_now(i,1:ke_snow) = 0._wp
          wliq_snow_now(i,1:ke_snow) = 0._wp
          t_snow_mult_now(i,1:ke_snow) = t_snow_now(i)
        ELSE IF (h_snow_fg(i) <= eps_soil .AND. h_snow(i) > eps_soil) THEN 

          ! snow analysis has either reestablished a snow cover that had been erroneously melted away by the model
          ! or generated a new snow cover that was missed by the model due to a lack of precipitation (or wrong phase of precip)
          ! ** snow density should be distinguished between new snow and reestablished old snow **
          rho_snow_now(i) = 250._wp ! needs to be improved
          w_snow_now(i)   = h_snow(i)*rho_snow_now(i)/rho_w
          t_snow_now(i)   = MIN(t_snow_now(i),t0_melt)

          rho_snow_mult_now(i,1:ke_snow) = rho_snow_now(i)
          dzh_snow_now (i,1:ke_snow) = h_snow(i)/REAL(ke_snow,wp)
          wtot_snow_now(i,1:ke_snow) = w_snow_now(i)/REAL(ke_snow,wp)
          wliq_snow_now(i,1:ke_snow) = 0._wp
          t_snow_mult_now(i,1:ke_snow) = t_snow_now(i)
          IF (h_snow(i) > REAL(ke_snow,wp)*max_toplaydepth) l_redist(i) = .TRUE.
        ELSE IF (h_snow_fg(i) > eps_soil .AND. h_snow_incr(i) < 0._wp) THEN

          ! a snow cover is present but reduced by the snow analysis increment
          ! in this case, we retain the snow densities but adjust the water equivalents
          l_redist(i) = .TRUE.
          dz_old(i,1) = MAX(dzh_snow_now (i,1),eps_soil/REAL(ke_snow,wp))
          dzh_snow_now (i,1) = dzh_snow_now (i,1) + h_snow_incr(i)
          wtot_snow_now(i,1) = wtot_snow_now(i,1)*dzh_snow_now (i,1)/dz_old(i,1)
        ELSE
          ! a snow cover is present and increased by the snow analysis
          ! ** snow density should be distinguished between new snow and reestablished old snow **
          ! IF (w_snow_incr(i) > 0._wp) THEN ! assume temperature-dependent fresh snow density - maybe air temp should be used here?!
          !   rho_new(i,1) = MAX(50._wp,125._wp-10._wp*(t0_melt-t_snow_mult_now(i,1)))
          ! ELSE
          !   rho_new(i,1) = rho_snow_mult_now(i,1)
          ! ENDIF
          l_redist(i) = .TRUE.
          rho_new(i,1) = rho_snow_mult_now(i,1)
          wtot_snow_now(i,1) = wtot_snow_now(i,1) + h_snow_incr(i)*rho_new(i,1)/rho_w
          dzh_snow_now (i,1) = dzh_snow_now (i,1) + h_snow_incr(i)
          rho_snow_mult_now(i,1) = rho_w*wtot_snow_now(i,1)/dzh_snow_now (i,1)
        ENDIF
      ENDDO

      ! Redistribute snow cover among layers where necessary
      sum_weight(:) = 0._wp
      DO ksn = ke_snow,1,-1
        DO i = ivstart, ivend
          IF (l_redist(i)) THEN
            dz_old(i,ksn) = dzh_snow_now(i,ksn)
            z_old(i,ksn) = -sum_weight(i) - dzh_snow_now(i,ksn)/2._wp
            sum_weight(i) = sum_weight(i) + dzh_snow_now(i,ksn)
          ENDIF
        END DO
      END DO

      k = MIN(2, ke_snow-1)
      DO ksn = 1,ke_snow
        DO i = ivstart, ivend
          IF (l_redist(i)) THEN

            IF (ksn == 1) THEN ! Limit top layer to max_toplaydepth
              zhh_snow(i,ksn) = -MAX( h_snow(i)-max_toplaydepth, h_snow(i)/ke_snow*(ke_snow-ksn) )
              zhm_snow(i,ksn) = (-h_snow(i) + zhh_snow(i,ksn))/2._wp   
              dzh_snow_now(i,ksn) = zhh_snow(i,ksn) + h_snow(i)    !layer thickness betw. half levels of uppermost snow layer
            ELSE IF (ksn == 2 .AND. ke_snow > 2) THEN ! Limit second layer to 8*max_toplaydepth
              zhh_snow(i,ksn) = MIN( 8._wp*max_toplaydepth+zhh_snow(i,1), zhh_snow(i,1)/(ke_snow-1)*(ke_snow-ksn) )
            ELSE ! distribute the remaining snow equally among the layers
              zhh_snow(i,ksn) = zhh_snow(i,k)/(ke_snow-k)*(ke_snow-ksn)
            ENDIF

            t_new  (i,ksn) = 0.0_wp
            rho_new(i,ksn) = 0.0_wp
            wl_new (i,ksn) = 0.0_wp

            IF(dz_old(i,ksn) > 0._wp .AND. rho_snow_mult_now(i,ksn) > 0._wp) THEN
              wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn)/dz_old(i,ksn)
            END IF

          ENDIF
        END DO
      END DO

      DO ksn = 2, ke_snow
        DO i = ivstart, ivend
          IF(l_redist(i)) THEN
            zhm_snow(i,ksn) = (zhh_snow(i,ksn) + zhh_snow(i,ksn-1))/2._wp
            dzh_snow_now(i,ksn) = zhh_snow(i,ksn) - zhh_snow(i,ksn-1) ! layer thickness betw. half levels
          END IF
        END DO
      END DO

      ! mass-weighted recomputation of temperatures and densities
      DO ksn = ke_snow,1,-1 
        DO k = ke_snow,1,-1
          DO i = ivstart, ivend
            IF(l_redist(i)) THEN
         
              weight = MAX(MIN(z_old(i,k)+dz_old(i,k)/2._wp,zhm_snow(i,ksn) + dzh_snow_now(i,ksn)/2._wp)-   &
                       MAX(z_old(i,k)-dz_old(i,k)/2._wp, zhm_snow(i,ksn)-dzh_snow_now(i,ksn)/2._wp),0._wp) &
                       &/dzh_snow_now(i,ksn)
      
              t_new  (i,ksn) = t_new  (i,ksn) + t_snow_mult_now  (i,k)*weight
              rho_new(i,ksn) = rho_new(i,ksn) + rho_snow_mult_now(i,k)*weight
              wl_new (i,ksn) = wl_new (i,ksn) + wliq_snow_now    (i,k)*weight
            END IF
          END DO
        END DO
      END DO

      DO ksn = ke_snow,1,-1
        DO i = ivstart, ivend
          IF(l_redist(i)) THEN
            t_snow_mult_now  (i,ksn) = t_new  (i,ksn)
            rho_snow_mult_now(i,ksn) = rho_new(i,ksn)
            wtot_snow_now    (i,ksn) = rho_new(i,ksn)*dzh_snow_now(i,ksn)/rho_w
            wliq_snow_now    (i,ksn) = wl_new (i,ksn)*dzh_snow_now(i,ksn)
          END IF
        END DO   
      END DO


      ! Re-diagnose integrated/averaged snow-fields, in order to 
      ! have a consistent state.
      !
      h_snow      (ivstart:ivend) = 0.0_wp
      w_snow_now  (ivstart:ivend) = 0.0_wp
      rho_snow_now(ivstart:ivend) = 0.0_wp
      t_snow_now  (ivstart:ivend) = 0.0_wp
      ! Re-diagnose h_snow, w_snow, rho_snow, t_snow
      DO ksn = ke_snow,1,-1
        DO i = ivstart, ivend
          h_snow(i)    = h_snow(i)     + dzh_snow_now(i,ksn)
          w_snow_now(i)= w_snow_now(i) + wtot_snow_now(i,ksn)
        END DO
      ENDDO

      DO ksn = ke_snow,1,-1
        DO i = ivstart, ivend
          IF (h_snow(i) > eps_soil) THEN
            rho_snow_now(i) = rho_snow_now(i) + rho_snow_mult_now(i,ksn) * dzh_snow_now(i,ksn)/h_snow(i)
          ELSE
            rho_snow_now(i) = 250._wp
          ENDIF
        END DO
      ENDDO

      DO i = ivstart, ivend
        t_snow_now(i) = t_snow_mult_now(i,1)
      ENDDO


    ELSE  ! single snow layer

      ! diagnose w_snow and rho_snow from analyzed h_snow
      ! h_snow is taken from analysis, while up to this point 
      ! rho_snow_now and w_snow_now contain the first guess.
      CALL get_wsnow(h_snow,         &  ! in
                     rho_snow_now,   &  ! inout
                     t_rhosnowini,   &  ! in
                     ivstart, ivend, &  ! in
                     soiltyp_subs,   &  ! in
                     w_snow_now      )  ! out

      DO i = ivstart, ivend
        IF (w_snow_now(i) <= eps_soil) THEN
          ! spurious snow is removed
          t_snow_now(i)= t_so_now(i,0)
          w_snow_now(i)= 0.0_wp
        ELSE
          t_snow_now(i) = MIN( t0_melt, t_snow_now(i))
        ENDIF
      ENDDO

    ENDIF ! multi snow

  ENDIF  ! init_mode

!$noacc end data

! End of timestep 0 preparations
! ==============================

  END SUBROUTINE terra_init

!==============================================================================

!>
!! Transform snow depth to snow water equivalent
!!
!! Description:
!!  Snow depth observations are reported in m and this is also the analysis
!!  variable. The prognostic model variables describing the evolution of
!!  snow are snow water eqivalent and snow density.
!!  This subroutine calculates snow water equivalent from snow depth.
!!  For non-permanent ice points, hsnow is limited to 4m. w_snow is limited 
!!  accordingly, using rho_snow from the first guess.
!!
!! Method:
!!  Simple linear relation between snow water and snow depth using snow density

SUBROUTINE get_wsnow(h_snow, rho_snow, t_snow, istart, iend, soiltyp, w_snow)

  REAL(wp), INTENT(INOUT) :: h_snow(:)    ! snow depth            [m]
  REAL(wp), INTENT(INOUT) :: rho_snow(:)  ! snow density          [kg/m**3]
  REAL(wp), INTENT(IN)    :: t_snow(:)    ! snow temperature      [K]
  REAL(wp), INTENT(INOUT) :: w_snow(:)    ! snow water equivalent [m H2O]
  INTEGER,  INTENT(IN)    :: istart       ! start index
  INTEGER,  INTENT(IN)    :: iend         ! end index
  INTEGER,  INTENT(IN)    :: soiltyp(:)   ! type of soil (1-9)

  ! local
  INTEGER :: jc         ! loop index

  REAL(wp), PARAMETER :: rho_snw_default=250._wp
  !--------------------------------------------------------------------

  ! calculate water equivalent from snow depth
  !
  DO jc = istart, iend
    IF (rho_snow(jc) < 1._wp) rho_snow(jc)=rho_snw_default      ! average initial density
    IF (h_snow(jc) > 0._wp) THEN
      ! multiply analysed snow depth [m] by (first-guess) density 
      ! to get water equivalent in [m H2O]
      IF (soiltyp(jc) /= 1) THEN       ! 1=ice
        ! limit snow depth to 4m for non-glacier points
        h_snow(jc) = MIN(h_snow(jc), 4._wp)
        w_snow(jc) = h_snow(jc) * rho_snow(jc)/rho_w
      ELSE
        ! limit snow depth to at least 1 m and set snow density to the equilibrium value for the current snow temperature
        h_snow(jc)   = MAX(h_snow(jc), 1._wp)
        rho_snow(jc) = crhosmax_tmin + MAX(0._wp,(MIN(t0_melt,t_snow(jc))-csnow_tmin)/(t0_melt-csnow_tmin)) &
                       * (crhosmax_ml-crhosmax_tmin)
        w_snow(jc)   = h_snow(jc) * rho_snow(jc)/rho_w
      ENDIF  ! soiltyp
    ENDIF
  ENDDO


END SUBROUTINE get_wsnow

!==============================================================================

!!$!
!!$!+ Transform snow depth to snow water equivalent and vice versa
!!$!
!!$      SUBROUTINE snw_tran (snw,rho,ie,je,imode)
!!$!
!!$! Description:
!!$!  Snow depth observations are reported in m and this is also the analysis
!!$!  variable. The prognostic model variables describing the evolution of
!!$!  snow are  snow water eqivalent and snow density.
!!$!  This subroutine calculates snow depth from water eqivalent and vice versa.
!!$!
!!$! Method:
!!$!  Simple linear relation between snow water and snow depth using snow density
!!$!
!!$! Input  files: none
!!$! Output files: none
!!$!
!!$!=======================================================================
!!$!
!!$! Declarations:
!!$!
!!$! Modules used:
!!$      USE data_all, ONLY:
!!$     & rhosnw
!!$
!!$      IMPLICIT NONE
!!$!=======================================================================
!!$!
!!$! Subroutine / Function arguments
!!$! Scalar arguments with intent(in):
!!$      INTEGER, INTENT(IN) ::  ie, je   ! Dimensions of snws
!!$      INTEGER, INTENT(IN) ::  imode  
!!$!     imode= 1: calculate water equivalent from snow depth
!!$!     imode=-1: calculate snow depth from water equivalent
!!$
!!$! Array arguments with intent(inout):
!!$      REAL, INTENT(INOUT) ::   snw (ie, je), rho(ie,je)
!!$!
!!$! Local Type Definitions:
!!$
!!$!
!!$! Local Parameters:
!!$
!!$!
!!$! Local Scalars:
!!$      INTEGER i,j
!!$      REAL zcrit, zwcrit
!!$!
!!$! Local Arrays:
!!$
!!$!=======================================================================
!!$! Include statements: none
!!$!
!!$!- End of header
!!$!------------------------------------------------------------------------
!!$      if (imode == 1) then
!!$!     calculate water equivalent from snow depth
!!$        do j = 1,je
!!$        do i = 1,ie
!!$          if (snw(i,j) > 0.) then
!!$             if (rho(i,j) < 1.) rho(i,j)=rhosnw
!!$!           scale analysed snow depth [mm] to get in m and multiply it by density 
!!$            snw(i,j) = snw(i,j)/1000.*rho(i,j)
!!$          endif
!!$        enddo
!!$        enddo
!!$
!!$      else if (imode ==-1) then
!!$!     calculate  snow depth from water eqivalent
!!$        do j = 1,je
!!$        do i = 1,ie
!!$          if (snw(i,j) >= 0.) then
!!$            if (rho(i,j) < 1.) rho(i,j)=rhosnw
!!$            snw(i,j) = snw(i,j)/rho(i,j)
!!$          else
!!$            rho(i,j) = 0.
!!$          endif
!!$!         scale snow depth to have it in mm 
!!$          snw(i,j) = snw(i,j)*1000.
!!$        enddo
!!$        enddo
!!$      endif
!!$
!!$      end

!==============================================================================

END MODULE sfc_terra_init
