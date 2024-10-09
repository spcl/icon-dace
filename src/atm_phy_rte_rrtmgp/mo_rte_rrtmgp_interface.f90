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

MODULE mo_rte_rrtmgp_interface
  USE mo_kind,                       ONLY: wp
  USE mo_math_constants,             ONLY: pi
  USE mo_physical_constants,         ONLY: rhoh2o, rd_o_cpd
  USE mo_exception,                  ONLY: finish
#ifdef _OPENACC
  USE mo_exception,                  ONLY: warning
#endif
  USE mo_parallel_config,            ONLY: nproma_sub
  USE mo_bc_aeropt_kinne,            ONLY: set_bc_aeropt_kinne
  USE mo_bc_aeropt_splumes,          ONLY: add_bc_aeropt_splumes

  USE mo_optical_props,              ONLY: ty_optical_props_1scl, &
                                           ty_optical_props_2str
  USE mo_gas_concentrations,         ONLY: ty_gas_concs
  USE mo_source_functions,           ONLY: ty_source_func_lw
  USE mo_rte_lw,                     ONLY: rte_lw
  USE mo_rte_sw,                     ONLY: rte_sw
  USE mo_fluxes,                     ONLY: ty_fluxes_broadband
  USE mo_icon_fluxes_sw,             ONLY: ty_icon_fluxes_sw, set_fractions
  USE mo_rte_rrtmgp_setup,           ONLY: k_dist_lw, k_dist_sw, &
                                           cloud_optics_lw, cloud_optics_sw, &
                                           stop_on_err, inhoml, inhomi, inhoms

  USE mo_rad_diag,                   ONLY: rad_aero_diag
  USE mo_timer,                      ONLY: ltimer, timer_start, timer_stop, &
   &                                       timer_rte_rrtmgp_int, &
   &                                       timer_rte_rrtmgp_int_onb, &
   &                                       timer_gas_concs, &
   &                                       timer_clamp_pr_temp, &
   &                                       timer_source_lw, &
   &                                       timer_atmos_lw, &
   &                                       timer_k_dist_lw, &
   &                                       timer_aerosol_lw, &
   &                                       timer_rte_lw_clrsky, &
   &                                       timer_clouds_bnd_lw, &
   &                                       timer_cloud_optics_lw, &
   &                                       timer_snow_bnd_lw, &
   &                                       timer_rte_lw_allsky, &
   &                                       timer_atmos_sw, &
   &                                       timer_k_dist_sw, &
   &                                       timer_aerosol_sw, &
   &                                       timer_rte_sw_clrsky, &
   &                                       timer_clouds_bnd_sw, &
   &                                       timer_cloud_optics_sw, &
   &                                       timer_snow_bnd_sw, &
   &                                       timer_rte_sw_allsky
  USE mo_radiation_general,          ONLY: wavenum1, wavenum2
  USE mo_aes_rad_config,             ONLY: aes_rad_config
  USE mtime,                         ONLY: datetime


#ifdef RRTMGP_MERGE_DEBUG
  USE mo_rte_rrtmgp_merge_debug, ONLY: write_record_interface_aes
#endif

! These need to be sent once in the init phae from the atmo to the ps_rad
!          zf(kbdim,klev),               & !< geometric height at full level in m
!          zh(kbdim,klev+1),             & !< geometric height at half level in m
!          dz(kbdim,klev),               & !< geometric height thickness in m

  IMPLICIT NONE

  PRIVATE
#ifdef _OPENACC
  LOGICAL, PARAMETER :: use_acc  = .TRUE.
#else
  LOGICAL, PARAMETER :: use_acc  = .FALSE.
#endif

  LOGICAL, PARAMETER :: top_at_1 = .true.
  LOGICAL            :: lneed_aerosols

  REAL(wp), PARAMETER :: pressure_scale = 100._wp,         &
                         inverse_pressure_scale = 1._wp/pressure_scale, &
                         droplet_scale = 1.0e2, &
                         nir_vis_boundary   = 14500._wp

  PUBLIC :: rte_rrtmgp_interface, pressure_scale, inverse_pressure_scale, &
            droplet_scale

CONTAINS

  !-------------------------------------------------------------------

  !>
  !! arranges input and calls rrtm sw and lw routines
  !!
  !! Remarks
  !!   Some cloud physical properties are prescribed, which are
  !!   required to derive cloud optical properties
  !!
  !! The gases are passed into RRTM via two multi-constituent arrays:


  ! TODO/BUG?
  !-------------------------------------------------------------------
  SUBROUTINE rte_rrtmgp_interface(                                          &
      & jg, jb, jcs, jce, nproma, klev                                     ,&
      & irad_aero       , lrad_aero_diag, lrad_coupled                     ,&
      & psctm           , ssi_factor                                       ,&
      & loland          ,loglac          ,this_datetime                    ,&
      & pcos_mu0        ,daylght_frc                                       ,&
      & alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif     ,&
      & emissivity                                                         ,&
      & zf              ,zh              ,dz                               ,&
      & pp_sfc          ,pp_fl           ,pp_hl                            ,&
      & tk_sfc          ,tk_fl           ,tk_hl                            ,&
      & rad_2d                                                             ,&
      & xvmr_vap        ,xm_liq          ,xm_ice                           ,&
      & reff_ice        ,tau_ice         ,reff_snow       ,tau_snow        ,&
      & cdnc            ,xc_frc          ,xm_snw                           ,&
      & xvmr_co2        ,xvmr_ch4        ,xvmr_n2o        ,xvmr_cfc        ,&
      & xvmr_o3         ,xvmr_o2                                           ,&
      & lw_upw          ,lw_upw_clr      ,lw_dnw          ,lw_dnw_clr      ,&
      & sw_upw          ,sw_upw_clr      ,sw_dnw          ,sw_dnw_clr      ,&
      & vis_dn_dir_sfc  ,par_dn_dir_sfc  ,nir_dn_dir_sfc                   ,&
      & vis_dn_dff_sfc  ,par_dn_dff_sfc  ,nir_dn_dff_sfc                   ,&
      & vis_up_sfc      ,par_up_sfc      ,nir_up_sfc                       ,&
      & aer_aod_533     ,aer_ssa_533     ,aer_asy_533                      ,&
      & aer_aod_2325    ,aer_ssa_2325    ,aer_asy_2325                     ,&
      & aer_aod_9731                                                       )
#ifdef __INTEL_COMPILER
!DIR$ OPTIMIZE:1
#endif
     !-------------------------------------------------------------------

    INTEGER,INTENT(IN) :: &
         jg,           & !< domain index
         jb,           & !< block index
         jcs, jce,     & !< starting and ending columns
         nproma, klev, & !< array dimensions(?)
         irad_aero       !< aerosol control

    LOGICAL, INTENT(IN) :: lrad_aero_diag                !< diagnose aerosol optical properties
    LOGICAL, INTENT(IN) :: lrad_coupled                  !< kinne aerosol from coupler (true) or file
    REAL(wp),INTENT(IN) :: psctm                         !< orbit and time dependent solar constant for radiation time step
    REAL(wp),INTENT(IN) :: ssi_factor(:)                 !< fraction of TSI in the 14 RRTM SW bands

    LOGICAL,INTENT(IN) ::              &
         loland(:),                & !< land sea mask, land=.true.
         loglac(:)                   !< glacier mask, glacier=.true.

    TYPE(datetime), POINTER ::  this_datetime !< actual time step

    REAL(WP),INTENT(IN)  :: &
         pcos_mu0(:),     & !< mu0 for solar zenith angle
         daylght_frc(:),  & !< daylight fraction; with diurnal cycle 0 or 1, with zonal mean in [0,1]
         alb_vis_dir(:),  & !< surface albedo for vis range and dir light
         alb_nir_dir(:),  & !< surface albedo for NIR range and dir light
         alb_vis_dif(:),  & !< surface albedo for vis range and dif light
         alb_nir_dif(:),  & !< surface albedo for NIR range and dif light
         emissivity(:),   & !< sufrace emissivity
         zf(:,:),         & !< geometric height at full level in m
         zh(:,:),         & !< geometric height at half level in m
         dz(:,:),         & !< geometric height thickness in m
         pp_sfc(:),       & !< surface pressure in Pa
         pp_fl(:,:),      & !< full level pressure in Pa
         pp_hl(:,:),      & !< half level pressure in Pa
         tk_sfc(:),       & !< surface temperature in K
         tk_fl(:,:),      & !< full level temperature in K
         tk_hl(:,:),      & !< half level temperature in K
         xvmr_vap(:,:),   & !< water vapor volume mixing ratio 
         xm_liq(:,:),     & !< cloud water mass in kg/m2
         xm_ice(:,:),     & !< cloud ice   mass in kg/m2
         cdnc(:,:),       & !< cloud nuclei concentration
         xc_frc(:,:),     & !< fractional cloud cover
         xm_snw(:,:),     & !< snow        mass in kg/m2
         xvmr_co2(:,:),   & !< co2 volume mixing ratio
         xvmr_ch4(:,:),   & !< ch4 volume mixing ratio
         xvmr_n2o(:,:),   & !< n2o volume mixing ratio
         xvmr_cfc(:,:,:), & !< cfc volume mixing ratio
         xvmr_o3(:,:),    & !< o3  volume mixing ratio
         xvmr_o2(:,:),    & !< o2  volume mixing ratio
         reff_ice(:,:),   & !< cloud ice effective radius in m
         reff_snow(:,:)     !< snow effective radius in m

    REAL(wp), INTENT(INOUT) :: &
         tau_ice(:,:),    & !< optical depth of cloud ice integrated over bands
         tau_snow(:,:)      !< optical depth of snow integrated over bands

    REAL(wp), INTENT(INOUT) :: &
         rad_2d(:)          !< arbitrary 2d field for output inside radiation

    REAL(wp), INTENT(OUT)   :: &
      & lw_dnw_clr(:,:),& !< Clear-sky downward longwave  at all levels
      & lw_upw_clr(:,:),& !< Clear-sky upward   longwave  at all levels
      & sw_dnw_clr(:,:),& !< Clear-sky downward shortwave at all levels
      & sw_upw_clr(:,:),& !< Clear-sky upward   shortwave at all levels
      & lw_dnw(:,:),    & !< All-sky   downward longwave  at all levels
      & lw_upw(:,:),    & !< All-sky   upward   longwave  at all levels
      & sw_dnw(:,:),    & !< All-sky   downward shortwave at all levels
      & sw_upw(:,:)       !< All-sky   upward   shortwave at all levels

    REAL(wp), INTENT(OUT) ::  &
         vis_dn_dir_sfc(:),   & !< Diffuse downward flux surface visible radiation
         par_dn_dir_sfc(:),   & !< Diffuse downward flux surface PAR
         nir_dn_dir_sfc(:),   & !< Diffuse downward flux surface near-infrared radiation
         vis_dn_dff_sfc(:),   & !< Direct  downward flux surface visible radiation
         par_dn_dff_sfc(:),   & !< Direct  downward flux surface PAR
         nir_dn_dff_sfc(:),   & !< Direct  downward flux surface near-infrared radiation
         vis_up_sfc    (:),   & !< Upward  flux surface visible radiation
         par_up_sfc    (:),   & !< Upward  flux surface PAR
         nir_up_sfc    (:),   & !< Upward  flux surface near-infrared radiation
         aer_aod_533   (:,:), & !< Aerosol optical density at 533 nm
         aer_ssa_533   (:,:), & !< Single scattering albedo at 533 nm
         aer_asy_533   (:,:), & !< Asymmetry factor at 533 nm
         aer_aod_2325  (:,:), & !< Aerosol optical density at 2325 nm
         aer_ssa_2325  (:,:), & !< Single scattering albedo at 2325 nm
         aer_asy_2325  (:,:), & !< Asymmetry factor at 2325 nm
         aer_aod_9731  (:,:)    !< Aerosol optical density at 9731 nm

    LOGICAL :: lclrsky_lw, lclrsky_sw
    LOGICAL :: inhom_lts
    REAL(wp) :: inhom_lts_max

    ! --------------------------------------------------------------------------
    INTEGER :: ncol_supplied, ncol_needed, jchunk_start, jchunk_end
    INTEGER :: nbndsw, nbndlw
    ! --- Aerosol optical properites - vertically reversed fields
    REAL(wp), ALLOCATABLE :: &
         aer_tau_lw(:,:,:),  & !< LW optical thickness of aerosols
         aer_tau_sw(:,:,:),  & !< aerosol optical thickness
         aer_ssa_sw(:,:,:),  & !< aerosol single scattering albedo
         aer_asy_sw(:,:,:)     !< aerosol asymmetry factor

    IF (ltimer) CALL timer_start(timer_rte_rrtmgp_int)

    ! --------------------------------------------------------------------------
    !
    ! Aerosol optical properties are computed at this level because they require
    !   geographic and temporal information; the geographic information is lost
    !   when the data provided to RTE+RRTMGP are extracted from larger arrarys
    !
    ! set all aerosols to zero first

    lneed_aerosols = (irad_aero /= 0)
    IF(lneed_aerosols) THEN
      nbndlw = k_dist_lw%get_nband()
      nbndsw = k_dist_sw%get_nband()
      
      ALLOCATE( aer_tau_lw(nproma,klev,nbndlw), &
                aer_tau_sw(nproma,klev,nbndsw), &
                aer_ssa_sw(nproma,klev,nbndsw), &
                aer_asy_sw(nproma,klev,nbndsw)  )
         
      !$ACC ENTER DATA CREATE(aer_tau_lw, aer_tau_sw, aer_ssa_sw, aer_asy_sw)

      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
      aer_tau_lw(:,:,:) = 0.0_wp
      aer_tau_sw(:,:,:) = 0.0_wp
      aer_ssa_sw(:,:,:) = 1.0_wp
      aer_asy_sw(:,:,:) = 0.0_wp
      !$ACC END KERNELS
      IF (irad_aero==12 .OR. irad_aero==13 .OR. irad_aero==19) THEN
      ! irad_aero=12 Kinne aerosols (natural background, data are read
      !      from a file without year in its name.
      ! irad_aero=13: only Kinne aerosols are used
      ! irad_aero=19: Kinne aerosols (background of natural origin,
      ! read from a file without year in its name!) + simple plumes
        CALL set_bc_aeropt_kinne(this_datetime,                        &
              & jg,                                                    &
              & jcs,            jce,                   nproma,         &
              & klev,           jb,                                    &
              & nbndsw,         nbndlw,                                &
              & zf,             dz,                                    &
              & aer_tau_sw,     aer_ssa_sw,            aer_asy_sw,     &
              & aer_tau_lw,                                            &
              & opt_from_coupler=lrad_coupled, opt_use_acc=use_acc     )
      END IF
      IF (irad_aero==19) THEN
      ! Simple plumes are added to ...
      ! iaero=19: ... Kinne background aerosols (of natural origin, 1850)
#ifdef _OPENACC
        CALL warning('mo_rte_rrtmgp_interface/rte_rrtmgp_interface','Plumes ACC not implemented')
#endif
        !$ACC UPDATE HOST(aer_tau_lw, aer_tau_sw, aer_ssa_sw, aer_asy_sw, zf, dz, zh(:,klev+1)) ASYNC(1)
        !$ACC WAIT(1)
        CALL add_bc_aeropt_splumes(                                      &
              & jg,          jcs,         jce,           nproma,         & 
              & klev,        jb,          nbndsw,        this_datetime,  &
              & zf,          dz,          zh(:,klev+1),  wavenum1,       &
              & wavenum2,    aer_tau_sw,  aer_ssa_sw,    aer_asy_sw      )
        !$ACC UPDATE DEVICE(aer_tau_sw, aer_ssa_sw, aer_asy_sw) ASYNC(1)
      END IF

      ! this should be decativated in the concurrent version and make the aer_* global variables for output
      IF (lrad_aero_diag) THEN
        CALL rad_aero_diag (                                  &
          & jcs,             jce,             nproma,         &
          & klev,            nbndlw,          nbndsw,         &
          & aer_tau_lw,      aer_tau_sw,      aer_ssa_sw,     &
          & aer_asy_sw,                                       &
          & aer_aod_533,     aer_ssa_533,     aer_asy_533,    &
          & aer_aod_2325,    aer_ssa_2325,    aer_asy_2325,   &
          & aer_aod_9731,                                     &
          & opt_use_acc = use_acc                             )
      ENDIF
      !
      ! Map solar bands from RRTMG to GP order
      !
      CALL rearrange_bands2rrtmgp(nproma,klev,nbndsw, aer_tau_sw)
      CALL rearrange_bands2rrtmgp(nproma,klev,nbndsw, aer_ssa_sw)
      CALL rearrange_bands2rrtmgp(nproma,klev,nbndsw, aer_asy_sw)

      ! Aerosol optical properties have reverse vertical orientation - reorient
      !DA TODO
      CALL reorient_3d_wrt2(aer_tau_lw)
      CALL reorient_3d_wrt2(aer_tau_sw)
      CALL reorient_3d_wrt2(aer_ssa_sw)
      CALL reorient_3d_wrt2(aer_asy_sw)
    ELSE
      ! allocate dummy zero-size arrays
      ALLOCATE( aer_tau_lw(1,1,0), &
                aer_tau_sw(1,1,0), &
                aer_ssa_sw(1,1,0), &
                aer_asy_sw(1,1,0)  )
    END IF
    !
    ! Turn geography and number concentrations into effective radii
    !
    ! --------------------------------------------------------------------------
    ! Set flag for the optional computation of clear-sky fluxes
    lclrsky_lw    = aes_rad_config(jg)%lclrsky_lw
    lclrsky_sw    = aes_rad_config(jg)%lclrsky_sw
    !
    inhom_lts     = aes_rad_config(jg)%inhom_lts
    inhom_lts_max = aes_rad_config(jg)%inhom_lts_max
    ! --------------------------------------------------------------------------
    !
    !
    ! Assumption: all arrays share the "column" dimension
    !
    ncol_supplied = size(pcos_mu0) ! baustelle - this should be = nproma?
    ncol_needed   = jce-jcs+1
    !
    ! RTE+RRTMGP process all columns supplied and assume a starting index of 1.
    !   If these conditions are satisfied we can call the interface directly...
    !
    IF (jcs==1 .and. ncol_needed == ncol_supplied .and. nproma_sub == ncol_needed) THEN

       CALL rte_rrtmgp_interface_onBlock(                              &
          & lclrsky_lw,        lclrsky_sw,                             &
          & inhom_lts,         inhom_lts_max,                          &
          & ncol_needed,       klev,                                   &
          & psctm,             ssi_factor,                             &
          & loland(:),         loglac(:),                              &
          & pcos_mu0(:),       daylght_frc(:),                         &
          & alb_vis_dir(:),    alb_nir_dir(:),                         &
          & alb_vis_dif(:),    alb_nir_dif(:),                         &
          & emissivity(:),                                             &
          & zf(:,:),           zh(:,:),           dz(:,:),             &
          & pp_sfc(:),         pp_fl(:,:),        pp_hl(:,:),          &
          & tk_sfc(:),         tk_fl(:,:),        tk_hl(:,:),          &
          & rad_2d(:),                                                 &
          & xvmr_vap(:,:),     xm_liq(:,:),                            &
          & xm_ice(:,:),       reff_ice(:,:),     tau_ice(:,:),        &
          & reff_snow(:,:),    tau_snow(:,:),                          &
          & cdnc(:,:),         xc_frc(:,:),       xm_snw(:,:),         &
          & xvmr_co2(:,:),     xvmr_ch4(:,:),     xvmr_n2o (:,:),      &
          & xvmr_cfc(:,:,:),   xvmr_o3(:,:),      xvmr_o2(:,:),        &
          & aer_tau_lw(:,:,:),                                         &
          & aer_tau_sw(:,:,:), aer_ssa_sw(:,:,:), aer_asy_sw(:,:,:),   &
          !
          & lw_upw(:,:),       lw_upw_clr (:,:),                       &
          & lw_dnw(:,:),       lw_dnw_clr (:,:),                       &
          & sw_upw(:,:),       sw_upw_clr(:,:),                        &
          & sw_dnw(:,:),       sw_dnw_clr(:,:),                        &
          & vis_dn_dir_sfc(:), par_dn_dir_sfc(:), nir_dn_dir_sfc(:),   &
          & vis_dn_dff_sfc(:), par_dn_dff_sfc(:), nir_dn_dff_sfc(:),   &
          & vis_up_sfc(:),     par_up_sfc(:),     nir_up_sfc(:)        )
       !
    ELSE
       !
       ! ... but if the starting column is > 1 and/or there are trailing columns
       !   which shouldn't be processed we make copies to supply contiguous memory
       !   to RTE+RRMTPG
       !

       DO jchunk_start = jcs,jce, nproma_sub
        jchunk_end = MIN(jchunk_start + nproma_sub - 1, jce)
        CALL shift_and_call_rte_rrtmgp_interface_onBlock(                &
            & lclrsky_lw,        lclrsky_sw,                             &
            & inhom_lts,         inhom_lts_max,                          &
            & jchunk_start,      jchunk_end,                             &
            & klev,                                                      &
            & psctm,             ssi_factor,                             &
            & loland(:),         loglac(:),                              &
            & pcos_mu0(:),       daylght_frc(:),                         &
            & alb_vis_dir(:),    alb_nir_dir(:),                         &
            & alb_vis_dif(:),    alb_nir_dif(:),                         &
            & emissivity(:),                                             &
            & zf(:,:),           zh(:,:),           dz(:,:),             &
            & pp_sfc(:),         pp_fl(:,:),        pp_hl(:,:),          &
            & tk_sfc(:),         tk_fl(:,:),        tk_hl(:,:),          &
            & rad_2d(:),                                                 &
            & xvmr_vap(:,:),     xm_liq(:,:),                            &
            & xm_ice(:,:),       reff_ice(:,:),     tau_ice(:,:),        &
            & reff_snow(:,:),    tau_snow(:,:),                          &
            & cdnc(:,:),         xc_frc(:,:),       xm_snw(:,:),         &
            & xvmr_co2(:,:),     xvmr_ch4(:,:),     xvmr_n2o (:,:),      &
            & xvmr_cfc(:,:,:),   xvmr_o3(:,:),      xvmr_o2(:,:),        &
            & aer_tau_lw(:,:,:),                                         &
            & aer_tau_sw(:,:,:), aer_ssa_sw(:,:,:), aer_asy_sw(:,:,:),   &
            !
            & lw_upw(:,:),       lw_upw_clr (:,:),                       &
            & lw_dnw(:,:),       lw_dnw_clr (:,:),                       &
            & sw_upw(:,:),       sw_upw_clr(:,:),                        &
            & sw_dnw(:,:),       sw_dnw_clr(:,:),                        &
            & vis_dn_dir_sfc(:), par_dn_dir_sfc(:), nir_dn_dir_sfc(:),   &
            & vis_dn_dff_sfc(:), par_dn_dff_sfc(:), nir_dn_dff_sfc(:),   &
            & vis_up_sfc(:),     par_up_sfc(:),     nir_up_sfc(:)        )
       END DO
       !
    END IF

  !$ACC WAIT
  !$ACC EXIT DATA DELETE(aer_tau_lw, aer_tau_sw, aer_ssa_sw, aer_asy_sw) IF(lneed_aerosols)

    IF (ltimer) CALL timer_stop(timer_rte_rrtmgp_int)

  END SUBROUTINE rte_rrtmgp_interface
 ! -------------------------------------------------------------------------------------
  SUBROUTINE clamp_pressure(src, tgt, low, high)
    REAL(wp), INTENT(IN) :: src(:,:)
    REAL(wp), INTENT(OUT) :: tgt(:,:)
    REAL(wp), INTENT(IN) :: low, high

    !tgt(:,:) = min(high, max(low, src))
    INTEGER :: i, j, m, n

    ! min and max are level-dependent
    REAL(wp), DIMENSION(SIZE(src,2)) :: tgt_min, tgt_max 
    
    !$ACC DATA CREATE(tgt_min, tgt_max) PRESENT(src, tgt)
    
    m = SIZE(src,1)
    n = SIZE(src,2)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO j = 1, n
      tgt_min(j)=low  + (j-1)*epsilon(tgt_min)
      tgt_max(j)=high - (n-j)*epsilon(tgt_max)
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO j = 1, n
      DO i = 1 , m
         tgt(i,j) = min(tgt_max(j), max(tgt_min(j), src(i,j)))
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC WAIT(1)
    !$ACC END DATA
  END SUBROUTINE clamp_pressure

  SUBROUTINE clamp_temperature(src, tgt, low, high)
    REAL(wp), INTENT(IN) :: src(:,:)
    REAL(wp), INTENT(OUT) :: tgt(:,:)
    REAL(wp), INTENT(IN) :: low, high

    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
    tgt(:,:) = min(high, max(low, src(:,:)))
    !$ACC END KERNELS
  END SUBROUTINE clamp_temperature
  !----------------------------------------------- !>
  !! arranges input and calls rrtm sw and lw routines
  !!
  !! Remarks
  !!   Because the RRTM indexes vertical levels differently than ECHAM a chief
  !!   function of thise routine is to reorder the input in the vertical.  In
  !!   addition some cloud physical properties are prescribed, which are
  !!   required to derive cloud optical properties
  !!

  SUBROUTINE rte_rrtmgp_interface_onBlock(                   &
       & lclrsky_lw,     lclrsky_sw,                         &
       & inhom_lts,      inhom_lts_max,                      &
       & ncol,           klev,                               &
       & psctm,          ssi_factor,                         &
       & laland,         laglac,                             &
       & pcos_mu0,       daylght_frc,                        &
       & alb_vis_dir,    alb_nir_dir,                        &
       & alb_vis_dif,    alb_nir_dif,                        &
       & emissivity,                                         &
       & zf,             zh,             dz,                 &
       & pp_sfc,         pp_fl,          pp_hl,              &
       & tk_sfc,         tk_fl,          tk_hl,              &
       & rad_2d,                                             &
       & xvmr_vap,       xm_liq,                             &
       & xm_ice,         reff_ice,       tau_ice,            &
       & reff_snow,      tau_snow,                           &
       & cdnc,           cld_frc,        xm_snw,             &
       & xvmr_co2,       xvmr_ch4,       xvmr_n2o ,          &
       & xvmr_cfc ,      xvmr_o3,        xvmr_o2,            &
       & aer_tau_lw,                                         &
       & aer_tau_sw,     aer_ssa_sw  ,   aer_asy_sw,         &
       & flx_uplw,       flx_uplw_clr,                       &
       & flx_dnlw,       flx_dnlw_clr,                       &
       & flx_upsw,       flx_upsw_clr,                       &
       & flx_dnsw,       flx_dnsw_clr,                       &
       & vis_dn_dir_sfc, par_dn_dir_sfc, nir_dn_dir_sfc,     &
       & vis_dn_dff_sfc, par_dn_dff_sfc, nir_dn_dff_sfc,     &
       & vis_up_sfc,     par_up_sfc,     nir_up_sfc          )

#ifdef __INTEL_COMPILER
!DIR$ OPTIMIZE:1
#endif

    LOGICAL,INTENT(IN)  :: lclrsky_lw                    !< flag for LW clear-sky computations
    LOGICAL,INTENT(IN)  :: lclrsky_sw                    !< flag for SW clear-sky computations
    LOGICAL,INTENT(IN)  :: inhom_lts
    REAL(wp),INTENT(IN) :: inhom_lts_max                 !< maximum value on inhoml

    INTEGER,INTENT(IN)  :: &
         ncol,             & !< number of columns
         klev                !< number of levels

    REAL(wp),INTENT(IN) :: psctm                         !< orbit and time dependent solar constant for radiation time step
    REAL(wp),INTENT(IN) :: ssi_factor(:)                 !< fraction of TSI in the 14 RRTM SW bands

    LOGICAL,INTENT(IN) :: &
         laland(:),   & !< land sea mask, land=.true.
         laglac(:)      !< glacier mask, glacier=.true.

    REAL(WP),INTENT(IN)  ::    &
         pcos_mu0(:),      & !< mu0 for solar zenith angle
         daylght_frc(:),   & !< daylight fraction; with diurnal cycle 0 or 1, with zonal mean in [0,1]
         alb_vis_dir(:),   & !< surface albedo for vis range and dir light
         alb_nir_dir(:),   & !< surface albedo for NIR range and dir light
         alb_vis_dif(:),   & !< surface albedo for vis range and dif light
         alb_nir_dif(:),   & !< surface albedo for NIR range and dif light
         emissivity(:),    & !< surface longwave emissivity
         zf(:,:),          & !< geometric height at full level in m
         zh(:,:),          & !< geometric height at half level in m
         dz(:,:),          & !< geometric height thickness in m
         
         pp_sfc(:),        & !< surface pressure in Pa
         pp_fl(:,:),       & !< full level pressure in Pa
         pp_hl(:,:),       & !< full level pressure in Pa
         tk_sfc(:),        & !< surface temperature in K
         tk_fl(:,:),       & !< full level temperature in K
         tk_hl(:,:),       & !< half level temperature in K
         xvmr_vap(:,:),    & !< water vapor volume mixing ratio
         xm_liq(:,:),      & !< cloud water mass in kg/m2
         xm_ice(:,:),      & !< cloud ice   mass in kg/m2
         xm_snw(:,:),      & !< snow        mass in kg/m2
         aer_tau_lw(:,:,:),& !< aerosol optical depth, longwave (ncol, nlay, nbndlw)
         aer_tau_sw(:,:,:),& !< aerosol optical depth,            shortwave (ncol, nlay, nbndlw)
         aer_ssa_sw(:,:,:),& !< aerosol single-scattering albedo, shortwave (ncol, nlay, nbndlw)
         aer_asy_sw(:,:,:),& !< aerosol asymetry parameter,       shortwave (ncol, nlay, nbndlw)
         cdnc(:,:),        & !< cloud nuclei concentration
         cld_frc(:,:),     & !< fractional cloud cover
         xvmr_co2(:,:),    & !< co2 volume mixing ratio
         xvmr_ch4(:,:),    & !< ch4 volume mixing ratio
         xvmr_n2o(:,:),    & !< n2o volume mixing ratio
         xvmr_cfc(:,:,:),  & !< cfc volume mixing ratio (kbdim,klev,2)
         xvmr_o3(:,:),     & !< o3  volume mixing ratio
         xvmr_o2(:,:),     & !< o2  volume mixing ratio
         reff_ice(:,:),    & !< cloud ice effective radius m
         reff_snow(:,:)      !< snow effectiv radius m

    REAL (wp), INTENT (INOUT) :: &
         tau_ice(:,:),     & !< optical depth of cloud ice integrated over bands
         tau_snow(:,:)       !< optical depth of snow integrated over bands

    REAL (wp), INTENT (INOUT) :: &
         rad_2d(:)           !< arbitrary 2d-field in radiation for output
         

    REAL (wp), TARGET, INTENT (INOUT) ::       &
         flx_uplw    (:,:), & !<   upward LW flux profile, all sky
         flx_uplw_clr(:,:), & !<   upward LW flux profile, clear sky
         flx_dnlw    (:,:), & !< downward LW flux profile, all sky
         flx_dnlw_clr(:,:), & !< downward LW flux profile, clear sky
         flx_upsw    (:,:), & !<   upward SW flux profile, all sky
         flx_upsw_clr(:,:), & !<   upward SW flux profile, clear sky
         flx_dnsw    (:,:), & !< downward SW flux profile, all sky
         flx_dnsw_clr(:,:)    !< downward SW flux profile, clear sky

    REAL (wp), TARGET, INTENT (INOUT) :: &
         vis_dn_dir_sfc(:) , & !< Diffuse downward flux surface visible radiation
         par_dn_dir_sfc(:) , & !< Diffuse downward flux surface PAR
         nir_dn_dir_sfc(:) , & !< Diffuse downward flux surface near-infrared radiation
         vis_dn_dff_sfc(:) , & !< Direct  downward flux surface visible radiation
         par_dn_dff_sfc(:) , & !< Direct  downward flux surface PAR
         nir_dn_dff_sfc(:) , & !< Direct  downward flux surface near-infrared radiation
         vis_up_sfc    (:) , & !< Upward  flux surface visible radiation
         par_up_sfc    (:) , & !< Upward  flux surface PAR
         nir_up_sfc    (:)     !< Upward  flux surface near-infrared radiation

    ! -----------------------------------------------------------------------
    ! At this stage, columns are stored in the arrays in a contiguous form, i.e.
    ! all loops run over (1:ncol, 1:klev) where ncol=jce-jcs+1 and arrays are copied
    ! into the arrays before accordingly.

    INTEGER  :: jk, jl !< loop indices
    INTEGER  :: nbndlw, nbndsw !, ngptsw
    REAL(wp) ::                      &
         zsemiss(k_dist_lw%get_nband(),ncol) !< LW surface emissivity by band
    ! --- local scaled variables
    REAL(wp) ::                   &
         cld_frc_loc, & !< secure cloud fraction
         ziwp       (ncol,klev), & !< in cloud ice water content       [g/m2]
         zlwp       (ncol,klev), & !< in cloud liquid water content    [g/m2]
         zlwc                      !< in cloud water concentration     [g/m3]

    REAL(wp) ::                &
         re_drop (ncol,klev), & !< effective radius of liquid
         re_cryst(ncol,klev)
    REAL(wp) ::                  &
         zswp       (ncol,klev), & !< snow water path [g/m2]
         zdwp       (ncol,klev), & !< dummy water path 
         re_snow(ncol,klev)        !< snow effective radius
    !
    ! Random seeds for sampling. Needs to get somewhere upstream
    !
    INTEGER :: band, i, j
    REAL(wp) :: low, high

    TYPE(ty_source_func_lw)     :: source_lw !check types regarding acc later
    TYPE(ty_optical_props_1scl) :: atmos_lw !check types regarding acc later
    TYPE(ty_optical_props_1scl) :: aerosol_lw !check types regarding acc later
    TYPE(ty_optical_props_1scl) :: clouds_bnd_lw !check types regarding acc later
    TYPE(ty_optical_props_1scl) :: snow_bnd_lw ! for snow optics
    TYPE(ty_optical_props_2str) :: atmos_sw !check types regarding acc later
    TYPE(ty_optical_props_2str) :: aerosol_sw !check types regarding acc later
    TYPE(ty_optical_props_2str) :: clouds_bnd_sw !check types regarding acc later
    TYPE(ty_optical_props_2str) :: snow_bnd_sw

    TYPE(ty_fluxes_broadband) :: fluxes_lw !check acc
    TYPE(ty_icon_fluxes_sw  ) :: fluxes_sw !check acc
    TYPE(ty_fluxes_broadband) :: fluxes_lwcs, fluxes_swcs !check acc

    TYPE(ty_gas_concs) :: gas_concs !check acc
    REAL(wp), DIMENSION(ncol       ) :: mu0
    REAL(wp), DIMENSION(ncol,klev  ) :: tlay, play
    REAL(wp), DIMENSION(ncol,klev+1) :: tlev, plev
    REAL(wp), DIMENSION(k_dist_sw%get_nband(),ncol) :: albdir, albdif
    REAL(wp) :: band_lims(2,k_dist_sw%get_nband()), delwave, frc_vis
    REAL(wp) :: toa_flux(ncol,k_dist_sw%get_ngpt())
    ! Threshold within which a cloud fraction is considered = 0 or 1.
    REAL(wp) :: cld_frc_thresh
    !--------------------------------
    INTEGER, PARAMETER :: n_gas_names = 8
    CHARACTER(len=5), PARAMETER :: gas_names(n_gas_names) = (/ &
       'h2o  ', 'co2  ', 'ch4  ', 'o2   ', 'o3   ', 'n2o  ','cfc11', 'cfc12'/)
    !--------------------------------
    ! Variables for effective radii computations
    REAL (wp), PARAMETER :: &
       ccwmin = 1.e-7_wp, &    ! min condensate for lw cloud opacity
       zkap_cont = 1.143_wp, & ! continental (Martin et al. ) breadth param
       zkap_mrtm = 1.077_wp, & ! maritime (Martin et al.) breadth parameter
       del1      = 2._wp,    & ! transition factor for inhomogeneity stability scaling
       del2      = 20._wp      ! cut-overpoint for inhomogeneity stability scaling
    REAL (wp) :: effective_radius
    REAL (wp) :: reimin, reimax, relmin, relmax, zkap
    REAL (wp) :: lts
    LOGICAL   :: lcldlyr
    !
    !DA TODO: rearrange the data section to reduce memory consumption
    !
    !$ACC DATA PRESENT(cld_frc, xm_ice, xm_liq, dz, pcos_mu0, emissivity) &
    !$ACC   PRESENT(alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif) &
    !$ACC   PRESENT(daylght_frc, laland, laglac, dz, cdnc, xm_snw) &
    !$ACC   PRESENT(reff_ice, tau_ice, reff_snow, tau_snow) &
    !$ACC   PRESENT(tk_sfc, pp_sfc, tk_fl, pp_fl) &
    !$ACC   CREATE(ziwp, zlwp, mu0, zsemiss, albdif, re_cryst, re_drop) &
    !$ACC   CREATE(zswp, zdwp, re_snow) &
    !$ACC   CREATE(albdir, toa_flux) &
    !$ACC   CREATE(plev, play, tlev, tlay)

    IF (ltimer) CALL timer_start(timer_rte_rrtmgp_int_onb)

    nbndlw = k_dist_lw%get_nband()
    nbndsw = k_dist_sw%get_nband()
!    ngptsw = k_dist_sw%get_ngpt()
    cld_frc_thresh = 4._wp*spacing(1._wp)

    ! 1.0 Constituent properties
    !--------------------------------
    !
    ! 1.1 Cloud condensate and cloud fraction
    !
    effective_radius = &
      1.0e6_wp * droplet_scale * (3.0e-9_wp / (4.0_wp * pi * rhoh2o))**(1.0_wp/3.0_wp) 

    reimin = MAX(cloud_optics_lw%get_min_radius_ice(), cloud_optics_sw%get_min_radius_ice()) ! 10.0_wp  !
    reimax = MIN(cloud_optics_lw%get_max_radius_ice(), cloud_optics_sw%get_max_radius_ice()) ! 124.0_wp !
    
    relmin = MAX(cloud_optics_lw%get_min_radius_liq(), cloud_optics_sw%get_min_radius_liq()) ! 2.5_wp  ! 
    relmax = MIN(cloud_optics_lw%get_max_radius_liq(), cloud_optics_sw%get_max_radius_liq()) ! 21.5_wp ! 

    IF (relmax <= relmin .OR. reimax <= reimin) THEN
      CALL finish('rte_rrtmgp_interface_onBlock (mo_rte_rrtmgp_interface.f90)', &
                  'Droplet minimun size required is bigger than maximum')
    END IF

    IF (inhom_lts) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO jl = 1, ncol
         lts = tk_fl(jl,min(73,klev))*(1e5_wp/pp_fl(jl,min(73,klev)))**(rd_o_cpd) - tk_sfc(jl)*(1e5_wp/pp_sfc(jl))**(rd_o_cpd)
         rad_2d(jl) = inhoml + (inhom_lts_max-inhoml)*(1._wp - atan2(del1,(lts - del2))/pi)
      END DO 
     !$ACC END PARALLEL LOOP
     ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO jl = 1, ncol
         rad_2d(jl) = inhoml
      END DO
      !$ACC END PARALLEL LOOP
    END IF
    !
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1, klev
      DO jl = 1, ncol
        !
        ! --- Cloud liquid and ice mass: [kg/m2 in cell] --> [g/m2 in cloud]
        !
        cld_frc_loc = MAX(EPSILON(1.0_wp),cld_frc(jl,jk))
        ziwp(jl,jk) = xm_ice(jl,jk)*1000.0_wp/cld_frc_loc
        zlwp(jl,jk) = xm_liq(jl,jk)*1000.0_wp/cld_frc_loc
        zswp(jl,jk) = xm_snw(jl,jk)*1000.0_wp

        !
        ! --- Initialize variables for snow optics
        !
        tau_ice(jl,jk) = 0._wp
        tau_snow(jl,jk) = 0._wp
        zdwp(jl,jk) = 0._wp

        ! Mask which tells cloud optics that this cell is clear
        lcldlyr = cld_frc(jl,jk) > cld_frc_thresh !!!
        IF (.NOT. lcldlyr) THEN
          ziwp(jl,jk) = 0.0_wp
          zlwp(jl,jk) = 0.0_wp
        END IF
        !
        ! --- cloud water concentration [g/m3]
        !
        zlwc = zlwp(jl,jk)/dz(jl,jk)
        !
        IF (lcldlyr .AND. (zlwp(jl,jk)+ziwp(jl,jk))>ccwmin) THEN

          zkap = zkap_mrtm
          IF ( laland(jl) .AND. .NOT.laglac(jl) ) zkap = zkap_cont
          re_cryst(jl,jk) = MAX(reimin, MIN(reimax, 1.e6_wp * reff_ice(jl,jk)))
          re_drop(jl,jk)  = MAX(relmin, MIN(relmax, &
            effective_radius * zkap * (zlwc / cdnc(jl,jk))**(1.0_wp/3.0_wp) ))
        ELSE
          re_cryst(jl,jk) = reimin
          re_drop (jl,jk) = relmin
        END IF

        ! we take the same minimum condition for "snow water path" as for cloud condensate
        ! Since we take the same interpolation tables as for cloud ice, we need the same bounds reimin, reimax
        IF (zswp(jl,jk)>ccwmin) THEN
          re_snow(jl,jk) = MAX(reimin, MIN(reimax, 1.e6_wp * reff_snow(jl,jk)))
        ELSE
          re_snow(jl,jk) = reimin
        ENDIF
      END DO
    END DO
    !$ACC END PARALLEL
    !
    ! RRTMGP cloud optics: here compute effective radius of liquid and ice from formulae in
    !   mo_radiation_cloud_optics.
    !
    !--------------------------------
    !
    ! 1.2 Gas concentrations
    !

    ! At this stage, columns are stored in the arrays in a contiguous form, i.e.
    ! all loops run over (1:ncol, 1:klev) where ncol=jce-jcs+1 and arrays are copied
    ! into the arrays before accordingly.
    ! The gas profile routine provides all gas concentrations in volume mixing ratios
    !
    ! RTE-RRTMGP ACC code is synchronous, so need to wait before calling it
    !$ACC WAIT
    !
    IF (ltimer) CALL timer_start(timer_gas_concs)
    !
    CALL stop_on_err(gas_concs%init(gas_names))
    CALL stop_on_err(gas_concs%set_vmr('h2o',   xvmr_vap))
    CALL stop_on_err(gas_concs%set_vmr('co2',   xvmr_co2))
    CALL stop_on_err(gas_concs%set_vmr('ch4',   xvmr_ch4))
    CALL stop_on_err(gas_concs%set_vmr('o2',    xvmr_o2))
    CALL stop_on_err(gas_concs%set_vmr('o3',    xvmr_o3))
    CALL stop_on_err(gas_concs%set_vmr('n2o',   xvmr_n2o))
    CALL stop_on_err(gas_concs%set_vmr('cfc11', xvmr_cfc(:,:,1)))
    CALL stop_on_err(gas_concs%set_vmr('cfc12', xvmr_cfc(:,:,2)))
    !
    IF (ltimer) CALL timer_stop (timer_gas_concs)

    !--------------------------------
    !
    ! Restrict out-of-bounds temperatures and pressures
    !
    ! The air pressure on levels plev, on the upper and lower boundaries of a layer,
    ! is used to determine the air mass transferred by radiation. For safety
    ! reasons plev is limited to values >= 0 Pa (and <=10**6 Pa so that high is defined).
    !
    IF (ltimer) CALL timer_start(timer_clamp_pr_temp)
    !
    low = 0._wp
    high = 1000000._wp
    CALL clamp_pressure(pp_hl, plev, low, high)
    !
    ! The pressure in the layer play is used for optical properties and
    ! is limited here to the range for which tables are defined in the
    ! RRTMGP data file and stored in k_dist_lw.
    ! Thus radiation can be computed for pressure values lower or higher
    ! than the defined range, albeit at lower precision.

    low =  k_dist_lw%get_press_min()
    high = k_dist_lw%get_press_max()
    CALL clamp_pressure(pp_fl, play, low, high)

    low =  k_dist_lw%get_temp_min()
    high = k_dist_lw%get_temp_max()
    CALL clamp_temperature(tk_hl, tlev, low, high)
    CALL clamp_temperature(tk_fl, tlay, low, high)
    !
    IF (ltimer) CALL timer_stop (timer_clamp_pr_temp)

    !--------------------------------
    !
    ! Boundary conditions
    !
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
    mu0(:) = MAX(1.e-10_wp,MIN(1.0_wp,pcos_mu0(:)))
    !$ACC END KERNELS

    ! 2.0 Surface Properties
    ! --------------------------------
   !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
   !$ACC LOOP GANG VECTOR COLLAPSE(2)
   DO j=1,ncol
      DO i=1,nbndlw
        zsemiss(i,j)=emissivity(j)
      END DO
    END DO
    !$ACC END PARALLEL

    !
    ! Surface albedo interpolation
    !
    !DA TODO: the next function has to run on GPU
    band_lims = k_dist_sw%get_band_lims_wavenumber()

    !$ACC PARALLEL DEFAULT(PRESENT) COPYIN(band_lims) ASYNC(1)
    !$ACC LOOP COLLAPSE(2)
    DO j=1,ncol
      DO band=1,nbndsw
        delwave = band_lims(2,band) - band_lims(1,band)
        frc_vis = MAX(0.0_wp, MIN(1.0_wp, &
          (band_lims(2,band) - nir_vis_boundary) / delwave))

        albdif(band,j) = alb_vis_dif(j) * frc_vis + &
                         alb_nir_dif(j) * (1.0_wp - frc_vis)
        albdir(band,j) = alb_vis_dir(j) * frc_vis + &
                         alb_nir_dir(j) * (1.0_wp - frc_vis)
      END DO
    END DO
    !$ACC END PARALLEL

    !
    ! 3.0 Particulate Optical Properties
    ! --------------------------------
    !
    ! 3.1 Aerosols - moved to a layer above
    !


    !
    ! 4.0 Radiative Transfer Routines
    ! --------------------------------
    !
    ! 4.1 Longwave radiative Transfer

    !
    ! 4.1.2 Gas optics
    !
    ! RTE-RRTMGP ACC code is synchronous, so need to wait before calling it
    !$ACC WAIT
    IF (ltimer) CALL timer_start(timer_source_lw)
    CALL stop_on_err(source_lw%alloc    (ncol, klev, k_dist_lw))
    IF (ltimer) CALL timer_stop (timer_source_lw)
    IF (ltimer) CALL timer_start(timer_atmos_lw)
    CALL stop_on_err(atmos_lw%alloc_1scl(ncol, klev, k_dist_lw))
    IF (ltimer) CALL timer_stop (timer_atmos_lw)
    !$ACC DATA CREATE(source_lw, atmos_lw)

    !$ACC DATA CREATE(source_lw%lay_source, source_lw%lev_source_inc) &
    !$ACC   CREATE(source_lw%lev_source_dec, source_lw%sfc_source) &
    !$ACC   CREATE(source_lw%sfc_source_Jac, atmos_lw%tau)

    IF (ltimer) CALL timer_start(timer_k_dist_lw)
    CALL stop_on_err( &
           k_dist_lw%gas_optics(play, plev, tlay, tk_sfc, &
                                gas_concs, atmos_lw, source_lw, &
                                tlev = tlev))
    IF (ltimer) CALL timer_stop (timer_k_dist_lw)
    !
    ! 4.1.2 Aerosol optical depth: add to clear-sky
    !  If irad_aero == 0, aer_tau_lw will not be allocated here
    !  and we need to skip this step
    !
    IF ( lneed_aerosols ) THEN
      IF (ltimer) CALL timer_start(timer_aerosol_lw)
      CALL stop_on_err(aerosol_lw%alloc_1scl(ncol, klev, &
                                             k_dist_lw%get_band_lims_wavenumber()))
      IF (ltimer) CALL timer_stop (timer_aerosol_lw)
      !$ACC DATA PRESENT(aer_tau_lw) CREATE(aerosol_lw)
      !$ACC DATA CREATE(aerosol_lw%tau)
      !
      !DA TODO: this can be just a pointer assignment
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(3)
      DO band = 1, nbndlw
        DO j = 1, klev
          DO i = 1, ncol
            aerosol_lw%tau(i,j,band) = aer_tau_lw(i,j,band)
          END DO
        END DO
      END DO
      !$ACC END PARALLEL
      !
      ! RTE-RRTMGP ACC code is synchronous, so need to wait before calling it
      !$ACC WAIT
      IF (ltimer) CALL timer_start(timer_aerosol_lw)
      CALL stop_on_err(aerosol_lw%increment(atmos_lw))
      IF (ltimer) CALL timer_stop (timer_aerosol_lw)
      ! aerosols
      !$ACC END DATA
      DEALLOCATE(aerosol_lw%tau)
      !$ACC END DATA
      IF (ltimer) CALL timer_start(timer_aerosol_lw)
      CALL aerosol_lw%finalize()
      IF (ltimer) CALL timer_stop (timer_aerosol_lw)
    END IF
    !
    !
    ! 4.1.3 Longwave clear-sky fluxes
    !
    IF (lclrsky_lw) THEN
       !
       fluxes_lwcs%flux_up => flx_uplw_clr
       fluxes_lwcs%flux_dn => flx_dnlw_clr
       !
       ! RTE-RRTMGP ACC code is synchronous, so need to wait before calling it
       !$ACC WAIT
       IF (ltimer) CALL timer_start(timer_rte_lw_clrsky)
       CALL stop_on_err(rte_lw(atmos_lw, top_at_1, source_lw, zsemiss, fluxes_lwcs))
       IF (ltimer) CALL timer_stop (timer_rte_lw_clrsky)
       !
    END IF

    ! new cloud optics: allocate memory for cloud optical properties:
    IF (ltimer) CALL timer_start(timer_clouds_bnd_lw)
    CALL stop_on_err(clouds_bnd_lw%alloc_1scl(ncol, klev, &
                     k_dist_lw%get_band_lims_wavenumber()))
    IF (ltimer) CALL timer_stop (timer_clouds_bnd_lw)
    !$ACC DATA CREATE(clouds_bnd_lw)
    !$ACC DATA CREATE(clouds_bnd_lw%tau)
    ! then compute cloud optics

    ! !$ACC update host(zlwp,     ziwp,    re_drop,    re_cryst)
    ! write (0,*) "newcloudsss", sum(zlwp),     sum(ziwp),    sum(re_drop),    sum(re_cryst)
!++jsr, first, detect cloud ice optical depth with zdwp=0,
!       then calculate cloud optical depth
    !$ACC WAIT(1)
    IF (ltimer) CALL timer_start(timer_cloud_optics_lw)
    CALL stop_on_err(cloud_optics_lw%cloud_optics( &
                     zdwp,     ziwp,    re_drop,    re_cryst,   clouds_bnd_lw ))
    IF (ltimer) CALL timer_stop (timer_cloud_optics_lw)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO band = 1, nbndlw
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO j = 1, klev
        DO i = 1, ncol
          tau_ice(i,j) = tau_ice(i,j) + clouds_bnd_lw%tau(i,j,band)
        END DO
      END DO
    END DO
    !$ACC END PARALLEL
!--jsr, calculate cloud optics including ice and water hydrometeors now
!       only these are used in the sequel.    
    !$ACC WAIT(1)
    IF (ltimer) CALL timer_start(timer_cloud_optics_lw)
    CALL stop_on_err(cloud_optics_lw%cloud_optics( &
                     zlwp,     ziwp,    re_drop,    re_cryst,   clouds_bnd_lw ))
    IF (ltimer) CALL timer_stop (timer_cloud_optics_lw)
    ! This will require computing logical masks for ice and liquid clouds
    !   nrghice (ice roughness) is 1, 2, or 3; probably any values is fine
    !
    IF (ltimer) CALL timer_start(timer_clouds_bnd_lw)
    CALL stop_on_err(clouds_bnd_lw%increment(atmos_lw))
    IF (ltimer) CALL timer_stop (timer_clouds_bnd_lw)

    !$ACC WAIT
    !$ACC END DATA
    DEALLOCATE(clouds_bnd_lw%tau)
    !$ACC END DATA
    IF (ltimer) CALL timer_start(timer_clouds_bnd_lw)
    CALL clouds_bnd_lw%finalize()
    IF (ltimer) CALL timer_stop (timer_clouds_bnd_lw)

    ! Snow optics
    ! snow optics using optical properties of cloud ice
    ! allocate memory for snow optical properties:
    IF (ltimer) CALL timer_start(timer_snow_bnd_lw)
    CALL stop_on_err(snow_bnd_lw%alloc_1scl(ncol, klev, &
                     k_dist_lw%get_band_lims_wavenumber()))
    IF (ltimer) CALL timer_stop (timer_snow_bnd_lw)
    !$ACC DATA CREATE(snow_bnd_lw)
    !$ACC DATA CREATE(snow_bnd_lw%tau)
    ! compute snow optics from table of cloud_optics
    !$ACC WAIT(1)
    IF (ltimer) CALL timer_start(timer_cloud_optics_lw)
    CALL stop_on_err(cloud_optics_lw%cloud_optics( &
         zdwp,     zswp,  re_snow,  re_snow,   snow_bnd_lw ))
    IF (ltimer) CALL timer_stop (timer_cloud_optics_lw)
    !++jsr scale tau with reimax/reff_snow for reff_snow > reimax
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO band = 1, nbndlw
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO j = 1, klev
        DO i = 1, ncol
          IF ((1.e6_wp * reff_snow(i, j)) > reimax) THEN
            snow_bnd_lw%tau(i, j, band) = snow_bnd_lw%tau(i, j, band) * reimax / (1.e6_wp * reff_snow(i, j))
          END IF
          tau_snow(i, j) = tau_snow(i, j) + snow_bnd_lw%tau(i, j, band)
        END DO
      END DO
    END DO
    !$ACC END PARALLEL
    !--jsr
    !$ACC WAIT(1)
    IF (ltimer) CALL timer_start(timer_snow_bnd_lw)
    CALL stop_on_err(snow_bnd_lw%increment(atmos_lw))
    IF (ltimer) CALL timer_stop (timer_snow_bnd_lw)
    !$ACC END DATA
    DEALLOCATE(snow_bnd_lw%tau)
    !$ACC END DATA
    IF (ltimer) CALL timer_start(timer_snow_bnd_lw)
    CALL snow_bnd_lw%finalize()
    IF (ltimer) CALL timer_stop (timer_snow_bnd_lw)
    
    !
    ! 4.1.5 Longwave all-sky fluxes
    !
    fluxes_lw%flux_up => flx_uplw
    fluxes_lw%flux_dn => flx_dnlw
    ! RTE-RRTMGP ACC code is synchronous, so need to wait before calling it
    !$ACC WAIT
    IF (ltimer) CALL timer_start(timer_rte_lw_allsky)
    CALL stop_on_err(rte_lw(atmos_lw, top_at_1, source_lw, &
                            zsemiss, fluxes_lw))
    IF (ltimer) CALL timer_stop (timer_rte_lw_allsky)
    !
    ! 4.1.6 End of longwave calculations - free memory
    !
    !$ACC END DATA
    DEALLOCATE(atmos_lw%tau)
    !$ACC END DATA
    IF (ltimer) CALL timer_start(timer_source_lw)
    CALL source_lw%finalize()
    IF (ltimer) CALL timer_stop (timer_source_lw)
    IF (ltimer) CALL timer_start(timer_atmos_lw)
    CALL atmos_lw%finalize()
    IF (ltimer) CALL timer_stop (timer_atmos_lw)
    !
    !-------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------
    !
    ! 4.2 Shortwave calculations
    !
    !-------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------
    !
    ! 4.2.1 Array and type allocation for shortwave
    !--------------------------------
    !
    ! Shortwave gas optical properties and source functions
    !
    IF (ltimer) CALL timer_start(timer_atmos_sw)
    CALL stop_on_err(atmos_sw%alloc_2str(ncol, klev, k_dist_sw))
    IF (ltimer) CALL timer_stop (timer_atmos_sw)
    !$ACC DATA CREATE(atmos_sw)
    !$ACC DATA CREATE(atmos_sw%tau, atmos_sw%ssa, atmos_sw%g) &
    !$ACC   CREATE(toa_flux)

    ! RTE-RRTMGP ACC code is synchronous, so need to wait before calling it
    !$ACC WAIT
    IF (ltimer) CALL timer_start(timer_k_dist_sw)
    CALL stop_on_err(&
       k_dist_sw%gas_optics(play, plev, tlay, &
                            gas_concs, atmos_sw, &
                            toa_flux))
    IF (ltimer) CALL timer_stop (timer_k_dist_sw)
    !toa_flux is output, some flux of rrtmgp, see mo_gas_optics_rrtmgp.F90
    !
    ! 4.2.2 Aerosol optical depth: add to clear-sky, reorder bands
    !
    IF ( lneed_aerosols ) THEN
      IF (ltimer) CALL timer_start(timer_aerosol_sw)
      CALL stop_on_err(aerosol_sw%alloc_2str(ncol, klev, &
                                            k_dist_sw%get_band_lims_wavenumber()))
      IF (ltimer) CALL timer_stop (timer_aerosol_sw)
      !$ACC DATA CREATE(aerosol_sw)
      !$ACC DATA CREATE(aerosol_sw%tau, aerosol_sw%ssa, aerosol_sw%g) &
      !$ACC   PRESENT(aer_tau_sw, aer_ssa_sw, aer_asy_sw)
      !DA TODO: this could be just a pointer assignment
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(3)
      DO band = 1, nbndsw
        DO j = 1, klev
          DO i = 1, ncol
            aerosol_sw%tau(i,j,band) = aer_tau_sw(i,j,band)
            aerosol_sw%ssa(i,j,band) = aer_ssa_sw(i,j,band)
            aerosol_sw%g  (i,j,band) = aer_asy_sw(i,j,band)
          END DO
        END DO
      END DO
      !$ACC END PARALLEL
      ! RTE-RRTMGP ACC code is synchronous, so need to wait before calling it
      !$ACC WAIT
      IF (ltimer) CALL timer_start(timer_aerosol_sw)
      CALL stop_on_err(aerosol_sw%increment(atmos_sw))
      IF (ltimer) CALL timer_stop (timer_aerosol_sw)
      ! aerosol_sw
      !$ACC END DATA
      !$ACC END DATA
      IF (ltimer) CALL timer_start(timer_aerosol_sw)
      CALL aerosol_sw%finalize()
      IF (ltimer) CALL timer_stop (timer_aerosol_sw)
    END IF
    !
    ! 4.2.3 Shortwave clear-sky fluxes
    !
    IF (lclrsky_sw) THEN
       !
       fluxes_swcs%flux_up => flx_upsw_clr
       fluxes_swcs%flux_dn => flx_dnsw_clr
       !
       ! RTE-RRTMGP ACC code is synchronous, so need to wait before calling it
       !$ACC WAIT
       IF (ltimer) CALL timer_start(timer_rte_sw_clrsky)
       CALL stop_on_err(rte_sw(atmos_sw, top_at_1, mu0, toa_flux, albdir, albdif, fluxes_swcs))
       IF (ltimer) CALL timer_stop (timer_rte_sw_clrsky)
       !
    END IF

    ! hack inhom implementation by scaling the condensate water paths
    ! it's important to run this AFTER the longwave
    !!$ACC DATA CREATE(zlwp,ziwp,zswp)
    !!$ACC DATA PRESENT(rad_2d)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO j = 1, klev
      DO i = 1, ncol
        zlwp(i,j) = zlwp(i,j) * rad_2d(i)
        ziwp(i,j) = ziwp(i,j) * inhomi
        zswp(i,j) = zswp(i,j) * inhoms
      END DO
    END DO
    !$ACC END PARALLEL
    
    ! new cloud optics: allocate memory for cloud optical properties:
    IF (ltimer) CALL timer_start(timer_clouds_bnd_sw)
    CALL stop_on_err(clouds_bnd_sw%alloc_2str(ncol, klev, &
                     k_dist_sw%get_band_lims_wavenumber()))
    IF (ltimer) CALL timer_stop (timer_clouds_bnd_sw)
    !$ACC DATA CREATE(clouds_bnd_sw)
    !$ACC DATA CREATE(clouds_bnd_sw%tau, clouds_bnd_sw%ssa, clouds_bnd_sw%g)
    ! then compute cloud optics
!++jsr, first, detect cloud ice optical depth with zdwp=0,
!       then calculate cloud optical depth
    !$ACC WAIT(1)
    IF (ltimer) CALL timer_start(timer_cloud_optics_sw)
    CALL stop_on_err(cloud_optics_sw%cloud_optics( &
                     zdwp,     ziwp,    re_drop,    re_cryst,   clouds_bnd_sw ))
    IF (ltimer) CALL timer_stop (timer_cloud_optics_sw)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO band = 1, nbndsw
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO j = 1, klev
        DO i = 1, ncol
          tau_ice(i,j) = tau_ice(i,j) + clouds_bnd_sw%tau(i,j,band)
        END DO
      END DO
    END DO
    !$ACC END PARALLEL
!--jsr, calculate cloud optics including ice and water hydrometeors now
!       only these are used in the sequel.    
    !$ACC WAIT(1)
    IF (ltimer) CALL timer_start(timer_cloud_optics_sw)
    CALL stop_on_err(cloud_optics_sw%cloud_optics( &
                     zlwp,     ziwp,    re_drop,    re_cryst,   clouds_bnd_sw ))
    IF (ltimer) CALL timer_stop (timer_cloud_optics_sw)
    !
    IF (ltimer) CALL timer_start(timer_clouds_bnd_sw)
    CALL stop_on_err(clouds_bnd_sw%delta_scale()) ! necessary for cases w=g near 1
    CALL stop_on_err(clouds_bnd_sw%increment(atmos_sw))
    IF (ltimer) CALL timer_stop (timer_clouds_bnd_sw)

    !$ACC END DATA
    !$ACC END DATA
    IF (ltimer) CALL timer_start(timer_clouds_bnd_sw)
    CALL clouds_bnd_sw%finalize()
    IF (ltimer) CALL timer_stop (timer_clouds_bnd_sw)
    !
    ! optics for snow
    IF (ltimer) CALL timer_start(timer_snow_bnd_sw)
    CALL stop_on_err(snow_bnd_sw%alloc_2str(ncol, klev, &
                     k_dist_sw%get_band_lims_wavenumber()))
    IF (ltimer) CALL timer_stop (timer_snow_bnd_sw)
    !$ACC DATA CREATE(snow_bnd_sw)
    !$ACC DATA CREATE(snow_bnd_sw%tau, snow_bnd_sw%ssa, snow_bnd_sw%g)
    ! then compute snow optics
    !$ACC WAIT(1)
    IF (ltimer) CALL timer_start(timer_cloud_optics_sw)
    CALL stop_on_err(cloud_optics_sw%cloud_optics( &
                     zdwp,     zswp,  re_snow,  re_snow,   snow_bnd_sw ))
    IF (ltimer) CALL timer_stop (timer_cloud_optics_sw)
    ! delta scale for the case ssa and g close to 1
    IF (ltimer) CALL timer_start(timer_snow_bnd_sw)
    CALL stop_on_err(snow_bnd_sw%delta_scale())
    IF (ltimer) CALL timer_stop (timer_snow_bnd_sw)
    !++jsr scale tau with reimax/reff_snow for reff_snow > reimax
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO band = 1, nbndsw
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO j = 1, klev
        DO i = 1, ncol
          IF ((1.e6_wp * reff_snow(i, j)) > reimax) THEN
            snow_bnd_sw%tau(i, j, band) = snow_bnd_sw%tau(i, j, band) * reimax / (1.e6_wp * reff_snow(i, j))
          END IF
          tau_snow(i, j) = tau_snow(i, j) + snow_bnd_sw%tau(i, j, band)
        END DO
      END DO
    END DO
    !$ACC END PARALLEL
    !--jsr
    ! increment the optcial properties of the atmosphere
    !$ACC WAIT(1)
    IF (ltimer) CALL timer_start(timer_snow_bnd_sw)
    CALL stop_on_err(snow_bnd_sw%increment(atmos_sw))
    IF (ltimer) CALL timer_stop(timer_snow_bnd_sw)
    !$ACC END DATA
    DEALLOCATE(snow_bnd_sw%tau)
    !$ACC END DATA
    IF (ltimer) CALL timer_start(timer_snow_bnd_sw)
    CALL snow_bnd_sw%finalize()
    IF (ltimer) CALL timer_stop (timer_snow_bnd_sw)

    !
    ! 4.2.5 Shortwave all-sky fluxes
    !
    fluxes_sw%flux_up => flx_upsw
    fluxes_sw%flux_dn => flx_dnsw
    fluxes_sw%vis_dn_dir_sfc => vis_dn_dir_sfc
    fluxes_sw%par_dn_dir_sfc => par_dn_dir_sfc
    fluxes_sw%nir_dn_dir_sfc => nir_dn_dir_sfc
    fluxes_sw%vis_dn_dff_sfc => vis_dn_dff_sfc
    fluxes_sw%par_dn_dff_sfc => par_dn_dff_sfc
    fluxes_sw%nir_dn_dff_sfc => nir_dn_dff_sfc
    fluxes_sw%vis_up_sfc => vis_up_sfc
    fluxes_sw%par_up_sfc => par_up_sfc
    fluxes_sw%nir_up_sfc => nir_up_sfc

    CALL set_fractions(fluxes_sw, atmos_sw, psctm, ssi_factor)
    ! RTE-RRTMGP ACC code is synchronous, so need to wait before calling it
    !$ACC WAIT
    IF (ltimer) CALL timer_start(timer_rte_sw_allsky)
    CALL stop_on_err(rte_sw(atmos_sw, top_at_1, &
                            mu0, toa_flux, albdir, albdif, &
                            fluxes_sw))
    IF (ltimer) CALL timer_stop (timer_rte_sw_allsky)

    !
    ! 4.2.6 End of shortwave calculations - free memory
    !
    !$ACC END DATA
    !$ACC END DATA
    CALL atmos_sw%finalize()
        
#ifdef RRTMGP_MERGE_DEBUG
!$OMP CRITICAL (write_record)
    CALL write_record_interface_aes(nproma, pcos_mu0, daylght_frc, &
      alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif, &
      tk_sfc, zf, zh, dz, pp_fl, pp_hl, tk_fl, tk_hl, &
      play, plev, tlay, tlev, &
      xvmr_vap, xvmr_co2, xvmr_ch4, xvmr_o2, xvmr_o3, xvmr_n2o, cdnc, &
      cld_frc, &
      flx_dnlw_clr, flx_uplw_clr, flx_dnsw_clr, flx_upsw_clr, &
      flx_dnlw, flx_uplw, flx_dnsw, flx_upsw, &
      vis_dn_dir_sfc, par_dn_dir_sfc, nir_dn_dir_sfc, &
      vis_dn_dff_sfc, par_dn_dff_sfc, nir_dn_dff_sfc, &
      vis_up_sfc,     par_up_sfc,     nir_up_sfc      )
!$OMP END CRITICAL (write_record)
#endif

  !$ACC END DATA

    IF (ltimer) CALL timer_stop(timer_rte_rrtmgp_int_onb)

  END SUBROUTINE rte_rrtmgp_interface_onBlock
  ! ----------------------------------------------------------------------------
  SUBROUTINE shift_and_call_rte_rrtmgp_interface_onBlock(    &
    & lclrsky_lw,     lclrsky_sw,                     &
    & inhom_lts,      inhom_lts_max,                  &
    & jcs,            jce,                            &
    &                 klev,                           &
    !
    & psctm,          ssi_factor,                     &
    & laland,         laglac,                         &
    & pcos_mu0,       daylght_frc,                    &
    & alb_vis_dir,    alb_nir_dir,                    &
    & alb_vis_dif,    alb_nir_dif,                    &
    & emissivity,                                     &
    & zf,             zh,             dz,             &
    & pp_sfc,         pp_fl,          pp_hl,          &
    & tk_sfc,         tk_fl,          tk_hl,          &
    & rad_2d,                                         &
    & xvmr_vap,       xm_liq,                         &
    & xm_ice,         reff_ice,       tau_ice,        &
    & reff_snow,      tau_snow,                       &
    & cdnc,           xc_frc,         xm_snw,         &
    & xvmr_co2,       xvmr_ch4,       xvmr_n2o,       &
    & xvmr_cfc,       xvmr_o3,        xvmr_o2,        &
    & aer_tau_lw,                                     &
    & aer_tau_sw,     aer_ssa_sw,     aer_asy_sw,     &
    !
    & lw_upw,         lw_upw_clr,                     &
    & lw_dnw,         lw_dnw_clr,                     &
    & sw_upw,         sw_upw_clr,                     &
    & sw_dnw,         sw_dnw_clr,                     &
    & vis_dn_dir_sfc, par_dn_dir_sfc, nir_dn_dir_sfc, &
    & vis_dn_dff_sfc, par_dn_dff_sfc, nir_dn_dff_sfc, &
    & vis_up_sfc,     par_up_sfc,     nir_up_sfc      )

 LOGICAL,INTENT(IN)  :: lclrsky_lw                    !< flag for LW clear-sky computations
 LOGICAL,INTENT(IN)  :: lclrsky_sw                    !< flag for SW clear-sky computations
 LOGICAL,INTENT(IN)  :: inhom_lts
 REAL(wp),INTENT(IN) :: inhom_lts_max

 INTEGER,INTENT(IN)  :: &
      & jcs,            & !< cell/column index, start
      & jce,            & !< cell/column index, end
      & klev              !< number of full levels

 REAL(wp),INTENT(IN) :: psctm                         !< orbit and time dependent solar constant for radiation time step
 REAL(wp),INTENT(IN) :: ssi_factor(:)                 !< fraction of TSI in the 14 RRTM SW bands

 LOGICAL,INTENT(IN) :: &
      & laland(:),   & !< land sea mask, land=.true.
      & laglac(:)      !< glacier mask, glacier=.true.

 REAL(WP),INTENT(IN)  ::    &
      & pcos_mu0(:),      & !< mu0 for solar zenith angle
      & daylght_frc(:),   & !< daylight fraction; with diurnal cycle 0 or 1, with zonal mean in [0,1]
      & alb_vis_dir(:),   & !< surface albedo for vis range and dir light
      & alb_nir_dir(:),   & !< surface albedo for NIR range and dir light
      & alb_vis_dif(:),   & !< surface albedo for vis range and dif light
      & alb_nir_dif(:),   & !< surface albedo for NIR range and dif light
      & emissivity(:),    & !< surface longwave emissivity
      & zf(:,:),          & !< geometric height at full level in m
      & zh(:,:),          & !< geometric height at half level in m
      & dz(:,:),          & !< geometric height thickness in m
      & pp_sfc(:),        & !< surface pressure in Pa
      & pp_fl(:,:),       & !< full level pressure in Pa
      & pp_hl(:,:),       & !< full level pressure in Pa
      & tk_sfc(:),        & !< surface temperature in K
      & tk_fl(:,:),       & !< full level temperature in K
      & tk_hl(:,:),       & !< half level temperature in K
      & xvmr_vap(:,:),    & !< water vapor volume mixing ratio
      & xm_liq(:,:),      & !< cloud water mass in kg/m2
      & xm_ice(:,:),      & !< cloud ice   mass in kg/m2
      & cdnc(:,:),        & !< cloud nuclei concentration
      & xc_frc(:,:),      & !< fractional cloud cover
      & xm_snw(:,:),      & !< snow        mass in kg/m2
      & xvmr_co2(:,:),    & !< co2 volume mixing ratio
      & xvmr_ch4(:,:),    & !< ch4 volume mixing ratio
      & xvmr_n2o(:,:),    & !< n2o volume mixing ratio
      & xvmr_cfc(:,:,:),  & !< cfc volume mixing ratio (kbdim,klev,2)
      & xvmr_o3(:,:),     & !< o3  volume mixing ratio
      & xvmr_o2(:,:),     & !< o2  volume mixing ratio
      & aer_tau_lw(:,:,:),& !< aerosol optical depth, longwave (ncol, nlay, nbndlw)
      & aer_tau_sw(:,:,:),& !< aerosol optical depth,            shortwave (ncol, nlay, nbndlw)
      & aer_ssa_sw(:,:,:),& !< aerosol single-scattering albedo, shortwave (ncol, nlay, nbndlw)
      & aer_asy_sw(:,:,:),& !< aerosol asymetry parameter,       shortwave (ncol, nlay, nbndlw)
      & reff_ice(:,:),    & !< effective radius of cloud ice m
      & reff_snow(:,:)      !< effective radius of snow m

 REAL (wp), INTENT (INOUT) :: &
      & tau_ice(:,:),     & !< optical depth of cloud ice integrated over bands
      & tau_snow(:,:)       !< optical depth of snow integrated over bands
 
 REAL (wp), INTENT (INOUT) :: &
      & rad_2d(:)           !< arbitrary 2d-field in radiation for output


 REAL (wp), TARGET, INTENT (INOUT) ::       &
      & lw_upw    (:,:), & !<   upward LW flux profile, all sky
      & lw_upw_clr(:,:), & !<   upward LW flux profile, clear sky
      & lw_dnw    (:,:), & !< downward LW flux profile, all sky
      & lw_dnw_clr(:,:), & !< downward LW flux profile, clear sky
      & sw_upw    (:,:), & !<   upward SW flux profile, all sky
      & sw_upw_clr(:,:), & !<   upward SW flux profile, clear sky
      & sw_dnw    (:,:), & !< downward SW flux profile, all sky
      & sw_dnw_clr(:,:)    !< downward SW flux profile, clear sky

 REAL (wp), TARGET, INTENT (INOUT) :: &
      & vis_dn_dir_sfc(:) , & !< Diffuse downward flux surface visible radiation
      & par_dn_dir_sfc(:) , & !< Diffuse downward flux surface PAR
      & nir_dn_dir_sfc(:) , & !< Diffuse downward flux surface near-infrared radiation
      & vis_dn_dff_sfc(:) , & !< Direct  downward flux surface visible radiation
      & par_dn_dff_sfc(:) , & !< Direct  downward flux surface PAR
      & nir_dn_dff_sfc(:) , & !< Direct  downward flux surface near-infrared radiation
      & vis_up_sfc    (:) , & !< Upward  flux surface visible radiation
      & par_up_sfc    (:) , & !< Upward  flux surface PAR
      & nir_up_sfc    (:)     !< Upward  flux surface near-infrared radiation

 INTEGER :: ncol !< number of columns needed

 ! Shifted input arguments
 !
 REAL(wp)  ::                                    &
      & s_zf             (jce-jcs+1,klev),       & !< geometric height at full level in m
      & s_zh             (jce-jcs+1,klev+1),     & !< geometric height at half level in m
      & s_dz             (jce-jcs+1,klev),       & !< geometric height thickness in m
      & s_pp_fl          (jce-jcs+1,klev),       & !< full level pressure in Pa
      & s_pp_hl          (jce-jcs+1,klev+1),     & !< full level pressure in Pa
      & s_tk_fl          (jce-jcs+1,klev),       & !< full level temperature in K
      & s_tk_hl          (jce-jcs+1,klev+1),     & !< half level temperature in K
      & s_xvmr_vap       (jce-jcs+1,klev),       & !< water vapor volume mixing ratio
      & s_xm_liq         (jce-jcs+1,klev),       & !< cloud water mass in kg/m2
      & s_xm_ice         (jce-jcs+1,klev),       & !< cloud ice   mass in kg/m2
      & s_reff_ice       (jce-jcs+1,klev),       & !< cloud ice effective radius
      & s_tau_ice        (jce-jcs+1,klev),       & !< optical depth of cloud ice integrated over bands
      & s_reff_snow      (jce-jcs+1,klev),       & !< snow effective radius
      & s_tau_snow       (jce-jcs+1,klev),       & !< optical depth of snow integrated over bands
      & s_cdnc           (jce-jcs+1,klev),       & !< cloud nuclei concentration
      & s_xc_frc         (jce-jcs+1,klev),       & !< fractional cloud cover
      & s_xm_snw         (jce-jcs+1,klev),       & !< snow        mass in kg/m2
      & s_xvmr_co2       (jce-jcs+1,klev),       & !< co2 volume mixing ratio
      & s_xvmr_ch4       (jce-jcs+1,klev),       & !< ch4 volume mixing ratio
      & s_xvmr_n2o       (jce-jcs+1,klev),       & !< n2o volume mixing ratio
      & s_xvmr_cfc       (jce-jcs+1,klev,2),     & !< cfc volume mixing ratio
      & s_xvmr_o3        (jce-jcs+1,klev),       & !< o3  volume mixing ratio
      & s_xvmr_o2        (jce-jcs+1,klev)          !< o2  volume mixing ratio

 REAL(wp), ALLOCATABLE ::    & 
      & s_aer_tau_lw(:,:,:), &
      & s_aer_tau_sw(:,:,:), &
      & s_aer_ssa_sw(:,:,:), &
      & s_aer_asy_sw(:,:,:)

 ! Shifted output arguments
 ! - Note: The size of the 2nd dimension of the "clr" fields can be different from klev+1
 !         Therefore SIZE(.,2) is used her for the 2nd dimension of all fields.
 !
 REAL(wp)  ::                              &
      & s_lw_upw       (jce-jcs+1,SIZE(lw_upw    ,2)), & !<   upward LW flux profile, all sky
      & s_lw_upw_clr   (jce-jcs+1,SIZE(lw_upw_clr,2)), & !<   upward LW flux profile, clear sky
      & s_lw_dnw       (jce-jcs+1,SIZE(lw_dnw    ,2)), & !< downward LW flux profile, all sky
      & s_lw_dnw_clr   (jce-jcs+1,SIZE(lw_dnw_clr,2)), & !< downward LW flux profile, clear sky
      & s_sw_upw       (jce-jcs+1,SIZE(sw_upw    ,2)), & !<   upward SW flux profile, all sky
      & s_sw_upw_clr   (jce-jcs+1,SIZE(sw_upw_clr,2)), & !<   upward SW flux profile, clear sky
      & s_sw_dnw       (jce-jcs+1,SIZE(sw_dnw    ,2)), & !< downward SW flux profile, all sky
      & s_sw_dnw_clr   (jce-jcs+1,SIZE(sw_dnw_clr,2))    !< downward SW flux profile, clear sky

  ! Shift input arguments that would be non-contiguous when sliced
  !
  ncol = jce-jcs+1

  !$ACC DATA CREATE(s_zf, s_zh, s_dz) &
  !$ACC   CREATE(s_pp_fl, s_pp_hl) &
  !$ACC   CREATE(s_tk_fl, s_tk_hl) &
  !$ACC   CREATE(s_xvmr_vap, s_xm_liq) &
  !$ACC   CREATE(s_xm_ice, s_reff_ice, s_tau_ice, s_reff_snow, s_tau_snow) &
  !$ACC   CREATE(s_cdnc, s_xc_frc, s_xm_snw) &
  !$ACC   CREATE(s_xvmr_co2, s_xvmr_ch4, s_xvmr_n2o) &
  !$ACC   CREATE(s_xvmr_cfc, s_xvmr_o3, s_xvmr_o2) &
  !$ACC   CREATE(s_lw_upw, s_lw_upw_clr) &
  !$ACC   CREATE(s_lw_dnw, s_lw_dnw_clr) &
  !$ACC   CREATE(s_sw_upw, s_sw_upw_clr) &
  !$ACC   CREATE(s_sw_dnw, s_sw_dnw_clr)

  ! (ncol, klev)
  !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
  s_zf          (1:ncol,:)   = zf          (jcs:jce,:)
  s_dz          (1:ncol,:)   = dz          (jcs:jce,:)
  s_pp_fl       (1:ncol,:)   = pp_fl       (jcs:jce,:)
  s_tk_fl       (1:ncol,:)   = tk_fl       (jcs:jce,:)
  s_xvmr_vap    (1:ncol,:)   = xvmr_vap    (jcs:jce,:)
  s_xm_liq      (1:ncol,:)   = xm_liq      (jcs:jce,:)
  s_xm_ice      (1:ncol,:)   = xm_ice      (jcs:jce,:)
  s_reff_ice    (1:ncol,:)   = reff_ice    (jcs:jce,:)
  s_tau_ice     (1:ncol,:)   = tau_ice     (jcs:jce,:)
  s_reff_snow   (1:ncol,:)   = reff_snow   (jcs:jce,:)
  s_tau_snow    (1:ncol,:)   = tau_snow    (jcs:jce,:)
  s_cdnc        (1:ncol,:)   = cdnc        (jcs:jce,:)
  s_xc_frc      (1:ncol,:)   = xc_frc      (jcs:jce,:)
  s_xm_snw      (1:ncol,:)   = xm_snw      (jcs:jce,:)
  s_xvmr_co2    (1:ncol,:)   = xvmr_co2    (jcs:jce,:)
  s_xvmr_ch4    (1:ncol,:)   = xvmr_ch4    (jcs:jce,:)
  s_xvmr_n2o    (1:ncol,:)   = xvmr_n2o    (jcs:jce,:)
  s_xvmr_o3     (1:ncol,:)   = xvmr_o3     (jcs:jce,:)
  s_xvmr_o2     (1:ncol,:)   = xvmr_o2     (jcs:jce,:)
  !$ACC END KERNELS

  ! (ncol, klev+1)
  !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
  s_zh          (1:ncol,:)   = zh          (jcs:jce,:)
  s_pp_hl       (1:ncol,:)   = pp_hl       (jcs:jce,:)
  s_tk_hl       (1:ncol,:)   = tk_hl       (jcs:jce,:)
  !$ACC END KERNELS

  ! (ncol, klev, 2)
  !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
  s_xvmr_cfc    (1:ncol,:,:) = xvmr_cfc    (jcs:jce,:,:)
  !$ACC END KERNELS

  IF ( lneed_aerosols ) THEN
    ! Aerosols are present, irad_aero /= 0
    !      
    ALLOCATE( s_aer_tau_lw(jce-jcs+1,klev,k_dist_lw%get_nband()), &
              s_aer_tau_sw(jce-jcs+1,klev,k_dist_sw%get_nband()), &
              s_aer_ssa_sw(jce-jcs+1,klev,k_dist_sw%get_nband()), &
              s_aer_asy_sw(jce-jcs+1,klev,k_dist_sw%get_nband())  )
    !
    !$ACC ENTER DATA CREATE(s_aer_tau_lw, s_aer_tau_sw, s_aer_ssa_sw, s_aer_asy_sw)
    !
    ! (ncol, klev, nbndlw)
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
    s_aer_tau_lw  (1:ncol,:,:) = aer_tau_lw(jcs:jce,:,:)
    !$ACC END KERNELS

    ! (ncol, klev, nbndsw)
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
    s_aer_tau_sw  (1:ncol,:,:) = aer_tau_sw(jcs:jce,:,:)
    s_aer_ssa_sw  (1:ncol,:,:) = aer_ssa_sw(jcs:jce,:,:)
    s_aer_asy_sw  (1:ncol,:,:) = aer_asy_sw(jcs:jce,:,:)
    !$ACC END KERNELS
  ELSE
    ! allocate dummy zero-size arrays
    ALLOCATE( s_aer_tau_lw(1,1,0), &
              s_aer_tau_sw(1,1,0), &
              s_aer_ssa_sw(1,1,0), &
              s_aer_asy_sw(1,1,0)  )
  END IF

  ! Call radiation with shifted input arguments and receive shifted output arguments
  !
  CALL rte_rrtmgp_interface_onBlock(                                                 &
      & lclrsky_lw,               lclrsky_sw,                                        &
      & inhom_lts,                inhom_lts_max,                                     &
      & ncol,                     klev,                                              &
      !
      & psctm,                    ssi_factor,                                        &
      & laland     (jcs:jce),     laglac     (jcs:jce),                              &
      & pcos_mu0   (jcs:jce),     daylght_frc(jcs:jce),                              &
      & alb_vis_dir(jcs:jce),     alb_nir_dir(jcs:jce),                              &
      & alb_vis_dif(jcs:jce),     alb_nir_dif(jcs:jce),                              &
      & emissivity (jcs:jce),                                                        &
      & s_zf(:,:),                s_zh(:,:),                s_dz(:,:),               &
      & pp_sfc     (jcs:jce),     s_pp_fl(:,:),             s_pp_hl(:,:),            &
      & tk_sfc     (jcs:jce),     s_tk_fl(:,:),             s_tk_hl(:,:),            &
      & rad_2d     (jcs:jce),                                                        &
      & s_xvmr_vap(:,:),          s_xm_liq(:,:),                                     &
      & s_xm_ice(:,:),            s_reff_ice(:,:),          s_tau_ice(:,:),          &
      & s_reff_snow(:,:),         s_tau_snow(:,:),                                   &
      & s_cdnc(:,:),              s_xc_frc(:,:),            s_xm_snw(:,:),           &
      & s_xvmr_co2(:,:),          s_xvmr_ch4(:,:),          s_xvmr_n2o(:,:),         &
      & s_xvmr_cfc(:,:,:),        s_xvmr_o3(:,:),           s_xvmr_o2(:,:),          &
      & s_aer_tau_lw(:,:,:),                                                         &
      & s_aer_tau_sw(:,:,:),      s_aer_ssa_sw(:,:,:),      s_aer_asy_sw(:,:,:),     &
      !     
      & s_lw_upw(:,:),            s_lw_upw_clr(:,:),                                 &
      & s_lw_dnw(:,:),            s_lw_dnw_clr(:,:),                                 &
      & s_sw_upw(:,:),            s_sw_upw_clr(:,:),                                 &
      & s_sw_dnw(:,:),            s_sw_dnw_clr(:,:),                                 &
      & vis_dn_dir_sfc(jcs:jce),  par_dn_dir_sfc(jcs:jce),  nir_dn_dir_sfc(jcs:jce), &
      & vis_dn_dff_sfc(jcs:jce),  par_dn_dff_sfc(jcs:jce),  nir_dn_dff_sfc(jcs:jce), &
      & vis_up_sfc    (jcs:jce),  par_up_sfc    (jcs:jce),  nir_up_sfc    (jcs:jce)  )

  ! Shift output arguments
  !
  ! (ncol, klev+1)
  !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
  tau_snow       (jcs:jce,:) = s_tau_snow       (1:ncol,:)
  tau_ice        (jcs:jce,:) = s_tau_ice        (1:ncol,:)
  lw_upw         (jcs:jce,:) = s_lw_upw         (1:ncol,:)
  lw_dnw         (jcs:jce,:) = s_lw_dnw         (1:ncol,:)
  sw_upw         (jcs:jce,:) = s_sw_upw         (1:ncol,:)
  sw_dnw         (jcs:jce,:) = s_sw_dnw         (1:ncol,:)
  !$ACC END KERNELS
  !
  IF (lclrsky_lw) THEN
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
    lw_upw_clr   (jcs:jce,:) = s_lw_upw_clr     (1:ncol,:)
    lw_dnw_clr   (jcs:jce,:) = s_lw_dnw_clr     (1:ncol,:)
    !$ACC END KERNELS
  END IF
  !
  IF (lclrsky_sw) THEN
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
    sw_upw_clr   (jcs:jce,:) = s_sw_upw_clr     (1:ncol,:)
    sw_dnw_clr   (jcs:jce,:) = s_sw_dnw_clr     (1:ncol,:)
    !$ACC END KERNELS
  END IF

  !$ACC WAIT(1)
  !$ACC EXIT DATA DELETE(s_aer_tau_lw, s_aer_tau_sw, s_aer_ssa_sw, s_aer_asy_sw) IF(lneed_aerosols)
  !$ACC END DATA
END SUBROUTINE shift_and_call_rte_rrtmgp_interface_onBlock

SUBROUTINE reorient_3d_wrt2 (field)
  REAL(wp), INTENT(INOUT)  :: field(:,:,:)
  INTEGER                  :: nl1, nl2, nl3, il1, il2, il3, idx1, idx2
  REAL(WP)                 :: tmp

  nl1=SIZE(field,1)
  nl2=SIZE(field,2)
  nl3=SIZE(field,3)

  !$ACC PARALLEL PRESENT(field) ASYNC(1)
  !$ACC LOOP GANG VECTOR COLLAPSE(3)
  DO il3 = 1, nl3
    DO il2 = 1, nl2/2
      DO il1 = 1, nl1
        idx1 = il2
        idx2 = nl2 - il2 + 1
        tmp = field(il1,idx1,il3)
        field(il1,idx1,il3) = field(il1,idx2,il3)
        field(il1,idx2,il3) = tmp
      END DO
    END DO
  END DO
  !$ACC END PARALLEL
END SUBROUTINE reorient_3d_wrt2

SUBROUTINE rearrange_bands2rrtmgp(nproma, klev, nbnd, field)
  INTEGER,  INTENT(IN)    :: nproma, klev, nbnd
  REAL(wp), INTENT(INOUT) :: field(nproma,klev,nbnd) 

#ifndef _OPENACC
  INTEGER  :: i
#else
  REAL(wp) :: last
  INTEGER  :: jk, jl, jband
#endif
  
#ifndef _OPENACC
  field(:,:,:) = field(:,:,[nbnd, (i, i = 1, nbnd-1)])
#else
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO jk=1,klev
    DO jl=1,nproma
      last = field(jl,jk,nbnd)

      !$ACC LOOP SEQ
      DO jband=nbnd, 2, -1
        field(jl,jk,jband) = field(jl,jk,jband-1)
      END DO

      field(jl,jk,1) = last
    END DO
  END DO
  !$ACC END PARALLEL
#endif
END SUBROUTINE rearrange_bands2rrtmgp

END MODULE mo_rte_rrtmgp_interface
