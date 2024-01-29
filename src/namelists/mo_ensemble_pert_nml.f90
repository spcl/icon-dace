! Namelist for perturbing nwp physics
!
! These subroutines are called by read_atmo_namelists and set the ranges
! for ensemble physics perturbations
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

MODULE mo_ensemble_pert_nml

  USE mo_kind,                ONLY: wp
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_ensemble_pert_config,ONLY: config_range_gkwake    => range_gkwake,    &
    &                               config_range_gkdrag    => range_gkdrag,    &
    &                               config_range_gfrcrit   => range_gfrcrit,   &
    &                               config_range_gfluxlaun => range_gfluxlaun, &
    &                               config_range_zvz0i     => range_zvz0i,     &
    &                               config_range_rain_n0fac => range_rain_n0fac, &
    &                               config_range_entrorg   => range_entrorg,   &  
    &                               config_range_rdepths   => range_rdepths,   &  
    &                               config_range_rprcon    => range_rprcon,    &  
    &                               config_range_capdcfac_et => range_capdcfac_et, &  
    &                               config_range_capdcfac_tr => range_capdcfac_tr, &  
    &                               config_range_lowcapefac  => range_lowcapefac,  &  
    &                               config_range_negpblcape  => range_negpblcape,  &  
    &                               config_range_minsnowfrac => range_minsnowfrac, &
    &                               config_range_c_soil    => range_c_soil,    &
    &                               config_range_cwimax_ml => range_cwimax_ml, &
    &                               config_range_rhebc     => range_rhebc,     &  
    &                               config_range_texc      => range_texc,      &
    &                               config_range_qexc      => range_qexc,      &
    &                               config_range_box_liq   => range_box_liq,   &
    &                               config_range_box_liq_asy => range_box_liq_asy, &
    &                               config_range_thicklayfac => range_thicklayfac, &
    &                               config_range_fac_ccqc  => range_fac_ccqc,  &
    &                               config_range_tkhmin    => range_tkhmin,    &  
    &                               config_range_tkmmin    => range_tkmmin,    &
    &                               config_range_turlen    => range_turlen,    &
    &                               config_range_a_hshr    => range_a_hshr,    &
    &                               config_range_a_stab    => range_a_stab,    &
    &                               config_range_c_diff    => range_c_diff,    &
    &                               config_range_q_crit    => range_q_crit,    &
    &                               config_range_tkred_sfc => range_tkred_sfc, &
    &                               config_range_rlam_heat => range_rlam_heat, &
    &                               config_range_charnock  => range_charnock,  &
    &                               config_range_lhn_coef  => range_lhn_coef,  &
    &                               config_range_lhn_artif_fac => range_lhn_artif_fac, &
    &                               config_range_fac_lhn_down => range_fac_lhn_down, &
    &                               config_range_fac_lhn_up => range_fac_lhn_up, &
    &                               config_range_z0_lcc    => range_z0_lcc,    &
    &                               config_range_rootdp    => range_rootdp,    &
    &                               config_range_rsmin     => range_rsmin,     &
    &                               config_range_laimax    => range_laimax,    &
    &                               config_stdev_sst_pert  => stdev_sst_pert,  &
    &                               config_itype_pert_gen  => itype_pert_gen,  &
    &                               config_timedep_pert    => timedep_pert,    &
    &                               config_use_ensemble_pert => use_ensemble_pert

  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_ensemble_pert_namelist


  !--------------------------------------!
  ! ensemble_pert_nml namelist variables !
  !--------------------------------------!

  REAL(wp) :: &                    !< low level wake drag constant
    &  range_gkwake

  REAL(wp) :: &                    !< gravity wave drag constant
    &  range_gkdrag

  REAL(wp) :: &                    !< critical Froude number used for computing blocking layer depth
    &  range_gfrcrit

  REAL(wp) :: &                    !< gravity wave flux emission
    &  range_gfluxlaun

  REAL(wp) :: &                    !< Terminal fall velocity of ice 
    &  range_zvz0i

  REAL(wp) :: &                    !< Tuning factor for intercept parameter of raindrop size distribution 
    &  range_rain_n0fac

  REAL(wp) :: &                    !< Entrainment parameter for deep convection valid at dx=20 km 
    &  range_entrorg

  REAL(wp) :: &                    !< Maximum shallow convection depth 
    &  range_rdepths

  REAL(wp) :: &                    !< Entrainment parameter for deep convection valid at dx=20 km 
    &  range_rprcon

  REAL(wp) :: &                    !< Fraction of CAPE diurnal cycle correction applied in the extratropics
    &  range_capdcfac_et            ! (relevant only if icapdcycl = 3)

  REAL(wp) :: &                    !< Fraction of CAPE diurnal cycle correction applied in the tropics
    &  range_capdcfac_tr            ! (relevant only if icapdcycl = 3)

  REAL(wp) :: &                    !< Tuning factor for reducing the diurnal cycle correction in low-cape situations
    &  range_lowcapefac            ! (relevant only if icapdcycl = 3)

  REAL(wp) :: &                    !< Minimum allowed negative PBL cape in diurnal cycle correction
    &  range_negpblcape            ! (relevant only if icapdcycl = 3)

  REAL(wp) :: &                    !< RH thresholds for evaporation below cloud base
    &  range_rhebc

  REAL(wp) :: &                    !< Excess value for temperature used in test parcel ascent
    &  range_texc

  REAL(wp) :: &                    !< Excess value for mixing ratio used in test parcel ascent
    &  range_qexc

  REAL(wp) :: &                    !< Minimum value to which the snow cover fraction is artificially reduced
    &  range_minsnowfrac           !  in case of melting show (in case of idiag_snowfrac = 20)

  REAL(wp) :: &                    !< Fraction of surface area available for bare soil evaporation
    &  range_c_soil

  REAL(wp) :: &                    !< Capacity of interception storage (multiplicative perturbation)
    &  range_cwimax_ml

  REAL(wp) :: &                    !< Box width for liquid clouds assumed in the cloud cover scheme
    &  range_box_liq               ! (in case of inwp_cldcover = 1)

  REAL(wp) :: &                    !< factor [1/m] for increasing the box with for layer thicknesses exceeding 150 m
    &  range_thicklayfac           ! (in case of inwp_cldcover = 1)

  REAL(wp) :: &                    !< Asymmetry factor for sub-grid scale liquid cloud distribution
    &  range_box_liq_asy           ! (in case of inwp_cldcover = 1)

  REAL(wp) :: &                    !< Factor for CLC-QC relationship in cloud cover scheme
    &  range_fac_ccqc              ! (in case of inwp_cldcover = 1)

  REAL(wp) :: &                    !< Minimum vertical diffusion for heat/moisture 
    &  range_tkhmin

  REAL(wp) :: &                    !< Minimum vertical diffusion for momentum 
    &  range_tkmmin

  REAL(wp) :: &                    !< Perturbation of reduction of minimum diffusion coefficients near the surface 
    &  range_tkred_sfc

  REAL(wp) :: &                    !< Laminar transport resistance parameter 
    &  range_rlam_heat

  REAL(wp) :: &                    !< Maximum turbulent mixing length scale 
    &  range_turlen

  REAL(wp) :: &                    !< Scaling factor for extended horizontal shear term in turbulence scheme 
    &  range_a_hshr

  REAL(wp) :: &                    !< Scaling factor for stability correction in turbulence scheme 
    &  range_a_stab

  REAL(wp) :: &                    !< Length scale factor for vertical diffusion in turbulence scheme 
    &  range_c_diff

  REAL(wp) :: &                    !< Critical value for normalized supersaturation in turbulent cloud scheme
    &  range_q_crit

  REAL(wp) :: &                    !< Upper and lower bound of wind-speed dependent Charnock parameter 
    &  range_charnock

  REAL(wp) :: &                    !< Scaling factor for latent heat nudging increments
    &  range_lhn_coef

  REAL(wp) :: &                    !< Scaling factor for artificial heating profile in latent heat nudging
    &  range_lhn_artif_fac

  REAL(wp) :: &                    !< Lower limit for reduction of existing latent heating in LHN
    &  range_fac_lhn_down

  REAL(wp) :: &                    !< Upper limit for increase of existing latent heating in LHN
    &  range_fac_lhn_up

  REAL(wp) :: &                    !< Roughness length attributed to land-cover class 
    &  range_z0_lcc

  REAL(wp) :: &                    !< Root depth related to land-cover class
    &  range_rootdp

  REAL(wp) :: &                    !< Minimum stomata resistance related to land-cover class
    &  range_rsmin

  REAL(wp) :: &                    !< Maximum leaf area index related to land-cover class
    &  range_laimax

  REAL(wp) :: &                    !< Standard deviation of SST perturbations specified in SST analysis (K)
    &  stdev_sst_pert              !  this switch controls a correction term compensating the systematic
                                   !  increase of evaporation related to the SST perturbations

  INTEGER :: itype_pert_gen        !< type of ensemble perturbation generation
  INTEGER :: timedep_pert          !< time dependence of ensemble perturbations
  LOGICAL :: use_ensemble_pert     !< main switch

  NAMELIST/ensemble_pert_nml/ use_ensemble_pert, range_gkwake, range_gkdrag, range_gfluxlaun, range_zvz0i, &
    &                         range_entrorg, range_capdcfac_et, range_box_liq, range_tkhmin, range_tkmmin, &
    &                         range_rlam_heat, range_rhebc, range_texc, range_minsnowfrac, range_z0_lcc,   &
    &                         range_rootdp, range_rsmin, range_laimax, range_charnock, range_tkred_sfc,    &
    &                         range_gfrcrit, range_c_soil, range_cwimax_ml, range_capdcfac_tr,             &
    &                         range_lowcapefac, range_negpblcape, stdev_sst_pert, itype_pert_gen,          &
    &                         timedep_pert, range_a_stab, range_c_diff, range_q_crit, range_box_liq_asy,   &
    &                         range_rdepths, range_turlen, range_rain_n0fac, range_a_hshr, range_qexc,     &
    &                         range_rprcon, range_thicklayfac, range_lhn_coef, range_lhn_artif_fac,        &
    &                         range_fac_lhn_down, range_fac_lhn_up, range_fac_ccqc


CONTAINS


  !-------------------------------------------------------------------------
  !
  !! Read Namelist for NWP ensemble perturbations. 
  !!
  !! This subroutine 
  !! - reads the Namelist for NWP ensemble perturbations
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)   
  !!
  SUBROUTINE read_ensemble_pert_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: iunit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_ensemble_pert_nml: read_ensemble_pert_namelist'

    !-----------------------
    ! 1. default settings   
    !-----------------------

    ! Ranges for ensemble perturbations:

    ! SSO tuning
    range_gkwake     = 1.5_wp       ! low-level blocking parameter
    range_gkdrag     = 0.04_wp      ! parameter for wave drag deposition in the atmosphere
    range_gfrcrit    = 0.1_wp       ! parameter for critical Froude number
    !
    ! GWD tuning
    range_gfluxlaun  = 0.75e-3_wp   ! scaling parameter for GWD flux production
    !
    ! grid scale microphysics
    range_zvz0i      = 0.25_wp      ! scaling for cloud ice sedimentation speed
    range_rain_n0fac = 4._wp        ! multiplicative change of intercept parameter of raindrop size distribution
    !
    ! convection
    range_entrorg    = 0.2e-3_wp    ! entrainment parameter for deep convection
    range_rdepths    = 5.e3_wp      ! maximum allowed shallow convection depth (Pa)
    range_rprcon     = 0.25e-3_wp   ! Coefficient for conversion of cloud water into precipitation
    range_capdcfac_et = 0.75_wp     ! fraction of CAPE diurnal cycle correction applied in the extratropics
    range_capdcfac_tr = 0.75_wp     ! fraction of CAPE diurnal cycle correction applied in the tropics
    range_lowcapefac = 0.5_wp       ! Tuning factor for reducing the diurnal cycle correction in low-cape situations
    range_negpblcape = 500._wp      ! Minimum allowed negative PBL cape in diurnal cycle correction
    range_rhebc      = 0.05_wp      ! RH thresholds for evaporation below cloud base
    range_texc       = 0.05_wp      ! Excess value for temperature used in test parcel ascent
    range_qexc       = 0.005_wp     ! Excess value for moisture (QV) used in test parcel ascent
    !
    ! cloud cover
    range_box_liq     = 0.01_wp      ! box width scale of liquid clouds
    range_box_liq_asy = 0.25_wp      ! Asymmetry factor for sub-grid scale liquid cloud distribution
    range_thicklayfac = 0.0025_wp    ! factor [1/m] for increasing the box with for layer thicknesses exceeding 150 m
    range_fac_ccqc    = 4.0_wp       ! Factor for CLC-QC relationship in cloud cover scheme
    !
    ! turbulence scheme
    range_tkhmin     = 0.2_wp       ! minimum vertical diffusion (m**2/s) for heat/moisture
    range_tkmmin     = 0.2_wp       ! minimum vertical diffusion (m**2/s) for momentum
    range_tkred_sfc  = 4.0_wp       ! multiplicative change of reduction of minimum diffusion coefficients near the surface
    range_rlam_heat  = 8.0_wp       ! additive change of laminar transport resistance parameter
                                    ! (compensated by an inverse change of rat_sea)
    range_turlen     = 150._wp      ! turbulent length scale (m)
    range_a_hshr     = 1._wp        ! scaling factor for extended horizontal shear term
    range_a_stab     = 1._wp        ! scaling factor for stability correction in turbulence scheme
    range_c_diff     = 2._wp        ! length scale factor for vertical diffusion in turbulence scheme (multiplicative)
    range_q_crit     = 1._wp        ! critical value for normalized super-saturation in turbulent cloud scheme
    range_charnock   = 1.5_wp       ! multiplicative change of upper and lower bound of wind-speed dependent
                                    ! Charnock parameter
    !
    ! snow cover diagnosis
    range_minsnowfrac = 0.1_wp      ! Minimum value to which the snow cover fraction is artificially reduced
                                    ! in case of melting show (in case of idiag_snowfrac = 20)
    !
    ! TERRA
    range_c_soil      = 0.25_wp     ! evaporative surface area
    range_cwimax_ml   = 2._wp       ! capacity of interception storage (multiplicative perturbation)

    !
    ! LHN
    range_lhn_coef      = 0.0_wp    ! Scaling factor for latent heat nudging increments
    range_lhn_artif_fac = 0.0_wp    ! Scaling factor for artificial heating profile in latent heat nudging
    range_fac_lhn_down  = 0.0_wp    ! Lower limit for reduction of pre-existing latent heating in LHN
    range_fac_lhn_up    = 0.0_wp    ! Upper limit for increase of pre-existing latent heating in LHN

    ! external parameters specified depending on land-cover class
    ! all subsequent ranges indicate relative changes of the respective parameter
    range_z0_lcc   = 0.25_wp        ! Roughness length
    range_rootdp   = 0.2_wp         ! Root depth
    range_rsmin    = 0.2_wp         ! Minimum stomata resistance
    range_laimax   = 0.15_wp        ! Leaf area index

    use_ensemble_pert = .FALSE.     ! Usage of ensemble perturbations must be turned on explicitly
    itype_pert_gen    = 1           ! Type of ensemble perturbation generation:
                                    ! 1: equal distribution within perturbation range
                                    ! 2: use either default or extrema of perturbation range
    timedep_pert      = 0           ! 0: no time-dependence of ensemble perturbations
                                    ! 1: perturbations depend on start date, but remain fixed during forecast 
                                    ! 2: time-dependent perturbations varying sinusoidally within their range

    stdev_sst_pert    = 0._wp       ! No compensation for SST perturbations is applied by default

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('ensemble_pert_nml')
      READ(funit,NML=ensemble_pert_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('ensemble_pert_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, ensemble_pert_nml)   ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, ensemble_pert_nml)             ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, ensemble_pert_nml)    ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml



    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------



    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    config_range_gkwake       = range_gkwake
    config_range_gkdrag       = range_gkdrag
    config_range_gfrcrit      = range_gfrcrit
    config_range_gfluxlaun    = range_gfluxlaun
    config_range_zvz0i        = range_zvz0i
    config_range_rain_n0fac   = range_rain_n0fac
    config_range_entrorg      = range_entrorg
    config_range_rdepths      = range_rdepths
    config_range_rprcon       = range_rprcon
    config_range_capdcfac_et  = range_capdcfac_et
    config_range_capdcfac_tr  = range_capdcfac_tr
    config_range_lowcapefac   = range_lowcapefac
    config_range_negpblcape   = range_negpblcape
    config_range_rhebc        = range_rhebc
    config_range_texc         = range_texc
    config_range_qexc         = range_qexc
    config_range_minsnowfrac  = range_minsnowfrac
    config_range_c_soil       = range_c_soil
    config_range_cwimax_ml    = range_cwimax_ml
    config_range_box_liq      = range_box_liq
    config_range_thicklayfac  = range_thicklayfac
    config_range_box_liq_asy  = range_box_liq_asy
    config_range_fac_ccqc     = range_fac_ccqc
    config_range_tkhmin       = range_tkhmin
    config_range_tkmmin       = range_tkmmin
    config_range_tkred_sfc    = range_tkred_sfc
    config_range_rlam_heat    = range_rlam_heat
    config_range_turlen       = range_turlen
    config_range_a_hshr       = range_a_hshr
    config_range_a_stab       = range_a_stab
    config_range_c_diff       = range_c_diff
    config_range_q_crit       = range_q_crit
    config_range_charnock     = range_charnock
    config_range_lhn_coef     = range_lhn_coef
    config_range_lhn_artif_fac = range_lhn_artif_fac
    config_range_fac_lhn_down = range_fac_lhn_down
    config_range_fac_lhn_up   = range_fac_lhn_up
    config_range_z0_lcc       = range_z0_lcc
    config_range_rootdp       = range_rootdp
    config_range_rsmin        = range_rsmin
    config_range_laimax       = range_laimax
    config_stdev_sst_pert     = stdev_sst_pert
    config_use_ensemble_pert  = use_ensemble_pert
    config_itype_pert_gen     = itype_pert_gen
    config_timedep_pert       = timedep_pert

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=ensemble_pert_nml)                    
      CALL store_and_close_namelist(funit, 'ensemble_pert_nml')             
    ENDIF

    !--------------------------------------------------------
    ! 7. write the contents of the namelist to an ASCII file
    !--------------------------------------------------------
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=ensemble_pert_nml)


  END SUBROUTINE read_ensemble_pert_namelist

END MODULE mo_ensemble_pert_nml
