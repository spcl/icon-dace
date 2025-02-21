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

MODULE mo_aerosol_util

  USE mo_impl_constants,         ONLY: min_rlcell_int, iss, iorg, ibc, iso4, idu, nclass_aero
  USE mo_impl_constants_grf,     ONLY: grf_bdywidth_c
  USE mo_math_constants,         ONLY: rad2deg
  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish
  USE mo_loopindices,            ONLY: get_indices_c
  USE mo_lrtm_par,               ONLY: jpband => nbndlw
  USE mo_model_domain,           ONLY: t_patch
  USE mo_intp_data_strc,         ONLY: t_int_state
  USE mo_srtm_config,            ONLY: jpsw
  USE mo_lnd_nwp_config,         ONLY: ntiles_lnd, dzsoil, isub_water
  USE mo_nwp_tuning_config,      ONLY: tune_dust_abs
  USE mo_advection_config,       ONLY: advection_config
  USE mo_nwp_phy_state,          ONLY: phy_params
  USE mo_aerosol_sources_types,  ONLY: p_dust_source_const
  USE mo_aerosol_sources,        ONLY: aerosol_dust_aod_source, aerosol_ssa_aod_source, &
                                   &   calc_anthro_aod, calc_so4_nucleation
  USE mo_math_laplace,           ONLY: nabla2_scalar
  USE mo_util_phys,              ONLY: rel_hum
  USE mtime,                     ONLY: datetime
  USE mo_fortran_tools,          ONLY: set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_aerosol_util'

  !RRTM
  REAL  (wp) ::             &
  zaea_rrtm(jpsw+jpband,5), &  ! ratio of optical thickness for the absorption in spectral
                               ! interval jpspec  and total optical thickness at 0.55m*1.E-06 
                               ! for an aerosoltyp specified by second array index
  zaes_rrtm(jpsw+jpband,5), &  ! analog for the optical thickness of scattering 
  zaeg_rrtm(jpsw+jpband,5)!, zaef_rrtm(jpsw+jpband,5)

  ! The (long wave) wavenumbers from rrtm are not set when not using RRTM
  ! So we copy the values from mo_lrtm_setup (long wave) and mo_srtm_config (short wave)
  REAL(wp), PARAMETER :: & !< Lower wavenumber bound long wave
    &  rrtm_wavenum1_lw(jpband) = (/  10._wp, 350._wp, 500._wp, 630._wp, 700._wp, 820._wp, &
                                &    980._wp,1080._wp,1180._wp,1390._wp,1480._wp,1800._wp, &
                                &   2080._wp,2250._wp,2380._wp,2600._wp/)
  REAL(wp), PARAMETER :: & !< Upper wavenumber bound long wave
    &  rrtm_wavenum2_lw(jpband) = (/ 350._wp, 500._wp, 630._wp, 700._wp, 820._wp, 980._wp, &
                                &   1080._wp,1180._wp,1390._wp,1480._wp,1800._wp,2080._wp, &
                                &   2250._wp,2380._wp,2600._wp,3250._wp/)
  REAL(wp), PARAMETER :: & !< Lower wavenumber bound short wave
    &  rrtm_wavenum1_sw(jpsw) = (/ 2600._wp, 3250._wp, 4000._wp, 4650._wp, 5150._wp, 6150._wp, &
                              &    7700._wp, 8050._wp,12850._wp,16000._wp,22650._wp,29000._wp, &
                              &   38000._wp, 820._wp  /)
  REAL(wp), PARAMETER :: & !< Upper wavenumber bound short wave
    &  rrtm_wavenum2_sw(jpsw) = (/ 3250._wp, 4000._wp, 4650._wp, 5150._wp, 6150._wp, 7700._wp, &
                              &    8050._wp,12850._wp,16000._wp,22650._wp,29000._wp,38000._wp, &
                              &   50000._wp, 2600._wp /)

  !ecRad
  TYPE t_tegen_scal_factors
    ! Total number of wavelength bands
    INTEGER :: n_bands
    ! Scaling factors from 550nm to wavelengths bands
    REAL(wp), ALLOCATABLE :: &
      &  absorption(:,:),            & !< Dim [n_bands, nspecies=5]
      &  scattering(:,:),            & !< Dim [n_bands, nspecies=5] 
      &  asymmetry(:,:)                !< Dim [n_bands, nspecies=5]
    CONTAINS
      PROCEDURE :: init     => init_tegen_scal_factors
      PROCEDURE :: finalize => finalize_tegen_scal_factors
  END TYPE t_tegen_scal_factors

  TYPE(t_tegen_scal_factors), TARGET :: &
    &  tegen_scal_factors_rrtm,         & !< Scaling factors original rrtm bands
    &  tegen_scal_factors_mod,          & !< Modified scaling factors
    &  tegen_scal_factors                 !< Used scaling factors

  PUBLIC :: zaea_rrtm, zaes_rrtm, zaeg_rrtm
  PUBLIC :: aerdis
  PUBLIC :: init_aerosol_props_tegen_rrtm, tune_dust
  PUBLIC :: prog_aerosol_2D, aerosol_2D_diffusion
  PUBLIC :: tegen_scal_factors

  !$ACC DECLARE CREATE(zaea_rrtm, zaes_rrtm, zaeg_rrtm)

CONTAINS

  !!  Subroutine aerdis is simplified version from COSMO model (version 4.16).
  !!
   
  SUBROUTINE aerdis ( klevp1, kbdim, jcs, jce, petah,  pvdaes, pvdael, pvdaeu, pvdaed, lacc )
    
    !------------------------------------------------------------------------------
    !
    ! Description:
    !
    ! The module procedure aerdis provides parameters for the vertical distribution
    ! of aerosols (based on the original code of J.F. Geleyn (ECMWF, 4.11.82).
    !
    ! The routine computes the values PVDAE* (* = s, l, u or d for sea, land
    ! urban or desert) of a surfach-normalised vertical distribution of aerosols'
    ! optical depth from the argument petah (vertical coordinate) at klevp1 levels.
    ! It also sets values for non-geograpically weighted total optical depths (at
    ! 55 micrometer wavelength) paeopn for the same four types and similar optical
    ! depths diveded by pressure for bachground well-mixed aerosols of three types
    ! p**bga (** = tr, vo or st for tropospheric, volcanic (stratosperic ashes) or
    ! stratosperic (sulfuric type)). It finally sets values for the power to be
    ! applied to a temperature ratio smaller than two in order to obtain an index
    ! one in the stratosphere and zero in the troposphere with a relatively smooth
    ! transistion (ptrpt), as well as for adsorption coefficients fo water to the
    ! three type of troposperic aerosols (paeadk) with a minimum value ( in the 
    ! whole atmosphere) for the sum of the products paeadk by the optical depths
    ! divided by pressure thickness: paeadm. 
    !
    ! Method:
    !
    ! Straightforward, equivalent heights are given in meters (8434 for the
    ! atmosphere) and tropospheric and stratospheric pressure boundary values
    ! are set at 101325 and 19330 Pascal. 
    !
    !------------------------------------------------------------------------------
    
    ! Subroutine arguments:
    ! --------------------
    
    ! Input data
    ! ----------
    INTEGER, INTENT (IN) ::  &
      & klevp1,         &           ! number of model layer interfaces
      & kbdim,          &
      & jcs,            &
      & jce

    REAL    (wp), INTENT (IN) ::  &
      petah(kbdim,klevp1)    ! normalized vertical coordinate at half levels

    LOGICAL, OPTIONAL, INTENT(IN) :: lacc

    ! Output data
    ! -----------
    REAL    (wp), INTENT (OUT) ::  &
      pvdaes(kbdim,klevp1), & ! normalized vertical distribution (sea)
      pvdael(kbdim,klevp1), & ! normalized vertical distribution (land)
      pvdaeu(kbdim,klevp1), & ! normalized vertical distribution (urban)
      pvdaed(kbdim,klevp1)    ! normalized vertical distrubution (desert)

    ! Local parameters:
    ! -------------
    REAL (wp), PARAMETER  ::  &
      zhss = 8434.0_wp/1000.0_wp ,  & !
      zhsd = 8434.0_wp/3000.0_wp      !

    INTEGER :: jc,jk
    REAL(wp) :: log_eta
    LOGICAL :: lzacc

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    ! Begin Subroutine aerdis              
    !------------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)
    !$ACC DATA PRESENT(petah, pvdaes, pvdael, pvdaeu, pvdaed) IF(lzacc)

    ! default data present
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO jc=jcs,jce
      pvdaes(jc,1) = 0.0_wp
      pvdael(jc,1) = 0.0_wp
      pvdaeu(jc,1) = 0.0_wp
      pvdaed(jc,1) = 0.0_wp
    ENDDO
    !$ACC END PARALLEL

!!$  IF(petah(1).NE.0._wp) THEN
!!$     pvdaes(1) = petah(1)**zhss
!!$     pvdael(1) = petah(1)**zhsl
!!$     pvdaeu(1) = petah(1)**zhsu
!!$     pvdaed(1) = petah(1)**zhsd
!!$  END IF

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(log_eta)
    DO jk=2,klevp1
      DO jc=jcs,jce
        log_eta       = LOG(petah(jc,jk))
        pvdaes(jc,jk) = EXP(zhss*log_eta) ! petah(jc,jk)**zhss
        pvdael(jc,jk) = pvdaes(jc,jk)     ! petah(jc,jk)**zhsl; zhsl is the same as zhss
        pvdaeu(jc,jk) = pvdaes(jc,jk)     ! petah(jc,jk)**zhsu; zhsu is the same as zhss
        pvdaed(jc,jk) = EXP(zhsd*log_eta) ! petah(jc,jk)**zhsd
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC END DATA

  END SUBROUTINE aerdis

  
  SUBROUTINE init_aerosol_props_tegen_rrtm

  ! the following aerosol types (second array index) are considered:
  ! 1. continental, 2. maritime, 3. desert, 4. urban, 5. stratospheric background (SB)

   !absorption
   zaea_rrtm=RESHAPE( (/ &
     &0.0304_wp,0.0367_wp,0.0462_wp,0.0566_wp,0.0496_wp,0.0336_wp,0.0355_wp,0.0456_wp,&
     &0.0272_wp,0.0264_wp,0.0290_wp,0.0156_wp,0.0165_wp,0.0157_wp,0.0138_wp,0.0401_wp,&
     &0.0401_wp,0.0760_wp,0.0214_wp,0.0227_wp,0.0295_wp,0.0394_wp,0.0431_wp,0.0519_wp,&
     &0.0611_wp,0.0774_wp,0.1012_wp,0.1412_wp,0.2632_wp,0.0324_wp,                    &
     &0.1096_wp,0.1614_wp,0.2294_wp,0.2506_wp,0.2242_wp,0.1190_wp,0.0680_wp,0.0664_wp,&
     &0.0656_wp,0.0749_wp,0.1250_wp,0.0425_wp,0.0498_wp,0.0425_wp,0.0259_wp,0.1619_wp,&
     &0.1619_wp,0.2152_wp,0.0139_wp,0.0119_wp,0.0046_wp,0.0036_wp,0.0020_wp,0.0016_wp,&
     &0.0012_wp,0.0013_wp,0.0016_wp,0.0035_wp,0.0147_wp,0.0882_wp,                    &
     &0.0974_wp,0.1529_wp,0.1643_wp,0.1373_wp,0.1753_wp,0.1923_wp,0.2804_wp,0.2426_wp,&
     &0.1263_wp,0.1321_wp,0.0979_wp,0.0664_wp,0.0360_wp,0.0311_wp,0.0325_wp,0.0833_wp,&
     &0.0833_wp,0.1170_wp,0.0739_wp,0.0631_wp,0.0604_wp,0.0628_wp,0.0645_wp,0.0677_wp,&
     &0.0843_wp,0.1328_wp,0.2224_wp,0.3022_wp,0.3579_wp,0.1820_wp,                    &
     &0.0267_wp,0.0329_wp,0.0420_wp,0.0515_wp,0.0461_wp,0.0332_wp,0.0354_wp,0.0447_wp,&
     &0.0303_wp,0.0306_wp,0.0342_wp,0.0248_wp,0.0274_wp,0.0276_wp,0.0271_wp,0.0526_wp,&
     &0.0526_wp,0.0903_wp,0.0450_wp,0.0492_wp,0.0596_wp,0.0754_wp,0.0842_wp,0.1082_wp,&
     &0.1429_wp,0.1926_wp,0.2595_wp,0.3379_wp,0.4761_wp,0.0340_wp,                    &
     &0.0060_wp,0.0117_wp,0.0269_wp,0.0222_wp,0.0195_wp,0.0398_wp,0.0733_wp,0.1091_wp,&     ! SB
     &0.1124_wp,0.0415_wp,0.0424_wp,0.0495_wp,0.0451_wp,0.0484_wp,0.0540_wp,0.0735_wp,&     ! SB
     &0.0735_wp,0.0188_wp,0.0021_wp,0.0014_wp,0.0007_wp,0.0002_wp,0.0000_wp,0.0000_wp,&     ! SB
     &0.0000_wp,0.0000_wp,0.0000_wp,0.0000_wp,0.0000_wp,0.0628_wp/),(/jpsw+jpband,5/))      ! SB

   !scattering
   zaes_rrtm=RESHAPE( (/ &
     &0.0060_wp,0.0107_wp,0.0134_wp,0.0150_wp,0.0152_wp,0.0200_wp,0.0232_wp,0.0211_wp,&
     &0.0112_wp,0.0186_wp,0.0128_wp,0.0260_wp,0.0339_wp,0.0368_wp,0.0409_wp,0.0527_wp,&
     &0.0527_wp,0.0621_wp,0.0715_wp,0.0929_wp,0.1276_wp,0.1895_wp,0.2350_wp,0.3930_wp,&
     &0.6641_wp,0.9834_wp,1.3737_wp,1.7160_wp,1.9115_wp,0.0198_wp,                    &
     &0.0188_wp,0.0421_wp,0.0576_wp,0.0547_wp,0.0430_wp,0.0367_wp,0.0806_wp,0.1209_wp,&
     &0.1681_wp,0.2257_wp,0.2440_wp,0.3622_wp,0.4540_wp,0.5026_wp,0.5765_wp,0.5986_wp,&
     &0.5986_wp,0.5225_wp,0.7420_wp,0.8311_wp,0.8970_wp,0.9444_wp,0.9637_wp,0.9763_wp,&
     &0.9855_wp,1.0034_wp,1.0337_wp,1.0640_wp,1.0795_wp,0.1312_wp,                    &
     &0.0458_wp,0.0823_wp,0.0667_wp,0.0642_wp,0.1080_wp,0.1471_wp,0.2422_wp,0.1216_wp,&
     &0.0717_wp,0.1616_wp,0.2027_wp,0.3042_wp,0.4045_wp,0.4369_wp,0.4685_wp,0.5043_wp,&
     &0.5043_wp,0.5782_wp,0.6898_wp,0.7477_wp,0.7926_wp,0.8320_wp,0.8503_wp,0.8736_wp,&
     &0.8874_wp,0.8737_wp,0.8278_wp,0.7857_wp,0.7571_wp,0.1714_wp,                    &
     &0.0048_wp,0.0085_wp,0.0107_wp,0.0119_wp,0.0121_wp,0.0160_wp,0.0185_wp,0.0170_wp,&
     &0.0090_wp,0.0150_wp,0.0103_wp,0.0210_wp,0.0274_wp,0.0298_wp,0.0332_wp,0.0430_wp,&
     &0.0430_wp,0.0485_wp,0.0593_wp,0.0776_wp,0.1073_wp,0.1610_wp,0.2008_wp,0.3398_wp,&
     &0.5809_wp,0.8701_wp,1.2309_wp,1.5535_wp,1.7368_wp,0.0159_wp,                    &
     &0.0000_wp,0.0000_wp,0.0000_wp,0.0000_wp,0.0001_wp,0.0003_wp,0.0006_wp,0.0008_wp,&     ! SB
     &0.0005_wp,0.0003_wp,0.0008_wp,0.0013_wp,0.0024_wp,0.0030_wp,0.0040_wp,0.0059_wp,&     ! SB
     &0.0059_wp,0.0123_wp,0.0236_wp,0.0384_wp,0.0651_wp,0.1246_wp,0.1801_wp,0.3807_wp,&     ! SB
     &0.7105_wp,1.0514_wp,1.3754_wp,1.5334_wp,1.5495_wp,0.0009_wp/),(/jpsw+jpband,5/))      ! SB

   !asymmetry factor
   zaeg_rrtm=RESHAPE( (/ &
     &0.4388_wp,0.5396_wp,0.6191_wp,0.6535_wp,0.6876_wp,0.6718_wp,0.6493_wp,0.6782_wp,&
     &0.7958_wp,0.7537_wp,0.7757_wp,0.7821_wp,0.7583_wp,0.7487_wp,0.7351_wp,0.6917_wp,&
     &0.6917_wp,0.6989_wp,0.6982_wp,0.6726_wp,0.6426_wp,0.6294_wp,0.6337_wp,0.6582_wp,&
     &0.6850_wp,0.7061_wp,0.7212_wp,0.7306_wp,0.7417_wp,0.6978_wp,                    &
     &0.4062_wp,0.4507_wp,0.4878_wp,0.5302_wp,0.5850_wp,0.6962_wp,0.7242_wp,0.7293_wp,&
     &0.7414_wp,0.7484_wp,0.7607_wp,0.7785_wp,0.7805_wp,0.7785_wp,0.7724_wp,0.7690_wp,&
     &0.7690_wp,0.8348_wp,0.8316_wp,0.8170_wp,0.8074_wp,0.7990_wp,0.7954_wp,0.7897_wp,&
     &0.7884_wp,0.7927_wp,0.8001_wp,0.8057_wp,0.8076_wp,0.7462_wp,                    &
     &0.4219_wp,0.3928_wp,0.5306_wp,0.6229_wp,0.5544_wp,0.5454_wp,0.4353_wp,0.5736_wp,&
     &0.7502_wp,0.6957_wp,0.7038_wp,0.6881_wp,0.6740_wp,0.6739_wp,0.6784_wp,0.6969_wp,&
     &0.6969_wp,0.7068_wp,0.6965_wp,0.6918_wp,0.6904_wp,0.6911_wp,0.6915_wp,0.6952_wp,&
     &0.7080_wp,0.7326_wp,0.7689_wp,0.8000_wp,0.8206_wp,0.5788_wp,                    &
     &0.4387_wp,0.5394_wp,0.6187_wp,0.6531_wp,0.6871_wp,0.6712_wp,0.6482_wp,0.6756_wp,&
     &0.7930_wp,0.7498_wp,0.7685_wp,0.7766_wp,0.7520_wp,0.7419_wp,0.7277_wp,0.6828_wp,&
     &0.6828_wp,0.6875_wp,0.6872_wp,0.6622_wp,0.6333_wp,0.6209_wp,0.6250_wp,0.6479_wp,&
     &0.6725_wp,0.6912_wp,0.7043_wp,0.7129_wp,0.7254_wp,0.6956_wp,                    &
     &0.0021_wp,0.0039_wp,0.0061_wp,0.0078_wp,0.0109_wp,0.0161_wp,0.0201_wp,0.0206_wp,&     ! SB
     &0.0217_wp,0.0320_wp,0.0428_wp,0.0583_wp,0.0773_wp,0.0856_wp,0.0985_wp,0.1310_wp,&     ! SB
     &0.1310_wp,0.1906_wp,0.2625_wp,0.3154_wp,0.3869_wp,0.4787_wp,0.5279_wp,0.6272_wp,&     ! SB
     &0.6941_wp,0.7286_wp,0.7358_wp,0.7177_wp,0.6955_wp,0.0616_wp/),(/jpsw+jpband,5/))      ! SB

    !$ACC UPDATE DEVICE(zaea_rrtm, zaes_rrtm, zaeg_rrtm) ASYNC(1)

  END SUBROUTINE init_aerosol_props_tegen_rrtm


  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE init_tegen_scal_factors
  !! Constructor for t_tegen_scal_factors
  !!
  !---------------------------------------------------------------------------------------
  SUBROUTINE init_tegen_scal_factors(this, n_bands)
    CLASS(t_tegen_scal_factors), INTENT(inout) :: &
      &  this      !< Scaling factor information
    INTEGER, INTENT(in) :: &
      &  n_bands   !< Total number of wavelength bands

    this%n_bands = n_bands

    ALLOCATE(this%absorption(this%n_bands,5))
    ALLOCATE(this%scattering(this%n_bands,5))
    ALLOCATE(this%asymmetry (this%n_bands,5))
    !$ACC ENTER DATA CREATE(this%absorption, this%scattering, this%asymmetry)
  END SUBROUTINE
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE finalize_tegen_scal_factors
  !! Destructor for t_tegen_scal_factors
  !!
  !---------------------------------------------------------------------------------------
  SUBROUTINE finalize_tegen_scal_factors(this)
    CLASS(t_tegen_scal_factors), INTENT(inout) :: &
      &  this      !< Scaling factor information

    this%n_bands = 0

    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(this%absorption) IF(ALLOCATED(this%absorption))
    IF (ALLOCATED(this%absorption)) &
      DEALLOCATE(this%absorption)
    !$ACC EXIT DATA DELETE(this%scattering) IF(ALLOCATED(this%scattering))
    IF (ALLOCATED(this%scattering)) &
      DEALLOCATE(this%scattering)
    !$ACC EXIT DATA DELETE(this%asymmetry) IF(ALLOCATED(this%asymmetry))
    IF (ALLOCATED(this%asymmetry)) &
      DEALLOCATE(this%asymmetry)
  END SUBROUTINE
  !---------------------------------------------------------------------------------------


  ! Very simple parameterization of source and sink terms for prognostic 2D aerosol fields
  !
  SUBROUTINE prog_aerosol_2D (jcs, jce, jg, nproma, nlev, dtime, iprog_aero, aerosol, &
    &                         aercl_ss,aercl_or,aercl_bc,aercl_su,aercl_du,           &
    &                         exner, temp, qv, cosmu0, rr_gsp,sr_gsp,rr_con,sr_con,   &
    &                         soiltype,plcov_t,frac_t,w_so_t, w_so_ice_t, h_snow_t,   &
    &                         t_seasfc, lc_class_t, rho, tcm_t, u, v, sp_10m, emi_bc, &
    &                         emi_oc, emi_so2, bcfire, ocfire, so2fire,               &
    &                         idx_lst_t, gp_count_t , i_count_sea, idx_sea)
    REAL(wp), INTENT(in)            :: &
      &  dtime,                        & !< Time step (s)
      &  aercl_ss(:), aercl_or(:),     & !< AOD climatology (sea salt, organic)
      &  aercl_bc(:), aercl_su(:),     & !< AOD climatology (black carbon, sulfate)
      &  aercl_du(:),                  & !< AOD climatology (dust)
      &  exner(:,:),                   & !< Exner pressure
      &  temp(:,:), qv(:,:),           & !< Air temperature, specific humidity
      &  cosmu0(:),                    & !< Cosine of solar zenith angle
      &  rr_gsp(:),sr_gsp(:),          & !< Grid-scale rain & snow rate
      &  rr_con(:),sr_con(:),          & !< Convective rain & snow rate
      &  plcov_t(:,:),                 & !< Plant cover (tiled)
      &  frac_t(:,:),                  & !< Tile fraction
      &  w_so_t(:,:), w_so_ice_t(:,:), & !< Soil water & ice (tiled)
      &  h_snow_t(:,:),                & !< Snow height (tiled)
      &  t_seasfc(:),                  & !< Sea surface temperature
      &  rho(:),                       & !< Air density
      &  tcm_t(:,:),                   & !< Transfer coefficient for momentum
      &  u(:), v(:),                   & !< Wind vector components
      &  sp_10m(:),                    & !< Wind speed in 10m
      &  emi_bc(:),                    & !< Precursor for anthropogenic emissions (black carbon)
      &  emi_oc(:),                    & !< Precursor for anthropogenic emissions (organic carbon)
      &  emi_so2(:),                   & !< Precursor for anthropogenic emissions (SO2)
      &  bcfire(:),                    & !< Precursor for wildfire emissions (black carbon)
      &  ocfire(:),                    & !< Precursor for wildfire emissions (organic carbon)
      &  so2fire(:)                      !< Precursor for wildfire emissions (SO2)
    INTEGER,  INTENT(in) :: &
      &  nproma, nlev,      & !< Array dimensions
      &  jcs, jce,          & !< Start and end index of nproma loop
      &  jg,                & !< Domain index
      &  iprog_aero,        & !< Prognostic aerosol mode: 1 only dust, 2 all
      &  soiltype(:),       & !< Soil type index (dim: nproma)
      &  lc_class_t(:,:),   & !< Land use class index (dim: nproma, ntiles)
      &  idx_lst_t(:,:),    & !< Tiled index list to loop over land points (dim: nproma,ntiles)
      &  gp_count_t(:),     & !< Returns number of local grid points per tile (dim: ntiles)
      &  i_count_sea,       & !< Number of open water points in current block
      &  idx_sea(:)           !< Indices of open water points in current block
    REAL(wp), INTENT(inout) :: &
      &  aerosol(:,:)         !< Aerosol Optical Depth (AOD)
    ! Local variables
    REAL(wp) ::                 &
      &  relhum(nproma,nlev)      !< Relative humidity (0 - 1)
    REAL(wp) ::                           &
      &  relax_bc,  relax_oc,             & !< Relaxation time scales black & organic carbon
      &  relax_so4, relax_du,             & !< Relaxation time scales sulphate & dust
      &  relax_ss,                        & !< Relaxation time scale sea salt
      &  minfrac,                         & !< minimum allowed fraction of climatological AOD
      &  washout, washout_scale,          & !< Washout and washout scale for dust
      &  aod_flux,                        & !< Source function for aerosol optical depth
      &  tunefac_bc_ant, tunefac_org_ant, & !< Conversion factor anthr. bc/oc emission to bc/oc AOD emission
      &  tunefac_so4_ant, tunefac_bc_wf,  & !< Conversion factor anthr. so2/wildfire bc emission to so4/bc AOD emission
      &  tunefac_org_wf,  tunefac_so4_wf    !< Conversion factor wildfire oc/so2 emission to oc/so4 AOD emission
    INTEGER ::              &
      &  jc, jt, jcl, jk,   & !< Loop indices
      &  i_count_lnd          !< Number of land grid points in current block

    relax_ss       = 1._wp/(3._wp*86400._wp)  ! 3 days
    relax_du       = 1._wp/(8._wp*86400._wp)  ! 8 days
    relax_so4      = 1._wp/(7._wp*86400._wp)  ! 7 days
    relax_bc       = 1._wp/(10._wp*86400._wp) ! 10 days
    relax_oc       = 1._wp/(14._wp*86400._wp) ! 14 days
    washout_scale  = 1._wp/5._wp             ! e-folding scale 7.5 mm WE precipitation
    minfrac        = 0.025_wp
    tunefac_bc_ant = 5.e3_wp
    tunefac_org_ant= 2.e4_wp
    tunefac_so4_ant= 3.e3_wp
    tunefac_bc_wf  = 3.e4_wp
    tunefac_org_wf = 2.e4_wp
    tunefac_so4_wf = 5.e3_wp

    ! Prediction of mineral dust; other aerosol classes are treated prognostically only if iprog_aero=2

    ! Relaxation to scaled climatology
    DO jc = jcs, jce
      aerosol(jc,idu)  = aerosol(jc,idu)  + dtime*relax_du*(aercl_du(jc)-aerosol(jc,idu))
    ENDDO

    DO jt = 1, ntiles_lnd
      i_count_lnd = gp_count_t(jt)
      IF (i_count_lnd == 0) CYCLE ! skip loop if the index list for the given tile is empty
!$NEC ivdep
      DO jcl = 1, i_count_lnd
        jc = idx_lst_t(jcl,jt)
        CALL aerosol_dust_aod_source (p_dust_source_const(jg), dzsoil(1), w_so_t(jc,jt), h_snow_t(jc,jt), &
          &                           w_so_ice_t(jc,jt), soiltype(jc), plcov_t(jc,jt), lc_class_t(jc,jt), &
          &                           rho(jc), tcm_t(jc,jt), u(jc), v(jc), aod_flux)
        ! Update AOD field with tendency from aod_flux
        aerosol(jc,idu) = aerosol(jc,idu) + aod_flux * frac_t(jc,jt) * dtime
      ENDDO ! jcl
    ENDDO !jt

    DO jc = jcs, jce
      ! Washout using scale-dependent convective area fraction rcucov from convection param.
      washout = dtime*washout_scale*(rr_gsp(jc)+sr_gsp(jc)+phy_params(jg)%rcucov*(rr_con(jc)+sr_con(jc)))*aerosol(jc,idu)
      aerosol(jc,idu)  = aerosol(jc,idu) - washout
      ! Ensure that the aerosol optical depth does not fall below 2.5% of the climatological value
      aerosol(jc,idu)  = MAX(aerosol(jc,idu),  minfrac*aercl_du(jc))
    ENDDO

    IF (iprog_aero >= 2) THEN

      ! Calculate relative humidity for nucleation
      DO jk = advection_config(jg)%kstart_aero(1), advection_config(jg)%kend_aero(1)
        DO jc = jcs, jce
          relhum(jc,jk)    = rel_hum(temp(jc,jk), qv(jc,jk), exner(jc,jk)) /100._wp
        ENDDO
      ENDDO

      DO jc = jcs, jce
        ! Relaxation to scaled climatology
        aerosol(jc,iss)  = aerosol(jc,iss)  + dtime*relax_ss *(aercl_ss(jc)-aerosol(jc,iss))
        aerosol(jc,iorg) = aerosol(jc,iorg) + dtime*relax_oc *(aercl_or(jc)-aerosol(jc,iorg))
        aerosol(jc,ibc)  = aerosol(jc,ibc)  + dtime*relax_bc *(aercl_bc(jc)-aerosol(jc,ibc))
        aerosol(jc,iso4) = aerosol(jc,iso4) + dtime*relax_so4*(aercl_su(jc)-aerosol(jc,iso4))
        ! Sources based on anthropogenic emission datasets
        aerosol(jc,ibc)  = aerosol(jc,ibc)  + dtime * calc_anthro_aod( emi_bc(jc),  tunefac_bc_ant )
        aerosol(jc,iorg) = aerosol(jc,iorg) + dtime * calc_anthro_aod( emi_oc(jc),  tunefac_org_ant )
        aerosol(jc,iso4) = aerosol(jc,iso4) + dtime * calc_anthro_aod( emi_so2(jc), tunefac_so4_ant )
      ENDDO

      IF (iprog_aero > 2) THEN
        ! Sources based on wildfire emission datasets
        DO jc = jcs, jce
          aerosol(jc,ibc)  = aerosol(jc,ibc)  + dtime * calc_anthro_aod( bcfire(jc),  tunefac_bc_wf )
          aerosol(jc,iorg) = aerosol(jc,iorg) + dtime * calc_anthro_aod( ocfire(jc),  tunefac_org_wf )
          aerosol(jc,iso4) = aerosol(jc,iso4) + dtime * calc_anthro_aod( so2fire(jc), tunefac_so4_wf )
        ENDDO
      ENDIF

      ! Source based on nucleation
      CALL calc_so4_nucleation(jcs, jce, advection_config(jg)%kstart_aero(1), advection_config(jg)%kend_aero(1), &
        &                      temp(:,:), relhum(:,:), cosmu0(:), aerosol(:,iso4))

      ! Sea salt aerosol source
!$NEC ivdep
      DO jcl = 1, i_count_sea
        jc = idx_sea(jcl)
        CALL aerosol_ssa_aod_source (t_seasfc(jc), sp_10m(jc), aod_flux)
        aerosol(jc,iss) = aerosol(jc,iss) + aod_flux * frac_t(jc,isub_water) * dtime
      ENDDO

      DO jc = jcs, jce
        ! Washout using scale-dependent convective area fraction rcucov from convection param.
        washout = dtime*washout_scale*(rr_gsp(jc)+sr_gsp(jc)+phy_params(jg)%rcucov*(rr_con(jc)+sr_con(jc)))*aerosol(jc,iss)
        aerosol(jc,iss)   = aerosol(jc,iss)  - washout
        washout = dtime*washout_scale*(rr_gsp(jc)+sr_gsp(jc)+phy_params(jg)%rcucov*(rr_con(jc)+sr_con(jc)))*aerosol(jc,ibc)
        aerosol(jc,ibc)   = aerosol(jc,ibc)  - washout
        washout = dtime*washout_scale*(rr_gsp(jc)+sr_gsp(jc)+phy_params(jg)%rcucov*(rr_con(jc)+sr_con(jc)))*aerosol(jc,iorg)
        aerosol(jc,iorg)  = aerosol(jc,iorg) - washout
        washout = dtime*washout_scale*(rr_gsp(jc)+sr_gsp(jc)+phy_params(jg)%rcucov*(rr_con(jc)+sr_con(jc)))*aerosol(jc,iso4)
        aerosol(jc,iso4)  = aerosol(jc,iso4) - washout
        ! Ensure that the aerosol optical depth does not fall below 2.5% of the climatological value
        aerosol(jc,iss)  = MAX(aerosol(jc,iss),  minfrac*aercl_ss(jc))
        aerosol(jc,iorg) = MAX(aerosol(jc,iorg), minfrac*aercl_or(jc))
        aerosol(jc,ibc)  = MAX(aerosol(jc,ibc),  minfrac*aercl_bc(jc))
        aerosol(jc,iso4) = MAX(aerosol(jc,iso4), minfrac*aercl_su(jc))
      ENDDO

    ENDIF

  END SUBROUTINE prog_aerosol_2D


  ! Tuning of longwave absorption coefficient of mineral dust in order to reduce cold bias in the Saharan region
  !
  SUBROUTINE tune_dust (lat,lon,iend,tunefac)

    REAL(wp), INTENT(in) :: lat(:), lon(:)
    INTEGER,  INTENT(in) :: iend

    REAL(wp), INTENT(out) :: tunefac(:,:)

    INTEGER :: jc, jb
    REAL(wp) :: maxfac

    DO jb = 1, jpband
      maxfac = tune_dust_abs*5._wp*(jpband-MAX(8,jb))/REAL(jpband-8,wp)
      DO jc = 1, iend
        tunefac(jc,jb) = 1._wp + maxfac*(1._wp - MIN(1._wp,((rad2deg*lat(jc)-15._wp)/20._wp)**4)) * &
         (1._wp - MIN(1._wp,((rad2deg*lon(jc)-20._wp)/50._wp)**4))
      ENDDO
    ENDDO


  END SUBROUTINE tune_dust


  SUBROUTINE aerosol_2D_diffusion( p_patch, p_int_state, nproma, aerosol )
    TYPE(t_patch), INTENT(in)     :: &
      &  p_patch                       !< Current patch
    TYPE(t_int_state), INTENT(in) :: &
      &  p_int_state                   !< interpolation state
    INTEGER,  INTENT(in)          :: &
      &  nproma
    REAL(wp), INTENT(inout)       :: &
      &  aerosol(:,:,:)                !< Aerosol container
    ! Local variables
    REAL(wp)                      :: &
      &  diff_coeff, diff_coeff_so4, & !< Diffusion coefficients general and so4)
      &  diff_coeff_dust,            & !< Diffusion coefficient dust
      &  nabla2_aero(nproma,nclass_aero,p_patch%nblks_c) !< Laplacian of aerosol(:,:,:)
    INTEGER                       :: &
      &  jb, jc,                     &
      &  i_rlstart, i_rlend,         &
      &  i_startblk, i_endblk,       & 
      &  i_startidx, i_endidx

    diff_coeff      = 0.1_wp
    diff_coeff_so4  = 0.05_wp
    diff_coeff_dust = 0.05_wp

    CALL nabla2_scalar(aerosol(:,:,:),          &
      &                p_patch, p_int_state,    &
      &                nabla2_aero(:,:,:),      &
      &                iss, idu, grf_bdywidth_c+1, min_rlcell_int)

    i_rlstart  = grf_bdywidth_c+1
    i_rlend    = min_rlcell_int
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jc = i_startidx, i_endidx
        aerosol(jc,idu,jb)  = MAX(0.0_wp, aerosol(jc,idu,jb)  + diff_coeff_dust*       &
                                  p_patch%cells%area(jc,jb) * nabla2_aero(jc,idu,jb))
        aerosol(jc,iss,jb)  = MAX(0.0_wp, aerosol(jc,iss,jb)  + diff_coeff     *       &
                                  p_patch%cells%area(jc,jb) * nabla2_aero(jc,iss,jb))
        aerosol(jc,iorg,jb) = MAX(0.0_wp, aerosol(jc,iorg,jb) + diff_coeff     *       &
                                  p_patch%cells%area(jc,jb) * nabla2_aero(jc,iorg,jb))
        aerosol(jc,ibc,jb)  = MAX(0.0_wp, aerosol(jc,ibc,jb)  + diff_coeff     *       &
                                  p_patch%cells%area(jc,jb) * nabla2_aero(jc,ibc,jb))
        aerosol(jc,iso4,jb) = MAX(0.0_wp, aerosol(jc,iso4,jb) + diff_coeff_so4 *       &
                                  p_patch%cells%area(jc,jb) * nabla2_aero(jc,iso4,jb))
      ENDDO !jc
    ENDDO !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE aerosol_2D_diffusion

END MODULE mo_aerosol_util

