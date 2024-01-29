! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Global Data Module File
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
! SPDX-License-Identifier: GPL-3.0-only  
! ---------------------------------------------------------------

MODULE messy_mecca_kpp_global 

  USE messy_cmn_photol_mem ! IP_MAX, ip_*, jname
  USE messy_mecca_kpp_parameters, ONLY:  NSPEC, NVAR, NFIX, NREACT 
  USE mo_kind,                    ONLY:  dp
  PUBLIC
  SAVE


! Declaration of global variables

! C - Concentration of all species
  REAL(kind=dp) :: C(NSPEC)
! VAR - Concentrations of variable species (global)
  REAL(kind=dp) :: VAR(NVAR)
! FIX - Concentrations of fixed species (global)
  REAL(kind=dp) :: FIX(NFIX)
! VAR, FIX are chunks of array C
      EQUIVALENCE( C(1),VAR(1) )
      EQUIVALENCE( C(13),FIX(1) )
! RCONST - Rate constants (global)
  REAL(kind=dp) :: RCONST(NREACT)
! TIME - Current integration time
  REAL(kind=dp) :: TIME
! SUN - Sunlight intensity between [0,1]
  REAL(kind=dp) :: SUN
! TEMP - Temperature
  REAL(kind=dp) :: TEMP
! TSTART - Integration start time
  REAL(kind=dp) :: TSTART
! TEND - Integration end time
  REAL(kind=dp) :: TEND
! DT - Integration step
  REAL(kind=dp) :: DT
! ATOL - Absolute tolerance
  REAL(kind=dp) :: ATOL(NVAR)
! RTOL - Relative tolerance
  REAL(kind=dp) :: RTOL(NVAR)
! STEPMIN - Lower bound for integration step
  REAL(kind=dp) :: STEPMIN
! STEPMAX - Upper bound for integration step
  REAL(kind=dp) :: STEPMAX
! CFACTOR - Conversion factor for concentration units
  REAL(kind=dp) :: CFACTOR

! INLINED global variable declarations

  ! MECCA info from xmecca:
  CHARACTER(LEN=*), PUBLIC, PARAMETER :: &
    timestamp            = 'xmecca was run on 2020-04-16 at 16:48:52 by my9453 on machine fh2n1995.localdomain', &
    gas_spc_file         = '-rw-r--r-- 1 my9453 fh2-project-iconart 41911 Dec  3 08:55 gas.spc', &
    aqueous_spc_file     = '-rw------- 1 my9453 fh2-project-iconart 8444 Dec  3 08:55 aqueous.spc', &
    gas_eqn_file         = '-rw-r-x--- 1 my9453 fh2-project-iconart 153236 Dec  3 08:55 gasjs.eqn', &
    aqueous_eqn_file     = '-rw------- 1 my9453 fh2-project-iconart 56776 Dec  3 08:55 aqueous.eqn', &
    gas_spc_file_sum     = '13492    41', &
    aqueous_spc_file_sum = '54945     9', &
    gas_eqn_file_sum     = '62550   150', &
    aqueous_eqn_file_sum = '12944    56', &
    rplfile              = 'radmka', &
    wanted               = 'Js', &
    diagtracfile         = '', &
    rxnrates             = 'n', &
    tagdbl               = 'n'
  LOGICAL, PARAMETER :: REQ_MCFCT = .FALSE.

  ! from xmecca for aerosol:
  INTEGER, PARAMETER, PUBLIC :: APN = 1
  ! from aerosol.awk:
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_O2_a        = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_O3_a        = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_OH_a        = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HO2_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_H2O_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_H2O2_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_NH3_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_NO_a        = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_NO2_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_NO3_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HONO_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HNO3_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HNO4_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_N2O5_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_CH3OH_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HCOOH_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HCHO_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_CH3O2_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_CH3OOH_a    = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_CO2_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_CH3CO2H_a   = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_PAN_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_C2H5O2_a    = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_CH3CHO_a    = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_CH3COCH3_a  = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_Cl_a        = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_Cl2_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HCl_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HOCl_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_Br_a        = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_Br2_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HBr_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HOBr_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_BrCl_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_I2_a        = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_IO_a        = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HI_a        = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HOI_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_ICl_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_IBr_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HIO3_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_SO2_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_H2SO4_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_DMS_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_DMSO_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_Hg_a        = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HgO_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HgOH_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HgOHOH_a    = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HgOHCl_a    = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HgCl2_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HgBr2_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HgSO3_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_ClHgBr_a    = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_BrHgOBr_a   = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_ClHgOBr_a   = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_O2m_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_OHm_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_Hp_a        = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_NH4p_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_NO2m_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_NO3m_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_NO4m_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_CO3m_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HCOOm_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HCO3m_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_CH3COOm_a   = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_Clm_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_Cl2m_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_ClOm_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_ClOHm_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_Brm_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_Br2m_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_BrOm_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_BrOHm_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_BrCl2m_a    = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_Br2Clm_a    = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_Im_a        = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_IO2m_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_IO3m_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_ICl2m_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_IClBrm_a    = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_IBr2m_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_SO3m_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_SO3mm_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_SO4m_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_SO4mm_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_SO5m_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HSO3m_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HSO4m_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HSO5m_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_CH3SO3m_a   = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_CH2OHSO3m_a = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_Hgp_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_Hgpp_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HgOHp_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HgClp_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HgCl3m_a    = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HgCl4mm_a   = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HgBrp_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HgBr3m_a    = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HgBr4mm_a   = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_HgSO32mm_a  = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_D1O_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_D2O_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_DAHp_a      = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_DA_a        = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_DAm_a       = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_DGtAi_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_DGtAs_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_PROD1_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_PROD2_a     = 0
  INTEGER, PUBLIC, DIMENSION(APN) :: ind_Nap_a       = 0

  ! from gas.eqn:
  REAL :: k_HO2_HO2, k_NO3_NO2, k_NO2_HO2, k_HNO3_OH, k_CH3O2, k_CH3OOH_OH, &
          k_CH3CO3_NO2, k_PAN_M, k_ClO_ClO, k_BrO_NO2, k_I_NO2, k_DMS_OH, &
          G7402a_yield
  ! mz_dt_20140617+ new methane chemistry
  REAL :: beta_null_CH3NO3, beta_inf_CH3NO3, beta_CH3NO3, &
          k_NO2_CH3O2, k_G4138, k_G9408
  ! mz_dt_20140617-
  REAL :: KRO2NO, KRO2HO2, KAPHO2, KAPNO, KRO2NO3, KNO3AL
  REAL :: J_IC3H7NO3, J_ACETOL
  REAL(dp) :: RO2      ! sum of peroxy radicals
  REAL :: k_O3s ! mz_pj_20080812
  ! mz_rs_20140904+ only for mim1
  REAL :: k_PrO2_HO2, k_PrO2_NO, k_PrO2_CH3O2
  REAL :: k0_NO_HO2, k1d_NO_HO2, k1w_NO_HO2, k2d_NO_HO2, k2w_NO_HO2, &
          alpha_NO_HO2, beta_NO_HO2
  ! mz_rs_20140904- only for mim1
  ! mz_dt_20140604+
  ! SAR for OH reactions (reference 3030)
  ! parameters for the SAR on H-abstraction OH
  REAL(kind=dp), PARAMETER :: &
                     k_s = 8.42E-13 , k_t = 1.75E-12, k_p = 1.24E-13 &
                   , k_rohro = 1.6E-13 & ! value for CH3CH2OH + OH
                   , f_soh = 3.44   , f_toh = 2.68 &
                   , f_sooh = 7. , f_tooh = 7. & !from Wang and Chen, AE 2008
                   , f_ono2 = 0.04  , f_ch2ono2 = 0.2 &
                   , f_cpan= .25 , f_allyl = 3.6 , f_alk= 1.23 &
                   , f_cho = 0.55 , f_co2h = 1.67 , f_co = 0.73 &
                   , f_o = 8.15 , f_pch2oh = 1.29, f_tch2oh = 0.53
 ! parameters for the SAR on OH-addtion to alkenes
 ! Note for ach2ooh calculated
 ! assuming RO2 distribution in Paulot et al(2009)
 ! B,E,D,F,AC are .41, .02, .22, .05, .30 respectively
  REAL(kind=dp), PARAMETER :: &
                     k_adp = 0.45E-11 , k_ads = 3.0E-11, k_adt = 5.5E-11 &
                   , k_adsecprim = 3.0E-11 , k_adtertprim = 5.7E-11 &
                   , a_pan = 0.56       , a_cho = 0.31    , a_coch3 = 0.76 &
                   , a_ch2ono2 = 0.47   , a_ch2oh = 1.7   , a_ch2ooh = 0.21 &
                   , a_coh = 2.2        , a_cooh = 2.2    , a_co2h = 0.25
  ! mz_dt_201405604-
  ! mz_dt_20140121+
  ! Branching ratios for RCO3 + HO2 reactions:
  REAL(dp), PARAMETER :: RCO3_OH = 0.69, RCO3_O3 = 0.10, RCO3_OOH = 0.21
  ! mz_dt_20140121-
  ! ALL PHASES:
!KPPPP_DIRECTIVE vector variable definition start
  ! -------------------------------------------------------
  ! IMPORTANT NOTES ABOUT temp, PRESS, AND CAIR:
  ! - The KPP variable "temp" is already defined automatically by KPP in
  !   messy_mecca_kpp_global.f90. The KPP variables "press" and "cair" are
  !   defined here.
  ! - The 3 variables temp, press, and cair are only used inside KPP.
  !   They are different from the variables with the same names in the base
  !   model (as used in the SMIL files *_si.f90 and *_box.f90)
  ! - Data transfer between the SMIL and the KPP variables is done via the
  !   fill subroutines in messy_mecca_kpp.f90:
  !   - fill_temp transfers temperature
  !   - fill_press transfers pressure
  !   - fill_cair transfers cair (this is redundant because cair could be
  !     calculated from temp and press; nevertheless, redundant transfer to
  !     KPP is preferred to avoid potential numerical differences when
  !     recalculating cair inside KPP)
  ! - qqq describe potential KP4 incompatibility here...
  ! -------------------------------------------------------
  REAL(dp) :: cair    ! c(air) (wet) [mcl/cm^3]
  REAL(dp) :: press   ! pressure [Pa]
  ! mz_ab_20101119+
  REAL(dp) :: temp_ion ! ion temperature [K]
  REAL(dp) :: temp_elec! electron temperature [K]
  ! mz_ab_20101119-
!KPPPP_DIRECTIVE vector variable definition end
  ! AEROSOL ONLY:
!KPPPP_DIRECTIVE vector variable definition start
  REAL(dp) :: xaer(APN)
  REAL(dp) :: cvfac(APN)    ! unit conversion factor
  REAL(dp) :: lwc(APN)      ! liquid water content
  REAL(dp) :: k_exf(APN,NSPEC) = 0.
  REAL(dp) :: k_exb(APN,NSPEC) = 0.
  REAL(dp) :: k_exf_N2O5(APN)  = 0.
  REAL(dp) :: k_exf_ClNO3(APN) = 0.
  REAL(dp) :: k_exf_BrNO3(APN) = 0.
!KPPPP_DIRECTIVE vector variable definition end
  INTEGER, PUBLIC  :: xnom7sulf = 1 ! = 1-xm7sulf
!KPPPP_DIRECTIVE vector variable definition start
  REAL(dp) :: jx(IP_MAX) = 0.
!KPPPP_DIRECTIVE vector variable definition end
  ! iht_ = index of troposheric heterogeneous reactions
  INTEGER, PARAMETER, PUBLIC :: &
    iht_N2O5      =  1, iht_HNO3      =  2, iht_Hg      =  3, &
    iht_RGM       = 4
  INTEGER, PARAMETER, PUBLIC :: IHT_MAX = 4
!KPPPP_DIRECTIVE vector variable definition start
  REAL(dp) :: khet_Tr(IHT_MAX) = 0.
!KPPPP_DIRECTIVE vector variable definition end
  ! ihs_ = index of stratospheric heterogeneous reactions
  ! (must be the same as in messy_msbm.f90!)
  INTEGER, PARAMETER :: &
    ihs_N2O5_H2O  =  1, ihs_HOCl_HCl  =  2, ihs_ClNO3_HCl =  3, &
    ihs_ClNO3_H2O =  4, ihs_N2O5_HCl  =  5, ihs_ClNO3_HBr =  6, &
    ihs_BrNO3_HCl =  7, ihs_HOCl_HBr  =  8, ihs_HOBr_HCl  =  9, &
    ihs_HOBr_HBr  = 10, ihs_BrNO3_H2O = 11, ihs_Hg        = 12, &
    ihs_RGM       = 13
  INTEGER, PARAMETER, PUBLIC :: IHS_MAX = 13
!KPPPP_DIRECTIVE vector variable definition start
  REAL(dp) :: khet_St(IHS_MAX) = 0.
!KPPPP_DIRECTIVE vector variable definition end
  ! Parameters included for acid-base equilibria calculation
  ! used to enable the double use of the aqueous.eqn for liquid
  ! and aerosol phase.
  REAL(dp), PARAMETER :: &
    testfac_HO2   = 1.e5_dp, testfac_HONO   = 1.e5_dp, &
    testfac_HNO3  = 1.e7_dp, testfac_HNO4   = 1.e5_dp, &
    testfac_HCOOH = 1.e5_dp, testfac_SO2    = 1.e9_dp, &
    testfac_HSO3m = 1.e9_dp, testfac_HSO4m  = 1.e7_dp, &
    testfac_NH3   = 1.e7_dp, testfac_H2O    = 1.e9_dp, &
    testfac_CO2   = 1.e5_dp, testfac_HCl    = 1.e2_dp, &
    testfac_HBr   = 1.e6_dp,                           &
    testfac_HOCl  = 1.e2_dp, testfac_HOBr   = 1.e2_dp, &
    testfac_ICl   = 1.e2_dp, testfac_IBr    = 1.e2_dp, &
    testfac_IClBr = 1.e2_dp, testfac_H2SO4  = 1.e7_dp

  ! from xmecca:
  LOGICAL, PARAMETER :: REQ_HET     = .FALSE.
  LOGICAL, PARAMETER :: REQ_PHOTRAT = .TRUE.
  LOGICAL, PARAMETER :: REQ_AEROSOL = .FALSE.

  ! from xmecca:
  INTEGER, PARAMETER, PUBLIC :: MAX_MCEXP = 1
!KPPPP_DIRECTIVE vector variable definition start
  REAL :: mcexp(MAX_MCEXP) ! dummy Monte-Carlo factor
!KPPPP_DIRECTIVE vector variable definition end

  ! KPP info from xmecca (via integr.kpp):
  CHARACTER(LEN=*), PUBLIC, PARAMETER :: &
    mecca_spc_file     = '-rw------- 1 my9453 fh2-project-iconart 50765 Apr 16 16:48 mecca.spc', &
    mecca_eqn_file     = '-rw-r--r-- 1 my9453 fh2-project-iconart 32401 Apr 16 16:48 mecca.eqn', &
    mecca_spc_file_sum = '30485    50', &
    mecca_eqn_file_sum = '12142    32', &
    kppoption          = 'k', &
    KPP_HOME           = '/pfs/data4/project/fh2-project-iconart/my9453/caaba_3.0/mecca/kpp', &
    KPP_version        = '2.2.1_rs5', &
    integr             = 'rosenbrock_mz'

! INLINED global variable declarations


 END MODULE messy_mecca_kpp_global 

