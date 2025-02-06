!
! mo_mecicon_init
! This module wraps the box model MECCA (Code based on
! Sander, Rolf, et al. "Technical note: The new comprehensive
!        atmospheric chemistry module MECCA."
!        Atmospheric Chemistry and Physics 5.2 (2005): 445-450.
!Sander,R. et al.: The atmospheric chemistry box model CAABA/MECCA-3.0,
!       Geosci. Model Dev., 4, 373-380, doi:10.5194/gmd-4-373-2011, 2011.
!
!
!
!   KPP - symbolic chemistry Kinetics PreProcessor
!         (http://www.cs.vt.edu/~asandu/Software/KPP)
!   KPP is distributed under GPL, the general public licence
!         (http://www.gnu.org/copyleft/gpl.html)
!     (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
!     (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!         with contributions from:
!         R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany

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

MODULE mo_mecicon_init

    USE mo_kind,                 ONLY: wp
    USE mo_art_mecicon_data,     ONLY: R_gas, N_A, TINY_DP, HLINE1, OneDay, &
                                   &   STRLEN_VLONG, t_art_mecicon_utils
    USE messy_mecca_kpp_parameters
    USE messy_cmn_photol_mem
    USE mo_util_string,          ONLY: toupper
    USE messy_mecca_kpp_global
    
    IMPLICIT NONE 

    PUBLIC :: mecicon_init, initialize_kpp_variables
             

    CONTAINS

    SUBROUTINE mecicon_init(pres_icon, temp_icon,vmr_kpp, kpp_ind)
        REAL(wp), intent(in)        :: pres_icon, temp_icon
        REAL(wp), intent(inout)     :: vmr_kpp
        INTEGER, intent(in)         :: kpp_ind


       ! IF (.NOT. ALLOCATED(c)) ALLOCATE(c(NSPEC))

        temp  = temp_icon
        press = pres_icon
        cair  = (N_A/1.E6) * press / (R_gas*temp)

        c(:)  = 0. ! default value unless explicitly initialized
        
        CALL initialize_kpp_variables

        CALL x0                         ! set initial number concentrations

        vmr_kpp = c(kpp_ind)/cair       ! convert #/cm3 -> mol/mol

    END SUBROUTINE mecicon_init




    SUBROUTINE initialize_kpp_variables
        USE messy_mecca_kpp_initialize, ONLY: initialize
        USE messy_mecca_kpp_parameters  ! ind_*
        USE messy_mecca_kpp_global
        

        IMPLICIT NONE

    ! initialize kpp variables
        CALL initialize

        rtol(:) = 1E-4_DP ! relative tolerance
        atol(:) = 1.E-1_DP   ! absolute tolerance
!    rtol(:) = 1E-2_dp ! relative tolerance
!    atol(:) = 1E1_dp   ! absolute tolerance

        IF ((ind_OH  >0).AND.(ind_OH  <=NVAR)) atol(ind_OH)  = 1._dp
        IF ((ind_NO3 >0).AND.(ind_NO3 <=NVAR)) atol(ind_NO3) = 1._dp
        IF ((ind_Cl  >0).AND.(ind_Cl  <=NVAR)) atol(ind_Cl)  = 1._dp
        IF ((ind_Br  >0).AND.(ind_Br  <=NVAR)) atol(ind_Br)  = 1._dp
        IF ((ind_O1D >0).AND.(ind_O1D <=NVAR)) atol(ind_O1D) = 1._dp

    END SUBROUTINE initialize_kpp_variables


    SUBROUTINE x0

    ! ------------------------------------------------------------------------

    ! initialize some mixing ratios
    ! values in mol/mol, cair converts to particles/cm3
  !  SELECT CASE (TRIM(init_scenario))
  !  CASE ('') ! default if no scenario selected in namelist:
  !    print *, "c0 chapman"
  
!      CALL x0_radmka
      CALL x0_ccmi
!      CALL x0_psc
      
   !   CALL x0_chapman
!      CALL x0_free_trop
 !     CALL x0_free_trop
 !     CALL x0_mbl
  !  CASE ('FF_ANTARCTIC')
  !    CALL x0_ff_antarctic
  !  CASE ('FF_ARCTIC')
  !    CALL x0_ff_arctic
  !  CASE ('FREE_TROP')
  !    CALL x0_free_trop
  !  CASE ('HOOVER')
  !    CALL x0_hoover
  !  CASE ('LAB','LAB_C15')
  !    CALL x0_lab
  !  CASE ('MBL')
  !    CALL x0_mbl
  !  CASE ('MIM2')
  !    CALL x0_mim2
  !  CASE ('MTCHEM')
  !    CALL x0_mtchem
  !  CASE ('OOMPH')
  !    CALL x0_oomph
  !  CASE ('ISOO')
  !    CALL x0_isoo
  !  CASE ('STRATO')
  !    ! choose one:
  !    CALL x0_strato10
  !    !CALL x0_strato20
  !  CASE DEFAULT
  !    PRINT *, 'ERROR, init_scenario '//TRIM(init_scenario)//' is not defined'
  !    STOP
  !  END SELECT

    ! ------------------------------------------------------------------------
  CONTAINS

    ! ------------------------------------------------------------------------
    SUBROUTINE x0_chapman
    !  IF (ind_O2       /= 0) c(ind_O2)      = 3.602181818181818e-3*cair 
    !  IF (ind_N2       /= 0) c(ind_N2)      = 0.002330909090909091 *cair
    !  IF (ind_O1D      /= 0) c(ind_O1D)     = 3.602181818181818e-18*cair
    !  IF (ind_O3P      /= 0) c(ind_O3P)     = 2.4088727272727273e-8*cair
    !  IF (ind_NO      /= 0) c(ind_NO)       = 2.7E-12*cair
    !  IF (ind_NO2      /= 0) c(ind_NO2)     = 4.2E-10*cair
    !  IF (ind_O3     /= 0) c(ind_O3)        = 1.9367272727272727e-06*cair
    !print *, "ind_O2", ind_O2
    !print*, c(ind_O2)
 !    IF (ind_O1D      /= 0) c(ind_O1D)    =   1.E-16 * cair
 !    IF (ind_O3P      /= 0) c(ind_O3P)    =   9.E-34 * cair
 !     IF (ind_O2       /= 0) c(ind_O2)     = 210.E-03 * cair
 !     IF (ind_N        /= 0) c(ind_N)      =   1.E-12 * cair
 !     IF (ind_NO2      /= 0) c(ind_NO2)    =   1.E-09 * cair
 !     IF (ind_O3       /= 0) c(ind_O3)     =   4.E-06 * cair

      IF (ind_O1D      /= 0) c(ind_O1D)    =   5.4E-17 * cair
      IF (ind_O3P      /= 0) c(ind_O3P)    =   8.9E-11 * cair
      IF (ind_O2       /= 0) c(ind_O2)     = 210.E-03 * cair
      IF (ind_N        /= 0) c(ind_N)      =   1.E-12 * cair
      IF (ind_NO2      /= 0) c(ind_NO2)    =   2.9E-11 * cair
    ! IF (ind_O3       /= 0) c(ind_O3)     =   10.0E12
      IF (ind_O3       /= 0) c(ind_O3)     =   8.3E-6*cair
      IF (ind_N2       /= 0) c(ind_N2)      = 780.E-03 * cair
      IF (ind_NO       /= 0) c(ind_NO)     = 2.7E-9* cair
      IF (ind_HNO3        /= 0) c(ind_HNO3)        =  1.33769E-10 * cair
      IF (ind_H2O         /= 0) c(ind_H2O)         =  0.0038177   * cair
      !IF (ind_OH          /= 0) c(ind_OH)          =  7.07567E-14 * cair
      IF (ind_OH          /= 0) c(ind_OH)          =  7.07567E-11 * cair
    ! IF (ind_HO2         /= 0) c(ind_HO2)         =  9.74215E-13 * cair
      IF (ind_HO2         /= 0) c(ind_HO2)         =  9.74215E-10 * cair
      IF (ind_NO3         /= 0) c(ind_NO3)         =  3.1458E-12  * cair
      IF (ind_N2O5        /= 0) c(ind_N2O5)        =  4.37511E-12 * cair
      IF (ind_CF2Cl2      /= 0) c(ind_CF2Cl2)      =  7.87179E-10 * cair
      IF (ind_Cl          /= 0) c(ind_Cl)          =  7.12196E-17 * cair
      IF (ind_ClO         /= 0) c(ind_ClO)         =  8.18961E-14 * cair
      IF (ind_CH4         /= 0) c(ind_CH4)         =  1.74813E-06 * cair
      IF (ind_HCl         /= 0) c(ind_HCl)         =  5.60562E-11 * cair
      IF (ind_ClNO3       /= 0) c(ind_ClNO3)       =  3.81274E-12 * cair
      IF (ind_N2O         /= 0) c(ind_N2O)         =  3.16661E-07 * cair
      IF (ind_CH3O2       /= 0) c(ind_CH3O2)       =  1.1287E-12  * cair

  END SUBROUTINE x0_chapman
    SUBROUTINE x0_simple
      IF (ind_O2       /= 0) c(ind_O2)      = 210.E-03 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 780.E-03 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)     =  1.8E-06 * cair
      IF (ind_CO       /= 0) c(ind_CO)      =  70.E-09 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)     = 350.E-06 * cair
    END SUBROUTINE x0_simple

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_ff_antarctic
      IF (ind_O3       /= 0) c(ind_O3)      =  30.E-09 * cair
      IF (ind_O2       /= 0) c(ind_O2)      = 210.E-03 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 780.E-03 * cair
      IF (ind_NO       /= 0) c(ind_NO)      =  2.0E-12 * cair
      IF (ind_NO2      /= 0) c(ind_NO2)     =  2.0E-12 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)     =  1.8E-06 * cair
      IF (ind_HCHO     /= 0) c(ind_HCHO)    =  50.E-12 * cair
      IF (ind_CH3CHO   /= 0) c(ind_CH3CHO)  =  10.E-12 * cair
      IF (ind_CO       /= 0) c(ind_CO)      = 170.E-09 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)     = 350.E-06 * cair
      IF (ind_DMS      /= 0) c(ind_DMS)     =  10.E-12 * cair
      IF (ind_SO2      /= 0) c(ind_SO2)     =  30.E-12 * cair
      IF (ind_CH3I     /= 0) c(ind_CH3I)    =  2.0E-12 * cair
      IF (ind_CH3Br    /= 0) c(ind_CH3Br)   =  5.0E-12 * cair
      IF (ind_C2H4     /= 0) c(ind_C2H4)    =  10.E-12 * cair
      IF (ind_C2H2     /= 0) c(ind_C2H2)    =  10.E-12 * cair
      IF (ind_C2H6     /= 0) c(ind_C2H6)    = 300.E-12 * cair
      IF (ind_Hg       /= 0) c(ind_Hg)      = 1.68E-13 * cair
    END SUBROUTINE x0_ff_antarctic

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_ff_arctic
      IF (ind_O3       /= 0) c(ind_O3)      =  40.E-09 * cair
      IF (ind_O2       /= 0) c(ind_O2)      = 210.E-03 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 780.E-03 * cair
      IF (ind_NO       /= 0) c(ind_NO)      =  10.E-12 * cair
      IF (ind_NO2      /= 0) c(ind_NO2)     =  10.E-12 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)     =  1.8E-06 * cair
      IF (ind_HCHO     /= 0) c(ind_HCHO)    = 200.E-12 * cair
      IF (ind_CH3CHO   /= 0) c(ind_CH3CHO)  = 100.E-12 * cair
      IF (ind_CO       /= 0) c(ind_CO)      = 170.E-09 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)     = 350.E-06 * cair
      IF (ind_DMS      /= 0) c(ind_DMS)     =  10.E-12 * cair
      IF (ind_SO2      /= 0) c(ind_SO2)     = 100.E-12 * cair
      IF (ind_CH3I     /= 0) c(ind_CH3I)    =  2.0E-12 * cair
      IF (ind_CH3Br    /= 0) c(ind_CH3Br)   =  5.0E-12 * cair
      IF (ind_C2H4     /= 0) c(ind_C2H4)    =   26E-12 * cair ! ref1737
      ! see also: C2H4 = 100E-12 ! ref0351, Tab.1, 3 Apr
      IF (ind_C2H2     /= 0) c(ind_C2H2)    =  329E-12 * cair ! ref1737
      ! see also: C2H2 = 840E-12 ! ref0351, Tab.1, 3 Apr
      IF (ind_C2H6     /= 0) c(ind_C2H6)    =  2.0E-09 * cair
      IF (ind_Hg       /= 0) c(ind_Hg)      = 1.68E-13 * cair
    END SUBROUTINE x0_ff_arctic

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_hoover
      IF (ind_O2       /= 0) c(ind_O2)      = 210.E-03 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 780.E-03 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)     = 350.E-06 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)     = 1.8E-06 * cair
      IF (ind_H2       /= 0) c(ind_H2)      = 0.6E-06 * cair
    END SUBROUTINE x0_hoover

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_free_trop
      ! average init from CARIBIC2 trajectories 5-day back init
      ! flight: 20060706_CAN_FRA_157_EMAC
      IF (ind_CH3OH       /= 0) c(ind_CH3OH)       =  7.23043E-10 * cair
      IF (ind_HCHO        /= 0) c(ind_HCHO)        =  2.53104E-10 * cair
      IF (ind_CO          /= 0) c(ind_CO)          =  7.29323E-08 * cair
      IF (ind_HCOOH       /= 0) c(ind_HCOOH)       =  4.13365E-11 * cair
      IF (ind_CH3CHO      /= 0) c(ind_CH3CHO)      =  2.00554E-11 * cair
      IF (ind_CH3CO2H     /= 0) c(ind_CH3CO2H)     =  6.93619E-11 * cair
      IF (ind_CH3COCH3    /= 0) c(ind_CH3COCH3)    =  3.34549E-10 * cair
      IF (ind_ACETOL      /= 0) c(ind_ACETOL)      =  5.64005E-11 * cair
      IF (ind_MGLYOX      /= 0) c(ind_MGLYOX)      =  2.03489E-11 * cair
      IF (ind_BIACET      /= 0) c(ind_BIACET)      =  4.19846E-13 * cair
      IF (ind_Br          /= 0) c(ind_Br)          =  3.38058E-14 * cair
      IF (ind_Br2         /= 0) c(ind_Br2)         =  3.87407E-16 * cair
      IF (ind_BrO         /= 0) c(ind_BrO)         =  2.39308E-13 * cair
      IF (ind_HBr         /= 0) c(ind_HBr)         =  5.70085E-13 * cair
      IF (ind_HOBr        /= 0) c(ind_HOBr)        =  2.86479E-13 * cair
      IF (ind_BrNO2       /= 0) c(ind_BrNO2)       =  0.          * cair
      IF (ind_BrNO3       /= 0) c(ind_BrNO3)       =  6.53232E-13 * cair
      IF (ind_BrCl        /= 0) c(ind_BrCl)        =  9.05211E-17 * cair
      IF (ind_Cl          /= 0) c(ind_Cl)          =  7.12196E-17 * cair
      IF (ind_Cl2         /= 0) c(ind_Cl2)         =  3.46324E-17 * cair
      IF (ind_ClO         /= 0) c(ind_ClO)         =  8.18961E-14 * cair
      IF (ind_HCl         /= 0) c(ind_HCl)         =  5.60562E-11 * cair
      IF (ind_HOCl        /= 0) c(ind_HOCl)        =  2.5891E-13  * cair
      IF (ind_Cl2O2       /= 0) c(ind_Cl2O2)       =  3.47767E-17 * cair
      IF (ind_OClO        /= 0) c(ind_OClO)        =  1.33148E-16 * cair
      IF (ind_ClNO3       /= 0) c(ind_ClNO3)       =  3.81274E-12 * cair
      IF (ind_CCl4        /= 0) c(ind_CCl4)        =  8.9102E-11  * cair
      IF (ind_CH3Cl       /= 0) c(ind_CH3Cl)       =  4.55034E-10 * cair
      IF (ind_CH3CCl3     /= 0) c(ind_CH3CCl3)     =  1.50035E-11 * cair
      IF (ind_CF2Cl2      /= 0) c(ind_CF2Cl2)      =  7.87179E-10 * cair
      IF (ind_CFCl3       /= 0) c(ind_CFCl3)       =  2.42555E-10 * cair
      IF (ind_CH2ClBr     /= 0) c(ind_CH2ClBr)     =  9.40689E-14 * cair
      IF (ind_CHCl2Br     /= 0) c(ind_CHCl2Br)     =  5.51093E-14 * cair
      IF (ind_CHClBr2     /= 0) c(ind_CHClBr2)     =  4.54419E-14 * cair
      IF (ind_CH2Br2      /= 0) c(ind_CH2Br2)      =  1.02541E-12 * cair
      IF (ind_CH3Br       /= 0) c(ind_CH3Br)       =  8.78744E-15 * cair
      IF (ind_CHBr3       /= 0) c(ind_CHBr3)       =  4.8533E-13  * cair
      IF (ind_CF3Br       /= 0) c(ind_CF3Br)       =  2.43175E-12 * cair
      IF (ind_CF2ClBr     /= 0) c(ind_CF2ClBr)     =  4.20382E-12 * cair
      IF (ind_CH4         /= 0) c(ind_CH4)         =  1.74813E-06 * cair
      IF (ind_C2H6        /= 0) c(ind_C2H6)        =  5.04266E-10 * cair
      IF (ind_C2H4        /= 0) c(ind_C2H4)        =  1.51334E-11 * cair
      IF (ind_C3H8        /= 0) c(ind_C3H8)        =  4.6215E-11  * cair
      IF (ind_C3H6        /= 0) c(ind_C3H6)        =  1.68475E-12 * cair
      IF (ind_NC4H10      /= 0) c(ind_NC4H10)      =  7.76293E-11 * cair
      IF (ind_MVK         /= 0) c(ind_MVK)         =  3.90326E-11 * cair
      IF (ind_MEK         /= 0) c(ind_MEK)         =  3.98487E-11 * cair
      IF (ind_C5H8        /= 0) c(ind_C5H8)        =  8.45986E-12 * cair
      IF (ind_HgO         /= 0) c(ind_HgO)         =  3.57349E-16 * cair
      IF (ind_HgBr2       /= 0) c(ind_HgBr2)       =  9.89811E-16 * cair
      IF (ind_ClHgBr      /= 0) c(ind_ClHgBr)      =  1.13395E-16 * cair
      IF (ind_BrHgOBr     /= 0) c(ind_BrHgOBr)     =  7.30999E-15 * cair
      IF (ind_ClHgOBr     /= 0) c(ind_ClHgOBr)     =  1.1631E-15  * cair
      IF (ind_HgCl        /= 0) c(ind_HgCl)        =  6.52374E-17 * cair
      IF (ind_HgBr        /= 0) c(ind_HgBr)        =  7.25411E-16 * cair
      IF (ind_Hg          /= 0) c(ind_Hg)          =  1.75258E-13 * cair
      IF (ind_NACA        /= 0) c(ind_NACA)        =  2.85441E-12 * cair
      IF (ind_MPAN        /= 0) c(ind_MPAN)        =  5.31725E-12 * cair
      IF (ind_IC3H7NO3    /= 0) c(ind_IC3H7NO3)    =  8.94082E-13 * cair
      IF (ind_LC4H9NO3    /= 0) c(ind_LC4H9NO3)    =  8.32727E-12 * cair
      IF (ind_ISON        /= 0) c(ind_ISON)        =  9.0462E-12  * cair
      IF (ind_N2          /= 0) c(ind_N2)          =  0.78        * cair
      IF (ind_NH3         /= 0) c(ind_NH3)         =  1.2068E-10  * cair
      IF (ind_N2O         /= 0) c(ind_N2O)         =  3.16661E-07 * cair
      IF (ind_NO          /= 0) c(ind_NO)          =  4.65026E-11 * cair
      IF (ind_NO2         /= 0) c(ind_NO2)         =  1.25056E-10 * cair
      IF (ind_NO3         /= 0) c(ind_NO3)         =  3.1458E-12  * cair
      IF (ind_N2O5        /= 0) c(ind_N2O5)        =  4.37511E-12 * cair
      IF (ind_HONO        /= 0) c(ind_HONO)        =  6.12167E-13 * cair
      IF (ind_HNO3        /= 0) c(ind_HNO3)        =  1.33769E-10 * cair
      IF (ind_HNO4        /= 0) c(ind_HNO4)        =  2.77194E-11 * cair
      IF (ind_PAN         /= 0) c(ind_PAN)         =  3.1475E-10  * cair
      IF (ind_NH2OH       /= 0) c(ind_NH2OH)       =  1.31453E-12 * cair
      IF (ind_NHOH        /= 0) c(ind_NHOH)        =  1.00251E-11 * cair
      IF (ind_HNO         /= 0) c(ind_HNO)         =  9.47677E-14 * cair
      IF (ind_NH2         /= 0) c(ind_NH2)         =  8.9126E-18  * cair
      IF (ind_O3P         /= 0) c(ind_O3P)         =  2.31275E-15 * cair
      IF (ind_O2          /= 0) c(ind_O2)          =  0.21        * cair
      ! also lower strat in average
      IF (ind_O3          /= 0) c(ind_O3)          =  1.60052E-07 * cair
      IF (ind_H2          /= 0) c(ind_H2)          =  5.35521E-07 * cair
      IF (ind_OH          /= 0) c(ind_OH)          =  7.07567E-14 * cair
      IF (ind_HO2         /= 0) c(ind_HO2)         =  9.74215E-13 * cair
      IF (ind_H2O2        /= 0) c(ind_H2O2)        =  6.10601E-10 * cair
      IF (ind_H2O         /= 0) c(ind_H2O)         =  0.0038177   * cair
      IF (ind_CH3O2       /= 0) c(ind_CH3O2)       =  1.1287E-12  * cair
      IF (ind_CH3OOH      /= 0) c(ind_CH3OOH)      =  2.37819E-10 * cair
      IF (ind_C2H5O2      /= 0) c(ind_C2H5O2)      =  1.21222E-14 * cair
      IF (ind_C2H5OOH     /= 0) c(ind_C2H5OOH)     =  9.14669E-12 * cair
      IF (ind_CH3CO3      /= 0) c(ind_CH3CO3)      =  9.65657E-14 * cair
      IF (ind_CH3CO3H     /= 0) c(ind_CH3CO3H)     =  5.72797E-11 * cair
      IF (ind_IC3H7O2     /= 0) c(ind_IC3H7O2)     =  3.568E-15   * cair
      IF (ind_IC3H7OOH    /= 0) c(ind_IC3H7OOH)    =  1.98726E-12 * cair
      IF (ind_LHOC3H6O2   /= 0) c(ind_LHOC3H6O2)   =  2.6452E-14  * cair
      IF (ind_LHOC3H6OOH  /= 0) c(ind_LHOC3H6OOH)  =  4.19966E-12 * cair
      IF (ind_CH3COCH2O2  /= 0) c(ind_CH3COCH2O2)  =  4.96396E-15 * cair
      IF (ind_HYPERACET   /= 0) c(ind_HYPERACET)   =  1.95516E-12 * cair
      IF (ind_LC4H9O2     /= 0) c(ind_LC4H9O2)     =  1.79502E-14 * cair
      IF (ind_LC4H9OOH    /= 0) c(ind_LC4H9OOH)    =  9.50258E-12 * cair
      IF (ind_MVKO2       /= 0) c(ind_MVKO2)       =  7.3364E-14  * cair
      IF (ind_MVKOOH      /= 0) c(ind_MVKOOH)      =  2.49326E-11 * cair
      IF (ind_LMEKO2      /= 0) c(ind_LMEKO2)      =  6.33844E-15 * cair
      IF (ind_LMEKOOH     /= 0) c(ind_LMEKOOH)     =  3.52573E-12 * cair
      IF (ind_ISO2        /= 0) c(ind_ISO2)        =  2.85094E-14 * cair
      IF (ind_ISOOH       /= 0) c(ind_ISOOH)       =  7.02518E-12 * cair
      IF (ind_SO2         /= 0) c(ind_SO2)         =  6.63355E-11 * cair
      IF (ind_H2SO4       /= 0) c(ind_H2SO4)       =  5.51045E-14 * cair
      IF (ind_CH3SO3H     /= 0) c(ind_CH3SO3H)     =  6.24971E-11 * cair
      IF (ind_DMS         /= 0) c(ind_DMS)         =  5.13875E-12 * cair
      IF (ind_DMSO        /= 0) c(ind_DMSO)        =  7.64715E-14 * cair
      IF (ind_CH3SO2      /= 0) c(ind_CH3SO2)      =  4.44594E-17 * cair
      IF (ind_CH3SO3      /= 0) c(ind_CH3SO3)      =  1.33623E-14 * cair
      IF (ind_CO2         /= 0) c(ind_CO2)         =  0.000382388 * cair
      IF (ind_SF6         /= 0) c(ind_SF6)         =  5.82353E-12 * cair
      IF (ind_CH3I        /= 0) c(ind_CH3I)        =  4.39702E-14 * cair
    END SUBROUTINE x0_free_trop

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_lab
      IF (ind_O2       /= 0) c(ind_O2)      = 210.E-03 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 780.E-03 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)     =  1.8E-06 * cair
      !qqq todo: more values from Sergej?
    END SUBROUTINE x0_lab

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_mbl
      IF (ind_H2       /= 0) c(ind_H2)      =   1.E-06 * cair
      IF (ind_O3       /= 0) c(ind_O3)      =  25.E-09 * cair
      IF (ind_O2       /= 0) c(ind_O2)      = 210.E-03 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 780.E-03 * cair
      !IF (ind_NO      /= 0) c(ind_NO)      =  10.E-12 * cair
      IF (ind_NO2      /= 0) c(ind_NO2)     =  20.E-12 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)     =  1.8E-06 * cair
      IF (ind_HCHO     /= 0) c(ind_HCHO)    = 300.E-12 * cair
      !IF (ind_CH3CHO  /= 0) c(ind_CH3CHO)  =  1.0E-10 * cair
      IF (ind_CO       /= 0) c(ind_CO)      =  70.E-09 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)     = 350.E-06 * cair
      IF (ind_DMS      /= 0) c(ind_DMS)     =  60.E-12 * cair
      IF (ind_SO2      /= 0) c(ind_SO2)     =  90.E-12 * cair
      IF (ind_CH3I     /= 0) c(ind_CH3I)    =  2.0E-12 * cair
      IF (ind_H2O2     /= 0) c(ind_H2O2)    = 600.E-12 * cair
      IF (ind_NH3      /= 0) c(ind_NH3)     = 200.E-12 * cair
      IF (ind_HNO3     /= 0) c(ind_HNO3)    =  5.0E-12 * cair
      IF (ind_HCl      /= 0) c(ind_HCl)     =  40.E-12 * cair
      IF (ind_C3H7I    /= 0) c(ind_C3H7I)   =  1.0E-12 * cair
      IF (ind_CH3OH    /= 0) c(ind_CH3OH)   = 300.E-12 * cair
      !IF (ind_C5H8    /= 0) c(ind_C5H8)    =  1.0E-09 * cair
      IF (ind_Hg       /= 0) c(ind_Hg)      = 1.68E-13 * cair

      ! examples for initializing aqueous-phase species:
      ! (qqq the index must not be greater than APN)
      ! IF (ind_NH4p_a(2)  /=0) c(ind_NH4p_a(2))  = 300.E-12 * cair
      ! 1 nmol/L DMS:
      ! IF (ind_DMS_a(3) /=0) c(ind_DMS_a(3)) = 1E-9 * lwc(3) * N_A / 1.E3 * cair

    END SUBROUTINE x0_mbl

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_mim2
      ! data from Tab. 2 of Taraborrelli et al., ACP, 9, 2751-2777 (2009):
      IF (ind_H2O2     /= 0) c(ind_H2O2)    =   7.E-09 * cair
      IF (ind_O3       /= 0) c(ind_O3)      =  30.E-09 * cair
      IF (ind_O2       /= 0) c(ind_O2)      = 210.E-03 * cair
      IF (ind_NH3      /= 0) c(ind_NH3)     = 100.E-12 * cair
      IF (ind_NO       /= 0) c(ind_NO)      =  10.E-12 * cair
      IF (ind_NO2      /= 0) c(ind_NO2)     = 100.E-12 * cair
      IF (ind_HONO     /= 0) c(ind_HONO)    =  40.E-14 * cair
      IF (ind_HNO3     /= 0) c(ind_HNO3)    =  5.0E-12 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 780.E-03 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)     =  1.8E-06 * cair
      IF (ind_HCHO     /= 0) c(ind_HCHO)    =  5.0E-09 * cair
      IF (ind_CH3OH    /= 0) c(ind_CH3OH)   = 500.E-12 * cair
      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH)  =  4.0E-09 * cair
      IF (ind_HCOOH    /= 0) c(ind_HCOOH)   = 350.E-12 * cair
      IF (ind_CO       /= 0) c(ind_CO)      = 100.E-09 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)     = 350.E-06 * cair
      IF (ind_CH3CO2H  /= 0) c(ind_CH3CO2H) =  2.0E-09 * cair
      IF (ind_CH3CO3H  /= 0) c(ind_CH3CO3H) =  1.5E-09 * cair
      IF (ind_ACETOL   /= 0) c(ind_ACETOL)  =  4.0E-09 * cair
      IF (ind_MGLYOX   /= 0) c(ind_MGLYOX)  = 500.E-12 * cair
      IF (ind_C5H8     /= 0) c(ind_C5H8)    =  2.0E-09 * cair
      IF (ind_PAN      /= 0) c(ind_PAN)     = 100.E-12 * cair
    END SUBROUTINE x0_mim2

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_oomph
      ! fixed species:
      IF (ind_CO2      /= 0) c(ind_CO2)     = 382.E-06 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 780.E-03 * cair
      IF (ind_O2       /= 0) c(ind_O2)      = 210.E-03 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)     = 1.75E-06 * cair ! JW

      ! def:   default as in usual mbl setup
      ! Heard: North Atlantic Campaign NAMBLEX
      ! JW:    J. Williams' guess
      ! RS:    Rolf Sander's guess
      ! Sa:    Sander's lit compilation
      ! Singh: Hanwant Singh's missions 15, 16, 18 (remote clean air)
      !        in tropical Pacific = ref0314
      ! Wa:    Warneck: Nat Atm + info (page) = ref0067
      ! mod<#> different modifications to reach steady state earlier

      ! REMOTE MARINE BACKGROUND
      PRINT *, 'OOMPH init: marine background'
      IF (ind_ACETOL   /= 0) c(ind_ACETOL)  = 205.E-12 * cair ! 3rdgen
      !IF (ind_ACETOL   /= 0) c(ind_ACETOL)  = 250.E-12 * cair ! new
      IF (ind_HYPERACET    /= 0) c(ind_HYPERACET)   =  0.74E-12 * cair ! 4Igen
      !IF (ind_HYPERACET    /= 0) c(ind_HYPERACET)   =  0.7E-12 * cair ! 4gen
      !IF (ind_CH3CHO   /= 0) c(ind_CH3CHO)  = 100.E-12 * cair ! def(off)
      IF (ind_C3H7I    /= 0) c(ind_C3H7I)   =  0.6E-12 * cair ! 4Igen
      !IF (ind_C3H7I    /= 0) c(ind_C3H7I)   =  0.7E-12 * cair ! mod5
      IF (ind_CH3CHO   /= 0) c(ind_CH3CHO)  = 100.E-12 * cair ! 4gen
      IF (ind_CH3COCH3 /= 0) c(ind_CH3COCH3)= 100.E-12 * cair ! 4gen
      IF (ind_MGLYOX /= 0) c(ind_MGLYOX)= 230.E-12 * cair ! 4Igen
      !IF (ind_MGLYOX /= 0) c(ind_MGLYOX)= 280.E-12 * cair ! 3rdgen
      !IF (ind_MGLYOX /= 0) c(ind_MGLYOX)= 340.E-12 * cair ! new
      IF (ind_CH3I     /= 0) c(ind_CH3I)    =  1.24E-12 * cair ! 4Igen
      !IF (ind_CH3I     /= 0) c(ind_CH3I)    =  1.3E-12 * cair ! mod8
      IF (ind_CH3OH    /= 0) c(ind_CH3OH)   = 500.E-12 * cair ! 4gen
      !IF (ind_CH3OH    /= 0) c(ind_CH3OH)   = 1.47E-09 * cair ! 3rdgen
      !IF (ind_CH3OH    /= 0) c(ind_CH3OH)   = 1.44E-09 * cair ! mod10
      !IF (ind_CH3OH    /= 0) c(ind_CH3OH)   = 800.E-12 * cair ! mod9+Sinha
      !IF (ind_CH3OH    /= 0) c(ind_CH3OH)   = 300.E-12 * cair ! def+JW
      !IF (ind_CH3OOH   /= 0) c(ind_CH3OOH)  = 0.E-12 * cair ! 4Igentest
      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH)  = 590.E-12 * cair ! 3rdgen
      !IF (ind_CH3OOH   /= 0) c(ind_CH3OOH)  = 530.E-12 * cair ! mod8
      !IF (ind_CH3CO2H  /= 0) c(ind_CH3CO2H) = 0.E-09 * cair ! 4Igentest
      IF (ind_CH3CO2H  /= 0) c(ind_CH3CO2H) = 1.34E-09 * cair ! 4Igen
      !IF (ind_CH3CO2H  /= 0) c(ind_CH3CO2H) = 1.54E-09 * cair ! 4gen
      !IF (ind_CH3CO2H  /= 0) c(ind_CH3CO2H) = 1.48E-09 * cair ! 3rdgen
      !IF (ind_CH3CO2H  /= 0) c(ind_CH3CO2H) = 1.50E-09 * cair ! new
      IF (ind_CO       /= 0) c(ind_CO)      =  35.E-09 * cair ! 4Igen
      !IF (ind_CO       /= 0) c(ind_CO)      = 115.E-09 * cair ! 3rdgen, 4gen
      !IF (ind_CO       /= 0) c(ind_CO)      = 113.E-09 * cair ! mod6
      IF (ind_DMS      /= 0) c(ind_DMS)     = 200.E-12 * cair ! 4gen
      !IF (ind_DMS      /= 0) c(ind_DMS)     = 130.E-12 * cair ! 3rdgen
      !IF (ind_DMS      /= 0) c(ind_DMS)     =  30.E-12 * cair ! mod8
      IF (ind_DMSO     /= 0) c(ind_DMSO)    =  18.E-12 * cair ! 4Igen
      IF (ind_H2       /= 0) c(ind_H2)      =   1.E-06 * cair ! def
      IF (ind_H2O2     /= 0) c(ind_H2O2)    = 220.E-12 * cair ! 4Igen
      !IF (ind_H2O2     /= 0) c(ind_H2O2)    = 260.E-12 * cair ! mod6
      !IF (ind_H2O2     /= 0) c(ind_H2O2)    = 552.E-12 * cair ! Singh
      IF (ind_HCHO     /= 0) c(ind_HCHO)    = 220.E-12 * cair ! 4Igen
      !IF (ind_HCHO     /= 0) c(ind_HCHO)    = 200.E-12 * cair ! 3rdgen
      !IF (ind_HCHO     /= 0) c(ind_HCHO)    = 180.E-12 * cair ! mod3
      !IF (ind_HCHO     /= 0) c(ind_HCHO)    = 200.E-12 * cair ! Sa
      IF (ind_HCOOH    /= 0) c(ind_HCOOH)   =   7.E-12 * cair ! 4Igen
      !IF (ind_HCOOH    /= 0) c(ind_HCOOH)   =   8.E-12 * cair ! 3rdgen
      !IF (ind_HCOOH    /= 0) c(ind_HCOOH)   =  11.E-12 * cair ! mod5
      IF (ind_HCl      /= 0) c(ind_HCl)     =  30.E-12 * cair ! 4gen
      !IF (ind_HCl      /= 0) c(ind_HCl)     =  40.E-12 * cair ! def,3gen
      !IF (ind_HCl      /= 0) c(ind_HCl)     = 100.E-12 * cair ! Sa
      IF (ind_HNO3     /= 0) c(ind_HNO3)    = 0.15E-12 * cair ! mod3
      IF (ind_ISON     /= 0) c(ind_ISON)    =  2.4E-12 * cair ! 4Igen
      !IF (ind_ISON     /= 0) c(ind_ISON)    =  3.4E-12 * cair ! 4gen
      !IF (ind_ISON     /= 0) c(ind_ISON)    =  2.2E-12 * cair ! 3rdgen
      !IF (ind_ISON     /= 0) c(ind_ISON)    =  2.8E-12 * cair ! new
      IF (ind_ISOOH    /= 0) c(ind_ISOOH)   =  30.E-12 * cair ! 4Igen
      !IF (ind_ISOOH    /= 0) c(ind_ISOOH)   =  40.E-12 * cair ! new
      !IF (ind_ISOOH    /= 0) c(ind_ISOOH)   =  40.E-12 * cair ! new
      IF (ind_C5H8     /= 0) c(ind_C5H8)    =  50.E-12 * cair ! new
      IF (ind_MVK      /= 0) c(ind_MVK)     = 130.E-12 * cair ! 4Igen
      !IF (ind_MVK      /= 0) c(ind_MVK)     = 160.E-12 * cair ! 3rdgen
      !IF (ind_MVK      /= 0) c(ind_MVK)     = 200.E-12 * cair ! new
      !IF (ind_MVKOOH   /= 0) c(ind_MVKOOH)  = 0.E-12 * cair ! 4Igentest
      IF (ind_MVKOOH   /= 0) c(ind_MVKOOH)  = 180.E-12 * cair ! 4Igen
      !IF (ind_MVKOOH   /= 0) c(ind_MVKOOH)  = 260.E-12 * cair ! 3rdgen
      !IF (ind_MVKOOH   /= 0) c(ind_MVKOOH)  = 295.E-12 * cair ! new
      IF (ind_NACA     /= 0) c(ind_NACA)    =  1.6E-12 * cair ! 4Igen
      !IF (ind_NACA     /= 0) c(ind_NACA)    =  2.1E-12 * cair ! 4gen
      !IF (ind_NACA     /= 0) c(ind_NACA)    = 1.35E-12 * cair ! 3rdgen
      !IF (ind_NACA     /= 0) c(ind_NACA)    = 1.65E-12 * cair ! new
      IF (ind_NH3      /= 0) c(ind_NH3)     = 170.E-12 * cair ! 4Igen
      !IF (ind_NH3      /= 0) c(ind_NH3)     = 200.E-12 * cair ! 4gen
      !IF (ind_NH3      /= 0) c(ind_NH3)     = 150.E-12 * cair ! 3rdgen
      !IF (ind_NH3      /= 0) c(ind_NH3)     = 250.E-12 * cair ! mod8
      !IF (ind_NO       /= 0) c(ind_NO)      =  10.E-12 * cair ! def(off)
      !IF (ind_NO       /= 0) c(ind_NO)      =  0.3E-12 * cair ! JW (<<Singh)
      !IF (ind_NO2      /= 0) c(ind_NO2)     = 1.0E-12 * cair ! JW (<<Singh)
      IF (ind_NO2      /= 0) c(ind_NO2)     =  2.E-12 * cair ! 4gen
      !IF (ind_NO2      /= 0) c(ind_NO2)     =  1.6E-12 * cair ! mod3
      IF (ind_O3       /= 0) c(ind_O3)      = 10.4E-09 * cair ! 4Igen
      !IF (ind_O3       /= 0) c(ind_O3)      = 10.0E-09 * cair ! 3rdgen
      !IF (ind_O3       /= 0) c(ind_O3)      = 10.6E-09 * cair ! mod6
      IF (ind_CH3CO3H      /= 0) c(ind_CH3CO3H)     = 300.E-12 * cair ! 4Igen
      !IF (ind_CH3CO3H      /= 0) c(ind_CH3CO3H)     = 460.E-12 * cair ! 4gen
      !IF (ind_CH3CO3H      /= 0) c(ind_CH3CO3H)     = 390.E-12 * cair ! 3rdgen
      !IF (ind_CH3CO3H      /= 0) c(ind_CH3CO3H)     = 420.E-12 * cair ! new
      IF (ind_PAN      /= 0) c(ind_PAN)     =   2.E-12 * cair ! 4gen
      !IF (ind_PAN      /= 0) c(ind_PAN)     =   1.E-12 * cair ! new
      IF (ind_SO2      /= 0) c(ind_SO2)     = 135.E-12 * cair ! 4Igen
      !IF (ind_SO2      /= 0) c(ind_SO2)     = 130.E-12 * cair ! 3rdgen
      !IF (ind_SO2      /= 0) c(ind_SO2)     =  35.E-12 * cair ! mod8

      ! NOT TAKEN INTO ACCOUNT EVEN THOUGH PRESENT:
      !IF (ind_C2H2     /= 0) c(ind_C2H2)     = 28.7E-12 * cair ! Singh
      !IF (ind_C2H4     /= 0) c(ind_C2H4)     =  21.E-12 * cair ! Singh
      !IF (ind_C2H6     /= 0) c(ind_C2H6)     = 28.7E-12 * cair ! Singh
      !IF (ind_C3H6     /= 0) c(ind_C3H6)     = 11.4E-12 * cair ! Singh
      !IF (ind_C3H8     /= 0) c(ind_C3H8)     =  16.E-12 * cair ! Singh
      !IF (ind_CH3COCH3 /= 0) c(ind_CH3COCH3) = 300.E-12 * cair ! JW
      !IF (ind_PAN      /= 0) c(ind_PAN)      =   2.E-12 * cair ! Singh
    END SUBROUTINE x0_oomph

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_isoo
      IF (ind_CO2      /= 0) c(ind_CO2)     = 382.E-06 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 780.E-03 * cair
      IF (ind_O2       /= 0) c(ind_O2)      = 210.E-03 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)     = 1.80E-06 * cair

      ! Photochem.smog parcel chemically degrading in the low troposphere
      ! + oxygen isotopes

      PRINT *, 'ISOO init: tropospheric photochem. smog parcel'
      IF (ind_ACETOL   /= 0) c(ind_ACETOL)    = 250.E-12 * cair
      IF (ind_HYPERACET/= 0) c(ind_HYPERACET) = 0.74E-12 * cair ! 4Igen
      IF (ind_CH3CHO   /= 0) c(ind_CH3CHO)    = 100.E-12 * cair ! 4gen
      IF (ind_CH3COCH3 /= 0) c(ind_CH3COCH3)  = 100.E-12 * cair ! 4gen
      IF (ind_MGLYOX   /= 0) c(ind_MGLYOX)    = 230.E-12 * cair ! 4Igen
      IF (ind_CH3OH    /= 0) c(ind_CH3OH)     = 500.E-12 * cair ! 4gen
      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH)    = 590.E-12 * cair ! 3rdgen
      IF (ind_CH3CO2H  /= 0) c(ind_CH3CO2H)   = 1.34E-09 * cair ! 4Igen
      IF (ind_CO       /= 0) c(ind_CO)        = 200.E-09 * cair ! 4Igen
      IF (ind_DMS      /= 0) c(ind_DMS)       = 200.E-12 * cair ! 4gen
      IF (ind_DMSO     /= 0) c(ind_DMSO)      =  18.E-12 * cair ! 4Igen
      IF (ind_H2       /= 0) c(ind_H2)        =   1.E-06 * cair ! def
      IF (ind_H2O2     /= 0) c(ind_H2O2)      = 220.E-12 * cair ! 4Igen
      IF (ind_HCHO     /= 0) c(ind_HCHO)      = 220.E-12 * cair ! 4Igen
      IF (ind_HCOOH    /= 0) c(ind_HCOOH)     =   7.E-12 * cair ! 4Igen
      IF (ind_HNO3     /= 0) c(ind_HNO3)      = 0.15E-12 * cair ! mod3
      IF (ind_NH3      /= 0) c(ind_NH3)       = 170.E-12 * cair ! 4Igen
      IF (ind_NO2      /= 0) c(ind_NO2)       =   2.E-12 * cair ! 4gen
      IF (ind_O3       /= 0) c(ind_O3)        = 10.4E-09 * cair ! 4Igen
      IF (ind_CH3CO3H  /= 0) c(ind_CH3CO3H)   = 300.E-12 * cair ! 4Igen
      IF (ind_PAN      /= 0) c(ind_PAN)       =   2.E-12 * cair ! 4gen
      IF (ind_SO2      /= 0) c(ind_SO2)       = 135.E-12 * cair ! 4Igen
    END SUBROUTINE x0_isoo

    ! ------------------------------------------------------------------------

    ! mz_ab_20091111+
    SUBROUTINE x0_strato20

      ! stratosphere, 20 hPa
      ! from scout02 (ProSECCO simulation)

      IF (ind_H        /= 0) c(ind_H)      =   1.E-12 * cair
      IF (ind_OH       /= 0) c(ind_OH)     =   1.E-16 * cair
      IF (ind_HO2      /= 0) c(ind_HO2)    =   1.E-15 * cair
      IF (ind_N        /= 0) c(ind_N)      =   1.E-12 * cair
      IF (ind_NO3      /= 0) c(ind_NO3)    =   1.E-12 * cair
      IF (ind_N2O5     /= 0) c(ind_N2O5)   =   1.E-10 * cair
      IF (ind_HNO4     /= 0) c(ind_HNO4)   =   1.E-10 * cair
      IF (ind_CL       /= 0) c(ind_CL)     =   1.E-21 * cair
      IF (ind_CLO      /= 0) c(ind_CLO)    =   1.E-15 * cair
      IF (ind_HOCl     /= 0) c(ind_HOCl)   =  40.E-12 * cair
      IF (ind_CL2O2    /= 0) c(ind_CL2O2)  =   1.E-13 * cair
      IF (ind_CL2      /= 0) c(ind_CL2)    =   1.E-13 * cair
      IF (ind_CH3O2    /= 0) c(ind_CH3O2)  =   1.E-12 * cair
      IF (ind_N2O      /= 0) c(ind_N2O)    =  1.3E-07 * cair
      IF (ind_CO       /= 0) c(ind_CO)     =  1.4E-08 * cair
      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH) =  12.E-12 * cair
      IF (ind_ClNO3    /= 0) c(ind_ClNO3)  =   6.E-10 * cair
      IF (ind_CFCl3    /= 0) c(ind_CFCl3)  =  1.4E-11 * cair
      IF (ind_CF2Cl2   /= 0) c(ind_CF2Cl2) =   1.E-12 * cair
      IF (ind_CH3CL    /= 0) c(ind_CH3CL)  =   1.E-12 * cair
      IF (ind_CCL4     /= 0) c(ind_CCL4)   =   1.E-12 * cair
      IF (ind_CH3CCL3  /= 0) c(ind_CH3CCL3)=   1.E-12 * cair
      IF (ind_HNO3     /= 0) c(ind_HNO3)   =   5.E-09 * cair
      IF (ind_H2O      /= 0) c(ind_H2O)    =   1.E-12 * cair
      IF (ind_O3P      /= 0) c(ind_O3P)    =   9.E-34 * cair
      IF (ind_O1D      /= 0) c(ind_O1D)    =   1.E-16 * cair
      IF (ind_H2       /= 0) c(ind_H2)     =   5.E-07 * cair
      IF (ind_O3       /= 0) c(ind_O3)     =   4.E-06 * cair
      IF (ind_NO       /= 0) c(ind_NO)     =   1.E-24 * cair
      IF (ind_NO2      /= 0) c(ind_NO2)    =   1.E-09 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)    =  1.8E-06 * cair
      IF (ind_HCHO     /= 0) c(ind_HCHO)   =   7.E-11 * cair
      IF (ind_CO       /= 0) c(ind_CO)     =  70.E-09 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)    = 350.E-06 * cair
      IF (ind_H2O2     /= 0) c(ind_H2O2)   = 450.E-12 * cair
      IF (ind_HCl      /= 0) c(ind_HCl)    = 400.E-12 * cair
      ! additional for mecca
      IF (ind_SO2      /= 0) c(ind_SO2)    =   0. * cair
      IF (ind_O2       /= 0) c(ind_O2)     = 210.E-03 * cair
      IF (ind_N2       /= 0) c(ind_N2)     = 780.E-03 * cair

    END SUBROUTINE x0_strato20
    ! mz_ab_20091111-

    ! ------------------------------------------------------------------------

    ! mz_ab_20091111+
    SUBROUTINE x0_strato10

      ! stratosphere, 10 hPa
      IF (ind_H        /= 0) c(ind_H)      =   1.E-16 * cair
      IF (ind_OH       /= 0) c(ind_OH)     =   1.E-16 * cair
      IF (ind_HO2      /= 0) c(ind_HO2)    =   1.E-15 * cair
      IF (ind_N        /= 0) c(ind_N)      =   1.E-12 * cair
      IF (ind_NO3      /= 0) c(ind_NO3)    =   1.E-12 * cair
      IF (ind_N2O5     /= 0) c(ind_N2O5)   =   1.E-10 * cair
      IF (ind_HNO4     /= 0) c(ind_HNO4)   =   1.5E-10 * cair
      IF (ind_CL       /= 0) c(ind_CL)     =   1.E-21 * cair
      IF (ind_CLO      /= 0) c(ind_CLO)    =   1.E-15 * cair
      IF (ind_HOCl     /= 0) c(ind_HOCl)   =  40.E-12 * cair
      IF (ind_CL2O2    /= 0) c(ind_CL2O2)  =   1.E-13 * cair
      IF (ind_CL2      /= 0) c(ind_CL2)    =   1.E-13 * cair
      IF (ind_CH3O2    /= 0) c(ind_CH3O2)  =   1.E-12 * cair
      IF (ind_N2O      /= 0) c(ind_N2O)    =  1.3E-07 * cair
      IF (ind_CO       /= 0) c(ind_CO)     =  1.4E-08 * cair
      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH) =  12.E-12 * cair
      IF (ind_ClNO3    /= 0) c(ind_ClNO3)  =   6.E-10 * cair
      IF (ind_CFCl3    /= 0) c(ind_CFCl3)  =  1.4E-11 * cair
      IF (ind_CF2Cl2   /= 0) c(ind_CF2Cl2) =   1.E-12 * cair
      IF (ind_CH3CL    /= 0) c(ind_CH3CL)  =   1.E-12 * cair
      IF (ind_CCL4     /= 0) c(ind_CCL4)   =   1.E-12 * cair
      IF (ind_CH3CCL3  /= 0) c(ind_CH3CCL3)=   1.E-12 * cair
      IF (ind_HNO3     /= 0) c(ind_HNO3)   =   5.E-09 * cair
      IF (ind_H2O      /= 0) c(ind_H2O)    = 4.251E-06 * cair
      IF (ind_O3P      /= 0) c(ind_O3P)    =   9.E-34 * cair
      IF (ind_O1D      /= 0) c(ind_O1D)    =   1.E-16 * cair
      IF (ind_H2       /= 0) c(ind_H2)     =   5.E-07 * cair
      IF (ind_O3       /= 0) c(ind_O3)     =   8.E-06 * cair
      IF (ind_NO       /= 0) c(ind_NO)     =   1.E-24 * cair
      IF (ind_NO2      /= 0) c(ind_NO2)    =   1.E-09 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)    =  1.8E-06 * cair
      IF (ind_HCHO     /= 0) c(ind_HCHO)   =   7.E-11 * cair
      IF (ind_CO       /= 0) c(ind_CO)     =  70.E-09 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)    = 350.E-06 * cair
      IF (ind_H2O2     /= 0) c(ind_H2O2)   = 180.E-12 * cair
      IF (ind_HCl      /= 0) c(ind_HCl)    = 400.E-12 * cair
      ! additional for mecca
      IF (ind_SO2      /= 0) c(ind_SO2)    =   0. * cair
      IF (ind_O2       /= 0) c(ind_O2)     = 210.E-03 * cair
      IF (ind_N2       /= 0) c(ind_N2)     = 780.E-03 * cair

    END SUBROUTINE x0_strato10
    ! mz_ab_20091111-

    ! ------------------------------------------------------------------------

    ! mz_ab_20091111+
    SUBROUTINE x0_mtchem

!!$      ! mesosphere, 0.01 hPa
!!$      IF (ind_H        /= 0) c(ind_H)      = 1.876E-07* cair
!!$      IF (ind_OH       /= 0) c(ind_OH)     = 1.240E-08* cair
!!$      IF (ind_HO2      /= 0) c(ind_HO2)    = 5.012E-09* cair
!!$      IF (ind_N        /= 0) c(ind_N)      = 2.114E-10* cair
!!$      IF (ind_NO3      /= 0) c(ind_NO3)    = 3.083E-21* cair
!!$      IF (ind_N2O5     /= 0) c(ind_N2O5)   = 9.072E-27* cair
!!$      IF (ind_HNO4     /= 0) c(ind_HNO4)   = 4.754E-18* cair
!!$      IF (ind_CL       /= 0) c(ind_CL)     = 8.001E-11* cair
!!$      IF (ind_CLO      /= 0) c(ind_CLO)    = 1.564E-13* cair
!!$      IF (ind_HOCl     /= 0) c(ind_HOCl)   = 7.015E-15* cair
!!$      IF (ind_CL2O2    /= 0) c(ind_CL2O2)  = 1.558E-25* cair
!!$      IF (ind_CL2      /= 0) c(ind_CL2)    = 1.E-13* cair ! ????
!!$      IF (ind_CH3O2    /= 0) c(ind_CH3O2)  = 1.E-12* cair! ????
!!$      IF (ind_N2O      /= 0) c(ind_N2O)    = 7.077E-11* cair
!!$      IF (ind_H2O2     /= 0) c(ind_H2O2)   = 1.685E-12* cair
!!$      IF (ind_CO       /= 0) c(ind_CO)     = 9.862E-07* cair
!!$      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH) = 3.558E-13* cair
!!$      IF (ind_ClNO3    /= 0) c(ind_ClNO3)  = 1.460E-22* cair
!!$      IF (ind_CFCl3    /= 0) c(ind_CFCl3)  = 2.854E-38* cair
!!$      IF (ind_CF2Cl2   /= 0) c(ind_CF2Cl2) = 2.776E-17* cair
!!$      IF (ind_CH3CL    /= 0) c(ind_CH3CL)  = 2.212E-14* cair
!!$      IF (ind_CCL4     /= 0) c(ind_CCL4)   = 5.605E-45* cair
!!$      IF (ind_HNO3     /= 0) c(ind_HNO3)   = 2.084E-13* cair
!!$      IF (ind_H2O      /= 0) c(ind_H2O)    = 4.567E-06* cair
!!$      IF (ind_O3P      /= 0) c(ind_O3P)    = 5.350E-06* cair
!!$      IF (ind_O1D      /= 0) c(ind_O1D)    = 3.143E-14* cair
!!$      IF (ind_H2       /= 0) c(ind_H2)     = 1.310E-06* cair
!!$      IF (ind_O3       /= 0) c(ind_O3)     = 7.823E-08* cair
!!$      IF (ind_NO       /= 0) c(ind_NO)     = 2.454E-09* cair
!!$      IF (ind_NO2      /= 0) c(ind_NO2)    = 1.685E-12* cair
!!$      IF (ind_CH4      /= 0) c(ind_CH4)    = 1.113E-07* cair
!!$      IF (ind_HCHO     /= 0) c(ind_HCHO)   = 1.417E-12* cair
!!$      IF (ind_CO       /= 0) c(ind_CO)     = 9.862E-07* cair
!!$      IF (ind_CO2      /= 0) c(ind_CO2)    = 3.641E-04* cair
!!$      IF (ind_H2O2     /= 0) c(ind_H2O2)   = 1.169E-10* cair
!!$      IF (ind_HCl      /= 0) c(ind_HCl)    = 3.342E-09* cair
!!$      ! additional for mecca
!!$      IF (ind_SO2      /= 0) c(ind_SO2)    =   0. * cair
!!$      IF (ind_O2       /= 0) c(ind_O2)     = 210.E-03 * cair
!!$      IF (ind_N2       /= 0) c(ind_N2)     = 780.E-03 * cair
!!$      IF (ind_em       /= 0) c(ind_em)     = 1.2E-14 * cair
!!$      IF (ind_NOp      /= 0) c(ind_NOp)     = 1.2E-10 * cair
!!$      IF (ind_O2p      /= 0) c(ind_O2p)     = 1.2E-10 * cair
!!$      IF (ind_Op      /= 0) c(ind_Op)     = 1.2E-10 * cair
!!$      IF (ind_Np      /= 0) c(ind_Np)     = 1.2E-10 * cair
!!$      IF (ind_H2O      /= 0) c(ind_H2O)     = 1.2E-7 * cair

      ! from Miriam Sinnhuber, 0.51 Pa
      IF (ind_H        /= 0) c(ind_H)      = 7.69140E-07* cair
      IF (ind_OH       /= 0) c(ind_OH)     = 6.53740E-10* cair
      IF (ind_HO2      /= 0) c(ind_HO2)    = 2.09240E-10* cair
      IF (ind_N        /= 0) c(ind_N)      = 2.114E-10* cair
      IF (ind_NO3      /= 0) c(ind_NO3)    = 2.07260E-24* cair
      IF (ind_N2O5     /= 0) c(ind_N2O5)   = 8.59140E-30* cair
      IF (ind_HNO4     /= 0) c(ind_HNO4)   = 4.754E-18* cair
      IF (ind_CL       /= 0) c(ind_CL)     = 8.09210E-10* cair
      IF (ind_CLO      /= 0) c(ind_CLO)    = 1.28810E-13* cair
      IF (ind_HOCl     /= 0) c(ind_HOCl)   = 8.07360E-18* cair
      IF (ind_CL2O2    /= 0) c(ind_CL2O2)  = 4.67980E-27* cair
      IF (ind_CL2      /= 0) c(ind_CL2)    = 1.E-13* cair ! ????
      IF (ind_CH3O2    /= 0) c(ind_CH3O2)  = 1.E-12* cair! ????
      IF (ind_N2O      /= 0) c(ind_N2O)    = 7.84100E-10* cair
      IF (ind_H2O2     /= 0) c(ind_H2O2)   = 4.64940E-13* cair
      IF (ind_CO       /= 0) c(ind_CO)     = 3.11580E-05* cair
      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH) = 3.558E-13* cair
      IF (ind_ClNO3    /= 0) c(ind_ClNO3)  = 3.84100E-26* cair
      IF (ind_CFCl3    /= 0) c(ind_CFCl3)  = 2.854E-38* cair
      IF (ind_CF2Cl2   /= 0) c(ind_CF2Cl2) = 2.776E-17* cair
      IF (ind_CH3CL    /= 0) c(ind_CH3CL)  = 2.212E-14* cair
      IF (ind_CCL4     /= 0) c(ind_CCL4)   = 5.605E-35* cair
      IF (ind_HNO3     /= 0) c(ind_HNO3)   = 1.14160E-18* cair
      IF (ind_O3P      /= 0) c(ind_O3P)    = 1.18100E-04* cair
      IF (ind_O1D      /= 0) c(ind_O1D)    = 9.87100E-14* cair
      IF (ind_H2       /= 0) c(ind_H2)     = 5.29290E-06* cair
      IF (ind_O3       /= 0) c(ind_O3)     = 8.45990E-08* cair
      IF (ind_NO       /= 0) c(ind_NO)     = 2.05670E-08* cair
      IF (ind_NO2      /= 0) c(ind_NO2)    = 4.20810E-14* cair
      IF (ind_CH4      /= 0) c(ind_CH4)    = 2.63900E-09* cair
      IF (ind_HCHO     /= 0) c(ind_HCHO)   = 1.417E-12* cair
      IF (ind_CO       /= 0) c(ind_CO)     = 3.11580E-05 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)    = 3.54360E-04* cair
      IF (ind_H2O2     /= 0) c(ind_H2O2)   = 1.169E-10* cair
      IF (ind_HCl      /= 0) c(ind_HCl)    = 2.39140E-09* cair
      IF (ind_SO2      /= 0) c(ind_SO2)    =   0. * cair
      IF (ind_O2       /= 0) c(ind_O2)     = 210.E-03 * cair
      IF (ind_N2       /= 0) c(ind_N2)     = 780.E-03 * cair
!!$      IF (ind_em       /= 0) c(ind_em)     = 1.2E-14 * cair
!!$      IF (ind_NOp      /= 0) c(ind_NOp)     = 1.2E-10 * cair
!!$      IF (ind_O2p      /= 0) c(ind_O2p)     = 1.2E-10 * cair
!!$      IF (ind_Op      /= 0) c(ind_Op)     = 1.2E-10 * cair
!!$      IF (ind_Np      /= 0) c(ind_Np)     = 1.2E-10 * cair
      IF (ind_H2O      /= 0) c(ind_H2O)     = 3.36210E-09 * cair
      IF (ind_N       /= 0) c(ind_N)     =9.27150E-10

    END SUBROUTINE x0_mtchem


	    SUBROUTINE x0_radmka
      ! ka_rr_20150807+ 
      IF (ind_H2O2     /= 0) c(ind_H2O2)    =   7.E-09 * cair
      IF (ind_O3       /= 0) c(ind_O3)      =  30.E-09 * cair
      IF (ind_O2       /= 0) c(ind_O2)      = 210.E-03 * cair
      IF (ind_NH3      /= 0) c(ind_NH3)     = 100.E-12 * cair
      IF (ind_NO       /= 0) c(ind_NO)      =  10.E-12 * cair
      IF (ind_NO2      /= 0) c(ind_NO2)     = 100.E-12 * cair
      IF (ind_HONO     /= 0) c(ind_HONO)    =  40.E-14 * cair
      IF (ind_HNO3     /= 0) c(ind_HNO3)    =  5.0E-12 * cair
      IF (ind_N2       /= 0) c(ind_N2)      = 780.E-03 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)     =  1.8E-06 * cair
      IF (ind_HCHO     /= 0) c(ind_HCHO)    =  5.0E-09 * cair
      IF (ind_CO       /= 0) c(ind_CO)      = 100.E-09 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)     = 350.E-06 * cair
     ! IF (ind_ISO      /= 0) c(ind_ISO)     =  2.0E-09 * cair
     ! IF (ind_PAN_Ka   /= 0) c(ind_PAN_Ka)  = 100.E-12 * cair
      !IF (ind_ISO      /= 0) c(ind_ISO)     =  2.0E-09 * cair
      !IF (ind_PAN_Ka   /= 0) c(ind_PAN_Ka)  = 100.E-12 * cair
      IF (ind_OH       /= 0) c(ind_OH)     =   1.E-16 * cair
      IF (ind_N2O         /= 0) c(ind_N2O)         =  3.16661E-07 * cair
      IF (ind_CH3O2       /= 0) c(ind_CH3O2)       =  1.1287E-12  * cair
      IF (ind_CH3COCH3    /= 0) c(ind_CH3COCH3)    =  3.34549E-10 * cair
      
      
      
      
      ! ka_rr_20150807- 
    END SUBROUTINE x0_radmka


        SUBROUTINE x0_ccmi

      IF (ind_CF3Br       /= 0) c(ind_CF3Br)       =  2.43175E-12 * cair
      IF (ind_CF2ClBr     /= 0) c(ind_CF2ClBr)     =  4.20382E-12 * cair
      
      IF (ind_H        /= 0) c(ind_H)      =   1.E-16 * cair
      !IF (ind_OH       /= 0) c(ind_OH)     =   1.E-16 * cair
      IF (ind_OH       /= 0) c(ind_OH)     =   1.E-15 * cair
      IF (ind_HO2      /= 0) c(ind_HO2)    =   1.E-15 * cair
      IF (ind_N        /= 0) c(ind_N)      =   1.E-12 * cair
      IF (ind_NO3      /= 0) c(ind_NO3)    =   1.E-12 * cair
      IF (ind_N2O5     /= 0) c(ind_N2O5)   =   1.E-10 * cair
      IF (ind_HNO4     /= 0) c(ind_HNO4)   =   1.5E-10 * cair
      IF (ind_CL       /= 0) c(ind_CL)     =   1.E-21 * cair
      IF (ind_CLO      /= 0) c(ind_CLO)    =   1.E-15 * cair
      IF (ind_HOCl     /= 0) c(ind_HOCl)   =  40.E-12 * cair
      IF (ind_CL2O2    /= 0) c(ind_CL2O2)  =   1.E-13 * cair
      IF (ind_CL2      /= 0) c(ind_CL2)    =   1.E-13 * cair
      IF (ind_CH3O2    /= 0) c(ind_CH3O2)  =   1.E-12 * cair
      IF (ind_N2O      /= 0) c(ind_N2O)    =  1.3E-07 * cair
      IF (ind_CO       /= 0) c(ind_CO)     =  1.4E-08 * cair
      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH) =  12.E-12 * cair
      IF (ind_ClNO3    /= 0) c(ind_ClNO3)  =   6.E-10 * cair
      IF (ind_CFCl3    /= 0) c(ind_CFCl3)  =  1.4E-11 * cair
      IF (ind_CF2Cl2   /= 0) c(ind_CF2Cl2) =   1.E-12 * cair
      IF (ind_CH3CL    /= 0) c(ind_CH3CL)  =   1.E-12 * cair
      IF (ind_CCL4     /= 0) c(ind_CCL4)   =   1.E-12 * cair
      IF (ind_CH3CCL3  /= 0) c(ind_CH3CCL3)=   1.E-12 * cair
      IF (ind_HNO3     /= 0) c(ind_HNO3)   =   5.E-09 * cair
      IF (ind_H2O      /= 0) c(ind_H2O)    = 4.251E-06 * cair
      IF (ind_O3P      /= 0) c(ind_O3P)    =   9.E-34 * cair
      IF (ind_O1D      /= 0) c(ind_O1D)    =   1.E-16 * cair
      IF (ind_H2       /= 0) c(ind_H2)     =   5.E-07 * cair
      IF (ind_O3       /= 0) c(ind_O3)     =   8.E-06 * cair
      IF (ind_NO       /= 0) c(ind_NO)     =   1.E-24 * cair
      !IF (ind_NO2      /= 0) c(ind_NO2)    =   1.E-09 * cair
      IF (ind_NO2      /= 0) c(ind_NO2)    =   3.E-09 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)    =  1.8E-06 * cair
      IF (ind_HCHO     /= 0) c(ind_HCHO)   =   7.E-11 * cair
      IF (ind_CO       /= 0) c(ind_CO)     =  70.E-09 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)    = 350.E-06 * cair
      IF (ind_H2O2     /= 0) c(ind_H2O2)   = 180.E-12 * cair
      IF (ind_HCl      /= 0) c(ind_HCl)    = 400.E-12 * cair
      ! additional for mecca
      IF (ind_SO2      /= 0) c(ind_SO2)    =   0. * cair
      IF (ind_O2       /= 0) c(ind_O2)     = 210.E-03 * cair
      IF (ind_N2       /= 0) c(ind_N2)     = 780.E-03 * cair 


      END SUBROUTINE
    ! mz_ab_20091111-
        SUBROUTINE x0_psc

      IF (ind_H        /= 0) c(ind_H)      =   1.E-16 * cair
      !IF (ind_OH       /= 0) c(ind_OH)     =   1.E-16 * cair
      IF (ind_OH       /= 0) c(ind_OH)     =   1.E-15 * cair
      IF (ind_HO2      /= 0) c(ind_HO2)    =   1.E-15 * cair
      IF (ind_N        /= 0) c(ind_N)      =   1.E-12 * cair
      IF (ind_N2O5     /= 0) c(ind_N2O5)   =   1.E-10 * cair
      IF (ind_CL       /= 0) c(ind_CL)     =      0.0 * cair
      IF (ind_CLO      /= 0) c(ind_CLO)    =      0.0 * cair
      IF (ind_HOCl     /= 0) c(ind_HOCl)   =      0.0 * cair
      IF (ind_CL2O2    /= 0) c(ind_CL2O2)  =      0.0 * cair
      IF (ind_CL2      /= 0) c(ind_CL2)    =      0.0 * cair
      IF (ind_N2O      /= 0) c(ind_N2O)    =  1.3E-07 * cair
      IF (ind_CO       /= 0) c(ind_CO)     =  1.4E-08 * cair
      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH) =  12.E-12 * cair
      IF (ind_ClNO3    /= 0) c(ind_ClNO3)  =1500.E-12 * cair
      IF (ind_CFCl3    /= 0) c(ind_CFCl3)  =  1.4E-11 * cair
      IF (ind_HNO3     /= 0) c(ind_HNO3)   =   5.E-09 * cair
      IF (ind_H2O      /= 0) c(ind_H2O)    = 4.251E-06 * cair
      IF (ind_O3P      /= 0) c(ind_O3P)    =   9.E-34 * cair
      IF (ind_O1D      /= 0) c(ind_O1D)    =   1.E-16 * cair
      IF (ind_H2       /= 0) c(ind_H2)     =   5.E-07 * cair
      IF (ind_O3       /= 0) c(ind_O3)     =   8.E-06 * cair
      IF (ind_NO       /= 0) c(ind_NO)     =   1.E-24 * cair
      !IF (ind_NO2      /= 0) c(ind_NO2)    =   1.E-09 * cair
      IF (ind_NO2      /= 0) c(ind_NO2)    =   3.E-09 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)    =  1.8E-06 * cair
      IF (ind_CO       /= 0) c(ind_CO)     =  70.E-09 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)    = 350.E-06 * cair
      IF (ind_HCl      /= 0) c(ind_HCl)    = 1500.E-12 * cair
      IF (ind_OCS      /= 0) c(ind_OCS)    = 5.E-11 * cair
      ! additional for mecca
      IF (ind_DMS      /= 0) c(ind_DMS)    =  10.E-12 * cair
      IF (ind_SO2      /= 0) c(ind_SO2)    =   0. * cair
      IF (ind_H2SO4    /= 0) c(ind_H2SO4)  =  70.E-12 * cair
      IF (ind_O2       /= 0) c(ind_O2)     = 210.E-03 * cair
      IF (ind_N2       /= 0) c(ind_N2)     = 780.E-03 * cair 
      IF (ind_NH3      /= 0) c(ind_NH3)    = 100.E-12 * cair


      END SUBROUTINE x0_psc

    ! ------------------------------------------------------------------------

  END SUBROUTINE x0
END MODULE mo_mecicon_init
