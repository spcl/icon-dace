! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! The Reaction Rates File
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

MODULE messy_mecca_kpp_rates 
  USE mo_kind,                 ONLY: dp

  USE messy_mecca_kpp_Parameters
  USE messy_mecca_kpp_global
  IMPLICIT NONE 
  PUBLIC :: Update_RCONST, Update_PHOTO

CONTAINS



! Begin Rate Law Functions from KPP_HOME/util/UserRateLaws

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  User-defined Rate Law functions
!  Note: the default argument type for rate laws, as read from the equations file, is single precision
!        but all the internal calculations are performed in double precision
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~>  Arrhenius
   REAL(kind=dp) FUNCTION ARR( A0,B0,C0 )
      REAL A0,B0,C0      
      ARR =  DBLE(A0) * EXP(-DBLE(B0)/TEMP) * (TEMP/300.0_dp)**DBLE(C0)
   END FUNCTION ARR        

!~~~> Simplified Arrhenius, with two arguments
!~~~> Note: The argument B0 has a changed sign when compared to ARR
   REAL(kind=dp) FUNCTION ARR2( A0,B0 )
      REAL A0,B0           
      ARR2 =  DBLE(A0) * EXP( DBLE(B0)/TEMP )              
   END FUNCTION ARR2          

   REAL(kind=dp) FUNCTION EP2(A0,C0,A2,C2,A3,C3)
      REAL A0,C0,A2,C2,A3,C3
      REAL(kind=dp) K0,K2,K3            
      K0 = DBLE(A0) * EXP(-DBLE(C0)/TEMP)
      K2 = DBLE(A2) * EXP(-DBLE(C2)/TEMP)
      K3 = DBLE(A3) * EXP(-DBLE(C3)/TEMP)
      K3 = K3*CFACTOR*1.0E6_dp
      EP2 = K0 + K3/(1.0_dp+K3/K2 )
   END FUNCTION EP2

   REAL(kind=dp) FUNCTION EP3(A1,C1,A2,C2) 
      REAL A1, C1, A2, C2
      REAL(kind=dp) K1, K2      
      K1 = DBLE(A1) * EXP(-DBLE(C1)/TEMP)
      K2 = DBLE(A2) * EXP(-DBLE(C2)/TEMP)
      EP3 = K1 + K2*(1.0E6_dp*CFACTOR)
   END FUNCTION EP3 

   REAL(kind=dp) FUNCTION FALL ( A0,B0,C0,A1,B1,C1,CF)
      REAL A0,B0,C0,A1,B1,C1,CF
      REAL(kind=dp) K0, K1     
      K0 = DBLE(A0) * EXP(-DBLE(B0)/TEMP)* (TEMP/300.0_dp)**DBLE(C0)
      K1 = DBLE(A1) * EXP(-DBLE(B1)/TEMP)* (TEMP/300.0_dp)**DBLE(C1)
      K0 = K0*CFACTOR*1.0E6_dp
      K1 = K0/K1
      FALL = (K0/(1.0_dp+K1))*   &
           DBLE(CF)**(1.0_dp/(1.0_dp+(LOG10(K1))**2))
   END FUNCTION FALL

  !---------------------------------------------------------------------------

  ELEMENTAL REAL(kind=dp) FUNCTION k_3rd(temp,cair,k0_300K,n,kinf_300K,m,fc)

    INTRINSIC LOG10

    REAL(kind=dp), INTENT(IN) :: temp      ! temperature [K]
    REAL(kind=dp), INTENT(IN) :: cair      ! air concentration [molecules/cm3]
    REAL, INTENT(IN) :: k0_300K   ! low pressure limit at 300 K
    REAL, INTENT(IN) :: n         ! exponent for low pressure limit
    REAL, INTENT(IN) :: kinf_300K ! high pressure limit at 300 K
    REAL, INTENT(IN) :: m         ! exponent for high pressure limit
    REAL, INTENT(IN) :: fc        ! broadening factor (usually fc=0.6)
    REAL(kind=dp) :: zt_help, k0_T, kinf_T, k_ratio

    zt_help = 300._dp/temp
    k0_T    = k0_300K   * zt_help**(n) * cair ! k_0   at current T
    kinf_T  = kinf_300K * zt_help**(m)        ! k_inf at current T
    k_ratio = k0_T/kinf_T
    k_3rd   = k0_T/(1._dp+k_ratio)*fc**(1._dp/(1._dp+LOG10(k_ratio)**2))

  END FUNCTION k_3rd

  !---------------------------------------------------------------------------

  ELEMENTAL REAL(kind=dp) FUNCTION k_arr (k_298,tdep,temp)
    ! Arrhenius function

    REAL,     INTENT(IN) :: k_298 ! k at T = 298.15K
    REAL,     INTENT(IN) :: tdep  ! temperature dependence
    REAL(kind=dp), INTENT(IN) :: temp  ! temperature

    INTRINSIC EXP

    k_arr = k_298 * EXP(tdep*(1._dp/temp-3.3540E-3_dp)) ! 1/298.15=3.3540e-3

  END FUNCTION k_arr

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  End of User-defined Rate Law functions
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! End Rate Law Functions from KPP_HOME/util/UserRateLaws


! Begin INLINED Rate Law Functions


  ELEMENTAL REAL(dp) FUNCTION k_SIV_H2O2 (k_298,tdep,cHp,temp)
    ! special rate function for S(IV) + H2O2
    REAL,     INTENT(IN) :: k_298 ! k at T = 298.15K
    REAL,     INTENT(IN) :: tdep  ! temperature dependence
    REAL(dp), INTENT(IN) :: cHp   ! c(H+)
    REAL(dp), INTENT(IN) :: temp  ! temperature
    INTRINSIC :: EXP
    k_SIV_H2O2 = k_298 * EXP(tdep*(1._dp/temp-3.3540E-3_dp)) &
      * cHp / (cHp+0.1_dp)
  END FUNCTION k_SIV_H2O2
  ELEMENTAL REAL(dp) FUNCTION k_3rd_iupac(temp,cair,k0_300K,n,kinf_300K,m,fc)
    ! IUPAC three body reaction formula (iupac.pole-ether.fr)
    INTRINSIC :: LOG10
    REAL(dp), INTENT(IN) :: temp      ! temperature [K]
    REAL(dp), INTENT(IN) :: cair      ! air concentration [molecules/cm3]
    REAL,     INTENT(IN) :: k0_300K   ! low pressure limit at 300 K
    REAL,     INTENT(IN) :: n         ! exponent for low pressure limit
    REAL,     INTENT(IN) :: kinf_300K ! high pressure limit at 300 K
    REAL,     INTENT(IN) :: m         ! exponent for high pressure limit
    REAL,     INTENT(IN) :: fc        ! broadening factor (e.g. 0.45 or 0.6...)
    REAL                 :: nu        ! N
    REAL                 :: zt_help, k0_T, kinf_T, k_ratio
    zt_help = 300._dp/temp
    k0_T    = k0_300K   * zt_help**(n) * cair ! k_0   at current T
    kinf_T  = kinf_300K * zt_help**(m)        ! k_inf at current T
    k_ratio = k0_T/kinf_T
    nu      = 0.75-1.27*LOG10(fc)
    k_3rd_iupac = k0_T/(1._dp+k_ratio)* &
      fc**(1._dp/(1._dp+(LOG10(k_ratio)/nu)**2))
  END FUNCTION k_3rd_iupac
! mz_dt_20141104+
  ELEMENTAL REAL(dp) FUNCTION alpha_AN(n,ro2type,temp,cair)
  ! Alkyl nitrate yields dependent on T and P according to
  ! Arey et al., doi:10.1021/jp003292z
  ! Teng et al., ACPD 2014 doi:10.5194/acpd-14-6721-2014
    INTRINSIC :: LOG10
    INTEGER,  INTENT(IN) :: n         ! number of heavy atoms (C, O, N) excluding the O-O moiety
    INTEGER,  INTENT(IN) :: ro2type   ! 1, 2 or 3 for primary, secondary and tertiary RO2
    REAL(dp), INTENT(IN) :: temp      ! temperature [K]
    REAL(dp), INTENT(IN) :: cair      ! air concentration [molecules/cm3]
    REAL, PARAMETER      :: alpha=2.E-22, beta=1.0, Yinf_298K=0.43, &
                            F=0.41, m0=0., minf=8.0
    REAL                 :: m ! factor for primary, secondary and tertiary RO2
    REAL                 :: Y0_298K, Y0_298K_tp, Yinf_298K_t, zeta, k_ratio
! According to Teng et al. the nature of the radical does not affect the nitrate yield
!!    IF (ro2type = 1) THEN
!!      m = 0.4              ! primary RO2
!!    ELSE IF (ro2type = 2) THEN
!!      m = 1.               ! secondary RO2
!!    ELSE IF (ro2type = 3) THEN
!!      m = 0.3              ! tertiary RO2
!!    ELSE
!!      m = 1.
!!    ENDIF
        m = 1.
    Y0_298K     = alpha*EXP(beta*n)
    Y0_298K_tp  = Y0_298K * cair * (temp/298)**(-m0)
    Yinf_298K_t = Yinf_298K * (temp/298)**(-minf)
    zeta        = 1/(1+LOG10(Y0_298K_tp/Yinf_298K_t)**2)
    k_ratio     = (Y0_298K_tp/(1+Y0_298K_tp/Yinf_298K_t))*F**zeta
    alpha_AN    = k_ratio/(1+k_ratio) * m
  END FUNCTION alpha_AN
! mz_dt_20141104+
  ELEMENTAL REAL(dp) FUNCTION k_limited (k3rd,cHp)
    ! diffusion limitation caps 3rd order rate coefficients
    REAL(dp), INTENT(IN) :: k3rd  ! 3rd order rate coefficient
    REAL(dp), INTENT(IN) :: cHp   ! c(H+)
    REAL(dp), PARAMETER  :: DiffLimit = 1E10 ! diffusion limitation [M-1s-1]
    INTRINSIC :: EXP
    k_limited = 1._dp / ( 1._dp/k3rd + cHp/DiffLimit)
  END FUNCTION k_limited
  ELEMENTAL REAL(dp) FUNCTION k_N2_O(temp,temp_ion)
    REAL(dp), INTENT(IN) :: temp      ! temperature [K]
    REAL(dp), INTENT(IN) :: temp_ion  ! ion temperature [K]
    REAL                 :: temp_mean
    temp_mean = (temp_ion + temp)/2
    k_N2_O = 1.4E-10*(300./temp_mean)**0.44
  END FUNCTION k_N2_O
  ELEMENTAL REAL(dp) FUNCTION k_Op_O2(temp,temp_ion)
    REAL(dp), INTENT(IN) :: temp      ! temperature [K]
    REAL(dp), INTENT(IN) :: temp_ion  ! ion temperature [K]
    REAL(dp)             :: temp_mean
    temp_mean  = .667*temp_ion + .333*temp
    k_Op_O2  = 2.82E-11 - 7.74E-12*(temp_mean/300.) + &
      1.073E-12*(temp_mean/300.)**2 - 5.17E-14*(temp_mean/300.)**3 + &
      9.65E-16*(temp_mean/300.)**4
  END FUNCTION k_Op_O2
  ELEMENTAL REAL(dp) FUNCTION k_Op_N2(temp,temp_ion)
    REAL(dp), INTENT(IN) :: temp      ! temperature [K]
    REAL(dp), INTENT(IN) :: temp_ion  ! ion temperature [K]
    REAL(dp)             :: temp_mean
    temp_mean = .6363*temp_ion + .3637*temp
    k_Op_N2 = 1.533E-12 - 5.92E-13*(temp_mean/300.) + 8.6E-14*(temp_mean/300.)**2
  END FUNCTION k_Op_N2
!C
!C*********************************************************************
!C Die Funktion TROE bestimmt die Werte der Troefunktion
!C nach Stockwell et al. [2],
!C die die Druckabhaengigkeit von Reaktionskonstanten durch eine
!C quasi-bimolekulare Reaktionskonstante beschreibt.
!C*********************************************************************
!C
 FUNCTION TROE( K0300, Q, KU300, R, M, T )
!C CALCULATION OF RATE CONSTANTS FOR TROE REACTIONS
!      IMPLICIT REAL(dp) (A-Z)
      REAL(dp) :: K0300, Q, KU300, R, M, T, TROE
      REAL(dp) :: tt, k0, ku, k0m, kk, lgkk, e, f
      TT= T / 3.D2
      K0= K0300 / TT**Q
      KU= KU300 / TT**R
      K0M= K0 * M
      KK= K0M / KU
      LGKK=0.434294481D0 * LOG(KK)
      E=1.D0 / ( 1.D0 + LGKK*LGKK )
      F=0.6D0 ** E
      TROE= F * K0M / ( 1.D0 + KK )
END FUNCTION TROE
!C

! End INLINED Rate Law Functions

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Update_SUN - update SUN light using TIME
!   Arguments :
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE Update_SUN()
      !USE messy_mecca_kpp_Parameters
      !USE messy_mecca_kpp_global

    IMPLICIT NONE

    REAL(kind=dp) :: SunRise, SunSet
    REAL(kind=dp) :: Thour, Tlocal, Ttmp 
    ! PI - Value of pi
    REAL(kind=dp), PARAMETER :: PI = 3.14159265358979d0
    
    SunRise = 4.5_dp 
    SunSet  = 19.5_dp 
    Thour = TIME/3600.0_dp 
    Tlocal = Thour - (INT(Thour)/24)*24

    IF ((Tlocal>=SunRise).AND.(Tlocal<=SunSet)) THEN
       Ttmp = (2.0*Tlocal-SunRise-SunSet)/(SunSet-SunRise)
       IF (Ttmp.GT.0) THEN
          Ttmp =  Ttmp*Ttmp
       ELSE
          Ttmp = -Ttmp*Ttmp
       END IF
       SUN = ( 1.0_dp + COS(PI*Ttmp) )/2.0_dp 
    ELSE
       SUN = 0.0_dp 
    END IF

 END SUBROUTINE Update_SUN

! End of Update_SUN function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Update_RCONST - function to update rate constants
!   Arguments :
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Update_RCONST ( )




! Begin INLINED RCONST


  USE mo_art_mecicon_data ! atm2Pa, N_A, R_gas
  USE messy_cmn_photol_mem     ! IP_MAX, ip_*, jname
  ! end of USE statements
  !
  ! start of executable statements
  ! define some rate constants here if the expressions are too long
  ! for kpp or if they are used more than once
  ! mz_rs_20140904+ only for mim1
  !op_kg_20110805+ for HO2 + NO -> HNO3 (butkovskaya_*.rpl)
  alpha_NO_HO2 = C(ind_H2O)*6.6E-27*temp*EXP(3700./temp)
  beta_NO_HO2  = MAX(((530./temp)+(press*4.8004E-6)-1.73)*0.01,0._dp) 
  k0_NO_HO2    = 3.5E-12*EXP(250./temp)
  !without humidity correction
  k2d_NO_HO2   = (beta_NO_HO2*k0_NO_HO2)/(1.+beta_NO_HO2)
  k1d_NO_HO2   = k0_NO_HO2 - k2d_NO_HO2
  !with humidity correction
  k2w_NO_HO2   = (beta_NO_HO2*k0_NO_HO2*(1.+42.*alpha_NO_HO2))/ &
                 ((1.+alpha_NO_HO2)*(1.+beta_NO_HO2))
  k1w_NO_HO2   = k0_NO_HO2 - k2w_NO_HO2
  !op_kg_20110805-
  k_PrO2_HO2   = 1.9E-13*EXP(1300./temp)
  k_PrO2_NO    = 2.7E-12*EXP(360./temp)
  k_PrO2_CH3O2 = 9.46E-14*EXP(431./temp)
  ! mz_rs_20140904- only for mim1
  k_HO2_HO2    = (1.5E-12*EXP(19./temp)+1.7E-33*EXP(1000./temp)*cair)* & ! {1.5+/-0.2}
                 (1.+1.4E-21*EXP(2200./temp)*C(ind_H2O))
  k_NO3_NO2    = k_3rd(temp,cair,2.0E-30,4.4,1.4E-12,0.7,0.6)
  k_NO2_HO2    = k_3rd(temp,cair,2.0E-31,3.4,2.9E-12,1.1,0.6)
  k_HNO3_OH    = 2.4E-14*EXP(460./temp) + 1./ &
                 ( 1./(6.5E-34*EXP(1335./temp)*cair) + &
                 1./(2.7E-17*EXP(2199./temp)) )
  k_CH3O2      = 1.03E-13*EXP(365./temp) ! self-reaction of CH3O2
  k_CH3OOH_OH  = 5.3E-12*EXP(190./temp)
  k_CH3CO3_NO2 = k_3rd(temp,cair,9.7E-29,5.6,9.3E-12,1.5,0.6)
  k_PAN_M      = k_CH3CO3_NO2/(9.0E-29*EXP(14000./temp))
  k_ClO_ClO    = k_3rd_iupac(temp,cair,2.0E-32, &
                 4.0,1.0E-11,0.0,0.45)
  ! JPL: k_ClO_ClO   = k_3rd(temp,cair,1.6E-32,4.5,2.0E-12,2.4,0.6)
  k_BrO_NO2    = k_3rd_iupac(temp,cair,4.7E-31,3.1,1.8E-11,0.0,0.4)
  ! JPL: k_BrO_NO2   = k_3rd(temp,cair,5.2E-31,3.2,6.9E-12,2.9,0.6)
  k_I_NO2      = k_3rd_iupac(temp,cair,3.0E-31,1.0,6.6E-11,0.0,0.63)
  ! for numerical reasons, the expression is multiplied by 1e30/1e30
  k_DMS_OH     = 1.E-9*EXP(5820./temp)*C(ind_O2)/ &
                 (1.E30+5.*EXP(6280./temp)*C(ind_O2))
  G7402a_yield = 0.8/1.1 ! 0.8+/-0.2
  ! mz_pj_20080812+
  k_O3s = (1.7E-12*EXP(-940./temp)) * C(ind_OH) &  ! <G2104>
        + (1.E-14*EXP(-490./temp)) * C(ind_HO2) &  ! <G2107>
        + jx(ip_O1D) * 2.2E-10 * C(ind_H2O) &      !
        / ( 3.2E-11*EXP(70./temp)*C(ind_O2)   &
           + 1.8E-11*EXP(110./temp)*C(ind_N2) &
           + 2.2E-10*C(ind_H2O) )
  ! mz_pj_20080812-
  ! mz_dt_20140617+
  ! branching ratio for methyl nitrate according to Butkovskaya et al., JPC A 2012
  ! fit of data according to Lindemann-Hinshelwood scheme by J. Crowley
  !http://iupac.pole-ether.fr/datasheets/pdf/ROO_1_CH3O2_NO.pdf
  beta_null_CH3NO3 = 0.00295 + 5.15E-22*cair*(temp/298)**7.4
  beta_inf_CH3NO3  = 0.022
  beta_CH3NO3      = (beta_null_CH3NO3 * beta_inf_CH3NO3) /&
                     (beta_null_CH3NO3 + beta_inf_CH3NO3)
  k_NO2_CH3O2      = k_3rd(temp,cair,1.0E-30,4.8,7.2E-12,2.1,0.6)
  k_G4138          = 4.25E-12 ! average of the two estimates for k(CH2OO + NO2), 7.E-12 and 1.5E-12 
  k_G9408          = 3.66E-11 ! average of the two estimates for k(CH2OO + SO2), 3.9E-11 + 3.42E-11
  ! mz_dt_20140617-
  ! rate coefficients from the MCM:
  KRO2NO  = 2.54E-12*EXP(360./temp)
  KRO2HO2 = 2.91E-13*EXP(1300./temp)
  KAPHO2  = 4.30E-13*EXP(1040./temp) ! CH3CO3 + HO2
  KAPNO   = 8.10E-12*EXP(270./temp)  ! CH3CO3 + NO
  KRO2NO3 = 2.50E-12
  KNO3AL  = 1.4E-12*EXP(-1900./temp)
  J_IC3H7NO3 = 3.7*jx(ip_PAN)
  J_ACETOL   = 0.65*0.11*jx(ip_CHOH)
  RO2 = 0.
  IF (ind_LISOPACO2>0)  RO2 = RO2 + C(ind_LISOPACO2)
  IF (ind_ISOPBO2>0)    RO2 = RO2 + C(ind_ISOPBO2)
  IF (ind_ISOPDO2>0)    RO2 = RO2 + C(ind_ISOPDO2)
  IF (ind_NISOPO2>0)    RO2 = RO2 + C(ind_NISOPO2)
  IF (ind_LHC4ACCO3>0)  RO2 = RO2 + C(ind_LHC4ACCO3)
  IF (ind_LC578O2>0)    RO2 = RO2 + C(ind_LC578O2)
  IF (ind_C59O2>0)      RO2 = RO2 + C(ind_C59O2)
  IF (ind_LNISO3>0)     RO2 = RO2 + C(ind_LNISO3)
  IF (ind_CH3O2>0)      RO2 = RO2 + C(ind_CH3O2)
  IF (ind_HOCH2O2>0)    RO2 = RO2 + C(ind_HOCH2O2)
  IF (ind_CH3CO3>0)     RO2 = RO2 + C(ind_CH3CO3)
  IF (ind_C2H5O2>0)     RO2 = RO2 + C(ind_C2H5O2)
  IF (ind_HOCH2CO3>0)   RO2 = RO2 + C(ind_HOCH2CO3)
  IF (ind_HYPROPO2>0)   RO2 = RO2 + C(ind_HYPROPO2)
  IF (ind_HCOCO3>0)     RO2 = RO2 + C(ind_HCOCO3)
  IF (ind_CO2H3CO3>0)   RO2 = RO2 + C(ind_CO2H3CO3)
  IF (ind_LHMVKABO2>0)  RO2 = RO2 + C(ind_LHMVKABO2)
  IF (ind_MACO3>0)      RO2 = RO2 + C(ind_MACO3)
  IF (ind_MACRO2>0)     RO2 = RO2 + C(ind_MACRO2)
  IF (ind_LMVKOHABO2>0) RO2 = RO2 + C(ind_LMVKOHABO2)
  IF (ind_PRONO3BO2>0)  RO2 = RO2 + C(ind_PRONO3BO2)
  IF (ind_HOCH2CH2O2>0) RO2 = RO2 + C(ind_HOCH2CH2O2)
  IF (ind_CH3COCH2O2>0) RO2 = RO2 + C(ind_CH3COCH2O2)
  IF (ind_IC3H7O2>0)    RO2 = RO2 + C(ind_IC3H7O2)
  IF (ind_LC4H9O2>0)    RO2 = RO2 + C(ind_LC4H9O2)
  IF (ind_LMEKO2>0)     RO2 = RO2 + C(ind_LMEKO2)
  ! APINENE-related:
  IF (ind_LAPINABO2> 0)        RO2 = RO2 + C(ind_LAPINABO2)
  IF (ind_C96O2> 0)            RO2 = RO2 + C(ind_C96O2)
  IF (ind_C97O2> 0)            RO2 = RO2 + C(ind_C97O2)
  IF (ind_C98O2> 0)            RO2 = RO2 + C(ind_C98O2)
  IF (ind_C85O2> 0)            RO2 = RO2 + C(ind_C85O2)
  IF (ind_C86O2> 0)            RO2 = RO2 + C(ind_C86O2)
  IF (ind_PINALO2> 0)          RO2 = RO2 + C(ind_PINALO2)
  IF (ind_C96CO3> 0)           RO2 = RO2 + C(ind_C96CO3)
  IF (ind_C89CO3> 0)           RO2 = RO2 + C(ind_C89CO3)
  IF (ind_C85CO3> 0)           RO2 = RO2 + C(ind_C85CO3)
  IF (ind_ROO6R1O2> 0)         RO2 = RO2 + C(ind_ROO6R1O2)
  IF (ind_RO6R1O2> 0)          RO2 = RO2 + C(ind_RO6R1O2)
  IF (ind_OHMENTHEN6ONEO2> 0)  RO2 = RO2 + C(ind_OHMENTHEN6ONEO2)
  IF (ind_C511O2> 0)           RO2 = RO2 + C(ind_C511O2)
  IF (ind_C106O2> 0)           RO2 = RO2 + C(ind_C106O2)
  IF (ind_CO235C6CO3> 0)       RO2 = RO2 + C(ind_CO235C6CO3)
  IF (ind_CHOC3COCO3> 0)       RO2 = RO2 + C(ind_CHOC3COCO3)
  IF (ind_CO235C6O2> 0)        RO2 = RO2 + C(ind_CO235C6O2)
  IF (ind_C716O2> 0)           RO2 = RO2 + C(ind_C716O2)
  IF (ind_C614O2> 0)           RO2 = RO2 + C(ind_C614O2)
  IF (ind_HCOCH2CO3> 0)        RO2 = RO2 + C(ind_HCOCH2CO3)
  IF (ind_BIACETO2> 0)         RO2 = RO2 + C(ind_BIACETO2)
  IF (ind_CO23C4CO3> 0)        RO2 = RO2 + C(ind_CO23C4CO3)
  IF (ind_C109O2> 0)           RO2 = RO2 + C(ind_C109O2)
  IF (ind_C811CO3> 0)          RO2 = RO2 + C(ind_C811CO3)
  IF (ind_C89O2> 0)            RO2 = RO2 + C(ind_C89O2)
  IF (ind_C812O2> 0)           RO2 = RO2 + C(ind_C812O2)
  IF (ind_C813O2> 0)           RO2 = RO2 + C(ind_C813O2)
  IF (ind_C721CO3> 0)          RO2 = RO2 + C(ind_C721CO3)
  IF (ind_C721O2> 0)           RO2 = RO2 + C(ind_C721O2)
  IF (ind_C722O2> 0)           RO2 = RO2 + C(ind_C722O2)
  IF (ind_C44O2> 0)            RO2 = RO2 + C(ind_C44O2)
  IF (ind_C512O2> 0)           RO2 = RO2 + C(ind_C512O2)
  IF (ind_C513O2> 0)           RO2 = RO2 + C(ind_C513O2)
  IF (ind_CHOC3COO2> 0)        RO2 = RO2 + C(ind_CHOC3COO2)
  IF (ind_C312COCO3> 0)        RO2 = RO2 + C(ind_C312COCO3)
  IF (ind_HOC2H4CO3> 0)        RO2 = RO2 + C(ind_HOC2H4CO3)
  IF (ind_LNAPINABO2> 0)       RO2 = RO2 + C(ind_LNAPINABO2)
  IF (ind_C810O2> 0)           RO2 = RO2 + C(ind_C810O2)
  IF (ind_C514O2> 0)           RO2 = RO2 + C(ind_C514O2)
  IF (ind_CHOCOCH2O2> 0)       RO2 = RO2 + C(ind_CHOCOCH2O2)
  ! BPINENE-related:
  IF (ind_ROO6R1O2> 0)         RO2 = RO2 + C(ind_ROO6R1O2)
  IF (ind_ROO6R3O2> 0)         RO2 = RO2 + C(ind_ROO6R3O2)
  IF (ind_RO6R1O2> 0)          RO2 = RO2 + C(ind_RO6R1O2)
  IF (ind_RO6R3O2> 0)          RO2 = RO2 + C(ind_RO6R3O2)
  IF (ind_BPINAO2> 0)          RO2 = RO2 + C(ind_BPINAO2)
  IF (ind_C8BCO2> 0)           RO2 = RO2 + C(ind_C8BCO2)
  IF (ind_NOPINDO2> 0)         RO2 = RO2 + C(ind_NOPINDO2)
  IF (ind_LNBPINABO2> 0)       RO2 = RO2 + C(ind_LNBPINABO2)

! End INLINED RCONST

  RCONST(1) = (3.3E-11*EXP(55./temp))
  RCONST(2) = (6.E-34*((temp/300.)**(-2.4))*cair)
  RCONST(3) = (8.E-12*EXP(-2060./temp))
  RCONST(4) = (1.7E-12*EXP(-940./temp))
  RCONST(5) = (1.E-14*EXP(-490./temp))
  RCONST(6) = (4.8E-11*EXP(250./temp))
  RCONST(7) = (1.63E-10*EXP(60./temp))
  RCONST(8) = (2.15E-11*EXP(110./temp))
  RCONST(9) = (7.25E-11*EXP(20./temp))
  RCONST(10) = (3.E-12*EXP(-1500./temp))
  RCONST(11) = (5.1E-12*EXP(210./temp))
  RCONST(12) = (1.2E-13*EXP(-2450./temp))
  RCONST(13) = (k_NO3_NO2)
  RCONST(14) = (k_3rd(temp,cair,1.8E-30,3.0,2.8E-11,0.,0.6))
  RCONST(15) = (k_HNO3_OH)
  RCONST(16) = (jx(ip_O2))
  RCONST(17) = (jx(ip_O1D))
  RCONST(18) = (jx(ip_O3P))
  RCONST(19) = (jx(ip_NO2))
  RCONST(20) = (jx(ip_NO2O))
  RCONST(21) = (jx(ip_N2O5))
  RCONST(22) = (jx(ip_HNO3))
      
END SUBROUTINE Update_RCONST

! End of Update_RCONST function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Update_PHOTO - function to update photolytical rate constants
!   Arguments :
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Update_PHOTO ( )


   USE messy_mecca_kpp_global

  RCONST(16) = (jx(ip_O2))
  RCONST(17) = (jx(ip_O1D))
  RCONST(18) = (jx(ip_O3P))
  RCONST(19) = (jx(ip_NO2))
  RCONST(20) = (jx(ip_NO2O))
  RCONST(21) = (jx(ip_N2O5))
  RCONST(22) = (jx(ip_HNO3))
      
END SUBROUTINE Update_PHOTO

! End of Update_PHOTO function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



 END MODULE messy_mecca_kpp_rates 

