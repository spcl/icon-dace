!<
!! mo_art_watercont_seas
!! Calculation of water content of seasalt aerosol
!!
!! Author: Heike Vogel, Sven Werchner, KIT
!! Initial Release: 2023-06-21
!!
!! Modifications:
!! YYYY-MM-DD: <name>, <institution>
!! - ...
!!
!! Copyright and License
!! This code is subject to the KIT-Software-License-Agreement in
!! its most recent form.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!! This software is provided for non-commercial use only.
!!
!! See the LICENSE conditions:
!! http://icon-art.imk-tro.kit.edu/
!! Institute of Meteorology and Climate Research, KIT, Karlsruhe
!>
MODULE mo_art_seas_watercont
! ICON
  USE mo_kind,                          ONLY: wp
! ART

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_watercont'

  PUBLIC :: art_watercont

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_watercont(relhum,temp,rho,vseas,vso4seas,vh2oseas,istart,iend,kstart,kend)
!<
! SUBROUTINE art_watercont
!********************************************************************************************
!Description of the method:
!The aerosol liquid content is calculated for an aerosol containing NaCl and
!H2SO4 with the use of the ZSR mixing rule [Zdanovskii,1948;Stokes and
!Robinson,1966].

!The ZSR mixing rule gives the liquid water content wh2o according to:
! wh2o=(n1/m1o(aw)+n2/m2o(aw)+n3/m3o(aw)+..)*1e9
!where
!--> wh2o:   total aerosol liquid water content     [ug/ m**3]
!--> ni:     concentration of electrolyte i         [mol/m**3]
!--> mio(aw):binary electrolyte molality, which
!            implies the amount of moles of solute
!            per kg solvent(water) if the electrolyte
!            were alone in solution at the
!            ambient RH(=aw).                        [mol/kg]
!--> aw:     water activity, which is the fractional
!            relative humidity RH if the Kelvin effect
!            is neglected.
!--> 1e9:    factor for transformation from kg to ug.

!The binary molality mi for each electrolyte i is calculated as function
! of aw with a polynomial of the form:
! mi0(aw)=Y0+Y1(aw)+Y2(aw**2)+Y3(aw**3)+Y4(aw**4)+..   (see Eq. 17.66 in Jacobson 2005)
!for which the coefficients Yi are for example found in
!Jacobson, 2005: "Fundamentals of atmospheric modeling", Appendix Table B.10, Page 749.
!(Can also be found in Tang et al., 1997.)

!The amount of moles of an electrolyte is most often not known directly if we have a mixture
!(i.e.,NaCl+H2SO4, and not only NaCl). Rather is the amount of each ion known and from
!this knowledge the amount of moles of each electrolyte is become by applying the Eq. 9 and 10
!in Zaveri et al., 2005 in JGR,110. doi:10.1029/1004JD004681 or equation 10 in Topping et al., 2009
!in JGR,114. doi:10.1029/2008JD10099. In this version the method applied by Topping and others is applied.
!No sulfate ratio is considered and NaCl and Na2SO4 are able to form.

!(When using the equation of Zaveri:
!According to e.g. Zaveri et al., 2005 two electrolyte formation domains are possible, depending
!on the sulfate ratio Xt=nNa/nSO4. In this code it is assumed that nNa >>nSO4. Thus, we are in
!the sulfate poor domain and sulfate is assumed to be completely neutralized by Na+. (Xt.ge.2).
!The sulfate rich domain is here neglected. In the sulfate poor domain the possible electrolytes
!Based on code of Kristina Lundgren COSMO-ART, 2009

! Part of Module: mo_art_watercont_seas
! Author: Heike Vogel, Sven Werchner, KIT
! Initial Release: 2023-06-25
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>


IMPLICIT NONE
!***main array of variables
  REAL(wp), INTENT(in)           :: &
    &  temp(:,:),                   &
    &  rho(:,:),                    &
    &  relhum(:,:),                 &
    &  vseas(:,:),                  & 
    &  vso4seas(:,:)                 

  INTEGER, INTENT(in)            :: &
    &  istart,iend,                 & !< Start and end index of nproma loop
    &  kstart,kend                    !< Start and end index of vertical loop

  REAL(wp), INTENT(inout)        :: &
    &  vh2oseas(:,:)                  !< water content of seasalt aerosol

  REAL(wp)                       :: &
    &  nNaCl,                       &  !electrolyte concentrations in each sea salt mode[mol/m**3]
    &  nNa2SO4,                     &  !-""-
    &  moNaCl,moNa2SO4,             &  !binary electrolyte molality [mol/kg]
    &  n0NaCl,                      &  !initial total mole concentrations in each mode (before ion pairing)[mol/m**3]
    &  n0H2SO4,                     &  !-"-
    &  nNa,                         &  !moles of each ion in each mode [mol/m**3]
    &  nCl,                         &  !-"-
    &  nSO4,                        &  !-"-
    &  nH,                          &  !-"-
    &  xt,                          &  !sulfate ratio [-]
    &  aw,                          &  !water activity [-]
    &  eNa,                         &  !equivalent ion fractions in each mode for ion pairing
    &  eH,                          &  !-"-
    &  eCl,                         &  !-"-
    &  eSO4,                        &  !-"-
    &  wh2o_1,                      &  !aerosol liquid water content with respect to NaCl for each mode [ug/m**3]
    &  wh2o_2,                      &  !aerosol liquid water content with respect to Na2SO4 for each mode [ug/m**3]
    &  growth                          !growth factor for each sea salt mode with respect to water uptake

  REAL(wp)                       :: &
    &  MwNa,MwCl,MwSO4,MwH,MwNaCl,  &  !Molecular weights         [g/mol]
    &  y00,y01,y02,y03,y04,         &  !parameters for NaCl molality
    &  y10,y11,y12,y13,y14                 !parameters for Na2SO4 molality

  INTEGER                        :: &
    &  jc, jk                  !< Loop counter


  PARAMETER (MwNa=22.990_wp,MwCl=35.453_wp,MwSO4=96.058_wp,MwH=1.008_wp,MwNaCl=58.443_wp)
  PARAMETER (y00=5.875248e1_wp, y01=-1.8781997e2_wp, y02=2.7211377e2_wp, y03=-1.8458287e2_wp, y04=4.153689e1_wp )
  PARAMETER (y10=5.5983158e2_wp,y11=-2.56942664e3_wp,y12=4.47450201e3_wp,y13=-3.45021842e3_wp,y14=9.8527913e2_wp)

!**************************
!Begin subroutine watercont
!**************************
!Set initial values
  nNaCl     =    1.0e-19_wp
  nNa2SO4   =    1.0e-19_wp
  moNaCl    =    1.0e-19_wp
  moNa2SO4  =    1.0e-19_wp
  n0NaCl    =    1.0e-19_wp
  n0H2SO4   =    1.0e-19_wp
  nNa       =    1.0e-19_wp
  nCl       =    1.0e-19_wp
  nSO4      =    1.0e-19_wp
  xt        =    1.0e-19_wp
  nH        =    1.0e-19_wp
  aw        =    1.0e-19_wp
  eNa       =    1.0e-19_wp
  eH        =    1.0e-19_wp
  eCl       =    1.0e-19_wp
  eSO4      =    1.0e-19_wp
  wh2o_1    =    1.0e-19_wp
  wh2o_2    =    1.0e-19_wp
  growth    =    1.0e-19_wp


  DO jk = kstart, kend
!NEC$ ivdep
    DO jc = istart, iend

!fetch the relative humidity.
      aw     =    relhum(jc,jk)*0.01_wp
!Calculate ion concentrations in moles/m3 from mass in ug/m3 for each mode
      n0NaCl  =  MAX(1.0e-19_wp,(vseas(jc, jk) * rho(jc,jk) / MwNaCl*1.0e-6_wp))
      nNa     =  n0NaCl
      nCl     =  nNa
      n0H2SO4 =  MAX(1.0e-19_wp,(vso4seas(jc,jk) * rho(jc,jk) / (2._wp*MwH+MwSO4)*1.e-06_wp))
      nSO4    =  n0H2SO4
      nH      =  2._wp*nSO4

!Determining the electrolyte concentration!
!Check that the sulfate ratio is within the interval xt.ge.2, which is assumed by [nacl]>>[h2so4].
!Calculate the modes separately and skip if sulfate ratio is too small.
      xt      =  MAX(1.0e-19_wp,n0NaCl/nSO4)
      IF(xt >= 2._wp) THEN
!equivalent ion fractions for ion-pairing
        eNa     =  MAX(1.0e-19_wp,(nNa   /(nNa +nH)))
        eH      =  MAX(1.0e-19_wp,(nH    /(nNa +nH)))
        eCl     =  MAX(1.0e-19_wp,(nCl   /(nCl +2._wp*nSO4)))
        eSO4    =  MAX(1.0e-19_wp,(2._wp*nSO4/(nCl +2._wp*nSO4)))
!calculate the electrolyte concentration in moles
        nNaCl   =  MAX(1.0e-19_wp,((eNa *nCl *MwCl   +eCl  *nNa *MwNa) / (MwNaCl)))
        nNa2SO4 =  MAX(1.0e-19_wp,((eNa *nSO4 *MwSO4 +eSO4 *nNa *MwNa) / (2._wp*MwNa+MwSO4)))
      ELSE !xt<2
        nNaCl   =  1.0e-19_wp
        nNa2SO4 =  1.0e-19_wp
      ENDIF !check the fine mode sulfate ratio

!Determine the binary molality with respect to the current humidity and the corresponding water content. Factor 1.d9 for kg-->ug
!IF(aw.LT.0.99_wp) THEN
      aw    =  MIN(0.99_wp,aw)
      IF(aw >= 0.47_wp) THEN
        moNaCl  =  y00+y01*aw+y02*aw**2+y03*aw**3+y04*aw**4
        wh2o_1  =  MAX(1.0e-19_wp,(nNaCl / moNaCl *1.0e+09_wp))
      ELSE !if rh<47% -->dry NaCl, no aerosol water content with respect to NaCl
        wh2o_1  =  1.0e-19_wp
      ENDIF

      IF(aw >= 0.58_wp) THEN
        moNa2SO4  =  y10+y11*aw+y12*aw**2+y13*aw**3+y14*aw**4
        wh2o_2    =  MAX(1.0e-19_wp,(nNa2SO4 / moNa2SO4 *1.0e+09_wp))
      ELSE !if rh<58% -->dry Na2SO4, no aerosol liquid water content with respect to Na2SO4
        wh2o_2    =  1.0e-19_wp
      ENDIF

!calculate the total aerosol water content for each mode as the sum of water with respect to the two electrolyte components [ug/m**3]
      vh2oseas(jc,jk) = MAX(1.0e-19_wp,(wh2o_1 + wh2o_2))
      vh2oseas(jc,jk) = vh2oseas(jc,jk) / (rho(jc,jk))
    ENDDO
  ENDDO


END SUBROUTINE art_watercont

!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_seas_watercont

