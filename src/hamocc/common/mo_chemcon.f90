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

! @brief T, S dependencies of chemical constants used in the inorganic carbon cycle.
!
! Computes chemical constants in the surface layer (aksurf)
! and in the water column (ak13, ak23, akb3, aksp, ak1p3, ak2p3, ak3p3
! aks3,akf3,aksi3,akw3)
! Parametrizations follow Dickson 2007, 2010 (Guides to best practise for CO2
! measurements, OA research)
! Constant parameter values can be found in mo_bgc_constants.f90

MODULE mo_chemcon


  USE mo_kind, ONLY        : wp
  USE mo_control_bgc, ONLY : bgc_nproma, bgc_zlevs
  USE mo_hamocc_nml, ONLY  : l_limit_sal

  USE mo_bgc_constants
  USE mo_bgc_memory_types, ONLY  : t_bgc_memory
  USE mo_fortran_tools, ONLY     : set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: chemcon
           

CONTAINS


SUBROUTINE CHEMCON (local_bgc_mem, start_idx, end_idx, klevs, psao, ptho,  &
     &               pddpo, ptiestu, kldtday, lacc)

  
  IMPLICIT NONE

  !! Arguments
  TYPE(t_bgc_memory), POINTER :: local_bgc_mem

  INTEGER, INTENT(in)  :: start_idx            !< start index for j loop (ICON cells, MPIOM lat dir)  
  INTEGER, INTENT(in)  :: end_idx              !< end index  for j loop  (ICON cells, MPIOM lat dir) 
  INTEGER, INTENT(in)  :: klevs(bgc_nproma)    !<  vertical levels

  REAL(wp) :: psao(bgc_nproma,bgc_zlevs)  !< salinity [psu]
  REAL(wp) :: ptho(bgc_nproma,bgc_zlevs)  !< potential temperature [deg C]
  REAL(wp) :: pddpo(bgc_nproma,bgc_zlevs) !< size of scalar grid cell (3rd REAL) [m]
  REAL(wp) :: ptiestu(bgc_nproma,bgc_zlevs)  !< depth of scalar grid cell [m]
  LOGICAL, INTENT(IN), OPTIONAL :: lacc

  !! Local variables

  INTEGER ::  jc, k, kldtday, kpke, js
  REAL(wp) :: t, q, s, log_t, log_q, sqrt_s, ti, qi, q2
  REAL(wp) :: cek0, ckb, ck1, ck2, ckw, oxy, ani
  REAL(wp) :: ak1, ak2, akb, akw, ak0, aksp0, log10ksp
  REAL(wp) :: p, cp, tc,tt,ts,ts2,ts3,ts4,ts5
  REAL(wp) :: pis, pis2, rs,s2,deltav,deltak,lnkpk0(11)
  REAL(wp) :: aksi, cksi, aks, cks, akf, ckf, free2sws,total2sws
  REAL(wp) :: ck1p,ck2p,ck3p, ak1p,ak2p,ak3p, total2free
  REAL(wp) :: sti, fti, total2free_0p, free2SWS_0p, total2SWS_0p,SWS2total
  LOGICAL :: lzacc

  CALL set_acc_host_or_device(lzacc, lacc)

  !     -----------------------------------------------------------------
  !*            SET MEAN TOTAL [CA++] IN SEAWATER (MOLES/KG)
  !             (SEE BROECKER A. PENG, 1982, P. 26)
  !             ([CA++](MOLES/KG)=1.026E-2*(S/35.) AFTER
  !             CULKIN(1965), CF. BROECKER ET AL. 1982)
  !             ------------- --- -------- -- --- -----
  ! Done in mo_ini_bgc, this routine assigns local-memory only
!   
!   calcon = 1.03e-2_wp
!   rrrcl = salchl * 1.025_wp * bor1 * bor2

  !
  !     -----------------------------------------------------------------
  !*        21. CHEMICAL CONSTANTS - SEA SURFACE
  !             --------------------------------
  !

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG VECTOR
  DO jc = start_idx, end_idx

      IF (pddpo(jc, 1) > EPSILON(0.5_wp)) THEN           ! wet cell
              !
              !*        21.1 SET ABSOLUTE TEMPERATURE
              !              ------------------------
              T = ptho(jc,1) + tmelt                  ! degC to K
              log_t = LOG(t)
              ti = 1.0_wp / t
              Q = T * PERC                             ! perc=0.01
              log_q = LOG(q)
              qi =  1.0_wp / q
              q2 = q**2
              IF (l_limit_sal) THEN
                s = MAX(25._wp, psao(jc, 1))           ! minimum salinity 25
              ELSE
                s = psao(jc, 1)
              ENDIF
              sqrt_s = SQRT(s)
              s2= s*s
              !
              !
              !*        21.2 CHLORINITY (WOOSTER ET AL., 1969)
              !              ---------------------------------

              !
              !*        21.3 SOLUBILITY OF GASES
              !              -------------------------------------------------
              !
              !      CO2 CARBON DIOXIDE, SEA SURFACE
              !      --------------------------------------
              !         LN(K0) OF SOLUBILITY OF CO2 (EQ. 12, WEISS, 1974)
              !

              cek0 = c00*ti + c01 + c02 * LOG_q + s * (c03 + c04*t + c05*q2)
              AK0 = EXP(CEK0)*SMICR
              local_bgc_mem%solco2(jc) = AK0
              

              pis = 19.924_wp * s / ( 1000._wp - 1.005_wp * s )
              pis2 = pis * pis

              !*         DISSOCIATION OF SILICIC ACID
              !
              !          aksi = [H][SiO(OH)3]/[Si(OH)4] 
              !          Millero p.671 (1995) using data from Yao and Millero
              !          (1995)

              CKSI =  cksi1/t + cksi2 + cksi3 * LOG(t) & 
                 &     + ( cksi4/t + cksi5 ) * sqrt(pis)  &
                 &     + ( cksi6/t + cksi7) * pis     &
                 &     + ( cksi8/t + cksi9) * pis2    &
                 &     + LOG(1.0_wp + cksi10*s)
              AKSI=EXP(CKSI)


              !
              !*         DISSOCIATION OF HYDROGEN SULFATE 
              !          aks = [H][SO4]/[HSO4]
              !          Dickson (1990, J. chem. Thermodynamics 22, 113)


              CKS =  cks1/t + cks2 + cks3 * LOG(t) &
                 &    + ( cks4/t + cks5 + cks6 * LOG(t) ) * sqrt(pis) &
                 &    + ( cks7/t + cks8 + cks9 * LOG(t) ) * pis &
                 &    + cks10/t * pis**1.5_wp & 
                 &    + cks11/t * pis2 &
                 &    + LOG(1.0_wp + cks12 * s ) 
              AKS=EXP(CKS)    

           
           

              !
              !*         DISSOCIATION OF HYDROGEN FLUORIDE 
              !
              !          akf = [H][F]/[HF]
              !          Dickson and Riley (1979) 

              CKF = ckf1*ti + ckf2  + ckf3*sqrt_s
              AKF=EXP(CKF)


                            !

              !*         DISSOCIATION OF PHOSPHORIC ACID 
              !
              !          DOE(1994) eq 7.2.20 with footnote using data from
              !          Millero (1974)
              !          AK1P, AK2P, AK3P
              !
              !          ak1p = [H][H2PO4]/[H3PO4] 


              CK1P =  ck1p1/t + ck1p2 + ck1p3 * LOG(t) &
                &  + ( ck1p4/t + ck1p5 ) *  sqrt(s) &
                &  + ( ck1p6/t + ck1p7 ) * s 
              AK1P=EXP(CK1P)

              !          ak2p = [H][HPO4]/[H2PO4] 
              CK2P = ck2p1/t + ck2p2 + ck2p3 * LOG(t)  & 
                  &   + ( ck2p4/t + ck2p5 ) *  sqrt(s) &
                  &   + ( ck2p6/t + ck2p7 ) *s 
              AK2P=EXP(CK2P)

              !          ak3p = [H][PO4]/[HPO4]

              CK3P = ck3p1/t + ck3p2                  &
                &     + ( ck3p3/t + ck3p4 ) * sqrt(s) &
                &     + ( ck3p5/t + ck3p6 ) * s 
              AK3P=EXP(CK3P)



              !      N2O LAUGHING GAS, SEA SURFACE  
              !      --------------------------------------
              rs = a1 + a2 * (100._wp *ti) + a3 * LOG(t*0.01_wp) &
             &    + a4*(t*0.01_wp)**2 + s * ( b1 +b2*(t*0.01_wp) + b3*(t*0.01_wp)**2)

              local_bgc_mem%satn2o(jc) = atn2o*EXP(rs)
              !
              !      O2 OXYGEN, SEA SURFACE
              ! Garcia & Gordon, 1992, EQ. 8, p 1310
              ! w/o A3*ts**2 --> see mocsy -> gasx.f90, OMIP paper 
              ! t in degC, s in permill

              !      -----------------------------------------------

               ! t in degC, s in permill

                tt  = 298.15_wp - ptho(jc,1)
                ts  = LOG(tt/t)
                ts2 = ts*ts
                ts3 = ts*ts2
                ts4 = ts*ts3
                ts5 = ts*ts4

              

              ! O2sat ml/L
                oxy  = oxya0 + oxya1*ts + oxya2*ts2 + oxya3*ts3 + oxya4*ts4 +oxya5*ts5  &
               + s*(oxyb0 + oxyb1*ts + oxyb2*ts2 + oxyb3*ts3)         &
               + oxyc0*(s*s)


              ! satoxy is a 3d-field
              local_bgc_mem%satoxy(jc,1) = EXP(OXY)*OXYCO


              !      N2 DINITROGEN, SEA SURFACE
              !      -----------------------------------------------
              !        WEISS, R. F. (1970), DEEP-SEA RESEARCH, VOL. 17, 721-735.

              ani = an0 + an1 * qi + an2 * log_q + an3 * q &
                  + s * (an4 + an5 * q + an6 * q2)
              local_bgc_mem%satn2(jc)  = EXP(ANI)*OXYCO

              !
              !
              !*        21.4 DISSOZIATION CONSTANTS
              !              -----------------------------------------
              !
              !*      PK1 and PK2  OF CARBONIC ACID
              !       -------------------------------

              ck1 =  c10*ti + c11 + c12*log_t + c13*s + c14* s2
              ck2 =  c20*ti + c21 + c22*log_t + c23*s + c24* s2



              !*      PKB OF BORIC ACID
              !       -------------------------------

              CKB= (cb0 + cb1*sqrt_s + cb2*s + cb3*(s*sqrt_s) + cb4*(s**2)) * ti &
                   + cb5 + cb6*sqrt_s + cb7*s                                &
                   -(cb8+cb9*sqrt_s + cb10*s) * log_t                       &
                   + cb11*sqrt_s*t
              !
              !*      PKW OF WATER (H2O)
              !       ------------------------------------

              CKW = cw0  + cw1*ti + cw2*log_t                                &
                    +SQRT_S*(cw3*ti + cw4 +cw5*log_t) + cw6*s

              !*        21.5 K1, K2 OF CARB. ACID, KB OF BORIC ACID, KW of WATER
              !              ----------------------------------------------------------------
              AK1 = 10._wp**CK1
              AK2 = 10._wp**CK2
              AKB = EXP(CKB)
              AKW = EXP(CKW)

             ! Conversion factors to total scale
             ! sulfate Morris & Riley (1966)
              sti   = 0.14_wp *  s*1.025_wp/1.80655_wp  / 96.062_wp

             ! fluoride Riley 1965
              fti    = 0.000067_wp * s*1.025_wp/1.80655_wp / 18.9984_wp
             ! conversion factor total to free scale 0p
              total2free_0p = 1._wp/(1._wp + sti/AKS)

             ! conversion factor free to sea water scale 0p
              free2SWS_0p  = 1._wp + sti/AKS + fti/(AKF*total2free_0p)  ! using Kf on free scale

             ! conversion factor total to sea water 0p
              total2SWS_0p = total2free_0p * free2SWS_0p             ! KSWS =Ktotal*total2SWS

              SWS2total = 1._wp/total2SWS_0p

              


              !*       21.11 STORE CHEMICAL CONSTANTS AT SEA SURFACE
              !              pressure independent
               local_bgc_mem%aksurf(jc,1) = ak1
               local_bgc_mem%aksurf(jc,2) = ak2
               local_bgc_mem%aksurf(jc,3) = akb
               local_bgc_mem%aksurf(jc,4) = AKW * SWS2total
               local_bgc_mem%aksurf(jc,5) = aksi * SWS2total
               local_bgc_mem%aksurf(jc,6) = akf ! already at total scale
               local_bgc_mem%aksurf(jc,7) = aks !  here free=total scale
               local_bgc_mem%aksurf(jc,8) = ak1p * SWS2total
               local_bgc_mem%aksurf(jc,9) = ak2p * SWS2total
               local_bgc_mem%aksurf(jc,10) = ak3p * SWS2total

      ENDIF

  END DO
  !$ACC END PARALLEL

  !
  !     -----------------------------------------------------------------
  !*        22. CHEMICAL CONSTANTS - DEEP OCEAN
  !     ----------------------------------------------------------------

  IF ( kldtday == 1 ) THEN
 
     !
     !*     22.1 APPROX. SEAWATER PRESSURE AT U-POINT DEPTH (BAR)
     !  ----------------------------------------------------------------

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG VECTOR PRIVATE(lnkpk0)
  DO jc = start_idx, end_idx
     kpke=klevs(jc)
     !$ACC LOOP SEQ
      DO k = 1, kpke
        p = 1.025e-1_wp * ptiestu(jc,k)   ! pressure
           IF(kpke.ne.0)THEN
           IF (pddpo(jc, k) > EPSILON(0.5_wp)) THEN           ! wet cell
           !
           !*    22.1 SET ABSOLUTE TEMPERATURE
           ! ----------------------------------------------------------------

           t = ptho(jc,k) + tmelt
           ti = 1.0_wp / t
           log_t = LOG(t)
           q = t * perc
           log_q = LOG(q)
           qi = 1.0_wp / q
           !
           !*    22.2 SET LIMITS FOR SALINITY
           ! ----------------------------------------------------------------
           !   (THIS IS DONE TO AVOID COMPUTATIONAL CRASH AT DRY
           !   POINTS DURING CALCULATION OF CHEMICAL CONSTANTS)

           IF (l_limit_sal) THEN
             s = MAX(25._wp, psao(jc, k))           ! minimum salinity 25
           ELSE
             s = psao(jc, k)
           ENDIF

           sqrt_s = SQRT(s)
           s2 = s*s


           !
           !*    22.4 SOLUBILITY OF O2 
           ! 2 saturation concentration at 1atm total pressure
           ! Garcia & Gordon, 1992, EQ. 8, p 1310
           ! w/o A3*ts**2 --> see mocsy -> gasx.f90, OMIP paper 
           ! t in degC, s in permill
           tt  = 298.15_wp - ptho(jc,k)
           ts  = LOG(tt/t)
           ts2 = ts*ts
           ts3 = ts*ts2
           ts4 = ts*ts3
           ts5 = ts*ts4

           ! O2sat ml/L
           oxy  = oxya0 + oxya1*ts + oxya2*ts2 + oxya3*ts3 + oxya4*ts4 +oxya5*ts5  &
               + s*(oxyb0 + oxyb1*ts + oxyb2*ts2 + oxyb3*ts3)         &
               + oxyc0*(s*s)


           if (k > 1) local_bgc_mem%satoxy(jc,k) = EXP(oxy) * oxyco

           !
           !*    22.5 PK1, PK2 OF CARBONIC ACID, PKB OF BORIC ACID, PKW OF WATER
           ! ----------------------------------------------------------------
           ck1 =  c10*ti + c11 + c12*log_t + c13*s + c14* s2
           ck2 =  c20*ti + c21 + c22*log_t + c23*s + c24* s2


           CKB= (cb0 + cb1*sqrt_s + cb2*s + cb3*(s*sqrt_s) + cb4*(s**2)) * ti &
                   + cb5 + cb6*sqrt_s + cb7*s                                &
                   -(cb8+cb9*sqrt_s + cb10*s) * log_t                       &
                   + cb11*sqrt_s*t

           CKW = cw0  + cw1*ti + cw2*log_t                                &
                    +sqrt_s*(cw3*ti + cw4 +cw5*log_t) + cw6*s

           !
           !*    22.6 K1, K2 OF CARBONIC ACID, KB OF BORIC ACID, KW (H2O)
           ! ----------------------------------------------------------------
           AK1 = 10._wp**CK1
           AK2 = 10._wp**CK2
           AKB = EXP(CKB)
           AKW = EXP(CKW)


           !*         DISSOCIATION OF SILICIC ACID
           !
           !          aksi = [H][SiO(OH)3]/[Si(OH)4] 
           !          Millero p.671 (1995) using data from Yao and Millero
           !          (1995)

           pis = 19.924_wp * s / ( 1000._wp - 1.005_wp * s )
           pis2 = pis * pis

           CKSI =  cksi1/t + cksi2 + cksi3 * LOG(t) & 
              &     + ( cksi4/t + cksi5 ) * sqrt(pis)  &
              &     + ( cksi6/t + cksi7) * pis     &
              &     + ( cksi8/t + cksi9) * pis2    & 
              &     + LOG(1.0_wp + cksi10*s)

           AKSI=EXP(CKSI)

           !
           !*         DISSOCIATION OF HYDROGEN SULFATE 
           !          aks = [H][SO4]/[HSO4]
           !          Dickson (1990, J. chem. Thermodynamics 22, 113)

           

           CKS =  cks1/t + cks2 + cks3 * LOG(t) &
              &    + ( cks4/t + cks5 + cks6 * LOG(t) ) * sqrt(pis) &
              &    + ( cks7/t + cks8 + cks9 * LOG(t) ) * pis &
              &    + cks10/t * pis**1.5_wp & 
              &    + cks11/t * pis2 &
              &    + LOG(1.0_wp + cks12 * s ) 

           AKS=EXP(CKS)    


          !
          !*         DISSOCIATION OF HYDROGEN FLUORIDE 
          !
          !          akf = [H][F]/[HF]
          !          Dickson 2007/2010 (lomip)
          !          Dickson and Riley (1979) 


          CKF = ckf1*ti + ckf2  + ckf3*sqrt_s
          AKF=EXP(CKF)

           !
           !*         DISSOCIATION OF PHOSPHORIC ACID 
           !
           !          DOE(1994) eq 7.2.20 with footnote using data from Millero
           !          (1974)
           !          AK1P, AK2P, AK3P
           !
           !          ak1p = [H][H2PO4]/[H3PO4] 

               

           CK1P =  ck1p1/t + ck1p2 + ck1p3 * LOG(t) &
             &  + ( ck1p4/t + ck1p5 ) *  sqrt(s) &
             &  + ( ck1p6/t + ck1p7 ) * s 


           AK1P=EXP(CK1P)

           !          ak2p = [H][HPO4]/[H2PO4] 
           CK2P = ck2p1/t + ck2p2 + ck2p3 * LOG(t)  & 
               &   + ( ck2p4/t + ck2p5 ) *  sqrt(s) &
               &   + ( ck2p6/t + ck2p7 ) *s 

               

           AK2P=EXP(CK2P)

           !          ak3p = [H][PO4]/[HPO4]
           CK3P = ck3p1/t + ck3p2                  &
             &     + ( ck3p3/t + ck3p4 ) * sqrt(s) &
             &     + ( ck3p5/t + ck3p6 ) * s 
           AK3P=EXP(CK3P)

           TC = ptho(jc,k)


           !
           !       FORMULA FOR CP
           !       ----------------------------------------------

           CP = P / (RGAS*T)

           !
           !*    22.7 APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE IN SEAWATER
           ! ----------------------------------------------------------------
           !     Mucci (1983)

           log10ksp = akcc1 + akcc2*t + akcc3*ti + akcc4*log_t*invlog10       &
                 + (akcc5 +akcc6*t + akcc7*ti) * sqrt_s                 &
                 + akcc8*s + akcc9*(s*sqrt_s)
           aksp0 = 10._wp**(log10ksp)


          ! Pressure effect for diss. const.s of  hydrogen fluoride,
          ! hydrogen sulf ate, phosphoric acid
          ! based on Millero 1995
          !   1) aks 2) akf 3) ak1p 4) ak2p 
          !   5) ak3p 6) aksi  7) ak1 8) ak2 
          !   9) akb 10) akw 11) aksp 

           DO js = 1,11

            deltav      = pa0(js) + pa1(js) * tc + pa2(js) * tc * tc
            deltak      = pb0(js) + pb1(js) * tc + pb2(js) * tc * tc
            lnkpk0(js) = - ( deltav * cp + 0.5_wp * deltak * cp * p )

           ENDDO

           ! sulfate Morris & Riley (1966)
           sti   = 0.14_wp *  s*1.025_wp/1.80655_wp  / 96.062_wp

           ! fluoride Riley 1965
           fti    = 0.000067_wp * s*1.025_wp/1.80655_wp / 18.9984_wp

          ! Conversion at pressure zero
          ! conversion factor total to free scale 0p
           total2free_0p = 1._wp/(1._wp + sti/AKS)

          ! conversion factor free to sea water scale 0p
          free2SWS_0p  = 1._wp + sti/AKS + fti/(AKF*total2free_0p)  ! using Kf on free scale

          ! conversion factor total to sea water 0p
          total2SWS_0p = total2free_0p * free2SWS_0p             ! KSWS =Ktotal*total2SWS

          ! Pressure correction on Ks (free scale)
           AKS  = AKS  * EXP(lnkpk0(1))

         ! conversion factor total to free scale
           total2free = 1._wp/(1._wp + sti/AKS)


          ! Pressure correction on Kf (free scale)
           AKF  = AKF  * EXP(lnkpk0(2))

          ! convert to total scale 
           AKF  = AKF / total2free

         ! Convert between seawater and total hydrogen (pH) scales
          free2SWS  = 1._wp + sti/aks + fti/(akf*total2free)  ! using Kf on free scale

          total2SWS = total2free * free2SWS                   ! KSWS =Ktotal*total2SWS

          SWS2total = 1._wp / total2SWS

         ! Convert K1, K2, Kb to seawater scale
           AK1 = AK1 * total2SWS_0p    
           AK2 = AK2 * total2SWS_0p    
           AKB = AKB * total2SWS_0p    

        ! Pressure correction 

           AK1  = AK1 * EXP(lnkpk0(7)) 
           AK2  = AK2 * EXP(lnkpk0(8)) 
           AKB = AKB * EXP(lnkpk0(9)) 
           AKW = AKW * EXP(lnkpk0(10)) 
           AKSP0 = ARACAL* AKSP0 * EXP(lnkpk0(11)) 
           AK1P = AK1P * EXP(lnkpk0(3))
           AK2P = AK2P * EXP(lnkpk0(4))
           AK3P = AK3P * EXP(lnkpk0(5))
           AKSI = AKSI * EXP(lnkpk0(6)) 

          ! Conversion to total scale and move to 3D var

           local_bgc_mem%aks3(jc,k)  = aks ! here free=total scale
           local_bgc_mem%akf3(jc,k)  = akf ! already at total scale  
           local_bgc_mem%ak1p3(jc,k) = ak1p * sws2total
           local_bgc_mem%ak2p3(jc,k) = ak2p * sws2total
           local_bgc_mem%ak3p3(jc,k) = ak3p * sws2total
           local_bgc_mem%aksi3(jc,k) = aksi * sws2total 
           local_bgc_mem%ak13(jc,k)  = ak1  * sws2total
           local_bgc_mem%ak23(jc,k)  = ak2  * sws2total
           local_bgc_mem%akb3(jc,k) = akb  * sws2total
           local_bgc_mem%akw3(jc,k) = akw  * sws2total
           local_bgc_mem%aksp(jc,k) = aksp0  ! independent


          END IF
          END IF
        END DO
  END DO
  !$ACC END PARALLEL



  ENDIF  ! kldtay == 1

END SUBROUTINE CHEMCON

END MODULE
