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

MODULE mo_gme_turbdiff

  USE mo_kind,               ONLY : wp
  USE mo_physical_constants, ONLY: grav, cpd, rcpd, vtmpc1, p0ref, rd_o_cpd, &
                                   tmelt, alvdcp, alv, als, rd, rdv, O_m_rdv
  USE mo_satad,              ONLY: sat_pres_water, sat_pres_ice, spec_humi, dqsatdT
  USE data_turbulence,       ONLY: Rkarman => akt, tkhmin, tkmmin
  USE mo_lnd_nwp_config,     ONLY: lseaice

  USE mo_exception,          ONLY: message
  USE mo_run_config,         ONLY: msg_level

  USE mo_fortran_tools,      ONLY: set_acc_host_or_device, assert_acc_device_only

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: partura, parturs, progimp_turb, nearsfc




  LOGICAL, PARAMETER :: using_cuda_graph = .FALSE.


  CONTAINS

SUBROUTINE partura( zh  , zf , u  , v  , t   ,  &
     &              qv  , qc , ph , pf ,        &
     &              ie  , ke , ke1,             &
     &              i_startidx, i_endidx,       &
     &              tkvm, tkvh)
!
 
!**** *partura*  calculates atmospheric turbulent exchange coefficients
 
!     Purpose:   Calculation of turbulent atmospheric exchange
!                coefficients for momentum and heat/moisture
!     METHOD
!     ------
!     -  second order closure on level 2 following MELLOR and YAMADA
!      
!==============================================================================

  IMPLICIT NONE
!
! array dimensions
  INTEGER, INTENT(IN) :: ie, ke, ke1, &        ! horizontal, vertical
                         i_startidx, i_endidx  ! start and end indices of loops in horizontal

  REAL(wp), INTENT(IN) :: zh (ie,ke1),&  ! height (m) of half levels
                          zf (ie,ke), &  ! height (m) of full levels
                          u  (ie,ke), &  ! zonal wind   (m/s)
                          v  (ie,ke), &  ! merid. wind  (m/s)
                          t  (ie,ke), &  ! layer temperature (K)
                          qv (ie,ke), &  ! water vapour (kg/kg)
                          qc (ie,ke), &  ! liquid water (kg/kg)
                          ph (ie,ke1),&  ! pressure at half levels (Pa)
                          pf (ie,ke)     ! pressure at full levels (Pa)

! turbulent diffusion coefficients (m^2/s)
  REAL(KIND=wp), INTENT(INOUT) :: tkvm(ie,2:ke) ! for momentum
  REAL(KIND=wp), INTENT(INOUT) :: tkvh(ie,2:ke) ! for heat
 
! Local arrays and variables

      REAL(KIND=wp) :: ztpm   (ie)
      REAL(KIND=wp) :: ztph   (ie)
      REAL(KIND=wp) :: zdzq   (ie)
      REAL(KIND=wp) :: zriza  (ie)
      REAL(KIND=wp) :: ztetlu (ie) 
      REAL(KIND=wp) :: ztetlo (ie)
      REAL(KIND=wp) :: zbetatu(ie)
      REAL(KIND=wp) :: zbetato(ie)
      REAL(KIND=wp) :: zbetawu(ie)
      REAL(KIND=wp) :: zbetawo(ie)
      REAL(KIND=wp) :: zdpu   (ie)
      REAL(KIND=wp) :: zdpo   (ie)
      REAL(KIND=wp) :: zdza   (ie)
      REAL(KIND=wp) :: ztphww (ie)
      REAL(KIND=wp) :: ztet   (ie)
      REAL(KIND=wp) :: zclcu  (ie)
      REAL(KIND=wp) :: zclco  (ie)
      REAL(KIND=wp) :: zclvcu (ie)
      REAL(KIND=wp) :: zclvco (ie)
      REAL(KIND=wp) :: zclwcu (ie)
      REAL(KIND=wp) :: zclwco (ie)

      INTEGER :: j1,j3 ! loop indices over spatial dimensions
      INTEGER :: ms       ! 
!    
!     local physical parameters
!      
!     minimum diffusion coefficient (fraction of neutral value) for
!     very stable conditions (RI >= RIK) in the free atmosphere
      REAL(KIND=wp) :: ztmmin_a  ! for momentum
      REAL(KIND=wp) :: zthmin_a  ! for heat
!     minimum diffusion coefficient (absolute value) for
!     very stable conditions (RI >= RIK) in the free atmosphere
      REAL(KIND=wp) :: ztkmmin_a  ! for momentum
      REAL(KIND=wp) :: ztkhmin_a  ! for heat

!     local utility variables
      REAL(KIND=wp) :: zgqd, ztl
      REAL(KIND=wp) :: zsw, zsi, zqdw, zedtetv, zuc, zsigma, zgewo, zgewu
      REAL(KIND=wp) :: zdpn2, ztn,ztvo,ztvu, zcs, zclwcmn, zdqdt
      REAL(KIND=wp) :: za, zb, zmzb, zbetatn, zbetawn, zbetatog, zbetawog
      REAL(KIND=wp) :: ztphwg, zrf, zgam, zgs, zsm, zalh, ztkvmom, ztkvhom
!
!=======================================================================

!     constants
      REAL(KIND=wp), PARAMETER :: zrik   = 0.380_wp  ! critical Richardson number
      REAL(KIND=wp), PARAMETER :: zalf   = 0.900_wp  ! weighting coefficients for temporal
      REAL(KIND=wp), PARAMETER :: z1malf = 1._wp - zalf != 0.100_wp ! smoothing of exchange coefficients
      REAL(KIND=wp), PARAMETER :: zalhn  = 1.00_wp
      REAL(KIND=wp), PARAMETER :: zdze   = 5.E2_wp   ! asymptotic limit for height in mixing length relation
      REAL(KIND=wp), PARAMETER :: zdvbmin= 0.01_wp
      REAL(KIND=wp), PARAMETER :: zclwfak= 0.005_wp  ! = 0.5%

!                                              stable    unstable
      REAL(KIND=wp), PARAMETER :: zaa (2) = (/ 3.7_wp    , 4.025_wp  /)
      REAL(KIND=wp), PARAMETER :: zab1(2) = (/ 2.5648_wp , 3.337_wp  /)
      REAL(KIND=wp), PARAMETER :: zab2(2) = (/ 1.1388_wp , 0.688_wp  /)
      REAL(KIND=wp), PARAMETER :: zac1(2) = (/ 0.8333_wp , 1.285_wp  /)
      REAL(KIND=wp), PARAMETER :: zac2(2) = (/ 0.2805_wp , 0.2305_wp /)
      REAL(KIND=wp), PARAMETER :: zac3(2) = (/ 0.1122_wp ,-0.1023_wp /)
 
!     minimum diffusion coefficient (fraction of neutral value) for
!     very stable conditions (RI >= RIK) in the ABL
      REAL(KIND=wp), PARAMETER :: ztmmin = 0.010_wp   ! for momentum
      REAL(KIND=wp), PARAMETER :: zthmin = 0.007_wp   ! for heat
 
!     minimum diffusion coefficient (absolute value) for
!     very stable conditions (RI >= RIK) in the ABL (specified via namelist now)
      REAL(KIND=wp) :: ztkmmin  ! old default 1.0_wp  ! for momentum
      REAL(KIND=wp) :: ztkhmin  ! old default 0.4_wp  ! for heat

!     relative humidities to calculate partial cloud cover
      REAL(KIND=wp), PARAMETER :: Rh_cr1 =  0.8_wp
      REAL(KIND=wp)            :: Rh_cr2
      REAL(KIND=wp), PARAMETER :: Rh_cr3 =  1.0_wp
!
!=======================================================================
!
      IF ( msg_level >= 15) CALL message( 'mo_gme_turbdiff:', 'partura')

      ! take minimum diffusion coefficients from namelist
      ztkmmin   = tkmmin
      ztkhmin   = tkhmin

!     Presettings
      ztmmin_a  = ztmmin *0.1_wp
      zthmin_a  = zthmin *0.1_wp
      ztkmmin_a = ztkmmin*1.0_wp
      ztkhmin_a = ztkhmin*1.0_wp

      Rh_cr2 =  SQRT(3._wp)

!     Prepare calculation of modified Richardson number in the presence
!     of clouds

!     calculations for lowest layer
      DO j1 = i_startidx, i_endidx
!       Cloud cover
        IF ( qc(j1,ke) > 0.0_wp ) THEN
          zclcu (j1) = 1.0_wp      ! cloud cover
          ztet  (j1) = t(j1,ke) * EXP( rd_o_cpd*LOG(p0ref/pf(j1,ke)) )
          ztetlu(j1) = ztet(j1)*( 1._wp-alvdcp*qc(j1,ke)/t(j1,ke) )
          ztl        = t(j1,ke) - qc(j1,ke)*alvdcp
          zgqd       = spec_humi( sat_pres_water(ztl), pf(j1,ke))
          zdqdt      = dqsatdT ( ztl, zgqd)
          za         = 1._wp/(1._wp+zdqdt*alvdcp)
          zb         = za*t(j1,ke)*zdqdt/ztet(j1)
          zedtetv    = 1._wp/( ztet(j1)*(1._wp+vtmpc1*qv(j1,ke) - qc(j1,ke)))

!         Store beta-parameters for lower layer
          zbetatu(j1)= 1._wp/ztet(j1) - zb*( alvdcp/t(j1,ke) -      &
     &                       (1._wp+vtmpc1)*ztet(j1)*zedtetv )
          zbetawu(j1)= ztet(j1)*zedtetv*( vtmpc1 - (1._wp+vtmpc1)*za)  &
     &                 + za*alvdcp/t(j1,ke)
        ELSE
          zclcu(j1) = 0.0_wp      ! cloud cover
        ENDIF

!       Store liquid water and water vapour content of lowest layer
        zclwcu(j1) = qc(j1,ke)
        zclvcu(j1) = qv(j1,ke) 

      END DO     ! lowest layer calculations

      DO j3 = ke, 2, -1     ! Vertical loop
 
!       Calculation of non-convective partial cloud cover
        DO j1 = i_startidx, i_endidx
 
!         Critical relative humdity as function of p/ps     
          zsigma = pf(j1,j3-1)/ph(j1,ke1)
          zuc    = 0.95_wp - Rh_cr1*zsigma*(1._wp-zsigma)    &
           &                 *(1._wp+Rh_cr2*(zsigma-0.5_wp)) 
 
!         total water content (water vapour + liquid) of top layer
          zqdw   = qv(j1,j3-1) + qc(j1,j3-1)
!         saturation specific humdity over water
          zsw    = spec_humi ( sat_pres_water( t(j1,j3-1)), pf(j1,j3-1))

!         cloud liquid water limit value: 0.5% of saturation humidity
          IF ( t(j1,j3-1) >= tmelt ) THEN
            zclwcmn = zsw*zclwfak  
          ELSE
!           saturation specific humdity over ice
            zsi     = spec_humi ( sat_pres_ice  ( t(j1,j3-1)), pf(j1,j3-1))
            zclwcmn = zsi*zclwfak             !  for T<tmelt
          ENDIF
 
!         partial cloud cover as function of relative humidity   
          zcs    = MAX( 0.0_wp, MIN( 1._wp,(zqdw/zsw-zuc)/(Rh_cr3-zuc)))**2
 
!         store top layer partial cloud cover and liquid water content
          zclco(j1) = zcs
          IF(zclco(j1) > 0.0_wp) THEN
            zclwco(j1) = zclwcmn
          ELSE
            zclwco(j1) = 0.0_wp
          END IF

!         grid-scale cloud existence has first priority
          IF ( qc(j1,j3-1) > 0.0_wp ) THEN
            zclco(j1)= 1.0_wp         
            zclwco(j1)= qc(j1,j3-1) 
          ENDIF
 
!         total water content must not be affected by partial cloud cover
          zclvco(j1)= qv(j1,j3-1) + qc(j1,j3-1) - zclwco(j1)
 
        END DO
 
!     compute source functions and Richardson number in cloud-free
!     situation: Zu is below (unten), Zo is above (oben)
        DO j1 = i_startidx, i_endidx
          zdpu(j1) = ph(j1,j3+1) - ph(j1,j3)
          zdpo(j1) = ph(j1,j3)   - ph(j1,j3-1)
          zdpn2    = ph(j1,j3+1) - ph(j1,j3-1)
          ztvu     = t(j1,j3  ) * (1._wp+vtmpc1*qv(j1,j3  ))
          ztvo     = t(j1,j3-1) * (1._wp+vtmpc1*qv(j1,j3-1))
          ztn      = (zdpu(j1)*ztvo + zdpo(j1)*ztvu)/zdpn2
          zdza(j1) = zf(j1,j3) - zf(j1,j3-1)
          ztpm(j1) = (  (u(j1,j3) - u(j1,j3-1))**2    &
           &          + (v(j1,j3) - v(j1,j3-1))**2 )  &
           &         / zdza(j1)**2
          ztpm(j1) = MAX( ztpm(j1) , (zdvbmin/zdza(j1))**2 )
          ztph(j1) = grav*  ((ztvu - ztvo)/zdza(j1) + grav*rcpd)/ztn
          zmzb        = zh(j1,j3) - zh(j1,ke1)
          zdzq(j1) = (Rkarman*zmzb/(1._wp+zmzb/zdze))**2

!       utility variable for modified Richardson number 
          ztet(j1) = t(j1,j3-1) * EXP( rd_o_cpd*LOG(p0ref/pf(j1,j3-1)) )
        END DO

!     modified Richardson number in the presence of clouds
        DO j1 = i_startidx, i_endidx
 
!         at cloud base, only upper layer values are computed 
          IF( zclco(j1) > 0.0_wp .AND. zclcu(j1) <= 0.0_wp ) THEN
            ztetlo(j1) = ztet(j1)*( 1._wp-alvdcp*zclwco(j1)/t(j1,j3-1) )
            ztl        = t(j1,j3-1)-zclwco(j1)*alvdcp
            zgqd       = spec_humi ( sat_pres_water( ztl), pf(j1,j3-1))
            zdqdt      = dqsatdT (ztl, zgqd)
            za         = 1._wp/(1._wp+zdqdt*alvdcp)
            zb         = za*t(j1,j3-1)*zdqdt/ztet(j1)
            zedtetv    = 1._wp/(ztet(j1)*(1._wp + vtmpc1*zclvco(j1) - zclwco(j1)))
            zbetato(j1)= 1._wp/ztet(j1) - zb*( alvdcp/t(j1,j3-1)      &
             &                  -(1._wp+vtmpc1)*ztet(j1)*zedtetv )
            zbetawo(j1)= ztet(j1)*zedtetv*(vtmpc1-(1._wp+vtmpc1)*za)  &
             &           + za*alvdcp/t(j1,j3-1)
          ELSE 
            ztetlo(j1) = 0._wp
            zbetato(j1)= 0._wp
            zbetawo(j1)= 0._wp
          ENDIF

        END DO
!
!     inside of cloud
        DO j1 = i_startidx, i_endidx
          IF( zclco(j1) > 0.0_wp .AND. zclcu(j1) > 0.0_wp ) THEN
            ztetlo(j1) = ztet(j1)*(1._wp-alvdcp*zclwco(j1)/t(j1,j3-1))
            ztl        = t(j1,j3-1)-zclwco(j1)*alvdcp
            zgqd       = spec_humi ( sat_pres_water(ztl), pf(j1,j3-1))
            zdqdt      = dqsatdT (ztl,zgqd)
            za         = 1._wp/(1._wp+zdqdt*alvdcp)
            zb         = za*t(j1,j3-1)*zdqdt/ztet(j1)
            zedtetv    = 1._wp/(ztet(j1)*(1._wp+vtmpc1*zclvco(j1) - zclwco(j1)))
            zbetato(j1)= 1._wp/ztet(j1)-zb*(alvdcp/t(j1,j3-1)        &
     &                    -(1._wp+vtmpc1)*ztet(j1)*zedtetv)
            zbetawo(j1)= ztet(j1)*zedtetv*(vtmpc1-(1._wp+vtmpc1)*za)  &
     &                   + za*alvdcp/t(j1,j3-1)
            zbetatn    = zbetato(j1)*zbetatu(j1)* ( zdpu(j1)+zdpo(j1) )/  &
     &                    ( zdpu(j1)*zbetatu(j1)+zdpo(j1)*zbetato(j1) )
            zbetawn    = zbetawo(j1)*zbetawu(j1)* ( zdpu(j1)+zdpo(j1) )/  &
     &                    ( zdpu(j1)*zbetawu(j1)+zdpo(j1)*zbetawo(j1) )
            ztphww(j1) = grav* ( zbetatn * (ztetlo(j1)-ztetlu(j1))   &
     &                    + zbetawn * (zclwco(j1)+zclvco(j1)-        &
     &                                 zclwcu(j1)-zclvcu(j1)))/      &
     &                              (-1._wp*zdza(j1))
          ELSE
            ztphww(j1) = 0._wp
          ENDIF
 
!     if cloud cover in top layer is greater equal lower layer  
          IF( zclcu(j1) > 0._wp .AND. zclco(j1) >= zclcu(j1) ) THEN
            ztph(j1) = (1._wp-zclcu(j1))*ztph(j1) + zclcu(j1)*ztphww(j1)
          ENDIF
!
        END DO

!     cloud top or cloud cover is smaller in top layer
        DO j1 = i_startidx, i_endidx
          IF ( zclcu(j1) > 0._wp .AND. zclco(j1) < zclcu(j1) ) THEN
            zgewu      = 0.7_wp
            zgewo      = 1._wp-zgewu
            ztetlo(j1) = ztet(j1)
            ztl        = t(j1,j3-1)
            zgqd       = spec_humi ( sat_pres_water(ztl), pf(j1,j3-1))
            zdqdt      = dqsatdT (ztl,zgqd)
            za         = 1._wp/(1._wp+zdqdt*alvdcp)
            zb         = za*t(j1,j3-1)*zdqdt/ztet(j1)
            zedtetv    = 1._wp/(ztet(j1)*(1._wp+vtmpc1*zclvco(j1))) 
            zbetatog   = 1._wp/ztet(j1)
            zbetawog   = ztet(j1)*zedtetv*vtmpc1
            zbetatn    = (zgewo*zbetatog + zgewu*zbetatu(j1))/         &
     &                   (zgewu+zgewo)
            zbetawn    = (zgewo*zbetawog + zgewu*zbetawu(j1))/         &
     &                   (zgewu+zgewo)
            ztphwg     = grav* ( zbetatn * (ztetlo(j1)-ztetlu(j1)) +   &
     &                     zbetawn * (zclwco(j1)+zclvco(j1)-        &
     &                                zclwcu(j1)-zclvcu(j1)))/      &
     &                  (-1._wp*zdza(j1))

            ztph(j1) = (1._wp-zclcu(j1))*ztph(j1)            &
     &                   + (zclcu(j1)-zclco(j1) )*ztphwg     &
     &                   +  zclco(j1)*ztphww(j1)
          ENDIF
        END DO

        DO j1 = i_startidx, i_endidx
          zriza(j1) = ztph(j1)/ztpm(j1)
 
!       swap values needed for next layer:
          ztetlu(j1) = ztetlo(j1)
          zbetatu(j1)= zbetato(j1)
          zbetawu(j1)= zbetawo(j1)
          zclcu(j1)  = zclco(j1)
          zclwcu(j1) = zclwco(j1) 
          zclvcu(j1) = zclvco(j1)
        END DO
!
!=======================================================================
! 
!     compute vertical exchange coefficients
        DO j1 = i_startidx, i_endidx
      
!     very stable case   (Ri > Rik)
          IF(zriza(j1) >= zrik) THEN   ! very stable case

            IF( pf(j1,j3) >= 0.8_wp*ph(j1,ke1) ) THEN   ! Within the ABL
              tkvm(j1,j3) = ztmmin  *zdzq(j1)*SQRT( ztpm(j1) )
              tkvh(j1,j3) = zthmin  *zdzq(j1)*SQRT( ztpm(j1) )
            ELSE                     ! In the free atmosphere
              tkvm(j1,j3) = ztmmin_a*zdzq(j1)*SQRT( ztpm(j1) )
              tkvh(j1,j3) = zthmin_a*zdzq(j1)*SQRT( ztpm(j1) )
            ENDIF

          ELSE                            ! neutral or unstable case

            IF(zriza(j1) <= 0._wp) THEN
              ms = 2
            ELSE
              ms = 1
            END IF
 
!     flux Richardson number as function of gradient Richardson number
            zrf = zac1(ms)*( zriza(j1) + zac2(ms)         &
     &           -SQRT(zriza(j1)**2-zac3(ms)*zriza(j1)    &
     &                             +zac2(ms)**2) )
            zrf = MIN ( zrf, 0.99999_wp )
 
            zgam = zrf/(1.0_wp-zrf)    ! stability parameter gamma
 
!         SM(zgam)**1.5 and ALH(zgam)    
            zgs  = (1.0_wp-zab1(ms)*zgam)/(1.0_wp-zab2(ms)*zgam)
            zsm  = SQRT(((1._wp-zaa(ms)*zgam)/zgs)**3)
            zalh = zalhn*zgs
 
!         exchange coefficients for momentum and heat
!         a) for actual time step
            ztkvmom   = zdzq(j1)*zsm* SQRT( ztpm(j1) - zalh*ztph(j1))
            ztkvhom   = ztkvmom*zalh
!         b) weighted average using previous values
            tkvm(j1,j3) = zalf*ztkvmom + z1malf*tkvm(j1,j3) 
            tkvh(j1,j3) = zalf*ztkvhom + z1malf*tkvh(j1,j3) 
          ENDIF

!       set lower limit for exchange coefficients
          IF( pf(j1,j3) >= 0.8_wp*ph(j1,ke1) ) THEN   ! Within the ABL
            tkvm (j1,j3) = MAX ( ztkmmin,   tkvm(j1,j3) )
            tkvh (j1,j3) = MAX ( ztkhmin,   tkvh(j1,j3) )
          ELSE                     ! In the free atmosphere
            tkvm (j1,j3) = MAX ( ztkmmin_a, tkvm(j1,j3) )
            tkvh (j1,j3) = MAX ( ztkhmin_a, tkvh(j1,j3) )
          ENDIF

        END DO

      END DO     ! end of vertical loop

END SUBROUTINE partura

!=======================================================================

SUBROUTINE parturs( zsurf, z1  , u1   , v1     , t1   , qv1  ,    &
                    t_g  , qv_s, ps   , fr_land, h_ice,           &
                    ie   , i_startidx , i_endidx,                 &
                    tcm  , tch , gz0  ,                           &
                    shfl_s, lhfl_s, qhfl_s, umfl_s, vmfl_s, lacc )
 
!**** *parturs*  calculates turbulent transfer coefficients
!=======================================================================
!
!     Purpose:   Calculation of transfer coefficients for momentum
!                and heat (tcm/tch) and of the surface roughness
!                length over the sea ( gz0)
!     METHOD
!     ------
!     -  DYER-BUSINGER relations (modified by J.F.LOUIS)
!      
!------------------------------------------------------------------------------

  IMPLICIT NONE
!
!     Input data
!
  INTEGER,       INTENT(IN) :: ie,      &   ! array dimensions
                               i_startidx, i_endidx  ! start and end indices of loops

  REAL(KIND=wp), INTENT(IN) :: zsurf  (ie)  ! height of surface (m)
  REAL(KIND=wp), INTENT(IN) :: z1     (ie)  ! height (m) of lowest full level
  REAL(KIND=wp), INTENT(IN) :: u1     (ie)  ! zonal wind   (m/s) at lowest full level
  REAL(KIND=wp), INTENT(IN) :: v1     (ie)  ! merid.wind   (m/s) at lowest full level
  REAL(KIND=wp), INTENT(IN) :: t1     (ie)  ! temperature  (K)   at lowest full level
  REAL(KIND=wp), INTENT(IN) :: qv1    (ie)  ! spec. humidity (kg/kg) at "    "    "
  REAL(KIND=wp), INTENT(IN) :: t_g    (ie)  ! surface temperature (K)
  REAL(KIND=wp), INTENT(IN) :: qv_s   (ie)  ! surface humidity (kg/kg)
  REAL(KIND=wp), INTENT(IN) :: ps     (ie)  ! surface pressure (Pa)
  REAL(KIND=wp), INTENT(IN) :: fr_land(ie)  ! surface land fraction
  REAL(KIND=wp), INTENT(IN) :: h_ice  (ie)  ! thickness of sea ice

!     Output data

  REAL(KIND=wp), INTENT(INOUT)   :: tcm (ie) ! transfer coefficient for momentum  
  REAL(KIND=wp), INTENT(INOUT)   :: tch (ie) ! transfer coefficient for heat/moisture 
  REAL(KIND=wp), INTENT(INOUT)   :: gz0 (ie) ! roughness length * g (m2/s2)

  REAL(KIND=wp), INTENT(INOUT)   :: shfl_s(ie) ! sensible   heat flux at surface (W/m2)
  REAL(KIND=wp), INTENT(INOUT)   :: lhfl_s(ie) ! latent     heat flux at surface (W/m2)
  REAL(KIND=wp), INTENT(INOUT)   :: qhfl_s(ie) ! moisture   flux at surface (kg/m2/s)
  REAL(KIND=wp), INTENT(INOUT)   :: umfl_s(ie) ! u-momentum flux at the surface (N/m2)
  REAL(KIND=wp), INTENT(INOUT)   :: vmfl_s(ie) ! v-momentum flux at the surface (N/m2)
 
!     Local arrays and variables

  REAL(KIND=wp) :: zvpb (ie) ! wind speed in Prandtl layer
  REAL(KIND=wp) :: zx   (ie) ! utitlity variable for computation of Richardson number
  REAL(KIND=wp) :: ztcm (ie) ! transfer coefficient for momentum
  REAL(KIND=wp) :: ztch (ie) ! transfer coefficient for heat/moisture
  REAL(KIND=wp) :: zdfip(ie) ! geopotential thickness of Prandtl layer
  REAL(KIND=wp) :: zris (ie) ! Richardson number
  REAL(KIND=wp) :: zgz0m(ie) ! roughness length for momentum
  REAL(KIND=wp) :: zgz0h(ie) ! roughness length for sensible and latent heat
  LOGICAL       :: lo_ice(ie) ! logical sea ice indicator
 
 
!     local physical parameters
  REAL(KIND=wp), PARAMETER :: zah = 5.3_wp
  REAL(KIND=wp), PARAMETER :: zgz0hh = 0.98_wp     ! upper limit for roughness length for heat
  REAL(KIND=wp), PARAMETER :: zalphaf = 1._wp
  REAL(KIND=wp), PARAMETER :: zalpha0 = 0.0150_wp  ! Charnock constant for roughness length computation
                                                   !  over sea for momentum (modifified by B.Ritter 29/3/2001)
                                                   ! In COSMO alpha0= 0.0123
  REAL(KIND=wp), PARAMETER :: zalphah = 0.60_wp    !  constant for roughness length computation 
                                                   !  over sea for scalar quantities (latent and sensible heat)
                                                   ! modifified by B.Ritter 13/6/2001 (from 0.50)
  REAL(KIND=wp), PARAMETER :: zviscos = 1.5E-05_wp !  kinematic viscosity constant (m**2/s)
  REAL(KIND=wp), PARAMETER :: zbeta10 = 0.042_wp
  REAL(KIND=wp), PARAMETER :: z10 = 10._wp
  REAL(KIND=wp), PARAMETER :: zvmin = 0.01_wp      !  minimum wind velocity
! minimum value for transfer coefficients (as fractions of the neutral value for stable conditions)
  REAL(KIND=wp), PARAMETER :: ztmmin = 0.140_wp    !  minimum transfer coefficient      
  REAL(KIND=wp), PARAMETER :: zthmin = 0.010_wp    !  minimum transfer coefficient      
  REAL(KIND=wp), PARAMETER :: z0tWMO = 3.0E-2_wp   ! "WMO rougness lengths" for temperature [m]
  REAL(KIND=wp) :: zustar   !

! local utility variables
  REAL(KIND=wp) :: ztvg, ztvs, zgz0d, zgz0dd, zxi, zxih, zy, rho_s
  INTEGER       :: j1          ! loop indices

  REAL(KIND=wp), PARAMETER :: z1d3  = 1.0_wp/3.0_wp
  REAL(KIND=wp), PARAMETER :: z2d3  = 2.0_wp/3.0_wp

! OpenACC variables
  LOGICAL, OPTIONAL, INTENT(in)   :: lacc  ! GPU flag

      CALL assert_acc_device_only("parturs", lacc)

      !$ACC DATA CREATE(zvpb, zx, ztcm, ztch, zdfip, zris, zgz0m, zgz0h, lo_ice)

      IF ( msg_level >= 25) CALL message( 'mo_gme_turbdiff:', 'parturs')
 
!     wind velocity in Prandtl layer
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO j1 = i_startidx, i_endidx
        zvpb(j1) = MAX( SQRT( u1(j1)**2 +v1(j1)**2), zvmin)
        IF ( lseaice) THEN
          lo_ice(j1) = h_ice(j1) > 0._wp
        ELSE
          lo_ice(j1) = t_g(j1) < tmelt - 1.7_wp
        END IF
      END DO
      !$ACC END PARALLEL
 
!     calculation of new transfer coefficients   
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(ztvg, ztvs, zgz0d, zgz0dd, zxi, zxih, zy, rho_s, zustar)
      DO j1 = i_startidx, i_endidx
        ztvg = t_g(j1)  * (1._wp + vtmpc1*qv_s(j1)   )
        ztvs = t1 (j1)*(1._wp + vtmpc1*qv1 (j1))
        zdfip(j1)= grav*( z1(j1) - zsurf(j1) )
        zx   (j1)= ( ztvs - ztvg + rcpd*zdfip(j1) )*zdfip(j1)/t_g(j1)
        zris (j1)= zx(j1)/zvpb(j1)**2
        zx   (j1)= ABS(zx(j1))
 
!     for sea points compute initial value of roughness length, if
!     necessary (e.g. model cold start)
!     note, that sea roughness lengths for timestep 0 are generally
!     defined via the first guess forecast in the data assimilation 
        IF ( fr_land(j1) < 0.5_wp .AND. gz0(j1) <= 0.0_wp ) THEN

!         use constant value of 0.001 m for roughness length over sea ice
          IF ( lo_ice(j1) ) THEN
            gz0(j1) = 0.001_wp*grav

!         define z0 over open water                        
          ELSE
            zgz0d   = zalpha0*( zvpb(j1)/(1.0_wp/zbeta10                  &
                               +LOG( zdfip(j1)/(grav*z10) )/Rkarman) )**2
            IF ( zris(j1) < 0.0_wp ) THEN
              zgz0dd  = zalphaf**3/( zah*SQRT(zdfip(j1)) )                &
                       *( zalpha0*zx(j1) )**1.5_wp
            ELSE
              zgz0dd  = 0.0_wp
            ENDIF

!         use ECMWF method, i.e. convective velocity scale added to
!         near surface grid scale velocity
            gz0(j1) = zgz0d + zgz0dd

          ENDIF
        ENDIF     ! need to define z0 over water at initial time
 
        zgz0m(j1)=MIN( gz0(j1),0.5*zdfip(j1)) ! limit z0 for momentum

        IF ( ( fr_land(j1) < 0.5_wp) .AND.  &      ! derive z0 for heat
         &  .NOT. lo_ice(j1) ) THEN                ! over open sea
          zustar      = SQRT(zgz0m(j1)/zalpha0)    ! friction velocity derived
                                                   ! from z0 for momentum (Charnock) 
          zgz0h(j1)= grav*zviscos*zalphah/MAX(1.E-8_wp,zustar)
        ELSE 
          zgz0h(j1)= zgz0m(j1)                  ! use z0m over land and sea ice
        ENDIF
        zgz0h(j1)=MIN( zgz0h(j1),zgz0hh)        ! limit z0 for heat

!       IF ( lz0_snow .AND.( fr_land(j1) >= 0.5_wp) ) THEN
!          IF( gz0s(j1) > 0._wp)   zgz0h(j1) = MIN( gz0s(j1), zgz0hh)
!       END IF

        zxi         = zdfip(j1)/zgz0m(j1)
        zxih        = zdfip(j1)/zgz0h(j1) 
        zy          = (Rkarman/LOG(zxi))**2

        IF ( zris(j1) >= 0.0_wp ) THEN
!         stable case (land and sea)

          ztcm(j1) = zy*zvpb(j1)*MAX ( ztmmin, 1.0_wp/               &
                    (1.0_wp + 10.0_wp*zris(j1)/                      &
                     SQRT( 1.0_wp + 5.0_wp*zris(j1) ) ) ) 
          ztch(j1) = Rkarman**2/(LOG(zxi)*LOG(zxih))*zvpb(j1)*          &
                     MAX ( zthmin, 1.0_wp/(1.0_wp + 15.0_wp*zris(j1)*   &
                     SQRT( 1.0_wp + 5.0_wp*zris(j1) ) ) ) 
 
!         new z0 (for momentum) over sea
          IF ( fr_land(j1) < 0.5_wp ) THEN
            IF ( lo_ice(j1) ) THEN
!             z0=0.001 m above sea ice              
              gz0(j1) = 0.001_wp*grav
            ELSE
!             Charnock formula over open water for z0
              gz0 (j1) = zalpha0*ztcm(j1)*zvpb(j1)
            ENDIF
          ENDIF
 
        ELSE
!         unstable case
 
!         Land points (only transfer coefficients)
          IF ( fr_land(j1) >= 0.5_wp ) THEN
            ztcm(j1) = zy*zvpb(j1)*(1._wp - 10.0_wp*zris(j1)/                      &
             & (1._wp + 75.0_wp*zy*(zxi**z1d3-1.0_wp)**1.5_wp*SQRT( -zris(j1) ) ))
            ztch(j1) = Rkarman**2/(LOG(zxi)*LOG(zxih))*zvpb(j1)*                   &
             & (1._wp-15.0_wp*zris(j1)/(1._wp+75.0_wp*SQRT(zy)*Rkarman/LOG(zxih)*  &
             &  (zxih**z1d3-1.0_wp)**1.5_wp*SQRT( -zris(j1) ) ))
 
!       Sea points (transfer coefficients and z0)
          ELSE
            ztcm(j1) = zy*zvpb(j1)*(1._wp - 10.0_wp*zris(j1)/       &
             &        (1._wp + 75.0_wp*zy*SQRT( -zris(j1)*zxi ) ))
            ztch(j1) = Rkarman**2/(LOG(zxi)*LOG(zxih))*zvpb(j1)*    &
             & (1._wp - 15.0_wp*zris(j1)/(1._wp + 75.0_wp*SQRT(zy)*Rkarman/LOG(zxih)* &
             &   SQRT( -zris(j1)*zxih ) ))

!         z0=0.001 m above sea ice              
            IF ( lo_ice(j1) ) THEN
              gz0(j1) = 0.001_wp*grav
            ELSE

!         Charnock formula over open water for z0(momentum)
!         use ECMWF method, i.e. convective velocity scale added to
!         near surface grid scale velocity
              gz0(j1) = zalpha0*( ztcm(j1)*zvpb(j1)                    &
                                  +zalphaf**2*(ztch(j1)*zx(j1))**z2d3)

            ENDIF
          ENDIF
 
        ENDIF
 
!     store final results in output arrays
!     tcm and tch are dimensionless quantities
        tcm  (j1) =  ztcm(j1)/zvpb(j1)
        tch  (j1) =  ztch(j1)/zvpb(j1)
        gz0  (j1) = MAX( gz0 (j1), 1.0E-10_wp )

        rho_s     = ps(j1)/( rd * ztvg )

        shfl_s(j1) = ztch(j1)*rho_s*( cpd*(t1(j1)-t_g(j1)) + zdfip(j1) )

        lhfl_s(j1) = ztch(j1)*rho_s*( qv1(j1) - qv_s(j1) )
        qhfl_s(j1) = lhfl_s(j1)
        IF ( .NOT. lseaice) THEN
          lhfl_s(j1) = alv*lhfl_s(j1)
        ELSE IF ( h_ice(j1) > 0._wp) THEN
          lhfl_s(j1) = als*lhfl_s(j1)
        ELSE
          lhfl_s(j1) = alv*lhfl_s(j1)
        END IF

        umfl_s(j1) = -ztcm(j1)*rho_s * u1(j1)
        vmfl_s(j1) = -ztcm(j1)*rho_s * v1(j1)
        
 
      END DO
      !$ACC END PARALLEL

      IF (.NOT. using_cuda_graph) THEN
        !$ACC WAIT(1)
      END IF
      !$ACC END DATA ! zvpb, zx, ztcm, ztch, zdfip, zris, zgz0m, zgz0h, lo_ice

  END SUBROUTINE parturs

!=======================================================================

  SUBROUTINE progimp_turb ( t      , qv     , qc     , u       , v       ,  &
!                           qi     ,  qr    , qs     , o3      ,            &
                            zh     , zf     , rho    , ps      ,            &
                            tkvm   , tkvh   , t_g    , qv_s   , h_ice  ,    &
                            tcm    , tch    , ie     , ke      , ke1     ,  &
                            i_startidx, i_endidx,                           &
                            dt     , du_turb, dv_turb, dt_turb, dqv_turb ,  &
                            dqc_turb )
!!                          shfl_s , lhfl_s , qhfl_s, umfl_s , vmfl_s )
!
!=======================================================================
!
!  *progimp_turb* implicitly calculates tendencies of t qv, qc, u, v
!  due to turbulent diffusion.
!     following variables
!     - temperature                         t
!     - specific humidity                   qv
!     - specific cloud liquid water content qc
!     - zonal wind component                u
!     - meridional wind component           v
!
!  For diagnostic purposes, the turbulent fluxes of heat, moisture
!  and momentum (--> shfl_s, lhfl_s, qhfl_s, umfl_s, vmfl_s) are computed.
!
!  The code was extracted from subroutine "progexp" of GME and modified.
!
!=======================================================================

  IMPLICIT NONE
!
  INTEGER      , INTENT(IN) :: ie, ke, ke1,  &       ! dimensions
                               i_startidx, i_endidx  ! start and end indices of loops
!
  REAL(KIND=wp), INTENT(IN) :: t    (ie,ke)   ! temperature at full levels
  REAL(KIND=wp), INTENT(IN) :: qv   (ie,ke)   ! specific humidity at full levels
  REAL(KIND=wp), INTENT(IN) :: qc   (ie,ke)   ! specific cloud liquid water at full levels
  REAL(KIND=wp), INTENT(IN) :: u    (ie,ke)   ! zonal wind component at full levels
  REAL(KIND=wp), INTENT(IN) :: v    (ie,ke)   ! meridional wind component at full levels
! REAL(KIND=wp), INTENT(IN) :: qi   (ie,ke)   ! cloud ice at full levels
! REAL(KIND=wp), INTENT(IN) :: qr   (ie,ke)   ! specific rain content; optional (ntype_gsp==3)
! REAL(KIND=wp), INTENT(IN) :: qs   (ie,ke)   ! specific snow content; optional (ntype_gsp==3)
! REAL(KIND=wp), INTENT(IN) :: o3   (ie,ke)   ! ozone mass mixing ratio at full levels; 
  REAL(KIND=wp), INTENT(IN) :: zh   (ie,ke1)  ! height of half levels (m)
  REAL(KIND=wp), INTENT(IN) :: zf   (ie,ke)   ! height of full levels (m)
  REAL(KIND=wp), INTENT(IN) :: rho  (ie,ke)   ! air density at full levels (kg/m3)
!
  REAL(KIND=wp), INTENT(IN) :: ps   (ie)      ! surface pressre (Pa)
  REAL(KIND=wp), INTENT(IN) :: tkvm (ie,2:ke) ! turbulent diffusion coefficient for momentum
  REAL(KIND=wp), INTENT(IN) :: tkvh (ie,2:ke) ! turbulent diffusion coefficient for heat/moisture
!
  REAL(KIND=wp), INTENT(IN) :: t_g  (ie)      ! temperature at the interface surface - atmosphere
  REAL(KIND=wp), INTENT(IN) :: qv_s (ie)      ! specific humidity at the surface
  REAL(KIND=wp), INTENT(IN) :: h_ice(ie)      ! thickness of sea ice
  REAL(KIND=wp), INTENT(IN) :: tcm  (ie)      ! transfer coefficient for momentum at the surface
  REAL(KIND=wp), INTENT(IN) :: tch  (ie)      ! transfer coefficient for heat/moisture at the surface

  REAL(KIND=wp), INTENT(IN) :: dt             ! time step (s)
!
!=======================================================================
!
  REAL(KIND=wp), INTENT(INOUT) :: du_turb (ie,ke) ! tendency of u (m/s2)
  REAL(KIND=wp), INTENT(INOUT) :: dv_turb (ie,ke) ! tendency of v (m/s2)
  REAL(KIND=wp), INTENT(INOUT) :: dt_turb (ie,ke) ! tendency of T (K/s)
  REAL(KIND=wp), INTENT(INOUT) :: dqv_turb(ie,ke) ! tendency of qv (1/s)
  REAL(KIND=wp), INTENT(INOUT) :: dqc_turb(ie,ke) ! tendency of qc (1/s)
! REAL(KIND=wp), INTENT(INOUT) :: shfl_s  (ie)    ! sensible heat flux at the surface (W/m2)
! REAL(KIND=wp), INTENT(INOUT) :: lhfl_s  (ie)    ! latent   heat flux at the surface (W/m2)
! REAL(KIND=wp), INTENT(INOUT) :: qhfl_s  (ie)    ! moisture      flux at the surface (kg/m2/s)
! REAL(KIND=wp), INTENT(INOUT) :: umfl_s  (ie)    ! u-momentum    flux at the surface (N/m2)
! REAL(KIND=wp), INTENT(INOUT) :: vmfl_s  (ie)    ! v-momentum    flux at the surface (N/m2)
!
!=======================================================================
!
!     Local arrays

  REAL(KIND=wp) :: a1t  (ke1)      ! weight of "t+dt" for implicit treatment
  REAL(KIND=wp) :: ztmcm(ie)       ! modified transfer coefficient for momentum and heat/moisture
  REAL(KIND=wp) :: ztmkv(ie, 2:ke) ! modified vertical diffusion coefficient for momentum and heat/moisture
  REAL(KIND=wp) :: rdzrho          ! reciprocal of the layer thickness "zdp"
!
!     The following zag* arrays are used in the solution of the
!     tri-diagonal equation system resulting from the implicit 
!     treatment of diffusion
!
  REAL(KIND=wp) :: zaga       ,  &  ! terms involving the layer above
                   zagb       ,  &  ! terms involving the center layer
                   zagc(ie,ke),  &  ! terms involving the layer below
                   zag1(ie,ke),  &  ! r.h.s. of the tri-diagonal linear equation system
                   zag2(ie,ke),  &
                   zag3(ie,ke)
!
!     Local variables
!
  REAL(KIND=wp) :: zvm(ie)   ! wind speed at lowest model layer ("t-dt")
  REAL(KIND=wp) :: rho_s(ie) ! air density at surface
  REAL(KIND=wp) :: rho_h     ! air density at half levels
  REAL(KIND=wp) :: rdt       ! 1./dt  where dt is the time step
  REAL(KIND=wp) :: zagat     ! term in "zaga" proportional to vertical dif.
  REAL(KIND=wp) :: zagct     ! term in "zagc" proportional to vertical dif.
  REAL(KIND=wp) :: zt3       ! scalar temporary for optimization

  REAL(KIND=wp) :: g_o_cp  ! grav / Cpd
!
  INTEGER ::       j1, j3      ! DO loop variables

  IF ( msg_level >= 15) CALL message( 'mo_gme_turbdiff:', 'progimp_turb')
!
!=======================================================================
!
! 1.  Preliminary settings
      rdt  = 1._wp/dt
      a1t(:)    = 0.75
      a1t(ke1 ) = 1.2_wp
      a1t(ke  ) = 1.2_wp
      a1t(ke-1) = 1.1_wp
      a1t(ke-2) = 1.0_wp
      a1t(ke-3) = 0.9_wp
      a1t(ke-4) = 0.8_wp
!a1t(:) = 0._wp   ! explicit test
!   
!=======================================================================
!
      zagc(:,ke) = 0._wp
!
!=======================================================================
!
!     The tri-diagonal linear equation system to be solved for du (and dv)
!     has the general form
!
!        A(j3)*du(j3-1) + B(j3)*du(j3) + C(j3)*du(j3+1) = D(j3)
!
!     where j3 is the index of the layers.
!     
!=======================================================================
!
! 3.3 Compute the modified vertical transfer and diffusion coefficients
!
      DO j1 = i_startidx, i_endidx
!
!     Wind speed at lowest model layer (j3 = ke)
!  
        zvm(j1) = SQRT( u(j1,ke)**2 + v(j1,ke)**2)
!
!     density at the surface
!
        rho_s(j1) = ps(j1)/( rd * t_g(j1)*( 1._wp + vtmpc1*qv_s(j1) ) )
!
!     Modified transfer coefficient for momentum ---> "ztmcm"
!
        ztmcm(j1) = tcm(j1)*zvm(j1)*rho_s(j1)
!
      ENDDO
!
!     Modified vertical diffusion coefficient
!
      DO j3 = 2, ke
        DO j1 = i_startidx, i_endidx
          rho_h = 0.5_wp*( rho(j1,j3) + rho(j1,j3-1) )
          ztmkv(j1,j3) = rho_h*tkvm(j1,j3)/( zf(j1,j3) - zf(j1,j3-1) )   ! < 0
        ENDDO
      ENDDO
!
!=======================================================================
!
! 3.4 Set up the tri-diagonal linear systems for du (zag1) and dv (zag2)
! !Thomas algorithm
!=======================================================================
!
!     Top layer: j3 = 1
!
      DO j1 = i_startidx, i_endidx
        rdzrho     = 1._wp/( rho(j1,1)*( zh(j1,2)-zh(j1,1)) )
        zagct      = ztmkv(j1,2) * rdzrho
        zagc(j1,1) = -zagct*a1t(2)

        zagb       = rdt - zagc(j1,1)
        zag1(j1,1) = zagct* ( u(j1,2) - u(j1,1))
        zag2(j1,1) = zagct* ( v(j1,2) - v(j1,1))
!
        zagc(j1,1) = zagc(j1,1)/zagb
        zag1(j1,1) = zag1(j1,1)/zagb
        zag2(j1,1) = zag2(j1,1)/zagb

      ENDDO
!
!       The layers between j3 = 2 and ke-1
!
      DO j3 = 2, ke-1
        DO j1 = i_startidx, i_endidx
          rdzrho      = 1._wp/( rho(j1,j3)*( zh(j1,j3+1)-zh(j1,j3)) )
          zagat       = ztmkv (j1,j3)  *rdzrho          
          zagct       = ztmkv (j1,j3+1)*rdzrho          
          zaga        = -zagat*a1t(j3)
          zagc(j1,j3) = -zagct*a1t(j3+1)
          zagb        = rdt - zaga - zagc(j1,j3)

          zag1(j1,j3) =  zagat*( u(j1,j3-1) - u(j1,j3)) &
            &          + zagct*( u(j1,j3+1) - u(j1,j3))
          zag2(j1,j3) =  zagat*( v(j1,j3-1) - v(j1,j3)) &
            &          + zagct*( v(j1,j3+1) - v(j1,j3))
          zt3         = 1._wp/( zagb - zaga*zagc(j1,j3-1))
!
          zagc(j1,j3) =  zagc(j1,j3)*zt3
          zag1(j1,j3) = ( zag1(j1,j3) - zaga*zag1(j1,j3-1) )*zt3
          zag2(j1,j3) = ( zag2(j1,j3) - zaga*zag2(j1,j3-1) )*zt3
        ENDDO
      ENDDO
!
!     Lowest model layer (j3 = ke)
!
      DO j1 = i_startidx, i_endidx
        rdzrho      = 1._wp/( rho(j1,ke)*( zh(j1,ke1)-zh(j1,ke)) )
        zagat       = ztmkv (j1,ke)*rdzrho
        zagct       = ztmcm (j1)   *rdzrho
        zaga        = - zagat*a1t(ke)
        zagb        = rdt - zaga - zagct*a1t(ke1)
        zag1(j1,ke) =  zagat*( u(j1,ke-1) - u(j1,ke)) &
          &          + zagct*  u(j1,ke)
        zag2(j1,ke) =  zagat*( v(j1,ke-1) - v(j1,ke)) &
          &          + zagct*  v(j1,ke)

        du_turb(j1,ke) = ( zag1(j1,ke) - zaga*zag1(j1,ke-1))/  &
          &              ( zagb        - zaga*zagc(j1,ke-1))
        dv_turb(j1,ke) = ( zag2(j1,ke) - zaga*zag2(j1,ke-1))/  &
          &              ( zagb        - zaga*zagc(j1,ke-1))
      ENDDO
!
!=======================================================================
!
! 3.5 The tri-diagonal linear system is solved by Gaussian elimination
!     The matrices are transformed to an upper triangular form with
!     unit diagonal
!
      DO j3 = ke-1, 1, -1
        DO j1 = i_startidx, i_endidx
          du_turb(j1,j3) = zag1(j1,j3) - zagc(j1,j3)* du_turb(j1,j3+1)
          dv_turb(j1,j3) = zag2(j1,j3) - zagc(j1,j3)* dv_turb(j1,j3+1)
        ENDDO
      ENDDO
!
!=======================================================================
!
! 3.6 Turbulent momentum fluxes
!
!     DO j1 = i_startidx, i_endidx
!       umfl_s(j1) = -ztmcm(j1) * ( du_turb(j1,ke) + u(j1,ke) )
!       vmfl_s(j1) = -ztmcm(j1) * ( dv_turb(j1,ke) + v(j1,ke) )
!     ENDDO
! 
!     Calculate tendencies du_turb = du/dt, dv_turb = dv/dt
      du_turb(:,:) = rdt*du_turb(:,:)
      dv_turb(:,:) = rdt*dv_turb(:,:)
!
!=======================================================================
!
! 4.  Calculation of the temperature (T), specific humidity (water
!     vapour, qv), liquid cloud water content.
!
!     Since the vertical diffusion is
!     solved implicitly for stability reasons, the complete solution
!     is the result of a tri-diagonal linear system connecting the
!     vertical layers.
!
!=======================================================================
!
!     The computation of the equations is performed row by row to save
!     local storage. The tri-diagonal linear equation system to be
!     solved for dT ( dqv, dqc) has the general form
!
!        A(j3)*dT(j3-1) + B(j3)*dT(j3) + C(j3)*dT(j3+1) = D(j3)
!
!     where j3 is the index of the layers. Uses Thomas Algorithm
!     
!=======================================================================
!
! 4.7 Compute the modified vertical transfer and diffusion coefficients
!
      DO j1 = i_startidx, i_endidx
!
!     Modified transfer coefficient for heat/moisture ---> "ztmcm"
!
        ztmcm(j1) = tch(j1)*zvm(j1)*rho_s(j1)
!
      ENDDO
!
!     Modified vertical diffusion coefficient
!
      DO j3 = 2, ke
        DO j1 = i_startidx, i_endidx
          rho_h = 0.5_wp*( rho(j1,j3) + rho(j1,j3-1) )
          ztmkv(j1,j3) = rho_h*tkvh(j1,j3)/( zf(j1,j3) - zf(j1,j3-1) )
        ENDDO
      ENDDO
      g_o_cp = grav*rcpd
!
!=======================================================================
!
! 4.8 Set up the tri-diagonal linear systems for t (zag1),
!     qv (zag2), qc (zag3), qi (zag4), o3 (zag5), qr(zag6), and qs(zag7)
!
!=======================================================================
!
!     Top layer (j3 = 1)
!
      DO j1 = i_startidx, i_endidx
!
        rdzrho     = 1._wp/( rho(j1,1)*( zh(j1,2)-zh(j1,1)) )
        zagct      = ztmkv(j1,2)*rdzrho
        zagc(j1,1) = -zagct*a1t(2)

        zagb       = rdt - zagc(j1,1)
        zag1(j1,1) = zagct*( t (j1,2) - t (j1,1) + g_o_cp*( zf(j1,2)-zf(j1,1)) )
        zag2(j1,1) = zagct*( qv(j1,2) - qv(j1,1))
        zag3(j1,1) = zagct*( qc(j1,2) - qc(j1,1))
!       zag4(j1,1) = zagct*( qi(j1,2) - qi(j1,1))
!
        zagc(j1,1) = zagc(j1,1)/zagb
        zag1(j1,1) = zag1(j1,1)/zagb         
        zag2(j1,1) = zag2(j1,1)/zagb         
        zag3(j1,1) = zag3(j1,1)/zagb         
!       zag4(j1,1) = zag4(j1,1)/zagb         
!
      ENDDO
!
!     The layers between j3 = 2 and ke-1
!
      DO j3 = 2, ke-1
        DO j1 = i_startidx, i_endidx
          rdzrho      = 1._wp/( rho(j1,j3)*( zh(j1,j3+1)-zh(j1,j3)) )
          zagat       = ztmkv(j1,j3)  *rdzrho          
          zagct       = ztmkv(j1,j3+1)*rdzrho          
          zaga        = -zagat*a1t(j3)
          zagc(j1,j3) = -zagct*a1t(j3+1)
          zagb        = rdt - zaga - zagc(j1,j3)

          zag1(j1,j3) =  zagat*( t (j1,j3-1) - t (j1,j3) + g_o_cp*( zf(j1,j3-1)-zf(j1,j3))) &
            &          + zagct*( t (j1,j3+1) - t (j1,j3) + g_o_cp*( zf(j1,j3+1)-zf(j1,j3)))
          zag2(j1,j3) =  zagat*( qv(j1,j3-1) - qv(j1,j3))  &
            &          + zagct*( qv(j1,j3+1) - qv(j1,j3))
          zag3(j1,j3) =  zagat*( qc(j1,j3-1) - qc(j1,j3))  &
            &          + zagct*( qc(j1,j3+1) - qc(j1,j3)) 
!         zag4(j1,j3) =  zagat*( qi(j1,j3-1) - qi(j1,j3))  &
!           &          + zagct*( qi(j1,j3+1) - qi(j1,j3)) 
!
          zt3         = 1._wp/( zagb - zaga*zagc(j1,j3-1))
!
          zagc(j1,j3) =  zagc(j1,j3)*zt3
          zag1(j1,j3) = ( zag1(j1,j3) - zaga*zag1(j1,j3-1) )*zt3
          zag2(j1,j3) = ( zag2(j1,j3) - zaga*zag2(j1,j3-1) )*zt3
          zag3(j1,j3) = ( zag3(j1,j3) - zaga*zag3(j1,j3-1) )*zt3
!         zag4(j1,j3) = ( zag4(j1,j3) - zaga*zag4(j1,j3-1) )*zt3
        ENDDO
      ENDDO
!
!       Lowest model layer (j3 = ke)
!

      DO j1 = i_startidx, i_endidx
!
        rdzrho      = 1._wp/( rho(j1,ke)*( zh(j1,ke1)-zh(j1,ke)) )
        zagat       = ztmkv(j1,ke)*rdzrho
        zagct       = ztmcm (j1)  *rdzrho
        zaga        = -zagat*a1t(ke)
        zagb        = rdt - zaga - zagct*a1t(ke1)           
        zag1(j1,ke) =  zagat*( t (j1,ke-1) - t(j1,ke) + g_o_cp*( zf(j1,ke-1)-zf(j1,ke)) ) &
          &          - zagct*( t_g(j1    ) - t(j1,ke) + g_o_cp*( zh(j1,ke1 )-zf(j1,ke)) )
        zag2(j1,ke) =  zagat*( qv(j1,ke-1) - qv(j1,ke)) &
          &          - zagct*( qv_s(j1   ) - qv(j1,ke)) 
        zag3(j1,ke) =  zagat*( qc(j1,ke-1) - qc(j1,ke)) &
          &          + zagct*  qc(j1,ke)
!       zag4(j1,ke) =  zagat*( qi(j1,ke-1) - qi(j1,ke)) &
!         &          - zagct*  qi(j1,ke)

        dt_turb (j1,ke) = ( zag1(j1,ke) - zaga*zag1(j1,ke-1))/  &
          &               ( zagb        - zaga*zagc(j1,ke-1))
        dqv_turb(j1,ke) = ( zag2(j1,ke) - zaga*zag2(j1,ke-1))/  &
          &               ( zagb        - zaga*zagc(j1,ke-1))
        dqc_turb(j1,ke) = ( zag3(j1,ke) - zaga*zag3(j1,ke-1))/  &
          &               ( zagb        - zaga*zagc(j1,ke-1))
!       dqi_turb(j1,ke) = ( zag4(j1,ke) - zaga*zag4(j1,ke-1))/  &
!         &               ( zagb        - zaga*zagc(j1,ke-1))
!
      ENDDO
!
!=======================================================================
!
! 4.9 The tri-diagonal linear system is solved by Gaussian elimination
!     The matrices are transformed to an upper triangular form with
!     unit diagonal
!
      DO j3 = ke-1, 1,-1
        DO j1 = i_startidx, i_endidx
          dt_turb (j1,j3) = zag1(j1,j3) - zagc(j1,j3)*dt_turb (j1,j3+1)
          dqv_turb(j1,j3) = zag2(j1,j3) - zagc(j1,j3)*dqv_turb(j1,j3+1)
          dqc_turb(j1,j3) = zag3(j1,j3) - zagc(j1,j3)*dqc_turb(j1,j3+1)
!         dqi_trub(j1,j3) = zag4(j1,j3) - zagc(j1,j3)*dqi_turb(j1,j3+1)
        ENDDO
      ENDDO
!
!=======================================================================
!
! 4.10 Sensible and latent heat fluxes
!
!     DO j1 = i_startidx, i_endidx
!
!       shfl_s(j1) = - ztmcm(j1)*cpd*                          &
!         &           ( t_g(j1) - t(j1,ke) -dt_turb(j1,ke)     &
!         &            + g_o_cp*(zh(j1,ke1)-zf(j1,ke) ) )
!
!       lhfl_s(j1) = - ztmcm(j1)*( qv_s(j1) - qv(j1,ke) - dqv_turb(j1,ke) )
!       qhfl_s(j1) = lhfl_s(j1)
!       IF ( .NOT. lseaice) THEN
!         lhfl_s(j1) = alv*lhfl_s(j1)
!       ELSE IF ( h_ice(j1) > 0._wp) THEN
!         lhfl_s(j1) = als*lhfl_s(j1)
!       ELSE
!         lhfl_s(j1) = alv*lhfl_s(j1)
!       END IF
!
!     ENDDO
!
!     Calculate tendencies dT_turb = dT/dt, dqv_turb = dqv/dt, dqc_turb = dqc/dt
      dt_turb (:,:) = rdt*dt_turb (:,:)
      dqv_turb(:,:) = rdt*dqv_turb(:,:)
      dqc_turb(:,:) = rdt*dqc_turb(:,:)
!
  END SUBROUTINE progimp_turb

!==============================================================================

  SUBROUTINE nearsfc( t     , qv    , u     , v     , zf   ,  &
               &      ps    , t_g   , tcm   , tch   , gz0  ,  &  ! gz0s, 
               &      shfl_s, lhfl_s, umfl_s, vmfl_s, zsurf,  &
               &      fr_land,pf1   , qv_s  , ie    , ke   ,  &
               &      i_startidx, i_endidx,                   &
               &      t_2m  , qv_2m , td_2m , rh_2m ,         &
               &      u_10m , v_10m, lacc )
!
!=======================================================================
!
!     *nearsfc* computes the near surface fields:
!     t_2m:      temperature               at  2 m above ground
!     qv_2m:     specific humidity         at  2 m above ground
!     td_2m:     dew point temperature     at  2 m above ground
!     rh_2m:     relative humidity         at  2 m above ground
!     u_10m:     zonal wind component      at 10 m above ground
!     v_10m:     meridional wind component at 10 m above ground
!     
!==============================================================================

  USE mo_convect_tables,     ONLY: b1  => c1es,     & ! variables for computing the saturation steam pressure
                               &   b2w => c3les,    & ! over water (w)
                               &   b4w => c4les

  IMPLICIT NONE
!
!=======================================================================
!
! array dimensions
  INTEGER, INTENT(IN) :: ie, ke,   &         ! horizontal, vertical
                         i_startidx, i_endidx  ! start and end indices of loops in horizontal

!=======================================================================
!
  REAL(KIND=wp), INTENT(IN) ::  &
               t    (ie, ke),  & ! temperature at full levels
               qv   (ie, ke),  & ! specific humidity at full levels
               u    (ie, ke),  & ! zonal wind component at full levels
               v    (ie, ke),  & ! meridional wind component at full levels
               zf   (ie, ke)     ! height of full levels (m)
!
  REAL(KIND=wp), INTENT(IN) :: &
               ps   (ie),  & ! surface pressure
               tcm  (ie),  & ! transfer coefficient for momentum at the surface
               tch  (ie),  & ! transfer coefficient for heat/moisture at the surface
               gz0  (ie),  & ! roughness length (*G)
!              gz0s (ie),  & ! roughness length (*G) for snow
               zsurf(ie),  & ! orography (m)
               fr_land(ie),& ! land fraction (m)
               pf1  (ie),  & ! pressure at first full level above ground (Pa)
               shfl_s(ie), & ! sensible heat flux (W/m2)
               lhfl_s(ie), & ! latent heat flux (W/m2)
               umfl_s(ie), & ! u-momentum flux at surface (N/m2)
               vmfl_s(ie)    ! v-momentum flux at surface (N/m2)
!
! LOGICAL, INTENT(IN) ::  lolp (ie) ! logical switch indicating land (.true.) or water point
!
!=======================================================================
!
      REAL(KIND=wp), INTENT(IN) ::  &
               t_g  (ie),  & ! temperature at the interface surface - atmosphere
               qv_s (ie)     ! specific humidity at the surface
!
!=======================================================================
!     Output:
!
  REAL(KIND=wp), INTENT(INOUT) :: &
               t_2m (ie),  & ! temperature at 2 m above ground
               qv_2m(ie),  & ! specific humidity at 2 m above ground
               td_2m(ie),  & ! dew point temperature at 2 m above ground
               rh_2m(ie),  & ! relative humidity     at 2 m above ground
               u_10m(ie),  & ! zonal wind component at 10 m above ground
               v_10m(ie)     ! meridional wind component at 10 m above ground
!
!=======================================================================
!
!     Local arrays
!
  REAL(KIND=wp) ::          &
               zv   (ie),  & ! wind speed in the Prandtl layer (j3 = i3e)
               zcm  (ie),  & ! turbulent transfer coefficient for momentum
               zch  (ie),  & ! turbulent transfer coefficient for heat/moisture
               zgz0 (ie)     ! roughness length (*G)
!              zgz0s(ie)     ! roughness length (*G) for snow
!
!     Local variables
!
      REAL(KIND=wp) :: &
               zvmin,      & ! minimum wind speed
               z10g,       & ! geopotential at a height of 10 m
               z2g,        & ! geopotential at a height of  2 m
               zh_s,       & ! dry static energy at the surface
               zh_p,       & ! dry static energy at the Prandtl layer
               zfi_p,      & ! height of the Prandtl layer above ground
               zfi_2,      & ! height of the 2. layer above ground
               zsqcm,      & ! square root of turbulent transfer coeff.
               zchdcm,     & ! stability dependent factor for 10 m wind
               zlnz,       & ! stability dependent factor for 10 m wind
               zh_05m,     & ! dry static energy at 0.5 m (empirical)
               zh_2m,      & ! dry static energy at 2 m
               zp_2m,      & ! pressure at 2 m
               zp_p,       & ! pressure at the top of the Prandtl layer
               zpsat_2m,   & ! saturation pressure with respect to water 
!                              at a height of 2 m above ground
               zpsat_p,    & ! saturation pressure with respect to water 
!                              at the height of the Prandtl layer
               zqvs_2m,    & ! specific humidity  at a height of 2 m at
!                              saturation
               zqvs_p,     & ! specific humidity  at the height of the
!                              Prandtl layer at saturation
               zqv_2m,     & ! specific humidity  at a height of 2 m
               zqv_2m_min, & ! minimum specific humidity  at a height of 2 m
!                              (dependent on the minimum relative humidity)
               zfrac         ! factor for the dew point temperature at 2 m
!

!_dm>
!  Local variables needed for the new diagnostics of T2m
      REAL(KIND=wp) :: &
               rho_s     , & ! Surface density [kg m^(-3)]
               ufr_s     , & ! Surface firction velocity [m s^(-1)]
               psiMO_tp  , & ! The Monin-Obukhov (MO) stability function for temperature profile 
                             ! relative to the first model layer above the surface [-]
               psiMO_t2  , & ! The MO stability function for temperature profile 
                             ! relative to the 2m above the surface [-]
               lOburecip , & ! Reciprocal of the Obukhov length [m^(-1)]
               z0tWMO_g  , & ! "WMO rougness lengths" for temperature times G [m^2 s^(-2)]
               rho_srec      ! 1/rho_s
!_dm<

!
      INTEGER ::  j1        ! DO loop variables
!

!_dm>
!  Local parameters needed for the new diagnostics of T2m
      REAL (KIND = wp), PARAMETER :: &
        z0uthr     = 3.0E-2_wp     , & ! Threshold value of the rougness length for wind velocity
                                       ! used to determine z_0 for temperature in the T2m calculations [m]
                                       ! (adopted from ECMWF IFS documentation CY31r1)
        z0tWMO     = 3.0E-2_wp     , & ! "WMO rougness lengths" for temperature [m]
                                       ! (values adopted from ECMWF IFS documentation CY31r1 times 10)
        cMOtstsfac = 5.0_wp        , & ! Factor in the MO flux-profile relationships
                                       ! (temperature, stable stratification) [-]
        cMOtconfac = 15.0_wp       , & ! Factor in the MO flux-profile relationships
                                       ! (temperature, convection) [-]
        cMOtconexp = 0.5_wp        , & ! Exponent in the MO flux-profile relationships
                                       ! (temperature, convection) [-]
        small_num  = 1.0E-04_wp        ! A small number (security constant)
!
!=======================================================================
!
  LOGICAL, OPTIONAL, INTENT(in)   :: lacc  ! GPU flag
  LOGICAL                         :: lzacc ! non-optional version of lacc

      CALL set_acc_host_or_device(lzacc, lacc)

      !$ACC DATA CREATE(zv, zcm, zch, zgz0) IF(lzacc)

!_dm<
!
!=======================================================================
!
! 1.  Preliminary settings
!   
!     Minimum wind speed
      zvmin   = 0.01_wp
!
!     Geopotential at a height of 2 and 10 m above ground
      z10g    = 10._wp*grav
      z2g     =  2._wp*grav
!
!=======================================================================
!
! 3.  Compute the wind speed in the Prandtl layer, i.e. the lowest model
!     layer (j3 = ke)     ---> zv
!     Set minimum values of the turbulent transfer coefficients for
!     momentum              ---> zcm
!     heat/moisture         ---> zch
!     roughness length      ---> zgz0
!
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO j1 = i_startidx, i_endidx
        zv  (j1) = MAX ( SQRT (u(j1,ke)**2+v(j1,ke)**2), zvmin )
        zcm (j1) = MAX ( tcm(j1), 5.E-4_wp )
        zch (j1) = MAX ( tch(j1), 4.E-5_wp )
        zgz0(j1) = MAX ( gz0(j1), 1.E-3_wp )
      ENDDO
      !$ACC END PARALLEL
!
!=======================================================================
!
! 4.  Compute the temperature and dew point temperature at 2 m, and the
!     wind components at 10 m
!
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR PRIVATE(zh_s, zh_p, zfi_p, zfi_2, zsqcm, zchdcm, zlnz, zh_05m) &
      !$ACC   PRIVATE(zh_2m, zp_2m, zp_p, zpsat_2m, zpsat_p, zqvs_2m, zqvs_p, zqv_2m, zqv_2m_min, zfrac) &
      !$ACC   PRIVATE(rho_s, ufr_s, psiMO_tp, psiMO_t2, lOburecip, z0tWMO_g, rho_srec)
      DO j1 = i_startidx, i_endidx
!
!     Dry static energy at the surface            ---> zh_s
!     Dry static energy at the Prandtl layer      ---> zh_p
!
        zh_s    = cpd*t_g(j1   ) + grav*zsurf(j1   )
        zh_p    = cpd*t  (j1,ke) + grav* zf  (j1,ke)
!
!     Height (*G) of the Prandtl layer            ---> zfi_p
!
        zfi_p   = grav*( zf(j1,ke) - zsurf(j1) )
        zfi_2   = grav*( zf(j1,ke-1)-zsurf(j1) )
!
!     Stability-dependent factors for the derivation of the 2 m and 10 m values
!
        zsqcm   = SQRT ( zcm(j1) )
        zchdcm  = zch(j1)/(Rkarman*zsqcm)
!
!     10 m values
!
        IF ( zfi_p > z10g ) THEN
!
!     Derive the 10 m wind by inverting the flux-
!     gradient relation in the Prandtl layer (assuming constant fluxes 
!     of momentum and heat)
!
!     Stable stratification
!
          IF ( zh_s <= zh_p ) THEN
!
            zlnz = LOG((z10g  + zgz0(j1))/zgz0(j1)) - z10g/zfi_p* &
                   LOG((zfi_p + zgz0(j1))/zgz0(j1))
!
!     Wind components at 10 m
!
            u_10m(j1) = u(j1,ke)*( z10g/zfi_p + zsqcm/Rkarman*zlnz )
            v_10m(j1) = v(j1,ke)*( z10g/zfi_p + zsqcm/Rkarman*zlnz )
!
          ELSE
!
!     Set minimum values of the stability-dependent factors
!
            zsqcm = MAX ( zsqcm , 0.01_wp )
!
            zlnz  = 1._wp - zsqcm/Rkarman*LOG( 1._wp +       &
              &         ( EXP(Rkarman/zsqcm) - 1._wp )*      &
              &                 zgz0(j1)*(zfi_p - z10g)/     &
              &               (zfi_p*(z10g + zgz0(j1))) )
!
!     Wind components at 10 m
!
            u_10m(j1) = u(j1,ke)*zlnz
            v_10m(j1) = v(j1,ke)*zlnz
!
          ENDIF

        ELSE
!
!     Prandl layer less than 10 m. Interpolate wind components at 10 m
!
          u_10m(j1) = ( u(j1,ke)*(zfi_2-z10g) + u(j1,ke-1)*(z10g-zfi_p) )/( zfi_2-zfi_p )
          v_10m(j1) = ( v(j1,ke)*(zfi_2-z10g) + v(j1,ke-1)*(z10g-zfi_p) )/( zfi_2-zfi_p )

        ENDIF
!
!     2 m values
!
        IF ( zfi_p > z2g ) THEN
!
!     Derive the 2 m temperature by inverting the flux-
!     gradient relation in the Prandtl layer (assuming constant fluxes 
!     of momentum and heat)
!

!_dm>
!  A "new" procedure to diagnose T2m over land is used.
!  A complete form of the MO surface-layer flux-profile relationships is applied
!  [not the approximations due to J.-F. Geleyn (1988), or the like]
!  in case the height of the first model level above the ground exceeds 2 m
!  (otherwise, a simple linear inerpolation is used; see below).
!  At the moment, this new procedure is only used to diagnose T2m.
!  It should later be applied to diagnose the 2m humidity and the 10m wind
!  as well as to compute the surface fluxes in "parturs.f90".
!  Note that the surface layer is considered to be stable/unstable depending on
!  whether dry static energy at the first model level above the surface, zh_p,
!  is larger/smaller than the dry static energy at the surface, zh_s.
!  This should be reformulated in terms of moist static energy,
!  perhaps with due regard for the water loading effect.
!
!  Over water, an old procedure to diagnose T2m is used.
!Water_or_Land: IF( .NOT. lolp(j1) ) THEN 
Water_or_Land: IF( fr_land(j1) < 0.5_wp ) THEN 
!_dm<

!     Stable stratification
!
            IF ( zh_s <= zh_p ) THEN
!
!     For the temperature, an empirical correction is applied to avoid
!     unrealistic minimum values in very stable situations
!
!     Dry static energy at 0.5 m (empirical value)
!
              zh_05m = zh_s + 0.25_wp*(zh_p - zh_s)
!
!     Dry static energy at 2 m with empirical correction
!
              zh_2m  = zh_05m + 1.5_wp*grav*(zh_p - zh_05m)/(zfi_p - 0.5_wp*grav)
!
!     Unstable stratification
!
            ELSE
!
!     Set minimum values of the stability-dependent factors
!
              zchdcm = MAX ( zchdcm, 0.05_wp )
!
!     Dry static energy at 2 m
!
              zh_2m = zh_s + (zh_p - zh_s)*(1._wp - zchdcm*LOG (1._wp + &
              &  (EXP (1._wp/zchdcm) - 1._wp)*(zgz0(j1)*(zfi_p - z2g)/  &
              &  (zfi_p*(z2g  + zgz0(j1))))) )
!
            ENDIF
!

!_dm>
!  Over land, a new  procedure to diagnose T2m is used.
          ELSE Water_or_Land
!  Compute quantities required in both stable and convective conditions
!  The "WMO rougness length" with respect to the temperature times G
            IF( zgz0(j1) < z0uthr*grav) THEN
              z0tWMO_g = zgz0(j1)
            ELSE
              z0tWMO_g = z0tWMO*grav
            ENDIF
!           IF ( lz0_snow) THEN
!             z0tWMO_g = MIN( zgz0s(j1), z0tWMO_g)
!           END IF
!  Air density at the surface
            rho_s = t_g(j1)*(1._wp+vtmpc1*qv_s(j1))
            rho_s = ps(j1)/(rd*MAX(rho_s,small_num))
            rho_srec = 1._wp/MAX( rho_s, small_num)
!  Surface friction velocity
            ufr_s = SQRT( umfl_s(j1)*umfl_s(j1) + vmfl_s(j1)*vmfl_s(j1) )
            ufr_s = SQRT( ufr_s * rho_srec)
!  Surface buoyancy flux times 1/G
            lOburecip = ( shfl_s(j1)*rcpd/MAX(t_g(j1),small_num) +   &
                          vtmpc1*lhfl_s(j1)/alv ) * rho_srec
!  Reciprocal of the Obukhov length times 1/G
            lOburecip = lOburecip*Rkarman/MAX(ufr_s,small_num)**3

            IF ( zh_s <= zh_p ) THEN
!  Stable stratification
!  The Moni-Obukhov stability functions
              IF(lOburecip <= 0._wp) THEN
!  The Obukhov length proves to be negative, indicating unstable stratification
!  (this is unlikely but it may happen as the stability criterium used above is based upon 
!  the dry static energy difference between the surface and the first model layer,
!  not upon the sign of the surface buoyancy flux).
!  Set the MO stability functions to zero and interpolate along the logarithmic curve.
                psiMO_tp = 0._wp
                psiMO_t2 = 0._wp
              ELSE
                psiMO_tp = cMOtstsfac*zfi_p*lOburecip*(1._wp-MIN(z0tWMO_g/zfi_p, 1._wp))
                psiMO_t2 = cMOtstsfac*z2g  *lOburecip*(1._wp-MIN(z0tWMO_g/z2g  , 1._wp))
              ENDIF

            ELSE
!  Unstable stratification
!  The Moni-Obukhov stability functions
              IF(lOburecip .GE. 0._wp) THEN
!  The Obukhov length proves to be positive, indicating stable stratification
!  (this is unlikely but it may happen as the stability criterium used above is based upon 
!  the dry static energy difference between the surface and the first model layer,
!  not upon the sign of the surface buoyancy flux).
!  Set the MO stability functions to zero and interpolate along the logarithmic curve.
                psiMO_tp = 0._wp
                psiMO_t2 = 0._wp
              ELSE
                psiMO_t2 = (1._wp-cMOtconfac*z0tWMO_g*lOburecip)**cMOtconexp
                psiMO_tp = (1._wp-cMOtconfac*zfi_p   *lOburecip)**cMOtconexp
                psiMO_tp = 2._wp*LOG((1._wp+psiMO_t2)/(1._wp+psiMO_tp))
                psiMO_t2 = 2._wp*LOG((1._wp+psiMO_t2)/   &
                  &       (1._wp+(1._wp-cMOtconfac*z2g*lOburecip)**cMOtconexp))
              ENDIF
            ENDIF

!  Dry static energy at 2m
            zh_2m = LOG(MAX(zfi_p/z0tWMO_g,1._wp)) + psiMO_tp
            zh_2m = MAX(zh_2m,small_num)
            zh_2m = ( LOG(MAX(z2g/z0tWMO_g,1._wp)) + psiMO_t2 ) / zh_2m
            zh_2m = MIN(1._wp,MAX(zh_2m,0._wp))
            zh_2m = zh_s+(zh_p-zh_s)*zh_2m

!
          ENDIF Water_or_Land
!_dm>

!     Temperature at 2 m based on the dry static energy at 2 m
!
          t_2m(j1) = (zh_2m - z2g - grav*zsurf(j1))*rcpd

        ELSE
!
!       First level less than 2 m: Interpolate temperature at 2 m
!
          t_2m(j1) = ( t(j1,ke)*(zfi_2-z2g )                  &
                       +t(j1,ke-1)*( z2g-zfi_p) )/( zfi_2-zfi_p )

        ENDIF
!
!=======================================================================
!
! 5.  The dew point temperature at 2m is derived under the assumption
!     of a constant relative humidity in the Prandtl layer
!
!     Pressure at a height of 2 m above ground
!
        zp_2m = ps(j1)*(1._wp - z2g/(rd*t_2m(j1)       &
                          *(1._wp + vtmpc1*qv(j1,ke))))
!
!     Pressure at the top of the Prandtl layer, i.e. the lowest model
!     layer
!
        zp_p  = pf1(j1)
!
!     Saturation water vapour pressure with respect to water at 2 m
!     above ground and at the lowest model layer
!
        zpsat_2m = sat_pres_water ( t_2m(j1) )
        zpsat_p  = sat_pres_water ( t   (j1,ke) )
!
!     Specific humidity at saturation at 2 m and at the lowest model
!     layer
!
        zqvs_2m = spec_humi( zpsat_2m, zp_2m)
        zqvs_p  = spec_humi( zpsat_p , zp_p )
!
!     Specific humidity at 2 m above ground under the assumption of
!     constant relative humidity in the Prandtl layer
!
        zqv_2m_min = 0.05_wp*zqvs_2m ! set the minimum specific humidity in 2 m 
                                     ! according to a minimum relative humidity of 5%
        zqv_2m     = MAX ( zqv_2m_min, qv(j1,ke)*zqvs_2m/zqvs_p )
        qv_2m(j1) = zqv_2m
!
!     Dew point temperature at 2 m above the ground
!
        zfrac     = LOG (zp_2m*zqv_2m/(b1*(rdv + O_m_rdv*zqv_2m)))
        td_2m(j1) = (b2w*tmelt - b4w*zfrac)/(b2w - zfrac)
        td_2m(j1) = MIN ( td_2m(j1), t_2m(j1) )

        rh_2m(j1) = sat_pres_water( td_2m(j1) )/zpsat_p
!
      ENDDO
      !$ACC END PARALLEL

      IF (.NOT. using_cuda_graph) THEN
        !$ACC WAIT(1)
      END IF
      !$ACC END DATA ! zv, zcm, zch, zgz0

  END SUBROUTINE nearsfc

END MODULE mo_gme_turbdiff
