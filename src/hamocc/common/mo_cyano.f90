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

! @brief N2 fixation and cyanobateria dynamics.
!
! Contains computation of diagnostic and prognostic N2 fixation and.
! cyanobacteria dyanmics (growth, decay, buouyancy)

MODULE mo_cyano

  USE mo_kind, ONLY           : wp
  USE mo_control_bgc, ONLY    : dtb, dtbgc, bgc_nproma, bgc_zlevs
  USE mo_bgc_memory_types, ONLY  : t_bgc_memory
  USE mo_fortran_tools, ONLY  : set_acc_host_or_device
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cyano, &
            cyadyn

CONTAINS

SUBROUTINE cyano (local_bgc_mem, start_idx,end_idx,pddpo, za, lacc)

!! @brief diagostic N2 fixation

  USE mo_memory_bgc, ONLY     : rnit, n2_fixation, rn2
  USE mo_param1_bgc, ONLY     : iano3, iphosph, igasnit, &
       &                        ioxygen, ialkali, knfixd, &
       &                        kn2b, kaou

  IMPLICIT NONE

  !! Arguments
  TYPE(t_bgc_memory), POINTER    :: local_bgc_mem


  INTEGER, INTENT(in) :: start_idx            !< start index for j loop (ICON cells, MPIOM lat dir)   
  INTEGER, INTENT(in) :: end_idx              !< end index  for j loop  (ICON cells, MPIOM lat dir) 

  REAL(wp), INTENT(in) :: pddpo(bgc_nproma,bgc_zlevs) !< size of scalar grid cell (3rd dimension) [m].
  REAL(wp), INTENT(in) :: za(bgc_nproma)              !< surface height [m].
  LOGICAL, INTENT(IN), OPTIONAL :: lacc

  !! Local variables

  REAL(wp) :: oldnitrate
  INTEGER  :: j
  LOGICAL :: lzacc

  CALL set_acc_host_or_device(lzacc, lacc)
  !
  ! --------------------------------------------------------------------
  !
     !$ACC DATA COPY(local_bgc_mem, local_bgc_mem%bgctra, local_bgc_mem%bgcflux) &
     !$ACC   COPY(local_bgc_mem%bgctend, local_bgc_mem%satoxy) IF(lzacc)
     !$ACC PARALLEL ASYNC(1) IF(lzacc)
     !$ACC LOOP GANG VECTOR
     DO j = start_idx, end_idx

        IF (pddpo(j,1) > EPSILON(0.5_wp)) THEN
           IF (local_bgc_mem%bgctra(j,1,iano3) < (rnit*local_bgc_mem%bgctra(j,1,iphosph))) THEN

              oldnitrate = local_bgc_mem%bgctra(j,1,iano3)

              local_bgc_mem%bgctra(j,1,iano3) = local_bgc_mem%bgctra(j,1,iano3)                 &
                                     * (1._wp - n2_fixation*dtb)              &
                                     + n2_fixation * dtb * rnit * local_bgc_mem%bgctra(j,1,iphosph)


              local_bgc_mem%bgctra(j,1,igasnit)  = local_bgc_mem%bgctra(j,1,igasnit)                  &
                                     - (local_bgc_mem%bgctra(j,1,iano3)-oldnitrate) /rn2  


         
              local_bgc_mem%bgctra(j,1,ioxygen)  = local_bgc_mem%bgctra(j,1,ioxygen)                  &
                                     - (local_bgc_mem%bgctra(j,1,iano3)-oldnitrate) * 1.5_wp  ! factor 1.5 results from: NO3 = N + 3O = N + 1.5*O2 
                                                                                 ! (O2 is used for the balancing)

              !alkalinity is decreased during n-fixation due to release of H+ 
              ! change in ALK is identical to change in NO3

              local_bgc_mem%bgctra(j,1,ialkali)  = local_bgc_mem%bgctra(j,1,ialkali)                  &
                   &                 - (local_bgc_mem%bgctra(j,1,iano3)-oldnitrate)
             
              local_bgc_mem%bgcflux(j,knfixd)    = (- oldnitrate + local_bgc_mem%bgctra(j,1,iano3))/rn2 /dtbgc
              local_bgc_mem%bgctend(j,1,kn2b)    = local_bgc_mem%bgctend(j,1,kn2b) - (local_bgc_mem%bgctra(j,1,iano3) - oldnitrate) * (pddpo(j,1) + za(j))
              local_bgc_mem%bgctend(j,1, kaou)   = local_bgc_mem%satoxy(j,1) - local_bgc_mem%bgctra(j,1,ioxygen)
           ENDIF
        ENDIF

     ENDDO
     !$ACC END PARALLEL
     !$ACC WAIT(1)
     !$ACC END DATA
 
 
END SUBROUTINE cyano



SUBROUTINE cyadyn(local_bgc_mem, klevs, start_idx, end_idx, pddpo, za, ptho, ptiestu, l_dynamic_pi, max_klevs, lacc)
!! @brief prognostic N2 fixation, cyanobacteria

      USE mo_memory_bgc, ONLY      : pi_alpha_cya,          &
       &                            Topt_cya,T1_cya,T2_cya,bkcya_N,      &
       &                            fPAR, ro2ut, ro2ut_cya,ralk,      &
       &                            doccya_fac, rnit, riron, rcar, rn2, &
       &                            wcya, rnoi, cyamin, &
       &                            ro2ammo, bknh4_cya, bkno3_cya

      USE mo_param1_bgc, ONLY     : iano3, iphosph, igasnit, &
           &                        ioxygen, ialkali, icya,  &
           &                        isco212, idoc, kaou, &
           &                        idet, iiron, knfix, &
           &                        kpho_cya, kcyaloss, kn2b, &
           &                        kcTlim, kcLlim, kcPlim, kcFlim, &
           &                        iammo, kcyapro

      USE mo_hamocc_nml,ONLY      : cycdec, cya_growth_max, bkcya_fe, bkcya_P, &
           &                        l_N_cycle

      IMPLICIT NONE
      TYPE(t_bgc_memory), POINTER    :: local_bgc_mem

      INTEGER, INTENT(in) ::  start_idx               !< 1st REAL of model grid
      INTEGER, INTENT(in) ::  end_idx                 !< 1st REAL of model grid
      INTEGER, INTENT(in) ::  klevs(bgc_nproma)       !< 3rd (vertical) REAL of model grid.
   
      REAL(wp), INTENT(in) :: pddpo(bgc_nproma, bgc_zlevs)  !< size of scalar grid cell (3rd dimension) [m].
      REAL(wp), INTENT(in) :: ptho(bgc_nproma, bgc_zlevs)   !< potential temperature [deg C]
      REAL(wp), INTENT(in) :: za(bgc_nproma)                !< potential temperature [deg C]
      REAL(wp), INTENT(in) :: ptiestu(bgc_nproma,bgc_zlevs) !< depth of scalar grid cell [m]

      LOGICAL, INTENT(in) :: l_dynamic_pi
      INTEGER, INTENT(IN) :: max_klevs
      LOGICAL, INTENT(IN), OPTIONAL :: lacc

      !! Local variables
   
      INTEGER  :: j,k, kpke
      REAL(wp) :: oldigasnit                                                                          
      REAL(wp) :: cyapro,cyaloss
      REAL(wp) :: avanut,avcyabac                                             
      REAL(wp) :: pho,avanfe 
      REAL(wp) :: l_I, l_T
      REAL(wp) :: T_min_Topt,sgnT
      REAL(wp) :: xa_P, xa_fe, avnit,l_P,l_fe
      REAL(wp) ::  phosy_cya
      REAL(wp) :: dyn_pi_alpha_cya
      REAL(wp) :: surface_height
      LOGICAL :: lzacc

      ! for N-cycle
      REAL(wp) :: no3cya, nh4cya, hib, xa_nh4, xa_no3, no3lim
      REAL(wp) :: xn_fe, xn_p

      CALL set_acc_host_or_device(lzacc, lacc)
  
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO k = 1, max_klevs
   DO j = start_idx, end_idx

            IF( pddpo(j,k) .GT. EPSILON(0.5_wp) .and. k <= klevs(j)) THEN
                   

              avcyabac = MAX(1.e-11_wp,local_bgc_mem%bgctra(j,k,icya))                !available cyanobacteria
              avanut = MAX(0._wp,local_bgc_mem%bgctra(j,k,iphosph))                   !available phosphate    
              avanfe = MAX(0._wp,(local_bgc_mem%bgctra(j,k,iiron)/riron))             !available iron
              avnit = MAX(0._wp,local_bgc_mem%bgctra(j,k,iano3)/rnit)                 !available nitrate
                        
 
              if (l_dynamic_pi)then
                 
                   dyn_pi_alpha_cya = pi_alpha_cya + 0.05_wp* ptiestu(j,k)/(ptiestu(j,k) + 90._wp) ! pi_alpha  
                   
                   l_I = (dyn_pi_alpha_cya*fPAR*local_bgc_mem%strahl(j))*local_bgc_mem%meanswr(j,k) &     ! light limitation
                        /SQRT(cya_growth_max**2 + (dyn_pi_alpha_cya**2)*(fPAR*local_bgc_mem%strahl(j)*local_bgc_mem%meanswr(j,k))**2) 

              else
                   l_I = (pi_alpha_cya*fPAR*local_bgc_mem%strahl(j))*local_bgc_mem%swr_frac(j,k) &     ! light limitation
                        /SQRT(cya_growth_max**2 + (pi_alpha_cya**2)*(fPAR*local_bgc_mem%strahl(j)*local_bgc_mem%swr_frac(j,k))**2) 
         
              endif
              local_bgc_mem%bgctend(j,k,kcLlim) = l_I 
 
              T_min_Topt = ptho(j,k)-Topt_cya                          
              sgnT = sign(1._wp,T_min_Topt)
              IF(T_min_Topt .eq. 0._wp) sgnT = 0._wp 

              l_T = exp(-((T_min_Topt**4)/(T1_cya-T2_cya*sgnT)**4))          !temperature limitation 
              local_bgc_mem%bgctend(j,k,kcTlim) = l_T 

              xa_p = avanut                                                 
              l_P = xa_P / (bkcya_P + xa_P)                      !phosphate limitation 
              local_bgc_mem%bgctend(j,k,kcPlim) = l_P 

              xa_fe = avanfe         
              l_fe = xa_fe / (bkcya_fe + xa_fe)                  !iron limitation
              local_bgc_mem%bgctend(j,k,kcFlim) = l_fe 

               
              IF (.not. l_N_cycle) THEN

                 pho=dtb*cya_growth_max*l_I*l_T*l_P*l_fe            !growth 
             
                 cyapro = MIN(avnit-1.e-11_wp*rnoi, &
   &                       avcyabac*pho*avnit**2/(bkcya_N**2 + avnit**2))


                 phosy_cya = pho*avcyabac  
              
                 ! limitation on DIC
                 if (local_bgc_mem%bgctra(j,k,isco212).le.rcar*phosy_cya) then
                     cyapro=0._wp
                     phosy_cya=0._wp    
!                     cyapro = min(avnit-1.e-11_wp*rnoi,max(0._wp,(local_bgc_mem%bgctra(j,k,isco212)-EPSILON(1.0_wp))/(rcar*phosy_cya)*cyapro))
!                    phosy_cya=max(0._wp,(local_bgc_mem%bgctra(j,k,isco212)-EPSILON(1.0_wp)))/rcar
                 endif
                 ! ---------- nutrient uptake

              
                 oldigasnit = local_bgc_mem%bgctra(j,k,igasnit)  
                 local_bgc_mem%bgctra(j,k,igasnit) = local_bgc_mem%bgctra(j,k,igasnit) - (phosy_cya - cyapro)*rnit/rn2  ! gasnit [N2]
                 local_bgc_mem%bgctend(j,k,knfix) =  -1._wp * rn2 *(local_bgc_mem%bgctra(j,k,igasnit) - oldigasnit)/dtbgc   ! N fixation
                 surface_height = MERGE(za(j), 0._wp, k==1)
                 local_bgc_mem%bgctend(j,k,kn2b) = local_bgc_mem%bgctend(j,k,kn2b) - (phosy_cya -cyapro) * rnit*(pddpo(j,k) +surface_height)
             
                 local_bgc_mem%bgctra(j,k,iano3) = local_bgc_mem%bgctra(j,k,iano3) - cyapro*rnit                                                           

                 ! ---------- change of alkalinity 
                 !(only for production on nitrate, no change of alkalinity 
                 !for N2 fixation Wolf-Gladrow et al.(2007))
                 ! P uptake due to N2fixation changes alkalinity by 1, as 
                 ! total phosphorous is considered in alk

                 local_bgc_mem%bgctra(j,k,ialkali) = local_bgc_mem%bgctra(j,k,ialkali) + cyapro*ralk + (phosy_cya-cyapro)

                 ! ---------- oxygen production
                 ! O2 from cyano growth using NO3: cyapro [NO3] --> ro2ut 
                 ! O2 from cyano growth fixing N2: (pho - cyapro) [NO3]
                 ! -> ro2ut - 3 * rnit/rno2 = 172 - 3*16/2 = 148
                 local_bgc_mem%bgctra(j,k,ioxygen) = local_bgc_mem%bgctra(j,k,ioxygen) + (phosy_cya - cyapro)*ro2ut_cya &
               &                       + cyapro*ro2ut                                            

              ELSE
              
                 xa_nh4 =  max(0._wp,local_bgc_mem%bgctra(j,k,iammo)-1.E-9_wp)
                 xa_no3 =  max(0._wp,local_bgc_mem%bgctra(j,k,iano3)-1.E-7_wp)

              
                 pho=dtb*cya_growth_max*l_I*l_T*avcyabac 

                 xn_P = avanut/ (1._wp + pho*l_fe/(bkcya_P + xa_P))  ! max photosynthesis on P

                 xn_fe = avanfe/ (1._wp + pho*l_P/(bkcya_fe + xa_fe))  ! ma
             
             
                 phosy_cya = min(xa_P - xn_P, xa_fe - xn_fe)     ! in P units
              
                 nh4cya = MIN(xa_nh4, phosy_cya*rnit * xa_nh4/(xa_nh4 + bknh4_cya))
                 no3cya = MIN(xa_no3, (phosy_cya*rnit - nh4cya) * xa_no3/(xa_no3 + bkno3_cya))

                 cyapro = no3cya + nh4cya          ! cyapro in N units

                 local_bgc_mem%bgctra(j,k,iano3) = local_bgc_mem%bgctra(j,k,iano3) - no3cya
                 local_bgc_mem%bgctra(j,k,iammo) = local_bgc_mem%bgctra(j,k,iammo) - nh4cya   
                                                       
                 ! ---------- change of alkalinity 
                 !(only for production on nitrate, phosphate, no change of alkalinity 
                 !for N2 fixation Wolf-Gladrow et al.(2007))
                 ! P uptake due to N2fixation changes alkalinity by 1, as 
                 ! total phosphorous is considered in alk
                 local_bgc_mem%bgctra(j,k,ialkali) = local_bgc_mem%bgctra(j,k,ialkali) + phosy_cya  - nh4cya + no3cya

                 ! ---------- oxygen production
                 ! O2 from cyano growth using NO3: cyapro [NO3] --> ro2ut 
                 ! O2 from cyano growth fixing N2: (phosy_cya - cyapro) [NO3]
                 ! -> ro2ut - 3 * rnit/rno2 = 172 - 3*16/2 = 148

                 local_bgc_mem%bgctra(j,k,ioxygen)= local_bgc_mem%bgctra(j,k,ioxygen) + (phosy_cya - cyapro*rnoi)*148._wp &
                  &                       + (no3cya*ro2ut + nh4cya*ro2ammo)*rnoi 


                 local_bgc_mem%bgctend(j,k,kcyapro) = cyapro/dtbgc

                 local_bgc_mem%bgctra(j,k,igasnit) = local_bgc_mem%bgctra(j,k,igasnit) - (phosy_cya*rnit - cyapro)*0.5_wp  ! gasnit [N2]
        
                 local_bgc_mem%bgctend(j,k,knfix) =  (phosy_cya*rnit - cyapro)/dtbgc ! output budgets [N]

                 local_bgc_mem%bgctend(j,k,kn2b) = local_bgc_mem%bgctend(j,k,kn2b) - (phosy_cya*rnit - cyapro)*(pddpo(j,k) +surface_height) 
          
              ENDIF ! l_N_cycle


              local_bgc_mem%bgctra(j,k,iphosph) = local_bgc_mem%bgctra(j,k,iphosph) - phosy_cya    
              local_bgc_mem%bgctra(j,k,iiron) = local_bgc_mem%bgctra(j,k,iiron) - phosy_cya * riron  

              local_bgc_mem%bgctend(j,k,kpho_cya) =  phosy_cya/dtbgc
              local_bgc_mem%bgctend(j,k, kaou)   = local_bgc_mem%satoxy(j,k) - local_bgc_mem%bgctra(j,k,ioxygen)
              

              ! --------- change of total CO2
              local_bgc_mem%bgctra(j,k,isco212) = local_bgc_mem%bgctra(j,k,isco212) - phosy_cya*rcar


              ! --------- decay of cyanobacteria
              cyaloss = cycdec * max(0._wp,local_bgc_mem%bgctra(j,k,icya) - 2._wp*cyamin)                                

              ! --------- change of cyanobacteria
              local_bgc_mem%bgctra(j,k,icya) = local_bgc_mem%bgctra(j,k,icya) + phosy_cya - cyaloss
              local_bgc_mem%bgctend(j,k,kcyaloss) =  cyaloss/dtbgc   

              ! --------- decaying cyanobacteria are distributed to DOCCYA and detritus  
              local_bgc_mem%bgctra(j,k,idoc) = local_bgc_mem%bgctra(j,k,idoc) + doccya_fac*cyaloss     
              local_bgc_mem%bgctra(j,k,idet) = local_bgc_mem%bgctra(j,k,idet) + (1.0_wp - doccya_fac)*cyaloss        




            ENDIF ! wet cells
      ENDDO ! 
  ENDDO ! 
  !$ACC END PARALLEL
 

! -------------- buoyancy of cyanobacteria----------------------------------------

  ! implicit method:
  ! C(k,T+dt)=C(k,T) + (w*dt/ddpo(k))*(C(k-1,T+1)-C(k,T+1))
  ! -->
  ! C(k,T+dt)=(ddpo(k)*C(k,T)+w*dt*C(k-1,T+dt))/(ddpo(k)+w*dt)
  ! sedimentation=w*dt*C(ks,T+dt)
  !
 !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
 !$ACC LOOP GANG VECTOR
 DO j=start_idx,end_idx
    kpke=klevs(j)
    IF (kpke > 0)THEN
    IF(pddpo(j,kpke).GT.EPSILON(0.5_wp)) THEN
           local_bgc_mem%bgctra(j,kpke,icya)  = (local_bgc_mem%bgctra(j,kpke,icya)*pddpo(j,kpke))      &
                   &              / (pddpo(j,kpke)+wcya)
             
    ENDIF
    !$ACC LOOP SEQ
    do k=(kpke-1),2,-1  
         ! water column
        if(pddpo(j,k+1).LE.EPSILON(0.5_wp))then ! last wet cell
              local_bgc_mem%bgctra(j,k,icya)  = (local_bgc_mem%bgctra(j,k,icya)*pddpo(j,k))      &
                   &              / (pddpo(j,k)+wcya)

         else
               local_bgc_mem%bgctra(j,k,icya)    = (local_bgc_mem%bgctra(j,k  ,icya)*pddpo(j,k)    &
                   &                +  local_bgc_mem%bgctra(j,k+1,icya)*wcya)/          &
                   &                          (pddpo(j,k)+wcya)
         endif 
    
   ENDDO
    k=1
    IF((pddpo(j,k).GT.EPSILON(0.5_wp)) .and. (pddpo(j,k+1).GT.EPSILON(0.5_wp)) )then ! only if next cell also wet
         local_bgc_mem%bgctra(j,k,icya)  =  local_bgc_mem%bgctra(j,k,icya) + (wcya*local_bgc_mem%bgctra(j,k+1,icya))/(pddpo(j,k)+za(j))
    endif   
   ENDIF
 ENDDO
 !$ACC END PARALLEL
                  
 
END SUBROUTINE  cyadyn

END MODULE


