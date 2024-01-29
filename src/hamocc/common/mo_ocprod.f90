! @brief Computes plankton dynamics, OM degradation
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

MODULE mo_ocprod

  USE mo_kind, ONLY           : wp

  USE mo_bgc_memory_types, ONLY  : t_bgc_memory

  USE mo_memory_bgc, ONLY     :phytomi, grami, rnoi, riron, pi_alpha, &
       &                        fpar, bkphy, bkzoo, epsher,         &
       &                        zinges, ro2ut, remido, dyphy, spemor,     &
       &                        gammaz, gammap, ecan, rnit, ropal, bkopal,         &
       &                        rcar, relaxfe, fesoly,            &
       &                        nitdem, dremn2o,         &
       &                        n2prod, sulfate_reduction,       &
       &                        thresh_aerob, thresh_o2, prodn2o, & 
       &                        thresh_sred, dmsp, &
       &                        ralk, bkh2sox, rh2sox,&
       &                        docmin, &
       &                        bkpo4, bkfe, bkno3, bknh4, ro2ammo, no2denit, &
       &                        anamoxra, bkno2, bkrad, nitriox, nitrira, bkfe, o2thresh, rno3no2, &
       &                        rno3nh4, rnh4no2, rno2no3, alk_nrn2, rno2n2, o2den_lim
     

  USE mo_control_bgc, ONLY    : dtb, bgc_nproma, bgc_zlevs, dtbgc , inv_dtbgc
  USE mo_param1_bgc, ONLY     : icalc, iopal, ian2o, igasnit, idms, &
       &                        iphy, izoo, isilica, iphosph, &
       &                        iano3, ioxygen, idet, idoc, isco212, &
       &                        ialkali, kphosy, kremin, ih2s,    &
       &                        iiron, ksred, kdenit, kgraz, kbacfra, &
       &                        kn2b, kh2ob, kdelcar, kdelsil, kbacfrac, &
       &                        kdmsprod, kdmsuv, kdmsbac,keuexp, knlim, kflim,&
       &                        kplim, kgraton, kexudp, kexudz, kpdy,kzdy,kaou,&
       &                        kh2sprod, kh2sloss, &
       &                        iammo, iano2, &
       &                        kgppnh, kdnrn, kdnra, kanam, kammox, knitox

    USE mo_hamocc_nml, ONLY    : grazra, calmax, dremopal, drempoc, denitrification, &
       &                         l_N_cycle, no3nh4red, no3no2red , &
       &                          l_opal_q10, opal_remin_q10, opal_remin_tref, &
       &                          l_doc_q10, doc_remin_q10, doc_remin_tref, &
       &                          l_poc_q10, poc_remin_q10, poc_remin_tref
    USE mo_fortran_tools, ONLY : set_acc_host_or_device
  PUBLIC :: ocprod
           

CONTAINS


SUBROUTINE ocprod (local_bgc_mem, klev,start_idx, end_idx, ptho, pddpo, za, ptiestu, l_dynamic_pi, max_klevs, lacc)
    

  IMPLICIT NONE

  ! Arguments
  TYPE(t_bgc_memory), POINTER    :: local_bgc_mem
  
  INTEGER, INTENT(in), TARGET    :: klev(bgc_nproma)       !<  vertical levels
  INTEGER, INTENT(in)            :: start_idx              !< start index for j loop (ICON cells, MPIOM lat dir)  
  INTEGER, INTENT(in)            :: end_idx                !< end index  for j loop  (ICON cells, MPIOM lat dir) 
  REAL(wp), INTENT(in), TARGET   :: ptho(bgc_nproma,bgc_zlevs)       !<  potential temperature (degC)
  REAL(wp), INTENT(in), TARGET   :: pddpo(bgc_nproma,bgc_zlevs)      !< size of scalar grid cell (3rd dimension) [m]

  REAL(wp), INTENT(in), TARGET   :: za(bgc_nproma)      !< surface height
  REAL(wp), INTENT(in) :: ptiestu(bgc_nproma,bgc_zlevs) !< depth of scalar grid cell [m]

  LOGICAL, INTENT(in) :: l_dynamic_pi
  INTEGER, INTENT(IN) :: max_klevs
  LOGICAL, INTENT(IN), OPTIONAL :: lacc

 !  Local variables

  INTEGER                :: j,k, kpke
  REAL(wp) :: avphy,avanut,avanfe,xa,xn,ya,yn,phosy,                   &
       &      pho, phofa, temfa,                                         &
       &      avgra,grazing,avsil,graton,                                &
       &      gratpoc,grawa,bacfra,phymor,zoothresh,zoomor,excdoc,exud,  &
       &      export, delsil, delcar, remin,             &
       &      opalrem, remin2o, aou, refra,                      &
       &      o2lim, actn2o, avdet,   &
       &      maxn2o, avoxy, oxid

  REAL(wp) :: surface_height

  REAL(wp) :: dms_prod, dms_uv, dms_bac 

  LOGICAL :: lzacc

  REAL(wp) :: reminfac

 ! N-cycle variables
  REAL(wp) :: limp, limf,limn, po4lim, felim, no3lim, nh4lim, hib
  REAL(wp) :: nfrac, detn, remin_nit, rdnrn, rdnra, fdnrn, fdnra, no3rmax
  REAL(wp) :: no3a, no3c_max, detc_max, detc_act, n2ormax
  REAL(wp) :: anam, nh4n, nh4a, annpot, fammox, fnitox, newammo
  REAL(wp) :: newnitr, oxpot, nitox, ammox, remsulf, detnew, ntotlim
  REAL(wp) :: dnrn, dnra, nlim, no2a, no2c_max, n2oa, nrn2, n2oc_max 
  REAL(wp) :: anamox, avo2, n2oprod, n2on2, no2rmax

  CALL set_acc_host_or_device(lzacc, lacc)

 !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
 !$ACC LOOP GANG VECTOR COLLAPSE(2)
 DO k = 1, max_klevs
   DO j = start_idx, end_idx

      IF(pddpo(j,k) > EPSILON(0.5_wp) .and. k <= klev(j)) then

       surface_height = MERGE(za(j), 0._wp, k==1) 

       !=====PLANKTON DYNAMICS: ========================
       !===GROWTH
       !===GRAZING
       !===DYING
       !===(EXUDATION)
       !===(EXCRETION)

       avphy  = MAX(phytomi,local_bgc_mem%bgctra(j,k,iphy))                    ! 'available' phytoplankton
       avgra  = MAX(grami, local_bgc_mem%bgctra(j,k,izoo))                     ! 'available' zooplankton
       avsil  = MAX(0._wp, local_bgc_mem%bgctra(j,k,isilica))                  ! available silicate
       avanut = MAX(0._wp, MIN(local_bgc_mem%bgctra(j,k,iphosph),        &     ! available nutrients (phosphate   [kmol P /m3]
 &          rnoi*local_bgc_mem%bgctra(j,k,iano3)))                 !                     + nitrate)
       avanfe = MAX(0._wp, MIN(avanut,local_bgc_mem%bgctra(j,k,iiron)/riron))  ! available iron

       ! Nutrient limitation diagnostic
       !   Iron: flim=merge(1,0,(P<=N).and.(Fe<P))
        local_bgc_mem%bgctend(j,k,kflim)=merge(1._wp,0._wp,avanfe.eq.local_bgc_mem%bgctra(j,k,iiron)/riron)
       !Phosphate: plim=merge(1.,0.,(flim==0.).and.(P<=N))
        local_bgc_mem%bgctend(j,k,kplim)=merge(1._wp,0._wp,(local_bgc_mem%bgctend(j,k,kflim).eq.0._wp).and.(local_bgc_mem%bgctra(j,k,iphosph).le.rnoi*local_bgc_mem%bgctra(j,k,iano3)))
        ! Nitrate: nlim=merge(1.,0.,(flim==0.).and.(N<P))
        local_bgc_mem%bgctend(j,k,knlim)=merge(1._wp,0._wp,(local_bgc_mem%bgctend(j,k,kflim).eq.0._wp).and.(rnoi*local_bgc_mem%bgctra(j,k,iano3).lt.local_bgc_mem%bgctra(j,k,iphosph)))


       ! phytoplankton growth
      if(l_dynamic_pi)then
         phofa=(pi_alpha + 0.05_wp*ptiestu(j,k)/(ptiestu(j,k)+90._wp)) &
     &            *fPAR*local_bgc_mem%strahl(j)*local_bgc_mem%meanswr(j,k)
   
      else
       phofa = pi_alpha*fPAR*local_bgc_mem%strahl(j)*local_bgc_mem%meanswr(j,k)
      endif
       temfa = 0.6_wp * 1.066_wp**ptho(j,k)
       pho   = dtb*phofa*temfa/SQRT(phofa**2 + temfa**2)

       !!!! N cycle !!!!!!!!
       IF (l_N_cycle) THEN
          xa = MAX(0._wp, (local_bgc_mem%bgctra(j,k,iphosph)))
          xn = xa/(1._wp + pho*avphy/(xa + bkpo4))
          po4lim = MAX(0._wp, xa - xn)
          limp = xn/(xa + bkpo4)             

          xa = MAX(0._wp, (local_bgc_mem%bgctra(j,k,iiron)))
          xn = xa/(1._wp + pho*avphy*riron/(xa + bkfe))
          felim = MAX(0._wp,xa - xn)
          limf = xn/(xa + bkfe)              

          hib = 1._wp/(1._wp + local_bgc_mem%bgctra(j,k,iammo)/bknh4) ! sigma_inhib

          xa = MAX(0._wp, (local_bgc_mem%bgctra(j,k,iano3)))
          xn = xa/(1._wp + pho*avphy*rnit*hib/(xa + bkno3))
          no3lim = MAX(0._wp, xa - xn)
          limn = xn/(xa + bkno3)             

          xa = MAX(0._wp, (local_bgc_mem%bgctra(j,k,iammo)))
          xn = xa/(1._wp + pho*avphy*rnit/(xa + bknh4))
          nh4lim = MAX(0._wp, xa - xn )       
          ntotlim = no3lim + nh4lim          
       
          phosy = MIN(po4lim, ntotlim/rnit, felim/riron)

          nfrac = 1._wp
          if(phosy .gt. 1.E-18_wp) nfrac= nh4lim/ntotlim     ! fraction of photosynthesis on NH4
          limn = limn*(1._wp - nfrac) + nfrac*xn/( xa + bknh4) 


          IF ( limf .le. limp .and. limf .le. limn) THEN
             local_bgc_mem%bgctend(j,k,kflim) = 1._wp
             local_bgc_mem%bgctend(j,k,knlim) = 0._wp
             local_bgc_mem%bgctend(j,k,kplim) = 0._wp
          ELSE IF (limn .le. limp) THEN
             local_bgc_mem%bgctend(j,k,kplim) = 0._wp
             local_bgc_mem%bgctend(j,k,knlim) = 1._wp
             local_bgc_mem%bgctend(j,k,kflim) = 0._wp
          ELSE
             local_bgc_mem%bgctend(j,k,kplim) = 1._wp
             local_bgc_mem%bgctend(j,k,knlim) = 0._wp
             local_bgc_mem%bgctend(j,k,kflim) = 0._wp
          END IF

       ELSE

          xa    = avanfe
          xn    = xa / (1._wp + pho*avphy / (xa+bkphy) )                ! bkphy = half saturation constant
          phosy = MAX(0._wp, xa-xn)                                     ! photo synthesis 

       ENDIF
       !!!! N cycle !!!!!!!!

       phosy=MERGE(local_bgc_mem%bgctra(j,k,isco212)/rcar,phosy,local_bgc_mem%bgctra(j,k,isco212).le.rcar*phosy) ! limit phosy by available DIC
       if(pho < 0.000001_wp)phosy=0._wp 

       ! zooplankton growth, phy grazing
       ya    = avphy + phosy                                         ! new phytoplankton concentration before grazing
       yn    = (ya+grazra*avgra*phytomi/(avphy+bkzoo))            &  ! grazing
&         / (1._wp + grazra * avgra / (avphy + bkzoo))
       grazing = MAX(0._wp, ya-yn)                                   
       graton  = epsher * (1._wp - zinges) * grazing                 ! "grazing to (re-dissolved) nutrients"
       gratpoc = (1._wp - epsher) * grazing                          ! epsher=0.8 "grazing to POC"
       grawa   = epsher*zinges*grazing                               ! grazer 'wachstum(?)'

       ! plankton dying
       phymor    = dyphy * MAX(0._wp, (local_bgc_mem%bgctra(j,k,iphy) - 2._wp * phytomi))   ! phytoplankton mortality dyphy=.008*dt
       zoothresh = MAX(0._wp, (local_bgc_mem%bgctra(j,k,izoo) - 2._wp * grami))
       zoomor    = spemor*zoothresh*zoothresh                             ! zooplankton mortality
       excdoc    = gammaz*zoothresh                                       ! excretion to DOC (zooplankton)
       exud      = gammap * MAX(0._wp, (local_bgc_mem%bgctra(j,k,iphy) - 2._wp*phytomi))    ! exudation to DOC (phytoplankton)
   
       local_bgc_mem%bgctra(j,k,iphosph) = local_bgc_mem%bgctra(j,k,iphosph) &
                   &    - phosy + graton  &
                   &    + ecan*zoomor 

       !!!! N cycle !!!!!!!!
       IF (l_N_cycle) THEN
          local_bgc_mem%bgctra(j,k,iammo) =  local_bgc_mem%bgctra(j,k,iammo)                            & 
                   &                + (graton + ecan*zoomor)*rnit                   & !  all remineralization products added to NH4
                   &                - nfrac*phosy*rnit
          
          local_bgc_mem%bgctra(j,k,iano3) = local_bgc_mem%bgctra(j,k,iano3) - (1._wp - nfrac)*phosy*rnit
       ELSE
          local_bgc_mem%bgctra(j,k,iano3) = local_bgc_mem%bgctra(j,k,iano3) &
        &        + (-phosy+graton+ecan*zoomor)*rnit
       ENDIF
       !!!! N cycle !!!!!!!!
        

       export = zoomor * (1._wp - ecan) + phymor + gratpoc               ! ecan=.95, gratpoc= .2*grazing [P-units]

       local_bgc_mem%bgctra(j,k,idet) = local_bgc_mem%bgctra(j,k,idet) + export

       !=====SHELL PRODUCTIOIN ========================

       delsil = MIN(ropal * export * avsil / (avsil + bkopal), 0.5_wp * avsil)
       delcar = calmax * rcar *export * bkopal/(avsil+bkopal)  ! 'detritus linked calcium carbonate '

       local_bgc_mem%bgctra(j,k,isco212) = local_bgc_mem%bgctra(j,k,isco212) - delcar         &                ! - CACO3 production
                   &  + rcar*(  - phosy + graton + ecan*zoomor)   ! + remineralization C-units

       local_bgc_mem%bgctra(j,k,iphy) = local_bgc_mem%bgctra(j,k,iphy) + phosy - grazing       &
                   &             - phymor - exud
       !!!! N cycle !!!!!!!!
       IF (l_N_cycle) THEN

          ! LR: Why is rnit used here and not ralk ?!
          local_bgc_mem%bgctra(j,k,ialkali) = local_bgc_mem%bgctra(j,k,ialkali)   &    ! ocean with NH4 - alkalinity change
                    &           - nfrac*phosy*rnit                          &    ! alk decrease if OM from NH4 
                    &           + (1.-nfrac)*phosy*rnit                     &    ! alk increase if OM from NO3
                    &           + rnit*(graton + ecan*zoomor)               &    ! remin all to NH4
                    &           - (graton - phosy + ecan*zoomor)            &    ! PO4 changes                      
                    &           - 2._wp*delcar

          local_bgc_mem%bgctra(j,k,ioxygen) = local_bgc_mem%bgctra(j,k,ioxygen)             &
                   &            + phosy*(ro2ut*(1._wp - nfrac) + ro2ammo*nfrac)       & ! phosy from NO3 produces ro2ut, from NH4 only ro2ammo
                   &            - (graton + ecan*zoomor)*ro2ammo                      ! since all re
 
       ELSE

          local_bgc_mem%bgctra(j,k,ialkali) = local_bgc_mem%bgctra(j,k,ialkali) - 2._wp * delcar &
                               - ralk * ( - phosy + graton + ecan*zoomor)

          local_bgc_mem%bgctra(j,k,ioxygen) = local_bgc_mem%bgctra(j,k,ioxygen)                   &
                            &   + ro2ut*phosy                    &
                            &   - (graton + ecan*zoomor)*ro2ut

       ENDIF
       !!!! N cycle !!!!!!!!     

       local_bgc_mem%bgctra(j,k,izoo)  = local_bgc_mem%bgctra(j,k,izoo) + grawa - excdoc - zoomor

       local_bgc_mem%bgctra(j,k,idoc)  = local_bgc_mem%bgctra(j,k,idoc) + excdoc + exud

       local_bgc_mem%bgctra(j,k,icalc) = local_bgc_mem%bgctra(j,k,icalc) + delcar

       IF(l_opal_q10)then  ! T-dependency (Maerz et al., 2020)
         opalrem =  MAX(dremopal * opal_remin_q10**((ptho(j,k)-opal_remin_Tref)/10._wp),0.001_wp*dtb)&
                  & * local_bgc_mem%bgctra(j, k, iopal)
       ELSE
         opalrem = dremopal * 0.1_wp * (ptho(j,k) + 3.0_wp) * local_bgc_mem%bgctra(j,k,iopal)
       ENDIF


       local_bgc_mem%bgctra(j,k,isilica) = local_bgc_mem%bgctra(j,k,isilica) - delsil + opalrem

       local_bgc_mem%bgctra(j,k,iopal) = local_bgc_mem%bgctra(j,k,iopal) + delsil - opalrem

       local_bgc_mem%bgctra(j,k,iiron) = local_bgc_mem%bgctra(j,k,iiron) + (- phosy + graton + ecan*zoomor)*riron  &
                    &   - relaxfe * MAX(local_bgc_mem%bgctra(j,k,iiron) - fesoly, 0._wp)

       local_bgc_mem%bgctend(j,k,kphosy) = phosy * inv_dtbgc 
       local_bgc_mem%bgctend(j,k,kgraz) = grazing * inv_dtbgc 
       local_bgc_mem%bgctend(j,k,kgraton) = graton * inv_dtbgc 
       local_bgc_mem%bgctend(j,k,kexudp) = exud * inv_dtbgc
       local_bgc_mem%bgctend(j,k,kexudz) = excdoc * inv_dtbgc
       local_bgc_mem%bgctend(j,k,kzdy) = zoomor * inv_dtbgc
       local_bgc_mem%bgctend(j,k,kpdy) =  phymor * inv_dtbgc
       local_bgc_mem%bgctend(j,k,kdelsil) = delsil * inv_dtbgc
       local_bgc_mem%bgctend(j,k,kdelcar) = delcar * inv_dtbgc
       local_bgc_mem%bgctend(j,k,keuexp) = export * inv_dtbgc
       if (l_N_cycle) local_bgc_mem%bgctend(j,k,kgppnh) = nfrac*phosy * inv_dtbgc 

      !===== DMS ===

       dms_prod = (dmsp(5)*delsil+dmsp(4)*delcar)                        &             ! production
         &  * (1._wp + 1._wp / (ptho(j,k) + dmsp(1))**2)

       dms_bac = dmsp(3) * dtb * ABS(ptho(j,k) + 3._wp) * local_bgc_mem%bgctra(j,k,idms) & ! bacterial consumption
                   &  * (local_bgc_mem%bgctra(j,k,idms)/(dmsp(6)+local_bgc_mem%bgctra(j,k,idms)))

       dms_uv  = dmsp(2) * 4._wp * dtb * phofa * local_bgc_mem%bgctra(j, k, idms)         ! decay due to UV-radiation


       local_bgc_mem%bgctra(j,k,idms) = local_bgc_mem%bgctra(j,k,idms)                      &
                   &             + dms_prod - dms_bac - dms_uv
  
       local_bgc_mem%bgctend(j,k,kdmsprod) = dms_prod * inv_dtbgc
       local_bgc_mem%bgctend(j,k,kdmsbac)  = dms_bac * inv_dtbgc
       local_bgc_mem%bgctend(j,k,kdmsuv)   = dms_uv * inv_dtbgc

       IF (local_bgc_mem%bgctra(j,k,ioxygen) > thresh_aerob) THEN                          

           !=====AEROB REMINERALIZATION ========================
           avoxy = local_bgc_mem%bgctra(j,k,ioxygen) - thresh_aerob      ! available O2                       
       
           o2lim = local_bgc_mem%bgctra(j,k,ioxygen)/(thresh_o2 + local_bgc_mem%bgctra(j,k,ioxygen))

           ! POC decomposition
            IF(l_poc_q10)then  ! T-dependency (Maerz et al., 2020)
            xn=local_bgc_mem%bgctra(j,k,idet)/(1._wp+o2lim*drempoc*poc_remin_q10**((ptho(j,k)-poc_remin_tref)/10._wp))
           ELSE
            xn=local_bgc_mem%bgctra(j,k,idet)/(1._wp+o2lim*drempoc)
           ENDIF
           
           remin=MAX(0._wp,local_bgc_mem%bgctra(j,k,idet)-xn)

           !!!! N cycle !!!!!!!!
           IF (l_N_cycle) THEN
              remin = MIN(remin, 0.5_wp*avoxy/ro2ammo)
           ELSE
              remin = MIN(remin, avoxy/ro2ut)      
           ENDIF
           !!!! N cycle !!!!!!!!

           ! DOC decomposition
           avoxy = local_bgc_mem%bgctra(j,k,ioxygen) -remin*ro2ut - thresh_aerob      ! available O2
           IF(l_doc_q10)then  ! T-dependency (Maerz et al., 2020)
            reminfac=remido*doc_remin_q10**((ptho(j,k)-doc_remin_tref)/10._wp)
           ELSE
            reminfac=remido
           ENDIF
           xn=(local_bgc_mem%bgctra(j,k,idoc)+reminfac*docmin)/(1._wp+reminfac)
           bacfra=MAX(0._wp,local_bgc_mem%bgctra(j,k,idoc) - xn)

           !!!! N cycle !!!!!!!!
           IF (l_N_cycle) THEN
              bacfra = MERGE(-0._wp,bacfra, avoxy-bacfra*ro2ammo.lt.thresh_aerob) 
           ELSE
              bacfra = MERGE(-0._wp,bacfra, avoxy-bacfra*ro2ut.lt.thresh_aerob)
           ENDIF
           !!!! N cycle !!!!!!!!

           local_bgc_mem%bgctra(j,k,idoc)  = local_bgc_mem%bgctra(j,k,idoc) - bacfra 
       

           local_bgc_mem%bgctra(j,k,idet)  = local_bgc_mem%bgctra(j,k,idet) - remin

           local_bgc_mem%bgctra(j,k,iiron)  = local_bgc_mem%bgctra(j,k,iiron) + (bacfra+remin)*riron  

           local_bgc_mem%bgctra(j,k,iphosph) =  local_bgc_mem%bgctra(j,k,iphosph) + bacfra +remin                    

           

           local_bgc_mem%bgctra(j,k,isco212) = local_bgc_mem%bgctra(j,k,isco212)             &  
        &                + rcar*( bacfra +  remin)                       ! + remineralization C-units


           !!!! N cycle !!!!!!!!
           IF (l_N_cycle) THEN

              local_bgc_mem%bgctra(j,k,ioxygen) = local_bgc_mem%bgctra(j,k,ioxygen)             &
                         &            - bacfra*ro2ammo            & ! since all remineralization products go to NH4
                         &            - remin*ro2ammo               ! there is less oxygen demand ro2ammo=(ro2ut - 2*rnit) 

              ! LR: Why is rnit used here and not ralk ?!
              ! This is strange: why no prefactor for "PO4 update" ?
              local_bgc_mem%bgctra(j,k,ialkali) = local_bgc_mem%bgctra(j,k,ialkali)            &
                         &            -    (remin + bacfra)       & ! PO4 update
                         &            + rnit*(bacfra + remin)       ! NH

              local_bgc_mem%bgctra(j,k,iammo)= local_bgc_mem%bgctra(j,k,iammo)                 & !  all remineralization products in EU added to NH4
                         &            + (bacfra + remin)*rnit

           ELSE
       
              local_bgc_mem%bgctra(j,k,ioxygen) = local_bgc_mem%bgctra(j,k,ioxygen)       &
                        &            -(bacfra + remin)*ro2ut     

          
              local_bgc_mem%bgctra(j,k,ialkali) = local_bgc_mem%bgctra(j,k,ialkali)       &
                        &            -(bacfra +  remin)*ralk

              local_bgc_mem%bgctra(j,k,iano3)= local_bgc_mem%bgctra(j,k,iano3)            &
                        &            +(bacfra + remin)*rnit
           ENDIF
           !!!! N cycle !!!!!!!!

           !***********************************************************************
           !   There is about 1.e4 O2 on 1 N2O molecule (Broecker&Peng)
           !    refra : Tim Rixen, pers. communication
           !***********************************************************************

           aou   = local_bgc_mem%satoxy(j,k) - local_bgc_mem%bgctra(j,k,ioxygen)
           refra = 1._wp + 3._wp * (0.5_wp + SIGN(0.5_wp, aou - 1.97e-4_wp))
       

           maxn2o = (remin+bacfra)*prodn2o*ro2ut*refra*0.5_wp
           avoxy = max(0._wp, local_bgc_mem%bgctra(j,k,ioxygen)-thresh_aerob)
           actn2o = min(avoxy,maxn2o)

           local_bgc_mem%bgctra(j,k,ian2o)   = local_bgc_mem%bgctra(j,k,ian2o)                  &
                   &    + 2._wp *actn2o

           local_bgc_mem%bgctra(j,k,igasnit) = local_bgc_mem%bgctra(j,k,igasnit)                &
                   &     - 2._wp *actn2o

           local_bgc_mem%bgctra(j,k,ioxygen) = local_bgc_mem%bgctra(j,k,ioxygen)               &
                   &     - actn2o
         
           local_bgc_mem%bgctend(j,k,kremin) = remin * inv_dtbgc
           local_bgc_mem%bgctend(j,k,kbacfra) = bacfra * inv_dtbgc
           local_bgc_mem%bgctend(j,k,kdenit) = 0._wp 
       ENDIF   ! O2 >= thresh_aerob


       !!!! N-cycle !!!!!!!!
       IF (l_N_cycle) THEN
       IF (local_bgc_mem%bgctra(j,k,idet) > 1.e-15_wp) THEN
       IF (local_bgc_mem%bgctra(j,k,ioxygen) <= o2thresh) THEN  ! o2thresh = 1e-5 currently

          ! o2 limitation identical for all suboxic processes
          ! LR: this seems to be there for avoiding cases where oxygen is
          !     negative?
          o2lim = min(1._wp,1._wp - local_bgc_mem%bgctra(j,k,ioxygen)/o2thresh)

          ! convert detritus in P-units to N-units for nitrogen cycle changes  
          detn = max(0._wp, local_bgc_mem%bgctra(j,k,idet)*rnit)
          remin_nit = 0._wp

          ! denitrification on NO3
          rdnrn = o2lim*no3no2red *detn/(local_bgc_mem%bgctra(j,k,iano3) + 3.E-5_wp)
          rdnra = o2lim*no3nh4red *detn/(local_bgc_mem%bgctra(j,k,iano3) + 3.E-5_wp)

          no3rmax = rdnrn + rdnra                       ! max pot loss of NO3

          ! LR: no3rmax < 0 only occurs, if no3 is strongly negative?!
          IF (no3rmax > 0._wp) THEN
          
          fdnrn = rdnrn/no3rmax                       ! fraction each process
          fdnra = rdnra/no3rmax

          !< implicit formulation to avoid neg. nitrate concentration
          no3a = local_bgc_mem%bgctra(j,k,iano3)/(1._wp + no3rmax)   ! max change in NO3  
          no3c_max = local_bgc_mem%bgctra(j,k,iano3) - no3a         ! corresponding max NO3 loss 
          detc_max =  no3c_max*(fdnrn/rno3no2 + fdnra/rno3nh4)  ! corresponding max change in det 
          detc_act = min (0.9_wp*local_bgc_mem%bgctra(j,k,idet), detc_max)

          dnrn = fdnrn*detc_act             ! in P units
          dnra = fdnra*detc_act             ! in P untis 

          remin_nit = dnrn + dnra      ! change for DIC and PO4

          local_bgc_mem%bgctend(j,k,kdnrn) =  dnrn*rno3no2 / dtbgc    ! in N units
          local_bgc_mem%bgctend(j,k,kdnra) =  dnra*rno3nh4 / dtbgc    ! in N units
          ! LR: rate diagnostics not yet implemented
          ! bgc_o_pro(i,j,k,kradnrn) = rdnrn*86400._wp*inv_dtbgc    ! in per day
          ! bgc_o_pro(i,j,k,kradnra) = rdnra*86400._wp*inv_dtbgc

          ! change from dnrn and dnra in other tracers 
          local_bgc_mem%bgctra(j,k,idet) = local_bgc_mem%bgctra(j,k,idet) - dnrn - dnra           ! change in detritus                        

          local_bgc_mem%bgctra(j,k,iano3) = local_bgc_mem%bgctra(j,k,iano3)    &    ! change in nitrate 
                        &     - rno3no2*dnrn       &    ! from DNRN
                        &     - rno3nh4*dnra            ! from DNRA

          local_bgc_mem%bgctra(j,k,iammo) = local_bgc_mem%bgctra(j,k,iammo)    &    ! change in ammonium 
                      &       + rnit*dnrn          &    ! from DNRN
                      &       + 86._wp*dnra             ! from DNRA

          local_bgc_mem%bgctra(j,k,iano2) = local_bgc_mem%bgctra(j,k,iano2)    &    ! change in nitrite 
                      &       + rno3no2*dnrn            ! from DNRN

          ! LR: Why rnit?!
          !     Where does the 86 come from and why is it used for changes in 
          !     ammo and alkali, even though units should be different?
          local_bgc_mem%bgctra(j,k,ialkali) = local_bgc_mem%bgctra(j,k,ialkali)  &   ! change from DNRN and DNRA
                                &          + rnit*dnrn             &   ! from DNRN
                                &          + 86._wp*dnra           &   ! from DNRA
                                &          - remin_nit                 ! PO4 update

          local_bgc_mem%bgctend(j,k,kn2b) = local_bgc_mem%bgctend(j,k,kn2b) - 70._wp*dnra*(pddpo(j,k) + surface_height)


          local_bgc_mem%bgctra(j,k,iphosph) = local_bgc_mem%bgctra(j,k,iphosph) + remin_nit
          local_bgc_mem%bgctra(j,k,isco212) = local_bgc_mem%bgctra(j,k,isco212) + rcar*remin_nit
          local_bgc_mem%bgctra(j,k,iiron)   = local_bgc_mem%bgctra(j,k,iiron)  + riron*remin_nit 

          ENDIF ! no3rmax > 0.0

          ! allow complete denitrification only at very low O2
          ! denitrification  on NO2
          if (local_bgc_mem%bgctra(j,k,ioxygen) < o2den_lim) then
             detn = max(0._wp,local_bgc_mem%bgctra(j,k,idet)*rnit)
             o2lim = 1._wp - local_bgc_mem%bgctra(j,k,ioxygen)/o2thresh
             no2rmax = o2lim*no2denit*detn/(local_bgc_mem%bgctra(j,k,iano2) + 0.1E-6_wp)

             ! implicit formulation to avoid neg. nitrite concentration
             no2a = local_bgc_mem%bgctra(j,k,iano2)/(1._wp + no2rmax) 
             no2c_max = max(0._wp,local_bgc_mem%bgctra(j,k,iano2) - no2a)         ! corresponding max NO2 loss 
             detc_max=  no2c_max/rno2n2                    ! corresponding max change in det, rno2n2 conversion to P 
             nrn2 = min (local_bgc_mem%bgctra(j,k,idet),detc_max)   ! in P units
             nrn2 = max(0._wp,nrn2)

             ! change from nrn2 in other tracers 
             local_bgc_mem%bgctra(j,k,idet)  = local_bgc_mem%bgctra(j,k,idet) - nrn2          ! change in detritus
             local_bgc_mem%bgctra(j,k,iano2) = local_bgc_mem%bgctra(j,k,iano2) - rno2n2*nrn2 ! change in nitrite      
             local_bgc_mem%bgctra(j,k,iammo) = local_bgc_mem%bgctra(j,k,iammo) + rnit*nrn2   ! change in ammonium 

             local_bgc_mem%bgctra(j,k,igasnit) = local_bgc_mem%bgctra(j,k,igasnit) & 
                           &       + rno2n2*nrn2/2._wp           ! from nitrite reduction      

             local_bgc_mem%bgctra(j,k,ialkali) = local_bgc_mem%bgctra(j,k,ialkali)  &  
                           &       + alk_nrn2*nrn2      &
                           &       - nrn2                   

             local_bgc_mem%bgctend(j,k,kn2b) = local_bgc_mem%bgctend(j,k,kn2b)       &
                           &       + (alk_nrn2 - rnit)*nrn2*(pddpo(j,k) + surface_height)


             local_bgc_mem%bgctend(j,k,kh2ob) = local_bgc_mem%bgctend(j,k,kh2ob) + rno2n2*nrn2*0.25_wp         &
                           &      * (pddpo(j,k) + surface_height)  ! produces losses 560/3 o2

             local_bgc_mem%bgctra(j,k,iphosph) = local_bgc_mem%bgctra(j,k,iphosph) + nrn2
             local_bgc_mem%bgctra(j,k,isco212) = local_bgc_mem%bgctra(j,k,isco212) + rcar*nrn2
             local_bgc_mem%bgctra(j,k,iiron)   = local_bgc_mem%bgctra(j,k,iiron)  + riron*nrn2

             ! LR: not yet implemented
             ! bgcprod(i,j,k,kremnit) = bgcprod(i,j,k,kremnit) + nrn2*inv_dtbgc ! in P -units
             !bgc_o_pro(i,j,k,kradeni) = no2rmax*86400._wp*inv_dtbgc

             local_bgc_mem%bgctend(j,k,kdenit)  = rno2n2*nrn2 / dtbgc  ! in N units

             ! denitrification using laughing gas
             ! with correct stoichometric
             detn = max(0._wp,local_bgc_mem%bgctra(j,k,idet)*rnit)

             n2ormax = o2lim*dremn2o*detn/(local_bgc_mem%bgctra(j,k,ian2o)+5.E-9_wp)

             n2oa = local_bgc_mem%bgctra(j,k,ian2o)/(1._wp + n2ormax) 
             n2oc_max = max(0._wp, local_bgc_mem%bgctra(j,k,ian2o) - n2oa)         ! corresponding max N2O loss 
             detc_max =  n2oc_max/280._wp                    ! corresponding max change in det, rno2n2 conversion to P 
             n2on2 = min (local_bgc_mem%bgctra(j,k,idet), detc_max)   ! in P units
             n2on2 = max(0._wp, n2on2)   

             ! change from nrn2 in other tracers excl. DIC and PO4, done later
             local_bgc_mem%bgctra(j,k,idet)  = local_bgc_mem%bgctra(j,k,idet) - n2on2              ! change in detritus
             local_bgc_mem%bgctra(j,k,ian2o) = local_bgc_mem%bgctra(j,k,ian2o) - 280._wp*n2on2    ! change in nitrous oxide
             local_bgc_mem%bgctra(j,k,iammo) = local_bgc_mem%bgctra(j,k,iammo) + rnit*n2on2       ! change in ammonium 

             local_bgc_mem%bgctra(j,k,igasnit) = local_bgc_mem%bgctra(j,k,igasnit) & 
                           &       + 280._wp*n2on2                    ! dinitrogen production

             ! LR: again, why rnit?
             local_bgc_mem%bgctra(j,k,ialkali) = local_bgc_mem%bgctra(j,k,ialkali)  &  
                           &       + rnit*n2on2                &
                           &       - n2on2                     

             local_bgc_mem%bgctra(j,k,iphosph) = local_bgc_mem%bgctra(j,k,iphosph) + n2on2
             local_bgc_mem%bgctra(j,k,isco212) = local_bgc_mem%bgctra(j,k,isco212) + rcar*n2on2
             local_bgc_mem%bgctra(j,k,iiron)   = local_bgc_mem%bgctra(j,k,iiron)  + riron*n2on2

             ! LR: not yet implemented
             !bgcprod(i,j,k,kremn2o) = n2on2*inv_dtbgc ! in P -units
          else
             local_bgc_mem%bgctend(j,k,kdenit) = 0._wp
          end if ! O2 < o2denlim   
       end if ! O2 < o2thresh   
       ENDIF ! det > 1.e-15

       ELSE ! no extendend N-cycle

       IF (local_bgc_mem%bgctra(j,k,ioxygen) <= o2den_lim) THEN                          
            !=====DENITRIFICATION ========================

           avdet = MAX(1.e-15_wp,local_bgc_mem%bgctra(j,k,idet))

           ! NO3 reduction
           o2lim = 1._wp - local_bgc_mem%bgctra(j,k,ioxygen)/thresh_o2

           remin  = denitrification*o2lim  &
                   & * MIN(avdet,0.5_wp*local_bgc_mem%bgctra(j,k,iano3)/nitdem)

           local_bgc_mem%bgctra(j,k,idet)      = local_bgc_mem%bgctra(j,k,idet)      &
                   &                -  remin
           ! N2O reduction

           remin2o = dremn2o * MIN(avdet,                    &  ! remineralization using N2O
                   &   0.003_wp * local_bgc_mem%bgctra(j,k,ian2o) / (2._wp*ro2ut))


           local_bgc_mem%bgctra(j,k,idet)      = local_bgc_mem%bgctra(j,k,idet)      &
                   &                -  remin2o
     
           local_bgc_mem%bgctra(j,k,iphosph) = local_bgc_mem%bgctra(j,k,iphosph) + remin + remin2o 

           local_bgc_mem%bgctra(j,k,isco212)  = local_bgc_mem%bgctra(j,k,isco212)  &
                   &                 + rcar*(remin + remin2o)



           local_bgc_mem%bgctra(j,k,iano3)   = local_bgc_mem%bgctra(j,k,iano3)     &
                   &          - nitdem*remin + rnit*(remin + remin2o)

           local_bgc_mem%bgctra(j,k,igasnit) = local_bgc_mem%bgctra(j,k,igasnit)   &
                   &                + n2prod*remin + 2._wp*ro2ut*remin2o
          

           local_bgc_mem%bgctra(j,k,iiron)   = local_bgc_mem%bgctra(j,k,iiron)     &
                   &          + riron*(remin+remin2o)

           local_bgc_mem%bgctra(j,k,ian2o)   = local_bgc_mem%bgctra(j,k,ian2o)     &
                   &         - 2._wp * ro2ut * remin2o

           !alkalinity is increased during denitrification due to consumption of H+ (see Wolf-Gladrow etal,2007)
           local_bgc_mem%bgctra(j,k,ialkali) = local_bgc_mem%bgctra(j,k,ialkali)   &
                   &          + nitdem*remin - ralk*(remin+remin2o)

          local_bgc_mem%bgctend(j,k,kn2b) = local_bgc_mem%bgctend(j,k,kn2b) + nitdem * remin * (pddpo(j,k) + surface_height) 
          !denitrification produces water (H2O), the corresponding O2 uptake is budgeted in h2obudget
          local_bgc_mem%bgctend(j,k,kh2ob) = local_bgc_mem%bgctend(j,k,kh2ob) + 0.5_wp * n2prod * remin * (pddpo(j,k) + surface_height) 

           local_bgc_mem%bgctend(j,k,kdenit) = remin * inv_dtbgc
           local_bgc_mem%bgctend(j,k,kbacfrac) = 0._wp 

       ENDIF ! oxygen < thresh_aerob
       ENDIF !!!! N-cycle !!!!!!!!


       !!!! N-cycle !!!!!!!!
       IF (l_N_cycle) THEN

          !!!! Anammox
          IF(local_bgc_mem%bgctra(j,k,ioxygen) < o2thresh) THEN       ! in suboxic water

             ! anamox on NO2 and NH4
             nh4a = max(0._wp,local_bgc_mem%bgctra(j,k,iammo) - 1.E-15_wp)
             no2a = max(0._wp,local_bgc_mem%bgctra(j,k,iano2) - 1.E-15_wp)

             o2lim = 1._wp -max(0._wp,local_bgc_mem%bgctra(j,k,ioxygen)/o2thresh)
             anam = o2lim*anamoxra*no2a/(no2a+bkno2)   

             nh4n = nh4a/(1._wp + anam)
             annpot = nh4a - nh4n            ! potential change in NH4 due to anamox
             anamox = min(annpot, no2a/1.3_wp)

             local_bgc_mem%bgctra(j,k,iammo) = local_bgc_mem%bgctra(j,k,iammo) - anamox           ! from anamox               
             local_bgc_mem%bgctra(j,k,iano2) = local_bgc_mem%bgctra(j,k,iano2) - 1.3_wp*anamox    ! from anamox               
             local_bgc_mem%bgctra(j,k,igasnit) = local_bgc_mem%bgctra(j,k,igasnit) + anamox       ! from anamox, gasnit in N2 
             local_bgc_mem%bgctra(j,k,iano3) = local_bgc_mem%bgctra(j,k,iano3) + 0.3_wp*anamox    ! change in nitrite 
             local_bgc_mem%bgctra(j,k,ialkali) = local_bgc_mem%bgctra(j,k,ialkali) - 0.3_wp*anamox

             ! loss of Os from NO2 - gain NO3 - 0.5 from NH4 =1.3 - 0.3*1.5- 0.5 = +0.35
             local_bgc_mem%bgctend(j,k,kh2ob) = local_bgc_mem%bgctend(j,k,kh2ob)          &
                           &       + anamox*0.35_wp * (pddpo(j,k) + surface_height)          

             local_bgc_mem%bgctend(j,k,kn2b) = local_bgc_mem%bgctend(j,k,kn2b)  + 1.7_wp*anamox  &
                          &      * (pddpo(j,k) + surface_height)  ! alk is changed for NO3, but for change 
                                                               ! in NH4 we need n2bugdet change 

             local_bgc_mem%bgctend(j,k,kanam)  = 2._wp*anamox / dtbgc !loss to N2 from NH4 and NO2, in kmolN /m3 s 

             ! LR: not yet implemented  
             ! bgc_o_pro(i,j,k,kraanam) = anam*86400._wp*inv_dtbgc
          ELSE
             local_bgc_mem%bgctend(j,k,kanam)  = 0._wp
             ! bgc_o_pro(i,j,k,kraanam)= 0._wp
          ENDIF ! O2 < 0.5


          !!!! Oxidation of NH4 and NO2
          fammox = nitrira*bkrad/(bkrad + local_bgc_mem%swr_frac(j,k)*local_bgc_mem%strahl(j)) ! local rate of ammox per time step
          fnitox = nitriox*bkrad/(bkrad + local_bgc_mem%swr_frac(j,k)*local_bgc_mem%strahl(j)) ! local rate of nitrite ox

          ! LR: rates not yet implemented
          !bgc_o_pro(i,j,k,kraammox) = fammox*86400._wp*inv_dtbgc
          !bgc_o_pro(i,j,k,kranitox) = fnitox*86400._wp*inv_dtbgc

          nh4a  =  local_bgc_mem%bgctra(j,k,iammo)                         
          newammo = nh4a/(1._wp + fammox)

          no2a =  max(0._wp, local_bgc_mem%bgctra(j,k,iano2))  
          newnitr = no2a/ (1._wp+fnitox)                        ! change of nitrite

          oxpot = max(0._wp, (local_bgc_mem%bgctra(j,k,ioxygen) - 0.5E-6_wp)/rno2no3)     ! max change in o2
          nitox = min(no2a - newnitr, oxpot)

          oxpot =  max(0._wp, (local_bgc_mem%bgctra(j,k,ioxygen) - 0.5E-6_wp - rno2no3*nitox)/rnh4no2)
          ammox = min(nh4a - newammo, oxpot)

          local_bgc_mem%bgctra(j,k,iammo) = nh4a - ammox
          local_bgc_mem%bgctra(j,k,iano2) = local_bgc_mem%bgctra(j,k,iano2) + ammox - nitox ! change of nitrite
          local_bgc_mem%bgctra(j,k,iano3) = local_bgc_mem%bgctra(j,k,iano3) + nitox
          local_bgc_mem%bgctra(j,k,ioxygen) = local_bgc_mem%bgctra(j,k,ioxygen)               &
                   &            - rno2no3*nitox - rnh4no2*ammox        ! O2 will be used during nitrification 

          local_bgc_mem%bgctra(j,k,ialkali) = local_bgc_mem%bgctra(j,k,ialkali)               &    ! ocean with NH4 - alkalinity change
                    &           - 2._wp*ammox              ! according to Wolf-Gladrow etal Mar. Chem. (2007)

          local_bgc_mem%bgctend(j,k,kammox) =  ammox/dtbgc                    
          local_bgc_mem%bgctend(j,k,knitox) =  nitox/dtbgc                    

          ! N2O production relies on nitrification on ammonimum
          aou   = local_bgc_mem%satoxy(j,k)-local_bgc_mem%bgctra(j,k,ioxygen)
          refra = 1._wp + 3._wp * (0.5_wp + SIGN(0.5_wp, aou - 1.97e-4_wp))
          maxn2o = rnh4no2*ammox * prodn2o * refra
          avo2 = max(0._wp,local_bgc_mem%bgctra(j,k,ioxygen)- 0.5e-6_wp)
          n2oprod =  min(avo2,maxn2o)
       
          local_bgc_mem%bgctra(j,k,ian2o)   = local_bgc_mem%bgctra(j,k,ian2o) + n2oprod
          local_bgc_mem%bgctra(j,k,igasnit) = local_bgc_mem%bgctra(j,k,igasnit) - n2oprod
          ! divide by 2 since it consumes 1 mol oxygen which is given in mol O2 
          local_bgc_mem%bgctra(j,k,ioxygen) = local_bgc_mem%bgctra(j,k,ioxygen) - n2oprod *0.5_wp

       ENDIF
       !!!! N-cycle !!!!!!!!


       local_bgc_mem%bgctend(j,k,ksred) = 0._wp 

       !=====SULFATE REDUCTION ========================
       IF (local_bgc_mem%bgctra(j,k, ioxygen) < thresh_sred) THEN

             IF (l_N_cycle) THEN

                o2lim = 1._wp - local_bgc_mem%bgctra(j,k,ioxygen)/o2thresh
                nlim = max(0._wp,1._wp - (local_bgc_mem%bgctra(j,k,iano3) + local_bgc_mem%bgctra(j,k,iano2))/10.E-6_wp)
        
                avdet = max(0._wp,local_bgc_mem%bgctra(j,k,idet))
                remsulf = o2lim*nlim*sulfate_reduction*dtb 
                detnew = avdet/(1._wp+remsulf)
                remin = avdet - detnew

                ! if we assume full oxdidation of ammo with sulfur according to Thamdrup (anammox tpye)
                ! we get N2 instead of ammonium and a smaller alkalinity change
                ! LR: Check use of rnit/ralk !!
                local_bgc_mem%bgctra(j,k,ialkali) = local_bgc_mem%bgctra(j,k,ialkali)  + 2._wp*rnit * remin + remin
                local_bgc_mem%bgctra(j,k,igasnit)   = local_bgc_mem%bgctra(j,k,igasnit)  + 0.5_wp *rnit  *remin

                local_bgc_mem%bgctra(j,k,iphosph) = local_bgc_mem%bgctra(j,k,iphosph) + remin
                local_bgc_mem%bgctra(j,k,iiron)   = local_bgc_mem%bgctra(j,k,iiron) + riron * remin
                local_bgc_mem%bgctra(j,k,ih2s)    = local_bgc_mem%bgctra(j,k,ih2s) + 2._wp*rnit * remin + remin!

                local_bgc_mem%bgctend(j,k,kn2b) = local_bgc_mem%bgctend(j,k,kn2b)  + 3._wp*rnit *remin*(pddpo(j,k) + surface_height)
                !sulphate reduction indirectly effects O2 bugdet, which is budgeted in h2obudget
                local_bgc_mem%bgctend(j,k,kh2ob) = local_bgc_mem%bgctend(j,k,kh2ob) - (ro2ammo + 0.5_wp*rnit)*remin*(pddpo(j,k) + surface_height)
       
             ELSE
                o2lim = 1._wp - local_bgc_mem%bgctra(j,k,ioxygen)/thresh_o2

                xa = max(0._wp,local_bgc_mem%bgctra(j,k,idet))
                xn = xa / (1._wp + sulfate_reduction*o2lim)
                remin = MAX(0._wp, xa-xn)

 
                local_bgc_mem%bgctra(j,k,idet)    = local_bgc_mem%bgctra(j,k,idet)    -        remin
                local_bgc_mem%bgctra(j,k,ialkali) = local_bgc_mem%bgctra(j,k,ialkali) + ralk * remin 
                local_bgc_mem%bgctra(j,k,ih2s)    = local_bgc_mem%bgctra(j,k,ih2s) + ralk * remin 
                local_bgc_mem%bgctra(j,k,isco212) = local_bgc_mem%bgctra(j,k,isco212) + rcar * remin
                local_bgc_mem%bgctra(j,k,iphosph) = local_bgc_mem%bgctra(j,k,iphosph) +        remin

                local_bgc_mem%bgctra(j,k,iano3)   = local_bgc_mem%bgctra(j,k,iano3)   + rnit  * remin
                local_bgc_mem%bgctra(j,k,iiron)   = local_bgc_mem%bgctra(j,k,iiron)   + riron * remin
     
                local_bgc_mem%bgctend(j,k,kn2b) = local_bgc_mem%bgctend(j,k,kn2b) + 2._wp * ralk * remin * (pddpo(j,k) + surface_height) 
                local_bgc_mem%bgctend(j,k,kh2ob) = local_bgc_mem%bgctend(j,k,kh2ob) - ro2ut * remin * (pddpo(j,k) + surface_height)
             ENDIF
             !!!! N-cycle !!!!!!!!

             local_bgc_mem%bgctend(j,k,ksred) = remin* inv_dtbgc
             local_bgc_mem%bgctend(j,k,kh2sprod) =  ralk * remin * inv_dtbgc
             local_bgc_mem%bgctend(j,k,kh2sloss) =  0._wp

       else
             ! HS oxidation 
               o2lim = local_bgc_mem%bgctra(j,k,ioxygen)/(local_bgc_mem%bgctra(j,k,ioxygen) + bkh2sox)
               xa = max(0._wp,local_bgc_mem%bgctra(j,k,ih2s))
               xn = xa / ( 1._wp + rh2sox*o2lim)
               oxid = max(0._wp, xa-xn)

               local_bgc_mem%bgctra(j,k,ih2s) = local_bgc_mem%bgctra(j,k,ih2s) - oxid

               local_bgc_mem%bgctra(j,k,ialkali) = local_bgc_mem%bgctra(j,k,ialkali) - 2._wp  * oxid
               local_bgc_mem%bgctend(j,k,kn2b) = local_bgc_mem%bgctend(j,k,kn2b) - 2._wp * oxid * (pddpo(j,k) + surface_height) 
   
               local_bgc_mem%bgctend(j,k,kh2sprod) =  0._wp
               local_bgc_mem%bgctend(j,k,kh2sloss) =  oxid * inv_dtbgc

       ENDIF ! O2 < thresh_sred

       local_bgc_mem%bgctend(j,k,kaou)   = local_bgc_mem%satoxy(j,k) - local_bgc_mem%bgctra(j,k,ioxygen)
  
      ENDIF ! wet cells
     ENDDO ! k=1,kpke
 ENDDO ! j=start_idx,end_idx 
 !$ACC END PARALLEL
 
 

END SUBROUTINE ocprod
END MODULE
