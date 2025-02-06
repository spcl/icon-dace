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

MODULE mo_settling


      USE mo_kind,    ONLY        : wp
      USE mo_bgc_memory_types, ONLY  : t_bgc_memory, t_sediment_memory
      USE mo_fortran_tools, ONLY  : set_acc_host_or_device


      IMPLICIT NONE

      PRIVATE

      PUBLIC:: settling, settling_pdm


CONTAINS

!>
! !! @file settling.f90
!! @brief compute settling of debris
      SUBROUTINE settling (local_bgc_mem, local_sediment_mem, klev,start_idx, end_idx, pddpo, za, lacc)

      USE mo_param1_bgc, ONLY     : icalc, iopal, kopex90,   &
       &                            idet, kcalex90, &
       &                            kprorca,kprcaca,   &
       &                            ksilpro,kcoex90, idust, &
       &                            kprodus, kopex1000, &
       &                            kopex2000,kcoex1000,&
       &                            kcoex2000, kcalex1000, &
       &                            kcalex2000,kwdust,kwpoc,&
       &                            kwopal,kwcal

      USE mo_memory_bgc,  ONLY    : n90depth,&
      &                             n1000depth,n2000depth

      USE mo_control_bgc, ONLY    : bgc_nproma, bgc_zlevs, dtbgc, inv_dtbgc

      IMPLICIT NONE


      ! Arguments
      TYPE(t_bgc_memory), POINTER    :: local_bgc_mem
      TYPE(t_sediment_memory), POINTER :: local_sediment_mem

      INTEGER, INTENT(in), TARGET    :: klev(bgc_nproma)       !<  vertical levels
      INTEGER, INTENT(in)            :: start_idx              !< start index for j loop (ICON cells, MPIOM lat dir)  
      INTEGER, INTENT(in)            :: end_idx                !< end index  for j loop  (ICON cells, MPIOM lat dir) 

      REAL(wp), INTENT(in), TARGET   :: pddpo(bgc_nproma,bgc_zlevs)      !< size of scalar grid cell (3rd dimension) [m]
      REAL(wp), INTENT(in), TARGET   :: za(bgc_nproma)      !< surface height
      LOGICAL, INTENT(IN), OPTIONAL :: lacc


      ! Local variables
      INTEGER,  POINTER  :: kbo(:)   !< k-index of bottom layer (2d)
      INTEGER :: k, kpke,j
      LOGICAL :: lzacc

      CALL set_acc_host_or_device(lzacc, lacc)


     ! implicit method:
     ! C(k,T+dt)=C(k,T) + (w*dt/ddpo(k))*(C(k-1,T+1)-C(k,T+1))
     ! -->
     ! C(k,T+dt)=(ddpo(k)*C(k,T)+w*dt*C(k-1,T+dt))/(ddpo(k)+w*dt)
     ! sedimentation=w*dt*C(ks,T+dt)
     !
     kbo => local_bgc_mem%kbo
     
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO j=start_idx,end_idx
        kpke=klev(j)        

        IF(kpke > 0)THEN
        IF(pddpo(j,1) > EPSILON(0.5_wp))THEN

         if(kpke>=n90depth)then
           local_bgc_mem%bgcflux(j,kcoex90)  = local_bgc_mem%bgctra(j,n90depth,idet)* &
             & local_bgc_mem%wpoc(j,n90depth)*inv_dtbgc
           local_bgc_mem%bgcflux(j,kcalex90) = local_bgc_mem%bgctra(j,n90depth,icalc)* &
             & local_bgc_mem%wcal(j,n90depth)*inv_dtbgc
           local_bgc_mem%bgcflux(j,kopex90)  = local_bgc_mem%bgctra(j,n90depth,iopal)* &
             & local_bgc_mem%wopal(j,n90depth)*inv_dtbgc
         endif

         if(kpke>=n1000depth)then
           local_bgc_mem%bgcflux(j,kcoex1000)  = local_bgc_mem%bgctra(j,n1000depth,idet)* &
             & local_bgc_mem%wpoc(j,n1000depth)*inv_dtbgc
           local_bgc_mem%bgcflux(j,kcalex1000) = local_bgc_mem%bgctra(j,n1000depth,icalc)* &
             & local_bgc_mem%wcal(j,n1000depth)*inv_dtbgc
           local_bgc_mem%bgcflux(j,kopex1000)  = local_bgc_mem%bgctra(j,n1000depth,iopal)* &
             & local_bgc_mem%wopal(j,n1000depth)*inv_dtbgc
         endif

         if(kpke>=n2000depth)then
           local_bgc_mem%bgcflux(j,kcoex2000)  = local_bgc_mem%bgctra(j,n2000depth,idet)* &
             & local_bgc_mem%wpoc(j,n2000depth)*inv_dtbgc
           local_bgc_mem%bgcflux(j,kcalex2000) = local_bgc_mem%bgctra(j,n2000depth,icalc)* & 
             & local_bgc_mem%wcal(j,n2000depth)*inv_dtbgc
           local_bgc_mem%bgcflux(j,kopex2000)  = local_bgc_mem%bgctra(j,n2000depth,iopal)* &
             & local_bgc_mem%wopal(j,n2000depth)*inv_dtbgc
         endif

                  ! -----------surface layer
         k=1 
           local_bgc_mem%bgctra(j,k,idet)  = (local_bgc_mem%bgctra(j,k,idet)*(pddpo(j,k)+za(j)))      &
                &              / (pddpo(j,k)+ za(j) +local_bgc_mem%wpoc(j,k))

           local_bgc_mem%bgctra(j,k,icalc)   = (local_bgc_mem%bgctra(j,k,icalc)*(pddpo(j,k)+za(j)))   &
                &                / (pddpo(j,k)+za(j)+local_bgc_mem%wcal(j,k))

           local_bgc_mem%bgctra(j,k,iopal)   = (local_bgc_mem%bgctra(j,k,iopal)*(pddpo(j,k)+za(j)))    &
                &                / (pddpo(j,k)+za(j)+local_bgc_mem%wopal(j,k))

           local_bgc_mem%bgctra(j,k,idust)   = (local_bgc_mem%bgctra(j,k,idust)*(pddpo(j,k)+za(j)))    &
                &                / (pddpo(j,k)+za(j)+local_bgc_mem%wdust(j,k))

           local_bgc_mem%bgctend(j,k,kwdust) = local_bgc_mem%wdust(j,k)*inv_dtbgc
           local_bgc_mem%bgctend(j,k,kwpoc)  = local_bgc_mem%wpoc(j,k)*inv_dtbgc
           local_bgc_mem%bgctend(j,k,kwopal) = local_bgc_mem%wopal(j,k)*inv_dtbgc
           local_bgc_mem%bgctend(j,k,kwcal)  = local_bgc_mem%wcal(j,k)*inv_dtbgc
          ENDIF

         !$ACC LOOP SEQ
         DO k=2,kpke
          IF(pddpo(j,k) > EPSILON(0.5_wp))THEN
          ! water column
              local_bgc_mem%bgctra(j,k,idet)  = (local_bgc_mem%bgctra(j,k  ,idet)*pddpo(j,k)    &
                   &        +  local_bgc_mem%bgctra(j,k-1,idet)*local_bgc_mem%wpoc(j,k-1))/  &
                   &                     (pddpo(j,k)+local_bgc_mem%wpoc(j,k))
              local_bgc_mem%bgctra(j,k,icalc)   = (local_bgc_mem%bgctra(j,k  ,icalc)*pddpo(j,k)   &
                   &                +  local_bgc_mem%bgctra(j,k-1,icalc)*local_bgc_mem%wcal(j,k-1))/         &
                   &                           (pddpo(j,k)+local_bgc_mem%wcal(j,k))

              local_bgc_mem%bgctra(j,k,iopal)   = (local_bgc_mem%bgctra(j,k  ,iopal)*pddpo(j,k)   &
                   &                +  local_bgc_mem%bgctra(j,k-1,iopal)*local_bgc_mem%wopal(j,k-1))/        &
                   &                           (pddpo(j,k)+local_bgc_mem%wopal(j,k))

              local_bgc_mem%bgctra(j,k,idust)   = (local_bgc_mem%bgctra(j,k  ,idust)*pddpo(j,k)   &
                   &                +  local_bgc_mem%bgctra(j,k-1,idust)*local_bgc_mem%wdust(j,k-1))/        &
                   &                           (pddpo(j,k)+local_bgc_mem%wdust(j,k))

           local_bgc_mem%bgctend(j,k,kwdust) = local_bgc_mem%wdust(j,k)*inv_dtbgc
           local_bgc_mem%bgctend(j,k,kwpoc)  = local_bgc_mem%wpoc(j,k)*inv_dtbgc
           local_bgc_mem%bgctend(j,k,kwopal) = local_bgc_mem%wopal(j,k)*inv_dtbgc
           local_bgc_mem%bgctend(j,k,kwcal)  = local_bgc_mem%wcal(j,k)*inv_dtbgc
          ENDIF
         ENDDO

        IF(pddpo(j,kbo(j)) > EPSILON(0.5_wp))THEN
 
           ! sediment fluxes at the bottom
            local_sediment_mem%prorca(j) = local_bgc_mem%bgctra(j,kbo(j),idet )*local_bgc_mem%wpoc(j,kbo(j))
            local_sediment_mem%prcaca(j) = local_bgc_mem%bgctra(j,kbo(j),icalc)*local_bgc_mem%wcal(j,kbo(j))
            local_sediment_mem%silpro(j) = local_bgc_mem%bgctra(j,kbo(j),iopal)*local_bgc_mem%wopal(j,kbo(j))
            local_sediment_mem%produs(j) = local_bgc_mem%bgctra(j,kbo(j),idust)*local_bgc_mem%wdust(j,kbo(j))
            local_bgc_mem%bgcflux(j,kprorca) = local_sediment_mem%prorca(j)*inv_dtbgc
            local_bgc_mem%bgcflux(j,kprcaca) = local_sediment_mem%prcaca(j)*inv_dtbgc
            local_bgc_mem%bgcflux(j,ksilpro) = local_sediment_mem%silpro(j)*inv_dtbgc
            local_bgc_mem%bgcflux(j,kprodus) = local_sediment_mem%produs(j)*inv_dtbgc

          ENDIF
         ENDIF
      END DO
      !$ACC END PARALLEL

      END SUBROUTINE settling 


!>
!! @file settling.f90
!! @brief compute settling of debris
      SUBROUTINE settling_pdm (local_bgc_mem, local_sediment_mem, klev,start_idx, end_idx,pddpo)
  ! using explicit P2-PDM limiter scheme (=ULTIMATE QUICKEST in Leonard 1991:
  ! The ULTIMATE conservative difference scheme applied to unsteady
  ! one-dimensional advection.
  ! see also Pietrzak 1998: The use of TVD limiters for forward-in-time
  ! upstream-biased advection schemes in ocean modeling, there described as
  ! P2-PDM)

      USE mo_param1_bgc, ONLY     : icalc, iopal, kopex90,   &
       &                            idet, kcalex90, &
       &                            kprorca,kprcaca,   &
       &                            ksilpro,kcoex90, idust, &
       &                            kprodus, kopex1000, &
       &                            kopex2000,kcoex1000,&
       &                            kcoex2000, kcalex1000, &
       &                            kcalex2000,kwdust,kwpoc,kwopal,&
       &                            kwcal

      USE mo_memory_bgc,  ONLY    : n90depth, n1000depth,n2000depth

      USE mo_control_bgc, ONLY    : bgc_nproma, bgc_zlevs, inv_dtbgc
  

      IMPLICIT NONE


      ! Arguments
      TYPE(t_bgc_memory), POINTER    :: local_bgc_mem
      TYPE(t_sediment_memory), POINTER :: local_sediment_mem

      INTEGER, INTENT(in), TARGET    :: klev(bgc_nproma)       !<  vertical levels
      INTEGER, INTENT(in)            :: start_idx              !< start index for j loop (ICON cells, MPIOM lat dir)  
      INTEGER, INTENT(in)            :: end_idx                !< end index  for j loop  (ICON cells, MPIOM lat dir) 

      REAL(wp), INTENT(in), TARGET   :: pddpo(bgc_nproma,bgc_zlevs)      !< size of scalar grid cell (3rd dimension) [m]


      ! Local variables
      INTEGER,  POINTER  :: kbo(:)   !< k-index of bottom layer (2d)
      INTEGER :: k, kpke,j

      REAL(wp) :: det_flux(bgc_zlevs+1), calc_flux(bgc_zlevs+1)
      REAL(wp) :: opal_flux(bgc_zlevs+1), dust_flux(bgc_zlevs+1)
      REAL(wp) :: cn_det, cn_calc, cn_opal, cn_dust, maxflux 
      REAL(wp) :: det_plim, calc_plim, opal_plim, dust_plim

      kbo => local_bgc_mem%kbo

 
       DO j=start_idx,end_idx
        kpke=klev(j)        

        IF(kpke > 0)THEN
        IF(pddpo(j,1) > EPSILON(0.5_wp))THEN
        
         DO k = 1,kpke + 1

           det_flux(k)  = 0._wp
           calc_flux(k) = 0._wp
           opal_flux(k) = 0._wp
           dust_flux(k) = 0._wp

           if (k < kpke+1)then
            local_bgc_mem%bgctend(j,k,kwdust) = local_bgc_mem%wdust(j,k)*inv_dtbgc
            local_bgc_mem%bgctend(j,k,kwpoc)  = local_bgc_mem%wpoc(j,k)*inv_dtbgc
            local_bgc_mem%bgctend(j,k,kwopal) = local_bgc_mem%wopal(j,k)*inv_dtbgc
            local_bgc_mem%bgctend(j,k,kwcal)  = local_bgc_mem%wcal(j,k)*inv_dtbgc
           endif
         ENDDO

        ! flux from atmosphere to surface layer euqals zero (here)
        ! hence, starting calculation of flux with k=2

         IF( kbo(j) > 1 ) THEN ! only, when water column is larger than one box 
                             ! otherwise go immediately to the calculation 
                             ! for flux to sediment

            k = 2 ! flux at the bottom of first level
            cn_det  = local_bgc_mem%wpoc(j,k-1)/pddpo(j,k-1)    ! assuming both positive
            cn_calc = local_bgc_mem%wcal(j,k-1)/pddpo(j,k-1)   ! assuming both positive
            cn_opal = local_bgc_mem%wopal(j,k-1)/pddpo(j,k-1)   ! assuming both positive
            cn_dust = local_bgc_mem%wdust(j,k-1)/pddpo(j,k-1)   ! assuming both positive

            det_plim    = 0._wp
            calc_plim   = 0._wp
            opal_plim   = 0._wp
            dust_plim   = 0._wp

            CALL plimiter(det_plim,  cn_det,  local_bgc_mem%bgctra(j,k,idet),  local_bgc_mem%bgctra(j,k-1,idet),  local_bgc_mem%bgctra(j,k-1,idet))
            CALL plimiter(calc_plim, cn_calc, local_bgc_mem%bgctra(j,k,icalc), local_bgc_mem%bgctra(j,k-1,icalc), local_bgc_mem%bgctra(j,k-1,icalc))
            CALL plimiter(opal_plim, cn_opal, local_bgc_mem%bgctra(j,k,iopal), local_bgc_mem%bgctra(j,k-1,iopal), local_bgc_mem%bgctra(j,k-1,iopal))
            CALL plimiter(dust_plim, cn_dust, local_bgc_mem%bgctra(j,k,idust), local_bgc_mem%bgctra(j,k-1,idust), local_bgc_mem%bgctra(j,k-1,idust))

            det_flux(k)  = local_bgc_mem%wpoc(j,k-1) * (local_bgc_mem%bgctra(j,k-1,idet)                    &
                              & + 0.5_wp * det_plim * (1._wp-cn_det)            &
                              & * (local_bgc_mem%bgctra(j,k,idet) - local_bgc_mem%bgctra(j,k-1,idet))) 
            maxflux      = max(local_bgc_mem%bgctra(j,k-1,idet) - EPSILON(1._wp),0._wp) * pddpo(j,k-1)  
            det_flux(k)  = min(det_flux(k), maxflux)

            calc_flux(k) = local_bgc_mem%wcal(j,k-1) * (local_bgc_mem%bgctra(j,k-1,icalc)                   &
                              & + 0.5_wp * calc_plim * (1._wp-cn_calc)          &
                              & * (local_bgc_mem%bgctra(j,k,icalc) - local_bgc_mem%bgctra(j,k-1,icalc)))
            maxflux      = max(local_bgc_mem%bgctra(j,k-1,icalc) - EPSILON(1._wp),0._wp) * pddpo(j,k-1)  
            calc_flux(k) = min(calc_flux(k), maxflux)

            opal_flux(k) = local_bgc_mem%wopal(j,k-1) * (local_bgc_mem%bgctra(j,k-1,iopal)                &
                              & + 0.5_wp * opal_plim * (1._wp-cn_opal)          &
                              & * (local_bgc_mem%bgctra(j,k,iopal) - local_bgc_mem%bgctra(j,k-1,iopal)))
            maxflux      = max(local_bgc_mem%bgctra(j,k-1,iopal) - EPSILON(1._wp),0._wp) * pddpo(j,k-1)
            opal_flux(k) = min(opal_flux(k), maxflux)
            

            dust_flux(k) = local_bgc_mem%wdust(j,k-1) * (local_bgc_mem%bgctra(j,k-1,idust)               &
                            & + 0.5_wp * dust_plim * (1._wp-cn_dust)          &
                            & * (local_bgc_mem%bgctra(j,k,idust) - local_bgc_mem%bgctra(j,k-1,idust)))
            maxflux      = max(local_bgc_mem%bgctra(j,k-1,idust) - EPSILON(1._wp),0._wp) * pddpo(j,k-1) 
            dust_flux(k) = min(dust_flux(k), maxflux)
            
        
            DO k = 3, kbo(j) ! water column

              cn_det  = local_bgc_mem%wpoc(j,k-1)/pddpo(j,k-1)    ! assuming both positive
              cn_calc = local_bgc_mem%wcal(j,k-1)/pddpo(j,k-1)   ! assuming both positive
              cn_opal = local_bgc_mem%wopal(j,k-1)/pddpo(j,k-1)   ! assuming both positive
              cn_dust = local_bgc_mem%wdust(j,k-1)/pddpo(j,k-1)   ! assuming both positive

              det_plim    = 0._wp
              calc_plim   = 0._wp
              opal_plim   = 0._wp
              dust_plim   = 0._wp

              CALL plimiter(det_plim,  cn_det,  local_bgc_mem%bgctra(j,k,idet),  local_bgc_mem%bgctra(j,k-1,idet),  local_bgc_mem%bgctra(j,k-2,idet))
              CALL plimiter(calc_plim, cn_calc, local_bgc_mem%bgctra(j,k,icalc), local_bgc_mem%bgctra(j,k-1,icalc), local_bgc_mem%bgctra(j,k-2,icalc))
              CALL plimiter(opal_plim, cn_opal, local_bgc_mem%bgctra(j,k,iopal), local_bgc_mem%bgctra(j,k-1,iopal), local_bgc_mem%bgctra(j,k-2,iopal))
              CALL plimiter(dust_plim, cn_dust, local_bgc_mem%bgctra(j,k,idust), local_bgc_mem%bgctra(j,k-1,idust), local_bgc_mem%bgctra(j,k-2,idust))

              det_flux(k)  = local_bgc_mem%wpoc(j,k-1) * (local_bgc_mem%bgctra(j,k-1,idet)                 &
                               & + 0.5_wp * det_plim  * (1._wp-cn_det)          &
                               & * (local_bgc_mem%bgctra(j,k,idet) - local_bgc_mem%bgctra(j,k-1,idet))) 
              maxflux      = max(local_bgc_mem%bgctra(j,k-1,idet) - EPSILON(1._wp),0._wp) * pddpo(j,k-1) ! 
              det_flux(k)  = min(det_flux(k), maxflux)

              calc_flux(k) = local_bgc_mem%wcal(j,k-1) * (local_bgc_mem%bgctra(j,k-1,icalc)                &
                               & + 0.5_wp * calc_plim * (1._wp-cn_calc)         &
                              & * (local_bgc_mem%bgctra(j,k,icalc) - local_bgc_mem%bgctra(j,k-1,icalc)))
              maxflux      = max(local_bgc_mem%bgctra(j,k-1,icalc) - EPSILON(1._wp),0._wp) * pddpo(j,k-1)
              calc_flux(k) = min(calc_flux(k), maxflux)

              opal_flux(k) = local_bgc_mem%wopal(j,k-1) * (local_bgc_mem%bgctra(j,k-1,iopal)                 &
                              & + 0.5_wp * opal_plim * (1._wp-cn_opal)         &
                              & * (local_bgc_mem%bgctra(j,k,iopal) - local_bgc_mem%bgctra(j,k-1,iopal)))
              maxflux      = max(local_bgc_mem%bgctra(j,k-1,iopal) - EPSILON(1._wp),0._wp) * pddpo(j,k-1) 
              opal_flux(k) = min(opal_flux(k), maxflux)

              dust_flux(k) = local_bgc_mem%wdust(j,k-1) * (local_bgc_mem%bgctra(j,k-1,idust)              &
                             & + 0.5_wp * dust_plim * (1._wp-cn_dust)          &
                             & * (local_bgc_mem%bgctra(j,k,idust) - local_bgc_mem%bgctra(j,k-1,idust)))
              maxflux      = max(local_bgc_mem%bgctra(j,k-1,idust) - EPSILON(1._wp),0._wp) * pddpo(j,k-1) 
              dust_flux(k) = min(dust_flux(k), maxflux)

            ENDDO
         ENDIF !kbo > 1?

         k = kbo(j) + 1 ! flux to sediment
         det_flux(k)  = local_bgc_mem%wpoc(j,k-1)  * local_bgc_mem%bgctra(j,k-1,idet)
         calc_flux(k) = local_bgc_mem%wcal(j,k-1)  * local_bgc_mem%bgctra(j,k-1,icalc)
         opal_flux(k) = local_bgc_mem%wopal(j,k-1) * local_bgc_mem%bgctra(j,k-1,iopal)
         dust_flux(k) = local_bgc_mem%wdust(j,k-1) * local_bgc_mem%bgctra(j,k-1,idust)

         ! calculating change of concentrations water column:
         DO k = 1,kbo(j)
             local_bgc_mem%bgctra(j,k,idet)  = local_bgc_mem%bgctra(j,k,idet)  + (det_flux(k)  - det_flux(k+1))/pddpo(j,k)
             local_bgc_mem%bgctra(j,k,icalc) = local_bgc_mem%bgctra(j,k,icalc) + (calc_flux(k) - calc_flux(k+1))/pddpo(j,k)
             local_bgc_mem%bgctra(j,k,iopal) = local_bgc_mem%bgctra(j,k,iopal) + (opal_flux(k) - opal_flux(k+1))/pddpo(j,k)
             local_bgc_mem%bgctra(j,k,idust) = local_bgc_mem%bgctra(j,k,idust) + (dust_flux(k) - dust_flux(k+1))/pddpo(j,k)
         ENDDO

         ! here n90depth+1, since the interface from box k to k+1 is w-level k+1
         ! in the previous implicit scheme, it's defined as k-flux
         IF (kbo(j) > n90depth) THEN
           local_bgc_mem%bgcflux(j,kcoex90) = det_flux(n90depth+1)*inv_dtbgc   
           local_bgc_mem%bgcflux(j,kopex90) = opal_flux(n90depth+1)*inv_dtbgc 
           local_bgc_mem%bgcflux(j,kcalex90) = calc_flux(n90depth+1)*inv_dtbgc 
         ENDIF

         ! write out fluxes at about 1000 m
         IF (kbo(j) > n1000depth) THEN
           local_bgc_mem%bgcflux(j,kcoex1000) = det_flux(n1000depth+1)*inv_dtbgc     
           local_bgc_mem%bgcflux(j,kopex1000) = opal_flux(n1000depth+1)*inv_dtbgc    
           local_bgc_mem%bgcflux(j,kcalex1000) = calc_flux(n1000depth+1)*inv_dtbgc    
         ENDIF

         ! write out fluxes at about 1950 m
         IF (kbo(j) > n2000depth) THEN
           local_bgc_mem%bgcflux(j,kcoex2000) = det_flux(n2000depth+1)*inv_dtbgc     
           local_bgc_mem%bgcflux(j,kopex2000) = opal_flux(n2000depth+1)*inv_dtbgc    
           local_bgc_mem%bgcflux(j,kcalex2000) = calc_flux(n2000depth+1)*inv_dtbgc    
         ENDIF

         IF (pddpo(j,1) > EPSILON(0.5_wp)) THEN

           local_sediment_mem%prorca(j) = det_flux(kbo(j)+1)
           local_sediment_mem%prcaca(j) = calc_flux(kbo(j)+1)
           local_sediment_mem%silpro(j) = opal_flux(kbo(j)+1)
           local_sediment_mem%produs(j) = dust_flux(kbo(j)+1)
           
           !
           !  write output
           !
           local_bgc_mem%bgcflux(j,kprorca) = local_sediment_mem%prorca(j)*inv_dtbgc     
           local_bgc_mem%bgcflux(j,kprcaca) = local_sediment_mem%prcaca(j)*inv_dtbgc     
           local_bgc_mem%bgcflux(j,ksilpro) = local_sediment_mem%silpro(j)*inv_dtbgc     
           local_bgc_mem%bgcflux(j,kprodus) = local_sediment_mem%produs(j)*inv_dtbgc     

         ENDIF 

        ENDIF ! ddpo > 0.5   
        ENDIF !    
       ENDDO ! j

 
 

      END SUBROUTINE settling_pdm 


      SUBROUTINE plimiter(p_lim, cnnr, conc_p1, conc_m, conc_m1)

      IMPLICIT NONE

      REAL(wp), INTENT(inout) :: p_lim
      REAL(wp), INTENT(in)    :: cnnr
      REAL(wp), INTENT(in)    :: conc_p1
      REAL(wp), INTENT(in)    :: conc_m
      REAL(wp), INTENT(in)    :: conc_m1

      REAL(wp)                :: ratio_p, limiter, delta_conc

      delta_conc = conc_p1 - conc_m
      IF(ABS(delta_conc) > 0._wp) THEN ! unfortunately required, as it can rarely happen that it's equal 0
         ratio_p    = (conc_m - conc_m1) / delta_conc
      ELSE 
         ! range, where limiter 2._wp/(1._wp-cnnr)) should jump in
         ! don't have to care about (1-cn), since model is close to unstable anyway, when cn=1 
         ratio_p    = 1.e37_wp ! fix me - is there a value for the largest real value defined?
      ENDIF
      

      limiter   = 0.5_wp+1._wp/6._wp*(1._wp-2._wp*cnnr)  &
                   & +  (0.5_wp-1._wp/6._wp*(1._wp-2._wp*cnnr)) * ratio_p
      p_lim  = MAX(0._wp,MIN(MIN(limiter,2._wp/(1._wp-cnnr)),2._wp * ratio_p / cnnr))
     
      END SUBROUTINE plimiter
END MODULE
