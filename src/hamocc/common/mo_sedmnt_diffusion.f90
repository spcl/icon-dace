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

! @brief sediment pore water diffusion
! routines: dipowa, powadi

MODULE mo_sedmnt_diffusion

 USE mo_kind, ONLY           : wp
 USE mo_hamocc_nml, ONLY     : ks,porwat, l_N_cycle
 USE mo_bgc_memory_types, ONLY  : t_bgc_memory, t_sediment_memory
 USE mo_sedmnt, ONLY : zcoefsu,zcoeflo
 USE mo_fortran_tools, ONLY  : set_acc_host_or_device
 
 IMPLICIT NONE

 PRIVATE

 PUBLIC:: dipowa, powadi


CONTAINS

SUBROUTINE DIPOWA (local_bgc_mem, local_sediment_mem, start_idx, end_idx, lacc)
!! @brief diffusion of pore water
!!
!! vertical diffusion of sediment pore water tracers
!! calculate vertical diffusion of sediment pore water properties
!! and diffusive flux through the ocean/sediment interface.
!! integration.
!!
!! implicit formulation;
!! constant diffusion coefficient : 1.e-9 set in mo_sedment.
!! diffusion coefficient : zcoefsu/zcoeflo for upper/lower
!! sediment layer boundary.

  USE mo_sedmnt, ONLY         : sedict, seddzi, seddw, &
       &                        porwah

  USE mo_control_bgc, ONLY    : dtbgc

  USE mo_param1_bgc, ONLY     : npowtra, ipowaox,  ioxygen,  &
  &                             ipowno3, ipowasi, iphosph, iano3, &
  &                             isilica, ipowafe, iiron, ialkali, &  
  &                             isco212, igasnit, ipowaph, ipowaal, &
  &                             ipown2, ipowaic, ipowh2s, ih2s, &
  &                             iammo, iano2, ipownh4, ipowno2 
  USE mo_ocean_nml, ONLY      : lsediment_only

  IMPLICIT NONE

  !! Arguments
  TYPE(t_bgc_memory), POINTER :: local_bgc_mem
  TYPE(t_sediment_memory), POINTER :: local_sediment_mem

  INTEGER, INTENT(in)  :: start_idx    !< start index for j loop (ICON cells, MPIOM lat dir)          
  INTEGER, INTENT(in)  :: end_idx      !< end index  for j loop  (ICON cells, MPIOM lat dir) 
  LOGICAL, INTENT(IN), OPTIONAL :: lacc

  !! Local variables
  INTEGER,  POINTER  :: kbo(:)   !< k-index of bottom layer (2d)
  INTEGER :: j,k,l,iv
  INTEGER :: iv_oc                         !< index of local_bgc_mem%bgctra in local_sediment_mem%powtra loop

  REAL(wp) :: sedb1(0:ks,npowtra)          !< 
  REAL(wp) :: tredsy(0:ks,3)               !< redsy for 'reduced system'

  REAL(wp) :: aprior                       !< start value of oceanic tracer in bottom layer
  LOGICAL :: lzacc

  CALL set_acc_host_or_device(lzacc, lacc)

  !
  ! --------------------------------------------------------------------
  !
  kbo => local_bgc_mem%kbo  

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG VECTOR PRIVATE(sedb1, tredsy)
  DO j = start_idx, end_idx
        
    if(local_bgc_mem%bolay(j) > EPSILON(0.5_wp))then
        k = 0

        tredsy(k,1) = zcoefsu(k)
        tredsy(k,3) = zcoeflo(k)
        ! dz(kbo) - diff upper - diff lower
        tredsy(k,2) =  local_bgc_mem%bolay(j) - tredsy(k,1) - tredsy(k,3)

        !$ACC LOOP SEQ
        DO iv = 1, npowtra      ! loop over pore water tracers

          iv_oc = iv

          if(iv == ipowaox) iv_oc = ioxygen
          if(iv == ipowno3) iv_oc = iano3
          if(iv == ipowasi) iv_oc = isilica
          if(iv == ipowafe) iv_oc = iiron
          if(iv == ipowaal) iv_oc = ialkali
          if(iv == ipowaph) iv_oc = iphosph
          if(iv == ipown2) iv_oc = igasnit
          if(iv == ipowaic) iv_oc = isco212
          if(iv == ipowh2s) iv_oc = ih2s

          if (l_N_cycle) then
            if(iv == ipownh4) iv_oc = iammo
            if(iv == ipowno2) iv_oc = iano2
          endif

          sedb1( k, iv) = 0._wp
           ! tracer_concentration(kbo) * dz(kbo)
          sedb1(k,iv) = local_bgc_mem%bgctra(j,kbo(j),iv_oc) * local_bgc_mem%bolay(j) 
   
           
        END DO

        !$ACC LOOP SEQ
        DO k = 1, ks
           tredsy(k,1) = zcoefsu(k)
           tredsy(k,3) = zcoeflo(k)
           tredsy(k,2) = seddw(k) * porwat(k) - tredsy(k,1) - tredsy(k,3)
       END DO

       !$ACC LOOP SEQ
       DO iv= 1, npowtra
        !$ACC LOOP SEQ
        DO k = 1, ks
              ! tracer_concentration(k[1:ks]) * porewater fraction(k) * dz(k)
              sedb1(k,iv) = local_sediment_mem%powtra(j,k,iv) * porwat(k) * seddw(k)
        END DO
       END DO

       !$ACC LOOP SEQ
       DO k = 1, ks
              ! this overwrites tredsy(k=0) for k=1
              tredsy(k-1,1) = tredsy(k,1) / tredsy(k-1,2)
              !                 diff upper    / conc (k-1)
              tredsy(k,2)   = tredsy(k,2)                      &
                   &          - tredsy(k-1,3) * tredsy(k,1) / tredsy(k-1,2)
              !   concentration -diff lower     * diff upper    / conc(k-1)
       END DO

     ! diffusion from above
     !$ACC LOOP SEQ
     DO iv = 1, npowtra
        !$ACC LOOP SEQ
        DO k = 1, ks
              sedb1(k,iv) = sedb1(k,iv)                        &
                   &        - tredsy(k-1,1) * sedb1(k-1,iv)
        END DO
     END DO

     ! sediment bottom layer
     k = ks
     !$ACC LOOP SEQ
     DO iv = 1, npowtra
        local_sediment_mem%powtra(j,k,iv) = sedb1(k,iv) / tredsy(k,2)
     END DO


     ! sediment column
     !$ACC LOOP SEQ
     DO iv = 1, npowtra
        !$ACC LOOP SEQ
        DO k = 1, ks-1
           l = ks-k
           local_sediment_mem%powtra(j,l,iv) = ( sedb1(l,iv)            &
                      &             - tredsy(l,3) * local_sediment_mem%powtra(j,l+1,iv) )    &
                      &             / tredsy(l,2)
        END DO
     END DO

     ! sediment ocean interface
     !$ACC LOOP SEQ
     DO iv = 1, npowtra
        ! 
        ! check mo_param1_bgc.f90 for consistency
        iv_oc = iv
        if(iv == ipowaox) iv_oc = ioxygen
        if(iv == ipowno3) iv_oc = iano3
        if(iv == ipowasi) iv_oc = isilica
        if(iv == ipowafe) iv_oc = iiron
        if(iv == ipowaal) iv_oc = ialkali
        if(iv == ipowaph) iv_oc = iphosph
        if(iv == ipown2) iv_oc = igasnit
        if(iv == ipowaic) iv_oc = isco212
        if(iv == ipowh2s) iv_oc = ih2s

        if (l_N_cycle) then
          if(iv == ipownh4) iv_oc = iammo
          if(iv == ipowno2) iv_oc = iano2
        endif

           l = 0

         aprior = local_bgc_mem%bgctra(j,kbo(j),iv_oc)
         local_bgc_mem%bgctra(j,kbo(j),iv_oc) =                                  &
              &         ( sedb1(l,iv) - tredsy(l,3) * local_sediment_mem%powtra(j,l+1,iv) ) &
              &         / tredsy(l,2)

         local_bgc_mem%sedfluxo(j,iv) = (local_bgc_mem%bgctra(j,kbo(j),iv_oc)-aprior)*local_bgc_mem%bolay(j)/dtbgc

         IF (lsediment_only) local_bgc_mem%bgctra(j,kbo(j),iv_oc) = aprior

      END DO
     endif

  END DO ! j loop
  !$ACC END PARALLEL

END SUBROUTINE DIPOWA

SUBROUTINE powadi (local_bgc_mem, j,  solrat, sedb1, sediso, bolven)
!! @file powadi.f90
!! @brief vertical diffusion with simultaneous dissolution,
!! implicit discretisation.
!!
!!

  USE mo_sedmnt, ONLY     : sedict, seddzi, seddw, porwah

  IMPLICIT NONE

  !! Arguments
  TYPE(t_bgc_memory), POINTER :: local_bgc_mem

  INTEGER, INTENT(in)     :: j             !< zonal grid index

  REAL(wp), INTENT(in)    :: solrat(ks)    !< dissolution rate
  REAL(wp), INTENT(inout) :: sedb1(0:ks)   !< tracer at entry
  REAL(wp), INTENT(inout) :: sediso(0:ks)  !<
  REAL(wp), INTENT(in)    :: bolven        !< bottom layer ventilation rate

  !! Local variables

  INTEGER,  SAVE  :: k,l   
  REAL(wp) :: tredsy(0:ks,3)
  REAL(wp), SAVE :: asu,alo
  !
  !----------------------------------------------------------------------
  !
  DO k = 1, ks

     asu=sedict*seddzi(k)*porwah(k)
     alo = 0._wp

     IF (k < ks) alo = sedict*seddzi(k+1)*porwah(k+1)

        tredsy(k,1) = -asu
        tredsy(k,3) = -alo
        tredsy(k,2) = seddw(k)*porwat(k) - tredsy(k,1)       &
             &        - tredsy(k,3) + solrat(k)*porwat(k)*seddw(k)

  END DO

  k = 0

  asu = 0._wp
  alo = sedict*seddzi(1)*porwah(1)

     IF(local_bgc_mem%bolay(j) > 0._wp) THEN
        tredsy(k,1) = -asu
        tredsy(k,3) = -alo
        tredsy(k,2) = bolven*local_bgc_mem%bolay(j) - tredsy(k,1) - tredsy(k,3)
     ELSE
        tredsy(k,1) = 0._wp
        tredsy(k,3) = 0._wp
        tredsy(k,2) = 0._wp
     ENDIF


  DO k = 1, ks
        IF (local_bgc_mem%bolay(j) > 0._wp) THEN
           tredsy(k-1,1) = tredsy(k,1) / tredsy(k-1,2)
           tredsy(k,2)   = tredsy(k,2)                       &
                &          - tredsy(k-1,3) * tredsy(k,1) / tredsy(k-1,2)
        ENDIF
  END DO

  DO k = 1, ks
        sedb1(k) = sedb1(k) - tredsy(k-1,1) * sedb1(k-1)
  END DO

  k = ks

     IF (local_bgc_mem%bolay(j) > 0._wp)sediso(k) = sedb1(k) / tredsy(k,2)

  DO k = 1, ks
     l = ks-k
        IF (local_bgc_mem%bolay(j) > 0._wp) sediso(l) =                           &
             &           ( sedb1(l) - tredsy(l,3) * sediso(l+1) )       &
             &           / tredsy(l,2)
  END DO

END SUBROUTINE powadi

END MODULE mo_sedmnt_diffusion
