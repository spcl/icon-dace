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


   MODULE mo_bgc_icon_comm

! icon specific routines for output, update etc

      USE mo_ocean_diagnostics_types,   ONLY:  t_ocean_regions

      USE mo_control_bgc,          ONLY: dtb,bgc_gin, bgc_arctic, bgc_lab, & 
       &                                 bgc_natl, bgc_atl, bgc_tatl, &
       &                                 bgc_tropac, &
       &                                 bgc_land, bgc_ind, &
       &                                 bgc_soce, bgc_npac, bgc_carb,inv_dtb

      USE mo_exception, ONLY      : message, finish,message_to_own_unit

      USE mo_model_domain,   ONLY: t_patch_3D, t_patch

      USE mo_ocean_nml,  ONLY: n_zlev, lsediment_only

      USE mo_kind,     ONLY: wp

      USE mo_impl_constants,      ONLY: max_char_length


      USE mo_grid_config,         ONLY: n_dom


      USE mo_hamocc_types,       ONLY: t_hamocc_diag, t_hamocc_state, &
    &                                  t_hamocc_sed, t_hamocc_tend,   &
    &                                  t_hamocc_monitor, t_hamocc_prog, t_hamocc_agg


      USE mo_bgc_constants,      ONLY: molw_co2

       
      USE mo_parallel_config,     ONLY: nproma

      USE mo_hamocc_nml,         ONLY: io_stdo_bgc, l_cpl_co2, ks, l_N_cycle,i_settling

      USE mo_bgc_memory_types,   ONLY: t_bgc_memory, t_sediment_memory, t_aggregates_memory

      USE mo_var_list_gpu,        ONLY: gpu_update_var_list

      USE mo_fortran_tools,      ONLY: set_acc_host_or_device
 
      IMPLICIT NONE

      PUBLIC

      TYPE(t_hamocc_state), TARGET                  :: hamocc_state


      INTERFACE to_bgcout
        module procedure to_bgcout_real
        module procedure to_bgcout_int
        module procedure to_bgcout_logical
      END INTERFACE



      CONTAINS

!================================================================================== 
       SUBROUTINE ini_bgc_regions

        IMPLICIT NONE
        TYPE(t_ocean_regions)         :: ocean_regions


         bgc_land=ocean_regions%land                                !0
         bgc_gin= ocean_regions%greenland_iceland_norwegian_sea     !1
         bgc_arctic=ocean_regions%arctic_ocean                      !2
         bgc_lab= ocean_regions%labrador_sea                        !3
         bgc_natl= ocean_regions%north_atlantic                     !4
         bgc_tatl= ocean_regions%tropical_atlantic                  !5 
         bgc_soce=ocean_regions%southern_ocean                      !6
         bgc_ind=ocean_regions%indian_ocean                         !7
         bgc_tropac= ocean_regions%tropical_pacific                 !8
         bgc_npac=ocean_regions%north_pacific                       !9
         bgc_carb=ocean_regions%caribbean                           !-33
 
       END SUBROUTINE
!================================================================================== 
    
      SUBROUTINE update_icon(local_bgc_mem, start_idx, end_idx, &
&             klevs, pddpo, jb, ptracer, pco2flx, lacc)

      USE mo_param1_bgc, ONLY: n_bgctra,kcflux_cpl

      TYPE(t_bgc_memory), POINTER :: local_bgc_mem
      REAL(wp)     :: ptracer(:,:,:,:)    
      INTEGER, INTENT(in)::klevs(nproma), jb
      REAL(wp),INTENT(in) :: pddpo(nproma,n_zlev) !< size of scalar grid cell (3rd REAL) [m]
      REAL(wp),INTENT(inout) :: pco2flx(nproma)
      LOGICAL, INTENT(IN), OPTIONAL :: lacc

      INTEGER :: jc, jk, kpke
      INTEGER :: start_idx, end_idx
      INTEGER :: itrac
      LOGICAL :: lzacc
      CHARACTER(LEN=max_char_length), PARAMETER :: &
                routine = 'update_icon'

     ! CALL message(TRIM(routine), 'start' )

     CALL set_acc_host_or_device(lzacc, lacc)
 
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc=start_idx,end_idx 
        kpke=klevs(jc)
        IF (pddpo(jc, 1) .GT. EPSILON(0.5_wp)) THEN
          pco2flx(jc)=local_bgc_mem%bgcflux(jc,kcflux_cpl) * molw_co2
        !$ACC LOOP SEQ
        DO jk =1,kpke
          !$ACC LOOP SEQ
          DO itrac=1,n_bgctra
             ptracer(jc,jk,jb,itrac) = local_bgc_mem%bgctra(jc,jk,itrac)
          ENDDO 
        ENDDO
        ENDIF
      ENDDO
      !$ACC END PARALLEL
 
      END SUBROUTINE

!================================================================================== 
      SUBROUTINE update_bgc(local_bgc_mem, local_sediment_mem, start_index, end_index, &
                          & klevs, pddpo, jb, ptracer, pco2mr, p_diag, p_sed, p_tend,  &
                          & max_klevs, lacc)

     USE mo_bgc_bcond,  ONLY: ext_data_bgc

      USE MO_PARAM1_BGC, ONLY: n_bgctra, issso12,         &
 &                             isssc12, issssil, issster, &
 &                             ipowaic, ipowaal, ipowaph, &
 &                             ipowaox, ipown2, ipowno3,  &
 &                             ipowasi, ipowafe, kn2b,    &
 &                             kh2ob,    &
&                              ipowh2s, &
&                              iatmco2, &
&                              ipownh4, ipowno2


      TYPE(t_bgc_memory), POINTER :: local_bgc_mem
      TYPE(t_sediment_memory), POINTER :: local_sediment_mem
      
      
      REAL(wp)     :: ptracer(:,:,:,:)    
      REAL(wp)     :: pco2mr(nproma)
      INTEGER, INTENT(in)::klevs(nproma), jb
      TYPE(t_hamocc_diag) :: p_diag
      TYPE(t_hamocc_sed) :: p_sed
      TYPE(t_hamocc_tend) :: p_tend
      REAL(wp),INTENT(in) :: pddpo(nproma,n_zlev) !< size of scalar grid cell (3rd REAL) [m]
      INTEGER, INTENT(IN) :: max_klevs
      LOGICAL, INTENT(IN), OPTIONAL :: lacc

      INTEGER :: jc, jk
      INTEGER :: start_index, end_index
      INTEGER :: itrac
      CHARACTER(LEN=max_char_length), PARAMETER :: &
                routine = 'update_icon'
      LOGICAL :: lzacc

      CALL set_acc_host_or_device(lzacc, lacc)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc=start_index,end_index
        IF (pddpo(jc, 1) .GT. EPSILON(0.5_wp)) THEN
          IF (l_cpl_co2) local_bgc_mem%atm(jc,iatmco2) = pco2mr(jc)
          local_bgc_mem%satn2(jc)  = p_tend%satn2(jc,jb)
          local_bgc_mem%satn2o(jc) = p_tend%satn2o(jc,jb)
          local_bgc_mem%solco2(jc) = p_tend%solco2(jc,jb)
          local_bgc_mem%kbo(jc) = p_sed%kbo(jc,jb)
          local_bgc_mem%bolay(jc) = p_sed%bolay(jc,jb)
          local_sediment_mem%burial(jc,issso12) = p_sed%bo12(jc,jb)
          local_sediment_mem%burial(jc,isssc12) = p_sed%bc12(jc,jb)
          local_sediment_mem%burial(jc,issssil) = p_sed%bsil(jc,jb)
          local_sediment_mem%burial(jc,issster) = p_sed%bter(jc,jb)
        END IF
      END DO

      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk =1,max_klevs
        DO jc=start_index,end_index
          IF (pddpo(jc, 1) .GT. EPSILON(0.5_wp) .and. jk <= klevs(jc)) THEN
            local_bgc_mem%hi(jc, jk)   = p_diag%hi(jc,jk,jb)
            local_bgc_mem%co3(jc, jk)  = p_diag%co3(jc,jk,jb)
            local_bgc_mem%aksp(jc, jk) = p_tend%aksp(jc,jk,jb)
            local_bgc_mem%ak13(jc, jk) = p_tend%ak1(jc,jk,jb)
            local_bgc_mem%ak23(jc, jk) = p_tend%ak2(jc,jk,jb)
            local_bgc_mem%akb3(jc, jk) = p_tend%akb(jc,jk,jb)
            local_bgc_mem%akw3(jc, jk) = p_tend%akw(jc,jk,jb)
            local_bgc_mem%aksi3(jc, jk) = p_tend%aksi(jc,jk,jb)
            local_bgc_mem%ak1p3(jc, jk) = p_tend%ak1p(jc,jk,jb)
            local_bgc_mem%ak2p3(jc, jk) = p_tend%ak2p(jc,jk,jb)
            local_bgc_mem%ak3p3(jc, jk) = p_tend%ak3p(jc,jk,jb)
            local_bgc_mem%aks3(jc, jk) = p_tend%aks(jc,jk,jb)
            local_bgc_mem%akf3(jc, jk) = p_tend%akf(jc,jk,jb)
            local_bgc_mem%satoxy(jc, jk) = p_tend%satoxy(jc,jk,jb)
            local_bgc_mem%bgctend(jc,jk,kh2ob) =  p_tend%h2obudget(jc,jk,jb)
            local_bgc_mem%bgctend(jc,jk,kn2b) =  p_tend%n2budget(jc,jk,jb)
          END IF
        END DO
      END DO

      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1,ks
        DO jc = start_index, end_index
          IF(pddpo(jc,1) .GT. EPSILON(0.5_wp)) THEN
            ! Solid sediment
            local_sediment_mem%sedlay(jc,jk,issso12) = p_sed%so12(jc,jk,jb)
            local_sediment_mem%sedlay(jc,jk,isssc12) = p_sed%sc12(jc,jk,jb)
            local_sediment_mem%sedlay(jc,jk,issssil) = p_sed%ssil(jc,jk,jb)
            local_sediment_mem%sedlay(jc,jk,issster) = p_sed%ster(jc,jk,jb)
            ! Pore water
            local_sediment_mem%powtra(jc,jk,ipowaic) = p_sed%pwic(jc,jk,jb)
            local_sediment_mem%powtra(jc,jk,ipowaal) = p_sed%pwal(jc,jk,jb)
            local_sediment_mem%powtra(jc,jk,ipowaph) = p_sed%pwph(jc,jk,jb)
            local_sediment_mem%powtra(jc,jk,ipowaox) = p_sed%pwox(jc,jk,jb)
            local_sediment_mem%powtra(jc,jk,ipowasi) = p_sed%pwsi(jc,jk,jb)
            local_sediment_mem%powtra(jc,jk,ipowafe) = p_sed%pwfe(jc,jk,jb)
            local_sediment_mem%powtra(jc,jk,ipown2)  = p_sed%pwn2(jc,jk,jb)
            local_sediment_mem%powtra(jc,jk,ipowno3) = p_sed%pwno3(jc,jk,jb)
            local_sediment_mem%powtra(jc,jk,ipowh2s) = p_sed%pwh2s(jc,jk,jb)

            local_sediment_mem%sedhpl(jc,jk)         = p_sed%sedhi(jc,jk,jb)
            local_sediment_mem%powh2obud(jc,jk)    = p_sed%pwh2ob(jc,jk,jb)
            local_sediment_mem%pown2bud(jc,jk)     = p_sed%pwn2b(jc,jk,jb)
          END IF
        END DO
      END DO

      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk =1,max_klevs
        DO jc=start_index,end_index
          IF (pddpo(jc, 1) .GT. EPSILON(0.5_wp) .and. jk <= klevs(jc)) THEN
            !$ACC LOOP SEQ
            DO itrac = 1,n_bgctra
              local_bgc_mem%bgctra(jc,jk,itrac)=ptracer(jc,jk,jb,itrac)
            END DO
          END IF
        END DO
      END DO

      IF (lsediment_only) THEN
        !$ACC LOOP GANG VECTOR
        DO jc=start_index,end_index
          IF (pddpo(jc, 1) .GT. EPSILON(0.5_wp)) THEN
            local_sediment_mem%silpro(jc) = ext_data_bgc%silpro(jc,jb)
            local_sediment_mem%produs(jc) = ext_data_bgc%produs(jc,jb)
            local_sediment_mem%prcaca(jc) = ext_data_bgc%prcaca(jc,jb)
            local_sediment_mem%prorca(jc) = ext_data_bgc%prorca(jc,jb)
          END IF
        END DO
      ENDIF

      IF (l_N_cycle) THEN
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1,ks
          DO jc = start_index, end_index
            IF( pddpo(jc,1) .GT. EPSILON(0.5_wp)) THEN
              local_sediment_mem%powtra(jc,jk,ipownh4) = p_sed%pwnh4(jc,jk,jb)
              local_sediment_mem%powtra(jc,jk,ipowno2) = p_sed%pwno2(jc,jk,jb)
            END IF
          END DO
        END DO
      END IF

      !$ACC END PARALLEL

      END SUBROUTINE

!================================================================================== 
  SUBROUTINE set_bgc_tendencies_output(local_bgc_mem, local_sediment_mem, local_aggregate_memory, start_idx, end_idx, &
&             klevs,pddpo,jb,p_tend, p_diag, p_sed, p_agg, max_klevs, lacc)
      
      USE mo_param1_bgc, ONLY: kphosy, ksred, kremin, kdenit, &
 &                             kcflux, koflux, knflux, knfixd, &
 &                             knfix, kgraz, ksilpro, kprorca, &
 &                             kn2oflux,korginp,ksilinp, kcalinp,  &
 &                             kprcaca, kcoex90,issso12, &
 &                             isssc12, issssil, issster, &
 &                             ipowaic, ipowaal, ipowaph, &
 &                             ipowaox, ipown2, ipowno3,  &
 &                             ipownh4, ipowno2, &
 &                             ipowasi, ipowafe, kpho_cya, &
 &                             kcyaloss, kn2b, kh2ob, kprodus, &
&                              kbacfra, kdelsil, kdelcar,  &
&                              kdmsflux, kdmsprod, kdmsbac, kdmsuv, &
&                              keuexp, kplim, kflim, knlim,kcalex90,&
&                              kopex90, kgraton, kexudp, kexudz, &
&                              kzdy, kpdy,kcoex1000,kcoex2000, &
&                              kopex1000,kopex2000,kcalex1000,&
&                              kcalex2000, kaou, kcTlim, kcLlim, &
&                              kcPlim, kcFlim, ipowh2s,kh2sprod, &   
&                              kh2sloss,iatmco2,kpco2,klysocl,knitinp, &
&                              kwdust, kwpoc, kwopal, kwcal, &
&                              isremino, isreminn, isremins, &
&                              knh3flux, kgppnh, kcyapro, kammox, knitox, &
&                              kdnrn, kdnra, kanam, ksammox, ksnitox, &
&                              ksdnrn, ksdnra, ksanam, ksnrn2

      USE mo_memory_agg, ONLY : kavdp, kavrhop, ksticka, klmaxagg, kdfagg, kavrhof
  
      TYPE(t_bgc_memory), POINTER :: local_bgc_mem
      TYPE(t_sediment_memory), POINTER :: local_sediment_mem
      TYPE(t_aggregates_memory), POINTER :: local_aggregate_memory

      TYPE(t_hamocc_tend) :: p_tend
      TYPE(t_hamocc_diag) :: p_diag
      TYPE(t_hamocc_sed) :: p_sed
      TYPE(t_hamocc_agg) :: p_agg


      INTEGER, INTENT(in) :: klevs(nproma)
      INTEGER, INTENT(in) :: start_idx, end_idx,jb
      INTEGER, INTENT(IN) :: max_klevs
      REAL(wp),INTENT(in) :: pddpo(nproma,n_zlev) !< size of scalar grid cell (3rd REAL) [m]
      LOGICAL, INTENT(IN), OPTIONAL :: lacc


      INTEGER :: jc, jk
      LOGICAL :: lzacc

      CALL set_acc_host_or_device(lzacc, lacc)
 
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc=start_idx,end_idx 
        IF (pddpo(jc, 1) .GT. EPSILON(0.5_wp)) THEN
          p_tend%co2mr(jc,jb) = local_bgc_mem%atm(jc,iatmco2)
          p_tend%cflux(jc,jb) = local_bgc_mem%bgcflux(jc,kcflux)
          p_tend%pco2(jc,jb)  = local_bgc_mem%bgcflux(jc,kpco2)
          p_tend%oflux(jc,jb) = local_bgc_mem%bgcflux(jc,koflux)
          p_tend%nflux(jc,jb) = local_bgc_mem%bgcflux(jc,knflux)
          p_tend%dmsflux(jc,jb) = local_bgc_mem%bgcflux(jc,kdmsflux)
          p_tend%orginp(jc,jb) = local_bgc_mem%bgcflux(jc,korginp)
          p_tend%silinp(jc,jb) = local_bgc_mem%bgcflux(jc,ksilinp)
          p_tend%calinp(jc,jb) = local_bgc_mem%bgcflux(jc,kcalinp)
          p_tend%n2oflux(jc,jb) = local_bgc_mem%bgcflux(jc,kn2oflux)
          p_tend%nfixd(jc,jb) = local_bgc_mem%bgcflux(jc,knfixd)
          p_tend%prcaca(jc,jb) = local_bgc_mem%bgcflux(jc,kprcaca)
          p_tend%prorca(jc,jb) = local_bgc_mem%bgcflux(jc,kprorca)
          p_tend%silpro(jc,jb) = local_bgc_mem%bgcflux(jc,ksilpro)
          p_tend%produs(jc,jb) = local_bgc_mem%bgcflux(jc,kprodus)
          p_tend%coex90(jc,jb) = local_bgc_mem%bgcflux(jc,kcoex90)
          p_tend%calex90(jc,jb) = local_bgc_mem%bgcflux(jc,kcalex90)
          p_tend%opex90(jc,jb) = local_bgc_mem%bgcflux(jc,kopex90)
          p_tend%coex1000(jc,jb) = local_bgc_mem%bgcflux(jc,kcoex1000)
          p_tend%calex1000(jc,jb) = local_bgc_mem%bgcflux(jc,kcalex1000)
          p_tend%opex1000(jc,jb) = local_bgc_mem%bgcflux(jc,kopex1000)
          p_tend%coex2000(jc,jb) = local_bgc_mem%bgcflux(jc,kcoex2000)
          p_tend%calex2000(jc,jb) = local_bgc_mem%bgcflux(jc,kcalex2000)
          p_tend%opex2000(jc,jb) = local_bgc_mem%bgcflux(jc,kopex2000)
          p_tend%lysocline(jc,jb)=local_bgc_mem%bgcflux(jc,klysocl)
          p_tend%nitrogeninp(jc,jb)=local_bgc_mem%bgcflux(jc,knitinp)
        END IF
      END DO

      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk =1,max_klevs
        DO jc=start_idx,end_idx
          IF (pddpo(jc, 1) .GT. EPSILON(0.5_wp) .and. jk <= klevs(jc)) THEN
             p_tend%npp(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kphosy)
             p_tend%graz(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kgraz)
             p_tend%zoomor(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kzdy)
             p_tend%phymor(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kpdy)
             p_tend%exudz(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kexudz)
             p_tend%graton(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kgraton)
             p_tend%exud(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kexudp)
             p_tend%nfix(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,knfix)
             p_tend%phoc(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kpho_cya)
             p_tend%cyloss(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kcyaloss)
             p_tend%bacfra(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kbacfra)
             p_tend%remina(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kremin)
             p_tend%remins(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,ksred)
             p_tend%reminn(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kdenit)
             p_tend%h2obudget(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kh2ob)
             p_tend%n2budget(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kn2b)
             p_tend%delsil(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kdelsil)
             p_tend%delcar(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kdelcar)
             p_tend%h2sprod(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kh2sprod)
             p_tend%h2sloss(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kh2sloss)
             p_tend%dmsprod(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kdmsprod)
             p_tend%dmsbac(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kdmsbac)
             p_tend%dmsuv(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kdmsuv)
             p_tend%euexp(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,keuexp)
             p_diag%hi(jc,jk,jb)     = local_bgc_mem%hi(jc, jk)
             p_diag%co3(jc,jk,jb)    = local_bgc_mem%co3(jc,jk) 
             p_tend%akb(jc,jk,jb)    = local_bgc_mem%akb3(jc,jk) 
             p_tend%akw(jc,jk,jb)    = local_bgc_mem%akw3(jc,jk) 
             p_tend%ak1(jc,jk,jb)    = local_bgc_mem%ak13(jc,jk) 
             p_tend%ak2(jc,jk,jb)    = local_bgc_mem%ak23(jc,jk) 
             p_tend%aks(jc,jk,jb)    = local_bgc_mem%aks3(jc,jk) 
             p_tend%akf(jc,jk,jb)    = local_bgc_mem%akf3(jc,jk) 
             p_tend%ak1p(jc,jk,jb)   = local_bgc_mem%ak1p3(jc,jk) 
             p_tend%ak2p(jc,jk,jb)   = local_bgc_mem%ak2p3(jc,jk)
             p_tend%ak3p(jc,jk,jb)   = local_bgc_mem%ak3p3(jc,jk)
             p_tend%aksi(jc,jk,jb)   = local_bgc_mem%aksi3(jc,jk)
             p_tend%aksp(jc,jk,jb)   = local_bgc_mem%aksp(jc,jk) 
             p_tend%flim(jc,jk,jb)   = local_bgc_mem%bgctend(jc,jk,kflim)
             p_tend%nlim(jc,jk,jb)   = local_bgc_mem%bgctend(jc,jk,knlim)
             p_tend%plim(jc,jk,jb)   = local_bgc_mem%bgctend(jc,jk,kplim)
             p_tend%cTlim(jc,jk,jb)  = local_bgc_mem%bgctend(jc,jk,kcTlim)
             p_tend%cLlim(jc,jk,jb)  = local_bgc_mem%bgctend(jc,jk,kcLlim)
             p_tend%cPlim(jc,jk,jb)  = local_bgc_mem%bgctend(jc,jk,kcPlim)
             p_tend%cFlim(jc,jk,jb)  = local_bgc_mem%bgctend(jc,jk,kcFlim)
             p_tend%satoxy(jc,jk,jb)    = local_bgc_mem%satoxy(jc,jk) 
             p_tend%aou(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kaou)
           END IF
        END DO
      END DO

      IF(i_settling==2)then
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1,max_klevs
          DO jc=start_idx,end_idx
            IF (pddpo(jc, 1) .GT. EPSILON(0.5_wp) .and. jk <= klevs(jc)) THEN
              p_tend%wdust(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kwdust)
              p_tend%wpoc(jc,jk,jb)  = local_bgc_mem%bgctend(jc,jk,kwpoc)
              p_tend%wopal(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kwopal)
              p_tend%wcal(jc,jk,jb)  = local_bgc_mem%bgctend(jc,jk,kwcal)
              p_agg%avdp(jc,jk,jb)   = local_aggregate_memory%aggdiag(jc,jk,kavdp)
              p_agg%avrhop(jc,jk,jb) = local_aggregate_memory%aggdiag(jc,jk,kavrhop)
              p_agg%sticka(jc,jk,jb) = local_aggregate_memory%aggdiag(jc,jk,ksticka)
              p_agg%lmaxagg(jc,jk,jb) = local_aggregate_memory%aggdiag(jc,jk,klmaxagg)
              p_agg%dfagg(jc,jk,jb)  = local_aggregate_memory%aggdiag(jc,jk,kdfagg)
            END IF
          END DO
        END DO
      END IF
                 

      IF (l_N_cycle) THEN
        !$ACC LOOP GANG VECTOR
        DO jc=start_idx,end_idx
          IF (pddpo(jc, 1) .GT. EPSILON(0.5_wp)) THEN
            p_tend%nh3flux(jc,jb) = local_bgc_mem%bgcflux(jc,knh3flux)
            p_tend%sedflnh4(jc,jb) = local_bgc_mem%sedfluxo(jc,ipownh4)
            p_tend%sedflno2(jc,jb) = local_bgc_mem%sedfluxo(jc,ipowno2)
          END IF
        END DO

        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1, max_klevs
          DO jc = start_idx, end_idx
            IF( pddpo(jc,1) .GT. EPSILON(0.5_wp) .and. jk <= klevs(jc)) THEN
              p_tend%gppnh4(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kgppnh)
              p_tend%cyapro(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kcyapro)
              p_tend%ammox(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kammox)
              p_tend%nitox(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,knitox)
              p_tend%dnrn(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kdnrn)
              p_tend%dnra(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kdnra)
              p_tend%anam(jc,jk,jb) = local_bgc_mem%bgctend(jc,jk,kanam)
            END IF
          END DO
        END DO

        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1,ks
          DO jc = start_idx, end_idx
            IF( pddpo(jc,1) .GT. EPSILON(0.5_wp)) THEN
              p_sed%pwnh4(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipownh4) 
              p_sed%pwno2(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowno2)
              p_tend%sedammox(jc,jk,jb) = local_sediment_mem%sedtend(jc,jk,ksammox)
              p_tend%sednitox(jc,jk,jb) = local_sediment_mem%sedtend(jc,jk,ksnitox)
              p_tend%seddnrn(jc,jk,jb) = local_sediment_mem%sedtend(jc,jk,ksdnrn)
              p_tend%seddnra(jc,jk,jb) = local_sediment_mem%sedtend(jc,jk,ksdnra)
              p_tend%sedanam(jc,jk,jb) = local_sediment_mem%sedtend(jc,jk,ksanam)
              p_tend%sednrn2(jc,jk,jb) = local_sediment_mem%sedtend(jc,jk,ksnrn2)
            END IF
          END DO
        END DO
      END IF

      !$ACC LOOP GANG VECTOR
      DO jc=start_idx,end_idx
        IF (pddpo(jc, 1) .GT. EPSILON(0.5_wp)) THEN
          p_tend%satn2(jc,jb)    = local_bgc_mem%satn2(jc) 
          p_tend%satn2o(jc,jb)   = local_bgc_mem%satn2o(jc) 
          p_tend%solco2(jc,jb)   = local_bgc_mem%solco2(jc) 
          ! Sediment
          ! local_sediment_mem%burial layers
          p_sed%bo12(jc,jb) = local_sediment_mem%burial(jc,issso12)
          p_sed%bc12(jc,jb) = local_sediment_mem%burial(jc,isssc12)
          p_sed%bsil(jc,jb) = local_sediment_mem%burial(jc,issssil)
          p_sed%bter(jc,jb) = local_sediment_mem%burial(jc,issster) 
          ! Sediment-ocean fluxes
          p_tend%sedflic(jc,jb) = local_bgc_mem%sedfluxo(jc,ipowaic) 
          p_tend%sedflal(jc,jb) = local_bgc_mem%sedfluxo(jc,ipowaal) 
          p_tend%sedflph(jc,jb) = local_bgc_mem%sedfluxo(jc,ipowaph) 
          p_tend%sedflox(jc,jb) = local_bgc_mem%sedfluxo(jc,ipowaox) 
          p_tend%sedflsi(jc,jb) = local_bgc_mem%sedfluxo(jc,ipowasi) 
          p_tend%sedflfe(jc,jb) = local_bgc_mem%sedfluxo(jc,ipowafe) 
          p_tend%sedfln2(jc,jb) = local_bgc_mem%sedfluxo(jc,ipown2) 
          p_tend%sedflno3(jc,jb) = local_bgc_mem%sedfluxo(jc,ipowno3) 
          p_tend%sedflh2s(jc,jb) = local_bgc_mem%sedfluxo(jc,ipowh2s)
        END IF
      END DO

      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk =1,ks
        DO jc=start_idx,end_idx
          IF (pddpo(jc, 1) .GT. EPSILON(0.5_wp)) THEN
            ! Solid sediment
            p_sed%so12(jc,jk,jb) = local_sediment_mem%sedlay(jc,jk,issso12)
            p_sed%sc12(jc,jk,jb) = local_sediment_mem%sedlay(jc,jk,isssc12)
            p_sed%ssil(jc,jk,jb) = local_sediment_mem%sedlay(jc,jk,issssil)
            p_sed%ster(jc,jk,jb) = local_sediment_mem%sedlay(jc,jk,issster)
            ! Pore water
            p_sed%pwic(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowaic)
            p_sed%pwal(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowaal)
            p_sed%pwph(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowaph)
            p_sed%pwox(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowaox)
            p_sed%pwsi(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowasi)
            p_sed%pwfe(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowafe)
            p_sed%pwn2(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipown2)
            p_sed%pwno3(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowno3)
            p_sed%pwh2s(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowh2s)
            p_sed%sedhi(jc,jk,jb) = local_sediment_mem%sedhpl(jc,jk)
            p_sed%pwh2ob(jc,jk,jb) = local_sediment_mem%powh2obud(jc,jk)
            p_sed%pwn2b(jc,jk,jb) = local_sediment_mem%pown2bud(jc,jk)
            ! tendencies
            p_tend%sedro2(jc,jk,jb) = local_sediment_mem%sedtend(jc,jk,isremino) 
            p_tend%sedrn(jc,jk,jb)  = local_sediment_mem%sedtend(jc,jk,isreminn) 
            p_tend%sedrs(jc,jk,jb)  = local_sediment_mem%sedtend(jc,jk,isremins) 
          END IF
        END DO
      END DO
      !$ACC END PARALLEL

  END SUBROUTINE 

!================================================================================== 
  SUBROUTINE set_bgc_tendencies_output_sedon(local_bgc_mem, local_sediment_mem, start_idx,end_idx,pddpo,jb, &
 &                                           p_tend, p_sed)

      USE mo_param1_bgc, ONLY: issso12, &
 &                             isssc12, issssil, issster, &
 &                             ipowaic, ipowaal, ipowaph, &
 &                             ipowaox, ipown2, ipowno3,  &
 &                             ipownh4, ipowno2, &
 &                             ipowasi, ipowafe, ipowh2s, &
 &                             isremino, isreminn, isremins, &
 &                             ksammox, ksnitox, ksdnrn, ksdnra, &
 &                             ksanam, ksnrn2

      USE mo_bgc_bcond, ONLY: ext_data_bgc
      TYPE(t_bgc_memory), POINTER :: local_bgc_mem
      TYPE(t_sediment_memory), POINTER :: local_sediment_mem
  
      TYPE(t_hamocc_sed) :: p_sed
      TYPE(t_hamocc_tend):: p_tend

      INTEGER, INTENT(in) :: start_idx, end_idx,jb
      REAL(wp),INTENT(in) :: pddpo(nproma,n_zlev) !< size of scalar grid cell (3rd REAL) [m]


      INTEGER :: jc, jk

      DO jc=start_idx,end_idx 
        IF (pddpo(jc, 1) .GT. EPSILON(0.5_wp)) THEN
        p_tend%prcaca(jc,jb) = ext_data_bgc%prcaca(jc,jb)
        p_tend%prorca(jc,jb) = ext_data_bgc%prorca(jc,jb)
        p_tend%silpro(jc,jb) = ext_data_bgc%silpro(jc,jb)
        p_tend%produs(jc,jb) = ext_data_bgc%produs(jc,jb)
        ! Sediment
        ! local_sediment_mem%burial layers
        p_sed%bo12(jc,jb) = local_sediment_mem%burial(jc,issso12)
        p_sed%bc12(jc,jb) = local_sediment_mem%burial(jc,isssc12)
        p_sed%bsil(jc,jb) = local_sediment_mem%burial(jc,issssil)
        p_sed%bter(jc,jb) = local_sediment_mem%burial(jc,issster) 
        ! Sediment-ocean fluxes
        p_tend%sedflic(jc,jb) = local_bgc_mem%sedfluxo(jc,ipowaic) 
        p_tend%sedflal(jc,jb) = local_bgc_mem%sedfluxo(jc,ipowaal) 
        p_tend%sedflph(jc,jb) = local_bgc_mem%sedfluxo(jc,ipowaph) 
        p_tend%sedflox(jc,jb) = local_bgc_mem%sedfluxo(jc,ipowaox) 
        p_tend%sedflsi(jc,jb) = local_bgc_mem%sedfluxo(jc,ipowasi) 
        p_tend%sedflfe(jc,jb) = local_bgc_mem%sedfluxo(jc,ipowafe) 
        p_tend%sedfln2(jc,jb) = local_bgc_mem%sedfluxo(jc,ipown2) 
        p_tend%sedflno3(jc,jb) = local_bgc_mem%sedfluxo(jc,ipowno3) 
        p_tend%sedflh2s(jc,jb) = local_bgc_mem%sedfluxo(jc,ipowh2s) 
        if (l_N_cycle) THEN
           p_tend%sedflnh4(jc,jb) = local_bgc_mem%sedfluxo(jc,ipownh4) 
           p_tend%sedflno2(jc,jb) = local_bgc_mem%sedfluxo(jc,ipowno2)
           DO jk = 1,ks
                p_sed%pwnh4(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipownh4) 
                p_sed%pwno2(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowno2)
                p_tend%sedammox(jc,jk,jb) = local_sediment_mem%sedtend(jc,jk,ksammox)
                p_tend%sednitox(jc,jk,jb) = local_sediment_mem%sedtend(jc,jk,ksnitox)
                p_tend%seddnrn(jc,jk,jb) = local_sediment_mem%sedtend(jc,jk,ksdnrn)
                p_tend%seddnra(jc,jk,jb) = local_sediment_mem%sedtend(jc,jk,ksdnra)
                p_tend%sedanam(jc,jk,jb) = local_sediment_mem%sedtend(jc,jk,ksanam)
                p_tend%sednrn2(jc,jk,jb) = local_sediment_mem%sedtend(jc,jk,ksnrn2)
           ENDDO
        ENDIF
        DO jk =1,ks
             ! Solid sediment
             p_sed%so12(jc,jk,jb) = local_sediment_mem%sedlay(jc,jk,issso12)
             p_sed%sc12(jc,jk,jb) = local_sediment_mem%sedlay(jc,jk,isssc12)
             p_sed%ssil(jc,jk,jb) = local_sediment_mem%sedlay(jc,jk,issssil)
             p_sed%ster(jc,jk,jb) = local_sediment_mem%sedlay(jc,jk,issster)
             ! Pore water
             p_sed%pwic(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowaic)
             p_sed%pwal(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowaal)
             p_sed%pwph(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowaph)
             p_sed%pwox(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowaox)
             p_sed%pwsi(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowasi)
             p_sed%pwfe(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowafe)
             p_sed%pwn2(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipown2)
             p_sed%pwno3(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowno3)
             p_sed%pwh2s(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowh2s)
             p_sed%sedhi(jc,jk,jb) = local_sediment_mem%sedhpl(jc,jk)
             p_sed%pwh2ob(jc,jk,jb) = local_sediment_mem%powh2obud(jc,jk)
             p_sed%pwn2b(jc,jk,jb) = local_sediment_mem%pown2bud(jc,jk)

             ! tendencies
             p_tend%sedro2(jc,jk,jb) = local_sediment_mem%sedtend(jc,jk,isremino) 
             p_tend%sedrn(jc,jk,jb)  = local_sediment_mem%sedtend(jc,jk,isreminn) 
             p_tend%sedrs(jc,jk,jb)  = local_sediment_mem%sedtend(jc,jk,isremins) 
             
        ENDDO
      ENDIF
      ENDDO

 

  END SUBROUTINE 

!================================================================================== 
    
      SUBROUTINE initial_update_icon(local_bgc_mem, local_sediment_mem, start_index, end_index, &
&             klevs,pddpo,jb,ptracer, p_sed,p_diag,pco2flux)

      USE mo_param1_bgc, ONLY: n_bgctra, issso12, &
 &                             isssc12, issssil, issster, &
 &                             ipowaic, ipowaal, ipowaph, &
 &                             ipowaox, ipown2, ipowno3,  &
 &                             ipownh4, ipowno2, &
 &                             ipowasi, ipowafe, ipowh2s, &
 &                             kcflux_cpl
  
      TYPE(t_bgc_memory), POINTER :: local_bgc_mem
      TYPE(t_sediment_memory), POINTER :: local_sediment_mem
      REAL(wp)     :: ptracer(nproma,n_zlev,n_bgctra)    

      INTEGER, INTENT(in)::klevs(nproma)
      TYPE(t_hamocc_sed) :: p_sed
      TYPE(t_hamocc_diag) :: p_diag
      REAL(wp),INTENT(in) :: pddpo(nproma,n_zlev) !< size of scalar grid cell (3rd REAL) [m]
      REAL(wp),INTENT(inout)  :: pco2flux(nproma)

      INTEGER :: jc, jk, kpke,jb
      INTEGER :: start_index, end_index
      INTEGER :: itrac
      CHARACTER(LEN=max_char_length), PARAMETER :: &
                routine = 'update_icon'

     ! CALL message(TRIM(routine), 'start' )
 
      DO jc=start_index,end_index 
        kpke=klevs(jc)
        IF (pddpo(jc, 1) .GT. EPSILON(0.5_wp)) THEN
         if(l_cpl_co2)pco2flux(jc)=local_bgc_mem%bgcflux(jc,kcflux_cpl) * molw_co2
        DO jk =1,kpke
          DO itrac=1,n_bgctra
             ptracer(jc,jk,itrac) = local_bgc_mem%bgctra(jc,jk,itrac)
          ENDDO

             p_diag%hi(jc,jk,jb)  = local_bgc_mem%hi(jc, jk)
             p_diag%co3(jc,jk,jb) = local_bgc_mem%co3(jc,jk) 
        ENDDO
        ! Sediment
        ! local_sediment_mem%burial layers
        p_sed%bo12(jc,jb) = local_sediment_mem%burial(jc,issso12)
        p_sed%bc12(jc,jb) = local_sediment_mem%burial(jc,isssc12)
        p_sed%bsil(jc,jb) = local_sediment_mem%burial(jc,issssil)
        p_sed%bter(jc,jb) = local_sediment_mem%burial(jc,issster) 
        p_sed%bolay(jc,jb) = local_bgc_mem%bolay(jc)
        p_sed%kbo(jc,jb) = local_bgc_mem%kbo(jc)
        DO jk =1,ks
             ! Solid sediment
             p_sed%so12(jc,jk,jb) = local_sediment_mem%sedlay(jc,jk,issso12)
             p_sed%sc12(jc,jk,jb) = local_sediment_mem%sedlay(jc,jk,isssc12)
             p_sed%ssil(jc,jk,jb) = local_sediment_mem%sedlay(jc,jk,issssil)
             p_sed%ster(jc,jk,jb) = local_sediment_mem%sedlay(jc,jk,issster)
             ! Pore water
             p_sed%pwic(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowaic)
             p_sed%pwal(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowaal)
             p_sed%pwph(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowaph)
             p_sed%pwox(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowaox)
             p_sed%pwsi(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowasi)
             p_sed%pwfe(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowafe)
             p_sed%pwn2(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipown2)
             p_sed%pwno3(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowno3)
             p_sed%pwh2s(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowh2s)
             if (l_N_cycle) THEN
                p_sed%pwnh4(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipownh4) 
                p_sed%pwno2(jc,jk,jb) = local_sediment_mem%powtra(jc,jk,ipowno2) 
             ENDIF
             p_sed%sedhi(jc,jk,jb) = local_sediment_mem%sedhpl(jc,jk)
        ENDDO
      ENDIF
 
      ENDDO
 
 
      END SUBROUTINE

!================================================================================== 


  SUBROUTINE print_bgc_parameters
  USE mo_memory_bgc, ONLY      : phytomi, grami, remido, dyphy, zinges,        &
       &                      bkphy, bkzoo, bkopal,                      &
       &                     sulfate_reduction,                &
       &                     n2_fixation,                             &
       &                     ropal, perc_diron, riron, fesoly, relaxfe,         &
       &                     pi_alpha_cya,                     &
       &                     Topt_cya,T1_cya,T2_cya,bkcya_N, &
       &                     buoyancyspeed_cya,                                 &
       &                     doccya_fac, thresh_aerob, thresh_sred, &
       &                     bkno3, bkfe, bkpo4, bkno3_cya, bknh4, bknh4_cya, &
       &                     ro2ammo, rmm, kg_denom, no2denit, anamoxra, &
       &                     nitriox, nitrira, bkno2, rno3nh4, rno3no2

   USE mo_hamocc_nml, ONLY: l_cyadyn, denit_sed, disso_po, &
      &                 sinkspeed_opal, sinkspeed_calc,grazra,cycdec,l_dynamic_pi, &
      &                 drempoc,dremopal,dremcalc,denitrification, bkcya_P,bkcya_fe, &
      &                 no3nh4red, no3no2red, calmax, &
      &                    l_opal_q10, opal_remin_q10, opal_remin_tref,&
      &                    l_doc_q10, doc_remin_q10, doc_remin_tref,&
      &                    l_poc_q10, poc_remin_q10, poc_remin_tref


   USE mo_sedmnt, ONLY: disso_op, disso_cal,sred_sed

  
  CHARACTER(LEN=max_char_length) :: &
                cpara_name,cpara_val

   cpara_name='========PARAMETER SETUP'
   cpara_val="============="
   CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )

   ! Phytoplankton
   cpara_name='PHYTOPLANKTON'
   cpara_val="========"
   CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   CALL to_bgcout("phytomi",phytomi)
   CALL to_bgcout("dyphy",dyphy)
   CALL to_bgcout("bkphy",bkphy)
   CALL to_bgcout("l_dynamic_pi",l_dynamic_pi)

   ! Zooplankton
   cpara_name='ZOOPLANKTON'
   cpara_val="========"
   CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   CALL to_bgcout("grami",grami)
   CALL to_bgcout("zinges",zinges)
   CALL to_bgcout("grazra",grazra*inv_dtb)
   CALL to_bgcout("bkzoo",bkzoo)

   ! Cyanobacteria
   cpara_name='CYANOBACTERIA'
   cpara_val="========"
   CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   CALL to_bgcout("l_cyadyn",l_cyadyn)
   CALL to_bgcout("n2_fixation",n2_fixation)
   CALL to_bgcout("n2_fixation",n2_fixation)
   CALL to_bgcout("buoyancyspeed_cya 1/d",buoyancyspeed_cya*inv_dtb)
   CALL to_bgcout("cycdec 1/d",cycdec*inv_dtb)
   CALL to_bgcout("pi_alpha_cya",pi_alpha_cya)
   CALL to_bgcout("Topt_cya",Topt_cya)
   CALL to_bgcout("T1",T1_cya)
   CALL to_bgcout("T2",T2_cya)
   CALL to_bgcout("bkcya_p",bkcya_p)
   CALL to_bgcout("bkcya_n",bkcya_n)
   CALL to_bgcout("bkcya_fe",bkcya_fe)
   CALL to_bgcout("doccya_fac",doccya_fac)
  
   ! Detritus
   cpara_name='DETRITUS'
   cpara_val="========"
   CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   CALL to_bgcout("i_settling",i_settling)
   CALL to_bgcout("drempoc 1/d",drempoc*inv_dtb)
   CALL to_bgcout("denitrification 1/d",denitrification*inv_dtb)
   CALL to_bgcout("sulfate_reduction 1/d",sulfate_reduction)
   CALL to_bgcout("thresh_aerob",thresh_aerob)
   CALL to_bgcout("thresh_sred",thresh_sred)
   
   IF (i_settling==2) THEN
       CALL to_bgcout("l_poc_q10",l_poc_q10)     
       CALL to_bgcout("poc_remin_tref",poc_remin_tref)
       CALL to_bgcout("poc_remin_q10",poc_remin_q10)
   END IF

   ! extended N-cycle
   if (l_N_cycle) then

   cpara_name='Ncycle'
   cpara_val="========"
   CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   CALL to_bgcout("bkno3",bkno3)
   CALL to_bgcout("bkfe",bkfe)
   CALL to_bgcout("bkpo4",bkpo4)
   CALL to_bgcout("bkno3_cya",bkno3_cya)
   CALL to_bgcout("bknh4",bknh4)
   CALL to_bgcout("bknh4_cya",bknh4_cya)
   CALL to_bgcout("ro2ammo",ro2ammo)
   CALL to_bgcout("rmm",rmm)
   CALL to_bgcout("kg_denom",kg_denom)
   CALL to_bgcout("no3no2red 1/d", no3no2red*inv_dtb)
   CALL to_bgcout("no3nh4red 1/d", no3nh4red*inv_dtb)
   CALL to_bgcout("no2denit 1/d", no2denit*inv_dtb)
   CALL to_bgcout("anamoxra 1/d", anamoxra*inv_dtb)
   CALL to_bgcout("nitriox 1/d", nitriox*inv_dtb)
   CALL to_bgcout("nitrira 1/d", nitrira*inv_dtb)
   CALL to_bgcout("bkno2",bkno2)
   CALL to_bgcout("rno3nh4",rno3nh4)
   CALL to_bgcout("rno3no2",rno3no2)

   endif


   ! DOC
   cpara_name='DOC'
   cpara_val="========"
   CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   CALL to_bgcout("remido",remido*inv_dtb)
   IF (i_settling==2) THEN
       CALL to_bgcout("l_doc_q10",l_doc_q10)     
       CALL to_bgcout("doc_remin_tref",doc_remin_tref)
       CALL to_bgcout("doc_remin_q10",doc_remin_q10)
   END IF
  

   ! Opal
   cpara_name='Si/Opal'
   cpara_val="========"
   CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   CALL to_bgcout("bkopal",bkopal)
   CALL to_bgcout("dremopal 1/d",dremopal*inv_dtb)
   CALL to_bgcout("ropal",ropal)
   CALL to_bgcout("sinkspeed_opal",sinkspeed_opal*inv_dtb)
   IF (i_settling==2) THEN
       CALL to_bgcout("l_opal_q10",l_opal_q10)     
       CALL to_bgcout("opal_remin_tref",opal_remin_tref)
       CALL to_bgcout("opal_remin_q10",opal_remin_q10)
   END IF
   CALL  to_bgcout("calmax",calmax)
   ! Iron
   cpara_name='Iron'
   cpara_val="========"
   CALL to_bgcout("perc_diron",perc_diron)
   CALL to_bgcout("fesoly",fesoly)
   CALL to_bgcout("riron",riron)
   CALL to_bgcout("relaxfe",relaxfe)

   ! Sediment
   cpara_name='Sediment'
   cpara_val="========"
   CALL to_bgcout("denit_sed",denit_sed)
   CALL to_bgcout("disso_po",disso_po)
   CALL to_bgcout("disso_op",disso_op)
   CALL to_bgcout("disso_cal",disso_cal)
   CALL to_bgcout("sred_sed",sred_sed)
   

   ! Calc
   cpara_name='Calc'
   cpara_val="========"
   CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   CALL to_bgcout("dremcalc 1/d",dremcalc*inv_dtb)
   CALL to_bgcout("sinkspeed_calc",sinkspeed_calc*inv_dtb)

   cpara_name='======================='
   cpara_val="==========="
   CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
  END SUBROUTINE print_bgc_parameters
!================================================================================== 

  SUBROUTINE print_wpoc(wpoc)
    REAL(wp), INTENT(in) :: wpoc(:,:)
     
    CHARACTER(LEN=max_char_length) :: &
                cpara_name,cpara_val

    INTEGER:: k

 ! CALL to_bgcout("wdust",wdust)
   cpara_name='========WPOC [m/d]'
   cpara_val="============="
   CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   
   DO k=1,n_zlev
    write(cpara_name,'(i2)')k
    write(cpara_val,'(f6.2)')wpoc(1,k)/dtb
    CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
   enddo

   cpara_name='======================='
   cpara_val="==========="
   CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
  
  END SUBROUTINE print_wpoc

!================================================================================== 

SUBROUTINE to_bgcout_real(cname,val)
  REAL(wp),INTENT(in) ::val
  CHARACTER( LEN = * ):: cname
  CHARACTER(LEN=max_char_length) :: cpara_name, cpara_val

  cpara_name=cname
  IF(abs(val)<100._wp.and.abs(val)>0.1_wp)then
    write(cpara_val, '(f9.2)') val
  ELSE
    write(cpara_val, '(ES22.15)') val
  ENDIF
  CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc)
END SUBROUTINE

SUBROUTINE to_bgcout_logical(cname,val)
  LOGICAL,INTENT(in) ::val
  CHARACTER( LEN = * ):: cname
  CHARACTER(LEN=max_char_length) :: cpara_name, cpara_val

  cpara_name=cname
  cpara_val=merge(".true. ",".false.",val)
  CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc)
END SUBROUTINE

SUBROUTINE to_bgcout_int(cname,val)
  INTEGER,INTENT(in) ::val
  CHARACTER( LEN = * ):: cname
  CHARACTER(LEN=max_char_length) :: cpara_name, cpara_val

  cpara_name=cname
  write(cpara_val, '(i0)') val
  CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc)
END SUBROUTINE


 END MODULE
