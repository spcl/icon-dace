!
! mo_art_photolysis
! This module provides the subroutines for photolysis rate calculation
! based on CloudJ
! Prather, M. J.: Photolysis rates in correlated overlapping cloud fields:
! Cloud-J 7.3c, Geosci. Model Dev., 8, 2587-2595, doi:10.5194/gmd-8-2587-2015, 2015
!
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
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE mo_art_photolysis
    

  USE mo_kind,                    ONLY: wp
  USE mo_physical_constants,      ONLY: ppmv2gg,rhoh2o, amd, avo, argas, grav, p0sl_bg, rd
  USE mo_math_constants,          ONLY: rad2deg, pi, pi2
  USE mo_util_phys,               ONLY: rel_hum
  
  USE mo_run_config,              ONLY: iqc, iqv, iqi
  USE mo_art_config,              ONLY: art_config
  
  USE CLD_SUB_MOD,                ONLY: CLOUD_JX
  USE FJX_CMN_MOD,                ONLY: JXL1_, AN_, JVN_, JFACTA, JIND, NRATJ, S_
  USE mo_atm_phy_nwp_config,      ONLY: atm_phy_nwp_config
   
  USE mo_art_data,                ONLY: p_art_data
  USE mo_art_chem_data,           ONLY: nphot,          &
                                    &   t_art_photolysis
  USE mo_art_atmo_data,           ONLY: t_art_atmo
  USE mo_art_wrapper_routines,    ONLY: art_get_indices_c
  
  USE mo_fjx_lyman,               ONLY: photlyman
#ifdef __ART_GPL
  USE messy_mecca_kpp_parameters, ONLY: ind_O3
#endif


  IMPLICIT NONE


  PRIVATE
  
  PUBLIC   ::             &
      &   art_photolysis  
         
CONTAINS

!!
!!-----------------------------------------------------------------------------
!!



SUBROUTINE art_photolysis(jg, p_tracer_now, vmr2Nconc)

  !<
  ! SUBROUTINE art_photolysis                  
  ! This subroutine calculates the photolysis rates using CloudJ
  ! Part of Module: mo_art_photolysis
  ! Author: Jennifer Schroeter, KIT
  ! Initial Release: 2014-10-15
  ! Modifications:
  !>

  INTEGER, INTENT(in)  ::  &
    &  jg                !< patch id on which computation is performed
  REAL(wp), INTENT(INOUT) ::    &
    &  p_tracer_now(:,:,:,:)  !< tracer mixing ratios (specific concentrations)
                              !< at current time level n (before transport)
                              !< [kg/kg]
                              !< dim: (art_atmo%nproma,art_atmo%nlev,art_atmo%nblks_c,ntracer)
  REAL(wp), INTENT(in) :: &
    &  vmr2Nconc(:,:,:)       !< volume mixing ratio (mol/mol) to number concentration (#/cm3)

  !Local variables
  TYPE(t_art_atmo),POINTER    :: &
    &  art_atmo                     !< Pointer to ART atmo fields 
  INTEGER, POINTER :: &
    &  mapping_indices_kpp(:)
  TYPE(t_art_photolysis),POINTER    :: &
    &  art_photo                   !< Pointer to ART photolysis fields

  REAL(wp)              ::  fjx_zfact
  REAL(wp)              ::  fjx_zkap_cont
  REAL(wp)              ::  fjx_zkap_mrtm 

  
  !===================================================================================

  !--- local fastj parameters as profiles
  INTEGER                                   :: L_, L1_
  REAL(wp)                                  :: U0, SZA, REFLB, SOLF, CLDCOR
  REAL(wp)                                  :: FG0
  LOGICAL                                   :: LPRTJ
  


  INTEGER                                   :: CLDFLAG, NRANDO, IRAN, LNRG
  INTEGER                                   :: NICA, JCOUNT

  !--- other locals
  REAL(wp) :: PMID, WLC, ICWC
  REAL(wp) :: H_scale, ZDEL, PDEL,  F1
  INTEGER :: k, L, ANU, NJXU

  
  REAL(wp), PARAMETER :: M_air_si=amd/1000._wp
  
  INTEGER :: J

  LOGICAL :: lart_lyman
  
  !==================================================================================
 
  INTEGER ::                     &
     &    jc, jk, jb,            &                   !< loop indizes
     &    i_startidx, i_endidx
  REAL (wp), PARAMETER ::     &
     &     ccwmin    = 1.e-7_wp


  ! ----------------------------------
  ! --- Start routine
  ! ----------------------------------

  art_atmo  => p_art_data(jg)%atmo
  art_photo => p_art_data(jg)%chem%photo
#ifdef __ART_GPL
  IF (ALLOCATED(p_art_data(jg)%chem%mecicon%utils%mapping_indices_kpp)) THEN
    mapping_indices_kpp => p_art_data(jg)%chem%mecicon%utils%mapping_indices_kpp
  ELSE
    mapping_indices_kpp => NULL()
  ENDIF
#endif
  
  lart_lyman = .FALSE.


  ! ----------------------------------
  ! --- Initialisation
  ! ----------------------------------

  L_     = art_atmo%nlev-1
  L1_    = L_+1

  art_photo%PPP=0.0
  art_photo%ZZZ=0.0
  art_photo%DDD=0.0
  art_photo%TTT=0.0
  art_photo%RRR=0.0
  art_photo%OOO=0.0
  art_photo%LWP=0.0
  art_photo%IWP=0.0
  art_photo%REFFL=0.0
  art_photo%REFFI=0.0
  art_photo%CLF=0.0
  art_photo%CLDIW=0.0
  art_photo%AERSP=0.0
  art_photo%NDXAER=0

  art_photo%rate(:,:,:,:) = 0.0_wp
    


!  ! ----------------------------------
!  ! --- Variables needed for CloudJ
!  ! ----------------------------------

  fjx_zfact = 1.0e6_wp*(3.0e-9_wp/(4.0_wp*pi*rhoh2o))**(1.0_wp/3.0_wp)
  fjx_zkap_cont = 1.143_wp
  fjx_zkap_mrtm = 1.077_wp
  
  
  art_photo%input_fastjx_lwc(:,:,:) = 0.0
  art_photo%input_fastjx_iwc(:,:,:) = 0.0
  art_photo%input_fastjx_reice(:,:,:) = 0.0
  art_photo%input_fastjx_reliq(:,:,:) = 0.0

!  ! ----------------------------------
!  ! --- 2d information
!  ! ----------------------------------


!  ! Convert land fraction into 1 / 0 values
!  ! Get fraction of glacier
          
  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

!NEC$ ivdep
    DO jc=i_startidx,i_endidx
      art_photo%fjx_zland(jc,jb)   =  art_atmo%fr_land(jc,jb)

      art_photo%fjx_zglac(jc,jb) = art_atmo%fr_glac(jc,jb)
    ENDDO
  ENDDO


!  ! Preparation for cloud radius calculation adapted from mo_newcld_optics

  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

!NEC$ ivdep
    DO jc=i_startidx,i_endidx

      art_photo%fjx_zkap(jc,jb) = fjx_zkap_cont * ( art_photo%fjx_zland(jc,jb) - art_photo%fjx_zglac(jc,jb))           &
               +  fjx_zkap_mrtm * (1.0_wp - art_photo%fjx_zland(jc,jb) + art_photo%fjx_zglac(jc,jb))

    ENDDO
  ENDDO

!  ! ----------------------------------
!  ! --- 3d information
!  ! ----------------------------------

  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

    DO jk=1,art_atmo%nlev
!NEC$ ivdep
      DO jc=i_startidx,i_endidx
        ! get values for temperature, etc. for calculation of relative humidity

        art_photo%input_fastjx_rh(jc,jk,jb)   = rel_hum(art_atmo%temp(jc,jk,jb),       &
                   &                          MAX(p_tracer_now(jc,jk,jb,iqv),0._wp),   &
                   &                          art_atmo%exner(jc,jk,jb)) / 100.

      ENDDO
    ENDDO
  ENDDO

  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

    DO jk=1,art_atmo%nlev
!NEC$ ivdep
      DO jc=i_startidx,i_endidx
        art_photo%input_fastjx_clf(jc, jk, jb) = art_atmo%clc(jc,jk,jb)
        

        art_photo%input_fastjx_iwc(jc,jk,jb) = MAX(art_atmo%tot_cld(jc,jk,jb,iqi),0._wp)   &
                &                    * art_atmo%rho(jc,jk,jb) * 1000.
        art_photo%input_fastjx_lwc(jc,jk,jb) = MAX(art_atmo%tot_cld(jc,jk,jb,iqc),0._wp)   &
                &                    * art_atmo%rho(jc,jk,jb) * 1000.
      

        art_photo%input_fastjx_lwp(jc,jk,jb)  = MAX(art_atmo%tot_cld(jc,jk,jb,iqc),0._wp)   &
                   &       * art_atmo%rho(jc,jk,jb) * 1000.      &
                   &       * art_atmo%dz(jc,jk,jb) 
        art_photo%input_fastjx_iwp(jc,jk,jb)  = MAX(art_atmo%tot_cld(jc,jk,jb,iqi),0._wp)   &
                    &      * art_atmo%rho(jc,jk,jb) * 1000.      &
                    &      * art_atmo%dz(jc,jk,jb) 

        art_photo%fjx_cdnc(jc,jk,jb) = art_atmo%acdnc(jc,jk,jb)
      ENDDO
    ENDDO
  ENDDO

  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

  ENDDO
! -------------------------------------------------------------------------------------------------

!
!  ! ----------------------------------
!  ! --- Main CloudJ call
!  ! ---------------------------------


  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

!NEC$ ivdep
    DO jc=i_startidx,i_endidx

     
      U0 = art_atmo%sza(jc,jb)


      IF (U0 .gt. -0.13) THEN
        SZA = rad2deg*acos( U0 )

        REFLB = art_atmo%albedo(jc,jb)

        SOLF  = 1.d0-(0.034d0*cos(dble(39-186)*pi2/365.d0))
        FG0   = 0.0_wp
        LPRTJ = .false.


        art_photo%PPP(L1_)   = exp(-1.) * art_atmo%pres(jc,1, jb) / 100.

        DO k=1,L_-1
           art_photo%PPP(L_-k+1) =( art_atmo%pres(jc,k+1,jb) )/ 100. ! -> hP
        ENDDO

        art_photo%PPP(1) = art_photo%PPP(2) * art_photo%PPP(2) / art_photo%PPP(3)

        DO K=1,L1_
           art_photo%DDD(K) = ( art_photo%PPP(K) - art_photo%PPP(K+1) )*100._wp/grav  & ! -> kg/m2
           !DDD(K) = ( PPP(K) - PPP(K+1) )/g  & ! -> kg/m2
                &    /M_air_si                    & ! -> -> mol/m2 (note SI)
                &    *avo                         & ! -> molec,
                &    /1.e4_wp                          ! -> cm-2
        ENDDO
  

        IF (art_config(jg)%lart_chemtracer) THEN
          DO K=1,L_
            art_photo%OOO(K) =  (art_atmo%o3_field_icon(jc,L_-K+1,jb) / ppmv2gg)  &
                      &          * 1e-6_wp * art_photo%DDD(K)
          ENDDO
        END IF

        IF (art_config(jg)%lart_mecca) THEN
#ifdef __ART_GPL
          IF (ind_O3 /= 0) THEN
            DO K=1,L_
              art_photo%OOO(K) =  (p_tracer_now(jc,L_-K+1,jb, mapping_indices_kpp(ind_O3))  &
                      &   / vmr2Nconc(jc,L_-K+1,jb)) * art_photo%DDD(K)
            END DO
          ELSE
            DO K=1,L_
              art_photo%OOO(K) =  (art_atmo%o3_field_icon(jc,L_-K+1,jb) / ppmv2gg) &
                                &  * 1e-6_wp * art_photo%DDD(K)
            END DO
          END IF
#endif
        END IF

        DO K=1,L_
          art_photo%TTT(K) = art_atmo%temp(jc,L_-K+1,jb)
          art_photo%RRR(K) = art_photo%input_fastjx_rh(jc,L_-K+1,jb)
          art_photo%CLF(K) =  art_photo%input_fastjx_clf(jc,L_-K+1,jb)
        ENDDO

        art_photo%TTT(L1_) = art_photo%TTT(L_)
        art_photo%RRR(L1_) = art_photo%RRR(L_)
        art_photo%OOO(L1_) = art_photo%OOO(L_)
        art_photo%CLF(L1_) = art_photo%CLF(L_)

        DO L=1,L1_
          IF ( art_photo%CLF(L) < 0.005_wp ) then
             !  CLDIW is an integer flag: 1 = water cld, 2 = ice cloud, 3 = both
             art_photo%CLDIW(L) = 0
             art_photo%CLF(L)   = 0._wp
             art_photo%LWP(L)   = 0._wp
             art_photo%IWP(L)   = 0._wp
          ELSE
             ! in PHOTO_JX rescaled
             art_photo%LWP(L) =   art_photo%input_fastjx_lwp(jc,L_-L+1,jb) / art_photo%CLF(L)
             art_photo%IWP(L) =   art_photo%input_fastjx_iwp(jc,L_-L+1,jb) / art_photo%CLF(L)
          ENDIF
        ENDDO

        art_photo%ZZZ(1)  = 16.d5*log10(p0sl_bg/art_photo%PPP(1))     ! zzz in cm

        DO L = 1,L_
          H_scale  = Rd * art_photo%TTT(L) / grav * 100.         ! in cm
          ZDEL     = -( log(art_photo%PPP(L+1)/art_photo%PPP(L)) * H_scale )
          PDEL     = art_photo%PPP(L) - art_photo%PPP(L+1)
          PMID     = ( art_photo%PPP(L) + art_photo%PPP(L+1) ) / 2
          art_photo%ZZZ(L+1) = art_photo%ZZZ(L) + ZDEL
          WLC      = art_photo%LWP(L) * grav / 100000. / PDEL

          IF (WLC > 1.d-11) THEN                ! note: cldj.f90 uses 1e-11 and 1e-12
             F1       = 0.005d0 * (PMID - 610.d0)
             F1       = min(1.d0, max(0.d0, F1))
             !art_photo%REFFL(L) = 9.6d0*F1 + 12.68d0*(1.d0-F1)
             art_photo%CLDIW(L) = 1
          ELSE
             art_photo%LWP(L)   = 0.d0
             art_photo%REFFL(L) = 0.d0
          ENDIF

          IF ( art_photo%IWP(L)*grav/10000./PDEL > 1.d-11 ) THEN
             !ICWC     = art_photo%IWP(L) / ZDEL          ! g/m3
             !art_photo%REFFI(L) = 164.d0 * (ICWC**0.23d0)
             art_photo%CLDIW(L) = art_photo%CLDIW(L) + 2
          ENDIF
        ENDDO

        art_photo%ZZZ(L1_+1) = art_photo%ZZZ(L1_) + H_scale ! use scale height from highest layer = last loop element
        
        CLDCOR  = 0._wp
        FG0     = 0._wp
        LPRTJ   = .FALSE.
        art_photo%AERSP   = 0.
        art_photo%NDXAER  = 0
        CLDFLAG = 2
        NRANDO  = 0
        IRAN    = 0
        LNRG    = 0
        NICA    = 0
        ANU     = AN_
        NJXU    = JVN_
        JCOUNT  = 0

        CALL CLOUD_JX ( U0       & ! U0 = cos (SZA)
                      , SZA      & ! SZA = solar zenith angle in degree
                      , REFLB    & ! REFLB = Lambertian reflectivity at the Lower Boundary
                      , SOLF     & ! SOLF = solar flux factor for sun-earth distance
                      , FG0      & ! FG0 = scale for asymmetry factor (CLDFLAG=3
                      , LPRTJ    & ! LPRTJ = .true. = turn on internal print
                      , art_photo%PPP      & ! P = top edge press (hPa)), #tr PPP(1) surface
                      , art_photo%ZZZ      & ! Z = edge alt (cm)))
                      , art_photo%TTT      & ! T = layer temp (K
                      , art_photo%DDD      & ! D = layer dens-path (# molec /cm2)
                      , art_photo%RRR      & ! R = layer rel.hum.(fraction)
                      , art_photo%OOO      & ! O = layer O3 path (# O3 /cm2
                      , art_photo%LWP      & ! LWP/IWP = Liquid/Ice water path (g/m2)
                      , art_photo%IWP      & ! (probably meant for cloudy part only as scaled for CLDFLG < 4)
                      , art_photo%REFFL    & ! REFFL/REFFI = R-effective(microns) in liquid/ice cloud
                      , art_photo%REFFI    &
                      , art_photo%CLF      & ! CLF = cloud fraction (0.0 to 1.0)
                      , CLDCOR   &
                      , art_photo%CLDIW    & ! CLDIW = denoting cloud in layer: 1=water, 2=ice, 3=both
                      , art_photo%AERSP    & ! AERSP = aerosol path (g/m2) & NDXAER = aerosol index type
                      , art_photo%NDXAER   & ! aerosols with up to AN_ different types in an ICA layer
                      , L1_      & ! L1_ = dim of profile variables, L_+1 top (non CTM) layer
                      , ANU      & ! AN_ = parameter, dim of number of aerosols being passed
                      , art_photo%VALJXX   & ! VALJXX = J-values from CLOUD_JX & PHOTO_JX
                      , NJXU     &
                      , CLDFLAG  & ! CLDFLAG = integer index for type of cloud overlap
                      , NRANDO   & ! NRANDO = number of random ICAs to do J's for (CLDFLAG=4)
                      , IRAN     & ! IRAN = starting index for random number selection
                      , LNRG     & ! LNRG = flag max-ran overlap groups: 0, 3, 6 
                      , NICA     & ! NICA = total number of ICAs
                      , JCOUNT   & ! 
                      , art_photo%SKPERD  &
                      )

        DO J = 1,S_+2
          DO K=1,L_
            art_photo%heating_rates(jc,K,jb,J) = art_photo%SKPERD(J,K)
          ENDDO
        ENDDO

        DO J = 1,NRATJ
          IF (JIND(J).gt.0) then
            art_photo%rate(jc,1:art_atmo%nlev-1,jb,J) = art_photo%VALJXX(:,JIND(J))  &
                                 &                      * JFACTA(J)
          ELSE
             art_photo%rate(jc,:,jb,J) = 0.0_wp
          ENDIF
        ENDDO

        IF (lart_lyman) THEN
          CALL photlyman(art_photo%TTT,              &
                  &      art_photo%DDD,              &
                  &      SZA,                        &
                  &      art_photo%rate(jc,:, jb,:), &
                  &      art_atmo%nlev)
        ENDIF
        !  Photolysis rates are upside down
      ELSE ! of if (U0 .gt. -0.13) 
        art_photo%rate(jc,:,jb,:)= 0.0_wp
      ENDIF

    ENDDO
  ENDDO

  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

    DO jk = 1,art_atmo%nlev-1
  
!NEC$ ivdep
      DO jc = i_startidx,i_endidx
 
        art_photo%output_fastjx_photo_new(jc,jk,jb,:) = art_photo%rate(jc,art_atmo%nlev-jk,jb,:)

      ENDDO
    ENDDO
   
    art_photo%rate(:,:,jb,:)    =  art_photo%output_fastjx_photo_new(:,:,jb,:)
    art_photo%rate(:,art_atmo%nlev,jb,:) =  art_photo%rate(:,art_atmo%nlev-1,jb,:)
  ENDDO


END SUBROUTINE art_photolysis

END MODULE mo_art_photolysis




