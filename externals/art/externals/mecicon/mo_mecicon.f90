!
! mo_mecicon
! This module wraps the box model MECCA (Code based on
! Sander, Rolf, et al. "Technical note: The new comprehensive
!        atmospheric chemistry module MECCA."
!        Atmospheric Chemistry and Physics 5.2 (2005): 445-450.
!Sander,R. et al.: The atmospheric chemistry box model CAABA/MECCA-3.0,
!       Geosci. Model Dev., 4, 373-380, doi:10.5194/gmd-4-373-2011, 2011.
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
! SPDX-License-Identifier: GPL-3.0-only  
! ---------------------------------------------------------------

MODULE mo_mecicon

  USE mo_kind,                          ONLY: wp
  USE mo_art_config,                    ONLY: art_config
  USE mo_art_psc_types,                 ONLY: t_art_psc
  USE messy_mecca_kpp_global
  USE messy_mecca_kpp_function,         ONLY: messy_A => A

  ! ART
  USE mo_art_impl_constants,            ONLY: IART_QV
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_chem_data,                 ONLY: t_art_photolysis
  USE mo_art_mecicon_data,              ONLY: t_art_mecicon
  USE messy_mecca_kpp_parameters
  USE messy_cmn_photol_mem

  

  IMPLICIT NONE 

  PUBLIC :: mecicon_call

  PRIVATE

CONTAINS





SUBROUTINE mecicon_call(jc, jk, jb, jg, p_dtime,pres_icon, &
                               temp_icon, p_tracer_now)
!<
! SUBROUTINE art_mecicon_call
! This subroutine calls the subroutine for
! the 1d box of mecca (CALL mecca_physc)
! Photolysis rates of CloudJ are mappend in the
! respective order 
! Part of Module: mo_art_mecicon
! Author: Jennifer Schroeter, KIT
! Initial Release: 2015-06-25                
! Modifications:
!>
  IMPLICIT NONE


  REAL(wp), INTENT(INOUT)      :: p_tracer_now(:,:,:,:)   !< tracer mixing ratios (specific conc.)
                                                          !< at current time level n 
                                                          !< (before transport) [kg/kg]
                                                          !< dim: (nproma,nlev,nblks_c,ntracer)
  REAL(wp), INTENT(IN)         :: p_dtime                 !<time step
  INTEGER, INTENT(IN)          :: jc,jk,jb

  INTEGER, INTENT(in)          :: jg                      !< patch on which computation is performed
  REAL(wp), intent(in)         :: pres_icon, temp_icon
  ! local variables
  INTEGER ::  &
    &  ind_tr      !< index of the tracer
  INTEGER ::  &
    &  ihs_reac    !< reaction rate index
  TYPE(t_art_photolysis),POINTER  :: &
    &  art_photo   !< Pointer to ART photolysis fields 
  TYPE(t_art_mecicon), POINTER :: &
    &  art_mecicon !< Pointer to ART MECCA fields 
  TYPE(t_art_psc), POINTER  :: &
    &  art_PSC     !< Pointer to ART PSC fields 
  INTEGER, POINTER :: &
    &  mapping_indices_kpp(:)


  !! Initialisation, getting meteorological variables temp and pres
  !! from ICON and calculate cair

  art_photo       => p_art_data(jg)%chem%photo
  art_mecicon     => p_art_data(jg)%chem%mecicon
  art_PSC         => p_art_data(jg)%chem%PSC_meta
  mapping_indices_kpp => p_art_data(jg)%chem%mecicon%utils%mapping_indices_kpp

  temp            = temp_icon
  press           = pres_icon
  cair            = p_art_data(jg)%chem%vmr2Nconc(jc,jk,jb)

  !! Loop for 3d (p_tracer_now) onto 1d (c) concentration

  DO ind_tr = 1,NSPEC
    c(ind_tr) = p_tracer_now(jc,jk,jb, mapping_indices_kpp(ind_tr))
  ENDDO



  !! This part remains from the original box model
  !! Not valid anymore, since we are calculating 
  !! chemical reation rates with the same time step
  !! as the main ICON timestep p_dtime

  art_mecicon%kpp_control%time_step = p_dtime

  !! Photolysis rates of CloudJ (art_photo%rate) in 3d
  !! are mappend on 1d fields (jx_fjx)

  art_mecicon%photo%jx_fjx(ip_NO)       = art_photo%rate(jc,jk,jb,ip_NO)
  art_mecicon%photo%jx_fjx(ip_O2)       = art_photo%rate(jc,jk,jb,ip_O2)
  art_mecicon%photo%jx_fjx(ip_O3P)      = art_photo%rate(jc,jk,jb,3) - art_photo%rate(jc,jk,jb,4)
  art_mecicon%photo%jx_fjx(ip_O1D)      = art_photo%rate(jc,jk,jb,ip_O1D)
  art_mecicon%photo%jx_fjx(ip_CHOH)     = art_photo%rate(jc,jk,jb,ip_CHOH)
  art_mecicon%photo%jx_fjx(ip_COH2)     = art_photo%rate(jc,jk,jb,ip_COH2)
  art_mecicon%photo%jx_fjx(ip_H2O2)     = art_photo%rate(jc,jk,jb,ip_H2O2)
  art_mecicon%photo%jx_fjx(ip_CH3OOH)   = art_photo%rate(jc,jk,jb,ip_CH3OOH)
  art_mecicon%photo%jx_fjx(ip_NO2)      = art_photo%rate(jc,jk,jb,ip_NO2)
  art_mecicon%photo%jx_fjx(ip_NOO2)     = art_photo%rate(jc,jk,jb,ip_NOO2)
  art_mecicon%photo%jx_fjx(ip_NO2O)     = art_photo%rate(jc,jk,jb,ip_NO2O)
  art_mecicon%photo%jx_fjx(ip_N2O5)     = art_photo%rate(jc,jk,jb,ip_N2O5)
  art_mecicon%photo%jx_fjx(ip_HONO)     = art_photo%rate(jc,jk,jb,ip_HONO)
  art_mecicon%photo%jx_fjx(ip_HNO3)     = art_photo%rate(jc,jk,jb,ip_HNO3)
  art_mecicon%photo%jx_fjx(ip_HNO4)     = art_photo%rate(jc,jk,jb,ip_HNO4)
  art_mecicon%photo%jx_fjx(ip_ClNO3)    = art_photo%rate(jc,jk,jb,ip_ClNO3)
  art_mecicon%photo%jx_fjx(ip_ClONO2)   = art_photo%rate(jc,jk,jb,ip_ClONO2)
  art_mecicon%photo%jx_fjx(ip_Cl2)      = art_photo%rate(jc,jk,jb,ip_Cl2)
  art_mecicon%photo%jx_fjx(ip_HOCl)     = art_photo%rate(jc,jk,jb,ip_HOCl)
  art_mecicon%photo%jx_fjx(ip_OClO)     = art_photo%rate(jc,jk,jb,ip_OClO)
  art_mecicon%photo%jx_fjx(ip_Cl2O2)    = art_photo%rate(jc,jk,jb,ip_Cl2O2)
  art_mecicon%photo%jx_fjx(ip_BrO)      = art_photo%rate(jc,jk,jb,ip_BrO)
  art_mecicon%photo%jx_fjx(ip_BrNO3)    = art_photo%rate(jc,jk,jb,ip_BrNO3)
  art_mecicon%photo%jx_fjx(ip_HOBr)     = art_photo%rate(jc,jk,jb,ip_HOBr)
  art_mecicon%photo%jx_fjx(ip_BrCl)     = art_photo%rate(jc,jk,jb,ip_BrCl)
  art_mecicon%photo%jx_fjx(ip_OCS)      = art_photo%rate(jc,jk,jb,ip_OCS)
  art_mecicon%photo%jx_fjx(ip_SO2)      = art_photo%rate(jc,jk,jb,ip_SO2)
  art_mecicon%photo%jx_fjx(ip_N2O)      = art_photo%rate(jc,jk,jb,ip_N2O)
  art_mecicon%photo%jx_fjx(ip_CFCl3)    = art_photo%rate(jc,jk,jb,ip_CFCl3)
  art_mecicon%photo%jx_fjx(ip_CF2Cl2)   = art_photo%rate(jc,jk,jb,ip_CF2Cl2)
  art_mecicon%photo%jx_fjx(ip_CCl4)     = art_photo%rate(jc,jk,jb,ip_CCl4)
  art_mecicon%photo%jx_fjx(ip_CH3Cl)    = art_photo%rate(jc,jk,jb,ip_CH3Cl)
  art_mecicon%photo%jx_fjx(ip_CH3CCl3)  = art_photo%rate(jc,jk,jb,ip_CH3CCl3)
  art_mecicon%photo%jx_fjx(ip_CH3Br)    = art_photo%rate(jc,jk,jb,ip_CH3Br)
  art_mecicon%photo%jx_fjx(ip_CF2ClBr)  = art_photo%rate(jc,jk,jb,ip_CF2ClBr)
  art_mecicon%photo%jx_fjx(ip_CH2Br2)   = art_photo%rate(jc,jk,jb,ip_CH2Br2)
  art_mecicon%photo%jx_fjx(ip_CHBr3)    = art_photo%rate(jc,jk,jb,ip_CHBr3)
  art_mecicon%photo%jx_fjx(ip_CH3I)     = art_photo%rate(jc,jk,jb,ip_CH3I)
  art_mecicon%photo%jx_fjx(ip_PAN)      = art_photo%rate(jc,jk,jb,ip_PAN)
  art_mecicon%photo%jx_fjx(ip_CH3NO3)   = art_photo%rate(jc,jk,jb,ip_MGLYOX)
  art_mecicon%photo%jx_fjx(ip_CH3CHO)   = art_photo%rate(jc,jk,jb,ip_CH3CHO)
  art_mecicon%photo%jx_fjx(ip_MVK)      = art_photo%rate(jc,jk,jb,55) + art_photo%rate(jc,jk,jb,56)
  art_mecicon%photo%jx_fjx(ip_CF3Br)    = art_photo%rate(jc,jk,jb,ip_CF3Br)
  art_mecicon%photo%jx_fjx(ip_MACR)     = art_photo%rate(jc,jk,jb,ip_MACR)
  art_mecicon%photo%jx_fjx(ip_HOCH2CHO) = art_photo%rate(jc,jk,jb,58)+  art_photo%rate(jc,jk,jb,59)
  art_mecicon%photo%jx_fjx(ip_KET)      = art_photo%rate(jc,jk,jb,61)+  art_photo%rate(jc,jk,jb,62)
  art_mecicon%photo%jx_fjx(ip_MGLYOX)   = art_photo%rate(jc,jk,jb,ip_MGLYOX)

  art_mecicon%photo%jx_fjx(ip_GLYOX)    = art_photo%rate(jc,jk,jb,65)         &
                                            &  +  art_photo%rate(jc,jk,jb,66) &
                                            &  +  art_photo%rate(jc,jk,jb,67)

  art_mecicon%photo%jx_fjx(ip_GLYOXb)    = art_photo%rate(jc,jk,jb,65)+  art_photo%rate(jc,jk,jb,66)
  art_mecicon%photo%jx_fjx(ip_GLYOXa)    = art_photo%rate(jc,jk,jb,67) 
  art_mecicon%photo%jx_fjx(ip_CH3COCH3)  = art_photo%rate(jc,jk,jb,68)+  art_photo%rate(jc,jk,jb,69)
  art_mecicon%photo%jx_fjx(ip_H2O)       = art_photo%rate(jc,jk,jb,ip_H2O)
  art_mecicon%photo%jx_fjx(ip_CO2)       = art_photo%rate(jc,jk,jb,ip_CO2)
 
  ! calculate ratelysis rate of Br2 from that of Cl2
  art_mecicon%photo%jx_fjx(ip_Br2)       = 2._wp * art_photo%rate(jc,jk,jb,ip_Cl2)

  ! calculate ratelysis of ClNO2 from ClONO2
  art_mecicon%photo%jx_fjx(ip_ClNO2)      =  art_photo%rate(jc,jk,jb,ip_ClONO2)


  ! Heterogeneous reactions on polar stratospheric clouds
  IF (art_config(jg)%lart_psc) THEN
    DO ihs_reac = 1,IHS_MAX
      khet_St(ihs_reac)  = art_PSC%k_het(jc,jk,jb,ihs_reac)
    END DO
  ELSE
    DO ihs_reac = 1,IHS_MAX
      khet_St(ihs_reac)  = 0._wp
    END DO
  END IF


!  !! Main MECCA call

!  !!  ! time loop:
  CALL mecca_physc(jc,jk,jb,jg)

  DO ind_tr=1,NSPEC
    p_tracer_now(jc,jk,jb, mapping_indices_kpp(ind_tr)) = c(ind_tr)
  ENDDO

  p_tracer_now(jc,jk,jb, mapping_indices_kpp(ind_H2O)) =   &
    &      p_art_data(jg)%chem%water_tracers(jc,jk,jb,IART_QV)

END SUBROUTINE mecicon_call



SUBROUTINE mecca_physc(jc, jk, jb,jg)

!<
! SUBROUTINE mecca_physc
! This subroutine is part of the original code
! of MECCA, see reference 
!    Sander, Rolf, et al. "Technical note: The new comprehensive
!        atmospheric chemistry module MECCA." 
!        Atmospheric Chemistry and Physics 5.2 (2005): 445-450.
!>
  USE messy_mecca_kpp_global,     ONLY: jx
  USE mo_art_data,                ONLY: p_art_data
  USE messy_mecca_kpp_rates,      ONLY: UPDATE_RCONST
  USE messy_mecca_kpp_integrator, ONLY: integrate
  
  INTEGER, INTENT(IN) :: jc, jk, jb, jg
  
  ! local variables
  INTEGER, PARAMETER :: NBL = 1 ! N_block_length

  INTEGER :: ip, i

  TYPE(t_art_mecicon),POINTER    :: &
    &  art_mecicon                     !< Pointer to ART chem fields

  ! declarations for all variables that are transported to the SMCL via fill
  ! subroutines (or via kpp_integrate for "C")
  !REAL(DP) :: jx(IP_MAX) = 0.
  REAL(DP) :: dummy_khet_Tr(IHT_MAX) ! dummy, not needed in box model
  ! (for C, cair, and temp, see caaba_mem.f90)
  ! (see also important notes about temp, press, and cair in gas.eqn!)

  INTEGER, DIMENSION(20) :: istatus_u
  INTEGER                :: ierr_u

  art_mecicon => p_art_data(jg)%chem%mecicon
    
 !   IF (l_ff) THEN
 !     ! at day 4 switch frost flowers on:
 !     IF (model_time >= (model_start_day+4.)*OneDay) THEN
 !       xaer(2) = 1.
 !     ENDIF
 !   ENDIF



  ! The following values are only needed for the mesosphere. Since
  ! CAABA has only one type of temperature, it is also used for
  ! temp_elec and temp_ion:
  !CALL fill_temp_elec(status, SPREAD(temp,1,NBL))
  !CALL fill_temp_ion(status, SPREAD(temp,1,NBL))

  ! dummy values:
  dummy_khet_Tr(:) = 0.
  khet_Tr(:) = dummy_khet_Tr

  DO ip=1, IP_MAX
    JX(ip) = art_mecicon%photo%jx_fjx(ip)
  ENDDO

!  IF (art_mecicon%gen_control%l_aero) THEN
!    CALL aerosol_exchng(jg) ! exchange with fresh aerosol
!    ! first, calculate exchange coefficients for aerosols:
!    CALL mecca_aero_calc_k_ex( &
!      radius, temp, press, l_het, xaer, lwc, c, &         ! in
!      k_exf, k_exb, k_exf_N2O5, k_exf_ClNO3, k_exf_BrNO3) ! out
!    ! next, calculate exchange coefficients for ocean surface:
!    CALL mecca_aero_calc_k_ex_ocean(xaer, radius, temp, art_mecicon%gen_control%zmix, k_exf, k_exb)
!  ENDIF
  

  IF (.NOT.art_mecicon%gen_control%l_skipkpp) THEN
    c(:) = MAX(c(:),1.e-30_DP) ! set negative values to zero

    art_mecicon%kpp_control%istatus = 0



    CALL UPDATE_RCONST

    ! integrate from t=0 to t=dt
    CALL integrate(TIN = 0._dp, &
      & TOUT = art_mecicon%kpp_control%time_step, &
      & ICNTRL_U = art_mecicon%kpp_control%icntrl,& 
      & RCNTRL_U =art_mecicon%kpp_control%rcntrl, &
      & istatus_u=istatus_u,ierr_u=ierr_u)

    ! Return Diagnostic Information

    art_mecicon%kpp_control%ierr = ierr_u
    art_mecicon%kpp_control%xNacc = istatus_u(4)
    art_mecicon%kpp_control%xNrej = istatus_u(5)

    art_mecicon%kpp_control%istatus(1:8) =  art_mecicon%kpp_control%istatus(1:8) + istatus_u(1:8)



    IF (art_config(jg)%lart_diag_out) THEN
      DO i = 1,NREACT
        ! rates of the reactions
        art_mecicon%utils%reac_rates(jc,jk,jb,i) = messy_A(i)
      END DO
    END IF

    !  CALL kpp_integrate(timesteplen,c)  ! main kpp call
  ENDIF


!-------------------------------------------------------------------------

END SUBROUTINE mecca_physc

!***************************************************************************


END MODULE mo_mecicon
