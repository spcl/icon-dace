!
! mo_art_potential_vorticity
! This module provides the subroutines to define variables
! needed for the module mo_define_polar_vortex
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

MODULE mo_art_potential_vorticity
    ! ICON
    USE mo_kind,                   ONLY: wp, vp
    USE mo_physical_constants,     ONLY: g => grav
    USE mo_model_domain,           ONLY: t_patch
    USE mo_physical_constants,     ONLY: earth_radius
    USE mo_math_gradients,         ONLY: grad_fe_cell
    USE mo_impl_constants,         ONLY: min_rlcell_int          
    USE mo_impl_constants_grf,     ONLY: grf_bdywidth_c          
    USE mo_intp_data_strc,         ONLY: p_int_state
    ! ART
    USE mo_art_data,               ONLY: p_art_data
    USE mo_art_atmo_data,          ONLY: t_art_atmo
    USE mo_art_wrapper_routines,   ONLY: art_get_indices_c
      
#ifndef NOMPI
#if !defined (__SUNPRO_F95)
  USE mpi
#endif
#endif
      
  IMPLICIT NONE


  PRIVATE
  
  PUBLIC   :: art_get_potential_vorticity       
         
CONTAINS

SUBROUTINE art_get_potential_vorticity(p_patch,vorticity,dpres,pv)
!<
! SUBROUTINE grt_et_potential_vorticity
! Routine to calculate the potential vortcity
! Part of Module: mo_art_potential_vorticity
! Author: Christian Stassen, KIT
! Initial Release: 2015-03-24
! Remarks: 
!         Jennifer Schroeter (2016-07-13)
!         This routine should be replaced by ICON's internal
!         PV calculation in the near future
!>

  TYPE(t_patch), TARGET, INTENT(IN) ::  &
    &  p_patch                !< patch on which computation is performed
  REAL(wp), INTENT(IN)   ::  &
    &  vorticity(:,:,:)       !< vorticity
  REAL(wp), INTENT(IN)   ::  &
    &  dpres(:,:,:)           !< pressure thickness
  REAL(wp),INTENT(INOUT) ::  &
    &  pv(:,:,:)              !< Potential Vorticity

  !Local variables
  INTEGER ::                     &
     &    jc, jk, jb, jg,        &                   !< loop indizes
     &    i_startidx, i_endidx

  REAL(vp),ALLOCATABLE   :: theta_grad(:,:,:,:)      !< theta_grad(1,:,:,:) zonal

  ! These parameters have to be provided for the gradient calculation which
  ! should be replaced by a wrapper in the future
  INTEGER :: i_rlstart, i_rlend
                                                     !< theta_grad(2,:,:,:) meridional
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  jg = p_patch%id

  art_atmo => p_art_data(jg)%atmo


  ! ----------------------------------
  ! --- Allocation
  ! ----------------------------------

  IF (.NOT. ALLOCATED(theta_grad))  &
    &  ALLOCATE( theta_grad(2,art_atmo%nproma,art_atmo%nlev,art_atmo%nblks) )

  ! ----------------------------------
  ! --- Init. Values
  ! ----------------------------------

  pv(:,:,:)           = 0.
  theta_grad(:,:,:,:) = 0.

  ! ----------------------------------
  ! --- Calculate Potential Vorticity (PV)
  ! ----------------------------------

   i_rlstart   =  grf_bdywidth_c+1
   i_rlend     =  min_rlcell_int
   CALL grad_fe_cell(art_atmo%theta(:,:,:), p_patch, p_int_state(jg), theta_grad(:,:,:,:), &
          & opt_slev=1, opt_elev=art_atmo%nlev, opt_rlstart=i_rlstart, opt_rlend=i_rlend)
              

   DO jb = art_atmo%i_startblk, art_atmo%i_endblk
     CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

     DO jk = 1,art_atmo%nlev-1
       DO jc = i_startidx,i_endidx

         !Potential Vorticity using Eq. 3.1.4 of Andrews "Middle Atmosphere Dynamics"
         pv(jc,jk,jb) = g * ( (p_patch%cells%f_c(jc,jb) + vorticity(jc,jk,jb))                    &
                    & * (art_atmo%theta(jc,jk,jb)-art_atmo%theta(jc,jk+1,jb)) / dpres(jc,jk,jb)   &
                    & - theta_grad(1,jc,jk,jb)*(art_atmo%v(jc,jk,jb)-art_atmo%v(jc,jk+1,jb))      &
                    & / dpres(jc,jk,jb) / (earth_radius*COS(art_atmo%lat(jc,jb)))                 &
                    & + theta_grad(2,jc,jk,jb)*(art_atmo%u(jc,jk,jb)-art_atmo%u(jc,jk+1,jb))      &
                    & / dpres(jc,jk,jb) / earth_radius )
       ENDDO
     ENDDO
   ENDDO

END SUBROUTINE art_get_potential_vorticity

END MODULE mo_art_potential_vorticity
