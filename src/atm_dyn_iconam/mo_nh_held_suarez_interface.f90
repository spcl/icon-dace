!
! This module contains the interface between ICONAM dynamics and Held-Suarez forcing
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nh_held_suarez_interface

  USE mo_kind,                  ONLY: wp, vp
  USE mo_parallel_config,       ONLY: nproma
  USE mo_model_domain,          ONLY: t_patch
  USE mo_nonhydro_types,        ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_intp_data_strc,        ONLY: t_int_state
  USE mo_intp,                  ONLY: cells2edges_scalar
  USE mo_intp_rbf,              ONLY: rbf_vec_interpol_cell
  USE mo_impl_constants_grf,    ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_loopindices,           ONLY: get_indices_c, get_indices_e
  USE mo_nh_hs_test,            ONLY: held_suarez_forcing_temp, &
                                      & held_suarez_forcing_vn
  USE mo_nh_diagnose_pres_temp, ONLY: diagnose_pres_temp
  USE mo_physical_constants,    ONLY: rd_o_cpd
  USE mo_nh_testcases_nml,      ONLY: lhs_fric_heat
  USE mo_timer,                 ONLY: ltimer, timer_start, timer_stop, timer_held_suarez_intr
  USE mo_fortran_tools,         ONLY: init

  IMPLICIT NONE

  PRIVATE


  PUBLIC  :: held_suarez_nh_interface

CONTAINS

  !>
  !! SUBROUTINE held_suarez_nh_interface -- interface between ICONAN dynamics
  !! and Held-Suarez forcing
  !!
  SUBROUTINE held_suarez_nh_interface (p_nh_prog,p_patch,p_int_state,p_metrics,p_nh_diag)

    TYPE(t_patch),TARGET,INTENT(in):: p_patch    !< single patch
    TYPE(t_int_state),INTENT(in)  :: p_int_state!< single interpolation state
    TYPE(t_nh_metrics),INTENT(in) :: p_metrics  !< single metrics state
    TYPE(t_nh_prog), INTENT(inout)   :: p_nh_prog  !< single nh prognostic state
    TYPE(t_nh_diag), TARGET, INTENT(inout):: p_nh_diag  !< single nh diagnostic state

    ! Local scalar

    INTEGER :: jk, jb, jbs, is, ie
    INTEGER :: nblks_c, nblks_e
    INTEGER :: nlev              !< number of full levels

    ! Local arrays

    REAL(wp) :: zlat(nproma)      ! latitude

    REAL(wp) ::   & !< pressure @ cells
      &  zsigma_mc(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp) ::   & !< pressure @ edges
      &  zsigma_me(nproma,p_patch%nlev,p_patch%nblks_e)

    REAL(wp) ::   & !< tendency of temp due to HS forcing
      &  ddt_temp (nproma,p_patch%nlev,p_patch%nblks_c)

    REAL(wp) ::   & !< kinetic energy @ cells
      &  z_ekin(nproma,p_patch%nlev)

    REAL(wp) ::   & !< forcing on vn
      &  zddt_vn(nproma,p_patch%nlev)

    INTEGER :: i

!!$    CHARACTER(len=*), PARAMETER :: routine =  &
!!$                   '(mo_nh_held_suarez_interface) held_suarez_nh_interface:'

    !-------------------------------------------------------------------------
    ! Dimension parameters related to refinement and MPI parallelisation
    IF (ltimer) CALL timer_start(timer_held_suarez_intr)

    nblks_e = p_patch%nblks_e
    nblks_c = p_patch%nblks_c

    ! number of vertical levels
    nlev = p_patch%nlev

    !$ACC DATA CREATE(zsigma_mc, zsigma_me, ddt_temp, z_ekin, zddt_vn, zlat) &
    !$ACC   PRESENT(p_nh_diag, p_int_state, p_metrics, p_nh_prog, p_patch%cells%center)
    
    !-------------------------------------------------------------------------
    ! First the surface pressure, pressure and temperature must be diagnosed

    CALL diagnose_pres_temp ( p_metrics, p_nh_prog,               &
      &                       p_nh_prog, p_nh_diag,               &
      &                       p_patch,                            &
      &                       opt_calc_temp=.TRUE.,               &
      &                       opt_calc_pres=.TRUE. )

    IF (lhs_fric_heat) CALL rbf_vec_interpol_cell(p_nh_prog%vn,p_patch,p_int_state,p_nh_diag%u,p_nh_diag%v)

    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    ! Newtonian cooling (and optionallly frictional heating due to Rayleigh friction)
    !-------------------------------------------------------------------------

    jbs = p_patch%cells%start_blk( grf_bdywidth_c+1,1 )
!$OMP PARALLEL
    CALL init(ddt_temp(:,:,:))
!$OMP BARRIER
!$OMP DO PRIVATE(jb,is,ie,jk,z_ekin,zlat) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs,nblks_c
       CALL get_indices_c( p_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1 )

       !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
       !$ACC LOOP GANG VECTOR COLLAPSE(2)
       DO jk=1,nlev
          DO i=is,ie
             zsigma_mc(i,jk,jb) = p_nh_diag%pres(i,jk,jb)/p_nh_diag%pres_sfc(i,jb)
          ENDDO
       ENDDO
       !$ACC END PARALLEL

       IF (lhs_fric_heat) THEN
         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
         !$ACC LOOP GANG VECTOR COLLAPSE(2)
         DO jk=1,nlev
            DO i=is,ie
               z_ekin(i,jk) = 0.5_wp*(p_nh_diag%u(i,jk,jb)**2+p_nh_diag%v(i,jk,jb)**2)
            ENDDO
         ENDDO
         !$ACC END PARALLEL
       ELSE
         !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
         z_ekin(:,:) = 0._wp
         !$ACC END KERNELS
       ENDIF

! WS: DEFAULT(NONE): Data clause required with default(none): p_patch   Why?
       !$ACC PARALLEL ASYNC(1)
       !$ACC LOOP GANG VECTOR
       DO i=is,ie
          zlat(i) = p_patch%cells%center(i,jb)%lat
       ENDDO
       !$ACC END PARALLEL
       ! last 2 inputs in case of additional computation of frictional heating
       CALL held_suarez_forcing_temp( p_nh_diag%temp(:,:,jb),     &! in
                                    & p_nh_diag%pres(:,:,jb),     &! in
                                    & zsigma_mc(:,:,jb), zlat(:), &! in
                                    & nlev, nproma, is, ie,       &! in
                                    & ddt_temp(:,:,jb),           &! out
                                    & z_ekin(:,:), lhs_fric_heat)  ! optional in

       ! the tendency in temp must be transfromed to a tendency in the exner function
       ! For this it is assumed that the density is constant
       
       !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
       !$ACC LOOP GANG VECTOR COLLAPSE(2)
       DO jk=1,nlev
          DO i=is,ie
             p_nh_diag%ddt_exner_phy(i,jk,jb)=rd_o_cpd/p_nh_prog%theta_v(i,jk,jb)*ddt_temp(i,jk,jb)
          ENDDO
       ENDDO
       !$ACC END PARALLEL
       
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !-------------------------------------------------------------------------
    ! Rayleigh friction
    !-------------------------------------------------------------------------
    ! First interpolate sigma values from cells to edges

    CALL cells2edges_scalar( zsigma_mc,                    &! in
                           & p_patch, p_int_state%c_lin_e, &! in
                           & zsigma_me, lacc=.TRUE.)        ! out
    ! Now compute the velocity tendency due to friction

    jbs = p_patch%edges%start_blk( grf_bdywidth_e+1,1 )
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,jk,zddt_vn) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs,nblks_e
       CALL get_indices_e( p_patch, jb,jbs,nblks_e, is,ie, grf_bdywidth_e+1 )

       CALL held_suarez_forcing_vn( p_nh_prog%vn(:,:,jb),  &! in
                                  & zsigma_me(:,:,jb),     &! in
                                  & nlev, nproma, is, ie,  &! in
                                  & zddt_vn )               ! inout
       !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
       !$ACC LOOP GANG VECTOR COLLAPSE(2)
       DO jk=1,nlev
          DO i=is,ie
             p_nh_diag%ddt_vn_phy(i,jk,jb) = zddt_vn(i,jk)
          ENDDO
       ENDDO
       !$ACC END PARALLEL
    ENDDO
    !$ACC WAIT(1)
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC END DATA

    !--------------------------------------------------------------------------
    IF (ltimer) CALL timer_stop(timer_held_suarez_intr)

  END SUBROUTINE held_suarez_nh_interface

END MODULE mo_nh_held_suarez_interface

