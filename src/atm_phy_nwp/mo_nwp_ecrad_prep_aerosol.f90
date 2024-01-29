!
! This module prepares aerosol climatologies in a format that can be used by ecRad
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

#if defined __xlC__
@PROCESS SPILL(1058)
#endif
MODULE mo_nwp_ecrad_prep_aerosol

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish
  USE mo_fortran_tools,          ONLY: assert_acc_host_only, assert_acc_device_only,t_ptr_2d
  USE mo_impl_constants,         ONLY: n_camsaermr
#ifdef __ECRAD
  USE mo_ecrad,                  ONLY: t_ecrad_aerosol_type, t_ecrad_conf, t_opt_ptrs
  USE mo_aerosol_util,           ONLY: get_nbands_lw_aerosol, get_nbands_sw_aerosol,   &
                                   &   tegen_scal_factors
#endif

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nwp_ecrad_prep_aerosol'

#ifdef __ECRAD
  PUBLIC :: nwp_ecrad_prep_aerosol

INTERFACE nwp_ecrad_prep_aerosol
  MODULE PROCEDURE nwp_ecrad_prep_aerosol_constant
  MODULE PROCEDURE nwp_ecrad_prep_aerosol_tegen
  MODULE PROCEDURE nwp_ecrad_prep_aerosol_td
  MODULE PROCEDURE nwp_ecrad_prep_aerosol_CAMS
END INTERFACE nwp_ecrad_prep_aerosol
  
  
CONTAINS

  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE nwp_ecrad_prep_aerosol_constant
  !! Prepare aerosol from constant values. If these optional values are not passed,
  !! the corresponding field in ecRad is set to zero.
  !!
  !---------------------------------------------------------------------------------------
  SUBROUTINE nwp_ecrad_prep_aerosol_constant ( slev, nlev, i_startidx, i_endidx,         &
    &                                          ecrad_conf, ecrad_aerosol,                &
    &                                          od_lw, ssa_lw, g_lw, od_sw, ssa_sw, g_sw, &
    &                                          lacc )
    INTEGER, INTENT(in)      :: &
      &  slev, nlev,            & !< Start and end index of vertical loop
      &  i_startidx, i_endidx     !< Start and end index of horizontal loop
    TYPE(t_ecrad_conf),        INTENT(in)    :: &
      &  ecrad_conf                        !< ecRad configuration object
    TYPE(t_ecrad_aerosol_type),INTENT(inout) :: &
      &  ecrad_aerosol                     !< ecRad aerosol information (input)
    REAL(wp), INTENT(in), OPTIONAL :: &
      &  od_lw, ssa_lw, g_lw,   & !< Optical depth, single scattering albedo, assymetry factor long wave
      &  od_sw, ssa_sw, g_sw      !< Optical depth, single scattering albedo, assymetry factor short wave
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc
! Local variables
    INTEGER                  :: &
      &  n_bands_sw,            & !< Number of ecrad shortwave bands
      &  n_bands_lw,            & !< Number of ecrad longwave bands
      &  jc, jk, jband            !< Loop indices

    CALL assert_acc_device_only("nwp_ecrad_prep_aerosol_constant", lacc)

    !$ACC DATA PRESENT(ecrad_conf, ecrad_aerosol, ssa_lw, od_lw, g_lw, ssa_sw, od_sw, g_sw)

    IF (ecrad_conf%do_lw) THEN
      n_bands_lw = get_nbands_lw_aerosol(ecrad_conf)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(3)
      DO jband = 1, n_bands_lw
        DO jk = slev, nlev
          DO jc = i_startidx, i_endidx
            ecrad_aerosol%od_lw(jband,jk,jc)  = 0._wp
            ecrad_aerosol%ssa_lw(jband,jk,jc) = 0._wp
            ecrad_aerosol%g_lw(jband,jk,jc)   = 0._wp
            IF ( PRESENT(od_lw) ) THEN
              ecrad_aerosol%od_lw(jband,jk,jc)  = od_lw
            ENDIF
            IF ( PRESENT(ssa_lw) ) THEN
              ecrad_aerosol%ssa_lw(jband,jk,jc) = ssa_lw
            ENDIF
            IF ( PRESENT(g_lw) ) THEN
              ecrad_aerosol%g_lw(jband,jk,jc)   = g_lw
            ENDIF
          ENDDO ! jc
        ENDDO ! jk
      ENDDO ! jband
      !$ACC END PARALLEL
    ENDIF

    IF (ecrad_conf%do_sw) THEN
      n_bands_sw = get_nbands_sw_aerosol(ecrad_conf)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(3)
      DO jband = 1, n_bands_sw
        DO jk = slev, nlev
          DO jc = i_startidx, i_endidx
            ecrad_aerosol%od_sw(jband,jk,jc)  = 0._wp
            ecrad_aerosol%ssa_sw(jband,jk,jc) = 0._wp
            ecrad_aerosol%g_sw(jband,jk,jc)   = 0._wp
            IF ( PRESENT(od_sw) ) THEN
              ecrad_aerosol%od_sw(jband,jk,jc)  = od_sw
            ENDIF
            IF ( PRESENT(ssa_sw) ) THEN
              ecrad_aerosol%ssa_sw(jband,jk,jc) = ssa_sw
            ENDIF
            IF ( PRESENT(g_sw) ) THEN
              ecrad_aerosol%g_sw(jband,jk,jc)   = g_sw
            ENDIF
          ENDDO ! jc
        ENDDO ! jk
      ENDDO ! jband
      !$ACC END PARALLEL
    ENDIF

    !$ACC END DATA

  END SUBROUTINE nwp_ecrad_prep_aerosol_constant
  !---------------------------------------------------------------------------------------



  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE nwp_ecrad_prep_aerosol_tegen
  !! Prepare aerosol from Tegen climatology for ecRad. All the necessary
  !! information on the vertical and spatial distribution of the aerosol is given by the
  !! preprocessed fields zaeq1 - zaeq5. Code taken and adapted from rrtm, module mo_radiation
  !!
  !---------------------------------------------------------------------------------------
  SUBROUTINE nwp_ecrad_prep_aerosol_tegen ( slev, nlev, i_startidx, i_endidx,       &
    &                                       zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,      &
    &                                       ecrad_conf, ecrad_aerosol, lacc )
    INTEGER, INTENT(in)      :: &
      &  slev, nlev,            & !< Start and end index of vertical loop
      &  i_startidx, i_endidx     !< Start and end index of horizontal loop
    REAL(wp), INTENT(in)     :: &
      &  zaeq1(:,:),            & !< Tegen optical thicknesses         1: continental
      &  zaeq2(:,:),            & !<   relative to 550 nm, including   2: maritime
      &  zaeq3(:,:),            & !<   a vertical profile              3: desert
      &  zaeq4(:,:),            & !<   for 5 different                 4: urban
      &  zaeq5(:,:)               !<   aerosol species.                5: stratospheric background
    TYPE(t_ecrad_conf),        INTENT(in)    :: &
      &  ecrad_conf               !< ecRad configuration object
    TYPE(t_ecrad_aerosol_type),INTENT(inout) :: &
      &  ecrad_aerosol            !< ecRad aerosol information (input)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc
! Local variables
    REAL(wp)                 :: &
      &  tau_abs, tau_sca         !< Absorption and scattering optical depth
    REAL(wp), POINTER        :: &
      &  scal_abs(:,:),         & !< Scaling factor absorption
      &  scal_sct(:,:),         & !< Scaling factor scattering
      &  scal_asy(:,:)            !< Scaling factor asymmetry
    CHARACTER(len=*), PARAMETER :: &
      &  routine = modname//'::nwp_ecrad_prep_aerosol_tegen' 
    INTEGER                  :: &
      &  jc, jk, jband,         & !< Loop indices
      &  n_bands_sw,            & !< Number of ecrad shortwave bands
      &  n_bands_lw,            & !< Number of ecrad longwave bands
      &  jband_shift              !< Band index in container (for shortwave: shifted by n_bands_lw)

    scal_abs => tegen_scal_factors%absorption
    scal_sct => tegen_scal_factors%scattering
    scal_asy => tegen_scal_factors%asymmetry

    n_bands_lw = get_nbands_lw_aerosol(ecrad_conf)
    n_bands_sw = get_nbands_sw_aerosol(ecrad_conf)

    CALL assert_acc_device_only("nwp_ecrad_prep_aerosol_tegen", lacc)

    !$ACC DATA PRESENT(ecrad_conf, ecrad_aerosol, zaeq1, zaeq2, zaeq3, zaeq4, zaeq5, scal_abs, scal_sct, scal_asy)

! LONGWAVE
    IF (ecrad_conf%do_lw) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(3)
      DO jk = slev, nlev
!NEC$ nointerchange
!NEC$ nounroll
        DO jband = 1, n_bands_lw
          DO jc = i_startidx, i_endidx
            ! LW optical thickness
            ecrad_aerosol%od_lw (jband,jk,jc) =  zaeq1(jc,jk) * scal_abs(jband,1) &
              &                                + zaeq2(jc,jk) * scal_abs(jband,2) &
              &                                + zaeq3(jc,jk) * scal_abs(jband,3) &
              &                                + zaeq4(jc,jk) * scal_abs(jband,4) &
              &                                + zaeq5(jc,jk) * scal_abs(jband,5)
            ! No scattering at aerosol in longwave
            ecrad_aerosol%ssa_lw(jband,jk,jc) = 0._wp
            ecrad_aerosol%g_lw  (jband,jk,jc) = 0._wp
          ENDDO ! jc
        ENDDO ! jband
      ENDDO ! jk
      !$ACC END PARALLEL
    ENDIF

! SHORTWAVE
    IF (ecrad_conf%do_sw) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(3) PRIVATE(jband_shift, tau_abs, tau_sca)
      DO jk = slev, nlev
!NEC$ nointerchange
        DO jband = 1, n_bands_sw
          DO jc = i_startidx, i_endidx

            jband_shift = n_bands_lw + jband

            ! SW absorption optical depth
            tau_abs =  zaeq1(jc,jk) * scal_abs(jband_shift,1) &
              &      + zaeq2(jc,jk) * scal_abs(jband_shift,2) &
              &      + zaeq3(jc,jk) * scal_abs(jband_shift,3) &
              &      + zaeq4(jc,jk) * scal_abs(jband_shift,4) &
              &      + zaeq5(jc,jk) * scal_abs(jband_shift,5)
            ! SW scattering optical depth
            tau_sca =  zaeq1(jc,jk) * scal_sct(jband_shift,1) &
              &      + zaeq2(jc,jk) * scal_sct(jband_shift,2) &
              &      + zaeq3(jc,jk) * scal_sct(jband_shift,3) &
              &      + zaeq4(jc,jk) * scal_sct(jband_shift,4) &
              &      + zaeq5(jc,jk) * scal_sct(jband_shift,5)
    
            ! Total optical depth for band jband
            ecrad_aerosol%od_sw (jband,jk,jc) = tau_abs + tau_sca
            
            ! Bulk SW single scattering albedo for band jband
            ecrad_aerosol%ssa_sw(jband,jk,jc) = tau_sca / ( tau_abs + tau_sca )
    
            ! Bulk SW asymmetry factor
            ecrad_aerosol%g_sw  (jband,jk,jc) =                                        &
              & (   zaeq1(jc,jk) * scal_sct(jband_shift,1) * scal_asy(jband_shift,1)   &
              &   + zaeq2(jc,jk) * scal_sct(jband_shift,2) * scal_asy(jband_shift,2)   &
              &   + zaeq3(jc,jk) * scal_sct(jband_shift,3) * scal_asy(jband_shift,3)   &
              &   + zaeq4(jc,jk) * scal_sct(jband_shift,4) * scal_asy(jband_shift,4)   &
              &   + zaeq5(jc,jk) * scal_sct(jband_shift,5) * scal_asy(jband_shift,5) ) / tau_sca
          ENDDO ! jc
        ENDDO ! jband
      ENDDO ! jk
      !$ACC END PARALLEL
    ENDIF

    !$ACC END DATA

  END SUBROUTINE nwp_ecrad_prep_aerosol_tegen

  !---------------------------------------------------------------------------------------


  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE nwp_ecrad_prep_aerosol_td
  !! Time-dependent aerosol
  !!
  !---------------------------------------------------------------------------------------
  SUBROUTINE nwp_ecrad_prep_aerosol_td (slev, nlev, i_startidx, i_endidx, &
    &                                   opt_ptrs_lw, opt_ptrs_sw,         &
    &                                   ecrad_conf, ecrad_aerosol, lacc)

    INTEGER, INTENT(in)      :: &
      &  slev, nlev,            & !< Start and end index of vertical loop
      &  i_startidx, i_endidx     !< Start and end index of horizontal loop

    TYPE(t_ecrad_conf),        INTENT(in)    :: &
      &  ecrad_conf               !< ecRad configuration object
    TYPE(t_ecrad_aerosol_type),INTENT(inout) :: &
      &  ecrad_aerosol            !< ecRad aerosol information (input)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    INTEGER                  :: &
      &  jc, jk, jband            !< Loop indices

    TYPE(t_opt_ptrs), DIMENSION(ecrad_conf%n_bands_lw), INTENT(in):: opt_ptrs_lw
    TYPE(t_opt_ptrs), DIMENSION(ecrad_conf%n_bands_sw), INTENT(in):: opt_ptrs_sw

    CALL assert_acc_host_only("nwp_ecrad_prep_aerosol_td", lacc)

   ! LONGWAVE
    IF (ecrad_conf%do_lw) THEN
      DO jband = 1, ecrad_conf%n_bands_lw
        DO jk = slev, nlev
          DO jc = i_startidx, i_endidx
            ! LW optical thickness
            ecrad_aerosol%od_lw  (jband,jk,jc) = opt_ptrs_lw(jband)%ptr_od(jc,jk)
            ! No scattering at aerosol in longwave
            ecrad_aerosol%ssa_lw (jband,jk,jc) = 0._wp
            ecrad_aerosol%g_lw   (jband,jk,jc) = 0._wp
          ENDDO ! jc
        ENDDO   ! jk
      ENDDO     ! jband
    ENDIF

   !SHORTWAVE
    IF (ecrad_conf%do_sw) THEN
      DO jband = 1, ecrad_conf%n_bands_sw
        DO jk = slev, nlev
          DO jc = i_startidx, i_endidx
            ecrad_aerosol%od_sw  (jband,jk,jc) = opt_ptrs_sw(jband)%ptr_od  (jc,jk)
            ecrad_aerosol%ssa_sw (jband,jk,jc) = opt_ptrs_sw(jband)%ptr_ssa (jc,jk)
            ecrad_aerosol%g_sw   (jband,jk,jc) = opt_ptrs_sw(jband)%ptr_g   (jc,jk)
          ENDDO ! jc
        ENDDO   ! jk
      ENDDO     ! jband
    ENDIF
  END SUBROUTINE nwp_ecrad_prep_aerosol_td
  !---------------------------------------------------------------------------------------

  SUBROUTINE nwp_ecrad_prep_aerosol_CAMS ( slev, nlev, i_startidx, i_endidx, &
    &                                      ecrad_aerosol, ptr_camsaermr)

    INTEGER, INTENT(in)                      :: &
      &  slev, nlev,            & !< Start and end index of vertical loop
      &  i_startidx, i_endidx     !< Start and end index of horizontal loop
    TYPE(t_ptr_2d), INTENT(in)               :: &
      &  ptr_camsaermr(11)
    TYPE(t_ecrad_aerosol_type),INTENT(inout) :: &
      &  ecrad_aerosol            !< ecRad aerosol information (input)
    ! Local variables
    CHARACTER(len=*), PARAMETER :: &
      &  routine = modname//'::nwp_ecrad_prep_aerosol_CAMS' 
    INTEGER                  :: &
      &  jc, jk, js              !< Loop indices

    DO js = 1, n_camsaermr
      DO jk = slev, nlev
        DO jc = i_startidx, i_endidx
          ecrad_aerosol%mixing_ratio(jc,jk,js)= ptr_camsaermr(js)%p(jc,jk)
        ENDDO
      ENDDO
    ENDDO    

  END SUBROUTINE nwp_ecrad_prep_aerosol_CAMS
#endif
END MODULE mo_nwp_ecrad_prep_aerosol
