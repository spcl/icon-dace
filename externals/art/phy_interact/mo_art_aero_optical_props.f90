!
! mo_art_aero_optical_props
! This module provides routines to calculate optical properties
! of aerosol at AERNONET and Ceilometer wavelengths
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

MODULE mo_art_aero_optical_props
! ICON
  USE mo_exception,                     ONLY: finish
  USE mo_kind,                          ONLY: wp
  USE mo_impl_constants,                ONLY: SUCCESS, MAX_CHAR_LENGTH
! ART
  USE mo_art_data,                      ONLY: t_art_data
  USE mo_art_modes_linked_list,         ONLY: t_mode
  USE mo_art_emiss_types,               ONLY: t_art_emiss2tracer
  USE mo_art_config,                    ONLY: art_config

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_calc_aod, art_calc_bsc, art_calc_aodvar

CONTAINS

SUBROUTINE art_calc_aod(rho, tracer, dz, istart, iend, nlev, jb, jg, p_art_data)

!<
! SUBROUTINE art_calc_aod
! This subroutine calculates the aerosol optical depth at 9 different
! Wavelengths used by AERONET
! Extinction coefficient values calculated by P. Gasch, 2015 (Dust)
! and C. Walter (2015) (Ash)
! Part of Module: mo_art_aero_optical_props
! Author: Daniel Rieger, KIT
! Initial Release: 2014-08-04
! Modifications:
! 2014-11-24: Daniel Rieger, KIT
! - Put block loop into interface and adapted values to 550 nm
! 2015-12-15: Philipp Gasch, Carolin Walter, KIT
! - Combined dust/seas AOD routines to a more
!   generic routine for every aerosol
!>

  REAL(wp), INTENT(in)   :: &
    &  rho(:,:),            & !< Air density
    &  tracer(:,:,:),       & !< Tracer mixing ratios [kg/kg]
    &  dz(:,:)                !< Layer height
  INTEGER, INTENT(in)    :: &
    &  istart, iend,        & !< Start and end index of nproma loop
    &  nlev, jb, jg           !< Number of verical levels, Block index
  TYPE(t_art_data),INTENT(inout) :: &
    &  p_art_data             !< Data container for ART
  ! Local variables
  REAL(wp),ALLOCATABLE   :: &
    &  ext_so4_sol(:,:,:),  & !< wavelength specific AODs
    &  ext_ash_insol(:,:,:),&
    &  ext_ash_mixed(:,:,:),&   
    &  ext_ash_giant(:,:,:)
  REAL(wp),ALLOCATABLE   :: &
    &  ext_dust(:,:),       &
    &  ext_seas(:,:),       &
    &  ext_volc(:,:),		&
	&  ext_soot(:,:),       &
    &  ext_dust_insol(:,:),       &
    &  ext_seas_sol(:,:),       &
	&  ext_soot_insol(:,:)
    
  REAL(wp)               :: &
    &  ext_so4_sol_ait(1),      & !< extinction coefficients from Mie calculations for 
                                  ! soluble Aitken mode for 550nm
    &  ext_so4_sol_acc(1),      & !< extinction coefficients from Mie calculations for
                                  ! soluble accumulation mode for 550nm
    &  ext_ash_insol_acc(1),    & !< extinction coefficients from Mie calculations for
                                  ! insoluble accumulation mode for 550nm
    &  ext_ash_insol_coa(1),    & !< extinction coefficients from Mie calculations for
                                  ! insoluble coarse mode for 550nm
    &  ext_ash_mixed_acc(1),    & !< extinction coefficients from Mie calculations for
                                  ! mixed accumulation mode for 550nm
    &  ext_ash_mixed_coa(1),    & !< extinction coefficients from Mie calculations for
                                  ! mixed coarse mode for 550nm
    &  ext_ash_giant_mode(1),   & !< extinction coefficients from Mie calculations for
                                  ! giant mode for 550nm

    &  ext_dust_m1(1:9),    & !< extinction coefficients from Mie calculations for mode 1
                              !  for each wavelength
    &  ext_dust_m2(1:9),    & !< extinction coefficients from Mie calculations for mode 2
                              !  for each wavelength
    &  ext_dust_m3(1:9),    & !< extinction coefficients from Mie calculations for mode 3
                              !  for each wavelength
    &  ext_seas_m1(1:9),    &
    &  ext_seas_m2(1:9),    &
    &  ext_seas_m3(1:9),    &
    &  ext_volc_m1(1:9),    &
    &  ext_volc_m2(1:9),    &
    &  ext_volc_m3(1:9),    &
    &  ext_soot_m1(1:9)
  INTEGER                :: &
    &  jc,jk,i_wavel,       & !< Loop indices (nproma, nlev)
    &  idusta, idustb,      & !< Tracer container indices
    &  idustc, iseasa,      & !< Tracer container indices
    &  iseasb, iseasc,      & !< Tracer container indices
    &  iasha, iashb,        & !< Tracer container indices
    &  iashc,               & !< Tracer container indices
    &  isoot,               & !< Tracer container indices
    &  ierror,              & !< error return value
    &  iso4_sol_ait,        & !< Tracer container indices
    &  iso4_sol_acc,        & !< Tracer container indices
    &  iash_insol_acc,      & !< Tracer container indices
    &  iash_insol_coa,      & !< Tracer container indices
    &  iash_mixed_acc,      & !< Tracer container indices
    &  iash_mixed_coa,      & !< Tracer container indices
    &  iash_giant,          & !< Tracer container indices
    &  isoot_insol_ait,     & !< Tracer container indices
    &  isoot_insol_acc,     & !< Tracer container indices
    &  idust_insol_acc,     & !< Tracer container indices
    &  idust_insol_coa,     & !< Tracer container indices
    &  ina_sol_acc,         & !< Tracer container indices
    &  ina_sol_coa,         & !< Tracer container indices
    &  icl_sol_acc,         & !< Tracer container indices
    &  icl_sol_coa,         & !< Tracer container indices
    &  idust_giant         !< Tracer container indices
    
  LOGICAL ::                &
      &  is_aerodyn                     !< determines if aerodyn is used
         
  CHARACTER(LEN=MAX_CHAR_LENGTH)  :: &
    &  thisroutine = "mo_art_aero_optical_props:art_calc_aod"

  TYPE(t_mode), POINTER  :: &
    &  current

  IF (art_config(jg)%iart_dust > 0) THEN
    CALL p_art_data%dict_tracer%get('dusta',idusta,ierror)
      IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "dusta" not found in dictionary.')
    CALL p_art_data%dict_tracer%get('dustb',idustb,ierror)
      IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "dustb" not found in dictionary.')
    CALL p_art_data%dict_tracer%get('dustc',idustc,ierror)
      IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "dustc" not found in dictionary.')
  ENDIF

  IF (art_config(jg)%iart_seasalt > 0) THEN
    CALL p_art_data%dict_tracer%get('seasa',iseasa,ierror)
      IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "seasa" not found in dictionary.')
    CALL p_art_data%dict_tracer%get('seasb',iseasb,ierror)
      IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "seasb" not found in dictionary.')
    CALL p_art_data%dict_tracer%get('seasc',iseasc,ierror)
      IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "seasc" not found in dictionary.')
  ENDIF

  IF (art_config(jg)%iart_volcano == 2) THEN
    CALL p_art_data%dict_tracer%get('asha',iasha,ierror)
      IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "asha" not found in dictionary.')
    CALL p_art_data%dict_tracer%get('ashb',iashb,ierror)
      IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "ashb" not found in dictionary.')
    CALL p_art_data%dict_tracer%get('ashc',iashc,ierror)
      IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "ashc" not found in dictionary.')
  ENDIF


  iso4_sol_ait = 0
  iso4_sol_acc = 0
  iash_insol_acc = 0
  iash_insol_coa = 0
  iash_mixed_acc = 0
  iash_mixed_coa = 0
  iash_giant = 0

  CALL p_art_data%dict_tracer%get('so4_sol_ait',iso4_sol_ait)
  CALL p_art_data%dict_tracer%get('so4_sol_acc',iso4_sol_acc)
  CALL p_art_data%dict_tracer%get('ash_insol_acc',iash_insol_acc)
  CALL p_art_data%dict_tracer%get('ash_insol_coa',iash_insol_coa)
  CALL p_art_data%dict_tracer%get('ash_mixed_acc',iash_mixed_acc)
  CALL p_art_data%dict_tracer%get('ash_mixed_coa',iash_mixed_coa)
  CALL p_art_data%dict_tracer%get('ash_giant',iash_giant)  
  CALL p_art_data%dict_tracer%get('soot_insol_ait',isoot_insol_ait)    
  CALL p_art_data%dict_tracer%get('soot_insol_acc',isoot_insol_acc)    
  CALL p_art_data%dict_tracer%get('dust_insol_acc',idust_insol_acc)   
  CALL p_art_data%dict_tracer%get('dust_insol_coa',idust_insol_coa)   
  CALL p_art_data%dict_tracer%get('na_sol_acc',ina_sol_acc)
  CALL p_art_data%dict_tracer%get('na_sol_coa',ina_sol_coa)
  CALL p_art_data%dict_tracer%get('cl_sol_acc',icl_sol_acc)
  CALL p_art_data%dict_tracer%get('cl_sol_coa',icl_sol_coa)
  CALL p_art_data%dict_tracer%get('dust_giant',idust_giant)

  ! Calculate AOD at specific wavelengths : SO4 - just for 550nm

  IF (iso4_sol_ait > 0 .OR. iso4_sol_acc > 0 ) THEN
    ALLOCATE(ext_so4_sol(1,iend,nlev))
    ext_so4_sol(:,:,:) = 0._wp
  ENDIF

  IF (iso4_sol_ait > 0) THEN
    ext_so4_sol_ait = (/ 0.1939_wp/)
    DO jk=1,nlev
      DO jc = istart, iend
!!        DO i_wavel=1
          i_wavel=1
          ext_so4_sol(i_wavel,jc,jk) = ext_so4_sol(i_wavel,jc,jk) &
            &     + ext_so4_sol_ait(i_wavel) * tracer(jc,jk,iso4_sol_ait) * rho(jc,jk) *1.e-6_wp
!!        ENDDO
      ENDDO
    ENDDO
  ENDIF
  IF (iso4_sol_acc > 0) THEN
    ext_so4_sol_acc = (/2.6912_wp/)
    DO jk=1,nlev
      DO jc = istart, iend
!!        DO i_wavel=1
          i_wavel=1
          ext_so4_sol(i_wavel,jc,jk) = ext_so4_sol(i_wavel,jc,jk) &
            &     + ext_so4_sol_acc(i_wavel) * tracer(jc,jk,iso4_sol_acc) * rho(jc,jk) *1.e-6_wp
!!        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF (iso4_sol_ait > 0 .OR. iso4_sol_acc > 0 ) THEN
    DO jk=1,nlev
      DO jc = istart, iend
!!        DO i_wavel=1
          i_wavel=1
          p_art_data%diag%so4_sol_aeronet(i_wavel)%tau(jc,jk,jb) = ext_so4_sol(i_wavel,jc,jk)  & 
            &                                                    * dz(jc,jk)
!!        ENDDO ! i_wavel
      ENDDO !jc
    ENDDO !jk
    DEALLOCATE(ext_so4_sol)
  ENDIF


 ! Calculate AOD at specific wavelengths : ash insol - just for 550nm
  IF (iash_insol_acc > 0 .OR. iash_insol_coa > 0 ) THEN
    ALLOCATE(ext_ash_insol(1,iend,nlev))
    ext_ash_insol(:,:,:) = 0._wp
  ENDIF

  IF (iash_insol_acc > 0) THEN
    ext_ash_insol_acc = (/ 1.60644_wp/)
    DO jk=1,nlev
      DO jc = istart, iend
!!        DO i_wavel=1
          i_wavel=1
          ext_ash_insol(i_wavel,jc,jk) = ext_ash_insol(i_wavel,jc,jk) &
            &    + ext_ash_insol_acc(i_wavel) * tracer(jc,jk,iash_insol_acc) * rho(jc,jk) *1.e-6_wp
!!        ENDDO
      ENDDO
    ENDDO
  ENDIF
  IF (iash_insol_coa > 0) THEN
    ext_ash_insol_coa = (/0.15129_wp/)
    DO jk=1,nlev
      DO jc = istart, iend
!!        DO i_wavel=1
          i_wavel=1
          ext_ash_insol(i_wavel,jc,jk) = ext_ash_insol(i_wavel,jc,jk) &
            &    + ext_ash_insol_coa(i_wavel) * tracer(jc,jk,iash_insol_coa) * rho(jc,jk) *1.e-6_wp
!!        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF (iash_insol_acc > 0 .OR. iash_insol_coa > 0 ) THEN
    DO jk=1,nlev
      DO jc = istart, iend
!!        DO i_wavel=1
          i_wavel=1
          p_art_data%diag%ash_insol_aeronet(i_wavel)%tau(jc,jk,jb) = ext_ash_insol(i_wavel,jc,jk) &
            &                                                      * dz(jc,jk)
!!        ENDDO ! i_wavel
      ENDDO !jc
    ENDDO !jk
    DEALLOCATE(ext_ash_insol)
  ENDIF

 ! Calculate AOD at specific wavelengths : ash mixed - just for 550nm
  IF (iash_mixed_acc > 0 .OR. iash_mixed_coa > 0 ) THEN
    ALLOCATE(ext_ash_mixed(1,iend,nlev))
    ext_ash_mixed(:,:,:) = 0._wp
  ENDIF

  IF (iash_mixed_acc > 0) THEN
    ext_ash_mixed_acc = (/ 4.70122_wp/)
    DO jk=1,nlev
      DO jc = istart, iend
!!        DO i_wavel=1
          i_wavel=1
          ext_ash_mixed(i_wavel,jc,jk) = ext_ash_mixed(i_wavel,jc,jk) &
            &    + ext_ash_mixed_acc(i_wavel) * tracer(jc,jk,iash_mixed_acc) * rho(jc,jk) *1.e-6_wp
!!        ENDDO
      ENDDO
    ENDDO
  ENDIF
  IF (iash_mixed_coa > 0) THEN
    ext_ash_mixed_coa = (/0.26733_wp/)
    DO jk=1,nlev
      DO jc = istart, iend
!!        DO i_wavel=1
          i_wavel=1
          ext_ash_mixed(i_wavel,jc,jk) = ext_ash_mixed(i_wavel,jc,jk) &
            &    + ext_ash_mixed_coa(i_wavel) * tracer(jc,jk,iash_mixed_coa) * rho(jc,jk) *1.e-6_wp
!!        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF (iash_mixed_acc > 0 .OR. iash_mixed_coa > 0 ) THEN
    DO jk=1,nlev
      DO jc = istart, iend
!!        DO i_wavel=1
          i_wavel=1
          p_art_data%diag%ash_mixed_aeronet(i_wavel)%tau(jc,jk,jb) = ext_ash_mixed(i_wavel,jc,jk) &
            &                                                      * dz(jc,jk)
!!        ENDDO ! i_wavel
      ENDDO !jc
    ENDDO !jk
    DEALLOCATE(ext_ash_mixed)
  ENDIF

 ! Calculate AOD at specific wavelengths : ash giant - just for 550nm
  IF (iash_giant > 0 ) THEN
    ALLOCATE(ext_ash_giant(1,iend,nlev))
    ext_ash_giant(:,:,:) = 0._wp
  ENDIF

  IF (iash_giant > 0) THEN
    ext_ash_giant_mode = (/ 0.06923_wp/)
    DO jk=1,nlev
      DO jc = istart, iend
!!        DO i_wavel=1
          i_wavel=1
          ext_ash_giant(i_wavel,jc,jk) = ext_ash_giant(i_wavel,jc,jk) &
            &     + ext_ash_giant_mode(i_wavel) * tracer(jc,jk,iash_giant) * rho(jc,jk) *1.e-6_wp
!!        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF (iash_giant > 0 ) THEN
    DO jk=1,nlev
      DO jc = istart, iend
!!        DO i_wavel=1
          i_wavel=1
          p_art_data%diag%ash_giant_aeronet(i_wavel)%tau(jc,jk,jb) = ext_ash_giant(i_wavel,jc,jk) &
            &                                                      * dz(jc,jk)
!!        ENDDO ! i_wavel
      ENDDO !jc
    ENDDO !jk
    DEALLOCATE(ext_ash_giant)
  ENDIF
 
  IF (art_config(jg)%lart_diag_out) THEN
    ! ----------------------------------
    ! Calculate AOD at specific wavelengths: dust
    ! ----------------------------------
    IF (art_config(jg)%iart_dust > 0) THEN
      SELECT CASE(art_config(jg)%iart_nonsph)
        CASE (1)! Ellipsoid Mixed 35 Shapes
          ext_dust_m1 = (/ 1.51852_wp, 1.54100_wp, 1.56019_wp, 1.55844_wp, &
            &              1.54720_wp, 1.48451_wp, 1.33041_wp, 1.18977_wp, &
            &              1.15010_wp/)
          ext_dust_m2 = (/ 0.24463_wp, 0.24786_wp, 0.25357_wp, 0.25830_wp, &
            &              0.26251_wp, 0.27245_wp, 0.28714_wp, 0.29419_wp, &
            &              0.29738_wp/)
          ext_dust_m3 = (/ 0.10144_wp, 0.10233_wp, 0.10394_wp, 0.10528_wp, &
            &              0.10651_wp, 0.10943_wp, 0.11356_wp, 0.11410_wp, &
            &              0.11486_wp/)
        CASE (0)! USE Mie values
          ext_dust_m1 = (/ 1.08098_wp, 1.10009_wp, 1.14874_wp, 1.20190_wp, &
            &              1.22450_wp, 1.27422_wp, 1.32228_wp, 1.30219_wp, &
            &              1.27715_wp/)
          ext_dust_m2 = (/ 0.20341_wp, 0.20418_wp, 0.20764_wp, 0.20869_wp, &
            &              0.21023_wp, 0.21165_wp, 0.21636_wp, 0.21884_wp, &
            &              0.21912_wp/)
          ext_dust_m3 = (/ 0.09030_wp, 0.09053_wp, 0.09128_wp, 0.09205_wp, &
            &              0.09276_wp, 0.09380_wp, 0.09282_wp, 0.09538_wp, &
            &              0.09409_wp/)
      END SELECT
      IF (idust_insol_acc > 0 .OR. idust_insol_coa > 0 ) THEN

      ALLOCATE(ext_dust_insol(iend,nlev))
      DO i_wavel = 1, 9
        DO jk = 1, nlev
          DO jc = istart, iend
            ext_dust_insol(jc,jk) = ( ext_dust_m1(i_wavel) * tracer(jc,jk,idust_insol_acc) &
                          &   + ext_dust_m2(i_wavel) * tracer(jc,jk,idust_insol_coa) &
                          &   + ext_dust_m3(i_wavel) * tracer(jc,jk,idust_giant) &
                          &   ) * rho(jc,jk) * 1.e-6_wp
            p_art_data%diag%dust_aeronet(i_wavel)%tau(jc,jk,jb) =           &
                          & ext_dust_insol(jc,jk) * dz(jc,jk)
          ENDDO !jc
        ENDDO !jk
      ENDDO ! i_wavel
      DEALLOCATE(ext_dust_insol)
      
      ELSE

      ALLOCATE(ext_dust(iend,nlev))
      DO i_wavel = 1, 9
        DO jk = 1, nlev
          DO jc = istart, iend
            ext_dust(jc,jk) = ( ext_dust_m1(i_wavel) * tracer(jc,jk,idusta) &
                          &   + ext_dust_m2(i_wavel) * tracer(jc,jk,idustb) &
                          &   + ext_dust_m3(i_wavel) * tracer(jc,jk,idustc) &
                          &   ) * rho(jc,jk) * 1.e-6_wp
            p_art_data%diag%dust_aeronet(i_wavel)%tau(jc,jk,jb) =           &
                          & ext_dust(jc,jk) * dz(jc,jk)
          ENDDO !jc
        ENDDO !jk
      ENDDO ! i_wavel
      DEALLOCATE(ext_dust)
      ENDIF
            
    ENDIF

    IF (art_config(jg)%iart_seasalt > 0) THEN
      ext_seas_m1 = (/ 5.63021_wp, 5.37511_wp, 4.97794_wp, 4.55793_wp, &
        &              4.18355_wp, 3.34615_wp, 2.36107_wp, 1.81212_wp, &
        &              1.67608_wp/) 
      ext_seas_m2 = (/ 0.39173_wp, 0.39347_wp, 0.39831_wp, 0.39967_wp, &
        &              0.40494_wp, 0.41040_wp, 0.42740_wp, 0.43641_wp, &
        &              0.43613_wp/) 
      ext_seas_m3 = (/ 0.10179_wp, 0.10175_wp, 0.10324_wp, 0.10367_wp, &
        &              0.10353_wp, 0.10317_wp, 0.10527_wp, 0.10344_wp, &
        &              0.10625_wp/) 
    IF (ina_sol_acc > 0 .OR. ina_sol_coa > 0 ) THEN

      ALLOCATE(ext_seas_sol(iend,nlev))
      DO i_wavel = 1, 9
        DO jk = 1, nlev
          DO jc = istart, iend
            ext_seas_sol(jc,jk) = ( ext_seas_m1(i_wavel) * tracer(jc,jk,ina_sol_acc) &
                          &   + ext_seas_m2(i_wavel) * tracer(jc,jk,ina_sol_coa) &
                          &   + ext_seas_m1(i_wavel) * tracer(jc,jk,icl_sol_acc) &
                          &   + ext_seas_m2(i_wavel) * tracer(jc,jk,icl_sol_coa) &
                          &   ) * rho(jc,jk) * 1.e-6_wp
            p_art_data%diag%seas_aeronet(i_wavel)%tau(jc,jk,jb) =           &
                          & ext_seas_sol(jc,jk) * dz(jc,jk)
          ENDDO !jc
        ENDDO !jk
      ENDDO ! i_wavel
      DEALLOCATE(ext_seas_sol)
      ELSE
      ALLOCATE(ext_seas(iend,nlev))
      DO i_wavel = 1, 9
        DO jk = 1, nlev
          DO jc = istart, iend
            ext_seas(jc,jk) = ( ext_seas_m1(i_wavel) * tracer(jc,jk,iseasa) &
                          &   + ext_seas_m2(i_wavel) * tracer(jc,jk,iseasb) &
                          &   + ext_seas_m3(i_wavel) * tracer(jc,jk,iseasc) &
                          &   ) * rho(jc,jk) * 1.e-6_wp
            p_art_data%diag%seas_aeronet(i_wavel)%tau(jc,jk,jb) =           &
                          & ext_seas(jc,jk) * dz(jc,jk)
          ENDDO !jc
        ENDDO !jk
      ENDDO ! i_wavel
      DEALLOCATE(ext_seas)
      ENDIF
    ENDIF

    IF (art_config(jg)%iart_volcano == 2) THEN
      !basalt
      ext_volc_m1 = (/ 0.84250_wp, 0.87836_wp, 0.85691_wp, 0.87497_wp, &
        &              0.90654_wp, 0.96593_wp, 1.12086_wp, 1.20607_wp, &
        &              1.22174_wp/)
      ext_volc_m2 = (/ 0.31586_wp, 0.31632_wp, 0.32207_wp, 0.32474_wp, &
        &              0.32193_wp, 0.33511_wp, 0.33976_wp, 0.34670_wp, &
        &              0.35050_wp/)
      ext_volc_m3 = (/ 0.14137_wp, 0.14240_wp, 0.14072_wp, 0.14193_wp, &
        &              0.14320_wp, 0.14594_wp, 0.14625_wp, 0.15059_wp, &
        &              0.15152_wp/)
      !!andesite
      !ext_volc_m1 = (/ 0.84890_wp, 0.87173_wp, 0.86897_wp, 0.89135_wp, &
      !  &              0.90668_wp, 1.00513_wp, 1.15530_wp, 1.21450_wp, &
      !  &              1.22160_wp/)
      !ext_volc_m2 = (/ 0.31302_wp, 0.31689_wp, 0.32308_wp, 0.32034_wp, &
      !  &              0.32914_wp, 0.33188_wp, 0.34809_wp, 0.35241_wp, &
      !  &              0.35591_wp/)
      !ext_volc_m3 = (/ 0.14190_wp, 0.14280_wp, 0.14226_wp, 0.14283_wp, &
      !  &              0.14450_wp, 0.14614_wp, 0.14820_wp, 0.15020_wp, &
      !  &              0.14951_wp/)
      !obsidian
      !ext_volc_m1 = (/ 0.85408_wp, 0.85202_wp, 0.88423_wp, 0.89398_wp, &
      !  &              0.90708_wp, 1.00552_wp, 1.16487_wp, 1.21921_wp, &
      !  &              1.22527_wp/)
      !ext_volc_m2 = (/ 0.31333_wp, 0.32130_wp, 0.32096_wp, 0.32424_wp, &
      !  &              0.32649_wp, 0.33153_wp, 0.34874_wp, 0.34994_wp, &
      !  &              0.35960_wp/)
      !ext_volc_m3 = (/ 0.13960_wp, 0.14066_wp, 0.14365_wp, 0.14353_wp, &
      !  &              0.14495_wp, 0.14594_wp, 0.14811_wp, 0.15086_wp, &
      !  &              0.14983_wp/)

      ALLOCATE(ext_volc(iend,nlev))
      DO i_wavel = 1, 9
        DO jk = 1, nlev
          DO jc = istart, iend
            ext_volc(jc,jk) = ( ext_volc_m1(i_wavel) * tracer(jc,jk,iasha) &
                          &   + ext_volc_m2(i_wavel) * tracer(jc,jk,iashb) &
                          &   + ext_volc_m3(i_wavel) * tracer(jc,jk,iashc) &
                          &   ) * rho(jc,jk) * 1.e-6_wp
            p_art_data%diag%volc_aeronet(i_wavel)%tau(jc,jk,jb) =          &
                          & ext_volc(jc,jk) * dz(jc,jk)
          ENDDO !jc
        ENDDO !jk
      ENDDO ! i_wavel
      DEALLOCATE(ext_volc)
    ENDIF  

    IF (art_config(jg)%iart_fire > 0) THEN
        ! coated soot wit cmd 150 nm 
      ext_soot_m1 = (/ 5.59227_wp, 5.37464_wp, 5.05216_wp, 4.68786_wp, &
        &              4.40029_wp, 3.71617_wp, 2.84327_wp, 2.31510_wp, &
        &              2.18563_wp/)
     
    IF(isoot_insol_acc > 0) THEN
      ALLOCATE(ext_soot_insol(iend,nlev))
      DO i_wavel = 1, 9
        DO jk = 1, nlev
          DO jc = istart, iend
            ext_soot_insol(jc,jk) = ( ext_soot_m1(i_wavel) * tracer(jc,jk,isoot_insol_acc) &
                          &   ) * rho(jc,jk) * 1.e-6_wp
            p_art_data%diag%soot_aeronet(i_wavel)%tau(jc,jk,jb) =           &
                          & ext_soot_insol(jc,jk) * dz(jc,jk)
          ENDDO !jc
        ENDDO !jk
      ENDDO ! i_wavel
      DEALLOCATE(ext_soot_insol)
      ELSE
      ALLOCATE(ext_soot(iend,nlev))
      DO i_wavel = 1, 9
        DO jk = 1, nlev
          DO jc = istart, iend
            ext_soot(jc,jk) = ext_soot_m1(i_wavel) * tracer(jc,jk,isoot) &
                          &   * rho(jc,jk) * 1.e-6_wp
            p_art_data%diag%soot_aeronet(i_wavel)%tau(jc,jk,jb) =           &
                          & ext_soot(jc,jk) * dz(jc,jk)
          ENDDO !jc
        ENDDO !jk
      ENDDO ! i_wavel
      DEALLOCATE(ext_soot)
      ENDIF
  
    ENDIF

  ENDIF !lart_diag_out

END SUBROUTINE art_calc_aod

SUBROUTINE art_calc_aodvar(rho, tracer, dz, istart, iend, nlev, jb, jg, p_art_data, &
  &                        cmd, ini_cmd, mode_name, l_init_aod)

!<
! SUBROUTINE art_calc_aodvar
! This subroutine calculates the aerosol optical depth at 9 different
! Wavelengths considering the variable cmd at 550 nm
! Part of Module: mo_art_aero_optical_props
! Author: Ali Hoshyaripour, KIT
! Initial Release: 2018-08-14
!>

  REAL(wp), INTENT(in)   :: &
    &  rho(:,:),            & !< Air density
    &  tracer(:,:,:),       & !< Tracer mixing ratios [kg/kg]
    &  dz(:,:),             & !< Layer height
    &  cmd(:,:),            &
    &  ini_cmd
  INTEGER, INTENT(in)    :: &
    &  istart, iend,        & !< Start and end index of nproma loop
    &  nlev, jb,            & !< Number of verical levels, Block index
    &  jg                     !< Patch id
  CHARACTER(*), INTENT(in) :: &
    &  mode_name
  LOGICAL, INTENT(in) :: &
    &  l_init_aod
  TYPE(t_art_data),INTENT(inout) :: &
    &  p_art_data             !< Data container for ART
  ! Local variables
  REAL(wp),ALLOCATABLE   :: &
    &  ext_seas(:,:),       &
    &  ext_volc(:,:),		&
	&  ext_soot(:,:),       &
    &  ext_dust_insol(:,:), &
    &  ext_seas_sol(:,:),   &
	&  ext_soot_insol(:,:)
  REAL(wp)               :: &
    &  ext_dust_m1(1:9),    & !< extinction coefficients from Mie calculations for mode 1
                              !  for each wavelength
    &  ext_dust_m2(1:9),    & !< extinction coefficients from Mie calculations for mode 2
                              !  for each wavelength
    &  ext_dust_m3(1:9),    & !< extinction coefficients from Mie calculations for mode 3
                              !  for each wavelength
    &  ext_seas_m1(1:9),    &
    &  ext_seas_m2(1:9),    &
    &  ext_seas_m3(1:9),    &
    &  ext_volc_m1(1:9),    &
    &  ext_volc_m2(1:9),    &
    &  ext_volc_m3(1:9),    &
	&  ext_soot_m1(1:9),    &
    &  fac_cmd
  INTEGER                :: &
    &  jc,jk,i_wavel,       & !< Loop indices (nproma, nlev)
    &  idusta, idustb,      & !< Tracer container indices
    &  idustc, iseasa,      & !< Tracer container indices
    &  iseasb, iseasc,      & !< Tracer container indices
    &  iasha, iashb,        & !< Tracer container indices
    &  iashc, isoot,        & !< Tracer container indices
    &  isoot_insol_ait,     & !< Tracer container indices
    &  isoot_insol_acc,     & !< Tracer container indices
    &  idust_insol_acc,     & !< Tracer container indices
    &  idust_insol_coa,     & !< Tracer container indices
    &  ina_sol_acc,         & !< Tracer container indices
    &  ina_sol_coa,         & !< Tracer container indices
    &  icl_sol_acc,         & !< Tracer container indices
    &  icl_sol_coa,         & !< Tracer container indices
    &  idust_giant,         &!< Tracer container indices
	&  ierror          		  !< Tracer container indices/ error return value
    
  LOGICAL ::                &
    &  is_aerodyn                     !< determines if aerodyn is used
  CHARACTER(LEN=MAX_CHAR_LENGTH)  :: &
    &  thisroutine = "mo_art_aero_optical_props:art_calc_aodvar"

    idusta = 0
    idustb = 0
    idustc = 0
    iseasa = 0
    iseasb = 0
    iseasc = 0
    iasha = 0
    iashb = 0
    iashc = 0
    isoot = 0 

     is_aerodyn = .FALSE.
     IF(.NOT.is_aerodyn) THEN
          IF (art_config(jg)%iart_dust > 0) THEN
            CALL p_art_data%dict_tracer%get('dusta',idusta,ierror)
              IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "dusta" not found in dictionary.')
            CALL p_art_data%dict_tracer%get('dustb',idustb,ierror)
              IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "dustb" not found in dictionary.')
            CALL p_art_data%dict_tracer%get('dustc',idustc,ierror)
              IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "dustc" not found in dictionary.')
          ENDIF

          IF (art_config(jg)%iart_seasalt > 0) THEN
            CALL p_art_data%dict_tracer%get('seasa',iseasa,ierror)
              IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "seasa" not found in dictionary.')
            CALL p_art_data%dict_tracer%get('seasb',iseasb,ierror)
              IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "seasb" not found in dictionary.')
            CALL p_art_data%dict_tracer%get('seasc',iseasc,ierror)
              IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "seasc" not found in dictionary.')
          ENDIF

          IF (art_config(jg)%iart_volcano == 2) THEN
            CALL p_art_data%dict_tracer%get('asha',iasha,ierror)
              IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "asha" not found in dictionary.')
            CALL p_art_data%dict_tracer%get('ashb',iashb,ierror)
              IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "ashb" not found in dictionary.')
            CALL p_art_data%dict_tracer%get('ashc',iashc,ierror)
              IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "ashc" not found in dictionary.')
          ENDIF
          IF (art_config(jg)%iart_fire > 0) THEN
            CALL p_art_data%dict_tracer%get('soot',isoot,ierror)
              IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "soot" not found in dictionary.')
          ENDIF
     ENDIF
     
      isoot_insol_acc = 0
      idust_insol_acc = 0
      idust_insol_coa = 0
      ina_sol_acc = 0
      ina_sol_coa = 0
      icl_sol_acc = 0
      icl_sol_coa = 0

        IF (art_config(jg)%iart_fire > 0 .AND. isoot == 0) THEN
          CALL p_art_data%dict_tracer%get('soot_insol_acc',isoot_insol_acc)
        ENDIF
        IF (art_config(jg)%iart_dust > 0 .AND. idusta == 0) THEN
          CALL p_art_data%dict_tracer%get('dust_insol_acc',idust_insol_acc)   
          CALL p_art_data%dict_tracer%get('dust_insol_coa',idust_insol_coa)
        ENDIF
        IF (art_config(jg)%iart_seasalt > 0 .AND. iseasa == 0) THEN
          CALL p_art_data%dict_tracer%get('na_sol_acc',ina_sol_acc)
          CALL p_art_data%dict_tracer%get('na_sol_coa',ina_sol_coa)
          CALL p_art_data%dict_tracer%get('cl_sol_acc',icl_sol_acc)
          CALL p_art_data%dict_tracer%get('cl_sol_coa',icl_sol_coa)
        ENDIF
  
  IF (art_config(jg)%lart_diag_out) THEN

    ! ----------------------------------
    ! --- Calculate AOD at specific wavelengths
    ! ----------------------------------
    IF (art_config(jg)%iart_dust > 0) THEN
      SELECT CASE(art_config(jg)%iart_nonsph)
      CASE (1)! Ellipsoid Mixed 35 Shapes
        ext_dust_m1 = (/ 1.51852_wp, 1.54100_wp, 1.56019_wp, 1.55844_wp, &
          &              1.54720_wp, 1.48451_wp, 1.33041_wp, 1.18977_wp, &
          &              1.15010_wp/)
        ext_dust_m2 = (/ 0.24463_wp, 0.24786_wp, 0.25357_wp, 0.25830_wp, &
          &              0.26251_wp, 0.27245_wp, 0.28714_wp, 0.29419_wp, &
          &              0.29738_wp/)
        ext_dust_m3 = (/ 0.10144_wp, 0.10233_wp, 0.10394_wp, 0.10528_wp, &
          &              0.10651_wp, 0.10943_wp, 0.11356_wp, 0.11410_wp, &
          &              0.11486_wp/)
      CASE (0)! USE Mie values
        ext_dust_m1 = (/ 1.08098_wp, 1.10009_wp, 1.14874_wp, 1.20190_wp, &
          &              1.22450_wp, 1.27422_wp, 1.32228_wp, 1.30219_wp, &
          &              1.27715_wp/)
        ext_dust_m2 = (/ 0.20341_wp, 0.20418_wp, 0.20764_wp, 0.20869_wp, &
          &              0.21023_wp, 0.21165_wp, 0.21636_wp, 0.21884_wp, &
          &              0.21912_wp/)
        ext_dust_m3 = (/ 0.09030_wp, 0.09053_wp, 0.09128_wp, 0.09205_wp, &
          &              0.09276_wp, 0.09380_wp, 0.09282_wp, 0.09538_wp, &
          &              0.09409_wp/)
      END SELECT

      IF (mode_name(:4) == 'dust' .AND. l_init_aod) THEN
        DO i_wavel = 1, 9
          DO jk = 1, nlev
            DO jc = istart, iend
              p_art_data%diag%dust_aeronet(i_wavel)%tau(jc,jk,jb) = 0.0_wp
            ENDDO !jc
          ENDDO !jk
        ENDDO ! i_wavel
      END IF

      SELECT CASE (mode_name)
      CASE ('dusta')
        DO i_wavel = 1, 9
          DO jk = 1, nlev
            DO jc = istart, iend
              ! for 550 nm only
              IF (i_wavel==5 .AND. cmd(jc,jk) > 0.25_wp * ini_cmd) THEN
                fac_cmd = ini_cmd / cmd(jc,jk)
              ELSE
                fac_cmd = 1.0_wp
              ENDIF
              p_art_data%diag%dust_aeronet(i_wavel)%tau(jc,jk,jb) =      &
                &   p_art_data%diag%dust_aeronet(i_wavel)%tau(jc,jk,jb)  &
                & + ext_dust_m1(i_wavel) * tracer(jc,jk,idusta)          &
                &   * fac_cmd * rho(jc,jk) * 1.e-6_wp * dz(jc,jk)
            ENDDO !jc
          ENDDO !jk
        ENDDO ! i_wavel
      CASE ('dustb')
        DO i_wavel = 1, 9
          DO jk = 1, nlev
            DO jc = istart, iend
              ! for 550 nm only
              IF (i_wavel==5 .AND. cmd(jc,jk) > 0.25_wp * ini_cmd) THEN
                 fac_cmd = ini_cmd / cmd(jc,jk)
              ELSE
                fac_cmd = 1.0_wp
              ENDIF
              p_art_data%diag%dust_aeronet(i_wavel)%tau(jc,jk,jb) =      &
                &   p_art_data%diag%dust_aeronet(i_wavel)%tau(jc,jk,jb)  &
                & + ext_dust_m2(i_wavel) * tracer(jc,jk,idustb)          &
                &   * fac_cmd * rho(jc,jk) * 1.e-6_wp * dz(jc,jk)
            ENDDO !jc
          ENDDO !jk
        ENDDO ! i_wavel
      CASE ('dustc')
        DO i_wavel = 1, 9
          DO jk = 1, nlev
            DO jc = istart, iend
              ! for 550 nm only
              IF (i_wavel==5 .AND. cmd(jc,jk) > 0.25_wp * ini_cmd) THEN
                 fac_cmd = ini_cmd / cmd(jc,jk)
              ELSE
                fac_cmd = 1.0_wp
              ENDIF
              p_art_data%diag%dust_aeronet(i_wavel)%tau(jc,jk,jb) =      &
                &   p_art_data%diag%dust_aeronet(i_wavel)%tau(jc,jk,jb)  &
                & + ext_dust_m3(i_wavel) * tracer(jc,jk,idustc)          &
                &   * fac_cmd * rho(jc,jk) * 1.e-6_wp * dz(jc,jk)
            ENDDO !jc
          ENDDO !jk
        ENDDO ! i_wavel
      CASE DEFAULT
        ! nothing!
      END SELECT
      !WRITE(*,*), p_art_data%diag%dusta_diam(jc,jk)
    ENDIF

    IF (art_config(jg)%iart_seasalt > 0) THEN
      ext_seas_m1 = (/ 5.63021_wp, 5.37511_wp, 4.97794_wp, 4.55793_wp, &
        &              4.18355_wp, 3.34615_wp, 2.36107_wp, 1.81212_wp, &
        &              1.67608_wp/) 
      ext_seas_m2 = (/ 0.39173_wp, 0.39347_wp, 0.39831_wp, 0.39967_wp, &
        &              0.40494_wp, 0.41040_wp, 0.42740_wp, 0.43641_wp, &
        &              0.43613_wp/) 
      ext_seas_m3 = (/ 0.10179_wp, 0.10175_wp, 0.10324_wp, 0.10367_wp, &
        &              0.10353_wp, 0.10317_wp, 0.10527_wp, 0.10344_wp, &
        &              0.10625_wp/)

      ALLOCATE(ext_seas(iend,nlev))
      DO i_wavel = 1, 9
        DO jk = 1, nlev
          DO jc = istart, iend
            ext_seas(jc,jk) = ( ext_seas_m1(i_wavel) * tracer(jc,jk,iseasa) &
                          &   + ext_seas_m2(i_wavel) * tracer(jc,jk,iseasb) &
                          &   + ext_seas_m3(i_wavel) * tracer(jc,jk,iseasc) &
                          &   ) * rho(jc,jk) * 1.e-6_wp
            p_art_data%diag%seas_aeronet(i_wavel)%tau(jc,jk,jb) =           &
                          & ext_seas(jc,jk) * dz(jc,jk)
          ENDDO !jc
        ENDDO !jk
      ENDDO ! i_wavel
      DEALLOCATE(ext_seas)
    ENDIF

    IF (art_config(jg)%iart_volcano == 2) THEN
      !basalt
      ext_volc_m1 = (/ 0.84250_wp, 0.87836_wp, 0.85691_wp, 0.87497_wp, &
        &              0.90654_wp, 0.96593_wp, 1.12086_wp, 1.20607_wp, &
        &              1.22174_wp/)
      ext_volc_m2 = (/ 0.31586_wp, 0.31632_wp, 0.32207_wp, 0.32474_wp, &
        &              0.32193_wp, 0.33511_wp, 0.33976_wp, 0.34670_wp, &
        &              0.35050_wp/)
      ext_volc_m3 = (/ 0.14137_wp, 0.14240_wp, 0.14072_wp, 0.14193_wp, &
        &              0.14320_wp, 0.14594_wp, 0.14625_wp, 0.15059_wp, &
        &              0.15152_wp/)
      !!andesite
      !ext_volc_m1 = (/ 0.84890_wp, 0.87173_wp, 0.86897_wp, 0.89135_wp, &
      !  &              0.90668_wp, 1.00513_wp, 1.15530_wp, 1.21450_wp, &
      !  &              1.22160_wp/)
      !ext_volc_m2 = (/ 0.31302_wp, 0.31689_wp, 0.32308_wp, 0.32034_wp, &
      !  &              0.32914_wp, 0.33188_wp, 0.34809_wp, 0.35241_wp, &
      !  &              0.35591_wp/)
      !ext_volc_m3 = (/ 0.14190_wp, 0.14280_wp, 0.14226_wp, 0.14283_wp, &
      !  &              0.14450_wp, 0.14614_wp, 0.14820_wp, 0.15020_wp, &
      !  &              0.14951_wp/)
      !obsidian
      !ext_volc_m1 = (/ 0.85408_wp, 0.85202_wp, 0.88423_wp, 0.89398_wp, &
      !  &              0.90708_wp, 1.00552_wp, 1.16487_wp, 1.21921_wp, &
      !  &              1.22527_wp/)
      !ext_volc_m2 = (/ 0.31333_wp, 0.32130_wp, 0.32096_wp, 0.32424_wp, &
      !  &              0.32649_wp, 0.33153_wp, 0.34874_wp, 0.34994_wp, &
      !  &              0.35960_wp/)
      !ext_volc_m3 = (/ 0.13960_wp, 0.14066_wp, 0.14365_wp, 0.14353_wp, &
      !  &              0.14495_wp, 0.14594_wp, 0.14811_wp, 0.15086_wp, &
      !  &              0.14983_wp/)

      ALLOCATE(ext_volc(iend,nlev))
      DO i_wavel = 1, 9
        DO jk = 1, nlev
          DO jc = istart, iend
            ext_volc(jc,jk) = ( ext_volc_m1(i_wavel) * tracer(jc,jk,iasha) &
                          &   + ext_volc_m2(i_wavel) * tracer(jc,jk,iashb) &
                          &   + ext_volc_m3(i_wavel) * tracer(jc,jk,iashc) &
                          &   ) * rho(jc,jk) * 1.e-6_wp
            p_art_data%diag%volc_aeronet(i_wavel)%tau(jc,jk,jb) =          &
                          & ext_volc(jc,jk) * dz(jc,jk)
          ENDDO !jc
        ENDDO !jk
      ENDDO ! i_wavel
      DEALLOCATE(ext_volc)
        ENDIF
        IF (art_config(jg)%iart_fire > 0) THEN
        ! coated soot wit cmd 150 nm 
      ext_soot_m1 = (/ 5.59227_wp, 5.37464_wp, 5.05216_wp, 4.68786_wp, &
        &              4.40029_wp, 3.71617_wp, 2.84327_wp, 2.31510_wp, &
        &              2.18563_wp/) 
            IF(isoot_insol_acc > 0) THEN
              ALLOCATE(ext_soot_insol(iend,nlev))
              DO i_wavel = 1, 9
                DO jk = 1, nlev
                  DO jc = istart, iend
                    ext_soot_insol(jc,jk) = ( ext_soot_m1(i_wavel) * tracer(jc,jk,isoot_insol_acc) &
                                  &   ) * rho(jc,jk) * 1.e-6_wp
                    p_art_data%diag%soot_aeronet(i_wavel)%tau(jc,jk,jb) =           &
                                  & ext_soot_insol(jc,jk) * dz(jc,jk)
                  ENDDO !jc
                ENDDO !jk
              ENDDO ! i_wavel
              DEALLOCATE(ext_soot_insol)
              ELSE
              ALLOCATE(ext_soot(iend,nlev))
              DO i_wavel = 1, 9
                DO jk = 1, nlev
                  DO jc = istart, iend
                    ext_soot(jc,jk) = ext_soot_m1(i_wavel) * tracer(jc,jk,isoot) &
                                  &   * rho(jc,jk) * 1.e-6_wp
                    p_art_data%diag%soot_aeronet(i_wavel)%tau(jc,jk,jb) =           &
                                  & ext_soot(jc,jk) * dz(jc,jk)
                  ENDDO !jc
                ENDDO !jk
              ENDDO ! i_wavel
              DEALLOCATE(ext_soot)
              ENDIF
        ENDIF
  ENDIF !lart_diag_out
END SUBROUTINE art_calc_aodvar

SUBROUTINE art_calc_bsc(rho,tracer, dz, istart, iend, nlev, jb, jg, p_art_data)
!<
! SUBROUTINE art_calc_bsc
! This subroutine calculates the aerosol backscatter and attenuated
! backscattter at 3 different Wavelengths used by lidars
! from ground and from satellite
! Backscatter coefficient values calculated by C. Walter (2015) (Ash)
! Part of Module: mo_art_aero_optical_props
! Author: Carolin Walter, KIT
! Initial Release: 2016-03-14
!>
  REAL(wp), INTENT(in)   :: &
    &  rho(:,:),            & !< Air density
    &  tracer(:,:,:),       & !< Tracer mixing ratios [kg/kg]
    &  dz(:,:)                !< Layer height
  INTEGER, INTENT(in)    :: &
    &  istart, iend,        & !< Start and end index of nproma loop
    &  nlev, jb, jg           !< Number of verical levels, Block index
  TYPE(t_art_data),INTENT(inout) :: &
    &  p_art_data             !< Data container for ART
  ! Local variables
  REAL(wp),ALLOCATABLE   :: &
    &  bsc_dust(:,:),       & !< wavelength specific BSCs
    &  bsc_volc(:,:),       &
    &  bsc_soot(:,:),       &
    &  ext_dust(:,:),       & !< wavelength specific EXT
    &  ext_volc(:,:),       &
    &  ext_soot(:,:),       &
    &  att_dust(:,:),       & !< wavelength specific ATT BSCs
    &  att_volc(:,:),       &
    &  att_sum_dust(:,:),   &
    &  att_sum_volc(:,:),   &
    &  att_seas(:,:),       &       
    &  att_soot(:,:)
  REAL(wp)               :: &
    &  bsc_dust_m1(1:3),    & !< backscatter coefficients from Mie calculations for mode 1
                              !  for each wavelength
    &  bsc_dust_m2(1:3),    & !< backscatter coefficients from Mie calculations for mode 2
                              !  for each wavelength
    &  bsc_dust_m3(1:3),    & !< backscatter coefficients from Mie calculations for mode 3
                              !  for each wavelength
    &  bsc_seas_m1(1:3),    &
    &  bsc_seas_m2(1:3),    &
    &  bsc_seas_m3(1:3),    &
    &  bsc_volc_m1(1:3),    &
    &  bsc_volc_m2(1:3),    &
    &  bsc_volc_m3(1:3),    &
    &  bsc_soot_m1(1:3),    &
    &  ext_dust_m1(1:3),    &
    &  ext_dust_m2(1:3),    &
    &  ext_dust_m3(1:3),    &
    &  ext_seas_m1(1:3),    &
    &  ext_seas_m2(1:3),    &
    &  ext_seas_m3(1:3),    &
    &  ext_volc_m1(1:3),    &
    &  ext_volc_m2(1:3),    &
    &  ext_volc_m3(1:3),    &
    &  ext_soot_m1(1:3)

  INTEGER                :: &
    &  jc, jk, i_wavel,     & !< Loop indices (nproma, nlev)
    &  jkm1, jkp1,          & !< layer above, layer below
    &  idusta, idustb,      & !< Tracer container indices
    &  idustc,iseasa,       & !< Tracer container indices
    &  iseasb, iseasc,      & !< Tracer container indices              & !< Tracer container indices
    &  iasha, iashb,        & !< Tracer container indices
    &  iashc, ierror,       & !< Tracer container indices/ error return value
    &  isoot,               &
    &  isoot_insol_ait,     & !< Tracer container indices
    &  isoot_insol_acc,     & !< Tracer container indices
    &  idust_insol_acc,     & !< Tracer container indices
    &  idust_insol_coa,     & !< Tracer container indices
    &  ina_sol_acc,         & !< Tracer container indices
    &  ina_sol_coa,         & !< Tracer container indices
    &  icl_sol_acc,         & !< Tracer container indices
    &  icl_sol_coa,         & !< Tracer container indices
    &  idust_giant            !< Tracer container indices
  LOGICAL ::                &
    &  is_aerodyn                     !< determines if aerodyn is used 
  CHARACTER(LEN=MAX_CHAR_LENGTH)  :: &
    &  thisroutine = "mo_art_aero_optical_props:art_calc_bsc"
  TYPE(t_mode), POINTER  :: &
    &  current

  idusta = 0
  idustb = 0
  idustc = 0
  iasha  = 0
  iashb  = 0
  iashc  = 0
  iseasa = 0
  iseasb = 0
  iseasc = 0
  isoot = 0
 
          IF (art_config(jg)%iart_dust > 0) THEN
            CALL p_art_data%dict_tracer%get('dusta',idusta,ierror)
              IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "dusta" not found in dictionary.')
            CALL p_art_data%dict_tracer%get('dustb',idustb,ierror)
              IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "dustb" not found in dictionary.')
            CALL p_art_data%dict_tracer%get('dustc',idustc,ierror)
              IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "dustc" not found in dictionary.')
          ENDIF

          IF (art_config(jg)%iart_seasalt > 0) THEN
            CALL p_art_data%dict_tracer%get('seasa',iseasa,ierror)
              IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "seasa" not found in dictionary.')
            CALL p_art_data%dict_tracer%get('seasb',iseasb,ierror)
              IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "seasb" not found in dictionary.')
            CALL p_art_data%dict_tracer%get('seasc',iseasc,ierror)
              IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "seasc" not found in dictionary.')
          ENDIF

          IF (art_config(jg)%iart_volcano == 2) THEN
            CALL p_art_data%dict_tracer%get('asha',iasha,ierror)
              IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "asha" not found in dictionary.')
            CALL p_art_data%dict_tracer%get('ashb',iashb,ierror)
              IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "ashb" not found in dictionary.')
            CALL p_art_data%dict_tracer%get('ashc',iashc,ierror)
              IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "ashc" not found in dictionary.')
          ENDIF
          IF (art_config(jg)%iart_fire > 0) THEN
            CALL p_art_data%dict_tracer%get('soot',isoot,ierror)
              IF(ierror /= SUCCESS) CALL finish (thisroutine, 'Tracer labeled "soot" not found in dictionary.')
          ENDIF
     
      isoot_insol_acc = 0
      idust_insol_acc = 0
      idust_insol_coa = 0
      ina_sol_acc = 0
      ina_sol_coa = 0
      icl_sol_acc = 0
      icl_sol_coa = 0
        IF (art_config(jg)%iart_fire > 0 .AND. isoot == 0) THEN
          CALL p_art_data%dict_tracer%get('soot_insol_acc',isoot_insol_acc)
        ENDIF
        IF (art_config(jg)%iart_dust > 0 .AND. idusta == 0) THEN
          CALL p_art_data%dict_tracer%get('dust_insol_acc',idust_insol_acc)   
          CALL p_art_data%dict_tracer%get('dust_insol_coa',idust_insol_coa)
        ENDIF
        IF (art_config(jg)%iart_seasalt > 0 .AND. iseasa == 0) THEN
          CALL p_art_data%dict_tracer%get('na_sol_acc',ina_sol_acc)
          CALL p_art_data%dict_tracer%get('na_sol_coa',ina_sol_coa)
          CALL p_art_data%dict_tracer%get('cl_sol_acc',icl_sol_acc)
          CALL p_art_data%dict_tracer%get('cl_sol_coa',icl_sol_coa)
        ENDIF
     

  IF (art_config(jg)%lart_diag_out) THEN

    ! ----------------------------------
    ! --- Calculate backscatter at specific wavelengths
    ! ----------------------------------
    IF (art_config(jg)%iart_dust > 0) THEN
      SELECT CASE(art_config(jg)%iart_nonsph)
        CASE (1) ! Ellipsoid Mixed 35 shapes
          bsc_dust_m1 = (/ 0.03650_wp,  0.02320_wp,  0.01355_wp /)
          bsc_dust_m2 = (/ 0.00114_wp,  0.00310_wp,  0.00991_wp /)
          bsc_dust_m3 = (/ 0.00023_wp,  0.00066_wp,  0.00162_wp /)
          ext_dust_m1 = (/ 1.52664_wp,  1.55218_wp,  1.15010_wp /)
          ext_dust_m2 = (/ 0.24533_wp,  0.26100_wp,  0.29738_wp /)
          ext_dust_m3 = (/ 0.10163_wp,  0.10607_wp,  0.11486_wp /)
        CASE (0)  ! Mie values (bsc crrected Meng 2010)
          bsc_dust_m1 = (/ 0.11602_wp,  0.02977_wp,  0.01305_wp /)
          bsc_dust_m2 = (/ 0.00326_wp,  0.01643_wp,  0.04098_wp /)
          bsc_dust_m3 = (/ 0.00054_wp,  0.00685_wp,  0.00740_wp /)
          ext_dust_m1 = (/ 1.10031_wp,  1.19550_wp,  1.27715_wp /)
          ext_dust_m2 = (/ 0.20446_wp,  0.20769_wp,  0.21912_wp /)
          ext_dust_m3 = (/ 0.09036_wp,  0.09139_wp,  0.09409_wp /)
      END SELECT

      ALLOCATE(bsc_dust(iend,nlev))
      DO i_wavel = 1, 3
        DO jk=1,nlev
          DO jc = istart, iend
            bsc_dust(jc,jk) = ( bsc_dust_m1(i_wavel) * tracer(jc,jk,idusta) &
                          &   + bsc_dust_m2(i_wavel) * tracer(jc,jk,idustb) &
                          &   + bsc_dust_m3(i_wavel) * tracer(jc,jk,idustc) &
                          &   ) * rho(jc,jk) * 1.e-6_wp
            p_art_data%diag%dust_ceilo(i_wavel)%bsc(jc,jk,jb) = bsc_dust(jc,jk)

          ENDDO !jc
        ENDDO !jk
      ENDDO ! i_wavel
      DEALLOCATE(bsc_dust)

      !calculate attenuated backscatter
      ALLOCATE(att_dust(iend,nlev))
      ALLOCATE(ext_dust(iend,nlev))
      
      ! ... for ceilometer
      DO i_wavel = 1, 3
        DO jk = nlev, 1, -1
          jkp1 = MIN(nlev, jk+1)
!NEC$ ivdep
          DO jc = istart, iend
            IF (jk == nlev) att_dust(jc,jk) = 0.0_wp
            ext_dust(jc,jk) = ( ext_dust_m1(i_wavel) * tracer(jc,jk,idusta) &
                          &   + ext_dust_m2(i_wavel) * tracer(jc,jk,idustb) &
                          &   + ext_dust_m3(i_wavel) * tracer(jc,jk,idustc) &
                          &   ) * rho(jc,jk) * 1.e-6_wp
            att_dust(jc,jk) = att_dust(jc,jkp1) + ext_dust(jc,jk) * dz(jc,jk)
            IF (p_art_data%diag%dust_ceilo(i_wavel)%bsc(jc,jk,jb) > 0.0_wp) THEN
              p_art_data%diag%dust_att(i_wavel)%ceil_bsc(jc,jk,jb) = &
                &  p_art_data%diag%dust_ceilo(i_wavel)%bsc(jc,jk,jb) * EXP(-2 * att_dust(jc,jk))
            ELSE
              p_art_data%diag%dust_att(i_wavel)%ceil_bsc(jc,jk,jb) = 0.0_wp
            END IF
          ENDDO !jc
        ENDDO ! jk
      ENDDO ! i_wavel
      
      ! ... for satellite
      DO i_wavel = 1, 3
        DO jk = 1, nlev
          jkm1 = MAX(1, jk-1)
!NEC$ ivdep
          DO jc = istart, iend
            IF (jk == 1) att_dust(jc,jk) = 0.0_wp
              ext_dust(jc,jk) = ( ext_dust_m1(i_wavel) * tracer(jc,jk,idusta) &
                          &   + ext_dust_m2(i_wavel) * tracer(jc,jk,idustb) &
                          &   + ext_dust_m3(i_wavel) * tracer(jc,jk,idustc) &
                          &   ) * rho(jc,jk) * 1.e-6_wp

            att_dust(jc,jk) = att_dust(jc,jkm1) + ext_dust(jc,jk) * dz(jc,jk)
            IF (p_art_data%diag%dust_ceilo(i_wavel)%bsc(jc,jk,jb) == 0.0_wp) THEN
              p_art_data%diag%dust_sat(i_wavel)%sat_bsc(jc,jk,jb) = 0.0_wp
            ELSE
              p_art_data%diag%dust_sat(i_wavel)%sat_bsc(jc,jk,jb) = &
                &  p_art_data%diag%dust_ceilo(i_wavel)%bsc(jc,jk,jb) * EXP(-2 * att_dust(jc,jk))
            END IF
          ENDDO !jc
        ENDDO !jk
      ENDDO ! i_wavel
      
      DEALLOCATE(ext_dust)
      DEALLOCATE(att_dust)
    ENDIF !iart_dust>0

!! MARKER: seasalt needs to be adjusted for aerodyn-treatment !!
!! ============================================================
!    IF (art_config(jg)%iart_seasalt > 0) THEN
!      bsc_seas_m1 = (/ 0.37411_wp,  0.27673_wp,  0.13381_wp /) 
!      bsc_seas_m2 = (/ 0.02118_wp,  0.02521_wp,  0.03019_wp /) 
!      bsc_seas_m3 = (/ 0.00423_wp,  0.00453_wp,  0.00469_wp /) 
!      ext_seas_m1 = (/ 5.60142_wp,  4.30244_wp,  1.67568_wp /) 
!      ext_seas_m2 = (/ 0.39031_wp,  0.40685_wp,  0.43607_wp /) 
!      ext_seas_m3 = (/ 0.10244_wp,  0.10235_wp,  0.10624_wp /)
!
!      ALLOCATE(bsc_seas(iend,nlev))
!      DO i_wavel = 1, 3
!        DO jk = 1, nlev
!          DO jc = istart, iend
!            bsc_seas(jc,jk) = ( bsc_seas_m1(i_wavel) * p_trac(jc,jk,iseasa) &
!                          &   + bsc_seas_m2(i_wavel) * p_trac(jc,jk,iseasb) &
!                          &   + bsc_seas_m3(i_wavel) * p_trac(jc,jk,iseasc) &
!                          &   ) * rho(jc,jk) * 1.e-6_wp
!            p_art_data%diag%seas_ceilo(i_wavel)%bsc(jc,jk,jb) = bsc_seas(jc,jk)
!          ENDDO !jc
!        ENDDO !jk
!      ENDDO ! i_wavel
!      DEALLOCATE(bsc_seas)
!
!      !calculate attenuated backscatter
!      ALLOCATE(att_seas(iend,nlev))
!      ALLOCATE(ext_seas(iend,nlev))
!
!      ! ... for ceilometer
!      DO i_wavel = 1, 3
!        DO jk = nlev, 1, -1
!          jkp1 = MIN(nlev, jk+1)
!!NEC$ ivdep
!          DO jc = istart, iend
!            IF (jk == nlev) att_seas(jc,jk) = 0.0_wp
!            ext_seas(jc,jk) = ( ext_seas_m1(i_wavel) * p_trac(jc,jk,iseasa) &
!                          &   + ext_seas_m2(i_wavel) * p_trac(jc,jk,iseasb) &
!                          &   + ext_seas_m3(i_wavel) * p_trac(jc,jk,iseasc) &
!                          &   ) * rho(jc,jk) * 1.e-6_wp
!            att_seas(jc,jk) = att_seas(jc,jkp1) + ext_seas(jc,jk) * dz(jc,jk)
!            IF (p_art_data%diag%seas_ceilo(i_wavel)%bsc(jc,jk,jb) > 0.0_wp) THEN
!              p_art_data%diag%seas_att(i_wavel)%ceil_bsc(jc,jk,jb) = &
!                &  p_art_data%diag%seas_ceilo(i_wavel)%bsc(jc,jk,jb) * EXP(-2 * att_seas(jc,jk))
!            ELSE
!              p_art_data%diag%seas_att(i_wavel)%ceil_bsc(jc,jk,jb) = 0.0_wp
!            END IF
!          ENDDO !jc
!        ENDDO ! jk
!      ENDDO ! i_wavel
!
!      ! ... for satellite
!      DO i_wavel = 1, 3
!        DO jk = 1, nlev
!          jkm1 = MAX(1, jk-1)
!!NEC$ ivdep
!          DO jc = istart, iend
!            IF (jk == 1) att_seas(jc,jk) = 0.0_wp
!            ext_seas(jc,jk) = ( ext_seas_m1(i_wavel) * p_trac(jc,jk,iseasa) &
!                          &   + ext_seas_m2(i_wavel) * p_trac(jc,jk,iseasb) &
!                          &   + ext_seas_m3(i_wavel) * p_trac(jc,jk,iseasc) &
!                          &   ) * rho(jc,jk) * 1.e-6_wp
!            att_seas(jc,jk) = att_seas(jc,jkm1) + ext_seas(jc,jk) * dz(jc,jk)
!            IF (p_art_data%diag%seas_ceilo(i_wavel)%bsc(jc,jk,jb) == 0.0_wp) THEN
!              p_art_data%diag%seas_sat(i_wavel)%sat_bsc(jc,jk,jb) = 0.0_wp
!            ELSE
!              p_art_data%diag%seas_sat(i_wavel)%sat_bsc(jc,jk,jb) = &
!                &  p_art_data%diag%seas_ceilo(i_wavel)%bsc(jc,jk,jb) * EXP(-2 * att_seas(jc,jk))
!            END IF
!          ENDDO !jc
!        ENDDO !jk
!      ENDDO ! i_wavel
!
!      DEALLOCATE(ext_seas)
!      DEALLOCATE(att_seas)
!    ENDIF !iart_seasalt >0

    ! Calculate backscatter at specific wavelengths: VOLC

    IF (iasha > 0 .AND. iashb > 0 .AND. iashc > 0) THEN
      !!basalt
      bsc_volc_m1 = (/ 0.06846_wp, 0.09828_wp, 0.11124_wp /)
      ext_volc_m1 = (/ 0.84682_wp,  0.89899_wp,  1.22174_wp /)
      bsc_volc_m2 = (/ 0.01736_wp, 0.02019_wp, 0.03105_wp /)
      ext_volc_m2 = (/ 0.31802_wp,  0.32237_wp,  0.35050_wp /)
      bsc_volc_m3 = (/ 0.00553_wp, 0.00699_wp, 0.00924_wp /)
      ext_volc_m3 = (/ 0.13917_wp,  0.14501_wp,  0.15152_wp /)
      !!andesite
      !bsc_volc_m1 = (/ 0.07159_wp, 0.09492_wp, 0.09302_wp /)
      !ext_volc_m1 = (/ 0.85448_wp,  0.91301_wp,  1.22160_wp /)
      !bsc_volc_m2 = (/ 0.01579_wp, 0.02010_wp, 0.02964_wp /)
      !ext_volc_m2 = (/ 0.31621_wp,  0.32427_wp,  0.35591_wp /)
      !bsc_volc_m3 = (/ 0.00506_wp, 0.00575_wp, 0.00818_wp /)
      !ext_volc_m3 = (/ 0.14055_wp,  0.14414_wp,  0.14951_wp /)
      !!obsidian
      !bsc_volc_m1 = (/ 0.08343_wp, 0.09757_wp, 0.09385_wp /)
      !ext_volc_m1 = (/ 0.86564_wp,  0.89756_wp,  1.22527_wp /)
      !bsc_volc_m2 = (/ 0.02054_wp, 0.02539_wp, 0.03332_wp /)
      !ext_volc_m2 = (/ 0.31696_wp,  0.32982_wp,  0.35960_wp /)
      !bsc_volc_m3 = (/ 0.00740_wp, 0.00826_wp, 0.01033_wp /)
      !ext_volc_m3 = (/ 0.14218_wp,  0.14429_wp,  0.14983_wp /)
    
      ALLOCATE(bsc_volc(iend,nlev))
      DO i_wavel = 1, 3
        DO jk = 1, nlev
          DO jc = istart, iend
            bsc_volc(jc,jk) = ( bsc_volc_m1(i_wavel) * tracer(jc,jk,iasha) &
                          &   + bsc_volc_m2(i_wavel) * tracer(jc,jk,iashb) &
                          &   + bsc_volc_m3(i_wavel) * tracer(jc,jk,iashc) &
                          &   ) * rho(jc,jk) * 1.e-6_wp
            p_art_data%diag%volc_ceilo(i_wavel)%bsc(jc,jk,jb) = bsc_volc(jc,jk)
          ENDDO !jc
        ENDDO !jk
      ENDDO ! i_wavel
      DEALLOCATE(bsc_volc)
    
      !calculate attenuated backscatter
      ALLOCATE(att_volc(iend,nlev))
      ALLOCATE(ext_volc(iend,nlev))
      
      ! ... for ceilometer
      DO i_wavel = 1, 3
        DO jk = nlev, 1, -1
          jkp1 = MIN(nlev, jk+1)
!NEC$ ivdep
          DO jc = istart, iend
            IF (jk == nlev) att_volc(jc,jk) = 0.0_wp
            ext_volc(jc,jk) = ( ext_volc_m1(i_wavel) * tracer(jc,jk,iasha) &
                          &   + ext_volc_m2(i_wavel) * tracer(jc,jk,iashb) &
                          &   + ext_volc_m3(i_wavel) * tracer(jc,jk,iashc) &
                          &   ) * rho(jc,jk) * 1.e-6_wp
            att_volc(jc,jk) = att_volc(jc,jkp1) + ext_volc(jc,jk) * dz(jc,jk)
            IF (p_art_data%diag%volc_ceilo(i_wavel)%bsc(jc,jk,jb) > 0.0_wp) THEN
              p_art_data%diag%volc_att(i_wavel)%ceil_bsc(jc,jk,jb) = &
                &  p_art_data%diag%volc_ceilo(i_wavel)%bsc(jc,jk,jb) * EXP(-2 * att_volc(jc,jk))
            ELSE
              p_art_data%diag%volc_att(i_wavel)%ceil_bsc(jc,jk,jb) = 0.0_wp
            END IF
          ENDDO !jc
        END DO ! jk
      ENDDO ! i_wavel
      
      ! ... for satellite
      DO i_wavel = 1, 3
        DO jk = 1, nlev
          jkm1 = MAX(1, jk-1)
!NEC$ ivdep
          DO jc = istart, iend
            IF (jk == 1) att_volc(jc,jk) = 0.0_wp
              ext_volc(jc,jk) = ( ext_volc_m1(i_wavel) * tracer(jc,jk,iasha) &
                          &   + ext_volc_m2(i_wavel) * tracer(jc,jk,iashb) &
                          &   + ext_volc_m3(i_wavel) * tracer(jc,jk,iashc) &
                          &   ) * rho(jc,jk) * 1.e-6_wp
            att_volc(jc,jk) = att_volc(jc,jkm1) + ext_volc(jc,jk) * dz(jc,jk)
            IF (p_art_data%diag%volc_ceilo(i_wavel)%bsc(jc,jk,jb) == 0.0_wp) THEN
              p_art_data%diag%volc_sat(i_wavel)%sat_bsc(jc,jk,jb) = 0.0_wp
            ELSE
              p_art_data%diag%volc_sat(i_wavel)%sat_bsc(jc,jk,jb) = &
                &  p_art_data%diag%volc_ceilo(i_wavel)%bsc(jc,jk,jb) * EXP(-2 * att_volc(jc,jk))
            END IF
          ENDDO !jc
        ENDDO !jk
      ENDDO ! i_wavel
      
      DEALLOCATE(ext_volc)
      DEALLOCATE(att_volc)
    ENDIF

    IF (art_config(jg)%iart_fire > 0) THEN
    ! coated soot with cmd 150 nm
      bsc_soot_m1 = (/ 0.32933_wp,  0.28423_wp,  0.17841_wp /) 
      ext_soot_m1 = (/ 5.42035_wp,  4.47650_wp,  2.17680_wp /) 
      
      ALLOCATE(bsc_soot(iend,nlev))
      IF(isoot_insol_acc > 0) THEN
         DO i_wavel = 1, 3
            DO jk = 1, nlev
              DO jc = istart, iend
                bsc_soot(jc,jk) = bsc_soot_m1(i_wavel) * tracer(jc,jk,isoot_insol_acc) &
                              &   * rho(jc,jk) * 1.e-6_wp
                p_art_data%diag%soot_ceilo(i_wavel)%bsc(jc,jk,jb) = bsc_soot(jc,jk)
              ENDDO !jc
            ENDDO !jk
          ENDDO ! i_wavel
      ELSE
          DO i_wavel = 1, 3
            DO jk = 1, nlev
              DO jc = istart, iend
                bsc_soot(jc,jk) = bsc_soot_m1(i_wavel) * tracer(jc,jk,isoot) &
                              &   * rho(jc,jk) * 1.e-6_wp
                p_art_data%diag%soot_ceilo(i_wavel)%bsc(jc,jk,jb) = bsc_soot(jc,jk)
              ENDDO !jc
            ENDDO !jk
          ENDDO ! i_wavel
       ENDIF
      DEALLOCATE(bsc_soot)

      !calculate attenuated backscatter
      ALLOCATE(att_soot(iend,nlev))
      ALLOCATE(ext_soot(iend,nlev))

      ! ... for ceilometer
      DO i_wavel = 1, 3
        DO jk = nlev, 1, -1
          jkp1 = MIN(nlev, jk+1)
!NEC$ ivdep
          DO jc = istart, iend
            IF (jk == nlev) att_soot(jc,jk) = 0.0_wp
                IF(isoot_insol_acc > 0) THEN
                  ext_soot(jc,jk) = ( ext_soot_m1(i_wavel) * tracer(jc,jk,isoot_insol_acc) &
                              &   ) * rho(jc,jk) * 1.e-6_wp
                ELSE
                  ext_soot(jc,jk) = ( ext_soot_m1(i_wavel) * tracer(jc,jk,isoot) &
                              &   ) * rho(jc,jk) * 1.e-6_wp
                ENDIF
            att_soot(jc,jk) = att_soot(jc,jkp1) + ext_soot(jc,jk) * dz(jc,jk)
            IF (p_art_data%diag%soot_ceilo(i_wavel)%bsc(jc,jk,jb) > 0.0_wp) THEN
              p_art_data%diag%soot_att(i_wavel)%ceil_bsc(jc,jk,jb) = &
                &  p_art_data%diag%soot_ceilo(i_wavel)%bsc(jc,jk,jb) * EXP(-2 * att_soot(jc,jk))
            ELSE
              p_art_data%diag%soot_att(i_wavel)%ceil_bsc(jc,jk,jb) = 0.0_wp
            END IF
          ENDDO !jc
        ENDDO ! jk
      ENDDO ! i_wavel

      ! ... for satellite
      DO i_wavel = 1, 3
        DO jk = 1, nlev
          jkm1 = MAX(1, jk-1)
!NEC$ ivdep
          DO jc = istart, iend
            IF (jk == 1) att_soot(jc,jk) = 0.0_wp
              IF(isoot_insol_acc > 0) THEN
                ext_soot(jc,jk) = ( ext_soot_m1(i_wavel) * tracer(jc,jk,isoot_insol_acc) &
                          &   ) * rho(jc,jk) * 1.e-6_wp
              ELSE
                ext_soot(jc,jk) = ( ext_soot_m1(i_wavel) * tracer(jc,jk,isoot) &
                          &   ) * rho(jc,jk) * 1.e-6_wp
              ENDIF
            att_soot(jc,jk) = att_soot(jc,jkm1) + ext_soot(jc,jk) * dz(jc,jk)
            IF (p_art_data%diag%soot_ceilo(i_wavel)%bsc(jc,jk,jb) == 0.0_wp) THEN
              p_art_data%diag%soot_sat(i_wavel)%sat_bsc(jc,jk,jb) = 0.0_wp
            ELSE
              p_art_data%diag%soot_sat(i_wavel)%sat_bsc(jc,jk,jb) = &
                &  p_art_data%diag%soot_ceilo(i_wavel)%bsc(jc,jk,jb) * EXP(-2 * att_soot(jc,jk))
            END IF
          ENDDO !jc
        ENDDO !jk
      ENDDO ! i_wavel

      DEALLOCATE(ext_soot)
      DEALLOCATE(att_soot)
    ENDIF !iart_fire >0

  ENDIF !lart_diag_out

END SUBROUTINE art_calc_bsc
END MODULE mo_art_aero_optical_props
