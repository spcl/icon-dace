!
! mo_art_external_init_soilhwsd
! This module provides initialization procedures for data from
! HWSD soil needed by the mineral dust emission routines
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

MODULE mo_art_external_init_soilhwsd
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_ext_data_types,                ONLY: t_external_data
  USE mo_impl_constants,                ONLY: min_rlcell
  USE mo_loopindices,                   ONLY: get_indices_c
! ART
  USE mo_art_external_types,            ONLY: t_art_soil_table, t_art_soil_properties

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_extinit_hwsdsoil_prepare
  PUBLIC :: art_extinit_hwsdsoil_finalize

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_extinit_hwsdsoil_prepare(dims2d, soil_prop)
!<
! SUBROUTINE art_extinit_hwsdsoil_prepare
! - Allocation of HWSD soil fields
! - Setting size distributions from Shao et al. (2010)
! Based on: Rieger (2016) - Der Einfluss von natuerlichem Aerosol auf 
!                           Wolken ueber Mitteleuropa, 
!                           Dissertation an der Fakultaet fuer Physik des 
!                           Karlsruher Instituts fuer Technologie (KIT)
! Part of Module: mo_art_external_init_soilhwsd
! Author: Daniel Rieger, KIT
! Initial Release: 2014-05-14
! Modifications:
! 2016-12-14: Daniel Rieger, KIT
! - Refactoring of code: Splitting into prepare and finalize SR
!>
  INTEGER, INTENT(IN)                       :: &
    &  dims2d(2)                                 !< 2D-Dimensions used to allocate fields
  TYPE(t_art_soil_properties),INTENT(inout) :: &
    &  soil_prop
! Local variables
  TYPE(t_art_soil_table),POINTER            :: &
    &  this_soil_table                           !< pointer to soil table
  INTEGER                                   :: &
    &  js                                        !< loop counter for soil types

  soil_prop%nsoil_types = 14
  ALLOCATE(soil_prop%soil_type(soil_prop%nsoil_types) )
  ALLOCATE(soil_prop%f_z0(dims2d(1),dims2d(2))        )
  ALLOCATE(soil_prop%wstrich(dims2d(1),dims2d(2))     )
  ALLOCATE(soil_prop%emiss_rate_a(dims2d(1),dims2d(2)))
  ALLOCATE(soil_prop%emiss_rate_b(dims2d(1),dims2d(2)))
  ALLOCATE(soil_prop%emiss_rate_c(dims2d(1),dims2d(2)))
  ALLOCATE(soil_prop%dust_mask(dims2d(1),dims2d(2))   )

  DO js = 1, soil_prop%nsoil_types

    this_soil_table => soil_prop%soil_type(js)

    SELECT CASE(js)
      CASE(1) ! Heavy Clay
        this_soil_table = soil_table('hcla', 4,                 & !< Name, nsoil_modes
          &        (/3.5542_wp,4.2239_wp,5.1638_wp,0.0000_wp/), & !< Median diameter min. disp.
          &        (/0.6931_wp,3.9323_wp,5.4486_wp,0.0000_wp/), & !< Median diameter ful. disp.
          &        (/1.0000_wp,0.2507_wp,0.4632_wp,0.0000_wp/), & !< Std. deviation min. disp.
          &        (/1.0000_wp,0.9181_wp,0.3916_wp,0.0000_wp/), & !< Std. deviation ful. disp.
          &        (/0.3902_wp,0.2813_wp,0.3286_wp,0.0000_wp/), & !< Massfraction min. disp.
          &        (/0.0872_wp,0.4464_wp,0.4665_wp,0.0000_wp/))   !< Massfraction ful. disp.
      CASE(2) ! Silty Clay
        this_soil_table = soil_table('silc', 4,                 & !< Name, nsoil_modes
          &        (/3.5542_wp,4.2239_wp,5.1638_wp,0.0000_wp/), & !< Median diameter min. disp.
          &        (/0.6931_wp,3.9323_wp,5.4486_wp,0.0000_wp/), & !< Median diameter ful. disp.
          &        (/1.0000_wp,0.2507_wp,0.4632_wp,0.0000_wp/), & !< Std. deviation min. disp.
          &        (/1.0000_wp,0.9181_wp,0.3916_wp,0.0000_wp/), & !< Std. deviation ful. disp.
          &        (/0.3902_wp,0.2813_wp,0.3286_wp,0.0000_wp/), & !< Massfraction min. disp.
          &        (/0.0872_wp,0.4464_wp,0.4665_wp,0.0000_wp/))   !< Massfraction ful. disp.
      CASE(3) ! Light Clay (In Shao this is called just 'clay')
        this_soil_table = soil_table('lcla', 4,                 & !< Name, nsoil_modes
          &        (/3.5542_wp,4.2239_wp,5.1638_wp,0.0000_wp/), & !< Median diameter min. disp.
          &        (/0.6931_wp,3.9323_wp,5.4486_wp,0.0000_wp/), & !< Median diameter ful. disp.
          &        (/1.0000_wp,0.2507_wp,0.4632_wp,0.0000_wp/), & !< Std. deviation min. disp.
          &        (/1.0000_wp,0.9181_wp,0.3916_wp,0.0000_wp/), & !< Std. deviation ful. disp.
          &        (/0.3902_wp,0.2813_wp,0.3286_wp,0.0000_wp/), & !< Massfraction min. disp.
          &        (/0.0872_wp,0.4464_wp,0.4665_wp,0.0000_wp/))   !< Massfraction ful. disp.
      CASE(4) ! Silty Clay Loam
        this_soil_table = soil_table('sicl', 4,                 & !< Name, nsoil_modes
          &        (/4.3565_wp,5.1674_wp,5.4092_wp,0.0000_wp/), & !< Median diameter min. disp.
          &        (/4.6079_wp,5.2050_wp,7.0553_wp,0.6931_wp/), & !< Median diameter ful. disp.
          &        (/0.4257_wp,0.3824_wp,1.0000_wp,0.0000_wp/), & !< Std. deviation min. disp.
          &        (/0.6141_wp,0.2897_wp,1.0000_wp,1.0000_wp/), & !< Std. deviation ful. disp.
          &        (/0.1114_wp,0.4554_wp,0.4331_wp,0.0000_wp/), & !< Massfraction min. disp.
          &        (/0.5844_wp,0.3304_wp,0.0522_wp,0.0330_wp/))   !< Massfraction ful. disp.
      CASE(5) ! Clay Loam
        this_soil_table = soil_table('cloa', 4,                 & !< Name, nsoil_modes
          &        (/4.3565_wp,5.1674_wp,5.4092_wp,0.0000_wp/), & !< Median diameter min. disp.
          &        (/4.6079_wp,5.2050_wp,7.0553_wp,0.6931_wp/), & !< Median diameter ful. disp.
          &        (/0.4257_wp,0.3824_wp,1.0000_wp,0.0000_wp/), & !< Std. deviation min. disp.
          &        (/0.6141_wp,0.2897_wp,1.0000_wp,1.0000_wp/), & !< Std. deviation ful. disp.
          &        (/0.1114_wp,0.4554_wp,0.4331_wp,0.0000_wp/), & !< Massfraction min. disp.
          &        (/0.5844_wp,0.3304_wp,0.0522_wp,0.0330_wp/))   !< Massfraction ful. disp.
      CASE(6) ! Silt --> equal to sand, as silt is not contained in shao data
        this_soil_table = soil_table('silt', 4,                 &
          &        (/0.0000_wp,4.3733_wp,5.7689_wp,0.0000_wp/), & !< Median diameter min. disp.
          &        (/0.0000_wp,0.6931_wp,5.6300_wp,0.0000_wp/), & !< Median diameter ful. disp.
          &        (/0.0000_wp,0.8590_wp,0.2526_wp,0.0000_wp/), & !< Std. deviation min. disp.
          &        (/0.0000_wp,1.0000_wp,0.2542_wp,0.0000_wp/), & !< Std. deviation ful. disp.
          &        (/0.0000_wp,0.0329_wp,0.9671_wp,0.0000_wp/), & !< Massfraction min. disp.
          &        (/0.0000_wp,0.0338_wp,0.9662_wp,0.0000_wp/))   !< Massfraction ful. disp.
      CASE(7) ! Silt Loam
        this_soil_table = soil_table('silo', 4,                 & !< Name, nsoil_modes
          &        (/3.9110_wp,4.56446_wp,0.000_wp,0.0000_wp/), & !< Median diameter min. disp.
          &        (/1.9996_wp,4.3235_wp,6.0420_wp,0.0000_wp/), & !< Median diameter ful. disp.
          &        (/0.1488_wp,0.2610_wp,0.0000_wp,0.0000_wp/), & !< Std. deviation min. disp.
          &        (/0.8649_wp,0.4658_wp,0.3464_wp,0.0000_wp/), & !< Std. deviation ful. disp.
          &        (/0.1112_wp,0.8888_wp,0.0000_wp,0.0000_wp/), & !< Massfraction min. disp.
          &        (/0.0748_wp,0.8900_wp,0.0352_wp,0.0000_wp/))   !< Massfraction ful. disp.
      CASE(8) ! Sandy Clay
        this_soil_table = soil_table('scla', 4,                 & !< Name, nsoil_modes
          &        (/4.3565_wp,5.1674_wp,5.4092_wp,0.0000_wp/), & !< Median diameter min. disp.
          &        (/0.6931_wp,3.9323_wp,5.4486_wp,0.0000_wp/), & !< Median diameter ful. disp.
          &        (/0.4257_wp,0.3824_wp,1.0000_wp,0.0000_wp/), & !< Std. deviation min. disp.
          &        (/1.0000_wp,0.9181_wp,0.3916_wp,0.0000_wp/), & !< Std. deviation ful. disp.
          &        (/0.1114_wp,0.4554_wp,0.4331_wp,0.0000_wp/), & !< Massfraction min. disp.
          &        (/0.0872_wp,0.4464_wp,0.4665_wp,0.0000_wp/))   !< Massfraction ful. disp.
      CASE(9) ! Loam
        this_soil_table = soil_table('loam', 4,                 & !< Name, nsoil_modes
          &        (/4.3565_wp,5.1674_wp,5.4092_wp,0.0000_wp/), & !< Median diameter min. disp.
          &        (/4.6079_wp,5.2050_wp,7.0553_wp,0.6931_wp/), & !< Median diameter ful. disp.
          &        (/0.4257_wp,0.3824_wp,1.0000_wp,0.0000_wp/), & !< Std. deviation min. disp.
          &        (/0.6141_wp,0.2897_wp,1.0000_wp,1.0000_wp/), & !< Std. deviation ful. disp.
          &        (/0.1114_wp,0.4554_wp,0.4331_wp,0.0000_wp/), & !< Massfraction min. disp.
          &        (/0.5844_wp,0.3304_wp,0.0522_wp,0.0330_wp/))   !< Massfraction ful. disp.
      CASE(10) ! Sandy Clay Loam
        this_soil_table = soil_table('sclo', 4,                 & !< Name, nsoil_modes
          &        (/3.9110_wp,4.56446_wp,0.000_wp,0.0000_wp/), & !< Median diameter min. disp.
          &        (/1.9996_wp,4.3235_wp,6.0420_wp,0.0000_wp/), & !< Median diameter ful. disp.
          &        (/0.1488_wp,0.2610_wp,0.0000_wp,0.0000_wp/), & !< Std. deviation min. disp.
          &        (/0.8649_wp,0.4658_wp,0.3464_wp,0.0000_wp/), & !< Std. deviation ful. disp.
          &        (/0.1112_wp,0.8888_wp,0.0000_wp,0.0000_wp/), & !< Massfraction min. disp.
          &        (/0.0748_wp,0.8900_wp,0.0352_wp,0.0000_wp/))   !< Massfraction ful. disp.
      CASE(11) ! Sandy Loam
        this_soil_table = soil_table('sloa', 4,                 & !< Name, nsoil_modes
          &        (/2.2675_wp,4.9654_wp,5.5819_wp,0.0000_wp/), & !< Median diameter min. disp.
          &        (/1.8079_wp,4.2050_wp,5.6553_wp,0.0000_wp/), & !< Median diameter ful. disp.
          &        (/1.0000_wp,0.3496_wp,0.5893_wp,0.0000_wp/), & !< Std. deviation min. disp.
          &        (/0.6141_wp,0.2897_wp,1.0000_wp,0.0000_wp/), & !< Std. deviation ful. disp.
          &        (/0.0722_wp,0.6266_wp,0.3012_wp,0.0000_wp/), & !< Massfraction min. disp.
          &        (/0.2344_wp,0.3634_wp,0.4022_wp,0.0000_wp/))   !< Massfraction ful. disp.
      CASE(12) ! Loamy Sand
        this_soil_table = soil_table('lsan', 4,                 & !< Name, nsoil_modes
          &        (/0.0000_wp,4.3733_wp,5.7689_wp,0.0000_wp/), & !< Median diameter min. disp.
          &        (/0.0000_wp,0.6931_wp,5.6300_wp,0.0000_wp/), & !< Median diameter ful. disp.
          &        (/0.0000_wp,0.8590_wp,0.2526_wp,0.0000_wp/), & !< Std. deviation min. disp.
          &        (/0.0000_wp,1.0000_wp,0.2542_wp,0.0000_wp/), & !< Std. deviation ful. disp.
          &        (/0.0000_wp,0.0329_wp,0.9671_wp,0.0000_wp/), & !< Massfraction min. disp.
          &        (/0.0000_wp,0.0338_wp,0.9662_wp,0.0000_wp/))   !< Massfraction ful. disp.
      CASE(13) ! Sand
        this_soil_table = soil_table('sand', 4,                 & !< Name, nsoil_modes
          &        (/0.0000_wp,4.3733_wp,5.7689_wp,0.0000_wp/), & !< Median diameter min. disp.
          &        (/0.0000_wp,0.6931_wp,5.6300_wp,0.0000_wp/), & !< Median diameter ful. disp.
          &        (/0.0000_wp,0.8590_wp,0.2526_wp,0.0000_wp/), & !< Std. deviation min. disp.
          &        (/0.0000_wp,1.0000_wp,0.2542_wp,0.0000_wp/), & !< Std. deviation ful. disp.
          &        (/0.0000_wp,0.0329_wp,0.9671_wp,0.0000_wp/), & !< Massfraction min. disp.
          &        (/0.0000_wp,0.0338_wp,0.9662_wp,0.0000_wp/))   !< Massfraction ful. disp.
      CASE(14) ! Undefined // Water
        this_soil_table = soil_table('udef', 4,                 & !< Name, nsoil_modes
          &        (/0.0000_wp,0.0000_wp,0.0000_wp,0.0000_wp/), & !< Median diameter min. disp.
          &        (/0.0000_wp,0.0000_wp,0.0000_wp,0.0000_wp/), & !< Median diameter ful. disp.
          &        (/0.0000_wp,0.0000_wp,0.0000_wp,0.0000_wp/), & !< Std. deviation min. disp.
          &        (/0.0000_wp,0.0000_wp,0.0000_wp,0.0000_wp/), & !< Std. deviation ful. disp.
          &        (/0.0000_wp,0.0000_wp,0.0000_wp,0.0000_wp/), & !< Massfraction min. disp.
          &        (/0.0000_wp,0.0000_wp,0.0000_wp,0.0000_wp/))   !< Massfraction ful. disp.
    END SELECT
  ENDDO

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
  TYPE(t_art_soil_table) FUNCTION soil_table(sname,nsoil_modes,                     &
    &                                        diam_m,diam_f,std_m,std_f,frm_m,frm_f)
    CHARACTER(LEN=*),INTENT(in)          :: &
      &  sname
    INTEGER,INTENT(in)                   :: &
      &  nsoil_modes
    REAL(wp),INTENT(in)                  :: &
      &  diam_m(:),                         &
      &  diam_f(:),                         &
      &  std_m(:),                          &
      &  std_f(:),                          &
      &  frm_m(:),                          &
      &  frm_f(:)

    soil_table%nsoil_modes = nsoil_modes
    ALLOCATE(soil_table%diam_med(nsoil_modes, 2))
    ALLOCATE(soil_table%std_dev (nsoil_modes, 2))
    ALLOCATE(soil_table%fr_mode (nsoil_modes, 2))
    ALLOCATE(soil_table%fr_soil(dims2d(1),dims2d(2)))

    soil_table%sname         = TRIM(sname)
    soil_table%diam_med(:,1) = EXP(diam_m(:)) * 1.e-6_wp !< conversion ln(dp) with dp in [mu m] to dp in [m]
    soil_table%diam_med(:,2) = EXP(diam_f(:)) * 1.e-6_wp !< conversion ln(dp) with dp in [mu m] to dp in [m]
    soil_table%std_dev(:,1)  = EXP(SQRT(std_m(:)))       !< the Shao et al. table includes ln(sigma)**2
    soil_table%std_dev(:,2)  = EXP(SQRT(std_f(:)))       !< the Shao et al. table includes ln(sigma)**2
    soil_table%fr_mode(:,1)  = frm_m(:)
    soil_table%fr_mode(:,2)  = frm_f(:)
  
  END FUNCTION soil_table
!!
!!-------------------------------------------------------------------------
!!
END SUBROUTINE art_extinit_hwsdsoil_prepare
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_extinit_hwsdsoil_finalize(p_patch,ext_data,soil_prop)
!<
! SUBROUTINE art_extinit_hwsdsoil_finalize
! The subroutine calculates the cross sectional area of the
! particles on the ground (Vogel et al, 2006: Eq. 3.8, Denominator
! combinded with 3.9 and 3.11). The integral in the deniminator of
! Eq. 3.8 is solved analytically. For a descritpion of the analyti-
! cal solution, contact Daniel Rieger. Additionally, a dust emission 
! mask is created.
! Based on: Rieger (2016) - Der Einfluss von natuerlichem Aerosol auf 
!                           Wolken ueber Mitteleuropa, 
!                           Dissertation an der Fakultaet fuer Physik des 
!                           Karlsruher Instituts fuer Technologie (KIT)
! Part of Module: mo_art_external_init_soilhwsd
! Author: Daniel Rieger, KIT
! Initial Release: 2014-05-27
! Modifications:
! 2014-06-06: Daniel Rieger, KIT
! - Added creation of dust emission mask
! 2016-11-16: Daniel Rieger, KIT
! - Refactoring as a part of the creation of mo_art_external state
!>
  TYPE(t_patch),INTENT(in)                  :: &
    &  p_patch                                   !< Domain information
  TYPE(t_external_data), INTENT(in)         :: &
    &  ext_data                                  !< Atmosphere external data (ICON)
  TYPE(t_art_soil_properties),INTENT(inout) :: &
    &  soil_prop                                 !< Storage type for all soil properties
! Local variables
  TYPE(t_art_soil_table),POINTER            :: &
    &  this_soil_table                           !< pointer to soil table
  REAL(wp)                                  :: &
    &  rho_p,                                  & !< Bulk density of particles
    &  bulkdens,                               & !< Bulk soil density
    &  lambda,                                 & !< factor for roughness correction
    &  claycont                                  !< clay content (%)
  INTEGER                                   :: &
    &  js, jm,                                 & !< loop counter for soil types/soil modes
    &  jc, jb,                                 & !< Loop indizes
    &  i_startblk, i_endblk,                   & !< Index of start and end block 
    &  i_rlstart, i_rlend,                     & !< Start and end values of relaxation zone
    &  istart, iend,                           & !< Start and end values of jc loop
    &  i_nchdom                                  !< Number of child domains

  rho_p    = 2.65E3_wp
  bulkdens = 1.5E3_wp

! Calculation of crosssectional areas
  DO js = 1, soil_prop%nsoil_types
    this_soil_table => soil_prop%soil_type(js)
    this_soil_table%stot_min = 0.0_wp
    this_soil_table%stot_ful = 0.0_wp
    DO jm = 1, this_soil_table%nsoil_modes
      IF(this_soil_table%fr_mode(jm,1) > 1.e-5_wp) THEN
! Area for minimal disperged (Vogel et al, 2006: Eq. 3.8, Denominator combined with 3.9 and 3.11)
        this_soil_table%stot_min = this_soil_table%stot_min                             &
          &                      + ( 3._wp / ( 2._wp * rho_p )                          &
          &                      * this_soil_table%fr_mode(jm,1) * bulkdens             &
          &                      / this_soil_table%diam_med(jm,1)                       &
          &                      * EXP(0.5_wp*LOG(this_soil_table%std_dev(jm,1))**2))
      ENDIF
      IF(this_soil_table%fr_mode(jm,2) > 1.e-5_wp) THEN
! Area for fully disperged (Vogel et al, 2006: Eq. 3.8, Denominator combined with 3.9 and 3.11)
        this_soil_table%stot_ful = this_soil_table%stot_ful                             &
          &                      + ( 3._wp / ( 2._wp * rho_p )                          &
          &                      * this_soil_table%fr_mode(jm,2) * bulkdens             &
          &                      / this_soil_table%diam_med(jm,2)                       &
          &                      * EXP(0.5_wp*LOG(this_soil_table%std_dev(jm,2))**2))
      ENDIF
    ENDDO
  ENDDO

  ! --- Get the loop indizes
  i_nchdom   = MAX(1,p_patch%n_childdom)
  i_rlstart  = 1
  i_rlend    = min_rlcell
  i_startblk = p_patch%cells%start_blk(i_rlstart,1)
  i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)
! Calculate roughness length for smooth conditions in [m]
  DO jb = i_startblk, i_endblk
    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
      &                istart, iend, i_rlstart, i_rlend)
    DO jc=istart, iend

! Roughness correction after Raupach et al. (1993)
      lambda = - 0.35_wp* LOG(1._wp - ext_data%atm%plcov(jc,jb))
      IF (lambda >= 2._wp) THEN
        soil_prop%f_z0(jc,jb) = 999._wp
      ELSE
        soil_prop%f_z0(jc,jb) = MAX(1._wp,(SQRT(1._wp - 0.5_wp * lambda) &
          &                   * SQRT(1._wp + 0.5_wp * 90._wp * lambda)))
      ENDIF

! Calculate Clay Content for Soil Water Content Correction
      claycont = 0._wp
      DO js = 1,soil_prop%nsoil_types
        this_soil_table => soil_prop%soil_type(js)
        ! Values obtained from USDA soil textural triangle (see Dissertation Rieger, 2016)
        IF (TRIM(this_soil_table%sname) == 'hcla') &
          &  claycont = claycont + this_soil_table%fr_soil(jc,jb) * 100._wp
        IF (TRIM(this_soil_table%sname) == 'lcla') &
          &  claycont = claycont + this_soil_table%fr_soil(jc,jb) * 80._wp
        IF (TRIM(this_soil_table%sname) == 'silc') &
          &  claycont = claycont + this_soil_table%fr_soil(jc,jb) * 50._wp
        IF (TRIM(this_soil_table%sname) == 'sicl') &
          &  claycont = claycont + this_soil_table%fr_soil(jc,jb) * 30._wp
        IF (TRIM(this_soil_table%sname) == 'cloa') &
          &  claycont = claycont + this_soil_table%fr_soil(jc,jb) * 30._wp
        IF (TRIM(this_soil_table%sname) == 'scla') &
          &  claycont = claycont + this_soil_table%fr_soil(jc,jb) * 45._wp
        IF (TRIM(this_soil_table%sname) == 'sclo') &
          &  claycont = claycont + this_soil_table%fr_soil(jc,jb) * 30._wp
        IF (TRIM(this_soil_table%sname) == 'loam') &
          &  claycont = claycont + this_soil_table%fr_soil(jc,jb) * 15._wp
        IF (TRIM(this_soil_table%sname) == 'silo') &
          &  claycont = claycont + this_soil_table%fr_soil(jc,jb) * 10._wp
        IF (TRIM(this_soil_table%sname) == 'sloa') &
          &  claycont = claycont + this_soil_table%fr_soil(jc,jb) * 10._wp 
        IF (TRIM(this_soil_table%sname) == 'silt') &
          &  claycont = claycont + this_soil_table%fr_soil(jc,jb) * 5._wp
        IF (TRIM(this_soil_table%sname) == 'lsan') &
          &  claycont = claycont + this_soil_table%fr_soil(jc,jb) * 5._wp
        IF (TRIM(this_soil_table%sname) == 'sand') &
          &  claycont = claycont + this_soil_table%fr_soil(jc,jb) * 5._wp
      ENDDO!js
      soil_prop%wstrich(jc,jb)  = 5._wp * (0.0014_wp*claycont*claycont+0.17_wp*claycont)

! Create dust mask
      soil_prop%dust_mask(jc,jb) = .FALSE.
      IF (ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_shrub_eg) > 0.1_wp) THEN
        ! Shrub cover Grassland/Forest // Evergreen
        soil_prop%dust_mask(jc,jb) = .TRUE.
      ENDIF
      IF (ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_shrub) > 0.1_wp) THEN
        ! Closed to open Shrubland (deciduous).
        soil_prop%dust_mask(jc,jb) = .TRUE.
      ENDIF
      IF (ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_grass) > 0.1_wp) THEN
        ! Grassland//herbaceous
        soil_prop%dust_mask(jc,jb) = .TRUE.
      ENDIF
      IF (ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_bare_soil) > 0.1_wp) THEN
        ! Bare soil
        soil_prop%dust_mask(jc,jb) = .TRUE.
      ENDIF
      IF (ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_sparse) > 0.1_wp) THEN
        ! Sparse shrub + sparse herbs
        soil_prop%dust_mask(jc,jb) = .TRUE.
      ENDIF
      IF (ext_data%atm%soiltyp(jc,jb) == 1) THEN ! ... no ice and glacier
        soil_prop%dust_mask(jc,jb) = .FALSE.
      ENDIF
      IF (ext_data%atm%soiltyp(jc,jb) == 2) THEN ! ... no rock, no lithosols
        soil_prop%dust_mask(jc,jb) = .FALSE.
      ENDIF
      IF (ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_urban) > 0.01_wp) THEN ! no urban fraction
        soil_prop%dust_mask(jc,jb) = .FALSE.
      ENDIF
      IF (ext_data%atm%fr_land(jc,jb) < 0.9_wp) THEN ! ... no water
        soil_prop%dust_mask(jc,jb) = .FALSE.
      ENDIF

    ENDDO!jc
  ENDDO!jb

END SUBROUTINE art_extinit_hwsdsoil_finalize
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_external_init_soilhwsd
