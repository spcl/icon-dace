!
! Provides interface to ART-routines dealing with washout
!
! This module provides an interface to the ART-routines.
! The interface is written in such a way, that ICON will compile and run
! properly, even if the ART-routines are not available at compile time.
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

MODULE mo_art_radiation_interface

  USE mo_kind,                          ONLY: wp
  USE mo_run_config,                    ONLY: lart
  USE mo_timer,                         ONLY: timers_level, timer_start, timer_stop,   &
                                          &   timer_art, timer_art_radInt
  ! ART Routines
  USE mo_art_radiation_aero,            ONLY: art_radiation_aero
  USE mo_art_config,                    ONLY: art_config

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_rad_aero_interface

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_rad_aero_interface(zaeq1,zaeq2,zaeq3,zaeq4,zaeq5,    &
        &                         zaea_rrtm,zaes_rrtm,zaeg_rrtm,    &
        &                         jg,jb,ks,ke,jcs,jce,nlong,nshort, &
        &                         aer_tau_lw_vr,                    &
        &                         aer_tau_sw_vr,                    &
        &                         aer_piz_sw_vr,                    &
        &                         aer_cg_sw_vr)
  REAL(wp),INTENT(in) ::     & !< Note that these values are only valid for Tegen climatology!
    &  zaeq1(:,:),           & !< continental aerosol
    &  zaeq2(:,:),           & !< maritime aerosol
    &  zaeq3(:,:),           & !< urban aerosol
    &  zaeq4(:,:),           & !< mineral dust aerosol
    &  zaeq5(:,:),           & !< stratospheric background aerosol
    &  zaea_rrtm(:,:),       & !< absorption coefficient
    &  zaes_rrtm(:,:),       & !< scatter coefficient
    &  zaeg_rrtm(:,:)          !< asymmetry coefficient
  INTEGER,INTENT(in) ::      &
    &  jg, jb,               & !< domain ID, block index
    &  ks, ke,               & !< loop index jk
    &  jcs, jce,             & !< loop index jc
    &  nlong,nshort            !< number of bands long/shortwave. Sorting of arrays: longwave 
                               !< first, then shortwave bands
  REAL(wp), OPTIONAL, INTENT(out) ::   &
    &  aer_tau_lw_vr(:,:,:), & !< longwave aerosol optical depth [layer-1], vertically reverse
    &  aer_tau_sw_vr(:,:,:), & !< shortwave aerosol optical depth [layer-1], vertically reverse
    &  aer_piz_sw_vr(:,:,:), & !< shortwave aerosol single scattering albedo [layer-1], 
                               !  vertically reverse
    &  aer_cg_sw_vr(:,:,:)     !< shortwave aerosol asymmetry factor [layer-1], vertically reverse
  
  if (lart) then
    IF (timers_level > 3) CALL timer_start(timer_art)
    IF (timers_level > 3) CALL timer_start(timer_art_radInt)

    if(art_config(jg)%iart_ari > 0) then
      CALL art_radiation_aero(zaeq1,zaeq2,zaeq3,zaeq4,zaeq5,    &
            &                 zaea_rrtm,zaes_rrtm,zaeg_rrtm,    &
            &                 jg,jb,ks,ke,jcs,jce,nlong,nshort, &
            &                 aer_tau_lw_vr,                    &
            &                 aer_tau_sw_vr,                    &
            &                 aer_piz_sw_vr,                    &
            &                 aer_cg_sw_vr)
    end if

    IF (timers_level > 3) CALL timer_stop(timer_art_radInt)
    IF (timers_level > 3) CALL timer_stop(timer_art)
  end if

END SUBROUTINE art_rad_aero_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_radiation_interface

