!
! mo_art_2mom_prepare
! This module is an adapted version of the mo_2mom_prepare by A. Seifert
!
! NOTE/DISCLAIMER: The ART version of the 2MOM scheme does not grant
! exact reproducability of ICON 2MOM scheme results. There are differences
! due to clipping etc.
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

MODULE mo_art_2mom_prepare
!ICON
  USE mo_kind,                          ONLY: wp
  USE mo_2mom_mcrph_types,              ONLY: particle, atmosphere
  USE mo_run_config,                    ONLY: iqr, iqi, iqs, iqg, iqh, iqc, iqv,         &
    &                                         iqnr,iqni,iqns,iqng,iqnh,iqnc,ininact
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: art_prepare_2mom
  PUBLIC :: art_post_2mom
  
  CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_prepare_2mom(atmo, cloud, rain, ice, snow, graupel, hail, &
    &                        rho, rhocorr, rhocld, pres, w, tk, p_trac, ik_slice)
!<
! SUBROUTINE art_prepare_2mom
! Preparations that are needed before microphysics calculation
! Based on: -
! Part of Module: mo_art_2mom_prepare
! Author: Daniel Rieger, KIT
! Initial Release: 2015-12-01
! Modifications:
! YYYY-MM-D: <name>, <institution>
! - ...
!>
  TYPE(atmosphere), INTENT(inout) :: &
    &  atmo                            !< Atmospheric state
  CLASS(particle), INTENT(inout)  :: &
    &  cloud, rain, ice,             & !< Hydrometeor characteristics
    &  snow, graupel, hail             !< Hydrometeor characteristics
  REAL(wp), TARGET, INTENT(in) :: &
    &  rho(:,:), rhocorr(:,:),    & !< Density, density dependency of particle fall speed
    &  rhocld(:,:),               & !< Density dependency of particle fall speed for cloud droplets
    &  pres(:,:), w(:,:), tk(:,:)   !< Pressure, vertical velocity, temperature
  REAL(wp), TARGET, INTENT(inout) :: &
    &  p_trac(:,:,:)                   !< Tracer fields
  INTEGER, INTENT(in) ::       &
    &  ik_slice(4)               !< Loop indices,ICON: (/jc_start,jc_end,jk_start,jk_end/)
  ! Local variables
  INTEGER :: ii, kk
    
    
  ! ... Transformation of microphysics variables to densities
  DO kk = ik_slice(3), ik_slice(4)
    DO ii = ik_slice(1), ik_slice(2)
      ! ... concentrations --> number densities
      p_trac(ii,kk,iqnc) = rho(ii,kk) * p_trac(ii,kk,iqnc)
      p_trac(ii,kk,iqnr) = rho(ii,kk) * p_trac(ii,kk,iqnr)
      p_trac(ii,kk,iqni) = rho(ii,kk) * p_trac(ii,kk,iqni)
      p_trac(ii,kk,iqns) = rho(ii,kk) * p_trac(ii,kk,iqns)
      p_trac(ii,kk,iqng) = rho(ii,kk) * p_trac(ii,kk,iqng)
      p_trac(ii,kk,iqnh) = rho(ii,kk) * p_trac(ii,kk,iqnh)
      
      ! ... mixing ratios -> mass densities
      p_trac(ii,kk,iqv) = rho(ii,kk) * p_trac(ii,kk,iqv)
      p_trac(ii,kk,iqc) = rho(ii,kk) * p_trac(ii,kk,iqc)
      p_trac(ii,kk,iqr) = rho(ii,kk) * p_trac(ii,kk,iqr)
      p_trac(ii,kk,iqi) = rho(ii,kk) * p_trac(ii,kk,iqi)
      p_trac(ii,kk,iqs) = rho(ii,kk) * p_trac(ii,kk,iqs)
      p_trac(ii,kk,iqg) = rho(ii,kk) * p_trac(ii,kk,iqg)
      p_trac(ii,kk,iqh) = rho(ii,kk) * p_trac(ii,kk,iqh)
      
      ! others
      p_trac(ii,kk,ininact) = rho(ii,kk) * p_trac(ii,kk,ininact)
    END DO
  END DO
   
  ! set pointers
  atmo%w    => w
  atmo%T    => tk
  atmo%p    => pres
  atmo%qv   => p_trac(:,:,iqv)
  atmo%rho  => rho
  
  cloud%rho_v   => rhocld
  rain%rho_v    => rhocorr
  ice%rho_v     => rhocorr
  graupel%rho_v => rhocorr
  snow%rho_v    => rhocorr
  hail%rho_v    => rhocorr
  
  cloud%q   => p_trac(:,:,iqc)
  cloud%n   => p_trac(:,:,iqnc)
  rain%q    => p_trac(:,:,iqr)
  rain%n    => p_trac(:,:,iqnr)
  ice%q     => p_trac(:,:,iqi)
  ice%n     => p_trac(:,:,iqni)
  snow%q    => p_trac(:,:,iqs)
  snow%n    => p_trac(:,:,iqns)
  graupel%q => p_trac(:,:,iqg)
  graupel%n => p_trac(:,:,iqng)
  hail%q    => p_trac(:,:,iqh)
  hail%n    => p_trac(:,:,iqnh)
  
END SUBROUTINE art_prepare_2mom
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_post_2mom(atmo, cloud, rain, ice, snow, graupel, hail, &
  &                      rho_r, p_trac, ik_slice)
!<
! SUBROUTINE art_post_2mom
! Postprocessing that is needed after microphysics calculation
! Based on: -
! Part of Module: mo_art_2mom_prepare
! Author: Daniel Rieger, KIT
! Initial Release: 2015-12-01
! Modifications:
! YYYY-MM-D: <name>, <institution>
! - ...
!>
  TYPE(atmosphere), INTENT(inout) :: &
    &  atmo                     !< Atmospheric state
  CLASS(particle), INTENT(inout)  :: &
    &  cloud, rain, ice,             & !< Hydrometeor characteristics
    &  snow, graupel, hail             !< Hydrometeor characteristics
  REAL(wp), TARGET, INTENT(in)    :: &
    &  rho_r(:,:)                   !< 1/Density
  REAL(wp), TARGET, INTENT(inout) :: &
    &  p_trac(:,:,:)            !< Tracer fields
  INTEGER, INTENT(in) ::       &
    &  ik_slice(4)              !< Loop indices,ICON: (/jc_start,jc_end,jk_start,jk_end/)
  ! Local variables
  INTEGER :: ii, kk
  
  ! nullify pointers
  atmo%w   => null()
  atmo%T   => null()
  atmo%p   => null()
  atmo%qv  => null()
  atmo%rho => null()
  
  cloud%rho_v   => null()
  rain%rho_v    => null()
  ice%rho_v     => null()
  graupel%rho_v => null()
  snow%rho_v    => null()
  hail%rho_v    => null()
  
  cloud%q   => null()
  cloud%n   => null()
  rain%q    => null()
  rain%n    => null()
  ice%q     => null()
  ice%n     => null()
  snow%q    => null()
  snow%n    => null()
  graupel%q => null()
  graupel%n => null()
  hail%q    => null()
  hail%n    => null()
  
  ! ... Transformation of variables back to ICON standard variables
  DO kk = ik_slice(3), ik_slice(4)
    DO ii = ik_slice(1), ik_slice(2)
      ! ... from mass densities back to mixing ratios
      p_trac(ii,kk,iqv) = rho_r(ii,kk) * p_trac(ii,kk,iqv)
      p_trac(ii,kk,iqc) = rho_r(ii,kk) * p_trac(ii,kk,iqc)
      p_trac(ii,kk,iqr) = rho_r(ii,kk) * p_trac(ii,kk,iqr)
      p_trac(ii,kk,iqi) = rho_r(ii,kk) * p_trac(ii,kk,iqi)
      p_trac(ii,kk,iqs) = rho_r(ii,kk) * p_trac(ii,kk,iqs)
      p_trac(ii,kk,iqg) = rho_r(ii,kk) * p_trac(ii,kk,iqg)
      p_trac(ii,kk,iqh) = rho_r(ii,kk) * p_trac(ii,kk,iqh)
      ! ... number concentrations
      p_trac(ii,kk,iqnc) = rho_r(ii,kk) * p_trac(ii,kk,iqnc)
      p_trac(ii,kk,iqnr) = rho_r(ii,kk) * p_trac(ii,kk,iqnr)
      p_trac(ii,kk,iqni) = rho_r(ii,kk) * p_trac(ii,kk,iqni)
      p_trac(ii,kk,iqns) = rho_r(ii,kk) * p_trac(ii,kk,iqns)
      p_trac(ii,kk,iqng) = rho_r(ii,kk) * p_trac(ii,kk,iqng)
      p_trac(ii,kk,iqnh) = rho_r(ii,kk) * p_trac(ii,kk,iqnh)
      ! Others
      p_trac(ii,kk,ininact) = rho_r(ii,kk) * p_trac(ii,kk,ininact)
    ENDDO
  ENDDO
END SUBROUTINE art_post_2mom
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_2mom_prepare
