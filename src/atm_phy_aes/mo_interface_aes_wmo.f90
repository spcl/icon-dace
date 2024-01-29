!
! Subroutine aes_phy_main calls all the parameterization schemes
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

MODULE mo_interface_aes_wmo

  USE mo_aes_phy_dims        ,ONLY: aes_phy_dims
  USE mo_aes_phy_memory      ,ONLY: t_aes_phy_field, prm_field

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_wmo

  USE mo_tropopause          ,ONLY: WMO_tropopause

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_aes_wmo

CONTAINS

  SUBROUTINE interface_aes_wmo(jg, jb, jcs, jce)

    INTEGER, INTENT(in)     :: jg, jb, jcs, jce

    ! Local variables
    !
    INTEGER  :: nlev
    INTEGER  :: nproma
    !
    TYPE(t_aes_phy_field), POINTER :: field

    IF (ltimer) call timer_start(timer_wmo)

    nlev   = aes_phy_dims(jg)%nlev
    nproma = aes_phy_dims(jg)%nproma

    ! associate pointers
    field  => prm_field(jg)

    CALL WMO_tropopause( jg,                       &! in
                       & jcs, jce, nproma, nlev,   &! in
                       & field% ta(:,:,jb),        &! in
                       & field% pfull(:,:,jb),     &! in
                       & field% ptp(:,jb),         &! inout for diagnostics
                       & lacc = .TRUE.)
    !

    IF (ltimer) call timer_stop(timer_wmo)

  END SUBROUTINE interface_aes_wmo
  !-------------------------------------------------------------------

END MODULE mo_interface_aes_wmo
