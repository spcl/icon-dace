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

! Subroutine interface_aes_art calls the art reaction interface.

MODULE mo_interface_aes_art

  USE mo_model_domain           ,ONLY: t_patch

  USE mo_kind                   ,ONLY: wp
  USE mtime                     ,ONLY: t_datetime => datetime

  USE mo_aes_phy_config         ,ONLY: aes_phy_tc
  USE mo_aes_phy_memory         ,ONLY: t_aes_phy_field, prm_field

  USE mo_nonhydro_state         ,ONLY: p_nh_state_lists
  USE mo_dynamics_config        ,ONLY: nnew_rcf

#ifdef __ICON_ART
  USE mo_art_reaction_interface ,ONLY: art_reaction_interface
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_aes_art

CONTAINS

  SUBROUTINE interface_aes_art(patch)

    ! Arguments
    !
    TYPE(t_patch), TARGET, INTENT(in) :: patch

    ! Pointers
    !
    TYPE(t_aes_phy_field), POINTER    :: field

    ! Local variables
    !
    TYPE(t_datetime), POINTER :: datetime
    REAL(wp) :: pdtime
    LOGICAL  :: is_in_sd_ed_interval
    LOGICAL  :: is_active
    !
    INTEGER  :: jg

    jg  = patch%id

    datetime             => aes_phy_tc(jg)%datetime
    pdtime               =  aes_phy_tc(jg)%dt_phy_sec
    is_in_sd_ed_interval =  aes_phy_tc(jg)%is_in_sd_ed_interval_art
    is_active            =  aes_phy_tc(jg)%is_active_art

    ! associate pointers
    field     => prm_field(jg)

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN
          !
#ifdef __ICON_ART
          CALL art_reaction_interface(jg,                                           & !> in
               &                      datetime,                                     & !> in
               &                      pdtime,                                       & !> in
               &                      p_nh_state_lists(jg)%prog_list(nnew_rcf(jg)), & !> in
               &                      field%qtrc_phy)
#endif
          !
       END IF
       !
    END IF

    ! disassociate pointers
    NULLIFY(datetime, field)

  END SUBROUTINE interface_aes_art

END MODULE mo_interface_aes_art
