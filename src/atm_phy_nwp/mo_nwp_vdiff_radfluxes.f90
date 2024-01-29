!
! Surface radiation fluxes for the NWP VDIFF interface.
!
! This module is separate because it is required by the interface and sea model. Due to the
! dependence on mo_nwp_phy_types it cannot be included with mo_nwp_vdiff_types because that would
! introduce a cycle.
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

MODULE mo_nwp_vdiff_radfluxes

  USE mo_impl_constants, ONLY: end_prog_cells, start_prog_cells
  USE mo_kind, ONLY: wp
  USE mo_loopindices, ONLY: get_indices_c
  USE mo_model_domain, ONLY: t_patch
  USE mo_nwp_phy_types, ONLY: t_nwp_phy_diag

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_nwp_vdiff_surface_rad_fluxes

  !> Collection of surface radiation fluxes and diffuse fractions.
  TYPE t_nwp_vdiff_surface_rad_fluxes
    !> Long-wave downward radiation flux at surface [W/m**2] (nproma,nblks_c).
    REAL(wp), ALLOCATABLE :: flx_lw_down(:,:)
    !> Visible downward radiation flux at surface [W/m**2] (nproma,nblks_c).
    REAL(wp), CONTIGUOUS, POINTER :: flx_vis_down(:,:)
    !> Near-IR downward radiation flux at surface [W/m**2] (nproma,nblks_c).
    REAL(wp), CONTIGUOUS, POINTER :: flx_nir_down(:,:)
    !> PAR downward radiation flux at surface [W/m**2] (nproma,nblks_c).
    REAL(wp), CONTIGUOUS, POINTER :: flx_par_down(:,:)
    !> Diffuse fraction of downward visible flux [1] (nproma,nblks_c).
    REAL(wp), CONTIGUOUS, POINTER :: fr_vis_diffuse(:,:)
    !> Diffuse fraction of downward near-IR flux [1] (nproma,nblks_c).
    REAL(wp), CONTIGUOUS, POINTER :: fr_nir_diffuse(:,:)
    !> Diffuse fraction of downward PAR flux [1] (nproma,nblks_c).
    REAL(wp), CONTIGUOUS, POINTER :: fr_par_diffuse(:,:)
  CONTAINS
    PROCEDURE :: init => nwp_vdiff_surface_rad_fluxes_init
    PROCEDURE :: get_albdif => nwp_vdiff_surface_rad_fluxes_get_albdif
    FINAL :: nwp_vdiff_surface_rad_fluxes_destroy
  END TYPE t_nwp_vdiff_surface_rad_fluxes

CONTAINS

  !> Initialize a `t_nwp_vdiff_surface_rad_fluxes` structure and compute downward radiation fluxes
  !! at surface from fluxes provided in `phy_diag`. This routine is OpenMP orphaned.
  SUBROUTINE nwp_vdiff_surface_rad_fluxes_init (self, patch, phy_diag)

    CLASS(t_nwp_vdiff_surface_rad_fluxes), INTENT(OUT) :: self !< Object to initialize.

    TYPE(t_patch), INTENT(IN) :: patch !< Current patch.

    !> Diagnostic NWP physics variables.
    TYPE(t_nwp_phy_diag), INTENT(IN) :: phy_diag

    INTEGER :: i_startblk, i_endblk, i_blk
    INTEGER :: ics, ice, ic

    !$OMP SINGLE
      ALLOCATE(self%flx_lw_down(SIZE(phy_diag%lwflxsfc, DIM=1), SIZE(phy_diag%lwflxsfc, DIM=2)))
      self%flx_nir_down => phy_diag%swflx_nir_sfc
      self%flx_vis_down => phy_diag%swflx_vis_sfc
      self%flx_par_down => phy_diag%swflx_par_sfc
      self%fr_nir_diffuse => phy_diag%fr_nir_sfc_diff
      self%fr_vis_diffuse => phy_diag%fr_vis_sfc_diff
      self%fr_par_diffuse => phy_diag%fr_par_sfc_diff
    !$OMP END SINGLE

    !$ACC UPDATE DEVICE(self) ASYNC(1)

    !$ACC ENTER DATA ASYNC(1) CREATE(self%flx_lw_down) &
    !$ACC   ATTACH(self%flx_nir_down, self%flx_vis_down, self%flx_par_down) &
    !$ACC   ATTACH(self%fr_nir_diffuse, self%fr_vis_diffuse, self%fr_par_diffuse)

    i_startblk = patch%cells%start_block(start_prog_cells)
    i_endblk = patch%cells%end_block(end_prog_cells)

    !$OMP DO
    DO i_blk = i_startblk, i_endblk
      CALL get_indices_c(patch, i_blk, i_startblk, i_endblk, ics, ice, start_prog_cells, &
          & end_prog_cells)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO ic = ics, ice
        self%flx_lw_down(ic,i_blk) = phy_diag%lwflxsfc(ic,i_blk) + phy_diag%lwflx_up_sfc(ic,i_blk)
      END DO
      !$ACC END PARALLEL
    END DO

  END SUBROUTINE nwp_vdiff_surface_rad_fluxes_init

  !> Destroy a `t_surface_rad_fluxes` structure.
  SUBROUTINE nwp_vdiff_surface_rad_fluxes_destroy (self)

    TYPE(t_nwp_vdiff_surface_rad_fluxes), INTENT(INOUT) :: self !< Object to destroy.

    IF (.NOT. ALLOCATED(self%flx_lw_down)) RETURN

    !$ACC WAIT

    !$ACC EXIT DATA DELETE(self%flx_lw_down) &
    !$ACC   DETACH(self%flx_nir_down, self%flx_vis_down, self%flx_par_down) &
    !$ACC   DETACH(self%fr_nir_diffuse, self%fr_vis_diffuse, self%fr_par_diffuse)

    DEALLOCATE(self%flx_lw_down)
    NULLIFY(self%flx_nir_down, self%flx_par_down, self%flx_vis_down, self%flx_vis_down, &
        & self%fr_nir_diffuse, self%fr_par_diffuse, self%fr_vis_diffuse)

  END SUBROUTINE nwp_vdiff_surface_rad_fluxes_destroy

  !> Compute average diffuse albedo from NIR and VIS albedos. The average albedo yields the same
  !! short-wave flux as the separate ones under the current downwelling NIR/VIS flux ratio.
  SUBROUTINE nwp_vdiff_surface_rad_fluxes_get_albdif (self, patch, albnirdif, albvisdif, albdif)
    CLASS(t_nwp_vdiff_surface_rad_fluxes), INTENT(IN) :: self !< Surface fluxes.
    TYPE(t_patch), INTENT(IN) :: patch !< Current patch.
    REAL(wp), INTENT(IN) :: albnirdif(:,:) !< diffuse NIR albedo.
    REAL(wp), INTENT(IN) :: albvisdif(:,:) !< diffuse VIS albedo.
    REAL(wp), INTENT(OUT) :: albdif(:,:) !< average diffuse albedo.

    !> Minimal reference flux for calculating albedos [W/m**2].
    REAL(wp), PARAMETER :: FLX_RAD_EPSILON = 1e-3_wp

    INTEGER :: i_startblk, i_endblk, i_blk
    INTEGER :: ics, ice, ic

    REAL(wp) :: nir_dif_down !< Diffuse downward NIR radiation.
    REAL(wp) :: vis_dif_down !< Diffuse downward VIS radiation.

    i_startblk = patch%cells%start_block(start_prog_cells)
    i_endblk = patch%cells%end_block(end_prog_cells)

    DO i_blk = i_startblk, i_endblk
      CALL get_indices_c(patch, i_blk, i_startblk, i_endblk, ics, ice, start_prog_cells, &
          & end_prog_cells)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR PRIVATE(nir_dif_down, vis_dif_down)
        DO ic = ics, ice
          nir_dif_down = self%flx_nir_down(ic,i_blk) * self%fr_nir_diffuse(ic,i_blk)
          vis_dif_down = self%flx_vis_down(ic,i_blk) * self%fr_vis_diffuse(ic,i_blk)

          albdif(ic,i_blk) = &
              & (nir_dif_down * albnirdif(ic,i_blk) + vis_dif_down * albvisdif(ic,i_blk)) &
              & / MAX(FLX_RAD_EPSILON, nir_dif_down + vis_dif_down)
        END DO
      !$ACC END PARALLEL
    END DO

  END SUBROUTINE nwp_vdiff_surface_rad_fluxes_get_albdif

END MODULE mo_nwp_vdiff_radfluxes
