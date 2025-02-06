!> QUINCY litter transport and partitioning
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### routines for litter transport and partitioning for the both JSM and SSM (soil models in QUINCY)
!>
MODULE mo_q_sb_litter_processes
#ifndef __NO_QUINCY__

  USE mo_kind,            ONLY: wp
  USE mo_lnd_bgcm_idx

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_litter_partitioning

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_sb_litter_processes'

CONTAINS

  ! ======================================================================================================= !
  !>Calculates the partitioning of litter fall from vegetation into metabolic, structural,
  !>and woody litter given
  !>  a) tissue specifc preferences (labile + reserve = metabolic; sapwood and heartwood = woody);
  !>  b) the tissue-specific lignin concentration and the actual tissue C:N ratio.
  !>All above-ground litter enters the first layer, only fine and coarse roots are distributed
  !>according to their root profile
  !>
  SUBROUTINE calc_litter_partitioning( &
    & nc, &                         ! in
    & nsoil_sb, &
    & num_sl_above_bedrock, &
    & lctlib_sla, &
    & lctlib_growthform, &
    & sb_model_scheme, &
    & soil_depth_sl, &
    & root_fraction_sl, &
    & veg_litterfall_mt, &             ! in
    & sb_formation_mt )                ! inout

    USE mo_veg_constants,           ONLY: ITREE, IGRASS
    USE mo_jsb_math_constants,      ONLY: eps8
    USE mo_sb_constants,            ONLY: fc_soluable_max, k_fc_soluable, k_fn_soluable, k_fp_soluable, &
                                          lc_leaf_max, lc_leaf2sla, lc_fine_root, lc_coarse_root, lc_fruit, lc_seed_bed
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                           INTENT(in)    :: nc                        !< dimension
    INTEGER,                           INTENT(in)    :: nsoil_sb                  !< nr soil layer
    REAL(wp), DIMENSION(nc),           INTENT(in)    :: num_sl_above_bedrock      !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    REAL(wp),                          INTENT(in)    :: lctlib_sla                !< lctlib parameter
    INTEGER,                           INTENT(in)    :: lctlib_growthform         !< lctlib parameter
    CHARACTER(len=*),                  INTENT(in)    :: sb_model_scheme           !< config option: soil model
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: soil_depth_sl             !< soil layer thickness
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: root_fraction_sl          !< root fraction per soil layer
    REAL(wp),                          INTENT(in)    :: veg_litterfall_mt(:,:,:)  !< bgcm flux: vegetation litterfall
    REAL(wp),                          INTENT(inout) :: sb_formation_mt(:,:,:,:)  !< bgcm flux: soil biogeochemistry formation flux
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), DIMENSION(nc) :: fc_soluable_leaf           !< fc: ...
    REAL(wp), DIMENSION(nc) :: fc_soluable_fine_root
    REAL(wp), DIMENSION(nc) :: fc_soluable_coarse_root
    REAL(wp), DIMENSION(nc) :: fc_soluable_fruit
    REAL(wp), DIMENSION(nc) :: fc_soluable_seed_bed
    REAL(wp), DIMENSION(nc) :: fc_soluable_labile
    REAL(wp), DIMENSION(nc) :: fc_soluable_reserve
    REAL(wp), DIMENSION(nc) :: fn_soluable_leaf           !< fn: ...
    REAL(wp), DIMENSION(nc) :: fn_soluable_fine_root
    REAL(wp), DIMENSION(nc) :: fn_soluable_coarse_root
    REAL(wp), DIMENSION(nc) :: fn_soluable_fruit
    REAL(wp), DIMENSION(nc) :: fn_soluable_seed_bed
    REAL(wp), DIMENSION(nc) :: fn_soluable_labile
    REAL(wp), DIMENSION(nc) :: fn_soluable_reserve
    REAL(wp), DIMENSION(nc) :: fp_soluable_leaf           !< fp: ...
    REAL(wp), DIMENSION(nc) :: fp_soluable_fine_root
    REAL(wp), DIMENSION(nc) :: fp_soluable_coarse_root
    REAL(wp), DIMENSION(nc) :: fp_soluable_fruit
    REAL(wp), DIMENSION(nc) :: fp_soluable_seed_bed
    REAL(wp), DIMENSION(nc) :: fp_soluable_labile
    REAL(wp), DIMENSION(nc) :: fp_soluable_reserve
    INTEGER                 :: ic, is                     !< loop over dimensions
    INTEGER                 :: isoil                      !< loop over soil layers
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_litter_partitioning'
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.9 init local variables
    !>  labile and reserve compartments are init with 1.0
    !>
    fc_soluable_leaf        = 0.0_wp
    fc_soluable_fine_root   = 0.0_wp
    fc_soluable_coarse_root = 0.0_wp
    fc_soluable_fruit       = 0.0_wp
    fc_soluable_seed_bed    = 0.0_wp
    fn_soluable_leaf        = 0.0_wp
    fn_soluable_fine_root   = 0.0_wp
    fn_soluable_coarse_root = 0.0_wp
    fn_soluable_fruit       = 0.0_wp
    fn_soluable_seed_bed    = 0.0_wp
    fp_soluable_leaf        = 0.0_wp
    fp_soluable_fine_root   = 0.0_wp
    fp_soluable_coarse_root = 0.0_wp
    fp_soluable_fruit       = 0.0_wp
    fp_soluable_seed_bed    = 0.0_wp
    ! labile & reserve
    fc_soluable_labile      = 1.0_wp
    fc_soluable_reserve     = 1.0_wp
    fn_soluable_labile      = 1.0_wp
    fn_soluable_reserve     = 1.0_wp
    fp_soluable_labile      = 1.0_wp
    fp_soluable_reserve     = 1.0_wp

    !>1.0 determine allocation of fresh litter fall to litter pools
    !>
    DO ic = 1,nc
      ! leaf
      IF(veg_litterfall_mt(ix_leaf, ixN, ic) > eps8) THEN
        fc_soluable_leaf(ic) = MAX(0.0_wp, fc_soluable_max - k_fc_soluable * (lc_leaf_max + lc_leaf2sla * lctlib_sla) &
          &                    * veg_litterfall_mt(ix_leaf, ixC, ic) / veg_litterfall_mt(ix_leaf, ixN, ic))
        IF(fc_soluable_leaf(ic) > eps8) THEN
          fn_soluable_leaf(ic) = 1._wp / (1._wp + (1._wp - fc_soluable_leaf(ic)) / (k_fn_soluable * fc_soluable_leaf(ic)))
          fp_soluable_leaf(ic) = 1._wp / (1._wp + (1._wp - fc_soluable_leaf(ic)) / (k_fp_soluable * fc_soluable_leaf(ic)))
        END IF
      END IF
      ! fine_root
      IF(veg_litterfall_mt(ix_fine_root, ixN, ic) > eps8) THEN
        fc_soluable_fine_root(ic) = MAX(0.0_wp, fc_soluable_max - k_fc_soluable * lc_fine_root &
          &                         * veg_litterfall_mt(ix_fine_root, ixC, ic) / veg_litterfall_mt(ix_fine_root, ixN, ic))
        IF(fc_soluable_fine_root(ic) > eps8) THEN
          fn_soluable_fine_root(ic) = 1._wp &
            &                         / (1._wp + (1._wp - fc_soluable_fine_root(ic)) / (k_fn_soluable * fc_soluable_fine_root(ic)))
          fp_soluable_fine_root(ic) = 1._wp &
            &                         / (1._wp + (1._wp - fc_soluable_fine_root(ic)) / (k_fp_soluable * fc_soluable_fine_root(ic)))
        END IF
      END IF
      ! coarse_root
      IF(veg_litterfall_mt(ix_coarse_root, ixN, ic) > eps8) THEN
        fc_soluable_coarse_root(ic) = MAX(0.0_wp, fc_soluable_max - k_fc_soluable * lc_coarse_root &
          &                           * veg_litterfall_mt(ix_coarse_root, ixC, ic) / veg_litterfall_mt(ix_coarse_root, ixN, ic))
        IF(fc_soluable_coarse_root(ic) > eps8) THEN
          fn_soluable_coarse_root(ic) = 1._wp / (1._wp + (1._wp - fc_soluable_coarse_root(ic)) &
            &                           / (k_fn_soluable*fc_soluable_coarse_root(ic)))
          fp_soluable_coarse_root(ic) = 1._wp / (1._wp + (1._wp - fc_soluable_coarse_root(ic)) &
            &                           / (k_fp_soluable*fc_soluable_coarse_root(ic)))
         END IF
      END IF
      ! fruit
      IF(veg_litterfall_mt(ix_fruit, ixN, ic) > eps8) THEN
        fc_soluable_fruit(ic) = MAX(0.0_wp, fc_soluable_max - k_fc_soluable * lc_fruit &
          &                     * veg_litterfall_mt(ix_fruit, ixC, ic) / veg_litterfall_mt(ix_fruit, ixN, ic))
        IF(fc_soluable_fruit(ic) > eps8) THEN
          fn_soluable_fruit(ic) = 1._wp / (1._wp + (1._wp - fc_soluable_fruit(ic)) / (k_fn_soluable * fc_soluable_fruit(ic)))
          fp_soluable_fruit(ic) = 1._wp / (1._wp + (1._wp - fc_soluable_fruit(ic)) / (k_fp_soluable * fc_soluable_fruit(ic)))
        END IF
      END IF
      ! seed_bed
      IF(veg_litterfall_mt(ix_seed_bed, ixN, ic) > eps8) THEN
        fc_soluable_seed_bed(ic) = MAX(0.0_wp, fc_soluable_max - k_fc_soluable * lc_seed_bed &
          &                        * veg_litterfall_mt(ix_seed_bed, ixC, ic) / veg_litterfall_mt(ix_seed_bed, ixN, ic))
        IF(fc_soluable_seed_bed(ic) > eps8) THEN
          fn_soluable_seed_bed(ic) = 1._wp &
            &                        / (1._wp + (1._wp - fc_soluable_seed_bed(ic)) / (k_fn_soluable * fc_soluable_seed_bed(ic)))
          fp_soluable_seed_bed(ic) = 1._wp &
            &                        / (1._wp + (1._wp - fc_soluable_seed_bed(ic)) / (k_fp_soluable * fc_soluable_seed_bed(ic)))
        END IF
      END IF
    END DO ! loop over nc

    !>2.0 Update surface metabolic and soluable litter pools with litter fall
    !>
    !>  2.1 soluable litter
    !>
    sb_formation_mt(ix_soluable_litter, ixC, :, 1) = &
      &                              sb_formation_mt(ix_soluable_litter, ixC, :, 1) &
      &                              + (fc_soluable_leaf(:)       * veg_litterfall_mt(ix_leaf, ixC, :)      &
      &                                 + fc_soluable_fruit(:)    * veg_litterfall_mt(ix_fruit, ixC, :)     &
      &                                 + fc_soluable_labile(:)   * veg_litterfall_mt(ix_labile, ixC, :)    &
      &                                 + fc_soluable_reserve(:)  * veg_litterfall_mt(ix_reserve, ixC, :)   &
      &                                 + fc_soluable_seed_bed(:) * veg_litterfall_mt(ix_seed_bed, ixC, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_soluable_litter, ixN, :, 1) = &
      &                              sb_formation_mt(ix_soluable_litter, ixN, :, 1) &
      &                              + (fn_soluable_leaf(:)       * veg_litterfall_mt(ix_leaf, ixN, :)      &
      &                                 + fn_soluable_fruit(:)    * veg_litterfall_mt(ix_fruit, ixN, :)     &
      &                                 + fn_soluable_labile(:)   * veg_litterfall_mt(ix_labile, ixN, :)    &
      &                                 + fn_soluable_reserve(:)  * veg_litterfall_mt(ix_reserve, ixN, :)   &
      &                                 + fn_soluable_seed_bed(:) * veg_litterfall_mt(ix_seed_bed, ixN, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_soluable_litter, ixP, :, 1) = &
      &                              sb_formation_mt(ix_soluable_litter, ixP, :, 1) &
      &                              + (fp_soluable_leaf(:)       * veg_litterfall_mt(ix_leaf, ixP, :)      &
      &                                 + fp_soluable_fruit(:)    * veg_litterfall_mt(ix_fruit, ixP, :)     &
      &                                 + fp_soluable_labile(:)   * veg_litterfall_mt(ix_labile, ixP, :)    &
      &                                 + fp_soluable_reserve(:)  * veg_litterfall_mt(ix_reserve, ixP, :)   &
      &                                 + fp_soluable_seed_bed(:) * veg_litterfall_mt(ix_seed_bed, ixP, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_soluable_litter, ixC13, :, 1) = &
      &                              sb_formation_mt(ix_soluable_litter, ixC13, :, 1) &
      &                              + (fc_soluable_leaf(:)       * veg_litterfall_mt(ix_leaf, ixC13, :)      &
      &                                 + fc_soluable_fruit(:)    * veg_litterfall_mt(ix_fruit, ixC13, :)     &
      &                                 + fc_soluable_labile(:)   * veg_litterfall_mt(ix_labile, ixC13, :)    &
      &                                 + fc_soluable_reserve(:)  * veg_litterfall_mt(ix_reserve, ixC13, :)   &
      &                                 + fc_soluable_seed_bed(:) * veg_litterfall_mt(ix_seed_bed, ixC13, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_soluable_litter, ixC14, :, 1) = &
      &                              sb_formation_mt(ix_soluable_litter, ixC14, :, 1) &
      &                              + (fc_soluable_leaf(:)       * veg_litterfall_mt(ix_leaf, ixC14, :)      &
      &                                 + fc_soluable_fruit(:)    * veg_litterfall_mt(ix_fruit, ixC14, :)     &
      &                                 + fc_soluable_labile(:)   * veg_litterfall_mt(ix_labile, ixC14, :)    &
      &                                 + fc_soluable_reserve(:)  * veg_litterfall_mt(ix_reserve, ixC14, :)   &
      &                                 + fc_soluable_seed_bed(:) * veg_litterfall_mt(ix_seed_bed, ixC14, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_soluable_litter, ixN15, :, 1) = &
      &                              sb_formation_mt(ix_soluable_litter, ixN15, :, 1) &
      &                              + (fn_soluable_leaf(:)       * veg_litterfall_mt(ix_leaf, ixN15, :)      &
      &                                 + fn_soluable_fruit(:)    * veg_litterfall_mt(ix_fruit, ixN15, :)     &
      &                                 + fn_soluable_labile(:)   * veg_litterfall_mt(ix_labile, ixN15, :)    &
      &                                 + fn_soluable_reserve(:)  * veg_litterfall_mt(ix_reserve, ixN15, :)   &
      &                                 + fn_soluable_seed_bed(:) * veg_litterfall_mt(ix_seed_bed, ixN15, :)) &
      &                              / soil_depth_sl(:,1)
    !>  2.2 polymeric litter
    !>
    sb_formation_mt(ix_polymeric_litter, ixC, :, 1) = &
      &                              sb_formation_mt(ix_polymeric_litter, ixC, :, 1) &
      &                              + ((1.0_wp - fc_soluable_leaf(:))       * veg_litterfall_mt(ix_leaf, ixC, :)      &
      &                                 + (1.0_wp - fc_soluable_fruit(:))    * veg_litterfall_mt(ix_fruit, ixC, :)     &
      &                                 + (1.0_wp - fc_soluable_labile(:))   * veg_litterfall_mt(ix_labile, ixC, :)    &
      &                                 + (1.0_wp - fc_soluable_reserve(:))  * veg_litterfall_mt(ix_reserve, ixC, :)   &
      &                                 + (1.0_wp - fc_soluable_seed_bed(:)) * veg_litterfall_mt(ix_seed_bed, ixC, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_polymeric_litter, ixN, :, 1) = &
      &                              sb_formation_mt(ix_polymeric_litter, ixN, :, 1) &
      &                              + ((1.0_wp - fn_soluable_leaf(:))       * veg_litterfall_mt(ix_leaf, ixN, :)      &
      &                                 + (1.0_wp - fn_soluable_fruit(:))    * veg_litterfall_mt(ix_fruit, ixN, :)     &
      &                                 + (1.0_wp - fn_soluable_labile(:))   * veg_litterfall_mt(ix_labile, ixN, :)    &
      &                                 + (1.0_wp - fn_soluable_reserve(:))  * veg_litterfall_mt(ix_reserve, ixN, :)   &
      &                                 + (1.0_wp - fn_soluable_seed_bed(:)) * veg_litterfall_mt(ix_seed_bed, ixN, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_polymeric_litter, ixP, :, 1) = &
      &                              sb_formation_mt(ix_polymeric_litter, ixP, :, 1) &
      &                              + ((1.0_wp - fp_soluable_leaf(:))       * veg_litterfall_mt(ix_leaf, ixP, :)      &
      &                                 + (1.0_wp - fp_soluable_fruit(:))    * veg_litterfall_mt(ix_fruit, ixP, :)     &
      &                                 + (1.0_wp - fp_soluable_labile(:))   * veg_litterfall_mt(ix_labile, ixP, :)    &
      &                                 + (1.0_wp - fp_soluable_reserve(:))  * veg_litterfall_mt(ix_reserve, ixP, :)   &
      &                                 + (1.0_wp - fp_soluable_seed_bed(:)) * veg_litterfall_mt(ix_seed_bed, ixP, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_polymeric_litter, ixC13, :, 1) = &
      &                              sb_formation_mt(ix_polymeric_litter, ixC13, :, 1) &
      &                              + ((1.0_wp - fc_soluable_leaf(:))       * veg_litterfall_mt(ix_leaf, ixC13, :)      &
      &                                 + (1.0_wp - fc_soluable_fruit(:))    * veg_litterfall_mt(ix_fruit, ixC13, :)     &
      &                                 + (1.0_wp - fc_soluable_labile(:))   * veg_litterfall_mt(ix_labile, ixC13, :)    &
      &                                 + (1.0_wp - fc_soluable_reserve(:))  * veg_litterfall_mt(ix_reserve, ixC13, :)   &
      &                                 + (1.0_wp - fc_soluable_seed_bed(:)) * veg_litterfall_mt(ix_seed_bed, ixC13, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_polymeric_litter, ixC14, :, 1) = &
      &                              sb_formation_mt(ix_polymeric_litter, ixC14, :, 1) &
      &                              + ((1.0_wp - fc_soluable_leaf(:))       * veg_litterfall_mt(ix_leaf, ixC14, :)      &
      &                                 + (1.0_wp - fc_soluable_fruit(:))    * veg_litterfall_mt(ix_fruit, ixC14, :)     &
      &                                 + (1.0_wp - fc_soluable_labile(:))   * veg_litterfall_mt(ix_labile, ixC14, :)    &
      &                                 + (1.0_wp - fc_soluable_reserve(:))  * veg_litterfall_mt(ix_reserve, ixC14, :)   &
      &                                 + (1.0_wp - fc_soluable_seed_bed(:)) * veg_litterfall_mt(ix_seed_bed, ixC14, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_polymeric_litter, ixN15, :, 1) = &
      &                              sb_formation_mt(ix_polymeric_litter, ixN15, :, 1) &
      &                              + ((1.0_wp - fn_soluable_leaf(:))       * veg_litterfall_mt(ix_leaf, ixN15, :)      &
      &                                 + (1.0_wp - fn_soluable_fruit(:))    * veg_litterfall_mt(ix_fruit, ixN15, :)     &
      &                                 + (1.0_wp - fn_soluable_labile(:))   * veg_litterfall_mt(ix_labile, ixN15, :)    &
      &                                 + (1.0_wp - fn_soluable_reserve(:))  * veg_litterfall_mt(ix_reserve, ixN15, :)   &
      &                                 + (1.0_wp - fn_soluable_seed_bed(:)) * veg_litterfall_mt(ix_seed_bed, ixN15, :)) &
      &                              / soil_depth_sl(:,1)

    !>3.0 Add surface litter fall from sapwood and heartwood pools to surface woody litter pool for trees
    !>    for grasses, stems and halms (=sapwood) are not actually wood, and therefore paritioned to
    !>    metabolic and structural litter as if their were fine roots
    SELECT CASE(lctlib_growthform)
    ! tree PFT
    CASE (ITREE)
      ! woody litter
      ! @TODO may be re-written to a one-liner using bgcm matrices
      sb_formation_mt(ix_woody_litter, ixC, :, 1) = &
        &           sb_formation_mt(ix_woody_litter, ixC, :, 1) &
        &           + (veg_litterfall_mt(ix_sap_wood, ixC, :) + veg_litterfall_mt(ix_heart_wood, ixC, :)) / soil_depth_sl(:,1)
      sb_formation_mt(ix_woody_litter, ixN, :, 1) = &
        &           sb_formation_mt(ix_woody_litter, ixN, :, 1) &
        &           + (veg_litterfall_mt(ix_sap_wood, ixN, :) + veg_litterfall_mt(ix_heart_wood, ixN, :)) / soil_depth_sl(:,1)
      sb_formation_mt(ix_woody_litter, ixP, :, 1) = &
        &           sb_formation_mt(ix_woody_litter, ixP, :, 1) &
        &           + (veg_litterfall_mt(ix_sap_wood, ixP, :) + veg_litterfall_mt(ix_heart_wood, ixP, :)) / soil_depth_sl(:,1)
      sb_formation_mt(ix_woody_litter, ixC13, :, 1) = &
        &           sb_formation_mt(ix_woody_litter, ixC13, :, 1) &
        &           + (veg_litterfall_mt(ix_sap_wood, ixC13, :) + veg_litterfall_mt(ix_heart_wood, ixC13, :)) / soil_depth_sl(:,1)
      sb_formation_mt(ix_woody_litter, ixC14, :, 1) = &
        &           sb_formation_mt(ix_woody_litter, ixC14, :, 1) &
        &           + (veg_litterfall_mt(ix_sap_wood, ixC14, :) + veg_litterfall_mt(ix_heart_wood, ixC14, :)) / soil_depth_sl(:,1)
      sb_formation_mt(ix_woody_litter, ixN15, :, 1) = &
        &           sb_formation_mt(ix_woody_litter, ixN15, :, 1) &
        &           + (veg_litterfall_mt(ix_sap_wood, ixN15, :) + veg_litterfall_mt(ix_heart_wood, ixN15, :)) / soil_depth_sl(:,1)
    ! grass PFT
    CASE (IGRASS)
      ! soluable litter
      sb_formation_mt(ix_soluable_litter, ixC, :, 1) = &
        &             sb_formation_mt(ix_soluable_litter, ixC, :, 1) &
        &             + fc_soluable_fine_root(:) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixC, :) + veg_litterfall_mt(ix_heart_wood, ixC, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_soluable_litter, ixN, :, 1) = &
        &             sb_formation_mt(ix_soluable_litter, ixN, :, 1) &
        &             + fn_soluable_fine_root(:) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixN, :) + veg_litterfall_mt(ix_heart_wood, ixN, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_soluable_litter, ixP, :, 1) = &
        &             sb_formation_mt(ix_soluable_litter, ixP, :, 1) &
        &             + fp_soluable_fine_root(:) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixP, :) + veg_litterfall_mt(ix_heart_wood, ixP, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_soluable_litter, ixC13, :, 1) = &
        &             sb_formation_mt(ix_soluable_litter, ixC13, :, 1) &
        &             + fc_soluable_fine_root(:) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixC13, :) + veg_litterfall_mt(ix_heart_wood, ixC13, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_soluable_litter, ixC14, :, 1) = &
        &             sb_formation_mt(ix_soluable_litter, ixC14, :, 1) &
        &             + fc_soluable_fine_root(:) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixC14, :) + veg_litterfall_mt(ix_heart_wood, ixC14, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_soluable_litter, ixN15, :, 1) = &
        &             sb_formation_mt(ix_soluable_litter, ixN15, :, 1) &
        &             + fn_soluable_fine_root(:) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixN15, :) + veg_litterfall_mt(ix_heart_wood, ixN15, :)) &
        &             / soil_depth_sl(:,1)
      ! polymeric litter
      sb_formation_mt(ix_polymeric_litter, ixC, :, 1) = &
        &             sb_formation_mt(ix_polymeric_litter, ixC, :, 1) &
        &             + (1.0_wp - fc_soluable_fine_root(:)) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixC, :) + veg_litterfall_mt(ix_heart_wood, ixC, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_polymeric_litter, ixN, :, 1) = &
        &             sb_formation_mt(ix_polymeric_litter, ixN, :, 1) &
        &             + (1.0_wp - fn_soluable_fine_root(:)) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixN, :) + veg_litterfall_mt(ix_heart_wood, ixN, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_polymeric_litter, ixP, :, 1) = &
        &             sb_formation_mt(ix_polymeric_litter, ixP, :, 1) &
        &             + (1.0_wp - fp_soluable_fine_root(:)) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixP, :) + veg_litterfall_mt(ix_heart_wood, ixP, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_polymeric_litter, ixC13, :, 1) = &
        &             sb_formation_mt(ix_polymeric_litter, ixC13, :, 1) &
        &             + (1.0_wp - fc_soluable_fine_root(:)) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixC13, :) + veg_litterfall_mt(ix_heart_wood, ixC13, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_polymeric_litter, ixC14, :, 1) = &
        &             sb_formation_mt(ix_polymeric_litter, ixC14, :, 1) &
        &             + (1.0_wp - fc_soluable_fine_root(:)) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixC14, :) + veg_litterfall_mt(ix_heart_wood, ixC14, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_polymeric_litter, ixN15, :, 1) = &
        &             sb_formation_mt(ix_polymeric_litter, ixN15, :, 1) &
        &             + (1.0_wp - fn_soluable_fine_root(:)) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixN15, :) + veg_litterfall_mt(ix_heart_wood, ixN15, :)) &
        &             / soil_depth_sl(:,1)
    END SELECT

    !>4.0 Distribute fine and coarse root litter to soil depth according to the root profile
    !>
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        ! soluable litter (fine root & coarse root)
        sb_formation_mt(ix_soluable_litter, ixC, ic, is) = &
          &             sb_formation_mt(ix_soluable_litter, ixC, ic, is) &
          &             + (fc_soluable_fine_root(ic)   * veg_litterfall_mt(ix_fine_root, ixC, ic) &
          &             +  fc_soluable_coarse_root(ic) * veg_litterfall_mt(ix_coarse_root, ixC, ic)) &
          &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        sb_formation_mt(ix_soluable_litter, ixN, ic, is) = &
          &             sb_formation_mt(ix_soluable_litter, ixN, ic, is) &
          &             + (fn_soluable_fine_root(ic)   * veg_litterfall_mt(ix_fine_root, ixN, ic) &
          &             +  fn_soluable_coarse_root(ic) * veg_litterfall_mt(ix_coarse_root, ixN, ic)) &
          &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        sb_formation_mt(ix_soluable_litter, ixP, ic, is) = &
          &             sb_formation_mt(ix_soluable_litter, ixP, ic, is) &
          &             + (fp_soluable_fine_root(ic)   * veg_litterfall_mt(ix_fine_root, ixP, ic) &
          &             +  fp_soluable_coarse_root(ic) * veg_litterfall_mt(ix_coarse_root, ixP, ic)) &
          &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        sb_formation_mt(ix_soluable_litter, ixC13, ic, is) = &
          &             sb_formation_mt(ix_soluable_litter, ixC13, ic, is) &
          &             + (fc_soluable_fine_root(ic)   * veg_litterfall_mt(ix_fine_root, ixC13, ic) &
          &             +  fc_soluable_coarse_root(ic) * veg_litterfall_mt(ix_coarse_root, ixC13, ic)) &
          &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        sb_formation_mt(ix_soluable_litter, ixC14, ic, is) = &
          &             sb_formation_mt(ix_soluable_litter, ixC14, ic, is) &
          &             + (fc_soluable_fine_root(ic)   * veg_litterfall_mt(ix_fine_root, ixC14, ic) &
          &             +  fc_soluable_coarse_root(ic) * veg_litterfall_mt(ix_coarse_root, ixC14, ic)) &
          &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        sb_formation_mt(ix_soluable_litter, ixN15, ic, is) = &
          &             sb_formation_mt(ix_soluable_litter, ixN15, ic, is) &
          &             + (fn_soluable_fine_root(ic)   * veg_litterfall_mt(ix_fine_root, ixN15, ic) &
          &             +  fn_soluable_coarse_root(ic) * veg_litterfall_mt(ix_coarse_root, ixN15, ic)) &
          &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        ! simple soil model / jsm
        SELECT CASE(TRIM(sb_model_scheme))
        CASE ("simple_1d")
          ! polymeric litter (fine root & coarse root)
          sb_formation_mt(ix_polymeric_litter, ixC, ic, is)    = &
            &           sb_formation_mt(ix_polymeric_litter, ixC, ic, is) &
            &           + ((1.0_wp - fc_soluable_fine_root(ic))   * veg_litterfall_mt(ix_fine_root, ixC, ic) &
            &           +  (1.0_wp - fc_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixC, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixN, ic, is)    = &
            &           sb_formation_mt(ix_polymeric_litter, ixN, ic, is) &
            &           + ((1.0_wp - fn_soluable_fine_root(ic))   * veg_litterfall_mt(ix_fine_root, ixN, ic) &
            &           +  (1.0_wp - fn_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixN, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixP, ic, is)    = &
            &           sb_formation_mt(ix_polymeric_litter, ixP, ic, is) &
            &           + ((1.0_wp - fp_soluable_fine_root(ic))   * veg_litterfall_mt(ix_fine_root, ixP, ic) &
            &           +  (1.0_wp - fp_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixP, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixC13, ic, is)    = &
            &           sb_formation_mt(ix_polymeric_litter, ixC13, ic, is) &
            &           + ((1.0_wp - fc_soluable_fine_root(ic))   * veg_litterfall_mt(ix_fine_root, ixC13, ic) &
            &           +  (1.0_wp - fc_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixC13, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixC14, ic, is)    = &
            &           sb_formation_mt(ix_polymeric_litter, ixC14, ic, is) &
            &           + ((1.0_wp - fc_soluable_fine_root(ic))   * veg_litterfall_mt(ix_fine_root, ixC14, ic) &
            &           +  (1.0_wp - fc_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixC14, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixN15, ic, is)    = &
            &           sb_formation_mt(ix_polymeric_litter, ixN15, ic, is) &
            &           + ((1.0_wp - fn_soluable_fine_root(ic))   * veg_litterfall_mt(ix_fine_root, ixN15, ic) &
            &           +  (1.0_wp - fn_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixN15, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
        CASE ("jsm")
          ! polymeric litter (fine root)
          sb_formation_mt(ix_polymeric_litter, ixC, ic, is) = &
            &           sb_formation_mt(ix_polymeric_litter, ixC, ic, is) &
            &           + ((1.0_wp - fc_soluable_fine_root(ic)) * veg_litterfall_mt(ix_fine_root, ixC, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixN, ic, is) = &
            &           sb_formation_mt(ix_polymeric_litter, ixN, ic, is) &
            &           + ((1.0_wp - fn_soluable_fine_root(ic)) * veg_litterfall_mt(ix_fine_root, ixN, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixP, ic, is) = &
            &           sb_formation_mt(ix_polymeric_litter, ixP, ic, is) &
            &           + ((1.0_wp - fp_soluable_fine_root(ic)) * veg_litterfall_mt(ix_fine_root, ixP, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixC13, ic, is) = &
            &           sb_formation_mt(ix_polymeric_litter, ixC13, ic, is) &
            &           + ((1.0_wp - fc_soluable_fine_root(ic)) * veg_litterfall_mt(ix_fine_root, ixC13, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixC14, ic, is) = &
            &           sb_formation_mt(ix_polymeric_litter, ixC14, ic, is) &
            &           + ((1.0_wp - fc_soluable_fine_root(ic)) * veg_litterfall_mt(ix_fine_root, ixC14, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixN15, ic, is) = &
            &           sb_formation_mt(ix_polymeric_litter, ixN15, ic, is) &
            &           + ((1.0_wp - fn_soluable_fine_root(ic)) * veg_litterfall_mt(ix_fine_root, ixN15, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          ! polymeric litter (coarse root)
          sb_formation_mt(ix_woody_litter, ixC, ic, is) = &
            &           sb_formation_mt(ix_woody_litter, ixC, ic, is) &
            &           + ((1.0_wp - fc_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixC, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_woody_litter, ixN, ic, is) = &
            &           sb_formation_mt(ix_woody_litter, ixN, ic, is) &
            &           + ((1.0_wp - fn_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixN, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_woody_litter, ixP, ic, is) = &
            &           sb_formation_mt(ix_woody_litter, ixP, ic, is) &
            &           + ((1.0_wp - fp_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixP, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_woody_litter, ixC13, ic, is) = &
            &           sb_formation_mt(ix_woody_litter, ixC13, ic, is) &
            &           + ((1.0_wp - fc_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixC13, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_woody_litter, ixC14, ic, is) = &
            &           sb_formation_mt(ix_woody_litter, ixC14, ic, is) &
            &           + ((1.0_wp - fc_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixC14, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_woody_litter, ixN15, ic, is) = &
            &           sb_formation_mt(ix_woody_litter, ixN15, ic, is) &
            &           + ((1.0_wp - fn_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixN15, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
        END SELECT
      END DO
    END DO
  END SUBROUTINE calc_litter_partitioning

#endif
END MODULE mo_q_sb_litter_processes
