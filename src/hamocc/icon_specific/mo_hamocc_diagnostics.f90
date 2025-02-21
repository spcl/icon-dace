! calculates HAMOCC diagnostics:
! monitoring variables, global inventories!
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

MODULE mo_hamocc_diagnostics
!

   USE mo_kind,     ONLY: wp
   USE mo_sync,     ONLY: global_sum_array
   USE mo_sedmnt,   ONLY: ks, porsol, porwat, seddw
   USE mo_exception, ONLY: message_to_own_unit
   USE mo_impl_constants, ONLY: max_char_length
   USE mo_hamocc_types, ONLY: t_hamocc_state
   USE mo_ocean_types, ONLY: t_hydro_ocean_state
   USE mo_model_domain, ONLY: t_patch_3d, t_patch
   USE mo_grid_subset, ONLY: t_subset_range, get_index_range
   USE mo_hamocc_nml, ONLY: io_stdo_bgc, l_cyadyn, l_N_cycle
   USE mo_dynamics_config,     ONLY: nold
   USE mo_ocean_nml,   ONLY: n_zlev,no_tracer
   USE mo_bgc_constants, ONLY:  n2tgn, c2gtc, kilo
   USE mo_memory_bgc, ONLY: p2gtc, totalarea
   USE mo_control_bgc, ONLY: dtbgc, bgc_nproma, bgc_zlevs
   USE mo_bgc_icon_comm, ONLY: to_bgcout
   USE mo_param1_bgc, ONLY: isco212, ialkali, iphosph,iano3, igasnit, &
&                           iphy, izoo, icya, ioxygen, isilica, idoc, &
&                           ian2o, idet, iiron, icalc, iopal,         &
&                           iammo, iano2

   USE mo_name_list_output_init, ONLY: isRegistered
   USE mo_statistics, ONLY: levels_horizontal_mean
   USE mo_fortran_tools, ONLY: set_acc_host_or_device





IMPLICIT NONE

PRIVATE

PUBLIC:: get_inventories, get_monitoring,get_omz


CONTAINS

SUBROUTINE get_omz(hamocc_state, p_patch_3d, pddpo, ssh, lacc)

TYPE(t_hamocc_state) :: hamocc_state
TYPE(t_patch_3d ),TARGET, INTENT(in)   :: p_patch_3d
REAL(wp), INTENT(IN) :: pddpo(:,:,:)
REAL(wp), INTENT(IN) :: ssh(:,:)
LOGICAL, INTENT(IN), OPTIONAL :: lacc

! Local variables
INTEGER :: jc, jb, jk
INTEGER :: start_index, end_index
TYPE(t_subset_range), POINTER :: all_cells
INTEGER:: i_time_stat
INTEGER:: max_lev, omz_depth_index
REAL(wp):: ref_o2
LOGICAL :: lzacc

CALL set_acc_host_or_device(lzacc, lacc)

all_cells => p_patch_3d%p_patch_2d(1)%cells%ALL
i_time_stat=nold(1)


DO jb = all_cells%start_block, all_cells%end_block

    CALL get_index_range(all_cells, jb, start_index, end_index)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO jc=start_index, end_index
     
      max_lev = p_patch_3D%p_patch_1d(1)%dolic_c(jc,jb) 

      IF (max_lev > 0)then

        !omz_depth_index = calc_omz_depth_index(max_levels,o2)
        omz_depth_index=max_lev
        ref_o2=100._wp
        !$ACC LOOP SEQ
        DO jk = 1, max_lev
          IF (hamocc_state%p_prog(i_time_stat)%tracer(jc,jk,jb,ioxygen) < ref_o2)then
            omz_depth_index=jk
            ref_o2 = hamocc_state%p_prog(i_time_stat)%tracer(jc,jk,jb,ioxygen)
          END IF
        END DO

        !o2min = o2(jc,omz_depth_index,jb) in mol m-3
        hamocc_state%p_tend%o2min(jc,jb)=kilo*hamocc_state%p_prog(i_time_stat)%tracer(jc,omz_depth_index,jb,ioxygen)
    
        !zo2min = sum(thickness(1:omz_depth_index)) + zeta
        hamocc_state%p_tend%zo2min(jc,jb)=SUM(pddpo(jc,1:omz_depth_index,jb)) + ssh(jc,jb)

       END IF

    END DO
    !$ACC END PARALLEL

END DO

END SUBROUTINE get_omz

  SUBROUTINE get_monitoring(hamocc_state, tracer, ssh, pddpo, p_patch_3d, lacc)

    USE mo_memory_bgc, ONLY: n2prod, doccya_fac, rcar
    TYPE(t_hamocc_state) :: hamocc_state
    REAL(wp), INTENT(IN) :: ssh(:,:)
    REAL(wp), INTENT(IN) :: pddpo(:,:,:)
    REAL(wp), INTENT(IN) :: tracer(:,:,:,:)
    TYPE(t_patch_3d),TARGET, INTENT(in)   :: p_patch_3d
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    REAL(wp) :: glob_n2b, glob_pwn2b, glob_srf_thick
    REAL(wp) :: glob_det, glob_doc, glob_phy, glob_zoo, glob_cya
    REAL(wp) :: glob_dic, glob_calc, glob_pwic, glob_sedo12
    REAL(wp) :: glob_bo12, glob_sedc12, glob_bc12 

    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_subset_range), POINTER :: owned_cells
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    patch_2d => p_patch_3d%p_patch_2d(1)
    owned_cells    => patch_2d%cells%owned

    glob_n2b = 0.0_wp
    glob_pwn2b = 0.0_wp
    glob_det = 0.0_wp
    glob_doc = 0.0_wp
    glob_phy = 0.0_wp
    glob_zoo = 0.0_wp
    glob_cya = 0.0_wp
    glob_dic = 0.0_wp
    glob_calc = 0.0_wp
    glob_pwic = 0.0_wp
    glob_sedo12 = 0.0_wp
    glob_bo12 = 0.0_wp
    glob_sedc12 = 0.0_wp
    glob_bc12 = 0.0_wp

    IF (isRegistered('global_carbon_inventory')) THEN
      CALL calc_inventory3d(p_patch_3d, &
         &                  ssh, &
         &                  pddpo, &
         &                  tracer(:,:,:,idoc), &
         &                  glob_doc)
      CALL calc_inventory3d(p_patch_3d, &
         &                  ssh, &
         &                  pddpo, &
         &                  tracer(:,:,:,idet), &
         &                  glob_det)
      CALL calc_inventory3d(p_patch_3d, &
         &                  ssh, &
         &                  pddpo, &
         &                  tracer(:,:,:,iphy), &
         &                  glob_phy)
      CALL calc_inventory3d(p_patch_3d, &
         &                  ssh, &
         &                  pddpo, &
         &                  tracer(:,:,:,izoo), &
         &                  glob_zoo)
      CALL calc_inventory3d(p_patch_3d, &
         &                  ssh, &
         &                  pddpo, &
         &                  tracer(:,:,:,icya), &
         &                  glob_cya)
      CALL calc_inventory3d(p_patch_3d, &
         &                  ssh, &
         &                  pddpo, &
         &                  tracer(:,:,:,isco212), &
         &                  glob_dic)
      CALL calc_inventory3d(p_patch_3d, &
         &                  ssh, &
         &                  pddpo, &
         &                  tracer(:,:,:,icalc), &
         &                  glob_calc)
      CALL calc_inventory_sed(p_patch_3d, &
         &                    hamocc_state%p_sed%pwic(:,:,:), &
         &                    porwat, &
         &                    glob_pwic)
      CALL calc_inventory_sed(p_patch_3d, &
         &                    hamocc_state%p_sed%so12(:,:,:), &
         &                    porsol, &
         &                    glob_sedo12)
      CALL calc_inventory_sed(p_patch_3d, &
         &                    hamocc_state%p_sed%sc12(:,:,:), &
         &                    porsol, &
         &                    glob_sedc12)
      CALL calc_inventory2d(p_patch_3d, &
         &                  hamocc_state%p_sed%bo12(:,:), &
         &                  glob_bo12,-2)
      CALL calc_inventory2d(p_patch_3d, &
         &                  hamocc_state%p_sed%bc12(:,:), &
         &                  glob_bc12,-2)
    ENDIF

    IF (isRegistered('global_primary_production')) THEN
      CALL calc_inventory3d(p_patch_3d, &
        &                   ssh, &
        &                   pddpo, &
        &                   hamocc_state%p_tend%npp(:,:,:), &
        &                   hamocc_state%p_tend%monitor%phosy(1))
    ENDIF
    IF (isRegistered('global_npp_cya')) THEN
      CALL calc_inventory3d(p_patch_3d, &
        &                   ssh, &
        &                   pddpo, &
       &                    hamocc_state%p_tend%phoc(:,:,:), &
       &                    hamocc_state%p_tend%monitor%phosy_cya(1))
    ENDIF
    IF (isRegistered('global_zooplankton_grazing')) THEN
      CALL calc_inventory3d(p_patch_3d, &
        &                   ssh, &
        &                   pddpo, &
        &                   hamocc_state%p_tend%graz(:,:,:), &
        &                   hamocc_state%p_tend%monitor%grazing(1))
    ENDIF
    IF (isRegistered('global_remin_via_grazer')) THEN
      CALL calc_inventory3d(p_patch_3d, &
        &                   ssh, &
        &                   pddpo, &
        &                   hamocc_state%p_tend%graton(:,:,:), &
        &                   hamocc_state%p_tend%monitor%graton(1))
    ENDIF
    IF (isRegistered('global_exudation_phytoplankton')) THEN
      CALL calc_inventory3d(p_patch_3d, &
        &                   ssh, &
        &                   pddpo, &
        &                   hamocc_state%p_tend%exud(:,:,:), &
        &                   hamocc_state%p_tend%monitor%exud(1))
    ENDIF
    IF (isRegistered('global_exudation_zooplankton')) THEN
      CALL calc_inventory3d(p_patch_3d, &
        &                   ssh, &
        &                   pddpo, &
        &                   hamocc_state%p_tend%exudz(:,:,:), &
        &                   hamocc_state%p_tend%monitor%exudz(1))
    ENDIF
    IF (isRegistered('global_zooplankton_dying')) THEN
      CALL calc_inventory3d(p_patch_3d, &
        &                   ssh, &
        &                   pddpo, &
        &                   hamocc_state%p_tend%zoomor(:,:,:), &
        &                   hamocc_state%p_tend%monitor%zoomor(1))
    ENDIF
    IF (isRegistered('global_phytoplankton_dying')) THEN
      CALL calc_inventory3d(p_patch_3d, &
        &                   ssh, &
        &                   pddpo, &
        &                   hamocc_state%p_tend%phymor(:,:,:), &
        &                   hamocc_state%p_tend%monitor%phymor(1))
    ENDIF
    IF (isRegistered('global_opal_production')) THEN
      CALL calc_inventory3d(p_patch_3d, &
        &                   ssh, &
        &                   pddpo, &
        &                   hamocc_state%p_tend%delsil(:,:,:), &
        &                   hamocc_state%p_tend%monitor%delsil(1))
    ENDIF
    IF (isRegistered('global_caco3_production')) THEN
      CALL calc_inventory3d(p_patch_3d, &
        &                   ssh, &
        &                   pddpo, &
        &                   hamocc_state%p_tend%delcar(:,:,:), &
        &                   hamocc_state%p_tend%monitor%delcar(1))
    ENDIF
    IF (isRegistered('bacterial_activity')) THEN
      CALL calc_inventory3d(p_patch_3d, &
        &                   ssh, &
        &                   pddpo, &
        &                   hamocc_state%p_tend%bacfra(:,:,:), &
        &                   hamocc_state%p_tend%monitor%bacfra(1))
    ENDIF
    IF (isRegistered('Aerob_remin_of_detritus')) THEN
      CALL calc_inventory3d(p_patch_3d, &
        &                   ssh, &
        &                   pddpo, &
        &                   hamocc_state%p_tend%remina(:,:,:), &
        &                   hamocc_state%p_tend%monitor%remina(1))
    ENDIF
    IF (isRegistered('remin_of_det_by_S')) THEN
      CALL calc_inventory3d(p_patch_3d, &
        &                   ssh, &
        &                   pddpo, &
        &                   hamocc_state%p_tend%remins(:,:,:), &
        &                   hamocc_state%p_tend%monitor%remins(1))
    ENDIF
    IF (l_cyadyn) THEN
      IF (isRegistered('N2_fixation')) THEN
      CALL calc_inventory3d(p_patch_3d, &
        &                   ssh, &
        &                   pddpo, &
        &                   hamocc_state%p_tend%nfix(:,:,:), &
        &                   hamocc_state%p_tend%monitor%n2fix(1))
      ENDIF
      IF (isRegistered('global_cya_loss_det')) THEN
        CALL calc_inventory3d(p_patch_3d, &
        &                     ssh, &
        &                     pddpo, &
        &                     hamocc_state%p_tend%cyloss(:,:,:), &
        &                     hamocc_state%p_tend%monitor%cyaldet(1))
      ENDIF
    ELSE
      IF (isRegistered('N2_fixation')) THEN
        CALL calc_inventory2d(p_patch_3d, &
          &                   hamocc_state%p_tend%nfixd(:,:), &
          &                   hamocc_state%p_tend%monitor%n2fix(1), &
          &                   1, &
          &                   ssh, &
          &                   pddpo)
      ENDIF
    ENDIF
    IF (isRegistered('WC_denit')) THEN
      CALL calc_inventory3d(p_patch_3d, &
        &                   ssh, &
        &                   pddpo, &
        &                   hamocc_state%p_tend%reminn(:,:,:), &
        &                   hamocc_state%p_tend%monitor%wcdenit(1))
    ENDIF
    IF (isRegistered('global_net_co2_flux')) THEN
      CALL calc_inventory2d(p_patch_3d, &
        &                   hamocc_state%p_tend%cflux(:,:), &
        &                   hamocc_state%p_tend%monitor%net_co2_flux(1), &
        &                   -2)
    ENDIF
    IF (isRegistered('global_OM_export_at_90m')) THEN
      CALL calc_inventory2d(p_patch_3d, &
        &                   hamocc_state%p_tend%coex90(:,:), &
        &                   hamocc_state%p_tend%monitor%omex90(1), &
        &                   -2)
    ENDIF
    IF (isRegistered('global_calc_export_at_90m')) THEN
      CALL calc_inventory2d(p_patch_3d, &
        &                   hamocc_state%p_tend%calex90(:,:), &
        &                   hamocc_state%p_tend%monitor%calex90(1), &
        &                   -2)
    ENDIF
    IF (isRegistered('global_opal_export_at_90m')) THEN
      CALL calc_inventory2d(p_patch_3d, &
        &                   hamocc_state%p_tend%opex90(:,:), &
        &                   hamocc_state%p_tend%monitor%opex90(1), &
        &                   -2)
    ENDIF
    IF (isRegistered('global_OM_export_at_1000m')) THEN
      CALL calc_inventory2d(p_patch_3d, &
        &                   hamocc_state%p_tend%coex1000(:,:), &
        &                   hamocc_state%p_tend%monitor%omex1000(1), &
        &                   -2)
    ENDIF
    IF (isRegistered('global_calc_export_at_1000m')) THEN
      CALL calc_inventory2d(p_patch_3d, &
        &                   hamocc_state%p_tend%calex1000(:,:), &
        &                   hamocc_state%p_tend%monitor%calex1000(1), &
        &                   -2)
    ENDIF
    IF (isRegistered('global_opal_export_at_1000m')) THEN
      CALL calc_inventory2d(p_patch_3d, &
        &                   hamocc_state%p_tend%opex1000(:,:), &
        &                   hamocc_state%p_tend%monitor%opex1000(1), &
        &                   -2)
    ENDIF
    IF (isRegistered('global_OM_export_at_2000m')) THEN
      CALL calc_inventory2d(p_patch_3d, &
        &                   hamocc_state%p_tend%coex2000(:,:), &
        &                   hamocc_state%p_tend%monitor%omex2000(1), &
        &                   -2)
    ENDIF
    IF (isRegistered('global_calc_export_at_2000m')) THEN
      CALL calc_inventory2d(p_patch_3d, &
        &                   hamocc_state%p_tend%calex2000(:,:), &
        &                   hamocc_state%p_tend%monitor%calex2000(1), &
        &                   -2)
    ENDIF
    IF (isRegistered('global_opal_export_at_2000m')) THEN
      CALL calc_inventory2d(p_patch_3d, &
        &                   hamocc_state%p_tend%opex2000(:,:), &
        &                   hamocc_state%p_tend%monitor%opex2000(1), &
        &                   -2)
    ENDIF
    IF (isRegistered('global_surface_alk')) THEN
      CALL calc_inventory2d(p_patch_3d, &
        &                   tracer(:,1,:,ialkali), &
        &                   hamocc_state%p_tend%monitor%sfalk(1), &
        &                   1, &
        &                   ssh, &
        &                   pddpo)
    ENDIF
    IF (isRegistered('global_surface_dic')) THEN
      CALL calc_inventory2d(p_patch_3d, &
        &                   tracer(:,1,:,isco212), &
        &                   hamocc_state%p_tend%monitor%sfdic(1), &
        &                   1, &
        &                   ssh, &
        &                   pddpo)
    ENDIF
    IF (isRegistered('global_surface_phosphate')) THEN
      CALL calc_inventory2d(p_patch_3d, &
        &                   tracer(:,1,:,iphosph), &
        &                   hamocc_state%p_tend%monitor%sfphos(1), &
        &                   1, &
        &                   ssh, &
        &                   pddpo)
    ENDIF
    IF (isRegistered('global_surface_silicate')) THEN
      CALL calc_inventory2d(p_patch_3d, &
        &                   tracer(:,1,:,isilica), &
        &                   hamocc_state%p_tend%monitor%sfsil(1), &
        &                   1, &
        &                   ssh, &
        &                   pddpo)
    ENDIF
    IF (isRegistered('global_surface_nitrate')) THEN
      CALL calc_inventory2d(p_patch_3d, &
        &                   tracer(:,1,:,iano3), &
        &                   hamocc_state%p_tend%monitor%sfnit(1), &
         &                   1, &
        &                   ssh, &
        &                   pddpo)
    ENDIF

    IF (isRegistered('global_zalkn2')) THEN
      CALL calc_inventory_sed(p_patch_3d, &
        &                     hamocc_state%p_sed%pwn2b(:,:,:), &
        &                     porwat, &
        &                     glob_pwn2b)
      CALL calc_inventory3d(p_patch_3d, &
        &                   ssh, &
        &                   pddpo, &
        &                   hamocc_state%p_tend%n2budget(:,:,:), &
        &                   glob_n2b, &
        &                   .TRUE.)
    ENDIF

    IF (isRegistered('SED_denit')) THEN
      CALL calc_inventory_sed(p_patch_3d, &
        &                     hamocc_state%p_tend%sedrn(:,:,:), &
        &                     porwat, &
        &                     hamocc_state%p_tend%monitor%seddenit(1))
    ENDIF


    IF (l_N_cycle) THEN
      IF (isRegistered('global_primary_production_nh4')) THEN
        CALL calc_inventory3d(p_patch_3d, &
          &                   ssh, &
          &                   pddpo, &
          &                   hamocc_state%p_tend%gppnh4(:,:,:), &
          &                   hamocc_state%p_tend%monitor%phosy_nh4(1))
      ENDIF

      IF (isRegistered('global_npp_cya_nh4')) THEN
        CALL calc_inventory3d(p_patch_3d, &
          &                   ssh, &
          &                   pddpo, &
          &                   hamocc_state%p_tend%cyapro(:,:,:), &
          &                   hamocc_state%p_tend%monitor%phosy_cya_nh4(1))
      ENDIF

      IF (isRegistered('global_net_nh3_flux')) THEN
        CALL calc_inventory2d(p_patch_3d, &
          &                   hamocc_state%p_tend%nh3flux(:,:), &
          &                   hamocc_state%p_tend%monitor%net_nh3_flux(1), &
          &                   -2)
      ENDIF

      IF (isRegistered('global_surface_nh4')) THEN
        CALL calc_inventory2d(p_patch_3d, &
          &                   tracer(:,1,:,iammo), &
          &                    hamocc_state%p_tend%monitor%sfnh4(1), &
          &                   1, &
          &                   ssh, &
          &                   pddpo)
      ENDIF

      IF (isRegistered('WC_nitri_no2')) THEN
        CALL calc_inventory3d(p_patch_3d, &
          &                   ssh, &
          &                   pddpo, &
          &                   hamocc_state%p_tend%nitox(:,:,:), &
          &                   hamocc_state%p_tend%monitor%wc_nitri_no2(1))
      ENDIF

      IF (isRegistered('WC_nitri_nh4')) THEN
        CALL calc_inventory3d(p_patch_3d, &
          &                   ssh, &
          &                   pddpo, &
          &                   hamocc_state%p_tend%ammox(:,:,:), &
          &                   hamocc_state%p_tend%monitor%wc_nitri_nh4(1))
      ENDIF

      IF (isRegistered('WC_dnrn')) THEN
        CALL calc_inventory3d(p_patch_3d, &
          &                   ssh, &
          &                   pddpo, &
          &                   hamocc_state%p_tend%dnrn(:,:,:), &
          &                   hamocc_state%p_tend%monitor%wc_dnrn(1))
      ENDIF

      IF (isRegistered('WC_dnra')) THEN
        CALL calc_inventory3d(p_patch_3d, &
          &                   ssh, &
          &                   pddpo, &
          &                   hamocc_state%p_tend%dnra(:,:,:), &
          &                   hamocc_state%p_tend%monitor%wc_dnra(1))
      ENDIF

      IF (isRegistered('WC_anammox')) THEN
        CALL calc_inventory3d(p_patch_3d, &
          &                   ssh, &
          &                   pddpo, &
          &                   hamocc_state%p_tend%anam(:,:,:), &
          &                   hamocc_state%p_tend%monitor%wc_anammox(1))
      ENDIF
    ENDIF

    ! Unit conversion
    hamocc_state%p_tend%monitor%phosy(1)        = hamocc_state%p_tend%monitor%phosy(1) * p2gtc
    hamocc_state%p_tend%monitor%phosy_cya(1)    = hamocc_state%p_tend%monitor%phosy_cya(1) * p2gtc
    hamocc_state%p_tend%monitor%grazing(1)      = hamocc_state%p_tend%monitor%grazing(1) * p2gtc
    hamocc_state%p_tend%monitor%remina(1)       = hamocc_state%p_tend%monitor%remina(1) * p2gtc
    hamocc_state%p_tend%monitor%remins(1)       = hamocc_state%p_tend%monitor%remins(1) * p2gtc
    hamocc_state%p_tend%monitor%exud(1)         = hamocc_state%p_tend%monitor%exud(1) * p2gtc
    hamocc_state%p_tend%monitor%exudz(1)        = hamocc_state%p_tend%monitor%exudz(1) * p2gtc
    hamocc_state%p_tend%monitor%zoomor(1)       = hamocc_state%p_tend%monitor%zoomor(1) * p2gtc
    hamocc_state%p_tend%monitor%phymor(1)       = hamocc_state%p_tend%monitor%phymor(1) * p2gtc
    hamocc_state%p_tend%monitor%graton(1)       = hamocc_state%p_tend%monitor%graton(1) * p2gtc
    hamocc_state%p_tend%monitor%bacfra(1)       = hamocc_state%p_tend%monitor%bacfra(1) * p2gtc
    hamocc_state%p_tend%monitor%net_co2_flux(1) = hamocc_state%p_tend%monitor%net_co2_flux(1) * c2gtc
    hamocc_state%p_tend%monitor%delcar(1)       = hamocc_state%p_tend%monitor%delcar(1) * c2gtc
    
    hamocc_state%p_tend%monitor%n2fix(1)        = hamocc_state%p_tend%monitor%n2fix(1) * n2tgn
    hamocc_state%p_tend%monitor%omex90(1)       = hamocc_state%p_tend%monitor%omex90(1) * p2gtc
    hamocc_state%p_tend%monitor%calex90(1)      = hamocc_state%p_tend%monitor%calex90(1) * c2gtc
    hamocc_state%p_tend%monitor%omex1000(1)     = hamocc_state%p_tend%monitor%omex1000(1) * p2gtc
    hamocc_state%p_tend%monitor%calex1000(1)    = hamocc_state%p_tend%monitor%calex1000(1) * c2gtc
    hamocc_state%p_tend%monitor%omex2000(1)     = hamocc_state%p_tend%monitor%omex2000(1) * p2gtc
    hamocc_state%p_tend%monitor%calex2000(1)    = hamocc_state%p_tend%monitor%calex2000(1) * c2gtc
    
    hamocc_state%p_tend%monitor%cyaldoc(1)      = hamocc_state%p_tend%monitor%cyaldet(1) * p2gtc * doccya_fac
    hamocc_state%p_tend%monitor%cyaldet(1)      = hamocc_state%p_tend%monitor%cyaldet(1) * p2gtc *(1._wp - doccya_fac)

    CALL levels_horizontal_mean(pddpo(:,1,:), patch_2d%cells%area(:,:), owned_cells, glob_srf_thick, lopenacc=lzacc)

    ! mean values of surface concentrations
    hamocc_state%p_tend%monitor%sfalk(1)        = hamocc_state%p_tend%monitor%sfalk(1)/totalarea/glob_srf_thick
    hamocc_state%p_tend%monitor%sfdic(1)        = hamocc_state%p_tend%monitor%sfdic(1)/totalarea/glob_srf_thick
    hamocc_state%p_tend%monitor%sfsil(1)        = hamocc_state%p_tend%monitor%sfsil(1)/totalarea/glob_srf_thick
    hamocc_state%p_tend%monitor%sfphos(1)       = hamocc_state%p_tend%monitor%sfphos(1)/totalarea/glob_srf_thick
    hamocc_state%p_tend%monitor%sfnit(1)        = hamocc_state%p_tend%monitor%sfnit(1)/totalarea/glob_srf_thick
    hamocc_state%p_tend%monitor%zalkn2(1)       = glob_n2b + glob_pwn2b

    ! total carbon
    hamocc_state%p_tend%monitor%carbinv(1)   = ((glob_det + glob_doc + glob_phy + glob_zoo + glob_cya) * rcar &
           &                                     + glob_dic + glob_calc + glob_pwic &
           &                                     + (glob_sedo12 + glob_bo12) * rcar + glob_sedc12 + glob_bc12 ) &
           &                                     * c2gtc

    IF (.not. l_N_cycle) THEN
      hamocc_state%p_tend%monitor%wcdenit(1)      = hamocc_state%p_tend%monitor%wcdenit(1) * 2._wp * n2prod* n2tgn
      hamocc_state%p_tend%monitor%seddenit(1)     = hamocc_state%p_tend%monitor%seddenit(1) * 2._wp* n2prod* n2tgn
    ELSE
      hamocc_state%p_tend%monitor%wcdenit(1)      = hamocc_state%p_tend%monitor%wcdenit(1) * n2tgn
      hamocc_state%p_tend%monitor%seddenit(1)     = hamocc_state%p_tend%monitor%seddenit(1) * n2tgn

      hamocc_state%p_tend%monitor%phosy_nh4(1)        = hamocc_state%p_tend%monitor%phosy_nh4(1) * p2gtc
      hamocc_state%p_tend%monitor%phosy_cya_nh4(1)    = hamocc_state%p_tend%monitor%phosy_cya_nh4(1) * p2gtc
      hamocc_state%p_tend%monitor%sfnh4(1)            = hamocc_state%p_tend%monitor%sfnh4(1)/totalarea/glob_srf_thick
      hamocc_state%p_tend%monitor%net_nh3_flux(1)     = hamocc_state%p_tend%monitor%net_nh3_flux(1) * n2tgn

      ! LR: should be in N units, but check
      hamocc_state%p_tend%monitor%wc_nitri_no2(1)      = hamocc_state%p_tend%monitor%wc_nitri_no2(1) * n2tgn
      hamocc_state%p_tend%monitor%wc_nitri_nh4(1)      = hamocc_state%p_tend%monitor%wc_nitri_nh4(1) * n2tgn
      hamocc_state%p_tend%monitor%wc_dnrn(1)           = hamocc_state%p_tend%monitor%wc_dnrn(1) * n2tgn
      hamocc_state%p_tend%monitor%wc_dnra(1)           = hamocc_state%p_tend%monitor%wc_dnra(1) * n2tgn
      hamocc_state%p_tend%monitor%wc_anammox(1)        = hamocc_state%p_tend%monitor%wc_anammox(1) * n2tgn
    ENDIF


  END SUBROUTINE get_monitoring

SUBROUTINE get_inventories(hamocc_state,ssh,pddpo, tracer, p_patch_3d, weathering_flag, flux_flag)

USE mo_memory_bgc,      ONLY: rnit,rn2, ro2bal,rcar,ralk

REAL(wp),INTENT(IN) :: ssh(:,:)
REAL(wp),INTENT(IN) :: pddpo(:,:,:)
REAL(wp),INTENT(IN) :: tracer(:,:,:,:)
REAL(wp),INTENT(IN) :: weathering_flag, flux_flag
TYPE(t_hamocc_state) :: hamocc_state
TYPE(t_patch_3d ),TARGET, INTENT(in)   :: p_patch_3d

! Local variables
REAL(wp) :: glob_det, glob_doc, glob_phy, glob_zoo
REAL(wp) :: glob_phos, glob_sedo12, glob_bo12
REAL(wp) :: glob_sedc12, glob_nit, glob_gnit, glob_pwn2b, glob_pwh2ob
REAL(wp) :: glob_n2o,glob_n2fl,glob_n2ofl, glob_orginp,glob_nitinp
REAL(wp) :: glob_calinp, glob_silinp, glob_alk, glob_calc
REAL(wp) :: glob_sil, glob_opal, glob_sedsi, glob_pwsi
REAL(wp) :: glob_bsil, glob_silpro, glob_n2b, glob_h2ob
REAL(wp) :: glob_prorca,  glob_cya, glob_produs
REAL(wp) :: glob_pwn2, glob_pwno3, glob_prcaca, glob_ofl
REAL(wp) :: glob_pwic, glob_pwal, glob_pwph, glob_cfl
REAL(wp) :: glob_dic, glob_o2, glob_fe, glob_co3, glob_hi
REAL(wp) :: glob_pwox, glob_pwfe, glob_bclay, glob_sedclay
REAL(wp) :: total_ocean, watersum, sedsum
REAL(wp) :: rcyano, glob_bc12
! extended N-cycle
REAL(wp) :: glob_nh4, glob_no2, glob_pwnh4, glob_pwno2
REAL(wp) :: glob_nh3fl

CHARACTER(LEN=max_char_length) :: cpara_name, cpara_val

cpara_name='======================='
cpara_val="==========="
CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )

!Check if cyanobacteria are calculcated
rcyano = MERGE(1._wp,0._wp, l_cyadyn)

! Calculate global inventories of individual tracers
! Water column

CALL calc_inventory3d(p_patch_3d, ssh, pddpo, tracer(:,:,:,idoc), &
&                      glob_doc)
CALL calc_inventory3d(p_patch_3d, ssh, pddpo, tracer(:,:,:,idet), &
&                      glob_det)
CALL calc_inventory3d(p_patch_3d, ssh, pddpo, tracer(:,:,:,iphy), &
&                      glob_phy)
CALL calc_inventory3d(p_patch_3d, ssh, pddpo, tracer(:,:,:,izoo), &
&                      glob_zoo)
CALL calc_inventory3d(p_patch_3d, ssh, pddpo, tracer(:,:,:,iphosph),&
&                      glob_phos)
CALL calc_inventory3d(p_patch_3d, ssh, pddpo, tracer(:,:,:,iano3),&
&                      glob_nit)
CALL calc_inventory3d(p_patch_3d, ssh, pddpo, tracer(:,:,:,igasnit),&
&                      glob_gnit)
CALL calc_inventory3d(p_patch_3d, ssh, pddpo, tracer(:,:,:,ian2o), &
&                      glob_n2o)
CALL calc_inventory3d(p_patch_3d, ssh, pddpo, tracer(:,:,:,isilica), &
&                      glob_sil)
CALL calc_inventory3d(p_patch_3d, ssh, pddpo, tracer(:,:,:,icalc), &
&                      glob_calc)
CALL calc_inventory3d(p_patch_3d, ssh, pddpo, tracer(:,:,:,iopal), &
&                      glob_opal)
CALL calc_inventory3d(p_patch_3d, ssh, pddpo, tracer(:,:,:,ialkali),&
&                      glob_alk)
CALL calc_inventory3d(p_patch_3d, ssh, pddpo, tracer(:,:,:,isco212),&
&                      glob_dic)
CALL calc_inventory3d(p_patch_3d, ssh, pddpo, tracer(:,:,:,ioxygen),&
&                      glob_o2)
CALL calc_inventory3d(p_patch_3d, ssh, pddpo, tracer(:,:,:,iiron),&
&                      glob_fe)
IF(l_cyadyn)THEN
 CALL calc_inventory3d(p_patch_3d, ssh, pddpo, tracer(:,:,:,icya),&
&                       glob_cya)
else
 glob_cya=0._wp
ENDIF
CALL calc_inventory3d(p_patch_3d, ssh, pddpo, hamocc_state%p_diag%co3(:,:,:), glob_co3)
CALL calc_inventory3d(p_patch_3d, ssh, pddpo, hamocc_state%p_diag%hi(:,:,:), glob_hi)

! Sediment
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%so12(:,:,:), porsol,  glob_sedo12)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%sc12(:,:,:), porsol,  glob_sedc12)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%ster(:,:,:), porsol,  glob_sedclay)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%ssil(:,:,:), porsol,  glob_sedsi)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwsi(:,:,:), porwat,  glob_pwsi)
!CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwsi(:,:,:), porsol, glob_pwsi)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwal(:,:,:), porwat,  glob_pwal)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwic(:,:,:), porwat,  glob_pwic)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwph(:,:,:), porwat,  glob_pwph)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwox(:,:,:), porwat,  glob_pwox)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwfe(:,:,:), porwat,  glob_pwfe)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwn2b(:,:,:), porwat, glob_pwn2b)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwn2(:,:,:), porwat,  glob_pwn2)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwno3(:,:,:), porwat,  glob_pwno3)
CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwh2ob(:,:,:), porwat, glob_pwh2ob)

CALL calc_inventory2d(p_patch_3d, hamocc_state%p_sed%bo12(:,:), glob_bo12,-2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_sed%bsil(:,:), glob_bsil,-2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_sed%bc12(:,:), glob_bc12,-2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_sed%bter(:,:), glob_bclay,-2)

! Tendencies
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_tend%prorca(:,:), glob_prorca,-2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_tend%prcaca(:,:), glob_prcaca,-2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_tend%silpro(:,:), glob_silpro,-2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_tend%produs(:,:), glob_produs,-2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_tend%cflux(:,:), glob_cfl, -2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_tend%oflux(:,:), glob_ofl, -2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_tend%nflux(:,:), glob_n2fl, -2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_tend%n2oflux(:,:), glob_n2ofl, -2)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_tend%orginp(:,:), glob_orginp, 1, ssh, pddpo)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_tend%silinp(:,:), glob_silinp, 1, ssh, pddpo)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_tend%calinp(:,:), glob_calinp, 1, ssh, pddpo)
CALL calc_inventory2d(p_patch_3d, hamocc_state%p_tend%nitrogeninp(:,:), glob_nitinp, 1, ssh, pddpo)
CALL calc_inventory3d(p_patch_3d, ssh, pddpo, hamocc_state%p_tend%h2obudget(:,:,:), glob_h2ob,.TRUE.)
CALL calc_inventory3d(p_patch_3d, ssh, pddpo, hamocc_state%p_tend%n2budget(:,:,:), glob_n2b,.TRUE.)

! Convert unit for fluxes
![kmol/s] to [kmol]

glob_cfl=glob_cfl*dtbgc
glob_ofl=glob_ofl*dtbgc
glob_n2fl=glob_n2fl*dtbgc
glob_n2ofl=glob_n2ofl*dtbgc

IF (l_N_cycle) THEN
   CALL calc_inventory3d(p_patch_3d, ssh, pddpo, tracer(:,:,:,iammo),glob_nh4)
   CALL calc_inventory3d(p_patch_3d, ssh, pddpo, tracer(:,:,:,iano2),glob_no2)
   CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwnh4(:,:,:), porwat,  glob_pwnh4)
   CALL calc_inventory_sed(p_patch_3d, hamocc_state%p_sed%pwno2(:,:,:), porwat,  glob_pwno2)
   CALL calc_inventory2d(p_patch_3d, hamocc_state%p_tend%nh3flux(:,:), glob_nh3fl, -2)
   glob_nh3fl = glob_nh3fl*dtbgc
ENDIF



! Print tracer output
CALL message_to_own_unit(' ', ' ', io_stdo_bgc)
CALL message_to_own_unit('Global inventory of', 'ocean tracers', io_stdo_bgc)
cpara_name='-----------------------'
cpara_val="-----------"
CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )

CALL to_bgcout('DIC',glob_dic)
CALL to_bgcout('Alkalinity',glob_alk)
CALL to_bgcout('Phosphate',glob_phos)
CALL to_bgcout('Oxygen',glob_o2)
CALL to_bgcout('N2',glob_gnit)
CALL to_bgcout('Nitrate',glob_nit)
CALL to_bgcout('Silicate',glob_sil)
CALL to_bgcout('DOC',glob_doc)
CALL to_bgcout('Phytoplankton',glob_phy)
CALL to_bgcout('Zooplankton',glob_zoo)
CALL to_bgcout('N2O',glob_n2o)
CALL to_bgcout('Iron',glob_fe)
CALL to_bgcout('Detritus',glob_det)
IF(l_cyadyn)THEN
 CALL to_bgcout('Cyanobacteria',glob_cya)
ENDIF
CALL to_bgcout('Calc',glob_calc)
CALL to_bgcout('Opal',glob_opal)
IF (l_N_cycle) THEN
   CALL to_bgcout('Ammonium',glob_nh4)
   CALL to_bgcout('Nitrite',glob_no2)
ENDIF

CALL message_to_own_unit('Global inventory of', 'additional tracers', io_stdo_bgc)
cpara_name='-----------------------'
cpara_val="-----------"
CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
CALL to_bgcout('hi',glob_hi)
CALL to_bgcout('co3',glob_co3)

CALL message_to_own_unit(' ', ' ', io_stdo_bgc)
CALL message_to_own_unit('Global inventory of', 'aqueous sediment tracers', io_stdo_bgc)
cpara_name='-----------------------'
cpara_val="-----------"
CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
CALL to_bgcout('DIC',glob_pwic)
CALL to_bgcout('Alkalinity',glob_pwal)
CALL to_bgcout('Phosphate',glob_pwph)
CALL to_bgcout('Oxygen',glob_pwox)
CALL to_bgcout('N2',glob_pwn2)
CALL to_bgcout('Nitrate',glob_pwno3)
CALL to_bgcout('Silicate',glob_pwsi)
CALL to_bgcout('Iron',glob_pwfe)
IF (l_N_cycle) THEN
   CALL to_bgcout('Ammonium',glob_pwnh4)
   CALL to_bgcout('Nitrite',glob_pwno2)
ENDIF

CALL message_to_own_unit(' ', ' ', io_stdo_bgc)
CALL message_to_own_unit('Global inventory of', 'solid sediment constituents', io_stdo_bgc)
cpara_name='-----------------------'
cpara_val="-----------"
CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
CALL to_bgcout('Solid C org',glob_sedo12)
CALL to_bgcout('Burial C org',glob_bo12)
CALL to_bgcout('Solid CaCO3',glob_sedc12)
CALL to_bgcout('Burial CaCO3',glob_bc12)
CALL to_bgcout('Solid opal',glob_sedsi)
CALL to_bgcout('Burial opal',glob_bsil)
CALL to_bgcout('Solid clay',glob_sedclay)
CALL to_bgcout('Burial clay',glob_bclay)

CALL message_to_own_unit(' ', ' ', io_stdo_bgc)
cpara_name='======================='
cpara_val="==========="
CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )

CALL to_bgcout('Redfield (global)',glob_nit/glob_phos)

CALL message_to_own_unit(' ', ' ', io_stdo_bgc)
CALL message_to_own_unit('Global fluxes into', 'atmosphere [kmol]', io_stdo_bgc)
cpara_name='-----------------------'
cpara_val="-----------"
CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
CALL to_bgcout('CO2 flux',glob_cfl)
CALL to_bgcout('O2 flux',glob_ofl)
CALL to_bgcout('N2 flux',glob_n2fl)
CALL to_bgcout('N2O flux',glob_n2ofl)
IF (l_N_cycle) THEN
   CALL to_bgcout('NH3 flux',glob_nh3fl)
ENDIF

CALL message_to_own_unit(' ', ' ', io_stdo_bgc)
CALL message_to_own_unit('Global fluxes into', 'sediment [kmol]', io_stdo_bgc)
cpara_name='-----------------------'
cpara_val="-----------"
CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )

CALL to_bgcout('prcaca',glob_prcaca)
CALL to_bgcout('prorca',glob_prorca)
CALL to_bgcout('silpro',glob_silpro)
CALL to_bgcout('produs',glob_produs)

CALL to_bgcout('zalkn2',glob_n2b+glob_pwn2b)

CALL message_to_own_unit(' ', ' ', io_stdo_bgc)
CALL message_to_own_unit('Global weathering fluxes', ' [kmol]', io_stdo_bgc)
cpara_name='-----------------------'
cpara_val="-----------"
CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )

CALL to_bgcout('orginp',glob_orginp)
CALL to_bgcout('silinp',glob_silinp)
CALL to_bgcout('calcinp',glob_calinp)

cpara_name='======================='
cpara_val="==========="
CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )
CALL message_to_own_unit(' ', ' ', io_stdo_bgc)


! Calculate total inventory diagnostics
! DO not consider prcaca, prorca, silpro and produs, as the call is after the loop
! and the fluxes are zeros after powach, the output in bgcflux however is not equal zero
!-------- Phosphate
watersum = glob_det + glob_doc + glob_phy + glob_zoo + glob_phos  &
     &     + rcyano*glob_cya

sedsum =  glob_sedo12 + glob_bo12 - weathering_flag *  glob_orginp + glob_pwph

total_ocean = watersum + sedsum


CALL to_bgcout('Global water phosphate [kmol]',watersum)
CALL to_bgcout('Global sediment phosphate [kmol]',sedsum)
CALL to_bgcout('Global total phosphate [kmol]',total_ocean)
CALL message_to_own_unit(' ', ' ', io_stdo_bgc)


!-------- Nitrate
watersum = rnit * (glob_det + glob_doc + glob_phy + glob_zoo  &
     &     + rcyano*glob_cya ) + glob_nit     &
     &     + rn2 * (glob_gnit + glob_n2o + flux_flag * (glob_n2fl + glob_n2ofl)) &
     &     + rn2*glob_pwn2 + glob_pwno3 - weathering_flag * glob_nitinp

sedsum =  rnit* (glob_sedo12 + glob_bo12)   &
     &     - weathering_flag * rnit * glob_orginp

IF (l_N_cycle) THEN
   watersum = watersum + glob_nh4 + glob_no2 + flux_flag * glob_nh3fl
   sedsum = sedsum + glob_pwnh4 + glob_pwno2
ENDIF

total_ocean = watersum + sedsum

CALL to_bgcout('Global water nitrate [kmol]',watersum)
CALL to_bgcout('Global sediment nitrate [kmol]',sedsum)
CALL to_bgcout('Global total nitrate [kmol]',total_ocean)
CALL message_to_own_unit(' ', ' ', io_stdo_bgc)

!-------- Silicate
watersum =  glob_sil + glob_opal + glob_pwsi
    
sedsum = glob_sedsi + glob_bsil  - weathering_flag * glob_silinp

total_ocean = watersum + sedsum

CALL to_bgcout('Global water silicate [kmol]',watersum)
CALL to_bgcout('Global sediment silicate [kmol]',sedsum)
CALL to_bgcout('Global total silicate [kmol]',total_ocean)
CALL message_to_own_unit(' ', ' ', io_stdo_bgc)

! Alkalinity

watersum = glob_alk - ralk* (glob_det + glob_doc + glob_phy + glob_zoo &
  &        + rcyano* glob_cya) - (glob_n2b+glob_pwn2b)          &
  &        + 2._wp * glob_calc + weathering_flag * (ralk * glob_orginp - 2._wp * glob_calinp)

sedsum =  2._wp * (glob_sedc12 + glob_bc12) + glob_pwal - ralk * (glob_sedo12 + glob_bo12)

IF (l_N_cycle) THEN
   watersum = watersum - 2._wp * (glob_nh4 + flux_flag * glob_nh3fl)
   sedsum = sedsum - 2._wp * glob_pwnh4
ENDIF

total_ocean = watersum + sedsum

CALL to_bgcout('Global water alkalinity [kmol]',watersum)
CALL to_bgcout('Global sediment alkalinity [kmol]',sedsum)
CALL to_bgcout('Global total alkalinity [kmol]',total_ocean)
CALL message_to_own_unit(' ', ' ', io_stdo_bgc)

! Oxygen

watersum = (glob_det + glob_doc + glob_phy + glob_zoo +         &
  &         rcyano*glob_cya )*(-ro2bal) + &
  &         glob_o2 + glob_phos*2._wp + glob_dic + glob_calc +  &
  &         glob_nit * 1.5_wp + glob_n2o* 0.5_wp + glob_pwno3* 1.5 + &
  &         glob_pwic + glob_pwox + glob_pwph*2._wp + flux_flag * (glob_ofl + &
  &         glob_n2ofl * 0.5_wp + glob_cfl) + glob_h2ob + glob_pwh2ob +    &
  &         weathering_flag * (glob_orginp*ro2bal - glob_calinp  &
  &         - glob_nitinp*1.5_wp)

sedsum =   (glob_sedo12 + glob_bo12)*(-ro2bal) + glob_sedc12 + glob_bc12

IF (l_N_cycle) THEN
   watersum = watersum - 0.5_wp * (glob_nh4 + flux_flag * glob_nh3fl) + glob_no2
   sedsum = sedsum - 0.5_wp * glob_pwnh4 + glob_pwno2
ENDIF

total_ocean = watersum + sedsum

CALL to_bgcout('Global water oxygen [kmol]',watersum)
CALL to_bgcout('Global sediment oxygen [kmol]',sedsum)
CALL to_bgcout('Global total oxygen [kmol]',total_ocean)
CALL message_to_own_unit(' ', ' ', io_stdo_bgc)


! Carbon

watersum = (glob_det + glob_doc + glob_phy + glob_zoo   &
     &     + rcyano*glob_cya ) *rcar + &
     &     glob_dic + glob_calc + flux_flag * glob_cfl - weathering_flag * (glob_calinp + &
     &     rcar * glob_orginp) + glob_pwic
    
sedsum =  (glob_sedo12 + glob_bo12 ) * rcar + glob_sedc12 + glob_bc12


total_ocean = watersum + sedsum


CALL to_bgcout('Global water carbon [kmol]',watersum)
CALL to_bgcout('Global sediment carbon [kmol]',sedsum)
CALL to_bgcout('Global total carbon [kmol]',total_ocean)
CALL message_to_own_unit(' ', ' ', io_stdo_bgc)


cpara_name='======================='
cpara_val="==========="
CALL message_to_own_unit(TRIM(cpara_name), TRIM(cpara_val), io_stdo_bgc )

END SUBROUTINE
!--------------------------------------------------------------------------------------------
SUBROUTINE calc_inventory3d(patch3D, ssh, pddpo, pfield3d, field_globsum, no_thick)
! Calculate inventory of the whole water column
! inv = sum (tracer * volume)
! if no_thick is given:
! inv = sum(tracer * area) --> thickness needs to be cons. elsewhere
REAL(wp), TARGET:: pfield3d(:,:,:)
TYPE(t_patch_3d), TARGET, INTENT(in) :: patch3D
REAL(wp), INTENT(IN), TARGET:: ssh(:,:)
REAL(wp), INTENT(IN), TARGET:: pddpo(:,:,:)
REAL(wp), INTENT(OUT):: field_globsum
LOGICAL, OPTIONAL:: no_thick

! Local
INTEGER:: jk
REAL(wp) :: ptmp
TYPE(t_patch), POINTER :: patch_2d

!$ACC UPDATE HOST(pfield3d) ASYNC(1) IF(i_am_accel_node .AND. acc_is_present(pfield3d))
!$ACC WAIT(1)

patch_2d => patch3D%p_patch_2d(1)

field_globsum=0._wp

IF(PRESENT(no_thick))then
! if multiplied elsewhere, do not consider cell thickness and surface area here
 DO jk=1,n_zlev
  ptmp = global_sum_array(&
&   patch_2d%cells%area(:,:) * &
&   patch3d%wet_halo_zero_c(:,jk,:) * &
&   pfield3d(:,jk,:) &
)
 field_globsum = field_globsum + ptmp
 ENDDO
ELSE
 DO jk=1,n_zlev
  ptmp = global_sum_array(&
&   patch_2d%cells%area(:,:) * &
&   patch3d%wet_halo_zero_c(:,jk,:) * &
&   MERGE(pddpo(:,jk,:)+ssh , pddpo(:,jk,:),jk==1) * &
&   pfield3d(:,jk,:) &
)
 field_globsum = field_globsum + ptmp
 ENDDO

ENDIF

END SUBROUTINE


SUBROUTINE calc_inventory_sed(patch3D, pfield3d, sed_state, field_globsum)
! calculated sediment inventory
! sed_state: porsol or porwat (for solid or pore water tracer)
! inv = sum (tracer * {porsol or porwat} * area * seddw)
REAL(wp), TARGET:: pfield3d(:,:,:)
REAL(wp) :: sed_state(ks)
TYPE(t_patch_3d), TARGET, INTENT(in) :: patch3D
REAL(wp), INTENT(OUT) :: field_globsum

! Local
INTEGER:: jk
TYPE(t_patch), POINTER :: patch_2d
REAL(wp) :: ptmp

!$ACC UPDATE HOST(pfield3d) ASYNC(1) IF(i_am_accel_node .AND. acc_is_present(pfield3d))
!$ACC WAIT(1)

patch_2d => patch3D%p_patch_2d(1)

field_globsum=0._wp
DO jk =1,ks
  ptmp = global_sum_array( &
&   patch_2d%cells%area(:,:) * &
&   patch3d%wet_halo_zero_c(:,1,:) * &
&   seddw(jk) * sed_state(jk) * &
&   pfield3d(:,jk,:) &
 )
 field_globsum = field_globsum + ptmp
ENDDO

END SUBROUTINE


SUBROUTINE calc_inventory2d(patch3D, pfield2d, field_globsum, jk, ssh, pddpo)
! calculate inventory of 2d fields (optionally at given level)
! inv = (tracer * volume)
! volume = area* (cell thickness  + surface_height) : if ocean_state given
! volume = area* (cell thickness(level) ) : if ocean_state not given
REAL(wp), TARGET:: pfield2d(:,:)
TYPE(t_patch_3d), TARGET, INTENT(in) :: patch3D
REAL(wp), TARGET, OPTIONAL:: ssh(:,:)
REAL(wp), TARGET, OPTIONAL:: pddpo(:,:,:)
REAL(wp), INTENT(OUT) :: field_globsum
INTEGER,  OPTIONAL :: jk
! Local
TYPE(t_patch), POINTER :: patch_2d
INTEGER :: ik

!$ACC UPDATE HOST(pfield2d) ASYNC(1) IF(i_am_accel_node .AND. acc_is_present(pfield2d))
!$ACC WAIT(1)

patch_2d => patch3D%p_patch_2d(1)

ik = 1
IF(PRESENT(jk)) ik=jk

IF(PRESENT(ssh))THEN ! if needed add surface_height to cell thickness
field_globsum=global_sum_array(&
&   patch_2d%cells%area(:,:) * &
&   patch3d%wet_halo_zero_c(:,1,:) * &
&   (pddpo(:,1,:)+ssh)*&
&   pfield2d(:,:))
ELSEIF(ik < 0) THEN ! do not consider thickness
field_globsum=global_sum_array(&
&   patch_2d%cells%area(:,:) * &
&   patch3d%wet_halo_zero_c(:,1,:) * &
&   pfield2d(:,:))
ELSE
field_globsum=global_sum_array(&
&   patch_2d%cells%area(:,:) * &
&   patch3d%wet_halo_zero_c(:,ik,:) * &
&   pddpo(:,ik,:)*&
&   pfield2d(:,:))
ENDIF


END SUBROUTINE

END MODULE mo_hamocc_diagnostics

