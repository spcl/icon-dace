!
! This module contains the I/O routines for initicon
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

MODULE mo_initicon_utils

  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: msg_level, ntracer, iqv, iqc, iqi, iqr, iqs, iqg, iforcing, &
                                    iqh, iqnc, iqni, iqnr, iqns, iqng, iqnh, iqgl, iqhl
  USE mo_dynamics_config,     ONLY: nnow, nnow_rcf, nnew, nnew_rcf
  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_state,      ONLY: p_nh_state
  USE mo_nonhydro_types,      ONLY: t_nh_state, t_nh_metrics, t_nh_diag, t_nh_prog
  USE mo_nonhydrostatic_config, ONLY: kstart_moist, kstart_tracer
  USE mo_nwp_lnd_types,       ONLY: t_lnd_state, t_lnd_prog, t_lnd_diag, t_wtr_prog
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_initicon_types,      ONLY: t_initicon_state, alb_snow_var, t_pi_atm_in, t_pi_sfc_in, t_pi_atm, &
    &                               t_pi_sfc, t_sfc_inc, ana_varnames_dict, t_init_state_const
  USE mo_initicon_config,     ONLY: init_mode, l_sst_in, qcana_mode, qiana_mode, qrsgana_mode, &
    &                               ana_varnames_map_file, lread_vn, fire2d_filename,          &
    &                               lvert_remap_fg, aerosol_fg_present, icpl_da_sfcevap
  USE mo_impl_constants,      ONLY: MODE_DWDANA, MODE_IAU,                              &
                                    MODE_IAU_OLD, MODE_IFSANA, MODE_COMBINED,           &
    &                               MODE_COSMO, MODE_ICONVREMAP, MODIS, LSS_JSBACH,     &
    &                               min_rlcell_int, grf_bdywidth_c, min_rlcell,         &
    &                               iss, iorg, ibc, iso4, idu, SUCCESS, iaes
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_radiation_config,    ONLY: albedo_type
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_grid_config,         ONLY: n_dom, nroot
  USE mo_mpi,                 ONLY: my_process_is_stdio, p_io, &
    &                               p_comm_work, my_process_is_mpi_workroot, &
    &                               p_min, p_max, p_sum, num_work_procs, my_process_is_work
  USE mo_util_string,         ONLY: tolower
  USE mo_lnd_nwp_config,      ONLY: nlev_soil, ntiles_total, lseaice, llake, lmulti_snow,         &
    &                               isub_lake, frlnd_thrhld,             &
    &                               frlake_thrhld, frsea_thrhld, nlev_snow, ntiles_lnd,           &
    &                               l2lay_rho_snow, lprog_albsi, dzsoil, frsi_min
  USE mo_nwp_sfc_utils,       ONLY: init_snowtile_lists
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config, iprog_aero
  USE mo_aes_phy_config,      ONLY: aes_phy_config
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
  USE sfc_terra_data,         ONLY: csalb_snow_min, csalb_snow_max, csalb_snow, crhosmin_ml, crhosmax_ml, &
    &                               cpwp, cfcap
  USE mo_physical_constants,  ONLY: tf_salt, tmelt, cpd, rd, cvd_o_rd, p0ref, vtmpc1, ci, clw, cpv, cvd, cvv
  USE mo_hydro_adjust,        ONLY: hydro_adjust
  USE sfc_seaice,             ONLY: seaice_coldinit_nwp
  USE mo_post_op,             ONLY: perform_post_op
  USE mo_var_metadata_types,  ONLY: t_var_metadata, POST_OP_NONE
  USE mo_var_metadata,        ONLY: get_var_name
  USE mo_var_list_register,   ONLY: t_vl_register_iter
  USE mo_var,                 ONLY: level_type_ml
  USE sfc_flake,              ONLY: flake_coldinit
  USE mtime,                  ONLY: datetime, newDatetime, deallocateDatetime, &
    &                               OPERATOR(==), OPERATOR(+) 
  USE mo_intp_data_strc,      ONLY: p_int_state
  USE mo_dictionary,          ONLY: t_dictionary
  USE mo_checksum,            ONLY: printChecksum
  USE mo_fortran_tools,       ONLY: init
  USE mo_time_config,         ONLY: time_config
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights,         &
    &                                  calculate_time_interpolation_weights
  USE mo_aerosol_sources,     ONLY: inquire_fire2d_data
  USE mo_aerosol_sources_types,  ONLY: p_dust_source_const, p_fire_source_info
  USE mo_upatmo_config,       ONLY: upatmo_config
  USE mo_2mom_mcrph_util,     ONLY: set_qnc, set_qnr, set_qni,   &
    &                               set_qns, set_qng, set_qnh_expPSD_N0const


  IMPLICIT NONE


  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_initicon_utils'


  ! construct/destruct
  PUBLIC :: construct_initicon
  PUBLIC :: deallocate_initicon
  PUBLIC :: allocate_extana_atm
  PUBLIC :: allocate_extana_sfc

  PUBLIC :: initicon_inverse_post_op
  PUBLIC :: copy_initicon2prog_atm
  PUBLIC :: copy_initicon2prog_sfc
  PUBLIC :: copy_fg2initicon
  PUBLIC :: initVarnamesDict
  PUBLIC :: fill_tile_points
  PUBLIC :: init_snowtiles
  PUBLIC :: printChecksums
  PUBLIC :: init_aerosol
  PUBLIC :: init_qnx_from_qx_twomom
  PUBLIC :: init_qnxinc_from_qxinc_twomom
  PUBLIC :: get_diag_stat_comm_work
  PUBLIC :: new_land_from_ocean
  PUBLIC :: prepare_thermo_src_term

  CONTAINS




  !-------------
  !>
  !! SUBROUTINE initicon_inverse_post_op
  !! Perform inverse post_op on input field, if necessary 
  !!
  SUBROUTINE initicon_inverse_post_op(varname, optvar_out2D, optvar_out3D)
    CHARACTER(len=*), INTENT(IN)      :: varname             !< var name of field to be read
    REAL(wp), OPTIONAL, INTENT(INOUT) :: optvar_out2D(:,:)   !< 3D output field
    REAL(wp), OPTIONAL, INTENT(INOUT) :: optvar_out3D(:,:,:) !< 2D output field
    INTEGER                         :: i              ! loop count
    TYPE(t_var_metadata), POINTER   :: info           ! variable metadata
    CHARACTER(*), PARAMETER     :: routine = 'initicon_inverse_post_op'
    CHARACTER(LEN=LEN_TRIM(varname)) :: lc_varname
    TYPE(t_vl_register_iter) :: vl_iter

    !-------------------------------------------------------------------------
    ! Check consistency of optional arguments
    IF (PRESENT( optvar_out2D ) .AND. PRESENT( optvar_out3D )) THEN
      CALL finish(routine, 'Only one optional argument must be present')
    ELSE IF (.NOT.PRESENT( optvar_out2D ) .AND. .NOT.PRESENT( optvar_out3D )) THEN
      CALL finish(routine, 'One of 2 optional arguments must be present')
    ENDIF
    lc_varname = tolower(varname)
    ! get metadata information for field to be read
    NULLIFY(info)
    DO WHILE(vl_iter%next() .AND. .NOT.ASSOCIATED(info))
      ! loop only over model level variables
      IF (vl_iter%cur%p%vlevel_type /= level_type_ml) CYCLE 
      DO i = 1, vl_iter%cur%p%nvars
        info => vl_iter%cur%p%vl(i)%p%info
        IF (lc_varname == tolower(get_var_name(info))) EXIT
        NULLIFY(info)
      END DO
    ENDDO
    IF (.NOT.ASSOCIATED(info)) THEN
      CALL message(TRIM(varname)//' not found',message_text)
      CALL finish(routine, 'Varname does not match any of the ICON variable names')
    ENDIF
    ! perform post_op
    IF (info%post_op%ipost_op_type /= POST_OP_NONE) THEN
      IF(my_process_is_stdio() .AND. msg_level>10) &
        & CALL message(routine, 'Inverse Post_op for: ' // varname)
      IF (PRESENT(optvar_out2D)) THEN
        CALL perform_post_op(info%post_op, optvar_out2D, opt_inverse=.TRUE.)
      ELSE IF (PRESENT(optvar_out3D)) THEN
        CALL perform_post_op(info%post_op, optvar_out3D, opt_inverse=.TRUE.)
      ENDIF
    ENDIF
  END SUBROUTINE initicon_inverse_post_op

  !>
  !! SUBROUTINE init_aersosol
  !! Initializes the aerosol field from the climatology if no first-guess data are available
  !!
  SUBROUTINE init_aerosol(p_patch, ext_data, prm_diag)

    TYPE(t_patch),          INTENT(in)    :: p_patch(:)
    TYPE(t_external_data),  INTENT(in)    :: ext_data(:)
    TYPE(t_nwp_phy_diag),   INTENT(inout) :: prm_diag(:)
    
    CHARACTER(*), PARAMETER     :: routine = 'init_aerosol'

    TYPE(t_time_interpolation_weights)  :: current_time_interpolation_weights

    TYPE(datetime), POINTER :: mtime_hour
    
    INTEGER  :: rl_start, rl_end, i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: jb, jc, jg
    
    INTEGER  :: mo1, mo2
    REAL(wp) :: zw1, zw2

    !    CALL month2hour (time_config%cur_datetime, mo1, mo2, wgt)
    mtime_hour => newDatetime(time_config%tc_current_date)
    mtime_hour%time%minute = 0
    mtime_hour%time%second = 0
    mtime_hour%time%ms     = 0        
    current_time_interpolation_weights = calculate_time_interpolation_weights(mtime_hour)
    call deallocateDatetime(mtime_hour)
    mo1 = current_time_interpolation_weights%month1
    mo2 = current_time_interpolation_weights%month2
    zw1 = current_time_interpolation_weights%weight1
    zw2 = current_time_interpolation_weights%weight2
    
!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)
    DO jg = 1, n_dom

      rl_start = 1
      rl_end   = min_rlcell

      i_startblk = p_patch(jg)%cells%start_block(rl_start)
      i_endblk   = p_patch(jg)%cells%end_block(rl_end)

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        DO jc = i_startidx, i_endidx

          IF (.NOT. aerosol_fg_present(jg,iss)) prm_diag(jg)%aerosol(jc,iss,jb) =   &
            &  zw1*ext_data(jg)%atm_td%aer_ss(jc,jb,mo1)  + zw2*ext_data(jg)%atm_td%aer_ss(jc,jb,mo2)
          IF (.NOT. aerosol_fg_present(jg,iorg)) prm_diag(jg)%aerosol(jc,iorg,jb) = &
            &  zw1*ext_data(jg)%atm_td%aer_org(jc,jb,mo1) + zw2*ext_data(jg)%atm_td%aer_org(jc,jb,mo2)
          IF (.NOT. aerosol_fg_present(jg,ibc)) prm_diag(jg)%aerosol(jc,ibc,jb) =   &
            &  zw1*ext_data(jg)%atm_td%aer_bc(jc,jb,mo1)  + zw2*ext_data(jg)%atm_td%aer_bc(jc,jb,mo2)
          IF (.NOT. aerosol_fg_present(jg,iso4)) prm_diag(jg)%aerosol(jc,iso4,jb) = &
            &  zw1*ext_data(jg)%atm_td%aer_so4(jc,jb,mo1) + zw2*ext_data(jg)%atm_td%aer_so4(jc,jb,mo2)
          IF (.NOT. aerosol_fg_present(jg,idu)) prm_diag(jg)%aerosol(jc,idu,jb) =   &
            &  zw1*ext_data(jg)%atm_td%aer_dust(jc,jb,mo1)+ zw2*ext_data(jg)%atm_td%aer_dust(jc,jb,mo2)

        ENDDO

      ENDDO
!$OMP END DO
!$OMP MASTER
      IF ( ANY(.NOT.aerosol_fg_present(jg,1:5))) THEN
        WRITE(message_text,'(a,i3,a,i3,a)') 'Aerosol initialized from climatology, domain ',jg,', for',&
          COUNT(.NOT.aerosol_fg_present(jg,1:5)),' of 5 types'
        CALL message(routine, TRIM(message_text))
      ENDIF

      CALL p_dust_source_const(jg)%init(ext_data(jg)%atm%i_lc_shrub_eg,  &
        &                               ext_data(jg)%atm%i_lc_shrub,     &
        &                               ext_data(jg)%atm%i_lc_grass,     &
        &                               ext_data(jg)%atm%i_lc_bare_soil, &
        &                               ext_data(jg)%atm%i_lc_sparse,    &
        &                               i_st_sand=3, i_st_sandyloam=4,   &
        &                               i_st_loam=5, i_st_clayloam=6,    &
        &                               i_st_clay=7, nlu_classes=23,     &
        &                               soiltype_sidx=0, soiltype_eidx=9)
        
      IF (iprog_aero > 2) THEN
        CALL p_fire_source_info(jg)%init( ext_data(jg)%atm%bcfire,  &
          &                               ext_data(jg)%atm%ocfire,  &
          &                               ext_data(jg)%atm%so2fire  )
        CALL inquire_fire2d_data(p_patch(jg), nroot, fire2d_filename, p_fire_source_info(jg), &
          &                      time_config%tc_current_date, .TRUE.)
      ENDIF ! iprog_aero > 2
!$OMP END MASTER

    ENDDO
!$OMP END PARALLEL

  END SUBROUTINE init_aerosol


  !>
  !! SUBROUTINE fill_tile_points
  !! Used in the case of a 'cold' tile initialization
  !!  i.e. initializing a run with tiles with first guess data not containing tiles. The first guess data 
  !!  orignate from a run without tiles.
  !! or tile coldstart
  !!  i.e. initializing a run with tiles with first guess data not containing tiles. The first guess data
  !!  orignate from a run without tiles (but tile-averaged variables).
  !!  In the latter case the filling routine is only applied to the ANA fields fr_seaice and t_seasfc.
  !! 
  !! Specifically, this routine fills sub-grid scale (previously nonexistent) land and water points
  !! with appropriate data from neighboring grid points where possible
  !!
  !!
  SUBROUTINE fill_tile_points(p_patch, p_lnd_state, ext_data, process_ana_vars)

    TYPE(t_patch),             INTENT(IN)    :: p_patch(:)
    TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_external_data),     INTENT(INOUT) :: ext_data(:)

    LOGICAL,                   INTENT(IN)    :: process_ana_vars  ! neighbour filling only for analysed fields
                                                                  ! fr_seaice and t_seasfc

    TYPE(t_lnd_prog),  POINTER :: lnd_prog
    TYPE(t_lnd_diag),  POINTER :: lnd_diag
    TYPE(t_wtr_prog),  POINTER :: wtr_prog

    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

    INTEGER :: jg, jb, jk, jc, jt, ji
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx

    REAL(wp), DIMENSION(nproma,10) :: lpmask, fpmask, spmask
    REAL(wp), DIMENSION(nproma)    :: lpcount, fpcount, spcount
    REAL(wp) :: wgt(10), fr_lnd, fr_lk, fr_sea, th_notile, aux_lk(nproma,12)
    REAL(wp), ALLOCATABLE :: frac_ml_lk(:,:)

    !-------------------------------------------------------------------------

    th_notile = 0.5_wp

    ! Weighting factors for neighbor filling
    DO ji = 2, 10
      SELECT CASE (ji)
      CASE(2,5,8) ! direct neighbors
        wgt(ji) = 1._wp
      CASE DEFAULT ! indirect neighbors
        wgt(ji) = 0.5_wp
      END SELECT
    ENDDO

    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      lnd_prog => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
      lnd_diag => p_lnd_state(jg)%diag_lnd
      wtr_prog => p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))

      iidx     => p_int_state(jg)%rbf_c2grad_idx
      iblk     => p_int_state(jg)%rbf_c2grad_blk

      i_startblk = p_patch(jg)%cells%start_block(grf_bdywidth_c+1)
      i_endblk   = p_patch(jg)%cells%end_block(min_rlcell_int)

      ! Compute ratio between mixed-layer depth and lake depth
      ALLOCATE(frac_ml_lk(nproma,p_patch(jg)%nblks_c))
      WHERE(ext_data(jg)%atm%depth_lk(:,:) > 0.0_wp)
        frac_ml_lk(:,:) = wtr_prog%h_ml_lk(:,:) / MAX(0.05_wp,ext_data(jg)%atm%depth_lk(:,:))
      ELSEWHERE
        frac_ml_lk(:,:) = 0._wp
      END WHERE

!$OMP PARALLEL
!$OMP DO PRIVATE(ji,jb,jk,jc,jt,i_startidx,i_endidx,lpmask,fpmask,spmask,lpcount,fpcount,spcount,fr_lnd,fr_lk,fr_sea,aux_lk)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell_int)


        ! call flake_coldinit and store results on a local auxiliary array; they are used as a backup if no
        ! appropriate neighbor points are found
        CALL flake_coldinit(                                     &
          &     nflkgb      = ext_data(jg)%atm%list_lake%ncount(jb),&  ! in
          &     idx_lst_fp  = ext_data(jg)%atm%list_lake%idx(:,jb), &  ! in
          &     depth_lk    = ext_data(jg)%atm%depth_lk  (:,jb),    &  ! in
          &     tskin       = lnd_prog%t_g_t(:,jb,isub_lake)   ,    &  ! in
          &     t_snow_lk_p = aux_lk(:,1),                          &
          &     h_snow_lk_p = aux_lk(:,2),                          &
          &     t_ice_p     = aux_lk(:,3),                          &
          &     h_ice_p     = aux_lk(:,4),                          &
          &     t_mnw_lk_p  = aux_lk(:,5),                          &
          &     t_wml_lk_p  = aux_lk(:,6),                          &
          &     t_bot_lk_p  = aux_lk(:,7),                          &
          &     c_t_lk_p    = aux_lk(:,8),                          &
          &     h_ml_lk_p   = aux_lk(:,9),                          &
          &     t_b1_lk_p   = aux_lk(:,10),                         &
          &     h_b1_lk_p   = aux_lk(:,11),                         &
          &     t_g_lk_p    = aux_lk(:,12)                          ) 


        lpmask(:,:) = 0._wp
        fpmask(:,:) = 0._wp
        spmask(:,:) = 0._wp

        lpcount(:) = 0._wp
        fpcount(:) = 0._wp
        spcount(:) = 0._wp

        ! Prepare control variables to determine for which grid points neighbor filling is possible
        DO jc = i_startidx, i_endidx

          fr_lnd = ext_data(jg)%atm%fr_land(jc,jb)
          IF (fr_lnd > frlnd_thrhld .AND. fr_lnd < th_notile) THEN
            lpmask(jc,1) = 1._wp
            DO ji = 2,10
              fr_lnd = ext_data(jg)%atm%fr_land(iidx(ji,jc,jb),iblk(ji,jc,jb))
              IF (fr_lnd > th_notile) THEN
                lpmask(jc,ji) = wgt(ji)
                lpcount(jc) = lpcount(jc)+wgt(ji)
              ENDIF
            ENDDO
            IF (lpcount(jc) == 0._wp) lpmask(jc,1) = 0._wp
          ENDIF

          fr_lk = ext_data(jg)%atm%fr_lake(jc,jb)
          IF (fr_lk > frlake_thrhld .AND. fr_lk < th_notile) THEN
            fpmask(jc,1) = 1._wp
            DO ji = 2,10
              fr_lk = ext_data(jg)%atm%fr_lake(iidx(ji,jc,jb),iblk(ji,jc,jb))
              IF (fr_lk > th_notile) THEN
                fpmask(jc,ji) = wgt(ji)
                fpcount(jc) = fpcount(jc)+wgt(ji)
              ENDIF
            ENDDO
            IF (fpcount(jc) == 0._wp) THEN
              fpmask(jc,1) = 0._wp
              ! use backup values from flake_coldinit if no neighbors are available
              IF (.NOT. process_ana_vars) THEN
                wtr_prog%t_mnw_lk(jc,jb) = aux_lk(jc,5)
                wtr_prog%t_wml_lk(jc,jb) = aux_lk(jc,6)
                wtr_prog%t_bot_lk(jc,jb) = aux_lk(jc,7)
                wtr_prog%c_t_lk  (jc,jb) = aux_lk(jc,8)
                wtr_prog%h_ml_lk (jc,jb) = aux_lk(jc,9)
                wtr_prog%t_b1_lk (jc,jb) = aux_lk(jc,10)
                wtr_prog%h_b1_lk (jc,jb) = aux_lk(jc,11)
              ENDIF
            ENDIF
          ENDIF

          fr_sea = 1._wp - (ext_data(jg)%atm%fr_lake(jc,jb)+ext_data(jg)%atm%fr_land(jc,jb))
          IF (fr_sea > frsea_thrhld .AND. fr_sea < th_notile) THEN
            spmask(jc,1) = 1._wp
            DO ji = 2,10
              fr_sea = 1._wp - (ext_data(jg)%atm%fr_lake(iidx(ji,jc,jb),iblk(ji,jc,jb)) + &
                                ext_data(jg)%atm%fr_land(iidx(ji,jc,jb),iblk(ji,jc,jb))   )
              IF (fr_sea > th_notile) THEN
                spmask(jc,ji) = wgt(ji)
                spcount(jc) = spcount(jc)+wgt(ji)
              ENDIF
            ENDDO
            IF (spcount(jc) == 0._wp) spmask(jc,1) = 0._wp
          ENDIF

        ENDDO

        ! Apply neighbor filling
        !
        IF (process_ana_vars) THEN

          DO jc = i_startidx, i_endidx

            ! a) ocean points
            IF (spmask(jc,1) == 1._wp) THEN
              CALL ngb_search(lnd_diag%fr_seaice, iidx, iblk, spmask, spcount, jc, jb)
              CALL ngb_search(lnd_diag%t_seasfc, iidx, iblk, spmask, spcount, jc, jb)
            ENDIF

          ENDDO  ! jc

        ELSE   ! .NOT. process_ana_vars

          DO jc = i_startidx, i_endidx
            ! a) ocean points
            IF (spmask(jc,1) == 1._wp) THEN
              CALL ngb_search(wtr_prog%t_ice, iidx, iblk, spmask, spcount, jc, jb)
              CALL ngb_search(wtr_prog%h_ice, iidx, iblk, spmask, spcount, jc, jb)
              IF ( lprog_albsi ) THEN
                CALL ngb_search(wtr_prog%alb_si, iidx, iblk, spmask, spcount, jc, jb)
              ENDIF
            ENDIF

            ! b) lake points
            IF (fpmask(jc,1) == 1._wp) THEN
              CALL ngb_search(wtr_prog%t_mnw_lk, iidx, iblk, fpmask, fpcount, jc, jb)
              CALL ngb_search(wtr_prog%t_wml_lk, iidx, iblk, fpmask, fpcount, jc, jb)
              CALL ngb_search(frac_ml_lk,        iidx, iblk, fpmask, fpcount, jc, jb)
              CALL ngb_search(wtr_prog%t_bot_lk, iidx, iblk, fpmask, fpcount, jc, jb)
              CALL ngb_search(wtr_prog%c_t_lk,   iidx, iblk, fpmask, fpcount, jc, jb)
              CALL ngb_search(wtr_prog%t_b1_lk,  iidx, iblk, fpmask, fpcount, jc, jb)
              CALL ngb_search(wtr_prog%h_b1_lk,  iidx, iblk, fpmask, fpcount, jc, jb)
              ! restore mixed-layer depth
              wtr_prog%h_ml_lk(jc,jb) = frac_ml_lk(jc,jb) * ext_data(jg)%atm%depth_lk(jc,jb)
            ENDIF
          ENDDO  ! jc

          ! c) land points
          DO jt = 1, ntiles_total

            ! single-layer fields
            DO jc = i_startidx, i_endidx
              IF (lpmask(jc,1) == 1._wp) THEN
                CALL ngb_search(lnd_diag%freshsnow_t(:,:,jt), iidx, iblk, lpmask, lpcount, jc, jb)
                CALL ngb_search(lnd_prog%w_snow_t   (:,:,jt), iidx, iblk, lpmask, lpcount, jc, jb)
                CALL ngb_search(lnd_prog%w_i_t      (:,:,jt), iidx, iblk, lpmask, lpcount, jc, jb)
                CALL ngb_search(lnd_diag%h_snow_t   (:,:,jt), iidx, iblk, lpmask, lpcount, jc, jb)
                CALL ngb_search(lnd_prog%t_snow_t   (:,:,jt), iidx, iblk, lpmask, lpcount, jc, jb)
                CALL ngb_search(lnd_prog%rho_snow_t (:,:,jt), iidx, iblk, lpmask, lpcount, jc, jb)

                CALL ngb_search(lnd_prog%t_so_t(:,nlev_soil+1,:,jt), iidx, iblk, lpmask, lpcount, jc, jb)

                IF (l2lay_rho_snow) THEN ! only rho_snow layer 1 is relevant in this case
                  CALL ngb_search(lnd_prog%rho_snow_mult_t(:,1,:,jt), iidx, iblk, lpmask, lpcount, jc, jb)
                ENDIF
              ENDIF
            ENDDO

            ! soil fields
            DO jk = 1, nlev_soil
              DO jc = i_startidx, i_endidx
                IF (lpmask(jc,1) == 1._wp) THEN
                  CALL ngb_search(lnd_prog%t_so_t(:,jk,:,jt),     iidx, iblk, lpmask, lpcount, jc, jb)
                  CALL ngb_search(lnd_prog%w_so_t(:,jk,:,jt),     iidx, iblk, lpmask, lpcount, jc, jb)
                  CALL ngb_search(lnd_prog%w_so_ice_t(:,jk,:,jt), iidx, iblk, lpmask, lpcount, jc, jb)
                ENDIF
              ENDDO
            ENDDO

            IF (lmulti_snow) THEN ! multi-layer snow fields
              DO jk = 1, nlev_snow
                DO jc = i_startidx, i_endidx
                  IF (lpmask(jc,1) == 1._wp) THEN
                    CALL ngb_search(lnd_prog%t_snow_mult_t(:,jk,:,jt),   iidx, iblk, lpmask, lpcount, jc, jb)
                    CALL ngb_search(lnd_prog%rho_snow_mult_t(:,jk,:,jt), iidx, iblk, lpmask, lpcount, jc, jb)
                    CALL ngb_search(lnd_prog%wtot_snow_t(:,jk,:,jt),     iidx, iblk, lpmask, lpcount, jc, jb)
                    CALL ngb_search(lnd_prog%wliq_snow_t(:,jk,:,jt),     iidx, iblk, lpmask, lpcount, jc, jb)
                    CALL ngb_search(lnd_prog%dzh_snow_t(:,jk,:,jt),      iidx, iblk, lpmask, lpcount, jc, jb)

                    IF (jk == 1) CALL ngb_search(lnd_prog%t_snow_mult_t(:,nlev_snow+1,:,jt), iidx, iblk, lpmask, lpcount, jc, jb)
                  ENDIF
                ENDDO  ! jc
              ENDDO
            ENDIF  ! lmulti_snow

          ENDDO  ! jt

        ENDIF  ! process_ana_vars

      ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      DEALLOCATE(frac_ml_lk)

    ENDDO  ! jg

    CONTAINS

    SUBROUTINE ngb_search (fld, iidx, iblk, mask, cnt, jc, jb)

      REAL(wp), INTENT(INOUT) :: fld(:,:)
      REAL(wp), INTENT(IN)    :: mask(:,:), cnt(:)
      INTEGER,  INTENT(IN)    :: iidx(:,:,:), iblk(:,:,:), jc, jb

      INTEGER :: ji

      fld(jc,jb) = 0._wp
      DO ji = 2, 10
        fld(jc,jb) = fld(jc,jb) + fld(iidx(ji,jc,jb),iblk(ji,jc,jb))*mask(jc,ji)
      ENDDO
      fld(jc,jb) = fld(jc,jb)/cnt(jc)

    END SUBROUTINE ngb_search

  END SUBROUTINE fill_tile_points

  !>
  !! SUBROUTINE init_snowtiles
  !! Active in the case of a tile warmstart in combination with snowtiles.
  !! In this case, the tile-based index lists and the tile fractions (frac_t) need to be restored
  !! from the landuse-class fractions and the snow-cover fractions
  !!
  SUBROUTINE init_snowtiles(p_patch, p_lnd_state, ext_data)

    TYPE(t_patch),             INTENT(IN)    :: p_patch(:)
    TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_external_data),     INTENT(INOUT) :: ext_data(:)

    TYPE(t_lnd_diag),  POINTER :: lnd_diag

    INTEGER :: jg, jt

    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      lnd_diag => p_lnd_state(jg)%diag_lnd

      ! initialize snowfrac_t with appropriate values
      DO jt = ntiles_lnd+1, ntiles_total
        WHERE (lnd_diag%snowfrac_lc_t(:,:,jt) > 0._wp) 
          lnd_diag%snowfrac_t(:,:,jt) = 1._wp
        ELSEWHERE
          lnd_diag%snowfrac_t(:,:,jt) = 0._wp
        END WHERE
      ENDDO

      CALL init_snowtile_lists(p_patch(jg), ext_data(jg), lnd_diag)

    ENDDO

  END SUBROUTINE init_snowtiles

  !>
  !! SUBROUTINE copy_initicon2prog_atm
  !! Copies atmospheric fields interpolated by init_icon to the
  !! prognostic model state variables 
  !!
  !! Required input: initicon state
  !! Output is written on fields of NH state
  !!
  SUBROUTINE copy_initicon2prog_atm(p_patch, initicon, p_nh_state)

    TYPE(t_patch),          INTENT(IN) :: p_patch(:)
    TYPE(t_initicon_state), INTENT(IN) :: initicon(:)

    TYPE(t_nh_state),      INTENT(INOUT) :: p_nh_state(:)

    INTEGER :: jg, jb, jk, jc, je, idx, itracer
    INTEGER :: nblks_c, npromz_c, nblks_e, npromz_e, nlen, nlev, nlevp1, ntl, ntlr

!$OMP PARALLEL PRIVATE(jg,nblks_c,npromz_c,nblks_e,npromz_e,nlev,nlevp1,ntl,ntlr)
    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      nblks_c   = p_patch(jg)%nblks_c
      npromz_c  = p_patch(jg)%npromz_c
      nblks_e   = p_patch(jg)%nblks_e
      npromz_e  = p_patch(jg)%npromz_e
      nlev      = p_patch(jg)%nlev
      nlevp1    = p_patch(jg)%nlevp1
      ntl       = nnow(jg)
      ntlr      = nnow_rcf(jg)

!$OMP DO PRIVATE(jb,jk,je,nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_e

        IF (jb /= nblks_e) THEN
          nlen = nproma
        ELSE
          nlen = npromz_e
        ENDIF

        ! Wind speed
        DO jk = 1, nlev
          DO je = 1, nlen
            p_nh_state(jg)%prog(ntl)%vn(je,jk,jb) = initicon(jg)%atm%vn(je,jk,jb)
          ENDDO
        ENDDO

      ENDDO
!$OMP END DO

!$OMP DO PRIVATE(jb,jk,jc,nlen,idx,itracer) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_c

        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF

        ! 3D fields
        DO jk = 1, nlev
          DO jc = 1, nlen
            !
            ! Diagnostic pressure and temperature
            p_nh_state(jg)%diag%pres(jc,jk,jb)         = initicon(jg)%atm%pres(jc,jk,jb)
            p_nh_state(jg)%diag%temp(jc,jk,jb)         = initicon(jg)%atm%temp(jc,jk,jb)
            !
            ! Dynamic prognostic variables on cell points
            p_nh_state(jg)%prog(ntl)%w(jc,jk,jb)       = initicon(jg)%atm%w(jc,jk,jb)
            p_nh_state(jg)%prog(ntl)%theta_v(jc,jk,jb) = initicon(jg)%atm%theta_v(jc,jk,jb)
            p_nh_state(jg)%prog(ntl)%exner(jc,jk,jb)   = initicon(jg)%atm%exner(jc,jk,jb)
            p_nh_state(jg)%prog(ntl)%rho(jc,jk,jb)     = initicon(jg)%atm%rho(jc,jk,jb)

            ! Water vapor
            p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqv) = initicon(jg)%atm%qv(jc,jk,jb)
          ENDDO
        ENDDO

        ! Cloud and precipitation hydrometeors - these are supposed to be zero in the region where
        ! moisture physics is turned off
        !
        ! above kstart_moist(jg): set to zero
        DO jk = 1, kstart_moist(jg)-1
          DO jc = 1, nlen
            p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqc) = 0.0_wp
            p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqi) = 0.0_wp
            IF ( iqr /= 0 .AND. iqr <= ntracer) THEN
              p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqr) = 0.0_wp
            END IF
            IF ( iqs /= 0 .AND. iqs <= ntracer) THEN
              p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqs) = 0.0_wp
            END IF
            IF ( (atm_phy_nwp_config(jg)%lhave_graupel) .OR. ( iqg /= 0 .AND. iqg <= ntracer) ) THEN
              p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqg) = 0.0_wp
            END IF
            IF ( atm_phy_nwp_config(jg)%l2moment .OR. aes_phy_config(jg)%l2moment) THEN
              p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqh)  = 0.0_wp
              p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqnc) = 0.0_wp
              p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqni) = 0.0_wp
              p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqnr) = 0.0_wp
              p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqns) = 0.0_wp
              p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqng) = 0.0_wp
              p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqnh) = 0.0_wp
            END IF
          ENDDO
        ENDDO
        !
        IF (iforcing == iaes) THEN
          ! at and below kstart_moist(jg): copy from initicon%atm or set to zero
          ! HAS TO BE CHECKED for possibility to initialize
          DO jk = kstart_moist(jg), nlev
            DO jc = 1, nlen
              p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqc) = initicon(jg)%atm%qc(jc,jk,jb)
              p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqi) = initicon(jg)%atm%qi(jc,jk,jb)
              IF ( iqr /= 0 .AND. iqr <= ntracer) THEN
                p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqr) = 0.0_wp
              END IF
              IF ( iqs /= 0 .AND. iqs <= ntracer) THEN
                p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqs) = 0.0_wp
              END IF
              IF ( iqg /= 0 .AND. iqg <= ntracer) THEN
                ! as qg is not in atm initialize with zero
                p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqg) = 0.0_wp
              END IF
              IF ( aes_phy_config(jg)%l2moment) THEN
                p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqh)  = 0.0_wp
                p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqnc) = 0.0_wp
                p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqni) = 0.0_wp
                p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqnr) = 0.0_wp
                p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqns) = 0.0_wp
                p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqng) = 0.0_wp
                p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqnh) = 0.0_wp
              END IF
            ENDDO
          ENDDO
        ELSE
          ! at and below kstart_moist(jg): copy from initicon%atm
          DO jk = kstart_moist(jg), nlev
            DO jc = 1, nlen
              p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqc) = initicon(jg)%atm%qc(jc,jk,jb)
              p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqi) = initicon(jg)%atm%qi(jc,jk,jb)
              IF ( iqr /= 0 .AND. iqr <= ntracer) THEN
                p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqr) = initicon(jg)%atm%qr(jc,jk,jb)
              END IF
              IF ( iqs /= 0 .AND. iqs <= ntracer) THEN
                p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqs) = initicon(jg)%atm%qs(jc,jk,jb)
              END IF
              IF ( atm_phy_nwp_config(jg)%lhave_graupel ) THEN
                p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqg) = initicon(jg)%atm%qg(jc,jk,jb)
              END IF
              IF ( atm_phy_nwp_config(jg)%l2moment ) THEN
                p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqh)  = initicon(jg)%atm%qh(jc,jk,jb)
                p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqnc) = initicon(jg)%atm%qnc(jc,jk,jb)
                p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqni) = initicon(jg)%atm%qni(jc,jk,jb)
                p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqnr) = initicon(jg)%atm%qnr(jc,jk,jb)
                p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqns) = initicon(jg)%atm%qns(jc,jk,jb)
                p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqng) = initicon(jg)%atm%qng(jc,jk,jb)
                p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqnh) = initicon(jg)%atm%qnh(jc,jk,jb)
              END IF
            ENDDO
          ENDDO
        END IF

        ! ... additional tracer variables
        DO idx = 1, ntracer
          IF ( ASSOCIATED(initicon(jg)%atm%tracer(idx)%field) ) THEN
            itracer = initicon(jg)%atm%tracer(idx)%var_element%info%ncontained
            ! above kstart_tracer(jg,itracer): set to zero
            DO jk = 1, kstart_tracer(jg,itracer)-1
              DO jc = 1, nlen
                p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,itracer) = 0.0_wp
              ENDDO
            ENDDO
            ! at and below kstart_tracer(jg,itracer): copy from initicon%atm
            DO jk = kstart_tracer(jg,itracer), nlev
              DO jc = 1, nlen
                p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,itracer) =  &
                  &  initicon(jg)%atm%tracer(idx)%field(jc,jk,jb)
              ENDDO
            ENDDO
          END IF
        ENDDO

        ! w at surface level
        DO jc = 1, nlen
          p_nh_state(jg)%prog(ntl)%w(jc,nlevp1,jb)      = initicon(jg)%atm%w(jc,nlevp1,jb)
          p_nh_state(jg)%prog(nnew(jg))%w(jc,nlevp1,jb) = initicon(jg)%atm%w(jc,nlevp1,jb)
        ENDDO

        IF (init_mode == MODE_ICONVREMAP .OR. lvert_remap_fg) THEN ! copy also TKE field
          DO jk = 1, nlevp1
            DO jc = 1, nlen
              p_nh_state(jg)%prog(ntlr)%tke(jc,jk,jb) = initicon(jg)%atm%tke(jc,jk,jb)
            ENDDO
          ENDDO
        ENDIF

      ENDDO  ! jb
!$OMP END DO NOWAIT

    ENDDO  ! jg
!$OMP END PARALLEL

    ! Finally, compute exact hydrostatic adjustment for thermodynamic fields
    ! (in case of an upper-atmosphere extrapolation, this is already part 
    ! of the pressure "extrapolation" in 'src/atm_dyn_iconam/mo_nh_vert_interp: vert_interp')
    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE
      ntl = nnow(jg)

      IF (.NOT. upatmo_config(jg)%exp%l_expol) THEN
        CALL hydro_adjust(p_patch(jg), p_nh_state(jg)%metrics,                                  &
                          p_nh_state(jg)%prog(ntl)%rho,     p_nh_state(jg)%prog(ntl)%exner,     &
                          p_nh_state(jg)%prog(ntl)%theta_v )
      ENDIF

    ENDDO

  END SUBROUTINE copy_initicon2prog_atm


  !>
  !! SUBROUTINE copy_fg2initicon
  !! Copies first-guess fields from the assimilation cycle to the initicon state
  !! in order to prepare subsequent vertical remapping
  !!
  SUBROUTINE copy_fg2initicon(p_patch, initicon, p_nh_state)

    TYPE(t_patch),          INTENT(IN) :: p_patch(:)
    TYPE(t_nh_state),       INTENT(IN) :: p_nh_state(:)

    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::copy_fg2initicon"
    INTEGER :: jg, jb, jk, jc, je, ierrstat, idx, itracer
    INTEGER :: nblks_c, npromz_c, nblks_e, npromz_e, nlen, nlev, nlevp1, ntl, ntlr
    REAL(wp), ALLOCATABLE :: w_ifc(:,:,:), tke_ifc(:,:,:)

    REAL(wp) :: exner, tempv

    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      nblks_c   = p_patch(jg)%nblks_c
      npromz_c  = p_patch(jg)%npromz_c
      nblks_e   = p_patch(jg)%nblks_e
      npromz_e  = p_patch(jg)%npromz_e
      nlev      = p_patch(jg)%nlev
      nlevp1    = p_patch(jg)%nlevp1
      ntl       = nnow(jg)
      ntlr      = nnow_rcf(jg)

      ! allocate temporary array:
      ALLOCATE(w_ifc(nproma,    (nlev+1), nblks_c), &
        &      tke_ifc(nproma,  (nlev+1), nblks_c), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_e

        nlen = MERGE( nproma, npromz_e, jb /= nblks_e)

        ! Wind speed
        DO jk = 1, nlev
          DO je = 1, nlen
            initicon(jg)%atm_in%vn(je,jk,jb) = p_nh_state(jg)%prog(ntl)%vn(je,jk,jb)
          ENDDO
        ENDDO

      ENDDO
!$OMP END DO

!$OMP DO PRIVATE(jb,jk,jc,nlen,exner,tempv,idx,itracer) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_c

        nlen = MERGE(nproma, npromz_c, jb /= nblks_c)

        ! 3D fields
        DO jk = 1, nlev
          DO jc = 1, nlen

            ! Dynamic prognostic variables on cell points
            w_ifc(jc,jk,jb)                       = p_nh_state(jg)%prog(ntl)%w(jc,jk,jb)
            tke_ifc(jc,jk,jb)                     = p_nh_state(jg)%prog(ntlr)%tke(jc,jk,jb)
            initicon(jg)%atm_in%theta_v(jc,jk,jb) = p_nh_state(jg)%prog(ntl)%theta_v(jc,jk,jb)
            initicon(jg)%atm_in%rho(jc,jk,jb)     = p_nh_state(jg)%prog(ntl)%rho(jc,jk,jb)

            ! Moisture variables
            initicon(jg)%atm_in%qv(jc,jk,jb) = p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqv)
            initicon(jg)%atm_in%qc(jc,jk,jb) = p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqc)
            initicon(jg)%atm_in%qi(jc,jk,jb) = p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqi)
            initicon(jg)%atm_in%qr(jc,jk,jb) = p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqr)
            initicon(jg)%atm_in%qs(jc,jk,jb) = p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqs)
          ENDDO
        ENDDO

        ! ... additional tracer variables
        DO idx = 1, ntracer
          IF ( ASSOCIATED(initicon(jg)%atm_in%tracer(idx)%field) ) THEN
            itracer = initicon(jg)%atm_in%tracer(idx)%var_element%info%ncontained
            DO jk = 1, nlev
              DO jc = 1, nlen
                initicon(jg)%atm_in%tracer(idx)%field(jc,jk,jb) =  &
                  &  p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,itracer)
              ENDDO
            ENDDO
          END IF
        ENDDO

        IF (atm_phy_nwp_config(jg)%lhave_graupel) THEN
          DO jk = 1, nlev
            DO jc = 1, nlen
              initicon(jg)%atm_in%qg(jc,jk,jb) = p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqg)
            ENDDO
          ENDDO
        ELSE
          ! Probably unnecessary due to previous initialization?
          DO jk = 1, nlev
            DO jc = 1, nlen
              initicon(jg)%atm_in%qg(jc,jk,jb) = 0.0_wp
            ENDDO
          ENDDO          
        END IF
        
        ! 2-moment hydrometeors
        IF (atm_phy_nwp_config(jg)%l2moment) THEN
          DO jk = 1, nlev
            DO jc = 1, nlen
              initicon(jg)%atm_in%qh(jc,jk,jb)  = p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqh)
              initicon(jg)%atm_in%qnc(jc,jk,jb) = p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqnc)
              initicon(jg)%atm_in%qni(jc,jk,jb) = p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqni)
              initicon(jg)%atm_in%qnr(jc,jk,jb) = p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqnr)
              initicon(jg)%atm_in%qns(jc,jk,jb) = p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqns)
              initicon(jg)%atm_in%qng(jc,jk,jb) = p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqng)
              initicon(jg)%atm_in%qnh(jc,jk,jb) = p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqnh)
            ENDDO
          ENDDO
        END IF
          
        ! w and TKE at surface level
        DO jc = 1, nlen
          w_ifc(jc,nlevp1,jb)   = p_nh_state(jg)%prog(ntl)%w(jc,nlevp1,jb)
          tke_ifc(jc,nlevp1,jb) = p_nh_state(jg)%prog(ntlr)%tke(jc,nlevp1,jb)
        ENDDO

        ! diagnose pressure and temperature 
        IF (atm_phy_nwp_config(jg)%l2moment) THEN
          DO jk = 1, nlev
            DO jc = 1, nlen

              initicon(jg)%atm_in%w(jc,jk,jb) = (w_ifc(jc,jk,jb) + w_ifc(jc,jk+1,jb)) * 0.5_wp
              initicon(jg)%atm_in%tke(jc,jk,jb) = (tke_ifc(jc,jk,jb) + tke_ifc(jc,jk+1,jb)) * 0.5_wp

              exner = (initicon(jg)%atm_in%rho(jc,jk,jb)*initicon(jg)%atm_in%theta_v(jc,jk,jb)*rd/p0ref)**(1._wp/cvd_o_rd)
              tempv = initicon(jg)%atm_in%theta_v(jc,jk,jb)*exner

              initicon(jg)%atm_in%pres(jc,jk,jb) = exner**(cpd/rd)*p0ref
              initicon(jg)%atm_in%temp(jc,jk,jb) = tempv / (1._wp + vtmpc1*initicon(jg)%atm_in%qv(jc,jk,jb) - &
                (initicon(jg)%atm_in%qc(jc,jk,jb) + initicon(jg)%atm_in%qi(jc,jk,jb) +                        &
                 initicon(jg)%atm_in%qr(jc,jk,jb) + initicon(jg)%atm_in%qs(jc,jk,jb) +                        &
                 initicon(jg)%atm_in%qg(jc,jk,jb) + initicon(jg)%atm_in%qh(jc,jk,jb)   ))

            ENDDO
          ENDDO
        ELSE IF (atm_phy_nwp_config(jg)%lhave_graupel) THEN
          DO jk = 1, nlev
            DO jc = 1, nlen

              initicon(jg)%atm_in%w(jc,jk,jb) = (w_ifc(jc,jk,jb) + w_ifc(jc,jk+1,jb)) * 0.5_wp
              initicon(jg)%atm_in%tke(jc,jk,jb) = (tke_ifc(jc,jk,jb) + tke_ifc(jc,jk+1,jb)) * 0.5_wp

              exner = (initicon(jg)%atm_in%rho(jc,jk,jb)*initicon(jg)%atm_in%theta_v(jc,jk,jb)*rd/p0ref)**(1._wp/cvd_o_rd)
              tempv = initicon(jg)%atm_in%theta_v(jc,jk,jb)*exner

              initicon(jg)%atm_in%pres(jc,jk,jb) = exner**(cpd/rd)*p0ref
              initicon(jg)%atm_in%temp(jc,jk,jb) = tempv / (1._wp + vtmpc1*initicon(jg)%atm_in%qv(jc,jk,jb) - &
                (initicon(jg)%atm_in%qc(jc,jk,jb) + initicon(jg)%atm_in%qi(jc,jk,jb) +                        &
                 initicon(jg)%atm_in%qr(jc,jk,jb) + initicon(jg)%atm_in%qs(jc,jk,jb) +                        &
                 initicon(jg)%atm_in%qg(jc,jk,jb)                                      ))

            ENDDO
          ENDDO
        ELSE
          DO jk = 1, nlev
            DO jc = 1, nlen

              initicon(jg)%atm_in%w(jc,jk,jb) = (w_ifc(jc,jk,jb) + w_ifc(jc,jk+1,jb)) * 0.5_wp
              initicon(jg)%atm_in%tke(jc,jk,jb) = (tke_ifc(jc,jk,jb) + tke_ifc(jc,jk+1,jb)) * 0.5_wp

              exner = (initicon(jg)%atm_in%rho(jc,jk,jb)*initicon(jg)%atm_in%theta_v(jc,jk,jb)*rd/p0ref)**(1._wp/cvd_o_rd)
              tempv = initicon(jg)%atm_in%theta_v(jc,jk,jb)*exner

              initicon(jg)%atm_in%pres(jc,jk,jb) = exner**(cpd/rd)*p0ref
              initicon(jg)%atm_in%temp(jc,jk,jb) = tempv / (1._wp + vtmpc1*initicon(jg)%atm_in%qv(jc,jk,jb) - &
                (initicon(jg)%atm_in%qc(jc,jk,jb) + initicon(jg)%atm_in%qi(jc,jk,jb) +                        &
                 initicon(jg)%atm_in%qr(jc,jk,jb) + initicon(jg)%atm_in%qs(jc,jk,jb)   ))

            ENDDO
          ENDDO

        END IF

      ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ! cleanup
      IF (ALLOCATED(w_ifc))    DEALLOCATE(w_ifc)
      IF (ALLOCATED(tke_ifc))  DEALLOCATE(tke_ifc)

    ENDDO  ! jg

    ! Tell the vertical interpolation routine that vn needs to be processed
    lread_vn = .TRUE.

  END SUBROUTINE copy_fg2initicon


  !-------------
  !>
  !! SUBROUTINE copy_initicon2prog_sfc
  !! Copies surface fields interpolated by init_icon to the prognostic model 
  !! state variables. 
  !!
  !! Required input: initicon state
  !! Output is written on fields of land state
  !!
  SUBROUTINE copy_initicon2prog_sfc(p_patch, initicon, p_lnd_state, ext_data)

    ! Procedure arguments

    TYPE(t_patch),          INTENT(IN) :: p_patch(:)
    TYPE(t_initicon_state), INTENT(IN) :: initicon(:)

    TYPE(t_lnd_state),     INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_external_data), INTENT(   IN) :: ext_data(:)

    ! Local variables

    INTEGER  :: jg, jb, jc, jt, js, jp, ic, ilu
    INTEGER  :: nblks_c, npromz_c, nlen
    REAL(wp) :: zfrice_thrhld, zminsnow_alb, zmaxsnow_alb, zsnowalb_lu, t_fac

    ! Local arrays

    REAL(wp), DIMENSION(nproma) ::             &
                                &  frsi_in   , &  !< sea-ice fraction [-]                                                         
                                &  temp_in   , &  !< meaningfull guess of sea-ice surface temperature [K] 
                                                  !< (e.g. tskin from IFS)                                                          
                                &  tice_now  , &  !< sea-ice temperature at previous time level [K] 
                                &  hice_now  , &  !< sea-ice thickness at previous time level [m]
                                &  tsnow_now , &  !< snow temperature at previous time level [K]
                                &  hsnow_now , &  !< snow thickness at previous time level [m]
                                &  albsi_now , &  !< sea-ice albedo at previous time level [-]
                                &  tice_new  , &  !< sea-ice temperature at new time level [K] 
                                &  hice_new  , &  !< sea-ice thickness at new time level [m]
                                &  tsnow_new , &  !< snow temperature at new time level [K]
                                &  hsnow_new , &  !< snow thickness at new time level [m]
                                &  albsi_new      !< sea-ice albedo at new time level [-]


    ! set frice_thrhld depending on tile usage
    IF ( ntiles_total == 1 .AND. ALL(atm_phy_nwp_config(:)%inwp_surface /= LSS_JSBACH) ) THEN  ! no tile approach
      zfrice_thrhld = 0.5_wp
    ELSE
      zfrice_thrhld = frsi_min
    ENDIF


!$OMP PARALLEL PRIVATE(jg,nblks_c,npromz_c)
    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      nblks_c   = p_patch(jg)%nblks_c
      npromz_c  = p_patch(jg)%npromz_c


!$OMP DO PRIVATE(jb,jc,nlen,jt,js,jp,ic,zminsnow_alb,zmaxsnow_alb,zsnowalb_lu,    &
!$OMP            t_fac,ilu,frsi_in,temp_in,tice_now,hice_now,tsnow_now,hsnow_now, &
!$OMP            albsi_now,tice_new,hice_new,tsnow_new,hsnow_new,albsi_new) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_c

        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF


        ! ground temperature
        DO jc = 1, nlen
          p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_g(jc,jb) = initicon(jg)%sfc%tskin(jc,jb)
          p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%t_g(jc,jb) = initicon(jg)%sfc%tskin(jc,jb)
        ENDDO
        ! In addition, write skin temperature to lake points, limited to 33 deg C. We stick 
        ! to that until something more reasonable becomes available
        DO ic = 1, ext_data(jg)%atm%list_lake%ncount(jb)
          jc = ext_data(jg)%atm%list_lake%idx(ic,jb)
          p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_g(jc,jb) = MIN(306.15_wp,initicon(jg)%sfc%tskin(jc,jb))
          p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%t_g(jc,jb) = MIN(306.15_wp,initicon(jg)%sfc%tskin(jc,jb))
        ENDDO

        ! Fill also SST and sea ice fraction fields over ocean points; SST is limited to 30 deg C
        ! Note: missing values of the sea ice fraction, which may occur due to differing land-sea masks, 
        ! are indicated with -999.9; non-ocean points are filled with zero for both fields
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, ext_data(jg)%atm%list_sea%ncount(jb)
          jc = ext_data(jg)%atm%list_sea%idx(ic,jb)
          IF ( l_sst_in .AND. initicon(jg)%sfc%sst(jc,jb) > 270._wp  ) THEN
            p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) = initicon(jg)%sfc%sst(jc,jb)
          ELSE
            p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) = MIN(303.15_wp,initicon(jg)%sfc%tskin(jc,jb))
          ENDIF
          !
          ! In case of missing sea ice fraction values, we make use of the sea 
          ! surface temperature (tskin over ocean points). For tskin<=tf_salt, 
          ! we set the sea ice fraction to one. For tskin>tf_salt, we set it to 0.
          ! Note: tf_salt=271.45K is the salt-water freezing point
          !
          IF ( initicon(jg)%sfc%seaice(jc,jb) > -999.0_wp ) THEN
            p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) = initicon(jg)%sfc%seaice(jc,jb) 
          ELSE    ! missing value
            IF ( initicon(jg)%sfc%tskin(jc,jb) <= tf_salt ) THEN
              p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) = 1._wp     ! sea ice point
            ELSE
              p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) = 0._wp     ! water point
            ENDIF
          ENDIF

        ENDDO


        IF ( atm_phy_nwp_config(jg)%inwp_surface > 0 ) THEN
          DO jt = 1, ntiles_total
            DO jc = 1, nlen
              p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_snow_t(jc,jb,jt) = initicon(jg)%sfc%tsnow(jc,jb)

              ilu = MAX(1,ext_data(jg)%atm%lc_class_t(jc,jb,jt))
              zsnowalb_lu = ABS(ext_data(jg)%atm%snowalb_lcc(ilu))

              IF (ntiles_total > 1 .AND. albedo_type == MODIS) THEN
                zmaxsnow_alb = MIN(csalb_snow_max,zsnowalb_lu)
              ELSE
                zmaxsnow_alb = csalb_snow_max
              ENDIF

               ! Initialize freshsnow
               ! for seapoints, freshsnow is set to 0
              IF(alb_snow_var == 'ALB_SNOW') THEN

                IF (albedo_type == MODIS) THEN
                  IF (ext_data(jg)%atm%alb_dif(jc,jb) > csalb_snow_min) THEN
                    t_fac = MIN(1._wp,0.1_wp*(tmelt-p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_snow_t(jc,jb,jt)))
                    zminsnow_alb = (1._wp-t_fac)*csalb_snow_min + t_fac*ext_data(jg)%atm%alb_dif(jc,jb)
                  ELSE
                    zminsnow_alb = MAX(0.4_wp*csalb_snow_min,MIN(csalb_snow_min,0.6_wp*zsnowalb_lu))
                  ENDIF
                ELSE
                  IF (ext_data(jg)%atm%lc_class_t(jc,jb,jt) == ext_data(jg)%atm%i_lc_snow_ice) THEN
                    t_fac = MIN(1._wp,0.1_wp*(tmelt-p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_snow_t(jc,jb,jt)))
                    zminsnow_alb = (1._wp-t_fac)*csalb_snow_min + t_fac*csalb_snow
                  ELSE
                    zminsnow_alb = csalb_snow_min
                  ENDIF
                ENDIF

                p_lnd_state(jg)%diag_lnd%freshsnow_t(jc,jb,jt)    =  MAX(0._wp,MIN(1._wp, &
            &                           (initicon(jg)%sfc%snowalb (jc,jb)-zminsnow_alb)   &
            &                          /(zmaxsnow_alb-zminsnow_alb)))                     &
            &                          * REAL(NINT(ext_data(jg)%atm%fr_land(jc,jb)),wp) 
              ELSE
                p_lnd_state(jg)%diag_lnd%freshsnow_t(jc,jb,jt)    =  MAX(0._wp,MIN(1._wp, &
            &                     1._wp - ((initicon(jg)%sfc%snowalb (jc,jb)-crhosmin_ml) &
            &                    /(crhosmax_ml-crhosmin_ml))))                            &
            &                    * REAL(NINT(ext_data(jg)%atm%fr_land(jc,jb)),wp)
              END IF


              p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_snow_t(jc,jb,jt)           = &
                &                                                initicon(jg)%sfc%snowweq (jc,jb)
              p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%rho_snow_t(jc,jb,jt)         = &
                &                                                initicon(jg)%sfc%snowdens(jc,jb) 
              p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_i_t(jc,jb,jt)              = &
                &                                                initicon(jg)%sfc%skinres (jc,jb)

              p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%t_snow_t(jc,jb,jt)   = &
                &                                                initicon(jg)%sfc%tsnow   (jc,jb)
              p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%w_snow_t(jc,jb,jt)   = &
                &                                                initicon(jg)%sfc%snowweq (jc,jb)
              p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%rho_snow_t(jc,jb,jt) = &
                &                                                initicon(jg)%sfc%snowdens(jc,jb)
              p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%w_i_t(jc,jb,jt)      = &
                &                                                initicon(jg)%sfc%skinres (jc,jb)
            ENDDO
          ENDDO

          ! Multi-layer surface fields
          DO jt = 1, ntiles_total

            DO js = 0, nlev_soil
              jp = js+1 ! indexing for the ICON state field starts at 1
              DO jc = 1, nlen
                p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_so_t(jc,jp,jb,jt)= &
                  &                                              initicon(jg)%sfc%tsoil(jc,js,jb)
                p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%t_so_t(jc,jp,jb,jt)= &
                  &                                              initicon(jg)%sfc%tsoil(jc,js,jb)
              ENDDO
            ENDDO

            ! For soil water, no comparable layer shift exists
            DO js = 1, nlev_soil
              DO jc = 1, nlen
                p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%w_so_t(jc,js,jb,jt)= &
                  &                                              initicon(jg)%sfc%wsoil(jc,js,jb)
                p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%w_so_t(jc,js,jb,jt)= &
                  &                                              initicon(jg)%sfc%wsoil(jc,js,jb)
              ENDDO
            ENDDO

            ! set t_s and t_sk for land tiles to t_so_t(1)
            DO ic = 1, ext_data(jg)%atm%lp_count_t(jb,jt)
              jc = ext_data(jg)%atm%idx_lst_lp_t(ic,jb,jt)
              p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_s_t(jc,jb,jt)= &
                &                                              p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_so_t(jc,1,jb,jt)
              p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%t_s_t(jc,jb,jt)= &
                &                                              p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%t_so_t(jc,1,jb,jt)
              p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_sk_t(jc,jb,jt)= &
                &                                              p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_so_t(jc,1,jb,jt)
              p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%t_sk_t(jc,jb,jt)= &
                &                                              p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%t_so_t(jc,1,jb,jt)
            ENDDO
          ENDDO


          ! Coldstart for sea-ice parameterization scheme
          ! Sea-ice surface temperature is initialized with tskin from IFS.
          ! Since the seaice index list is not yet available at this stage, we loop over 
          ! all sea points and initialize points with fr_seaice>threshold. 
          ! The threshold is 0.5 without tiles and frsi_min with tiles.
          ! Note that exactly the same threshold values must be used as in init_sea_lists. 
          ! If not, you will see what you get.
          !@Pilar: This should still work out for you, since the non-sea-ice points are 
          !        now initialized during warmstart initialization 
          !        in mo_nwp_sfc_utils:nwp_surface_init
          !

          IF (lseaice) THEN

            DO ic = 1, ext_data(jg)%atm%list_sea%ncount(jb)
              jc = ext_data(jg)%atm%list_sea%idx(ic,jb)
              frsi_in(ic)   = p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb)             
              temp_in(ic)   = initicon(jg)%sfc%tskin(jc,jb)                        
              tice_now(ic)  = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_ice    (jc,jb)
              hice_now(ic)  = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_ice    (jc,jb)
              tsnow_now(ic) = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_snow_si(jc,jb)
              hsnow_now(ic) = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_snow_si(jc,jb)
              albsi_now(ic) = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%alb_si   (jc,jb)
              tice_new(ic)  = p_lnd_state(jg)%prog_wtr(nnew_rcf(jg))%t_ice    (jc,jb)
              hice_new(ic)  = p_lnd_state(jg)%prog_wtr(nnew_rcf(jg))%h_ice    (jc,jb)
              tsnow_new(ic) = p_lnd_state(jg)%prog_wtr(nnew_rcf(jg))%t_snow_si(jc,jb)
              hsnow_new(ic) = p_lnd_state(jg)%prog_wtr(nnew_rcf(jg))%h_snow_si(jc,jb)
              albsi_new(ic) = p_lnd_state(jg)%prog_wtr(nnew_rcf(jg))%alb_si   (jc,jb)
            ENDDO  ! ic


            CALL seaice_coldinit_nwp(                                        &
              &         nswgb        = ext_data(jg)%atm%list_sea%ncount(jb), &
              &         frice_thrhld = zfrice_thrhld,                        &
              &         frsi         = frsi_in(:),                           &
              &         temp_in      = temp_in(:),                           &
              &         tice_p       = tice_now(:),                          &
              &         hice_p       = hice_now(:),                          &
              &         tsnow_p      = tsnow_now(:),                         &
              &         hsnow_p      = hsnow_now(:),                         &
              &         albsi_p      = albsi_now(:),                         &
              &         tice_n       = tice_new(:),                          &
              &         hice_n       = hice_new(:),                          &
              &         tsnow_n      = tsnow_new(:),                         &
              &         hsnow_n      = hsnow_new(:),                         &
              &         albsi_n      = albsi_new(:)                          )

            DO ic = 1, ext_data(jg)%atm%list_sea%ncount(jb)
              jc = ext_data(jg)%atm%list_sea%idx(ic,jb)
              p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_ice    (jc,jb) = tice_now(ic)
              p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_ice    (jc,jb) = hice_now(ic)
              p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_snow_si(jc,jb) = tsnow_now(ic)
              p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_snow_si(jc,jb) = hsnow_now(ic)
              p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%alb_si   (jc,jb) = albsi_now(ic)
              p_lnd_state(jg)%prog_wtr(nnew_rcf(jg))%t_ice    (jc,jb) = tice_new(ic)
              p_lnd_state(jg)%prog_wtr(nnew_rcf(jg))%h_ice    (jc,jb) = hice_new(ic)
              p_lnd_state(jg)%prog_wtr(nnew_rcf(jg))%t_snow_si(jc,jb) = tsnow_new(ic)
              p_lnd_state(jg)%prog_wtr(nnew_rcf(jg))%h_snow_si(jc,jb) = hsnow_new(ic)
              p_lnd_state(jg)%prog_wtr(nnew_rcf(jg))%alb_si   (jc,jb) = albsi_new(ic)
            ENDDO  ! ic

          ENDIF  ! lseaice

          ! Cold-start initialization of the fresh-water lake model FLake.
          ! The procedure is the same as in "int2lm".
          ! Note that no lake ice is assumed at the cold start.

          ! Make use of sfc%ls_mask in order to identify potentially problematic points, 
          ! where depth_lk>0 (lake point in ICON) but ls_mask >0.5 (land point in IFS).
          ! At these points, tskin should not be used to initialize the water temperature.

          IF (llake) THEN
            CALL flake_coldinit(                                        &
              &     nflkgb      = ext_data(jg)%atm%list_lake%ncount(jb),&  ! in
              &     idx_lst_fp  = ext_data(jg)%atm%list_lake%idx(:,jb), &  ! in
              &     depth_lk    = ext_data(jg)%atm%depth_lk     (:,jb), &  ! in
              &     tskin       = initicon(jg)%sfc%tskin        (:,jb), &  ! in
              &     t_snow_lk_p = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_snow_lk(:,jb), &
              &     h_snow_lk_p = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_snow_lk(:,jb), &
              &     t_ice_p     = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_ice    (:,jb), &
              &     h_ice_p     = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_ice    (:,jb), &
              &     t_mnw_lk_p  = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_mnw_lk (:,jb), &
              &     t_wml_lk_p  = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_wml_lk (:,jb), & 
              &     t_bot_lk_p  = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_bot_lk (:,jb), &
              &     c_t_lk_p    = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%c_t_lk   (:,jb), &
              &     h_ml_lk_p   = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_ml_lk  (:,jb), &
              &     t_b1_lk_p   = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_b1_lk  (:,jb), &
              &     h_b1_lk_p   = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%h_b1_lk  (:,jb), &
              &     t_g_lk_p    = p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_g_t    (:,jb,isub_lake) )

            ! t_s for lake tile
            DO ic = 1, ext_data(jg)%atm%list_lake%ncount(jb)
              jc = ext_data(jg)%atm%list_lake%idx(ic,jb)
              p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_s_t(jc,jb,isub_lake) = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_wml_lk(jc,jb)
              p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%t_s_t(jc,jb,isub_lake) = p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))%t_wml_lk(jc,jb)
            ENDDO

          ELSE

            DO ic = 1, ext_data(jg)%atm%list_lake%ncount(jb)
              jc = ext_data(jg)%atm%list_lake%idx(ic,jb)
              p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))%t_s_t(jc,jb,isub_lake) = MIN(306.15_wp,initicon(jg)%sfc%tskin(jc,jb))
              p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%t_s_t(jc,jb,isub_lake) = MIN(306.15_wp,initicon(jg)%sfc%tskin(jc,jb))
            ENDDO
          ENDIF  ! llake

        ENDIF   ! inwp_surface > 0
      ENDDO
!$OMP END DO NOWAIT

    ENDDO
!$OMP END PARALLEL

    ! NOTE: Initialization of sea-water and sea-ice tiles 
    ! for t_s_t is done later in mo_nwp_sfc_utils:nwp_surface_init, 
    ! because
    ! I)  index lists for sea-ice and sea-water are 
    !     not yet available at this point.
  END SUBROUTINE copy_initicon2prog_sfc


  SUBROUTINE initVarnamesDict(dictionary)
    TYPE(t_dictionary), INTENT(INOUT) :: dictionary

    ! read the map file into dictionary data structure:
    CALL dictionary%init(.FALSE.)
    IF(ana_varnames_map_file /= ' ') THEN
      IF (my_process_is_mpi_workroot()) &
        CALL dictionary%loadfile(TRIM(ana_varnames_map_file))
      CALL dictionary%bcast(p_io, p_comm_work)
    END IF
  END SUBROUTINE initVarnamesDict


  !-------------
  !>
  !! SUBROUTINE construct_initicon
  !! Ensures that all fields have a defined VALUE.
  !!   * resets all linitialized flags
  !!   * copies topography AND coordinate surfaces
  !!   * allocates the fields we USE
  !!       * zeros OUT these fields to ensure deteministic checksums
  !!   * nullificates all other pointers
  !!
  !! This initalizes all ALLOCATED memory to avoid nondeterministic
  !! checksums when ONLY a part of a field IS READ from file due to
  !! nonfull blocks.
  SUBROUTINE construct_initicon(initicon, p_patch, topography_c, metrics)
    TYPE(t_initicon_state), INTENT(INOUT) :: initicon
    TYPE(t_patch), INTENT(IN) :: p_patch
    REAL(wp), POINTER :: topography_c(:,:)
    TYPE(t_nh_metrics), INTENT(IN) :: metrics

    ! Local variables: loop control and dimensions
    INTEGER :: nlev, nlevp1, nblks_c, nblks_e, jg

    nlev = p_patch%nlev
    nlevp1 = nlev + 1
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e
    jg      = p_patch%id

    ! basic init_icon data
    initicon%const%topography_c => topography_c
    initicon%const%z_ifc        => metrics%z_ifc
    initicon%const%z_mc         => metrics%z_mc

!WS 2017-04-12:  ORDERED was added here to work around a CCE 8.5.5 bug
!$OMP ORDERED
    CALL construct_atm_in(initicon%atm_in, initicon%const)
    CALL construct_sfc_in(initicon%sfc_in)
    CALL construct_atm(initicon%atm)
    CALL construct_atm_inc(initicon%atm_inc)
    CALL construct_sfc(initicon%sfc)
    CALL construct_sfc_inc(initicon%sfc_inc)
!$OMP END ORDERED

!-------------------------------------------------------------------------

  CONTAINS

    SUBROUTINE construct_atm_in(atm_in, const)
        TYPE(t_pi_atm_in),        INTENT(INOUT) :: atm_in
        TYPE(t_init_state_const), INTENT(INOUT) :: const

        NULLIFY(atm_in%temp, &
        &       atm_in%pres, &
        &       atm_in%u, &
        &       atm_in%v, &
        &       atm_in%w, &
        &       atm_in%vn, &
        &       atm_in%qv, &
        &       atm_in%qc, &
        &       atm_in%qi, &
        &       atm_in%qr, &
        &       atm_in%qs, &
        &       atm_in%qg, &
        &       atm_in%rho, &
        &       atm_in%theta_v, &
        &       atm_in%tke, &
        !
        &       const%z_mc_in)
        atm_in%nlev         = 0
        atm_in%linitialized = .FALSE.
    END SUBROUTINE construct_atm_in

    SUBROUTINE construct_sfc_in(sfc_in)
        TYPE(t_pi_sfc_in), INTENT(INOUT) :: sfc_in

        sfc_in%nlevsoil     = 0
        sfc_in%linitialized = .FALSE.
    END SUBROUTINE construct_sfc_in

    ! Allocate atmospheric output data
    SUBROUTINE construct_atm(atm)
        TYPE(t_pi_atm), INTENT(INOUT) :: atm

        IF(lvert_remap_fg .OR. ANY((/MODE_IFSANA, MODE_DWDANA, MODE_COSMO, MODE_COMBINED, MODE_ICONVREMAP/)==init_mode)) THEN
            ALLOCATE(atm%vn     (nproma,nlev  ,nblks_e), &
            &        atm%u      (nproma,nlev  ,nblks_c), &
            &        atm%v      (nproma,nlev  ,nblks_c), &
            &        atm%w      (nproma,nlevp1,nblks_c), &
            &        atm%temp   (nproma,nlev  ,nblks_c), &
            &        atm%exner  (nproma,nlev  ,nblks_c), &
            &        atm%pres   (nproma,nlev  ,nblks_c), &
            &        atm%rho    (nproma,nlev  ,nblks_c), &
            &        atm%theta_v(nproma,nlev  ,nblks_c), &
            &        atm%qv     (nproma,nlev  ,nblks_c), &
            &        atm%qc     (nproma,nlev  ,nblks_c), &
            &        atm%qi     (nproma,nlev  ,nblks_c), &
            &        atm%qr     (nproma,nlev  ,nblks_c), &
            &        atm%qs     (nproma,nlev  ,nblks_c)  )
!$OMP PARALLEL 
            CALL init(atm%vn(:,:,:))
            CALL init(atm%u(:,:,:))
            CALL init(atm%v(:,:,:))
            CALL init(atm%w(:,:,:))
            CALL init(atm%temp(:,:,:))
            CALL init(atm%exner(:,:,:))
            CALL init(atm%pres(:,:,:))
            CALL init(atm%rho(:,:,:))
            CALL init(atm%theta_v(:,:,:))
            CALL init(atm%qv(:,:,:))
            CALL init(atm%qc(:,:,:))
            CALL init(atm%qi(:,:,:))
            CALL init(atm%qr(:,:,:))
            CALL init(atm%qs(:,:,:))
!$OMP END PARALLEL

            IF(lvert_remap_fg .OR. init_mode == MODE_ICONVREMAP) THEN
                ALLOCATE(atm%tke(nproma,nlevp1,nblks_c))
!$OMP PARALLEL 
                CALL init(atm%tke(:,:,:))
!$OMP END PARALLEL
            END IF

            IF (atm_phy_nwp_config(jg)%lhave_graupel) THEN
                ALLOCATE(atm%qg(nproma,nlev,nblks_c))
!$OMP PARALLEL 
                CALL init(atm%qg(:,:,:))
!$OMP END PARALLEL
            END IF

            IF (atm_phy_nwp_config(jg)%l2moment) THEN
                ALLOCATE(atm%qh(nproma,nlev,nblks_c))
                ALLOCATE(atm%qnc(nproma,nlev,nblks_c))
                ALLOCATE(atm%qni(nproma,nlev,nblks_c))
                ALLOCATE(atm%qnr(nproma,nlev,nblks_c))
                ALLOCATE(atm%qns(nproma,nlev,nblks_c))
                ALLOCATE(atm%qng(nproma,nlev,nblks_c))
                ALLOCATE(atm%qnh(nproma,nlev,nblks_c))
!$OMP PARALLEL 
                CALL init(atm%qh(:,:,:))
                CALL init(atm%qnc(:,:,:))
                CALL init(atm%qni(:,:,:))
                CALL init(atm%qnr(:,:,:))
                CALL init(atm%qns(:,:,:))
                CALL init(atm%qng(:,:,:))
                CALL init(atm%qnh(:,:,:))
!$OMP END PARALLEL
            END IF

            atm%nlev         = nlev
            atm%linitialized = .TRUE.
        ELSE
            atm%nlev         = 0
            atm%linitialized = .FALSE.
        END IF
    END SUBROUTINE construct_atm

    ! atmospheric assimilation increments
    SUBROUTINE construct_atm_inc(atm_inc)
        TYPE(t_pi_atm), INTENT(INOUT) :: atm_inc

        IF ( ANY((/MODE_IAU, MODE_IAU_OLD/) == init_mode) ) THEN
            ALLOCATE(atm_inc%temp(nproma,nlev,nblks_c), &
            &        atm_inc%pres(nproma,nlev,nblks_c), &
            &        atm_inc%u   (nproma,nlev,nblks_c), &
            &        atm_inc%v   (nproma,nlev,nblks_c), &
            &        atm_inc%vn  (nproma,nlev,nblks_e), &
            &        atm_inc%qv  (nproma,nlev,nblks_c)  )
!$OMP PARALLEL 
            CALL init(atm_inc%temp(:,:,:))
            CALL init(atm_inc%pres(:,:,:))
            CALL init(atm_inc%u(:,:,:))
            CALL init(atm_inc%v(:,:,:))
            CALL init(atm_inc%vn(:,:,:))
            CALL init(atm_inc%qv(:,:,:))
!$OMP END PARALLEL 

            IF (init_mode == MODE_IAU) THEN
              IF (qcana_mode > 0) THEN
                ALLOCATE(atm_inc%qc(nproma,nlev,nblks_c))
!$OMP PARALLEL 
                CALL init(atm_inc%qc(:,:,:))
!$OMP END PARALLEL 
              ENDIF
              IF (qiana_mode > 0) THEN
                ALLOCATE(atm_inc%qi(nproma,nlev,nblks_c))
!$OMP PARALLEL 
                CALL init(atm_inc%qi(:,:,:))
!$OMP END PARALLEL 
              ENDIF
              IF (qrsgana_mode > 0) THEN
                ALLOCATE(atm_inc%qr(nproma,nlev,nblks_c))
                ALLOCATE(atm_inc%qs(nproma,nlev,nblks_c))
!$OMP PARALLEL 
                CALL init(atm_inc%qr(:,:,:))
                CALL init(atm_inc%qs(:,:,:))
!$OMP END PARALLEL 
              ENDIF
              IF (qrsgana_mode > 0 .AND. atm_phy_nwp_config(jg)%lhave_graupel) THEN
                ALLOCATE(atm_inc%qg(nproma,nlev,nblks_c))
!$OMP PARALLEL 
                CALL init(atm_inc%qg(:,:,:))
!$OMP END PARALLEL 
              END IF
              IF (atm_phy_nwp_config(jg)%l2moment) THEN
                IF (qcana_mode > 0) THEN
                  ALLOCATE(atm_inc%qnc(nproma,nlev,nblks_c))
!$OMP PARALLEL 
                  CALL init(atm_inc%qnc(:,:,:))
!$OMP END PARALLEL
                END IF
                IF (qiana_mode > 0) THEN
                  ALLOCATE(atm_inc%qni(nproma,nlev,nblks_c))
!$OMP PARALLEL 
                  CALL init(atm_inc%qni(:,:,:))
!$OMP END PARALLEL
                END IF
                IF (qrsgana_mode > 0) THEN
                  ALLOCATE(atm_inc%qh(nproma,nlev,nblks_c))
!$OMP PARALLEL 
                  CALL init(atm_inc%qh(:,:,:))
!$OMP END PARALLEL 
                  ALLOCATE(atm_inc%qnr(nproma,nlev,nblks_c))
                  ALLOCATE(atm_inc%qns(nproma,nlev,nblks_c))
                  ALLOCATE(atm_inc%qng(nproma,nlev,nblks_c))
                  ALLOCATE(atm_inc%qnh(nproma,nlev,nblks_c))
!$OMP PARALLEL 
                  CALL init(atm_inc%qnr(:,:,:))
                  CALL init(atm_inc%qns(:,:,:))
                  CALL init(atm_inc%qng(:,:,:))
                  CALL init(atm_inc%qnh(:,:,:))
!$OMP END PARALLEL 
                END IF
              END IF
            ENDIF

            atm_inc%nlev         = nlev
            atm_inc%linitialized = .TRUE.
        ELSE
            atm_inc%nlev         = 0
            atm_inc%linitialized = .FALSE.
        ENDIF
    END SUBROUTINE construct_atm_inc

    ! Allocate surface output data
    SUBROUTINE construct_sfc(sfc)
        TYPE(t_pi_sfc), INTENT(INOUT) :: sfc

        ! always allocate sst (to be on the safe side)
        ALLOCATE(sfc%sst(nproma,nblks_c))
!$OMP PARALLEL 
        CALL init(sfc%sst(:,:))
!$OMP END PARALLEL

        IF(init_mode == MODE_IFSANA) THEN
            ALLOCATE(sfc%tskin   (nproma,nblks_c            ), &
            &        sfc%tsnow   (nproma,nblks_c            ), &
            &        sfc%snowalb (nproma,nblks_c            ), &
            &        sfc%snowweq (nproma,nblks_c            ), &
            &        sfc%snowdens(nproma,nblks_c            ), &
            &        sfc%skinres (nproma,nblks_c            ), &
            &        sfc%ls_mask (nproma,nblks_c            ), &
            &        sfc%seaice  (nproma,nblks_c            ), &
            &        sfc%tsoil   (nproma,0:nlev_soil,nblks_c)  )
            IF (  nlev_soil == 0 ) THEN
               ALLOCATE(sfc%wsoil   (nproma,0:nlev_soil,nblks_c)  )
            ELSE
               ALLOCATE(sfc%wsoil   (nproma,1:nlev_soil,nblks_c)  )
            ENDIF
            
!$OMP PARALLEL 
            CALL init(sfc%tskin(:,:))
            CALL init(sfc%tsnow(:,:))
            CALL init(sfc%snowalb(:,:))
            CALL init(sfc%snowweq(:,:))
            CALL init(sfc%snowdens(:,:))
            CALL init(sfc%skinres(:,:))
            CALL init(sfc%ls_mask(:,:))
            CALL init(sfc%seaice(:,:))
            CALL init(sfc%tsoil(:,:,:))
            CALL init(sfc%wsoil(:,:,:))
!$OMP END PARALLEL 
            ! note the flipped dimensions with respect to sfc_in!

            sfc%nlevsoil     = nlev_soil
            sfc%linitialized = .TRUE.
        ELSE
            sfc%nlevsoil     = 0
            sfc%linitialized = .FALSE.
        ENDIF
    END SUBROUTINE construct_sfc

    ! surface assimilation increments
    SUBROUTINE construct_sfc_inc(sfc_inc)
        TYPE(t_sfc_inc), INTENT(INOUT) :: sfc_inc

        IF ( (init_mode == MODE_IAU) .OR. (init_mode == MODE_IAU_OLD) ) THEN
            ALLOCATE(sfc_inc%w_so (nproma,nlev_soil,nblks_c ) )
!$OMP PARALLEL 
            CALL init(sfc_inc%w_so(:,:,:))
!$OMP END PARALLEL

            ! allocate additional fields for MODE_IAU
            IF (init_mode == MODE_IAU) THEN
                ALLOCATE(sfc_inc%h_snow   (nproma,nblks_c), &
                &        sfc_inc%freshsnow(nproma,nblks_c) )
                IF (icpl_da_sfcevap == 1 .OR. icpl_da_sfcevap == 2) ALLOCATE(sfc_inc%t_2m(nproma,nblks_c))

                ! initialize with 0, since some increments are only read
                ! for specific times
!$OMP PARALLEL 
                CALL init(sfc_inc%h_snow   (:,:))
                CALL init(sfc_inc%freshsnow(:,:))
                IF (icpl_da_sfcevap == 1 .OR. icpl_da_sfcevap == 2) CALL init(sfc_inc%t_2m(:,:))
!$OMP END PARALLEL
            ENDIF  ! MODE_IAU

            sfc_inc%nlevsoil     = nlev_soil
            sfc_inc%linitialized = .TRUE.
        ELSE
            sfc_inc%nlevsoil     = 0
            sfc_inc%linitialized = .FALSE.
        ENDIF
    END SUBROUTINE construct_sfc_inc

  END SUBROUTINE construct_initicon


  !-------------
  !>
  !! SUBROUTINE allocate_extana_atm
  !! Allocates fields for reading in external analysis data
  !!
  !! This initalizes all ALLOCATED memory to avoid nondeterministic
  !! checksums when ONLY a part of a field IS READ from file due to
  !! nonfull blocks.
  SUBROUTINE allocate_extana_atm (nblks_c, nblks_e, nlev_in, atm_in, const)
    INTEGER,                  INTENT(IN)    :: nblks_c, nblks_e
    INTEGER,                  INTENT(IN)    :: nlev_in           ! number of vertical levels
                                                                 ! of external analysis data
    TYPE(t_pi_atm_in),        INTENT(INOUT) :: atm_in
    TYPE(t_init_state_const), INTENT(INOUT) :: const
    ! Local variables: loop control and dimensions
    CHARACTER(len=*), PARAMETER :: routine = modname//':allocate_extana_atm'


    IF (nlev_in == 0) THEN
      CALL finish(routine, "Number of input levels <nlev_in> not yet initialized.")
    END IF

    atm_in%nlev = nlev_in

    ! Allocate atmospheric input data
    ALLOCATE(  &
      atm_in%pres    (nproma,nlev_in,nblks_c),   &
      atm_in%temp    (nproma,nlev_in,nblks_c),   &
      atm_in%u       (nproma,nlev_in,nblks_c),   &
      atm_in%v       (nproma,nlev_in,nblks_c),   &
      atm_in%vn      (nproma,nlev_in,nblks_e),   &
      atm_in%w       (nproma,nlev_in,nblks_c),   &
      atm_in%qv      (nproma,nlev_in,nblks_c),   &
      atm_in%qc      (nproma,nlev_in,nblks_c),   &
      atm_in%qi      (nproma,nlev_in,nblks_c),   &
      atm_in%qr      (nproma,nlev_in,nblks_c),   &
      atm_in%qs      (nproma,nlev_in,nblks_c),   &
      atm_in%qg      (nproma,nlev_in,nblks_c),   &
      const%z_mc_in  (nproma,nlev_in,nblks_c) )
!$OMP PARALLEL 
    CALL init(atm_in%pres(:,:,:))
    CALL init(const%z_mc_in(:,:,:))
    CALL init(atm_in%temp(:,:,:))
    CALL init(atm_in%u(:,:,:))
    CALL init(atm_in%v(:,:,:))
    CALL init(atm_in%vn(:,:,:))
    CALL init(atm_in%w(:,:,:))
    CALL init(atm_in%qv(:,:,:))
    CALL init(atm_in%qc(:,:,:))
    CALL init(atm_in%qi(:,:,:))
    CALL init(atm_in%qr(:,:,:))
    CALL init(atm_in%qs(:,:,:))
    CALL init(atm_in%qg(:,:,:))
!$OMP END PARALLEL

    IF (init_mode == MODE_ICONVREMAP .OR. lvert_remap_fg) THEN
      ALLOCATE( &
        atm_in%rho     (nproma,nlev_in  ,nblks_c), &
        atm_in%theta_v (nproma,nlev_in  ,nblks_c), &
        atm_in%tke     (nproma,nlev_in  ,nblks_c) )
!$OMP PARALLEL
      CALL init(atm_in%rho(:,:,:))
      CALL init(atm_in%theta_v(:,:,:))
      CALL init(atm_in%tke(:,:,:))
!$OMP END PARALLEL
    ENDIF

    atm_in%linitialized = .TRUE.
  END SUBROUTINE allocate_extana_atm


  !-------------
  !>
  !! SUBROUTINE allocate_extana_sfc
  !! Allocates fields for reading in external analysis data
  !!
  SUBROUTINE allocate_extana_sfc (nblks_c, nlevsoil_in, sfc_in)
    INTEGER,                INTENT(IN)    :: nblks_c
    TYPE(t_pi_sfc_in),      INTENT(INOUT) :: sfc_in
    INTEGER,                INTENT(IN)    :: nlevsoil_in    ! number of soil levels
                                                            ! of external analysis data
    ! Local variables: loop control and dimensions
    CHARACTER(len=*), PARAMETER :: routine = modname//':allocate_extana_sfc'

    IF (nlevsoil_in == 0) THEN
      CALL finish(routine, "Number of input soil levels <nlevsoil_in> not yet initialized.")
    END IF

    sfc_in%nlevsoil = nlevsoil_in

    ! Allocate surface input data
    ! The extra soil temperature levels are not read in; they are only used to simplify vertical interpolation
    ALLOCATE(  &
      sfc_in%phi      (nproma,nblks_c                ), &
      sfc_in%tskin    (nproma,nblks_c                ), &
      sfc_in%sst      (nproma,nblks_c                ), &
      sfc_in%tsnow    (nproma,nblks_c                ), &
      sfc_in%snowalb  (nproma,nblks_c                ), &
      sfc_in%snowweq  (nproma,nblks_c                ), &
      sfc_in%snowdens (nproma,nblks_c                ), &
      sfc_in%skinres  (nproma,nblks_c                ), &
      sfc_in%ls_mask  (nproma,nblks_c                ), &
      sfc_in%seaice   (nproma,nblks_c                ), &
      sfc_in%tsoil    (nproma,nblks_c,0:nlevsoil_in+1), &
      sfc_in%wsoil    (nproma,nblks_c,0:nlevsoil_in+1)  )
!$OMP PARALLEL 
    CALL init(sfc_in%phi(:,:))
    CALL init(sfc_in%tskin(:,:))
    CALL init(sfc_in%sst(:,:))
    CALL init(sfc_in%tsnow(:,:))
    CALL init(sfc_in%snowalb(:,:))
    CALL init(sfc_in%snowweq(:,:))
    CALL init(sfc_in%snowdens(:,:))
    CALL init(sfc_in%skinres(:,:))
    CALL init(sfc_in%ls_mask(:,:))
    CALL init(sfc_in%seaice(:,:))
    CALL init(sfc_in%tsoil(:,:,:))
    CALL init(sfc_in%wsoil(:,:,:))
!$OMP END PARALLEL

    sfc_in%linitialized = .TRUE.
  END SUBROUTINE allocate_extana_sfc


  !-------------
  !>
  !! SUBROUTINE deallocate_initicon
  !! Deallocates the components of the initicon data type
  !!
  SUBROUTINE deallocate_initicon (initicon)

    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)

    ! Local variables: loop control
    INTEGER :: jg
    TYPE(t_init_state_const), POINTER :: const

!-------------------------------------------------------------------------

    ! Loop over model domains
    DO jg = 1, n_dom

      const => initicon(jg)%const

      ! call destructor
      CALL initicon(jg)%finalize()
      CALL const%finalize()

    ENDDO ! loop over model domains

    ! destroy variable name dictionaries:
    CALL ana_varnames_dict%finalize()

  END SUBROUTINE deallocate_initicon


  !> output checksums of all possible input fields
  !!
  !! XXX: This FUNCTION should have been written using a few
  !! preprocessor macros taking little more than the NAME of the
  !! respective variable.
  !!
  !!      Alas, such macros would have generated code violating the
  !!      fortran line limit, so we are stuck with the expanded
  !!      version.
  SUBROUTINE printChecksums(initicon, p_nh_state, p_lnd_state)
    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)
    TYPE(t_nh_state), INTENT(INOUT) :: p_nh_state(:)
    TYPE(t_lnd_state), INTENT(INOUT), OPTIONAL :: p_lnd_state(:)

    INTEGER :: jg, i, pfx_tlen
    CHARACTER(LEN = 256) :: prefix

    IF(msg_level < 13) RETURN

    DO jg = 1, n_dom
      WRITE (prefix, '(a,i0,a)') "checksum of initicon(", jg, ")%"
      pfx_tlen = LEN_TRIM(prefix)
      IF(ASSOCIATED(initicon(jg)%atm_in%temp)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_in%temp: ", &
        & initicon(jg)%atm_in%temp)
      IF(ASSOCIATED(initicon(jg)%atm_in%pres)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_in%pres: ", &
        & initicon(jg)%atm_in%pres)
      IF(ASSOCIATED(initicon(jg)%const%z_mc_in)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"const%z_mc_in: ", &
        & initicon(jg)%const%z_mc_in)
      IF(ASSOCIATED(initicon(jg)%atm_in%u)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_in%u: ", &
        & initicon(jg)%atm_in%u)
      IF(ASSOCIATED(initicon(jg)%atm_in%v)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_in%v: ", &
        & initicon(jg)%atm_in%v)
      IF(ASSOCIATED(initicon(jg)%atm_in%w)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_in%w: ", &
        & initicon(jg)%atm_in%w)
      IF(ASSOCIATED(initicon(jg)%atm_in%vn)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_in%vn: ", &
        & initicon(jg)%atm_in%vn)
      IF(ASSOCIATED(initicon(jg)%atm_in%qv)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_in%qv: ", &
        & initicon(jg)%atm_in%qv)
      IF(ASSOCIATED(initicon(jg)%atm_in%qc)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_in%qc: ", &
        & initicon(jg)%atm_in%qc)
      IF(ASSOCIATED(initicon(jg)%atm_in%qi)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_in%qi: ", &
        & initicon(jg)%atm_in%qi)
      IF(ASSOCIATED(initicon(jg)%atm_in%qr)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_in%qr: ", &
        & initicon(jg)%atm_in%qr)
      IF(ASSOCIATED(initicon(jg)%atm_in%qs)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_in%qs: ", &
        & initicon(jg)%atm_in%qs)
      IF(ASSOCIATED(initicon(jg)%atm_in%qg)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_in%qg: ", &
        & initicon(jg)%atm_in%qg)
      IF(ASSOCIATED(initicon(jg)%atm_in%qh)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_in%qh: ", &
        & initicon(jg)%atm_in%qh)
      IF(ASSOCIATED(initicon(jg)%atm_in%qnc)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_in%qnc: ", &
        & initicon(jg)%atm_in%qnc)
      IF(ASSOCIATED(initicon(jg)%atm_in%qni)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_in%qni: ", &
        & initicon(jg)%atm_in%qni)
      IF(ASSOCIATED(initicon(jg)%atm_in%qnr)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_in%qnr: ", &
        & initicon(jg)%atm_in%qnr)
      IF(ASSOCIATED(initicon(jg)%atm_in%qns)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_in%qns: ", &
        & initicon(jg)%atm_in%qns)
      IF(ASSOCIATED(initicon(jg)%atm_in%qng)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_in%qng: ", &
        & initicon(jg)%atm_in%qng)
      IF(ASSOCIATED(initicon(jg)%atm_in%qnh)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_in%qnh: ", &
        & initicon(jg)%atm_in%qnh)
      IF(ASSOCIATED(initicon(jg)%atm_in%rho)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_in%rho: ", &
        & initicon(jg)%atm_in%rho)
      IF(ASSOCIATED(initicon(jg)%atm_in%theta_v)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_in%theta_v: ", &
        & initicon(jg)%atm_in%theta_v)
      IF(ASSOCIATED(initicon(jg)%atm_in%tke)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_in%tke: ", &
        & initicon(jg)%atm_in%tke)
      IF(ALLOCATED(initicon(jg)%sfc_in%tsnow)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc_in%tsnow: ", &
        & initicon(jg)%sfc_in%tsnow)
      IF(ALLOCATED(initicon(jg)%sfc_in%tskin)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc_in%tskin: ", &
        & initicon(jg)%sfc_in%tskin)
      IF(ALLOCATED(initicon(jg)%sfc_in%sst)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc_in%sst: ", &
        & initicon(jg)%sfc_in%sst)
      IF(ALLOCATED(initicon(jg)%sfc_in%snowalb)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc_in%snowalb: ", &
        & initicon(jg)%sfc_in%snowalb)
      IF(ALLOCATED(initicon(jg)%sfc_in%snowweq)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc_in%snowweq: ", &
        & initicon(jg)%sfc_in%snowweq)
      IF(ALLOCATED(initicon(jg)%sfc_in%snowdens)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc_in%snowdens: ", &
        & initicon(jg)%sfc_in%snowdens)
      IF(ALLOCATED(initicon(jg)%sfc_in%skinres)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc_in%skinres: ", &
        & initicon(jg)%sfc_in%skinres)
      IF(ALLOCATED(initicon(jg)%sfc_in%ls_mask)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc_in%ls_mask: ", &
        & initicon(jg)%sfc_in%ls_mask)
      IF(ALLOCATED(initicon(jg)%sfc_in%seaice)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc_in%seaice: ", &
        & initicon(jg)%sfc_in%seaice)
      IF(ALLOCATED(initicon(jg)%sfc_in%phi)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc_in%phi: ", &
        & initicon(jg)%sfc_in%phi)
      IF(ALLOCATED(initicon(jg)%sfc_in%tsoil)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc_in%tsoil: ", &
        & initicon(jg)%sfc_in%tsoil)
      IF(ALLOCATED(initicon(jg)%sfc_in%wsoil)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc_in%wsoil: ", &
        & initicon(jg)%sfc_in%wsoil)
      IF(ALLOCATED(initicon(jg)%atm%vn)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%vn: ", &
        & initicon(jg)%atm%vn)
      IF(ALLOCATED(initicon(jg)%atm%u)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%u: ", &
        & initicon(jg)%atm%u)
      IF(ALLOCATED(initicon(jg)%atm%v)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%v: ", &
        & initicon(jg)%atm%v)
      IF(ALLOCATED(initicon(jg)%atm%w)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%w: ", &
        & initicon(jg)%atm%w)
      IF(ALLOCATED(initicon(jg)%atm%temp)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%temp: ", &
        & initicon(jg)%atm%temp)
      IF(ALLOCATED(initicon(jg)%atm%theta_v)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%theta_v: ", &
        & initicon(jg)%atm%theta_v)
      IF(ALLOCATED(initicon(jg)%atm%exner)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%exner: ", &
        & initicon(jg)%atm%exner)
      IF(ALLOCATED(initicon(jg)%atm%rho)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%rho: ", &
        & initicon(jg)%atm%rho)
      IF(ALLOCATED(initicon(jg)%atm%pres)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%pres: ", &
        & initicon(jg)%atm%pres)
      IF(ALLOCATED(initicon(jg)%atm%qv)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%qv: ", &
        & initicon(jg)%atm%qv)
      IF(ALLOCATED(initicon(jg)%atm%qc)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%qc: ", &
        & initicon(jg)%atm%qc)
      IF(ALLOCATED(initicon(jg)%atm%qi)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%qi: ", &
        & initicon(jg)%atm%qi)
      IF(ALLOCATED(initicon(jg)%atm%qr)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%qr: ", &
        & initicon(jg)%atm%qr)
      IF(ALLOCATED(initicon(jg)%atm%qs)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%qs: ", &
        & initicon(jg)%atm%qs)
      IF(ALLOCATED(initicon(jg)%atm%qg)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%qg: ", &
        & initicon(jg)%atm%qg)
      IF(ALLOCATED(initicon(jg)%atm%qh)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%qh: ", &
        & initicon(jg)%atm%qh)
      IF(ALLOCATED(initicon(jg)%atm%qnc)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%qnc: ", &
        & initicon(jg)%atm%qnc)
      IF(ALLOCATED(initicon(jg)%atm%qni)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%qni: ", &
        & initicon(jg)%atm%qni)
      IF(ALLOCATED(initicon(jg)%atm%qnr)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%qnr: ", &
        & initicon(jg)%atm%qnr)
      IF(ALLOCATED(initicon(jg)%atm%qns)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%qns: ", &
        & initicon(jg)%atm%qns)
      IF(ALLOCATED(initicon(jg)%atm%qng)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%qng: ", &
        & initicon(jg)%atm%qng)
      IF(ALLOCATED(initicon(jg)%atm%qnh)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%qnh: ", &
        & initicon(jg)%atm%qnh)
      IF(ALLOCATED(initicon(jg)%atm%tke)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm%tke: ", &
        & initicon(jg)%atm%tke)
      IF(ALLOCATED(initicon(jg)%atm_inc%vn)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%vn: ", &
        & initicon(jg)%atm_inc%vn)
      IF(ALLOCATED(initicon(jg)%atm_inc%u)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%u: ", &
        & initicon(jg)%atm_inc%u)
      IF(ALLOCATED(initicon(jg)%atm_inc%v)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%v: ", &
        & initicon(jg)%atm_inc%v)
      IF(ALLOCATED(initicon(jg)%atm_inc%w)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%w: ", &
        & initicon(jg)%atm_inc%w)
      IF(ALLOCATED(initicon(jg)%atm_inc%temp)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%temp: ", &
        & initicon(jg)%atm_inc%temp)
      IF(ALLOCATED(initicon(jg)%atm_inc%theta_v)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%theta_v: ", &
        & initicon(jg)%atm_inc%theta_v)
      IF(ALLOCATED(initicon(jg)%atm_inc%exner)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%exner: ", &
        & initicon(jg)%atm_inc%exner)
      IF(ALLOCATED(initicon(jg)%atm_inc%rho)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%rho: ", &
        & initicon(jg)%atm_inc%rho)
      IF(ALLOCATED(initicon(jg)%atm_inc%pres)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%pres: ", &
        & initicon(jg)%atm_inc%pres)
      IF(ALLOCATED(initicon(jg)%atm_inc%qv)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%qv: ", &
        & initicon(jg)%atm_inc%qv)
      IF(ALLOCATED(initicon(jg)%atm_inc%qc)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%qc: ", &
        & initicon(jg)%atm_inc%qc)
      IF(ALLOCATED(initicon(jg)%atm_inc%qi)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%qi: ", &
        & initicon(jg)%atm_inc%qi)
      IF(ALLOCATED(initicon(jg)%atm_inc%qr)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%qr: ", &
        & initicon(jg)%atm_inc%qr)
      IF(ALLOCATED(initicon(jg)%atm_inc%qs)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%qs: ", &
        & initicon(jg)%atm_inc%qs)
      IF(ALLOCATED(initicon(jg)%atm_inc%qg)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%qg: ", &
        & initicon(jg)%atm_inc%qg)
      IF(ALLOCATED(initicon(jg)%atm_inc%qh)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%qh: ", &
        & initicon(jg)%atm_inc%qh)
      IF(ALLOCATED(initicon(jg)%atm_inc%qnc)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%qnc: ", &
        & initicon(jg)%atm_inc%qnc)
      IF(ALLOCATED(initicon(jg)%atm_inc%qni)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%qni: ", &
        & initicon(jg)%atm_inc%qni)
      IF(ALLOCATED(initicon(jg)%atm_inc%qnr)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%qnr: ", &
        & initicon(jg)%atm_inc%qnr)
      IF(ALLOCATED(initicon(jg)%atm_inc%qns)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%qns: ", &
        & initicon(jg)%atm_inc%qns)
      IF(ALLOCATED(initicon(jg)%atm_inc%qng)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%qng: ", &
        & initicon(jg)%atm_inc%qng)
      IF(ALLOCATED(initicon(jg)%atm_inc%qnh)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%qnh: ", &
        & initicon(jg)%atm_inc%qnh)
      IF(ALLOCATED(initicon(jg)%atm_inc%tke)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"atm_inc%tke: ", &
        & initicon(jg)%atm_inc%tke)
      IF(ALLOCATED(initicon(jg)%sfc%tsnow)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc%tsnow: ", &
        & initicon(jg)%sfc%tsnow)
      IF(ALLOCATED(initicon(jg)%sfc%tskin)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc%tskin: ", &
        & initicon(jg)%sfc%tskin)
      IF(ALLOCATED(initicon(jg)%sfc%sst)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc%sst: ", &
        & initicon(jg)%sfc%sst)
      IF(ALLOCATED(initicon(jg)%sfc%snowalb)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc%snowalb: ", &
        & initicon(jg)%sfc%snowalb)
      IF(ALLOCATED(initicon(jg)%sfc%snowweq)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc%snowweq: ", &
        & initicon(jg)%sfc%snowweq)
      IF(ALLOCATED(initicon(jg)%sfc%snowdens)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc%snowdens: ", &
        & initicon(jg)%sfc%snowdens)
      IF(ALLOCATED(initicon(jg)%sfc%skinres)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc%skinres: ", &
        & initicon(jg)%sfc%skinres)
      IF(ALLOCATED(initicon(jg)%sfc%ls_mask)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc%ls_mask: ", &
        & initicon(jg)%sfc%ls_mask)
      IF(ALLOCATED(initicon(jg)%sfc%seaice)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc%seaice: ", &
        & initicon(jg)%sfc%seaice)
      IF(ALLOCATED(initicon(jg)%sfc%tsoil)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc%tsoil: ", &
        & initicon(jg)%sfc%tsoil)
      IF(ALLOCATED(initicon(jg)%sfc%wsoil)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc%wsoil: ", &
        & initicon(jg)%sfc%wsoil)
      IF(ALLOCATED(initicon(jg)%sfc%w_so)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc%w_so: ", &
        & initicon(jg)%sfc%w_so)
      IF(ALLOCATED(initicon(jg)%sfc_inc%w_so)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc_inc%w_so: ", &
        & initicon(jg)%sfc_inc%w_so)
      IF(ALLOCATED(initicon(jg)%sfc_inc%h_snow)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc_inc%h_snow: ", &
        & initicon(jg)%sfc_inc%h_snow)
      IF(ALLOCATED(initicon(jg)%sfc_inc%freshsnow)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"sfc_inc%freshsnow: ", &
        & initicon(jg)%sfc_inc%freshsnow)

      IF(ALLOCATED(p_nh_state(jg)%prog)) THEN
        DO i = 1, SIZE(p_nh_state(jg)%prog, 1)
          WRITE (prefix, '(a,2(i0,a))') &
            "checksum of p_nh_state(", jg, ")%prog(", i, ")%"
          pfx_tlen = LEN_TRIM(prefix)
          IF(ASSOCIATED(p_nh_state(jg)%prog(i)%w)) &
            & CALL printChecksum(prefix(1:pfx_tlen)//"w: ", &
            & p_nh_state(jg)%prog(i)%w)
          IF(ASSOCIATED(p_nh_state(jg)%prog(i)%vn)) &
            & CALL printChecksum(prefix(1:pfx_tlen)//"vn: ", &
            & p_nh_state(jg)%prog(i)%vn)
          IF(ASSOCIATED(p_nh_state(jg)%prog(i)%rho)) &
            & CALL printChecksum(prefix(1:pfx_tlen)//"rho: ", &
            & p_nh_state(jg)%prog(i)%rho)
          IF(ASSOCIATED(p_nh_state(jg)%prog(i)%exner)) &
            & CALL printChecksum(prefix(1:pfx_tlen)//"exner: ", &
            & p_nh_state(jg)%prog(i)%exner)
          IF(ASSOCIATED(p_nh_state(jg)%prog(i)%theta_v)) &
            & CALL printChecksum(prefix(1:pfx_tlen)//"theta_v: ", &
            & p_nh_state(jg)%prog(i)%theta_v)
          IF(ASSOCIATED(p_nh_state(jg)%prog(i)%tracer)) &
            & CALL printChecksum(prefix(1:pfx_tlen)//"tracer: ", &
            & p_nh_state(jg)%prog(i)%tracer)
          IF(ASSOCIATED(p_nh_state(jg)%prog(i)%tke)) &
            & CALL printChecksum(prefix(1:pfx_tlen)//"tke: ", &
            & p_nh_state(jg)%prog(i)%tke)
        END DO
      END IF
      WRITE (prefix, '(a,i0,a)') "checksum of p_nh_state(", jg, ")%diag"
      pfx_tlen = LEN_TRIM(prefix)
      IF(ASSOCIATED(p_nh_state(jg)%diag%u)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"u: ", &
        & p_nh_state(jg)%diag%u)
      IF(ASSOCIATED(p_nh_state(jg)%diag%v)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"v: ", &
        & p_nh_state(jg)%diag%v)
      IF(ASSOCIATED(p_nh_state(jg)%diag%vt)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"vt: ", &
        & p_nh_state(jg)%diag%vt)
      IF(ASSOCIATED(p_nh_state(jg)%diag%omega_z)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"omega_z: ", &
        & p_nh_state(jg)%diag%omega_z)
      IF(ASSOCIATED(p_nh_state(jg)%diag%vor)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"vor: ", &
        & p_nh_state(jg)%diag%vor)
      IF(ASSOCIATED(p_nh_state(jg)%diag%ddt_vn_phy)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"ddt_vn_phy: ", &
        & p_nh_state(jg)%diag%ddt_vn_phy)
      IF(ASSOCIATED(p_nh_state(jg)%diag%ddt_exner_phy)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"ddt_exner_phy: ", &
        & p_nh_state(jg)%diag%ddt_exner_phy)
      IF(ASSOCIATED(p_nh_state(jg)%diag%ddt_temp_dyn)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"ddt_temp_dyn: ", &
        & p_nh_state(jg)%diag%ddt_temp_dyn)
      IF(ASSOCIATED(p_nh_state(jg)%diag%ddt_tracer_adv)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"ddt_tracer_adv: ", &
        & p_nh_state(jg)%diag%ddt_tracer_adv)
      IF(ASSOCIATED(p_nh_state(jg)%diag%tracer_vi)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"tracer_vi: ", &
        & p_nh_state(jg)%diag%tracer_vi)
      IF(ASSOCIATED(p_nh_state(jg)%diag%exner_pr)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"exner_pr: ", &
        & p_nh_state(jg)%diag%exner_pr)
      IF(ASSOCIATED(p_nh_state(jg)%diag%exner_dyn_incr)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"exner_dyn_incr: ", &
        & p_nh_state(jg)%diag%exner_dyn_incr)
      IF(ASSOCIATED(p_nh_state(jg)%diag%temp)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"temp: ", &
        & p_nh_state(jg)%diag%temp)
      IF(ASSOCIATED(p_nh_state(jg)%diag%tempv)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"tempv: ", &
        & p_nh_state(jg)%diag%tempv)
      IF(ASSOCIATED(p_nh_state(jg)%diag%temp_ifc)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"temp_ifc: ", &
        & p_nh_state(jg)%diag%temp_ifc)
      IF(ASSOCIATED(p_nh_state(jg)%diag%pres)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"pres: ", &
        & p_nh_state(jg)%diag%pres)
      IF(ASSOCIATED(p_nh_state(jg)%diag%pres_ifc)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"pres_ifc: ", &
        & p_nh_state(jg)%diag%pres_ifc)
      IF(ASSOCIATED(p_nh_state(jg)%diag%pres_sfc)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"pres_sfc: ", &
        & p_nh_state(jg)%diag%pres_sfc)
      IF(ASSOCIATED(p_nh_state(jg)%diag%pres_sfc_old)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"pres_sfc_old: ", &
        & p_nh_state(jg)%diag%pres_sfc_old)
      IF(ASSOCIATED(p_nh_state(jg)%diag%pres_msl)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"pres_msl: ", &
        & p_nh_state(jg)%diag%pres_msl)
      IF(ASSOCIATED(p_nh_state(jg)%diag%dpres_mc)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"dpres_mc: ", &
        & p_nh_state(jg)%diag%dpres_mc)
      IF(ASSOCIATED(p_nh_state(jg)%diag%hfl_tracer)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"hfl_tracer: ", &
        & p_nh_state(jg)%diag%hfl_tracer)
      IF(ASSOCIATED(p_nh_state(jg)%diag%vfl_tracer)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"vfl_tracer: ", &
        & p_nh_state(jg)%diag%vfl_tracer)
      IF(ASSOCIATED(p_nh_state(jg)%diag%div)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"div: ", &
        & p_nh_state(jg)%diag%div)
      ! For some reason, these provide nondeterministic results.
      ! IF(ASSOCIATED(p_nh_state(jg)%diag%div_ic)) &
      !   & CALL printChecksum(prefix(1:pfx_tlen)//"div_ic: ", &
      !   & p_nh_state(jg)%diag%div_ic)
      ! IF(ASSOCIATED(p_nh_state(jg)%diag%hdef_ic)) &
      !   & CALL printChecksum(prefix(1:pfx_tlen)//"hdef_ic: ", &
      !   & p_nh_state(jg)%diag%hdef_ic)
      ! IF(ASSOCIATED(p_nh_state(jg)%diag%dwdx)) &
      !   & CALL printChecksum(prefix(1:pfx_tlen)//"dwdx: ", &
      !   & p_nh_state(jg)%diag%dwdx)
      ! IF(ASSOCIATED(p_nh_state(jg)%diag%dwdy)) &
      !   & CALL printChecksum(prefix(1:pfx_tlen)//"dwdy: ", &
      !   & p_nh_state(jg)%diag%dwdy)
      ! IF(ASSOCIATED(p_nh_state(jg)%diag%mass_fl_e)) &
      !   & CALL printChecksum(prefix(1:pfx_tlen)//"mass_fl_e: ", &
      !   & p_nh_state(jg)%diag%mass_fl_e)
      IF(ASSOCIATED(p_nh_state(jg)%diag%mass_fl_e_sv)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"mass_fl_e_sv: ", &
        & p_nh_state(jg)%diag%mass_fl_e_sv)
      IF(ASSOCIATED(p_nh_state(jg)%diag%rho_ic)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"rho_ic: ", &
        & p_nh_state(jg)%diag%rho_ic)
      IF(ASSOCIATED(p_nh_state(jg)%diag%theta_v_ic)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"theta_v_ic: ", &
        & p_nh_state(jg)%diag%theta_v_ic)
      IF(ASSOCIATED(p_nh_state(jg)%diag%w_concorr_c)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"w_concorr_c: ", &
        & p_nh_state(jg)%diag%w_concorr_c)
      IF(ASSOCIATED(p_nh_state(jg)%diag%vn_ie)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"vn_ie: ", &
        & p_nh_state(jg)%diag%vn_ie)
      IF(ASSOCIATED(p_nh_state(jg)%diag%ddt_vn_apc_pc)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"ddt_vn_apc_pc: ", &
        & p_nh_state(jg)%diag%ddt_vn_apc_pc)
      IF(ASSOCIATED(p_nh_state(jg)%diag%ddt_vn_cor_pc)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"ddt_vn_cor_pc: ", &
        & p_nh_state(jg)%diag%ddt_vn_cor_pc)
      IF(ASSOCIATED(p_nh_state(jg)%diag%ddt_w_adv_pc)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"ddt_w_adv_pc: ", &
        & p_nh_state(jg)%diag%ddt_w_adv_pc)
      IF(ASSOCIATED(p_nh_state(jg)%diag%airmass_now)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"airmass_now: ", &
        & p_nh_state(jg)%diag%airmass_now)
      IF(ASSOCIATED(p_nh_state(jg)%diag%airmass_new)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"airmass_new: ", &
        & p_nh_state(jg)%diag%airmass_new)
      IF(ASSOCIATED(p_nh_state(jg)%diag%grf_tend_vn)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"grf_tend_vn: ", &
        & p_nh_state(jg)%diag%grf_tend_vn)
      IF(ASSOCIATED(p_nh_state(jg)%diag%grf_tend_w)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"grf_tend_w: ", &
        & p_nh_state(jg)%diag%grf_tend_w)
      IF(ASSOCIATED(p_nh_state(jg)%diag%grf_tend_rho)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"grf_tend_rho: ", &
        & p_nh_state(jg)%diag%grf_tend_rho)
      IF(ASSOCIATED(p_nh_state(jg)%diag%grf_tend_mflx)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"grf_tend_mflx: ", &
        & p_nh_state(jg)%diag%grf_tend_mflx)
      IF(ASSOCIATED(p_nh_state(jg)%diag%grf_bdy_mflx)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"grf_bdy_mflx: ", &
        & p_nh_state(jg)%diag%grf_bdy_mflx)
      IF(ASSOCIATED(p_nh_state(jg)%diag%grf_tend_thv)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"grf_tend_thv: ", &
        & p_nh_state(jg)%diag%grf_tend_thv)
      IF(ASSOCIATED(p_nh_state(jg)%diag%grf_tend_tracer)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"grf_tend_tracer: ", &
        & p_nh_state(jg)%diag%grf_tend_tracer)
      IF(ASSOCIATED(p_nh_state(jg)%diag%vn_ie_int)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"vn_ie_int: ", &
        & p_nh_state(jg)%diag%vn_ie_int)
      IF(ASSOCIATED(p_nh_state(jg)%diag%vn_ie_ubc)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"vn_ie_ubc: ", &
        & p_nh_state(jg)%diag%vn_ie_ubc)
      IF(ASSOCIATED(p_nh_state(jg)%diag%mflx_ic_int)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"mflx_ic_int: ", &
        & p_nh_state(jg)%diag%mflx_ic_int)
      IF(ASSOCIATED(p_nh_state(jg)%diag%mflx_ic_ubc)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"mflx_ic_ubc: ", &
        & p_nh_state(jg)%diag%mflx_ic_ubc)
      IF(ASSOCIATED(p_nh_state(jg)%diag%vn_incr)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"vn_incr: ", &
        & p_nh_state(jg)%diag%vn_incr)
      IF(ASSOCIATED(p_nh_state(jg)%diag%exner_incr)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"exner_incr: ", &
        & p_nh_state(jg)%diag%exner_incr)
      IF(ASSOCIATED(p_nh_state(jg)%diag%rho_incr)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"rho_incr: ", &
        & p_nh_state(jg)%diag%rho_incr)
      IF(ASSOCIATED(p_nh_state(jg)%diag%rhov_incr)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"rhov_incr: ", &
        & p_nh_state(jg)%diag%rhov_incr)
      WRITE (prefix, '(a,i0,a)') "checksum of p_nh_state(", jg, ")%ref%"
      pfx_tlen = LEN_TRIM(prefix)
      IF(ASSOCIATED(p_nh_state(jg)%ref%vn_ref)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"vn_ref: ", &
        & p_nh_state(jg)%ref%vn_ref)
      IF(ASSOCIATED(p_nh_state(jg)%ref%w_ref)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"w_ref: ", &
        & p_nh_state(jg)%ref%w_ref)
      WRITE (prefix, '(a,i0,a)') "checksum of p_nh_state(", jg, ")%metrics%"
      pfx_tlen = LEN_TRIM(prefix)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%ddqz_z_full)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"ddqz_z_full: ", &
        & p_nh_state(jg)%metrics%ddqz_z_full)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%geopot)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"geopot: ", &
        & p_nh_state(jg)%metrics%geopot)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%geopot_agl)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"geopot_agl: ", &
        & p_nh_state(jg)%metrics%geopot_agl)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%geopot_agl_ifc)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"geopot_agl_ifc: ", &
        & p_nh_state(jg)%metrics%geopot_agl_ifc)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%dgeopot_mc)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"dgeopot_mc: ", &
        & p_nh_state(jg)%metrics%dgeopot_mc)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%rayleigh_w)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"rayleigh_w: ", &
        & p_nh_state(jg)%metrics%rayleigh_w)
      ! For some reason, this provides nondeterministic results.
      !       IF(ASSOCIATED(p_nh_state(jg)%metrics%rayleigh_vn)) &
      !         & CALL printChecksum(prefix(1:pfx_tlen)//"rayleigh_vn: ", &
      !         & p_nh_state(jg)%metrics%rayleigh_vn)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%enhfac_diffu)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"enhfac_diffu: ", &
        & p_nh_state(jg)%metrics%enhfac_diffu)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%scalfac_dd3d)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"scalfac_dd3d: ", &
        & p_nh_state(jg)%metrics%scalfac_dd3d)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%vwind_expl_wgt)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"vwind_expl_wgt: ", &
        & p_nh_state(jg)%metrics%vwind_expl_wgt)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%vwind_impl_wgt)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"vwind_impl_wgt: ", &
        & p_nh_state(jg)%metrics%vwind_impl_wgt)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%theta_ref_mc)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"theta_ref_mc: ", &
        & p_nh_state(jg)%metrics%theta_ref_mc)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%theta_ref_me)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"theta_ref_me: ", &
        & p_nh_state(jg)%metrics%theta_ref_me)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%theta_ref_ic)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"theta_ref_ic: ", &
        & p_nh_state(jg)%metrics%theta_ref_ic)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%tsfc_ref)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"tsfc_ref: ", &
        & p_nh_state(jg)%metrics%tsfc_ref)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%exner_ref_mc)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"exner_ref_mc: ", &
        & p_nh_state(jg)%metrics%exner_ref_mc)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%rho_ref_mc)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"rho_ref_mc: ", &
        & p_nh_state(jg)%metrics%rho_ref_mc)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%rho_ref_me)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"rho_ref_me: ", &
        & p_nh_state(jg)%metrics%rho_ref_me)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%zd_intcoef)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"zd_intcoef: ", &
        & p_nh_state(jg)%metrics%zd_intcoef)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%zd_geofac)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"zd_geofac: ", &
        & p_nh_state(jg)%metrics%zd_geofac)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%zd_e2cell)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"zd_e2cell: ", &
        & p_nh_state(jg)%metrics%zd_e2cell)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%zd_diffcoef)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"zd_diffcoef: ", &
        & p_nh_state(jg)%metrics%zd_diffcoef)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%inv_ddqz_z_half_e)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"inv_ddqz_z_half_e: ", &
        & p_nh_state(jg)%metrics%inv_ddqz_z_half_e)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%inv_ddqz_z_full_e)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"inv_ddqz_z_full_e: ", &
        & p_nh_state(jg)%metrics%inv_ddqz_z_full_e)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%inv_ddqz_z_half)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"inv_ddqz_z_half: ", &
        & p_nh_state(jg)%metrics%inv_ddqz_z_half)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%inv_ddqz_z_half_v)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"inv_ddqz_z_half_v: ", &
        & p_nh_state(jg)%metrics%inv_ddqz_z_half_v)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%wgtfac_v)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"wgtfac_v: ", &
        & p_nh_state(jg)%metrics%wgtfac_v)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%mixing_length_sq)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"mixing_length_sq: ", &
        & p_nh_state(jg)%metrics%mixing_length_sq)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%rho_ref_corr)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"rho_ref_corr: ", &
        & p_nh_state(jg)%metrics%rho_ref_corr)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%ddxn_z_full)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"ddxn_z_full: ", &
        & p_nh_state(jg)%metrics%ddxn_z_full)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%ddxt_z_full)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"ddxt_z_full: ", &
        & p_nh_state(jg)%metrics%ddxt_z_full)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%ddqz_z_full_e)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"ddqz_z_full_e: ", &
        & p_nh_state(jg)%metrics%ddqz_z_full_e)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%ddqz_z_half)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"ddqz_z_half: ", &
        & p_nh_state(jg)%metrics%ddqz_z_half)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%inv_ddqz_z_full)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"inv_ddqz_z_full: ", &
        & p_nh_state(jg)%metrics%inv_ddqz_z_full)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%wgtfac_c)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"wgtfac_c: ", &
        & p_nh_state(jg)%metrics%wgtfac_c)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%wgtfac_e)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"wgtfac_e: ", &
        & p_nh_state(jg)%metrics%wgtfac_e)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%wgtfacq_c)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"wgtfacq_c: ", &
        & p_nh_state(jg)%metrics%wgtfacq_c)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%wgtfacq_e)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"wgtfacq_e: ", &
        & p_nh_state(jg)%metrics%wgtfacq_e)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%wgtfacq1_c)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"wgtfacq1_c: ", &
        & p_nh_state(jg)%metrics%wgtfacq1_c)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%wgtfacq1_e)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"wgtfacq1_e: ", &
        & p_nh_state(jg)%metrics%wgtfacq1_e)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%coeff_gradekin)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"coeff_gradekin: ", &
        & p_nh_state(jg)%metrics%coeff_gradekin)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%coeff1_dwdz)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"coeff1_dwdz: ", &
        & p_nh_state(jg)%metrics%coeff1_dwdz)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%coeff2_dwdz)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"coeff2_dwdz: ", &
        & p_nh_state(jg)%metrics%coeff2_dwdz)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%zdiff_gradp)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"zdiff_gradp: ", &
        & p_nh_state(jg)%metrics%zdiff_gradp)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%coeff_gradp)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"coeff_gradp: ", &
        & p_nh_state(jg)%metrics%coeff_gradp)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%exner_exfac)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"exner_exfac: ", &
        & p_nh_state(jg)%metrics%exner_exfac)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%d_exner_dz_ref_ic)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"d_exner_dz_ref_ic: ", &
        & p_nh_state(jg)%metrics%d_exner_dz_ref_ic)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%d2dexdz2_fac1_mc)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"d2dexdz2_fac1_mc: ", &
        & p_nh_state(jg)%metrics%d2dexdz2_fac1_mc)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%d2dexdz2_fac2_mc)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"d2dexdz2_fac2_mc: ", &
        & p_nh_state(jg)%metrics%d2dexdz2_fac2_mc)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%pg_exdist)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"pg_exdist: ", &
        & p_nh_state(jg)%metrics%pg_exdist)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%vertidx_gradp)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"vertidx_gradp: ", &
        & p_nh_state(jg)%metrics%vertidx_gradp)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%zd_indlist)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"zd_indlist: ", &
        & p_nh_state(jg)%metrics%zd_indlist)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%zd_blklist)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"zd_blklist: ", &
        & p_nh_state(jg)%metrics%zd_blklist)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%zd_edgeidx)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"zd_edgeidx: ", &
        & p_nh_state(jg)%metrics%zd_edgeidx)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%zd_edgeblk)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"zd_edgeblk: ", &
        & p_nh_state(jg)%metrics%zd_edgeblk)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%zd_vertidx)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"zd_vertidx: ", &
        & p_nh_state(jg)%metrics%zd_vertidx)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%pg_edgeidx)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"pg_edgeidx: ", &
        & p_nh_state(jg)%metrics%pg_edgeidx)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%pg_edgeblk)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"pg_edgeblk: ", &
        & p_nh_state(jg)%metrics%pg_edgeblk)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%pg_vertidx)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"pg_vertidx: ", &
        & p_nh_state(jg)%metrics%pg_vertidx)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%nudge_c_idx)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"nudge_c_idx: ", &
        & p_nh_state(jg)%metrics%nudge_c_idx)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%nudge_e_idx)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"nudge_e_idx: ", &
        & p_nh_state(jg)%metrics%nudge_e_idx)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%nudge_c_blk)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"nudge_c_blk: ", &
        & p_nh_state(jg)%metrics%nudge_c_blk)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%nudge_e_blk)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"nudge_e_blk: ", &
        & p_nh_state(jg)%metrics%nudge_e_blk)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%bdy_halo_c_idx)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"bdy_halo_c_idx: ", &
        & p_nh_state(jg)%metrics%bdy_halo_c_idx)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%bdy_halo_c_blk)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"bdy_halo_c_blk: ", &
        & p_nh_state(jg)%metrics%bdy_halo_c_blk)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%ovlp_halo_c_dim)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"ovlp_halo_c_dim: ", &
        & p_nh_state(jg)%metrics%ovlp_halo_c_dim)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%ovlp_halo_c_idx)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"ovlp_halo_c_idx: ", &
        & p_nh_state(jg)%metrics%ovlp_halo_c_idx)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%ovlp_halo_c_blk)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"ovlp_halo_c_blk: ", &
        & p_nh_state(jg)%metrics%ovlp_halo_c_blk)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%bdy_mflx_e_idx)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"bdy_mflx_e_idx: ", &
        & p_nh_state(jg)%metrics%bdy_mflx_e_idx)
      IF(ASSOCIATED(p_nh_state(jg)%metrics%bdy_mflx_e_blk)) &
        & CALL printChecksum(prefix(1:pfx_tlen)//"bdy_mflx_e_blk: ", &
        & p_nh_state(jg)%metrics%bdy_mflx_e_blk)

      IF(PRESENT(p_lnd_state)) THEN
        IF(ALLOCATED(p_lnd_state(jg)%prog_lnd)) THEN
          DO i = 1, SIZE(p_lnd_state(jg)%prog_lnd(:), 1)
            WRITE (prefix, '(a,2(i0,a))') &
              "checksum of p_lnd_state(", jg, ")%prog_lnd(", i, ")%"
            pfx_tlen = LEN_TRIM(prefix)
            IF(ASSOCIATED(p_lnd_state(jg)%prog_lnd(i)%t_s_t)) &
              & CALL printChecksum(prefix(1:pfx_tlen)//"t_s_t: ", &
              & p_lnd_state(jg)%prog_lnd(i)%t_s_t)
            IF(ASSOCIATED(p_lnd_state(jg)%prog_lnd(i)%t_g)) &
              & CALL printChecksum(prefix(1:pfx_tlen)//"t_g: ", &
              & p_lnd_state(jg)%prog_lnd(i)%t_g)
            IF(ASSOCIATED(p_lnd_state(jg)%prog_lnd(i)%t_g_t)) &
              & CALL printChecksum(prefix(1:pfx_tlen)//"t_g_t: ", &
              & p_lnd_state(jg)%prog_lnd(i)%t_g_t)
            IF(ASSOCIATED(p_lnd_state(jg)%prog_lnd(i)%w_i_t)) &
              & CALL printChecksum(prefix(1:pfx_tlen)//"w_i_t: ", &
              & p_lnd_state(jg)%prog_lnd(i)%w_i_t)
            IF(ASSOCIATED(p_lnd_state(jg)%prog_lnd(i)%w_p_t)) &
              & CALL printChecksum(prefix(1:pfx_tlen)//"w_p_t: ", &
              & p_lnd_state(jg)%prog_lnd(i)%w_p_t)
            IF(ASSOCIATED(p_lnd_state(jg)%prog_lnd(i)%w_s_t)) &
              & CALL printChecksum(prefix(1:pfx_tlen)//"w_s_t: ", &
              & p_lnd_state(jg)%prog_lnd(i)%w_s_t)
            IF(ASSOCIATED(p_lnd_state(jg)%prog_lnd(i)%t_so_t)) &
              & CALL printChecksum(prefix(1:pfx_tlen)//"t_so_t: ", &
              & p_lnd_state(jg)%prog_lnd(i)%t_so_t)
            IF(ASSOCIATED(p_lnd_state(jg)%prog_lnd(i)%w_so_t)) &
              & CALL printChecksum(prefix(1:pfx_tlen)//"w_so_t: ", &
              & p_lnd_state(jg)%prog_lnd(i)%w_so_t)
            IF(ASSOCIATED(p_lnd_state(jg)%prog_lnd(i)%w_so_ice_t)) &
              & CALL printChecksum(prefix(1:pfx_tlen)//"w_so_ice_t: ", &
              & p_lnd_state(jg)%prog_lnd(i)%w_so_ice_t)
            IF(ASSOCIATED(p_lnd_state(jg)%prog_lnd(i)%t_snow_t)) &
              & CALL printChecksum(prefix(1:pfx_tlen)//"t_snow_t: ", &
              & p_lnd_state(jg)%prog_lnd(i)%t_snow_t)
            IF(ASSOCIATED(p_lnd_state(jg)%prog_lnd(i)%w_snow_t)) &
              & CALL printChecksum(prefix(1:pfx_tlen)//"w_snow_t: ", &
              & p_lnd_state(jg)%prog_lnd(i)%w_snow_t)
            IF(ASSOCIATED(p_lnd_state(jg)%prog_lnd(i)%rho_snow_t)) &
              & CALL printChecksum(prefix(1:pfx_tlen)//"rho_snow_t: ", &
              & p_lnd_state(jg)%prog_lnd(i)%rho_snow_t)
            IF(ASSOCIATED(p_lnd_state(jg)%prog_lnd(i)%t_snow_mult_t)) &
              & CALL printChecksum(prefix(1:pfx_tlen)//"t_snow_mult_t: ", &
              & p_lnd_state(jg)%prog_lnd(i)%t_snow_mult_t)
            IF(ASSOCIATED(p_lnd_state(jg)%prog_lnd(i)%wtot_snow_t)) &
              & CALL printChecksum(prefix(1:pfx_tlen)//"wtot_snow_t: ", &
              & p_lnd_state(jg)%prog_lnd(i)%wtot_snow_t)
            IF(ASSOCIATED(p_lnd_state(jg)%prog_lnd(i)%wliq_snow_t)) &
              & CALL printChecksum(prefix(1:pfx_tlen)//"wliq_snow_t: ", &
              & p_lnd_state(jg)%prog_lnd(i)%wliq_snow_t)
            IF(ASSOCIATED(p_lnd_state(jg)%prog_lnd(i)%rho_snow_mult_t)) &
              & CALL printChecksum(prefix(1:pfx_tlen)//"rho_snow_mult_t: ", &
              & p_lnd_state(jg)%prog_lnd(i)%rho_snow_mult_t)
            IF(ASSOCIATED(p_lnd_state(jg)%prog_lnd(i)%dzh_snow_t)) &
              & CALL printChecksum(prefix(1:pfx_tlen)//"dzh_snow_t: ", &
              & p_lnd_state(jg)%prog_lnd(i)%dzh_snow_t)
            !Can't checksum the t_ptr_2d3d fields because they
            !generally have ONLY one of the two pointers
            !initialized. Thus, checking the association
            !status of the pointers would RESULT IN undefined
            !behavior.
          END DO
        END IF
        IF(ALLOCATED(p_lnd_state(jg)%prog_wtr)) THEN
          DO i = 1, SIZE(p_lnd_state(jg)%prog_wtr(:), 1)
            WRITE (prefix, '(a,2(i0,a))') &
              "checksum of p_lnd_state(", jg, ")%prog_wtr(", i, ")%"
            pfx_tlen = LEN_TRIM(prefix)
            ! For some reason, these checksums explode with floating point exception.
            !   IF(ASSOCIATED(p_lnd_state(jg)%prog_wtr(i)%t_ice)) &
            !     &CALL printChecksum(prefix(1:pfx_tlen)//"t_ice: ", &
            !     & p_lnd_state(jg)%prog_wtr(i)%t_ice)
            !   IF(ASSOCIATED(p_lnd_state(jg)%prog_wtr(i)%h_ice)) &
            !     & CALL printChecksum(prefix(1:pfx_tlen)//"h_ice: ", &
            !     & p_lnd_state(jg)%prog_wtr(i)%h_ice)
            !   IF(ASSOCIATED(p_lnd_state(jg)%prog_wtr(i)%t_snow_si)) &
            !     & CALL printChecksum(prefix(1:pfx_tlen)//"t_snow_si: ", &
            !     & p_lnd_state(jg)%prog_wtr(i)%t_snow_si)
            !   IF(ASSOCIATED(p_lnd_state(jg)%prog_wtr(i)%h_snow_si)) &
            !     & CALL printChecksum(prefix(1:pfx_tlen)//"h_snow_si: ", &
            !     & p_lnd_state(jg)%prog_wtr(i)%h_snow_si)
            !   IF(ASSOCIATED(p_lnd_state(jg)%prog_wtr(i)%t_snow_lk)) &
            !     & CALL printChecksum(prefix(1:pfx_tlen)//"t_snow_lk: ", &
            !     & p_lnd_state(jg)%prog_wtr(i)%t_snow_lk)
            !   IF(ASSOCIATED(p_lnd_state(jg)%prog_wtr(i)%h_snow_lk)) &
            !     & CALL printChecksum(prefix(1:pfx_tlen)//"h_snow_lk: ", &
            !     & p_lnd_state(jg)%prog_wtr(i)%h_snow_lk)
            !   IF(ASSOCIATED(p_lnd_state(jg)%prog_wtr(i)%t_mnw_lk)) &
            !     & CALL printChecksum(prefix(1:pfx_tlen)//"t_mnw_lk: ", &
            !     & p_lnd_state(jg)%prog_wtr(i)%t_mnw_lk)
            !   IF(ASSOCIATED(p_lnd_state(jg)%prog_wtr(i)%t_wml_lk)) &
            !     & CALL printChecksum(prefix(1:pfx_tlen)//"t_wml_lk: ", &
            !     & p_lnd_state(jg)%prog_wtr(i)%t_wml_lk)
            !   IF(ASSOCIATED(p_lnd_state(jg)%prog_wtr(i)%h_ml_lk)) &
            !     & CALL printChecksum(prefix(1:pfx_tlen)//"h_ml_lk: ", &
            !     & p_lnd_state(jg)%prog_wtr(i)%h_ml_lk)
            !   IF(ASSOCIATED(p_lnd_state(jg)%prog_wtr(i)%t_bot_lk)) &
            !     & CALL printChecksum(prefix(1:pfx_tlen)//"t_bot_lk: ", &
            !     & p_lnd_state(jg)%prog_wtr(i)%t_bot_lk)
            !   IF(ASSOCIATED(p_lnd_state(jg)%prog_wtr(i)%c_t_lk)) &
            !     & CALL printChecksum(prefix(1:pfx_tlen)//"c_t_lk: ", &
            !     & p_lnd_state(jg)%prog_wtr(i)%c_t_lk)
            !   IF(ASSOCIATED(p_lnd_state(jg)%prog_wtr(i)%t_b1_lk)) &
            !     & CALL printChecksum(prefix(1:pfx_tlen)//"t_b1_lk: ", &
            !     & p_lnd_state(jg)%prog_wtr(i)%t_b1_lk)
            !   IF(ASSOCIATED(p_lnd_state(jg)%prog_wtr(i)%h_b1_lk)) &
            !     & CALL printChecksum(prefix(1:pfx_tlen)//"h_b1_lk: ", &
            !     & p_lnd_state(jg)%prog_wtr(i)%h_b1_lk)
          END DO
        END IF
        WRITE (prefix, '(a,i0,a)') &
          "checksum of p_lnd_state(", jg, ")%diag_lnd%"
        pfx_tlen = LEN_TRIM(prefix)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%qv_s)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"qv_s: ", &
          & p_lnd_state(jg)%diag_lnd%qv_s)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%t_s)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"t_s: ", &
          & p_lnd_state(jg)%diag_lnd%t_s)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%t_seasfc)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"t_seasfc: ", &
          & p_lnd_state(jg)%diag_lnd%t_seasfc)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%w_i)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"w_i: ", &
          & p_lnd_state(jg)%diag_lnd%w_i)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%w_p)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"w_p: ", &
          & p_lnd_state(jg)%diag_lnd%w_p)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%w_s)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"w_s: ", &
          & p_lnd_state(jg)%diag_lnd%w_s)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%t_so)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"t_so: ", &
          & p_lnd_state(jg)%diag_lnd%t_so)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%w_so)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"w_so: ", &
          & p_lnd_state(jg)%diag_lnd%w_so)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%w_so_ice)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"w_so_ice: ", &
          & p_lnd_state(jg)%diag_lnd%w_so_ice)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%runoff_s)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"runoff_s: ", &
          & p_lnd_state(jg)%diag_lnd%runoff_s)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%runoff_g)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"runoff_g: ", &
          & p_lnd_state(jg)%diag_lnd%runoff_g)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%resid_wso)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"resid_wso: ", &
          & p_lnd_state(jg)%diag_lnd%resid_wso)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%fr_seaice)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"fr_seaice: ", &
          & p_lnd_state(jg)%diag_lnd%fr_seaice)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%qv_s_t)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"qv_s_t: ", &
          & p_lnd_state(jg)%diag_lnd%qv_s_t)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%runoff_s_t)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"runoff_s_t: ", &
          & p_lnd_state(jg)%diag_lnd%runoff_s_t)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%runoff_g_t)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"runoff_g_t: ", &
          & p_lnd_state(jg)%diag_lnd%runoff_g_t)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%resid_wso_t)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"resid_wso_t: ", &
          & p_lnd_state(jg)%diag_lnd%resid_wso_t)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%rstom)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"rstom: ", &
          & p_lnd_state(jg)%diag_lnd%rstom)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%rstom_t)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"rstom_t: ", &
          & p_lnd_state(jg)%diag_lnd%rstom_t)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%t_snow)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"t_snow: ", &
          & p_lnd_state(jg)%diag_lnd%t_snow)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%rho_snow)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"rho_snow: ", &
          & p_lnd_state(jg)%diag_lnd%rho_snow)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%w_snow)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"w_snow: ", &
          & p_lnd_state(jg)%diag_lnd%w_snow)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%h_snow)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"h_snow: ", &
          & p_lnd_state(jg)%diag_lnd%h_snow)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%h_snow_t)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"h_snow_t: ", &
          & p_lnd_state(jg)%diag_lnd%h_snow_t)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%freshsnow)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"freshsnow: ", &
          & p_lnd_state(jg)%diag_lnd%freshsnow)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%freshsnow_t)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"freshsnow_t: ", &
          & p_lnd_state(jg)%diag_lnd%freshsnow_t)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%snowfrac)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"snowfrac: ", &
          & p_lnd_state(jg)%diag_lnd%snowfrac)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%snowfrac_lc)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"snowfrac_lc: ", &
          & p_lnd_state(jg)%diag_lnd%snowfrac_lc)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%snowfrac_t)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"snowfrac_t: ", &
          & p_lnd_state(jg)%diag_lnd%snowfrac_t)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%snowfrac_lc_t)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"snowfrac_lc_t: ", &
          & p_lnd_state(jg)%diag_lnd%snowfrac_lc_t)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%t_snow_mult)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"t_snow_mult: ", &
          & p_lnd_state(jg)%diag_lnd%t_snow_mult)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%rho_snow_mult)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"rho_snow_mult: ", &
          & p_lnd_state(jg)%diag_lnd%rho_snow_mult)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%wliq_snow)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"wliq_snow: ", &
          & p_lnd_state(jg)%diag_lnd%wliq_snow)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%wtot_snow)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"wtot_snow: ", &
          & p_lnd_state(jg)%diag_lnd%wtot_snow)
        IF(ASSOCIATED(p_lnd_state(jg)%diag_lnd%dzh_snow)) &
          & CALL printChecksum(prefix(1:pfx_tlen)//"dzh_snow: ", &
          & p_lnd_state(jg)%diag_lnd%dzh_snow)
        !Can't checksum the t_ptr_2d3d fields because they
        !generally have ONLY one of the two pointers
        !initialized. Thus, checking the association status of the
        !pointers would RESULT IN undefined behavior.
      END IF
    END DO
  END SUBROUTINE printChecksums

  SUBROUTINE init_qnx_from_qx_twomom (caller, p_patch, p_prog, lqnx_init)

    CHARACTER(len=*), INTENT(in)    :: caller   ! Name of calling routine for messages
    TYPE(t_patch)   , INTENT(in)    :: p_patch
    TYPE(t_nh_prog) , INTENT(inout) :: p_prog
    LOGICAL         , INTENT(in)    :: lqnx_init(:)  ! List of switches for each hydrometeor, if its
                                                     ! number conc. should be initialized from its mass conc.
                                                     ! The corresponding indices into this list have to be iqnc, iqni, iqnr, etc.

    INTEGER                             :: jb, jk, jc, nlen
    REAL(wp), POINTER, DIMENSION(:,:,:) :: my_qc,  my_qi,  my_qr,  my_qs,  my_qg,  my_qh, my_rho, &
                                           my_qnc, my_qni, my_qnr, my_qns, my_qng, my_qnh
    CHARACTER(len=110)                  :: ncmaxstr
    REAL(wp)                            :: rholoc

    my_rho => p_prog%rho(:,:,:)
    my_qc  => p_prog%tracer(:,:,:,iqc)
    my_qi  => p_prog%tracer(:,:,:,iqi)
    my_qr  => p_prog%tracer(:,:,:,iqr)
    my_qs  => p_prog%tracer(:,:,:,iqs)
    my_qg  => p_prog%tracer(:,:,:,iqg)
    my_qh  => p_prog%tracer(:,:,:,iqh)

    my_qnc => p_prog%tracer(:,:,:,iqnc)
    my_qni => p_prog%tracer(:,:,:,iqni)
    my_qnr => p_prog%tracer(:,:,:,iqnr)
    my_qns => p_prog%tracer(:,:,:,iqns)
    my_qng => p_prog%tracer(:,:,:,iqng)
    my_qnh => p_prog%tracer(:,:,:,iqnh)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,nlen,rholoc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, p_patch%nblks_c
      nlen = MERGE(nproma, p_patch%npromz_c, jb /= p_patch%nblks_c)

      IF (lqnx_init(iqnc)) THEN
        DO jk = 1, p_patch%nlev
          DO jc = 1, nlen
            rholoc = MAX(my_rho(jc,jk,jb), 1e-20_wp)
            my_qnc(jc,jk,jb) = set_qnc( MAX(my_qc(jc,jk,jb), 0.0_wp)*rholoc ) / rholoc
          END DO
        END DO
      END IF
      IF (lqnx_init(iqni)) THEN
        DO jk = 1, p_patch%nlev
          DO jc = 1, nlen
            rholoc = MAX(my_rho(jc,jk,jb), 1e-20_wp)
            my_qni(jc,jk,jb) = set_qni( MAX(my_qi(jc,jk,jb), 0.0_wp)*rholoc ) / rholoc
          END DO
        END DO
      END IF
      IF (lqnx_init(iqnr)) THEN
        DO jk = 1, p_patch%nlev
          DO jc = 1, nlen
            rholoc = MAX(my_rho(jc,jk,jb), 1e-20_wp)
            my_qnr(jc,jk,jb) = set_qnr( MAX(my_qr(jc,jk,jb), 0.0_wp)*rholoc ) / rholoc
          END DO
        END DO
      END IF
      IF (lqnx_init(iqns)) THEN
        DO jk = 1, p_patch%nlev
          DO jc = 1, nlen
            rholoc = MAX(my_rho(jc,jk,jb), 1e-20_wp)
            my_qns(jc,jk,jb) = set_qns( MAX(my_qs(jc,jk,jb), 0.0_wp)*rholoc ) / rholoc
          END DO
        END DO
      END IF
      IF (lqnx_init(iqng)) THEN
        DO jk = 1, p_patch%nlev
          DO jc = 1, nlen
            rholoc = MAX(my_rho(jc,jk,jb), 1e-20_wp)
            my_qng(jc,jk,jb) = set_qng( MAX(my_qg(jc,jk,jb), 0.0_wp)*rholoc ) / rholoc
          END DO
        END DO
      END IF
      IF (lqnx_init(iqnh)) THEN
        DO jk = 1, p_patch%nlev
          DO jc = 1, nlen
            rholoc = MAX(my_rho(jc,jk,jb), 1e-20_wp)
            ! init qnh assuming hail bulk density of 750 kg/m^3 and an
            !  exponential PSD w.r.t. diameter with const. intercept N0 of 1e6 m-4:
            my_qnh(jc,jk,jb) = set_qnh_expPSD_N0const( &
              MAX(my_qh(jc,jk,jb), 0.0_wp)*rholoc, 750.0_wp, 1.0e6_wp ) / rholoc
          END DO
        END DO
      END IF
    END DO
!$OMP END DO
!$OMP END PARALLEL

    IF (lqnx_init(iqnc)) THEN
      ncmaxstr = get_diag_stat_str_3d ( p_patch, my_qnc )
      CALL message(caller, 'init_qnx_from_qx_twomom: set_qnc() from qc, '//TRIM(ncmaxstr))
    END IF
    IF (lqnx_init(iqni)) THEN
      ncmaxstr = get_diag_stat_str_3d ( p_patch, my_qni )
      CALL message(caller, 'init_qnx_from_qx_twomom: set_qni() from qi, '//TRIM(ncmaxstr))
    END IF
    IF (lqnx_init(iqnr)) THEN
      ncmaxstr = get_diag_stat_str_3d ( p_patch, my_qnr )
      CALL message(caller, 'init_qnx_from_qx_twomom: set_qnr() from qr, '//TRIM(ncmaxstr))
    END IF
    IF (lqnx_init(iqns)) THEN
      ncmaxstr = get_diag_stat_str_3d ( p_patch, my_qns )
      CALL message(caller, 'init_qnx_from_qx_twomom: set_qns() from qs, '//TRIM(ncmaxstr))
    END IF
    IF (lqnx_init(iqng)) THEN
      ncmaxstr = get_diag_stat_str_3d ( p_patch, my_qng )
      CALL message(caller, 'init_qnx_from_qx_twomom: set_qng() from qg, '//TRIM(ncmaxstr))
    END IF
    IF (lqnx_init(iqnh)) THEN
      ncmaxstr = get_diag_stat_str_3d ( p_patch, my_qnh )
      CALL message(caller, 'init_qnx_from_qx_twomom: set_qnh() from qh, '//TRIM(ncmaxstr))
    END IF

  END SUBROUTINE init_qnx_from_qx_twomom

  SUBROUTINE init_qnxinc_from_qxinc_twomom (caller, p_patch, p_prog, initicon, lqx_avail, lqxinc_avail, lqnxinc_init)

    CHARACTER(len=*), INTENT(in)    :: caller   ! Name of calling routine for messages
    TYPE(t_patch)   , INTENT(in)    :: p_patch
    TYPE(t_nh_prog) , INTENT(in)    :: p_prog
    TYPE(t_initicon_state), INTENT(inout), TARGET    :: initicon
    LOGICAL         , INTENT(in)    :: lqx_avail(:), lqxinc_avail(:), lqnxinc_init(:)

    INTEGER                             :: jb, jk, jc, nlen
    REAL(wp), POINTER, DIMENSION(:,:,:) :: my_qc,  my_qi,  my_qr,  my_qs,  my_qg,  my_qh, my_rho, &
                                           my_qc_inc, my_qi_inc, my_qr_inc, my_qs_inc, my_qg_inc, my_qh_inc, &
                                           my_qnc_inc, my_qni_inc, my_qnr_inc, my_qns_inc, my_qng_inc, my_qnh_inc, &
                                           my_qnc,  my_qni,  my_qnr,  my_qns,  my_qng,  my_qnh
    REAL(wp)                            :: qtmp0, qtmp1, rholoc, meanmass
    CHARACTER(len=110)                  :: ncmaxstr

    REAL(wp), PARAMETER                     :: qc_xmax = 2.60e-10_wp 
    REAL(wp), PARAMETER                     :: qc_xmin = 4.20e-15_wp
    REAL(wp), PARAMETER                     :: qr_xmax = 3.00e-06_wp
    REAL(wp), PARAMETER                     :: qr_xmin = 2.60e-10_wp
    REAL(wp), PARAMETER                     :: qi_xmax = 1.00e-05_wp
    REAL(wp), PARAMETER                     :: qi_xmin = 1.00e-12_wp
    REAL(wp), PARAMETER                     :: qs_xmax = 2.00e-05_wp
    REAL(wp), PARAMETER                     :: qs_xmin = 1.00e-10_wp
    REAL(wp), PARAMETER                     :: qg_xmax = 5.30e-04_wp
    REAL(wp), PARAMETER                     :: qg_xmin = 4.19e-09_wp
    REAL(wp), PARAMETER                     :: qh_xmax = 5.00e-03_wp
    REAL(wp), PARAMETER                     :: qh_xmin = 2.60e-09_wp
    REAL(wp), PARAMETER                     :: myeps_q = 1.0e-6_wp
    REAL(wp), PARAMETER                     :: myeps_n = 0.1_wp

    my_rho => p_prog%rho(:,:,:)
    my_qc  => p_prog%tracer(:,:,:,iqc)
    my_qi  => p_prog%tracer(:,:,:,iqi)
    my_qr  => p_prog%tracer(:,:,:,iqr)
    my_qs  => p_prog%tracer(:,:,:,iqs)
    my_qg  => p_prog%tracer(:,:,:,iqg)
    my_qh  => p_prog%tracer(:,:,:,iqh)

    my_qnc  => p_prog%tracer(:,:,:,iqnc)
    my_qni  => p_prog%tracer(:,:,:,iqni)
    my_qnr  => p_prog%tracer(:,:,:,iqnr)
    my_qns  => p_prog%tracer(:,:,:,iqns)
    my_qng  => p_prog%tracer(:,:,:,iqng)
    my_qnh  => p_prog%tracer(:,:,:,iqnh)

    my_qc_inc => initicon%atm_inc%qc
    my_qi_inc => initicon%atm_inc%qi
    my_qr_inc => initicon%atm_inc%qr
    my_qs_inc => initicon%atm_inc%qs
    my_qg_inc => initicon%atm_inc%qg
    my_qh_inc => initicon%atm_inc%qh

    my_qnc_inc => initicon%atm_inc%qnc
    my_qni_inc => initicon%atm_inc%qni
    my_qnr_inc => initicon%atm_inc%qnr
    my_qns_inc => initicon%atm_inc%qns
    my_qng_inc => initicon%atm_inc%qng
    my_qnh_inc => initicon%atm_inc%qnh

            ! Alberto: performance could be improved by defining 1ovmeanmas

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,nlen,qtmp0,qtmp1,rholoc,meanmass) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, p_patch%nblks_c

      nlen = MERGE(nproma, p_patch%npromz_c, jb /= p_patch%nblks_c)

      IF (lqnxinc_init(iqnc) .AND. lqx_avail(iqc) .AND. lqxinc_avail(iqc) .AND. qcana_mode > 0) THEN
        DO jk = 1, p_patch%nlev
          DO jc = 1, nlen
            IF ( my_qc(jc,jk,jb) > myeps_q .AND. my_qnc(jc,jk,jb) > myeps_n) THEN
              meanmass = MAX(MIN(my_qc(jc,jk,jb)/(my_qnc(jc,jk,jb)+myeps_n),qc_xmax),qc_xmin ) 
              my_qnc_inc(jc,jk,jb) = my_qc_inc(jc,jk,jb) / meanmass 
            ELSE
              qtmp1 = MAX( my_qc(jc,jk,jb) + my_qc_inc(jc,jk,jb) , 0.0_wp)
              qtmp0 = MAX( my_qc(jc,jk,jb) , 0.0_wp)
              rholoc = MAX(my_rho(jc,jk,jb), 1e-20_wp)
              my_qnc_inc(jc,jk,jb) = ( set_qnc( qtmp1*rholoc ) - set_qnc( qtmp0*rholoc ) ) / rholoc 
            END IF
          END DO
        END DO
      END IF

      IF (lqnxinc_init(iqni) .AND. lqx_avail(iqi) .AND. lqxinc_avail(iqi) .AND. qiana_mode > 0) THEN
        DO jk = 1, p_patch%nlev
          DO jc = 1, nlen
            IF ( my_qi(jc,jk,jb) > myeps_q .AND. my_qni(jc,jk,jb) > myeps_n) THEN
              meanmass = MAX(MIN(my_qi(jc,jk,jb)/(my_qni(jc,jk,jb)+myeps_n),qi_xmax),qi_xmin ) 
              my_qni_inc(jc,jk,jb) = my_qi_inc(jc,jk,jb) / meanmass 
            ELSE
              qtmp1 = MAX( my_qi(jc,jk,jb) + my_qi_inc(jc,jk,jb) , 0.0_wp)
              qtmp0 = MAX( my_qi(jc,jk,jb) , 0.0_wp)
              rholoc = MAX(my_rho(jc,jk,jb), 1e-20_wp)
              my_qni_inc(jc,jk,jb) = ( set_qni( qtmp1*rholoc ) - set_qni( qtmp0*rholoc ) ) / rholoc 
            END IF
          END DO
        END DO
      END IF

      IF (lqnxinc_init(iqnr) .AND. lqx_avail(iqr) .AND. lqxinc_avail(iqr) .AND. qrsgana_mode > 0) THEN
        DO jk = 1, p_patch%nlev
          DO jc = 1, nlen
            IF ( my_qr(jc,jk,jb) > myeps_q .AND. my_qnr(jc,jk,jb) > myeps_n) THEN
              meanmass = MAX(MIN(my_qr(jc,jk,jb)/(my_qnr(jc,jk,jb)+myeps_n),qr_xmax),qr_xmin ) 
              my_qnr_inc(jc,jk,jb) = my_qr_inc(jc,jk,jb) / meanmass 
            ELSE
              qtmp1 = MAX( my_qr(jc,jk,jb) + my_qr_inc(jc,jk,jb) , 0.0_wp)
              qtmp0 = MAX( my_qr(jc,jk,jb) , 0.0_wp)
              rholoc = MAX(my_rho(jc,jk,jb), 1e-20_wp)
              my_qnr_inc(jc,jk,jb) = ( set_qnr( qtmp1*rholoc ) - set_qnr( qtmp0*rholoc ) ) / rholoc 
            END IF
          END DO
        END DO
      END IF

      IF (lqnxinc_init(iqns) .AND. lqx_avail(iqs) .AND. lqxinc_avail(iqs) .AND. qrsgana_mode > 0) THEN
        DO jk = 1, p_patch%nlev
          DO jc = 1, nlen
            IF ( my_qs(jc,jk,jb) > myeps_q .AND. my_qns(jc,jk,jb) > myeps_n) THEN
              meanmass = MAX(MIN(my_qs(jc,jk,jb)/(my_qns(jc,jk,jb)+myeps_n),qs_xmax),qs_xmin ) 
              my_qns_inc(jc,jk,jb) = my_qs_inc(jc,jk,jb) / meanmass 
            ELSE
              qtmp1 = MAX( my_qs(jc,jk,jb) + my_qs_inc(jc,jk,jb) , 0.0_wp)
              qtmp0 = MAX( my_qs(jc,jk,jb) , 0.0_wp)
              rholoc = MAX(my_rho(jc,jk,jb), 1e-20_wp)
              my_qns_inc(jc,jk,jb) = ( set_qns( qtmp1*rholoc ) - set_qns( qtmp0*rholoc ) ) / rholoc 
            END IF
          END DO
        END DO
      END IF

      IF (lqnxinc_init(iqng) .AND. lqx_avail(iqg) .AND. lqxinc_avail(iqg) .AND. qrsgana_mode > 0) THEN
        DO jk = 1, p_patch%nlev
          DO jc = 1, nlen
            IF ( my_qg(jc,jk,jb) > myeps_q .AND. my_qng(jc,jk,jb) > myeps_n) THEN
              meanmass = MAX(MIN(my_qg(jc,jk,jb)/(my_qng(jc,jk,jb)+myeps_n),qg_xmax),qg_xmin ) 
              my_qng_inc(jc,jk,jb) = my_qg_inc(jc,jk,jb) / meanmass 
            ELSE
              qtmp1 = MAX( my_qg(jc,jk,jb) + my_qg_inc(jc,jk,jb) , 0.0_wp)
              qtmp0 = MAX( my_qg(jc,jk,jb) , 0.0_wp)
              rholoc = MAX(my_rho(jc,jk,jb), 1e-20_wp)
              my_qng_inc(jc,jk,jb) = ( set_qng( qtmp1*rholoc ) - set_qng( qtmp0*rholoc ) ) / rholoc 
            END IF
          END DO
        END DO
      END IF

      IF (lqnxinc_init(iqnh) .AND. lqx_avail(iqh) .AND. lqxinc_avail(iqh) .AND. qrsgana_mode > 0) THEN
        DO jk = 1, p_patch%nlev
          DO jc = 1, nlen
            IF ( my_qh(jc,jk,jb) > myeps_q .AND. my_qnh(jc,jk,jb) > myeps_n) THEN
              meanmass = MAX(MIN(my_qh(jc,jk,jb)/(my_qnh(jc,jk,jb)+myeps_n),qh_xmax),qh_xmin ) 
              my_qnh_inc(jc,jk,jb) = my_qh_inc(jc,jk,jb) / meanmass 
            ELSE
              qtmp1 = MAX( my_qh(jc,jk,jb) + my_qh_inc(jc,jk,jb) , 0.0_wp)
              qtmp0 = MAX( my_qh(jc,jk,jb) , 0.0_wp)
              rholoc = MAX(my_rho(jc,jk,jb), 1e-20_wp)
!              my_qnh_inc(jc,jk,jb) = ( set_qnh( qtmp1*rholoc ) - set_qnh( qtmp0*rholoc ) ) / rholoc 
              ! init qnh increment by assuming hail bulk density of 750 kg/m^3 and an
              !  exponential PSD w.r.t. diameter with const. intercept N0 of 1e6 m-4:
              my_qnh_inc(jc,jk,jb) = ( set_qnh_expPSD_N0const( qtmp1*rholoc, 750.0_wp, 1.0e6_wp ) - &
                                       set_qnh_expPSD_N0const( qtmp0*rholoc, 750.0_wp, 1.0e6_wp ) ) / rholoc
            END IF
          END DO
        END DO
      END IF

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    
    IF (qcana_mode > 0 .AND. lqnxinc_init(iqnc)) THEN
      IF (lqx_avail(iqc) .AND. lqxinc_avail(iqc)) THEN
        ncmaxstr = get_diag_stat_str_3d ( p_patch, my_qnc_inc )
        CALL message(caller, 'init_qnxinc_from_qxinc_twomom for IAU: set_qnc_inc() from qc and qcinc, '//TRIM(ncmaxstr))
      ELSE
        CALL message(caller, 'init_qnxinc_from_qxinc_twomom for IAU: set_qnc_inc() failed '// &
                                   'due to missing qc (FG) and/or qcinc (ANA)')        
      END IF
    END IF
    IF (qiana_mode > 0 .AND. lqnxinc_init(iqni)) THEN
      IF (lqx_avail(iqi) .AND. lqxinc_avail(iqi)) THEN
        ncmaxstr = get_diag_stat_str_3d ( p_patch, my_qni_inc )
        CALL message(caller, 'init_qnxinc_from_qxinc_twomom for IAU: set_qni_inc() from qi and qiinc, '//TRIM(ncmaxstr))
      ELSE
        CALL message(caller, 'init_qnxinc_from_qxinc_twomom for IAU: set_qni_inc() failed '// &
                                   'due to missing qi (FG) and/or qiinc (ANA)')        
      END IF
    END IF
    IF (qrsgana_mode > 0) THEN
      IF (lqnxinc_init(iqnr)) THEN
        IF (lqx_avail(iqr) .AND. lqxinc_avail(iqr)) THEN
          ncmaxstr = get_diag_stat_str_3d ( p_patch, my_qnr_inc )
          CALL message(caller, 'init_qnxinc_from_qxinc_twomom for IAU: set_qnr_inc() from qr and qrinc, '//TRIM(ncmaxstr))
        ELSE
          CALL message(caller, 'init_qnxinc_from_qxinc_twomom for IAU: set_qnr_inc() failed '// &
                                     'due to missing qr (FG) and/or qrinc (ANA)')
        END IF
      END IF
      IF (lqnxinc_init(iqns)) THEN
        IF (lqx_avail(iqs) .AND. lqxinc_avail(iqs)) THEN
          ncmaxstr = get_diag_stat_str_3d ( p_patch, my_qns_inc )
          CALL message(caller, 'init_qnxinc_from_qxinc_twomom for IAU: set_qns_inc() from qs and qsinc, '//TRIM(ncmaxstr))
        ELSE
          CALL message(caller, 'init_qnxinc_from_qxinc_twomom for IAU: set_qns_inc() failed '// &
                                     'due to missing qs (FG) and/or qsinc (ANA)')        
        END IF
      END IF
      IF (lqnxinc_init(iqng)) THEN
        IF (lqx_avail(iqg) .AND. lqxinc_avail(iqg)) THEN
          ncmaxstr = get_diag_stat_str_3d ( p_patch, my_qng_inc )
          CALL message(caller, 'init_qnxinc_from_qxinc_twomom for IAU: set_qng_inc() from qg and qginc, '//TRIM(ncmaxstr))
        ELSE
          CALL message(caller, 'init_qnxinc_from_qxinc_twomom for IAU: set_qng_inc() failed '// &
                                     'due to missing qg (FG) and/or qginc (ANA)')        
        END IF
      END IF
      IF (lqnxinc_init(iqnh)) THEN
        IF (lqx_avail(iqh) .AND. lqxinc_avail(iqh)) THEN
          ncmaxstr = get_diag_stat_str_3d ( p_patch, my_qnh_inc )
          CALL message(caller, 'init_qnxinc_from_qxinc_twomom for IAU: set_qnh_inc() from qh and qhinc, '//TRIM(ncmaxstr))
        ELSE
          CALL message(caller, 'init_qnxinc_from_qxinc_twomom for IAU: set_qnh_inc() failed '// &
                                     'due to missing qh (FG) and/or qhinc (ANA)')        
        END IF
      END IF
    END IF

  END SUBROUTINE init_qnxinc_from_qxinc_twomom

  !> Diagnostic min/mean/max in a character representation of a 3D-field in the interior domain (no boundaries and halos)
  FUNCTION get_diag_stat_str_3d ( p_patch_jg, field3d ) RESULT(statstr)

    REAL(wp), INTENT(in) :: field3d(:,:,:)
    TYPE(t_patch)        :: p_patch_jg

    REAL(wp)             :: stats(3)

    CHARACTER(len=110)   :: statstr

    stats = get_diag_stat_comm_work ( p_patch_jg, field3d )

    IF (my_process_is_work()) THEN

      WRITE(statstr, '("min/mean/max interior domain: ",2(es12.5," / "),es12.5)') stats

    ELSE

      statstr = modname//': get_diag_stat_str_3d(): not a worker PE, no min/mean/max available!'

    END IF

  END FUNCTION get_diag_stat_str_3d

  !> Diagnostic min/mean/max values of a 3D-field in the interior domain (no boundaries and halos)
  FUNCTION get_diag_stat_comm_work ( p_patch_jg, field3d ) RESULT(stats)

    REAL(wp), INTENT(in) :: field3d(:,:,:)
    TYPE(t_patch)        :: p_patch_jg

    REAL(wp)             :: stats(3)

    REAL(wp)             :: mn, mx, mm
    INTEGER              :: size_field3d
    INTEGER              :: jb, jc, jk, is, ie, i_startblk, i_endblk

    IF (my_process_is_work()) THEN

      mn = HUGE(1.0_wp)
      mx = -HUGE(1.0_wp)
      mm = 0.0_wp
      size_field3d = 0

      i_startblk = p_patch_jg % cells % start_block(grf_bdywidth_c+1)   ! interior cells
      i_endblk   = p_patch_jg % cells % end_block(min_rlcell_int)       ! excluding halo cells

!$OMP PARALLEL PRIVATE(jb, jc, jk, is, ie) REDUCTION(+:size_field3d,mm) REDUCTION(min:mn) REDUCTION(max:mx)
!$OMP DO
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch_jg, jb, i_startblk, i_endblk, is, ie, grf_bdywidth_c+1, min_rlcell_int)
        DO jk = 1, p_patch_jg%nlev
          DO jc = is, ie
            mn           = MIN(field3d(jc,jk,jb), mn)
            mx           = MAX(field3d(jc,jk,jb), mx)
            mm           = field3d(jc,jk,jb) + mm
            size_field3d = size_field3d + 1
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL




      mm = mm / MAX(size_field3d, 1)

      stats = (/mn, mm, mx/)

    ELSE

      stats(:) = -HUGE(1.0_wp)

    END IF

  END FUNCTION get_diag_stat_comm_work


  !-------------------------------------------------------------------------
  !>
  !! for coupled simulation with new land points forced from the ocean: initialize soil
  !!
  !-------------------------------------------------------------------------

  SUBROUTINE new_land_from_ocean (p_patch, p_nh_state, p_lnd_state, ext_data)

    TYPE(t_patch)                 ,INTENT(IN)    :: p_patch(:)
    TYPE(t_nh_state) , TARGET     ,INTENT(IN)    :: p_nh_state(:)
    TYPE(t_lnd_state), TARGET     ,INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_external_data)         ,INTENT(IN)    :: ext_data(:)

    INTEGER :: jg, jb, jt, jk, jc, nlev            ! loop indices
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk 
    INTEGER :: i_startidx, i_endidx
    INTEGER :: ist

    TYPE(t_nh_diag) , POINTER :: p_diag            ! shortcut to diag state
    TYPE(t_lnd_prog), POINTER :: lnd_prog_now      ! shortcut to prognostic land state

    !-----------------------------------------------------------------------

    CALL message('','Initialize soil w_so and t_so for coupled simulation new land points.')

    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      rl_start  = 1
      rl_end    = min_rlcell

      i_startblk = p_patch(jg)%cells%start_block(rl_start)
      i_endblk   = p_patch(jg)%cells%end_block(rl_end)

      p_diag       => p_nh_state (jg)%diag
      lnd_prog_now => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))

      nlev         =  p_patch(jg)%nlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jt, jk, jc, ist,i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        DO jt = 1, ntiles_total
          DO jk = 1, nlev_soil+1
            DO jc = i_startidx, i_endidx
 
              IF ( ext_data(jg)%atm%lsm_switch(jc,jb) == 1 ) THEN
                ist = ext_data(jg)%atm%soiltyp(jc,jb)

                IF ( jk <= nlev_soil ) THEN  ! t_so_t has one level more than w_so_t
                  ! set default soil moisture between field capacity and wilting point
                  lnd_prog_now%w_so_t(jc,jk,jb,jt) = 0.5_wp * ( cfcap(ist) + cpwp(ist) ) * dzsoil(jk)
                  !          ( cporv(ist) + cadp(ist)                ! pore capacity + field capacity
                ENDIF

                lnd_prog_now%t_so_t(jc,jk,jb,jt) = p_diag%temp(jc,nlev,jb)
 
              ENDIF

            ENDDO
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
    ENDDO

  END SUBROUTINE new_land_from_ocean

  !>
  !! Interface for calling diagnose condensate
  !!
  SUBROUTINE prepare_thermo_src_term(p_patch)

    TYPE(t_patch), TARGET, INTENT(INOUT) :: p_patch

    INTEGER :: i_startidx, i_endidx, jc, jk, jt, jb, nlev
    INTEGER :: i_rlstart , i_rlend , i_startblk, i_endblk
    INTEGER :: nn, jg

    REAL(wp), POINTER :: p_csum(:,:)

    REAL(wp), TARGET :: qsum_liq(nproma,p_patch%nlev), qsum_ice(nproma,p_patch%nlev)
    REAL(wp) :: z_a, z_b


    jg = p_patch%id
    nn = nnow_rcf(jg)

    ! number of vertical levels
    nlev = p_patch%nlev

    ! boundary zone and halo points are not needed
    i_rlstart  = grf_bdywidth_c + 1
    i_rlend    = min_rlcell_int
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    !$ACC DATA CREATE(qsum_liq, qsum_ice) &
    !$ACC   COPYIN(kstart_moist) &
    !$ACC   IF(i_am_accel_node)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jk,jt,jb,i_startidx,i_endidx,p_csum,z_a,z_b,qsum_liq,qsum_ice) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
      &                   i_startidx, i_endidx, i_rlstart, i_rlend )

      !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx

          qsum_liq(jc,jk) = 0.0_wp
          qsum_ice(jc,jk) = 0.0_wp

        END DO ! jc
      END DO ! jk
      !$ACC END PARALLEL

      DO jt = 1, ntracer

        IF (ANY((/iqc,iqr,iqhl,iqgl/)==jt)) THEN
          p_csum => qsum_liq
        ELSE IF (ANY((/iqi,iqs,iqg,iqh/)==jt)) THEN
          p_csum => qsum_ice
        ELSE
          CYCLE
        END IF

        !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
         DO jk = kstart_moist(jg),nlev
           DO jc = i_startidx, i_endidx

             p_csum(jc,jk) = p_csum(jc,jk) + p_nh_state(jg)%prog(nn)%tracer(jc,jk,jb,jt)

           END DO ! jc
         END DO ! jk
        !$ACC END PARALLEL

      END DO ! jt

      !$ACC PARALLEL IF(i_am_accel_node) DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(z_a, z_b)
      DO jk = kstart_moist(jg),nlev
        DO jc = i_startidx, i_endidx

          z_a = 1._wp - ( p_nh_state(jg)%prog(nn)%tracer(jc,jk,jb,iqv) &
                          + qsum_liq(jc,jk) + qsum_ice(jc,jk) )

          z_b = clw * qsum_liq(jc,jk) + ci * qsum_ice(jc,jk)

          p_nh_state(jg)%diag%chi_q(jc,jk,jb) = 1._wp - (                                   &
              ( z_a + ( z_b  + cpv * p_nh_state(jg)%prog(nn)%tracer(jc,jk,jb,iqv) ) / cpd ) &
            / ( z_a + ( z_b  + cvv * p_nh_state(jg)%prog(nn)%tracer(jc,jk,jb,iqv) ) / cvd ) )

        END DO ! jc
      END DO ! jk
      !$ACC END PARALLEL

    END DO ! jb
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE prepare_thermo_src_term

END MODULE mo_initicon_utils
!
! Local Variables:
! f90-continuation-indent: 2
! End:
!
