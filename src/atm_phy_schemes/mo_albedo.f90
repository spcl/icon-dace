!
! Calculation of surface albedo
!
! Calculation of surface albedo taking soil type, vegetation
! and snow/ice conditions into account.
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

MODULE mo_albedo

  USE mo_exception,            ONLY: finish
  USE mo_kind,                 ONLY: wp
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_model_domain,         ONLY: t_patch
  USE mo_nwp_lnd_types,        ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_parallel_config,      ONLY: nproma
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_radiation_config,     ONLY: rad_csalbw, direct_albedo, direct_albedo_water, albedo_whitecap
  USE mo_lnd_nwp_config,       ONLY: ntiles_total, ntiles_water, ntiles_lnd, lterra_urb, lurbalb,  &
    &                                lseaice, lprog_albsi, llake, isub_water, isub_lake,           &
    &                                isub_seaice
  USE mo_extpar_config,        ONLY: itype_vegetation_cycle
  USE sfc_terra_data,          ONLY: csalb, csalb_snow_fe, csalb_snow_fd,     &
    &                                csalb_snow_min, csalb_snow_max, csalb_p, csalb_snow, &
    &                                ist_seawtr, ist_seaice
  USE mo_physical_constants,   ONLY: tmelt, tf_salt
  USE sfc_flake_data,          ONLY: albedo_whiteice_ref, albedo_blueice_ref, &
    &                                c_albice_MR, tpl_T_f, h_Ice_min_flk
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_impl_constants,       ONLY: min_rlcell_int, LSS_JSBACH, LSS_TERRA
  USE sfc_seaice,              ONLY: alb_seaice_equil
  USE mo_initicon_config,      ONLY: icpl_da_snowalb

  IMPLICIT NONE

  PRIVATE



  PUBLIC  :: sfc_albedo
  PUBLIC  :: sfc_albedo_modis
  PUBLIC  :: sfc_albedo_scm

  PUBLIC :: sfc_albedo_dir_rg
  PUBLIC :: sfc_albedo_dir_yang
  PUBLIC :: sfc_albedo_dir_taylor
  PUBLIC :: sfc_albedo_whitecap

CONTAINS


  !>
  !! Calculation of surface albedo
  !!
  !! Calculation of surface albedo based on tabulated shortwave bare soil 
  !! albedo data. In addition, soil type, vegetation and snow/ice conditions 
  !! are taken into account
  !!
  SUBROUTINE sfc_albedo(pt_patch, ext_data, lnd_prog, wtr_prog, lnd_diag, prm_diag, lacc)

    TYPE(t_patch),          INTENT(   in):: pt_patch  !< grid/patch info.

    TYPE(t_external_data),  INTENT(   in):: ext_data  !< external data

    TYPE(t_lnd_prog),       INTENT(   in):: lnd_prog  !< land prognostic state (new)

    TYPE(t_wtr_prog),       INTENT(   in):: wtr_prog  !< water prognostic state (new)

    TYPE(t_lnd_diag),       INTENT(   in):: lnd_diag  !< land diagnostic state

    TYPE(t_nwp_phy_diag),   INTENT(inout):: prm_diag

    LOGICAL,                INTENT(   in):: lacc      !< accelerator flag

    ! Local scalars:
    REAL(wp):: zvege                   !< plant cover fraction
    REAL(wp):: zsnow                   !< snow cover fraction
    REAL(wp):: zminsnow_alb            !< temperature-dependent minimum snow albedo
    REAL(wp):: t_fac                   !< factor for temperature dependency of zminsnow_alb over glaciers
    REAL(wp):: zsalb_snow              !< snow albedo (predictor)
    REAL(wp):: zsnow_alb               !< snow albedo (corrector)
    REAL(wp):: wc_fraction             !< whitecap fraction
    REAL(wp):: wc_albedo               !< whitecap albedo

    INTEGER :: jg                      !< patch ID
    INTEGER :: jb, jc, ic, jt          !< loop indices
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index
    INTEGER :: ist
    INTEGER :: ilu                     !< land cover class
    INTEGER :: i_count_lnd             !< number of land points
    INTEGER :: i_count_sea             !< number of sea points
    INTEGER :: i_count_flk             !< number of lake points
    INTEGER :: i_count_seaice          !< number of seaice points

    !-----------------------------------------------------------------------

    jg = pt_patch%id

    IF (atm_phy_nwp_config(jg)%inwp_surface == LSS_JSBACH) THEN
      ! Albedo is already up-to-date.
      RETURN
    END IF


    i_nchdom  = MAX(1,pt_patch%n_childdom)

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jt,jc,ic,i_startidx,i_endidx,ist,zvege,zsnow,  &
!$OMP            zsalb_snow,zsnow_alb,ilu,i_count_lnd,i_count_sea, &
!$OMP            i_count_flk,i_count_seaice,zminsnow_alb,t_fac,    &
!$OMP            wc_fraction, wc_albedo) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = i_startblk, i_endblk


      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        &                 i_startidx, i_endidx, rl_start, rl_end)


      !------------------------------------------------------------------------------
      ! Calculation of surface albedo taking soil type,              
      ! vegetation and snow/ice conditions into account
      !------------------------------------------------------------------------------
      
      IF ( atm_phy_nwp_config(jg)%inwp_surface == LSS_TERRA ) THEN


        !
        ! 1. Consider land points (may have tiles)
        !
        ! - loop over surface tiles
        ! - note that different grid points may have different numbers 
        !   of active tiles (1<=ntiles<=ntiles_total). Therefore each tile has a 
        !   separate index list.
        ! 
        DO jt = 1, ntiles_total

          i_count_lnd = ext_data%atm%gp_count_t(jb,jt)

          IF (i_count_lnd == 0) CYCLE ! skip loop if the index list for the given tile is empty

!$NEC ivdep
          DO ic = 1, i_count_lnd

            jc = ext_data%atm%idx_lst_t(ic,jb,jt)

            ist = ext_data%atm%soiltyp(jc,jb) ! water (ist=9) and sea ice (ist=10) included

            ! surface albedo including moisture correction
            prm_diag%albdif_t(jc,jb,jt) = csalb(ist)&
              &                         - rad_csalbw(ist)*lnd_prog%w_so_t(jc,1,jb,jt)

          ENDDO  ! ic


          ! Account for Snow cover and vegetation
          !---------------------------------------
          IF (ntiles_lnd > 1) THEN     ! tile approach


!$NEC ivdep
            DO ic = 1, i_count_lnd

              jc = ext_data%atm%idx_lst_t(ic,jb,jt)

              zsnow= 0.0_wp

              ! temperature-dependent minimum snow albedo over glaciers
              ! this is needed to prevent unrealistically low albedos over Antarctica
              IF (ext_data%atm%lc_class_t(jc,jb,jt) == ext_data%atm%i_lc_snow_ice) THEN
                t_fac = MIN(1._wp,0.1_wp*(tmelt-lnd_prog%t_snow_t(jc,jb,jt)))
                zminsnow_alb = (1._wp-t_fac)*csalb_snow_min + t_fac*csalb_snow
              ELSE
                zminsnow_alb = csalb_snow_min
              ENDIF
              !
              ! 1. Consider effects of aging on solar snow albedo
              !
              zsalb_snow = zminsnow_alb + lnd_diag%freshsnow_t(jc,jb,jt)*(csalb_snow_max-zminsnow_alb)

              ! special treatment for forests and artificial surfaces
              ! - aging of snow not considered
              ! - instead, snow albedo is limited to land-class specific value

              ! get land cover class
              ilu = ext_data%atm%lc_class_t(jc,jb,jt) 
              zsnow_alb = MIN(zsalb_snow, ABS(ext_data%atm%snowalb_lcc(ilu)))


              ! 2. account for plant cover and snow cover
              !
              ! plant cover
              zvege = ext_data%atm%plcov_t(jc,jb,jt)
              ! snow cover fraction
              zsnow = lnd_diag%snowfrac_t(jc,jb,jt)

              ! 3. compute final solar snow albedo
              !
              prm_diag%albdif_t(jc,jb,jt) = zsnow * zsnow_alb           &  ! snow covered
                &  + (1.0_wp - zsnow) * (zvege * csalb_p                &  ! snow-free with vege
                &  + (1.0_wp - zvege) * prm_diag%albdif_t(jc,jb,jt))       ! snow-free bare

            ENDDO  ! ic


          ELSE  ! without tiles

!$NEC ivdep
            DO ic = 1, i_count_lnd

              jc = ext_data%atm%idx_lst_t(ic,jb,jt)

              zsnow= 0.0_wp

              ! temperature-dependent minimum snow albedo over glaciers
              ! this is needed to prevent unrealistically low albedos over Antarctica
              IF (ext_data%atm%lc_class_t(jc,jb,jt) == ext_data%atm%i_lc_snow_ice) THEN
                t_fac = MIN(1._wp,0.1_wp*(tmelt-lnd_prog%t_snow_t(jc,jb,jt)))
                zminsnow_alb = (1._wp-t_fac)*csalb_snow_min + t_fac*csalb_snow
              ELSE
                zminsnow_alb = csalb_snow_min
              ENDIF
              !
              ! 1. Consider effects of aging on solar snow albedo
              !
              zsalb_snow = zminsnow_alb + lnd_diag%freshsnow_t(jc,jb,jt)*(csalb_snow_max-zminsnow_alb)

              ! special treatment for forests
              zsnow_alb = zsalb_snow*(1._wp-ext_data%atm%for_e(jc,jb)-ext_data%atm%for_d(jc,jb)) &
                + csalb_snow_fe * ext_data%atm%for_e(jc,jb)                       &
                + csalb_snow_fd * ext_data%atm%for_d(jc,jb)

              ! 2. account for plant cover and snow cover
              !
              ! plant cover
              zvege = ext_data%atm%plcov_t(jc,jb,jt)
              ! snow cover fraction
              zsnow = lnd_diag%snowfrac_t(jc,jb,jt)

              ! 3. compute final solar snow albedo
              !
              prm_diag%albdif_t(jc,jb,jt) = zsnow * zsnow_alb           &  ! snow covered
                &  + (1.0_wp - zsnow) * (zvege * csalb_p                &  ! snow-free with vege
                &  + (1.0_wp - zvege) * prm_diag%albdif_t(jc,jb,jt))       ! snow-free bare

            ENDDO  ! ic

          ENDIF ! ntiles_lnd

        ENDDO  !ntiles





        ! 2. Consider water points with/without seaice model
        !
        IF ( (lseaice) ) THEN  ! seaice model switched on

          !
          ! Open sea points
          !
          ! - loop over sea points (points which are at least partly ice-free) 
          !
          i_count_sea = ext_data%atm%list_seawtr%ncount(jb)

          DO ic = 1, i_count_sea
            jc = ext_data%atm%list_seawtr%idx(ic,jb)

            prm_diag%albdif_t(jc,jb,isub_water) = csalb(ist_seawtr)
          ENDDO

          ! whitecap albedo by breaking ocean waves

          IF (albedo_whitecap == 1) THEN
            DO ic = 1, i_count_sea
              jc = ext_data%atm%list_seawtr%idx(ic,jb)

              CALL sfc_albedo_whitecap(prm_diag%sp_10m(jc,jb), wc_fraction, wc_albedo)
  
              prm_diag%albdif_t(jc,jb,isub_water) = wc_fraction * wc_albedo + &
            & prm_diag%albdif_t(jc,jb,isub_water) * (1.0_wp - wc_fraction)

            ENDDO
          ENDIF

          !
          ! Sea-ice points
          !
          ! - loop over sea-ice points (points which are at least partly ice-covered)
          !
          i_count_seaice = ext_data%atm%list_seaice%ncount(jb)

          PrognosticSeaIceAlbedo: IF ( lprog_albsi ) THEN
            ! Use prognostic sea-ice albedo (computed within the routines of the sea-ice scheme)
            DO ic = 1, i_count_seaice
              jc = ext_data%atm%list_seaice%idx(ic,jb)
              prm_diag%albdif_t(jc,jb,isub_seaice) = MAX(0._wp,wtr_prog%alb_si(jc,jb))
              ! Note that missing values for alb_si are chosen to be negative (usually -1)
              ! If something goes wrong with the sea-ice index list, or if the initialization 
              ! fails, there is the danger that we make use negative albedo values. To force a
              ! model crash in this case, we set negative albedo values to 0.
              ! 
            ENDDO
          ELSE 
            ! Use diagnostic sea-ice albedo (computed here)
            DO ic = 1, i_count_seaice
              jc = ext_data%atm%list_seaice%idx(ic,jb)
              ! In case the sea ice model is used, compute ice albedo for seaice 
              ! points with an empirical formula taken from GME.
              ! The ice albedo is the lower the warmer, and therefore wetter 
              ! the ice is. Use ice temperature at time level nnew 
              ! (2-time level scheme in sea ice model).
              prm_diag%albdif_t(jc,jb,isub_seaice) = alb_seaice_equil( wtr_prog%t_ice(jc,jb) )
              ! gives alb_max = 0.70
              !       alb_min = 0.43
              ! compare with Mironov et. al (2012), Tellus
              !       alb_max = 0.65
              !       alb_min = 0.40
            ENDDO
          ENDIF PrognosticSeaIceAlbedo

        ELSE   ! no seaice model

          !
          ! Sea points
          !
          ! - loop over sea points (including sea-ice points)
          !
          i_count_sea = ext_data%atm%list_seawtr%ncount(jb)

          DO ic = 1, i_count_sea
            jc = ext_data%atm%list_seawtr%idx(ic,jb)

            ! special handling of sea ice points
            IF (lnd_prog%t_g_t(jc,jb,isub_water) < tf_salt) THEN 
              ist = ist_seaice  ! sea ice
            ELSE
              ist = ist_seawtr  ! sea water
            ENDIF

            prm_diag%albdif_t(jc,jb,isub_water) = csalb(ist)

          ENDDO

          ! whitecap albedo by breaking ocean waves

          IF (albedo_whitecap == 1) THEN
            DO ic = 1, i_count_sea
              jc = ext_data%atm%list_seawtr%idx(ic,jb)

              CALL sfc_albedo_whitecap(prm_diag%sp_10m(jc,jb), wc_fraction, wc_albedo)
  
              prm_diag%albdif_t(jc,jb,isub_water) = wc_fraction * wc_albedo + &
            & prm_diag%albdif_t(jc,jb,isub_water) * (1.0_wp - wc_fraction)

            ENDDO
          ENDIF

        ENDIF  ! seaice



        ! 3. Consider lake points with/without lake model
        !
        IF (llake) THEN  ! Lake model switched on

          !
          ! Lake points
          !
          ! - loop over lake points
          !
          i_count_flk = ext_data%atm%list_lake%ncount(jb)

!$NEC ivdep
          DO ic = 1, i_count_flk
            jc = ext_data%atm%list_lake%idx(ic,jb)

            !  In case ice is present at lake points, compute ice albedo 
            !  for lake points with an empirical formulation 
            !  proposed by Mironov and Ritter (2004) for use in GME 
            !  [ice_albedo=function(ice_surface_temperature)].
            !  Use surface temperature at time level "nnew".

            ! special handling for ice-covered lake points
            IF (wtr_prog%h_ice(jc,jb) > h_Ice_min_flk) THEN 

              prm_diag%albdif_t(jc,jb,isub_lake) = albedo_whiteice_ref                      &
                &              - (albedo_whiteice_ref - albedo_blueice_ref)                 &
                &              * EXP(-c_albice_MR*(tpl_T_f-lnd_prog%t_g_t(jc,jb,isub_lake)) &
                &              /tpl_T_f)
            ELSE
              prm_diag%albdif_t(jc,jb,isub_lake) = csalb(ist_seawtr)
            ENDIF
          ENDDO

        ELSE  ! Lake model switched off

          !
          ! Lake points
          !
          ! - loop over lake points
          !
          i_count_flk = ext_data%atm%list_lake%ncount(jb)

          DO ic = 1, i_count_flk
            jc = ext_data%atm%list_lake%idx(ic,jb)

            ! special handling of sea ice points
            IF (lnd_prog%t_g_t(jc,jb,isub_lake) < tpl_T_f) THEN 
              ist = ist_seaice ! sea ice
            ELSE
              ist = ist_seawtr   ! water
            ENDIF

            prm_diag%albdif_t(jc,jb,isub_lake) = csalb(ist)
          ENDDO

        ENDIF  ! llake





        !*****************************!
        !                             !
        !  Aggregate surface albedo   !
        !                             !
        !*****************************!

        ! Loop over ALL grid points


        !
        ! Aggregate surface albedo on all points
        !
        IF (ntiles_total == 1) THEN
 
          DO jc = i_startidx, i_endidx
            prm_diag%albdif(jc,jb) = prm_diag%albdif_t(jc,jb,1)
            ! For RRTM copy albdif to albnirdif and albvisdif 
            ! i.e. no distiction is made between shortwave, vis and nir albedo
            prm_diag%albvisdif(jc,jb) = prm_diag%albdif_t(jc,jb,1)
            prm_diag%albnirdif(jc,jb) = prm_diag%albdif_t(jc,jb,1)
          ENDDO

        ELSE ! aggregate fields over tiles

          prm_diag%albdif(i_startidx:i_endidx,jb)    = 0._wp

          DO jt = 1, ntiles_total+ntiles_water
            DO jc = i_startidx, i_endidx
              prm_diag%albdif(jc,jb) = prm_diag%albdif(jc,jb)      &
                &                    + prm_diag%albdif_t(jc,jb,jt) &
                &                    * ext_data%atm%frac_t(jc,jb,jt)
            ENDDO
          ENDDO

          DO jc = i_startidx, i_endidx
            ! For RRTM copy albdif to albnirdif and albvisdif 
            ! i.e. no distiction is made between shortwave, vis and nir albedo
            prm_diag%albvisdif(jc,jb) = prm_diag%albdif(jc,jb)
            prm_diag%albnirdif(jc,jb) = prm_diag%albdif(jc,jb)
          ENDDO
        ENDIF  ! ntiles_total = 1

      ELSE  ! surface model switched OFF


        DO jc = i_startidx, i_endidx
          
          ist = ist_seaice

          IF (ext_data%atm%llsm_atm_c(jc,jb) .OR. lnd_prog%t_g(jc,jb) >= tf_salt ) THEN
            ist = ext_data%atm%soiltyp(jc,jb) ! water (ist=9) and sea ice (ist=10) included
          ENDIF

          prm_diag%albdif(jc,jb) = csalb(ist)
          ! For RRTM copy albdif to albnirdif and albvisdif
          ! i.e. no distiction is made between shortwave, vis and nir albedo
          prm_diag%albvisdif(jc,jb) = csalb(ist)
          prm_diag%albnirdif(jc,jb) = csalb(ist)
          
        ENDDO

      ENDIF ! inwp_surface=1


      ! albvisdir, albnirdir only needed for RRTM 
      !
      ! compute black sky albedo from white sky albedo and solar zenith angle formula 
      ! as in Ritter-Geleyn's fesft. So far we do not distinguish between  
      ! visible and NIR spectral bands.
      DO jc = i_startidx, i_endidx

        prm_diag%albvisdir(jc,jb) = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb),prm_diag%albvisdif(jc,jb))

        ! no need to do the computation twice, since albvisdif=albnirdif=albdif
        ! Thus: just copy
        prm_diag%albnirdir(jc,jb) = prm_diag%albvisdir(jc,jb)
      ENDDO

    ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  END SUBROUTINE sfc_albedo



  !>
  !! Calculation of surface albedo based on MODIS data
  !!
  !! Calculation of surface albedo based on snow-free MODIS data. In addition, 
  !! snow/ice conditions are taken into account. The snow-free MODIS albedo is updated 
  !! on a daily basis.
  !! We distinguish between
  !! - shortwave broadband albedo  (diffuse, 0.3-5.0um): albdif
  !! - UV visible broadband albedo (diffuse, 0.3-0.7um): albvisdif
  !! - near IR broadband albedo    (diffuse, 0.7-5.0um): albnirdif
  !! - UV visible broadband albedo (direct , 0.3-0.7um): albvisdir
  !! - near IR broadband albedo    (direct , 0.7-5.0um): albnirdir
  !!
  !! albvisdif/albvisdir and albnirdif/albnirdir are exclusively used by the RRTM scheme
  !!
  !! Diffuse albedo:
  !!=================
  !! land points (snow-free): MODIS albedo is used with separate values for visible and 
  !!                          nera-IR spectral bands. No distinction for tiles
  !! land points (snow covered): based on land-class specific tabulated values with
  !!                             consideration of aging snow. No distinction between 
  !!                             visible and near-IR spectral bands, yet
  !! sea points: tabulated value. No distinction yet between visible and near-IR spectral 
  !!             bands.
  !! lake points: tabulated value. No distinction yet between visible and near-IR spectral 
  !!              bands.
  !! sea-ice points: based on an empirical formula taken from GME. No distinction yet 
  !!                 between visible and near-IR spectral bands.
  !! lake-ice points: based on an empirical formula taken from GME. No distinction yet 
  !!                 between visible and near-IR spectral bands.
  !!
  !!
  !! Direct albedo:
  !!===============
  !! missing documentation
  !!
  !! Possible improvements:
  !!=======================
  !! snow-covered land points: separate values for visible and near-IR spectral bands.
  !! sea-ice points: separate values for visible and near-IR spectral bands.
  !! lake-ice points: separate values for visible and near-IR spectral bands.
  !! sea/lake points: separate values for direct and diffuse radiation (see IFS)
  !!
  !! Note that when using separate values for visible and near-IR spectral bands, the 
  !! shortwave albedo albdif_t must be derived by spectral integration of the visible 
  !! and near-IR albedo.                         
  !! 
  SUBROUTINE sfc_albedo_modis(pt_patch, ext_data, lnd_prog, wtr_prog, lnd_diag, prm_diag, lacc)

    TYPE(t_patch),          INTENT(   in):: pt_patch  !< grid/patch info.

    TYPE(t_external_data),  INTENT(   in):: ext_data  !< external data

    TYPE(t_lnd_prog),       INTENT(   in):: lnd_prog  !< land prognostic state (new)

    TYPE(t_wtr_prog),       INTENT(   in):: wtr_prog  !< water prognostic state (new)

    TYPE(t_lnd_diag),       INTENT(   in):: lnd_diag  !< land diagnostic state

    TYPE(t_nwp_phy_diag),   INTENT(inout):: prm_diag

    LOGICAL,                INTENT(   in):: lacc      !< accelerator flag

    ! Local scalars:
    REAL(wp):: snow_frac               !< snow cover fraction
    REAL(wp):: zminsnow_alb            !< temperature-dependent minimum snow albedo
    REAL(wp):: zmaxsnow_alb            !< maximum snow albedo depending on landuse class
    REAL(wp):: zlimsnow_alb            !< upper limit snow albedo depending on snow depth and roughness length
    REAL(wp):: zsnowalb_lu             !< maximum snow albedo specified in landuse table
    REAL(wp):: zurb_isa                !< impervious surface area fraction of the urban canopy
    REAL(wp):: zsnowfree_albdif        !< shortwave broadband surface albedo (white sky)
    REAL(wp):: zsnowfree_albvisdif     !< UV visible broadband surface albedo (white sky)
    REAL(wp):: zsnowfree_albnirdif     !< near IR broadband surface albedo (white sky)
    REAL(wp):: t_fac                   !< factor for temperature dependency of zminsnow_alb over glaciers
    REAL(wp):: zsnowfrac(nproma)       !< aggregated snow-cover fraction
    REAL(wp):: wc_fraction             !< whitecap fraction
    REAL(wp):: wc_albedo               !< whitecap albedo

    REAL(wp):: zsnow_alb(nproma,ntiles_total) !< snow albedo

    ! Auxiliaries for tile-specific calculation of direct beam albedo
    REAL(wp):: zalbvisdir_t(nproma,ntiles_total+ntiles_water)
    REAL(wp):: zalbnirdir_t(nproma,ntiles_total+ntiles_water)

    LOGICAL :: lfrozenwater

    INTEGER :: jg                      !< patch ID
    INTEGER :: jb, jc, ic, jt          !< loop indices
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index
    INTEGER :: ist
    INTEGER :: ilu                     !< land cover class
    INTEGER :: i_count_lnd             !< number of land points
    INTEGER :: i_count_sea             !< number of sea points
    INTEGER :: i_count_flk             !< number of lake points
    INTEGER :: i_count_seaice          !< number of seaice points

    !-----------------------------------------------------------------------
    jg = pt_patch%id

    IF (atm_phy_nwp_config(jg)%inwp_surface == LSS_JSBACH) THEN
      ! Albedo is already up-to-date.
      RETURN
    END IF

    i_nchdom  = MAX(1,pt_patch%n_childdom)

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)
    

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jt,jc,ic,i_startidx,i_endidx,ist,snow_frac,t_fac,               &
!$OMP            zurb_isa,zsnowfree_albdif,zsnowfree_albvisdif,zsnowfree_albnirdif, &
!$OMP            zsnow_alb,ilu,i_count_lnd,i_count_sea,i_count_flk,                 &
!$OMP            i_count_seaice,zminsnow_alb,zmaxsnow_alb,zlimsnow_alb,zsnowalb_lu, &
!$OMP            zalbvisdir_t,zalbnirdir_t,zsnowfrac, wc_fraction, wc_albedo,       &
!$OMP            lfrozenwater)                                                      &
!$OMP            ICON_OMP_DEFAULT_SCHEDULE

    DO jb = i_startblk, i_endblk


      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        &                 i_startidx, i_endidx, rl_start, rl_end)


      zalbvisdir_t(:,:) = 0._wp
      zalbnirdir_t(:,:) = 0._wp

      !------------------------------------------------------------------------------
      ! Calculation of land surface albedo based on MODIS input data
      !------------------------------------------------------------------------------

      IF ( atm_phy_nwp_config(jg)%inwp_surface == LSS_TERRA ) THEN


        !
        ! 1. Consider land points (may have tiles)
        !
        ! - loop over surface tiles
        ! - note that different grid points may have different numbers 
        !   of active tiles (1<=ntiles<=ntiles_total). Therefore each tile has a 
        !   separate index list.
        ! 
        !$ACC DATA COPYIN(zalbvisdir_t, zalbnirdir_t) CREATE(zsnowfrac, zsnow_alb) IF(lacc)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
        !$ACC LOOP SEQ
        DO jt = 1, ntiles_total

          i_count_lnd = ext_data%atm%gp_count_t(jb,jt)

          IF (i_count_lnd == 0) CYCLE ! skip loop if the index list for the given tile is empty

!$NEC ivdep 
          !$ACC LOOP GANG VECTOR &
          !$ACC   PRIVATE(jc, snow_frac, t_fac, ilu) &
          !$ACC   PRIVATE(zminsnow_alb, zmaxsnow_alb, zlimsnow_alb, zsnowalb_lu) &
          !$ACC   PRIVATE(zurb_isa, zsnowfree_albdif, zsnowfree_albvisdif, zsnowfree_albnirdif)
          DO ic = 1, i_count_lnd

            jc = ext_data%atm%idx_lst_t(ic,jb,jt)

            ! current land-cover class
            ilu = ext_data%atm%lc_class_t(jc,jb,jt)

            ! maximum snow albedo specified in landuse table
            zsnowalb_lu = ABS(ext_data%atm%snowalb_lcc(ilu))

            IF (ilu == ext_data%atm%i_lc_snow_ice) THEN
              ! temperature-dependent minimum snow albedo over glaciers, ranging between 0.5 and 0.7
              ! this is needed to prevent unrealistically low albedos over Antarctica
              IF (itype_vegetation_cycle > 1) THEN ! use climatological T2M to avoid unphysical diurnal cycle of albedo
                t_fac = MIN(1._wp,0.1_wp*MAX(0._wp,tmelt-ext_data%atm%t2m_clim_hc(jc,jb)))
              ELSE
                t_fac = MIN(1._wp,0.1_wp*(tmelt-lnd_prog%t_snow_t(jc,jb,jt)))
              ENDIF
              zminsnow_alb = (1._wp-t_fac)*csalb_snow_min + t_fac*0.7_wp
            ELSE
              ! otherwise take minimum snow albedo as 60% of the landuse-class specific maximum snow albedo
              zminsnow_alb = MAX(0.4_wp*csalb_snow_min,MIN(csalb_snow_min,0.6_wp*zsnowalb_lu))
            ENDIF

            ! maximum snow albedo depending on landuse class if tile approach is used
            ! (without tiles, this would conflict with the snow albedo limit depending on the forest fraction)
            IF (ntiles_lnd > 1) THEN
              zmaxsnow_alb = MIN(csalb_snow_max,zsnowalb_lu)
              zlimsnow_alb = MIN(csalb_snow_max,SQRT(zsnowalb_lu))
            ELSE
              zmaxsnow_alb = csalb_snow_max
              zlimsnow_alb = csalb_snow_max
            ENDIF

            ! Further limitation of maximum snow albedo in case of thin snow cover depending on landuse 
            ! roughness length and SSO standard deviation; without tiles, it is assumed that non-forest
            ! vegetation is fully covered by snow if the snow depth exceeds three times the roughness length.
            ! The snow depth is taken to be at least 5 cm here because the effect of partial snow
            ! cover is considered below.
            !
            zmaxsnow_alb = MIN(zmaxsnow_alb, zlimsnow_alb *                                         &
              MIN(1._wp,SQRT(0.25_wp+0.25_wp*MAX(0.05_wp,lnd_diag%h_snow_t(jc,jb,jt))/              &
              MAX(MIN(0.5_wp,ext_data%atm%z0_lcc(ilu)),1.e-3_wp*ext_data%atm%sso_stdh_raw(jc,jb)) )))

            ! Consider effects of aging on solar snow albedo
            !
            zsnow_alb(jc,jt) = zminsnow_alb + lnd_diag%freshsnow_t(jc,jb,jt)*(zmaxsnow_alb-zminsnow_alb)

            ! Adaptive tuning of snow albedo if T2M is assimilated
            ! The upper limit is set to csalb_snow_max rather than zmaxsnow_alb in order to allow larger corrections
            ! in sparsely forested regions where the LU-based maximum albedo might be too low
            IF (icpl_da_snowalb >= 1) THEN
              zsnow_alb(jc,jt) = MAX(zminsnow_alb,MIN(zsnow_alb(jc,jt)*prm_diag%snowalb_fac(jc,jb),csalb_snow_max))
            ENDIF

            IF (ntiles_lnd == 1) THEN
              ! special treatment for forests
              ! - no landuse-class specific limitation of snow albedo
              ! - instead, snow albedo is limited as a function of the forest fraction
              zsnow_alb(jc,jt) = zsnow_alb(jc,jt)*(1._wp-ext_data%atm%for_e(jc,jb)-ext_data%atm%for_d(jc,jb))  &
                + csalb_snow_fe * ext_data%atm%for_e(jc,jb) + csalb_snow_fd * ext_data%atm%for_d(jc,jb)
            ENDIF

            ! snow cover fraction
            snow_frac = lnd_diag%snowfrac_t(jc,jb,jt)

            !
            ! TERRA_URB: Set urban albedo by modifying background albedo (of snow-free fraction)
            !
            IF (lterra_urb .AND. lurbalb) THEN
              ! impervious surface area fraction of the urban canopy
              zurb_isa = ext_data%atm%urb_isa_t(jc,jb,jt)

              zsnowfree_albdif    = zurb_isa * ext_data%atm%urb_alb_so_t(jc,jb,jt)      &
                                  ! * ext_data%atm%urb_alb_red_t(jc,jb,jt)  & Multiplication already made in mo_ext_data_init.f90
                &                 + (1._wp - zurb_isa) * ext_data%atm%alb_dif(jc,jb)

              zsnowfree_albvisdif = zurb_isa * ext_data%atm%urb_alb_so_t(jc,jb,jt)      &
                                  ! * ext_data%atm%urb_alb_red_t(jc,jb,jt)  & Multiplication already made in mo_ext_data_init.f90
                &                 + (1._wp - zurb_isa) * ext_data%atm%albuv_dif(jc,jb)

              zsnowfree_albnirdif = zurb_isa * ext_data%atm%urb_alb_so_t(jc,jb,jt)      &
                                  ! * ext_data%atm%urb_alb_red_t(jc,jb,jt)  & Multiplication already made in mo_ext_data_init.f90
                &                 + (1._wp - zurb_isa) * ext_data%atm%albni_dif(jc,jb)
            ELSE
              zsnowfree_albdif    = ext_data%atm%alb_dif(jc,jb)
              zsnowfree_albvisdif = ext_data%atm%albuv_dif(jc,jb)
              zsnowfree_albnirdif = ext_data%atm%albni_dif(jc,jb)
            END IF

            !
            ! Set albedo (background albedo plus modification due to snow)
            !
            ! shortwave broadband surface albedo (white sky)
            prm_diag%albdif_t(jc,jb,jt)    = snow_frac * zsnow_alb(jc,jt)               &
              &                            + (1._wp - snow_frac) * zsnowfree_albdif

            ! UV visible broadband surface albedo (white sky)
            prm_diag%albvisdif_t(jc,jb,jt) = snow_frac * zsnow_alb(jc,jt)               &
              &                            + (1._wp - snow_frac) * zsnowfree_albvisdif

            ! near IR broadband surface albedo (white sky)
            prm_diag%albnirdif_t(jc,jb,jt) = snow_frac * zsnow_alb(jc,jt)               &
              &                            + (1._wp - snow_frac) * zsnowfree_albnirdif

          ENDDO  ! ic

        ENDDO  !ntiles
        !$ACC END PARALLEL


        ! Albedo correction for artificially reduced snow-cover fractions (melting-rate parameterization):
        ! The albedo of the snow-free tile is adjusted such that the tile-averaged albedo equals that obtained
        ! without the artificial snow-cover reduction
        IF (ntiles_lnd+1 <= ntiles_total) THEN ! this IF is necessary to circumvent an ACC not-present bug with ntiles=1 (NVfortran 21.3)
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
          !$ACC LOOP SEQ
          DO jt = ntiles_lnd+1, ntiles_total

            i_count_lnd = ext_data%atm%gp_count_t(jb,jt)

            IF (i_count_lnd == 0) CYCLE ! skip loop if the index list for the given tile is empty

!$NEC ivdep
            !$ACC LOOP GANG VECTOR PRIVATE(jc)
            DO ic = 1, i_count_lnd

              jc = ext_data%atm%idx_lst_t(ic,jb,jt)

              IF (lnd_diag%snowfrac_lcu_t(jc,jb,jt) > lnd_diag%snowfrac_lc_t(jc,jb,jt) .AND. &
                  lnd_diag%snowfrac_lc_t(jc,jb,jt) > 0._wp) THEN

                prm_diag%albdif_t(jc,jb,jt-ntiles_lnd) = 1._wp / (1._wp-lnd_diag%snowfrac_lc_t(jc,jb,jt)) *     (    &
                  prm_diag%albdif_t(jc,jb,jt)*(lnd_diag%snowfrac_lcu_t(jc,jb,jt)-lnd_diag%snowfrac_lc_t(jc,jb,jt)) + &
                  prm_diag%albdif_t(jc,jb,jt-ntiles_lnd)*(1._wp-lnd_diag%snowfrac_lcu_t(jc,jb,jt))              )

                prm_diag%albnirdif_t(jc,jb,jt-ntiles_lnd) = 1._wp / (1._wp-lnd_diag%snowfrac_lc_t(jc,jb,jt)) *     (    &
                  prm_diag%albnirdif_t(jc,jb,jt)*(lnd_diag%snowfrac_lcu_t(jc,jb,jt)-lnd_diag%snowfrac_lc_t(jc,jb,jt)) + &
                  prm_diag%albnirdif_t(jc,jb,jt-ntiles_lnd)*(1._wp-lnd_diag%snowfrac_lcu_t(jc,jb,jt))              )

                prm_diag%albvisdif_t(jc,jb,jt-ntiles_lnd) = 1._wp / (1._wp-lnd_diag%snowfrac_lc_t(jc,jb,jt)) *     (    &
                  prm_diag%albvisdif_t(jc,jb,jt)*(lnd_diag%snowfrac_lcu_t(jc,jb,jt)-lnd_diag%snowfrac_lc_t(jc,jb,jt)) + &
                  prm_diag%albvisdif_t(jc,jb,jt-ntiles_lnd)*(1._wp-lnd_diag%snowfrac_lcu_t(jc,jb,jt))              )

              ENDIF

            ENDDO  ! ic

          ENDDO  !ntiles
          !$ACC END PARALLEL
        ENDIF


        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
        !$ACC LOOP SEQ
        DO jt = 1, ntiles_total

          i_count_lnd = ext_data%atm%gp_count_t(jb,jt)

          IF (i_count_lnd == 0) CYCLE ! skip loop if the index list for the given tile is empty

!$NEC ivdep
          !$ACC LOOP GANG VECTOR PRIVATE(jc, ilu, snow_frac)
          DO ic = 1, i_count_lnd

            jc = ext_data%atm%idx_lst_t(ic,jb,jt)

            ! current land-cover class
            ilu = ext_data%atm%lc_class_t(jc,jb,jt)

            ! snow cover fraction
            snow_frac = lnd_diag%snowfrac_t(jc,jb,jt)

            ! direct albedo (black sky) (vis and nir)
            !
            IF ( direct_albedo == 1 ) THEN       ! Ritter and Geleyn (1992)
              zalbvisdir_t(jc,jt) = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb), &
                &                                        prm_diag%albvisdif_t(jc,jb,jt))
              zalbnirdir_t(jc,jt) = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb), &
                &                                        prm_diag%albnirdif_t(jc,jb,jt)) 

            ELSE IF ( direct_albedo == 2 ) THEN  ! Zaengl
              zalbvisdir_t(jc,jt) = sfc_albedo_dir_zaengl (prm_diag%cosmu0(jc,jb),         &
                &                                          prm_diag%albvisdif_t(jc,jb,jt), &
                &                                          ext_data%atm%z0_lcc(ilu),       &
                &                                          ext_data%atm%sso_stdh_raw(jc,jb))

              zalbnirdir_t(jc,jt) = sfc_albedo_dir_zaengl (prm_diag%cosmu0(jc,jb),         &
                &                                          prm_diag%albnirdif_t(jc,jb,jt), &
                &                                          ext_data%atm%z0_lcc(ilu),       &
                &                                          ext_data%atm%sso_stdh_raw(jc,jb))

            ELSE IF ( direct_albedo == 3 ) THEN  ! Yang et al (2008)
              zalbvisdir_t(jc,jt) = snow_frac                                                  &
                &                 * sfc_albedo_dir_zaengl (prm_diag%cosmu0(jc,jb),             &
                &                                          zsnow_alb(jc,jt),                   &
                &                                          ext_data%atm%z0_lcc(ilu),           &
                &                                          ext_data%atm%sso_stdh_raw(jc,jb))   &
                &                 + (1._wp-snow_frac)                                          &
                &                 * sfc_albedo_dir_yang(prm_diag%cosmu0(jc,jb),                &
                &                                       ext_data%atm%albuv_dif(jc,jb))

              zalbnirdir_t(jc,jt) = snow_frac                                                  &
                &                 * sfc_albedo_dir_zaengl (prm_diag%cosmu0(jc,jb),             &
                &                                          zsnow_alb(jc,jt),                   &
                &                                          ext_data%atm%z0_lcc(ilu),           &
                &                                          ext_data%atm%sso_stdh_raw(jc,jb))   &
                &                 + (1._wp-snow_frac)                                          &
                &                 * sfc_albedo_dir_yang(prm_diag%cosmu0(jc,jb),                &
                &                                       ext_data%atm%albni_dif(jc,jb))

            ELSE IF ( direct_albedo == 4 ) THEN  ! Briegleb and Ramanathan (1992)
              zalbvisdir_t(jc,jt) = snow_frac                                                  &
                &                 * sfc_albedo_dir_zaengl (prm_diag%cosmu0(jc,jb),             &
                &                                          zsnow_alb(jc,jt),                   &
                &                                          ext_data%atm%z0_lcc(ilu),           &
                &                                          ext_data%atm%sso_stdh_raw(jc,jb))   &
                &                 + (1._wp-snow_frac)                                          &
                &                 * sfc_albedo_dir_briegleb(prm_diag%cosmu0(jc,jb),            &
                &                                       ext_data%atm%albuv_dif(jc,jb),         &
                &                                       ext_data%atm%z0_lcc(ilu))

              zalbnirdir_t(jc,jt) = snow_frac                                                  &
                &                 * sfc_albedo_dir_zaengl (prm_diag%cosmu0(jc,jb),             &
                &                                          zsnow_alb(jc,jt),                   &
                &                                          ext_data%atm%z0_lcc(ilu),           &
                &                                          ext_data%atm%sso_stdh_raw(jc,jb))   &
                &                 + (1._wp-snow_frac)                                          &
                &                 * sfc_albedo_dir_briegleb(prm_diag%cosmu0(jc,jb),            &
                &                                       ext_data%atm%albni_dif(jc,jb),         &
                &                                       ext_data%atm%z0_lcc(ilu))
            ENDIF


          ENDDO  ! ic

        ENDDO  !ntiles
        !$ACC END PARALLEL

        ! 2. Consider water points with/without seaice model
        !
        IF ( (lseaice) ) THEN  ! seaice model switched on

          !
          ! Open sea points
          !
          ! - loop over sea points  (points which are at least partly ice-free) 
          !
          i_count_sea = ext_data%atm%list_seawtr%ncount(jb)
!$NEC ivdep
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
          !$ACC LOOP GANG VECTOR PRIVATE(jc)
          DO ic = 1, i_count_sea
            jc = ext_data%atm%list_seawtr%idx(ic,jb)

            ! diffuse albedo

            prm_diag%albdif_t   (jc,jb,isub_water) = csalb(ist_seawtr)
            prm_diag%albvisdif_t(jc,jb,isub_water) = csalb(ist_seawtr)
            prm_diag%albnirdif_t(jc,jb,isub_water) = csalb(ist_seawtr)

            ! direct albedo

            SELECT CASE (direct_albedo_water)
              CASE (1)                      ! Ritter and Geleyn (1992)
                zalbvisdir_t(jc,isub_water) = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb), &
                  &                          prm_diag%albvisdif_t(jc,jb,isub_water))
                zalbnirdir_t(jc,isub_water) = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb), &
                  &                          prm_diag%albnirdif_t(jc,jb,isub_water))
              CASE (2,4)                    ! Yang et al (2008)
                zalbvisdir_t(jc,isub_water) = sfc_albedo_dir_yang(prm_diag%cosmu0(jc,jb), &
                  &                          prm_diag%albvisdif_t(jc,jb,isub_water))
                zalbnirdir_t(jc,isub_water) = sfc_albedo_dir_yang(prm_diag%cosmu0(jc,jb), &
                  &                          prm_diag%albnirdif_t(jc,jb,isub_water))
              CASE (3)                      ! Taylor et al (1996)
                zalbvisdir_t(jc,isub_water) = sfc_albedo_dir_taylor(prm_diag%cosmu0(jc,jb))
                zalbnirdir_t(jc,isub_water) = sfc_albedo_dir_taylor(prm_diag%cosmu0(jc,jb))
            END SELECT

          ENDDO
          !$ACC END PARALLEL

          ! whitecap albedo by breaking ocean waves

          IF (albedo_whitecap == 1) THEN
!$NEC ivdep
            !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
            !$ACC LOOP GANG VECTOR PRIVATE(jc, wc_fraction, wc_albedo)
            DO ic = 1, i_count_sea
              jc = ext_data%atm%list_seawtr%idx(ic,jb)

              CALL sfc_albedo_whitecap(prm_diag%sp_10m(jc,jb), wc_fraction, wc_albedo)
  
              prm_diag%albdif_t   (jc,jb,isub_water) = wc_fraction * wc_albedo + &
            & prm_diag%albdif_t   (jc,jb,isub_water) * (1.0_wp - wc_fraction)
              prm_diag%albvisdif_t(jc,jb,isub_water) = wc_fraction * wc_albedo + &
            & prm_diag%albvisdif_t(jc,jb,isub_water) * (1.0_wp - wc_fraction)
              prm_diag%albnirdif_t(jc,jb,isub_water) = wc_fraction * wc_albedo + &
            & prm_diag%albnirdif_t(jc,jb,isub_water) * (1.0_wp - wc_fraction)
              zalbvisdir_t           (jc,isub_water) = wc_fraction * wc_albedo + &
            & zalbvisdir_t           (jc,isub_water) * (1.0_wp - wc_fraction)
              zalbnirdir_t           (jc,isub_water) = wc_fraction * wc_albedo + &
            & zalbnirdir_t           (jc,isub_water) * (1.0_wp - wc_fraction)

            ENDDO
            !$ACC END PARALLEL
          ENDIF

          !
          ! Sea-ice points
          !
          ! - loop over sea-ice points (points which are at least partly ice-covered)
          !
          !
          i_count_seaice = ext_data%atm%list_seaice%ncount(jb)

          PrognosticSeaIceAlbedo_modis: IF ( lprog_albsi ) THEN 
            ! Use prognostic diffuse sea-ice albedo (computed within the routines of the sea-ice scheme)
!$NEC ivdep
            !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
            !$ACC LOOP GANG VECTOR PRIVATE(jc, t_fac)
            DO ic = 1, i_count_seaice
              jc = ext_data%atm%list_seaice%idx(ic,jb)

              prm_diag%albdif_t(jc,jb,isub_seaice) = MAX(0._wp,wtr_prog%alb_si(jc,jb))
              ! Note that missing values for alb_si are chosen to be negative (usually -1)
              ! If something goes wrong with the sea-ice index list, or if the initialization 
              ! fails, there is the danger that we make use negative albedo values. To force a
              ! model crash in this case, we set negative albedo values to 0.
              !
              IF (icpl_da_snowalb >= 2) THEN
                ! Adaptive tuning of sea-ice albedo analogous to snow albedo:
                ! the constant limits reflect the albedo bounds assumed in the sea-ice scheme;
                ! in addition, albedo reduction is limited to 10% to avoid strong albedo reduction if a cold bias occurs in the polar night
                prm_diag%albdif_t(jc,jb,isub_seaice) =                                                &
                  MAX(0.685_wp*csalb(ist_seaice), 0.9_wp*prm_diag%albdif_t(jc,jb,isub_seaice),        &
                  MIN(prm_diag%albdif_t(jc,jb,isub_seaice)*prm_diag%snowalb_fac(jc,jb),csalb_snow_max))
                ! Moreover, we reset the sea ice albedo to the original value close to the melting point
                ! in order to avoid potential unwanted impacts on sea ice melt during summer
                t_fac = MIN(1._wp, MAX(0._wp, 0.5_wp*(tmelt - wtr_prog%t_ice(jc,jb)) ))
                prm_diag%albdif_t(jc,jb,isub_seaice) = t_fac*prm_diag%albdif_t(jc,jb,isub_seaice) + &
                  (1._wp-t_fac)*wtr_prog%alb_si(jc,jb)
              ENDIF

            ENDDO
            !$ACC END PARALLEL
          ELSE 
            ! Use diagnostic diffuse sea-ice albedo (computed here)
!$NEC ivdep
            !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
            !$ACC LOOP GANG VECTOR PRIVATE(jc)
            DO ic = 1, i_count_seaice
              jc = ext_data%atm%list_seaice%idx(ic,jb)
              prm_diag%albdif_t(jc,jb,isub_seaice) = alb_seaice_equil( wtr_prog%t_ice(jc,jb) )
            ENDDO
            !$ACC END PARALLEL
          ENDIF PrognosticSeaIceAlbedo_modis
!$NEC ivdep
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
          !$ACC LOOP GANG VECTOR PRIVATE(jc)
          DO ic = 1, i_count_seaice
            jc = ext_data%atm%list_seaice%idx(ic,jb)

            ! diffuse albedo             
            prm_diag%albvisdif_t(jc,jb,isub_seaice) = prm_diag%albdif_t(jc,jb,isub_seaice)
            prm_diag%albnirdif_t(jc,jb,isub_seaice) = prm_diag%albdif_t(jc,jb,isub_seaice)

            ! direct albedo
            zalbvisdir_t(jc,isub_seaice) = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb),                &
              &                                              prm_diag%albvisdif_t(jc,jb,isub_seaice))
            zalbnirdir_t(jc,isub_seaice) = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb),                &
              &                                              prm_diag%albnirdif_t(jc,jb,isub_seaice))
          ENDDO
          !$ACC END PARALLEL

        ELSE

          !
          ! Consider sea points (including sea-ice points)
          !
          ! - loop over sea points
          !
          i_count_sea = ext_data%atm%list_seawtr%ncount(jb)
!$NEC ivdep
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
          !$ACC LOOP GANG VECTOR PRIVATE(jc, ist, lfrozenwater, wc_fraction, wc_albedo)
          DO ic = 1, i_count_sea
            jc = ext_data%atm%list_seawtr%idx(ic,jb)

            ! special handling of sea ice points
            IF (lnd_prog%t_g_t(jc,jb,isub_water) < tf_salt) THEN 
              lfrozenwater = .TRUE.
              ist = ist_seaice  ! sea ice
            ELSE
              lfrozenwater = .FALSE.
              ist = ist_seawtr  ! sea water
            ENDIF

            ! diffuse albedo
            prm_diag%albdif_t   (jc,jb,isub_water) = csalb(ist)
            prm_diag%albvisdif_t(jc,jb,isub_water) = csalb(ist)
            prm_diag%albnirdif_t(jc,jb,isub_water) = csalb(ist)

            ! direct albedo
            SELECT CASE (direct_albedo_water)
              CASE (1)                        ! Ritter and Geleyn (1992)
                zalbvisdir_t(jc,isub_water)   = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb), &
                  &                             prm_diag%albvisdif_t(jc,jb,isub_water))
                zalbnirdir_t(jc,isub_water)   = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb), &
                  &                             prm_diag%albnirdif_t(jc,jb,isub_water))
              CASE (2,4)                      ! Yang et al (2008)
                zalbvisdir_t(jc,isub_water)   = sfc_albedo_dir_yang(prm_diag%cosmu0(jc,jb), &
                  &                             prm_diag%albvisdif_t(jc,jb,isub_water))
                zalbnirdir_t(jc,isub_water)   = sfc_albedo_dir_yang(prm_diag%cosmu0(jc,jb), &
                  &                             prm_diag%albnirdif_t(jc,jb,isub_water))
              CASE (3)
                IF (lfrozenwater) THEN        ! Ritter and Geleyn (1992)
                  zalbvisdir_t(jc,isub_water) = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb), &
                    &                           prm_diag%albvisdif_t(jc,jb,isub_water))
                  zalbnirdir_t(jc,isub_water) = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb), &
                    &                           prm_diag%albnirdif_t(jc,jb,isub_water))
                ELSE                          ! Taylor et al (1996)
                  zalbvisdir_t(jc,isub_water) = sfc_albedo_dir_taylor(prm_diag%cosmu0(jc,jb))
                  zalbnirdir_t(jc,isub_water) = sfc_albedo_dir_taylor(prm_diag%cosmu0(jc,jb))
                ENDIF
            END SELECT

            ! whitecap albedo by breaking ocean waves
            IF ( albedo_whitecap == 1 .AND. .NOT. lfrozenwater ) THEN
              CALL sfc_albedo_whitecap(prm_diag%sp_10m(jc,jb), wc_fraction, wc_albedo)

              prm_diag%albdif_t   (jc,jb,isub_water) = wc_fraction * wc_albedo + &
            & prm_diag%albdif_t   (jc,jb,isub_water) * (1.0_wp - wc_fraction)
              prm_diag%albvisdif_t(jc,jb,isub_water) = wc_fraction * wc_albedo + &
            & prm_diag%albvisdif_t(jc,jb,isub_water) * (1.0_wp - wc_fraction)
              prm_diag%albnirdif_t(jc,jb,isub_water) = wc_fraction * wc_albedo + &
            & prm_diag%albnirdif_t(jc,jb,isub_water) * (1.0_wp - wc_fraction)
              zalbvisdir_t           (jc,isub_water) = wc_fraction * wc_albedo + &
            & zalbvisdir_t           (jc,isub_water) * (1.0_wp - wc_fraction)
              zalbnirdir_t           (jc,isub_water) = wc_fraction * wc_albedo + &
            & zalbnirdir_t           (jc,isub_water) * (1.0_wp - wc_fraction)
            ENDIF
          ENDDO
          !$ACC END PARALLEL

        ENDIF



        ! 3. Consider lake points with/without lake model
        !
        IF (llake) THEN

          !
          ! Lake points
          !
          ! - loop over lake points
          !
          i_count_flk = ext_data%atm%list_lake%ncount(jb)

!$NEC ivdep
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
          !$ACC LOOP GANG VECTOR PRIVATE(jc, lfrozenwater)
          DO ic = 1, i_count_flk
            jc = ext_data%atm%list_lake%idx(ic,jb)

            !  In case ice is present at lake points, compute ice albedo 
            !  for lake points with an empirical formulation 
            !  proposed by Mironov and Ritter (2004) for use in GME 
            !  [ice_albedo=function(ice_surface_temperature)].
            !  Use surface temperature at time level "nnew".

            ! special handling for ice-covered lake points
            !
            ! diffuse albedo 
            IF (wtr_prog%h_ice(jc,jb) > h_Ice_min_flk) THEN 
              lfrozenwater = .TRUE.
              prm_diag%albdif_t(jc,jb,isub_lake) = albedo_whiteice_ref                      &
                &              - (albedo_whiteice_ref - albedo_blueice_ref)                 &
                &              * EXP(-c_albice_MR*(tpl_T_f-lnd_prog%t_g_t(jc,jb,isub_lake)) &
                &              /tpl_T_f)
            ELSE
              lfrozenwater = .FALSE.
              prm_diag%albdif_t(jc,jb,isub_lake) = csalb(ist_seawtr)
            ENDIF

            prm_diag%albvisdif_t(jc,jb,isub_lake) = prm_diag%albdif_t(jc,jb,isub_lake)
            prm_diag%albnirdif_t(jc,jb,isub_lake) = prm_diag%albdif_t(jc,jb,isub_lake)

            ! direct albedo
            SELECT CASE (direct_albedo_water)
              CASE (1,4)                     ! Ritter and Geleyn (1992)
                zalbvisdir_t(jc,isub_lake)   = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb), &
                  &                            prm_diag%albvisdif_t(jc,jb,isub_lake))
                zalbnirdir_t(jc,isub_lake)   = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb), &
                  &                            prm_diag%albnirdif_t(jc,jb,isub_lake))
              CASE (2)                       ! Yang et al (2008)
                zalbvisdir_t(jc,isub_lake)   = sfc_albedo_dir_yang(prm_diag%cosmu0(jc,jb), &
                  &                            prm_diag%albvisdif_t(jc,jb,isub_lake))
                zalbnirdir_t(jc,isub_lake)   = sfc_albedo_dir_yang(prm_diag%cosmu0(jc,jb), &
                  &                            prm_diag%albnirdif_t(jc,jb,isub_lake))
              CASE (3)
                IF (lfrozenwater) THEN       ! Mironov and Ritter (2004) & Ritter and Geleyn (1992)
                  zalbvisdir_t(jc,isub_lake) = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb), &
                    &                          prm_diag%albvisdif_t(jc,jb,isub_lake))
                  zalbnirdir_t(jc,isub_lake) = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb), &
                    &                          prm_diag%albnirdif_t(jc,jb,isub_lake))
                ELSE                         ! Taylor et al (1996)
                  zalbvisdir_t(jc,isub_lake) = sfc_albedo_dir_taylor(prm_diag%cosmu0(jc,jb))
                  zalbnirdir_t(jc,isub_lake) = sfc_albedo_dir_taylor(prm_diag%cosmu0(jc,jb))
                ENDIF
            END SELECT

          ENDDO
          !$ACC END PARALLEL

        ELSE

          !
          ! Lake points
          !
          ! - loop over lake points
          !
          i_count_flk = ext_data%atm%list_lake%ncount(jb)
!$NEC ivdep
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
          !$ACC LOOP GANG VECTOR PRIVATE(jc, ist, lfrozenwater)
          DO ic = 1, i_count_flk
            jc = ext_data%atm%list_lake%idx(ic,jb)

            ! special handling of sea ice points
            IF (lnd_prog%t_g_t(jc,jb,isub_lake) < tpl_T_f) THEN 
              lfrozenwater = .TRUE.
              ist = ist_seaice  ! sea ice
            ELSE
              lfrozenwater = .FALSE.
              ist = ist_seawtr  ! water
            ENDIF

            ! diffuse albedo 
            prm_diag%albdif_t   (jc,jb,isub_lake) = csalb(ist)
            prm_diag%albvisdif_t(jc,jb,isub_lake) = csalb(ist)
            prm_diag%albnirdif_t(jc,jb,isub_lake) = csalb(ist)

            ! direct albedo
            SELECT CASE (direct_albedo_water)
              CASE (1,4)                     ! Ritter and Geleyn (1992)
                zalbvisdir_t(jc,isub_lake)   = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb), &
                  &                            prm_diag%albvisdif_t(jc,jb,isub_lake))
                zalbnirdir_t(jc,isub_lake)   = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb), &
                  &                            prm_diag%albnirdif_t(jc,jb,isub_lake))
              CASE (2)                       ! Yang et al (2008)
                zalbvisdir_t(jc,isub_lake)   = sfc_albedo_dir_yang(prm_diag%cosmu0(jc,jb), &
                  &                            prm_diag%albvisdif_t(jc,jb,isub_lake))
                zalbnirdir_t(jc,isub_lake)   = sfc_albedo_dir_yang(prm_diag%cosmu0(jc,jb), &
                  &                            prm_diag%albnirdif_t(jc,jb,isub_lake))
              CASE (3)
                IF (lfrozenwater) THEN       ! Ritter and Geleyn (1992)
                  zalbvisdir_t(jc,isub_lake) = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb), &
                    &                          prm_diag%albvisdif_t(jc,jb,isub_lake))
                  zalbnirdir_t(jc,isub_lake) = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb), &
                    &                          prm_diag%albnirdif_t(jc,jb,isub_lake))
                ELSE                         ! Taylor et al (1996)
                  zalbvisdir_t(jc,isub_lake) = sfc_albedo_dir_taylor(prm_diag%cosmu0(jc,jb))
                  zalbnirdir_t(jc,isub_lake) = sfc_albedo_dir_taylor(prm_diag%cosmu0(jc,jb))
                ENDIF
            END SELECT

          ENDDO
          !$ACC END PARALLEL


        ENDIF





        !*****************************!
        !                             !
        !  Aggregate surface albedo   !
        !                             !
        !*****************************!

        ! Loop over ALL grid points


        !
        ! Aggregate surface albedo on all points
        !
        IF (ntiles_total == 1) THEN

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx, i_endidx
            prm_diag%albdif(jc,jb) = prm_diag%albdif_t(jc,jb,1)
            ! albvisdif, albnirdif only needed for RRTM 
            prm_diag%albvisdif(jc,jb) = prm_diag%albvisdif_t(jc,jb,1)
            prm_diag%albnirdif(jc,jb) = prm_diag%albnirdif_t(jc,jb,1)

            ! albvisdir, albnirdir only needed for RRTM
            prm_diag%albvisdir(jc,jb) = zalbvisdir_t(jc,1)
            prm_diag%albnirdir(jc,jb) = zalbnirdir_t(jc,1)

            zsnowfrac(jc) = lnd_diag%snowfrac_t(jc,jb,1)
          ENDDO
          !$ACC END PARALLEL

        ELSE ! aggregate fields over tiles

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = i_startidx, i_endidx
            prm_diag%albdif   (jc,jb) = 0._wp
            prm_diag%albvisdif(jc,jb) = 0._wp
            prm_diag%albnirdif(jc,jb) = 0._wp
            prm_diag%albvisdir(jc,jb) = 0._wp
            prm_diag%albnirdir(jc,jb) = 0._wp
          ENDDO

          !$ACC LOOP SEQ
          DO jt = 1, ntiles_total+ntiles_water
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO jc = i_startidx, i_endidx

              prm_diag%albdif(jc,jb) = prm_diag%albdif(jc,jb)         &
                &                    + prm_diag%albdif_t(jc,jb,jt)    &
                &                    * ext_data%atm%frac_t(jc,jb,jt)

              ! albvisdif, albnirdif only needed for RRTM 
              prm_diag%albvisdif(jc,jb) = prm_diag%albvisdif(jc,jb)   &
                &                    + prm_diag%albvisdif_t(jc,jb,jt) &
                &                    * ext_data%atm%frac_t(jc,jb,jt)

              prm_diag%albnirdif(jc,jb) = prm_diag%albnirdif(jc,jb)   &
                &                    + prm_diag%albnirdif_t(jc,jb,jt) &
                &                    * ext_data%atm%frac_t(jc,jb,jt)

              ! albvisdir, albnirdir only needed for RRTM 
              prm_diag%albvisdir(jc,jb) = prm_diag%albvisdir(jc,jb)   &
                &                       + zalbvisdir_t(jc,jt) * ext_data%atm%frac_t(jc,jb,jt)

              prm_diag%albnirdir(jc,jb) = prm_diag%albnirdir(jc,jb)   &
                &                       + zalbnirdir_t(jc,jt) * ext_data%atm%frac_t(jc,jb,jt)

            ENDDO
          ENDDO

          ! aggregated snow-cover fraction for LW emissivity
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = i_startidx, i_endidx
            zsnowfrac(jc) = 0._wp
          ENDDO

          !$ACC LOOP SEQ
          DO jt = 1, ntiles_total
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO jc = i_startidx, i_endidx
              IF (ext_data%atm%fr_land(jc,jb) > 0._wp) THEN
              zsnowfrac(jc) = zsnowfrac(jc) + ext_data%atm%frac_t(jc,jb,jt) * &
                              lnd_diag%snowfrac_t(jc,jb,jt)/ext_data%atm%fr_land(jc,jb)
              ENDIF
            ENDDO
          ENDDO
          !$ACC END PARALLEL

        ENDIF  ! ntiles_total = 1

        ! Account for snow effect on LW emissivity, using values consistent with the CAMEL climatology
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          prm_diag%lw_emiss(jc,jb) = (1._wp-zsnowfrac(jc))*ext_data%atm%emis_rad(jc,jb) + 0.978_wp*zsnowfrac(jc)
          ! analogously use lower values for sea ice than for open water
          IF (lnd_diag%fr_seaice(jc,jb) > 0._wp) prm_diag%lw_emiss(jc,jb) = 0.975_wp
        ENDDO
        !$ACC END PARALLEL

        !$ACC WAIT
        !$ACC END DATA

      ELSE IF (atm_phy_nwp_config(jg)%inwp_surface == LSS_JSBACH) THEN

        ! Albedo is already up-to-date.

      ELSE  ! surface model switched OFF
        DO jc = i_startidx, i_endidx
          
          ist = ist_seaice

          IF (ext_data%atm%llsm_atm_c(jc,jb) .OR. lnd_prog%t_g(jc,jb) >= tf_salt ) THEN
            ist = ext_data%atm%soiltyp(jc,jb) ! water (ist=9) and sea ice (ist=10) included
          ENDIF

          prm_diag%albdif(jc,jb) = csalb(ist)
          ! albvisdif, albnirdif only needed for RRTM 
          prm_diag%albvisdif(jc,jb) = csalb(ist)
          prm_diag%albnirdif(jc,jb) = csalb(ist)

          ! albvisdir, albnirdir only needed for RRTM 
          prm_diag%albvisdir(jc,jb) = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb),  &
              &                                         prm_diag%albvisdif(jc,jb))
          prm_diag%albnirdir(jc,jb) = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb),  &
              &                                         prm_diag%albnirdif(jc,jb))
        ENDDO

      ENDIF ! inwp_surface=1

    ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  END SUBROUTINE sfc_albedo_modis


  !>
  !! Surface albedo set to particular value (for single column model)
  !!
  !! We distinguish between
  !! - shortwave broadband albedo  (diffuse, 0.3-5.0um): albdif
  !! - UV visible broadband albedo (diffuse, 0.3-0.7um): albvisdif
  !! - near IR broadband albedo    (diffuse, 0.7-5.0um): albnirdif
  !! - UV visible broadband albedo (direct , 0.3-0.7um): albvisdir
  !! - near IR broadband albedo    (direct , 0.7-5.0um): albnirdir
  !!
  !! albvisdif/albvisdir and albnirdif/albnirdir are exclusively used by the RRTM scheme
  !! 
  SUBROUTINE sfc_albedo_scm(pt_patch, albedo_fixed, prm_diag, lacc)

    TYPE(t_patch),          INTENT(   in):: pt_patch     !< grid/patch info.

    REAL(wp),               INTENT(   in):: albedo_fixed ! surface albedo value that is used fpor albedo_type ==3
                                                         ! (for single column model)
    LOGICAL,                INTENT(   in):: lacc         !< accelerator flag
    TYPE(t_nwp_phy_diag),   INTENT(inout):: prm_diag

    INTEGER :: jb, jc                  !< loop indices
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index

    !-----------------------------------------------------------------------

    i_nchdom  = MAX(1,pt_patch%n_childdom)

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx)      &
!$OMP            ICON_OMP_DEFAULT_SCHEDULE

    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        prm_diag%albdif   (jc,jb) = albedo_fixed
        prm_diag%albvisdif(jc,jb) = albedo_fixed
        prm_diag%albnirdif(jc,jb) = albedo_fixed
        prm_diag%albvisdir(jc,jb) = albedo_fixed
        prm_diag%albnirdir(jc,jb) = albedo_fixed
      END DO
      !$ACC END PARALLEL

    ENDDO

!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE sfc_albedo_scm


  !>
  !! Surface albedo for direct beam
  !!
  !! Surface albedo for direct beam, according to Ritter-Geleyn 
  !! radiation scheme.
  !!
  !!
  FUNCTION sfc_albedo_dir_rg (cosmu0, alb_dif)  RESULT (alb_dir)
    !
    REAL(wp), INTENT(IN) :: cosmu0           !< cosine of solar zenith angle (SZA)

    REAL(wp), INTENT(IN) :: alb_dif          !< diffuse albedo (NIR or VIS or broadband)

    REAL(wp) :: alb_dir

  !------------------------------------------------------------------------------

     alb_dir = MIN(0.999_wp,( 1.0_wp                             &
       &     +  0.5_wp * (cosmu0 * (1.0_wp/alb_dif - 1.0_wp)))   &
       &     / (1.0_wp + (cosmu0 * (1.0_wp/alb_dif - 1.0_wp)))**2)

  END FUNCTION sfc_albedo_dir_rg


  !>
  !! Surface albedo for direct beam - water surface (ocean or lake)
  !! accoring to IFS (43r1)
  !!
  !! Reference: Taylor et al. (1996, QJRMS), based on aircraft data
  !!
  FUNCTION sfc_albedo_dir_taylor (cosmu0)  RESULT (alb_dir)

    REAL(wp), INTENT(IN) :: cosmu0           !< cosine of solar zenith angle (SZA)

    REAL(wp) :: alb_dir



  !------------------------------------------------------------------------------

    alb_dir = MAX( 1.E-10_wp,                                    &
      &  0.037_wp / (1.1_wp * MIN(MAX(cosmu0,0.0_wp),1.0_wp) ** 1.4_wp + 0.15_wp) )

  END FUNCTION sfc_albedo_dir_taylor


  !>
  !! Surface albedo (SA) for direct beam over land after Zaengl
  !!
  !! Surface albedo (SA) for direct beam over land after Zaengl.
  !! Builds upon the parameterization used in the 
  !! Ritter-Geleyn radiation scheme. However, the direct beam albedo 
  !! is limited to diffuse albedo for landuse classes with "rough" vegetation
  !! and in mountainous regions. Surfaces are regarded as 'rough', 
  !! if z0>= 15 cm or SSO standard deviation >= 150 m. 
  !!
  !!
  FUNCTION sfc_albedo_dir_zaengl (cosmu0, alb_dif, z0, sso_stdh)  RESULT (alb_dir)
    !
    REAL(wp), INTENT(IN) :: cosmu0           !< cosine of solar zenith angle (SZA)

    REAL(wp), INTENT(IN) :: alb_dif          !< diffuse albedo (NIR or VIS or broadband)

    REAL(wp), INTENT(IN) :: z0               !< surface roughness length

    REAL(wp), INTENT(IN) :: sso_stdh         !< Standard deviation of sub-grid scale orography

    REAL(wp) :: alb_dir

    REAL(wp) :: zalb_dir                     !< unlimited direct albedo
    REAL(wp) :: zalbdirlim                   !< limit, which the computed albedo must not exceed
    REAL(wp) :: zalbdirfac                   !< factor for limit computation




  !------------------------------------------------------------------------------

    ! unlimited direct beam albedo according to Ritter-Geleyn's formulation
    zalb_dir = sfc_albedo_dir_rg(cosmu0, alb_dif)


    ! Full limitation is applied for a roughness length >= 15 cm or an SSO standard deviation >= 150 m
    zalbdirfac = MAX(0.01_wp*(sso_stdh - 50._wp), 10._wp *(z0 - 0.05_wp) )
    zalbdirfac = MIN(1._wp, MAX(0._wp, zalbdirfac))

    zalbdirlim = zalbdirfac*alb_dif + (1._wp-zalbdirfac)*zalb_dir
    alb_dir     = MIN(zalb_dir, zalbdirlim)

  END FUNCTION sfc_albedo_dir_zaengl



  !>
  !! Surface albedo for direct beam
  !!
  !! Surface albedo for direct beam, according to Yang (2008).
  !!
  !!
  !! Literature
  !! - Yang, F. et al. (2008), Dependence of Land Surface Albedo on Solar Zenith Angle:
  !!   Observations and Model Parameterization. J. App. Meteorology and Climatology, 47, 
  !!   2963-2982
  !!
  FUNCTION sfc_albedo_dir_yang (cosmu0, alb_dif)  RESULT (alb_dir)
    !
    REAL(wp), INTENT(IN) :: cosmu0           !< cosine of solar zenith angle (SZA)

    REAL(wp), INTENT(IN) :: alb_dif          !< diffuse albedo (NIR or VIS or broadband)

    REAL(wp) :: alb_dir




  !------------------------------------------------------------------------------

     alb_dir = MIN(0.999_wp, alb_dif*( 1.0_wp + 1.14_wp)/(1._wp + 1.48*cosmu0) )

  END FUNCTION sfc_albedo_dir_yang


  !>
  !! Surface albedo for direct beam
  !!
  !! Surface albedo for direct beam, according to Briegleb (1992).
  !!
  !!
  !! Literature
  !! - Yang, F. et al. (2008), Dependence of Land Surface Albedo on Solar Zenith Angle:
  !!   Observations and Model Parameterization. J. App. Meteorology and Climatology, 47, 
  !!   2963-2982
  !! - Briegleb, B. and V. Ramanathan (1982), Spectral and Diurnal Variations in Clear Sky 
  !!   Planetary Albedo, J. App. Met., 21, 1160-1171
  !! - Briegleb, B. (1992), Delta-Eddington approximation for solar radiation 
  !!   in the NCAR Community Climate Model. J. Geophys. Res., 97, 7603-7612
  !!
  FUNCTION sfc_albedo_dir_briegleb (cosmu0, alb_dif, z0)  RESULT (alb_dir)
    !
    REAL(wp), INTENT(IN) :: cosmu0           !< cosine of solar zenith angle (SZA)

    REAL(wp), INTENT(IN) :: alb_dif          !< diffuse albedo (NIR or VIS or broadband)

    REAL(wp), INTENT(IN) :: z0               !< roughness length

    REAL(wp) :: d                !< tuning constant
                                 ! 0.4 (0.1) for vegetation classes where the albedo 
                                 ! has a strong (weak) dependence on solar zenith angle.
                                 ! weak dependence is assumed for 'rough' surfaces
                                 ! Therefore we use 0.4, if z0<=0.15 and 0.1 otherwise

    REAL(wp) :: alb_dir




  !------------------------------------------------------------------------------

     d       = MERGE(0.4_wp,0.1_wp,z0<=0.15_wp)
     alb_dir = MIN(0.999_wp, alb_dif*( 1.0_wp + d)/(1._wp + 2._wp*d*cosmu0) )

  END FUNCTION sfc_albedo_dir_briegleb


  !>
  !! Compute the fraction and albedo of whitecaps.
  SUBROUTINE sfc_albedo_whitecap (wsp_10m, wc_fraction, wc_albedo)

    !$ACC ROUTINE SEQ

    REAL(wp), INTENT(IN) :: wsp_10m !< 10m wind speed [m/s].
    REAL(wp), INTENT(OUT) :: wc_fraction !< Whitecap fraction as fraction of water surface [1].
    REAL(wp), INTENT(OUT) :: wc_albedo !< Whitecap albedo [1].

    !> Exponent for whitecap formation [1].
    REAL(wp), PARAMETER :: frac_exp = 1.59_wp
    !> Scaling factor for whitecap formation [(m/s)**-frac_exp].
    REAL(wp), PARAMETER :: frac_scale = 0.000397_wp
    !> Whitecap albedo [1].
    REAL(wp), PARAMETER :: albedo = 0.174_wp
    !> Maximum wind speed [m/s]. Value results in a maximum whitecap fraction of 4.6%.
    REAL(wp), PARAMETER :: max_wind = 20._wp

    wc_fraction = frac_scale * min(wsp_10m, max_wind) ** frac_exp
    wc_albedo = albedo

  END SUBROUTINE sfc_albedo_whitecap

END MODULE mo_albedo

