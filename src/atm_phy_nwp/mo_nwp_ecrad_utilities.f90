!
! The radiation scheme ecRad expects information from the host model (i.e., ICON) to be
!   copied to different ecRad data structures: t_ecrad_single_level_type, t_ecrad_gas_type,
!   t_ecrad_thermodynamics_type and t_ecrad_cloud_type. Similarly, the output of ecRad, i.e.
!   the radiative fluxes, are stored in a data structure t_ecrad_flux_type.
! This module offers subroutines that transfer the data from ICON variables into the
!   correct ecRad data structure. The subroutines are written in a way that they can be used
!   on the reduced radiation grid as well as on the full radiation grid. This ensures
!   consistency between the two modes.
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

MODULE mo_nwp_ecrad_utilities

  USE mo_kind,                   ONLY: wp
  USE mo_math_constants,         ONLY: rad2deg, pi
  USE mo_exception,              ONLY: finish
  USE mo_fortran_tools,          ONLY: set_acc_host_or_device, assert_acc_device_only
  USE mo_math_types,             ONLY: t_geographical_coordinates
  USE mo_atm_phy_nwp_config,     ONLY: atm_phy_nwp_config
  USE mo_physical_constants,     ONLY: rd, grav
  USE mo_radiation_config,       ONLY: vmr_co2, vmr_n2o, vmr_o2, vmr_ch4,        &
                                   &   vmr_cfc11, vmr_cfc12,                     &
                                   &   irad_h2o, irad_o3, irad_co2,              &
                                   &   irad_n2o, irad_ch4,                       &
                                   &   irad_o2, irad_cfc11, irad_cfc12,          &
                                   &   vpp_ch4, vpp_n2o, decorr_pole, decorr_equator
  USE mo_nwp_tuning_config,      ONLY: tune_difrad_3dcont
  USE mtime,                     ONLY: datetime
  USE mo_bc_greenhouse_gases,    ONLY: ghg_co2mmr, ghg_ch4mmr, ghg_n2ommr, ghg_cfcmmr
  USE mo_parallel_config,        ONLY: nproma, nproma_sub
#ifdef __ECRAD
  USE mo_ecrad,                  ONLY: ecrad_set_gas_units,                      &
                                   &   t_ecrad_conf,                             &
                                   &   t_ecrad_single_level_type,                &
                                   &   t_ecrad_thermodynamics_type,              &
                                   &   t_ecrad_gas_type, t_ecrad_flux_type,      &
                                   &   t_ecrad_cloud_type,                       &
                                   &   IMassMixingRatio, IVolumeMixingRatio,     &
                                   &   ecRad_IH2O, ecRad_ICO2, ecRad_IO3,        &
                                   &   ecRad_IN2O, ecRad_ICH4,                   &
                                   &   ecRad_IO2, ecRad_ICFC11, ecRad_ICFC12,    &
                                   &   ecRad_IHCFC22, ecRad_ICCl4,               &
                                   &   nweight_nir_ecrad, iband_nir_ecrad,       &
                                   &   weight_nir_ecrad, nweight_vis_ecrad,      &
                                   &   iband_vis_ecrad, weight_vis_ecrad,        &
                                   &   nweight_par_ecrad, iband_par_ecrad,       &
                                   &   weight_par_ecrad,                         &
                                   &   ecrad_hyd_list,                           &   
                                   &   ecrad_iqc, ecrad_iqi, ecrad_iqr,          &
                                   &   ecrad_iqs, ecrad_iqg
#endif

  USE mo_exception,              ONLY: message
  USE mo_grid_config,            ONLY: l_scm_mode
  USE mo_scm_nml,                ONLY: lon_scm, lat_scm 

  IMPLICIT NONE

  PRIVATE
#ifdef __ECRAD
  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nwp_ecrad_utilities'


  PUBLIC :: ecrad_set_single_level
  PUBLIC :: ecrad_set_thermodynamics
  PUBLIC :: ecrad_set_clouds
  PUBLIC :: ecrad_set_gas
  PUBLIC :: ecrad_store_fluxes
  PUBLIC :: add_3D_diffuse_rad
  PUBLIC :: get_indices_rad_subblock

  ! helper functions to be removed once acc is merged into libecrad

CONTAINS


  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE ecrad_set_single_level:
  !! Set ecRad single level information, i.e. fill t_ecrad_single_level_type with information
  !! from ICON.
  !!
  SUBROUTINE ecrad_set_single_level(ecrad_single_level, current_datetime, cell_center, cosmu0, tsfc, &
    &                               albvisdif, albnirdif, albvisdir, albnirdir, emis_rad, i_startidx, i_endidx, &
    &                               lacc)

    TYPE(t_ecrad_single_level_type), INTENT(inout) :: &
      &  ecrad_single_level       !< ecRad single level information
    TYPE(datetime), POINTER, INTENT(in) :: &
      &  current_datetime         !< Current date and time
    TYPE(t_geographical_coordinates), INTENT(in), TARGET :: &
      &  cell_center(:)           !< lon/lat information of cell centers
    REAL(wp), INTENT(in)     :: &
      &  cosmu0(:),             & !< Cosine of solar zenith angle
      &  tsfc(:),               & !< Surface temperature
      &  albvisdif(:),          & !< Surface albedo for visible range (diffuse)
      &  albnirdif(:),          & !< Surface albedo for near IR range (diffuse)
      &  albvisdir(:),          & !< Surface albedo for visible range (direct)
      &  albnirdir(:),          & !< Surface albedo for near IR range (direct)
      &  emis_rad(:)              !< Longwave surface emissivity
    INTEGER, INTENT(in)      :: &
      &  i_startidx, i_endidx     !< Start and end index of nproma loop in current block
    LOGICAL, OPTIONAL, INTENT(IN) :: lacc
! Local variables
    INTEGER                  :: &
      &  jc,                    & !< loop index
      &  seed_in_time             !< temporal contribution to the seed

    TYPE(t_geographical_coordinates), TARGET, ALLOCATABLE :: scm_center(:)
    TYPE(t_geographical_coordinates), POINTER             :: ptr_center(:)

    CALL assert_acc_device_only("ecrad_set_single_level", lacc)

    ! SCM: read lat/lon for horizontally uniform zenith angle
    IF ( l_scm_mode ) THEN
#ifdef _OPENACC
      CALL finish('ecrad_set_single_level','l_scm_mode not ported to gpu')
#endif
      ALLOCATE(scm_center(SIZE(cell_center)))
      DO jc = i_startidx, i_endidx
        scm_center(jc)%lat = lat_scm * pi/180.
        scm_center(jc)%lon = lon_scm * pi/180.
      ENDDO
      ptr_center => scm_center
    ELSE
      ptr_center => cell_center
    ENDIF

    seed_in_time = create_rdm_seed_in_time(current_datetime)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = i_startidx, i_endidx
        ecrad_single_level%cos_sza(jc)            = cosmu0(jc)
        ecrad_single_level%skin_temperature(jc)   = tsfc(jc)
        ecrad_single_level%sw_albedo(jc,1)        = albvisdif(jc)
        ecrad_single_level%sw_albedo(jc,2)        = albnirdif(jc)
        ecrad_single_level%sw_albedo_direct(jc,1) = albvisdir(jc)
        ecrad_single_level%sw_albedo_direct(jc,2) = albnirdir(jc)
        ecrad_single_level%lw_emissivity(jc,1)    = emis_rad(jc)
        ecrad_single_level%iseed(jc)              = create_rdm_seed(ptr_center(jc)%lon, &
          &                                                         ptr_center(jc)%lat, &
          &                                                         seed_in_time)
    ENDDO !jc
    !$ACC END PARALLEL

    IF ( l_scm_mode ) THEN
      DEALLOCATE(scm_center)
    ENDIF

  END SUBROUTINE ecrad_set_single_level
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE ecrad_set_thermodynamics:
  !! Set ecRad thermodynamics information, i.e. fill t_ecrad_thermodynamics_type with information
  !! from ICON.
  !!
  SUBROUTINE ecrad_set_thermodynamics(ecrad_thermodynamics, temp, pres, pres_ifc, &
    &                                 nlev, nlevp1, i_startidx, i_endidx, lacc)

    TYPE(t_ecrad_thermodynamics_type), INTENT(inout) :: &
      &  ecrad_thermodynamics     !< ecRad thermodynamics information
    REAL(wp), INTENT(in)     :: &
      &  temp(:,:),             & !< Full level temperature field
      &  pres(:,:),             & !< Full level pressure field
      &  pres_ifc(:,:)            !< Half level pressure field
    INTEGER, INTENT(in)      :: &
      &  nlev, nlevp1,          & !< Number of vertical full and half levels
      &  i_startidx, i_endidx     !< Start and end index of nproma loop in current block
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

! Local variables
    INTEGER                  :: &
      &  jc, jk                   !< loop indices

      CALL assert_acc_device_only('ecrad_set_thermodynamics', lacc)

      !$ACC DATA PRESENT(ecrad_thermodynamics, temp, pres, pres_ifc)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk=1,nlevp1
        DO jc = i_startidx, i_endidx
          ecrad_thermodynamics%pressure_hl(jc,jk)    = pres_ifc(jc,jk)
        ENDDO !jc
      ENDDO !jk
      !$ACC END PARALLEL

      ! Temperature at half levels is interpolated in the same way as in rrtm so far.
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk=2,nlev
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          ecrad_thermodynamics%temperature_hl(jc,jk) =                                             &
            &                  (temp(jc,jk-1) * pres(jc,jk-1)  * ( pres(jc,jk) - pres_ifc(jc,jk) ) &
            &                + temp(jc,jk) * pres(jc,jk) * ( pres_ifc(jc,jk) - pres(jc,jk-1)))     &
            &                / ( pres_ifc(jc,jk) * (pres(jc,jk) - pres(jc,jk-1) ) )
        ENDDO !jc
      ENDDO !jk
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        ecrad_thermodynamics%temperature_hl(jc,nlevp1) = temp(jc,nlev) + (pres_ifc(jc,nlevp1) - pres(jc,nlev)) * &
                               (temp(jc,nlev-1) - temp(jc,nlev))/(pres(jc,nlev-1) - pres(jc,nlev))
        ecrad_thermodynamics%temperature_hl(jc,1)         = temp(jc,1)                        &
          &                   + ( pres_ifc(jc,1) - pres(jc,1) )                               &
          &                   * (temp(jc,1)      - ecrad_thermodynamics%temperature_hl(jc,2)) &
          &                   / (pres(jc,1)      - pres_ifc(jc,2) )
      ENDDO !jc
      !$ACC END PARALLEL

      ! Directly provide full level temperature and pressure to rrtm gas_optics in ecrad (see rrtm_pass_temppres_fl).
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          ecrad_thermodynamics%pressure_fl(jc,jk)    = pres(jc,jk)
          ecrad_thermodynamics%temperature_fl(jc,jk) = temp(jc,jk)
        ENDDO !jc
      ENDDO !jk
      !$ACC END PARALLEL

      CALL ecrad_thermodynamics%calc_saturation_wrt_liquid(istartcol=i_startidx, iendcol=i_endidx)

      !$ACC END DATA

  END SUBROUTINE ecrad_set_thermodynamics
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE ecrad_set_clouds:
  !! Set ecRad clouds information, i.e. fill t_ecrad_cloud_type with information
  !! from ICON.
  !!
  SUBROUTINE ecrad_set_clouds(ecrad_cloud, ecrad_thermodynamics, qc, qi, clc, temp, pres, acdnc,                  &
    &                         fr_glac, fr_land, qr,qs,qg,reff_liq, reff_frz, reff_rain, reff_snow, reff_graupel,  &
    &                         icpl_reff, fact_reffc, clc_min, use_general_cloud_optics, cell_center,              &
    &                         nlev, i_startidx, i_endidx, lacc)

    TYPE(t_ecrad_cloud_type), INTENT(inout) :: &
      &  ecrad_cloud              !< ecRad cloud information
    TYPE(t_ecrad_thermodynamics_type), INTENT(inout) :: &
      &  ecrad_thermodynamics     !< ecRad thermodynamics information
    TYPE(t_geographical_coordinates), INTENT(in), TARGET :: &
      &  cell_center(:)           !< lon/lat information of cell centers 
    REAL(wp), TARGET, INTENT(in) :: &
      &  qc(:,:),               & !< Total cloud water (gridscale + subgridscale)
      &  qi(:,:)                  !< Total cloud ice   (gridscale + subgridscale)
    REAL(wp), INTENT(in)     :: &
      &  clc(:,:),              & !< Cloud cover
      &  temp(:,:),             & !< Full level temperature field
      &  pres(:,:),             & !< Full level pressure field
      &  fact_reffc,            & !< Factor in the calculation of cloud droplet effective radius
      &  clc_min                  !< Minimum cloud cover value to be considered as partly cloudy
    REAL(wp), POINTER, INTENT(in)     :: &
      &  acdnc(:,:),            & !< Cloud droplet numb. conc. (m-3)
      &  fr_glac(:),            & !< fraction of land covered by glaciers
      &  fr_land(:),            & !< land-sea mask. (1. = land, 0. = sea/lakes)
      &  qr(:,:),               & !< rain
      &  qs(:,:),               & !< snow
      &  qg(:,:),               & !< graupel
      &  reff_liq(:,:),         & !< effective radius of the liquid phase (external)
      &  reff_frz(:,:),         & !< effective radius of the frozen phase (external)
      &  reff_rain(:,:),        & !< effective radius of the rain phase (external)
      &  reff_snow(:,:),        & !< effective radius of the snow phase (external)
      &  reff_graupel(:,:)        !< effective radius of the graupel phase (external)
    INTEGER, INTENT(in)      :: &
      &  icpl_reff,             & !< Option for effective radius
      &  nlev,                  & !< Number of vertical full levels
      &  i_startidx, i_endidx     !< Start and end index of nproma loop in current block
    LOGICAL, INTENT(in)      :: &
      &  use_general_cloud_optics !< Use general cloud optics
    LOGICAL, OPTIONAL, INTENT(IN) :: lacc
! Local variables
    REAL(wp), POINTER, CONTIGUOUS :: &
      &  ptr_qx(:,:),           & !< Generic pointer to mass concentratio
      &  ptr_reff_x(:,:)          !< Generic pointer to effective radius
    REAL(wp)                 :: &
      &  lwc, iwc,              & !< Cloud liquid and ice water content
      &  liwcfac,               & !< Factor to calculate cloud liquid and ice water content
      &  zcos_lat,              & !< latitude factor cos(lat)
      &  zdecorr(i_startidx:i_endidx), &!< decorrelation length scale
      &  reff_min, reff_max       !< Limits for reff needed by ecrad (hardcoded)
    REAL (wp), PARAMETER     :: &
      &  qcrit_rad = 5E-5_wp      !< Limit for when to consider large hydrometeors for radiation
    INTEGER                  :: &
      &  ntypes_cloud_opt,      & !< Number of hydrometeors in generalized cloud optics
      &  jc, jk, iqx              !< loop indices
    INTEGER                  :: &
      &  iqc_loc, iqi_loc         !< indices for water and ice in ecrad_hyd_list
    LOGICAL                  :: &
      &  l_large_hyd              !< large hydrometeors change cloud fraction / small have max/min limits 

    TYPE(t_geographical_coordinates), TARGET, ALLOCATABLE :: scm_center(:)
    TYPE(t_geographical_coordinates), POINTER             :: ptr_center(:)

    NULLIFY(ptr_qx,ptr_reff_x)

    ! SCM: read lat for latitude-dependent decorrelation length scale, in radians
    IF ( l_scm_mode ) THEN
      ALLOCATE(scm_center(SIZE(cell_center)))
      DO jc = i_startidx, i_endidx
        scm_center(jc)%lat = lat_scm * pi/180.
        scm_center(jc)%lon = lon_scm * pi/180.
      ENDDO
      ptr_center => scm_center
    ELSE
      ptr_center => cell_center
    ENDIF

    ! Calculate latitude-dependent decorrelation length scale based on
    ! Shonk et al. 2010, but using COS function to be smoother over equator
    ! as is implemented in IFS
    !$ACC DATA CREATE(zdecorr, ptr_center)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
    !$ACC LOOP GANG VECTOR PRIVATE(zcos_lat)
    DO jc = i_startidx, i_endidx
      ! Shonk et al. (2010) but smoothed over the equator
      zcos_lat     = COS(ptr_center(jc)%lat)
      zdecorr(jc)  = decorr_pole + (decorr_equator-decorr_pole) * zcos_lat*zcos_lat
    ENDDO 
    !$ACC END PARALLEL

    CALL assert_acc_device_only("ecrad_set_clouds", lacc)

    ! Use latitude-dependent array of decorrelation length-scale values
    CALL ecrad_cloud%set_overlap_param(ecrad_thermodynamics, zdecorr, istartcol=i_startidx, iendcol=i_endidx)

    !$ACC DATA PRESENT(ecrad_cloud, qc, qi, clc, temp, pres, acdnc, fr_land, fr_glac, reff_frz, reff_liq)
    !$ACC DATA PRESENT(ecrad_cloud%q_liq, ecrad_cloud%q_ice, ecrad_cloud%re_liq, ecrad_cloud%re_ice)

    ! Generalized Hydrometeors
    IF ( use_general_cloud_optics ) THEN

      ntypes_cloud_opt = SIZE(ecrad_hyd_list)

      ! Set cloud fraction
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          IF ( clc(jc,jk) > clc_min ) THEN
            ecrad_cloud%fraction(jc,jk) = clc(jc,jk)
          ELSE
            ecrad_cloud%fraction(jc,jk) = 0._wp
          ENDIF
        END DO
      END DO

      ! Use effective radius from the reff module
      IF ( icpl_reff > 0 ) THEN
        ! Set up hydrometeors
        DO iqx = 1,ntypes_cloud_opt

          SELECT CASE( ecrad_hyd_list(iqx) )
            CASE( ecrad_iqc )
              ptr_qx => qc
              ptr_reff_x => reff_liq
              reff_min = 2.0e-6_wp
              reff_max = 32.0e-6_wp
              l_large_hyd = .false.

            CASE( ecrad_iqi )
              ptr_qx => qi
              ptr_reff_x => reff_frz
              reff_min = 5.0e-6_wp
              reff_max = 99.0e-6_wp
              l_large_hyd = .false.

            CASE( ecrad_iqr )
              ptr_qx => qr
              ptr_reff_x => reff_rain
              l_large_hyd = .true.

            CASE( ecrad_iqs )
              ptr_qx => qs
              ptr_reff_x => reff_snow
              l_large_hyd = .true.

            CASE( ecrad_iqg )
              ptr_qx => qg
              ptr_reff_x => reff_graupel
              l_large_hyd = .true.           
          END SELECT
          
          IF ( l_large_hyd) THEN ! Large hydrometeors update
            DO jk = 1, nlev
              DO jc = i_startidx, i_endidx
                IF ( ptr_qx(jc,jk) > qcrit_rad ) THEN
                  ecrad_cloud%mixing_ratio(jc,jk,iqx) = ptr_qx(jc,jk)
                  ecrad_cloud%fraction(jc,jk) = 1.0    ! Set cloud cover to one if large hydrometors are present
                ELSE
                  ecrad_cloud%mixing_ratio(jc,jk,iqx) = 0.0_wp
                ENDIF
                ecrad_cloud%effective_radius(jc,jk,iqx) = ptr_reff_x(jc,jk) 
              END DO
            END DO
          ELSE ! Small hydrometeors update
            DO jk = 1, nlev
              DO jc = i_startidx, i_endidx
                ecrad_cloud%mixing_ratio(jc,jk,iqx) = ptr_qx(jc,jk)
                ecrad_cloud%effective_radius(jc,jk,iqx) = MAX(MIN(ptr_reff_x(jc,jk),reff_max),reff_min)
              END DO
            END DO
          END IF

        END DO

      ELSE ! Caclulation of effective radius for water and ice when it is not coupled to the parameterization.
        ! Find indices for water and ice
        iqc_loc = 0
        iqi_loc = 0
        DO iqx = 1,ntypes_cloud_opt
          SELECT CASE( ecrad_hyd_list(iqx) )
            CASE( ecrad_iqc )
              iqc_loc = iqx
            CASE( ecrad_iqi )
              iqi_loc = iqx
          END SELECT
        END DO

 
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            IF (clc(jc,jk) > clc_min ) THEN
              liwcfac = 1000.0_wp / clc(jc,jk) * pres(jc,jk) / temp(jc,jk) / rd
            ELSE
              liwcfac = 0._wp
            END IF

            ! Careful with acdnc input: A division is performed and it is not checked for 0 as the function used
            ! to create acdnc returns always positive values
            IF ( iqc_loc >  0 ) THEN
              ecrad_cloud%mixing_ratio(jc,jk,iqc_loc) = qc(jc,jk)
              lwc                         = qc(jc,jk) * liwcfac
              ecrad_cloud%effective_radius(jc,jk,iqc_loc)   = reff_droplet(lwc, acdnc(jc,jk), fr_land(jc), fr_glac(jc), fact_reffc)
            END IF
            IF ( iqi_loc >  0 ) THEN
              ecrad_cloud%mixing_ratio(jc,jk,iqi_loc) = qi(jc,jk)
              iwc                         = qi(jc,jk) * liwcfac
              ecrad_cloud%effective_radius(jc,jk,iqi_loc)   = reff_crystal(iwc)
            END IF

          END DO
        END DO
      END IF

    ELSE ! No generalized hydrometeors
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(liwcfac, lwc, iwc)
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          ecrad_cloud%q_liq(jc,jk) = qc(jc,jk)
          ecrad_cloud%q_ice(jc,jk) = qi(jc,jk)
          IF ( clc(jc,jk) > clc_min ) THEN
            liwcfac = 1000.0_wp / clc(jc,jk) * pres(jc,jk) / temp(jc,jk) / rd
            ecrad_cloud%fraction(jc,jk) = clc(jc,jk)
          ELSE
            liwcfac = 0._wp
            ecrad_cloud%fraction(jc,jk) = 0._wp
          ENDIF
          IF ( icpl_reff == 0 ) THEN ! No external calculationcof reff.
            lwc                         = qc(jc,jk) * liwcfac
            iwc                         = qi(jc,jk) * liwcfac
            ! Careful with acdnc input: A division is performed and it is not checked for 0 as the function used
            ! to create acdnc returns always positive values
            ecrad_cloud%re_liq(jc,jk)   = reff_droplet(lwc, acdnc(jc,jk), fr_land(jc), fr_glac(jc), fact_reffc)
            ecrad_cloud%re_ice(jc,jk)   = reff_crystal(iwc)
          END IF
        ENDDO
      ENDDO
      !$ACC END PARALLEL
      IF ( icpl_reff > 0 ) THEN
        IF (.NOT. ASSOCIATED(reff_liq) .OR. .NOT. ASSOCIATED(reff_frz)) &
          & CALL finish('ecrad_set_clouds','effective radius fields not associated')
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            ecrad_cloud%re_liq(jc,jk) = MAX(MIN(reff_liq(jc,jk),32.0e-6_wp),2.0e-6_wp)  
            ecrad_cloud%re_ice(jc,jk) = MAX(MIN(reff_frz(jc,jk),99.0e-6_wp),5.0e-6_wp) 
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDIF
    END IF
    !$ACC WAIT
    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA
  END SUBROUTINE ecrad_set_clouds
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE ecrad_set_gas:
  !! Set ecRad gas information, i.e. fill t_ecrad_gas_type with information
  !! from ICON namelist parameters. The parameters irad_xyz and vmr_xyz are used
  !! (xyz=h2o, o3, co2, o2, cfc11, cfc12, n20 and ch4). Not all options are available for
  !! each of the gases. irad_xyz can have the following values:
  !!   0         : Set the gas globally constant to a concentration of zero
  !!   1         : Use information from a tracer variable
  !!   2         : Set the concentration to a globally constant value taken from vmr_xyz
  !!   3         : Use vmr_xyz at the surface, tanh-decay with height
  !!   7/9/79/97 : Use climatologies (only implemented for ozone)
  !!   11        : Read ozone from SCM input instead of here.
  !! The finish calls in case default should never trigger as the values for irad_xyz
  !! were already checked in mo_nml_crosscheck.
  !! 
  !! Care has to be taken for the indices that are passed to ecRad:
  !! - Values outside i_startidx:i_endidx are not initialized and must not be accessed
  !! - This interval might be smaller than jcs:jce
  !! - There is a documentation of the indices in the header of get_indices_rad_subblock
  !!
  SUBROUTINE ecrad_set_gas(ecrad_gas, ecrad_conf, o3, qv, pres, i_startidx, i_endidx, nlev, lacc)

    TYPE(t_ecrad_gas_type), INTENT(inout) :: &
      &  ecrad_gas                !< ecRad cloud information
    TYPE(t_ecrad_conf), INTENT(in) :: &
      &  ecrad_conf               !< ecRad configuration type
    REAL(wp), INTENT(in)     :: &
      &  o3(:,:),               & !< Ozone concentration
      &  qv(:,:),               & !< Water vapor
      &  pres(:,:)                !< Full level pressure
    INTEGER, INTENT(in)      :: &
      &  i_startidx, i_endidx,  & !< Start and end index of nproma loop in current block
      &  nlev                     !< Number of vertical full levels
    LOGICAL, INTENT(IN), OPTIONAL :: &
      &  lacc ! If true, use openacc
    ! Local variables
    REAL(wp), ALLOCATABLE    :: &
      &  ch4(:,:),              & !< Methane volume mixing ratio
      &  n2o(:,:)                 !< N2O volume mixing ratio
    INTEGER                  :: ncol
    CHARACTER(len=*), PARAMETER :: &
      &  routine = modname//'::ecrad_set_gas'

    CALL assert_acc_device_only(routine, lacc)

    ncol = SIZE(pres, 1)

    ! Water Vapor
    SELECT CASE(irad_h2o)
      CASE(0) ! No water vapor
        CALL ecrad_gas%put_well_mixed(ecRad_IH2O, IVolumeMixingRatio, 0._wp,  istartcol=i_startidx, iendcol=i_endidx)
      CASE(1) ! Use values from diagnosed water vapor content
        CALL ecrad_gas%put(ecRad_IH2O, IMassMixingRatio, qv(i_startidx:i_endidx,:), istartcol=i_startidx)
      CASE DEFAULT
        CALL finish(routine, 'Current implementation only supports irad_h2o = 0, 1')
    END SELECT

    ! Ozone
    SELECT CASE(irad_o3)
      CASE(0) ! No Ozone
        CALL ecrad_gas%put_well_mixed(ecRad_IO3,IVolumeMixingRatio, 0._wp,  istartcol=i_startidx,iendcol=i_endidx)
      CASE(10) ! Use values from interactive ozone
        CALL ecrad_gas%put(ecRad_IO3, IMassMixingRatio, o3(i_startidx:i_endidx,:), istartcol=i_startidx)
      CASE(5,7,9,79,97) ! Use values from GEMS/MACC (different profiles)
                        ! or time dependent concentration from external file
        CALL ecrad_gas%put(ecRad_IO3,  IMassMixingRatio, o3(i_startidx:i_endidx,:), istartcol=i_startidx)
      CASE(11) ! Ozone is read from SCM input file
        CALL message('mo_nwp_ecrad_utilities: irad_o3=11', &
          &          'Ozone used for radiation is read from SCM input file')
      CASE DEFAULT
        CALL finish(routine, 'Current implementation only supports irad_o3 = 0, 5, 7, 9, 10, 79, 97, 11')
    END SELECT

    !CO2
    SELECT CASE(irad_co2)
      CASE(0) ! No CO2
        CALL ecrad_gas%put_well_mixed(ecRad_ICO2,IVolumeMixingRatio, 0._wp,    istartcol=i_startidx,iendcol=i_endidx)
      CASE(2) ! Constant value derived from namelist parameter vmr_co2
        CALL ecrad_gas%put_well_mixed(ecRad_ICO2,IVolumeMixingRatio, vmr_co2,  istartcol=i_startidx,iendcol=i_endidx)
      CASE(4) ! time dependent concentration from external file
        CALL ecrad_gas%put_well_mixed(ecRad_ICO2,IMassMixingRatio, ghg_co2mmr,  istartcol=i_startidx,iendcol=i_endidx)
      CASE DEFAULT
        CALL finish(routine, 'Current implementation only supports irad_co2 = 0, 2, 4')
    END SELECT

    !O2
    SELECT CASE(irad_o2)
      CASE(0) ! No O2
        CALL ecrad_gas%put_well_mixed(ecRad_IO2,IVolumeMixingRatio, 0._wp,   istartcol=i_startidx,iendcol=i_endidx)
      CASE(2) ! Constant value derived from namelist parameter vmr_o2
        ! O2 is hardcoded within rrtm gas solver of ecRad, option for ecRad_IO2 only in case of psrad gas solver
        ! We still put it in ecRad, because this bug should be fixed with the next ecRad release
        CALL ecrad_gas%put_well_mixed(ecRad_IO2,IVolumeMixingRatio, vmr_o2,  istartcol=i_startidx,iendcol=i_endidx)
      CASE DEFAULT
        CALL finish(routine, 'Current implementation only supports irad_o2 = 0, 2')
    END SELECT

    !CFC11
    SELECT CASE(irad_cfc11)
      CASE(0) ! No CFC11
        CALL ecrad_gas%put_well_mixed(ecRad_ICFC11,IVolumeMixingRatio, 0._wp,    istartcol=i_startidx,iendcol=i_endidx)
      CASE(2) ! Constant value derived from namelist parameter vmr_cfc11
        CALL ecrad_gas%put_well_mixed(ecRad_ICFC11,IVolumeMixingRatio, vmr_cfc11,istartcol=i_startidx,iendcol=i_endidx)
      CASE(4) ! time dependent concentration from external file
        CALL ecrad_gas%put_well_mixed(ecRad_ICFC11,IMassMixingRatio, ghg_cfcmmr(1),istartcol=i_startidx, &
          &                           iendcol=i_endidx)
      CASE DEFAULT
        CALL finish(routine, 'Current implementation only supports irad_cfc11 = 0, 2, 4')
    END SELECT

    !CFC12
    SELECT CASE(irad_cfc12)
      CASE(0) ! No CFC12
        CALL ecrad_gas%put_well_mixed(ecRad_ICFC12,IVolumeMixingRatio, 0._wp,    istartcol=i_startidx,iendcol=i_endidx)
      CASE(2) ! Constant value derived from namelist parameter vmr_cfc12
        CALL ecrad_gas%put_well_mixed(ecRad_ICFC12,IVolumeMixingRatio, vmr_cfc12,istartcol=i_startidx,iendcol=i_endidx)
      CASE(4) ! time dependent concentration from external file
        CALL ecrad_gas%put_well_mixed(ecRad_ICFC12,IMassMixingRatio, ghg_cfcmmr(2),istartcol=i_startidx, &
          &                           iendcol=i_endidx)
      CASE DEFAULT
        CALL finish(routine, 'Current implementation only supports irad_cfc12 = 0, 2, 4')
    END SELECT

    !N2O
    SELECT CASE(irad_n2o)
      CASE(0) ! No N2O
        CALL ecrad_gas%put_well_mixed(ecRad_IN2O,IVolumeMixingRatio, 0._wp,  istartcol=i_startidx,iendcol=i_endidx)
      CASE(2) ! Constant value derived fromecrad_set_gas namelist parameter vmr_n2o
        CALL ecrad_gas%put_well_mixed(ecRad_IN2O,IVolumeMixingRatio, vmr_n2o,istartcol=i_startidx,iendcol=i_endidx)
      CASE(3) ! Tanh profile
        ALLOCATE(n2o(ncol,nlev))
        !$ACC ENTER DATA CREATE(n2o) ASYNC(1)
        CALL gas_profile(vmr_n2o, pres, vpp_n2o, i_startidx, i_endidx, nlev, n2o(:,:), lacc=.TRUE.)
        CALL ecrad_gas%put(ecRad_IN2O,  IVolumeMixingRatio, n2o(i_startidx:i_endidx,:), istartcol=i_startidx)
        !$ACC WAIT
        !$ACC EXIT DATA DELETE(n2o)
        DEALLOCATE(n2o)
      CASE(4) ! time dependent concentration from external file
        CALL ecrad_gas%put_well_mixed(ecRad_IN2O,IMassMixingRatio, ghg_n2ommr,istartcol=i_startidx,iendcol=i_endidx)
      CASE DEFAULT
        CALL finish(routine, 'Current implementation only supports irad_n2o = 0, 2, 3, 4')
    END SELECT

    !CH4
    SELECT CASE(irad_ch4)
      CASE(0) ! No CH4
        CALL ecrad_gas%put_well_mixed(ecRad_ICH4,IVolumeMixingRatio, 0._wp,  istartcol=i_startidx,iendcol=i_endidx)
      CASE(2) ! Constant value derived from namelist parameter vmr_ch4
        CALL ecrad_gas%put_well_mixed(ecRad_ICH4,IVolumeMixingRatio, vmr_ch4,istartcol=i_startidx,iendcol=i_endidx)
      CASE(3) ! Tanh profile
        ALLOCATE(ch4(ncol,nlev))
        !$ACC ENTER DATA CREATE(ch4) ASYNC(1)
        CALL gas_profile(vmr_ch4, pres, vpp_ch4, i_startidx, i_endidx, nlev, ch4(:,:), lacc=.TRUE.)
        CALL ecrad_gas%put(ecRad_ICH4,  IVolumeMixingRatio, ch4(i_startidx:i_endidx,:), istartcol=i_startidx)
        !$ACC WAIT
        !$ACC EXIT DATA DELETE(ch4)
        DEALLOCATE(ch4)
      CASE(4) ! time dependent concentration from external file
        CALL ecrad_gas%put_well_mixed(ecRad_ICH4,IMassMixingRatio, ghg_ch4mmr,istartcol=i_startidx,iendcol=i_endidx)
      CASE DEFAULT
        CALL finish(routine, 'Current implementation only supports irad_ch4 = 0, 2, 3, 4')
    END SELECT

    ! The following gases are currently not filled from the ICON side. Although they are set to 0 inside ecrad, 
    ! they are set to 0 here for completeness. 
    CALL ecrad_gas%put_well_mixed(ecRad_IHCFC22,IVolumeMixingRatio, 0._wp, istartcol=i_startidx, iendcol=i_endidx)
    CALL ecrad_gas%put_well_mixed(ecRad_ICCl4,  IVolumeMixingRatio, 0._wp, istartcol=i_startidx, iendcol=i_endidx)

    CALL ecrad_set_gas_units(ecrad_conf, ecrad_gas)

  END SUBROUTINE ecrad_set_gas
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE ecrad_store_fluxes:
  !! Stores radiative fluxes calculated by ecRad from the ecRad type t_ecrad_flux_type
  !! in the corresponding ICON data structure.
  !!
  SUBROUTINE ecrad_store_fluxes(jg, ecrad_flux, cosmu0, trsolall, trsol_up_toa, trsol_up_sfc, trsol_nir_sfc,   &
    &                           trsol_vis_sfc, trsol_par_sfc, fr_nir_sfc_diff, fr_vis_sfc_diff,                &
    &                           fr_par_sfc_diff,trsol_dn_sfc_diff, trsolclr_sfc, lwflxall, lwflx_up_sfc_rs,    &
    &                           lwflxclr_sfc, lwflx_up    , lwflx_dn    , swflx_up    , swflx_dn,              &
    &                           lwflx_up_clr, lwflx_dn_clr, swflx_up_clr, swflx_dn_clr,                        &
    &                           cosmu0mask, zsct, i_startidx, i_endidx, nlevp1, lacc)

    INTEGER, INTENT(in)   :: &
      &  jg                       !< domain index
    TYPE(t_ecrad_flux_type), INTENT(inout) :: &
      &  ecrad_flux               !< ecRad cloud information

    REAL(wp), INTENT(inout)  :: &
      &  cosmu0(:),             & !< Cosine of solar zenith angle
      &  trsolall(:,:),         & !< solar transmissivity, all sky, net down
      &  trsol_up_toa(:),       & !< upward solar transmissivity at TOA
      &  trsol_up_sfc(:),       & !< upward solar transmissivity at surface
      &  trsol_nir_sfc(:),      & !< downward transmissivity for near-IR rad. at surface
      &  trsol_vis_sfc(:),      & !< downward transmissivity for visible rad. at surface
      &  trsol_par_sfc(:),      & !< downward transmissivity for photosynthetically active rad. at surface
      &  fr_nir_sfc_diff(:),    & !< diffuse fraction for downward near-IR rad. at surface
      &  fr_vis_sfc_diff(:),    & !< diffuse fraction for downward visible rad. at surface
      &  fr_par_sfc_diff(:),    & !< diffuse fraction for downward photosynthetically active rad. at surface
      &  trsol_dn_sfc_diff(:),  & !< downward diffuse solar transmissivity at surface
      &  trsolclr_sfc(:),       & !< clear-sky net transmissivity at surface
      &  lwflxall(:,:),         & !< terrestrial flux, all sky, net down
      &  lwflx_up_sfc_rs(:),    & !< longwave upward flux at surface
      &  lwflxclr_sfc(:),       & !< longwave clear-sky flux at surface
      &  lwflx_up(:,:),         & !< longwave  3D upward   flux            [W/m2]
      &  lwflx_dn(:,:),         & !< longwave  3D downward flux            [W/m2]
      &  swflx_up(:,:),         & !< shortwave 3D upward   flux            [W/m2]
      &  swflx_dn(:,:),         & !< shortwave 3D downward flux            [W/m2]
      &  lwflx_up_clr(:,:),     & !< longwave  3D upward   flux clear-sky  [W/m2]
      &  lwflx_dn_clr(:,:),     & !< longwave  3D downward flux clear-sky  [W/m2]
      &  swflx_up_clr(:,:),     & !< shortwave 3D upward   flux clear-sky  [W/m2]
      &  swflx_dn_clr(:,:)        !< shortwave 3D downward flux clear-sky  [W/m2]

    LOGICAL, INTENT(in)      :: &
      &  cosmu0mask(:)            !< Mask if cosmu0 > 0
    REAL(wp),                INTENT(in)    :: zsct        !< Time-dependent solar constant

    INTEGER, INTENT(in)      :: &
      &  i_startidx, i_endidx,  & !< Start and end index of nproma loop in current block
      &  nlevp1                   !< Number of vertical half levels

    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    ! Local Variables
    INTEGER                  :: &
      &  jband, jc, jk            !< Loop indices

      CALL assert_acc_device_only("ecrad_store_fluxes", lacc)

      !$ACC DATA PRESENT(ecrad_flux, cosmu0, trsolall, trsol_up_toa, trsol_up_sfc, trsol_nir_sfc) &
      !$ACC   PRESENT(trsol_vis_sfc, trsol_par_sfc, fr_nir_sfc_diff, fr_vis_sfc_diff) &
      !$ACC   PRESENT(fr_par_sfc_diff, trsol_dn_sfc_diff, trsolclr_sfc, lwflxall, lwflx_up_sfc_rs) &
      !$ACC   PRESENT(lwflxclr_sfc, lwflx_up, lwflx_dn, swflx_up, swflx_dn, lwflx_up_clr) &
      !$ACC   PRESENT(lwflx_dn_clr, swflx_up_clr, swflx_dn_clr, cosmu0mask)

      ! Initialize output fields
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk = 1, nlevp1
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1,SIZE(trsolall,1)
          trsolall(jc,jk)      = 0._wp
        ENDDO
      ENDDO
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = 1,SIZE(trsolall,1)
        trsol_up_toa(jc)      = 0._wp
        trsol_up_sfc(jc)      = 0._wp
        trsol_nir_sfc(jc)     = 0._wp
        trsol_vis_sfc(jc)     = 0._wp
        trsol_par_sfc(jc)     = 0._wp
        fr_nir_sfc_diff(jc)   = 0._wp
        fr_vis_sfc_diff(jc)   = 0._wp
        fr_par_sfc_diff(jc)   = 0._wp
        trsol_dn_sfc_diff(jc) = 0._wp
        trsolclr_sfc(jc)      = 0._wp
      ENDDO
      !$ACC END PARALLEL

      ! Store output of 3-D Fluxes
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlevp1
        DO jc = i_startidx, i_endidx
          IF ( cosmu0mask(jc) ) THEN 
            ! solar transmissivity, all sky, net down
            trsolall(jc,jk)     = (ecrad_flux%sw_dn(jc,jk)-ecrad_flux%sw_up(jc,jk))/cosmu0(jc)
          ENDIF
          ! terrestrial flux, all sky, net down
          lwflxall(jc,jk)       = ecrad_flux%lw_dn(jc,jk)-ecrad_flux%lw_up(jc,jk)
        ENDDO
      ENDDO
      !$ACC END PARALLEL

      IF (atm_phy_nwp_config(jg)%l_3d_rad_fluxes) THEN    
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1, nlevp1
          DO jc = i_startidx, i_endidx
            ! LW/SW, up/down, all/clear 3D fluxes
            lwflx_up    (jc,jk)   = ecrad_flux%lw_up(jc,jk)
            lwflx_dn    (jc,jk)   = ecrad_flux%lw_dn(jc,jk)
  
            swflx_up    (jc,jk)   = ecrad_flux%sw_up(jc,jk)       *zsct
            swflx_dn    (jc,jk)   = ecrad_flux%sw_dn(jc,jk)       *zsct
            lwflx_up_clr(jc,jk)   = ecrad_flux%lw_up_clear(jc,jk)
            lwflx_dn_clr(jc,jk)   = ecrad_flux%lw_dn_clear(jc,jk)
            swflx_up_clr(jc,jk)   = ecrad_flux%sw_up_clear(jc,jk) *zsct
            swflx_dn_clr(jc,jk)   = ecrad_flux%sw_dn_clear(jc,jk) *zsct
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      END IF

      ! Store output of 2-D Fluxes
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        lwflx_up_sfc_rs(jc) = ecrad_flux%lw_up(jc,nlevp1)
        lwflxclr_sfc(jc)    = ecrad_flux%lw_dn_clear(jc,nlevp1) - ecrad_flux%lw_up_clear(jc,nlevp1)

        IF ( cosmu0mask(jc) ) THEN
          trsol_up_toa(jc)  = ecrad_flux%sw_up(jc,1)/cosmu0(jc)
          trsol_up_sfc(jc)  = ecrad_flux%sw_up(jc,nlevp1)/cosmu0(jc)

          trsol_dn_sfc_diff(jc) = (ecrad_flux%sw_dn(jc,nlevp1) - ecrad_flux%sw_dn_direct(jc,nlevp1))      &
            &                     / cosmu0(jc)
          trsolclr_sfc(jc)      = (ecrad_flux%sw_dn_clear(jc,nlevp1) - ecrad_flux%sw_up_clear(jc,nlevp1)) &
            &                     / cosmu0(jc)
        ENDIF
      ENDDO
      !$ACC END PARALLEL

      ! Total transmissivities in nir, vis, and par bands.
      CALL transmissivity_in_band(weight_nir_ecrad(:), iband_nir_ecrad(:), &
          & ecrad_flux%sw_dn_surf_band(:,i_startidx:i_endidx), cosmu0(i_startidx:i_endidx), &
          & cosmu0mask(i_startidx:i_endidx), trsol_nir_sfc(i_startidx:i_endidx), nbands=nweight_nir_ecrad)
      CALL transmissivity_in_band(weight_vis_ecrad(:), iband_vis_ecrad(:), &
          & ecrad_flux%sw_dn_surf_band(:,i_startidx:i_endidx), cosmu0(i_startidx:i_endidx), &
          & cosmu0mask(i_startidx:i_endidx), trsol_vis_sfc(i_startidx:i_endidx), nbands=nweight_vis_ecrad)
      CALL transmissivity_in_band(weight_par_ecrad(:), iband_par_ecrad(:), &
          & ecrad_flux%sw_dn_surf_band(:,i_startidx:i_endidx), cosmu0(i_startidx:i_endidx), &
          & cosmu0mask(i_startidx:i_endidx), trsol_par_sfc(i_startidx:i_endidx), nbands=nweight_par_ecrad)

      ! Transmissivities for direct radiation. We want the diffuse fraction, so this gets subtracted from and
      ! divided by the total below.
      CALL transmissivity_in_band(weight_nir_ecrad(:), iband_nir_ecrad(:), &
          & ecrad_flux%sw_dn_direct_surf_band(:,i_startidx:i_endidx), cosmu0(i_startidx:i_endidx), &
          & cosmu0mask(i_startidx:i_endidx), fr_nir_sfc_diff(i_startidx:i_endidx), nbands=nweight_nir_ecrad)
      CALL transmissivity_in_band(weight_vis_ecrad(:), iband_vis_ecrad(:), &
          & ecrad_flux%sw_dn_direct_surf_band(:,i_startidx:i_endidx), cosmu0(i_startidx:i_endidx), &
          & cosmu0mask(i_startidx:i_endidx), fr_vis_sfc_diff(i_startidx:i_endidx), nbands=nweight_vis_ecrad)
      CALL transmissivity_in_band(weight_par_ecrad(:), iband_par_ecrad(:), &
          & ecrad_flux%sw_dn_direct_surf_band(:,i_startidx:i_endidx), cosmu0(i_startidx:i_endidx), &
          & cosmu0mask(i_startidx:i_endidx), fr_par_sfc_diff(i_startidx:i_endidx), nbands=nweight_par_ecrad)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        fr_nir_sfc_diff(jc) = 1._wp - MIN(1._wp, fr_nir_sfc_diff(jc) / MAX(trsol_nir_sfc(jc), EPSILON(1._wp)))
        fr_vis_sfc_diff(jc) = 1._wp - MIN(1._wp, fr_vis_sfc_diff(jc) / MAX(trsol_vis_sfc(jc), EPSILON(1._wp)))
        fr_par_sfc_diff(jc) = 1._wp - MIN(1._wp, fr_par_sfc_diff(jc) / MAX(trsol_par_sfc(jc), EPSILON(1._wp)))
      END DO
      !$ACC END PARALLEL

      !$ACC END DATA

  END SUBROUTINE ecrad_store_fluxes
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! FUNCTION transmissivity_in_band:
  !! Computes a total transmissivity given a set of band-based transmissivities with their
  !! accompanying weights. NOTE: `tr_wgt` must be zero-initialized.
  !!
  SUBROUTINE transmissivity_in_band (weights, bands, tr_band, cosmu0, mask, tr_wgt, nbands)

    REAL(wp), INTENT(IN) :: weights(:) !< Weight for each band index present in `band`.
    INTEGER, INTENT(IN) :: bands(:) !< List of band indices.
    REAL(wp), INTENT(IN) :: tr_band(:,:) !< Transmissivities for each band (nbands,ncells).
    REAL(wp), INTENT(IN) :: cosmu0(:) !< Cosine of the solar zenith angle.
    LOGICAL, INTENT(IN) :: mask(:) !< Set to .TRUE. for elements that should be processed.
    REAL(wp), INTENT(INOUT) :: tr_wgt(:) !< Total (band-weighted) transmissivity (ncells).

    INTEGER, INTENT(IN), OPTIONAL :: nbands !< Number of bands (default: SIZE(bands))

    INTEGER :: jband
    INTEGER :: jc
    INTEGER :: nb
    INTEGER :: ncells

    IF (PRESENT(nbands)) THEN
      nb = nbands
    ELSE
      nb = SIZE(bands)
    END IF

    ncells = SIZE(tr_band, 2)

    !$ACC DATA PRESENT(weights, bands, tr_band, cosmu0, mask, tr_wgt)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
      DO jband = 1, nb
      !$ACC LOOP GANG VECTOR
      DO jc = 1, ncells
        IF (mask(jc)) THEN
          tr_wgt(jc) = tr_wgt(jc) + &
              & weights(jband) * tr_band(bands(jband),jc) / cosmu0(jc)
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC END DATA

  END SUBROUTINE transmissivity_in_band

  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE add_3D_diffuse_rad:
  !! Adds 3D contribution to diffuse radiation by reflection of direct solar radiation on scattered low clouds
  !!
  SUBROUTINE add_3D_diffuse_rad( &
        & ecrad_flux, clc, pres, temp, cosmu0, fr_nir_sfc_diff, fr_vis_sfc_diff, fr_par_sfc_diff, &
        & trsol_dn_sfc_diff, i_startidx, i_endidx, nlev, lacc &
      )

    TYPE(t_ecrad_flux_type), INTENT(in) :: ecrad_flux !< ecRad cloud information

    REAL(wp), INTENT(in)  :: &
      &  cosmu0(:),             & !< Cosine of solar zenith angle
      &  clc(:,:),              & !< cloud cover fraction
      &  pres(:,:),             & !< pressure
      &  temp(:,:)                !< temperature

    REAL(wp), INTENT(inout)  :: fr_nir_sfc_diff(:) !< diffuse fraction of downward near-IR transmissivity at surface
    REAL(wp), INTENT(inout)  :: fr_vis_sfc_diff(:) !< diffuse fraction of downward visible transmissivity at surface
    REAL(wp), INTENT(inout)  :: fr_par_sfc_diff(:) !< diffuse fraction of downward PAR transmissivity at surface
    REAL(wp), INTENT(inout)  :: trsol_dn_sfc_diff(:) !< downward diffuse solar transmissivity at surface

    INTEGER, INTENT(in)      :: &
      &  i_startidx, i_endidx,  & !< Start and end index of nproma loop in current block
      &  nlev                     !< Number of vertical levels

    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    ! Local Variables
    INTEGER                   ::  jc, jk  !< Loop indices

    REAL(wp), PARAMETER :: zdecorr = 2000.0_wp, & ! decorrelation length scale for cloud overlap scheme
                           epsi    = 1.e-20_wp

    REAL(wp) :: zcloud(i_startidx:i_endidx), ccmax, ccran, deltaz, alpha
    REAL(wp) :: diff_frac_corr !< Correction to the diffuse fraction.

    CALL assert_acc_device_only("add_3D_diffuse_rad", lacc)

    !$ACC DATA PRESENT(ecrad_flux, pres, clc, temp, cosmu0, fr_nir_sfc_diff, fr_vis_sfc_diff, fr_par_sfc_diff) &
    !$ACC   PRESENT(trsol_dn_sfc_diff) CREATE(zcloud)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = i_startidx, i_endidx
      zcloud(jc)     = 0.0_wp
    ENDDO
    !$ACC END PARALLEL

    ! Calculate low-level cloud cover fraction
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO jk = 2, nlev
      !$ACC LOOP GANG VECTOR PRIVATE(ccmax, ccran, deltaz, alpha)
      DO jc = i_startidx, i_endidx
        IF (pres(jc,jk)/pres(jc,nlev) > 0.75_wp) THEN
          ccmax = MAX(clc(jc,jk),  zcloud(jc))
          ccran = clc(jc,jk) + zcloud(jc) - clc(jc,jk)*zcloud(jc)

          ! layer thickness [m] between level jk and next upper level jk-1
          deltaz = (pres(jc,jk)-pres(jc,jk-1))/(pres(jc,jk-1)+pres(jc,jk)) * &
                   (temp(jc,jk-1)+temp(jc,jk))*rd/grav

          alpha  = MIN(EXP(-deltaz/zdecorr), clc(jc,jk-1)/MAX(epsi,clc(jc,jk)) )

          zcloud(jc) = alpha * ccmax + (1-alpha) * ccran
        ENDIF
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR PRIVATE(diff_frac_corr)
    DO jc = i_startidx, i_endidx
      IF (cosmu0(jc) > 0.05_wp) THEN
        diff_frac_corr = tune_difrad_3dcont * zcloud(jc) * (1._wp-zcloud(jc))**2

        fr_nir_sfc_diff(jc) = MIN(1._wp, fr_nir_sfc_diff(jc) + diff_frac_corr)
        fr_vis_sfc_diff(jc) = MIN(1._wp, fr_vis_sfc_diff(jc) + diff_frac_corr)
        fr_vis_sfc_diff(jc) = MIN(1._wp, fr_par_sfc_diff(jc) + diff_frac_corr)

        trsol_dn_sfc_diff(jc) = MIN(ecrad_flux%sw_dn(jc,nlev+1)/cosmu0(jc), trsol_dn_sfc_diff(jc) + &
            & diff_frac_corr * ecrad_flux%sw_dn(jc,nlev+1)/cosmu0(jc))
      ENDIF
    ENDDO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE add_3D_diffuse_rad
  !---------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------
  !>
  !! Function create_rdm_seed_in_time:
  !! Create a unique but reproducable random seed in time for the McICA solver
  !!
  FUNCTION create_rdm_seed_in_time(current_datetime)
    !$ACC ROUTINE SEQ
    ! In:
    TYPE(datetime), POINTER, INTENT(in) :: &
      &  current_datetime     !< Current date and time
    ! Out:
    INTEGER              :: &
      &  create_rdm_seed_in_time
    ! Local:
    INTEGER              :: &
      &  time, day            !< time of the day in minutes, the day of the month

    day  = INT(current_datetime%date%day)
    time = (INT(current_datetime%time%hour) * 60) + INT(current_datetime%time%minute)

    create_rdm_seed_in_time = time + day

  END FUNCTION create_rdm_seed_in_time
  !---------------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------------
  !>
  !! Function create_rdm_seed:
  !! Create a unique but reproducable random seed for the McICA solver
  !! 
  !! -------Algorithm taken from IFS, courtesy of R.J. Hogan-------
  !! This method gives a unique value for roughly every 1-km square
  !! on the globe and every minute.  (lat * rad2deg) gives rough
  !! latitude in degrees, which we multiply by 100 to give a unique
  !! value for roughly every km. lon*60*100 gives a unique number
  !! for roughly every km of longitude around the equator, which we
  !! multiply by 180*100 so there is no overlap with the latitude
  !! values.  The result can be contained in a 32-byte integer (but
  !! since random numbers are generated with the help of integer
  !! overflow, it should not matter if the number did overflow).
  !! 
  !! A more simple algorithm using the sum of int(lat), int(lon) and
  !! int(simtime) creates stripe patterns in the instantaneous fluxes.
  !! --------------------------------------------------------------
  !!
  FUNCTION create_rdm_seed(lon,lat,seed_in_time)
    !$ACC ROUTINE SEQ
    ! In:
    REAL(wp),       INTENT(in) :: &
      &  lon, lat             !< Longitude and Latitude value (radian)
    INTEGER, INTENT(in) :: &
      &  seed_in_time     !< contribution to the seed in time
    ! Out:
    INTEGER              :: &
      &  create_rdm_seed

    create_rdm_seed = seed_in_time +  NINT(lon*108000000._wp + (lat * rad2deg * 100._wp) )

  END FUNCTION create_rdm_seed
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! Function reff_crystal:
  !! Function to calculate effective radius of ice crystals (extracted for the use in ecRad
  !! from mo_newcld_optics.f90). Author of the original code: Bjorn Stevens, MPI-M, Hamburg
  !! see ECHAM5 documentation (Roeckner et al, MPI report 349)
  !!
  REAL(wp) FUNCTION reff_crystal(ziwc)  ![m]
    !$ACC ROUTINE SEQ
    REAL(wp), INTENT (IN)  :: &
      &  ziwc                    !< ice water content (g/m3)
    REAL(wp)                :: &
      &  reff_crystal_min,     & !< Minimum value ice crystal effective radius
      &  reff_crystal_max        !< Maximum value ice crystal effective radius
    
    ! Minimum and maximum value (derived from file ECHAM6_CldOptProps.nc)
    reff_crystal_min = 4._wp
    reff_crystal_max = 99._wp ! 124._wp < modified as values > 100 mu m lead to crashes in ecRad
                              ! (if delta_eddington_scat_od = .false.)

    reff_crystal = MAX(reff_crystal_min ,MIN(reff_crystal_max  ,83.8_wp*ziwc**0.216_wp)) * 1.e-6_wp
  END FUNCTION
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! Function reff_droplet:
  !! Function to calculate effective radius of water droplets (extracted for the use in ecRad
  !! from mo_newcld_optics.f90). Author of the original code: Bjorn Stevens, MPI-M, Hamburg
  !!
  REAL(wp) FUNCTION reff_droplet(zlwc,zcdnc,zland,zglac,zfact)  ![m]
    !$ACC ROUTINE SEQ
    REAL (wp), INTENT (IN)  :: &
      &  zlwc,                 & !< liquid water content (g/m3)
      &  zcdnc,                & !< cloud drop number concentration
      &  zglac,                & !< fraction of land covered by glaciers
      &  zland,                & !< land-sea mask. (1. = land, 0. = sea/lakes)
      &  zfact                   !< factor
    REAL(wp)                :: &
      &  zkap,                 & !< Factor
      &  reff_droplet_min,     & !< Minimum value ice crystal effective radius
      &  reff_droplet_max        !< Maximum value ice crystal effective radius
    REAL (wp), PARAMETER ::    &
      &  zkap_cont = 1.143_wp, & !< continental (Martin et al. ) breadth param
      &  zkap_mrtm = 1.077_wp    !< maritime (Martin et al.) breadth parameter

    ! Minimum and maximum value (derived from file ECHAM6_CldOptProps.nc)
    reff_droplet_min = 2.e-6_wp
    reff_droplet_max = 32.e-6_wp

    zkap         = zkap_cont*(zland-zglac) + zkap_mrtm*(1.0_wp-zland+zglac)
    reff_droplet = MAX(reff_droplet_min,MIN(reff_droplet_max,zfact*zkap*(zlwc / (zcdnc * 1.e-6_wp))**(1.0_wp/3.0_wp)))
  END FUNCTION
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! Function gas_profile:
  !! Function to calculate gas profile decaying with height with a tanh function.
  !! Extracted for the use in ecRad from mo_radiation.f90.
  !! Author of the original code: Bjorn Stevens, MPI-M, Hamburg
  !!
  SUBROUTINE gas_profile(vmr_gas, pres, xp, i_startidx, i_endidx, nlev, profile, lacc)

    REAL(wp), INTENT (IN)         :: &
      &  vmr_gas,                    & !< Constant volume mixing ratio of gas specified via namelist
      &  pres(:,:),                  & !< Full level pressure
      &  xp(3)                         !< Gas-specific coefficient
    INTEGER,  INTENT (IN)         :: &
      &  i_startidx, i_endidx,       & !< Start and end index of nproma loop in current block
      &  nlev                          !< Number of vertical full levels
    REAL(wp), INTENT (OUT)        :: &
      &  profile(:,:) !< Profile to be calculated
    LOGICAL, INTENT(IN), OPTIONAL :: &
      &  lacc ! If true, use openacc
    REAL(wp)                      :: &
      &  zx_d, zx_m
    INTEGER                       :: &
      &  jc, jk                        !< loop indices
    LOGICAL                       :: &
      &  lzacc ! non-optional version of lacc

    CALL set_acc_host_or_device(lzacc, lacc)

    zx_m = (vmr_gas+xp(1)*vmr_gas)*0.5_wp
    zx_d = (vmr_gas-xp(1)*vmr_gas)*0.5_wp

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk=1,nlev
      DO jc = i_startidx, i_endidx
        profile(jc,jk) = (1._wp-(zx_d/zx_m)*TANH(LOG(pres(jc,jk) /xp(2)) /xp(3))) * zx_m
      ENDDO !jc
    ENDDO !jk
    !$ACC END PARALLEL

  END SUBROUTINE

  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE get_indices_rad_subblock:
  !! This subroutine gets the index used for subblocking the ecRad radiation.
  !! There are the indices that point to an original array: jcs:jce,      and i_startidx_sub:i_endidx_sub.
  !! There are the indices that point to a  subblock array: 1:nproma_sub, and i_startidx_rad:i_endidx_rad.
  !! In each step of the subblocking, values are copied from an original array to a subblock array and vice versa.
  !! Here, is an illustration for three subblocks and how the index are set.
  !!
  !! jb_rad = 1
  !!
  !! original array:
  !!    !1   i_startidx                                                                   i_endidx    nproma|
  !!    |-------^----------------------------------------------------------------------------^--------------|
  !!    |^      ^                        ^|
  !!    |jcs    |                      jce|
  !!    |       |                        ^|
  !!    |  i_startidx_sub     i_endidx_sub|
  !!
  !! subblock array:
  !!    |1                      nproma_sub|
  !!    |---------------------------------|
  !!    |       ^                        ^|
  !!    |  i_startidx_rad     i_endidx_rad|
  !!
  !! jb_rad = 2
  !!
  !! original array:
  !!    !1   i_startidx                                                                   i_endidx    nproma|
  !!    |-------^----------------------------------------------------------------------------^--------------|
  !!                                      |^                               ^|
  !!                                      |jcs                           jce|
  !!                                      |^                               ^|
  !!                                      |i_startidx_sub       i_endidx_sub|
  !!
  !! subblock array:
  !!                                      |1                      nproma_sub|
  !!                                      |---------------------------------|
  !!                                      |^                               ^|
  !!                                      |i_startidx_rad       i_endidx_rad|
  !!
  !! jb_rad = 3
  !!
  !! original array:
  !!    !1   i_startidx                                                                   i_endidx    nproma|
  !!    |-------^----------------------------------------------------------------------------^--------------|
  !!                                                                        |^               ^             ^  |
  !!                                                                        |jcs             |           jce  |
  !!                                                                        |^               |                |
  !!                                                                        |i_startidx_sub  i_endidx_sub     |
  !!
  !! subblock array:
  !!                                                                        |1                      nproma_sub|
  !!                                                                        |---------------------------------|
  !!                                                                        |^               ^                |
  !!                                                                        |i_startidx_rad  i_endidx_rad     |
  !!
  SUBROUTINE get_indices_rad_subblock(i_startidx, i_endidx, jb_rad, jcs, jce, i_startidx_rad, i_endidx_rad, &
        &  l_3d_rad_fluxes, jnps, jnpe, i_startidx_sub, i_endidx_sub)

    INTEGER, INTENT(IN)  :: i_startidx, i_endidx
    INTEGER, INTENT(IN)  :: jb_rad
    
    INTEGER, INTENT(OUT) :: jcs, jce
    INTEGER, INTENT(OUT) :: i_startidx_rad, i_endidx_rad

    LOGICAL, INTENT(IN), OPTIONAL :: l_3d_rad_fluxes

    INTEGER, INTENT(OUT), OPTIONAL :: jnps, jnpe
    INTEGER, INTENT(OUT), OPTIONAL :: i_startidx_sub, i_endidx_sub

    jcs = nproma_sub*(jb_rad-1) + 1 
    jce = MIN(nproma_sub*jb_rad, nproma)

    i_startidx_rad = MAX(i_startidx-jcs+1, 1)
    i_endidx_rad   = MIN(i_endidx  -jcs+1, nproma_sub)

    IF (PRESENT(l_3d_rad_fluxes) .AND. PRESENT(jnps) .AND. PRESENT(jnpe)) THEN
      IF (l_3d_rad_fluxes) THEN
        jnps = jcs
        jnpe = jce
      ELSE
        jnps = 1
        jnpe = 1
      ENDIF
    ENDIF

    IF (PRESENT(i_startidx_sub) .AND. PRESENT(i_endidx_sub)) THEN
      i_startidx_sub = MAX(jcs, i_startidx)
      i_endidx_sub   = MIN(jce, i_endidx)
    ENDIF

  END SUBROUTINE get_indices_rad_subblock
  !---------------------------------------------------------------------------------------

#endif
END MODULE mo_nwp_ecrad_utilities
