!
! This module prepares aerosol for the use in radiation
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
MODULE mo_nwp_aerosol

! ICON infrastructure
  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: finish, message, message_text
  USE mo_model_domain,            ONLY: t_patch
  USE mo_grid_config,             ONLY: nroot 
  USE mo_ext_data_types,          ONLY: t_external_data
  USE mo_nonhydro_types,          ONLY: t_nh_diag
  USE mo_nwp_phy_types,           ONLY: t_nwp_phy_diag
  USE mo_parallel_config,         ONLY: nproma
  USE mo_loopindices,             ONLY: get_indices_c
  USE mo_impl_constants,          ONLY: min_rlcell_int, SUCCESS, &
                                    &   iss, iorg, ibc, iso4, idu, n_camsaermr
  USE mo_impl_constants_grf,      ONLY: grf_bdywidth_c
  USE mo_physical_constants,      ONLY: rd, grav, cpd
  USE mo_reader_cams,             ONLY: t_cams_reader
  USE mo_interpolate_time,        ONLY: t_time_intp
  USE mo_io_units,                ONLY: filename_max
  USE mo_fortran_tools,           ONLY: set_acc_host_or_device
  USE mo_util_string,             ONLY: int2string, associate_keyword, t_keyword_list, with_keywords
! ICON configuration
  USE mo_atm_phy_nwp_config,      ONLY: atm_phy_nwp_config, iprog_aero, icpl_aero_conv
  USE mo_radiation_config,        ONLY: irad_aero, iRadAeroConstKinne, iRadAeroKinne, iRadAeroCAMSclim, &
                                    &   iRadAeroVolc, iRadAeroKinneVolc, iRadAeroART,                   &
                                    &   iRadAeroKinneVolcSP, iRadAeroKinneSP, iRadAeroTegen,            &
                                    &   cams_clim_filename
! External infrastruture
  USE mtime,                      ONLY: datetime, timedelta, newDatetime, newTimedelta,       &
                                    &   operator(+), deallocateTimedelta, deallocateDatetime
! Aerosol-specific
  USE mo_aerosol_util,            ONLY: aerdis
  USE mo_bc_aeropt_kinne,         ONLY: read_bc_aeropt_kinne, set_bc_aeropt_kinne
  USE mo_bc_aeropt_cmip6_volc,    ONLY: read_bc_aeropt_cmip6_volc, add_bc_aeropt_cmip6_volc
  USE mo_bc_aeropt_splumes,       ONLY: add_bc_aeropt_splumes
  USE mo_bcs_time_interpolation,  ONLY: t_time_interpolation_weights,         &
    &                                   calculate_time_interpolation_weights

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nwp_aerosol'

  PUBLIC :: nwp_aerosol_interface
  PUBLIC :: nwp_aerosol_cleanup
  PUBLIC :: nwp_aerosol_init
  PUBLIC :: cams_reader
  PUBLIC :: cams_intp

  TYPE(t_cams_reader),      ALLOCATABLE, TARGET :: cams_reader(:)
  TYPE(t_time_intp),        ALLOCATABLE         :: cams_intp(:)

CONTAINS

  !---------------------------------------------------------------------------------------
  !! This subroutine uploads CAMS aerosols mixing ratios 3D climatology and updates them once a day
  SUBROUTINE nwp_aerosol_init(mtime_datetime, p_patch)

    TYPE(datetime), POINTER, INTENT(in) :: &
      &  mtime_datetime                      !< Current datetime
    TYPE(t_patch), INTENT(in)           :: &
      &  p_patch
    ! Local variables
    INTEGER ::                 &
      &  jg                                  !< Domain index
    CHARACTER(LEN=filename_max) :: &
      &  cams_clim_td_file                   !< CAMS climatology file name

    IF (irad_aero == iRadAeroCAMSclim) THEN

      jg     = p_patch%id

      cams_clim_td_file = generate_cams_clim_filename(TRIM(cams_clim_filename), nroot, p_patch%level, p_patch%id)
      CALL message  ('nwp_aerosol_init opening CAMS 3D climatology file: ', TRIM(cams_clim_td_file))

      CALL cams_reader(jg)%init(p_patch, TRIM(cams_clim_td_file))
      CALL cams_intp(jg)%init(cams_reader(jg), mtime_datetime, '', 2)

    ENDIF

    CONTAINS

      FUNCTION generate_cams_clim_filename(filename_in, nroot, jlev, idom) RESULT(result_str)
        CHARACTER(filename_max)                 :: result_str
        CHARACTER(LEN=*), INTENT(in)            :: filename_in
        INTEGER,                     INTENT(in) :: nroot, jlev, idom
        ! Local variables
        TYPE (t_keyword_list), POINTER      :: keywords => NULL()

        CALL associate_keyword("<nroot>",    TRIM(int2string(nroot,"(i0)")),   keywords)
        CALL associate_keyword("<nroot0>",   TRIM(int2string(nroot,"(i2.2)")), keywords)
        CALL associate_keyword("<jlev>",     TRIM(int2string(jlev, "(i2.2)")), keywords)
        CALL associate_keyword("<idom>",     TRIM(int2string(idom, "(i2.2)")), keywords)

        result_str = TRIM(with_keywords(keywords, TRIM(filename_in)))
      END FUNCTION generate_cams_clim_filename

  END SUBROUTINE nwp_aerosol_init

  !---------------------------------------------------------------------------------------
  SUBROUTINE nwp_aerosol_interface(mtime_datetime, pt_patch, ext_data, pt_diag, prm_diag,          &
    &                              zf, zh, dz, dt_rad,                                             &
    &                              inwp_radiation, nbands_lw, nbands_sw, wavenum1_sw, wavenum2_sw, &
    &                              zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,                              &
    &                              od_lw, od_sw, ssa_sw, g_sw, lacc)
    CHARACTER(len=*), PARAMETER :: &
      &  routine = modname//':nwp_aerosol_interface'

    TYPE(datetime), POINTER, INTENT(in) :: &
      &  mtime_datetime          !< Current datetime
    TYPE(t_patch), TARGET, INTENT(in) :: &
      &  pt_patch                !< Grid/patch info
    TYPE(t_external_data), INTENT(inout) :: &
      &  ext_data                !< External data
    TYPE(t_nh_diag), TARGET, INTENT(inout) :: &
      &  pt_diag                 !< the diagnostic variables
    TYPE(t_nwp_phy_diag), INTENT(inout) :: &
      &  prm_diag                !< Physics diagnostics
    REAL(wp), INTENT(in) ::    &
      &  zf(:,:,:), zh(:,:,:), & !< model full/half layer height
      &  dz(:,:,:),            & !< Layer thickness
      &  dt_rad                  !< Radiation time step
    REAL(wp), POINTER, INTENT(in) :: &
      &  wavenum1_sw(:),       & !< Shortwave wavenumber lower band bounds
      &  wavenum2_sw(:)          !< Shortwave wavenumber upper band bounds
    REAL(wp), INTENT(inout) :: &
      &  zaeq1(:,:,:),         & !< Tegen optical thicknesses       1: continental
      &  zaeq2(:,:,:),         & !< relative to 550 nm, including   2: maritime
      &  zaeq3(:,:,:),         & !< a vertical profile              3: desert
      &  zaeq4(:,:,:),         & !< for 5 different                 4: urban
      &  zaeq5(:,:,:)            !< aerosol species.                5: stratospheric background
    INTEGER, INTENT(in) ::     &
      &  inwp_radiation,       & !< Radiation scheme (1=rrtmg, 4=ecrad)
      &  nbands_lw, nbands_sw    !< Number of short and long wave bands
    REAL(wp), ALLOCATABLE, INTENT(out) :: &
      &  od_lw(:,:,:,:),       & !< Longwave optical thickness
      &  od_sw(:,:,:,:),       & !< Shortwave optical thickness
      &  ssa_sw(:,:,:,:),      & !< Shortwave asymmetry factor
      &  g_sw(:,:,:,:)           !< Shortwave single scattering albedo
    LOGICAL, OPTIONAL, INTENT(in) :: lacc ! If true, use openacc
! Local variables



    REAL(wp) ::                &
      &  latitude(nproma),     & !< Geographical latitude
      &  time_weight             !< Weihting for temporal interpolation
    REAL(wp),  ALLOCATABLE ::  &
      &  cams(:,:,:,:)           !< CAMS climatology fields taken from external file 
    INTEGER ::                 &
      &  jk, jc, jb, jt,       &
      &  rl_start, rl_end,     &
      &  i_startblk, i_endblk, &
      &  i_startidx, i_endidx, &
      &  istat,                & !< Error code
      &  imo1 , imo2             !< Month index (current and next month)
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)
  
    SELECT CASE(irad_aero)
!---------------------------------------------------------------------------------------
! Tegen aerosol (+ART if chosen)
!---------------------------------------------------------------------------------------
      CASE(iRadAeroTegen, iRadAeroART)

        !$ACC DATA CREATE(latitude) IF(lzacc)

        ! Outer two rows need dummy values as RRTM always starts at 1
        rl_start   = 1
        rl_end     = 2
        i_startblk = pt_patch%cells%start_block(rl_start)
        i_endblk   = pt_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)  ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(pt_patch,jb,i_startblk,i_endblk,i_startidx,i_endidx,rl_start,rl_end)
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 1, pt_patch%nlev
            DO jc = i_startidx,i_endidx
              zaeq1(jc,jk,jb) = 0._wp
              zaeq2(jc,jk,jb) = 0._wp
              zaeq3(jc,jk,jb) = 0._wp
              zaeq4(jc,jk,jb) = 0._wp
              zaeq5(jc,jk,jb) = 0._wp
              ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

        ! Start at third row instead of fifth as two rows are needed by the reduced grid aggregation
        rl_start   = grf_bdywidth_c-1
        rl_end     = min_rlcell_int
        i_startblk = pt_patch%cells%start_block(rl_start)
        i_endblk   = pt_patch%cells%end_block(rl_end)

        ! Calculate the weighting factor and month indices for temporal interpolation
        CALL get_time_intp_weights(mtime_datetime, imo1 , imo2, time_weight)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,latitude)  ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(pt_patch,jb,i_startblk,i_endblk,i_startidx,i_endidx,rl_start,rl_end)

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx,i_endidx
            latitude(jc) = pt_patch%cells%center(jc,jb)%lat
          ENDDO
          !$ACC END PARALLEL

          CALL nwp_aerosol_tegen(i_startidx, i_endidx, pt_patch%nlev, pt_patch%nlevp1, &
            &                    prm_diag%k850(:,jb), pt_diag%temp(:,:,jb), pt_diag%pres(:,:,jb), pt_diag%pres_ifc(:,:,jb), &
            &                    ext_data%atm_td%aer_ss(:,jb,imo1),   ext_data%atm_td%aer_org(:,jb,imo1), &
            &                    ext_data%atm_td%aer_bc(:,jb,imo1),   ext_data%atm_td%aer_so4(:,jb,imo1), &
            &                    ext_data%atm_td%aer_dust(:,jb,imo1), ext_data%atm_td%aer_ss(:,jb,imo2),  &
            &                    ext_data%atm_td%aer_org(:,jb,imo2),  ext_data%atm_td%aer_bc(:,jb,imo2),  &
            &                    ext_data%atm_td%aer_so4(:,jb,imo2),  ext_data%atm_td%aer_dust(:,jb,imo2),&
            &                    prm_diag%pref_aerdis(:,jb),latitude(:), pt_diag%dpres_mc(:,:,jb), time_weight, &
            &                    prm_diag%aerosol(:,:,jb), prm_diag%aercl_ss(:,jb), prm_diag%aercl_or(:,jb), &
            &                    prm_diag%aercl_bc(:,jb), prm_diag%aercl_su(:,jb), prm_diag%aercl_du(:,jb), &
            &                    zaeq1(:,:,jb),zaeq2(:,:,jb),zaeq3(:,:,jb),zaeq4(:,:,jb),zaeq5(:,:,jb),lacc )

          ! This is where ART should be placed

          ! Compute cloud number concentration depending on aerosol climatology if 
          ! aerosol-microphysics or aerosol-convection coupling is turned on
          IF (atm_phy_nwp_config(pt_patch%id)%icpl_aero_gscp == 1 .OR. icpl_aero_conv == 1) THEN
            CALL nwp_cpl_aero_gscp_conv(i_startidx, i_endidx, pt_patch%nlev, pt_diag%pres_sfc(:,jb), pt_diag%pres(:,:,jb), &
              &                         prm_diag%acdnc(:,:,jb), prm_diag%cloud_num(:,jb), lacc)
          ENDIF

        ENDDO !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL



        !$ACC WAIT
        !$ACC END DATA

!---------------------------------------------------------------------------------------
! Kinne aerosol
!---------------------------------------------------------------------------------------
      CASE(iRadAeroConstKinne, iRadAeroKinne, iRadAeroVolc, iRadAeroKinneVolc, iRadAeroKinneVolcSP, iRadAeroKinneSP)

        rl_start   = grf_bdywidth_c-1
        rl_end     = min_rlcell_int
        i_startblk = pt_patch%cells%start_block(rl_start)
        i_endblk   = pt_patch%cells%end_block(rl_end)
      
        ! Compatibility checks
        WRITE(message_text,'(a,i2,a)') 'irad_aero = ', irad_aero,' requires to compile with --enable-ecrad.'
        CALL finish(routine, message_text)

        ! Update Kinne aerosol from files once per day
        CALL nwp_aerosol_daily_update_kinne(mtime_datetime, pt_patch, dt_rad, inwp_radiation, &
          &                                 nbands_lw, nbands_sw)

        ! Allocations
        ALLOCATE(od_lw (nproma,pt_patch%nlev,pt_patch%nblks_c,nbands_lw)  , &
          &      od_sw (nproma,pt_patch%nlev,pt_patch%nblks_c,nbands_sw)  , &
          &      ssa_sw(nproma,pt_patch%nlev,pt_patch%nblks_c,nbands_sw)  , &
          &      g_sw  (nproma,pt_patch%nlev,pt_patch%nblks_c,nbands_sw)  , &
          &      STAT=istat)
        IF(istat /= SUCCESS) &
          &  CALL finish(routine, 'Allocation of od_lw, od_sw, ssa_sw, g_sw failed')

        IF ( .NOT. ASSOCIATED(wavenum1_sw) .OR. .NOT. ASSOCIATED(wavenum2_sw) ) &
          &  CALL finish(routine, 'wavenum1 or wavenum2 not associated')
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(pt_patch,jb,i_startblk,i_endblk,i_startidx,i_endidx,rl_start,rl_end)
          IF (i_startidx>i_endidx) CYCLE

          CALL nwp_aerosol_kinne(mtime_datetime, zf(:,:,jb), zh(:,:,jb), dz(:,:,jb),   &
            &                    pt_patch%id, jb, i_endidx, pt_patch%nlev,             &
            &                    nbands_lw, nbands_sw, wavenum1_sw(:), wavenum2_sw(:), &
            &                    od_lw(:,:,jb,:), od_sw(:,:,jb,:),                     &
            &                    ssa_sw(:,:,jb,:), g_sw(:,:,jb,:)                      )

        END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ! CAMS climatology aerosol
      CASE(iRadAeroCAMSclim)

        rl_start   = grf_bdywidth_c+1
        rl_end     = min_rlcell_int
        i_startblk = pt_patch%cells%start_block(rl_start)
        i_endblk   = pt_patch%cells%end_block(rl_end)

        ! Compatibility checks
        WRITE(message_text,'(a,i2,a)') 'irad_aero = ', irad_aero,' requires to compile with --enable-ecrad.'
        CALL finish(routine, message_text)

        CALL nwp_aerosol_update_cams(mtime_datetime, pt_patch%id, cams)

        DO jt=1, n_camsaermr
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
          DO jb = i_startblk,i_endblk
            CALL get_indices_c(pt_patch,jb,i_startblk,i_endblk,i_startidx,i_endidx,rl_start,rl_end)
            IF (i_startidx>i_endidx) CYCLE
            pt_diag%camsaermr(:,:,jb,jt) = 0._wp
            CALL vinterp_cams_aerosols(i_startidx, i_endidx, pt_diag%pres_ifc(:,:,jb),    &
              &                   cams_pres_in=cams(:,:,jb,n_camsaermr+1), cams=cams(:,:,jb,jt), &
              &                   nlev=pt_patch%nlev, camsaermr=pt_diag%camsaermr(:,:,jb,jt))
          END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
        END DO ! jt aerosol loop

      CASE DEFAULT
        ! Currently continue as not all cases are ported to nwp_aerosol_interface yet
    END SELECT

  END SUBROUTINE nwp_aerosol_interface

  !---------------------------------------------------------------------------------------
  SUBROUTINE nwp_aerosol_daily_update_kinne(mtime_datetime, pt_patch, dt_rad, inwp_radiation, nbands_lw, nbands_sw)
    TYPE(datetime), POINTER, INTENT(in) :: &
      &  mtime_datetime                    !< Current datetime
    TYPE(t_patch), TARGET, INTENT(in) :: &
      &  pt_patch                          !< Grid/patch info
    REAL(wp), INTENT(in) ::              &
      &  dt_rad                            !< Radiation time step
    INTEGER, INTENT(in) ::               &
      &  inwp_radiation,                 & !< Radiation scheme (1=rrtmg, 4=ecrad)
      &  nbands_lw, nbands_sw              !< Number of short and long wave bands
    ! Local variables
    TYPE(datetime), POINTER ::           &
      &  prev_radtime                      !< Datetime of previous radiation time step
    TYPE(timedelta), POINTER ::          &
      &  td_dt_rad                         !< Radiation time step

    td_dt_rad => newTimedelta('-',0,0,0,0,0, second=NINT(dt_rad), ms=0)
    prev_radtime => newDatetime(mtime_datetime + td_dt_rad)

    IF (prev_radtime%date%day /= mtime_datetime%date%day) THEN
      IF (inwp_radiation == 4) THEN
        IF (ANY(irad_aero == [iRadAeroKinne, iRadAeroKinneVolc])) &
            & CALL read_bc_aeropt_kinne(mtime_datetime, pt_patch, .TRUE., nbands_lw, nbands_sw)
        IF (ANY(irad_aero == [iRadAeroVolc, iRadAeroKinneVolc, iRadAeroKinneVolcSP])) &
            & CALL read_bc_aeropt_cmip6_volc(mtime_datetime, nbands_lw, nbands_sw)
      ENDIF
    ENDIF

    CALL deallocateTimedelta(td_dt_rad)
    CALL deallocateDatetime(prev_radtime)

  END SUBROUTINE nwp_aerosol_daily_update_kinne

  !---------------------------------------------------------------------------------------
  SUBROUTINE nwp_aerosol_kinne(mtime_datetime, zf, zh, dz, jg, jb, i_endidx, nlev, &
    &                          nbands_lw, nbands_sw, wavenum1_sw, wavenum2_sw,     &
    &                          od_lw, od_sw, ssa_sw, g_sw)
    TYPE(datetime), POINTER, INTENT(in) :: &
      &  mtime_datetime                      !< Current datetime
    REAL(wp), INTENT(in) ::                &
      &  zf(:,:), zh(:,:), dz(:,:),        & !< model full/half layer height, layer thickness
      &  wavenum1_sw(:),                   & !< Shortwave wavenumber lower band bounds
      &  wavenum2_sw(:)                      !< Shortwave wavenumber upper band bounds
    INTEGER, INTENT(in) ::                 &
      &  jg, jb,                           & !< Domain and block index
      &  i_endidx,                         & !< Loop bound
      &  nlev,                             & !< Number of vertical levels
      &  nbands_lw, nbands_sw                !< Number of short and long wave bands
    REAL(wp), INTENT(out) ::               &
      &  od_lw(:,:,:), od_sw(:,:,:),       & !< LW/SW optical thickness
      &  ssa_sw(:,:,:), g_sw(:,:,:)          !< SW asymmetry factor, SW single scattering albedo
    ! Local variables
    REAL(wp) ::                            &
      &  od_lw_vr (nproma,nlev,nbands_lw), & !< LW optical thickness of aerosols    (vertically reversed)
      &  od_sw_vr (nproma,nlev,nbands_sw), & !< SW aerosol optical thickness        (vertically reversed)
      &  g_sw_vr  (nproma,nlev,nbands_sw), & !< SW aerosol asymmetry factor         (vertically reversed)
      &  ssa_sw_vr(nproma,nlev,nbands_sw)    !< SW aerosol single scattering albedo (vertically reversed)
    REAL(wp) ::                            &
      &  x_cdnc(nproma)                      !< Scale factor for Cloud Droplet Number Concentration (currently not used)
    INTEGER ::                             &
      &  jk                                  !< Loop index

    od_lw_vr(:,:,:)  = 0.0_wp
    od_sw_vr(:,:,:)  = 0.0_wp
    ssa_sw_vr(:,:,:) = 1.0_wp
    g_sw_vr (:,:,:)  = 0.0_wp

    ! Tropospheric Kinne aerosol
    IF (ANY( irad_aero == (/iRadAeroConstKinne,iRadAeroKinne,iRadAeroKinneVolc, &
      &                     iRadAeroKinneVolcSP,iRadAeroKinneSP/) )) THEN
      CALL set_bc_aeropt_kinne(mtime_datetime, jg, 1, i_endidx, nproma, nlev, jb, &
        &                      nbands_sw, nbands_lw, zf(:,:), dz(:,:),            &
        &                      od_sw_vr(:,:,:), ssa_sw_vr(:,:,:),                 &
        &                      g_sw_vr (:,:,:), od_lw_vr(:,:,:)                   )
    ENDIF

    ! Volcanic stratospheric aerosols for CMIP6
    IF (ANY( irad_aero == (/iRadAeroVolc,iRadAeroKinneVolc,iRadAeroKinneVolcSP/) )) THEN 
     CALL add_bc_aeropt_cmip6_volc(mtime_datetime, jg, 1, i_endidx, nproma, nlev, jb, &
       &                           nbands_sw, nbands_lw, zf(:,:), dz(:,:),            &
       &                           od_sw_vr(:,:,:), ssa_sw_vr(:,:,:),                 &
       &                           g_sw_vr (:,:,:), od_lw_vr(:,:,:)                   )
    END IF

    ! Simple plumes
    IF (ANY( irad_aero == (/iRadAeroKinneVolcSP,iRadAeroKinneSP/) )) THEN
      CALL add_bc_aeropt_splumes(jg, 1, i_endidx, nproma, nlev, jb,  &
        &                        nbands_sw, mtime_datetime,          &
        &                        zf(:,:), dz(:,:), zh(:,nlev+1),     &
        &                        wavenum1_sw(:), wavenum2_sw(:),     &
        &                        od_sw_vr(:,:,:), ssa_sw_vr(:,:,:),  &
        &                        g_sw_vr (:,:,:), x_cdnc(:)          )
    END IF

    ! Vertically reverse the fields:
    DO jk = 1, nlev
      od_lw (:,jk,:) = od_lw_vr (:,nlev-jk+1,:)
      od_sw (:,jk,:) = od_sw_vr (:,nlev-jk+1,:)
      ssa_sw(:,jk,:) = ssa_sw_vr(:,nlev-jk+1,:)
      g_sw  (:,jk,:) = g_sw_vr  (:,nlev-jk+1,:)
    ENDDO

  END SUBROUTINE nwp_aerosol_kinne


  !---------------------------------------------------------------------------------------
  !! This subroutine uploads CAMS aerosols mixing ratios 3D climatology and updates them
  SUBROUTINE nwp_aerosol_update_cams(mtime_datetime, jg, cams)

    TYPE(datetime), POINTER, INTENT(in)  :: &
      &  mtime_datetime                    !< Current datetime
    REAL(wp), ALLOCATABLE, INTENT(inout) :: &
      &  cams(:,:,:,:)                     !< CAMS fields taken from external file
    INTEGER, INTENT(in)                  :: &
      &  jg                                !< Domain index
    ! Local variables
    REAL(wp), ALLOCATABLE                :: &
      &  camscl_dat(:,:,:,:)

    ALLOCATE(cams( nproma, cams_reader(jg)%nlev_cams, cams_reader(jg)%p_patch%nblks_c, n_camsaermr+1 ))
    cams(:,:,:,:) = 0.0_wp

    CALL cams_intp(jg)%intp(mtime_datetime, camscl_dat)

    cams(:,:,:,:) = camscl_dat(:,:,:,:)

  END SUBROUTINE nwp_aerosol_update_cams

  !---------------------------------------------------------------------------------------
  !! Vertical interpolation of CAMS aerosols mixing ratios 3D climatology
  !! This routine uses the original data of CAMS climatology which has different
  !! number of pressure levels from ICON and also different lowest/highest pressure values. 
  !! Therefore, we first move to sigma coordinates in both models. In both ICON and CAMS
  !! half level pressure are used. The layer directly above surface (nk1+1) pressure thickness 
  !! is absent from the original CAMS file and was set to 240 Pa for all grid points (evaluated from CAMS data).
  !! The code loops on each ICON sigma level and calculates the amount of layer-integrated
  !! mass that is confined within it. For each icon sigma layer we find lmax which is the closest
  !! CAMS level below icon_sigma2 (jk+1) and lmin which is the closest level above 
  !! icon_sigma1 (jk). There are 3 cases:
  !! CASE 1: Current ICON layer completely within one CAMS layer
  !!                       -----lmin------
  !!   --------jk-------
  !!   -------jk+1------
  !!                       -----lmax------  
  !! CASE 2: Current ICON layer covers 2 CAMS layer
  !!                       -----lmin------
  !!   --------jk-------
  !!                       ----lmin+1-----  
  !!   -------jk+1------
  !!                       -----lmax------
  !! CASE 3: Current ICON layer covers a few CAMS layer
  !!                       -----lmin------
  !!   --------jk-------
  !!                       ----lmin+1-----  
  !!                       ---- ... ------ 
  !!                       ----lmax-1----- 
  !!   -------jk+1------
  !!                       -----lmax------
  !!
  SUBROUTINE vinterp_cams_aerosols(i_startidx, i_endidx, pres, cams_pres_in, cams, nlev, camsaermr )

    REAL(wp),      INTENT(in)      :: &
      &  cams(:,:),                      & !< CAMS fields taken from external file; layer integrated mass [kg/m^2]
      &  pres(:,:),                      & !< ICON diagnosed pressure
      &  cams_pres_in(:,:)                 !< CAMS climatology pressure taken from external file
    REAL(wp),      INTENT(inout)   :: &
      &  camsaermr(:,:)                    !< CAMS aerosols mixing ratios [kg/kg]
    INTEGER,       INTENT(in)      :: &
      &  nlev,                           & !< ICON number of vertical levels
      &  i_startidx,                     & !< Loop indices
      &  i_endidx                          !< Loop indices

    ! local variables
    REAL(wp)                       :: &
      &  g_dp,                           & !< graviational constant divided by pressure thickness
      &  icon_sigma1, icon_sigma2,       & !< ICON pressure levels
      &  delta1, delta, delt1, delt,     & !< variables for looping
      &  dp_cams_surface,                & !< cams near surface pressure thickness
      &  layer_mass                        !< mass at specific layer
    REAL(wp) , ALLOCATABLE         :: &
      &  cams_sigma(:,:),                & !< CAMS sigma coordinate
      &  icon_sigma(:,:)                   !< ICON sigma coordinate
    INTEGER                        :: &
      &  jc, jk, jk1, k, lmin, lmax,     & !< Loop indices
      &  nk1                               !< number of vertical levels in original CAMS climatology data

    CHARACTER(len=*), PARAMETER    :: &
      &  routine = modname//':vinterp_cams_aerosols'

    nk1 = size(cams_pres_in,2)  
 
    ALLOCATE(cams_sigma(size(cams_pres_in,1),size(cams_pres_in,2)+1))
    ALLOCATE(icon_sigma(size(pres,1),size(pres,2)))

    dp_cams_surface = 240.0_wp

    ! move to sigma coordinate  
    DO jc = i_startidx, i_endidx 
      cams_sigma(jc,nk1+1) = 1.0_wp
      DO jk1 = 1, nk1
        cams_sigma(jc,jk1) = cams_pres_in(jc,jk1)/(cams_pres_in(jc,nk1) + dp_cams_surface)
      ENDDO
      DO jk = 1, nlev+1 ! loop on icon vertical pressure ifc levels
        icon_sigma(jc,jk) = pres(jc,jk)/pres(jc,nlev+1)
      ENDDO
    ENDDO

    DO jc = i_startidx, i_endidx ! loop on icon horizontal index
      DO jk = 1, nlev ! loop on icon vertical levels

        camsaermr(jc,jk) = 0.0_wp
        layer_mass = 0.0_wp

        icon_sigma1 = icon_sigma(jc,jk)
        icon_sigma2 = icon_sigma(jc,jk+1)
        g_dp = grav/(pres(jc,nlev+1)*(icon_sigma2-icon_sigma1))

        ! Exclude all ICON levels which have lower pressure than CAMS lowest pressure
        IF (icon_sigma1 > cams_sigma(jc,1)) THEN

          ! loop to find lmin & lmax that are above/below jk and jk+1 respectively
          delta = 1.0_wp
          delt = 1.0_wp

          DO jk1 = 1,nk1+1
            delta1 = icon_sigma1 - cams_sigma(jc,jk1) 
            IF (delta1 >= 0.0_wp .AND. ABS(delta1) < delta) THEN
              lmin = jk1 
              delta = ABS(delta1)
            END IF

            delt1 = cams_sigma(jc,jk1)-icon_sigma2 
            IF (delt1 >= 0.0_wp .AND. ABS(delt1) < delt) THEN
              lmax = jk1 
              delt = ABS(delt1)
            END IF
          ENDDO

          ! accumulate all mass in between jk and jk+1
          ! Current ICON layer completely within one cams layer (CASE 1)
          IF ( lmax == (lmin+1) ) THEN
          layer_mass = layer_mass + cams(jc,lmin) *            &
            &         ( (icon_sigma2-icon_sigma1) &
            &         / (cams_sigma(jc,lmax)-cams_sigma(jc,lmin)) )

          ! Current ICON layer covers more than one CAMS layer (CASE 2,3)
          ELSE          
            ! this IF is for CASE 3 only
            IF (lmax > lmin + 2) THEN
              DO k = lmin+1, lmax-2
                layer_mass = layer_mass + cams(jc,k) 
              END DO           
            END IF

            layer_mass = layer_mass + cams(jc,lmin) *         &
              &         ( (cams_sigma(jc,lmin+1)-icon_sigma1) &
              &         / (cams_sigma(jc,lmin+1)-cams_sigma(jc,lmin)) )

            layer_mass = layer_mass + cams(jc,lmax-1) *       &
              &         ( (icon_sigma2-cams_sigma(jc,lmax-1)) &
              &         / (cams_sigma(jc,lmax)-cams_sigma(jc,lmax-1)) )

          END IF

          ! add upper (lowest pressure) CAMS levels mass 
          ! which does not have any corresponding ICON level
          ! add this mass to ICON level jk = 1 (lowest pressure level) 
          IF (jk == 1 .AND. lmin > 1) THEN
            DO k = 1, lmin-1
              layer_mass = layer_mass + cams(jc,k)
            END DO
              layer_mass = layer_mass + cams(jc,lmin) *       &
                &         ( (icon_sigma1-cams_sigma(jc,lmin)) &
                &         / (cams_sigma(jc,lmin+1)-cams_sigma(jc,lmin)) )
          END IF

        END IF

        ! Include the case where ICON level partialy covers the first CAMS level
        IF (icon_sigma1 < cams_sigma(jc,1) .AND. icon_sigma2 > cams_sigma(jc,1) ) THEN
          layer_mass = layer_mass + cams(jc,1) *       &
            &         ( (icon_sigma2-cams_sigma(jc,1)) &
            &         / (cams_sigma(jc,2)-cams_sigma(jc,1)) )
        END IF

        ! move from integrated mass (kg/m^2) to mixing ratios
        camsaermr(jc,jk) = layer_mass * g_dp

        ! for checking if all gridpoints have meaningful values
        IF (camsaermr(jc,jk) < 0.0_wp) THEN
          CALL finish(routine,'mo_nwp_aerosol: vinterp_cams_aerosols failed')
        END IF

      ENDDO !jk
    ENDDO !jc

    DEALLOCATE(cams_sigma)
    DEALLOCATE(icon_sigma)

  END SUBROUTINE vinterp_cams_aerosols

  !---------------------------------------------------------------------------------------
  SUBROUTINE nwp_aerosol_tegen ( istart, iend, nlev, nlevp1, k850, temp, pres, pres_ifc,                 &
    &                            aer_ss_mo1, aer_org_mo1, aer_bc_mo1, aer_so4_mo1, aer_dust_mo1,         &
    &                            aer_ss_mo2, aer_org_mo2, aer_bc_mo2, aer_so4_mo2, aer_dust_mo2,         &
    &                            pref_aerdis,latitude, dpres_mc, time_weight,                            &
    &                            aerosol, aercl_ss, aercl_or, aercl_bc, aercl_su, aercl_du, &
    &                            zaeq1,zaeq2,zaeq3,zaeq4,zaeq5,lacc )

    INTEGER,  INTENT(in)                :: &
      &  nlev, nlevp1,                     & !< Number of vertical full/half levels
      &  istart, iend,                     & !< Start and end index of jc loop
      &  k850(:)                             !< Index of 850 hPa layer
    REAL(wp), INTENT(in)                :: &
      &  temp(:,:), pres(:,:),             & !< temperature and pressure at full level
      &  pres_ifc(:,:),                    & !< pressure at half level
      &  aer_ss_mo1(:), aer_org_mo1(:),    & !< Month 1 climatology from extpar file (sea salt, organic)
      &  aer_bc_mo1(:), aer_so4_mo1(:),    & !< Month 1 climatology from extpar file (blck carbon, sulphate)
      &  aer_dust_mo1(:),                  & !< Month 1 climatology from extpar file (dust)
      &  aer_ss_mo2(:), aer_org_mo2(:),    & !< Month 2 climatology from extpar file (sea salt, organic)
      &  aer_bc_mo2(:), aer_so4_mo2(:),    & !< Month 2 climatology from extpar file (blck carbon, sulphate)
      &  aer_dust_mo2(:),                  & !< Month 2 climatology from extpar file (dust)
      &  pref_aerdis(:),                   & !< Reference pressure for vertical distribution of aerosol
      &  latitude(:),                      & !< geographical latitude
      &  dpres_mc(:,:),                    & !< pressure thickness
      &  time_weight                         !< Temporal weighting factor
    REAL(wp), TARGET, INTENT(inout)     :: &
      &  aerosol(:,:),                     & !< Aerosol field incl. temporal interpolation
      &  aercl_ss(:), aercl_or(:),         & !< Climatological fields for relaxation (iprog_aero > 0)
      &  aercl_bc(:), aercl_su(:),         & !< Climatological fields for relaxation (iprog_aero > 0)
      &  aercl_du(:)                         !< Climatological fields for relaxation (iprog_aero > 0)
    REAL(wp), INTENT(inout)             :: &
      &  zaeq1(:,:), zaeq2(:,:),           & !< organics, sea salt
      &  zaeq3(:,:), zaeq4(:,:),           & !< dust, black carbon
      &  zaeq5(:,:)                          !< sulphate (incl. stratospheric background zstbga)
    LOGICAL, INTENT(in), OPTIONAL       :: &
      &  lacc                                !< If true, use openacc
! Local variables
    CHARACTER(len=*), PARAMETER         :: &
      &  routine = modname//':nwp_aerosol_tegen'
    REAL(wp)                            :: &
      &  zsign (nproma,nlevp1),            &
      &  zvdaes(nproma,nlevp1),            &
      &  zvdael(nproma,nlevp1),            &
      &  zvdaeu(nproma,nlevp1),            &
      &  zvdaed(nproma,nlevp1),            &
      &  zaeqdo (nproma), zaeqdn,          &
      &  zaequo (nproma), zaequn,          &
      &  zaeqlo (nproma), zaeqln,          &
      &  zaeqsuo(nproma), zaeqsun,         &
      &  zaeqso (nproma), zaeqsn,          &
      &  zptrop (nproma), zdtdz(nproma),   &
      &  zlatfac(nproma), zstrfac,         &
      &  zpblfac, zslatq
    REAL(wp), PARAMETER                 :: &
      & ztrbga = 0.03_wp  / (101325.0_wp - 19330.0_wp), &
      ! original value for zstbga of 0.045 is much higher than recently published climatologies
      & zstbga = 0.015_wp  / 19330.0_wp
    INTEGER                             :: &
      &  jc,jk                               !< Loop indices
    LOGICAL                             :: &
      &  lzacc                               !< non-optional version of lacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA CREATE(zsign, zvdaes, zvdael, zvdaeu, zvdaed) &
    !$ACC   CREATE(zaeqdo, zaequo, zaeqlo, zaeqsuo, zaeqso, zptrop) &
    !$ACC   CREATE(zdtdz, zlatfac) IF(lzacc)

    SELECT CASE(iprog_aero)
      CASE(0) ! Purely climatological aerosol
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO jc = istart, iend
          aerosol(jc,iss)  = aer_ss_mo1  (jc) + ( aer_ss_mo2  (jc) - aer_ss_mo1  (jc) ) * time_weight
          aerosol(jc,iorg) = aer_org_mo1 (jc) + ( aer_org_mo2 (jc) - aer_org_mo1 (jc) ) * time_weight
          aerosol(jc,ibc)  = aer_bc_mo1  (jc) + ( aer_bc_mo2  (jc) - aer_bc_mo1  (jc) ) * time_weight
          aerosol(jc,iso4) = aer_so4_mo1 (jc) + ( aer_so4_mo2 (jc) - aer_so4_mo1 (jc) ) * time_weight
          aerosol(jc,idu)  = aer_dust_mo1(jc) + ( aer_dust_mo2(jc) - aer_dust_mo1(jc) ) * time_weight
        ENDDO
        !$ACC END PARALLEL
      CASE(1) ! Simple prognostic for dust
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO jc = istart, iend
          aerosol(jc,iss)  = aer_ss_mo1  (jc) + ( aer_ss_mo2  (jc) - aer_ss_mo1  (jc) ) * time_weight
          aerosol(jc,iorg) = aer_org_mo1 (jc) + ( aer_org_mo2 (jc) - aer_org_mo1 (jc) ) * time_weight
          aerosol(jc,ibc)  = aer_bc_mo1  (jc) + ( aer_bc_mo2  (jc) - aer_bc_mo1  (jc) ) * time_weight
          aerosol(jc,iso4) = aer_so4_mo1 (jc) + ( aer_so4_mo2 (jc) - aer_so4_mo1 (jc) ) * time_weight
          aercl_du(jc)     = aer_dust_mo1(jc) + ( aer_dust_mo2(jc) - aer_dust_mo1(jc) ) * time_weight
        ENDDO
        !$ACC END PARALLEL
      CASE(2,3) ! Simple prognostic for all species
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO jc = istart, iend
          aercl_ss(jc)     = aer_ss_mo1  (jc) + ( aer_ss_mo2  (jc) - aer_ss_mo1  (jc) ) * time_weight
          aercl_or(jc)     = aer_org_mo1 (jc) + ( aer_org_mo2 (jc) - aer_org_mo1 (jc) ) * time_weight
          aercl_bc(jc)     = aer_bc_mo1  (jc) + ( aer_bc_mo2  (jc) - aer_bc_mo1  (jc) ) * time_weight
          aercl_su(jc)     = aer_so4_mo1 (jc) + ( aer_so4_mo2 (jc) - aer_so4_mo1 (jc) ) * time_weight
          aercl_du(jc)     = aer_dust_mo1(jc) + ( aer_dust_mo2(jc) - aer_dust_mo1(jc) ) * time_weight
        ENDDO
        !$ACC END PARALLEL
      CASE DEFAULT
        CALL finish(routine,'iprog_aero setting not implemented')
    END SELECT

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 2, nlevp1
      DO jc = istart, iend
        zsign(jc,jk) = pres_ifc(jc,jk) / MAX(pref_aerdis(jc),0.95_wp*pres_ifc(jc,nlevp1))
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    CALL aerdis ( &
      & kbdim  = nproma,      & !in
      & jcs    = istart,      & !in
      & jce    = iend,        & !in
      & klevp1 = nlevp1,      & !in
      & petah  = zsign(1,1),  & !in
      & pvdaes = zvdaes(1,1), & !out
      & pvdael = zvdael(1,1), & !out
      & pvdaeu = zvdaeu(1,1), & !out
      & pvdaed = zvdaed(1,1), & !out
      & lacc = lzacc)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR PRIVATE(jk, zslatq)
    DO jc = istart, iend
      ! top level
      zaeqso (jc) = zvdaes(jc,1) * aerosol(jc,iss)
      zaeqlo (jc) = zvdael(jc,1) * aerosol(jc,iorg)
      zaeqsuo(jc) = zvdael(jc,1) * aerosol(jc,iso4)
      zaequo (jc) = zvdaeu(jc,1) * aerosol(jc,ibc) 
      zaeqdo (jc) = zvdaed(jc,1) * aerosol(jc,idu)
    
      ! tropopause pressure and PBL stability
      jk          = k850(jc)
      zslatq      = SIN(latitude(jc))**2
      zptrop(jc)  = 1.e4_wp + 2.e4_wp*zslatq ! 100 hPa at the equator, 300 hPa at the poles
      zdtdz(jc)   = (temp(jc,jk)-temp(jc,nlev-1))/(-rd/grav*                     &
        &           (temp(jc,jk)+temp(jc,nlev-1))*(pres(jc,jk)-pres(jc,nlev-1))/ &
        &           (pres(jc,jk)+pres(jc,nlev-1)))
      ! latitude-dependence of tropospheric background
      zlatfac(jc) = MAX(0.1_wp, 1._wp-MERGE(zslatq**3, zslatq, latitude(jc) > 0._wp))
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP SEQ
    DO jk = 1,nlev
      !$ACC LOOP GANG VECTOR PRIVATE(zaeqsn, zaeqln, zaeqsun, zaequn, zaeqdn, zstrfac, zpblfac)
      DO jc = istart, iend
        zaeqsn  = zvdaes(jc,jk+1) * aerosol(jc,iss)
        zaeqln  = zvdael(jc,jk+1) * aerosol(jc,iorg)
        zaeqsun = zvdael(jc,jk+1) * aerosol(jc,iso4)
        zaequn  = zvdaeu(jc,jk+1) * aerosol(jc,ibc)
        zaeqdn  = zvdaed(jc,jk+1) * aerosol(jc,idu)

        ! stratosphere factor: 1 in stratosphere, 0 in troposphere, width of transition zone 0.1*p_TP
        zstrfac = MIN(1._wp,MAX(0._wp,10._wp*(zptrop(jc)-pres(jc,jk))/zptrop(jc)))
        ! PBL stability factor; enhance organic, sulfate and black carbon aerosol for stable stratification
        zpblfac = 1._wp + MIN(1.5_wp,1.e2_wp*MAX(0._wp, zdtdz(jc) + grav/cpd))

        zaeq1(jc,jk) = (1._wp-zstrfac)*MAX(zpblfac*(zaeqln-zaeqlo(jc)), ztrbga*zlatfac(jc)*dpres_mc(jc,jk))
        zaeq2(jc,jk) = (1._wp-zstrfac)*(zaeqsn-zaeqso(jc))
        zaeq3(jc,jk) = (1._wp-zstrfac)*(zaeqdn-zaeqdo(jc))
        zaeq4(jc,jk) = (1._wp-zstrfac)*zpblfac*(zaequn-zaequo(jc))
        zaeq5(jc,jk) = (1._wp-zstrfac)*zpblfac*(zaeqsun-zaeqsuo(jc)) + zstrfac*zstbga*dpres_mc(jc,jk)

        zaeqso(jc)  = zaeqsn
        zaeqlo(jc)  = zaeqln
        zaeqsuo(jc) = zaeqsun
        zaequo(jc)  = zaequn
        zaeqdo(jc)  = zaeqdn

      ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$ACC END DATA

  END SUBROUTINE nwp_aerosol_tegen
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  SUBROUTINE nwp_cpl_aero_gscp_conv(istart, iend, nlev, pres_sfc, pres, acdnc, cloud_num, lacc)
  INTEGER, INTENT(in)                 :: &
    &  istart, iend, nlev                  !< loop start and end indices (nproma, vertical)
  REAL(wp), INTENT(in)                :: &
    &  pres_sfc(:), pres(:,:)              !< Surface and atmospheric pressure
  REAL(wp), INTENT(inout)                :: &
    &  acdnc(:,:),                       & !< cloud droplet number concentration
    &  cloud_num(:)                        !< cloud droplet number concentration
  LOGICAL, INTENT(in), OPTIONAL       :: &
    &  lacc                                !< If true, use openacc
  ! Local variables
  REAL(wp)                            :: &
    &  wfac, ncn_bg
  INTEGER                             :: &
    &  jc, jk                              !< Loop indices
  LOGICAL                             :: &
    &  lzacc                               !< non-optional version of lacc

  CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(wfac, ncn_bg)
    DO jk = 1,nlev
      DO jc = istart, iend
        wfac         = MAX(1._wp,MIN(8._wp,0.8_wp*pres_sfc(jc)/pres(jc,jk)))**2
        ncn_bg       = MIN(cloud_num(jc),50.e6_wp)
        acdnc(jc,jk) = (ncn_bg+(cloud_num(jc)-ncn_bg)*(EXP(1._wp-wfac)))
      END DO
    END DO
    !$ACC END PARALLEL

  END SUBROUTINE nwp_cpl_aero_gscp_conv
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  SUBROUTINE get_time_intp_weights(mtime_datetime, imo1 , imo2, time_weight)
    TYPE(datetime), POINTER, INTENT(in) :: &
      &  mtime_datetime                      !< Current datetime
    INTEGER, INTENT(out)                :: &
      &  imo1 , imo2                         !< Month indices for temporal interpolation
    REAL(wp), INTENT(out)               :: &
      &  time_weight
    ! Local variables
    TYPE(datetime), POINTER             :: &
      &  current_time_hours
    TYPE(t_time_interpolation_weights)  :: &
      &  current_time_interpolation_weights

    current_time_hours => newDatetime(mtime_datetime)
    current_time_hours%time%minute = 0
    current_time_hours%time%second = 0
    current_time_hours%time%ms = 0

    current_time_interpolation_weights = calculate_time_interpolation_weights(current_time_hours)

    imo1        = current_time_interpolation_weights%month1
    imo2        = current_time_interpolation_weights%month2
    time_weight = current_time_interpolation_weights%weight2

    CALL deallocateDatetime(current_time_hours)
    
  END SUBROUTINE get_time_intp_weights

  !---------------------------------------------------------------------------------------
  SUBROUTINE nwp_aerosol_cleanup(od_lw, od_sw, ssa_sw, g_sw)
    CHARACTER(len=*), PARAMETER :: &
      &  routine = modname//':nwp_aerosol_cleanup'

    REAL(wp), ALLOCATABLE, INTENT(inout) :: &
      &  od_lw(:,:,:,:),       & !< Longwave optical thickness
      &  od_sw(:,:,:,:),       & !< Shortwave optical thickness
      &  ssa_sw(:,:,:,:),      & !< Shortwave asymmetry factor
      &  g_sw(:,:,:,:)           !< Shortwave single scattering albedo
    ! Local variables
    INTEGER :: istat

    IF( ALLOCATED(od_lw) ) THEN
      DEALLOCATE(od_lw, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of od_lw failed.')
    ENDIF
    IF( ALLOCATED(od_sw) ) THEN
      DEALLOCATE(od_sw, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of od_sw failed.')
    ENDIF
    IF( ALLOCATED(ssa_sw) ) THEN
      DEALLOCATE(ssa_sw, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of ssa_sw failed.')
    ENDIF
    IF( ALLOCATED(g_sw) ) THEN
      DEALLOCATE(g_sw, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of g_sw failed.')
    ENDIF

  END SUBROUTINE nwp_aerosol_cleanup

END MODULE mo_nwp_aerosol

