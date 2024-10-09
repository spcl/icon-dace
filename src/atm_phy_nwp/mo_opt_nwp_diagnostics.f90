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

! Routines for optional diagnostic output variables in NWP
! (formerly located in mo_util_phys)

!NEC$ options "-finline-max-depth=3 -finline-max-function-size=1000"

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_opt_nwp_diagnostics

  USE mo_kind,                  ONLY: vp, wp
  USE mo_parallel_config,       ONLY: nproma, proc0_shift
  USE mo_math_constants,        ONLY: pi
  USE mo_physical_constants,    ONLY: o_m_rdv        , & !! 1 - r_d/r_v &
    &                                 rdv,             & !! r_d / r_v
    &                                 vtmpc1,          &
    &                                 grav,            &
    &                                 tmelt, earth_radius, &
    &                                 alvdcp, rd_o_cpd, &
    &                                 rhoh2o, rhoice, K_w_0, K_i_0
  USE mo_lookup_tables_constants,ONLY: c1es, c3les, c4les
  USE mo_nh_diagnose_pres_temp, ONLY: calc_qsum
  USE mo_opt_nwp_reflectivity,  ONLY: compute_field_dbz_1mom, compute_field_dbz_2mom
  USE mo_exception,             ONLY: finish, message, warning
  USE mo_fortran_tools,         ONLY: assign_if_present, set_acc_host_or_device, assert_acc_host_only, &
    &                                 assert_acc_device_only, init, copy
  USE mo_impl_constants,        ONLY: min_rlcell_int, min_rledge_int, &
    &                                 min_rlcell, grf_bdywidth_c
  USE mo_impl_constants_grf,    ONLY: grf_bdyintp_start_c,  &
    &                                 grf_ovlparea_start_c, grf_fbk_start_c
  USE mo_model_domain,          ONLY: t_patch
  USE mo_nonhydro_types,        ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,         ONLY: t_nwp_phy_diag
  USE mo_run_config,            ONLY: iqv, iqc, iqi, iqr, iqs, iqg, iqni, &
       &                              iqh, iqnc, iqnr, iqns, iqng, iqnh, iqgl, iqhl, msg_level, dtime
  USE mo_loopindices,           ONLY: get_indices_c, get_indices_e
  USE mo_atm_phy_nwp_config,    ONLY: atm_phy_nwp_config
  USE mo_nonhydrostatic_config, ONLY: kstart_moist
  USE mo_io_config,             ONLY: echotop_meta, &
       &                              wdur_min_hailcast
  USE mo_lnd_nwp_config,        ONLY: nlev_soil, dzsoil
  USE mo_nwp_lnd_types,         ONLY: t_lnd_diag
  USE mo_ext_data_types,        ONLY: t_external_data
  USE mo_intp_data_strc,        ONLY: t_int_state
  USE mo_nwp_sfc_interp,        ONLY: wsoil2smi
  USE mo_icon_interpolation_scalar,                     &
    &                           ONLY: edges2cells_scalar, cells2edges_scalar, &
    &                                 cells2verts_scalar, verts2edges_scalar, &
    &                                 cells2edges_scalar
  USE mo_math_gradients,        ONLY: grad_fd_norm, grad_fd_tang
  USE mo_intp_rbf,              ONLY: rbf_vec_interpol_edge
  USE mo_sync,                  ONLY: sync_patch_array, SYNC_C
  USE mo_grf_intp_data_strc,    ONLY: p_grf_state_local_parent
  USE mo_communication,         ONLY: exchange_data
  USE mo_grid_config,           ONLY: l_limited_area
  USE mo_mpi,                   ONLY: my_process_is_mpi_workroot, get_my_mpi_work_id
  USE mo_timer,                 ONLY: timer_start, timer_stop, timers_level
  USE mo_diag_hailcast,         ONLY: hailstone_driver
  USE mo_util_phys,             ONLY: inversion_height_index  
  USE mo_nwp_tuning_config,     ONLY: tune_dursun_scaling
  USE microphysics_1mom_schemes,ONLY: get_cloud_number, get_snow_temperature
#ifdef HAVE_RADARFWO
  USE radar_data_mie,             ONLY: ldebug_dbz, T0C_emvorado => T0C_fwo
  USE radar_interface,            ONLY: initialize_tmax_atomic_1mom, &
    &                                   initialize_tmax_atomic_2mom, &
    &                                   initialize_tmin_atomic_1mom, &
    &                                   initialize_tmin_atomic_2mom, &
    &                                   init_1mom_types, init_2mom_types      
  USE radar_mie_iface_cosmo_1mom, ONLY: radar_mie_1mom_vec, &
    &                                   radar_rayleigh_oguchi_1mom_vec
  USE radar_mie_iface_cosmo_2mom, ONLY: radar_mie_2mom_vec, &
    &                                   radar_rayleigh_oguchi_2mom_vec
  USE mo_synradar_config,         ONLY: synradar_meta, ydir_mielookup_read, ydir_mielookup_write
  USE mo_mpi,                     ONLY: get_my_mpi_work_comm_size
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: calsnowlmt
  PUBLIC :: compute_field_omega
  PUBLIC :: compute_field_pv
  PUBLIC :: compute_field_sdi
  PUBLIC :: compute_field_lpi
  PUBLIC :: maximize_field_lpi
  PUBLIC :: compute_field_ceiling
  PUBLIC :: compute_field_hbas_sc
  PUBLIC :: compute_field_htop_sc
  PUBLIC :: compute_field_twater
  PUBLIC :: compute_field_q_sedim
  PUBLIC :: compute_field_dursun
  PUBLIC :: compute_field_tcond_max
  PUBLIC :: compute_field_uh_max
  PUBLIC :: compute_field_vorw_ctmax
  PUBLIC :: compute_field_w_ctmax
  PUBLIC :: compute_field_smi
  PUBLIC :: cal_cloudtop
  PUBLIC :: cal_cape_cin
  PUBLIC :: cal_cape_cin_mu
  PUBLIC :: cal_cape_cin_mu_COSMO
  PUBLIC :: cal_si_sli_swiss
  PUBLIC :: compute_field_dbz3d_lin
  PUBLIC :: compute_field_dbzcmax
  PUBLIC :: compute_field_dbz850
  PUBLIC :: compute_field_dbzlmx
  PUBLIC :: maximize_field_dbzctmax
  PUBLIC :: compute_field_echotop
  PUBLIC :: compute_field_echotopinm
  PUBLIC :: compute_field_wshear
  PUBLIC :: compute_field_lapserate  
  PUBLIC :: compute_field_mconv
  PUBLIC :: compute_field_srh
  PUBLIC :: compute_field_visibility
  PUBLIC :: compute_field_inversion_height
  PUBLIC :: compute_updraft_duration
  PUBLIC :: compute_hail_statistics
  
  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_opt_nwp_diagnostics'

CONTAINS



  !------------------------------------------------------------------------------
  !>
  !! Description:
  !!   This subroutine calculates height of the snowfall limit (snowlmt).
  !!
  !! Method:
  !!   In a first step the wet bulb temperature is derived from pres, t and qv.
  !!   In a second step the snowfall limit is evaluated from 8000m down to the
  !!   the lowest model level (ke) and linearly interpolated to the height where
  !!   the wet bulb temperature is >= wbl (=+1.3C after P. Haechler, MeteoSwiss).
  !!   A flag (-999) is set to indicate that no snowlmt was found.
  !!
  SUBROUTINE calsnowlmt ( snowlmt, temp, pres, qv, hhl, hhlr, istart, iend, wbl, lacc)

    ! Parameter list:

    INTEGER, INTENT (IN)     ::  &
      istart, iend           ! loop start/end indices

    REAL(wp), INTENT (INOUT)   ::  &
      snowlmt(:)    ! height of the snowfall limit in m above sea level

    REAL(wp), INTENT (IN)    ::  &
      temp  (:,:), & ! temperature
      pres  (:,:), & ! pressure at full levels
      qv    (:,:), & ! specific humidity
      hhl   (:,:), & ! height of model half levels
      hhlr  (:)      ! height of model half levels resp. sea level

    REAL (wp), INTENT (IN)    ::  &
      wbl               ! (empirical) wet bulb temperature at snowfall limit (1.3C)

    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc
    !------------------------------------------------------------------------------
    ! Local variables

    INTEGER ::     i, k, ktopmin, nlev

    LOGICAL                  ::    &
      lfound(SIZE(temp,1))     ! Logical flag : =.TRUE when wet bulb temp corresponding to
                     !                  parameter "wbl" is found

    REAL (wp)       ::    &
      za = 0.78588481_wp,      & ! local storage
      zb = 7.567_wp,           &
      zc = 2066.92605_wp,      &
      zd = 33.45_wp,           &
      ze = 0.622_wp,           &
      zf = 0.378_wp,           &
      zg = 0.5_wp,             &
      zh = 0.6_wp,             &
      zi = 700._wp,            &
      zl = 0.1_wp,             &
      zm = 6400._wp,           &
      zn = 11.564_wp,          &
      zo = 1742._wp,           &
      td,tl,tp,                &
      zp,                      &  ! pressure in hPa
      ppp,                     &  ! pressure in dPa
      deltat,zt,               &
      ep,const,                &
      zh_bot, zh_top,          &
      zdt

    REAL(wp) :: wetblb(SIZE(temp,1),SIZE(temp,2))  ! wet-bulb temperature in Celsius

    LOGICAL :: lzacc ! non-optional version of lacc

  !------------------------------------------------------------------------------
    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA CREATE(lfound, wetblb) PRESENT(snowlmt, temp, pres, qv, hhl, hhlr) IF(lzacc)

    ! Begin subroutine calsnowlmt

    ! number of vertical full levels
    nlev = SIZE(temp,2)

    ! Set the uppermost model level for the occurence of a wet bulb temperature (wbl)
    ! to about 8000m above surface
    ktopmin = nlev+2
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) REDUCTION(MIN: ktopmin) IF(lzacc)
    DO k = nlev+1, 1, -1
      IF ( hhlr(k) < 8000.0_wp ) THEN
        ktopmin = k
      ENDIF
    ENDDO
    !$ACC WAIT(1)
    if( ktopmin>nlev+1 ) ktopmin = 2

    ! Initialize the definition mask and the output array snowlmt
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    DO i = 1, SIZE(temp,1)
      lfound (i) = .FALSE.
      snowlmt(i) = -999.0_wp
    ENDDO

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP SEQ
    DO k = ktopmin, nlev
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zp, ep, CONST, td, tl, tp, ppp) &
      !$ACC   PRIVATE(deltat, zt)
      DO i = istart, iend
        zp     = (pres(i,k))/100._wp     ! in hPa
        ep     = MAX(1.0E-10_wp,qv(i,k))*zp /      &
                 (ze + zf*MAX(1.0E-10_wp,qv(i,k)))
        ep     = MAX(ep,1.0E-10_wp)
        CONST  = LOG10(ep) - za
        td     = (zd*CONST-zc) / (CONST-zb)              ! in Kelvin
        ! Wet bulb temperature after Egger/Joss
        tl     = (temp(i,k) - tmelt) *10._wp
        tp     = (td-tmelt) *10._wp
        ppp    = zp * 10._wp
        deltat = tl-tp
        zt     = tp + zg*deltat*(zh-tp/zi)
        wetblb(i,k) = zl * ( tp +                      & ! in Celsius
                      (deltat / (1._wp + zm*EXP(zn*zt/(zo+zt))/ppp)))

        IF ( wetblb(i,k) >= wbl ) THEN
          ! definition of snowlmt can be made in this column
          lfound (i) = .TRUE.
        ENDIF
      ENDDO
    ENDDO

    !$ACC LOOP SEQ
    DO k = ktopmin+1, nlev
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zh_bot, zh_top, zdt)
      DO i = istart, iend
        IF ( lfound(i) .AND. wetblb(i,k) >= wbl ) THEN
          ! definition of snowlmt is now made once
          lfound (i) = .FALSE.
          zh_bot     = 0.5_wp * ( hhl(i,k) + hhl(i,k+1) )
          zh_top     = 0.5_wp * ( hhl(i,k) + hhl(i,k-1) )
          zdt        = ( wbl - wetblb(i,k) ) /                 &
                       ( wetblb(i,k-1) - wetblb(i,k) )
          snowlmt(i) = zh_bot + (zh_top-zh_bot)*zdt
        ENDIF
      ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$ACC WAIT(1)

    !$ACC END DATA

  END SUBROUTINE calsnowlmt






  !> computation of vertical velocity (dp/dt)
  !!
  SUBROUTINE compute_field_omega(ptr_patch, p_prog, out_var, &
    &                            opt_slev, opt_elev, opt_rlstart, opt_rlend, lacc)

    TYPE(t_patch)        , INTENT(IN)    :: ptr_patch              !< patch on which computation is performed
    TYPE(t_nh_prog)      , INTENT(IN)    :: p_prog                 !< nonhydrostatic state
    REAL(wp)             , INTENT(INOUT) :: out_var(:,:,:)         !< output variable, dim: (nproma,nlev,nblks_c)
    INTEGER, INTENT(IN), OPTIONAL        :: opt_slev, opt_elev     !< optional vertical start/end level
    INTEGER, INTENT(IN), OPTIONAL        :: opt_rlstart, opt_rlend !< start and end values of refin_ctrl flag
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc


    ! local
    REAL(wp):: w_avg               ! vertical velocity averaged to full level
    INTEGER :: slev, elev          ! vertical start and end index
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_nchdom
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc, jk, jb  
    LOGICAL :: lzacc ! non-optional version of lacc

    ! default values
    CALL set_acc_host_or_device(lzacc, lacc)
    slev     = 1
    elev     = UBOUND(out_var,2)
    rl_start = 2
    rl_end   = min_rlcell_int-1
    ! check optional arguments
    CALL assign_if_present(slev,     opt_slev)
    CALL assign_if_present(elev,     opt_elev)
    CALL assign_if_present(rl_start, opt_rlstart)
    CALL assign_if_present(rl_end,   opt_rlend)

    ! values for the blocking
    i_nchdom   = MAX(1,ptr_patch%n_childdom)
    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL    
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx,w_avg), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        i_startidx, i_endidx, rl_start, rl_end)
      
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(w_avg)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif
          ! half level to full level interpolation
          w_avg = 0.5_wp * (p_prog%w(jc,jk,jb) + p_prog%w(jc,jk+1,jb))

          out_var(jc,jk,jb) = -p_prog%rho(jc,jk,jb)*grav*w_avg

        ENDDO
      ENDDO
      !$ACC END PARALLEL

    ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_omega



  !> computation of soil mositure index (smi)
  !!
  !! Conversion of soil moisture into soil moisture index
  !! smi = (soil moisture - wilting point) / (field capacity - wilting point)
  !!
  SUBROUTINE compute_field_smi(ptr_patch, diag_lnd, ext_data, out_var, &
    &                            opt_rlstart, opt_rlend, lacc)

    TYPE(t_patch)        , INTENT(IN)    :: ptr_patch              !< patch on which computation is performed
    TYPE(t_lnd_diag)     , INTENT(IN)    :: diag_lnd               !< nwp diag land state
    TYPE(t_external_data), INTENT(IN)    :: ext_data               !< ext_data state
    REAL(wp)             , INTENT(INOUT) :: out_var(:,:,:)         !< output variable, dim: (nproma,nlev,nblks_c)
    INTEGER, INTENT(IN), OPTIONAL        :: opt_rlstart, opt_rlend !< start and end values of refin_ctrl flag
    LOGICAL              , INTENT(IN), OPTIONAL :: lacc            !< If true, use openacc

    ! local
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk
    INTEGER :: jc, jk, jb, ic  
    INTEGER :: i_count
    INTEGER :: ierr, ierr_wsoil2smi
    LOGICAL :: lzacc ! non-optional version of lacc

    CALL set_acc_host_or_device(lzacc, lacc)

   !$ACC DATA &
   !$ACC   PRESENT(ptr_patch, diag_lnd%w_so, dzsoil, ext_data%atm%list_land%ncount) &
   !$ACC   PRESENT(ext_data%atm%soiltyp, out_var) IF(lzacc)
    ! default values
    rl_start = 2
    rl_end   = min_rlcell_int-1
    ! check optional arguments
    CALL assign_if_present(rl_start, opt_rlstart)
    CALL assign_if_present(rl_end,   opt_rlend)

    ! values for the blocking
    i_startblk = ptr_patch%cells%start_block(rl_start)
    i_endblk   = ptr_patch%cells%end_block(rl_end)

!$OMP PARALLEL    
!$OMP DO PRIVATE(jc,jk,jb,ic,i_count,ierr,ierr_wsoil2smi), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      ierr = 0

      ! loop over target (ICON) land points only
      i_count = ext_data%atm%list_land%ncount(jb)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG PRIVATE(jc)
      DO ic = 1, i_count
        jc = ext_data%atm%list_land%idx(ic,jb)
        !$ACC LOOP VECTOR PRIVATE(ierr_wsoil2smi) REDUCTION(MIN: ierr)
        DO jk = 1, nlev_soil-1

          CALL wsoil2smi(wsoil   = diag_lnd%w_so(jc,jk,jb),     & !in
            &            dzsoil  = dzsoil(jk),                  & !in
            &            soiltyp = ext_data%atm%soiltyp(jc,jb), & !in
            &            smi     = out_var(jc,jk,jb),           & !out
            &            ierr    = ierr_wsoil2smi               ) !out
          !
          ierr = MIN(ierr, ierr_wsoil2smi)
        ENDDO
        ! assume no-gradient condition for soil moisture reservoir layer
        out_var(jc,nlev_soil,jb) = out_var(jc,nlev_soil-1,jb)
      ENDDO
      !$ACC END PARALLEL
      IF (ierr < 0) THEN
        CALL finish("compute_field_smi", "Landpoint has invalid soiltype (sea water or sea ice)")
      ENDIF

    ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  !$ACC END DATA
  END SUBROUTINE compute_field_smi



  !> Computation of potential vorticity
  !! The full 3D-Ertel PV is calculated at the edges and interpolated to cells.
  !! The shallow atmosphere approximations are used.
  !!
  !! Implemented by Tobias Selz, LMU
  
  SUBROUTINE compute_field_pv(p_patch, p_int_state, p_metrics, p_prog, p_diag, out_var, lacc )

    TYPE(t_patch)        , INTENT(INOUT) :: p_patch              !< patch on which computation is performed
    TYPE(t_int_state)    , INTENT(IN)    :: p_int_state
    TYPE(t_nh_metrics)   , INTENT(IN)    :: p_metrics
    TYPE(t_nh_prog)      , INTENT(IN)    :: p_prog                 !< nonhydrostatic state
    TYPE(t_nh_diag)      , INTENT(IN)    :: p_diag
    REAL(wp)             , INTENT(INOUT) :: out_var(:,:,:)         !< output variable, dim: (nproma,nlev,nblks_c)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    !Local variables
    !Indices
    INTEGER  :: slev, elev, rl_start, rl_end, i_nchdom,     &
      &         i_startblk, i_endblk, i_startidx, i_endidx, &
      &         jc, je, jk, jb, ivd1, ivd2
    LOGICAL :: lzacc ! non-optional version of lacc

    REAL(wp) :: vdfac

    !temporary fields
    REAL(wp) :: pv_ef    (nproma,p_patch%nlev  ,p_patch%nblks_e),  &
                vt       (nproma,p_patch%nlev  ,p_patch%nblks_e),  &
                theta_cf (nproma,p_patch%nlev  ,p_patch%nblks_c),  &
                theta_vf (nproma,p_patch%nlev  ,p_patch%nblks_v),  &
                theta_ef (nproma,p_patch%nlev  ,p_patch%nblks_e),  &
                w_vh     (nproma,p_patch%nlev+1,p_patch%nblks_v),  & 
                w_eh     (nproma,p_patch%nlev+1,p_patch%nblks_e),  &
                ddtw_eh  (nproma,p_patch%nlev+1,p_patch%nblks_e),  &
                ddnw_eh  (nproma,p_patch%nlev+1,p_patch%nblks_e),  &
                ddtth_ef (nproma,p_patch%nlev  ,p_patch%nblks_e),  &
                ddnth_ef (nproma,p_patch%nlev  ,p_patch%nblks_e),  &
                vor_ef   (nproma,p_patch%nlev  ,p_patch%nblks_e)
                
    !Pointers to metric terms
    REAL(vp), POINTER :: ddnz(:,:,:), ddtz(:,:,:), gamma(:,:,:)

    CALL set_acc_host_or_device(lzacc, lacc)

    ddnz  => p_metrics%ddxn_z_full
    ddtz  => p_metrics%ddxt_z_full
    gamma => p_metrics%ddqz_z_full_e


    ! Index bounds
    slev     = 1
    elev     = UBOUND(out_var,2)
    rl_start = 1
    rl_end   = min_rlcell

    ! values for the blocking
    i_nchdom   = MAX(1,p_patch%n_childdom)
    i_startblk = p_patch%cells%start_blk (rl_start,1)
    i_endblk   = p_patch%cells%end_blk   (rl_end,i_nchdom)

    !$ACC DATA CREATE(pv_ef, vt, theta_cf, theta_vf, theta_ef, w_vh, w_eh, ddtw_eh, ddnw_eh) &
    !$ACC   CREATE(ddtth_ef, ddnth_ef, vor_ef) &
    !$ACC   PRESENT(ddnz, ddtz, gamma, p_prog, p_diag, p_patch, out_var) IF(lzacc)
    
!$OMP PARALLEL    
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    !compute theta on cells
    DO jb = i_startblk, i_endblk
    
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx

          theta_cf(jc,jk,jb) = p_diag%temp(jc,jk,jb) / p_prog%exner(jc,jk,jb)
          
        ENDDO
      ENDDO
      !$ACC END PARALLEL

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! synchronize theta
    CALL sync_patch_array(SYNC_C, p_patch, theta_cf)

    !Get vt at edges (p_diag%vt is not up to date)
    CALL rbf_vec_interpol_edge( p_prog%vn, p_patch, p_int_state, vt)
    
    !Interpolate theta to vertices
    CALL cells2verts_scalar( theta_cf, p_patch, p_int_state%cells_aw_verts, theta_vf )
    
    !Interpolate theta to edges
    CALL cells2edges_scalar( theta_cf, p_patch, p_int_state%c_lin_e, theta_ef, lacc=lzacc )
    
    !Interpolate w to vertices
    CALL cells2verts_scalar( p_prog%w, p_patch, p_int_state%cells_aw_verts, w_vh )
    
    !Interpolate w to edges
    CALL cells2edges_scalar( p_prog%w, p_patch, p_int_state%c_lin_e, w_eh, lacc=lzacc )
    
    !Interpolate vorticity to edges
    CALL verts2edges_scalar( p_diag%omega_z, p_patch, p_int_state%v_1o2_e, vor_ef )
    
    !Calculate horizontal derivatives of w and theta
    CALL grad_fd_norm ( p_prog%w, p_patch, ddnw_eh  )
    CALL grad_fd_tang ( w_vh,     p_patch, ddtw_eh  )
    CALL grad_fd_norm ( theta_cf, p_patch, ddnth_ef )
    CALL grad_fd_tang ( theta_vf, p_patch, ddtth_ef )
    
    !Recompute loop indices for edges
    rl_start   = 3
    rl_end     = min_rledge_int-1
    i_startblk = p_patch%edges%start_blk (rl_start,1)
    i_endblk   = p_patch%edges%end_blk   (rl_end,i_nchdom)
    
!$OMP PARALLEL    
!$OMP DO PRIVATE(je,jk,jb,i_startidx,i_endidx,ivd1,ivd2,vdfac), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk
      
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(ivd1, ivd2, vdfac)
      DO jk = slev, elev
        DO je = i_startidx, i_endidx

          !Get indices for vertical derivatives of full level variables
          IF ( jk == slev ) THEN
            ivd1=slev
            ivd2=slev+1
            vdfac=1_wp
          ELSE IF ( jk == elev ) THEN
            ivd1=elev-1
            ivd2=elev
            vdfac=1_wp
          ELSE
            ivd1=jk-1
            ivd2=jk+1
            vdfac=2_wp
          END IF
          
          !Ertel-PV calculation on edges
          pv_ef(je,jk,jb) =                                                                                     &
            &     (   0.5_wp*(ddnw_eh(je,jk,jb)+ddnw_eh(je,jk+1,jb))                                            &
            &       + ddnz(je,jk,jb)/gamma(je,jk,jb)*(w_eh(je,jk+1,jb)-w_eh(je,jk,jb))                          &
            &       + (p_prog%vn(je,ivd2,jb)-p_prog%vn(je,ivd1,jb))/vdfac/gamma(je,jk,jb)                       &
            &       - p_prog%vn(je,jk,jb)/earth_radius                                                          &
            &     )                                                                                             &
            &   * (   ddtth_ef(je,jk,jb)                                                                        &
            &       + ddtz(je,jk,jb)/gamma(je,jk,jb) * (theta_ef(je,ivd2,jb)-theta_ef(je,ivd1,jb))/vdfac        &
            &     )                                                                                             &
            &   + ( - (vt(je,ivd2,jb)-vt(je,ivd1,jb))/vdfac/gamma(je,jk,jb)                                     &
            &       - 0.5_wp*(ddtw_eh(je,jk,jb)+ddtw_eh(je,jk+1,jb))                                            &
            &       - ddtz(je,jk,jb)/gamma(je,jk,jb) * (w_eh(je,jk+1,jb)-w_eh(je,jk,jb))                        &
            &       + vt(je,jk,jb)/earth_radius                                                                 &
            &      )                                                                                            &
            &   * (   ddnth_ef(je,jk,jb)                                                                        &
            &       + ddnz(je,jk,jb)/gamma(je,jk,jb) * (theta_ef(je,ivd2,jb)-theta_ef(je,ivd1,jb))/vdfac        &
            &     )                                                                                             &
            &   + (   vor_ef(je,jk,jb)                                                                          &
            &       + ddtz(je,jk,jb)/gamma(je,jk,jb) * (p_prog%vn(je,ivd2,jb)-p_prog%vn(je,ivd1,jb))/vdfac      &
            &       - ddnz(je,jk,jb)/gamma(je,jk,jb) * (vt(je,ivd2,jb)-vt(je,ivd1,jb))/vdfac                    &
            &       + p_patch%edges%f_e(je,jb)                                                                  &
            &     )                                                                                             &
            &   * ( -(theta_ef(je,ivd2,jb)-theta_ef(je,ivd1,jb))/vdfac/gamma(je,jk,jb) )
                   
        ENDDO
      ENDDO
      !$ACC END PARALLEL

    ENDDO 
    !$ACC WAIT(1)
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    
    !Interpolate to cells
    CALL edges2cells_scalar( pv_ef, p_patch, p_int_state%e_bln_c_s, out_var, opt_rlstart=2 )
    

    rl_start = 2
    rl_end   = min_rlcell_int

    ! values for the blocking
    i_nchdom   = MAX(1,p_patch%n_childdom)
    i_startblk = p_patch%cells%start_blk (rl_start,1)
    i_endblk   = p_patch%cells%end_blk   (rl_end,i_nchdom)

    !Normalize with density
    !
!$OMP PARALLEL    
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk
    
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
          out_var(jc,jk,jb) = out_var(jc,jk,jb) / p_prog%rho(jc,jk,jb)
        ENDDO
      ENDDO
      !$ACC END PARALLEL

    ENDDO  ! jb
    !$ACC WAIT(1)
!$OMP END DO NOWAIT
!$OMP END PARALLEL    

    !$ACC END DATA

  END SUBROUTINE compute_field_pv


  !>
  !! compute_field_sdi
  !!
  !! Description:
  !!   calculation of the supercell detection indices (SDI1, SDI2)
  !!
  !! Method:
  !!   defined in:
  !!   Wicker L, J. Kain, S. Weiss and D. Bright, A Brief Description of the
  !!          Supercell Detection Index, (available from
  !!   http://www.spc.noaa.gov/exper/Spring_2005/SDI-docs.pdf)
  !!
  SUBROUTINE compute_field_sdi( ptr_patch, jg, ptr_patch_local_parent, p_int,    &
                                p_metrics, p_prog, p_diag,                 &
                                sdi_2, lacc)

    IMPLICIT NONE

    TYPE(t_patch),      INTENT(IN)    :: ptr_patch           !< patch on which computation is performed
    INTEGER,            INTENT(IN)    :: jg    ! domain ID of main grid
    TYPE(t_patch), TARGET, INTENT(IN) :: ptr_patch_local_parent    !< parent grid for larger exchange halo
    TYPE(t_int_state),  INTENT(IN)    :: p_int   ! for reduced grid

    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_prog),    INTENT(IN)    :: p_prog
    TYPE(t_nh_diag),    INTENT(IN)    :: p_diag
 
    REAL(wp),           INTENT(OUT)   :: sdi_2(:,:)    !< output variable, dim: (nproma,nblks_c)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    INTEGER  :: i_rlstart,  i_rlend
    INTEGER  :: i_startblk, i_endblk
    INTEGER  :: i_startidx, i_endidx

    REAL(wp) :: z_min, z_max

    REAL(wp) :: hsurf( nproma )
    REAL(wp) :: delta_z
    REAL(wp) :: vol( nproma )
    REAL(wp) :: vol_inv
    REAL(wp) :: area_norm

    TYPE(t_patch),  POINTER  :: p_pp

    REAL(wp) :: w_c

    ! vertical averages:
    REAL(wp) :: w_vmean        ( nproma, ptr_patch%nblks_c)
    REAL(wp) :: zeta_vmean     ( nproma, ptr_patch%nblks_c)
    REAL(wp) :: w_w_vmean      ( nproma, ptr_patch%nblks_c)
    REAL(wp) :: zeta_zeta_vmean( nproma, ptr_patch%nblks_c)
    REAL(wp) :: w_zeta_vmean   ( nproma, ptr_patch%nblks_c)

    ! volume averages:
    REAL(wp) :: w_mean        ( nproma, ptr_patch%nblks_c)
    REAL(wp) :: zeta_mean     ( nproma, ptr_patch%nblks_c)
    REAL(wp) :: w_w_mean      ( nproma, ptr_patch%nblks_c)
    REAL(wp) :: zeta_zeta_mean( nproma, ptr_patch%nblks_c)
    REAL(wp) :: w_zeta_mean   ( nproma, ptr_patch%nblks_c)

    INTEGER, PARAMETER :: idx_w         = 1
    INTEGER, PARAMETER :: idx_zeta      = 2
    INTEGER, PARAMETER :: idx_w_w       = 3
    INTEGER, PARAMETER :: idx_zeta_zeta = 4
    INTEGER, PARAMETER :: idx_w_zeta    = 5

    ! the corresponding fields on the parent grid:
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: p_vmean  ! vertical means on the parent grid
    REAL(wp) :: p_mean( nproma, 5)                      ! means on the parent grid

    REAL(wp) :: w2_mean
    REAL(wp) :: zeta2_mean
    REAL(wp) :: helic_mean
    REAL(wp) :: helic_w_corr

    INTEGER :: jb, jc, jk
    INTEGER :: jb2, jc2
    INTEGER :: i, l, idx
    INTEGER :: ist
    INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk

    REAL(wp), POINTER :: p_fbkwgt(:,:,:)

    INTEGER :: nblks_c_lp

    REAL(wp) :: EPS = 1.0E-20_wp
    LOGICAL :: lzacc ! non-optional version of lacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA &
    !$ACC   PRESENT(ptr_patch, p_diag%vor, p_int%cell_environ%idx) &
    !$ACC   PRESENT(p_int%cell_environ%blk, p_int%cell_environ%area_norm) &
    !$ACC   PRESENT(p_metrics%z_ifc, p_prog%w, sdi_2) &
    !$ACC   CREATE(hsurf, vol) &
    !$ACC   CREATE(w_vmean, zeta_vmean, w_w_vmean, zeta_zeta_vmean, w_zeta_vmean) &
    !$ACC   CREATE(w_mean, zeta_mean, w_w_mean, zeta_zeta_mean, w_zeta_mean) &
    !$ACC   CREATE(p_mean) IF(lzacc)


    ! definition of the vertical integration limits:
    z_min = 1500.0  ! in m
    z_max = 5500.0  ! in m

    ! --- to prevent errors at the boundaries, set some field(s) to 0:

    i_rlstart = 1
    i_rlend   = grf_bdywidth_c

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,           &
                          i_startidx, i_endidx, i_rlstart, i_rlend)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        !sdi_1   ( jc, jb ) = 0.0_wp
        sdi_2 ( jc, jb ) = 0.0_wp

        w_vmean        (jc,jb) = 0.0_wp
        zeta_vmean     (jc,jb) = 0.0_wp
        w_w_vmean      (jc,jb) = 0.0_wp
        zeta_zeta_vmean(jc,jb) = 0.0_wp
        w_zeta_vmean   (jc,jb) = 0.0_wp

      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    ! --- calculate vertical averages ---

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,delta_z,vol,vol_inv,hsurf,w_c) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = i_startidx, i_endidx

        w_vmean        (jc,jb) = 0.0_wp
        zeta_vmean     (jc,jb) = 0.0_wp
        w_w_vmean      (jc,jb) = 0.0_wp
        zeta_zeta_vmean(jc,jb) = 0.0_wp
        w_zeta_vmean   (jc,jb) = 0.0_wp

        vol  (jc) = 0.0_wp
        hsurf(jc) = p_metrics%z_ifc(jc,ptr_patch%nlev+1,jb)
      END DO
      !$ACC LOOP SEQ
      DO jk = 1, ptr_patch%nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(delta_z, w_c)
        DO jc = i_startidx, i_endidx

          IF ( ( p_metrics%z_ifc(jc,jk+1,jb) >= hsurf(jc) + z_min ) .AND.      &
            &  ( p_metrics%z_ifc(jc,jk  ,jb) <= hsurf(jc) + z_max ) ) THEN
            ! integrate only between z_min and z_max

            delta_z = p_metrics%z_ifc(jc,jk,jb) - p_metrics%z_ifc(jc,jk+1,jb)

            w_c = 0.5_wp * ( p_prog%w(jc,jk  ,jb)    &
              &            + p_prog%w(jc,jk+1,jb) )

            w_vmean        (jc,jb) = w_vmean        (jc,jb) + delta_z * w_c
            zeta_vmean     (jc,jb) = zeta_vmean     (jc,jb) + delta_z * p_diag%vor(jc,jk,jb)
            w_w_vmean      (jc,jb) = w_w_vmean      (jc,jb) + delta_z * w_c * w_c
            zeta_zeta_vmean(jc,jb) = zeta_zeta_vmean(jc,jb) + delta_z * p_diag%vor(jc,jk,jb) * p_diag%vor(jc,jk,jb)
            w_zeta_vmean   (jc,jb) = w_zeta_vmean   (jc,jb) + delta_z * w_c * p_diag%vor(jc,jk,jb) 

            vol(jc) = vol(jc) + delta_z
          END IF

        END DO
      END DO
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(vol_inv)
      DO jc = i_startidx, i_endidx
        vol_inv = 1.0_wp / vol(jc)
        w_vmean        (jc,jb) = w_vmean        (jc,jb) * vol_inv
        zeta_vmean     (jc,jb) = zeta_vmean     (jc,jb) * vol_inv
        w_w_vmean      (jc,jb) = w_w_vmean      (jc,jb) * vol_inv
        zeta_zeta_vmean(jc,jb) = zeta_zeta_vmean(jc,jb) * vol_inv
        w_zeta_vmean   (jc,jb) = w_zeta_vmean   (jc,jb) * vol_inv
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    ! --- average these values to the parent grid cells
    !
    !     motivation: enhance the horiz. averaging area, which is otherwise 
    !     too strongly limited by the limited halo of the domain decomposition.

    p_pp => ptr_patch_local_parent

    p_fbkwgt => p_grf_state_local_parent(jg)%fbk_wgt_bln

    ! Set pointers to index and coefficient fields for cell-based variables
    iidx => p_pp%cells%child_idx
    iblk => p_pp%cells%child_blk


    nblks_c_lp = p_pp%cells%end_block( min_rlcell )

    ALLOCATE( p_vmean ( nproma, 5, nblks_c_lp), STAT=ist )
    IF ( ist /= 0 ) THEN
      CALL finish( modname//':compute_field_sdi', "allocate failed" )
    END IF
    !$ACC DATA CREATE(p_vmean) IF(lzacc)
    ! first nullify all lateral grid points

    ! Start/End block in the parent domain
    IF (jg == 1 .AND. l_limited_area) THEN
      i_rlstart = grf_bdyintp_start_c
    ELSE
      i_rlstart = grf_ovlparea_start_c
    ENDIF
    i_rlend   = grf_fbk_start_c + 1

    i_startblk = p_pp%cells%start_block( i_rlstart )
    i_endblk   = p_pp%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,idx,jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_pp, jb, i_startblk, i_endblk,           &
                          i_startidx, i_endidx, i_rlstart, i_rlend)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO idx=1, 5
        DO jc = i_startidx, i_endidx
          p_vmean ( jc, idx, jb ) = 0.0_wp
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    
    ! Start/End block in the parent domain
    IF (jg == 1 .AND. l_limited_area) THEN
      i_rlstart = grf_fbk_start_c
    ELSE
      i_rlstart = grf_ovlparea_start_c
    ENDIF
    i_rlend   = min_rlcell_int

    i_startblk = p_pp%cells%start_block( i_rlstart )
    i_endblk   = p_pp%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_pp, jb, i_startblk, i_endblk,           &
                          i_startidx, i_endidx, i_rlstart, i_rlend)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx

        p_vmean(jc, idx_w, jb) =                                                 &
            w_vmean        ( iidx(jc,jb,1), iblk(jc,jb,1) ) * p_fbkwgt(jc,jb,1)  &
          + w_vmean        ( iidx(jc,jb,2), iblk(jc,jb,2) ) * p_fbkwgt(jc,jb,2)  &
          + w_vmean        ( iidx(jc,jb,3), iblk(jc,jb,3) ) * p_fbkwgt(jc,jb,3)  &
          + w_vmean        ( iidx(jc,jb,4), iblk(jc,jb,4) ) * p_fbkwgt(jc,jb,4)

        p_vmean(jc, idx_zeta, jb) =                                              &
            zeta_vmean     ( iidx(jc,jb,1), iblk(jc,jb,1) ) * p_fbkwgt(jc,jb,1)  &
          + zeta_vmean     ( iidx(jc,jb,2), iblk(jc,jb,2) ) * p_fbkwgt(jc,jb,2)  &
          + zeta_vmean     ( iidx(jc,jb,3), iblk(jc,jb,3) ) * p_fbkwgt(jc,jb,3)  &
          + zeta_vmean     ( iidx(jc,jb,4), iblk(jc,jb,4) ) * p_fbkwgt(jc,jb,4)

        p_vmean(jc, idx_w_w, jb) =                                               &
            w_w_vmean      ( iidx(jc,jb,1), iblk(jc,jb,1) ) * p_fbkwgt(jc,jb,1)  &
          + w_w_vmean      ( iidx(jc,jb,2), iblk(jc,jb,2) ) * p_fbkwgt(jc,jb,2)  &
          + w_w_vmean      ( iidx(jc,jb,3), iblk(jc,jb,3) ) * p_fbkwgt(jc,jb,3)  &
          + w_w_vmean      ( iidx(jc,jb,4), iblk(jc,jb,4) ) * p_fbkwgt(jc,jb,4)

        p_vmean(jc, idx_zeta_zeta, jb) =                                         &
            zeta_zeta_vmean( iidx(jc,jb,1), iblk(jc,jb,1) ) * p_fbkwgt(jc,jb,1)  &
          + zeta_zeta_vmean( iidx(jc,jb,2), iblk(jc,jb,2) ) * p_fbkwgt(jc,jb,2)  &
          + zeta_zeta_vmean( iidx(jc,jb,3), iblk(jc,jb,3) ) * p_fbkwgt(jc,jb,3)  &
          + zeta_zeta_vmean( iidx(jc,jb,4), iblk(jc,jb,4) ) * p_fbkwgt(jc,jb,4)

        p_vmean(jc, idx_w_zeta, jb) =                                            &
            w_zeta_vmean   ( iidx(jc,jb,1), iblk(jc,jb,1) ) * p_fbkwgt(jc,jb,1)  &
          + w_zeta_vmean   ( iidx(jc,jb,2), iblk(jc,jb,2) ) * p_fbkwgt(jc,jb,2)  &
          + w_zeta_vmean   ( iidx(jc,jb,3), iblk(jc,jb,3) ) * p_fbkwgt(jc,jb,3)  &
          + w_zeta_vmean   ( iidx(jc,jb,4), iblk(jc,jb,4) ) * p_fbkwgt(jc,jb,4)

       END DO
       !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  
    ! --- Exchange of these fields on the parent grid

    CALL exchange_data(p_pat=p_pp%comm_pat_c, lacc=lzacc, recv=p_vmean )

    ! --- Average over the neighbouring parent grid cells

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,p_mean,  &
!$OMP            l,jc2,jb2,area_norm,i) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_pp, jb, i_startblk, i_endblk,             &
                          i_startidx, i_endidx, i_rlstart, i_rlend)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO idx=1, 5
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = i_startidx, i_endidx
          p_mean( jc, idx ) = 0.0_wp
        END DO
      ENDDO
      !$ACC LOOP SEQ
      DO l=1, p_int%cell_environ%max_nmbr_nghbr_cells
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jc2, jb2, area_norm)
        DO jc = i_startidx, i_endidx

          jc2       = p_int%cell_environ%idx      ( jc, jb, l)
          jb2       = p_int%cell_environ%blk      ( jc, jb, l)
          area_norm = p_int%cell_environ%area_norm( jc, jb, l)

          p_mean( jc, idx_w        ) = p_mean( jc, idx_w        ) + area_norm * p_vmean( jc2, idx_w,         jb2)
          p_mean( jc, idx_zeta     ) = p_mean( jc, idx_zeta     ) + area_norm * p_vmean( jc2, idx_zeta,      jb2)
          p_mean( jc, idx_w_w      ) = p_mean( jc, idx_w_w      ) + area_norm * p_vmean( jc2, idx_w_w,       jb2)
          p_mean( jc, idx_zeta_zeta) = p_mean( jc, idx_zeta_zeta) + area_norm * p_vmean( jc2, idx_zeta_zeta, jb2)
          p_mean( jc, idx_w_zeta   ) = p_mean( jc, idx_w_zeta   ) + area_norm * p_vmean( jc2, idx_w_zeta,    jb2)

        END DO
      END DO

      ! write back to the 4 child grid cells:
      !$ACC LOOP SEQ
      DO i = 1, 4
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = i_startidx, i_endidx
          w_mean        ( iidx(jc,jb,i), iblk(jc,jb,i) ) = p_mean( jc, idx_w )
          zeta_mean     ( iidx(jc,jb,i), iblk(jc,jb,i) ) = p_mean( jc, idx_zeta )
          w_w_mean      ( iidx(jc,jb,i), iblk(jc,jb,i) ) = p_mean( jc, idx_w_w )
          zeta_zeta_mean( iidx(jc,jb,i), iblk(jc,jb,i) ) = p_mean( jc, idx_zeta_zeta )
          w_zeta_mean   ( iidx(jc,jb,i), iblk(jc,jb,i) ) = p_mean( jc, idx_w_zeta )
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL
    !$ACC END DATA ! p_vmean
    DEALLOCATE( p_vmean, STAT=ist ) 
    IF ( ist /= 0 ) THEN
      CALL finish( modname//':compute_field_sdi', "deallocate failed" )
    END IF

    ! --- calculate SDI_1, SDI_2 (now again on the child-grid)

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,helic_mean,w2_mean,zeta2_mean,helic_w_corr) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,      &
                          i_startidx, i_endidx, i_rlstart, i_rlend)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR PRIVATE(helic_mean, zeta2_mean, w2_mean, helic_w_corr)
      DO jc = i_startidx, i_endidx

        helic_mean  = w_zeta_mean(jc,jb)    - w_mean(jc,jb) * zeta_mean(jc,jb)
        zeta2_mean  = zeta_zeta_mean(jc,jb) - zeta_mean(jc,jb) * zeta_mean(jc,jb)
        w2_mean     = w_w_mean(jc,jb)       - w_mean(jc,jb) * w_mean(jc,jb)

        IF ( ( w2_mean > EPS ) .AND. ( zeta2_mean > EPS ) ) THEN

          helic_w_corr = helic_mean / SQRT( w2_mean * zeta2_mean )

          !sdi_1(jc,jb) = helic_w_corr * zeta_vmean(jc,jb)

          IF ( w_vmean(jc,jb) > 0.0_wp ) THEN
            sdi_2(jc,jb) = helic_w_corr * ABS( zeta_vmean(jc,jb) )
          ELSE
            sdi_2(jc,jb) = 0.0_wp
          END IF

        ELSE
         !sdi_1(jc,jb) = 0.0_wp
          sdi_2(jc,jb) = 0.0_wp
        END IF

      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  !$ACC END DATA
  END SUBROUTINE compute_field_sdi


  !>
  !! Calculate the Lightning Potential Index (LPI)
  !!
  !! Literature
  !!       - B. Lynn, Y. Yair, 2010: Prediction of lightning flash density with the WRF model,
  !!           Adv. Geosci., 23, 11-16
  !! adapted from the COSMO-implementation by Uli Blahak.
  !!
  SUBROUTINE compute_field_lpi( ptr_patch, jg, ptr_patch_local_parent, p_int,   &
                                p_metrics, p_prog, p_prog_rcf, p_diag,          &
                                lpi, lacc )

    IMPLICIT NONE

    TYPE(t_patch),      INTENT(IN)    :: ptr_patch         !< patch on which computation is performed
    INTEGER,            INTENT(IN)    :: jg                ! domain ID of main grid
    TYPE(t_patch), TARGET, INTENT(IN) :: ptr_patch_local_parent  !< parent grid for larger exchange halo
    TYPE(t_int_state),  INTENT(IN)    :: p_int
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_prog),    INTENT(IN)    :: p_prog, p_prog_rcf
    TYPE(t_nh_diag),    INTENT(INOUT) :: p_diag

    REAL(wp),           INTENT(OUT)   :: lpi(:,:)          !< output variable, dim: (nproma,nblks_c)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    TYPE(t_patch),  POINTER  :: p_pp

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx

    REAL(wp) :: delta_z
    REAL(wp) :: q_liqu, q_solid, epsw
    REAL(wp) :: lpi_incr
    REAL(wp) :: w_c
    REAL(wp) :: Tmelt_m_20K

    REAL(wp) :: w_updraft_crit
    REAL(wp) :: qg_low_limit, qg_upp_limit, qg_scale_inv, q_i, q_s, q_g
    REAL(wp) :: w_thresh, w_frac_tresh

    REAL(wp) :: vol( nproma )

    INTEGER :: jc, jk, jb
    INTEGER :: jc2, jb2
    INTEGER :: l, i
    INTEGER :: ist

    INTEGER :: nmbr_w  ( nproma, ptr_patch%nblks_c )
    !INTEGER :: nmbr_all( nproma, ptr_patch%nblks_c )  ! only necessary for 'variant 1 of updraft in environment crit.'

    INTEGER, ALLOCATABLE :: p_nmbr_w(:,:)
    INTEGER  :: p_nmbr_w_sum  (nproma)
    INTEGER  :: p_nmbr_all_sum(nproma)
    REAL(wp) :: p_frac_w      (nproma)

    REAL(wp) :: frac_w( nproma, ptr_patch%nblks_c )

    INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk
    INTEGER :: nblks_c_lp
    LOGICAL :: lzacc ! non-optional version of lacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA CREATE(vol, nmbr_w, p_nmbr_w_sum, p_nmbr_all_sum, p_frac_w, frac_w) IF(lzacc)

    !$ACC DATA PRESENT(p_int%cell_environ%idx, p_int%cell_environ%blk, p_int%cell_environ%area_norm) &
    !$ACC   PRESENT(ptr_patch, p_diag, p_metrics, p_prog, p_prog_rcf, atm_phy_nwp_config(jg)) &
    !$ACC   PRESENT(p_diag%temp, p_metrics%ddqz_z_full, p_prog%w, p_prog_rcf%tracer) &
    !$ACC   PRESENT(lpi, vol, nmbr_w, p_nmbr_w_sum, p_nmbr_all_sum, p_frac_w, frac_w) &
    !$ACC   PRESENT(kstart_moist(jg)) IF(lzacc)

    IF (.not.atm_phy_nwp_config(jg)%lhave_graupel) THEN
      CALL finish( modname//'compute_field_lpi',  &
        &     "no graupel available! Either switch off LPI output or change the microphysics scheme" )
    END IF

    Tmelt_m_20K = Tmelt - 20.0_wp

    ! --- Thresholds, limits, ... for several criteria: ---

    ! for the 'updraft criterion' in the LPI integral:
    w_updraft_crit = 0.5_wp  ! in m/s; threshold for w

    ! for the 'Graupel-criterion' in the LPI integral:
    qg_low_limit = 0.0002_wp  ! in kg/kg; lower limit for the 'Graupel-criterion'
    qg_upp_limit = 0.001_wp   ! in kg/kg; upper limit for the 'Graupel-criterion'

    qg_scale_inv = 1.0_wp / ( qg_upp_limit - qg_low_limit )


!!$ UB: need to check this threshold for ICON-D2!
    ! for the 'updraft in environment'-criterion:
    w_thresh = 0.5_wp  ! in m/s; threshold for w
                       ! see Lynn, Yair (2010)
                       ! however, Uli Blahak determined 1.1 m/s for COSMO

    w_frac_tresh = 0.5_wp   ! in 100%; threshold for the fraction of points with 'w>w_thresh'
                            ! to identify growing thunderstorms
                            ! Lynn, Yair (2010) propose 0.5 (however, the whole criterion is unclear)


    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    ! nullify every grid point (lateral boundary, too)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jb=1,ptr_patch%nblks_c
      DO jc=1,nproma
        lpi(jc,jb) = 0.0_wp
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    ! --- calculation of the LPI integral ---

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,vol,delta_z,w_c,q_liqu,q_i,q_s,q_g,q_solid,epsw,lpi_incr), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        vol     (jc)     = 0.0_wp
        nmbr_w  (jc, jb) = 0
       !nmbr_all(jc, jb) = 0  ! only necessary for 'variant 1 of updraft in environment crit.'
      END DO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO jk = kstart_moist(jg), ptr_patch%nlev
        !$ACC LOOP GANG VECTOR PRIVATE(delta_z, w_c, q_liqu, q_i, q_s, q_g, epsw, q_solid, lpi_incr)
        DO jc = i_startidx, i_endidx

          IF (      (p_diag%temp(jc,jk,jb) <= Tmelt)     &
              .AND. (p_diag%temp(jc,jk,jb) >= Tmelt_m_20K) ) THEN

            delta_z = p_metrics%ddqz_z_full(jc,jk,jb)
            w_c     = 0.5_wp * ( p_prog%w(jc,jk,jb) + p_prog%w(jc,jk+1,jb) )

            q_liqu  = p_prog_rcf%tracer(jc,jk,jb,iqc)   &
                    + p_prog_rcf%tracer(jc,jk,jb,iqr)

            IF (atm_phy_nwp_config(jg)%l2moment) THEN
              ! sometimes during IAU slightly negative values of the hydrometeors were encountered, which
              ! lead to a crash in the sqrt() below. Therefore we clip all hydrometeors at 0.0 for security:
              q_i = MAX(p_prog_rcf%tracer(jc,jk,jb,iqi), 0.0_wp)
              q_s = MAX(p_prog_rcf%tracer(jc,jk,jb,iqs), 0.0_wp)
              q_g = MAX(p_prog_rcf%tracer(jc,jk,jb,iqg) + p_prog_rcf%tracer(jc,jk,jb,iqh), 0.0_wp)
            ELSE
              q_i = MAX(p_prog_rcf%tracer(jc,jk,jb,iqi), 0.0_wp)
              q_s = MAX(p_prog_rcf%tracer(jc,jk,jb,iqs), 0.0_wp)
              q_g = MAX(p_prog_rcf%tracer(jc,jk,jb,iqg), 0.0_wp)
            END IF

            IF (atm_phy_nwp_config(jg)%inwp_gscp == 7) THEN
              q_g = q_g + p_prog_rcf%tracer(jc,jk,jb,iqgl) + p_prog_rcf%tracer(jc,jk,jb,iqhl)
            END IF
            
            q_solid = q_g *                                                 &
                 &    ( SQRT( q_i * q_g  ) / MAX( q_i + q_g, 1.0e-20_wp) +  &
                 &      SQRT( q_s * q_g  ) / MAX( q_s + q_g, 1.0e-20_wp) )

            epsw = 2.0_wp * SQRT( q_liqu * q_solid ) / MAX( q_liqu + q_solid, 1.0e-20_wp)

            ! 'updraft-criterion' in the LPI integral
            IF ( w_c >= w_updraft_crit ) THEN
              ! only (strong enough) updrafts should be counted, no downdrafts
              lpi_incr = w_c * w_c * epsw * delta_z

              ! additional 'Graupel-criterion' in the LPI integral
              lpi_incr = lpi_incr * MAX( MIN( (q_g - qg_low_limit) * qg_scale_inv, 1.0_wp ), 0.0_wp )

            ELSE
              lpi_incr = 0.0_wp
            END IF

            lpi(jc,jb) = lpi(jc,jb) + lpi_incr
            vol(jc)    = vol(jc)    + delta_z

            ! sum up for the later test on growth of thunderstorm ('updraft in environment crit.'):

            ! 'variant 1 of updraft in environment crit.' (3D-criterion, MB)
            !IF (w_c >= w_thresh) THEN
            !  nmbr_w(jc,jb) = nmbr_w(jc,jb) + 1
            !END IF
            !nmbr_all(jc,jb) = nmbr_all(jc,jb) + 1

            ! 'variant 2 of updraft in environment crit.' (2D-criterion, UB)
            IF (w_c >= w_thresh) THEN
              nmbr_w(jc,jb) = 1
            END IF

          END IF

        END DO
      END DO
      !$ACC END PARALLEL

      ! normalization
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        IF ( vol(jc) > 1.0e-30_wp ) THEN
          lpi(jc,jb) = lpi(jc,jb) / vol(jc)
        ELSE
          lpi(jc,jb) = 0.0_wp
        END IF
      END DO
      !$ACC END PARALLEL

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    ! --- 'updraft in environment criterion' ---


    ! spatial filtering, i.e. 
    ! test for growth phase of the thunderstorm
    ! (remark: currently only 'variant 2' (2D-criterion) is implemented here)

    ! --- average nmbr_w to the parent grid cells

    p_pp => ptr_patch_local_parent

    ! Set pointers to index and coefficient fields for cell-based variables
    iidx => p_pp%cells%child_idx
    iblk => p_pp%cells%child_blk

    !$ACC DATA PRESENT(iidx, iblk) IF(lzacc)

    nblks_c_lp = p_pp%cells%end_block( min_rlcell )

    ALLOCATE( p_nmbr_w ( nproma, nblks_c_lp), STAT=ist )
    IF ( ist /= 0 ) THEN
      CALL finish( modname//':compute_field_lpi', "allocate failed" )
    END IF

    !$ACC DATA CREATE(p_nmbr_w) IF(lzacc)

    ! first nullify all lateral grid points

    ! Start/End block in the parent domain
    IF (jg == 1 .AND. l_limited_area) THEN
      i_rlstart = grf_bdyintp_start_c
    ELSE
      i_rlstart = grf_ovlparea_start_c
    ENDIF
    i_rlend   = grf_fbk_start_c + 1

    i_startblk = p_pp%cells%start_block( i_rlstart )
    i_endblk   = p_pp%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_pp, jb, i_startblk, i_endblk,           &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        p_nmbr_w( jc, jb) = 0.0_wp
      END DO
      !$ACC END PARALLEL

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    ! Start/End block in the parent domain
    IF (jg == 1 .AND. l_limited_area) THEN
      i_rlstart = grf_fbk_start_c
    ELSE
      i_rlstart = grf_ovlparea_start_c
    ENDIF
    i_rlend   = min_rlcell_int

    i_startblk = p_pp%cells%start_block( i_rlstart )
    i_endblk   = p_pp%cells%end_block  ( i_rlend   )


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_pp, jb, i_startblk, i_endblk,           &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        p_nmbr_w(jc, jb) =                            &
              nmbr_w( iidx(jc,jb,1), iblk(jc,jb,1) )  &
            + nmbr_w( iidx(jc,jb,2), iblk(jc,jb,2) )  &
            + nmbr_w( iidx(jc,jb,3), iblk(jc,jb,3) )  &
            + nmbr_w( iidx(jc,jb,4), iblk(jc,jb,4) ) 
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! --- Exchange of these fields on the parent grid
    CALL exchange_data(p_pat=p_pp%comm_pat_c, lacc=lzacc, recv=p_nmbr_w )

    ! --- Average over the neighbouring parent grid cells

    ! consistency check
    IF ( p_int%cell_environ%max_nmbr_iter /= 1 ) THEN
      CALL finish( modname//':compute_field_lpi', "cell_environ is not built with the right number of iterations" )
    END IF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,  &
!$OMP            p_nmbr_w_sum, p_nmbr_all_sum,  &
!$OMP            l,jc2,jb2,p_frac_w,i) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_pp, jb, i_startblk, i_endblk,           &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        p_nmbr_w_sum  (jc)  = 0
        p_nmbr_all_sum(jc) = 0
      END DO
      !$ACC END PARALLEL

      !$ACC PARALLEL ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO l=1, p_int%cell_environ%max_nmbr_nghbr_cells
        
        !$ACC LOOP GANG VECTOR PRIVATE(jc2, jb2)
        DO jc = i_startidx, i_endidx

          jc2 = p_int%cell_environ%idx( jc, jb, l)
          jb2 = p_int%cell_environ%blk( jc, jb, l)

          IF ( p_int%cell_environ%area_norm( jc, jb, l) > 1.0e-7_wp ) THEN
            p_nmbr_w_sum  (jc) = p_nmbr_w_sum  (jc) + p_nmbr_w(jc2, jb2)
            p_nmbr_all_sum(jc) = p_nmbr_all_sum(jc) + 4
          END IF
        END DO
      END DO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        p_frac_w(jc) = DBLE( p_nmbr_w_sum(jc) ) / DBLE( p_nmbr_all_sum(jc) )
      END DO
      !$ACC END PARALLEL

      ! write back to the 4 child grid cells:
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO i = 1, 4
        DO jc = i_startidx, i_endidx
          frac_w( iidx(jc,jb,i), iblk(jc,jb,i) ) = p_frac_w(jc)
        END DO
      END DO
      !$ACC END PARALLEL

    END DO
    !$ACC WAIT(1)
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    !$ACC END DATA ! p_nmbr_w

    DEALLOCATE( p_nmbr_w, STAT=ist )
    IF ( ist /= 0 ) THEN
      CALL finish( modname//':compute_field_lpi', "deallocate failed" )
    END IF

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        ! finally, this is the 'updraft in environment criterion':
        IF ( frac_w(jc,jb) < w_frac_tresh ) THEN
          lpi(jc, jb) = 0.0_wp
        END IF
      END DO
      !$ACC END PARALLEL

    END DO
    !$ACC WAIT(1)
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  !$ACC END DATA ! iidx,iblk

  !$ACC END DATA ! header
  !$ACC END DATA ! create data
  END SUBROUTINE compute_field_lpi


  !>
  !! Do a maximization step for the calculation of the LPI_MAX.
  !!
  !!
  SUBROUTINE maximize_field_lpi( ptr_patch, jg, ptr_patch_local_parent, p_int,   &
                                p_metrics, p_prog, p_prog_rcf, p_diag,           &
                                lpi_max, lacc )

    IMPLICIT NONE

    TYPE(t_patch),      INTENT(IN)    :: ptr_patch         !< patch on which computation is performed
    INTEGER,            INTENT(IN)    :: jg                ! domain ID of main grid
    TYPE(t_patch), TARGET, INTENT(IN) :: ptr_patch_local_parent  !< parent grid for larger exchange halo
    TYPE(t_int_state),  INTENT(IN)    :: p_int
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_prog),    INTENT(IN)    :: p_prog, p_prog_rcf
    TYPE(t_nh_diag),    INTENT(INOUT) :: p_diag

    REAL(wp), INTENT(INOUT) :: lpi_max( nproma, ptr_patch%nblks_c )
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jc
    LOGICAL :: lzacc ! non-optional version of lacc

    REAL(wp) :: lpi( nproma, ptr_patch%nblks_c )
    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA CREATE(lpi) PRESENT(lpi_max) IF(lzacc)
    CALL compute_field_lpi( ptr_patch, jg, ptr_patch_local_parent, p_int,   &
                            p_metrics, p_prog, p_prog_rcf, p_diag,          &
                            lpi, lacc=lzacc )

    ! Maximize

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int
    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE

    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        lpi_max(jc,jb) = MAX( lpi_max(jc,jb), lpi(jc,jb) )
      END DO
      !$ACC END PARALLEL
    END DO
    !$ACC WAIT(1)
    !$ACC END DATA
!$OMP END PARALLEL

  END SUBROUTINE maximize_field_lpi


  !>
  !! Calculate the ceiling height
  !! = height above MSL, for which cloud coverage > 4/8
  !!
  SUBROUTINE compute_field_ceiling( ptr_patch, jg,    &
                                p_metrics, prm_diag,  &
                                ceiling_height, lacc )

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch     !< patch on which computation is performed
    INTEGER,              INTENT(IN)  :: jg            ! domain ID of main grid
    TYPE(t_nh_metrics),   INTENT(IN)  :: p_metrics
    TYPE(t_nwp_phy_diag), INTENT(IN)  :: prm_diag

    REAL(wp),             INTENT(OUT) :: ceiling_height(:,:)    !< output variable, dim: (nproma,nblks_c)
    LOGICAL, INTENT(IN), OPTIONAL     :: lacc ! If true, use openacc
    ! Local scalar and arrays
    ! -----------------------

    LOGICAL ::  cld_base_found( nproma ) 

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jk, jc
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)
    !$ACC DATA &
    !$ACC   PRESENT(ceiling_height, ptr_patch) &
    !$ACC   PRESENT(p_metrics, prm_diag, kstart_moist(jg)) &
    !$ACC   CREATE(cld_base_found) &
    !$ACC   IF(lzacc)
    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,cld_base_found), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(i_startidx, i_endidx, jb, jg) IF(lzacc)
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = i_startidx, i_endidx
        cld_base_found(jc) = .FALSE.
        ceiling_height(jc,jb) = p_metrics%z_mc(jc,1,jb)  ! arbitrary default value
      END DO
      !$ACC LOOP SEQ
      DO jk = ptr_patch%nlev, kstart_moist(jg), -1
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = i_startidx, i_endidx
          IF ( .NOT.(cld_base_found(jc)) .AND. (prm_diag%clc(jc,jk,jb) > 0.5_wp) ) THEN
            ceiling_height(jc,jb) = p_metrics%z_mc(jc,jk,jb)
            cld_base_found(jc) = .TRUE.
          ENDIF
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ENDDO
    !$ACC WAIT(1)
    !$ACC END DATA
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_ceiling


  !>
  !! Calculate the height of base over MSL from the shallow convection parameterization.
  !!
  !! This subroutine is quite similar to compute_field_htop_sc.
  !!
  SUBROUTINE compute_field_hbas_sc( ptr_patch,        &
                                p_metrics, prm_diag,  &
                                hbas_sc, lacc)

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch     !< patch on which computation is performed
    TYPE(t_nh_metrics),   INTENT(IN)  :: p_metrics
    TYPE(t_nwp_phy_diag), INTENT(IN)  :: prm_diag

    REAL(wp),             INTENT(OUT) :: hbas_sc(:,:)    !< output variable, dim: (nproma,nblks_c)
    LOGICAL,              INTENT(IN), OPTIONAL :: lacc   !< If true, use openacc
    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jc, idx
    

    REAL(wp), PARAMETER :: undefValue = 0.0_wp
    LOGICAL :: lzacc ! non-optional version of lacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA &
    !$ACC   PRESENT(ptr_patch, p_metrics%z_mc, prm_diag%ktype) &
    !$ACC   PRESENT(prm_diag%mbas_con, hbas_sc) IF(lzacc)

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,idx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR PRIVATE(idx)
      DO jc = i_startidx, i_endidx

        IF ( prm_diag%ktype(jc,jb) == 2 ) THEN
          ! the column is identified as shallow convection
          idx = prm_diag%mbas_con( jc, jb)
          hbas_sc(jc,jb) = p_metrics%z_mc( jc, idx, jb)
       ELSE
          hbas_sc(jc,jb) = undefValue
        END IF

      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  !$ACC END DATA
  END SUBROUTINE compute_field_hbas_sc


  !>
  !! Calculate the height of top over MSL from the shallow convection parameterization.
  !!
  !! This subroutine is quite similar to compute_field_hbas_sc.
  !!
  SUBROUTINE compute_field_htop_sc( ptr_patch,        &
                                p_metrics, prm_diag,  &
                                htop_sc, lacc)

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch     !< patch on which computation is performed
    TYPE(t_nh_metrics),   INTENT(IN)  :: p_metrics
    TYPE(t_nwp_phy_diag), INTENT(IN)  :: prm_diag

    REAL(wp),             INTENT(OUT) :: htop_sc(:,:)    !< output variable, dim: (nproma,nblks_c)
    LOGICAL,              INTENT(IN), OPTIONAL :: lacc   !< If true, use openacc
    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jc, idx

    REAL(wp), PARAMETER :: undefValue = 0.0_wp
    LOGICAL :: lzacc ! non-optional version of lacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA &
    !$ACC   PRESENT(ptr_patch, p_metrics%z_mc, prm_diag%ktype) &
    !$ACC   PRESENT(prm_diag%mbas_con, htop_sc) IF(lzacc)

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,idx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR PRIVATE(idx)
      DO jc = i_startidx, i_endidx

        IF ( prm_diag%ktype(jc,jb) == 2 ) THEN
          ! the column is identified as shallow convection
          idx = prm_diag%mtop_con( jc, jb)
          htop_sc(jc,jb) = p_metrics%z_mc( jc, idx, jb)
       ELSE
          htop_sc(jc,jb) = undefValue
        END IF

      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  !$ACC END DATA
  END SUBROUTINE compute_field_htop_sc


  !>
  !! Calculate total column integrated water in kg m-2 (twater)
  !!
  SUBROUTINE compute_field_twater( p_patch, ddqz_z_full, rho, tracer,  &
                                   idx_list_condensate, twater, opt_slev, lacc )

    TYPE(t_patch),     INTENT(IN)  :: p_patch                !< patch on which computation is performed
    REAL(wp),          INTENT(IN)  :: ddqz_z_full(:,:,:)     !< cell height [m]
    REAL(wp),          INTENT(IN)  :: rho(:,:,:)             !< total air density [kg m-3]
    REAL(wp),          INTENT(IN)  :: tracer(:,:,:,:)        !< tracer mass fractions [kg kg-1]
    INTEGER,           INTENT(IN)  :: idx_list_condensate(:) !< IDs of all water tracers
                                                             !< excluding qv
    REAL(wp),          INTENT(OUT) :: twater(:,:)            !< total column integrated water [kg m-2]
                                                             !< dim: (nproma,nblks_c)
    INTEGER, OPTIONAL, INTENT(IN)  :: opt_slev               !< vertical start level
    LOGICAL, OPTIONAL, INTENT(IN)  :: lacc                   !< if true, use OpenACC

    ! Local variables
    !----------------
    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jk, jc
    INTEGER :: jg
    INTEGER :: slev              ! start level
    INTEGER :: slev_moist        ! start level for moisture variables other than qv

    REAL(wp):: q_water( nproma, p_patch%nlev )
    LOGICAL :: lzacc             ! OpenACC flag

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_field_twater'


    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA &
    !$ACC   CREATE(q_water) &
    !$ACC   IF(lzacc)

    jg = p_patch%id

    ! sanity check
    !
    IF (ANY(idx_list_condensate > SIZE(tracer,4))) THEN
      CALL finish( routine, "tracer ID exceeds size of tracer container" )
    END IF

    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    ENDIF
    ! start index for moisture variables other than qv
    slev_moist = MAX(kstart_moist(jg),slev)

    ! without halo or boundary points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = p_patch%cells%start_block( i_rlstart )
    i_endblk   = p_patch%cells%end_block  ( i_rlend   )


!$OMP PARALLEL
    CALL init(twater(:,:), lacc=lzacc)

!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx,q_water), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      ! get total mass fraction of condensates (excluding qv)
      ! calc_sum ensures that q_water is initialized with zero for slev<=jk<=slev_moist
      !
      CALL calc_qsum (tracer(:,:,:,:), q_water(:,:), idx_list_condensate, &
        &             jb, i_startidx, i_endidx, slev, slev_moist, p_patch%nlev)

      ! add contribution by qv
      !
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = slev, p_patch%nlev
        DO jc = i_startidx, i_endidx
          q_water(jc,jk) = q_water(jc,jk) + tracer(jc,jk,jb,iqv)
        END DO
      END DO
      !$ACC END PARALLEL

      ! integrate over column
      !
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO jk = slev, p_patch%nlev
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          twater(jc,jb) = twater(jc,jb)       &
            &            + rho(jc,jk,jb) * q_water(jc,jk) * ddqz_z_full(jc,jk,jb)
        END DO
      END DO
      !$ACC END PARALLEL

    ENDDO  !jb
    !$ACC WAIT(1)
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC END DATA
  END SUBROUTINE compute_field_twater


  !>
  !! Calculate specific content of precipitation particles
  !!
  SUBROUTINE compute_field_q_sedim( ptr_patch, jg, p_prog_rcf, q_sedim, lacc )

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch     !< patch on which computation is performed
    INTEGER,              INTENT(IN)  :: jg            ! domain ID of main grid
    TYPE(t_nh_prog),      INTENT(IN)  :: p_prog_rcf

    REAL(wp),             INTENT(OUT) :: q_sedim(:,:,:)  !< output variable, dim: (nproma,nlev,nblks_c)
    LOGICAL, OPTIONAL,    INTENT(IN)  :: lacc            !< if true, use OpenACC

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jk, jc
    LOGICAL :: lzacc ! non-optional version of lacc

    CALL set_acc_host_or_device(lzacc, lacc)

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )
!$OMP PARALLEL
    CALL init(q_sedim( :, :, 1:i_startblk-1 ), lacc=lzacc)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      ! it is assumed that at least a warm rain Kessler scheme is used, i.e. qr is available too
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO jk = kstart_moist(jg), ptr_patch%nlev
        DO jc = i_startidx, i_endidx
          q_sedim(jc,jk,jb) = p_prog_rcf%tracer(jc,jk,jb,iqr)
        END DO
      END DO

      IF ( ASSOCIATED( p_prog_rcf%tracer_ptr(iqs)%p_3d ) ) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            q_sedim(jc,jk,jb) = q_sedim(jc,jk,jb) + p_prog_rcf%tracer(jc,jk,jb,iqs)
          END DO
        END DO
      END IF

      IF ( atm_phy_nwp_config(jg)%lhave_graupel ) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            q_sedim(jc,jk,jb) = q_sedim(jc,jk,jb) + p_prog_rcf%tracer(jc,jk,jb,iqg)
          END DO
        END DO
      END IF

      IF ( atm_phy_nwp_config(jg)%l2moment ) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            q_sedim(jc,jk,jb) = q_sedim(jc,jk,jb) + p_prog_rcf%tracer(jc,jk,jb,iqh)
          END DO
        END DO
      END IF
      !$ACC END PARALLEL

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_q_sedim


  !>
  !! Calculate the low level (mean over 0-1000 m AGL) moisture convergence,
  !! assuming flat surface for simplicity. Any metrical terms due to
  !! terrain following height coordinates are neglected.
  !!
  !! Adapted from the COSMO-implementation by Uli Blahak.
  !!
  !!
  SUBROUTINE compute_field_mconv( ptr_patch, p_int,   &
                                  p_metrics, p_prog, p_prog_rcf, &
                                  z_low, z_up, mconv )

    ! Input/output variables:
    TYPE(t_patch),      INTENT(IN), TARGET  :: ptr_patch     !< patch on which computation is performed
    TYPE(t_int_state),  INTENT(IN), TARGET  :: p_int
    TYPE(t_nh_metrics), INTENT(IN)          :: p_metrics
    TYPE(t_nh_prog),    INTENT(IN), TARGET  :: p_prog, p_prog_rcf
    REAL(wp),           INTENT(IN)          :: z_low         !< Lower height AGL for mconv averaging [m AGL]
    REAL(wp),           INTENT(IN)          :: z_up          !< Upper height AGL for mconv averaging [m AGL]

    REAL(wp),           INTENT(OUT)         :: mconv(:,:)    !> output variable, dim: (nproma,nblks_c)

    ! Local variables:
    CHARACTER(len=*), PARAMETER :: routine = modname//': compute_field_mconv'
    INTEGER, POINTER            :: ieidx(:,:,:), ieblk(:,:,:)
    REAL(wp)                    :: qv_e(nproma,ptr_patch%nlev,ptr_patch%nblks_e) !< qv interpolated to edges
    REAL(wp), POINTER           :: vn_e(:,:,:), geofac_e(:,:,:)
    REAL(wp)                    :: div_qvv_layer(nproma,ptr_patch%nlev)
    REAL(wp), DIMENSION(nproma) :: div_qvv_mean, p_conv_sum, p_conv_wgt

    INTEGER               :: i_rlstart,  i_rlend
    INTEGER               :: i_startblk, i_endblk
    INTEGER               :: i_startidx, i_endidx
    INTEGER               :: jb, jc, jk, k_start, k_start_vec(nproma), &
                              iex(3), ieb(3), l, jc2, jb2, iter

    REAL(wp)              :: wgt_loc, area_norm
    REAL(wp)              :: mconv_smth(SIZE(mconv,dim=1),SIZE(mconv,dim=2))

    ! Parameters for smoothing filter:
    INTEGER, PARAMETER :: niter_smooth = 1 ! number of successive filter applications

    ! Linear interpolation of qv to cell edges:
    CALL cells2edges_scalar(p_prog_rcf%tracer(:,:,:,iqv), ptr_patch, p_int%c_lin_e, qv_e, lacc=.FALSE.)

    ieidx    => ptr_patch%cells%edge_idx
    ieblk    => ptr_patch%cells%edge_blk
    vn_e     => p_prog%vn
    geofac_e => p_int%geofac_div
    
    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
    CALL init(mconv, 0.0_wp)
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,k_start,k_start_vec, &
!$OMP            div_qvv_layer, &
!$OMP            div_qvv_mean,iex,ieb), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      ! Determine the lowermost height level (highest index) which is
      !  everywhere just above the maximum height considered
      !  for the following vertical averages and integrals,
      !  and use this as starting level for all vertical loops to save computing time.
      ! We do a bottom-up search for the level just below z_up and subtract 1 to get one level above:
      k_start_vec(:) = ptr_patch%nlev
      DO jk = ptr_patch%nlev, 2, -1
        DO jc = i_startidx, i_endidx
          IF (p_metrics%z_ifc(jc,jk,jb) - p_metrics%z_ifc(jc,ptr_patch%nlev+1,jb) < z_up) THEN
            k_start_vec(jc) = jk - 1 
          END IF
        END DO
      END DO
      ! k_start_vec contains the max height index in each column. The overall k_start index
      ! is the smallest among them:
      k_start = MINVAL(k_start_vec(i_startidx:i_endidx))

      div_qvv_layer(:,:) = 0.0_wp
      DO jk = k_start, ptr_patch%nlev
        DO jc = i_startidx, i_endidx

          iex(1:3) = ieidx(jc,jb,1:3)
          ieb(1:3) = ieblk(jc,jb,1:3)

          ! Horizontal divergence by Gauss' theorem:
          ! sum of oriented outward edge-normal qv fluxes throuth the side faces of the cell
          ! divided by the cell volume, assuming flat orography. We neglect any metrical terms
          ! due to the terrain-following vertical coordinates.
          ! The geofac_e for each edge is the edge length times orientation factor (+ or -1) divided by cell area
          div_qvv_layer(jc,jk) = &
               vn_e(iex(1),jk,ieb(1)) * qv_e(iex(1),jk,ieb(1)) * geofac_e(jc,1,jb) + &
               vn_e(iex(2),jk,ieb(2)) * qv_e(iex(2),jk,ieb(2)) * geofac_e(jc,2,jb) + &
               vn_e(iex(3),jk,ieb(3)) * qv_e(iex(3),jk,ieb(3)) * geofac_e(jc,3,jb)

        END DO
      END DO

      ! average div_qvv of z_low-z_up AGL:
      CALL vert_integral_vec_1d ( i_startidx, i_endidx, k_start,    &
           &                      hhl  = p_metrics % z_ifc(:,:,jb), &
           &                      f    = div_qvv_layer(:,:),        &
           &                      zlow = z_low,                     &
           &                      zup  = z_up,                      &
           &                      fint = div_qvv_mean(:),           &
           &                      l_agl= .TRUE.,                    &
           &                      l_calc_mean = .TRUE.,             &
           &                      l_rescale_to_full_thickness = .FALSE. &
           &                      )

      DO jc = i_startidx, i_endidx
        mconv(jc,jb) = -div_qvv_mean(jc)
      END DO
      
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !-------------------------------------------------------------------
    !
    ! Spatial smoothing over neighbouring points
    !
    !-------------------------------------------------------------------


    ! Apply approximate binomial smoother niter_smooth times:
    iterloop: DO iter=1, niter_smooth
      
      ! --- Exchange of mconv for reproducible results:
      CALL sync_patch_array(SYNC_C, ptr_patch, mconv)

      ! --- Weighted average over the neighbouring grid cells:
!$OMP PARALLEL
      CALL init(mconv_smth, 0.0_wp)    
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,  &
!$OMP            p_conv_sum,p_conv_wgt,wgt_loc,area_norm, &
!$OMP            l,jc2,jb2), ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                            i_startidx, i_endidx, i_rlstart, i_rlend)

        DO jc = i_startidx, i_endidx
          p_conv_sum(jc)  = 0.0_wp
          p_conv_wgt(jc)  = 0.0_wp
        END DO

        DO l=1, p_int%cell_environ%max_nmbr_nghbr_cells

          IF (l == 1) THEN
            ! This is the center cell, which gets the most weight in the average
            wgt_loc = 1.0_wp
          ELSE
            ! These are the direct neighbours of order 1, which get a lower weight
            wgt_loc = 0.333_wp
          END IF
          
          DO jc = i_startidx, i_endidx
            
            jc2 = p_int%cell_environ%idx( jc, jb, l)
            jb2 = p_int%cell_environ%blk( jc, jb, l)
            area_norm = p_int%cell_environ%area_norm( jc, jb, l)
            IF ( area_norm > 1.0e-7_wp ) THEN
              p_conv_sum(jc) = p_conv_sum(jc) + mconv(jc2,jb2)*wgt_loc*area_norm
              p_conv_wgt(jc) = p_conv_wgt(jc) + wgt_loc*area_norm
            END IF
          END DO
        END DO
        
        DO jc = i_startidx, i_endidx
          IF (p_conv_wgt(jc) > 1e-20_wp) THEN
            mconv_smth(jc,jb) = p_conv_sum(jc) / p_conv_wgt(jc)
          END IF
        END DO
        
      END DO
!$OMP END DO
      ! Copy back the smoothed field to the output variable:
      CALL copy(mconv_smth, mconv)
!$OMP END PARALLEL


    END DO iterloop
    
  END SUBROUTINE compute_field_mconv
  
  !>
  !! Calculate 
  !!     TCOND_MAX   (total column-integrated condensate, max. during the last hour)
  !! and TCOND10_MAX (total column-integrated condensate above z(T=-10 degC), max. during the last hour)
  !! Here, compute columnwise amximum of these input fields and the newly computed fields.
  !! 
  !! Implementation analogous to those of Uli Blahak in COSMO.
  !!
  SUBROUTINE compute_field_tcond_max( ptr_patch, jg,                      &
                                   p_metrics, p_prog, p_prog_rcf, p_diag, &
                                   flag_tcond_max, flag_tcond10_max,      &
                                   tcond_max, tcond10_max, lacc )

    TYPE(t_patch),      INTENT(IN)    :: ptr_patch         !< patch on which computation is performed
    INTEGER,            INTENT(IN)    :: jg                ! domain ID of main grid
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_prog),    INTENT(IN)    :: p_prog, p_prog_rcf
    TYPE(t_nh_diag),    INTENT(INOUT) :: p_diag

    LOGICAL,            INTENT(IN)    :: flag_tcond_max    ! if true, then calculate tcond_max
    LOGICAL,            INTENT(IN)    :: flag_tcond10_max  ! if true, then calculate tcond10_max

    REAL(wp),           INTENT(INOUT) :: tcond_max  (:,:)    !< output variable, dim: (nproma,nblks_c)
    REAL(wp),           INTENT(INOUT) :: tcond10_max(:,:)    !< output variable, dim: (nproma,nblks_c)
    LOGICAL, OPTIONAL,  INTENT(IN)  :: lacc                  !< if true, use OpenACC

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jk, jc

    REAL(wp) :: q_cond( nproma, ptr_patch%nlev )
    REAL(wp) :: tcond( nproma )

    REAL(wp), PARAMETER :: Tzero_m_10K = 263.15    ! in K

    ! Consistency check:
    IF ( (.NOT. flag_tcond_max) .AND. (.NOT. flag_tcond10_max) ) THEN
      CALL finish( "compute_field_tcond_max", "at least one of the two flags must be set to .TRUE." )
    END IF

    CALL assert_acc_device_only("compute_field_tcond_max", lacc)

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    !$ACC DATA CREATE(q_cond, tcond)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,q_cond,tcond), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      ! it is assumed that at least a warm rain Kessler scheme is used, i.e. qr is available too
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO jk = kstart_moist(jg), ptr_patch%nlev
        DO jc = i_startidx, i_endidx
          q_cond(jc,jk) = p_prog_rcf%tracer(jc,jk,jb,iqc)  &
            &           + p_prog_rcf%tracer(jc,jk,jb,iqr)
        END DO
      END DO

      IF ( ASSOCIATED( p_prog_rcf%tracer_ptr(iqi)%p_3d ) ) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            q_cond(jc,jk) = q_cond(jc,jk) + p_prog_rcf%tracer(jc,jk,jb,iqi)
          END DO
        END DO
      END IF

      IF ( ASSOCIATED( p_prog_rcf%tracer_ptr(iqs)%p_3d ) ) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            q_cond(jc,jk) = q_cond(jc,jk) + p_prog_rcf%tracer(jc,jk,jb,iqs)
          END DO
        END DO
      END IF

      IF ( atm_phy_nwp_config(jg)%lhave_graupel ) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            q_cond(jc,jk) = q_cond(jc,jk) + p_prog_rcf%tracer(jc,jk,jb,iqg)
          END DO
        END DO
      END IF

      IF ( atm_phy_nwp_config(jg)%l2moment ) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            q_cond(jc,jk) = q_cond(jc,jk) + p_prog_rcf%tracer(jc,jk,jb,iqh)
          END DO
        END DO
      END IF
      !$ACC END PARALLEL
      ! calculate vertically integrated mass

      IF ( flag_tcond_max ) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = i_startidx, i_endidx
          tcond(jc) = 0.0_wp
        END DO
        !$ACC LOOP SEQ
        DO jk = kstart_moist(jg), ptr_patch%nlev
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = i_startidx, i_endidx
            tcond(jc) = tcond(jc)       &
              &       + p_prog%rho(jc,jk,jb) * q_cond(jc,jk) * p_metrics%ddqz_z_full(jc,jk,jb)
          END DO
        END DO

        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = i_startidx, i_endidx
          tcond_max( jc,jb ) = MAX( tcond_max( jc,jb ), tcond(jc) )
        END DO
        !$ACC END PARALLEL
      END IF

      IF ( flag_tcond10_max ) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = i_startidx, i_endidx
          tcond(jc) = 0.0_wp
        END DO
        !$ACC LOOP SEQ
        DO jk = kstart_moist(jg), ptr_patch%nlev
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = i_startidx, i_endidx
            IF ( p_diag%temp( jc,jk,jb) <= Tzero_m_10K ) THEN
              tcond(jc) = tcond(jc)       &
                &       + p_prog%rho(jc,jk,jb) * q_cond(jc,jk) * p_metrics%ddqz_z_full(jc,jk,jb)
            END IF
          END DO
        END DO

        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = i_startidx, i_endidx
          tcond10_max( jc,jb ) = MAX( tcond10_max( jc,jb ), tcond(jc) )
        END DO
        !$ACC END PARALLEL
      END IF

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE compute_field_tcond_max


  !>
  !! Calculate UH_MAX (updraft helicity, max.  during the last hour)
  !! For the definition see: Kain et al. (2008) Wea. Forecasting
  !!
  !! Implementation analogous to those of Uli Blahak in COSMO.
  !!
  SUBROUTINE compute_field_uh_max( ptr_patch,                 &
                                   p_metrics, p_prog, p_diag, &
                                   zmin_in, zmax_in,          &
                                   uh_max, lacc )

    TYPE(t_patch),      INTENT(IN)    :: ptr_patch         !< patch on which computation is performed
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_prog),    INTENT(IN)    :: p_prog
    TYPE(t_nh_diag),    INTENT(IN)    :: p_diag

    REAL(wp),           INTENT(IN)    :: zmin_in        !< lower boundary for vertical integration
    REAL(wp),           INTENT(IN)    :: zmax_in        !< upper boundary for vertical integration

    REAL(wp),           INTENT(INOUT) :: uh_max(:,:)    !< input/output variable, dim: (nproma,nblks_c)
    LOGICAL, OPTIONAL,  INTENT(IN)  :: lacc             !< if true, use OpenACC

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jk, jc

    REAL(wp) :: zmin( nproma )
    REAL(wp) :: zmax( nproma )
    REAL(wp) :: uhel( nproma )
    REAL(wp) :: w_c

    CALL assert_acc_device_only("compute_field_uh_max", lacc)

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    !$ACC DATA CREATE(zmin, zmax, uhel)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,zmin,zmax,uhel,w_c), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = i_startidx, i_endidx
        IF (zmin_in < 500._wp) THEN
          zmin(jc) = MAX( p_metrics%z_ifc( jc, ptr_patch%nlev+1, jb), zmin_in )
        ELSE
          zmin(jc) = MAX( p_metrics%z_ifc( jc, ptr_patch%nlev+1, jb) + 500.0_wp, zmin_in )
        END IF
        zmax(jc) = zmin(jc) + zmax_in - zmin_in
        uhel(jc) = 0.0_wp
      END DO

      !$ACC LOOP SEQ
      DO jk = 1, ptr_patch%nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(w_c)
        DO jc = i_startidx, i_endidx

          ! Parts of the grid boy are within the bounds, integrate over the exact bounds [zmin,zmax]:
          !  (It also works if the integration layer is so narrow that the bounds are in the same grid box)
          IF ( ( p_metrics%z_ifc( jc, jk+1, jb) <= zmax(jc) ) .AND.     &
            &  ( p_metrics%z_ifc( jc, jk, jb)   >= zmin(jc) ) ) THEN

            w_c = 0.5_wp * ( p_prog%w(jc,jk,jb) + p_prog%w(jc,jk+1,jb) )
            
            ! a simple box-integration in the vertical, but honouring the exact integration bounds zmin, zmax;
            ! only updrafts are counted:
            uhel(jc) = uhel(jc) + MAX( w_c, 0.0_wp) * p_diag%vor(jc,jk,jb) * &
                 ( MIN(p_metrics%z_ifc(jc,jk,jb), zmax(jc)) - MAX(p_metrics%z_ifc( jc, jk+1, jb), zmin(jc)) )
            
          END IF

        END DO
      END DO

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = i_startidx, i_endidx
        uh_max(jc,jb) = MERGE(uhel(jc), uh_max(jc,jb), ABS(uhel(jc)) > ABS(uh_max(jc,jb)) )
      END DO
      !$ACC END PARALLEL

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE compute_field_uh_max


  !>
  !! Calculate VORW_CTMAX (Maximum rotation amplitude during the last hour)
  !!
  !! Implementation analogous to those of Uli Blahak in COSMO.
  !!
  SUBROUTINE compute_field_vorw_ctmax( ptr_patch,          &
                                       p_metrics, p_diag,  &
                                       vorw_ctmax, lacc )

    TYPE(t_patch),      INTENT(IN)    :: ptr_patch         !< patch on which computation is performed
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_diag),    INTENT(IN)    :: p_diag

    REAL(wp),           INTENT(INOUT) :: vorw_ctmax(:,:)  !< input/output variable, dim: (nproma,nblks_c)
    LOGICAL, OPTIONAL,  INTENT(IN)    :: lacc             !< if true, use OpenACC

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jk, jc

    REAL(wp) :: zmin( nproma )
    REAL(wp) :: zmax( nproma )
    REAL(wp) :: vort( nproma )

    CALL assert_acc_device_only("compute_field_vorw_ctmax", lacc)

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    !$ACC DATA CREATE(zmin, zmax, vort)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,zmin,zmax,vort), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = i_startidx, i_endidx
        zmin(jc) = p_metrics%z_ifc( jc, ptr_patch%nlev+1, jb)
        zmax(jc) = MAX( 3000.0_wp,  zmin(jc) + 1500.0_wp )
        vort(jc) = 0.0_wp
      END DO

      DO jk = 1, ptr_patch%nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = i_startidx, i_endidx

          ! Parts of the grid box are within the bounds, integrate over the exact bounds [zmin,zmax]:
          !  (It also works if the integration layer is so narrow that the bounds are in the same grid box)
          IF ( ( p_metrics%z_ifc( jc, jk+1, jb) <= zmax(jc) ) .AND.     &
            &  ( p_metrics%z_ifc( jc, jk, jb)   >= zmin(jc) ) ) THEN

            ! simple box-integration in the vertical, but honouring the exact integration bounds zmin, zmax:
            vort(jc) = vort(jc) + p_diag%vor(jc,jk,jb) * &
                 ( MIN(p_metrics%z_ifc(jc,jk,jb), zmax(jc)) - MAX(p_metrics%z_ifc( jc, jk+1, jb), zmin(jc)) )
            
          END IF

        END DO
      END DO

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = i_startidx, i_endidx
        vort(jc) = vort(jc) / (zmax(jc) - zmin(jc))
        vorw_ctmax(jc,jb) = MERGE(vort(jc), vorw_ctmax(jc,jb), ABS(vort(jc)) > ABS(vorw_ctmax(jc,jb)) )
      END DO
      !$ACC END PARALLEL

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE compute_field_vorw_ctmax


  !>
  !! Calculate W_CTMAX (Maximum updraft track during the last hour)
  !!
  !! Implementation analogous to those of Uli Blahak in COSMO.
  !!
  SUBROUTINE compute_field_w_ctmax( ptr_patch,             &
                                    p_metrics, p_prog,     &
                                    w_ctmax, lacc )

    TYPE(t_patch),      INTENT(IN)    :: ptr_patch         !< patch on which computation is performed
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_prog),    INTENT(IN)    :: p_prog

    REAL(wp),           INTENT(INOUT) :: w_ctmax(:,:)    !< input/output variable, dim: (nproma,nblks_c)
    LOGICAL, OPTIONAL,  INTENT(IN)    :: lacc            !< if true, use OpenACC

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jk, jc

    CALL assert_acc_device_only("compute_field_w_ctmax", lacc)

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk = 1, ptr_patch%nlev
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx

          IF ( p_metrics%z_mc( jc, jk, jb) <= 10000.0_wp ) THEN
            w_ctmax(jc,jb) = MAX( w_ctmax(jc,jb), p_prog%w(jc,jk,jb) )
          END IF

        END DO
      END DO
      !$ACC END PARALLEL

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_w_ctmax


  !------------------------------------------------------------------------------
  !
  !>
  !! Description:
  !!  Computation of Convective Available Potential Energy CAPE_ML,
  !!  Convective Inhibition CIN_ML based on parcel theory.
  !!  This subroutine is based on COSMO code.
  !!        Helmut Frank
  !! 
  !! Input:  
  !!         - Temperature, specific humidity and pressure of environment
  !!
  !! Output: 
  !!         - cape_ml/cin_ml: CAPE/CIN based on a parcel with thermodynamical 
  !!                           properties of the lowest mean layer in the PBL (50hPa)
  !!         - cape_3km/cin_3km: CAPE/CIN based on a parcel with thermodynamical 
  !!                           properties of the lowest mean layer in the PBL (50hPa)
  !!                           with end of ascent at 3 km HAG.
  !!         - lcl_ml/lfc_ml: Lifting Condensation Level/Level of Free Convection
  !!                           based on a parcel with thermodynamical properties of the lowest
  !!                           mean layer in the PBL (50hPa)(ABOVE GROUND level) 
  !!
  !!----------------------------------------------------------------------------
  
  SUBROUTINE cal_cape_cin ( i_startidx, i_endidx, kmoist, te, qve, prs, hhl, & ! in
                            cape_ml, cin_ml, cape_3km, lcl_ml, lfc_ml,       & ! out 
                            lacc )                                             ! in (optional)

    ! Input data
    !----------- 
    INTEGER, INTENT (IN) ::  &
         i_startidx, i_endidx,  &  ! start and end indices of loops in horizontal patch
         kmoist                    ! start index for moist processes

    REAL    (wp),    INTENT (IN) ::  &
         te  (:,:),         & ! environment temperature
         qve (:,:),         & ! environment specific humidity
         prs (:,:),         & ! full level pressure
         hhl (:,:)            ! height of half levels

    ! Output data
    !------------ 
    REAL (wp), INTENT (OUT) :: &
         cape_ml   (:),  & ! mixed layer CAPE_ML
         cin_ml    (:)     ! mixed layer CIN_ML

    REAL (wp), INTENT (OUT), OPTIONAL :: &
         cape_3km  (:),  & ! mixed layer CAPE_ML, with endpoint 3km
         lcl_ml    (:),  & ! mixed layer Lifting Condensation Level ABOVE GROUND level
         lfc_ml    (:)     ! mixed layer Level of Free Convection ABOVE GROUND level

    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    ! Local scalars and automatic arrays
    !-----------------------------------

    REAL(wp) :: &
         qvp_start(SIZE(te,1)), & ! parcel initial specific humidity in mixed layer
         tp_start (SIZE(te,1)), & ! parcel initial potential temperature in mixed layer
         te_start (SIZE(te,1)), & ! parcel initial temperature at center height of mixed layer
         hl_l,hl_u                ! height of full levels in order to find the closest model level

    INTEGER :: &
         jc, k, nlev,             &
         kstart(SIZE(te,1)),     & ! Model level corresponding to start height of parcel
         klcl  (SIZE(te,1)),     & ! Indices for Lifting Condensation Level LCL,
         klfc  (SIZE(te,1)),     & ! Level of Free Convection LFC and
         kel   (SIZE(te,1)),     & ! Equilibrium level EL
         k_ml  (SIZE(te,1))        ! Index for calculation of mixed layer averages 
                                   ! (potential temperature, moisture)

#ifndef _OPENACC
    LOGICAL :: lexit(SIZE(te,1))
#endif

    LOGICAL :: &
         lzacc,         & ! non-optional version of lacc
         l3km,          & ! true if cape_3km is requested
         llev             ! true if either lfc_ml or lcl_ml is requested

    ! Local parameters:
    !------------------
    
    REAL (wp), PARAMETER :: p0 = 1.e5_wp   ! reference pressure for calculation of potential temperature
    REAL (wp), PARAMETER :: missing_value  = -999.9_wp   ! Missing value for CIN (if no LFC/CAPE was found),

    ! Depth of mixed surface layer: 50hPa following Huntrieser, 1997.
    ! Other frequently used value is 100hPa.
    REAL (wp), PARAMETER :: ml_depth = 5000._wp


    !------------------------------------------------------------------------------
    ! 
    ! A well mixed near surface layer is assumed (its depth is specified with 
    ! parameter ml_depth) Potential temperature and specific humidity are constant
    ! in this layer, they are calculated as arithmetical means of the corresponding
    ! variables of the environment (model) profile. The parcel starts from a level 
    ! approximately in the middle of this well mixed layer, with the average spec. 
    ! humidity and potential temperature as start values. 
    !
    !------------------------------------------------------------------------------
    CALL set_acc_host_or_device(lzacc, lacc)

    IF (PRESENT(cape_3km)) THEN
      l3km = .TRUE.
    ELSE
      l3km = .FALSE.
    ENDIF

    IF (PRESENT(lfc_ml) .OR. PRESENT(lcl_ml)) THEN
      llev = .TRUE.
    ELSE
      llev = .FALSE.
    ENDIF

    nlev = SIZE( te,2)

    !$ACC DATA &
    !$ACC   PRESENT(te, qve, prs, hhl) &
    !$ACC   PRESENT(cape_ml, cin_ml, cape_3km, lcl_ml, lfc_ml) &
    !$ACC   CREATE(qvp_start, tp_start, te_start) &
    !$ACC   CREATE(kstart, klcl, klfc, k_ml) &
    !$ACC   IF(lzacc)

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(i_startidx, i_endidx, nlev, kmoist) IF(lzacc)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jc = i_startidx, i_endidx
      k_ml  (jc)  = nlev  ! index used to step through the well mixed layer
      kstart(jc)  = nlev  ! index of model level corresponding to average 
      klfc  (jc)  = nlev  !
      klcl  (jc)  = nlev  !
      ! mixed layer pressure
      qvp_start(jc) = 0.0_wp ! specific humidities in well mixed layer
      tp_start (jc) = 0.0_wp ! potential temperatures in well mixed layer
      ! outputs
      cape_ml  (jc) = 0.0_wp 
      cin_ml   (jc) = missing_value
#ifndef _OPENACC
      lexit(jc)     = .FALSE.
#endif
    ENDDO

    
    ! now calculate the mixed layer average potential temperature and 
    ! specific humidity
    !$ACC LOOP SEQ
    DO k = nlev, kmoist, -1
#ifndef _OPENACC
      IF (ALL(lexit(i_startidx:i_endidx))) EXIT
#endif
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = i_startidx, i_endidx
        IF ( prs(jc,k) > (prs(jc,nlev) - ml_depth)) THEN
          qvp_start(jc) = qvp_start(jc) + qve(jc,k)
          tp_start (jc) = tp_start (jc) + te (jc,k)*(p0/prs(jc,k))**rd_o_cpd

          ! Find the level, where pressure approximately corresponds to the 
          ! average pressure of the well mixed layer. Simply assume a threshold
          ! of ml_depth/2 as average pressure in the layer, if this threshold 
          ! is surpassed the level with approximate mean pressure is found
          IF (prs(jc,k) > prs(jc,nlev) - ml_depth*0.5_wp) THEN
            kstart(jc) = k
          ENDIF

          k_ml(jc) = k - 1
#ifndef _OPENACC
        ELSE
          lexit(jc) = .TRUE.
#endif
        ENDIF

      ENDDO
    ENDDO
    ! Calculate the start values for the parcel ascent, 
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jc = i_startidx, i_endidx
      IF (k_ml(jc) < nlev) THEN
        qvp_start(jc) =  qvp_start(jc) / (nlev-k_ml(jc))
        tp_start (jc) =  tp_start (jc) / (nlev-k_ml(jc))
        te_start (jc) =  tp_start (jc)*(prs(jc,kstart(jc))/p0)**rd_o_cpd
      ELSE
        qvp_start(jc) =  qve(jc,nlev)
        te_start (jc) =  te (jc,nlev)
      END IF
    ENDDO
    !$ACC END PARALLEL

    ! The pseudoadiabatic ascent of the test parcel:
    IF (l3km) THEN
      CALL ascent ( i_startidx=i_startidx, i_endidx=i_endidx, kmoist=kmoist, te=te, qve=qve, prs=prs, hhl=hhl,  &
                    kstart=kstart, qvp_start=qvp_start, te_start=te_start,                                      &
                    acape=cape_ml, acape3km=cape_3km, acin=cin_ml, alcl=klcl, alfc=klfc,                        &
                    lacc= lacc )
    ELSE
      ! cape_3km is not computed
      CALL ascent ( i_startidx=i_startidx, i_endidx=i_endidx, kmoist=kmoist, te=te, qve=qve, prs=prs, hhl=hhl,  &
                    kstart=kstart, qvp_start=qvp_start, te_start=te_start,                                      &
                    acape=cape_ml, acin=cin_ml, alcl=klcl, alfc=klfc,                        &
                    lacc= lacc )
    ENDIF
    IF (llev) THEN
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(i_startidx, i_endidx, nlev) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        IF ((klcl(jc) .LE. 0) .OR. (klcl(jc) .GT. nlev)) THEN ! if index not defined
          lcl_ml(jc) = missing_value                       ! lcl=missing_val
        ELSE
          ! Lifting condensation level computations ABOVE GROUND level
          lcl_ml (jc) = 0.5_wp * ( hhl(jc,klcl(jc)) + hhl(jc,klcl(jc)+1) )-hhl(jc,nlev+1)
          ! Level of free condensation ABOVE GROUND level
        ENDIF
        IF ((klfc(jc) .LE. 0) .OR. (klfc(jc) .GT. nlev)) THEN ! if index not defined
          lfc_ml(jc) = missing_value                       ! lfc=missing_val
        ELSE
          lfc_ml (jc) = 0.5_wp * ( hhl(jc,klfc(jc)) + hhl(jc,klfc(jc)+1) )-hhl(jc,nlev+1)
          ! Equilibrium level ABOVE GROUND level
        ENDIF
      ENDDO
      !$ACC END PARALLEL
    ENDIF
    !$ACC WAIT(1)
    !$ACC END DATA
  END SUBROUTINE cal_cape_cin


  !------------------------------------------------------------------------------
  !
  !>
  !! Description:
  !! Compute the most unstable CAPE as an approximation: there is no search loop over the lowest 3000 m AGL,
  !!  but we rather search for the maximum equivalent potential temperature from the ground up to a
  !!  certain height z_limit AGL (usually 3000.0 m) and use this model
  !!  level as starting height for the test parcel.
  !!
  SUBROUTINE cal_cape_cin_mu(i_startidx, i_endidx, kmoist, z_limit, te, qve, prs, hhl,  &
                             cape_mu, cin_mu, lacc )

    IMPLICIT NONE

    INTEGER, INTENT (IN) ::  &
         i_startidx, i_endidx,  &  !> start and end indices of loops in horizontal patch
         kmoist                    !> start index for moist processes

    REAL    (wp),    INTENT (IN) ::  &
         z_limit,     & !> Max. height AGL for starting parcels to determine cape_mu/cin_mu
         te  (:,:),   & !> environment temperature,        dim: (nproma,nlev)
         qve (:,:),   & !> environment specific humidity,  dim: (nproma,nlev)
         prs (:,:),   & !> full level pressure,            dim: (nproma,nlev)
         hhl (:,:)      !> height of half levels,          dim: (nproma,nlev)

    REAL(wp),             INTENT(OUT) :: cape_mu(:), cin_mu(:)    !< output variables, dim: (nproma)

    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    INTEGER  :: jk, jc, nk
    REAL(wp) :: t_lcl    !> T at lifting condensation level approximation of Bolton (1980)
    REAL(wp) :: t_dew, zml, p_lcl
    INTEGER,  DIMENSION(SIZE(te,1))            :: kstart
    REAL(wp), DIMENSION(SIZE(te,1))            :: qvp_start, te_start
    REAL(wp), DIMENSION(SIZE(te,1),SIZE(te,2)) :: tequiv

    LOGICAL :: lzacc

    REAL (wp), PARAMETER :: p0 = 1.e5_wp   ! reference pressure for calculation of potential temperature
    REAL (wp), PARAMETER :: missing_value  = -999.9_wp   ! Missing value for CIN (if no LFC/CAPE was found),

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA &
    !$ACC   CREATE(tequiv, kstart, qvp_start, te_start) &
    !$ACC   PRESENT(te, qve, prs, hhl, cape_mu, cin_mu) &
    !$ACC   IF(lzacc)

    nk = SIZE(te,2)
    
    ! initialize outputs (good practice)
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(i_startidx, i_endidx, kmoist, nk, z_limit) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO jc = i_startidx, i_endidx
      cape_mu  (jc) = 0.0_wp 
      cin_mu   (jc) = missing_value
    ENDDO
    
    ! Compute equivalent potential temperature T_equiv approximation after Bolton (1980), Eq. 28:
    
    !$ACC LOOP SEQ
    DO jk = kmoist, nk
      !$ACC LOOP GANG VECTOR PRIVATE(zml, t_dew, t_lcl, p_lcl)
      DO jc = i_startidx, i_endidx
        zml = 0.5_wp * (hhl(jc,jk)+hhl(jc,jk+1)) - hhl(jc,nk+1)  ! m AGL
        IF (zml <= z_limit) THEN
          ! Dew point:
          t_dew = dewpoint_water (qve(jc,jk), prs(jc,jk))
          ! T at LCL approximation after Bolton (1980), Eq. 15. If te == t_dew, we are already at saturation:
          t_dew = MAX(t_dew,57.0_wp) ! for security
          t_lcl = 56.0_wp + 1.0_wp / ( 1.0_wp/(t_dew-56.0_wp) + LOG(MAX(te(jc,jk)/t_dew,1.0_wp))/800.0_wp )
          ! p at LCL:
          p_lcl = prs(jc,jk) * (t_lcl/te(jc,jk))**(1.0_wp/rd_o_cpd)
          ! Equivalent potential temperature using Bolton (1980), Eq. 28 at saturation:
          tequiv(jc,jk) = fthetae( t_lcl, p_lcl, qve(jc,jk) )
        ELSE
          tequiv(jc,jk) = -HUGE(1.0_wp)
        END IF
      END DO
    END DO

    ! The layer with the maximum T_equiv defines the starting values of the test parcel for the pseudo-adiabatic ascent
    !  to compute an approximation of cape_mu and cin_mu:
    !$ACC LOOP GANG VECTOR
    DO jc = i_startidx, i_endidx
      kstart(jc) = kmoist
    END DO
    !$ACC LOOP SEQ
    DO jk = kmoist+1, nk
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        IF (tequiv(jc,jk) > tequiv(jc,jk-1)) THEN
          kstart(jc) = jk
        END IF
      END DO
    END DO
    !$ACC LOOP GANG VECTOR
    DO jc = i_startidx, i_endidx
      qvp_start(jc) = qve(jc,kstart(jc))
      te_start (jc) = te (jc,kstart(jc))
    END DO
    !$ACC END PARALLEL

    ! The pseudoadiabatic ascent of the test parcel:
    CALL ascent ( i_startidx = i_startidx,  &
                  i_endidx   = i_endidx,    &
                  kmoist     = kmoist,      &
                  te         = te,          &
                  qve        = qve,         &
                  prs        = prs,         &
                  hhl        = hhl,         &
                  kstart     = kstart,      &
                  qvp_start  = qvp_start,   &
                  te_start   = te_start,    &
                  acape      = cape_mu,     &
                  acin       = cin_mu,      &
                  lacc       = lzacc )
    !$ACC END DATA

  END SUBROUTINE cal_cape_cin_mu

  SUBROUTINE cal_cape_cin_mu_COSMO(i_startidx, i_endidx, kmoist, te, qve, prs, hhl,  &
    cape_mu_COSMO, cin_mu_COSMO, lacc )

      ! Input data
      !----------- 
      INTEGER, INTENT (IN) ::  &
          i_startidx, i_endidx,  &  ! start and end indices of loops in horizontal patch
          kmoist                    ! start index for moist processes

      REAL    (wp),    INTENT (IN) ::  &
          te  (:,:),   & ! environment temperature
          qve (:,:),   & ! environment specific humidity
          prs (:,:),   & ! full level pressure
          hhl (:,:)      ! height of half levels

    ! Output data
      !------------ 
      REAL (wp), INTENT (OUT) :: &
        cape_mu_COSMO  (:),   & ! CAPE 
        cin_mu_COSMO   (:)      ! CIN  with respect to the starting values qvp_start, tp_start at level kstart

      LOGICAL, INTENT(IN), OPTIONAL :: lacc

    ! Local scalars and automatic arrays
      !-----------------------------------
        REAL (wp)             :: &
          acape       (SIZE(te,1)),   & ! CAPE Helper variable for the output from the routine ASCENT
          acin        (SIZE(te,1)),   & ! CIN Helper variable for the output from the routine ASCENT
          mup_lay_thck         ! thickness of layer in which most unstable parcel 
        
        INTEGER       :: &
        jc, k, nlev,           & !
        kstart(SIZE(te,1))       ! Model level corresponding to start height of parcel

        REAL (wp), PARAMETER :: missing_value  = -999.9_wp   ! Missing value for CIN (if no LFC/CAPE was found),

        LOGICAL,  DIMENSION(SIZE(te,1))            :: lcomp

        LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    nlev = SIZE(te,2)
    ! Thickness of the layer within which the most unstable parcel 
    ! is searched for -> smaller values of this parameter result in
    ! less computational cost, adapt if the code is running very slowly!
    mup_lay_thck = 30000._wp
    !$ACC DATA &
    !$ACC   PRESENT(te, qve, prs, hhl) &
    !$ACC   CREATE(cape_mu_COSMO, cin_mu_COSMO, acape, acin, lcomp, kstart) &
    !$ACC   IF(lzacc)

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(i_startidx, i_endidx, nlev) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO jc = i_startidx, i_endidx
      cape_mu_COSMO(jc)  = 0.0_wp
      cin_mu_COSMO (jc)  = missing_value
      acape        (jc)  = 0.0_wp
      acin         (jc)  = missing_value
      lcomp        (jc)  = .FALSE.
      kstart       (jc)  = nlev
    ENDDO
    !$ACC END PARALLEL
  

      

    !------------------------------------------------------------------------------
    ! Start computation. Overview:
    ! 
    ! parcelloop varies the start level
    ! of the parcel in most unstable calculation method.
    ! One single ascent of a parcel is calculated within the 
    ! SUBROUTINE ascent further below. 
    ! Here we make sure that "ascent" is being called with the 
    ! correct initial temperature, moisture, model level and 
    ! gridpoint indices. 
    !
    ! The initial model level, from where the parcel starts ascending is varied 
    ! until the pressure at this level is lower than the threshold mup_lay_thck, 
    ! defined above.
    !------------------------------------------------------------------------------
    parcelloop:  DO k = kmoist, nlev
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(k, nlev, i_startidx, i_endidx, mup_lay_thck) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        kstart(jc) = k
        lcomp (jc) = (prs(jc,k) > (prs(jc,nlev)-mup_lay_thck))
        acape (jc) = 0.0_wp
        acin  (jc) = missing_value
      ENDDO
      !$ACC END PARALLEL

      ! ! Take temperature and moisture of environment profile at current 
      ! ! level as initial values for the ascending parcel, call "ascent"
      ! ! to perform the dry/moist adiabatic parcel ascent and get back
      ! ! the calculated CAPE/CIN
      ! 
      ! 
      CALL ascent ( i_startidx=i_startidx, i_endidx=i_endidx, kmoist=kmoist,   & ! in (indeces)
                     te=te, qve=qve, prs=prs, hhl=hhl,                         & ! in (environment)
                     kstart=kstart, qvp_start=qve(:,k), te_start=te(:,k),      & ! in (initial conditions)
                     acape=acape, acin=acin, lcomp = lcomp,                    & ! out (cape, cin) in (logical array for computation)
                     lacc= lacc )
      
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(i_startidx, i_endidx) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        IF ( acape(jc) > cape_mu_COSMO(jc) ) THEN
          cape_mu_COSMO(jc)  = acape(jc)
          cin_mu_COSMO (jc)  = acin (jc)
        ENDIF
      ENDDO
      !$ACC END PARALLEL
    ENDDO parcelloop
    !$ACC WAIT(1)
    !$ACC END DATA
  END SUBROUTINE

  SUBROUTINE cal_si_sli_swiss ( i_startidx, i_endidx, kmoist,                    & ! in
                                  te, qve, prs, hhl, u, v,                         & ! in
                                  si, sli, swiss12, swiss00,                       & ! out 
                                  lacc )                                          ! in (optional)


    ! Input data
    !----------- 
    INTEGER, INTENT (IN) ::  &
         i_startidx, i_endidx,  &  ! start and end indices of loops in horizontal patch
         kmoist                    ! start index for moist processes

    REAL    (wp),    INTENT (IN) ::  &
         te  (:,:),   & ! environment temperature
         qve (:,:),   & ! environment specific humidity
         prs (:,:),   & ! full level pressure
         u   (:,:),   & ! environment zonal wind speed
         v   (:,:),   & ! environment meridional wind speed
         hhl (:,:)      ! height of half levels

    ! Output data
    !------------ 
    REAL (wp), INTENT (OUT) :: &
         si       (:),    & ! Showalter Index SI
         sli      (:),    & ! Surface Lifted Index SLI
         swiss12  (:),    & ! SWISS12 index
         swiss00  (:)       ! SWISS00 index


    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    ! Local scalars and automatic arrays
    !-----------------------------------

    REAL(wp) :: &
         qvp_start(SIZE(te,1)), & ! parcel initial specific humidity
         te_start (SIZE(te,1)), & ! parcel initial temperature
         hl_l, hl_u,            & ! height of full levels in order to find the closest model level
         e,                     & ! water vapor pressure
         td,                    & ! dew point temperature
         vh3000,                & ! norm of the horiz. wind vectors at approx. 3000m
         vh6000,                & ! norm of the horiz. wind vectors at approx. 3000m
         vhfirstlev               ! norm of the horiz. wind vectors at the first model level

    INTEGER :: &
         jc, k, nlev,             &
         kstart(SIZE(te,1)),     & ! Model level corresponding to start height of parcel
         k3000m(SIZE(te,1)),     & ! Model level corresponding to 3 km a.s.l.
         k6000m(SIZE(te,1)),     & ! Model level corresponding to 6 km a.s.l.
         k600  (SIZE(te,1)),     & ! Model level corresponding to 600 hPa
         k650  (SIZE(te,1))        ! Model level corresponding to 650 hPa

    REAL (wp), PARAMETER :: sistartprs = 85000.0_wp ! lower limit pressure for SI calculation, per definition 850hPa
    REAL (wp), PARAMETER :: missing_value  = -999.9_wp   ! Missing value for CIN (if no LFC/CAPE was found),

    LOGICAL :: lzacc

  ! initialization
  nlev = SIZE( te,2)

  CALL set_acc_host_or_device(lzacc, lacc)
  !$ACC DATA &
  !$ACC   PRESENT(te, qve, prs, u, v, hhl) &
  !$ACC   PRESENT(si, sli, swiss12, swiss00) &
  !$ACC   CREATE(kstart, k3000m, k6000m, k600, k650) &
  !$ACC   CREATE(qvp_start, te_start) &
  !$ACC   IF(lzacc)
  !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(i_startidx, i_endidx, nlev, kmoist) IF(lzacc)
  !$ACC LOOP GANG VECTOR
  DO jc = i_startidx, i_endidx
    kstart(jc) = nlev  
    si    (jc) = missing_value         
  ENDDO
  !------------------------------------------------------------------------------
  ! Section 1: Showalter Index SI calculation
  !
  ! Definition: Tp - Te at 500hPa, where Tp is temperature of parcel, ascending 
  ! from start level 850hPa and Te is environment temperature.
  ! Implementation here is done straightforward based on this definition.
  !------------------------------------------------------------------------------
    ! Loop through the levels, from highest pressure to lowest (ground -> up)
    ! in order to find the first level >= 850 hPa to initialize the parcel ascent
    !$ACC LOOP SEQ
    siloop: DO k = nlev, kmoist, -1
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        IF (prs(jc,k) >=  sistartprs) THEN
          kstart(jc) = k           
        ENDIF
      ENDDO
    ENDDO siloop
    ! set the initial conditions for the parcel ascent
    !$ACC LOOP GANG VECTOR
    DO jc = i_startidx, i_endidx
        te_start(jc) = te  (jc,kstart(jc))
        qvp_start(jc) = qve (jc,kstart(jc))
    ENDDO
    !$ACC END PARALLEL

    CALL ascent ( i_startidx=i_startidx, i_endidx=i_endidx,                & ! in (grid)
                  kmoist=kmoist, te=te, qve=qve, prs=prs, hhl=hhl,         & ! in (environment properties)
                  kstart=kstart, qvp_start=qvp_start, te_start=te_start,   & ! in (initial parcel conditions)
                  asi=si,                                                  & ! out (Showalter Index)
                  lacc= lacc )                                            ! in (GPU flag)
    
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(i_startidx, i_endidx, nlev) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO jc = i_startidx, i_endidx
        IF (prs(jc,nlev) < sistartprs) THEN
          si(jc) = missing_value
        ENDIF
    ENDDO
  !------------------------------------------------------------------------------
  ! Section 2: surface lifed index SLI calculation
  !
  ! Definition: Tp - Te at 500hPa, where Tp is temperature of parcel, ascending 
  ! from lowest level and Te is environment temperature.
  ! Implementation here is done straightforward based on this definition.
  !------------------------------------------------------------------------------
  
  ! initialization
  !$ACC LOOP GANG VECTOR
  DO jc = i_startidx, i_endidx
    kstart   (jc) = nlev  
    sli      (jc) = missing_value     
    te_start (jc) = te  (jc,kstart(jc))
    qvp_start(jc) = qve (jc,kstart(jc))    
  ENDDO
  !$ACC END PARALLEL

  CALL ascent ( i_startidx=i_startidx, i_endidx=i_endidx,                  & ! in (grid)
                  kmoist=kmoist, te=te, qve=qve, prs=prs, hhl=hhl,         & ! in (environment properties)
                  kstart=kstart, qvp_start=qvp_start, te_start=te_start,   & ! in (initial parcel conditions)
                  asi=sli,                                                 & ! out (Surface Lifted Index)
                  lacc= lacc )                                            ! in (GPU flag)

  !------------------------------------------------------------------------------
  ! Section 3: SWISS indices
  !
  ! Definition: SWISS00 and SWISS12 are two statistically based indices developped
  !             using data originating from the meteorological station of
  !             Payerne in Switzerland. SWISS00 has been developed with observations at
  !             00:00UTC whereas SWISS12 has been develped with observations at 12:00UTC.
  !------------------------------------------------------------------------------
  ! For these indeces, we need to obtain the levels corresponding to 600 hPa, 650 hPa,
  ! 3000m and 6000m (NOTE that this is not Height Above Ground (HAG) but Above mean Sea Level (a.s.l.))
  !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(i_startidx, i_endidx, nlev, kmoist) IF(lzacc)
  ! inizialisation
  !$ACC LOOP GANG VECTOR PRIVATE(hl_l, hl_u)
  DO jc = i_startidx, i_endidx
    k600  (jc) = -1
    k650  (jc) = -1
    k3000m(jc) = -1
    k6000m(jc) = -1  
  ENDDO       
  
  ! Find these indeces
  !$ACC LOOP SEQ
  kloop: DO k = nlev-1, kmoist+1, -1
    !$ACC LOOP GANG VECTOR PRIVATE(hl_l, hl_u)
    DO jc = i_startidx, i_endidx
      !Find the k index coressponding to the model level closest to 600 hPa
      IF ((prs (jc,k) >= 60000._wp).AND.(prs (jc,k-1) <= 60000._wp))THEN
          IF (abs(prs (jc,k)- 60000._wp) <= abs(prs(jc,k-1)- 60000._wp))THEN
            k600(jc) = k
          ELSE
            k600(jc) = k-1
          ENDIF
      ENDIF

      !Find the k index coressponding to the model level closest to 650 hPa
      IF ((prs (jc,k) >= 65000._wp).AND.(prs (jc,k-1) <= 65000._wp))THEN
          IF (abs(prs (jc,k)- 65000._wp) <= abs(prs(jc,k-1)- 65000._wp))THEN
            k650(jc) = k
          ELSE
            k650(jc) = k-1
          ENDIF
      ENDIF

      hl_l = 0.5_wp * (hhl(jc,k) + hhl(jc,k+1)) ! height a.s.l.
      hl_u = 0.5_wp * (hhl(jc,k-1) + hhl(jc,k)) ! height a.s.l.
      !Find the k index corresponding to the model level closest to 3000 m
      IF ((hl_l  <=  3000._wp) .AND. (hl_u  >=  3000._wp)) THEN
          IF (abs(hl_l - 3000._wp) <= abs(hl_u - 3000._wp)) THEN
            k3000m(jc) = k
          ELSE
            k3000m(jc) = k-1
          ENDIF
      ENDIF

      !Find the k index corresponding to the model level closest to 6000 m
      IF ((hl_l  <=  6000._wp) .AND. (hl_u  >=  6000._wp)) THEN
          IF (abs(hl_l - 6000._wp) <= abs(hl_u - 6000._wp)) THEN
            k6000m(jc) = k
          ELSE
            k6000m(jc) = k-1
          ENDIF
      ENDIF
    ENDDO
  ENDDO kloop
  !------------------------------------------------------------------------------
  ! Section 3a: SWISS00 index
  !
  !             Its components are:
  !               - the showalter index (see above)
  !               - the wind shear between 3000m and 6000m
  !               - the dew point depression at 600hPa
  !
  !             A SWISS00 value less than 5.1 means "likely thunderstorms".
  !------------------------------------------------------------------------------
  !$ACC LOOP GANG VECTOR PRIVATE(e, td, vh6000, vh3000)
  DO jc = i_startidx, i_endidx
    IF ( (k600(jc)  <  0) .OR. (k3000m(jc)  <  0  .OR.  k6000m(jc) < 0 .OR. si(jc) < -900.0_wp) ) THEN
      swiss00(jc) = missing_value
    ELSE
      !Compute vapor pressure at approximately 600 hPa using index k600
      e = prs(jc,k600(jc))*qve(jc,k600(jc))/(qve(jc,k600(jc))+rdv*(1 - qve(jc,k600(jc))))
      !Compute dew point temperature at approximately 650 hPa
      !Obtained using this identity: e = fesatw(td) at approximately 650 hPa
      !It means that the vapor pressure is equal to the saturation vapor
      !pressure at the dew point temperature

      IF (e <= 0.0_wp) THEN
          swiss00(jc) = missing_value
      ELSE
        td = ( 273.16_wp * 17.2693882_wp -log(e/610.78_wp)* 35.86_wp ) / ( 17.2693882_wp -log(e/ 610.78_wp) )
        !Compute norm of horizontal wind vectors at approximately 6000 m
        !using k6000 index
        vh6000 = sqrt ( u (jc,k6000m(jc)) **2 + v (jc,k6000m(jc)) **2)
        !Compute norm of horizontal wind vectors at approximately 3000 m
        !using k3000 index
        vh3000 = sqrt ( u (jc,k3000m(jc)) **2 + v (jc,k3000m(jc)) **2)
        swiss00(jc) = si(jc) + 0.4_wp*(vh6000-vh3000) + 0.1_wp*MAX(0.0_wp,(te(jc,k600(jc))-td))
      ENDIF
    ENDIF
  ENDDO
  !------------------------------------------------------------------------------
  ! Section 3b: SWISS12 index
  !
  !             Its components are:
  !               - the surface lifted index (see above)
  !               - the wind shear between the ground and 3000m
  !               - the dew point depression at 650hPa
  !
  !             The concept of "wind at the ground" is ambiguous: the wind at the
  !             surface is equal to zero but we could also consider the wind at
  !             10m for example as "ground value". We use the wind at the first
  !             model level as an approximation for the surface wind.
  !
  !             A SWISS12 value less than 0.6 means "likely thunderstorms".
  !------------------------------------------------------------------------------
  !$ACC LOOP GANG VECTOR PRIVATE(e, td, vh3000, vhfirstlev)
  DO jc = i_startidx, i_endidx
    IF ( (k650(jc)  <  0) .OR. (k3000m(jc)  <  0) ) THEN
      swiss12(jc) = missing_value
    ELSE
      !Compute vapor pressure at approximately 650 hPa using index k650
      e = prs(jc,k650(jc))*qve(jc,k650(jc))/(qve(jc,k650(jc))+rdv*(1 - qve(jc,k650(jc))))

      !Compute dew point temperature at approximately 650 hPa
      !Obtained using this identity: e = fesatw(td) at approximately 650 hPa
      !It means that the vapor pressure is equal to the saturation vapor
      !pressure at the dew point temperature

      IF (e <= 0.0_wp) THEN
          swiss12(jc) = missing_value
      ELSE
        ! dewpoint: e = vapor pressure
        ! td = dewpoint_water(qve(i, k650(i)), prs(i, k650(i)))
        ! td = ftd (e)
        ! td = ( b3*b2w -log(e/b1)*b4w ) / ( b2w -log(e/b1) )
        td = ( 273.16_wp * 17.2693882_wp -log(e/610.78_wp)* 35.86_wp ) / ( 17.2693882_wp -log(e/ 610.78_wp) )
        !Compute norm of horizontal wind vectors  at approximately 3000 m
        vh3000 = sqrt ( u(jc,k3000m(jc)) **2 + v(jc,k3000m(jc)) **2)
        !Compute norm of horizontal wind vectors  at the first model level
        ! WARNING: Currently the first model level is at 10m a.s.l. but this
        !          height could change!!!
        vhfirstlev = sqrt ( u(jc,nlev) **2 + v(jc,nlev) **2)
        !Compute SWISS12 Index with the above calculated parameters and the
        !surface lifted index
        swiss12(jc) = sli(jc) - 0.3_wp*(vh3000-vhfirstlev) + 0.3_wp*(MAX(0.0_wp,(te(jc,k650(jc))-td)))
      ENDIF
    ENDIF
  ENDDO
  !$ACC END PARALLEL
  !$ACC WAIT(1)
  !$ACC END DATA
  END SUBROUTINE

  SUBROUTINE cal_cloudtop(i_startidx, i_endidx, kmoist,          & ! in
                          clc, h,                                & ! in
                          cloudtop,                              & ! out
                          lacc)                                 ! in (optional)
  
    ! Input data
    !----------- 
    INTEGER, INTENT (IN) ::  &
         i_startidx, i_endidx,  &  ! start and end indices of loops in horizontal patch
         kmoist                    ! start index for moist processes

    REAL    (wp),    INTENT (IN) ::  &
         clc          (:,:),             & ! cloud coverage
         h            (:,:)                ! Geometric height at full level center

    ! Output data
    !------------ 
    REAL (wp), INTENT (OUT) :: &
         cloudtop    (:)           ! CLOUDTOP: 

    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    ! Local scalars and automatic arrays
    !-----------------------------------
    INTEGER :: &
         jc, k, nlev             

#ifndef _OPENACC
    LOGICAL :: lexit(SIZE(clc,1))
#endif

    LOGICAL :: lzacc

    ! Local parameters:
    !------------------
    LOGICAL :: &
    lcloudtop(SIZE(clc,1))         ! logical array to stop computation after findind cloudtop

    REAL (wp), PARAMETER :: missing_value  = -999.9_wp   ! Missing value


    IF (PRESENT(lacc)) THEN
      lzacc = lacc
    ELSE
      lzacc = .FALSE.
    END IF
    ! Initialize
    nlev = SIZE( clc,2)
    !$ACC DATA &
    !$ACC   PRESENT(clc, h, cloudtop) &
    !$ACC   CREATE(lcloudtop) &
    !$ACC   IF(lzacc)

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(i_startidx, i_endidx, kmoist, nlev) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO jc = i_startidx, i_endidx
      cloudtop     (jc) = missing_value
      lcloudtop    (jc) = .FALSE.
    ENDDO
    !$ACC LOOP SEQ
    DO k = kmoist, nlev ! the only difference with CEILING computation, looping from top to surface (nlev is surface)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        IF((clc(jc,k) > 0.5) .AND. (.NOT. lcloudtop(jc))) THEN
          cloudtop  (jc) = h(jc,k)
          lcloudtop (jc) = .TRUE.
        ENDIF
      ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$ACC WAIT(1)
    !$ACC END DATA
  END SUBROUTINE

  SUBROUTINE ascent ( i_startidx, i_endidx, kmoist, te, qve, prs, hhl,  & ! in
                      kstart, qvp_start, te_start,                      & ! in
                      acape, acape3km, acin,                            & ! out (optional)
                      alcl, alfc,                                       & ! out (optional)
                      asi,                                              & ! out (optional)
                      lcomp, lacc )                                    ! in (optional)

    !------------------------------------------------------------------------------
    !
    !>
    !! Description:
    !!  Computation of Convective Available Potential Energy CAPE,
    !!  Convective Inhibition CIN based on parcel theory with respect to a test
    !!  parcel having T = tp_start, QV = qvp_start and a height corresponding to
    !!  model level kstart. Based on this, different versions of cape and cin
    !!  may be computed, e.g., mixed layer cape/cin, most unstable cape/cin or the like.
    !!
    !!  This subroutine is based on COSMO code.
    !!        Helmut Frank
    !! 
    !! Input:  
    !!         - Temperature, specific humidity and pressure of environment
    !!
    !! Output: 
    !!         - cape/cin: CAPE/CIN based on a parcel with thermodynamical 
    !!                           properties of the lowest mean layer in the PBL (50hPa)
    !!      
    !! Motivation: 
    !!  Current parameter CAPE_CON is calculated in LM in the framework of the 
    !!  convective parametrisation scheme. Therefore this parameter is only available
    !!  at those gridpoints, where the scheme is called, but not continuously on the 
    !!  whole domain. This subroutine, on the other hand, provides continuous fields. 
    !!
    !! Method:
    !!  A dry/moist parcel ascent is performed following classic parcel theory.
    !!  Moist adiabatic ascent is calculated iteratively with an appropriate scheme.
    !!  Based on the temperature and moisture of the ascending parcel, CAPE and CIN
    !!  are computed, closely following the recommendations of Doswell and Rasmussen 
    !!  (1994), including a virtual temperature correction and searching for the 
    !!  most unstable parcel in the lower troposphere. Additionally, a mixed layer 
    !!  CAPE as well as the traditional Showalter Index and the surface lifted 
    !!  index are computed as further variables. 
    !!
    !!  References used during development: 
    !!  - C. A. Doswell and Rasmussen, E. N.: The Effect of Neglecting the 
    !!    Virtual Temperature Correction on CAPE Calculations. 
    !!    Weather and Forecasting, 9, 625-629.
    !!
    !!  - K. A. Emanuel (1994): Atmospheric Convection. Oxford University Press.
    !!
    !!  - H. Huntrieser et al. (1997): Comparison of Traditional and Newly Developed 
    !!    Thunderstorm Indices for Switzerland. Weather and Forecasting, 12, 
    !!    108-125.
    !!
    !!  - D. Bolton (1980): The Computation of Equivalent Potential Temperature. 
    !!    Monthly Weather Review, 108, 1046-1053
    !!
    !!  - Davies, J.M.,2002: On low-level thermodynamic parameters
    !!    associated with tornadic and nontornadic supercells.
    !!    Preprints, 21st Conf. On Severe Local Storms, San Antonio, Amer. Meteor. Soc.
    !!    http://members.cox.net/jondavies1/LLthermo.PDF
    !!

    ! Input data
    !----------- 
    INTEGER, INTENT (IN) ::  &
         i_startidx, i_endidx,  &  ! start and end indices of loops in horizontal patch
         kmoist                    ! start index for moist processes

    REAL    (wp),    INTENT (IN) ::  &
         te  (:,:),   & ! environment temperature
         qve (:,:),   & ! environment specific humidity
         prs (:,:),   & ! full level pressure
         hhl (:,:)      ! height of half levels

    REAL(wp), INTENT(IN) :: &
         qvp_start(:), & ! parcel initial specific humidity
         te_start (:)    ! parcel initial temperature

    INTEGER, INTENT(IN) :: &
         kstart(:)  ! Model level corresponding to start height of parcel

    ! Output data
    !------------ 
    REAL (wp), INTENT (OUT), OPTIONAL :: &
         acape   (:),   & ! CAPE 
         acape3km(:),   & ! CAPE_3KM
         acin    (:),   & ! CIN
         asi     (:)      ! SI

    INTEGER, INTENT (OUT), OPTIONAL :: & ! these are also in the local scalars so that they can be optional outputs
         alcl(:),     & ! Indices for Lifting Condensation Level LCL,
         alfc(:)        ! Level of Free Convection LFC

    LOGICAL, INTENT(IN), OPTIONAL :: &
      lacc   ,      & ! If true, use openacc
      lcomp  (:)      ! Most Unstable Computation

    ! Local scalars and automatic arrays
    !-----------------------------------
    INTEGER :: &
      nlev, jc, k,             & ! Indices of input/output fields
      k3000m(SIZE(te,1)),     & ! Indices for endpoint of 3KM ascent
      lcllev(SIZE(te,1)),     & ! Indices for Lifting Condensation Level LCL,
      lfclev(SIZE(te,1)),     & ! Level of Free Convection LFC
      ellev (SIZE(te,1))        ! Equilibrium Level EL

    REAL (wp), PARAMETER :: p0 = 1.e5_wp   ! reference pressure for calculation of potential temperature
    REAL (wp), PARAMETER :: missing_value  = -999.9_wp   ! Missing value for CIN (if no LFC/CAPE was found),
   

    ! The following parameters are help values for the iterative calculation 
    ! of the parcel temperature during the moist adiabatic ascent
    REAL    (wp)             :: esat,tguess1,tguess2,thetae1,thetae2
#ifdef __SX__
    REAL    (wp)             :: tguess1v(SIZE(te,1))
    LOGICAL                  :: lcalc(SIZE(te,1))
#endif
    ! REAL    (wp)             :: rp, r1,r2
    REAL    (wp)             :: q1, q2
    REAL    (wp), PARAMETER  :: eps=0.03
    REAL    (wp), PARAMETER  :: sistopprs = 50000.0_wp ! upper limit pressure for SI calculation, per definition 500hPa


    ! this parameter helps to find the LFC above a capping inversion in cases, 
    ! where a LFC already was found in an unstable layer in the convective 
    ! boundary layer below. 
    REAL (wp), PARAMETER :: cc_comp    = 2.0_wp                      

    INTEGER ::    icount              ! counter for the iterative process

    REAL (wp) ::             &
         cape    (SIZE(te,1)),  & ! CAPE computed in this parcel ascent
         cape3km (SIZE(te,1)),  & ! CAPE_3KM computed in this parcel ascent
         cin     (SIZE(te,1)),  & ! CIN computed in this parcel ascent
         si      (SIZE(te,1)),  & ! SI computed in this parcel ascent
         cin_help(SIZE(te,1)),  & ! help variable, the CIN above the LFC
         buo     (SIZE(te,1)),  & ! parcel buoyancy at level k
         tp      (SIZE(te,1)),  & ! temperature profile of ascending parcel
         tp_start(SIZE(te,1)),  & ! parcel initial potential temperature
         qvp     (SIZE(te,1)),  & ! specific moisture profile of ascending parcel
         thp     (SIZE(te,1)),  & ! 1st guess theta_e of parcel for iterative 
         tvp,                   & ! virtual temperature of parcel at level k
         tve,                   & ! virtual temperature of environment at level k
         buo_belo,              & ! parcel buoyancy of level k+1 below
         esatp,                 & ! saturation vapour pressure at level k
         qvsp,                  & ! saturation specific humidity at level k
         hl_l,hl_u                ! height of full levels in order to find the closest model level
    ! calculation of moist adiabatic ascent

    INTEGER :: lfcfound(SIZE(te,1))   ! flag indicating if a LFC has already been found
    ! below, in cases where several EL and LFC's occur
    LOGICAL :: &
    lzacc,     & ! GPU flag
    lacape,    & ! check if ACAPE is present
    lacape3km, & ! check if ACAPE3KM is present
    lacin,     & ! check if ACIN is present
    lalcl,     & ! check if ALCL is present
    lalfc,     & ! check if ALFC is present
    lasi,      & ! check if ASI is present
    lmu(SIZE(te,1)) ! local non-optional version of lcomp for
                    ! COSMO-like Most Unstable computation

    CALL set_acc_host_or_device(lzacc, lacc)



    IF (PRESENT(acape)) THEN
      lacape = .TRUE.
    ELSE
      lacape = .FALSE.
    END IF

    IF (PRESENT(acape3km)) THEN
      lacape3km = .TRUE.
    ELSE
      lacape3km = .FALSE.
    END IF

    IF (PRESENT(acin)) THEN
      lacin = .TRUE.
    ELSE
      lacin = .FALSE.
    END IF

    IF (PRESENT(alcl)) THEN
      lalcl = .TRUE.
    ELSE
      lalcl = .FALSE.
    END IF

    IF (PRESENT(alfc)) THEN
      lalfc = .TRUE.
    ELSE
      lalfc = .FALSE.
    END IF

    IF (PRESENT(asi)) THEN
      lasi = .TRUE.
    ELSE
      lasi = .FALSE.
    END IF
    !$ACC DATA &
    !$ACC   PRESENT(te, qve, prs, hhl, qvp_start, te_start, kstart) &
    !$ACC   PRESENT(acape, acape3km, acin) &
    !$ACC   PRESENT(alfc, alcl, asi, lcomp) &
    !$ACC   CREATE(k3000m, lcllev, lfclev, ellev) &
    !$ACC   CREATE(cape, cape3km, cin, si, cin_help, buo, tp, tp_start) &
    !$ACC   CREATE(qvp, thp, lfcfound, lmu) &
    !$ACC   IF(lzacc)
    nlev = SIZE( te,2)

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(i_startidx, i_endidx) IF(lzacc)
    IF (PRESENT(lcomp)) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = i_startidx, i_endidx  
        lmu(jc) = lcomp(jc)
      ENDDO
    ELSE
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = i_startidx, i_endidx  
        lmu(jc) = .TRUE.
      ENDDO
    END IF

    !------------------------------------------------------------------------------
    !
    ! Description:
    !   A single parcel ascent is performed, based on the given start 
    !   values kstart (level), te_start (initial parcel temperature) and
    !   qvp_start (initial parcel specific humidity). 
    !
    !------------------------------------------------------------------------------

    ! Initialization
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jc = i_startidx, i_endidx  
      tp_start(jc) = te_start(jc) * (p0/prs(jc,kstart(jc)))**rd_o_cpd
      cape    (jc) = 0.0_wp
      cape3km (jc) = 0.0_wp
      cin     (jc) = 0.0_wp ! not missing value because we add to this 0_wp
      si      (jc) = missing_value
      lcllev  (jc) = 0
      lfclev  (jc) = 0
      ellev   (jc) = 0
      lfcfound(jc) = 0
      cin_help(jc) = 0.0_wp
      tp      (jc) = 0.0_wp
      qvp     (jc) = 0.0_wp               
      buo     (jc) = 0.0_wp
      k3000m  (jc) = -1 ! this is done also in the COSMO subroutine so it is reproduced here
    ENDDO
    !$ACC END PARALLEL

    IF(lacape3km) THEN
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(nlev, kmoist, i_startidx, i_endidx) IF(lzacc)
      !$ACC LOOP SEQ
      DO k = nlev-1, kmoist+1, -1
        !$ACC LOOP GANG VECTOR PRIVATE(hl_l, hl_u)
        DO jc = i_startidx, i_endidx
          hl_l = 0.5_wp * (hhl(jc,k) + hhl(jc,k+1)) - hhl(jc,nlev+1)
          hl_u = 0.5_wp * (hhl(jc,k-1) + hhl(jc,k)) - hhl(jc,nlev+1)
          IF ((hl_l  <=  3000._wp) .AND. (hl_u  >=  3000._wp)) THEN
            IF (abs(hl_l - 3000._wp) <= abs(hl_u - 3000._wp)) THEN
              k3000m(jc) = k
            ELSE
              k3000m(jc) = k-1
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ENDIF

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    ! Loop over all model levels above kstart
    !$ACC LOOP SEQ
    kloop: DO k = nlev, kmoist, -1
      !$ACC LOOP GANG VECTOR PRIVATE(esatp, qvsp)
      DO jc = i_startidx, i_endidx
        IF ( k > kstart(jc) .OR. (.NOT. lmu(jc)) ) CYCLE 
        ! Dry ascent if below cloud base, assume first level is not saturated 
        ! (first approximation)
        IF (k > lcllev(jc)) THEN
          tp (jc)   = tp_start(jc)*( prs(jc,k)/p0)**rd_o_cpd   ! dry adiabatic process
          qvp(jc)   = qvp_start(jc)                           ! spec humidity conserved

          ! Calculate parcel saturation vapour pressure and saturation 
          ! specific humidity
          esatp = esat_water( tp(jc))
          qvsp  = fqvs( esatp, prs(jc,k))

          ! Check whether parcel is saturated or not and 
          ! no LCL was already found below
          IF ( (qvp(jc) >= qvsp) .AND. (lcllev(jc) == 0) ) THEN  
            lcllev(jc) = k                                    ! LCL is reached

            ! Moist ascent above LCL, first calculate an approximate thetae to hold 
            ! constant during the remaining ascent
            !         rp      = qvp(jc)/( 1._wp - qvp(jc) )
            !         thp(jc)  = fthetae( tp(jc),prs(jc,k),rp )
            thp(jc)  = fthetae( tp(jc),prs(jc,k), qvp(jc) )
          ENDIF
        ENDIF

#ifdef __SX__
      ENDDO ! i = i_startidx, i_endidx
      !$ACC LOOP GANG VECTOR
      ! Vectorized version
      DO jc = i_startidx, i_endidx
        lcalc(jc) = .FALSE.
        IF ( k > kstart(jc) .OR. (.NOT. lmu(jc)) ) CYCLE 
        IF ( k <= lcllev(jc) ) THEN ! If we are above the LCL
          ! The scheme uses a first guess temperature, which is the parcel
          ! temperature at the level below. If it happens that the initial
          ! parcel is already saturated, the environmental temperature
          ! is taken as first guess instead
          IF (  k == kstart(jc) ) THEN
            tguess1v(jc) = te(jc,kstart(jc))
          ELSE
            tguess1v(jc) = tp(jc)
          END IF
          lcalc(jc) = .TRUE.
        ENDIF
      ENDDO ! jc = i_startidx, i_endidx

      !$ACC LOOP SEQ
      ! Calculate iteratively parcel temperature from thp, prs and 1st guess tguess1
      DO icount = 1, 21
        IF (COUNT(lcalc(i_startidx:i_endidx)) > 0) THEN
          !$ACC LOOP GANG VECTOR PRIVATE(esat, q1, thetae1, tguess2, esat, q2, thetae2)
          DO jc = i_startidx, i_endidx
            IF ( lcalc(jc) ) THEN
              esat     = esat_water( tguess1v(jc))
              q1       = fqvs( esat, prs(jc,k) )
              thetae1  = fthetae( tguess1v(jc),prs(jc,k),q1)

              tguess2  = tguess1v(jc) - 1.0_wp
              esat     = esat_water( tguess2)
              q2       = fqvs( esat, prs(jc,k) )
              thetae2  = fthetae( tguess2,prs(jc,k),q2)

              tguess1v(jc)  = tguess1v(jc)+(thetae1-thp(jc))/(thetae2-thetae1)

              IF ( ABS( thetae1-thp(jc)) < eps .OR. icount > 20) THEN
                tp(jc) = tguess1v(jc)
                lcalc(jc) = .false.
              END IF
            END IF
          ENDDO ! jc = i_startidx, i_endidx
        ELSE
          EXIT
        ENDIF
      END DO

      ! update specific humidity of the saturated parcel for new temperature
      !$ACC LOOP GANG VECTOR PRIVATE(esatp)
      DO jc = i_startidx, i_endidx
        IF ( k > kstart(jc) .OR. (.NOT. lmu(jc)) ) CYCLE 
        IF ( k <= lcllev(jc) ) THEN
          esatp  = esat_water( tp(jc))
          qvp(jc) = fqvs( esatp,prs(jc,k))
        END IF
      ENDDO ! jc = i_startidx, i_endidx
      !$ACC LOOP GANG VECTOR PRIVATE(tguess1, icount, esat, q1, thetae1, tguess2, q2, thetae2, esatp, tvp, tve, buo_belo)
      DO jc = i_startidx, i_endidx
        IF ( k > kstart(jc) .OR. (.NOT. lmu(jc)) ) CYCLE 
#else

        ! Moist adiabatic process: the parcel temperature during this part of 
        ! the ascent is calculated iteratively using the iterative newton
        ! scheme, assuming the equivalent potential temperature of the parcel 
        ! at the LCL (thp) is held constant. The scheme converges usually within
        ! few (less than 10) iterations, its accuracy can be tuned with the 
        ! parameter "eps", a value of 0.03 is tested and recommended. 

        IF ( k <= lcllev(jc) ) THEN                                
          ! The scheme uses a first guess temperature, which is the parcel 
          ! temperature at the level below. If it happens that the initial 
          ! parcel is already saturated, the environmental temperature 
          ! is taken as first guess instead
          IF (  k == kstart(jc) ) THEN
            tguess1 = te(jc,kstart(jc))            
          ELSE
            tguess1 = tp(jc)
          END IF
          icount = 0       ! iterations counter

          ! Calculate iteratively parcel temperature from 
          ! thp, prs and 1st guess tguess1
          DO
            esat     = esat_water( tguess1)
            !         r1       = rdv*esat/(prs(jc,k)-esat)
            !         thetae1  = fthetae( tguess1,prs(jc,k),r1)
            q1       = fqvs( esat, prs(jc,k))
            thetae1  = fthetae( tguess1,prs(jc,k),q1)

            tguess2  = tguess1 - 1.0_wp
            esat     = esat_water( tguess2)
            !         r2       = rdv*esat/(prs(jc,k)-esat)
            !         thetae2  = fthetae( tguess2,prs(jc,k),r2)
            q2       = fqvs( esat, prs(jc,k))
            thetae2  = fthetae( tguess2,prs(jc,k),q2)

            tguess1  = tguess1+(thetae1-thp(jc))/(thetae2-thetae1)
            icount   = icount    + 1   

            IF ( ABS( thetae1-thp(jc)) < eps .OR. icount > 20 ) THEN
              tp(jc) = tguess1
              EXIT
            END IF
          END DO

          ! update specific humidity of the saturated parcel for new temperature
          esatp  = esat_water( tp(jc))
          qvp(jc) = fqvs( esatp,prs(jc,k))
        END IF
#endif       

        ! Calculate virtual temperatures of parcel and environment
        tvp    = tp(jc  ) * (1.0_wp + vtmpc1*qvp(jc  )/(1.0_wp - qvp(jc  )) )  
        tve    = te(jc,k) * (1.0_wp + vtmpc1*qve(jc,k)/(1.0_wp - qve(jc,k)) ) 

        ! Calculate the buoyancy of the parcel at current level k, 
        ! save buoyancy from level k+1 below (buo_belo) to check if LFC or EL have been passed
        buo_belo = buo(jc)
        buo(jc)   = tvp - tve

        ! Check for level of free convection (LFC) and set flag accordingly. 
        ! Basic LFC condition is that parcel buoyancy changes from negative to 
        ! positive (comparison of buo with buo_belo). Tests showed that very 
        ! often the LFC is already found within the boundary layer below even if 
        ! significant capping inversions are present above (and since CIN is only
        ! defined below the LFC no CIN was accumulated in these cases.)
        ! To handle these situations in a meteorologically meaningful way an 
        ! additional flag "lfcfound" was introduced which is initially zero but 
        ! set to 1 if a second LFC was found, under the condition that the CIN 
        ! within the capping inversion is greater than the CAPE in the convective
        ! boundary layer below times the factor cc_comp (cc_comp = 1 - 2 
        ! recommended.)
        ! Help variable CIN_HELP saves all contributions to the total cin above 
        ! the LFC and has to be subtracted at the end from the final CIN in order
        ! to get the CIN only below the LFC (this is necessary since we do not 
        ! know yet where exactly we will find an LFC when performing the ascent
        ! from bottom to top in a stepwise manner.)

        ! Find the first LFC
        IF ( (buo(jc) > 0.0_wp) .AND. (buo_belo <= 0.0_wp)            &
             .AND. ( lfcfound(jc)==0) ) THEN

          ! Check whether it is an LFC at one of the lowest model levels 
          ! (indicated by CAPE=0)
          IF ( (cape(jc) > 0.0_wp) .AND. ( lfcfound(jc) == 0 ) ) THEN
            ! Check if there is a major capping inversion below, defined as 
            ! having CIN with an absolute value larger than the CAPE accumulated
            ! below times some arbitrary factor cc_comp - if this is the case the
            ! LFC index "lfclev" is updated to the current level k and 
            ! "lfcfound"-flag is now set to 1 assuming that we have found the 
            ! level of free convection finally. 
            IF ( cc_comp * ABS(cin_help(jc)) > cape(jc) ) THEN
              lfclev   (jc) = k
              cape     (jc) = 0.0_wp
              cape3km  (jc) = 0.0_wp
              cin_help (jc) = 0.0_wp
              lfcfound (jc) = 1
            ENDIF
          ELSE
            ! the LFC found is near the surface, set the LFC index to the current
            ! level k (lfclev) but do not set the flag "lfcfound" to zero to 
            ! indicate that a further LFC may be present above the boundary layer
            ! and an eventual capping inversion. Reset the CIN_HELP to zero to 
            ! store the contribution of CIN above this LFC.
            lfclev(jc)   = k
            cin_help(jc) = 0.0_wp
          ENDIF
        ENDIF
        
        IF ( (buo(jc) < 0_wp) .AND. (buo_belo >= 0_wp) .AND. (lfclev(jc) /= 0) ) THEN
          ellev(jc) = k
        ENDIF
        ! Accumulation of CAPE and CIN according to definition given in Doswell 
        ! and Rasmussen (1994), 
        IF ( (buo(jc) >= 0.0_wp) .AND. (k <= lfclev(jc)) ) THEN   
          cape(jc)  = cape(jc)  + (buo(jc)/tve)*grav*(hhl(jc,k) - hhl(jc,k+1))
          IF (lacape3km) THEN
             IF ( k3000m(jc) > 0 .AND. k >= k3000m(jc) ) THEN
                cape3km(jc) = cape3km(jc) + (buo(jc)/tve)*grav*(hhl(jc,k) - hhl(jc,k+1))
             ENDIF
          ENDIF
        ELSEIF ( (buo(jc) < 0.0) .AND. (k < kstart(jc)) ) THEN  
          ! buo is negative, hhl(jc,k) > hhl(jc,k+1) so we add a negative contribution
          cin(jc)      = cin(jc)      + (buo(jc)/tve)*grav*(hhl(jc,k) - hhl(jc,k+1))
          cin_help(jc) = cin_help(jc) + (buo(jc)/tve)*grav*(hhl(jc,k) - hhl(jc,k+1))
        ENDIF
        
        ! If 500hPa level approximately reached, save parcel temperature for 
        ! calculation of Showalter Index. Do not check at lowest level of parcel
        ! ascent since pressure at level below is not defined if kstart=kdim 
        ! (parcel starting at lowest model level.)
        IF (lasi .AND. (k < nlev)) THEN
          ! Assume 500hPa level is approximately reached if model level pressure
          ! changes from >500hPa to <=500hPa
          IF ( (prs(jc,k) <= sistopprs) .AND. (prs(jc,k+1) > sistopprs) ) THEN
            asi(jc) = te(jc,k)-tp(jc)
          ENDIF
        ENDIF
      ENDDO ! jc = i_startidx, i_endidx
    ENDDO  kloop       ! End k-loop over levels
    !$ACC END PARALLEL

    ! Subtract the CIN above the LFC from the total accumulated CIN to 
    ! get only contriubtions from below the LFC as the definition demands.
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(i_startidx, i_endidx) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO jc = i_startidx, i_endidx
      ! make CIN positive
      cin(jc) = ABS (cin(jc) - cin_help(jc))
      ! set the CIN to missing value if no LFC was found or no CAPE exists
      IF ( (lfclev(jc) == 0) .OR. (ABS(cape(jc)) < 1.0E-8_wp)) cin(jc) = missing_value 
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(i_startidx, i_endidx, lacape, lacape3km, lacin, lalcl, lalfc) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO jc = i_startidx, i_endidx
      IF (lacape)    acape(jc)    = cape      (jc)
      IF (lacape3km) acape3km(jc) = cape3km   (jc)
      IF (lacin)     acin (jc)    = cin       (jc)
      IF (lalcl)     alcl (jc)    = lcllev    (jc)
      IF (lalfc)     alfc (jc)    = lfclev    (jc) 
    ENDDO
    !$ACC END PARALLEL
    !$ACC WAIT(1)
    !$ACC END DATA
  END SUBROUTINE ascent

  !!>
  !! Specific humidity at saturation as function of water vapor pressure zex,
  !! air pressure zpx, and specific humidity zqx.
  !! Initial version: Helmut Frank
  !! Corrected: Ulrich Blahak, 13.4.2022
  ELEMENTAL FUNCTION fqvs( zex, zpx)
!!!  ELEMENTAL FUNCTION fqvs( zex, zpx, zqx)

    REAL(wp), INTENT(IN) :: zex   ! vapor pressure        [Pa]
    REAL(wp), INTENT(IN) :: zpx   ! atmospheric pressure  [Pa]
!!!      REAL(wp), INTENT(IN) :: zqx   ! specific humidity     [kg/kg]
    REAL(wp)             :: fqvs  ! Equivalent potential temperature
    REAL(wp)             :: zex_lim

    !$ACC ROUTINE SEQ

    ! limit zex by its maximum permissible value, which is the total pressure:
    zex_lim = MIN(zex, zpx)
    
    fqvs = rdv * zex_lim/ (zpx - o_m_rdv*zex_lim )        
!!!    fqvs = zex/zpx *( rdv + o_m_rdv*zqx )        

  END FUNCTION fqvs

  !!>
  !! Equivalent potential temperature to hold constant during ascent, assuming saturation.
  !!   Bolton (1980), the simple approx. Eq. 28, not the full-fledged Eq. 43
  !! Initial version: Helmut Frank
  ELEMENTAL FUNCTION fthetae( ztx,zpx,zqx)

    REAL(wp), INTENT(IN) :: ztx     ! air temperature       [K]
    REAL(wp), INTENT(IN) :: zpx     ! atmospheric pressure  [Pa]
    !   REAL(wp), INTENT(IN) :: zrx     ! mixing ratio          [kg/kg]
    REAL(wp), INTENT(IN) :: zqx     ! specific humidity     [kg/kg]
    REAL(wp)             :: fthetae ! Equivalent potential temperature [K]

    REAL (wp), PARAMETER :: p0 = 1.e5_wp   ! reference pressure for calculation of potential temperature

    !$ACC ROUTINE SEQ

    !   fthetae = (p0/zpx)**rd_o_cpd *ztx*exp( alvdcp*zrx/ztx)  
    fthetae = (p0/zpx)**rd_o_cpd *ztx*EXP( alvdcp*zqx/(ztx*(1._wp-zqx)) )

  END FUNCTION fthetae

  !!>
  !! Saturation vapor pressure over flat water surface as function of temperature with a
  !!  safety measure for very low temperatures:
  !! Initial version: Ulrich Blahak, 13.4.2022
  ELEMENTAL FUNCTION esat_water(temp)
    IMPLICIT NONE

    REAL (wp)              :: esat_water
    REAL (wp), INTENT(IN)  :: temp

    !$ACC ROUTINE SEQ

    IF (temp > c4les + 1e-6_wp) THEN
      esat_water = c1es*EXP( c3les*(temp-tmelt)/(temp-c4les) )
    ELSE
      esat_water = 0.0_wp
    END IF
    
  END FUNCTION esat_water
  
  !!>
  !! Dew point with respect to water as function of specific hum. r and total pressure p
  !! Initial version: Ulrich Blahak, 13.4.2022
  ELEMENTAL FUNCTION dewpoint_water (r, p) RESULT(td)

    REAL(wp), INTENT(in) :: r  !< specific humidity
    REAL(wp), INTENT(in) :: p  !< total pressure

    REAL(wp)             :: td !< RESULT: dewpoint

    REAL(wp)             :: esatt, z

    !$ACC ROUTINE SEQ

    IF (r > 1e-16_wp) THEN
      esatt = p / ( rdv/r + o_m_rdv )
      z = LOG(esatt / c1es)
      td = (z*c4les-c3les*tmelt) / (z-c3les)
    ELSE
      td = c4les
    END IF
    
  END FUNCTION dewpoint_water

  
  !!-------------------------------------------------------------------------------------------------
  !!
  !> Wrapper routine to get the 3D radar reflectivity field depending on the microphysics scheme
  !!  and store it in p_diag%dbz3d(:,:,:)
  !!
  SUBROUTINE compute_field_dbz3d_lin(jg, ptr_patch, p_prog,  p_prog_rcf, p_diag, prm_diag, dbz3d_lin, lacc)

    INTEGER, INTENT(in)  :: jg
    ! patch on which computation is performed:
    TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
    ! nonhydrostatic state
    TYPE(t_nh_prog), INTENT(IN)       :: p_prog            !< at timelevel nnow(jg)
    TYPE(t_nh_prog), INTENT(IN)       :: p_prog_rcf        !< at timelevel nnow_rcf(jg)
    TYPE(t_nh_diag), INTENT(IN)       :: p_diag
    TYPE(t_nwp_phy_diag), INTENT(IN)  :: prm_diag
    REAL(wp),        INTENT(OUT)      :: dbz3d_lin(:,:,:)  !< reflectivity in mm^6/m^3
    LOGICAL,    OPTIONAL, INTENT(IN)  :: lacc              !< initialization flag

    ! local variables
    CHARACTER(len=*), PARAMETER :: routine = modname//': compute_field_dbz3d_lin'
    REAL(wp) :: rho, qnc_s(nproma,ptr_patch%nblks_c), cloud_num
    INTEGER  :: i_rlstart, i_rlend, i_startblk, i_endblk, i_startidx, i_endidx, i_startidx_1, i_endidx_2, &
         &      jc, jk, jb, ilow, iup, jlow, jup, klow, kup, itype_gscp_emvo, isnow_n0temp

    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: Tmax_i, Tmax_s, Tmax_g, Tmax_h, Tmin_g, Tmin_h
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:), TARGET :: dummy0
    REAL(wp), POINTER, DIMENSION(:,:,:)   :: t, p, rho_tot, qc, qr, qi, qs, qg, qh, qnc, qnr, qni, qns, qng, qnh, qgl, qhl

    LOGICAL :: lzacc             ! OpenACC flag
    CALL set_acc_host_or_device(lzacc, lacc)

#ifdef HAVE_RADARFWO
    IF ( synradar_meta%itype_refl == 4 ) THEN
#endif

      !=========================================================================
      ! Simple Rayleigh approximation where melting particles just get the refractive index of water
      ! instead of a proper EMA. This is the default.
      !=========================================================================

      ! without halo or boundary points:
      i_rlstart = grf_bdywidth_c + 1
      i_rlend   = min_rlcell_int

      i_startblk = ptr_patch%cells%start_block( i_rlstart )
      i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

      ! Determine first and last index for first and last block, respectively
      CALL get_indices_c( ptr_patch, i_startblk, i_startblk, i_endblk, i_startidx_1, i_endidx, i_rlstart, i_rlend)
      CALL get_indices_c( ptr_patch, i_endblk, i_startblk, i_endblk, i_startidx, i_endidx_2, i_rlstart, i_rlend)

      SELECT CASE ( atm_phy_nwp_config(jg)%inwp_gscp )
      CASE ( 1,3 )

        CALL get_cloud_number(cloud_num)
        IF (atm_phy_nwp_config(jg)%icpl_aero_gscp == 2) THEN
          ! Not yet implemented in microphysics! We give a dummy value here.
          qnc_s(:,:) = cloud_num               ! 1/kg
        ELSE IF ( ANY ( atm_phy_nwp_config(jg)%icpl_aero_gscp == (/1, 3/) ) ) THEN
          qnc_s(:,:) = prm_diag%cloud_num(:,:) ! neglect difference of 1/m^3 and 1/kg for this near-surface value
        ELSE
          qnc_s(:,:) = cloud_num               ! 1/kg
        END IF

        CALL compute_field_dbz_1mom( npr       = nproma,                           &
             nlev      = ptr_patch%nlev,                   &
             nblks     = ptr_patch%nblks_c,                &
             startblk  = i_startblk,                       &
             endblk    = i_endblk,                         &
             jk_start  = kstart_moist(jg),                 &
             startidx1 = i_startidx_1,                     &
             endidx2   = i_endidx_2,                       &
             lmessage_light = (msg_level > 12 .AND. my_process_is_mpi_workroot()), &
             lmessage_full  = (msg_level > 15),            &
             my_id_for_message = get_my_mpi_work_id(),     &
             rho_w     = rhoh2o,                           &
             rho_ice   = rhoice,                           &
             K_w       = K_w_0,                            &
             K_ice     = K_i_0,                            &
             T_melt    = Tmelt,                            &
             igscp     = atm_phy_nwp_config(jg)%inwp_gscp, &
             q_crit_radar = 1e-8_wp,                       &
             T         = p_diag%temp(:,:,:),               &
             rho       = p_prog%rho(:,:,:),                &
             q_cloud   = p_prog_rcf%tracer(:,:,:,iqc),     &
             q_ice     = p_prog_rcf%tracer(:,:,:,iqi),     &
             q_rain    = p_prog_rcf%tracer(:,:,:,iqr),     &
             q_snow    = p_prog_rcf%tracer(:,:,:,iqs),     &
             n_cloud_s = qnc_s(:,:),                       &  ! 1/kg
             z_radar   = dbz3d_lin(:,:,:),                 &
             lacc      = lzacc                             )

      CASE ( 2 )

        CALL get_cloud_number(cloud_num)
        IF (atm_phy_nwp_config(jg)%icpl_aero_gscp == 2) THEN
          ! Not yet implemented in microphysics! We give a dummy value here.
          qnc_s(:,:) = cloud_num               ! 1/kg
        ELSE IF (atm_phy_nwp_config(jg)%icpl_aero_gscp == 1) THEN
          qnc_s(:,:) = prm_diag%cloud_num(:,:) ! neglect difference of 1/m^3 and 1/kg for this near-surface value
        ELSE
          qnc_s(:,:) = cloud_num               ! 1/kg
        END IF

        CALL compute_field_dbz_1mom( npr       = nproma,                           &
             nlev      = ptr_patch%nlev,                   &
             nblks     = ptr_patch%nblks_c,                &
             startblk  = i_startblk,                       &
             endblk    = i_endblk,                         &
             jk_start  = kstart_moist(jg),                 &
             startidx1 = i_startidx_1,                     &
             endidx2   = i_endidx_2,                       &
             lmessage_light = (msg_level > 12 .AND. my_process_is_mpi_workroot()), &
             lmessage_full  = (msg_level > 15),            &
             my_id_for_message = get_my_mpi_work_id(),     &
             rho_w     = rhoh2o,                           &
             rho_ice   = rhoice,                           &
             K_w       = K_w_0,                            &
             K_ice     = K_i_0,                            &
             T_melt    = Tmelt,                            &
             igscp     = atm_phy_nwp_config(jg)%inwp_gscp, &
             q_crit_radar = 1e-8_wp,                       &
             T         = p_diag%temp(:,:,:),               &
             rho       = p_prog%rho(:,:,:),                &
             q_cloud   = p_prog_rcf%tracer(:,:,:,iqc),     &
             q_ice     = p_prog_rcf%tracer(:,:,:,iqi),     &
             q_rain    = p_prog_rcf%tracer(:,:,:,iqr),     &
             q_snow    = p_prog_rcf%tracer(:,:,:,iqs),     &
             q_graupel = p_prog_rcf%tracer(:,:,:,iqg),     &
             n_cloud_s = qnc_s(:,:),                       &  ! 1/kg
             z_radar   = dbz3d_lin(:,:,:),                 &
             lacc      = lzacc                             )

      CASE ( 4, 5, 6, 8 )

        CALL compute_field_dbz_2mom( npr       = nproma,                           &
             nlev      = ptr_patch%nlev,                   &
             nblks     = ptr_patch%nblks_c,                &
             startblk  = i_startblk,                       &
             endblk    = i_endblk,                         &
             jk_start  = kstart_moist(jg),                 &
             startidx1 = i_startidx_1,                     &
             endidx2   = i_endidx_2,                       &
             lmessage_light = (msg_level > 12 .AND. my_process_is_mpi_workroot()), &
             lmessage_full  = (msg_level > 15),            &
             my_id_for_message = get_my_mpi_work_id(),     &
             rho_w     = rhoh2o,                           &
             rho_ice   = rhoice,                           &
             K_w       = K_w_0,                            &
             K_ice     = K_i_0,                            &
             T_melt    = Tmelt,                            &
             q_crit_radar = 1e-8_wp,                       &
             luse_mu_Dm_rain = atm_phy_nwp_config(jg)%cfg_2mom%luse_mu_Dm_rain, &
             T         = p_diag%temp(:,:,:),               &
             rho       = p_prog%rho(:,:,:),                &
             q_cloud   = p_prog_rcf%tracer(:,:,:,iqc),     &
             q_ice     = p_prog_rcf%tracer(:,:,:,iqi),     &
             q_rain    = p_prog_rcf%tracer(:,:,:,iqr),     &
             q_snow    = p_prog_rcf%tracer(:,:,:,iqs),     &
             q_graupel = p_prog_rcf%tracer(:,:,:,iqg),     &
             q_hail    = p_prog_rcf%tracer(:,:,:,iqh),     &
             n_cloud   = p_prog_rcf%tracer(:,:,:,iqnc),    &
             n_ice     = p_prog_rcf%tracer(:,:,:,iqni),    &
             n_rain    = p_prog_rcf%tracer(:,:,:,iqnr),    &
             n_snow    = p_prog_rcf%tracer(:,:,:,iqns),    &
             n_graupel = p_prog_rcf%tracer(:,:,:,iqng),    &
             n_hail    = p_prog_rcf%tracer(:,:,:,iqnh),    &
             z_radar   = dbz3d_lin(:,:,:),                 &
             lacc      = lzacc                             )


      CASE ( 7 )
#ifdef _OPENACC
        CALL finish(routine, 'compute_field_dbz_2mom is supported by OpenACC, but never tested.')
#endif
        CALL compute_field_dbz_2mom( npr       = nproma,                           &
             nlev      = ptr_patch%nlev,                   &
             nblks     = ptr_patch%nblks_c,                &
             startblk  = i_startblk,                       &
             endblk    = i_endblk,                         &
             jk_start  = kstart_moist(jg),                 &
             startidx1 = i_startidx_1,                     &
             endidx2   = i_endidx_2,                       &
             lmessage_light = (msg_level > 12 .AND. my_process_is_mpi_workroot()), &
             lmessage_full  = (msg_level > 15),            &
             my_id_for_message = get_my_mpi_work_id(),     &
             rho_w     = rhoh2o,                           &
             rho_ice   = rhoice,                           &
             K_w       = K_w_0,                            &
             K_ice     = K_i_0,                            &
             T_melt    = Tmelt,                            &
             q_crit_radar = 1e-8_wp,                       &
             luse_mu_Dm_rain = atm_phy_nwp_config(jg)%cfg_2mom%luse_mu_Dm_rain, &
             T         = p_diag%temp(:,:,:),               &
             rho       = p_prog%rho(:,:,:),                &
             q_cloud   = p_prog_rcf%tracer(:,:,:,iqc),     &
             q_ice     = p_prog_rcf%tracer(:,:,:,iqi),     &
             q_rain    = p_prog_rcf%tracer(:,:,:,iqr),     &
             q_snow    = p_prog_rcf%tracer(:,:,:,iqs),     &
             q_graupel = p_prog_rcf%tracer(:,:,:,iqg),     &
             q_hail    = p_prog_rcf%tracer(:,:,:,iqh),     &
             n_cloud   = p_prog_rcf%tracer(:,:,:,iqnc),    &
             n_ice     = p_prog_rcf%tracer(:,:,:,iqni),    &
             n_rain    = p_prog_rcf%tracer(:,:,:,iqnr),    &
             n_snow    = p_prog_rcf%tracer(:,:,:,iqns),    &
             n_graupel = p_prog_rcf%tracer(:,:,:,iqng),    &
             n_hail    = p_prog_rcf%tracer(:,:,:,iqnh),    &
             ql_graupel= p_prog_rcf%tracer(:,:,:,iqgl),    &
             ql_hail   = p_prog_rcf%tracer(:,:,:,iqhl),    &
             z_radar   = dbz3d_lin(:,:,:),                 &
             lacc      = lzacc                             )


      CASE DEFAULT

        CALL finish( routine,  &
             &     "dbz3d-computation not available for this microphysics scheme! Available for inwp_gscp=1,2,4,5,6 or 7" )

      END SELECT


#ifdef HAVE_RADARFWO
    ELSE   ! synradar_meta%itype_refl == 1, 3, 5, 6

      !=========================================================================
      ! Mie-scattering or Rayleigh-Oguchi-Approximation from EMVORADO:
      !=========================================================================

      ! .. Reduced domain of this PE, i.e., excluding interior
      !    boundary lines but including boudary halo at the outer
      !    model domain boundaries:
      ilow = 1
      iup  = nproma   ! the few halo/boundary points in the first and last block are not explicitly excluded here for simplicity
      jlow = kstart_moist(jg)
      jup  = ptr_patch % nlev
      klow = ptr_patch % cells % start_block(grf_bdywidth_c+1)
      kup  = ptr_patch % cells % end_block(min_rlcell_int)

      ! .. Set module switch ldebug_dbz from EMVORADO according to ICONs msg_level. Has to be the same on all workers,
      !     so there is no "light" version of the debug output, only all or nothing:
      ldebug_dbz = (msg_level > 15)

      SELECT CASE ( atm_phy_nwp_config(jg)%inwp_gscp )
      CASE ( 1, 2, 3 )

        t   => p_diag%temp(:,:,:)
        p   => p_diag%pres(:,:,:)
        rho_tot => p_prog%rho(:,:,:)
        qc  => p_prog_rcf%tracer(:,:,:,iqc)
        qr  => p_prog_rcf%tracer(:,:,:,iqr)
        qi  => p_prog_rcf%tracer(:,:,:,iqi)
        qs  => p_prog_rcf%tracer(:,:,:,iqs)
        IF (atm_phy_nwp_config(jg)%lhave_graupel) THEN 
          qg => p_prog_rcf%tracer(:,:,:,iqg)
        ELSE
          ALLOCATE(dummy0(nproma,ptr_patch%nlev,ptr_patch%nblks_c))
          dummy0 = 0.0_wp
          qg => dummy0(:,:,:)
        END IF

        IF (atm_phy_nwp_config(jg)%inwp_gscp == 1) THEN
          itype_gscp_emvo = 140 ! "140" is the corresponding itype_gscp in COSMO and EMVORADO
        ELSE
          itype_gscp_emvo = 150 ! "150" is the corresponding itype_gscp in COSMO and EMVORADO
        END IF
        CALL init_1mom_types(itype_gscp_loc=itype_gscp_emvo, rho_w=rhoh2o)

        ALLOCATE (Tmax_i(nproma,ptr_patch%nblks_c), Tmax_s(nproma,ptr_patch%nblks_c), Tmax_g(nproma,ptr_patch%nblks_c))
        CALL initialize_tmax_atomic_1mom(qx=qi, t=t, neigh=0.0_wp, qthresh=synradar_meta%qthresh_i, &
             &                           Tmax_min=synradar_meta%Tmax_min_i, Tmax_max=synradar_meta%Tmax_max_i, Tmax_x=Tmax_i)
        CALL initialize_tmax_atomic_1mom(qx=qs, t=t, neigh=0.0_wp, qthresh=synradar_meta%qthresh_s, &
             &                           Tmax_min=synradar_meta%Tmax_min_s, Tmax_max=synradar_meta%Tmax_max_s, Tmax_x=Tmax_s)
        CALL initialize_tmax_atomic_1mom(qx=qg, t=t, neigh=0.0_wp, qthresh=synradar_meta%qthresh_g, &
             &                           Tmax_min=synradar_meta%Tmax_min_g, Tmax_max=synradar_meta%Tmax_max_g, Tmax_x=Tmax_g)

        ALLOCATE (Tmin_g(nproma,ptr_patch%nblks_c))
        IF (synradar_meta%ldynamic_wetgrowth_gh .AND. synradar_meta%Tmeltbegin_g < T0C_emvorado) THEN
          CALL initialize_tmin_atomic_1mom(hydrotype='graupel', qx=qg, t=t, p=p, ql=qc+qr, qf=qi+qs, rho=rho_tot, &
                                           Tmin_min=synradar_meta%Tmeltbegin_g, Tmin_x=Tmin_g)
        ELSE
          Tmin_g = synradar_meta%Tmeltbegin_g
        END IF

        CALL get_snow_temperature(isnow_n0temp)
        SELECT CASE ( synradar_meta%itype_refl )
        CASE ( 1, 5, 6 )
          ! Mie- or T-matrix scattering from EMVORADO:
#ifdef _OPENACC
          CALL finish(routine, 'radar_mie_1mom_vec is not supported by OpenACC.')
#endif
          CALL radar_mie_1mom_vec( &
               myproc               = get_my_mpi_work_id(), &
               lambda_radar         = synradar_meta%lambda_radar, &
               itype_gscp_fwo       = itype_gscp_emvo, &
               itype_refl           = synradar_meta%itype_refl, &
               luse_tmatrix         = (synradar_meta%itype_refl >= 5), &
               ldo_nonsphere        = (synradar_meta%itype_refl == 5), &
               isnow_n0temp         = isnow_n0temp, &
               igraupel_type        = synradar_meta%igraupel_type, &
               itype_Dref_fmelt     = synradar_meta%itype_Dref_fmelt, &
               ctype_dryice         = synradar_meta%ctype_dryice_mie, &
               ctype_wetice         = synradar_meta%ctype_wetice_mie, &
               ctype_drysnow        = synradar_meta%ctype_drysnow_mie, &
               ctype_wetsnow        = synradar_meta%ctype_wetsnow_mie, &
               ctype_drygraupel     = synradar_meta%ctype_drygraupel_mie, &
               ctype_wetgraupel     = synradar_meta%ctype_wetgraupel_mie, &
               ldynamic_wetgrowth_gh= synradar_meta%ldynamic_wetgrowth_gh, &
               Tmeltbegin_i         = synradar_meta%Tmeltbegin_i, &
               meltdegTmin_i        = synradar_meta%meltdegTmin_i, &
               Tmax_min_i           = synradar_meta%Tmax_min_i, &
               Tmax_max_i           = synradar_meta%Tmax_max_i, &
               Tmeltbegin_s         = synradar_meta%Tmeltbegin_s, &
               meltdegTmin_s        = synradar_meta%meltdegTmin_s, &
               Tmax_min_s           = synradar_meta%Tmax_min_s, &
               Tmax_max_s           = synradar_meta%Tmax_max_s, &
               Tmeltbegin_g         = synradar_meta%Tmeltbegin_g, &
               meltdegTmin_g        = synradar_meta%meltdegTmin_g, &
               Tmax_min_g           = synradar_meta%Tmax_min_g, &
               Tmax_max_g           = synradar_meta%Tmax_max_g, &
               pMPr                 = synradar_meta%polMP_r, &
               pMPi                 = synradar_meta%polMP_i, &
               pMPs                 = synradar_meta%polMP_s, &
               pMPg                 = synradar_meta%polMP_g, &
               rho                  = rho_tot(:,:,:), &
               t                    = t(:,:,:), &
               qc                   = qc(:,:,:), &
               qr                   = qr(:,:,:), &
               qi                   = qi(:,:,:), &
               qs                   = qs(:,:,:), &
               qg                   = qg(:,:,:), &
               Tmax_i               = Tmax_i(:,:), &
               Tmax_s               = Tmax_s(:,:), &
               Tmax_g               = Tmax_g(:,:), &
               Tmin_g               = Tmin_g(:,:), &
               ilow=ilow, iup=iup, jlow=jlow, jup=jup, klow=klow, kup=kup, &
               lalloc_qi            = .TRUE., &
               lalloc_qs            = .TRUE., &
               lalloc_qg            = atm_phy_nwp_config(jg)%lhave_graupel, &
               llookup              = synradar_meta%llookup_mie, &
               impipar_lookupgen    = 2, &
               pe_start             = proc0_shift,   & ! Start-PE of the gang which computes the lookup tables, numbering within the work-communicator
               pe_end               = get_my_mpi_work_comm_size()-1, &  ! End-PE of the gang. Can be at most the number of work PEs minus 1
               linterp_mode_dualpol = (synradar_meta%itype_refl >= 5), &
               ydir_lookup_read     = TRIM(ydir_mielookup_read), &
               ydir_lookup_write    = TRIM(ydir_mielookup_write), &
               ext_tune_fac_pure    = synradar_meta%ext_tune_fac_pure, &
               ext_tune_fac_melt    = synradar_meta%ext_tune_fac_melt, &
               zh_radar             = dbz3d_lin(:,:,:), &
               lhydrom_choice_testing = synradar_meta%lhydrom_choice_testing &
               )

        CASE ( 3 )
#ifdef _OPENACC
          CALL finish(routine, 'radar_rayleigh_oguchi_1mom_vec is not supported by OpenACC.')
#endif
          CALL radar_rayleigh_oguchi_1mom_vec( &
               myproc         = get_my_mpi_work_id(), &
               lambda_radar   = synradar_meta%lambda_radar, &
               itype_gscp_fwo = itype_gscp_emvo, &
               isnow_n0temp   = isnow_n0temp, &
               Tmeltbegin_i   = synradar_meta%Tmeltbegin_i, &
               meltdegTmin_i  = synradar_meta%meltdegTmin_i, &
               Tmeltbegin_s   = synradar_meta%Tmeltbegin_s, &
               meltdegTmin_s  = synradar_meta%meltdegTmin_s, &
               Tmeltbegin_g   = synradar_meta%Tmeltbegin_g, &
               meltdegTmin_g  = synradar_meta%meltdegTmin_g, &
               rho            = rho_tot(:,:,:), &
               t              = t(:,:,:), &
               qc             = qc(:,:,:), &
               qr             = qr(:,:,:), &
               qi             = qi(:,:,:), &
               qs             = qs(:,:,:), &
               qg             = qg(:,:,:), &
               Tmax_i         = Tmax_i(:,:), &
               Tmax_s         = Tmax_s(:,:), &
               Tmax_g         = Tmax_g(:,:), &
               Tmin_g         = Tmin_g(:,:), &
               ilow=ilow, iup=iup, jlow=jlow, jup=jup, klow=klow, kup=kup, &
               lalloc_qi      = .TRUE., &
               lalloc_qs      = .TRUE., &
               lalloc_qg      = atm_phy_nwp_config(jg)%lhave_graupel, &
               zh_radar       = dbz3d_lin(:,:,:), &
               lhydrom_choice_testing = synradar_meta%lhydrom_choice_testing &
               )

        CASE default

          CALL finish( routine,  &
               &     "dbz3d-computation not available for this synradar_meta%itype_refl! Available are options 1, 3, 4, 5 and 6" )

        END SELECT

        DEALLOCATE (Tmax_i, Tmax_s, Tmax_g, Tmin_g)
        IF (ALLOCATED(dummy0)) DEALLOCATE(dummy0)

      CASE ( 4, 5, 6, 7, 8)

        t   => p_diag%temp(:,:,:)
        p   => p_diag%pres(:,:,:)
        rho_tot => p_prog%rho(:,:,:)
        qc  => p_prog_rcf%tracer(:,:,:,iqc)
        qr  => p_prog_rcf%tracer(:,:,:,iqr)
        qi  => p_prog_rcf%tracer(:,:,:,iqi)
        qs  => p_prog_rcf%tracer(:,:,:,iqs)
        qg  => p_prog_rcf%tracer(:,:,:,iqg)
        qh  => p_prog_rcf%tracer(:,:,:,iqh)
        qnc => p_prog_rcf%tracer(:,:,:,iqnc)
        qnr => p_prog_rcf%tracer(:,:,:,iqnr)
        qni => p_prog_rcf%tracer(:,:,:,iqni)
        qns => p_prog_rcf%tracer(:,:,:,iqns)
        qng => p_prog_rcf%tracer(:,:,:,iqng)
        qnh => p_prog_rcf%tracer(:,:,:,iqnh)
        IF (atm_phy_nwp_config(jg)%inwp_gscp == 7) THEN
          qgl => p_prog_rcf%tracer(:,:,:,iqgl)
          qhl => p_prog_rcf%tracer(:,:,:,iqhl)
        ELSE
          ALLOCATE(dummy0(nproma,ptr_patch%nlev,ptr_patch%nblks_c))
          dummy0 = 0.0_wp
          qgl => dummy0(:,:,:)
          qhl => dummy0(:,:,:)
        END IF

        CALL init_2mom_types()

        ALLOCATE ( Tmax_i(nproma,ptr_patch%nblks_c), Tmax_s(nproma,ptr_patch%nblks_c), &
             Tmax_g(nproma,ptr_patch%nblks_c), Tmax_h(nproma,ptr_patch%nblks_c) )
        CALL initialize_tmax_atomic_2mom( qx=qi, qnx=qni, t=t, neigh=0.0_wp, &
             &                            qthresh=synradar_meta%qthresh_i, qnthresh=synradar_meta%qnthresh_i, &
             &                            Tmax_min=synradar_meta%Tmax_min_i, Tmax_max=synradar_meta%Tmax_max_i, Tmax_x=Tmax_i )
        CALL initialize_tmax_atomic_2mom( qx=qs, qnx=qns, t=t, neigh=0.0_wp, &
             &                            qthresh=synradar_meta%qthresh_s, qnthresh=synradar_meta%qnthresh_s, &
             &                            Tmax_min=synradar_meta%Tmax_min_s, Tmax_max=synradar_meta%Tmax_max_s, Tmax_x=Tmax_s )
        CALL initialize_tmax_atomic_2mom( qx=qg+qgl, qnx=qng, t=t, neigh=0.0_wp, &
             &                            qthresh=synradar_meta%qthresh_g, qnthresh=synradar_meta%qnthresh_g, &
             &                            Tmax_min=synradar_meta%Tmax_min_g, Tmax_max=synradar_meta%Tmax_max_g, Tmax_x=Tmax_g )
        CALL initialize_tmax_atomic_2mom( qx=qh+qhl, qnx=qnh, t=t, neigh=0.0_wp, &
             &                            qthresh=synradar_meta%qthresh_h, qnthresh=synradar_meta%qnthresh_h, &
             &                            Tmax_min=synradar_meta%Tmax_min_h, Tmax_max=synradar_meta%Tmax_max_h, Tmax_x=Tmax_h )

        ALLOCATE (Tmin_g(nproma,ptr_patch%nblks_c), Tmin_h(nproma,ptr_patch%nblks_c))
        IF (synradar_meta%ldynamic_wetgrowth_gh .AND. synradar_meta%Tmeltbegin_g < T0C_emvorado) THEN
          CALL initialize_tmin_atomic_2mom(hydrotype='graupel', qx=qg+qgl, qnx=qng, t=t, p=p, ql=qc+qr, qf=qi+qs, rho=rho_tot, &
                                           Tmin_min=synradar_meta%Tmeltbegin_g, Tmin_x=Tmin_g)
        ELSE
          Tmin_g = synradar_meta%Tmeltbegin_g
        END IF
        IF (synradar_meta%ldynamic_wetgrowth_gh .AND. synradar_meta%Tmeltbegin_h < T0C_emvorado) THEN
          CALL initialize_tmin_atomic_2mom(hydrotype='hail', qx=qh+qhl, qnx=qnh, t=t, p=p, ql=qc+qr, qf=qi+qs, rho=rho_tot, &
                                           Tmin_min=synradar_meta%Tmeltbegin_h, Tmin_x=Tmin_h)
        ELSE
          Tmin_h = synradar_meta%Tmeltbegin_h
        END IF

        SELECT CASE ( synradar_meta%itype_refl )
        CASE ( 1, 5, 6 )
          ! Mie-scattering from EMVORADO:
#ifdef _OPENACC
          CALL finish(routine, 'radar_mie_2mom_vec is not supported by OpenACC.')
#endif
          CALL radar_mie_2mom_vec( &
               myproc            = get_my_mpi_work_id(), &
               lambda_radar      = synradar_meta%lambda_radar, &
               itype_gscp_fwo    = 260, &
               itype_refl        = synradar_meta%itype_refl, &
               luse_tmatrix      = (synradar_meta%itype_refl >= 5), &
               ldo_nonsphere     = (synradar_meta%itype_refl == 5), &
               igraupel_type     = synradar_meta%igraupel_type, &
               itype_Dref_fmelt  = synradar_meta%itype_Dref_fmelt, &
               ctype_dryice      = synradar_meta%ctype_dryice_mie, &
               ctype_wetice      = synradar_meta%ctype_wetice_mie, &
               ctype_drysnow     = synradar_meta%ctype_drysnow_mie, &
               ctype_wetsnow     = synradar_meta%ctype_wetsnow_mie, &
               ctype_drygraupel  = synradar_meta%ctype_drygraupel_mie, &
               ctype_wetgraupel  = synradar_meta%ctype_wetgraupel_mie, &
               ctype_dryhail     = synradar_meta%ctype_dryhail_mie, &
               ctype_wethail     = synradar_meta%ctype_wethail_mie, &
               ldynamic_wetgrowth_gh= synradar_meta%ldynamic_wetgrowth_gh, &
               Tmeltbegin_i      = synradar_meta%Tmeltbegin_i, &
               meltdegTmin_i     = synradar_meta%meltdegTmin_i, &
               Tmax_min_i        = synradar_meta%Tmax_min_i, &
               Tmax_max_i        = synradar_meta%Tmax_max_i, &
               Tmeltbegin_s      = synradar_meta%Tmeltbegin_s, &
               meltdegTmin_s     = synradar_meta%meltdegTmin_s, &
               Tmax_min_s        = synradar_meta%Tmax_min_s, &
               Tmax_max_s        = synradar_meta%Tmax_max_s, &
               Tmeltbegin_g      = synradar_meta%Tmeltbegin_g, &
               meltdegTmin_g     = synradar_meta%meltdegTmin_g, &
               Tmax_min_g        = synradar_meta%Tmax_min_g, &
               Tmax_max_g        = synradar_meta%Tmax_max_g, &
               Tmeltbegin_h      = synradar_meta%Tmeltbegin_h, &
               meltdegTmin_h     = synradar_meta%meltdegTmin_h, &
               Tmax_min_h        = synradar_meta%Tmax_min_h, &
               Tmax_max_h        = synradar_meta%Tmax_max_h, &
               pMPr              = synradar_meta%polMP_r, &
               pMPi              = synradar_meta%polMP_i, &
               pMPs              = synradar_meta%polMP_s, &
               pMPg              = synradar_meta%polMP_g, &
               pMPh              = synradar_meta%polMP_g, &
               rho               = rho_tot(:,:,:), &
               t                 = t(:,:,:), &
               qc                = qc(:,:,:), &
               qr                = qr(:,:,:), &
               qi                = qi(:,:,:), &
               qs                = qs(:,:,:), &
               qg                = qg(:,:,:), &
               qh                = qh(:,:,:), &
               qnc               = qnc(:,:,:), &
               qnr               = qnr(:,:,:), &
               qni               = qni(:,:,:), &
               qns               = qns(:,:,:), &
               qng               = qng(:,:,:), &
               qnh               = qnh(:,:,:), &
               qgl               = qgl(:,:,:), &
               qhl               = qhl(:,:,:), &
               Tmax_i            = Tmax_i(:,:), &
               Tmax_s            = Tmax_s(:,:), &
               Tmax_g            = Tmax_g(:,:), &
               Tmax_h            = Tmax_h(:,:), &
               Tmin_g            = Tmin_g(:,:), &
               Tmin_h            = Tmin_h(:,:), &
               ilow=ilow, iup=iup, jlow=jlow, jup=jup, klow=klow, kup=kup, &
               lalloc_qi         = .TRUE., &
               lalloc_qs         = .TRUE., &
               lalloc_qg         = .TRUE., &
               lalloc_qh         = .TRUE., &
               llookup           = synradar_meta%llookup_mie, &
               impipar_lookupgen = 2, &
               pe_start          = proc0_shift,   & ! Start-PE of the gang which computes the lookup tables, numbering within the work-communicator
               pe_end            = get_my_mpi_work_comm_size()-1, &  ! End-PE of the gang. Can be at most the number of work PEs minus 1
               linterp_mode_dualpol = (synradar_meta%itype_refl >= 5), &
               luse_muD_relation_rain  = atm_phy_nwp_config(jg)%cfg_2mom%luse_mu_Dm_rain, &
               ydir_lookup_read  = TRIM(ydir_mielookup_read), &
               ydir_lookup_write = TRIM(ydir_mielookup_write), &
               ext_tune_fac_pure    = synradar_meta%ext_tune_fac_pure, &
               ext_tune_fac_melt    = synradar_meta%ext_tune_fac_melt, &
               zh_radar          = dbz3d_lin(:,:,:), &
               lhydrom_choice_testing = synradar_meta%lhydrom_choice_testing &
               )

        CASE ( 3 )
#ifdef _OPENACC
          CALL finish(routine, 'radar_rayleigh_oguchi_2mom_vec is not supported by OpenACC.')
#endif
          CALL radar_rayleigh_oguchi_2mom_vec( &
               myproc         = get_my_mpi_work_id(), &
               lambda_radar   = synradar_meta%lambda_radar, &
               Tmeltbegin_i   = synradar_meta%Tmeltbegin_i, &
               meltdegTmin_i  = synradar_meta%meltdegTmin_i, &
               Tmeltbegin_s   = synradar_meta%Tmeltbegin_s, &
               meltdegTmin_s  = synradar_meta%meltdegTmin_s, &
               Tmeltbegin_g   = synradar_meta%Tmeltbegin_g, &
               meltdegTmin_g  = synradar_meta%meltdegTmin_g, &
               Tmeltbegin_h   = synradar_meta%Tmeltbegin_h, &
               meltdegTmin_h  = synradar_meta%meltdegTmin_h, &
               rho            = rho_tot(:,:,:), &
               t              = t(:,:,:), &
               qc             = qc(:,:,:), &
               qr             = qr(:,:,:), &
               qi             = qi(:,:,:), &
               qs             = qs(:,:,:), &
               qg             = qg(:,:,:), &
               qh             = qh(:,:,:), &
               qnc            = qnc(:,:,:), &
               qnr            = qnr(:,:,:), &
               qni            = qni(:,:,:), &
               qns            = qns(:,:,:), &
               qng            = qng(:,:,:), &
               qnh            = qnh(:,:,:), &
               qgl            = qgl(:,:,:), &
               qhl            = qhl(:,:,:), &
               Tmax_i         = Tmax_i(:,:), &
               Tmax_s         = Tmax_s(:,:), &
               Tmax_g         = Tmax_g(:,:), &
               Tmax_h         = Tmax_h(:,:), &
               Tmin_g         = Tmin_g(:,:), &
               Tmin_h         = Tmin_h(:,:), &
               ilow=ilow, iup=iup, jlow=jlow, jup=jup, klow=klow, kup=kup, &
               lalloc_qi      = .TRUE., &
               lalloc_qs      = .TRUE., &
               lalloc_qg      = .TRUE., &
               lalloc_qh      = .TRUE., &
               luse_muD_relation_rain  = atm_phy_nwp_config(jg)%cfg_2mom%luse_mu_Dm_rain, &
               zh_radar       = dbz3d_lin(:,:,:), &
               lhydrom_choice_testing = synradar_meta%lhydrom_choice_testing &
               )

        CASE default

          CALL finish( routine,  &
               &     "dbz3d-computation not available for this synradar_meta%itype_refl! Available are options 1, 3, 4, 5, and 6" )

        END SELECT

        DEALLOCATE (Tmax_i, Tmax_s, Tmax_g, Tmax_h, Tmin_g, Tmin_h)
        IF (ALLOCATED(dummy0)) DEALLOCATE(dummy0)

      CASE DEFAULT

        CALL finish( routine,  &
             &     "dbz3d-computation not available for this microphysics scheme! Available for inwp_gscp=1,2,4,5,6 or 7" )

      END SELECT

    END IF   ! synradar_meta%itype_refl == 4
#endif

  END SUBROUTINE compute_field_dbz3d_lin


  !>
  !! Compute column maximum reflectivity from dbz3d_lin
  !!
  SUBROUTINE compute_field_dbzcmax( ptr_patch, jg, dbz3d_lin, dbz_cmax, lacc )

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch        !< patch on which computation is performed
    INTEGER,              INTENT(IN)  :: jg               ! domain ID of main grid
    REAL(wp),             INTENT(IN)  :: dbz3d_lin(:,:,:) !< reflectivity in mm^6/m^3

    REAL(wp),             INTENT(OUT) :: dbz_cmax(:,:)  !< output variable, dim: (nproma,nblks_c)
    
    LOGICAL,    OPTIONAL, INTENT(IN)  :: lacc           !< initialization flag

    REAL(wp) :: most_negative_value

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jk, jc

    LOGICAL :: lzacc             ! OpenACC flag
    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA PRESENT(dbz3d_lin, dbz_cmax, kstart_moist(jg:jg), ptr_patch) IF(lzacc)

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    most_negative_value = -HUGE(1.0_wp)

!$OMP PARALLEL
    CALL init(dbz_cmax(:,i_startblk:i_endblk), most_negative_value, lacc=lzacc)
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO jk = kstart_moist(jg), ptr_patch%nlev
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx

          dbz_cmax(jc,jb) = MAX (dbz_cmax(jc,jb), dbz3d_lin(jc,jk,jb))

        END DO
      END DO
      !$ACC END PARALLEL

    END DO
    !$ACC WAIT(1)
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC END DATA

  END SUBROUTINE compute_field_dbzcmax

  !>
  !! Compute column maximum radar reflectivity from dbz3d_lin and maximize over time
  !!
  SUBROUTINE maximize_field_dbzctmax( ptr_patch, jg, dbz3d_lin, dbz_ctmax, lacc )

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch        !< patch on which computation is performed
    INTEGER,              INTENT(IN)  :: jg               ! domain ID of main grid
    REAL(wp),             INTENT(IN)  :: dbz3d_lin(:,:,:) !< reflectivity in mm^6/m^3

    REAL(wp),             INTENT(INOUT) :: dbz_ctmax(:,:)  !< input/output variable, dim: (nproma,nblks_c)

    LOGICAL,    OPTIONAL, INTENT(IN)  :: lacc             !< initialization flag

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jk, jc

    LOGICAL :: lzacc             ! OpenACC flag
    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA PRESENT(dbz3d_lin, dbz_ctmax, kstart_moist(jg:jg), ptr_patch) IF(lzacc)

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO jk = kstart_moist(jg), ptr_patch%nlev
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx

          dbz_ctmax(jc,jb) = MAX (dbz_ctmax(jc,jb), dbz3d_lin(jc,jk,jb))

        END DO
      END DO
      !$ACC END PARALLEL

    END DO
    !$ACC WAIT(1)
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC END DATA

  END SUBROUTINE maximize_field_dbzctmax

  !>
  !! Compute radar reflectivity around approx 850 hPa from dbz3d_lin.
  !!
  SUBROUTINE compute_field_dbz850( ptr_patch, k850, dbz3d_lin, dbz_850, lacc )

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch        !< patch on which computation is performed
    INTEGER,              INTENT(IN)  :: k850(:,:)        !< level index field indicating 850 hPa
    REAL(wp),             INTENT(IN)  :: dbz3d_lin(:,:,:) !< reflectivity in mm^6/m^3

    REAL(wp),             INTENT(OUT) :: dbz_850(:,:)  !< output variable, dim: (nproma,nblks_c)

    LOGICAL,    OPTIONAL, INTENT(IN)  :: lacc             !< initialization flag

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jk, jc

    LOGICAL :: lzacc             ! OpenACC flag
    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA PRESENT(dbz3d_lin, dbz_850, k850) IF(lzacc)

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
    CALL init(dbz_850(:,i_startblk:i_endblk), lacc=lzacc)
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)
      
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR PRIVATE(jk)
      DO jc = i_startidx, i_endidx

        jk = k850(jc,jb)

        ! Just take over the values from dbz3d_lin(:,:,:) in linear space. The conversion to dBZ
        ! will be done by a post_op ("post operation") right before output to file. See the
        ! corresponding add_var(..., post_op=(...) ) in mo_nwp_phy_state.f90!
        dbz_850(jc,jb) = dbz3d_lin(jc,jk,jb)

      END DO
      !$ACC END PARALLEL
      
    END DO
    !$ACC WAIT(1)
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC END DATA

  END SUBROUTINE compute_field_dbz850


  !>
  !! Compute layer maximum of radar reflectivity. This is for example to mimick the
  !! typical radar composites from observations. For this purpose, the maximum reflectivity
  !! in a layer from 500 to 2500 m AGL would be appropriate, to mimick the typical max composite method, where in spatially
  !! overlapping radar measuring circles the largest value measured by any station at a
  !! specific geographic location is taken. In today's radar networks in Europe, usually
  !! each location is sampled by at least 2 radar stations, and measuring heights vary between about
  !! 500 m AGL and 2500 m AGL, depending on the station height, the elevation angle used for
  !! the composite, and the local orography.
  !!
  SUBROUTINE compute_field_dbzlmx( ptr_patch, jg, z_agl_low, z_agl_up, p_metrics, dbz3d_lin, dbzlmx, lacc )

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch        !< patch on which computation is performed
    REAL(wp),             INTENT(in)  :: z_agl_low        !< lower height bound AGL for maximisation of reflectivity
    REAL(wp),             INTENT(in)  :: z_agl_up         !< upper height bound AGL for maximisation of reflectivity
    INTEGER,              INTENT(IN)  :: jg
    TYPE(t_nh_metrics),   INTENT(IN)  :: p_metrics
    REAL(wp),             INTENT(IN)  :: dbz3d_lin(:,:,:) !< reflectivity in mm^6/m^3

    REAL(wp),             INTENT(OUT) :: dbzlmx(:,:)  !< output variable, dim: (nproma,nblks_c)

    LOGICAL,    OPTIONAL, INTENT(IN)  :: lacc             !< initialization flag

    INTEGER  :: i_rlstart,  i_rlend
    INTEGER  :: i_startblk, i_endblk
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: jb, jk, jc, nlevp1
    REAL(wp) :: zml

    LOGICAL :: lzacc             ! OpenACC flag
    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA PRESENT(dbzlmx, dbz3d_lin, kstart_moist(jg:jg), ptr_patch, p_metrics, p_metrics%z_mc, p_metrics%z_ifc) &
    !$ACC   IF(lzacc)

    nlevp1 = ptr_patch%nlev+1
    
    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
    CALL init(dbzlmx(:,i_startblk:i_endblk), lacc=lzacc)
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,zml), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO jk = ptr_patch%nlev, kstart_moist(jg), -1
        !$ACC LOOP GANG VECTOR PRIVATE(zml)
        DO jc = i_startidx, i_endidx

          zml = p_metrics%z_mc(jc,jk,jb) - p_metrics%z_ifc(jc,nlevp1,jb)
          IF (z_agl_low <= zml .AND. zml <= z_agl_up) THEN
            ! Just take over the values from dbz3d_lin(:,:,:) in linear space. The conversion to dBZ
            ! will be done by a post_op ("post operation") right before output to file. See the
            ! corresponding add_var(..., post_op=(...) ) in mo_nwp_phy_state.f90!
            dbzlmx(jc,jb) = MAX(dbz3d_lin(jc,jk,jb), dbzlmx(jc,jb))
          END IF

        END DO
      END DO
      !$ACC END PARALLEL

    END DO
    !$ACC WAIT(1)
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC END DATA

  END SUBROUTINE compute_field_dbzlmx

  !>
  !! Compute ECHOTOPs in Pa from linear dbz3d_lin
  !!
  SUBROUTINE compute_field_echotop( ptr_patch, jg, p_diag, dbz3d_lin, echotop_p, lacc )

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch        !< patch on which computation is performed
    INTEGER,              INTENT(IN)  :: jg               !< domain ID of main grid
    TYPE(t_nh_diag),      INTENT(IN)  :: p_diag           !< type which contains the dbz3d_lin(:,:,:) field
    REAL(wp),             INTENT(IN)  :: dbz3d_lin(:,:,:) !< reflectivity in mm^6/m^3

    REAL(wp),             INTENT(INOUT) :: echotop_p(:,:,:)  !< input/output variable, dim: (nproma,nechotop,nblks_c)

    LOGICAL,    OPTIONAL, INTENT(IN)  :: lacc             !< initialization flag

    INTEGER               :: i_rlstart,  i_rlend
    INTEGER               :: i_startblk, i_endblk
    INTEGER               :: i_startidx, i_endidx
    INTEGER               :: jb, jc, jk, lev_etop
    INTEGER               :: jk_echotop(1:nproma)
    REAL(wp)              :: pechotop
    REAL(wp)              :: zthresh, zzthresh, zpA, zpB, zzdbzA,  zzdbzB

    REAL(wp), PARAMETER   :: repsilon = 1.0E8_wp*TINY(1.0_wp)  ! To prevent numerical division by 0 below

    LOGICAL :: lzacc             ! OpenACC flag
    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA CREATE(jk_echotop) PRESENT(dbz3d_lin, echotop_p, kstart_moist(jg:jg), ptr_patch, p_diag, p_diag%pres) &
    !$ACC   IF(lzacc)
    
    ! NOTE: pressure does not have to be recomputed/diagnosed here because this was already done when computing
    !       dbz3d_lin in the call to compute_field_dbz3d_lin() in mo_nh_stepping().
    
    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    DO lev_etop = 1, echotop_meta(jg)%nechotop

      zthresh  = 10.0_wp ** ( 0.1_wp*echotop_meta(jg)%dbzthresh(lev_etop))
      zzthresh = (zthresh+repsilon) ** 0.66666_wp  ! convert to something ~rain rate for interpolation,
                                                   ! take into account repsilon in the same way than in the EXP(LOG()) below

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,zzdbzA,zzdbzB,zpA,zpB,pechotop,jk_echotop), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                            i_startidx, i_endidx, i_rlstart, i_rlend)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        
        ! Find the model level just below the echotop:
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = i_startidx, i_endidx
          jk_echotop(jc)  = -999
        END DO

        !$ACC LOOP SEQ
        DO jk = MAX(kstart_moist(jg),2), ptr_patch%nlev
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = i_startidx, i_endidx
            IF ( jk_echotop(jc) < -900 .AND. dbz3d_lin(jc,jk,jb) >= zthresh ) THEN
              jk_echotop(jc) = jk
            END IF
          END DO
        END DO

        ! Interpolate the exact echotop pressure log-linearily and take the min to the pre-existing "old" value:
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jk, pechotop, zpA, zpB, zzdbzA, zzdbzB)
        DO jc = i_startidx, i_endidx
          IF (jk_echotop(jc) >= -900) THEN
            jk = jk_echotop(jc)
            ! Interpolation to the pressure where the reflectivity-equivalent normalized rain rate Z^(2/3)
            !  equals the equivalently transformed threshold:
            ! The reflectivity in dbz3d_lin(:,:,:) is linear:
            zzdbzA = EXP( 0.66666_wp * LOG(dbz3d_lin(jc,jk  ,jb)+repsilon) )
            zzdbzB = EXP( 0.66666_wp * LOG(dbz3d_lin(jc,jk-1,jb)+repsilon) ) ! Should be < zzdbzA, according to the logic above
            zzdbzB = MIN(zzdbzB - zzdbzA, -repsilon) ! To prevent numerical division by "almost" 0 in the next line
            ! This is the logarithmically interpolated pressure:
            zpA  = LOG( p_diag%pres(jc,jk  ,jb) )
            zpB  = LOG( p_diag%pres(jc,jk-1,jb) )
            pechotop = EXP(zpA + (zpB-zpA) / zzdbzB * (zzthresh-zzdbzA))
            IF ( echotop_p(jc,lev_etop,jb) < 0.0_wp ) THEN
              echotop_p (jc,lev_etop,jb) = pechotop
            ELSE
              echotop_p (jc,lev_etop,jb) = MIN(echotop_p(jc,lev_etop,jb), pechotop)
            END IF
          END IF
        END DO

        !$ACC END PARALLEL

      END DO
      !$ACC WAIT(1)
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    END DO

    !$ACC END DATA

  END SUBROUTINE compute_field_echotop

  !>
  !! Compute ECHOTOPs in m MSL from linear dbz3d_lin
  !!
  SUBROUTINE compute_field_echotopinm( ptr_patch, jg, p_metrics, dbz3d_lin, echotop_z, lacc )

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch        !< patch on which computation is performed
    INTEGER,              INTENT(IN)  :: jg               !< domain ID of grid
    TYPE(t_nh_metrics),   INTENT(IN)  :: p_metrics 
    REAL(wp),             INTENT(IN)  :: dbz3d_lin(:,:,:) !< reflectivity in mm^6/m^3

    REAL(wp),             INTENT(INOUT) :: echotop_z(:,:,:)  !< input/output variable, dim: (nproma,nechotop,nblks_c)

    LOGICAL,    OPTIONAL, INTENT(IN)  :: lacc             !< initialization flag

    INTEGER               :: i_rlstart,  i_rlend
    INTEGER               :: i_startblk, i_endblk
    INTEGER               :: i_startidx, i_endidx
    INTEGER               :: jb, jc, jk, lev_etop
    INTEGER               :: jk_echotop(1:nproma)
    REAL(wp)              :: zechotop
    REAL(wp)              :: zthresh, zzthresh, zA, zB, zzdbzA,  zzdbzB

    REAL(wp), PARAMETER   :: repsilon = 1.0E8_wp*TINY(1.0_wp)  ! To prevent numerical division by 0 below

    LOGICAL :: lzacc             ! OpenACC flag
    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA CREATE(jk_echotop) &
    !$ACC   PRESENT(dbz3d_lin, echotop_z, kstart_moist(jg:jg), ptr_patch, p_metrics, p_metrics%z_mc) &
    !$ACC   IF(lzacc)

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    DO lev_etop = 1, echotop_meta(jg)%nechotop

      zthresh  = 10.0_wp ** ( 0.1_wp*echotop_meta(jg)%dbzthresh(lev_etop))
      zzthresh = (zthresh+repsilon) ** 0.66666_wp  ! convert to something ~rain rate for interpolation,
                                                   ! take into account repsilon in the same way than in the EXP(LOG()) below

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,zzdbzA,zzdbzB,zA,zB,zechotop,jk_echotop), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                            i_startidx, i_endidx, i_rlstart, i_rlend)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)

        ! Find the model level just below the echotop:
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = i_startidx, i_endidx
          jk_echotop(jc)  = -999
        END DO

        !$ACC LOOP SEQ
        DO jk = MAX(kstart_moist(jg),2), ptr_patch%nlev
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = i_startidx, i_endidx
            IF ( jk_echotop(jc) < -900 .AND. dbz3d_lin(jc,jk,jb) >= zthresh ) THEN
              jk_echotop(jc) = jk
            END IF
          END DO
        END DO

        ! Interpolate the exact echotop height linearily and take the max to the pre-existing "old" value:
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jk, zA, zB, zechotop, zzdbzA, zzdbzB)
        DO jc = i_startidx, i_endidx
          IF (jk_echotop(jc) >= -900) THEN
            jk = jk_echotop(jc)
            ! Interpolation to the height where the reflectivity-equivalent normalized rain rate Z^(2/3)
            !  equals the equivalently transformed threshold:
            ! The reflectivity in dbz3d_lin(:,:,:) is linear:
            zzdbzA = EXP( 0.66666_wp * LOG(dbz3d_lin(jc,jk  ,jb)+repsilon) )
            zzdbzB = EXP( 0.66666_wp * LOG(dbz3d_lin(jc,jk-1,jb)+repsilon) ) ! Should be < zzdbzA, according to the logic above
            zzdbzB = MIN(zzdbzB - zzdbzA, -repsilon) ! To prevent numerical division by "almost" 0 in the next line
            ! This is the interpolated height:
            zA  = p_metrics%z_mc( jc, jk  , jb)  ! lower bound for linear interpolation
            zB  = p_metrics%z_mc( jc, jk-1, jb)  ! upper bound
            zechotop = zA + (zB-zA) / zzdbzB * (zzthresh-zzdbzA)
            echotop_z (jc,lev_etop,jb) = MAX(echotop_z(jc,lev_etop,jb), zechotop)
          END IF
        END DO

        !$ACC END PARALLEL

      END DO
      !$ACC WAIT(1)
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    END DO

    !$ACC END DATA

  END SUBROUTINE compute_field_echotopinm

  !>
  !! Calculate sunshine duration
  !!
  !! According to WMO (2003),2 sunshine duration during a given period is defined as
  !! the sum of that sub-period for which the perpendicular direct solar irradiance exceeds 120 W m-2
  !! WMO-No. 8 Guide to Meteorological Instruments and Methods of Observation
  !!
  !! The direct solar irradiance at the surface is calculated from the shortwave net flux at surface, 
  !! the shortwave upward flux and the shortwave diffuse downward radiative flux. It is divided
  !! by the cosine of solar zenith angle to get the perpendicular solar irradiance.
  !! If the direct solar irradiance exeeds 120 Wm-2 the sunshine duration is extended by the fast physics timestep
  !!
  !! The direct solar irradiance at the surface can be scaled by the namelist parameter tune_dursun_scaling
  !! (default is 1) to reduce the sunshine duration bias. This might be needed to account for the delta-Eddington 
  !! scaling in ecRad and other biases (e.g. ice water path)
  !!
  !! settings for sunshine duration: 
  !!
  !! dursun_thresh is the threshold for solar direct irradiance in W/m2
  !!       above which the sunshine duration is increased (default 120 W/m2)
  !!
  !! dursun_thresh_width is the smoothness / width of the threshold
  !!       function (e.g. if equal to 60 W/m2 the sunshine duration will
  !!       increase from zero at 170 W/m2 to dt at 230 W/m2)
  !!
  !! Note: MeteoSwiss uses 200 W/m2 instead of WMO value of 120 W/m2 and a width of 60 W/m2
  !!
  !! Note: use dursun_thresh=120.0_wp and dursun_thresh_width=0.01_wp to
  !!       reproduce original behaviour, according to WMO
  !!
  SUBROUTINE compute_field_dursun( pt_patch, dt_phy, dursun,                    &
    &                              swflxsfc, swflx_up_sfc, swflx_dn_sfc_diff,   &
    &                              cosmu0, dursun_thresh, dursun_thresh_width,  &
    &                              dursun_m, dursun_r, zsct, pres, twater, lacc)

    TYPE(t_patch),      INTENT(IN)    :: pt_patch              !< patch on which computation is performed
    REAL(wp),           INTENT(IN)    :: dt_phy                !< time interval for fast physics
    REAL(wp), INTENT(INOUT)           :: dursun(:,:)           !< sunshine duration (s)
    REAL(wp),           INTENT(IN)    :: swflxsfc(:,:)         !< shortwave net flux at surface [W/m2]
    REAL(wp),           INTENT(IN)    :: swflx_up_sfc(:,:)     !< shortwave upward flux at the surface [W/m2]
    REAL(wp),           INTENT(IN)    :: swflx_dn_sfc_diff(:,:)!< shortwave diffuse downward radiative flux at the surface [W/m2]
    REAL(wp),           INTENT(IN)    :: cosmu0(:,:)           !< cosine of solar zenith angle
    REAL(wp),           INTENT(IN)    :: dursun_thresh         !< threshold for solar direct irradiance in W/m2
    REAL(wp),           INTENT(IN)    :: dursun_thresh_width   !< smoothness / width of the threshold
    REAL(wp), INTENT(INOUT), POINTER  :: dursun_m(:,:)         !< maximum sunshine duration (s)
    REAL(wp), INTENT(INOUT), POINTER  :: dursun_r(:,:)         !< relative sunshine duration (s)
    REAL(wp), INTENT(IN)              :: zsct                  !< solar constant (at time of year) [W/m2]
    REAL(wp), INTENT(IN), OPTIONAL    :: pres(:,:)             !< pressure
    REAL(wp), INTENT(IN), OPTIONAL    :: twater(:,:)           !< total column water
    LOGICAL,  INTENT(IN), OPTIONAL    :: lacc                  !< initialization flag

    ! Use a minimum value to avoid div0: The exact value does not make a difference
    ! as the dursun_thresh [W/m2] will not be hit for such small values anyway.
    REAL(wp) :: cosmu0_dark = 1.e-9_wp

    ! local variables
    LOGICAL  :: l_present_dursun_m, l_present_dursun_r
    REAL(wp) :: sun_el, swrad_dir, theta_sun, xval
    INTEGER  :: i_rlstart,  i_rlend
    INTEGER  :: i_startblk, i_endblk
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: jb, jc
    LOGICAL  :: lzacc             ! OpenACC flag 
    CALL set_acc_host_or_device(lzacc, lacc)
    !$ACC DATA &
    !$ACC   PRESENT(cosmu0, dursun, dursun_m, dursun_r, pres, pt_patch) &
    !$ACC   PRESENT(swflxsfc, swflx_up_sfc, swflx_dn_sfc_diff, twater) &
    !$ACC   IF(lzacc)

    l_present_dursun_m = .FALSE.
    l_present_dursun_r = .FALSE.
    IF (ASSOCIATED(dursun_m)) l_present_dursun_m=.TRUE.
    IF (ASSOCIATED(dursun_r)) l_present_dursun_r=.TRUE.

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = pt_patch%cells%start_block( i_rlstart )
    i_endblk   = pt_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,sun_el,swrad_dir,theta_sun,xval), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( pt_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR PRIVATE(sun_el, swrad_dir, theta_sun, xval)
      DO jc = i_startidx, i_endidx
        IF(cosmu0(jc,jb)>cosmu0_dark) THEN

          ! compute sunshine duration
          ! direct solar irradiance is scaled by the tuning factor tune_dursun_scaling (default is 1) to reduce 
          ! a possible bias
          xval = (swflxsfc(jc,jb) + swflx_up_sfc(jc,jb) - swflx_dn_sfc_diff(jc,jb))*tune_dursun_scaling/cosmu0(jc,jb)
          xval = (xval - dursun_thresh)/(dursun_thresh_width/pi)
          IF (xval > 0.5_wp*pi) THEN
            dursun(jc,jb) = dursun(jc,jb) + dt_phy
          ELSEIF (xval > -0.5_wp*pi) THEN
            dursun(jc,jb) = dursun(jc,jb) + dt_phy* 0.5_wp*(SIN(xval) + 1.0_wp)
          ENDIF
        ENDIF

        IF (l_present_dursun_m) THEN
          ! estimate direct solar radiation for cloud free conditions 
          ! (after R. G. Allen et al. 2006, Agricultural and Forest Meteorology 
          !  doi:10.1016/j.agrformet.2006.05.012                               )
          ! The prefactor in eq. (17) is 0.94 instead of 0.98 in order to better
          ! fit the COSMO clear sky radiation. The turbidity factor $K_t$ is set
          ! to 0.8 (between 0.5 for extremely turbid air and 1.0 for clean air)

          ! sun elevation angle
          sun_el = ASIN(cosmu0(jc,jb))
          ! sun elevation angle in radians
          theta_sun = MAX(0.0_wp, sun_el)

          IF ( swflxsfc(jc,jb) > 0.0001_wp .AND. SIN(theta_sun) > 0.0001_wp ) THEN
            swrad_dir = zsct * 0.94_wp * EXP(                                    &
                 - 0.00146_wp * pres(jc,jb) / 1.0E3_wp / 0.8_wp / SIN(theta_sun) &
                 - 0.075_wp * (twater(jc,jb)/SIN(theta_sun))**0.4_wp             )
          ELSE
            swrad_dir = 0.0_wp
          ENDIF
  
          ! maximum possible sunshine duration (same formula as for SSD above)
          xval = (swrad_dir-dursun_thresh)/(dursun_thresh_width/pi)
          IF (xval > 0.5_wp*pi) THEN
            dursun_m(jc,jb) = dursun_m(jc,jb) + dt_phy
          ELSEIF (xval > -0.5_wp*pi) THEN
            dursun_m(jc,jb) = dursun_m(jc,jb) + dt_phy* 0.5_wp*(SIN(xval) + 1.0_wp)
          ENDIF
        ENDIF
  
        IF (l_present_dursun_r) THEN
          ! relative sunshine duration (%)
          IF (dursun_m(jc,jb) > 0.0_wp) THEN
            dursun_r(jc,jb) = 100.0_wp*dursun(jc,jb)/dursun_m(jc,jb)
          ENDIF
          dursun_r(jc,jb) = MIN(100.0_wp, MAX(0.0_wp, dursun_r(jc,jb)))
        ENDIF

      END DO ! jc loop
      !$ACC END PARALLEL

    END DO ! jb loop
    !$ACC WAIT(1)
    !$ACC END DATA
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE compute_field_dursun

  SUBROUTINE compute_updraft_duration(ptr_patch, w, wdur, wup_mask)
    TYPE(t_patch),      INTENT(IN)    :: ptr_patch              !< patch on which computation is performed
    REAL(wp), INTENT(IN) :: w(:,:,:)
    REAL(wp), INTENT(INOUT) :: wdur(:,:)
    INTEGER, INTENT(INOUT) :: wup_mask(:,:)
    REAL(wp), DIMENSION(nproma,ptr_patch%nblks_c) :: wdur_prev
    INTEGER, DIMENSION(nproma,ptr_patch%nblks_c) :: wup_mask_prev

    ! Vertex blocks of current cell
    INTEGER :: my_neighbor_blk ( 3 )

    ! Vertex indices of current cell
    INTEGER :: my_neighbor_idx ( 3 )

    INTEGER :: i,k,jb,jc,i_rlstart,i_rlend,i_startblk,i_endblk,i_startidx,i_endidx,dt

    !$ACC DATA CREATE(wdur_prev, wup_mask_prev, my_neighbor_idx, my_neighbor_blk)

    CALL init(wup_mask_prev, lacc=.TRUE.)
    CALL init(wdur_prev, lacc=.TRUE.)

    dt = dtime
    !------------------------------------------------------------------------------ 

    ! without halo or boundary points: 
    i_rlstart = grf_bdywidth_c + 1 
    i_rlend   = min_rlcell_int 
    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,           &
                          i_startidx, i_endidx, i_rlstart, i_rlend)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        wdur_prev(jc,jb) = wdur(jc,jb)
        wup_mask_prev(jc,jb) = wup_mask(jc,jb)
        wup_mask(jc,jb) = 0
        wdur(jc,jb) = 0.0_wp
        DO k = 1, ptr_patch%nlev
          IF ( w(jc,k,jb) >= 10.0_wp) THEN
            wup_mask(jc,jb) = 1
          END IF
        END DO
      END DO
      !$ACC END PARALLEL
    END DO

    CALL sync_patch_array(SYNC_C,ptr_patch, wdur_prev)
    CALL sync_patch_array(SYNC_C,ptr_patch, wup_mask_prev)

    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,           &
                          i_startidx, i_endidx, i_rlstart, i_rlend)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(my_neighbor_idx, my_neighbor_blk)
      DO jc = i_startidx, i_endidx

        ! ------------------------------------- Find neighbors ------------------------------------

        DO i = 1,3
            my_neighbor_idx(i) = ptr_patch%cells%neighbor_idx(jc, jb, i)
            my_neighbor_blk(i) = ptr_patch%cells%neighbor_blk(jc, jb, i)
        END DO

        ! ------------------------------------------------------------------------
        IF ( wup_mask(jc,jb) == 1                                        .OR. &
             wup_mask_prev(jc,jb) == 1                                   .OR. &
             wup_mask_prev(my_neighbor_idx(1), my_neighbor_blk(1)) == 1  .OR. &
             wup_mask_prev(my_neighbor_idx(2), my_neighbor_blk(2)) == 1  .OR. &
             wup_mask_prev(my_neighbor_idx(3), my_neighbor_blk(3)) == 1 ) THEN
          wdur(jc,jb) = MAX(                                    &
             wdur_prev(jc,jb),                                  &
             wdur_prev(my_neighbor_idx(1), my_neighbor_blk(1)), &
             wdur_prev(my_neighbor_idx(2), my_neighbor_blk(2)), &
             wdur_prev(my_neighbor_idx(3), my_neighbor_blk(3)) ) + REAL(dt,kind=wp)

        END IF
      END DO
      !$ACC END PARALLEL
    END DO
    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE compute_updraft_duration

  !>
  !! Calculate dhail
  !!
  SUBROUTINE compute_field_dhail(ptr_patch, p_metrics, p_prog, p_prog_rcf, p_diag,nwp_p_diag, topography_c, dhail)

    TYPE(t_patch),      INTENT(IN)    :: ptr_patch              !< patch on which computation is performed
    TYPE(t_nh_prog),    INTENT(IN)    :: p_prog                 !< prognostic variables
    TYPE(t_nh_prog),    INTENT(IN)    :: p_prog_rcf             !< prognostic variables (reduced calling frequency)
    TYPE(t_nh_diag),    INTENT(IN)    :: p_diag                 !< diagnostic variables
    TYPE(t_nwp_phy_diag),  INTENT(IN) :: nwp_p_diag             !< diagnostic variables
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics              !< metric variables
    REAL(wp),           INTENT(IN)    :: topography_c(:,:)      !< height of ground surface above sea level

    REAL(wp), INTENT(INOUT) :: dhail(:,:,:) !< expected hail diameter

#ifdef _OPENMP
    CALL warning('compute_field_dhail','No OpenMP pragmas implemented for hailcast!')
#endif

    CALL hailstone_driver(ptr_patch, p_metrics, p_prog, p_prog_rcf, p_diag,nwp_p_diag%wdur, topography_c, wdur_min_hailcast , dhail)

  END SUBROUTINE compute_field_dhail

  SUBROUTINE compute_hail_statistics(ptr_patch, p_metrics, p_prog, p_prog_rcf, p_diag,nwp_p_diag, topography_c)
    TYPE(t_patch),      INTENT(IN)    :: ptr_patch              !< patch on which computation is performed
    TYPE(t_nh_prog),    INTENT(IN)    :: p_prog                 !< prognostic variables
    TYPE(t_nh_prog),    INTENT(IN)    :: p_prog_rcf             !< prognostic variables (reduced calling frequency)
    TYPE(t_nh_diag),    INTENT(IN) :: p_diag                 !< diagnostic variables
    TYPE(t_nwp_phy_diag),    INTENT(INOUT) :: nwp_p_diag                 !< diagnostic variables
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics              !< metric variables
    REAL(wp),           INTENT(IN)    :: topography_c(:,:)      !< height of ground surface above sea level

    ! Helper variable for computation of standard deviation
    REAL(wp) :: sum_sd

    INTEGER :: k,jb,jc,i_rlstart,i_rlend,i_startblk,i_endblk,i_startidx,i_endidx


    ! without halo or boundary points: 
    i_rlstart = grf_bdywidth_c + 1 
    i_rlend   = min_rlcell_int 
    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    CALL compute_field_dhail ( ptr_patch, p_metrics, p_prog, &
                               p_prog_rcf, p_diag, nwp_p_diag, topography_c, nwp_p_diag%dhail )


    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,           &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(sum_sd)
      DO jc = i_startidx, i_endidx


        ! maximize hail size
        nwp_p_diag%dhail_mx(jc,jb) = MAX(nwp_p_diag%dhail(jc,1,jb),nwp_p_diag%dhail(jc,2,jb),nwp_p_diag%dhail(jc,3,jb),nwp_p_diag%dhail(jc,4,jb),nwp_p_diag%dhail(jc,5,jb))

        ! Compute average hail size
        nwp_p_diag%dhail_av(jc,jb) = SUM(nwp_p_diag%dhail(jc,:,jb)) / REAL(SIZE(nwp_p_diag%dhail,2),wp)

        sum_sd = 0.0_wp
        DO k = 1, SIZE(nwp_p_diag%dhail,2)
            sum_sd = sum_sd + ( nwp_p_diag%dhail(jc,k,jb) - nwp_p_diag%dhail_av(jc,jb) )**2
        END DO
        nwp_p_diag%dhail_sd(jc,jb) = SQRT( sum_sd / REAL(SIZE(nwp_p_diag%dhail, 2)-1,wp) ) 

      END DO
      !$ACC END PARALLEL

    END DO

  END SUBROUTINE compute_hail_statistics

  !! Compute vertical wind shear of either u or v component as the difference between a height AGL and the lowest model level.
  !! This is done for a number of heights and is stored in a pseudo 3D field.
  !!
  SUBROUTINE compute_field_wshear( ptr_patch, p_metrics, u_or_v, wshear_heights, wshear )

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch         !< Patch on which computation is performed
    TYPE(t_nh_metrics),   INTENT(IN)  :: p_metrics
    REAL(wp),             INTENT(IN)  :: u_or_v(:,:,:)     !< Either U or V field on model levels in m/s
    REAL(wp),             INTENT(IN)  :: wshear_heights(:) !< List of height levels for which to compute the shear

    REAL(wp),             INTENT(OUT) :: wshear(:,:,:)  !< output variable, dim: (nproma,SIZE(wshear_heights),nblks_c)

    INTEGER  :: i_rlstart,  i_rlend
    INTEGER  :: i_startblk, i_endblk
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: jb, jk, jc, lev_wshear
    REAL(wp) :: zl, zu, zint, zsurf

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    wshear(:,:,:) = 0.0_wp
    
    DO lev_wshear = 1, SIZE(wshear_heights)
      
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,zsurf,zint,zl,zu), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                            i_startidx, i_endidx, i_rlstart, i_rlend)

        ! Linear interpolate u_or_v to height zint AGL:
        DO jk = 2, ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            zsurf = p_metrics%z_ifc(jc,ptr_patch%nlev+1,jb)
            ! Height AGL for vertical interpolation with a correction if it is below the lowest level:
            zint = MAX( wshear_heights(lev_wshear), p_metrics%z_mc(jc,ptr_patch%nlev,jb) - zsurf )
            zl = p_metrics%z_mc(jc,jk  ,jb) - zsurf  ! lower bound for linear interpolation
            zu = p_metrics%z_mc(jc,jk-1,jb) - zsurf  ! upper bound
            IF ( zl <= zint .AND. zint < zu ) THEN
              wshear(jc,lev_wshear,jb) = u_or_v(jc,jk,jb) &
                   &                     + (u_or_v(jc,jk-1,jb)-u_or_v(jc,jk,jb)) / (zu-zl) * (zint-zl) &
                   &                     - u_or_v(jc,ptr_patch%nlev,jb)
            END IF
          END DO
        END DO

      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    END DO
      
  END SUBROUTINE compute_field_wshear

  !>
  !! Compute the temperature lapse rate between two pressure levels:
  !!
  !!   dTdz = T(pu) - T(pl)
  !!
  SUBROUTINE compute_field_lapserate( ptr_patch, p_metrics, p_diag, pu, pl, lapserate )

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch         !< Patch on which computation is performed
    TYPE(t_nh_metrics),   INTENT(IN)  :: p_metrics
    TYPE(t_nh_diag),      INTENT(IN)  :: p_diag            !< Diag. state which contains temp and pres
    REAL(wp),             INTENT(IN)  :: pu, pl            !< Upper and lower pressure level for computing the lapse rate

    REAL(wp),             INTENT(OUT) :: lapserate(:,:)    !< output variable, dim: (nproma,nblks_c)

    INTEGER  :: i_rlstart,  i_rlend
    INTEGER  :: i_startblk, i_endblk
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: jb, jk, jc, n_limit
    REAL(wp) :: zu(1:nproma), zl(1:nproma), tu(1:nproma), tl(1:nproma), pu_loc(1:nproma), pl_loc(1:nproma), &
                pA, pB, zpA, zpB, tpA, tpB
#ifndef _OPENACC
    LOGICAL  :: lfound_pu(1:nproma), lfound_pl(1:nproma)
#endif

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    lapserate(:,:) = 0.0_wp

    ! To avoid effects of extreme near-surface T-profiles, the pressure levels should at least
    !  be some levels above ground. We choose the third-lowest level:
    n_limit = ptr_patch%nlev - 2
    
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,pu_loc,pl_loc,&
#ifndef _OPENACC
!$OMP            lfound_pu,lfound_pl,&
#endif
!$OMP            pA,pB,zu,zl,tu,tl,zpA,zpB,tpA,tpB), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      ! limit the pressure level values by the pressure at the third-lowest level:
      DO jc = i_startidx, i_endidx
        pu_loc(jc) = MIN(pu, p_diag%pres(jc,n_limit,jb))
        pl_loc(jc) = MIN(pl, p_diag%pres(jc,n_limit,jb))
      END DO

      zu(:) = -HUGE(1.0_wp)
      zl(:) = -HUGE(1.0_wp)

#ifndef _OPENACC
      lfound_pu(:) = .FALSE.
      lfound_pl(:) = .FALSE.
#endif

      DO jk = n_limit, 2, -1

#ifndef _OPENACC
        IF (ALL( lfound_pu(i_startidx:i_endidx) .AND. lfound_pl(i_startidx:i_endidx) )) EXIT
#endif

        DO jc = i_startidx, i_endidx

          pA = p_diag%pres(jc,jk-1,jb)
          pB = p_diag%pres(jc,jk  ,jb)
          
          IF (pA < pu_loc(jc) .AND. pu_loc(jc) <= pB) THEN
            ! Logarithmically interpolated height equivalent to pressure pu:
            zpA  = p_metrics%z_mc(jc,jk-1,jb)
            zpB  = p_metrics%z_mc(jc,jk  ,jb)
            zu(jc) = zpA + (zpB-zpA) / (LOG(pB)-LOG(pA)) * (LOG(pu_loc(jc))-LOG(pA))
            ! Interpolated temperature at that height:
            tpA  = p_diag%temp(jc,jk-1,jb)
            tpB  = p_diag%temp(jc,jk  ,jb)
            tu(jc) = tpA + (tpB-tpA) / (zpB-zpA) * (zu(jc)-zpA)
#ifndef _OPENACC
            lfound_pu(jc) = .TRUE.
#endif
          END IF
          
          IF (pA < pl_loc(jc) .AND. pl_loc(jc) <= pB) THEN
            ! Logarithmically interpolated height equivalent to pressure pl:
            zpA  = p_metrics%z_mc(jc,jk-1,jb)
            zpB  = p_metrics%z_mc(jc,jk  ,jb)
            zl(jc) = zpA + (zpB-zpA) / (LOG(pB)-LOG(pA)) * (LOG(pl_loc(jc))-LOG(pA))
            ! Interpolated temperature at that height:
            tpA  = p_diag%temp(jc,jk-1,jb)
            tpB  = p_diag%temp(jc,jk  ,jb)
            tl(jc) = tpA + (tpB-tpA) / (zpB-zpA) * (zl(jc)-zpA)
#ifndef _OPENACC
            lfound_pl(jc) = .TRUE.
#endif
          END IF
          
        END DO
      END DO

      DO jc = i_startidx, i_endidx
        IF (ABS(zu(jc)-zl(jc)) >= 1e-6_wp .AND. zu(jc) > -HUGE(1.0_wp) .AND. zl(jc) > -HUGE(1.0_wp)) THEN
          lapserate(jc,jb) = (tu(jc)-tl(jc)) / (zu(jc)-zl(jc))
        END IF
      END DO
      
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_lapserate

  !>
  !! Compute the storm relative helicity with the storm motion vector estimate following Bunkers et al. (2000).
  !!  Up to now, we consider right-movers.
  !!  The height layer for computing the mean wind is taken as [0, z_up_meanwind] (m AGL).
  !!  The height layers for computing the Bunkers shear vector are taken as:
  !!    lower layer = [z_low_shear-dz_shear/2 , z_low_shear+dz_shear/2 ] (m AGL)
  !!    upper layer = [z_up_shear -dz_shear/2 , z_up_shear +dz_shear/2 ] (m AGL)
  !!  The upper limit for SRH integration is z_up_srh (m AGL).
  !!
  SUBROUTINE compute_field_srh( ptr_patch, p_metrics, p_diag, &
       z_up_srh, z_up_meanwind, z_low_shear, z_up_shear, dz_shear, &
       srh )

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch         !< Patch on which computation is performed
    TYPE(t_nh_metrics),   INTENT(IN)  :: p_metrics
    TYPE(t_nh_diag),      INTENT(IN)  :: p_diag            !< Diag. state which contains temp and pres
    REAL(wp),             INTENT(IN)  :: z_up_srh(:)       !< Vector of upper height levels AGL for vertical integration of SRH
    REAL(wp),             INTENT(IN)  :: z_up_meanwind     !< Upper height AGL for mean wind vector computation
    REAL(wp),             INTENT(IN)  :: z_low_shear       !< Lower height AGL for shear vector computation
    REAL(wp),             INTENT(IN)  :: z_up_shear        !< Upper height AGL for shear vector computation
    REAL(wp),             INTENT(IN)  :: dz_shear          !< level thickness for the lower and upper layers in shear vector computation

    REAL(wp),             INTENT(OUT) :: srh(:,:,:)          !< output variable, dim: (nproma,SIZE(z_up_srh),nblks_c)

    INTEGER  :: i_rlstart,  i_rlend
    INTEGER  :: i_startblk, i_endblk
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: lev_srh, jb, jk, jc, jku, jkl, k_start, k_start_vec(nproma)
    REAL(wp) :: speed_shear, r_or_left_fac, du, dv, dz, max_height
    REAL(wp), DIMENSION(nproma) :: u_mean, v_mean, u_shear, v_shear, u_storm, v_storm, &
         u_shear_up, u_shear_low, v_shear_up, v_shear_low
    REAL(wp), DIMENSION(nproma,ptr_patch%nlev) :: srh_layer

    REAL(wp), PARAMETER :: ref_speed_bunkers = 7.5_wp  ! m/s

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    max_height = MAXVAL([MAXVAL(z_up_srh(:)), z_up_meanwind, z_up_shear+dz_shear*0.5_wp, z_low_shear+dz_shear*0.5_wp]) ! m AGL

!$OMP PARALLEL
    CALL init(srh(:,:,:), 0.0_wp)
!$OMP DO PRIVATE(jb,jc,lev_srh,i_startidx,i_endidx,k_start,k_start_vec, &
!$OMP            speed_shear,u_mean,v_mean,u_shear,v_shear,u_storm,v_storm, &
!$OMP            u_shear_up,u_shear_low,v_shear_up,v_shear_low,r_or_left_fac, &
!$OMP            srh_layer,du,dv,jku,jkl,dz), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
           i_startidx, i_endidx, i_rlstart, i_rlend)

      ! Determine the lowermost height level (highest index) which is
      !  everywhere just above the maximum height considered
      !  for the following vertical averages and integrals,
      !  and use this as starting level for all vertical loops to save computing time.
      ! We do a bottom-up search for the level just below z_up and subtract 1 to get one level above:
      k_start_vec(:) = ptr_patch%nlev
      DO jk = ptr_patch%nlev, 2, -1
        DO jc = i_startidx, i_endidx
          IF (p_metrics%z_ifc(jc,jk,jb) - p_metrics%z_ifc(jc,ptr_patch%nlev+1,jb) < max_height) THEN
            k_start_vec(jc) = jk - 1 
          END IF
        END DO
      END DO
      ! k_start_vec contains the max height index in each column. The overall k_start index
      ! is the smallest among them:
      k_start = MINVAL(k_start_vec(i_startidx:i_endidx))

      ! mean U-component of 0-z_up_meanwind AGL:
      CALL vert_integral_vec_1d ( i_startidx, i_endidx, k_start,    &
           &                      hhl  = p_metrics % z_ifc(:,:,jb), &
           &                      f    = p_diag    % u(:,:,jb),     &
           &                      zlow = 0.0_wp,                    &
           &                      zup  = z_up_meanwind,             &
           &                      fint = u_mean(:),                 &
           &                      l_agl= .TRUE.,                    &
           &                      l_calc_mean = .TRUE.,             &
           &                      l_rescale_to_full_thickness = .FALSE. &
           &                      )

      ! mean V-component of 0-z_up_meanwind AGL:
      CALL vert_integral_vec_1d ( i_startidx, i_endidx, k_start,    &
           &                      hhl  = p_metrics % z_ifc(:,:,jb), &
           &                      f    = p_diag    % v(:,:,jb),     &
           &                      zlow = 0.0_wp,                    &
           &                      zup  = z_up_meanwind,             &
           &                      fint = v_mean(:),                 &
           &                      l_agl= .TRUE.,                    &
           &                      l_calc_mean = .TRUE.,             &
           &                      l_rescale_to_full_thickness = .FALSE. &
           &                      )

      ! Upper U-component for shear vector, averaged over the layer
      !  [ z_up_shear-dz_shear/2 ; z_up_shear+dz_shear/2] AGL:
      CALL vert_integral_vec_1d ( i_startidx, i_endidx, k_start,      &
           &                      hhl  = p_metrics % z_ifc(:,:,jb),   &
           &                      f    = p_diag    % u(:,:,jb),       &
           &                      zlow = z_up_shear-0.5_wp*dz_shear,  &
           &                      zup  = z_up_shear+0.5_wp*dz_shear,  &
           &                      fint = u_shear_up(:),               &
           &                      l_agl= .TRUE.,                      &
           &                      l_calc_mean = .TRUE.,               &
           &                      l_rescale_to_full_thickness = .FALSE. &
           &                      )

      ! Upper V-component for shear vector, averaged over the layer
      !  [ z_up_shear-dz_shear/2 ; z_up_shear+dz_shear/2] AGL:
      CALL vert_integral_vec_1d ( i_startidx, i_endidx, k_start,      &
           &                      hhl  = p_metrics % z_ifc(:,:,jb),   &
           &                      f    = p_diag    % v(:,:,jb),       &
           &                      zlow = z_up_shear-0.5_wp*dz_shear,  &
           &                      zup  = z_up_shear+0.5_wp*dz_shear,  &
           &                      fint = v_shear_up(:),               &
           &                      l_agl= .TRUE.,                      &
           &                      l_calc_mean = .TRUE.,               &
           &                      l_rescale_to_full_thickness = .FALSE. &
           &                      )

      ! Lower U-component for shear vector, averaged over the layer
      !  [ z_low_shear-dz_shear/2 ; z_low_shear+dz_shear/2] AGL:
      CALL vert_integral_vec_1d ( i_startidx, i_endidx, k_start,      &
           &                      hhl  = p_metrics % z_ifc(:,:,jb),   &
           &                      f    = p_diag    % u(:,:,jb),       &
           &                      zlow = z_low_shear-0.5_wp*dz_shear, &
           &                      zup  = z_low_shear+0.5_wp*dz_shear, &
           &                      fint = u_shear_low(:),              &
           &                      l_agl= .TRUE.,                      &
           &                      l_calc_mean = .TRUE.,               &
           &                      l_rescale_to_full_thickness = .FALSE. &
           &                      )

      ! Lower V-component for shear vector, averaged over the layer
      !  [ z_low_shear-dz_shear/2 ; z_low_shear+dz_shear/2] AGL:
      CALL vert_integral_vec_1d ( i_startidx, i_endidx, k_start,      &
           &                      hhl  = p_metrics % z_ifc(:,:,jb),   &
           &                      f    = p_diag    % v(:,:,jb),       &
           &                      zlow = z_low_shear-0.5_wp*dz_shear, &
           &                      zup  = z_low_shear+0.5_wp*dz_shear, &
           &                      fint = v_shear_low(:),              &
           &                      l_agl= .TRUE.,                      &
           &                      l_calc_mean = .TRUE.,               &
           &                      l_rescale_to_full_thickness = .FALSE. &
           &                      )

      ! Storm motion vector following Bunkers et al. (2000):
      r_or_left_fac = 1.0_wp   ! 1.0 for rightmovers, -1.0 for leftmovers
      DO jc = i_startidx, i_endidx
        u_shear(jc) = u_shear_up(jc) - u_shear_low(jc)
        v_shear(jc) = v_shear_up(jc) - v_shear_low(jc)
        speed_shear = SQRT(u_shear(jc)**2 + v_shear(jc)**2)
        u_storm(jc) = u_mean(jc) + r_or_left_fac * ref_speed_bunkers * v_shear(jc) / speed_shear
        v_storm(jc) = v_mean(jc) - r_or_left_fac * ref_speed_bunkers * u_shear(jc) / speed_shear
      END DO

      ! SRH-contribution in each model layer:
      DO jk = k_start, ptr_patch%nlev
        jku = MAX(jk-1,1)
        jkl = MIN(jk+1,ptr_patch%nlev)
        DO jc = i_startidx, i_endidx
          du = p_diag%u(jc,jku,jb) - p_diag%u(jc,jkl,jb)
          dv = p_diag%v(jc,jku,jb) - p_diag%v(jc,jkl,jb)
          dz = p_metrics%z_mc(jc,jku,jb) - p_metrics%z_mc(jc,jkl,jb)
          srh_layer(jc,jk) = (p_diag%v(jc,jk,jb)-v_storm(jc)) * du - (p_diag%u(jc,jk,jb)-u_storm(jc)) * dv
          srh_layer(jc,jk) = srh_layer(jc,jk) / dz
        END DO
      END DO

      ! SRH integrated from the ground up to height z_up_srh:
      DO lev_srh = 1, SIZE(z_up_srh)

        CALL vert_integral_vec_1d ( i_startidx, i_endidx, k_start,      &
             &                      hhl  = p_metrics % z_ifc(:,:,jb),   &
             &                      f    = srh_layer(:,:),              &
             &                      zlow = 0.0_wp,                      &
             &                      zup  = z_up_srh(lev_srh),           &
             &                      fint = srh(:,lev_srh,jb),           &
             &                      l_agl= .TRUE.,                      &
             &                      l_calc_mean = .FALSE.,              &
             &                      l_rescale_to_full_thickness = .FALSE. &
             &                      )

      END DO
      
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_srh


  !!>
  !! Vertical integration or averaging of a quantity f over a height interval [zlow,zup].
  !! If l_agl=.true., the heights zlow and zup are with respect to AGL. Otherwise MSL.
  !! If l_calc_mean=.true., the vertical integral is divided by the layer thickness.
  !! Otherwise, if l_rescale_to_full_thickness=.true., the vertical integral is rescaled to the
  !!  full layer thickness in case one or both layers are below the surface or above the model top.
  !!  Re-scaling means that voids are filled with the average of the present heights in the integral.
  !!
  !! NOTE: if kstart is > 1, it must be small enough so that hhl(kstart) is above zup everywhere!
  
  SUBROUTINE vert_integral_vec_1d (istart, iend, kstart, hhl, f, zlow, zup, fint, &
                                   l_agl, l_calc_mean, l_rescale_to_full_thickness)

    INTEGER,  INTENT(in)  :: istart, iend !< start- and end-indices for computations along the horizontal (nproma) dimension
    INTEGER,  INTENT(in)  :: kstart       !< start index in the vertical, may be used to shorten the vertical loops if appropriate, to save time
    REAL(wp), INTENT(in)  :: hhl(:,:)     !< one vertical slice of half level heights,            dim: (nproma,nlev+1)
    REAL(wp), INTENT(in)  :: f  (:,:)     !< one vertical slice of data to integrate over height, dim: (nproma,nlev)
    REAL(wp), INTENT(in)  :: zlow, zup    !< lower and upper height bounds for vertical integration (either MSL or AGL, depending on switch l_agl
    LOGICAL,  INTENT(in)  :: l_agl        !< if .TRUE., zlow and zup are taken as height AGL. Otherwise, MSL.
    LOGICAL,  INTENT(in)  :: l_calc_mean  !< if .TRUE., compute the vertical MEAN, not the vertical INTEGRAL
    LOGICAL,  INTENT(in)  :: l_rescale_to_full_thickness !< if .TRUE. and if the lower bound zlow is below the surface,
                                                         !  fill integral of part below surface WITH the mean VALUE from above
    
    REAL(wp), INTENT(out) :: fint(:)      !> integral of f over height, dim: (nproma)

    INTEGER  :: i, k, nlev
    REAL(wp) :: h_offset(SIZE(f, DIM=1)) ! offset for height to discriminate AGL and MSL
    REAL(wp) :: dz_layer(SIZE(f, DIM=1)) ! total layer thickness for integration
    REAL(wp) :: dz_loc                   ! contribution of the actual layer to dz_layer

    nlev = SIZE(f, DIM=2)

    IF (l_agl) THEN
      h_offset(:) = hhl(:,nlev+1)
    ELSE
      h_offset(:) = 0.0_wp
    END IF
    
    fint    (:) = 0.0_wp
    dz_layer(:) = 0.0_wp
    DO k = kstart, nlev
      DO i = istart, iend
    
        ! Parts of the grid box are within the bounds, integrate over the exact bounds [zlow,zup]:
        !  (It also works if the integration layer is so narrow that the bounds are in the same grid box)
        IF ( ( hhl(i,k+1)-h_offset(i) <= zup ) .AND. ( hhl(i,k)-h_offset(i)   >= zlow ) ) THEN

          dz_loc = MIN(hhl(i,k)-h_offset(i), zup) - MAX(hhl(i,k+1)-h_offset(i), zlow)
          dz_layer(i) = dz_layer(i) + dz_loc
          
          ! a simple box-integration in the vertical, but honouring the exact integration bounds zlow, zup;
          fint(i) = fint(i) + f(i,k) * dz_loc
            
        END IF
      END DO
    END DO

    IF (l_calc_mean) THEN
      ! Compute vertical average:
      DO i = istart, iend
        IF (dz_layer(i) > 1e-20_wp) THEN
          fint(i) = fint(i) / dz_layer(i)
        END IF
      END DO
    ELSE IF (l_rescale_to_full_thickness) THEN
      ! Fill the voids below the surface and/or above the top in the integral with it's average value:
      DO i = istart, iend
        IF (dz_layer(i) > 1e-20_wp) THEN
          fint(i) = fint(i) / dz_layer(i) * (zup-zlow)
        END IF
      END DO
    END IF

  END SUBROUTINE vert_integral_vec_1d


  !>
  !! compute near surface visibility
  !!
  !! Description:
  !!  This routine computes horizontal visibility [km] at the
  !!   surface or lowest model layer from qv, qc, qr, qi, qs, and qg.
  !! 
  !!------------------------------------------------------------------------------
  !!
  !! SUBPROGRAM:    CALVIS      CALCULATE HORIZONTAL VISIBILITY
  !!
  !!   PRGMMR:  BENJAMIN, STAN ORG: NOAA/FSL       DATE: 99-09-07
  !!
  !! ABSTRACT:
  !!
  !!   Started with Stoelinga-Warner algorithm for hydrometeors only.
  !!    Added coefficients for graupel.
  !!    Added algorithm for clear-air RH-based visibility.
  !!
  !!
  !! HISTORY
  !! PROGRAM HISTORY LOG:
  !!    99-05-                        Version from Eta model and from
  !!                                    Mark Stoelinga and Tom Warner
  !!    99-09-07      S. Benjamin     Modified for MM5 microphysics variables
  !!                                    include graupel mixing ratio
  !!    99-09         S. Benjamin     Added algorithm for RH-based clear-air
  !!                                    visibility
  !!    00-08         S. Benjamin     Added mods for base of 60km instead of 90km,
  !!                                    max RH from lowest 2 levels instead of
  !!                                    lev 2 only, max hydrometeor mix ratio
  !!                                    from lowest 5 levs instead of lev 1 only
  !!                                  Based on Schwartz stats and Smirnova et al
  !!                                    paper, and on METAR verif started this week
  !!    Dec 03        S. Benjamin     - updates
  !!                              - day/night distinction for vis constants
  !!                                  from Roy Rasmussen
  !!                              - low-level wind shear term
  !!                                  - recommended by Evan Kuchera
  !!
  !!   Mar 22        T. Goecke     - RH dependence now after Gultepe etal (2009)
  !!  !!
  !!

  SUBROUTINE compute_field_visibility(ptr_patch, p_prog, p_diag, prm_diag, jg, vis_out, lacc)

    TYPE(t_patch),        INTENT(IN)            :: ptr_patch
    TYPE(t_nh_prog),      INTENT(IN)            :: p_prog       ! nonhydrostatic state
    TYPE(t_nh_diag),      INTENT(IN)            :: p_diag
    TYPE(t_nwp_phy_diag), INTENT(INOUT)         :: prm_diag     ! physics variables
    INTEGER,              INTENT(IN)            :: jg           ! domain ID of the grid
    REAL(wp),             INTENT(OUT)           :: vis_out(:,:) ! output variable
    LOGICAL,              INTENT(IN), OPTIONAL  :: lacc         ! if true, use openacc

    !local variables
    REAL(wp), POINTER ::   qv(:,:,:) ! specific humidity (subgrid)
    REAL(wp), POINTER ::   qc(:,:,:) ! cloud water (subgrid)
    REAL(wp), POINTER ::   qi(:,:,:) ! cloud ice (subgrid)
    REAL(wp), POINTER ::   qr(:,:,:) ! cloud water (subgrid)
    REAL(wp), POINTER ::   qs(:,:,:) ! cloud water (subgrid)
    REAL(wp), POINTER ::   qg(:,:,:) ! cloud water (subgrid)
    REAL(wp), POINTER ::  rho(:,:,:) ! total density (inkluding hydrometeors)

    ! specific tracer concentrations
    REAL(wp) :: Ccmax, Cimax, Crmax, Csmax, Cgmax

    ! parameters for vis parametrization
    REAL(wp), PARAMETER :: a_c = 144.7_wp, b_c = 0.88_wp
    REAL(wp), PARAMETER :: a_i = 327.8_wp, b_i = 1.0_wp
    REAL(wp), PARAMETER :: a_r = 2.24_wp, b_r = 0.75_wp
    REAL(wp), PARAMETER :: a_s_wet = 6.0_wp, a_s_dry = 10.0_wp, b_s = 1.0_wp
    REAL(wp), PARAMETER :: a_g = 4.0_wp, b_g = 0.75_wp
    REAL(wp)            :: a_s, temp_fac, beta
    
    INTEGER, PARAMETER :: top_lev = 3  ! number of levels from ground used for vis diagnostic

    REAL(wp) :: vis, vis_night, visrh, qrh
    REAL(wp) :: pvsat, pv, rhmax, visrh_clip, &
	     &  shear, shear_fac, czen, zen_fac
    REAL(wp), ALLOCATABLE :: rh(:,:)

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc, jk, jb, nlev

    LOGICAL :: lzacc ! non-optional version of lacc

    CALL set_acc_host_or_device(lzacc, lacc)

    ! local pointers
    rho => p_prog%rho
    qv  => prm_diag%tot_ptr(iqv)%p_3d
    qc  => prm_diag%tot_ptr(iqc)%p_3d
    qi  => prm_diag%tot_ptr(iqi)%p_3d
    qr  => p_prog%tracer_ptr(iqr)%p_3d
    qs  => p_prog%tracer_ptr(iqs)%p_3d
    IF(atm_phy_nwp_config(jg)%lhave_graupel) qg => p_prog%tracer_ptr(iqg)%p_3d

    ! some parameters
    nlev = size(qv,2)
    allocate(rh(nproma,nlev-top_lev+1:nlev))
    visrh_clip = 1.0_wp ! clip RH-VIS at this km

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    !$ACC DATA PRESENT(qc, qg, qi, qr, qs, qv, rho, vis_out) &
    !$ACC   CREATE(rh) IF(lzacc)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,pvsat,pv,rh,rhmax,visrh,qrh,vis,vis_night,shear,shear_fac,&
!$OMP     Ccmax,Cimax,Crmax,Csmax,Cgmax,temp_fac,a_s,beta,czen,zen_fac), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!  I) VIS due to relative humidity   !!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(pvsat, pv)
      DO jk = nlev-top_lev+1, nlev
        DO jc = i_startidx, i_endidx
          ! compute RH in lowest 'top_level' levels
          pvsat     = esat_water(p_diag%temp(jc,jk,jb))
          pv        = p_diag%pres(jc,jk,jb)*qv(jc,jk,jb)/(rdv + o_m_rdv * qv(jc,jk,jb))
          rh(jc,jk) = 100.0_wp * MIN(1.0_wp, pv/pvsat)
        END DO
      END DO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR PRIVATE(a_s, beta) &
      !$ACC   PRIVATE(Ccmax, Cgmax, Cimax, Crmax, Csmax) &
      !$ACC   PRIVATE(czen, rhmax, shear, shear_fac, temp_fac) &
      !$ACC   PRIVATE(vis, visrh, vis_night, zen_fac)
      DO jc = i_startidx, i_endidx
        ! maximum lower two levels
        rhmax = MAX(rh(jc,nlev), rh(jc,nlev-1))

        ! vis due to haze parametrized as function of rh only, form found via fit to
        ! SYNOP station data over Germany from 11/2021
        IF (rhmax <= 40.0_wp) THEN
          visrh = 88950.37269485_wp - 327.73380915_wp*rhmax
        ELSE IF (rhmax <= 98.2_wp) THEN
          visrh = 2.74158753e-04_wp*rhmax**5 - 8.04508715e-02_wp*rhmax**4 &
            &   + 9.4148139_wp*rhmax**3 - 5.78127237e+02_wp*rhmax**2      &
            &   + 1.82682914e+04_wp*rhmax - 1.54588988e+05_wp
        ELSE
          visrh = -1171.04931497_wp*rhmax + 117111.72541055_wp
        END IF
        visrh = visrh/1000.0_wp ! convert to units of km

        ! clip below X km
        visrh = MAX(visrh, visrh_clip)

        ! add term to increase RH vis term for
        ! low-level wind shear increasing from 4 to 6 ms-1
        ! (using Evan Kuchera's paper as a guideline)
        ! calculate term for shear in the lowest about 15 mb
        ! 15mb about 120 meters, so try lev nlev-2
        shear     = SQRT(  (p_diag%u(jc,nlev-2,jb) - p_diag%u(jc,nlev,jb))**2 &
          &       + (p_diag%v(jc,nlev-2,jb) - p_diag%v(jc,nlev,jb))**2  )
        shear_fac = MIN( 1.0_wp, MAX(0.0_wp, (shear-4.0_wp)/2.0_wp) )
        IF (visrh < 10.0_wp) visrh = visrh + (10.0_wp-visrh) * shear_fac

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!  II) VIS due to hydrometeors       !!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!   The routine uses the following
        !!   expressions for extinction coefficient, beta (in km**-1),
        !!   with C being the mass concentration (in g/m**3):
        !!
        !!      cloud water:  beta = 144.7 * C ** (0.8800)
        !!      rain water:   beta =  2.24 * C ** (0.7500)
        !!      cloud ice:    beta = 327.8 * C ** (1.0000)
        !!      snow:         beta = 10.36 * C ** (0.7776)
        !!      graupel:      beta =  8.0  * C ** (0.7500)
        !!
        !!   These expressions were obtained from the following sources:
        !!
        !!      for cloud water: from Kunkel (1984)
        !!      for rainwater: from M-P dist'n, with N0=8e6 m**-4 and
        !!         rho_w=1000 kg/m**3
        !!      for cloud ice: assume randomly oriented plates which follow
        !!         mass-diameter relationship from Rutledge and Hobbs (1983)
        !!      for snow: from Stallabrass (1985), assuming beta = -ln(.02)/vis
        !!      for graupel: guestimate by John Brown and Stan Benjamin,
        !!         similar to snow, but a smaller extinction coef seemed
        !!         reasonable.  27 Aug 99
        !!
        !!   The extinction coefficient for each water species present is
        !!   calculated, and then all applicable betas are summed to yield
        !!   a single beta. Then the following relationship is used to
        !!   determine visibility (in km), where epsilon is the threshhold
        !!   of contrast, usually taken to be .02:
        !!
        !!      vis = -ln(epsilon)/beta      [found in Kunkel (1984)]
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! conversion to volumetric conentration.
        ! rho = V/m_tot (V is given by the grid anyway
        ! specific quantities q_k = m_k/m_tot -> C_k = q_k * rho  	
        ! maximize hydrometeors over lowest 'top_lev' levels, final unit= g/m^3
        Ccmax = MAXVAL(qc(jc,nlev-top_lev+1:nlev,jb)*rho(jc,nlev-top_lev+1:nlev,jb))*1000.0_wp
        Cimax = MAXVAL(qi(jc,nlev-top_lev+1:nlev,jb)*rho(jc,nlev-top_lev+1:nlev,jb))*1000.0_wp
        Crmax = MAXVAL(qr(jc,nlev-top_lev+1:nlev,jb)*rho(jc,nlev-top_lev+1:nlev,jb))*1000.0_wp
        Csmax = MAXVAL(qs(jc,nlev-top_lev+1:nlev,jb)*rho(jc,nlev-top_lev+1:nlev,jb))*1000.0_wp
        IF (atm_phy_nwp_config(jg)%lhave_graupel) THEN
          Cgmax = MAXVAL(qg(jc,nlev-top_lev+1:nlev,jb)*rho(jc,nlev-top_lev+1:nlev,jb))*1000.0_wp
        ELSE
          Cgmax = 0.0_wp
        END IF

        ! snow coefficient temperature dependent
        temp_fac  = MIN( 1.0_wp, MAX((p_diag%temp(jc,nlev,jb)-271.15_wp), 0.0_wp) )
        a_s       = a_s_dry * (1.0_wp-temp_fac) + a_s_wet * temp_fac

        ! calculate extinction coefficient  
        beta = a_c * Ccmax**b_c & ! cloud water	
          &  + a_i * Cimax**b_i & ! cloud ice
          &  + a_r * Crmax**b_r & ! rain
          &  + a_s * Csmax**b_s & ! snow
          &  + a_g * Cgmax**b_g & ! graupel   
          &  + 1.0e-10_wp         ! small offsett to prevent zero division

        ! vis after koschmieder formula with 2 percent of initial beam intensity
        vis = MIN(90.0_wp, -log(0.02)/beta)

        ! zenith angle
        czen = prm_diag%cosmu0(jc,jb)

        ! Dec 2003 - Roy Rasmussen (NCAR) expression for night vs. day vis
        ! 1.609 factor is number of km in mile.
        vis_night = 1.69_wp * ( (vis/1.609_wp)**0.86_wp ) * 1.609_wp
        zen_fac   = MIN( 0.1_wp, MAX(czen, 0.0_wp) ) / 0.1_wp
        vis       = zen_fac * vis + (1.0_wp-zen_fac) * vis_night

        ! take minumum from vis and visrh
        vis = MIN(vis, visrh)

        ! convert to meter
        vis = vis * 1000.0_wp

        ! write to diagnostic field
        vis_out(jc,jb) = vis
      END DO ! jc
      !$ACC END PARALLEL

    END DO ! jb
!$OMP END DO NOWAIT    
!$OMP END PARALLEL

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE compute_field_visibility

  !! Find the lowest inversion and provide its inversion height and lowest point of the entrainment zone 
  !! It follows Van Wevweberg et al. Month Weath. Rev. 2021
  !! This function just produces the variable for output
  !! The calculations are done in compute_field_inversion_height

  SUBROUTINE compute_field_inversion_height(ptr_patch,jg,p_metrics,p_prog,p_diag,prm_diag,inv_height)

    TYPE(t_patch),        INTENT(IN)    :: ptr_patch
    INTEGER,              INTENT(IN)    :: jg       !< domain ID of grid
    TYPE(t_nh_metrics),   INTENT(IN)    :: p_metrics
    TYPE(t_nh_prog),      INTENT(IN)    :: p_prog   ! nonhydrostatic state
    TYPE(t_nh_diag),      INTENT(IN)    :: p_diag   ! diagnostic variables  
    TYPE(t_nwp_phy_diag), INTENT(INOUT)  :: prm_diag ! physics variables
    
    REAL(WP), INTENT(OUT)   :: inv_height(:,:) ! output variable


    ! Parameters
    REAL(wp), PARAMETER   ::   no_inversion_value = -99.99_wp ! Output value when no inversion is found

    
    ! Local variables
    INTEGER  ::   i_inversion(nproma)  ! k-idex for inversion
    INTEGER  ::   i_ent_zone(nproma)   ! k-idex for entrainment zone
    LOGICAL  ::   lfound_inversion(nproma)  ! To stop loop when inversion is found



    ! Pointers
    REAL(wp),POINTER      ::   z(:,:,:)        ! height at model levels
    REAL(wp),POINTER      ::   z_ifc(:,:,:)    ! height at interface levels
    REAL(wp),POINTER      ::   te(:,:,:)       ! temperature
    REAL(wp),POINTER      ::   qc(:,:,:)       ! cloud water
    REAL(wp),POINTER      ::   prs(:,:,:)      ! pressure (from physics)
    REAL(WP),POINTER      ::   low_ent_zone(:,:) ! output variable 

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jktop,jkbot,nlev
    
    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    ! Set pointers
    prs => p_diag%pres
    te  => p_diag%temp
    qc  => p_prog%tracer_ptr(iqc)%p_3d
    z   => p_metrics%z_mc
    z_ifc => p_metrics%z_ifc
    low_ent_zone => prm_diag%low_ent_zone

    
    ! Integration limits
    jktop = kstart_moist(jg)
    jkbot = ptr_patch%nlev
    nlev  = ptr_patch%nlev


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,lfound_inversion,i_inversion, &
!$OMP            i_ent_zone), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk
      
        CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                            i_startidx, i_endidx, i_rlstart, i_rlend)

        CALL inversion_height_index(z(:,:,jb),z_ifc(:,nlev+1,jb),qc(:,:,jb),te(:,:,jb),prs(:,:,jb), &
                      &            i_startidx,i_endidx,jktop,jkbot,nlev, &
                      &            i_inversion(:),i_ent_zone(:),lfound_inversion(:))

! Calculate the inversion height in meters and set non-values
        DO jc = i_startidx, i_endidx
          IF (lfound_inversion(jc)) THEN
            inv_height(jc,jb)   = z(jc,i_inversion(jc),jb)
            low_ent_zone(jc,jb) = z(jc,i_ent_zone(jc),jb)
          ELSE
            inv_height(jc,jb)   = no_inversion_value
            low_ent_zone(jc,jb) = no_inversion_value  
          END IF
        END DO

      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE compute_field_inversion_height
END MODULE mo_opt_nwp_diagnostics

