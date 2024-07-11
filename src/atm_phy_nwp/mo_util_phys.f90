!
! Implementation of physics utility routines.
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
#include "omp_definitions.inc"
!----------------------------

MODULE mo_util_phys

  USE mo_kind,                  ONLY: wp
  USE mo_parallel_config,       ONLY: nproma
  USE mo_physical_constants,    ONLY: o_m_rdv        , & !! 1 - r_d/r_v &
    &                                 rdv,             & !! r_d / r_v
    &                                 cpd, p0ref, rd,  &
    &                                 vtmpc1, t3, grav,&
    &                                 alv,alvdcp, rd_o_cpd
  USE mo_exception,             ONLY: finish
  USE mo_satad,                 ONLY: sat_pres_water, sat_pres_ice
  USE mo_fortran_tools,         ONLY: assign_if_present, set_acc_host_or_device, assert_acc_host_only
  USE mo_impl_constants,        ONLY: min_rlcell_int
  USE mo_model_domain,          ONLY: t_patch
  USE mo_nonhydro_types,        ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,         ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_run_config,            ONLY: iqv, iqc, iqi, iqr, iqs, iqni, ininact, &
       &                              iqm_max, nqtendphy, lart, iqnc, iqnr, iqns, &
       &                              iqb_i, iqb_e
#ifndef __NO_ICON_LES__
  USE mo_ls_forcing_nml,        ONLY: is_ls_forcing, is_nudging_tq, &
       &                              nudge_start_height, nudge_full_height, dt_relax
#endif
  USE mo_loopindices,           ONLY: get_indices_c
  USE mo_atm_phy_nwp_config,    ONLY: atm_phy_nwp_config
  USE mo_nwp_tuning_config,     ONLY: tune_gust_factor, itune_gust_diag, tune_gustsso_lim
  USE mo_art_config,            ONLY: art_config
  USE mo_nonhydrostatic_config, ONLY: kstart_moist
  USE mo_2mom_mcrph_util,       ONLY: set_qnc, set_qnr, set_qni, set_qns
  USE microphysics_1mom_schemes,ONLY: get_mean_crystal_mass

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: nwp_dyn_gust
  PUBLIC :: nwp_con_gust
  PUBLIC :: virtual_temp
  PUBLIC :: vap_pres
  PUBLIC :: swdir_s
  PUBLIC :: rel_hum
  PUBLIC :: compute_field_rel_hum_wmo
  PUBLIC :: compute_field_rel_hum_ifs
  PUBLIC :: tracer_add_phytend
  PUBLIC :: calc_ustar
  PUBLIC :: inversion_height_index
  
  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_util_phys'

CONTAINS


  !-------------------------------------------------------------------------
  !!
  !! Calculate dynamic gusts as in near_surface of COSMO code
  !!     gust = ff10m + gust_factor * ustar
  !! where ff10m is the 10 m wind and the friction velocity ustar = SQRT(tcm)*ff1
  !!
  ELEMENTAL FUNCTION nwp_dyn_gust( u_10m, v_10m, tcm, u1, v1, u_env, v_env, fr_oce, mtnmask) RESULT( vgust_dyn)
    !$ACC ROUTINE SEQ

    REAL(wp), INTENT(IN) :: u_10m, &    ! zonal wind component at 10 m above ground [m/s]
      &                     v_10m, &    ! meridional wind component at 10 m above ground [m/s]
      &                     tcm  , &    ! transfer coefficient for momentum at surface
      &                     u1   , &    ! zonal wind at lowest model layer above ground [m/s]
      &                     v1   , &    ! meridional wind at lowest model layer above ground [m/s]
      &                     u_env, &    ! zonal wind at top of SSO envelope layer [m/s]
      &                     v_env, &    ! meridional wind at top of SSO envelope layer [m/s]
      &                     fr_oce,&    ! ocean fraction
      &                     mtnmask     ! mask field for weighting SSO enhancement

    REAL(wp) :: vgust_dyn               ! dynamic gust at 10 m above ground [m/s]

    REAL(wp) :: ff10m, ustar, uadd_sso, gust_nonlin, offset, base_gust, mtn_lim, oce_shift

    uadd_sso = MAX(0._wp, SQRT(u_env**2 + v_env**2) - SQRT(u1**2 + v1**2))
    SELECT CASE (itune_gust_diag)
    CASE (4)     ! gust param based on 10-min averaged wind
      offset = 6._wp+6._wp*fr_oce
    CASE (3)     ! ICON-D2 with subgrid-scale condensation
      offset = 10._wp+6._wp*fr_oce
    CASE (2)     ! ICON global with MERIT/REMA orography data
      offset = MAX(0._wp,6._wp-uadd_sso)
    CASE default ! actually (1), but code does not vectorize without default branch
      offset = 10._wp
    END SELECT
    oce_shift = MERGE(fr_oce,0._wp,itune_gust_diag>=3)

    ff10m = SQRT( u_10m**2 + v_10m**2)
    ustar = calc_ustar(tcm, u1, v1)
    gust_nonlin = MAX(0._wp,MIN(2._wp,0.2_wp*(ff10m-offset)))
    base_gust = ff10m + (tune_gust_factor-oce_shift+gust_nonlin)*ustar
    mtn_lim = MAX(0._wp, MIN(1._wp, 2._wp - base_gust/tune_gustsso_lim) )
    vgust_dyn = base_gust + mtnmask*mtn_lim*(uadd_sso + (2._wp+gust_nonlin)*ustar )

  END FUNCTION nwp_dyn_gust

  !-------------------------------------------------------------------------
  !!
  !! Calculate friction velocity ustar = SQRT(tcm)*ff1.
  !! Taken from the original implementation by H. Frank from function 
  !! nwp_dyn_gust to be also usable for other purposes.
  !!
  ELEMENTAL FUNCTION calc_ustar(tcm, u1, v1) RESULT (ustar)
    !$ACC ROUTINE SEQ

    REAL(wp), INTENT(in)  :: &
      &  tcm,                & !< Transfer coefficient for momentum at surface
      &  u1, v1                !< Horizontal wind components at lowest model layer (m/s)
    REAL(wp)              :: &
      &  ustar                 !< Friction velocity

    ustar = SQRT( MAX( tcm, 5.e-4_wp) * ( u1**2 + v1**2) )

  END FUNCTION calc_ustar

  !-------------------------------------------------------------------------
  !!
  !! Calculate convective contribution to the wind gusts
  !!     gust_conv = \alpha MAX(0,U_850 - U_950)
  !! where \alpha=0.6 is a tunable constant and U_850-U_950 is the difference between
  !! the 850 hPa and 950 hPa wind speeds, which represents the low-level wind shear.
  !!
  !! Literature
  !! Bechthold, P. and J. Bidlot (2009): Parameterization of convective gusts. 
  !! ECMWF Newsletter No. 119
  !!
  ELEMENTAL FUNCTION nwp_con_gust( u_850, u_950, v_850, v_950) RESULT(vgust_con)
    !$ACC ROUTINE SEQ

    REAL(wp), INTENT(IN) :: u_850, &    ! zonal wind component at 850 hPa [m/s]
      &                     u_950, &    ! zonal wind component at 950 hPa [m/s]
      &                     v_850, &    ! meridional wind component at 850 hPa [m/s]
      &                     v_950       ! meridional wind component at 950 hPa [m/s]

    REAL(wp) :: vgust_con               ! convective contribution to the wind gusts [m/s]

    REAL(wp), PARAMETER :: alpha = 0.6_wp ! convective mixing parameter

    vgust_con = alpha * MAX(0._wp, SQRT((u_850**2 + v_850**2)) - SQRT((u_950**2 + v_950**2)))

  END FUNCTION nwp_con_gust



  !-------------
  !>
  !! SUBROUTINE virtual_temp
  !! Computes virtual temperature
  !!
  !! Required input fields: temperature, specific humidity, cloud and precipitation variables
  !! Output: virtual temperature
  !!
  SUBROUTINE virtual_temp(p_patch, temp, qv, qc, qi, qr, qs, qg, qh, temp_v)


    TYPE(t_patch), INTENT(IN) :: p_patch

    ! Input fields - all defined at full model levels
    REAL(wp), INTENT(IN)                   :: temp(:,:,:) ! temperature (K)
    REAL(wp), INTENT(IN)                   :: qv  (:,:,:) ! specific humidity
    REAL(wp), INTENT(IN), OPTIONAL, TARGET :: qc  (:,:,:) ! specific cloud water
    REAL(wp), INTENT(IN), OPTIONAL, TARGET :: qi  (:,:,:) ! specific cloud ice
    REAL(wp), INTENT(IN), OPTIONAL, TARGET :: qr  (:,:,:) ! specific rain water
    REAL(wp), INTENT(IN), OPTIONAL, TARGET :: qs  (:,:,:) ! specific snow
    REAL(wp), INTENT(IN), OPTIONAL, TARGET :: qg  (:,:,:) ! specific graupel
    REAL(wp), INTENT(IN), OPTIONAL, TARGET :: qh  (:,:,:) ! specific hail

    REAL(wp), INTENT(OUT) :: temp_v(:,:,:) ! virtual temperature (K)

    INTEGER :: jb, jk, jc, jt
    INTEGER :: nlen, nlev
    INTEGER :: num_qcpvars ! number of cloud or precipitation variables
    REAL(wp):: z_qsum(nproma,SIZE(temp,2))

    TYPE t_fieldptr
      REAL(wp), POINTER :: fld(:,:,:)
    END TYPE t_fieldptr
    TYPE(t_fieldptr) :: qptr(6)

    nlev = SIZE(temp,2) ! in order to be usable for input and output data

    num_qcpvars = 0
    IF (PRESENT(qc)) THEN
      num_qcpvars = num_qcpvars + 1
      qptr(num_qcpvars)%fld => qc
    ENDIF
    IF (PRESENT(qr)) THEN
      num_qcpvars = num_qcpvars + 1
      qptr(num_qcpvars)%fld => qr
    ENDIF
    IF (PRESENT(qi)) THEN
      num_qcpvars = num_qcpvars + 1
      qptr(num_qcpvars)%fld => qi
    ENDIF
    IF (PRESENT(qs)) THEN
      num_qcpvars = num_qcpvars + 1
      qptr(num_qcpvars)%fld => qs
    ENDIF
    IF (PRESENT(qg)) THEN
      num_qcpvars = num_qcpvars + 1
      qptr(num_qcpvars)%fld => qg
    ENDIF
    IF (PRESENT(qh)) THEN
      num_qcpvars = num_qcpvars + 1
      qptr(num_qcpvars)%fld => qh
    ENDIF


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,jt,z_qsum) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, p_patch%nblks_c
      IF (jb /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      ENDIF
      
      z_qsum(:,:) = 0._wp
      IF (num_qcpvars > 0) THEN
        DO jt = 1, num_qcpvars
          DO jk = 1, nlev
            DO jc = 1, nlen
              z_qsum(jc,jk) = z_qsum(jc,jk) + qptr(jt)%fld(jc,jk,jb)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      DO jk = 1, nlev
        DO jc = 1, nlen
          temp_v(jc,jk,jb) = temp(jc,jk,jb) * (1._wp + vtmpc1*qv(jc,jk,jb) - z_qsum(jc,jk))
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE virtual_temp



  !> POINTWISE computation of shortwave direct downward flux as
  !! swdir_s = (1+ (\alpha/(1-\alpha))) * sobs - swdifd_s
  !! in W m**-2
  !!
  !! (domain independent and elemental)
  !!
  ELEMENTAL FUNCTION swdir_s(albedo, swdifd_s, sobs)
    !$ACC ROUTINE SEQ
    REAL(wp)             :: swdir_s
    REAL(wp), INTENT(IN) :: albedo      ! shortwave broadband albedo
    REAL(wp), INTENT(IN) :: swdifd_s    ! shortwave diffuse downward flux (sfc)
    REAL(wp), INTENT(IN) :: sobs        ! shortwave net flux (sfc)

    swdir_s = (1._wp + albedo/(1._wp - albedo)) * sobs - swdifd_s

  END FUNCTION swdir_s 


  !> POINTWISE computation of relative humidity as r=100. * e/e_sat,
  !! according to WMO standard
  !!
  !! (domain independent and elemental)
  !!
  ELEMENTAL FUNCTION rel_hum(temp, qv, p_ex)
    !$ACC ROUTINE SEQ

    REAL(wp) :: rel_hum
    REAL(wp), INTENT(IN) :: temp, &  ! temperature
      &                     qv,   &  ! spec. water vapor content
      &                     p_ex     ! exner pressure
    ! local variables
    REAL(wp) :: pres, e_s, e

    ! compute dynamic pressure from Exner pressure:
    pres = p0ref * EXP((cpd/rd)*LOG(p_ex))
    ! approx. saturation vapor pressure:
    e_s = sat_pres_water(temp)
    ! compute vapor pressure from formula for specific humidity:
    e   = pres*qv / (rdv + o_m_rdv*qv)

    rel_hum = 100._wp * e/e_s

  END FUNCTION rel_hum


  !> POINTWISE computation of relative humidity as r=100. * e/e_sat, 
  !! according to IFS documentation
  !! I.e. For the temperature range 250.16<=T<=273.16, the saturation 
  !! vapour pressure is computed as a combination of the values over 
  !! water e_s_water and over ice e_s_ice.
  !!
  !! (domain independent and elemental)
  !!
  ELEMENTAL FUNCTION rel_hum_ifs(temp, qv, p_ex)
    !$ACC ROUTINE SEQ

    REAL(wp) :: rel_hum_ifs
    REAL(wp), INTENT(IN) :: temp, &  ! temperature
      &                     qv,   &  ! spec. water vapor content
      &                     p_ex     ! exner pressure
    ! local variables
    REAL(wp) :: pres, e_s, e_s_water, e_s_ice, e
    REAL(wp), PARAMETER:: t_i = 250.16_wp  ! threshold value for mixed-phase clouds

    ! compute dynamic pressure from Exner pressure:
    pres = p0ref * EXP((cpd/rd)*LOG(p_ex))
    ! approx. saturation vapor pressure:
    IF (temp > t3) THEN
      e_s       = sat_pres_water(temp)
    ELSE IF (temp < (t3-23._wp)) THEN
      e_s       = sat_pres_ice(temp)
    ELSE
      e_s_water = sat_pres_water(temp)
      e_s_ice   = sat_pres_ice(temp)

      e_s       = e_s_ice + (e_s_water - e_s_ice) * ((temp - t_i)/(t3 - t_i))**2
    ENDIF

    ! compute vapor pressure from formula for specific humidity:
    e   = pres*qv / (rdv + o_m_rdv*qv)

    rel_hum_ifs = 100._wp * e/e_s

  END FUNCTION rel_hum_ifs



  !> computation of relative humidity as r=e/e_sat, according to WMO standard
  !!
  SUBROUTINE compute_field_rel_hum_wmo(ptr_patch, p_prog, p_diag, out_var, &
    &                              opt_slev, opt_elev, opt_rlstart, opt_rlend, lacc)

    ! patch on which computation is performed:
    TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
    ! nonhydrostatic state
    TYPE(t_nh_prog), INTENT(IN) :: p_prog
    TYPE(t_nh_diag), INTENT(IN) :: p_diag
    ! output variable, dim: (nproma,nlev,nblks_c):
    REAL(wp), TARGET, INTENT(INOUT)   :: out_var(:,:,:)
    ! optional vertical start/end level:
    INTEGER, INTENT(in), OPTIONAL     :: opt_slev, opt_elev
    ! start and end values of refin_ctrl flag:
    INTEGER, INTENT(in), OPTIONAL     :: opt_rlstart, opt_rlend
    LOGICAL, INTENT(IN), OPTIONAL     :: lacc ! If true, use openacc

    ! local variables
    REAL(wp) :: temp, qv, p_ex
    INTEGER  :: slev, elev, rl_start, rl_end, i_nchdom,     &
      &         i_startblk, i_endblk, i_startidx, i_endidx, &
      &         jc, jk, jb
    LOGICAL :: lzacc ! non-optional version of lacc

    CALL set_acc_host_or_device(lzacc, lacc)

    ! default values
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
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc,temp,qv,p_ex), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        i_startidx, i_endidx, rl_start, rl_end)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(temp, qv, p_ex)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

!!$ UB: do we need p_prog_rcf%tracer_ptr(iqv) instead of p_prog%tracer_ptr(iqv)?
          ! get values for temperature, etc.:
          temp = p_diag%temp(jc,jk,jb)
          qv   = p_prog%tracer_ptr(iqv)%p_3d(jc,jk,jb)
          p_ex = p_prog%exner(jc,jk,jb)
          !-- compute relative humidity as r = e/e_s:
          out_var(jc,jk,jb) = rel_hum(temp, qv, p_ex)

        END DO
      END DO
      !$ACC END PARALLEL

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_rel_hum_wmo



  !> computation of relative humidity as r=e/e_sat, according to IFS
  !!
  SUBROUTINE compute_field_rel_hum_ifs(ptr_patch, p_prog, p_diag, out_var, &
    &                              opt_lclip, opt_slev, opt_elev,          &
    &                              opt_rlstart, opt_rlend)

    ! patch on which computation is performed:
    TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
    ! nonhydrostatic state
    TYPE(t_nh_prog), INTENT(IN) :: p_prog
    TYPE(t_nh_diag), INTENT(IN) :: p_diag
    ! output variable, dim: (nproma,nlev,nblks_c):
    REAL(wp),INTENT(INOUT) :: out_var(:,:,:)
    ! optional clipping to rh<=100%
    LOGICAL, INTENT(IN), OPTIONAL     :: opt_lclip
    ! optional vertical start/end level:
    INTEGER, INTENT(in), OPTIONAL     :: opt_slev, opt_elev
    ! start and end values of refin_ctrl flag:
    INTEGER, INTENT(in), OPTIONAL     :: opt_rlstart, opt_rlend
   
    ! local variables
    REAL(wp) :: temp, qv, p_ex
    INTEGER  :: slev, elev, rl_start, rl_end, i_nchdom,     &
      &         i_startblk, i_endblk, i_startidx, i_endidx, &
      &         jc, jk, jb
    LOGICAL  :: lclip       ! clip rel. hum. to values <=100% 

#ifdef _OPENACC
    CALL finish ('mo_util_phys:compute_field_rel_hum_ifs', 'OpenACC version currently not implemented')
#endif

    IF (PRESENT(opt_lclip)) THEN
      lclip = opt_lclip
    ELSE
      lclip = .FALSE.
    ENDIF


    ! default values
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
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc,temp,qv,p_ex), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        i_startidx, i_endidx, rl_start, rl_end)
      
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

!!$ UB: do we need p_prog_rcf%tracer_ptr(iqv) instead of p_prog%tracer_ptr(iqv)?
          ! get values for temperature, etc.:
          temp = p_diag%temp(jc,jk,jb)
          qv   = p_prog%tracer_ptr(iqv)%p_3d(jc,jk,jb)
          p_ex = p_prog%exner(jc,jk,jb)
          !-- compute relative humidity as r = e/e_s:
          out_var(jc,jk,jb) = rel_hum_ifs(temp, qv, p_ex)

          ! optional clipping, if lclip=.TRUE.
          out_var(jc,jk,jb) = MERGE(MIN(100._wp,out_var(jc,jk,jb)), out_var(jc,jk,jb), lclip)

        END DO
      END DO

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_rel_hum_ifs



  !> computation of water vapour pressure
  !!
  !! water vapour pressure is computed as a function of specific humidity 
  !! qv and atmospheric pressure pres.
  !!
  ELEMENTAL FUNCTION vap_pres(qv,pres)
  !$ACC ROUTINE SEQ

  IMPLICIT NONE

    REAL(wp), INTENT(IN)  :: qv   ! specific humidity         [kg/kg]
    REAL(wp), INTENT(IN)  :: pres ! atmospheric pressure      [Pa]
    REAL(wp) :: vap_pres          ! water vapour pressure     [Pa]

    vap_pres = (qv * pres) / (rdv + O_m_rdv*qv)

  END FUNCTION vap_pres


  !
  ! Add slow-physics tendencies to tracer fields
  !
  ! Add slow-physics tendencies to tracer fields. Currently, 
  ! convection is the only slow-physics routine which provides tracer 
  ! tendencies.
  ! In addition, this routine
  ! - makes sure that tendencies from advection and/or convection 
  !   do not result in negative mass fractions. If negative values in qx  
  !   occur, these are clipped. The moisture which is spuriously created by this 
  !   clipping is substracted from qv.
  ! - Diagnoses amount of convective rain and snow (rain_con, snow_con), 
  !   as well as the total convective precipitation (prec_con).
  ! - applies large-scale-forcing tendencies, if ICON is run in single-column-mode.
  ! 
  SUBROUTINE tracer_add_phytend( p_rho_now, prm_nwp_tend, pdtime, prm_diag, &
    &                            pt_prog_rcf, p_metrics, dt_loc, jg, jb, i_startidx, i_endidx, kend, lacc)

    REAL(wp), CONTIGUOUS,INTENT(IN)   :: p_rho_now(:,:)  !< total air density
    TYPE(t_nwp_phy_tend),INTENT(IN)   :: prm_nwp_tend    !< atm tend vars
    REAL(wp)            ,INTENT(IN)   :: pdtime          !< time step
    TYPE(t_nwp_phy_diag),INTENT(INOUT):: prm_diag        !< the physics variables
    TYPE(t_nh_prog)     ,INTENT(INOUT):: pt_prog_rcf     !< the tracer field at
                                                         !< reduced calling frequency
    TYPE(t_nh_metrics)  ,INTENT(IN)   :: p_metrics       !< NH metrics variables
    REAL(wp)            ,INTENT(IN)   :: dt_loc          !< (advective) time step applicable to local grid level
    INTEGER             ,INTENT(IN)   :: jg              !< domain ID
    INTEGER             ,INTENT(IN)   :: jb              !< block index
    INTEGER             ,INTENT(IN)   :: i_startidx, i_endidx
    INTEGER             ,INTENT(IN)   :: kend            !< vertical end index                             
    LOGICAL, OPTIONAL   ,INTENT(IN)   :: lacc            ! If true, use openacc

    ! Local variables
    INTEGER  :: jt          ! tracer loop index
    INTEGER  :: idx         ! tracer position in container
    INTEGER  :: pos_qv      ! position of qv in local array zrhox
    INTEGER  :: jk,jc
    INTEGER  :: iq_start
    REAL(wp) :: zrhox(nproma,kend,5)
    REAL(wp) :: zrhox_clip(nproma,kend) !< Negative sum of clipped water tracer density [kg/m3].
    REAL(wp) :: zwtr_clip_rate(nproma) !< Total negative sum of clipped vapor mass rate [kg/m2/s].
    REAL(wp) :: zrhoqv
#ifndef __NO_ICON_LES__
    REAL(wp) :: nudgecoeff  ! SCM Nudging
    REAL(wp) :: z_ddt_q_nudge
#endif
    !
    INTEGER, DIMENSION(5) :: conv_list
    LOGICAL :: lzacc ! non-optional version of lacc
    REAL(wp) :: zxiconv

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA PRESENT(p_rho_now, prm_nwp_tend, prm_diag, pt_prog_rcf) &
    !$ACC   CREATE(zrhox, zrhox_clip, zwtr_clip_rate) &
    !$ACC   COPYIN(kstart_moist) &
    !$ACC   IF(lzacc)


    ! get list of water tracers which are affected by convection
    IF (atm_phy_nwp_config(jg)%ldetrain_conv_prec) THEN
      conv_list = (/iqv,iqc,iqi,iqr,iqs/)
    ELSE
      conv_list = (/iqv,iqc,iqi,-1,-1/)
    ENDIF
    !$ACC DATA COPYIN(conv_list) IF(lzacc)

    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    zrhox_clip(:,:) = 0._wp
    zwtr_clip_rate(:) = 0._wp
    !$ACC END KERNELS

    ! add tendency due to convection
    !$ACC PARALLEL PRIVATE(pos_qv) DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP SEQ
    DO jt=1,SIZE(conv_list)
      idx = conv_list(jt)
      IF (idx <= 0) CYCLE

      !$ACC LOOP SEQ
      DO jk = kstart_moist(jg), kend
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = i_startidx, i_endidx
          zrhox(jc,jk,jt) = p_rho_now(jc,jk)*pt_prog_rcf%tracer(jc,jk,jb,idx)  &
            &             + pdtime*prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,idx)

          ! keep mass that is created due to artificial clipping
          zrhox_clip(jc,jk) = zrhox_clip(jc,jk) + MIN(0._wp,zrhox(jc,jk,jt))

          ! clip
          zrhox(jc,jk,jt) = MAX(0._wp, zrhox(jc,jk,jt))
        ENDDO
      ENDDO
      !
      ! Re-diagnose tracer mass fraction from partial mass
      IF (idx == iqv) THEN
        pos_qv = jt   ! store local qv-position for later use
        CYCLE         ! special treatment see below
      ENDIF
      !
      !$ACC LOOP SEQ
      DO jk = kstart_moist(jg), kend
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = i_startidx, i_endidx
          pt_prog_rcf%tracer(jc,jk,jb,idx) = zrhox(jc,jk,jt)/p_rho_now(jc,jk)
        ENDDO
      ENDDO
    ENDDO ! jt
    !
    ! Special treatment for qv.
    ! Rediagnose tracer mass fraction and substract mass created by artificial clipping.
    !$ACC LOOP SEQ
    DO jk = kstart_moist(jg), kend
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zrhoqv)
      DO jc = i_startidx, i_endidx
        zrhoqv = zrhox(jc,jk,pos_qv) + zrhox_clip(jc,jk)
        ! Keep total mass of clipped water vapor. Used to reduce convective rain and snow.
        zwtr_clip_rate(jc) = zwtr_clip_rate(jc) &
          & + p_metrics%ddqz_z_full(jc,jk,jb) / pdtime * MIN(0._wp, zrhoqv)
        pt_prog_rcf%tracer(jc,jk,jb,iqv) = MAX(0._wp, zrhoqv/p_rho_now(jc,jk))
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO jc = i_startidx, i_endidx
      ! First remove clipped mass from rain, then from snow.
      prm_diag%rain_con_rate_corr(jc,jb) = MAX(0._wp, prm_diag%rain_con_rate(jc,jb) + zwtr_clip_rate(jc))
      zwtr_clip_rate(jc) = MIN(0._wp, zwtr_clip_rate(jc) + prm_diag%rain_con_rate(jc,jb))
      prm_diag%snow_con_rate_corr(jc,jb) = MAX(0._wp, prm_diag%snow_con_rate(jc,jb) + zwtr_clip_rate(jc))
    END DO
    !$ACC END PARALLEL

!!  Update of two-moment number densities using the updates from the convective parameterization
    IF (atm_phy_nwp_config(jg)%l2moment) THEN
      DO jt=1,SIZE(conv_list)
        idx = conv_list(jt)
        IF ( idx > 0 .AND. idx /= iqv .AND. idx /= iqc .AND. idx /= iqi .AND. idx /= iqr .AND. idx /= iqs ) THEN
          CALL finish("mo_util_phys", "Not valid update from convective parameterization for two-moment scheme.")
        ENDIF
      ENDDO

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO jt=1,SIZE(conv_list)
        idx = conv_list(jt)
        IF (idx <= 0 .OR. idx == iqv) CYCLE
     
        IF ( idx == iqc ) THEN
          !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
          DO jk = kstart_moist(jg), kend
            DO jc = i_startidx, i_endidx              
              pt_prog_rcf%tracer(jc,jk,jb,iqnc) =  pt_prog_rcf%tracer(jc,jk,jb,iqnc)    &
                 + set_qnc( pdtime *  prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqc) )/p_rho_now(jc,jk)   
            END DO
          END DO
        ELSEIF ( idx == iqi ) THEN
          !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
          DO jk = kstart_moist(jg), kend
            DO jc = i_startidx, i_endidx              
              pt_prog_rcf%tracer(jc,jk,jb,iqni) =  pt_prog_rcf%tracer(jc,jk,jb,iqni)    &
                 + set_qni( pdtime *  prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqi) )/p_rho_now(jc,jk)   
            END DO
          END DO
        ELSEIF ( idx == iqr ) THEN
          !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
          DO jk = kstart_moist(jg), kend
            DO jc = i_startidx, i_endidx              
              pt_prog_rcf%tracer(jc,jk,jb,iqnr) =  pt_prog_rcf%tracer(jc,jk,jb,iqnr)    &
                 + set_qnr( pdtime *  prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqr) )/p_rho_now(jc,jk)   
            END DO
          END DO
        ELSEIF ( idx == iqs ) THEN
          !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
          DO jk = kstart_moist(jg), kend
            DO jc = i_startidx, i_endidx              
              pt_prog_rcf%tracer(jc,jk,jb,iqns) =  pt_prog_rcf%tracer(jc,jk,jb,iqns)    &
                 + set_qns( pdtime *  prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqs) )/p_rho_now(jc,jk)   
            END DO
          END DO
        END IF
      END DO
      !$ACC END PARALLEL
    ELSEIF (atm_phy_nwp_config(jg)%inwp_gscp == 3) THEN
      CALL get_mean_crystal_mass(zxiconv)
      !CALL assert_acc_host_only("tracer_add_phytend l2moment", lacc)  ! some GPU stuff?
      DO jt=1,SIZE(conv_list)
        idx = conv_list(jt)
        IF ( idx == iqi ) THEN
          DO jk = kstart_moist(jg), kend
            DO jc = i_startidx, i_endidx
              pt_prog_rcf%tracer(jc,jk,jb,iqni) = pt_prog_rcf%tracer(jc,jk,jb,iqni)    &
                 + pdtime * prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqi)/(zxiconv*p_rho_now(jc,jk))
            END DO
          END DO
        END IF
      END DO
    END IF

    IF(lart .AND. art_config(jg)%lart_conv) THEN
      ! add convective tendency and fix to positive values
      DO jt=1,art_config(jg)%nconv_tracer  ! ASH
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1, kend
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            pt_prog_rcf%conv_tracer(jb,jt)%ptr(jc,jk)=MAX(0._wp,pt_prog_rcf%conv_tracer(jb,jt)%ptr(jc,jk) &
               +pdtime*prm_nwp_tend%conv_tracer_tend(jb,jt)%ptr(jc,jk)/p_rho_now(jc,jk))
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDDO
    ENDIF !lart

    ! additional clipping for qr, qs, ... up to iqm_max
    ! (very small negative values may occur during the transport process (order 10E-15))
    iq_start = MAXVAL(conv_list(:)) + 1  ! all others have already been clipped above
    !
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(3)
    DO jt=iq_start, iqm_max  ! qr,qs,etc. 
      DO jk = kstart_moist(jg), kend
        DO jc = i_startidx, i_endidx
          pt_prog_rcf%tracer(jc,jk,jb,jt) = MAX(0._wp, pt_prog_rcf%tracer(jc,jk,jb,jt))
        ENDDO
      ENDDO
    ENDDO
    !$ACC END PARALLEL
    
    ! clipping for number concentrations
    IF(atm_phy_nwp_config(jg)%l2moment)THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(3)
      DO jt=iqni, ininact  ! qni,qnr,qns,qng,qnh,qnc and ninact (but not yet ninpot)
        DO jk = kstart_moist(jg), kend
          DO jc = i_startidx, i_endidx
            pt_prog_rcf%tracer(jc,jk,jb,jt) = MAX(0._wp, pt_prog_rcf%tracer(jc,jk,jb,jt))
          ENDDO          
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    END IF

    ! clipping for mass-bins
    IF(atm_phy_nwp_config(jg)%lsbm)THEN
      CALL assert_acc_host_only("tracer_add_phytend lsbm", lacc)
      DO jt = iqb_i, iqb_e
        DO jk = kstart_moist(jg), kend
          DO jc = i_startidx, i_endidx
            pt_prog_rcf%tracer(jc,jk,jb,jt) = MAX(0._wp, pt_prog_rcf%tracer(jc,jk,jb,jt))
          ENDDO
        ENDDO
      ENDDO
    END IF

    ! Diagnose convective precipitation amount
    IF (atm_phy_nwp_config(jg)%lcalc_acc_avg) THEN
!DIR$ IVDEP
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx

        prm_diag%rain_con(jc,jb) = prm_diag%rain_con(jc,jb)    &
          &                      + pdtime * prm_diag%rain_con_rate(jc,jb)

        prm_diag%snow_con(jc,jb) = prm_diag%snow_con(jc,jb)    &
          &                      + pdtime * prm_diag%snow_con_rate(jc,jb)

        prm_diag%prec_con(jc,jb) = prm_diag%rain_con(jc,jb) + prm_diag%snow_con(jc,jb)

        ! to compute tot_prec_d lateron:
        prm_diag%prec_con_d(jc,jb) = prm_diag%prec_con_d(jc,jb) + pdtime * ( &
          &                          prm_diag%rain_con_rate(jc,jb) + prm_diag%snow_con_rate(jc,jb) )

      ENDDO
      !$ACC END PARALLEL
    ENDIF

#ifndef __NO_ICON_LES__
    ! Add LS forcing to moisture variable including nudging
    IF(is_ls_forcing)THEN
      CALL assert_acc_host_only("tracer_add_phytend is_ls_forcing", lacc)
      DO jt=1, nqtendphy  ! qv,qc,qi
        DO jk = kstart_moist(jg), kend
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx

            ! add q nudging (T, U, V nudging is in mo_nh_interface_nwp)
            IF ( is_nudging_tq ) THEN

              ! linear nudging profile between "start" and "full" heights - prevent sfc layer instability
              IF ( nudge_full_height == nudge_start_height ) THEN
                nudgecoeff = 1.0_wp
              ELSE
                nudgecoeff = ( p_metrics%geopot_agl(jc,jk,jb)/grav - nudge_start_height ) / &
                           & ( nudge_full_height                   - nudge_start_height )
                nudgecoeff = MAX( MIN( nudgecoeff, 1.0_wp ), 0.0_wp )
              END IF
              ! analytic implicit: (q,n+1 - q,n) / dt = (q,nudge - q,n) / dt_relax * exp(-dt/dt_relax)
              z_ddt_q_nudge =                                                          &
                &  - ( pt_prog_rcf%tracer(jc,jk,jb,jt) - prm_nwp_tend%q_nudge(jk,jt) ) &
                &  / dt_relax * exp(-dt_loc/dt_relax) * nudgecoeff
            ELSE
              z_ddt_q_nudge = 0.0_wp
            END IF

            pt_prog_rcf%tracer(jc,jk,jb,jt) = MAX(0._wp, pt_prog_rcf%tracer(jc,jk,jb,jt)   &
              &                                 + pdtime*prm_nwp_tend%ddt_tracer_ls(jk,jt) &
              &                                 + pdtime*z_ddt_q_nudge)

          ENDDO
        ENDDO
      END DO
    ENDIF  ! is_ls_forcing
#endif

    !$ACC WAIT
    !$ACC END DATA ! COPYIN(conv_list)
    !$ACC END DATA ! DATA PRESENT

  END SUBROUTINE tracer_add_phytend

  !>
  !! Find the lowest inversion and provide its inversion height and lowest point of the entrainment zone 
  !! It follows Van Wevweberg et al. Month Weath. Rev. 2021
  !!
  !! The inversion height is identified as the maximum gradient of liquid potential temperature
  
  SUBROUTINE inversion_height_index(z,zsurf,qc,te,prs,i_startidx,i_endidx,jktop,jkbot,nlev, &
                      &             i_inversion,i_ent_zone,lfound_inversion)
    REAL(wp),      INTENT(IN)  ::  z(:,:)     ! Height above sea level
    REAL(wp),      INTENT(IN)  ::  zsurf(:)   ! Surface height above sea level
    REAL(wp),      INTENT(IN)  ::  qc(:,:)  ! Liquid water
    REAL(wp),      INTENT(IN)  ::  te(:,:)  ! Temperature
    REAL(wp),      INTENT(IN)  ::  prs(:,:) ! Pressure
    INTEGER,       INTENT(IN)  ::  i_startidx,i_endidx,jktop,jkbot,nlev ! loop indices

    INTEGER,       INTENT(OUT) ::  i_inversion(nproma) ! Inversion index
    INTEGER,       INTENT(OUT) ::  i_ent_zone(nproma)  ! Lowest inversion index
    LOGICAL,       INTENT(OUT) ::  lfound_inversion(nproma) ! Inversion found (true/false) 

    REAL (wp),      PARAMETER  ::   p0 = 1.e5_wp    ! reference pressure for calculation of potential temperature
    REAL (wp),     PARAMETER   ::  zmin_inv =  400.0_wp  ! Lowest possible inversion (in m Above Surface)
    REAL (wp),     PARAMETER   ::  zmax_inv = 3000.0_wp  ! Highest possible inversion(in m Above Surface)

    ! Local variables
    REAL(wp) ::   theta_l(nproma,nlev)    ! Liquid potential temperature
    REAL(WP) ::   dthetadz(nproma,3)      ! Gradiente of liquid potential temperature
    LOGICAL  ::   lbelow_zmax(nproma)     ! Height below zmax_inv       

    REAL    ::     lapse_lim              ! Stratification limit to be considered as a inversion    
    INTEGER ::     jc,jk


    ! Limit to be in the entrainment zone (Van Wevweberg et al. Month Weath. Rev. 2021)
    lapse_lim = grav/cpd*0.1_wp            

    ! Start arrays  
    theta_l(:,:)   = 0.0_wp
    i_ent_zone(:)  = jkbot-3
    i_inversion(:) = jkbot-3
    lfound_inversion(:) = .false.
    lbelow_zmax(:)      = .true.
        
    ! Calculate the liquid potential temperature (constant latent heat approximation)
    DO jk = jktop,jkbot
      DO jc = i_startidx, i_endidx            
        theta_l(jc,jk) = (te(jc,jk) - alvdcp *qc(jc,jk))*(prs(jc,jk)/p0)**rd_o_cpd
      END DO
    END DO


! Lowest two levels
    DO jc = i_startidx, i_endidx
      dthetadz(jc,1) = (theta_l(jc,jkbot-2) - theta_l(jc,jkbot  ) ) / (z(jc,jkbot-2) - z(jc,jkbot  )) 
      dthetadz(jc,2) = (theta_l(jc,jkbot-3) - theta_l(jc,jkbot-1) ) / (z(jc,jkbot-3) - z(jc,jkbot-1)) 
    END DO

! Loop from bottom to top
    DO jk = jkbot-3, jktop+1 ,-1
      DO jc = i_startidx, i_endidx
        ! Calculate when the inversion has not been found and below max z level
        IF ( lbelow_zmax(jc) .AND. .NOT. lfound_inversion(jc) ) THEN 
          dthetadz(jc,3) = (theta_l(jc,jk-1) - theta_l(jc,jk+1) ) / (z(jc,jk-1) - z(jc,jk+1))
          ! Criteria for entrainment zone
          IF ( dthetadz(jc,2) > lapse_lim .AND. z(jc,jk+1) > (zmin_inv + zsurf(jc)) ) THEN
            ! Maximum: criteria for inversion height
            IF ( dthetadz(jc,2) > dthetadz(jc,3) .AND. dthetadz(jc,2) > dthetadz(jc,1) ) THEN
              lfound_inversion(jc) = .true.
              i_inversion(jc) = MAX(MIN(jk + 1,jkbot),1)
            END IF
          ELSE
            ! We are not in the entrainment zone, we shift the limit upwards
            i_ent_zone(jc) = MAX(MIN(jk,jkbot),1)
          END IF
          ! Shift the array
          dthetadz(jc,1) = dthetadz(jc,2)
          dthetadz(jc,2) = dthetadz(jc,3)
          IF ( z(jc,jk) > (zmax_inv + zsurf(jc)) ) THEN
            lbelow_zmax(jc) = .false.
          END IF
        END IF
      END DO
    END DO


  END SUBROUTINE inversion_height_index
  
END MODULE mo_util_phys
