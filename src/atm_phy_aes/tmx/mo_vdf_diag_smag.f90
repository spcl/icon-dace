!
! Classes and functions for the turbulent mixing package (tmx)
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

MODULE mo_vdf_diag_smag

  USE mo_kind,              ONLY: wp, vp, i1
  USE mo_exception,         ONLY: message, finish
  USE mo_tmx_process_class, ONLY: t_tmx_process
  USE mo_tmx_field_class,   ONLY: t_tmx_field, t_domain, isfc_oce, isfc_ice, isfc_lnd
  USE mo_util_string,       ONLY: real2string
  USE mo_model_domain,      ONLY: t_patch
  USE mo_intp_data_strc,    ONLY: t_int_state, p_int_state
  USE mo_nonhydro_types,    ONLY: t_nh_metrics
  USE mo_nonhydro_state,    ONLY: p_nh_state
  USE mo_aes_sfc_indices,   ONLY: nsfc_type, iwtr, iice, ilnd
  ! USE mo_physical_constants,ONLY: grav, rd, vtmpc1, p0ref, rgrav, alv, als
  USE mo_physical_constants,ONLY: grav, rd, vtmpc1, rgrav
  USE mo_turb_vdiff_params, ONLY: ckap
  USE mo_math_constants,    ONLY: pi_2, ln2
  USE mo_impl_constants,    ONLY: min_rlcell, min_rledge_int, min_rlcell_int, min_rlvert_int
  USE mo_nh_vert_interp_les,ONLY: brunt_vaisala_freq, vert_intp_full2half_cell_3d
  USE mo_intp,              ONLY: cells2verts_scalar, cells2edges_scalar
  USE mo_sync,              ONLY: SYNC_E, SYNC_C, SYNC_V, sync_patch_array,     &
    &                             sync_patch_array_mult
  USE mo_loopindices,       ONLY: get_indices_e, get_indices_c
  USE mo_impl_constants_grf,ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_intp_rbf,          ONLY: rbf_vec_interpol_edge
  USE mo_fortran_tools,     ONLY: init

  USE mo_aes_thermo,        ONLY: potential_temperature, sat_pres_water, sat_pres_ice, specific_humidity
  ! USE mo_jsb_interface,     ONLY: jsbach_get_var

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: sfc_exchange_coefficients, compute_sfc_exchange_coefficients, &
    &       compute_wind_speed, &
    &       compute_atm_potential_temperature, compute_moist_richardson, &
    ! &       compute_sfc_sat_spec_humidity, compute_sfc_fluxes, compute_sfc_roughness &
    &       compute_sfc_potential_temperature, compute_sfc_density

  !Parameters for surface layer parameterizations: From Zeng_etal 1997 J. Clim
  REAL(wp), PARAMETER :: bsm = 5.0_wp  !Businger Stable Momentum
  REAL(wp), PARAMETER :: bum = 16._wp  !Businger Untable Momentum
  REAL(wp), PARAMETER :: bsh = 5.0_wp  !Businger Stable Heat
  REAL(wp), PARAMETER :: buh = 16._wp  !Businger Untable Heat

  CHARACTER(len=*), PARAMETER :: modname = 'mo_vdf_diag_smag'

CONTAINS

  !
  !=================================================================
  !
  ! Subroutines for diagnostics
  !
  ! Note: will be re-factored
  !
  SUBROUTINE compute_sfc_density(                 &
    & domain,                                            &
    & nvalid, indices,                                   &
    & qsat_sfc, ppsfc, ptsfc,  &
    & rhos)

    ! Domain information
    TYPE(t_domain),        INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:)
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      & qsat_sfc, &
      ppsfc,     &
      ptsfc
    !
    ! Output variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: rhos

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_sfc_density'

    INTEGER  :: jb, jls, js
    
!$OMP PARALLEL
    CALL init(rhos, 1._wp)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c, domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1) PRIVATE(js)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)
        rhos(js,jb) = ppsfc(js,jb) / (rd * ptsfc(js,jb) * (1._wp + vtmpc1 * qsat_sfc(js,jb)))
        ! TODO: in vdiff sfc_exchange_coeff 1._wp+vtmpc1*pqm1_b(jl)-pqxm1_b(jl) is used
        !       and temperature at lowest model level instead of ptsfc !?

      END DO !jls
      !$ACC END PARALLEL LOOP
    END DO !jb
!$OMP END PARALLEL DO

    !$ACC WAIT(1)

  END SUBROUTINE compute_sfc_density  !
!
!=================================================================
!
  SUBROUTINE compute_wind_speed( &
    & domain,                  &
    & min_sfc_wind,            &
    & pum1, pvm1,              &
    & mwind                    &
    )

    ! Domain information
    TYPE(t_domain),        INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    REAL(wp), INTENT(in) :: min_sfc_wind
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      pum1 ,     &
      pvm1
    !
    ! Output variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: mwind

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_wind_speed'

    INTEGER  :: jb, jc
    
!$OMP PARALLEL
    CALL init(mwind)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c, domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
        mwind(jc,jb) = MAX( min_sfc_wind, SQRT(pum1(jc,jb)**2 + pvm1(jc,jb)**2) )
      END DO !jc
      !$ACC END PARALLEL LOOP
    END DO !jb
!$OMP END PARALLEL DO

  END SUBROUTINE compute_wind_speed
  !
  !=================================================================
  !
  SUBROUTINE compute_atm_potential_temperature( &
    & domain,                  &
    & ptm1, ptvm1, papm1,      &
    & theta, thetav            &
    )

    ! Domain information
    TYPE(t_domain),        INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
    & ptm1 , &
    & ptvm1, &
    & papm1
    !
    ! Output variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: theta, thetav

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_atm_potential_temperature'

    INTEGER  :: jb, jc

!$OMP PARALLEL
    CALL init(theta)
    CALL init(thetav)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c, domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
        theta (jc,jb) = potential_temperature(ptm1(jc,jb),  papm1(jc,jb))
        thetav(jc,jb) = potential_temperature(ptvm1(jc,jb), papm1(jc,jb))
      END DO !jc
      !$ACC END PARALLEL LOOP
    END DO !jb
!$OMP END PARALLEL DO

  END SUBROUTINE compute_atm_potential_temperature
  !
  !=================================================================
  !
  SUBROUTINE compute_sfc_potential_temperature( &
    & domain,                  &
    & nvalid, indices,         &
    & ppsfc, ptsfc, qsat,      &
    & theta, thetav            &
    )

    ! Domain information
    TYPE(t_domain),        INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:)
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      ptsfc,     &
      qsat,      &
      ppsfc
    !
    ! Output variables
    !
    REAL(wp), DIMENSION(:,:), INTENT(out) :: theta, thetav

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_sfc_potential_temperature'

    INTEGER  :: jb, jls, js

!$OMP PARALLEL
    CALL init(theta)
    CALL init(thetav)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jls, js) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)
        theta(js,jb)  = potential_temperature(ptsfc(js,jb), ppsfc(js,jb))
        thetav(js,jb) = theta(js,jb) * (1._wp + vtmpc1 * qsat(js,jb))
      END DO !jls
      !$ACC END PARALLEL LOOP
    END DO !jb
!$OMP END PARALLEL DO

  END SUBROUTINE compute_sfc_potential_temperature
  !
  !=================================================================
  !
  SUBROUTINE compute_moist_richardson(     &
    & domain,                    &
    & nvalid, indices,           &
    & fsl, zf, &
    & thetav_atm, thetav_sfc, wind, &
    & richardson_number)

    ! Domain information
    TYPE(t_domain),  INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:)
    REAL(wp), INTENT(in) :: &
      & fsl
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      & zf, thetav_atm, thetav_sfc, wind
    REAL(wp), DIMENSION(:,:), INTENT(out) :: &
      & richardson_number

    REAL(wp) :: zthetav_mid, w1, w2
    INTEGER :: jb, jls, js

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_moist_richardson'

    w1 = fsl
    w2 = 1._wp - fsl

!$OMP PARALLEL
    CALL init(richardson_number)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jls, js, zthetav_mid) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1) PRIVATE(js)
      DO jls = 1, nvalid(jb)
        js = indices(jls,jb)
        zthetav_mid = w1 * thetav_atm(js,jb) + w2 * thetav_sfc(js,jb)
        richardson_number(js,jb) = zf(js,jb) * grav * (thetav_atm(js,jb) - thetav_sfc(js,jb)) / (zthetav_mid * wind(js,jb))
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!$OMP END PARALLEL DO

  END SUBROUTINE compute_moist_richardson
  !
  !=================================================================
  !
! #ifndef _OPENACC
!   ELEMENTAL &
! #endif
  PURE SUBROUTINE sfc_exchange_coefficients(                   &
    & dz,                                                 &
    & pqm1, thetam1, mwind, rough_m, theta_sfc, qsat_sfc, &
    & km, kh, km_neutral, kh_neutral                      &
    & )

    REAL(wp), INTENT(in) :: &
      dz,        & ! height to be used as a reference height in surface layer
      thetam1,   &
      pqm1 ,     &
      mwind,     &
      rough_m,   &
      theta_sfc, &
      qsat_sfc
      !
    REAL(wp), INTENT(out) :: &
      kh,         &
      km,         &
      kh_neutral, &
      km_neutral

    !$ACC ROUTINE SEQ
  
    REAL(wp) :: z_mc, RIB, tcn_mom, tcn_heat, &
      &         shfl_local, lhfl_local, bflx1, ustar, obukhov_length, inv_bus_mom
    REAL(wp) :: tch
    REAL(wp) :: tcm

    INTEGER  :: itr

    ! to prevent floating-point arithmetic inconsistencies later in
    ! the interpolation to u 10m and 2m T/T_d: has been 0.01 before
    ! (Cray FP instead of IEEE 754 FP format)
    REAL(wp),parameter :: zepsec = 0.028_wp
    REAL(wp),parameter :: zcons17 = 1._wp / ckap**2

    z_mc = dz
    !First guess for tch and tcm using bulk approach
    RIB = grav * (thetam1-theta_sfc) * (z_mc-rough_m) / (theta_sfc * mwind**2)
    tcn_mom = (ckap/LOG(z_mc/rough_m))**2
    tcm     = tcn_mom * stability_function_mom(RIB,z_mc/rough_m,tcn_mom)

    tcn_heat        = ckap**2/(LOG(z_mc/rough_m)*LOG(z_mc/rough_m))
    tch             = tcn_heat * stability_function_heat(RIB,z_mc/rough_m,tcn_heat)

    !now iterate
    DO itr = 1, 5
      shfl_local = tch*mwind * (theta_sfc - thetam1)
      lhfl_local = tch*mwind * (qsat_sfc - pqm1)
      bflx1= shfl_local + vtmpc1 * theta_sfc * lhfl_local
      ustar= SQRT(tcm)*mwind

      obukhov_length = -ustar**3 * theta_sfc * rgrav / (ckap * bflx1)

      inv_bus_mom = 1._wp / businger_mom(rough_m,z_mc,obukhov_length)
      tch         = inv_bus_mom / businger_heat(rough_m,z_mc,obukhov_length)
      tcm         = inv_bus_mom * inv_bus_mom
    END DO

    ! pcfm = tcm*mwind
    ! pcfh = tch*mwind
    kh  = tch
    km  = tcm
    kh_neutral = ckap / MAX(zepsec, SQRT(tcn_heat))
    km_neutral = ckap / MAX(zepsec, SQRT(tcn_mom))

    ! pcfm = pcfm + pfrc*pcfm
    ! pcfh = pcfh + pfrc*pcfh

    ! pbn = ckap / MAX( zepsec, sqrt(tcn_mom) )
    ! pbhn= ckap / MAX( zepsec, sqrt(tcn_heat) )
    ! pbm = MAX( zepsec, sqrt(pcfm * tch*zcons17/ (tcn_mom*mwind)) )
    ! pbh = MAX( zepsec, tch/pbm*zcons17)
    ! pbm = 1._wp / pbm
    ! pbh = 1._wp / pbh
  
  END SUBROUTINE sfc_exchange_coefficients
  !
  !=================================================================
  !
  ! Initial release by Junhong Lee; MPI-M (2020-02)
  ! Rewritten for the new framework of ICON-Sapphire by Junhong Lee, Mira Shevchenko; MPI-M (2022-10)
  !
  SUBROUTINE compute_sfc_exchange_coefficients(          &
    & domain,                                            &
    & nvalid, indices, dz,                               &
    & pqm1,                                              &
    & thetam1, mwind, rough_m, theta_sfc, qsat_sfc,      &
    & km, kh, km_neutral, kh_neutral                     &
    )

    ! Domain information
    TYPE(t_domain),        INTENT(in), POINTER :: domain
    !
    ! Input variables
    !
    INTEGER,  INTENT(in)  :: &
      & nvalid(:),           &
      & indices(:,:)
    REAL(vp), DIMENSION(:,:), INTENT(in) :: &
      dz
    REAL(wp), DIMENSION(:,:), INTENT(in) :: &
      thetam1,   &
      pqm1 ,     &
      ! pxim1,     &
      ! ppsfc,     &
      ! ptsfc,     &
      ! pcsat,     &
      ! pcair,     &
      mwind,     &
      rough_m,   &
      theta_sfc, &
      qsat_sfc
      ! pch,       &
      ! pbn,       &
      ! pbhn,      &
      ! pbm,       &
      ! pbh,       &
      !
      ! Output variables
      !
      REAL(wp), DIMENSION(:,:), INTENT(out) :: &
        kh,         &
        km,         &
        kh_neutral, &
        km_neutral

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_sfc_exchange_coefficients'

    INTEGER  :: jb, jls, js
    REAL(wp) :: dz_temp
    
!$OMP PARALLEL
    CALL init(km)
    CALL init(kh)
    CALL init(km_neutral)
    CALL init(kh_neutral)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb,jls,js, dz_temp) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c, domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1) PRIVATE(js, dz_temp)
      DO jls = 1, nvalid(jb)
        js=indices(jls,jb)

        ! Note: dz=p_nh_metrics%ddqz_z_half(:,nlevp1,:) from atmosphere (src/atm_dyn_iconam/mo_vertical_grid.f90)
        ! is twice of what we want
        dz_temp = REAL(dz(js,jb), KIND=wp) * 0.5_wp
        CALL sfc_exchange_coefficients(                                                      &
          & dz_temp,                                                                         &
          & pqm1(js,jb),                                                                     &
          & thetam1(js,jb), mwind(js,jb), rough_m(js,jb), theta_sfc(js,jb), qsat_sfc(js,jb), &
          & km(js,jb), kh(js,jb), km_neutral(js,jb), kh_neutral(js,jb)                       &
          & )
      END DO !jls
      !$ACC END PARALLEL LOOP
    END DO !jb
!$OMP END PARALLEL DO

    !$ACC WAIT(1)

  END SUBROUTINE compute_sfc_exchange_coefficients  !
  !
  !=================================================================
  !
  !! stability_function_mom
  !! Taken from COSMO docs and Holstag & Boville 1992
  !!------------------------------------------------------------------------
#ifndef _OPENACC
  ELEMENTAL &
#endif
  PURE FUNCTION stability_function_mom(RIB, hz0, tc) RESULT(stab_fun)
     REAL(wp), INTENT(IN) :: RIB, hz0, tc

     REAL(wp) :: stab_fun, hz0_fac

     !$ACC ROUTINE SEQ

     IF(RIB.GE.0._wp)THEN
       !Cosmo
       !stab_fun = 1._wp / ( 1._wp + 10._wp*RIB/SQRT(1._wp+5*RIB) )

       !H&B
       stab_fun = 1._wp / ( 1._wp + 10._wp*RIB*(1._wp+8._wp*RIB) )
    ELSE
       hz0_fac = ( max(hz0, 1._wp)**(1._wp/3._wp) - 1._wp )**1.5_wp ! FLO - hz0 can be < 1 then, the **1.5 is invalid.
       !for water surface (z0/h)**(1/3)<<1 giving hz0_fac=SQRT(h/z0)
       !Generally it is explicitly written for water surface but i don't
       !see any reason to do that.
       stab_fun = 1._wp + 10._wp*ABS(RIB)/(1._wp + 75._wp*tc*hz0_fac*SQRT(ABS(RIB)))
     END IF

  END FUNCTION stability_function_mom
  !
  ! stability_function_heat
  !------------------------------------------------------------------------
#ifndef _OPENACC
  ELEMENTAL &
#endif
  PURE FUNCTION stability_function_heat(RIB, hzh, tc) RESULT(stab_fun)
     REAL(wp), INTENT(IN) :: RIB, hzh, tc

     REAL(wp) :: stab_fun, hzh_fac

     !$ACC ROUTINE SEQ

     IF(RIB.GE.0._wp)THEN
       !Cosmo
       !stab_fun = 1._wp / ( 1._wp + 15._wp*RIB*SQRT(1._wp+5*RIB) )

       !H&B
       stab_fun = 1._wp / ( 1._wp + 10._wp*RIB*(1._wp+8._wp*RIB) )
     ELSE
       hzh_fac = ( max(hzh, 1._wp)**(1._wp/3._wp) - 1._wp )**1.5_wp
       !for water surface (zh/h)**(1/3)<<1 giving hzh_fac=SQRT(h/zh)
       !Generally it is explicitly written for water surface but i don't
       !see any reason to do that.
       stab_fun = 1._wp + 15._wp*ABS(RIB)/(1._wp + 75._wp*tc*hzh_fac*SQRT(ABS(RIB)))
     END IF
  END FUNCTION stability_function_heat

  !
  ! factor_mom
  !------------------------------------------------------------------------
  ! Businger Dyer similarity profile:
  ! Louis (1979) A Parametirc model of vertical eddy fluxes in the atmosphere
  ! and R. B. Stull's book
  !------------------------------------------------------------------------
#ifndef _OPENACC
  ELEMENTAL &
#endif
  PURE FUNCTION businger_mom(z0, z1, L) RESULT(factor)
    ! Use VALUE attribute in order to prevent 'Reference argument passing prevents parallelization' message by nvhpc
    REAL(wp), VALUE, INTENT(IN) :: z0, z1, L
    REAL(wp) :: factor, zeta, psi, lamda
    REAL(wp) :: zeta0, psi0, lamda0

    !$ACC ROUTINE SEQ

    IF(L > 0._wp)THEN !Stable
      zeta  = z1/L
      zeta0 = z0/L
      IF(zeta > 1._wp)THEN !Zeng etal 1997 J. Clim
        psi    = -bsm*LOG(zeta) - zeta + 1
        psi0   = -bsm*LOG(zeta0) - zeta0 + 1
        factor = (LOG(L/z0) + bsh - psi + psi0  ) / ckap
      ELSE
        psi  = -bsm*zeta
        psi0 = -bsm*zeta0
        factor = ( LOG(z1/z0) - psi + psi0 ) / ckap
      END IF
    ELSEIF(L < 0._wp)THEN !unstable
      zeta   = z1/L
      zeta0  = z0/L
      lamda  = SQRT(SQRT(1._wp - bum*zeta))
      lamda0 = SQRT(SQRT(1._wp - bum*zeta0))

      psi    = 2._wp * LOG(1._wp+lamda) + LOG(1._wp+lamda*lamda) - &
              2._wp * ATAN(lamda) + pi_2 - 3._wp*ln2

      psi0   = 2._wp * LOG(1._wp+lamda0) + LOG(1._wp+lamda0*lamda0) - &
              2._wp * ATAN(lamda0) + pi_2 - 3._wp*ln2

      factor = ( LOG(z1/z0) - psi + psi0 ) / ckap
    ELSE !neutral
      factor = LOG(z1/z0) / ckap
    END IF

  END FUNCTION businger_mom
  !
  ! factor_heat
  !------------------------------------------------------------------------
  ! Businger Dyer similarity profile:
  ! Louis (1979) A Parametirc model of vertical eddy fluxes in the atmosphere
  ! and R. B. Stull's book
  !------------------------------------------------------------------------
#ifndef _OPENACC
  ELEMENTAL &
#endif
  PURE FUNCTION businger_heat(z0, z1, L) RESULT(factor)
    ! Use VALUE attribute in order to prevent 'Reference argument passing prevents parallelization' message by nvhpc
    REAL(wp), VALUE, INTENT(IN) :: z0, z1, L
    REAL(wp) :: factor, zeta, lamda, psi
    REAL(wp) :: zeta0, lamda0, psi0

    !$ACC ROUTINE SEQ

    IF(L > 0._wp)THEN !Stable
      zeta   = z1/L
      zeta0  = z0/L
      IF(zeta > 1._wp)THEN !Zeng etal 1997 J. Clim
        psi    = -bsh*LOG(zeta) - zeta + 1
        psi0   = -bsh*LOG(zeta0) - zeta0 + 1
        factor = (LOG(L/z0) + bsh - psi + psi0  ) / ckap
      ELSE
        psi    = -bsh*zeta
        psi0   = -bsh*zeta0
        factor = (LOG(z1/z0) - psi + psi0) / ckap
      END IF
    ELSEIF(L < 0._wp)THEN !unstable
      zeta   = z1/L
      zeta0  = z0/L
      lamda  = SQRT(1._wp - buh*zeta)
      lamda0 = SQRT(1._wp - buh*zeta0)
      psi    = 2._wp * ( LOG(1._wp+lamda) - ln2 )
      psi0   = 2._wp * ( LOG(1._wp+lamda0) - ln2 )
      factor = (LOG(z1/z0) - psi + psi0) / ckap
    ELSE !Neutral
      factor = LOG(z1/z0) / ckap
    END IF

  END FUNCTION businger_heat

END MODULE mo_vdf_diag_smag
