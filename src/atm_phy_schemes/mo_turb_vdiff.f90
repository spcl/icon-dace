!
! VDIFF turbulent mixing scheme
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





MODULE mo_turb_vdiff

  USE mo_exception,          ONLY: finish
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c
  USE mo_impl_constants,     ONLY: min_rlcell_int, SUCCESS
  USE mo_kind,               ONLY: wp
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_model_domain,       ONLY: t_patch
  USE mo_nh_testcases_nml,   ONLY: isrfc_type, shflx, lhflx
  USE mo_physical_constants, ONLY: rgrav, cpd
  USE mo_turb_vdiff_sma,     ONLY: atm_exchange_coeff3d, diffuse_hori_velocity, &
                                 & diffuse_vert_velocity,diffuse_scalar
  USE mo_turb_vdiff_config,  ONLY: t_vdiff_config
  USE mo_turb_vdiff_params,  ONLY: VDIFF_TURB_3DSMAGORINSKY, VDIFF_TURB_TTE, cchar, totte_min, &
                                 & tpfac1, tpfac2, tpfac3
  USE mo_turb_vdiff_diag,    ONLY: atm_exchange_coeff, sfc_exchange_coeff

  IMPLICIT NONE
  PRIVATE

  !-------------------
  ! Module parameters
  !-------------------

  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_turb_vdiff'

  !------------------
  ! Module variables
  !------------------

  ! general
  INTEGER :: vdiff_initialized_cnt = 0 !< Number of times vdiff_init was called.
  INTEGER :: initial_khydromet !< Initial init call: number of hydrometeors.
  INTEGER :: initial_ntracers !< Initial init call: number of additional tracers.

  LOGICAL :: ljsb !< Use JSBACH LSS.

  ! solver
  INTEGER :: nvar_vdiff    !< total number of variables affected by turbulent mixing
  INTEGER :: iu !< Right-hand-side index of zonal wind.
  INTEGER :: iv !< Right-hand-side index of meridional wind.
  INTEGER :: ih !< Right-hand-side index of heat / dry static energy.
  INTEGER :: iqv !< Right-hand-side index of water vapor.
  INTEGER :: ixl !< Right-hand-side index of liquid water.
  INTEGER :: ixi !< Right-hand-side index of ice.
  INTEGER :: ixv !< Right-hand-side index of water variance.
  INTEGER :: itotte !< Right-hand-side index of total turbulent energy.
  INTEGER :: ithv !< Right-hand-side index of virtual temperature variance.
  INTEGER :: itrc_start !< Start index in right-hand side of additional tracers.
  INTEGER :: nmatrix !< Number of tridiagonal matrices produced by `vdiff_down`.
  INTEGER :: imh !< Matrix index for heat / dry static energy.
  INTEGER :: imqv !< Matrix index for water vapor.
  INTEGER :: imuv !< Matrix index for wind.

  !> Index of the matrix for the given quantity. Index with iu,iv,ih,... (nvar_vdiff).
  INTEGER, ALLOCATABLE :: matrix_idx(:)
  !> Offset from klev of the bottom level for the given quantity. Index with iu,iv,ih,... (nvar_vdiff).
  INTEGER, ALLOCATABLE :: ibtmoffset_var  (:)
  !> Offset from klev for the given matrix index (nmatrix).
  INTEGER, ALLOCATABLE :: ibtmoffset_mtrx (:)

  !------------------
  ! Module interface
  !------------------

  PUBLIC :: vdiff_init
  PUBLIC :: vdiff_cleanup
  PUBLIC :: vdiff_down
  PUBLIC :: vdiff_up

  PUBLIC :: matrix_to_richtmyer_coeff
  PUBLIC :: vdiff_update_boundary
  PUBLIC :: vdiff_surface_flux
  PUBLIC :: vdiff_mixed_time_value
  PUBLIC :: vdiff_new_time_value
  PUBLIC :: vdiff_get_richtmyer_coeff_momentum
  PUBLIC :: vdiff_get_tke

  PUBLIC    :: ih
  PROTECTED :: ih
  PUBLIC    :: imh
  PROTECTED :: imh
  PUBLIC    :: imqv
  PROTECTED :: imqv
  PUBLIC    :: imuv
  PROTECTED :: imuv
  PUBLIC    :: iqv
  PROTECTED :: iqv
  PUBLIC    :: iu
  PROTECTED :: iu
  PUBLIC    :: iv
  PROTECTED :: iv
  PUBLIC    :: ixl
  PROTECTED :: ixl
  PUBLIC    :: nmatrix
  PROTECTED :: nmatrix
  PUBLIC    :: nvar_vdiff
  PROTECTED :: nvar_vdiff

CONTAINS

  SUBROUTINE vdiff_init( khydromet, ktrac )

    INTEGER, INTENT(IN) :: khydromet !< Number of hydrometeors to consider.
    INTEGER, INTENT(IN) :: ktrac !< Number of additional tracers.

    ! FIXME: hard-coded JSBACH use
    ljsb = .TRUE.

    IF (vdiff_initialized_cnt > 0) THEN
      IF (initial_khydromet /= khydromet .OR. initial_ntracers /= ktrac) THEN

        CALL finish(thismodule//':vdiff_init', 'vdiff_init called twice with different arguments')

      END IF

      vdiff_initialized_cnt = vdiff_initialized_cnt + 1
    ELSE
      vdiff_initialized_cnt = 1
      initial_khydromet = khydromet
      initial_ntracers = ktrac

      CALL init_vdiff_solver( khydromet, ktrac )
    END IF

  END SUBROUTINE vdiff_init
  !-------------

  SUBROUTINE vdiff_cleanup
    vdiff_initialized_cnt = vdiff_initialized_cnt - 1

    IF (vdiff_initialized_cnt == 0) THEN
      CALL cleanup_vdiff_solver
    ELSE IF (vdiff_initialized_cnt < 0) THEN
      CALL finish(thismodule//':vdiff_cleanup', 'vdiff_cleanup called more often than vdiff_init')
    END IF
  END SUBROUTINE vdiff_cleanup

  !>
  !! First half of the driver routine for turbulent mixing.
  !!
  SUBROUTINE vdiff_down( kbdim, nblks_c, nblks_v, nblks_e,              &! in
                       & klev, klevm1, klevp1, ktrac,                   &! in
                       & ksfc_type, idx_wtr, idx_ice, idx_lnd,          &! in
                       & pdtime,  pcoriol,                              &! in
                       & patch,                                         &! in
                       & l2moment,                                      &! in
                       & pzf, pzh,                                      &! in
                       & pfrc,                                          &! in
                       & ptsfc_tile, pocu,      pocv,       ppsfc,      &! in
                       & pum1,       pvm1,      pwm1,                   &! in
                       & ptm1,       pqm1,                              &! in
                       & pxlm1,      pxim1,     pxm1,       pxtm1,      &! in
                       & pmair,      rho,                               &! in
                       & paphm1,     papm1,                             &! in
                       & ptvm1,      paclc,     pxt_emis,   pthvvar,    &! in
                       & pxvar,      pz0m_tile,                         &! in
                       & ptottem1,                                      &! in
                       & vdiff_config,                                  &! in
                       & pustar,     pwstar,    pwstar_tile,            &! inout, out, inout
                       & pqsat_tile, phdtcbl,                           &! out
                       & pri,        pri_tile,  pmixlen,                &! out
                       & pcfm,       pcfm_tile, pcfh,       pcfh_tile,  &! out
                       & pcfv,       pcftotte,  pcfthv,                 &! out
                       & aa,         aa_btm,    bb,         bb_btm,     &! out
                       & ddt_u, ddt_v, ddt_w,                           &! out
                       & ta_hori_tend,                                  &! out
                       & qv_hori_tend,                                  &! out
                       & ql_hori_tend,                                  &! out
                       & qi_hori_tend,                                  &! out
                       & qnc_hori_tend,                                 &! out
                       & qni_hori_tend,                                 &! out
                       & pfactor_sfc, pcpt_tile,                        &! out
                       & pcptgz,                                        &! out
                       & pzthvvar,   pztottevn,                         &! out
                       & pch_tile,                                      &! out, for "nsurf_diag"
                       & pbn_tile,   pbhn_tile,                         &! out, for "nsurf_diag"
                       & pbm_tile,   pbh_tile,                          &! out, for "nsurf_diag"
                       & pcsat,                                         &! in
                       & pcair,                                         &! in
                       & paz0lh)                                         ! in

    INTEGER, INTENT(IN) :: kbdim, klev, klevm1, klevp1, ktrac, nblks_c, nblks_v, nblks_e
    INTEGER, INTENT(IN) :: ksfc_type, idx_wtr, idx_ice, idx_lnd
    LOGICAL, INTENT(IN),OPTIONAL :: l2moment
    REAL(wp),INTENT(IN) :: pdtime

    TYPE(t_patch)   ,TARGET ,INTENT(IN) :: patch

    REAL(wp),INTENT(IN) ::          &
      & pcoriol   (:,:)   ,&!< (kbdim) Coriolis parameter: 2*omega*sin(lat)
      & pzf       (:,:,:) ,&!< (kbdim,klev) geopotential height above sea level, full level
      & pzh       (:,:,:) ,&!< (kbdim,klevp1) geopotential height above sea level, half level
      & pfrc      (:,:,:) ,&!< (kbdim,ksfc_type) area fraction of each surface type
      & ptsfc_tile(:,:,:) ,&!< (kbdim,ksfc_type) surface temperature
      & pocu      (:,:)   ,&!< (kbdim) eastward  velocity of ocean sfc current
      & pocv      (:,:)   ,&!< (kbdim) northward velocity of ocean sfc current
      & ppsfc     (:,:)     !< (kbdim) surface pressure

    REAL(wp),INTENT(IN) ::        &
      & ptm1    (:,:,:)   ,&!< (kbdim,klev) temperature at step t-dt
      & pqm1    (:,:,:)   ,&!< (kbdim,klev) specific humidity at step t-dt
      & pxlm1   (:,:,:)   ,&!< (kbdim,klev) cloud water concentration at step t-dt
      & pxim1   (:,:,:)   ,&!< (kbdim,klev) cloud ice   concentration at step t-dt
      & pxm1    (:,:,:)   ,&!< (kbdim,klev) cloud water + cloud ice at step t-dt
      & pxtm1   (:,:,:,:)   !< (kbdim,klev,ktrac) specific density of other tracers at t-dt

    REAL(wp),INTENT(IN) ::        &
      & pmair   (:,:,:)   ,&!< (kbdim,klev)     air mass [kg/m2]
      & paphm1  (:,:,:)   ,&!< (kbdim,klevp1) half level pressure [Pa]
      & papm1   (:,:,:)   ,&!< (kbdim,klev) full level pressure [Pa]
      & ptvm1   (:,:,:)   ,&!< (kbdim,klev) virtual temperature
      & paclc   (:,:,:)   ,&!< (kbdim,klev) cloud fraction
      & pxt_emis(:,:,:)     !< (kbdim,ktrac) tracer tendency due to surface emission
                          !< and dry deposition

    REAL(wp),INTENT(IN) ::         &
      & pthvvar  (:,:,:)  ,&!< (kbdim,klev) variance of virtual pot. temp. at step t-dt
      & pxvar    (:,:,:)  ,&!< (kbdim,klev) step t-dt
      & pz0m_tile(:,:,:)    !< (kbdim,ksfc_type) roughness length at step t-dt

    REAL(wp),INTENT(IN) ::        &
      & rho     (:,:,:)   ,&!< (kbdim,klev) air density [kg/m3]
      & pum1    (:,:,:)   ,&!< (kbdim,klev) u-wind at step t-dt
      & pvm1    (:,:,:)     !< (kbdim,klev) q-wind at step t-dt

    REAL(wp),INTENT(IN)  :: ptottem1(:,:,:)    !< (kbdim,klev) TTE at step t-dt

    TYPE(t_vdiff_config), INTENT(IN) :: vdiff_config !< tuning parameters

    ! Grid-box mean friction velocity.
    ! In: value at step t-2dt computed in the previous time step,
    ! used in the computation of PBL height (then mixing length);
    ! Out: computed in sfc_exchange_coeff at step t-dt.

    REAL(wp),INTENT(INOUT) :: pustar (:,:)         !< (kbdim)
    REAL(wp),INTENT(OUT)   :: pwstar (:,:)         !< (kbdim)
    REAL(wp),INTENT(INOUT) :: pwstar_tile(:,:,:)   !< (kbdim,ksfc_type)
    REAL(wp),INTENT(IN)    :: pwm1    (:,:,:)      !< (kbdim,klevp1) vertical wind in m/s
    REAL(wp),INTENT(OUT)   :: ddt_u (:,:,:),      &
                            & ddt_v (:,:,:),      &
                            & ddt_w (:,:,:)
    REAL(wp),INTENT(OUT)   :: ta_hori_tend (:,:,:)
    REAL(wp),INTENT(OUT)   :: qv_hori_tend (:,:,:)
    REAL(wp),INTENT(OUT)   :: ql_hori_tend (:,:,:)
    REAL(wp),INTENT(OUT)   :: qi_hori_tend (:,:,:)
    ! Only needed for 2 moment scheme
    REAL(wp),INTENT(OUT),OPTIONAL :: qnc_hori_tend (:,:,:)
    REAL(wp),INTENT(OUT),OPTIONAL :: qni_hori_tend (:,:,:)

    ! Variables with intent(out)

    REAL(wp),INTENT(OUT) :: pqsat_tile(:,:,:)   !< (kbdim,ksfc_type) saturation specific
                                              !< humidity at sfc.
                                              !< (step t-dt)

    REAL(wp),INTENT(OUT) :: phdtcbl(:,:)    !< (kbdim) height of the top of the atmospheric dry
                                          !< convective boundary layer

    REAL(wp),INTENT(OUT) ::      &   ! out
      & pri      (:,:,:)  ,&!< (kbdim,klev) Richardson number
      & pri_tile (:,:,:)  ,&!< (kbdim,ksfc_type) Richardson number
      & pmixlen  (:,:,:)  ,&!< (kbdim,klev) mixing length
      & pcfm     (:,:,:)  ,&!< (kbdim,klev) exchange coeff. for u, v
      & pcfm_tile(:,:,:)  ,&!< (kbdim,ksfc_type) exchange coeff. for u, v
      & pcfh     (:,:,:)  ,&!< (kbdim,klev) exchange coeff. for heat and tracers
      & pcfh_tile(:,:,:)  ,&!< (kbdim,ksfc_type) exchange coeff. for heat and tracers
      & pcfv     (:,:,:)  ,&!< (kbdim,klev) exchange coeff. for variance of qx
      & pcftotte (:,:,:)  ,&!< (kbdim,klev) exchange coeff. for TTE
      & pcfthv   (:,:,:)    !< (kbdim,klev) exchange coeff. for variance of theta_v

    ! Coefficient matrices and right-hand-side vectors.
    ! _btm refers to the lowest model level (i.e., full level "klev", not the surface)

    REAL(wp),INTENT(OUT) ::           &  ! out
      & aa     (:,:,:,:,:)    ,&!< (kbdim,klev,3,nmatrix) coeff. matrices, all variables
      & aa_btm (:,:,:,imh:,:) ,&!< (kbdim,3,ksfc_type,imh:imqv) last row of coeff. matrix of heat and moisture
      & bb     (:,:,:,:)      ,&!< (kbdim,klev,nvar_vdiff) r.h.s., all variables
      & bb_btm (:,:,ih:,:)      !< (kbdim,ksfc_type,ih:iqv) last row of r.h.s. of heat and moisture

    ! Other variables to be passed on to the second part of turbulence solver

    REAL(wp),INTENT(OUT) ::       &  ! out
      & pfactor_sfc(:,:)  ,&!< (kbdim) prefactor for the exchange coeff.
      & pcpt_tile (:,:,:) ,&!< (kbdim,ksfc_type) dry static energy at surface
      & pcptgz    (:,:,:) ,&!< (kbdim,klev) dry static energy
      & pzthvvar  (:,:,:) ,&!< (kbdim,klev)
      & pztottevn (:,:,:)   !< (kbdim,klev) intermediate value of TTE
    REAL(wp) :: jztottevn(kbdim,nblks_c)

    REAL(wp), INTENT(OUT) :: pch_tile(:,:,:)    !< (kbdim,ksfc_type) out, for "nsurf_diag"
    REAL(wp), INTENT(OUT) :: pbn_tile(:,:,:)    !< (kbdim,ksfc_type) out, for "nsurf_diag"
    REAL(wp), INTENT(OUT) :: pbhn_tile(:,:,:)   !< (kbdim,ksfc_type) out, for "nsurf_diag"
    REAL(wp), INTENT(OUT) :: pbm_tile(:,:,:)    !< (kbdim,ksfc_type) out, for "nsurf_diag"
    REAL(wp), INTENT(OUT) :: pbh_tile(:,:,:)    !< (kbdim,ksfc_type) out, for "nsurf_diag"

    REAL(wp), OPTIONAL, INTENT(IN) ::          &
      & pcsat     (:,:)          ,&!< (kbdim) area fraction with wet land surface
      & pcair     (:,:)          ,&!< (kbdim) area fraction with wet land surface
      & paz0lh    (:,:)            !< (kbdim) surface roughness length over land for heat

    ! Local variables

    REAL(wp) :: zghf   (kbdim,klev,nblks_c)   !< geopotential height above ground, full level
    REAL(wp) :: zghh   (kbdim,klevp1,nblks_c) !< geopotential height above ground, full level

    REAL(wp) :: zfactor(kbdim,klev,nblks_c)   !< prefactor for the exchange coefficients
    REAL(wp) :: zrmairm(kbdim,klev,nblks_c)
    REAL(wp) :: zrmairh(kbdim,klevm1,nblks_c)

    REAL(wp), DIMENSION(kbdim,klev,nblks_c)   :: km_c
    REAL(wp), DIMENSION(kbdim,klevp1,nblks_v) :: km_iv
    REAL(wp), DIMENSION(kbdim,klevp1,nblks_e) :: km_ie
    REAL(wp), DIMENSION(kbdim,klevp1,nblks_c) :: kh_ic, km_ic
    REAL(wp), DIMENSION(kbdim,klev,nblks_e)   :: vn

    REAL(wp) :: u_vert(kbdim,klev,nblks_v), v_vert(kbdim,klev,nblks_v)
    REAL(wp) :: div_c(kbdim,klev,nblks_c)
    REAL(wp) :: rho_ic(kbdim,klevp1,nblks_c)
    REAL(wp) :: w_vert(kbdim,klevp1,nblks_v)
    REAL(wp) :: w_ie(kbdim,klevp1,nblks_e)

    ! _b denotes value at the bottom level (the klev-th full level)

    REAL(wp) :: ztheta_b (kbdim,nblks_c)  !< potential temperature
    REAL(wp) :: zthetav_b(kbdim,nblks_c)  !< virtual potential temperature
    REAL(wp) :: zthetal_b(kbdim,nblks_c)  !< liquid (and ice?) pot. temp.
    REAL(wp) :: zqsat_b  (kbdim,nblks_c)  !< specific humidity at saturation
    REAL(wp) :: zlh_b    (kbdim,nblks_c)  !< latent heat

    REAL(wp) :: zconst

    INTEGER  :: jl, jk
    INTEGER, PARAMETER :: tracer_dry_static = 1
    INTEGER, PARAMETER :: tracer_water = 2

    ! OMP variables
    INTEGER                             :: jb,jbs,jbe,jcs,jce,ncd,rls,rle

    !---- Local variables
    !$ACC DATA &
    !$ACC   CREATE(zghf, zghh, zfactor, zrmairm, zrmairh, jztottevn) &
    !$ACC   CREATE(ztheta_b, zthetav_b, zthetal_b, zqsat_b, zlh_b)


    rls = grf_bdywidth_c+1
    rle = min_rlcell_int

    ncd = MAX(1,patch%n_childdom)
    jbs = patch%cells%start_blk(rls,  1)
    jbe = patch%cells%end_blk  (rle,ncd)

    !##############################################################################
    !## jb loop1
    !$OMP PARALLEL DO PRIVATE(jb,jcs,jce)
    DO jb = jbs, jbe
      CALL get_indices_c(patch, jb, jbs, jbe, jcs, jce, rls, rle)
      IF (jcs>jce) CYCLE
    !##############################################################################
      !----------------------------------------------------------------------
      ! 0. Compute useful local fields
      !----------------------------------------------------------------------

      ! geopotential height above ground

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO jk = 1,klev
        DO jl = jcs,jce
          zghf(jl,jk,jb) = pzf(jl,jk,jb) - pzh(jl,klevp1,jb)
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO jk = 1,klevp1
        DO jl = jcs,jce
          zghh (jl,jk,jb) = pzh(jl,jk,jb) - pzh(jl,klevp1,jb)
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      ! reciprocal layer mass
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO jk = 1,klev
        DO jl = jcs,jce
          zrmairm(jl,jk,jb) = 1._wp / pmair(jl,jk,jb)
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO jk = 1,klevm1
        DO jl = jcs,jce
          zrmairh(jl,jk,jb) = 2._wp / (pmair(jl,jk,jb) + pmair(jl,jk+1,jb))
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO jk = 1,klev
        DO jl = jcs,jce
          ddt_u(jl,jk,jb) = 0._wp
          ddt_v(jl,jk,jb) = 0._wp
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO jk = 1,klev+1
         DO jl = jcs,jce
           ddt_w(jl,jk,jb) = 0._wp
         END DO
      END DO
      !$ACC END PARALLEL LOOP

    !##############################################################################
    !## jb end loop1
    END DO
    !$OMP END PARALLEL DO
    !##############################################################################

    SELECT CASE ( vdiff_config%turb ) ! select turbulent scheme
    CASE ( VDIFF_TURB_TTE )
      !----------------------------------------------------------------------
      ! 1. Compute various thermodynamic variables; Diagnose PBL extension;
      !    Compute exchange coefficients for half levels [1+1/2, klev-1/2];
      !    Get TTE and variance of theta_v at intermediate time step.
      !----------------------------------------------------------------------

      rls = grf_bdywidth_c+1
      rle = min_rlcell_int

      ncd = MAX(1,patch%n_childdom)
      jbs = patch%cells%start_blk(rls,  1)
      jbe = patch%cells%end_blk  (rle,ncd)
      !##############################################################################
      !## jb loop2
      !$OMP PARALLEL DO PRIVATE(jb,jcs,jce)
      DO jb = jbs, jbe
        CALL get_indices_c(patch, jb, jbs, jbe, jcs, jce, rls, rle)
        IF (jcs>jce) CYCLE
      !##############################################################################

        ! DA: this routine is async aware, so it's safe not not wait here
        CALL atm_exchange_coeff( jb,                                                         &! in, for debugging only
                              & jcs, jce, kbdim, klev, klevm1,                              &! in
                              & pdtime, pcoriol(:,jb),                                      &! in
                              & zghf(:,:,jb), zghh(:,:,jb),                                 &! in
                              & pum1(:,:,jb), pvm1(:,:,jb), ptm1(:,:,jb), ptvm1(:,:,jb),    &! in
                              & pqm1(:,:,jb), pxm1(:,:,jb),                                 &! in
                              & papm1(:,:,jb), paphm1(:,:,jb), paclc(:,:,jb), pustar(:,jb), &! in
                              & pthvvar(:,:,jb), ptottem1(:,:,jb),                          &! in
                              & vdiff_config,                                           &! in
                              & pcptgz(:,:,jb), phdtcbl(:,jb),                              &! out
                              & pzthvvar(:,1:klevm1,jb),                                    &! out
                              & pztottevn(:,1:klevm1,jb),                                   &! out
                              & pcfm    (:,1:klevm1,jb), pcfh  (:,1:klevm1,jb),             &! out
                              & pcfv    (:,1:klevm1,jb),                                    &! out
                              & pcftotte(:,1:klevm1,jb),                                    &! out
                              & pcfthv  (:,1:klevm1,jb),zfactor(:,1:klevm1,jb),             &! out
                              & ztheta_b(:,jb), zthetav_b(:,jb), zthetal_b(:,jb),           &! out, for "sfc_exchange_coeff"
                              & zqsat_b(:,jb),  zlh_b(:,jb),                                &! out, for "sfc_exchange_coeff"
                              & pri(:,1:klevm1,jb), pmixlen(:,1:klevm1,jb)                  )! out, for output

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
        DO jl = jcs,jce
          pmixlen(jl,klev,jb) = -999._wp
        END DO
        !$ACC END PARALLEL LOOP

        !-----------------------------------------------------------------------
        ! 2. Compute exchange coefficients at the air-sea/ice/land interface.
        !    Get boundary condition for TTE and variance of theta_v.
        !-----------------------------------------------------------------------

        ! DA: this routine is async, no need to wait
        CALL sfc_exchange_coeff( jb, jcs, jce, kbdim, ksfc_type, patch,        &! in
                              & idx_wtr, idx_ice, idx_lnd,                     &! in
                              & pz0m_tile(:,jb,:),  ptsfc_tile(:,jb,:),        &! in
                              & pfrc(:,jb,:),       phdtcbl(:,jb),             &! in
                              & pocu(:,jb),         pocv(:,jb),   ppsfc(:,jb), &! in
                              & zghf(:,klev,jb),                               &! in
                              & pum1(:,klev,jb),    pvm1  (:,klev,jb),         &! in
                              & ptm1(:,klev,jb),                               &! in
                              & pqm1(:,klev,jb),    pxm1  (:,klev,jb),         &! in
                              & zqsat_b  (:,jb),    zlh_b    (:,jb),           &! in
                              & ztheta_b (:,jb),    zthetav_b(:,jb),           &! in
                              & zthetal_b(:,jb),    paclc (:,klev,jb),         &! in
                              & ptottem1(:,klev,jb),pzthvvar(:,klevm1,jb),     &! in
                              & vdiff_config,                                  &! in
                              & pwstar(:,jb),       pwstar_tile(:,jb,:),       &! out, inout
                              & pqsat_tile(:,jb,:), pcpt_tile(:,jb,:),         &! out
                              & pri    (:,klev,jb), pri_tile(:,jb,:),          &! out
                              & pcfm   (:,klev,jb), pcfm_tile(:,jb,:),         &! out
                              & pcfh   (:,klev,jb), pcfh_tile(:,jb,:),         &! out
                              & pcfv   (:,klev,jb),                            &! out
                              & pcftotte(:,klev,jb),pcfthv  (:,klev,jb),       &! out
                              & zfactor(:,klev,jb),                            &! out
                              & pztottevn(:,klev,jb),                          &! out
                              & pzthvvar(:,klev,jb),                           &! out
                              & jztottevn(:,jb),                               &! out
                              & pustar(:,jb),                                  &! out, for "atm_exchange_coeff" at next time step
                              & pch_tile(:,jb,:),                              &! out, for "nsurf_diag"
                              & pbn_tile(:,jb,:),   pbhn_tile(:,jb,:),         &! out, for "nsurf_diag"
                              & pbm_tile(:,jb,:),   pbh_tile(:,jb,:),          &! out, for "nsurf_diag"
                              & paz0lh(:,jb),                                  &! in, optional
                              & pcsat(:,jb),                                   &! in, optional
                              & pcair(:,jb))                                    ! in, optional

        IF ( isrfc_type == 1 ) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
          DO jl = jcs,jce
            pztottevn(jl,klev,jb) = jztottevn(jl,jb)
          END DO
          !$ACC END PARALLEL LOOP
        END IF

      !##############################################################################
      !## jb end loop2
      END DO
      !$OMP END PARALLEL DO
      !##############################################################################

    CASE ( VDIFF_TURB_3DSMAGORINSKY )

      !$ACC DATA &
      !$ACC   CREATE(kh_ic, km_ic, km_c, km_iv, km_ie, vn, u_vert, v_vert, w_vert, rho_ic, div_c, w_ie)
      CALL atm_exchange_coeff3d ( kbdim, nblks_c,                                     &! in
                            & klev, klevm1, klevp1,                                   &! in
                            & ksfc_type, idx_lnd,                                     &! in
                            & patch,                                                  &! in
                            & pz0m_tile(:,:,:), ptsfc_tile(:,:,:), pfrc(:,:,:),       &! in
                            & ppsfc(:,:),                                             &! in
                            & zghf(:,:,:),                                            &! in
                            & pum1(:,:,:), pvm1(:,:,:), pwm1(:,:,:),                  &! in
                            & ptm1(:,:,:), ptvm1(:,:,:),                              &! in
                            & pqm1(:,:,:), pxm1(:,:,:),                               &! in
                            & rho(:,:,:),                                             &! in
                            & papm1(:,:,:), paphm1(:,:,:),                            &! in
                            & vdiff_config,                                           &! in
                            & pri_tile(:,:,:),                                        &! out
                            & pcfm_tile(:,:,:),                                       &! out
                            & pcfh_tile(:,:,:),                                       &! out
                            & pqsat_tile(:,:,:), pcpt_tile(:,:,:),                    &! out
                            & pcptgz(:,:,:),                                          &! out
                            & pzthvvar(:,:,:),                                        &! out
                            & pztottevn(:,:,:), pmixlen(:,:,:),                       &! out
                            & pcfm    (:,:,:), pcfh  (:,:,:),                         &! out
                            & pcfv    (:,:,:),                                        &! out
                            & pcftotte(:,:,:),                                        &! out
                            & pcfthv  (:,:,:),                                        &! out
                            & km_c(:,:,:), km_iv(:,:,:), km_ie(:,:,:),                &! out
                            & km_ic(:,:,:), kh_ic(:,:,:),                             &! out,
                            & zfactor(:,:,:),                                         &! out
                            & u_vert(:,:,:), v_vert(:,:,:), div_c(:,:,:),             &! out, for "sfc_exchange_coeff"
                            & rho_ic(:,:,:), w_vert(:,:,:), w_ie(:,:,:),              &! out
                            & vn(:,:,:),                                              &! out,
                            & pch_tile(:,:,:),                                        &! out, for "nsurf_diag"
                            & pbn_tile(:,:,:), pbhn_tile(:,:,:),                      &! out
                            & pbm_tile(:,:,:), pbh_tile(:,:,:),                       &! out
                            & pcsat=pcsat(:,:), pcair=pcair(:,:)                      )! in, optional


      CALL diffuse_hori_velocity( kbdim,                                            &
                                & patch,                                            &
                                & km_c(:,:,:), km_iv(:,:,:),                        &
                                & u_vert(:,:,:), v_vert(:,:,:), div_c(:,:,:),       &
                                & rho(:,:,:), vn(:,:,:),                            &
                                & ddt_u(:,:,:), ddt_v(:,:,:) )

      CALL diffuse_vert_velocity( kbdim,                                            &
                                & patch,                                            &
                                & rho_ic(:,:,:), w_vert(:,:,:), w_ie(:,:,:),        &
                                & km_c(:,:,:), km_iv(:,:,:), km_ic(:,:,:),          &
                                & u_vert(:,:,:), v_vert(:,:,:), div_c(:,:,:),       &
                                & pum1(:,:,:), pvm1(:,:,:), pwm1(:,:,:), vn(:,:,:), &
                                & ddt_w(:,:,:), pdtime)


      ! dry static energy
      call diffuse_scalar( kbdim, ptm1(:,:,:),               &
                        & patch,                            &
                        & km_ie(:,:,:),                     &
                        & ta_hori_tend(:,:,:),              &
                        & rho,                              &
                        & tracer_dry_static,                &
                        & vdiff_config%rturb_prandtl)

      call diffuse_scalar( kbdim, pqm1(:,:,:),               &
                        & patch,                            &
                        & km_ie(:,:,:),                     &
                        & qv_hori_tend(:,:,:),              &
                        & rho,                              &
                        & tracer_water,                     &
                        & vdiff_config%rturb_prandtl)

      call diffuse_scalar( kbdim, pxlm1(:,:,:),              &
                        & patch,                            &
                        & km_ie(:,:,:),                     &
                        & ql_hori_tend(:,:,:),              &
                        & rho,                              &
                        & tracer_water,                     &
                        & vdiff_config%rturb_prandtl)

      call diffuse_scalar( kbdim, pxim1(:,:,:),              &
                        & patch,                            &
                        & km_ie(:,:,:),                     &
                        & qi_hori_tend(:,:,:),              &
                        & rho,                              &
                        & tracer_water,                     &
                        & vdiff_config%rturb_prandtl)

      IF (PRESENT(l2moment)) THEN
        IF (l2moment) THEN
          CALL diffuse_scalar( kbdim, pxtm1(:,:,:,1),          & !pxtm1(:,:,:,1) = qtrc_phy (:,:,:,iqnc)
                             & patch,                          &
                             & km_ie(:,:,:),                   &
                             & qnc_hori_tend(:,:,:),           &
                             & rho,                            &
                             & tracer_water,                   &
                             & vdiff_config%rturb_prandtl)

          CALL diffuse_scalar( kbdim, pxtm1(:,:,:,2),          & !pxtm1(:,:,:,2) = qtrc_phy (:,:,:,iqni)
                             & patch,                          &
                             & km_ie(:,:,:),                   &
                             & qni_hori_tend(:,:,:),           &
                             & rho,                            &
                             & tracer_water,                   &
                             & vdiff_config%rturb_prandtl)
        END IF
      END IF
      !$ACC WAIT
      !$ACC END DATA

    END SELECT    !select turbulent scheme

    !-----------------------------------------------------------------------
    ! 3. Set up coefficient matrix of the tri-diagonal system, then perform
    !    Gauss elimination for it. The matrix is built from
    !    - the exchange coefficients;
    !    - the prefactor "zfactor" and some additional constants ("zconst")
    !      which are determined by the spatial and temporal discretization
    !      employed for vertical diffusion;
    !    - the assumption about upper and lower boundaries, especially
    !      whether there is turbulent flux at the lower boundary for each
    !      quantity subject to turbulent mixing.
    !-----------------------------------------------------------------------

    rls = grf_bdywidth_c+1
    rle = min_rlcell_int

    ncd = MAX(1,patch%n_childdom)
    jbs = patch%cells%start_blk(rls,  1)
    jbe = patch%cells%end_blk  (rle,ncd)

    !##############################################################################
    !## jb loop3
    !$OMP PARALLEL DO PRIVATE(jb,jcs,jce)
    DO jb = jbs, jbe
      CALL get_indices_c(patch, jb, jbs, jbe, jcs, jce, rls, rle)
      IF (jcs>jce) CYCLE
    !##############################################################################

      zconst = tpfac1*pdtime

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1,klevm1
        DO jl = jcs,jce
          zfactor(jl,jk,jb) = zfactor(jl,jk,jb)*zconst
        END DO
      END DO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jl = jcs,jce
        zfactor(jl,klev,jb) = zfactor(jl, klev,jb)*zconst
      END DO
      !$ACC END PARALLEL

      CALL matrix_setup_elim( jcs, jce, kbdim, klev, klevm1, ksfc_type,           &! in
                            & pcfm     (:,:,jb),   pcfh  (:,1:klevm1,jb),         &! in
                            & pcfh_tile(:,jb,:),   pcfv  (:,:,jb),                &! in
                            & pcftotte (:,:,jb),   pcfthv(:,:,jb),                &! in
                            & zfactor  (:,:,jb),                                  &! in
                            & zrmairm(:,:,jb), zrmairh(:,:,jb),                   &! in
                            & aa(:,:,:,:,jb), aa_btm(:,:,:,:,jb)                  )! out

      ! Save for output, to be used in "update_surface"
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jl = jcs,jce
        pfactor_sfc(jl,jb) = zfactor(jl,klev,jb)
      END DO
      !$ACC END PARALLEL

      !-----------------------------------------------------------------------
      ! 4. Set up right-hand side of the tri-diagonal system and perform
      !    Gauss elimination. Factors that determine the r.h.s. include
      !    - time stepping scheme used for vertical diffusion
      !    - whether there is any other process (e.g., tracer emission)
      !      solved together with vertical diffusion.
      !-----------------------------------------------------------------------

      CALL rhs_setup( jcs, jce, kbdim, klev, klevm1,                                                  &! in
                    & ksfc_type, ktrac, pdtime,                                                       &! in
                    & pum1(:,:,jb), pvm1(:,:,jb), pcptgz(:,:,jb), pqm1(:,:,jb),                       &! in
                    & pxlm1(:,:,jb), pxim1(:,:,jb), pxvar(:,:,jb), pxtm1(:,:,jb,:), pxt_emis(:,:,jb), &! in
                    & zrmairm(:,:,jb), pztottevn(:,:,jb), pzthvvar(:,:,jb), aa(:,:,:,:,jb),           &! in
                    & bb(:,:,:,jb), bb_btm(:,:,:,jb)                                                  )! out

      CALL rhs_elim ( jcs, jce, klev,              &! in
                    & aa(:,:,:,:,jb), bb(:,:,:,jb) )! in, inout

    !##############################################################################
    !## jb end loop3
    END DO
    !$OMP END PARALLEL DO
    !##############################################################################

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE vdiff_down
  !-------------


  !> Compute surface flux of the given quantity from its surface value and Richtmyer coefficients.
  !! Flux is positive if direction is into the surface.
  !!
  !! Uses the Richtmyer coefficients to calculate the value in the lowest layer given the surface
  !! value of the requested quantity, \f$ x_a = E x_{sfc} + F \f$. As `pfactor_sfc` is \f$ \alpha
  !! \Delta t \rho_{sfc} \f$, and the flux is given by \f$ \rho C v (x_a - x_{sfc}) \f$, we have
  !! to divide the prefactor by \f$ \alpha \Delta t \f$.
  !!
  SUBROUTINE vdiff_surface_flux (jcs, kproma, dtime, pfactor_sfc, pcf, pen, pfn, x_sfc, flx)

    INTEGER, INTENT(IN) :: jcs !< Start cell index.
    INTEGER, INTENT(IN) :: kproma !< End cell index.
    REAL(wp), INTENT(IN) :: dtime !< Time step.
    REAL(wp), INTENT(IN) :: pfactor_sfc(:) !< Prefactor for exchange coefficients [kg s/m**3].
    REAL(wp), INTENT(IN) :: pcf(:) !< Exchange coefficient [m/s].
    REAL(wp), INTENT(IN) :: pen(:) !< Richtmyer E coefficient for lowest layer [1].
    REAL(wp), INTENT(IN) :: pfn(:) !< Richtmyer F coefficient for lowest layer [X/kg].
    !> Surface value of quantity for which the flux is calculated [X/kg].
    REAL(wp), INTENT(IN) :: x_sfc(:)
    REAL(wp), INTENT(INOUT) :: flx(:) !< Surface flux. Down is positive [X/m**2/s].

    INTEGER :: jc

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs, kproma
      flx(jc) = pfactor_sfc(jc) / (tpfac1 * dtime) * pcf(jc) &
          & * ((pen(jc) - 1._wp) * x_sfc(jc) + pfn(jc))
    END DO
    !$ACC END PARALLEL

  END SUBROUTINE


  !> Calculate a mixed-time value from the old and new ones.
  REAL(wp) FUNCTION vdiff_mixed_time_value (x_new, x_old) RESULT(x_hat)

    !$ACC ROUTINE SEQ

    !> Value at time `t+1`.
    REAL(wp), INTENT(IN) :: x_new
    !> Value at time `t`.
    REAL(wp), INTENT(IN) :: x_old

    x_hat = tpfac1 * x_new + (1._wp - tpfac1) * x_old

  END FUNCTION vdiff_mixed_time_value


  !> Calculate a new-time (`t+1`) value from the old and mixed-time ones.
  REAL(wp) FUNCTION vdiff_new_time_value (x_hat, x_old, scaled) RESULT(x_new)

    !$ACC ROUTINE SEQ

    !> Value at time `t+1`.
    REAL(wp), INTENT(IN) :: x_hat
    !> Value at time `t`.
    REAL(wp), INTENT(IN) :: x_old
    !> Flag indicating that `x_hat` is already divided by the mixing parameter (default: false).
    LOGICAL, OPTIONAL, INTENT(IN) :: scaled

    LOGICAL :: lscaled

    IF (PRESENT(scaled)) THEN
      lscaled = scaled
    ELSE
      lscaled = .FALSE.
    END IF

    IF (lscaled) THEN
      x_new = x_hat + (1._wp - 1._wp/tpfac1) * x_old
    ELSE
      x_new = x_hat/tpfac1 + (1._wp - 1._wp/tpfac1) * x_old
    END IF

  END FUNCTION vdiff_new_time_value


  !> Compute Richtmyer-Morton coefficients for zonal and meridional momentum.
  !! Assumes that `aa` and `bb` are in the state right after `vdiff_down`, i.e., before
  !! `vdiff_update_boundary` gets called.
  !!
  SUBROUTINE vdiff_get_richtmyer_coeff_momentum( &
        & jcs, kproma, klev, ksfc_type, aa, bb, pfactor_sfc, pcfm_tile, pmair, pen_uv, pfn_u, &
        & pfn_v &
      )

    INTEGER, INTENT(IN) :: jcs !< Start cell index.
    INTEGER, INTENT(IN) :: kproma !< End cell index.
    INTEGER, INTENT(IN) :: klev !< Number of vertical levels.
    INTEGER, INTENT(IN) :: ksfc_type !< Number of surface types.

    !> Matrices (jcs:kproma, klev, 3, nmatrix).
    REAL(wp), INTENT(IN) :: aa(:,:,:,:)
    !> Right-hand sides (jcs:kproma, klev, nvar_vdiff).
    REAL(wp), INTENT(IN) :: bb(:,:,:)
    !> Common prefactor for exchange coefficients [kg*s/m**3] (jcs:kproma).
    REAL(wp), INTENT(IN) :: pfactor_sfc(:)
    !> Exchange coefficient for momentum [m/s] (jcs:kproma,ksfc_type).
    REAL(wp), INTENT(IN) :: pcfm_tile(:,:)
    !> Layer air mass [kg/m**2] (jcs:kproma,klev)
    REAL(wp), INTENT(IN) :: pmair(:,:)

    !> Richtmyer E for u and v [1] (jcs:kproma,ksfc_type).
    REAL(wp), INTENT(INOUT) :: pen_uv(:,:)

    !> Richtmyer F for u [m/s] (jcs:kproma,ksfc_type).
    REAL(wp), INTENT(INOUT) :: pfn_u(:,:)
    !> Richtmyer F for v [m/s] (jcs:kproma,ksfc_type).
    REAL(wp), INTENT(INOUT) :: pfn_v(:,:)

    INTEGER :: jc, jsfc
    REAL(wp) :: a2, a3, a2sub, bu, bv

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR PRIVATE(a2, a3, a2sub, bu, bv) COLLAPSE(2)
    DO jsfc = 1, ksfc_type
      DO jc = jcs, kproma
        a3 = -pcfm_tile(jc, jsfc) * pfactor_sfc(jc) / pmair(jc, klev)
        a2 = 1._wp - aa(jc, klev, 1, imuv) - a3

        a2sub = a2 - aa(jc, klev, 1, imuv) * aa(jc, klev-1, 3, imuv)

        bu = bb(jc, klev, iu) - aa(jc, klev, 1, imuv) * bb(jc, klev-1, iu)
        bu = bu / a2sub

        bv = bb(jc, klev, iv) - aa(jc, klev, 1, imuv) * bb(jc, klev-1, iv)
        bv = bv / a2sub

        pen_uv(jc,jsfc) = -a3 / a2sub
        pfn_u(jc,jsfc) = bu * tpfac1
        pfn_v(jc,jsfc) = bv * tpfac1
      END DO
    END DO
    !$ACC END PARALLEL

  END SUBROUTINE vdiff_get_richtmyer_coeff_momentum


  !> Update the lower boundary of the linear system for dry static energy and humidity with the
  !! given fluxes. Optionally also incorporates momentum fluxes from the surface. The fluxes
  !! should be fluxes of mixed-time quantities \f$ \hat{x} = \alpha x^{t+1} + (1-\alpha) x^t \f$.
  !!
  !! This routine adds fluxes from the surface for dry static energy, humidity, and, optionally,
  !! momentum to the right-hand side of the linear system. Then, it performs Gaussian elimination
  !! on the bottom row of the right-hand side (the last row is considered to consist only of a
  !! single 1 on the diagonal at this point). On exit, the last row of `bb` contains the solution
  !! of the system for all variables (Gaussian elimination for the other variables has already
  !! been performed by `rhs_elim`, called by `vdiff_down`).
  !!
  !! Note When omitting the momentum fluxes, the routine implicitly (through the way the matrices
  !!   are set up initially) assumes a no-slip boundary condition `u_sfc=0`, not a zero-flux one.
  !!
  SUBROUTINE vdiff_update_boundary( &
        & jcs, kproma, klev, dtime, pmair, shflx, qflx, aa, aa_btm, s_btm, q_btm, bb, &
        & uflx, vflx &
      )

    INTEGER, INTENT(IN) :: jcs !< Start cell index.
    INTEGER, INTENT(IN) :: kproma !< End cell index.
    INTEGER, INTENT(IN) :: klev !< Number of vertical levels.
    REAL(wp), INTENT(IN) :: dtime !< Time step [s].

    !> Air mass in lowest layer [kg/m**2] (jcs:kproma).
    REAL(wp), INTENT(IN) :: pmair(:)

    !> Sensible heat flux into surface [W/m**2] (jcs:kproma).
    REAL(wp), INTENT(IN) :: shflx(:)
    !> Humidity flux into surface [kg/m**2/s] (jcs:kproma).
    REAL(wp), INTENT(IN) :: qflx(:)

    !> Matrices (jcs:kproma, klev, 3, nmatrix).
    REAL(wp), INTENT(IN) :: aa(:,:,:,:)
    !> Bottom row of matrix for heat and humidity (jcs:kproma, 3, SFC_NUM, imh:imqv).
    REAL(wp), INTENT(IN) :: aa_btm(:,:,:,imh:)
    !> Dry static energy at lowest layer [J/kg] (jcs:kproma).
    REAL(wp), INTENT(IN) :: s_btm(:)
    !> Specific humidity at lowest layer [kg/kg] (jcs:kproma).
    REAL(wp), INTENT(IN) :: q_btm(:)

    !> Right-hand sides (jcs:kproma, klev, nvar_vdiff).
    REAL(wp), INTENT(INOUT) :: bb(:,:,:)

    !> Zonal momentum flux into surface [kg*m/s/m**2/s] (jcs:kproma).
    REAL(wp), INTENT(IN), OPTIONAL :: uflx(:)
    !> Meridional momentum flux into surface [kg*m/s/m**2/s] (jcs:kproma).
    REAL(wp), INTENT(IN), OPTIONAL :: vflx(:)

    INTEGER :: jc
    REAL(wp) :: bs, bqv, bu, bv

    ! The way things are currently set up, the bottom layer of `aa` and `bb` for heat and humidity
    ! is unset. We have to construct the bottom row of the matrix and perform Gaussian elimination
    ! ourselves.

    ! The `aa_btm(:, 1, isfc, im)` element is independent of the surface type, we use the first
    ! one, whatever it is. Since we are implementing a Neumann boundary condition, the diagonal
    ! element before elimination is just `a_i,i = 1 - a_i,i-1`.

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR PRIVATE(bs, bqv)
    DO jc = jcs, kproma
      bs = s_btm(jc) / tpfac1 - dtime * shflx(jc) / pmair(jc)
      bs = bs - aa_btm(jc, 1, 1, imh) * bb(jc, klev-1, ih)
      bs = bs / (1._wp - aa_btm(jc, 1, 1, imh) - aa_btm(jc, 1, 1, imh) * aa(jc, klev-1, 3, imh))
      bb(jc, klev, ih) = bs

      bqv = q_btm(jc) / tpfac1 - dtime * qflx(jc) / pmair(jc)
      bqv = bqv - aa_btm(jc, 1, 1, imqv) * bb(jc, klev-1, iqv)
      bqv = bqv / (1._wp - aa_btm(jc,1,1,imqv) - aa_btm(jc,1,1,imqv) * aa(jc, klev-1, 3, imqv))
      bb(jc, klev, iqv) = bqv
    END DO
    !$ACC END PARALLEL

    ! For momentum, the last row is set up but not eliminated, and set up for a Dirichlet "no
    ! slip" boundary condition where the momentum at the surface is kept at zero. This is fine
    ! until ocean currents are included.
    IF (PRESENT(uflx)) THEN
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR PRIVATE(bu)
      DO jc = jcs, kproma
        bu = bb(jc, klev, iu) - dtime * uflx(jc) / pmair(jc)
        bu = bu - aa(jc, klev, 1, imuv) * bb(jc, klev-1, iu)
        bu = bu / (1._wp - aa(jc, klev, 1, imuv) - aa(jc, klev, 1, imuv) * aa(jc, klev-1, 3, imuv))
        bb(jc, klev, iu) = bu
      END DO
      !$ACC END PARALLEL
    ELSE
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR PRIVATE(bu)
      DO jc = jcs, kproma
        bu = bb(jc, klev, iu) - aa(jc, klev, 1, imuv) * bb(jc, klev-1, iu)
        bu = bu / (aa(jc, klev, 2, imuv) - aa(jc, klev, 1, imuv) * aa(jc, klev-1, 3, imuv))
        bb(jc, klev, iu) = bu
      END DO
      !$ACC END PARALLEL
    END IF

    IF (PRESENT(vflx)) THEN
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR PRIVATE(bv)
      DO jc = jcs, kproma
        bv = bb(jc, klev, iv) - dtime * vflx(jc) / pmair(jc)
        bv = bv - aa(jc, klev, 1, imuv) * bb(jc, klev-1, iv)
        bv = bv / (1._wp - aa(jc, klev, 1, imuv) - aa(jc, klev, 1, imuv) * aa(jc, klev-1, 3, imuv))
        bb(jc, klev, iv) = bv
      END DO
      !$ACC END PARALLEL
    ELSE
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR PRIVATE(bv)
      DO jc = jcs, kproma
        bv = bb(jc, klev, iv) - aa(jc, klev, 1, imuv) * bb(jc, klev-1, iv)
        bv = bv / (aa(jc, klev, 2, imuv) - aa(jc, klev, 1, imuv) * aa(jc, klev-1, 3, imuv))
        bb(jc, klev, iv) = bv
      END DO
      !$ACC END PARALLEL
    END IF

  END SUBROUTINE vdiff_update_boundary


  SUBROUTINE vdiff_get_tke (jcs, kproma, klev, vdiff_config, ptotte, pri, tke)

    INTEGER, INTENT(IN) :: jcs !< Start cell index.
    INTEGER, INTENT(IN) :: kproma !< End cell index.
    INTEGER, INTENT(IN) :: klev !< Number of vertical levels.

    !> VDIFF configuration for the current domain.
    TYPE(t_vdiff_config), INTENT(IN) :: vdiff_config

    !> Total turbulence energy (jcs:kproma,nlev) [m**2/s**2].
    REAL(wp), INTENT(IN) :: ptotte(:,:)
    !> Richardson number (jcs:kproma,nlev) [1].
    REAL(wp), INTENT(IN) :: pri(:,:)

    !> Turbulence kinetic energy as diagnosed by the TTE scheme (jcs:kproma,nlev+1) [m**2/s**2].
    REAL(wp), INTENT(INOUT) :: tke(:,:)

    INTEGER :: jc, jl
    REAL(wp) :: tte, ri
    REAL(wp) :: pr0, ek_ep_ratio_stable, ek_ep_ratio_unstable

    pr0 = vdiff_config%pr0
    ek_ep_ratio_stable = vdiff_config%ek_ep_ratio_stable
    ek_ep_ratio_unstable = vdiff_config%ek_ep_ratio_unstable

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs, kproma
      tke(jc,1) = 0._wp
    END DO

    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(tte, ri)
    DO jl = 1, klev
      DO jc = jcs, kproma
        tte = ptotte(jc,jl)
        ri = pri(jc,jl)

        IF(ri > 0._wp) THEN
          tke(jc,jl+1) = tte / (1._wp + ri/(pr0 + ek_ep_ratio_stable * ri))
        ELSE
          tke(jc,jl+1) = tte / (1._wp + ri/(ek_ep_ratio_unstable * ri - pr0))
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

  END SUBROUTINE vdiff_get_tke


  SUBROUTINE vdiff_up( jcs, kproma, kbdim, klev, klevm1, &! in
    ktrac,      ksfc_type,   idx_wtr,                    &! in
    pdtime, pfrc,                                        &! in
    pcfm_tile,                                           &! in
    aa,         pcptgz,                                  &! in
    pum1,       pvm1,        ptm1,                       &! in
    pmair,                                               &! in
    pqm1,       pxlm1,       pxim1,       pxtm1,         &! in
    pgeom1,      pztottevn,                              &! in
    vdiff_config,                                        &! in
    bb,                                                  &! inout
    pzthvvar,   pxvar,       pz0m_tile,                  &! in, inout, inout
    pkedisp,                                             &! out
    pute_vdf,   pvte_vdf,    pq_vdf,                     &! out
    pqte_vdf,   pxlte_vdf,   pxite_vdf,   pxtte_vdf,     &! out
    pz0m,                                                &! out
    pthvvar,                                             &! out
    ptotte                                               )! out

    INTEGER, INTENT(IN) :: jcs, kproma, kbdim, klev, klevm1, ktrac
    INTEGER, INTENT(IN) :: ksfc_type, idx_wtr
    REAL(wp),INTENT(IN) :: pdtime

    REAL(wp),INTENT(IN) ::    &
    & pfrc      (:,:)     , & !< (kbdim,ksfc_type) area fraction of each surface type
    & pcfm_tile (:,:)         !< (kbdim,ksfc_type) exchange coeff

    REAL(wp),INTENT(IN) :: aa    (:,:,:,:) !< (kbdim,klev,3,nmatrix) for all variables


    ! The input variables below are needed only by "vdiff_tendencies"

    REAL(wp),INTENT(IN) :: pcptgz  (:,:)   !< (kbdim,klev) dry static energy

    REAL(wp),INTENT(IN) :: pum1    (:,:)   !< (kbdim,klev) u-wind at step t-dt
    REAL(wp),INTENT(IN) :: pvm1    (:,:)   !< (kbdim,klev) q-wind at step t-dt
    REAL(wp),INTENT(IN) :: ptm1    (:,:)   !< (kbdim,klev) temperature at step t-dt
    REAL(wp),INTENT(IN) :: pmair   (:,:)   !< (kbdim,klev) moist air mass [kg/m2]
    REAL(wp),INTENT(IN) :: pqm1    (:,:)   !< (kbdim,klev) specific humidity at step t-dt
    REAL(wp),INTENT(IN) :: pxlm1   (:,:)   !< (kbdim,klev) cloud water concentration at step t-dt
    REAL(wp),INTENT(IN) :: pxim1   (:,:)   !< (kbdim,klev) cloud ice   concentration at step t-dt
    REAL(wp),INTENT(IN) :: pxtm1   (:,:,:) !< (kbdim,klev,ktrac) specific density of other tracers at step t-dt

    REAL(wp),INTENT(IN) :: pgeom1 (:,:)   !< (kbdim,klev) geopotential above ground
    REAL(wp),INTENT(IN) :: pztottevn(:,:) !< (kbdim,klev) intermediate value of TTE

    TYPE(t_vdiff_config), INTENT(IN) :: vdiff_config !< VDIFF configuration for current domain.

    REAL(wp),INTENT(INOUT) :: bb    (:,:,:)  !< (kbdim,klev,nvar_vdiff)

    REAL(wp),INTENT(IN)    :: pzthvvar (:,:) !< (kbdim,klev) intermediate value of thvvar
    REAL(wp),INTENT(INOUT) :: pxvar    (:,:) !< (kbdim,klev) distribution width (b-a)
                              !< in: step t-dt, out: modified
                              !< due to vertical diffusion

    ! Roughness length
    ! In: values over each surface type.
    ! Out: z0m_tile over the ocean is updated using the time average ("hat" value)
    ! of u and v and the t-dt value of some other variables.

    REAL(wp),INTENT(INOUT) :: pz0m_tile (:,:) !< (kbdim,ksfc_type)

    ! Vertically integrated dissipation of kinetic energy [W/m2]

    REAL(wp),INTENT(OUT) :: pkedisp(:) !< (kbdim)

    ! Tendencies

    REAL(wp),INTENT(OUT) :: pute_vdf (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pvte_vdf (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pq_vdf   (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pqte_vdf (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pxlte_vdf(:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pxite_vdf(:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pxtte_vdf(:,:,:) !< (kbdim,klev,ktrac)

    ! Some other diagnostics

    REAL(wp),INTENT(OUT) :: pz0m      (:)     !< (kbdim) grid-box mean roughness height
    REAL(wp),INTENT(OUT) :: pthvvar   (:,:)   !< (kbdim,klev) variance of virtual potential temperature
                              !< at the new time step t
    REAL(wp),INTENT(OUT) :: ptotte    (:,:)   !< (kbdim,klev)

    !-----------------------------------------------------------------------
    ! 6. Obtain solution of the tri-diagonal system by back-substitution.
    !    Then compute tendencies and diagnose moisture flux etc.
    !-----------------------------------------------------------------------
    CALL rhs_bksub( jcs, kproma, klev, aa, bb ) ! in,...,in, inout

    CALL vdiff_tendencies( jcs, kproma, kbdim, klev, klevm1, &! in
      & ktrac, ksfc_type, idx_wtr,                   &! in
      & pdtime,                                      &! in
      & pum1, pvm1, ptm1,                            &! in
      & pmair,                                       &! in
      & pqm1, pxlm1, pxim1, pxtm1,                   &! in
      & pgeom1, pcptgz,                              &! in
      & pztottevn, pzthvvar,                         &! in
      & pcfm_tile, pfrc, bb,                         &! in
      & vdiff_config,                                &! in
      & pkedisp,                                     &! out
      & pxvar,                                       &! inout
      & pz0m_tile,                                   &! inout (out for iwtr)
      & pute_vdf, pvte_vdf, pq_vdf,                  &! out
      & pqte_vdf, pxlte_vdf, pxite_vdf, pxtte_vdf,   &! out
      & pz0m, ptotte, pthvvar                        )! out

    ! Note: computation of additional diagnostics, e.g., surface sensible heat flux,
    !       wind stress, 10m wind, 2m temperature etc., has not been implemented yet.

  END SUBROUTINE vdiff_up
  !-------------


  !>
  !!
  !! In this prototype it is assumed that the following variables are subject to
  !! turbulent mixing:
  !!
  !!    variables                              |  # of variables
  !! -----------------------------------------------------------
  !!   u, v, T, qv                             |  4
  !!   all hydrometeors                        |  khydromet
  !!   variance of cloud droplet concentration |  1
  !!   TTE                                     |  1
  !!   variance of theta_v                     |  1
  !!   additional tracers                      |  ktrac
  !! -----------------------------------------------------------
  !!
  SUBROUTINE init_vdiff_solver( khydromet, ktrac )

    INTEGER,INTENT(IN) :: khydromet, ktrac
    INTEGER :: ist

    !------------------------------------------
    ! Set up index for prognostic variables
    !------------------------------------------

    nvar_vdiff = 7 + khydromet + ktrac

    iu   = 1;   iv   = 2
    ixl  = 3;   ixi  = 4;  ixv = 5
    itotte= 6;  ithv = 7
    ih   = 8;   iqv  = 9

    !>KF suggestion
    IF((7 + khydromet) > iqv ) &
    CALL finish( TRIM(thismodule),'matrix for vdiff is not properly defined')

    IF(ktrac > 0)  THEN
      itrc_start = 7 + khydromet +1
    ELSE
      itrc_start = 7 + khydromet
    ENDIF
    !<KF

    !-------------------------------------------------------------------
    ! # of vertical levels on which the prognostic equations are solved
    !-------------------------------------------------------------------

    ALLOCATE( ibtmoffset_var(nvar_vdiff),STAT=ist )
    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),&
      & 'Allocation of ibtmoffset_var failed')

    ! momentum, heat, water substances and tracers are solved on
    ! klev full levels

    ibtmoffset_var(:)    = 0

    ! TTE and the variance of $\theta_v$ are solved at klev-1 half levels.
    ! The upper and lower boundaries of the atmosphere are excluded.

    ibtmoffset_var(itotte) = -1
    ibtmoffset_var(ithv) = -1

    !------------------------------------------
    ! Set up matrix indices
    !------------------------------------------

    ALLOCATE( matrix_idx(nvar_vdiff),STAT=ist )
    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),&
      & 'Allocation of matrix_idx failed')

    matrix_idx(iu)   = 1  ; imuv = 1
    matrix_idx(iv)   = 1  ! u and v share the same exchange coeff.
    matrix_idx(ixl)  = 2
    matrix_idx(ixi)  = 2  ! cloud water and ice share the same exchange coeff.
    matrix_idx(ixv)  = 3
    matrix_idx(itotte) = 4
    matrix_idx(ithv) = 5
    matrix_idx(ih)   = 6 ; imh  = 6
    matrix_idx(iqv)  = 7 ; imqv = 7

    IF (ktrac>0) matrix_idx(nvar_vdiff-ktrac+1:nvar_vdiff) = 2

    nmatrix = 7    ! total number of matrices

    !---------------------------------------------------------------------------
    ! # of vertical levels on which elimination of the coefficient will be done
    !---------------------------------------------------------------------------

    ALLOCATE( ibtmoffset_mtrx(nmatrix),STAT=ist )
    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule),&
    & 'Allocation of ibtmoffset_mtrx failed')

    ibtmoffset_mtrx(:)                = 0
    ibtmoffset_mtrx(matrix_idx(itotte)) = -1
    ibtmoffset_mtrx(matrix_idx(ithv)) = -1

    !$ACC ENTER DATA COPYIN(matrix_idx, ibtmoffset_mtrx, ibtmoffset_var)

  END SUBROUTINE init_vdiff_solver
  !-------------
  !>
  !!
  SUBROUTINE cleanup_vdiff_solver

    INTEGER :: ist

    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(matrix_idx, ibtmoffset_mtrx, ibtmoffset_var)
    DEALLOCATE( matrix_idx,ibtmoffset_mtrx,ibtmoffset_var, STAT=ist)
    IF (ist/=SUCCESS) CALL finish('cleanup_vdiff_solver','Deallocation failed')

  END SUBROUTINE cleanup_vdiff_solver
  !-------------
  !>
  !!
  !! Set up coeffient matrix of the linear algebraic system and
  !! perform Gauss elimination. For moisture, the last row of the
  !! matrix (aa_btm) can not be finished yet because the evapotranspiration
  !! coefficients "cair" and "csat" are not yet available. Thus for this variable
  !! elimination is performed only till level klev-1.
  !!
  SUBROUTINE matrix_setup_elim( jcs, kproma, kbdim, klev, klevm1, &! in
                              & ksfc_type,                    &! in
                              & pcfm, pcfh, pcfh_tile, pcfv,  &! in
                              & pcftotte, pcfthv,             &! in
                              & pprfac,                       &! in
                              & prmairm, prmairh,             &! in
                              & aa, aa_btm                    )! out
    ! Arguments

    INTEGER, INTENT(IN) :: jcs, kproma, kbdim, klev, klevm1, ksfc_type

    REAL(wp),INTENT(IN) :: pcfm     (:,:)   !< (kbdim,klev) exchange coeff. for u, v
    REAL(wp),INTENT(IN) :: pcfh     (:,:)   !< (kbdim,klevm1) exchange coeff. for heat and tracers
    REAL(wp),INTENT(IN) :: pcfh_tile(:,:)   !< (kbdim,ksfc_type) exchange coeff. for heat and qv, at surface
    REAL(wp),INTENT(IN) :: pcfv     (:,:)   !< (kbdim,klev) exchange coeff. for total water variance
    REAL(wp),INTENT(IN) :: pcftotte (:,:)   !< (kbdim,klev) exchange coeff. for TTE
    REAL(wp),INTENT(IN) :: pcfthv   (:,:)   !< (kbdim,klev) exchange coeff. for variance of theta_v
    REAL(wp),INTENT(IN) :: pprfac   (:,:)   !< (kbdim,klev) prefactor for the exchange coefficients
    REAL(wp),INTENT(IN) :: prmairm  (:,:)   !< (kbdim,klev) reciprocal of layer air mass, full levels
    REAL(wp),INTENT(IN) :: prmairh  (:,:)   !< (kbdim,klevm1) reciprocal of layer air mass, half levels

    REAL(wp),INTENT(OUT) :: aa    (:,:,:,:)     !< (kbdim,klev,3,nmatrix) exchange coeff. matrices    out
    REAL(wp),INTENT(OUT) :: aa_btm(:,:,:,imh:)  !< (kbdim,3,ksfc_type,imh:imqv) out
                                                !< last (the klev-th) row of the coeff. matrices
                                                !< of dry static energy and moisture

    ! Local variables

    REAL(wp) :: zkstar (kbdim,0:klev)     !< scaled exchange coeff on half-levels
    REAL(wp) :: zkh    (kbdim,0:klevm1)   !< scaled exchange coeff on full-levels,
                                          !< for TTE and variance of theta_v
    INTEGER  :: im             !< index of coefficient matrix
    INTEGER  :: jc, jk, jsfc   !< loop indices
    INTEGER  :: jkm1, jmax

    !---- Local Variables
    !$ACC DATA CREATE(zkstar, zkh)

    !-----------------------------------------------------------------------
    ! For all prognostic variables: no turbulent flux at the upper boundary
    !-----------------------------------------------------------------------

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs,kproma
      zkstar(jc,0) = 0._wp
    ENDDO
    !$ACC END PARALLEL

    !-----------------------------------------------------------------------
    ! For momentum: surface flux is considered
    !-----------------------------------------------------------------------
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,klev
      DO jc = jcs,kproma
          zkstar(jc,jk) = pprfac(jc,jk)*pcfm(jc,jk)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    im = matrix_idx(iu)    ! also = matrix_idx(iv)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,klev
      DO jc = jcs,kproma
        aa(jc,jk,1,im) = -zkstar(jc,jk-1)*prmairm(jc,jk)  ! -K*_{k-1/2}/dm_k
        aa(jc,jk,3,im) = -zkstar(jc,jk  )*prmairm(jc,jk)  ! -K*_{k+1/2}/dm_k
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO
    !$ACC END PARALLEL


    !---------------------------------------------------------------------
    ! Dry static energy: surface fluxes on different surface types
    ! are handled separately.
    !---------------------------------------------------------------------
    im = imh
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,klevm1
      DO jc = jcs,kproma
        zkstar(jc,jk) =  pprfac(jc,jk) &
                                  &   *pcfh(jc,jk)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,klevm1
      DO jc = jcs,kproma
        aa(jc,jk,1,im) = -zkstar(jc,jk-1)*prmairm(jc,jk)  ! -K*_{k-1/2}/dm_k
        aa(jc,jk,3,im) = -zkstar(jc,jk  )*prmairm(jc,jk)  ! -K*_{k+1/2}/dm_k
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    ! Set the bottom row of the coeff matrix. The same formula applies
    ! for all surface types (land, water, ice).

    jk = klev
    IF ( isrfc_type == 1 ) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jsfc = 1,ksfc_type
        DO jc = jcs,kproma
          aa_btm(jc,1,jsfc,im) = -zkstar(jc,jk-1)*prmairm(jc,jk)    ! -K*_{k-1/2}/dm_k
          aa_btm(jc,3,jsfc,im) = -pcfh_tile(jc,jsfc)*pprfac(jc,jk)*prmairm(jc,jk)
          aa_btm(jc,2,jsfc,im) = 1._wp - aa_btm(jc,1,jsfc,im)
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ELSE
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jsfc = 1,ksfc_type
        DO jc = jcs,kproma
          aa_btm(jc,1,jsfc,im) = -zkstar(jc,jk-1)*prmairm(jc,jk)    ! -K*_{k-1/2}/dm_k
          aa_btm(jc,3,jsfc,im) = -pcfh_tile(jc,jsfc)*pprfac(jc,jk)*prmairm(jc,jk)
          aa_btm(jc,2,jsfc,im) = 1._wp - aa_btm(jc,1,jsfc,im) - aa_btm(jc,3,jsfc,im)
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    END IF

    !---------------------------------------------------------------------
    ! Moisture: different surface types are handled separately.
    !---------------------------------------------------------------------
    im = imqv
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,klevm1
      DO jc = jcs,kproma
      zkstar(jc,jk) =  pprfac(jc,jk) * pcfh(jc,jk)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,klevm1
      DO jc = jcs,kproma
        aa(jc,jk,1,im) = -zkstar(jc,jk-1)*prmairm(jc,jk)  ! -K*_{k-1/2}/dm_k
        aa(jc,jk,3,im) = -zkstar(jc,jk  )*prmairm(jc,jk)  ! -K*_{k+1/2}/dm_k
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    ! Bottom row of the matrix: finish the setup over water and ice;
    ! do part of the computation for land. Later in subroutine
    ! matrix_to_richtmyer_coeff, aa_btm(:,3,idx_land,imqv) will be
    ! modified, and aa_btm(:,2,idx_land,imqv) re-computed.

    jk = klev
    IF ( isrfc_type == 1 ) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jsfc = 1,ksfc_type
        DO jc = jcs,kproma
          aa_btm(jc,1,jsfc,im) = -zkstar(jc,jk-1)*prmairm(jc,jk)    ! -K*_{k-1/2}/dm_k
          aa_btm(jc,3,jsfc,im) = -pcfh_tile(jc,jsfc)*pprfac(jc,jk)*prmairm(jc,jk)
          aa_btm(jc,2,jsfc,im) = 1._wp - aa_btm(jc,1,jsfc,im)
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ELSE
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jsfc = 1,ksfc_type
        DO jc = jcs,kproma
          aa_btm(jc,1,jsfc,im) = -zkstar(jc,jk-1)*prmairm(jc,jk)    ! -K*_{k-1/2}/dm_k
          aa_btm(jc,3,jsfc,im) = -pcfh_tile(jc,jsfc)*pprfac(jc,jk)*prmairm(jc,jk)
          aa_btm(jc,2,jsfc,im) = 1._wp - aa_btm(jc,1,jsfc,im) - aa_btm(jc,3,jsfc,im)
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    END IF

    !----------------------------------------------------------------------
    ! For all advected tracers except water vapour: no turbulent flux at
    ! the surface.
    !----------------------------------------------------------------------
    !im = matrix_idx(ixl)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs,kproma
      zkstar(jc,klev) = 0._wp  ! lower boundary, no turbulent flux
    END DO
    !$ACC END PARALLEL

    im = matrix_idx(ixl)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,klev
      DO jc = jcs,kproma
        aa(jc,jk,1,im) = -zkstar(jc,jk-1)*prmairm(jc,jk)  ! -K*_{k-1/2}/dm_k
        aa(jc,jk,3,im) = -zkstar(jc,jk  )*prmairm(jc,jk)  ! -K*_{k+1/2}/dm_k
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !----------------------------------------------------------------------
    ! For total water variance: no surface flux. The exchange coefficient
    ! pcfv has been set to to zero in subroutine sfc_exchange_coeff, which
    ! automatically leads to zkstar(:,klev) = 0._wp, thus no additional
    ! attention is needed here.
    !----------------------------------------------------------------------

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,klev
      DO jc = jcs,kproma
      zkstar(jc,jk) =  pprfac(jc,jk) &
                                  &   *pcfv(jc,jk)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !zkstar(1:kproma,1:klev) = pprfac(1:kproma,1:klev) &
    !                           &  *pcfv(1:kproma,1:klev)

    im = matrix_idx(ixv)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,klev
      DO jc = jcs,kproma
        aa(jc,jk,1,im) = -zkstar(jc,jk-1)*prmairm(jc,jk)
        aa(jc,jk,3,im) = -zkstar(jc,jk  )*prmairm(jc,jk)
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !------------------------------------------------------------------------
    ! For TTE: Note that
    ! - Vertical averaging is needed to convert exchange coefficient from
    !   half to full levels, because TTE equation is solved on half levels.
    ! - TTE equation is solved only till array subscript klevm1, which
    !   corresponds to half level (klev - 1/2), i.e., the lowest
    !   interface above surface. Surface value of TTE is (already)
    !   computed in subroutine "sfc_exchange_coeff".
    !------------------------------------------------------------------------

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,klev
      DO jc = jcs,kproma
      zkstar(jc,jk) =  pprfac(jc,jk) &
                                  &   *pcftotte(jc,jk)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,klevm1
      DO jc = jcs,kproma
        zkh(jc,jk) = 0.5_wp*(zkstar(jc,jk)+zkstar(jc,jk+1))
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs,kproma
      zkh(jc,0) = 0._wp  ! upper boundary, no flux
    ENDDO
    !$ACC END PARALLEL

    im = matrix_idx(itotte)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,klevm1
      DO jc = jcs,kproma
        aa(jc,jk,1,im) = -zkh(jc,jk-1)*prmairh(jc,jk)
        aa(jc,jk,3,im) = -zkh(jc,jk  )*prmairh(jc,jk)
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !------------------------------------------------
    ! For the variance of theta_v (similar to TTE)
    !------------------------------------------------
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,klev
      DO jc = jcs,kproma
      zkstar(jc,jk) =  pprfac(jc,jk) &
                                  &   *pcfthv(jc,jk)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,klevm1
      DO jc = jcs,kproma
        zkh(jc,jk) = 0.5_wp*(zkstar(jc,jk)+zkstar(jc,jk+1))
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    im = matrix_idx(ithv)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,klevm1
      DO jc = jcs,kproma
        aa(jc,jk,1,im) = -zkh(jc,jk-1)*prmairh(jc,jk)
        aa(jc,jk,3,im) = -zkh(jc,jk  )*prmairh(jc,jk)
        aa(jc,jk,2,im) = 1._wp - aa(jc,jk,1,im) - aa(jc,jk,3,im)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !-----------------------------------------------------------------------------
    ! Gauss elimination for the coefficient matrices at
    ! - vertical levels [1,klev-2], for TTE and variance of theta_v;
    ! - vertical levels [1,klev-1], for all the other variables.
    !-----------------------------------------------------------------------------


    DO im = 1,nmatrix
      DO jc = jcs,kproma
        aa(jc,1,3,im) = aa(jc,1,3,im)/aa(jc,1,2,im)
      ENDDO

      jmax = klev + ibtmoffset_mtrx(im) - 1
      DO jk = 2,jmax
        jkm1 = jk - 1
        DO jc = jcs,kproma
          aa(jc,jk,2,im) =  aa(jc,jk,2,im)                       &
                            & -aa(jc,jk,1,im)*aa(jc,jkm1,3,im)
          aa(jc,jk,3,im) =  aa(jc,jk,3,im)/aa(jc,jk,2,im)
        ENDDO
      ENDDO
    END DO


    ! Translation for developers who prefer to think in terms of
    ! the Richtmyer-Morthon formula and are familiar with the paper by
    ! Polcher et al (1998): after this elimination,
    !  aa(:,1:klev+ibtmoffset_mtrx(im)-1,2,:) becomes C  (Eqn. 17),
    !  aa(:,1:klev+ibtmoffset_mtrx(im)-1,3,:) becomes -A (Eqn. 19).
    ! See subroutine matrix_to_richtmyer_coeff.

  !$ACC WAIT
  !$ACC END DATA

  END SUBROUTINE matrix_setup_elim

  !--------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------
  !>
  SUBROUTINE rhs_setup( jcs, kproma, kbdim, klev, klevm1,    &! in
                      & ksfc_type, ktrac, pdtime,            &! in
                      & pum1, pvm1, pcptgz, pqm1,            &! in
                      & pxlm1, pxim1, pxvar, pxtm1, pxt_emis,&! in
                      & prmairm, ptottevn, pzthvvar, aa,     &! in
                      & bb, bb_btm                           )! out

    ! Arguments

    INTEGER, INTENT(IN) :: jcs, kproma, kbdim, klev, klevm1
    INTEGER, INTENT(IN) :: ksfc_type, ktrac
    REAL(wp),INTENT(IN) :: pdtime

    REAL(wp),INTENT(IN) :: pum1     (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN) :: pvm1     (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN) :: pcptgz   (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN) :: pqm1     (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN) :: pxlm1    (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN) :: pxim1    (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN) :: pxvar    (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN) :: pxtm1    (:,:,:) !< (kbdim,klev,ktrac)
    REAL(wp),INTENT(IN) :: pxt_emis (:,:)   !< (kbdim,ktrac)
    !REAL(wp),INTENT(IN) :: pxt_emis (:,:,:) ! (kbdim,klev,ktrac) backup for later use
    REAL(wp),INTENT(IN) :: ptottevn (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN) :: pzthvvar (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN) :: prmairm  (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN) :: aa       (:,:,:,:) !< (kbdim,klev,3,nmatrix)

    REAL(wp),INTENT(OUT) :: bb    (:,:,:)   !< (kbdim,klev,nvar_vdiff) OUT
    REAL(wp),INTENT(OUT) :: bb_btm(:,:,ih:) !< (kbdim,ksfc_type,ih:iqv) OUT

    ! Local variables

    REAL(wp) :: ztmp(kbdim,klev)
    INTEGER  :: jsfc, jt, irhs, im, jk, jc

    !---- Local Variables
    !$ACC DATA &
    !$ACC   CREATE(ztmp)

    !-------------------------------------------------------------------
    ! First handle variables that are defined on full levels
    !-------------------------------------------------------------------
    ! u and v
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,klev
      DO jc = jcs,kproma
        bb(jc,jk,iu) = pum1(jc,jk)
        bb(jc,jk,iv) = pvm1(jc,jk)

    ! Hydrometeors and the variance of cloud droplets

        bb(jc,jk,ixl) = pxlm1(jc,jk)
        bb(jc,jk,ixi) = pxim1(jc,jk)
        bb(jc,jk,ixv) = pxvar(jc,jk)
      END DO
    END DO
    !$ACC END PARALLEL

    ! Other tracers

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(ktrac > 0)
    !$ACC LOOP SEQ
    DO jt = 1,ktrac
      irhs = jt - 1 + itrc_start
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1,klev
          DO jc = jcs,kproma
            bb(jc,jk,irhs) =  pxtm1(jc,jk,jt)
          END DO
        END DO
        !bb(1:kproma,1:klev,irhs) =  pxtm1(1:kproma,1:klev,jt)
    ENDDO
    !$ACC END PARALLEL

    ! Heat and moisture

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,klevm1
      DO jc = jcs,kproma
        bb(jc,jk,ih ) = pcptgz(jc,jk)
        bb(jc,jk,iqv) = pqm1  (jc,jk)
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jsfc = 1,ksfc_type
      DO jc = jcs,kproma
        bb_btm(jc,jsfc,ih)  = pcptgz(jc,klev)
        bb_btm(jc,jsfc,iqv) =   pqm1(jc,klev)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !-------------------------------------------------------------------
    ! TTE and the variance of theta_v:
    ! These variables are defined at half levels. Array index jk
    ! correspond to half level k+1/2. Thus klev correspond to the
    ! lower boundary. The linear solver only solves till index klevm1.
    !-------------------------------------------------------------------

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,klevm1
      DO jc = jcs,kproma
        bb(jc,jk,itotte) =  ptottevn(jc,jk)
      END DO
    END DO
    !$ACC END PARALLEL

    im = matrix_idx(itotte)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs,kproma
      bb(jc,     klevm1,itotte) =  bb(jc,klevm1,itotte)   &
                              & -aa(jc,klevm1,3,im)   &
                              & *ptottevn(jc,klev)
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,klevm1
      DO jc = jcs,kproma
        bb(jc,jk,ithv) =  pzthvvar(jc,jk)
      END DO
    END DO
    !$ACC END PARALLEL

    im = matrix_idx(ithv)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs,kproma
      bb(jc,     klevm1,ithv) =  bb(jc,klevm1,ithv)   &
                              & -aa(jc,klevm1,3,im)   &
                              & *pzthvvar(jc,klev)
    ENDDO
    !$ACC END PARALLEL

    !--------------------------------------------------------------------
    ! Apply the implicitness factor
    !--------------------------------------------------------------------
    !bb     = tpfac2*bb
    !bb_btm = tpfac2*bb_btm

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(3)
    DO jt = 1, itotte-1
      DO jk = 1,klev
        DO jc = jcs,kproma
          bb(jc,jk,jt)  = tpfac2*bb(jc,jk,jt)
        ENDDO
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(3)
    DO jt = itotte, iqv
      DO jk = 1,klevm1
        DO jc = jcs,kproma
          bb(jc,jk,jt)  = tpfac2*bb(jc,jk,jt)
        ENDDO
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    IF (ktrac>0) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(3)
      DO jt = itrc_start, nvar_vdiff
        DO jk = 1,klev
          DO jc = jcs,kproma
            bb(jc,jk,jt)  = tpfac2*bb(jc,jk,jt)
          ENDDO
        ENDDO
      ENDDO
    !$ACC END PARALLEL

    ENDIF

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(3)
    DO jt = ih,iqv
      DO jk = 1,ksfc_type
        DO jc = jcs,kproma
          bb_btm(jc,jk,jt)  = tpfac2*bb_btm(jc,jk,jt)
        ENDDO
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !--------------------------------------------------------------------
    ! Add tracer emissions
    !--------------------------------------------------------------------
    ! Currently we follow ECHAM in which only the surface emission
    ! is treated in "vdiff".
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs,kproma
      ztmp(jc,klev) = prmairm(jc,klev)*pdtime
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(ktrac > 0)
    !$ACC LOOP SEQ
    DO jt = 1,ktrac
      irhs = jt - 1 + itrc_start
        !$ACC LOOP GANG VECTOR
        DO jc = jcs,kproma
          bb(jc,klev,irhs) =         bb(jc,klev,irhs) &
                              & + pxt_emis(jc,jt)        &
                              &      *ztmp(jc,klev)
        END DO
    ENDDO
    !$ACC END PARALLEL


    !DO jt = 1,ktrac
    !   irhs = jt - 1 + itrc_start
    !   bb(1:kproma,klev,irhs) =         bb(1:kproma,klev,irhs) &
    !                          & + pxt_emis(1:kproma,jt)        &
    !                          &      *ztmp(1:kproma,klev)
    !ENDDO

    ! Later we may consider treating emission on all vertical levels
    ! in the same way.
    !
    !ztmp(jcs:kproma,1:klev) = prmairm(jcs:kproma,1:klev)*pdtime
    !
    !DO jt = 1,ktrac
    !   irhs = jt - 1 + itrc_start
    !   bb(jcs:kproma,1:klev,irhs) =         bb(jcs:kproma,1:klev,irhs) &
    !                               & + pxt_emis(jcs:kproma,1:klev,jt)   &
    !                               &      *ztmp(jcs:kproma,1:klev)
    !ENDDO
  !$ACC WAIT
  !$ACC END DATA


  END SUBROUTINE rhs_setup

  !--------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------
  !>
  !!
  !! Gauss elimination of the right-hand-side vector
  !! using coefficients obtained in subroutine "matrix_setup_elim".
  !!
  !! Translation for developers who prefer to think in terms of
  !! the Richtmyer-Morthon formula and are familiar with the paper by
  !! Polcher et al (1998): after the elimination at the end of
  !! subroutine matrix_setup_elim, aa(:,1:klev+ibtmoffset_mtrx(im)-1,2,:)
  !! became the coeff C defined by Eqn. 17 of Polcher et al (1998).
  !! It is used in this subroutine to convert the variable bb
  !! into the Richtmyer coeff B (cf Eqn. 19 of Polcher et al 1998).
  !!
  SUBROUTINE rhs_elim( jcs, kproma, klev, aa,   &! in
                      & bb                      )! in, inout

    INTEGER, INTENT(IN)    :: jcs, kproma, klev
    REAL(wp),INTENT(IN)    :: aa(:,:,:,:) !< (jcs:kproma,klev,3,nmatrix)
    REAL(wp),INTENT(INOUT) :: bb(:,:,:)   !< (jcs:kproma,klev,nvar_vdiff)

    REAL(wp) :: znum, zden
    INTEGER  :: jvar, im, jk, jkm1, jmax, jc

    ! 1. Vertical levels [2,klev-2] for TTE and variance of theta_v;
    !    [2,klev-1] for all the other variables.

    DO jvar = 1,nvar_vdiff
      im = matrix_idx(jvar)  ! Index of coefficient matrix
      DO jc = jcs, kproma
        bb(jc,1,jvar) =  bb(jc,1,jvar)/aa(jc,1,2,im)
      ENDDO

      jmax = klev + ibtmoffset_var(jvar) - 1

      DO jk = 2,jmax
        jkm1 = jk - 1
        DO jc = jcs,kproma
          znum =  bb(jc,jk  ,jvar)                     &
                    & -bb(jc,jkm1,jvar)*aa(jc,jk,1,im)
          bb(jc,jk,jvar) = znum/aa(jc,jk,2,im)
        ENDDO
      ENDDO
    ENDDO !jvar: variable loop

    ! 2. Bottom level for all variables except u, v, dry static energy
    !    and moisture. After this step the array bb contains the
    !    solution of the linear system.

    DO jvar = 1,nvar_vdiff

      IF (jvar==iu.OR.jvar==iv.OR.jvar==ih.OR.jvar==iqv ) THEN
          CYCLE
      ELSE

      im   = matrix_idx(jvar)  ! Index of coefficient matrix
      jk   = klev + ibtmoffset_var(jvar)  ! Bottom level index
      jkm1 = jk - 1

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jc = jcs,kproma
        zden =  aa(jc,jk,2,im)                      &
                      & -aa(jc,jk,1,im)*aa(jc,jkm1,3,im)
        znum =  bb(jc,jk,jvar)                      &
                      & -aa(jc,jk,1,im)*bb(jc,jkm1,jvar)
        bb(jc,jk,jvar) = znum/zden
      ENDDO
      !$ACC END PARALLEL

      END IF
    ENDDO !jvar: variable loop

    ! Note that for TTE and the variance of theta_v, klev-1 is the lowest
    ! level above surface. Now set boundary condition for the variance
    ! of theta_v.

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = jcs,kproma
      bb(jc,klev,ithv) = bb(jc,klev-1,ithv)
    ENDDO
    !$ACC END PARALLEL

    !$ACC WAIT

  END SUBROUTINE rhs_elim

  !--------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------
  !>
  !!
  !! Prepare the Richtmyer-Morton coeffcients for dry static energy and
  !! moisture, to be used by the surface models (ocean, sea-ice, land).
  !!
  SUBROUTINE matrix_to_richtmyer_coeff( jcs, kproma, klev, ksfc_type, idx_lnd, &! in
                                      & aa, bb,                                      &! in
                                      & pdtime, delz,                                &! in
                                      & aa_btm, bb_btm,                              &! inout
                                      & pen_h, pfn_h, pen_qv, pfn_qv,                &! out
                                      & pcair,                                       &! in
                                      & pcsat)                                        ! in

    INTEGER,INTENT(IN)     :: jcs, kproma, klev, ksfc_type, idx_lnd
    REAL(wp),INTENT(IN)    :: aa    (:,:,:,imh:) !< (jcs:kproma,klev,3,imh:imqv)
    REAL(wp),INTENT(IN)    :: bb    (:,:,ih:)    !< (jcs:kproma,klev,ih:iqv)
    REAL(wp),INTENT(IN)    :: pdtime
    REAL(wp),INTENT(IN)    :: delz(:)            !< (jcs:kproma)
    REAL(wp),INTENT(INOUT) :: aa_btm(:,:,:,imh:) !< (jcs:kproma,3,ksfc_type,imh:imqv)
    REAL(wp),INTENT(INOUT) :: bb_btm(:,:,ih:)    !< (jcs:kproma,ksfc_type,ih:iqv)

    REAL(wp),INTENT(OUT) :: pen_h (:,:)  !< (jcs:kproma,ksfc_type) OUT
    REAL(wp),INTENT(OUT) :: pfn_h (:,:)  !< (jcs:kproma,ksfc_type) OUT
    REAL(wp),INTENT(OUT) :: pen_qv(:,:)  !< (jcs:kproma,ksfc_type) OUT
    REAL(wp),INTENT(OUT) :: pfn_qv(:,:)  !< (jcs:kproma,ksfc_type) OUT

    REAL(wp),OPTIONAL,INTENT(IN)    :: pcair(:) !< (jcs:kproma)
    REAL(wp),OPTIONAL,INTENT(IN)    :: pcsat(:) !< (jcs:kproma)

    INTEGER  :: jk, jsfc

    !---------------------------------------------------------
    ! Matrix setup and bottom level elimination for moisture
    !---------------------------------------------------------
    ! Evapotranspiration has to be considered over land

    IF (ljsb .AND. idx_lnd<=ksfc_type) THEN

      jsfc = idx_lnd

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jk = jcs, kproma
        aa_btm(jk,2,jsfc,imqv) =           1._wp - aa_btm(jk,1,jsfc,imqv) &
                                    & - pcair(jk)*aa_btm(jk,3,jsfc,imqv)
        aa_btm(jk,3,jsfc,imqv) =   pcsat(jk)*aa_btm(jk,3,jsfc,imqv)
      END DO
      !$ACC END PARALLEL

    END IF ! ljsbach

    ! Bottom level elimination for all surface types

    IF ( isrfc_type == 1) THEN ! 1 = fixed surface heat fluxes
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jsfc = 1,ksfc_type
        DO jk = jcs, kproma
          aa_btm(jk,2,jsfc,imqv) =  aa_btm(jk,2,jsfc,imqv)  &
                                  & -aa_btm(jk,1,jsfc,imqv)  &
                                  & *aa    (jk,klev-1,3,imqv)

          aa_btm(jk,3,jsfc,imqv) =  -lhflx*pdtime/delz(jk) &
                                  & /aa_btm(jk,2,jsfc,imqv)

          bb_btm(jk,jsfc,iqv)    = (bb_btm(jk,jsfc,iqv)    &
                                  & -aa_btm(jk,1,jsfc,imqv) &
                                  & *bb    (jk,klev-1,iqv) )&
                                  & /aa_btm(jk,2,jsfc,imqv)

        END DO
      END DO
      !$ACC END PARALLEL
    ELSE
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jsfc = 1,ksfc_type
        DO jk = jcs, kproma
          aa_btm(jk,2,jsfc,imqv) =  aa_btm(jk,2,jsfc,imqv)  &
                                  & -aa_btm(jk,1,jsfc,imqv)  &
                                  & *aa    (jk,klev-1,3,imqv)

          aa_btm(jk,3,jsfc,imqv) =  aa_btm(jk,3,jsfc,imqv)  &
                                  & /aa_btm(jk,2,jsfc,imqv)

          bb_btm(jk,jsfc,iqv)    = (bb_btm(jk,jsfc,iqv)    &
                                  & -aa_btm(jk,1,jsfc,imqv) &
                                  & *bb    (jk,klev-1,iqv) )&
                                  & /aa_btm(jk,2,jsfc,imqv)
        END DO
      END DO
      !$ACC END PARALLEL
    END IF

    !---------------------------------------------------------
    ! Bottom level elimination for dry static energy
    !---------------------------------------------------------
    IF ( isrfc_type == 1 ) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jsfc = 1,ksfc_type
        DO jk = jcs, kproma
          aa_btm(jk,2,jsfc,imh) =  aa_btm(jk,2,jsfc,imh) &
                                      & -aa_btm(jk,1,jsfc,imh) &
                                      & *aa    (jk,klev-1,3,imh)

          aa_btm(jk,3,jsfc,imh) =  -shflx*cpd*pdtime/delz(jk) &
                                      & /aa_btm(jk,2,jsfc,imh)

          bb_btm(jk,jsfc,ih)    = (bb_btm(jk,jsfc,ih)    &
                                      & -aa_btm(jk,1,jsfc,imh) &
                                      & *bb    (jk,klev-1,ih) )&
                                      & /aa_btm(jk,2,jsfc,imh)
        END DO
      END DO
      !$ACC END PARALLEL
    ELSE
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jsfc = 1,ksfc_type
        DO jk = jcs, kproma

          aa_btm(jk,2,jsfc,imh) =  aa_btm(jk,2,jsfc,imh) &
                                      & -aa_btm(jk,1,jsfc,imh) &
                                      & *aa    (jk,klev-1,3,imh)

          aa_btm(jk,3,jsfc,imh) =  aa_btm(jk,3,jsfc,imh) &
                                      & /aa_btm(jk,2,jsfc,imh)

          bb_btm(jk,jsfc,ih)    = (bb_btm(jk,jsfc,ih)    &
                                      & -aa_btm(jk,1,jsfc,imh) &
                                      & *bb    (jk,klev-1,ih) )&
                                      & /aa_btm(jk,2,jsfc,imh)
        END DO
      END DO
      !$ACC END PARALLEL
    END IF

    !---------------------------------------------------------
    ! Convert matrix entries to Richtmyer-Morton coefficients
    !---------------------------------------------------------

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jsfc = 1,ksfc_type
      DO jk = jcs, kproma
        pen_h (jk,jsfc) = -aa_btm(jk,3,jsfc,imh)
        pen_qv(jk,jsfc) = -aa_btm(jk,3,jsfc,imqv)

        pfn_h (jk,jsfc) =  bb_btm(jk,jsfc,ih )*tpfac1
        pfn_qv(jk,jsfc) =  bb_btm(jk,jsfc,iqv)*tpfac1
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC WAIT

  END SUBROUTINE matrix_to_richtmyer_coeff
  !--------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------
  !>
  !!
  !! Do Back-substitution to get the solution of the linear system.
  !!
  !! Translation for developers who prefer to think in terms of
  !! the Richtmyer-Morthon formula and are familiar with the paper by
  !! Polcher et al (1998): on entry bb contains the solution
  !! at the bottom level and the coeff B on upper levels;
  !! On exit it becomes the solution of the linear system.
  !! aa(:,:,3,:) used here corresponds to -A in the Appendix of
  !! Polcher et al (1998).
  !! Note that VDIFF uses the implicit time stepping as in IFS
  !! in contrast to Polcher et al (1998). Thus the solution is
  !! not yet the new value at time step t+dt.
  !!
  SUBROUTINE rhs_bksub( jcs, kproma, klev, aa, bb )

    INTEGER, INTENT(IN)   :: jcs, kproma, klev
    REAL(wp),INTENT(IN)   :: aa(:,:,:,:) !< (jcs:kproma,klev,3,nmatrix)
    REAL(wp),INTENT(INOUT):: bb(:,:,:)   !< (jcs:kproma,klev,nvar_vdiff)

    INTEGER  :: jvar, im, jk, jl, jmax

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jvar = 1,nvar_vdiff
      DO jl = jcs,kproma

        im   = matrix_idx(jvar)
        jmax = klev + ibtmoffset_var(jvar) - 1

        !$ACC LOOP SEQ
        DO jk = jmax,1,-1
          bb(jl,jk,jvar) =  bb(jl,jk ,jvar) &
                                & -bb(jl,jk+1,jvar) &
                                & *aa(jl,jk  ,3,im)
        ENDDO
      ENDDO
    END DO
    !$ACC END PARALLEL

    !$ACC WAIT

  END SUBROUTINE rhs_bksub
  !-------------
  !>
  !!
  SUBROUTINE vdiff_tendencies( jcs, kproma, kbdim, klev, klevm1,            &! in
                              & ktrac, ksfc_type, idx_wtr,                  &! in
                              & pdtime,                                     &! in
                              & pum1, pvm1, ptm1,                           &! in
                              & pmair,                                      &! in
                              & pqm1, pxlm1, pxim1, pxtm1,                  &! in
                              & pgeom1, pcptgz,                             &! in
                              & pztottevn, pzthvvar,                        &! in
                              & pcfm_tile, pfrc, bb,                        &! in
                              & vdiff_config,                               &! in
                              & pkedisp,                                    &! out
                              & pxvar, pz0m_tile,                           &! inout
                              & pute_vdf, pvte_vdf, pq_vdf,                 &! out
                              & pqte_vdf, pxlte_vdf, pxite_vdf, pxtte_vdf,  &! out
                              & pz0m, ptotte, pthvvar                       )! out

    INTEGER, INTENT(IN) :: jcs, kproma, kbdim, klev, klevm1, ktrac !!$, klevp1
    INTEGER, INTENT(IN) :: ksfc_type, idx_wtr
    REAL(wp),INTENT(IN) :: pdtime

    REAL(wp),INTENT(IN)  :: pum1   (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pvm1   (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN)  :: ptm1   (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pmair  (:,:)   !< (kbdim,klev) moist air mass [kg/m2]
    REAL(wp),INTENT(IN)  :: pqm1   (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pxlm1  (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pxim1  (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pxtm1  (:,:,:) !< (kbdim,klev,ktrac)
    REAL(wp),INTENT(IN)  :: pgeom1 (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pcptgz (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pztottevn(:,:) !< (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pzthvvar(:,:) !< (kbdim,klev)
    REAL(wp),INTENT(IN)  :: pcfm_tile     (:,:) !< (kbdim,ksfc_type)
    REAL(wp),INTENT(IN)  :: pfrc          (:,:) !< (kbdim,ksfc_type)
    REAL(wp),INTENT(IN)  :: bb            (:,:,:) !<(kbdim,klev,nvar_vdiff)

    TYPE(t_vdiff_config), INTENT(IN) :: vdiff_config

    REAL(wp),INTENT(OUT) :: pkedisp(:) !< (kbdim) vertically integrated dissipation
                                        !  of kinetic energy [W/m2]

    REAL(wp),INTENT(INOUT) :: pxvar    (:,:) !< (kbdim,klev)
    REAL(wp),INTENT(INOUT) :: pz0m_tile(:,:) !< (kbdim,ksfc_type)

    REAL(wp),INTENT(OUT) :: pute_vdf (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pvte_vdf (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pq_vdf   (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pqte_vdf (:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pxlte_vdf(:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pxite_vdf(:,:)   !< (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pxtte_vdf(:,:,:) !< (kbdim,klev,ktrac)

    REAL(wp),INTENT(OUT) :: pz0m     (:)   !< (kbdim)
    REAL(wp),INTENT(OUT) :: ptotte   (:,:) !< (kbdim,klev)
    REAL(wp),INTENT(OUT) :: pthvvar  (:,:) !< (kbdim,klev)
!!$    REAL(wp),INTENT(OUT) :: psh_vdiff(:)   !< (kbdim)
!!$    REAL(wp),INTENT(OUT) :: pqv_vdiff(:)   !< (kbdim)

    REAL(wp) :: ztest, zrdt
    REAL(wp) :: zunew, zvnew, zqnew, zsnew, zhnew
    REAL(wp) :: zcp
    REAL(wp) :: zdis  (kbdim,klev)
    REAL(wp) :: z0m_min

    INTEGER  :: jk, jl, jt, irhs, jsfc

    !-------------------------------------------------------------------
    ! Start GPU data region
    !-------------------------------------------------------------------
    !$ACC DATA CREATE(zdis)

    zrdt   = 1._wp/pdtime

    z0m_min = vdiff_config%z0m_min

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1, klev
      DO jl = 1, kbdim
        pute_vdf (jl,jk)   = 0._wp
        pvte_vdf (jl,jk)   = 0._wp
        pq_vdf   (jl,jk)   = 0._wp
        pqte_vdf (jl,jk)   = 0._wp
        pxlte_vdf(jl,jk)   = 0._wp
        pxite_vdf(jl,jk)   = 0._wp
        ptotte     (jl,jk)   = 0._wp
        pthvvar  (jl,jk)   = 0._wp
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(ktrac > 0)
    !$ACC LOOP GANG VECTOR COLLAPSE(3)
    DO jt = 1, ktrac
      DO jk = 1, klev
        DO jl = 1, kbdim
          pxtte_vdf(jl,jk,jt) = 0._wp
        END DO
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jl = 1, kbdim
      pz0m     (jl)     = 0._wp
    END DO
    !$ACC END PARALLEL
    !-------------------------------------------------------------------
    ! Compute TTE at the new time step.
    !-------------------------------------------------------------------
    ztest = 0._wp
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) REDUCTION(+: ztest) ASYNC(1)
    DO jk = 1,klevm1
      DO jl = jcs,kproma
        ptotte(jl,jk) = bb(jl,jk,itotte) + tpfac3*pztottevn(jl,jk)
        ztest = ztest+MERGE(1._wp,0._wp,ptotte(jl,jk)<0._wp)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    IF( vdiff_config%turb == VDIFF_TURB_3DSMAGORINSKY ) THEN
      ztest = 1._wp
    ELSE
      IF(ztest.NE.0._wp) THEN
        CALL finish('vdiff_tendencies','TTE IS NEGATIVE')
      ENDIF
    ENDIF

    !ptotte(jcs:kproma,klev) = pztottevn(jcs:kproma,klev)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jl = jcs,kproma
      ptotte(jl,klev) = pztottevn(jl,klev)
    END DO
    !$ACC END PARALLEL

    !-------------------------------------------------------------
    ! Variance of virtual potential temperature
    !-------------------------------------------------------------
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 1,klev
      DO jl = jcs,kproma
        pthvvar(jl,jk) = bb(jl,jk,ithv) + tpfac3*pzthvvar(jl,jk)
        pthvvar(jl,jk) = MAX(totte_min,pthvvar(jl,jk))
      END DO
    END DO
    !$ACC END PARALLEL
    !-------------------------------------------------------------
    ! Tendency of velocity; kinetic energy dissipation
    !-------------------------------------------------------------

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jk = 1,kbdim
      pkedisp(jk) = 0._wp   ! initilize the vertical integral
    END DO
    !$ACC END PARALLEL


    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO jk = 1,klev
      !$ACC LOOP GANG VECTOR PRIVATE(zunew, zvnew)
      DO jl = jcs,kproma
        pute_vdf(jl,jk) = (bb(jl,jk,iu)-tpfac2*pum1(jl,jk))*zrdt
        pvte_vdf(jl,jk) = (bb(jl,jk,iv)-tpfac2*pvm1(jl,jk))*zrdt

        zunew = bb(jl,jk,iu) + tpfac3*pum1(jl,jk)
        zvnew = bb(jl,jk,iv) + tpfac3*pvm1(jl,jk)

        zdis(jl,jk) = 0.5_wp*( pum1(jl,jk)**2 - zunew**2 &
                    &         +pvm1(jl,jk)**2 - zvnew**2 )
        pkedisp(jl)  = pkedisp(jl) + zdis(jl,jk)*pmair(jl,jk)*zrdt
      END DO
    END DO
    !$ACC END PARALLEL
    !-------------------------------------------------------------
    ! Tendency of T and qv, ql, qi; xvar at the new time step
    !-------------------------------------------------------------

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(zqnew, zsnew, zhnew, zcp)
    DO jk=1,klev
      DO jl=jcs,kproma

        zqnew = bb(jl,jk,iqv) + tpfac3*pqm1(jl,jk)
        pqte_vdf(jl,jk) = (zqnew-pqm1(jl,jk))*zrdt


        ! The computation of the new temperature must be consistent with the computation of
        ! the static energy pcptgz in the subroutine mo_turbulence_diag:atm_exchange_coeff.
        ! The same specific heat must be used.
        !
        zsnew = bb(jl,jk,ih) + tpfac3*pcptgz(jl,jk)
        zhnew = (zsnew + zdis(jl,jk) - pgeom1(jl,jk))
        zcp   = cpd!+(cpv-cpd)*pqm1(jl,jk) ! cp of moist air
        !
        ! Now derive the heating for constant pressure conditions
        ! as needed for provisional updating in the physics.
        !
        pq_vdf(jl,jk)   = (zhnew - ptm1(jl,jk)*zcp)*zrdt*pmair(jl,jk)

        pxlte_vdf(jl,jk) = (bb(jl,jk,ixl) - tpfac2*pxlm1(jl,jk))*zrdt
        pxite_vdf(jl,jk) = (bb(jl,jk,ixi) - tpfac2*pxim1(jl,jk))*zrdt

        pxvar(jl,jk) = bb(jl,jk,ixv) + tpfac3*pxvar(jl,jk)
      END DO
    END DO
    !$ACC END PARALLEL

!!$    IF ( get_lebudget() ) THEN
!!$      psh_vdiff(:) = 0._wp
!!$      pqv_vdiff(:) = 0._wp
!!$      DO jk=1,klev
!!$        ! compute heat budget diagnostic
!!$        psh_vdiff(jcs:kproma) = psh_vdiff(jcs:kproma) + pmair(jcs:kproma,jk) * &
!!$        & (bb(jcs:kproma,jk,ih)  + (tpfac3 - 1._wp)*pcptgz(jcs:kproma,jk)) * zrdt
!!$        ! compute moisture budget diagnostic
!!$        ! ? zdis appears to be dissipation, probably we don't need this for qv??
!!$        pqv_vdiff(jcs:kproma) = pqv_vdiff(jcs:kproma) + pmair(jcs:kproma,jk)* &
!!$        & (bb(jcs:kproma,jk,iqv) + (tpfac3 - 1._wp)*pqm1(jcs:kproma,jk)) * zrdt
!!$      END DO
!!$    END IF
    !-------------------------------------------------------------
    ! Tendency of tracers
    !-------------------------------------------------------------
!   IF (trlist% anyvdiff /= 0) THEN   ! ECHAM
!     DO 577 jt=1,trlist% ntrac       ! ECHAM
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(ktrac > 0)
        !$ACC LOOP GANG
        DO jt = 1,ktrac
          irhs = itrc_start + jt - 1
!         IF (trlist% ti(jt)% nvdiff /= 1) CYCLE  ! ECHAM
          !$ACC LOOP
          DO jk = 1,klev
            !$ACC LOOP VECTOR
            DO jl = jcs,kproma
              pxtte_vdf(jl,jk,jt) = (bb(jl,jk,irhs)-tpfac2*pxtm1(jl,jk,jt))*zrdt
            ENDDO
          ENDDO
        ENDDO
        !$ACC END PARALLEL

!577  ENDDO
!     END IF

    !----------------------------------------------------------------------------
    ! Update roughness height over open water, then update the grid-box mean
    !----------------------------------------------------------------------------
    IF (idx_wtr<=ksfc_type) THEN  ! water surface exists in the simulation
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jl = 1,kbdim
        pz0m_tile(jl,idx_wtr) = vdiff_config%z0m_oce
      ENDDO
      !$ACC END PARALLEL
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jl = jcs,kproma
        IF(pfrc(jl,idx_wtr).GT.0._wp) THEN
          pz0m_tile(jl,idx_wtr) = tpfac1*SQRT( bb(jl,klev,iu)**2+bb(jl,klev,iv)**2 ) &
                                & *pcfm_tile(jl,idx_wtr)*cchar*rgrav
          pz0m_tile(jl,idx_wtr) = MAX(z0m_min,pz0m_tile(jl,idx_wtr))
        ENDIF
      ENDDO
      !$ACC END PARALLEL
    ENDIF

    ! Compute grid-box mean

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jl = 1,kbdim
      pz0m(jl) = 0._wp
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO jsfc = 1,ksfc_type
      !$ACC LOOP GANG VECTOR
      DO jl = jcs,kproma
        pz0m(jl) = pz0m(jl) + pfrc(jl,jsfc)*pz0m_tile(jl,jsfc)
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !-------------------------------------------------------------------
    ! End GPU data region
    !-------------------------------------------------------------------
    !$ACC WAIT
    !$ACC END DATA


  END SUBROUTINE vdiff_tendencies
  !-------------

END MODULE mo_turb_vdiff
