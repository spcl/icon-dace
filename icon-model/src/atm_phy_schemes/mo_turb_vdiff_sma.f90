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

! Subroutines for computing turbulent exchange coefficients of
! 3D Smagorinsky turbulent scheme.

!------------------
#include "fsel.inc"
!------------------


MODULE mo_turb_vdiff_sma

  USE mo_kind              ,ONLY: wp, i1
  USE mo_convect_tables    ,ONLY: compute_qsat
  USE mo_turb_vdiff_config ,ONLY: t_vdiff_config
  USE mo_turb_vdiff_params ,ONLY: ckap
  USE mo_physical_constants,ONLY: grav, rd, cpd, rd_o_cpd,             &
    &                             vtmpc1, p0ref, rgrav
  USE mo_model_domain      ,ONLY: t_patch
  USE mo_nonhydro_types    ,ONLY: t_nh_metrics
  USE mo_intp_data_strc    ,ONLY: t_int_state
  USE mo_nonhydro_state    ,ONLY: p_nh_state
  USE mo_intp_data_strc    ,ONLY: p_int_state
  USE mo_fortran_tools     ,ONLY: init
  USE mo_impl_constants    ,ONLY: min_rlcell, min_rledge_int, min_rlcell_int, &
    &                             min_rlvert_int
  USE mo_parallel_config   ,ONLY: p_test_run
  USE mo_loopindices       ,ONLY: get_indices_e, get_indices_c
  USE mo_nh_vert_interp_les,ONLY: brunt_vaisala_freq, vert_intp_full2half_cell_3d
  USE mo_index_list        ,ONLY: generate_index_list_batched
  USE mo_intp              ,ONLY: cells2verts_scalar, cells2edges_scalar
  USE mo_sync              ,ONLY: SYNC_E, SYNC_C, SYNC_V, sync_patch_array, &
    &                             sync_patch_array_mult
  USE mo_intp_rbf          ,ONLY: rbf_vec_interpol_vertex, rbf_vec_interpol_edge
  USE mo_impl_constants_grf,ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_nh_testcases_nml  ,ONLY: is_dry_cbl
  USE mo_satad             ,ONLY: spec_humi, sat_pres_water
  USE mo_math_constants    ,ONLY: pi_2, ln2
  USE mo_math_utilities    ,ONLY: tdma_solver_vec
  USE mo_intp_rbf          ,ONLY: rbf_vec_interpol_cell

  IMPLICIT NONE
  PRIVATE
  REAL(wp),         PARAMETER :: z_1by3  = 1._wp/3._wp

  !Parameter for vertical scheme type
  INTEGER, PARAMETER :: iexplicit = 1
  INTEGER, PARAMETER :: iimplicit = 2

  CHARACTER(len=*), PARAMETER :: inmodule = 'mo_sgs_turbulence:'

  PUBLIC :: atm_exchange_coeff3d, diffuse_hori_velocity, diffuse_vert_velocity, &
          & diffuse_scalar

  !Parameters for surface layer parameterizations: From Zeng_etal 1997 J. Clim
  REAL(wp), PARAMETER :: bsm = 5.0_wp  !Businger Stable Momentum
  REAL(wp), PARAMETER :: bum = 16._wp  !Businger Untable Momentum
  REAL(wp), PARAMETER :: bsh = 5.0_wp  !Businger Stable Heat
  REAL(wp), PARAMETER :: buh = 16._wp  !Businger Untable Heat

CONTAINS
  !>
  !! Compute various thermodynamic variables for all (full) vertical levels;
  !! Diagnose PBL extension;
  !! Diagnose wind shear, buoyancy, Ri-number, mixing length, then compute
  !! the turbulent exchange coefficients of momentum, dry static energy,
  !! tracers, TTE, variance of virtual optential temperation at half levels
  !! [1+1/2, klev-1/2].
  !!
  !! Note that
  !! - for all coeffcient arrays, vertical index k in this subroutine
  !!   correspond to interface (half level) k+1/2;
  !! - the exchange coefficient at model top (level 1/2) is zero, thus does
  !!   not need computing;
  !! Therefore a number of variables are defined on what shall be referred to
  !! as "mid-levels", which are defined on levels with indices xmid(k) corre-
  !! sponding to xh(k+1), where x is some quantity at half-level k+1/2.
  !!
  SUBROUTINE atm_exchange_coeff3d( kbdim, nblks_c,                        &! in
                               & klev, klevm1, klevp1,                    &! in
                               & ksfc_type, idx_lnd,                      &! in
                               & p_patch,                                 &! in
                               & pz0m, ptsfc, pfrc,                       &! in
                               & ppsfc,                                   &! in
                               & pghf, pgeof,                             &! in
                               & pum1, pvm1, pwm1,                        &! in
                               & ptm1, ptvm1,                             &! in
                               & pqm1, pxm1,                              &! in
                               & rho,                                     &! in
                               & papm1, paphm1,                           &! in
                               & vdiff_config,                            &! in
                               & pri_tile,                                &! out
                               & pcfm_tile,                               &! out
                               & pcfh_tile,                               &! out
                               & pqsat_tile, pcpt_tile,                   &! out
                               & pcptgz,                                  &! out
                               & pzthvvar, ptottevn,                      &! out
                               & pmixlen,                                 &! out
                               & pcfm, pcfh, pcfv, pcftotte, pcfthv,      &! out
                               & km_c, km_iv, km_ie, km_ic, kh_ic,        &! out
                               & pprfac,                                  &! out
                               & u_vert, v_vert, div_c,                   &! out
                               & inv_rho_ic, w_vert, w_ie,                &! out
                               & vn,                                      &! out
                               & pch_tile,                                &! out
                               & pbn_tile, pbhn_tile, pbm_tile, pbh_tile, &! out
                               & pcsat, pcair                             &! in, optional
                               & )

    ! Arguments

    INTEGER, INTENT(IN) :: nblks_c
    INTEGER, INTENT(IN) :: kbdim
    INTEGER :: nproma
    INTEGER, INTENT(IN) :: klev, klevm1, klevp1
    INTEGER :: nlev, nlevp1
    REAL(wp),INTENT(IN) :: pghf(:,:,:)
    REAL(wp),INTENT(IN) :: pgeof(:,:,:) !< Geopotential above ground (full level) [m2/s2]
    REAL(wp),INTENT(IN) :: pxm1(:,:,:)
    REAL(wp),INTENT(IN) :: ptvm1(:,:,:)
    REAL(wp),INTENT(IN) :: pqm1(:,:,:)
    REAL(wp),INTENT(IN) :: ptm1(:,:,:), rho(:,:,:)
    REAL(wp),INTENT(IN) :: papm1(:,:,:),  paphm1(:,:,:)
    REAL(wp),INTENT(IN) :: pum1(:,:,:),  pvm1(:,:,:), pwm1(:,:,:)

    TYPE(t_vdiff_config), INTENT(IN) :: vdiff_config

    REAL(wp),INTENT(OUT) :: ptottevn(:,:,:) !< TTE at intermediate time step
    REAL(wp),INTENT(OUT) :: pcftotte(:,:,:) !< exchange coeff. for TTE
    REAL(wp),INTENT(OUT) :: pcfthv  (:,:,:) !< exchange coeff. for var. of theta_v
    REAL(wp),INTENT(OUT) :: pcfm    (:,:,:) !< exchange coeff. for u, v
    REAL(wp),INTENT(OUT) :: pcfh    (:,:,:) !< exchange coeff. for cptgz and tracers
    REAL(wp),INTENT(OUT) :: pcfv    (:,:,:) !< exchange coeff. for variance of qx
    REAL(wp),INTENT(OUT) :: pzthvvar(:,:,:) !< variance of theta_v at interm. step
    REAL(wp),INTENT(OUT) :: pcptgz  (:,:,:) !< dry static energy
    REAL(wp),INTENT(OUT) :: pprfac  (:,:,:) !< prefactor for the exchange coeff.
    REAL(wp),INTENT(OUT) :: pmixlen (:,:,:) !< prefactor for the exchange coeff.

    ! Local variables
    ! - Variables defined at full levels

    REAL(wp) :: ztheta (kbdim,klev,nblks_c)  !< potential temperature

    ! - Variables defined at mid-levels

    REAL(wp) :: zdgmid !< geopotential height difference between two full levels
    REAL(wp) :: ztvmid

    TYPE(t_patch)   ,TARGET ,INTENT(IN)   :: p_patch
    TYPE(t_nh_metrics) ,POINTER :: p_nh_metrics
    TYPE(t_int_state)  ,POINTER :: p_int         !< interpolation state

    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end
    INTEGER :: jc, jb, je

  !Variables for the module
    REAL(wp), INTENT(OUT), DIMENSION(:,:,:)   :: km_c
    REAL(wp), INTENT(OUT), DIMENSION(:,:,:)   :: km_iv
    REAL(wp), INTENT(OUT), DIMENSION(:,:,:)   :: km_ie
    REAL(wp), INTENT(OUT), DIMENSION(:,:,:)   :: kh_ic, km_ic
    REAL(wp), INTENT(OUT), DIMENSION(:,:,:)   :: u_vert, v_vert
    REAL(wp), INTENT(OUT), DIMENSION(:,:,:)   :: div_c
    REAL(wp), INTENT(OUT), DIMENSION(:,:,:)   :: w_vert
    REAL(wp), INTENT(OUT), DIMENSION(:,:,:)   :: w_ie
    REAL(wp), INTENT(OUT), DIMENSION(:,:,:)   :: inv_rho_ic
    REAL(wp), INTENT(OUT), DIMENSION(:,:,:)   :: vn !< normal wind vector

    REAL(wp), DIMENSION(kbdim,klev,nblks_c)   :: theta_v
    REAL(wp), DIMENSION(kbdim,klevp1,nblks_c) :: bruvais

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: vn_ie, vt_ie, mech_prod, &
                                                 shear, div_of_stress

    INTEGER, DIMENSION(:,:,:), POINTER :: ividx, ivblk, iecidx, iecblk, ieidx, ieblk

    REAL(wp) :: vn_vert1, vn_vert2, vn_vert3, vn_vert4
    REAL(wp) :: vt_vert1, vt_vert2, vt_vert3, vt_vert4
    REAL(wp) :: w_full_c1, w_full_c2, w_full_v1, w_full_v2
    REAL(wp) :: D_11, D_12, D_13, D_22, D_23, D_33

    ! - surface variables
    INTEGER, INTENT(IN) :: ksfc_type, idx_lnd

    REAL(wp),INTENT(IN) :: pz0m  (:,:,:) !< aerodynamic roughness length
    REAL(wp),INTENT(IN) :: ptsfc (:,:,:) !< temp. at surface
    REAL(wp),INTENT(IN) :: pfrc  (:,:,:) !< fraction of the grid box occupied
    REAL(wp),INTENT(IN) :: ppsfc (:,:)  !< surface pressure

    ! optional arguments for use with jsbach
    REAL(wp),OPTIONAL,INTENT(IN) :: pcsat  (:,:)  !< area fraction with wet land surface
    REAL(wp),OPTIONAL,INTENT(IN) :: pcair  (:,:)  !< area fraction with wet land surface (air)

    REAL(wp),INTENT(OUT) :: pqsat_tile (:,:,:)!< saturation specific humidity
    REAL(wp),INTENT(OUT) :: pcpt_tile (:,:,:) !< dry static energy
    REAL(wp),INTENT(OUT) :: pcfm_tile (:,:,:) !< exchange coeff. of momentum,
                                                                !< for each type of surface
    REAL(wp),INTENT(OUT) :: pcfh_tile (:,:,:) !< exchange coeff. of heat and
                                                                !<  vapor for each surface type
    REAL(wp),INTENT(OUT) :: pbn_tile  (:,:,:) !< for diagnostics
    REAL(wp),INTENT(OUT) :: pbhn_tile (:,:,:) !< for diagnostics
    REAL(wp),INTENT(OUT) :: pbm_tile  (:,:,:) !< for diagnostics
    REAL(wp),INTENT(OUT) :: pbh_tile  (:,:,:) !< for diagnostics
    REAL(wp),INTENT(OUT) :: pch_tile  (:,:,:) !< for TTE boundary condition
    REAL(wp),INTENT(OUT) :: pri_tile  (:,:,:) !< Richardson number for diagnostics

    REAL(wp) :: zqts
    REAL(wp) :: zvn1, zvn2
    INTEGER(i1)::pfrc_test(kbdim,ksfc_type) !< integer mask to pass to CUB (can be removed later)
    INTEGER  :: loidx (kbdim,ksfc_type)     !< counter for masks
    INTEGER  :: is    (ksfc_type)           !< counter for masks
    INTEGER  :: jsfc, jls, js
    INTEGER  :: jcn,jbn                     !< jc and jb of neighbor cells sharing an edge je
    REAL(wp),parameter :: zcons17 = 1._wp / ckap**2

    INTEGER  :: itr
    REAL(wp) :: zrough, theta_sfc, qv_s, mwind, z_mc, RIB, tcn_mom, tcn_heat, &
                shfl, lhfl, bflx1, ustar, obukhov_length, inv_bus_mom
    REAL(wp) :: tch
    REAL(wp) :: tcm
    REAL(wp) :: zthetavmid

    ! Required for intermediate usage of the density at interface levels
    REAL(wp) :: rho_ic(kbdim,klevp1,nblks_c)

    ! - 1D variables and scalars

    INTEGER  :: jg, jk, jl
    REAL(wp) :: zrdp
    REAL(wp) :: zsdep1
    REAL(wp) :: zsdep2

    ! to prevent floating-point arithmetic inconsistencies later in
    ! the interpolation to u 10m and 2m T/T_d: has been 0.01 before
    ! (Cray FP instead of IEEE 754 FP format)
    REAL(wp) :: zepsec = 0.028_wp

    ! Shortcuts to components of aes_vdf_config
    !
    REAL(wp) :: fsl, min_sfc_wind, km_min, turb_prandtl, rturb_prandtl
    !

    jg = p_patch%id
    p_nh_metrics => p_nh_state(jg)%metrics
    p_int        => p_int_state(jg)

    nproma = kbdim
    nlev = klev; nlevp1 = klevp1

    fsl           = vdiff_config%fsl
    min_sfc_wind  = vdiff_config%min_sfc_wind
    km_min        = vdiff_config%km_min
    turb_prandtl  = vdiff_config%turb_prandtl
    rturb_prandtl = vdiff_config%rturb_prandtl

    rl_start   = 1
    rl_end     = min_rlcell
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    !$ACC DATA &
    !---- Argument arrays - intent(in)
    !$ACC   PRESENT(pghf, pgeof, pxm1, ptvm1, pqm1, pwm1, ptm1, rho, papm1, paphm1) &
    !$ACC   PRESENT(pz0m, ptsfc, pfrc, ppsfc, pcsat, pcair) &
    !---- Argument arrays - intent(inout)
    !$ACC   PRESENT(pum1, pvm1, p_patch) &
    !---- Argument arrays - intent(out)
    !$ACC   PRESENT(ptottevn, pcftotte, pcfthv, pcfm, pcfh, pcfv, pzthvvar, pcptgz, pprfac, pmixlen) &
    !$ACC   PRESENT(pqsat_tile, pcpt_tile, pcfm_tile, pcfh_tile, pbn_tile, pbhn_tile, pbm_tile) &
    !$ACC   PRESENT(pbh_tile, pch_tile, pri_tile) &
    !$ACC   PRESENT(km_c, km_iv, km_ie, kh_ic, km_ic, vn) &
    !$ACC   PRESENT(u_vert, v_vert, w_vert, inv_rho_ic, div_c, w_ie) &
    !---- Argument arrays - Module Variables
    !$ACC   PRESENT(p_nh_metrics, p_int) &
    !$ACC   CREATE(loidx, pfrc_test, ztheta, is) &
    !$ACC   CREATE(theta_v, bruvais, rho_ic)



!#########################################################################
!## initialize
!#########################################################################

!$OMP PARALLEL
    CALL init(km_iv)
    CALL init(km_c)
    CALL init(km_ie)
    CALL init(kh_ic)
    CALL init(km_ic)
    CALL init(vn)

    IF(p_test_run)THEN
      CALL init(u_vert(:,:,:))
      CALL init(v_vert(:,:,:))
      CALL init(w_vert(:,:,:))
    END IF
!$OMP END PARALLEL

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        ptottevn(jc,:,jb) = 10._wp  ! for vdiff_up
        pzthvvar(jc,:,jb) = 1._wp   ! for vdiff_up

        pcftotte(jc,:,jb) = 10._wp  ! K TEE
        pcfthv(jc,:,jb)   = 0._wp   ! K theta_v
        pcfv(jc,:,jb)     = 0._wp   ! K qx
        pmixlen(jc,:,jb)  = 0._wp   ! K qx

        pcfm_tile(jc,jb,:) = 0._wp
        pcfh_tile(jc,jb,:) = 0._wp

        pqsat_tile(jc,jb,:)= 0._wp
        pri_tile(jc,jb,:)  = 0._wp
        pch_tile(jc,jb,:)  = 0._wp
        pcpt_tile(jc,jb,:) = 0._wp
      END DO
      !$ACC END PARALLEL

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   ! Here pum1 and pvm1 are assumed to be synchronized already, see interface_iconam_aes.

    rl_start   = 2
    rl_end     = min_rledge_int
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,jcn,jbn,zvn1,zvn2)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
          jcn  =   p_patch%edges%cell_idx(je,jb,1)
          jbn  =   p_patch%edges%cell_blk(je,jb,1)
          zvn1 =   pum1(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(je,jb,1)%v1 &
            &    + pvm1(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(je,jb,1)%v2
          !
          jcn  =   p_patch%edges%cell_idx(je,jb,2)
          jbn  =   p_patch%edges%cell_blk(je,jb,2)
          zvn2 =   pum1(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(je,jb,2)%v1 &
            &    + pvm1(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(je,jb,2)%v2
          !
          vn(je,jk,jb) = p_int%c_lin_e(je,1,jb)*zvn1 &
            &          + p_int%c_lin_e(je,2,jb)*zvn2
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO
!$OMP END PARALLEL

   !$ACC WAIT
   CALL sync_patch_array(SYNC_E, p_patch, vn)

!#########################################################################
!## variables for TTE scheme and JSBACH LSM
!#########################################################################

    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jl,jls,js,i_startidx,i_endidx,is,jsfc,zqts,zrough, theta_sfc,   &
!$OMP            qv_s, mwind, z_mc, RIB, tcn_mom, tcn_heat, itr, shfl, lhfl,        &
!$OMP            bflx1, ustar, obukhov_length, inv_bus_mom, loidx, pfrc_test, zrdp, &
!$OMP            zsdep1, zsdep2, ztvmid, zdgmid, zthetavmid, tch, tcm)

  DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)


    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jl = i_startidx, i_endidx
      DO jk = 1, klev
       IF (jk .le. klevm1) THEN
        zdgmid = pghf(jl,jk,jb) - pghf(jl,jk+1,jb)

        ! interpolation coefficients
        zrdp   = 1._wp/(papm1(jl,jk,jb) - papm1(jl,jk+1,jb))
        zsdep1 = (paphm1(jl,jk+1,jb)- papm1(jl,jk+1,jb))*zrdp
        zsdep2 = (papm1(jl,jk,jb)  - paphm1(jl,jk+1,jb))*zrdp

        ztvmid = zsdep1*ptvm1(jl,jk,jb)+zsdep2*ptvm1(jl,jk+1,jb)

        ! Virtual dry static energy
        pcptgz(jl,jk,jb) = pgeof(jl,jk,jb)+ptm1(jl,jk,jb)*cpd !+(cpv-cpd)*pqm1(jl,jk,jb))

        ! Potential temperature
        ztheta(jl,jk,jb) = ptm1(jl,jk,jb)*(p0ref/papm1(jl,jk,jb))**rd_o_cpd

        ! Air density at mid levels, p/(Tv*R)/dz = air density/dz, and the prefactor
        ! that will be multiplied later to the exchange coeffcients to build a linear
        ! algebraic equation set.
        pprfac(jl,jk,jb) = paphm1(jl,jk+1,jb)/(ztvmid*rd*zdgmid)

       ELSE

      ! dry static energy pcpt_tile
      !
      pcptgz(jl,klev,jb) = pgeof(jl,klev,jb)+ptm1(jl,klev,jb)*cpd!+(cpv-cpd)*pqm1(jl,klev,jb))

      ! Potential temperature
      ztheta(jl,klev,jb) = ptm1(jl,klev,jb)*(p0ref/papm1(jl,klev,jb))**rd_o_cpd

      ! The prefactor (= air density) that will be multiplied to the exchange
      ! coefficients when building the linear algebraic equations. The density here is
      ! computed using air temperature of the lowest model level at time step n-1.
      pprfac(jl,klev,jb) = ppsfc(jl,jb)                                       &
                           &  /( rd*ptm1(jl,klev,jb)                             &
                           &  *(1._wp+vtmpc1*pqm1(jl,klev,jb)-pxm1(jl,klev,jb)) )
       ENDIF
      ENDDO
    END DO
    !$ACC END PARALLEL


    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
    pcfm(:,klev,jb) = 0._wp
    pcfh(:,klev,jb) = 0._wp
    !$ACC END KERNELS

    !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1)
    DO jsfc = 1,ksfc_type
      DO jl = i_startidx, i_endidx
        pfrc_test(jl,jsfc) = MERGE(1_i1, 0_i1, pfrc(jl,jb,jsfc) > 0.0_wp)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    CALL generate_index_list_batched(pfrc_test, loidx, i_startidx, i_endidx, is, 1)
    !$ACC UPDATE HOST(is) ASYNC(1)
    !$ACC WAIT(1)

    DO jsfc = 1, ksfc_type
      CALL compute_qsat( kbdim, is(jsfc), loidx(:,jsfc), ppsfc(:,jb), ptsfc(:,jb,jsfc), pqsat_tile(:,jb,jsfc) )

     ! loop over mask only
     !

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jls = 1, is(jsfc)
        js=loidx(jls,jsfc)

        ! dry static energy pcpt_tile
        !
        IF(jsfc == idx_lnd ) THEN
          zqts = pcsat(js,jb) * pqsat_tile(js,jb,jsfc) + (1._wp - pcair(js,jb))*pqm1(js,klev,jb) ! q_total at land surface
        ELSE
          zqts = pqsat_tile(js,jb,jsfc)                                              ! q_total at non-land surface
        END IF

        pcpt_tile(js,jb,jsfc) = ptsfc(js,jb,jsfc) * cpd ! (cpd + (cpv - cpd) * zqts)

        !Surface roughness length
        IF ( pz0m(js,jb,jsfc) .GT. 0.0_wp ) THEN
          zrough = pz0m(js,jb,jsfc)
        ELSE
          zrough = 1.E-3_wp
        END IF

        !Get surface pot. temperature and humidity
        theta_sfc = ptsfc(js,jb,jsfc) / EXP( rd_o_cpd*LOG(ppsfc(js,jb)/p0ref) )
        qv_s    = spec_humi(sat_pres_water(ptsfc(js,jb,jsfc)),ppsfc(js,jb))


        mwind = MAX( min_sfc_wind, SQRT(pum1(js,klev,jb)**2+pvm1(js,klev,jb)**2) )

        !Z height to be used as a reference height in surface layer
        z_mc   = p_nh_metrics%z_mc(js,klev,jb) - p_nh_metrics%z_ifc(js,klevp1,jb)

        !First guess for tch and tcm using bulk approach
        RIB = grav * (ztheta(js,klev,jb)-theta_sfc) * (z_mc-zrough) / (theta_sfc * mwind**2)
        tcn_mom = (ckap/LOG(z_mc/zrough))**2
        tcm     = tcn_mom * stability_function_mom(RIB,z_mc/zrough,tcn_mom)

        tcn_heat        = ckap**2/(LOG(z_mc/zrough)*LOG(z_mc/zrough))
        tch             = tcn_heat * stability_function_heat(RIB,z_mc/zrough,tcn_heat)

        !now iterate
        DO itr = 1 , 5
           shfl = tch*mwind*(theta_sfc-ztheta(js,klev,jb))
           lhfl = tch*mwind*(qv_s-pqm1(js,klev,jb))
           bflx1= shfl + vtmpc1 * theta_sfc * lhfl
           ustar= SQRT(tcm)*mwind

           obukhov_length = -ustar**3 * theta_sfc * rgrav / (ckap * bflx1)

           inv_bus_mom = 1._wp / businger_mom(zrough,z_mc,obukhov_length)
           tch         = inv_bus_mom / businger_heat(zrough,z_mc,obukhov_length)
           tcm         = inv_bus_mom * inv_bus_mom
        END DO

        pcfm_tile(js,jb,jsfc) = tcm*mwind
        pcfh_tile(js,jb,jsfc) = tch*mwind
        pch_tile (js,jb,jsfc) = tch

        pcfm(js,klev,jb) = pcfm(js,klev,jb) + pfrc(js,jb,jsfc)*pcfm_tile(js,jb,jsfc)
        pcfh(js,klev,jb) = pcfh(js,klev,jb) + pfrc(js,jb,jsfc)*pcfh_tile(js,jb,jsfc)


        pbn_tile(js,jb,jsfc) = ckap / MAX( zepsec, sqrt(tcn_mom) )
        pbhn_tile(js,jb,jsfc)= ckap / MAX( zepsec, sqrt(tcn_heat) )
        pbm_tile(js,jb,jsfc) = MAX( zepsec, sqrt(pcfm_tile(js,jb,jsfc) * tch*zcons17/ (tcn_mom*mwind)) )
        pbh_tile(js,jb,jsfc) = MAX( zepsec, tch/pbm_tile(js,jb,jsfc)*zcons17)
        pbm_tile(js,jb,jsfc) = 1._wp / pbm_tile(js,jb,jsfc)
        pbh_tile(js,jb,jsfc) = 1._wp / pbh_tile(js,jb,jsfc)

        zthetavmid           = fsl*ptvm1(js,nlev,jb)*(p0ref/papm1(js,nlev,jb))**rd_o_cpd + &
                               & (1._wp-fsl) * theta_sfc * (1._wp+vtmpc1*zqts)

        pri_tile(js,jb,jsfc) = pghf(js,klev,jb) * grav *                                   &
                             & ( ptvm1(js,nlev,jb)*(p0ref/papm1(js,nlev,jb))**rd_o_cpd -   &
                             &   theta_sfc * (1._wp+vtmpc1*zqts)                       ) / &
                             & ( zthetavmid * mwind )
      END DO !jls
      !$ACC END PARALLEL
    END DO !jsfc


  END DO !jb
!$OMP END DO
!$OMP END PARALLEL



!#########################################################################
!## Convert temperature to potential temperature: all routines within
!## use theta.
!#########################################################################


    rl_start   = 3
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        theta_v(jc,1:nlev,jb) = ptvm1(jc,1:nlev,jb)*(p0ref/papm1(jc,1:nlev,jb))**rd_o_cpd
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    !Get rho at interfaces to be used later
    CALL vert_intp_full2half_cell_3d(p_patch, p_nh_metrics, rho, rho_ic, &
                                     1, min_rlcell_int-2, lacc=.TRUE.)
    
    ! Compute the Brunt Vaisala frequency where theta_v was defined
                                     
    CALL brunt_vaisala_freq(p_patch, p_nh_metrics, kbdim, theta_v, bruvais, &
                            opt_rlstart=3, lacc=.TRUE.)
    
    
!#########################################################################
!## Smagorinsky_model
  !!------------------------------------------------------------------------
  !! Computes the sgs viscosity and diffusivity using Smagorinsky model
  !! \tau_ij = KD_ij where D_ij = du_i/dx_j + du_j/dx_i
  !!
  !! and  K = cs * \Delta**2 * D / sqrt(2), where D = sqrt(D_ijD_ij)
  !!
  !! and D**2 = D_11**2 + D_22**2 + D_33**2 + 2D_12**2 + 2D_13**2 + 2D_23**2
  !!
  !! where, D_11 = 2 * du_1/dx_1
  !!        D_22 = 2 * d_u2/dx_2
  !!        D_33 = 2 * d_u3/dx_3
  !!        D_12 = du_1/dx_2 + du_2/dx_1
  !!        D_13 = du_1/dx_3 + du_3/dx_1
  !!        D_23 = du_2/dx_3 + du_3/dx_2
  !! For triangles: 1=normal, 2=tangential, and 3 = z directions
  !!------------------------------------------------------------------------

    ALLOCATE( vn_ie(nproma,p_patch%nlevp1,p_patch%nblks_e)      &
             ,vt_ie(nproma,p_patch%nlevp1,p_patch%nblks_e)      &
             ,shear(nproma,p_patch%nlev,p_patch%nblks_e)        &
             ,div_of_stress(nproma,p_patch%nlev,p_patch%nblks_e)&
             ,mech_prod(nproma,nlevp1,p_patch%nblks_c)          &
            )
    !$ACC DATA CREATE(vn_ie, vt_ie, shear, div_of_stress, mech_prod)

    IF(p_test_run)THEN
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
      kh_ic(:,:,:) = 0._wp
      km_ic(:,:,:) = 0._wp
      !$ACC END KERNELS
    END IF

    !--------------------------------------------------------------------------
    !1) Interpolate velocities at desired locations- mostly around the quadrilateral
    !
    !<It assumes that prog values are all synced at this stage while diag values might not>
    !--------------------------------------------------------------------------

    ! Here pwm1 is assumed to be synchronized already, see interface_iconam_aes.

    CALL cells2verts_scalar(pwm1, p_patch, p_int%cells_aw_verts, w_vert,                   &
                            opt_rlend=min_rlvert_int, opt_acc_async=.TRUE.)
    CALL cells2edges_scalar(pwm1, p_patch, p_int%c_lin_e, w_ie, opt_rlend=min_rledge_int-2,&
                            lacc=.TRUE.)

    ! RBF reconstruction of velocity at vertices: include halos
    CALL rbf_vec_interpol_vertex(vn, p_patch, p_int, u_vert, v_vert, &
                                 opt_rlend=min_rlvert_int, opt_acc_async=.TRUE. )

    !$ACC WAIT

    !sync them
    CALL sync_patch_array_mult(SYNC_V, p_patch, 3, w_vert, u_vert, v_vert)

    !Get vn at interfaces and then get vt at interfaces
    !Boundary values are extrapolated like dynamics although
    !they are not required in current implementation

!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)
    rl_start   = 2
    rl_end     = min_rledge_int-3
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 2, nlev
#else
      DO jk = 2, nlev
        DO je = i_startidx, i_endidx
#endif
          vn_ie(je,jk,jb) = p_nh_metrics%wgtfac_e(je,jk,jb) * vn(je,jk,jb) +            &
                            ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * vn(je,jk-1,jb)
        END DO
      END DO
      !$ACC END PARALLEL
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO je = i_startidx, i_endidx
        vn_ie(je,1,jb)      = p_nh_metrics%wgtfacq1_e(je,1,jb) * vn(je,1,jb) +          &
                              p_nh_metrics%wgtfacq1_e(je,2,jb) * vn(je,2,jb) +          &
                              p_nh_metrics%wgtfacq1_e(je,3,jb) * vn(je,3,jb)

        vn_ie(je,nlevp1,jb) = p_nh_metrics%wgtfacq_e(je,1,jb) * vn(je,nlev,jb)   +      &
                              p_nh_metrics%wgtfacq_e(je,2,jb) * vn(je,nlev-1,jb) +      &
                              p_nh_metrics%wgtfacq_e(je,3,jb) * vn(je,nlev-2,jb)
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO
!$OMP END PARALLEL

    CALL rbf_vec_interpol_edge(vn_ie, p_patch, p_int, vt_ie, opt_rlstart=3, &
                               opt_rlend=min_rledge_int-2, opt_acc_async=.TRUE.)

    !--------------------------------------------------------------------------
    !2) Compute horizontal strain rate tensor at full levels
    !--------------------------------------------------------------------------
    ividx => p_patch%edges%vertex_idx
    ivblk => p_patch%edges%vertex_blk

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk

    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    rl_start   = 4
    rl_end     = min_rledge_int-2
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,vn_vert1,vn_vert2,vn_vert3,vn_vert4,  &
!$OMP            vt_vert1,vt_vert2,vt_vert3,vt_vert4,w_full_c1,w_full_c2,w_full_v1,  &
!$OMP            w_full_v2,D_11,D_12,D_13,D_22,D_23,D_33)
    DO jb = i_startblk,i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR TILE(32, 4)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
#endif

          vn_vert1 = u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,1)%v1 +   &
                     v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,1)%v2

          vn_vert2 = u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,2)%v1 +   &
                     v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,2)%v2

          vn_vert3 = u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,3)%v1 +   &
                     v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,3)%v2

          vn_vert4 = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,4)%v1 +   &
                     v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,4)%v2

          vt_vert1 = u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1))     *   &
                     p_patch%edges%dual_normal_vert(je,jb,1)%v1   +   &
                     v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1))     *   &
                     p_patch%edges%dual_normal_vert(je,jb,1)%v2

          vt_vert2 = u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2))     *   &
                     p_patch%edges%dual_normal_vert(je,jb,2)%v1   +   &
                     v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2))     *   &
                     p_patch%edges%dual_normal_vert(je,jb,2)%v2

          vt_vert3 = u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))     *   &
                     p_patch%edges%dual_normal_vert(je,jb,3)%v1   +   &
                     v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))     *   &
                     p_patch%edges%dual_normal_vert(je,jb,3)%v2

          vt_vert4 = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))     *   &
                     p_patch%edges%dual_normal_vert(je,jb,4)%v1   +   &
                     v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))     *   &
                     p_patch%edges%dual_normal_vert(je,jb,4)%v2

          ! W at full levels
          w_full_c1  = 0.5_wp *                                              &
                       ( pwm1(iecidx(je,jb,1),jk,iecblk(je,jb,1)) +   &
                         pwm1(iecidx(je,jb,1),jk+1,iecblk(je,jb,1)) )

          w_full_c2  = 0.5_wp *                                              &
                       ( pwm1(iecidx(je,jb,2),jk,iecblk(je,jb,2)) +   &
                         pwm1(iecidx(je,jb,2),jk+1,iecblk(je,jb,2)) )

          ! W at full levels vertices from w at vertices at interface levels
          w_full_v1  = 0.5_wp *                                              &
                       ( w_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) +          &
                         w_vert(ividx(je,jb,1),jk+1,ivblk(je,jb,1)) )

          w_full_v2  = 0.5_wp *                                              &
                       ( w_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) +          &
                         w_vert(ividx(je,jb,2),jk+1,ivblk(je,jb,2)) )


          ! Strain rates at edge center
          D_11       = 2._wp * ( vn_vert4 - vn_vert3 ) *                     &
                       p_patch%edges%inv_vert_vert_length(je,jb)

          D_12       = p_patch%edges%tangent_orientation(je,jb) *            &
                       ( vn_vert2 - vn_vert1 ) *                             &
                       p_patch%edges%inv_primal_edge_length(je,jb) +         &
                       ( vt_vert4-vt_vert3 ) *                               &
                       p_patch%edges%inv_vert_vert_length(je,jb)

          D_13       = ( vn_ie(je,jk,jb) - vn_ie(je,jk+1,jb) ) *             &
                       p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)  +           &
                       ( w_full_c2 - w_full_c1 ) *                           &
                       p_patch%edges%inv_dual_edge_length(je,jb)

          D_22       = 2._wp * ( vt_vert2-vt_vert1 ) *                       &
                       p_patch%edges%tangent_orientation(je,jb) *            &
                       p_patch%edges%inv_primal_edge_length(je,jb)

          D_23       = ( vt_ie(je,jk,jb) - vt_ie(je,jk+1,jb) ) *             &
                       p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)  +           &
                       p_patch%edges%tangent_orientation(je,jb) *            &
                       ( w_full_v2 - w_full_v1 ) *                           &
                       p_patch%edges%inv_primal_edge_length(je,jb)

          D_33       = 2._wp * ( w_ie(je,jk,jb) - w_ie(je,jk+1,jb) ) *       &
                       p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)

          ! Mechanical prod is half of this value divided by km
          shear(je,jk,jb) = D_11**2 + D_22**2 + D_33**2 +                    &
                            2._wp * ( D_12**2 + D_13**2 + D_23**2 )

          ! calculate divergence to get the deviatoric part of stress tensor in
          ! diffusion: D_11 - 1/3 * (D_11 + D_22 + D_33)
          div_of_stress(je,jk,jb) = 0.5_wp * ( D_11 + D_22 + D_33 )
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ENDDO
!$OMP END DO
!$OMP END PARALLEL


    !Interpolate div(stress) from edge to cell-scalar, incl. halo for its use in hor diffusion
    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int-1
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,      &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
#endif
          div_c(jc,jk,jb) =                                                                       &
                  ( div_of_stress(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%e_bln_c_s(jc,1,jb) +  &
                    div_of_stress(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%e_bln_c_s(jc,2,jb) +  &
                    div_of_stress(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%e_bln_c_s(jc,3,jb) )
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO
!$OMP END PARALLEL


    ! Interpolate mech. production term from mid level edge to interface level cell
    ! except top and bottom boundaries.
    rl_start   = 3
    rl_end     = min_rlcell_int-1
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,      &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 2, nlev
#else
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
#endif

          mech_prod(jc,jk,jb) = p_nh_metrics%wgtfac_c(jc,jk,jb) * (                      &
                      shear(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%e_bln_c_s(jc,1,jb)   +      &
                      shear(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%e_bln_c_s(jc,2,jb)   +      &
                      shear(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%e_bln_c_s(jc,3,jb) ) +      &
                      ( 1._wp - p_nh_metrics%wgtfac_c(jc,jk,jb) ) * (                             &
                      shear(ieidx(jc,jb,1),jk-1,ieblk(jc,jb,1)) * p_int%e_bln_c_s(jc,1,jb) +      &
                      shear(ieidx(jc,jb,2),jk-1,ieblk(jc,jb,2)) * p_int%e_bln_c_s(jc,2,jb) +      &
                      shear(ieidx(jc,jb,3),jk-1,ieblk(jc,jb,3)) * p_int%e_bln_c_s(jc,3,jb) )
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO
!$OMP END PARALLEL

    !--------------------------------------------------------------------------
    ! 3) Classical Smagorinsky model with stability correction due to Lilly 1962
    !    at interface cell centers. At this point mech_prod is twice the actual
    !    mechanical production term.
    !--------------------------------------------------------------------------
    ! MP = Mechanical prod term calculated above
    ! visc = mixing_length_sq * SQRT(MP/2) * SQRT(1-Ri/Pr) where
    ! Ri = (g / theta) * d_theta_dz / (MP/2), where Brunt_vaisala_freq (byncy prod term/kh)
    !    = (g / theta) * d_theta_dz.
    ! After simplification: visc = mixing_length_sq/SQRT(2) * SQRT[MP/2 - (Brunt_vaisala_frq/Pr)]
    ! Note that the factor SQRT(2) with mixing_length_sq is considered into the Smag constant
    rl_start   = 3
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 2 , nlev
#else
      DO jk = 2 , nlev
        DO jc = i_startidx, i_endidx
#endif
          kh_ic(jc,jk,jb) = rho_ic(jc,jk,jb) * rturb_prandtl *                           &
                            p_nh_metrics%mixing_length_sq(jc,jk,jb) *                    &
                            SQRT( MAX( 0._wp, mech_prod(jc,jk,jb) * 0.5_wp -             &
                            rturb_prandtl * bruvais(jc,jk,jb) ) )
          km_ic(jc,jk,jb) = kh_ic(jc,jk,jb) * turb_prandtl
        END DO
      END DO
      !$ACC END PARALLEL
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        kh_ic(jc,1,jb)      = kh_ic(jc,2,jb)
        kh_ic(jc,nlevp1,jb) = kh_ic(jc,nlev,jb)
        km_ic(jc,1,jb)      = km_ic(jc,2,jb)
        km_ic(jc,nlevp1,jb) = km_ic(jc,nlev,jb)
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO
!$OMP END PARALLEL

    !$ACC WAIT

    CALL sync_patch_array(SYNC_C, p_patch, kh_ic)
    CALL sync_patch_array(SYNC_C, p_patch, km_ic)

    !--------------------------------------------------------------------------
    !4) Interpolate difusivity (viscosity) to different locations: calculate them for
    !   halos also because they will be used later in diffusion
    !--------------------------------------------------------------------------

    !4a) visc at cell center
    rl_start = grf_bdywidth_c
    rl_end   = min_rlcell_int-1
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1 , nlev
#else
      DO jk = 1 , nlev
        DO jc = i_startidx, i_endidx
#endif
          km_c(jc,jk,jb) = MAX( km_min,                                   &
                                ( kh_ic(jc,jk,jb) + kh_ic(jc,jk+1,jb) ) * &
                                0.5_wp * turb_prandtl )
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !4b) visc at vertices
    CALL cells2verts_scalar(kh_ic, p_patch, p_int%cells_aw_verts, km_iv, &
                            opt_rlstart=5, opt_rlend=min_rlvert_int-1,   &
                            opt_acc_async=.TRUE.)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(3)
    DO jb = 1, SIZE(km_iv, 3)
      DO jk = 1, SIZE(km_iv, 2)
        DO jc = 1, SIZE(km_iv, 1)
          km_iv(jc,jk,jb) = MAX( km_min,  km_iv(jc,jk,jb) * turb_prandtl )
        END DO
      END DO
    END DO
    !$ACC END PARALLEL

    !4c) Now calculate visc at half levels at edge
    CALL cells2edges_scalar(kh_ic, p_patch, p_int%c_lin_e, km_ie,                   &
                            opt_rlstart=grf_bdywidth_e, opt_rlend=min_rledge_int-1, &
                            lacc=.TRUE.)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(3)
    DO jb = 1,  SIZE(km_ie, 3)
      DO jk = 1, SIZE(km_ie, 2)
        DO jc = 1, SIZE(km_ie, 1)
          km_ie(jc,jk,jb) = MAX( km_min,  km_ie(jc,jk,jb) * turb_prandtl )
        END DO
      END DO
    END DO
    !$ACC END PARALLEL

    !4d)Get visc at the center on interface level
!    prm_diag%tkvm = MAX( km_min, prm_diag%tkvh * turb_prandtl )

    rl_start   = 1
    rl_end     = min_rlcell
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, nlev-1
#else
      DO jk = 1, nlev-1
        DO jc = i_startidx, i_endidx
#endif
          ! Since the TTE scheme provides the turbulent diffusion coefficients
          ! excluding the density, the same should apply for the Smagorinsky
          ! scheme. Hence, km and kh are devided by rho at this point. The
          ! inverse of the density is saved, as it is required in
          ! diffuse_vert_velocity.
          inv_rho_ic(jc,jk+1,jb) = 1._wp / rho_ic(jc,jk+1,jb)
          pcfh(jc,jk,jb) = kh_ic(jc,jk+1,jb) * inv_rho_ic(jc,jk+1,jb)
          pcfm(jc,jk,jb) = km_ic(jc,jk+1,jb) * inv_rho_ic(jc,jk+1,jb)
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO
!$OMP END PARALLEL

    !$ACC WAIT

    !$ACC END DATA
    !$ACC END DATA

  NULLIFY(p_nh_metrics)
  NULLIFY(p_int)

  END SUBROUTINE atm_exchange_coeff3d
  !-------------
  !>
  !!

  !-------------------------------------------------------------------------------------
  !>
  !! diffuse_hori_velocity
  !!------------------------------------------------------------------------
  !! Calculate the SGS diffusion term for normal velocity component
  !! - Uses the forward Euler time scheme in time split (sequential) manner adopted
  !!   in the NH version.
  !! - Option to switch on implicit scheme in vertical
  !!
  !! d_vn/d_t =  d_tau_11/d_x1 + d_tau_12/d_x2 + d_tau_13/d_x3
  !!------------------------------------------------------------------------

  SUBROUTINE diffuse_hori_velocity( nproma,                &
                                  & p_patch,               &
                                  & km_c, km_iv,           &
                                  & u_vert, v_vert, div_c, &
                                  & rho, vn,               &
                                  & ddt_u, ddt_v)

    INTEGER,INTENT(in) :: nproma
    TYPE(t_patch), TARGET, INTENT(in) :: p_patch      !< single patch
    TYPE(t_nh_metrics) ,POINTER :: p_nh_metrics
    TYPE(t_int_state)  ,POINTER :: p_int            !< interpolation state
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: km_c
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev+1,p_patch%nblks_v) :: km_iv
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_v) :: u_vert, v_vert
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: div_c
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: rho

    REAL(wp), INTENT(out) :: ddt_u(nproma,p_patch%nlev,p_patch%nblks_c) !< u tendency
    REAL(wp), INTENT(out) :: ddt_v(nproma,p_patch%nlev,p_patch%nblks_c) !< v tendency

    REAL(wp) :: flux_up_v, flux_dn_v, flux_up_c, flux_dn_c
    REAL(wp) :: vn_vert1, vn_vert2, vn_vert3, vn_vert4, dvt
    REAL(wp) :: inv_rhoe(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp), INTENT(IN) :: vn(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) :: tot_tend(nproma,p_patch%nlev,p_patch%nblks_e)

    INTEGER,  DIMENSION(:,:,:), POINTER :: ividx, ivblk, iecidx, iecblk
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, je, jcn, jbn, jvn
    INTEGER :: nlev, jg

    jg = p_patch%id

    p_int        => p_int_state(jg)

    ! number of vertical levels
    nlev     = p_patch%nlev

    ividx  => p_patch%edges%vertex_idx
    ivblk  => p_patch%edges%vertex_blk

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk

    !$ACC DATA CREATE(inv_rhoe, tot_tend)

    !total tendency
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
    tot_tend(:,:,:) = 0._wp
    !$ACC END KERNELS

    ! Here rho is assumed to be synchronized already, see interface_iconam_aes.

    !density at edge
    CALL cells2edges_scalar(rho, p_patch, p_int%c_lin_e, inv_rhoe,                  &
                            opt_rlstart=grf_bdywidth_e+1, opt_rlend=min_rledge_int, &
                            lacc=.TRUE.)

    rl_start   = grf_bdywidth_e+1
    rl_end     = min_rledge_int
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
#endif
          inv_rhoe(je,jk,jb) = 1._wp / inv_rhoe(je,jk,jb)
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO

    ! 1) First get the horizontal tendencies

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,vn_vert1,vn_vert2,vn_vert3,vn_vert4,&
!$OMP            dvt,jcn,jbn,flux_up_c,flux_dn_c,jvn,flux_up_v,flux_dn_v)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR TILE(32, 4)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
#endif
          vn_vert1 = u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,1)%v1 +   &
                     v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,1)%v2

          vn_vert2 = u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,2)%v1 +   &
                     v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,2)%v2

          vn_vert3 = u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,3)%v1 +   &
                     v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,3)%v2

          vn_vert4 = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,4)%v1 +   &
                     v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))     *   &
                     p_patch%edges%primal_normal_vert(je,jb,4)%v2

          dvt      = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))     *   &
                     p_patch%edges%dual_normal_vert(je,jb,4)%v1   +   &
                     v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))     *   &
                     p_patch%edges%dual_normal_vert(je,jb,4)%v2   -   &
                     (u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))    *   &
                     p_patch%edges%dual_normal_vert(je,jb,3)%v1   +   &
                     v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))     *   &
                     p_patch%edges%dual_normal_vert(je,jb,3)%v2)

          ! tendency in normal direction:
          ! flux = visc*(D_11-2/3DIV) = visc*(2*delta_v/(vert_vert_len/2)-2/3*div_of_stress)

          jcn       = iecidx(je,jb,2)
          jbn       = iecblk(je,jb,2)
          flux_up_c = km_c(jcn,jk,jbn) * ( 4._wp * ( vn_vert4 - vn(je,jk,jb) ) *     &
                      p_patch%edges%inv_vert_vert_length(je,jb) - 2._wp * z_1by3 *   &
                      div_c(jcn,jk,jbn) )


          jcn       = iecidx(je,jb,1)
          jbn       = iecblk(je,jb,1)
          flux_dn_c = km_c(jcn,jk,jbn) * ( 4._wp * ( vn(je,jk,jb) - vn_vert3 ) *     &
                      p_patch%edges%inv_vert_vert_length(je,jb) -  2._wp * z_1by3 *  &
                      div_c(jcn,jk,jbn) )

          ! tendency in tangential direction

          ! D_12 between edge center and the vertex: delta_v/(primal_edge_len/2) +
          ! ((vt4+vt2)/2-(vt3+vt2)/2)/(distance_opp_edges)
          ! flux = D_12 * visc

          ! Note that the tangential velocities at vertices are used in D_12 is an
          ! approximation for speed. Better way is to use vt reconstructed from vn at
          ! each edges. Also, visc at somewhere between edge mid point and the vertex
          ! should be used but this is a good approximation

          jvn       = ividx(je,jb,2)
          jbn       = ivblk(je,jb,2)
          flux_up_v = 0.5_wp * ( km_iv(jvn,jk,jbn) + km_iv(jvn,jk+1,jbn) ) *                      &
                      ( p_patch%edges%tangent_orientation(je,jb) *                                &
                        ( vn_vert2 - vn(je,jk,jb) )    *                                          &
                        p_patch%edges%inv_primal_edge_length(je,jb) * 2._wp +                     &
                        dvt * p_patch%edges%inv_vert_vert_length(je,jb) )

          jvn       = ividx(je,jb,1)
          jbn       = ivblk(je,jb,1)
          flux_dn_v = 0.5_wp * ( km_iv(jvn,jk,jbn) + km_iv(jvn,jk+1,jbn) ) *                      &
                      ( p_patch%edges%tangent_orientation(je,jb) *                                &
                        ( vn(je,jk,jb) - vn_vert1 )    *                                          &
                        p_patch%edges%inv_primal_edge_length(je,jb) * 2._wp +                     &
                        dvt * p_patch%edges%inv_vert_vert_length(je,jb) )

          tot_tend(je,jk,jb) = ( ( flux_up_c - flux_dn_c ) *                                      &
                                 p_patch%edges%inv_dual_edge_length(je,jb) +                      &
                                 p_patch%edges%tangent_orientation(je,jb) *                       &
                                 ( flux_up_v - flux_dn_v ) *                                      &
                                 p_patch%edges%inv_primal_edge_length(je,jb) * 2._wp ) *          &
                               inv_rhoe(je,jk,jb)

        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO
!$OMP END PARALLEL
    !$ACC WAIT

    CALL sync_patch_array(SYNC_E, p_patch, tot_tend)
    CALL rbf_vec_interpol_cell(tot_tend, p_patch, p_int, ddt_u, ddt_v, opt_rlend=min_rlcell_int)

  NULLIFY(p_int)

  !$ACC WAIT
  !$ACC END DATA

  END SUBROUTINE diffuse_hori_velocity
  !-------------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------------
  !>
  !! diffuse_vert_velocity
  !!------------------------------------------------------------------------
  !! Calculate the SGS diffusion term for vertical velocity component
  !! - Uses the forward Euler time scheme in time split (sequential) manner adopted
  !!   in the NH version.
  !! - Option to switch on implicit scheme in vertical
  !! - only solves for jk=2 to nlev. The bottom and top boundaries are left untouched
  !!------------------------------------------------------------------------
  SUBROUTINE diffuse_vert_velocity( nproma,                   &
                                  & p_patch,                  &
                                  & inv_rho_ic, w_vert, w_ie, &
                                  & km_c, km_iv, km_ic,       &
                                  & u_vert, v_vert, div_c,    &
                                  & pum1, pvm1, pwm1, vn,     &
                                  & ddt_w, dt)

    INTEGER,INTENT(in) :: nproma
    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch      !< single patch
    TYPE(t_nh_metrics) ,POINTER :: p_nh_metrics
    TYPE(t_int_state)  ,POINTER :: p_int         !< interpolation state
    REAL(wp),          INTENT(in)        :: km_ic(nproma,p_patch%nlev+1,p_patch%nblks_c)
    REAL(wp),          INTENT(in)        :: dt

    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: inv_rho_ic
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev+1,p_patch%nblks_v) :: w_vert
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev+1,p_patch%nblks_e) :: w_ie
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: km_c
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev+1,p_patch%nblks_v) :: km_iv
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_v) :: u_vert, v_vert
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: div_c
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: pum1, pvm1
    REAL(wp), INTENT(IN), DIMENSION(nproma,p_patch%nlevp1,p_patch%nblks_c) :: pwm1
    REAL(wp), INTENT(OUT),DIMENSION(nproma,p_patch%nlevp1,p_patch%nblks_c) :: ddt_w

    REAL(wp) :: flux_up_c, flux_dn_c, dvn1, dvn2, dvt1, dvt2, flux_up_v, flux_dn_v
    REAL(wp) :: vt_e(nproma,p_patch%nlev,p_patch%nblks_e), inv_dt
    REAL(wp), INTENT(IN) :: vn(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp), DIMENSION(nproma,p_patch%nlev) :: a, b, c, rhs, var_new

    !interface level variables but only nlev quantities are needed
    REAL(wp) :: hor_tend(nproma,p_patch%nlev,p_patch%nblks_e)

    INTEGER,  DIMENSION(:,:,:), POINTER :: ividx, ivblk, iecidx, iecblk, ieidx, ieblk
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end, jg
    INTEGER :: jk, jb, je, jc, jcn, jbn, jvn
    INTEGER  :: nlev

    !patch id
    jg = p_patch%id

    p_nh_metrics => p_nh_state(jg)%metrics
    p_int        => p_int_state(jg)

    ! number of vertical levels
    nlev     = p_patch%nlev

    inv_dt  = 1._wp / dt

    ividx  => p_patch%edges%vertex_idx
    ivblk  => p_patch%edges%vertex_blk

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk

    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    !$ACC DATA &
    !---- Argument arrays - intent(out)
    !$ACC   CREATE(vt_e, hor_tend) &
    !$ACC   CREATE(a, b, c, rhs, var_new) &
    !$ACC   PRESENT(km_c, km_ic, km_iv, u_vert, v_vert, w_vert, w_ie) &
    !$ACC   PRESENT(p_int, div_c, pum1, pvm1, pwm1, ddt_w, inv_rho_ic) &
    !$ACC   PRESENT(p_nh_metrics, p_patch, ividx, ivblk, iecidx, iecblk, ieblk, ieidx)

    !Some initializations
    !total tendency
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
    a = 0._wp; c = 0._wp
    ddt_w(:,:,:) = 0._wp

    IF(p_test_run)THEN
      hor_tend(:,:,:) = 0._wp
    END IF
    !$ACC END KERNELS

    CALL rbf_vec_interpol_edge( vn, p_patch, p_int, vt_e, opt_rlend=min_rledge_int-1, &
                                opt_acc_async=.TRUE.)

    ! 1) Get horizontal tendencies at half level edges
    rl_start   = grf_bdywidth_e
    rl_end     = min_rledge_int-1
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx,jcn,jbn,dvn1,dvn2,flux_up_c,flux_dn_c,&
!$OMP            jvn,dvt1,dvt2,flux_up_v,flux_dn_v)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR TILE(32, 4)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 2, nlev
#else
      DO jk = 2, nlev
        DO je = i_startidx, i_endidx
#endif

          ! tendency in normal direction
          ! flux = visc_c * D_31_c where D_31(=D_13) is calculated at half level
          ! cell center

          jcn   = iecidx(je,jb,2)
          jbn   = iecblk(je,jb,2)

          dvn2  = pum1(jcn,jk-1,jbn) * p_patch%edges%primal_normal_cell(je,jb,2)%v1 +  &
                  pvm1(jcn,jk-1,jbn) * p_patch%edges%primal_normal_cell(je,jb,2)%v2 -  &
                  pum1(jcn,jk,jbn) * p_patch%edges%primal_normal_cell(je,jb,2)%v1   -  &
                  pvm1(jcn,jk,jbn) * p_patch%edges%primal_normal_cell(je,jb,2)%v2

          flux_up_c = km_ic(jcn,jk,jbn) * (                                            &
                      dvn2 * p_nh_metrics%inv_ddqz_z_half(jcn,jk,jbn) +                &
                      ( w_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) - w_ie(je,jk,jb) ) *  &
                      p_patch%edges%inv_vert_vert_length(je,jb) * 2.0_wp )

          jcn   = iecidx(je,jb,1)
          jbn   = iecblk(je,jb,1)

          dvn1  = pum1(jcn,jk-1,jbn) * p_patch%edges%primal_normal_cell(je,jb,1)%v1 +  &
                  pvm1(jcn,jk-1,jbn) * p_patch%edges%primal_normal_cell(je,jb,1)%v2 -  &
                  pum1(jcn,jk,jbn) * p_patch%edges%primal_normal_cell(je,jb,1)%v1   -  &
                  pvm1(jcn,jk,jbn) * p_patch%edges%primal_normal_cell(je,jb,1)%v2


          flux_dn_c = km_ic(jcn,jk,jbn) * (                                            &
                      dvn1 * p_nh_metrics%inv_ddqz_z_half(jcn,jk,jbn) +                &
                      ( w_ie(je,jk,jb) - w_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) ) *  &
                      p_patch%edges%inv_vert_vert_length(je,jb) * 2.0_wp )


         ! tendency in tangential direction
         ! flux = visc_v * D_32_v where D_32(= D_23) is calculated at half level
         ! between vertex and edge center

          jvn  = ividx(je,jb,2)
          jbn  = ivblk(je,jb,2)

          dvt2 = ( u_vert(jvn,jk-1,jbn) * p_patch%edges%dual_normal_vert(je,jb,2)%v1 +            &
                   v_vert(jvn,jk-1,jbn) * p_patch%edges%dual_normal_vert(je,jb,2)%v2 +            &
                   vt_e(je,jk-1,jb) ) * 0.5_wp  -                                                 &
                 ( u_vert(jvn,jk,jbn) * p_patch%edges%dual_normal_vert(je,jb,2)%v1 +              &
                   v_vert(jvn,jk,jbn) * p_patch%edges%dual_normal_vert(je,jb,2)%v2 +              &
                   vt_e(je,jk,jb) ) * 0.5_wp

          flux_up_v = km_iv(jvn,jk,jbn) * ( dvt2 * p_nh_metrics%inv_ddqz_z_half_v(jvn,jk,jbn) +   &
                      p_patch%edges%tangent_orientation(je,jb) * ( w_vert(jvn,jk,jbn) -           &
                      w_ie(je,jk,jb) ) / p_patch%edges%edge_cell_length(je,jb,2) )


          jvn  = ividx(je,jb,1)
          jbn  = ivblk(je,jb,1)

          dvt1 = ( u_vert(jvn,jk-1,jbn) * p_patch%edges%dual_normal_vert(je,jb,1)%v1 +            &
                   v_vert(jvn,jk-1,jbn) * p_patch%edges%dual_normal_vert(je,jb,1)%v2 +            &
                   vt_e(je,jk-1,jb) ) * 0.5_wp           - &
                 ( u_vert(jvn,jk,jbn) * p_patch%edges%dual_normal_vert(je,jb,1)%v1 +              &
                   v_vert(jvn,jk,jbn) * p_patch%edges%dual_normal_vert(je,jb,1)%v2 +              &
                   vt_e(je,jk,jb) ) * 0.5_wp


          flux_dn_v = km_iv(jvn,jk,jbn) * ( dvt1 * p_nh_metrics%inv_ddqz_z_half_v(jvn,jk,jbn) +   &
                      p_patch%edges%tangent_orientation(je,jb) * ( w_ie(je,jk,jb) -               &
                      w_vert(jvn,jk,jbn) ) / p_patch%edges%edge_cell_length(je,jb,1) )

          hor_tend(je,jk,jb) = ( flux_up_c - flux_dn_c ) *                                        &
                                 p_patch%edges%inv_dual_edge_length(je,jb) +                      &
                                 p_patch%edges%tangent_orientation(je,jb)  *                      &
                               ( flux_up_v - flux_dn_v ) *                                        &
                                 p_patch%edges%inv_primal_edge_length(je,jb) * 2._wp

        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO
!$OMP END PARALLEL

    ! Interpolate horizontal tendencies to w point: except top and bottom boundaries
    ! w==0 at these boundaries
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 2, nlev
#else
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
#endif
          ddt_w(jc,jk,jb)    = inv_rho_ic(jc,jk,jb)       *                                       &
                               ( hor_tend(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) *                     &
                                 p_int%e_bln_c_s(jc,1,jb) +                                       &
                                 hor_tend(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) *                     &
                                 p_int%e_bln_c_s(jc,2,jb) +                                       &
                                 hor_tend(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) *                     &
                                 p_int%e_bln_c_s(jc,3,jb) )
         END DO
       END DO
       !$ACC END PARALLEL
    END DO
!$OMP END DO
!$OMP END PARALLEL

    ! 2) Vertical tendency: evaluated at w point

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx,a,b,c,rhs,var_new)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,      &
                          i_startidx, i_endidx, rl_start, rl_end)
     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
     !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
     DO jc = i_startidx, i_endidx
       DO jk = 3, nlev-1
#else
     DO jk = 3, nlev-1
       DO jc = i_startidx, i_endidx
#endif

           a(jc,jk)   = - 2._wp * km_c(jc,jk-1,jb) * p_nh_metrics%inv_ddqz_z_full(jc,jk-1,jb) * &
                          p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) * inv_rho_ic(jc,jk,jb)

           c(jc,jk)   = - 2._wp * km_c(jc,jk,jb) * p_nh_metrics%inv_ddqz_z_full(jc,jk,jb)     * &
                          p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) * inv_rho_ic(jc,jk,jb)

           b(jc,jk)   =  inv_dt - a(jc,jk) - c(jc,jk)

           rhs(jc,jk) =  pwm1(jc,jk,jb) * inv_dt +                                       &
                         2._wp * ( km_c(jc,jk,jb)   * z_1by3 * div_c(jc,jk,jb)     -            &
                                   km_c(jc,jk-1,jb) * z_1by3 * div_c(jc,jk-1,jb) ) *            &
                         p_nh_metrics%inv_ddqz_z_half(jc,jk,jb) * inv_rho_ic(jc,jk,jb)
       END DO
     END DO
     !$ACC END PARALLEL

        ! Boundary treatment
        !--------------------------------------------------------
        ! jk = 2 (w == 0)
        !--------------------------------------------------------
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          c(jc,2)   = - 2._wp * km_c(jc,2,jb) * p_nh_metrics%inv_ddqz_z_full(jc,2,jb) *         &
                        p_nh_metrics%inv_ddqz_z_half(jc,2,jb) * inv_rho_ic(jc,2,jb)

          b(jc,2)   = inv_dt - c(jc,2) + 2._wp * km_c(jc,1,jb) *                                &
                      p_nh_metrics%inv_ddqz_z_full(jc,1,jb) *                                   &
                      p_nh_metrics%inv_ddqz_z_half(jc,2,jb) * inv_rho_ic(jc,2,jb)

          rhs(jc,2) = pwm1(jc,2,jb) * inv_dt +                                           &
                      2._wp * ( km_c(jc,2,jb) * z_1by3 * div_c(jc,2,jb) -                       &
                                km_c(jc,1,jb) * z_1by3 * div_c(jc,1,jb) ) *                     &
                      p_nh_metrics%inv_ddqz_z_half(jc,2,jb) * inv_rho_ic(jc,2,jb)
        END DO
        !$ACC END PARALLEL
        !--------------------------------------------------------
        ! jk = nlev (w == 0)
        !--------------------------------------------------------
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          a(jc,nlev)   = - km_c(jc,nlev-1,jb) * p_nh_metrics%inv_ddqz_z_full(jc,nlev-1,jb) *    &
                           p_nh_metrics%inv_ddqz_z_half(jc,nlev,jb) * 2._wp *                   &
                           inv_rho_ic(jc,nlev,jb)

          b(jc,nlev)   =   inv_dt - a(jc,nlev) + 2._wp * km_c(jc,nlev,jb) *                     &
                           p_nh_metrics%inv_ddqz_z_full(jc,nlev,jb) *                           &
                           p_nh_metrics%inv_ddqz_z_half(jc,nlev,jb) * inv_rho_ic(jc,nlev,jb)

          rhs(jc,nlev) =   pwm1(jc,nlev,jb) * inv_dt +                                   &
                           2._wp * ( km_c(jc,nlev,jb) * z_1by3 * div_c(jc,nlev,jb) -            &
                                     km_c(jc,nlev-1,jb) * z_1by3 * div_c(jc,nlev-1,jb) ) *      &
                           p_nh_metrics%inv_ddqz_z_half(jc,nlev,jb) * inv_rho_ic(jc,nlev,jb)
        END DO
        !$ACC END PARALLEL

        !$ACC WAIT
        CALL tdma_solver_vec(a,b,c,rhs,2,nlev,i_startidx,i_endidx,var_new)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = 2, nlev
#else
        DO jk = 2, nlev
          DO jc = i_startidx, i_endidx
#endif
            ddt_w(jc,jk,jb) = ddt_w(jc,jk,jb) +                                  &
                                     ( var_new(jc,jk) - pwm1(jc,jk,jb) ) * inv_dt
          END DO
        END DO
        !$ACC END PARALLEL

    END DO !jb
!$OMP END DO
!$OMP END PARALLEL

  NULLIFY(p_nh_metrics)
  NULLIFY(p_int)

  !$ACC WAIT
  !$ACC END DATA

  END SUBROUTINE diffuse_vert_velocity
  !-------------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------------
  !>
  !! diffuse_scalar
  !!------------------------------------------------------------------------
  !! Calculate the SGS diffusion term for cell based scalars
  !! - Uses the forward Euler time scheme in time split (sequential) manner adopted
  !!   in the NH version.
  !! - Option to switch on implicit scheme in vertical
  !!------------------------------------------------------------------------
  SUBROUTINE diffuse_scalar( nproma, var_temp, &
                           & p_patch,          &
                           & km_ie,            &
                           & hori_tend,        &
                           & rho,              &
                           & scalar_name,      &
                           & rturb_prandtl)

    INTEGER,INTENT(in) :: nproma
    REAL(wp),          INTENT(in)           :: var_temp(:,:,:)      ! input scalar
    TYPE(t_patch),     INTENT(in),TARGET    :: p_patch         !< single patch
    TYPE(t_int_state)  ,POINTER :: p_int         !< interpolation state
    REAL(wp),          INTENT(in)           :: rho(:,:,:)      ! density at cell center
    INTEGER,           INTENT(in)           :: scalar_name
    REAL(wp), INTENT(IN) :: rturb_prandtl !< inverse turbulent prandtl number
    REAL(wp), INTENT(IN)                    :: km_ie(nproma,p_patch%nlev+1,p_patch%nblks_e)
    REAL(wp), INTENT(OUT)                   :: hori_tend(nproma,p_patch%nlev,p_patch%nblks_c) !< total tendency
    REAL(wp)                                :: var(nproma,p_patch%nlev,p_patch%nblks_c)      ! input scalar

    !Local variables
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end
    INTEGER :: jk, jb, je, jc, jg
    INTEGER :: nlev

    INTEGER,  DIMENSION(:,:,:), POINTER :: iecidx, iecblk, ieidx, ieblk

    REAL(wp) :: nabla2_e(nproma,p_patch%nlev,p_patch%nblks_e)

    INTEGER, PARAMETER :: tracer_water = 2

    !patch id
    jg = p_patch%id

    p_int        => p_int_state(jg)

    ! number of vertical levels
    nlev = p_patch%nlev

    iecidx => p_patch%edges%cell_idx
    iecblk => p_patch%edges%cell_blk

    ieidx => p_patch%cells%edge_idx
    ieblk => p_patch%cells%edge_blk

    !$ACC DATA &
    !---- Argument arrays - intent(out)
    !$ACC   CREATE(nabla2_e, var) &
    !$ACC   PRESENT(p_patch, km_ie, rho, p_int, hori_tend) &
    !$ACC   PRESENT(iecidx, iecblk, ieidx, ieblk, var_temp)

    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) DEFAULT(PRESENT) ASYNC(1)
    DO jc = 1, p_patch%nblks_c
       DO je = 1, p_patch%nlev
          DO jk = 1, nproma
             hori_tend(jk,je,jc) = 0._wp
             var(jk,je,jc) = var_temp(jk,je,jc)
          END DO
       END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC WAIT
    CALL sync_patch_array(SYNC_C, p_patch, var)

    !1) First set local vars to 1 for other scalars
    !   Soon get different routines for different scalars

!$OMP PARALLEL PRIVATE(rl_start, rl_end, i_startblk, i_endblk)

    !---------------------------------------------------------------
    ! Horizontal diffusion (conservative; following mo_nh_diffusion)
    !---------------------------------------------------------------

    !include halo points and boundary points because these values will be
    !used in next loop
    rl_start   = grf_bdywidth_e
    rl_end     = min_rledge_int-1
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(jk,je,jb,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      ! compute kh_ie * grad_horiz(var)
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
#endif
          nabla2_e(je,jk,jb) = 0.5_wp * ( km_ie(je,jk,jb) + km_ie(je,jk+1,jb) ) *               &
                               rturb_prandtl *                                                  &
                               p_patch%edges%inv_dual_edge_length(je,jb) *                      &
                               ( var(iecidx(je,jb,2),jk,iecblk(je,jb,2)) -                      &
                                 var(iecidx(je,jb,1),jk,iecblk(je,jb,1)) )
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ENDDO
!$OMP END DO

    ! now compute the divergence of the quantity above
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
#endif
          ! horizontal tendency
          hori_tend(jc,jk,jb) = (                                                  &
                        nabla2_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%geofac_div(jc,1,jb) +  &
                        nabla2_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%geofac_div(jc,2,jb) +  &
                        nabla2_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%geofac_div(jc,3,jb) ) / rho(jc,jk,jb)
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    IF (is_dry_cbl .AND. scalar_name==tracer_water) THEN
!$OMP PARALLEL
      CALL init(hori_tend(:,:,:), lacc=.TRUE.)
!$OMP END PARALLEL
    END IF

  NULLIFY(p_int)

  !$ACC WAIT
  !$ACC END DATA

  END SUBROUTINE diffuse_scalar
  !-------------------------------------------------------------------------------------


  !
  !! stability_function_mom
  !! Taken from COSMO docs and Holstag & Boville 1992
  !!------------------------------------------------------------------------
  FUNCTION stability_function_mom(RIB, hz0, tc) RESULT(stab_fun)
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
  !>
  !! stability_function_heat
  !!------------------------------------------------------------------------
  FUNCTION stability_function_heat(RIB, hzh, tc) RESULT(stab_fun)
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

  !>
  !! factor_mom
  !!------------------------------------------------------------------------
  !! Businger Dyer similarity profile:
  !! Louis (1979) A Parametirc model of vertical eddy fluxes in the atmosphere
  !! and R. B. Stull's book
  !!------------------------------------------------------------------------
  FUNCTION businger_mom(z0, z1, L) RESULT(factor)
     REAL(wp), INTENT(IN) :: z0, z1, L
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
  !>
  !! factor_heat
  !!------------------------------------------------------------------------
  !! Businger Dyer similarity profile:
  !! Louis (1979) A Parametirc model of vertical eddy fluxes in the atmosphere
  !! and R. B. Stull's book
  !!------------------------------------------------------------------------
  FUNCTION businger_heat(z0, z1, L) RESULT(factor)
     REAL(wp), INTENT(IN) :: z0, z1, L
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


  !-------------
END MODULE mo_turb_vdiff_sma
