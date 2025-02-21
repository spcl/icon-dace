! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! This file has been modified for the use in ICON
!-------------------------------------------------------------------------------
! SPDX-License-Identifier: Apache-2.0
!-------------------------------------------------------------------------------
!
! Source module for the sub-grid scale orography
!
!-----------------------------------------------------------------------------
!
! Description:
!   The module "mo_sso_cosmo" performs calculations related to the
!   parameterization of sub-grid scale orographic (SSO) effects. The present
!   SSO scheme is based on Lott and Miller (1997).
!
! For ICON:
!   The subroutine "sso" is directly called in the interface routine
!   "mo_nh_interface", the variables are transferred via a variable list.
!
! Reference  The parameterization package has been extracted from the DWD
!   limited area model COSMO (V4_16 of the DWD lm_f90 library) and is
!   originally based on an earlier version by M. Miller and F. Lott (ECMWF).
!
! Lott, F. and M. J. Miller, 1997: A new subgrid-scale
!   orographic drag parametrization: Its formulation and testing.
!   Q. J. R. Meteor. Soc., 123, 101-127.
!
!-----------------------------------------------------------------------------

MODULE mo_sso_cosmo

! Declarations:
!
! Modules used:


USE mo_kind   , ONLY :   &
    wp, vp,      & ! KIND-type parameter for real variables
    i4             ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------

USE mo_physical_constants , ONLY :   &
    r_d   => rd   , & ! gas constant for dry air
    cp_d  => cpd  , & ! specific heat of dry air at constant pressure
    g     => grav     ! acceleration due to gravity

USE mo_nwp_parameters,  ONLY: t_phy_params
USE mo_exception,       ONLY: message
USE mo_fortran_tools,   ONLY: set_acc_host_or_device
! end of mo_physical_constants


!==============================================================================

IMPLICIT NONE

PRIVATE


!------------------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------------------

PUBLIC :: sso

!==============================================================================

! Declarations

! The following variables are tunable parameters for the sub-grid scale
! orography scheme. They are global in this module.

REAL (KIND = wp) ::      &
!
! Tunable parameters
! ------------------
  Gkdrag, gkdrag_add      , &   ! gw drag constants (set in mo_nwp_tuning_nml)
  Gkwake                  , &   ! gw drag constant (set in mo_nwp_tuning_nml)
  Grcrit, grcrit_add      , &   ! critical Richardson numbers (set in mo_nwp_tuning_nml)
  Gfrcrit                 , &   ! critical Froude number (determines depth of blocking layer; set in mo_nwp_tuning_nml)
  minsso, minsso_gwd      , &   ! minimum SSO standard deviation (m) for which SSO information is used (blocking/GWD part)
  blockred                , &   ! multiple of SSO standard deviation above which blocking tendency is reduced

! Security constants
! ------------------
  Gvsec   = 0.10_wp   , &   ! to secure the projection calculation
  Gssec   = 1.E-12_wp , &   ! to secure stability
  Gtsec   = 1.E-07_wp       ! to secure the stress calculation

!==============================================================================
! Module procedures in "mo_sso_cosmo"
!==============================================================================

CONTAINS


!==============================================================================
!+ Module procedure in "mo_sso_cosmo"
!------------------------------------------------------------------------------

SUBROUTINE sso (                                                       &
           ie     , ke     , ke1    ,  istart  , iend   ,              &
           ppf    , pph    , pfif   , pt       , pu , pv  , pfis     , &
           psso_stdh, psso_gamma, psso_theta, psso_sigma, sfcfric_fac, &
           pdt    , mkenvh, params, sso_lat_mask, ldebug ,             &
           pdu_sso, pdv_sso, pustr_sso, pvstr_sso, pvdis_sso, lacc )

!------------------------------------------------------------------------------
!
! Purpose:   The module procedure sso performs the parameterisation of
!            sub-grid scale orographic effects according to Lott and
!            Miller (1997). It computes the physical tendencies of the
!            prognostic variables u, v and T due to vertical transports
!            by sub-grid scale orographically excited gravity waves.
!
! Method:    The scheme consists of three parts, i.e.
!             - the calculation of lowlevel drag due to wake effects of
!               the sub-grid scale orography
!             - the gravity wave stress and
!             - the stress profile in the vertical
!            The stress is computed based on the low level wind, the
!            static stability and subgrid orographic parameters (i.e.
!            standard deviation of height, angle of the main orographic
!            axis relative to east, anisotropy and slope of the orography).
!            At each level a wave Richardson number is computed and by
!            requiring that its value is never less than a critical one,
!            a value of stress can be determined for each model level.
!            The critical Froude number of the flow determines the depth
!            of the low level drag.
!
! Reference: Paper on SSO scheme by Lott and Miller (1997).
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------

!     Input variables
!     ---------------

!     Scalars
!     -------
      INTEGER, INTENT(IN) ::  &
      ie        ,    & ! number of grid points in first (zonal) direction
      ke        ,    & ! number of grid points in vertical direction
      ke1              ! ke + 1

      INTEGER, INTENT(IN) ::  &
      istart    ,    & ! start index for first (zonal) direction
      iend             ! end index for first (zonal) direction


      LOGICAL, OPTIONAL, INTENT(in)    :: lacc ! If true, use openacc

      ! Tuning parameters
      TYPE(t_phy_params),INTENT(in)    :: params

!     Grid scale variables
!     --------------------
      REAL(KIND=wp), INTENT(IN) :: pt   (:,:)  ! (ie,ke)
!                            temperature at full levels         (K)
      REAL(KIND=wp), INTENT(IN) :: pu   (:,:)  ! (ie,ke)
!                            zonal wind component               (m/s)
      REAL(KIND=wp), INTENT(IN) :: pv   (:,:)  ! (ie,ke)
!                            meridional wind component          (m/s)
      REAL(KIND=wp), INTENT(IN) :: pfif (:,:)  ! (ie,ke)
!                            geopotential at full levels        (m**2/s**2)
      REAL(KIND=wp), INTENT(IN) :: pfis (:)  ! (ie)
!                            geopotential at surface            (m**2/s**2)
      REAL(KIND=wp), INTENT(IN) :: pph  (:,:)  ! (ie,ke1)
!                            pressure at half levels            (Pa)
      REAL(KIND=wp), INTENT(IN) :: ppf  (:,:)  ! (ie,ke)
!                            pressure at full levels            (Pa)

!     GWD-Parameters for each model grid point
!     (standard deviation, anisotropy, angle and slope)
!     -------------------------------------------------
      REAL(KIND=wp), INTENT(IN) :: psso_stdh  (:)  !  (ie)
      REAL(KIND=wp), INTENT(inout) :: psso_gamma (:)  !  (ie)
      REAL(KIND=wp), INTENT(IN) :: psso_theta (:)  !  (ie)
      REAL(KIND=wp), INTENT(IN) :: psso_sigma (:)  !  (ie)

      ! Latitude mask for improved GWD tuning
      REAL(KIND=wp), INTENT(IN) :: sso_lat_mask (:)  !  (ie)

      ! Factor for adaptive tuning of surface friction
      REAL(KIND=wp), INTENT(IN) :: sfcfric_fac(:)  !  (ie)

      REAL(KIND=wp) :: pdt           ! time step

      LOGICAL ldebug ! debug control switch

!     Output variables
!     ----------------

!     Tendencies of T, u and v
!
      REAL(KIND=vp) :: pdv_sso(:,:)  ! (ie,ke)
      REAL(KIND=vp) :: pdu_sso(:,:)  ! (ie,ke)

!     Surface (u,v) momentum fluxes and vertically integrated dissipation
!
      REAL(KIND=wp),INTENT(OUT), OPTIONAL :: pustr_sso(:)  !  (ie)
      REAL(KIND=wp),INTENT(OUT), OPTIONAL :: pvstr_sso(:)  !  (ie)
      REAL(KIND=wp),INTENT(OUT), OPTIONAL :: pvdis_sso(:)  !  (ie)

!     Arrays and variables local to *sso* or used for communication
!     with higher level subroutines
!     -----------------------------

!     Indicators of significant levels for the operation of the scheme
      INTEGER, INTENT(OUT) ::  mkenvh (ie) ! level index denoting top of envelope layer
      INTEGER  mcrit  (ie)
      INTEGER  mkcrith(ie)
      INTEGER  mknu   (ie)
      INTEGER  mknu2  (ie)

      LOGICAL  lo_sso (ie), lo_sso_gwd(ie)

      LOGICAL :: lzacc ! non-optional version of lacc

      REAL(KIND=wp) :: zfi    (ie,ke)
      ! geopotential minus surface geopotential (m**2/s**2)
      REAL(KIND=wp) :: ztau   (ie,ke1)
      ! gravity wave stress               (Pa)
      REAL(KIND=wp) :: zstrdu (ie,ke1)
      ! flux of u-momentum (GWD and Blocking)
      REAL(KIND=wp) :: zstrdv (ie,ke1)
      ! flux of v-momentum (GWD and Blocking)
      REAL(KIND=wp) :: zstab  (ie,ke1)
      ! squared Brunt-vaisala frequency   (1/s**2)
      REAL(KIND=wp) :: zvph   (ie,ke1)
      ! wind profile projected onto plane of low level wind
      REAL(KIND=wp) :: zrho   (ie,ke1)
      ! density at half levels            (kg/m**3)
      REAL(KIND=wp) :: zri    (ie,ke1)
      ! mean flow Richardson number       (-)
      REAL(KIND=wp) :: zpsi   (ie,ke1)
      ! angle between flow and main axis of orography
      !
      REAL(KIND=wp) :: zzdep  (ie,ke)
      !
      REAL(KIND=wp) :: zdudt  (ie)
      ! sso-tendency for zonal wind    (m/s**2)
      REAL(KIND=wp) :: zdvdt  (ie)
      ! sso-tendency for merid.wind    (m/s**2)
      REAL(KIND=wp) :: zdtdt  (ie)
      ! sso-tendency for temperature   (K/s)
      REAL(KIND=wp) :: zulow  (ie)
      ! u-component of low level wind  (m/s)
      REAL(KIND=wp) :: zvlow  (ie)
      ! v-component of low level wind  (m/s)
      ! directional parameters (see *sso_setup*)
      REAL(KIND=wp) :: zvidis (ie)
      REAL(KIND=wp) :: zd1    (ie)
      REAL(KIND=wp) :: zd2    (ie)
      REAL(KIND=wp) :: zdmod  (ie)
!
!     Utility variables
!     -----------------
      REAL(KIND=wp) :: zgdph,zcons1
      REAL(KIND=wp) :: zdelp,ztemp,zb,zc,zcs,zss,zconb,zabsv,zzd1
      REAL(KIND=wp) :: zratio,zbet,zdt2,fac_sfc
!     REAL(KIND=wp) :: zust,zvst

      INTEGER j1,j3      ! loop indices

!     Set tuning parameters
      Gkdrag  = params%Gkdrag
      Gkdrag_add  = MAX(0._wp,params%Gkdrag_enh-params%Gkdrag)
      Gkwake  = params%Gkwake
      Grcrit  = params%Grcrit
      Grcrit_add  = MAX(0._wp,params%Grcrit_enh-params%Grcrit)
      Gfrcrit = params%Gfrcrit
      minsso  = params%minsso
      minsso_gwd  = params%minsso_gwd
      blockred= params%blockred

      CALL set_acc_host_or_device(lzacc, lacc)

      !Declaration of GPU arrays  
      !$ACC DATA PRESENT(pt, pu, pv, pfif, pfis, pph, ppf, psso_stdh, psso_gamma, sfcfric_fac) &
      !$ACC   PRESENT(psso_theta, psso_sigma, pdv_sso, pdu_sso, pustr_sso, pvstr_sso, pvdis_sso, sso_lat_mask) &
      !$ACC   CREATE(mcrit, mkcrith, mknu, mknu2, lo_sso, lo_sso_gwd) &
      !$ACC   CREATE(zfi, ztau, zstrdu, zstrdv, zstab, zvph, zrho, zri, zpsi, zzdep) &
      !$ACC   CREATE(zdudt, zdvdt, zdtdt, zulow, zvlow, zvidis, zd1, zd2, zdmod, mkenvh) IF(lzacc)

!     Timestep is already set for 2TL or 3TL scheme, respectively,
!     in calling routine organize_sso.
!     zdt2  = 2._wp* pdt
      zdt2  = pdt

!     zcons1=1._wp/(G*pdt*2._wp)
      zcons1=1._wp/(G*zdt2)

!     Initialize tendencies and compute geopotential above ground
!     ===========================================================
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO j3=1,ke
        DO j1=istart,iend
          pdu_sso(j1,j3) = 0.0_vp
          pdv_sso(j1,j3) = 0.0_vp
          zfi    (j1,j3) = pfif(j1,j3)-pfis(j1)
        END DO
      END DO

!     Control operation of scheme by selection of points with standard
!     deviation of sub-grid scale orography > 10 m only (default value of minsso)
!     =================================================
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO j1=istart,iend
        lo_sso(j1)     = MERGE(.TRUE., .FALSE., psso_stdh(j1) > minsso)
        lo_sso_gwd(j1) = MERGE(.TRUE., .FALSE., psso_stdh(j1) > MAX(minsso,minsso_gwd))
      END DO
      !$ACC END PARALLEL
      !$ACC WAIT IF(lzacc)

! ========================================================
!     Computation of basic state variables in *sso_setup*
! ========================================================

      CALL sso_setup (                                            &
         ie     , ke   , ke1 ,   istart , iend   ,                &
         pph   , ppf   , pu    , pv  , pt  ,zfi  ,                &
         psso_stdh, sso_lat_mask, psso_theta, psso_gamma, lo_sso, &
         zrho  , zri   , zstab, ztau, zvph, zpsi, zzdep,          &
         zulow , zvlow , zd1  , zd2 ,zdmod, mkcrith ,             &
         mcrit, mkenvh,mknu,mknu2, lacc=lzacc )

! ========================================================
!     Surface gravity wave stress amplitude
! ========================================================

      CALL gw_stress (                                            &
         ie     , ke1 , istart , iend   ,                         &
         zrho,zstab,zvph,psso_stdh,sso_lat_mask,psso_sigma,zdmod, &
         lo_sso_gwd, zfi, mkenvh, ztau, lacc=lzacc )

! ========================================================
!     Gravity wave stress profile
! ========================================================

      CALL gw_profil(                                          &
         ie     , ke     , ke1 , istart , iend ,               &
         pph   , zrho  , zstab , zvph    , zri  ,              &
         ztau  , zdmod , psso_sigma, psso_stdh , sso_lat_mask, &
         mkcrith, mcrit, mkenvh, mknu    , mknu2,              &
         lo_sso_gwd , lacc=lzacc )

! ========================================================
!     Computation of SSO effects' tendencies
! ========================================================

!     Initialisation of tendencies for ALL grid points
!     ------------------------------------------------
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO j1=istart,iend
        zvidis (j1)=0.0_wp
        zdudt  (j1)=0.0_wp
        zdvdt  (j1)=0.0_wp
        zdtdt  (j1)=0.0_wp
      END DO
      !$ACC END PARALLEL

!     Compute and add low level drag tendencies to the GWD ones
!     ---------------------------------------------------------
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO j3=1,ke
        !$ACC LOOP GANG VECTOR
        DO j1=istart,iend

        IF (lo_sso(j1)) THEN

!         Gravity wave drag (cf. documentation EQ.4.13)
!         ---------------------------------------------
          zdelp = pph(j1,j3+1)-pph(j1,j3)
          ztemp = -G*(ztau(j1,j3+1)-ztau(j1,j3))           &
                    /(zvph(j1,ke1)*zdelp)
          zdudt(j1)=(zulow(j1)*zd1(j1)-zvlow(j1)*zd2(j1))  &
                                      *ztemp/zdmod(j1)
          zdvdt(j1)=(zvlow(j1)*zd1(j1)+zulow(j1)*zd2(j1))  &
                                      *ztemp/zdmod(j1)
          IF (j3 < 4) THEN
           zdudt(j1)= SIGN(MIN(ABS(zdudt(j1)),20._wp/3600._wp),zdudt(j1))
           zdvdt(j1)= SIGN(MIN(ABS(zdvdt(j1)),20._wp/3600._wp),zdvdt(j1))
          ENDIF

          IF (j3==ke) THEN
            fac_sfc = sfcfric_fac(j1)
          ELSE
            fac_sfc = 1._wp
          ENDIF

!         Low level drag ('blocking') (cf. documentation EQ.4.14 ff.)
!         -----------------------------------------------------------
          IF (j3.GE.mkenvh(j1)) THEN
           zb  = 1.0_wp-0.18_wp*psso_gamma(j1)-0.04_wp*psso_gamma(j1)**2
           zc  = 0.48_wp*psso_gamma(j1)+0.3_wp*psso_gamma(j1)**2
           zcs = COS(zpsi(j1,j3))**2
           zss = 1.0_wp-zcs
           zzd1  =zb*zcs+zc*zss
           zconb =zdt2*Gkwake*psso_sigma(j1)/(2._wp*psso_stdh(j1))
           zabsv =0.5_wp*SQRT(pu(j1,j3)**2+pv(j1,j3)**2)
           zratio=(zcs+psso_gamma(j1)*zss)/(psso_gamma(j1)*zcs+zss)
           zbet  =MAX(0._wp,2._wp-1._wp/zratio)*zconb*zzdep(j1,j3)*zzd1*zabsv
           zbet = zbet * fac_sfc*MIN(1._wp,blockred*psso_stdh(j1)/(zfi(j1,j3)/G))
!          Partially implicit tendency calculation
!          ---------------------------------------
           zdudt(j1)=-pu(j1,j3)/zdt2*(zbet/(1._wp+zbet))
           zdvdt(j1)=-pv(j1,j3)/zdt2*(zbet/(1._wp+zbet))
          END IF

          pdu_sso(j1,j3)=zdudt(j1)
          pdv_sso(j1,j3)=zdvdt(j1)
        ENDIF

        END DO
      END DO     ! loop over vertical layers
      !$ACC END PARALLEL

! ======================================================================
!     Flux computations of original code of *GWDRAG* are not required in
!     COSMO model, but may be reactivated by uncommenting the following
!     statements and declaration of the appropriate arrays
! ======================================================================

!     Stress components and dissipation
!     ---------------------------------

      IF(PRESENT(pvdis_sso))THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO j1=istart,iend
          IF(lo_sso(j1)) THEN
!           pvdis_sso(j1)=zcons1*zvidis(j1)
            pvdis_sso(j1)=zvidis(j1)/G
          ELSE
            pvdis_sso(j1)=0._wp
          ENDIF
        END DO
        !$ACC END PARALLEL
      ENDIF

!     Initialize flux at top
!     ----------------------
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO j1=istart,iend
        zstrdu(j1,1)=0._wp
        zstrdv(j1,1)=0._wp
      END DO
      !$ACC END PARALLEL

!     Increment flux based on tendency in each layer
!     ----------------------------------------------
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO j3=1,ke
        !$ACC LOOP GANG VECTOR
        DO j1=istart,iend
          zgdph=-G  /(pph(j1,j3+1)-pph(j1,j3))
          zstrdu(j1,j3+1)=pdu_sso(j1,j3)/zgdph + zstrdu(j1,j3)
          zstrdv(j1,j3+1)=pdv_sso(j1,j3)/zgdph + zstrdv(j1,j3)
        END DO
      END DO
      !$ACC END PARALLEL

      IF(PRESENT(pustr_sso))THEN
!     Store flux at surface
!     ---------------------
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO j1=istart,iend
          pustr_sso(j1)=zstrdu(j1,ke1)
          pvstr_sso(j1)=zstrdv(j1,ke1)
        END DO
        !$ACC END PARALLEL
      ENDIF
!
!     Control printout
!     ----------------

      IF (ldebug) THEN
      !$ACC UPDATE HOST(pfis, psso_stdh, psso_gamma, psso_theta, psso_sigma) ASYNC(1) IF(lzacc)
      !$ACC UPDATE HOST(pph, ppf, pu, pv, pt, pfif, pdu_sso, pdv_sso) ASYNC(1) IF(lzacc)
      !$ACC WAIT(1) IF(lzacc)
        DO j1=istart,iend
          IF (j1.EQ.55) THEN
            PRINT *, ' '
            PRINT *, ' '
            PRINT *, ' Diagnosis SSO scheme j1=55'
            PRINT *, ' '
            PRINT *, ' fis      : ', pfis       (j1)
            PRINT *, ' sso_stdh : ', psso_stdh  (j1)
            PRINT *, ' sso_gamma: ', psso_gamma (j1)
            PRINT *, ' sso_theta: ', psso_theta (j1)
            PRINT *, ' sso_sigma: ', psso_sigma (j1)
            PRINT *, ' '
            PRINT *, '  j3    ph (Pa)     pf (Pa)     u(m/s)    ',  &
           'v (m/s)     T (K)      fif (m^2/s^2)     du_sso     ',  &
           'dv_sss     dt_sso'
          ENDIF
        ENDDO

        DO j1=istart,iend
          IF (j1.EQ.55) THEN
            DO j3=1,ke
              WRITE (*,'(i3, 8E13.6)') j3, pph(j1,j3+1), ppf(j1,j3),  &
              pu(j1,j3), pv(j1,j3), pt(j1,j3), pfif(j1,j3),     &
              pdu_sso(j1,j3), pdv_sso(j1,j3)
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      !$ACC WAIT IF(lzacc)
      !$ACC END DATA
!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE sso

!==============================================================================
!+ Module procedure in "mo_sso_cosmo"
!------------------------------------------------------------------------------

SUBROUTINE sso_setup (                                      &
           ie     , ke   , ke1 ,   istart , iend   ,        &
           pph ,ppf ,pu ,pv ,pt ,pfi,                       &
           psso_stdh, sso_lat_mask, psso_theta, psso_gamma , lo_sso,      &
           prho  , pri   , pstab, ptau, pvph , ppsi, pzdep, &
           pulow , pvlow , pd1  , pd2 , pdmod,              &
           kkcrith, kcrit, kkenvh,kknu,kknu2, lacc)


!------------------------------------------------------------------------------
!
! Purpose: The module procedure sso_setup sets up parameters for the
!          parameterization of sub-grid scale orographic effects:
!
!          - definition of various reference model levels
!          - Brunt-Vaisala frequency on half levels
!          - mean wind components in the layer between one and two
!            standard deviations of sso above ground
!          - geometry factors, relating the orientation of sso and wind
!          - Phillips parameters
!          - vertical wind profile in plane of gravity wave stress
!          - basic flow Richardson number
!          - calculation of depth of 'blocked' layer
!          - calculation of layer in which low level wave braking occurs
!          - calculation of assumed vertical profile of sso
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------

!     Input arrays and variables
!     ==========================
!
      INTEGER, INTENT(IN) ::  &
      ie        ,    & ! number of grid points in first (zonal) direction
      ke        ,    & ! number of grid points in vertical direction
      ke1              ! ke + 1

      INTEGER, INTENT(IN) ::  &
      istart    ,    & ! start index for first (zonal) direction
      iend             ! end index for first (zonal) direction

      LOGICAL, OPTIONAL, INTENT(IN)  :: lacc ! If true, use openacc

      REAL(KIND=wp) :: pph (:,:) ! (ie,ke1)
      REAL(KIND=wp) :: ppf (:,:) ! (ie,ke)
!
      REAL(KIND=wp) :: pu  (:,:) ! (ie,ke)
      REAL(KIND=wp) :: pv  (:,:) ! (ie,ke)
      REAL(KIND=wp) :: pt  (:,:) ! (ie,ke)
      REAL(KIND=wp) :: pfi (:,:) ! (ie,ke)

!     subgrid scale orography parameters
      REAL(KIND=wp) :: psso_stdh (:) ! (ie)
      REAL(KIND=wp) :: sso_lat_mask (:) ! (ie)
      REAL(KIND=wp) :: psso_theta(:) ! (ie)
      REAL(KIND=wp) :: psso_gamma(:) ! (ie)

      LOGICAL lo_sso(ie)

!     Output arrays
!     =============

      REAL(KIND=wp) :: prho (:,:) ! (ie,ke1)
!     density on half levels          (kg/m**3)
      REAL(KIND=wp) :: pri  (:,:) ! (ie,ke1)
!     mean flow Richardson number     (-)
      REAL(KIND=wp) :: pstab(:,:) ! (ie,ke1)
!     squared Brunt-Vaisala frequency (1/s**2)
      REAL(KIND=wp) :: ptau (:,:) ! (ie,ke1)
!     gravity wave stress profile     (Pa)
      REAL(KIND=wp) :: pvph (:,:) ! (ie,ke1)
!     projected flow on half levels   (m/s)
      REAL(KIND=wp) :: ppsi (:,:) ! (ie,ke1)
!     angle between orography and blocked flow (1:ke)
!                            or low level flow (ke1)
      REAL(KIND=wp) :: pzdep(:,:) ! (ie,ke)
!     height dependency factor for 'blocking' tendency
      REAL(KIND=wp) :: pulow(:) ! (ie)
!     low level zonal wind            (m/s)
      REAL(KIND=wp) :: pvlow(:) ! (ie)
!     low level meridional wind       (m/s)

!     directional parameters
      REAL(KIND=wp) :: pd1  (:) ! (ie)
      REAL(KIND=wp) :: pd2  (:) ! (ie)
      REAL(KIND=wp) :: pdmod(:) ! (ie)

      INTEGER kkcrith(:) ! (ie)
!     maximum level for wave breaking

      INTEGER kcrit  (:) ! (ie)
!     critical level

      INTEGER kkenvh (:) ! (ie)
!     index of top of 'envelope' / blocking level

      INTEGER kknu   (:) ! (ie)
!     model level at 4*stdh

      INTEGER kknu2  (:) ! (ie)
!     model level at 3*stdh

!     local arrays and variables
!     ==========================

      REAL(KIND=wp) :: zvpf   (ie,ke)
      ! projected flow on full levels (m/s)
      REAL(KIND=wp) :: zdp    (ie,ke)
      ! pressure difference between layers
      REAL(KIND=wp) :: zsqst  (ie,ke)
      REAL(KIND=wp) :: znorm  (ie)
      REAL(KIND=wp) :: znup   (ie)
      REAL(KIND=wp) :: znum   (ie)
      REAL(KIND=wp) :: znu    (ie)

      REAL(KIND=wp) :: zcons1, zcons2 ! utility constants
      REAL(KIND=wp) :: zu             ! security for low level zonal wind
      REAL(KIND=wp) :: zb             ! Phillips parameter B
      REAL(KIND=wp) :: zc             ! Phillips parameter C
      REAL(KIND=wp) :: zdelp          ! pressure thickness of layers
      REAL(KIND=wp) :: zvt1,zvt2      ! utility variables for flow projection
      REAL(KIND=wp) :: zst            ! utility variable for stability calculation
      REAL(KIND=wp) :: zdwind         ! utility variable for wind shear calculation
      REAL(KIND=wp) :: zwind          ! utility variable for proj. wind calculation
      REAL(KIND=wp) :: zggeenv,zggeo,zgvar ! geopotential utility variables
      REAL(KIND=wp) :: zhcrit

      INTEGER mknub(ie)
      INTEGER mknul(ie)

      INTEGER mi3h              ! vertical loop limit
      INTEGER j1,j3             ! loop variables

      LOGICAL lo1  (ie,ke1)
      LOGICAL llo               ! utility switch
      LOGICAL :: lzacc

!     The following parameter is a tunable constant for the sub-grid scale
!     orography scheme
!     ================

      INTEGER (KIND=i4) :: Nktopg
                                ! number of topmost layer used to define low level
                                ! flow in case of high vertical model resolution

!-------------------------------------------------------------------------------------

! Begin subroutine

      CALL set_acc_host_or_device(lzacc, lacc)

      Nktopg = ke               ! number of topmost layer used to defined low level

!     computational constants
!     =======================

!     mi3h =(ki3e-ki3s+1)/3
      mi3h =ke/3

      zcons1=1._wp/R_d
      zcons2=G**2/Cp_d

!C*         2.1     DEFINE LOW LEVEL WIND, PROJECT WINDS IN PLANE OF
!C*                 LOW LEVEL WIND, DETERMINE SECTOR IN WHICH TO TAKE
!C*                 THE VARIANCE AND SET INDICATOR FOR CRITICAL LEVELS.
!C
!     security on anisotropy factor and presetting of critical levels

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc) &
      !$ACC   CREATE(lo1, mknul, mknub, znu, znum, znup, znorm, zsqst, zdp, zvpf)
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO j1=istart,iend
        psso_gamma(j1) = MAX(psso_gamma(j1),Gtsec)
        kknu      (j1) = ke
        kknu2     (j1) = ke
        mknub     (j1) = ke
        mknul     (j1) = ke
        lo1(j1,ke1)    =.FALSE.    ! Initialize control variable
      END DO

!
!!!!  define top of low level drag calculation layer (*kkcrit*)
!     and other critical levels
!     ============================================================

      !$ACC LOOP SEQ
      DO j3=ke,mi3h,-1     ! vertical loop
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zhcrit)
        DO j1=istart,iend
          zhcrit = 4._wp*psso_stdh(j1)
          lo1(j1,j3)=((pfi(j1,j3)/G).GT.zhcrit)
          IF(lo1(j1,j3).NEQV.lo1(j1,j3+1)) THEN
            kknu(j1)=j3  ! first layer with height > 4*stdh
          ENDIF
        END DO
      END DO                ! end of vertical loop


      !$ACC LOOP SEQ
      DO j3=ke,mi3h,-1    ! vertical loop
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zhcrit)
        DO j1=istart,iend
          zhcrit          =3._wp*psso_stdh(j1)
          lo1(j1,j3)=((pfi(j1,j3)/G).GT.zhcrit          )
          IF(lo1(j1,j3).NEQV.lo1(j1,j3+1)) THEN
            kknu2(j1)=j3 ! first layer with height > 3*stdh
          ENDIF
        END DO
      END DO                ! end of vertical loop


      !$ACC LOOP SEQ
      DO j3=ke,mi3h,-1    ! vertical loop
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zhcrit)
        DO j1=istart,iend
          zhcrit          =2._wp*psso_stdh(j1)
          lo1(j1,j3)=((pfi(j1,j3)/G).GT.zhcrit          )
          IF(lo1(j1,j3).NEQV.lo1(j1,j3+1)) THEN
            mknub(j1)=j3  ! first layer with height > 2*stdh
          ENDIF
        END DO
      END DO                ! end of vertical loop


      !$ACC LOOP SEQ
      DO j3=ke,mi3h,-1    ! vertical loop
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zhcrit)
        DO j1=istart,iend
          zhcrit          =psso_stdh(j1)
          lo1(j1,j3)=((pfi(j1,j3)/G).GT.zhcrit          )
          IF(lo1(j1,j3).NEQV.lo1(j1,j3+1)) THEN
            mknul(j1)=j3 ! first layer with height > 1*stdh
          ENDIF
        END DO
      END DO                ! end of vertical loop


!     Confine critical level indices to be less or equal to Nktopg
!     ============================================================
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO j1=istart,iend
        kknu(j1) =MIN(kknu(j1),Nktopg)
        mknub(j1)=MIN(mknub(j1),Nktopg)
        IF(mknub(j1).EQ.Nktopg) mknul(j1)=ke
        IF(mknub(j1).EQ.mknul(j1)) mknub(j1) = mknub(j1) - 1
      END DO



      mi3h = MIN(ke-2,MINVAL(kknu(istart:iend)))


!     Initialize various arrays
!     =========================
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO j1=istart,iend
        prho (j1,ke1) = 0.0_wp
        pstab(j1,1  ) = 0.0_wp
        pstab(j1,ke1) = 0.0_wp
        pri  (j1,1  ) = 0.0_wp
        pri  (j1,ke1) = 9999.0_wp
        ppsi (j1,ke1) = 0.0_wp
        pvph (j1,1)   = 0.0_wp
        pulow(j1)     = 0.0_wp
        pvlow(j1)     = 0.0_wp
        kkcrith(j1)   = ke
        kkenvh(j1)    = ke ! default for top of envelope layer
        kcrit(j1)     = 1  ! default for critical level
        znu  (j1)     = 0.0_wp
        znum (j1)     = 0.0_wp
        lo1  (j1,ke1) = .FALSE.
      END DO


!     pressure thickness, density and Brunt-Vaisala frequency (squared)
!     =================================================================

      !$ACC LOOP SEQ
      DO j3=ke,2,-1        ! vertical loop
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO j1=istart,iend
          IF(lo_sso(j1)) THEN
            zdp (j1,j3) = ppf(j1,j3)-ppf(j1,j3-1)
!           density on half levels
            prho (j1,j3) = 2._wp*pph(j1,j3)*zcons1              &
     &                       /(pt(j1,j3)+pt(j1,j3-1))
!           squared Brunt-Vaisala frequency on half levels
            pstab(j1,j3)= 2._wp*zcons2/(pt(j1,j3)+pt(j1,j3-1))  &
     &                      *( 1._wp-Cp_d*prho(j1,j3)           &
     &                           *(pt(j1,j3)-pt(j1,j3-1))       &
     &                               / zdp(j1,j3) )
!           security on Brunt-Vaisala frequency
            pstab(j1,j3)=MAX(pstab(j1,j3),gssec)
            zsqst(j1,j3)=SQRT(pstab(j1,j3))
          ENDIF
        END DO
      END DO                ! end of vertical loop


!     Definition of blocked flow
!     ==========================
      !$ACC LOOP SEQ
      DO j3=ke,mi3h,-1          ! vertical loop
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO j1=istart,iend
          IF(lo_sso(j1)) THEN
            IF(j3.GE.mknub(j1).AND.j3.LE.mknul(j1)) THEN
!             only layers with height between one and two stdh contribute
              pulow(j1) = pulow(j1) +pu(j1,j3)*(pph(j1,j3+1)-pph(j1,j3))
              pvlow(j1) = pvlow(j1) +pv(j1,j3)*(pph(j1,j3+1)-pph(j1,j3))
            END IF
          END IF
        END DO
      END DO                ! end of vertical loop


!     Division by pressure thickness of contributing layers and
!     determination of wind speed of blocked flow

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO j1=istart,iend
        IF(lo_sso(j1)) THEN
          pulow(j1) = pulow(j1) / (pph(j1,mknul(j1)+1)-pph(j1,mknub(j1)))
          pvlow(j1) = pvlow(j1) / (pph(j1,mknul(j1)+1)-pph(j1,mknub(j1)))
          znorm(j1)=MAX(SQRT(pulow(j1)**2+pvlow(j1)**2),Gvsec)
          pvph(j1,ke1)=znorm(j1)  ! projected flow at lower boundary
        END IF
      END DO


!     Axes of subgrid scale orography and plane of profiles
!     =====================================================

      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(llo, ZU, zb, zc)
      DO j1=istart,iend
        IF (lo_sso(j1)) THEN
          llo=(pulow(j1).LT.Gvsec).AND.(pulow(j1).GE.-Gvsec)
          IF(llo) THEN
            ZU=pulow(j1)+2._wp*Gvsec
          ELSE
            ZU=pulow(j1)
          ENDIF
!         angle between principal axis of orog. and low-level wind direction
!         ------------------------------------------------------------------
          ppsi (j1,ke1) = psso_theta(j1)-ATAN(pvlow(j1)/ZU)
!         Phillips parameters B and C
!         ---------------------------
          zb = 1._wp - 0.18_wp*psso_gamma(j1) - 0.04_wp*psso_gamma(j1)**2
          zc = 0.48_wp*psso_gamma(j1) + 0.3_wp*psso_gamma(j1)**2
!         projection parameters D1 and D2 (see documentation)
          pd1  (j1) = zb-(zb-zc)*(SIN(ppsi(j1,ke1))**2)
          pd2  (j1) = (zb-zc)*SIN(ppsi(j1,ke1))*COS(ppsi(j1,ke1))
          pdmod(j1) = SQRT(pd1(j1)**2+pd2(j1)**2)
        END IF
      END DO


!     projection of flow into plane of low level stress  (eq.4.7)
!     ===========================================================
      !$ACC LOOP SEQ
      DO j3=1,ke                  ! vertical loop
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zvt1, zvt2)
        DO j1=istart,iend
          IF(lo_sso(j1)) THEN
            zvt1 = pulow(j1)*pu(j1,j3)+pvlow(j1)*pv(j1,j3)
            zvt2 =-pvlow(j1)*pu(j1,j3)+pulow(j1)*pv(j1,j3)
            zvpf(j1,j3)=(zvt1*pd1(j1)+zvt2*pd2(j1)) / (znorm(j1)*pdmod(j1))
          ENDIF

!      initialize stress array *ptau*, depth of blocked layer *pzdep*
!      and angle array *ppsi* (not for j3=ke1) and reset control
!      variable *llo1*
!      -----------------------------------------------------------------
          ptau(j1,j3)  =0.0_wp
          pzdep(j1,j3) =0.0_wp
          ppsi(j1,j3)  =0.0_wp
          lo1(j1,j3)   =.FALSE.
        END DO
      END DO                ! end of vertical loop

!
!     linear interpolation of projected flow to half levels
!     ========================================================
      !$ACC LOOP SEQ
      DO j3=2,ke     ! vertical loop
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO j1=istart,iend
          IF(lo_sso(j1)) THEN
            pvph(j1,j3)=                                  &
     &     ((pph(j1,j3)-ppf(j1,j3-1))*zvpf(j1,j3) +     &
     &      (ppf(j1,j3)-pph(j1,j3  ))*zvpf(j1,j3-1))    &
     &                         /zdp(j1,j3)
            IF(pvph(j1,j3).LT.Gvsec) THEN  ! critical layer
              pvph(j1,j3)=Gvsec
              if (j3.lt.mknub(j1)) kcrit(j1)=j3
              ENDIF
          ENDIF
        END DO
      END DO                ! end of vertical loop


!     Brunt-Vaisala frequency and density for lowest level
!     ====================================================

      !$ACC LOOP SEQ
      DO j3=mi3h,ke   ! vertical loop
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zst)
        DO j1=istart,iend
          IF(lo_sso(j1)) THEN
            IF(j3.GE.(mknub(j1)+1).AND.j3.LE.mknul(j1)) THEN
              zst=zcons2/pt(j1,j3)*(1._wp-Cp_d*prho(j1,j3)*     &
     &                      (pt(j1,j3)-pt(j1,j3-1))/zdp(j1,j3))
              pstab(j1,ke1)=pstab(j1,ke1)+zst*zdp(j1,j3)
              pstab(j1,ke1)=MAX(pstab(j1,ke1),Gssec)
              prho (j1,ke1)= prho(j1,ke1)                     &
     &                      +pph(j1,j3)*2._wp*zdp(j1,j3)  &
     &                      *zcons1/(pt(j1,j3)+pt(j1,j3-1))
            ENDIF
          ENDIF
        END DO
      END DO                ! end of vertical loop


!     normalization
!     -------------

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO j1=istart,iend
        IF(lo_sso(j1)) THEN
          pstab(j1,ke1)=pstab(j1,ke1) / (ppf(j1,mknul(j1))-ppf(j1,mknub(j1)))
          pstab(j1,ke1)=MAX(pstab(j1,ke1),gssec)
          prho (j1,ke1)=prho(j1,ke1) / (ppf(j1,mknul(j1))-ppf(j1,mknub(j1)))
        END IF
      END DO


!     mean flow Richardson number on half levels
!     ==========================================
      !$ACC LOOP SEQ
      DO j3=2,ke     ! vertical loop
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zdwind)
        DO j1=istart,iend
          IF(lo_sso(j1)) THEN
            zdwind=MAX(ABS(zvpf(j1,j3)-zvpf(j1,j3-1)),Gvsec)
            pri(j1,j3)=pstab(j1,j3)*(zdp(j1,j3) / (G*prho(j1,j3)*zdwind))**2
            pri(j1,j3)=MAX(pri(j1,j3),Grcrit+Grcrit_add*sso_lat_mask(j1))
          ENDIF
        END DO
      END DO                ! end of vertical loop


!     define top of 'envelope' layer (cf. eq.4.8)
!     ===========================================
      !$ACC LOOP SEQ
      DO j3=mi3h,ke-1     ! vertical loop
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zwind, zdelp)
        DO j1=istart,iend
          IF(lo_sso(j1) .AND. j3.GE.kknu2(j1)) THEN
            znum(j1)=znu(j1)
            zwind= (pulow(j1)*pu(j1,j3)                           &
     &             +pvlow(j1)*pv(j1,j3))/                         &
     &         MAX(SQRT(pulow(j1)**2+pvlow(j1)**2),Gvsec)
            zwind=MAX(ABS(zwind),Gvsec)
            zdelp=pph(j1,j3+1)-pph(j1,j3)
!           vertically integrated left side of eq.4.8
            znu(j1) = znu(j1) + (zdelp/G)*                             &
     &              ( (zsqst(j1,j3+1)/prho(j1,j3+1)                    &
     &                +zsqst(j1,j3  )/prho(j1,j3  ) )/2._wp)/zwind
            IF((znum(j1).LE.gfrcrit).AND.(znu(j1).GT.gfrcrit)          &
     &                          .AND.(kkenvh(j1).EQ.ke))               &
     &      kkenvh(j1)=j3
          ENDIF
        END DO
      END DO                ! end of vertical loop


!     dynamical mixing height for the breaking of gravity waves
!     =========================================================
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO j1=istart,iend
        znup(j1)=0.0_wp
        znum(j1)=0.0_wp
      END DO


      !$ACC LOOP SEQ
      DO j3=ke-1,2,-1  ! vertical loop
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zwind, zdelp)
        DO j1=istart,iend
          IF(lo_sso(j1)) THEN
            IF (j3.LT.kkenvh(j1)) THEN    ! only above envelope height
            znum(j1)=znup(j1)
            zwind=(pulow(j1)*pu(j1,j3)+pvlow(j1)*pv(j1,j3))/    &
     &            MAX(SQRT(pulow(j1)**2+pvlow(j1)**2),Gvsec)
            zwind=MAX(ABS(zwind),Gvsec)
            zdelp=pph(j1,j3+1)-pph(j1,j3)
            znup(j1) = znup(j1) + (zdelp/G)*                          &
     &            ( (zsqst(j1,j3+1)/prho(j1,j3+1)                     &
     &              +zsqst(j1,j3  )/prho(j1,j3  ) )/2._wp)/zwind
            IF((znum(j1).LE.1.5_wp).AND.(znup(j1).GT.1.5_wp)           &
     &                          .AND.(kkcrith(j1).EQ.ke))              &
     &      kkcrith(j1)=j3
            ENDIF
          ENDIF
        END DO
      END DO                ! end of vertical loop


!     allow low level wave breaking only above height of 4*stdh
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO j1=istart,iend
        kkcrith(j1)=MIN(kkcrith(j1),kknu(j1))
      END DO


!     directional information for flow blocking
!     =========================================
      !$ACC LOOP SEQ
      DO j3=mi3h,ke     ! vertical loop
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(llo, ZU)
        DO j1=istart,iend
          IF(lo_sso(j1)) THEN
            IF(j3.GE.kkenvh(j1)) THEN  ! only within envelope layer
            llo=(pu(j1,j3).LT.Gvsec).AND.(pu(j1,j3).GE.-Gvsec)
            IF(llo) THEN
              ZU=pu(j1,j3)+2._wp*Gvsec
            ELSE
              ZU=pu(j1,j3)
            ENDIF
          ppsi(j1,j3)=psso_theta(j1)-ATAN(pv(j1,j3)/ZU)
          ENDIF
          ENDIF
        END DO
      END DO                ! end of vertical loop

!     assumed vertical profile of sso for blocking calculations
!     =========================================================
      !$ACC LOOP SEQ
      DO j3=mi3h,ke     ! vertical loop
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zggeenv, zggeo, zgvar)
        DO j1=istart,iend
          IF(lo_sso(j1)) THEN
            IF(j3.GE.kkenvh(j1)) THEN   ! only within envelope layer
              zggeenv= MAX(1._wp, &
     &         (pfi(j1,kkenvh(j1))+pfi(j1,kkenvh(j1)-1))/2._wp)
              zggeo  = MAX(pfi(j1,j3),1._wp)
              zgvar  = MAX(psso_stdh(j1)*G,1._wp)
              pzdep(j1,j3)=SQRT((zggeenv-zggeo)/(zggeo+zgvar))
            END IF
          END IF
        END DO
      END DO                ! end of vertical loop

      !$ACC END PARALLEL

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE sso_setup

!==============================================================================
!+ Module procedure in "mo_sso_cosmo"
!------------------------------------------------------------------------------

SUBROUTINE gw_stress (                                  &
           ie     , ke1 , istart , iend   ,             &
           prho,pstab,pvph,psso_stdh,sso_lat_mask,psso_sigma,pdmod,  &
           lo_sso, pfi,mkenvh,ptau ,lacc )

!------------------------------------------------------------------------------
!
! Purpose: The module procedure gw_stress computes the gravity stress
!          amplitude following eq.4.11 of the documentation.
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------

!     Input
!     =====
!
      INTEGER, INTENT(IN) ::  &
      ie        ,    & ! number of grid points in first (zonal) direction
      ke1              ! ke + 1

      LOGICAL, OPTIONAL, INTENT(IN)  :: lacc ! If true, use openacc

      INTEGER, INTENT(IN) ::  &
      istart    ,    & ! start index for first (zonal) direction
      iend             ! end index for first (zonal) direction

      REAL(KIND=wp) :: prho (:,:) ! (ie,ke1)
      ! density on half levels    (kg/m**3)
      REAL(KIND=wp) :: pstab(:,:) ! (ie,ke1)
      ! squared Brunt-Vaisala frequency  (1/s**2)
      REAL(KIND=wp) :: pvph (:,:) ! (ie,ke1)
      ! wind on half levels projected on plane of low level wind (m/s)
      REAL(KIND=wp) :: psso_stdh (:) ! (ie)
      ! standard deviation of sso-height (m)
      REAL(KIND=wp) :: psso_sigma(:) ! (ie)
      ! mean sso-slope                   (-)
      REAL(KIND=wp) :: pdmod(:) ! (ie)
      ! projection parameter = SQRT(D1**2+D2**2)    cf. eq.4.7
      REAL(KIND=wp) :: pfi (:,:) ! (ie,ke) full-level geopotential
      INTEGER mkenvh (:) ! (ie) index of top of envelope layer
      REAL(KIND=wp) :: sso_lat_mask (:) ! (ie)
      LOGICAL lo_sso(:) ! (ie)
      LOGICAL :: lzacc ! non-optional version of lacc
      !

!     Output
!     ======

      REAL(KIND=wp) :: ptau(:,:) ! (ie,ke1)
      ! gravity wave stress amplitude   (Pa)

!     local variables
!     ===============
                              ! utility variables, which may be used to modify
      REAL(KIND=wp) :: zblock ! the magnitude of the subgrid scale standard
      REAL(KIND=wp) :: zeff   ! deviation which enters the stress amplitude
                              ! calculation
      REAL(KIND=wp) :: zramp  ! ramp function to avoid sudden jumps at sso stdev threshold (minsso_gwd)
      INTEGER j1  ! loop variable

      CALL set_acc_host_or_device(lzacc, lacc)

      !Declaration of GPU arrays
      !$ACC DATA PRESENT(ptau, lo_sso, prho, pstab, pvph, psso_stdh, psso_sigma, sso_lat_mask, pdmod, pfi, mkenvh) IF(lzacc)

!     gravity wave stress amplitude (eq.4.11)
!     =======================================
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO j1=istart,iend
        IF(lo_sso(j1)) THEN
          zblock=pfi(j1,mkenvh(j1))/g
          zeff=MAX(0._wp,3._wp*psso_stdh(j1)-zblock)
          zramp = MAX(0._wp, MIN(1._wp, (psso_stdh(j1)-minsso_gwd)/MAX(0.1_wp,minsso_gwd)))

          ptau(j1,ke1)=(Gkdrag+Gkdrag_add*sso_lat_mask(j1))*prho(j1,ke1)*psso_sigma(j1)*zeff**2*0.5_wp      &
     &                /psso_stdh(j1)*pvph(j1,ke1)*pdmod(j1)*SQRT(pstab(j1,ke1))*zramp
        ELSE
          ptau(j1,ke1)=0.0_wp
        ENDIF
      END DO
      !$ACC END PARALLEL

      !$ACC WAIT IF(lzacc)
      !$ACC END DATA

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE gw_stress

!============================================================================
!+ Module procedure in "mo_sso_cosmo"
!------------------------------------------------------------------------------

SUBROUTINE gw_profil(                                    &
           ie     , ke     , ke1 , istart , iend ,       &
           pph    , prho   , pstab  , pvph     , pri ,   &
           ptau   , pdmod  , psso_sigma, psso_stdh, sso_lat_mask   ,   &
           kkcrith, kcrit, kkenvh, kknu, kknu2 ,         &
           lo_sso , lacc)

!------------------------------------------------------------------------------
!
! Purpose: The module procedure gw_profil computes the vertical profile of
!          gravity wave stress.
!
! Method:  The stress profile for gravity waves is computed as follows:
!
!          - it is constant (no gravity wave drag) at all levels
!            between the ground and the top of the blocked layer kkenvh
!          - it decreases linearly with height from the top of the
!            blocked layer to 3 * stdh (kknu), to simulate lee waves or
!            non-linear gravity wave breaking
!          - at levels above (kknu) it is constant, except when the
!            wave encounters a critical level (kcrit) or when it breaks
!
!------------------------------------------------------------------------------
!
! Declarations:
!
!------------------------------------------------------------------------------

!     Input
!     =====
!
      INTEGER, INTENT(IN) ::  &
      ie        ,    & ! number of grid points in first (zonal) direction
      ke        ,    & ! number of grid points in vertical direction
      ke1              ! ke + 1

      LOGICAL, OPTIONAL, INTENT(IN)  :: lacc ! If true, use openacc

      INTEGER, INTENT(IN) ::  &
      istart    ,    & ! start index for first (zonal) direction
      iend             ! end index for first (zonal) direction

      REAL(KIND=wp) :: pph (:,:) !  (ie,ke1)
      ! half level pressure           (Pa)
      REAL(KIND=wp) :: prho (:,:) ! (ie,ke1)
      ! density on half levels    (kg/m**3)
      REAL(KIND=wp) :: pstab(:,:) ! (ie,ke1)
      ! squared Brunt-Vaisala frequency  (1/s**2)
      REAL(KIND=wp) :: pri  (:,:) ! (ie,ke1)
      ! mean flow Richardson number      ( - )
      REAL(KIND=wp) :: pvph (:,:) ! (ie,ke1)
      ! wind on half levels projected on plane of low level wind (m/s)
      REAL(KIND=wp) :: ptau (:,:) ! (ie,ke1)
      !gravity wave stress profile
      REAL(KIND=wp) :: pdmod(:) ! (ie)
      ! projection parameter = SQRT(D1**2+D2**2)    cf. eq.4.7
      REAL(KIND=wp) :: psso_stdh (:) ! (ie)
      ! standard deviation of sso-height (m)
      REAL(KIND=wp) :: psso_sigma(:) ! (ie)
      ! mean sso-slope                   (-)
      REAL(KIND=wp) :: sso_lat_mask (:) ! (ie)

!     various significant levels
      INTEGER kkcrith(:) ! (ie)
      INTEGER kcrit  (:) ! (ie)
      INTEGER kkenvh (:) ! (ie)
      INTEGER kknu   (:) ! (ie)
      INTEGER kknu2  (:) ! (ie)

      LOGICAL lo_sso (:) ! (ie)

!     local arrays and variables
!     ==========================

      REAL(KIND=wp) :: zdz2   (ie,ke)
      REAL(KIND=wp) :: ztau1(ie), ztau2(ie), ztau3(ie)
      REAL(KIND=wp) :: znorm  (ie)
      REAL(KIND=wp) :: zoro   (ie)

      REAL(KIND=wp) :: zb,zdelp,zdelpt,zalpha,zalfa,zdel,zriw    ! utitility variables
      REAL(KIND=wp) :: zsqri,zdz2n,rcrit                               ! utitility variables

      INTEGER j1,j3                 ! loop indices
      LOGICAL :: lzacc ! non-optional version of lacc

      CALL set_acc_host_or_device(lzacc, lacc)

      !Declaration of GPU arrays
      !$ACC DATA PRESENT(pph, prho, pstab, pri, pvph, ptau, pdmod, psso_stdh, psso_sigma) &
      !$ACC   PRESENT(kkcrith, kcrit, kkenvh, kknu, kknu2, lo_sso, sso_lat_mask) &
      !$ACC   CREATE(zdz2, ztau1, ztau2, ztau3, znorm, zoro) IF(lzacc)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO j1=istart,iend
        IF(lo_sso(j1)) THEN
          zoro(j1) = psso_sigma(j1)*pdmod(j1)           &
     &                 /(4._wp*MAX(psso_stdh(j1),1.0_wp))
          ztau1(j1) = ptau(j1,ke1)
        ENDIF
      END DO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO j3=ke,2,-1     ! vertical loop
        !$ACC LOOP GANG VECTOR
        DO j1=istart,iend
!     constant stress up to top of blocking layer
!     ===========================================  
          IF(lo_sso(j1)) THEN
            IF(j3.GE.kknu2(j1)) THEN
              ptau(j1,j3)=ztau1(j1)
            ENDIF
          ENDIF

!     wave displacement at next level
!     ===============================
          IF(lo_sso(j1)) THEN
            IF(j3.LT.kknu2(j1)) THEN
            znorm(j1)=(Gkdrag+Gkdrag_add*sso_lat_mask(j1))*prho(j1,j3)*SQRT(pstab(j1,j3))  &
                      *pvph(j1,j3)*zoro(j1)
            zdz2(j1,j3)=ptau(j1,j3+1)/MAX(znorm(j1),Gssec)
            ENDIF
          ENDIF

!     wave Richardson number, new wave displacement and stress
!     breaking evaluation and critical level
!     ========================================================

          IF(lo_sso(j1)) THEN
          IF(j3.LT.kknu2(j1)) THEN    ! only above blocking layer
            IF((ptau(j1,j3+1).LT.Gtsec).OR.(j3.LE.kcrit(j1))) THEN
            ptau(j1,j3)=0.0_wp    ! above critical levels
            ELSE
            zsqri=SQRT(pri(j1,j3))
            zalfa=SQRT(pstab(j1,j3)*zdz2(j1,j3))/pvph(j1,j3)
            zriw=pri(j1,j3)*(1._wp-zalfa)/(1._wp+zalfa*zsqri)**2
            rcrit=Grcrit+Grcrit_add*sso_lat_mask(j1)
            IF(zriw.LT.rcrit) THEN      ! breaking occurs
              zdel=4._wp/zsqri/rcrit+1._wp/rcrit**2+4._wp/rcrit
              zb=1._wp/rcrit+2._wp/zsqri
              zalpha=0.5_wp*(-zb+SQRT(zdel))
              zdz2n=(pvph(j1,j3)*zalpha)**2/pstab(j1,j3)
              ptau(j1,j3)=znorm(j1)*zdz2n
            ELSE
              ptau(j1,j3)=znorm(j1)*zdz2(j1,j3)
            ENDIF
            ptau(j1,j3)=MIN(ptau(j1,j3),ptau(j1,j3+1))
            ENDIF
          ENDIF
          ENDIF
        END DO

      END DO       ! end of vertical loop
      !$ACC END PARALLEL

!     reorganisation of stress profile, if breaking occurs at low levels
!     ==================================================================
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO j1=istart,iend
        IF(lo_sso(j1)) THEN
          ztau2(j1) = ptau(j1,kkenvh(j1))
          ztau3(j1) = ptau(j1,kkcrith(j1))
        ENDIF
      END DO
      !$ACC END PARALLEL

!     linear decrease between kkenvh and kkcrith
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO j3=1,ke      ! vertical loop
        DO j1=istart,iend
          IF(lo_sso(j1)) THEN
            IF(j3.GT.kkcrith(j1).AND.j3.LT.kkenvh(j1))THEN
              zdelp=pph(j1,j3)-pph(j1,kkenvh(j1))
              zdelpt=pph(j1,kkcrith(j1))-pph(j1,kkenvh(j1))
              ptau(j1,j3)=ztau2(j1) + (ztau3(j1)-ztau2(j1) ) * zdelp/zdelpt
            ENDIF
          ENDIF
        END DO
      END DO       ! end of vertical loop

      !$ACC WAIT IF(lzacc)
      !$ACC END DATA
!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE gw_profil


!------------------------------------------------------------------------------
! End of module mo_sso_cosmo
!------------------------------------------------------------------------------

END MODULE mo_sso_cosmo
