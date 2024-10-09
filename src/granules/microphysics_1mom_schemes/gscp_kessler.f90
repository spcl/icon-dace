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

! Description of *gscp_kessler*:
!   This module procedure calculates the rates of change of temperature, cloud
!   water, water vapor and rain due to cloud microphysical processes related
!   to the formation of grid scale precipitation. This includes the sedimentation
!   of rain and snow.
!   The precipitation fluxes at the surface are also calculated here.
!
! Method:
!   Prognostic one-moment bulk microphysical parameterization.
!   The sedimentation of rain and snow is computed implicitly.
!
! Reference   This is an adaption of subroutine kessler in file src_gscp.f90
!  of the COSMO-Model. Equation numbers refer to
!  Doms, Foerstner, Heise, Herzog, Raschendorfer, Schrodin, Reinhardt, Vogel
!    (September 2005): "A Description of the Nonhydrostatic Regional Model LM",
!
!------------------------------------------------------------------------------

MODULE gscp_kessler

!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:

!------------------------------------------------------------------------------
! Microphysical constants and variables
!------------------------------------------------------------------------------

USE, INTRINSIC :: iso_fortran_env, ONLY: wp => real64, &
                                         i4 => int32
USE mo_physical_constants, ONLY: r_v   => rv    , & !> gas constant for water vapour
                                 o_m_rdv        , & !! 1 - r_d/r_v
                                 rdv            , & !! r_d / r_v
                                 lh_v  => alv   , & !! latent heat of vapourization
                                 cpdr  => rcpd  , & !! (spec. heat of dry air at constant press)^-1
                                 cvdr  => rcvd  , & !! (spec. heat of dry air at const vol)^-1
                                 b3    => tmelt , & !! melting temperature of ice/snow
                                 t0    => tmelt     !! melting temperature of ice/snow

USE mo_lookup_tables_constants, ONLY: &
                                 b1    => c1es  , & !! constants for computing the sat. vapour
                                 b2w   => c3les , & !! pressure over water (l) and ice (i)
                                 b4w   => c4les     !!               -- " --
USE mo_satad,              ONLY: sat_pres_water     !! saturation vapor pressure w.r.t. water
USE mo_exception,          ONLY: message, message_text

!------------------------------------------------------------------------------

USE gscp_data, ONLY: &          ! all variables are used here

    zconst,   zbev,    zvz0r,                 &
    x7o8,                                     &
    zkcac,   zkphi1,    zkphi2,    zkphi3,    &
    x1o8,      x3o16,   iautocon

!==============================================================================

IMPLICIT NONE
PRIVATE

!------------------------------------------------------------------------------
!! Public subroutines
!------------------------------------------------------------------------------

PUBLIC :: kessler

!==============================================================================

CONTAINS

!==============================================================================
!> Module procedure "kessler" in "gscp_kessler" for computing effects of 
!!  grid scale precipitation including cloud water, cloud ice, rain and snow
!------------------------------------------------------------------------------

SUBROUTINE kessler  (             &
  nvec,ke,                           & !> array dimensions
  ivstart,ivend, kstart,             & !! optional start/end indicies
  idbg,                              & !! optional debug level
  zdt, dz,                           & !! numerics parameters
  t,p,rho,qv,qc,qr,                  & !! prognostic variables
  qc0,                               & !! cloud ice/water threshold for autoconversion
  prr_gsp,                           & !! surface precipitation rates
  qrsflux,                           & !  precipitation flux
  l_cv,                              &
  ldass_lhn,                         &
  ldiag_ttend,     ldiag_qtend     , &
  ddt_tend_t     , ddt_tend_qv     , &
  ddt_tend_qc    ,                   & !> ddt_tend_xx are tendencies
  ddt_tend_qr)!!    necessary for dynamics

!------------------------------------------------------------------------------
!>
!! Description:
!!   This module procedure calculates the rates of change of temperature, cloud
!!   water, cloud ice, water vapor, rain and snow due to cloud microphysical
!!   processes related to the formation of grid scale precipitation. The
!!   variables are updated in this subroutine. Rain and snow are prognostic
!!   variables. The precipitation fluxes at the surface are stored on the
!!   corresponding global fields.
!!
!! Method:
!!
!------------------------------------------------------------------------------
!! Declarations:
!!
!------------------------------------------------------------------------------
!! Modules used: These are declared in the module declaration section
!! -------------

!! Subroutine arguments:
!! --------------------

  INTEGER, INTENT(IN) ::  &
    nvec          ,    & !> number of horizontal points
    ke                     !! number of grid points in vertical direction

  INTEGER, INTENT(IN), OPTIONAL ::  &
    ivstart   ,    & !> optional start index for horizontal direction
    ivend     ,    & !! optional end index   for horizontal direction
    kstart    ,    & !! optional start index for the vertical index
    idbg             !! optional debug level

  REAL(KIND=wp), INTENT(IN) :: &
    zdt                    !> time step for integration of microphysics     (  s  )

  REAL(KIND=wp), INTENT(IN) :: &
    qc0              !> cloud ice/water threshold for autoconversion

  REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) ::      &   ! (ie,ke)
    dz              ,    & !> layer thickness of full levels                (  m  )
    rho             ,    & !! density of moist air                          (kg/m3)
    p                      !! pressure                                      ( Pa  )

    LOGICAL, INTENT(IN):: l_cv, &                   !! if true, cv is used instead of cp
      ldass_lhn

  LOGICAL, INTENT(IN), OPTIONAL :: &
    ldiag_ttend,         & ! if true, temperature tendency shall be diagnosed
    ldiag_qtend            ! if true, moisture tendencies shall be diagnosed

  REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT) ::   &   ! dim (ie,ke)
    t               ,    & !> temperature                                   (  K  )
    qv              ,    & !! specific water vapor content                  (kg/kg)
    qc              ,    & !! specific cloud water content                  (kg/kg)
    qr                     !! specific rain content                         (kg/kg)

  REAL(KIND=wp), INTENT(INOUT) :: &
       qrsflux(:,:)        !  precipitation flux

  REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) ::   &   ! dim (ie)
    prr_gsp                !> precipitation rate of rain, grid-scale        (kg/(m2*s))

  REAL(KIND=wp), DIMENSION(:,:), INTENT(OUT), OPTIONAL ::   &     ! dim (ie,ke)
    ddt_tend_t      , & !> tendency T                                       ( 1/s )
    ddt_tend_qv     , & !! tendency qv                                      ( 1/s )
    ddt_tend_qc     , & !! tendency qc                                      ( 1/s )
    ddt_tend_qr         !! tendency qr                                      ( 1/s )

  !! Local parameters: None, parameters are in module header, gscp_data or data_constants
  !! ----------------
  
  REAL    (KIND=wp), PARAMETER ::  &
    ! basic constants of the parameterization scheme
    zaau   = 1.0_wp/1000.0_wp,            & ! coef. for autoconversion
    zaac   = 1.72_wp,                     & ! coef. for accretion (neu)
    zbev   = 9.0_wp,                      & ! coef. for drop ventilation,
    ! to avoid illegal operations in power expressions
    znull = 1.E-20_wp

  !> Local scalars:
  !! -------------
  
  INTEGER (KIND=i4) :: &
    iv, k             !> loop indices

  REAL    (KIND=wp   ) :: z_heat_cap_r !! reciprocal of cpdr or cvdr (depending on l_cv)

  LOGICAL :: lldiag_ttend, lldiag_qtend

  INTEGER ::  &
    iv_start     ,    & !> start index for horizontal direction
    iv_end       ,    & !! end index for horizontal direction
    k_start      ,    & !! model level where computations start
    izdebug             !! debug level

  REAL (KIND=wp)   ::  &
    ztx   ,            & ! 
    zpx   ,            & !
    fsa3  ,            & !
    zspw  ,            & ! equilibrium vapour pressure over water
    zsa3  ,            & !
    zc3   ,            & !
    zc1c  ,            & !
    zc1   ,            & !
    zx    ,            & !
    zsrmax

  REAL (KIND=wp)   ::  &
    zsqvw ,            & !  specific humidity at water saturation
    zqvts ,            & !  qv-tendency in a layer
    zqcts ,            & !  qc-tendency in a layer
    zqrts ,            & !  qr-tendency in a layer
    ztts  ,            & !  t -tendency in a layer
    zswra ,            & !  autoconversion rate
    zswrk ,            & !  accretion rate
    zsrd                 !  evaporation rate

  REAL(KIND=wp), DIMENSION(nvec,ke) ::   &
    t_in          ,    & !> temperature                                   (  K  )
    qv_in         ,    & !! specific water vapor content                  (kg/kg)
    qc_in         ,    & !! specific cloud water content                  (kg/kg)
    qr_in                !! specific rain content                         (kg/kg)

  REAL (KIND=wp)   ::  &
    zpv          ,     & !
    zdtdh, zphi  ,     & !
    zqrk, lnzqrk ,     & !
    zdtr, ztau   ,     & ! 
    zimr, zzar   ,     & !
    qvg          ,     & !  specific water vapor content: local grid cell value
    qcg          ,     & !  specfic cloud water content:  - "" -
    tg           ,     & !  temperature:                  - "" -
    qrg          ,     & !  specific rain content:        - "" -
    rhog         ,     & !  density                       - "" -
    rhogr        ,     & !  reciprocal density            - "" -
    ppg                  !  full level pressure           - "" -


!! Local (automatic) arrays:
!! -------------------------

  REAL    (KIND=wp   ) ::  &
    zvzr        (nvec),     & !
    zpkr        (nvec),     & !
    zprvr       (nvec),     & !
    zpkm1r      (nvec)        !

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
!! Begin Subroutine kessler
!------------------------------------------------------------------------------

!> Statement functions
! -------------------

  fsa3(ztx) = 3.86E-3_wp - 9.41E-5_wp*(ztx-t0)

!------------------------------------------------------------------------------
!  Section 1: Initial setting of local variables
!------------------------------------------------------------------------------

! Define reciprocal of heat capacity of dry air (at constant pressure vs at constant volume)


    IF (l_cv) THEN
      z_heat_cap_r = cvdr
    ELSE
      z_heat_cap_r = cpdr
    ENDIF

  ! Delete precipitation fluxes from previous timestep
  prr_gsp (:) = 0.0_wp
  zpkr    (:) = 0.0_wp
  zprvr   (:) = 0.0_wp
  zvzr    (:) = 0.0_wp

  ! Optional arguments

  IF (PRESENT(ivstart)) THEN
    iv_start = ivstart
  ELSE
    iv_start = 1
  END IF
  IF (PRESENT(ivend)) THEN
    iv_end = ivend
  ELSE
    iv_end = nvec
  END IF
  IF (PRESENT(kstart)) THEN
    k_start = kstart
  ELSE
    k_start = 1
  END IF
  IF (PRESENT(idbg)) THEN
    izdebug = idbg
  ELSE
    izdebug = 0
  END IF
  IF (PRESENT(ldiag_ttend)) THEN
    lldiag_ttend = ldiag_ttend
  ELSE
    lldiag_ttend = .FALSE.
  ENDIF
  IF (PRESENT(ldiag_qtend)) THEN
    lldiag_qtend = ldiag_qtend
  ELSE
    lldiag_qtend = .FALSE.
  ENDIF

  ! save input arrays for final tendency calculation
  IF (lldiag_ttend) THEN
    t_in  = t
  ENDIF
  IF (lldiag_qtend) THEN
    qv_in = qv
    qc_in = qc
    qr_in = qr
  END IF

! timestep for calculations
  zdtr  = 1.0_wp / zdt

  ! add part of latent heating calculated in subroutine kessler to model latent
  ! heating field: subtract temperature from model latent heating field
  IF (ldass_lhn) THEN
    ! CALL get_gs_lheating ('add',1,ke) !XL :should not be called from block physics
    qrsflux(:,:) = 0.0_wp
  ENDIF

! output for various debug levels
  IF (izdebug > 15) CALL message('','gscp_kessler:  Start of kessler')
  IF (izdebug > 20) THEN
    WRITE (message_text,*) '   nvec = ',nvec       ; CALL message('',message_text)
    WRITE (message_text,*) '   ke = ',ke           ; CALL message('',message_text)
    WRITE (message_text,*) '   ivstart = ',ivstart ; CALL message('',message_text)
    WRITE (message_text,*) '   ivend   = ',ivend   ; CALL message('',message_text)
  END IF

  IF (izdebug > 50) THEN
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN dz  = ',MAXVAL(dz),MINVAL(dz)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN T   = ',MAXVAL(t),MINVAL(t)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN p   = ',MAXVAL(p),MINVAL(p)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN rho = ',MAXVAL(rho),MINVAL(rho)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qv  = ',MAXVAL(qv),MINVAL(qv)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qc  = ',MAXVAL(qc),MINVAL(qc)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qr  = ',MAXVAL(qr),MINVAL(qr)
    CALL message('',message_text)
  ENDIF

! ----------------------------------------------------------------------
! Loop from the top of the model domain to the surface to calculate the
! transfer rates  and sedimentation terms
! ----------------------------------------------------------------------

loop_over_levels: DO k = 1, ke

  !----------------------------------------------------------------------------
  !  Section 2: Test for clouds and precipitation in present layer.
  !             If no cloud or precipitation points have been found
  !             go to the next layer.
  !----------------------------------------------------------------------------

  loop_over_i: DO iv = iv_start, iv_end

    IF (qr(iv,k) < 1.0E-15_wp) qr(iv,k) = 0.0_wp
    IF (qc(iv,k) < 1.0E-15_wp) qc(iv,k) = 0.0_wp

    qcg   = qc(iv,k)
    qvg   = qv(iv,k)
    qrg   = qr(iv,k)
    tg    = t(iv,k)
    ppg   = p(iv,k)
    rhog  = rho(iv,k)
    rhogr = 1.0_wp / rhog

    zqrk  = qrg * rhog
    zdtdh = 0.5_wp * zdt / dz(iv,k)

    zpkm1r(iv) = zpkr(iv)
    zzar      = zqrk/zdtdh + zprvr(iv) + zpkm1r(iv)
    IF (zqrk > znull) THEN
      zpkr(iv) = zqrk * zvz0r * EXP(x1o8*LOG(zqrk))
    ELSE
      zpkr(iv) = 0.0_wp
    ENDIF
    zpkr(iv) = MIN( zpkr(iv), zzar )
    zzar    = zdtdh * (zzar-zpkr(iv))

    IF (zvzr(iv) == 0.0_wp) zvzr(iv) = zvz0r * EXP(x1o8*LOG(MAX(0.5_wp*zqrk,znull)))

    zimr    = 1.0_wp / (1.0_wp + zvzr(iv) * zdtdh)
    zqrk    = zzar*zimr

    IF (zqrk > znull) THEN
      lnzqrk = LOG(zqrk)
    ELSE
      lnzqrk = 0.0_wp
    ENDIF

    ztts  = 0.0_wp
    zqvts = 0.0_wp
    zqcts = 0.0_wp

    !------------------------------------------------------------------------
    !  Section 5: Calculation of cloud microphysics for cloud case
    !             ( qc > 0)
    !------------------------------------------------------------------------

    IF (qcg > 0.0_wp) THEN

      ! Calculate conversion rates

      IF (iautocon == 0) THEN
        ! Coefficients
        zc1c  = zaac *  EXP(lnzqrk*x7o8)
        zc1   = zaau + zc1c
        zx    = qcg / (1.0_wp + zc1*zdt)

        ! Conversion rates
        zswra  = zaau * zx
        zswrk  = zc1c * zx
      ELSEIF (iautocon == 1) THEN
        ! Seifert and Beheng (2001) autoconversion rate
        ! with constant cloud droplet number concentration cloud_num

        IF (qcg > 1.0E-6_wp) THEN
          ztau  = MIN(1.0_wp-qcg/(qcg+qrg),0.9_wp)
          zphi  = zkphi1 * ztau**zkphi2 * (1.0_wp - ztau**zkphi2)**3
          zswra = zconst * qcg*qcg*qcg*qcg &                      ! S_au   ! This formula is wrong. Missing cloud_num
                 * (1.0_wp + zphi/(1.0_wp - ztau)**2)                      ! which was removed from definition of zconst
          zphi  = (ztau/(ztau+zkphi3))**4
          zswrk = zkcac * qcg * qrg * zphi !* zrho1o2(i)          ! S_ac
        ELSE
          zswra = 0.0_wp  ! S_au                    
          zswrk = 0.0_wp  ! S_ac
        ENDIF
      ENDIF

      ! Store tendencies
      zqcts  = - zswra - zswrk
      zqrts  =   zswra + zswrk

      ! Update values
      qrg = MAX(0.0_wp,(zzar*rhogr + zqrts*zdt)*zimr)

    !------------------------------------------------------------------------
    !  Section 7: Calculation of cloud microphysics for
    !             precipitation case without cloud ( qc = 0 )
    !------------------------------------------------------------------------

    ELSEIF ( (zqrk) > 0.0_wp .AND. qcg <= 0.0_wp )    THEN

      zsqvw = sat_pres_water(tg)/(rhog * r_v *tg)

      ! Coefficients
      zsrmax = zzar*rhogr * zdtr
      zsa3   = fsa3(tg)
      zc3    = zsa3 * SQRT(zqrk) * (1.0_wp + zbev * EXP(lnzqrk*x3o16))

      ! Conversion rates
      zsrd   = -zc3 * (qvg - zsqvw)
      zsrd   = MIN (zsrmax, zsrd)

      ! Store tendencies
      zqvts =   zsrd
      ztts  = - lh_v*zsrd*z_heat_cap_r
      zqrts = - zsrd

      ! Update values
      qrg = MAX(0.0_wp,(zzar*rhogr + zqrts*zdt)*zimr)

    ENDIF

    !------------------------------------------------------------------------
    ! Section 8: Complete time step
    !------------------------------------------------------------------------

    IF ( k /= ke ) THEN
      ! Store precipitation fluxes and sedimentation velocities for the next level
      zprvr(iv) = qrg*rhog*zvzr(iv)
      zvzr(iv)  = zvz0r * EXP(x1o8 * LOG(MAX((qrg+qr(iv,k+1))*0.5_wp*rhog,znull)))
          ! for the latent heat nudging
          IF (ldass_lhn) THEN
            qrsflux(iv,k) = zprvr(iv)
            qrsflux(iv,k) = 0.5_wp*(qrsflux(iv,k)+zpkr(iv))
          ENDIF
    ELSE
      ! Precipitation flux at the ground
      prr_gsp(iv) = 0.5_wp * (qrg*rhog*zvzr(iv) + zpkr(iv))
      ! for the latent heat nudging
        IF (ldass_lhn) &
           qrsflux(iv,k) = prr_gsp(iv)
    ENDIF

    ! Update of prognostic variables or tendencies
    qr (iv,k) = qrg
!   qrs(iv,k) = qrg
    t  (iv,k) = t (iv,k) + ztts*zdt
    qv (iv,k) = MAX ( 0.0_wp, qv(iv,k) + zqvts*zdt )
    qc (iv,k) = MAX ( 0.0_wp, qc(iv,k) + zqcts*zdt )

  ENDDO loop_over_i

  IF (izdebug > 25) THEN
    ! Check for negative values
    DO iv = iv_start, iv_end
      IF (qr(iv,k) < 0.0_wp) THEN
        WRITE(message_text,'(a)') ' WARNING: kessler_pp, negative value in qr'
        CALL message('',message_text)
      ENDIF
      IF (qc(iv,k) < 0.0_wp) THEN
        WRITE(message_text,'(a)') ' WARNING: kessler_pp, negative value in qc'
        CALL message('',message_text)
      ENDIF
      IF (qv(iv,k) < 0.0_wp) THEN
        WRITE(message_text,'(a)') ' WARNING: kessler_pp, negative value in qv'
        CALL message('',message_text)
      ENDIF
    ENDDO
  ENDIF

ENDDO loop_over_levels

!------------------------------------------------------------------------------
! final tendency calculation for ICON
!
! Note: as soon as we have a new satad subroutine in ICON, this tendency
! calculation will be done in the k-loop and the original 3D variables wont
! be used to store the new values. Then we wont need the _in variables anymore.
!------------------------------------------------------------------------------

! calculated pseudo-tendencies

  IF ( lldiag_ttend ) THEN
    DO k=k_start,ke
      DO iv=iv_start,iv_end
        ddt_tend_t (iv,k) = (t (iv,k) - t_in (iv,k))*zdtr
     END DO
    END DO
  ENDIF
 
  IF ( lldiag_qtend ) THEN
    DO k=k_start,ke
      DO iv=iv_start,iv_end
        ddt_tend_qv(iv,k) = MAX(-qv_in(iv,k)*zdtr,(qv(iv,k) - qv_in(iv,k))*zdtr)
        ddt_tend_qc(iv,k) = MAX(-qc_in(iv,k)*zdtr,(qc(iv,k) - qc_in(iv,k))*zdtr)
        ddt_tend_qr(iv,k) = MAX(-qr_in(iv,k)*zdtr,(qr(iv,k) - qr_in(iv,k))*zdtr)
      END DO
    END DO
  ENDIF

  IF (izdebug > 25) THEN
    CALL message('mo_gscp', 'UPDATED VARIABLES')
   WRITE(message_text,'(a,2E20.9)') 'kessler_pp T= ' ,&
    MAXVAL( t(:,:)), MINVAL(t(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'kessler_pp qv= ',&
    MAXVAL( qv(:,:)), MINVAL(qv(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'kessler_pp qc= ',&
    MAXVAL( qc(:,:)), MINVAL(qc(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(a,2E20.9)') 'kessler_pp qr= ',&
    MAXVAL( qr(:,:)), MINVAL(qr(:,:) )
    CALL message('', TRIM(message_text))
  ENDIF

!------------------------------------------------------------------------------
! End of subroutine kessler
!------------------------------------------------------------------------------

END SUBROUTINE kessler

!==============================================================================

END MODULE gscp_kessler
