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

MODULE mo_cloud_diag

!------------------------------------------------------------------------------

USE mo_kind, ONLY: &
    ireals   =>wp       , &
    iintegers=>i4

!USE mo_exception,          ONLY: message, finish ,message_text
USE mo_physical_constants, ONLY: &
!    rvd_m_o  =>vtmpc1   , & !! r_v/r_d - 1
!    o_m_rdv             , & !! 1 - r_d/r_v
!    rdv                 , & !! r_d / r_v
    lhocp    => alvdcp  !, & !! lh_v/cp_d
!    b3       => tmelt   , & !! melting temperature of ice/snow

USE data_turbulence, ONLY: &
    clc_diag            , & !! cloud cover at saturation in statistical cloud diagnostic
    q_crit                  !! critical value for normalized over-saturation

USE mo_physical_constants, ONLY: &
    uc1                 , & !! variables for computing the rate of cloud cover in
    ucl                     !! the unsaturated case
USE mo_math_constants, ONLY: &
    uc2 => sqrt3            !!               -- " --

!USE mo_convect_tables, ONLY: &
!    b1       => c1es    , & !! constants for computing the sat. vapour
!    b2w      => c3les   , & !! pressure over water (l) and ice (i)
!    b2i      => c3ies   , & !!               -- " --
!    b4w      => c4les   , & !!               -- " --
!    b4i      => c4ies   , & !!               -- " --
!    b234w    => c5les       !!               -- " --

USE mo_satad, ONLY: zpsat_w => sat_pres_water, & !! saturation vapor pressure w.r.t. water
!                   zpsat_i => sat_pres_ice  , & !! saturation vapor pressure w.r.t. ice
                    zqvap   => spec_humi     , & !! Specific humidity
                    zdqsdt  => dqsatdT           !! Derivation of qsat w.r.t. temperature


!------------------------------------------------------------------------------

IMPLICIT NONE

PRIVATE


PUBLIC :: cloud_diag, turb_cloud

!------------------------------------------------------------------------------

CONTAINS


SUBROUTINE cloud_diag ( clc, clwc,                         &
                        iis, iie, ijs, ije, iks, ike,      &
                        ie , je , ke,                      &
                        t, qv, qc, pp, p0, rcld, ps,       &
                        itype_wcld )

!------------------------------------------------------------------------------
!
! Description:
!
!     This routine calculates the area fraction of a grid box covered
!     by stratiform (non-convective) clouds.
!     If subgrid-scale condensation is required, an additional
!     saturation adjustment is done.
!
! Method:
!
!     itype_wcld = 1 :
!     The fractional cloud cover clc is determined empirically from
!     relative humidity. Also, an in-cloud water content of sugrid-scale
!     clouds is determined as a fraction of the saturation specific
!     humidity. Both quantities define the grid-volume mean cloud water
!     content.
!     itype_wcld=2:
!     A Gaussion distribution is assumed for the saturation deficit
!     dq = qt - qs where qt = qv + ql is the total water content and
!     qs is the saturation specific humidity. Using the standard deviation
!     rcld of this distribution (on input) and the conservative grid-scale
!     quantities qt and tl (liquid water temperature), a corrected liquid
!     water content is determined which contains also the contributions from
!     subgrid-scale clouds. A corresponding cloudiness is also calculated.
!
!------------------------------------------------------------------------------

! Subroutine arguments
!----------------------

! Scalar arguments with intent(in):

INTEGER (KIND=iintegers), INTENT (IN) :: &  ! dimensions and run indices
  ie              ,    & ! number of grid points in zonal direction
  je              ,    & ! number of grid points in meridional direction
  ke              ,    & ! number of grid points in vertical direction
  iis             ,    & ! ie start index
  ijs             ,    & ! je start index
  iks             ,    & ! ke start index
  iie             ,    & ! ie end index
  ije             ,    & ! je end index
  ike             ,    & ! ke end index
  itype_wcld             ! cloud cover method, 1: RH based, 2: Gaussian satur. deficit

! Array arguments with intent(in):

REAL (KIND=ireals), INTENT (IN)        :: & !
  t   (ie,je,ke ) ,    & ! temperature (main levels)
  qv  (ie,je,ke ) ,    & ! water vapour (")
  qc  (ie,je,ke ) ,    & ! cloud water  (")
  pp  (ie,je,ke ) ,    & ! perturbation pressure (") /for ICON total pressure
  p0  (ie,je,ke ) ,    & ! base state pressure (")   /for ICON 0
  rcld(ie,je,ke+1),    & ! standard deviation of saturation deficit
  ps  (ie,je)            ! surface pressure

! Array arguments with intent(out):

REAL (KIND=ireals), INTENT (OUT)        :: &
  clc (ie,je,ke)  ,    & ! stratiform subgrid-scale cloud cover
  clwc(ie,je,ke)         ! liquid water content of ""


! Local variables and constants
! -----------------------------

INTEGER (KIND=iintegers) :: &
  i,j,k                             ! loop indices

REAL (KIND=ireals), PARAMETER :: &
  zsig_max = 1.0E-3_ireals,       & ! max. standard deviation of saturation deficit
  zclwfak  = 0.005_ireals,        & ! fraction of saturation specific humidity
  zuc      = 0.95_ireals            ! constant for critical relative humidity

REAL (KIND=ireals)       :: &
  temp, pres, ql, qt, qs, tl, dq, & !
  gam, q, sig, uc,                & !
  zsigma, zclc1, zq_max             !

!REAL (KIND=ireals)       ::      &
! zpsat_w, zqvap, zdqsdt,         & !statement functions and
! zpvap, zqsat, ztemp, zpres        !their formal arguments

!------------ End of header ---------------------------------------------------

! Definition of statement functions:

! saturation vapour pressure over water (zpsat_w) and over ice (zpsat_i):
! zpsat_w(ztemp) = b1 * exp( b2w*(ztemp-b3)/(ztemp-b4w) )
! zpsat_i(ztemp) = b1 * exp( b2i*(ztemp-b3)/(ztemp-b4i) )

! specific humidity:
! zqvap(zpvap,zpres) = rdv * zpvap / ( zpres - o_m_rdv*zpvap )

! Derivation of zqsat with respect to temperature:
! zdqsdt(ztemp,zqsat) = b234w * ( 1.0_ireals + rvd_m_o*zqsat ) * zqsat &
!                            / (ztemp-b4w)**2

! Begin Subroutine cloud_diag
! ---------------------------

  zq_max   = q_crit*(1.0_ireals/clc_diag - 1.0_ireals)

  DO k = iks, ike
    DO j = ijs, ije
      DO i = iis, iie

        ql   = qc(i,j,k)               ! cloud water content
        qt   = ql + qv(i,j,k)          ! total water content
        pres = p0(i,j,k) + pp(i,j,k)   ! pressure
        temp = t(i,j,k)                ! temperature
        tl   = temp - lhocp*ql         ! liquid water temperature
        qs   = zqvap(zpsat_w(tl),pres) ! saturation mixing ratio
        dq   = qt - qs                 ! saturation deficit
        gam  = 1.0_ireals / ( 1.0_ireals + lhocp*zdqsdt(tl,qs) )

        IF ( itype_wcld == 1 ) THEN

        ! Calculation of cloud cover and cloud water content
        ! using an empirical relative humidity criterion

          zsigma = pres / ps(i,j)

          ! critical relative humidity
          uc     = zuc - uc1 * zsigma * ( 1.0_ireals - zsigma )  &
                   * ( 1.0_ireals + uc2*(zsigma-0.5_ireals) )

          ! cloud cover
          clc(i,j,k) = MAX( 0.0_ireals,  &
                       MIN( 1.0_ireals, clc_diag * ((qt/qs-uc)/(ucl-uc))) )**2

          ! in-cloud water content
          ql = qs * zclwfak

          ! grid-volume water content
          IF ( dq > 0.0_ireals ) THEN
            zclc1 = clc_diag * ( (1.0_ireals-uc)/(ucl-uc) )**2
            ql    = ql + (gam*dq-ql)*(clc(i,j,k)-zclc1)/(1.0_ireals-zclc1)
          END IF
          ql = clc(i,j,k) * ql

        ELSEIF ( itype_wcld == 2 ) THEN

        ! Statistical calculation of cloud cover and cloud water content
        ! using the standard deviation of the saturation deficit

          sig = MIN ( zsig_max, rcld(i,j,k) )

          ! in case of sig=0, the method is similar to grid-scale
          ! saturation adjustment. Otherwise, a fractional cloud cover
          ! is diagnosed.
          IF ( sig <= 0.0_ireals ) THEN
            clc(i,j,k) = ABS ( (SIGN(1.0_ireals,dq)+1.0_ireals)*0.5_ireals )
            ql         = clc(i,j,k) * gam * dq
          ELSE
            q          = dq / sig
            clc(i,j,k) = MIN ( 1.0_ireals, MAX ( 0.0_ireals, &
                                        clc_diag * (1.0_ireals+q/q_crit) ) )
            IF ( q <= - q_crit ) THEN
              ql = 0.0_ireals
            ELSEIF ( q >= zq_max ) THEN
              ql = gam * dq
            ELSE
              ql = gam * sig * (q+q_crit) * (q+zq_max) / (2._ireals*(q_crit+zq_max))
            ENDIF
          ENDIF

        ENDIF

        clwc(i,j,k) = ql

     ENDDO
   ENDDO
 ENDDO

END SUBROUTINE cloud_diag


SUBROUTINE turb_cloud (                      &
!
   ie, je, ke, ke1,            kcs,    kce,  &
   istart, iend, jstart, jend, kstart, kend, &
!
   t, qv, qc, prs, ps, rcld,                 &
!
   clc, clwc ,itype_wcld                       )
   

!------------------------------------------------------------------------------
!
! Description:
!
!     This routine calculates the area fraction of a grid box covered
!     by stratiform (non-convective) clouds.
!     If subgrid-scale condensation is required, an additional
!     saturation adjustment is done.
!
! Method:
!
!     itype_wcld = 1 :
!     The fractional cloud cover clc is determined empirically from
!     relative humidity. Also, an in-cloud water content of sugrid-scale
!     clouds is determined as a fraction of the saturation specific
!     humidity. Both quantities define the grid-volume mean cloud water
!     content.
!     itype_wcld=2:
!     A Gaussion distribution is assumed for the saturation deficit
!     dq = qt - qs where qt = qv + ql is the total water content and
!     qs is the saturation specific humidity. Using the standard deviation
!     rcld of this distribution (on input) and the conservative grid-scale
!     quantities qt and tl (liquid water temperature), a corrected liquid
!     water content is determined which contains alse the contributions from
!     subgrid-scale clouds. A corresponding cloudiness is also calculated.
!
!------------------------------------------------------------------------------

! Subroutine arguments
!----------------------

! Scalar arguments with intent(in):

INTEGER (KIND=iintegers), INTENT(IN) :: &  ! indices used for allocation of arrays
  ie, je,            & ! number of grib points of horizontal dimensions
  ke, ke1,           & ! number of full and half level dimension
  kcs, kce,          & ! lower and upper level index bound for output variabls
  istart, iend, & ! zonal      start and end index
  jstart, jend, & ! meridional start and end index
  kstart, kend, & ! vertical   start and end index
  itype_wcld      ! cloud cover method, 1: RH based, 2: Gaussian satur. deficit

! Array arguments with intent(in):

REAL (KIND=ireals), INTENT(IN) :: & !
  t   (ie,je,ke ),    &    ! temperature (main levels)
  qv  (ie,je,ke ),    &    ! water vapour (")
  qc  (ie,je,ke ),    &    ! cloud water  (")
  prs (ie,je,ke ),    &    ! base state pressure (")
  rcld(ie,je,ke1),    &    ! standard deviation of saturation deficit
  ps  (ie,je)              ! surface pressure

! Array arguments with intent(out):

REAL (KIND=ireals), INTENT(OUT) :: &
  clc (ie,je,kcs:kce),  & ! stratiform subgrid-scale cloud cover
  clwc(ie,je,kcs:kce)     ! liquid water content of ""

! Local variables and constants
! -----------------------------

INTEGER (KIND=iintegers) :: &
  i,j,k                   ! loop indices

REAL (KIND=ireals), PARAMETER :: &
  zsig_max = 1.0E-3_ireals,  & ! max. standard deviation of saturation deficit
  zclwfak  = 0.005_ireals,   & ! fraction of saturation specific humidity
  zuc      = 0.95_ireals       ! constant for critical relative humidity

REAL (KIND=ireals) :: &
  tord, pres, ql, qt, qs, tl, dq, & !
  gam, q, sig, uc,                & !
  zsigma, zclc1, zq_max             !

!------------ End of header ----------------------------------------

! Begin Subroutine turb_cloud
! ---------------------------

  zq_max   = q_crit*(1.0_ireals/clc_diag - 1.0_ireals)

  DO k = kstart, kend
     DO j = jstart, jend
     DO i = istart, iend

        ql   = qc(i,j,k)               ! cloud water content
        qt   = ql + qv(i,j,k)          ! total water content
        pres = prs(i,j,k)              ! pressure
        tord = t(i,j,k)                ! ordinary temperature
        tl   = tord - lhocp*ql         ! liquid water temperature
        qs   = zqvap(zpsat_w(tl),pres) ! saturation mixing ratio
        dq   = qt - qs                 ! saturation deficit
        gam  = 1.0_ireals / ( 1.0_ireals + lhocp*zdqsdt(tl,qs) )

        IF ( itype_wcld == 1 ) THEN

        ! Calculation of cloud cover and cloud water content
        ! using an empirical relative humidity criterion

          zsigma = pres / ps(i,j)

          ! critical relative humidity
          uc     = zuc - uc1 * zsigma * ( 1.0_ireals - zsigma )  &
                   * ( 1.0_ireals + uc2*(zsigma-0.5_ireals) )

          ! cloud cover
          clc(i,j,k) = MAX( 0.0_ireals,  &
                       MIN( 1.0_ireals, clc_diag * ((qt/qs-uc)/(ucl-uc))) )**2

          ! in-cloud water content
          ql = qs * zclwfak

          ! grid-volume water content
          IF ( dq > 0.0_ireals ) THEN
            zclc1 = clc_diag * ( (1.0_ireals-uc)/(ucl-uc) )**2
            ql    = ql + (gam*dq-ql)*(clc(i,j,k)-zclc1)/(1.0_ireals-zclc1)
          END IF
          ql = clc(i,j,k) * ql

        ELSEIF ( itype_wcld == 2 ) THEN

        ! Statistical calculation of cloud cover and cloud water content
        ! using the standard deviation of the saturation deficit

          sig = MIN ( zsig_max, rcld(i,j,k) )

          ! in case of sig=0, the method is similar to grid-scale
          ! saturation adjustment. Otherwise, a fractional cloud cover
          ! is diagnosed.
          IF ( sig <= 0.0_ireals ) THEN
            clc(i,j,k) = ABS ( (SIGN(1.0_ireals,dq)+1.0_ireals)*0.5_ireals )
            ql         = clc(i,j,k) * gam * dq
          ELSE
            q          = dq / sig
            clc(i,j,k) = MIN ( 1.0_ireals, MAX ( 0.0_ireals, &
                                        clc_diag * (1.0_ireals+q/q_crit) ) )
            IF ( q <= - q_crit ) THEN
              ql = 0.0_ireals
            ELSEIF ( q >= zq_max ) THEN
              ql = gam * dq
            ELSE
              ql = gam * sig * (q+q_crit) * (q+zq_max) / (2.0_ireals*(q_crit+zq_max))
            ENDIF
          ENDIF

        ENDIF

        clwc(i,j,k) = ql

     ENDDO
     ENDDO
  ENDDO

END SUBROUTINE turb_cloud


!==============================================================================

END MODULE mo_cloud_diag
