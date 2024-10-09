! mo_gausgrid.f90 - Calculate quantities related to Gaussian grids
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!_________________________________________

MODULE mo_gaussgrid

  USE mo_kind,      ONLY: dp
  USE mo_control,   ONLY: ngl, nhgl, nlon
  USE mo_constants, ONLY: a, api, omega


  IMPLICIT NONE

  PUBLIC

  ! Calculate quantities related to the Gaussian grid.
  !
  !     References:
  !  
  !     S.L. Belousov, Tables of normalized associated Legendre Polynomials, 
  !                    Pergamon Press (1962)
  !     P.N. Swarztrauber, On computing the points and weights for 
  !                        Gauss-Legendre quadrature,
  !                    SIAM J. Sci. Comput. Vol. 24 (3) pp. 945-954 (2002)
  !
  ! Luis Kornblueh, MPI, 2010-02-09, updated for more precise calculation of
  !                                  Gaussian latitudes and weights for 
  !                                  higher resolution (>= T106).
  !                                  The same update is scheduled for
  !                                  afterburner and cdo.
  !
  !  This routine originates (year 2014) from MPI-ESM, the Earth System Model of the 
  !  Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
  !  Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
  !  version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
  !  doi: 10.1029/2018MS001400.

  REAL(dp), ALLOCATABLE :: gl_gw(:)       ! Gaussian weights
  REAL(dp), ALLOCATABLE :: gl_gmu(:)      ! mu = sin(Gaussian latitudes)
  REAL(dp), ALLOCATABLE :: gl_coriol(:)   ! coriolis parameter, 2*omega*mu
  REAL(dp), ALLOCATABLE :: gl_twomu(:)    ! 2*mu
  REAL(dp), ALLOCATABLE :: gl_cst(:)      ! square of cos(latitude):(1-mu**2)
  REAL(dp), ALLOCATABLE :: gl_sqcst(:)    ! sqrt(1-mu**2)
  REAL(dp), ALLOCATABLE :: gl_rcsth(:)    ! half reciprocal of *cst*:1/(2*cst)
  REAL(dp), ALLOCATABLE :: gl_racst(:)    ! 1/(a*(1-mu**2)
  REAL(dp), ALLOCATABLE :: gl_rsqcst(:)   ! 1/sqrt(1-mu**2)
  REAL(dp), ALLOCATABLE :: gl_budw(:)     ! weights for global budgets,
                                          ! budw = gw/nlon
  REAL(dp), ALLOCATABLE :: gridarea(:)    ! area of a grid cell [m**2]
  REAL(dp), ALLOCATABLE :: philat(:)      ! gaussian latitude   [degrees]
  REAL(dp), ALLOCATABLE :: philon(:)      ! longitudes          [degrees]
  REAL(dp), ALLOCATABLE :: sinlon(:)      ! sin(longitude).
  REAL(dp), ALLOCATABLE :: coslon(:)      ! cos(longitude).

CONTAINS

  SUBROUTINE inigau

    ! Description:
    !
    ! Preset constants in *mo_gaussgrid*.
    !
    ! Method:
    !
    ! *inigau* is called from *setdyn*.
    !
    ! Authors:
    !
    ! M. Jarraud, ECMWF, December 1982, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! A. Rhodin, MPI, January 1999, subroutine inigau -> module mo_gaussgrid
    ! U. Schlese, MPI, January 2001, grid area added
    ! A. Rhodin, MPI, June 2001, philon, philat added
    ! L. Kornblueh, MPI, October 2001, provide non 'ping-pong' properties
    ! U. Schulzweida, MPI, May 2002, change 'ping-pong' to N->S
    ! 
    ! for more details see file AUTHORS
    !

    !  Local scalars: 
    REAL(dp) :: zcst, zl, zsqcst
    INTEGER :: jgl, jlon

    !  Local arrays: 
    REAL(dp) :: zgmu(ngl), zgw(ngl)

    !  Executable statements 

    !-- 0. Allocate module provided fields

    IF (.NOT. ALLOCATED(gl_gw)) THEN
      ALLOCATE(gl_gw(ngl))      
      ALLOCATE(gl_gmu(ngl))     
      ALLOCATE(gl_coriol(ngl))  
      ALLOCATE(gl_twomu(ngl))   
      ALLOCATE(gl_cst(ngl))     
      ALLOCATE(gl_sqcst(ngl))   
      ALLOCATE(gl_rcsth(ngl))   
      ALLOCATE(gl_racst(ngl))   
      ALLOCATE(gl_rsqcst(ngl))  
      ALLOCATE(gl_budw(ngl))    
      ALLOCATE(gridarea(ngl))    
      ALLOCATE(philat(ngl))    
      ALLOCATE(philon(nlon))    
      ALLOCATE(sinlon(2*nlon))    
      ALLOCATE(coslon(2*nlon))    
    END IF

    !-- 1. Compute Gaussian latitudes and weights

!lk old call - may allow for checking.
!lk CALL gauaw_old(zgmu,zgw,ngl)
    CALL gauaw(zgmu, zgw, ngl)

    DO jgl = 1, nhgl
      gl_gw(jgl)        = zgw(jgl)*0.5_dp
      gl_gw(ngl-jgl+1)  = zgw(jgl)*0.5_dp
      gl_gmu(jgl)       = zgmu(jgl)
      gl_gmu(ngl-jgl+1) = zgmu(jgl)
      
      !-- 2. Derive some other constants
      
      gl_coriol(jgl)       =  2*omega*zgmu(jgl)
      gl_coriol(ngl-jgl+1) = -2*omega*zgmu(jgl)
      gl_twomu(jgl)        =  2*zgmu(jgl)
      gl_twomu(ngl-jgl+1)  = -2*zgmu(jgl)
      gl_budw(jgl)         = gl_gw(jgl)/nlon
      gl_budw(ngl-jgl+1)   = gl_gw(jgl)/nlon
      zcst                 = 1.0_dp-zgmu(jgl)**2
      zsqcst               = SQRT(zcst)
      gl_cst(jgl)          = zcst
      gl_cst(ngl-jgl+1)    = zcst
      gl_sqcst(jgl)        = zsqcst
      gl_sqcst(ngl-jgl+1)  = zsqcst
      gl_rsqcst(jgl)       = 1.0_dp/zsqcst
      gl_rsqcst(ngl-jgl+1) = 1.0_dp/zsqcst
      gl_rcsth(jgl)        = 0.5_dp/zcst
      gl_rcsth(ngl-jgl+1)  = 0.5_dp/zcst
      gl_racst(jgl)        = 1.0_dp/(a*zcst)
      gl_racst(ngl-jgl+1)  = 1.0_dp/(a*zcst)
    END DO


    DO jlon = 1, nlon               ! double size for rotated domains 
      zl = 2*api*(jlon-1.0_dp)/nlon    ! on decomposed grid
      sinlon(jlon) = SIN(zl)
      sinlon(jlon+nlon) = sinlon(jlon)
      coslon(jlon) = COS(zl)
      coslon(jlon+nlon) = coslon(jlon)
      philon(jlon) = 360._dp*(jlon-1.0_dp)/nlon
    END DO

    !  Grid area stored from N - > S !

    DO jgl=1,nhgl      
      gridarea(jgl)       = gl_budw(jgl)*4*api*a**2
      gridarea(ngl+1-jgl) = gridarea(jgl)
      philat  (jgl)       = 180._dp/api*ASIN(zgmu(jgl))
      philat  (ngl-jgl+1) = -philat(jgl)
    END DO

  END SUBROUTINE inigau
!------------------------------------------------------------------------------
  SUBROUTINE gauaw_old (pa, pw, nlat)

    ! Description:
    !
    ! Compute abscissas and weights for gaussian integration.
    !
    ! Method:
    !

    !  Scalar arguments 
    INTEGER :: nlat

    !  Array arguments 
    REAL(dp) :: pa(nlat), pw(nlat)
    ! *pa*  - array, length at least *k,* to receive abscis abscissas.
    ! *pw*  - array, length at least *k,* to receive weights.


    !  Local scalars: 
    REAL(dp), PARAMETER :: epsil = EPSILON(0.0_dp)
    INTEGER, PARAMETER :: itemax = 20

    INTEGER :: iter, ins2, isym, jn, jgl
    REAL(dp):: za, zw, z, zan
    REAL(dp):: zk, zkm1, zkm2, zx, zxn, zldn, zmod

    !  Intrinsic functions 
    INTRINSIC ABS, COS, MOD, TAN

    !  Executable statements 

    ins2 = nlat/2+MOD(nlat,2)

    ! Find first approximation of the roots of the
    ! Legendre polynomial of degree nlat
    
    DO jgl = 1, ins2
       z = REAL(4*jgl-1,dp)*api/REAL(4*nlat+2,dp)
       pa(jgl) = COS(z+1.0_dp/(TAN(z)*REAL(8*nlat**2,dp)))
    END DO

    ! Computes roots and weights
    ! Perform the Newton loop
    ! Find 0 of Legendre polynomial with Newton loop

    DO jgl = 1, ins2

       za = pa(jgl)
    
       DO iter = 1, itemax+1
          zk = 0.0_dp

          ! Newton iteration step
    
          zkm2 = 1.0_dp
          zkm1 = za
          zx = za
          DO jn = 2, nlat
             zk = (REAL(2*jn-1,dp)*zx*zkm1-REAL(jn-1,dp)*zkm2)/REAL(jn,dp)
             zkm2 = zkm1
             zkm1 = zk
          END DO
          zkm1 = zkm2
          zldn = (REAL(nlat,dp)*(zkm1-zx*zk))/(1.0_dp-zx*zx)
          zmod = -zk/zldn
          zxn = zx+zmod
          zan = zxn
    
          ! computes weight
    
          zkm2 = 1.0_dp
          zkm1 = zxn
          zx = zxn
          DO jn = 2,nlat
             zk = (REAL(2*jn-1,dp)*zx*zkm1-REAL(jn-1,dp)*zkm2)/REAL(jn,dp)
             zkm2 = zkm1
             zkm1 = zk
          END DO
          zkm1 = zkm2
          zw = (1.0_dp-zx*zx)/(REAL(nlat*nlat,dp)*zkm1*zkm1)
          za = zan
          IF (ABS(zmod) <= epsil) EXIT
       END DO

       pa(jgl) = zan
       pw(jgl) = 2*zw
    
    ENDDO

!DIR$ IVDEP
!OCL NOVREC

    DO jgl = 1, nlat/2
       isym = nlat-jgl+1
       pa(isym) = -pa(jgl)
       pw(isym) = pw(jgl)
    ENDDO

  END SUBROUTINE gauaw_old
  !----------------------------------------------------------------------------
  SUBROUTINE gauaw(pl,pw, kn)

    REAL(dp), INTENT(out) :: pl(:) !< Abscissas of Gauss integration
    REAL(dp), INTENT(out) :: pw(:) !< Weights of the Gaussian integration
    INTEGER,  INTENT(in)  :: kn    !< Number of Gauss abscissas     

    REAL(dp) :: zfn(0:kn,0:kn), zfnlat(0:kn/2)
    INTEGER :: iter(kn)
    
    REAL(dp) :: z, zfnn
    INTEGER :: jgl, iodd, ik, jn, ins2, isym
    
    !
    ! 1.0 Initialize Fourier coefficients for ordinary Legendre polynomials
    !
    ! Belousov, Swarztrauber, and ECHAM use zfn(0,0) = sqrt(2)
    ! IFS normalisation chosen to be 0.5*Integral(Pnm**2) = 1 (zfn(0,0) = 2.0)
    !
    
    zfn(0,0) = SQRT(2.0_dp)
    DO jn = 1, kn
      zfnn = zfn(0,0)
      DO jgl = 1, jn
        zfnn = zfnn*SQRT(1.0_dp-0.25_dp/REAL(jgl**2,dp))
      ENDDO
      
      zfn(jn,jn) = zfnn
      
      iodd = MOD(jn, 2)
      DO jgl = 2, jn-iodd, 2
        zfn(jn,jn-jgl) = zfn(jn,jn-jgl+2) &
             *REAL((jgl-1)*(2*jn-jgl+2),dp)/REAL(jgl*(2*jn-jgl+1),dp)
      ENDDO
    ENDDO
    
    !
    ! 2.0 Gaussian latitudes and weights
    !
    
    iodd = MOD(kn, 2)
    ik = iodd
    DO jgl = iodd,kn, 2
      zfnlat(ik) = zfn(kn,jgl)
      ik = ik+1
    ENDDO
    
    ins2 = kn/2+MOD(kn,2)
    
    !
    ! 2.1 Find first approximation of the roots of the
    !     Legendre polynomial of degree kn.
    !
    
    DO jgl = 1, ins2
      z = REAL(4*jgl-1,dp)*api/REAL(4*kn+2,dp)
      pl(jgl) = z+1.0_dp/(TAN(z)*REAL(8*kn**2,dp))
    ENDDO
    
    ! 2.2 Computes roots and weights for transformed theta
    
    DO jgl = ins2, 1, -1
      CALL gawl(zfnlat, pl(jgl), pw(jgl), kn, iter(jgl))
    ENDDO
    
    ! convert to physical latitude
    
    DO jgl = 1, ins2
      pl(jgl)  = COS(pl(jgl))
    ENDDO
    
    DO jgl = 1, kn/2
      isym = kn-jgl+1
      pl(isym) = -pl(jgl)
      pw(isym) = pw(jgl)
    ENDDO
    
  CONTAINS
    
    SUBROUTINE gawl(pfn, pl, pw, kn, kiter)
      
      INTEGER,  INTENT(in)    :: kn
      REAL(dp), INTENT(in)    :: pfn(0:kn/2)

      REAL(dp), INTENT(out)   :: pw
      INTEGER,  INTENT(out)   :: kiter

      REAL(dp), INTENT(inout) :: pl
      
      REAL(dp) :: pmod
      
      INTEGER :: iflag, itemax, jter, iodd
      REAL(dp) :: zw
      REAL(dp) :: zdlx,zdlxn
      
      ! 1.0 Initialization.
      
      iflag  =  0
      itemax = 20
      
      iodd   = MOD(kn, 2)
      
      zdlx   = pl
      
      ! 2.0 Newton iteration
      
      DO jter = 1, itemax+1
        kiter = jter
        CALL cpledn(kn, iodd, pfn, zdlx, iflag, zw, zdlxn, pmod)
        zdlx = zdlxn
        IF (iflag == 1) EXIT
        IF (ABS(pmod) <= EPSILON(zw)*1000.0_dp) iflag = 1
      ENDDO
      
      pl  = zdlxn
      pw  = zw
      
    END SUBROUTINE gawl
    
    SUBROUTINE cpledn(kn, kodd, pfn, pdx, kflag, pw, pdxn, pxmod)
      
      INTEGER,  INTENT(in) :: kn
      INTEGER,  INTENT(in) :: kodd
      REAL(dp), INTENT(in) :: pfn(0:kn/2)
      REAL(dp), INTENT(in) :: pdx
      INTEGER,  INTENT(in) :: kflag
      
      REAL(dp), INTENT(out) :: pw
      REAL(dp), INTENT(out) :: pxmod

      REAL(dp), INTENT(inout) :: pdxn
      
      REAL(dp) :: zdlx, zdlk, zdlldn, zdlxn, zdlmod
      
      INTEGER :: ik, jn
      
      ! 1.0 Newton iteration step
      
      zdlx = pdx
      zdlk = 0.0_dp
      IF (kodd == 0) zdlk = 0.5_dp*pfn(0)
      zdlxn  = 0.0_dp 
      zdlldn = 0.0_dp
      
      ik = 1
      
      IF (kflag == 0) THEN
        DO jn = 2-kodd, kn, 2
          ! normalised ordinary Legendre polynomial == \overbar{p_n}^0
          zdlk = zdlk + pfn(ik)*COS(REAL(jn,dp)*zdlx)
          ! normalised derivative == d/d\theta(\overbar{p_n}^0)
          zdlldn = zdlldn - pfn(ik)*REAL(jn,dp)*SIN(REAL(jn,dp)*zdlx)
          ik = ik+1
        ENDDO
        ! Newton method
        zdlmod = -zdlk/zdlldn
        zdlxn  = zdlx+zdlmod
        pdxn   = zdlxn
        pxmod  = zdlmod
      ENDIF
      
      ! 2.0 Computes weight
      
      IF (kflag == 1) THEN
        DO jn = 2-kodd, kn, 2
          ! normalised derivative
          zdlldn = zdlldn - pfn(ik)*REAL(jn,dp)*SIN(REAL(jn,dp)*zdlx)
          ik = ik+1
        ENDDO
        pw = REAL(2*kn+1,dp)/zdlldn**2
      ENDIF
      
    END SUBROUTINE cpledn
    
  END SUBROUTINE gauaw
  !----------------------------------------------------------------------------
  SUBROUTINE cleanup_gaussgrid
    IF (ALLOCATED(gl_gw)) THEN
      DEALLOCATE (gl_gw)
      DEALLOCATE (gl_gmu)
      DEALLOCATE (gl_coriol)
      DEALLOCATE (gl_twomu)
      DEALLOCATE (gl_cst)
      DEALLOCATE (gl_sqcst)
      DEALLOCATE (gl_rcsth)
      DEALLOCATE (gl_racst)
      DEALLOCATE (gl_rsqcst)
      DEALLOCATE (gl_budw)
      DEALLOCATE (gridarea)
      DEALLOCATE (philat)
      DEALLOCATE (philon)
      DEALLOCATE (sinlon)
      DEALLOCATE (coslon)
    END IF
  END SUBROUTINE cleanup_gaussgrid
!------------------------------------------------------------------------------
END MODULE mo_gaussgrid
