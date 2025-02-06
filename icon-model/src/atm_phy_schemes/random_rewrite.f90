! ---------------------------------------------------------------
! Copyright (C) Alan Miller, ACM
!
! Author: Alan Miller
! Contact: https://jblevins.org/mirror/amiller/
!          https://www.acm.org/publications
! This software is not licensed for commercial use.
! SPDX-License-Identifier: LicenseRef-ACM
! ---------------------------------------------------------------
! Code has been modified
! ---------------------------------------------------------------
MODULE random_rewrite

  USE mo_kind, ONLY : JPRB=>wp ,JPIM=>i4
  IMPLICIT NONE

! A module for random number generation from the following distributions:
!
!     Distribution                    Function/subroutine name
!
!     Normal (Gaussian)               random_normal
!     Gamma                           random_gamma
!     Chi-squared                     random_chisq
!     Exponential                     random_exponential
!     Weibull                         random_Weibull
!     Beta                            random_beta
!     t                               random_t
!     Multivariate normal             random_mvnorm
!     Generalized inverse Gaussian    random_inv_gauss
!     Poisson                         random_Poisson
!     Binomial                        random_binomial1   *
!                                     random_binomial2   *
!     Negative binomial               random_neg_binomial
!     von Mises                       random_von_Mises
!     Cauchy                          random_Cauchy
!
!  Generate a random ordering of the integers 1 .. N
!                                     random_order
!     Initialize (seed) the uniform random number generator for ANY compiler
!                                     seed_random_number

!     Lognormal - see note below.

!  ** Two functions are provided for the binomial distribution.
!  If the parameter values remain constant, it is recommended that the
!  first function is used (random_binomial1).   If one or both of the
!  parameters change, use the second function (random_binomial2).

! The compilers own random number generator, SUBROUTINE RANDOM_NUMBER(r),
! is used to provide a source of uniformly distributed random numbers.

! N.B. At this stage, only one random number is generated at each call to
!      one of the functions above.

! The module uses the following functions which are included here:
! bin_prob to calculate a single binomial probability
! lngamma  to calculate the logarithm to base e of the gamma function

! Some of the code is adapted from Dagpunar's book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
!
! In most of Dagpunar's routines, there is a test to see whether the value
! of one or two floating-point parameters has changed since the last call.
! These tests have been replaced by using a logical variable FIRST.
! This should be set to .TRUE. on the first call using new values of the
! parameters, and .FALSE. if the parameter values are the same as for the
! previous call.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Lognormal distribution
! If X has a lognormal distribution, then log(X) is normally distributed.
! Here the logarithm is the natural logarithm, that is to base e, sometimes
! denoted as ln.  To generate random variates from this distribution, generate
! a random deviate from the normal distribution with mean and variance equal
! to the mean and variance of the logarithms of X, then take its exponential.

! Relationship between the mean & variance of log(X) and the mean & variance
! of X, when X has a lognormal distribution.
! Let m = mean of log(X), and s^2 = variance of log(X)
! Then
! mean of X     = exp(m + 0.5s^2)
! variance of X = (mean(X))^2.[exp(s^2) - 1]

! In the reverse direction (rarely used)
! variance of log(X) = log[1 + var(X)/(mean(X))^2]
! mean of log(X)     = log(mean(X) - 0.5var(log(X))

! N.B. The above formulae relate to population parameters; they will only be
!      approximate if applied to sample values.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Version 1.13, 2 October 2000
! Changes from version 1.01
! 1. The random_order, random_Poisson & random_binomial routines have been
!    replaced with more efficient routines.
! 2. A routine, seed_random_number, has been added to seed the uniform random
!    number generator.   This requires input of the required number of seeds
!    for the particular compiler from a specified I/O unit such as a keyboard.
! 3. Made compatible with Lahey's ELF90.
! 4. Marsaglia & Tsang algorithm used for random_gamma when shape parameter > 1.
! 5. INTENT for array f corrected in random_mvnorm.
!     Author: Alan Miller
!             CSIRO Division of Mathematical & Information Sciences
!             Private Bag 10, Clayton South MDC
!             Clayton 3169, Victoria, Australia
!     Phone: (+61) 3 9545-8016      Fax: (+61) 3 9545-8080

REAL(KIND=JPRB), PRIVATE      :: zero = 0.0_JPRB, half = 0.5_JPRB, one = 1.0_JPRB
                      
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)

PUBLIC             :: random_Poisson, random_normal

CONTAINS


SUBROUTINE random_normal(streammax,rnd,idx,fn_val)

! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.
INTEGER(KIND=JPIM), INTENT(IN)    :: streammax
REAL(KIND=JPRB)   , INTENT(IN)    :: rnd(streammax)
INTEGER(KIND=JPIM), INTENT(INOUT) :: idx
REAL(KIND=JPRB),    INTENT(OUT)   :: fn_val

!     Local variables
REAL(KIND=JPRB) :: s = 0.449871_JPRB, t = -0.386595_JPRB, a = 0.19600_JPRB, b = 0.25472_JPRB,  &
            r1 = 0.27597_JPRB, r2 = 0.27846_JPRB, u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

DO
! Replace random_number call with random numbers passed through as input
! Random numbers to use in generation of new clouds
!   CALL RANDOM_NUMBER(u)
   u=rnd(idx)
   !safety! occasionally random numbers are exactly zero -> log function lower down can't cope
   if (u .eq. 0.) then
      idx=idx+1
      u=rnd(idx)
   endif

  idx=idx+1
!   CALL RANDOM_NUMBER(v)
  v=rnd(idx)
  if (v .eq. 0.) then
     idx=idx+1
     v=rnd(idx)
  endif
  idx=idx+1
!  if (idx.ge.5000) write(7,*) 'random_Poisson, ran out of idx, random_normal',idx
!  flush(7)

  v = 1.7156_JPRB * (v - half)

!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.0_JPRB*LOG(u)*u**2) EXIT
END DO

!     Return ratio of P's coordinates as the normal deviate
fn_val = v/u
END SUBROUTINE random_normal

SUBROUTINE random_exponential(rnd,idx,streammax,fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A NEGATIVE EXPONENTIAL DlSTRIBUTION WlTH DENSITY PROPORTIONAL
! TO EXP(-random_exponential), USING INVERSION.
INTEGER(KIND=JPIM), INTENT(IN)    :: streammax
REAL(KIND=JPRB),    INTENT(IN)    :: rnd(streammax)
INTEGER(KIND=JPIM), INTENT(INOUT) :: idx
REAL(KIND=JPRB),    INTENT(OUT)   :: fn_val

!     Local variable
REAL(KIND=JPRB) :: r

DO
! Replace random number call with random number passed into routine
  r=rnd(idx)
  if (r .eq. 0.) then
     idx=idx+1
     r=rnd(idx)
  endif
  idx=idx+1
!  if (idx.ge.5000) write(7,*) 'random_Poisson, ran out of idx, random_exponential',idx
!  flush(7)
!  CALL RANDOM_NUMBER(r)
  IF (r > zero) EXIT
END DO

fn_val = -LOG(r)

END SUBROUTINE random_exponential

FUNCTION random_Poisson(mu,idx,streammax,rnd) RESULT(ival)
!**********************************************************************
!     Translated to Fortran 90 by Alan Miller from:
!                           RANLIB
!
!     Library of Fortran Routines for Random Number Generation
!
!                    Compiled and Written by:
!
!                         Barry W. Brown
!                          James Lovato
!
!             Department of Biomathematics, Box 237
!             The University of Texas, M.D. Anderson Cancer Center
!             1515 Holcombe Boulevard
!             Houston, TX      77030
!
! This work was supported by grant CA-16672 from the National Cancer Institute.

!                    GENerate POIsson random deviate

!                            Function

! Generates a single random deviate from a Poisson distribution with mean mu.

!                            Arguments

!     mu --> The mean of the Poisson distribution from which
!            a random deviate is to be generated.
!                              REAL(KIND=JPRB)mu

!                              Method

!     For details see:

!               Ahrens, J.H. and Dieter, U.
!               Computer Generation of Poisson Deviates
!               From Modified Normal Distributions.
!               ACM Trans. Math. Software, 8, 2
!               (June 1982),163-179

!     TABLES: COEFFICIENTS A0-A7 FOR STEP F. FACTORIALS FACT
!     COEFFICIENTS A(K) - FOR PX = FK*V*V*SUM(A(K)*V**K)-DEL

!     SEPARATION OF CASES A AND B

!     .. Scalar Arguments ..
REAL(KIND=JPRB)   , INTENT(IN) :: mu
INTEGER(KIND=JPIM), INTENT(IN) :: streammax
REAL(KIND=JPRB)   , INTENT(IN) :: rnd(streammax)
INTEGER(KIND=JPIM), INTENT(INOUT):: idx
INTEGER             :: ival
!     ..
!     .. Local Scalars ..
REAL(KIND=JPRB)          :: b1, b2, c, c0, c1, c2, c3, del, difmuk, e, fk, fx, fy, g,  &
                            omega, px, py, t, u, v, x, xx
REAL(KIND=JPRB)   :: s, d, p, q, p0
INTEGER :: j, k, kflag
LOGICAL :: full_init
INTEGER:: l, m
!INTEGER:: idx
!     ..
!     .. Local Arrays ..
REAL(KIND=JPRB)    :: pp(35)
!     ..
!     .. Data statements ..
REAL(KIND=JPRB), PARAMETER :: a0 = -.5_JPRB, a1 = .3333333_JPRB, a2 = -.2500068_JPRB, a3 = .2000118_JPRB,  &
                   a4 = -.1661269_JPRB, a5 = .1421878_JPRB, a6 = -.1384794_JPRB,   &
                   a7 = .1250060_JPRB

REAL(KIND=JPRB), PARAMETER :: fact(10) = (/ 1._JPRB, 1._JPRB, 2._JPRB, 6._JPRB, 24._JPRB, 120._JPRB, 720._JPRB, 5040._JPRB,  &
                                 40320._JPRB, 362880._JPRB /)
REAL(KIND=JPRB)    :: nn

!  PRINT *,'Poisson process ..............'
!     ..
!     .. Executable Statements ..
! Random numbers to use in generation of new clouds

IF (mu > 10.0_JPRB) THEN
!     C A S E  A. (RECALCULATION OF S, D, L IF MU HAS CHANGED)

  s = SQRT(mu)
  d = 6.0_JPRB*mu*mu

!             THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
!             PROBABILITIES FK WHENEVER K >= M(MU). L=IFIX(MU-1.1484)
!             IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .
  
  l = INT(mu - 1.1484_JPRB)
  full_init = .FALSE.


!     STEP N. NORMAL SAMPLE - random_normal() FOR STANDARD NORMAL DEVIATE

  CALL random_normal(streammax,rnd,idx,nn)
  g = mu + s*nn
!  g = mu + s*random_normal(rnd,idx)
  IF (g > 0.0_JPRB) THEN
    ival = INT(g)

!     STEP I. IMMEDIATE ACCEPTANCE IF ival IS LARGE ENOUGH

    IF (ival>=l) RETURN

!     STEP S. SQUEEZE ACCEPTANCE - SAMPLE U

    fk = ival
    difmuk = mu - fk
!    CALL RANDOM_NUMBER(u)
    u=rnd(idx)
    if (u .eq. 0.) then
       idx=idx+1
       u=rnd(idx)
    endif
    idx=idx+1
!    if (idx.ge.5000) write(7,*) 'random_Poisson, ran out of idx, step N',idx
!    flush(7)
    IF (d*u >= difmuk*difmuk*difmuk) RETURN
  END IF

!     STEP P. PREPARATIONS FOR STEPS Q AND H.
!             (RECALCULATIONS OF PARAMETERS IF NECESSARY)
!             .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
!             THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
!             APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
!             C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.

  IF (.NOT. full_init) THEN
    omega = .3989423_JPRB/s
    b1 = .4166667E-1_JPRB/mu
    b2 = .3_JPRB*b1*b1
    c3 = .1428571_JPRB*b1*b2
    c2 = b2 - 15._JPRB*c3
    c1 = b1 - 6._JPRB*b2 + 45._JPRB*c3
    c0 = 1._JPRB - b1 + 3._JPRB*b2 - 15._JPRB*c3
    c = .1069_JPRB/mu
    full_init = .true.
  END IF

  IF (g < 0.0_JPRB) GO TO 50

!             'SUBROUTINE' F IS CALLED (KFLAG=0 FOR CORRECT RETURN)

  kflag = 0
  GO TO 70

!     STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)

  40 IF (fy-u*fy <= py*EXP(px-fx)) RETURN

!     STEP E. EXPONENTIAL SAMPLE - random_exponential() FOR STANDARD EXPONENTIAL
!             DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
!             (IF T <= -.6744 THEN PK < FK FOR ALL MU >= 10.)

!  50 e = random_exponential(rnd,idx)
  50 CALL random_exponential(rnd,idx,streammax,e)
!  CALL RANDOM_NUMBER(u)
  u=rnd(idx)
  if (u .eq. 0.) then
      idx=idx+1
      u=rnd(idx)
  endif
  idx=idx+1
!  if (idx.ge.5000) write(7,*) 'random_Poisson, ran out of idx, step E',idx
!  flush(7)
  u = u + u - one
  t = 1.8_JPRB + SIGN(e, u)
  IF (t <= (-.6744_JPRB)) GO TO 50
  ival = INT(mu + s*t)
  fk = ival
  difmuk = mu - fk

!             'SUBROUTINE' F IS CALLED (KFLAG=1 FOR CORRECT RETURN)

  kflag = 1
  GO TO 70

!     STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)

  60 IF (c*ABS(u) > py*EXP(px+e) - fy*EXP(fx+e)) GO TO 50
  RETURN

!     STEP F. 'SUBROUTINE' F. CALCULATION OF PX, PY, FX, FY.
!             CASE ival < 10 USES FACTORIALS FROM TABLE FACT

  70 IF (ival>=10) GO TO 80
  px = -mu
  py = mu**ival/fact(ival+1)
  GO TO 110

!             CASE ival >= 10 USES POLYNOMIAL APPROXIMATION
!             A0-A7 FOR ACCURACY WHEN ADVISABLE
!             .8333333E-1=1./12.  .3989423=(2*PI)**(-.5)

  80 del = .8333333E-1_JPRB/fk
  del = del - 4.8_JPRB*del*del*del
  v = difmuk/fk
  IF (ABS(v)>0.25_JPRB) THEN
    px = fk*LOG(one + v) - difmuk - del
  ELSE
    px = fk*v*v* (((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0) - del
  END IF
  py = .3989423_JPRB/SQRT(fk)
  110 x = (half - difmuk)/s
  xx = x*x
  fx = -half*xx
  fy = omega* (((c3*xx + c2)*xx + c1)*xx + c0)
  IF (kflag <= 0) GO TO 40
  GO TO 60

!---------------------------------------------------------------------------
!     C A S E  B.    mu < 10
!     START NEW TABLE AND CALCULATE P0 IF NECESSARY

ELSE

  m = MAX(1, INT(mu))
  l = 0
  p = EXP(-mu)
  q = p
  p0 = p

!     STEP U. UNIFORM SAMPLE FOR INVERSION METHOD

  DO

!    CALL RANDOM_NUMBER(u)
    u=rnd(idx)
    if (u .eq. 0.) then
      idx=idx+1
      u=rnd(idx)
    endif
    idx=idx+1
!    if (idx.ge.5000) write(7,*) 'random_Poisson, ran out of idx, step U',idx
!    flush(7)
    ival = 0
    IF (u <= p0) RETURN
!  PRINT *,'STEP U...' , u , p0

!     STEP T. TABLE COMPARISON UNTIL THE END PP(L) OF THE
!             PP-TABLE OF CUMULATIVE POISSON PROBABILITIES
!             (0.458=PP(9) FOR MU=10)

    IF (l == 0) GO TO 150
    j = 1
    IF (u > 0.458_JPRB) j = MIN(l, m)
    DO k = j, l
      IF (u <= pp(k)) GO TO 180
    END DO
    IF (l == 35) CYCLE

!     STEP C. CREATION OF NEW POISSON PROBABILITIES P
!             AND THEIR CUMULATIVES Q=PP(K)

    150 l = l + 1
    DO k = l, 35
      p = p*mu / k
      q = q + p
      pp(k) = q
      IF (u <= q) GO TO 170
    END DO
    l = 35
  END DO

  170 l = k
  180 ival = k

  RETURN
END IF
RETURN
END FUNCTION random_Poisson

!SUBROUTINE seed_random_number(iounit)

!INTEGER, INTENT(IN)  :: iounit

! Local variables

!INTEGER              :: k
!INTEGER, ALLOCATABLE :: seed(:)

!CALL RANDOM_SEED(SIZE=k)
!ALLOCATE( seed(k) )

!WRITE(*, '(a, i2, a)')' Enter ', k, ' integers for random no. seeds: '
!READ(*, *) seed
!WRITE(iounit, '(a, (7i10))') ' Random no. seeds: ', seed
!CALL RANDOM_SEED(PUT=seed)

!DEALLOCATE( seed )

!RETURN
!END SUBROUTINE seed_random_number


END MODULE random_rewrite
