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

! Contains the source function routines for
! - wind input S_in
! - nonlinear energy transfer S_nl
! - dissipation S_dis (bottom friction and whitecapping)
!
! as well as the implicit time integration scheme for the
! total source function S = S_in + S_dis + S_nl

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_wave_source

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_parallel_config,     ONLY: nproma
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, min_rlcell
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_run_config,          ONLY: dtime
  USE mo_physical_constants,  ONLY: grav
  USE mo_math_constants,      ONLY: deg2rad, pi2
  USE mo_wave_types,          ONLY: t_wave_source, t_wave_diag
  USE mo_wave_config,         ONLY: t_wave_config
  USE mo_wave_constants,      ONLY: DELTA, CONSS

  IMPLICIT NONE

  PRIVATE


  ! source functions
  PUBLIC :: src_wind_input
  PUBLIC :: src_dissipation
  PUBLIC :: src_bottom_friction
  PUBLIC :: src_nonlinear_transfer
  PUBLIC :: src_wave_breaking
  !
  ! implicit time integration scheme
  PUBLIC :: integrate_in_time_src


  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_source'

CONTAINS

  !>
  !! Implicit scheme for integration of source functions
  !!
  !! S.D.HASSELMANN.  MPI
  !! H. GUENTHER AND L. ZAMBRESKY  OPTIMIZATION PERFORMED.
  !! H. GUENTHER      GKSS/ECMWF   OCTOBER  1989  NEW WIND FIELD
  !!                                              INTERFACE AND
  !!                                              TIME COUNTING
  !! P.A.E.M. JANSSEN KNMI         AUGUST   1990  COUPLED MODEL
  !! H. GUENTHER      GKSS/ECMWF   JUNE     1991  NEW SEPARATION OF
  !!                                              DIAG- AND PROGNOSTIC
  !!                                              PART OF SPECTRUM.
  !! P.A.E.M. JANSSEN ECMWF        FEBRUARY 1995  ADD MINIMUM VALUE (FMIN).
  !! H. GUENTHER      GKSS         JANUARY  2002  FT90
  !! ROOP LALBEHARRY  MSC/ARMN     APRIL    2003  PHILLIPS SOURCE
  !! ERIK MYKLEBUST                FEBRUARY 2005  OPTIMIZATION
  !! H. GUENTHER      GKSS         JANUARY  2014  Wave breaking
  !!
  !! PURPOSE:
  !! --------
  !! THE IMPLICIT SCHEME ENABLES THE USE OF A TIMESTEP WHICH IS
  !! LARGE COMPARED WITH THE CHARACTERISTIC DYNAMIC TIME SCALE.
  !! THE SCHEME IS REQUIRED FOR THE HIGH FREQUENCIES WHICH
  !! RAPIDLY ADJUST TO A QUASI-EQUILIBRIUM.
  !!
  !! METHOD:
  !! -------
  !! THE SPECTRUM AT TIME (TN+1) IS COMPUTED AS
  !! FN+1=FN+DELT*(\alpha*SN + (1-\alpha)*SN+1),
  !! WHERE SN IS THE TOTAL SOURCE FUNCTION AT TIME TN, SN+1=SN+(DS/DF)*DF
  !! - ONLY THE DIAGONAL TERMS OF THE FUNCTIONAL MATRIX DS/DF ARE COMPUTED,
  !!   THE NONDIAGONAL TERMS ARE NEGLIGIBLE.
  !! INCREASE OF SPECTRUM IN A TIME STEP IS LIMITED TO A FINITE
  !! FRACTION OF A TYPICAL F**(-4) EQUILIBRIUM SPECTRUM.
  !!
  !! REFERENCE:
  !! ----------
  !! The WAMDI Group, The WAM Model - A Third Generation Wave Prediction Model,
  !! Journal of Physical Oceanography, 1988
  !!
  !! Adaptation of WAM 4.5 code.
  !!
  SUBROUTINE integrate_in_time_src(p_patch, wave_config, p_diag, p_source, dir10m, tracer)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = 'integrate_in_time_src'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    TYPE(t_wave_diag),           INTENT(IN)    :: p_diag
    TYPE(t_wave_source),         INTENT(IN)    :: p_source
    REAL(wp),                    INTENT(IN)    :: dir10m(:,:)
    REAL(wp),                    INTENT(INOUT) :: tracer(:,:,:,:)

    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jc,jf,jd,jt,jk

    REAL(wp) :: temp_1, temp_2, temp_3
    REAL(wp) :: sprd
    REAL(wp) :: delfl

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev

    ! save some paperwork
    wc => wave_config

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jf,jd,jt,i_startidx,i_endidx,delfl,sprd,temp_1,temp_2,temp_3) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1,wc%nfreqs
        DO jd = 1,wc%ndirs
          !
          jt = wc%tracer_ind(jd,jf)
          !
          DO jc = i_startidx, i_endidx
            delfl = 5.0E-07_wp * grav / wc%freqs(jf)**4 * dtime
            temp_2 = p_diag%ustar(jc,jb) * delfl &
                 &  * MAX(p_diag%femeanws(jc,jb),p_diag%femean(jc,jb))
            sprd = MAX(0._wp, COS(wc%dirs(jd)-dir10m(jc,jb)*deg2rad) )**2

            temp_1 = dtime * p_source%sl(jc,jb,jt) &
                 / MAX(1._wp, 1._wp -  dtime * wc%impl_fac * p_source%fl(jc,jb,jt))

            temp_3 = MIN(ABS(temp_1),temp_2)

            tracer(jc,jk,jb,jt) = tracer(jc,jk,jb,jt)  + SIGN(temp_3,temp_1)
            tracer(jc,jk,jb,jt) = MAX(tracer(jc,jk,jb,jt), p_diag%FLMINFR(jc,jb,jf)*sprd)
          END DO
        END DO
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE integrate_in_time_src


  !>
  !! Calculation of input source function
  !!
  !! Compute input source function and the functional derivative of input
  !! source function.
  !!
  !! Adaptation of WAM 4.5 code.
  !! SINPUT
  !!
  !! Reference
  !! P. Janssen, JPO, 1989.
  !! P. Janssen, JPO., 1991.
  !!
  SUBROUTINE src_wind_input(p_patch, wave_config, dir10m, tracer, p_diag, p_source)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
         & routine =  modname//'src_wind_input'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: dir10m(:,:)
    REAL(wp),                    INTENT(IN)    :: tracer(:,:,:,:)
    TYPE(t_wave_diag),           INTENT(IN)    :: p_diag
    TYPE(t_wave_source),         INTENT(INOUT) :: p_source

    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jc,jf,jd,jt,jk

    REAL(wp) :: fac, const, const3, xk, cm, ucn, zcn, sh, cnsn
    REAL(wp) :: xv1d, temp, zbeta1, x, zlog, zlog2x, ufac

    wc => wave_config

    const3   = 2.0_wp * wc%xkappa / wc%betamax * wc%xkappa**2

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jf,jd,jt,i_startidx,i_endidx,fac,const,xk,cm,  &
!$OMP            ucn,zcn,sh,cnsn,xv1d,temp,zbeta1,x,zlog,zlog2x,ufac) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      FRE:DO jf = 1,wc%nfreqs
        fac = pi2 * wc%freqs(jf)
        const = fac * wc%xeps * wc%betamax / (wc%xkappa*wc%xkappa)

        DIR:DO jd = 1,wc%ndirs
          jt = wc%tracer_ind(jd,jf)

          DO jc = i_startidx, i_endidx
            xk = p_diag%wave_num_c(jc,jb,jf)
            cm = xk / fac
            ucn = p_diag%ustar(jc,jb) * cm + wc%zalp
            zcn = LOG(xk * p_diag%z0(jc,jb))
            sh = fac*fac / (grav * xk)
            cnsn = const * sh
            xv1d = -1.0_wp / (p_diag%ustar(jc,jb) / wc%xkappa * zcn * cm)
            temp = COS(wc%dirs(jd) - dir10m(jc,jb)*deg2rad)
            zbeta1 = const3 * (temp - xv1d) * ucn*ucn

            IF (temp > 0.01_wp) THEN
              x = temp * ucn
              zlog = zcn + wc%xkappa / x
              IF (zlog < 0._wp) THEN
                zlog2x = zlog*zlog * x
                ufac = EXP(zlog) * zlog2x*zlog2x + zbeta1
                p_source%llws(jc,jb,jt) = 1
              ELSE
                ufac = zbeta1
                p_source%llws(jc,jb,jt) = 0
              END IF
            ELSE
              ufac = zbeta1
              p_source%llws(jc,jb,jt) = 0
            END IF

            p_source%fl(jc,jb,jt) = cnsn * ufac
            p_source%sl(jc,jb,jt) = tracer(jc,jk,jb,jt) * p_source%fl(jc,jb,jt) !SL
          END DO
        END DO DIR
      END DO FRE
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE src_wind_input


  !>
  !! Calculation of dissipation source function
  !!
  !! Compute dissipation source function and store additively into
  !! net source function array. Also compute functional derivative
  !! of dissipation source function.
  !!
  !! Adaptation of WAM 4.5 code SDISSIP
  !!     S.D.HASSELMANN.
  !!     MODIFIED TO SHALLOW WATER : G. KOMEN , P. JANSSEN
  !!     OPTIMIZATION : L. ZAMBRESKY
  !!     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
  !!     H. GUENTHER GKSS  FEBRUARY 2002       FT 90
  !!     J. BIDLOT   ECMWF  NOVEMBER 2004  REFORMULATION BASED ON AKMEAN
  !!                                       AND FMEAN.
  !! Reference
  !! G.Komen, S. Hasselmann And K. Hasselmann, On The Existence
  !!          Of A Fully Developed Windsea Spectrum, JGR, 1984.
  !!
  SUBROUTINE src_dissipation(p_patch, wave_config, wave_num_c, tracer, p_diag, p_source)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
         & routine =  modname//'src_dissipation'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: wave_num_c(:,:,:) !< wave number (1/m)
    REAL(wp),                    INTENT(IN)    :: tracer(:,:,:,:)
    TYPE(t_wave_diag),           INTENT(IN)    :: p_diag
    TYPE(t_wave_source),         INTENT(INOUT) :: p_source

    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jk,jc,jf,jd,jt

    REAL(wp) :: sds, temp, sdiss

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev

    ! save some paperwork
    wc => wave_config

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jd,jt,jc,sds,temp,sdiss,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jf = 1,wc%nfreqs
        DO jd = 1, wc%ndirs
          jt = wc%tracer_ind(jd,jf)
          DO jc = i_startidx, i_endidx

            sds = CONSS * p_diag%f1mean(jc,jb) * p_diag%emean(jc,jb)**2 * p_diag%xkmean(jc,jb)**4
            temp = wave_num_c(jc,jb,jf) / p_diag%xkmean(jc,jb)
            temp = sds * ((1.0_wp - DELTA) * temp +  DELTA * temp**2)
            sdiss = temp * tracer(jc,jk,jb,jt)

            p_source%sl(jc,jb,jt) = p_source%sl(jc,jb,jt) + sdiss
            p_source%fl(jc,jb,jt) = p_source%fl(jc,jb,jt) + temp
          END DO
        END DO
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE src_dissipation

  !>
  !! Calculation of dissipation due to depth-induced wave breaking
  !!
  !! Adaptation of WAM 4.6 code, SUBROUTINE SFBRK
  !!
  !! output:
  !! hrms_frac
  !! wbr_frac
  !! sl
  !! fl
  !!
  !! Reference
  !! Battjes & Janssen (Coastal Engineering, 1978)
  !!
  SUBROUTINE src_wave_breaking(p_patch, wave_config, depth_c, tracer, p_diag, p_source)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
         & routine =  modname//'src_wave_breaking'

    TYPE(t_patch),       INTENT(IN)         :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN) :: wave_config
    REAL(wp),            INTENT(IN)         :: depth_c(:,:)
    REAL(wp),            INTENT(IN)         :: tracer(:,:,:,:)
    TYPE(t_wave_diag),   INTENT(INOUT)      :: p_diag
    TYPE(t_wave_source), INTENT(INOUT)      :: p_source

    TYPE(t_wave_config), POINTER :: wc => NULL()

    REAL(wp), PARAMETER :: alpha = 1.0_wp

    REAL(wp) :: qb, sbr(nproma), dsbr(nproma)

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jc,jf,jd,jt,jk

    wc => wave_config

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev

    CALL breaking_waves_frac(p_patch, p_diag%emean, depth_c, & ! IN
      &                      p_diag%hrms_frac, & !OUT
      &                      p_diag%wbr_frac)    !OUT

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jf,jd,jt,i_startidx,i_endidx,qb,sbr,dsbr) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx

        qb = MIN(1.0_wp,p_diag%wbr_frac(jc,jb))

        sbr(jc) = -alpha*2.0_wp * p_diag%f1mean(jc,jb)

        IF (p_diag%hrms_frac(jc,jb) <= 1.0_wp) THEN
          sbr(jc) = sbr(jc)*qb/p_diag%hrms_frac(jc,jb)
        END IF

        IF ( (p_diag%hrms_frac(jc,jb) < 1.0_wp             ) .AND. &
          &  (ABS(p_diag%hrms_frac(jc,jb)-qb) > 0.0_wp ) ) THEN
          dsbr(jc) = sbr(jc) * (1.0_wp - qb) / (p_diag%hrms_frac(jc,jb) - qb)
        ELSE
          dsbr(jc) = 0.0_wp
        END IF
      END DO

      DO jf = 1,wc%nfreqs
        DO jd = 1,wc%ndirs
          !
          jt = wc%tracer_ind(jd,jf)
          !
          DO jc = i_startidx, i_endidx
            p_source%sl(jc,jb,jt) = p_source%sl(jc,jb,jt) + sbr(jc) * tracer(jc,jk,jb,jt)
            p_source%fl(jc,jb,jt) = p_source%fl(jc,jb,jt) + dsbr(jc)
          END DO
        END DO
      END DO

    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE src_wave_breaking

  !>
  !! Calculation of fraction of breaking waves
  !!
  !! Adaptation of WAM 4.6 code, SUBROUTINE CMPQB
  !! ADDED BY WEIMIN LUO, POL, MAY 1996
  !! BASED ON THE CODE OF G. Ph. van Vledder, Delft Hydraulics
  !!
  !! Method
  !!
  !!     Newton-Raphson implementation of
  !!
  !!     1 - QB
  !!     ------ = - (HRMS/HMAX)^2
  !!     ln(QB)
  !!
  !!
  !!     If HRMS > HMAX --> QB = 1
  !!
  SUBROUTINE breaking_waves_frac(p_patch, emean, depth_c, hrms_frac, wbr_frac)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      & routine =  modname//'breaking_waves_frac'

    TYPE(t_patch), INTENT(IN)    :: p_patch
    REAL(wp),      INTENT(IN)    :: emean(:,:)     ! total energy (nproma,nblks_c)
    REAL(wp),      INTENT(IN)    :: depth_c(:,:)
    REAL(wp),      INTENT(INOUT) :: hrms_frac(:,:) ! square ratio (Hrms / Hmax)**2 BB
    REAL(wp),      INTENT(INOUT) :: wbr_frac(:,:)  ! fraction of breaking waves QB

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb

    REAL(wp), PARAMETER :: gamd  = 0.8  !! Parameter of depth limited wave height
    REAL(wp) :: frac_0(nproma) !Q0

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,frac_0) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      ! calculation of BB (Hrms / Hmax)**2 BB = 8.*EMEAN/(GAMD*DEPTH)**2 
      DO jc = i_startidx, i_endidx
        hrms_frac(jc,jb) = 8.0_wp * emean(jc,jb)/(gamd*depth_c(jc,jb))**2
      END DO

      ! initialisation of frac_0
      DO jc = i_startidx, i_endidx
        IF (hrms_frac(jc,jb)>=0.25_wp) THEN
          frac_0(jc) = (2.0_wp * SQRT(hrms_frac(jc,jb))-1.0_wp )**2
        ELSE
          frac_0(jc) = 0.0_wp
        END IF
      END DO

      DO jc = i_startidx, i_endidx
        IF (hrms_frac(jc,jb) < 1.0_wp) THEN
          wbr_frac(jc,jb) = frac_0(jc) - &
            & hrms_frac(jc,jb) * (frac_0(jc)-EXP((frac_0(jc)-1.0_wp)/hrms_frac(jc,jb))) / &
            & (hrms_frac(jc,jb) - EXP((frac_0(jc)-1.0_wp)/hrms_frac(jc,jb)))
        ELSE
          wbr_frac(jc,jb) = 1.0_wp
        END IF
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE breaking_waves_frac


  !>
  !! Calculation of bottom friction source function
  !!
  !! Compute dissipation of wave energy due to bottom friction
  !!
  !! Adaptation of WAM 4.5 code
  !!     G.J.KOMEN AND Q.D.GAO
  !!     OPTIMIZED BY L.F. ZAMBRESKY
  !!     H. GUENTHER   GKSS  FEBRUARY 2002       FT 90
  !!     E. MYKLEBUST        FEBRUARY 2005       OPTIMIZATION
  !!
  !! Reference
  !!  HASSELMANN ET AL, D. HYDR. Z SUPPL A12(1973) (JONSWAP)
  !!  BOUWS AND KOMEN, JPO 13(1983)1653-1658
  !!
  SUBROUTINE src_bottom_friction(p_patch, wave_config, wave_num_c, depth, tracer, p_source)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
         & routine =  modname//'src_bottom_friction'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: wave_num_c(:,:,:) !< wave number (1/m)
    REAL(wp),                    INTENT(IN)    :: depth(:,:)
    REAL(wp),                    INTENT(IN)    :: tracer(:,:,:,:)
    TYPE(t_wave_source),         INTENT(INOUT) :: p_source

    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jk,jc,jf,jd,jt

    REAL(wp) :: const, sbo

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev

    ! save some paperwork
    wc => wave_config

    const = -2.0_wp*0.038_wp/grav

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jd,jt,jc,sbo,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jf = 1,wc%nfreqs
        DO jd = 1, wc%ndirs
          jt = wc%tracer_ind(jd,jf)
          DO jc = i_startidx, i_endidx
            sbo = MIN(2.0_wp * depth(jc,jb) * wave_num_c(jc,jb,jf),50.0_wp)
            sbo = const * wave_num_c(jc,jb,jf) / SINH(sbo)
            p_source%sl(jc,jb,jt) = p_source%sl(jc,jb,jt) + sbo*tracer(jc,jk,jb,jt)
            p_source%fl(jc,jb,jt) = p_source%fl(jc,jb,jt) + sbo
          END DO
        END DO
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE src_bottom_friction


  !>
  !!
  !! Computation of nonlinear transfer rate and its
  !! functional derivative (diagonal terms only) and
  !! addition to corresponding net expressions.
  !!
  !! Adaptation of WAM 4.5 SNONLIN
  !! S.D. Hasselmann.  MPI
  !! G. Komen, P. Janssen   KNMI        modified to shallow water
  !! H. Guenther, L. Zambresky          optimized
  !! H. Guenther  GKSS/ECMWF  June 1991 interactions between diag-
  !!                                    and prognostic part.
  !! H. Guenther  GKSS  February 2002   FT 90
  !! E. Myklebust       February 2005   optimization
  !!
  SUBROUTINE src_nonlinear_transfer(p_patch, wave_config, depth, tracer, p_diag, p_source)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
         & routine =  modname//'src_nonlinear_transfer'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: depth(:,:)
    REAL(wp),                    INTENT(IN)    :: tracer(:,:,:,:)
    TYPE(t_wave_diag),           INTENT(IN)    :: p_diag
    TYPE(t_wave_source),         INTENT(INOUT) :: p_source

    ! local
    TYPE(t_wave_config), POINTER :: wc             => NULL()
    INTEGER,             POINTER:: tr_idx(:,:,:,:) => NULL()

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jk,jb,jf

    INTEGER :: nfreqs, ndirs
    INTEGER :: MP, MP1, MM, MM1, IC, IP, IP1, IM, IM1, KH, K

    REAL(wp) :: FFACP, FFACP1, FFACM1, FTAIL, FKLAMP, FKLAMP1, GW1, GW2, GW3, GW4
    REAL(wp) :: FKLAMPA, FKLAMPB, FKLAMP2, FKLAPA2, FKLAPB2, FKLAP12, FKLAP22
    REAL(wp) :: FKLAMM, FKLAMM1, GW5, GW6, GW7, GW8, FKLAMMA, FKLAMMB, FKLAMM2
    REAL(wp) :: FKLAMA2, FKLAMB2, FKLAM12, FKLAM22
    REAL(wp) :: SAP, SAM, FIJ, FAD1, FAD2, FCEN

    REAL(wp) :: AD
    REAL(wp) :: DELAD, DELAP, DELAM
    REAL(wp) :: FTEMP, ENHFR


    ! convenience pointers
    wc     => wave_config
    tr_idx => p_diag%non_lin_tr_ind(:,:,:,:)
    nfreqs = wc%nfreqs
    ndirs  = wc%ndirs

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jf,i_startidx,i_endidx,ENHFR,FTEMP,                              &
!$OMP           MP,MP1,MM,MM1,IC,IP,IP1,IM,IM1,KH,K,FFACP,                               &
!$OMP           FFACP1,FFACM1,FTAIL,FKLAMP,FKLAMP1, GW1, GW2, GW3, GW4,FKLAMPA, FKLAMPB, &
!$OMP           FKLAMP2, FKLAPA2, FKLAPB2, FKLAP12, FKLAP22, FKLAMM, FKLAMM1, GW5, GW6, &
!$OMP           GW7, GW8, FKLAMMA, FKLAMMB, FKLAMM2,FKLAMA2, FKLAMB2, FKLAM12, FKLAM22, &
!$OMP           SAP, SAM, FIJ, FAD1, FAD2, FCEN, AD, DELAD, DELAP, DELAM ) ICON_OMP_DEFAULT_SCHEDULE
    ljb: DO jb = i_startblk, i_endblk
       CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
            &                 i_startidx, i_endidx, i_rlstart, i_rlend)
       DO jc = i_startidx, i_endidx
          ENHFR = MAX(0.75_wp*depth(jc,jb)*p_diag%AKMEAN(jc,jb), 0.5_wp)
          ENHFR = 1.0_wp + (5.5_wp/ENHFR) * (1.0_wp-0.833_wp*ENHFR) &
               &                * EXP(-1.25_wp*ENHFR)
          p_diag%ENH(jc,jb) = ENHFR
       END DO

      FRE4: DO jf = 1,nfreqs+4
        MP  = p_diag%IKP (jf)
        MP1 = p_diag%IKP1(jf)
        MM  = p_diag%IKM (jf)
        MM1 = p_diag%IKM1(jf)
        FFACP  = 1._wp
        FFACP1 = 1._wp
        FFACM1 = 1._wp
        FTAIL  = 1._wp
        IC  = jf
        IP  = MP
        IP1 = MP1
        IM  = MM
        IM1 = MM1
        IF (IP1.GT.nfreqs) THEN
          FFACP1 = wc%FRH(IP1-nfreqs+1)
          IP1 = nfreqs
          IF (IP.GT.nfreqs) THEN
            FFACP  = wc%FRH(IP-nfreqs+1)
            IP  = nfreqs
            IF (IC.GT.nfreqs) THEN
              FTAIL  = wc%FRH(IC-nfreqs+1)
              IC  = nfreqs
              IF (IM1.GT.nfreqs) THEN
                FFACM1 = wc%FRH(IM1-nfreqs+1)
                IM1 = nfreqs
              END IF
            END IF
          END IF
        END IF

        FKLAMP  = p_diag%FKLAP(jf)
        FKLAMP1 = p_diag%FKLAP1(jf)
        GW2 = FKLAMP1*FFACP*wc%DAL1
        GW1 = GW2*wc%CL11
        GW2 = GW2*wc%ACL1
        GW4 = FKLAMP*FFACP1*wc%DAL1
        GW3 = GW4*wc%CL11
        GW4 = GW4*wc%ACL1
        FKLAMPA = FKLAMP*wc%CL11
        FKLAMPB = FKLAMP*wc%ACL1
        FKLAMP2 = FKLAMP1*wc%ACL1
        FKLAMP1 = FKLAMP1*wc%CL11
        FKLAPA2 = FKLAMPA**2
        FKLAPB2 = FKLAMPB**2
        FKLAP12 = FKLAMP1**2
        FKLAP22 = FKLAMP2**2
        FKLAMM  = p_diag%FKLAM(jf)
        FKLAMM1 = p_diag%FKLAM1(jf)
        GW6 = FKLAMM1*wc%DAL2
        GW5 = GW6*wc%CL21
        GW6 = GW6*wc%ACL2
        GW8 = FKLAMM*FFACM1*wc%DAL2
        GW7 = GW8*wc%CL21
        GW8 = GW8*wc%ACL2
        FKLAMMA = FKLAMM*wc%CL21
        FKLAMMB = FKLAMM*wc%ACL2
        FKLAMM2 = FKLAMM1*wc%ACL2
        FKLAMM1 = FKLAMM1*wc%CL21
        FKLAMA2 = FKLAMMA**2
        FKLAMB2 = FKLAMMB**2
        FKLAM12 = FKLAMM1**2
        FKLAM22 = FKLAMM2**2

        ! shallow water case
        IF (jf.GT.4) THEN! UNTIL 7.
          IF (MM1.LE.nfreqs) THEN! UNTIL 6.
            IF (jf .LE.nfreqs) THEN! UNTIL 5.
              IF (MP .LE.nfreqs) THEN! UNTIL 4.
                IF (MP1.LE.nfreqs) THEN! UNTIL 3.
                  !     2.1.1   ANGULAR LOOP.                                     !
                  DIR2: DO K = 1,ndirs !DIR2

                    DO jc = i_startidx, i_endidx
                      !     2.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND     !
                      !             DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.         !
                      !             ----------------------------------------------    !

                      !     2.1 LOOP FOR ANLULAR SYMMETRY.                            !
                      MIR2: DO KH = 1,2

                         SAP = &
                              GW1*tracer(jc,jk,jb,tr_idx(1,jf,KH,K)) + &
                              GW2*tracer(jc,jk,jb,tr_idx(2,jf,KH,K)) + &
                              GW3*tracer(jc,jk,jb,tr_idx(3,jf,KH,K)) + &
                              GW4*tracer(jc,jk,jb,tr_idx(4,jf,KH,K))
                         SAM = &
                              GW5*tracer(jc,jk,jb,tr_idx(5,jf,KH,K)) + &
                              GW6*tracer(jc,jk,jb,tr_idx(6,jf,KH,K)) + &
                              GW7*tracer(jc,jk,jb,tr_idx(7,jf,KH,K)) + &
                              GW8*tracer(jc,jk,jb,tr_idx(8,jf,KH,K))

                         FTEMP = p_diag%AF11(jf) * p_diag%ENH(jc,jb)
                         FIJ = tracer(jc,jk,jb,tr_idx(9,jf,KH,K))*FTAIL
                         FAD1 = FIJ*(SAP+SAM)
                         FAD2 = FAD1-2._wp*SAP*SAM
                         FAD1 = FAD1+FAD2
                         FCEN = FTEMP*FIJ
                         AD = FAD2*FCEN
                         DELAD = FAD1*FTEMP
                         DELAP = (FIJ-2._wp*SAM)*wc%DAL1*FCEN
                         DELAM = (FIJ-2._wp*SAP)*wc%DAL2*FCEN

                         p_source%sl(jc,jb,tr_idx(10,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(10,jf,KH,K)) + AD*FKLAMM1 !SL
                         p_source%sl(jc,jb,tr_idx(11,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(11,jf,KH,K)) + AD*FKLAMM2 !SL
                         p_source%fl(jc,jb,tr_idx(10,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(10,jf,KH,K)) + DELAM*FKLAM12 !FL
                         p_source%fl(jc,jb,tr_idx(11,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(11,jf,KH,K)) + DELAM*FKLAM22 !FL

                         p_source%sl(jc,jb,tr_idx(12,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(12,jf,KH,K)) + AD*FKLAMMA
                         p_source%sl(jc,jb,tr_idx(13,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(13,jf,KH,K)) + AD*FKLAMMB
                         p_source%fl(jc,jb,tr_idx(12,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(12,jf,KH,K)) + DELAM*FKLAMA2
                         p_source%fl(jc,jb,tr_idx(13,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(13,jf,KH,K)) + DELAM*FKLAMB2

                         p_source%sl(jc,jb,tr_idx(14,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(14,jf,KH,K)) - 2._wp*AD
                         p_source%fl(jc,jb,tr_idx(14,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(14,jf,KH,K)) - 2._wp*DELAD

                         p_source%sl(jc,jb,tr_idx(15,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(15,jf,KH,K)) + AD*FKLAMP1
                         p_source%sl(jc,jb,tr_idx(16,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(16,jf,KH,K)) + AD*FKLAMP2
                         p_source%fl(jc,jb,tr_idx(15,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(15,jf,KH,K)) + DELAP*FKLAP12
                         p_source%fl(jc,jb,tr_idx(16,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(16,jf,KH,K)) + DELAP*FKLAP22

                         p_source%sl(jc,jb,tr_idx(17,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(17,jf,KH,K)) + AD*FKLAMPA
                         p_source%sl(jc,jb,tr_idx(18,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(18,jf,KH,K)) + AD*FKLAMPB
                         p_source%fl(jc,jb,tr_idx(17,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(17,jf,KH,K)) + DELAP*FKLAPA2
                         p_source%fl(jc,jb,tr_idx(18,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(18,jf,KH,K)) + DELAP*FKLAPB2
                      END DO MIR2
                    END DO  ! jc
                  END DO DIR2
                ELSE!IF (MP1.LE.ML) THEN
                  !     3.1.1   ANGULAR LOOP.                                   !
                  DIR3: DO K = 1, ndirs

                    DO jc = i_startidx, i_endidx
                      !     3.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND   !
                      !             DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.       !
                      !             ----------------------------------------------  !

                      !     3.1 LOOP FOR ANGULAR SYMMETRY.                          !
                      MIR3: DO KH = 1,2

                        SAP = &
                             GW1*tracer(jc,jk,jb,tr_idx(1,jf,KH,K)) + &
                             GW2*tracer(jc,jk,jb,tr_idx(2,jf,KH,K)) + &
                             GW3*tracer(jc,jk,jb,tr_idx(3,jf,KH,K)) + &
                             GW4*tracer(jc,jk,jb,tr_idx(4,jf,KH,K))
                        SAM = &
                             GW5*tracer(jc,jk,jb,tr_idx(5,jf,KH,K)) + &
                             GW6*tracer(jc,jk,jb,tr_idx(6,jf,KH,K)) + &
                             GW7*tracer(jc,jk,jb,tr_idx(7,jf,KH,K)) + &
                             GW8*tracer(jc,jk,jb,tr_idx(8,jf,KH,K))

                        FTEMP = p_diag%AF11(jf) * p_diag%ENH(jc,jb)
                        FIJ = tracer(jc,jk,jb,tr_idx(9,jf,KH,K))*FTAIL
                        FAD1 = FIJ*(SAP+SAM)
                        FAD2 = FAD1-2._wp*SAP*SAM
                        FAD1 = FAD1+FAD2
                        FCEN = FTEMP*FIJ
                        AD = FAD2*FCEN
                        DELAD = FAD1*FTEMP
                        DELAP = (FIJ-2._wp*SAM)*wc%DAL1*FCEN
                        DELAM = (FIJ-2._wp*SAP)*wc%DAL2*FCEN

                        p_source%sl(jc,jb,tr_idx(10,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(10,jf,KH,K)) + AD*FKLAMM1
                        p_source%sl(jc,jb,tr_idx(11,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(11,jf,KH,K)) + AD*FKLAMM2
                        p_source%fl(jc,jb,tr_idx(10,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(10,jf,KH,K)) + DELAM*FKLAM12
                        p_source%fl(jc,jb,tr_idx(11,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(11,jf,KH,K)) + DELAM*FKLAM22

                        p_source%sl(jc,jb,tr_idx(12,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(12,jf,KH,K)) + AD*FKLAMMA
                        p_source%sl(jc,jb,tr_idx(13,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(13,jf,KH,K)) + AD*FKLAMMB
                        p_source%fl(jc,jb,tr_idx(12,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(12,jf,KH,K)) + DELAM*FKLAMA2
                        p_source%fl(jc,jb,tr_idx(13,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(13,jf,KH,K)) + DELAM*FKLAMB2

                        p_source%sl(jc,jb,tr_idx(14,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(14,jf,KH,K)) - 2._wp*AD
                        p_source%fl(jc,jb,tr_idx(14,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(14,jf,KH,K)) - 2._wp*DELAD

                        p_source%sl(jc,jb,tr_idx(15,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(15,jf,KH,K)) + AD*FKLAMP1
                        p_source%sl(jc,jb,tr_idx(16,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(16,jf,KH,K)) + AD*FKLAMP2
                        p_source%fl(jc,jb,tr_idx(15,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(15,jf,KH,K)) + DELAP*FKLAP12
                        p_source%fl(jc,jb,tr_idx(16,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(16,jf,KH,K)) + DELAP*FKLAP22
                      END DO MIR3
                    END DO !jc
                  END DO DIR3!  BRANCH BACK TO 3.1.1 FOR NEXT DIRECTION.
                END IF
              ELSE!IF (MP .LE.ML) THEN
                !     4.1.1   ANGULAR LOOP.                                   !
                DIR4: DO K = 1, ndirs!

                  DO jc = i_startidx, i_endidx
                    !     4.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND   !
                    !             DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.       !
                    !             ----------------------------------------------  !

                    !     4.1 LOOP FOR ANGULAR SYMMETRY.
                    MIR4: DO KH = 1,2

                      SAP = &
                           GW1*tracer(jc,jk,jb,tr_idx(1,jf,KH,K)) + &
                           GW2*tracer(jc,jk,jb,tr_idx(2,jf,KH,K)) + &
                           GW3*tracer(jc,jk,jb,tr_idx(3,jf,KH,K)) + &
                           GW4*tracer(jc,jk,jb,tr_idx(4,jf,KH,K))
                      SAM = &
                           GW5*tracer(jc,jk,jb,tr_idx(5,jf,KH,K)) + &
                           GW6*tracer(jc,jk,jb,tr_idx(6,jf,KH,K)) + &
                           GW7*tracer(jc,jk,jb,tr_idx(7,jf,KH,K)) + &
                           GW8*tracer(jc,jk,jb,tr_idx(8,jf,KH,K))

                      FTEMP = p_diag%AF11(jf) * p_diag%ENH(jc,jb)
                      FIJ = tracer(jc,jk,jb,tr_idx(9,jf,KH,K))*FTAIL
                      FAD1 = FIJ*(SAP+SAM)
                      FAD2 = FAD1-2._wp*SAP*SAM
                      FAD1 = FAD1+FAD2
                      FCEN = FTEMP*FIJ
                      AD = FAD2*FCEN
                      DELAD = FAD1*FTEMP
                      DELAP = (FIJ-2._wp*SAM)*wc%DAL1*FCEN
                      DELAM = (FIJ-2._wp*SAP)*wc%DAL2*FCEN

                      p_source%sl(jc,jb,tr_idx(10,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(10,jf,KH,K)) + AD*FKLAMM1
                      p_source%sl(jc,jb,tr_idx(11,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(11,jf,KH,K)) + AD*FKLAMM2
                      p_source%fl(jc,jb,tr_idx(10,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(10,jf,KH,K)) + DELAM*FKLAM12
                      p_source%fl(jc,jb,tr_idx(11,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(11,jf,KH,K)) + DELAM*FKLAM22

                      p_source%sl(jc,jb,tr_idx(12,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(12,jf,KH,K)) + AD*FKLAMMA
                      p_source%sl(jc,jb,tr_idx(13,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(13,jf,KH,K)) + AD*FKLAMMB
                      p_source%fl(jc,jb,tr_idx(12,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(12,jf,KH,K)) + DELAM*FKLAMA2
                      p_source%fl(jc,jb,tr_idx(13,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(13,jf,KH,K)) + DELAM*FKLAMB2

                      p_source%sl(jc,jb,tr_idx(14,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(14,jf,KH,K)) - 2._wp*AD
                      p_source%fl(jc,jb,tr_idx(14,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(14,jf,KH,K)) - 2._wp*DELAD
                    END DO MIR4
                  END DO !jc
                END DO DIR4
              END IF
            ELSE
              !     5.1.1   ANGULAR LOOP.                                         !
              DIR5: DO K = 1, ndirs

                DO jc = i_startidx, i_endidx

                  !     5.1 LOOP FOR ANLULAR SYMMETRY.                                !
                  MIR5: DO KH = 1,2

                    SAP = &
                         GW1*tracer(jc,jk,jb,tr_idx(1,jf,KH,K)) + &
                         GW2*tracer(jc,jk,jb,tr_idx(2,jf,KH,K)) + &
                         GW3*tracer(jc,jk,jb,tr_idx(3,jf,KH,K)) + &
                         GW4*tracer(jc,jk,jb,tr_idx(4,jf,KH,K))
                    SAM = &
                         GW5*tracer(jc,jk,jb,tr_idx(5,jf,KH,K)) + &
                         GW6*tracer(jc,jk,jb,tr_idx(6,jf,KH,K)) + &
                         GW7*tracer(jc,jk,jb,tr_idx(7,jf,KH,K)) + &
                         GW8*tracer(jc,jk,jb,tr_idx(8,jf,KH,K))

                    FTEMP = p_diag%AF11(jf) * p_diag%ENH(jc,jb)
                    FIJ = tracer(jc,jk,jb,tr_idx(9,jf,KH,K))*FTAIL
                    FAD1 = FIJ*(SAP+SAM)
                    FAD2 = FAD1-2._wp*SAP*SAM
                    FAD1 = FAD1+FAD2
                    FCEN = FTEMP*FIJ
                    AD = FAD2*FCEN
                    DELAD = FAD1*FTEMP
                    DELAP = (FIJ-2._wp*SAM)*wc%DAL1*FCEN
                    DELAM = (FIJ-2._wp*SAP)*wc%DAL2*FCEN

                    p_source%sl(jc,jb,tr_idx(10,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(10,jf,KH,K)) + AD*FKLAMM1
                    p_source%sl(jc,jb,tr_idx(11,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(11,jf,KH,K)) + AD*FKLAMM2
                    p_source%fl(jc,jb,tr_idx(10,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(10,jf,KH,K)) + DELAM*FKLAM12
                    p_source%fl(jc,jb,tr_idx(11,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(11,jf,KH,K)) + DELAM*FKLAM22

                    p_source%sl(jc,jb,tr_idx(12,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(12,jf,KH,K)) + AD*FKLAMMA
                    p_source%sl(jc,jb,tr_idx(13,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(13,jf,KH,K)) + AD*FKLAMMB
                    p_source%fl(jc,jb,tr_idx(12,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(12,jf,KH,K)) + DELAM*FKLAMA2
                    p_source%fl(jc,jb,tr_idx(13,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(13,jf,KH,K)) + DELAM*FKLAMB2
                  END DO MIR5
                END DO !jc
              END DO DIR5
            END IF
          ELSE
            !     6.1.1   ANGULAR LOOP.                                        !
            DIR6: DO K = 1,ndirs

              DO jc = i_startidx, i_endidx

                !     6.1 LOOP FOR ANGULAR SYMMETRY.                                !
                MIR6: DO KH = 1,2

                  SAP = &
                       GW1*tracer(jc,jk,jb,tr_idx(1,jf,KH,K)) + &
                       GW2*tracer(jc,jk,jb,tr_idx(2,jf,KH,K)) + &
                       GW3*tracer(jc,jk,jb,tr_idx(3,jf,KH,K)) + &
                       GW4*tracer(jc,jk,jb,tr_idx(4,jf,KH,K))
                  SAM = &
                       GW5*tracer(jc,jk,jb,tr_idx(5,jf,KH,K)) + &
                       GW6*tracer(jc,jk,jb,tr_idx(6,jf,KH,K)) + &
                       GW7*tracer(jc,jk,jb,tr_idx(7,jf,KH,K)) + &
                       GW8*tracer(jc,jk,jb,tr_idx(8,jf,KH,K))

                  FTEMP = p_diag%AF11(jf) * p_diag%ENH(jc,jb)
                  FIJ = tracer(jc,jk,jb,tr_idx(9,jf,KH,K))*FTAIL
                  FAD1 = FIJ*(SAP+SAM)
                  FAD2 = FAD1-2._wp*SAP*SAM
                  FAD1 = FAD1+FAD2
                  FCEN = FTEMP*FIJ
                  AD = FAD2*FCEN
                  DELAD = FAD1*FTEMP
                  DELAP = (FIJ-2._wp*SAM)*wc%DAL1*FCEN
                  DELAM = (FIJ-2._wp*SAP)*wc%DAL2*FCEN

                  p_source%sl(jc,jb,tr_idx(10,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(10,jf,KH,K)) + AD*FKLAMM1
                  p_source%sl(jc,jb,tr_idx(11,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(11,jf,KH,K)) + AD*FKLAMM2
                  p_source%fl(jc,jb,tr_idx(10,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(10,jf,KH,K)) + DELAM*FKLAM12
                  p_source%fl(jc,jb,tr_idx(11,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(11,jf,KH,K)) + DELAM*FKLAM22
                END DO MIR6
              END DO !jc
            END DO DIR6
          END IF
        ELSE
          !     7.1.1   ANGULAR LOOP.                                           !
          DIR7: DO K = 1,ndirs

            DO jc = i_startidx, i_endidx

              !     7.1 LOOP FOR ANGULAR SYMMETRY.                                     !
              MIR7: DO KH = 1,2

                SAP = &
                     GW1*tracer(jc,jk,jb,tr_idx(1,jf,KH,K)) + &
                     GW2*tracer(jc,jk,jb,tr_idx(2,jf,KH,K)) + &
                     GW3*tracer(jc,jk,jb,tr_idx(3,jf,KH,K)) + &
                     GW4*tracer(jc,jk,jb,tr_idx(4,jf,KH,K))

                FTEMP = p_diag%AF11(jf) * p_diag%ENH(jc,jb)
                FIJ = tracer(jc,jk,jb,tr_idx(9,jf,KH,K))
                FAD2 = FIJ*SAP
                FAD1 = 2._wp*FAD2
                FCEN = FTEMP*FIJ
                AD = FAD2*FCEN
                DELAD = FAD1*FTEMP
                DELAP = FIJ*wc%DAL1*FCEN

                p_source%sl(jc,jb,tr_idx(14,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(14,jf,KH,K)) - 2._wp*AD
                p_source%sl(jc,jb,tr_idx(15,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(15,jf,KH,K)) + AD*FKLAMP1
                p_source%sl(jc,jb,tr_idx(16,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(16,jf,KH,K)) + AD*FKLAMP2
                p_source%sl(jc,jb,tr_idx(17,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(17,jf,KH,K)) + AD*FKLAMPA
                p_source%sl(jc,jb,tr_idx(18,jf,KH,K)) = p_source%sl(jc,jb,tr_idx(18,jf,KH,K)) + AD*FKLAMPB

                p_source%fl(jc,jb,tr_idx(14,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(14,jf,KH,K)) - 2._wp*DELAD
                p_source%fl(jc,jb,tr_idx(15,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(15,jf,KH,K)) + DELAP*FKLAP12
                p_source%fl(jc,jb,tr_idx(16,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(16,jf,KH,K)) + DELAP*FKLAP22
                p_source%fl(jc,jb,tr_idx(17,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(17,jf,KH,K)) + DELAP*FKLAPA2
                p_source%fl(jc,jb,tr_idx(18,jf,KH,K)) = p_source%fl(jc,jb,tr_idx(18,jf,KH,K)) + DELAP*FKLAPB2
              END DO MIR7
            END DO !jc
          END DO DIR7
        END IF
      END DO FRE4
    END DO ljb
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE src_nonlinear_transfer

END MODULE mo_wave_source
