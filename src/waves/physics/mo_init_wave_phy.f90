! Description:  Contains the data structures
! for initialisation of the physical model state and other auxiliary variables
! in order to run wave physics.
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

MODULE mo_init_wave_physics

  USE mo_kind,                 ONLY: wp
  USE mo_mpi,                  ONLY: my_process_is_stdio
  USE mo_exception,            ONLY: message, message_text, finish
  USE mo_model_domain,         ONLY: t_patch
  USE mo_impl_constants,       ONLY: MAX_CHAR_LENGTH, min_rlcell, SUCCESS
  USE mo_physical_constants,   ONLY: grav
  USE mo_math_constants,       ONLY: pi2, rpi_2, deg2rad, rad2deg
  USE mo_loopindices,          ONLY: get_indices_c

  USE mo_wave_types,           ONLY: t_wave_prog, t_wave_diag
  USE mo_wave_forcing_types,   ONLY: t_wave_forcing
  USE mo_wave_config,          ONLY: t_wave_config
  USE mo_wave_ext_data_types,  ONLY: t_external_wave
  USE mo_wave_constants,       ONLY: EMIN
  USE mo_wave_physics,         ONLY: wave_group_velocity_c,   &
    &                                wave_group_velocity_e,   &
    &                                wave_group_velocity_nt,  &
    &                                wave_group_velocity_bnd, &
    &                                wave_number_c, wave_number_e

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_wave_phy

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_init_wave_physics'

CONTAINS

  !>
  !! Initialisation of the wave physics
  !!
  !! 1. Calculation of initial spectrum from the fetch law
  !! and 1D JONSWAP spectum. The minimum of wave energy
  !! is limited to FLMIN.
  !! 2. tba
  !!
  SUBROUTINE init_wave_phy(p_patch, wave_config, p_prog, p_diag, wave_ext_data, p_forcing)

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    TYPE(t_wave_prog),           INTENT(INOUT) :: p_prog
    TYPE(t_wave_diag),           INTENT(INOUT) :: p_diag
    TYPE(t_external_wave),       INTENT(IN)    :: wave_ext_data
    TYPE(t_wave_forcing),        INTENT(IN)    :: p_forcing

    TYPE(t_wave_config), POINTER :: wc => NULL()

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':init_wave_phy'

    ! save some paperwork
    wc => wave_config

    CALL FETCH_LAW(p_patch, wave_config, p_diag, p_forcing)

    ! Set minimum values of energy allowed in the spectrum
    CALL JONSWAP(p_patch,           &
         wc%freqs,                  &
         p_diag%ALPHAJ*0.01_wp,     &
         wc%GAMMA_wave, wc%SIGMA_A, wc%SIGMA_B, &
         p_diag%FP,                 &
         p_diag%FLMINFR)           !OUT

    ! Set JONSWAP spectrum
    CALL JONSWAP(p_patch,           &
         wc%freqs,                  &
         p_diag%ALPHAJ,             &
         wc%GAMMA_wave, wc%SIGMA_A, wc%SIGMA_B, &
         p_diag%FP,                 &
         p_diag%ET)           !OUT

    CALL init_wave_spectrum(p_patch, wc, p_diag, p_forcing, p_prog%tracer)

    ! get wave number as a function of circular frequency and water depth
    ! at cell center
    CALL wave_number_c(p_patch     = p_patch,                    & !IN
      &                wave_config = wave_config,                & !IN
      &                depth       = wave_ext_data%bathymetry_c, & !IN
      &                wave_num_c  = p_diag%wave_num_c)            !OUT

    ! get wave number as a function of circular frequency and water depth
    ! at edge midpoint
    CALL wave_number_e(p_patch     = p_patch,                    & !IN
      &                wave_config = wave_config,                & !IN
      &                depth       = wave_ext_data%bathymetry_e, & !IN
      &                wave_num_e  = p_diag%wave_num_e)            !OUT


    ! compute absolute value of group velocity at cell centers
    !
    CALL wave_group_velocity_c(p_patch, wc, &
         p_diag%wave_num_c,          &  ! IN
         wave_ext_data%bathymetry_c, &  ! IN
         p_diag%gv_c)                   !INOUT

    ! compute absolute value of group velocity at edge midpoints
    !
    CALL wave_group_velocity_e(p_patch, wc, &
         p_diag%wave_num_e,          &  ! IN
         wave_ext_data%bathymetry_e, &  ! IN
         p_diag%gv_e)                   !INOUT

    ! compute normal and tangential components of group velociy vector
    ! at edge midpoints
    CALL wave_group_velocity_nt(p_patch, wc, &
         p_diag%gv_e, &
         p_diag%gvn_e, & !INOUT
         p_diag%gvt_e)   !INOUT

    CALL wave_group_velocity_bnd(p_patch, wc, &
         p_diag%gvn_e)  !INOUT

    ! initialisation of the nonlinear transfer
    CALL init_wave_nonlinear(wave_config = wave_config,                & !IN
         &                   p_diag      = p_diag) !INOUT



    CALL message(TRIM(routine),'finished')

  END SUBROUTINE init_wave_phy


  !>
  !! Initialisation of the wave spectrum
  !!
  !! Calculation of wind dependent initial spectrum from
  !! the fetch law and from the 1D JONSWAP spectum. The minimum
  !! of wave energy is limited to FLMIN.
  !!
  SUBROUTINE init_wave_spectrum(p_patch, wave_config, p_diag, p_forcing, tracer)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//':init_wave_spectrum'

    TYPE(t_patch),        INTENT(IN)    :: p_patch
    TYPE(t_wave_config),  INTENT(IN)    :: wave_config
    TYPE(t_wave_diag),    INTENT(IN)    :: p_diag
    TYPE(t_wave_forcing), INTENT(IN)    :: p_forcing
    REAL(wp),             INTENT(INOUT) :: tracer(:,:,:,:)

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jd,jf,jk,jt
    REAL(wp):: st

    REAL(wp), PARAMETER :: FLMIN = 0.000001_wp !! absolute minimum energy in spectral bins

    ! halo points must be included
    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)
    jk         = p_patch%nlev


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jd,jt,jc,i_startidx,i_endidx,st)
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1,wave_config%nfreqs
        DO jd = 1,wave_config%ndirs
          !
          jt = wave_config%tracer_ind(jd,jf)
          !
        DO jc = i_startidx, i_endidx
            st = rpi_2*MAX(0._wp, COS(wave_config%dirs(jd)-p_forcing%dir10m(jc,jb)*deg2rad) )**2
            IF (st < 0.1E-08_wp) st = 0._wp

            ! Avoid too small numbers of p_diag%ET
            p_diag%ET(jc,jb,jf) = MAX(p_diag%ET(jc,jb,jf),FLMIN)

            ! WAM initialisation
            tracer(jc,jk,jb,jt) = p_diag%ET(jc,jb,jf) * st
            tracer(jc,jk,jb,jt) = MAX(tracer(jc,jk,jb,jt),EMIN)
          END DO  !jc
        END DO  !jd
      END DO  !jf
    END DO  !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CALL message(TRIM(routine),'finished')

  END SUBROUTINE init_wave_spectrum

  !>
  !! Calculation of the JONSWAP spectrum according to
  !! Hasselmann et al. 1973. Adaptation of WAM 4.5
  !! subroutine JONSWAP.
  !!
  SUBROUTINE JONSWAP (p_patch, freqs, ALPHAJ, GAMMA, SA, SB, FP, ET)
    TYPE(t_patch), INTENT(IN)  :: p_patch
    REAL(wp),      INTENT(IN)  :: freqs(:)      !! FREQUENCiIES.
    REAL(wp),      INTENT(IN)  :: ALPHAJ(:,:)   !! OVERALL ENERGY LEVEL OF JONSWAP SPECTRA.
    REAL(wp),      INTENT(IN)  :: GAMMA         !! OVERSHOOT FACTOR.
    REAL(wp),      INTENT(IN)  :: SA            !! LEFT PEAK WIDTH.
    REAL(wp),      INTENT(IN)  :: SB            !! RIGHT PEAK WIDTH.
    REAL(wp),      INTENT(IN)  :: FP(:,:)       !! PEAK FREQUENCIES.
    REAL(wp),      INTENT(OUT) :: ET(:,:,:)     !! JONSWAP SPECTRA.

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//':JONSWAP'

    REAL(wp) :: ARG, sigma, G2ZPI4FRH5M

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jc,i_startidx,i_endidx,G2ZPI4FRH5M,sigma,ARG) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1,SIZE(freqs)

        G2ZPI4FRH5M = grav**2 / pi2**4 * freqs(jf)**(-5)

        DO jc = i_startidx, i_endidx

          sigma = MERGE(sb,sa, freqs(jf)>fp(jc,jb))
          ET(jc,jb,jf) = 0._wp

          ARG = 1.25_wp*(FP(jc,jb)/freqs(jf))**4
          IF (ARG.LT.50.0_wp) THEN
            ET(jc,jb,jf) = ALPHAJ(jc,jb) * G2ZPI4FRH5M * EXP(-ARG)
          END IF

          ARG = 0.5_wp*((freqs(jf)-FP(jc,jb)) / (sigma*FP(jc,jb)))**2
          IF (ARG.LT.99._wp) THEN
            ET(jc,jb,jf) = ET(jc,jb,jf)*exp(log(GAMMA)*EXP(-ARG))
          END IF

        END DO  !jc
      END DO  !jf
    END DO  !jb
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

    CALL message(TRIM(routine),'finished')

  END SUBROUTINE JONSWAP

  !>
  !! Calculation of JONSWAP parameters.
  !!
  !! Calculate the peak frequency from a fetch law
  !! and the JONSWAP alpha.
  !!
  !! Developted by S. Hasselmann (July 1990) and H. Guenther (December 1990).
  !! K.HASSELMAN,D.B.ROOS,P.MUELLER AND W.SWELL. A parametric wave prediction
  !! model. Journal of physical oceanography, Vol. 6, No. 2, March 1976.
  !!
  !! Adopted from WAM 4.5.
  !!
  SUBROUTINE FETCH_LAW (p_patch, wave_config, p_diag, p_forcing)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//':FETCH_LAW'
    !
    TYPE(t_patch),        INTENT(IN)    :: p_patch
    TYPE(t_wave_config),  INTENT(IN)    :: wave_config
    TYPE(t_wave_diag),    INTENT(INOUT) :: p_diag
    TYPE(t_wave_forcing), INTENT(IN)    :: p_forcing

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc, jb

    !@waves move to nml?
    REAL(wp), PARAMETER :: A = 2.84_wp,  D = -(3._wp/10._wp) !! PEAK FREQUENCY FETCH LAW CONSTANTS
    REAL(wp), PARAMETER :: B = 0.033_wp, E = 2._wp/3._wp     !! ALPHA-PEAK FREQUENCY LAW CONSTANTS

    REAL(wp) :: UG
    REAL(wp) :: fetch
    REAL(wp) :: fm

    ! halo points must be included !
    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    fetch = wave_config%fetch
    fm    = wave_config%fm

    ! ---------------------------------------------------------------------------- !
    !                                                                              !
    !     1. COMPUTE VALUES FROM FETCH LAWS.                                       !
    !        -------------------------------
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,UG)
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jc = i_startidx, i_endidx
        IF (p_forcing%sp10m(jc,jb).GT.0.1E-08_wp) THEN
          UG = grav / p_forcing%sp10m(jc,jb)
          p_diag%FP(jc,jb) = MAX(0.13_wp, A*((grav*fetch)/(p_forcing%sp10m(jc,jb)**2))**D)

          p_diag%FP(jc,jb) = MIN(p_diag%FP(jc,jb), fm/UG)
          p_diag%ALPHAJ(jc,jb) = MAX(0.0081_wp, B * p_diag%FP(jc,jb)**E)
          p_diag%FP(jc,jb) = p_diag%FP(jc,jb) * UG
        ELSE
          p_diag%ALPHAJ(jc,jb) = 0.0081_wp
          p_diag%FP(jc,jb) = fm
        END IF

      END DO

    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

    CALL message(TRIM(routine),'finished')

  END SUBROUTINE FETCH_LAW

  !>
  !! Calculation of index arrays and weights for the computation of
  !! the nonlinear transfer rate for shallow water.
  !!
  !! Computation of parameters used in discrete interaction
  !! parameterization of nonlinear transfer.
  !!
  !! Adoptation of NLWEIGT from WAM 4.5.
  !!
  !! SUSANNE HASSELMANN JUNE 86.
  !! H. GUNTHER   ECMWF/GKSS  DECEMBER 90 - CYCLE_4 MODIFICATIONS.
  !! P. Janssen   ECMWF June 2005                                         !
  !! H. Gunther   HZG   January 2015  cycle_4.5.4
  !!
  !! Reference
  !! S. Hasselmann and K. Hasselmann, JPO, 1985
  !!
  SUBROUTINE init_wave_nonlinear(wave_config, p_diag)!, ext_data)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//':init_wave_nonlinear'

    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    TYPE(t_wave_diag),           INTENT(INOUT) :: p_diag

    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: jf,jd
    INTEGER :: error

    INTEGER :: nfreqs, ndirs
    INTEGER :: klp1, ic, kh, klh, k, ks, icl1, icl2, isg, k1, k11, k2, k21
    INTEGER :: m, ikn, i, ie
    INTEGER :: mc,im,im1,ip,ip1,mm,mm1,mp,mp1,mct

    REAL(wp) :: alamd, con, delphi1, delphi2
    REAL(wp) :: deltha, cl1, cl2, al11, al12, ch, cl1h, cl2h
    REAL(wp) :: f1p1, frg, flp, flm, fkp, fkm

    REAL(WP), ALLOCATABLE, DIMENSION(:) :: frlon


    wc => wave_config
    nfreqs = wc%nfreqs
    ndirs = wc%ndirs

    ALLOCATE(frlon(2*nfreqs+2), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

    ! Parameters for discrete approximation of nonlinear transfer
    alamd   = 0.25_wp   ! lambda
    con     = 3000.0_wp ! weight for discrete approximation of nonlinear transfer
    delphi1 = -11.48_wp
    delphi2 = 33.56_wp

    ! 1. Computation for angular grid
    deltha = wc%delth * rad2deg

    cl1 = delphi1/deltha
    cl2 = delphi2/deltha

    ! 1.1 computation of indices of angular cell.
    klp1 = ndirs+1
    ic = 1

    DO kh = 1,2
       klh = ndirs
       IF (kh.eq.2) klh=klp1
       DO k = 1,klh
          ks = k
          IF (kh.gt.1) ks=klp1-k+1
          IF (ks.gt.ndirs) CYCLE
          ch = ic*cl1
          p_diag%ja1(ks,kh) = jafu(ch,k,ndirs)
          ch = ic*cl2
          p_diag%ja2(ks,kh) = jafu(ch,k,ndirs)
       END DO
       ic = -1
    END DO

    ! 1.2 computation of angular weights
    icl1 = cl1
    cl1  = cl1 - icl1
    icl2 = cl2
    cl2  = cl2 - icl2
    wc%acl1 = ABS(cl1)
    wc%acl2 = ABS(cl2)
    wc%cl11 = 1._wp - wc%acl1
    wc%cl21 = 1._wp - wc%acl2
    al11 = (1._wp + alamd)**4
    al12 = (1._wp - alamd)**4
    wc%dal1 = 1._wp / al11
    wc%dal2 = 1._wp / al12

    ! 1.3 computation of angular indices
    isg = 1
    DO kh = 1,2
       cl1h = isg*cl1
       cl2h = isg*cl2
       DO k = 1,ndirs
          ks = k
          IF (kh.eq.2) ks = ndirs-k+2
          IF(k.eq.1) ks = 1
          k1 = p_diag%ja1(k,kh)
          p_diag%k1w(ks,kh) = k1
          IF (cl1h.lt.0.) THEN
             k11 = k1-1
             IF (k11.lt.1) k11 = ndirs
          ELSE
             k11 = k1+1
             IF (k11.gt.ndirs) k11 = 1
          END IF
          p_diag%k11w(ks,kh) = k11
          k2 = p_diag%ja2(k,kh)
          p_diag%k2w(ks,kh) = k2
          IF (cl2h.lt.0) THEN
             k21 = k2-1
             IF(k21.lt.1) k21 = ndirs
          ELSE
             k21 = k2+1
             IF (k21.gt.ndirs) k21 = 1
          END IF
          p_diag%k21w(ks,kh) = k21
       END DO
       isg = -1
    END DO

    ! 2. computation for frequency grid
    frlon(1:nfreqs) = wc%freqs(1:nfreqs)

    DO m = nfreqs+1,2*nfreqs+2
       frlon(m) = wc%co*frlon(m-1)
    END DO

    f1p1 = LOG10(wc%co)
    DO m = 1,nfreqs+4
       frg = frlon(m)
       p_diag%af11(m) = con * frg**11
       flp = frg*(1.+alamd)
       flm = frg*(1.-alamd)
       ikn = INT(LOG10(1._wp+alamd)/f1p1+.000001_wp)
       ikn = m+ikn
       p_diag%ikp(m) = ikn
       fkp = frlon(p_diag%ikp(m))
       p_diag%ikp1(m) = p_diag%ikp(m)+1
       p_diag%fklap(m) = (flp-fkp)/(frlon(p_diag%ikp1(m))-fkp)

       p_diag%fklap1(m) = 1._wp-p_diag%fklap(m)
       IF (frlon(1).ge.flm) THEN
          p_diag%ikm(m) = 1
          p_diag%ikm1(m) = 1
          p_diag%fklam(m) = 0._wp
          p_diag%fklam1(m) = 0._wp
       ELSE
          ikn = INT(LOG10(1._wp-alamd)/f1p1+.0000001_wp)
          ikn = m+ikn-1
          IF (ikn.lt.1) ikn = 1
          p_diag%ikm(m) = ikn
          fkm = frlon(p_diag%ikm(m))
          p_diag%ikm1(m) = p_diag%ikm(m)+1
          p_diag%fklam(m) = (flm-fkm)/(frlon(p_diag%ikm1(m))-fkm)

          p_diag%fklam1(m) = 1._wp-p_diag%fklam(m)
       END IF
    END DO


    ! 3. calculate nonlinear tracer index p_diag%non_lin_tr_ind(18,nfreqs+4,2,ndirs)
    FRE4: DO MC = 1,nfreqs+4
      MP  = p_diag%IKP (MC)
      MP1 = p_diag%IKP1(MC)
      MM  = p_diag%IKM (MC)
      MM1 = p_diag%IKM1(MC)
      IC  = MC
      IP  = MP
      IP1 = MP1
      IM  = MM
      IM1 = MM1
      IF (IP1.GT.nfreqs) THEN
        IP1 = nfreqs
        IF (IP.GT.nfreqs) THEN
          IP  = nfreqs
          IF (IC.GT.nfreqs) THEN
            IC  = nfreqs
            IF (IM1.GT.nfreqs) THEN
              IM1 = nfreqs
            END IF
          END IF
        END IF
      END IF

      MCT = MC
      IF (MCT.GT.nfreqs) MCT  = nfreqs
      IF (MM.GT.nfreqs)  MM  = nfreqs
      IF (MM1.GT.nfreqs) MM1 = nfreqs
      IF (MP.GT.nfreqs)  MP  = nfreqs
      IF (MP1.GT.nfreqs) MP1 = nfreqs

      !     2.1.1   ANGULAR LOOP.                                     !
      DIR2: DO K = 1,ndirs !DIR2
        MIR2: DO KH = 1,2 !MIR2

          K1  = p_diag%K1W(K,KH)
          K2  = p_diag%K2W(K,KH)
          K11 = p_diag%K11W(K,KH)
          K21 = p_diag%K21W(K,KH)

          p_diag%non_lin_tr_ind( 1,MC,KH,K) = wc%tracer_ind(K1,IP)
          p_diag%non_lin_tr_ind( 2,MC,KH,K) = wc%tracer_ind(K11,IP)
          p_diag%non_lin_tr_ind( 3,MC,KH,K) = wc%tracer_ind(K1,IP1)
          p_diag%non_lin_tr_ind( 4,MC,KH,K) = wc%tracer_ind(K11,IP1)
          p_diag%non_lin_tr_ind( 5,MC,KH,K) = wc%tracer_ind(K2,IM)
          p_diag%non_lin_tr_ind( 6,MC,KH,K) = wc%tracer_ind(K21,IM)
          p_diag%non_lin_tr_ind( 7,MC,KH,K) = wc%tracer_ind(K2,IM1)
          p_diag%non_lin_tr_ind( 8,MC,KH,K) = wc%tracer_ind(K21,IM1)
          p_diag%non_lin_tr_ind( 9,MC,KH,K) = wc%tracer_ind(K,IC)
          p_diag%non_lin_tr_ind(10,MC,KH,K) = wc%tracer_ind(K2 ,MM)
          p_diag%non_lin_tr_ind(11,MC,KH,K) = wc%tracer_ind(K21,MM)
          p_diag%non_lin_tr_ind(12,MC,KH,K) = wc%tracer_ind(K2 ,MM1)
          p_diag%non_lin_tr_ind(13,MC,KH,K) = wc%tracer_ind(K21,MM1)
          p_diag%non_lin_tr_ind(14,MC,KH,K) = wc%tracer_ind(K  ,MCT)
          p_diag%non_lin_tr_ind(15,MC,KH,K) = wc%tracer_ind(K1 ,MP)
          p_diag%non_lin_tr_ind(16,MC,KH,K) = wc%tracer_ind(K11,MP)
          p_diag%non_lin_tr_ind(17,MC,KH,K) = wc%tracer_ind(K1 ,MP1)
          p_diag%non_lin_tr_ind(18,MC,KH,K) = wc%tracer_ind(K11,MP1)
        END DO MIR2
      END DO DIR2
    END DO FRE4


    ! 3. compute tail frequency ratios
    wc%frh(:) = 0._wp ! initialisation

    ie = MIN(30,nfreqs+3)
    DO i = 1,ie
       m = nfreqs+i-1
       wc%frh(i) = (frlon(nfreqs)/frlon(m))**5
    END DO

    IF (ALLOCATED(frlon))          DEALLOCATE(frlon)

    !print nonlinear status
    IF (my_process_is_stdio()) THEN
       WRITE(0,'(/,'' -------------------------------------------------'')')
       WRITE(0,*)'        non linear interaction parameters'
       WRITE(0,'(  '' -------------------------------------------------'')')
       WRITE(0,'(/,''  frequency arrays'')')
       WRITE(0,'(''     acl1       acl2       cl11       cl21   '',                &
            &            ''    dal1       dal2'')')
       WRITE(0,'(1x,6f11.8)') wc%acl1, wc%acl2, wc%cl11, wc%cl21, wc%dal1, wc%dal2
       WRITE(0,*) ' '
       WRITE(0,'(''  m   ikp ikp1  ikm ikm1   fklap       fklap1 '',               &
            &            ''   fklam       fklam1     af11'')')

       DO jf = 1,size(p_diag%ikp)
          WRITE(0,'(1x,i2,4i5,4f11.8,e11.3)') jf, p_diag%ikp(jf), p_diag%ikp1(jf), p_diag%ikm(jf), p_diag%ikm1(jf), &
               &            p_diag%fklap(jf), p_diag%fklap1(jf), p_diag%fklam(jf), p_diag%fklam1(jf), p_diag%af11(jf)
       END DO

       WRITE(0,'(/,''  angular arrays'')')
       WRITE(0,'(''   |--------kh = 1----------||--------kh = 2----------|'')')
       WRITE(0,'(''  k   k1w   k2w  k11w  k21w   k1w   k2w  k11w  k21w'')')
       DO jd = 1,size(p_diag%k1w,1)
          WRITE(0,'(1x,i2,8i6)') jd,(p_diag%k1w(jd,kh), p_diag%k2w(jd,kh), p_diag%k11w(jd,kh),              &
               &                            p_diag%k21w(jd,kh),kh=1,2)
       END DO

       WRITE(0,'(/,''  tail array frh'')')
       WRITE(0,'(1x,8f10.7)') wc%frh(1:30)
    END IF


    CALL message(TRIM(routine),'finished')

  END SUBROUTINE init_wave_nonlinear


  !>
  !! Function to compute the index array for the angles of the
  !! interacting wavenumbers.
  !!
  !! Adopted from WAM 4.5 JAFU
  !!
  !!  S. Hasselmann        MPIFM        01/12/1985
  !!
  !! Indices defining bins in frequency and direction plane into
  !! which nonlinear energy transfer increments are stored. Needed
  !! for computation of the nonlinear energy transfer.
  !!
  !! Reference
  !! S. Hasselmann and K. Hasselmann, JPO, 1985 B
  INTEGER FUNCTION JAFU (CL, J, IAN)

    REAL(wp),    INTENT(IN) :: CL !! weights.
    INTEGER, INTENT(IN) :: J      !! index in angular array.
    INTEGER, INTENT(IN) :: IAN    !! number of angles in array.

    JAFU = J + INT(CL)
    IF (JAFU.LE.0)   JAFU = JAFU+IAN
    IF (JAFU.GT.IAN) JAFU = JAFU-IAN

  END FUNCTION JAFU


END MODULE mo_init_wave_physics
