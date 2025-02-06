! Module for hailcast diagnostic
!
! ---------------------------------------------------------------
! Copyright: Copyright (C) 2014-2018, AER, Rebecca Adams-Selin 
!
! Authors: Rebecca Adams-Selin
! Contact: RSelin@aer.com
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------
! 
! Code has been modified
! A full description of HAILCAST is provided in 
! Adams-Selin R.D. and C.L. Ziegler, 2016:
! Forecasting hail using a  one-dimensional hail growth model 
! within WRF. Mon. Wea. Rev., 144, 4919-4939.
! 

MODULE mo_diag_hailcast

!  Inputs:
!    1-d (nz)
!     TCA          temperature (K) 
!     h1d          height above sea level (m) 
!     PA           total pressure (Pa)
!     rho1d        density (kg/m3)
!     RA           vapor mixing ratio (kg/kg)
!     qi1d         cloud ice mixing ratio (kg/kg)
!     qc1d         cloud water mixing ratio (kg/kg)
!     qr1d         rain water mixing ratio (kg/kg)
!     qg1d         graupel mixing ratio (kg/kg)
!     qs1d         snow mixing ratio (kg/kg)
!     VUU          updraft speed at each level (m/s)
!    Float
!     ht         terrain height (m)
!     wdur       duration of any updraft > 10 m/s within 1 surrounding 
!                 gridpoint 
!    Integer
!     nz         number of vertical levels
!
!  Output:
!     dhail      hail diameter in mm 
!                1st-5th rank-ordered hail diameters returned

  USE mo_kind,                  ONLY: dp, wp
  USE mo_impl_constants,        ONLY: min_rlcell_int, grf_bdywidth_c
  USE mo_model_domain,          ONLY: t_patch
  USE mo_nonhydro_types,        ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_intp_data_strc,        ONLY: t_int_state
  USE mo_nwp_phy_types,         ONLY: t_nwp_phy_diag
  USE mo_parallel_config,       ONLY: nproma
  USE mo_physical_constants,    ONLY: grav, rd, rv, alv, als, alf
  USE mo_math_constants,        ONLY: pi
  USE mo_run_config,            ONLY: iqv, iqi, iqc, iqr, iqg, iqs
  USE mo_nh_pzlev_config,       ONLY: t_nh_pzlev_config
  USE mo_loopindices,           ONLY: get_indices_c
  USE mo_fortran_tools,         ONLY: init 

!==============================================================================

IMPLICIT NONE

PRIVATE

PUBLIC :: hailstone_driver

!==============================================================================


CONTAINS

  SUBROUTINE hailstone_driver ( ptr_patch, p_metrics, p_prog, p_prog_rcf, p_diag, wdur,topography_c, wdur_min, &
                                dhail )

    TYPE(t_patch),      INTENT(IN)    :: ptr_patch             !< patch on which computation is performed
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics             !< metric variables
    TYPE(t_nh_prog),    INTENT(IN)    :: p_prog, p_prog_rcf    !< prognostic variables
    TYPE(t_nh_diag),    INTENT(IN)    :: p_diag                !< diagnostic variables
    
    REAL(wp),           INTENT(IN)    :: topography_c(:,:)     !< height of ground surface above sea level, dim: (nproma, nblks_c)
    REAL(wp),           INTENT(IN)    :: wdur(:,:)             !< updraft duration
    REAL(wp),           INTENT(IN)    :: wdur_min              !< minimal time an updraft has to persist to activate hailcast

    REAL(wp),           INTENT(INOUT)   :: dhail(:,:,:)          !< output variable, dim: (nproma,nhailstone,nblks_c)

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx


    !Local variables
    REAL (KIND=wp), DIMENSION( nproma,ptr_patch%nlev,ptr_patch%nblks_c ) :: &
      RIA,             &                         ! frozen content mixing ratio (kg/kg)
      RWA_new                                    ! adiabatic liquid content mixing ratio (kg/kg)

    INTEGER, DIMENSION(nproma,ptr_patch%nblks_c) ::             &
      KBAS                                       ! level of cloud base

    INTEGER, DIMENSION(nproma,5,ptr_patch%nblks_c) ::           &
      LEV                                        ! actual level of hailstone

    REAL (KIND=wp), DIMENSION(nproma,ptr_patch%nblks_c) ::      &
      ZBAS                                       ! height of cloudbase
    REAL (KIND=wp) RBAS                          ! mix ratio of cloud base
    REAL (KIND=wp) cwitot                        ! total cloud water, ice mix ratio
    REAL (KIND=wp) ZFZL, TFZL, WFZLP             ! height, temp, pressure of embryo start point
    REAL (KIND=wp) RFZL                          ! mix ratio of embryo start point
    REAL (KIND=wp) VUFZL, DENSAFZL               ! updraft speed, density of embryo start point
    INTEGER CLOUDON                              ! use to zero out cloud water, ice once past updraft duration
    REAL (KIND=wp) RTIME                         ! real updraft duration (sec)
    REAL (KIND=wp) TAU                           ! upper time limit of simulation (sec)
    !hailstone parameters
    REAL (KIND=wp) D, D_ICE                      ! hail diameter and ice diameter of it (m)
    REAL (KIND=wp) VT                            ! terminal velocity (m/s)
    REAL (KIND=wp) V                             ! actual stone velocity (m/s)
    !HAILSTONE temperature differencing
    REAL (KIND=wp), DIMENSION(nproma,5,ptr_patch%nblks_c) ::      &
      TS, TSm1, TSm2                             ! hailstone temperature of actual and two passed time steps (K)
    REAL (KIND=wp) FW                            ! fraction of stone that is liquid
    REAL (KIND=wp) WATER                         ! mass of stone that is liquid
    REAL (KIND=wp) CRIT                          ! critical water mass allowed on stone surface
    REAL (KIND=wp) DENSE                         ! hailstone density (kg/m3)
    INTEGER ITYPE                                ! wet (2) or dry (1) growth regime

    !variables of updraft characteristics for adiabatic water profile
    REAL (KIND=wp) RWA_adiabat                   ! adiabatic liquid content mixing ratio (kg/kg)
    REAL (KIND=wp) RSA                           ! saturation mixing ratio
    REAL (KIND=wp) ESVA                          ! saturation vapor pressure

    !in-cloud updraft parameters at location of hailstone
    REAL (KIND=wp) P                             ! in-cloud pressure (Pa)
    REAL (KIND=wp) RS                            ! in-cloud saturation mixing ratio 
    REAL (KIND=wp) RI, RW                        ! ice, liquid water mix. ratio (kg/kg)
    REAL (KIND=wp) XI, XW                        ! ice, liquid water content (kg/m3 air)
    REAL (KIND=wp) PC                            ! in-cloud fraction of frozen water
    REAL (KIND=wp) TC                            ! in-cloud temperature (K)
    REAL (KIND=wp) VU                            ! in-cloud updraft speed (m/s)
    REAL (KIND=wp) VUMAX                         ! in-cloud updraft speed read from WRF (max allowed)
    REAL (KIND=wp) VUCORE                        ! perturbed in-cloud updraft speed
    REAL (KIND=wp) DENSA                         ! in-cloud updraft density (kg/m3)
    REAL (KIND=wp) Z                             ! height of hailstone (m)
    REAL (KIND=wp) DELRW                         ! diff in sat vap. dens. between hail and air (kg/m3)
    !mean sub-cloud layer variables
    REAL (KIND=wp) TLAYER,RLAYER,PLAYER          ! mean sub-cloud temp, mix ratio, pres
    REAL (KIND=wp) TSUM,RSUM,PSUM                ! sub-cloud layer T, R, P sums
    REAL (KIND=wp) LDEPTH                        ! layer depth
    !internal function variables
    REAL (KIND=wp) GM,GM1,GMW,GMI,           &   ! mass variables of different parts of hailstone growth
      DGM,DGMW,DGMI,DGMV,DI,                 &   ! mass diff of different kind of growth
      ANU,RE,AE                                  ! some other parameter
    REAL (KIND=wp) dum,                      &   ! a dummy variable
      icefactor,                             &   ! factor of ice in the cloud
      VERH
    INTEGER dum_temp                             ! run time variable of internal timestep
    REAL (KIND=wp) sec, secdel                   ! time step, increment in seconds
    INTEGER i, j, k, IFOUT,                  &   ! some step variables 
      jc, jb, kk,                            &   ! horizontal and vertical step variables
      hs                                         ! hailstone number

    REAL (KIND=wp) PDIFF                         ! Interpollation variable for pressure
    REAL (KIND=wp) VDIFF                         ! Interpollation variable for arrays 
    REAL (KIND=wp) FRACT                         ! Interpollation variable of temperature

    REAL (KIND=wp) dhails                        ! hailstone at the end
    REAL (KIND=wp) tk_embryo                     ! initial embryo temperature - reflects the height (K)
    REAL (KIND=wp) DD                            ! initial embryo diameter (m)

    !$ACC DATA CREATE(ZBAS, KBAS, RWA_new, RIA, LEV, TS, TSm1, TSm2)

    !!!!!!!!!!!!!!!! 1. INITIAL EMBRYO !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!       FIND PROPERTIES OF INITIAL EMBRYO LOCATION       !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Find the cloud base for end-of-algorithm purposes.
    ! Start parallelizing


    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(RBAS, kk, cwitot, ESVA, RSA, icefactor, RWA_adiabat)
      DO jc = i_startidx, i_endidx 
        IF (wdur(jc,jb) .gt. wdur_min) THEN
          RWA_new(jc,:,jb) = 0.0_wp
          RIA(jc,:,jb) = 0.0_wp
          KBAS(jc,jb) = 1
          RSA = 0.0_wp
          ESVA = 0.0_wp

          !$ACC LOOP SEQ
          DO kk=ptr_patch%nlev,1,-1
            cwitot = p_prog_rcf%tracer(jc,kk,jb,iqi) + p_prog_rcf%tracer(jc,kk,jb,iqc)
            RIA(jc,kk,jb) = p_prog_rcf%tracer(jc,kk,jb,iqi) + p_prog_rcf%tracer(jc,kk,jb,iqs)
            IF ((cwitot .ge. 1.E-12_wp) .and. (kk .gt. KBAS(jc,jb))) THEN
              KBAS(jc,jb) = kk 
            ENDIF
          ENDDO

          !At coarser resolutions WRF underestimates liquid water aloft.
          !A fairer estimate is of the adiabatic liquid water profile, or 
          !the difference between the saturation mixing ratio at each level
          !and the cloud base water vapor mixing ratio
          ZBAS(jc,jb) = 0.5_wp * (p_metrics%z_ifc(jc,KBAS(jc,jb),jb) + p_metrics%z_ifc(jc,KBAS(jc,jb)+1,jb))
          RBAS = p_prog_rcf%tracer(jc,KBAS(jc,jb),jb,iqv)

          !$ACC LOOP SEQ
          DO kk=KBAS(jc,jb),1,-1
            !saturation vapor pressure
            ESVA = 611.2_wp*exp(17.67_wp*(p_diag%temp(jc,kk,jb)-273.155_wp)/(p_diag%temp(jc,kk,jb)-29.655_wp)) !Pa
            !saturation vapor mixing ratio
            RSA = 0.62197_wp * ESVA / (p_diag%pres(jc,kk,jb) - ESVA)
            !Up above -31, start converting to ice, entirely by -38
            !(Rosenfeld and Woodley 2000)
            IF (p_diag%temp(jc,kk,jb).gt.242.155_wp) THEN
              icefactor = 1._wp
            ELSE IF ((p_diag%temp(jc,kk,jb).LE.242.155_wp).AND.(p_diag%temp(jc,kk,jb).GT.235.155_wp)) THEN
              icefactor = (1._wp-(242.155_wp-p_diag%temp(jc,kk,jb))/5._wp)
            ELSE
              icefactor = 0._wp
            ENDIF
            !Don't want any negative liquid water values
            IF (RBAS.GT.RSA) THEN
              RWA_adiabat = (RBAS - RSA)*icefactor
            ELSE
              RWA_adiabat = p_prog_rcf%tracer(jc,kk,jb,iqc)
            ENDIF
            !Remove cloud liquid water outputted at previous lower levels
            IF (kk.eq.KBAS(jc,jb)) THEN
              RWA_new(jc,kk,jb) = RWA_adiabat
            ELSE IF ((kk.le.KBAS(jc,jb)-1).AND.(RWA_adiabat.ge.1.E-12_wp)) THEN
              RWA_new(jc,kk,jb) = RWA_adiabat* &
                                 (0.5_wp * (p_metrics%z_ifc(jc,kk,jb) + p_metrics%z_ifc(jc,kk+1,jb)) &
                                - 0.5_wp * (p_metrics%z_ifc(jc,kk+1,jb) + p_metrics%z_ifc(jc,kk+2,jb))) &
                                - RWA_new(jc,kk+1,jb)
              IF (RWA_new(jc,kk,jb).LT.0._wp) RWA_new(jc,kk,jb) = 0._wp
            ELSE
              RWA_new(jc,kk,jb) = p_prog_rcf%tracer(jc,kk,jb,iqc) 
            ENDIF
          ENDDO

          !remove the height factor from RWA_new
          !$ACC LOOP SEQ
          DO kk=KBAS(jc,jb)-1,1,-1
            RWA_new(jc,kk,jb) = RWA_new(jc,kk,jb) / (0.5_wp * (p_metrics%z_ifc(jc,kk,jb) + p_metrics%z_ifc(jc,kk+1,jb)) &
                               -0.5_wp * (p_metrics%z_ifc(jc,kk+1,jb) + p_metrics%z_ifc(jc,kk+2,jb)))
          ENDDO
        END IF
      END DO
      !$ACC END PARALLEL
    END DO

    !!!!!!!!!!!!!!!! 2. INITIAL EMBRYO SIZE  !!!!!!!!!!!!!!!!!!!!!
    !!!      SET CONSTANT RANGE OF INITIAL EMBRYO SIZES        !!!
    !!!                     START LOOP                         !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! See Adams-Selin and Ziegler 2016 MWR for explanation of why
    ! these sizes were picked.
    !Run each initial embryo size perturbation

    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(dhails, kk, dum_temp, RTIME, TFZL, IFOUT) &
      !$ACC   PRIVATE(FRACT, WFZLP, PDIFF, VDIFF, VERH, ZFZL, RFZL, VUFZL) &
      !$ACC   PRIVATE(DENSAFZL, sec, P, RS, TC, VU, Z, LDEPTH, DENSA) &
      !$ACC   PRIVATE(D, DENSE, ITYPE, CLOUDON, FW, DELRW) &
      !$ACC   PRIVATE(TAU, VUMAX, VUCORE, secdel, VT, V, RI, RW, XI, XW, PC) &
      !$ACC   PRIVATE(DGMV, GM1, DGM, GM, GMW, GMI, DGMW, DGMI, DI, ANU, RE, AE) &
      !$ACC   PRIVATE(D_ICE, WATER, CRIT, TSUM, RSUM, PSUM, TLAYER, PLAYER, RLAYER) &
      !$ACC   PRIVATE(DD, tk_embryo)
      DO hs = 1,5
        DO jc = i_startidx, i_endidx 
          IF (wdur(jc,jb) .gt. wdur_min) THEN
            IF (hs == 1) THEN
              DD = 5.E-3_wp
              tk_embryo = 265.155_wp
            ELSE IF (hs == 2) THEN
              DD = 7.5E-3_wp
              tk_embryo = 265.155_wp
            ELSE IF (hs == 3) THEN
              DD = 5.E-3_wp
              tk_embryo = 260.155_wp
            ELSE IF (hs == 4) THEN
              DD = 7.5E-3_wp
              tk_embryo = 260.155_wp
            ELSE IF (hs == 5) THEN
              DD = 1.E-2_wp
              tk_embryo = 260.155_wp
            END IF
            !!!!!!!!!!!!!!!! 3. UPDRAFT PROPERTIES !!!!!!!!!!!!!!!!!!!!!!!
            !!!              DEFINE UPDRAFT PROPERTIES                 !!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Set some initial default values
            ! Upper limit of simulation in seconds
            TAU = 7200.0_wp
            secdel = 5.0_wp

            ! Initial height of levels 
            LEV(jc,hs,jb) = ptr_patch%nlev - 1
            dhails = 0.0_wp

            ! Cap updraft lifetime at 2000 sec.
            RTIME = 2000._wp
            IF (wdur(jc,jb) .LT. RTIME) RTIME = wdur(jc,jb)
            kk = ptr_patch%nlev
            TFZL = tk_embryo
            IFOUT = 1

            ! Do interpolation for the hailstone between two levels
            ! for the pressure and the ratio for other arrays
            DO WHILE ((IFOUT == 1) .AND. (kk .GT. 1) .AND. (kk .LE. ptr_patch%nlev))
              IF ( (tk_embryo .LT. p_diag%temp(jc,kk,jb) .AND. tk_embryo .GE. p_diag%temp(jc,kk-1,jb)) .or.  &   ! dT/dz < 0
                   (tk_embryo .GT. p_diag%temp(jc,kk,jb) .AND. tk_embryo .LE. p_diag%temp(jc,kk-1,jb)) ) THEN    ! dT/dz > 0

                FRACT = (p_diag%temp(jc,kk,jb) - tk_embryo) / (p_diag%temp(jc,kk,jb) - p_diag%temp(jc,kk-1,jb))
                !.... compute the pressure value pval at temperature tval
                WFZLP = ((1.0_wp - FRACT) * p_diag%pres(jc,kk,jb)) + &
                        (FRACT * p_diag%pres(jc,kk-1,jb))
                IF ((WFZLP.LE. p_diag%pres(jc,kk,jb)) .AND. (WFZLP.GT.p_diag%pres(jc,kk-1,jb))) THEN
                  LEV(jc,hs,jb)   = kk
                  PDIFF = p_diag%pres(jc,kk,jb) - p_diag%pres(jc,kk-1,jb)
                  VDIFF = p_diag%pres(jc,kk,jb) - WFZLP
                  VERH  = VDIFF/PDIFF     
                  IFOUT = 0
                ENDIF
              ENDIF
              kk = kk - 1
            ENDDO

            ! Do interpolation for the hailstone between two levels
            ! for height(ZFZL), humidity(RFZL), Updraft(VUFZL), and density (DENSAFZL)
            ZFZL = 0.5_wp * (p_metrics%z_ifc(jc,LEV(jc,hs,jb),jb)+p_metrics%z_ifc(jc,LEV(jc,hs,jb)+1,jb)) + &
                  (0.5_wp * (p_metrics%z_ifc(jc,LEV(jc,hs,jb)-1,jb)+p_metrics%z_ifc(jc,LEV(jc,hs,jb),jb)) - &
                   0.5_wp * (p_metrics%z_ifc(jc,LEV(jc,hs,jb),jb)+p_metrics%z_ifc(jc,LEV(jc,hs,jb)+1,jb))) * VERH
            RFZL = p_prog_rcf%tracer(jc,LEV(jc,hs,jb),jb,iqv) + (p_prog_rcf%tracer(jc,LEV(jc,hs,jb)-1,jb,iqv) - &
                   p_prog_rcf%tracer(jc,LEV(jc,hs,jb),jb,iqv)) * VERH
            VUFZL = 0.5_wp * (p_prog%w(jc,LEV(jc,hs,jb),jb)+p_prog%w(jc,LEV(jc,hs,jb)+1,jb)) + &
                   (0.5_wp * (p_prog%w(jc,LEV(jc,hs,jb)-1,jb)+p_prog%w(jc,LEV(jc,hs,jb),jb)) - &
                    0.5_wp * (p_prog%w(jc,LEV(jc,hs,jb),jb)+p_prog%w(jc,LEV(jc,hs,jb)+1,jb))) * VERH
            DENSAFZL = p_prog%rho(jc,LEV(jc,hs,jb),jb) + (p_prog%rho(jc,LEV(jc,hs,jb)-1,jb) - &
                       p_prog%rho(jc,LEV(jc,hs,jb),jb)) * VERH

            !Begin hail simulation time (seconds)
            sec = 0.0_wp

            !!!!!!!!!!!!!!!!!!  4. INITIAL PARAMETERS    !!!!!!!!!!!!!!!!!
            !!!          PUT INITIAL PARAMETERS IN VARIABLES           !!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !Set initial values for parameters at freezing level
            P = WFZLP
            RS = RFZL
            TC = TFZL
            VU = VUFZL  
            Z = ZFZL - topography_c(jc,jb)
            LDEPTH = Z
            DENSA = DENSAFZL

            !Set initial hailstone parameters
            TS(jc,hs,jb) = TC
            TSm1(jc,hs,jb) = TS(jc,hs,jb)
            TSm2(jc,hs,jb) = TS(jc,hs,jb)      
            D = DD   !hailstone diameter in m
            FW = 0.0_wp
            DENSE = 500._wp  !kg/m3
            ITYPE = 1  !Assume starts in dry growth.
            CLOUDON = 1  !we'll eventually turn cloud "off" once updraft past time limit

            
            !Start time loop.
            !$ACC LOOP SEQ
            DO dum_temp = 0, nint(TAU), nint(secdel)
              IF (sec .lt. TAU) THEN
                sec = sec + secdel
                !!!!!!!!!!!!!!!!!!  5. CALCULATE PARAMETERS  !!!!!!!!!!!!!!!!!
                !!!              CALCULATE UPDRAFT PROPERTIES              !!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !Intepolate vertical velocity to our new pressure level
                VUMAX = 0.5_wp * (p_prog%w(jc,LEV(jc,hs,jb),jb)+p_prog%w(jc,LEV(jc,hs,jb)+1,jb)) + &
                       (0.5_wp * (p_prog%w(jc,LEV(jc,hs,jb)-1,jb)+p_prog%w(jc,LEV(jc,hs,jb),jb)) - &
                        0.5_wp * (p_prog%w(jc,LEV(jc,hs,jb),jb)+p_prog%w(jc,LEV(jc,hs,jb)+1,jb))) * VERH
         
                !Outside pressure levels?  If so, exit loop
                IF (IFOUT.EQ.1) sec = TAU
                  
                !Sine wave multiplier on updraft strength
                IF (SEC .GT. 0.0_wp .AND. SEC .LT. RTIME) THEN
                  VUCORE = VUMAX * (0.5_wp*SIN((pi*SEC)/(RTIME))+0.5_wp)*1.2_wp
                  VU = VUCORE      
                ELSEIF (SEC .GE. RTIME) THEN
                  VU = 0.0_wp
                  CLOUDON = 0
                ENDIF
         
                !Calculate terminal velocity of the hailstone 
                ! (use previous values)
                CALL TERMINL(DENSA,DENSE,D,VT,TC)
         
                !Actual velocity of hailstone (upwards positive)
                V = VU - VT
         
                !Use hydrostatic eq'n to calc height of next level
                P = P - DENSA*grav*V*secdel
                Z = Z + V*secdel

                ! Do interpolation for the hailstone between two levels
                ! for the ratio of other arrays
                IFOUT = 1
                IF ( V .GT. 0.0_wp ) THEN
                  kk = LEV(jc,hs,jb)
                  DO WHILE ((IFOUT == 1) .AND. (kk .GT. 1) .AND. (kk .LE. ptr_patch%nlev))
                    IF (P.LE. p_diag%pres(jc,kk,jb) .AND. P .GT. p_diag%pres(jc,kk-1,jb)) THEN
                      PDIFF = p_diag%pres(jc,kk,jb) - p_diag%pres(jc,kk-1,jb)
                      VDIFF = p_diag%pres(jc,kk,jb) - P
                      VERH  = VDIFF/PDIFF
                      LEV(jc,hs,jb)= kk
                      !Calculate ratio between vdiff and pdiff
                      IFOUT = 0
                    ENDIF
                    kk = kk - 1 
                  ENDDO
                ELSE
                  kk = LEV(jc,hs,jb)-1
                  DO WHILE ((IFOUT == 1) .AND. (kk .GT. 1) .AND. (kk .LE. ptr_patch%nlev))
                    IF (P.LE. p_diag%pres(jc,kk,jb) .AND. P .GT. p_diag%pres(jc,kk-1,jb)) THEN
                      PDIFF = p_diag%pres(jc,kk,jb) - p_diag%pres(jc,kk-1,jb)
                      VDIFF = p_diag%pres(jc,kk,jb) - P
                      VERH  = VDIFF/PDIFF
                      LEV(jc,hs,jb)= kk
                      !Calculate ratio between vdiff and pdiff
                      IFOUT = 0
                    ENDIF
                    kk = kk + 1
                  ENDDO
                END IF

                !Interpolate cloud temp, qvapor at new p-level

                TC = p_diag%temp(jc,LEV(jc,hs,jb),jb) + (p_diag%temp(jc,LEV(jc,hs,jb)-1,jb) - &
                     p_diag%temp(jc,LEV(jc,hs,jb),jb)) * VERH
                RS = p_prog_rcf%tracer(jc,LEV(jc,hs,jb),jb,iqv) + (p_prog_rcf%tracer(jc,LEV(jc,hs,jb) - 1,jb,iqv) - &
                     p_prog_rcf%tracer(jc,LEV(jc,hs,jb),jb,iqv)) * VERH

                !New density of in-cloud air
                DENSA=P/(rd*(1._wp+0.609_wp*RS/(1._wp+RS))*TC)

                !Interpolate liquid, frozen water mix ratio at new level
                RI = RIA(jc,LEV(jc,hs,jb),jb) + (RIA(jc,LEV(jc,hs,jb)-1,jb) - &
                     RIA(jc,LEV(jc,hs,jb),jb)) * VERH
                RW = RWA_new(jc,LEV(jc,hs,jb),jb) + (RWA_new(jc,LEV(jc,hs,jb)-1,jb) - &
                     RWA_new(jc,LEV(jc,hs,jb),jb)) * VERH
                XI = RI * DENSA * CLOUDON
                XW = RW * DENSA * CLOUDON
                IF( (XW+XI).GT.0._wp) THEN
                  PC = XI / (XW+XI)
                ELSE
                  PC = 1._wp
                ENDIF

                ! SATURATION VAPOUR DENSITY DIFFERENCE BETWTEEN STONE AND CLOUD
                CALL VAPORCLOSE(DELRW,TS(jc,hs,jb),TC,ITYPE)
        
                !!!!!!!!!!!!!!!!!!  6. STONE'S MASS GROWTH !!!!!!!!!!!!!!!!!!!!
                CALL MASSAGR(D,GM,GM1,GMW,GMI,DGM,DGMW,DGMI,DGMV,DI,ANU,RE,AE,&
                         TC,TS(jc,hs,jb),P,DENSE,DENSA,FW,VT,XW,XI,secdel,ITYPE,DELRW)

                !!!!!!!!!!!!!!!!!!  7. HEAT BUDGET OF HAILSTONE !!!!!!!!!!!!!!!
                CALL HEATBUD(TS(jc,hs,jb),TSm1(jc,hs,jb),TSm2(jc,hs,jb),FW,TC,DELRW,D,DENSA,GM1,GM,DGM,DGMW,  & 
                         DGMI,DI,RE,AE,secdel,ITYPE)

                !!!!! 8. TEST DIAMETER OF STONE AND HEIGHT ABOVE GROUND !!!!!!!
                !!!  TEST IF DIAMETER OF STONE IS GREATER THAN CRITICAL LIMIT, IF SO  
                !!!  BREAK UP 
                WATER=FW*GM  !KG
                ! CRTICAL MASS CAPABLE OF BEING "SUPPORTED" ON THE STONE'S SURFACE 
                CRIT = 2.0E-4_wp
                IF (WATER.GT.CRIT)THEN
                  CALL BREAKUP(DENSE,D,GM,FW,CRIT)
                ENDIF
        
                !!! Has stone reached below cloud base?
                IF (Z .LE. ZBAS(jc,jb)) sec = TAU

                !calculate ice-only diameter size
                D_ICE = ( (6._wp*GM*(1._wp-FW)) / (pi*DENSE) )**0.33333333_wp 

                !Has the stone entirely melted and it's below the freezing level?  
                IF ((D_ICE .LT. 1.E-8_wp) .AND. (TC.GT.273.155_wp)) sec = TAU

                !Did HAILCAST exceed the time limit?
                IF (sec .GE. TAU) THEN
                    LDEPTH = Z
                END IF

                !move values to previous timestep value
                TSm2(jc,hs,jb) = TSm1(jc,hs,jb)
                TSm1(jc,hs,jb) = TS(jc,hs,jb)

              END IF
            END DO  !end time loop


            !!!!!!!!!!!!!!!!!! 9. MELT STONE BELOW CLOUD !!!!!!!!!!!!!!!!!!!!
            !Did the stone shoot out the top of the storm? 
            !Then let's assume it's lost in the murky "outside storm" world.
            IF (P.lt. p_diag%pres(jc,1,jb)) THEN
              D=0.0_wp
            !Is the stone entirely water? Then set D=0 and exit.
            ELSE IF(ABS(FW - 1.0_wp).LT.0.001_wp) THEN
              D=0.0_wp
            ELSE IF (Z.GT.0_wp) THEN
              !If still frozen, then use melt routine to melt below cloud
              !based on mean below-cloud conditions.
       
              !Calculate mean sub-cloud layer conditions
              TSUM = 0._wp
              RSUM = 0._wp
              PSUM = 0._wp
              TLAYER = 0._wp
              PLAYER = 0._wp
              RLAYER = 0._wp
              DO kk = KBAS(jc,jb),ptr_patch%nlev
                TSUM = TSUM + p_diag%temp(jc,kk,jb)
                PSUM = PSUM + p_diag%pres(jc,kk,jb)
                RSUM = RSUM + p_prog_rcf%tracer(jc,kk,jb,iqv)
              ENDDO
              TLAYER = TSUM / (ptr_patch%nlev-KBAS(jc,jb)+1)
              PLAYER = PSUM / (ptr_patch%nlev-KBAS(jc,jb)+1)
              RLAYER = RSUM / (ptr_patch%nlev-KBAS(jc,jb)+1)

              !MELT is expecting a hailstone of only ice.  At the surface
              !we're only interested in the actual ice diameter of the hailstone,
              !so let's shed any excess water now.
              D_ICE = ( (6._wp*GM*(1._wp-FW)) / (pi*DENSE) )**0.33333333_wp 
              D = D_ICE  
              CALL MELT(D,TLAYER,PLAYER,RLAYER,LDEPTH,VT)

            ENDIF !end check for melting call

            !Check to make sure something hasn't gone horribly wrong
            IF (D.GT.0.254_wp) D = 0._wp  !just consider missing for now if > 10 in
      
            !assign hail size in mm for output
            dhails = D * 1000._wp

            IF (dhails .gt. dhail(jc,hs,jb)) THEN
              dhail(jc,hs,jb) = dhails
            END IF
          END IF
        END DO
      END DO
    !$ACC END PARALLEL
    END DO
    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE hailstone_driver



  SUBROUTINE TERMINL(DENSA,DENSE,D,VT,TC)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!
  !!!! INTERP: Calculate terminal velocity of the hailstone
  !!!!
  !!!! INPUT: DENSA  density of updraft air (kg/m3)
  !!!!        DENSE  density of hailstone
  !!!!        D      diameter of hailstone (m)
  !!!!        TC     updraft temperature (K)
  !!!! OUTPUT:VT     hailstone terminal velocity (m/s)
  !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !$ACC ROUTINE SEQ
      IMPLICIT NONE
      
      REAL(wp), INTENT(IN)  :: D, DENSA, DENSE, TC
      REAL(wp), INTENT(OUT)  :: VT

      !Local variables
      REAL(wp) :: GMASS, GX, RE, W, Y
      REAL(wp) :: ANU
      
      !Mass of stone in kg
      GMASS = (DENSE * PI * (D**3._wp)) / 6._wp
      
      !Dynamic viscosity
      ANU = (0.00001718_wp)*(273.155_wp+120._wp)/(TC+120._wp)*(TC/273.155_wp)**(1.5_wp)
      
      !CALC THE BEST NUMBER, X AND REYNOLDS NUMBER, RE 
      GX=(8.0_wp*GMASS*grav*DENSA)/(PI*(ANU*ANU))
      RE=(GX/0.6_wp)**0.5_wp

      !SELECT APPROPRIATE EQUATIONS FOR TERMINAL VELOCITY DEPENDING ON 
      !THE BEST NUMBER
      IF (GX.LT.550_wp) THEN
        W=LOG10(GX)
        Y= -1.7095_wp + 1.33438_wp*W - 0.11591_wp*(W**2.0_wp)      
        RE=10._wp**Y
        VT=ANU*RE/(D*DENSA)
      ELSE IF (GX.GE.550._wp.AND.GX.LT.1800._wp) THEN
        W=LOG10(GX)
        Y= -1.81391_wp + 1.34671_wp*W - 0.12427_wp*(W**2.0_wp) + 0.0063_wp*(W**3.0_wp)
        RE=10._wp**Y
        VT=ANU*RE/(D*DENSA)
      ELSE IF (GX.GE.1800._wp.AND.GX.LT.3.45E08_wp) THEN
        RE=0.4487_wp*(GX**0.5536_wp)
        VT=ANU*RE/(D*DENSA)
      ELSE
        RE=(GX/0.6_wp)**0.5_wp
        VT=ANU*RE/(D*DENSA)
      ENDIF
      
  END SUBROUTINE TERMINL   

   
   
  SUBROUTINE VAPORCLOSE(DELRW,TS,TC,ITYPE)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  VAPORCLOSE: CALC THE DIFFERENCE IN SATURATION VAPOUR DENSITY 
  !!!  BETWEEN THAT OVER THE HAILSTONE'S SURFACE AND THE IN-CLOUD 
  !!!  AIR, DEPENDS ON THE WATER/ICE RATIO OF THE UPDRAFT, 
  !!!  AND IF THE STONE IS IN WET OR DRY GROWTH REGIME
  !!!
  !!!  INPUT:  TS    temperature of hailstone (K)
  !!!          TC    temperature of updraft air (K)
  !!!          ITYPE wet (2) or dry (1) growth regime
  !!!  OUTPUT: DELRW difference in sat vap. dens. between hail and air
  !!!          (kg/m3)
  !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !$ACC ROUTINE SEQ

      IMPLICIT NONE
      REAL(wp), INTENT(OUT) :: DELRW
      REAL(wp), INTENT(IN)  :: TS, TC
      INTEGER,  INTENT(IN)  :: ITYPE
      !local variables
      REAL(wp) :: RATIO
      REAL(wp) :: ESAT, RHOKOR, ESATW, RHOOMGW, ESATI, RHOOMGI, RHOOMG

      !!!  FOR HAILSTONE:  FIRST TEST IF STONE IS IN WET OR DRY GROWTH
      RATIO = 1._wp/273.155_wp
      IF(ITYPE.EQ.2) THEN !!WET GROWTH
        ESAT=611._wp*EXP(alv/rv*(RATIO-1._wp/TS))
      ELSE  !!DRY GROWTH
        ESAT=611._wp*EXP(als/rv*(RATIO-1._wp/TS))
      ENDIF
      RHOKOR=ESAT/(rv*TS)
      
      !!!  NOW FOR THE AMBIENT/IN-CLOUD CONDITIONS 
      ESATW=611._wp*EXP(alv/rv*(RATIO-1._wp/TC))
      RHOOMGW=ESATW/(rv*TC)
      ESATI=611._wp*EXP(als/rv*(RATIO-1._wp/TC))
      RHOOMGI=ESATI/(rv*TC)
      RHOOMG = RHOOMGI  !done as in hailtraj.f

      !!!  CALC THE DIFFERENCE(KG/M3): <0 FOR CONDENSATION, 
      !!!  >0 FOR EVAPORATION
      DELRW=(RHOKOR-RHOOMG) 

  END SUBROUTINE VAPORCLOSE



  SUBROUTINE MASSAGR(D,GM,GM1,GMW,GMI,DGM,DGMW,DGMI,DGMV,DI,ANU,RE,AE,& 
                 TC,TS,P,DENSE,DENSA,FW,VT,XW,XI,SEKDEL,ITYPE,DELRW)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! CALC THE STONE'S INCREASE IN MASS 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !$ACC ROUTINE SEQ
      IMPLICIT NONE
      REAL(wp), INTENT(IN)    ::   P, TC, TS, DENSA, VT, &
                                   FW,XW,XI,SEKDEL,DELRW
      REAL(wp), INTENT(OUT)   ::   DGMV, GM1, DGM
      REAL(wp), INTENT(INOUT) ::   D
      REAL(wp), INTENT(INOUT) ::   DENSE
      REAL(wp), INTENT(INOUT) ::   GM,GMW,GMI,DGMW,DGMI,DI,ANU,RE,AE
      INTEGER,  INTENT(IN)    ::   ITYPE

      !local variables
      REAL(wp) :: D0, GMW2, GMI2, EW, EI
      REAL(wp) :: DENSEL, DENSELI, DENSELW 
      REAL(wp) :: DC !MEAN CLOUD DROPLET DIAMETER (MICRONS, 1E-6M)
      REAL(wp) :: VOLL, VOLT !VOLUME OF NEW LAYER, TOTAL (M3)
      REAL(wp) :: VOL1, DGMW_NOSOAK, SOAK, SOAKM
      REAL(wp) :: DENSAC, E, AFACTOR, NS, TSCELSIUS, VIMP, WIMP
      

      !!!  CALCULATE THE DIFFUSIVITY DI (m2/s)
      D0=0.226_wp*1.E-4_wp  ! change to m2/s, not cm2/s
      DI=D0*(TC/273.155_wp)**1.81_wp*(100000._wp/P)
  
      !!!  COLLECTION EFFICIENCY FOR WATER AND ICE 
      EW=1.0_wp
      !!!  Linear function for ice accretion efficiency
      IF (TC .GE. 273.155_wp) THEN
         EI=1.00_wp
      ELSE IF (TC.GE.233.155_wp) THEN
         EI=1.0_wp - ( (273.155_wp - TS) / 40._wp )
      ELSE  !cooler than -40C
         EI=0.0_wp
      ENDIF

      !!!  CALCULATE THE VENTILATION COEFFICIENT - NEEDED FOR GROWTH FROM VAPOR
      !The coefficients in the ventilation coefficient equations have been 
      !experimentally derived, and are expecting cal-C-g units.  Do some conversions.
      DENSAC = DENSA * (1.E3_wp) * (1.E-6_wp)
      !dynamic viscosity 
      ANU=1.717E-4_wp*(393.0_wp/(TC+120.0_wp))*(TC/273.155_wp)**1.5_wp
      !!!  CALCULATE THE REYNOLDS NUMBER - unitless
      RE=D*VT*DENSAC/ANU   
      E=(0.60_wp)**(0.333333333_wp)*(RE**0.50_wp) !ventilation coefficient vapor (fv)
      !!!   SELECT APPROPRIATE VALUES OF AE ACCORDING TO RE
      IF(RE.LT.6000.0_wp)THEN
         AE=0.78_wp+0.308_wp*E
      ELSEIF(RE.GE.6000.0_wp.AND.RE.LT.20000.0_wp)THEN
         AE=0.76_wp*E
      ELSEIF(RE.GE.20000.0_wp) THEN
         AE=(0.57_wp+9.0E-6_wp*RE)*E
      ENDIF

      !!!  CALC HAILSTONE'S MASS (GM), MASS OF WATER (GMW) AND THE  
      !!!  MASS OF ICE IN THE STONE (GMI)
      GM=PI/6._wp*(D**3._wp)*DENSE
      GMW=FW*GM
      GMI=GM-GMW
  
      !!!  STORE THE MASS
      GM1=GM
      
      !!! NEW MASS GROWTH CALCULATIONS WITH VARIABLE RIME 
      !!! LAYER DENSITY BASED ON ZIEGLER ET AL. (1983)
      
      !!! CALCULATE INCREASE IN MASS DUE INTERCEPTED CLD WATER, USE
      !!! ORIGINAL DIAMETER
      GMW2=GMW+SEKDEL*(PI/4._wp*D**2._wp*VT*XW*EW)
      DGMW=GMW2-GMW 
      GMW=GMW2

      !!! CALCULATE THE INCREASE IN MASS DUE INTERCEPTED CLOUD ICE
      GMI2=GMI+SEKDEL*(PI/4._wp*D**2._wp*VT*XI*EI)
      DGMI=GMI2-GMI 
      GMI=GMI2
  
      !!! CALCULATE INCREASE IN MASS DUE TO SUBLIMATION/CONDENSATION OF 
      !!! WATER VAPOR
      DGMV = SEKDEL*2._wp*PI*D*AE*DI*DELRW
      IF (DGMV .LT. 0._wp) DGMV=0._wp

      !!! CALCULATE THE TOTAL MASS CHANGE 
      DGM=DGMW+DGMI+DGMV

      !!! CALCULATE DENSITY OF NEW LAYER, DEPENDS ON FW AND ITYPE
      IF (ITYPE.EQ.1) THEN !DRY GROWTH
          !If hailstone encountered supercooled water, calculate new layer density 
          ! using Macklin form
          IF ((DGMW.GT.0._wp).OR.(DGMV.GT.0._wp)) THEN

            !MEAN CLOUD DROPLET RADIUS, ASSUME CLOUD DROPLET CONC OF 3E8 M-3 (300 CM-3)
            DC = (0.74_wp*XW / (PI*1000._wp*3.E8_wp))**0.33333333_wp * 1.E6_wp !MICRONS
            !!! FIND THE STOKES NUMBER  (rasmussen heymsfield 1985)
            NS = 2._wp*VT*100._wp*(DC*1.E-4_wp)**2._wp / (9._wp*ANU*D*50._wp)  !need hail radius in cm
            !!! FIND IMPACT VELOCITY (rasmussen heymsfield 1985)
            IF (NS .LT. 1.E-12_wp) THEN
              WIMP = 0._wp
            ELSE
              WIMP = LOG10(NS)
            ENDIF

            IF (RE.GT.200._wp) THEN
              IF (NS.LT.0.1_wp) THEN
                VIMP = 0.0_wp
              ELSEIF ((NS.GE.0.1_wp).AND.(NS.LE.10._wp)) THEN
                VIMP = (0.356_wp + 0.4738_wp*WIMP - 0.1233_wp*WIMP**2._wp &
                       -0.1618_wp*WIMP**3._wp + 0.0807_wp*WIMP**4._wp)*VT
              ELSEIF (NS.GT.10._wp) THEN
                VIMP = 0.63_wp*VT
              ENDIF
            ELSEIF ((RE.GT.65._wp).AND.(RE.LE.200._wp)) THEN
              IF (NS.LT.0.1_wp) THEN
                VIMP = 0.0_wp
              ELSEIF ((NS.GE.0.1_wp).AND.(NS.LE.10._wp)) THEN
                VIMP = (0.3272_wp + 0.4907_wp*WIMP - 0.09452_wp*WIMP**2._wp &
                       -0.1906_wp*WIMP**3._wp + 0.07105_wp*WIMP**4._wp)*VT
              ELSEIF (NS.GT.10._wp) THEN
                VIMP = 0.61_wp*VT
              ENDIF
            ELSEIF ((RE.GT.20._wp).AND.(RE.LE.65._wp)) THEN
              IF (NS.LT.0.1_wp) THEN
                VIMP = 0.0_wp
              ELSEIF ((NS.GE.0.1_wp).AND.(NS.LE.10._wp)) THEN
                VIMP = (0.2927_wp + 0.5085_wp*WIMP - 0.03453_wp*WIMP**2._wp &
                       -0.2184_wp*WIMP**3._wp + 0.03595_wp*WIMP**4._wp)*VT
              ELSEIF (NS.GT.10) THEN
                VIMP = 0.59_wp*VT
              ENDIF
            ELSEIF (RE.LE.20._wp) THEN
              IF (NS.LT.0.4_wp) THEN
                VIMP = 0.0_wp
              ELSEIF ((NS.GE.0.4_wp).AND.(NS.LE.10._wp)) THEN
                VIMP = (0.1701_wp + 0.7246_wp*WIMP + 0.2257_wp*WIMP**2._wp &
                       -1.13_wp*WIMP**3._wp + 0.5756_wp*WIMP**4._wp)*VT
              ELSEIF (NS.GT.10._wp) THEN
                VIMP = 0.57_wp*VT
              ENDIF
            ENDIF
              
            !RIME LAYER DENSITY, HEYMSFIELD AND PFLAUM 1985 FORM
            TSCELSIUS = TS - 273.155_wp 
            AFACTOR = -DC*VIMP/MIN(TSCELSIUS, -1.E-6_wp)
            IF ((TSCELSIUS.LE.-5._wp).OR.(AFACTOR.GE.-1.60_wp)) THEN
              DENSELW = 0.30_wp*(AFACTOR)**0.44_wp
            ELSEIF (TSCELSIUS.GT.-5._wp) THEN
              DENSELW = EXP(-0.03115_wp - 1.7030_wp*AFACTOR + &
                             0.9116_wp*AFACTOR**2._wp - 0.1224_wp*AFACTOR**3._wp)
            ENDIF

            DENSELW = DENSELW * 1000._wp !KG M-3
            !BOUND POSSIBLE DENSITIES
            IF (DENSELW.LT.100._wp) DENSELW=100._wp
            IF (DENSELW.GT.900._wp) DENSELW=900._wp
          ELSE
            !In case condition ((DGMW.GT.0._wp).OR.(DGMV.GT.0._wp)) is not met, 
            !DENSELW would be un-initialized.
            !We set it in this case to 500 which the middle value of the 
            !allowed range which is 100-900. "
            DENSELW = 500._wp
          ENDIF

          IF (DGMI.GT.0._wp) THEN
            !Ice collection main source of growth, so set new density layer
            DENSELI = 700._wp
          ENDIF
          
          !All liquid water contributes to growth, none is soaked into center.
          DGMW_NOSOAK = DGMW  !All liquid water contributes to growth,
                              ! none of it is soaked into center.

      ELSE !WET GROWTH
          !Collected liquid water can soak into the stone before freezing,
          ! increasing mass and density but leaving volume constant.
          !Volume of current drop, before growth 
          VOL1 = GM/DENSE
          !Difference b/w mass of stone if density is 900 kg/m3, and
          ! current mass
          SOAK = 900._wp*VOL1 - GM
          !Liquid mass available
          SOAKM = DGMW
          !Soak up as much liquid as we can, up to a density of 900 kg/m3
          IF (SOAKM.GT.SOAK) SOAKM=SOAK
          GM = GM+SOAKM  !Mass of current drop, plus soaking
          !New density of current drop, including soaking but before growth
          DENSE = GM/VOL1 
          !Mass increment of liquid water growth that doesn't
          ! include the liquid water we just soaked into the stone.
          DGMW_NOSOAK = DGMW - SOAKM
          
          !Whatever growth does occur has high density
          DENSELW = 900._wp  !KG M-3
          DENSELI = 900._wp
         
      ENDIF

      !!!VOLUME OF NEW LAYER
      IF (DGMI.LE.0._wp) THEN
         VOLL = (DGMW_NOSOAK+DGMV) / DENSELW
      ELSE IF (DGMW.LE.0._wp) THEN
         VOLL = (DGMI) / DENSELI
      ELSE
           VOLL = (DGMI) / DENSELI + (DGMW_NOSOAK+DGMV) / DENSELW
      ENDIF

      !!!NEW TOTAL VOLUME, DENSITY, DIAMETER
      VOLT = VOLL + GM/DENSE
      DENSE = (GM+DGMI+DGMV+DGMW_NOSOAK) / VOLT
      GM = GM+DGMI+DGMW_NOSOAK+DGMV
      D = ( (6._wp*GM) / (PI*DENSE) )**0.33333333_wp 

  END SUBROUTINE MASSAGR



  SUBROUTINE HEATBUD(TS,TSm1,TSm2,FW,TC,DELRW,D,DENSA,GM1,GM,DGM,DGMW,     &
                     DGMI,DI,RE,AE,SEKDEL,ITYPE)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! CALCULATE HAILSTONE'S HEAT BUDGET 
  !!! See Rasmussen and Heymsfield 1987; JAS
  !!! Original Hailcast's variable units
  !!! TS - Celsius
  !!! FW - unitless, between 0 and 1
  !!! TC - Celsius
  !!! VT - m/s
  !!! D  - m
  !!! DELRW - g/cm3 (per comment)
  !!! DENSA - g/cm3 (per comment)
  !!! GM1, DMG, DGMW, DGMV, DGMI, GMW, GMI - should all be kg
  !!! DI - cm2 / sec
  !!! P  - hPa
  !!! Original HAILCAST HEATBUD subroutine uses c-g-s units, so do some conversions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !$ACC ROUTINE SEQ
      IMPLICIT NONE
      REAL(wp), INTENT(IN)   :: D
      REAL(wp), INTENT(IN)   :: TSm1,TSm2,TC,DELRW,DENSA,GM1,GM,DGM,DGMW,  &
                     DGMI,DI,RE,AE,SEKDEL
      REAL(wp), INTENT(INOUT)  :: TS, FW
      INTEGER,  INTENT(INOUT):: ITYPE
      
      REAL(wp) ALF, ALV, ALS, CI, CW, AK
      REAL(wp) H, AH, TCC, TSC, DELRWC, DENSAC, TDIFF
      REAL(wp) DMLT
      REAL(wp) TSCm1, TSCm2

      ALF = 79.7_wp
      ALV = 597.3_wp
      ALS = 677.0_wp
      CI = 0.5_wp
      CW = 1._wp
      
      !Convert values to non-SI units here
      TSC = TS - 273.155_wp
      TSCm1 = TSm1 - 273.155_wp
      TSCm2 = TSm2 - 273.155_wp
      TCC = TC - 273.155_wp
      DELRWC = DELRW * (1.E3_wp) * (1.E-6_wp)
      DENSAC = DENSA * (1.E3_wp) * (1.E-6_wp)
      !DI still in cm2/sec

      !!!  CALCULATE THE CONSTANTS 
      AK=(5.8_wp+0.0184_wp*TCC)*1.E-5_wp  !thermal conductivity - cal/(cm*sec*K)

      H=(0.71_wp)**(0.333333333_wp)*(RE**0.50_wp) !ventilation coefficient heat (fh)

      !!!   SELECT APPROPRIATE VALUES OF AH AND AE ACCORDING TO RE
      IF(RE.LT.6000.0_wp)THEN
         AH=0.78_wp+0.308_wp*H
      ELSEIF(RE.GE.6000.0_wp.AND.RE.LT.20000.0_wp)THEN
         AH=0.76_wp*H
      ELSEIF(RE.GE.20000.0_wp) THEN
         AH=(0.57_wp+9.0E-6_wp*RE)*H
      ENDIF

      !!!  FOR DRY GROWTH FW=0, CALCULATE NEW TS, ITYPE=1 
      !!!  FOR WET GROWTH TS=0, CALCULATE NEW FW, ITYPE=2

      IF(ITYPE.EQ.1) THEN
      !!!  DRY GROWTH; CALC NEW TEMP OF THE STONE 

         TSC=0.6_wp*(TSC-TSC*DGM/GM1+SEKDEL/(GM1*CI)*       &
            (2._wp*PI*D*(AH*AK*(TCC-TSC)-AE*ALS*DI*DELRWC)+ &
            DGMW/SEKDEL*(ALF+CW*TCC)+DGMI/SEKDEL*CI*TCC)) + &
            0.2_wp*TSCm1 + 0.2_wp*TSCm2
         
         TS = TSC+273.155_wp
         IF (TS.GE.273.155_wp) THEN 
            TS=273.155_wp
         ENDIF
         TDIFF = ABS(TS-273.155_wp)         
         IF (TDIFF.LE.1.E-6_wp) ITYPE=2  !NOW IN WET GROWTH
     
      ELSE IF (ITYPE.EQ.2) THEN
      !!!  WET GROWTH; CALC NEW FW          
         
         IF (TCC.LT.0._wp) THEN
            !Original Hailcast algorithm
            FW=FW-FW*DGM/GM1+SEKDEL/(GM1*ALF)*               &
                (2._wp*PI*D*(AH*AK*TCC-AE*ALV*DI*DELRWC)+    &
                DGMW/SEKDEL*(ALF+CW*TCC)+DGMI/SEKDEL*CI*TCC)
         ELSE
            !Calculate decrease in ice mass due to melting
            DMLT = (2._wp*PI*D*AH*AK*TCC + 2._wp*PI*D*AE*ALV*DI*DELRWC + &
                    DGMW/SEKDEL*CW*TCC) / ALF
            FW = (FW*GM + DMLT) / GM
         ENDIF
         
         IF(FW.GT.1._wp)FW=1._wp
         IF(FW.LT.0._wp)FW=0._wp

         !IF ALL OUR ACCRETED WATER WAS FROZEN, WE ARE BACK IN DRY GROWTH
         IF(FW.LE.1.E-6_wp) THEN
            ITYPE=1  
         ENDIF
         
      ENDIF

  END SUBROUTINE HEATBUD


  
  SUBROUTINE BREAKUP(DENSE,D,GM,FW,CRIT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  TEST IF AMOUNT OF WATER ON SURFACE EXCEEDS CRTICAL LIMIT- 
  !!!  IF SO INVOKE SHEDDING SCHEME 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !$ACC ROUTINE SEQ
      IMPLICIT NONE
      REAL(wp), INTENT(OUT)   :: D
      REAL(wp), INTENT(INOUT) :: FW, GM
      REAL(wp), INTENT(IN)    :: DENSE, CRIT
      !local variables
      REAL(wp) WATER, WAT

      WATER=FW*GM  !KG

      WAT=WATER-CRIT
      GM=GM-WAT
      FW=(CRIT)/GM
    
      IF(FW.GT.1.0_wp) FW=1.0_wp
      IF(FW.LT.0.0_wp) FW=0.0_wp

      ! RECALCULATE DIAMETER AFTER SHEDDING 
      ! Assume density remains the same
      D=(6._wp*GM/(PI*DENSE))**(0.333333333_wp)
  END SUBROUTINE BREAKUP
  
  
  SUBROUTINE MELT(D,TLAYER,PLAYER,RLAYER,LDEPTH,VT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  This is a spherical hail melting estimate based on the Goyer 
  !!!  et al. (1969) eqn (3).  The depth of the warm layer, estimated 
  !!!  terminal velocity, and mean temperature of the warm layer are 
  !!!  used.  DRB.  11/17/2003.
  !!!
  !!!  INPUT:  TLAYER   mean sub-cloud layer temperature (K)
  !!!          PLAYER   mean sub-cloud layer pressure (Pa)
  !!!          RLAYER   mean sub-cloud layer mixing ratio (kg/kg)
  !!!          VT       terminal velocity of stone (m/s)
  !!!  OUTPUT: D        diameter (m)
  !!!          
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !$ACC ROUTINE SEQ
      IMPLICIT NONE

      REAL(wp), INTENT(INOUT) :: D
      REAL(wp), INTENT(IN)    :: TLAYER, PLAYER, RLAYER, LDEPTH, VT
      REAL(wp) eenv, delta, ewet, de, der, wetold, wetbulb, wetbulbk
      REAL(wp) tdclayer, tclayer, eps, b, hplayer
      REAL(wp) a
      REAL(wp) ka, t0, dv, rhoice, &
           tres, re, delt, esenv, rhosenv, essfc, rhosfc, dsig, &
           dmdt, mass, massorg, newmass, gamma, r, rho
      INTEGER wcnt
      REAL(wp), PARAMETER       ::           &
        b1       =   610.78_wp              ,&
        b2w      =    17.2693882_wp
      
      !Convert temp to Celsius, calculate dewpoint in celsius
      eps = 0.622_wp
      tclayer = TLAYER - 273.155_wp
      a = 2.53E11_wp
      b = 5.42E3_wp
      tdclayer = b / LOG(a*eps / (rlayer*player))
      hplayer = player / 100._wp
      
      !Calculate partial vapor pressure
      eenv = (player*rlayer) / (rlayer+eps)
      eenv = eenv / 100._wp  !convert to mb
      
      !Estimate wet bulb temperature (C)
      gamma = 6.6E-4_wp*player
      delta = (4098.0_wp*eenv)/((tdclayer+237.3_wp)*(tdclayer+237.3_wp))
      wetbulb = ((gamma*tclayer)+(delta*tdclayer))/(gamma+delta)
      
      !Iterate to get exact wet bulb
      wcnt = 0
      DO WHILE (wcnt .lt. 11)
        ewet = b1/100.0_wp*(exp((b2w*wetbulb)/(237.3_wp + wetbulb))) 
        de = (0.0006355_wp*hplayer*(tclayer-wetbulb))-(ewet-eenv)
        der= (ewet*(.0091379024_wp - (6106.396_wp/(273.155_wp+wetbulb)**2._wp))) &
             - (0.0006355_wp*hplayer)
        wetold = wetbulb
        wetbulb = wetbulb - de/der
        wcnt = wcnt + 1
        IF ((abs(wetbulb-wetold)/wetbulb.gt.0.0001_wp)) THEN
           wcnt = 12 
        ENDIF
      ENDDO
      
      wetbulbk = wetbulb + 273.155_wp  !convert to K
      ka = .02_wp ! thermal conductivity of air
      t0 = 273.155_wp ! temp of ice/water melting interface
      dv = 0.25e-4_wp ! diffusivity of water vapor (m2/s)
      rhoice = 917.0_wp ! density of ice (kg/m**3)
      r = D/2._wp ! radius of stone (m)
      
      !Compute residence time in warm layer
      tres = LDEPTH / VT
        
      !Calculate dmdt based on eqn (3) of Goyer et al. (1969)
      !Reynolds number...from pg 317 of Atmo Physics (Salby 1996)
      !Just use the density of air at 850 mb...close enough.
      rho = 85000._wp/(rd*TLAYER)
      re = rho*r*VT*.01_wp/1.7e-5_wp
      
      !Temperature difference between environment and hailstone surface
      delt = wetbulb !- 0.0 !assume stone surface is at 0C
                            !wetbulb is in Celsius

      !Difference in vapor density of air stream and equil vapor
      !density at the sfc of the hailstone
      esenv = b1*(exp((b2w*wetbulb)/  &
               (237.3_wp + wetbulb))) ! es environment in Pa
      rhosenv = esenv/(rv*wetbulbk)
      essfc = b1*(exp((b2w*(t0-273.155_wp))/  &
               (237.3_wp + (t0-273.155_wp)))) ! es environment in Pa
      rhosfc = essfc/(rv*t0)
      dsig = rhosenv - rhosfc

      !Calculate new mass growth
      dmdt = (-1.7_wp*pi*r*(re**0.5_wp)/alf)*((ka*delt)+((alv-alf)*dv*dsig))
      IF (dmdt.gt.0.0_wp) dmdt = 0._wp
      mass = dmdt*tres
      
      !Find the new hailstone diameter
      massorg = 1.33333333_wp*pi*r*r*r*rhoice
      newmass = massorg + mass
      if (newmass.lt.0.0_wp) newmass = 0.0_wp

      D = 2._wp*(0.75_wp*newmass/(pi*rhoice))**0.333333333_wp
  END SUBROUTINE MELT

END MODULE mo_diag_hailcast
