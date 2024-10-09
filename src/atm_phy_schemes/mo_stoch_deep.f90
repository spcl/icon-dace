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

! Stochastic differential equation (SDE) deep convection subroutine.
!
! Description:
! This routine is not a stand-alone convection parameterization, but is
! used in conjunction with a conventional mass flux parameterization
! scheme, currently the Tiedtke-Bechtold scheme. It replaces the cloud
! base mass flux calculated with the conventional scheme with a
! stochastically perturbed mass flux value for each grid cell.
! The SDE version of the scheme evolves two prognostic variables
! describing the state of the cloud ensemble (grid box mass flux and cloud
! numbers).
! Required input for the construction of the mass flux distribution is
! the "bulk" or first-guess mass flux calculated by the conventional
! convection scheme (applied to profiles averaged over a representative
! neighbourhood).
!
! A full description of the scheme can be found in the following
! publications:
! Plant and Craig (2008)
! Machulskaya and Seifert (2019)
!
! Compile with FFLAG -fall-intrinsics
!
!-----------------------------------------------------------------------------

MODULE mo_stoch_deep

  USE netcdf
#ifdef HAVE_ACM_LICENSE  
  USE random_rewrite,        ONLY: random_Poisson
#endif
  USE mo_kind,               ONLY: wp, i4
  USE mo_cuparameters, ONLY :                                    &
       & deep_k_wei, deep_alpha_mf, deep_beta_mf, deep_mean_mf, deep_mean_tau
  USE cloud_random_numbers,  ONLY: rng_type, IRngMinstdVector, IRngNative

  IMPLICIT NONE

  PUBLIC :: deep_stoch_sde

  CONTAINS

  !! Main subroutine calculating mass flux (MF) distribution, birth and death rates, determining
  !! the number of newborn/dying clouds, assigning each cloud a mass flux and incrementing the
  !! prognostic variables for grid box mass flux and cloud numbers.  
  SUBROUTINE deep_stoch_sde(i_startidx,i_endidx,klon,ptsphy,&
                           & mfb, mfp, ktype,clnum, &
                           & pclnum_d, pclmf_d, &
                           & lseed,extra_3d,luse3d,cell_area)

  !Interface variables
  INTEGER(i4),INTENT(IN)    :: i_startidx         ! block index start            
  INTEGER(i4),INTENT(IN)    :: i_endidx           ! block index end
  INTEGER(i4),INTENT(IN)    :: klon               ! number of grid points in block
  REAL(wp)   ,INTENT(IN)    :: ptsphy             ! length of physics time step (s)
  REAL(wp)   ,INTENT(IN)    :: mfb(klon)          ! bulk/first-guess mass flux derived from T-B scheme (normalised by grid cell area)
  REAL(wp)   ,INTENT(OUT)   :: mfp(klon)          ! stochastic perturbed mass flux (diagnostic)
  INTEGER(i4),INTENT(IN)    :: ktype(klon)        ! convection type (0:off, 1:deep, 2:shallow, 3:mid-level)
  REAL(wp)   ,INTENT(inout) :: pclnum_d(klon)     ! prognostic cloud number
  REAL(wp)   ,INTENT(inout) :: pclmf_d(klon)      ! prognostic mass flux
  REAL(wp)   ,INTENT(OUT)   :: clnum(klon)        ! diagnostic number of clouds per grid cell
  INTEGER(i4),INTENT(in)    :: lseed(klon)        ! seed for random number generator
  REAL(wp)   ,INTENT(INOUT), DIMENSION(:,:,:), OPTIONAL :: extra_3d! 3D extra diagnostics
  LOGICAL    ,INTENT(IN)    :: luse3d             ! logical to switch on output of extra 3D diagnostics
  REAL(wp)   ,INTENT(in)    :: cell_area(klon)    ! grid cell area

  !Local variables
  REAL(wp) :: M(klon)                             ! total bulk mass flux in each cell (not normalised by grid cell area)
  REAL(wp) :: mf_new(klon)                        ! increase of MF due to newly born clouds
  REAL(wp) :: mf_die(klon)                        ! decrease of MF due to dying clouds
  REAL(wp) :: dxy                                 ! grid cell area
  REAL(wp) :: rdxy                                ! inverse of grid cell area
  REAL(wp) :: eps                                 ! equivalent of 0.9 clouds in units of "clouds per grid cell area"
  REAL(wp) :: excess                              ! number of clouds exceeding upper limit (5000)
  REAL(wp) :: birth_rate(klon)                    ! birth rate
  REAL(wp) :: death_rate(klon)                    ! death rate
  REAL(wp) :: mean_m(klon)                        ! avg mass flux of current cloud ensemble
  REAL(wp) :: pclmf_now_d(klon),pclnum_now_d(klon)
  REAL(wp) :: z0   

  INTEGER(i4) :: i, k, kstart,kstop,idx,j           ! (loop) indices
  INTEGER(i4) :: rn_poisson1, rn_poisson2, rn1    ! cloud numbers drawn from Poisson distribution
  INTEGER(i4) :: iseed                            ! local copy of seed for random number generator
  INTEGER(i4) :: streammax                        ! length of random number sequence to be generated at once
  REAL(wp), DIMENSION(5120) :: rn_u               ! sequence of random numbers generated by call to generator
  ! The max number of random numbers that can be retrieved from the generator is hardcoded, currently to 5120
  
  ! Type for random number generator and stream for producing random
  ! numbers
  TYPE(rng_type) :: random_number_generator

  REAL(wp), PARAMETER :: pi = 3.14159265358

  ! Loop over grid points
  DO i = i_startidx, i_endidx

     pclmf_now_d(i)=pclmf_d(i)
     pclnum_now_d(i)=pclnum_d(i)
     
     !Initialize
     birth_rate(i)=0._wp
     death_rate(i)=0._wp
     mean_m(i)    =0._wp

     ! specify grid cell area
     dxy     = cell_area(i)
     ! inverse of grid cell area
     rdxy    = 1._wp/dxy
     ! equivalent of .9 clouds in units of clouds per unit area
     ! needed for checks of cloud numbers > zero
     eps     = .9*rdxy

     ! Only consider grid cells with active deep convection
     IF (ktype(i) == 1) THEN
       ! Bulk mass flux from traditional T-B closure, convert into units of kg/s
       ! by multiplying with grid box area
       M(i)  = mfb(i)*dxy ! in kg/s; rho =1 kg/m3

       ! Plant and Craig 2008 assume a fixed mean mass flux per cloud
       ! deep_mean_mf=2.0E+07_wp
       birth_rate(i) = M(i)/(deep_mean_mf*deep_mean_tau)
       IF(pclnum_now_d(i) < eps .AND. pclmf_now_d(i) == 0._wp) THEN
         ! If no clouds exist yet in cloud ensemble, we must make an assumption
         ! about the cloud number and their average mass flux in order to calculate
         ! a death rate. Here, we assume that the ensemble mass flux is the bulk
         ! mass flux mfb, and the cloud number the bulk mass flux divided by the
         ! fixed mean mass flux per cloud from Plant and Craig.
         ! The values in the prognostic variables here are overwritten later, so only
         ! used for the calculation of death_rate and mean_m.
         pclnum_d(i) = mfb(i)/deep_mean_mf
         pclmf_d (i) = mfb(i)
       ELSE
         pclnum_d(i) = pclnum_now_d(i)
         pclmf_d (i) = pclmf_now_d (i)
       END IF

       ! Calculate actual mean mass flux per cloud currently existing
       ! in the ensemble.
       mean_m(i) = pclmf_d(i)*dxy/MAX(pclnum_d(i)*dxy, 1.E-10_wp)
       ! Death rate see Machulskaya and Seifert 2019, Eqn 3
       ! Assuming exponentially distributed cloud life times with
       ! mean cloud life time deep_mean_tau=45min
       death_rate(i) = pclnum_d(i)*dxy/deep_mean_tau
    ENDIF
  ENDDO

  ! RANDOM NUMBER GENERATION AND DETERMINATION OF NUMBER OF NEWBORN CLOUDS BY DRAW FROM
  ! POISSON DISTRIBUTION
  ! Separate loop for random_Poisson call because it doesn't vectorise on NEC
  
  !Loop over grid points
  DO i = i_startidx, i_endidx
     
    ! The random seed for each grid point is calculated based on lat, lon, forecast time and ensemble number,
    ! such that the seed is reproducable for subsequent runs. NOTE: the seed is not updated in this routine!
    iseed=lseed(i)

    ! Get sufficient random numbers for two calls to the random_Poisson routine. This routine iterates, and
    ! 50 random numbers per call is a very conservative estimate of the random numbers needed.
    streammax=100
    ! Initialise the random number generator
    CALL random_number_generator%initialize(IRngMinstdVector, iseed=iseed, &
         &                                  nmaxstreams=streammax)
    ! Retrieve streammax new random numbers to use in draw from Poisson distributions
    CALL random_number_generator%uniform_distribution(rn_u(1:streammax))

#ifdef HAVE_ACM_LICENSE
    ! Determine number of newborn clouds by drawing from Poisson distribution with given birth rate,
    ! and add mass flux to grid cell total mass flux
    idx=1
    rn_poisson1 = random_Poisson(birth_rate(i)*ptsphy,idx,streammax,rn_u(1:streammax))
    ! Safety check in case random_Poisson routine does not converge and runs out of random numbers to use
    IF (idx > streammax) THEN
       WRITE(6,*) 'idx for deep cld birth out of range',idx,birth_rate(i)*ptsphy, &
            & rn_u(1),rn_u(idx),MINVAL(rn_u(1:streammax)),MAXVAL(rn_u(1:streammax))
       idx=2
    endif
    
    ! Determine number of dying clouds by drawing from Poisson distribution with given death rate,
    ! and subtract mass flux from grid cell total mass flux
    idx=50! Start at index 50, to make sure to use "fresh" random numbers
    rn_poisson2 = random_Poisson(death_rate(i)*ptsphy,idx,streammax,rn_u(1:streammax))
    IF (idx > streammax) THEN
       write(6,*) 'idx for deep cloud death out of range',idx,death_rate(i)*ptsphy, &
            & rn_u(1),rn_u(idx),MINVAL(rn_u(1:streammax)),MAXVAL(rn_u(1:streammax))
       idx=3
    ENDIF
    idx=100
#else
    ! Determine number of newborn active clouds by sampling a normal distribution
    ! (instead of Poisson) using the Box-Muller method, with given birth rate.
    idx=1 ! This index keeps track of which random numbers out of the 200 have already been used
    z0 = SQRT(-2._wp * LOG(rn_u(idx))) * COS(pi * rn_u(idx+25))
    rn_poisson1 = MAX(0, INT(z0 * SQRT(birth_rate(i)*ptsphy) + birth_rate(i)*ptsphy))
    
    ! Determine number of dying active clouds by sampling a normal distribution
    ! (instead of Poisson) using the Box-Muller method, with given death rate.
    idx=50 
    z0 = SQRT(-2._wp * LOG(rn_u(idx))) * COS(pi * rn_u(idx+25))
    rn_poisson2 = MAX(0, INT(z0 * SQRT(death_rate(i)*ptsphy) + death_rate(i)*ptsphy))
    
    idx=100
#endif

    ! Initialise diagnostic fields that keep track of how much mass flux is added by newborn clouds,
    ! and how much is removed by dying clouds
    mf_new(i)=0._wp
    mf_die(i)=0._wp

    ! specify grid cell area (must be repeated because it's grid point dependent, and we're in a new loop)
    dxy     = cell_area(i)
    ! inverse of grid cell area
    rdxy    = 1._wp/dxy
    ! equivalent of .9 clouds in units of clouds per unit area
    ! needed for checks of cloud numbers > zero
    eps     = .9*rdxy

    ! Only consider grid cells with active deep convection
    IF (ktype(i) == 1) THEN

      ! Update stochastic mass flux and cloud numbers from last time step
      ! The placeholder-values from the previous loop are overwritten here
      pclmf_d(i)=pclmf_now_d(i)
      pclnum_d(i)=pclnum_now_d(i)

      ! GENERATE RANDOM NUMBERS, AND DRAW EACH NEW CLOUD'S MASS FLUX FROM
      ! THE DISTRIBUTION

      ! Get a random number for each newly born/dying cloud, plus 100 extra
      ! so we can throw out the first 100 that were already used for the Poisson draws.
      
      streammax=rn_poisson1+rn_poisson2+idx
      ! Safety check - we can draw only 5120 random numbers at a time. This
      ! limit is hard-coded into the random number generator.
      ! In case this code is run at very coarse resolutions, this will make sure
      ! the number is capped at the limit. 
      ! Reduce cloud numbers equally in both categories of born/dying clouds.
      IF (streammax > 5120) THEN
        excess=streammax/5120._wp-1._wp
        rn_poisson1=rn_poisson1-CEILING(excess*(rn_poisson1+idx*.5))
        rn_poisson2=rn_poisson2-CEILING(excess*(rn_poisson2+idx*.5))
        streammax=rn_poisson1+rn_poisson2+idx
        write(6,*) 'Warning! The limit of 5000 clouds in deep stoch convection was exceeded!'
     ENDIF!

      ! Initialize and retrieve only as many random numbers as will be needed - speeds up code
      CALL random_number_generator%initialize(IRngMinstdVector, iseed=iseed, &
           &                                  nmaxstreams=streammax)
      CALL random_number_generator%uniform_distribution(rn_u(1:streammax))
      
      kstart=idx+1! Start using first random number at 101
      kstop=idx+rn_poisson1
      DO k = kstart, kstop
        ! Randomly draw mass flux from exponential distribution. 
        ! Note: for the birth rate, a fixed value for the average cloud
        ! mass flux is used.
        pclmf_d(i) = pclmf_d(i) - LOG(1.-rn_u(k))*deep_mean_mf*rdxy
        ! Keep track of active mass flux added in diagnostic
        mf_new(i)  = mf_new(i)  - LOG(1.-rn_u(k))*deep_mean_mf*rdxy
      END DO
     
      kstart=kstop+1
      kstop=kstop+rn_poisson2
      DO k = kstart,kstop
        ! Randomly draw mass flux from exponential distribution. Note:
        ! for the death rate, the actual mean mass flux per cloud
        ! of the existing cloud ensemble is used.
        pclmf_d(i) = pclmf_d(i) + LOG(1.-rn_u(k))*mean_m(i)*rdxy
        mf_die(i)  = mf_die(i)  - LOG(1.-rn_u(k))*mean_m(i)*rdxy
      END DO

      ! update prognostic cloud numbers, and check for negative values
      ! Note: prognostic cloud numbers have units of "clouds per grid box area"
      pclnum_d(i) = pclnum_d(i) + float(rn_poisson1 - rn_poisson2)*rdxy
      pclnum_d(i) = MAX(pclnum_d(i),0._wp)
      pclmf_d (i) = MAX(pclmf_d (i),0._wp)

      ! Enforce consistency between mass flux and cloud numbers
      ! Check if cloud numbers have dropped below 0.9, and
      ! remove any remaining MF if that is the case
      IF (pclnum_d(i) <= eps) THEN
         pclmf_d (i) = 0._wp
      ELSE
         pclmf_d (i) = MAX(pclmf_d (i),0._wp)
      ENDIF

      ! update diagnostic variables for MF and cloud numbers
      mfp(i)=pclmf_d(i)
      ! Diagnostic uses absolute cloud number, not number per unit area
      clnum(i)=(pclnum_d(i))*dxy !Diagnostic uses absolute cloud number, not number per unit area
      ! get rid of very small values 1e-16 left by unit conversion
      IF (clnum(i) < .1_wp) clnum(i)= 0.0_wp 

      IF (luse3d) THEN 
        ! write diagnostics
        extra_3d(i,1,3)=birth_rate(i)*ptsphy     ! birth rate
        extra_3d(i,2,3)=pclnum_d(i)*dxy          ! current number of clouds
        extra_3d(i,3,3)=float(rn_poisson1)       ! newborn clouds
        extra_3d(i,4,3)=float(rn_poisson2)       ! dying clouds
        extra_3d(i,5,3)=pclmf_d(i)               ! mass flux
        extra_3d(i,6,3)=mf_new(i)                ! mf of newborn clouds
        extra_3d(i,7,3)=mf_die(i)                ! mf of dying clouds
        extra_3d(i,8,3)=mean_m(i)                ! avg mass flux of cloud ensemble distribution
        extra_3d(i,9,3)=death_rate(i)*ptsphy     ! death rate
      ENDIF
    ELSE ! not ktype=1

      ! If grid cell has no active deep convection, continue decaying number and massflux according to death rate.
      ! While the mass flux calculated for this time step is not actively used by the convection scheme,
      ! the state of the decaying cloud ensemble may be picked up again if convection switches back on during
      ! subsequent time steps.

      ! We only need to generate random numbers for dying clouds +100 used for prior Poisson draw 
      streammax=rn_poisson2+idx
      ! Safety check again that random numbers to be drawn don't exceed the hard-coded max limit,
      ! and reduce numbers if necessary.
      IF (streammax > 5120) THEN
        rn_poisson2=5120
        streammax=rn_poisson2
      ENDIF
      ! Initialise random number generator
      CALL random_number_generator%initialize(IRngMinstdVector, iseed=iseed, &
           &                                  nmaxstreams=streammax)
      ! Random numbers to use in decay of dying clouds
      CALL random_number_generator%uniform_distribution(rn_u(1:streammax))
      
      kstart=idx+1
      kstop=idx+rn_poisson2
      DO k = kstart,kstop
        ! Randomly draw mass flux from exponential distribution. Note:
        ! for the death rate, the actual mean mass flux per cloud
        ! of the existing cloud ensemble is used.
        pclmf_d(i) = pclmf_d(i) + LOG(1.-rn_u(k))*mean_m(i)*rdxy
        mf_die(i)  = mf_die(i)  - LOG(1.-rn_u(k))*mean_m(i)*rdxy
      END DO
      
      ! Update prognostics cloud numbers
      pclnum_d(i) = pclnum_now_d(i) - float(rn_poisson2)*rdxy
      pclnum_d(i) = MAX(pclnum_d(i),0._wp)

      ! Enforce consistency between MF and cloud numbers
      IF (pclnum_d(i) <= eps) THEN
         pclmf_d (i) = 0._wp
      ELSE
         pclmf_d (i) = MAX(pclmf_d (i),0._wp)
      ENDIF

      ! update diagnostic variables for MF and cloud numbers
      ! NOTE: This means a perturbed MF and cloud number will be put out even
      ! for grid points that are not deep, but may be e.g. shallow at this particular
      ! time step. The perturbed mass flux is not used in cumaster if ktype=2, but
      ! it may be confusing to see perturbed MF output at grid points that have ktype=2.
      ! This update of the diagnostics should possibly be commented out.
      !mfp(i)=pclmf_d(i)
      ! Diagnostic uses absolute cloud number, not number per unit area
      !clnum(i)=(pclnum_d(i))*dxy !Diagnostic uses absolute cloud number, not number per unit area
      ! get rid of very small values 1e-16 left by unit conversion
      !IF (clnum(i) < .1_wp) clnum(i)= 0.0_wp

      ! Side note:
      ! The stochastic deep and shallow schemes are not meant to be run together.
      ! If they are, then updating the diagnostic variables here may overwrite
      ! entries from the shallow stochastic scheme.
      
      IF (luse3d) then 
         ! write diagnostics
         extra_3d(i,1,3)=0._wp                  ! birth rate
         extra_3d(i,2,3)=pclnum_d(i)*dxy            ! current number of clouds
         extra_3d(i,3,3)=0._wp                  ! newborn clouds
         extra_3d(i,4,3)=float(rn_poisson2)       ! dying clouds
         extra_3d(i,5,3)=pclmf_d(i)               ! mass flux
         extra_3d(i,6,3)=0._wp                  ! mf of newborn clouds
         extra_3d(i,7,3)=mf_die(i)                ! mf of dying clouds
         extra_3d(i,8,3)=mean_m(i)                ! avg mass flux of cloud ensemble distribution
         extra_3d(i,9,3)=death_rate(i)*ptsphy     ! death rate
      ENDIF
    ENDIF
 ENDDO
END SUBROUTINE deep_stoch_sde

END MODULE mo_stoch_deep
