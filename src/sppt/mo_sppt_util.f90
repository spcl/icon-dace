!
! Utility subroutines used for SPPT
! (Stochastic Perturbation of Physics Tendencies)
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

MODULE mo_sppt_util 

  USE mtime,                      ONLY: datetime, OPERATOR(+), &
    &                                   newDatetime, deallocateDatetime,  &
    &                                   datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_util_mtime,              ONLY: is_event_active
  USE mo_kind,                    ONLY: wp
  USE mo_parallel_config,         ONLY: proc0_offloading
  USE mo_model_domain,            ONLY: t_patch
  USE mo_sppt_config,             ONLY: t_sppt_config
  USE mo_sppt_core,               ONLY: time_interpol_rn
  USE mo_loopindices,             ONLY: get_indices_c
  USE mo_impl_constants,          ONLY: min_rlcell_int
  USE mo_gribout_config,          ONLY: gribout_config
  USE mo_exception,               ONLY: message, message_text, finish
  USE mo_fortran_tools,           ONLY: set_acc_host_or_device, assert_acc_device_only, copy


  IMPLICIT NONE

  PRIVATE


  PUBLIC :: construct_rn
  PUBLIC :: init_rn

  CONTAINS

  !---------------------------------------------------------------------
  ! Initiation of random pattern. Initially two random pattern need 
  ! to be generated to be able to interpolate in time. After initiation 
  ! and after hinc_rn amount of time has passed construct_rn is used.
  !---------------------------------------------------------------------

  SUBROUTINE init_rn (pt_patch, mtime_current, sppt_config, rn_2d_now, rn_2d_new)

    ! Subroutine arguments
    TYPE(t_patch),            INTENT(IN)        :: pt_patch         !< patch variables
    TYPE(datetime),  POINTER, INTENT(IN)        :: mtime_current    !< current_datetime
    TYPE(t_sppt_config),      INTENT(INOUT)     :: sppt_config      !< sppt config state

    REAL(wp),                 INTENT(INOUT)     :: rn_2d_now(:,:)   !< 2d utility field
    REAL(wp),                 INTENT(INOUT)     :: rn_2d_new(:,:)   !< 2d utility field

    ! Local variables
    TYPE(datetime), POINTER                     :: mtime_current_new
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)         :: valid_time_string


    !-----------------
    ! Initial read of random numbers
    !---------------

    ! ... since generate_rn expects a pointer
    mtime_current_new => newDatetime('0001-01-01T00:00:00')

    ! Fill rn_2d_now
    CALL generate_rn(pt_patch, sppt_config, mtime_current, rn_2d_now)

    CALL datetimeToString(mtime_current, valid_time_string)
    WRITE(message_text, '(a,a)') 'rn_2d_now generated for : ', valid_time_string
    CALL message('init_rn', message_text)

    ! mtime_current plus hinc_rn
    mtime_current_new = mtime_current + sppt_config%mtime_hinc_rn


    ! Fill rn_2d_new
    CALL generate_rn(pt_patch, sppt_config, mtime_current_new, rn_2d_new)

    CALL datetimeToString(mtime_current_new, valid_time_string)
    WRITE(message_text, '(a,a)') 'rn_2d_new generated for : ', valid_time_string
    CALL message('init_rn', message_text)

    ! store validity date of rn_2d_new
    sppt_config%validity_date_rn_2d_new = mtime_current_new


    CALL deallocateDatetime(mtime_current_new)

  END SUBROUTINE init_rn


  !---------------------------------------------------------------------
  ! This subroutine construct_rn() called within the time loop
  ! takes care of generating random number fields every time increment 
  ! as defined by hinc_rn.
  !
  ! Event timers are used depending on hinc_rn (namelist switch).
  !
  !---------------------------------------------------------------------

  SUBROUTINE construct_rn (pt_patch, mtime_current, sppt_config, rn_3d, rn_2d_now, rn_2d_new, lacc)

    ! Subroutine arguments
    TYPE(t_patch),            INTENT(IN)        :: pt_patch         !< patch variables
    TYPE(datetime),  POINTER, INTENT(IN)        :: mtime_current    !< current_datetime
    TYPE(t_sppt_config),      INTENT(INOUT)     :: sppt_config      !< sppt config state

    REAL(wp),                 INTENT(INOUT)     :: rn_3d(:,:,:)     !< 3d field of random numbers
    REAL(wp),                 INTENT(INOUT)     :: rn_2d_now(:,:)   !< 2d utility field
    REAL(wp),                 INTENT(INOUT)     :: rn_2d_new(:,:)   !< 2d utility field
    LOGICAL,     OPTIONAL,    INTENT(IN)        :: lacc             !< OpenACC flag

    ! Local
    LOGICAL                                     :: rapa_event_active      !< local logical for events
    TYPE(datetime), POINTER                     :: mtime_current_new
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)         :: valid_time_string

    CALL assert_acc_device_only("construct_rn", lacc)

    ! Check if event is active and assign
    rapa_event_active = is_event_active(sppt_config%read_rapa_Event,  mtime_current, proc0_offloading)


    !-----------------
    ! Read random new random number fields - event triggered
    !----------------

    IF(rapa_event_active) THEN ! event active

      CALL copy(rn_2d_new,rn_2d_now)

      ! ... since generate_rn expects a pointer
      mtime_current_new => newDatetime('0001-01-01T00:00:00')

      mtime_current_new = mtime_current + sppt_config%mtime_hinc_rn

      ! Generate random pattern - and write message to output
      CALL generate_rn(pt_patch, sppt_config, mtime_current_new, rn_2d_new, lacc=.TRUE.)

      CALL datetimeToString(mtime_current_new, valid_time_string)
      WRITE(message_text, '(a,a)') 'rn_2d_new generated for : ', valid_time_string
      CALL message('construct_rn', message_text)

      ! store validity date of rn_2d_new
      sppt_config%validity_date_rn_2d_new = mtime_current_new

      CALL deallocateDatetime(mtime_current_new)

    ENDIF

    ! time interpolation/tapering of random numbers
    !
    CALL time_interpol_rn(pt_patch, sppt_config, rn_2d_now(:,:), rn_2d_new(:,:), &
                          mtime_current, rn_3d(:,:,:), lacc=.TRUE.)

  END SUBROUTINE construct_rn



  !---------------------------------------------------------------------
  !
  ! Random number generator
  !
  ! This function returns a normally distributed pseudo-random number with zero
  ! mean and unit variance.
  !
  ! The algorithm uses the ratio of uniforms method of A.J. Kinderman
  ! and J.F. Monahan augmented with quadratic bounding curves.
  !
  ! Adapted from the following Fortran 77 code
  !   ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
  !   THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
  !   VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
  !
  !---------------------------------------------------------------------

#ifdef HAVE_ACM_LICENSE

  SUBROUTINE random_normal_values(seed, values_range, values)

    ! Subroutine arguments (in/out/inout)
    INTEGER,            INTENT(IN)                           :: seed          ! seed to be used
    REAL(wp),           INTENT(IN)                           :: values_range  ! allowed range of random numbers
    REAL(wp),           DIMENSION(:), INTENT(INOUT)          :: values        ! array of value to be filled with random numbers

    ! Local parameters
    CHARACTER(LEN=*), PARAMETER                :: nm = 'random_normal_values: ' ! Tag

    ! For random number generator
    !     Constants used in the random number generators - note particularly
    !     the values for A0, A and M. The type jpim is 4 bytes and jpib is 8 bytes.
    !     This is because in the operation mod(A*X,M) we don't want A*X to overflow
    !     or go negative.  The problem is avoided in C by the use of
    !     unsigned integers where the appropriate value of M basically
    !     represents overflow. If this operation can be done entirely with
    !     4-byte integers performance would presumably be improved.
    ! ... kind

    INTEGER, PARAMETER            :: jprb = SELECTED_REAL_KIND(12,307)  ! corresponds to double precision
    INTEGER, PARAMETER            :: jpim = SELECTED_INT_KIND(9)        ! at least 4 byte integer
    INTEGER, PARAMETER            :: jpib = SELECTED_INT_KIND(14)       ! at least 8 byte integer

    ! ... integer parameters
    INTEGER(KIND=jpib), PARAMETER :: IMinstdA0    = 16807
    INTEGER(KIND=jpib), PARAMETER :: IMinstdA     = 48271

    ! ... real parameters
    REAL(KIND=jprb), PARAMETER    :: IMinstdM     = 2147483647._jprb
    REAL(KIND=jprb), PARAMETER    :: IMinstdScale = 1.0_jprb / IMinstdM

    ! ... size of random numbers block
    INTEGER(KIND=jpim), PARAMETER              :: nmaxstreams = 512

    ! ... for normal distribution
    REAL(KIND=jprb), PARAMETER                 :: ss = 0.449871_jprb
    REAL(KIND=jprb), PARAMETER                 :: tt = -0.386595_jprb
    REAL(KIND=jprb), PARAMETER                 :: aa = 0.19600_jprb
    REAL(KIND=jprb), PARAMETER                 :: bb = 0.25472_jprb
    REAL(KIND=jprb), PARAMETER                 :: r1 = 0.27597_jprb
    REAL(KIND=jprb), PARAMETER                 :: r2 = 0.449871_jprb

    ! Local variables
    INTEGER(KIND=jpim)                         :: ii, jj
    INTEGER(KIND=jpim)                         :: idx
    REAL(KIND=jprb)                            :: rseed, rnd_init
    REAL(KIND=jprb)                            :: uu, vv, xx, yy, qq
    REAL(KIND=jprb), DIMENSION(nmaxstreams)    :: randnum, istate

    ! -----------------------------------------------------------------------
    ! Initialize a random number generator, by using the MINSTD linear
    ! congruential generator (LCG). "nmaxstreams" indicates that random
    ! numbers will be requested in blocks of this length. The generator
    ! is seeded with "seed".
    !-------------------------------------------------------------------------

    rseed = REAL(ABS(seed), KIND=jprb)

    DO ii = 1, nmaxstreams

      rnd_init = NINT(                                                                                &
                  MOD(rseed * ii * (1._jprb - 0.05_jprb*ii + 0.005_jprb*ii*ii) * IMinstdA0, IMinstdM), &
                                                                                            KIND=jpib)
      istate(ii) = MOD(IMinstdA * rnd_init, IMinstdM)

    END DO

    ! -------------------------------------------------------------------------
    ! For each element of values, draw a random number from a normal distribution
    ! -------------------------------------------------------------------------

    values_loop:          &
    DO ii = 1, SIZE(values)

    ! Retrieve nmaxstreams new random numbers to use in draw from normal distribution
    
      DO jj = 1, nmaxstreams
      istate(jj) = MOD(IMinstdA * istate(jj), IMinstdM)
      randnum(jj) = IMinstdScale * istate(jj)
      ENDDO

    ! --------------------------------------------------------------------------  
    ! Draw from normal distribution:
    !
    ! Generate normal distribution from uniform distribution:
    ! P = (u,v) uniform in rectangle enclosing acceptance region
    ! Occasionally random numbers are exactly zero >
    !    log function lower down can't cope and a new draw is done
    !---------------------------------------------------------------------------

      idx = 1
      DO
        IF ( idx > nmaxstreams - 4 ) EXIT values_loop

        uu = randnum(idx)
        IF ( uu == 0. ) THEN
          idx = idx + 1
          uu = randnum(idx)
        ENDIF
        idx = idx + 1

        vv = randnum(idx)
        IF ( vv == 0. ) THEN
          idx = idx + 1
          vv = randnum(idx)
        ENDIF
        vv = 1.7156_jprb * (vv - 0.5_jprb)
        idx = idx + 1

        ! Evaluate the quadratic form
        xx = uu - ss
        yy = ABS(vv) - tt
        qq = xx * xx + yy * (aa*yy - bb*xx)
        ! Accept P if inside inner ellipse
        IF ( qq < r1 ) EXIT
        ! Reject P if outside outer ellipse
        IF ( qq > r2 ) CYCLE
        ! Reject P if outside acceptance region
        IF ( vv * vv < -4.0_jprb * LOG(uu) * uu * uu ) EXIT
      END DO

      ! Assign value
      values(ii) = REAL( vv / uu )

    END DO values_loop

  ! Limit random number range
  WHERE     (values(:) > values_range)
    values(:)  = values_range
  ELSEWHERE (values(:) < -values_range)
    values(:) = -values_range
  END WHERE


  END SUBROUTINE random_normal_values

#endif  

  !---------------------------------------------------------------------
  ! Bilinear interpolation, c.f.
  !
  ! https://en.wikipedia.org/wiki/Bilinear_interpolation#Repeated_linear_interpolation
  !
  !---------------------------------------------------------------------

  FUNCTION lerp2(x, y, x1, x2, y1, y2, fq11, fq21, fq12, fq22) result(fxy)
 
    REAL(wp), INTENT(IN) :: x, y, x1, x2, y1, y2, fq11, fq12, fq21, fq22
    REAL(wp)             :: fxy

    REAL(wp)             :: w, QdotVec_1, QdotVec_2

    !$ACC ROUTINE SEQ

    w = 1.0_wp/((x2-x1)*(y2-y1))     

    QdotVec_1 = fq11*(y2 - y) + fq12*(y - y1) 
    QdotVec_2 = fq21*(y2 - y) + fq22*(y - y1)

    fxy = w*(QdotVec_1*(x2-x) + QdotVec_2*(x-x1))

  END FUNCTION lerp2


  !---------------------------------------------------------------------
  !
  ! Genration of random number filed
  !
  ! This subroutine returns a field of normally distributed pseudo-random numbers with zero
  !   mean and unit variance generated by random_normal_values(). Fields are generated on a regular 
  !   lat/lon grid with predefined grid spacing (namelist switch dlat_rn and dlon_rn) and interpolated 
  !   onto ICON native grid.
  !
  !---------------------------------------------------------------------

  SUBROUTINE generate_rn(p_patch, sppt_config, mtime_current, rn_2d, lacc)

    TYPE(t_patch),            INTENT(IN )    :: p_patch
    TYPE(t_sppt_config),      INTENT(IN )    :: sppt_config
    TYPE(datetime),  POINTER, INTENT(IN )    :: mtime_current     !< current_datetime 
    REAL(wp),                 INTENT(OUT)    :: rn_2d(:,:)        !< gaussian random number field (triangular grid)
    LOGICAL, OPTIONAL,        intent(in)     :: lacc              !< flag to run on GPU

 
    REAL(wp):: lat_cell, lon_cell

    REAL(wp):: lat_lo, lon_lo, lon_hi, lat_hi 
    REAL(wp):: fq11, fq21, fq12, fq22   

    INTEGER :: coarse_nlat, coarse_nlon
    INTEGER :: jb, jc
    INTEGER :: i_startidx, i_endidx, i_startblk, i_endblk, rl_start, rl_end    
    INTEGER :: ilat_lo, ilon_lo, ilon_hi, ilat_hi

    INTEGER :: seed_rn, kconseed ! seed part
    REAL(wp):: rn_2d_coarse(sppt_config%coarse_nlon, sppt_config%coarse_nlat)
    REAL(wp):: rn_1d_coarse(sppt_config%coarse_nlon * sppt_config%coarse_nlat)

    LOGICAL :: lzacc ! non-optional version of lacc

    CALL set_acc_host_or_device(lzacc, lacc)

    ! Initializations
    coarse_nlon = sppt_config%coarse_nlon
    coarse_nlat = sppt_config%coarse_nlat

    ! Calculate seed for random number generation
    kconseed = gribout_config(1)%perturbationNumber
    CALL set_seed_rand_numb(mtime_current, kconseed, seed_rn)

    ! Generate random numbers on coarse grid
#ifdef HAVE_ACM_LICENSE
    CALL random_normal_values(seed_rn, sppt_config%range_rn, rn_1d_coarse)
#else
    CALL finish('mo_sppt_util.f90:', "SPPT requires the ACM Software License Agreement when configuring ICON")
#endif
    rn_2d_coarse(:,:) = reshape(rn_1d_coarse, (/SIZE(rn_2d_coarse,1), SIZE(rn_2d_coarse,2)/))

   ! Interpolate form coarse to fine

    rl_start = 1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    !$ACC DATA COPYIN(rn_2d_coarse)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,lon_cell,lat_cell,ilon_lo,ilat_lo,ilon_hi,ilat_hi, &
!$OMP            lon_lo,lat_lo,lon_hi,lat_hi,fq11,fq21,fq12,fq22) ICON_OMP_DEFAULT_SCHEDULE
    ! Interpolate form coarse to fine
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR PRIVATE(lon_cell, lat_cell, ilon_lo, ilat_lo, ilon_hi, ilat_hi, lon_lo, lat_lo, lon_hi, lat_hi, fq11, fq21, fq12, fq22)
      DO jc = i_startidx,i_endidx

        lon_cell = p_patch%cells%center(jc,jb)%lon
        lat_cell = p_patch%cells%center(jc,jb)%lat

        ! indices of coarse cell        
        !   - lon_cell - bbmin_lon is mathematically guarantueed to be >= 0, however, 
        !     in floating point math this might be < 0. hence the call to max()
        !   - dividing by [dlon|dlat]_rn and flooring gives the index of the lower vertex
        !     on the coarse grid. add 1 because of 1 based indexing in fortran
        ilon_lo = floor(max(lon_cell - sppt_config%bbmin_lon, 0._wp)/sppt_config%dlon_rn)+1
        ilat_lo = floor(max(lat_cell - sppt_config%bbmin_lat, 0._wp)/sppt_config%dlat_rn)+1

        ! this is save, the above gives results in the range [1,rn_2d_coarse-1]
        ilon_hi = ilon_lo+1
        ilat_hi = ilat_lo+1

        ! spatial coordinates of coarse cell
        lon_lo = sppt_config%bbmin_lon + (ilon_lo-1)*sppt_config%dlon_rn
        lat_lo = sppt_config%bbmin_lat + (ilat_lo-1)*sppt_config%dlat_rn

        lon_hi = sppt_config%bbmin_lon + (ilon_hi-1)*sppt_config%dlon_rn
        lat_hi = sppt_config%bbmin_lat + (ilat_hi-1)*sppt_config%dlat_rn

        ! function value / random number at coarse cell
        fq11 = rn_2d_coarse(ilon_lo,ilat_lo)
        fq21 = rn_2d_coarse(ilon_hi,ilat_lo)
        fq12 = rn_2d_coarse(ilon_lo,ilat_hi)
        fq22 = rn_2d_coarse(ilon_hi,ilat_hi)        

        rn_2d(jc,jb) = REAL((-1)**kconseed, wp)* &
                       lerp2(x=lon_cell, y=lat_cell, x1=lon_lo, x2=lon_hi, y1=lat_lo, y2=lat_hi, &
                             fq11=fq11, fq21=fq21, fq12=fq12, fq22=fq22)

      END DO
      !$ACC END PARALLEL

    END DO
    !$ACC WAIT(1)

!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC END DATA
    
  END SUBROUTINE generate_rn

  !---------------------------------------------------------------------
  !     
  ! Sets the seed for a random number stream and member number
  !
  ! This subroutine is inherited from the COSMO code. A seed 'seed_rn' is set as a function of 'ydate' 
  !   (which can be e.g. the initial time of the model run), of integer 'kconseed' (which can depend
  !   e.g. on the ensemble member number), and of the input value of 'iseed' itself (which may depend
  !   e.g. on a namelist value). For dates chosen at random, the seeds are approximately uniformly
  !   distributed between 1 and HUGE(0). A highly nonlinear function is used to reduce the possibility 
  !   of correlations between random sequences generated for different initial dates.    
  !
  !---------------------------------------------------------------------

  SUBROUTINE set_seed_rand_numb(mtime_current, kconseed, seed_rn)

    TYPE(datetime),  POINTER, INTENT(IN)        :: mtime_current    !< current_datetime
    INTEGER                 , INTENT(IN)        :: kconseed         !< ID of the member in the ensemble
    INTEGER                 , INTENT(OUT)       :: seed_rn          !< output seed

    !< Local arrays
    !--------------   
    INTEGER, PARAMETER     ::  &
      ibtshf = 11            ! number of bits to shift 'kconseed'  
    INTEGER                ::  &
      iscale,              & ! odd integer between 0 and HUGE-100000000 used for seed scaling
      is,                  & ! utility, integer part of zs
      jdigit, idigits,     & ! loop index and number of digits of ztim
      iradix,              & ! RADIX of ztim
      ndd, nmm, nccaa,     & ! Intermediate integers for the
      nzzaa, nzzmm,        & !    function of the date
      keseed,              & ! function of kconseed (ID of member in the ensemble)
      nindat, nsssss         ! input date and input time

    REAL(wp)               ::  &
      zs, zt, ztim,        & ! utility reals for seed 
      zirr1,               & ! irrational number to randomize bits
      rjudat                 ! Intermediate real for the function of the date

    !   initialize output
    seed_rn = 0
    !   produce a number from'kconseed'
    keseed = (kconseed+1)/2 ! mapping of even/odd enemble ids to the same seed
    ! doing an xorshift32 see https://en.wikipedia.org/wiki/Xorshift
    keseed = IEOR(keseed*2**13, keseed)
    keseed = IEOR(keseed/2**17, keseed)
    keseed = IEOR(keseed*2** 5, keseed)

    !--- generate a unique number from the date and the input keseed
    ! The following transformations are done:
    ! mtime_current -> (nindat, nsssss) -> ... -> rjudat -> ztim -> zs -> seed_rn (output)
    ! the '...' represent intermediate variables ndd, nmm etc.

    ! Part 1: mtime_current -> (nindat, nsssss)
    ! nindat = year * 10 000 + month * 100 + days
    ! nsssss = 3600 * hours + 60 * minutes
    nindat = mtime_current%date%year * 10000 &
           + mtime_current%date%month * 100 &
           + mtime_current%date%day
    nsssss = mtime_current%time%hour * 3600 &
           + mtime_current%time%minute * 60
    
    ! Part 2: (nindat, nsssss) -> rjudat
    nccaa = nindat/10000
    ndd = MOD(nindat,100)
    nmm = MOD((nindat-ndd)/100,100)
    nzzaa = nccaa-( (1-SIGN(1,nmm-3))/2 )
    nzzmm = nmm+6*(1-SIGN(1,nmm-3))
    rjudat = 1720994.5_wp + REAL(2-nzzaa/100  &
                        + (nzzaa/100)/4 &
                        + INT(365.25_wp * REAL(nzzaa,wp))&
                        + INT(30.601_wp * REAL(nzzmm+1,wp))    &
                        + ndd,wp)

    ! Part 3:  (rjudat, nsssss, keseed) -> ztim                
    ztim =   rjudat     &
         - 1720994.5_wp +  REAL(nsssss,wp)/86400.0_wp             &
         - 2581470.3_wp * keseed  

    !--- multiply by an irrational number to randomize the bits and scale
    !--- to between 0 and 1.
    zirr1 = 0.5_wp*(SQRT(5.0_wp)-1.0_wp)
    ztim = FRACTION(zirr1*ABS(ztim))

    ! Part 4:  ztim  -> zs
    !--- reverse the bits
    zs = 0.0_wp
    zt = ztim
    idigits = DIGITS(ztim)
    iradix  = RADIX(ztim)
    DO jdigit = 1, idigits
      zt = zt*iradix
      is = int(zt)
      zt = zt-is
      zs = (zs+is)/iradix
    ENDDO

    ! Part 5:  zs  -> seed_rn
    !--- Scale to an odd number between 0 and HUGE-100000000
    !--- (Allow some headroom in order to use set_seed_rand_numb to set an initial 
    !---  seed and then generate new seeds by incrementing.)
    iscale = (HUGE(seed_rn)-100000000)/2
    seed_rn = 1 + 2*INT( iscale*zs )

  END SUBROUTINE set_seed_rand_numb 

END MODULE mo_sppt_util



