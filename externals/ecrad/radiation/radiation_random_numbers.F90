! # 1 "radiation/radiation_random_numbers.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_random_numbers.f90"
! this file has been modified for the use in icon

! radiation_random_numbers.f90 - generate random numbers for mcica solver
!
! (c) copyright 2020- ecmwf.
!
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.
!
! in applying this licence, ecmwf does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! author:  robin hogan
! email:   r.j.hogan@ecmwf.int
!
! the derived type "rng_type" is a random number generator that uses
! either (1) fortran's built-in random_number function, or (2) a
! vectorized version of the minstd linear congruential generator.  in
! the case of (2), an rng_type object is initialized with a seed that
! is used to fill up a state of "nmaxstreams" elements using the c++
! minstd_rand0 version of the minstd linear congruential generator
! (lng), which has the form istate[i+1] = mod(istate[i]*a0, m) from
! i=1 to i=nmaxstreams. subsequent requests for blocks of nmaxstreams
! of random numbers use the c++ minstd_ran algorithm in a vectorizable
! form, which modifies the state elements via istate[i] <-
! mod(istate[i]*a, m). uniform deviates are returned that normalize
! the state elements to the range 0-1.
!
! the minstd generator was coded because the random_numbers_mix
! generator in the ifs was found not to vectorize well on some
! hardware.  i am no expert on random number generators, so my
! implementation should really be looked at and improved by someone
! who knows what they are doing.
!
! reference for minstd: park, stephen k.; miller, keith
! w. (1988). "random number generators: good ones are hard to find"
! (pdf). communications of the acm. 31 (10):
! 1192-1201. doi:10.1145/63039.63042
!
! modifications
!   2022-12-01  r. hogan  fixed zeroed state in single precision

module radiation_random_numbers

  use parkind1, only : jprb, jprd, jpim, jpib

  implicit none

  public :: rng_type, irngminstdvector, irngnative, initialize_acc, &
    &  uniform_distribution_acc, iminstda0, iminstda, iminstdm

  enum, bind(c) 
    enumerator irngnative, &    ! built-in fortran-90 rng
         &     irngminstdvector ! vector minstd algorithm
  end enum
  
  ! maximum number of random numbers that can be computed in one call
  ! - this can be increased
  integer(kind=jpim), parameter :: nmaxstreams = 512
  
  ! a requirement of the generator is that the operation mod(a*x,m) is
  ! performed with no loss of precision, so type used for a and x must
  ! be able to hold the largest possible value of a*x without
  ! overflowing, going negative or losing precision. the largest
  ! possible value is 48271*2147483647 = 103661183124337. this number
  ! can be held in either a double-precision real number, or an 8-byte
  ! integer. either may be used, but on some hardwares it has been
  ! found that operations on double-precision reals are faster. select
  ! which you prefer by defining use_real_rng_state for double
  ! precision, or undefining it for an 8-byte integer.


  ! define rng_state_type based on 1, where jprd
  ! refers to a double-precision number regardless of the working
  ! precision described by jprb, while jpib describes an 8-byte
  ! integer






  ! the constants used in the main random number generator
  real(kind=jprd) , parameter :: iminstda  = 48271
  real(kind=jprd) , parameter :: iminstdm  = 2147483647

  ! an alternative value of a that can be used to initialize the
  ! members of the state from a single seed
  real(kind=jprd) , parameter :: iminstda0 = 16807
  
  ! scaling to convert the state to a uniform deviate in the range 0
  ! to 1 in working precision
  real(kind=jprb), parameter :: iminstdscale = 1.0_jprb / real(iminstdm,jprb)

  !---------------------------------------------------------------------
  ! a random number generator type: after being initialized with a
  ! seed, type and optionally a number of vector streams, subsequent
  ! calls to "uniform_distribution" are used to fill 1d or 2d arrays
  ! with random numbers in a way that ought to be fast.
  type rng_type

    integer(kind=jpim) :: itype = irngnative
    real(kind=jprd)     :: istate(nmaxstreams)
    integer(kind=jpim) :: nmaxstreams = nmaxstreams
    integer(kind=jpim) :: iseed

  contains
    procedure :: initialize
    procedure :: uniform_distribution_1d
    procedure :: uniform_distribution_2d
    procedure :: uniform_distribution_2d_masked
    generic   :: uniform_distribution => uniform_distribution_1d, &
         &                               uniform_distribution_2d, &
         &                               uniform_distribution_2d_masked

  end type rng_type

contains

  !---------------------------------------------------------------------
  ! initialize a random number generator, where "itype" may be either
  ! irngnative, indicating to use fortran's native random_number
  ! subroutine, or irngminstdvector, indicating to use the minstd
  ! linear congruential generator (lcg).  in the latter case
  ! "nmaxstreams" should be provided indicating that random numbers
  ! will be requested in blocks of this length. the generator is
  ! seeded with "iseed".
  subroutine initialize(this, itype, iseed, nmaxstreams)

    class(rng_type), intent(inout) :: this
    integer(kind=jpim), intent(in), optional :: itype
    integer(kind=jpim), intent(in), optional :: iseed
    integer(kind=jpim), intent(in), optional :: nmaxstreams

    integer, allocatable :: iseednative(:)
    integer :: nseed, jseed, jstr
    real(jprd) :: rseed ! note this must be in double precision

    if (present(itype)) then
      this%itype = itype
    else
      this%itype = irngnative
    end if
    
    if (present(iseed)) then
      this%iseed = iseed
    else
      this%iseed = 1
    end if

    if (present(nmaxstreams)) then
      this%nmaxstreams = nmaxstreams
    else
      this%nmaxstreams = nmaxstreams
    end if
    
    if (this%itype == irngminstdvector) then
      ! ! option 1: use the c++ minstd_rand0 algorithm to populate the
      ! ! state: this loop is not vectorizable because the state in
      ! ! one stream depends on the one in the previous stream.
      ! this%istate(1) = this%iseed
      ! do jseed = 2,this%nmaxstreams
      !   this%istate(jseed) = mod(iminstda0 * this%istate(jseed-1), iminstdm)
      ! end do

      ! option 2: use a modified (and vectorized) c++ minstd_rand0 algorithm to
      ! populate the state
      rseed = real(abs(this%iseed),jprd)
      do jstr = 1,this%nmaxstreams
        ! note that nint returns an integer of type jpib (8-byte)
        ! which may be converted to double if that is the type of
        ! istate
        this%istate(jstr) = nint(mod(rseed*jstr*(1.0_jprd-0.05_jprd*jstr &
             &      +0.005_jprd*jstr**2)*iminstda0, real(iminstdm,jprd)),kind=jpib)
      end do

      ! one warmup of the c++ minstd_rand algorithm
      do jstr = 1,this%nmaxstreams
        this%istate(jstr) = mod(iminstda * this%istate(jstr), iminstdm)
      end do
      
    else
      ! native generator by default
      call random_seed(size=nseed)
      allocate(iseednative(nseed))
      do jseed = 1,nseed
        iseednative(jseed) = this%iseed + jseed - 1
      end do
      call random_seed(put=iseednative)
      deallocate(iseednative)
    end if

  end subroutine initialize

  !---------------------------------------------------------------------
  ! populate vector "randnum" with pseudo-random numbers; if rannum is
  ! of length greater than nmaxstreams (specified when the generator
  ! was initialized) then only the first nmaxstreams elements will be
  ! assigned.
  subroutine uniform_distribution_1d(this, randnum)

    class(rng_type), intent(inout) :: this
    real(kind=jprb), intent(out)   :: randnum(:)

    integer :: imax, i

    if (this%itype == irngminstdvector) then
      
      imax = min(this%nmaxstreams, size(randnum))

      ! c++ minstd_rand algorithm
      do i = 1, imax
        ! the following calculation is computed entirely with 8-byte
        ! numbers (whether real or integer)
        this%istate(i) = mod(iminstda * this%istate(i), iminstdm)
        ! scale the current state to a number in working precision
        ! (jprb) between 0 and 1
        randnum(i) = iminstdscale * this%istate(i)
      end do

    else

      call random_number(randnum)

    end if

  end subroutine uniform_distribution_1d


  !---------------------------------------------------------------------
  ! populate matrix "randnum" with pseudo-random numbers; if the inner
  ! dimension of rannum is of length greater than nmaxstreams
  ! (specified when the generator was initialized) then only the first
  ! nmaxstreams elements along this dimension will be assigned.
  subroutine uniform_distribution_2d(this, randnum)

    class(rng_type), intent(inout) :: this
    real(kind=jprb), intent(out)   :: randnum(:,:)

    integer :: imax, jblock, i

    if (this%itype == irngminstdvector) then
      
      imax = min(this%nmaxstreams, size(randnum,1))

      ! c++ minstd_ran algorithm
      do jblock = 1,size(randnum,2)
        ! these lines should be vectorizable
        do i = 1, imax
          this%istate(i) = mod(iminstda * this%istate(i), iminstdm)
          randnum(i,jblock) = iminstdscale * this%istate(i)
        end do
      end do

    else

      call random_number(randnum)

    end if

  end subroutine uniform_distribution_2d

  !---------------------------------------------------------------------
  ! populate matrix "randnum" with pseudo-random numbers; if the inner
  ! dimension of rannum is of length greater than nmaxstreams
  ! (specified when the generator was initialized) then only the first
  ! nmaxstreams elements along this dimension will be assigned. this
  ! version only operates on outer dimensions for which "mask" is true.
  subroutine uniform_distribution_2d_masked(this, randnum, mask)

    class(rng_type), intent(inout) :: this
    real(kind=jprb), intent(inout) :: randnum(:,:)
    logical,         intent(in)    :: mask(:)

    integer :: imax, jblock, i

    if (this%itype == irngminstdvector) then
      
      imax = min(this%nmaxstreams, size(randnum,1))

      ! c++ minstd_ran algorithm
      do jblock = 1,size(randnum,2)
        if (mask(jblock)) then
          ! these lines should be vectorizable
          do i = 1, imax
            this%istate(i) = mod(iminstda * this%istate(i), iminstdm)
            randnum(i,jblock) = iminstdscale * this%istate(i)
          end do
        end if
      end do

    else

      do jblock = 1,size(randnum,2)
        if (mask(jblock)) then
          call random_number(randnum(:,jblock))
        end if
      end do

    end if

  end subroutine uniform_distribution_2d_masked

  !---------------------------------------------------------------------
  ! initialize a random number generator, using the minstd
  ! linear congruential generator (lcg). the generator is
  ! seeded with "iseed" and "jseed".
  ! note that this function is not used but manually inlined as the compiler didnot succed.
  ! the seperate function stays in the code, so that hopefully, when the
  ! compiler issue is fixed, it can be used instead of the manual inline.
  pure function initialize_acc(iseed, jseed) result(istate)

    integer(kind=jpim), intent(in)      :: iseed
    integer,            intent(in)      :: jseed

    real(kind=jprb) :: istate

    !$acc routine seq
    
    istate = real(abs(iseed),jprb)
    ! use a modified (and vectorized) c++ minstd_rand0 algorithm to populate the state
    istate = nint(mod( istate*jseed*(1._jprb-0.05_jprb*jseed+0.005_jprb*jseed**2)*iminstda0, iminstdm))
    
    ! one warmup of the c++ minstd_rand algorithm
    istate = mod(iminstda * istate, iminstdm)

  end function initialize_acc

  !---------------------------------------------------------------------
  ! populate vector "randnum" with pseudo-random numbers; if rannum is
  ! of length greater than nmaxstreams (specified when the generator
  ! was initialized) then only the first nmaxstreams elements will be
  ! assigned.
  function uniform_distribution_acc(istate) result(randnum)

    real(kind=jprb), intent(inout) :: istate
    real(kind=jprb)   :: randnum

    !$acc routine seq

    ! c++ minstd_rand algorithm
    istate = mod(iminstda * istate, iminstdm)
    randnum = iminstdscale * istate

  end function uniform_distribution_acc


end module radiation_random_numbers

! #define use_real_rng_state 1
! #define __atomic_acquire 2
! #define __char_bit__ 8
! #define __float_word_order__ __order_little_endian__
! #define __order_little_endian__ 1234
! #define __order_pdp_endian__ 3412
! #define __gfc_real_10__ 1
! #define __finite_math_only__ 0
! #define __gnuc_patchlevel__ 0
! #define __gfc_int_2__ 1
! #define __sizeof_int__ 4
! #define __sizeof_pointer__ 8
! #define __gfortran__ 1
! #define __gfc_real_16__ 1
! #define __stdc_hosted__ 0
! #define __no_math_errno__ 1
! #define __sizeof_float__ 4
! #define __pic__ 2
! #define _language_fortran 1
! #define __sizeof_long__ 8
! #define __gfc_int_8__ 1
! #define __dynamic__ 1
! #define __sizeof_short__ 2
! #define __gnuc__ 13
! #define __sizeof_long_double__ 16
! #define __biggest_alignment__ 16
! #define __atomic_relaxed 0
! #define _lp64 1
! #define __ecrad_little_endian 1
! #define __gfc_int_1__ 1
! #define __order_big_endian__ 4321
! #define __byte_order__ __order_little_endian__
! #define __sizeof_size_t__ 8
! #define __pic__ 2
! #define __sizeof_double__ 8
! #define __atomic_consume 1
! #define __gnuc_minor__ 3
! #define __gfc_int_16__ 1
! #define __lp64__ 1
! #define __atomic_seq_cst 5
! #define __sizeof_long_long__ 8
! #define rng_state_type real(kind=jprd)
! #define __atomic_acq_rel 4
! #define __atomic_release 3
! #define __version__ "13.3.0"

