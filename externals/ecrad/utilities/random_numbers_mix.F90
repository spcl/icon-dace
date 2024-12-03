! # 1 "utilities/random_numbers_mix.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "utilities/random_numbers_mix.f90"
!**** *random_numbers_mix*  - portable random number generator

! (c) copyright 2002- ecmwf.
!
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.
!
! in applying this licence, ecmwf does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!     purpose.
!     --------
!           generate machine-independent pseudo-random numbers

!**   interface.
!     ----------
!        call initialize_random_numbers (kseed, yd_stream)
!        call uniform_distribution      (px   , yd_stream)
!        call gaussian_distribution     (px   , yd_stream)
!        call random_number_restartfile (fname, action)
!        call wr_rangen_state           (nunit)

!        explicit arguments :
!        --------------------
!        kseed  (input)    : integer seed in the range [0,huge(kseed)]
!        yd_stream (optional) : the state of the random number generator
!        px     (output)   : array to receive random numbers in the range

!        in the case of uniform_distribution, px has values in the range [0.0,1.0)

!        implicit arguments :
!        --------------------
!        none

!     method.
!     -------
!        based loosly on zufall (petersen, 1994).

!        the main difference between this generator and zufall is that integer arithmetic
!        is used. this ensures portability to vector machines that implement different
!        real arithmetic. in particular, vector machines often implement non-ieee
!        arithmetic for their vector pipes. this routine will give identical results for
!        any integer type with at least 32 bits.

!        the generator is a lagged-fibonacci generator: x(i) = x(i-p) + x(i-q) mod 2**m.
!        lagged-fibonacci generators have very long repeat periods: (2**q -1) * 2**(m-1)
!        (i.e about 2.85e191 for q=607, m=30). they pass most tests for randomness.

!        p and q must be chosen carefully. values from the following table are ok.
!        larger values give better random numbers, but smaller values are more
!        cache-friendly.

!          q         p
!        9689      4187
!        4423      2098
!        2281      1029
!        1279       418
!         607       273
!         521       168
!         250       103
!         127        63
!          97        33
!          55        24
!          43        22
!          31        13
!          24        10

!        the initial q values of x are set using the binary shirt register method of
!        burns and pryor 1999.

!        mascagni et al (1995) show how to choose different sets of initial values that
!        are guaranteed to be drawn from different maximal-length cycles. this requires
!        the initial values of x(1)...x(q) to be in "canonical form". specifically,
!        x(1) must be zero and all but a particular one or two values of x must be
!        even. for q=607 and p=273, only one element (jpq-jps) must be odd.

!     externals.
!     ----------
!        none

!     reference.
!     ----------
!        burns p.j. and pryor d.v. 1999,
!                             surface radiative transport at large scale via monte carlo.
!                             annual review of heat transfer, vol 9.
!
!        petersen w.p., 1994, lagged fibonacci series random number generator
!                             for the nec sx-3. international journal of high speed computing
!                             vol. 6, no. 3, pp387-398.
!
!        mascagni m., cuccaro s.a., pryor d.v., robinson m.l., 1995,
!                             a fast, high quality and reproducible parallel lagged-fibonacci
!                             pseudorandom number generator. journal of computational physics
!                             vol 119. pp211-219.

!     author.
!     -------
!        mike fisher *ecmwf*

!     modifications.
!     --------------
!        original : 2002-09-25
!        made parallel friendly: 2003-08-11 robert pincus
!        m leutbecher: 2004-05-10 restart capability
!        m fisher:     2005-03-30 replaced lcg initialization with shift register
!     ------------------------------------------------------------------




module random_numbers_mix
use ecradhook,  only : lhook, dr_hook, jphook
use parkind1, only : jpim, jprb

implicit none

save

private
public randomnumberstream,wr_rangen_state, &
     & initialize_random_numbers, uniform_distribution, gaussian_distribution ! ,&
!     & random_number_restartfile, 

integer(kind=jpim), parameter      :: jpp=273, jpq=607, jps=105
integer(kind=jpim), parameter      :: jpmm=30
integer(kind=jpim), parameter      :: jpm=2**jpmm
integer(kind=jpim), parameter      :: jpnumsplit=(jpq-2)/(jpp-1)
integer(kind=jpim), parameter      :: jplensplit=(jpq-jpp+jpnumsplit-1)/jpnumsplit
integer(kind=jpim), parameter      :: initvalue = 12345678

type randomnumberstream
  private
  integer(kind=jpim)                 :: iused
  integer(kind=jpim)                 :: inittest ! should initialize to zero, but can't in f90
  integer(kind=jpim), dimension(jpq) :: ix 
  real(kind=jprb)                    :: zrm
end type randomnumberstream

contains
!-------------------------------------------------------------------------------
subroutine initialize_random_numbers (kseed, yd_stream) 
  !-------------------------------------------------------------------------------
  ! initialize fibgen
  !-------------------------------------------------------------------------------
  integer(kind=jpim),                intent(in   ) :: kseed
  type(randomnumberstream), intent(inout) :: yd_stream
  
  integer, parameter :: jpmask=123459876
  integer(kind=jpim), parameter     :: jpwarmup_shft=64, jpwarmup_lfg=999
  integer(kind=jpim)                :: idum,jk,jj,jbit
  real(kind=jprb), dimension(jpwarmup_lfg)   :: zwarmup

  !-------------------------------------------------------------------------------
  ! initialize the buffer using a binary shift register (burns and pryor, 1999).
  ! the galois representation is used for the shift register as it is more
  ! efficient than the fibonacci representation. the magic numbers 31 and 87
  ! define the shift register primitive polynomial=(32,7,5,3,2,1,0).
  !
  ! to ensure that different seeds produce distinct initial buffer states in
  ! canonical form, bits 0...jpmm-2 of the initial seed (after xoring with jpmask
  ! and spinning up using the linear congruential generator) are used to construct
  ! x(2), and the remaining bits are used to construct x(jpq).
  !-------------------------------------------------------------------------------
  
  real(kind=jphook) :: zhook_handle
  if (lhook) call dr_hook('random_numbers_mix:initialize_random_numbers',0,zhook_handle)
  idum = abs(ieor(kseed,jpmask))
  if (idum==0) idum=jpmask

  do jj=1,jpwarmup_shft
    if (btest(idum,31)) then
      idum=ibset(ishft(ieor(idum,87),1),0)
    else
      idum=ibclr(ishft(idum,1),0)
    endif
  enddo

  yd_stream%ix(1:jpq-1)= 0
  yd_stream%ix(2)      = ishft(ibits(idum,0,jpmm-1),1)
  yd_stream%ix(jpq)    = ibits(idum,jpmm-1,bit_size(idum)+1-jpmm)

  do jbit=1,jpmm-1
    do jj=3,jpq-1
      if (btest(idum,31)) then
        idum=ibset(ishft(ieor(idum,87),1),0)
        yd_stream%ix(jj)=ibset(yd_stream%ix(jj),jbit)
      else
        idum=ibclr(ishft(idum,1),0)
      endif
    enddo
  enddo

  yd_stream%ix(jpq-jps) = ibset(yd_stream%ix(jpq-jps),0)
  
  !-------------------------------------------------------------------------------
  ! initialize some constants
  !-------------------------------------------------------------------------------
  
  yd_stream%iused=jpq
  yd_stream%zrm=1.0_jprb/real(jpm,jprb)
  
  !-------------------------------------------------------------------------------
  ! check the calculation of jpnumsplit and jplensplit.
  !-------------------------------------------------------------------------------
  
  if (jpp+jpnumsplit*jplensplit < jpq) then
    call abor1 ('initialize_random_numbers: upper limit of last loop < jpq')
  endif
  
  if (jplensplit >=jpp) then
    call abor1 ('initialize_random_numbers: loop length > jpp')
  endif
  
  if (jpnumsplit>1) then
    if ((jpq-jpp+jpnumsplit-2)/(jpnumsplit-1) < jpp) then
      call abor1 ('initialize_random_numbers: jpnumsplit is bigger than necessary')
    endif
  endif

  !-------------------------------------------------------------------------------
  ! set inittest to show that the stream is initialized.
  !-------------------------------------------------------------------------------

  yd_stream%inittest = initvalue
  
  !-------------------------------------------------------------------------------
  ! warm up the generator.
  !-------------------------------------------------------------------------------

  call uniform_distribution (zwarmup, yd_stream)

if (lhook) call dr_hook('random_numbers_mix:initialize_random_numbers',1,zhook_handle)
end subroutine initialize_random_numbers

!@process hot nostrict
subroutine uniform_distribution (px,yd_stream)
  !--------------------------------------------------------------------------------
  ! generate uniformly distributed random numbers in the range 0.0<= px < 1.0
  !--------------------------------------------------------------------------------
  integer(kind=jpim), parameter :: ivar = int(z"3fffffff",jpim)
  type(randomnumberstream), intent(inout) :: yd_stream
  real(kind=jprb), dimension(:),     intent(  out) :: px

  integer(kind=jpim)                :: jj, jk, in, ifilled
  
  ! this test is a little dirty but fortran 90 doesn't allow for the initialization
  !   of components of derived types. 
  real(kind=jphook) :: zhook_handle
! dr_hook removed to reduce overhead
! if (lhook) call dr_hook('random_numbers_mix:uniform_distribution',0,zhook_handle)
  if(yd_stream%inittest /= initvalue) &
    & call abor1 ('uniform_distribution called before initialize_random_numbers')

  !--------------------------------------------------------------------------------
  ! copy numbers that were generated during the last call, but not used.
  !--------------------------------------------------------------------------------
  
  in=size(px)
  ifilled=0
  
  do jj=yd_stream%iused+1,min(jpq,in+yd_stream%iused)
    px(jj-yd_stream%iused) = yd_stream%ix(jj)*yd_stream%zrm
    ifilled=ifilled+1
  enddo
  
  yd_stream%iused=yd_stream%iused+ifilled
  
  if (ifilled==in)  then 
!   if (lhook) call dr_hook('random_numbers_mix:uniform_distribution',1,zhook_handle)
! dr_hook removed to reduce overhead
    return
  endif
  
  !--------------------------------------------------------------------------------
  ! generate batches of jpq numbers until px has been filled
  !--------------------------------------------------------------------------------
  
  do while (ifilled<in)
  
  !--------------------------------------------------------------------------------
  ! generate jpq numbers in vectorizable loops. the first loop is length jpp. the
  ! remaining jpq-jpp elements are calculated in loops of length shorter than jpp.
  !--------------------------------------------------------------------------------
  
  !ocl novrec
    do jj=1,jpp
!     yd_stream%ix(jj) = yd_stream%ix(jj) + yd_stream%ix(jj-jpp+jpq)
!     if (yd_stream%ix(jj)>=jpm) yd_stream%ix(jj) = yd_stream%ix(jj)-jpm
      yd_stream%ix(jj) = iand(ivar,yd_stream%ix(jj) + yd_stream%ix(jj-jpp+jpq))
    enddo
  
    do jk=1,jpnumsplit
  !ocl novrec
      do jj=1+jpp+(jk-1)*jplensplit,min(jpq,jpp+jk*jplensplit)
!       yd_stream%ix(jj) = yd_stream%ix(jj) + yd_stream%ix(jj-jpp)
!       if (yd_stream%ix(jj)>=jpm) yd_stream%ix(jj) = yd_stream%ix(jj)-jpm
        yd_stream%ix(jj) = iand(ivar,yd_stream%ix(jj) + yd_stream%ix(jj-jpp))
      enddo
    enddo
  
    yd_stream%iused = min(jpq,in-ifilled)
    px(ifilled+1:ifilled+yd_stream%iused) = yd_stream%ix(1:yd_stream%iused)*yd_stream%zrm
    ifilled = ifilled+yd_stream%iused
  enddo
  
!if (lhook) call dr_hook('random_numbers_mix:uniform_distribution',1,zhook_handle)
! dr_hook removed to reduce overhead
end subroutine uniform_distribution
!-------------------------------------------------------------------------------
subroutine gaussian_distribution (px, yd_stream)
  type(randomnumberstream), intent(inout) :: yd_stream
  real(kind=jprb),                   intent(  out) :: px(:)
  !--------------------------------------------------------------------------------
  ! generate normally-distributed random numbers using the box-muller method.
  !
  ! nb: this routine does not use buffering. this means that the following calls:
  !     call gaussian_distribution (zx(1:k))
  !     call gaussian_distribution (zx(k+1:n))
  ! will produce different numbers for elements k+1 onwards than the single call:
  !     call gaussian_distribution (zx(1:n))
  !--------------------------------------------------------------------------------
  
  integer(kind=jpim) :: ilen, j
  real(kind=jprb) :: zfac, ztwopi
  real(kind=jprb) :: zx(size(px)+1)
  
  !--------------------------------------------------------------------------------
  ! generate uniform random points in the range [0,1)
  !--------------------------------------------------------------------------------

    real(kind=jphook) :: zhook_handle
    if (lhook) call dr_hook('random_numbers_mix:gaussian_distribution',0,zhook_handle)
    call uniform_distribution (zx, yd_stream)

  !--------------------------------------------------------------------------------
  ! generate gaussian deviates using box-muller method
  !--------------------------------------------------------------------------------
  
  ztwopi = 8.0_jprb*atan(1.0_jprb)
  ilen=size(px)
  
  do j=1,ilen-1,2
    zfac = sqrt(-2.0_jprb*log(1.0_jprb-zx(j)))
    px(j  ) = zfac*cos(ztwopi*zx(j+1))
    px(j+1) = zfac*sin(ztwopi*zx(j+1))
  enddo
  
  !--------------------------------------------------------------------------------
  ! generate the last point if ilen is odd
  !--------------------------------------------------------------------------------
  
  if (mod(ilen,2) /= 0) then
    zfac = sqrt(-2.0_jprb*log(1.0_jprb-zx(ilen)))
    px(ilen) = zfac*cos(ztwopi*zx(ilen+1))
  endif
  
if (lhook) call dr_hook('random_numbers_mix:gaussian_distribution',1,zhook_handle)
end subroutine gaussian_distribution
!-------------------------------------------------------------------------------
!!$subroutine random_number_restartfile( cdfname, cdaction,yd_stream )
!!$!--------------------------------------------------------------------------------
!!$!
!!$! read (cdaction='r') or write (cdaction='w') restart file
!!$! for random number generator
!!$!
!!$!--------------------------------------------------------------------------------
!!$character (len=*),   intent(in) :: cdfname
!!$character (len=1  ), intent(in) :: cdaction
!!$type(randomnumberstream), intent(inout) :: yd_stream
!!$  
!!$integer(kind=jpim) :: iunit, iret, ibytes_in_jpim
!!$
!!$real(kind=jprb) :: zhook_handle
!!$if (lhook) call dr_hook('random_numbers_mix:random_number_restartfile',0,zhook_handle)
!!$ibytes_in_jpim= ceiling(real(bit_size(yd_stream%iused))/8.0_jprb - tiny(1.0_jprb))
!!$
!!$if (ibytes_in_jpim /= 4) then
!!$  call abor1('random_number_restartfile: number of bytes for jpim is not 4 ')        
!!$endif
!!$
!!$call pbopen(iunit, cdfname, cdaction, iret)
!!$if (iret /= 0) then
!!$  call abor1('random_number_restartfile: pbopen failed opening '//cdfname)    
!!$endif
!!$
!!$
!!$if (cdaction=='r' .or. cdaction=='r') then
!!$  call pbread(iunit, yd_stream%ix,    ibytes_in_jpim*jpq, iret)
!!$  if (iret < 0) then
!!$    call abor1('random_number_restartfile: pbread could not read ix from '//cdfname)    
!!$  endif
!!$  call pbread(iunit, yd_stream%iused, ibytes_in_jpim    , iret)
!!$  if (iret < 0) then
!!$    call abor1('random_number_restartfile: pbread could not read iused from '//cdfname)    
!!$  endif
!!$
!!$!  l_initialized = .true.
!!$  yd_stream%inittest = initvalue
!!$  yd_stream%zrm=1.0_jprb/real(jpm,jprb)
!!$elseif(cdaction=='w' .or. cdaction=='w') then
!!$  call pbwrite(iunit, yd_stream%ix, ibytes_in_jpim*jpq, iret)
!!$  if (iret < 0) then
!!$    call abor1('random_number_restartfile: pbwrite could not write ix on '//cdfname)    
!!$  endif
!!$  call pbwrite(iunit, yd_stream%iused, ibytes_in_jpim , iret)
!!$  if (iret < 0) then
!!$    call abor1('random_number_restartfile: pbwrite could not write iused on '//cdfname)    
!!$  endif
!!$
!!$else
!!$  call abor1 ('random_number_restartfile: cdaction = '//cdaction//' is undefined.')
!!$endif
!!$
!!$call pbclose(iunit, iret)
!!$if (iret /= 0) then
!!$  call abor1('random_number_restartfile: pbclose failed closing '//cdfname)    
!!$endif
!!$
!!$if (lhook) call dr_hook('random_numbers_mix:random_number_restartfile',1,zhook_handle)
!!$end subroutine random_number_restartfile


subroutine wr_rangen_state( kunit, yd_stream )
!--------------------------------------------------------------------------------
! write state of random number generator to unit kunit
!--------------------------------------------------------------------------------
integer(kind=jpim), intent(in) :: kunit
type(randomnumberstream), intent(in) :: yd_stream

real(kind=jphook) :: zhook_handle
if (lhook) call dr_hook('random_numbers_mix:wr_rangen_state',0,zhook_handle)
write( kunit, * ) 'module random_numbers_mix, generator state is'
write( kunit, '(8i10)') yd_stream%ix
write( kunit, '(i10)')  yd_stream%iused

if (lhook) call dr_hook('random_numbers_mix:wr_rangen_state',1,zhook_handle)
end subroutine wr_rangen_state

end module random_numbers_mix
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
! #define __atomic_acq_rel 4
! #define __atomic_release 3
! #define __version__ "13.3.0"

