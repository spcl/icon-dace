! # 1 "ifsrrtm/modify_wv_continuum.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/modify_wv_continuum.f90"
subroutine modify_wv_continuum(nwvcontinuum)

! modify_wv_continuum - adjust the shortwave continuum coefficients
!
! purpose
! -------
!   the default water vapour continuum model in srtm is mt_ckd 2.5,
!   but some measurement programmes, notably from the caviar project
!   (shine et al., j. mol. spectrosc., 2016) suggest a much stronger
!   absorption in the near infrared. this routine provides the option
!   to implement an approximate scaling of the shortwave continuum
!   coefficients to match the caviar continuum. further details on the
!   impact were provided by hogan et al. (2017, ecmwf tech. memo. 816).
!
! interface
! ---------
!   this routine is called from suecrad. if its argument is 0, it does
!   nothing so that the default srtm continuum is used. if its
!   argument is 1 then it implements the caviar continuum by scaling
!   coefficients within the relevant srtm modules.
!
! author
! ------
!   robin hogan, ecmwf
!   original: 2018-02-21
!
! modifications
! -------------
!
! -----------------------------------------------------------------------

use parkind1  , only : jpim, jprb
use ecradhook   , only : lhook, dr_hook, jphook
! load the coefficients for each relevant shortwave band
use yoesrta16, only : selfref16 => selfref, forref16 => forref
use yoesrta17, only : selfref17 => selfref, forref17 => forref
use yoesrta18, only : selfref18 => selfref, forref18 => forref
use yoesrta19, only : selfref19 => selfref, forref19 => forref
use yoesrta20, only : selfref20 => selfref, forref20 => forref
use yoesrta21, only : selfref21 => selfref, forref21 => forref
use yoesrta22, only : selfref22 => selfref, forref22 => forref
use yoesrta23, only : selfref23 => selfref, forref23 => forref
use yoesrta29, only : selfref29 => selfref, forref29 => forref

implicit none

! caviar continuum enhancements
real(kind=jprb), parameter :: self_enh16(16) = (/ 2.42,  2.42,  2.91,  2.91,  2.52, &
     &  2.52,  2.53,  2.53,  2.51,  2.51,  2.51,  2.51,  2.58,  2.58,  2.58,  2.58 /)
real(kind=jprb), parameter :: fore_enh16(16) = (/ 3.38,  3.38,  3.19,  3.19,  1.21,  &
     &  1.21,  1.09,  1.09,  1.07,  1.07,  1.07,  1.07,  1.12,  1.12,  1.12,  1.12 /)
real(kind=jprb), parameter :: self_enh17(16) = (/ 2.18,  1.40,  1.09,  1.19,  1.02,  1.00, &
     &  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00 /)
real(kind=jprb), parameter :: fore_enh17(16) = (/ 3.17,  3.40,  1.66,  1.00,  1.00,  1.00, &
     &  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00 /)
real(kind=jprb), parameter :: self_enh18(16) = (/ 9.67, 12.36,  9.22,  3.71,  1.12,  1.12, &
     &  0.53,  0.53,  0.49,  0.49,  0.49,  0.49,  0.35,  0.35,  0.35,  0.35 /)
real(kind=jprb), parameter :: fore_enh18(16) = (/ 38.90, 15.37, 16.55, 14.81,  4.91,  4.91, &
     &  2.59,  2.59,  2.21,  2.21,  2.21,  2.21,  1.77,  1.77,  1.77,  1.77 /)
real(kind=jprb), parameter :: self_enh19(16) = (/ 28.53, 26.12, 19.14, 10.12,  3.69, &
     &  3.69,  1.63,  1.63,  2.52,  2.52,  2.52,  2.52,  2.40,  2.40,  2.40,  2.40 /)
real(kind=jprb), parameter :: fore_enh19(16) = (/ 11.66,  9.78,  9.57,  9.55,  4.96, &
     &  4.96,  2.68,  2.68,  2.61,  2.61,  2.61,  2.61,  2.37,  2.37,  2.37,  2.37 /)
real(kind=jprb), parameter :: self_enh20(16) = (/ 4.93,  2.76,  1.23,  0.66,  1.41, &
     &  1.11,  1.07,  1.03,  1.03,  1.03,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00 /)
real(kind=jprb), parameter :: fore_enh20(16) = (/ 24.16,  9.04,  2.73,  2.17,  1.05, &
     &  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00 /)
real(kind=jprb), parameter :: self_enh21(16) = (/ 9.70,  4.56,  0.99,  1.21,  1.37, &
     &  1.25,  0.94,  0.99,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00 /)
real(kind=jprb), parameter :: fore_enh21(16) = (/ 50.84, 19.27,  1.49,  1.16,  0.97, &
     &  1.64,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00 /)
real(kind=jprb), parameter :: self_enh22(16) = (/ 3.37,  3.37,  3.37,  3.37,  3.37, &
     &  3.37,  3.37,  3.37,  1.42,  1.42,  1.42,  1.42,  1.42,  1.42,  1.42,  1.42 /)
real(kind=jprb), parameter :: fore_enh22(16) = (/ 12.31, 12.31, 12.31, 12.31, 12.31, &
     &  12.31, 12.31, 12.31,  3.20,  3.20,  3.20,  3.20,  3.20,  3.20,  3.20,  3.20 /)
real(kind=jprb), parameter :: self_enh23(16) = (/ 1.00,  1.00,  1.19,  1.19,  1.65, &
     &  1.46,  1.32,  1.07,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00 /)
real(kind=jprb), parameter :: fore_enh23(16) = (/ 1.04,  1.04,  1.08,  1.08,  1.12, &
     &  1.10,  1.18,  1.06,  1.01,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00 /)
real(kind=jprb), parameter :: self_enh29(16) = (/ 1.70,  1.00,  1.00,  1.03,  1.19,  &
     &  1.19, 1.43,  1.43,  1.30,  1.30,  1.33,  1.33,  1.28,  1.28,  1.08,  1.23 /)
real(kind=jprb), parameter :: fore_enh29(16) = (/ 107.42,  5.87,  3.26,  2.42,  1.39, &
     &  1.39,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00 /)

integer(kind=jpim), intent(in) :: nwvcontinuum

integer(kind=jpim) :: jg

real(kind=jphook) :: zhook_handle

if (lhook) call dr_hook('modify_wv_continuum',0,zhook_handle)
! -----------------------------------------------------------------------

if (nwvcontinuum == 1) then
  ! apply caviar continuum enhancements
  do jg = 1,16
    forref16(:,jg)  = forref16(:,jg)  * fore_enh16(jg)
    selfref16(:,jg) = selfref16(:,jg) * self_enh16(jg)
  enddo
  do jg = 1,16
    forref17(:,jg)  = forref17(:,jg)  * fore_enh17(jg)
    selfref17(:,jg) = selfref17(:,jg) * self_enh17(jg)
  enddo
  do jg = 1,16
    forref18(:,jg)  = forref18(:,jg)  * fore_enh18(jg)
    selfref18(:,jg) = selfref18(:,jg) * self_enh18(jg)
  enddo
  do jg = 1,16
    forref19(:,jg)  = forref19(:,jg)  * fore_enh19(jg)
    selfref19(:,jg) = selfref19(:,jg) * self_enh19(jg)
  enddo
  do jg = 1,16
    forref20(:,jg)  = forref20(:,jg)  * fore_enh20(jg)
    selfref20(:,jg) = selfref20(:,jg) * self_enh20(jg)
  enddo
  do jg = 1,16
    forref21(:,jg)  = forref21(:,jg)  * fore_enh21(jg)
    selfref21(:,jg) = selfref21(:,jg) * self_enh21(jg)
  enddo
  do jg = 1,16
    forref22(:,jg)  = forref22(:,jg)  * fore_enh22(jg)
    selfref22(:,jg) = selfref22(:,jg) * self_enh22(jg)
  enddo
  do jg = 1,16
    forref23(:,jg)  = forref23(:,jg)  * fore_enh23(jg)
    selfref23(:,jg) = selfref23(:,jg) * self_enh23(jg)
  enddo
  do jg = 1,16
    forref29(:,jg)  = forref29(:,jg)  * fore_enh29(jg)
    selfref29(:,jg) = selfref29(:,jg) * self_enh29(jg)
  enddo
endif
  
! -----------------------------------------------------------------------

if (lhook) call dr_hook('modify_wv_continuum',1,zhook_handle)

end subroutine modify_wv_continuum
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

