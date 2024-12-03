! # 1 "ifsrrtm/yoesrtwn.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/yoesrtwn.f90"
! this file has been modified for the use in icon

module yoesrtwn

use parkind1  ,only : jpim     ,jprb

implicit none

public

save

!    -------------------------------------------------------------------

integer(kind=jpim) , dimension(16:29) :: ng
integer(kind=jpim) , dimension(16:29) :: nspa
integer(kind=jpim) , dimension(16:29) :: nspb
integer(kind=jpim) , dimension(14)    :: nmpsrtm

real(kind=jprb), dimension(59) :: pref
real(kind=jprb), dimension(59) :: preflog
real(kind=jprb), dimension(59) :: tref

integer(kind=jpim), dimension(224) :: ngm
integer(kind=jpim), dimension(14)  :: ngc, ngs

real(kind=jprb), dimension(16)  :: wt, wtsm
real(kind=jprb), dimension(224) :: rwgt

!$acc declare create(nspa, nspb, preflog, tref, ngc)

! use for 56 g-points
!integer(kind=jpim), dimension(56) :: ngbsw, ngn
! use for 112 g-points
!integer(kind=jpim), dimension(112) :: ngbsw, ngn
! use for 224 g-points
!integer(kind=jpim), dimension(224) :: ngbsw, ngn

!     -----------------------------------------------------------------
!        * e.c.m.w.f. physics package ** rrtm sw radiation **

!     j.-j. morcrette       e.c.m.w.f.      03-03-07
!     m. j. iacono          aer             12/09/03

!  name     type     purpose
!  ----   : ----    : -------
!  ng     : integer : number of k-coefficients in spectral intervals
!  nspa   : integer :
!  nspb   : integer :
! nmpsrtm : integer : mapping indices for 6-spectral int. surface albedo 
! wavenum1: real    : lower wavenumber spectral limit
! wavenum2: real    : higher wavenumber spectral limit
! delwave : real    : spectral interval width
! pref    : real    : reference pressure
! preflog : real    : log reference pressure
! tref    : real    : reference temperature
!  ngc    : integer : the number of new g-points in each band
!  ngs    : integer : the cumulative sum of new g-points for each band
!  ngm    : integer : the index of each new g-point relative to the
!                     original 16 g-points for each band.
!
!  wt     : real    : rrtm weights for 16 g-points.
!  wtsum  : real    : sum of the weights
!  rwgt   : real    :
! 
!  ngn    : integer : the number of original g-points that are combined 
!                     to make each new g-point in each band.
!  ngbsw  : integer : the band index for each new g-point.
!     -----------------------------------------------------------------
end module yoesrtwn

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

