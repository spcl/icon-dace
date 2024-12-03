! # 1 "ifsrrtm/srtm_kgb26.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/srtm_kgb26.f90"
subroutine srtm_kgb26

!     originally by j.delamere, atmospheric & environmental research.
!     revision: 2.4
!     band 26:  22650-29000 cm-1 (low - nothing; high - nothing)
!     reformatted for f90 by jjmorcrette, ecmwf
!     g.mozdzynski march 2011 read constants from files

!     ------------------------------------------------------------------

use parkind1  , only : jprb
use ecradhook   , only : lhook, dr_hook, jphook
use yoesrta26 , only : sfluxref, rayl 

!     ------------------------------------------------------------------

implicit none

! kurucz
real(kind=jphook) :: zhook_handle
if (lhook) call dr_hook('srtm_kgb26',0,zhook_handle)

sfluxref = (/ &
 !  &     129.462_jprb, 15*_zero_ /)
 & 29.0079_jprb,  28.4088_jprb,     20.3099_jprb,  13.0283_jprb &
 & ,  11.8619_jprb,  9.95840_jprb,     6.68696_jprb,  5.38987_jprb &
 & ,  3.49829_jprb, 0.407693_jprb,    0.299027_jprb, 0.236827_jprb &
 & , 0.188502_jprb, 0.163489_jprb, 4.64335e-02_jprb, 2.72662e-03_jprb /)  

!     rayleigh extinction coefficient at all v 
rayl = (/ &
 & 1.21263e-06_jprb,1.43428e-06_jprb,1.67677e-06_jprb,1.93255e-06_jprb &
 & , 2.19177e-06_jprb,2.44195e-06_jprb,2.66926e-06_jprb,2.85990e-06_jprb &
 & , 3.00380e-06_jprb,3.06996e-06_jprb,3.08184e-06_jprb,3.09172e-06_jprb &
 & , 3.09938e-06_jprb,3.10456e-06_jprb,3.10727e-06_jprb,3.10818e-06_jprb /)  

!     ------------------------------------------------------------------
if (lhook) call dr_hook('srtm_kgb26',1,zhook_handle)
end subroutine srtm_kgb26
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

