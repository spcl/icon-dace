! # 1 "ifsrrtm/srtm_init.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/srtm_init.f90"
subroutine srtm_init(cdirectory, nwvcontinuum)

!-- read in the basic coefficients to configure rrtm_sw
!- creates module yoesrtwn with bg, nspa, nspb, wavenum1, wavenum2,
!  delwave, pref, preflog, tref

use parkind1  ,only : jpim , jprb
use ecradhook   ,only : lhook, dr_hook, jphook

use parsrtm  , only : jpg, jpsw
use yoesrtm  , only : ngn
use yoesrtwn , only : ng, ngm, wt, ngc, rwgt, wtsm
!use yoesrtwn , only : ng, ngm, wt, ngc, ngn, rwgt, wtsm
!use yomlun   , only : nulout

implicit none

character(len=*), intent(in) :: cdirectory

! water vapour continuum model (0=default mt_ckd2.5, 1=caviar)
integer(kind=jpim), intent(in), optional :: nwvcontinuum

! local variables
integer(kind=jpim) :: igc, igcsm, ibnd, ig, ind, ipr, iprsm
real(kind=jprb)    :: zwtsum

real(kind=jphook) :: zhook_handle

!#include "susrtmcf.intfb.h"

! # 1 "./include/susrtm.intfb.h" 1
interface
subroutine susrtm
end subroutine susrtm
end interface
! # 31 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_kgb16.intfb.h" 1
interface
subroutine srtm_kgb16(directory)
character(len=*), intent(in) :: directory
end subroutine srtm_kgb16
end interface
! # 32 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_kgb17.intfb.h" 1
interface
subroutine srtm_kgb17
end subroutine srtm_kgb17
end interface
! # 33 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_kgb18.intfb.h" 1
interface
subroutine srtm_kgb18
end subroutine srtm_kgb18
end interface
! # 34 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_kgb19.intfb.h" 1
interface
subroutine srtm_kgb19
end subroutine srtm_kgb19
end interface
! # 35 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_kgb20.intfb.h" 1
interface
subroutine srtm_kgb20
end subroutine srtm_kgb20
end interface
! # 36 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_kgb21.intfb.h" 1
interface
subroutine srtm_kgb21
end subroutine srtm_kgb21
end interface
! # 37 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_kgb22.intfb.h" 1
interface
subroutine srtm_kgb22
end subroutine srtm_kgb22
end interface
! # 38 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_kgb23.intfb.h" 1
interface
subroutine srtm_kgb23
end subroutine srtm_kgb23
end interface
! # 39 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_kgb24.intfb.h" 1
interface
subroutine srtm_kgb24
end subroutine srtm_kgb24
end interface
! # 40 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_kgb25.intfb.h" 1
interface
subroutine srtm_kgb25
end subroutine srtm_kgb25
end interface
! # 41 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_kgb26.intfb.h" 1
interface
subroutine srtm_kgb26
end subroutine srtm_kgb26
end interface
! # 42 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_kgb27.intfb.h" 1
interface
subroutine srtm_kgb27
end subroutine srtm_kgb27
end interface
! # 43 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_kgb28.intfb.h" 1
interface
subroutine srtm_kgb28
end subroutine srtm_kgb28
end interface
! # 44 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_kgb29.intfb.h" 1
interface
subroutine srtm_kgb29
end subroutine srtm_kgb29
end interface
! # 45 "ifsrrtm/srtm_init.f90" 2
!#include "susrtop.intfb.h"


! # 1 "./include/modify_wv_continuum.intfb.h" 1
interface
subroutine modify_wv_continuum(nwvcontinuum)
use parkind1 , only:&
 & jpim
integer(kind=jpim), intent(in) :: nwvcontinuum
end subroutine modify_wv_continuum
end interface
! # 48 "ifsrrtm/srtm_init.f90" 2


! # 1 "./include/srtm_cmbgb16.intfb.h" 1
interface
subroutine srtm_cmbgb16
end subroutine srtm_cmbgb16
end interface
! # 50 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_cmbgb17.intfb.h" 1
interface
subroutine srtm_cmbgb17
end subroutine srtm_cmbgb17
end interface
! # 51 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_cmbgb18.intfb.h" 1
interface
subroutine srtm_cmbgb18
end subroutine srtm_cmbgb18
end interface
! # 52 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_cmbgb19.intfb.h" 1
interface
subroutine srtm_cmbgb19
end subroutine srtm_cmbgb19
end interface
! # 53 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_cmbgb20.intfb.h" 1
interface
subroutine srtm_cmbgb20
end subroutine srtm_cmbgb20
end interface
! # 54 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_cmbgb21.intfb.h" 1
interface
subroutine srtm_cmbgb21
end subroutine srtm_cmbgb21
end interface
! # 55 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_cmbgb22.intfb.h" 1
interface
subroutine srtm_cmbgb22
end subroutine srtm_cmbgb22
end interface
! # 56 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_cmbgb23.intfb.h" 1
interface
subroutine srtm_cmbgb23
end subroutine srtm_cmbgb23
end interface
! # 57 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_cmbgb24.intfb.h" 1
interface
subroutine srtm_cmbgb24
end subroutine srtm_cmbgb24
end interface
! # 58 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_cmbgb25.intfb.h" 1
interface
subroutine srtm_cmbgb25
end subroutine srtm_cmbgb25
end interface
! # 59 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_cmbgb26.intfb.h" 1
interface
subroutine srtm_cmbgb26
end subroutine srtm_cmbgb26
end interface
! # 60 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_cmbgb27.intfb.h" 1
interface
subroutine srtm_cmbgb27
end subroutine srtm_cmbgb27
end interface
! # 61 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_cmbgb28.intfb.h" 1
interface
subroutine srtm_cmbgb28
end subroutine srtm_cmbgb28
end interface
! # 62 "ifsrrtm/srtm_init.f90" 2

! # 1 "./include/srtm_cmbgb29.intfb.h" 1
interface
subroutine srtm_cmbgb29
end subroutine srtm_cmbgb29
end interface
! # 63 "ifsrrtm/srtm_init.f90" 2

if (lhook) call dr_hook('srtm_init',0,zhook_handle)

!call susrtmcf
call susrtm

!-- read in the molecular absorption coefficients

call srtm_kgb16(cdirectory)
call srtm_kgb17
call srtm_kgb18
call srtm_kgb19
call srtm_kgb20
call srtm_kgb21
call srtm_kgb22
call srtm_kgb23
call srtm_kgb24
call srtm_kgb25
call srtm_kgb26
call srtm_kgb27
call srtm_kgb28
call srtm_kgb29

if (present(nwvcontinuum)) then
  ! modify the shortwave water vapour continuum, if requested
  call modify_wv_continuum(nwvcontinuum)
end if

!-- read in the cloud optical properties
!- creates module yoesrtop with extliq1, ssaliq1, asyliq1, 
!  extice3, ssaice3, asyice3, fdlice3  

!-- rrtm_sw cloud optical properties are not used
!   srtm_cldprop is not called
!   no need to call susrtop

!call susrtop ( -1 )


!mike iacono 20050804
!-- perform g-point reduction from 16 per band (224 total points) to
!-- a band dependent number (112 total points) for all absorption
!-- coefficient input data and planck fraction input data.
!-- compute relative weighting for new g-point combinations.

igcsm = 0
do ibnd = 1,jpsw
  iprsm = 0
  if (ngc(ibnd) < jpg) then
    do igc = 1,ngc(ibnd)
      igcsm = igcsm + 1
      zwtsum = 0.
      do ipr = 1, ngn(igcsm)
        iprsm = iprsm + 1
        zwtsum = zwtsum + wt(iprsm)
      enddo
      wtsm(igc) = zwtsum
    enddo

    do ig = 1,ng(ibnd+15)
      ind = (ibnd-1)*jpg + ig
      rwgt(ind) = wt(ig)/wtsm(ngm(ind))
    enddo
  else
    do ig = 1,ng(ibnd+15)
      igcsm = igcsm + 1
      ind = (ibnd-1)*jpg + ig
      rwgt(ind) = 1.0
    enddo
  endif
enddo

call srtm_cmbgb16
call srtm_cmbgb17
call srtm_cmbgb18
call srtm_cmbgb19
call srtm_cmbgb20
call srtm_cmbgb21
call srtm_cmbgb22
call srtm_cmbgb23
call srtm_cmbgb24
call srtm_cmbgb25
call srtm_cmbgb26
call srtm_cmbgb27
call srtm_cmbgb28
call srtm_cmbgb29

!-----------------------------------------------------------------------
if (lhook) call dr_hook('srtm_init',1,zhook_handle)
end subroutine srtm_init

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

