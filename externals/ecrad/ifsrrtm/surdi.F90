! # 1 "ifsrrtm/surdi.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/surdi.f90"
! this file has been modified for the use in icon

subroutine surdi

!**** *surdi*   - initialize common yoerdi controlling radint

!     purpose.
!     --------
!           initialize yoerdi, the common that controls the
!           radiation interface

!**   interface.
!     ----------
!        call *surdi* from *surad*
!              -----        -----

!        explicit arguments :
!        --------------------
!        none

!        implicit arguments :
!        --------------------
!        common yoerdi

!     method.
!     -------
!        see documentation

!     externals.
!     ----------
!        none

!     reference.
!     ----------
!        ecmwf research department documentation of the ifs model

!     author.
!     -------
!      original  jean-jacques morcrette  *ecmwf*
!      original : 88-12-15

!     modifications.
!     --------------
!      m.hamrud      01-oct-2003 cy28 cleaning
!      modified   p. viterbo   24-05-2004  surf library
!      jjmorcrette   2004-10-07 gas concentrations
!     ------------------------------------------------------------------

use parkind1  ,only : jprb
use ecradhook   ,only : lhook, dr_hook, jphook

use yoerdi   , only : rrae     ,&
 & rcardi   ,rch4     ,rn2o     ,rno2     ,ro3      ,&
 & rcfc11   ,rcfc12   ,rcfc22   ,rccl4    ,&
 & repclc   ,reph2o   ,rsundur  ,&
 & rcco2    ,rcch4    ,rcn2o    ,rcno2    ,rccfc11  ,&
 & rccfc12  ,rccfc22  ,rcccl4
use yomdyncore, only : laqua

implicit none

real(kind=jprb) :: zairmwg, zc11mwg, zc12mwg, zch4mwg, zco2mwg,&
 & zn2omwg, zno2mwg, zo3mwg, zc22mwg, zcl4mwg
real(kind=jphook) :: zhook_handle

!      ----------------------------------------------------------------

!*       1.    set default values.
!              -------------------

if (lhook) call dr_hook('surdi',0,zhook_handle)
rrae = 0.1277e-02_jprb

!* threshold for computing sunshine duration (w/m2)
rsundur=120._jprb

!*  for sea ice, monthly values are based on ebert and curry, 1993, table 2.
!   we take dry snow albedo as the representative value for non-summer
!   months, and bare sea-ice as the representative value for summer
!   months. the values for antarctic are shifted six-months.
! all computations brought back to *suswn*

!*  concentration of the various trace gases (ipcc/sacc values for 1990)
!        co2         ch4        n2o        cfc11       cfc12
!      353ppmv     1.72ppmv   310ppbv     280pptv     484pptv

zairmwg = 28.970_jprb
zco2mwg = 44.011_jprb
zch4mwg = 16.043_jprb
zn2omwg = 44.013_jprb
zno2mwg = 46.006_jprb
zo3mwg  = 47.9982_jprb
zc11mwg = 137.3686_jprb
zc12mwg = 120.9140_jprb
zc22mwg =  86.4690_jprb
zcl4mwg = 153.8230_jprb

!rcardi  = 353.e-06_jprb*zco2mwg/zairmwg
!rch4    = 1.72e-06_jprb*zch4mwg/zairmwg
!rn2o    = 310.e-09_jprb*zn2omwg/zairmwg
!rno2    = 500.e-13_jprb*zno2mwg/zairmwg
!ro3     =   1.e-06_jprb*zo3mwg /zairmwg
!rcfc11  = 280.e-12_jprb*zc11mwg/zairmwg
!rcfc12  = 484.e-12_jprb*zc12mwg/zairmwg
!rcfc22  =   1.e-12_jprb*zc22mwg/zairmwg
!rccl4   =   1.e-12_jprb*zcl4mwg/zairmwg

!> remove this block as the variables not used and accesses undefined values
!if( laqua ) then
!  rcardi  = 348.e-06_jprb*zco2mwg/zairmwg
!  rch4    = 1.65e-06_jprb*zch4mwg/zairmwg
!  rn2o    = 306.e-09_jprb*zn2omwg/zairmwg
!else
!  rcardi  = rcco2   * zco2mwg/zairmwg
!  rch4    = rcch4   * zch4mwg/zairmwg
!  rn2o    = rcn2o   * zn2omwg/zairmwg
!endif
!
!rno2    = rcno2   * zno2mwg/zairmwg 
!ro3     = 1.e-06_jprb*zo3mwg /zairmwg
!rcfc11  = rccfc11 * zc11mwg/zairmwg
!rcfc12  = rccfc12 * zc12mwg/zairmwg
!rcfc22  = rccfc22 * zc22mwg/zairmwg
!rccl4   = rcccl4  * zcl4mwg/zairmwg
!<end

repclc=1.e-12_jprb
reph2o=1.e-12_jprb

!     -----------------------------------------------------------------

if (lhook) call dr_hook('surdi',1,zhook_handle)
end subroutine surdi
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

