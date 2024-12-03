! # 1 "ifsaux/yomcst.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsaux/yomcst.f90"
! (c) copyright 2014- ecmwf.
!
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.
!
! in applying this licence, ecmwf does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

module yomcst

use parkind1  ,only : jprb

implicit none

public

save

! * rpi          : number pi
real(kind=jprb), parameter :: rpi = 3.14159265358979323846_jprb
! * rclum        : light velocity
real(kind=jprb), parameter :: rclum = 299792458._jprb
! * rhpla        : planck constant
real(kind=jprb), parameter :: rhpla = 6.6260755e-34_jprb
! * rkbol        : bolzmann constant
real(kind=jprb), parameter :: rkbol = 1.380658e-23_jprb
! * rnavo        : avogadro number
real(kind=jprb), parameter :: rnavo = 6.0221367e+23_jprb
! * rsigma       : stefan-bolzman constant
real(kind=jprb), parameter :: rsigma = 5.67037321e-8_jprb ! w m-2 k-4
! * rg           : gravity constant
real(kind=jprb), parameter :: rg = 9.80665_jprb ! m s-2
! * rmd          : dry air molar mass
real(kind=jprb), parameter :: rmd = 28.9644_jprb
! * rmv          : vapour water molar mass
real(kind=jprb), parameter :: rmv = 18.0153_jprb
! * r            : perfect gas constant
real(kind=jprb), parameter :: r = rnavo*rkbol
! * rd           : r_dry (dry air constant)
real(kind=jprb), parameter :: rd = 287.058_jprb! j kg-1 k-1
! * rv           : r_vap (vapour water constant)
real(kind=jprb), parameter :: rv = 1000._jprb*r/rmv
! * rmo3         : ozone molar mass
real(kind=jprb), parameter :: rmo3 = 47.9942_jprb
! * rtt          : tt = temperature of water fusion at "pre_n" 
real(kind=jprb), parameter :: rtt = 273.16_jprb
! * rlvtt        : rlvtt = vaporisation latent heat at t=tt
real(kind=jprb), parameter :: rlvtt = 2.5008e+6_jprb
! * rlstt        : rlstt = sublimation latent heat at t=tt
real(kind=jprb), parameter :: rlstt = 2.8345e+6_jprb
! * ri0          : solar constant
real(kind=jprb), parameter :: ri0 = 1366.0_jprb
! * retv         : r_vap/r_dry - 1
real(kind=jprb), parameter :: retv = rv/rd-1.0_jprb
! * rmco2        : co2 (carbon dioxide) molar mass
real(kind=jprb), parameter :: rmco2 = 44.0095_jprb
! * rmch4        : ch4 (methane) molar mass
real(kind=jprb), parameter :: rmch4 = 16.04_jprb
! * rmn2o        : n2o molar mass
real(kind=jprb), parameter :: rmn2o = 44.013_jprb
! * rmno2        : no2 (nitrogen dioxide) molar mass
real(kind=jprb), parameter :: rmno2 = 46.01_jprb
! * rmcfc11      : cfc11 molar mass
real(kind=jprb), parameter :: rmcfc11 = 137.3686_jprb
! * rmcfc12      : cfc12 molar mass
real(kind=jprb), parameter :: rmcfc12 = 120.914_jprb
! * rmhcfc12     : hcfc22 molar mass
real(kind=jprb), parameter :: rmhcfc22 = 86.469_jprb
! * rmccl4       : ccl4 molar mass
real(kind=jprb), parameter :: rmccl4 = 153.823_jprb

real(kind=jprb), parameter :: rcpd  = 3.5_jprb*rd
real(kind=jprb), parameter :: rlmlt = rlstt-rlvtt

end module yomcst
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

