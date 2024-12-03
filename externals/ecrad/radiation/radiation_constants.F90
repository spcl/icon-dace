! # 1 "radiation/radiation_constants.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_constants.f90"
! radiation_constants.f90 - constants used in radiation calculations
!
! (c) copyright 2014- ecmwf.
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

module radiation_constants

  use parkind1, only : jprb
!  use yomcst,   only : rpi, rsigma, rg, rd, rmv, rmo3, rhpla

  implicit none
  public

  ! rename some constants from their cryptic ifs names
  real(jprb), parameter :: pi                 = 3.14159265358979323846_jprb
  real(jprb), parameter :: accelduetogravity  = 9.80665_jprb! m s-2
  real(jprb), parameter :: stefanboltzmann    = 5.67037321e-8_jprb ! w m-2 k-4
  real(jprb), parameter :: densityliquidwater = 1000.0_jprb ! kg m-3
  real(jprb), parameter :: densitysolidice    = 916.7_jprb  ! kg m-3
  real(jprb), parameter :: gasconstantdryair  = 287.058_jprb! j kg-1 k-1
  real(jprb), parameter :: planckconstant     = 6.6260695729e-34_jprb ! j s
  real(jprb), parameter :: boltzmannconstant  = 1.380648813e-23_jprb ! j k-1
  real(jprb), parameter :: speedoflight       = 299792458.0_jprb ! m s-1

end module radiation_constants
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

