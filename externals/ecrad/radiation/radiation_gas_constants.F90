! # 1 "radiation/radiation_gas_constants.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_gas_constants.f90"
! this file has been modified for the use in icon

! radiation_gas_constants.f90 - molar mases and id codes of the various gases
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
! license: see the copying file for details
!

module radiation_gas_constants

  use parkind1, only : jprb

  implicit none

  public

  ! gas codes; these indices match those of rrtm-lw up to 7
  integer, parameter :: igasnotpresent = 0
  integer, parameter :: ih2o   = 1
  integer, parameter :: ico2   = 2
  integer, parameter :: io3    = 3
  integer, parameter :: in2o   = 4
  integer, parameter :: ico    = 5
  integer, parameter :: ich4   = 6
  integer, parameter :: io2    = 7
  integer, parameter :: icfc11 = 8
  integer, parameter :: icfc12 = 9
  integer, parameter :: ihcfc22= 10
  integer, parameter :: iccl4  = 11 
  integer, parameter :: ino2   = 12
  integer, parameter :: nmaxgases = 12
  !$acc declare copyin(nmaxgases)

  ! molar masses (g mol-1) of dry air and the various gases above
  real(jprb), parameter :: airmolarmass = 28.970_jprb
  real(jprb), parameter, dimension(0:nmaxgases) :: gasmolarmass = (/ &
       & 0.0_jprb,        & ! gas not present
       & 18.0152833_jprb, & ! h2o
       & 44.011_jprb,     & ! co2
       & 47.9982_jprb,    & ! o3
       & 44.013_jprb,     & ! n2o
       & 28.0101_jprb,    & ! co
       & 16.043_jprb,     & ! ch4
       & 31.9988_jprb,    & ! o2
       & 137.3686_jprb,   & ! cfc11
       & 120.914_jprb,    & ! cfc12
       & 86.469_jprb,     & ! hcfc22
       & 153.823_jprb,    & ! ccl4    
       & 46.0055_jprb /)    ! no2

  ! the corresponding names of the gases in upper and lower case, used
  ! for reading variables from the input file
  character*6, dimension(nmaxgases), parameter :: gasname &
       &  = (/'h2o   ','co2   ','o3    ','n2o   ','co    ','ch4   ', &
       &      'o2    ','cfc11 ','cfc12 ','hcfc22','ccl4  ','no2   '/)
  character*6, dimension(nmaxgases), parameter :: gaslowercasename &
       &  = (/'h2o   ','co2   ','o3    ','n2o   ','co    ','ch4   ', &
       &      'o2    ','cfc11 ','cfc12 ','hcfc22','ccl4  ','no2   '/)

end module radiation_gas_constants
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

