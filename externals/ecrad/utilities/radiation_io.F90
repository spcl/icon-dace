! # 1 "utilities/radiation_io.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "utilities/radiation_io.f90"
! radiation_io.f90 - provides logging and abort functionality
!
! (c) copyright 2015- ecmwf.
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
!  this file provides an interface to the provision of file units used
!  for logging (nulout and nulerr) and for reading data files
!  (nulrad), as well as an abort routine that should do clean-up
!  appropriate for the environment in which the radiation scheme is
!  embedded.
!
!  rewrite this file as appropriate if the radiation scheme is to be
!  embedded into a model other than the ecmwf integrated forecasting
!  system.

module radiation_io

  ! in the ifs, nulout is equivalent to standard output but only
  ! output from the primary node will be logged, while nulerr is
  ! equivalent to standard error, and text sent to this unit from any
  ! node will be logged. normally, nulerr should only be used before
  ! calling radiation_abort.
  use yomlun_ifsaux, only : nulout, nulerr

  implicit none
  public

  ! this unit may be used for reading radiation configuration files,
  ! but should be closed as soon as the file is read
  integer :: nulrad = 25

contains

  ! abort the program with optional error message. normally you would
  ! log details of the error to nulerr before calling this subroutine.
  subroutine radiation_abort(text)
    character(len=*), intent(in), optional :: text
    if (present(text)) then
      write(nulerr,'(a)') text



      error stop 1

    else



      error stop 'error in radiation scheme'

    end if
  end subroutine radiation_abort

end module radiation_io
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

