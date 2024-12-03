! # 1 "radiation/radiation_check.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_check.f90"
! radiation_check.f90 - checking routines
!
! (c) copyright 2020- ecmwf.
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

module radiation_check

  use parkind1, only : jprb

  implicit none
  public

contains

  !---------------------------------------------------------------------
  ! return .true. if 1d allocatable array "var" is out of physical
  ! range specified by boundmin and boundmax, and issue a warning.
  ! "do_fix" determines whether erroneous values are fixed to lie
  ! within the physical range. to check only a subset of the array,
  ! specify i1 and i2 for the range.
  function out_of_bounds_1d(var, var_name, boundmin, boundmax, do_fix, i1, i2) result (is_bad)

    use radiation_io,     only : nulout

    real(jprb), allocatable, intent(inout) :: var(:)
    character(len=*),        intent(in) :: var_name
    real(jprb),              intent(in) :: boundmin, boundmax
    logical,                 intent(in) :: do_fix
    integer,       optional, intent(in) :: i1, i2

    logical                       :: is_bad

    real(jprb) :: varmin, varmax

    is_bad = .false.

    if (allocated(var)) then

      if (present(i1) .and. present(i2)) then
        varmin = minval(var(i1:i2))
        varmax = maxval(var(i1:i2))
      else
        varmin = minval(var)
        varmax = maxval(var)
      end if

      if (varmin < boundmin .or. varmax > boundmax) then
        write(nulout,'(a,a,a,g12.4,a,g12.4,a,g12.4,a,g12.4)',advance='no') &
             &  '*** warning: ', var_name, ' range', varmin, ' to', varmax, &
             &  ' is out of physical range', boundmin, 'to', boundmax
        is_bad = .true.
        if (do_fix) then
          if (present(i1) .and. present(i2)) then
            var(i1:i2) = max(boundmin, min(boundmax, var(i1:i2)))
          else
            var = max(boundmin, min(boundmax, var))
          end if
          write(nulout,'(a)') ': corrected'
        else
          write(nulout,'(1x)')
        end if
      end if

    end if
    
  end function out_of_bounds_1d


  !---------------------------------------------------------------------
  ! return .true. if 2d allocatable array "var" is out of physical
  ! range specified by boundmin and boundmax, and issue a warning.  to
  ! check only a subset of the array, specify i1 and i2 for the range
  ! of the first dimension and j1 and j2 for the range of the second.
  function out_of_bounds_2d(var, var_name, boundmin, boundmax, do_fix, &
       &                    i1, i2, j1, j2) result (is_bad)

    use radiation_io,     only : nulout

    real(jprb), allocatable, intent(inout) :: var(:,:)
    character(len=*),        intent(in) :: var_name
    real(jprb),              intent(in) :: boundmin, boundmax
    logical,                 intent(in) :: do_fix
    integer,       optional, intent(in) :: i1, i2, j1, j2

    ! local copies of indices
    integer :: ii1, ii2, jj1, jj2

    logical                       :: is_bad

    real(jprb) :: varmin, varmax

    is_bad = .false.

    if (allocated(var)) then

      if (present(i1) .and. present(i2)) then
        ii1 = i1
        ii2 = i2
      else
        ii1 = lbound(var,1)
        ii2 = ubound(var,1)
      end if
      if (present(j1) .and. present(j2)) then
        jj1 = j1
        jj2 = j2
      else
        jj1 = lbound(var,2)
        jj2 = ubound(var,2)
      end if
      varmin = minval(var(ii1:ii2,jj1:jj2))
      varmax = maxval(var(ii1:ii2,jj1:jj2))

      if (varmin < boundmin .or. varmax > boundmax) then
        write(nulout,'(a,a,a,g12.4,a,g12.4,a,g12.4,a,g12.4)',advance='no') &
             &  '*** warning: ', var_name, ' range', varmin, ' to', varmax,&
             &  ' is out of physical range', boundmin, 'to', boundmax
        is_bad = .true.
        if (do_fix) then
          var(ii1:ii2,jj1:jj2) = max(boundmin, min(boundmax, var(ii1:ii2,jj1:jj2)))
          write(nulout,'(a)') ': corrected'
        else
          write(nulout,'(1x)')
        end if
      end if

    end if
    
  end function out_of_bounds_2d


  !---------------------------------------------------------------------
  ! return .true. if 3d allocatable array "var" is out of physical
  ! range specified by boundmin and boundmax, and issue a warning.  to
  ! check only a subset of the array, specify i1 and i2 for the range
  ! of the first dimension, j1 and j2 for the second and k1 and k2 for
  ! the third.
  function out_of_bounds_3d(var, var_name, boundmin, boundmax, do_fix, &
       &                    i1, i2, j1, j2, k1, k2) result (is_bad)

    use radiation_io,     only : nulout

    real(jprb), allocatable, intent(inout) :: var(:,:,:)
    character(len=*),        intent(in) :: var_name
    real(jprb),              intent(in) :: boundmin, boundmax
    logical,                 intent(in) :: do_fix
    integer,       optional, intent(in) :: i1, i2, j1, j2, k1, k2

    ! local copies of indices
    integer :: ii1, ii2, jj1, jj2, kk1, kk2

    logical                       :: is_bad

    real(jprb) :: varmin, varmax

    is_bad = .false.

    if (allocated(var)) then

      if (present(i1) .and. present(i2)) then
        ii1 = i1
        ii2 = i2
      else
        ii1 = lbound(var,1)
        ii2 = ubound(var,1)
      end if
      if (present(j1) .and. present(j2)) then
        jj1 = j1
        jj2 = j2
      else
        jj1 = lbound(var,2)
        jj2 = ubound(var,2)
      end if
      if (present(k1) .and. present(k2)) then
        kk1 = k1
        kk2 = k2
      else
        kk1 = lbound(var,3)
        kk2 = ubound(var,3)
      end if
      varmin = minval(var(ii1:ii2,jj1:jj2,kk1:kk2))
      varmax = maxval(var(ii1:ii2,jj1:jj2,kk1:kk2))

      if (varmin < boundmin .or. varmax > boundmax) then
        write(nulout,'(a,a,a,g12.4,a,g12.4,a,g12.4,a,g12.4)',advance='no') &
             &  '*** warning: ', var_name, ' range', varmin, ' to', varmax,&
             &  ' is out of physical range', boundmin, 'to', boundmax
        is_bad = .true.
        if (do_fix) then
          var(ii1:ii2,jj1:jj2,kk1:kk2) = max(boundmin, min(boundmax, &
               &                             var(ii1:ii2,jj1:jj2,kk1:kk2)))
          write(nulout,'(a)') ': corrected'
        else
          write(nulout,'(1x)')
        end if
      end if

    end if
    
  end function out_of_bounds_3d

end module radiation_check
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

