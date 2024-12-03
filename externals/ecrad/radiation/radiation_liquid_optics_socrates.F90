! # 1 "radiation/radiation_liquid_optics_socrates.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_liquid_optics_socrates.f90"
! this file has been modified for the use in icon

! radiation_liquid_optics_socrates.f90 - socrates method for parameterizing liquid droplet optics
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
! modifications
!   2020-08-10  r. hogan  bounded re to be >=1.2um and <=50um

module radiation_liquid_optics_socrates

  use parkind1, only : jprb

  implicit none
  public

  ! socrates (edwards-slingo) parameterizes info on the dependence of
  ! the scattering properties in each band on effective radius in
  ! terms of 16 coefficients
  integer, parameter :: nliqopticscoeffssocrates = 16

  ! range of valid input effective radius, in microns
  real(jprb), parameter :: mineffectiveradius = 1.2e-6
  real(jprb), parameter :: maxeffectiveradius = 50.0e-6

contains

  !---------------------------------------------------------------------
  ! compute liquid-droplet scattering properties using a
  ! parameterization consisting of pade approximants from the
  ! socrates (edwards-slingo) code
  subroutine calc_liq_optics_socrates(nb, coeff, lwp, re_in, od, scat_od, g)

    use parkind1, only : jprb
    !use ecradhook,  only : lhook, dr_hook, jphook

    ! number of bands
    integer, intent(in)  :: nb
    ! coefficients read from a data file
    real(jprb), intent(in) :: coeff(:,:)
    ! liquid water path (kg m-2) and effective radius (m)
    real(jprb), intent(in) :: lwp, re_in
    ! total optical depth, scattering optical depth and asymmetry factor
    real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)

    integer    :: jb
    ! local effective radius (m), after applying bounds
    real(jprb) :: re

    !real(jphook) :: hook_handle

    !$acc routine seq

    !if (lhook) call dr_hook('radiation_liquid_optics_socrates:calc_liq_optics_socrates',0,hook_handle)

    ! apply the bounds of validity to effective radius
    re = max(mineffectiveradius, min(re_in, maxeffectiveradius))

! added for dwd (2020)
!nec$ shortloop
    do jb = 1, nb
      od(jb) = lwp * (coeff(jb,1) + re*(coeff(jb,2) + re*coeff(jb,3))) &
         &  / (1.0_jprb + re*(coeff(jb,4) + re*(coeff(jb,5) &
         &  + re*coeff(jb,6))))
      scat_od(jb) = od(jb) * (1.0_jprb &
         &  - (coeff(jb,7) + re*(coeff(jb,8) + re*coeff(jb,9))) &
         &  / (1.0_jprb + re*(coeff(jb,10) + re*coeff(jb,11))))
      g(jb) = (coeff(jb,12) + re*(coeff(jb,13) + re*coeff(jb,14))) &
         &  / (1.0_jprb + re*(coeff(jb,15) + re*coeff(jb,16)))
    end do

    !if (lhook) call dr_hook('radiation_liquid_optics_socrates:calc_liq_optics_socrates',1,hook_handle)

  end subroutine calc_liq_optics_socrates

end module radiation_liquid_optics_socrates
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

