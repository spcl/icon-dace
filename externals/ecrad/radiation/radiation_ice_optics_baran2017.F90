! # 1 "radiation/radiation_ice_optics_baran2017.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_ice_optics_baran2017.f90"
! radiation_ice_optics_baran2017.f90 - 2017 parameterization of baran's ice optical properties
!
! (c) copyright 2017- ecmwf.
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

module radiation_ice_optics_baran2017

  implicit none
  public

  ! the number of ice coefficients depends on the parameterization
  integer, parameter :: niceopticscoeffsbaran2017 = 9
  integer, parameter :: niceopticsgeneralcoeffsbaran2017 = 5

contains

  
  !---------------------------------------------------------------------
  ! compute ice-particle scattering properties using a
  ! parameterization as a function of ice water mixing ratio and
  ! temperature
  subroutine calc_ice_optics_baran2017(nb, coeff_gen, coeff, ice_wp, &
       &  qi, temperature, od, scat_od, g)

    use parkind1, only : jprb
    !use ecradhook,  only : lhook, dr_hook, jphook

    ! number of bands
    integer, intent(in)  :: nb
    ! general coefficients read from a data file
    real(jprb), intent(in) :: coeff_gen(:)
    ! band-specific coefficients read from a data file
    real(jprb), intent(in) :: coeff(:,:)
    ! ice water path (kg m-2) and mixing ratio (kg kg-1)
    real(jprb), intent(in) :: ice_wp, qi
    ! temperature (k)
    real(jprb), intent(in) :: temperature
    ! total optical depth, scattering optical depth and asymmetry factor
    real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)
    
    ! modified ice mixing ratio, and the same raised to an appropriate power
    real(jprb) :: qi_mod, qi_mod_od, qi_mod_ssa, qi_mod_g
    
    !real(jphook) :: hook_handle

    !if (lhook) call dr_hook('radiation_ice_optics:calc_ice_optics_baran2017',0,hook_handle)

    qi_mod     = qi * exp(coeff_gen(1)*(temperature-coeff_gen(2)))
    qi_mod_od  = qi_mod ** coeff_gen(3)
    qi_mod_ssa = qi_mod ** coeff_gen(4)
    qi_mod_g   = qi_mod ** coeff_gen(5)

    od      = ice_wp * (coeff(1:nb,1) + coeff(1:nb,2)/(1.0_jprb+qi_mod_od *coeff(1:nb,3)))
    scat_od = od     * (coeff(1:nb,4) + coeff(1:nb,5)/(1.0_jprb+qi_mod_ssa*coeff(1:nb,6)))
    g       =           coeff(1:nb,7) + coeff(1:nb,8)/(1.0_jprb+qi_mod_g  *coeff(1:nb,9))

    !if (lhook) call dr_hook('radiation_ice_optics:calc_ice_optics_baran2017',1,hook_handle)

  end subroutine calc_ice_optics_baran2017

end module radiation_ice_optics_baran2017
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

