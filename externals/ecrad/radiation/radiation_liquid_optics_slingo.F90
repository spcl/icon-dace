! # 1 "radiation/radiation_liquid_optics_slingo.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_liquid_optics_slingo.f90"
! radiation_liquid_optics_slingo.f90 - slingo sw & lindner-li lw parameterization of liquid droplet optics
!
! (c) copyright 2016- ecmwf.
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

module radiation_liquid_optics_slingo

  implicit none
  public

  integer, parameter :: nliqopticscoeffsslingosw = 6
  integer, parameter :: nliqopticscoeffslindnerlilw = 13

contains

  !---------------------------------------------------------------------
  ! compute liquid-droplet scattering properties in the shortwave from
  ! slingo (1989). warning: this parameterization is known not to be
  ! very accurate: see nielsen et al. (gmd 2014).
  subroutine calc_liq_optics_slingo(nb, coeff, lwp, re, od, scat_od, g)

    use parkind1, only : jprb
    !use ecradhook,  only : lhook, dr_hook, jphook

    ! number of bands
    integer, intent(in)  :: nb
    ! coefficients read from a data file
    real(jprb), intent(in) :: coeff(:,:)
    ! liquid water path (kg m-2) and effective radius (m)
    real(jprb), intent(in) :: lwp, re
    ! total optical depth, scattering optical depth and asymmetry factor
    real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)

    ! liquid water path in g m-2
    real(jprb) :: lwp_gm_2
    ! effective radius in microns, and its inverse
    real(jprb) :: re_um, inv_re_um

    !real(jphook) :: hook_handle

    !if (lhook) call dr_hook('radiation_liquid_optics_slingo:calc_liq_optics_slingo',0,hook_handle)

    lwp_gm_2 = lwp * 1000.0_jprb
    ! range of validity reported by slingo (1989): 4.2-16.6 microns
    re_um = min(max(4.2_jprb, re * 1.0e6_jprb), 16.6_jprb)
    inv_re_um = 1.0_jprb / re_um

    od = lwp_gm_2 * (coeff(1:nb,1) + inv_re_um*coeff(1:nb,2))
    scat_od = od * (1.0_jprb - coeff(1:nb,3) - re_um*coeff(1:nb,4))
    g = coeff(1:nb,5) + re_um*coeff(1:nb,6)

    !if (lhook) call dr_hook('radiation_liquid_optics_slingo:calc_liq_optics_slingo',1,hook_handle)

  end subroutine calc_liq_optics_slingo


  !---------------------------------------------------------------------
  ! compute liquid-droplet scattering properties in the longwave from
  ! lindner & li (2000)
  subroutine calc_liq_optics_lindner_li(nb, coeff, lwp, re, od, scat_od, g)

    use parkind1, only : jprb
    use ecradhook,  only : lhook, dr_hook, jphook

    ! number of bands
    integer, intent(in)  :: nb
    ! coefficients read from a data file
    real(jprb), intent(in) :: coeff(:,:)
    ! liquid water path (kg m-2) and effective radius (m)
    real(jprb), intent(in) :: lwp, re
    ! total optical depth, scattering optical depth and asymmetry factor
    real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)

    ! liquid water path in g m-2
    real(jprb) :: lwp_gm_2
    ! effective radius in microns, and its inverse
    real(jprb) :: re_um, inv_re_um

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_liquid_optics_slingo:calc_liq_optics_lindner_li',0,hook_handle)

    lwp_gm_2 = lwp * 1000.0_jprb
    ! range of validity reported by lindner and li (2000): 2-40 microns
    re_um = min(max(2.0_jprb, re * 1.0e6_jprb), 40.0_jprb)
    inv_re_um = 1.0_jprb / re_um

    od = lwp_gm_2 * (coeff(1:nb,1) + re_um*coeff(1:nb,2) + inv_re_um*(coeff(1:nb,3) &
         &  + inv_re_um*(coeff(1:nb,4) + inv_re_um*coeff(1:nb,5))))
    scat_od = od * (1.0_jprb - (coeff(1:nb,6) + inv_re_um*coeff(1:nb,7) &
         &                      + re_um*(coeff(1:nb,8) + re_um*coeff(1:nb,9))))
    g = coeff(1:nb,10) + inv_re_um*coeff(1:nb,11) &
         &  + re_um*(coeff(1:nb,12) + re_um*coeff(1:nb,13))

    if (lhook) call dr_hook('radiation_liquid_optics_slingo:calc_liq_optics_lindner_li',1,hook_handle)

  end subroutine calc_liq_optics_lindner_li

end module radiation_liquid_optics_slingo
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

