! # 1 "radiation/radiation_ice_optics_fu.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_ice_optics_fu.f90"
! this file has been modified for the use in icon

! radiation_ice_optics_fu.f90 - fu's scheme for ice optical properties
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
!   2020-08-10  r. hogan  bounded re to be <= 100um and g to be < 1.0

module radiation_ice_optics_fu

  use parkind1, only : jprb

  implicit none
  public

  ! the number of ice coefficients depends on the parameterization
  integer, parameter :: niceopticscoeffsfusw  = 10
  integer, parameter :: niceopticscoeffsfulw  = 11

  ! limits based on the range of validity of the parameterizations
  real(jprb), parameter :: maxasymmetryfactor = 1.0_jprb - 10.0_jprb*epsilon(1.0_jprb)
  real(jprb), parameter :: maxeffectiveradius = 100.0e-6_jprb ! metres

contains

  !---------------------------------------------------------------------
  ! compute shortwave ice-particle scattering properties using fu
  ! (1996) parameterization.  the asymmetry factor in band 14 goes
  ! larger than one for re > 100.8 um, so we cap re at 100 um.
  ! asymmetry factor is capped at just less than 1 because if it is
  ! exactly 1 then delta-eddington scaling leads to a zero scattering
  ! optical depth and then division by zero.
  subroutine calc_ice_optics_fu_sw(nb, coeff, ice_wp, &
       &  re, od, scat_od, g)

    !use ecradhook,  only : lhook, dr_hook, jphook

    ! number of bands
    integer, intent(in)  :: nb
    ! coefficients read from a data file
    real(jprb), intent(in) :: coeff(:,:)
    ! ice water path (kg m-2)
    real(jprb), intent(in) :: ice_wp
    ! effective radius (m)
    real(jprb), intent(in) :: re
    ! total optical depth, scattering optical depth and asymmetry factor
    real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)

    ! fu's effective diameter (microns) and its inverse
    real(jprb) :: de_um, inv_de_um
    ! ice water path in g m-2
    real (jprb) :: iwp_gm_2

    integer :: jb
    !real(jphook) :: hook_handle

    !$acc routine seq

    !if (lhook) call dr_hook('radiation_ice_optics:calc_ice_optics_fu_sw',0,hook_handle)

    ! convert to effective diameter using the relationship in the ifs
    de_um     = min(re, maxeffectiveradius) * (1.0e6_jprb / 0.64952_jprb)
    inv_de_um = 1.0_jprb / de_um
    iwp_gm_2  = ice_wp * 1000.0_jprb

! added for dwd (2020)
!nec$ shortloop
    do jb = 1, nb
      od(jb) = iwp_gm_2 * (coeff(jb,1) + coeff(jb,2) * inv_de_um)
      scat_od(jb) = od(jb) * (1.0_jprb - (coeff(jb,3) + de_um*(coeff(jb,4) &
         &  + de_um*(coeff(jb,5) + de_um*coeff(jb,6)))))
      g(jb) = min(coeff(jb,7) + de_um*(coeff(jb,8) &
         &  + de_um*(coeff(jb,9) + de_um*coeff(jb,10))), &
         &  maxasymmetryfactor)
    end do

    !if (lhook) call dr_hook('radiation_ice_optics:calc_ice_optics_fu_sw',1,hook_handle)

  end subroutine calc_ice_optics_fu_sw


  !---------------------------------------------------------------------
  ! compute longwave ice-particle scattering properties using fu et
  ! al. (1998) parameterization
  subroutine calc_ice_optics_fu_lw(nb, coeff, ice_wp, &
       &  re, od, scat_od, g)

    !use ecradhook,  only : lhook, dr_hook, jphook

    ! number of bands
    integer, intent(in)  :: nb
    ! coefficients read from a data file
    real(jprb), intent(in) :: coeff(:,:)
    ! ice water path (kg m-2)
    real(jprb), intent(in) :: ice_wp
    ! effective radius (m)
    real(jprb), intent(in) :: re
    ! total optical depth, scattering optical depth and asymmetry factor
    real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)

    ! fu's effective diameter (microns) and its inverse
    real(jprb) :: de_um, inv_de_um
    ! ice water path in g m-2
    real (jprb) :: iwp_gm_2

    integer :: jb
    !real(jphook) :: hook_handle

    !$acc routine seq

    !if (lhook) call dr_hook('radiation_ice_optics:calc_ice_optics_fu_lw',0,hook_handle)

    ! convert to effective diameter using the relationship in the ifs
    de_um = min(re, maxeffectiveradius) * (1.0e6_jprb / 0.64952_jprb)

    inv_de_um = 1.0_jprb / de_um
    iwp_gm_2  = ice_wp * 1000.0_jprb

! added for dwd (2020)
!nec$ shortloop
    do jb = 1, nb
      od(jb) = iwp_gm_2 * (coeff(jb,1) + inv_de_um*(coeff(jb,2) &
         &  + inv_de_um*coeff(jb,3)))
      scat_od(jb) = od(jb) - iwp_gm_2*inv_de_um*(coeff(jb,4) + de_um*(coeff(jb,5) &
         &  + de_um*(coeff(jb,6) + de_um*coeff(jb,7))))
      g(jb) = min(coeff(jb,8) + de_um*(coeff(jb,9) &
         &  + de_um*(coeff(jb,10) + de_um*coeff(jb,11))), &
         &  maxasymmetryfactor)
    end do

    !if (lhook) call dr_hook('radiation_ice_optics:calc_ice_optics_fu_lw',1,hook_handle)

  end subroutine calc_ice_optics_fu_lw

end module radiation_ice_optics_fu
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

