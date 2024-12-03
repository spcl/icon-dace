! # 1 "radiation/radiation_ice_optics_yi.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_ice_optics_yi.f90"
! radiation_ice_optics_yi.f90 - yi et al. (2013) ice optical properties
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
! authors: mark fielding and robin hogan
! email:   r.j.hogan@ecmwf.int
!
! the reference for this ice optics parameterization is yi, b.,
! p. yang, b.a. baum, t. l'ecuyer, l. oreopoulos, e.j. mlawer,
! a.j. heymsfield, and k. liou, 2013: influence of ice particle
! surface roughening on the global cloud radiative
! effect. j. atmos. sci., 70, 2794-2807,
! https://doi.org/10.1175/jas-d-13-020.1

module radiation_ice_optics_yi

  implicit none
  public

  ! the number of ice coefficients depends on the parameterization
  integer, parameter :: niceopticscoeffsyisw  = 69
  integer, parameter :: niceopticscoeffsyilw  = 69

  integer, parameter :: nsinglecoeffs = 23

contains

  !---------------------------------------------------------------------
  ! compute shortwave ice-particle scattering properties using yi et
  ! al. (2013) parameterization
  subroutine calc_ice_optics_yi_sw(nb, coeff, ice_wp, &
       &  re, od, scat_od, g)

    use parkind1, only : jprb, jpim
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

    ! yi's effective diameter (microns)
    real(jprb) :: de_um
    ! ice water path in g m-2
    real (jprb) :: iwp_gm_2
    ! lut temp variables
    real(jprb) :: wts_1, wts_2
    integer(jpim) :: lu_idx
    real(kind=jprb), parameter    :: lu_scale  = 0.2_jprb
    real(kind=jprb), parameter    :: lu_offset = 1.0_jprb
    !real(jphook) :: hook_handle

    !if (lhook) call dr_hook('radiation_ice_optics:calc_ice_optics_yi_sw',0,hook_handle)

    ! convert to effective diameter using the relationship in the ifs
    !de_um     = re * (1.0e6_jprb / 0.64952_jprb)
    de_um     = re * 2.0e6_jprb

    ! limit de_um to validity of lut
    de_um = max(de_um,10.0_jprb)
    de_um = min(de_um,119.99_jprb) !avoid greater than or equal to 120 um

    iwp_gm_2  = ice_wp * 1000.0_jprb

    lu_idx = floor(de_um * lu_scale - lu_offset)
    wts_2  = (de_um * lu_scale - lu_offset) - lu_idx
    wts_1  = 1.0_jprb - wts_2
    od     = 0.001_jprb * iwp_gm_2 * & 
             & ( wts_1 * coeff(1:nb,lu_idx) + wts_2 * coeff(1:nb,lu_idx+1) )
    scat_od = od * & 
             & ( wts_1 * coeff(1:nb,lu_idx+nsinglecoeffs) + wts_2 * coeff(1:nb,lu_idx+nsinglecoeffs+1) )
    g = wts_1 * coeff(1:nb,lu_idx+2*nsinglecoeffs) + wts_2 * coeff(1:nb,lu_idx+2*nsinglecoeffs+1)

    !if (lhook) call dr_hook('radiation_ice_optics:calc_ice_optics_yi_sw',1,hook_handle)

  end subroutine calc_ice_optics_yi_sw


  !---------------------------------------------------------------------
  ! compute longwave ice-particle scattering properties using yi et
  ! al. (2013) parameterization
  subroutine calc_ice_optics_yi_lw(nb, coeff, ice_wp, &
       &  re, od, scat_od, g)

    use parkind1, only : jprb, jpim
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

    ! yi's effective diameter (microns)
    real(jprb) :: de_um
    ! ice water path in g m-2
    real (jprb) :: iwp_gm_2
    ! lut temp variables
    real(jprb) :: wts_1, wts_2
    integer(jpim) :: lu_idx
    real(kind=jprb), parameter    :: lu_scale  = 0.2_jprb
    real(kind=jprb), parameter    :: lu_offset = 1.0_jprb
    !real(jphook) :: hook_handle

    !if (lhook) call dr_hook('radiation_ice_optics:calc_ice_optics_yi_sw',0,hook_handle)

    ! convert to effective diameter using the relationship in the ifs
    !de_um     = re * (1.0e6_jprb / 0.64952_jprb)
    de_um     = re * 2.0e6_jprb

    ! limit de_um to validity of lut
    de_um = max(de_um,10.0_jprb)
    de_um = min(de_um,119.99_jprb) !avoid greater than or equal to 120 um

    iwp_gm_2  = ice_wp * 1000.0_jprb

    lu_idx = floor(de_um * lu_scale - lu_offset)
    wts_2  = (de_um * lu_scale - lu_offset) - lu_idx
    wts_1  = 1.0_jprb - wts_2
    od     = 0.001_jprb * iwp_gm_2 * & 
             & ( wts_1 * coeff(1:nb,lu_idx) + wts_2 * coeff(1:nb,lu_idx+1) )
    scat_od = od * & 
             & ( wts_1 * coeff(1:nb,lu_idx+nsinglecoeffs) + wts_2 * coeff(1:nb,lu_idx+nsinglecoeffs+1) )
    g = wts_1 * coeff(1:nb,lu_idx+2*nsinglecoeffs) + wts_2 * coeff(1:nb,lu_idx+2*nsinglecoeffs+1)

     !if (lhook) call dr_hook('radiation_ice_optics:calc_ice_optics_yi_lw',1,hook_handle)

  end subroutine calc_ice_optics_yi_lw

end module radiation_ice_optics_yi
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

