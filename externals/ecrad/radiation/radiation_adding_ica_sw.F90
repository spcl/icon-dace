! # 1 "radiation/radiation_adding_ica_sw.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_adding_ica_sw.f90"
! this file has been modified for the use in icon

! radiation_adding_ica_sw.f90 - shortwave adding method in independent column approximation
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
! modifications
!   2017-10-23  r. hogan  renamed single-character variables

module radiation_adding_ica_sw

  public

contains

  subroutine adding_ica_sw(ncol, nlev, incoming_toa, &
       &  albedo_surf_diffuse, albedo_surf_direct, cos_sza, &
       &  reflectance, transmittance, ref_dir, trans_dir_diff, trans_dir_dir, &
       &  flux_up, flux_dn_diffuse, flux_dn_direct, &
       &  albedo, source, inv_denominator)

    use parkind1, only           : jprb
    use ecradhook,  only           : lhook, dr_hook, jphook

    implicit none

    ! inputs
    integer, intent(in) :: ncol ! number of columns (may be spectral intervals)
    integer, intent(in) :: nlev ! number of levels

    ! incoming downwelling solar radiation at top-of-atmosphere (w m-2)
    real(jprb), intent(in),  dimension(ncol)         :: incoming_toa

    ! surface albedo to diffuse and direct radiation
    real(jprb), intent(in),  dimension(ncol)         :: albedo_surf_diffuse, &
         &                                              albedo_surf_direct

    ! cosine of the solar zenith angle
    real(jprb), intent(in)                           :: cos_sza

    ! diffuse reflectance and transmittance of each layer
    real(jprb), intent(in),  dimension(ncol, nlev)   :: reflectance, transmittance

    ! fraction of direct-beam solar radiation entering the top of a
    ! layer that is reflected back up or scattered forward into the
    ! diffuse stream at the base of the layer
    real(jprb), intent(in),  dimension(ncol, nlev)   :: ref_dir, trans_dir_diff

    ! direct transmittance, i.e. fraction of direct beam that
    ! penetrates a layer without being scattered or absorbed
    real(jprb), intent(in),  dimension(ncol, nlev)   :: trans_dir_dir

    ! resulting fluxes (w m-2) at half-levels: diffuse upwelling,
    ! diffuse downwelling and direct downwelling
    real(jprb), intent(out), dimension(ncol, nlev+1) :: flux_up, flux_dn_diffuse, &
         &                                              flux_dn_direct
    
    ! albedo of the entire earth/atmosphere system below each half
    ! level
    real(jprb), intent(out), dimension(ncol, nlev+1) :: albedo

    ! upwelling radiation at each half-level due to scattering of the
    ! direct beam below that half-level (w m-2)
    real(jprb), intent(out), dimension(ncol, nlev+1) :: source

    ! equal to 1/(1-albedo*reflectance)
    real(jprb), intent(out), dimension(ncol, nlev)   :: inv_denominator

    ! loop index for model level and column
    integer :: jlev, jcol

    real(jphook) :: hook_handle


    if (lhook) call dr_hook('radiation_adding_ica_sw:adding_ica_sw',0,hook_handle)


    !$acc routine worker

    ! compute profile of direct (unscattered) solar fluxes at each
    ! half-level by working down through the atmosphere
    flux_dn_direct(:,1) = incoming_toa
    !$acc loop seq
    do jlev = 1,nlev
      flux_dn_direct(:,jlev+1) = flux_dn_direct(:,jlev)*trans_dir_dir(:,jlev)
    end do

    albedo(:,nlev+1) = albedo_surf_diffuse

    ! at the surface, the direct solar beam is reflected back into the
    ! diffuse stream
    source(:,nlev+1) = albedo_surf_direct * flux_dn_direct(:,nlev+1) * cos_sza

    ! work back up through the atmosphere and compute the albedo of
    ! the entire earth/atmosphere system below that half-level, and
    ! also the "source", which is the upwelling flux due to direct
    ! radiation that is scattered below that level
!$acc loop seq
! added for dwd (2020)
!nec$ outerloop_unroll(8)
    do jlev = nlev,1,-1
      ! next loop over columns. we could do this by indexing the
      ! entire inner dimension as follows, e.g. for the first line:
      !   inv_denominator(:,jlev) = 1.0_jprb / (1.0_jprb-albedo(:,jlev+1)*reflectance(:,jlev))
      ! and similarly for subsequent lines, but this slows down the
      ! routine by a factor of 2!  rather, we do it with an explicit
      ! loop.
      !$acc loop worker vector
      do jcol = 1,ncol
        ! lacis and hansen (1974) eq 33, shonk & hogan (2008) eq 10:
        inv_denominator(jcol,jlev) = 1.0_jprb / (1.0_jprb-albedo(jcol,jlev+1)*reflectance(jcol,jlev))
        ! shonk & hogan (2008) eq 9, petty (2006) eq 13.81:
        albedo(jcol,jlev) = reflectance(jcol,jlev) + transmittance(jcol,jlev) * transmittance(jcol,jlev) &
             &                                     * albedo(jcol,jlev+1) * inv_denominator(jcol,jlev)
        ! shonk & hogan (2008) eq 11:
        source(jcol,jlev) = ref_dir(jcol,jlev)*flux_dn_direct(jcol,jlev) &
             &  + transmittance(jcol,jlev)*(source(jcol,jlev+1) &
             &        + albedo(jcol,jlev+1)*trans_dir_diff(jcol,jlev)*flux_dn_direct(jcol,jlev)) &
             &  * inv_denominator(jcol,jlev)
      end do
    end do

    ! at top-of-atmosphere there is no diffuse downwelling radiation
    flux_dn_diffuse(:,1) = 0.0_jprb

    ! at top-of-atmosphere, all upwelling radiation is due to
    ! scattering by the direct beam below that level
    flux_up(:,1) = source(:,1)

    ! work back down through the atmosphere computing the fluxes at
    ! each half-level
!$acc loop seq
! added for dwd (2020)
!nec$ outerloop_unroll(8)
    do jlev = 1,nlev
      !$acc loop worker vector
      do jcol = 1,ncol
        ! shonk & hogan (2008) eq 14 (after simplification):
        flux_dn_diffuse(jcol,jlev+1) &
             &  = (transmittance(jcol,jlev)*flux_dn_diffuse(jcol,jlev) &
             &     + reflectance(jcol,jlev)*source(jcol,jlev+1) &
             &     + trans_dir_diff(jcol,jlev)*flux_dn_direct(jcol,jlev)) * inv_denominator(jcol,jlev)
        ! shonk & hogan (2008) eq 12:
        flux_up(jcol,jlev+1) = albedo(jcol,jlev+1)*flux_dn_diffuse(jcol,jlev+1) &
             &            + source(jcol,jlev+1)
        flux_dn_direct(jcol,jlev) = flux_dn_direct(jcol,jlev)*cos_sza
      end do
    end do
    flux_dn_direct(:,nlev+1) = flux_dn_direct(:,nlev+1)*cos_sza


    if (lhook) call dr_hook('radiation_adding_ica_sw:adding_ica_sw',1,hook_handle)


  end subroutine adding_ica_sw

end module radiation_adding_ica_sw
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

