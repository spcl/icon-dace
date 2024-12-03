! # 1 "radiation/radiation_cloudless_lw.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_cloudless_lw.f90"
! radiation_cloudless_lw.f90 - longwave homogeneous cloudless solver
!
! (c) copyright 2019- ecmwf.
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

module radiation_cloudless_lw

public :: solver_cloudless_lw

contains

  !---------------------------------------------------------------------
  ! longwave homogeneous solver containing no clouds
  subroutine solver_cloudless_lw(nlev,istartcol,iendcol, &
       &  config, od, ssa, g, planck_hl, emission, albedo, flux)

    use parkind1, only           : jprb
    use ecradhook,  only           : lhook, dr_hook, jphook

    use radiation_config, only         : config_type
    use radiation_flux, only           : flux_type, indexed_sum_profile
    use radiation_two_stream, only     : calc_two_stream_gammas_lw, &
         &                               calc_reflectance_transmittance_lw, &
         &                               calc_no_scattering_transmittance_lw
    use radiation_adding_ica_lw, only  : adding_ica_lw, calc_fluxes_no_scattering_lw
    use radiation_lw_derivatives, only : calc_lw_derivatives_ica

    implicit none

    ! inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config

    ! gas and aerosol optical depth, single-scattering albedo and
    ! asymmetry factor at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: &
         &  od
    real(jprb), intent(in), dimension(config%n_g_lw_if_scattering, nlev, istartcol:iendcol) :: &
         &  ssa, g

    ! planck function at each half-level and the surface
    real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: &
         &  planck_hl
  
    ! emission (planck*emissivity) and albedo (1-emissivity) at the
    ! surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) &
         &  :: emission, albedo

    ! output
    type(flux_type), intent(inout):: flux

    ! local variables

    ! diffuse reflectance and transmittance for each layer in clear
    ! and all skies
    real(jprb), dimension(config%n_g_lw, nlev) :: reflectance, transmittance

    ! emission by a layer into the upwelling or downwelling diffuse
    ! streams, in clear and all skies
    real(jprb), dimension(config%n_g_lw, nlev) :: source_up, source_dn

    ! fluxes per g point
    real(jprb), dimension(config%n_g_lw, nlev+1) :: flux_up, flux_dn

    ! combined optical depth, single scattering albedo and asymmetry
    ! factor
    real(jprb), dimension(config%n_g_lw) :: ssa_total, g_total

    ! two-stream coefficients
    real(jprb), dimension(config%n_g_lw) :: gamma1, gamma2

    ! number of g points
    integer :: ng

    ! loop indices for level and column
    integer :: jlev, jcol

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloudless_lw:solver_cloudless_lw',0,hook_handle)

    ng = config%n_g_lw

    ! loop through columns
    do jcol = istartcol,iendcol

      ! compute the reflectance and transmittance of all layers,
      ! neglecting clouds
      do jlev = 1,nlev
        if (config%do_lw_aerosol_scattering) then
          ! scattering case: first compute clear-sky reflectance,
          ! transmittance etc at each model level
          ssa_total = ssa(:,jlev,jcol)
          g_total   = g(:,jlev,jcol)
          call calc_two_stream_gammas_lw(ng, ssa_total, g_total, &
               &  gamma1, gamma2)
          call calc_reflectance_transmittance_lw(ng, &
               &  od(:,jlev,jcol), gamma1, gamma2, &
               &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1,jcol), &
               &  reflectance(:,jlev), transmittance(:,jlev), &
               &  source_up(:,jlev), source_dn(:,jlev))
        else
          ! non-scattering case: use simpler functions for
          ! transmission and emission
          call calc_no_scattering_transmittance_lw(ng, od(:,jlev,jcol), &
               &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1, jcol), &
               &  transmittance(:,jlev), source_up(:,jlev), source_dn(:,jlev))          
          ! ensure that clear-sky reflectance is zero
          reflectance(:,jlev) = 0.0_jprb
        end if
      end do

      if (config%do_lw_aerosol_scattering) then
        ! then use adding method to compute fluxes
        call adding_ica_lw(ng, nlev, &
             &  reflectance, transmittance, source_up, source_dn, &
             &  emission(:,jcol), albedo(:,jcol), &
             &  flux_up, flux_dn)
      else
        ! simpler down-then-up method to compute fluxes
        call calc_fluxes_no_scattering_lw(ng, nlev, &
             &  transmittance, source_up, source_dn, &
             &  emission(:,jcol), albedo(:,jcol), &
             &  flux_up, flux_dn)
          
      end if

      ! sum over g-points to compute broadband fluxes
      flux%lw_up(jcol,:) = sum(flux_up,1)
      flux%lw_dn(jcol,:) = sum(flux_dn,1)
      ! store surface spectral downwelling fluxes
      flux%lw_dn_surf_g(:,jcol) = flux_dn(:,nlev+1)

      ! save the spectral fluxes if required
      if (config%do_save_spectral_flux) then
        call indexed_sum_profile(flux_up, config%i_spec_from_reordered_g_lw, &
             &                   flux%lw_up_band(:,jcol,:))
        call indexed_sum_profile(flux_dn, config%i_spec_from_reordered_g_lw, &
             &                   flux%lw_dn_band(:,jcol,:))
      end if

      if (config%do_clear) then
        ! clear-sky calculations are equal to all-sky for this solver:
        ! copy fluxes over
        flux%lw_up_clear(jcol,:) = flux%lw_up(jcol,:)
        flux%lw_dn_clear(jcol,:) = flux%lw_dn(jcol,:)
        flux%lw_dn_surf_clear_g(:,jcol) = flux%lw_dn_surf_g(:,jcol)
        if (config%do_save_spectral_flux) then
          flux%lw_up_clear_band(:,jcol,:) = flux%lw_up_band(:,jcol,:)
          flux%lw_dn_clear_band(:,jcol,:) = flux%lw_dn_band(:,jcol,:)
        end if
      end if

      ! compute the longwave derivatives needed by hogan and bozzo
      ! (2015) approximate radiation update scheme
      if (config%do_lw_derivatives) then
        call calc_lw_derivatives_ica(ng, nlev, jcol, transmittance, flux_up(:,nlev+1), &
             &                       flux%lw_derivatives)
       end if

    end do

    if (lhook) call dr_hook('radiation_cloudless_lw:solver_cloudless_lw',1,hook_handle)
    
  end subroutine solver_cloudless_lw

end module radiation_cloudless_lw
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

