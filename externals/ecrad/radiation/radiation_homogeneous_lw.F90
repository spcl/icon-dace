! # 1 "radiation/radiation_homogeneous_lw.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_homogeneous_lw.f90"
! radiation_homogeneous_lw.f90 - longwave homogeneous-column (no cloud fraction) solver
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
! modifications
!   2017-04-11  r. hogan  receive emission/albedo rather than planck/emissivity
!   2017-04-22  r. hogan  store surface fluxes at all g-points
!   2017-10-23  r. hogan  renamed single-character variables
!   2019-01-14  r. hogan  save spectral flux profile if required

module radiation_homogeneous_lw

  public

contains

  !---------------------------------------------------------------------
  ! longwave homogeneous solver, in which clouds are assumed to fill
  ! the gridbox horizontally
  subroutine solver_homogeneous_lw(nlev,istartcol,iendcol, &
       &  config, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
       &  emission, albedo, &
       &  flux)

    use parkind1, only           : jprb
    use ecradhook,  only           : lhook, dr_hook, jphook

    use radiation_config, only         : config_type
    use radiation_cloud, only          : cloud_type
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
    type(cloud_type),         intent(in) :: cloud

    ! gas and aerosol optical depth, single-scattering albedo and
    ! asymmetry factor at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: &
         &  od
    real(jprb), intent(in), dimension(config%n_g_lw_if_scattering, nlev, istartcol:iendcol) :: &
         &  ssa, g

    ! cloud and precipitation optical depth, single-scattering albedo and
    ! asymmetry factor in each longwave band
    real(jprb), intent(in), dimension(config%n_bands_lw,nlev,istartcol:iendcol)   :: &
         &  od_cloud
    real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, &
         &  nlev,istartcol:iendcol) :: ssa_cloud, g_cloud

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

    ! combined gas+aerosol+cloud optical depth, single scattering
    ! albedo and asymmetry factor
    real(jprb), dimension(config%n_g_lw) :: od_total, ssa_total, g_total

    ! two-stream coefficients
    real(jprb), dimension(config%n_g_lw) :: gamma1, gamma2

    ! optical depth of cloud in g-point space
    real(jprb), dimension(config%n_g_lw) :: od_cloud_g

    ! is there any cloud in the profile?
    logical :: is_cloudy_profile

    ! number of g points
    integer :: ng

    ! loop indices for level and column
    integer :: jlev, jcol

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_homogeneous_lw:solver_homogeneous_lw',0,hook_handle)

    ng = config%n_g_lw

    ! loop through columns
    do jcol = istartcol,iendcol

      ! is there any cloud in the profile?
      is_cloudy_profile = .false.
      do jlev = 1,nlev
        if (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold) then
          is_cloudy_profile = .true.
          exit
        end if
      end do

      ! if clear-sky fluxes need to be computed then we first compute
      ! the reflectance and transmittance of all layers, neglecting
      ! clouds. if clear-sky fluxes are not required then we only do
      ! the clear-sky layers since these will be needed when we come
      ! to do the total-sky fluxes.
      do jlev = 1,nlev
        if (config%do_clear .or. cloud%fraction(jcol,jlev) &
             &                 < config%cloud_fraction_threshold) then
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
            ! ensure that clear-sky reflectance is zero since it may be
            ! used in cloudy-sky case
            reflectance(:,jlev) = 0.0_jprb
          end if
         
        end if
      end do

      if (config%do_clear) then
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
        flux%lw_up_clear(jcol,:) = sum(flux_up,1)
        flux%lw_dn_clear(jcol,:) = sum(flux_dn,1)
        ! store surface spectral downwelling fluxes
        flux%lw_dn_surf_clear_g(:,jcol) = flux_dn(:,nlev+1)

        ! save the spectral fluxes if required
        if (config%do_save_spectral_flux) then
          call indexed_sum_profile(flux_up, config%i_spec_from_reordered_g_lw, &
               &                   flux%lw_up_clear_band(:,jcol,:))
          call indexed_sum_profile(flux_dn, config%i_spec_from_reordered_g_lw, &
               &                   flux%lw_dn_clear_band(:,jcol,:))
        end if

      end if ! do clear-sky calculations

      ! now the total-sky calculation.  if this is a clear profile and
      ! clear-sky fluxes have been calculated then we can simply copy
      ! over the clear-sky fluxes, otherwise we need to compute fluxes
      ! now.
      if (is_cloudy_profile .or. .not. config%do_clear) then
        do jlev = 1,nlev
          ! compute combined gas+aerosol+cloud optical properties;
          ! note that for clear layers, the reflectance and
          ! transmittance have already been calculated
          if (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold) then
            od_cloud_g = od_cloud(config%i_band_from_reordered_g_lw,jlev,jcol)
            od_total = od(:,jlev,jcol) + od_cloud_g
            ssa_total = 0.0_jprb
            g_total   = 0.0_jprb

            if (config%do_lw_cloud_scattering) then
              ! scattering case: calculate reflectance and
              ! transmittance at each model level
              if (config%do_lw_aerosol_scattering) then
                where (od_total > 0.0_jprb)
                  ssa_total = (ssa(:,jlev,jcol)*od(:,jlev,jcol) &
                       &     + ssa_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     *  od_cloud_g) & 
                       &     / od_total
                end where
                where (ssa_total > 0.0_jprb .and. od_total > 0.0_jprb)
                  g_total = (g(:,jlev,jcol)*ssa(:,jlev,jcol)*od(:,jlev,jcol) &
                       &     +   g_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     * ssa_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     *  od_cloud_g) &
                       &     / (ssa_total*od_total)
                end where
              else
                where (od_total > 0.0_jprb)
                  ssa_total = ssa_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     * od_cloud_g / od_total
                end where
                where (ssa_total > 0.0_jprb .and. od_total > 0.0_jprb)
                  g_total = g_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     * ssa_cloud(config%i_band_from_reordered_g_lw,jlev,jcol) &
                       &     *  od_cloud_g / (ssa_total*od_total)
                end where
              end if
            
              ! compute cloudy-sky reflectance, transmittance etc at
              ! each model level
              call calc_two_stream_gammas_lw(ng, ssa_total, g_total, &
                   &  gamma1, gamma2)
              call calc_reflectance_transmittance_lw(ng, &
                   &  od_total, gamma1, gamma2, &
                   &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1,jcol), &
                   &  reflectance(:,jlev), transmittance(:,jlev), &
                   &  source_up(:,jlev), source_dn(:,jlev))
            else
              ! no-scattering case: use simpler functions for
              ! transmission and emission
              call calc_no_scattering_transmittance_lw(ng, od_total, &
                   &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1, jcol), &
                   &  transmittance(:,jlev), source_up(:,jlev), source_dn(:,jlev))
            end if
          end if ! is cloudy layer
        end do
        
        if (config%do_lw_cloud_scattering) then
          ! use adding method to compute fluxes for an overcast sky
          call adding_ica_lw(ng, nlev, reflectance, transmittance, source_up, source_dn, &
               &  emission(:,jcol), albedo(:,jcol), &
               &  flux_up, flux_dn)
        else
          ! simpler down-then-up method to compute fluxes
          call calc_fluxes_no_scattering_lw(ng, nlev, &
               &  transmittance, source_up, source_dn, emission(:,jcol), albedo(:,jcol), &
               &  flux_up, flux_dn)
        end if
        
        ! store overcast broadband fluxes
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

      else
        ! no cloud in profile and clear-sky fluxes already
        ! calculated: copy them over
        flux%lw_up(jcol,:) = flux%lw_up_clear(jcol,:)
        flux%lw_dn(jcol,:) = flux%lw_dn_clear(jcol,:)
        flux%lw_dn_surf_g(:,jcol) = flux%lw_dn_surf_clear_g(:,jcol)

        if (config%do_save_spectral_flux) then
          flux%lw_up_band(:,jcol,:) = flux%lw_up_clear_band(:,jcol,:)
          flux%lw_dn_band(:,jcol,:) = flux%lw_dn_clear_band(:,jcol,:)
        end if

     end if

      ! compute the longwave derivatives needed by hogan and bozzo
      ! (2015) approximate radiation update scheme, using clear-sky
      ! transmittance if no clouds were present in the profile,
      ! all-sky transmittance otherwise
      if (config%do_lw_derivatives) then
        call calc_lw_derivatives_ica(ng, nlev, jcol, transmittance, flux_up(:,nlev+1), &
             &                       flux%lw_derivatives)
 
      end if

    end do

    if (lhook) call dr_hook('radiation_homogeneous_lw:solver_homogeneous_lw',1,hook_handle)
    
  end subroutine solver_homogeneous_lw

end module radiation_homogeneous_lw
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

