! # 1 "radiation/radiation_cloudless_sw.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_cloudless_sw.f90"
! this file has been modified for the use in icon

! radiation_cloudless_sw.f90 - shortwave homogeneous cloudless solver
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

module radiation_cloudless_sw

public :: solver_cloudless_sw

contains

  ! provides elemental function "delta_eddington"

! # 1 "radiation/radiation_delta_eddington.h" 1
! radiation_delta_eddington.h - delta-eddington scaling -*- f90 -*-
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
! this file is intended to be included inside a module to ensure that
! these simple functions may be inlined

!---------------------------------------------------------------------
! perform in-place delta-eddington scaling of the phase function
elemental subroutine delta_eddington(od, ssa, g)

  use parkind1, only : jprb
  
  ! total optical depth, single scattering albedo and asymmetry
  ! factor
  real(jprb), intent(inout) :: od, ssa, g
  
  ! fraction of the phase function deemed to be in the forward lobe
  ! and therefore treated as if it is not scattered at all
  real(jprb) :: f
  
  f   = g*g
  od  = od * (1.0_jprb - ssa*f)
  ssa = ssa * (1.0_jprb - f) / (1.0_jprb - ssa*f)
  g   = g / (1.0_jprb + g)
  
end subroutine delta_eddington


!---------------------------------------------------------------------
! perform in-place delta-eddington scaling of the phase function, but
! using extensive variables (i.e. the scattering optical depth,
! scat_od, rather than the single-scattering albedo, and the
! scattering-optical-depth-multiplied-by-asymmetry-factor, scat_od_g,
! rather than the asymmetry factor.
elemental subroutine delta_eddington_extensive(od, scat_od, scat_od_g)

  !$acc routine seq

  use parkind1, only : jprb

  ! total optical depth, scattering optical depth and asymmetry factor
  ! multiplied by the scattering optical depth
  real(jprb), intent(inout) :: od, scat_od, scat_od_g

  ! fraction of the phase function deemed to be in the forward lobe
  ! and therefore treated as if it is not scattered at all
  real(jprb) :: f, g

  if (scat_od > 0.0_jprb) then
    g = scat_od_g / scat_od
  else
    g = 0.0
  end if

  f         = g*g
  od        = od - scat_od * f
  scat_od   = scat_od * (1.0_jprb - f)
  scat_od_g = scat_od * g / (1.0_jprb + g)
  
end subroutine delta_eddington_extensive


!---------------------------------------------------------------------
! array version of delta_eddington_extensive, more likely to vectorize
 subroutine delta_eddington_extensive_vec(ng, od, scat_od, scat_od_g)

  use parkind1, only : jprb

  ! total optical depth, scattering optical depth and asymmetry factor
  ! multiplied by the scattering optical depth
  integer,                   intent(in)    :: ng
  real(jprb), dimension(ng), intent(inout) :: od, scat_od, scat_od_g

  ! fraction of the phase function deemed to be in the forward lobe
  ! and therefore treated as if it is not scattered at all
  real(jprb) :: f, g
  integer :: j

  do j = 1,ng
    g            = scat_od_g(j) / max(scat_od(j), 1.0e-24)
    f            = g*g
    od(j)        = od(j) - scat_od(j) * f
    scat_od(j)   = scat_od(j) * (1.0_jprb - f)
    scat_od_g(j) = scat_od(j) * g / (1.0_jprb + g)
  end do
  
end subroutine delta_eddington_extensive_vec


!---------------------------------------------------------------------
! perform in-place delta-eddington scaling of the phase function,
! using the scattering optical depth rather than the single scattering
! albedo
elemental subroutine delta_eddington_scat_od(od, scat_od, g)

  use parkind1, only : jprb
  
  ! total optical depth, scattering optical depth and asymmetry factor
  real(jprb), intent(inout) :: od, scat_od, g

  ! fraction of the phase function deemed to be in the forward lobe
  ! and therefore treated as if it is not scattered at all
  real(jprb) :: f

  !$acc routine seq

  f       = g*g
  od      = od - scat_od * f
  scat_od = scat_od * (1.0_jprb - f)
  g       = g / (1.0_jprb + g)

end subroutine delta_eddington_scat_od


!---------------------------------------------------------------------
! revert delta-eddington-scaled quantities in-place, back to their
! original state
elemental subroutine revert_delta_eddington(od, ssa, g)

  use parkind1, only : jprb
  
  ! total optical depth, single scattering albedo and asymmetry
  ! factor
  real(jprb), intent(inout) :: od, ssa, g
  
  ! fraction of the phase function deemed to be in the forward lobe
  ! and therefore treated as if it is not scattered at all
  real(jprb) :: f
  
  g   = g / (1.0_jprb - g)
  f   = g*g
  ssa = ssa / (1.0_jprb - f + f*ssa);
  od  = od / (1.0_jprb - ssa*f)
  
end subroutine revert_delta_eddington
! # 26 "radiation/radiation_cloudless_sw.f90" 2

  !---------------------------------------------------------------------
  ! shortwave homogeneous solver containing no clouds
  subroutine solver_cloudless_sw(nlev,istartcol,iendcol, &
       &  config, single_level, & 
       &  od, ssa, g, albedo_direct, albedo_diffuse, incoming_sw, &
       &  flux)

    use parkind1, only           : jprb
    use ecradhook,  only           : lhook, dr_hook, jphook

    use radiation_config, only         : config_type
    use radiation_single_level, only   : single_level_type
    use radiation_flux, only           : flux_type, indexed_sum_profile, &
         &                               add_indexed_sum_profile
    use radiation_two_stream, only     : calc_two_stream_gammas_sw, &
         &                       calc_reflectance_transmittance_sw
    use radiation_constants, only      : pi, gasconstantdryair, &
         &                               accelduetogravity
    use radiation_adding_ica_sw, only  : adding_ica_sw

    implicit none

    ! inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level

    ! gas and aerosol optical depth, single-scattering albedo and
    ! asymmetry factor at each shortwave g-point
    real(jprb), intent(in), dimension(config%n_g_sw, nlev, istartcol:iendcol) :: &
         &  od, ssa, g

    ! direct and diffuse surface albedos, and the incoming shortwave
    ! flux into a plane perpendicular to the incoming radiation at
    ! top-of-atmosphere in each of the shortwave g points
    real(jprb), intent(in), dimension(config%n_g_sw,istartcol:iendcol) :: &
         &  albedo_direct, albedo_diffuse, incoming_sw

    ! output
    type(flux_type), intent(inout):: flux

    ! local variables

    ! cosine of solar zenith angle
    real(jprb)                                 :: cos_sza

    ! diffuse reflectance and transmittance for each layer
    real(jprb), dimension(config%n_g_sw, nlev) :: reflectance, transmittance

    ! fraction of direct beam scattered by a layer into the upwelling
    ! or downwelling diffuse streams
    real(jprb), dimension(config%n_g_sw, nlev) :: ref_dir, trans_dir_diff

    ! transmittance for the direct beam in clear and all skies
    real(jprb), dimension(config%n_g_sw, nlev) :: trans_dir_dir

    ! fluxes per g point
    real(jprb), dimension(config%n_g_sw, nlev+1) :: flux_up, flux_dn_diffuse, flux_dn_direct

    ! combined optical depth, single scattering albedo and asymmetry
    ! factor
    real(jprb), dimension(config%n_g_sw) :: od_total, ssa_total, g_total

    ! two-stream coefficients
    real(jprb), dimension(config%n_g_sw) :: gamma1, gamma2, gamma3

    ! temporary working array
    real(jprb), dimension(config%n_g_sw,nlev+1) :: tmp_work_albedo, tmp_work_source
    real(jprb), dimension(config%n_g_sw,nlev) :: tmp_work_inv_denominator

    ! number of g points
    integer :: ng

    ! loop indices for level and column
    integer :: jlev, jcol

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloudless_sw:solver_cloudless_sw',0,hook_handle)

    ng = config%n_g_sw

    ! loop through columns
    do jcol = istartcol,iendcol
      ! only perform calculation if sun above the horizon
      if (single_level%cos_sza(jcol) > 0.0_jprb) then

        cos_sza = single_level%cos_sza(jcol)
        
        ! the following is the same as the clear-sky part of
        ! solver_homogeneous_sw
        if (.not. config%do_sw_delta_scaling_with_gases) then
          ! delta-eddington scaling has already been performed to the
          ! aerosol part of od, ssa and g
          do jlev = 1,nlev
            call calc_two_stream_gammas_sw(ng, cos_sza, &
                 &  ssa(:,jlev,jcol), g(:,jlev,jcol), &
                 &  gamma1, gamma2, gamma3)
            call calc_reflectance_transmittance_sw(ng, &
                 &  cos_sza, &
                 &  od(:,jlev,jcol), ssa(:,jlev,jcol), &
                 &  gamma1, gamma2, gamma3, &
                 &  reflectance(:,jlev), transmittance(:,jlev), &
                 &  ref_dir(:,jlev), trans_dir_diff(:,jlev), &
                 &  trans_dir_dir(:,jlev) )
          end do
        else
          ! apply delta-eddington scaling to the aerosol-gas mixture
          do jlev = 1,nlev
            od_total  =  od(:,jlev,jcol)
            ssa_total = ssa(:,jlev,jcol)
            g_total   =   g(:,jlev,jcol)
            call delta_eddington(od_total, ssa_total, g_total)
            call calc_two_stream_gammas_sw(ng, &
                 &  cos_sza, ssa_total, g_total, &
                 &  gamma1, gamma2, gamma3)
            call calc_reflectance_transmittance_sw(ng, &
                 &  cos_sza, od_total, ssa_total, &
                 &  gamma1, gamma2, gamma3, &
                 &  reflectance(:,jlev), transmittance(:,jlev), &
                 &  ref_dir(:,jlev), trans_dir_diff(:,jlev), &
                 &  trans_dir_dir(:,jlev) )
          end do
        end if
          
        ! use adding method to compute fluxes
        call adding_ica_sw(ng, nlev, incoming_sw(:,jcol), &
             &  albedo_diffuse(:,jcol), albedo_direct(:,jcol), &
             &  cos_sza, reflectance, transmittance, ref_dir, trans_dir_diff, &
             &  trans_dir_dir, flux_up, flux_dn_diffuse, flux_dn_direct, &
             &  albedo=tmp_work_albedo, &
             &  source=tmp_work_source, &
             &  inv_denominator=tmp_work_inv_denominator)
        
        ! sum over g-points to compute and save clear-sky broadband
        ! fluxes
        flux%sw_up(jcol,:) = sum(flux_up,1)
        if (allocated(flux%sw_dn_direct)) then
          flux%sw_dn_direct(jcol,:) = sum(flux_dn_direct,1)
          flux%sw_dn(jcol,:) = sum(flux_dn_diffuse,1) &
               &  + flux%sw_dn_direct(jcol,:)
        else
          flux%sw_dn(jcol,:) = sum(flux_dn_diffuse,1) + sum(flux_dn_direct,1)
        end if
        ! store spectral downwelling fluxes at surface
        flux%sw_dn_diffuse_surf_g(:,jcol) = flux_dn_diffuse(:,nlev+1)
        flux%sw_dn_direct_surf_g(:,jcol)  = flux_dn_direct(:,nlev+1)

        ! save the spectral fluxes if required
        if (config%do_save_spectral_flux) then
          call indexed_sum_profile(flux_up, config%i_spec_from_reordered_g_sw, &
               &                   flux%sw_up_band(:,jcol,:))
          call indexed_sum_profile(flux_dn_direct, config%i_spec_from_reordered_g_sw, &
               &                   flux%sw_dn_band(:,jcol,:))
          if (allocated(flux%sw_dn_direct_band)) then
            flux%sw_dn_direct_band(:,jcol,:) &
                 &  = flux%sw_dn_band(:,jcol,:)
          end if
          call add_indexed_sum_profile(flux_dn_diffuse, &
               &                       config%i_spec_from_reordered_g_sw, &
               &                       flux%sw_dn_band(:,jcol,:))
        end if

        if (config%do_clear) then
          ! clear-sky calculations are equal to all-sky for this
          ! solver: copy fluxes over
          flux%sw_up_clear(jcol,:) = flux%sw_up(jcol,:)
          flux%sw_dn_clear(jcol,:) = flux%sw_dn(jcol,:)
          if (allocated(flux%sw_dn_direct_clear)) then
            flux%sw_dn_direct_clear(jcol,:) = flux%sw_dn_direct(jcol,:)
          end if
          flux%sw_dn_diffuse_surf_clear_g(:,jcol) = flux%sw_dn_diffuse_surf_g(:,jcol)
          flux%sw_dn_direct_surf_clear_g(:,jcol)  = flux%sw_dn_direct_surf_g(:,jcol)

          if (config%do_save_spectral_flux) then
            flux%sw_up_clear_band(:,jcol,:) = flux%sw_up_band(:,jcol,:)
            flux%sw_dn_clear_band(:,jcol,:) = flux%sw_dn_band(:,jcol,:)
            if (allocated(flux%sw_dn_direct_clear_band)) then
              flux%sw_dn_direct_clear_band(:,jcol,:) = flux%sw_dn_direct_band(:,jcol,:)
            end if
          end if

        end if ! do_clear

      else
        ! set fluxes to zero if sun is below the horizon
        flux%sw_up(jcol,:) = 0.0_jprb
        flux%sw_dn(jcol,:) = 0.0_jprb
        if (allocated(flux%sw_dn_direct)) then
          flux%sw_dn_direct(jcol,:) = 0.0_jprb
        end if
        flux%sw_dn_diffuse_surf_g(:,jcol) = 0.0_jprb
        flux%sw_dn_direct_surf_g(:,jcol)  = 0.0_jprb

        if (config%do_clear) then
          flux%sw_up_clear(jcol,:) = 0.0_jprb
          flux%sw_dn_clear(jcol,:) = 0.0_jprb
          if (allocated(flux%sw_dn_direct_clear)) then
            flux%sw_dn_direct_clear(jcol,:) = 0.0_jprb
          end if
          flux%sw_dn_diffuse_surf_clear_g(:,jcol) = 0.0_jprb
          flux%sw_dn_direct_surf_clear_g(:,jcol)  = 0.0_jprb
        end if

        if (config%do_save_spectral_flux) then
          flux%sw_dn_band(:,jcol,:) = 0.0_jprb
          flux%sw_up_band(:,jcol,:) = 0.0_jprb
          if (allocated(flux%sw_dn_direct_band)) then
            flux%sw_dn_direct_band(:,jcol,:) = 0.0_jprb
          end if
          if (config%do_clear) then
            flux%sw_dn_clear_band(:,jcol,:) = 0.0_jprb
            flux%sw_up_clear_band(:,jcol,:) = 0.0_jprb
            if (allocated(flux%sw_dn_direct_clear_band)) then
              flux%sw_dn_direct_clear_band(:,jcol,:) = 0.0_jprb
            end if
          end if
        end if

      end if ! sun above horizon
    end do

    if (lhook) call dr_hook('radiation_cloudless_sw:solver_cloudless_sw',1,hook_handle)

  end subroutine solver_cloudless_sw

end module radiation_cloudless_sw
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

