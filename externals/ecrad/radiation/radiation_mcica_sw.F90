! # 1 "radiation/radiation_mcica_sw.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_mcica_sw.f90"
! this file has been modified for the use in icon

! radiation_mcica_sw.f90 - monte-carlo independent column approximation shortwave solver
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
!   2017-04-11  r. hogan  receive albedos at g-points
!   2017-04-22  r. hogan  store surface fluxes at all g-points
!   2017-10-23  r. hogan  renamed single-character variables


! # 1 "radiation/ecrad_config.h" 1
! ecrad_config.h - preprocessor definitions to configure compilation ecrad -*- f90 -*-
!
! (c) copyright 2023- ecmwf.
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
! this file should be included in fortran source files that require
! different optimizations or settings for different architectures and
! platforms.  feel free to maintain a site-specific version of it.

! the following settings turn on optimizations specific to the
! long-vector nec sx (the short-vector x86-64 architecture is assumed
! otherwise). 




  



  




! in the ifs, an mpi version of easy_netcdf capability is used so that
! only one mpi task reads the data files and shares with the other
! tasks. the mpi version is not used for writing files.

!#define easy_netcdf_read_mpi 1
! # 23 "radiation/radiation_mcica_sw.f90" 2

module radiation_mcica_sw

  public

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
! # 32 "radiation/radiation_mcica_sw.f90" 2

  !---------------------------------------------------------------------
  ! shortwave monte carlo independent column approximation
  ! (mcica). this implementation performs a clear-sky and a cloudy-sky
  ! calculation, and then weights the two to get the all-sky fluxes
  ! according to the total cloud cover. this method reduces noise for
  ! low cloud cover situations, and exploits the clear-sky
  ! calculations that are usually performed for diagnostic purposes
  ! simultaneously. the cloud generator has been carefully written
  ! such that the stochastic cloud field satisfies the prescribed
  ! overlap parameter accounting for this weighting.
  subroutine solver_mcica_sw(nlev,istartcol,iendcol, &
       &  config, single_level, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, &
       &  albedo_direct, albedo_diffuse, incoming_sw, &
       &  flux)

    use parkind1, only           : jprb
    use ecradhook,  only           : lhook, dr_hook, jphook

    use radiation_io,   only           : nulerr, radiation_abort
    use radiation_config, only         : config_type
    use radiation_single_level, only   : single_level_type
    use radiation_cloud, only          : cloud_type
    use radiation_flux, only           : flux_type
    use radiation_two_stream, only     : calc_ref_trans_sw
    use radiation_adding_ica_sw, only  : adding_ica_sw
    use radiation_cloud_generator, only: cloud_generator

    implicit none

    ! inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(cloud_type),         intent(in) :: cloud

    ! gas and aerosol optical depth, single-scattering albedo and
    ! asymmetry factor at each shortwave g-point
    real(jprb), intent(in), dimension(config%n_g_sw, nlev, istartcol:iendcol) :: &
         &  od, ssa, g

    ! cloud and precipitation optical depth, single-scattering albedo and
    ! asymmetry factor in each shortwave band
    real(jprb), intent(in), dimension(config%n_bands_sw,nlev,istartcol:iendcol)   :: &
         &  od_cloud, ssa_cloud, g_cloud

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

    ! diffuse reflectance and transmittance for each layer in clear
    ! and all skies
    real(jprb), dimension(config%n_g_sw, nlev) :: ref_clear, trans_clear, reflectance, transmittance

    ! fraction of direct beam scattered by a layer into the upwelling
    ! or downwelling diffuse streams, in clear and all skies
    real(jprb), dimension(config%n_g_sw, nlev) :: ref_dir_clear, trans_dir_diff_clear, ref_dir, trans_dir_diff

    ! transmittance for the direct beam in clear and all skies
    real(jprb), dimension(config%n_g_sw, nlev) :: trans_dir_dir_clear, trans_dir_dir

    ! fluxes per g point
    real(jprb), dimension(config%n_g_sw, nlev+1) :: flux_up, flux_dn_diffuse, flux_dn_direct

    ! combined gas+aerosol+cloud optical depth, single scattering
    ! albedo and asymmetry factor
    real(jprb), dimension(config%n_g_sw) :: od_total, ssa_total, g_total

    ! combined scattering optical depth
    real(jprb) :: scat_od

    ! optical depth scaling from the cloud generator, zero indicating
    ! clear skies
    real(jprb), dimension(config%n_g_sw,nlev) :: od_scaling

    ! modified optical depth after mcica scaling to represent cloud
    ! inhomogeneity
    real(jprb), dimension(config%n_g_sw) :: od_cloud_new

    ! temporary working array
    real(jprb), dimension(config%n_g_sw,nlev+1) :: tmp_work_albedo, &
      &                                            tmp_work_source
    real(jprb), dimension(config%n_g_sw,nlev) :: tmp_work_inv_denominator

    ! total cloud cover output from the cloud generator
    real(jprb) :: total_cloud_cover

    ! temporary storage for more efficient summation



    real(jprb) :: sum_up, sum_dn_diff, sum_dn_dir


    ! number of g points
    integer :: ng

    ! loop indices for level, column and g point
    integer :: jlev, jcol, jg

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_mcica_sw:solver_mcica_sw',0,hook_handle)

    if (.not. config%do_clear) then
      write(nulerr,'(a)') '*** error: shortwave mcica requires clear-sky calculation to be performed'
      call radiation_abort()      
    end if

    ng = config%n_g_sw

    ! loop through columns
    do jcol = istartcol,iendcol
      ! only perform calculation if sun above the horizon
      if (single_level%cos_sza(jcol) > 0.0_jprb) then
        cos_sza = single_level%cos_sza(jcol)

        ! clear-sky calculation - first compute clear-sky reflectance,
        ! transmittance etc at each model level
        if (.not. config%do_sw_delta_scaling_with_gases) then
          ! delta-eddington scaling has already been performed to the
          ! aerosol part of od, ssa and g
          call calc_ref_trans_sw(ng*nlev, &
               &  cos_sza, od(:,:,jcol), ssa(:,:,jcol), g(:,:,jcol), &
               &  ref_clear, trans_clear, &
               &  ref_dir_clear, trans_dir_diff_clear, &
               &  trans_dir_dir_clear)
        else
          ! apply delta-eddington scaling to the aerosol-gas mixture
          do jlev = 1,nlev
            od_total  =  od(:,jlev,jcol)
            ssa_total = ssa(:,jlev,jcol)
            g_total   =   g(:,jlev,jcol)
            call delta_eddington(od_total, ssa_total, g_total)
            call calc_ref_trans_sw(ng, &
                 &  cos_sza, od_total, ssa_total, g_total, &
                 &  ref_clear(:,jlev), trans_clear(:,jlev), &
                 &  ref_dir_clear(:,jlev), trans_dir_diff_clear(:,jlev), &
                 &  trans_dir_dir_clear(:,jlev) )
          end do
        end if

        ! use adding method to compute fluxes
        call adding_ica_sw(ng, nlev, incoming_sw(:,jcol), &
             &  albedo_diffuse(:,jcol), albedo_direct(:,jcol), cos_sza, &
             &  ref_clear, trans_clear, ref_dir_clear, trans_dir_diff_clear, &
             &  trans_dir_dir_clear, flux_up, flux_dn_diffuse, flux_dn_direct, &
             &  albedo=tmp_work_albedo, &
             &  source=tmp_work_source, &
             &  inv_denominator=tmp_work_inv_denominator)
        
        ! sum over g-points to compute and save clear-sky broadband
        ! fluxes. note that the built-in "sum" function is very slow,
        ! and before being replaced by the alternatives below
        ! accounted for around 40% of the total cost of this routine.
! # 215 "radiation/radiation_mcica_sw.f90"
        ! optimized summation for the x86-64 architecture
        do jlev = 1,nlev+1
          sum_up      = 0.0_jprb
          sum_dn_diff = 0.0_jprb
          sum_dn_dir  = 0.0_jprb
          !$omp simd reduction(+:sum_up, sum_dn_diff, sum_dn_dir)
          do jg = 1,ng
            sum_up      = sum_up      + flux_up(jg,jlev)
            sum_dn_diff = sum_dn_diff + flux_dn_diffuse(jg,jlev)
            sum_dn_dir  = sum_dn_dir  + flux_dn_direct(jg,jlev)
          end do
          flux%sw_up_clear(jcol,jlev) = sum_up
          flux%sw_dn_clear(jcol,jlev) = sum_dn_diff + sum_dn_dir
          if (allocated(flux%sw_dn_direct_clear)) then
            flux%sw_dn_direct_clear(jcol,jlev) = sum_dn_dir
          end if
        end do

        
        ! store spectral downwelling fluxes at surface
        do jg = 1,ng
          flux%sw_dn_diffuse_surf_clear_g(jg,jcol) = flux_dn_diffuse(jg,nlev+1)
          flux%sw_dn_direct_surf_clear_g(jg,jcol)  = flux_dn_direct(jg,nlev+1)
        end do

        ! do cloudy-sky calculation
        call cloud_generator(ng, nlev, config%i_overlap_scheme, &
             &  single_level%iseed(jcol), &
             &  config%cloud_fraction_threshold, &
             &  cloud%fraction(jcol,:), cloud%overlap_param(jcol,:), &
             &  config%cloud_inhom_decorr_scaling, cloud%fractional_std(jcol,:), &
             &  config%pdf_sampler, od_scaling, total_cloud_cover, &
             &  use_beta_overlap=config%use_beta_overlap, &
             &  use_vectorizable_generator=config%use_vectorizable_generator)

        ! store total cloud cover
        flux%cloud_cover_sw(jcol) = total_cloud_cover
        
        if (total_cloud_cover >= config%cloud_fraction_threshold) then
          ! total-sky calculation
          do jlev = 1,nlev
            ! compute combined gas+aerosol+cloud optical properties
            if (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold) then
              do jg = 1,ng
                od_cloud_new(jg) = od_scaling(jg,jlev) &
                   &  * od_cloud(config%i_band_from_reordered_g_sw(jg),jlev,jcol)
                od_total(jg)  = od(jg,jlev,jcol) + od_cloud_new(jg)
                ssa_total(jg) = 0.0_jprb
                g_total(jg)   = 0.0_jprb

                ! in single precision we need to protect against the
                ! case that od_total > 0.0 and ssa_total > 0.0 but
                ! od_total*ssa_total == 0 due to underflow
                if (od_total(jg) > 0.0_jprb) then
                  scat_od = ssa(jg,jlev,jcol)*od(jg,jlev,jcol) &
                       &     + ssa_cloud(config%i_band_from_reordered_g_sw(jg),jlev,jcol) &
                       &     *  od_cloud_new(jg)
                  ssa_total(jg) = scat_od / od_total(jg)
                  if (scat_od > 0.0_jprb) then
                    g_total(jg) = (g(jg,jlev,jcol)*ssa(jg,jlev,jcol)*od(jg,jlev,jcol) &
                         &     +   g_cloud(config%i_band_from_reordered_g_sw(jg),jlev,jcol) &
                         &     * ssa_cloud(config%i_band_from_reordered_g_sw(jg),jlev,jcol) &
                         &     *  od_cloud_new(jg)) &
                         &     / scat_od
                  end if
                end if
              end do

              ! apply delta-eddington scaling to the cloud-aerosol-gas
              ! mixture
              if (config%do_sw_delta_scaling_with_gases) then
                call delta_eddington(od_total, ssa_total, g_total)
              end if

              ! compute cloudy-sky reflectance, transmittance etc at
              ! each model level
              call calc_ref_trans_sw(ng, &
                   &  cos_sza, od_total, ssa_total, g_total, &
                   &  reflectance(:,jlev), transmittance(:,jlev), &
                   &  ref_dir(:,jlev), trans_dir_diff(:,jlev), &
                   &  trans_dir_dir(:,jlev))
              
            else
              ! clear-sky layer: copy over clear-sky values
              do jg = 1,ng
                reflectance(jg,jlev) = ref_clear(jg,jlev)
                transmittance(jg,jlev) = trans_clear(jg,jlev)
                ref_dir(jg,jlev) = ref_dir_clear(jg,jlev)
                trans_dir_diff(jg,jlev) = trans_dir_diff_clear(jg,jlev)
                trans_dir_dir(jg,jlev) = trans_dir_dir_clear(jg,jlev)
              end do
            end if
          end do
            
          ! use adding method to compute fluxes for an overcast sky
          call adding_ica_sw(ng, nlev, incoming_sw(:,jcol), &
               &  albedo_diffuse(:,jcol), albedo_direct(:,jcol), cos_sza, &
               &  reflectance, transmittance, ref_dir, trans_dir_diff, &
               &  trans_dir_dir, flux_up, flux_dn_diffuse, flux_dn_direct, &
               &  albedo=tmp_work_albedo, &
               &  source=tmp_work_source, &
               &  inv_denominator=tmp_work_inv_denominator)
          
          ! store overcast broadband fluxes
! # 334 "radiation/radiation_mcica_sw.f90"
          do jlev = 1,nlev+1
            sum_up      = 0.0_jprb
            sum_dn_diff = 0.0_jprb
            sum_dn_dir  = 0.0_jprb
            !$omp simd reduction(+:sum_up, sum_dn_diff, sum_dn_dir)
            do jg = 1,ng
              sum_up      = sum_up      + flux_up(jg,jlev)
              sum_dn_diff = sum_dn_diff + flux_dn_diffuse(jg,jlev)
              sum_dn_dir  = sum_dn_dir  + flux_dn_direct(jg,jlev)
            end do
            flux%sw_up(jcol,jlev) = sum_up
            flux%sw_dn(jcol,jlev) = sum_dn_diff + sum_dn_dir
            if (allocated(flux%sw_dn_direct)) then
              flux%sw_dn_direct(jcol,jlev) = sum_dn_dir
            end if
          end do

          
          ! cloudy flux profiles currently assume completely overcast
          ! skies; perform weighted average with clear-sky profile
          do jlev = 1, nlev+1
            flux%sw_up(jcol,jlev) =  total_cloud_cover *flux%sw_up(jcol,jlev) &
                 &     + (1.0_jprb - total_cloud_cover)*flux%sw_up_clear(jcol,jlev)
            flux%sw_dn(jcol,jlev) =  total_cloud_cover *flux%sw_dn(jcol,jlev) &
                 &     + (1.0_jprb - total_cloud_cover)*flux%sw_dn_clear(jcol,jlev)
            if (allocated(flux%sw_dn_direct)) then
              flux%sw_dn_direct(jcol,jlev) = total_cloud_cover *flux%sw_dn_direct(jcol,jlev) &
                   &  + (1.0_jprb - total_cloud_cover)*flux%sw_dn_direct_clear(jcol,jlev)
            end if
          end do
          ! likewise for surface spectral fluxes
          do jg = 1,ng
            flux%sw_dn_diffuse_surf_g(jg,jcol) = flux_dn_diffuse(jg,nlev+1)
            flux%sw_dn_direct_surf_g(jg,jcol)  = flux_dn_direct(jg,nlev+1)
            flux%sw_dn_diffuse_surf_g(jg,jcol) = total_cloud_cover *flux%sw_dn_diffuse_surf_g(jg,jcol) &
                 &                 + (1.0_jprb - total_cloud_cover)*flux%sw_dn_diffuse_surf_clear_g(jg,jcol)
            flux%sw_dn_direct_surf_g(jg,jcol)  = total_cloud_cover *flux%sw_dn_direct_surf_g(jg,jcol) &
                 &                 + (1.0_jprb - total_cloud_cover)*flux%sw_dn_direct_surf_clear_g(jg,jcol)
          end do

        else
          ! no cloud in profile and clear-sky fluxes already
          ! calculated: copy them over
          do jlev = 1, nlev+1
            flux%sw_up(jcol,jlev) = flux%sw_up_clear(jcol,jlev)
            flux%sw_dn(jcol,jlev) = flux%sw_dn_clear(jcol,jlev)
            if (allocated(flux%sw_dn_direct)) then
              flux%sw_dn_direct(jcol,jlev) = flux%sw_dn_direct_clear(jcol,jlev)
            end if
          end do
          do jg = 1,ng
            flux%sw_dn_diffuse_surf_g(jg,jcol) = flux%sw_dn_diffuse_surf_clear_g(jg,jcol)
            flux%sw_dn_direct_surf_g(jg,jcol)  = flux%sw_dn_direct_surf_clear_g(jg,jcol)
          end do

        end if ! cloud is present in profile

      else
        ! set fluxes to zero if sun is below the horizon
        do jlev = 1, nlev+1
          flux%sw_up(jcol,jlev) = 0.0_jprb
          flux%sw_dn(jcol,jlev) = 0.0_jprb
          if (allocated(flux%sw_dn_direct)) then
            flux%sw_dn_direct(jcol,jlev) = 0.0_jprb
          end if
          flux%sw_up_clear(jcol,jlev) = 0.0_jprb
          flux%sw_dn_clear(jcol,jlev) = 0.0_jprb
          if (allocated(flux%sw_dn_direct_clear)) then
            flux%sw_dn_direct_clear(jcol,jlev) = 0.0_jprb
          end if
        end do
        do jg = 1,ng
          flux%sw_dn_diffuse_surf_g(jg,jcol) = 0.0_jprb
          flux%sw_dn_direct_surf_g(jg,jcol)  = 0.0_jprb
          flux%sw_dn_diffuse_surf_clear_g(jg,jcol) = 0.0_jprb
          flux%sw_dn_direct_surf_clear_g(jg,jcol)  = 0.0_jprb
        end do
      end if ! sun above horizon

    end do ! loop over columns

    if (lhook) call dr_hook('radiation_mcica_sw:solver_mcica_sw',1,hook_handle)
    
  end subroutine solver_mcica_sw

end module radiation_mcica_sw
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

