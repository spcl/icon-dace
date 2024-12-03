! # 1 "radiation/radiation_mcica_acc_lw.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_mcica_acc_lw.f90"
! this file has been modified for the use in icon

! radiation_mcica_acc_lw.f90 - monte-carlo independent column approximation longtwave solver
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
!   2017-04-11  r. hogan  receive emission/albedo rather than planck/emissivity
!   2017-04-22  r. hogan  store surface fluxes at all g-points
!   2017-07-12  r. hogan  call fast adding method if only clouds scatter
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
! # 24 "radiation/radiation_mcica_acc_lw.f90" 2

module radiation_mcica_acc_lw

  public

contains

  !---------------------------------------------------------------------
  ! longwave monte carlo independent column approximation
  ! (mcica). this implementation performs a clear-sky and a cloudy-sky
  ! calculation, and then weights the two to get the all-sky fluxes
  ! according to the total cloud cover. this method reduces noise for
  ! low cloud cover situations, and exploits the clear-sky
  ! calculations that are usually performed for diagnostic purposes
  ! simultaneously. the cloud generator has been carefully written
  ! such that the stochastic cloud field satisfies the prescribed
  ! overlap parameter accounting for this weighting.
  subroutine solver_mcica_acc_lw(nlev,istartcol,iendcol, &
       &  config, single_level, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
       &  emission, albedo, &
       &  flux)

    use parkind1, only           : jprb
    use ecradhook,  only           : lhook, dr_hook

    use radiation_io,   only           : nulerr, radiation_abort
    use radiation_config, only         : config_type
    use radiation_single_level, only   : single_level_type
    use radiation_cloud, only          : cloud_type
    use radiation_flux, only           : flux_type
    use radiation_two_stream, only     : calc_ref_trans_lw, &
         &                               calc_no_scattering_transmittance_lw
    use radiation_adding_ica_lw, only  : adding_ica_lw, fast_adding_ica_lw, &
         &                               calc_fluxes_no_scattering_lw
    use radiation_lw_derivatives, only : calc_lw_derivatives_ica, modify_lw_derivatives_ica
    use radiation_cloud_generator_acc, only: cloud_generator_acc
    use radiation_cloud_cover, only    : beta2alpha, maxcloudfrac

    implicit none

    ! inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
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
    real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) :: emission, albedo

    ! output
    type(flux_type), intent(inout):: flux

    ! local variables

    ! diffuse reflectance and transmittance for each layer in clear
    ! and all skies
    real(jprb), dimension(config%n_g_lw, nlev) :: ref_clear, trans_clear, reflectance, transmittance

    ! emission by a layer into the upwelling or downwelling diffuse
    ! streams, in clear and all skies
    real(jprb), dimension(config%n_g_lw, nlev) :: source_up_clear, source_dn_clear, source_up, source_dn

    ! fluxes per g point
    real(jprb), dimension(config%n_g_lw, nlev+1, istartcol:iendcol) :: flux_up, flux_dn
    real(jprb), dimension(config%n_g_lw, nlev+1, istartcol:iendcol) :: flux_up_clear, flux_dn_clear

    ! combined gas+aerosol+cloud optical depth, single scattering
    ! albedo and asymmetry factor
    real(jprb), dimension(config%n_g_lw) :: od_total, ssa_total, g_total

    ! combined scattering optical depth
    real(jprb) :: scat_od, scat_od_total(config%n_g_lw)

    ! optical depth scaling from the cloud generator, zero indicating
    ! clear skies
    real(jprb), dimension(config%n_g_lw,nlev) :: od_scaling

    ! modified optical depth after mcica scaling to represent cloud
    ! inhomogeneity
    real(jprb), dimension(config%n_g_lw) :: od_cloud_new

    ! total cloud cover output from the cloud generator
    real(jprb) :: total_cloud_cover

    ! identify clear-sky layers
    logical :: is_clear_sky_layer(nlev, istartcol:iendcol)

    ! workaround that allows inling of cloud generator
    real(jprb), dimension(config%pdf_sampler%ncdf, config%pdf_sampler%nfsd)  :: sample_val

    ! temporary arrays to increase performance
    real(jprb), dimension(nlev, istartcol:iendcol) :: frac, frac_std
    real(jprb), dimension(nlev-1, istartcol:iendcol) :: overlap_param

    ! temporary arrays
    real(jprb), dimension(nlev, istartcol:iendcol) :: cum_cloud_cover
    real(jprb), dimension(nlev-1, istartcol:iendcol) :: pair_cloud_cover

    ! "alpha" overlap parameter
    real(jprb) :: overlap_alpha

    ! cumulative product needed in computation of total_cloud_cover
    real(jprb) :: cum_product(istartcol:iendcol)

    ! first and last cloudy layers
    integer :: ibegin(istartcol:iendcol), iend(istartcol:iendcol)

    ! temporary working array
    real(jprb), dimension(config%n_g_lw,nlev+1) :: tmp_work_albedo, &
      &                                            tmp_work_source
    real(jprb), dimension(config%n_g_lw,nlev) :: tmp_work_inv_denominator
    ! temporary storage for more efficient summation
    real(jprb) :: sum_up, sum_dn, sum_up_clr, sum_dn_clr

    ! index of the highest cloudy layer
    integer :: i_cloud_top

    ! number of g points
    integer :: ng

    ! loop indices for level, column and g point
    integer :: jlev, jcol, jg

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_mcica_acc_lw:solver_mcica_acc_lw',0,hook_handle)

    if (.not. config%do_clear) then
      write(nulerr,'(a)') '*** error: longwave mcica acc requires clear-sky calculation to be performed'
      call radiation_abort()      
    end if

    ng = config%n_g_lw

    !$acc data create(flux_up, flux_dn, flux_up_clear, flux_dn_clear, &
    !$acc             is_clear_sky_layer, &
    !$acc             sample_val, frac, frac_std, overlap_param, cum_cloud_cover, &
    !$acc             pair_cloud_cover, cum_product, ibegin, iend) &
    !$acc     present(config, single_level, cloud, od, ssa, g, od_cloud, ssa_cloud, &
    !$acc             g_cloud, planck_hl, emission, albedo, flux)

    !$acc parallel default(none) async(1)
    !$acc loop gang vector collapse(2)
    do jlev = 1,config%pdf_sampler%nfsd
      do jcol = 1,config%pdf_sampler%ncdf
        sample_val(jcol,jlev) = config%pdf_sampler%val(jcol,jlev)
      end do
    end do
    !$acc end parallel

    !$acc parallel default(none) async(1)
    !$acc loop gang vector collapse(2)
    do jcol = istartcol,iendcol
      do jlev = 1, nlev
        frac(jlev, jcol) = cloud%fraction(jcol,jlev)
        frac_std(jlev, jcol) = cloud%fractional_std(jcol,jlev)
        is_clear_sky_layer(jlev,jcol) = .true.
      end do
    end do
    !$acc end parallel

    !$acc parallel default(none) async(1)
    !$acc loop gang vector collapse(2)
    do jcol = istartcol,iendcol
      do jlev = 1, nlev-1
        overlap_param(jlev, jcol) = cloud%overlap_param(jcol,jlev)
      end do
    end do
    !$acc end parallel


    !$acc parallel default(none) async(1)
    !$acc loop gang vector collapse(2) private(overlap_alpha)
    do jcol = istartcol,iendcol       
      !---------------------------------------------------------------------
      ! manual inline from cum_cloud_cover_exp_ran >>>>>>>>>>>>>>>>>>>>>>>>
      ! loop to compute total cloud cover and the cumulative cloud cover
      ! down to the base of each layer
      do jlev = 1,nlev-1      
        ! convert to "alpha" overlap parameter if necessary
        if (config%use_beta_overlap) then
          overlap_alpha = beta2alpha(overlap_param(jlev,jcol), &
                &                     frac(jlev,jcol), frac(jlev+1,jcol))
        else
          overlap_alpha = overlap_param(jlev,jcol)
        end if

        ! compute the combined cloud cover of layers jlev and jlev+1
        pair_cloud_cover(jlev, jcol) = overlap_alpha*max(frac(jlev,jcol),frac(jlev+1,jcol)) &
              &  + (1.0_jprb - overlap_alpha) &
              &  * (frac(jlev,jcol)+frac(jlev+1,jcol)-frac(jlev,jcol)*frac(jlev+1,jcol))
      end do
    end do
    !$acc end parallel
    
    !$acc parallel default(none) async(1)
    !$acc loop gang vector
    do jcol = istartcol,iendcol
      cum_cloud_cover(1, jcol) = frac(1,jcol)        
      cum_product(jcol) = 1.0_jprb - frac(1,jcol)
      !$acc loop seq 
      do jlev = 1,nlev-1
        if (frac(jlev,jcol) >= maxcloudfrac) then
          ! cloud cover has reached one
          cum_product(jcol) = 0.0_jprb
        else
          cum_product(jcol) = cum_product(jcol) * (1.0_jprb - pair_cloud_cover(jlev, jcol)) &
                &  / (1.0_jprb - frac(jlev,jcol))
        end if
        cum_cloud_cover(jlev+1, jcol) = 1.0_jprb - cum_product(jcol) 
      end do
      flux%cloud_cover_lw(jcol) = cum_cloud_cover(nlev,jcol);
      if (flux%cloud_cover_lw(jcol) < config%cloud_fraction_threshold) then
        ! treat column as clear sky: calling function therefore will not
        ! use od_scaling so we don't need to calculate it
        flux%cloud_cover_lw(jcol) = 0.0_jprb
      end if
    end do
    !$acc end parallel

    !$acc parallel default(none) async(1)
    !$acc loop gang vector
    do jcol = istartcol,iendcol
      if (flux%cloud_cover_lw(jcol) >= config%cloud_fraction_threshold) then
        ! cloud is present: need to calculate od_scaling
        ! find range of cloudy layers
        ibegin(jcol) = nlev
        !$acc loop seq
        do jlev = 1, nlev
          if( frac(jlev,jcol) > 0.0_jprb ) then
            ibegin(jcol) = min(jlev, ibegin(jcol))
          end if
        end do

        iend(jcol) = ibegin(jcol)
        !$acc loop seq
        do jlev = ibegin(jcol)+1,nlev
          if (frac(jlev,jcol) > 0.0_jprb) then
            iend(jcol) = max(jlev, iend(jcol))
          end if
        end do
      end if
    end do
    !$acc end parallel


    ! loop through columns
    !$acc parallel default(none) async(1)  &
    !$acc   num_gangs(iendcol-istartcol+1) num_workers((config%n_g_lw-1)/32+1) vector_length(32)
    !$acc loop gang private(g_total, i_cloud_top, od_cloud_new, od_scaling, od_total, ref_clear, &
    !$acc   reflectance, tmp_work_inv_denominator, tmp_work_albedo, tmp_work_source, &
    !$acc   scat_od_total, source_dn, source_dn_clear, source_up, source_up_clear, ssa_total, &
    !$acc   trans_clear, transmittance)
    do jcol = istartcol,iendcol

      ! clear-sky calculation

      if (config%do_lw_aerosol_scattering) then
        ! scattering case: first compute clear-sky reflectance,
        ! transmittance etc at each model level
        call calc_ref_trans_lw(ng*nlev, &
             &  od(:,:,jcol), ssa(:,:,jcol), g(:,:,jcol), &
             &  planck_hl(:,1:jlev,jcol), planck_hl(:,2:jlev+1,jcol), &
             &  ref_clear, trans_clear, &
             &  source_up_clear, source_dn_clear)
        ! then use adding method to compute fluxes
        call adding_ica_lw(ng, nlev, &
             &  ref_clear, trans_clear, source_up_clear, source_dn_clear, &
             &  emission(:,jcol), albedo(:,jcol), &
             &  flux_up_clear(:,:,jcol), flux_dn_clear(:,:,jcol))
      else

        ! non-scattering case: use simpler functions for
        ! transmission and emission
        !$acc loop seq
        do jlev = 1,nlev
          call calc_no_scattering_transmittance_lw(ng, od(:,jlev,jcol), &
               &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1, jcol), &
               &  trans_clear(:,jlev), source_up_clear(:,jlev), &
               &  source_dn_clear(:,jlev))
        end do
        ! ensure that clear-sky reflectance is zero since it may be
        ! used in cloudy-sky case
        !$acc loop seq
        do jlev = 1,nlev
          !$acc loop worker vector
          do jg = 1,ng
            ref_clear(jg,jlev) = 0.0_jprb
          end do
        end do
        ! simpler down-then-up method to compute fluxes
        call calc_fluxes_no_scattering_lw(ng, nlev, &
             &  trans_clear, source_up_clear, source_dn_clear, &
             &  emission(:,jcol), albedo(:,jcol), &
             &  flux_up_clear(:,:,jcol), flux_dn_clear(:,:,jcol))

      end if


      ! store surface spectral downwelling fluxes
      !$acc loop worker vector
      do jg = 1,ng
        flux%lw_dn_surf_clear_g(jg,jcol) = flux_dn_clear(jg,nlev+1,jcol)
      end do

      ! do cloudy-sky calculation; add a prime number to the seed in
      ! the longwave
      call cloud_generator_acc(ng, nlev, &
           &  single_level%iseed(jcol) + 997, &
           &  config%cloud_fraction_threshold, &
           &  frac(:,jcol), overlap_param(:,jcol), &
           &  config%cloud_inhom_decorr_scaling, frac_std(:,jcol), &
           &  config%pdf_sampler%ncdf, config%pdf_sampler%nfsd, &
           &  config%pdf_sampler%fsd1, config%pdf_sampler%inv_fsd_interval, &
           &  sample_val, &
           &  od_scaling, flux%cloud_cover_lw(jcol)+0.0_jprb, & ! workaround for nvhpc-24.1
           &  ibegin(jcol), iend(jcol), &
           &  cum_cloud_cover=cum_cloud_cover(:,jcol), &
           &  pair_cloud_cover=pair_cloud_cover(:,jcol))
      
      if (flux%cloud_cover_lw(jcol) >= config%cloud_fraction_threshold) then
        ! total-sky calculation

        i_cloud_top = nlev+1
        !$acc loop seq
        do jlev = 1,nlev
          ! compute combined gas+aerosol+cloud optical properties
          if (frac(jlev,jcol) >= config%cloud_fraction_threshold) then
            is_clear_sky_layer(jlev,jcol) = .false.
            ! get index to the first cloudy layer from the top
            if (i_cloud_top > jlev) then
              i_cloud_top = jlev
            end if

            !$acc loop worker vector
            do jg = 1,ng
              od_cloud_new(jg) = od_scaling(jg,jlev) &
                 &  * od_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol)
              od_total(jg)  = od(jg,jlev,jcol) + od_cloud_new(jg)
              ssa_total(jg) = 0.0_jprb
              g_total(jg)   = 0.0_jprb
            end do

            if (config%do_lw_cloud_scattering) then
              ! scattering case: calculate reflectance and
              ! transmittance at each model level


              if (config%do_lw_aerosol_scattering) then
                ! in single precision we need to protect against the
                ! case that od_total > 0.0 and ssa_total > 0.0 but
                ! od_total*ssa_total == 0 due to underflow
                do jg = 1,ng
                  if (od_total(jg) > 0.0_jprb) then
                    scat_od_total(jg) = ssa(jg,jlev,jcol)*od(jg,jlev,jcol) &
                     &     + ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                     &     *  od_cloud_new(jg)
                    ssa_total(jg) = scat_od_total(jg) / od_total(jg)

                    if (scat_od_total(jg) > 0.0_jprb) then
                      g_total(jg) = (g(jg,jlev,jcol)*ssa(jg,jlev,jcol)*od(jg,jlev,jcol) &
                         &     +   g_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                         &     * ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                         &     *  od_cloud_new(jg)) &
                         &     / scat_od_total(jg)
                    end if
                  end if
                end do

              else


                !$acc loop worker vector private(scat_od)
                do jg = 1,ng
                  if (od_total(jg) > 0.0_jprb) then
                    scat_od = ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                         &     * od_cloud_new(jg)
                    ssa_total(jg) = scat_od / od_total(jg)
                    if (scat_od > 0.0_jprb) then
                      g_total(jg) = g_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                           &     * ssa_cloud(config%i_band_from_reordered_g_lw(jg),jlev,jcol) &
                           &     *  od_cloud_new(jg) / scat_od
                    end if
                  end if
                end do


              end if

            
              ! compute cloudy-sky reflectance, transmittance etc at
              ! each model level
              call calc_ref_trans_lw(ng, &
                   &  od_total, ssa_total, g_total, &
                   &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1,jcol), &
                   &  reflectance(:,jlev), transmittance(:,jlev), &
                   &  source_up(:,jlev), source_dn(:,jlev))
            else
              ! no-scattering case: use simpler functions for
              ! transmission and emission
              call calc_no_scattering_transmittance_lw(ng, od_total, &
                   &  planck_hl(:,jlev,jcol), planck_hl(:,jlev+1, jcol), &
                   &  transmittance(:,jlev), source_up(:,jlev), &
                   &  source_dn(:,jlev))
            end if

          else
            ! clear-sky layer: copy over clear-sky values
            !$acc loop worker vector
            do jg = 1,ng
              reflectance(jg,jlev) = ref_clear(jg,jlev)
              transmittance(jg,jlev) = trans_clear(jg,jlev)
              source_up(jg,jlev) = source_up_clear(jg,jlev)
              source_dn(jg,jlev) = source_dn_clear(jg,jlev)
            end do
          end if
        end do
        

        if (config%do_lw_aerosol_scattering) then
          ! use adding method to compute fluxes for an overcast sky,
          ! allowing for scattering in all layers
          call adding_ica_lw(ng, nlev, reflectance, transmittance, source_up, source_dn, &
               &  emission(:,jcol), albedo(:,jcol), &
               &  flux_up(:,:,jcol), flux_dn(:,:,jcol))
        else if (config%do_lw_cloud_scattering) then



          ! use adding method to compute fluxes but optimize for the
          ! presence of clear-sky layers
          call fast_adding_ica_lw(ng, nlev, reflectance, transmittance, source_up, source_dn, &
               &  emission(:,jcol), albedo(:,jcol), &
               &  is_clear_sky_layer(:,jcol), i_cloud_top, flux_dn_clear(:,:,jcol), &
               &  flux_up(:,:,jcol), flux_dn(:,:,jcol), &
               &  albedo=tmp_work_albedo, &
               &  source=tmp_work_source, &
               &  inv_denominator=tmp_work_inv_denominator)
        else
          ! simpler down-then-up method to compute fluxes
          call calc_fluxes_no_scattering_lw(ng, nlev, &
               &  transmittance, source_up, source_dn, emission(:,jcol), albedo(:,jcol), &
               &  flux_up(:,:,jcol), flux_dn(:,:,jcol))
        end if

        ! cloudy flux profiles currently assume completely overcast
        ! skies; perform weighted average with clear-sky profile
        ! store surface spectral downwelling fluxes
        !$acc loop worker vector
        do jg = 1,ng
          flux%lw_dn_surf_g(jg,jcol) = flux%cloud_cover_lw(jcol)*flux_dn(jg,nlev+1,jcol) &
              &  + (1.0_jprb - flux%cloud_cover_lw(jcol))*flux%lw_dn_surf_clear_g(jg,jcol)
        end do

        ! compute the longwave derivatives needed by hogan and bozzo
        ! (2015) approximate radiation update scheme

        if (config%do_lw_derivatives) then
          call calc_lw_derivatives_ica(ng, nlev, jcol, transmittance, flux_up(:,nlev+1,jcol), &
               &                       flux%lw_derivatives)
          if (flux%cloud_cover_lw(jcol) < 1.0_jprb - config%cloud_fraction_threshold) then
            ! modify the existing derivative with the contribution from the clear sky
            call modify_lw_derivatives_ica(ng, nlev, jcol, trans_clear, flux_up_clear(:,nlev+1,jcol), &
                 &                         1.0_jprb-flux%cloud_cover_lw(jcol), flux%lw_derivatives)
          end if
        end if


      else
        ! no cloud in profile and clear-sky fluxes already
        ! calculated: copy them over
        !$acc loop worker vector
        do jg = 1,ng
          flux%lw_dn_surf_g(jg,jcol) = flux%lw_dn_surf_clear_g(jg,jcol)
        end do

        if (config%do_lw_derivatives) then
          call calc_lw_derivatives_ica(ng, nlev, jcol, trans_clear, flux_up_clear(:,nlev+1,jcol), &
               &                       flux%lw_derivatives)
 
        end if

      end if ! cloud is present in profile
    end do
    !$acc end parallel

    ! loop through columns
    !$acc parallel default(none) async(1) &
    !$acc   num_gangs((iendcol-istartcol+1)*(nlev+1)) num_workers(1) vector_length(32*((config%n_g_lw-1)/32+1))
    !$acc loop gang collapse(2) private(total_cloud_cover, sum_up, sum_dn, sum_up_clr, sum_dn_clr)
    do jcol = istartcol,iendcol
      do jlev = 1,nlev+1

        sum_up_clr = 0._jprb
        sum_dn_clr = 0._jprb
        !$acc loop vector reduction(+:sum_up_clr, sum_dn_clr)
        do jg = 1,ng
          sum_up_clr = sum_up_clr + flux_up_clear(jg,jlev,jcol)
          sum_dn_clr = sum_dn_clr + flux_dn_clear(jg,jlev,jcol)
        end do
        flux%lw_up_clear(jcol,jlev) = sum_up_clr
        flux%lw_dn_clear(jcol,jlev) = sum_dn_clr

        total_cloud_cover = flux%cloud_cover_lw(jcol)
        if (total_cloud_cover >= config%cloud_fraction_threshold) then

          ! store overcast broadband fluxes
          sum_up = 0._jprb
          sum_dn = 0._jprb
          !$acc loop vector reduction(+:sum_up, sum_dn)
          do jg = 1,ng
            sum_up = sum_up + flux_up(jg,jlev,jcol)
            sum_dn = sum_dn + flux_dn(jg,jlev,jcol)
          end do
          flux%lw_up(jcol,jlev) = total_cloud_cover*sum_up + (1.0_jprb - total_cloud_cover)*sum_up_clr
          flux%lw_dn(jcol,jlev) = total_cloud_cover*sum_dn + (1.0_jprb - total_cloud_cover)*sum_dn_clr

        else

          flux%lw_up(jcol,jlev) = sum_up_clr
          flux%lw_dn(jcol,jlev) = sum_dn_clr

        end if
      end do
    end do
    !$acc end parallel

    !$acc wait
    !$acc end data

    if (lhook) call dr_hook('radiation_mcica_acc_lw:solver_mcica_acc_lw',1,hook_handle)
    
  end subroutine solver_mcica_acc_lw

end module radiation_mcica_acc_lw
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

