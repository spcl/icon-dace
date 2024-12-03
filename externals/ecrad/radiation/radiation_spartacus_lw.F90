! # 1 "radiation/radiation_spartacus_lw.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_spartacus_lw.f90"
! radiation_spartacus_lw.f90 - spartacus longwave solver
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
!   2017-04-11  r. hogan  receive emission/albedo rather than planck/emissivity
!   2017-04-22  r. hogan  store surface fluxes at all g-points
!   2018-09-03  r. hogan  security via min_cloud_effective_size
!   2018-10-08  r. hogan  call calc_region_properties
!   2019-01-12  r. hogan  use inv_inhom_effective_size if allocated

module radiation_spartacus_lw

  public

contains

  ! small routine for scaling cloud optical depth in the cloudy
  ! regions

! # 1 "radiation/radiation_optical_depth_scaling.h" 1
! radiation_optical_depth_scaling.h - cloud optical-depth scaling for tripleclouds 
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
!
! author:  robin hogan
! email:   r.j.hogan@ecmwf.int
!
! modifications
!   2017-07-14  r. hogan  incorporate gamma distribution option
!
! this file is intended to be included inside a module to ensure that
! this simple routine may be inlined

!---------------------------------------------------------------------
! compute the optical depth scalings for the optically "thick" and
! "thin" regions of a tripleclouds representation of a sub-grid pdf of
! cloud optical depth. following shonk and hogan (2008), the 16th
! percentile is used for the thin region, and the formulas estimate
! this for both lognormal and gamma distributions.
pure subroutine optical_depth_scaling(nreg, frac_std, do_gamma, od_scaling)

  use parkind1, only : jprb

  ! number of regions
  integer, intent(in)     :: nreg

  ! fractional standard deviation of in-cloud water content
  real(jprb), intent(in)  :: frac_std

  ! do we do a lognormal or gamma distribution?
  logical, intent(in) :: do_gamma

  ! optical depth scaling for the cloudy regions
  real(jprb), intent(out) :: od_scaling(2:nreg)

  if (nreg == 2) then
    ! only one clear-sky and one cloudy region: cloudy region is
    ! homogeneous
    od_scaling(2) = 1.0_jprb
  else
    ! two cloudy regions with optical depth scaled by 1-x and
    ! 1+x.
    ! simple version which fails when fractional_std >= 1:
    !od_scaling(2) = 1.0_jprb-cloud%fractional_std(jcol,jlev)
    ! according to shonk and hogan (2008), 1-x should correspond to
    ! the 16th percentile. 
    if (.not. do_gamma) then
      ! if we treat the distribution as a lognormal such that the
      ! equivalent normal has a mean mu and standard deviation sigma,
      ! then the 16th percentile of the lognormal is very close to
      ! exp(mu-sigma).
      od_scaling(2) &
           &  = exp(-sqrt(log(frac_std**2+1))) / sqrt(frac_std**2+1)
    else
      ! if we treat the distribution as a gamma then the 16th
      ! percentile is close to the following
      od_scaling(2) = exp(-frac_std*(1.0_jprb + 0.5_jprb*frac_std &
           &                                   *(1.0_jprb+0.5_jprb*frac_std)))
    end if

    ! ensure mean optical depth is conserved
    od_scaling(3) = 2.0_jprb-od_scaling(2)
  end if

end subroutine optical_depth_scaling
! # 31 "radiation/radiation_spartacus_lw.f90" 2

  ! this module contains just one subroutine, the longwave solver
  ! using the speedy algorithm for radiative transfer through cloud
  ! sides (spartacus), which can represent 3d effects using a matrix
  ! form of the two-stream equations.
  !
  ! sections:
  !   1: prepare general variables and arrays
  !   2: prepare column-specific variables and arrays
  !   3: first loop over layers
  !     3.1: layer-specific initialization
  !     3.2: compute gamma variables
  !       3.2a: clear-sky case
  !       3.2b: cloudy case
  !     3.3: compute reflection, transmission and emission
  !       3.3a: g-points with 3d effects
  !       3.3b: g-points without 3d effects
  !   4: compute total sources and albedos
  !   5: compute fluxes
  subroutine solver_spartacus_lw(nlev,istartcol,iendcol, &
       &  config, thermodynamics, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
       &  emission, albedo, &
       &  flux)

    use parkind1,                 only : jprb
    use ecradhook,                  only : lhook, dr_hook, jphook

    use radiation_io,             only : nulout
    use radiation_config,         only : config_type, ipdfshapegamma
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_cloud,          only : cloud_type
    use radiation_regions,        only : calc_region_properties
    use radiation_overlap,        only : calc_overlap_matrices
    use radiation_flux,           only : flux_type, indexed_sum
    use radiation_matrix
    use radiation_two_stream,     only : calc_two_stream_gammas_lw, &
         calc_reflectance_transmittance_lw, lwdiffusivitywp
    use radiation_lw_derivatives, only : calc_lw_derivatives_matrix
    use radiation_constants,      only : pi, gasconstantdryair, &
         accelduetogravity

    implicit none

    ! inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in)        :: config
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(cloud_type),   intent(in)       :: cloud

    ! gas and aerosol optical depth of each layer at each longwave
    ! g-point
    real(jprb), intent(in), dimension(config%n_g_lw,nlev,istartcol:iendcol) :: od

    ! gas and aerosol single-scattering albedo and asymmetry factor,
    ! only if longwave scattering by aerosols is to be represented
    real(jprb), intent(in), &
         &  dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol) :: ssa, g

    ! cloud and precipitation optical depth of each layer in each
    ! longwave band
    real(jprb), intent(in) :: od_cloud(config%n_bands_lw,nlev,istartcol:iendcol)

    ! cloud and precipitation single-scattering albedo and asymmetry
    ! factor, only if longwave scattering by clouds is to be
    ! represented
    real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, &
         &                            nlev,istartcol:iendcol) :: ssa_cloud, g_cloud

    ! planck function (emitted flux from a black body) at half levels
    real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) &
         &  :: planck_hl

    ! emission (planck*emissivity) and albedo (1-emissivity) at the
    ! surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) &
         &  :: emission, albedo

    ! output
    type(flux_type),          intent(inout):: flux

    integer :: nreg, ng
    integer :: nregactive ! =1 in clear layer, =nreg in a cloudy layer
    integer :: jcol, jlev, jg, jreg, iband, jreg2
    integer :: ng3d ! number of g-points with small enough gas optical
                    ! depth that 3d effects need to be represented

    ! ratio of gas constant for dry air to acceleration due to gravity
    real(jprb), parameter &
         &  :: r_over_g = gasconstantdryair / accelduetogravity

    ! used in computing rates of lateral radiation transfer
    real(jprb), parameter :: four_over_pi = 4.0_jprb / pi

    ! the tangent of the effective zenith angle for diffuse radiation
    ! that is appropriate for 3d transport
    real(jprb), parameter :: tan_diffuse_angle_3d = pi * 0.5_jprb
    ! incorrect value of tand(53) used by hogan & shonk (2013)
    !    real(jprb), parameter :: tan_diffuse_angle_3d = 1.32704482162041

    ! optical depth, single scattering albedo and asymmetry factor in
    ! each region at each g-point
    real(jprb), dimension(config%n_g_lw, config%nregions) &
         &  :: od_region, ssa_region, g_region

    ! scattering optical depths of gases and clouds
    real(jprb) :: scat_od, scat_od_cloud

    ! the area fractions of each region
    real(jprb) :: region_fracs(1:config%nregions,nlev,istartcol:iendcol)

    ! the scaling used for the optical depth in the cloudy regions
    real(jprb) :: od_scaling(2:config%nregions,nlev,istartcol:iendcol)

    ! the length of the interface between regions per unit area of
    ! gridbox, equal to l_diff^ab in hogan and shonk (2013). this is
    ! actually the effective length oriented to a photon with random
    ! azimuth angle. the three values are the edge length between
    ! regions a-b, b-c and a-c.
    real(jprb) :: edge_length(3)

    ! element i,j gives the rate of 3d transfer of diffuse
    ! radiation from region i to region j, multiplied by the thickness
    ! of the layer in m
    real(jprb) :: transfer_rate(config%nregions,config%nregions)

    ! directional overlap matrices defined at all layer interfaces
    ! including top-of-atmosphere and the surface
    real(jprb), dimension(config%nregions,config%nregions,nlev+1,istartcol:iendcol) &
         &  :: u_matrix, v_matrix

    ! two-stream variables
    real(jprb), dimension(config%n_g_lw, config%nregions) &
         &  :: gamma1, gamma2

    ! matrix gamma multiplied by the layer thickness z1, so units
    ! metres.  after calling expm, this array contains the matrix
    ! exponential of the original.
    real(jprb), dimension(config%n_g_lw, 2*config%nregions, &
         &  2*config%nregions) :: gamma_z1

    ! diffuse reflection and transmission matrices of each layer
    real(jprb), dimension(config%n_g_lw, config%nregions, &
         &  config%nregions, nlev) :: reflectance, transmittance

    ! clear-sky diffuse reflection and transmission matrices of each
    ! layer
    real(jprb), dimension(config%n_g_lw, nlev) :: ref_clear, trans_clear

    ! if the planck term at the top of a layer in each region is the
    ! vector b0 then the following is vector [-b0; b0]*dz
    real(jprb), dimension(config%n_g_lw,2*config%nregions) :: planck_top

    ! the difference between the planck terms at the bottom and top of
    ! a layer, doubled as with planck_top; in terms of b' in the
    ! paper, planck_diff is [-b'; b']*dz*dz
    real(jprb), dimension(config%n_g_lw,2*config%nregions) :: planck_diff

    ! parts of the particular solution associated with the
    ! inhomogeneous ode. in terms of quantities from the paper,
    ! solution0 is [c0;d0] and solution_diff is [c';d']*dz.
    real(jprb), dimension(config%n_g_lw,2*config%nregions) &
         &  :: solution0, solution_diff

    ! used for computing the planck emission per layer
    real(jprb), dimension(config%n_g_lw,config%nregions) &
         &  :: tmp_vectors

    ! the fluxes upwelling from the top of a layer (source_up) and
    ! downwelling from the bottom of the layer (source_dn) due to
    ! emission
    real(jprb), dimension(config%n_g_lw, config%nregions, nlev) &
         &  :: source_up, source_dn
    ! ...clear-sky equivalents
    real(jprb), dimension(config%n_g_lw, nlev) &
         &  :: source_up_clear, source_dn_clear

    ! upwelling radiation just above a layer interface due to emission
    ! below that interface, where level index = 1 corresponds to the
    ! top-of-atmosphere
    real(jprb), dimension(config%n_g_lw, config%nregions, nlev+1) &
         &  :: total_source
    ! ...clear-sky equivalent
    real(jprb), dimension(config%n_g_lw, nlev+1) :: total_source_clear

    ! as total_source, but just below a layer interface
    real(jprb), dimension(config%n_g_lw, config%nregions) &
         &  :: total_source_below

    ! total albedo of the atmosphere/surface just above a layer
    ! interface with respect to downwelling diffuse radiation at that
    ! interface, where level index = 1 corresponds to the
    ! top-of-atmosphere
    real(jprb), dimension(config%n_g_lw, config%nregions, &
         &  config%nregions, nlev+1) :: total_albedo
    ! ...clear-sky equivalent
    real(jprb), dimension(config%n_g_lw, nlev+1) :: total_albedo_clear

    ! as total_albedo, but just below a layer interface
    real(jprb), dimension(config%n_g_lw, config%nregions, config%nregions) &
         &  :: total_albedo_below

    ! the following is used to store matrices of the form i-a*b that
    ! are used on the denominator of some expressions
    real(jprb), dimension(config%n_g_lw, config%nregions, config%nregions) &
         &  :: denominator
    ! clear-sky equivalent, but actually its reciprocal to replace
    ! some divisions by multiplications
    real(jprb), dimension(config%n_g_lw) :: inv_denom_scalar

    ! layer depth (m)
    real(jprb) :: dz

    ! upwelling and downwelling fluxes above/below layer interfaces
    real(jprb), dimension(config%n_g_lw, config%nregions) &
         &  :: flux_up_above, flux_dn_above, flux_dn_below
    ! clear-sky upwelling and downwelling fluxes (which don't
    ! distinguish between whether they are above/below a layer
    ! interface)
    real(jprb), dimension(config%n_g_lw) :: flux_up_clear, flux_dn_clear

    ! parameterization of internal flux distribution in a cloud that
    ! can lead to the flux just about to exit a cloud being different
    ! from the mean in-cloud value.  "lateral_od" is the typical
    ! absorption optical depth through the cloud in a horizontal
    ! direction.
    real(jprb), dimension(config%n_g_lw) :: lateral_od, sqrt_1_minus_ssa
    real(jprb), dimension(config%n_g_lw) :: side_emiss_thick, side_emiss

    ! in the optically thin limit, the flux near the edge is greater
    ! than that in the interior
    real(jprb), parameter :: side_emiss_thin = 1.4107

    real(jprb) :: aspect_ratio

    ! keep a count of the number of calls to the two ways of computing
    ! reflectance/transmittance matrices
    integer :: n_calls_expm, n_calls_meador_weaver

    ! identify clear-sky layers, with pseudo layers for outer space
    ! and below the ground, both treated as single-region clear skies
    logical :: is_clear_sky_layer(0:nlev+1)

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_spartacus_lw:solver_spartacus_lw',0,hook_handle)

    ! --------------------------------------------------------
    ! section 1: prepare general variables and arrays
    ! --------------------------------------------------------

    ! copy array dimensions to local variables for convenience
    nreg = config%nregions
    ng = config%n_g_lw

    ! reset count of number of calls to the two ways to compute
    ! reflectance/transmission matrices
    n_calls_expm = 0
    n_calls_meador_weaver = 0

    ! initialize dz to avoid compiler warning
    dz = 1.0_jprb

    ! compute the wavelength-independent region fractions and
    ! optical-depth scalings
    call calc_region_properties(nlev,nreg,istartcol,iendcol, &
         &  config%i_cloud_pdf_shape == ipdfshapegamma, &
         &  cloud%fraction, cloud%fractional_std, region_fracs, &
         &  od_scaling, config%cloud_fraction_threshold)

    ! compute wavelength-independent overlap matrices u_matrix and v_matrix
    call calc_overlap_matrices(nlev, nreg, istartcol, iendcol, &
         &  region_fracs, cloud%overlap_param, &
         &  u_matrix, v_matrix, decorrelation_scaling=config%cloud_inhom_decorr_scaling, &
         &  cloud_fraction_threshold=config%cloud_fraction_threshold, &
         &  use_beta_overlap=config%use_beta_overlap, &
         &  cloud_cover=flux%cloud_cover_lw)

    if (config%iverbose >= 3) then
      write(nulout,'(a)',advance='no') '  processing columns'
    end if

    ! main loop over columns
    do jcol = istartcol, iendcol
      ! --------------------------------------------------------
      ! section 2: prepare column-specific variables and arrays
      ! --------------------------------------------------------

      if (config%iverbose >= 3) then
        write(nulout,'(a)',advance='no') '.'
      end if

      ! define which layers contain cloud; assume that
      ! cloud%crop_cloud_fraction has already been called
      is_clear_sky_layer = .true.
      do jlev = 1,nlev
        if (cloud%fraction(jcol,jlev) > 0.0_jprb) then
          is_clear_sky_layer(jlev) = .false.
        end if
      end do

      ! --------------------------------------------------------
      ! section 3: first loop over layers
      ! --------------------------------------------------------
      ! in this section the reflectance, transmittance and sources
      ! are computed for each layer
      do jlev = 1,nlev ! start at top-of-atmosphere
        ! --------------------------------------------------------
        ! section 3.1: layer-specific initialization
        ! --------------------------------------------------------
        ! array-wise assignments
        gamma1 = 0.0_jprb
        gamma2 = 0.0_jprb
        gamma_z1= 0.0_jprb
        planck_top = 0.0_jprb
        planck_diff = 0.0_jprb
        transfer_rate(:,:) = 0.0_jprb
        solution0 = 0.0_jprb
        solution_diff = 0.0_jprb

        ! reset clear-sky absorption/scattering properties in each
        ! region
        od_region = 0.0_jprb
        ssa_region = 0.0_jprb ! no rayleigh scattering in longwave by default
        g_region = 0.0_jprb

        ! set optical depth of clear-sky region (region 1) to the
        ! gas/aerosol optical depth
        od_region(1:ng,1) = od(1:ng, jlev, jcol)

        ! can't deal with zero gas optical depths, but this is now
        ! dealt with in radiation_ifsrrtm
        ! where (od_region(:,1) < 1.0e-15_jprb) od_region(:,1) = 1.0e-15_jprb

        ! set clear-sky single-scattering albedo and asymmetry
        ! factor, if non-zero
        if (config%do_lw_aerosol_scattering) then
          ssa_region(1:ng,1) = ssa(1:ng, jlev, jcol)
          g_region(1:ng,1)   = g(1:ng, jlev, jcol)
        end if

        ! --------------------------------------------------------
        ! section 3.2: compute gamma variables
        ! --------------------------------------------------------
        if (is_clear_sky_layer(jlev)) then
          ! --- section 3.2a: clear-sky case --------------------

          nregactive = 1   ! only region 1 (clear-sky) is active

          ! compute two-stream variables gamma1 and gamma2
          call calc_two_stream_gammas_lw(ng, &
               &  ssa_region(1:ng,1), g_region(1:ng,1), &
               &  gamma1(1:ng,1), gamma2(1:ng,1))

          if (config%use_expm_everywhere) then
            ! treat cloud-free layer as 3d: the matrix exponential
            ! "expm" is used as much as possible (provided the
            ! optical depth is not too high). ng3d is initially
            ! set to the total number of g-points, then the
            ! clear-sky optical depths are searched until a value
            ! exceeds the threshold for treating 3d effects, and
            ! if found it is reset to that.  note that we are
            ! assuming that the g-points have been reordered in
            ! approximate order of gas optical depth.
            ng3d = ng
            do jg = 1, ng
              if (od_region(jg,1) > config%max_gas_od_3d) then
                ng3d = jg-1
                exit
              end if
            end do
          else
            ! otherwise treat cloud-free layer using the classical
            ! meador & weaver formulae for all g-points
            ng3d = 0
          end if

        else
          ! --- section 3.2b: cloudy case -----------------------

          ! default number of g-points to treat with
          ! matrix-exponential scheme
          if (config%use_expm_everywhere) then
            ng3d = ng   ! all g-points
          else
            ng3d = 0    ! no g-points
          end if

          ! do we compute 3d effects; note that if we only have one
          ! region and the sky is overcast then 3d calculations must
          ! be turned off as there will be only one region
          if (config%do_3d_effects .and. &
               &  allocated(cloud%inv_cloud_effective_size) .and. &
               &  .not. (nreg == 2 .and. cloud%fraction(jcol,jlev) &
               &  > 1.0_jprb-config%cloud_fraction_threshold)) then
            if (cloud%inv_cloud_effective_size(jcol,jlev) &
                 &  > 0.0_jprb) then
              ! 3d effects are only simulated if
              ! inv_cloud_effective_size is defined and greater
              ! than zero
              ng3d = ng

              ! the following is from the hydrostatic equation
              ! and ideal gas law: dz = dp * r * t / (p * g)
              dz = r_over_g &
                   &  * (thermodynamics%pressure_hl(jcol,jlev+1) &
                   &   - thermodynamics%pressure_hl(jcol,jlev)) &
                   &  * (thermodynamics%temperature_hl(jcol,jlev) &
                   &   + thermodynamics%temperature_hl(jcol,jlev+1)) &
                   &  / (thermodynamics%pressure_hl(jcol,jlev) &
                   &   + thermodynamics%pressure_hl(jcol,jlev+1))

              ! compute cloud edge length per unit area of gridbox
              ! from rearranging hogan & shonk (2013) eq. 45, but
              ! adding a factor of (1-frac) so that a region that
              ! fully occupies the gridbox (frac=1) has an edge of
              ! zero. we can therefore use the fraction of clear sky,
              ! region_fracs(1,jlev,jcol) for convenience instead. the
              ! pi on the denominator means that this is actually edge
              ! length with respect to a light ray with a random
              ! azimuthal direction.
              edge_length(1) = four_over_pi &
                   &  * region_fracs(1,jlev,jcol)*(1.0_jprb-region_fracs(1,jlev,jcol)) &
                   &  * min(cloud%inv_cloud_effective_size(jcol,jlev), &
                   &        1.0_jprb / config%min_cloud_effective_size)
              if (nreg > 2) then
                ! the corresponding edge length between the two cloudy
                ! regions is computed in the same way but treating the
                ! optically denser of the two regions (region 3) as
                ! the cloud; note that the fraction of this region may
                ! be less than that of the optically less dense of the
                ! two regions (region 2).  for increased flexibility,
                ! the user may specify the effective size of
                ! inhomogeneities separately from the cloud effective
                ! size.
                if (allocated(cloud%inv_inhom_effective_size)) then
                  edge_length(2) = four_over_pi &
                       &  * region_fracs(3,jlev,jcol)*(1.0_jprb-region_fracs(3,jlev,jcol)) &
                       &  * min(cloud%inv_inhom_effective_size(jcol,jlev), &
                       &        1.0_jprb / config%min_cloud_effective_size)
                else
                  edge_length(2) = four_over_pi &
                       &  * region_fracs(3,jlev,jcol)*(1.0_jprb-region_fracs(3,jlev,jcol)) &
                       &  * min(cloud%inv_cloud_effective_size(jcol,jlev), &
                       &        1.0_jprb / config%min_cloud_effective_size)
                end if
                ! in the case of three regions, some of the cloud
                ! boundary may go directly to the "thick" region
                if (config%clear_to_thick_fraction > 0.0_jprb) then
                  edge_length(3) = config%clear_to_thick_fraction &
                       &  * min(edge_length(1), edge_length(2))
                  edge_length(1) = edge_length(1) - edge_length(3)
                  edge_length(2) = edge_length(2) - edge_length(3)
                else 
                  edge_length(3) = 0.0_jprb
                end if
              end if

              do jreg = 1, nreg-1
                ! compute lateral transfer rate from region jreg to
                ! jreg+1 following hogan & shonk (2013) eq. 47, but
                ! multiplied by dz because the transfer rate is
                ! vertically integrated across the depth of the layer
                if (region_fracs(jreg,jlev,jcol) > epsilon(1.0_jprb)) then
                  transfer_rate(jreg,jreg+1) = dz &
                       &  * edge_length(jreg) &
                       &  * tan_diffuse_angle_3d / region_fracs(jreg,jlev,jcol)
                end if
                ! compute transfer rate from region jreg+1 to jreg
                if (region_fracs(jreg+1,jlev,jcol) > epsilon(1.0_jprb)) then
                  transfer_rate(jreg+1,jreg) = dz &
                       &  * edge_length(jreg) &
                       &  * tan_diffuse_angle_3d / region_fracs(jreg+1,jlev,jcol)
                end if
              end do

              ! compute transfer rates directly between regions 1 and
              ! 3
              if (edge_length(3) > 0.0_jprb) then
                if (region_fracs(1,jlev,jcol) > epsilon(1.0_jprb)) then
                  transfer_rate(1,3) = dz &
                       &  * edge_length(3) &
                       &  * tan_diffuse_angle_3d / region_fracs(1,jlev,jcol)
                end if
                if (region_fracs(3,jlev,jcol) > epsilon(1.0_jprb)) then
                  transfer_rate(3,1) = dz &
                       &  * edge_length(3) &
                       &  * tan_diffuse_angle_3d / region_fracs(3,jlev,jcol)
                end if
              end if

              ! don't allow the transfer rate out of a region to be
              ! equivalent to a loss of exp(-10) through the layer
              where (transfer_rate > config%max_3d_transfer_rate) 
                transfer_rate = config%max_3d_transfer_rate
              end where
            end if ! cloud has edge length required for 3d effects
          end if ! include 3d effects

          ! in a cloudy layer the number of active regions equals
          ! the number of regions
          nregactive = nreg

          ! compute scattering properties of the regions at each
          ! g-point, mapping from the cloud properties
          ! defined in each band.
          do jg = 1,ng
            ! mapping from g-point to band
            iband = config%i_band_from_reordered_g_lw(jg)

            ! scattering optical depth of clear-sky region
            scat_od = od_region(jg,1)*ssa_region(jg,1)

            ! loop over cloudy regions
            do jreg = 2,nreg
              ! add scaled cloud optical depth to clear-sky value
              od_region(jg,jreg) = od_region(jg,1) &
                   &  + od_cloud(iband,jlev,jcol)*od_scaling(jreg,jlev,jcol)
              if (config%do_lw_cloud_scattering) then
                ! compute single-scattering albedo of gas-cloud
                ! combination
                scat_od_cloud = od_cloud(iband,jlev,jcol) &
                     &  * ssa_cloud(iband,jlev,jcol)*od_scaling(jreg,jlev,jcol)
                ssa_region(jg,jreg) = (scat_od+scat_od_cloud) &
                     &  / od_region(jg,jreg)
                if (scat_od + scat_od_cloud > 0.0_jprb) then
                  ! compute asymmetry factor of gas-cloud
                  ! combination
                  g_region(jg,jreg) &
                       &  = (scat_od*g_region(jg,1) &
                       &  + scat_od_cloud*g_cloud(iband,jlev,jcol)) &
                       &  / (scat_od + scat_od_cloud)
                  ! else g_region is already set to zero
                end if
                ! else ssa_region is already set to zero
              end if

              ! apply maximum cloud optical depth for stability in the
              ! 3d case
              if (od_region(jg,jreg) > config%max_cloud_od) then
                od_region(jg,jreg) = config%max_cloud_od
              end if

            end do

            ! calculate two-stream variables gamma1 and gamma2 of
            ! all regions at once
            call calc_two_stream_gammas_lw(nreg, &
                 &  ssa_region(jg,:), g_region(jg,:), &
                 &  gamma1(jg,:), gamma2(jg,:))

            ! loop is in order of g-points with typically
            ! increasing optical depth: if optical depth of
            ! clear-sky region exceeds a threshold then turn off
            ! 3d effects for any further g-points
            if (ng3d == ng &
                 &  .and. od_region(jg,1) > config%max_gas_od_3d) then
              ng3d = jg-1
            end if
          end do ! loop over g points
        end if ! cloudy level

        ! --------------------------------------------------------
        ! section 3.3: compute reflection, transmission and emission
        ! --------------------------------------------------------
        if (ng3d > 0) then
          ! --- section 3.3a: g-points with 3d effects ----------

          ! 3d effects need to be represented in "ng3d" of the g
          ! points.  this is done by creating ng3d square matrices
          ! each of dimension 2*nreg by 2*nreg, computing the matrix
          ! exponential, then computing the various
          ! transmission/reflectance matrices from that.
          do jreg = 1,nregactive
            ! write the diagonal elements of -gamma1*z1
            gamma_z1(1:ng3d,jreg,jreg) &
                 &  = od_region(1:ng3d,jreg)*gamma1(1:ng3d,jreg)
            ! write the diagonal elements of +gamma2*z1
            gamma_z1(1:ng3d,jreg+nreg,jreg) &
                 &  = od_region(1:ng3d,jreg)*gamma2(1:ng3d,jreg)

            ! write the vectors corresponding to the inhomogeneous
            ! parts of the matrix ode
            planck_top(1:ng3d,nreg+jreg) = od_region(1:ng3d,jreg) &
                 &  *(1.0_jprb-ssa_region(1:ng3d,jreg))*region_fracs(jreg,jlev,jcol) &
                 &  *planck_hl(1:ng3d,jlev,jcol)*lwdiffusivitywp
            planck_top(1:ng3d,jreg) = -planck_top(1:ng3d,nreg+jreg)
            planck_diff(1:ng3d,nreg+jreg) = od_region(1:ng3d,jreg) &
                 &  * (1.0_jprb-ssa_region(1:ng3d,jreg))*region_fracs(jreg,jlev,jcol) &
                 &  * (planck_hl(1:ng3d,jlev+1,jcol) &
                 &  -planck_hl(1:ng3d,jlev,jcol))*lwdiffusivitywp
            planck_diff(1:ng3d,jreg) = -planck_diff(1:ng3d,nreg+jreg)
          end do

          if (nregactive < nreg) then
            ! to avoid nans in some situations we need to put
            ! non-zeros in the cloudy-sky regions even if the cloud
            ! fraction is zero
            do jreg = nregactive+1,nreg
              gamma_z1(1:ng3d,jreg,jreg) = gamma_z1(1:ng3d,1,1)
              gamma_z1(1:ng3d,nreg+jreg,jreg) = gamma_z1(1:ng3d,nreg+1,1)
            end do
          end if

          ! parameterization for the effective emissivity of the side
          ! of the cloud
          if (config%do_lw_side_emissivity &
             & .and. region_fracs(1,jlev,jcol) > 0.0_jprb .and. region_fracs(2,jlev,jcol) > 0.0_jprb &
             & .and. config%do_3d_effects &
             & .and. cloud%inv_cloud_effective_size(jcol,jlev) > 0.0_jprb) then
            aspect_ratio = 1.0_jprb / (min(cloud%inv_cloud_effective_size(jcol,jlev), &
                 &                         1.0_jprb / config%min_cloud_effective_size) &
                 &                     * region_fracs(1,jlev,jcol) * dz)
            lateral_od(1:ng3d) = (aspect_ratio / (nreg-1.0_jprb)) &
                 &  * sum(od_region(1:ng3d,2:nreg)*(1.0_jprb-ssa_region(1:ng3d,2:nreg)),2)
            sqrt_1_minus_ssa(1:ng3d) = sqrt(1.0_jprb - ssa_region(1:ng3d,2))
            side_emiss_thick(1:ng3d) = 2.0_jprb * sqrt_1_minus_ssa(1:ng3d) / &
                 & (sqrt_1_minus_ssa(1:ng3d) &
                 &  + sqrt(1.0_jprb-ssa_region(1:ng3d,2)*g_region(1:ng3d,2)))
            side_emiss(1:ng3d) = (side_emiss_thin - side_emiss_thick(1:ng3d)) &
                 &  / (lateral_od(1:ng3d) + 1.0_jprb) & 
                 &  + side_emiss_thick(1:ng3d)
          else
            side_emiss(1:ng3d) = 1.0_jprb
          end if

          do jreg = 1,nregactive-1
            ! add the terms assocated with 3d transport to gamma1*z1.
            ! first the rate from region jreg to jreg+1
            gamma_z1(1:ng3d,jreg,jreg) = gamma_z1(1:ng3d,jreg,jreg) &
                 &  + transfer_rate(jreg,jreg+1)
            gamma_z1(1:ng3d,jreg+1,jreg) = -transfer_rate(jreg,jreg+1)

            ! next the rate from region jreg+1 to jreg
            if (jreg > 1) then
              ! flow between one cloudy region and another
              gamma_z1(1:ng3d,jreg+1,jreg+1) = gamma_z1(1:ng3d,jreg+1,jreg+1) &
                   &  + transfer_rate(jreg+1,jreg)
              gamma_z1(1:ng3d,jreg,jreg+1) = -transfer_rate(jreg+1,jreg)
            else
              ! only for the lateral transfer between cloud and clear
              ! skies do we account for the effective emissivity of
              ! the side of the cloud
              gamma_z1(1:ng3d,jreg+1,jreg+1) = gamma_z1(1:ng3d,jreg+1,jreg+1) &
                   &  + side_emiss(1:ng3d) * transfer_rate(jreg+1,jreg)
              gamma_z1(1:ng3d,jreg,jreg+1) = -side_emiss(1:ng3d)*transfer_rate(jreg+1,jreg)
            end if
          end do

          ! possible flow between regions a and c
          if (edge_length(3) > 0.0_jprb) then
            gamma_z1(1:ng3d,1,1) = gamma_z1(1:ng3d,1,1) &
                 &  + transfer_rate(1,3)
            gamma_z1(1:ng3d,3,1) = -transfer_rate(1,3)
            gamma_z1(1:ng3d,3,3) = gamma_z1(1:ng3d,3,3) &
                 &  + side_emiss(1:ng3d) * transfer_rate(3,1)
            gamma_z1(1:ng3d,1,3) = -side_emiss(1:ng3d)*transfer_rate(3,1)
          end if

          ! copy gamma1*z1
          gamma_z1(1:ng3d,nreg+1:2*nreg,nreg+1:2*nreg) &
               &  = -gamma_z1(1:ng3d,1:nreg,1:nreg)
          ! copy gamma2*z1
          gamma_z1(1:ng3d,1:nreg,nreg+1:2*nreg) &
               &  = -gamma_z1(1:ng3d,nreg+1:2*nreg,1:nreg)

          ! compute the parts of the particular solution
          solution_diff(1:ng3d,1:2*nreg) &
               &  = solve_vec(ng,ng3d,2*nreg,gamma_z1,planck_diff)
          solution_diff(1:ng3d,1:2*nreg) = - solution_diff(1:ng3d,1:2*nreg) 
          solution0(1:ng3d,1:2*nreg) = solve_vec(ng,ng3d,2*nreg,gamma_z1, &
               &  solution_diff-planck_top)

          ! compute the matrix exponential of gamma_z1, returning the
          ! result in-place
          call expm(ng, ng3d, 2*nreg, gamma_z1, imatrixpatterndense)

          ! update count of expm calls
          n_calls_expm = n_calls_expm + ng3d

          ! diffuse reflectance matrix
          reflectance(1:ng3d,:,:,jlev) = -solve_mat(ng,ng3d,nreg, &
               &  gamma_z1(1:ng3d,1:nreg,1:nreg), &
               &  gamma_z1(1:ng3d,1:nreg,nreg+1:2*nreg))
          ! diffuse transmission matrix
          transmittance(1:ng3d,:,:,jlev) = mat_x_mat(ng,ng3d,nreg, &
               &  gamma_z1(1:ng3d,nreg+1:2*nreg,1:nreg), &
               &  reflectance(1:ng3d,:,:,jlev)) &
               &  + gamma_z1(1:ng3d,nreg+1:2*nreg,nreg+1:2*nreg)

          ! upward and downward sources due to emission within the
          ! layer
          tmp_vectors(1:ng3d,:) = solution0(1:ng3d,1:nreg) &
               &  + solution_diff(1:ng3d,1:nreg) &
               &  - mat_x_vec(ng,ng3d,nreg,gamma_z1(:,1:nreg,nreg+1:2*nreg), &
               &  solution0(:,nreg+1:2*nreg))
          source_up(1:ng3d,:,jlev) = solution0(1:ng3d,1:nreg) &
               &  - solve_vec(ng,ng3d,nreg,gamma_z1(:,1:nreg,1:nreg), &
               &  tmp_vectors)
          tmp_vectors(1:ng3d,:) = source_up(1:ng3d,:,jlev) &
               &  - solution0(1:ng3d,1:nreg)
          source_dn(1:ng3d,:,jlev) = mat_x_vec(ng,ng3d,nreg,&
               &  gamma_z1(:,nreg+1:2*nreg,1:nreg),tmp_vectors) &
               &  + solution0(1:ng3d,nreg+1:2*nreg) &
               &  - mat_x_vec(ng,ng3d,nreg, &
               &  gamma_z1(:,nreg+1:2*nreg,nreg+1:2*nreg), &
               &  solution0(:,nreg+1:2*nreg)) &
               &  + solution_diff(1:ng3d,nreg+1:2*nreg)
        end if ! we are treating 3d effects for some g points

        ! --- section 3.3b: g-points without 3d effects ----------

        ! compute reflectance, transmittance and upward and downward
        ! sources for clear skies, using the meador-weaver formulas
        ! for reflectance and transmittance, and equivalent solutions
        ! to the coupled odes for the sources.
        call calc_reflectance_transmittance_lw(ng, &
             &  od_region(1:ng,1), &
             &  gamma1(1:ng,1), gamma2(1:ng,1), &
             &  planck_hl(1:ng,jlev,jcol), planck_hl(1:ng,jlev+1,jcol), &
             &  ref_clear(1:ng,jlev), trans_clear(1:ng,jlev), &
             &  source_up_clear(1:ng,jlev), source_dn_clear(1:ng,jlev))

        n_calls_meador_weaver = n_calls_meador_weaver + ng

        if (ng3d < ng) then
          ! some of the g points are to be treated using the
          ! conventional plane-parallel method.  first zero the
          ! relevant parts of the arrays
          reflectance  (ng3d+1:ng,:,:,jlev) = 0.0_jprb
          transmittance(ng3d+1:ng,:,:,jlev) = 0.0_jprb
          source_up    (ng3d+1:ng,:,jlev) = 0.0_jprb
          source_dn    (ng3d+1:ng,:,jlev) = 0.0_jprb

          ! since there is no lateral transport, the clear-sky parts
          ! of the arrays can be copied from the clear-sky arrays
          reflectance(ng3d+1:ng,1,1,jlev) = ref_clear(ng3d+1:ng,jlev)
          transmittance(ng3d+1:ng,1,1,jlev) = trans_clear(ng3d+1:ng,jlev)
          source_up(ng3d+1:ng,1,jlev) = region_fracs(1,jlev,jcol)*source_up_clear(ng3d+1:ng,jlev)
          source_dn(ng3d+1:ng,1,jlev) = region_fracs(1,jlev,jcol)*source_dn_clear(ng3d+1:ng,jlev)

          ! compute reflectance, transmittance and upward and downward
          ! sources, for each cloudy region, using the meador-weaver
          ! formulas for reflectance and transmittance, and equivalent
          ! solutions to the coupled odes for the sources.
          do jreg = 2, nregactive
            call calc_reflectance_transmittance_lw(ng-ng3d, &
                 &  od_region(ng3d+1:ng,jreg), &
                 &  gamma1(ng3d+1:ng,jreg), gamma2(ng3d+1:ng,jreg), &
                 &  region_fracs(jreg,jlev,jcol)*planck_hl(ng3d+1:ng,jlev,jcol), &
                 &  region_fracs(jreg,jlev,jcol)*planck_hl(ng3d+1:ng,jlev+1,jcol), &
                 &  reflectance(ng3d+1:ng,jreg,jreg,jlev), &
                 &  transmittance(ng3d+1:ng,jreg,jreg,jlev), &
                 &  source_up(ng3d+1:ng,jreg,jlev), &
                 &  source_dn(ng3d+1:ng,jreg,jlev))
          end do
          n_calls_meador_weaver &
               &  = n_calls_meador_weaver + (ng-ng3d)*(nregactive-1)
        end if

      end do ! loop over levels

      ! --------------------------------------------------------
      ! section 4: compute total sources and albedos
      ! --------------------------------------------------------

      total_albedo(:,:,:,:) = 0.0_jprb
      total_source(:,:,:) = 0.0_jprb

      if (config%do_clear) then
        total_albedo_clear(:,:) = 0.0_jprb
        total_source_clear(:,:) = 0.0_jprb
      end if

      ! calculate the upwelling radiation emitted by the surface, and
      ! copy the surface albedo into total_albedo 
      do jreg = 1,nreg
        do jg = 1,ng
          ! region_fracs(jreg,nlev,jcol) is the fraction of each
          ! region in the lowest model level
          total_source(jg,jreg,nlev+1) = region_fracs(jreg,nlev,jcol)*emission(jg,jcol)
          total_albedo(jg,jreg,jreg,nlev+1) = albedo(jg,jcol)
        end do
      end do
      ! equivalent surface values for computing clear-sky fluxes 
      if (config%do_clear) then
        do jg = 1,ng
          total_source_clear(jg,nlev+1) = emission(jg,jcol)
        end do
        ! in the case of surface albedo there is no dependence on
        ! cloud fraction so we can copy the all-sky value
        total_albedo_clear(1:ng,nlev+1) = total_albedo(1:ng,1,1,nlev+1)
      end if

      ! loop back up through the atmosphere computing the total albedo
      ! and the total upwelling due to emission below each level
      do jlev = nlev,1,-1
        if (config%do_clear) then
          ! use adding method for clear-sky arrays; note that there
          ! is no need to consider "above" and "below" quantities
          ! since with no cloud overlap to worry about, these are
          ! the same
          inv_denom_scalar(:) = 1.0_jprb &
               &  / (1.0_jprb - total_albedo_clear(:,jlev+1)*ref_clear(:,jlev))
          total_albedo_clear(:,jlev) = ref_clear(:,jlev) &
               &  + trans_clear(:,jlev)*trans_clear(:,jlev)*total_albedo_clear(:,jlev+1) &
               &  * inv_denom_scalar(:)
          total_source_clear(:,jlev) = source_up_clear(:,jlev) &
               &  + trans_clear(:,jlev)*(total_source_clear(:,jlev+1) &
               &  + total_albedo_clear(:,jlev+1)*source_dn_clear(:,jlev)) &
               &  * inv_denom_scalar(:)
        end if

        if (is_clear_sky_layer(jlev)) then
          ! clear-sky layer: use scalar adding method
          inv_denom_scalar(:) = 1.0_jprb &
               &  / (1.0_jprb - total_albedo(:,1,1,jlev+1)*reflectance(:,1,1,jlev))
          total_albedo_below = 0.0_jprb
          total_albedo_below(:,1,1) = reflectance(:,1,1,jlev) &
               &  + transmittance(:,1,1,jlev)  * transmittance(:,1,1,jlev) &
               &  * total_albedo(:,1,1,jlev+1) * inv_denom_scalar(:)
          total_source_below = 0.0_jprb
          total_source_below(:,1) = source_up(:,1,jlev) &
               &  + transmittance(:,1,1,jlev)*(total_source(:,1,jlev+1) &
               &  + total_albedo(:,1,1,jlev+1)*source_dn(:,1,jlev)) &
               &  * inv_denom_scalar(:)
        else if (config%do_3d_effects .or. &
             &   config%do_3d_lw_multilayer_effects) then
          ! cloudy layer: use matrix adding method
          denominator = identity_minus_mat_x_mat(ng,ng,nreg,&
               &  total_albedo(:,:,:,jlev+1), reflectance(:,:,:,jlev))
          total_albedo_below = reflectance(:,:,:,jlev) &
               &  + mat_x_mat(ng,ng,nreg,transmittance(:,:,:,jlev), &
               &  solve_mat(ng,ng,nreg,denominator, &
               &  mat_x_mat(ng,ng,nreg,total_albedo(:,:,:,jlev+1), &
               &  transmittance(:,:,:,jlev))))
          total_source_below = source_up(:,:,jlev) & 
               &  + mat_x_vec(ng,ng,nreg,transmittance(:,:,:,jlev), &
               &  solve_vec(ng,ng,nreg,denominator, &
               &  total_source(:,:,jlev+1) + mat_x_vec(ng,ng,nreg, &
               &  total_albedo(:,:,:,jlev+1),source_dn(:,:,jlev))))
        else 
          ! cloudy layer for which reflectance, transmittance and
          ! total_albedo matrices are diagonal
          total_albedo_below = 0.0_jprb
          total_source_below = 0.0_jprb
          do jreg = 1,nreg
            inv_denom_scalar(:) = 1.0_jprb / (1.0_jprb &
                 &  - total_albedo(:,jreg,jreg,jlev+1)*reflectance(:,jreg,jreg,jlev))
            total_albedo_below(:,jreg,jreg) = reflectance(:,jreg,jreg,jlev) &
                 &  + transmittance(:,jreg,jreg,jlev)*transmittance(:,jreg,jreg,jlev) &
                 &  * total_albedo(:,jreg,jreg,jlev+1) &
                 &  * inv_denom_scalar(:)
            total_source_below(:,jreg) = source_up(:,jreg,jlev) &
                 &  + transmittance(:,jreg,jreg,jlev)*(total_source(:,jreg,jlev+1) &
                 &  + total_albedo(:,jreg,jreg,jlev+1)*source_dn(:,jreg,jlev)) &
                 &  * inv_denom_scalar(:)
          end do

        end if

        ! account for cloud overlap when converting albedo and
        ! source below a layer interface to the equivalent values
        ! just above
        if (is_clear_sky_layer(jlev) .and. is_clear_sky_layer(jlev-1)) then
          ! if both layers are cloud free, this is trivial...
          total_albedo(:,:,:,jlev) = 0.0_jprb
          total_albedo(:,1,1,jlev) = total_albedo_below(:,1,1)
          total_source(:,:,jlev) = 0.0_jprb
          total_source(:,1,jlev) = total_source_below(:,1)
        else 
          total_source(:,:,jlev) = singlemat_x_vec(ng,ng,nreg,&
               &  u_matrix(:,:,jlev,jcol), total_source_below)

!          if (config%do_3d_effects .or. config%do_3d_lw_multilayer_effects) then
          if (config%do_3d_lw_multilayer_effects) then
            ! use the overlap matrices u_matrix and v_matrix
            total_albedo(:,:,:,jlev) = singlemat_x_mat(ng,ng,nreg,&
                 &  u_matrix(:,:,jlev,jcol), &
                 &  mat_x_singlemat(ng,ng,nreg,total_albedo_below,&
                 &  v_matrix(:,:,jlev,jcol)))
          else
            total_albedo(:,:,:,jlev) = 0.0_jprb
            ! "total_albedo" is diagonal and we wish to exclude
            ! anomalous horizontal transport described by shonk &
            ! hogan (2008).  therefore, the operation we perform is
            ! essentially diag(total_albedo) = matmul(transpose(v_matrix),
            ! diag(total_albedo_below)).
            do jreg = 1,nreg
              do jreg2 = 1,nreg
                total_albedo(:,jreg,jreg,jlev) &
                     &  = total_albedo(:,jreg,jreg,jlev) &
                     &  + total_albedo_below(:,jreg2,jreg2) &
                     &  * v_matrix(jreg2,jreg,jlev,jcol)
              end do
            end do
          end if
        end if

      end do ! reverse loop over levels

      ! --------------------------------------------------------
      ! section 5: compute fluxes
      ! --------------------------------------------------------

      ! top-of-atmosphere fluxes into the regions of the top-most
      ! layer: zero since we assume no extraterrestrial longwave
      ! radiation (the shortwave scheme accounts for solar radiation
      ! in the "longwave" part of the spectrum)
      flux_dn_below = 0.0_jprb
      flux%lw_dn(jcol,1) = 0.0_jprb
      if (config%do_clear) then
        flux_dn_clear = 0.0_jprb
        flux%lw_dn_clear(jcol,1) = 0.0_jprb
      end if

      ! store the outgoing longwave radiation at top-of-atmosphere
      flux%lw_up(jcol,1) = sum(sum(total_source(:,:,1),1))
      if (config%do_clear) then
        flux%lw_up_clear(jcol,1) = sum(total_source_clear(:,1))
      end if

      if (config%do_save_spectral_flux) then
        call indexed_sum(sum(total_source(:,:,1),2), &
             &           config%i_spec_from_reordered_g_lw, &
             &           flux%lw_up_band(:,jcol,1))
        flux%lw_dn_band(:,jcol,1) = 0.0_jprb
        if (config%do_clear) then
          call indexed_sum(total_source_clear(:,1), &
               &           config%i_spec_from_reordered_g_lw, &
               &           flux%lw_up_clear_band(:,jcol,1))
          flux%lw_dn_clear_band(:,jcol,1) = 0.0_jprb
        end if
      end if

      ! final loop back down through the atmosphere to compute fluxes
      do jlev = 1,nlev
        if (config%do_clear) then
          ! scalar operations for clear-sky fluxes
          flux_dn_clear(:) = (trans_clear(:,jlev)*flux_dn_clear(:) &
               + ref_clear(:,jlev)*total_source_clear(:,jlev+1) &
               + source_dn_clear(:,jlev)) &
               / (1.0_jprb - ref_clear(:,jlev)*total_albedo_clear(:,jlev+1))
          flux_up_clear(:) = total_source_clear(:,jlev+1) &
               + total_albedo_clear(:,jlev+1)*flux_dn_clear(:)
        end if

        if (is_clear_sky_layer(jlev)) then
          ! scalar operations for clear-sky layer
          flux_dn_above(:,1) = (transmittance(:,1,1,jlev)*flux_dn_below(:,1) &
               &  + reflectance(:,1,1,jlev)*total_source(:,1,jlev+1) &
               &  + source_dn(:,1,jlev)) &
               &  / (1.0_jprb - reflectance(:,1,1,jlev)*total_albedo(:,1,1,jlev+1))
          flux_dn_above(:,2:nreg) = 0.0_jprb
          flux_up_above(:,1) = total_source(:,1,jlev+1) &
               &  + total_albedo(:,1,1,jlev+1)*flux_dn_above(:,1)
          flux_up_above(:,2:nreg) = 0.0_jprb
        else if (config%do_3d_effects .or. config%do_3d_lw_multilayer_effects) then
          ! matrix operations for cloudy layer
          denominator = identity_minus_mat_x_mat(ng,ng,nreg,reflectance(:,:,:,jlev), &
               &  total_albedo(:,:,:,jlev+1))
          flux_dn_above = solve_vec(ng,ng,nreg,denominator, &
               &  mat_x_vec(ng,ng,nreg,transmittance(:,:,:,jlev),flux_dn_below) &
               &  + mat_x_vec(ng,ng,nreg,reflectance(:,:,:,jlev), &
               &  total_source(:,:,jlev+1)) &
               &  + source_dn(:,:,jlev))
          flux_up_above = mat_x_vec(ng,ng,nreg,total_albedo(:,:,:,jlev+1), &
               &  flux_dn_above) + total_source(:,:,jlev+1)
        else
          do jreg = 1,nreg
            ! scalar operations for all regions, requiring that
            ! reflectance, transmittance and total_albedo are diagonal
            flux_dn_above(:,jreg) = (transmittance(:,jreg,jreg,jlev)*flux_dn_below(:,jreg) &
                 &  + reflectance(:,jreg,jreg,jlev)*total_source(:,jreg,jlev+1) &
                 &  + source_dn(:,jreg,jlev)) &
                 &  / (1.0_jprb - reflectance(:,jreg,jreg,jlev) &
                 &              * total_albedo(:,jreg,jreg,jlev+1))
            flux_up_above(:,jreg) = total_source(:,jreg,jlev+1) &
                 &  + total_albedo(:,jreg,jreg,jlev+1)*flux_dn_above(:,jreg)
          end do
        end if

        ! account for overlap rules in translating fluxes just above
        ! a layer interface to the values just below
        if (is_clear_sky_layer(jlev) .and. is_clear_sky_layer(jlev+1)) then
          flux_dn_below = flux_dn_above
        else
          flux_dn_below = singlemat_x_vec(ng,ng,nreg,v_matrix(:,:,jlev+1,jcol), &
               &    flux_dn_above)
        end if

        ! store the broadband fluxes
        flux%lw_up(jcol,jlev+1) = sum(sum(flux_up_above,1))
        flux%lw_dn(jcol,jlev+1) = sum(sum(flux_dn_above,1))
        if (config%do_clear) then
          flux%lw_up_clear(jcol,jlev+1) = sum(flux_up_clear)
          flux%lw_dn_clear(jcol,jlev+1) = sum(flux_dn_clear)
        end if

        if (config%do_save_spectral_flux) then
          call indexed_sum(sum(flux_up_above,2), &
               &           config%i_spec_from_reordered_g_lw, &
               &           flux%lw_up_band(:,jcol,jlev+1))
          call indexed_sum(sum(flux_dn_above,2), &
               &           config%i_spec_from_reordered_g_lw, &
               &           flux%lw_dn_band(:,jcol,jlev+1))
          if (config%do_clear) then
            call indexed_sum(flux_up_clear, &
                 &           config%i_spec_from_reordered_g_lw, &
                 &           flux%lw_up_clear_band(:,jcol,jlev+1))
            call indexed_sum(flux_dn_clear, &
                 &           config%i_spec_from_reordered_g_lw, &
                 &           flux%lw_dn_clear_band(:,jcol,jlev+1))
          end if
        end if

      end do ! final loop over levels

      ! store surface spectral downwelling fluxes, which at this point
      ! are at the surface
      flux%lw_dn_surf_g(:,jcol) = sum(flux_dn_above,2)
      if (config%do_clear) then
        flux%lw_dn_surf_clear_g(:,jcol) = flux_dn_clear
      end if

      ! compute the longwave derivatives needed by hogan and bozzo
      ! (2015) approximate radiation update scheme
      if (config%do_lw_derivatives) then
        ! note that at this point flux_up_above contains the spectral
        ! fluxes into the regions of the lowest layer; we sum over
        ! regions first to provide a simple spectral flux upwelling
        ! from the surface
        call calc_lw_derivatives_matrix(ng, nlev, nreg, jcol, transmittance, &
             &  u_matrix(:,:,:,jcol), sum(flux_up_above,2), flux%lw_derivatives)
      end if
        
    end do ! loop over columns

    if (config%iverbose >= 3) then
      write(nulout,*)
    end if

    ! report number of calls to each method of solving single-layer
    ! two-stream equations
    if (config%iverbose >= 4) then
      write(nulout,'(a,i0)') '  matrix-exponential calls: ', n_calls_expm
      write(nulout,'(a,i0)') '  meador-weaver calls: ', n_calls_meador_weaver
    end if

    if (lhook) call dr_hook('radiation_spartacus_lw:solver_spartacus_lw',1,hook_handle)

  end subroutine solver_spartacus_lw

end module radiation_spartacus_lw
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

