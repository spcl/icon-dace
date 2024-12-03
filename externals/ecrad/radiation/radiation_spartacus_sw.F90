! # 1 "radiation/radiation_spartacus_sw.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_spartacus_sw.f90"
! radiation_spartacus_sw.f90 - spartacus shortwave solver
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
!   2017-04-11  r. hogan  receive albedos at g-points
!   2017-04-22  r. hogan  store surface fluxes at all g-points
!   2017-06-30  r. hogan  reformulate to use total_albedo_direct not total_source
!   2017-07-03  r. hogan  explicit calculation of encroachment
!   2017-10-23  r. hogan  renamed single-character variables
!   2018-02-20  r. hogan  corrected "computed" encroachment
!   2018-03-09  r. hogan  security on computed encroachment, transmittance and reflectance
!   2018-08-29  r. hogan  reformulate horizontal migration distances in step_migrations
!   2018-09-02  r. hogan  bug fix in x_direct in step_migrations
!   2018-09-03  r. hogan  security via min_cloud_effective_size
!   2018-09-04  r. hogan  use encroachment_scaling and read encroachment edge_length from upper layer
!   2018-09-13  r. hogan  added "fractal" encroachment option
!   2018-09-14  r. hogan  added "zero" encroachment option
!   2018-10-08  r. hogan  call calc_region_properties
!   2018-10-15  r. hogan  added call to fast_expm_exchange instead of expm
!   2019-01-12  r. hogan  use inv_inhom_effective_size if allocated
!   2019-02-10  r. hogan  renamed "encroachment" to "entrapment"

module radiation_spartacus_sw

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
! # 43 "radiation/radiation_spartacus_sw.f90" 2

  ! this module contains just one exported subroutine, the shortwave
  ! solver using the speedy algorithm for radiative transfer through
  ! cloud sides (spartacus), which can represent 3d effects using a
  ! matrix form of the two-stream equations.
  !
  ! sections:
  !   1: prepare general variables and arrays
  !   2: prepare column-specific variables and arrays
  !   3: first loop over layers
  !     3.1: layer-specific initialization
  !     3.2: compute gamma variables
  !       3.2a: clear-sky case
  !       3.2b: cloudy case
  !     3.3: compute reflection, transmission and sources
  !       3.3a: g-points with 3d effects
  !       3.3b: g-points without 3d effects
  !   4: compute total albedos
  !     4.1: adding method
  !     4.2: overlap and entrapment
  !   5: compute fluxes
  subroutine solver_spartacus_sw(nlev,istartcol,iendcol, &
       &  config, single_level, thermodynamics, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, &
       &  albedo_direct, albedo_diffuse, incoming_sw, &
       &  flux)

    use parkind1, only           : jprb
    use ecradhook,  only           : lhook, dr_hook, jphook

    use radiation_io, only             : nulout
    use radiation_config, only         : config_type, ipdfshapegamma, &
         &  ientrapmentzero, ientrapmentedgeonly, ientrapmentexplicit, &
         &  ientrapmentexplicitnonfractal, ientrapmentmaximum
    use radiation_single_level, only   : single_level_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_cloud, only          : cloud_type
    use radiation_regions, only        : calc_region_properties
    use radiation_overlap, only        : calc_overlap_matrices
    use radiation_flux, only           : flux_type, &
         &                               indexed_sum, add_indexed_sum
    use radiation_matrix
    use radiation_two_stream, only     : calc_two_stream_gammas_sw, &
         &  calc_reflectance_transmittance_sw, calc_frac_scattered_diffuse_sw
    use radiation_constants, only      : pi, gasconstantdryair, &
         &                               accelduetogravity

    implicit none

    ! inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(cloud_type),         intent(in) :: cloud

    ! gas and aerosol optical depth, single-scattering albedo and
    ! asymmetry factor at each shortwave g-point
    real(jprb), intent(in), dimension(config%n_g_sw,nlev,istartcol:iendcol) :: &
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

    integer :: nreg, ng
    integer :: nregactive ! =1 in clear layer, =nreg in a cloudy layer
    integer :: jcol, jlev, jg, jreg, iband, jreg2,jreg3



    integer :: ng3d ! number of g-points with small enough gas optical
                    ! depth that 3d effects need to be represented

    ! ratio of gas constant for dry air to acceleration due to gravity
    real(jprb), parameter :: r_over_g = gasconstantdryair / accelduetogravity

    real(jprb) :: mu0, one_over_mu0, tan_sza, transfer_scaling

    ! the tangent of the effective zenith angle for diffuse radiation
    ! that is appropriate for 3d transport
    real(jprb), parameter :: tan_diffuse_angle_3d = pi * 0.5_jprb

    ! the minimum cosine of solar zenith angle to consider for 3d
    ! effects, equivalent to one solar radius (0.2615 degrees)
    real(jprb), parameter :: min_mu0_3d = 0.004625_jprb

    ! optical depth, single scattering albedo and asymmetry factor in
    ! each region and (except for asymmetry) at each g-point
    real(jprb), dimension(config%n_g_sw, config%nregions) &
         &  :: od_region, ssa_region
    real(jprb), dimension(config%nregions) :: g_region

    ! scattering optical depths of gases and clouds
    real(jprb) :: scat_od, scat_od_cloud

    ! the area fractions of each region
    real(jprb) :: region_fracs(1:config%nregions,nlev,istartcol:iendcol)

    ! the scaling used for the optical depth in the cloudy regions
    real(jprb) :: od_scaling(2:config%nregions,nlev,istartcol:iendcol)

    ! the length of the interface between regions jreg and jreg+1 per
    ! unit area of gridbox, equal to l_diff^ab in hogan and shonk
    ! (2013). when the first index is 3, this is the length of the
    ! interface between regions 3 and 1. this is actually the
    ! effective length oriented to a photon with random azimuth angle,
    ! so is the true edge length divided by pi.
    real(jprb) :: edge_length(3,nlev)

    ! element i,j gives the rate of 3d transfer of diffuse/direct
    ! radiation from region i to region j, multiplied by the thickness
    ! of the layer in m
    real(jprb) :: transfer_rate_diffuse(config%nregions,config%nregions)
    real(jprb) :: transfer_rate_direct(config%nregions,config%nregions)

    ! directional overlap matrices defined at all layer interfaces
    ! including top-of-atmosphere and the surface
    real(jprb), dimension(config%nregions,config%nregions,nlev+1, &
         &                istartcol:iendcol) :: u_matrix, v_matrix

    ! two-stream variables
    real(jprb), dimension(config%n_g_sw, config%nregions) &
         &  :: gamma1, gamma2, gamma3

    ! matrix gamma multiplied by the layer thickness z1, so units
    ! metres.  after calling expm, this array contains the matrix
    ! exponential of the original.
    real(jprb) :: gamma_z1(config%n_g_sw,3*config%nregions,3*config%nregions)

    ! diffuse reflection and transmission matrices of each layer
    real(jprb), dimension(config%n_g_sw, config%nregions, &
         &  config%nregions, nlev) :: reflectance, transmittance

    ! clear-sky diffuse reflection and transmission matrices of each
    ! layer
    real(jprb), dimension(config%n_g_sw, nlev) :: ref_clear, trans_clear

    ! matrices translating the direct flux entering the layer from
    ! above to the reflected radiation exiting upwards (ref_dir) and
    ! the scattered radiation exiting downwards (trans_dir_diff),
    ! along with the direct unscattered transmission matrix
    ! (trans_dir_dir).
    real(jprb), dimension(config%n_g_sw, config%nregions, config%nregions, nlev) &
         &  :: ref_dir, trans_dir_diff, trans_dir_dir
    ! ...clear-sky equivalents
    real(jprb), dimension(config%n_g_sw, nlev) &
         &  :: ref_dir_clear, trans_dir_diff_clear, trans_dir_dir_clear

    ! the fluxes downwelling from the bottom of the layer due to
    ! scattering by the direct beam within the layer
    real(jprb), dimension(config%n_g_sw, config%nregions) :: source_dn
    ! ...clear-sky equivalent
    real(jprb), dimension(config%n_g_sw) :: source_dn_clear

    ! the fluxes upwelling just above the base of a layer due to
    ! reflection of the direct downwelling beam; this is just used as
    ! a temporary variable
    real(jprb), dimension(config%n_g_sw, config%nregions) :: total_source

    ! direct downwelling flux below and above an interface between
    ! layers into a plane perpendicular to the direction of the sun
    real(jprb), dimension(config%n_g_sw, config%nregions) &
         &  :: direct_dn_below, direct_dn_above
    ! ...clear-sky equivalent (no distinction between "above/below")
    real(jprb), dimension(config%n_g_sw) :: direct_dn_clear

    ! total albedo of the atmosphere/surface just above a layer
    ! interface with respect to downwelling diffuse and direct
    ! radiation at that interface, where level index = 1 corresponds
    ! to the top-of-atmosphere
    real(jprb), dimension(config%n_g_sw, config%nregions, &
         &  config%nregions, nlev+1) :: total_albedo, total_albedo_direct
    ! ...clear-sky equivalent
    real(jprb), dimension(config%n_g_sw, nlev+1) &
         &  :: total_albedo_clear, total_albedo_clear_direct

    ! as total_albedo, but just below a layer interface
    real(jprb), dimension(config%n_g_sw, config%nregions, config%nregions) &
         &  :: total_albedo_below, total_albedo_below_direct

    ! temporary array for applying adding method with entrapment to
    ! albedo matrices
    real(jprb), dimension(config%n_g_sw, config%nregions, config%nregions) &
         &  :: albedo_part

    ! horizontal migration distance (m) of reflected light
    real(jprb), dimension(config%n_g_sw, config%nregions) &
         &  :: x_diffuse, x_direct
    ! temporary variables when applying overlap rules
    real(jprb), dimension(config%n_g_sw, config%nregions) &
         &  :: x_diffuse_above, x_direct_above

    real(jprb), dimension(config%n_g_sw, config%nregions, config%nregions) &
         &  :: entrapment

    ! the following is used to store matrices of the form i-a*b that
    ! are used on the denominator of some expressions
    real(jprb) :: denominator(config%n_g_sw,config%nregions,config%nregions)

    ! clear-sky equivalent, but actually its reciprocal to replace
    ! some divisions by multiplications
    real(jprb), dimension(config%n_g_sw) :: inv_denom_scalar

    ! final step in working out how much transport between regions
    ! above occurs
    real(jprb), dimension(config%n_g_sw) :: fractal_factor

    ! inverse of cloud effective size (m^-1)
    real(jprb) :: inv_effective_size

    ! layer depth (m)
    real(jprb) :: dz, layer_depth(nlev)

    ! upwelling and downwelling fluxes above/below layer interfaces
    real(jprb), dimension(config%n_g_sw, config%nregions) &
         &  :: flux_up_above, flux_dn_above, flux_dn_below
    ! clear-sky upwelling and downwelling fluxes (which don't
    ! distinguish between whether they are above/below a layer
    ! interface)
    real(jprb), dimension(config%n_g_sw) :: flux_up_clear, flux_dn_clear

    ! index of top-most cloudy layer, or nlev+1 if no cloud
    integer :: i_cloud_top

    ! keep a count of the number of calls to the two ways of computing
    ! reflectance/transmittance matrices
    integer :: n_calls_expm, n_calls_meador_weaver

    ! identify clear-sky layers, with pseudo layers for outer space
    ! and below the ground, both treated as single-region clear skies
    logical :: is_clear_sky_layer(0:nlev+1)

    ! used in computing rates of lateral radiation transfer
    real(jprb), parameter :: four_over_pi = 4.0_jprb / pi

    ! maximum entrapment coefficient
    real(jprb) :: max_entr

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_spartacus_sw:solver_spartacus_sw',0,hook_handle)

    ! --------------------------------------------------------
    ! section 1: prepare general variables and arrays
    ! --------------------------------------------------------

    ! copy array dimensions to local variables for convenience
    nreg = config%nregions
    ng   = config%n_g_sw

    ! reset count of number of calls to the two ways to compute
    ! reflectance/transmission matrices
    n_calls_expm          = 0
    n_calls_meador_weaver = 0

    ! compute the wavelength-independent region fractions and
    ! optical-depth scalings
    call calc_region_properties(nlev,nreg,istartcol,iendcol, &
         &  config%i_cloud_pdf_shape == ipdfshapegamma, &
         &  cloud%fraction, cloud%fractional_std, region_fracs, &
         &  od_scaling, config%cloud_fraction_threshold)

    ! compute wavelength-independent overlap matrices u_matrix and v_matrix
    call calc_overlap_matrices(nlev,nreg,istartcol,iendcol, &
         &  region_fracs, cloud%overlap_param, &
         &  u_matrix, v_matrix, decorrelation_scaling=config%cloud_inhom_decorr_scaling, &
         &  cloud_fraction_threshold=config%cloud_fraction_threshold, &
         &  use_beta_overlap=config%use_beta_overlap, &
         &  cloud_cover=flux%cloud_cover_sw)

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

      ! copy local cosine of the solar zenith angle
      mu0 = single_level%cos_sza(jcol)

      ! skip profile if sun is too low in the sky
      if (mu0 < 1.0e-10_jprb) then
        flux%sw_dn(jcol,:) = 0.0_jprb
        flux%sw_up(jcol,:) = 0.0_jprb
        if (allocated(flux%sw_dn_direct)) then
          flux%sw_dn_direct(jcol,:) = 0.0_jprb
        end if
        if (config%do_clear) then
          flux%sw_dn_clear(jcol,:) = 0.0_jprb
          flux%sw_up_clear(jcol,:) = 0.0_jprb
          if (allocated(flux%sw_dn_direct_clear)) then
            flux%sw_dn_direct_clear(jcol,:) = 0.0_jprb
          end if
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

        flux%sw_dn_diffuse_surf_g(:,jcol) = 0.0_jprb
        flux%sw_dn_direct_surf_g(:,jcol)  = 0.0_jprb
        if (config%do_clear) then
          flux%sw_dn_diffuse_surf_clear_g(:,jcol) = 0.0_jprb
          flux%sw_dn_direct_surf_clear_g(:,jcol)  = 0.0_jprb
        end if

        cycle
      end if ! sun is below the horizon

      ! at this point mu0 >= 1.0e-10

      ! used to compute rate of attenuation of direct solar beam
      ! through the atmosphere
      one_over_mu0 = 1.0_jprb / mu0

      ! the rate at which direct radiation enters cloud sides is
      ! proportional to the tangent of the solar zenith angle, but
      ! this gets very large when the sun is low in the sky, in which
      ! case we limit it to one solar radius above the horizon
      if (mu0 < min_mu0_3d) then
        tan_sza = sqrt(1.0_jprb/(min_mu0_3d*min_mu0_3d) - 1.0_jprb)
      else if (one_over_mu0 > 1.0_jprb) then
        tan_sza = sqrt(one_over_mu0*one_over_mu0 - 1.0_jprb &
             &         + config%overhead_sun_factor)
      else
        ! just in case we get mu0 > 1...
        tan_sza = sqrt(config%overhead_sun_factor)
      end if

      ! define which layers contain cloud; assume that
      ! cloud%crop_cloud_fraction has already been called
      is_clear_sky_layer = .true.
      i_cloud_top = nlev+1
      do jlev = nlev,1,-1
        if (cloud%fraction(jcol,jlev) > 0.0_jprb) then
          is_clear_sky_layer(jlev) = .false.
          i_cloud_top = jlev
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
        gamma3 = 0.0_jprb
        gamma_z1= 0.0_jprb
        transfer_rate_direct(:,:)  = 0.0_jprb
        transfer_rate_diffuse(:,:) = 0.0_jprb
        edge_length(:,jlev) = 0.0_jprb

        ! the following is from the hydrostatic equation
        ! and ideal gas law: dz = dp * r * t / (p * g)
        layer_depth(jlev) = r_over_g &
             &  * (thermodynamics%pressure_hl(jcol,jlev+1) &
             &     - thermodynamics%pressure_hl(jcol,jlev)) &
             &  * (thermodynamics%temperature_hl(jcol,jlev) &
             &     + thermodynamics%temperature_hl(jcol,jlev+1)) &
             &  / (thermodynamics%pressure_hl(jcol,jlev) &
             &     + thermodynamics%pressure_hl(jcol,jlev+1))

        ! --------------------------------------------------------
        ! section 3.2: compute gamma variables
        ! --------------------------------------------------------
        if (is_clear_sky_layer(jlev)) then
          ! --- section 3.2a: clear-sky case --------------------

          nregactive = 1   ! only region 1 (clear-sky) is active

          ! copy optical depth and single-scattering albedo of
          ! clear-sky region
          od_region(1:ng,1) = od(1:ng,jlev,jcol)
          ssa_region(1:ng,1) = ssa(1:ng,jlev,jcol)

          ! compute two-stream variables gamma1-gamma3 (gamma4=1-gamma3)
          call calc_two_stream_gammas_sw(ng, &
               &  mu0, ssa(1:ng,jlev,jcol), g(1:ng,jlev,jcol), &
               &  gamma1(1:ng,1), gamma2(1:ng,1), gamma3(1:ng,1))

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

          if (config%do_3d_effects .and. &
               &  allocated(cloud%inv_cloud_effective_size) .and. &
               &  .not. (nreg == 2 .and. cloud%fraction(jcol,jlev) &
               &  > 1.0-config%cloud_fraction_threshold)) then
            if (cloud%inv_cloud_effective_size(jcol,jlev) > 0.0_jprb) then
              ! 3d effects are only simulated if
              ! inv_cloud_effective_size is defined and greater
              ! than zero
              ng3d = ng

              ! depth of current layer
              dz = layer_depth(jlev)

              ! compute cloud edge length per unit area of gridbox
              ! from rearranging hogan & shonk (2013) eq. 45, but
              ! adding a factor of (1-frac) so that a region that
              ! fully occupies the gridbox (frac=1) has an edge of
              ! zero. we can therefore use the fraction of clear sky,
              ! region_fracs(1,jlev,jcol) for convenience instead. the
              ! pi on the denominator means that this is actually edge
              ! length with respect to a light ray with a random
              ! azimuthal direction.
              edge_length(1,jlev) = four_over_pi &
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
                  edge_length(2,jlev) = four_over_pi &
                       &  * region_fracs(3,jlev,jcol)*(1.0_jprb-region_fracs(3,jlev,jcol)) &
                       &  * min(cloud%inv_inhom_effective_size(jcol,jlev), &
                       &        1.0_jprb / config%min_cloud_effective_size)
                else
                  edge_length(2,jlev) = four_over_pi &
                       &  * region_fracs(3,jlev,jcol)*(1.0_jprb-region_fracs(3,jlev,jcol)) &
                       &  * min(cloud%inv_cloud_effective_size(jcol,jlev), &
                       &        1.0_jprb / config%min_cloud_effective_size)
                end if

                ! in the case of three regions, some of the cloud
                ! boundary may go directly to the "thick" region
                if (config%clear_to_thick_fraction > 0.0_jprb) then
                  edge_length(3,jlev) = config%clear_to_thick_fraction &
                       &  * min(edge_length(1,jlev), edge_length(2,jlev))
                  edge_length(1,jlev) = edge_length(1,jlev) - edge_length(3,jlev)
                  edge_length(2,jlev) = edge_length(2,jlev) - edge_length(3,jlev)
                else 
                  edge_length(3,jlev) = 0.0_jprb
                end if
              end if

              do jreg = 1, nreg-1
                ! compute lateral transfer rates from region jreg to
                ! jreg+1 following hogan & shonk (2013) eq. 47, but
                ! multiplied by dz because the transfer rate is
                ! vertically integrated across the depth of the layer
                if (region_fracs(jreg,jlev,jcol) > epsilon(1.0_jprb)) then
                  transfer_rate_direct(jreg,jreg+1) = dz &
                       &  * edge_length(jreg,jlev) * tan_sza / region_fracs(jreg,jlev,jcol)
                  transfer_rate_diffuse(jreg,jreg+1) = dz &
                       &  * edge_length(jreg,jlev) &
                       &  * tan_diffuse_angle_3d / region_fracs(jreg,jlev,jcol)
                end if
                ! compute transfer rates from region jreg+1 to
                ! jreg
                if (region_fracs(jreg+1,jlev,jcol) > epsilon(1.0_jprb)) then
                  transfer_rate_direct(jreg+1,jreg) = dz &
                       &  * edge_length(jreg,jlev) &
                       &  * tan_sza / region_fracs(jreg+1,jlev,jcol)
                  transfer_rate_diffuse(jreg+1,jreg) = dz &
                       &  * edge_length(jreg,jlev) &
                       &  * tan_diffuse_angle_3d / region_fracs(jreg+1,jlev,jcol)
                end if
              end do

              ! compute transfer rates directly between regions 1 and
              ! 3
              if (edge_length(3,jlev) > 0.0_jprb) then
                if (region_fracs(1,jlev,jcol) > epsilon(1.0_jprb)) then
                  transfer_rate_direct(1,3) = dz &
                       &  * edge_length(3,jlev) * tan_sza / region_fracs(1,jlev,jcol)
                  transfer_rate_diffuse(1,3) = dz &
                       &  * edge_length(3,jlev) &
                       &  * tan_diffuse_angle_3d / region_fracs(1,jlev,jcol)
                end if
                if (region_fracs(3,jlev,jcol) > epsilon(1.0_jprb)) then
                  transfer_rate_direct(3,1) = dz &
                       &  * edge_length(3,jlev) * tan_sza / region_fracs(3,jlev,jcol)
                  transfer_rate_diffuse(3,1) = dz &
                       &  * edge_length(3,jlev) &
                       &  * tan_diffuse_angle_3d / region_fracs(3,jlev,jcol)
                end if
              end if

              ! don't allow the transfer rate out of a region to be
              ! equivalent to a loss of exp(-10) through the layer
              where (transfer_rate_direct > config%max_3d_transfer_rate) 
                transfer_rate_direct = config%max_3d_transfer_rate
              end where
              where (transfer_rate_diffuse > config%max_3d_transfer_rate) 
                transfer_rate_diffuse = config%max_3d_transfer_rate
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
            iband = config%i_band_from_reordered_g_sw(jg)

            ! scattering optical depth of clear-sky region
            scat_od = od(jg,jlev,jcol)*ssa(jg,jlev,jcol)

            ! scattering properties of clear-sky regions copied
            ! over
            od_region(jg,1)  = od(jg, jlev, jcol)
            ssa_region(jg,1) = ssa(jg, jlev, jcol)
            g_region(1)      = g(jg, jlev, jcol)

            ! loop over cloudy regions
            do jreg = 2,nreg
              scat_od_cloud = od_cloud(iband,jlev,jcol) &
                   &  * ssa_cloud(iband,jlev,jcol)*od_scaling(jreg,jlev,jcol)
              ! add scaled cloud optical depth to clear-sky value
              od_region(jg,jreg) = od(jg,jlev,jcol) &
                   &  + od_cloud(iband,jlev,jcol)*od_scaling(jreg,jlev,jcol)
              ! compute single-scattering albedo and asymmetry
              ! factor of gas-cloud combination
              ssa_region(jg,jreg) = (scat_od+scat_od_cloud) &
                   &  / od_region(jg,jreg)
              g_region(jreg) = (scat_od*g(jg,jlev,jcol) &
                   &  + scat_od_cloud * g_cloud(iband,jlev,jcol)) &
                   &  / (scat_od + scat_od_cloud)

              ! apply maximum cloud optical depth for stability in the
              ! 3d case
              if (od_region(jg,jreg) > config%max_cloud_od) then
                od_region(jg,jreg) = config%max_cloud_od
              end if

            end do

            ! calculate two-stream variables gamma1-gamma3 of all
            ! regions at once
            call calc_two_stream_gammas_sw(nreg, &
                 &  mu0, ssa_region(jg,:), g_region, &
                 &  gamma1(jg,:), gamma2(jg,:), gamma3(jg,:))

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
          ! each of dimension 3*nreg by 3*nreg, computing the matrix
          ! exponential, then computing the various
          ! transmission/reflectance matrices from that.
          do jreg = 1,nregactive
            ! write the diagonal elements of -gamma1*z1
            gamma_z1(1:ng3d,jreg,jreg) &
                 &  = od_region(1:ng3d,jreg)*gamma1(1:ng3d,jreg)
            ! write the diagonal elements of +gamma2*z1
            gamma_z1(1:ng3d,jreg+nreg,jreg) &
                 &  = od_region(1:ng3d,jreg)*gamma2(1:ng3d,jreg)
            ! write the diagonal elements of -gamma3*z1
            gamma_z1(1:ng3d,jreg,jreg+2*nreg) &
                 &  = -od_region(1:ng3d,jreg)*ssa_region(1:ng3d,jreg) &
                 &  * gamma3(1:ng3d,jreg)

            ! write the diagonal elements of +gamma4*z1
            gamma_z1(1:ng3d,jreg+nreg,jreg+2*nreg) &
                 &  = od_region(1:ng3d,jreg)*ssa_region(1:ng3d,jreg) &
                 &  * (1.0_jprb - gamma3(1:ng3d,jreg))

            ! write the diagonal elements of +gamma0*z1
            gamma_z1(1:ng3d,jreg+2*nreg,jreg+2*nreg) &
                 &  = -od_region(1:ng3d,jreg)*one_over_mu0
          end do

          do jreg = 1,nregactive-1
            ! write the elements of -gamma1*z1 concerned with 3d
            ! transport
            gamma_z1(1:ng3d,jreg,jreg) = gamma_z1(1:ng3d,jreg,jreg) &
                 &  + transfer_rate_diffuse(jreg,jreg+1)
            gamma_z1(1:ng3d,jreg+1,jreg+1) = gamma_z1(1:ng3d,jreg+1,jreg+1) &
                 &  + transfer_rate_diffuse(jreg+1,jreg)
            gamma_z1(1:ng3d,jreg+1,jreg) = -transfer_rate_diffuse(jreg,jreg+1)
            gamma_z1(1:ng3d,jreg,jreg+1) = -transfer_rate_diffuse(jreg+1,jreg)
            ! write the elements of +gamma0*z1 concerned with 3d
            ! transport
            gamma_z1(1:ng3d,jreg+2*nreg,jreg+2*nreg) &
                 &  = gamma_z1(1:ng3d,jreg+2*nreg,jreg+2*nreg) &
                 &  - transfer_rate_direct(jreg,jreg+1)
            gamma_z1(1:ng3d,jreg+2*nreg+1,jreg+2*nreg+1) &
                 &  = gamma_z1(1:ng3d,jreg+2*nreg+1,jreg+2*nreg+1) &
                 &  - transfer_rate_direct(jreg+1,jreg)
            gamma_z1(1:ng3d,jreg+2*nreg+1,jreg+2*nreg) &
                 &  = transfer_rate_direct(jreg,jreg+1)
            gamma_z1(1:ng3d,jreg+2*nreg,jreg+2*nreg+1) &
                 &  = transfer_rate_direct(jreg+1,jreg)
          end do

          ! possible flow between regions a and c
          if (edge_length(3,jlev) > 0.0_jprb) then
            ! diffuse transport
            gamma_z1(1:ng3d,1,1) = gamma_z1(1:ng3d,1,1) &
                 &  + transfer_rate_diffuse(1,3)
            gamma_z1(1:ng3d,3,3) = gamma_z1(1:ng3d,3,3) &
                 &  + transfer_rate_diffuse(3,1)
            gamma_z1(1:ng3d,3,1) = -transfer_rate_diffuse(1,3)
            gamma_z1(1:ng3d,1,3) = -transfer_rate_diffuse(3,1)
            ! direct transport
            gamma_z1(1:ng3d,1+2*nreg,1+2*nreg) = gamma_z1(1:ng3d,1+2*nreg,1+2*nreg) &
                 &  - transfer_rate_direct(1,3)
            gamma_z1(1:ng3d,3+2*nreg,3+2*nreg) = gamma_z1(1:ng3d,3+2*nreg,3+2*nreg) &
                 &  - transfer_rate_direct(3,1)
            gamma_z1(1:ng3d,3+2*nreg,1+2*nreg) = transfer_rate_direct(1,3)
            gamma_z1(1:ng3d,1+2*nreg,3+2*nreg) = transfer_rate_direct(3,1)
          end if

          ! copy gamma1*z1
          gamma_z1(1:ng3d,nreg+1:nreg+nregactive,nreg+1:nreg+nregactive) &
               &  = -gamma_z1(1:ng3d,1:nregactive,1:nregactive)
          ! copy gamma2*z1
          gamma_z1(1:ng3d,1:nregactive,nreg+1:nreg+nregactive) &
               &  = -gamma_z1(1:ng3d,nreg+1:nreg+nregactive,1:nregactive)

          ! compute the matrix exponential of gamma_z1, returning the
          ! result in-place
          call expm(ng, ng3d, 3*nreg, gamma_z1, imatrixpatternshortwave)

          ! update count of expm calls
          n_calls_expm = n_calls_expm + ng3d

          ! direct transmission matrix
          trans_dir_dir(1:ng3d,:,:,jlev) = min(1.0_jprb,max(0.0_jprb, &
               &  gamma_z1(1:ng3d,2*nreg+1:3*nreg, 2*nreg+1:3*nreg)))
          ! diffuse reflectance matrix; security on negative values
          ! necessary occasionally for very low cloud fraction and very high
          ! in-cloud optical depth
          reflectance(1:ng3d,:,:,jlev) = min(1.0_jprb,max(0.0_jprb, &
               &  -solve_mat(ng,ng3d,nreg,gamma_z1(1:ng3d,1:nreg,1:nreg), &
               &             gamma_z1(1:ng3d,1:nreg,nreg+1:2*nreg))))
          ! diffuse transmission matrix
          transmittance(1:ng3d,:,:,jlev) = min(1.0_jprb,max(0.0_jprb, &
               &  mat_x_mat(ng,ng3d,nreg,gamma_z1(1:ng3d,nreg+1:2*nreg,1:nreg), &
               &            reflectance(1:ng3d,:,:,jlev)) &
               &  + gamma_z1(1:ng3d,nreg+1:2*nreg,nreg+1:2*nreg)))
          ! transfer matrix between downward direct and upward
          ! diffuse
          ref_dir(1:ng3d,:,:,jlev) = min(mu0,max(0.0_jprb, &
               &  -solve_mat(ng,ng3d,nreg,gamma_z1(1:ng3d,1:nreg,1:nreg), &
               &             gamma_z1(1:ng3d,1:nreg,2*nreg+1:3*nreg))))
          ! transfer matrix between downward direct and downward
          ! diffuse in layer interface below.  include correction for
          ! trans_dir_diff out of plausible bounds (note that meador &
          ! weaver has the same correction in radiation_two_stream.f90
          ! - this is not just an expm thing)
          trans_dir_diff(1:ng3d,:,:,jlev) = min(mu0,max(0.0_jprb, &
               &  mat_x_mat(ng,ng3d,nreg,gamma_z1(1:ng3d,nreg+1:2*nreg,1:nreg), &
               &            ref_dir(1:ng3d,:,:,jlev)) &
               &  + gamma_z1(1:ng3d,nreg+1:2*nreg,2*nreg+1:3*nreg)))

        end if ! we are treating 3d effects for some g points

        ! --- section 3.3b: g-points without 3d effects ----------

        ! compute reflectance, transmittance and associated terms for
        ! clear skies, using the meador-weaver formulas
        call calc_reflectance_transmittance_sw(ng, &
             &  mu0, od_region(1:ng,1), ssa_region(1:ng,1), &
             &  gamma1(1:ng,1), gamma2(1:ng,1), gamma3(1:ng,1), &
             &  ref_clear(1:ng,jlev), trans_clear(1:ng,jlev), &
             &  ref_dir_clear(1:ng,jlev), trans_dir_diff_clear(1:ng,jlev), &
             &  trans_dir_dir_clear(1:ng,jlev) )

        n_calls_meador_weaver = n_calls_meador_weaver + ng

        if (ng3d < ng) then
          ! some of the g points are to be treated using the
          ! conventional plane-parallel method.  first zero the
          ! relevant parts of the matrices
          trans_dir_dir (ng3d+1:ng,:,:,jlev) = 0.0_jprb
          reflectance   (ng3d+1:ng,:,:,jlev) = 0.0_jprb
          transmittance (ng3d+1:ng,:,:,jlev) = 0.0_jprb
          ref_dir       (ng3d+1:ng,:,:,jlev) = 0.0_jprb
          trans_dir_diff(ng3d+1:ng,:,:,jlev) = 0.0_jprb

          ! since there is no lateral transport, the clear-sky parts
          ! of the arrays can be copied from the clear-sky arrays
          trans_dir_dir (ng3d+1:ng,1,1,jlev) = trans_dir_dir_clear (ng3d+1:ng,jlev)
          reflectance   (ng3d+1:ng,1,1,jlev) = ref_clear           (ng3d+1:ng,jlev)
          transmittance (ng3d+1:ng,1,1,jlev) = trans_clear         (ng3d+1:ng,jlev)
          ref_dir       (ng3d+1:ng,1,1,jlev) = ref_dir_clear       (ng3d+1:ng,jlev)
          trans_dir_diff(ng3d+1:ng,1,1,jlev) = trans_dir_diff_clear(ng3d+1:ng,jlev)

          ! compute reflectance, transmittance and associated terms
          ! for each cloudy region, using the meador-weaver formulas
          do jreg = 2, nregactive
            call calc_reflectance_transmittance_sw(ng-ng3d, &
                 &  mu0, &
                 &  od_region(ng3d+1:ng,jreg), ssa_region(ng3d+1:ng,jreg), &
                 &  gamma1(ng3d+1:ng,jreg), gamma2(ng3d+1:ng,jreg), &
                 &  gamma3(ng3d+1:ng,jreg), &
                 &  reflectance(ng3d+1:ng,jreg,jreg,jlev), &
                 &  transmittance(ng3d+1:ng,jreg,jreg,jlev), &
                 &  ref_dir(ng3d+1:ng,jreg,jreg,jlev), &
                 &  trans_dir_diff(ng3d+1:ng,jreg,jreg,jlev), &
                 &  trans_dir_dir(ng3d+1:ng,jreg,jreg,jlev) )
          end do
          n_calls_meador_weaver &
               &  = n_calls_meador_weaver + (ng-ng3d)*(nregactive-1)
        end if

      end do ! loop over levels

      ! --------------------------------------------------------
      ! section 4: compute total albedos
      ! --------------------------------------------------------

      total_albedo(:,:,:,:)        = 0.0_jprb
      total_albedo_direct(:,:,:,:) = 0.0_jprb

      if (config%do_clear) then
        total_albedo_clear(:,:)        = 0.0_jprb
        total_albedo_clear_direct(:,:) = 0.0_jprb
      end if

      ! calculate the upwelling radiation scattered from the direct
      ! beam incident on the surface, and copy the surface albedo
      ! into total_albedo
      do jreg = 1,nreg
        do jg = 1,ng
          total_albedo(jg,jreg,jreg,nlev+1) = albedo_diffuse(jg,jcol)
          total_albedo_direct(jg,jreg,jreg,nlev+1) &
               &  = mu0 * albedo_direct(jg,jcol)
        end do
      end do

      if (config%do_clear) then
        ! surface albedo is the same
        total_albedo_clear(1:ng,nlev+1) = total_albedo(1:ng,1,1,nlev+1)
        total_albedo_clear_direct(1:ng,nlev+1) &
             &  = total_albedo_direct(1:ng,1,1,nlev+1)
      end if

      ! horizontal migration distances of reflected radiation at the
      ! surface are zero
      x_diffuse = 0.0_jprb
      x_direct  = 0.0_jprb

      ! work back up through the atmosphere computing the total albedo
      ! of the atmosphere below that point using the adding method
      do jlev = nlev,1,-1

        ! --------------------------------------------------------
        ! section 4.1: adding method
        ! --------------------------------------------------------

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
          total_albedo_clear_direct(:,jlev) = ref_dir_clear(:,jlev) &
               &  + (trans_dir_dir_clear(:,jlev) * total_albedo_clear_direct(:,jlev+1) &
               &    +trans_dir_diff_clear(:,jlev) * total_albedo_clear(:,jlev+1)) &
               &  * trans_clear(:,jlev) * inv_denom_scalar(:)
        end if

        if (is_clear_sky_layer(jlev)) then
          ! clear-sky layer: use scalar adding method
          inv_denom_scalar(:) = 1.0_jprb &
               &  / (1.0_jprb - total_albedo(:,1,1,jlev+1)*reflectance(:,1,1,jlev))
          total_albedo_below = 0.0_jprb
          total_albedo_below(:,1,1) = reflectance(:,1,1,jlev) &
               &  + transmittance(:,1,1,jlev)  * transmittance(:,1,1,jlev) &
               &  * total_albedo(:,1,1,jlev+1) * inv_denom_scalar(:)
          total_albedo_below_direct = 0.0_jprb
          total_albedo_below_direct(:,1,1) = ref_dir(:,1,1,jlev) &
               &  + (trans_dir_dir(:,1,1,jlev)*total_albedo_direct(:,1,1,jlev+1) &
               &    +trans_dir_diff(:,1,1,jlev)*total_albedo(:,1,1,jlev+1)) &
               &  * transmittance(:,1,1,jlev) * inv_denom_scalar(:)
        else 
          ! cloudy layer: use matrix adding method
          denominator = identity_minus_mat_x_mat(ng,ng,nreg, &
               &  total_albedo(:,:,:,jlev+1), reflectance(:,:,:,jlev))
          total_albedo_below = reflectance(:,:,:,jlev) &
               &  + mat_x_mat(ng,ng,nreg,transmittance(:,:,:,jlev), &
               &  solve_mat(ng,ng,nreg,denominator, &
               &  mat_x_mat(ng,ng,nreg,total_albedo(:,:,:,jlev+1), &
               &  transmittance(:,:,:,jlev))))
          total_albedo_below_direct = ref_dir(:,:,:,jlev) &
               &  + mat_x_mat(ng,ng,nreg,transmittance(:,:,:,jlev), &
               &  solve_mat(ng,ng,nreg,denominator, &
               &    mat_x_mat(ng,ng,nreg,total_albedo_direct(:,:,:,jlev+1), &
               &                       trans_dir_dir(:,:,:,jlev)) &
               &   +mat_x_mat(ng,ng,nreg,total_albedo(:,:,:,jlev+1), &
               &                       trans_dir_diff(:,:,:,jlev))))
        end if

        ! --------------------------------------------------------
        ! section 4.2: overlap and entrapment
        ! --------------------------------------------------------


        if ((config%i_3d_sw_entrapment == ientrapmentexplicitnonfractal &
             &  .or. config%i_3d_sw_entrapment == ientrapmentexplicit) &
             &  .and. jlev >= i_cloud_top) then




          !  "explicit entrapment": we have the horizontal migration
          !  distances just above the base of the layer, and need to
          !  step them to just below the top of the same layer
          call step_migrations(ng, nreg, cloud%fraction(jcol,jlev), &
               & layer_depth(jlev), tan_diffuse_angle_3d, tan_sza, &
               &  reflectance(:,:,:,jlev), transmittance(:,:,:,jlev), &
               &  ref_dir(:,:,:,jlev), trans_dir_dir(:,:,:,jlev), &
               &  trans_dir_diff(:,:,:,jlev), total_albedo(:,:,:,jlev+1), &
               &  total_albedo_direct(:,:,:,jlev+1), &
               &  x_diffuse, x_direct)

! # 966 "radiation/radiation_spartacus_sw.f90"

        end if

        ! account for cloud overlap when converting albedo and source
        ! below a layer interface to the equivalent values just above
        if (is_clear_sky_layer(jlev) .and. is_clear_sky_layer(jlev-1)) then
          ! if both layers are cloud free, this is trivial...
          total_albedo(:,:,:,jlev) = 0.0_jprb
          total_albedo(:,1,1,jlev) = total_albedo_below(:,1,1)
          total_albedo_direct(:,:,:,jlev) = 0.0_jprb
          total_albedo_direct(:,1,1,jlev) = total_albedo_below_direct(:,1,1)

        else if (config%i_3d_sw_entrapment == ientrapmentmaximum &
             &  .or. is_clear_sky_layer(jlev-1)) then
          ! "maximum entrapment": use the overlap matrices u_matrix and v_matrix
          ! (this is the original spartacus method)
          total_albedo(:,:,:,jlev) = singlemat_x_mat(ng,ng,nreg,&
               &  u_matrix(:,:,jlev,jcol), &
               &  mat_x_singlemat(ng,ng,nreg,total_albedo_below,&
               &  v_matrix(:,:,jlev,jcol)))
          total_albedo_direct(:,:,:,jlev) = singlemat_x_mat(ng,ng,nreg,&
               &  u_matrix(:,:,jlev,jcol), &
               &  mat_x_singlemat(ng,ng,nreg,total_albedo_below_direct,&
               &  v_matrix(:,:,jlev,jcol)))

        else if (config%i_3d_sw_entrapment == ientrapmentzero) then
          ! "zero entrapment": even radiation transported
          ! laterally between regions in the layers below is
          ! reflected back up into the same region. first diffuse
          ! radiation:
          total_albedo(:,:,:,jlev) = 0.0_jprb
          do jreg = 1,nreg    ! target layer (jlev-1)
            do jreg2 = 1,nreg ! current layer (jlev)
              total_albedo(:,jreg,jreg,jlev) = total_albedo(:,jreg,jreg,jlev) &
                   &  + sum(total_albedo_below(:,:,jreg2),2) &
                   &  * v_matrix(jreg2,jreg,jlev,jcol)
            end do
          end do
          ! ...then direct radiation:
          total_albedo_direct(:,:,:,jlev) = 0.0_jprb
          do jreg = 1,nreg    ! target layer (jlev-1)
            do jreg2 = 1,nreg ! current layer (jlev)
              total_albedo_direct(:,jreg,jreg,jlev) = total_albedo_direct(:,jreg,jreg,jlev) &
                   &  + sum(total_albedo_below_direct(:,:,jreg2),2) &
                   &  * v_matrix(jreg2,jreg,jlev,jcol)
            end do
          end do

        else
          ! controlled entrapment

! # 1027 "radiation/radiation_spartacus_sw.f90"
          ! add the contribution from off-diagonal elements of the
          ! albedo matrix in the lower layer, i.e. radiation that
          ! flows between regions...

          ! first diffuse radiation:
          albedo_part = total_albedo_below
          do jreg = 1,nreg
            albedo_part(:,jreg,jreg) = 0.0_jprb
          end do
          total_albedo(:,:,:,jlev) = singlemat_x_mat(ng,ng,nreg,&
               &  u_matrix(:,:,jlev,jcol), &
               &  mat_x_singlemat(ng,ng,nreg,albedo_part,&
               &  v_matrix(:,:,jlev,jcol)))
          ! ...then direct radiation:
          albedo_part = total_albedo_below_direct
          do jreg = 1,nreg
            albedo_part(:,jreg,jreg) = 0.0_jprb
          end do
          total_albedo_direct(:,:,:,jlev) = singlemat_x_mat(ng,ng,nreg,&
               &  u_matrix(:,:,jlev,jcol), &
               &  mat_x_singlemat(ng,ng,nreg,albedo_part,&
               &  v_matrix(:,:,jlev,jcol)))




          
          ! now the contribution from the diagonals of the albedo
          ! matrix in the lower layer
          if (config%i_3d_sw_entrapment == ientrapmentedgeonly &
               &  .or. (.not. config%do_3d_effects)) then
            ! "edge-only entrapment": the operation we perform is
            ! essentially diag(total_albedo) += matmul(transpose(v_matrix),
            ! diag(total_albedo_below)).
            do jreg = 1,nreg
              do jreg2 = 1,nreg
                total_albedo(:,jreg,jreg,jlev) &
                     &  = total_albedo(:,jreg,jreg,jlev) &
                     &  + total_albedo_below(:,jreg2,jreg2) &
                     &  * v_matrix(jreg2,jreg,jlev,jcol)
                total_albedo_direct(:,jreg,jreg,jlev) &
                     &  = total_albedo_direct(:,jreg,jreg,jlev) &
                     &  + total_albedo_below_direct(:,jreg2,jreg2) &
                     &  * v_matrix(jreg2,jreg,jlev,jcol)
              end do
            end do

          else
            ! "explicit entrapment"

            do jreg2 = 1,nreg
              ! loop through each region in the lower layer. for one
              ! of the regions in the lower layer, we are imagining it
              ! to be divided into "nreg" subregions that map onto the
              ! regions in the upper layer. the rate of exchange
              ! between these subregions is computed via a coupled
              ! differential equation written in terms of a singular
              ! exchange matrix (there are only terms for the exchange
              ! between subregions, but no propagation effects since
              ! we already know the albedo of this region). this is
              ! solved using the matrix-exponential method.

              ! use the following array for the transfer of either
              ! diffuse or direct radiation (despite the name), per
              ! unit horizontal distance travelled
              transfer_rate_diffuse = 0.0_jprb

              ! as we need to reference the layer above the interface,
              ! don't do the following on the highest layer
              if (jlev > 1) then

                ! given a horizontal migration distance, there is
                ! still uncertainty about how much entrapment occurs
                ! associated with how one assumes cloud boundaries
                ! line up in adjacent layers. "overhang_factor"
                ! can be varied between 0.0 (the boundaries line up to
                ! the greatest extent possible given the overlap
                ! parameter) and 1.0 (the boundaries line up to the
                ! minimum extent possible); here this is used to
                ! produce a scaling factor for the transfer rate.
                transfer_scaling = 1.0_jprb - (1.0_jprb - config%overhang_factor) & 
                     &  * cloud%overlap_param(jcol,jlev-1) &
                     &  * min(region_fracs(jreg2,jlev,jcol),region_fracs(jreg2,jlev-1,jcol)) &
                     &  / max(config%cloud_fraction_threshold, region_fracs(jreg2,jlev,jcol))

                do jreg = 1, nreg-1
                  ! compute lateral transfer rates from region jreg to
                  ! jreg+1 as before, but without length scale which
                  ! is wavelength dependent.

                  ! recall that overlap indexing is
                  ! u_matrix(upper_region, lower_region, level,
                  ! column).
                  transfer_rate_diffuse(jreg,jreg+1) = transfer_scaling &
                       &  * edge_length(jreg,jlev-1) / max(u_matrix(jreg,jreg2,jlev,jcol),1.0e-5_jprb)
                  ! compute transfer rates from region jreg+1 to jreg
                  transfer_rate_diffuse(jreg+1,jreg) = transfer_scaling &
                       &  * edge_length(jreg,jlev-1) / max(u_matrix(jreg+1,jreg2,jlev,jcol),1.0e-5_jprb)
                end do
              
                ! compute transfer rates directly between regions 1
                ! and 3 (not used below)
                if (edge_length(3,jlev) > 0.0_jprb) then
                  transfer_rate_diffuse(1,3) = transfer_scaling &
                       &  * edge_length(3,jlev-1) / max(u_matrix(1,jreg2,jlev,jcol),1.0e-5_jprb)
                  transfer_rate_diffuse(3,1) = transfer_scaling &
                       &  * edge_length(3,jlev-1) / max(u_matrix(3,jreg2,jlev,jcol),1.0e-5_jprb)
                end if
              end if

              ! compute matrix of exchange coefficients
              entrapment = 0.0_jprb
              inv_effective_size = min(cloud%inv_cloud_effective_size(jcol,jlev-1), &
                   &                   1.0_jprb/config%min_cloud_effective_size)
              do jreg = 1,nreg-1
                ! diffuse transport down and up with one random
                ! scattering event
                if (config%i_3d_sw_entrapment == ientrapmentexplicit) then
                  fractal_factor = 1.0_jprb / sqrt(max(1.0_jprb, 2.5_jprb*x_diffuse(:,jreg2) &
                       &                                         * inv_effective_size))
                  entrapment(:,jreg+1,jreg) = entrapment(:,jreg+1,jreg) &
                       &  + transfer_rate_diffuse(jreg,jreg+1)*x_diffuse(:,jreg2) &
                       &  * fractal_factor
                  entrapment(:,jreg,jreg+1) = entrapment(:,jreg,jreg+1) &
                       &  + transfer_rate_diffuse(jreg+1,jreg)*x_diffuse(:,jreg2) &
                       &  * fractal_factor
                else
                  entrapment(:,jreg+1,jreg) = entrapment(:,jreg+1,jreg) &
                       &  + transfer_rate_diffuse(jreg,jreg+1)*x_diffuse(:,jreg2)
                  entrapment(:,jreg,jreg+1) = entrapment(:,jreg,jreg+1) &
                       &  + transfer_rate_diffuse(jreg+1,jreg)*x_diffuse(:,jreg2)                  
                end if
                entrapment(:,jreg,jreg) = entrapment(:,jreg,jreg) &
                     &  - entrapment(:,jreg+1,jreg)
                entrapment(:,jreg+1,jreg+1) = entrapment(:,jreg+1,jreg+1) &
                     &  - entrapment(:,jreg,jreg+1)
              end do

              ! if rate of exchange is excessive the expm can throw a
              ! floating point exception, even if it tends towards a
              ! trival limit, so we cap the maximum input to expm by
              ! scaling down if necessary
              do jg = 1,ng
                max_entr = -min(entrapment(jg,1,1),entrapment(jg,2,2))
                if (max_entr > config%max_cloud_od) then
                  ! scale down all inputs for this g point
                  entrapment(jg,:,:) = entrapment(jg,:,:) * (config%max_cloud_od/max_entr)
                end if
              end do

              ! since the matrix to be exponentiated has a simple
              ! structure we may use a faster method described in the
              ! appendix of hogan et al. (gmd 2018)


              if (nreg == 2) then
                call fast_expm_exchange(ng, ng, entrapment(:,2,1), entrapment(:,1,2), &
                     &                  albedo_part)
              else
                call fast_expm_exchange(ng, ng, entrapment(:,2,1), entrapment(:,1,2), &
                     &                          entrapment(:,3,2), entrapment(:,2,3), &
                     &                  albedo_part)
              end if








              ! scale to get the contribution to the diffuse albedo
              do jreg3 = 1,nreg
                do jreg = 1,nreg
                  albedo_part(:,jreg3,jreg) = albedo_part(:,jreg3,jreg) &
                       &  * v_matrix(jreg2,jreg,jlev,jcol) * total_albedo_below(:,jreg2,jreg2)
                end do
              end do
! # 1232 "radiation/radiation_spartacus_sw.f90"

              ! increment diffuse albedo
              total_albedo(:,:,:,jlev) = total_albedo(:,:,:,jlev) + albedo_part

              ! now do the same for the direct albedo
              entrapment = 0.0_jprb
              do jreg = 1,nreg-1
                ! direct transport down and diffuse up with one random
                ! scattering event
                if (config%i_3d_sw_entrapment == ientrapmentexplicit) then
                  fractal_factor = 1.0_jprb / sqrt(max(1.0_jprb, 2.5_jprb*x_direct(:,jreg2) &
                       &                                         * inv_effective_size))
                  entrapment(:,jreg+1,jreg) = entrapment(:,jreg+1,jreg) &
                       &  + transfer_rate_diffuse(jreg,jreg+1)*x_direct(:,jreg2) &
                       &  * fractal_factor
                  entrapment(:,jreg,jreg+1) = entrapment(:,jreg,jreg+1) &
                       &  + transfer_rate_diffuse(jreg+1,jreg)*x_direct(:,jreg2) &
                       &  * fractal_factor
                else
                  entrapment(:,jreg+1,jreg) = entrapment(:,jreg+1,jreg) &
                       &  + transfer_rate_diffuse(jreg,jreg+1)*x_direct(:,jreg2)
                  entrapment(:,jreg,jreg+1) = entrapment(:,jreg,jreg+1) &
                       &  + transfer_rate_diffuse(jreg+1,jreg)*x_direct(:,jreg2)
                end if
                entrapment(:,jreg,jreg) = entrapment(:,jreg,jreg) &
                     &  - entrapment(:,jreg+1,jreg)
                entrapment(:,jreg+1,jreg+1) = entrapment(:,jreg+1,jreg+1) &
                     &  - entrapment(:,jreg,jreg+1)
              end do

              ! if rate of exchange is excessive the expm can throw a
              ! floating point exception, even if it tends towards a
              ! trival limit, so we cap the maximum input to expm by
              ! scaling down if necessary
              do jg = 1,ng
                max_entr = -min(entrapment(jg,1,1),entrapment(jg,2,2))
                if (max_entr > config%max_cloud_od) then
                  ! scale down all inputs for this g point
                  entrapment(jg,:,:) = entrapment(jg,:,:) * (config%max_cloud_od/max_entr)
                end if
              end do



              if (nreg == 2) then
                call fast_expm_exchange(ng, ng, entrapment(:,2,1), entrapment(:,1,2), &
                     &                  albedo_part)
              else
                call fast_expm_exchange(ng, ng, entrapment(:,2,1), entrapment(:,1,2), &
                     &                          entrapment(:,3,2), entrapment(:,2,3), &
                     &                  albedo_part)
              end if







              do jreg3 = 1,nreg
                do jreg = 1,nreg
                  albedo_part(:,jreg3,jreg) = albedo_part(:,jreg3,jreg) &
                       &  * v_matrix(jreg2,jreg,jlev,jcol) * total_albedo_below_direct(:,jreg2,jreg2)
                end do
              end do
! # 1321 "radiation/radiation_spartacus_sw.f90"
              ! increment direct albedo
              total_albedo_direct(:,:,:,jlev) = total_albedo_direct(:,:,:,jlev) + albedo_part

            end do

          end if
        end if

        if ((config%i_3d_sw_entrapment == ientrapmentexplicitnonfractal &
             &  .or. config%i_3d_sw_entrapment == ientrapmentexplicit) &
             &  .and. .not. (is_clear_sky_layer(jlev) .and. is_clear_sky_layer(jlev-1))) then
          ! horizontal migration distances are averaged when
          ! applying overlap rules, so equation is
          ! x_above=matmul(transpose(v_matrix),x_below)
          
          ! we do this into temporary arrays...
          x_direct_above = 0.0_jprb
          x_diffuse_above = 0.0_jprb
          
          nregactive = nreg
          if (is_clear_sky_layer(jlev)) then
            nregactive = 1
          end if

          do jreg = 1,nreg          ! target layer (jlev-1)
            do jreg2 = 1,nregactive ! current layer (jlev)
              x_direct_above(:,jreg) = x_direct_above(:,jreg) &
                   &  + x_direct(:,jreg2) * v_matrix(jreg2,jreg,jlev,jcol)
              x_diffuse_above(:,jreg) = x_diffuse_above(:,jreg) &
                   &  + x_diffuse(:,jreg2) * v_matrix(jreg2,jreg,jlev,jcol)
            end do
          end do
          
          !... then copy out of the temporary arrays
          x_direct = x_direct_above
          x_diffuse = x_diffuse_above
        end if

      end do ! reverse loop over levels

      ! --------------------------------------------------------
      ! section 5: compute fluxes
      ! --------------------------------------------------------

      ! top-of-atmosphere fluxes into the regions of the top-most
      ! layer, zero since we assume no diffuse downwelling
      flux_dn_below = 0.0_jprb
      ! direct downwelling flux (into a plane perpendicular to the
      ! sun) entering the top of each region in the top-most layer
      do jreg = 1,nreg
        direct_dn_below(:,jreg) = incoming_sw(:,jcol)*region_fracs(jreg,1,jcol)
      end do
      ! we're using flux_up_above as a container; actually its
      ! interpretation at top of atmosphere here is just 'below' the
      ! toa interface, so using the regions of the first model layer
      flux_up_above = mat_x_vec(ng,ng,nreg,total_albedo_direct(:,:,:,1),direct_dn_below)

      if (config%do_clear) then
        flux_dn_clear = 0.0_jprb
        direct_dn_clear(:) = incoming_sw(:,jcol)
        flux_up_clear = direct_dn_clear*total_albedo_clear_direct(:,1)
      end if

      ! store the toa broadband fluxes
      flux%sw_up(jcol,1) = sum(sum(flux_up_above,1))
      flux%sw_dn(jcol,1) = mu0 * sum(direct_dn_clear(:))
      if (allocated(flux%sw_dn_direct)) then
        flux%sw_dn_direct(jcol,1) = flux%sw_dn(jcol,1)
      end if
      if (config%do_clear) then
        flux%sw_up_clear(jcol,1) = sum(flux_up_clear)
        flux%sw_dn_clear(jcol,1) = flux%sw_dn(jcol,1)
        if (allocated(flux%sw_dn_direct_clear)) then
          flux%sw_dn_direct_clear(jcol,1) = flux%sw_dn_clear(jcol,1)
        end if
      end if

      ! save the spectral fluxes if required
      if (config%do_save_spectral_flux) then
        call indexed_sum(sum(flux_up_above(:,:),2), &
             &           config%i_spec_from_reordered_g_sw, &
             &           flux%sw_up_band(:,jcol,1))
        call indexed_sum(sum(direct_dn_below(:,:),2), &
             &           config%i_spec_from_reordered_g_sw, &
             &           flux%sw_dn_band(:,jcol,1))
        flux%sw_dn_band(:,jcol,1) = mu0 * flux%sw_dn_band(:,jcol,1)
        if (allocated(flux%sw_dn_direct_band)) then
          flux%sw_dn_direct_band(:,jcol,1) = flux%sw_dn_band(:,jcol,1)
        end if
        if (config%do_clear) then
          flux%sw_dn_clear_band(:,jcol,1) = flux%sw_dn_band(:,jcol,1)
          call indexed_sum(flux_up_clear, &
               &           config%i_spec_from_reordered_g_sw, &
               &           flux%sw_up_clear_band(:,jcol,1))
          if (allocated(flux%sw_dn_direct_clear_band)) then
            flux%sw_dn_direct_clear_band(:,jcol,1) &
                 &   = flux%sw_dn_clear_band(:,jcol,1)
          end if
        end if
      end if

      ! final loop back down through the atmosphere to compute fluxes
      do jlev = 1,nlev

! # 1437 "radiation/radiation_spartacus_sw.f90"

        ! compute the solar downwelling "source" at the base of the
        ! layer due to scattering of the direct beam within it
        if (config%do_clear) then
          source_dn_clear = trans_dir_diff_clear(:,jlev)*direct_dn_clear
        end if
        source_dn(:,:) = mat_x_vec(ng,ng,nreg,trans_dir_diff(:,:,:,jlev),direct_dn_below, &
             &  is_clear_sky_layer(jlev))

        ! compute direct downwelling flux in each region at base of
        ! current layer
        if (config%do_clear) then
          direct_dn_clear = trans_dir_dir_clear(:,jlev)*direct_dn_clear
        end if
        direct_dn_above = mat_x_vec(ng,ng,nreg,trans_dir_dir(:,:,:,jlev),direct_dn_below, &
             &  is_clear_sky_layer(jlev))

        ! integrate downwelling direct flux across spectrum and
        ! regions, and store (the diffuse part will be added later)
        flux%sw_dn(jcol,jlev+1) = mu0 * sum(sum(direct_dn_above,1))
        if (allocated(flux%sw_dn_direct)) then
          flux%sw_dn_direct(jcol,jlev+1) = flux%sw_dn(jcol,jlev+1)
        end if
        if (config%do_clear) then
          flux%sw_dn_clear(jcol,jlev+1) = mu0 * sum(direct_dn_clear)
          if (allocated(flux%sw_dn_direct_clear)) then
            flux%sw_dn_direct_clear(jcol,jlev+1) &
                 &  = flux%sw_dn_clear(jcol,jlev+1)
          end if
        end if

        if (config%do_save_spectral_flux) then
          call indexed_sum(sum(direct_dn_above,2), &
               &           config%i_spec_from_reordered_g_sw, &
               &           flux%sw_dn_band(:,jcol,jlev+1))
          flux%sw_dn_band(:,jcol,jlev+1) = mu0 * flux%sw_dn_band(:,jcol,jlev+1)

          if (allocated(flux%sw_dn_direct_band)) then
            flux%sw_dn_direct_band(:,jcol,jlev+1) &
                 &   = flux%sw_dn_band(:,jcol,jlev+1)
          end if
          if (config%do_clear) then
            call indexed_sum(direct_dn_clear, &
                 &           config%i_spec_from_reordered_g_sw, &
                 &           flux%sw_dn_clear_band(:,jcol,jlev+1))
            flux%sw_dn_clear_band(:,jcol,jlev+1) = mu0 &
                 &   * flux%sw_dn_clear_band(:,jcol,jlev+1)
            if (allocated(flux%sw_dn_direct_clear_band)) then
              flux%sw_dn_direct_clear_band(:,jcol,jlev+1) &
                   &  = flux%sw_dn_clear_band(:,jcol,jlev+1)
            end if
          end if
        end if

        if (config%do_clear) then
          ! scalar operations for clear-sky fluxes
          flux_dn_clear(:) = (trans_clear(:,jlev)*flux_dn_clear(:) &
               &  + ref_clear(:,jlev)*total_albedo_clear_direct(:,jlev+1)*direct_dn_clear &
               &  + source_dn_clear) &
               &  / (1.0_jprb - ref_clear(:,jlev)*total_albedo_clear(:,jlev+1))
          flux_up_clear(:) = total_albedo_clear_direct(:,jlev+1)*direct_dn_clear &
               &  + total_albedo_clear(:,jlev+1)*flux_dn_clear
        end if

        if (is_clear_sky_layer(jlev)) then
          ! scalar operations for clear-sky layer
          flux_dn_above(:,1) = (transmittance(:,1,1,jlev)*flux_dn_below(:,1) &
               &  + reflectance(:,1,1,jlev)*total_albedo_direct(:,1,1,jlev+1)*direct_dn_above(:,1) &
               &  + source_dn(:,1)) &
               &  / (1.0_jprb - reflectance(:,1,1,jlev)*total_albedo(:,1,1,jlev+1))
          flux_dn_above(:,2:nreg) = 0.0_jprb
          flux_up_above(:,1) = total_albedo_direct(:,1,1,jlev+1)*direct_dn_above(:,1) &
               &  + total_albedo(:,1,1,jlev+1)*flux_dn_above(:,1)
          flux_up_above(:,2:nreg) = 0.0_jprb
        else
          ! matrix operations for cloudy layer
          denominator = identity_minus_mat_x_mat(ng,ng,nreg,reflectance(:,:,:,jlev), &
               &  total_albedo(:,:,:,jlev+1))
          total_source = mat_x_vec(ng,ng,nreg,total_albedo_direct(:,:,:,jlev+1),direct_dn_above)

          flux_dn_above = solve_vec(ng,ng,nreg,denominator, &
               &  mat_x_vec(ng,ng,nreg,transmittance(:,:,:,jlev),flux_dn_below) &
               &  + mat_x_vec(ng,ng,nreg,reflectance(:,:,:,jlev), total_source(:,:)) &
               &  + source_dn(:,:))
          flux_up_above = mat_x_vec(ng,ng,nreg,total_albedo(:,:,:,jlev+1), &
               &  flux_dn_above) + total_source(:,:)
        end if

        ! account for overlap rules in translating fluxes just above
        ! a layer interface to the values just below
        if (is_clear_sky_layer(jlev) .and. is_clear_sky_layer(jlev+1)) then
          ! regions in current layer map directly on to regions in
          ! layer below
          flux_dn_below = flux_dn_above
          direct_dn_below = direct_dn_above
        else
          ! apply downward overlap matrix to compute direct
          ! downwelling flux entering the top of each region in the
          ! layer below
          flux_dn_below = singlemat_x_vec(ng,ng,nreg,v_matrix(:,:,jlev+1,jcol), &
               &  flux_dn_above)
          direct_dn_below = singlemat_x_vec(ng,ng,nreg,v_matrix(:,:,jlev+1,jcol), &
               &  direct_dn_above)
        end if

        ! store the broadband fluxes
        flux%sw_up(jcol,jlev+1) = sum(sum(flux_up_above,1))
        flux%sw_dn(jcol,jlev+1) &
             &  = flux%sw_dn(jcol,jlev+1) + sum(sum(flux_dn_above,1))
        if (config%do_clear) then
          flux%sw_up_clear(jcol,jlev+1) = sum(flux_up_clear)
          flux%sw_dn_clear(jcol,jlev+1) &
               &  = flux%sw_dn_clear(jcol,jlev+1) + sum(flux_dn_clear)
        end if

        ! save the spectral fluxes if required
        if (config%do_save_spectral_flux) then
          call indexed_sum(sum(flux_up_above,2), &
               &           config%i_spec_from_reordered_g_sw, &
               &           flux%sw_up_band(:,jcol,jlev+1))
          call add_indexed_sum(sum(flux_dn_above,2), &
               &           config%i_spec_from_reordered_g_sw, &
               &           flux%sw_dn_band(:,jcol,jlev+1))
          if (config%do_clear) then
            call indexed_sum(flux_up_clear, &
                 &           config%i_spec_from_reordered_g_sw, &
                 &           flux%sw_up_clear_band(:,jcol,jlev+1))
            call add_indexed_sum(flux_dn_clear, &
                 &           config%i_spec_from_reordered_g_sw, &
                 &           flux%sw_dn_clear_band(:,jcol,jlev+1))
          end if
        end if

      end do ! final loop over levels

      ! store surface spectral fluxes, if required (after the end of
      ! the final loop over levels, the current values of these arrays
      ! will be the surface values)
      flux%sw_dn_diffuse_surf_g(:,jcol) = sum(flux_dn_above,2)
      flux%sw_dn_direct_surf_g(:,jcol)  = mu0 * sum(direct_dn_above,2)
      if (config%do_clear) then
        flux%sw_dn_diffuse_surf_clear_g(:,jcol) = flux_dn_clear
        flux%sw_dn_direct_surf_clear_g(:,jcol)  = mu0 * direct_dn_clear
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

    if (lhook) call dr_hook('radiation_spartacus_sw:solver_spartacus_sw',1,hook_handle)

  end subroutine solver_spartacus_sw


  ! step the horizontal migration distances from the base of a layer
  ! to the top, accounting for the extra distance travelled within the
  ! layer
  subroutine step_migrations(ng, nreg, cloud_frac, &
       &  layer_depth, tan_diffuse_angle_3d, tan_sza, &
       &  reflectance, transmittance, ref_dir, trans_dir_dir, &
       &  trans_dir_diff, total_albedo_diff, total_albedo_dir, &
       &  x_diffuse, x_direct)
    
    use parkind1, only : jprb

    implicit none

    ! inputs

    ! number of g points and regions
    integer, intent(in) :: ng, nreg
    ! cloud fraction
    real(jprb), intent(in) :: cloud_frac
    ! layer depth (m), tangent of diffuse zenith angle and tangent of
    ! solar zenith angle
    real(jprb), intent(in) :: layer_depth, tan_diffuse_angle_3d, tan_sza
    ! reflectance and transmittance to diffuse downwelling radiation
    real(jprb), intent(in), dimension(ng, nreg, nreg) :: reflectance, transmittance
    ! reflectance and transmittance to direct downwelling radiation
    real(jprb), intent(in), dimension(ng, nreg, nreg) :: ref_dir, trans_dir_dir
    ! transmittance involving direct entering a layer from the top and
    ! diffuse leaving from the bottom
    real(jprb), intent(in), dimension(ng, nreg, nreg) :: trans_dir_diff

    ! total albedo of direct and diffuse radiation of the atmosphere
    ! below the layer in question
    real(jprb), intent(in), dimension(ng, nreg, nreg) &
         &  :: total_albedo_diff, total_albedo_dir

    ! inputs/outputs

    ! horizontal migration distance (m) of reflected light
    real(jprb), intent(inout), dimension(ng, nreg) :: x_diffuse, x_direct

    ! local variables

    ! top albedo, i.e. the albedo of the top of the layer assuming no
    ! lateral transport
    real(jprb), dimension(ng) :: top_albedo

    ! multiple-scattering amplitude enhancement
    real(jprb), dimension(ng) :: ms_enhancement

    ! multiple-scattering distance enhancement
    real(jprb), dimension(ng) :: x_enhancement


    real(jprb) :: x_layer_diffuse, x_layer_direct
    integer :: jreg, istartreg, iendreg

    istartreg = 1
    iendreg   = nreg

    if (cloud_frac <= 0.0_jprb) then
      ! clear-sky layer: don't waste time on cloudy regions
      iendreg = 1
    else if (cloud_frac >= 1.0_jprb) then
      ! overcast layer: don't waste time on clear region
      istartreg = 2
    end if

    ! this is the mean horizontal distance travelled by diffuse
    ! radiation that travels from the top of a layer to the centre and
    ! is then scattered back up and out
    x_layer_diffuse = layer_depth * tan_diffuse_angle_3d/sqrt(2.0_jprb) 

    ! this is the mean horizontal distance travelled by direct
    ! radiation that travels from the top of a layer to the centre and
    ! is then scattered back up and out
    x_layer_direct  = layer_depth * sqrt(tan_sza*tan_sza &
         &                             + tan_diffuse_angle_3d*tan_diffuse_angle_3d) * 0.5_jprb

    do jreg = istartreg,iendreg
      ! geometric series enhancement due to multiple scattering: the
      ! amplitude enhancement is equal to the limit of
      ! t*[1+ra+(ra)^2+(ra)^3+...]
      ms_enhancement = transmittance(:,jreg,jreg) &
           &  / (1.0_jprb - reflectance(:,jreg,jreg)*total_albedo_diff(:,jreg,jreg))
      ! ...and the distance enhancement is approximately equal to the
      ! limit of t*[1+sqrt(2)*ra+sqrt(3)*(ra)^2+sqrt(4)*(ra)^3+...]
      x_enhancement = (1.0_jprb - reflectance(:,jreg,jreg)*total_albedo_diff(:,jreg,jreg))**(-1.5_jprb)

      ! horizontal migration of direct downwelling radiation
      top_albedo = max(1.0e-8_jprb, ref_dir(:,jreg,jreg) + ms_enhancement &
           &  * (trans_dir_diff(:,jreg,jreg)*total_albedo_diff(:,jreg,jreg) &
           &     +trans_dir_dir(:,jreg,jreg)*total_albedo_dir(:,jreg,jreg)))
      ! the following is approximate and has been found to
      ! occasionally go negative
      x_direct(:,jreg) = max(0.0_jprb, x_layer_direct &
           &  + ((trans_dir_diff(:,jreg,jreg)*total_albedo_diff(:,jreg,jreg)*x_enhancement &
           &      +trans_dir_dir(:,jreg,jreg)*total_albedo_dir(:,jreg,jreg)*(x_enhancement-1.0_jprb)) &
           &     *(x_diffuse(:,jreg)+x_layer_diffuse) &
           &    +trans_dir_dir(:,jreg,jreg)*total_albedo_dir(:,jreg,jreg) &
           &     *(x_direct(:,jreg)+x_layer_direct)) &
           &    * transmittance(:,jreg,jreg) / top_albedo)

      ! horizontal migration of diffuse downwelling radiation
      top_albedo = max(1.0e-8_jprb, reflectance(:,jreg,jreg) &
           &  + ms_enhancement*transmittance(:,jreg,jreg)*total_albedo_diff(:,jreg,jreg))
      x_diffuse(:,jreg) = x_layer_diffuse + x_enhancement*total_albedo_diff(:,jreg,jreg) &
           &  *(transmittance(:,jreg,jreg)*transmittance(:,jreg,jreg)) &
           &  * (x_diffuse(:,jreg) + x_layer_diffuse) / top_albedo

    end do
    if (iendreg < nreg) then
      x_diffuse(:,iendreg+1:nreg)      = 0.0_jprb
      x_direct(:,iendreg+1:nreg)       = 0.0_jprb
    else if (istartreg == 2) then
      x_diffuse(:,1)      = 0.0_jprb
      x_direct(:,1)       = 0.0_jprb
    end if

  end subroutine step_migrations

end module radiation_spartacus_sw
! #define __atomic_acquire 2
! #define __char_bit__ 8
! #define __float_word_order__ __order_little_endian__
! #define __order_little_endian__ 1234
! #define __order_pdp_endian__ 3412
! #define __gfc_real_10__ 1
! #define __finite_math_only__ 0
! #define __gnuc_patchlevel__ 0
! #define __gfc_int_2__ 1
! #define use_fast_expm_exchange 1
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

