! # 1 "radiation/radiation_cloud_optics.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_cloud_optics.f90"
! this file has been modified for the use in icon

! radiation_cloud_optics.f90 - computing cloud optical properties
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
!   2017-07-22  r. hogan  added yi et al. ice optics model

module radiation_cloud_optics

  implicit none

  public

contains

  ! provides elemental function "delta_eddington_scat_od"

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
! # 30 "radiation/radiation_cloud_optics.f90" 2

  !---------------------------------------------------------------------
  ! load cloud scattering data; this subroutine delegates to one
  ! in radiation_cloud_optics_data.f90, but checks the size of
  ! what is returned
  subroutine setup_cloud_optics(config)

    use parkind1,         only : jprb
    use ecradhook,          only : lhook, dr_hook, jphook

    use radiation_io,     only : nulerr, radiation_abort
    use radiation_config, only : config_type, iicemodelfu, iicemodelbaran, &
         &                       iicemodelbaran2016, iicemodelbaran2017, &
         &                       iicemodelyi, &
         &                       iliquidmodelsocrates, iliquidmodelslingo
    use radiation_ice_optics_fu, only    : niceopticscoeffsfusw, &
         &                                 niceopticscoeffsfulw
    use radiation_ice_optics_baran, only : niceopticscoeffsbaran, &
         &                                 niceopticscoeffsbaran2016
    use radiation_ice_optics_baran2017, only : niceopticscoeffsbaran2017, &
         &                                 niceopticsgeneralcoeffsbaran2017
    use radiation_ice_optics_yi, only    : niceopticscoeffsyisw, &
         &                                 niceopticscoeffsyilw
    use radiation_liquid_optics_socrates, only : nliqopticscoeffssocrates
    use radiation_liquid_optics_slingo, only : nliqopticscoeffsslingosw, &
         &                                     nliqopticscoeffslindnerlilw

    type(config_type), intent(inout) :: config

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloud_optics:setup_cloud_optics',0,hook_handle)

    call config%cloud_optics%setup(trim(config%liq_optics_file_name), &
         &     trim(config%ice_optics_file_name), iverbose=config%iverbosesetup)

    ! check liquid coefficients
    if (size(config%cloud_optics%liq_coeff_lw, 1) /= config%n_bands_lw) then
      write(nulerr,'(a,i0,a,i0,a)') &
           &  '*** error: number of longwave bands for droplets (', &
           &  size(config%cloud_optics%liq_coeff_lw, 1), &
           &  ') does not match number for gases (', config%n_bands_lw, ')'
      call radiation_abort()
    end if
    if (size(config%cloud_optics%liq_coeff_sw, 1) /= config%n_bands_sw) then
      write(nulerr,'(a,i0,a,i0,a)') &
           &  '*** error: number of shortwave bands for droplets (', &
           &  size(config%cloud_optics%liq_coeff_sw, 1), &
           &  ') does not match number for gases (', config%n_bands_sw, ')'
      call radiation_abort()
    end if

    if (config%i_liq_model == iliquidmodelsocrates) then
      if (size(config%cloud_optics%liq_coeff_lw, 2) /= nliqopticscoeffssocrates) then
        write(nulerr,'(a,i0,a,i0,a,i0,a)') &
             &  '*** error: number of liquid cloud optical coefficients (', &
             &  size(config%cloud_optics%liq_coeff_lw, 2), &
             &  ') does not match number expected (', nliqopticscoeffssocrates,')'
        call radiation_abort()
      end if
    else if (config%i_liq_model == iliquidmodelslingo) then
      if (size(config%cloud_optics%liq_coeff_sw, 2) /= nliqopticscoeffsslingosw) then
        write(nulerr,'(a,i0,a,i0,a,i0,a)') &
             &  '*** error: number of shortwave liquid cloud optical coefficients (', &
             &  size(config%cloud_optics%liq_coeff_sw, 2), &
             &  ') does not match number expected (', nliqopticscoeffsslingosw,')'
        call radiation_abort()
      end if
      if (size(config%cloud_optics%liq_coeff_lw, 2) /= nliqopticscoeffslindnerlilw) then
        write(nulerr,'(a,i0,a,i0,a,i0,a)') &
             &  '*** error: number of longwave liquid cloud optical coefficients (', &
             &  size(config%cloud_optics%liq_coeff_lw, 2), &
             &  ') does not match number expected (', nliqopticscoeffslindnerlilw,')'
        call radiation_abort()
      end if
    end if

    ! check ice coefficients
    if (size(config%cloud_optics%ice_coeff_lw, 1) /= config%n_bands_lw) then
      write(nulerr,'(a,i0,a,i0,a)') &
           &  '*** error: number of longwave bands for ice particles (', &
           &  size(config%cloud_optics%ice_coeff_lw, 1), &
           &  ') does not match number for gases (', config%n_bands_lw, ')'
      call radiation_abort()
    end if
    if (size(config%cloud_optics%ice_coeff_sw, 1) /= config%n_bands_sw) then
      write(nulerr,'(a,i0,a,i0,a)') &
           &  '*** error: number of shortwave bands for ice particles (', &
           &  size(config%cloud_optics%ice_coeff_sw, 1), &
           &  ') does not match number for gases (', config%n_bands_sw, ')'
      call radiation_abort()
    end if

    if (config%i_ice_model == iicemodelfu) then
      if (size(config%cloud_optics%ice_coeff_lw, 2) &
           &  /= niceopticscoeffsfulw) then
        write(nulerr,'(a,i0,a,i0,a,i0,a)') &
             &  '*** error: number of lw ice-particle optical coefficients (', &
             &  size(config%cloud_optics%ice_coeff_lw, 2), &
             &  ') does not match number expected (', niceopticscoeffsfulw,')'
        call radiation_abort()
      end if
      if (size(config%cloud_optics%ice_coeff_sw, 2) &
           &  /= niceopticscoeffsfusw) then
        write(nulerr,'(a,i0,a,i0,a,i0,a)') &
             &  '*** error: number of sw ice-particle optical coefficients (', &
             &  size(config%cloud_optics%ice_coeff_sw, 2), &
             &  ') does not match number expected (', niceopticscoeffsfusw,')'
        call radiation_abort()
      end if
    else if (config%i_ice_model == iicemodelbaran &
         &  .and. size(config%cloud_optics%ice_coeff_lw, 2) &
         &  /= niceopticscoeffsbaran) then
      write(nulerr,'(a,i0,a,i0,a,i0,a)') &
           &  '*** error: number of ice-particle optical coefficients (', &
           &  size(config%cloud_optics%ice_coeff_lw, 2), &
           &  ') does not match number expected (', niceopticscoeffsbaran,')'
      call radiation_abort()
    else if (config%i_ice_model == iicemodelbaran2016 &
         &  .and. size(config%cloud_optics%ice_coeff_lw, 2) &
         &  /= niceopticscoeffsbaran2016) then
      write(nulerr,'(a,i0,a,i0,a,i0,a)') &
           &  '*** error: number of ice-particle optical coefficients (', &
           &  size(config%cloud_optics%ice_coeff_lw, 2), &
           &  ') does not match number expected (', niceopticscoeffsbaran2016,')'
      call radiation_abort()
    else if (config%i_ice_model == iicemodelbaran2017) then
      if (size(config%cloud_optics%ice_coeff_lw, 2) &
           &  /= niceopticscoeffsbaran2017) then
        write(nulerr,'(a,i0,a,i0,a,i0,a)') &
             &  '*** error: number of ice-particle optical coefficients (', &
             &  size(config%cloud_optics%ice_coeff_lw, 2), &
             &  ') does not match number expected (', niceopticscoeffsbaran2017,')'
        call radiation_abort()
      else if (.not. allocated(config%cloud_optics%ice_coeff_gen)) then
        write(nulerr,'(a)') &
             &  '*** error: coeff_gen needed for baran-2017 ice optics parameterization'
        call radiation_abort()
      else if (size(config%cloud_optics%ice_coeff_gen) &
           &  /= niceopticsgeneralcoeffsbaran2017) then
        write(nulerr,'(a,i0,a,i0,a,i0,a)') &
             &  '*** error: number of general ice-particle optical coefficients (', &
             &  size(config%cloud_optics%ice_coeff_gen), &
             &  ') does not match number expected (', niceopticsgeneralcoeffsbaran2017,')'
        call radiation_abort()
      end if
    else if (config%i_ice_model == iicemodelyi) then
      if (size(config%cloud_optics%ice_coeff_lw, 2) &
           &  /= niceopticscoeffsyilw) then
        write(nulerr,'(a,i0,a,i0,a,i0,a)') &
             &  '*** error: number of lw ice-particle optical coefficients (', &
             &  size(config%cloud_optics%ice_coeff_lw, 2), &
             &  ') does not match number expected (', niceopticscoeffsyilw,')'
        call radiation_abort()
      end if
      if (size(config%cloud_optics%ice_coeff_sw, 2) &
           &  /= niceopticscoeffsyisw) then
        write(nulerr,'(a,i0,a,i0,a,i0,a)') &
             &  '*** error: number of sw ice-particle optical coefficients (', &
             &  size(config%cloud_optics%ice_coeff_sw, 2), &
             &  ') does not match number expected (', niceopticscoeffsyisw,')'
        call radiation_abort()
      end if
    end if

    if (lhook) call dr_hook('radiation_cloud_optics:setup_cloud_optics',1,hook_handle)

  end subroutine setup_cloud_optics


  !---------------------------------------------------------------------
  ! compute cloud optical properties
  subroutine cloud_optics(nlev,istartcol,iendcol, &
       &  config, thermodynamics, cloud, &
       &  od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
       &  od_sw_cloud, ssa_sw_cloud, g_sw_cloud)

    use parkind1, only           : jprb
    use ecradhook,  only           : lhook, dr_hook, jphook

    use radiation_io,     only : nulout, nulerr, radiation_abort
    use radiation_config, only : config_type, iicemodelfu, iicemodelbaran, &
         &                       iicemodelbaran2016, iicemodelbaran2017, &
         &                       iicemodelyi, &
         &                       iliquidmodelsocrates, iliquidmodelslingo
    use radiation_thermodynamics, only    : thermodynamics_type
    use radiation_cloud, only             : cloud_type
    use radiation_constants, only         : accelduetogravity
    use radiation_ice_optics_fu, only     : calc_ice_optics_fu_sw, &
         &                                  calc_ice_optics_fu_lw
    use radiation_ice_optics_baran, only  : calc_ice_optics_baran, &
         &                                  calc_ice_optics_baran2016
    use radiation_ice_optics_baran2017, only  : calc_ice_optics_baran2017
    use radiation_ice_optics_yi, only     : calc_ice_optics_yi_sw, &
         &                                  calc_ice_optics_yi_lw
    use radiation_liquid_optics_socrates, only:calc_liq_optics_socrates
    use radiation_liquid_optics_slingo, only:calc_liq_optics_slingo, &
         &                                   calc_liq_optics_lindner_li

    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in), target :: config
    type(thermodynamics_type),intent(in)  :: thermodynamics
    type(cloud_type),   intent(in)        :: cloud

    ! layer optical depth, single scattering albedo and g factor of
    ! clouds in each longwave band, where the latter two
    ! variables are only defined if cloud longwave scattering is
    ! enabled (otherwise both are treated as zero).
    real(jprb), dimension(config%n_bands_lw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_lw_cloud
    real(jprb), dimension(config%n_bands_lw_if_scattering,nlev,istartcol:iendcol), &
         &   intent(out) :: ssa_lw_cloud, g_lw_cloud

    ! layer optical depth, single scattering albedo and g factor of
    ! clouds in each shortwave band
    real(jprb), dimension(config%n_bands_sw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_sw_cloud, ssa_sw_cloud, g_sw_cloud

    ! longwave and shortwave optical depth, scattering optical depth
    ! and asymmetry factor, for liquid and ice in all bands but a
    ! single cloud layer
    real(jprb), dimension(config%n_bands_lw) :: &
         &  od_lw_liq, scat_od_lw_liq, g_lw_liq, &
         &  od_lw_ice, scat_od_lw_ice, g_lw_ice
    real(jprb), dimension(config%n_bands_sw) :: &
         &  od_sw_liq, scat_od_sw_liq, g_sw_liq, &
         &  od_sw_ice, scat_od_sw_ice, g_sw_ice

    ! in-cloud water path of cloud liquid or ice (i.e. liquid or ice
    ! gridbox-mean water path divided by cloud fraction); kg m-2
    real(jprb) :: lwp_in_cloud, iwp_in_cloud

    ! full-level temperature (k)
    real(jprb) :: temperature

    ! factor to convert gridbox-mean mixing ratio to in-cloud water
    ! path
    real(jprb) :: factor

    integer    :: jcol, jlev, jb

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloud_optics:cloud_optics',0,hook_handle)

    if (config%iverbose >= 2) then
      write(nulout,'(a)') 'computing cloud absorption/scattering properties'
    end if

    associate(ho => config%cloud_optics)

      !$acc data present (config, thermodynamics, thermodynamics%pressure_hl, &
      !$acc              cloud, cloud%fraction, cloud%q_liq, cloud%q_ice, &
      !$acc              cloud%re_liq, cloud%re_ice) &
      !$acc      present(od_lw_cloud, g_lw_cloud, ssa_lw_cloud, od_sw_cloud, g_sw_cloud, ssa_sw_cloud)

      ! array-wise assignment
      !$acc parallel default(none) async(1)
      !$acc loop gang collapse(2)
      do jcol=istartcol, iendcol
        do jlev=1, nlev
          !$acc loop vector
          do jb=1, config%n_bands_sw
            od_sw_cloud(jb,jlev,jcol) = 0.0_jprb
            ssa_sw_cloud(jb,jlev,jcol) = 0.0_jprb
            g_sw_cloud(jb,jlev,jcol) = 0.0_jprb
          end do
          !$acc loop vector
          do jb=1, config%n_bands_lw
            od_lw_cloud(jb,jlev,jcol) = 0.0_jprb
          end do
          !$acc loop vector
          do jb=1, config%n_bands_lw_if_scattering
            ssa_lw_cloud(jb,jlev,jcol) = 0.0_jprb
            g_lw_cloud(jb,jlev,jcol) = 0.0_jprb
          end do
        end do
      end do
      !$acc end parallel

      !$acc parallel default(none) async(1)
      !$acc loop seq
      do jlev = 1,nlev
        !$acc loop gang vector &
        !$acc private(od_lw_liq, scat_od_lw_liq, g_lw_liq, od_lw_ice, scat_od_lw_ice, g_lw_ice)  &
        !$acc private(od_sw_liq, scat_od_sw_liq, g_sw_liq, od_sw_ice, scat_od_sw_ice, g_sw_ice)
        do jcol = istartcol,iendcol
          ! only do anything if cloud is present (assume that
          ! cloud%crop_cloud_fraction has already been called)
          if (cloud%fraction(jcol,jlev) > 0.0_jprb) then

            ! compute in-cloud liquid and ice water path
            if (config%is_homogeneous) then
              ! homogeneous solvers assume cloud fills the box
              ! horizontally, so we don't divide by cloud fraction
              factor = ( thermodynamics%pressure_hl(jcol,jlev+1)    &
                   &  -thermodynamics%pressure_hl(jcol,jlev  )  ) &
                   &  / accelduetogravity
            else
              factor = ( thermodynamics%pressure_hl(jcol,jlev+1)    &
                   &  -thermodynamics%pressure_hl(jcol,jlev  )  ) &
                   &  / (accelduetogravity * cloud%fraction(jcol,jlev))
            end if
            lwp_in_cloud = factor * cloud%q_liq(jcol,jlev)
            iwp_in_cloud = factor * cloud%q_ice(jcol,jlev)

            ! only compute liquid properties if liquid cloud is
            ! present
            if (lwp_in_cloud > 0.0_jprb) then

              if (config%i_liq_model == iliquidmodelsocrates) then

                ! compute longwave properties
                call calc_liq_optics_socrates(config%n_bands_lw, &
                    &  config%cloud_optics%liq_coeff_lw, &
                    &  lwp_in_cloud, cloud%re_liq(jcol,jlev), &
                    &  od_lw_liq, scat_od_lw_liq, g_lw_liq)
                ! compute shortwave properties
                call calc_liq_optics_socrates(config%n_bands_sw, &
                    &  config%cloud_optics%liq_coeff_sw, &
                    &  lwp_in_cloud, cloud%re_liq(jcol,jlev), &
                    &  od_sw_liq, scat_od_sw_liq, g_sw_liq)

              else if (config%i_liq_model == iliquidmodelslingo) then
                ! compute longwave properties
                call calc_liq_optics_lindner_li(config%n_bands_lw, &
                    &  config%cloud_optics%liq_coeff_lw, &
                    &  lwp_in_cloud, cloud%re_liq(jcol,jlev), &
                    &  od_lw_liq, scat_od_lw_liq, g_lw_liq)
                ! compute shortwave properties
                call calc_liq_optics_slingo(config%n_bands_sw, &
                    &  config%cloud_optics%liq_coeff_sw, &
                    &  lwp_in_cloud, cloud%re_liq(jcol,jlev), &
                    &  od_sw_liq, scat_od_sw_liq, g_sw_liq)
              else
                write(nulerr,*) '*** error: unknown liquid model with code', &
                    &          config%i_liq_model
                call radiation_abort()
              end if


              ! delta-eddington scaling in the shortwave only
              if (.not. config%do_sw_delta_scaling_with_gases) then
                call delta_eddington_scat_od(od_sw_liq, scat_od_sw_liq, g_sw_liq)
              end if
              ! originally delta-eddington has been off in ecrad for
              ! liquid clouds in the longwave, but it should be on
              !call delta_eddington_scat_od(od_lw_liq, scat_od_lw_liq, g_lw_liq)

            else
              ! liquid not present: set properties to zero
              od_lw_liq = 0.0_jprb
              scat_od_lw_liq = 0.0_jprb
              g_lw_liq = 0.0_jprb

              od_sw_liq = 0.0_jprb
              scat_od_sw_liq = 0.0_jprb
              g_sw_liq = 0.0_jprb
            end if ! liquid present

            ! only compute ice properties if ice cloud is present
            if (iwp_in_cloud > 0.0_jprb) then

              if (config%i_ice_model == iicemodelbaran) then
                ! compute longwave properties
                call calc_ice_optics_baran(config%n_bands_lw, &
                    &  config%cloud_optics%ice_coeff_lw, &
                    &  iwp_in_cloud, cloud%q_ice(jcol,jlev), &
                    &  od_lw_ice, scat_od_lw_ice, g_lw_ice)
                ! compute shortwave properties
                call calc_ice_optics_baran(config%n_bands_sw, &
                    &  config%cloud_optics%ice_coeff_sw, &
                    &  iwp_in_cloud, cloud%q_ice(jcol,jlev), &
                    &  od_sw_ice, scat_od_sw_ice, g_sw_ice)
              else if (config%i_ice_model == iicemodelbaran2016) then
                temperature = 0.5_jprb * (thermodynamics%temperature_hl(jcol,jlev) &
                    &                   +thermodynamics%temperature_hl(jcol,jlev+1))
                ! compute longwave properties
                call calc_ice_optics_baran2016(config%n_bands_lw, &
                    &  config%cloud_optics%ice_coeff_lw, &
                    &  iwp_in_cloud, cloud%q_ice(jcol,jlev), &
                    &  temperature, &
                    &  od_lw_ice, scat_od_lw_ice, g_lw_ice)
                ! compute shortwave properties
                call calc_ice_optics_baran2016(config%n_bands_sw, &
                    &  config%cloud_optics%ice_coeff_sw, &
                    &  iwp_in_cloud, cloud%q_ice(jcol,jlev), &
                    &  temperature, &
                    &  od_sw_ice, scat_od_sw_ice, g_sw_ice)
              else if (config%i_ice_model == iicemodelbaran2017) then
                temperature = 0.5_jprb * (thermodynamics%temperature_hl(jcol,jlev) &
                    &                   +thermodynamics%temperature_hl(jcol,jlev+1))
                ! compute longwave properties
                call calc_ice_optics_baran2017(config%n_bands_lw, &
                    &  config%cloud_optics%ice_coeff_gen, &
                    &  config%cloud_optics%ice_coeff_lw, &
                    &  iwp_in_cloud, cloud%q_ice(jcol,jlev), &
                    &  temperature, &
                    &  od_lw_ice, scat_od_lw_ice, g_lw_ice)
                ! compute shortwave properties
                call calc_ice_optics_baran2017(config%n_bands_sw, &
                    &  config%cloud_optics%ice_coeff_gen, &
                    &  config%cloud_optics%ice_coeff_sw, &
                    &  iwp_in_cloud, cloud%q_ice(jcol,jlev), &
                    &  temperature, &
                    &  od_sw_ice, scat_od_sw_ice, g_sw_ice)
              else if (config%i_ice_model == iicemodelfu) then

                ! compute longwave properties
                call calc_ice_optics_fu_lw(config%n_bands_lw, &
                    &  config%cloud_optics%ice_coeff_lw, &
                    &  iwp_in_cloud, cloud%re_ice(jcol,jlev), &
                    &  od_lw_ice, scat_od_lw_ice, g_lw_ice)
                if (config%do_fu_lw_ice_optics_bug) then
                  ! reproduce bug in old ifs scheme
                  scat_od_lw_ice = od_lw_ice - scat_od_lw_ice
                end if
                ! compute shortwave properties
                call calc_ice_optics_fu_sw(config%n_bands_sw, &
                    &  config%cloud_optics%ice_coeff_sw, &
                    &  iwp_in_cloud, cloud%re_ice(jcol,jlev), &
                    &  od_sw_ice, scat_od_sw_ice, g_sw_ice)

              else if (config%i_ice_model == iicemodelyi) then
                ! compute longwave properties
                call calc_ice_optics_yi_lw(config%n_bands_lw, &
                    &  config%cloud_optics%ice_coeff_lw, &
                    &  iwp_in_cloud, cloud%re_ice(jcol,jlev), &
                    &  od_lw_ice, scat_od_lw_ice, g_lw_ice)
                ! compute shortwave properties
                call calc_ice_optics_yi_sw(config%n_bands_sw, &
                    &  config%cloud_optics%ice_coeff_sw, &
                    &  iwp_in_cloud, cloud%re_ice(jcol,jlev), &
                    &  od_sw_ice, scat_od_sw_ice, g_sw_ice)
              else
                write(nulerr,*) '*** error: unknown ice model with code', &
                    &          config%i_ice_model
                call radiation_abort()
              end if


              ! delta-eddington scaling in both longwave and shortwave
              ! (assume that particles are larger than wavelength even
              ! in longwave)
              if (.not. config%do_sw_delta_scaling_with_gases) then
                call delta_eddington_scat_od(od_sw_ice, scat_od_sw_ice, g_sw_ice)
              end if
              call delta_eddington_scat_od(od_lw_ice, scat_od_lw_ice, g_lw_ice)

            else
              ! ice not present: set properties to zero
              od_lw_ice = 0.0_jprb
              scat_od_lw_ice = 0.0_jprb
              g_lw_ice = 0.0_jprb

              od_sw_ice = 0.0_jprb
              scat_od_sw_ice = 0.0_jprb
              g_sw_ice = 0.0_jprb
            end if ! ice present

            ! combine liquid and ice
            if (config%do_lw_cloud_scattering) then
! added for dwd (2020)
!nec$ shortloop
              !$acc loop seq
              do jb = 1, config%n_bands_lw
                od_lw_cloud(jb,jlev,jcol) = od_lw_liq(jb) + od_lw_ice(jb)
                if (scat_od_lw_liq(jb)+scat_od_lw_ice(jb) > 0.0_jprb) then
                  g_lw_cloud(jb,jlev,jcol) = (g_lw_liq(jb) * scat_od_lw_liq(jb) &
                    &  + g_lw_ice(jb) * scat_od_lw_ice(jb)) &
                    &  / (scat_od_lw_liq(jb)+scat_od_lw_ice(jb))
                else
                  g_lw_cloud(jb,jlev,jcol) = 0.0_jprb
                end if
                ssa_lw_cloud(jb,jlev,jcol) = (scat_od_lw_liq(jb) + scat_od_lw_ice(jb)) &
                  &                    / (od_lw_liq(jb) + od_lw_ice(jb))
              end do
            else
              ! if longwave scattering is to be neglected then the
              ! best approximation is to set the optical depth equal
              ! to the absorption optical depth
! added for dwd (2020)
!nec$ shortloop
              !$acc loop seq
              do jb = 1, config%n_bands_lw
                od_lw_cloud(jb,jlev,jcol) = od_lw_liq(jb) - scat_od_lw_liq(jb) &
                      &                   + od_lw_ice(jb) - scat_od_lw_ice(jb)
              end do
            end if
! added for dwd (2020)
!nec$ shortloop
            !$acc loop seq
            do jb = 1, config%n_bands_sw
              od_sw_cloud(jb,jlev,jcol) = od_sw_liq(jb) + od_sw_ice(jb)
              g_sw_cloud(jb,jlev,jcol) = (g_sw_liq(jb) * scat_od_sw_liq(jb) &
                &  + g_sw_ice(jb) * scat_od_sw_ice(jb)) &
                &  / (scat_od_sw_liq(jb) + scat_od_sw_ice(jb))
              ssa_sw_cloud(jb,jlev,jcol) &
                &  = (scat_od_sw_liq(jb) + scat_od_sw_ice(jb)) / (od_sw_liq(jb) + od_sw_ice(jb))
            end do
          end if ! cloud present
        end do ! loop over column
      end do ! loop over level

      !$acc end parallel
      !$acc wait
      !$acc end data ! present

    end associate

    if (lhook) call dr_hook('radiation_cloud_optics:cloud_optics',1,hook_handle)

  end subroutine cloud_optics

end module radiation_cloud_optics
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

