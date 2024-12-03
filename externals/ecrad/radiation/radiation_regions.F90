! # 1 "radiation/radiation_regions.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_regions.f90"
! radiation_regions.f90 -- properties of horizontal regions in tripleclouds & spartacus
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
!   2017-07-14  r. hogan  incorporate gamma distribution option
!   2018-10-06  r. hogan  merged from radiation_optical_depth_scaling.h and radiation_overlap.f90

module radiation_regions

  implicit none

  public :: calc_region_properties

contains

  !---------------------------------------------------------------------
  ! compute the optical depth scalings for the optically "thick" and
  ! "thin" regions of a tripleclouds representation of a sub-grid pdf
  ! of cloud optical depth. following shonk and hogan (2008), the 16th
  ! percentile is used for the thin region, and the formulas estimate
  ! this for both lognormal and gamma distributions. however, an
  ! adjustment is needed for the gamma distribution at large
  ! fractional standard deviations.
  subroutine calc_region_properties(nlev, nreg, istartcol, iendcol, do_gamma, &
       &  cloud_fraction, frac_std, reg_fracs, od_scaling, cloud_fraction_threshold)

    use parkind1,     only : jprb
    use ecradhook,      only : lhook, dr_hook, jphook
    use radiation_io, only : nulerr, radiation_abort

    ! minimum od_scaling in the case of a gamma distribution
    real(jprb), parameter :: mingammaodscaling = 0.025_jprb

    ! at large fractional standard deviations (fsds), we cannot
    ! capture the behaviour of a gamma distribution with two equally
    ! weighted points; we need to weight the first ("lower") point
    ! more.  the weight of the first point is normally 0.5, but for
    ! fsds between 1.5 and 3.725 it rises linearly to 0.9, and for
    ! higher fsd it is capped at this value.  the weight of the second
    ! point is one minus the first point.
    real(jprb), parameter :: minlowerfrac      = 0.5_jprb
    real(jprb), parameter :: maxlowerfrac      = 0.9_jprb
    real(jprb), parameter :: fsdatminlowerfrac = 1.5_jprb
    real(jprb), parameter :: fsdatmaxlowerfrac = 3.725_jprb
    ! between fsdatminlowerfrac and fsdatmaxlowerfrac,
    ! lowerfrac=lowerfracfsdintercept+fsd*lowerfracfsdgradient
    real(jprb), parameter :: lowerfracfsdgradient &
         &  = (maxlowerfrac-minlowerfrac) / (fsdatmaxlowerfrac-fsdatminlowerfrac)
    real(jprb), parameter :: lowerfracfsdintercept &
         &  = minlowerfrac - fsdatminlowerfrac*lowerfracfsdgradient

    ! number of levels and regions
    integer, intent(in) :: nlev, nreg

    ! range of columns to process
    integer, intent(in) :: istartcol, iendcol

    ! do we do a lognormal or gamma distribution?
    logical, intent(in) :: do_gamma

    ! cloud fraction, i.e. the fraction of the gridbox assigned to all
    ! regions numbered 2 and above (region 1 is clear sky)
    real(jprb), intent(in), dimension(:,:)  :: cloud_fraction ! (ncol,nlev)

    ! fractional standard deviation of in-cloud water content
    real(jprb), intent(in), dimension(:,:)  :: frac_std       ! (ncol,nlev)

    ! fractional area coverage of each region
    real(jprb), intent(out) :: reg_fracs(1:nreg,nlev,istartcol:iendcol)

    ! optical depth scaling for the cloudy regions
    real(jprb), intent(out) :: od_scaling(2:nreg,nlev,istartcol:iendcol)

    ! regions smaller than this are ignored
    real(jprb), intent(in), optional :: cloud_fraction_threshold

    ! in case the user doesn't supply cloud_fraction_threshold we use
    ! a default value
    real(jprb) :: frac_threshold

    ! loop indices
    integer :: jcol, jlev

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_region_properties:calc_region_properties',0,hook_handle)

    if (present(cloud_fraction_threshold)) then
      frac_threshold = cloud_fraction_threshold
    else
      frac_threshold = 1.0e-20_jprb
    end if

    if (nreg == 2) then
      ! only one clear-sky and one cloudy region: cloudy region is
      ! homogeneous
      reg_fracs(2,1:nlev,istartcol:iendcol)  = transpose(cloud_fraction(istartcol:iendcol,1:nlev))
      reg_fracs(1,1:nlev,istartcol:iendcol)  = 1.0_jprb - reg_fracs(2,1:nlev,istartcol:iendcol)
      od_scaling(2,1:nlev,istartcol:iendcol) = 1.0_jprb

    else if (nreg == 3) then
      ! two cloudy regions with optical depth scaled by 1-x and
      ! 1+x.
      ! simple version which fails when fractional_std >= 1:
      !od_scaling(2) = 1.0_jprb-cloud%fractional_std(jcol,jlev)
      ! according to shonk and hogan (2008), 1-fsd should correspond to
      ! the 16th percentile. 
      if (.not. do_gamma) then
        ! if we treat the distribution as a lognormal such that the
        ! equivalent normal has a mean mu and standard deviation
        ! sigma, then the 16th percentile of the lognormal is very
        ! close to exp(mu-sigma).
        do jcol = istartcol,iendcol
          do jlev = 1,nlev
            if (cloud_fraction(jcol,jlev) < frac_threshold) then
              reg_fracs(1,jlev,jcol)       = 1.0_jprb
              reg_fracs(2:3,jlev,jcol)  = 0.0_jprb
              od_scaling(2:3,jlev,jcol) = 1.0_jprb
            else
              reg_fracs(1,jlev,jcol) = 1.0_jprb - cloud_fraction(jcol,jlev)
              reg_fracs(2:3,jlev,jcol) = cloud_fraction(jcol,jlev)*0.5_jprb
              od_scaling(2,jlev,jcol) &
                   &  = exp(-sqrt(log(frac_std(jcol,jlev)**2+1))) &
                   &  / sqrt(frac_std(jcol,jlev)**2+1)
              od_scaling(3,jlev,jcol) = 2.0_jprb - od_scaling(2,jlev,jcol)
            end if
          end do
        end do
      else
        ! if we treat the distribution as a gamma then the 16th
        ! percentile is close to the following.  note that because it
        ! becomes vanishingly small for fsd >~ 2, we have a lower
        ! limit of 1/40, and for higher fsds reduce the fractional
        ! cover of the denser region - see region_fractions routine
        ! below
        do jcol = istartcol,iendcol
          do jlev = 1,nlev
            if (cloud_fraction(jcol,jlev) < frac_threshold) then
              reg_fracs(1,jlev,jcol)    = 1.0_jprb
              reg_fracs(2:3,jlev,jcol)  = 0.0_jprb
              od_scaling(2:3,jlev,jcol) = 1.0_jprb
            else
              ! fraction of the clear-sky region
              reg_fracs(1,jlev,jcol) = 1.0_jprb - cloud_fraction(jcol,jlev)
!#define old_gamma_region_behaviour 1
! # 172 "radiation/radiation_regions.f90"
              ! improved behaviour.
              ! fraction and optical-depth scaling of the lower of the
              ! two cloudy regions
              reg_fracs(2,jlev,jcol) = cloud_fraction(jcol,jlev) &
                   &  * max(minlowerfrac, min(maxlowerfrac, &
                   &  lowerfracfsdintercept + frac_std(jcol,jlev)*lowerfracfsdgradient))
              od_scaling(2,jlev,jcol) = mingammaodscaling &
                   &  + (1.0_jprb - mingammaodscaling) &
                   &    * exp(-frac_std(jcol,jlev)*(1.0_jprb + 0.5_jprb*frac_std(jcol,jlev) &
                   &                     *(1.0_jprb+0.5_jprb*frac_std(jcol,jlev))))

              ! fraction of the upper of the two cloudy regions
              reg_fracs(3,jlev,jcol) = 1.0_jprb-reg_fracs(1,jlev,jcol)-reg_fracs(2,jlev,jcol)
              ! ensure conservation of the mean optical depth
              od_scaling(3,jlev,jcol) = (cloud_fraction(jcol,jlev) &
                   &  -reg_fracs(2,jlev,jcol)*od_scaling(2,jlev,jcol)) / reg_fracs(3,jlev,jcol)
            end if
          end do
        end do
      end if
    else ! nreg > 3
      write(nulerr,'(a)') '*** error: only 2 or 3 regions may be specified'
      call radiation_abort()
    end if

    if (lhook) call dr_hook('radiation_region_properties:calc_region_properties',1,hook_handle)

  end subroutine calc_region_properties

end module radiation_regions

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

