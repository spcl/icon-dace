! # 1 "radiation/radiation_two_stream.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_two_stream.f90"
! this file has been modified for the use in icon

! radiation_two_stream.f90 - compute two-stream coefficients
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
!   2017-05-04  p dueben/r hogan  use jprd where double precision essential
!   2017-07-12  r hogan  optimized lw coeffs in low optical depth case
!   2017-07-26  r hogan  added calc_frac_scattered_diffuse_sw routine
!   2017-10-23  r hogan  renamed single-character variables
!   2021-02-19  r hogan  security for shortwave singularity
!   2022-11-22  p ukkonen/r hogan  single precision uses no double precision
!   2023-09-28  r hogan  increased security for single-precision sw "k"


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
! # 27 "radiation/radiation_two_stream.f90" 2

module radiation_two_stream

  use parkind1, only : jprb, jprd

  implicit none
  public

  ! elsasser's factor: the effective factor by which the zenith
  ! optical depth needs to be multiplied to account for longwave
  ! transmission at all angles through the atmosphere.  alternatively
  ! think of acos(1/lw_diffusivity) to be the effective zenith angle
  ! of longwave radiation.
  real(jprd), parameter :: lwdiffusivity   = 1.66_jprd
  real(jprb), parameter :: lwdiffusivitywp = 1.66_jprb ! working precision version

  ! the routines in this module can be called millions of times, so
  ! calling dr hook for each one may be a significant overhead.
  ! uncomment the following to turn dr hook on.
!#define do_dr_hook_two_stream

contains

  !---------------------------------------------------------------------
  ! calculate the two-stream coefficients gamma1 and gamma2 for the
  ! longwave
  subroutine calc_two_stream_gammas_lw(ng, ssa, g, &
       &                               gamma1, gamma2)





    integer, intent(in) :: ng
    ! sngle scattering albedo and asymmetry factor:
    real(jprb), intent(in),  dimension(:) :: ssa, g
    real(jprb), intent(out), dimension(:) :: gamma1, gamma2

    real(jprb) :: factor

    integer    :: jg

    !$acc routine worker







!$acc loop worker vector 
! added for dwd (2020)
!nec$ shortloop
    do jg = 1, ng
      ! fu et al. (1997), eq 2.9 and 2.10:
      !      gamma1(jg) = lwdiffusivity * (1.0_jprb - 0.5_jprb*ssa(jg) &
      !           &                    * (1.0_jprb + g(jg)))
      !      gamma2(jg) = lwdiffusivity * 0.5_jprb * ssa(jg) &
      !           &                    * (1.0_jprb - g(jg))
      ! reduce number of multiplications
      factor = (lwdiffusivity * 0.5_jprb) * ssa(jg)
      gamma1(jg) = lwdiffusivity - factor*(1.0_jprb + g(jg))
      gamma2(jg) = factor * (1.0_jprb - g(jg))
    end do





  end subroutine calc_two_stream_gammas_lw


  !---------------------------------------------------------------------
  ! calculate the two-stream coefficients gamma1-gamma4 in the
  ! shortwave
  subroutine calc_two_stream_gammas_sw(ng, mu0, ssa, g, &
       &                               gamma1, gamma2, gamma3)





    integer, intent(in) :: ng
    ! cosine of solar zenith angle, single scattering albedo and
    ! asymmetry factor:
    real(jprb), intent(in)                :: mu0
    real(jprb), intent(in),  dimension(:) :: ssa, g
    real(jprb), intent(out), dimension(:) :: gamma1, gamma2, gamma3

    real(jprb) :: factor

    integer    :: jg







    !$acc routine worker

    ! zdunkowski "pifm" (zdunkowski et al., 1980; contributions to
    ! atmospheric physics 53, 147-66)
!$acc loop worker vector private(factor)
! added for dwd (2020)
!nec$ shortloop
    do jg = 1, ng
      !      gamma1(jg) = 2.0_jprb  - ssa(jg) * (1.25_jprb + 0.75_jprb*g(jg))
      !      gamma2(jg) = 0.75_jprb *(ssa(jg) * (1.0_jprb - g(jg)))
      !      gamma3(jg) = 0.5_jprb  - (0.75_jprb*mu0)*g(jg)
      ! optimized version:
      factor = 0.75_jprb*g(jg)
      gamma1(jg) = 2.0_jprb  - ssa(jg) * (1.25_jprb + factor)
      gamma2(jg) = ssa(jg) * (0.75_jprb - factor)
      gamma3(jg) = 0.5_jprb  - mu0*factor
    end do





  end subroutine calc_two_stream_gammas_sw


  !---------------------------------------------------------------------
  ! compute the longwave reflectance and transmittance to diffuse
  ! radiation using the meador & weaver formulas, as well as the
  ! upward flux at the top and the downward flux at the base of the
  ! layer due to emission from within the layer assuming a linear
  ! variation of planck function within the layer.
  subroutine calc_reflectance_transmittance_lw(ng, &
       &    od, gamma1, gamma2, planck_top, planck_bot, &
       &    reflectance, transmittance, source_up, source_dn)





    implicit none
    
    integer, intent(in) :: ng

    ! optical depth and single scattering albedo
    real(jprb), intent(in), dimension(ng) :: od

    ! the two transfer coefficients from the two-stream
    ! differentiatial equations (computed by
    ! calc_two_stream_gammas_lw)
    real(jprb), intent(in), dimension(ng) :: gamma1, gamma2

    ! the planck terms (functions of temperature) at the top and
    ! bottom of the layer
    real(jprb), intent(in), dimension(ng) :: planck_top, planck_bot

    ! the diffuse reflectance and transmittance, i.e. the fraction of
    ! diffuse radiation incident on a layer from either top or bottom
    ! that is reflected back or transmitted through
    real(jprb), intent(out), dimension(ng) :: reflectance, transmittance

    ! the upward emission at the top of the layer and the downward
    ! emission at its base, due to emission from within the layer
    real(jprb), intent(out), dimension(ng) :: source_up, source_dn

    real(jprd) :: k_exponent, reftrans_factor
    real(jprd) :: exponential  ! = exp(-k_exponent*od)
    real(jprd) :: exponential2 ! = exp(-2*k_exponent*od)

    real(jprd) :: coeff, coeff_up_top, coeff_up_bot, coeff_dn_top, coeff_dn_bot

    integer :: jg







    !$acc routine worker


!$acc loop worker vector
! added for dwd (2020)
!nec$ shortloop
    do jg = 1, ng
      if (od(jg) > 1.0e-3_jprd) then
        k_exponent = sqrt(max((gamma1(jg) - gamma2(jg)) * (gamma1(jg) + gamma2(jg)), &
             1.0e-12_jprd)) ! eq 18 of meador & weaver (1980)
        exponential = exp(-k_exponent*od(jg))
        exponential2 = exponential*exponential
        reftrans_factor = 1.0 / (k_exponent + gamma1(jg) + (k_exponent - gamma1(jg))*exponential2)
        ! meador & weaver (1980) eq. 25
        reflectance(jg) = gamma2(jg) * (1.0_jprd - exponential2) * reftrans_factor
        ! meador & weaver (1980) eq. 26
        transmittance(jg) = 2.0_jprd * k_exponent * exponential * reftrans_factor
      
        ! compute upward and downward emission assuming the planck
        ! function to vary linearly with optical depth within the layer
        ! (e.g. wiscombe , jqsrt 1976).

        ! stackhouse and stephens (jas 1991) eqs 5 & 12
        coeff = (planck_bot(jg)-planck_top(jg)) / (od(jg)*(gamma1(jg)+gamma2(jg)))
        coeff_up_top  =  coeff + planck_top(jg)
        coeff_up_bot  =  coeff + planck_bot(jg)
        coeff_dn_top  = -coeff + planck_top(jg)
        coeff_dn_bot  = -coeff + planck_bot(jg)
        source_up(jg) =  coeff_up_top - reflectance(jg) * coeff_dn_top - transmittance(jg) * coeff_up_bot
        source_dn(jg) =  coeff_dn_bot - reflectance(jg) * coeff_up_bot - transmittance(jg) * coeff_dn_top
      else
        k_exponent = sqrt(max((gamma1(jg) - gamma2(jg)) * (gamma1(jg) + gamma2(jg)), &
             1.0e-12_jprd)) ! eq 18 of meador & weaver (1980)
        reflectance(jg) = gamma2(jg) * od(jg)
        transmittance(jg) = (1.0_jprb - k_exponent*od(jg)) / (1.0_jprb + od(jg)*(gamma1(jg)-k_exponent))
        source_up(jg) = (1.0_jprb - reflectance(jg) - transmittance(jg)) &
             &       * 0.5 * (planck_top(jg) + planck_bot(jg))
        source_dn(jg) = source_up(jg)
      end if
    end do
    



  
  end subroutine calc_reflectance_transmittance_lw
  

  !---------------------------------------------------------------------
  ! compute the longwave reflectance and transmittance to diffuse
  ! radiation using the meador & weaver formulas, as well as the
  ! upward flux at the top and the downward flux at the base of the
  ! layer due to emission from within the layer assuming a linear
  ! variation of planck function within the layer; this version
  ! computes gamma1 and gamma2 within the same routine.
  subroutine calc_ref_trans_lw(ng, &
       &    od, ssa, asymmetry, planck_top, planck_bot, &
       &    reflectance, transmittance, source_up, source_dn)





    integer, intent(in) :: ng

    ! optical depth and single scattering albedo
    real(jprb), intent(in), dimension(ng) :: od

    ! single scattering albedo and asymmetry factor
    real(jprb), intent(in), dimension(ng) :: ssa, asymmetry

    ! the planck terms (functions of temperature) at the top and
    ! bottom of the layer
    real(jprb), intent(in), dimension(ng) :: planck_top, planck_bot

    ! the diffuse reflectance and transmittance, i.e. the fraction of
    ! diffuse radiation incident on a layer from either top or bottom
    ! that is reflected back or transmitted through
    real(jprb), intent(out), dimension(ng) :: reflectance, transmittance

    ! the upward emission at the top of the layer and the downward
    ! emission at its base, due to emission from within the layer
    real(jprb), intent(out), dimension(ng) :: source_up, source_dn

    ! the two transfer coefficients from the two-stream
    ! differentiatial equations
    real(jprb) :: gamma1, gamma2

    real(jprb) :: k_exponent, reftrans_factor, factor
    real(jprb) :: exponential  ! = exp(-k_exponent*od)
    real(jprb) :: exponential2 ! = exp(-2*k_exponent*od)

    real(jprb) :: coeff, coeff_up_top, coeff_up_bot, coeff_dn_top, coeff_dn_bot

    integer :: jg







    !$acc routine worker

    !$acc loop worker vector private(factor, gamma1, gamma2, k_exponent, &
    !$acc   reftrans_factor, exponential, exponential2, &
    !$acc   coeff, coeff_up_top, coeff_up_bot, coeff_dn_top, coeff_dn_bot)
    do jg = 1, ng
      factor = (lwdiffusivitywp * 0.5_jprb) * ssa(jg)
      gamma1 = lwdiffusivitywp - factor*(1.0_jprb + asymmetry(jg))
      gamma2 = factor * (1.0_jprb - asymmetry(jg))
      k_exponent = sqrt(max((gamma1 - gamma2) * (gamma1 + gamma2), &
           1.0e-12_jprb)) ! eq 18 of meador & weaver (1980)
      if (od(jg) > 1.0e-3_jprb) then
        exponential = exp(-k_exponent*od(jg))
        exponential2 = exponential*exponential
        reftrans_factor = 1.0_jprb / (k_exponent + gamma1 + (k_exponent - gamma1)*exponential2)
        ! meador & weaver (1980) eq. 25
        reflectance(jg) = gamma2 * (1.0_jprb - exponential2) * reftrans_factor
        ! meador & weaver (1980) eq. 26
        transmittance(jg) = 2.0_jprb * k_exponent * exponential * reftrans_factor
      
        ! compute upward and downward emission assuming the planck
        ! function to vary linearly with optical depth within the layer
        ! (e.g. wiscombe , jqsrt 1976).

        ! stackhouse and stephens (jas 1991) eqs 5 & 12
        coeff = (planck_bot(jg)-planck_top(jg)) / (od(jg)*(gamma1+gamma2))
        coeff_up_top  =  coeff + planck_top(jg)
        coeff_up_bot  =  coeff + planck_bot(jg)
        coeff_dn_top  = -coeff + planck_top(jg)
        coeff_dn_bot  = -coeff + planck_bot(jg)
        source_up(jg) =  coeff_up_top - reflectance(jg) * coeff_dn_top - transmittance(jg) * coeff_up_bot
        source_dn(jg) =  coeff_dn_bot - reflectance(jg) * coeff_up_bot - transmittance(jg) * coeff_dn_top
      else
        reflectance(jg) = gamma2 * od(jg)
        transmittance(jg) = (1.0_jprb - k_exponent*od(jg)) / (1.0_jprb + od(jg)*(gamma1-k_exponent))
        source_up(jg) = (1.0_jprb - reflectance(jg) - transmittance(jg)) &
             &       * 0.5 * (planck_top(jg) + planck_bot(jg))
        source_dn(jg) = source_up(jg)
      end if
    end do
    



  
  end subroutine calc_ref_trans_lw
  
  
  !---------------------------------------------------------------------
  ! compute the longwave transmittance to diffuse radiation in the
  ! no-scattering case, as well as the upward flux at the top and the
  ! downward flux at the base of the layer due to emission from within
  ! the layer assuming a linear variation of planck function within
  ! the layer.
  subroutine calc_no_scattering_transmittance_lw(ng, &
       &    od, planck_top, planck_bot, transmittance, source_up, source_dn)





    integer, intent(in) :: ng

    ! optical depth and single scattering albedo
    real(jprb), intent(in), dimension(ng) :: od

    ! the planck terms (functions of temperature) at the top and
    ! bottom of the layer
    real(jprb), intent(in), dimension(ng) :: planck_top, planck_bot

    ! the diffuse transmittance, i.e. the fraction of diffuse
    ! radiation incident on a layer from either top or bottom that is
    ! reflected back or transmitted through
    real(jprb), intent(out), dimension(ng) :: transmittance

    ! the upward emission at the top of the layer and the downward
    ! emission at its base, due to emission from within the layer
    real(jprb), intent(out), dimension(ng) :: source_up, source_dn

    real(jprb) :: coeff, coeff_up_top, coeff_up_bot, coeff_dn_top, coeff_dn_bot !, planck_mean

    integer :: jg

    !$acc routine worker








    transmittance = exp(-lwdiffusivitywp*od)


!$acc loop worker vector
    do jg = 1, ng
      ! compute upward and downward emission assuming the planck
      ! function to vary linearly with optical depth within the layer
      ! (e.g. wiscombe , jqsrt 1976).
      coeff = lwdiffusivitywp*od(jg)



      if (od(jg) > 1.0e-3_jprb) then
        ! simplified from calc_reflectance_transmittance_lw above
        coeff = (planck_bot(jg)-planck_top(jg)) / coeff
        coeff_up_top  =  coeff + planck_top(jg)
        coeff_up_bot  =  coeff + planck_bot(jg)
        coeff_dn_top  = -coeff + planck_top(jg)
        coeff_dn_bot  = -coeff + planck_bot(jg)
        source_up(jg) =  coeff_up_top - transmittance(jg) * coeff_up_bot
        source_dn(jg) =  coeff_dn_bot - transmittance(jg) * coeff_dn_top
      else
        ! linear limit at low optical depth
        source_up(jg) = coeff * 0.5_jprb * (planck_top(jg)+planck_bot(jg))
        source_dn(jg) = source_up(jg)
      end if
    end do





  end subroutine calc_no_scattering_transmittance_lw
   
   
  !---------------------------------------------------------------------
  ! compute the shortwave reflectance and transmittance to diffuse
  ! radiation using the meador & weaver formulas, as well as the
  ! "direct" reflection and transmission, which really means the rate
  ! of transfer of direct solar radiation (into a plane perpendicular
  ! to the direct beam) into diffuse upward and downward streams at
  ! the top and bottom of the layer, respectively.  finally,
  ! trans_dir_dir is the transmittance of the atmosphere to direct
  ! radiation with no scattering.
  subroutine calc_reflectance_transmittance_sw(ng, mu0, od, ssa, &
       &      gamma1, gamma2, gamma3, ref_diff, trans_diff, &
       &      ref_dir, trans_dir_diff, trans_dir_dir)
    




    integer, intent(in) :: ng

    ! cosine of solar zenith angle
    real(jprb), intent(in) :: mu0

    ! optical depth and single scattering albedo
    real(jprb), intent(in), dimension(ng) :: od, ssa

    ! the three transfer coefficients from the two-stream
    ! differentiatial equations (computed by calc_two_stream_gammas)
    real(jprb), intent(in), dimension(ng) :: gamma1, gamma2, gamma3

    ! the direct reflectance and transmittance, i.e. the fraction of
    ! incoming direct solar radiation incident at the top of a layer
    ! that is either reflected back (ref_dir) or scattered but
    ! transmitted through the layer to the base (trans_dir_diff)
    real(jprb), intent(out), dimension(ng) :: ref_dir, trans_dir_diff

    ! the diffuse reflectance and transmittance, i.e. the fraction of
    ! diffuse radiation incident on a layer from either top or bottom
    ! that is reflected back or transmitted through
    real(jprb), intent(out), dimension(ng) :: ref_diff, trans_diff

    ! transmittance of the direct been with no scattering
    real(jprb), intent(out), dimension(ng) :: trans_dir_dir

    real(jprd) :: gamma4, alpha1, alpha2, k_exponent, reftrans_factor
    real(jprb) :: exponential0 ! = exp(-od/mu0)
    real(jprd) :: exponential  ! = exp(-k_exponent*od)
    real(jprd) :: exponential2 ! = exp(-2*k_exponent*od)
    real(jprd) :: k_mu0, k_gamma3, k_gamma4
    real(jprd) :: k_2_exponential, od_over_mu0
    integer    :: jg

    ! local value of cosine of solar zenith angle, in case it needs to be
    ! tweaked to avoid near division by zero. this is intentionally in working
    ! precision (jprb) rather than fixing at double precision (jprd).
    real(jprb) :: mu0_local







    !$acc routine worker

!$acc loop worker vector private(gamma4, alpha1, alpha2, k_exponent, &
!$acc   reftrans_factor, exponential0, exponential, exponential2, k_mu0, &
!$acc   k_gamma3, k_gamma4, k_2_exponential, od_over_mu0)
! added for dwd (2020)
!nec$ shortloop
    do jg = 1, ng

        gamma4 = 1.0_jprd - gamma3(jg)
        alpha1 = gamma1(jg)*gamma4     + gamma2(jg)*gamma3(jg) ! eq. 16
        alpha2 = gamma1(jg)*gamma3(jg) + gamma2(jg)*gamma4    ! eq. 17

        k_exponent = sqrt(max((gamma1(jg) - gamma2(jg)) * (gamma1(jg) + gamma2(jg)), &
             &       1.0e-12_jprd)) ! eq 18

        ! we had a rare crash where k*mu0 was within around 1e-13 of 1,
        ! leading to ref_dir and trans_dir_diff being well outside the range
        ! 0-1. the following approach is appropriate when k_exponent is double
        ! precision and mu0_local is single precision, although work is needed
        ! to make this entire routine secure in single precision.
        mu0_local = mu0
        if (abs(1.0_jprd - k_exponent*mu0) < 1000.0_jprd * epsilon(1.0_jprd)) then
          mu0_local = mu0 * (1.0_jprb + sign(1._jprd,k_exponent*mu0-1._jprd)*10.0_jprb*epsilon(1.0_jprb))
        end if

        od_over_mu0 = max(od(jg) / mu0_local, 0.0_jprd)

        ! note that if the minimum value is reduced (e.g. to 1.0e-24)
        ! then noise starts to appear as a function of solar zenith
        ! angle
        k_mu0 = k_exponent*mu0_local
        k_gamma3 = k_exponent*gamma3(jg)
        k_gamma4 = k_exponent*gamma4
        ! check for mu0 <= 0!
        exponential0 = exp(-od_over_mu0)
        trans_dir_dir(jg) = exponential0
        exponential = exp(-k_exponent*od(jg))
        
        exponential2 = exponential*exponential
        k_2_exponential = 2.0_jprd * k_exponent * exponential
        
        reftrans_factor = 1.0_jprd / (k_exponent + gamma1(jg) + (k_exponent - gamma1(jg))*exponential2)
        
        ! meador & weaver (1980) eq. 25
        ref_diff(jg) = gamma2(jg) * (1.0_jprd - exponential2) * reftrans_factor
        
        ! meador & weaver (1980) eq. 26
        trans_diff(jg) = k_2_exponential * reftrans_factor
        
        ! here we need mu0 even though it wasn't in meador and weaver
        ! because we are assuming the incoming direct flux is defined
        ! to be the flux into a plane perpendicular to the direction of
        ! the sun, not into a horizontal plane
        reftrans_factor = mu0_local * ssa(jg) * reftrans_factor / (1.0_jprd - k_mu0*k_mu0)
        
        ! meador & weaver (1980) eq. 14, multiplying top & bottom by
        ! exp(-k_exponent*od) in case of very high optical depths
        ref_dir(jg) = reftrans_factor &
             &  * ( (1.0_jprd - k_mu0) * (alpha2 + k_gamma3) &
             &     -(1.0_jprd + k_mu0) * (alpha2 - k_gamma3)*exponential2 &
             &     -k_2_exponential*(gamma3(jg) - alpha2*mu0_local)*exponential0)
        
        ! meador & weaver (1980) eq. 15, multiplying top & bottom by
        ! exp(-k_exponent*od), minus the 1*exp(-od/mu0) term representing direct
        ! unscattered transmittance.  
        trans_dir_diff(jg) = reftrans_factor * ( k_2_exponential*(gamma4 + alpha1*mu0_local) &
            & - exponential0 &
            & * ( (1.0_jprd + k_mu0) * (alpha1 + k_gamma4) &
            &    -(1.0_jprd - k_mu0) * (alpha1 - k_gamma4) * exponential2) )

        ! final check that ref_dir + trans_dir_diff <= 1
        ref_dir(jg) = max(0.0_jprb, min(ref_dir(jg), 1.0_jprb))
        trans_dir_diff(jg) = max(0.0_jprb, min(trans_dir_diff(jg), 1.0_jprb-ref_dir(jg)))

    end do
    



 
  end subroutine calc_reflectance_transmittance_sw


  !---------------------------------------------------------------------
  ! compute the shortwave reflectance and transmittance to diffuse
  ! radiation using the meador & weaver formulas, as well as the
  ! "direct" reflection and transmission, which really means the rate
  ! of transfer of direct solar radiation (into a plane perpendicular
  ! to the direct beam) into diffuse upward and downward streams at
  ! the top and bottom of the layer, respectively.  finally,
  ! trans_dir_dir is the transmittance of the atmosphere to direct
  ! radiation with no scattering. this version incorporates the
  ! calculation of the gamma terms.
  subroutine calc_ref_trans_sw(ng, mu0, od, ssa, &
       &      asymmetry, ref_diff, trans_diff, &
       &      ref_dir, trans_dir_diff, trans_dir_dir)
    




    implicit none
    
    integer, intent(in) :: ng

    ! cosine of solar zenith angle
    real(jprb), intent(in) :: mu0

    ! optical depth and single scattering albedo
    real(jprb), intent(in), dimension(ng) :: od, ssa, asymmetry

    ! the direct reflectance and transmittance, i.e. the fraction of
    ! incoming direct solar radiation incident at the top of a layer
    ! that is either reflected back (ref_dir) or scattered but
    ! transmitted through the layer to the base (trans_dir_diff)
    real(jprb), intent(out), dimension(ng) :: ref_dir, trans_dir_diff

    ! the diffuse reflectance and transmittance, i.e. the fraction of
    ! diffuse radiation incident on a layer from either top or bottom
    ! that is reflected back or transmitted through
    real(jprb), intent(out), dimension(ng) :: ref_diff, trans_diff

    ! transmittance of the direct been with no scattering
    real(jprb), intent(out), dimension(ng) :: trans_dir_dir

    ! the three transfer coefficients from the two-stream
    ! differentiatial equations 

    real(jprb), dimension(ng) :: gamma1, gamma2, gamma3, gamma4 
    real(jprb), dimension(ng) :: alpha1, alpha2, k_exponent
    real(jprb), dimension(ng) :: exponential ! = exp(-k_exponent*od)





    
    real(jprb) :: reftrans_factor, factor
    real(jprb) :: exponential2 ! = exp(-2*k_exponent*od)
    real(jprb) :: k_mu0, k_gamma3, k_gamma4
    real(jprb) :: k_2_exponential, one_minus_kmu0_sqr
    integer    :: jg








    ! gcc 9.3 strange error: intermediate values of ~ -8000 cause a
    ! fpe when vectorizing exp(), but not in non-vectorized loop, nor
    ! with larger negative values!
    trans_dir_dir = max(-max(od * (1.0_jprb/mu0), 0.0_jprb),-1000.0_jprb)
    trans_dir_dir = exp(trans_dir_dir)

    do jg = 1, ng

      ! zdunkowski "pifm" (zdunkowski et al., 1980; contributions to
      ! atmospheric physics 53, 147-66)
      factor = 0.75_jprb*asymmetry(jg)

      gamma1(jg) = 2.0_jprb  - ssa(jg) * (1.25_jprb + factor)
      gamma2(jg) = ssa(jg) * (0.75_jprb - factor)
      gamma3(jg) = 0.5_jprb  - mu0*factor
      gamma4(jg) = 1.0_jprb - gamma3(jg)

      alpha1(jg) = gamma1(jg)*gamma4(jg) + gamma2(jg)*gamma3(jg) ! eq. 16
      alpha2(jg) = gamma1(jg)*gamma3(jg) + gamma2(jg)*gamma4(jg) ! eq. 17
      ! the following line crashes inexplicably with gfortran 8.5.0 in
      ! single precision - try a later version. note that the minimum
      ! value is needed to produce correct results for single
      ! scattering albedos very close to or equal to one.




      k_exponent(jg) = sqrt(max((gamma1(jg) - gamma2(jg)) * (gamma1(jg) + gamma2(jg)), &
           &       1.0e-12_jprb)) ! eq 18

    end do

    exponential = exp(-k_exponent*od)

    do jg = 1, ng
      k_mu0 = k_exponent(jg)*mu0
      one_minus_kmu0_sqr = 1.0_jprb - k_mu0*k_mu0
      k_gamma3 = k_exponent(jg)*gamma3(jg)
      k_gamma4 = k_exponent(jg)*gamma4(jg)
      exponential2 = exponential(jg)*exponential(jg)
      k_2_exponential = 2.0_jprb * k_exponent(jg) * exponential(jg)
      reftrans_factor = 1.0_jprb / (k_exponent(jg) + gamma1(jg) + (k_exponent(jg) - gamma1(jg))*exponential2)
        
      ! meador & weaver (1980) eq. 25
      ref_diff(jg) = gamma2(jg) * (1.0_jprb - exponential2) * reftrans_factor
      !ref_diff(jg)       = max(0.0_jprb, min(ref_diff(jg)), 1.0_jprb)

      ! meador & weaver (1980) eq. 26, with security (which is
      ! sometimes needed, but apparently not on ref_diff)
      trans_diff(jg) = max(0.0_jprb, min(k_2_exponential * reftrans_factor, 1.0_jprb-ref_diff(jg)))

      ! here we need mu0 even though it wasn't in meador and weaver
      ! because we are assuming the incoming direct flux is defined to
      ! be the flux into a plane perpendicular to the direction of the
      ! sun, not into a horizontal plane
      reftrans_factor = mu0 * ssa(jg) * reftrans_factor &
            &  / merge(one_minus_kmu0_sqr, epsilon(1.0_jprb), abs(one_minus_kmu0_sqr) > epsilon(1.0_jprb))
      
      ! meador & weaver (1980) eq. 14, multiplying top & bottom by
      ! exp(-k_exponent*od) in case of very high optical depths
      ref_dir(jg) = reftrans_factor &
           &  * ( (1.0_jprb - k_mu0) * (alpha2(jg) + k_gamma3) &
           &     -(1.0_jprb + k_mu0) * (alpha2(jg) - k_gamma3)*exponential2 &
           &     -k_2_exponential*(gamma3(jg) - alpha2(jg)*mu0)*trans_dir_dir(jg) )
        
      ! meador & weaver (1980) eq. 15, multiplying top & bottom by
      ! exp(-k_exponent*od), minus the 1*exp(-od/mu0) term
      ! representing direct unscattered transmittance.
      trans_dir_diff(jg) = reftrans_factor * ( k_2_exponential*(gamma4(jg) + alpha1(jg)*mu0) &
           & - trans_dir_dir(jg) &
           & * ( (1.0_jprb + k_mu0) * (alpha1(jg) + k_gamma4) &
           &    -(1.0_jprb - k_mu0) * (alpha1(jg) - k_gamma4) * exponential2) )

      ! final check that ref_dir + trans_dir_diff <= 1
      ref_dir(jg)        = max(0.0_jprb, min(ref_dir(jg), mu0*(1.0_jprb-trans_dir_dir(jg))))
      trans_dir_diff(jg) = max(0.0_jprb, min(trans_dir_diff(jg), mu0*(1.0_jprb-trans_dir_dir(jg))-ref_dir(jg)))
    end do

! # 798 "radiation/radiation_two_stream.f90"




 
  end subroutine calc_ref_trans_sw

  
  !---------------------------------------------------------------------
  ! compute the fraction of shortwave transmitted diffuse radiation
  ! that is scattered during its transmission, used to compute
  ! entrapment in spartacus
  subroutine calc_frac_scattered_diffuse_sw(ng, od, &
       &      gamma1, gamma2, frac_scat_diffuse)
    




    integer, intent(in) :: ng

    ! optical depth
    real(jprb), intent(in), dimension(ng) :: od

    ! the first two transfer coefficients from the two-stream
    ! differentiatial equations (computed by calc_two_stream_gammas)
    real(jprb), intent(in), dimension(ng) :: gamma1, gamma2

    ! the fraction of shortwave transmitted diffuse radiation that is
    ! scattered during its transmission
    real(jprb), intent(out), dimension(ng) :: frac_scat_diffuse

    real(jprd) :: k_exponent, reftrans_factor
    real(jprd) :: exponential  ! = exp(-k_exponent*od)
    real(jprd) :: exponential2 ! = exp(-2*k_exponent*od)
    real(jprd) :: k_2_exponential
    integer    :: jg







! added for dwd (2020)
!nec$ shortloop
    do jg = 1, ng
      ! note that if the minimum value is reduced (e.g. to 1.0e-24)
      ! then noise starts to appear as a function of solar zenith
      ! angle
      k_exponent = sqrt(max((gamma1(jg) - gamma2(jg)) * (gamma1(jg) + gamma2(jg)), &
           &       1.0e-12_jprd)) ! eq 18
      exponential = exp(-k_exponent*od(jg))
      exponential2 = exponential*exponential
      k_2_exponential = 2.0_jprd * k_exponent * exponential
        
      reftrans_factor = 1.0_jprd / (k_exponent + gamma1(jg) + (k_exponent - gamma1(jg))*exponential2)
        
      ! meador & weaver (1980) eq. 26.
      ! until 1.1.8, used lwdiffusivity instead of 2.0, although the
      ! effect is very small
      !      frac_scat_diffuse(jg) = 1.0_jprb - min(1.0_jprb,exp(-lwdiffusivity*od(jg)) &
      !           &  / max(1.0e-8_jprb, k_2_exponential * reftrans_factor))
      frac_scat_diffuse(jg) = 1.0_jprb &
           &  - min(1.0_jprb,exp(-2.0_jprb*od(jg)) &
           &  / max(1.0e-8_jprb, k_2_exponential * reftrans_factor))
    end do
    



 
  end subroutine calc_frac_scattered_diffuse_sw

end module radiation_two_stream
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

