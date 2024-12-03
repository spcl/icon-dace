! # 1 "radiation/radiation_cloud_generator_acc.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_cloud_generator_acc.f90"
! this file has been modified for the use in icon

! radiation_cloud_generator_acc.f90 - generate water-content or optical-depth scalings for mcica
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
! generate clouds for mcica using a method modified from raisanen et
! al. (2002)
! this is a copy of the original cloud_generator, that is better suited for openacc
!
! modifications
!   2018-02-22  r. hogan  call masked version of pdf sampler for speed
!   2020-03-31  r. hogan  more vectorizable version of exp-ran
!   2022-11-07  d. hupp adaptation for acc

module radiation_cloud_generator_acc

  public

contains

  !---------------------------------------------------------------------
  ! generate scaling factors for the cloud optical depth to represent
  ! cloud overlap, the overlap of internal cloud inhomogeneities and
  ! the fractional standard deviation of these inhomogeneities, for
  ! use in a monte carlo independent column approximation radiation
  ! scheme. all returned profiles contain cloud, and the total cloud
  ! cover is also returned, so the calling function can then do a
  ! weighted average of clear and cloudy skies; this is a way to
  ! reduce the monte carlo noise in profiles with low cloud cover.
  subroutine cloud_generator_acc(ng, nlev, &
    &  iseed, frac_threshold, &
    &  frac, overlap_param, decorrelation_scaling, &
    &  fractional_std, &
    &  sample_ncdf, sample_nfsd, sample_fsd1, &
    &  sample_inv_fsd_interval, sample_val, &
    &  od_scaling, total_cloud_cover, &
    &  ibegin, iend, &
    &  cum_cloud_cover, &
    &  pair_cloud_cover)

    use parkind1,                 only : jprb
    use radiation_random_numbers, only : initialize_acc, uniform_distribution_acc, iminstda0, iminstda, iminstdm

    implicit none

    ! inputs
    integer, intent(in)     :: ng    ! number of g points
    integer, intent(in)     :: nlev  ! number of model levels
    integer, intent(in)     :: iseed ! seed for random number generator

    ! only cloud fractions above this threshold are considered to be
    ! clouds
    real(jprb), intent(in)  :: frac_threshold

    ! cloud fraction on full levels
    real(jprb), intent(in)  :: frac(nlev)

    ! cloud overlap parameter for interfaces between model layers,
    ! where 0 indicates random overlap and 1 indicates maximum-random
    ! overlap
    real(jprb), intent(in)  :: overlap_param(nlev-1)

    ! overlap parameter for internal inhomogeneities
    real(jprb), intent(in)  :: decorrelation_scaling

    ! fractional standard deviation at each layer
    real(jprb), intent(in)  :: fractional_std(nlev)

    ! object for sampling from a lognormal or gamma distribution
    integer, intent(in)  :: sample_ncdf, sample_nfsd
    real(jprb), intent(in)  :: sample_fsd1, sample_inv_fsd_interval
    real(jprb), intent(in), dimension(:,:)  :: sample_val

    ! outputs

    ! cloud optical depth scaling factor, with 0 indicating clear sky
    real(jprb), intent(out) :: od_scaling(ng,nlev)

    ! total cloud cover using cloud fraction and overlap parameter
    real(jprb), intent(in) :: total_cloud_cover

    ! local variables

    ! cumulative cloud cover from toa to the base of each layer
    real(jprb), intent(in) :: cum_cloud_cover(nlev)

    ! first and last cloudy layers
    integer, intent(in) :: ibegin, iend

    ! scaled random number for finding cloud
    real(jprb) :: trigger

    ! uniform deviates between 0 and 1
    real(jprb) :: rand_top

    ! overlap parameter of inhomogeneities
    real(jprb) :: overlap_param_inhom

    real(jprb) :: rand_cloud, rand_inhom1, rand_inhom2

    integer :: itrigger

    ! loop index for model level and g-point
    integer :: jlev, jg

    ! cloud cover of a pair of layers, and amount by which cloud at
    ! next level increases total cloud cover as seen from above
    real(jprb), intent(inout), dimension(nlev-1) :: pair_cloud_cover
    real(jprb) :: overhang

    integer          :: jcloud
    ! number of contiguous cloudy layers for which to compute optical
    ! depth scaling
    integer :: n_layers_to_scale

    ! is it time to fill the od_scaling variable?
    logical :: do_fill_od_scaling

    ! variables from manual inlining of sample_vec_from_pdf
    ! index to look-up table
    integer    :: ifsd, icdf
    ! weights in bilinear interpolation
    real(jprb) :: wfsd, wcdf

    ! local variable for the acc rng
    real(kind=jprb) :: istate

    !$acc routine worker

    if (total_cloud_cover >= frac_threshold) then
      ! loop over ng columns (this loop should be paralized as soon as acc rng is available)
      !$acc loop worker vector private(istate, rand_top, trigger, itrigger, n_layers_to_scale)
      do jg = 1,ng

        !$acc loop seq 
        do jlev = 1,nlev
          od_scaling(jg,jlev) = 0.0_jprb
        end do

        ! begin manuel inline of istate = initialize_acc(iseed, jg)
        istate = real(abs(iseed),jprb)
        ! use a modified (and vectorized) c++ minstd_rand0 algorithm to populate the state
        istate = nint(mod( istate*jg*(1._jprb-0.05_jprb*jg+0.005_jprb*jg**2)*iminstda0, iminstdm))

        ! one warmup of the c++ minstd_rand algorithm
        istate = mod(iminstda * istate, iminstdm)
        ! end manuel inline of istate = initialize_acc(iseed, jg)
        rand_top = uniform_distribution_acc(istate)

        ! find the cloud top height corresponding to the current
        ! random number, and store in itrigger
        trigger = rand_top * total_cloud_cover
        itrigger = iend
        !$acc loop seq
        do jlev = ibegin,iend
          if (trigger <= cum_cloud_cover(jlev)) then
            itrigger = min(jlev, itrigger)
          end if
        end do

        ! manual inline: call generate_column_exp_ran >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ! so far our vertically contiguous cloud contains only one layer
        n_layers_to_scale = 1

        ! locate the clouds below this layer: first generate some more
        ! random numbers
        ! loop from the layer below the local cloud top down to the
        ! bottom-most cloudy layer
        !$acc loop seq private(do_fill_od_scaling, rand_cloud, rand_inhom1, overhang)
        do jlev = itrigger+1,iend+1
          do_fill_od_scaling = .false.
          if (jlev <= iend) then
            rand_cloud = uniform_distribution_acc(istate)
            if (n_layers_to_scale > 0) then
              ! there is a cloud above, in which case the probability
              ! of cloud in the layer below is as follows
              if (rand_cloud*frac(jlev-1) &
                  &  < frac(jlev) + frac(jlev-1) - pair_cloud_cover(jlev-1)) then
                ! add another cloudy layer
                n_layers_to_scale = n_layers_to_scale + 1
              else
                ! reached the end of a contiguous set of cloudy layers and
                ! will compute the optical depth scaling immediately.
                do_fill_od_scaling = .true.
              end if
            else
              overhang = cum_cloud_cover(jlev)-cum_cloud_cover(jlev-1)
              ! there is clear-sky above, in which case the
              ! probability of cloud in the layer below is as follows
              if (rand_cloud*(cum_cloud_cover(jlev-1) - frac(jlev-1)) &
                  &  < pair_cloud_cover(jlev-1) - overhang - frac(jlev-1)) then
                ! a new cloud top
                n_layers_to_scale = 1
              end if
            end if
          else
            ! we are at the bottom of the cloudy layers in the model,
            ! so in a moment need to populate the od_scaling array
            do_fill_od_scaling = .true.
          end if ! (jlev <= iend)

          if (do_fill_od_scaling) then

            rand_inhom1 = uniform_distribution_acc(istate)
            ! loop through the sequence of cloudy layers
            !$acc loop seq private(rand_inhom2, overlap_param_inhom, wcdf, icdf, &
            !$acc   wfsd, ifsd)
            do jcloud = max(2,jlev-n_layers_to_scale),jlev-1

              rand_inhom2 = uniform_distribution_acc(istate)

              ! set overlap parameter of inhomogeneities
              overlap_param_inhom = overlap_param(jcloud-1)
              if ( ibegin<=jcloud-1 .and. jcloud-1< iend .and. &
                & overlap_param(jcloud-1) > 0.0_jprb) then
                overlap_param_inhom = &
                  & overlap_param(jcloud-1)**(1.0_jprb/decorrelation_scaling)
              end if

              ! use second random number, and inhomogeneity overlap
              ! parameter, to decide whether the first random number
              ! should be repeated (corresponding to maximum overlap)
              ! or not (corresponding to random overlap)
              if ( jcloud > jlev-n_layers_to_scale .and. &
                & rand_inhom2 >= overlap_param_inhom) then
                rand_inhom1 = uniform_distribution_acc(istate)
              end if
              ! manual inline >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
              ! bilinear interpolation with bounds
              wcdf = rand_inhom1 * (sample_ncdf-1) + 1.0_jprb
              icdf = max(1, min(int(wcdf), sample_ncdf-1))
              wcdf = max(0.0_jprb, min(wcdf - icdf, 1.0_jprb))

              wfsd = (fractional_std(jcloud)-sample_fsd1) * sample_inv_fsd_interval + 1.0_jprb
              ifsd = max(1, min(int(wfsd), sample_nfsd-1))
              wfsd = max(0.0_jprb, min(wfsd - ifsd, 1.0_jprb))

              od_scaling(jg,jcloud) =    (1.0_jprb-wcdf)*(1.0_jprb-wfsd) * sample_val(icdf  ,ifsd)   &
                  & + (1.0_jprb-wcdf)*          wfsd  * sample_val(icdf  ,ifsd+1) &
                  & +           wcdf *(1.0_jprb-wfsd) * sample_val(icdf+1,ifsd)   &
                  & +           wcdf *          wfsd  * sample_val(icdf+1,ifsd+1)
              ! manual inline <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            end do

            n_layers_to_scale = 0

          end if ! (do_fill_od_scaling)
        end do ! jlev = itrigger+1,iend+1
      end do ! jg = 1,ng
    end if ! (total_cloud_cover < frac_threshold)

  end subroutine cloud_generator_acc

end module radiation_cloud_generator_acc
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

