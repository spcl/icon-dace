! # 1 "radiation/radiation_adding_ica_lw.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_adding_ica_lw.f90"
! this file has been modified for the use in icon

! radiation_adding_ica_lw.f90 - longwave adding method in independent column approximation
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
!   2017-07-12  r. hogan  fast adding method for if only clouds scatter
!   2017-10-23  r. hogan  renamed single-character variables

module radiation_adding_ica_lw

  public

contains

  !---------------------------------------------------------------------
  ! use the scalar "adding" method to compute longwave flux profiles,
  ! including scattering, by successively adding the contribution of
  ! layers starting from the surface to compute the total albedo and
  ! total upward emission of the increasingly larger block of
  ! atmospheric layers.
  subroutine adding_ica_lw(ncol, nlev, &
       &  reflectance, transmittance, source_up, source_dn, emission_surf, albedo_surf, &
       &  flux_up, flux_dn)

    use parkind1, only           : jprb
    use ecradhook,  only           : lhook, dr_hook, jphook

    implicit none

    ! inputs
    integer, intent(in) :: ncol ! number of columns (may be spectral intervals)
    integer, intent(in) :: nlev ! number of levels

    ! surface emission (w m-2) and albedo
    real(jprb), intent(in),  dimension(ncol) :: emission_surf, albedo_surf

    ! diffuse reflectance and transmittance of each layer
    real(jprb), intent(in),  dimension(ncol, nlev)   :: reflectance, transmittance

    ! emission from each layer in an upward and downward direction
    real(jprb), intent(in),  dimension(ncol, nlev)   :: source_up, source_dn

    ! resulting fluxes (w m-2) at half-levels: diffuse upwelling and
    ! downwelling
    real(jprb), intent(out), dimension(ncol, nlev+1) :: flux_up, flux_dn
    
    ! albedo of the entire earth/atmosphere system below each half
    ! level
    real(jprb), dimension(ncol, nlev+1) :: albedo

    ! upwelling radiation at each half-level due to emission below
    ! that half-level (w m-2)
    real(jprb), dimension(ncol, nlev+1) :: source

    ! equal to 1/(1-albedo*reflectance)
    real(jprb), dimension(ncol, nlev)   :: inv_denominator

    ! loop index for model level and column
    integer :: jlev, jcol

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_adding_ica_lw:adding_ica_lw',0,hook_handle)

    albedo(:,nlev+1) = albedo_surf

    ! at the surface, the source is thermal emission
    source(:,nlev+1) = emission_surf

    ! work back up through the atmosphere and compute the albedo of
    ! the entire earth/atmosphere system below that half-level, and
    ! also the "source", which is the upwelling flux due to emission
    ! below that level
    do jlev = nlev,1,-1
      ! next loop over columns. we could do this by indexing the
      ! entire inner dimension as follows, e.g. for the first line:
      !   inv_denominator(:,jlev) = 1.0_jprb / (1.0_jprb-albedo(:,jlev+1)*reflectance(:,jlev))
      ! and similarly for subsequent lines, but this slows down the
      ! routine by a factor of 2!  rather, we do it with an explicit
      ! loop.
      do jcol = 1,ncol
        ! lacis and hansen (1974) eq 33, shonk & hogan (2008) eq 10:
        inv_denominator(jcol,jlev) = 1.0_jprb &
             &  / (1.0_jprb-albedo(jcol,jlev+1)*reflectance(jcol,jlev))
        ! shonk & hogan (2008) eq 9, petty (2006) eq 13.81:
        albedo(jcol,jlev) = reflectance(jcol,jlev) + transmittance(jcol,jlev)*transmittance(jcol,jlev) &
             &  * albedo(jcol,jlev+1) * inv_denominator(jcol,jlev)
        ! shonk & hogan (2008) eq 11:
        source(jcol,jlev) = source_up(jcol,jlev) &
             &  + transmittance(jcol,jlev) * (source(jcol,jlev+1) &
             &                    + albedo(jcol,jlev+1)*source_dn(jcol,jlev)) &
             &                   * inv_denominator(jcol,jlev)
      end do
    end do

    ! at top-of-atmosphere there is no diffuse downwelling radiation
    flux_dn(:,1) = 0.0_jprb

    ! at top-of-atmosphere, all upwelling radiation is due to emission
    ! below that level
    flux_up(:,1) = source(:,1)

    ! work back down through the atmosphere computing the fluxes at
    ! each half-level
    do jlev = 1,nlev
      do jcol = 1,ncol
        ! shonk & hogan (2008) eq 14 (after simplification):
        flux_dn(jcol,jlev+1) &
             &  = (transmittance(jcol,jlev)*flux_dn(jcol,jlev) &
             &     + reflectance(jcol,jlev)*source(jcol,jlev+1) &
             &     + source_dn(jcol,jlev)) * inv_denominator(jcol,jlev)
        ! shonk & hogan (2008) eq 12:
        flux_up(jcol,jlev+1) = albedo(jcol,jlev+1)*flux_dn(jcol,jlev+1) &
             &            + source(jcol,jlev+1)
      end do
    end do

    if (lhook) call dr_hook('radiation_adding_ica_lw:adding_ica_lw',1,hook_handle)

  end subroutine adding_ica_lw


  !---------------------------------------------------------------------
  ! use the scalar "adding" method to compute longwave flux profiles,
  ! including scattering in cloudy layers only.
  subroutine fast_adding_ica_lw(ncol, nlev, &
       &  reflectance, transmittance, source_up, source_dn, emission_surf, albedo_surf, &
       &  is_clear_sky_layer, i_cloud_top, flux_dn_clear, &
       &  flux_up, flux_dn, albedo, source, inv_denominator)

    use parkind1, only           : jprb
    use ecradhook,  only           : lhook, dr_hook, jphook

    implicit none

    ! inputs
    integer, intent(in) :: ncol ! number of columns (may be spectral intervals)
    integer, intent(in) :: nlev ! number of levels

    ! surface emission (w m-2) and albedo
    real(jprb), intent(in),  dimension(ncol) :: emission_surf, albedo_surf

    ! diffuse reflectance and transmittance of each layer
    real(jprb), intent(in),  dimension(ncol, nlev)   :: reflectance, transmittance

    ! emission from each layer in an upward and downward direction
    real(jprb), intent(in),  dimension(ncol, nlev)   :: source_up, source_dn

    ! determine which layers are cloud-free
    logical, intent(in) :: is_clear_sky_layer(nlev)

    ! index to highest cloudy layer
    integer, intent(in) :: i_cloud_top

    ! pre-computed clear-sky downwelling fluxes (w m-2) at half-levels
    real(jprb), intent(in), dimension(ncol, nlev+1)  :: flux_dn_clear

    ! resulting fluxes (w m-2) at half-levels: diffuse upwelling and
    ! downwelling
    real(jprb), intent(out), dimension(ncol, nlev+1) :: flux_up, flux_dn
    
    ! albedo of the entire earth/atmosphere system below each half
    ! level
    real(jprb), intent(out), dimension(ncol, nlev+1) :: albedo

    ! upwelling radiation at each half-level due to emission below
    ! that half-level (w m-2)
    real(jprb), intent(out), dimension(ncol, nlev+1) :: source

    ! equal to 1/(1-albedo*reflectance)
    real(jprb), intent(out), dimension(ncol, nlev)   :: inv_denominator

    ! loop index for model level and column
    integer :: jlev, jcol

    real(jphook) :: hook_handle

    !$acc routine worker 


    if (lhook) call dr_hook('radiation_adding_ica_lw:fast_adding_ica_lw',0,hook_handle)


    ! copy over downwelling fluxes above cloud from clear sky
    flux_dn(:,1:i_cloud_top) = flux_dn_clear(:,1:i_cloud_top)

    albedo(:,nlev+1) = albedo_surf
    
    ! at the surface, the source is thermal emission
    source(:,nlev+1) = emission_surf

    ! work back up through the atmosphere and compute the albedo of
    ! the entire earth/atmosphere system below that half-level, and
    ! also the "source", which is the upwelling flux due to emission
    ! below that level
    !$acc loop seq
    do jlev = nlev,i_cloud_top,-1
      if (is_clear_sky_layer(jlev)) then
        ! reflectance of this layer is zero, simplifying the expression
        !$acc loop worker vector
        do jcol = 1,ncol
          albedo(jcol,jlev) = transmittance(jcol,jlev)*transmittance(jcol,jlev)*albedo(jcol,jlev+1)
          source(jcol,jlev) = source_up(jcol,jlev) &
               &  + transmittance(jcol,jlev) * (source(jcol,jlev+1) &
               &                    + albedo(jcol,jlev+1)*source_dn(jcol,jlev))
        end do
      else
        ! loop over columns; explicit loop seems to be faster
        !$acc loop worker vector
        do jcol = 1,ncol
          ! lacis and hansen (1974) eq 33, shonk & hogan (2008) eq 10:
          inv_denominator(jcol,jlev) = 1.0_jprb &
               &  / (1.0_jprb-albedo(jcol,jlev+1)*reflectance(jcol,jlev))
          ! shonk & hogan (2008) eq 9, petty (2006) eq 13.81:
          albedo(jcol,jlev) = reflectance(jcol,jlev) + transmittance(jcol,jlev)*transmittance(jcol,jlev) &
               &  * albedo(jcol,jlev+1) * inv_denominator(jcol,jlev)
          ! shonk & hogan (2008) eq 11:
          source(jcol,jlev) = source_up(jcol,jlev) &
               &  + transmittance(jcol,jlev) * (source(jcol,jlev+1) &
               &                    + albedo(jcol,jlev+1)*source_dn(jcol,jlev)) &
               &                   * inv_denominator(jcol,jlev)
        end do
      end if
    end do

    ! compute the fluxes above the highest cloud
    flux_up(:,i_cloud_top) = source(:,i_cloud_top) &
         &                 + albedo(:,i_cloud_top)*flux_dn(:,i_cloud_top)
    !$acc loop seq
    do jlev = i_cloud_top-1,1,-1
      flux_up(:,jlev) = transmittance(:,jlev)*flux_up(:,jlev+1) + source_up(:,jlev)
    end do

    ! work back down through the atmosphere from cloud top computing
    ! the fluxes at each half-level
    !$acc loop seq
    do jlev = i_cloud_top,nlev
      if (is_clear_sky_layer(jlev)) then
        !$acc loop worker vector
        do jcol = 1,ncol
          flux_dn(jcol,jlev+1) = transmittance(jcol,jlev)*flux_dn(jcol,jlev) &
               &               + source_dn(jcol,jlev)
          flux_up(jcol,jlev+1) = albedo(jcol,jlev+1)*flux_dn(jcol,jlev+1) &
               &               + source(jcol,jlev+1)
        end do
      else
        !$acc loop worker vector
        do jcol = 1,ncol
          ! shonk & hogan (2008) eq 14 (after simplification):
          flux_dn(jcol,jlev+1) &
               &  = (transmittance(jcol,jlev)*flux_dn(jcol,jlev) &
               &     + reflectance(jcol,jlev)*source(jcol,jlev+1) &
               &     + source_dn(jcol,jlev)) * inv_denominator(jcol,jlev)
          ! shonk & hogan (2008) eq 12:
          flux_up(jcol,jlev+1) = albedo(jcol,jlev+1)*flux_dn(jcol,jlev+1) &
               &               + source(jcol,jlev+1)
        end do
      end if
    end do


    if (lhook) call dr_hook('radiation_adding_ica_lw:fast_adding_ica_lw',1,hook_handle)


  end subroutine fast_adding_ica_lw


  !---------------------------------------------------------------------
  ! if there is no scattering then fluxes may be computed simply by
  ! passing down through the atmosphere computing the downwelling
  ! fluxes from the transmission and emission of each layer, and then
  ! passing back up through the atmosphere to compute the upwelling
  ! fluxes in the same way.
  subroutine calc_fluxes_no_scattering_lw(ncol, nlev, &
       &  transmittance, source_up, source_dn, emission_surf, albedo_surf, flux_up, flux_dn)

    use parkind1, only           : jprb

    use ecradhook,  only           : lhook, dr_hook, jphook


    implicit none

    ! inputs
    integer, intent(in) :: ncol ! number of columns (may be spectral intervals)
    integer, intent(in) :: nlev ! number of levels

    ! surface emission (w m-2) and albedo
    real(jprb), intent(in),  dimension(ncol) :: emission_surf, albedo_surf

    ! diffuse reflectance and transmittance of each layer
    real(jprb), intent(in),  dimension(ncol, nlev)   :: transmittance

    ! emission from each layer in an upward and downward direction
    real(jprb), intent(in),  dimension(ncol, nlev)   :: source_up, source_dn

    ! resulting fluxes (w m-2) at half-levels: diffuse upwelling and
    ! downwelling
    real(jprb), intent(out), dimension(ncol, nlev+1) :: flux_up, flux_dn
    
    ! loop index for model level
    integer :: jlev, jcol


    real(jphook) :: hook_handle


    !$acc routine worker


    if (lhook) call dr_hook('radiation_adding_ica_lw:calc_fluxes_no_scattering_lw',0,hook_handle)


    ! at top-of-atmosphere there is no diffuse downwelling radiation
    flux_dn(:,1) = 0.0_jprb

    ! work down through the atmosphere computing the downward fluxes
    ! at each half-level
!$acc loop seq
! added for dwd (2020)
!nec$ outerloop_unroll(8)
    do jlev = 1,nlev
      !$acc loop worker vector
      do jcol = 1,ncol
        flux_dn(jcol,jlev+1) = transmittance(jcol,jlev)*flux_dn(jcol,jlev) + source_dn(jcol,jlev)
      end do
    end do

    ! surface reflection and emission
    flux_up(:,nlev+1) = emission_surf + albedo_surf * flux_dn(:,nlev+1)

    ! work back up through the atmosphere computing the upward fluxes
    ! at each half-level
!$acc loop seq
! added for dwd (2020)
!nec$ outerloop_unroll(8)
    do jlev = nlev,1,-1
      !$acc loop worker vector
      do jcol = 1,ncol
        flux_up(jcol,jlev) = transmittance(jcol,jlev)*flux_up(jcol,jlev+1) + source_up(jcol,jlev)
      end do
    end do
    

    if (lhook) call dr_hook('radiation_adding_ica_lw:calc_fluxes_no_scattering_lw',1,hook_handle)


  end subroutine calc_fluxes_no_scattering_lw

end module radiation_adding_ica_lw
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

