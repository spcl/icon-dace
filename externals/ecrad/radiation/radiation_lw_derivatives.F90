! # 1 "radiation/radiation_lw_derivatives.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_lw_derivatives.f90"
! radiation_lw_derivatives.f90 - compute longwave derivatives for hogan and bozzo (2015) method
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
! this module provides routines to compute the rate of change of
! broadband upwelling longwave flux at each half level with respect to
! the surface broadband upwelling flux.  this is done from the surface
! spectral fluxes and the spectral transmittance of each atmospheric
! layer, assuming no longwave scattering. the result may be used to
! perform approximate updates to the longwave flux profile in between
! calls to the full radiation scheme, accounting for the change in
! skin temperature, following the method of hogan and bozzo (james
! 2015).  separate routines are provided for each solver.
!
! note that currently a more approximate calculation is performed from
! the exact one in hogan and bozzo (2015); here we assume that a
! change in temperature increases the spectral fluxes in proportion,
! when in reality there is a change in shape of the planck function in
! addition to an overall increase in the total emission.
!
! modifications
!   2017-10-23  r. hogan  renamed single-character variables
!   2022-11-22  p. ukkonen / r. hogan  optimized calc_lw_derivatives_region

module radiation_lw_derivatives

  public

contains

  !---------------------------------------------------------------------
  ! calculation for the independent column approximation
  subroutine calc_lw_derivatives_ica(ng, nlev, icol, transmittance, flux_up_surf, lw_derivatives)

    use parkind1, only           : jprb
    use ecradhook,  only           : lhook, dr_hook, jphook

    implicit none

    ! inputs
    integer,    intent(in) :: ng   ! number of spectral intervals
    integer,    intent(in) :: nlev ! number of levels
    integer,    intent(in) :: icol ! index of column for output
    real(jprb), intent(in) :: transmittance(ng,nlev)
    real(jprb), intent(in) :: flux_up_surf(ng) ! upwelling surface spectral flux (w m-2)
    
    ! output
    real(jprb), intent(out) :: lw_derivatives(:,:) ! dimensioned (ncol,nlev+1)

    ! rate of change of spectral flux at a given height with respect
    ! to the surface value
    real(jprb) :: lw_derivatives_g(ng)

    integer    :: jlev

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_lw_derivatives:calc_lw_derivatives_ica',0,hook_handle)

    ! initialize the derivatives at the surface
    lw_derivatives_g = flux_up_surf / sum(flux_up_surf)
    lw_derivatives(icol, nlev+1) = 1.0_jprb

    ! move up through the atmosphere computing the derivatives at each
    ! half-level
    do jlev = nlev,1,-1
      lw_derivatives_g = lw_derivatives_g * transmittance(:,jlev)
      lw_derivatives(icol,jlev) = sum(lw_derivatives_g)
    end do

    if (lhook) call dr_hook('radiation_lw_derivatives:calc_lw_derivatives_ica',1,hook_handle)

  end subroutine calc_lw_derivatives_ica


  !---------------------------------------------------------------------
  ! calculation for the independent column approximation
  subroutine modify_lw_derivatives_ica(ng, nlev, icol, transmittance, &
       &                               flux_up_surf, weight, lw_derivatives)

    use parkind1, only           : jprb
    use ecradhook,  only           : lhook, dr_hook, jphook

    implicit none

    ! inputs
    integer,    intent(in) :: ng   ! number of spectral intervals
    integer,    intent(in) :: nlev ! number of levels
    integer,    intent(in) :: icol ! index of column for output
    real(jprb), intent(in) :: transmittance(ng,nlev)
    real(jprb), intent(in) :: flux_up_surf(ng) ! upwelling surface spectral flux (w m-2)
    real(jprb), intent(in) :: weight ! weight new values against existing
    
    ! output
    real(jprb), intent(inout) :: lw_derivatives(:,:) ! dimensioned (ncol,nlev+1)

    ! rate of change of spectral flux at a given height with respect
    ! to the surface value
    real(jprb) :: lw_derivatives_g(ng)

    integer    :: jlev

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_lw_derivatives:modify_lw_derivatives_ica',0,hook_handle)

    ! initialize the derivatives at the surface
    lw_derivatives_g = flux_up_surf / sum(flux_up_surf)
    ! this value must be 1 so no weighting applied
    lw_derivatives(icol, nlev+1) = 1.0_jprb

    ! move up through the atmosphere computing the derivatives at each
    ! half-level
    do jlev = nlev,1,-1
      lw_derivatives_g = lw_derivatives_g * transmittance(:,jlev)
      lw_derivatives(icol,jlev) = (1.0_jprb - weight) * lw_derivatives(icol,jlev) &
           &                    + weight * sum(lw_derivatives_g)
    end do

    if (lhook) call dr_hook('radiation_lw_derivatives:modify_lw_derivatives_ica',1,hook_handle)

  end subroutine modify_lw_derivatives_ica



  !---------------------------------------------------------------------
  ! calculation for solvers involving multiple regions and matrices
  subroutine calc_lw_derivatives_matrix(ng, nlev, nreg, icol, transmittance, &
       &                                u_matrix, flux_up_surf, lw_derivatives)

    use parkind1, only           : jprb
    use ecradhook,  only           : lhook, dr_hook, jphook

    use radiation_matrix

    implicit none

    ! inputs
    integer,    intent(in) :: ng   ! number of spectral intervals
    integer,    intent(in) :: nlev ! number of levels
    integer,    intent(in) :: nreg ! number of regions
    integer,    intent(in) :: icol ! index of column for output
    real(jprb), intent(in) :: transmittance(ng,nreg,nreg,nlev)
    real(jprb), intent(in) :: u_matrix(nreg,nreg,nlev+1) ! upward overlap matrix
    real(jprb), intent(in) :: flux_up_surf(ng) ! upwelling surface spectral flux (w m-2)
    
    ! output
    real(jprb), intent(out) :: lw_derivatives(:,:) ! dimensioned (ncol,nlev+1)

    ! rate of change of spectral flux at a given height with respect
    ! to the surface value
    real(jprb) :: lw_derivatives_g_reg(ng,nreg)

    integer    :: jlev

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_lw_derivatives:calc_lw_derivatives_matrix',0,hook_handle)

    ! initialize the derivatives at the surface; the surface is
    ! treated as a single clear-sky layer so we only need to put
    ! values in region 1.
    lw_derivatives_g_reg = 0.0_jprb
    lw_derivatives_g_reg(:,1) = flux_up_surf / sum(flux_up_surf)
    lw_derivatives(icol, nlev+1) = 1.0_jprb

    ! move up through the atmosphere computing the derivatives at each
    ! half-level
    do jlev = nlev,1,-1
      ! compute effect of overlap at half-level jlev+1, yielding
      ! derivatives just above that half-level
      lw_derivatives_g_reg = singlemat_x_vec(ng,ng,nreg,u_matrix(:,:,jlev+1),lw_derivatives_g_reg)

      ! compute effect of transmittance of layer jlev, yielding
      ! derivatives just below the half-level above (jlev)
      lw_derivatives_g_reg = mat_x_vec(ng,ng,nreg,transmittance(:,:,:,jlev),lw_derivatives_g_reg)

      lw_derivatives(icol, jlev) = sum(lw_derivatives_g_reg)
    end do

    if (lhook) call dr_hook('radiation_lw_derivatives:calc_lw_derivatives_matrix',1,hook_handle)

  end subroutine calc_lw_derivatives_matrix


  !---------------------------------------------------------------------
  ! calculation for solvers involving multiple regions but no 3d
  ! effects: the difference from calc_lw_derivatives_matrix is that transmittance
  ! has one fewer dimensions
  subroutine calc_lw_derivatives_region(ng, nlev, nreg, icol, transmittance, &
       &                                u_matrix, flux_up_surf, lw_derivatives)

    use parkind1, only           : jprb
    use ecradhook,  only           : lhook, dr_hook, jphook

    use radiation_matrix

    implicit none

    ! inputs
    integer,    intent(in) :: ng   ! number of spectral intervals
    integer,    intent(in) :: nlev ! number of levels
    integer,    intent(in) :: nreg ! number of regions
    integer,    intent(in) :: icol ! index of column for output
    real(jprb), intent(in) :: transmittance(ng,nreg,nlev)
    real(jprb), intent(in) :: u_matrix(nreg,nreg,nlev+1) ! upward overlap matrix
    real(jprb), intent(in) :: flux_up_surf(ng) ! upwelling surface spectral flux (w m-2)
    
    ! output
    real(jprb), intent(out) :: lw_derivatives(:,:) ! dimensioned (ncol,nlev+1)

    ! rate of change of spectral flux at a given height with respect
    ! to the surface value
    real(jprb) :: lw_deriv(ng,nreg), lw_deriv_below(ng,nreg)
    real(jprb) :: partial_sum(ng)

    integer    :: jlev, jg

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_lw_derivatives:calc_lw_derivatives_region',0,hook_handle)

    ! initialize the derivatives at the surface; the surface is
    ! treated as a single clear-sky layer so we only need to put
    ! values in region 1.
    lw_deriv = 0.0_jprb
    lw_deriv(:,1) = flux_up_surf / sum(flux_up_surf)
    lw_derivatives(icol, nlev+1) = 1.0_jprb

    if (nreg == 3) then 
      ! optimize the most common case of 3 regions by removing the
      ! nested call to singlemat_x_vec and unrolling the matrix
      ! multiplication inline
      
      do jlev = nlev,1,-1
        ! compute effect of overlap at half-level jlev+1, yielding
        ! derivatives just above that half-level
        lw_deriv_below = lw_deriv
        
        associate(a=>u_matrix(:,:,jlev+1), b=>lw_deriv_below)
          do jg = 1,ng   
            ! both inner and outer loop of the matrix loops j1 and j2 unrolled
            ! inner loop:        j2=1             j2=2             j2=3 
            lw_deriv(jg,1) = a(1,1)*b(jg,1) + a(1,2)*b(jg,2) + a(1,3)*b(jg,3) 
            lw_deriv(jg,2) = a(2,1)*b(jg,1) + a(2,2)*b(jg,2) + a(2,3)*b(jg,3) 
            lw_deriv(jg,3) = a(3,1)*b(jg,1) + a(3,2)*b(jg,2) + a(3,3)*b(jg,3) 

            ! compute effect of transmittance of layer jlev, yielding
            ! derivatives just below the half-level above (jlev)
            lw_deriv(jg,1) = lw_deriv(jg,1) * transmittance(jg,1,jlev)
            lw_deriv(jg,2) = lw_deriv(jg,2) * transmittance(jg,2,jlev)
            lw_deriv(jg,3) = lw_deriv(jg,3) * transmittance(jg,3,jlev)

            partial_sum(jg) = lw_deriv(jg,1) + lw_deriv(jg,2) + lw_deriv(jg,3)
          end do
        end associate

        lw_derivatives(icol, jlev) = sum(partial_sum)
      end do
    else
      ! general case when number of regions is not 3
      
      ! move up through the atmosphere computing the derivatives at each
      ! half-level
      do jlev = nlev,1,-1
        ! compute effect of overlap at half-level jlev+1, yielding
        ! derivatives just above that half-level
        lw_deriv = singlemat_x_vec(ng,ng,nreg,u_matrix(:,:,jlev+1),lw_deriv)
        
        ! compute effect of transmittance of layer jlev, yielding
        ! derivatives just below the half-level above (jlev)
        lw_deriv = transmittance(:,:,jlev) * lw_deriv
        
        lw_derivatives(icol, jlev) = sum(lw_deriv)
      end do
    end if
    
    if (lhook) call dr_hook('radiation_lw_derivatives:calc_lw_derivatives_region',1,hook_handle)

  end subroutine calc_lw_derivatives_region


end module radiation_lw_derivatives
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

