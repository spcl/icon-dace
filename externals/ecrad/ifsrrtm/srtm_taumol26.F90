! # 1 "ifsrrtm/srtm_taumol26.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/srtm_taumol26.f90"
! this file has been modified for the use in icon

subroutine srtm_taumol26 &
 & ( kidia   , kfdia    , klev,&
 & p_colmol  ,k_laytrop,&
 & p_sfluxzen, p_taug   , p_taur    , prmu0   &
 & )  

!     written by eli j. mlawer, atmospheric & environmental research.

!     band 26:  22650-29000 cm-1 (low - nothing; high - nothing)

!      parameter (mg=16, mxlay=203, nbands=14)

! modifications
!        m.hamrud      01-oct-2003 cy28 cleaning

!     jjmorcrette 2003-02-24 adapted to ecmwf environment
!        d.salmond  31-oct-2007 vector version in the style of rrtm from meteo france & nec
!     jjmorcrette 20110610 flexible configuration for number of g-points

use parkind1 , only : jpim, jprb
use ecradhook  , only : lhook, dr_hook
use parsrtm  , only : jpg
use yoesrtm  , only : ng26
use yoesrta26, only : sfluxrefc, raylc

implicit none

!-- output
integer(kind=jpim),intent(in)    :: kidia, kfdia 
integer(kind=jpim),intent(in)    :: klev 
real(kind=jprb)   ,intent(in)    :: p_colmol(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: k_laytrop(kidia:kfdia) 

real(kind=jprb)   ,intent(inout) :: p_sfluxzen(kidia:kfdia,jpg) 
real(kind=jprb)   ,intent(inout) :: p_taug(kidia:kfdia,klev,jpg) 
real(kind=jprb)   ,intent(inout) :: p_taur(kidia:kfdia,klev,jpg) 
real(kind=jprb)   ,intent(in)    :: prmu0(kidia:kfdia)
!- from aer
!- from intfac      
!- from intind
!- from precise             
!- from profdata             
!- from self             
integer(kind=jpim) :: ig, i_lay, i_laysolfr(kidia:kfdia), i_nlayers, iplon
integer(kind=jpim) :: laytrop_min, laytrop_max
real(kind=jprb) :: zhook_handle

    !$acc data create(i_laysolfr) &
    !$acc     present(p_colmol, k_laytrop, p_sfluxzen, p_taug, p_taur, prmu0)


    laytrop_min = minval(k_laytrop(kidia:kfdia))
    laytrop_max = maxval(k_laytrop(kidia:kfdia))
! # 67 "ifsrrtm/srtm_taumol26.f90"

    i_nlayers = klev
    !$acc parallel default(none) async(1)
    !$acc loop gang vector
    do iplon = kidia, kfdia
      i_laysolfr(iplon) = k_laytrop(iplon)
    enddo
    !$acc end parallel

    !$acc wait
    !$acc parallel default(none) async(1)
    !$acc loop gang vector collapse(3)
    do i_lay = 1, laytrop_min
       do iplon = kidia, kfdia
!$nec unroll(ng26)
         do ig = 1 , ng26
           if(i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig)=sfluxrefc(ig)
           p_taug(iplon,i_lay,ig) = 0.0_jprb
           p_taur(iplon,i_lay,ig) = p_colmol(iplon,i_lay) * raylc(ig)
         enddo
       enddo
    enddo
    !$acc end parallel

    !$acc parallel default(none) async(1)
    !$acc loop gang vector collapse(2)
    do i_lay = laytrop_min+1, laytrop_max
       do iplon = kidia, kfdia
          if (i_lay <= k_laytrop(iplon)) then
!$nec unroll(ng26)
            !$acc loop seq
            do ig = 1 , ng26
              if(i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig)=sfluxrefc(ig)
              p_taug(iplon,i_lay,ig) = 0.0_jprb
              p_taur(iplon,i_lay,ig) = p_colmol(iplon,i_lay) * raylc(ig)
            enddo
          else
!$nec unroll(ng26)
            !$acc loop seq
            do ig = 1 , ng26
              p_taug(iplon,i_lay,ig) = 0.0_jprb
              p_taur(iplon,i_lay,ig) = p_colmol(iplon,i_lay) * raylc(ig)
            enddo
          endif
       enddo
    enddo
    !$acc end parallel

    !$acc parallel default(none) async(1)
    !$acc loop gang vector collapse(3)
    do ig = 1 , ng26
       do i_lay = laytrop_max+1, i_nlayers
         do iplon = kidia, kfdia
           p_taug(iplon,i_lay,ig) = 0.0_jprb
           p_taur(iplon,i_lay,ig) = p_colmol(iplon,i_lay) * raylc(ig)
         enddo
       enddo
    enddo
    !$acc end parallel

    !$acc wait
    !$acc end data

end subroutine srtm_taumol26
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

