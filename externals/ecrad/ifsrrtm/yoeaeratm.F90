! # 1 "ifsrrtm/yoeaeratm.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/yoeaeratm.f90"
module yoeaeratm

use parkind1  ,only : jpim     ,jprb

implicit none

save

!     ------------------------------------------------------------------
!*    ** *yoeaeratm* - control parameters for aerosols in the atmosphere
!     ------------------------------------------------------------------

integer(kind=jpim) :: naerconf
integer(kind=jpim) :: niniday   
integer(kind=jpim) :: nxt3daer
integer(kind=jpim) :: ndd1, nss1
integer(kind=jpim) :: nbcoptp, nddoptp, nomoptp, nssoptp, nsuoptp
integer(kind=jpim) :: nviswl
integer(kind=jpim) :: ndrydep
integer(kind=jpim) :: nwhtddep, ntddep, nindddep(15)
integer(kind=jpim) :: nwhtscav, ntscav, nindscav(15)
integer(kind=jpim) :: nwhtsedm, ntsedm, nindsedm(15)
integer(kind=jpim) :: nwhtwdep, ntwdep, nindwdep(15)

real(kind=jprb) :: rgrate
real(kind=jprb) :: rmasse(15)

real(kind=jprb) :: repscaer

logical :: laerclimg, laerclimz, laerclist, laerdrydp, laerhygro, laerlisi 
logical :: laerngat , laerscav , laersedim, laersurf , laerelvs , laer6sdia
logical :: laergtop , laerrad  , laerccn  , laeropt(8),laerinit , laervol
logical :: laercstr , laerdiag1, laerdiag2, laerrrtm , laeruvp  , luvindx   
logical :: laerextr , laergbud , laerprnt , laercalip

!     ------------------------------------------------------------------
! ndd1       : location of first bin for desert dust
! nss1       : location of first bin for sea-salt
! nbcoptp    : index for choosing the black carbon lw and sw optic.prop. (1: boucher, 2: bond, bergstrom, 3: stier et al.)
! nddoptp    : index for choosing the dust lw and sw optic.prop. (1: boucher, 2: highwood, 3: woodward)
! nomoptp    : index for choosing the organic carbon optic.prop.
! nssoptp    : index for choosing the sea salt optic.prop.
! nsuoptp    : index for choosing the sulphate optic.prop.
! nviswl     : index of wavelength for visibility computations
! rmfmin     : minimum mass flux for convective aerosol transport
! rgrate     : transformation rate from hygrophopic to hygrophilic for bc and om aerosols
! rmasse     : molar mass: n.b.: either g/mol or avogadro number

! repscaer   : security on aerosol concentration: always >= 1.e-15

! laerclimg  : .t. to start prognostic aerosols with geographical monthly 
!                  mean climatology
! laerclimz  : .t. to start prognostic aerosols with zonal annual mean 
!                  climatology
! laerclist  : .t. to start prognostic aerosols with geographical monthly 
!                  mean climatology for background stratospheric only
! laerdrydp  : .t. dry deposition is active
! laerhydro  : .t. hygroscopic effects on bc and om aerosols
! laerngat   : .t. prevents negative aerosol concentrations
! laerscav   : .t. in-cloud and below cloud scavenging is active
! laersedim  : .t. sedimentation is active
! laersurf   : .t. if surface emissions
! laerelvs   : .t. if "elevated" source
! laer6sdia  : .t. if radiance diagnostics with 6s
! laergtop   : .t. if gas-to-particle conversion for so2/so4
! laerrad    : .t. if there is any prognostic aerosols used for rt
! laeropt(.) : .t. if a given aerosol type is radiatively interactive
! laerccn    : .t. if prognostic aerosols are used to define the re of liq.wat.clds
! laeruvp    : .t. if prognostic aerosols are used in uv-processor
! laercstr   : .t. if climatological stratospheric aerosols are used in radiation 
!                  schemes with the prognostic tropospheric aerosols.
! laerrrtm   : .t. if rrtm schemes get the information from the prognostic aerosols
! laerinit   : .t. if analysed prognostic aerosols are only used as "climatological" aerosols
! laervol    : .t. if volcanic aerosol is considered
! ndrydep    : dry deposition 1: a la gems; 2: "exp" (a la d_grg_4.6)
! nwhtddep   : =0 only ss+du;  =1 ss+du+vol;  =2 ss+du+om+bc+su+vol
! ntddep     : total number of aerosol species to be in- and below-cloud scavenged
! nindddep   : actual set of aerosol indices to be dry deposited
! nwhtsedm   : =0 only ss+du;  =1 ss+du+vol;  =2 ss+du+om+bc+su+vol
! ntsedm     : total number of aerosol species to be in- and below-cloud scavenged
! nindsedm   : actual set of aerosol indices to be sedimented
! nwhtwdep   : =0 only ss+du;  =1 ss+du+vol;  =2 ss+du+om+bc+su+vol
! ntwdep     : total number of aerosol species to be in- and below-cloud scavenged
! nindwdep   : actual set of aerosol indices to be scavenged
! nwhtscav   : =0 only ss+du;  =1 ss+du+vol;  =2 ss+du+om+bc+su+vol
! ntscav     : total number of aerosol species to be in- and below-cloud scavenged
! nindscav   : actual set of aerosol indices to be scavenged
! laercalip  : .t. works with laerlisi=t, store only the caliop-type profile at 532 nm
!     ------------------------------------------------------------------
end module yoeaeratm

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

