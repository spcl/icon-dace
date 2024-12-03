! # 1 "radiation/radiation_cloud_optics_data.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_cloud_optics_data.f90"
! radiation_cloud_optics_data.f90 - type to store cloud optical properties
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
! # 17 "radiation/radiation_cloud_optics_data.f90" 2

module radiation_cloud_optics_data

  use parkind1, only : jprb

  implicit none
  public

  !---------------------------------------------------------------------
  ! this type holds the configuration information to compute
  ! cloud optical properties
  type cloud_optics_type
     ! band-specific coefficients are provided separately in the
     ! shortwave and longwave, and are dimensioned (nband,ncoeff),
     ! where ncoeff depends on the nature of the parameterization
     real(jprb), allocatable, dimension(:,:) :: &
          &  liq_coeff_lw, liq_coeff_sw, &
          &  ice_coeff_lw, ice_coeff_sw
     ! general coefficients are vectors of length ncoeffgen, which
     ! depends on the nature of the parameterization; note that most
     ! parameterizations use only band-specific coefficients
     real(jprb), allocatable, dimension(:) :: &
          &  liq_coeff_gen, ice_coeff_gen

   contains
     procedure :: setup => setup_cloud_optics

  end type cloud_optics_type

contains

  !---------------------------------------------------------------------
  ! setup cloud optics coefficients by reading them from a file
  subroutine setup_cloud_optics(this, liq_file_name, ice_file_name, iverbose)
    
    use ecradhook,              only : lhook, dr_hook, jphook



    use easy_netcdf,          only : netcdf_file


    class(cloud_optics_type), intent(inout) :: this
    character(len=*), intent(in)            :: liq_file_name, ice_file_name
    integer, intent(in), optional           :: iverbose

    ! the netcdf file containing the coefficients
    type(netcdf_file)  :: file
    integer            :: iverb
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_cloud_optics_data:setup',0,hook_handle)

    if (present(iverbose)) then
      iverb = iverbose
    else
      iverb = 2
    end if

    ! open the droplet scattering file and configure the way it is
    ! read
    call file%open(trim(liq_file_name), iverbose=iverb)
    call file%transpose_matrices()

    ! read the band-specific coefficients
    call file%get('coeff_lw',this%liq_coeff_lw)
    call file%get('coeff_sw',this%liq_coeff_sw)

    ! read the general  coefficients
    if (file%exists('coeff_gen')) then
      call file%get('coeff_gen',this%liq_coeff_gen)
    end if

    ! close droplet scattering file
    call file%close()

    ! open the ice scattering file and configure the way it is read
    call file%open(trim(ice_file_name), iverbose=iverb)
    call file%transpose_matrices()

    ! read the band-specific  coefficients
    call file%get('coeff_lw',this%ice_coeff_lw)
    call file%get('coeff_sw',this%ice_coeff_sw)

    ! read the general  coefficients
    if (file%exists('coeff_gen')) then
      call file%get('coeff_gen',this%ice_coeff_gen)
    end if

    ! close ice scattering file
    call file%close()

    if (lhook) call dr_hook('radiation_cloud_optics_data:setup',1,hook_handle)

  end subroutine setup_cloud_optics

end module radiation_cloud_optics_data
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

