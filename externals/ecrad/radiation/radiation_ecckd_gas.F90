! # 1 "radiation/radiation_ecckd_gas.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_ecckd_gas.f90"
! radiation_ecckd_gas.f90 - type representing a single ecckd gas
!
! (c) copyright 2020- ecmwf.
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
! license: see the copying file for details
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
! # 18 "radiation/radiation_ecckd_gas.f90" 2

module radiation_ecckd_gas

  use parkind1, only : jprb
  use radiation_gas_constants

  implicit none

  public

  ! concentration dependence of individual gases
  enum, bind(c)
    enumerator :: iconcdependencenone = 0, &
         &        iconcdependencelinear, &
         &        iconcdependencelut, &
         &        iconcdependencerelativelinear
  end enum

  !---------------------------------------------------------------------
  ! this derived type describes a correlated k-distribution
  ! representation of an individual gas (including composite gases)
  type ckd_gas_type

    ! code identifying the gas, from the codes in the
    ! radiation_gas_constants module
    integer :: i_gas_code = -1

    ! one of the iconcdependence* enumerators
    integer :: i_conc_dependence

    ! molar absorption coefficient in m2 mol-1. if
    ! i_conc_dependence==iconcdependencenone then it is the absorption
    ! cross section per mole of dry air.  if
    ! conc_dependence==iconcdependencelinear|iconcdependencerelativelinear,
    ! it is the absorption cross section per mole of the gas in
    ! question. it is dimensioned (g_point,pressure,temperature).
    real(jprb), allocatable :: molar_abs(:,:,:)
    
    ! if i_conc_dependence==iconcdependencelut then we have an
    ! additional dimension for concentration. it is dimensioned
    ! (g_point,pressure,temperature,conc)
    real(jprb), allocatable :: molar_abs_conc(:,:,:,:)

    ! if i_conc_dependence==iconcdependencerelativelinear then the
    ! following reference concentration is subtracted from the actual
    ! concentration before the result is multiplied by the mass
    ! absorption coefficient
    real(jprb) :: reference_mole_frac = 0.0_jprb

    ! mole fraction coordinate variable if
    ! i_conc_dependence==iconcdependencelut
    real(jprb) :: log_mole_frac1 = 0.0_jprb, d_log_mole_frac = 1.0_jprb
    integer    :: n_mole_frac = 0

  contains

    procedure :: read => read_ckd_gas
!    procedure :: deallocate => deallocate_ckd_gas

  end type ckd_gas_type

contains

  !---------------------------------------------------------------------
  ! read information about the representation of a single gas from a
  ! netcdf file, identifying it with code i_gas_code
  subroutine read_ckd_gas(this, file, gas_name, i_gas_code)




    use easy_netcdf,          only : netcdf_file


    class(ckd_gas_type), intent(inout) :: this
    type(netcdf_file),   intent(inout) :: file
    character(len=*),    intent(in)    :: gas_name
    integer,             intent(in)    :: i_gas_code
    
    ! local storage for mole fraction coordinate variable
    real(jprb), allocatable :: mole_fraction(:)

    this%i_gas_code = i_gas_code

    call file%get(gas_name // "_conc_dependence_code", this%i_conc_dependence)
    if (this%i_conc_dependence == iconcdependencelut) then
      call file%get(gas_name // "_molar_absorption_coeff", &
           &        this%molar_abs_conc)
      call file%get(gas_name // "_mole_fraction", mole_fraction)
      this%log_mole_frac1  = log(mole_fraction(1))
      this%n_mole_frac     = size(mole_fraction)
      this%d_log_mole_frac = (log(mole_fraction(size(mole_fraction))) &
           &                  - this%log_mole_frac1) / (this%n_mole_frac-1)
      deallocate(mole_fraction)
    else
      call file%get(gas_name // "_molar_absorption_coeff", &
           &        this%molar_abs)
    end if

    if (this%i_conc_dependence == iconcdependencerelativelinear) then
      call file%get(gas_name // "_reference_mole_fraction", &
           &        this%reference_mole_frac)
    end if

  end subroutine read_ckd_gas

end module radiation_ecckd_gas
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

