! # 1 "radiation/radiation_pdf_sampler.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_pdf_sampler.f90"
! radiation_pdf_sampler.f90 - get samples from a pdf for mcica
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

module radiation_pdf_sampler

  use parkind1, only : jprb

  implicit none
  public

  !---------------------------------------------------------------------
  ! derived type for sampling from a lognormal or gamma distribution,
  ! or other pdf, used to generate water content or optical depth
  ! scalings for use in the monte carlo independent column
  ! approximation (mcica)
  type pdf_sampler_type
    ! number of points in look-up table for cumulative distribution
    ! function (cdf) and fractional standard deviation (fsd)
    ! dimensions
    integer :: ncdf, nfsd

    ! first value of fsd and the reciprocal of the interval between
    ! fsd values (which are assumed to be uniformly distributed)
    real(jprb) :: fsd1, inv_fsd_interval

    ! value of the distribution for each cdf and fsd bin
    real(jprb), allocatable, dimension(:,:) :: val

  contains

    procedure :: setup => setup_pdf_sampler
    procedure :: sample => sample_from_pdf
    procedure :: masked_sample => sample_from_pdf_masked
    procedure :: block_sample => sample_from_pdf_block
    procedure :: masked_block_sample => sample_from_pdf_masked_block
    procedure :: deallocate => deallocate_pdf_sampler

  end type pdf_sampler_type

contains

  !---------------------------------------------------------------------
  ! load look-up table from a file 
  subroutine setup_pdf_sampler(this, file_name, iverbose)
    
    use ecradhook,     only : lhook, dr_hook, jphook
    use easy_netcdf, only : netcdf_file

    class(pdf_sampler_type), intent(inout) :: this
    character(len=*),        intent(in)    :: file_name
    integer, optional,       intent(in)    :: iverbose

    type(netcdf_file)  :: file
    integer            :: iverb
    real(jprb), allocatable :: fsd(:)

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_pdf_sampler:setup',0,hook_handle)

    if (present(iverbose)) then
      iverb = iverbose
    else
      iverb = 2
    end if

    if (allocated(this%val)) then
      deallocate(this%val)
    end if

    call file%open(trim(file_name), iverbose=iverb)

    call file%get('fsd',fsd)
    call file%get('x',  this%val)

    call file%close()

    this%ncdf = size(this%val,1)
    this%nfsd = size(this%val,2)
    this%fsd1 = fsd(1)
    this%inv_fsd_interval = 1.0_jprb / (fsd(2)-fsd(1))

    deallocate(fsd)

    if (lhook) call dr_hook('radiation_pdf_sampler:setup',1,hook_handle)

  end subroutine setup_pdf_sampler

  !---------------------------------------------------------------------
  ! deallocate data in pdf_sampler_type derived type
  subroutine deallocate_pdf_sampler(this)

    use ecradhook,     only : lhook, dr_hook, jphook

    class(pdf_sampler_type), intent(inout) :: this
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_pdf_sampler:deallocate',0,hook_handle)

    if (allocated(this%val)) then
      deallocate(this%val)
    end if

    if (lhook) call dr_hook('radiation_pdf_sampler:deallocate',1,hook_handle)
    
  end subroutine deallocate_pdf_sampler


  !---------------------------------------------------------------------
  ! extract the value from a pdf with fractional standard deviation
  ! "fsd" corresponding to the cumulative distribution function value
  ! "cdf", and return it in val. since this is an elemental
  ! subroutine, fsd, cdf and val may be arrays.
  elemental subroutine sample_from_pdf(this, fsd, cdf, val)
    
    class(pdf_sampler_type), intent(in)  :: this

    ! fractional standard deviation (0 to 4) and cumulative
    ! distribution function (0 to 1)
    real(jprb),              intent(in)  :: fsd, cdf

    ! sample from distribution
    real(jprb),              intent(out) :: val

    ! index to look-up table
    integer    :: ifsd, icdf

    ! weights in bilinear interpolation
    real(jprb) :: wfsd, wcdf

    ! bilinear interpolation with bounds
    wcdf = cdf * (this%ncdf-1) + 1.0_jprb
    icdf = max(1, min(int(wcdf), this%ncdf-1))
    wcdf = max(0.0_jprb, min(wcdf - icdf, 1.0_jprb))

    wfsd = (fsd-this%fsd1) * this%inv_fsd_interval + 1.0_jprb
    ifsd = max(1, min(int(wfsd), this%nfsd-1))
    wfsd = max(0.0_jprb, min(wfsd - ifsd, 1.0_jprb))

    val =    (1.0_jprb-wcdf)*(1.0_jprb-wfsd) * this%val(icdf  ,ifsd)   &
         & + (1.0_jprb-wcdf)*          wfsd  * this%val(icdf  ,ifsd+1) &
         & +           wcdf *(1.0_jprb-wfsd) * this%val(icdf+1,ifsd)   &
         & +           wcdf *          wfsd  * this%val(icdf+1,ifsd+1)

  end subroutine sample_from_pdf


  !---------------------------------------------------------------------
  ! for true elements of mask, extract the values of a pdf with
  ! fractional standard deviation "fsd" corresponding to the
  ! cumulative distribution function values "cdf", and return in
  ! val. for false elements of mask, return zero in val.
  subroutine sample_from_pdf_masked(this, nsamp, fsd, cdf, val, mask)
    
    class(pdf_sampler_type), intent(in)  :: this

    ! number of samples
    integer,    intent(in) :: nsamp

    ! fractional standard deviation (0 to 4) and cumulative
    ! distribution function (0 to 1)
    real(jprb), intent(in)  :: fsd(nsamp), cdf(nsamp)

    ! sample from distribution
    real(jprb), intent(out) :: val(:)

    ! mask
    logical,    intent(in) :: mask(nsamp)

    ! loop index
    integer    :: jsamp

    ! index to look-up table
    integer    :: ifsd, icdf

    ! weights in bilinear interpolation
    real(jprb) :: wfsd, wcdf

    do jsamp = 1,nsamp
      if (mask(jsamp)) then
        ! bilinear interpolation with bounds
        wcdf = cdf(jsamp) * (this%ncdf-1) + 1.0_jprb
        icdf = max(1, min(int(wcdf), this%ncdf-1))
        wcdf = max(0.0_jprb, min(wcdf - icdf, 1.0_jprb))
        
        wfsd = (fsd(jsamp)-this%fsd1) * this%inv_fsd_interval + 1.0_jprb
        ifsd = max(1, min(int(wfsd), this%nfsd-1))
        wfsd = max(0.0_jprb, min(wfsd - ifsd, 1.0_jprb))
        
        val(jsamp)=(1.0_jprb-wcdf)*(1.0_jprb-wfsd) * this%val(icdf  ,ifsd)   &
             &    +(1.0_jprb-wcdf)*          wfsd  * this%val(icdf  ,ifsd+1) &
             &    +          wcdf *(1.0_jprb-wfsd) * this%val(icdf+1,ifsd)   &
             &    +          wcdf *          wfsd  * this%val(icdf+1,ifsd+1)
      else
        val(jsamp) = 0.0_jprb
      end if
    end do
  end subroutine sample_from_pdf_masked

  !---------------------------------------------------------------------
  ! extract the values of a pdf with fractional standard deviation
  ! "fsd" corresponding to the cumulative distribution function values
  ! "cdf", and return in val. this version works on 2d blocks of data.
  subroutine sample_from_pdf_block(this, nz, ng, fsd, cdf, val)
    
    class(pdf_sampler_type), intent(in)  :: this

    ! number of samples
    integer,    intent(in) :: nz, ng

    ! fractional standard deviation (0 to 4) and cumulative
    ! distribution function (0 to 1)
    real(jprb), intent(in)  :: fsd(nz), cdf(ng, nz)

    ! sample from distribution
    real(jprb), intent(out) :: val(:,:)

    ! loop index
    integer    :: jz, jg

    ! index to look-up table
    integer    :: ifsd, icdf

    ! weights in bilinear interpolation
    real(jprb) :: wfsd, wcdf

    do jz = 1,nz
      do jg = 1,ng
        if (cdf(jg, jz) > 0.0_jprb) then
          ! bilinear interpolation with bounds
          wcdf = cdf(jg,jz) * (this%ncdf-1) + 1.0_jprb
          icdf = max(1, min(int(wcdf), this%ncdf-1))
          wcdf = max(0.0_jprb, min(wcdf - icdf, 1.0_jprb))
          
          wfsd = (fsd(jz)-this%fsd1) * this%inv_fsd_interval + 1.0_jprb
          ifsd = max(1, min(int(wfsd), this%nfsd-1))
          wfsd = max(0.0_jprb, min(wfsd - ifsd, 1.0_jprb))
          
          val(jg,jz)=(1.0_jprb-wcdf)*(1.0_jprb-wfsd) * this%val(icdf  ,ifsd)   &
               &    +(1.0_jprb-wcdf)*          wfsd  * this%val(icdf  ,ifsd+1) &
               &    +          wcdf *(1.0_jprb-wfsd) * this%val(icdf+1,ifsd)   &
               &    +          wcdf *          wfsd  * this%val(icdf+1,ifsd+1)
        else
          val(jg,jz) = 0.0_jprb
        end if
      end do
    end do

  end subroutine sample_from_pdf_block

  !---------------------------------------------------------------------
  ! extract the values of a pdf with fractional standard deviation
  ! "fsd" corresponding to the cumulative distribution function values
  ! "cdf", and return in val. this version works on 2d blocks of data.
  subroutine sample_from_pdf_masked_block(this, nz, ng, fsd, cdf, val, mask)
    
    class(pdf_sampler_type), intent(in)  :: this

    ! number of samples
    integer,    intent(in) :: nz, ng

    ! fractional standard deviation (0 to 4) and cumulative
    ! distribution function (0 to 1)
    real(jprb), intent(in)  :: fsd(nz), cdf(ng, nz)

    ! sample from distribution
    real(jprb), intent(out) :: val(:,:)

    ! mask
    logical,    intent(in), optional :: mask(nz)

    ! loop index
    integer    :: jz, jg

    ! index to look-up table
    integer    :: ifsd, icdf

    ! weights in bilinear interpolation
    real(jprb) :: wfsd, wcdf

    do jz = 1,nz

      if (mask(jz)) then
        
        do jg = 1,ng
          if (cdf(jg, jz) > 0.0_jprb) then
            ! bilinear interpolation with bounds
            wcdf = cdf(jg,jz) * (this%ncdf-1) + 1.0_jprb
            icdf = max(1, min(int(wcdf), this%ncdf-1))
            wcdf = max(0.0_jprb, min(wcdf - icdf, 1.0_jprb))
          
            wfsd = (fsd(jz)-this%fsd1) * this%inv_fsd_interval + 1.0_jprb
            ifsd = max(1, min(int(wfsd), this%nfsd-1))
            wfsd = max(0.0_jprb, min(wfsd - ifsd, 1.0_jprb))
            
            val(jg,jz)=(1.0_jprb-wcdf)*(1.0_jprb-wfsd) * this%val(icdf  ,ifsd)   &
                 &    +(1.0_jprb-wcdf)*          wfsd  * this%val(icdf  ,ifsd+1) &
                 &    +          wcdf *(1.0_jprb-wfsd) * this%val(icdf+1,ifsd)   &
                 &    +          wcdf *          wfsd  * this%val(icdf+1,ifsd+1)
          else
            val(jg,jz) = 0.0_jprb
          end if
        end do

      end if

    end do

  end subroutine sample_from_pdf_masked_block

end module radiation_pdf_sampler
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

