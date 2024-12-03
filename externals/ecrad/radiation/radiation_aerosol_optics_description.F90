! # 1 "radiation/radiation_aerosol_optics_description.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_aerosol_optics_description.f90"
! radiation_aerosol_optics_description.f90 - type to store aerosol optics metadata
!
! (c) copyright 2022- ecmwf.
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
! # 17 "radiation/radiation_aerosol_optics_description.f90" 2

module radiation_aerosol_optics_description

  use parkind1,      only : jprb

  implicit none
  public

  !---------------------------------------------------------------------
  ! this type holds the metadata from an aerosol optical property
  ! file, enabling the user to request the index to the aerosol type
  ! with particular properties.  note that string information is held
  ! in the form of 2d arrays of single characters, so comparison to
  ! character strings requires the to_string helper function at the
  ! end of this file.
  type aerosol_optics_description_type

    ! two-character code describing the aerosol family, dimensioned
    ! (2,naer), e.g.
    !   ss: sea salt
    !   om: organic matter
    !   su: sulfate
    !   ob: secondary organic biogenic
    !   oa: secondary organic anthropogenic
    !   am: fine-mode ammonium sulfate
    !   ni: nitrate
    !   dd: desert dust
    !   bc: black carbon
    character(len=1), allocatable :: code_phobic(:,:)
    character(len=1), allocatable :: code_philic(:,:)

    ! size bin, typically 1-2 or 1-3 in order fine to coarse, or zero
    ! if no division by size is used, dimensioned (naer)
    integer, allocatable :: bin_phobic(:)
    integer, allocatable :: bin_philic(:)

    ! character string characterizing the optical model, e.g. opac,
    ! gacp, glomap, dubovik2002 etc.
    character(len=1), allocatable :: optical_model_phobic(:,:)
    character(len=1), allocatable :: optical_model_philic(:,:)

    ! the user can call preferred_optical_model to specify that a
    ! certain optical model for a certain aerosol family is to be
    ! preferred when get_index is called
    logical, allocatable :: is_preferred_phobic(:)
    logical, allocatable :: is_preferred_philic(:)

    ! verbosity level
    integer :: iverbose
    
  contains
    procedure :: read
    procedure :: preferred_optical_model
    procedure :: get_index

  end type aerosol_optics_description_type

contains

  !---------------------------------------------------------------------
  ! read optical property file file_name into an
  ! aerosol_optics_description_type object
  subroutine read(this, file_name, iverbose)

    use ecradhook,              only : lhook, dr_hook, jphook



    use easy_netcdf,          only : netcdf_file


    class(aerosol_optics_description_type), intent(inout) :: this
    character(len=*), intent(in)              :: file_name
    integer, intent(in), optional             :: iverbose
    
    ! the netcdf file containing the aerosol optics data
    type(netcdf_file)  :: file

    real(jphook)       :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics_description:load',0,hook_handle)

    ! open the aerosol scattering file and configure the way it is
    ! read
    call file%open(trim(file_name), iverbose=iverbose)

    ! read metadata variables
    call file%get('code_hydrophilic', this%code_philic)
    call file%get('code_hydrophobic', this%code_phobic)
    call file%get('bin_hydrophilic',  this%bin_philic)
    call file%get('bin_hydrophobic',  this%bin_phobic)
    call file%get('optical_model_hydrophilic', this%optical_model_philic)
    call file%get('optical_model_hydrophobic', this%optical_model_phobic)

    ! allocate logical arrays of the appropriate size and set to false
    allocate(this%is_preferred_philic(size(this%bin_philic)))
    allocate(this%is_preferred_phobic(size(this%bin_phobic)))
    this%is_preferred_philic = .false.
    this%is_preferred_phobic = .false.

    call file%close()

    if (present(iverbose)) then
      this%iverbose = iverbose
    else
      this%iverbose = 3
    end if
    
    if (lhook) call dr_hook('radiation_aerosol_optics_description:load',1,hook_handle)

  end subroutine read

  !---------------------------------------------------------------------
  ! specify the preferred optical model for a particular aerosol
  ! family, e.g. "call
  ! aer_desc%preferred_optical_model('dd','woodward2001')" would mean
  ! that subsequent calls to get_index in which the optical model is
  ! not specified would return the woodward model rather than the
  ! first matching model in the file.  the check is only done on the
  ! first len(optical_model_str) characters, so "woodward" and
  ! "woodward2001" would both match the woodward2001 model.
  subroutine preferred_optical_model(this, code_str, optical_model_str)

    use ecradhook,              only : lhook, dr_hook, jphook
    use radiation_io,         only : nulout, nulerr, radiation_abort
    
    class(aerosol_optics_description_type), intent(inout) :: this
    character(len=2), intent(in) :: code_str
    character(len=*), intent(in) :: optical_model_str

    ! aerosol loop counter
    integer :: ja

    logical :: is_found, is_philic, is_phobic

    real(jphook)         :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics_description:preferred_optical_model',0,hook_handle)

    ! check for empty string
    if (optical_model_str == ' ') then
      if (lhook) call dr_hook('radiation_aerosol_optics_description:preferred_optical_model',1,hook_handle)
      return
    end if

    is_found  = .false.
    is_philic = .false.
    is_phobic = .false.
    
    ! loop over hydrophilic types
    do ja = 1,size(this%bin_philic)
      ! check if we have a match
      if (to_string(this%code_philic(:,ja)) == code_str &
           &  .and. trim(to_string(this%optical_model_philic(:,ja))) &
           &          == optical_model_str) then
        this%is_preferred_philic(ja) = .true.
        is_found  = .true.
        is_philic = .true.
      end if
    end do
    ! repeat for the hydrophobic types
    do ja = 1,size(this%bin_phobic)
      if (to_string(this%code_phobic(:,ja)) == code_str &
           &  .and. trim(to_string(this%optical_model_phobic(:,ja))) &
           &          == optical_model_str) then
        this%is_preferred_phobic(ja) = .true.
        is_found  = .true.
        is_phobic = .true.
      end if
    end do

    if (.not. is_found) then
      write(nulerr,'(a,a2,a,a,a)') '*** error: preferred "', code_str ,'" aerosol optical model "', &
           &  trim(optical_model_str), '" not found in file'
      call radiation_abort()
    else if (this%iverbose > 2) then
      write(nulout,'(a,a2,a,a,a)',advance='no') 'preferred "', code_str, '" aerosol optical model set to "', &
           &  trim(optical_model_str), '" ('
      if (is_phobic) then
        write(nulout,'(a)',advance='no') ' hydrophobic'
      end if
      if (is_philic) then
        write(nulout,'(a)',advance='no') ' hydrophilic'
      end if
      write(nulout,'(a)') ' )'
    end if
    
    if (lhook) call dr_hook('radiation_aerosol_optics_description:preferred_optical_model',1,hook_handle)

  end subroutine preferred_optical_model

  
  !---------------------------------------------------------------------
  ! return the index to the aerosol optical properties corresponding
  ! to the aerosol family in code_str (e.g. ss, dd etc), whether or
  ! not the requested aerosol is hydrophilic in the logical
  ! lhydrophilic, and optionally the size bin ibin and optical model
  ! in optical_model_str. the return value may be used to populate the
  ! radiation_config%i_aerosol_map vector, where a positive number is
  ! a hydrophobic index, a negative number is a hydrophilic index and
  ! zero indicates that the aerosol type was not found in the file.
  ! this is a valid entry in i_aerosol_map meaning the aerosol is
  ! ignored, but the calling routine to get_index might wish to throw
  ! a warning or error. this routine works by assigning a score based
  ! on the closeness of the match.
  function get_index(this, code_str, lhydrophilic, ibin, optical_model_str)
    
    use ecradhook,              only : lhook, dr_hook, jphook
    use radiation_io,         only : nulout

    class(aerosol_optics_description_type), intent(in) :: this
    character(len=2), intent(in) :: code_str
    logical, intent(in) :: lhydrophilic
    integer, intent(in), optional :: ibin
    character(len=*), intent(in), optional :: optical_model_str

    ! score of the currently selected aerosol index, and the score of
    ! the current one under consideration
    integer :: score, current_score

    ! loop index for aerosol type
    integer :: ja

    ! return value
    integer :: get_index

    ! issue a warning if there is more than one equal match
    logical :: is_ambiguous

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_aerosol_optics_description:get_index',0,hook_handle)

    ! initial values
    get_index = 0
    score = 0
    is_ambiguous = .false.

    if (lhydrophilic) then
      ! loop over hydrophilic aerosol types
      do ja = 1,size(this%bin_philic)
        current_score = 0
        if (to_string(this%code_philic(:,ja)) == code_str) then
          ! aerosol code matches
          if (present(ibin) .and. this%bin_philic(ja) > 0) then
            if (ibin > 0) then
              if (ibin == this%bin_philic(ja)) then
                ! requested bin number matches
                current_score = 4
              else
                ! requested bin number does not match
                current_score = -1
              end if
            else
              ! bin number is zero: no request
              current_score = 2
            end if
          else
            ! no bin number present
            current_score = 2
          end if
          if (present(optical_model_str)) then
            if (trim(to_string(this%optical_model_philic(:,ja))) &
                 &  == optical_model_str) then
              ! requested optical model matches
              if (current_score >= 0) then
                current_score = current_score + 4
              end if
            else
              ! requested optical model does not match
              current_score = -1
            end if
          else if (current_score >= 0) then
            ! no requested optical model
            current_score = current_score + 2
          end if
          if (current_score > 0 .and. this%is_preferred_philic(ja)) then
            current_score = current_score + 1
          end if
          if (current_score > score) then
            ! better score than any existing aerosol type
            get_index = -ja
            score = current_score
            is_ambiguous = .false.
          else if (current_score > 0 .and. current_score == score) then
            is_ambiguous = .true.
          end if
        end if
      end do
    else
      ! loop over hydrophobic aerosol types
      do ja = 1,size(this%bin_phobic)
        current_score = 0
        if (to_string(this%code_phobic(:,ja)) == code_str) then
          ! aerosol code matches
          if (present(ibin) .and. this%bin_phobic(ja) > 0) then
            if (ibin > 0) then
              if (ibin == this%bin_phobic(ja)) then
                ! requested bin number matches
                current_score = 4
              else
                ! requested bin number does not match
                current_score = -1
              end if
            else
              ! bin number is zero: no request
              current_score = 2
            end if
          else
            ! no bin number requested or present
            current_score = 2
          end if
          if (present(optical_model_str)) then
            if (trim(to_string(this%optical_model_phobic(:,ja))) &
                 &  == optical_model_str) then
              ! requested optical model matches
              if (current_score >= 0) then
                current_score = current_score + 4
              end if
            else
              ! requested optical model does not match
              current_score = -1
            end if
          else if (current_score >= 0) then
            ! no requested optical model
            current_score = current_score + 2
          end if
          if (current_score > 0 .and. this%is_preferred_phobic(ja)) then
            current_score = current_score + 1
          end if
          if (current_score > score) then
            ! better score than any existing aerosol type
            get_index = ja
            score = current_score
            is_ambiguous = .false.
          else if (current_score > 0 .and. current_score == score) then
            is_ambiguous = .true.
          end if          
        end if
      end do
    end if

    if (is_ambiguous) then
      write(nulout,'(a,a2,a,l1,a)') 'warning: radiation_aerosol_optics_description:get_index("', &
           &  code_str, '",', lhydrophilic, &
           &  ',...) does not unambiguously identify an aerosol optical property index'
    end if
    
    if (lhook) call dr_hook('radiation_aerosol_optics_description:get_index',1,hook_handle)

  end function get_index

  !---------------------------------------------------------------------
  ! utility function to convert an array of single characters to a
  ! character string (yes fortran's string handling is a bit
  ! rubbish). we set null characters (ascii code 0) returned from the
  ! netcdf library to spaces, so that trim can remove them.
  pure function to_string(arr) result(str)
    character, intent(in)  :: arr(:)
    character(len=size(arr)) :: str
    integer :: jc
    do jc = 1,size(arr)
      if (ichar(arr(jc)) == 0) then
        ! replace null character with a space
        str(jc:jc) = ' '
      else
        str(jc:jc) = arr(jc)
      end if
    end do
  end function to_string

end module radiation_aerosol_optics_description
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

