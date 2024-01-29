!
! mo_art_create_filenames
! This module provides an automatic filename creation for:
! - Input files for initialization of aerosols and gases
! - External parameter files
! - Emission Files (time-dependent)
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE mo_art_create_filenames
! ICON
  USE mo_exception,                ONLY: finish
  USE mo_art_config,               ONLY: art_config
! ART
  USE mo_art_io_constants,         ONLY: IINIT_AERO, IINIT_CHEM,                &
                                     &   IEXT_POCON, IEXT_POVAR, IEXT_POAMB,    &
                                     &   art_get_abbr_string, IRES_LEN
  USE mo_art_data,                 ONLY: p_art_data
  USE mo_art_atmo_data,            ONLY: t_art_atmo

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_create_filenames'

  PUBLIC :: art_create_filenames, art_check_filename
  PUBLIC :: art_get_res_string
  PUBLIC :: art_set_io_suffix

  CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_create_filenames(jg,cart_input_folder,type_of_data, &
  &                             cinput_filestring,ltime_dep,datetime)
!<
! SUBROUTINE art_create_filenames
! This subroutine automatically creates the file names according to the
! ICON-ART name convention: ART_<XXX>_iconR<n>B<kk>-grid-<yyyy-mm-dd-hh>_<grid_string>.nc
! where <XXX> stands for abbreviation of type, <n> and <kk> indicate ICON resolution,
! <yyyy-mm-dd-hh> is the date (if available) and <grid_string> is the number of grid used
! or a user specified string via the namelist parameter art_nml:cart_io_suffix(p_patch%id)
! 
! For datasets without time dependence file names are created according to
! ART_<XXX>_iconR<n>B<kk>-grid_<grid_string>.nc
!
! Part of Module: mo_art_create_filenames
! Author: Daniel Rieger, KIT
! Initial Release: 2015-02-10
! Modifications:
! 2015-02-26: Michael Weimer, KIT
! - Added ANT,BIO,BBE and change name creation with grid_filename
! 2016-04-11: Michael Weimer, KIT
! - included number_of_grid_used and adaption to new namelist parameter, create filenames 
!   now with nroot and p_patch%level
!>
  INTEGER, INTENT(in)                 :: &
    &  jg                                  !< patch on which computation is performed
  CHARACTER(LEN=*), INTENT(in)        :: &
    &  cart_input_folder                   !< Input folder of ART initialization and external files
  INTEGER, INTENT(in)                 :: &
    &  type_of_data                        !< Type of dataset (see mo_art_io_constants.f90)
  CHARACTER(LEN=*), INTENT(out)       :: &
    &  cinput_filestring                   !< Path and Filename of dataset that will be read
  LOGICAL, INTENT(in)                 :: &
    &  ltime_dep                           !< Time-dependent data? -> Timestring has to be attached
                                           !    to filename
  CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
    &  datetime                            !< Datetime string for time-dependent data
! Local variables
  CHARACTER(LEN=3)    :: &
    &  typestring          !< Label of dataset 
  CHARACTER(LEN=120)  :: &
    &  datasetname         !< Name of dataset (without path)
  CHARACTER(LEN=20)   :: &
    &  grid_part           !< p_patch%grid_filename without suffix if it exists

  typestring = art_get_abbr_string(type_of_data)
  
  grid_part =  'icon'//TRIM(art_get_res_string(jg))//'-grid'

  IF (ltime_dep) THEN
    IF (.NOT. PRESENT(datetime)) CALL finish('mo_art_create_filenames:art_create_filenames', &
      &          'Passed time-dependent data without datetime')
    datasetname = 'ART_'//TRIM(typestring)//'_'//TRIM(grid_part)//'-'//TRIM(datetime)
  ELSE
    datasetname = 'ART_'//TRIM(typestring)//'_'//TRIM(grid_part)
  ENDIF

  datasetname = TRIM(datasetname)//'_'//TRIM(art_config(jg)%cart_io_suffix)

  ! In case of initial data, always return the grib filename. The read routine will try to open
  ! the grib file. If no grib file is available it will try to open the netcdf file.
  ! Non-initialdata filenames are always returned as netcdf filenames
  IF (ANY(type_of_data == (/IINIT_AERO, IINIT_CHEM,                 &
    &                       IEXT_POCON, IEXT_POVAR, IEXT_POAMB/) )) THEN
    cinput_filestring = TRIM(cart_input_folder)//'/'//TRIM(datasetname)
  ELSE    
    cinput_filestring = TRIM(cart_input_folder)//'/'//TRIM(datasetname)//'.nc'
  END IF

END SUBROUTINE art_create_filenames
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_check_filename(filename, lrequired, lexist)
!<
! SUBROUTINE art_check_filename
! This subroutine checks whether a file really exists
! Part of Module: mo_art_create_filenames
! Author: Michael Weimer, KIT, Daniel Rieger, KIT
! Initial Release: 2016-11-16
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  CHARACTER(LEN=*),INTENT(IN)    :: &
    &  filename                       !< Path and name of file to check
  LOGICAL, INTENT(in), OPTIONAL  :: &
    &  lrequired                      !< Is this file required?
  LOGICAL, INTENT(out), OPTIONAL :: &
    &  lexist                         !< Does this filename exist (return value)?
!Local variables
  LOGICAL                        :: &
    &  lrequired_d,                 & !< Default value for lrequired
    &  lexist_loc                     !< Does this filename exist (local)?

  lrequired_d = .TRUE.

  IF(PRESENT(lrequired)) lrequired_d = lrequired

  INQUIRE(FILE = TRIM(filename), EXIST = lexist_loc)
  IF (lrequired_d) THEN ! stop only if file is required
    IF (.NOT. lexist_loc) THEN
      CALL finish(TRIM(routine)//':art_check_filename',  &
             &    TRIM(filename)//' does not exist.')
    ENDIF
  ENDIF

  IF(PRESENT(lexist)) lexist = lexist_loc

END SUBROUTINE art_check_filename
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_set_io_suffix(cart_io_suffix, number_of_grid_used)
!<
! SUBROUTINE art_set_io_suffix
! This routine sets the string of the io suffix
!
! Part of Module: mo_art_create_filenames
! Author: Michael Weimer, KIT
! Initial Release: 2018-06-21
! Modifications:
!>
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(inout) :: &
    &  cart_io_suffix
  INTEGER, INTENT(in) :: &
    &  number_of_grid_used

  IF (LEN_TRIM(cart_io_suffix) > 0) THEN
    IF (TRIM(cart_io_suffix) == 'grid-number') THEN
      WRITE(cart_io_suffix,'(I4.4)')  number_of_grid_used
    END IF
  ELSE
    cart_io_suffix = ''
  END IF

END SUBROUTINE art_set_io_suffix
!!
!!-------------------------------------------------------------------------
!!
CHARACTER(LEN=IRES_LEN) FUNCTION art_get_res_string(jg)
!<
! FUNCTION art_get_res_string
! This function creates a string with the current ICON resolution (e.g. "R2B06")
!
! Part of Module: mo_art_create_filenames
! Author: Michael Weimer, KIT
! Initial Release: 2016-10-27
! Modifications:
!>
  IMPLICIT NONE
  INTEGER, INTENT(in)  :: jg
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  art_atmo => p_art_data(jg)%atmo

  IF ( art_atmo%nroot < 10 ) THEN
    WRITE(art_get_res_string, '(A1,I1.1,A1,I2.2)') 'R',art_atmo%nroot,'B',art_atmo%nbisect
  ELSE
    WRITE(art_get_res_string, '(A1,I2.2,A1,I2.2)') 'R',art_atmo%nroot,'B',art_atmo%nbisect
  END IF
END FUNCTION art_get_res_string
END MODULE mo_art_create_filenames
