! Utility module: Print namelist to file, annotating all changed values.
!
! Every ICON run generates annotated lists of namelist parameters
! during the setup. These lists are written to text files
! "nml.atmo.log", "nml.cpl.log", "nml.ocean.log" and have the
! following form:
!
!     NAMELIST IO_NML
!         OUT_EXPNAME          'case4                                   [...]' (truncated)
!                              >> DEFAULT: 'IIIEEEETTTT                             [...]' (truncated)
!         OUT_FILETYPE         2
!         LKEEP_IN_SYNC        F
!         DT_DATA              43200.000000000000
!                              >> DEFAULT: 21600.000000000000
!         DT_DIAG              1728000.0000000000
!                              >> DEFAULT: 86400.000000000000
!
! and so on.
! The "DEFAULT" annotation denotes all those parameters that have
! been modified by the user; in this case the default value of the
! namelist parameter is stated together with the modified value.  All
! other namelist parameters are listed only with their default value.
!
! Note that it is not necessary to update this module whenever a namelist is
! extended: Lists are generated automatically.
!
! The implementation is more or less a "poor man's approach". This
! module writes and re-reads the namelists in text form, since
! Fortran is missing the language features for looping over namelist
! elements.
!
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

MODULE mo_nml_annotate

  USE mo_io_units,    ONLY: find_next_free_unit, filename_max
  USE mo_exception,   ONLY: finish
  USE mo_util_file,   ONLY: util_tmpnam, util_unlink
  USE mo_util_nml,    ONLY: util_annotate_nml
  
  IMPLICIT NONE
  PRIVATE
  

  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_nml_annotate')

  INTEGER :: tmpnml = -1   !< file handle for temporary text file with defaults and settings
  
  PUBLIC :: temp_defaults
  PUBLIC :: temp_settings
  PUBLIC :: log_nml_settings 
  
CONTAINS
  
  !> Opens a new temporary text file.
  !
  FUNCTION temp_defaults()
    INTEGER :: temp_defaults
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::temp_defaults'
    LOGICAL                     :: lopen
    CHARACTER(len=filename_max) :: filename
    INTEGER                     :: flen

    IF (tmpnml == -1) THEN
      lopen = .FALSE.
    ELSE
      INQUIRE(tmpnml, OPENED=lopen)
    END IF
    IF (.NOT. lopen) THEN
      flen   = util_tmpnam(filename, filename_max)
      tmpnml = find_next_free_unit(10,100)
      IF (tmpnml < 0) THEN
        CALL finish(routine, "Failed call to find_next_free_unit.")
      END IF
      OPEN(UNIT=tmpnml, FILE=filename(1:flen), &
        &  ACCESS='sequential',   action='write', delim='apostrophe')
    END IF
    temp_defaults = tmpnml
  END FUNCTION temp_defaults
  
  
  !> Opens a new temporary text file.
  !
  FUNCTION temp_settings()
    INTEGER :: temp_settings
    ! Local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::temp_settings'
    LOGICAL :: lopen

    INQUIRE(tmpnml, OPENED=lopen)
    IF (.NOT. lopen) CALL finish(routine, "Internal error!")
    temp_settings = tmpnml
  END FUNCTION temp_settings
  
  
  !> Read defaults and settings into string buffer, compare them and do the
  !  print-out.
  !
  SUBROUTINE log_nml_settings(dst_filename)
    CHARACTER(LEN=*), INTENT(IN) ::dst_filename !< log file name (optional)
    ! Local variables
    CHARACTER(len=filename_max) :: tmp_filename
    INTEGER                     :: iret
    LOGICAL                     :: lopen

    INQUIRE(tmpnml, NAME=tmp_filename)
    INQUIRE(tmpnml, OPENED=lopen)
    IF (lopen) THEN
      CLOSE (tmpnml)
      tmpnml = -1
      ! call a C routine which reads and parses the temporary
      ! namelists with a finite state machine:
      iret = util_annotate_nml(TRIM(tmp_filename), TRIM(dst_filename))
      iret = util_unlink(TRIM(tmp_filename))
    END IF
  END SUBROUTINE log_nml_settings
  
END MODULE mo_nml_annotate
