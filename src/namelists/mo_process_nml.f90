! Process Fortran namelists containing model configuration parameters from an external file.
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

MODULE mo_process_nml

  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_control,     ONLY: use_restart_namelists
  USE mo_restart_nml_and_att,ONLY: open_tmpfile, store_and_close_namelist, &
                                 & open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,       ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: process_nml

CONTAINS
  !>
  !!
  SUBROUTINE process_nml( filename, nml_name, nml_read, nml_write )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: nml_name

    INTERFACE
       !
       SUBROUTINE nml_read (funit)
         INTEGER, INTENT(in) :: funit
       END SUBROUTINE nml_read
       !
       SUBROUTINE nml_write(funit)
         INTEGER, INTENT(in) :: funit
       END SUBROUTINE nml_write
       !
    END INTERFACE
    
    INTEGER :: istat
    INTEGER :: funit
    INTEGER :: iunit

    !------------------------------------------------------------------
    ! 1. If this is a resumed integration, overwrite the initialization
    !    values with values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
       funit = open_and_restore_namelist(nml_name)
       !
       CALL nml_read(funit)       ! read old settings
       !
       CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------
    ! 2. Read user's (new) specifications (done by all MPI processes)
    !------------------------------------------------------------------
    CALL open_nml(filename)
    CALL position_nml (nml_name, STATUS=istat)
    IF (my_process_is_stdio()) THEN
       iunit = temp_defaults()
       !
       CALL nml_write(iunit)      ! write defaults to temporary text file
       !
    END IF
    SELECT CASE (istat)
    CASE (positioned)
       !
       CALL nml_read (nnml)       ! read new settings
       !
       IF (my_process_is_stdio()) THEN
          iunit = temp_settings()
          !
          CALL nml_write(iunit)   ! write settings to temporary text file
          !
       END IF
    END SELECT
    CALL close_nml

    !------------------------------------------------------------------
    ! 3. Store the namelist for restart
    !------------------------------------------------------------------
    IF (my_process_is_stdio()) THEN
       funit = open_tmpfile()
       !
       CALL nml_write(funit)
       !
       CALL store_and_close_namelist(funit, nml_name)
    END IF

    !------------------------------------------------------------------
    ! 4. Write the namelist to an ASCII file
    !------------------------------------------------------------------
    IF ( my_process_is_stdio() ) THEN
       !
       CALL nml_write(nnml_output)
       !
    END IF

  END SUBROUTINE process_nml

END MODULE mo_process_nml
