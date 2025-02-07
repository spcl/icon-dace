! mo_namelist.f90 - Routines for handling of namelists
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!_________________________________________

MODULE mo_namelist
  ! position_nml - position namelist file for reading
  !
  ! Author:
  !
  ! L. Kornblueh, MPI, March 2001, original source
  ! L. Kornblueh, MPI, January 2003, added error message for open
  ! L. Kornblueh, MPI, April 2005, added close_nml
  ! L. Kornblueh, MPI, February 2011, make it more usable, when 
  !               having multiple files open
  !
  !
  !  This routine originates (year 2014) from MPI-ESM, the Earth System Model of the 
  !  Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
  !  Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
  !  version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
  !  doi: 10.1029/2018MS001400.

  USE mo_util_string,   ONLY: tolower
  USE mo_filename,      ONLY: find_next_free_unit  
  USE mo_exception,     ONLY: message, finish, message_text

  IMPLICIT NONE

  PRIVATE

  ! accessible functions

  PUBLIC :: position_nml
  PUBLIC :: open_nml
  PUBLIC :: close_nml

  ! return values from position_nml

  PUBLIC :: positioned, missing, closed, length_error, read_error

  TYPE t_namelist
    CHARACTER(len=64) :: filename
    INTEGER :: fileunit
  END TYPE t_namelist

  TYPE(t_namelist), ALLOCATABLE :: namelist_file_list(:) 
  
  INTEGER, PARAMETER :: chunk = 32
  INTEGER, SAVE :: top = 0

  ! return values of function 'position_nml':

  INTEGER, PARAMETER :: positioned   =  0 ! file pointer set to namelist group
  INTEGER, PARAMETER :: missing      =  1 ! namelist group is missing
  INTEGER, PARAMETER :: closed       =  2 
  INTEGER, PARAMETER :: length_error =  3 
  INTEGER, PARAMETER :: read_error   =  4

CONTAINS

  !================================================================================================
  FUNCTION open_nml (file) RESULT(unit)
    INTEGER :: unit                         ! file handler
    CHARACTER(len=*), INTENT(in)  :: file   ! namelist filename

    TYPE(t_namelist), ALLOCATABLE :: tmp_list(:) 
    INTEGER :: ihandler, istat, iunit

    unit = -1

    ! Create or extend table of namelist files

    IF (.NOT. ALLOCATED(namelist_file_list)) THEN
      ALLOCATE(namelist_file_list(chunk))
    ELSE
      IF (MOD(top,chunk) == 0) THEN
        CALL move_alloc(namelist_file_list, tmp_list)
        ALLOCATE(namelist_file_list(top+chunk))
        namelist_file_list(:top) = tmp_list(:)
        DEALLOCATE(tmp_list)
      ENDIF
    ENDIF

    ! Look up file name in table, re-open if close, otherwise go to begin.
    ! Return from procedure if an entry was found.

    IF (top > 0) THEN
      DO ihandler = 1, top
        IF (TRIM(namelist_file_list(ihandler)%filename) == file) THEN
          IF (namelist_file_list(ihandler)%fileunit == -1) THEN
            iunit = find_next_free_unit(100,200)
            OPEN (iunit, file=file, iostat=istat, status='old', action='read', delim='apostrophe')
            namelist_file_list(ihandler)%fileunit = iunit
            CALL message('', 'Namelist file '//TRIM(file)//' reopened')
          ELSE
            CALL message('', 'Namelist file '//TRIM(file)//' already opened, rewind file to beginning')
            REWIND (namelist_file_list(ihandler)%fileunit)
          ENDIF
          unit = ihandler
          RETURN
        ENDIF  
      ENDDO
    ENDIF

    ! If no entry was found, open file and create table entry.

    iunit = find_next_free_unit(100,200)
    OPEN (iunit, file=file, iostat=istat, status='old', action='read', delim='apostrophe')

    IF (istat == 0) THEN
      top = top+1
      namelist_file_list(top)%filename = TRIM(file)
      namelist_file_list(top)%fileunit = iunit      
      unit = top
      CALL message('', 'Opened namelist file '//TRIM(file))
    ELSE
      CALL finish('open_nml','Could not open namelist file '//TRIM(file))
    ENDIF

  END FUNCTION open_nml
  !================================================================================================
  FUNCTION position_nml (name, unit, rewind, status) RESULT (iunit)
    
    INTEGER :: iunit ! Fortran unit to read namelist from

    ! position namelist file for reading on first occurrence
    ! namelist /name/ (case independent). 

    CHARACTER(len=*), INTENT(in)            :: name   ! namelist group name
    INTEGER,          INTENT(in)            :: unit   ! file handler
    LOGICAL,          INTENT(in)  ,OPTIONAL :: rewind ! default: true
    INTEGER,          INTENT(out) ,OPTIONAL :: status ! error return value

    CHARACTER(len=256) :: yline     ! line read
    CHARACTER(len=256) :: test      ! potentially uppercase namelist group name
    INTEGER            :: stat      ! local copy of status variable
    INTEGER            :: ios       ! status variable from read operation
    LOGICAL            :: lrew      ! local copy of rewind flag
    CHARACTER(len=64)  :: yfilename ! local copy of namelist filename 
    INTEGER            :: len_name  ! length of requested namelist group name
    CHARACTER          :: ytest     ! character to test for delimiter
    CHARACTER(len=12)  :: action    ! current action 
    INTEGER            :: ind       ! index from index routine
    INTEGER            :: indc      ! index of comment character (!)
    CHARACTER(len=256) :: iomsg

    IF (namelist_file_list(unit)%fileunit == -1) THEN
      CALL finish('position_nml','Namelist file '//TRIM(namelist_file_list(unit)%filename)//' already closed')
    ELSEIF (unit > top) THEN
      WRITE (message_text,'(a,i0,a)') 'Namelist file for handler ', unit, ' has not been opened'
      CALL finish('position_nml',message_text)
    ELSE
      iunit     = namelist_file_list(unit)%fileunit
      yfilename = TRIM(namelist_file_list(unit)%filename)
    ENDIF

    lrew  = .TRUE.   
    IF (PRESENT(rewind)) lrew  = rewind
    stat  =  missing
    action  = 'missing'

    len_name = LEN_TRIM (name)

    IF (len_name > LEN(test)) THEN
      stat =  length_error
      action = 'length_error'
    ENDIF

    test = tolower (name)
    
    ! Reposition file at beginning:

    IF (lrew) THEN
      REWIND (iunit)
    ENDIF

    ! Search start of namelist
    
    DO
      IF (stat /= MISSING) EXIT

      yline = ''
    
      READ (iunit,'(a)',IOSTAT=ios,IOMSG=iomsg) yline

      IF (ios < 0) THEN
        EXIT  ! missing
      ELSE IF (ios > 0) THEN
        stat =  read_error
        action = 'read_error'
        EXIT
      END IF

      yline = tolower (yline)
    
      ind = INDEX(yline,'&'//TRIM(test))
    
      IF (ind == 0) CYCLE
      
      indc = INDEX(yline,'!')

      IF (indc > 0 .AND. indc < ind) CYCLE

      ! test for delimiter
    
      ytest = yline(ind+len_name+1:ind+len_name+1)

      IF ( (LGE(ytest,'0') .AND. LLE(ytest,'9')) .OR. &
           (LGE(ytest,'a') .AND. LLE(ytest,'z')) .OR. &
            ytest == '_'                         .OR. &
           (LGE(ytest,'A') .AND. LLE(ytest,'Z'))) THEN
        CYCLE
      ELSE 
        stat = positioned
        BACKSPACE (iunit)
        EXIT
      END IF
    ENDDO
    
    IF (PRESENT (status)) status = stat
      SELECT CASE (stat)
      CASE (positioned)
        RETURN
      CASE (missing)
        IF (PRESENT (status)) RETURN
      END SELECT

     CALL message ('position_nml','File '//TRIM(yfilename)// 'namelist /'//TRIM(test)//'/ '//TRIM(action))

   END FUNCTION position_nml
  !================================================================================================
  SUBROUTINE close_nml (unit)

    INTEGER, INTENT(in) :: unit ! file handler
    INTEGER :: istat

    IF (unit == -1) THEN
      CALL finish('close_nml','Namelist file not opened')
    ELSE
      IF (namelist_file_list(unit)%fileunit == -1) THEN
        CALL finish ('close_nml','Namelist file '//TRIM(namelist_file_list(unit)%filename)//' already closed')
      ENDIF
    ENDIF

    CLOSE (namelist_file_list(unit)%fileunit, IOSTAT=istat)

    IF (istat == 0) THEN
      namelist_file_list(unit)%fileunit = -1
    ELSE
      CALL finish('close_nml','Could not close namelist file '//TRIM(namelist_file_list(unit)%filename))
    ENDIF

  END SUBROUTINE close_nml
  !================================================================================================
END MODULE mo_namelist
