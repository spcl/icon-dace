! mo_filename.f90- Quantities for file names etc.
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!_________________________________________

MODULE mo_filename

  ! -----------------------------------------------------------------
  !
  ! module *mo_filename* - quantities needed for file names etc.
  !
  ! -----------------------------------------------------------------
  !
  ! L. Kornblueh, MPI, April 1998, added NWP forecast mode
  ! I. Kirchner,  MPI, April 2001, revision
  ! A. Rhodin,    MPI, June  2001, adapted to stream interface
  !
  !  This routine originates (year 2014) from MPI-ESM, the Earth System Model of the 
  !  Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
  !  Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
  !  version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
  !  doi: 10.1029/2018MS001400.

  USE mo_time_control,  ONLY: next_date, start_date, str_date,     &
                              timelabel_type, get_date_components, &
                              get_forecast_hours,                  &
                              FLT_YMD, FLT_NWP_YMD
  USE mo_control,       ONLY: lnwp       ! prediction mode flag
  USE mo_exception,     ONLY: message_text, message, finish

  IMPLICIT NONE

  PRIVATE                           ! used by:
  PUBLIC :: yomdn                   ! mo_io.f90 labrun.f90
  PUBLIC :: find_next_free_unit     ! mo_nudging_init.f90 mo_nudging_io.f90
  PUBLIC :: out_expname             ! mo_nudging_init.f90 inictl.f90
  PUBLIC :: str_filter              ! mo_nudging_init.f90 mo_nudging_io.f90 ..
  PUBLIC :: out_datapath            ! inictl.f90
  PUBLIC :: rerun_filetype          ! filetype for rerun (NETCDF)
  PUBLIC :: out_filetype            ! filetype for echam output(GRIB or NETCDF)
  PUBLIC :: trac_filetype           ! filetype for tracers     (GRIB or NETCDF)
  PUBLIC :: out_ztype               ! compress type for echam GRIB output(NONE or SZIP)
                                    ! compress type for echam netCDF4 output(NONE or ZIP)
  PUBLIC :: compose_filenames       ! stepon.f90, mo_output.f90
  PUBLIC :: yexp, ypath, syr, smo, &! subjob.f90
            sdy, shr, smn, sse      !
  PUBLIC :: standard_output_file    ! mo_output.f90
  PUBLIC :: standard_rerun_file     !
  PUBLIC :: GRIB, NETCDF, NETCDF2, NETCDF4  ! Valid values for file types
  PUBLIC :: NONE, SZIP, ZIP         ! Valid values for compress types
  PUBLIC :: name_limit, path_limit

  INTEGER ,PARAMETER :: GRIB      = 1       ! Valid values for 'filetype'
  INTEGER ,PARAMETER :: NETCDF    = 2       !
  INTEGER ,PARAMETER :: NETCDF2   = 4       !
  INTEGER ,PARAMETER :: NETCDF4   = 6       !

  INTEGER ,PARAMETER :: NONE      = 0       ! Valid values for 'compresstype'
  INTEGER ,PARAMETER :: SZIP      = 1       ! only available for GRIB output
  INTEGER ,PARAMETER :: ZIP       = 2       ! only available for NETCDF4 output

#ifdef __linux
  ! Linux
  INTEGER, PARAMETER :: path_limit = 4096   ! PATH_MAX (Posix, limits.h)  
  INTEGER, PARAMETER :: name_limit =  255   ! NAME_PATH (Posix, limits.h)
#else
  ! MacOSX, AIX, Solaris
  INTEGER, PARAMETER :: path_limit = 1024   ! PATH_MAX (Posix, limits.h)  
  INTEGER, PARAMETER :: name_limit =  255   ! NAME_PATH (Posix, limits.h)
#endif
  INTEGER, PARAMETER :: exp_limit  = name_limit
  INTEGER, PARAMETER :: dir_limit  = path_limit-name_limit

  INTEGER :: syr, smo, sdy, shr, smn, sse   ! used for subjob procedure

  CHARACTER (len=exp_limit)  :: yexp

  CHARACTER (len=dir_limit)  :: yomdn
  CHARACTER (len=dir_limit)  :: ypath

  CHARACTER (len=path_limit) :: out_datapath  = ' '
  CHARACTER (len=name_limit) :: out_expname   = ' '

  INTEGER                    :: rerun_filetype = NETCDF
  INTEGER                    :: out_filetype   = GRIB
  INTEGER                    :: trac_filetype  = GRIB
  INTEGER                    :: out_ztype      = NONE

  CHARACTER (len=path_limit) :: standard_output_file
  CHARACTER (len=path_limit) :: standard_rerun_file

CONTAINS
!------------------------------------------------------------------------------
  SUBROUTINE compose_filenames

    INTEGER           :: i, iypt, ioff, fchour
    CHARACTER(len=19) :: time_label, suffix

    !-- Compose file names:

    !-- copy directory path first 

    standard_output_file = TRIM(out_datapath)
    iypt = MIN(LEN_TRIM(standard_output_file),dir_limit)

    !-- merge experiment label

    ioff = MAX(LEN_TRIM(out_expname), exp_limit)
    DO i = 1, ioff
      standard_output_file(iypt+i:iypt+i) = '_'
    ENDDO
    yexp = TRIM(out_expname)

    standard_output_file(iypt+1:iypt+ioff) = TRIM(yexp)

    !-- compose time label

    iypt = LEN_TRIM(standard_output_file)

    ypath = standard_output_file(1:INDEX(standard_output_file,'/',.TRUE.))
    WRITE (message_text, '(a,a,a,a,a)') &
         'Experiment: ', TRIM(yexp), ' - ',' data path: ', TRIM(ypath)
    CALL message('',message_text)

    ! open standard grib file
    
    IF (lnwp) THEN

       time_label = TRIM(str_date(FLT_NWP_YMD,start_date))
WRITE (0,*) '0) .... ', TRIM(time_label)
       CALL get_date_components(next_date,syr,smo,sdy,shr,smn,sse)

       standard_output_file(iypt+1:iypt+LEN_TRIM(time_label)+1) &
         = '_' // TRIM(time_label)

       fchour = get_forecast_hours()

       WRITE(suffix,'(i3.3)') fchour
       standard_output_file = TRIM(standard_output_file) // '+' // TRIM(suffix)

       WRITE (message_text,'(a,i3,a,a)') &
            'NWP forecast (', fchour, &
            'hr) GRIB output: ', &
            TRIM(standard_output_file(INDEX(standard_output_file,'/',.TRUE.)+1:))
       CALL message('',message_text)

    ELSE

       time_label = TRIM(str_date(timelabel_type,next_date))
       CALL get_date_components(next_date,syr,smo,sdy,shr,smn,sse)

       standard_output_file(iypt+1:iypt+LEN_TRIM(time_label)+1) &
         = '_' // TRIM(time_label)

       WRITE (message_text,'(a,a,a)') &
            'Standard GRIB output: ', &
            TRIM(standard_output_file(INDEX(standard_output_file,'/',.TRUE.)+1:))
       CALL message('',message_text)
    END IF

    CALL message('','')
    
  END SUBROUTINE compose_filenames
!------------------------------------------------------------------------------
!++mgs: added tt and ll arguments
  FUNCTION str_filter(vorlage,yr,mo,dy,hr,mi,se,nn,cc,tt,ll) RESULT(str1)

    CHARACTER(len=*)   ,INTENT(in)           :: vorlage
    INTEGER            ,INTENT(in) ,OPTIONAL :: yr, mo, dy, hr, mi, se, nn
    CHARACTER(len=*)   ,INTENT(in) ,OPTIONAL :: cc
    CHARACTER(len=*)   ,INTENT(in) ,OPTIONAL :: tt, ll   ! resolution as 'Tnnn', 'Lnn'
    CHARACTER(len=256)                       :: str1
    !---------------------------------------------------
    ! filter the string until all token are solved
    ! same token can be multiple used
    ! 
    ! Token    Replacement                       format
    !
    !  %Yx     yr (year)                          Ix.x
    !  %Mx     mo (month as number)               Ix.x
    !  %Ax     mo (month as short name 'JAN'...)  A3
    !  %Dx     dy (day)                           Ix.x
    !  %Hx     hr (hour)                          Ix.x
    !  %Ix     mi (minute)                        Ix.x
    !  %Sx     se (second)                        Ix.x
    !  %Nx     nn (a number)                      Ix.x
    !  %Cx     cc (a character string)            A
    !  %Tx     tt (spectral truncation as char)   A
    !  %Lx     ll (number of levels as char)      A
    !
    !  x is a digit (0...9). If x is 0 the minimum length
    !    required to write the integer number is used.
    !----------------------------------------------------

    CHARACTER(256) :: str2
    CHARACTER(1), PARAMETER :: token = '%'

    CHARACTER(3)   :: mon(12) = (/&
         'JAN','FEB','MAR','APR','MAY','JUN',&
         'JUL','AUG','SEP','OCT','NOV','DEC' /)

    CHARACTER(1)  :: str
    CHARACTER(64) :: transfer
    CHARACTER(20) :: form, cval
    INTEGER       :: pos_token, str_len, value

    str1 = TRIM(vorlage)
    DO
       ! find token position

       pos_token = INDEX(str1, token)
       IF (pos_token == 0) EXIT  ! all token expanded

       str_len  = LEN_TRIM(str1)
       IF (pos_token+2 > str_len) &
         CALL finish('str_filter','String length too large.')

       ! determine parameter to insert

       SELECT CASE(str1(pos_token+1:pos_token+1))    ! expand token and insert
       CASE('y','Y'); value = yr        !   year
       CASE('m','M'); value = mo        !   month as number
       CASE('d','D'); value = dy        !   day
       CASE('h','H'); value = hr        !   hour
       CASE('i','I'); value = mi        !   minute
       CASE('s','S'); value = se        !   second
       CASE('n','N'); value = nn        !   a number
       CASE('a','A'); value = mo        !   month as short name
       CASE('c','C'); value = 0         !   character string
       CASE('t','T'); value = 0         !   spectral truncation (character string)
       CASE('l','L'); value = 0         !   number of levels (character string)
       CASE default
          CALL finish('str_filter','Token not defined.')
       END SELECT

       ! determine format

       str = str1(pos_token+2:pos_token+2)

       ! determine length automatically if zero length is given

       IF (str == '0') THEN
         cval = ''
         WRITE (cval, *    ) value
         do
           if(cval(1:1)/=' ') exit
           cval = cval(2:)
         end do
         WRITE (str ,'(i1)') LEN_TRIM (cval)
       ENDIF

       SELECT CASE(str)
       CASE('1','2','3','4','5','6','7','8','9')
          WRITE(form,'(a2,a1,a1,a1,a1)') '(i',str,'.',str,')'
       CASE default
          CALL finish('str_filter','Token format not correct.')
       END SELECT
       SELECT CASE(str1(pos_token+1:pos_token+1))  ! expand token and insert
         CASE('a','A')                             !   month as short name
           str2 = str1(1:pos_token-1) // mon(mo) // str1(pos_token+3:str_len)
         CASE('c','C')                             ! string
           str2 = str1(1:pos_token-1) // TRIM(cc)// str1(pos_token+3:str_len)
         CASE('t','T')                             ! spectral truncation string
           str2 = str1(1:pos_token-1) // TRIM(tt)// str1(pos_token+3:str_len)
         CASE('l','L')                             ! nlev string
           str2 = str1(1:pos_token-1) // TRIM(ll)// str1(pos_token+3:str_len)
         CASE default                              !   integer value
           WRITE(transfer,form)   value
           str2 = str1(1:pos_token-1) // TRIM(transfer) //&
                  str1(pos_token+3:str_len)
       END SELECT
       str1 = TRIM(str2)

    END DO

  END FUNCTION str_filter
!------------------------------------------------------------------------------
  FUNCTION find_next_free_unit(istart,istop) RESULT(unit)
    INTEGER :: istart, istop, unit
    LOGICAL :: found, opened
    INTEGER :: i

    found = .FALSE.
    DO i=istart,istop
       INQUIRE(unit=i,opened=opened)
       IF (.NOT.opened) THEN
          unit = i
          found = .TRUE.
          EXIT
       END IF
    END DO

    IF (.NOT. found) THEN
       WRITE(message_text,'(a,i2.2,a,i2.2,a)') &
         'No unit in range <',istart,':',istop,'> free.'
       CALL finish('find_next_free_unit',message_text)
    END IF

  END FUNCTION find_next_free_unit
!------------------------------------------------------------------------------
END MODULE mo_filename
