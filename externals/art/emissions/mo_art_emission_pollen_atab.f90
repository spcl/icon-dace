! Description:!
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

MODULE mo_art_emission_pollen_atab
!=============================================================================
!
! Description:
!   Module containing routines to handle ATAB formatted data.
!
! Author: Pirmin Kaufmann, MeteoSwiss, Zurich, Switzerland.
!
! Current Code Owner: ditto.
! 
! Language: Fortran 95.
! 
! See svn log for history of modifications
!
!-----------------------------------------------------------------------------

! ICON
USE mo_kind,                          ONLY: wp 
USE mo_grid_config,                   ONLY: grid_sphere_radius
USE mo_art_config,                    ONLY: IART_PATH_LEN
USE mo_exception,                     ONLY: message, finish, message_text
USE mo_parallel_config,               ONLY: nproma, p_test_run
USE mo_math_constants,                ONLY: pi
USE mo_model_domain,                  ONLY: t_patch
USE mo_gnat_gridsearch,               ONLY: gnat_init_grid, gnat_destroy, t_gnat_tree, &
                                        &   gnat_query_containing_triangles,           &
                                        &   gnat_merge_distributed_queries, gk

! ART
USE mo_art_impl_constants,            ONLY: IART_VARNAMELEN

! Module variables
IMPLICIT NONE
PRIVATE ! Make all module variables inaccessible outside module by default
SAVE    ! Save module variables in case module is not used by main program
        ! (default in Fortran 2008+)

CHARACTER(LEN=*), PARAMETER :: routine = 'mo_art_emission_pollen_atab'  

PUBLIC :: t_art_all_stns, t_art_stn, t_art_pol_atab



!+----------------------------------------------------------------------------
! Public parameters
PUBLIC :: maxlen_head, maxlen_label, maxlen_data, i_undefined, r_undefined

! Public types
PUBLIC :: t_art_pollen_atab_header_info

! Public subroutines
PUBLIC :: art_pollen_read_atab
!-----------------------------------------------------------------------------

! Public functions 
PUBLIC:: ntoken, tolower


! Declarations for atab
! ============

! Parameter
INTEGER,PARAMETER   :: maxlen_head  =  32760    ! Maximal length of header line
INTEGER,PARAMETER   :: maxlen_label =     80    ! Maximal length of label
INTEGER,PARAMETER   :: maxlen_data  =  32760    ! Maximal length of data line
INTEGER,PARAMETER   :: i_undefined  = -HUGE(1)  ! Missing value, integer
REAL(wp),   PARAMETER   :: r_undefined  = -HUGE(1._wp) ! Missing value, real

! Type for header information
TYPE t_art_pollen_atab_header_info
   CHARACTER(LEN=maxlen_head) :: prod_type     = ''
   CHARACTER(LEN=maxlen_head) :: experiment    = ''
   CHARACTER(LEN=maxlen_head) :: description   = ''
   CHARACTER(LEN=maxlen_head) :: start_time    = ''
   CHARACTER(LEN=maxlen_head) :: end_time      = ''
   CHARACTER(LEN=maxlen_head) :: min_lead_time = ''
   CHARACTER(LEN=maxlen_head) :: max_lead_time = ''
   CHARACTER(LEN=maxlen_head) :: gen_process   = ''
   CHARACTER(LEN=maxlen_head) :: model_name    = ''
   CHARACTER(LEN=maxlen_head) :: model_version = ''
   REAL(wp)                   :: level         = r_undefined
   CHARACTER(LEN=maxlen_head) :: level_type    = ''
   REAL(wp)                   :: missing_value = r_undefined
   INTEGER                    :: txt_col_width = i_undefined
   INTEGER                    :: int_cols      = i_undefined
   INTEGER                    :: real_cols     = i_undefined
   INTEGER                    :: data_cols     = i_undefined
   INTEGER                    :: data_rows     = i_undefined
! Parameter-related keys
   CHARACTER(LEN=maxlen_label),POINTER :: descriptor(:) => NULL()
   CHARACTER(LEN=maxlen_label),POINTER :: PARAMETER(:)  => NULL()
   CHARACTER(LEN=maxlen_label),POINTER :: unit(:)       => NULL()
! Data column-related keys
   CHARACTER(LEN=maxlen_label),POINTER :: indicator(:)  => NULL()
   REAL,POINTER               :: latitude(:)    => NULL()
   REAL,POINTER               :: longitude(:)   => NULL()
   REAL,POINTER               :: elevation(:)   => NULL()
   REAL,POINTER               :: grid_i(:)      => NULL()
   REAL,POINTER               :: grid_j(:)      => NULL()
   REAL,POINTER               :: height(:)      => NULL()

   INTEGER,ALLOCATABLE        :: grid_idx(:)
   INTEGER,ALLOCATABLE        :: grid_blk(:)
END TYPE t_art_pollen_atab_header_info

TYPE t_art_stn
  CHARACTER(LEN=IART_VARNAMELEN) :: &
    &  id                      !< Name of source
  REAL(wp)                :: &
    &  lon,                  & !< longitude of the emission source
    &  lat,                  & !< latitude of the emission source
    &  saisl                   !< length of season at station
  INTEGER                 :: &
    &  ithis_nlocal_pts,     & !< local number of sources on 'subdomain' in parallel mode  
    &  tri_iidx_loc,         & !< location idx of source
    &  tri_iblk_loc            !< location blk of source
  LOGICAL                 :: &
    &  is_init
  CONTAINS
    PROCEDURE, PUBLIC :: init     => init_stn
    !PROCEDURE, PUBLIC :: print    => print_pntSrc_summary
END TYPE t_art_stn

TYPE t_art_all_stns
  TYPE(t_art_stn), ALLOCATABLE :: &
    &  p(:)                        !< List of all point sources
  INTEGER                     :: &
    &  n_stns                    !< number of active source scenarios
END TYPE t_art_all_stns

TYPE t_art_pol_atab
  TYPE(t_art_all_stns)                  :: &
    &  all_stns
  TYPE(t_art_pollen_atab_header_info)   :: &
    &  header_info
  INTEGER                               :: &
    &  n_days
  INTEGER,ALLOCATABLE                   :: &
    &  idays(:)
  REAL(wp),ALLOCATABLE                      :: &
    &  temp_array(:,:)
END TYPE t_art_pol_atab


! Settings
! ========
! Module variables
CHARACTER(LEN=*),PARAMETER    :: atab_type_string1 = 'ATAB' ! Normal type string
CHARACTER(LEN=*),PARAMETER    :: atab_type_string2 = 'XLS_TABLE'    ! Fieldextra
CHARACTER(LEN=*),PARAMETER    :: separator         = ACHAR(9)//' ,;!'
CHARACTER(LEN=*),PARAMETER    :: numeric           = separator//'0123456789.Ee+-'
CHARACTER(LEN=*),PARAMETER    :: comment_char      = '#!'
!-----------------------------------------------------------------------------

! Module procedures
! =================

CONTAINS

! Public functions
! ================

!+****************************************************************************
  INTEGER FUNCTION ntoken(string)
!=============================================================================
!
! Counts number of tokens in a string.
!
! ____________________________________________________________________________
! I/O   Name            Type    Description
! ____________________________________________________________________________
! I     i               C(*)    String
! O     ntoken          i       Number of tokens in string
!
! Copyright (c) 1997-1998 Pirmin Kaufmann
!
! Version history:
! 1.1   Use VERIFY, SCAN intrinsic instead of tabinln (faster, more reliable)
! 1.0   Basic version
!-----------------------------------------------------------------------------
IMPLICIT NONE
! Dummy arguments
CHARACTER(LEN=*),INTENT(IN) :: string

! Local parameters
CHARACTER(LEN=*), PARAMETER :: tab = ACHAR(9)
CHARACTER(LEN=*), PARAMETER :: whitespace = ACHAR(0)//ACHAR(9)//' '
CHARACTER(LEN=*), PARAMETER :: quotes = "'"//'"'
CHARACTER(LEN=*), PARAMETER :: separators = whitespace//';,'
! Lokal variables
INTEGER :: ilast,ii

ntoken = 0
ii = 1

DO
   ilast = ii
   ii = VERIFY(string(ii:),separators) + ii - 1
   IF (ii >= ilast) THEN
      ntoken = ntoken + 1
   ELSE
      ! all remaining characters are separators or no remaining characters
      EXIT
   END IF
   ilast = ii
   IF (VERIFY(string(ii:ii),quotes) <= 0 ) THEN
      ii = SCAN(string(ii+1:),quotes) + ii
   END IF
   ii = SCAN(string(ii:),separators) + ii - 1
   IF (ii < ilast) THEN
      ! no end quotes or no more separators found
      EXIT
   END IF
END DO

END FUNCTION ntoken


!+****************************************************************************
FUNCTION tolower(string)
!=============================================================================
!
! Returns UPPERCASE characters as lowercase.
!
! Requires explicit interface in calling program unit.
!
! ____________________________________________________________________________
! I/O   Name            Type    Description
! ____________________________________________________________________________
! I     string          C(*)    String
!
! Copyright (c) 1994-1998 Pirmin Kaufmann
!
! Version history:
! 1.0   Basic version
!-----------------------------------------------------------------------------
IMPLICIT NONE
! Dummy arguments
CHARACTER(LEN=*),INTENT(IN)  :: string
CHARACTER(LEN=LEN(string))   :: tolower
! Local variables
INTEGER   :: ii
CHARACTER :: ch

tolower = string
DO ii = 1, LEN_TRIM(string)
    ch = string(ii:ii)
    IF ('A' <= ch .AND. ch <= 'Z') THEN
       tolower(ii:ii) = ACHAR(IACHAR(ch) + 32)
    END IF
END DO

END FUNCTION tolower


! Public subroutines
! ==================

!+****************************************************************************
SUBROUTINE art_pollen_read_atab(pollname,cart_input_folder,pol_atab,p_patch)
!=============================================================================
!
! Read a data table from a file.
!
! Variables:
! ____________________________________________________________________________
! I/O   Name            Type            Description
! ____________________________________________________________________________
! I     file            char(*)         input filename
! I     pollname        char(*)         pollen name for select case (alt to file)
! I     cart_input_fld  char(*)         input base folder (alt to file)
! O     nrow            int             number of data rows
! O     lablen          int             Width of text label column
! O     icol            int             number of integer label columns
! O     fcol            int             number of real label columns
! O     ncol            int             number of data columns
! O     rlabel          char(*)         array of row labels (nrow)
! O     iarray          int(*,*)        array of integer labels (icol, nrow)
! O     farray          flt(*,*)        array of real labels (fcol, nrow)
! O     array           flt(*,*)        array of data (ncol, nrow)
! O     rclabel         char(*)         label of text label column
! O     ilabel          char(*)         labels of integer label columns (icol)
! O     flabel          char(*)         labels of real label columns (fcol)
! O     alabel          char(*)         labels of data columns (ncol)
! O     header_info     type            extract header information
! O     iostat          int             error flag (<0: warning, >0: error)
! O     iomsg           char(*)         error message
!-----------------------------------------------------------------------------

  ! Dummy arguments
  CHARACTER(LEN=*),                    INTENT(in)    :: pollname
  CHARACTER(LEN=*),                    INTENT(in)    :: cart_input_folder
  TYPE(t_art_pol_atab),                INTENT(inout) :: pol_atab
  INTEGER,                             ALLOCATABLE   :: iarray(:,:)
  REAL(wp),                            ALLOCATABLE   :: array(:,:)
  TYPE(t_art_pollen_atab_header_info)                :: header_info
  TYPE(t_patch),                       INTENT(in)    :: p_patch                
                                                        !< Patch on which computation is performed

  INTEGER                                            :: iostat
  CHARACTER(LEN=100)                                 :: iomsg

  ! Local variables
  CHARACTER(LEN=IART_PATH_LEN)                       :: file   !len is adjustable

  CHARACTER(LEN=*),                    PARAMETER     :: tag = 'art_read_read_atab'
  CHARACTER(LEN=*),                    PARAMETER     :: errmsg = ' not defined in file header.'

  CHARACTER(LEN=maxlen_label),         POINTER       :: rlabel(:)
  CHARACTER(LEN=maxlen_label),         POINTER       :: ilabel(:)
  CHARACTER(LEN=maxlen_label),         POINTER       :: flabel(:)
  CHARACTER(LEN=maxlen_label),         POINTER       :: alabel(:)

  REAL(wp),                            POINTER       :: farray(:,:)
  REAL(wp),                            ALLOCATABLE   :: arraytmp(:,:)

  INTEGER :: lablen
  INTEGER :: icol
  INTEGER :: fcol
  INTEGER :: ii,stn
  INTEGER :: icolon,nrow,ncol
  INTEGER :: headerlines
  INTEGER :: lun
  INTEGER :: i_comment_char
  INTEGER :: descr_count, param_count, unit_count

  CHARACTER(LEN=maxlen_label) :: rclabel
  CHARACTER(LEN=maxlen_head)  :: key
  CHARACTER(LEN=maxlen_head)  :: keyval
  CHARACTER(LEN=maxlen_head)  :: value
  CHARACTER(LEN=maxlen_data)  :: headline
  CHARACTER(LEN=maxlen_data)  :: cleanline
  CHARACTER(LEN=maxlen_data)  :: dataline

  NULLIFY(ilabel)
  NULLIFY(flabel)
  NULLIFY(alabel)
  NULLIFY(rlabel)
  NULLIFY(farray)

  ! Initialize
  ! ----------
  CALL dealloc_header_info(header_info)
  headline = ''
  ! Number of header lines
  HeaderLines  = 0
  ! Generate file name (if pollname is given)
  SELECT CASE(pollname)
    CASE('pollbetu')
      file = TRIM(cart_input_folder)//'/norm.T.stns.betu.atab'
    CASE('pollalnu')
      file = TRIM(cart_input_folder)//'/norm.T.stns.alnu.atab'
    CASE('pollcory')
      file = TRIM(cart_input_folder)//'/norm.T.stns.cory.atab'
    CASE DEFAULT
      ! no data table for other pollen species
      RETURN
  END SELECT

  WRITE (message_text,*) 'ART: Opening pollen ATAB file ', TRIM(file)
  CALL message(TRIM(routine)//':art_pollen_read_atab', message_text)

  ! Open file
  lun = 17
  OPEN(UNIT = lun, &
       FILE = file, &
       STATUS = 'OLD', &
       ACTION = 'READ', &
       POSITION = 'REWIND', &
       RECL = MAX(maxlen_data,255), &
       IOSTAT = iostat, &
       IOMSG  = iomsg)

  IF (iostat /= 0) THEN
    CALL finish(TRIM(routine)//':art_pollen_read_atab','Could not open pollen ATAB file '//TRIM(file)//'.')
  ENDIF


  ! Read header
  read_header_line: DO
     ! Read with '(a)' format, * format only reads until first separator
     READ(lun,'(a)') headline
     HeaderLines = HeaderLines + 1
     IF (HeaderLines == 1) THEN 
        IF (headline(:LEN_TRIM(atab_type_string1)) == atab_type_string1 .OR. &
            headline(:LEN_TRIM(atab_type_string2)) == atab_type_string2) THEN
           ! Skip to next line
           CYCLE
        ELSE
           ! Skip header processing
           EXIT
        END IF
     END IF
     ! Remove comments
     !i_comment_char = SCAN(headline,comment_char)
     i_comment_char = SCAN(headline,'#!')
     IF (i_comment_char > 0) THEN
        cleanline = headline(:i_comment_char-1)
     ELSE
        cleanline = headline
     END IF
     ! Check if file type indicator or empty line (after comment removal)
     IF (.NOT. (HeaderLines == 1 .OR. TRIM(cleanline) == '')) THEN
        ! Regular line
        ! Split line at first colon
        icolon = INDEX(cleanline,':')
        IF (icolon > 0) THEN
           ! extract key and value
           key = cleanline(:icolon-1)
           ! key: substitute underlines by whitespace, use lower case
           keyval = tolower(key)
           DO
              ii = INDEX(keyval,'_')
              IF (ii <= 0) EXIT
              keyval(ii:ii) = ' '
           END DO
           ! extract value
           value = TRIM(ADJUSTL(cleanline(icolon+1:)))
           IF (LEN_TRIM(value) > 0) THEN
              ! Read value into header_info structure
              SELECT CASE (keyval)
              CASE('type of product') 
                 header_info%prod_type = value
              CASE('experiment') 
                 header_info%experiment = value
              CASE('description')
                 header_info%description = value
              CASE('start time')
                 header_info%start_time = value
              CASE('observation start')  ! obsolete
                 header_info%start_time = value
              CASE('end time')
                 header_info%end_time = value
              CASE('observation end')    ! obsolete
                 header_info%end_time = value
              CASE('minimum lead time')
                 header_info%min_lead_time = value
              CASE('maximum lead time')
                 header_info%max_lead_time = value
              CASE('generating process')
                 header_info%gen_process = value
              CASE('model name')
                 header_info%model_name = value
              CASE('model version')
                 header_info%model_version = value
              CASE('level')
                 READ(value,*) header_info%level
              CASE('level type')
                 header_info%level_type = value
              CASE('missing value code')
                 READ(value,*) header_info%missing_value
              CASE('width of text label column')
                 READ(value,*) lablen
              CASE('number of integer label columns')
                 READ(value,*) icol
              CASE('number of real label columns')
                 READ(value,*) fcol
              CASE('number of data columns')
                 READ(value,*) ncol
              CASE('number of stations') ! obsolete
                 READ(value,*) ncol
              CASE('number of data rows')
                 READ(value,*) nrow
              ! Parameter-related keys
              CASE('descriptor')
                 descr_count = ntoken(value)
                 ALLOCATE(header_info%descriptor(descr_count))
                 READ(value,*) header_info%descriptor
                 !strtoken(value,/extract,count=descr_count)
              CASE('parameter')
                 param_count = ntoken(value)
                 ALLOCATE(header_info%parameter(param_count))
                 READ(value,*) header_info%parameter
              CASE('unit')
                 unit_count = ntoken(value)
                 IF (unit_count > 1 .AND. MAX(descr_count,param_count) == 1) THEN
                    ! Backward compatibility (unit with spaces
                    ! but without quotes) for single param files
                    ALLOCATE(header_info%unit(1))
                    header_info%unit = value
                 ELSE
                    ALLOCATE(header_info%unit(unit_count))
                    READ(value,*) header_info%unit
                 END IF
              ! Data-column-related keys
              CASE('indicator')
                 ALLOCATE(header_info%indicator(ncol))
                 READ(value,*) header_info%indicator
              CASE('latitude')
                 ALLOCATE(header_info%latitude(ncol))
                 READ(value,*) header_info%latitude
              CASE('longitude')
                 ALLOCATE(header_info%longitude(ncol))
                 READ(value,*) header_info%longitude
              CASE('elevation')
                 ALLOCATE(header_info%elevation(ncol))
                 READ(value,*) header_info%elevation
              CASE('grid i')
                 ALLOCATE(header_info%grid_i(ncol))
                 READ(value,*) header_info%grid_i
              CASE('grid j')
                 ALLOCATE(header_info%grid_j(ncol))
                 READ(value,*) header_info%grid_j
              CASE('height')
                 ALLOCATE(header_info%height(ncol))
                 READ(value,*) header_info%height
              CASE default
                 iostat = 1
                 iomsg  = 'Unknown ATAB keyword: '//TRIM(key)
                 RETURN
              END SELECT
           ELSE
              ! No value after keyword
              iostat = -1
              iomsg  = 'Value missing after keyword: '//TRIM(key)
           END IF
        ELSE
           ! Non-empty line contains no colon, not a valid header line
           IF (LEN_TRIM(cleanline) > 0) EXIT
        END IF
     END IF ! Regular line
  END DO read_header_line

  ! Fill in header_info argument
  header_info%txt_col_width =  lablen
  header_info%int_cols      =  icol
  header_info%real_cols     =  fcol
  header_info%data_cols     =  ncol
  header_info%data_rows     =  nrow

  ! Check if table dimensions have been defined
  IF (lablen < 0 .OR. icol < 0 .OR. ncol < 0 .OR. nrow < 0) THEN
     ! Error
     iostat = 1
     ! At least one dimension is missing or text label width is missing
     ! File must be ATAB format, dimensions must be present
     IF (lablen < 0) THEN
        iomsg = 'Width of text label column'//errmsg
     END IF
     IF (icol   < 0) THEN
        iomsg = 'Number of integer label columns'//errmsg
     END IF
     IF (ncol   < 0) THEN
        iomsg = 'Number of data columns'//errmsg
     END IF
     IF (nrow   < 0) THEN
        iomsg = 'Number of data rows'//errmsg
     END IF
     ! Close file
     CLOSE(lun)
     RETURN
  END IF

  ! Allocate storage
  ! If allocated with incorrect size/shape, deallocate first
  ! to avoid memory leaks (allocate on existing pointer creates new
  ! pointer, breaking link to existing storage without deallocating it)
  ! Pointer deallocation also nullifies pointer

  ! ilabel
  IF (ASSOCIATED(ilabel)) THEN
     IF (SIZE(ilabel) /= icol) THEN
        DEALLOCATE(ilabel)
        ALLOCATE(ilabel(icol))
     END IF
  ELSE
     ALLOCATE(ilabel(icol))
  END IF
  ! flabel
  IF (ASSOCIATED(flabel)) THEN
     IF (SIZE(flabel) /= fcol) THEN
         DEALLOCATE(flabel)
        ALLOCATE(flabel(fcol))
     END IF
  ELSE
     ALLOCATE(flabel(fcol))
  END IF
  ! alabel
  IF (ASSOCIATED(alabel)) THEN
     IF (SIZE(alabel) /= ncol) THEN
        DEALLOCATE(alabel)
        ALLOCATE(alabel(ncol))
     END IF
  ELSE
     ALLOCATE(alabel(ncol))
  END IF
  ! rlabel
  IF (ASSOCIATED(rlabel)) THEN
     IF (SIZE(rlabel) /= nrow) THEN
        DEALLOCATE(rlabel)
        ALLOCATE(rlabel(nrow))
     END IF
  ELSE
     ALLOCATE(rlabel(nrow))
  END IF
  ! iarray
  IF (ALLOCATED(iarray)) THEN
     IF (ANY(SHAPE(iarray) /= (/ icol,nrow /))) THEN
        DEALLOCATE(iarray)
        ALLOCATE(iarray(icol,nrow))
     END IF
  ELSE
     ALLOCATE(iarray(icol,nrow))
  END IF
  ! farray
  IF (ASSOCIATED(farray)) THEN
     IF (ANY(SHAPE(farray) /= (/ fcol,nrow /))) THEN
        DEALLOCATE(farray)
        ALLOCATE(farray(fcol,nrow))
     END IF
  ELSE
     ALLOCATE(farray(fcol,nrow))
  END IF
  ! array
  IF (ALLOCATED(array)) THEN
     IF (ANY(SHAPE(array) /= (/ ncol,nrow /))) THEN
        DEALLOCATE(array)
        ALLOCATE(array(ncol,nrow))
     END IF
  ELSE
     ALLOCATE(array(ncol,nrow))
  END IF
  IF (ALLOCATED(arraytmp)) THEN
     IF (ANY(SHAPE(arraytmp) /= (/ ncol,nrow /))) THEN
        DEALLOCATE(arraytmp)
        ALLOCATE(arraytmp(ncol,nrow))
     END IF
  ELSE
     ALLOCATE(arraytmp(ncol,nrow))
  END IF

  ! Extract labels
  IF (lablen > 0) THEN
     rclabel = TRIM(ADJUSTL(headline(1:lablen)))
  ELSE
     rclabel = ''
  END IF
! READ(headline(lablen+1:),*) ilabel, flabel, alabel

  ! Read data table
  IF (lablen <= 0) THEN
     ! No text labels
     DO ii = 1, nrow
        READ(lun,*) iarray(:,ii), farray(:,ii), arraytmp(:,ii)
     END DO
  ELSE
     ! Text labels present
     DO ii = 1, nrow
        READ(lun,'(a)') dataline
        rlabel(ii) = dataline(:lablen)
        READ(dataline(lablen+1:),*) iarray(:,ii), farray(:,ii), arraytmp(:,ii)
     END DO
  END IF

  ! Close file
  CLOSE(lun)
  
  ALLOCATE(pol_atab%temp_array(ncol,nrow))
  ALLOCATE(pol_atab%idays(nrow))
  pol_atab%temp_array=REAL(arraytmp,wp)
  pol_atab%n_days=nrow
  pol_atab%idays=iarray(1,:)

  ! Fill t_art_all_stns
  pol_atab%all_stns%n_stns=ncol
  IF(.NOT.ALLOCATED(pol_atab%all_stns%p)) THEN
    ALLOCATE(pol_atab%all_stns%p(ncol))
    DO stn=1,pol_atab%all_stns%n_stns
      pol_atab%all_stns%p(stn)%is_init=.FALSE.
    END DO
  END IF
! ALLOCATE(all_stns%p(ncol))
  DO stn=1,pol_atab%all_stns%n_stns
    pol_atab%all_stns%p(stn)%id=TRIM(header_info%indicator(stn))
    pol_atab%all_stns%p(stn)%lon=REAL(header_info%longitude(stn),wp)
    pol_atab%all_stns%p(stn)%lat=REAL(header_info%latitude(stn),wp)
    IF(.NOT.pol_atab%all_stns%p(stn)%is_init) then
!WRITE(0,*)'vorher',stn,pol_atab%all_stns%p(stn)%is_init
    CALL pol_atab%all_stns%p(stn)%init(p_patch)
!WRITE(0,*)'nachher',stn,pol_atab%all_stns%p(stn)%is_init
    ENDIF
  END DO

  pol_atab%header_info=header_info

END SUBROUTINE art_pollen_read_atab

!+****************************************************************************
SUBROUTINE dealloc_header_info(header_info)
!=============================================================================
!
! Deallocate header_info if allocated, to avoid memory leaks:
! Allocate on existing pointer creates new
! pointer, breaking link to existing storage without deallocating it.
! Pointer deallocation also nullifies pointer
!
! Variables:
! ____________________________________________________________________________
! I/O   Name            Type                            Description
! ____________________________________________________________________________
! O     header_info     t_art_pollen_atab_header_info   header information
!
! Copyright (c) 2006 Pirmin Kaufmann
!-----------------------------------------------------------------------------
  ! Dummy arguments
  TYPE(t_art_pollen_atab_header_info) :: header_info

  ! Pointer deallocation also nullifies pointer
  IF (ASSOCIATED(header_info%indicator)) DEALLOCATE(header_info%indicator)
  IF (ASSOCIATED(header_info%latitude))  DEALLOCATE(header_info%latitude)
  IF (ASSOCIATED(header_info%longitude)) DEALLOCATE(header_info%longitude)
  IF (ASSOCIATED(header_info%elevation)) DEALLOCATE(header_info%elevation)
  IF (ASSOCIATED(header_info%grid_i))    DEALLOCATE(header_info%grid_i)
  IF (ASSOCIATED(header_info%grid_j))    DEALLOCATE(header_info%grid_j)
  IF (ASSOCIATED(header_info%height))    DEALLOCATE(header_info%height)

END SUBROUTINE dealloc_header_info

!+****************************************************************************
SUBROUTINE init_stn ( this_stn, p_patch)
  CLASS(t_art_stn), INTENT(inout) :: &
    &  this_stn
  TYPE(t_patch), INTENT(in)  :: &
    &  p_patch                 !< Patch on which computation is performed
  !TYPE(t_art_all_stns),INTENT(INOUT)      :: &
  !  &  all_stns             !< header information in atab-file
! local variables                    
  TYPE(t_gnat_tree)          :: &
    &  gnat                   !< Grid search object
  INTEGER                    ::   &
    &  stn,                       & !< loop index for stations
    &  gnat_nblks, gnat_npromz,   & !< calculated for call to GNAT
    &  gnat_jc,gnat_jb              !< coordinates for GNAT

  REAL(wp),ALLOCATABLE         :: &
    &  lat_idx(:,:),lon_idx(:,:)

  REAL(wp)                     :: &
    &  stns_lon_deg,              & !< coordinate of synop station (in degree)
    &  stns_lat_deg                 !< coordinate of synop station (in degree)

  REAL(gk),ALLOCATABLE         :: &
    &  in_points(:,:,:),          & !< geographical locations (GNAT)
    &  min_dist(:,:)             !& !< minimal distance (GNAT)
   ! &  lat_idx(:,:),              &
   ! &  lon_idx(:,:)

  INTEGER, ALLOCATABLE         :: &
    &  tri_idx(:,:,:),            &
    &  owner(:)                 !< rank of sender PE for each station (GNAT)
  

  ! Initialize local fields with input data 
  ! in atab: lon values are -180>lon>180, maybe icon: 0>lon>360 ?
  stns_lon_deg = REAL(this_stn%lon,wp) + 360
  stns_lon_deg = MOD(stns_lon_deg,360._wp)
  stns_lat_deg = REAL(this_stn%lat,wp)

! WRITE(0,*) 'Start: gnat'
  gnat_nblks  = 1/nproma+1
  gnat_npromz = 1-nproma*(gnat_nblks-1)

  ! --- Allocate strcutures demanded by GNAT
  ALLOCATE( tri_idx(2,nproma,gnat_nblks)   )
  ALLOCATE( owner(1)                  )
  ALLOCATE( in_points(nproma,gnat_nblks,2) )
  ALLOCATE( min_dist(nproma,gnat_nblks)    )
  ALLOCATE( lat_idx(nproma,gnat_nblks)     )
  ALLOCATE( lon_idx(nproma,gnat_nblks)     )

  tri_idx   = 0
  owner     = 0
  in_points = 0._wp
  min_dist  = 0._wp
  lat_idx   = 0._wp
  lon_idx   = 0._wp

  ! --- Save lon/lat into structures demanded by GNAT
  gnat_jc = 1
  gnat_jb = 1

  lat_idx(gnat_jc,gnat_jb) = stns_lat_deg
  lon_idx(gnat_jc,gnat_jb) = stns_lon_deg

  in_points(gnat_jc,gnat_jb,1) = lon_idx(gnat_jc,gnat_jb) * pi/180._wp
  in_points(gnat_jc,gnat_jb,2) = lat_idx(gnat_jc,gnat_jb) * pi/180._wp

  ! --- Build GNAT data structure
  CALL gnat_init_grid(gnat, p_patch)

  ! --- Perform proximity query
  CALL gnat_query_containing_triangles(gnat, p_patch, in_points(:,:,:),               &
    &                                  nproma, gnat_nblks, gnat_npromz,               &
    &                                  grid_sphere_radius,p_test_run,tri_idx(:,:,:),  &
    &                                  min_dist(:,:))
  CALL gnat_merge_distributed_queries(p_patch, 1, nproma, gnat_nblks, min_dist,       &
    &                                 tri_idx(:,:,:), in_points(:,:,:),               &
    &                                 owner(:), this_stn%ithis_nlocal_pts)

  ! --- Cleanup GNAT
  CALL gnat_destroy(gnat)

  ! --- Get gathered variables at station 
  gnat_jc = 0
  gnat_jb = 1
  DO stn=1, this_stn%ithis_nlocal_pts
    gnat_jc = gnat_jc + 1
    IF (gnat_jc>nproma) THEN
      gnat_jc = 1
      gnat_jb = gnat_jb + 1
    ENDIF

    this_stn%tri_iidx_loc = tri_idx(1,gnat_jc,gnat_jb)
    this_stn%tri_iblk_loc = tri_idx(2,gnat_jc,gnat_jb)
    this_stn%is_init = .TRUE.
!   WRITE(0,*)'im ini',  this_stn%is_init,tri_idx(1,gnat_jc,gnat_jb), tri_idx(2,gnat_jc,gnat_jb)

  END DO

  DEALLOCATE( tri_idx    )
  DEALLOCATE( owner      )
  DEALLOCATE( in_points  )
  DEALLOCATE( min_dist   )
  DEALLOCATE( lat_idx    )
  DEALLOCATE( lon_idx    )

END SUBROUTINE init_stn
!=============================================================================

END MODULE mo_art_emission_pollen_atab
