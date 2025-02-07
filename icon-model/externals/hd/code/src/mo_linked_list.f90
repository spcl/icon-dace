! mo_linked_list.f90 - Managing of linked lists
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!_________________________________________

MODULE mo_linked_list

  ! This module holds the definitions of the data type and basic
  ! operators used for the memory buffer linked lists. The linked lists
  ! serve two purposes:
  !
  ! 1) Maintain a list of all 3D and 2D atmospheric and surface fields
  ! allocated permanently in the model.
  !
  ! 2) Provide information for the postprocessing.
  !
  ! Authors:
  !
  ! Luis Kornblueh, MPI,             original code
  ! Andreas Rhodin, MPI, April 2001, extended and documented
  ! I. Kirchner,    MPI, August 2002  lpout flag
  ! R. Schnur,      MPI, January 2003 addes grid representation LAND
  !
  !  This routine originates (year 2014) from MPI-ESM, the Earth System Model of the 
  !  Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
  !  Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
  !  version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
  !  doi: 10.1029/2018MS001400.
  !------------------------------------------------------------------------------
  !
  ! Modules used
  !

  USE mo_kind,         ONLY: dp, i8
  USE mo_exception,    ONLY: message_text, finish, message
  USE mo_util_string,  ONLY: separator
  USE mo_filename,     ONLY: NETCDF, GRIB, NETCDF2, NETCDF4, out_filetype, &
                             NONE, SZIP, ZIP, out_ztype
  USE mo_netcdf,       ONLY: BELOWSUR, SURFACE, ABOVESUR2, ABOVESUR10, &
                             HYBRID, HYBRID_H, TILES, SOILLEV, ROOTZONES, CANOPY

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  !
  ! Public entities
  !

  PRIVATE

  ! Valid values for the representation flag 'repr' are as follows:

  PUBLIC :: GAUSSIAN, FOURIER, SPECTRAL, HEXAGONAL, GRIDPOINT, LAND

  ! Values for 'level_type'

  PUBLIC :: BELOWSUR, SURFACE, ABOVESUR2, ABOVESUR10, HYBRID, HYBRID_H, TILES, SOILLEV, ROOTZONES, CANOPY

  ! Valid values for file types

  PUBLIC :: GRIB, NETCDF, NETCDF2, NETCDF4

  ! Valid values for compress types

  PUBLIC :: NONE, SZIP, ZIP

  ! Value for any unknown (not yet initialised) parameters

  PUBLIC :: UNKNOWN

  PUBLIC :: memory_type         ! data type to hold field list entry

  PUBLIC :: memory_info         ! data type to hold meta data on list entry
  PUBLIC :: default_info        ! empty (default) meta data list entry

  PUBLIC :: t_stream            ! anchor for a whole list
  PUBLIC :: default_stream      ! default stream header

  PUBLIC :: list_element        ! entry in the linked list

  PUBLIC :: construct_list      ! construct an (empty) list
  PUBLIC :: destruct_list       ! clean up the list

  PUBLIC :: add_list_element    ! add an element to the list
  PUBLIC :: find_list_element   ! find an element in the list
  PUBLIC :: remove_list_element ! remove one element from the list

  PUBLIC :: print_linked_list
  PUBLIC :: print_sinfo_list
  PUBLIC :: print_stream        ! print output stream
  PUBLIC :: print_memory_info   ! print stream element info
  !------------------------------------------------------------------------------
  ! some constants, may be changed

  INTEGER ,PARAMETER :: maxlev = 5 ! maximum number of special level values 

  !------------------------------------------------------------------------------
  ! TYPE memory_info holds the metainformation for an entry in the memory
  ! buffer list:
  ! 
  ! Each entry is identified by a unique name 'name'. For convenience the
  ! 'units', the 'long name' and the 'molecular weight' may be given.
  ! 
  ! If 'alloc' is .true., the list entry has been used to allocate memory
  ! for the field. In this case the pointer 'ptr' in type 'memory_type'
  ! must be deallocated if the list is deleted. If 'alloc' is
  ! .false. 'ptr' is merely a reference to a field allocated elsethere.
  ! 'ptr' always references a 4 dimensional array allthough less
  ! dimensions may be actually used as indicated by 'ndim'.  'dim(1)' to
  ! 'dim(4)' provide the dimensins allocated (or referenced) on the local
  ! processor. 'gdim' are the dimensions of the global field.
  ! 'repr' indicates the representation of the entry (gridpoint spectral..)
  ! 
  ! If 'laccu' is .true. the field is divided by the output interval
  ! before output, and set to 'reset' afterwards. 'reset' is zero
  ! for accumulated fields and large or small for minimum and maximum
  ! values respectively.
  ! 
  ! If 'lmiss' is .true. the field is using missing value and is initialzed
  ! by 'missval'.
  ! 
  ! 'lrerun' indicates if the field must go into the restart
  ! file. 'restart_read' indicates if it was actually found in the restart
  ! file. A prefix 'restart_pref' may be added to the variable name in the
  ! restart file to avoid naming conflicts (used for tracer names).
  ! 
  ! Some additional information is provided for GRIB and NetCDF processing.

  TYPE memory_info                      ! meta data type
    
    !
    SEQUENCE
    !
    ! Content of the field
    !
    CHARACTER(len= 64) :: name          ! variable name
    CHARACTER(len= 64) :: units         ! units
    CHARACTER(len=128) :: longname      ! long name
    !
    ! Memory buffer information
    !
    INTEGER            :: gdim(4)       ! global dimensions of variable
    INTEGER            :: dim (4)       ! local dimensions (lon,lat)
    INTEGER            :: dima(4)       ! allocated dimensions (nproma,ngpblks)
    LOGICAL            :: lreg          ! true for dim==dima
    INTEGER            :: ndim          ! rank of variable
    INTEGER            :: klev          ! number of vertical levels
    INTEGER            :: iklev         ! vertical level index
    LOGICAL            :: alloc         ! pointer must be dealloc. by destruct
    INTEGER            :: repr          ! representation (gridpoint,spectral..)
    !
    ! General postprocessing information
    !
    LOGICAL            :: lpost         ! field goes into output stream
    LOGICAL            :: laccu         ! accumulation flag
    LOGICAL            :: lmiss         ! missing value flag
    REAL(dp)           :: missval       ! missing value
    REAL(dp)           :: reset         ! reset value for accumulated fields
    LOGICAL            :: lrerun        ! restart file flag
    LOGICAL            :: contnorest    ! continue if not in restart file     
    LOGICAL            :: restart_read  ! field has been set from restart file
    CHARACTER(len=3)   :: restart_pref  ! prefix for name in restart file
    !
    ! GRIB output information
    !
    INTEGER            :: gribtable     ! gribcode table number
    INTEGER            :: gribcode      ! gribcode number
    INTEGER            :: gribbits      ! number of bits used for GRIB encoding
    INTEGER            :: levelindx     ! HYBRID, ABOVESUR, SURFACE, BELOWSUR,
                                       ! TILES, SOILLEV, ROOTZONES, CANOPY
    !
    ! NETCDF output information
    !
    INTEGER            :: tracidx       ! tracer index
    INTEGER            :: IO_var_indx(4)! index to dimension table
    INTEGER            :: IO_var_id     ! NETCDF id for internal use (restart)
    INTEGER            :: IO_var_stid   ! NETCDF id for internal use (stream)
    CHARACTER(len=128) :: IO_name       ! name specifier for NETCDF file
    CHARACTER(len= 32) :: IO_unit       ! unit specifier for NETCDF file
    !
    ! Output information
    !
    INTEGER            :: gridID
    INTEGER            :: zaxisID
    !
    ! I/O server to work with
    !
    INTEGER           :: IO_comm_indx    ! index of mapped comm. in p_io_comm. 
    INTEGER           :: IO_stream_indx  ! index of stream this is member
    ! 
    ! Mark last element of derived type to allow proper usage later 
    ! in derived MPI datatypes 
    INTEGER            :: last          ! just dummy

  END TYPE memory_info
  !------------------------------------------------------------------------------
  ! Valid values for the representation flag 'repr' are as follows:

  INTEGER ,PARAMETER :: UNKNOWN   = -HUGE(0)
  INTEGER ,PARAMETER :: GAUSSIAN  = 1
  INTEGER ,PARAMETER :: FOURIER   = 2
  INTEGER ,PARAMETER :: SPECTRAL  = 3
  INTEGER ,PARAMETER :: HEXAGONAL = 4
  INTEGER ,PARAMETER :: LAND      = 5
  INTEGER ,PARAMETER :: GRIDPOINT = GAUSSIAN ! subject to change
  !------------------------------------------------------------------------------
  ! default (empty) meta data entry:

  TYPE(memory_info), PARAMETER :: default_info = &
       memory_info( &
        ''        , &     ! name
        ''        , &     ! unit
        ''        , &     ! long name
       UNKNOWN    , &     ! global dimensions
       UNKNOWN    , &     ! local dimensions (lon,lat)
       UNKNOWN    , &     ! allocated dimensions (nproma,ngpblks)
       .FALSE.    , &     ! true for dim==dima
       UNKNOWN    , &     ! rank
       UNKNOWN    , &     ! number of vertical levels
       UNKNOWN    , &     ! vertical level index
       .FALSE.    , &     ! pointer allocated flag
       GAUSSIAN   , &     ! representation
       .FALSE.    , &     ! print flag
       .FALSE.    , &     ! accumulation flag
       .FALSE.    , &     ! missing value flag
       -9.e+33_dp , &     ! missing value
       0._dp      , &     ! reset value for accumulated fields
       .FALSE.    , &     ! lrerun file flag
       .FALSE.    , &     ! continue if not in restart file
       .FALSE.    , &     ! read from restart flag
       ''         , &     ! prefix to name in restart file
       0          , &     ! GRIB table number
       0          , &     ! GRIB code number
       16         , &     ! bits used for GRIB encoding
       UNKNOWN    , &     ! level type
       0          , &     ! tracer index
       UNKNOWN,     &     ! dimension table entry indices
       0,           &     ! NETCDF id for internal use (restart)
       0,           &     ! NETCDF id for internal use (stream)
       'undefined', &     ! name specifier for NETCDF file
       'undefined', &     ! unit specifier for NETCDF file
       -1,          &     ! gridID
       -1,          &     ! zaxisID
        0,          &     ! I/O server to work with
        0,          &     ! I/O stream this belongs to
       -1)                ! last - dummy

  !------------------------------------------------------------------------------
  ! Type 'memory_type' holds the pointer to the field as well as the meta
  ! associated information.

  TYPE memory_type                         ! linked list entry type
    !
    SEQUENCE
    ! 
    REAL(dp), POINTER   :: ptr (:,:,:,:)  ! pointer to 3D-field
    TYPE (memory_info)  :: info           ! meta data for this entry
  END TYPE memory_type
  !------------------------------------------------------------------------------
  ! Type 'list_element' provides the entry to the actual information 
  ! and a reference to the next element in the list.

  TYPE list_element
    !
    SEQUENCE
    !
    TYPE (memory_type)          :: field
    TYPE(list_element), POINTER :: next_list_element
  END TYPE list_element
  !------------------------------------------------------------------------------
  ! Type 'list' provides the anchor for a full list. The memory allocated
  ! 'memory_used' as well an the number of elements in the list
  ! 'list_elements'is traced. A default value for the meta information is
  ! used to initialize new list entries.
  ! 
  ! As the starting point for the postprocessing additional information is
  ! provided: 'lpost' indicates that the buffer should be used in te
  ! standard postprocessing of gridpoint fields (called from scan1).  
  ! Otherwise the postprocessing routine must be called from elsewhere. If
  ! 'suffix' is provided a file different from the standard output file is
  ! used (with the given suffix). 'filetype' determines the type of the
  ! file and 'ztype' the file compression. 'unit' is used to store the unit 
  ! number (or file ID) as long as the file is open.

  TYPE t_stream
    !
    SEQUENCE
    !
    ! Memory buffer information
    !
    CHARACTER (len=16)           :: name               ! name of the buffer
    TYPE (list_element), POINTER :: first_list_element ! reference to 1.l.e.
    INTEGER(i8)                  :: memory_used        ! memory allocated
    INTEGER                      :: list_elements      ! no. elements alloc.
    TYPE (memory_info)           :: default_info       ! default meta data
    !
    ! Postprocessing information
    !
    LOGICAL                      :: lpost           ! standard postproc.stream
    LOGICAL                      :: lpout           ! switch postprocessing
    LOGICAL                      :: lrerun          ! standard restart  stream
    LOGICAL                      :: lcontnorest     ! run if restart not available
    LOGICAL                      :: linit           ! standard initial  stream
    CHARACTER(len=128)           :: filename        ! name of file
    CHARACTER(len=  8)           :: post_suf        ! suffix of output  file
    CHARACTER(len=  8)           :: rest_suf        ! suffix of restart file
    CHARACTER(len=  8)           :: init_suf        ! suffix of initial file
    !
    ! temporaries used by NetCDF, GRIB routines
    !
    LOGICAL                      :: first           ! first stream in file
    INTEGER                      :: post_idx        ! postproc. interval index
    INTEGER                      :: filetype        ! GRIB,NETCDF,NETCDF2,NETCDF4
    INTEGER                      :: ztype           ! NONE,SZIP,ZIP
    INTEGER                      :: fileID          ! file ID
    INTEGER                      :: vlistID         ! vlist ID
    INTEGER                      :: timestep        ! Output timestep
    !
    ! Mark last element of derived type to allow proper usage later 
    ! in derived MPI datatypes 
    INTEGER                      :: last          ! just dummy
  END TYPE t_stream

  !------------------------------------------------------------------------------
  ! default (empty) stream entry:

! Fix for PGI compiler that doesn't handle named constants (default_info) in
! type declaration with attribute PARAMETER
#if defined (__PGI)
  TYPE(t_stream), SAVE :: default_stream = &
#else
  TYPE(t_stream), PARAMETER :: default_stream = &
#endif
       t_stream(     &
       ' ',          & ! name of the buffer
       NULL(),       & ! reference to 1.l.e.
       0_i8,         & ! memory allocated
       0,            & ! no. elements alloc.
       default_info, &
       .false.,      & ! standard postproc.stream
       .false.,      & ! switch postprocessing
       .false.,      & ! standard restart  stream
       .false.,      & ! run if restart not available
       .false.,      & ! standard initial  stream
       ' ',          & ! name of file
       ' ',          & ! suffix of output  file
       ' ',          & ! suffix of restart file
       ' ',          & ! suffix of initial file
       .true.,       & ! first stream in file
       UNKNOWN,      & ! postproc. interval index
       UNKNOWN,      & ! 'GRIB','NetCDF'
       UNKNOWN,      & ! 'NONE','SZIP'
       UNKNOWN,      & ! file ID
       UNKNOWN,      & ! vlist ID
       UNKNOWN,      & ! Output timestep
       -1)

  !------------------------------------------------------------------------------
CONTAINS
  !------------------------------------------------------------------------------
  !
  ! initialize a variable of type memnuf with default values
  ! nullify anchor to linked list
  !
  SUBROUTINE construct_list (this_list)

    TYPE (t_stream) ,INTENT(out) :: this_list

    NULLIFY (this_list%first_list_element)

    this_list% name            = ' '
    this_list% memory_used     = 0
    this_list% list_elements   = 0
    this_list% default_info    = default_info

    this_list% lpost           = .TRUE.
    this_list% lpout           = .TRUE.
    this_list% lrerun          = .TRUE.
    this_list% linit           = .TRUE.
    this_list% filename        = ' '
    this_list% post_suf        = ' '
    this_list% rest_suf        = ' '
    this_list% init_suf        = ' '
    this_list% filetype        = out_filetype
    this_list% ztype           = out_ztype
    this_list% post_idx        = 1
    this_list% fileID          = -1
    this_list% vlistID         = -1
    this_list% timestep        = -1

  END SUBROUTINE construct_list
  !------------------------------------------------------------------------------
  !
  ! remove all elements of a linked list
  ! check if all elements are removed
  !
  SUBROUTINE destruct_list (this_list)

    TYPE (t_stream) ,INTENT(inout) :: this_list

    CALL destruct_list_element (this_list, this_list%first_list_element)

    NULLIFY (this_list%first_list_element)
    IF (this_list%memory_used /= 0) THEN
      CALL message('', &
           'List destructor didn''t work proper (memory counter) ...')
      CALL finish ('destruct_list',&
           'List destructor didnt work proper (memory counter)')
    ENDIF

    IF (this_list%list_elements /= 0) THEN
      CALL message('', &
           'List destructor didn''t work proper (element counter) ...')
      CALL finish ('destruct_list',&
           'List destructor didnt work proper (element counter)')
    ENDIF

    !ik_bugfix_R1.02a    CALL construct_list (this_list)

  END SUBROUTINE destruct_list
  !------------------------------------------------------------------------------
  !
  ! deallocate a list element and all its sucessors
  !
  SUBROUTINE destruct_list_element (this_list, this_list_element)

    TYPE (t_stream)     ,INTENT(inout) :: this_list
    TYPE (list_element) ,POINTER       :: this_list_element

    TYPE (list_element) ,POINTER       :: this, next

    next  => this_list_element
    NULLIFY (this_list_element)
    DO
      IF (.NOT. ASSOCIATED(next)) EXIT
      this => next
      next => this% next_list_element
      !
      ! 8 as constant has to be adjusted with information from mo_machine
      ! the variable to be used is mp_real8
      !
      IF (this% field% info% alloc) THEN
        this_list%memory_used = this_list%memory_used -8*SIZE(this%field%ptr)
        DEALLOCATE (this%field%ptr)
        this% field% info% alloc = .FALSE.
      ENDIF
      this_list%list_elements = this_list%list_elements-1
      DEALLOCATE (this)
    END DO

  END SUBROUTINE destruct_list_element
  !------------------------------------------------------------------------------
  subroutine remove_list_element (this_list, this_list_element)
    TYPE (t_stream)     ,INTENT(inout) :: this_list
    TYPE (list_element) ,POINTER :: this_list_element
    !
    ! remove one element from the list
    ! the element is identified by its pointer
    ! the pointer is associated with the next element
    !
    TYPE (list_element) ,POINTER :: p
    p => this_list_element
    this_list_element => this_list_element% next_list_element

    IF (p% field% info% alloc) THEN
      this_list%memory_used = this_list%memory_used -8*SIZE(p%field%ptr)
      DEALLOCATE (p%field%ptr)
      p% field% info% alloc = .FALSE.
    ENDIF
    this_list%list_elements = this_list%list_elements-1
    DEALLOCATE (p)

  end subroutine remove_list_element
  !------------------------------------------------------------------------------
  SUBROUTINE create_list_element (this_list, current_list_element)

    TYPE (t_stream)     ,INTENT(inout) :: this_list
    TYPE (list_element) ,POINTER       :: current_list_element
    INTEGER                            :: istat

    ALLOCATE (current_list_element, STAT=istat)
    IF (istat /= 0) THEN
      CALL message('', 'Cannot add element to linked list ...')
      CALL finish('create_list_element','Cannot add element to linked list')
    ENDIF

    this_list%list_elements = this_list%list_elements+1

    NULLIFY (current_list_element%next_list_element)
    NULLIFY (current_list_element%field%ptr)
    current_list_element%field%info = this_list% default_info

  END SUBROUTINE create_list_element
  !------------------------------------------------------------------------------
  SUBROUTINE add_list_element (this_list, new_list_element, code)
    TYPE (t_stream)     ,INTENT(inout) :: this_list
    TYPE (list_element) ,POINTER       :: new_list_element
    INTEGER   ,OPTIONAL ,INTENT(in)    :: code
    !--------------------------------------------------------
    ! Add a list element to the linked list
    ! if 'code' is present order with respect to the gribcode
    ! else search for end of the list
    !-------------------------------------------------------- 
    TYPE (list_element), POINTER :: current_list_element
    INTEGER                      :: ic
    ic = 1000; IF(PRESENT(code)) ic = code
    !-----------------------------------------
    ! insert as first element if list is empty
    !-----------------------------------------
    IF (.NOT. ASSOCIATED (this_list% first_list_element)) THEN
      CALL create_list_element (this_list, this_list% first_list_element)
      new_list_element => this_list% first_list_element
      RETURN
    ENDIF
    !----------------------------
    ! insert before first element
    !----------------------------
    IF (this_list% first_list_element% field% info% gribcode > ic) THEN
      CALL create_list_element (this_list, new_list_element)
      new_list_element% next_list_element => this_list% first_list_element
      this_list%       first_list_element => new_list_element
      RETURN
    ENDIF
    !-----------------------------------------
    ! loop over list elements to find position
    !-----------------------------------------
    current_list_element => this_list% first_list_element
    DO WHILE (ASSOCIATED(current_list_element% next_list_element)) 
      IF (current_list_element% next_list_element% field% info% gribcode > ic)&
           EXIT
      current_list_element => current_list_element% next_list_element
    ENDDO
    !---------------
    ! insert element
    !---------------
    CALL create_list_element (this_list, new_list_element)
    new_list_element%     next_list_element => current_list_element% &
         next_list_element
    current_list_element% next_list_element => new_list_element

  END SUBROUTINE add_list_element
  !------------------------------------------------------------------------------
  SUBROUTINE delete_list_element (this_list, delete_this_list_element)

    TYPE (t_stream)     ,INTENT(inout) :: this_list
    TYPE (list_element) ,POINTER       :: delete_this_list_element

    TYPE (list_element), POINTER :: current_list_element

    IF (ASSOCIATED(delete_this_list_element, &
         this_list%first_list_element)) THEN
      this_list%first_list_element &
           => delete_this_list_element%next_list_element
    ELSE
      current_list_element => this_list%first_list_element
      DO WHILE ((ASSOCIATED(current_list_element)) &
           .AND. (.NOT. ASSOCIATED(current_list_element%next_list_element, &
           delete_this_list_element)))
        current_list_element => current_list_element%next_list_element
      ENDDO
      IF (.NOT. ASSOCIATED(current_list_element)) THEN
        CALL message('', 'Cannot find element to be deleted ...')
        RETURN
      ENDIF
      current_list_element%next_list_element &
           => current_list_element%next_list_element%next_list_element
    ENDIF

    IF (delete_this_list_element% field% info% alloc) THEN

      this_list%memory_used = this_list% memory_used &
           -8*SIZE(delete_this_list_element% field% ptr)

      DEALLOCATE (delete_this_list_element% field% ptr)

      delete_this_list_element% field% info% alloc = .FALSE.

    ENDIF

    this_list%list_elements = this_list%list_elements-1

    DEALLOCATE (delete_this_list_element)

  END SUBROUTINE delete_list_element
  !------------------------------------------------------------------------------
  ! Should be overloaded to be able to search for the different information 
  ! In the proposed structure for the linked list, in the example only
  ! A character string is used so it is straight forward only one find

  FUNCTION find_list_element (this_list, name) RESULT (this_list_element)

    TYPE (t_stream)   ,INTENT(in) :: this_list
    CHARACTER (len=*) ,INTENT(in) :: name

    TYPE (list_element), POINTER :: this_list_element

    this_list_element => this_list%first_list_element
    DO WHILE (ASSOCIATED(this_list_element))
      IF (name == this_list_element%field%info%name) THEN
        RETURN
      ENDIF
      this_list_element => this_list_element%next_list_element
    ENDDO

    NULLIFY (this_list_element)

  END FUNCTION find_list_element
  !------------------------------------------------------------------------------
  SUBROUTINE print_linked_list (this_list)

    TYPE (t_stream) ,INTENT(in) :: this_list

    TYPE (list_element) ,POINTER :: this_list_element

    this_list_element => this_list%first_list_element
    DO WHILE (ASSOCIATED(this_list_element))

      ! print ....
      IF (this_list_element%field%info%name /= '') THEN

        WRITE (message_text,'(a,a)')       &
             'Table entry name      : ', &
             TRIM(this_list_element%field%info%name)
        CALL message('', message_text)
#ifdef DEBUG
        WRITE (message_text,'(a,i10)') 'Address of data field : ', &
             111111111
!            LOC(this_list_element%field%ptr)
        CALL message('', message_text)
#endif
        IF (ASSOCIATED(this_list_element%field%ptr)) THEN
          CALL message ('','Pointer status        : in use.')

          IF (SIZE(this_list_element%field%ptr,4) == 1 & 
               .AND. SIZE(this_list_element%field%ptr,3) == 1 &
               .AND. SIZE(this_list_element%field%ptr,2) == 1) THEN
            WRITE (message_text,'(a,1(i4,a))') &
                 'Local field dimensions      : (',  &
                 SIZE(this_list_element%field%ptr,1), ')' 
            CALL message('', message_text)
          ELSE IF (SIZE(this_list_element%field%ptr,4) == 1 &
               .AND. SIZE(this_list_element%field%ptr,3) == 1) THEN
            WRITE (message_text,'(a,2(i4,a))') &
                 'Local field dimensions      : (',  &
                 SIZE(this_list_element%field%ptr,1), ',', &
                 SIZE(this_list_element%field%ptr,2), ')' 
            CALL message('', message_text)
          ELSE IF (SIZE(this_list_element%field%ptr,4) == 1) THEN
            WRITE (message_text,'(a,3(i4,a))') &
                 'Local field dimensions      : (',  &
                 SIZE(this_list_element%field%ptr,1), ',', &
                 SIZE(this_list_element%field%ptr,2), ',', &
                 SIZE(this_list_element%field%ptr,3), ')'       
            CALL message('', message_text)
          ELSE
            WRITE (message_text,'(a,4(i4,a))') &
                 'Local field dimensions      : (',  &
                 SIZE(this_list_element%field%ptr,1), ',', &
                 SIZE(this_list_element%field%ptr,2), ',', &
                 SIZE(this_list_element%field%ptr,3), ',', &
                 SIZE(this_list_element%field%ptr,4), ')'       
            CALL message('', message_text)
          ENDIF
        ELSE
          CALL message('', 'Pointer status       : not in use.')
        ENDIF

        IF (this_list_element%field%info%gdim(4) /= 0 & 
             .AND. this_list_element%field%info%gdim(3) /= 0 &
             .AND. this_list_element%field%info%gdim(2) /= 0) THEN
          WRITE (message_text,'(a,1(i4,a))') &
               'Global field dimensions      : (',  &
               this_list_element%field%info%gdim(1), ')' 
          CALL message('', message_text)
        ELSE IF (this_list_element%field%info%gdim(4) /= 0 &
             .AND. this_list_element%field%info%gdim(3) /= 0) THEN
          WRITE (message_text,'(a,2(i4,a))') &
               'Global field dimensions      : (',  &
               this_list_element%field%info%gdim(1), ',', &
               this_list_element%field%info%gdim(2), ')' 
          CALL message('', message_text)
        ELSE IF (this_list_element%field%info%gdim(4) /= 0) THEN
          WRITE (message_text,'(a,3(i4,a))') &
               'Global field dimensions      : (',  &
               this_list_element%field%info%gdim(1), ',', &
               this_list_element%field%info%gdim(2), ',', &
               this_list_element%field%info%gdim(3), ')'       
          CALL message('', message_text)
        ELSE
          WRITE (message_text,'(a,4(i4,a))') &
               'Global field dimensions      : (',  &
               this_list_element%field%info%gdim(1), ',', &
               this_list_element%field%info%gdim(2), ',', &
               this_list_element%field%info%gdim(3), ',', &
               this_list_element%field%info%gdim(4), ')'       
          CALL message('', message_text)
        ENDIF

        WRITE (message_text,'(a,i3,/,a,i3)') &
             'Assigned GRIB table   : ', &
             this_list_element%field%info%gribtable, &
             '         GRIB code    : ', &
             this_list_element%field%info%gribcode
        CALL message('', message_text)
        WRITE (message_text,'(a,i6,/,a,a,/,a,a)') &
             'IO id                 : ', &
             this_list_element%field%info%IO_var_id, &
             '   name               : ', &
             TRIM(this_list_element%field%info%IO_name), &
             '   unit               : ', &
             TRIM(this_list_element%field%info%IO_unit)
        CALL message('', message_text)
        IF (this_list_element%field%info%laccu) THEN
          CALL message('', 'Accumulation          : on.')
        ELSE
          CALL message('', 'Accumulation          : off.')
        ENDIF
        IF (this_list_element%field%info%lmiss) THEN
          WRITE (message_text,'(a,e20.12)')      &
               'Missing value         : ', &
               this_list_element%field%info%missval
          CALL message('', message_text)
        ELSE
          CALL message('', 'Missing values        : off.')
        ENDIF
        IF (this_list_element%field%info%lrerun) THEN
          CALL message('', 'Restart table         : added.')
        ELSE
          CALL message('', 'Restart table         : unused.')
        ENDIF
        CALL message('', '')
      ENDIF

      ! select next element in linked list 

      this_list_element => this_list_element%next_list_element
    ENDDO

  END SUBROUTINE print_linked_list
  !------------------------------------------------------------------------------
  SUBROUTINE print_sinfo_list (this_list)

    TYPE (t_stream)       ,INTENT(in) :: this_list

    TYPE (list_element) ,POINTER :: this_list_element
    CHARACTER (len=80)           :: cout

    cout = ''

    CALL message('','   Name   Local dimension   Tab Code Outint Accu  Restart')
    this_list_element => this_list%first_list_element
    DO WHILE (ASSOCIATED(this_list_element))

      ! print ....
      IF (this_list_element%field%info%name /= '') THEN

        WRITE(cout(1:10),'(a10)') TRIM(this_list_element%field%info%name)

        IF (ASSOCIATED(this_list_element%field%ptr)) THEN
          IF (SIZE(this_list_element%field%ptr,4) == 1 & 
               .AND. SIZE(this_list_element%field%ptr,3) == 1 &
               .AND. SIZE(this_list_element%field%ptr,2) == 1) THEN
            WRITE (cout(12:30),'(a,1(i3,a))') &
                 ' (',  &
                 SIZE(this_list_element%field%ptr,1), ')' 
          ELSE IF (SIZE(this_list_element%field%ptr,4) == 1 &
               .AND. SIZE(this_list_element%field%ptr,3) == 1) THEN
            WRITE (cout(12:30),'(a,2(i3,a))') &
                 ' (',  &
                 SIZE(this_list_element%field%ptr,1), ',', &
                 SIZE(this_list_element%field%ptr,2), ')' 
          ELSE IF (SIZE(this_list_element%field%ptr,4) == 1) THEN
            WRITE (cout(12:30),'(a,3(i3,a))') &
                 ' (',  &
                 SIZE(this_list_element%field%ptr,1), ',', &
                 SIZE(this_list_element%field%ptr,2), ',', &
                 SIZE(this_list_element%field%ptr,3), ')'       
          ELSE
            WRITE (cout(12:30),'(a,4(i3,a))') &
                 ' (',  &
                 SIZE(this_list_element%field%ptr,1), ',', &
                 SIZE(this_list_element%field%ptr,2), ',', &
                 SIZE(this_list_element%field%ptr,3), ',', &
                 SIZE(this_list_element%field%ptr,4), ')'       
          ENDIF
        ELSE
          WRITE (cout(12:30),'(a)')      &
               '    not in use '
        ENDIF
        WRITE (cout(31:33),'(i3)') this_list_element%field%info%gribtable
        WRITE (cout(34:38),'(i5)') this_list_element%field%info%gribcode

        IF (this_list_element%field%info%laccu) THEN
          WRITE (cout(47:52),'(a)') ' on  '
        ELSE
          WRITE (cout(47:52),'(a)') ' off '
        ENDIF
        IF (this_list_element%field%info%lrerun) THEN
          WRITE (cout(53:60),'(a)') ' added '
        ELSE
          WRITE (cout(53:60),'(a)') ' unused'
        ENDIF
      ENDIF

      CALL message ('', cout)

      ! select next element in linked list 

      this_list_element => this_list_element%next_list_element
    ENDDO

  END SUBROUTINE print_sinfo_list
  !------------------------------------------------------------------------------
  !
  ! printout of memory_info content
  !
  SUBROUTINE print_memory_info (info)
    TYPE (memory_info) ,INTENT(in) :: info

!OSBUTLS    CALL message('',separator)
    CALL message('','')
    CALL message('','  Stream Element Info :')
    CALL message('','')

    WRITE(message_text,'(a,a )')  '  name               : ',info% name
    CALL message('', message_text)
    WRITE(message_text,'(a,a )')  '  units              : ',info% units
    CALL message('', message_text)
    WRITE(message_text,'(a,a )')  '  long name          : ',TRIM(info% longname)

    CALL message('','')

    WRITE(message_text,'(a,5i4)') '  local  dimensions  : ',info% DIM(:)
    CALL message('', message_text)
    WRITE(message_text,'(a,5i4)') '  global dimensions  : ',info%gdim(:)
    CALL message('', message_text)
    WRITE(message_text,'(a, i4)') '  rank               : ',info% ndim
    CALL message('', message_text)
    WRITE(message_text,'(a, i4)') '  number of levels   : ',info% klev
    CALL message('', message_text)
    WRITE(message_text,'(a, l4)') '  allocated flag     : ',info% alloc
    CALL message('', message_text)
    WRITE(message_text,'(a, i4)') '  grid type          : ',info% repr
    CALL message('', message_text)

    CALL message('','')

    WRITE(message_text,'(a, l4)') '  postprocessing flag: ',info% lpost
    CALL message('', message_text)
    WRITE(message_text,'(a, l4)') '  accumulate flag    : ',info% laccu
    CALL message('', message_text)
    WRITE(message_text,'(a,f9.2)')'  reset value        : ',info% reset
    CALL message('', message_text)
    WRITE(message_text,'(a, l4)') '  restart flag       : ',info% lrerun
    CALL message('', message_text)
    WRITE(message_text,'(a, l4)') '  restart feedback   : ',info% restart_read
    CALL message('', message_text)
    WRITE(message_text,'(a,a )')  '  prefix for restart : ',info% restart_pref
    CALL message('', message_text)

    CALL message('','')

    WRITE(message_text,'(a, i4)') '  GRIB table         : ',info% gribtable
    CALL message('', message_text)
    WRITE(message_text,'(a, i4)') '  GRIB code          : ',info% gribcode
    CALL message('', message_text)
    WRITE(message_text,'(a, i4)') '  GRIB bits          : ',info% gribbits
    CALL message('', message_text)

    CALL message('','')

    WRITE(message_text,'(a, i4)') '  levelindx          : ',info% levelindx
    CALL message('', message_text)

    CALL message('','')

    WRITE(message_text,'(a,4i4)') '  IO_var_indx (int.) : ',info% IO_var_indx
    CALL message('', message_text)
    WRITE(message_text,'(a, i4)') '  IO_var_id (int.use): ',info% IO_var_id
    CALL message('', message_text)
    WRITE(message_text,'(a,a )')  '  NETCDF file        : ',TRIM(info% IO_name)
    CALL message('', message_text)
    WRITE(message_text,'(a,a )')  '  NETCDF unit        : ',TRIM(info% IO_unit)
    CALL message('', message_text)
    WRITE(message_text,'(a, i4)') '  gridID             : ',info% gridID
    CALL message('', message_text)
    WRITE(message_text,'(a, i4)') '  zaxisID            : ',info% zaxisID
    CALL message('', message_text)

    CALL message('','')
    
  END SUBROUTINE print_memory_info
  !------------------------------------------------------------------------------
  !
  ! print information on an output stream
  ! print one line for each element of the list in tabular form
  !

  SUBROUTINE print_stream (stream)

    ! Print output stream

    TYPE (t_stream) ,INTENT(in) :: stream ! output stream to print

    TYPE (list_element) ,POINTER   :: next
    TYPE (memory_info)  ,POINTER   :: info
    CHARACTER(len=10)              :: grid
    CHARACTER(len=10)              :: ltype

!OSBUTLS    CALL message('',separator)
    CALL message('','')
    CALL message('','  Output stream :')
    CALL message('','')

    WRITE(message_text,'(a,a  )') '  name               : ',stream% name
    CALL message('', message_text)
    WRITE(message_text,'(a,i8 )') '  memory used        : ',stream% memory_used
    CALL message('', message_text)
    WRITE(message_text,'(a,i8 )') '  list elements      : ',stream% list_elements
    CALL message('', message_text)

    CALL message('','')

    WRITE(message_text,'(a,l1 )') '  lpost flag         : ',stream% lpost
    CALL message('', message_text)
    WRITE(message_text,'(a,l1 )') '  lpout flag         : ',stream% lpout
    CALL message('', message_text)
    WRITE(message_text,'(a,l1 )') '  linit flag         : ',stream% linit
    CALL message('', message_text)
    WRITE(message_text,'(a,l1 )') '  lrerun flag        : ',stream% lrerun
    CALL message('', message_text)
    WRITE(message_text,'(a,a  )') '  output  file suffix: ',stream% post_suf
    CALL message('', message_text)
    WRITE(message_text,'(a,a  )') '  restart file suffix: ',stream% rest_suf
    CALL message('', message_text)
    WRITE(message_text,'(a,a  )') '  initial file suffix: ',stream% init_suf
    CALL message('', message_text)

    SELECT CASE (stream% filetype)
    CASE (GRIB)
      CALL message('','  filetype           : GRIB1')
    CASE (NETCDF)
      CALL message('','  filetype           : NETCDF')
    CASE (NETCDF2)
      CALL message('','  filetype           : NETCDF2')
    CASE (NETCDF4)
      CALL message('','  filetype           : NETCDF4')
    CASE default
      CALL message('','  filetype           : UNKNOWN')
    END SELECT

    SELECT CASE (stream% ztype)
    CASE (NONE)
      CALL message('','  output compression : NONE')
    CASE (SZIP)
      CALL message('','  output compression : SZIP')
    CASE (ZIP)
      CALL message('','  output compression : ZIP')
    CASE default
      CALL message('','  output compression : UNKNOWN')
    END SELECT
    WRITE(message_text,'(a,i8 )') '  interval index     : ',stream% post_idx
    CALL message('', message_text)
    WRITE(message_text,'(a,i20)') '  fileID             : ',stream% fileID
    CALL message('', message_text)

    next => stream% first_list_element

    CALL message('','')
    CALL message('','  name            units rank  ke alloc.  grid prt acc mis rst tbl cde bit lev_type')
    CALL message('','')

    DO
      IF (.NOT.ASSOCIATED(next)) EXIT
      info => next% field% info

      SELECT CASE (info% repr)
      CASE     (GAUSSIAN)
        grid = 'GAUSSIAN'
      CASE     (FOURIER)
        grid = 'FOURIER'
      CASE     (SPECTRAL)
        grid = 'SPECTRAL'
      CASE     (HEXAGONAL)
        grid = 'HEXAGONAL'
      CASE     (LAND)
        grid = 'LAND'
      CASE default
        grid = 'UNKNOWN'
      END SELECT

      ltype                                    = 'UNKNOWN'
      IF (info% levelindx == BELOWSUR)   ltype = 'BELOWSUR'
      IF (info% levelindx == SURFACE)    ltype = 'SURFACE'
      IF (info% levelindx == ABOVESUR2)  ltype = 'ABOVESUR2'
      IF (info% levelindx == ABOVESUR10) ltype = 'ABOVESUR10'
      IF (info% levelindx == HYBRID)     ltype = 'HYBRID'
      IF (info% levelindx == HYBRID_H)   ltype = 'HYBRID_H'
      IF (info% levelindx == TILES)      ltype = 'TILES'
      IF (info% levelindx == SOILLEV)    ltype = 'SOILLEV'
      IF (info% levelindx == ROOTZONES)  ltype = 'ROOTZONES'
      IF (info% levelindx == CANOPY)     ltype = 'CANOPY'

      WRITE(message_text,'(2x,a16,a8,i2,i4,l2,1x,a10,l3,l4,l4,l4,i5,i4,i4,1x,a8)') &
           info% name       ,&
           info% units      ,&
           info% ndim       ,&
           info% klev       ,&
           info% alloc      ,&
           grid       ,&
           info% lpost      ,&
           info% laccu      ,&
           info% lmiss      ,&
           info% lrerun     ,&
           info% gribtable  ,&
           info% gribcode   ,&
           info% gribbits   ,&
           ltype

      CALL message ('', message_text)
      next => next% next_list_element
    END DO

    CALL message('','')

  END SUBROUTINE print_stream
  !------------------------------------------------------------------------------
END MODULE mo_linked_list
