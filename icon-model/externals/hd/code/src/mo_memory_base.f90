! mo_memory_base.f90 - Memory manager routines for output or memory streams
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!_________________________________________

MODULE mo_memory_base

  ! Routines to insert new buffers into an output or memory stream 
  ! (new_stream_element) or get the reference to an existing entry 
  ! (get_stream_element)
  !
  ! Authors:
  !
  ! Luis Kornblueh, MPI,             original code
  ! Andreas Rhodin, MPI, April 2001, extended
  ! I. Kirchner,    MPI, August 2002 lpout flag
  !
  !  This routine originates (year 2014) from MPI-ESM, the Earth System Model of the 
  !  Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
  !  Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
  !  version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
  !  doi: 10.1029/2018MS001400.
  !----------------------------------------------------------------------------
  !
  ! Modules used
  !

  USE mo_kind,         ONLY: dp 
  USE mo_exception,    ONLY: finish, message, message_text
  USE mo_linked_list,  ONLY: memory_info, memory_type, t_stream, list_element,&
                             add_list_element, find_list_element,             &
                             remove_list_element,                             &
                             construct_list, destruct_list,                   &
                             print_linked_list, print_sinfo_list,             &
                             print_memory_info, default_info,                 &
                             GRIB, NETCDF, NETCDF2, NETCDF4,                  &
                             BELOWSUR, SURFACE, ABOVESUR2, ABOVESUR10,        &
                             HYBRID, HYBRID_H, UNKNOWN,                       &
                             GAUSSIAN, FOURIER, SPECTRAL, HEXAGONAL,          &
                             GRIDPOINT, LAND, TILES, SOILLEV, CANOPY
  USE mo_netCDF,       ONLY: IO_dim_ids, IO_get_varindx, &
                             add_dim, add_unknown_dim
  USE mo_decomposition,ONLY: ldc => local_decomposition
  USE mo_time_event,   ONLY: io_time_event
  USE mo_time_control, ONLY: get_new_ev_putdata
  USE mo_control,      ONLY: nsp
  USE mo_jsbach_comm_to_echam5mods,  ONLY: nland, domain_nland

  IMPLICIT NONE
  !----------------------------------------------------------------------------
  !
  ! Public entities
  !
  
  PRIVATE
  
  PUBLIC :: new_stream             ! get a pointer to a new output stream
  PUBLIC :: get_stream             ! get a pointer to an old output stream
  PUBLIC :: set_stream             ! set parameters of the output stream
  PUBLIC :: delete_stream          ! delete an output stream
  PUBLIC :: delete_streams         ! delete all output streams
  PUBLIC :: t_stream               ! output stream data type

  PUBLIC :: ostreams               ! list of output streams
  PUBLIC :: nstreams               ! number of output streams defined so far

  PUBLIC :: GRIB, NETCDF, NETCDF2, NETCDF4  ! valid file types
  PUBLIC :: GAUSSIAN, FOURIER, SPECTRAL, HEXAGONAL, GRIDPOINT, LAND ! flag 'repr'
  PUBLIC :: BELOWSUR, SURFACE, ABOVESUR2, ABOVESUR10, &       ! 'level_type'
            HYBRID, HYBRID_H, TILES, SOILLEV, CANOPY
  PUBLIC :: UNKNOWN
  PUBLIC :: AUTO                   ! determe code automatically

  PUBLIC :: add_stream_element     ! create/allocate a new stream list entry
  PUBLIC :: get_stream_element     ! obtain reference to existing list entry
  PUBLIC :: remove_stream_element  ! delete a list element
  PUBLIC :: add_stream_reference   ! add a reference to another element
  PUBLIC :: default_stream_setting ! change default output stream settings

  PUBLIC :: get_stream_element_info ! obtain meta data of list entry
  PUBLIC :: set_stream_element_info ! set meta data
  PUBLIC :: print_memory_table      ! print list/memory information
  PUBLIC :: print_memory_use        ! print used memory
  PUBLIC :: print_sinfo             ! print short information
  
  PUBLIC :: memory_info            ! meta data structure
  PUBLIC :: suffix                 ! preliminary: convert unit number to suffix
  PUBLIC :: maxstr                 ! max number of output streams

  PUBLIC :: nofiles                ! number of closed files
  PUBLIC :: cfilenames             ! cmd with names of closed files

!------------------------------------------------------------------------------
  !
  ! Interfaces
  !

  INTERFACE add_stream_element
     MODULE PROCEDURE add_stream_element_4d ! create a new list entry
     MODULE PROCEDURE add_stream_element_3d 
     MODULE PROCEDURE add_stream_element_2d 
     MODULE PROCEDURE add_stream_element_1d 
  END INTERFACE
  
  INTERFACE get_stream_element
     MODULE PROCEDURE get_stream_element_4d ! obtain reference to a list entry
     MODULE PROCEDURE get_stream_element_3d
     MODULE PROCEDURE get_stream_element_2d
     MODULE PROCEDURE get_stream_element_1d
  END INTERFACE

  INTERFACE assign_if_present
     MODULE PROCEDURE assign_if_present_character
     MODULE PROCEDURE assign_if_present_logical
     MODULE PROCEDURE assign_if_present_integer
     MODULE PROCEDURE assign_if_present_integers
     MODULE PROCEDURE assign_if_present_real
  END INTERFACE
!------------------------------------------------------------------------------
  !
  ! Module variables
  !

  INTEGER      ,PARAMETER    :: AUTO     = -1    ! determine GRIB code autom.
  INTEGER      ,PARAMETER    :: maxstr   = 50    ! max number of output streams
  LOGICAL                    :: lstream_init = .TRUE.
                               ! true during the initialization loop of streams
  INTEGER                     :: nstreams =  0    ! streams allocated so far
  TYPE(t_stream),TARGET, SAVE :: ostreams(maxstr) ! memory buffer array

  ! subjob interface
  INTEGER, SAVE           :: nofiles = 0        ! number of closed files    
  CHARACTER(len=256),SAVE :: cfilenames(maxstr) ! cmd with names of closed files

!==============================================================================
CONTAINS
!==============================================================================
  !
  ! Get a reference to a memory buffer / output stream
  !
  SUBROUTINE get_stream (stream, name)
  TYPE (t_stream)  ,POINTER              :: stream ! pointer
  CHARACTER(len=*) ,INTENT(in)           :: name   ! name of output stream

    INTEGER :: i

    NULLIFY (stream)
    DO i=1, nstreams
      IF (ostreams(i)% name == name) THEN
        stream => ostreams(i)
        EXIT
      ENDIF
    END DO

  END SUBROUTINE get_stream
!==============================================================================
  !
  ! Create a new memory buffer / output stream
  ! Get a pointer to the new stream
  !

  SUBROUTINE new_stream (stream, name, filetype, ztype,     &
                         post_suf, rest_suf, init_suf,      &
                         lpost, lpout, lrerun, lcontnorest, &
                         linit, interval)

#ifdef __SUNPRO_F95
  TYPE (t_stream)    ,POINTER              :: stream    ! pointer
#else
  TYPE (t_stream)    ,POINTER, INTENT(out) :: stream    ! pointer
#endif
  CHARACTER(len=*)   ,INTENT(in)           :: name      ! name of output stream
  INTEGER            ,INTENT(in) ,OPTIONAL :: filetype  ! GRIB, NETCDF, NETCDF2, NETCDF4
  INTEGER            ,INTENT(in) ,OPTIONAL :: ztype     ! NONE, SZIP, ZIP
  CHARACTER(len=*)   ,INTENT(in) ,OPTIONAL :: post_suf  ! suffix of outputfile
  CHARACTER(len=*)   ,INTENT(in) ,OPTIONAL :: rest_suf  ! suffix of restartfile
  CHARACTER(len=*)   ,INTENT(in) ,OPTIONAL :: init_suf  ! suffix of initialfile
  LOGICAL            ,INTENT(in) ,OPTIONAL :: lpost     ! open outputfile
  LOGICAL            ,INTENT(in) ,OPTIONAL :: lpout     ! write to outputfile
  LOGICAL            ,INTENT(in) ,OPTIONAL :: lrerun    ! write to restartfile
  LOGICAL            ,INTENT(in) ,OPTIONAL :: lcontnorest
  LOGICAL            ,INTENT(in) ,OPTIONAL :: linit     ! write to initialfile
  TYPE(io_time_event),INTENT(in) ,OPTIONAL :: interval  ! output interval
    !
    ! Local variables
    !
    INTEGER             :: i
    TYPE(io_time_event) :: int
    !
    CALL message('new_stream','Adding new stream '//name)
    !
    ! name must be unique
    !
    NULLIFY (stream)

    IF (lstream_init) THEN
      !
      ! look for old name in list
      !
      IF (ANY(ostreams(1:nstreams)% name == name)) CALL finish('new_list',&
                      'output stream '//TRIM(name)//' already used.')
      !
      ! find free place
      !
      DO i=1, nstreams    
        IF (ostreams(i)% name == '') THEN
          stream => ostreams(i)
          EXIT
        ENDIF
      END DO
      IF(.NOT.ASSOCIATED(stream)) THEN
        nstreams = nstreams + 1
        IF (nstreams > maxstr) CALL finish('new_list',&
                      'output streams exhausted, increase "maxstr"')
        stream => ostreams(nstreams)
      ENDIF
      CALL construct_list (stream)
      !
      ! set default list characteristics
      !
      stream% name     = name
      stream% post_suf = '_'//name
      stream% rest_suf = stream% post_suf
      stream% init_suf = stream% post_suf
      !
      ! set non-default list characteristics
      !
      IF (PRESENT(interval   )) int = interval
      IF (PRESENT(interval   )) stream% post_idx    = get_new_ev_putdata (int)
      IF (PRESENT(lpost      )) stream% lpost       = lpost
      IF (PRESENT(lpout      )) stream% lpout       = lpout
      IF (PRESENT(lrerun     )) stream% lrerun      = lrerun
      IF (PRESENT(lcontnorest)) stream% lcontnorest = lcontnorest
      IF (PRESENT(linit      )) stream% linit       = linit
      IF (PRESENT(post_suf   )) stream% post_suf    = post_suf 
      IF (PRESENT(rest_suf   )) stream% rest_suf    = rest_suf 
      IF (PRESENT(init_suf   )) stream% init_suf    = init_suf 
      IF (PRESENT(filetype   )) stream% filetype    = filetype
      IF (PRESENT(ztype      )) stream% ztype       = ztype

    ELSE

      ! definition of stream structure is available
      CALL get_stream(stream, name)
 
    END IF
  END SUBROUTINE new_stream
!==============================================================================
  !
  ! Change parameters of an already existent output stream
  !

  SUBROUTINE set_stream (stream, filetype, ztype,      &
                         post_suf, rest_suf, init_suf, &
                         lpost, lpout, lrerun, lcontnorest,linit, interval)

  TYPE (t_stream)    ,INTENT(inout)       :: stream   ! output stream to change
  INTEGER            ,INTENT(in),OPTIONAL :: filetype ! GRIB, NETCDF, NETCDF2, NETCDF4
  INTEGER            ,INTENT(in),OPTIONAL :: ztype    ! NONE, SZIP, ZIP
  CHARACTER(len=*)   ,INTENT(in),OPTIONAL :: post_suf ! suffix of output  file
  CHARACTER(len=*)   ,INTENT(in),OPTIONAL :: rest_suf ! suffix of restart file
  CHARACTER(len=*)   ,INTENT(in),OPTIONAL :: init_suf ! suffix of initial file
  LOGICAL            ,INTENT(in),OPTIONAL :: lpost    ! in standard output file
  LOGICAL            ,INTENT(in),OPTIONAL :: lpout    ! in standard output file
  LOGICAL            ,INTENT(in),OPTIONAL :: lrerun   ! in standard restartfile
  LOGICAL            ,INTENT(in) ,OPTIONAL :: lcontnorest
  LOGICAL            ,INTENT(in),OPTIONAL :: linit    ! in standard initialfile
  TYPE(io_time_event),INTENT(in),OPTIONAL :: interval ! output interval

    TYPE(io_time_event) :: int

    IF (PRESENT(interval)) int = interval
    IF (PRESENT(interval)) stream% post_idx = get_new_ev_putdata (int)
    IF (PRESENT(lpost   )) stream% lpost    = lpost
    IF (PRESENT(lpout   )) stream% lpout    = lpout
    IF (PRESENT(lrerun  )) stream% lrerun   = lrerun
    IF (PRESENT(lcontnorest)) stream% lcontnorest = lcontnorest
    IF (PRESENT(linit   )) stream% linit    = linit
    IF (PRESENT(post_suf)) stream% post_suf = post_suf 
    IF (PRESENT(rest_suf)) stream% rest_suf = rest_suf 
    IF (PRESENT(init_suf)) stream% init_suf = init_suf 
    IF (PRESENT(filetype)) stream% filetype = filetype  
    IF (PRESENT(ztype   )) stream% ztype    = ztype

  END SUBROUTINE set_stream
!==============================================================================
  !
  ! Delete an output stream, nullify the associated pointer
  !

  SUBROUTINE delete_stream (stream)

    TYPE (t_stream) ,POINTER :: stream

    IF (ASSOCIATED (stream)) THEN
      CALL destruct_list (stream)
      NULLIFY (stream)
    ENDIF

  END SUBROUTINE delete_stream
!==============================================================================
  !
  ! Delete all output streams 
  !

  SUBROUTINE delete_streams
  TYPE (t_stream) ,POINTER :: buffer

    INTEGER :: i

    DO i=1,nstreams
      buffer => ostreams(i)
      CALL delete_stream (buffer)
    END DO

    lstream_init = .FALSE.

  END SUBROUTINE delete_streams
!==============================================================================
  !
  ! Get a copy of the metadata concerning a stream element
  !

  SUBROUTINE get_stream_element_info (stream, name, info)
    TYPE (t_stream)    ,INTENT(in)  :: stream ! list
    CHARACTER (*)      ,INTENT(in)  :: name   ! name of variable
    TYPE (memory_info) ,INTENT(out) :: info   ! variable meta data

    TYPE (list_element), POINTER :: p

    p => find_list_element (stream, name)
  
    IF (ASSOCIATED (p)) THEN
      info = p% field% info
    ELSE
      info = default_info
    ENDIF

  END SUBROUTINE get_stream_element_info
!==============================================================================
  !
  ! Set default meta data of output stream
  !

  SUBROUTINE default_stream_setting (stream ,units                            &
                           ,ldims ,gdims ,repr                                &
                           ,lpost ,laccu, reset, lrerun ,contnorest           &
                           ,table ,code ,bits, leveltype                      &
                           ,dimnames, no_default)

  TYPE (t_stream)  ,INTENT(inout)        :: stream        ! output stream

  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: units         ! units

  INTEGER          ,INTENT(in) ,OPTIONAL :: ldims (:)     ! local dimensions
  INTEGER          ,INTENT(in) ,OPTIONAL :: gdims (:)     ! global dimensions
  INTEGER          ,INTENT(in) ,OPTIONAL :: repr          ! representation

  LOGICAL          ,INTENT(in) ,OPTIONAL :: lpost         ! open output stream
  LOGICAL          ,INTENT(in) ,OPTIONAL :: laccu         ! accumulation flag
  REAL(dp)         ,INTENT(in) ,OPTIONAL :: reset         ! reset value
  LOGICAL          ,INTENT(in) ,OPTIONAL :: lrerun        ! rerun file flag
  LOGICAL          ,INTENT(in) ,OPTIONAL :: contnorest    ! continue on restart

  INTEGER          ,INTENT(in) ,OPTIONAL :: table         ! gribtable number
  INTEGER          ,INTENT(in) ,OPTIONAL :: code          ! gribcode number
  INTEGER          ,INTENT(in) ,OPTIONAL :: bits          ! bits used for GRIB
  INTEGER          ,INTENT(in) ,OPTIONAL :: leveltype     ! grib level type

  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: dimnames(:)   ! dimension names

  LOGICAL          ,INTENT(in) ,OPTIONAL :: no_default    ! use no defaults

    CALL set_memory_info (stream% default_info                                &
     ,units=units ,ldims=ldims ,gdims=gdims ,repr=repr                        &
     ,ndim=99                                                                 &
     ,lpost=lpost ,laccu=laccu ,reset=reset ,lrerun=lrerun                    &
     ,contnorest=contnorest                                                   &
     ,table=table ,code=code ,bits=bits ,leveltype=leveltype                  &
     ,dimnames=dimnames, no_default=no_default, autodet=.FALSE.)

  END SUBROUTINE default_stream_setting
!------------------------------------------------------------------------------
  !
  ! Change metadata of a stream element
  !

  SUBROUTINE set_stream_element_info (stream, name  ,longname ,units         &
                          ,ldims ,gdims ,ndim ,klev ,ktrac ,alloc ,repr      &
                          ,lpost ,laccu, lmiss, missval, reset, lrerun       &
                          ,contnorest ,table ,code ,bits, leveltype          &
                          ,dimnames, no_default)

  TYPE (t_stream)  ,INTENT(in)           :: stream        ! i/o stream
  CHARACTER(len=*) ,INTENT(in)           :: name          ! variable name
  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: units         ! units
  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: longname      ! long name

  INTEGER          ,INTENT(in) ,OPTIONAL :: ldims (:)     ! local dimensions
  INTEGER          ,INTENT(in) ,OPTIONAL :: gdims (:)     ! global dimensions
  INTEGER          ,INTENT(in) ,OPTIONAL :: ndim          ! rank
  INTEGER          ,INTENT(in) ,OPTIONAL :: klev          ! number of levels
  INTEGER          ,INTENT(in) ,OPTIONAL :: ktrac         ! number of 'tracers'
  LOGICAL          ,INTENT(in) ,OPTIONAL :: alloc         ! pointer allocated
  INTEGER          ,INTENT(in) ,OPTIONAL :: repr          ! representation

  LOGICAL          ,INTENT(in) ,OPTIONAL :: lpost         ! into output stream
  LOGICAL          ,INTENT(in) ,OPTIONAL :: laccu         ! accumulation flag
  LOGICAL          ,INTENT(in) ,OPTIONAL :: lmiss         ! missing value flag
  REAL(dp)         ,INTENT(in) ,OPTIONAL :: missval       ! missing value
  REAL(dp)         ,INTENT(in) ,OPTIONAL :: reset         ! reset value
  LOGICAL          ,INTENT(in) ,OPTIONAL :: lrerun        ! rerun file flag
  LOGICAL          ,INTENT(in) ,OPTIONAL :: contnorest    ! continue on restart

  INTEGER          ,INTENT(in) ,OPTIONAL :: table         ! gribtable number
  INTEGER          ,INTENT(in) ,OPTIONAL :: code          ! gribcode number
  INTEGER          ,INTENT(in) ,OPTIONAL :: bits          ! bits used for GRIB
  INTEGER          ,INTENT(in) ,OPTIONAL :: leveltype     ! grib level type

  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: dimnames(:)   ! dimension names
  LOGICAL          ,INTENT(in) ,OPTIONAL :: no_default    ! use no defaults

    TYPE (list_element), POINTER :: element

    element => find_list_element (stream, name)

    IF (ASSOCIATED (element)) &
      CALL set_memory_info (element% field% info                              &
                           ,name, longname ,units                             &
                           ,ldims ,gdims ,ndim ,klev ,ktrac ,alloc ,repr      &
                           ,lpost ,laccu, lmiss, missval, reset, lrerun       &
                           ,contnorest ,table ,code ,bits, leveltype          &
                           ,dimnames, no_default)

  END SUBROUTINE set_stream_element_info
!==============================================================================
  !
  ! Set parameters of list element already created
  ! (private routine within this module)
  !

  SUBROUTINE set_memory_info (info  ,name  ,longname ,units                   &
                           ,ldims ,gdims ,ndim ,klev ,ktrac ,alloc ,repr      &
                           ,lpost ,laccu, lmiss, missval, reset, lrerun       &
                           ,contnorest, table ,code ,bits, leveltype          &
                           ,dimnames, no_default, autodet, verbose,tracidx)
  !--------------------------------------------------------------
  ! Set each parameter in data type memory_info if the respective 
  ! optional parameter is present.
  ! If AUTODET is not false perform consistency check 
  !--------------------------------------------------------------

  TYPE(memory_info),INTENT(inout)        :: info          ! memory info struct.

  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: name          ! variable name
  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: units         ! units
  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: longname      ! long name

  INTEGER          ,INTENT(in) ,OPTIONAL :: ldims (:)     ! local dimensions
  INTEGER          ,INTENT(in) ,OPTIONAL :: gdims (:)     ! global dimensions
  INTEGER          ,INTENT(in) ,OPTIONAL :: ndim          ! rank
  INTEGER          ,INTENT(in) ,OPTIONAL :: klev          ! number of levels
  INTEGER          ,INTENT(in) ,OPTIONAL :: ktrac         ! number of 'tracers'
  LOGICAL          ,INTENT(in) ,OPTIONAL :: alloc         ! pointer allocated
  INTEGER          ,INTENT(in) ,OPTIONAL :: repr          ! representation

  LOGICAL          ,INTENT(in) ,OPTIONAL :: lpost         ! into output stream
  LOGICAL          ,INTENT(in) ,OPTIONAL :: laccu         ! accumulation flag
  LOGICAL          ,INTENT(in) ,OPTIONAL :: lmiss         ! missing value flag
  REAL(dp)         ,INTENT(in) ,OPTIONAL :: missval       ! missing value
  REAL(dp)         ,INTENT(in) ,OPTIONAL :: reset         ! reset value
  LOGICAL          ,INTENT(in) ,OPTIONAL :: lrerun        ! rerun file flag
  LOGICAL          ,INTENT(in) ,OPTIONAL :: contnorest    ! continue on restart

  INTEGER          ,INTENT(in) ,OPTIONAL :: table         ! gribtable number
  INTEGER          ,INTENT(in) ,OPTIONAL :: code          ! gribcode number
  INTEGER          ,INTENT(in) ,OPTIONAL :: bits          ! bits used for GRIB
  INTEGER          ,INTENT(in) ,OPTIONAL :: leveltype     ! grib level type

  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: dimnames(:)   ! dimension names
  LOGICAL          ,INTENT(in) ,OPTIONAL :: no_default    ! use no defaults
  LOGICAL          ,INTENT(in) ,OPTIONAL :: autodet       ! no auto detection 
                                                          !   of properties

  LOGICAL          ,INTENT(in) ,OPTIONAL :: verbose

  INTEGER          ,INTENT(in) ,OPTIONAL :: tracidx

    INTEGER :: isize
    INTEGER :: i
    LOGICAL :: verb
    LOGICAL :: a
    INTEGER :: iklev, iklon, iklat, ikland, iktrac, iknsp, ikcpl
    LOGICAL :: linvalid

    !===================================
    ! set flags from optional parameters
    !===================================
    verb = .FALSE.; IF(PRESENT(verbose)) verb = verbose
    a    = .TRUE.;  IF(PRESENT(autodet)) a    = autodet
    !-------------------------------------------------------
    ! reset to original settings if 'no_default' flag is set
    !-------------------------------------------------------
    IF (PRESENT(no_default)) THEN
      IF (no_default) info = default_info 
    ENDIF
    !-----------------------------------------------------
    ! set components describing the 'Content of the field'
    !-----------------------------------------------------
    CALL assign_if_present (info% name       ,name)
    CALL assign_if_present (info% longname   ,longname)
    CALL assign_if_present (info% units      ,units)
    CALL assign_if_present (info% alloc      ,alloc)
    !---------------------------------
    ! set dimensions on this processor
    !---------------------------------
    IF (PRESENT(ldims)) THEN
      info% DIM(:) = 1
      info% ndim = SIZE (ldims)
      CALL assign_if_present (info% dim     ,ldims)
    ENDIF
    !-----------------------------------
    ! set global dimensions of the field
    !-----------------------------------
    IF (PRESENT(gdims)) THEN
      info% gdim(:) = 1
      info% ndim = SIZE (gdims)
      CALL assign_if_present (info% gdim    ,gdims)
    ENDIF
    !-------------------------------------
    ! explicitely set number of dimensions
    !-------------------------------------
    IF (PRESENT(ndim)) info% ndim = ndim
    isize = info% ndim
    info% gdim(isize+1:) = 1
    info%  DIM(isize+1:) = 1
    !------------------------------
    ! set number of vertical levels
    !------------------------------
    CALL assign_if_present (info% klev   ,klev)
    !--------------
    ! set grid type
    !--------------
    CALL assign_if_present (info% repr   ,repr)
    !-------------------------
    ! set flags concerning I/O
    !-------------------------
    CALL assign_if_present (info% lpost      ,lpost)
    CALL assign_if_present (info% reset      ,reset)
    CALL assign_if_present (info% laccu      ,laccu)
    CALL assign_if_present (info% lmiss      ,lmiss)
    CALL assign_if_present (info% missval    ,missval)
    CALL assign_if_present (info% lrerun     ,lrerun)
    CALL assign_if_present (info% contnorest ,contnorest)

    CALL assign_if_present (info% gribtable  ,table)
    CALL assign_if_present (info% gribcode   ,code)
    CALL assign_if_present (info% gribbits   ,bits)
    CALL assign_if_present (info% levelindx  ,leveltype)
    CALL assign_if_present (info% tracidx    ,tracidx)

    IF (PRESENT(dimnames)) THEN
      DO i = 1, SIZE (dimnames)
       info% IO_var_indx(i) = IO_get_varindx(dimnames(i))
       if (info% IO_var_indx(i)<0) then
         if(present(gdims)) then
           call add_dim (dimnames(i), gdims(i), indx=info% IO_var_indx(i))
         else
           CALL finish ('set_memory_info',&
             trim(info%name)//': unknown dimension:'//dimnames(i))
         endif
       endif
      END DO
    END IF

    !===================
    ! consistency checks
    !===================

    IF (a) THEN
      !-----
      ! rank
      !-----
      IF (info% ndim < 1 .OR. info% ndim > 4) &
        CALL finish ('set_memory_info','ndims == 0')
      !---------------------------
      ! identify lon,lat,lev index
      !---------------------------
      iklev = 0; iklon = 0; iklat = 0; ikland = 0; iktrac = 0; iknsp = 0; ikcpl = 0
      SELECT CASE (info% repr)
      CASE (GAUSSIAN)
        SELECT CASE (info% ndim)
        CASE (2)
          iklon  = 1
          iklat  = 2
        CASE (3)
          iklon  = 1
          iklev  = 2
          iklat  = 3
        CASE (4)
          iklon  = 1
          iklev  = 2
          iktrac = 3
          iklat  = 4
        END SELECT
      CASE (LAND)
          ikland = 1
          IF (info%ndim > 1) iklev = 2
      CASE (SPECTRAL)
        iklev  = 1
        SELECT CASE (info% ndim)
        CASE (3)
          iknsp = 3
          ikcpl = 2
        END SELECT
      CASE (FOURIER)
        iklev  = 1
      END SELECT
      info% iklev = iklev

      !-----------------------------------------
      ! set default values of lon,lat dimensions
      !-----------------------------------------       
      IF (iklon > 0) THEN
        IF (info%  DIM       (iklon) < 0) info%  DIM (iklon) = ldc% nglon
        IF (info% gdim       (iklon) < 0) info% gdim (iklon) = ldc% nlon
        IF (info% IO_var_indx(iklon) < 0) info% IO_var_indx(iklon) = &
          IO_get_varindx('lon')
      ENDIF

      IF (iklat > 0) THEN
        IF (info%  DIM       (iklat) < 0) info%  DIM (iklat) = ldc% nglat
        IF (info% gdim       (iklat) < 0) info% gdim (iklat) = ldc% nlat
        IF (info% IO_var_indx(iklat) < 0) info% IO_var_indx(iklat) = &
          IO_get_varindx('lat')
      ENDIF

      IF (iknsp > 0) THEN
        IF (info% gdim       (iknsp) < 0) info% gdim (iknsp) = nsp
        IF (info% IO_var_indx(iknsp) < 0) info% IO_var_indx(iknsp) = &
          IO_get_varindx('spc')
      ENDIF

      IF (ikland > 0) THEN
        IF (info%  DIM     (ikland) < 0) info%  DIM(ikland) = domain_nland
        IF (info% gdim     (ikland) < 0) info% gdim(ikland) = nland
        IF (info% IO_var_indx(ikland) < 0) info%IO_var_indx(ikland) = &
          IO_get_varindx('landpoint')
      ENDIF

      IF (ikcpl > 0) THEN
        IF (info% gdim       (ikcpl) < 0) info% gdim (ikcpl) = 2
        IF (info% IO_var_indx(ikcpl) < 0) info% IO_var_indx(ikcpl) = &
          IO_get_varindx('complex')
      ENDIF
      IF (iklev > 0) THEN
        !------------------------
        ! If not defined, derive:
        ! leveltype <- io_dim_ids
        !------------------------
        IF (info% levelindx < 0) info% levelindx = info% IO_var_indx (iklev)
        !------------------------
        ! If not defined, derive:
        ! klev <- dim <- gdim
        ! klev <- leveltype
        ! klev <- nlev
        !------------------------
        IF (info%  DIM (iklev) < 0) info%  DIM (iklev) = info% gdim (iklev)
        IF (info% klev         < 0) info% klev         = info% gdim (iklev)
        IF (info% klev < 0 .AND. info% levelindx > 0) &
            info% klev = IO_dim_ids (info% levelindx)% dim_len
        IF (info% klev < 0        ) info% klev         = ldc% nlev
        !-------------------------------
        ! check consistency of leveltype
        !-------------------------------
        IF (info% levelindx <= 0) THEN
          linvalid = .TRUE.
        ELSE
          linvalid = (info% klev /= IO_dim_ids (info% levelindx)% dim_len)
        ENDIF
        !--------------------------------
        ! set leveltype if not set so far
        !--------------------------------
        IF (linvalid) THEN
          IF (info% klev == ldc% nlev) THEN
            info% levelindx = HYBRID
          ELSE IF (info% klev == ldc% nlev+1) THEN
            info% levelindx = HYBRID_H
          ELSE IF (info% klev == 1) THEN
            info% levelindx = SURFACE
          ELSE
            CALL add_unknown_dim (info% klev, info% levelindx)
          ENDIF
        ENDIF
        !-----------------------------------------------
        ! set dim, gdim, io_dim_ids from klev, levelindx
        !-----------------------------------------------
        if (info% DIM(iklev) < 0) &
          info% DIM       (iklev) = info% klev
        info% gdim        (iklev) = info% klev
        info% IO_var_indx (iklev) = info% levelindx
      ELSE ! iklev == 0
        !-----------------------------
        ! surface field, set levelindx
        !-----------------------------
        info% klev = 1
        IF (info% levelindx <= 0) THEN
          linvalid = .TRUE.
        ELSE
          linvalid = (info% klev /= IO_dim_ids (info% levelindx)% dim_len)
        ENDIF
        IF (linvalid) info% levelindx = SURFACE
      ENDIF
      if (iktrac > 0) then
        CALL assign_if_present (info% DIM(iktrac) ,ktrac)      
        CALL assign_if_present (info%gdim(iktrac) ,ktrac)
        IF(info%gdim(iktrac)<0) info%gdim(iktrac) = info% DIM(iktrac) 
      endif
      !-----------------------------------------------
      ! set dimension from netcdf dimension if defined
      !-----------------------------------------------
      DO i=1,info% ndim
       if (info% IO_var_indx (i) > 0) then
        IF(info% DIM(i)<0) info% DIM(i)=io_dim_ids(info%IO_var_indx(i))%dim_len
        IF(info%gdim(i)<0) info%gdim(i)=io_dim_ids(info%IO_var_indx(i))%dim_len
       end if
      end do
      !------------------------------------
      ! finally all dimensions must be set.
      !------------------------------------
      DO i=1,info% ndim
        IF (i == info% iklev) CYCLE
        IF (info%  DIM (i) < 0) CALL finish ('set_memory_info',' dim not set')
        IF (info% gdim (i) < 0) CALL finish ('set_memory_info','gdim not set')
        IF (info% lpost .OR. info% lrerun) THEN
          !-------------------------------------------------------------
          ! add dummy-levelindex for postprocessing and rerun if not set
          !-------------------------------------------------------------
          IF (info% IO_var_indx(i) <= 0) THEN
            linvalid = .TRUE.
          ELSE
            linvalid = info%gdim(i) /= IO_dim_ids(info%IO_var_indx(i))%dim_len
          ENDIF
          IF(linvalid) CALL add_unknown_dim(info% gdim(i),info% IO_var_indx(i))
        ENDIF
      END DO

    ENDIF
    !----------------------------------
    ! check for nproma/ngpblks blocking
    !----------------------------------
    info% lreg = .false.
    info% dima = info% dim
    if (info% repr /= GAUSSIAN .AND. info%repr /= LAND) info% lreg = .true.
    if (ldc % lreg)             info% lreg = .true.
    if (.not.info% lreg .AND. info%repr == GAUSSIAN) then
      select case (info% ndim)
      case (2)
        info% dima(1) = ldc% nproma
        info% dima(2) = ldc% ngpblks
      case (3)
        info% dima(1) = ldc% nproma
        info% dima(3) = ldc% ngpblks
      case (4)
        info% dima(1) = ldc% nproma
        info% dima(4) = ldc% ngpblks
      case default
        info% lreg = .true.
      end select
    endif
    IF (info%repr == LAND) THEN
       info%dima(1) = domain_nland
    ENDIF
    !====================
    ! printout (optional)
    !====================
    IF (verb) CALL print_memory_info (info)

  END SUBROUTINE set_memory_info
!==============================================================================
  !
  ! Create a list new entry
  !
  ! Specific routines for pointers of different rank
  !
  SUBROUTINE add_stream_element_4d (stream, name, ptr, ldims, gdims,       &
             klev, ktrac, units, longname, repr,                           &
             lpost, laccu, lmiss, missval, reset, lrerun, contnorest,      &
             table, code, bits, leveltype,                                 &
             dimnames, mem_info, p4, no_default, verbose)
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 4d-field
  ! optionally overwrite default meta data 
  !
  TYPE (t_stream)  ,INTENT(inout)        :: stream        ! list

  CHARACTER (*)    ,INTENT(in)           :: name          ! name of variable
  REAL(dp)         ,POINTER              :: ptr(:,:,:,:)  ! reference to field
  INTEGER          ,INTENT(in) ,OPTIONAL :: ldims (4)     ! local dimensions
  INTEGER          ,INTENT(in) ,OPTIONAL :: gdims (4)     ! global dimensions
  INTEGER          ,INTENT(in) ,OPTIONAL :: klev          ! number of levels
  INTEGER          ,INTENT(in) ,OPTIONAL :: ktrac         ! number of 'tracers'
  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: units         ! units
  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: longname      ! long name
  INTEGER          ,INTENT(in) ,OPTIONAL :: repr          ! representation
  LOGICAL          ,INTENT(in) ,OPTIONAL :: lpost         ! into output stream
  LOGICAL          ,INTENT(in) ,OPTIONAL :: laccu         ! accumulation flag
  LOGICAL          ,INTENT(in) ,OPTIONAL :: lmiss         ! missing value flag
  REAL(dp)         ,INTENT(in) ,OPTIONAL :: missval       ! missing value
  REAL(dp)         ,INTENT(in) ,OPTIONAL :: reset         ! reset value
  LOGICAL          ,INTENT(in) ,OPTIONAL :: lrerun        ! rerun file flag
  LOGICAL          ,INTENT(in) ,OPTIONAL :: contnorest    ! continue
  INTEGER          ,INTENT(in) ,OPTIONAL :: table         ! gribtable number
  INTEGER          ,INTENT(in) ,OPTIONAL :: code          ! gribcode number
  INTEGER          ,INTENT(in) ,OPTIONAL :: bits          ! bits used for GRIB
  INTEGER          ,INTENT(in) ,OPTIONAL :: leveltype     ! grib level type
  CHARACTER (*)    ,INTENT(in) ,OPTIONAL :: dimnames(4)   ! dimension names
  TYPE(memory_info),POINTER    ,OPTIONAL :: mem_info      ! returned reference
  REAL(dp)         ,POINTER    ,OPTIONAL :: p4(:,:,:,:)   ! provided pointer
  LOGICAL          ,INTENT(in) ,OPTIONAL :: no_default    ! use no defaults
  LOGICAL          ,INTENT(in) ,OPTIONAL :: verbose       ! produce printout

    TYPE (list_element), POINTER :: new_list_element
    LOGICAL                      :: alloc
    INTEGER                      :: istat
    !
    ! add list entry
    !
    CALL add_list_element (stream, new_list_element, code=code)
    !
    ! and set meta data
    !
    alloc = .NOT. PRESENT(p4)
    new_list_element%field%info%alloc = alloc

    CALL set_memory_info (new_list_element% field% info                       &
     ,name=name  ,longname=longname ,units=units                              &
     ,ldims=ldims ,gdims=gdims ,ndim=4 ,klev=klev, ktrac=ktrac ,alloc=alloc   &
     ,repr=repr ,lpost=lpost ,laccu=laccu, lmiss=lmiss, missval=missval       &
     ,reset=reset, lrerun=lrerun                                              &
     ,contnorest=contnorest ,table=table ,code=code ,bits=bits                &
     ,leveltype=leveltype ,dimnames=dimnames                                  &
     ,no_default=no_default, verbose=verbose)

    IF (alloc) THEN
      ALLOCATE (new_list_element%field%      ptr(      &
                new_list_element%field%info% dima(1),  &
                new_list_element%field%info% dima(2),  &
                new_list_element%field%info% dima(3),  &
                new_list_element%field%info% dima(4)), &
                STAT=istat)
      IF (istat /= 0) THEN
        WRITE(message_text,'(a,a)') 'Could not allocate ', TRIM(name)
        CALL finish('add_stream_element(4d)', message_text)
      ENDIF
      stream%memory_used = stream%memory_used &
                 +8*SIZE(new_list_element%field%ptr)
    ELSE
      new_list_element%field%ptr => p4
    ENDIF
    ptr => new_list_element%field%ptr(:,:,:,:)
    IF(PRESENT(mem_info)) mem_info => new_list_element%field%info

    IF ( PRESENT(lmiss) ) THEN
      new_list_element%field%ptr=new_list_element%field%info%missval
    ELSE
      new_list_element%field%ptr=0._dp
    END IF

  END SUBROUTINE add_stream_element_4d
!-----------------------------------------------------------------------------
  SUBROUTINE add_stream_element_3d (stream, name, ptr, ldims, gdims,       &
             klev, units, longname, repr,                                  &
             lpost, laccu, lmiss, missval, reset, lrerun, contnorest,      &
             table, code, bits, leveltype,                                 &
             dimnames, mem_info, p4, no_default, verbose, tracidx)
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 3d-field
  ! optionally overwrite default meta data 
  !
  TYPE (t_stream)  ,INTENT(inout)        :: stream        ! list

  CHARACTER (*)    ,INTENT(in)           :: name          ! name of variable
  REAL(dp)         ,POINTER              :: ptr(:,:,:)    ! reference to field
  INTEGER          ,INTENT(in) ,OPTIONAL :: ldims (3)     ! local dimensions
  INTEGER          ,INTENT(in) ,OPTIONAL :: gdims (3)     ! global dimensions
  INTEGER          ,INTENT(in) ,OPTIONAL :: klev          ! number of levels
  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: units         ! units
  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: longname      ! long name
  INTEGER          ,INTENT(in) ,OPTIONAL :: repr          ! representation
  LOGICAL          ,INTENT(in) ,OPTIONAL :: lpost         ! into output stream
  LOGICAL          ,INTENT(in) ,OPTIONAL :: laccu         ! accumulation flag
  LOGICAL          ,INTENT(in) ,OPTIONAL :: lmiss         ! missing value flag
  REAL(dp)         ,INTENT(in) ,OPTIONAL :: missval       ! missing value
  REAL(dp)         ,INTENT(in) ,OPTIONAL :: reset         ! reset value
  LOGICAL          ,INTENT(in) ,OPTIONAL :: lrerun        ! rerun file flag
  LOGICAL          ,INTENT(in) ,OPTIONAL :: contnorest    ! continue
  INTEGER          ,INTENT(in) ,OPTIONAL :: table         ! gribtable number
  INTEGER          ,INTENT(in) ,OPTIONAL :: code          ! gribcode number
  INTEGER          ,INTENT(in) ,OPTIONAL :: bits          ! bits used for GRIB
  INTEGER          ,INTENT(in) ,OPTIONAL :: leveltype     ! grib level type
  CHARACTER (*)    ,INTENT(in) ,OPTIONAL :: dimnames(3)   ! dimension names
  TYPE(memory_info),POINTER    ,OPTIONAL :: mem_info      ! returned reference
  REAL(dp)         ,POINTER    ,OPTIONAL :: p4(:,:,:,:)   ! provided pointer
  LOGICAL          ,INTENT(in) ,OPTIONAL :: no_default    ! use no defaults
  LOGICAL          ,INTENT(in) ,OPTIONAL :: verbose       ! produce printout
  INTEGER          ,INTENT(in) ,OPTIONAL :: tracidx

    TYPE (list_element), POINTER :: new_list_element
    LOGICAL                      :: alloc
    INTEGER                      :: istat
    !
    ! add list entry
    !
    CALL add_list_element (stream, new_list_element, code=code)
    !
    ! and set meta data
    !

    alloc = .NOT. PRESENT(p4)
    new_list_element%field%info%alloc = alloc

    CALL set_memory_info (new_list_element% field% info                       &
     ,name=name  ,longname=longname ,units=units                              &
     ,ldims=ldims ,gdims=gdims ,ndim=3 ,klev=klev ,alloc=alloc ,repr=repr     &
     ,lpost=lpost ,laccu=laccu, lmiss=lmiss, missval=missval                  &
     ,reset=reset, lrerun=lrerun                                              &
     ,contnorest=contnorest ,table=table ,code=code ,bits=bits                &
     ,leveltype=leveltype ,dimnames=dimnames                                  &
     ,no_default=no_default, verbose=verbose, tracidx=tracidx)

    IF (alloc) THEN
      ALLOCATE (new_list_element%field%      ptr(        &
                new_list_element%field%info% dima(1),    &
                new_list_element%field%info% dima(2),    &
                new_list_element%field%info% dima(3),1), &
                STAT=istat)
      IF (istat /= 0) THEN
        WRITE(message_text,'(a,a)') 'Could not allocate ', TRIM(name)
        CALL finish('add_stream_element(3d)', message_text)
      ENDIF
      stream%memory_used = stream%memory_used &
                 +8*SIZE(new_list_element%field%ptr)
    ELSE
      new_list_element%field%ptr => p4
    ENDIF
    ptr => new_list_element%field%ptr(:,:,:,1)
    IF(PRESENT(mem_info)) mem_info => new_list_element%field%info

    IF ( PRESENT(lmiss) ) THEN
      new_list_element%field%ptr=new_list_element%field%info%missval
    ELSE
      new_list_element%field%ptr=0._dp
    END IF

  END SUBROUTINE add_stream_element_3d
!-----------------------------------------------------------------------------
  SUBROUTINE add_stream_element_2d (stream, name, ptr, ldims, gdims,       &
             klev, units, longname, repr,                                  &
             lpost, laccu, lmiss, missval, reset, lrerun, contnorest,      &
             table, code, bits, leveltype,                                 &
             dimnames, no_default, verbose)
  !
  ! create (allocate) a new table entry
  ! optionally obtain pointer to 2d-field
  ! optionally overwrite default meta data 
  !
  TYPE (t_stream)  ,INTENT(inout)        :: stream        ! list

  CHARACTER (*)    ,INTENT(in)           :: name          ! name of variable
  REAL(dp)         ,POINTER              :: ptr(:,:)      ! reference to field
  INTEGER          ,INTENT(in) ,OPTIONAL :: ldims (2)     ! local dimensions
  INTEGER          ,INTENT(in) ,OPTIONAL :: gdims (2)     ! global dimensions

  INTEGER          ,INTENT(in) ,OPTIONAL :: klev          ! number of levels
  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: units         ! units
  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: longname      ! long name

  INTEGER          ,INTENT(in) ,OPTIONAL :: repr          ! representation

  LOGICAL          ,INTENT(in) ,OPTIONAL :: lpost         ! into output stream
  LOGICAL          ,INTENT(in) ,OPTIONAL :: laccu         ! accumulation flag
  LOGICAL          ,INTENT(in) ,OPTIONAL :: lmiss         ! missing value flag
  REAL(dp)         ,INTENT(in) ,OPTIONAL :: missval       ! missing value
  REAL(dp)         ,INTENT(in) ,OPTIONAL :: reset         ! reset value
  LOGICAL          ,INTENT(in) ,OPTIONAL :: lrerun        ! rerun file flag
  LOGICAL          ,INTENT(in) ,OPTIONAL :: contnorest    ! continue

  INTEGER          ,INTENT(in) ,OPTIONAL :: table         ! gribtable number
  INTEGER          ,INTENT(in) ,OPTIONAL :: code          ! gribcode number
  INTEGER          ,INTENT(in) ,OPTIONAL :: bits          ! bits used for GRIB
  INTEGER          ,INTENT(in) ,OPTIONAL :: leveltype     ! grib level type

  CHARACTER (*)    ,INTENT(in) ,OPTIONAL :: dimnames(2)   ! dimension names
  LOGICAL          ,INTENT(in) ,OPTIONAL :: no_default    ! use no defaults

  LOGICAL          ,INTENT(in) ,OPTIONAL :: verbose       ! produce printout

    TYPE (list_element), POINTER :: element
    LOGICAL                      :: alloc
    INTEGER                      :: istat

    CALL add_list_element (stream, element, code=code)

    alloc = .TRUE.

    CALL set_memory_info (element% field% info                                &
     ,name=name  ,longname=longname ,units=units                              &
     ,ldims=ldims ,gdims=gdims ,ndim=2 ,klev=klev, alloc=alloc ,repr=repr     &
     ,lpost=lpost ,laccu=laccu, lmiss=lmiss, missval=missval                  &
     ,reset=reset, lrerun=lrerun                                              &
     ,contnorest=contnorest ,table=table ,code=code ,bits=bits                &
     ,leveltype=leveltype ,dimnames=dimnames                                  &
     ,no_default=no_default, verbose=verbose)

    ALLOCATE (element%field%ptr (               &
              element%field%info% dima(1),      &
              element%field%info% dima(2),1,1), &
              STAT=istat)
    IF (istat /= 0) THEN
      WRITE(message_text,'(a,a)') 'Could not allocate ', TRIM(name)
      CALL finish('add_stream_element(4d)', message_text)
    ENDIF
    stream%memory_used = stream%memory_used &
               +8*SIZE(element%field%ptr)
    ptr => element%field%ptr(:,:,1,1)

    IF ( PRESENT(lmiss) ) THEN
      ! Note: element%field%ptr = element%field%info%missval crashes for unknown reasons for T63 when called from
      !       add_stream_element(...,soil_diag%soil_temperature,lmiss=T,missing_value=missval) in mo_soil.f90
      element%field%ptr(1:element%field%info% dima(1),1:element%field%info% dima(2),1,1) = element%field%info%missval
    ELSE
      element%field%ptr=0._dp
    END IF

  END SUBROUTINE add_stream_element_2d
!-----------------------------------------------------------------------------
  SUBROUTINE add_stream_element_1d (stream, name, ptr, ldims, gdims, &
             klev, units, longname, repr,                                  &
             lpost, laccu, lmiss, missval, reset, lrerun, contnorest,      &
             table, code, bits, leveltype,                                 &
             dimnames, no_default, verbose)

    TYPE (t_stream)       ,INTENT(inout):: stream      ! list
    CHARACTER (*)         ,INTENT(in)   :: name        ! name of variable
    REAL(dp)              ,POINTER      :: ptr(:)      ! reference to field
    INTEGER      ,OPTIONAL,INTENT(in)   :: ldims(1)    ! shape of array 
    INTEGER      ,OPTIONAL,INTENT(in)   :: gdims(1)    ! global size of field

    INTEGER      ,OPTIONAL,INTENT(in)   :: klev        ! number of levels
    CHARACTER(*) ,OPTIONAL,INTENT(in)   :: units       ! units
    CHARACTER(*) ,OPTIONAL,INTENT(in)   :: longname    ! long name

    INTEGER      ,OPTIONAL,INTENT(in)   :: repr        ! representation

    LOGICAL      ,OPTIONAL,INTENT(in)   :: lpost       ! into output stream
    LOGICAL      ,OPTIONAL,INTENT(in)   :: lmiss       ! missing value flag
    REAL(dp)     ,OPTIONAL,INTENT(in)   :: missval     ! missing value
    REAL(dp)     ,OPTIONAL,INTENT(in)   :: reset       ! reset value

    INTEGER      ,OPTIONAL,INTENT(in)   :: bits        ! bits used for GRIB
    INTEGER      ,OPTIONAL,INTENT(in)   :: leveltype   ! grib level type

    LOGICAL      ,OPTIONAL,INTENT(in)   :: no_default  ! use no defaults

    LOGICAL      ,OPTIONAL,INTENT(in)   :: verbose     ! produce printout
    CHARACTER(*) ,OPTIONAL,INTENT(in)   :: dimnames(1) ! dimension names
    INTEGER      ,OPTIONAL,INTENT(in)   :: code        ! gribcode number
    INTEGER      ,OPTIONAL,INTENT(in)   :: table       ! gribcode table number
    LOGICAL      ,OPTIONAL,INTENT(in)   :: laccu       ! accumulation flag
    LOGICAL      ,OPTIONAL,INTENT(in)   :: lrerun      ! rerun file flag
    LOGICAL      ,OPTIONAL,INTENT(in)   :: contnorest  ! continue on restart
    !
    ! create (allocate) a new table entry
    ! optionally obtain pointer to 1d-field
    ! optionally overwrite default meta data 
    !
!    TYPE (list_element), POINTER :: new_list_element
!    !
!    ! add list entry
!    !
!    CALL add_list_element (stream, new_list_element, code=code)
!    !
!    ! and set meta data
!    !
!    new_list_element%field%info%name = name
!    ALLOCATE (new_list_element%field%ptr(ldims(1),1,1,1))
!
!    new_list_element%field%ptr=0._dp
!
!    stream%memory_used = stream%memory_used &
!               +8*SIZE(new_list_element%field%ptr)
!
!    new_list_element%field%info%alloc = .TRUE.
!
!    new_list_element%field%info%dim(1) = ldims(1)
!    new_list_element%field%info%dim(2) = 1
!    new_list_element%field%info%dim(3) = 1
!    new_list_element%field%info%dim(4) = 1
!    new_list_element%field%info%dima   = new_list_element%field%info%dim
!    new_list_element%field%info%lreg   = .true.
!
!    new_list_element%field%info%gdim(1) = gdims(1)
!    new_list_element%field%info%gdim(2) = 1
!    new_list_element%field%info%gdim(3) = 1
!    new_list_element%field%info%gdim(4) = 1
!
!    new_list_element%field%info%ndim = 1
!
!    ptr => new_list_element%field%ptr(:,1,1,1)
!    !
!    ! pass optional arguments
!    !
!    IF (PRESENT(dimnames)) THEN
!       new_list_element%field%info%IO_var_indx(1) = IO_get_varindx(dimnames(1))
!    END IF
!
!    IF (PRESENT(table))     new_list_element%field%info%gribtable  = table 
!    IF (PRESENT(code))      new_list_element%field%info%gribcode   = code
!    IF (PRESENT(laccu))     new_list_element%field%info%laccu      = laccu
!    IF (PRESENT(lrerun))    new_list_element%field%info%lrerun     = lrerun
!    IF (PRESENT(contnorest))new_list_element%field%info%contnorest = contnorest

    TYPE (list_element), POINTER :: element
    LOGICAL                      :: alloc
    INTEGER                      :: istat

    CALL add_list_element (stream, element, code=code)

    alloc = .TRUE.

    CALL set_memory_info (element% field% info                                &
     ,name=name  ,longname=longname ,units=units                              &
     ,ldims=ldims ,gdims=gdims ,ndim=1 ,klev=klev, alloc=alloc ,repr=repr     &
     ,lpost=lpost ,laccu=laccu, lmiss=lmiss, missval=missval                  &
     ,reset=reset, lrerun=lrerun                                              &
     ,contnorest=contnorest ,table=table ,code=code ,bits=bits                &
     ,leveltype=leveltype ,dimnames=dimnames                                  &
     ,no_default=no_default, verbose=verbose)

    IF (element%field%info%repr /= LAND) THEN
       element%field%info%levelindx = SURFACE
    ENDIF

    ALLOCATE (element%field%ptr (                 &
              element%field%info% dima(1),1,1,1), &                
              STAT=istat)
    IF (istat /= 0) THEN
      WRITE(message_text,'(a,a)') 'Could not allocate ', TRIM(name)
      CALL finish('add_stream_element(1d)', message_text)
    ENDIF
    stream%memory_used = stream%memory_used &
               +8*SIZE(element%field%ptr)
    ptr => element%field%ptr(:,1,1,1)

    IF ( PRESENT(lmiss) ) THEN
      element%field%ptr=element%field%info%missval
    ELSE
      element%field%ptr=0._dp
    END IF

  END SUBROUTINE add_stream_element_1d
!==============================================================================
  subroutine remove_stream_element (this_list, name)
  TYPE (t_stream)  ,INTENT(inout) :: this_list
  character(len=*) ,INTENT(in)    :: name   
  !
  ! remove one element from the list
  ! the element is identified by its name
  !
    TYPE (list_element) ,POINTER :: p

    if (this_list% first_list_element% field% info% name == name) then
      call remove_list_element (this_list, this_list% first_list_element)
      return
    else
      p => this_list% first_list_element
      do
        if (.not.associated (p% next_list_element)) exit
        if (p% next_list_element% field% info% name == name) then
          call remove_list_element (this_list, p% next_list_element)
          exit
        endif
        p => p% next_list_element
      end do
    endif
  end subroutine remove_stream_element
!==============================================================================

  !
  ! Get element of a list, specific routines for pointers of different rank
  !

  SUBROUTINE get_stream_element_4d (stream, name, ptr)

    TYPE (t_stream) ,INTENT(in) :: stream       ! list
    CHARACTER (*)   ,INTENT(in) :: name         ! name of variable
    REAL(dp)        ,POINTER    :: ptr(:,:,:,:) ! reference to allocated field

    TYPE (list_element), POINTER :: element

    ! obtain pointer to 4d-field

    element => find_list_element (stream, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%ptr

  END SUBROUTINE get_stream_element_4d
!-----------------------------------------------------------------------------
  SUBROUTINE get_stream_element_3d (stream, name, ptr)

    TYPE (t_stream) ,INTENT(in) :: stream       ! list
    CHARACTER (*)   ,INTENT(in) :: name         ! name of variable
    REAL(dp)        ,POINTER    :: ptr(:,:,:)   ! reference to allocated field

    TYPE (list_element), POINTER :: element
    !
    ! obtain pointer to 3d-field
    !
    element => find_list_element (stream, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%ptr(:,:,:,1)

  END SUBROUTINE get_stream_element_3d
!-----------------------------------------------------------------------------
  SUBROUTINE get_stream_element_2d (stream, name, ptr)

    TYPE (t_stream)   ,INTENT(in) :: stream     ! list
    CHARACTER (*)     ,INTENT(in) :: name       ! name of variable
    REAL(dp)          ,POINTER    :: ptr(:,:)   ! reference to allocated field

    TYPE (list_element), POINTER :: element
    !
    ! obtain pointer to 2d-field
    !
    element => find_list_element (stream, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%ptr(:,:,1,1)

  END SUBROUTINE get_stream_element_2d
!-----------------------------------------------------------------------------
  SUBROUTINE get_stream_element_1d (stream, name, ptr)

    TYPE (t_stream) ,INTENT(in) :: stream    ! list
    CHARACTER (*)   ,INTENT(in) :: name      ! name of variable
    REAL(dp)        ,POINTER    :: ptr(:)    ! reference to allocated field

    TYPE (list_element), POINTER :: element
    !
    ! obtain pointer to 1d-field
    !
    element => find_list_element (stream, name)
    NULLIFY (ptr)
    IF (ASSOCIATED (element)) ptr => element%field%ptr(:,1,1,1)

  END SUBROUTINE get_stream_element_1d
!==============================================================================
  !
  ! Print routines for control output and debuggung
  !

  SUBROUTINE print_memory_use (stream)

    TYPE (t_stream) ,INTENT(in) :: stream ! list

    WRITE (message_text,'(a16,a,a,i10,a,i4,a)')                 &
         ADJUSTR(stream%name), '-buffer: ',                     &
         'Memory in use: ', stream%memory_used/1024, ' kb in ', &
         stream%list_elements, ' fields.'
    CALL message('',message_text,adjust_right=.TRUE.)

  END SUBROUTINE print_memory_use
!-----------------------------------------------------------------------------
  SUBROUTINE print_memory_table (stream)

    TYPE (t_stream),  INTENT(in) :: stream ! list
    !
    ! print current memory table 
    !
    CALL message('','')
    CALL message('','')
    CALL message('','Status of base memory:')    
    CALL message('','')

    CALL print_linked_list (stream)
    
  END SUBROUTINE print_memory_table
!-----------------------------------------------------------------------------
  SUBROUTINE print_sinfo (stream)

    TYPE (t_stream),  INTENT(in) :: stream ! list
    !
    ! print current stat table 
    !
    WRITE (message_text,'(a16,a)') TRIM(stream%name), '-buffer: '
    CALL message('',message_text)
    CALL message('','')    
    CALL message('','')
    CALL message('','Statistic of base memory:')
    CALL message('','')

    CALL print_sinfo_list (stream)
    
  END SUBROUTINE print_sinfo
!==============================================================================
  SUBROUTINE add_stream_reference (stream, name, fromstream, lpost, kprec)
  !
  ! add supplementary fields (eg. geopotential, surface pressure, gridbox area)
  !
  TYPE (t_stream)  ,INTENT(inout)           :: stream
  CHARACTER(len=*) ,INTENT(in)              :: name
  CHARACTER(len=*) ,INTENT(in)    ,OPTIONAL :: fromstream
  LOGICAL          ,INTENT(in)    ,OPTIONAL :: lpost
  INTEGER          ,INTENT(in)    ,OPTIONAL :: kprec

    TYPE (memory_type)  ,POINTER :: source
    TYPE (list_element) ,POINTER :: new_list_element

    CALL locate (source ,name ,fromstream)
    IF(ASSOCIATED(source)) THEN
      IF (source% info% gribcode > 0) THEN
        CALL add_list_element (stream, new_list_element, &
                               code=source% info% gribcode)
      ELSE
        CALL add_list_element (stream, new_list_element)
      ENDIF
      new_list_element% field                = source
      new_list_element% field% info% alloc   = .FALSE.
      new_list_element% field% info% lrerun  = .FALSE.
      IF (PRESENT(lpost)) &
        new_list_element% field% info% lpost = lpost
      IF (PRESENT(kprec)) &
        new_list_element% field% info% gribbits = kprec
    ENDIF

  END SUBROUTINE add_stream_reference
!------------------------------------------------------------------------------
  SUBROUTINE locate (ENTRY, name, stream)
  !
  ! find an entry
  !
  TYPE(memory_type) ,POINTER              :: ENTRY
  CHARACTER(len=*)  ,INTENT(in)           :: name
  CHARACTER(len=*)  ,INTENT(in) ,OPTIONAL :: stream

    INTEGER                     :: i
    TYPE(list_element) ,POINTER :: link

    NULLIFY (ENTRY)

    DO i=1,nstreams
      IF (PRESENT(stream)) THEN
        IF (stream /= ostreams(i)% name) CYCLE
      ENDIF
      link => find_list_element (ostreams(i), name)
      IF (ASSOCIATED(link)) THEN
        ENTRY => link% field
        EXIT
      ENDIF
    END DO

  END SUBROUTINE locate
!==============================================================================
  !
  ! private routines to assign values if actual parameters are present
  !
!------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_character (y,x)
  CHARACTER(len=*) ,INTENT(inout)        :: y
  CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: x
    IF (.NOT. PRESENT(x)) RETURN
    IF ( x == ' ' )       RETURN      
    y = x
  END SUBROUTINE assign_if_present_character
!------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_logical (y,x)
  LOGICAL ,INTENT(inout)        :: y
  LOGICAL ,INTENT(in) ,OPTIONAL :: x
    IF (PRESENT(x)) y = x
  END SUBROUTINE assign_if_present_logical
!------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_integer (y,x)
  INTEGER ,INTENT(inout)        :: y
  INTEGER ,INTENT(in) ,OPTIONAL :: x
    IF (.NOT. PRESENT(x)) RETURN
    IF ( x == -HUGE(x)  ) RETURN
    y = x
  END SUBROUTINE assign_if_present_integer
!------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_integers (y,x)
  INTEGER ,INTENT(inout)        :: y (:)
  INTEGER ,INTENT(in) ,OPTIONAL :: x (:)
    INTEGER :: n
    IF (PRESENT(x)) THEN
      n = MIN(SIZE(x), SIZE(y))
      y(1:n) = x(1:n)
    ENDIF
  END SUBROUTINE assign_if_present_integers
!------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_real (y,x)
  REAL(dp), INTENT(inout)        :: y
  REAL(dp) ,INTENT(in) ,OPTIONAL :: x
    IF (.NOT.PRESENT(x)) RETURN
    IF ( x == -HUGE(x) ) RETURN
    y = x
  END SUBROUTINE assign_if_present_real
!==============================================================================
  CHARACTER(len=3) FUNCTION suffix (unit)
  INTEGER, INTENT(in) :: unit
    !
    ! preliminary: convert a unit number to a file suffix ('.ii')
    !
    suffix(1:1) = '.'
    suffix(2:2) = CHAR (    unit/10 + ICHAR('0'))
    suffix(3:3) = CHAR (MOD(unit,10)+ ICHAR('0'))
  END FUNCTION suffix
!==============================================================================
END MODULE mo_memory_base
