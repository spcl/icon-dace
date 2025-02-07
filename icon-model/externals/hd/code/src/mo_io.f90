! mo_io.f90 - IO routines and related variables
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!_________________________________________

MODULE mo_io
  !
  !  This routine originates (year 2014) from MPI-ESM, the Earth System Model of the 
  !  Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
  !  Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
  !  version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
  !  doi: 10.1029/2018MS001400.

  USE mo_kind,          ONLY: dp
  USE mo_netcdf
  USE mo_filename,      ONLY: standard_rerun_file, rerun_filetype, &
                              NONE, NETCDF, NETCDF2
!OTBHDF                              NONE, NETCDF, NETCDF2, NETCDF4
  USE mo_exception,     ONLY: message_text, message, finish
  USE mo_control,       ONLY: ldebugio, lindependent_read, lcollective_write
  USE mo_mpi,           ONLY: p_io, p_pe, p_bcast, p_all_comm, p_info_null
  USE mo_util_string,   ONLY: separator, toupper
  USE mo_linked_list,   ONLY: list_element, FOURIER, GAUSSIAN, SPECTRAL, LAND
  USE mo_decomposition, ONLY: dcg => global_decomposition
!OSBTRGT  USE mo_tr_gather,     ONLY: gather_field, gather_spectral, gather_fourier
!OSBTRSC  USE mo_tr_scatter,    ONLY: scatter_spectral, scatter_field, scatter_fourier
  USE mo_util_string,   ONLY: toupper
  USE mo_time_control,  ONLY: delta_time, next_date, resume_date,     &
                              start_date, init_step, lresume,         &
                              set_delta_time, get_time_step,          &
                              out_convert_date, inp_convert_date,     &
                              ec_manager_init, write_date, lfirst_cycle

  IMPLICIT NONE

  PUBLIC

  TYPE (FILE_INFO), SAVE         :: yearnc
  TYPE (FILE_INFO), SAVE         :: ini_ozon
  TYPE (FILE_INFO), SAVE         :: ini_field
  TYPE (FILE_INFO), SAVE         :: sstnc0, sstnc1, sstnc2
  TYPE (FILE_INFO), SAVE         :: icenc0, icenc1, icenc2
  TYPE (FILE_INFO), SAVE         :: flunc1
  TYPE (FILE_INFO), SAVE, TARGET :: header
  TYPE (FILE_INFO), SAVE         :: ini_surf
  TYPE (FILE_INFO), SAVE         :: ini_spec
  TYPE (FILE_INFO), SAVE         :: restart(31:38)
  TYPE (FILE_INFO), SAVE         :: ini_vol

  REAL(dp), ALLOCATABLE :: vlon(:), vlat(:)

  ! IO variables

  INTEGER, PARAMETER :: IO_READ = 1, IO_WRITE = 2

  INTEGER :: IO_file_id, IO_var_id, IO_dim_id
  INTEGER :: IO_dims(4)
  INTEGER :: IO_timestep = 0
  INTEGER :: istep

  ! time variables

  INTEGER :: forecast_date     ! YYYYMMDD
  INTEGER :: forecast_time     ! HHMMSS
  INTEGER :: verification_date ! YYYYMMDD
  INTEGER :: verification_time ! HHMMSS

  ! global land surface fraction (formerly defined in mo_couple.f90)

  REAL(dp) :: slm_glob

  ! independent read of file from each processor

  LOGICAL :: l_read = .FALSE.

  CHARACTER (len= 80) :: label(8)

  !-----------------------------------------------------------------------------
CONTAINS
  !-----------------------------------------------------------------------------
  INTEGER FUNCTION inquire_rerun_file(filename)

    CHARACTER(len=*), INTENT(IN) :: filename
    
    INTEGER :: fileID, status
    
    status = NF_OPEN(filename, NF_NOWRITE, fileID)

    IF ( status == NF_NOERR ) THEN
      ! file is NETCDF
      status = NF_CLOSE(fileID)
      inquire_rerun_file = NETCDF
    ELSE
      CALL message  ('inquire_rerun_file', &
           TRIM(filename)//' - '//TRIM(NF_STRERROR(status)) )
      inquire_rerun_file = NONE
    ENDIF

  END FUNCTION inquire_rerun_file
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_close(fileinfo)

    TYPE (FILE_INFO), INTENT(INOUT)  :: fileinfo
    INTEGER :: status

    IF (ldebugio) THEN
      WRITE (message_text,'(a,a,i0)') &
           fileinfo%file_name(1:30), ' Id = ',fileinfo%file_id
      CALL message('IO_close', message_text)
    ENDIF
    
    IF (fileinfo%opened) THEN
      status = NF_CLOSE(fileinfo%file_id)
      
      IF (status == NF_NOERR) THEN
        fileinfo%opened = .FALSE.
      ELSE
        IF ( fileinfo%format == NETCDF ) THEN
          CALL message ('IO_close', NF_STRERROR(status))
        ENDIF
        CALL finish  ('IO_close', 'Run terminated.')
      ENDIF
    ELSE
      CALL message ('IO_close','file is not open: '//TRIM(fileinfo%file_name))
    ENDIF
    
    IO_dim_ids(:)%dim_id = -1

  END SUBROUTINE IO_close
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_open(filename, fileinfo, mode, iostat, lparallel)

    INTEGER,           INTENT(IN)    :: mode
    CHARACTER(len=*),  INTENT(IN)    :: filename
    TYPE (FILE_INFO),  INTENT(INOUT) :: fileinfo
    INTEGER, OPTIONAL, INTENT(OUT)   :: iostat
    LOGICAL, OPTIONAL, INTENT(in)    :: lparallel

    INTEGER :: status, fill_mode, write_mode

    LOGICAL :: lparallel_open

    IF (PRESENT (iostat)) iostat = 0

    IF (PRESENT(lparallel)) THEN
      lparallel_open = lparallel
    ELSE
      lparallel_open = .FALSE.
    ENDIF

    IF (fileinfo%opened) THEN
      CALL message ('IO_open', 'file '//TRIM(fileinfo%file_name)//' already open')
      CALL finish  ('IO_open', 'Run terminated.')
    ENDIF
    
    IF (mode /= IO_READ .AND. mode /= IO_WRITE) THEN
      CALL message ('IO_open', 'unexpected mode')
      CALL finish  ('IO_open', 'Run terminated.')
    ENDIF
    
    ! clean up all elements of netCDF structure
    
    CALL IO_info_construct(fileinfo)
    
    fileinfo%file_name = filename
    
    IF (mode == IO_READ) THEN
print*,'OTBIO_open in ', trim(filename), NF_NOWRITE, fileinfo%file_name
!      status = NF_OPEN (TRIM(filename), NF_NOWRITE, fileinfo%file_id)
       status = NF__OPEN (TRIM(filename), NF_NOWRITE, chunksize, fileinfo%file_id)
!      fileinfo%file_id = ncopn(TRIM(filename), NF_NOWRITE,status)
print*,'OTBIO_open out ', trim(filename), status, fileinfo%file_id
    ELSE
      write_mode = NF_CLOBBER
      IF (lparallel_open) THEN
!OTBHDF        IF ( rerun_filetype /= NETCDF4 ) THEN 
          CALL finish('IO_open', 'Parallel write requires NetCDF4 format')
!OTBHDF        ENDIF
#ifdef PARALLEL_NC4
        status = NF_CREATE_PAR(filename, write_mode,    &
                               p_all_comm, p_info_null, &
                               fileinfo%file_id)
        fileinfo%parallel = .TRUE.
#else
        CALL finish('IO_open', 'Parallel write requires -DPARALLEL_NC4 enabled during compilation')
#endif
      ELSE
        IF ( rerun_filetype == NETCDF2 ) THEN
          write_mode = write_mode + NF_64BIT_OFFSET
        ENDIF
!OTBHDF        IF ( rerun_filetype == NETCDF4 ) THEN
!OTBHDF          write_mode = write_mode + NF_NETCDF4
!OTBHDF        ENDIF
        status = NF__CREATE (filename, write_mode,   &
                             initialsize, chunksize, & 
                             fileinfo%file_id)
      ENDIF
    ENDIF
    
    IF (ldebugio) THEN
      WRITE (message_text,'(a,a,i0,a,i0)') &
           TRIM(filename), ' mode = ', mode, ' Id = ',fileinfo%file_id
      CALL message('IO_open', message_text)
    ENDIF
    
    IF (status == NF_NOERR) THEN
      fileinfo%opened = .TRUE.
    ELSE
      CALL message ('IO_open failed on', TRIM(filename))
      IF ( fileinfo%format == NETCDF ) THEN
        CALL message ('IO_open', NF_STRERROR(status))
      ELSE
        WRITE(message_text,'(a,i0)') 'IO_open  : error status = ', status
        CALL message('', message_text)
      ENDIF
      IF (PRESENT (iostat)) THEN
        iostat = status
      ELSE
        CALL finish  ('IO_open', 'Run terminated.')
      ENDIF
    ENDIF
    ! successfull open, switch fill mode off.
 
    status = NF_SET_FILL (fileinfo%file_id, NF_NOFILL, fill_mode)    

    ! only for detailed debugging ...
    ! CALL IO_info_print(fileinfo)

  END SUBROUTINE IO_open
!==============================================================================
  SUBROUTINE IO_read_header(fileinfo)
  !
  ! read a NetCDF header information
  !        read into argument        'fileinfo' 
  !        and check for initial or restart file
  !
!OSBVRSN    USE mo_version, ONLY: label
    
    TYPE (FILE_INFO), INTENT(INOUT) :: fileinfo
    INTEGER :: fileID
    INTEGER :: i
    CHARACTER(len=80) :: clabel(SIZE(label))
    CHARACTER(len=16) :: attribute

#ifdef OLD_READ
    IF (p_pe == p_io) THEN
#endif

    fileID = fileinfo%file_id
    CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, 'file_type', &
         fileinfo%file_type)
    !
    ! check for initial or restart file
    !
    IF ( fileinfo%file_type(1:3) /= 'Ini' .AND. &
         fileinfo%file_type(1:3) /= 'Res' ) THEN
      CALL message ('IO_read_header', fileinfo%file_type)
      CALL message ('IO_read_header', 'No ECHAM initial or restart file.')
      CALL finish  ('IO_read_header', 'Run terminated.')
    ENDIF

 !   CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, 'title',       fileinfo% title)
    CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, 'file_type',   fileinfo% file_type)
    CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, 'source_type', fileinfo% binary_source)
    CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, 'history',     fileinfo% creation_program)
    CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, 'user',        fileinfo% creation_user)
    CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, 'created',     fileinfo% creation_date)
    
    DO i = 1, SIZE(clabel)
      clabel(i) = ' '
    ENDDO
    
    DO i = 1, SIZE(clabel)
      WRITE (attribute,'(a,i0)') 'label_', i
      CALL IO_GET_ATT_TEXT (fileID, NF_GLOBAL, TRIM(attribute), clabel(i))
    ENDDO
    
    IF (ldebugio) THEN
      CALL message ('','')
      CALL message('',TRIM(fileinfo%file_type))
      CALL message('',TRIM(fileinfo%binary_source))
      CALL message('',TRIM(fileinfo%creation_program))
      CALL message('',TRIM(fileinfo%creation_user))
      CALL message('',TRIM(fileinfo%creation_date))
      DO i = 1, SIZE(clabel)
        CALL message ('',TRIM(clabel(i)))
      ENDDO
      CALL message ('','')
    ENDIF

#ifdef OLD_READ
    ENDIF
#endif

  END SUBROUTINE IO_read_header
!==============================================================================
  SUBROUTINE IO_init

    USE mo_control,       ONLY: ngl, nhgl, nlon, nlp2, nlev, nlevp1, nsp,   &
                                nvclev, nmp1, vct, nm, nn, nk, n2sp, n4mp1, &
                                n2mp1, nnp1, nkp1, nisp
    USE mo_filename,      ONLY: out_expname

    INTEGER :: status
    INTEGER :: fileID

    IF (p_pe == p_io .OR. lindependent_read) THEN
      l_read = .TRUE.
    ENDIF

    yearnc%opened = .FALSE.
    sstnc0%opened = .FALSE.
    sstnc1%opened = .FALSE.
    sstnc2%opened = .FALSE.
    icenc0%opened = .FALSE.
    icenc1%opened = .FALSE.
    icenc2%opened = .FALSE.
    flunc1%opened = .FALSE.
    header%opened = .FALSE.
    ini_surf%opened = .FALSE.
    ini_spec%opened = .FALSE.
    ini_ozon%opened = .FALSE.
    ini_field%opened = .FALSE.
    restart(31:38)%opened = .FALSE.
    ini_vol%opened = .FALSE.

    ! open rerun or initial file

    IF(out_expname == ' ') &
         CALL finish(' ','specify out_expname in namelist runctl !')
    IF (lresume) THEN
      header%format = inquire_rerun_file('rerun_'//TRIM(out_expname)//'_echam')
      IF (header% format == NONE) &
           CALL finish(' ','restart file not available !')        
      IF (l_read) &
           CALL io_open ('rerun_'//TRIM(out_expname)//'_echam',header, IO_READ)
    ELSE
      header%format = NETCDF
      IF (l_read) &
      CALL IO_open_unit(nisp, header, IO_READ)
    ENDIF

    IF (l_read) &
      CALL IO_read_header(header)

    IF (l_read) THEN
      fileID = header%file_id

      ! get data reference times

        CALL IO_GET_ATT_INT (fileID, NF_GLOBAL, 'fdate', forecast_date)
        CALL IO_GET_ATT_INT (fileID, NF_GLOBAL, 'ftime', forecast_time)
        CALL IO_GET_ATT_INT (fileID, NF_GLOBAL, 'vdate', verification_date)
        CALL IO_GET_ATT_INT (fileID, NF_GLOBAL, 'vtime', verification_time)
    ENDIF

    IF (.NOT. lindependent_read) THEN
      CALL p_bcast (forecast_date, p_io)
      CALL p_bcast (forecast_time, p_io)
      CALL p_bcast (verification_date, p_io)
      CALL p_bcast (verification_time, p_io)
    ENDIF

    !*** date/time setting
    CALL inp_convert_date(forecast_date, forecast_time, start_date)
    CALL inp_convert_date(verification_date, verification_time, resume_date)

    IF (l_read) THEN

      CALL write_date(start_date,'Read date (initial/restart): ')

      ! get spherical truncations ...

      CALL IO_GET_ATT_INT (fileID, NF_GLOBAL, 'spherical_truncation_n', nn)
      CALL IO_GET_ATT_INT (fileID, NF_GLOBAL, 'spherical_truncation_m', nm)
      CALL IO_GET_ATT_INT (fileID, NF_GLOBAL, 'spherical_truncation_k', nk)

      IF (lresume) THEN
        !
        ! get nstep
        !
        CALL IO_GET_ATT_INT (fileID, NF_GLOBAL, 'nstep', istep)
        !
        ! get timestep
        !
        CALL IO_GET_ATT_INT (fileID, NF_GLOBAL, 'timestep', IO_timestep)
        
        CALL IO_GET_ATT_DOUBLE (fileID, NF_GLOBAL, 'land_fraction_global', slm_glob)
      ELSE
        istep = INIT_STEP
      ENDIF

      !
      ! inquire for dimensions and get values
      !
      CALL IO_INQ_DIMID  (fileID, 'lat', IO_dim_id)
      CALL IO_INQ_DIMLEN (fileID, IO_dim_id, ngl)
      CALL IO_INQ_DIMID  (fileID, 'lon', IO_dim_id)
      CALL IO_INQ_DIMLEN (fileID, IO_dim_id, nlon)
      !
      ! levels are named 'lev' or 'nlev'   (old style)
      ! spectral coeff. are 'spc' or 'nsp' (old style)
      !
      status = NF_INQ_DIMID (fileID, 'lev', IO_dim_id)
      IF (status /= NF_NOERR) &
           CALL IO_INQ_DIMID   (fileID, 'nlev', IO_dim_id)
      CALL   IO_INQ_DIMLEN  (fileID, IO_dim_id, nlev)
      status = NF_INQ_DIMID (fileID, 'spc', IO_dim_id)
      IF (status /= NF_NOERR) &
           CALL IO_INQ_DIMID   (fileID, 'nsp', IO_dim_id)
      CALL   IO_INQ_DIMLEN  (fileID, IO_dim_id, nsp)
      CALL   IO_INQ_DIMID   (fileID, 'nvclev', IO_dim_id)
      CALL   IO_INQ_DIMLEN  (fileID, IO_dim_id, nvclev)
      
      !        IF (lresume) CALL IO_INQ_DIMID  (fileID, 'nhtrac', IO_dim_id)
      !        IF (lresume) CALL IO_INQ_DIMLEN (fileID, IO_dim_id, nhtrac)
      
    ENDIF

    IF (.NOT. lindependent_read) THEN
      CALL p_bcast (nn, p_io)
      CALL p_bcast (nm, p_io)
      CALL p_bcast (nk, p_io)
      CALL p_bcast (istep, p_io)
      CALL p_bcast (IO_timestep, p_io)
      CALL p_bcast (slm_glob, p_io)
      CALL p_bcast (ngl, p_io)
      CALL p_bcast (nlon, p_io)
      CALL p_bcast (nlev, p_io)
      CALL p_bcast (nsp, p_io)
      CALL p_bcast (nvclev, p_io)  
    ENDIF

    ! derive dependend dimensions 

    nkp1   = nk+1
    nmp1   = nm+1
    nnp1   = nn+1
    n2mp1  = nmp1+nmp1
    n4mp1  = n2mp1+n2mp1
    nlevp1 = nlev+1
    nhgl   = ngl/2
    nlp2   = nlon+2
    n2sp   = nsp+nsp

    IF (.NOT.lresume) IO_timestep = set_delta_time(nn,nlev)

    ! read lon

    IF (.NOT. ALLOCATED(vlon)) ALLOCATE (vlon(nlon))
    IF (l_read) THEN
      CALL IO_INQ_VARID (fileID, 'lon', IO_var_id)
      CALL IO_GET_VAR_DOUBLE (fileID, IO_var_id, vlon)
    ENDIF
    IF (.NOT. lindependent_read) THEN
      CALL p_bcast (vlon, p_io)
    ENDIF

    ! read lat

    IF (.NOT. ALLOCATED(vlat)) ALLOCATE (vlat(ngl))
    IF (l_read) THEN
      CALL IO_INQ_VARID (fileID, 'lat', IO_var_id)
      CALL IO_GET_VAR_DOUBLE (fileID, IO_var_id, vlat)
    ENDIF
    IF (.NOT. lindependent_read) THEN
      CALL p_bcast (vlat, p_io)
    ENDIF

    ! read vct 

    IF (lfirst_cycle) NULLIFY(vct)
    IF (.NOT. ASSOCIATED(vct)) ALLOCATE (vct(nvclev*2))
    IF (l_read) THEN
      CALL IO_INQ_VARID (fileID, 'vct_a', IO_var_id)
      CALL IO_GET_VAR_DOUBLE (fileID, IO_var_id, vct(1:nvclev))
      
      CALL IO_INQ_VARID (fileID, 'vct_b', IO_var_id)   
      CALL IO_GET_VAR_DOUBLE (fileID, IO_var_id, vct(nvclev+1:2*nvclev))
    ENDIF
    IF (.NOT. lindependent_read) THEN
      CALL p_bcast (vct, p_io)
    ENDIF
    
    ! initialize echam time manager
    
    CALL ec_manager_init(IO_timestep,istep)
    
    IF (l_read) &
         CALL IO_close(header)

    ! Allocate arrays which depend on nlev
    
print*,'OTBPR IO_init_dims'
    CALL IO_init_dims
    
  END SUBROUTINE IO_init
!==============================================================================
  !
  ! Write (restart) files
  !
  SUBROUTINE write_streams

    USE mo_memory_base, ONLY: ostreams, nstreams
  
    LOGICAL           :: out_stream  (SIZE(ostreams)) ! streams to write
    LOGICAL           :: in_this_file(SIZE(ostreams)) ! streams to this file
    INTEGER           :: i, j                         ! loop indices
    TYPE(FILE_INFO)   :: fileinfo                     ! 

    !
    ! mark all streams to be written
    !
    out_stream = ostreams% lrerun 
    !
!OSBUTLS    CALL message ('',separator)
    CALL message ('','')
    CALL message ('','writing restart files :')
    CALL message ('','')
    SELECT CASE (rerun_filetype)
!OTBHDF    CASE (NETCDF4)
!OTBHDF      CALL message ('','          file format : NETCDF4')
    CASE (NETCDF2)
      CALL message ('','          file format : NETCDF2')
    CASE (NETCDF)
      CALL message ('','          file format : NETCDF')
    END SELECT
    CALL message ('','')
    CALL message ('','file suffix   buffer')
    !
    ! loop over files to process
    !
    DO i = 1,nstreams
      IF (.NOT. out_stream(i)) CYCLE
      !
      ! identify all streams to be written to this file
      !
      in_this_file = out_stream .AND. ostreams(i)%rest_suf == ostreams%rest_suf
      !
      ! open file
      !
      fileinfo%opened = .FALSE.
      fileinfo%format = rerun_filetype
      IF (p_pe == p_io .OR. lcollective_write) &
           CALL io_open (filename = TRIM(standard_rerun_file)      &
                                  //TRIM(ostreams(i)% rest_suf),   &
                         fileinfo = fileinfo, mode     = IO_WRITE, &
                         lparallel = lcollective_write)
      !
      ! write header
      !
!OSBIOWR      CALL io_write_header1 (fileinfo, lcollective_write)
      !
      DO j=i,nstreams
        IF (in_this_file(j)) THEN
!OSBIOWR            CALL io_write_header2 (fileinfo, ostreams(j)% first_list_element, lcollective_write)
        ENDIF
      END DO
      !
!OSBIOWR      CALL io_write_header3(fileinfo, lcollective_write)
      !
      ! write variables
      !
      DO j=i,nstreams
        IF (in_this_file(j)) THEN
          IF (p_pe == p_io) &
            WRITE(message_text,'(5x,a,1x,a)') ostreams(i)% rest_suf, ostreams(j)% name
            CALL message ('',message_text)
!OSBIOWR            CALL io_write_stream (fileinfo, ostreams(j)% first_list_element)
        ENDIF
      END DO
      !
      ! close file
      !
      IF (p_pe == p_io .OR. lcollective_write) &
           CALL io_close (fileinfo)
      !
      ! mark streams written to this file
      !
      out_stream = out_stream .AND. .NOT. in_this_file
    END DO
    IF (p_pe == p_io) THEN
      CALL message ('','')
!OSBUTLS      CALL message ('',separator)
    ENDIF
  END SUBROUTINE write_streams
!------------------------------------------------------------------------------
  SUBROUTINE IO_read_streams

  USE mo_memory_base, ONLY: ostreams, nstreams
  USE mo_filename,    ONLY: standard_rerun_file
  USE mo_control,     ONLY: ngl, nlon
  
    LOGICAL           :: inp_stream  (SIZE(ostreams)) ! streams to read
    LOGICAL           :: in_this_file(SIZE(ostreams)) ! streams to this file
    INTEGER           :: i, j                         ! loop indices
    TYPE(FILE_INFO)   :: fileinfo                     ! 
    INTEGER           :: ios                          ! status from open
    INTEGER           :: dimlen                       ! dimension
    !
    IF (p_pe == p_io .OR. lindependent_read) THEN
      l_read = .TRUE.
    ENDIF
    !
    ! mark all streams to be written
    !
    inp_stream = ostreams% lrerun 
    !
    ! loop over files to process
    !
    DO i = 1,nstreams
      IF (.NOT. inp_stream(i)) CYCLE
      !
      ! identify all streams to be written to this file
      !
      in_this_file = inp_stream .AND. ostreams(i)%rest_suf == ostreams%rest_suf
      !
      ! open file
      !
      fileinfo%opened = .FALSE.
      fileinfo%format = inquire_rerun_file(TRIM(standard_rerun_file) &
                                         //TRIM(ostreams(i)% rest_suf))
      !
      IF (fileinfo%format == NONE .AND. ostreams(i)% lcontnorest) CYCLE
      !
      IF (l_read) THEN
        CALL io_open (filename = TRIM(standard_rerun_file)&
                               //TRIM(ostreams(i)% rest_suf), &
                      fileinfo = fileinfo,                &
                      mode     = IO_READ,                 &
                      iostat   = ios)
        !
        ! check dimensions of the restart files
        !    (model dimensions are read from the _echam restart file in 
        !    restarted runs, but the other restart files do not necessarily 
        !    match)
        !

        CALL IO_INQ_DIMID  (fileinfo%file_id, 'lat', IO_dim_id)
        CALL IO_INQ_DIMLEN (fileinfo%file_id, IO_dim_id, dimlen)
        IF (dimlen /= ngl) THEN
          WRITE (message_text,*) 'restart file dimensions do not match: ngl=', &
               ngl, '; dimlen(lat): ', dimlen
          CALL finish ('IO_read_streams',message_text)
        ENDIF
        CALL IO_INQ_DIMID  (fileinfo%file_id, 'lon', IO_dim_id)
        CALL IO_INQ_DIMLEN (fileinfo%file_id, IO_dim_id, dimlen)
        IF (dimlen /= nlon) THEN
          WRITE (message_text,*) 'restart file dimensions do not match: nlon=', &
               nlon, '; dimlen(lon): ', dimlen
          CALL finish ('IO_read_streams',message_text)
        ENDIF
        SELECT CASE (TRIM(ostreams(i)% rest_suf))
        CASE ('_jsbach','_veg')
          CALL IO_INQ_DIMID  (fileinfo%file_id, 'tiles', IO_dim_id)
          CALL IO_INQ_DIMLEN (fileinfo%file_id, IO_dim_id, dimlen)
          IF (dimlen /= io_dim_ids(TILES)%dim_len) THEN
            WRITE (message_text,*) 'restart file dimensions do not match: ntiles=', &
                 io_dim_ids(TILES)%dim_len, '; dimlen(tiles) in ', &
                 TRIM(standard_rerun_file)//TRIM(ostreams(i)% rest_suf), ': ', dimlen
            CALL finish ('IO_read_streams',message_text)
          ENDIF
        END SELECT
      ENDIF
      !
      ! read variables
      !
      DO j=i,nstreams
        IF (in_this_file(j)) THEN
            CALL io_read_stream (fileinfo, ostreams(j)% first_list_element,ios)
        ENDIF
      END DO
      !
      ! close file
      !
      IF (l_read) &
           CALL io_close (fileinfo)
      !
      ! mark streams written to this file
      !
      inp_stream = inp_stream .AND. .NOT. in_this_file
    END DO

  END SUBROUTINE IO_read_streams
!------------------------------------------------------------------------------
  SUBROUTINE IO_read_stream (fileinfo, io_list, ios)

    USE mo_jsbach_comm_to_echam5mods, ONLY: mask, domain_mask

    TYPE (FILE_INFO), INTENT(IN)    :: fileinfo
    TYPE (list_element), POINTER    :: io_list
    INTEGER,             INTENT(in) :: ios     ! status from call to io_open

    REAL(dp), POINTER :: zin(:,:,:,:), zptr(:,:,:,:), zin2(:,:)
    REAL(dp), ALLOCATABLE :: zptr2(:,:)

    TYPE (list_element), POINTER :: next
    INTEGER :: fileID, status, i, j
    !
    IF (p_pe == p_io .OR. lindependent_read) THEN
      l_read = .TRUE.
    ENDIF
    !
    ! loop over list elements
    !
    fileID = fileinfo%file_id
    next => io_list
    DO WHILE (ASSOCIATED(next))
      IF (next%field%info%lrerun) THEN
        !
        ! get var_id, broadcast status, set restart_read flag on success
        ! on error try uppercase for old restart files
        !
        IF (l_read) THEN
          IF (ios /= 0) THEN
            status = NF_NOERR - 1
          ELSE
            status = NF_INQ_VARID(fileID, next%field%info%name, IO_var_id)
          ENDIF
          IF (status /= NF_NOERR) &
            status = NF_INQ_VARID(fileID, &
              toupper(next%field%info%name), IO_var_id)
          IF (status /= NF_NOERR) &
               CALL message('IO_read_stream','failed to read from rerun file: '&
                            //TRIM(next%field%info%name))
        ENDIF
        IF (.NOT. lindependent_read) THEN
          CALL p_bcast (status, p_io)
        ENDIF
        next%field%info%restart_read = status == NF_NOERR

        IF (status == NF_NOERR) THEN
          !
          ! read field if no error occured
          !
          IF (l_read) THEN
            ALLOCATE(zin(next%field%info%gdim(1), &
                     next%field%info%gdim(2), &
                     next%field%info%gdim(3), &
                     next%field%info%gdim(4)))
            CALL IO_GET_VAR_DOUBLE (fileID, IO_var_id, zin)
          ELSE
            NULLIFY (zin)
          ENDIF ! p_io
          !
          ! scatter the field over the processors
          !
          zptr => next%field%ptr(:,:,:,:)
          SELECT CASE (next%field%info%repr)
!OSBTRSC          CASE (FOURIER)
!OSBTRSC            CALL scatter_fourier (zin, next%field%info%gdim, zptr, lindependent_read)
!OSBTRSC          CASE (GAUSSIAN)
!OSBTRSC            CALL scatter_field (zin, next%field%info%gdim, zptr, lindependent_read)
          CASE (LAND)
            IF (next%field%info%gdim(4)>1) &
                 CALL finish('IO_read_streams','Only 3 dimensions in LAND streams allowed')
            IF (l_read) ALLOCATE(zin2(SIZE(mask,1),SIZE(mask,2)))
            ALLOCATE(zptr2(SIZE(domain_mask,1),SIZE(domain_mask,2)))
            DO j=1,next%field%info%gdim(3)
               DO i=1,next%field%info%gdim(2)
                  IF (l_read) zin2 = UNPACK(zin(:,i,j,1), MASK=mask, FIELD=0._dp)
!OSBTRSC                  CALL scatter_field (zin2, zptr2, lindependent_read)
                  next%field%ptr(:,i,j,1) = PACK(zptr2, MASK=domain_mask)
               END DO
            END DO
            DEALLOCATE(zptr2)
            IF (l_read) DEALLOCATE(zin2)
!OSBTRSC          CASE (SPECTRAL)
!OSBTRSC            CALL scatter_spectral (zin, next%field%info%gdim, zptr, lindependent_read)
          CASE default
            CALL finish('IO_read_stream', &
                        'wrong representation of: '//next%field%info%name)
          END SELECT
          IF(ASSOCIATED(zin)) DEALLOCATE(zin)
        ELSE
          !
          ! Error handling
          !
          IF (.NOT. next% field% info% contnorest) THEN
              CALL finish ('IO_read_stream',                                 &
                           'variable expected in rerun file is not present: '&
                           //TRIM(next%field%info%name))
          ENDIF
        ENDIF ! nf_noerr
      ENDIF   ! restart
      next => next%next_list_element
    END DO     ! next

  END SUBROUTINE IO_read_stream
!------------------------------------------------------------------------------
  SUBROUTINE cleanup_io
    !
    ! deallocate module variables
    !
    USE mo_control,       ONLY: vct

    IF (ALLOCATED  (vlon)) DEALLOCATE (vlon)
    IF (ALLOCATED  (vlat)) DEALLOCATE (vlat)
    IF (ASSOCIATED (vct))  DEALLOCATE (vct)
  END SUBROUTINE cleanup_io
!------------------------------------------------------------------------------
  SUBROUTINE IO_open_unit(unit, fileinfo, mode)

    INTEGER,          INTENT(IN)    :: unit
    INTEGER,          INTENT(IN)    :: mode
    TYPE (FILE_INFO), INTENT(INOUT) :: fileinfo

    CHARACTER(len=13) :: filename

    IF (fileinfo%opened) THEN
      WRITE(message_text,'(a,i0,a,a)') 'IO_open_unit: unit ',unit, &
           ' allready assigned to ', fileinfo%file_name
      CALL message ('',message_text)
      CALL finish  ('IO_open_unit', 'Run terminated.')
    ENDIF
    
    IF (unit < 10 .OR. unit > 99) THEN
      WRITE(message_text,'(a,i0,a)') &
           'IO_open_unit: unit ', unit, ' out of range'
      CALL message ('',message_text)
      CALL finish  ('IO_open_unit', 'Run terminated.')
    ENDIF

    WRITE (filename,'(A5,I2.2)') 'unit.',unit 
    
    CALL IO_open(filename, fileinfo, mode)
    
  END SUBROUTINE IO_open_unit
  !==============================================================================
END MODULE mo_io
