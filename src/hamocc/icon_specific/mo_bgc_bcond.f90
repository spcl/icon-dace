! Allocation/deallocation and reading of HAMOCC boundary conditions
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

MODULE mo_bgc_bcond

USE mo_master_control,       ONLY: get_my_process_name
  USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: filename_max
  USE mo_parallel_config,    ONLY: nproma
  USE mo_impl_constants,     ONLY: max_char_length
  USE mo_model_domain,       ONLY: t_patch, t_patch_3D
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_grid_config,        ONLY: n_dom
  USE mo_grid_subset,        ONLY: get_index_range, t_subset_range
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_hamocc_types,       ONLY: t_hamocc_bcond
  USE mo_var_list,           ONLY: add_var, t_var_list_ptr
  USE mo_var_list_register,  ONLY: vlr_add, vlr_del
  USE mo_master_config,      ONLY: getModelBaseDir
  USE mo_cf_convention,      ONLY: t_cf_var
  USE mo_grib2,              ONLY: t_grib2_var, grib2_var
  USE mo_read_interface,     ONLY: openInputFile, closeFile, on_cells, &
    &                              t_stream_id, read_3D
  USE mo_util_string,        ONLY: t_keyword_list,  &
    &                              associate_keyword
  USE mo_cdi,                ONLY: DATATYPE_FLT32, GRID_UNSTRUCTURED
  USE mo_zaxis_type,         ONLY: ZA_SURFACE
  USE mo_hamocc_nml,         ONLY: io_stdo_bgc
  USE mo_ext_data_types,     ONLY: t_external_data, t_external_bgc
  USE mo_ocean_ext_data,     ONLY: ext_data
  USE mtime,                 ONLY: datetime
  USE mo_cdi_constants,      ONLY: GRID_UNSTRUCTURED_CELL, GRID_CELL
  USE mo_ocean_nml,          ONLY: lsediment_only
  USE mo_run_config,         ONLY: dtime
  USE mo_netcdf_errhandler,  ONLY: nf
  USE mo_fortran_tools,      ONLY: set_acc_host_or_device
  USE mo_netcdf

#include "add_var_acc_macro.inc"

  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_bgc_ext_data
  PUBLIC :: destruct_bgc_ext_data
  PUBLIC :: update_bgc_bcond
  PUBLIC :: ext_data_bgc

  TYPE(t_hamocc_bcond),TARGET :: ext_data_bgc


!-------------------------------------------------------------------------

CONTAINS



  !-------------------------------------------------------------------------
  !>
  !!
!<Optimize:inUse>
  SUBROUTINE construct_bgc_ext_data (p_patch, ext_data,pext_data_bgc)

    TYPE(t_patch), INTENT(IN)            :: p_patch
    TYPE(t_external_data), INTENT(INOUT) :: ext_data
    TYPE(t_hamocc_bcond), INTENT(INOUT) :: pext_data_bgc

    INTEGER :: jg
    CHARACTER(len=MAX_CHAR_LENGTH) :: listname

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_bgc_bcond:construct_bgc_data'

    !-------------------------------------------------------------------------
    CALL message (TRIM(routine), 'Start')


    ! top-level procedure for building data structures for
    ! external data.
    CALL message (TRIM(routine), 'Construction of data structure for ' // &
      &                          'external data started')


    WRITE(listname,'(a,i2.2)') 'ext_data_bgc_D',jg
    ! At this point ext_data is already copied to the device, because of an early copy of ext_data%oce
    !$ACC ENTER DATA COPYIN(ext_data%bgc)
    CALL new_ext_data_bgc_list(p_patch, ext_data%bgc, pext_data_bgc, ext_data%bgc_list, TRIM(listname))


    CALL message (TRIM(routine), 'Construction of data structure for ' // &
      &                          'external data finished')


    IF (lsediment_only) THEN
     CALL read_ext_data_sedon(p_patch, ext_data)
    ELSE
     CALL read_ext_data_bgc(p_patch, ext_data)
    ENDIF

    !$ACC ENTER DATA COPYIN(ext_data_bgc, ext_data_bgc%dusty, ext_data_bgc%nitro)

  END SUBROUTINE construct_bgc_ext_data
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !
  !!
!<Optimize:inUse>
  SUBROUTINE new_ext_data_bgc_list ( p_patch, p_ext_bgc, p_ext_data_bgc, p_ext_bgc_list, &
    &                                listname)
!
    TYPE(t_patch), TARGET, INTENT(IN)   :: & !< current patch
      &  p_patch

    TYPE(t_external_bgc), INTENT(INOUT) :: & !< current external data structure
      &  p_ext_bgc
    TYPE(t_hamocc_bcond), INTENT(INOUT) :: & !< current external data structure
      &  p_ext_data_bgc

    TYPE(t_var_list_ptr) :: p_ext_bgc_list !< current external data list

    CHARACTER(len=*), INTENT(IN)  :: & !< list name
      &  listname

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c    !< number of cell blocks to allocate

    INTEGER :: shape2d_c(2),  shape3d_c(3)

    INTEGER :: ibits         !< "entropy" of horizontal slice


    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%alloc_cell_blocks


    ibits = 16   ! "entropy" of horizontal slice

    ! predefined array shapes
    shape2d_c = (/ nproma, nblks_c /)


    shape3d_c = (/ nproma, 12, nblks_c/)

    !
    ! Register a field list and apply default settings
    !
    CALL vlr_add(p_ext_bgc_list, TRIM(listname), patch_id=p_patch%id, &
      &               lrestart=.FALSE., model_type=TRIM(get_my_process_name()))
    IF (.not. lsediment_only) THEN
     cf_desc    = t_cf_var('Dust cell center', 'kg m-2 yr-1', &
       &                   'DUST', DATATYPE_FLT32)
     grib2_desc = grib2_var( 192, 140, 219, ibits, GRID_UNSTRUCTURED, GRID_CELL)
     CALL add_var( p_ext_bgc_list, 'DUSTin', p_ext_bgc%dust,      &
       &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3d_c, lopenacc = .TRUE. )
     __acc_attach(p_ext_bgc%dust)
     CALL add_var( p_ext_bgc_list, 'DUSTY', p_ext_data_bgc%dusty,      &
       &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )
     cf_desc    = t_cf_var('Nitrogen cell center', 'kg m-2 yr-1', &
       &                   'NDEP', DATATYPE_FLT32)
     grib2_desc = grib2_var( 192, 140, 239, ibits, GRID_UNSTRUCTURED, GRID_CELL)
     CALL add_var( p_ext_bgc_list, 'NDEP', p_ext_bgc%nitro,      &
       &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3d_c, lopenacc = .TRUE. )
     __acc_attach(p_ext_bgc%nitro)
     CALL add_var( p_ext_bgc_list, 'NITRO', p_ext_data_bgc%nitro,      &
       &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )
    ELSE
     cf_desc    = t_cf_var('dust sediment flx', 'kmol m-2 s', &
       &                   'produs', DATATYPE_FLT32)
     grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
     CALL add_var( p_ext_bgc_list, 'PRO_DUS', p_ext_bgc%produs,      &
       &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )
     CALL add_var( p_ext_bgc_list, 'produs', p_ext_data_bgc%produs,      &
       &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )
     cf_desc    = t_cf_var('organic sediment flx', 'kmol m-2 s', &
       &                   'prorca', DATATYPE_FLT32)
     CALL add_var( p_ext_bgc_list, 'pr_orca', p_ext_bgc%prorca,      &
       &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )
     CALL add_var( p_ext_bgc_list, 'prorca', p_ext_data_bgc%prorca,      &
       &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )
     cf_desc    = t_cf_var('calc sediment flx', 'kmol m-2 s', &
       &                   'prcaca', DATATYPE_FLT32)
     CALL add_var( p_ext_bgc_list, 'pr_caca', p_ext_bgc%prcaca,      &
       &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )
     CALL add_var( p_ext_bgc_list, 'prcaca', p_ext_data_bgc%prcaca,      &
       &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )
     cf_desc    = t_cf_var('opal sediment flx', 'kmol m-2 s', &
       &                   'silpro', DATATYPE_FLT32)
     CALL add_var( p_ext_bgc_list, 'sil_pro', p_ext_bgc%silpro,      &
       &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )
     CALL add_var( p_ext_bgc_list, 'silpro', p_ext_data_bgc%silpro,      &
       &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )
     CALL message('new_ext_data_bgc_list','add_var finished for particle fluxes')
    ENDIF

  END SUBROUTINE new_ext_data_bgc_list
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Destruct external data data structure and lists
  !!
!<Optimize:inUse>
  SUBROUTINE destruct_bgc_ext_data

    INTEGER :: jg
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_bgc_bcond:destruct_bgc_ext_data'
    !-------------------------------------------------------------------------

    CALL message (TRIM(routine), 'Destruction of data structure for ' // &
      &                          'external data started')

    DO jg = 1,n_dom
      ! Delete list of constant in time oceanic elements
      CALL vlr_del(ext_data(jg)%bgc_list)
    ENDDO

    CALL message (TRIM(routine), 'Destruction of data structure for ' // &
      &                          'external data finished')

  END SUBROUTINE destruct_bgc_ext_data
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !!
  !! Read HAMOCC external data from netcdf
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE read_ext_data_bgc (p_patch, ext_data)

    TYPE(t_patch), INTENT(IN)            :: p_patch
    TYPE(t_external_data), INTENT(INOUT) :: ext_data

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'read_ext_data_bgc'

    CHARACTER(filename_max) :: dust_file   !< file name for reading in

    LOGICAL :: l_exist
    INTEGER ::  no_cells,  no_tst
    INTEGER :: ncid, dimid
    TYPE(t_stream_id) :: stream_id

    REAL(wp):: z_flux(nproma,12,p_patch%alloc_cell_blocks)
    TYPE (t_keyword_list), POINTER :: keywords => NULL()

    CALL message (TRIM(routine), 'start')

    !-------------------------------------------------------------------------
    !  READ DUST
    !-------------------------------------------------------------------------

!     jg = 1

    z_flux(:,:,:) = 0.0_wp

    CALL associate_keyword("<path>", TRIM(getModelBaseDir()), keywords)

      dust_file='dust.nc'

      CALL message( TRIM(routine),'HAMOCC dust file is: '//TRIM(dust_file) )

      IF(my_process_is_stdio()) THEN
        !
        INQUIRE (FILE=dust_file, EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          write(io_stdo_bgc,*)'FORCING FILE: ',TRIM(dust_file)
          CALL finish(TRIM(routine),'DUST file is not found - ABORT')
        ENDIF

        !
        ! open file
        !
        CALL nf(nf_open(TRIM(dust_file), NF_NOWRITE, ncid), routine)
        CALL message( TRIM(routine),'HAMOCC dust file opened for read' )

        !
        !
        CALL nf(nf_inq_dimid (ncid, 'ncells', dimid), routine)
        CALL nf(nf_inq_dimlen(ncid, dimid, no_cells), routine)

        IF(p_patch%n_patch_cells_g /= no_cells) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of patch cells and cells in HAMOCC dust file do not match - ABORT')
        ENDIF

        !
        ! get number of timesteps
        !
        CALL nf(nf_inq_dimid (ncid, 'time', dimid), routine)
        CALL nf(nf_inq_dimlen(ncid, dimid, no_tst), routine)
        !
        ! check
        !
        WRITE(message_text,'(A,I6,A)')  'HAMOCC dust file contains',no_tst,' data sets'
        CALL message( TRIM(routine), TRIM(message_text) )
        IF(no_tst /= 12 ) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of forcing timesteps is not equal 12 specified in namelist - ABORT')
        ENDIF

        CALL nf(nf_close(ncid), routine)
      ENDIF

      CALL openinputfile(stream_id, dust_file, p_patch)

      no_tst = 12
      !-------------------------------------------------------
      !
      ! Read dust for triangle centers
      !
      !-------------------------------------------------------

        CALL read_3D(stream_id, on_cells, 'DUST', z_flux)
        ext_data%bgc%dust(:,:,:) = z_flux(:,:,:)


      !
      ! close file
      !
      CALL closeFile(stream_id)


      CALL message( TRIM(routine),'HAMOCC dust file read' )


      dust_file='nitrogen.nc'

      CALL message( TRIM(routine),'HAMOCC nitrogen input file is: '//TRIM(dust_file) )

      IF(my_process_is_stdio()) THEN
        !
        INQUIRE (FILE=dust_file, EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          write(io_stdo_bgc,*)'FORCING FILE: ',TRIM(dust_file)
          CALL finish(TRIM(routine),'Nitrogen input file is not found - ABORT')
        ENDIF

        !
        ! open file
        !
        CALL nf(nf_open(TRIM(dust_file), NF_NOWRITE, ncid), routine)
        CALL message( TRIM(routine),'HAMOCC nitrogen input file opened for read' )

        !
        !
        CALL nf(nf_inq_dimid (ncid, 'ncells', dimid), routine)
        CALL nf(nf_inq_dimlen(ncid, dimid, no_cells), routine)

        IF(p_patch%n_patch_cells_g /= no_cells) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of patch cells and cells in HAMOCC input file do not match - ABORT')
        ENDIF

        !
        ! get number of timesteps
        !
        CALL nf(nf_inq_dimid (ncid, 'time', dimid), routine)
        CALL nf(nf_inq_dimlen(ncid, dimid, no_tst), routine)
        !
        ! check
        !
        WRITE(message_text,'(A,I6,A)')  'HAMOCC nitrogen input file contains',no_tst,' data sets'
        CALL message( TRIM(routine), TRIM(message_text) )
        IF(no_tst /= 12 ) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of forcing timesteps is not equal 12 specified in namelist - ABORT')
        ENDIF

        CALL nf(nf_close(ncid), routine)
      ENDIF

      CALL openinputfile(stream_id, dust_file, p_patch)

      no_tst = 12
      !-------------------------------------------------------
      !
      ! Read nitrogen for triangle centers
      !
      !-------------------------------------------------------
        CALL read_3D(stream_id, on_cells, 'ndepo', z_flux)
        ext_data%bgc%nitro(:,:,:) = z_flux(:,:,:)


      !
      ! close file
      !
      CALL closeFile(stream_id)




  END SUBROUTINE read_ext_data_bgc
  !--------------------------------------------------
  !>
  !!
  !! Read HAMOCC external data from netcdf
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE read_ext_data_sedon (p_patch, ext_data)

    TYPE(t_patch), INTENT(IN)            :: p_patch
    TYPE(t_external_data), INTENT(INOUT) :: ext_data

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'read_ext_data_bgc'

    CHARACTER(filename_max) :: dust_file   !< file name for reading in

    LOGICAL :: l_exist
    INTEGER ::  no_cells,  no_tst
    INTEGER :: ncid, dimid
    TYPE(t_stream_id) :: stream_id

    REAL(wp):: z_flux(nproma,1,p_patch%alloc_cell_blocks)
    TYPE (t_keyword_list), POINTER :: keywords => NULL()

    CALL message (TRIM(routine), 'start')

    !-------------------------------------------------------------------------
    !  READ DUST
    !-------------------------------------------------------------------------

!     jg = 1

    z_flux(:,:,:) = 0.0_wp

    CALL associate_keyword("<path>", TRIM(getModelBaseDir()), keywords)

      dust_file='particle_fluxes.nc'

      CALL message( TRIM(routine),'HAMOCC sediment flux file is: '//TRIM(dust_file) )

      IF(my_process_is_stdio()) THEN
        !
        INQUIRE (FILE=dust_file, EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          write(io_stdo_bgc,*)'FORCING FILE: ',TRIM(dust_file)
          CALL finish(TRIM(routine),'Particle flux file is not found - ABORT')
        ENDIF

        !
        ! open file
        !
        CALL nf(nf_open(TRIM(dust_file), NF_NOWRITE, ncid), routine)
        CALL message( TRIM(routine),'HAMOCC particle flux file opened for read' )

        !
        !
        CALL nf(nf_inq_dimid (ncid, 'ncells', dimid), routine)
        CALL nf(nf_inq_dimlen(ncid, dimid, no_cells), routine)

        IF(p_patch%n_patch_cells_g /= no_cells) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of patch cells and cells in HAMOCC particle_flux file do not match - ABORT')
        ENDIF

        !
        ! get number of timesteps
        !
        CALL nf(nf_inq_dimid (ncid, 'time', dimid), routine)
        CALL nf(nf_inq_dimlen(ncid, dimid, no_tst), routine)
        !
        ! check
        !
        WRITE(message_text,'(A,I6,A)')  'HAMOCC particle flux file contains',no_tst,' data sets'
        CALL message( TRIM(routine), TRIM(message_text) )

        CALL nf(nf_close(ncid), routine)
      ENDIF

      CALL openInputFile(stream_id, dust_file, p_patch)

      no_tst = 12
      !-------------------------------------------------------
      !
      ! Read sediment fluxes for triangle centers
      !
      !-------------------------------------------------------

        CALL read_3D(stream_id, on_cells, 'prorca', z_flux)
        ext_data%bgc%prorca(:,:) = z_flux(:,1,:)*dtime
        CALL read_3D(stream_id, on_cells, 'prcaca', z_flux)
        ext_data%bgc%prcaca(:,:) = z_flux(:,1,:)*dtime
        CALL read_3D(stream_id, on_cells, 'silpro', z_flux)
        ext_data%bgc%silpro(:,:) = z_flux(:,1,:)*dtime
        CALL read_3D(stream_id, on_cells, 'produs', z_flux)
        ext_data%bgc%produs(:,:) = z_flux(:,1,:)*dtime

      !
      ! close file
      !
      CALL closeFile(stream_id)


      CALL message( TRIM(routine),'HAMOCC particle-flux file read' )

  END SUBROUTINE read_ext_data_sedon
  !--------------------------------------------------

!<Optimize:inUse>

   SUBROUTINE update_bgc_bcond(p_patch_3D, bgc_ext, this_datetime, lacc)
    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_hamocc_bcond)                        :: bgc_ext
    TYPE(datetime), INTENT(IN)                  :: this_datetime
    LOGICAL, INTENT(IN), OPTIONAL               :: lacc


 ! local variables
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_bgc_bcond:update_bgc_bcond'
    INTEGER  :: jmon, jdmon, jmon1, jmon2
    REAL(wp) :: rday1, rday2
    INTEGER  :: jc, jb, startidx, endidx
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: p_patch
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)


    !  calculate day and month
    jmon  = this_datetime%date%month         ! integer current month
    jdmon = this_datetime%date%day           ! integer day in month


      jmon1=jmon-1
      jmon2=jmon
      rday1=REAL(15-jdmon,wp)/30.0_wp
      rday2=REAL(15+jdmon,wp)/30.0_wp
      IF (jdmon > 15)  THEN
        jmon1=jmon
        jmon2=jmon+1
        rday1=REAL(45-jdmon,wp)/30.0_wp
        rday2=REAL(jdmon-15,wp)/30.0_wp
      END IF

      IF (jmon1 ==  0) jmon1=12
      IF (jmon1 == 13) jmon1=1
      IF (jmon2 ==  0) jmon2=12
      IF (jmon2 == 13) jmon2=1

      IF (lsediment_only) THEN
       bgc_ext%prorca(:,:) = ext_data(1)%bgc%prorca(:,:)
       bgc_ext%prcaca(:,:) = ext_data(1)%bgc%prcaca(:,:)
       bgc_ext%produs(:,:) = ext_data(1)%bgc%produs(:,:)
       bgc_ext%silpro(:,:) = ext_data(1)%bgc%silpro(:,:)
      ELSE
       p_patch => p_patch_3D%p_patch_2D(1)
       all_cells => p_patch%cells%all
       DO jb = all_cells%start_block, all_cells%end_block
         call get_index_range(all_cells, jb, startidx, endidx)
         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
         !$ACC LOOP GANG VECTOR
         DO jc = startidx, endidx
           bgc_ext%dusty(jc,jb) = rday1*ext_data(1)%bgc%dust(jc,jmon1,jb) + &
           &                                   rday2*ext_data(1)%bgc%dust(jc,jmon2,jb)
           bgc_ext%nitro(jc,jb) = rday1*ext_data(1)%bgc%nitro(jc,jmon1,jb) + &
           &                                   rday2*ext_data(1)%bgc%nitro(jc,jmon2,jb)
         END DO
         !$ACC END PARALLEL
       END DO
      ENDIF

  END SUBROUTINE update_bgc_bcond
END MODULE

