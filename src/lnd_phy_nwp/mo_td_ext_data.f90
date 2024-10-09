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

! Subroutines needed for updating the climatological boundary conditions
! in NWP physics.
!
! includes, e.g. updating time dependent SST and Sea Ice fraction

MODULE mo_td_ext_data

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_io_units,            ONLY: filename_max
  USE mo_master_config,       ONLY: getModelBaseDir
  USE mo_io_config,           ONLY: default_read_method
  USE mo_grid_config,         ONLY: n_dom
  USE mo_extpar_config,       ONLY: generate_td_filename, itopo, itype_vegetation_cycle, itype_lwemiss
  USE mo_lnd_nwp_config,      ONLY: sst_td_filename, ci_td_filename, sstice_mode, frsi_min
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_radiation_config,    ONLY: albedo_type
  USE mo_run_config,          ONLY: msg_level
  USE mo_dynamics_config,     ONLY: nnow_rcf, nnew_rcf
  USE mo_read_interface,      ONLY: openInputFile, closeFile, on_cells, &
    &                               t_stream_id, read_2D_1time
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_nonhydro_types,      ONLY: t_nh_state
  USE mo_nwp_lnd_types,       ONLY: t_lnd_state
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, min_rlcell_int, SSTICE_CLIM, &
    &                               SSTICE_AVG_MONTHLY, SSTICE_AVG_DAILY, SSTICE_ANA_CLINC, &
    &                               MODIS, ivdiff
  USE mo_ext_data_init,       ONLY: interpol_monthly_mean
  USE mo_nwp_sfc_utils,       ONLY: process_sst_and_seaice, sst_add_climatological_incr, &
    &                               update_ndvi_dependent_fields
  USE mo_nwp_vdiff_interface, ONLY: nwp_vdiff_update_seaice
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mtime,                  ONLY: datetime, newDatetime, deallocateDatetime, &
    &                               datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights,         &
    &                                  calculate_time_interpolation_weights
  USE mo_dynamics_config,     ONLY: nnow_rcf
  USE mo_nwp_phy_init,        ONLY: clim_cdnc

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: update_nwp_phy_bcs
  PUBLIC  :: set_sst_and_seaice
  PUBLIC  :: read_td_ext_data_file

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_td_extdata'

CONTAINS

  !>
  !! This subroutine updates the lower boundary conditions for NWP physics
  !!
  !! This subroutine updates the lower boundary conditions for NWP physics.
  !! The following fields are updated:
  !! sea surface temperature: t_seasfc
  !! seaice fraction        : fr_seaice
  !! normalized difference vegetation index ratio: ndviratio
  !! modis albedo           : alb_dif, albuv_dif, albni_dif
  !! climatological cloud droplet number         : cdnc
  !!
  !! Update of aerosol fields, which is based on a monthly climatology 
  !! happens in mo_nwp_rrtm_interface:nwp_aerosol.
  !!
  !! In addition, some fields which directly depend on the above ones
  !! are updated as well, in order to have a consistent state.
  !! Prominent examples are the sea index lists, t_g, t_s, plcov, sai, tai.
  !!
  SUBROUTINE update_nwp_phy_bcs (p_patch, ext_data, p_lnd_state, p_nh_state, prm_diag, &
    &                                ref_datetime, target_datetime, mtime_old )

    TYPE(t_patch)         , INTENT(IN)    :: p_patch
    TYPE(t_external_data) , INTENT(INOUT) :: ext_data
    TYPE(t_lnd_state)     , INTENT(INOUT) :: p_lnd_state
    TYPE(t_nh_state)      , INTENT(IN)    :: p_nh_state
    TYPE(t_nwp_phy_diag)  , INTENT(INOUT) :: prm_diag        
    TYPE(datetime)        , INTENT(IN)    :: ref_datetime     !< reference date for computation 
                                                              ! of climatological SST increments
    TYPE(datetime)        , INTENT(IN)    :: target_datetime  !< target date for temporal interpolation
    TYPE(datetime)        , INTENT(IN)    :: mtime_old

    ! local
    INTEGER :: jg
    INTEGER :: ierr
    CHARACTER(len=MAX_DATETIME_STR_LEN) :: target_datetime_PTString
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':update_nwp_phy_bcs'
  !--------------------------------------------------------------------------

    jg = p_patch%id

    IF (msg_level >= 12) THEN
      CALL datetimeToString(target_datetime, target_datetime_PTString, ierr)
      WRITE(message_text,'(a,i3,a,a)') 'Update clim. LBCs on DOM ',jg,' for target date ', &
        &                                 TRIM(target_datetime_PTString)
      CALL message(TRIM(routine),message_text)
    ENDIF

    ! Update ndvi normalized differential vegetation index
    IF (itopo == 1 .AND. atm_phy_nwp_config(jg)%inwp_surface >= 1) THEN
      CALL interpol_monthly_mean(p_patch,                     &! in
        &                        target_datetime,             &! in
        &                        ext_data%atm_td%ndvi_mrat,   &! in
        &                        ext_data%atm%ndviratio       )! out
      IF (itype_vegetation_cycle > 1) THEN
        CALL interpol_monthly_mean(p_patch,                   &! in
          &                        target_datetime,           &! in
          &                        ext_data%atm_td%t2m_m,     &! in
          &                        ext_data%atm%t2m_clim,     &! out
          &                        ext_data%atm%t2m_climgrad  )! optional out
      ENDIF

      ! updates plcov, tai and sai, which depend on ndvi
      CALL update_ndvi_dependent_fields(p_patch, ext_data, p_nh_state%diag)
    END IF


    SELECT CASE(sstice_mode)
    CASE (SSTICE_ANA_CLINC)
      !
      ! add climatological SST increment (updates t_g_t, t_s_t, t_sk_t, qv_s_t)
      !
      CALL sst_add_climatological_incr (p_patch        = p_patch,                       &
        &                               ext_data       = ext_data,                      &
        &                               prog_lnd       = p_lnd_state%prog_lnd(:),       &
        &                               diag_lnd       = p_lnd_state%diag_lnd,          &
        &                               pres_sfc       = p_nh_state%diag%pres_sfc(:,:), &
        &                               ref_datetime   = ref_datetime,                  &
        &                               target_datetime= target_datetime )

      ! VDIFF updates its SST internally. The routine above cannot be used because there
      ! are no ocean tiles in `p_lnd_state` with JSBACH/VDIFF.

    CASE(SSTICE_CLIM, SSTICE_AVG_MONTHLY, SSTICE_AVG_DAILY)
      !
      ! update SST (t_seasfc) and sea ice fraction (fr_seaice)
      !
      CALL set_sst_and_seaice (.FALSE., target_datetime,                   &
                                &  mtime_old, sstice_mode,                 &
                                &  p_patch, ext_data, p_lnd_state)

      IF (atm_phy_nwp_config(jg)%inwp_turb == ivdiff) THEN
        ! lacc fixed to TRUE here because current sea-ice variables are on the GPU.
        ! `set_sst_and_seaice` pushes the new SST and sea-ice fraction to the GPU after
        ! interpolating.
        CALL nwp_vdiff_update_seaice ( &
            & p_patch, .TRUE., p_lnd_state%diag_lnd%fr_seaice(:,:), &
            & ext_data%atm%list_sea, ext_data%atm%list_seaice, p_lnd_state%prog_wtr(nnow_rcf(jg)), &
            & lacc=.TRUE. &
          )
      ELSE
        ! rebuild index lists for water and seaice based on fr_seaice,
        ! and update tiled surface temperatures
        !
        CALL process_sst_and_seaice (p_patch      = p_patch,                             & !in
          &                          fr_seaice    = p_lnd_state%diag_lnd%fr_seaice(:,:), & !in(out)
          &                          t_seasfc     = p_lnd_state%diag_lnd%t_seasfc(:,:),  & !in
          &                          pres_sfc     = p_nh_state%diag%pres_sfc(:,:),       & !in
          &                          ext_data     = ext_data,                            & !inout
          &                          prog_lnd_now = p_lnd_state%prog_lnd(nnow_rcf(jg)),  & !inout
          &                          prog_lnd_new = p_lnd_state%prog_lnd(nnew_rcf(jg)),  & !inout
          &                          prog_wtr_now = p_lnd_state%prog_wtr(nnow_rcf(jg)),  & !inout
          &                          prog_wtr_new = p_lnd_state%prog_wtr(nnew_rcf(jg)),  & !inout
          &                          diag_lnd     = p_lnd_state%diag_lnd  )                !inout
      END IF

    CASE DEFAULT
      !
      ! do nothing (sstice_mode = SSTICE_ana, SSTICE_INST) 
    END SELECT


    ! Check if MODIS albedo needs to be updated
    IF (albedo_type == MODIS) THEN
      ! Note that here only an update of the external parameter fields is
      ! performed. The actual update happens in mo_albedo.
      !
      CALL interpol_monthly_mean(p_patch,                    &! in
        &                        target_datetime,            &! in
        &                        ext_data%atm_td%alb_dif,    &! in
        &                        ext_data%atm%alb_dif        )! out

      CALL interpol_monthly_mean(p_patch,                    &! in
        &                        target_datetime,            &! in
        &                        ext_data%atm_td%albuv_dif,  &! in
        &                        ext_data%atm%albuv_dif      )! out

      CALL interpol_monthly_mean(p_patch,                    &! in
        &                        target_datetime,            &! in
        &                        ext_data%atm_td%albni_dif,  &! in
        &                        ext_data%atm%albni_dif      )! out

    ENDIF

    ! Interpolate also longwave emissivity data if monthly climatology is available
    IF (itype_lwemiss == 2) THEN

      CALL interpol_monthly_mean(p_patch,                   &! in
        &                        target_datetime,           &! in
        &                        ext_data%atm_td%lw_emiss,  &! in
        &                        ext_data%atm%emis_rad      )! out

    ENDIF
    ! Interpolate monthly cdnc climatology
    IF (atm_phy_nwp_config(jg)%icpl_aero_gscp == 3) THEN
      CALL clim_cdnc(target_datetime, p_patch, ext_data, prm_diag)
    ENDIF

  END SUBROUTINE update_nwp_phy_bcs


  !>
  !! Set the SST and sea ice fraction to the values corresponding to datetime
  !!
  !! The SST and sea ice fraction are set to the values of the corresponding
  !!  day and month. Values are interpolated from  climatological monthly
  !!  means (ext_data_mode = 2) or from the actual monthly means (ext_data_mode = 3).
  !!  Another option not yet implemented (ext_data_mode = 4) will set SST and
  !!  sea ice cover to the actual (day, month, year) daily mean
  !!
  SUBROUTINE set_sst_and_seaice (lread, target_datetime, datetime_old, ext_data_mode,  &
                                  &  p_patch, ext_data, p_lnd_state)

    LOGICAL               , INTENT(IN)    :: lread !force the read of the ext dat files for sstice_mode=3
    TYPE(datetime), TARGET, INTENT(IN)    :: target_datetime
    TYPE(datetime), TARGET, INTENT(IN)    :: datetime_old
    INTEGER               , INTENT(IN)    :: ext_data_mode
    TYPE(t_patch)         , INTENT(IN)    :: p_patch
    TYPE(t_external_data) , INTENT(INOUT) :: ext_data
    TYPE(t_lnd_state)     , INTENT(INOUT) :: p_lnd_state

    CHARACTER(len=*), PARAMETER :: &
      routine = modname//':set_actual_td_ext_data  '

    INTEGER                            :: month1, month2, year1, year2
    INTEGER                            :: m1, m2
    REAL (wp)                          :: pw1, pw2, pw2_old
    INTEGER                            :: jb, jc, i_startidx, i_endidx
    INTEGER                            :: i_rlstart, i_rlend
    INTEGER                            :: i_startblk, i_endblk

    TYPE(datetime), POINTER            :: mtime_date
    TYPE(datetime), POINTER            :: mtime_date_old
    TYPE(datetime), POINTER            :: mtime_hour
    
    TYPE(t_time_interpolation_weights) :: current_time_interpolation_weights
    TYPE(t_time_interpolation_weights) :: old_time_interpolation_weights    

    ! this is because newDatetime requires a pointer as input argument
    mtime_date => target_datetime   
    mtime_date_old => datetime_old

    SELECT CASE (ext_data_mode)
      
    CASE (SSTICE_CLIM) 
      ! SST and sea ice fraction updated based on the climatological
      ! monthly values
      
      mtime_hour => newDatetime(mtime_date)
      mtime_hour%time%minute = 0
      mtime_hour%time%second = 0
      mtime_hour%time%ms     = 0                  
      current_time_interpolation_weights = calculate_time_interpolation_weights(mtime_hour)
      CALL deallocateDatetime(mtime_hour)

      month1 = current_time_interpolation_weights%month1
      month2 = current_time_interpolation_weights%month2
      pw2    = current_time_interpolation_weights%weight2
      pw1    = 1._wp - pw2


!$OMP PARALLEL PRIVATE(i_startblk,i_endblk, i_rlstart, i_rlend)
      i_rlstart = grf_bdywidth_c+1
      i_rlend   = min_rlcell_int

      i_startblk = p_patch%cells%start_block(i_rlstart)
      i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
      DO jb=i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, i_rlstart, i_rlend)

        DO jc = i_startidx, i_endidx
          IF (ext_data%atm_td%sst_m(jc,jb,1) >= 0._wp .AND.   &
            & ext_data%atm_td%sst_m(jc,jb,2) >= 0._wp ) THEN
            p_lnd_state%diag_lnd%t_seasfc(jc,jb) =                   &
              pw1*ext_data%atm_td%sst_m(jc,jb,month1)    &
              &         + pw2*ext_data%atm_td%sst_m(jc,jb,month2)
            p_lnd_state%diag_lnd%fr_seaice(jc,jb) =                  &
              pw1*ext_data%atm_td%fr_ice_m(jc,jb,month1) &
              &         + pw2*ext_data%atm_td%fr_ice_m(jc,jb,month2)

            IF (p_lnd_state%diag_lnd%fr_seaice(jc,jb) < frsi_min ) THEN
              p_lnd_state%diag_lnd%fr_seaice(jc,jb)= 0._wp
            ELSE IF (p_lnd_state%diag_lnd%fr_seaice(jc,jb) > 1._wp-frsi_min ) THEN
              p_lnd_state%diag_lnd%fr_seaice(jc,jb)= 1._wp
            END IF
          END IF
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL



    CASE (SSTICE_AVG_MONTHLY) 
      ! SST and sea ice fraction updated based on the actual monthly
      ! values

      mtime_hour => newDatetime(mtime_date_old)
      mtime_hour%time%minute = 0
      mtime_hour%time%second = 0
      mtime_hour%time%ms     = 0                  
      old_time_interpolation_weights = calculate_time_interpolation_weights(mtime_hour)
      call deallocateDatetime(mtime_hour)
      m1 = old_time_interpolation_weights%month1
      m2 = old_time_interpolation_weights%month2

      pw2_old = old_time_interpolation_weights%weight2

      mtime_hour => newDatetime(mtime_date)
      mtime_hour%time%minute = 0
      mtime_hour%time%second = 0
      mtime_hour%time%ms     = 0                  
      current_time_interpolation_weights = calculate_time_interpolation_weights(mtime_hour)
      call deallocateDatetime(mtime_hour)
      month1 = current_time_interpolation_weights%month1
      month2 = current_time_interpolation_weights%month2
      year1  = current_time_interpolation_weights%year1
      year2  = current_time_interpolation_weights%year2
      pw2    = current_time_interpolation_weights%weight2

      IF (m1 /= month1 .OR. lread ) THEN
        CALL read_td_ext_data_file (month1,month2,year1,year2,p_patch,ext_data)
      END IF
      pw1 = 1._wp - pw2



!$OMP PARALLEL PRIVATE(i_startblk,i_endblk, i_rlstart, i_rlend)
      i_rlstart = grf_bdywidth_c+1
      i_rlend   = min_rlcell_int

      i_startblk = p_patch%cells%start_block(i_rlstart)
      i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
      DO jb=i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, i_rlstart, i_rlend)

        DO jc = i_startidx, i_endidx
          IF (ext_data%atm_td%sst_m(jc,jb,1) < 0._wp .AND.   &
            & ext_data%atm_td%sst_m(jc,jb,2) >= 0._wp )  THEN
            WRITE( message_text,'(a,2g10.5)') 'something is wrong,', &
              & ext_data%atm_td%sst_m(jc,jb,1), ext_data%atm_td%sst_m(jc,jb,2)
            CALL message  (routine, TRIM(message_text))
          ELSE IF (ext_data%atm_td%sst_m(jc,jb,1) >= 0._wp .AND.   &
            & ext_data%atm_td%sst_m(jc,jb,2) >= 0._wp ) THEN
            p_lnd_state%diag_lnd%t_seasfc(jc,jb) =              &
              pw1*ext_data%atm_td%sst_m(jc,jb,1)    &
              &         + pw2*ext_data%atm_td%sst_m(jc,jb,2)
            p_lnd_state%diag_lnd%fr_seaice(jc,jb) =            &
              pw1*ext_data%atm_td%fr_ice_m(jc,jb,1) &
              &         + pw2*ext_data%atm_td%fr_ice_m(jc,jb,2)

            IF (p_lnd_state%diag_lnd%fr_seaice(jc,jb) < frsi_min ) THEN
              p_lnd_state%diag_lnd%fr_seaice(jc,jb)= 0._wp
            ELSE IF (p_lnd_state%diag_lnd%fr_seaice(jc,jb) > 1._wp-frsi_min ) THEN
              p_lnd_state%diag_lnd%fr_seaice(jc,jb)= 1._wp
            END IF
          END IF
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL



    CASE (SSTICE_AVG_DAILY) !SST and sea ice fraction updated based
      !  on the actual daily values
      !Not implemented
      WRITE( message_text,'(a)') 'ext_data_mode == 5 not yet implemented '
      CALL finish  (routine, TRIM(message_text))

    CASE DEFAULT
      ! do nothing
    END SELECT

    !$ACC UPDATE DEVICE(p_lnd_state%diag_lnd%t_seasfc, p_lnd_state%diag_lnd%fr_seaice) ASYNC(1)

  END SUBROUTINE set_sst_and_seaice
!-----------------------------------------------------------------------

  SUBROUTINE read_td_ext_data_file (m1,m2,y1,y2,p_patch,ext_data)

    INTEGER, INTENT(IN)                  :: m1,m2,y1,y2 ! month and year of the
                                               ! nearest months to datetime
    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch
    TYPE(t_external_data), INTENT(INOUT) :: ext_data

    CHARACTER(LEN=filename_max) :: extpar_file
    TYPE(t_stream_id)           :: stream_id
    INTEGER :: tlen
    CHARACTER(len=*), PARAMETER :: &
    &  routine = modname//':read_td_ext_data_file:'

!-----------------------------------------------------------------------
    ! extpar_td_filename = "<path>extpar_<year>_<month>_<gridfile>"

    ! Set the months needed to interpolate the ext_para to the actual day

    IF (p_patch%geometry_info%cell_type == 6) THEN ! hexagonal grid
      CALL finish(routine,&
        & 'Hexagonal grid is not supported, yet.')
    ENDIF
    tlen = LEN_TRIM(p_patch%grid_filename)
    WRITE (message_text,'(a,2i6,a)') 'n_dom, jg, grid_filename ',n_dom, &
      &                              p_patch%id, p_patch%grid_filename(1:tlen)
    CALL message(routine, message_text)

    !! READ SST files
    extpar_file = generate_td_filename(sst_td_filename,                   &
      &                                getModelBaseDir(),                 &
      &                                p_patch%grid_filename(1:tlen),   &
      &                                m1,y1                   )

    CALL message(routine, extpar_file)
    call openInputFile(stream_id, extpar_file, p_patch, default_read_method)
    CALL read_2D_1time(stream_id, on_cells, 'SST', &
      &          ext_data%atm_td%sst_m(:,:,1))
    CALL closeFile(stream_id)

    extpar_file = generate_td_filename(sst_td_filename,                   &
      &                                getModelBaseDir(),                 &
      &                                p_patch%grid_filename(1:tlen),   &
      &                                m2, y2                   )
    CALL message(routine, extpar_file)
    call openInputFile(stream_id, extpar_file, p_patch, default_read_method)
    CALL read_2D_1time(stream_id, on_cells, 'SST', &
      &          ext_data%atm_td%sst_m(:,:,2))
    CALL closeFile(stream_id)

    !! READ CI files

    extpar_file = generate_td_filename(ci_td_filename,                    &
      &                                getModelBaseDir(),                 &
      &                                p_patch%grid_filename(1:tlen),   &
      &                                m1,y1                   )
    CALL message(routine, extpar_file)
    call openInputFile(stream_id, extpar_file, p_patch, default_read_method)
    CALL read_2D_1time(stream_id, on_cells, 'CI', &
      &          ext_data%atm_td%fr_ice_m(:,:,1))
    CALL closeFile(stream_id)

    extpar_file = generate_td_filename(ci_td_filename,                    &
      &                                getModelBaseDir(),                 &
      &                             p_patch%grid_filename(1:tlen),      &
      &                             m2,y2                   )
    CALL message(routine, extpar_file)
    call openInputFile(stream_id, extpar_file, p_patch, default_read_method)
    CALL read_2D_1time(stream_id, on_cells, 'CI', &
      &          ext_data%atm_td%fr_ice_m(:,:,2))
    CALL closeFile(stream_id)

  END SUBROUTINE read_td_ext_data_file

!-----------------------------------------------------------------------

END MODULE mo_td_ext_data

