!
! Processing of the external data for the upper atmosphere.
!
!----------------------------
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_upatmo_extdat

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish, message
  USE mo_impl_constants,         ONLY: SUCCESS, MAX_CHAR_LENGTH, &
    &                                  min_rlcell_int, max_dom_dig10
  USE mo_impl_constants_grf,     ONLY: grf_bdywidth_c
  USE mo_upatmo_impl_const,      ONLY: iUpatmoExtdatId, iUpatmoPrcId
  USE mo_model_domain,           ONLY: t_patch
  USE mo_upatmo_types,           ONLY: t_upatmo_extdat, t_upatmo_tend
  USE mo_upatmo_config,          ONLY: t_upatmo_config
  USE mo_loopindices,            ONLY: get_indices_c
  USE mtime,                     ONLY: datetime, newDatetime, deallocateDatetime
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights, &
    &                                  calculate_time_interpolation_weights
  USE mo_util_string,            ONLY: int2string

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: update_upatmo_extdat_nwp

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_upatmo_extdat'

CONTAINS

  !>
  !! Update external data for the upper atmosphere 
  !! under NWP forcing.
  !!
  SUBROUTINE update_upatmo_extdat_nwp( mtime_datetime,    &  !in
    &                                  p_patch,           &  !in
    &                                  prm_upatmo_extdat, &  !inout
    &                                  prm_upatmo_tend,   &  !inout
    &                                  upatmo_config,     &  !in
    &                                  lupdate,           &  !in
    &                                  lmessage           )  !in

    ! In/out variables
    TYPE(datetime),                    INTENT(IN)    :: mtime_datetime     ! Date/time information
    TYPE(t_patch),            TARGET,  INTENT(IN)    :: p_patch            ! Grid/patch info
    TYPE(t_upatmo_extdat),    TARGET,  INTENT(INOUT) :: prm_upatmo_extdat  ! Upper-atmosphere external data
    TYPE(t_upatmo_tend),      TARGET,  INTENT(INOUT) :: prm_upatmo_tend    ! Upper-atmosphere physics tendencies
    TYPE(t_upatmo_config),             INTENT(IN)    :: upatmo_config      ! General upper-atmosphere configuration
    LOGICAL,                           INTENT(IN)    :: lupdate(:)         ! Update switches
    LOGICAL,                           INTENT(IN)    :: lmessage           ! Switch for message output

    ! Local variables
    REAL(wp), ALLOCATABLE :: ext_intrpl_time(:,:), ext_intrpl_lev(:,:)
    REAL(wp), POINTER     :: ext_data_time1(:,:), ext_data_time2(:,:), &
      &                      wgt_lev1(:), wgt_lev2(:),                 &
      &                      wgt_lat1(:,:), wgt_lat2(:,:),             &
      &                      intrpl_rslt(:,:,:)

    INTEGER, POINTER :: ilev1(:), ilev2(:), ilat1(:,:), ilat2(:,:)
    INTEGER  :: jg
    INTEGER  :: nlat_ext, nlev_ext, nlev, ngas_ext
    INTEGER  :: istartlat, iendlat, isteplat
    INTEGER  :: istartlev, iendlev, isteplev
    INTEGER  :: istartlev_chemheat, iendlev_chemheat
    INTEGER  :: jlat, jlev, jk, jb, jc, jgas
    INTEGER  :: rl_start, rl_end
    INTEGER  :: i_startblk, i_endblk 
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: istat

    CHARACTER(LEN=max_dom_dig10) :: dom_str

    TYPE(datetime) :: mtime_hour
    TYPE(t_time_interpolation_weights) :: time_intrpl

    CHARACTER(LEN=*), PARAMETER ::  &
      &  routine = modname//':update_upatmo_extdat_nwp'
    
    !--------------------------------------------------------------

    ! Some notes:
    !
    ! * For the time being, all external data for the upper atmosphere 
    !   under NWP forcing are monthly mean values, which we have to interpolate.
    !   The interpolation procedure is the same for chemical heating tendencies 
    !   and the various radiatively active gases, so that its outsourcing 
    !   into a separate subroutine might stand to reason. 
    !   Nevertheless, we think that the accompanied tightening up of the code 
    !   does not yet outweigh the loss in performance due to potentially several 
    !   invocations of that subroutine.
    ! 
    ! * The time resolution of the external data is monthly, 
    !   so the only reason to choose a shorter time scale
    !   for the update period of their time interpolation is 
    !   to avoid too strong jumps. 
    !   An update period on the order of a day may be recommendable. 
    !   Any shorter period is probably a wast of computational resources.
    !
    ! * The interpolation in time is a conceptual copy of 
    !   'src/shr_horizontal/mo_ext_data_init: interpol_monthly_mean'
    !
    ! * If we denote the number of horizontal and vertical grid points 
    !   of the external data by m and n, 
    !   and the number of horizontal and vertical grid points of ICON 
    !   by M and N, and assume that m > n, M > N, M > m and N > n, 
    !   the following interpolation order should come along with 
    !   the least number of operations:
    !   
    !   time  ->  vertical  ->  horizontal.
    !
    !   But this depends crucially on the above assumptions.
    !
    ! * Currently, external chemical heating tendencies are provided 
    !   on geometric heights, so we can do the vertical interpolation 
    !   onto the geometric heights of the ICON grid once here 
    !   and that interpolation is valid until the next call of this subroutine.
    ! 
    ! * Currently, external gas data are provided on pressure levels. 
    !   Their geometric vertical position may change 
    !   more or less strongly from time t to t + dt_fastphy. 
    !   Typically, the update period for the time interpolation 
    !   of the external data is much too long 
    !   to resolve this variation appropriately. 
    !   This is the reason, why the vertical interpolation 
    !   of the external gas data is not done here, 
    !   but in 'src/upper_atmosphere/mo_upatmo_phy_diag: update_diagnostic_variables'.

    ! Domain index
    jg = p_patch%id

    ! Number of vertical levels
    nlev = p_patch%nlev

    IF (lmessage) WRITE (dom_str, '(i0)') jg

    IF (lmessage) CALL message(routine, &
      & 'Start update of external data on domain '//dom_str)

    !---------------------------------------------------------------------
    !                            Preparation
    !---------------------------------------------------------------------

    !---------------------------
    ! For interpolation in time
    !---------------------------

    mtime_hour = mtime_datetime
    mtime_hour%time%minute = 0
    mtime_hour%time%second = 0
    mtime_hour%time%ms     = 0     
    time_intrpl = calculate_time_interpolation_weights(mtime_hour)
    
    !---------------------------------------------------------------------
    !                        Update external data:
    !                     Chemical heating tendencies
    !---------------------------------------------------------------------

    IF (lupdate( iUpatmoExtdatId%chemheat )) THEN

      nlat_ext  = prm_upatmo_extdat%chemheat%nlat
      nlev_ext  = prm_upatmo_extdat%chemheat%nlev
      istartlat = prm_upatmo_extdat%chemheat%istartlat
      iendlat   = prm_upatmo_extdat%chemheat%iendlat
      isteplat  = prm_upatmo_extdat%chemheat%isteplat
      istartlev = prm_upatmo_extdat%chemheat%istartlev
      iendlev   = prm_upatmo_extdat%chemheat%iendlev 
      isteplev  = prm_upatmo_extdat%chemheat%isteplev

      ! Set some convenience pointers:
      ! * For time interpolation
      ext_data_time1 => prm_upatmo_extdat%chemheat%data(:,:,time_intrpl%month1)
      ext_data_time2 => prm_upatmo_extdat%chemheat%data(:,:,time_intrpl%month2)
      ! * For vertical interpolation
      ilev1    => prm_upatmo_extdat%chemheat%intrpl%lev%idx(1,:)
      ilev2    => prm_upatmo_extdat%chemheat%intrpl%lev%idx(2,:)
      wgt_lev1 => prm_upatmo_extdat%chemheat%intrpl%lev%wgt(1,:)
      wgt_lev2 => prm_upatmo_extdat%chemheat%intrpl%lev%wgt(2,:)
      ! * For horizontal interpolation
      ilat1    => prm_upatmo_extdat%chemheat%intrpl%lat%idx(1,:,:)
      ilat2    => prm_upatmo_extdat%chemheat%intrpl%lat%idx(2,:,:)
      wgt_lat1 => prm_upatmo_extdat%chemheat%intrpl%lat%wgt(1,:,:)
      wgt_lat2 => prm_upatmo_extdat%chemheat%intrpl%lat%wgt(2,:,:)
      ! * For the final interpolation result
      intrpl_rslt => prm_upatmo_tend%ddt_temp_chemheat

      !-----------------------
      ! Interpolation in time
      !-----------------------

      ALLOCATE( ext_intrpl_time( nlat_ext, nlev_ext ), &
        &       ext_intrpl_lev( nlat_ext, nlev ),      &
        &       STAT=istat                             )
      IF(istat /= SUCCESS) CALL finish(routine, 'Allocation of ext_intrpl_time/lev failed.')

      ! Currently, the value of 'nlev_ext' is not very large, 
      ! so we might not yet profit from parallelizing the vertical loop
      DO jlev = istartlev, iendlev, isteplev
        DO jlat = istartlat, iendlat, isteplat
          ext_intrpl_time(jlat,jlev) = time_intrpl%weight1 * ext_data_time1(jlat,jlev) &
            &                        + time_intrpl%weight2 * ext_data_time2(jlat,jlev)
        ENDDO  !jlat
      ENDDO  !jlev

      NULLIFY(ext_data_time1, ext_data_time2)

      !------------------------
      ! Vertical interpolation
      !------------------------

      ! Please note that the vertical interpolation is onto 
      ! the nominal grid layer heights, not onto the actual heights 
      ! of the grid cells for reasons mentioned in 
      ! 'src/upper_atmosphere/mo_upatmo_extdat_state: construct_upatmo_extdat_nwp'.

      DO jk = 1, nlev
        DO jlat = istartlat, iendlat, isteplat
          ext_intrpl_lev(jlat,jk) = wgt_lev1(jk) * ext_intrpl_time(jlat,ilev1(jk)) &
            &                     + wgt_lev2(jk) * ext_intrpl_time(jlat,ilev2(jk))
        ENDDO  !jlat
      ENDDO  !jlev

      DEALLOCATE(ext_intrpl_time, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of ext_intrpl_time failed.')

      !--------------------------
      ! Horizontal interpolation
      !--------------------------

      ! Start level above which and end level 
      ! below which temperature tendencies 
      ! from chemical heating are set to zero.
      istartlev_chemheat = upatmo_config%nwp_phy%prc( iUpatmoPrcId%chemheat )%istartlev
      iendlev_chemheat   = upatmo_config%nwp_phy%prc( iUpatmoPrcId%chemheat )%iendlev

      ! Loop boundaries for prognostic domain.
      rl_start   = grf_bdywidth_c + 1
      rl_end     = min_rlcell_int
      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_GUIDED_SCHEDULE
      DO jb = i_startblk, i_endblk
        
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
        
        DO jk = 1, istartlev_chemheat - 1
          DO jc = i_startidx, i_endidx      
            intrpl_rslt(jc,jk,jb) = 0._wp
          ENDDO  !jc
        ENDDO  !jk
        
        DO jk = istartlev_chemheat, iendlev_chemheat
          DO jc = i_startidx, i_endidx      
            intrpl_rslt(jc,jk,jb) = wgt_lat1(jc,jb) * ext_intrpl_lev(ilat1(jc,jb),jk) &
              &                   + wgt_lat2(jc,jb) * ext_intrpl_lev(ilat2(jc,jb),jk)
          ENDDO  !jc
        ENDDO  !jk
        
        DO jk = iendlev_chemheat + 1, nlev
          DO jc = i_startidx, i_endidx      
            intrpl_rslt(jc,jk,jb) = 0._wp
          ENDDO  !jc
        ENDDO  !jk        
      ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL      
      
      !----------
      ! Clean-up
      !----------

      DEALLOCATE(ext_intrpl_lev, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of ext_intrpl_lev failed.')

      NULLIFY( ilev1, ilev2, wgt_lev1, wgt_lev2, ilat1, ilat2, wgt_lat1, wgt_lat2, &
        &      ext_data_time1, ext_data_time2, intrpl_rslt                         )

    ENDIF  !Update of external chemical heating tendencies due?

    !---------------------------------------------------------------------
    !                        Update external data:
    !                      Radiatively active gases
    !---------------------------------------------------------------------

    IF (lupdate( iUpatmoExtdatId%gases )) THEN

      ! Number of gases, whose concentrations are provided as external data
      ngas_ext = prm_upatmo_extdat%ngas

      IF (ngas_ext > 0) THEN
        nlat_ext  = prm_upatmo_extdat%gas(1)%nlat
        nlev_ext  = prm_upatmo_extdat%gas(1)%nlev
        ALLOCATE(ext_intrpl_time( nlat_ext, nlev_ext ), STAT=istat)
        IF(istat /= SUCCESS) CALL finish(routine, 'Allocation of ext_intrpl_time failed.')
      END IF
      ! Loop over these gases
      DO jgas = 1, ngas_ext

        nlat_ext  = prm_upatmo_extdat%gas( jgas )%nlat
        nlev_ext  = prm_upatmo_extdat%gas( jgas )%nlev
        istartlat = prm_upatmo_extdat%gas( jgas )%istartlat
        iendlat   = prm_upatmo_extdat%gas( jgas )%iendlat
        isteplat  = prm_upatmo_extdat%gas( jgas )%isteplat
        istartlev = prm_upatmo_extdat%gas( jgas )%istartlev
        iendlev   = prm_upatmo_extdat%gas( jgas )%iendlev 
        isteplev  = prm_upatmo_extdat%gas( jgas )%isteplev

        ext_data_time1 => prm_upatmo_extdat%gas( jgas )%data(:,:,time_intrpl%month1)
        ext_data_time2 => prm_upatmo_extdat%gas( jgas )%data(:,:,time_intrpl%month2)
        ilat1          => prm_upatmo_extdat%gas( jgas )%intrpl%lat%idx(1,:,:)
        ilat2          => prm_upatmo_extdat%gas( jgas )%intrpl%lat%idx(2,:,:)
        wgt_lat1       => prm_upatmo_extdat%gas( jgas )%intrpl%lat%wgt(1,:,:)
        wgt_lat2       => prm_upatmo_extdat%gas( jgas )%intrpl%lat%wgt(2,:,:)
        intrpl_rslt    => prm_upatmo_extdat%gas_interm( jgas )%p

        !-----------------------
        ! Interpolation in time
        !-----------------------

        IF (nlat_ext /= SIZE(ext_intrpl_time, 1) &
          & .OR. nlev_ext /= SIZE(ext_intrpl_time, 2)) THEN
          DEALLOCATE(ext_intrpl_time, STAT=istat)
          IF(istat == SUCCESS) ALLOCATE(ext_intrpl_time( nlat_ext, nlev_ext ), STAT=istat)
          IF(istat /= SUCCESS) CALL finish(routine, 'Reallocation of ext_intrpl_time failed.')
        END IF

!$OMP PARALLEL PRIVATE(rl_start, rl_end, i_startblk, i_endblk)
!$OMP DO PRIVATE(JLEV, JLAT)
        DO jlev = istartlev, iendlev, isteplev
          DO jlat = istartlat, iendlat, isteplat
            ext_intrpl_time(jlat,jlev) = time_intrpl%weight1 * ext_data_time1(jlat,jlev) &
              &                        + time_intrpl%weight2 * ext_data_time2(jlat,jlev)
          ENDDO  !jlat
        ENDDO  !jlev
!$OMP END DO

        !--------------------------
        ! Horizontal interpolation
        !--------------------------
        
        ! Loop boundaries for prognostic domain.
        rl_start   = grf_bdywidth_c + 1
        rl_end     = min_rlcell_int
        i_startblk = p_patch%cells%start_block(rl_start)
        i_endblk   = p_patch%cells%end_block(rl_end)
        
!$OMP DO PRIVATE(jb, jlev, jc, i_startidx, i_endidx) ICON_OMP_GUIDED_SCHEDULE
        DO jb = i_startblk, i_endblk
          
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
          
          DO jlev = istartlev, iendlev, isteplev
            DO jc = i_startidx, i_endidx      
              intrpl_rslt(jc,jlev,jb) = wgt_lat1(jc,jb) * ext_intrpl_time(ilat1(jc,jb),jlev) &
                &                     + wgt_lat2(jc,jb) * ext_intrpl_time(ilat2(jc,jb),jlev)
            ENDDO  !jc
          ENDDO  !jlev
        ENDDO  !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL      
        
        !----------
        ! Clean-up
        !----------


        NULLIFY(ilat1, ilat2, wgt_lat1, wgt_lat2, ext_data_time1, ext_data_time2, intrpl_rslt)
        
      ENDDO  !jgas

      DEALLOCATE(ext_intrpl_time, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of ext_intrpl_time failed.')

    ENDIF  !Update of external gas data due?

    IF (lmessage) CALL message(routine, &
      & 'Finish update of external data on domain '//dom_str)

  END SUBROUTINE update_upatmo_extdat_nwp

END MODULE mo_upatmo_extdat
