! Initialization routines for the terminator 'toy' chemistry as used in DCMIP2016.
!
! Literature
! Dynamical Core Model Intercomparison Project (DCMIP2016) Test Case Document
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

MODULE mo_nh_dcmip_terminator

  USE mo_kind,                ONLY: wp, i8
  USE mo_impl_constants,      ONLY: min_rlcell, SUCCESS
  USE mo_math_constants,      ONLY: deg2rad, rad2deg
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_parallel_config,     ONLY: nproma
  USE mo_dynamics_config,     ONLY: nnow_rcf, nnew_rcf, nnow, nnew
  USE mo_run_config,          ONLY: msg_level
  USE mo_advection_config,    ONLY: advection_config
  USE mo_exception,           ONLY: finish, message_text, message
  USE mo_nh_testcases_nml,    ONLY: toy_chem
  USE mo_action_types,        ONLY: t_var_action_element, new_action
  USE mtime,                  ONLY: event, newEvent, isCurrentEventActive, &
    &                               newDatetime, newTimedelta, timedelta,  &
    &                               datetime, getPTStringFromMS,           &
    &                               getPTStringFromSeconds,                &
    &                               MAX_EVENTNAME_STR_LEN, MAX_DATETIME_STR_LEN

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nh_dcmip_terminator'

  TYPE(event), POINTER       :: chem_event, cpl_event
  TYPE(t_var_action_element) :: chem_action, cpl_action 

  ! time rate of change of cl and cl2
  REAL(wp), ALLOCATABLE, SAVE :: ddtcl (:,:,:), &
    &                            ddtcl2(:,:,:)

  !=======================================================================
  !    Test case parameters
  !=======================================================================
  REAL(wp), PARAMETER :: cly_constant  = 4.e-6_wp
  REAL(wp), PARAMETER :: k1_lat_center = 20._wp*deg2rad
  REAL(wp), PARAMETER :: k1_lon_center = 300._wp*deg2rad



  PUBLIC :: init_nh_dcmip_terminator
  PUBLIC :: dcmip_terminator_interface

CONTAINS


  !>
  !! Setup initial conditions for DCMIP2016 terminator 'toy' chemistry.
  !!
  SUBROUTINE init_nh_dcmip_terminator (p_patch, p_metrics, p_nh_prog, p_nh_diag)

    TYPE(t_patch),        INTENT(INOUT) :: &  !< patch on which computation is performed
      &  p_patch

    TYPE(t_nh_metrics),   INTENT(IN)    :: &  !< NH metrics state
      &  p_metrics

    TYPE(t_nh_prog),      INTENT(INOUT) :: &  !< prognostic state vector (now)
      &  p_nh_prog(:)

    TYPE(t_nh_diag),      INTENT(INOUT) :: &  !< diagnostic state vector
      &  p_nh_diag

    ! local
    CHARACTER(LEN = *), PARAMETER :: routine = modname//':init_nh_dcmip_terminator'
    INTEGER :: jg                      ! patch ID
    INTEGER :: jc, jb                  ! loop indices for cell, edge, block
    INTEGER :: i_rlstart, i_rlend, i_nchdom
    INTEGER :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER :: nlev                        ! number of vertical (full) levels
    INTEGER :: ist                         ! status flag

    REAL(wp) :: zlon, zlat                 ! lat/lon of cell circumcenter
    REAL(wp) :: zcl, zcl2                  ! passive tracer mixing ratios

    CHARACTER(LEN=MAX_EVENTNAME_STR_LEN):: event_name
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: ptstr       ! PT string

    INTEGER, PARAMETER :: ACTION_CHEM = 98 ! arbitrary action id
    INTEGER, PARAMETER :: ACTION_CPL  = 99 ! arbitrary action id


    jg = p_patch%id

    nlev = p_patch%nlev

    ! Sanity checks
    !
    ! make sure that a sufficient number of tracer fields is allocated
    IF (.NOT. ASSOCIATED(p_nh_prog(nnow_rcf(jg))%tracer)) THEN
      CALL finish (routine, 'Tracer field not allocated')
    ENDIF
    IF (advection_config(jg)%npassive_tracer < 2) THEN
      CALL finish(routine, 'Testcase requires allocation of 2 passive tracers')
    ENDIF

    ! sanity checks for diagnostics 
    ! make sure that a minimum number of extra_2d fields are allocated
    !
    IF (.NOT. ASSOCIATED(p_nh_diag%extra_2d)) THEN
      CALL finish (routine, 'No extra_2d fields allocated. A minimum of 3 is required!')
    ENDIF
    IF (SIZE(p_nh_diag%extra_2d,3) < 3) THEN
      WRITE(message_text,'(a,i2,a)') 'Numer of extra_2d fields (',SIZE(p_nh_diag%extra_2d,3),') smaller than minimum (3)'
      CALL finish (routine, message_text)
    ENDIF


    ! Print some information about tracer IDs assigned to CL and CL2
    WRITE(message_text,'(a,i2)') 'Tracer ID assigned for CL ', toy_chem%id_cl
    CALL message(TRIM(routine),message_text)
    WRITE(message_text,'(a,i2)') 'Tracer ID assigned for CL2 ', toy_chem%id_cl2
    CALL message(TRIM(routine),message_text)


    ! allocate time rate of change for cl and cl2
    ALLOCATE(ddtcl (nproma,p_patch%nlev,p_patch%nblks_c), &
      &      ddtcl2(nproma,p_patch%nlev,p_patch%nblks_c), stat=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(routine,'allocation for ddtcl, ddtcl2 failed')
    ENDIF
    ddtcl (:,:,:) = 0._wp
    ddtcl2(:,:,:) = 0._wp



    i_rlstart = 1
    i_rlend   = min_rlcell

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

    !
    ! Initialize two chemically reactive, though otherwise passive tracers
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,zlon,zlat,zcl,zcl2)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)


      DO jc = i_startidx, i_endidx

          ! get geographical coordinates of cell circumcenter
          !
          zlon = rad2deg*p_patch%cells%center(jc,jb)%lon
          zlat = rad2deg*p_patch%cells%center(jc,jb)%lat
          CALL initial_value_Terminator( lat = zlat, & !in
            &                            lon = zlon, & !in
            &                            cl  = zcl , & !out
            &                            cl2 = zcl2  ) !out

          p_nh_prog(nnow_rcf(jg))%tracer(jc,1:nlev,jb,toy_chem%id_cl)  = zcl
          p_nh_prog(nnow_rcf(jg))%tracer(jc,1:nlev,jb,toy_chem%id_cl2) = zcl2 
      ENDDO  ! jc
    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL


    ! Setup events for
    ! * chemistry tendency update
    ! * transport-chemistry coupling
    CALL getPTStringFromSeconds(toy_chem%dt_chem, ptstr)
    chem_action = new_action(ACTION_CHEM, ptstr)  !"PT10M")
    event_name = 'compute chem tendency'
    chem_event =>newEvent(TRIM(event_name),   &
              &           chem_action%ref,    &
              &           chem_action%start,  &
              &           chem_action%end  ,  &
              &           chem_action%intvl)
    !
    CALL getPTStringFromSeconds(toy_chem%dt_cpl, ptstr)
    cpl_action = new_action(ACTION_CPL, ptstr)  !"PT10M")
    event_name = 'Transport-chemistry coupling'
    cpl_event =>newEvent(TRIM(event_name),    &
              &          cpl_action%ref,      &
              &          cpl_action%start,    &
              &          cpl_action%end  ,    &
              &          cpl_action%intvl)



    ! Call diagnostics routine for initial conditions
    CALL dcmip_terminator_diag_vint (p_patch   = p_patch,                        & !in
      &                              tracer    = p_nh_prog(nnow_rcf(jg))%tracer, & !in
      &                              rho       = p_nh_prog(nnow(jg))%rho,        & !in
      &                              dz        = p_metrics%ddqz_z_full,          & !in
      &                              p_nh_diag = p_nh_diag                       ) !inout

  END SUBROUTINE init_nh_dcmip_terminator



  !>
  !! Interface between transport and terminator 'toy' chemistry forcing.
  !! New chemistry tendencies are computed every dt_chem.
  !! The chemistry-transport coupling is performed every dt_cpl.
  !!
  SUBROUTINE dcmip_terminator_interface (p_patch, p_metrics, p_nh_prog, p_nh_diag, &
    &                                    cur_datetime, dtime)

    TYPE(t_patch),        INTENT(IN)    :: &  !< patch on which computation is performed
      &  p_patch

    TYPE(t_nh_metrics),   INTENT(IN)    :: &  !< NH metrics state
      &  p_metrics

    TYPE(t_nh_prog),      INTENT(INOUT) :: &  !< prognostic state vector
      &  p_nh_prog(:)

    TYPE(t_nh_diag),      INTENT(INOUT) :: &  !< diagnostic state vector
      &  p_nh_diag

    TYPE(datetime), POINTER, INTENT(INOUT)  :: &  !< current datetime
      &  cur_datetime

    REAL(wp),             INTENT(IN)    :: &  !< physics time step [s]
      &  dtime

    ! local
    CHARACTER(LEN = *), PARAMETER :: routine = modname//':dcmip_terminator_interface'
    INTEGER :: jg                      ! patch ID
    INTEGER :: jc, jk, jb              ! loop indices for cell, edge, level, block
    INTEGER :: i_rlstart, i_rlend, i_nchdom
    INTEGER :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER :: nlev                        ! number of vertical (full) levels

    REAL(wp) :: zlon, zlat                 ! lat/lon of cell circumcenter
    LOGICAL :: chem_event_isactive, cpl_event_isactive

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: str_slack       ! slack as string

    TYPE(timedelta), POINTER :: p_slack                    ! slack in 'timedelta'-Format

    TYPE(datetime),   POINTER :: &         !< Current date in mtime format
      &  mtime_date

    ! check whether the chemistry and/or the coupling event is active
    !
    mtime_date => newDatetime(cur_datetime)
    !
    ! compute allowed slack in PT-Format
    ! Use factor 999 instead of 1000, since no open interval is available
    ! needed [trigger_date, trigger_date + slack[
    ! used   [trigger_date, trigger_date + slack]
    CALL getPTStringFromMS(INT(999.0_wp*dtime,i8),str_slack)
    ! get slack in 'timedelta'-format appropriate for isCurrentEventActive
    p_slack => newTimedelta(str_slack)
    !
    chem_event_isactive = LOGICAL(isCurrentEventActive(chem_event,mtime_date, plus_slack=p_slack))
    cpl_event_isactive  = LOGICAL(isCurrentEventActive(cpl_event ,mtime_date, plus_slack=p_slack))

    IF (msg_level >= 8) THEN
      WRITE(message_text,'(a,l1)') 'chem_event_isactive: ', chem_event_isactive
      CALL message("toy chemistry: ",message_text)
      WRITE(message_text,'(a,l1)') 'cpl_event_isactive: ', cpl_event_isactive
      CALL message("toy chemistry: ",message_text)
    ENDIF


    jg = p_patch%id

    nlev = p_patch%nlev

    i_rlstart = 1
    i_rlend   = min_rlcell

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

    !
    ! compute updated tracer tendencies
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,zlon,zlat)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      IF (chem_event_isactive) THEN
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            ! get geographical coordinates of cell circumcenter
            !
            zlon = rad2deg*p_patch%cells%center(jc,jb)%lon
            zlat = rad2deg*p_patch%cells%center(jc,jb)%lat
            CALL tendency_Terminator( lat   = zlat,             &  !in
              &                       lon   = zlon,             &  !in
              &                       cl    = p_nh_prog(nnew_rcf(jg))%tracer(jc,jk,jb,toy_chem%id_cl),  & !in
              &                       cl2   = p_nh_prog(nnew_rcf(jg))%tracer(jc,jk,jb,toy_chem%id_cl2), & !in 
              &                       dt    = toy_chem%dt_chem, &  !in
              &                       cl_f  = ddtcl(jc,jk,jb),  &  !out
              &                       cl2_f = ddtcl2(jc,jk,jb)  )  !out
          ENDDO  ! jc
        ENDDO  ! jk
      ENDIF  ! chem_event_isactive
    
      ! apply tracer tendencies
      IF (cpl_event_isactive) THEN
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            p_nh_prog(nnew_rcf(jg))%tracer(jc,jk,jb,toy_chem%id_cl) =                 &
              &    MAX(0._wp,p_nh_prog(nnew_rcf(jg))%tracer(jc,jk,jb,toy_chem%id_cl)  &
              &  + toy_chem%dt_cpl*ddtcl(jc,jk,jb) )
            p_nh_prog(nnew_rcf(jg))%tracer(jc,jk,jb,toy_chem%id_cl2) =                &
              &    MAX(0._wp,p_nh_prog(nnew_rcf(jg))%tracer(jc,jk,jb,toy_chem%id_cl2) &
              &  + toy_chem%dt_cpl*ddtcl2(jc,jk,jb))
          ENDDO  ! jc
        ENDDO  ! jk
      ENDIF  ! cpl_event_isactive

    END DO  ! jb
!$OMP END DO
!$OMP END PARALLEL


    ! toy chemistry-specific diagnostics
    CALL dcmip_terminator_diag_vint (p_patch   = p_patch,                        & !in
      &                              tracer    = p_nh_prog(nnew_rcf(jg))%tracer, & !in
      &                              rho       = p_nh_prog(nnew(jg))%rho,        & !in
      &                              dz        = p_metrics%ddqz_z_full,          & !in
      &                              p_nh_diag = p_nh_diag                       ) !inout

  END SUBROUTINE dcmip_terminator_interface


  !>
  !! Diagnose vertical integrals of Cl, Cl2, and Cly = Cl + 2Cl2
  !! makes use of 3 extra_2d fields
  !! extra_2d(:,:,1) : Cl
  !! extra_2d(:,:,2) : Cl2
  !! extra_2d(:,:,3) : Cl + 2Cl2
  !!
  SUBROUTINE dcmip_terminator_diag_vint (p_patch, tracer, rho, dz, p_nh_diag)

    TYPE(t_patch),        INTENT(IN)    :: &  !< patch on which computation is performed
      &  p_patch

    REAL(wp),             INTENT(IN)    :: &  !< total (moist) density (updated state)
      &  rho(:,:,:)

    REAL(wp),             INTENT(IN)    :: &  !< tracer container (updated state)
     &  tracer(:,:,:,:)

    REAL(wp),             INTENT(IN)    :: &  !< full level thickness
     &  dz(:,:,:)

    TYPE(t_nh_diag),      INTENT(INOUT) :: &  !< diagnostic state vector
      &  p_nh_diag

    ! local
    CHARACTER(LEN = *), PARAMETER :: routine = modname//':dcmip_terminator_diag_vint'
    INTEGER :: jc, jk, jb                  ! loop indices for cell, edge, level, block
    INTEGER :: i_rlstart, i_rlend, i_nchdom
    INTEGER :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER :: nlev                        ! number of vertical (full) levels
    REAL(wp):: mass                        ! column integrated mass [kg m-2]
    REAL(wp):: rhodz(nproma,p_patch%nlev)  ! rho times delta z


    nlev = p_patch%nlev

    i_rlstart = 1
    i_rlend   = min_rlcell

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,rhodz,mass)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart, i_rlend)

      ! pre-computation of rho * \Delta z
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx 
          rhodz(jc,jk) = dz(jc,jk,jb) * rho(jc,jk,jb)  
        ENDDO
      ENDDO

      p_nh_diag%extra_2d(i_startidx:i_endidx,jb,:) = 0.0_wp

      DO jk = 1, nlev

        DO jc = i_startidx, i_endidx 

          ! vertical intergal of Cl
          p_nh_diag%extra_2d(jc,jb,1) = p_nh_diag%extra_2d(jc,jb,1)  &
            &                + rhodz(jc,jk) * tracer(jc,jk,jb,toy_chem%id_cl)
 
          ! vertical intergal of Cl2
          p_nh_diag%extra_2d(jc,jb,2) = p_nh_diag%extra_2d(jc,jb,2)  &
            &                + rhodz(jc,jk) * tracer(jc,jk,jb,toy_chem%id_cl2)

          ! vertical intergal of Cly
          p_nh_diag%extra_2d(jc,jb,3) = p_nh_diag%extra_2d(jc,jb,3)  &
            &  + rhodz(jc,jk) * (tracer(jc,jk,jb,toy_chem%id_cl) &
            &  + 2._wp*tracer(jc,jk,jb,toy_chem%id_cl2))
        ENDDO  ! jc
      ENDDO  ! jk

      ! compute vertical average
      DO jc = i_startidx, i_endidx
        mass = SUM(rhodz(jc,1:nlev))
        p_nh_diag%extra_2d(jc,jb,1:3) = p_nh_diag%extra_2d(jc,jb,1:3) / mass
      ENDDO
    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE dcmip_terminator_diag_vint



  ! The following routines were provided by the DCMIP2016 organizers.
  !=======================================================================

  !===============================================================================
  !  Compute initial values
  !===============================================================================

  SUBROUTINE initial_value_Terminator( lat, lon, cl, cl2 )

    !-----------------------------------------------------------------------
    ! Arguments:
    !-----------------------------------------------------------------------

    REAL(wp), intent(in)  :: lat, lon  ! latitude and longitude, degrees
    REAL(wp), intent(out) :: cl, cl2   ! molar mixing ratio of cl and cl2

    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------

    REAL(wp) :: r, det  ! useful algebraic forms
    REAL(wp) :: k1, k2  ! reaction rates

    CALL k_vals( lat*deg2rad, lon*deg2rad, k1, k2 )

    r = k1 / (4._wp*k2)
    det = sqrt(r*r + 2._wp*cly_constant*r)

    cl  = (det-r)
    cl2 = cly_constant/2._wp - (det-r)/2._wp

    RETURN

  END SUBROUTINE initial_value_Terminator


  !===============================================================================
  !  Solar photolysis rate and recombination rate
  !===============================================================================

  SUBROUTINE k_vals( lat, lon, k1, k2 )

    !-----------------------------------------------------------------------
    ! Arguments:
    !-----------------------------------------------------------------------
    REAL(wp), intent(in)    :: lat, lon  ! latitude and longitude, radians
    REAL(wp), intent(out)   :: k1, k2    ! reaction rates

    k1 = 1.0_wp*max(0._wp,sin(lat)*sin(k1_lat_center) + cos(lat)*cos(k1_lat_center)*cos(lon-k1_lon_center))
    k2 = 1._wp

  RETURN

  END SUBROUTINE k_vals


  !===============================================================================
  !  Tendencies of cl and cl2
  !===============================================================================

  SUBROUTINE tendency_Terminator( lat, lon, cl, cl2, dt, cl_f, cl2_f )

    !-----------------------------------------------------------------------
    ! Arguments:
    !-----------------------------------------------------------------------

    REAL(wp), intent(in)    :: lat, lon  ! latitude and longitude, degrees
    REAL(wp), intent(in)    :: cl, cl2   ! molar mixing ratio of cl and cl2
    REAL(wp), intent(in)    :: dt        ! size of physics time step

    REAL(wp), intent(out)   :: cl_f, cl2_f  ! time rate of change of cl and cl2

    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------

    REAL(wp) :: r, det, expdt, el ! useful algebaic quantities used in the computation
    REAL(wp) :: k1, k2            ! reaction rates
    REAL(wp) :: cly               ! quantity that should be conseved

    CALL k_vals( lat*deg2rad, lon*deg2rad, k1, k2 )

    r = k1 / (4._wp*k2)
    cly = cl + 2._wp* cl2

    det = sqrt( r*r + 2._wp*r*cly )
    expdt = exp( -4._wp*k2*det*dt )

    IF ( abs(det * k2 * dt) .gt. 1e-16 ) THEN
      el = (1._wp - expdt) /det /dt
    ELSE
      el = 4._wp*k2
    ENDIF

    cl_f  = -el * (cl - det + r)*(cl + det + r) / (1._wp + expdt + dt*el*(cl + r))
    cl2_f = -cl_f / 2._wp

    RETURN

  END SUBROUTINE tendency_Terminator


END MODULE mo_nh_dcmip_terminator

