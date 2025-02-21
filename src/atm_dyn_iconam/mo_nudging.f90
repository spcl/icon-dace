!
! Nudging.
!
! This module contains procedures related to the nudging
! of the atmospheric state simulated by ICON towards driving data.
! Covered nudging types:
! - Global nudging
! For the nudging types:
! - Lateral boundary nudging
! - Upper boundary nudging
! please see:
! - src/atm_dyn_iconam/mo_nh_nest_utilities: limarea_bdy_nudging
! The nudging we refer to in this module is unrelated
! to the nudging in the context of:
! - assimilation (assimilation_nml)
! - large-scale forcing (ls_forcing_nml)
!
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

!----------------------------
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

MODULE mo_nudging

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish, message
  USE mo_impl_constants,        ONLY: SUCCESS, min_rlcell, min_rledge
  USE mo_impl_constants_grf,    ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_physical_constants,    ONLY: rd, cvd_o_rd, p0ref, vtmpc1, rcpd
  USE mo_run_config,            ONLY: iqv, iqc
  USE mo_parallel_config,       ONLY: nproma
  USE mo_nudging_config,        ONLY: t_nudging_config, indg_type, indg_var, ithermdyn_type
  USE mo_model_domain,          ONLY: t_patch
  USE mo_nonhydro_types,        ONLY: t_nh_state, t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_initicon_types,        ONLY: t_pi_atm
  USE mo_async_latbc_types,     ONLY: t_latbc_data
  USE mtime,                    ONLY: datetime
  USE mo_loopindices,           ONLY: get_indices_c, get_indices_e
  USE mo_nh_diagnose_pres_temp, ONLY: diagnose_pres_temp
  USE mo_timer,                 ONLY: timer_start, timer_stop, timer_global_nudging

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nudging_interface

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nudging'

CONTAINS !..................................................................

  !>
  !! Nudging interface (for global nudging).
  !!
  !! This subroutine is meant as some kind of wrapper 
  !! for the procedures, which make up one nudging cycle 
  !! (except for the procedures related to the I/O of the driving data).
  !!
  SUBROUTINE nudging_interface( p_patch,          & !in
    &                           p_nh_state,       & !inout
    &                           latbc,            & !in
    &                           mtime_datetime,   & !in
    &                           nnew,             & !in
    &                           nnew_rcf,         & !in
    &                           nudging_config    ) !inout

    ! In/out variables
    TYPE(t_patch),            TARGET,  INTENT(IN)    :: p_patch           !< Grid/patch info
    TYPE(t_nh_state),         TARGET,  INTENT(INOUT) :: p_nh_state        !< Prognostic and diagnostic variables etc.
    TYPE(t_latbc_data),       TARGET,  INTENT(INOUT) :: latbc             !< Data structure for async latbc prefetching
    TYPE(datetime),           POINTER, INTENT(IN)    :: mtime_datetime    !< Date/time information
    INTEGER,                           INTENT(IN)    :: nnew, nnew_rcf    !< Time level indices
    TYPE(t_nudging_config),            INTENT(INOUT) :: nudging_config    !< Nudging switches

    ! Local variables
    TYPE(t_pi_atm), POINTER :: p_latbc_old, p_latbc_new
    REAL(wp) :: wfac_old, wfac_new
    INTEGER  :: jg
    LOGICAL  :: l_thermdyn, l_hydrostatic, l_message
    CHARACTER(LEN=*), PARAMETER :: &
      & routine = modname//":nudging_interface"
    
    !----------------------------------------

    ! Some notes may be in order: 
    !
    ! - Although it comes at the expense of computational efficiency   
    !   (due to the interposition of this subroutine 
    !   before the actual call of the nudging procedures), 
    !   we introduced this interface, to reduce the optical overload 
    !   of 'src/atm_dyn_iconam/mo_nh_stepping'
    !
    ! - This subroutine is called in 'mo_nh_stepping: integrate_nh' 
    !   only if global nudging is switched on (NOT nudging in general!) 
    !   and if we are on the primary domain ('jg = 1')
    !
    ! - The content of this subroutine follows closely the infrastructure 
    !   around the call of 'limarea_bdy_nudging' in 'mo_nh_stepping: integrate_nh'
    !
    ! - However, there is a lot of nudging infrastructure in 'mo_nh_stepping', 
    !   which cannot be included in this interface: 
    !   * Initialization and finalization of the nudging driving data processing 
    !   * Triggering of the read-in and postprocessing of the driving data
    !   In most cases a modification of this infrastructure for our purposes 
    !   is relatively moderate (search for "l_global_nudging"). 
    !   The prefetching of boundary data in case of asynchronous read-in 
    !   (num_prefetch_proc > 0) is triggered, if 'latbc_config%itype_latbc > 0'. 
    !   This should also work for global nudging, since 'latbc_config%itype_latbc = 1' 
    !   (i.e. time-dependent driving data) is enforced 
    !   in 'src/namelists/mo_nudging_nml: check_nudging'.

    ! Domain index
    jg = p_patch%id

    ! Better check
    IF (jg /= 1) RETURN

    !...............................................................
    !                       Global nudging
    !...............................................................

    IF (nudging_config%ltype(indg_type%globn)) THEN 

      IF (nudging_config%ltimer) CALL timer_start(timer_global_nudging)

      l_message = nudging_config%lmessage
      IF(l_message) CALL message(routine, 'Start global nudging.')

      !---------------------------------------------------------------
      !                        Preparation
      !---------------------------------------------------------------

      !
      ! Asynchronous read-in of driving data:
      !
      ! The following subroutine updates the weights: 
      ! * latbc%lc1 
      ! * latbc%lc2 
      ! for interpolation of the driving data in time
      CALL latbc%update_intp_wgt(mtime_datetime)


      ! Set pointer to past and future state of driving data,
      ! from which their current state is estimated by 
      ! linear interpolation in time
      p_latbc_old => latbc%latbc_data( latbc%prev_latbc_tlev() )%atm
      p_latbc_new => latbc%latbc_data( latbc%new_latbc_tlev    )%atm

      ! Get weights for time interpolation of driving data
      wfac_old = latbc%lc1
      wfac_new = latbc%lc2
      
      !---------------------------------------------------------------
      !                 Update diagnostic variables
      !---------------------------------------------------------------
      
      ! If we nudge the hydrostatic thermodynamic variables, 
      ! the diagnostic variables pressure and (virtual) temperature are required 
      ! (for the levels 'nudging_config%ilev_start' to 'nudging_config%ilev_end')

      l_thermdyn    = nudging_config%lvar(indg_var%thermdyn) 
      l_hydrostatic = nudging_config%thermdyn_type == ithermdyn_type%hydrostatic
      
      IF (l_thermdyn .AND. l_hydrostatic) THEN
        ! ('diagnose_pres_temp' updates the prognostic, halo and lateral boundary cells)
        CALL diagnose_pres_temp( p_metrics         = p_nh_state%metrics,           & !in
          &                      pt_prog           = p_nh_state%prog( nnew     ),  & !in 
          &                      pt_prog_rcf       = p_nh_state%prog( nnew_rcf ),  & !in
          &                      pt_diag           = p_nh_state%diag,              & !out
          &                      pt_patch          = p_patch,                      & !in
          &                      opt_calc_temp     = .TRUE.,                       & !optin
          &                      opt_calc_pres     = .TRUE.,                       & !optin 
          &                      opt_calc_temp_ifc = .FALSE.,                      & !optin
          &                      opt_slev          = nudging_config%ilev_start     ) !optin        
      ENDIF  !Update press and temp(v)?

      !---------------------------------------------------------------
      !                        Apply nudging
      !---------------------------------------------------------------
      
      CALL global_nudging( p_patch        = p_patch,                     & !in
        &                  p_prog         = p_nh_state%prog( nnew     ), & !inout
        &                  p_prog_rcf     = p_nh_state%prog( nnew_rcf ), & !inout
        &                  p_metrics      = p_nh_state%metrics,          & !in
        &                  p_diag         = p_nh_state%diag,             & !in
        &                  p_latbc_old    = p_latbc_old,                 & !in
        &                  p_latbc_new    = p_latbc_new,                 & !in
        &                  wfac_old       = wfac_old,                    & !in
        &                  wfac_new       = wfac_new,                    & !in
        &                  nudging_config = nudging_config               ) !in

      !---------------------------------------------------------------
      !                         Clean-up
      !---------------------------------------------------------------

      p_latbc_old => NULL()
      p_latbc_new => NULL()

      IF(l_message) CALL message(routine, 'End of global nudging.')
      
      IF (nudging_config%ltimer) CALL timer_stop(timer_global_nudging)

    ENDIF  !IF (nudging_config%ltype(indg_type%globn))
    
  END SUBROUTINE nudging_interface

!---------------------------------------------------------------------------

  !>
  !! This routine executes global nudging.
  !!
  SUBROUTINE global_nudging( p_patch,       & !in
    &                        p_prog,        & !inout
    &                        p_prog_rcf,    & !inout
    &                        p_metrics,     & !in
    &                        p_diag,        & !in
    &                        p_latbc_old,   & !in
    &                        p_latbc_new,   & !in
    &                        wfac_old,      & !in
    &                        wfac_new,      & !in
    &                        nudging_config ) !in

    ! In/out variables
    TYPE(t_patch),          INTENT(IN)    :: p_patch
    TYPE(t_nh_prog),        INTENT(INOUT) :: p_prog, p_prog_rcf
    TYPE(t_nh_metrics),     INTENT(IN)    :: p_metrics
    TYPE(t_nh_diag),        INTENT(IN)    :: p_diag
    TYPE(t_pi_atm),         INTENT(IN)    :: p_latbc_old, p_latbc_new  !< Past and future state of driving data
    REAL(wp),               INTENT(IN)    :: wfac_old, wfac_new        !< Weights for time interpolation of driving data
    TYPE(t_nudging_config), INTENT(IN)    :: nudging_config            !< Nudging switches

    ! Local variables
    REAL(wp), ALLOCATABLE :: qv_tend(:,:,:)
    REAL(wp), ALLOCATABLE :: nudge_coeff_thermdyn(:), &
      &                      nudge_coeff_vn(:),       &
      &                      nudge_coeff_qv(:)
    REAL(wp) :: rho_tend, thv_tend, vn_tend
    REAL(wp) :: pres_tend, temp_tend, tempv_tend, qv_update
    REAL(wp) :: nudge_coeff
    INTEGER  :: jg, jc, je, jb, jk
    INTEGER  :: istart, iend
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: rl_start, rl_end
    INTEGER  :: istat
    LOGICAL  :: l_thermdyn, l_vn, l_qv, l_hydrostatic, l_qv_tend
    REAL(wp), PARAMETER :: rd_o_cvd    = 1._wp / cvd_o_rd
    REAL(wp), PARAMETER :: rd_o_p0ref  = rd / p0ref
    REAL(wp), PARAMETER :: rrd         = 1._wp / rd
    REAL(wp), PARAMETER :: eps_qc      = 1.e-10_wp
    CHARACTER(LEN=*), PARAMETER :: &
      & routine = modname//":global_nudging"

    !----------------------------------------

    !........................................................................................
    ! Some notes:
    ! 
    ! - Global nudging is applied to the prognostic and the halo partitions of the domain, 
    !   so we assume that the prognostic variables 'p_prog' as well as the driving data 
    !   on 'p_latbc_old/new' have been synchronized before they enter this subroutine
    !
    ! - Some diagnostic variables ('p_diag%...') are required, 
    !   so they should have been updated before they enter this subroutine
    ! 
    ! - The following variables are potentially (directly) updated by nudging increments: 
    ! 
    !   Const. driving data         | Time-dependent driving data
    !   ----------------------------------------------------------
    !                 { * rho       | * rho            * press
    !   Thermodynamic {             |             or  
    !                 { * theta_v   | * theta_v        * temp
    !                               |
    !                   * vn        | * vn
    !                               |
    !                               | * qv
    !   ----------------------------------------------------------
    !
    ! - The nudging formula reads:
    !       X(t) = X'(t) + nudge_coeff_X(z) * [ X0(t) - X'(t) ],     (0)
    !                                         |_______________|
    !                                                 |
    !                                                = dX
    !   where X' denotes the value of the variable before the nudging step, 
    !   nudge_coeff_X is the variable-specific, height-dependent nudging coefficient, 
    !   and X0 is the value of the variable from the nudging driving data. 
    !   Note: the nudging update X - X' does not depend on the time step, 
    !   so it is not a typical state transition process, 
    !   but rather a state displacement or shift.
    !
    ! - In case of time-dependent driving data, there are two possibilities 
    !   for the provided thermodynamic nudging increments: 
    ! 
    !   Non-hydrostatic variables | Hydrostatic variables
    !   --------------------------------------------------
    !   * drho                    | * dpress
    !                             | 
    !   * dtheta_v                | * dtemp
    !   --------------------------------------------------
    !
    !   The increments drho and dtheta_v are multiplied with the nudging coefficient 
    !   and are directly added to the corresponding prognostic variables of ICON. 
    ! 
    !   The increments dpress and dtemp are transformed into drho and dtheta_v 
    !   using the linear map: 
    !                       drho  = A * dtemp_v + B * dpress,                        (1)
    !                    dtheta_v = C * dtemp_v + D * dpress,                        (2)    
    !   where the factors A, B, C and D are completely determined by the ICON-state 
    !   before the nudging. In order to find expressions for (1) and (2), 
    !   we proceed in the following way:
    !   First, given the definition of the virtual temperature 
    !   for moist (not cloudy) air: 
    !                      temp_v = temp * ( 1 + vtmpc1 * qv ),                      (3)
    !   
    !   we use its differential: 
    !          dtemp_v = dtemp * ( 1 + vtmpc1 * qv ) + temp * vtmpc1 * dqv,          (4)
    !   to compute the virtual temperature increment dtemp_v.
    !   Next, given the equation of state: 
    !                        rho = press / ( rd * temp_v ),                          (5)
    !   we use its differential: 
    !       drho = dpress / ( rd * temp_v ) - press * dtemp_v / ( rd * temp_v**2 )
    !            = [ dpress / rd - rho * dtemp_v ] / temp_v,                         (6)
    !   to compute the density increment drho. 
    !   And finally, given the definition of the virtual potential temperature: 
    !          theta_v = temp_v * ( p0ref / press )**(rd/cpd) = temp_v / exner,      (7)
    !   we use its differential: 
    !     dtheta_v = dtemp_v / exner - temp_v / exner**2 * ( dexner / dpress ) * dpress  
    !              = dtemp_v / exner 
    !              - ( rd / cpd ) * ( temp_v / exner ) * ( 1 / press ) * dpress  
    !              = [ dtemp_v - dpress / ( cpd * rho ) ] / exner,                   (8)
    !   where: 
    !        dexner / dpress = ( rd / cpd ) * ( press / p0ref )**(rd/cp-1) / p0ref
    !                        = ( rd / cpd ) * ( exner / press )                      (9)
    !   and (5) have been used.
    !
    !   After the nudging increments (multiplied with the nudging coefficient) 
    !   have been added to rho and theta_v, the Exner pressure is updated by
    !   the following equivalent to (5):
    !              exner = ( rd * rho * theta_v / p0ref )**(rd/cvd).                 (10)
    !
    ! - Please note that a nudging of the thermodynamic variables 
    !   based on the hydrostatic set of variables (pres and temp) 
    !   does not imply that the nudging increments themselves 
    !   are hydrostatically balanced!
    !........................................................................................

    ! Domain index
    jg = p_patch%id
    
    ! Although this query should have been likely covered before calling this subroutine, 
    ! we might feel better, if we make sure that:
    IF (jg == 1 .AND. nudging_config%ltype(indg_type%globn)) THEN

      ! Which variables should be nudged?
      ! Thermodynamic variables
      l_thermdyn = nudging_config%lvar(indg_var%thermdyn) 
      ! Horizontal wind
      l_vn       = nudging_config%lvar(indg_var%vn)  
      ! Water vapor 
      ! (includes 'ltransport=.true. implicitly)
      l_qv       = nudging_config%lvar(indg_var%qv)

      ! Thermodynamic variables: 
      ! should we nudge the hydrostatic pressure & temperature, 
      ! or the density and virtual potential temperature?
      l_hydrostatic = nudging_config%thermdyn_type == ithermdyn_type%hydrostatic

      ! Do we have to compute water vapor tendencies? 
      l_qv_tend = l_qv .OR. ( l_thermdyn .AND. l_hydrostatic )

      ! Start and end indices for vertical loops
      istart = nudging_config%ilev_start
      iend   = nudging_config%ilev_end

      ! Historically, the nudging coefficient was dependent on the number of substeps 
      ! 'tsrat = ndyn_substeps' which may change during runtime. To this reason, 
      ! it was necessary to re-compute the nudging coefficients each time step.
      ! As the dependency on tsrat no longer exists, we might think about computing the 
      ! nudging coefficients only once.
      !
      IF (l_thermdyn) THEN 
        ALLOCATE(nudge_coeff_thermdyn(istart:iend), STAT=istat)
        IF (istat /= SUCCESS) CALL finish (routine, 'Allocation of nudge_coeff_thermdyn failed!')
      ENDIF
      IF (l_vn) THEN 
        ALLOCATE(nudge_coeff_vn(istart:iend), STAT=istat)
        IF (istat /= SUCCESS) CALL finish (routine, 'Allocation of nudge_coeff_vn failed!')
      ENDIF
      IF (l_qv) THEN 
        ALLOCATE(nudge_coeff_qv(istart:iend), STAT=istat)
        IF (istat /= SUCCESS) CALL finish (routine, 'Allocation of nudge_coeff_qv failed!')
      ENDIF
      DO jk = istart, iend
        nudge_coeff = p_metrics%nudgecoeff_vert(jk)
        IF (l_thermdyn) nudge_coeff_thermdyn(jk) = nudging_config%max_nudge_coeff_thermdyn * nudge_coeff
        IF (l_vn)       nudge_coeff_vn(jk)       = nudging_config%max_nudge_coeff_vn       * nudge_coeff
        IF (l_qv)       nudge_coeff_qv(jk)       = nudging_config%max_nudge_coeff_qv       * nudge_coeff
      ENDDO  !jk

      ! Allocate the field for the water vapor increments
      IF (l_qv_tend) THEN 
        ALLOCATE(qv_tend(nproma,istart:iend,p_patch%nblks_c), STAT=istat)
        IF (istat /= SUCCESS) CALL finish (routine, 'Allocation of qv_tend failed!')
      ENDIF

!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)
 
      !...............................................................
      !               Nudging of cell-based variables
      !...............................................................
      
      ! The horizontal loop runs over the entire prognostic domain + halo cells, 
      ! in order to avoid a synchronization afterwards.
      ! (In case of the limited-area mode, the following setting should 
      ! guarantee that the lateral boundary interpolation zone and 
      ! its halo cells are excluded from nudging.)
      rl_start   = grf_bdywidth_c + 1
      rl_end     = min_rlcell
      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)
      
      !---------------------------------------------------------------
      !                  Get water vapor increments
      !---------------------------------------------------------------
      
      ! For efficiency reasons, we put the query for the nudging variables 
      ! outside the jb-loop, but have to pay for it with some code overhead. 
      ! 
      ! Water vapor increments are only required (and allocated!), 
      ! if we nudge water vapor and/or if we nudge the hydrostatic set 
      ! of thermodynamic variables.
      IF (l_qv_tend) THEN
        
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk
          
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
          
          DO jk = istart, iend
            DO jc = i_startidx, i_endidx
              
              qv_tend(jc,jk,jb) = wfac_old * p_latbc_old%qv(jc,jk,jb) + wfac_new * p_latbc_new%qv(jc,jk,jb) &
                &               - p_prog_rcf%tracer(jc,jk,jb,iqv)
              
              ! The actual nudging step of qv has to be postponed to after the nudging 
              ! of the thermodynamic variables, because in its computation, 
              ! we prefer to access the unnudged 'p_prog_rcf%tracer(jc,jk,jb,iqv)'.
              
            ENDDO  !jc
          ENDDO  !jk
        ENDDO  !jb
!$OMP END DO
        
      ENDIF  !IF (l_qv_tend)
      
      !---------------------------------------------------------------
      !              Nudging of thermodynamic variables
      !---------------------------------------------------------------
      
      IF (l_thermdyn .AND. l_hydrostatic) THEN
        
        !---------------------------------------------------------------
        !      CASE 1: Nudge hydrostatic pressure and temperature
        !---------------------------------------------------------------
        
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,pres_tend,temp_tend,tempv_tend,thv_tend,rho_tend) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk
          
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
          
          DO jk = istart, iend
            DO jc = i_startidx, i_endidx
              
              ! Given increments:
              ! (We nudge the hydrostatic pressure, instead of the non-hydrostatic pressure = rho * rd * tempv)
              pres_tend = wfac_old * p_latbc_old%pres(jc,jk,jb) + wfac_new * p_latbc_new%pres(jc,jk,jb) &
                &       - p_diag%pres(jc,jk,jb)
              temp_tend = wfac_old * p_latbc_old%temp(jc,jk,jb) + wfac_new * p_latbc_new%temp(jc,jk,jb) &
                &       - p_diag%temp(jc,jk,jb)
              
              ! Transform water vapor and temperature increments into a virtual temperature increment
              ! (see eq. (4) above)
              tempv_tend = ( 1._wp + vtmpc1 * p_prog_rcf%tracer(jc,jk,jb,iqv) ) * temp_tend & 
                &        + p_diag%temp(jc,jk,jb) * vtmpc1 * qv_tend(jc,jk,jb)
              
              ! Transform increments of virtual temperature and pressure into increments 
              ! of density and virtual potential temperature (see eqs. (6) and (8) above)
              rho_tend = ( rrd * pres_tend - p_prog%rho(jc,jk,jb) * tempv_tend ) / p_diag%tempv(jc,jk,jb)
              thv_tend = ( tempv_tend - rcpd * pres_tend / p_prog%rho(jc,jk,jb) ) / p_prog%exner(jc,jk,jb)
              
              ! Nudging step
              ! (Please note: we assume that the nudging update is well-behaved, 
              ! i.e. it will never lead to negative values for the density and the virtual potential temperature. 
              ! For reasons of computational efficiency, we have to refrain from countermeasures (e.g. ... = MAX(..., ...)).
              ! However, if only one of the two, rho or theta_v, would happen to become negative, 
              ! the subsequent computation of the Exner pressure would at least point us to this fact, 
              ! since the program would crash due to a negative argument of the LOG function.)
              p_prog%rho(jc,jk,jb)     = p_prog%rho(jc,jk,jb)     + nudge_coeff_thermdyn(jk) * rho_tend
              p_prog%theta_v(jc,jk,jb) = p_prog%theta_v(jc,jk,jb) + nudge_coeff_thermdyn(jk) * thv_tend
              
              ! Update Exner pressure
              p_prog%exner(jc,jk,jb) = EXP( rd_o_cvd * &
                & LOG( rd_o_p0ref * p_prog%rho(jc,jk,jb) * p_prog%theta_v(jc,jk,jb) ) )
              
            ENDDO  !jc
          ENDDO  !jk
        ENDDO  !jb
!$OMP END DO
        
      ELSEIF (l_thermdyn .AND. .NOT. l_hydrostatic) THEN
        
        !---------------------------------------------------------------
        !    CASE 2: Nudge density and virtual potential temperature 
        !---------------------------------------------------------------
        
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,thv_tend,rho_tend) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk
          
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
          
          DO jk = istart, iend
            DO jc = i_startidx, i_endidx
              
              thv_tend = wfac_old * p_latbc_old%theta_v(jc,jk,jb) + wfac_new * p_latbc_new%theta_v(jc,jk,jb) &
                &      - p_prog%theta_v(jc,jk,jb)
              rho_tend = wfac_old * p_latbc_old%rho(jc,jk,jb) + wfac_new * p_latbc_new%rho(jc,jk,jb)         &
                &      - p_prog%rho(jc,jk,jb)
              
              p_prog%rho(jc,jk,jb)     = p_prog%rho(jc,jk,jb)     + nudge_coeff_thermdyn(jk) * rho_tend
              p_prog%theta_v(jc,jk,jb) = p_prog%theta_v(jc,jk,jb) + nudge_coeff_thermdyn(jk) * thv_tend
              p_prog%exner(jc,jk,jb)   = EXP( rd_o_cvd * &
                & LOG( rd_o_p0ref * p_prog%rho(jc,jk,jb) * p_prog%theta_v(jc,jk,jb) ) )
              
            ENDDO  !jc
          ENDDO  !jk
        ENDDO  !jb
!$OMP END DO
        
      ENDIF  !What thermodynamic variables are to be nudged?
      
      !---------------------------------------------------------------
      !                    Nudging of water vapor
      !---------------------------------------------------------------
      
      IF (l_qv) THEN
        
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,qv_update) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk
          
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
          
          DO jk = istart, iend
            DO jc = i_startidx, i_endidx
              
              ! Suppress positive nudging tendencies in saturated (cloudy) regions,
              ! in order to avoid runaway effects
              qv_update = MERGE( nudge_coeff_qv(jk) * qv_tend(jc,jk,jb),        &  ! Yes 
                &                MERGE( nudge_coeff_qv(jk) * qv_tend(jc,jk,jb), &  ! {Yes                            }
                &                       0._wp,                                  &  ! {No                             } No
                &                       qv_tend(jc,jk,jb) < 0._wp ),            &  ! {Negative water vapor increment?}
                &                p_prog_rcf%tracer(jc,jk,jb,iqc) < eps_qc       )  ! Cell cloud-free?
              
              p_prog_rcf%tracer(jc,jk,jb,iqv) = p_prog_rcf%tracer(jc,jk,jb,iqv) + qv_update
              
            ENDDO  !jc
          ENDDO  !jk
        ENDDO  !jb
!$OMP END DO
        
      ENDIF  !IF (l_qv)
      
      !...............................................................
      !                Nudging of edge-based variables
      !...............................................................
      
      ! The horizontal loop runs over the entire prognostic domain + halo edges, 
      ! in order to avoid a synchronization afterwards
      rl_start   = grf_bdywidth_e + 1
      rl_end     = min_rledge
      i_startblk = p_patch%edges%start_block(rl_start)
      i_endblk   = p_patch%edges%end_block(rl_end)
      
      !---------------------------------------------------------------
      !                   Nudging of horizontal wind
      !---------------------------------------------------------------
      
      IF (l_vn) THEN
        
!$OMP DO PRIVATE(jk,je,jb,i_startidx,i_endidx,vn_tend) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk
          
          CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
          
          DO jk = istart, iend
            DO je = i_startidx, i_endidx
              
              vn_tend = wfac_old * p_latbc_old%vn(je,jk,jb) + wfac_new * p_latbc_new%vn(je,jk,jb) &
                &     - p_prog%vn(je,jk,jb)
              
              p_prog%vn(je,jk,jb) = p_prog%vn(je,jk,jb) + nudge_coeff_vn(jk) * vn_tend
              
            ENDDO  !je
          ENDDO  !jk
        ENDDO  !jb
!$OMP END DO
        
      ENDIF  !IF (l_vn)
      
!$OMP END PARALLEL
      
      ! Clean-up
      IF (ALLOCATED(nudge_coeff_thermdyn)) THEN 
        DEALLOCATE(nudge_coeff_thermdyn, STAT=istat)
        IF (istat /= SUCCESS) CALL finish (routine, 'Deallocation of nudge_coeff_thermdyn failed!')
      ENDIF
      IF (ALLOCATED(nudge_coeff_vn)) THEN 
        DEALLOCATE(nudge_coeff_vn, STAT=istat)
        IF (istat /= SUCCESS) CALL finish (routine, 'Deallocation of nudge_coeff_vn failed!')
      ENDIF
      IF (ALLOCATED(nudge_coeff_qv)) THEN 
        DEALLOCATE(nudge_coeff_qv, STAT=istat)
        IF (istat /= SUCCESS) CALL finish (routine, 'Deallocation of nudge_coeff_qv failed!')
      ENDIF
      IF (ALLOCATED(qv_tend)) THEN 
        DEALLOCATE(qv_tend, STAT=istat)
        IF (istat /= SUCCESS) CALL finish (routine, 'Deallocation of qv_tend failed!')
      ENDIF
      
    ENDIF  !IF (jg == 1 .AND. nudging_config%ltype(indg_type%globn))

  END SUBROUTINE global_nudging

END MODULE mo_nudging
