! Interface for non-hydrostatic testcases,
! which require some kind of update during the time integration.
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

MODULE mo_nh_testcase_interface

  USE mo_kind,                   ONLY: wp
  USE mo_impl_constants,         ONLY: TRACER_ONLY
  USE mo_model_domain,           ONLY: t_patch
  USE mo_nonhydro_types,         ONLY: t_nh_state
  USE mo_intp_data_strc,         ONLY: t_int_state
  USE mo_dynamics_config,        ONLY: nnow, nnew
  USE mo_nonhydrostatic_config,  ONLY: itime_scheme
  USE mo_nh_testcases_nml,       ONLY: nh_test_name, rotate_axis_deg, &
    &                                  lcoupled_rho
  USE mo_nh_pa_test,             ONLY: set_nh_w_rho
  USE mo_nh_df_test,             ONLY: get_nh_df_velocity
  USE mo_integrate_density_pa,   ONLY: integrate_density_pa
  USE mo_nh_dcmip_hadley,        ONLY: set_nh_velocity_hadley
  USE mo_exception,              ONLY: finish

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nh_testcase_interface

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nh_testcase_interface'

CONTAINS

  !>
  !! Interface for non-hydrostatic testcases 
  !!
  !! Interface for non-hydrostatic testcases, 
  !! which require some kind of update during the time integration.
  !!
  SUBROUTINE nh_testcase_interface( dt_loc,                  &  !in
    &                               sim_time,                &  !in
    &                               p_patch,                 &  !in 
    &                               p_nh_state,              &  !inout
    &                               p_int_state,             &  !in
    &                               jstep_adv_marchuk_order  )  !in

    ! In/out variables
    REAL(wp),                   INTENT(IN)    :: dt_loc                   !< advective time step on this grid level
    REAL(wp),                   INTENT(IN)    :: sim_time                 !< elapsed simulation time on this grid level
    TYPE(t_patch),     TARGET,  INTENT(INOUT) :: p_patch                  !< grid/patch info
    TYPE(t_nh_state),  TARGET,  INTENT(INOUT) :: p_nh_state               !< prognostic and diagnostic variables etc.
    TYPE(t_int_state), TARGET,  INTENT(IN)    :: p_int_state              !< interpolation state
    INTEGER,                    INTENT(IN)    :: jstep_adv_marchuk_order  !< Marchuk order

    ! Local variables
    INTEGER :: jg

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':nh_testcase_interface'
    
    !--------------------------------------------------------------

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! NOTE: a query for 'ltestcase_update' encompasses the calls of this subroutine. 
    ! Please, update the determination of 'ltestcase_update' 
    ! in 'src/testcases/mo_nh_testcases: init_nh_testcase', 
    ! if you add another testcase here. Otherwise no update will take place!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Domain index
    jg = p_patch%id

    ! Which time stepping scheme has been selected?
    IF (itime_scheme == TRACER_ONLY) THEN

      !----------------------------------------------------
      !                  Pure advection
      !----------------------------------------------------

      SELECT CASE ( TRIM(nh_test_name) )
        
      CASE ('PA') ! Solid body rotation

#ifdef _OPENACC
        CALL finish (routine, 'Test PA - Solid body rotation: OpenACC version currently not implemented')
#endif

        ! Set time-variant vertical velocity
        CALL set_nh_w_rho( p_patch,                      &  !in
          &                p_nh_state%metrics,           &  !in
          &                jstep_adv_marchuk_order,      &  !in
          &                dt_loc,                       &  !in
          &                sim_time-dt_loc,              &  !in
          &                p_nh_state%prog(nnew(jg))%w,  &  !inout
          &                p_nh_state%diag%pres,         &  !inout
          &                p_nh_state%diag%rho_ic        )  !inout
        
      CASE ('DF1', 'DF2', 'DF3', 'DF4') ! deformational flow

#ifdef _OPENACC
       CALL finish (routine, 'Tests DF1, DF2, DF3, DF4 - deformational flow: OpenACC version currently not implemented')
#endif

        ! Get velocity field
       CALL get_nh_df_velocity( p_patch,                    &  !in
         &                      p_nh_state%prog(nnew(jg)),  &  !inout
         &                      nh_test_name,               &  !in
         &                      rotate_axis_deg,            &  !in
         &                      sim_time-dt_loc+dt_loc      )  !in
        
        
        ! Get mass flux and new \rho. The latter one is only computed,
        ! if the density equation is re-integrated.
        CALL integrate_density_pa( p_patch,                    &  !in
          &                        p_int_state,                &  !in
          &                        p_nh_state%prog(nnow(jg)),  &  !in
          &                        p_nh_state%prog(nnew(jg)),  &  !in
          &                        p_nh_state%metrics,         &  !in
          &                        p_nh_state%diag,            &  !inout
          &                        dt_loc,                     &  !in
          &                        jstep_adv_marchuk_order,    &  !in
          &                        lcoupled_rho                )  !in

        
        
      CASE ('DCMIP_PA_12', 'dcmip_pa_12')

!#ifdef _OPENACC
!        CALL finish (routine, 'Test DCMIP_PA_12 - Hadley-like meridional circulation: OpenACC version currently not implemented')
!#endif
              
        ! Get velocity field for the DCMIP Hadley-like meridional circulation test
        !
        CALL set_nh_velocity_hadley( p_patch,                    &  !in
          &                          p_nh_state%prog(nnew(jg)),  &  !inout
          &                          p_nh_state%diag,            &  !in
          &                          p_int_state,                &  !in
          &                          p_nh_state%metrics,         &  !in   
          &                          sim_time-dt_loc+dt_loc,     &  !in
          &                          lacc=.TRUE.                 )  !in

!#ifdef _OPENACC
!        CALL finish (routine, 'Test DCMIP_PA_12 - Hadley-like meridional circulation: integrate_density_pa - OpenACC version currently not implemented')
!#endif

        ! Get mass flux and updated density for the DCMIP Hadley-like
        ! meridional circulation test
        !
        CALL integrate_density_pa( p_patch,                    &  !in
          &                        p_int_state,                &  !in
          &                        p_nh_state%prog(nnow(jg)),  &  !in
          &                        p_nh_state%prog(nnew(jg)),  &  !in
          &                        p_nh_state%metrics,         &  !in
          &                        p_nh_state%diag,            &  !inout
          &                        dt_loc,                     &  !in
          &                        jstep_adv_marchuk_order,    &  !in
          &                        lcoupled_rho                )  !in

      END SELECT

    ENDIF  !itime_scheme

  END SUBROUTINE nh_testcase_interface

END MODULE mo_nh_testcase_interface
