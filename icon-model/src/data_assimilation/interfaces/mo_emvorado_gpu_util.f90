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

MODULE mo_emvorado_gpu_util

  USE mo_run_config,              ONLY: ntracer, iqc, iqi, iqr, iqs, iqg, iqh, iqv, &
       &                              iqnc, iqni, iqnr, iqns, iqng, iqnh, iqgl, iqhl
  USE mo_nonhydro_state,        ONLY: p_nh_state
  USE mo_fortran_tools,           ONLY: assert_acc_device_only

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: radar_d2h_hydrometeors, radar_d2h_model_variables

  CONTAINS 

  !============================================================================
  ! Subroutine for copying all model hydrometeors required for EMVORADO. Amount of 
  ! data copies depends on the flags set for EMVORADO. Copies all data referenced in
  ! get_model_hydrometeors in module radar_interface
  !============================================================================
  SUBROUTINE radar_d2h_hydrometeors(ntlev, idom, lacc)

    INTEGER, INTENT(IN)           :: ntlev, idom
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    CALL assert_acc_device_only("radar_d2h_hydrometeors", lacc)

    !$ACC UPDATE &
    !$ACC   HOST(p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqv:iqv)) &
    !$ACC   HOST(p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqc:iqc)) &
    !$ACC   HOST(p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqr:iqr)) &
    !$ACC   HOST(p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqi:iqi)) &
    !$ACC   HOST(p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqnh:iqnh)) &
    !$ACC   HOST(p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqs:iqs)) &
    !$ACC   ASYNC(1)

    IF (iqg > 0 .AND. iqg <= ntracer) THEN
      !$ACC UPDATE HOST(p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqg:iqg)) ASYNC(1)
    END IF
    
    IF (iqh > 0 .AND. iqh <= ntracer) THEN
      !$ACC UPDATE HOST(p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqh:iqh)) ASYNC(1)
    END IF
    IF (iqnc > 0 .AND. iqnc <= ntracer) THEN
      !$ACC UPDATE HOST(p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqnc:iqnc)) ASYNC(1)
    END IF
    IF (iqnr > 0 .AND. iqnr <= ntracer) THEN
      !$ACC UPDATE HOST(p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqnr:iqnr)) ASYNC(1)
    END IF
    IF (iqni > 0 .AND. iqni <= ntracer) THEN
      !$ACC UPDATE HOST(p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqni:iqni)) ASYNC(1)
    END IF
    IF (iqns > 0 .AND. iqns <= ntracer) THEN
      !$ACC UPDATE HOST(p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqns:iqns)) ASYNC(1)
    END IF
    IF (iqng > 0 .AND. iqng <= ntracer) THEN
      !$ACC UPDATE HOST(p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqng:iqng)) ASYNC(1)
    END IF
    IF (iqnh > 0 .AND. iqnh <= ntracer) THEN
      !$ACC UPDATE HOST(p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqnh:iqnh)) ASYNC(1)
    END IF
    IF (iqgl > 0 .AND. iqgl <= ntracer) THEN
      ! 2mom scheme with liquid water fraction of graupel qgl:
      !$ACC UPDATE HOST(p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqgl:iqgl)) ASYNC(1)
    END IF
    IF (iqhl > 0 .AND. iqhl <= ntracer) THEN
      ! 2mom scheme with liquid water fraction of hail qhl:
      !$ACC UPDATE HOST(p_nh_state(idom)%prog(ntlev)%tracer(:,:,:,iqhl:iqhl)) ASYNC(1)
    END IF

    !$ACC WAIT(1)

  END SUBROUTINE radar_d2h_hydrometeors

  !============================================================================
  ! Subroutine for copying all model variables required for EMVORADO. Amount of 
  ! data copies depends on the flags set for EMVORADO. Copies all data referenced
  ! in get_model_variables in the module radar_interface
  !============================================================================
  SUBROUTINE radar_d2h_model_variables(ntlev_dyn, idom, lacc)

    INTEGER, INTENT(IN)           :: ntlev_dyn, idom
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    CALL assert_acc_device_only("radar_d2h_model_variables", lacc)

    !$ACC UPDATE &
    !$ACC   HOST(p_nh_state(idom)%metrics%z_ifc) &
    !$ACC   HOST(p_nh_state(idom)%metrics%z_mc) &
    !$ACC   HOST(p_nh_state(idom)%diag%u) &
    !$ACC   HOST(p_nh_state(idom)%diag%v) &
    !$ACC   HOST(p_nh_state(idom)%diag%temp) &
    !$ACC   HOST(p_nh_state(idom)%diag%pres) &
    !$ACC   HOST(p_nh_state(idom)%prog(ntlev_dyn)%w) &
    !$ACC   HOST(p_nh_state(idom)%prog(ntlev_dyn)%rho) &
    !$ACC   ASYNC(1)

    !$ACC WAIT(1)

  END SUBROUTINE radar_d2h_model_variables

END MODULE mo_emvorado_gpu_util
