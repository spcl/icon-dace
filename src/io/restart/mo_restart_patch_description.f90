! A CLASS that bundles all a patch's metadata together which IS relevant for restart purposes.
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

MODULE mo_restart_patch_description
  USE mo_communication, ONLY: t_comm_gather_pattern
  USE mo_dynamics_config, ONLY: nold, nnow, nnew, nnew_rcf, nnow_rcf
  USE mo_exception, ONLY: finish
  USE mo_fortran_tools, ONLY: assign_if_present_allocatable
  USE mo_cdi_constants, ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_EDGE
  USE mo_io_units, ONLY: filename_max
  USE mo_kind, ONLY: wp
  USE mo_model_domain, ONLY: p_patch, t_patch
  USE mo_mpi, ONLY: p_pe_work, my_process_is_work, process_mpi_all_comm, p_pe_work, process_mpi_root_id, &
                  & my_process_is_mpi_workroot, p_comm_work_2_restart
  USE mo_packed_message, ONLY: t_PackedMessage, kPackOp, kUnpackOp

  USE mo_upatmo_flowevent_utils, ONLY: t_upatmoRestartAttributes, upatmoRestartAttributesAssign, upatmoRestartAttributesPack

  USE mo_restart_util, ONLY: restartBcastRoot 

  IMPLICIT NONE
  PRIVATE

  TYPE optional_integer
    INTEGER :: v !< value if present
    LOGICAL :: present = .FALSE.
  END TYPE optional_integer

  TYPE t_gpat_ptr
    TYPE(t_comm_gather_pattern), POINTER :: p => NULL()
  END TYPE t_gpat_ptr

  PUBLIC :: t_restart_patch_description

  INTEGER, PARAMETER :: hmax = MAX(GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_EDGE)
  ! TYPE t_restart_patch_description contains all the DATA that
  ! describes a patch for restart purposes
  TYPE t_restart_patch_description
    ! logical patch id
    INTEGER :: id, nlev, cell_type, n_patch_elem_g(3), work_pe0_id, nold, &
      & nnow, nnew, nnew_rcf, nnow_rcf, opt_ndom = 1, opt_depth_lnd = -1
    INTEGER :: hmap(hmax)
    LOGICAL :: l_dom_active
    CHARACTER(LEN = filename_max) :: base_filename
    TYPE(optional_integer) :: opt_nlev_snow, opt_nice_class, opt_ndyn_substeps, &
         opt_jstep_adv_marchuk_order, opt_ocean_zlevels
    REAL(wp), ALLOCATABLE :: opt_pvct(:), opt_t_elapsed_phy(:), &
      & opt_ocean_zheight_cellMiddle(:), opt_ocean_zheight_cellInterfaces(:)

    TYPE(t_upatmoRestartAttributes) :: opt_upatmo_restart_atts

    TYPE(t_gpat_ptr) :: gpat(3)
  CONTAINS
    PROCEDURE :: init => restartPatchDescription_init
    ! called to set the DATA that may change from restart to restart
    PROCEDURE :: update => restartPatchDescription_update
    ! set the time level fields (nold, ...) to match the
    ! respective global variables
    PROCEDURE :: setTimeLevels => restartPatchDescription_setTimeLevels
    PROCEDURE :: updateOnMaster => restartPatchDescription_updateOnMaster
    PROCEDURE :: transferToRestart => restartPatchDescription_transferToRestart
    PROCEDURE :: packer => restartPatchDescription_packer
  END TYPE t_restart_patch_description

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_restart_patch_description"

CONTAINS

  SUBROUTINE restartPatchDescription_init(me, jg)
    CLASS(t_restart_patch_description), INTENT(INOUT) :: me
    INTEGER, INTENT(in) :: jg

    me%id = jg
    me%l_dom_active = .FALSE.
    me%hmap(:) = -1
    me%hmap(GRID_UNSTRUCTURED_CELL) = 1
    me%hmap(GRID_UNSTRUCTURED_VERT) = 2
    me%hmap(GRID_UNSTRUCTURED_EDGE) = 3
    CALL me%setTimeLevels()
    ! patch dependent info, p_patch IS NOT available on restart PEs
    IF(my_process_is_work()) THEN
      me%work_pe0_id = p_patch(jg)%proc0
      me%nlev = p_patch(jg)%nlev
      me%cell_type = p_patch(jg)%geometry_info%cell_type
      me%base_filename = TRIM(p_patch(jg)%grid_filename)
      me%n_patch_elem_g(1) = p_patch(jg)%n_patch_cells_g
      me%n_patch_elem_g(2) = p_patch(jg)%n_patch_verts_g
      me%n_patch_elem_g(3) = p_patch(jg)%n_patch_edges_g
      me%gpat(1)%p => p_patch(jg)%comm_pat_gather_c
      me%gpat(2)%p => p_patch(jg)%comm_pat_gather_v
      me%gpat(3)%p => p_patch(jg)%comm_pat_gather_e
    ELSE
      me%work_pe0_id = -1
      me%nlev = -1
      me%cell_type = -1
      me%base_filename = ''
      me%n_patch_elem_g(:) = -1
      NULLIFY(me%gpat(1)%p, me%gpat(2)%p, me%gpat(3)%p)
    END IF
  END SUBROUTINE restartPatchDescription_init

  SUBROUTINE restartPatchDescription_update(me, patch, opt_pvct, opt_t_elapsed_phy, &
                                           &opt_ndyn_substeps, opt_jstep_adv_marchuk_order, opt_depth_lnd, &
                                           &opt_nlev_snow, opt_nice_class, opt_ndom, & 

                                           &opt_upatmo_restart_atts,                                    &

                                           &opt_ocean_zlevels, &
                                           &opt_ocean_zheight_cellMiddle, opt_ocean_zheight_cellInterfaces )

    CLASS(t_restart_patch_description), INTENT(INOUT) :: me
    TYPE(t_patch), INTENT(IN) :: patch
    INTEGER, INTENT(IN), OPTIONAL :: opt_depth_lnd, opt_ndyn_substeps, opt_jstep_adv_marchuk_order, &
                                   & opt_nlev_snow, opt_nice_class, opt_ndom, opt_ocean_zlevels
    REAL(wp), INTENT(IN), OPTIONAL :: opt_pvct(:), opt_t_elapsed_phy(:), opt_ocean_zheight_cellMiddle(:), &
         & opt_ocean_zheight_cellInterfaces(:)

    TYPE(t_upatmoRestartAttributes), INTENT(IN), OPTIONAL :: opt_upatmo_restart_atts

    CHARACTER(*), PARAMETER :: routine = modname//":restartPatchDescription_update"

    IF(me%id /= patch%id) CALL finish(routine, "assertion failed: wrong patch passed to update()")
    ! update activity flag - this needs to be done on all compute
    ! PEs, otherwise starting a nest during runtime would lead to
    ! incomplete restart files
    me%l_dom_active = patch%ldom_active
    ! otherwise, only the patch master process needs the dynamic
    ! restart arguments, these will be communicated to the other
    ! processes anyway
    IF(p_pe_work == 0 .OR. p_pe_work == me%work_pe0_id) THEN
      ! Patch-dependent attributes
      CALL assign_if_present_allocatable(me%opt_pvct, opt_pvct)
      CALL assign_if_present_allocatable(me%opt_t_elapsed_phy, opt_t_elapsed_phy)
      CALL assign_if_present_allocatable(me%opt_ocean_zheight_cellMiddle, opt_ocean_zheight_cellMiddle)
      CALL assign_if_present_allocatable(me%opt_ocean_zheight_cellInterfaces, opt_ocean_zheight_cellInterfaces)
      CALL set_opt_int(me%opt_ndyn_substeps, opt_ndyn_substeps)
      CALL set_opt_int(me%opt_jstep_adv_marchuk_order, opt_jstep_adv_marchuk_order)
      IF (PRESENT(opt_depth_lnd)) me%opt_depth_lnd = opt_depth_lnd
      CALL set_opt_int(me%opt_nlev_snow, opt_nlev_snow)
      CALL set_opt_int(me%opt_nice_class, opt_nice_class)
      IF (PRESENT(opt_ndom)) me%opt_ndom = opt_ndom
      CALL set_opt_int(me%opt_ocean_zlevels, opt_ocean_zlevels)

      CALL upatmoRestartAttributesAssign(me%id, me%opt_upatmo_restart_atts, opt_upatmo_restart_atts)

      ! consistency check for OPTIONAL ocean variables
      IF(ALLOCATED(me%opt_ocean_zheight_cellMiddle)) THEN
        IF(.NOT. ALLOCATED(me%opt_ocean_Zheight_CellInterfaces) .OR. .NOT. me%opt_ocean_Zlevels%present) THEN
          CALL finish(routine, 'Ocean level parameteres not complete')
        END IF
      END IF
    END IF
  CONTAINS

    SUBROUTINE set_opt_int(o, opt_arg)
      TYPE(optional_integer), INTENT(out) :: o
      INTEGER, OPTIONAL, INTENT(in) :: opt_arg
      o%present = PRESENT(opt_arg)
      IF (o%present) THEN
        o%v = opt_arg
      ELSE
        o%v = -HUGE(o%v)
      END IF
    END SUBROUTINE set_opt_int
  END SUBROUTINE restartPatchDescription_update

  SUBROUTINE restartPatchDescription_setTimeLevels(me)
    CLASS(t_restart_patch_description), INTENT(INOUT) :: me

    me%nold = nold(me%id)
    me%nnow = nnow(me%id)
    me%nnow_rcf = nnow_rcf(me%id)
    me%nnew = nnew(me%id)
    me%nnew_rcf = nnew_rcf(me%id)
  END SUBROUTINE restartPatchDescription_setTimeLevels

  ! This ensures that the work master has complete up-to-date
  ! knowledge of the patch description.  To this END, this calls
  ! setTimeLevels() on the patchDescription, AND THEN sends the
  ! description from the subset master to the work master.
  SUBROUTINE restartPatchDescription_updateOnMaster(me)
    CLASS(t_restart_patch_description), INTENT(INOUT) :: me
  END SUBROUTINE restartPatchDescription_updateOnMaster

  SUBROUTINE restartPatchDescription_packer(me, op, pmsg)
    CLASS(t_restart_patch_description), INTENT(INOUT) :: me
    INTEGER, INTENT(in) :: op
    TYPE(t_PackedMessage), INTENT(INOUT) :: pmsg

    ! patch id AND activity flag
    CALL pmsg%packer(op, me%id)
    CALL pmsg%packer(op, me%l_dom_active)
    CALL pmsg%packer(op, me%nlev)
    CALL pmsg%packer(op, me%cell_type)
    CALL pmsg%packer(op, me%n_patch_elem_g(1))
    CALL pmsg%packer(op, me%n_patch_elem_g(2))
    CALL pmsg%packer(op, me%n_patch_elem_g(3))
    CALL pmsg%packer(op, me%work_pe0_id)
    CALL pmsg%packer(op, me%base_filename)
    ! time levels
    CALL pmsg%packer(op, me%nold)
    CALL pmsg%packer(op, me%nnow)
    CALL pmsg%packer(op, me%nnow_rcf)
    CALL pmsg%packer(op, me%nnew)
    CALL pmsg%packer(op, me%nnew_rcf)
    ! optional parameter values
    CALL pmsg%packer(op, me%opt_depth_lnd)
    CALL packer_opt_int(me%opt_nlev_snow)
    CALL packer_opt_int(me%opt_nice_class)
    CALL packer_opt_int(me%opt_ndyn_substeps)
    CALL packer_opt_int(me%opt_jstep_adv_marchuk_order)
    CALL pmsg%packer(op, me%opt_ndom)
    CALL packer_opt_int(me%opt_ocean_zlevels)
    ! optional parameter arrays
    CALL pmsg%packer(op, me%opt_pvct)
    CALL pmsg%packer(op, me%opt_t_elapsed_phy)
    CALL upatmoRestartAttributesPack(me%id, me%opt_upatmo_restart_atts, pmsg, op)
  CONTAINS

    SUBROUTINE packer_opt_int(oi)
      TYPE(optional_integer) :: oi

      CALL pmsg%packer(op, oi%present)
      IF (oi%present) CALL pmsg%packer(op, oi%v)
    END SUBROUTINE packer_opt_int
  END SUBROUTINE restartPatchDescription_packer

  SUBROUTINE restartPatchDescription_transferToRestart(me)
    CLASS(t_restart_patch_description), INTENT(INOUT) :: me
  END SUBROUTINE restartPatchDescription_transferToRestart

END MODULE mo_restart_patch_description
