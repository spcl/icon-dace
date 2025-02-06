!> Contains process factory
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_jsb_process_factory

  USE mo_jsb_memory_class,         ONLY: t_jsb_memory

  USE mo_jsb_process_factory_core, ONLY: &
    & Create_process_memory_core => Create_process_memory, &
    & Create_process, max_no_of_processes

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Create_process_memory
  PUBLIC :: Create_process, max_no_of_processes

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_process_factory'

CONTAINS

  FUNCTION Create_process_memory(iproc) RESULT(return_ptr)

    INTEGER, INTENT(in) :: iproc
    CLASS(t_jsb_memory), POINTER   :: return_ptr

    return_ptr => Create_process_memory_core(iproc)
    !$ACC ENTER DATA CREATE(return_ptr)
    !$ACC ENTER DATA CREATE(return_ptr%vars)

  END FUNCTION Create_process_memory

END MODULE mo_jsb_process_factory
