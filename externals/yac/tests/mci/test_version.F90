! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

#include "yac_config.h"

#include "test_macros.inc"

program test_version

  use utest
  use yac, only : yac_fget_version

  implicit none

  character(len=1024) :: macro_version
  character(len=:), allocatable :: function_version

  call start_test("yac_version")

  allocate(character(32) :: function_version)
  function_version = yac_fget_version()

  call test(allocated(function_version))

  if (allocated(function_version)) then
    write(macro_version, '(A,I0,A,I0,A,I0,A)') &
      "v", YAC_VERSION_MAJOR, ".", &
          YAC_VERSION_MINOR, ".", &
          YAC_VERSION_PATCH, &
          YAC_VERSION_TWEAK
    call test(TRIM(macro_version) == TRIM(function_version))

    write(macro_version, '(A,A)') "v", YAC_VERSION
    call test(TRIM(macro_version) == TRIM(function_version))

    print *, 'YAC Version: ' // TRIM(function_version)
  end if

  call stop_test

  call exit_tests

end program test_version

