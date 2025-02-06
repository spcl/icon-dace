! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

#include "test_macros.inc"

program test_def_comps

  use utest
  use yac
  implicit none

  integer :: comp_ids(2)

  character(len=YAC_MAX_CHARLEN) :: comp_names(2)

  logical  :: result

  call start_test("def_comps")

  call yac_finit ( )

  comp_names(1) = 'ICON-ocean'
  comp_names(2) = 'ICON-atmosphere'
  call yac_fdef_comps ( comp_names, 2, comp_ids )

  print *, ' def_comps returned comp_ids ', comp_ids

  result = ALL( comp_ids(:) /= -99 )

  call test ( result )

  call yac_ffinalize ( );

  call stop_test

  call exit_tests

end program test_def_comps
