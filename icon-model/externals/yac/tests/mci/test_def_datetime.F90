! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

#include "test_macros.inc"

PROGRAM test_def_datetime

  USE utest
  USE yac

  IMPLICIT NONE

  INTEGER :: instance_id, comp_id

  INTEGER :: arg_len
  INTEGER, PARAMETER :: max_opt_arg_len = 1024
  CHARACTER(max_opt_arg_len) :: config_dir

  CALL start_test("yac_def_datetime")

  CALL test(COMMAND_ARGUMENT_COUNT() == 1)
  CALL GET_COMMAND_ARGUMENT(1, config_dir, arg_len)

  ! test with default instance

  CALL yac_finit()
  CALL yac_fread_config_yaml(TRIM(config_dir) // 'test_def_datetime.yaml')

  CALL yac_fdef_comp("dummy", comp_id)
  CALL yac_fsync_def()

  CALL test(yac_fget_start_datetime() == '2008-03-09T16:05:07.000')
  CALL test(yac_fget_end_datetime() ==   '2008-03-10T16:05:07.000')

  CALL yac_fcleanup()
  CALL yac_finit()

  CALL yac_fdef_datetime(start_datetime = '2008-03-09T00:00:00.000')

  CALL yac_fdef_comp("dummy", comp_id)
  CALL yac_fsync_def()

  CALL test(yac_fget_start_datetime() == '2008-03-09T00:00:00.000')

  CALL yac_fcleanup()
  CALL yac_finit()

  CALL yac_fdef_datetime(end_datetime = '2008-03-10T00:00:00.000')

  CALL yac_fdef_comp("dummy", comp_id)
  CALL yac_fsync_def()

  CALL test(yac_fget_end_datetime() ==   '2008-03-10T00:00:00.000')

  CALL yac_fcleanup()
  CALL yac_finit()

  CALL yac_fdef_datetime( &
    start_datetime = '2008-03-09T16:05:07', &
    end_datetime = '2008-03-10T16:05:07')

  CALL yac_fdef_comp("dummy", comp_id)
  CALL yac_fsync_def()

  CALL test(yac_fget_start_datetime() == '2008-03-09T16:05:07.000')
  CALL test(yac_fget_end_datetime() ==   '2008-03-10T16:05:07.000')

  CALL yac_fcleanup()

  ! test with explicit instance

  CALL yac_finit_instance(instance_id)
  CALL yac_fread_config_yaml( &
    instance_id, TRIM(config_dir) // 'test_def_datetime.yaml')
  CALL yac_fdef_comp(instance_id, "dummy", comp_id)

  CALL yac_fsync_def(instance_id)

  CALL test(yac_fget_start_datetime_instance(instance_id) == '2008-03-09T16:05:07.000')
  CALL test(yac_fget_end_datetime_instance(instance_id) == '2008-03-10T16:05:07.000')

  CALL yac_fcleanup_instance(instance_id)
  CALL yac_finit_instance(instance_id)

  CALL yac_fdef_datetime_instance( &
       instance_id, start_datetime = '2008-03-09T00:00:00.000')

  CALL yac_fdef_comp(instance_id, "dummy", comp_id)
  CALL yac_fsync_def(instance_id)

  CALL test(yac_fget_start_datetime_instance(instance_id) == '2008-03-09T00:00:00.000')

  CALL yac_fcleanup_instance(instance_id)
  CALL yac_finit_instance(instance_id)

  CALL yac_fdef_datetime_instance( &
       instance_id, end_datetime = '2008-03-10T00:00:00.000')

  CALL yac_fdef_comp(instance_id, "dummy", comp_id)
  CALL yac_fsync_def(instance_id)

  CALL test(yac_fget_end_datetime_instance(instance_id) == '2008-03-10T00:00:00.000')

  CALL yac_fcleanup_instance(instance_id)
  CALL yac_finit_instance(instance_id)

  CALL yac_fdef_datetime_instance(                            &
    instance_id,  start_datetime = '2008-03-09T16:05:07.000', &
    end_datetime = '2008-03-10T16:05:07.000')

  CALL yac_fdef_comp(instance_id, "dummy", comp_id)
  CALL yac_fsync_def(instance_id)

  CALL test(yac_fget_start_datetime_instance(instance_id) == '2008-03-09T16:05:07.000')
  CALL test(yac_fget_end_datetime_instance(instance_id) == '2008-03-10T16:05:07.000')

  CALL yac_fcleanup_instance(instance_id)

  CALL yac_ffinalize()

  CALL stop_test

  CALL exit_tests

END PROGRAM test_def_datetime

