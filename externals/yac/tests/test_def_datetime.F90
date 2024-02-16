!!
!! @file test_def_datetime.F90
!!
!! @copyright Copyright  (C)  2013 DKRZ, MPI-M
!!
!! @author Moritz Hanke <hanke@dkrz.de>
!!
!
! Keywords:
! Maintainer: Moritz Hanke <hanke@dkrz.de>
!             Rene Redler <rene.redler@mpimet.mpg.de>
! URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
!
! This file is part of YAC.
!
! Redistribution and use in source and binary forms, with or without
! modification, are  permitted provided that the following conditions are
! met:
!
! Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! Neither the name of the DKRZ GmbH nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!


#include "test_macros.inc"

PROGRAM test_def_datetime

  USE utest
  USE mo_yac_finterface

  IMPLICIT NONE

  INTEGER :: instance_id, comp_id

  CALL start_test("yac_def_datetime")

  ! test with default instance

  CALL yac_finit()
  CALL yac_fread_config_yaml('test_def_datetime.yaml')

  CALL yac_fdef_comp("dummy", comp_id)
  CALL yac_fsync_def()

  CALL test(yac_fget_start_datetime() == '2008-03-09T16:05:07')
  CALL test(yac_fget_end_datetime() ==   '2008-03-10T16:05:07')

  CALL yac_fcleanup()
  CALL yac_finit()

  CALL yac_fdef_datetime(start_datetime = '2008-03-09T00:00:00')

  CALL yac_fdef_comp("dummy", comp_id)
  CALL yac_fsync_def()

  CALL test(yac_fget_start_datetime() == '2008-03-09T00:00:00')

  CALL yac_fcleanup()
  CALL yac_finit()

  CALL yac_fdef_datetime(end_datetime = '2008-03-10T00:00:00')

  CALL yac_fdef_comp("dummy", comp_id)
  CALL yac_fsync_def()

  CALL test(yac_fget_end_datetime() ==   '2008-03-10T00:00:00')

  CALL yac_fcleanup()
  CALL yac_finit()

  CALL yac_fdef_datetime( &
    start_datetime = '2008-03-09T16:05:07', &
    end_datetime = '2008-03-10T16:05:07')

  CALL yac_fdef_comp("dummy", comp_id)
  CALL yac_fsync_def()

  CALL test(yac_fget_start_datetime() == '2008-03-09T16:05:07')
  CALL test(yac_fget_end_datetime() ==   '2008-03-10T16:05:07')

  CALL yac_fcleanup()

  ! test with explicit instance

  CALL yac_finit_instance(instance_id)
  CALL yac_fread_config_yaml(instance_id, 'test_def_datetime.yaml')
  CALL yac_fdef_comp(instance_id, "dummy", comp_id)

  CALL yac_fsync_def(instance_id)

  CALL test(yac_fget_start_datetime_instance(instance_id) == '2008-03-09T16:05:07')
  CALL test(yac_fget_end_datetime_instance(instance_id) == '2008-03-10T16:05:07')

  CALL yac_fcleanup_instance(instance_id)
  CALL yac_finit_instance(instance_id)

  CALL yac_fdef_datetime_instance( &
       instance_id, start_datetime = '2008-03-09T00:00:00')

  CALL yac_fdef_comp(instance_id, "dummy", comp_id)
  CALL yac_fsync_def(instance_id)

  CALL test(yac_fget_start_datetime_instance(instance_id) == '2008-03-09T00:00:00')

  CALL yac_fcleanup_instance(instance_id)
  CALL yac_finit_instance(instance_id)

  CALL yac_fdef_datetime_instance( &
       instance_id, end_datetime = '2008-03-10T00:00:00')

  CALL yac_fdef_comp(instance_id, "dummy", comp_id)
  CALL yac_fsync_def(instance_id)

  CALL test(yac_fget_end_datetime_instance(instance_id) == '2008-03-10T00:00:00')

  CALL yac_fcleanup_instance(instance_id)
  CALL yac_finit_instance(instance_id)

  CALL yac_fdef_datetime_instance(                        &
    instance_id,  start_datetime = '2008-03-09T16:05:07', &
    end_datetime = '2008-03-10T16:05:07')

  CALL yac_fdef_comp(instance_id, "dummy", comp_id)
  CALL yac_fsync_def(instance_id)

  CALL test(yac_fget_start_datetime_instance(instance_id) == '2008-03-09T16:05:07')
  CALL test(yac_fget_end_datetime_instance(instance_id) == '2008-03-10T16:05:07')

  CALL yac_fcleanup_instance(instance_id)

  CALL yac_ffinalize()

  CALL stop_test

  CALL exit_tests

END PROGRAM test_def_datetime

