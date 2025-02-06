!>
!! @file test_xmap_dist_dir_intercomm_parallel_f.f90
!!
!! @copyright Copyright  (C)  2013 Jörg Behrens <behrens@dkrz.de>
!!                                 Moritz Hanke <hanke@dkrz.de>
!!                                 Thomas Jahns <jahns@dkrz.de>
!!
!! @author Jörg Behrens <behrens@dkrz.de>
!!         Moritz Hanke <hanke@dkrz.de>
!!         Thomas Jahns <jahns@dkrz.de>
!!
!
! Keywords:
! Maintainer: Jörg Behrens <behrens@dkrz.de>
!             Moritz Hanke <hanke@dkrz.de>
!             Thomas Jahns <jahns@dkrz.de>
! URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
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
#include "fc_feature_defs.inc"
PROGRAM test_xmap_dist_dir_intercomm_parallel
  USE xt_xmap_abstract, ONLY: xt_xmap_dist_dir_new
  USE test_xmap_common_intercomm_parallel, ONLY: &
       xmap_intercomm_parallel_test_main
  IMPLICIT NONE
  INTERFACE
    FUNCTION xmap_dist_dir_intercomm_wrap(src_idxlist, dst_idxlist, &
         inter_comm) RESULT(xmap)
      USE yaxt, ONLY: xt_xmap, xt_idxlist
      TYPE(xt_idxlist), INTENT(in) :: src_idxlist, dst_idxlist
      INTEGER, INTENT(in) :: inter_comm
      TYPE(xt_xmap) :: xmap
    END FUNCTION  xmap_dist_dir_intercomm_wrap
  END INTERFACE
  CALL xmap_intercomm_parallel_test_main(xt_xmap_dist_dir_new, &
       call_finalize=.FALSE.)
  CALL xmap_intercomm_parallel_test_main(xmap_dist_dir_intercomm_wrap, &
       call_initialize=.FALSE.)
END PROGRAM test_xmap_dist_dir_intercomm_parallel

FUNCTION xmap_dist_dir_intercomm_wrap(src_idxlist, dst_idxlist, inter_comm) &
     RESULT(xmap)
  USE yaxt, ONLY: xt_xmap_dist_dir_intercomm_new, xt_xmap, xt_idxlist
  USE test_xmap_common_intercomm_parallel, ONLY: intra_group_comm
  TYPE(xt_idxlist), INTENT(in) :: src_idxlist, dst_idxlist
  INTEGER, INTENT(in) :: inter_comm
  TYPE(xt_xmap) :: xmap
  xmap = xt_xmap_dist_dir_intercomm_new(src_idxlist, dst_idxlist, &
       inter_comm, intra_group_comm)
END FUNCTION  xmap_dist_dir_intercomm_wrap

!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
