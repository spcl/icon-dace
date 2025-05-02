changequote(`{',`}')dnl
include({f03-redist-gen.m4})dnl
!>
!! @file fname
!! @brief convenience wrappers of xt_redist exchanges for Fortran data
!!
!! @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
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
#include "fc_feature_defs.inc"
MODULE modname
  USE xt_redist_base, ONLY: xt_redist, xt_redist_s_exchange1, &
       xt_redist_a_exchange1
  USE xt_requests, ONLY: xt_request
#ifdef HAVE_FC_IS_CONTIGUOUS
  USE xt_core, ONLY: xt_abort
#endif
ifelse(typedecl,LOGICAL,
{#ifdef HAVE_FC_LOGICAL_INTEROP
})dnl
  USE iso_c_binding, ONLY: c_ptr, c_loc
ifelse(typedecl,{LOGICAL},
{#else
  USE iso_c_binding, ONLY: c_ptr
#endif
})dnl
  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER :: filename = 'modname.f90'
ifdef({prologue},{include(prologue)})dnl
  INTERFACE xt_redist_s_exchange
interface_gen_nd(typetag,{xt_redist_s_exchange})dnl
  END INTERFACE xt_redist_s_exchange
  PUBLIC :: xt_redist_s_exchange
  INTERFACE xt_redist_a_exchange
interface_gen_nd(typetag,{xt_redist_a_exchange})dnl
  END INTERFACE xt_redist_a_exchange
  PUBLIC :: xt_redist_a_exchange
CONTAINS
xt_redist_s_exchange_impl_gen(typetag,typedecl)dnl
xt_redist_a_exchange_impl_gen(typetag,typedecl)dnl
END MODULE modname
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! mode: f90
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
