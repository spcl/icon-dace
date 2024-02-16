!>
!! @file xt_config_f.f90
!! @brief Fortran interface to yaxt configuration object
!!
!! @copyright Copyright  (C)  2020 Jörg Behrens <behrens@dkrz.de>
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
MODULE xt_config_f
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_null_ptr, c_int, &
       c_char, c_null_char
  USE xt_core, ONLY: xt_abort
  IMPLICIT NONE
  PRIVATE
  ! note: this type must not be extended to contain any other
  ! components, its memory pattern has to match void * exactly, which
  ! it does because of C constraints
  TYPE, BIND(C), PUBLIC :: xt_config
#ifndef __G95__
    PRIVATE
#endif
    TYPE(c_ptr) :: cptr = c_null_ptr
  END TYPE xt_config

  INTERFACE
    ! this function must not be implemented in Fortran because
    ! PGI 11.x chokes on that
    FUNCTION xt_config_f2c(config) BIND(c, name='xt_config_f2c') RESULT(p)
      IMPORT :: c_ptr, xt_config
      IMPLICIT NONE
      TYPE(xt_config), INTENT(in) :: config
      TYPE(c_ptr) :: p
    END FUNCTION xt_config_f2c
  END INTERFACE

  PUBLIC :: xt_config_new, xt_config_delete
  PUBLIC :: xt_config_f2c
  PUBLIC :: xt_exchanger_id_by_name
  PUBLIC :: xt_config_get_exchange_method, xt_config_set_exchange_method
  PUBLIC :: xt_config_get_idxvec_autoconvert_size, &
       xt_config_set_idxvec_autoconvert_size
  INTEGER, PUBLIC, PARAMETER :: &
       xt_exchanger_irecv_send = 0, &
       xt_exchanger_irecv_isend = 1, &
       xt_exchanger_irecv_isend_packed = 2, &
       xt_exchanger_mix_isend_irecv = 3, &
       xt_exchanger_neigh_alltoall = 4

  CHARACTER(len=*), PARAMETER :: filename = 'xt_config_f.f90'

CONTAINS

  FUNCTION xt_config_new() RESULT(config)
    TYPE(xt_config) :: config
    INTERFACE
      FUNCTION xt_config_new_c() RESULT(config) &
           BIND(c, name='xt_config_new')
        IMPORT :: c_ptr
        IMPLICIT NONE
        TYPE(c_ptr) :: config
      END FUNCTION xt_config_new_c
    END INTERFACE
    config%cptr = xt_config_new_c()
  END FUNCTION xt_config_new

  SUBROUTINE xt_config_delete(config)
    TYPE(xt_config), INTENT(in) :: config
    INTERFACE
      SUBROUTINE xt_config_delete_c(config) BIND(c, name='xt_config_delete')
        IMPORT :: c_ptr
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in) :: config
      END SUBROUTINE xt_config_delete_c
    END INTERFACE
    CALL xt_config_delete_c(config%cptr)
  END SUBROUTINE xt_config_delete

  SUBROUTINE xt_config_set_exchange_method(config, method)
    TYPE(xt_config), INTENT(inout) :: config
    INTEGER, INTENT(in) :: method
    INTEGER(c_int) :: method_c
    INTERFACE
      SUBROUTINE xt_config_set_exchange_method_c(config, method) &
           BIND(c, name='xt_config_set_exchange_method')
        IMPORT :: c_int, c_ptr
        TYPE(c_ptr), VALUE :: config
        INTEGER(c_int), VALUE :: method
      END SUBROUTINE xt_config_set_exchange_method_c
    END INTERFACE
    method_c = INT(method, c_int)
    CALL xt_config_set_exchange_method_c(config%cptr, method_c)
  END SUBROUTINE xt_config_set_exchange_method

  FUNCTION xt_config_get_exchange_method(config) RESULT(method)
    TYPE(xt_config), INTENT(in) :: config
    INTEGER :: method
    INTERFACE
      FUNCTION xt_config_get_exchange_method_c(config) RESULT(method) &
           BIND(c, name='xt_config_get_exchange_method')
        IMPORT :: c_int, c_ptr
        TYPE(c_ptr), VALUE :: config
        INTEGER(c_int) :: method
      END FUNCTION xt_config_get_exchange_method_c
    END INTERFACE
    method = INT(xt_config_get_exchange_method_c(config%cptr))
  END FUNCTION xt_config_get_exchange_method

  FUNCTION xt_exchanger_id_by_name(name) RESULT(exchanger_id)
    CHARACTER(len=*), INTENT(in) :: name
    INTEGER :: exchanger_id
    INTERFACE
      FUNCTION xt_exchanger_id_by_name_c(name) RESULT(exchanger_id) &
           BIND(c, name='xt_exchanger_id_by_name')
        IMPORT :: c_char, c_int
        CHARACTER(len=1, kind=c_char), INTENT(in) :: name(*)
        INTEGER(c_int) :: exchanger_id
      END FUNCTION xt_exchanger_id_by_name_c
    END INTERFACE
    INTEGER(c_int) :: c_id
    CHARACTER(len=1) :: name_c(LEN(name)+1)
    INTEGER :: i, nlen
    nlen = LEN(name)
    DO i = 1, nlen
      name_c(i) = name(i:i)
    END DO
    name_c(nlen+1) = c_null_char
    c_id = xt_exchanger_id_by_name_c(name_c)
    exchanger_id = INT(c_id)
  END FUNCTION xt_exchanger_id_by_name

  FUNCTION xt_config_get_idxvec_autoconvert_size(config) RESULT(cnvsize)
    TYPE(xt_config), INTENT(in) :: config
    INTEGER :: cnvsize
    INTERFACE
      FUNCTION xt_config_get_idxvec_autoconvert_size_c(config) RESULT(cnvsize) &
           BIND(c, name='xt_config_get_idxvec_autoconvert_size')
        IMPORT :: c_int, c_ptr
        TYPE(c_ptr), VALUE :: config
        INTEGER(c_int) :: cnvsize
      END FUNCTION xt_config_get_idxvec_autoconvert_size_c
    END INTERFACE
    cnvsize = INT(xt_config_get_idxvec_autoconvert_size_c(config%cptr))
  END FUNCTION xt_config_get_idxvec_autoconvert_size

  SUBROUTINE xt_config_set_idxvec_autoconvert_size(config, cnvsize)
    TYPE(xt_config), INTENT(inout) :: config
    INTEGER, INTENT(in) :: cnvsize
    INTEGER(c_int) :: cnvsize_c
    INTERFACE
      SUBROUTINE xt_config_set_idxvec_autoconvert_size_c(config, cnvsize) &
           BIND(c, name='xt_config_set_idxvec_autoconvert_size')
        IMPORT :: c_int, c_ptr
        TYPE(c_ptr), VALUE :: config
        INTEGER(c_int), VALUE :: cnvsize
      END SUBROUTINE xt_config_set_idxvec_autoconvert_size_c
    END INTERFACE
    IF (cnvsize > HUGE(1_c_int) .OR. cnvsize < 1) &
      CALL xt_abort("invalid conversion size", filename, __LINE__)

    cnvsize_c = INT(cnvsize, c_int)
    CALL xt_config_set_idxvec_autoconvert_size_c(config%cptr, cnvsize_c)
  END SUBROUTINE xt_config_set_idxvec_autoconvert_size

END MODULE xt_config_f
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
