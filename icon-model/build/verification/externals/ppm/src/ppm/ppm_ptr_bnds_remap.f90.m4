changequote(`{',`}')dnl
include({forloop2.m4})dnl
!> @file ppm_ptr_bnds_remap.f90
!! @brief work-around for compilers lacking remapping of pointer bounds
!!
!!
!! @copyright Copyright  (C)  2016  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
!!         Moritz Hanke <hanke@dkrz.de>
!
! Maintainer: Thomas Jahns <jahns@dkrz.de>
! URL: https://www.dkrz.de/redmine/projects/scales-ppm
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
define({interface_gen},
{forloop({dim},{1},{7},{{    MODULE PROCEDURE $2_$1_}dim{d}
})})
MODULE ppm_ptr_bnds_remap
  USE ppm_std_type_kinds, ONLY: i4, i8, sp, dp
  IMPLICIT NONE
  PRIVATE
  INTERFACE ptr_bnds_remap
interface_gen(i4, ptr_bnds_remap)dnl
interface_gen(i8, ptr_bnds_remap)dnl
interface_gen(l, ptr_bnds_remap)dnl
interface_gen(sp, ptr_bnds_remap)dnl
interface_gen(dp, ptr_bnds_remap)dnl
  END INTERFACE ptr_bnds_remap
  PUBLIC :: ptr_bnds_remap
CONTAINS
define({ptr_bnds_remap_impl},
{forloop({dim},{1},{7},
{{  ! replacement for missing pointer bounds remapping in some compilers
  SUBROUTINE ptr_bnds_remap_}$1{_}dim{d(ptr, pte, lb)
    $2, POINTER :: ptr(:}forloop({pdim},{2},dim,{,:}){)
    INTEGER, INTENT(in) :: lb(}dim{)
    $2, TARGET, INTENT(in) :: &
         pte(lb(1):}forloop({pdim},{2},dim,{,lb(pdim):}){)
    ptr => pte
  END SUBROUTINE ptr_bnds_remap_$1_}dim{d
}})})
ptr_bnds_remap_impl({i4},{INTEGER(i4)})dnl
ptr_bnds_remap_impl({i8},{INTEGER(i8)})dnl
ptr_bnds_remap_impl({l},{LOGICAL})dnl
ptr_bnds_remap_impl({sp},{REAL(sp)})dnl
ptr_bnds_remap_impl({dp},{REAL(dp)})dnl
END MODULE ppm_ptr_bnds_remap
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! mode: f90
! license-default: "bsd"
! license-markup: "doxygen"
! End:
!
