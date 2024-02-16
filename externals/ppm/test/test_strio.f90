!> @file test_strio.f90
!! @brief test string parsing routines
!!
!! @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
! Keywords:
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
#include "fc_feature_defs.inc"
PROGRAM test_strio
  USE ppm_base, ONLY: abort_ppm, ppm_default_comm
  USE ppm_extents, ONLY: iinterval
  USE ppm_strio, ONLY: get_address, sscana, fmt_elem, &
       arg_i4, arg_i8, arg_sp, arg_dp
  USE ppm_math_extensions, ONLY: initialize_math_extensions, de_g_sp, de_g_dp, &
       de_g_sp_width, de_g_dp_width
  USE ppm_std_type_kinds, ONLY: i4, i8, sp, dp
  USE ppm_random, ONLY: initialize_irand, irand, irand8, irandr, frand, drand
  USE ieee_arithmetic, ONLY: ieee_is_nan
#ifdef _OPENMP
  USE omp_lib, ONLY: omp_get_thread_num
#endif
  IMPLICIT NONE
  INTERFACE
    SUBROUTINE ppm_test_de_coherence(de_g_sp_width, de_g_dp_width)
      INTEGER, INTENT(in) :: de_g_sp_width, de_g_dp_width
    END SUBROUTINE ppm_test_de_coherence
  END INTERFACE
  INTEGER, PARAMETER :: num_vars = 6
  TYPE(fmt_elem) :: var_fmt(num_vars)
  INTEGER(i4) :: count, i, ierror, i_ref, j, j_ref, rand_seed
  INTEGER(i8) :: k, l, k_ref, l_ref
  REAL(sp) :: f, f_ref
  REAL(dp) :: d, d_ref
#ifdef HAVE_VOLATILE
  VOLATILE :: i, j, k, l, f, d
#endif
  CHARACTER(len=1024) :: str
  CHARACTER(len=*), PARAMETER :: &
       out_fmt = '(4(i0,1x),' // de_g_sp // ',1x,' // de_g_dp // ')'
  INTEGER :: tid, t
  TYPE(iinterval) :: erange_f, erange_d
  CHARACTER(len=*), PARAMETER :: filename = 'test_strio.f90'
#ifdef _OPENMP
  tid = omp_get_thread_num()
#else
  tid = 0
#endif
  CALL initialize_irand(ppm_default_comm, 0, rand_seed)
  PRINT '(2(a,i0))', 'thread id=', tid, ', random seed=', rand_seed
  CALL initialize_math_extensions
  CALL ppm_test_de_coherence(de_g_sp_width, de_g_dp_width)
  erange_f = iinterval(INT(MINEXPONENT(f), i4), INT(MAXEXPONENT(f), i4))
  erange_d = iinterval(INT(MINEXPONENT(d), i4), INT(MAXEXPONENT(d), i4))
  var_fmt(1)%argtype = arg_i4
  CALL get_address(i, var_fmt(1)%addr)
  var_fmt(2)%argtype = arg_i4
  CALL get_address(j, var_fmt(2)%addr)
  var_fmt(3)%argtype = arg_i8
  CALL get_address(k, var_fmt(3)%addr)
  var_fmt(4)%argtype = arg_i8
  CALL get_address(l, var_fmt(4)%addr)
  var_fmt(5)%argtype = arg_sp
  CALL get_address(f, var_fmt(5)%addr)
  var_fmt(6)%argtype = arg_dp
  CALL get_address(d, var_fmt(6)%addr)
  DO t = 1, 100000
    i_ref = irand()
    j_ref = irand()
    k_ref = irand8()
    l_ref = irand8()
#ifdef UNRELIABLE_DENORMAL
    ! we need to prevent creation of denormal numbers
    f_ref = SET_EXPONENT(frand(), irandr(erange_f))
    d_ref = SET_EXPONENT(drand(), irandr(erange_d))
#else
    f_ref = SCALE(frand(), irandr(erange_f))
    d_ref = SCALE(drand(), irandr(erange_d))
#endif
    WRITE (str, out_fmt) i_ref, j_ref, k_ref, l_ref, f_ref, d_ref
    CALL sscana(TRIM(str), var_fmt, count, ierror)
    IF (i /= i_ref .OR. j /= j_ref &
         .OR. k /= k_ref .OR. l_ref /= l &
#ifdef __GNUC__
         .OR. ((f > f_ref .OR. f < f_ref) .AND. .NOT. ieee_is_nan(f) &
         &     .AND. .NOT. ieee_is_nan(f_ref)) &
         .OR. ((d > d_ref .OR. d < d_ref) .AND. .NOT. ieee_is_nan(d) &
         &     .AND. .NOT. ieee_is_nan(d_ref)) &
#else
         .OR. (f /= f_ref .AND. .NOT. ieee_is_nan(f) &
         &     .AND. .NOT. ieee_is_nan(f_ref)) &
         .OR. (d /= d_ref .AND. .NOT. ieee_is_nan(d) &
         &     .AND. .NOT. ieee_is_nan(d_ref)) &
#endif
         .OR. count /= num_vars) THEN
      PRINT '(a,i0,2a)', 't=', t, ', str=', TRIM(str)
      IF (count /= num_vars) PRINT '(2(a,i0))', 'num_vars=', num_vars, &
           ', count=', count
      IF (i_ref /= i) PRINT '(2(a,i0))', 'i_ref=', i_ref, ', i=', i
      IF (j_ref /= j) PRINT '(2(a,i0))', 'j_ref=', j_ref, ', j=', j
      IF (k_ref /= k) PRINT '(2(a,i0))', 'k_ref=', k_ref, ', k=', k
      IF (l_ref /= l) PRINT '(2(a,i0))', 'l_ref=', l_ref, ', l=', l
#ifdef __GNUC__
      IF (f_ref > f .OR. f_ref < f) &
           PRINT '(2(a,'//de_g_sp//'))', 'f_ref=', f_ref, ', f=', f
      IF (d_ref > d .OR. d_ref < d) &
           PRINT '(2(a,'//de_g_dp//'))', 'd_ref=', d_ref, ', d=', d
#else
      IF (f_ref /= f) PRINT '(2(a,'//de_g_sp//'))', 'f_ref=', f_ref, ', f=', f
      IF (d_ref /= d) PRINT '(2(a,'//de_g_dp//'))', 'd_ref=', d_ref, ', d=', d
#endif
      CALL abort_ppm(source=filename, line=__LINE__, &
           msg='parsed values do not equal written values, &
           &or parse failed')
    END IF
  END DO
END PROGRAM test_strio
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! license-markup: "doxygen"
! End:
!
