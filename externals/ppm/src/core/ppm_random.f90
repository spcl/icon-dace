!! @file ppm_random.f90
!! @brief gather interfaces to random numbers
!
! Copyright  (C)  2011  Thomas Jahns <jahns@dkrz.de>
!
! @version 1.0
! Keywords:
! @author Thomas Jahns <jahns@dkrz.de>
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

MODULE ppm_random
  USE ppm_std_type_kinds, ONLY: dp, sp
  USE ppm_math_extensions, ONLY: m_pi_dp, m_pi_sp
  USE ppm_irand_internal, ONLY: irand => ppm_irand, irandp => ppm_irandp, &
       irandr => ppm_irandr, &
       irand8 => ppm_irand8, irandp8 => ppm_irandp8, irand8_min, irand8_max, &
       irandr8 => ppm_irandr8, &
       drand => ppm_drand, drandp => ppm_drandp, drandr => ppm_drandr, &
       frand => ppm_frand, frandp => ppm_frandp, frandr => ppm_frandr, &
       irand_min, irand_max, &
       initialize_irand, finalize_irand, &
       a_rand, a_randp, a_randr, a_rand_mt, a_randp_mt, &
       a_randr_mt
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: irand, lrand, drand, frand, irandp, drandp, frandp
  PUBLIC :: irandr, drandr, frandr
  PUBLIC :: initialize_irand, finalize_irand
  PUBLIC :: a_rand, a_randp, a_randr, a_rand_mt, a_randp_mt
  PUBLIC :: a_randr_mt, drand_normal, frand_normal
  PUBLIC :: irand_min, irand_max
  PUBLIC :: irand8, irand8_min, irand8_max, irandp8, irandr8
CONTAINS
  !> return evenly distributed random logical value
  FUNCTION lrand() RESULT(p)
    LOGICAL :: p
    p = IAND(irandp(), 1) == 1
  END FUNCTION lrand

  !> normal distribution random number, double precision
  !! @param mean mean value of distribution
  !! @param variance square of standard deviation of distribution
  !! @return random value from normal distribution
  FUNCTION drand_normal(mean, variance) RESULT(r)
    REAL(dp), INTENT(in) :: mean, variance
    REAL(dp) :: r
    REAL(dp) :: u, v

    DO WHILE (.TRUE.)
      u = drandp()
      IF (u > 0.0_dp) EXIT
    END DO
    DO WHILE (.TRUE.)
      v = drandp()
      IF (v > 0.0_dp) EXIT
    END DO
    r = mean + variance * SQRT(-2.0_dp * LOG(u)) * COS(2.0_dp * m_pi_dp * v)
  END FUNCTION drand_normal

  !> normal distribution random number, single precision
  !! @param mean mean value of distribution
  !! @param variance square of standard deviation of distribution
  !! @return random value from normal distribution
  FUNCTION frand_normal(mean, variance) RESULT(r)
    REAL(sp), INTENT(in) :: mean, variance
    REAL(sp) :: r
    REAL(sp) :: u, v

    DO WHILE (.TRUE.)
      u = frandp()
      IF (u > 0.0_sp) EXIT
    END DO
    v = 0.0_sp
    DO WHILE (.TRUE.)
      v = frandp()
      IF (v > 0.0_sp) EXIT
    END DO
    r = mean + variance * SQRT(-2.0_sp * LOG(u)) * COS(2.0_sp * m_pi_sp * v)
  END FUNCTION frand_normal

  !> Chi distribution random number, double precision
  !! @param k degrees of freedom
  !! @return random value from chi distribution
  FUNCTION drand_chi(k) RESULT(y)
    INTEGER, INTENT(in) :: k
    REAL(dp) :: y

    REAL(dp) :: x
    INTEGER :: i

    y = 0.0_dp
    DO i = 1, k
      x = drand_normal(0.0_dp, 1.0_dp)
      y = y + x**2
    END DO
    y = SQRT(y)
  END FUNCTION drand_chi

  !> Chi distribution random number, single precision
  !! @param k degrees of freedom
  !! @return random value from chi distribution
  FUNCTION frand_chi(k) RESULT(y)
    INTEGER, INTENT(in) :: k
    REAL(sp) :: y

    REAL(sp) :: x
    INTEGER :: i

    y = 0.0_sp
    DO i = 1, k
      x = frand_normal(0.0_sp, 1.0_sp)
      y = y + x**2
    END DO
    y = SQRT(y)
  END FUNCTION frand_chi

  !> Maxwell-Boltzmann distribution random number, double precision
  !! @param a scale parameter
  !! @return random value from Maxwell-Boltzmann distribution
  FUNCTION drand_maxwell(a) RESULT(y)
    REAL(dp), INTENT(in) :: a
    REAL(dp) :: y

    REAL(dp) :: x
    INTEGER :: i

    y = 0.0_dp
    DO i = 1, 3
      x = drand_normal(0.0_dp, a**2)
      y = y + x**2
    END DO
    y = SQRT(y)
  END FUNCTION drand_maxwell

  !> Maxwell-Boltzmann distribution random number, single precision
  !! @param a scale parameter
  !! @return random value from Maxwell-Boltzmann distribution
  FUNCTION frand_maxwell(a) RESULT(y)
    REAL(sp), INTENT(in) :: a
    REAL(sp) :: y

    REAL(sp) :: x
    INTEGER :: i

    y = 0.0_sp
    DO i = 1, 3
      x = frand_normal(0.0_sp, a**2)
      y = y + x**2
    END DO
    y = SQRT(y)
  END FUNCTION frand_maxwell

END MODULE ppm_random
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
