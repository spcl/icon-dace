!> @file test_m1d.f90
!! @brief Test multi-level 1-dim partitioner
!!
!! @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>,
!!                                  Joerg Behrens <behrens@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
!! @author Joerg Behrens <behrens@dkrz.de>


!
! Keywords:
! Maintainer: Thomas Jahns <jahns@dkrz.de>, Joerg Behrens <behrens@dkrz.de>
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
PROGRAM test_m1d
  USE ppm_std_type_kinds, ONLY: dp
  USE ppm_base, ONLY: abort => abort_ppm, assertion
  USE ppm_extents, ONLY: iinterval, extent, ASSIGNMENT(=)
  USE ppm_uniform_partition, ONLY: uniform_partition_start
  USE ppm_m1d, ONLY: d2_deco, xy_bounds_t, zero_xy_bounds
  USE ppm_math_extensions, ONLY: m_pi_dp
  IMPLICIT NONE

  TYPE(xy_bounds_t),ALLOCATABLE :: pbounds(:,:)

  REAL(dp), ALLOCATABLE :: wload(:,:), pload(:,:), alpha(:,:)

  INTEGER, ALLOCATABLE :: xstart(:), ystart(:)

  REAL(dp) :: t_reg, e_reg, t_irreg, e_irreg

  CHARACTER(len=*), PARAMETER :: filename = 'test_m1d.f90'

  ALLOCATE(wload(800,400))
  CALL gen_workload(wload)


  CALL my_deco(regular=.TRUE., nprocx=32, nprocy=16, max_pload=t_reg, efficiency=e_reg)
  WRITE(0,*) '(test_m1d)   e_reg =', e_reg
  IF (e_reg<0.4_dp) CALL abort('unexpected bad regular efficiency',&
       & filename, __LINE__)

  CALL my_deco(regular=.FALSE., nprocx=32, nprocy=16, max_pload=t_irreg, efficiency=e_irreg)
  WRITE(0,*) '(test_m1d) e_irreg =',e_irreg
  IF (e_irreg<0.9_dp) CALL abort('unexpected bad irregular efficiency', &
       filename, __LINE__)


CONTAINS

  SUBROUTINE my_deco(regular, nprocx, nprocy, max_pload, efficiency)
    LOGICAL, INTENT(in) :: regular
    INTEGER, INTENT(in) :: nprocx, nprocy
    REAL(dp), INTENT(out) :: max_pload, efficiency

    IF (ALLOCATED(pload))  DEALLOCATE(pload)
    ALLOCATE(pload(nprocx, nprocy))

    IF (ALLOCATED(alpha))  DEALLOCATE(alpha)
    ALLOCATE(alpha(nprocx, nprocy))

    IF (ALLOCATED(xstart))  DEALLOCATE(xstart)
    ALLOCATE(xstart(nprocx))

    IF (ALLOCATED(ystart))  DEALLOCATE(ystart)
    ALLOCATE(ystart(nprocy))

    IF (ALLOCATED(pbounds))  DEALLOCATE(pbounds)
    ALLOCATE(pbounds(nprocx,nprocy))

    pbounds=zero_xy_bounds
    IF (regular) THEN
      CALL regular_d2_deco(wload, pbounds)
    ELSE
      CALL d2_deco(wload, pbounds, 1, refine=.FALSE.)
    ENDIF

    CALL pbounds2pload(wload, pbounds, pload)

    CALL pload2alpha(pload, alpha)

    efficiency=efficiency_d2(pload)

    max_pload=MAXVAL(pload)


  END SUBROUTINE my_deco

  SUBROUTINE regular_d2_deco(wload, pbounds)
    REAL(dp) :: wload(:,:)
    TYPE(xy_bounds_t), ALLOCATABLE :: pbounds(:,:)

    INTEGER :: i, j, m, n, npx, npy, &
         s_pos
    TYPE(extent) :: set_interval_x, set_interval_y

    npx = SIZE(pbounds, 1)
    npy = SIZE(pbounds, 2)
    set_interval_x = iinterval(1, SIZE(wload, 1))
    set_interval_y = iinterval(1, SIZE(wload, 2))
    m = SIZE(wload, 1)
    n = SIZE(wload, 2)
    pbounds(1, :)%xs = 1
    pbounds(:, 1)%ys = 1
    DO i = 2, npx
      s_pos = uniform_partition_start(set_interval_x, npx, i, .FALSE.)
      pbounds(i, :)%xs = s_pos
      pbounds(i - 1, :)%xe = s_pos - 1
    END DO
    pbounds(npx, :)%xe = m
    DO j = 2, npy
      s_pos = uniform_partition_start(set_interval_y, npy, j, .FALSE.)
      pbounds(:, j)%ys = s_pos
      pbounds(:, j - 1)%ye = s_pos - 1
    END DO
    pbounds(:, npy)%ye = n
  END SUBROUTINE regular_d2_deco


  CHARACTER(len=16) FUNCTION int2string(i)
    INTEGER, INTENT(in):: i

    CHARACTER(len=16)::ctmp

    WRITE(ctmp,*) i

    int2string=TRIM(ADJUSTL(ctmp))

  END FUNCTION int2string


  SUBROUTINE pload2alpha(pload, alpha)
    REAL(dp), INTENT(in)  :: pload(:,:) ! process space load
    REAL(dp), INTENT(out) :: alpha(:,:) ! load balance parameter alpha

    INTEGER  :: nx, ny, itx, ity
    REAL(dp) :: s, avg

    nx=SIZE(pload,1)
    ny=SIZE(pload,2)

    IF ( (SIZE(alpha,1) /= nx) .OR. (SIZE(alpha,2) /= ny) ) &
         & CALL abort('pload2alpha: size mismatch', filename, __LINE__)

    IF (nx*ny<1) CALL abort('zero process space', filename, __LINE__)

    s=0.0_dp
    DO ity=1,ny
      DO itx=1,nx
        s=s+pload(itx, ity)
      ENDDO
    ENDDO

    avg=s/REAL(nx*ny, dp)

    DO ity=1,ny
      DO itx=1,nx
        alpha(itx, ity)=(pload(itx, ity)-avg)/avg
      ENDDO
    ENDDO

  END SUBROUTINE pload2alpha


  SUBROUTINE pbounds2pload(wload, pb, pload)
    REAL(dp),INTENT(in)::wload(:,:)
    TYPE(xy_bounds_t),INTENT(in) :: pb(:,:)
    REAL(dp),INTENT(out)::pload(:,:)

    ! wload(igx, igy) = work load of grid space point (igx,igy)
    !
    ! pload(ipx, ipy) = summed load of all grid space points belonging
    !    to task WITH coords( ipx ,ipy) in process space
    !

    INTEGER :: igx, igy, ipx, ipy, nx, ny, nprocx, nprocy
    INTEGER :: xs, xe, ys, ye
    REAL(dp):: s

    nx=SIZE(wload,1)
    ny=SIZE(wload,2)

    nprocx=SIZE(pb,1)
    nprocy=SIZE(pb,2)


    DO ipy=1,nprocy
      DO ipx=1,nprocx

        xs=pb(ipx,ipy)%xs
        xe=pb(ipx,ipy)%xe
        ys=pb(ipx,ipy)%ys
        ye=pb(ipx,ipy)%ye
        IF (xs<1 .OR. xe>nx) CALL abort('gen_pload: bad x bounds', &
             filename, __LINE__)
        IF (ys<1 .OR. ye>ny) CALL abort('gen_pload: bad y bounds', &
             filename, __LINE__)

        s=0.0_dp
        DO igy=ys, ye
          DO igx=xs, xe
            s=s+wload(igx,igy)
          ENDDO
        ENDDO
        pload(ipx,ipy)=s

      ENDDO
    ENDDO

  END SUBROUTINE pbounds2pload


  REAL(dp) FUNCTION efficiency_d2(pload)
    REAL(dp),INTENT(in) :: pload(:,:)

    INTEGER  :: nx, ny, itx, ity
    REAL(dp) :: s, avg, alpha, e, max_pload

    ! estimates the quality of a given 2-dim. work distribution
    ! 1=best, 0=worst

    nx=SIZE(pload,1)
    ny=SIZE(pload,2)
    IF (nx<1 .OR. ny<1) THEN
      efficiency_d2=0.0_dp
      RETURN
    ENDIF

    max_pload=0.0_dp
    s=0.0_dp
    DO ity=1,ny
      DO itx=1,nx
        CALL assertion( pload(itx,ity) >= 0.0_dp , filename, __LINE__ , &
             'efficiency_d2: negative pload')
        max_pload=MAX(max_pload, pload(itx,ity))
        s=s+pload(itx, ity)
      ENDDO
    ENDDO

    avg=s/REAL(nx*ny, dp)
    alpha = ABS(max_pload-avg)/avg

    e=1.0_dp/(1.0_dp+alpha)

    CALL assertion( e>=0.0_dp .AND. e<=1.0_dp, &
         filename, __LINE__ ,'efficiency_d2: internal error')

    efficiency_d2=e

  END FUNCTION efficiency_d2

  SUBROUTINE gen_workload(w)
    REAL(dp), INTENT(out) :: w(:,:)

    REAL(dp), PARAMETER :: eps = 1.0e-10_dp
    INTEGER :: i,j, nx, ny
    REAL(dp) :: f0, f, ax, ay, s, qx, qy, qs

    nx = SIZE(w,1)
    ny = SIZE(w,2)

    qx = 1.0_dp/REAL(nx, dp)
    qy = 1.0_dp/REAL(ny, dp)

    f0 = 1.0_dp+eps
    s = 0.0_dp
    DO j=1,ny
      ay = 2.0_dp * m_pi_dp * REAL(j-1, dp) * qy + m_pi_dp
      DO i=1,nx
        ax = 2.0_dp * m_pi_dp * REAL(i-1, dp) * qx
        f =      f0 &
             & + 1.0_dp  * SIN(ax) * COS(ay) &
             & + 0.5_dp  * SIN(2.0_dp * ax) * COS(2.0_dp * ay) &
             & + 0.25_dp * SIN(3.0_dp * ax + 5.0_dp * ay)
        w(i,j) = abs(f)
        s = s + w(i,j)
      ENDDO
    ENDDO

    qs = 1.0_dp/s

    DO j=1,ny
      DO i=1,nx
        w(i,j) = w(i,j)*qs
      ENDDO
    ENDDO

  END SUBROUTINE gen_workload


END PROGRAM test_m1d
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
