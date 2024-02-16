!>
!! @file solver_internal_multi.f90
!! @author Florian Wilhelm
!!
!! @copyright Copyright  (C)  2010  Florian Wilhelm <Florian.Wilhelm@kit.edu>
!!
!! @version 1.0
!
! Keywords: scales ppm solver internal
! Maintainer: Florian Wilhelm <Florian.Wilhelm@kit.edu>
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
!
#define filename 'solver_internal_multi.f90'
!
SUBROUTINE ABORT_UNLESS_NORMAL0(x, x_str)
  REAL(PREC), INTENT(IN) :: x
  CHARACTER(len=*), INTENT(IN) :: x_str ! string representation of x
  INTEGER :: rank
#ifdef USE_MPI
  INTEGER :: error
#endif

#ifdef USE_MPI
  CALL mpi_comm_rank(mpi_comm_world, rank, error)
#else
  rank = 0
#endif

  IF ( .NOT. is_normal(x) ) THEN
    CALL abort_ppm("Not normal number found in " // x_str // &
         " on rank " // int2str(rank) // ".", filename, __LINE__)
  ENDIF

END SUBROUTINE ABORT_UNLESS_NORMAL0

! Aborts if abnormal (NaN, Inf etc.) entries are present
SUBROUTINE ABORT_UNLESS_NORMAL1(x, x_str)
  REAL(PREC), INTENT(IN) :: x(:)
  CHARACTER(len=*), INTENT(IN) :: x_str ! string representation of x
  INTEGER :: i, ie, rank
#ifdef USE_MPI
  INTEGER :: error
#endif

  ie = SIZE(x,1)

#ifdef USE_MPI
  CALL mpi_comm_rank(mpi_comm_world, rank, error)
#else
  rank = 0
#endif

  DO i=1,ie
    IF ( .NOT. is_normal(x(i)) ) THEN
      CALL abort_ppm("Abnormal number found in " // x_str // &
           ", i=" // int2str(i) // &
           " on rank " // int2str(rank) // ".", &
           filename, __LINE__)
    ENDIF
  ENDDO

END SUBROUTINE ABORT_UNLESS_NORMAL1

! Aborts if abnormal (NaN, Inf etc.) entries are present
SUBROUTINE ABORT_UNLESS_NORMAL2(x, x_str)
  REAL(PREC), INTENT(IN) :: x(:,:)
  CHARACTER(len=*), INTENT(IN) :: x_str ! string representation of x
  INTEGER :: i, j, ie, je, rank
#ifdef USE_MPI
  INTEGER :: error
#endif

  ie = SIZE(x,1)
  je = SIZE(x,2)

#ifdef USE_MPI
  CALL mpi_comm_rank(mpi_comm_world, rank, error)
#else
  rank = 0
#endif

  DO j=1,je
    DO i=1,ie
      IF ( .NOT. is_normal(x(i,j)) ) THEN
        CALL abort_ppm("Not normal number found in " // x_str // &
             ", i=" // int2str(i) // ", j=" // int2str(j) // &
             " on rank " // int2str(rank) // ".", &
             filename, __LINE__)
      ENDIF
    ENDDO
  ENDDO

END SUBROUTINE ABORT_UNLESS_NORMAL2

! Aborts if abnormal (NaN, Inf etc.) entries are present
SUBROUTINE ABORT_UNLESS_NORMAL3(x, x_str)
  REAL(PREC), INTENT(IN) :: x(:,:,:)
  CHARACTER(len=*), INTENT(IN) :: x_str ! string representation of x
  INTEGER :: i, j, k, ie, je, ke, rank
#ifdef USE_MPI
  INTEGER :: error
#endif

  ie = SIZE(x,1)
  je = SIZE(x,2)
  ke = SIZE(x,3)

#ifdef USE_MPI
  CALL mpi_comm_rank(mpi_comm_world, rank, error)
#else
  rank = 0
#endif

  DO k=1,ke
    DO j=1,je
      DO i=1,ie
        IF ( .NOT. is_normal(x(i,j,k)) ) THEN
          CALL abort_ppm("Not normal number found in " // x_str &
               // ", i=" // int2str(i) // ", j=" // int2str(j) // &
               ", k=" // int2str(k) // " on rank " &
               // int2str(rank) // ".", &
               filename, __LINE__)
        ENDIF
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE ABORT_UNLESS_NORMAL3

! sets halos to zero
SUBROUTINE CLEAR_HALOS(x, ext_x)
  REAL(PREC), INTENT(INOUT) :: x(:,:)
  TYPE(extent), INTENT(IN) :: ext_x(:)
  INTEGER :: jb, ib, il, jl, ie, je

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = PREC

  ie = SIZE(x,1)
  je = SIZE(x,2)
  ib = extent_start(ext_x(1))
  jb = extent_start(ext_x(2))
  il = extent_end(ext_x(1))
  jl = extent_end(ext_x(2))

  x(1:ib-1,:) = 0.0_precs
  x(il+1:ie,:) = 0.0_precs
  x(ib:il,1:jb-1) = 0.0_precs
  x(ib:il,jl+1:je) = 0.0_precs

END SUBROUTINE CLEAR_HALOS

#undef filename
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
