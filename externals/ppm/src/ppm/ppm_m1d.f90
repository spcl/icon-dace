!>
!! @file ppm_m1d.f90
!! @brief multi-level 1-dim partitioner
!!
!! @copyright Copyright  (C)  2012  Joerg Behrens <behrens@dkrz.de>
!!
!! @version 1.0
!! @author Joerg Behrens <behrens@dkrz.de>
!
! Maintainer: Joerg Behrens <behrens@dkrz.de>
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

MODULE ppm_m1d
  USE ppm_std_type_kinds, ONLY: dp
  USE ppm_base, ONLY: assertion




  IMPLICIT NONE
  PRIVATE


  TYPE :: xy_bounds_t
    ! rectangular region bounds
    ! end < start : zero size
    INTEGER :: xs=1 ! xstart
    INTEGER :: xe=0 ! xend
    INTEGER :: ys=1 ! ystart
    INTEGER :: ye=0 ! yend
  END TYPE xy_bounds_t
  TYPE(xy_bounds_t), PARAMETER:: zero_xy_bounds=xy_bounds_t(1, 0, 1, 0)


  CHARACTER(len=*),PARAMETER:: mcontext='ppm_m1d'

  PUBLIC:: xy_bounds_t, zero_xy_bounds

  PUBLIC :: d1_deco, d2_deco

  LOGICAL, PARAMETER:: debug=.FALSE.

CONTAINS

  !> private subroutine to calulate the solution of the 1-dim decomposition problem
  RECURSIVE SUBROUTINE vec_deco(wload, dstart, dend, nproc, refine)
    REAL(dp) :: wload(:) ! workload of grid points
    INTEGER, INTENT(out) :: dstart(:), dend(:) ! partitition borders
    INTEGER, INTENT(in) :: nproc ! number of procs
    LOGICAL, INTENT(in) :: refine

    REAL(dp) :: q(0:SIZE(wload))
    REAL(dp) :: avg, soll, d1, d2, pload(nproc)
    INTEGER  :: i, ix, nx, ip, ixs(nproc), ixe(nproc)

    nx=SIZE(wload)

    ! q =  integral of wload
    q(0)=0.0_dp
    DO i=1,nx
      q(i)=q(i-1)+wload(i)
    ENDDO

    ! gen slices with preferably uniform |q(i+1)-q(i)|
    ixs(:)=1
    ixe(:)=0
    avg=q(nx)/REAL(nproc, dp)
    ip=1
    ixs(1)=1
    soll=avg
    DO ix=1,nx
      IF (q(ix)>soll) THEN
        d1=ABS(soll-q(ix-1))
        d2=ABS(soll-q(ix))
        IF (d1<d2 .AND. ix>ixs(ip)) THEN
          ixe(ip)=ix-1
        ELSE
          ixe(ip)=ix
        ENDIF
        IF (ip==nproc) EXIT
        ip=ip+1
        ixs(ip)=ixe(ip-1)+1
        soll=REAL(ip, dp) * avg
      ENDIF
    ENDDO
    ixe(ip)=nx

    ! if we have zero occupation then we have an extreme peak load
    ! which we reconsider now:
    IF ( ixs(nproc) > ixe(nproc) ) THEN
      ! this is probably a bad solution but will do for now
      PRINT*,'WARNING: using simple workaround in vec_deco'
      ixs(nproc)=nx
      ixe(nproc)=nx
      CALL vec_deco(wload(1:nx-1), ixs(1:nproc-1), ixe(1:nproc-1), nproc-1, refine)
    ENDIF

    DO ip=1, nproc
      pload(ip)=0.0_dp
      DO ix=ixs(ip),ixe(ip)
        pload(ip)=pload(ip)+wload(ix)
      ENDDO
    ENDDO
    IF (debug) PRINT*,'pload error=', SUM(pload)-SUM(wload)
    IF (SIZE(dstart)<nproc) STOP 'vec_deco: bad size of start'

    DO ip=1,nproc
      dstart(ip)=ixs(ip)
      dend(ip)=ixe(ip)
      IF (dstart(ip)<1 .OR. dend(ip)>nx) THEN
        STOP 'vec_deco: bad bounds'
      ENDIF
    ENDDO
    IF (dend(nproc) /= nx) stop 'vec_deco: dend(nproc) /= nx'

    IF (refine) THEN
      CALL refine_c1d(wload, dend, pload)
      CALL gen_dstart(dstart=dstart, dend=dend)
    ENDIF

  END SUBROUTINE vec_deco

  !> private subroutine to calculate the start marks from end marks of connected intervals
  SUBROUTINE gen_dstart(dstart, dend)
    CHARACTER(len=*),PARAMETER:: context=mcontext//':gen_dstart'

    INTEGER, INTENT(in) :: dend(:)
    INTEGER, INTENT(out) :: dstart(:)

    INTEGER :: np, i

    np=SIZE(dend)
    CALL assertion(np == SIZE(dstart), msg=context//': np /= SIZE(dstart)')
    IF ( np == 0 ) RETURN

    dstart(1)=1
    DO i = 2, np
      dstart(i)=dend(i-1)+1
    ENDDO

  END SUBROUTINE gen_dstart

  !> Partially solve a two dimensional decomposition problem along one dimension only
  !!
  !! Along the selected dimension we get the optimal solution for elemental workloads
  !! in the limit of infinite resolution.
  !! @param wload workload for each element of the 2-dim domain
  !! @param dstart start indices of parts
  !! @param dend end indices of parts
  !! @param nproc number of intervals
  !! @param idim selected dimension
  !! @param refine try to improve the decomposition for finite resolution
  SUBROUTINE d1_deco(wload, dstart, dend, nproc, idim, refine)
    REAL(dp) :: wload(:,:) !wload(ix,iy)
    INTEGER, INTENT(inout) :: dstart(:), dend(:)
    INTEGER, INTENT(in) :: nproc, idim
    LOGICAL, INTENT(in) :: refine

    REAL(dp) :: xload(SIZE(wload,1)), yload(SIZE(wload,2))
    INTEGER :: ix,iy, nx, ny

    IF (SIZE(dstart) < nproc) STOP 'd1_deco: bad size of start'
    IF (idim<1 .OR. idim>2) STOP 'de_deco: bad idim'

    nx=SIZE(wload,1)
    ny=SIZE(wload,2)

    IF (idim==1) THEN

      DO ix=1,nx
        xload(ix)=0.0_dp
        DO iy=1,ny
          xload(ix)=xload(ix)+wload(ix,iy)
        ENDDO
      ENDDO
      IF (debug) PRINT*,'xload error=',SUM(xload)-SUM(wload)
      CALL vec_deco(xload, dstart, dend, nproc, refine)
      IF (dend(nproc) /= nx) STOP 'd1_deco: dend(nproc) /= nx'

    ELSE

      DO iy=1,ny
        yload(iy)=0.0_dp
        DO ix=1,nx
          yload(iy)=yload(iy)+wload(ix,iy)
        ENDDO
      ENDDO
      IF (debug) PRINT*,'yload error=',SUM(yload)-SUM(wload)
      CALL vec_deco(yload, dstart, dend, nproc, refine)
      IF (dend(nproc) /= ny) STOP 'd1_deco: dend(nproc) /= ny'

    ENDIF

  END SUBROUTINE d1_deco

  !> Solve a two dimensional decomposition problem
  !!
  !! Finds the optimal solution for elemental workloads
  !! in the limit of infinite resolution for the given process space
  !! @param wload workload for each element of the 2-dim domain
  !! @param pbounds bounds for each process, array shape == process space shape
  !! @param major_dim the dimension that is decomposed first
  !! @param refine try to improve the decomposition for finite resolution
  SUBROUTINE d2_deco(wload, pbounds, major_dim, refine)
    REAL(dp), INTENT(in) :: wload(:,:) !wload(ix,iy)
    TYPE(xy_bounds_t), INTENT(out) :: pbounds(:,:)
    INTEGER, INTENT(in) :: major_dim
    LOGICAL, INTENT(in) :: refine

    INTEGER :: nx, ny, nprocx, nprocy, xs, xe, ys, ye, ipx, ipy
    INTEGER :: pxs(SIZE(pbounds,1)),  pxe(SIZE(pbounds,1))
    INTEGER :: pys(SIZE(pbounds,2)),  pye(SIZE(pbounds,2))


    nprocx=SIZE(pbounds,1)
    nprocy=SIZE(pbounds,2)

    nx=SIZE(wload,1)
    ny=SIZE(wload,2)

    IF (major_dim == 1) THEN

      CALL d1_deco(wload, pxs, pxe, nprocx, 1, refine)
      IF (debug) PRINT*,'d2_deco: first d1_deco pxs=',pxs
      IF (debug) PRINT*,'d2_deco: first d1_deco pxe=',pxe
      DO ipx=1,nprocx
        xs=pxs(ipx)
        xe=pxe(ipx)
        IF (xs<1 .OR. xe>nx) STOP 'd2_deco: bad x bounds'
        IF (debug) PRINT*,'*** sub domain: ipx, xs, xe=',ipx, xs, xe
        CALL d1_deco(wload(xs:xe,:), pys, pye, nprocy, 2, refine)
        IF (debug) PRINT*,'d2_deco: second d1_deco pys=',pys
        IF (debug) PRINT*,'d2_deco: second d1_deco pye=',pye
        IF (pye(nprocy) /= ny) STOP 'd2_deco: pye(nprocy) /= ny'
        DO ipy=1,nprocy
          IF (pys(ipy)<1 .OR. pye(ipy)>ny) STOP 'd2_deco: bad y bounds'
          pbounds(ipx, ipy)%xs=pxs(ipx)
          pbounds(ipx, ipy)%xe=pxe(ipx)
          pbounds(ipx, ipy)%ys=pys(ipy)
          pbounds(ipx, ipy)%ye=pye(ipy)
        ENDDO
      END DO

    ELSEIF (major_dim == 2) THEN

      CALL d1_deco(wload, pys, pye, nprocy, 2, refine)
      IF (debug) PRINT*,'d2_deco: first d1_deco pys=',pys
      IF (debug) PRINT*,'d2_deco: first d1_deco pye=',pye
      DO ipy=1,nprocy
        ys=pys(ipy)
        ye=pye(ipy)
        IF (ys<1 .OR. ye>ny) STOP 'd2_deco: bad y bounds'
        IF (debug) PRINT*,'*** sub domain: ipx, ys, ye=',ipy, ys, ye
        CALL d1_deco(wload(:,ys:ye), pxs, pxe, nprocx, 1, refine)
        IF (debug) PRINT*,'d2_deco: second d1_deco pxs=',pxs
        IF (debug) PRINT*,'d2_deco: second d1_deco pxe=',pxe
        IF (pxe(nprocx) /= nx) STOP 'd2_deco: pxe(nprocx) /= nx'
        DO ipx=1,nprocx
          IF (pxs(ipx)<1 .OR. pxe(ipx)>nx) STOP 'd2_deco: bad x bounds'
          pbounds(ipx, ipy)%xs=pxs(ipx)
          pbounds(ipx, ipy)%xe=pxe(ipx)
          pbounds(ipx, ipy)%ys=pys(ipy)
          pbounds(ipx, ipy)%ye=pye(ipy)
        ENDDO
      END DO


    ELSE
      stop 'd2_deco: bad case'
    ENDIF

  END SUBROUTINE d2_deco

  !> refine a contiguous 1dim-decomposition
  !!
  !! @param wload workload for each element of the 1-dim domain
  !! @param em end marks of each part
  !! @param pload process load, must be valid or pload(1)<0
  !!        (in this case pload is recalculated from wload and em)
  !! @param nchange number of changes made by this subroutine
  SUBROUTINE refine_c1d(wload, em, pload, nchange)
    CHARACTER(len=*),PARAMETER:: context=mcontext//':refine_c1d'
    REAL(dp), INTENT(in) :: wload(:)
    INTEGER, INTENT(inout) :: em(:)
    REAL(dp), INTENT(inout) :: pload(:)
    INTEGER, OPTIONAL, INTENT(out) :: nchange

    REAL(dp) :: pmin, pmax, w
    INTEGER :: np, nc, ip, ipmin, ipmax

    np=size(pload)
    nc=0
    IF (np <= 0) RETURN ! nothing to do
    CALL assertion(np == SIZE(em), msg=context//': np /= SIZE(em)')

    IF (pload(1) < 0.0_dp) THEN
      pload(1) = SUM(wload(1:em(1)))
      DO ip=2,np
        pload(ip)=SUM(wload(em(ip-1)+1:em(ip)))
      ENDDO
    ENDIF

    refine: DO

      ! search max:
      ipmax=1
      pmax=pload(1)
      DO ip=2,np
        !PRINT*,'ip, pload=',ip, pload(ip)
        IF (pload(ip) > pmax) THEN
          pmax=pload(ip)
          ipmax=ip
        ENDIF
      ENDDO
      !DO ip=1,np
      !  PRINT*,'ip, delta=',ip, pload(ip)-pmax
      !ENDDO

      !PRINT*,'ipmax, pmax=',ipmax, pmax
      ! search right side for suitable hollow:

      ! can right neighbor proc carry my border weight?
      IF (ipmax<np) THEN
        w=wload(em(ipmax))
        !PRINT*,'w, nn cost=',w, pload(ipmax+1)+w - pmax
        ip=ipmax+1
        IF (pload(ip)+w < pmax) THEN
          !PRINT*,'w-success'
          em(ipmax)=em(ipmax)-1
          pload(ipmax)=pload(ipmax)-w
          pload(ip)=pload(ip)+w
          nc=nc+1
          CYCLE refine
        ENDIF
      ENDIF

      ! can left neighbor proc carry my border weight?
      IF (ipmax>1) THEN
        w=wload(em(ipmax-1)+1)
        !PRINT*,'w, nn cost=',w, pload(ipmax-1)+w - pmax
        ip=ipmax-1
        IF (pload(ip)+w < pmax) THEN
          !PRINT*,'max:w-success'
          em(ipmax-1)=em(ipmax-1)+1
          pload(ipmax)=pload(ipmax)-w
          pload(ip)=pload(ip)+w
          nc=nc+1
          CYCLE refine
        ENDIF
      ENDIF

      ! search min:
      ipmin=1
      pmin=pload(1)
      DO ip=2,np
        !PRINT*,'ip, pload=',ip, pload(ip)
        IF (pload(ip) < pmin) THEN
          pmin=pload(ip)
          ipmin=ip
        ENDIF
      ENDDO

      ! can I carry the border weight of my right neighbor?
      IF (ipmin<np) THEN
        w=wload(em(ipmin)+1)
        ip=ipmin+1
        IF (pload(ipmin)+w < pload(ip)) THEN
          !PRINT*,'min-success: right'
          em(ipmin)=em(ipmin)+1
          pload(ipmin)=pload(ipmin)+w
          pload(ip)=pload(ip)-w
          nc=nc+1
          CYCLE refine
        ENDIF
      ENDIF

      ! can I carry the border weight of my left neighbor?
      IF (ipmin>1) THEN
        w=wload(em(ipmin-1))
        ip=ipmin-1
        IF (pload(ipmin)+w < pload(ip)) THEN
          !PRINT*,'min-success: left'
          em(ipmin-1)=em(ipmin-1)-1
          pload(ipmin)=pload(ipmin)+w
          pload(ip)=pload(ip)-w
          nc=nc+1
          CYCLE refine
        ENDIF
      ENDIF

      EXIT refine
    ENDDO refine

    IF (PRESENT(nchange)) nchange=nc
  END SUBROUTINE refine_c1d

END MODULE ppm_m1d
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
!
