! This module provides definition of the FE mesh and P1-P1 elements.
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, AWI
!
! Authors: Sergey Danilov, Qiang Wang
! See LICENSES/Third_party_license_terms/BSD-3-Clause.fesom_seaice.txt for license information
! SPDX-License-Identifier: BSD-3-Clause
!
! This file has been modified for the use in ICON
! ---------------------------------------------------------------

module mo_ice_fem_mesh
  !
  ! Arrays defined here are used to keep mesh information
  !
  USE mo_kind,                  ONLY: wp
  USE mo_math_constants,        ONLY: deg2rad
  USE mo_physical_constants,    ONLY: earth_radius

  IMPLICIT NONE

  PUBLIC :: set_fem_mesh

  PRIVATE :: mesh_scaling
  PRIVATE :: basisfunctions
  PRIVATE :: build_nghbr_arrays

  PUBLIC
  SAVE
  !
  INTEGER                                   :: nod2D, elem2D
  REAL(wp), ALLOCATABLE, DIMENSION(:,:)     :: coord_nod2D
  INTEGER, ALLOCATABLE, DIMENSION(:)        :: index_nod2D
  INTEGER, ALLOCATABLE, DIMENSION(:,:)      :: elem2D_nodes, nod2D_elems, elem2D_nghbrs
  REAL(wp), ALLOCATABLE, DIMENSION(:)       :: cos_elem2D, sin_elem2D, metrics_elem2D, coriolis_nod2D
  REAL(wp), ALLOCATABLE, DIMENSION(:,:)     :: bafux, bafuy
  REAL(wp), ALLOCATABLE, DIMENSION(:)       :: voltriangle
  ! Mesh partition (trivial, retained for historical reasons). Equal to myDim_nod2D=nod2D, myDim_elem2D=elem2D
  INTEGER                                   :: myDim_nod2D, eDim_nod2D, myDim_elem2D
  INTEGER, ALLOCATABLE, DIMENSION(:)        :: myList_nod2D, myList_elem2D
  REAL(wp), ALLOCATABLE, DIMENSION(:,:)     :: bafux_nod, bafuy_nod

  ! vladimir: indices for elements/nodes with actual ice
  ! allocated/deallocated each time-step
  INTEGER                                   :: si_elem2D, si_nod2D
  INTEGER, ALLOCATABLE, DIMENSION(:)        :: si_idx_elem, si_idx_nodes
  !
  INTEGER, ALLOCATABLE, DIMENSION(:)        :: col_pos

  type addresstype
     INTEGER                                    :: nmb
     INTEGER(KIND=4), DIMENSION(:), pointer     :: addresses
  end type addresstype

  ! used in mo_ice_fem_advection
  type(addresstype),ALLOCATABLE,DIMENSION(:)    :: nod_in_elem2D, nghbr_nod2D
  !

CONTAINS

!----------------------------------------------------------------------------
!
! Prepare FE mesh and P1-P1 elements.
!
subroutine set_fem_mesh

  call mesh_scaling
  call basisfunctions
  call build_nghbr_arrays

end subroutine set_fem_mesh

!----------------------------------------------------------------------------
!
! Scales coord_nod2D from deg to rad, calculates element-wise metrics arrays.
!
subroutine mesh_scaling

  INTEGER         :: i

  coord_nod2D(1,:)=coord_nod2D(1,:)*deg2rad
  coord_nod2D(2,:)=coord_nod2D(2,:)*deg2rad

  allocate(cos_elem2D(elem2D), sin_elem2D(elem2D), metrics_elem2D(elem2D))

  do i=1, elem2D
     cos_elem2D(i)=sum(cos(coord_nod2D(2,elem2D_nodes(:,i))))/3.0_wp
     sin_elem2D(i)=sum(sin(coord_nod2D(2,elem2D_nodes(:,i))))/3.0_wp
     metrics_elem2D(i)=sin_elem2D(i)/cos_elem2D(i)/earth_radius
  end do
  !
end subroutine mesh_scaling

!----------------------------------------------------------------------------
!
! Calculate basisfunctions of finite elements: bafux, bafuy, voltriangle
!
subroutine basisfunctions

  REAL(wp)                          :: DET2D, Vol2D
  REAL(wp), DIMENSION(2, 3)         :: derivative_stdbf
  REAL(wp), DIMENSION(2, 3)         :: derivative_loczeta
  REAL(wp), DIMENSION(2,2)          :: jacobian2D, jacobian2D_inv
  INTEGER                           :: elem,i

  Vol2D=1._wp/2._wp

  allocate(bafux(3,elem2d), bafuy(3,elem2d))
  allocate(voltriangle(elem2d))

  bafux    =0._wp
  bafuy    =0._wp
  voltriangle   =0._wp

  ! BASISFUNCTIONS ON STANDARD 2D ELEMENTS (stdbafu):
  ! stdbafu(1)=1-x-y;   stdbafu(2)=x;   stdbafu(3)=y
  ! DERIVATIVES OF STANDARD BASISFUNCTIONS
  derivative_stdbf      = 0._wp
  derivative_stdbf(:,1) =-1._wp
  derivative_stdbf(1,2) = 1._wp
  derivative_stdbf(2,3) = 1._wp

  ! BASISFUNCTIONS on LOCAL 2D ELEMENTs (bafux, bafuy):
  do elem=1,elem2d
     call local_element_def(elem, derivative_stdbf, &
          jacobian2D, jacobian2D_inv, DET2D,  derivative_loczeta)
     do i=1,3
        bafux(i,elem) = derivative_loczeta(1,i)
        bafuy(i,elem) = derivative_loczeta(2,i)
     enddo
     voltriangle(elem) = abs(DET2D) * Vol2D
  enddo

end subroutine basisfunctions
!
!----------------------------------------------------------------------------
!
subroutine local_element_def(element, derivative_stdbf, &
  &                     jacobian2D, jacobian2D_inv, DET, derivative_loczeta)
  !--------------------------------------------------------------------------
  !                     2D LOCAL ELEMENT DEFINITIONS
  !     USES:
  !       coord_nod2D(2,nod2D)
  !       elem2D_nodes(3, elem2D)
  !       derivative_stabafu_x_2D(2, 3)
  !
  !     OUTPUT:
  !       jacobian2D(2,2), jacobian2D_inv(2,2), DET
  !       derivative_locbafu_x_2D(2, 3)
  !--------------------------------------------------------------------------
  !
  LOGICAL, PARAMETER    :: cyclic_geometry=.true.           ! periodic geometry
  REAL(wp), PARAMETER   :: domain_length=360.0_wp*deg2rad   ! in rad

  INTEGER, INTENT(IN)                   :: element
  REAL(wp), DIMENSION(2,3), intent(IN)  :: derivative_stdbf
  REAL(wp), DIMENSION(2,2), intent(OUT) :: jacobian2D
  REAL(wp), DIMENSION(2,2), intent(OUT) :: jacobian2D_inv
  REAL(wp),                 intent(OUT) :: DET
  REAL(wp), DIMENSION(2,3), intent(OUT) :: derivative_loczeta
  !
  REAL(wp), DIMENSION(2, 3)             :: local_cart
  REAL(wp), DIMENSION(3, 2)             :: der_transp
  INTEGER                               :: i, node
  REAL(wp)                              :: meancos
  REAL(wp), DIMENSION(3,2)              :: trans_temp
  !
  meancos=cos_elem2D(element)
  do i=1, 3
     node=elem2D_nodes(i,element)
     ! cartesian coordinates
     local_cart(1,i)=coord_nod2D(1,node)
     local_cart(2,i)=earth_radius * coord_nod2D(2,node)
  end do
  !
  !  T R A N S F O R M A T I O N - M A T R I X   "jacobian"
  do i=1,2
     jacobian2D(:,i)= local_cart(:,i+1)-local_cart(:,1)
     if (cyclic_geometry) then
        if (jacobian2D(1,i)> domain_length/2.0_wp) jacobian2D(1,i)=jacobian2D(1,i)-domain_length
        if (jacobian2D(1,i)<-domain_length/2.0_wp) jacobian2D(1,i)=jacobian2D(1,i)+domain_length
     end if
  end do
  jacobian2D(1,:)=jacobian2D(1,:)*meancos*earth_radius
  !
  !  I N V E R S E   O F   jacobian
  call matrix_inverse_2x2(jacobian2D, jacobian2D_inv, DET)
  !
  trans_temp=transpose(derivative_stdbf)
  der_transp=matmul(trans_temp, jacobian2D_inv)
  derivative_loczeta=transpose(der_transp)
  !
end subroutine local_element_def
!
!----------------------------------------------------------------------------
!
subroutine  matrix_inverse_2x2 (A, AINV, DET)
  !
  REAL(wp), DIMENSION(2,2), intent(IN)  :: A
  REAL(wp), DIMENSION(2,2), intent(OUT) :: AINV
  REAL(wp), intent(OUT)                 :: DET
  !
  INTEGER                                   :: i,j
  !
  DET  = A(1,1)*A(2,2) - A(1,2)*A(2,1)
  if ( DET .eq. 0.0_wp )  then
     do j=1,2
        write(*,*) (A(i,j),i=1,2)
     end do
     stop 'SINGULAR 2X2 MATRIX'
  else
     AINV(1,1) =  A(2,2)/DET
     AINV(1,2) = -A(1,2)/DET
     AINV(2,1) = -A(2,1)/DET
     AINV(2,2) =  A(1,1)/DET
  endif
end subroutine matrix_inverse_2x2
!
!---------------------------------------------------------------------------
!
subroutine build_nghbr_arrays
  !
  INTEGER                            :: j, k, m, a, tr(3), flag, el,ml!, pr(6)
  INTEGER                            :: i!, eltype, dgnod(2), pos
  INTEGER, DIMENSION(100)            :: AUX=0
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ind, check
  !
  !--------------- 2D mesh-------------------------------
  ! Builds nod_in_elem2D
  allocate(ind(nod2D), check(elem2D))
  ind=0
  do j=1,elem2D
     tr=elem2D_nodes(:,j)
     ind(tr)=ind(tr)+1
  end do
  allocate(nod_in_elem2D(nod2D))

  nod_in_elem2D%nmb=ind(1:nod2D)
  do j=1,nod2D
     allocate(nod_in_elem2D(j)%addresses(ind(j)))

  end do
  ind=0
  do i=1,elem2D
     tr=elem2D_nodes(:,i)
     ind(tr)=ind(tr)+1
     do j=1, 3
        nod_in_elem2D(tr(j))%addresses(ind(tr(j)))=i
     end do
  end do
  ! The list of elements is ordered, and no sorting is needed
  !
  ! Builds nghbr_nod2D
  allocate(nghbr_nod2D(nod2D))

  check=0
  do j=1, nod2D
     flag=0
     do m=1,nod_in_elem2D(j)%nmb
        el=nod_in_elem2D(j)%addresses(m)
        do k=1, 3
           a=elem2D_nodes(k,el)
           if (check(a)==0) then
              check(a)=1
              flag=flag+1
              aux(flag)=a
           end if
        end do
     end do
     nghbr_nod2D(j)%nmb=flag
     allocate(nghbr_nod2D(j)%addresses(flag))

     ! We need to sort array aux(1:flag)
     do m=flag,1,-1
        ml=maxloc(aux(1:flag), 1)
        nghbr_nod2D(j)%addresses(m)=aux(ml)
        check(aux(ml))=0
        aux(ml)=-999
     end do
  end do
  !
  deallocate(ind, check)
end subroutine build_nghbr_arrays

end module mo_ice_fem_mesh
!
!----------------------------------------------------------------------------
