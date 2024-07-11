! array allocation + matrix setup routines+field initialization+io.
! Initialization and IO are just very simple versions (when coupled to
! FESOM, ice uses FESOM machinery) that allow one to run test cases.
!
! Matrices are needed (i) if VP solver is used, (ii) in advection
! routines and (iii) for partitioning
!
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, AWI
!
! Authors: Sergey Danilov, Qiang Wang
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
!
! This file has been modified for the use in ICON
! ---------------------------------------------------------------

module mo_ice_fem_init

  USE mo_ice_fem_types
  USE mo_ice_fem_mesh

  USE mo_kind,          ONLY: wp
  USE mo_run_config,    ONLY: dtime
  USE mo_sea_ice_nml,   ONLY: Tevp_inv, i_ice_advec
!  USE mo_ice_fem_evp_old,   ONLY: init_evp_solver_coeffs_old
  USE mo_ice_fem_evp,       ONLY: init_evp_solver_coeffs
  USE mo_ice_fem_advection, ONLY: fem_fct_ice_init

  IMPLICIT NONE

  PUBLIC :: ice_init_fem

  PRIVATE :: set_par_support
  PRIVATE :: array_setup_ice
  PRIVATE :: icestiff_matrix
  PRIVATE :: ice_mass_matrix_fill

CONTAINS

!=================================================================
subroutine ice_init_fem

  Tevp_inv=3.0_wp/dtime ! inverse of the relaxation time for EVP, Tevp_inv=3/dtime is default

  ! Setting up FEM mesh
  call set_fem_mesh
  ! vla: trivial FEM mesh partition, retained for historical reasons
  call set_par_support

  ! allocate and inialize FEM sea ice fields
  call array_setup_ice

  ! allocate and inialize for FCT advection scheme
  if (i_ice_advec == 1) then
    call fem_fct_ice_init
  end if

  ! inialize coefficients for EVP solver
!  call init_evp_solver_coeffs_old
  call init_evp_solver_coeffs

end subroutine ice_init_fem

!=======================================================================
subroutine set_par_support

  integer     n

  !======= NODES
  ! Owned nodes + external nodes which I need:
  ! note: parallelization is treated by ICON
  myDim_nod2D=nod2D
  eDim_nod2D=0
  allocate(myList_nod2D(myDim_nod2D+eDim_nod2D))
  myList_nod2D(1:myDim_nod2D) = (/ (n, n=1,myDim_nod2D) /)

  !======= ELEMENTS
  myDim_elem2D=elem2D
  allocate(myList_elem2D(myDim_elem2D))
  myList_elem2D(1:myDim_elem2D)= (/ (n, n=1,myDim_elem2D) /)


end subroutine set_par_support

!=======================================================================

subroutine array_setup_ice
!
! inializing sea ice model
!

integer   :: k

!---------------------------------------------------------------------------------------------------
! The Coriolis term is now calculated in ice_fem_grid_init (modified by Vladimir Lapin 11/08/2015)
!---------------------------------------------------------------------------------------------------
! allocate(coriolis_nod2D(nod2D))
! coriolis_nod2D=2._wp*omega*sin(coord_nod2D(2,:))

 allocate(col_pos(nod2d))

! Allocate memory for variables of ice model
 allocate(m_ice(nod2D), a_ice(nod2D), m_snow(nod2D))
 allocate(u_ice(nod2D), v_ice(nod2D))
 allocate(sigma11(elem2D), sigma12(elem2D), sigma22(elem2D))
 allocate(rhs_m(nod2D), rhs_a(nod2D), rhs_u(nod2D), rhs_v(nod2D))
 allocate(rhs_mis(nod2D))

! initialize
 u_ice  = 0._wp
 v_ice  = 0._wp
 rhs_m  = 0._wp
 rhs_mis= 0._wp
 rhs_a  = 0._wp
 rhs_u  = 0._wp
 rhs_v  = 0._wp
 sigma11= 0._wp
 sigma22= 0._wp
 sigma12= 0._wp

! Allocate memory used for coupling (Partly used for input, and partly
! for output of information)
 allocate(stress_atmice_x(nod2D), stress_atmice_y(nod2D))
 allocate(elevation(nod2D))         ! =ssh  of ocean
 allocate(u_w(nod2D), v_w(nod2D))   ! =uf and vf of ocean at surface nodes

! Sets the structure of the stiffness matrix
 call icestiff_matrix
! Fill in  the mass matrix
 call ice_mass_matrix_fill
! Create lumped mass matrix
 allocate(lmass_matrix(nod2D))
DO k=1,nod2D
 lmass_matrix(k)=sum(mass_matrix(icestiff%rowptr(k):icestiff%rowptr(k+1)-1))
END DO

end subroutine array_setup_ice
!==========================================================================

SUBROUTINE icestiff_matrix()
!
! Sets the structure of the stiffness matrix for BE ice advection/VP problems
!

INTEGER                       :: k


! a) We must compute the number of potentially non-zero entries
! b) Allocate space
! c) Fill in the column index array with proper values


! a)
icestiff%dim=nod2D
allocate(icestiff%rowptr(icestiff%dim+1))
icestiff%rowptr(1)=1
DO k=1,nod2D
  icestiff%rowptr(k+1)=icestiff%rowptr(k)+nghbr_nod2D(k)%nmb
END DO
icestiff%nza=icestiff%rowptr(nod2D+1)-1
! b)
allocate(icestiff%colind(icestiff%nza))
allocate(icestiff%values(icestiff%nza))
allocate(mass_matrix(icestiff%nza))
icestiff%values=0.0_wp

! c)
        ! ===== FIND colind ======
DO k=1,nod2D
  icestiff%colind(icestiff%rowptr(k):icestiff%rowptr(k+1)-1)= &
    nghbr_nod2D(k)%addresses
END DO
END SUBROUTINE icestiff_matrix

!=========================================================================
SUBROUTINE ice_mass_matrix_fill
! Computes entries of mass matrix used to speed up stress2rhs and in the
! implementation of the FCT scheme

INTEGER                       :: col, elem, elnodes(3), offset
INTEGER                       :: ipos, row, q, n, i

  mass_matrix =0.0_wp
  DO i=1, myDim_elem2D
     elem=myList_elem2D(i)
     elnodes=elem2D_nodes(:,elem)
     DO n=1,3
        row=elnodes(n)
        ! Global-to-local neighbourhood correspondence
        DO q=1,nghbr_nod2D(row)%nmb
           col_pos(nghbr_nod2D(row)%addresses(q))=q
        END DO
        offset=icestiff%rowptr(row)-1
        DO q=1,3
           col=elnodes(q)
           ipos=offset+col_pos(col)
           mass_matrix(ipos)=mass_matrix(ipos)+voltriangle(elem)/12.0_wp
           IF(q==n) THEN
             mass_matrix(ipos)=mass_matrix(ipos)+voltriangle(elem)/12.0_wp
           END IF
        END DO
     END DO
  END DO

END SUBROUTINE ice_mass_matrix_fill

end module mo_ice_fem_init
