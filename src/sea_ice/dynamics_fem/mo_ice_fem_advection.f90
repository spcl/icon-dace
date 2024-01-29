! Contains the FESIM implementation of FE-FCT advection scheme due to Loehner et al.
! Driver routine fct_ice_solve calls other routines that do low-order,
! higher order solutions and then combine them in a flux corrected way.
! Taylor-Galerkin type of advection is used.
!
! There is a tunable paremeter gamma in ts_solve_low_order and fem_fct.
! Increasing it leads to positivity preserving solution.
!
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_ice_fem_advection
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp

  USE mo_ice_fem_mesh
  USE mo_ice_fem_types
  USE mo_ice_fem_icon_init,         ONLY: exchange_nod2D

  IMPLICIT NONE

  PUBLIC  :: fct_ice_solve
  PUBLIC  :: ice_TG_rhs
  PUBLIC  :: fem_fct_ice_init

  PRIVATE :: ice_solve_high_order
  PRIVATE :: ice_solve_low_order
  PRIVATE :: ice_fem_fct

  PRIVATE
  ! FCT advection parameters:
  INTEGER :: num_iter_solve = 3 ! Iterations of mass matrix solver
  REAL(wp):: ice_gamma_fct=0.25 ! smoothing parameter in ice fct advection

CONTAINS

  !-------------------------------------------------------------------------
  !
  !> Driving routine
  !!
  subroutine fct_ice_solve
      implicit none

      call ice_solve_high_order   ! uses arrays of low-order solutions as temp
                                  ! storage. It should preceed the call of low
                      ! order solution.
      call ice_solve_low_order

      call ice_fem_fct(1)    ! m_ice
      call ice_fem_fct(2)    ! a_ice
      call ice_fem_fct(3)    ! m_snow
  end subroutine fct_ice_solve
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !> Computes the rhs in a Taylor-Galerkin way (with upwind type of
  !! correction for the advection operator)
  !!
  subroutine ice_TG_rhs

    USE mo_run_config,          ONLY: dtime

      implicit none
       real(wp)     :: diff, entries(3),  um, vm, vol, dx(3), dy(3)
       integer          :: n, q, row, elem, elnodes(3), m
       real(wp)     :: ice_diff

      DO m=1, myDim_nod2D
         row=myList_nod2D(m)
         rhs_m(row)=0.0_wp
         rhs_a(row)=0.0_wp
         rhs_mis(row)=0.0_wp ! rhs_mis in the original AWI code is called rhs_ms
      END DO
      ice_diff=0.0_wp                 ! No diffusion is  needed, but the option
                                   ! to add it is still preserved.
                                   !
      do m=1,myDim_elem2D          ! assembling rhs over elements
         elem=myList_elem2D(m)
         elnodes=elem2D_nodes(:,elem)
          !derivatives
         dx=bafux(:, elem)
         dy=bafuy(:, elem)
         vol=voltriangle(elem)
         um=sum(u_ice(elnodes))
         vm=sum(v_ice(elnodes))
          !diffusivity

         diff=ice_diff*max(1.0_wp,2*voltriangle(elem)/1.0e+8)
         DO n=1,3
            row=elnodes(n)
        DO q = 1,3
           entries(q)= vol*dtime*((dx(n)*(um+u_ice(elnodes(q)))+ &
                                dy(n)*(vm+v_ice(elnodes(q))))/12.0_wp - &
                           diff*(dx(n)*dx(q)+ dy(n)*dy(q))- &
                   0.5*dtime*(um*dx(n)+vm*dy(n))*(um*dx(q)+vm*dy(q))/9.0_wp)
            END DO
        rhs_m(row)=rhs_m(row)+sum(entries*m_ice(elnodes))
            rhs_a(row)=rhs_a(row)+sum(entries*a_ice(elnodes))
            rhs_mis(row)=rhs_mis(row)+sum(entries*m_snow(elnodes))
         END DO
      end do

  end subroutine ice_TG_rhs
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !> Initialization of arrays necessary to implement FCT algorithm
  !!
  subroutine fem_fct_ice_init
      
      implicit none

      allocate(m_icel(nod2D), a_icel(nod2D), m_snowl(nod2D))  ! low-order solutions
      allocate(dm_ice(nod2D), da_ice(nod2D), dm_snow(nod2D)) ! increments of high
                                                              ! order solutions
      allocate(icefluxes(myDim_elem2D,3))
      allocate(icepplus(nod2D), icepminus(nod2D))
      m_icel=0.0_wp
      a_icel=0.0_wp
      m_snowl=0.0_wp
  end subroutine fem_fct_ice_init
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !> Does Taylor-Galerkin solution the first approximation
  !!
  subroutine ice_solve_high_order
      implicit none
      !
      integer                           :: n,clo,clo2,cn,location(100),row,m
      real(wp)                          :: rhs_new

      do m=1,myDim_nod2D
         row=myList_nod2D(m)
         dm_ice(row)=rhs_m(row)/lmass_matrix(row)
         da_ice(row)=rhs_a(row)/lmass_matrix(row)
         dm_snow(row)=rhs_mis(row)/lmass_matrix(row)
      end do
      CALL exchange_nod2D(dm_ice, da_ice, dm_snow)
      !iterate
      do n=1,num_iter_solve-1
         do m=1,myDim_nod2D
            row=myList_nod2D(m)
            clo=icestiff%rowptr(row)
            clo2=icestiff%rowptr(row+1)-1
            cn=clo2-clo+1
            location(1:cn)=icestiff%colind(clo:clo2)
            rhs_new=rhs_m(row) - sum(mass_matrix(clo:clo2)*dm_ice(location(1:cn)))
            m_icel(row)=dm_ice(row)+rhs_new/lmass_matrix(row)
            rhs_new=rhs_a(row) - sum(mass_matrix(clo:clo2)*da_ice(location(1:cn)))
            a_icel(row)=da_ice(row)+rhs_new/lmass_matrix(row)
            rhs_new=rhs_mis(row) - sum(mass_matrix(clo:clo2)*dm_snow(location(1:cn)))
            m_snowl(row)=dm_snow(row)+rhs_new/lmass_matrix(row)
       end do
         do m=1,myDim_nod2D
            row=myList_nod2D(m)
            dm_ice(row)=m_icel(row)
            da_ice(row)=a_icel(row)
        dm_snow(row)=m_snowl(row)
         end do
         CALL exchange_nod2D(dm_ice, da_ice, dm_snow)
      end do

  end subroutine ice_solve_high_order
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !> Low-order solution
  !! It is assumed that m_ice, a_ice and m_snow from the previous time step
  !! are known at 1:myDim_nod2D+eDim_nod2D.
  !! One adds diffusive contribution to the rhs. It is realized as
  !! difference between the  consistent and lumped mass matrices
  !! acting on the field from the previous time step. The mass matrix on the
  !! lhs is replaced with lumped one.
  !!
  subroutine ice_solve_low_order

      implicit none
      integer       :: m, row, clo, clo2, cn, location(100)
      real(wp)      :: gamma

      gamma=ice_gamma_fct       ! Added diffusivity parameter
                                ! Adjust it to ensure posivity of solution

     do m=1,myDim_nod2D
            row=myList_nod2D(m)
            clo=icestiff%rowptr(row)
            clo2=icestiff%rowptr(row+1)-1
            cn=clo2-clo+1
            location(1:cn)=icestiff%colind(clo:clo2)
            m_icel(row)=(rhs_m(row)+gamma*sum(mass_matrix(clo:clo2)* &
                  m_ice(location(1:cn))))/lmass_matrix(row) + &
              (1.-gamma)*m_ice(row)
            a_icel(row)=(rhs_a(row)+gamma*sum(mass_matrix(clo:clo2)* &
                  a_ice(location(1:cn))))/lmass_matrix(row) + &
              (1.-gamma)*a_ice(row)
            m_snowl(row)=(rhs_mis(row)+gamma*sum(mass_matrix(clo:clo2)* &
                  m_snow(location(1:cn))))/lmass_matrix(row) + &
              (1.-gamma)*m_snow(row)
     end do
     CALL exchange_nod2D(m_icel, a_icel, m_snowl)
            ! Low-order solution must be known to neighbours

  end subroutine ice_solve_low_order
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !> Flux corrected transport algorithm for tracer advection
  !! It is based on Loehner et al. (Finite-element flux-corrected
  !! transport (FEM-FCT) for the Euler and Navier-Stokes equation,
  !! Int. J. Numer. Meth. Fluids, 7 (1987), 1093--1109) as described by Kuzmin and
  !! Turek. (kuzmin@math.uni-dortmund.de)
  !!
  subroutine ice_fem_fct(tr_array_id)
      implicit none

      integer   :: tr_array_id
      integer   :: icoef(3,3),m,n,q, elem,elnodes(3),row
      real(wp), dimension(myDim_nod2D) :: tmax, tmin
      real(wp)  :: vol, flux, ae, gamma


      gamma=ice_gamma_fct        ! It should coinside with gamma in
                                 ! ts_solve_low_order

     !==========================
     ! Compute elemental antidiffusive fluxes to nodes
     !==========================
     ! This is the most unpleasant part ---
     ! it takes memory and time. For every element
     ! we need its antidiffusive contribution to
     ! each of its 3 nodes

     ! Auxiliary elemental operator (mass matrix- lumped mass matrix)
       icoef=1
       do n=1,3   ! three upper nodes
      ! Cycle over rows  row=elnodes(n)
            icoef(n,n)=-2
       end do

     do m=1, myDim_elem2D
        elem=myList_elem2D(m)
        elnodes=elem2D_nodes(:,elem)
        vol=voltriangle(elem)
        if (tr_array_id==1) then
        do q=1,3
        icefluxes(m,q)=-sum(icoef(:,q)*(gamma*m_ice(elnodes) + &
                      dm_ice(elnodes)))*(vol/lmass_matrix(elnodes(q)))/12.0_wp
        end do
        end if

        if (tr_array_id==2) then
        do q=1,3
        icefluxes(m,q)=-sum(icoef(:,q)*(gamma*a_ice(elnodes) + &
                      da_ice(elnodes)))*(vol/lmass_matrix(elnodes(q)))/12.0_wp
        end do
        end if

        if (tr_array_id==3) then
        do q=1,3
        icefluxes(m,q)=-sum(icoef(:,q)*(gamma*m_snow(elnodes) + &
                      dm_snow(elnodes)))*(vol/lmass_matrix(elnodes(q)))/12.0_wp
        end do
        end if
     end do

     !==========================
     ! Screening the low-order solution
     !==========================
     ! TO BE ADDED IF FOUND NECESSARY
     ! Screening means comparing low-order solutions with the
     ! solution on the previous time step and using whichever
     ! is greater/smaller in computations of max/min below

     !==========================
     ! Cluster min/max
     !==========================
     if (tr_array_id==1) then
     do m=1, myDim_nod2D
        row=myList_nod2D(m)
             n=nghbr_nod2D(row)%nmb
             tmax(m)=maxval(m_icel(nghbr_nod2D(row)%addresses(1:n)))
         tmin(m)=minval(m_icel(nghbr_nod2D(row)%addresses(1:n)))
             ! Admissible increments
         tmax(m)=tmax(m)-m_icel(row)
         tmin(m)=tmin(m)-m_icel(row)
     end do
     end if

     if (tr_array_id==2) then
     do m=1, myDim_nod2D
        row=myList_nod2D(m)
             n=nghbr_nod2D(row)%nmb
             tmax(m)=maxval(a_icel(nghbr_nod2D(row)%addresses(1:n)))
         tmin(m)=minval(a_icel(nghbr_nod2D(row)%addresses(1:n)))
             ! Admissible increments
         tmax(m)=tmax(m)-a_icel(row)
         tmin(m)=tmin(m)-a_icel(row)
     end do
     end if

     if (tr_array_id==3) then
     do m=1, myDim_nod2D
        row=myList_nod2D(m)
             n=nghbr_nod2D(row)%nmb
             tmax(m)=maxval(m_snowl(nghbr_nod2D(row)%addresses(1:n)))
         tmin(m)=minval(m_snowl(nghbr_nod2D(row)%addresses(1:n)))
             ! Admissible increments
         tmax(m)=tmax(m)-m_snowl(row)
         tmin(m)=tmin(m)-m_snowl(row)
     end do
     end if

     !=========================
     ! Sums of positive/negative fluxes to node row
     !=========================

         icepplus=0.0_wp
         icepminus=0.0_wp
         do m=1, myDim_elem2D
            elem=myList_elem2D(m)
            elnodes=elem2D_nodes(:,elem)
            do q=1,3
               n=elnodes(q)
               flux=icefluxes(m,q)
               if (flux>0) then
                  icepplus(n)=icepplus(n)+flux
               else
                  icepminus(n)=icepminus(n)+flux
               end if
            end do
         end do

     !========================
     ! The least upper bound for the correction factors
     !========================
         do m=1,myDim_nod2D
         n=myList_nod2D(m)
         flux=icepplus(n)
         if (abs(flux)>0) then
         icepplus(n)=min(1.0_wp,tmax(m)/flux)
         else
         icepplus(n)=0.0_wp
         end if

         flux=icepminus(n)
         if (abs(flux)>0) then
         icepminus(n)=min(1.0_wp,tmin(m)/flux)
         else
         icepminus(n)=0.0_wp
         end if
         end do
      ! pminus and pplus are to be known to neighbouting PE
         CALL exchange_nod2D(icepminus, icepplus)
      !========================
      ! Limiting
      !========================
         do m=1, myDim_elem2D
            elem=myList_elem2D(m)
            elnodes=elem2D_nodes(:,elem)
            ae=1.0_wp
            do q=1,3
            n=elnodes(q)
            flux=icefluxes(m,q)
            if(flux>=0.0_wp) ae=min(ae,icepplus(n))
            if(flux<0.0_wp) ae=min(ae,icepminus(n))
            end do
            icefluxes(m,:)=ae*icefluxes(m,:)
            !if (ae.le.0.0_wp) write (*,*) 'ae is too large', ae
         end do

       !==========================
       ! Update the solution
       !==========================
             if(tr_array_id==1) then
         do m=1,myDim_nod2D
            n=myList_nod2D(m)
            m_ice(n)=m_icel(n)
         end do
         do m=1, myDim_elem2D
            elem=myList_elem2D(m)
            elnodes=elem2D_nodes(:,elem)
            do q=1,3
            n=elnodes(q)
            m_ice(n)=m_ice(n)+icefluxes(m,q)
            end do
         end do
             end if

         if(tr_array_id==2) then
         do m=1,myDim_nod2D
            n=myList_nod2D(m)
            a_ice(n)=a_icel(n)
         end do
         do m=1, myDim_elem2D
            elem=myList_elem2D(m)
            elnodes=elem2D_nodes(:,elem)
            do q=1,3
            n=elnodes(q)
            a_ice(n)=a_ice(n)+icefluxes(m,q)
            end do
         end do
             end if
         if(tr_array_id==3) then
         do m=1,myDim_nod2D
            n=myList_nod2D(m)
            m_snow(n)=m_snowl(n)
         end do
         do m=1, myDim_elem2D
            elem=myList_elem2D(m)
            elnodes=elem2D_nodes(:,elem)
            do q=1,3
            n=elnodes(q)
            m_snow(n)=m_snow(n)+icefluxes(m,q)
            end do
         end do
       END IF

       CALL exchange_nod2D(m_ice, a_ice, m_snow)
  end subroutine ice_fem_fct
  !-------------------------------------------------------------------------

END MODULE mo_ice_fem_advection
