! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! Description:  Data type definition for initicon

MODULE mo_initicon_types

  USE mo_kind,                 ONLY: wp
  USE mo_impl_constants,       ONLY: max_ntracer, vname_len
  USE mo_run_config,           ONLY: ntracer
  USE mo_var,                  ONLY: t_var
  USE mo_dictionary,           ONLY: t_dictionary
  USE mo_ifs_coord,            ONLY: t_vct
  USE mo_fortran_tools,        ONLY: DO_DEALLOCATE, DO_PTR_DEALLOCATE
  USE mtime,                   ONLY: datetime


  IMPLICIT NONE
  PRIVATE


  !
  !variables
  PUBLIC :: t_init_state_const
  PUBLIC :: t_init_state      !> state vector for latbc data
  PUBLIC :: t_initicon_state  !> state vector for initicon
  PUBLIC :: t_pi_atm_in
  PUBLIC :: t_pi_sfc_in
  PUBLIC :: t_pi_atm
  PUBLIC :: t_pi_tracer
  PUBLIC :: t_pi_sfc
  PUBLIC :: t_sfc_inc
  PUBLIC :: geop_ml_var, alb_snow_var
  PUBLIC :: ana_varnames_dict
  ! functions
  PUBLIC :: t_init_state_finalize


 !
  TYPE :: t_pi_tracer

    REAL(wp), POINTER, DIMENSION(:,:,:) :: field => NULL()
    TYPE(t_var), POINTER :: var_element => NULL()
  CONTAINS
    PROCEDURE :: finalize => t_pi_tracer_finalize   !< destructor
  END TYPE t_pi_tracer


  ! atmospheric input variables
  TYPE :: t_pi_atm_in ! surface geopotential is regarded as
                      ! atmospheric variable here because the atmospheric fields cannot be processed without it

    ! Flag. True, if this data structure has been allocated
    LOGICAL :: linitialized

    ! vertical dimension of 3D input fields
    INTEGER :: nlev


    REAL(wp), POINTER, DIMENSION(:,:,:) :: temp    => NULL(), &
      &                                    pres    => NULL(), &
      &                                    u       => NULL(), &
      &                                    v       => NULL(), &
      &                                    w       => NULL(), &
      &                                    vn      => NULL(), &
      &                                    qv      => NULL(), &
      &                                    qc      => NULL(), &
      &                                    qi      => NULL(), &
      &                                    qr      => NULL(), &
      &                                    qs      => NULL(), &
      &                                    qg      => NULL(), &
      &                                    qh      => NULL(), &
      &                                    qnc     => NULL(), &
      &                                    qni     => NULL(), &
      &                                    qnr     => NULL(), &
      &                                    qns     => NULL(), &
      &                                    qng     => NULL(), &
      &                                    qnh     => NULL(), &
      &                                    rho     => NULL(), &
      &                                    theta_v => NULL(), &
      &                                    tke     => NULL()

    TYPE (t_pi_tracer), DIMENSION(max_ntracer) :: tracer

  CONTAINS
    PROCEDURE :: finalize => t_pi_atm_in_finalize   !< destructor
  END TYPE t_pi_atm_in


  ! vertical level height of input data
  TYPE :: t_init_state_const

    ! Flag. True, if this data structure has been allocated
    LOGICAL :: linitialized

    ! (half) level heights of input data
    REAL(wp), POINTER, DIMENSION(:,:,:) :: z_mc_in  => NULL()

    ! (half) level heights of the model
    REAL(wp), POINTER :: z_ifc(:,:,:), z_mc(:,:,:)

    REAL(wp), POINTER :: topography_c(:,:)

    ! vertical coordinate table for IFS coordinates
    TYPE(t_vct) :: vct

  CONTAINS
    PROCEDURE :: finalize => t_init_state_const_finalize   !< destructor
  END TYPE t_init_state_const


  ! surface input variables
  TYPE :: t_pi_sfc_in

    ! Flag. True, if this data structure has been allocated
    LOGICAL :: linitialized

    ! number of soil levels
    INTEGER :: nlevsoil

    REAL(wp), ALLOCATABLE, DIMENSION (:,:) :: tsnow, tskin, sst, snowalb,snowweq, snowdens, &
                                              skinres, ls_mask, seaice, phi
    REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: tsoil, wsoil

  CONTAINS
    PROCEDURE :: finalize => t_pi_sfc_in_finalize !< destructor
  END TYPE t_pi_sfc_in


  !
  TYPE :: t_pi_atm

    ! Flag. True, if this data structure has been allocated
    LOGICAL :: linitialized

    ! vertical dimension of 3D input fields
    INTEGER :: nlev

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: vn, u, v, w, temp, theta_v, exner, rho, &
                                               pres, qv, qc, qi, qr, qs, qg, tke, &
                                               qh, qnc, qni, qnr, qns, qng, qnh

    TYPE (t_pi_tracer), DIMENSION(max_ntracer) :: tracer

  CONTAINS
    PROCEDURE :: finalize => t_pi_atm_finalize   !< destructor
  END TYPE t_pi_atm


  !
  TYPE :: t_pi_sfc

    ! Flag. True, if this data structure has been allocated
    LOGICAL :: linitialized

    ! number of soil levels
    INTEGER :: nlevsoil

    REAL(wp), ALLOCATABLE, DIMENSION (:,:) :: tsnow, tskin, sst,  snowalb, snowweq, snowdens, &
                                              skinres, ls_mask, seaice
    REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: tsoil, wsoil

    REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: w_so

  CONTAINS
    PROCEDURE :: finalize => t_pi_sfc_finalize   !< destructor
  END TYPE t_pi_sfc


  ! surface field increments
  TYPE :: t_sfc_inc
    !
    ! Flag. True, if this data structure has been allocated
    LOGICAL :: linitialized

    ! number of soil levels
    INTEGER :: nlevsoil

    REAL(wp), ALLOCATABLE, DIMENSION (:,:,:) :: w_so
    REAL(wp), ALLOCATABLE, DIMENSION (:,:)   :: h_snow
    REAL(wp), ALLOCATABLE, DIMENSION (:,:)   :: freshsnow
    REAL(wp), ALLOCATABLE, DIMENSION (:,:)   :: t_2m

  CONTAINS
    PROCEDURE :: finalize => t_sfc_inc_finalize   !< destructor
  END TYPE t_sfc_inc


  ! state vector type: base class
  !
  TYPE :: t_init_state

    TYPE (t_pi_atm_in)     :: atm_in
    TYPE (t_pi_sfc_in)     :: sfc_in
    TYPE (t_pi_atm)        :: atm
    TYPE (t_pi_sfc)        :: sfc

    TYPE (t_init_state_const), POINTER :: const => NULL()

    TYPE (datetime)        :: vDateTime  ! validity DateTime
  CONTAINS
    PROCEDURE, PUBLIC :: finalize => t_init_state_finalize
  END TYPE t_init_state


  ! complete state vector type
  !
  TYPE, EXTENDS(t_init_state) :: t_initicon_state

    TYPE (t_pi_atm)        :: atm_inc
    TYPE (t_sfc_inc)       :: sfc_inc

  CONTAINS
    PROCEDURE, PUBLIC :: finalize => t_initicon_state_finalize
  END TYPE t_initicon_state


  CHARACTER(LEN=vname_len) :: geop_ml_var  ! model level surface geopotential
  CHARACTER(LEN=vname_len) :: alb_snow_var ! snow albedo


  ! dictionary which maps internal variable names onto
  ! GRIB2 shortnames or NetCDF var names.
  TYPE (t_dictionary) :: ana_varnames_dict


CONTAINS

  ! FINALIZE ROUTINES

  SUBROUTINE t_pi_tracer_finalize(tracer)
    CLASS(t_pi_tracer), INTENT(INOUT) :: tracer

    CALL DO_PTR_DEALLOCATE(tracer%field)
  END SUBROUTINE t_pi_tracer_finalize


  SUBROUTINE t_pi_atm_in_finalize(atm_in)
    CLASS(t_pi_atm_in), INTENT(INOUT) :: atm_in

    INTEGER :: idx

    atm_in%linitialized = .FALSE.
    CALL DO_PTR_DEALLOCATE(atm_in%temp)
    CALL DO_PTR_DEALLOCATE(atm_in%pres)
    CALL DO_PTR_DEALLOCATE(atm_in%u)
    CALL DO_PTR_DEALLOCATE(atm_in%v)
    CALL DO_PTR_DEALLOCATE(atm_in%w)
    CALL DO_PTR_DEALLOCATE(atm_in%vn)
    CALL DO_PTR_DEALLOCATE(atm_in%qv)
    CALL DO_PTR_DEALLOCATE(atm_in%qc)
    CALL DO_PTR_DEALLOCATE(atm_in%qi)
    CALL DO_PTR_DEALLOCATE(atm_in%qr)
    CALL DO_PTR_DEALLOCATE(atm_in%qs)
    CALL DO_PTR_DEALLOCATE(atm_in%qg)
    CALL DO_PTR_DEALLOCATE(atm_in%qh)
    CALL DO_PTR_DEALLOCATE(atm_in%qnc)
    CALL DO_PTR_DEALLOCATE(atm_in%qni)
    CALL DO_PTR_DEALLOCATE(atm_in%qnr)
    CALL DO_PTR_DEALLOCATE(atm_in%qns)
    CALL DO_PTR_DEALLOCATE(atm_in%qng)
    CALL DO_PTR_DEALLOCATE(atm_in%qnh)
    CALL DO_PTR_DEALLOCATE(atm_in%rho)
    CALL DO_PTR_DEALLOCATE(atm_in%theta_v)
    CALL DO_PTR_DEALLOCATE(atm_in%tke)

    DO idx=1, ntracer
      CALL atm_in%tracer(idx)%finalize()
    END DO
  END SUBROUTINE t_pi_atm_in_finalize


  SUBROUTINE t_init_state_const_finalize(const)
    CLASS(t_init_state_const), INTENT(INOUT) :: const

    const%linitialized = .FALSE.

    CALL const%vct%finalize()

    CALL DO_PTR_DEALLOCATE(const%z_mc_in)

    ! note: these pointers are not owned by this object
    NULLIFY(const%topography_c)
    NULLIFY(const%z_ifc)
    NULLIFY(const%z_mc)

  END SUBROUTINE t_init_state_const_finalize


  SUBROUTINE t_pi_sfc_in_finalize(sfc_in)
    CLASS(t_pi_sfc_in), INTENT(INOUT) :: sfc_in

    sfc_in%linitialized = .FALSE.
    CALL DO_DEALLOCATE(sfc_in%tsnow)
    CALL DO_DEALLOCATE(sfc_in%tskin)
    CALL DO_DEALLOCATE(sfc_in%sst)
    CALL DO_DEALLOCATE(sfc_in%snowalb)
    CALL DO_DEALLOCATE(sfc_in%snowweq)
    CALL DO_DEALLOCATE(sfc_in%snowdens)
    CALL DO_DEALLOCATE(sfc_in%skinres)
    CALL DO_DEALLOCATE(sfc_in%ls_mask)
    CALL DO_DEALLOCATE(sfc_in%seaice)
    CALL DO_DEALLOCATE(sfc_in%phi)
    CALL DO_DEALLOCATE(sfc_in%tsoil)
    CALL DO_DEALLOCATE(sfc_in%wsoil)
  END SUBROUTINE t_pi_sfc_in_finalize


  SUBROUTINE t_pi_atm_finalize(atm)
    CLASS(t_pi_atm), INTENT(INOUT) :: atm

    INTEGER :: idx

    atm%linitialized = .FALSE.
    CALL DO_DEALLOCATE(atm%vn)
    CALL DO_DEALLOCATE(atm%u)
    CALL DO_DEALLOCATE(atm%v)
    CALL DO_DEALLOCATE(atm%w)
    CALL DO_DEALLOCATE(atm%temp)
    CALL DO_DEALLOCATE(atm%theta_v)
    CALL DO_DEALLOCATE(atm%exner)
    CALL DO_DEALLOCATE(atm%rho)
    CALL DO_DEALLOCATE(atm%pres)
    CALL DO_DEALLOCATE(atm%qv)
    CALL DO_DEALLOCATE(atm%qc)
    CALL DO_DEALLOCATE(atm%qi)
    CALL DO_DEALLOCATE(atm%qr)
    CALL DO_DEALLOCATE(atm%qs)
    CALL DO_DEALLOCATE(atm%qg)
    CALL DO_DEALLOCATE(atm%qh)
    CALL DO_DEALLOCATE(atm%qnc)
    CALL DO_DEALLOCATE(atm%qni)
    CALL DO_DEALLOCATE(atm%qnr)
    CALL DO_DEALLOCATE(atm%qns)
    CALL DO_DEALLOCATE(atm%qng)
    CALL DO_DEALLOCATE(atm%qnh)
    CALL DO_DEALLOCATE(atm%tke)

    DO idx=1, ntracer
      CALL atm%tracer(idx)%finalize()
    END DO
  END SUBROUTINE t_pi_atm_finalize


  SUBROUTINE t_pi_sfc_finalize(sfc)
    CLASS(t_pi_sfc), INTENT(INOUT) :: sfc

    sfc%linitialized = .FALSE.
    CALL DO_DEALLOCATE(sfc%tsnow)
    CALL DO_DEALLOCATE(sfc%tskin)
    CALL DO_DEALLOCATE(sfc%sst)
    CALL DO_DEALLOCATE(sfc%snowalb)
    CALL DO_DEALLOCATE(sfc%snowweq)
    CALL DO_DEALLOCATE(sfc%snowdens)
    CALL DO_DEALLOCATE(sfc%skinres)
    CALL DO_DEALLOCATE(sfc%ls_mask)
    CALL DO_DEALLOCATE(sfc%seaice)
    CALL DO_DEALLOCATE(sfc%tsoil)
    CALL DO_DEALLOCATE(sfc%wsoil)
    CALL DO_DEALLOCATE(sfc%w_so)
  END SUBROUTINE t_pi_sfc_finalize


  SUBROUTINE t_sfc_inc_finalize(sfc_inc)
    CLASS(t_sfc_inc), INTENT(INOUT) :: sfc_inc

    sfc_inc%linitialized = .FALSE.
    CALL DO_DEALLOCATE(sfc_inc%w_so)
    CALL DO_DEALLOCATE(sfc_inc%h_snow)
    CALL DO_DEALLOCATE(sfc_inc%freshsnow)
    CALL DO_DEALLOCATE(sfc_inc%t_2m)
  END SUBROUTINE t_sfc_inc_finalize


  SUBROUTINE t_init_state_finalize(init_data)
    CLASS(t_init_state), INTENT(INOUT) :: init_data

    CALL init_data%atm_in%finalize()
    CALL init_data%sfc_in%finalize()
    CALL init_data%atm%finalize()
    CALL init_data%sfc%finalize()

    ! note: we do not call "initicon_data%const%finalize()", since
    ! this is a pointer to an external object which may be the target
    ! of another pointer.
    NULLIFY(init_data%const)
  END SUBROUTINE t_init_state_finalize


  SUBROUTINE t_initicon_state_finalize(init_data)
    CLASS(t_initicon_state), INTENT(INOUT) :: init_data

    ! call base class destructor
    CALL t_init_state_finalize(init_data)

    CALL init_data%atm_inc%finalize()
    CALL init_data%sfc_inc%finalize()
  END SUBROUTINE t_initicon_state_finalize

END MODULE mo_initicon_types
