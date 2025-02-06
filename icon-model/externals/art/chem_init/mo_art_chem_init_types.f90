!
! mo_chem_init_coord
! This module is based on mo_vertical_coord_table
!
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

MODULE mo_art_chem_init_types

  USE mo_kind,               ONLY: wp
  USE mo_art_impl_constants, ONLY: IART_VARNAMELEN

  IMPLICIT NONE

  PRIVATE

  ! Types
  PUBLIC :: t_art_chem_init_coord
  PUBLIC :: t_chem_init_state

  TYPE :: t_art_chem_init_coord
    INTEGER :: nlevm1         ! (number of levels)-1.
    INTEGER :: nplev          ! *number of pressure levels.
    INTEGER :: nplvp1         ! *nplev+1.*
    INTEGER :: nplvp2         ! *nplev+2.*
    INTEGER :: nplvpa         ! *nplvp1,* or 2 if *nplev=0.*
    INTEGER :: nlmsgl         ! *nlev* - (number of sigma levels).
    INTEGER :: nlmslp         ! *nlmsgl+1.*
    INTEGER :: nlmsla         ! *nlmslp,* or 2 if *nlmslp=1.*
  
    REAL(wp) ::apzero        ! *reference pressure for computation of the
    !                        !  hybrid vertical levels.
    REAL(wp) ::t0icao        ! *surface temperatur of reference atmosphere
    REAL(wp) ::tsticao       ! *stratospheric temperature of reference atmosphere
    REAL(wp) ::rdtstic       ! *rd*tsticao
    REAL(wp) ::rdlnp0i       ! *rd*ln(surface pressure) of reference atmosphere
    REAL(wp) ::alrrdic       ! *lapse-rate parameter of reference atmosphere
    REAL(wp) ::rdt0ral       ! *rd*t0icao/alphaic
    REAL(wp) ::ptricao       ! *tropopause pressure of reference atmosphere
    REAL(wp) ::rdlnpti       ! *rd*ln(ptricao)
    REAL(wp) ::gsticao       ! *constant used in geopotential calculation
  
    REAL(wp), ALLOCATABLE :: vct_chem_init(:) ! param. A and B of the vertical coordinate
  
    REAL(wp), ALLOCATABLE :: ralpha_chem_init(:) ! rd*alpha at pressure and sigma levels.
    REAL(wp), ALLOCATABLE :: rlnpr_chem_init(:)  ! rd*ln(p(k+.5)/p(k-.5))
    REAL(wp), ALLOCATABLE :: delpr_chem_init(:)  ! p(k+.5)-p(k-.5) of
                                       ! the refrence surface pressure.
    REAL(wp), ALLOCATABLE :: rdelpr_chem_init(:) ! *reciprocal of *delpr.*
  END TYPE t_art_chem_init_coord

  ! ------------------------------------------------------
  ! -     Structure for external data read in for chemical substances
  ! ------------------------------------------------------


  ! ----------------------------------
  ! --- t_chem_init_in   : fields needed for interpolation
  ! --- t_chem_init_chem : raw (spec) and interpolated chemical (spec_interp) species
  ! ----------------------------------

  TYPE :: t_chem_init_in

    LOGICAL :: linitialized  = .FALSE.
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: temp,q,temp_v,pres
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: z3d
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)     :: phi_sfc, ps
    INTEGER :: nlev_chem_init
    CHARACTER(LEN=IART_VARNAMELEN) :: model_name
  
    REAL(wp), ALLOCATABLE :: &
      &  pres_pl(:),         & !< pressure values at levels
      &  pres_pl3d(:,:,:),   & !< same as 3d field (for z_at_plevels)
      &  wfacpbl1(:,:),      & !< arrays for prepare_extrap
      &  wfacpbl2(:,:)         !< ...
    INTEGER, ALLOCATABLE  :: &
      &  kpbl1(:,:),  kpbl2(:,:) !< arrays for prepare_extrap
  END TYPE t_chem_init_in
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  TYPE :: t_chem_init_chem

    LOGICAL :: linitialized = .FALSE.
    INTEGER,POINTER :: n_spec_readin
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: spec
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: spec_interp

  END TYPE t_chem_init_chem
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  TYPE :: t_chem_init_state

    TYPE (t_art_chem_init_coord)  :: chem_init_coord
    TYPE (t_chem_init_in)         :: chem_init_in
    TYPE (t_chem_init_chem)       :: chem_init_chem

  END TYPE t_chem_init_state
  
  
END MODULE mo_art_chem_init_types


