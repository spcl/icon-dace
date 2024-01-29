!
! mo_art_modes
! This module provides types, parameters and initialization structures for modes
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

MODULE mo_art_modes
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_math_constants,                ONLY: pi
  USE mo_exception,                     ONLY: message, finish, message_text
  USE mo_fortran_tools,                 ONLY: init
  USE mo_impl_constants,                ONLY: SUCCESS
  USE mo_key_value_store,               ONLY: t_key_value_store
! ART
  USE mo_art_impl_constants,            ONLY: UNDEF_INT_ART, UNDEF_REAL_ART, &
                                          &   IART_VARNAMELEN
  USE mo_art_mode_fields,               ONLY: t_mode_fields
! ART aerosol dynamics processes routines
  USE mo_art_aerosol_utilities,         ONLY: art_modal_parameters, &
                                          &   art_modeshift
  USE mo_art_shift2mixed,               ONLY: art_shift2mixed
  USE mo_art_string_tools,              ONLY: key_value_storage_as_string
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_modes'

  TYPE t_optical_properties
    REAL(wp),POINTER  :: ext_coeff(:)       !< Mass specific extinction coefficient [m2 g-1], 
                                            !  dim:nbands
    REAL(wp),POINTER  :: ssa_coeff(:)       !< Single scattering albedo coefficient [-], dim:nbands
    REAL(wp),POINTER  :: asy_coeff(:)       !< Asymmetry coefficient [-], dim:nbands
    REAL(wp),POINTER  :: ext_param(:,:)     !< Mass specific extinction parameter for polynom [ ], 
                                            !  dim:nbands,npoly
    REAL(wp),POINTER  :: ssa_param(:,:)     !< Single scattering albedo parameter for polynom [ ], 
                                            !  dim:nbands,npoly
    REAL(wp),POINTER  :: asy_param(:,:)     !< Asymmetry coefficient parameter for polynom [ ], 
                                            !  dim:nbands,npoly
    REAL(wp),POINTER  :: ext_default(:,:)   !< Default value for extinction coefficient in case 
                                            !  diameter leaves boundaries [m2 g-1], dim:nbands,2 
                                            !  (dp > dp_max,dp < dp_min)
    REAL(wp),POINTER  :: ssa_default(:,:)   !< Default value for single scatter albedo in case 
                                            !  diameter leaves boundaries [-], dim:nbands,2
    REAL(wp),POINTER  :: asy_default(:,:)   !< Default value for asymmetry parameter in case 
                                            !  diameter leaves boundaries [-], dim:nbands,2
    REAL(wp),POINTER  :: dia_min_factor(:)  !< Lower border of parametrization validity as a 
                                            !  fraction of initial diameter [-]
    REAL(wp),POINTER  :: dia_max_factor(:)  !< Upper border of parametrization validity as a 
                                            !  fraction of initial diameter [-]
  END TYPE t_optical_properties

  TYPE t_mode_metadata
    REAL(wp)          :: diameter_ini_nmb   !< Initial Diameter number conc.
    REAL(wp)          :: diameter_ini_mass  !< Initial Diameter mass conc.
    REAL(wp)          :: sg_ini             !< Initial standard deviation
    REAL(wp)          :: init_nmb_conc      !< Initial number concentration of mode
    REAL(wp),POINTER  :: init_mass_conc(:)  !< Initial mass concentration of species (dim: ntr-1)
    REAL(wp)          :: exp_aero           !< Aerosol exponent
    REAL(wp)          :: mode_fac           !< from 3rd to 0th moment
    INTEGER           :: icondensation      !< Condensation scheme used for this mode
    INTEGER           :: inucleation        !< Nucleation scheme used for this mode
    INTEGER           :: itrnucl            !< Tracer index of nucleation tracer
    INTEGER           :: itrcond            !< Tracer index of condensation tracer
  END TYPE t_mode_metadata

  TYPE t_mode_connect
    LOGICAL              :: l_do_shift      !< Do this kind of shift for this mode?
    INTEGER              :: njsp            !< Number of species for mass mixing ratios
    INTEGER, ALLOCATABLE :: itr0(:)         !< Index pair number (dim1=2: shift from/to)
    INTEGER, ALLOCATABLE :: itr3(:,:)       !< Index pairs mass (dim1=2: shift from/to, dim2=njsp3)
    CHARACTER(LEN=IART_VARNAMELEN) :: shift2name !< Name of mode to shift to
    REAL(wp)             :: shift_diam      !< Optional: Diameter at which shift is performed
  END TYPE t_mode_connect
  
  TYPE coagPtr
    CLASS(t_mode_fields),POINTER :: p_coagulateWith   !< partner Mode to coagulate with
    CLASS(t_mode_fields),POINTER :: p_coagulateTo     !< resulting Mode by coagulating with 
                                                      !  other mode
    REAL(wp), POINTER           :: coagcoeff0(:,:,:)  !< coagulation coefficients 0th moment 
                                                      !  (dim1:dim2:dim3) [# m-3 s-1]
    REAL(wp), POINTER           :: coagcoeff3(:,:,:)  !< coagulation coefficients 3rd moment 
                                                      !  (dim1:dim2:dim3) [m3 m-3 s-1]
  END TYPE coagPtr
  
  TYPE t_coag_util
    INTEGER                     :: n_modes              !< How many coagulation-partner-modes exist
    TYPE(coagPtr), ALLOCATABLE  :: p_coag(:)            !< dimension: n_modes
  END TYPE t_coag_util

!-------------------------------------------------------------------------------
!-------       Hierarchy tree for t_mode_fields polymorphic objects      -------
!-------------------------------------------------------------------------------
!                                  t_mode_fields
!                t_fields_2mom <|                  |> t_fields_1mom
!                                   t_fields_pollen | t_fields_volc | t_fields_radio
!-------------------------------------------------------------------------------
  
  
  TYPE, extends(t_mode_fields) :: t_fields_1mom
    INTEGER               :: itr                     !< Index of species in tracer container
    REAL(wp),POINTER      :: flx_contra_vsed(:,:,:)  !< Flux due to sedimentation (mass or number)
    !
    CONTAINS
      PROCEDURE :: create   => create_t_fields_1mom
      PROCEDURE :: set_meta => set_meta_t_fields_1mom
      PROCEDURE :: print    => print_meta_t_fields_1mom
  END TYPE t_fields_1mom
  
  TYPE, extends(t_fields_1mom) :: t_fields_radio
    INTEGER    :: imis                               !< IMIS number of tracer
    REAL(wp)   :: halflife                           !< Half life of radioactive tracer
    REAL(wp)   :: vdep_const                         !< Constant deposition velocity
    REAL(wp)   :: fac_wetdep                         !< Wet deposition factor
    REAL(wp)   :: exp_wetdep                         !< Wet deposition exponent
    !
    CONTAINS
      PROCEDURE :: create   => create_t_fields_radio
      PROCEDURE :: set_meta => set_meta_t_fields_radio
      PROCEDURE :: print    => print_meta_t_fields_radio
  END TYPE t_fields_radio
  
  TYPE, extends(t_fields_1mom) :: t_fields_pollen
    REAL(wp)   :: diam                        !< Diameter of pollen
    REAL(wp)   :: rho                         !< Density of pollen
    INTEGER    :: doy_start_season            !< days since 1st December for start of pollen season
    INTEGER    :: doy_end_season              !< days since 1st December for end   of pollen season
    !
    CONTAINS
      PROCEDURE :: create   => create_t_fields_pollen
      PROCEDURE :: set_meta => set_meta_t_fields_pollen
      PROCEDURE :: print    => print_meta_t_fields_pollen
  END TYPE t_fields_pollen
  
  TYPE, extends(t_fields_1mom) :: t_fields_volc
    REAL(wp)   :: diam                               !< Diameter of volcanic ash
    REAL(wp)   :: rho                                !< Density of volcanic ash
    !
    CONTAINS
      PROCEDURE :: create   => create_t_fields_volc
      PROCEDURE :: set_meta => set_meta_t_fields_volc
      PROCEDURE :: print    => print_meta_t_fields_volc
  END TYPE t_fields_volc
  
  TYPE, extends(t_mode_fields) :: t_fields_2mom
    INTEGER,POINTER       :: itr3(:)                 !< Indices of mass species in tracer container
                                                     !  (dim: ntr-1)
    INTEGER               :: itr0                    !< Index of number species in tracer container
    LOGICAL,POINTER       :: linsol(:)               !< Mask for tracers which are insoluble 
                                                     !  (dim: ntr-1)
    TYPE(t_mode_metadata) :: info                    !< meta data for this entry
    TYPE(t_optical_properties) :: opt_props          !< Optical properties of mode
    REAL(wp),POINTER      :: third_moment(:,:,:)     !< 3rd moment
    REAL(wp),POINTER      :: mass(:,:,:)             !< total mass
    REAL(wp),POINTER      :: density(:,:,:)          !< Aerosol density
    REAL(wp),POINTER      :: diameter(:,:,:)         !< Diameter with respect to number conc.
    REAL(wp),POINTER      :: flx_contra_vsed0(:,:,:) !< flux due to sedimentation of zeroth moment
    REAL(wp),POINTER      :: flx_contra_vsed3(:,:,:) !< flux due to sedimentation of third moment
    REAL(wp),POINTER      :: dmdt_condso4(:,:,:)     !< Condensated SO4 mass
    REAL(wp),POINTER      :: knudsen_nr(:,:,:)       !< Knudsen Number
    REAL(wp),POINTER      :: rho(:)                  !< Density of the species (dim: ntr-1)
    INTEGER               :: ntr                     !< Number of species in mode 
                                                     !  (including number)
    TYPE(t_mode_connect)  :: shift2mixed             !< Connection of insoluble and mixed mode
    TYPE(t_mode_connect)  :: shift2larger            !< Connection of modes for shifting
    CLASS(t_mode_fields),POINTER :: p_shift2larger   !< Mode to shift to when threshold diam is 
                                                     !  exceeded
    CLASS(t_mode_fields),POINTER :: p_shift2mixed    !< Mode to shift to when threshold diam is 
                                                     !  exceeded
    INTEGER               :: do_coag                 !< 1 = inter+intra ; 2 = intra 
                                                     !  (maybe 0 = no coag)
    TYPE(t_coag_util)     :: coag_util               !< Type to handle coagulation
    LOGICAL               :: e2t_flag                !< Flag needed to realize emission2tracer
                                                     !<   mapping in initialization
    !
    CONTAINS
      PROCEDURE, PUBLIC  :: modal_param    => art_aerotbp_modpar
      PROCEDURE, PUBLIC  :: coagulation    => art_aerotbp_coag
      PROCEDURE, PUBLIC  :: mode_shift     => art_aerotbp_modeshift
      PROCEDURE, PUBLIC  :: create         => create_t_fields_2mom
      PROCEDURE, PUBLIC  :: set_meta       => set_meta_t_fields_2mom
      PROCEDURE, PUBLIC  :: print          => print_meta_t_fields_2mom
      PROCEDURE, PRIVATE :: update_number_1d
      PROCEDURE, PRIVATE :: update_number_2d
      PROCEDURE, PRIVATE :: update_number_2d_coag
      GENERIC, PUBLIC    :: update_number  => update_number_1d, update_number_2d, &
        &                                     update_number_2d_coag
      PROCEDURE, PRIVATE :: update_mass_1d
      PROCEDURE, PRIVATE :: update_mass_2d
      PROCEDURE, PRIVATE :: update_mass_2d_coag
      GENERIC, PUBLIC    :: update_mass    => update_mass_1d, update_mass_2d, update_mass_2d_coag
  END TYPE t_fields_2mom
  
  
  PUBLIC :: t_mode_fields
  PUBLIC :: t_mode_metadata
  PUBLIC :: t_fields_1mom
  PUBLIC :: t_fields_2mom
  PUBLIC :: t_fields_radio
  PUBLIC :: t_fields_pollen
  PUBLIC :: t_fields_volc
  
  CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_aerotbp_modpar(this_mode_fields,free_path, istart,iend,kstart,nlev,jb,tracer)
!<
! SUBROUTINE art_aerotbp_modpar
! Type bound procedure serving as interface to the routine
! calculating modal parameters art_modal_parameters
! Author: Daniel Rieger, KIT
! Initial Release: 2017-02-28
! Modifications:
! 2018-07-30: Lukas Muser, KIT
! - hand over initial values for nmb and mass concentration to modal_parameters 
!>
  CLASS(t_fields_2mom),INTENT(inout):: &
    &  this_mode_fields        !< Contains meta-information
  REAL(wp),INTENT(in)     :: &
    &  free_path(:,:)          !< Mean free path (Local)
  INTEGER, INTENT(in)     :: &
    &  istart, iend,         & !< Start and end indices of nproma loop
    &  kstart, nlev, jb        !< Number of vertical (full) levels and block index
  REAL(wp), INTENT(inout) :: &
    &  tracer(:,:,:)           !< Tracer field (jc,jk,ntracer)

  CALL art_modal_parameters(this_mode_fields%rho(:),                &
    &                       this_mode_fields%third_moment(:,:,jb),  &
    &                       this_mode_fields%mass(:,:,jb),          &
    &                       this_mode_fields%density(:,:,jb),       &
    &                       this_mode_fields%diameter(:,:,jb),      &
    &                       this_mode_fields%knudsen_nr(:,:,jb),    &
    &                       this_mode_fields%info%exp_aero,         &
    &                       this_mode_fields%info%diameter_ini_nmb, &
    &                       free_path(:,:),                         &
    &                       this_mode_fields%info%init_nmb_conc,    &
    &                       this_mode_fields%info%init_mass_conc(:),&
    &                       this_mode_fields%itr0,                  &
    &                       this_mode_fields%itr3(:),               &
    &                       this_mode_fields%ntr,                   &
    &                       istart,                                 &
    &                       iend,                                   &
    &                       kstart,                                 &
    &                       nlev,                                   &
    &                       tracer(:,:,:))

END SUBROUTINE art_aerotbp_modpar
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_aerotbp_coag(this_mode_fields,istart,iend,nlev,jb,tracer)
!<
! SUBROUTINE art_aerotbp_coag
! Type bound procedure serving as interface to the routine
! calculating coagulation art_coagulation
! Based on: -
! Author: Daniel Rieger, KIT
! Initial Release: 2017-02-28
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  CLASS(t_fields_2mom),INTENT(inout):: &
    &  this_mode_fields        !< Contains meta-information
  INTEGER, INTENT(in)     :: &
    &  istart, iend,         & !< Start and end indices of nproma loop
    &  nlev, jb                !< Number of vertical (full) levels and block index
  REAL(wp), INTENT(inout) :: &
    &  tracer(:,:,:)           !< Tracer field (jc,jk,ntracer)


  ! CALL art_coagulation()

END SUBROUTINE art_aerotbp_coag
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_aerotbp_modeshift(this_fields,istart,iend,kstart,nlev,jb,tracer)
!<
! SUBROUTINE art_aerotbp_modeshift
! Type bound procedure serving as interface to the routine
! calculating mode shifting art_modeshift
! Based on: -
! Author: Daniel Rieger, KIT
! Initial Release: 2017-02-28
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  CLASS(t_fields_2mom),INTENT(inout):: &
    &  this_fields             !< Contains meta-information
  CLASS(t_mode_fields), POINTER :: &
    &  that                    !< Contains meta-information mode to be shifted to
  INTEGER, INTENT(in)     :: &
    &  istart, iend,         & !< Start and end indices of nproma loop
    &  kstart, nlev, jb        !< Number of vertical (full) levels and block index
  REAL(wp), INTENT(inout) :: &
    &  tracer(:,:,:)           !< Tracer field (jc,jk,ntracer)
  ! Local variables
  INTEGER                 :: &
    &  jk, jc, jsp             !< Loop indices
  REAL(wp) :: &
    &  sol_mass(istart:iend,1:nlev) !< mass of soluble fraction

  IF (this_fields%shift2larger%l_do_shift) THEN
    that => this_fields%p_shift2larger
    SELECT TYPE(that_fields => that)
      CLASS IS(t_fields_2mom)
        CALL art_modeshift(this_fields%info%sg_ini,             &
          &                that_fields%info%sg_ini,             &
          &                this_fields%diameter(:,:,jb),        &
          &                that_fields%diameter(:,:,jb),        &
          &                this_fields%shift2larger%shift_diam, &
          &                this_fields%shift2larger%itr0,       &
          &                this_fields%shift2larger%itr3,       &
          &                this_fields%shift2larger%njsp,       &
          &                kstart,                              &
          &                nlev,                                &
          &                istart,                              &
          &                iend,                                &
          &                tracer(:,:,:))
      CLASS DEFAULT
        CALL finish(TRIM(routine)//':art_aerotbp_modeshift',          &
          &         'p_shift2larger is not of type t_fields_2mom for'//TRIM(this_fields%name)//'.')
    END SELECT
  ENDIF

  ! Convert insoluble to mixed mode
  IF (this_fields%shift2mixed%l_do_shift) THEN
    CALL init(sol_mass) 
    ! Determine soluble mass
    DO jk = kstart, nlev
      DO jc = istart, iend
        DO jsp = 1, this_fields%ntr - 1
          IF (.NOT. this_fields%linsol(jsp)) THEN
            sol_mass(jc,jk) = sol_mass(jc,jk) + tracer(jc,jk,this_fields%itr3(jsp))
          ENDIF
        ENDDO !jsp
      ENDDO !jc
    ENDDO !jk

    CALL art_shift2mixed(sol_mass(:,:),     &
      &                  this_fields%mass(:,:,jb),              &
      &                  tracer(:,:,:),                         &
      &                  istart,                                &
      &                  iend,                                  &
      &                  kstart,                                &
      &                  nlev,                                  &
      &                  this_fields%shift2mixed%njsp,          &
      &                  this_fields%shift2mixed%itr0,          &
      &                  this_fields%shift2mixed%itr3)
  ENDIF

END SUBROUTINE art_aerotbp_modeshift
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE create_t_fields_1mom(this_fields, modename, idims)
!<
! SUBROUTINE create_t_fields_1mom
! This subroutine performs a setup for the fields of type t_fields_1mom
! Part of Module: mo_art_modes
! Author: Daniel Rieger, KIT
! Initial Release: 2016-08-22
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_fields_1mom),INTENT(inout) :: &
    &  this_fields                        !< Container with fields
  CHARACTER(LEN=*),INTENT(in)        :: &
    &  modename                           !< Name of current mode
  INTEGER,INTENT(in)                 :: &
    &  idims(3)                           !< Dimensions to allocate fields
! Local variables
  INTEGER                            :: &
    &  istat                              !< Error stat

  istat = 0

  ! Calling create-SR from parent type
  CALL this_fields%t_mode_fields%create(modename,idims)
  this_fields%itr     = UNDEF_INT_ART

  ALLOCATE(this_fields%flx_contra_vsed(idims(1),(idims(2)+1),idims(3)),STAT=istat)
  ! idims(2)+1 as we need nlevp1 instead of nlev for sedimentative fluxes
  IF (istat /= 0) THEN
    CALL finish(TRIM(routine)//':create_t_fields_1mom',          &
      &         'Allocation failed for'//TRIM(this_fields%name))
  ENDIF

END SUBROUTINE create_t_fields_1mom
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE create_t_fields_radio(this_fields, modename, idims)
!<
! SUBROUTINE create_t_fields_radio
! This subroutine performs a setup for the fields of type t_fields_radio
! Part of Module: mo_art_modes
! Author: Daniel Rieger, KIT
! Initial Release: 2016-08-22
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_fields_radio),INTENT(inout) :: &
    &  this_fields                         !< Container with fields
  CHARACTER(LEN=*),INTENT(in)         :: &
    &  modename                            !< Name of current mode
  INTEGER,INTENT(in)                  :: &
    &  idims(3)                            !< Dimensions to allocate fields
! Local variables
  INTEGER                             :: &
    &  istat                               !< Error stat

  istat = 0

  ! Calling create-SR from parent type
  CALL this_fields%t_mode_fields%create(modename,idims)
  this_fields%itr     = UNDEF_INT_ART

  ALLOCATE(this_fields%flx_contra_vsed(idims(1),(idims(2)+1),idims(3)),STAT=istat)
  ! idims(2)+1 as we need nlevp1 instead of nlev for sedimentative fluxes
  IF (istat /= 0) THEN
    CALL finish(TRIM(routine)//':create_t_fields_radio',          &
      &         'Allocation failed for'//TRIM(this_fields%name))
  ENDIF

  this_fields%imis       = UNDEF_INT_ART
  this_fields%halflife   = UNDEF_REAL_ART
  this_fields%vdep_const = UNDEF_REAL_ART
  this_fields%fac_wetdep = UNDEF_REAL_ART
  this_fields%exp_wetdep = UNDEF_REAL_ART

END SUBROUTINE create_t_fields_radio
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE create_t_fields_pollen(this_fields, modename, idims)
!<
! SUBROUTINE create_t_fields_pollen
! This subroutine performs a setup for the fields of type t_fields_pollen
! Part of Module: mo_art_modes
! Author: Daniel Rieger, KIT
! Initial Release: 2016-08-22
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_fields_pollen),INTENT(inout) :: &
    &  this_fields                        !< Container with fields
  CHARACTER(LEN=*),INTENT(in)        :: &
    &  modename                           !< Name of current mode
  INTEGER,INTENT(in)                 :: &
    &  idims(3)                           !< Dimensions to allocate fields
! Local variables
  INTEGER                            :: &
    &  istat                              !< Error stat

  istat = 0

  ! Calling create-SR from parent type
  CALL this_fields%t_mode_fields%create(modename,idims)
  this_fields%itr     = UNDEF_INT_ART

  ALLOCATE(this_fields%flx_contra_vsed(idims(1),(idims(2)+1),idims(3)),STAT=istat)
  !$ACC ENTER DATA CREATE(this_fields%flx_contra_vsed)
  ! idims(2)+1 as we need nlevp1 instead of nlev for sedimentative fluxes
  IF (istat /= 0) THEN
    CALL finish(TRIM(routine)//':create_t_fields_pollen',          &
      &         'Allocation failed for'//TRIM(this_fields%name))
  ENDIF

  this_fields%diam = UNDEF_REAL_ART
  this_fields%rho  = UNDEF_REAL_ART

END SUBROUTINE create_t_fields_pollen
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE create_t_fields_volc(this_fields, modename, idims)
!<
! SUBROUTINE create_t_fields_volc
! This subroutine performs a setup for the fields of type t_fields_volc
! Part of Module: mo_art_modes
! Author: Daniel Rieger, KIT
! Initial Release: 2016-08-22
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_fields_volc),INTENT(inout) :: &
    &  this_fields                        !< Container with fields
  CHARACTER(LEN=*),INTENT(in)        :: &
    &  modename                           !< Name of current mode
  INTEGER,INTENT(in)                 :: &
    &  idims(3)                           !< Dimensions to allocate fields
! Local variables
  INTEGER                            :: &
    &  istat                              !< Error stat

  istat = 0

  ! Calling create-SR from parent type
  CALL this_fields%t_mode_fields%create(modename,idims)
  this_fields%itr     = UNDEF_INT_ART

  ALLOCATE(this_fields%flx_contra_vsed(idims(1),(idims(2)+1),idims(3)),STAT=istat)
  ! idims(2)+1 as we need nlevp1 instead of nlev for sedimentative fluxes
  IF (istat /= 0) THEN
    CALL finish(TRIM(routine)//':create_t_fields_volc',          &
      &         'Allocation failed for'//TRIM(this_fields%name))
  ENDIF

  this_fields%diam = UNDEF_REAL_ART
  this_fields%rho  = UNDEF_REAL_ART

END SUBROUTINE create_t_fields_volc
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE create_t_fields_2mom(this_fields, modename, idims)
!<
! SUBROUTINE create_t_mode_fields
! This subroutine performs a setup for the fields of type t_fields_2mom
! Part of Module: mo_art_modes
! Author: Daniel Rieger, KIT
! Initial Release: 2016-08-22
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_fields_2mom),INTENT(inout) :: &
    &  this_fields                        !< Container with fields
  CHARACTER(LEN=*),INTENT(in)        :: &
    &  modename                           !< Name of current mode
  INTEGER,INTENT(in)                 :: &
    &  idims(3)                           !< Dimensions to allocate fields
! Local variables
  INTEGER                            :: &
    &  istat, istat_all                   !< Error stat

  istat     = 0
  istat_all = 0

  ! Calling create-SR from parent type
  CALL this_fields%t_mode_fields%create(modename,idims)

  ! ALLOC
  ALLOCATE(this_fields%third_moment(idims(1),idims(2),idims(3)),    STAT=istat)
  istat_all = istat
  ALLOCATE(this_fields%mass(idims(1),idims(2),idims(3)),            STAT=istat)
  istat_all = istat_all + istat
  ALLOCATE(this_fields%density(idims(1),idims(2),idims(3)),         STAT=istat)
  istat_all = istat_all + istat
  ALLOCATE(this_fields%diameter(idims(1),idims(2),idims(3)),        STAT=istat)
  istat_all = istat_all + istat
  ALLOCATE(this_fields%flx_contra_vsed0(idims(1),(idims(2)+1),idims(3)),STAT=istat)
  istat_all = istat_all + istat
  ! idims(2)+1 as we need nlevp1 instead of nlev for sedimentative fluxes
  ALLOCATE(this_fields%flx_contra_vsed3(idims(1),(idims(2)+1),idims(3)),STAT=istat)
  istat_all = istat_all + istat
  ! idims(2)+1 as we need nlevp1 instead of nlev for sedimentative fluxes
  ALLOCATE(this_fields%dmdt_condso4(idims(1),idims(2),idims(3)),    STAT=istat)
  istat_all = istat_all + istat
  ALLOCATE(this_fields%knudsen_nr(idims(1),idims(2),idims(3)),      STAT=istat)
  istat_all = istat_all + istat
  IF (istat_all /= 0) THEN
    CALL finish(TRIM(routine)//':create_t_fields_2mom',          &
      &         'Allocation failed for'//TRIM(this_fields%name))
  ENDIF

  this_fields%info%diameter_ini_nmb  = UNDEF_REAL_ART
  this_fields%info%diameter_ini_mass = UNDEF_REAL_ART
  this_fields%info%sg_ini            = UNDEF_REAL_ART
  this_fields%info%init_nmb_conc     = UNDEF_REAL_ART
  this_fields%info%exp_aero          = UNDEF_REAL_ART
  this_fields%info%inucleation       = UNDEF_INT_ART
  this_fields%info%itrnucl           = UNDEF_INT_ART
  this_fields%info%icondensation     = UNDEF_INT_ART
  this_fields%info%itrcond           = UNDEF_INT_ART
  this_fields%ntr                    = 0
  this_fields%itr0                   = UNDEF_INT_ART

  this_fields%shift2larger%l_do_shift = .FALSE.
  this_fields%shift2larger%njsp       = 0
  this_fields%shift2larger%shift2name = ' '
  this_fields%shift2larger%shift_diam = UNDEF_REAL_ART

  this_fields%shift2mixed%l_do_shift  = .FALSE.
  this_fields%shift2mixed%njsp        = 0
  this_fields%shift2mixed%shift2name  = ' '
  this_fields%shift2mixed%shift_diam  = UNDEF_REAL_ART
  
  this_fields%p_shift2larger          => NULL()
  this_fields%p_shift2mixed           => NULL()

  this_fields%do_coag                 = 0
  this_fields%coag_util%n_modes       = 0
  this_fields%e2t_flag                = .FALSE.

END SUBROUTINE create_t_fields_2mom
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE set_meta_t_mode_fields(this_fields, meta_key_value_store)
!<
! SUBROUTINE set_meta_t_mode_fields
! This subroutine sets the metadata associated to type t_mode_fields
! Part of Module: mo_art_modes
! Author: Daniel Rieger, KIT
! Initial Release: 2016-08-22
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_mode_fields),INTENT(inout) :: &
    &  this_fields                        !< Container with fields
  TYPE(t_key_value_store),INTENT(in) :: &
    &  meta_key_value_store               !< Metadata container

  this_fields%linit = .TRUE.

END SUBROUTINE set_meta_t_mode_fields
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE set_meta_t_fields_1mom(this_fields, meta_key_value_store)
!<
! SUBROUTINE set_meta_t_fields_1mom
! This subroutine sets the metadata associated to type t_fields_1mom
! Part of Module: mo_art_modes
! Author: Daniel Rieger, KIT
! Initial Release: 2016-08-22
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_fields_1mom),INTENT(inout) :: &
    &  this_fields                        !< Container with fields
  TYPE(t_key_value_store),INTENT(in) :: &
    &  meta_key_value_store               !< Metadata container

  this_fields%linit = .TRUE.

END SUBROUTINE set_meta_t_fields_1mom
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE set_meta_t_fields_radio(this_fields, meta_key_value_store)
!<
! SUBROUTINE set_meta_t_fields_radio
! This subroutine sets the metadata associated to type t_fields_radio
! Part of Module: mo_art_modes
! Author: Daniel Rieger, KIT
! Initial Release: 2016-08-22
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_fields_radio),INTENT(inout) :: &
    &  this_fields                         !< Container with fields
  TYPE(t_key_value_store),INTENT(in)  :: &
    &  meta_key_value_store                !< Metadata container
  INTEGER                             :: &
    &  ierror                              !< Error check for key_value_store container

  this_fields%linit = .TRUE.
  ierror            = SUCCESS

  CALL meta_key_value_store%get('IMIS',this_fields%imis, ierror)
    IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':set_meta_t_fields_radio', &
      &                         'Metadata IMIS not available in key_value_store container.')
  CALL meta_key_value_store%get('halflife',this_fields%halflife, ierror)
    IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':set_meta_t_fields_radio', &
      &                         'Metadata halflife not available in key_value_store container.')
  CALL meta_key_value_store%get('vdep_const',this_fields%vdep_const, ierror)
    IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':set_meta_t_fields_radio', &
      &                         'Metadata vdep_const not available in key_value_store container.')
  CALL meta_key_value_store%get('fac_wetdep',this_fields%fac_wetdep, ierror)
    IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':set_meta_t_fields_radio', &
      &                         'Metadata fac_wetdep not available in key_value_store container.')
  CALL meta_key_value_store%get('exp_wetdep',this_fields%exp_wetdep, ierror)
    IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':set_meta_t_fields_radio', &
      &                         'Metadata exp_wetdep not available in key_value_store container.')

END SUBROUTINE set_meta_t_fields_radio
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE set_meta_t_fields_pollen(this_fields, meta_key_value_store)
!<
! SUBROUTINE set_meta_t_fields_pollen
! This subroutine sets the metadata associated to type t_fields_pollen
! Part of Module: mo_art_modes
! Author: Daniel Rieger, KIT
! Initial Release: 2016-08-22
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_fields_pollen),INTENT(inout) :: &
    &  this_fields                        !< Container with fields
  TYPE(t_key_value_store),INTENT(in)   :: &
    &  meta_key_value_store               !< Metadata container
  INTEGER                              :: &
    &  ierror                             !< Error check for key_value_store container

  this_fields%linit = .TRUE.
  ierror            = SUCCESS

  CALL meta_key_value_store%get('rho', this_fields%rho,  ierror)
    IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':set_meta_t_fields_pollen', &
      &                         'Metadata rho not available in key_value_store container.')
  CALL meta_key_value_store%get('diam',this_fields%diam, ierror)
    IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':set_meta_t_fields_pollen', &
      &                         'Metadata diam not available in key_value_store container.')
  CALL meta_key_value_store%get('doy_start_season',this_fields%doy_start_season, ierror)
    IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':set_meta_t_fields_pollen', &
      &                         'Metadata doy_start_season not available in storage container.')
  CALL meta_key_value_store%get('doy_end_season',this_fields%doy_end_season, ierror)
    IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':set_meta_t_fields_pollen', &
      &                         'Metadata doy_end_season not available in storage container.')

END SUBROUTINE set_meta_t_fields_pollen
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE set_meta_t_fields_volc(this_fields, meta_key_value_store)
!<
! SUBROUTINE set_meta_t_fields_volc
! This subroutine sets the metadata associated to type t_fields_volc
! Part of Module: mo_art_modes
! Author: Daniel Rieger, KIT
! Initial Release: 2016-08-22
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_fields_volc),INTENT(inout) :: &
    &  this_fields                        !< Container with fields
  TYPE(t_key_value_store),INTENT(in) :: &
    &  meta_key_value_store               !< Metadata container
  INTEGER                            :: &
    &  ierror                             !< Error check for key_value_store container

  this_fields%linit = .TRUE.
  ierror            = SUCCESS

  CALL meta_key_value_store%get('rho', this_fields%rho,  ierror)
    IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':set_meta_t_fields_volc', &
      &                                'Metadata rho not available in key_value_store container.')
  CALL meta_key_value_store%get('diam',this_fields%diam, ierror)
    IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':set_meta_t_fields_volc', &
      &                                'Metadata diam not available in key_value_store container.')

END SUBROUTINE set_meta_t_fields_volc
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE set_meta_t_fields_2mom(this_fields, meta_key_value_store)
!<
! SUBROUTINE set_meta_t_fields_2mom
! This subroutine sets the metadata associated to type t_fields_2mom
! Part of Module: mo_art_modes
! Author: Daniel Rieger, KIT
! Initial Release: 2016-08-22
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_fields_2mom),INTENT(inout) :: &
    &  this_fields                        !< Container with fields
  TYPE(t_key_value_store),INTENT(in) :: &
    &  meta_key_value_store               !< Metadata container
  INTEGER                            :: &
    &  ierror                             !< Error check for key_value_store container
  CHARACTER(LEN=:),ALLOCATABLE       :: &
    &  c_tmp

  this_fields%linit = .TRUE.
  ierror            = SUCCESS

! Geometric quantities
  CALL meta_key_value_store%get('d_gn',    this_fields%info%diameter_ini_nmb,  ierror)
  IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':set_meta_t_fields_2mom', &
    &                         'Metadata d_gn not available in key_value_store container.')
  CALL meta_key_value_store%get('sigma_g', this_fields%info%sg_ini,            ierror)
  IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':set_meta_t_fields_2mom', &
    &                         'Metadata sigma_g not available in key_value_store container.')
  CALL meta_key_value_store%get('d_gm',    this_fields%info%diameter_ini_mass, ierror)
  IF (ierror /= SUCCESS) THEN ! calculate it from d_gn and sigma_g
    this_fields%info%diameter_ini_mass = EXP( LOG(this_fields%info%diameter_ini_nmb)   &
      &                                +      (3._wp * LOG(this_fields%info%sg_ini)    &
      &                                *       LOG(this_fields%info%sg_ini) ) )
  ENDIF
  CALL meta_key_value_store%get('icoag',   this_fields%do_coag,                ierror)
  IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':set_meta_t_fields_2mom', &
    &                         'Metadata icoag not available in key_value_store container.')

  this_fields%info%exp_aero = EXP( 0.125_wp * (LOG(this_fields%info%sg_ini)**2))

! Mode shifting
  ! shift2larger
  CALL key_value_storage_as_string(meta_key_value_store,'shift2larger',c_tmp,ierror)
  IF (ierror /= SUCCESS) THEN ! no shifting
    this_fields%shift2larger%l_do_shift = .false.
  ELSE
    WRITE(this_fields%shift2larger%shift2name,'(A)') c_tmp   
    this_fields%shift2larger%l_do_shift = .true.
  ENDIF
  IF(TRIM(this_fields%shift2larger%shift2name) == TRIM(this_fields%name)) THEN
    CALL finish(TRIM(routine)//':set_meta_t_fields_2mom', &
      &         'Shifting mode '//TRIM(this_fields%name)//' to itself not allowed.')
  ENDIF
  IF (this_fields%shift2larger%l_do_shift) THEN
    CALL meta_key_value_store%get('shift_diam',this_fields%shift2larger%shift_diam, ierror)
    IF (ierror /= SUCCESS) THEN
      CALL finish(TRIM(routine)//':set_meta_t_fields_2mom', &
        &         'No diameter for shifting mode '//TRIM(this_fields%name)//' specified.')
    ENDIF
  ENDIF
  ! shift2mixed
  CALL key_value_storage_as_string(meta_key_value_store,'shift2mixed',c_tmp,ierror)
  IF (ierror /= SUCCESS) THEN ! no shifting
    this_fields%shift2mixed%l_do_shift = .false.
  ELSE
    WRITE(this_fields%shift2mixed%shift2name,'(A)') c_tmp   
    this_fields%shift2mixed%l_do_shift = .true.
  ENDIF
  IF(TRIM(this_fields%shift2mixed%shift2name) == TRIM(this_fields%name)) THEN
    CALL finish(TRIM(routine)//':set_meta_t_fields_2mom', &
      &         'Shifting mode '//TRIM(this_fields%name)//' to itself not allowed.')
  ENDIF

! Process control
  CALL meta_key_value_store%get('condensation',this_fields%info%icondensation,ierror)
  IF (ierror /= SUCCESS) this_fields%info%icondensation = 0

END SUBROUTINE set_meta_t_fields_2mom
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE print_meta_t_fields_1mom(this_fields)
!<
! SUBROUTINE print_meta_t_fields_1mom
! This subroutine prints the metadata contained in t_fields_1mom
! Part of Module: mo_art_modes
! Author: Daniel Rieger, KIT
! Initial Release: 2016-09-15
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_fields_1mom),INTENT(in)  :: &
    &  this_fields                      !< Fields of current mode

  IF(this_fields%linit) THEN
    WRITE (message_text,*) '==========ART: PRINTOUT MODEL METADATA FOR MODE '// &
      &                    TRIM(this_fields%name)//'=========='
    CALL message ('', message_text)
    WRITE (message_text,*) 'NAME            DATA'
    CALL message ('', message_text)
    WRITE (message_text,'(A16,I3)')    'Tracer index:   ',this_fields%itr
    CALL message ('', message_text)
  ELSE
    CALL finish(TRIM(routine)//':print_meta_t_fields_1mom',          &
      &         'Print not possible: Mode '//TRIM(this_fields%name)//' was not initialized.')
  ENDIF

END SUBROUTINE print_meta_t_fields_1mom
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE print_meta_t_fields_radio(this_fields)
!<
! SUBROUTINE print_meta_t_fields_radio
! This subroutine prints the metadata contained in t_fields_radio
! Part of Module: mo_art_modes
! Author: Daniel Rieger, KIT
! Initial Release: 2016-09-15
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_fields_radio),INTENT(in) :: &
    &  this_fields                      !< Fields of current mode

  IF(this_fields%linit) THEN
    WRITE (message_text,*) '==========ART: PRINTOUT MODEL METADATA FOR MODE '// &
      &                    TRIM(this_fields%name)//'=========='
    CALL message ('', message_text)
    WRITE (message_text,*) 'NAME            DATA'
    CALL message ('', message_text)
    WRITE (message_text,'(A15,I4)')    'IMIS:           ',this_fields%imis
    CALL message ('', message_text)
    WRITE (message_text,'(A15,E13.6)') 'Halflife:       ',this_fields%halflife
    CALL message ('', message_text)
    WRITE (message_text,'(A15,E13.6)') 'Dep. velocity:  ',this_fields%vdep_const
    CALL message ('', message_text)
    WRITE (message_text,'(A15,E13.6)') 'Fac. wet dep.:  ',this_fields%fac_wetdep
    CALL message ('', message_text)
    WRITE (message_text,'(A15,E13.6)') 'Exp. wet dep.:  ',this_fields%exp_wetdep
    CALL message ('', message_text)
    WRITE (message_text,'(A16,I3)')    'Tracer index:   ',this_fields%itr
    CALL message ('', message_text)
  ELSE
    CALL finish(TRIM(routine)//':print_meta_t_fields_radio',          &
      &         'Print not possible: Mode '//TRIM(this_fields%name)//' was not initialized.')
  ENDIF

END SUBROUTINE print_meta_t_fields_radio
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE print_meta_t_fields_pollen(this_fields)
!<
! SUBROUTINE print_meta_t_fields_pollen
! This subroutine prints the metadata contained in t_fields_pollen
! Part of Module: mo_art_modes
! Author: Daniel Rieger, KIT
! Initial Release: 2016-09-15
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_fields_pollen),INTENT(in)  :: &
    &  this_fields                      !< Fields of current mode

  IF(this_fields%linit) THEN
    WRITE (message_text,*) '==========ART: PRINTOUT MODEL METADATA FOR MODE '//&
      &                    TRIM(this_fields%name)//'=========='
    CALL message ('', message_text)
    WRITE (message_text,*) 'NAME                 DATA'
    CALL message ('', message_text)
    WRITE (message_text,'(A21,E13.6)') 'Diameter:            ',this_fields%diam
    CALL message ('', message_text)
    WRITE (message_text,'(A21,E13.6)') 'Density:             ',this_fields%rho
    CALL message ('', message_text)
    WRITE (message_text,'(A21,I3)')    'DoY Start of season: ',this_fields%doy_start_season
    CALL message ('', message_text)
    WRITE (message_text,'(A21,I3)')    'DoY End   of season: ',this_fields%doy_end_season
    CALL message ('', message_text)
    WRITE (message_text,'(A21,I3)')    'Tracer index:        ',this_fields%itr
    CALL message ('', message_text)
  ELSE
    CALL finish(TRIM(routine)//':print_meta_t_fields_pollen',          &
      &         'Print not possible: Mode '//TRIM(this_fields%name)//' was not initialized.')
  ENDIF

END SUBROUTINE print_meta_t_fields_pollen
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE print_meta_t_fields_volc(this_fields)
!<
! SUBROUTINE print_meta_t_fields_volc
! This subroutine prints the metadata contained in t_fields_volc
! Part of Module: mo_art_modes
! Author: Daniel Rieger, KIT
! Initial Release: 2016-09-15
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_fields_volc),INTENT(in)  :: &
    &  this_fields                      !< Fields of current mode

  IF(this_fields%linit) THEN
    WRITE (message_text,*) '==========ART: PRINTOUT MODEL METADATA FOR MODE '// &
      &                    TRIM(this_fields%name)//'=========='
    CALL message ('', message_text)
    WRITE (message_text,*) 'NAME            DATA'
    CALL message ('', message_text)
    WRITE (message_text,'(A15,E13.6)') 'Diameter:       ',this_fields%diam
    CALL message ('', message_text)
    WRITE (message_text,'(A15,E13.6)') 'Density:        ',this_fields%rho
    CALL message ('', message_text)
    WRITE (message_text,'(A16,I3)')    'Tracer index:   ',this_fields%itr
    CALL message ('', message_text)
  ELSE
    CALL finish(TRIM(routine)//':print_meta_t_fields_volc',          &
      &         'Print not possible: Mode '//TRIM(this_fields%name)//' was not initialized.')
  ENDIF

END SUBROUTINE print_meta_t_fields_volc
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE print_meta_t_fields_2mom(this_fields)
!<
! SUBROUTINE print_meta_t_fields_2mom
! This subroutine prints the metadata contained in t_fields_2mom
! Part of Module: mo_art_modes
! Author: Daniel Rieger, KIT
! Initial Release: 2016-09-15
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_fields_2mom),INTENT(in)  :: &
    &  this_fields                      !< Fields of current mode
!Local variables
  INTEGER                          :: &
    &  jsp
  CHARACTER(LEN=80)                :: &
    &  jspstring

  jspstring = '                                                                                '

  IF(this_fields%linit) THEN
    WRITE (message_text,*) '==========ART: PRINTOUT MODEL METADATA FOR MODE '// &
      &                    TRIM(this_fields%name)//'=========='
    CALL message ('', message_text)
    WRITE (message_text,*) 'NAME            DATA'
    CALL message ('', message_text)
    WRITE (message_text,'(A16,E13.6)') 'DG0:            ',this_fields%info%diameter_ini_nmb
    CALL message ('', message_text)
    WRITE (message_text,'(A16,E13.6)') 'DG3:            ',this_fields%info%diameter_ini_mass
    CALL message ('', message_text)
    WRITE (message_text,'(A16,E13.6)') 'SIGMAG:         ',this_fields%info%sg_ini
    CALL message ('', message_text)
    WRITE (message_text,'(A16,E13.6)') 'Aerosol Exp.:   ',this_fields%info%exp_aero
    CALL message ('', message_text)
    WRITE (message_text,'(A16,I3)')    'Nmb. of Species:',this_fields%ntr
    CALL message ('', message_text)

    IF (this_fields%info%inucleation > 0) THEN
      WRITE (message_text,'(A16,I3)')    'Nucleation:     ',this_fields%info%inucleation
      CALL message ('', message_text)
      WRITE (message_text,'(A16,I3)')    'Nucl. Species:  ',this_fields%info%itrnucl
      CALL message ('', message_text)
    ELSE
      WRITE (message_text,'(A16)')       'No Nucleation.  '
      CALL message ('', message_text)
    ENDIF

    IF (this_fields%info%icondensation > 0) THEN
      WRITE (message_text,'(A16,I3)')    'Condensation:   ',this_fields%info%icondensation
      CALL message ('', message_text)
      WRITE (message_text,'(A16,I3)')    'Cond. Species:  ',this_fields%info%itrcond
      CALL message ('', message_text)
    ELSE
      WRITE (message_text,'(A16)')       'No Condensation.'
      CALL message ('', message_text)
    ENDIF

    DO jsp = 1,this_fields%ntr-1
      WRITE (message_text,'(A16,E13.6)') 'Rho:            ',this_fields%rho(jsp)
      CALL message ('', message_text)
      WRITE(jspstring((4*jsp-3):(4*jsp)),'(I3,A1)') this_fields%itr3(jsp),','
    ENDDO
    
    WRITE(jspstring((4*this_fields%ntr-3):(4*this_fields%ntr)),'(I3,A1)') this_fields%itr0,'.'
    WRITE (message_text,*) 'Species contained: '//jspstring
    CALL message ('', message_text)

    IF (this_fields%shift2larger%l_do_shift) THEN
      WRITE (message_text,'(A16,A)') 'shift2larger:   ',TRIM(this_fields%shift2larger%shift2name)
      CALL message ('', message_text)
      WRITE (message_text,'(A16,E13.6)') 'Shift Diam.:    ',this_fields%shift2larger%shift_diam
      CALL message ('', message_text)
    ENDIF

    IF (this_fields%shift2mixed%l_do_shift) THEN
      WRITE (message_text,'(A16,A)') 'shift2mixed:    ',TRIM(this_fields%shift2mixed%shift2name)
      CALL message ('', message_text)
    ENDIF

  ELSE
    CALL finish(TRIM(routine)//':print_meta_t_fields_2mom',          &
      &         'Print not possible: Mode '//TRIM(this_fields%name)//' was not initialized.')
  ENDIF

END SUBROUTINE print_meta_t_fields_2mom
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE update_number_1d(this_fields, tracer, update_rate, dtime, istart, iend, opt_rho)
!<
! SUBROUTINE update_number_1d
! This subroutine updates the specific number
! of tracers associated to 2mom mode for 2d update rates (i.e. horizontal)
! Part of Module: mo_art_modes
! Author: Daniel Rieger, KIT
! Initial Release: 2017-04-20
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_fields_2mom),INTENT(in) :: &
    &  this_fields                     !< Fields of current mode
  REAL(wp),INTENT(inout)          :: &
    &  tracer(:)                     !< Field to be updated (kg-1)
  REAL(wp),INTENT(in)             :: &
    &  update_rate(:),               & !< Update rate (kg-1 s-1) or (m-3 s-1)
    &  dtime                           !< Integration time step (s)
  INTEGER,INTENT(in)              :: &
    &  istart,iend                     !< Loop indizes (start and end)
  REAL(wp),OPTIONAL               :: &
    &  opt_rho(:)                      !< Air density, needed if update_rate is in (m-3 s-1)
! Local variables
  INTEGER                         :: &
    &  jc                              !< Loop index

  IF(PRESENT(opt_rho)) THEN
!NEC$ ivdep
    DO jc = istart, iend
      tracer(jc) = tracer(jc) + update_rate(jc) * dtime * (1._wp/opt_rho(jc))
    ENDDO
  ELSE
!NEC$ ivdep
    DO jc = istart, iend
      tracer(jc) = tracer(jc) + update_rate(jc) * dtime
    ENDDO
  ENDIF

END SUBROUTINE update_number_1d
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE update_number_2d(this_fields, tracer, update_rate, dtime,           &
  &                         istart, iend, kstart, nlev, rho)
!<
! SUBROUTINE update_number_2d
! This subroutine updates the specific number
! of tracers associated to 2mom mode for 2d update rates (i.e. horizontal and vertical)
! Part of Module: mo_art_modes
! Author: Daniel Rieger, KIT
! Initial Release: 2017-04-20
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_fields_2mom),INTENT(in) :: &
    &  this_fields                     !< Fields of current mode
  REAL(wp),INTENT(inout)          :: &
    &  tracer(:,:)                   !< Field to be updated (kg-1)
  REAL(wp),INTENT(in)             :: &
    &  update_rate(:,:),             & !< Update rate (# m-3 s-1)
    &  dtime,                        & !< Integration time step (s)
    &  rho(:,:)                        !< Air density
  INTEGER,INTENT(in)              :: &
    &  istart, iend, kstart, nlev      !< Loop indizes (start and end), number of vertical levels
! Local variables
  INTEGER                         :: &
    &  jc, jk                          !< Loop index

  DO jk = kstart, nlev
!NEC$ ivdep
    DO jc = istart, iend
      tracer(jc,jk) = tracer(jc,jk) + update_rate(jc,jk) * dtime * (1._wp/rho(jc,jk))
    ENDDO
  ENDDO

END SUBROUTINE update_number_2d
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE update_number_2d_coag(this_fields, tracer, jb, dtime,               &
  &                              istart, iend, kstart, nlev, rho)
!<
! SUBROUTINE update_number_2d_coag
! This subroutine updates the specific number
! of tracers associated to 2mom mode for 2d update rates (i.e. horizontal and vertical)
! for the coagulation case - update rates are stored in type itself
! 
! Method: (rethink this - other variants are possible)
!     - remove the rate from this field
!     - add half the rate to target field
!     (other half will be considered when handling partner field)
!     - exception: self-coagulation (directly consider full rate - there is no second handling)
! Part of Module: mo_art_modes
! Author: Daniel Rieger, KIT
! Initial Release: 2017-04-20
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_fields_2mom),INTENT(in) :: &
    &  this_fields                     !< Fields of current mode
  REAL(wp),INTENT(inout)          :: &
    &  tracer(:,:,:,:)                   !< Field to be updated (kg-1)
  INTEGER, INTENT(in)             :: &
    &  jb                              !< Current block index
  REAL(wp),INTENT(in)             :: &
    &  dtime                           !< Integration time step (s)
  INTEGER,INTENT(in)              :: &
    &  istart, iend, kstart, nlev      !< Loop indizes (start and end), number of vertical levels
  REAL(wp),INTENT(in)             :: &
    &  rho(:,:)                        !< Air density, needed if update_rate is in (m-3 s-1)
! Local variables
  INTEGER                         :: &
    &  jn                              !< Loop index
    
  DO jn = 1, this_fields%coag_util%n_modes
    IF (TRIM(ADJUSTL(this_fields%name)) ==                                                        &
      & TRIM(ADJUSTL(this_fields%coag_util%p_coag(jn)%p_coagulateWith%name))) THEN
      CALL this_fields%update_number(tracer(:,:,jb,this_fields%itr0),                                                      &
        &                            this_fields%coag_util%p_coag(jn)%coagcoeff0(:,:,jb)*(-1.0_wp),  &
        &                            dtime,istart,iend,kstart,nlev,rho(:,:))
    ELSE
      CALL this_fields%update_number(tracer(:,:,jb,this_fields%itr0),                                                      &
        &                            this_fields%coag_util%p_coag(jn)%coagcoeff0(:,:,jb)*(-1.0_wp),  &
        &                            dtime,istart,iend,kstart,nlev,rho(:,:))
      SELECT TYPE(coagTo => this_fields%coag_util%p_coag(jn)%p_coagulateTo)
        CLASS IS(t_fields_2mom)
          CALL coagTo%update_number(tracer(:,:,jb,coagTo%itr0),                                                       &
            &                       this_fields%coag_util%p_coag(jn)%coagcoeff0(:,:,jb)*0.5_wp,      &
            &                       dtime,istart,iend,kstart,nlev,rho(:,:))
        CLASS DEFAULT
          CALL finish(TRIM(routine)//':update_number_2d_coag',                                    &
            &         'not t_fields_2mom, though it should be')
      END SELECT
    END IF
  
  END DO
  
END SUBROUTINE update_number_2d_coag
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE update_mass_1d(this_fields, jb, tracer, update_rate, dtime,          &
  &                         istart, iend, jk, opt_rho)
!<
! SUBROUTINE update_mass_1d
! This subroutine updates the specific mass
! of tracers associated to 2mom mode for 1d update rates (i.e. horizontal )
! Part of Module: mo_art_modes
! Author: Sven Werchner, KIT
! Initial Release: 2020-11-09
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_fields_2mom),INTENT(in) :: &
    &  this_fields                     !< Fields of current mode
  REAL(wp),INTENT(inout)          :: &
    &  tracer(:,:)                   !< Field to be updated (kg-1)
  REAL(wp),INTENT(in)             :: &
    &  update_rate(:),             & !< Update rate (kg-1 s-1) or (m-3 s-1)
    &  dtime                           !< Integration time step (s)
  INTEGER,INTENT(in)              :: &
    &  istart, iend,                 &!< Loop indizes (start and end)
    &  jb, jk                        !< Block ID and vertical level
  REAL(wp),OPTIONAL               :: &
    &  opt_rho(:)                    !< Air density, needed if update_rate is in (m-3 s-1)
! Local variables
  INTEGER                         :: &
    &  jc, jsp                     !< Loop index

  IF(PRESENT(opt_rho)) THEN
    DO jsp = 1, this_fields%ntr-1
!NEC$ ivdep
      DO jc = istart, iend
        tracer(jc,this_fields%itr3(jsp)) = tracer(jc,this_fields%itr3(jsp))                  &
          &                                 * (1._wp + update_rate(jc) * dtime               &
          &                                 / this_fields%third_moment(jc,jk,jb) / opt_rho(jc))
      ENDDO
    ENDDO
  ELSE
    DO jsp = 1, this_fields%ntr-1
!NEC$ ivdep
      DO jc = istart, iend
        tracer(jc,this_fields%itr3(jsp)) = tracer(jc,this_fields%itr3(jsp))      &
          &                                 * (1._wp + update_rate(jc) * dtime      &
          &                                 / this_fields%third_moment(jc,jk,jb)) 
      ENDDO
    ENDDO
  ENDIF

END SUBROUTINE update_mass_1d
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE update_mass_2d(this_fields, jb, tracer, update_rate, dtime,          &
  &                         istart, iend, kstart, nlev, rho)
!<
! SUBROUTINE update_mass_2d
! This subroutine updates the specific mass
! of tracers associated to 2mom mode for 2d update rates (i.e. horizontal and vertical)
! Part of Module: mo_art_modes
! Author: Sven Werchner, KIT
! Initial Release: 2020-11-09
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_fields_2mom),INTENT(in) :: &
    &  this_fields                     !< Fields of current mode
  REAL(wp),INTENT(inout)          :: &
    &  tracer(:,:,:,:)                 !< Field to be updated (kg-1)
  REAL(wp),INTENT(in)             :: &
    &  update_rate(:,:),             & !< Update rate (m3 m-3 s-1)
    &  dtime,                        & !< Integration time step (s)
    &  rho(:,:)                        !< Air density
  INTEGER,INTENT(in)              :: &
    &  istart, iend, kstart, nlev,   & !< Loop indizes (start and end), number of vertical levels
    &  jb                              !< block ID
! Local variables
  INTEGER                         :: &
    &  jc, jk, jsp                     !< Loop index

  DO jsp = 1, this_fields%ntr-1
    DO jk = kstart, nlev
!NEC$ ivdep
      DO jc = istart, iend
        tracer(jc,jk,jb,this_fields%itr3(jsp)) = tracer(jc,jk,jb,this_fields%itr3(jsp))               &
          &                                 * (1._wp + update_rate(jc,jk) * dtime               &
          &                                 / this_fields%third_moment(jc,jk,jb) /rho(jc,jk))
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE update_mass_2d
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE update_mass_2d_coag(this_fields, tracer, jb, dtime, istart, iend, kstart, nlev, rho)
!<
! SUBROUTINE update_number_2d_coag
! This subroutine updates the specific mass
! of tracers associated to 2mom mode for 2d update rates (i.e. horizontal and vertical)
! for the coagulation case - update rates are stored in type itself
! 
! Method: (rethink this - other variants are possible)
!     - remove corresponding rate from this field
!     - add corresponding rate to target field
!     (partner field rates are considered when partner calls)
! << SYMMETRICAL COAG-MATRIX >>
!
! Part of Module: mo_art_modes
! Author: Daniel Rieger, KIT
! Initial Release: 2017-04-20
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_fields_2mom),INTENT(in) :: &
    &  this_fields                     !< Fields of current mode
  REAL(wp),INTENT(inout)          :: &
    &  tracer(:,:,:,:)                   !< Field to be updated (kg-1)
  INTEGER, INTENT(in)             :: &
    &  jb                              !< Current block index
  REAL(wp),INTENT(in)             :: &
    &  dtime                           !< Integration time step (s)
  INTEGER,INTENT(in)              :: &
    &  istart, iend, kstart, nlev      !< Loop indizes (start and end), number of vertical levels
  REAL(wp), INTENT(in)               :: &
    &  rho(:,:)                    !< Air density, needed if update_rate is in (m-3 s-1)
! Local variables
  INTEGER                         :: &
    &  jk,jc,jn,ijsp                   !< Loop indices
  REAL(wp)                                     :: &
    &  rate            (istart:iend,kstart:nlev), &
    &  new_total_this  (istart:iend,kstart:nlev), & !< new total_mass of this_fields
    &  new_total_target(istart:iend,kstart:nlev)    !< new total_mass of target_fields
  !CLASS(t_mode_fields),POINTER     :: &
  !  &  target_fields

  DO jn = 1,this_fields%coag_util%n_modes
    IF (TRIM(ADJUSTL(this_fields%name)) ==                                                    &
      & TRIM(ADJUSTL(this_fields%coag_util%p_coag(jn)%p_coagulateTo%name))) THEN
      CYCLE
    END IF
    
    DO jk = kstart, nlev
!NEC$ ivdep
      DO jc = istart, iend
        IF (this_fields%coag_util%p_coag(jn)%coagcoeff3(jc,jk,jb) == 0.0_wp) THEN
          CYCLE
        END IF
        rate(jc,jk) = this_fields%coag_util%p_coag(jn)%coagcoeff3(jc,jk,jb) * dtime               &
          &         * this_fields%density(jc,jk,jb) / rho(jc,jk)
        new_total_this(jc,jk) = 0.0_wp
      ENDDO !jc
    ENDDO !jk

    ! this_fields -> loss
    DO ijsp = 1, this_fields%ntr-1
      DO jk = kstart, nlev
!NEC$ ivdep
        DO jc = istart, iend
          tracer(jc,jk,jb,this_fields%itr3(ijsp)) = tracer(jc,jk,jb,this_fields%itr3(ijsp))             &
            &                                  * ( 1.0_wp - rate(jc,jk) / this_fields%mass(jc,jk,jb))
          new_total_this(jc,jk) = new_total_this(jc,jk) + tracer(jc,jk,jb,this_fields%itr3(ijsp))
        ENDDO !jc
      ENDDO !jk
    ENDDO !ijsp

    DO jk = kstart, nlev
!NEC$ ivdep
      DO jc = istart, iend
        this_fields%mass(jc,jk,jb) = new_total_this(jc,jk)
        
        ! target_fields -> gain
        new_total_target(jc,jk) = 0.0_wp
      ENDDO !jc
    ENDDO !jk

    SELECT TYPE(target_fields => this_fields%coag_util%p_coag(jn)%p_coagulateTo)
      CLASS IS(t_fields_2mom)
        DO ijsp = 1, target_fields%ntr-1
          DO jk = kstart, nlev
!NEC$ ivdep
            DO jc = istart, iend
              tracer(jc,jk,jb,target_fields%itr3(ijsp)) = tracer(jc,jk,jb,target_fields%itr3(ijsp))     &
                &                                    * ( 1._wp + rate(jc,jk)                          &
                &                                    / target_fields%mass(jc,jk,jb))
              new_total_target(jc,jk) = new_total_target(jc,jk)                                   &
                &                     + tracer(jc,jk,jb,target_fields%itr3(ijsp))
            ENDDO !jc
          ENDDO !jk
        ENDDO !ijsp

        DO jk = kstart, nlev
!NEC$ ivdep
          DO jc = istart, iend
            target_fields%mass(jc,jk,jb) = new_total_target(jc,jk)
          ENDDO !jc
        ENDDO !jk

      CLASS DEFAULT
        CALL finish(TRIM(routine)//':update_mass_2d_coag',          &
          &         'not t_fields_2mom, though it should be')
      
    END SELECT
        
  ENDDO !jn
 
END SUBROUTINE update_mass_2d_coag
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_modes
