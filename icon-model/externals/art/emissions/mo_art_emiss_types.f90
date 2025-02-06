!
! mo_art_emiss_types
! This module provides data storage structures for emission datasets
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

MODULE mo_art_emiss_types
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_math_constants,                ONLY: pi
  USE mo_exception,                     ONLY: finish
  USE mo_impl_constants,                ONLY: SUCCESS
  USE mo_key_value_store,               ONLY: t_key_value_store
  USE mo_util_string,                   ONLY: split_string

  USE mtime,                            ONLY: deallocateDatetime
! ART
  USE mo_art_impl_constants,            ONLY: UNDEF_INT_ART, IART_VARNAMELEN   
  USE mo_art_prescribed_types,          ONLY: t_art_emiss_prescribed
  USE mo_art_mode_fields,               ONLY: t_mode_fields, t_art_map2tracer
  USE mo_art_modes_linked_list,         ONLY: t_mode_list, append_mode

  IMPLICIT NONE

  PUBLIC :: t_art_emiss_type_container
  PUBLIC :: t_art_emiss_storage
  PUBLIC :: t_art_aero_emiss
  PUBLIC :: t_art_emiss2tracer

  PRIVATE

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  TYPE t_art_emiss_standard
    INTEGER   :: &
      &  mode,   &            !< mode describing unit of the standard emission: 
                              !  1 = kg m-2 s-1; 2 = mol mol-1 s-1
      &  num_emiss_lev        !< number of lowest model levels into which emission 
                              !  shall be included
    REAL(wp)  ::  &
      &  val                  !< value of the standard emission
  END TYPE t_art_emiss_standard
!!
!!-------------------------------------------------------------------------
!!
  TYPE t_art_emiss_bioonl
    INTEGER ::  &
      &  num_emiss_lev,  &    !< number of lowest model levels into which emission 
                              !  shall be included
      &  idx                  !< internal tracer index in mo_art_bvoc_guenther2012
    REAL(wp), POINTER  :: &
      &  pft(:,:,:)           !< pointer to the pft data (has only to be saved once)
    REAL(wp), ALLOCATABLE :: &
      &  cbio(:,:)            !< online biogenic emission (in kg m-2 s-1)
    REAL(wp) :: &
      &  scaling_factor       !< scaling factor with which the emission is multipled (usually 1)
  END TYPE t_art_emiss_bioonl
!!
!!-------------------------------------------------------------------------
!!
  TYPE t_art_emiss_type_container
    TYPE(t_art_emiss_prescribed), ALLOCATABLE  ::  &
      &  types(:)             !< prescribed emission to be considered for the tracer
    INTEGER  :: &
      &  num_types_prescribed !< number of types of emission
    TYPE(t_art_emiss_standard)  :: &
      &  std                  !< standard emission
    TYPE(t_art_emiss_bioonl)  :: &
      &  bioonl               !< online biogenic emission
  END TYPE t_art_emiss_type_container
!!
!!-------------------------------------------------------------------------
!!
  TYPE t_art_emisstorage_element
    INTEGER            ::  &
     &   tracer_idx           !< ICON index of the tracer
    TYPE(t_art_emiss_type_container) ::  &
     &  emiss                 !< emission metadata for the tracer
    LOGICAL            ::  &
     &  loccupied             !< flag if the element in the storage is occupied
  END TYPE t_art_emisstorage_element
!!
!!-------------------------------------------------------------------------
!!
  TYPE t_art_emiss_storage
    TYPE(t_art_emisstorage_element),PRIVATE, ALLOCATABLE :: &
      &  elem(:)              !< "list" of the elements
    INTEGER :: nelements = 0  !< number of elements in the storage
    LOGICAL :: is_init = .FALSE. !< flag if the storage is initialised
    CONTAINS
      PROCEDURE :: init => init_table
      PROCEDURE :: next_free_element => get_next_free_element
      PROCEDURE :: free => free_table
      PROCEDURE, PRIVATE :: put_element_emiss
      GENERIC,   PUBLIC  :: put => put_element_emiss
      PROCEDURE, PRIVATE :: update_element_emiss
      GENERIC,   PUBLIC  :: update => update_element_emiss
      PROCEDURE, PRIVATE :: get_element_emiss
      GENERIC,   PUBLIC  :: get => get_element_emiss
      PROCEDURE, PRIVATE :: check_element_emiss
      GENERIC,   PUBLIC  :: exists => check_element_emiss
  END TYPE t_art_emiss_storage
  !!
  !!-------------------------------------------------------------------------
  !!
  TYPE, extends(t_art_map2tracer) :: t_art_emiss2tracer 
    !CHARACTER(LEN=100)    :: &
    !  &  name     = ''         !< Name of emission routine <- part of base class
    LOGICAL               :: &
      &  lcalcemiss = .FALSE.  !< Calculate emissions for this parameterization
    REAL(wp)              :: &
      &  rho                   !< Density
    REAL(wp), ALLOCATABLE :: &
      &  dg0(:,:),           & !< Emission median diameter with respect to 0th moment (nmodes,1)
      &  dg3(:,:),           & !< Emission median diameter with respect to 3th moment (nmodes,1)
      &  sigmag(:,:),        & !< Standard deviation (nmodes,ntr)
      &  molweight(:,:),     & !< Molar mass mass (nmodes,ntr)
      &  weight(:,:)           !< Scaling factor for sharing emission fluxes (nmodes,ntr)
    CONTAINS
      PROCEDURE :: create => create_t_art_emiss2tracer
      PROCEDURE :: set_meta => set_meta_t_art_emiss2tracer
      PROCEDURE :: distribute_emissions => distribute_emissions
      PROCEDURE :: calc_weights => calc_weights
      PROCEDURE :: calc_number_from_mass => calc_number_from_mass
      PROCEDURE :: calc_tot_ext => calc_tot_ext
  END TYPE t_art_emiss2tracer
!!
!!-------------------------------------------------------------------------
!!
  TYPE t_art_aero_emiss
    LOGICAL                   :: &
      &  lisinit = .FALSE.                  !< set TRUE if cart_aero_emiss_xml is given
    TYPE(t_mode_list),POINTER :: e2t_list   !< emission2tracer maps
  END TYPE t_art_aero_emiss

  CHARACTER(LEN=*), PARAMETER :: routine = 'mo_art_emiss_types'

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE create_t_art_emiss2tracer(this_fields,modename,idims)
!<
! SUBROUTINE create_t_art_emiss2tracer
! This subroutine performs a setup for the fields of type t_art_emiss2tracer
! Part of Module: mo_art_emiss_types
! Author: Sven Werchner, KIT (based on: create_t_mode_fields)
! Initial Release: 2018-11-23
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_art_emiss2tracer),INTENT(INOUT) :: &
    &  this_fields                           !< emiss2tracer-structre
  CHARACTER(LEN=*),INTENT(IN)        :: &
    &  modename                           !< Name of emission routine
  INTEGER,INTENT(IN)                 :: &
    &  idims(3)

  this_fields%linit = .FALSE.
  IF(LEN_TRIM(modename) > IART_VARNAMELEN) THEN
    CALL finish(TRIM(routine)//':create_t_art_emiss2tracer',          &
      &         'Length of name for '//TRIM(modename)//' to long.')
  ENDIF
  this_fields%name  = TRIM(modename)

END SUBROUTINE create_t_art_emiss2tracer
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE set_meta_t_art_emiss2tracer(this_fields, meta_key_value_store)
!<
! SUBROUTINE set_meta_t_art_emiss2tracer
! This subroutine sets the metadata associated to type t_art_emiss2tracer
! Part of Module: mo_art_emiss_types
! Author: Sven Werchner (KIT)
! Initial Release: 2018-11-23
! Modifications:
! YYYY-MM-DD: <Name>, <Institution>
! - <Description>
!>
  CLASS(t_art_emiss2tracer),INTENT(INOUT) :: &
    &  this_fields                        !< Container with fields
  TYPE(t_key_value_store),INTENT(IN)      :: &
    &  meta_key_value_store               !< Metadata container
  INTEGER                            :: &
    &  ierror,                          &  !< Error check for key_value_store
    &  imodes,                          &  !< loop variable
    &  isub,                            &  !< loop variable
    &  nsub                                !< number of substances
  REAL(wp)                           :: &
    &  dg0,dg3,sigmag
  CHARACTER(LEN=IART_VARNAMELEN) :: &
    &  imodes_str,                  &
    &  substances                         !< comma separated string of substances
  CHARACTER(:),ALLOCATABLE       :: &
    &  c_tmp
  INTEGER, ALLOCATABLE ::               &
    &  starts(:),                       & !< start-indices of individual substances in 
                                          !   substancesstring
    &  lens(:)                            !< lengths of individual substances in substancesstring

  this_fields%linit = .TRUE.
  ierror         = SUCCESS
  
  CALL meta_key_value_store%get('nmodes',this_fields%nmodes,ierror)
  IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':set_meta_t_art_emiss2tracer', &
      &                                'Metadata nmodes not available in key_value_store.')
      
  CALL meta_key_value_store%get('rho',this_fields%rho,ierror)
  IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':set_meta_t_art_emiss2tracer',&
                          &         'rho not found in XML file.')
      
  !ALLOCATE
  IF(.NOT.ALLOCATED(this_fields%modenames)) ALLOCATE(this_fields%modenames(this_fields%nmodes))
  IF(.NOT.ALLOCATED(this_fields%ntrpermode)) ALLOCATE(this_fields%ntrpermode(this_fields%nmodes))
  IF(.NOT.ALLOCATED(this_fields%dg0   )) ALLOCATE(this_fields%dg0   (this_fields%nmodes,1))
  IF(.NOT.ALLOCATED(this_fields%dg3   )) ALLOCATE(this_fields%dg3   (this_fields%nmodes,1))
  IF(.NOT.ALLOCATED(this_fields%sigmag)) ALLOCATE(this_fields%sigmag(this_fields%nmodes,1))
  
  DO imodes = 1, this_fields%nmodes
    WRITE(imodes_str,'(I3)') imodes
    CALL meta_key_value_store%get('d_g0_'//TRIM(ADJUSTL(imodes_str)),dg0,ierror)
    IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':set_meta_t_art_emiss2tracer',&
                           &         'd_g0_'//TRIM(ADJUSTL(imodes_str))//' not found in XML file.')
    CALL meta_key_value_store%get('d_g3_'//TRIM(ADJUSTL(imodes_str)),dg3,ierror)
    IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':set_meta_t_art_emiss2tracer',&
                           &         'd_g3_'//TRIM(ADJUSTL(imodes_str))//' not found in XML file.')
    CALL meta_key_value_store%get('sigma_g_'//TRIM(ADJUSTL(imodes_str)),sigmag,ierror)
    IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':set_meta_t_art_emiss2tracer',&
                           &      'sigma_g_'//TRIM(ADJUSTL(imodes_str))//' not found in XML file.')
    this_fields%dg0   (imodes,1) = dg0
    this_fields%dg3   (imodes,1) = dg3
    this_fields%sigmag(imodes,1) = sigmag
    this_fields%modenames(imodes)= ""
    this_fields%ntrpermode(imodes)=0
  ENDDO ! imodes

!  ! reset this_fields%nmodes to 0 for dynamic determination in mapping process
!  this_fields%nmodes = 0
  
  ! Substances
  CALL meta_key_value_store%get('substances',c_tmp,ierror)
  IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':set_meta_t_art_emiss2tracer',&
                          &         'substances not found in XML file.')
  WRITE(substances,'(A)') TRIM(c_tmp)
  !split substance-string into array of substances
  !maximum number of separated strings in string of length N is
  !int((N+1)/2)
  nsub=FLOOR((LEN_TRIM(substances)+1)/2._wp)
  ALLOCATE(starts(nsub),lens(nsub))
  nsub=1
  CALL split_string(substances,nsub,starts,lens)
  this_fields%nsub=nsub
  IF(.NOT.ALLOCATED(this_fields%substance)) ALLOCATE(this_fields%substance(nsub))
  
  DO isub = 1,nsub
    this_fields%substance(isub) = substances(starts(isub):starts(isub)+lens(isub)-1)
  END DO
  DEALLOCATE(starts,lens)

END SUBROUTINE set_meta_t_art_emiss2tracer
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE distribute_emissions(this, emiss_rate_mass, emiss_rate_numb, tracer, &
             &                  rho, dtime, istart, iend, kstart, nlev, jb)
!<
! SUBROUTINE distribute_emissions
! This subroutine updates mass- and number concentrations 
! with emission fluxes based on weights of molar masses
! Part of Module: mo_art_emiss_types
! Based on: COSMO-ART code
! Author: Simon Gruber, KIT
! Initial Release: 2018-10-08
!>
  CLASS(t_art_emiss2tracer), INTENT(INOUT) :: &
    &  this  
  INTEGER, INTENT(IN)   :: &
    &  istart, iend,       & !< Start and end indices of loops
    &  kstart, nlev,       & !<
    &  jb
  REAL(wp), INTENT(IN)  :: &
    &  emiss_rate_mass(istart:iend,1:nlev,1:this%nmodes), & !< Mass emission rate   [kg m-3 s-1]
    &  emiss_rate_numb(istart:iend,1:nlev,1:this%nmodes), & !< Number emission rate [m-3 s-1]
    &  rho(:,:),           & !< Density of air [kg/m3]
    &  dtime                 !< Integration time step [s]
  REAL(wp), INTENT(INOUT) :: &
    &  tracer(:,:,:,:)       !< Tracer mixing ratios [kg kg-1]
! Local variables
  REAL(wp)              :: &
    &  conv_fac              !< factor for conversion
  INTEGER               :: &
    &  jc, jk, itr, imod     !< loop indices

  DO imod = 1, this%nmodes
    DO itr = 1, this%ntr
      DO jc = istart, iend
        DO jk = kstart, nlev
            conv_fac = dtime/rho(jc,jk)*this%weight(imod,itr)
            tracer(jc,jk,jb,this%itr3(imod,itr)) = tracer(jc,jk,jb,this%itr3(imod,itr)) & 
              &                                  + emiss_rate_mass(jc,jk,imod)*conv_fac
            tracer(jc,jk,jb,this%itr0(imod)) = tracer(jc,jk,jb,this%itr0(imod)) & 
              &                              + emiss_rate_numb(jc,jk,imod)*conv_fac
        ENDDO ! jk
      ENDDO ! jc
    ENDDO ! itr
  ENDDO ! imod

END SUBROUTINE distribute_emissions
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE calc_weights(this)
!<    
! SUBROUTINE calc_weights
! This subroutine calculates weights for distribution of the emission fluxes
! based on molar mass of species considered 
! Part of Module: mo_art_emiss_types
! Based on: COSMO-ART code
! Author: Simon Gruber, KIT
! Initial Release: 2018-10-08
!>
  CLASS(t_art_emiss2tracer), INTENT(INOUT) :: &
    &  this  
! Local variables            
  REAL(wp)              :: &
    &  sum_weights           !< sum of molar masses
  INTEGER               :: &
    &  imod,               & !< loop index
    &  itr                   !< loop index

  DO imod = 1, this%nmodes
    sum_weights = SUM(this%molweight(imod,:))
    DO itr = 1, this%ntr
      this%weight(imod,itr) = this%molweight(imod,itr)/MAX(sum_weights,1.0e-20_wp)
    ENDDO
  ENDDO

END SUBROUTINE calc_weights
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE calc_number_from_mass(this, mass, number, istart, iend, kstart, nlev)
!<
! SUBROUTINE calc_number_from_mass
! This subroutine calculates a number number rate 
! from a mass rate
! Part of Module: mo_art_emiss_types
! Based on: COSMO-ART code
! Author: Daniel Rieger, KIT
! Initial Release: 2017-04-20
!>
  CLASS(t_art_emiss2tracer), INTENT(INOUT) :: &
    &  this  
  INTEGER, INTENT(IN)   :: &
    &  istart, iend,       & !< Start and end indices of loops
    &  kstart, nlev          !<
  REAL(wp), INTENT(IN)  :: &
    &  mass(istart:iend,1:nlev,1:this%nmodes)   !< Mass to be converted from
  REAL(wp), INTENT(INOUT) :: &
    &  number(istart:iend,1:nlev,1:this%nmodes) !< Number to be converted to
! Local variables
  REAL(wp)              :: &
    &  factnum,            & !< Factor to convert 3rd moment to 0th moment
    &  fac                   !< Factor to convert 3rd moment to mass
  INTEGER               :: &
    &  jc, jk,             & !< loop index
    &  imod                  !< loop index

  DO imod = 1, this%nmodes
    factnum = EXP(4.5_wp*(LOG(this%sigmag(imod,1))**2))/(this%dg3(imod,1)**3)
    fac     = 1.0e-9_wp*6.0_wp/pi/this%rho
    DO jc = istart, iend
      DO jk = kstart, nlev
        number(jc,jk,imod) = factnum*fac*mass(jc,jk,imod)
      ENDDO ! jk
    ENDDO ! jc
  ENDDO ! imod

END SUBROUTINE calc_number_from_mass
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE calc_tot_ext(this, istart, iend, kstart, nlev, ext, tracer, rho, ext_m1, ext_m2, ext_m3)
!<
! SUBROUTINE calc_tot_ext
! This subroutine calculates the total extinction 
! coefficient of tracers contained within this species
! Part of Module: mo_art_emiss_types
! Based on: 
! Author: Simon Gruber, KIT
! Initial Release: 2018-10-26
!>
  CLASS(t_art_emiss2tracer), INTENT(INOUT) :: &
    &  this  
  REAL(wp), INTENT(INOUT) :: &
    &  ext(:,:)            !< Total extinction coefficient
  REAL(wp), INTENT(IN)  :: &
    &  rho(:,:),           & !< Density of air
    &  tracer(:,:,:)         !< Tracer mass mixing ratios
  REAL(wp), TARGET, INTENT(IN)  :: &
    &  ext_m1,          & !< Extinction coefficients
    &  ext_m2,          & !  from Mie calculations
    &  ext_m3             !  for one wavelength
  INTEGER, INTENT(IN)   :: &
    &  istart, iend,       & !< Start and end indices of loops
    &  kstart, nlev          !<
! Local variables
  REAL(wp), POINTER     :: &
    &  ext_ptr            !< extinction coefficients pointer
  INTEGER               :: &
    &  jc, jk,             & !< loop index
    &  imod, itr             !< loop index

  ext(:,:) = 0._wp

  DO imod = 1, this%nmodes
    SELECT CASE(imod)
      CASE (1)
        ext_ptr => ext_m1
      CASE (2)
        ext_ptr => ext_m2
      CASE (3)
        ext_ptr => ext_m3
    END SELECT
    DO itr  = 1, this%ntr
      DO jk = kstart, nlev
        DO jc = istart, iend
          ext(jc,jk) = ext(jc,jk) + ext_ptr                &
            &          * tracer(jc,jk,this%itr3(imod,itr)) &
            &          * rho(jc,jk)*1.0e-6_wp
        ENDDO ! jc
      ENDDO ! jk
    ENDDO ! itr
  ENDDO ! imod

!!        DO i_wavel = 1, 9
!!          ext_dust(i_wavel,jc,jk) = ext_dust(i_wavel,jc,jk) &
!!            &                     + ext_dust_ptr(i_wavel)*tracer(jc,jk,this%itr3(imod,itr)) & 
!!            &                     * rho(jc,jk)*1.0e-6_wp
!!        ENDDO ! i_wavel

END SUBROUTINE calc_tot_ext
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE init_table(this_storage, nelements)
!<
! SUBROUTINE init_table
! This subroutine initializes a t_art_emiss_storage table
! Part of Module: mo_art_emiss_types
! Author: Michael Weimer, KIT
! Initial Release: 2016-11-02
! Modifications:
! YYYY-MM-DD: <name>,<institution>
! - <description>
!>
  CLASS(t_art_emiss_storage),INTENT(INOUT) :: &
    &  this_storage
  INTEGER,INTENT(IN) ::                       &
    &  nelements
  INTEGER ::                                  &
    &  nelem

  IF (.NOT.this_storage%is_init .AND. .NOT.ALLOCATED(this_storage%elem)) THEN
    ALLOCATE(this_storage%elem(nelements))
    this_storage%is_init   = .TRUE.
    this_storage%nelements = nelements
    DO nelem = 1, this_storage%nelements
      this_storage%elem(nelem)%tracer_idx = UNDEF_INT_ART
      this_storage%elem(nelem)%loccupied = .FALSE.
      this_storage%elem(nelem)%emiss%num_types_prescribed = 0
    ENDDO
  ELSE
    CALL finish (TRIM(routine)//':init_table','Storage is already initialized.')
  ENDIF
END SUBROUTINE init_table
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE get_next_free_element(this_storage, tracer_idx, next_free_element)
!<
! SUBROUTINE get_next_free_element
! This subroutine finds the next free element in a t_art_emiss_storage table
! Part of Module: mo_art_emiss_types
! Author: Michael Weimer, KIT
! Initial Release: 2016-11-02
! Modifications:
! YYYY-MM-DD: <name>,<institution>
! - <description>
!>
  CLASS(t_art_emiss_storage),INTENT(INOUT) :: &
    &  this_storage
  INTEGER,INTENT(IN)              :: &
    &  tracer_idx
  INTEGER,INTENT(INOUT)                      :: &
    &  next_free_element
! Local variables
  INTEGER                                  :: &
    &  nelem
  CHARACTER(LEN = 5)   :: &
    &  tracer_idx_str

  next_free_element = -1

  IF (.NOT.this_storage%is_init .OR. .NOT.ALLOCATED(this_storage%elem)) THEN
    CALL finish (TRIM(routine)//':get_next_free_element','Storage is not initialized.')
  ENDIF

  IF (this_storage%nelements < 1) THEN
    CALL finish (TRIM(routine)//':get_next_free_element','Storage does not contain elements.')
  ENDIF


  DO nelem = 1, this_storage%nelements
    ! Check for already existing entry
    IF(this_storage%elem(nelem)%tracer_idx == tracer_idx) THEN
      WRITE(tracer_idx_str,'(I5)') tracer_idx
      CALL finish (TRIM(routine)//':get_next_free_element',  &
        &  'Storage does already contain element for tracer index '//TRIM(tracer_idx_str)//'.')
    ENDIF

    IF (.NOT.this_storage%elem(nelem)%loccupied) THEN
      next_free_element = nelem
      EXIT
    ENDIF
  ENDDO

  IF (next_free_element == -1) THEN
    CALL finish (TRIM(routine)//':get_next_free_element','Storage full.')
  ENDIF
END SUBROUTINE get_next_free_element
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE free_table(this_storage)
!<
! SUBROUTINE free_table
! This subroutine sets a t_art_emiss_storage table free
! Part of Module: mo_art_emiss_types
! Author: Michael Weimer, KIT
! Initial Release: 2016-11-02
! Modifications:
! YYYY-MM-DD: <name>,<institution>
! - <description>
!>
  CLASS(t_art_emiss_storage),INTENT(INOUT), TARGET :: &
    &  this_storage
! Local variables
  INTEGER                                  :: &
    &  nelem, j, i_dt
  TYPE(t_art_emiss_type_container), POINTER ::  &
    & emiss

  IF (ALLOCATED(this_storage%elem)) THEN

    DO nelem = 1, this_storage%nelements
      this_storage%elem(nelem)%tracer_idx = -1
      emiss => this_storage%elem(nelem)%emiss
      emiss%std%val = 0.0_wp
      emiss%std%num_emiss_lev = UNDEF_INT_ART
      emiss%std%mode = UNDEF_INT_ART
      IF (emiss%num_types_prescribed > 0) THEN
        DO j = 1,emiss%num_types_prescribed
          IF (ALLOCATED(emiss%types(j)%emiss_2d)) DEALLOCATE(emiss%types(j)%emiss_2d)
          
          DO i_dt = 1,3,2
            CALL deallocateDatetime(emiss%types(j)%datetime(i_dt)%ptr)
          END DO
          CALL deallocateDatetime(emiss%types(j)%first_date)
          CALL deallocateDatetime(emiss%types(j)%last_date)
        END DO
      END IF

      IF (ALLOCATED(emiss%types)) DEALLOCATE(emiss%types)
      emiss%num_types_prescribed = UNDEF_INT_ART

      IF (ALLOCATED(emiss%bioonl%cbio)) DEALLOCATE(emiss%bioonl%cbio)
      
      this_storage%elem(nelem)%loccupied = .FALSE.
    ENDDO

    DEALLOCATE(this_storage%elem)
    this_storage%is_init   = .FALSE.
    this_storage%nelements = 0
  ENDIF
END SUBROUTINE free_table
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE put_element_emiss(this_storage, tracer_idx, emiss)
!<
! SUBROUTINE put_element_real
! This subroutine puts an element into a t_art_emiss_storage table
! Part of Module: mo_art_emiss_types
! Author: Michael Weimer, KIT
! Initial Release: 2016-11-02
! Modifications:
! YYYY-MM-DD: <name>,<institution>
! - <description>
!>
  CLASS(t_art_emiss_storage),INTENT(INOUT) :: &
    &  this_storage
  INTEGER ,INTENT(IN)              :: &
    &  tracer_idx
  TYPE(t_art_emiss_type_container),INTENT(IN)                      :: &
    &  emiss
! Local variables
  INTEGER                                  :: &
    &  id_next_free

  IF (.NOT.this_storage%is_init) THEN
    CALL finish (TRIM(routine)//':put_element_emiss','Storage was not initialized.')
  ENDIF

  CALL this_storage%next_free_element(tracer_idx,id_next_free)

  this_storage%elem(id_next_free)%tracer_idx = tracer_idx
  this_storage%elem(id_next_free)%emiss      = emiss
  this_storage%elem(id_next_free)%loccupied  = .TRUE.
END SUBROUTINE put_element_emiss
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE get_element_emiss(this_storage,tracer_idx , emiss, ierror)
!<
! SUBROUTINE get_element_real
! This subroutine gets an element from a t_art_emiss_storage table
! Part of Module: mo_art_emiss_types
! Author: Michael Weimer, KIT
! Initial Release: 2016-11-02
! Modifications:
! YYYY-MM-DD: <name>,<institution>
! - <description>
!>
  CLASS(t_art_emiss_storage),INTENT(IN)    :: &
    &  this_storage
  INTEGER,INTENT(IN)              :: &
    &  tracer_idx
  TYPE(t_art_emiss_type_container),INTENT(INOUT)                     :: &
    &  emiss
  INTEGER,INTENT(INOUT), OPTIONAL            :: &
    &  ierror
! Local variables
  INTEGER                                  :: &
    &  nelem

  IF (PRESENT(ierror)) ierror = 0

  IF (.NOT.this_storage%is_init) THEN
    CALL finish (TRIM(routine)//':get_element_emiss','Storage was not initialized.')
  ENDIF

  IF (this_storage%nelements < 1) THEN
    CALL finish (TRIM(routine)//':get_element_emiss','Storage does not contain elements.')
  ENDIF

  DO nelem = 1, this_storage%nelements
    IF (tracer_idx == this_storage%elem(nelem)%tracer_idx) THEN
      emiss = this_storage%elem(nelem)%emiss
      IF (.NOT. this_storage%elem(nelem)%loccupied) THEN
        CALL finish (TRIM(routine)//':get_element_emiss','Storage inconsistency.')
      ENDIF
      RETURN
    ENDIF
  ENDDO

  IF (PRESENT(ierror)) ierror = 1
END SUBROUTINE get_element_emiss
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE update_element_emiss(this_storage,tracer_idx , emiss, ierror)
!<
! SUBROUTINE update_element_emiss
! This subroutine updates an element from a t_art_emiss_storage table
! Part of Module: mo_art_emiss_types
! Author: Michael Weimer, KIT
! Initial Release: 2016-11-25
! Modifications:
! YYYY-MM-DD: <name>,<institution>
! - <description>
!>
  CLASS(t_art_emiss_storage),INTENT(INOUT)    :: &
    &  this_storage
  INTEGER,INTENT(IN)              :: &
    &  tracer_idx
  TYPE(t_art_emiss_type_container),INTENT(IN)                     :: &
    &  emiss
  INTEGER,INTENT(INOUT), OPTIONAL            :: &
    &  ierror
! Local variables
  INTEGER                                  :: &
    &  nelem

  IF (PRESENT(ierror)) ierror = 0

  IF (.NOT.this_storage%is_init) THEN
    CALL finish (TRIM(routine)//':update_element_emiss','Storage was not initialized.')
  ENDIF

  IF (this_storage%nelements < 1) THEN
    CALL finish (TRIM(routine)//':update_element_emiss','Storage does not contain elements.')
  ENDIF

  DO nelem = 1, this_storage%nelements
    IF (tracer_idx == this_storage%elem(nelem)%tracer_idx) THEN
      this_storage%elem(nelem)%emiss = emiss
      IF (.NOT. this_storage%elem(nelem)%loccupied) THEN
        CALL finish (TRIM(routine)//':update_element_emiss','Storage inconsistency.')
      ENDIF
      RETURN
    ENDIF
  ENDDO

  IF (PRESENT(ierror)) ierror = 1
END SUBROUTINE update_element_emiss
!!
!!-------------------------------------------------------------------------
!!
LOGICAL FUNCTION check_element_emiss(this_storage,tracer_idx)
!<
! FUNCITON check_element_emiss
! This subroutine check if the element with the given tracer index exists
! Part of Module: mo_art_emiss_types
! Author: Michael Weimer, KIT
! Initial Release: 2016-11-29
! Modifications:
! YYYY-MM-DD: <name>,<institution>
! - <description>
!>
  CLASS(t_art_emiss_storage),INTENT(IN)    :: &
    &  this_storage
  INTEGER,INTENT(IN)              :: &
    &  tracer_idx
! Local variables
  INTEGER                                  :: &
    &  nelem

  IF (.NOT.this_storage%is_init) THEN
    check_element_emiss = .FALSE.
    RETURN
  ENDIF

  IF (this_storage%nelements < 1) THEN
    check_element_emiss = .FALSE.
    RETURN
  ENDIF

  check_element_emiss = .FALSE.

  DO nelem = 1, this_storage%nelements
    IF (tracer_idx == this_storage%elem(nelem)%tracer_idx) THEN
      check_element_emiss = .TRUE.
      IF (.NOT. this_storage%elem(nelem)%loccupied) THEN
        CALL finish (TRIM(routine)//':check_element_emiss','Storage inconsistency.')
      ENDIF
      RETURN
    ENDIF
  ENDDO

END FUNCTION check_element_emiss
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_emiss_types
