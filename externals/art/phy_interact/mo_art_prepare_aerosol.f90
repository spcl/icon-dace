!
! mo_art_prepare_aerosol
! Prepares mineral dust aerosol concentration for the KL06 scheme
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

MODULE mo_art_prepare_aerosol
  USE mo_kind,                          ONLY: wp
  USE mo_math_constants,                ONLY: pi
  USE mo_impl_constants,                ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_exception,                     ONLY: finish
! ART ROUTINES
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_emiss_types,               ONLY: t_art_emiss2tracer
  USE mo_art_diag_types,                ONLY: t_art_diag
  USE mo_art_modes,                     ONLY: t_fields_2mom, t_fields_1mom
  USE mo_art_modes_linked_list,         ONLY: p_mode_state,t_mode
  USE mo_art_nucleation_interface,      ONLY: t_nuc_mode_cold
  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_prepare_dust_kl06, art_prepare_dust_inas

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_prepare_dust_kl06(jg, jb, istart, iend, kstart, kend, rho, tracer)
!<
! SUBROUTINE art_prepare_dust_kl06
! Prepares mineral dust aerosol concentration for the KL06 scheme
! Based on: -
! Part of Module: mo_art_prepare_aerosol
! Author: Daniel Rieger, KIT
! Initial Release: 2014-11-10
! Modifications:
!! YYYY-MM-DD: <name>, <institution>
!! - ...
!>
  INTEGER,  INTENT(IN)            :: &
    &  jg, jb,                       & !< Domain and block index
    &  istart, iend,                 & !< Start/end index jc loop
    &  kstart, kend                    !< Start/end index jk loop
  REAL(wp), INTENT(in), TARGET    :: &
    &  rho(:,:)                        !< Density
  REAL(wp), INTENT(INOUT), TARGET :: &
    &  tracer(:,:,:)                   !< Tracer fields
! Local Variables
  INTEGER                         :: &
    &  jc, jk,                       & !< Loop indices
    &  imod, itr                       !< Loop indices
  CHARACTER(LEN=MAX_CHAR_LENGTH)  :: &
    &  thisroutine = "mo_art_prepare_aerosol:art_prepare_dust_kl06"
  TYPE(t_mode), POINTER  :: &
    & current

  IF (p_art_data(jg)%tracer2aeroemiss%lisinit) THEN
    current=>p_art_data(jg)%tracer2aeroemiss%e2t_list%p%first_mode
    DO WHILE(ASSOCIATED(current))
      SELECT TYPE(this=>current%fields)
        TYPE IS(t_art_emiss2tracer)
          SELECT CASE(this%name)
            CASE('dust')
              IF(this%lcalcemiss) THEN
                p_art_data(jg)%diag%ndust_tot(:,:,jb) = 0._wp
                DO imod = 1, this%nmodes
                  DO jk = kstart, kend
                    DO jc = istart, iend
                      p_art_data(jg)%diag%ndust_tot(jc,jk,jb) =                                     &
                        &                                   p_art_data(jg)%diag%ndust_tot(jc,jk,jb) &
                        &                                 + tracer(jc,jk,this%itr0(imod))*rho(jc,jk)
                    ENDDO !jc
                  ENDDO !jk
                ENDDO ! imod
              ELSE
                CALL finish (thisroutine, 'no dust emissions found.')
              ENDIF ! this%lcalcemiss
            CASE DEFAULT
              !nothing to do
          END SELECT
      END SELECT
      current=>current%next_mode
    END DO
  END IF
  
END SUBROUTINE art_prepare_dust_kl06
!!
!!-------------------------------------------------------------------------
!!
!! Preparation for INAS-based ice nucleation
!!
SUBROUTINE art_prepare_dust_inas(jg, jb, istart, iend, kstart, kend, rho, p_trac, &
                                 numdust, sfcdust, aod_crit)
!>
  INTEGER,  INTENT(IN)            :: &
    &  jg, jb,                       & !< Domain and block index
    &  istart, iend,                 & !< Start/end index jc loop
    &  kstart, kend                    !< Start/end index jk loop
  REAL(wp), INTENT(in), TARGET    :: &
    &  rho(:,:)                        !< Density
  REAL(wp), INTENT(inout), TARGET :: &
    &  p_trac(:,:,:)                   !< Tracer fields
  REAL(wp), INTENT(inout), TARGET :: &
    &  numdust(:,:),                 & !< Number density of dust
    &  sfcdust(:,:)                    !< Surface area of dust
  REAL(wp), INTENT(in)            :: &
    &  aod_crit                        !< Threshold for dust AOD

  ! Local Variables
  TYPE(t_mode), POINTER           :: this_mode
  TYPE(t_art_diag),POINTER        :: &
    &  art_diag                      !< Pointer to ART diagnostic fields
  INTEGER                 ::       &
    &  i, k,                       & !< Loop indizes
    &  imodes,                     & !< Counter of modes
    &  idx                           !< Indexing variable
  REAL(wp) ::               &
    &  ndust,               &  ! Number density of dust mode
    &  ddust,               &  ! Mean diameter of dust mode
    &  sdust,               &  ! Surface area of dust mode
    &  sigdust                 ! Standard deviations of mineral dust distributions

  CHARACTER(len=*), PARAMETER :: routine = 'art_prepare_dust_inas'

  TYPE(t_nuc_mode_cold) :: &
       &  nuc_modes_cold(p_mode_state(jg)%p_mode_list%p%nmodes)

  ! Associate Pointer for short references
  art_diag => p_art_data(jg)%diag

  ! ----------------------------------
  ! --- Start routine
  ! ----------------------------------

  ! Loop through modes

  imodes  = 0
  this_mode => p_mode_state(jg)%p_mode_list%p%first_mode
  DO WHILE(ASSOCIATED(this_mode))
    imodes = imodes + 1
    ! Select type of mode // This may be not optimal to do a select type within a jb loop.
    ! But the object nuc_modes_cold cant be allocated outside the two-mom scheme and we would
    ! then suffer the problem that it needs the dimension jg
    ! ->Maybe a better solution can be found
    SELECT TYPE (fields=>this_mode%fields)
      CLASS is (t_fields_2mom)
        ! Calculate modal parameters (i.e. diameter, 3rd moment)
        CALL fields%modal_param(p_art_data(jg)%air_prop%art_free_path(:,:,jb),                &
          &                     istart, iend, kstart, kend, jb, p_trac(:,:,:))
        ! Drieg: do_homfreez/do_hetnuc needs to be set via Metadata in the future!!!
        nuc_modes_cold(imodes)%do_homfreez = .TRUE.
        nuc_modes_cold(imodes)%is_dust     = .FALSE.
        nuc_modes_cold(imodes)%is_soot     = .FALSE.
        nuc_modes_cold(imodes)%is_org      = .FALSE.
        ! Exceptions (Not necessary if using Metadata in the future)
        IF (TRIM(fields%name) == 'dusta') THEN
          nuc_modes_cold(imodes)%do_homfreez = .FALSE.
          nuc_modes_cold(imodes)%is_dust     = .TRUE.
        ENDIF
        IF (TRIM(fields%name) == 'dustb') THEN
          nuc_modes_cold(imodes)%do_homfreez = .FALSE.
          nuc_modes_cold(imodes)%is_dust     = .TRUE.
        ENDIF
        IF (TRIM(fields%name) == 'dustc') THEN
          nuc_modes_cold(imodes)%do_homfreez = .FALSE.
          nuc_modes_cold(imodes)%is_dust     = .TRUE.
        ENDIF
        nuc_modes_cold(imodes)%i_numb   = fields%itr0
        nuc_modes_cold(imodes)%sigma    = fields%info%sg_ini
        nuc_modes_cold(imodes)%diam_pointer => fields%diameter(:,:,jb)
        nuc_modes_cold(imodes)%mom3_pointer => fields%third_moment(:,:,jb)
        nuc_modes_cold(imodes)%jsp(:)   = fields%itr3(:)
        nuc_modes_cold(imodes)%njsp     = fields%ntr-1
      CLASS is (t_fields_1mom)
        nuc_modes_cold(imodes)%do_homfreez = .FALSE.
        nuc_modes_cold(imodes)%is_dust     = .FALSE.
        nuc_modes_cold(imodes)%is_soot     = .FALSE.
        nuc_modes_cold(imodes)%is_org      = .FALSE.
      CLASS DEFAULT
        ! Should not happen...
        CALL finish('mo_art_prepare_aerosol: art_prepare_dust_inas', &
             &      'ART: Unknown class')
    END SELECT
    this_mode => this_mode%next_mode
  ENDDO

  ! .. calculate dust for INAS parameterization (loop over dust modes)
  idx = 1
  numdust(:,:) = 0.0_wp
  sfcdust(:,:) = 0.0_wp
  DO imodes = 1, p_mode_state(jg)%p_mode_list%p%nmodes
    IF (nuc_modes_cold(imodes)%is_dust) THEN

      IF (idx > 3) THEN
        CALL finish('mo_art_prepare_aerosols:', &
             &      'ART: More dust species than 3 found')
      ENDIF

      ! with Eq. (3.16) of Riemer (2002)
      sigdust = pi * EXP( 2._wp * LOG( nuc_modes_cold(imodes)%sigma )**2 )
      DO k = kstart, kend
        DO i = istart, iend

          IF ( p_art_data(jg)%diag%dust_aeronet(5)%tau_vi(i,jb) > aod_crit ) THEN

            ! number, diameter and surface area of dust mode
            ndust = p_trac(i,k,nuc_modes_cold(imodes)%i_numb) * rho(i,k)
            ddust = nuc_modes_cold(imodes)%diam_pointer(i,k)
            sdust = sigdust * ddust**2 * ndust

            ! INAS dust input are number and surface area
            numdust(i,k) = numdust(i,k) + ndust
            sfcdust(i,k) = sfcdust(i,k) + sdust

          END IF

        ENDDO
      ENDDO
      idx = idx+1
    ENDIF
  ENDDO
  DO k = kstart, kend
    DO i = istart, iend
      ! sfcdust is the average surface area of the dust
      sfcdust(i,k) = sfcdust(i,k)/numdust(i,k)
    ENDDO
  ENDDO

  IF (idx == 1) CALL finish('mo_art_prepare_aerosols:', 'No dust found!')

END SUBROUTINE art_prepare_dust_inas
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_prepare_aerosol
