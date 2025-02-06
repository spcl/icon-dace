!> fage (forest age) process
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>#### Contains routines info for the fage (forest age) proc
!>
MODULE mo_fage_process
#ifndef __NO_JSBACH__

   USE mo_kind,      ONLY: wp
   USE mo_exception, ONLY: message, message_text, finish

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: recalc_fract_per_age_upon_ageing, recalc_fract_per_age_proportionally, &
     & synchronise_fracts, weighted_avg_var_upon_area_movement, weighted_avg_per_canopy_var_upon_area_movement

   CHARACTER(len=*), PARAMETER :: modname = 'mo_fage_process'

 CONTAINS

  ! ====================================================================================================== !
  !
  !> recalculate the tracked age fractions [y] upon ageing
  !> i.e. shifting fractions of maximum age from one age class to the next
  !
  SUBROUTINE recalc_fract_per_age_upon_ageing( nidx, nage, fract_per_age, mean_age)

    USE mo_fage_util,             ONLY: get_mean_age
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,  INTENT(in)    ::  nidx !< vector length
    INTEGER,  INTENT(in)    ::  nage !< number of tracked ages, i.e. years
    REAL(wp), INTENT(inout), DIMENSION(nidx, nage) :: fract_per_age !< vector to track exact age fractions
    REAL(wp), INTENT(inout), DIMENSION(nidx) :: mean_age            !< vector with the mean age of the pft
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER :: i
    ! -------------------------------------------------------------------------------------------------- !

    ! --- update age
    IF (nage > 1) THEN
      ! oldest age accumulates
      fract_per_age (:,nage) = fract_per_age (:,nage) + fract_per_age (:,nage-1)

      DO i=nage-1,2,-1
        fract_per_age (:,i) = fract_per_age (:,i-1)
      END DO

      ! youngest is emptied
      fract_per_age (:,1) = 0._wp

      !recalculate mean age
      mean_age (:) = get_mean_age(fract_per_age)
    END IF

  END SUBROUTINE recalc_fract_per_age_upon_ageing


  ! ====================================================================================================== !
  !
  !> synchronise fracts to avoid numerical issues
  !
  SUBROUTINE synchronise_fracts(nidx, nage, nacs, i_ac, pftid, age_ubounds, fract, fract_per_age)

    USE mo_jsb_lcc_class,     ONLY: min_tolerated_fract_mismatch
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,  INTENT(in)    ::  nidx     !< vector length
    INTEGER,  INTENT(in)    ::  nage     !< number of tracked ages, i.e. years
    INTEGER,  INTENT(in)    ::  nacs     !< number of age classes (ac)
    INTEGER,  INTENT(in)    ::  i_ac     !< index of the considered age class
    INTEGER,  INTENT(in)    ::  pftid    !< id of this pft
    REAL(wp), INTENT(inout), DIMENSION(nacs)       :: age_ubounds        !< max age of all acs
    REAL(wp), INTENT(inout), DIMENSION(nidx)       :: fract              !< area fract of this ac
    REAL(wp), INTENT(inout), DIMENSION(nidx, nage) :: fract_per_age      !< tracked age fractions all acs
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER :: j, thisMinAge, thisMaxAge
    REAL(wp), DIMENSION(nidx) :: thisAreaSum
    CHARACTER(len=*), PARAMETER :: routine = modname//':synchronise_fracts'
    ! -------------------------------------------------------------------------------------------------- !

    IF (i_ac .GT. 1) THEN
      thisMinAge = int(age_ubounds(i_ac-1)+1)
      thisMaxAge = int(age_ubounds(i_ac))

      IF (ANY(fract_per_age(:,thisMinAge:thisMaxAge) < - EPSILON(1._wp))) THEN
        WRITE (message_text,*) 'Violation of assertion: tracked age fraction got negative for age class ', &
          & i_ac, ' of pft ', pftid
        CALL finish (routine, message_text)
      END IF
      fract_per_age(:,thisMinAge:thisMaxAge) = MAX(fract_per_age(:,thisMinAge:thisMaxAge), 0.0_wp)

      thisAreaSum = SUM(fract_per_age(:,thisMinAge:thisMaxAge), DIM=2)

      IF (ANY(ABS(thisAreaSum(:) - fract(:)) .GT. min_tolerated_fract_mismatch)) THEN
        DO j = 1,nidx
          IF (ABS(thisAreaSum(j) - fract(j)) .GT. min_tolerated_fract_mismatch) THEN
            WRITE (message_text,*) ' fract and fage sum mismatch for age class ', i_ac, &
              & ' of pft ', pftid, ' found e.g. fract ', fract(j), ' and ', thisAreaSum(j)
            CALL finish (routine, message_text)
          END IF
        END DO
      END IF

      fract(:) = thisAreaSum(:)
    ELSE
      ! First age class has only 0-1 age
      fract_per_age(:,1) = fract(:)
    END IF

  END SUBROUTINE synchronise_fracts

  ! ====================================================================================================== !
  !
  !> recalculate the tracked age fractions [y] upon area loss
  !
  SUBROUTINE recalc_fract_per_age_proportionally( nidx, nage, nacs, i_ac, this_cf_loss,    &
    &                                                     age_ubounds, fract_per_age)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,  INTENT(in)    ::  nidx    !< vector length
    INTEGER,  INTENT(in)    ::  nage    !< number of tracked ages, i.e. years
    INTEGER,  INTENT(in)    ::  nacs    !< number of age classes
    INTEGER,  INTENT(in)    ::  i_ac    !< index of the age class with lost area
    REAL(wp), INTENT(in), DIMENSION(nacs) :: age_ubounds  !< max age of all age classes
    REAL(wp), INTENT(in), DIMENSION(nidx) :: this_cf_loss !< area fraction lost for this ac
    REAL(wp), INTENT(inout), DIMENSION(nidx, nage) :: fract_per_age !< tracked age fractions
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER :: j, k, thisMinAge, thisMaxAge
    REAL(wp), DIMENSION(nidx) :: thisAreaSum

    CHARACTER(len=*), PARAMETER :: routine = modname//':recalc_fract_per_age_proportionally'
    ! -------------------------------------------------------------------------------------------------- !

    IF (i_ac .EQ. 1) RETURN ! area is only re-distributed to the first child - not from

    ! lost area is assigned distributed proportionally over all ages in this age class
    thisMinAge = int(age_ubounds(i_ac-1)+1)
    thisMaxAge = int(age_ubounds(i_ac))

    thisAreaSum(:) = SUM(fract_per_age(:,thisMinAge:thisMaxAge), DIM=2)

    IF (ANY(thisAreaSum(:) .GT. 0._wp)) THEN
      ! Assertion: sum of available area needs to be larger than the disturbed area
      IF (ANY(this_cf_loss(:) .GT. thisAreaSum(:))) THEN
          WRITE (message_text,*) 'Violation of assertion: more area disturbed than available ' &
            & // 'for age class ', i_ac
          CALL finish (routine, message_text)
      END IF

      ! Subtract proportionally from all ages within this age class
      DO k = thisMinAge,thisMaxAge
        WHERE(thisAreaSum(:) > 0.0_wp)
          fract_per_age(:,k) = fract_per_age(:,k) - (this_cf_loss(:) * (fract_per_age(:,k) / thisAreaSum(:)))
        END WHERE
      END DO
    END IF

  END SUBROUTINE recalc_fract_per_age_proportionally


  ! ====================================================================================================== !
  !
  !> recalculate a given variable according to changes in tile area
  !
  ELEMENTAL SUBROUTINE weighted_avg_var_upon_area_movement(old_target_fract, delta_fract, source_var, target_var)

    ! -------------------------------------------------------------------------------------------------- !
    REAL(wp), intent(in)    :: old_target_fract !< area fraction of the target before changing area
    REAL(wp), intent(in)    :: delta_fract      !< change in area
    REAL(wp), intent(in)    :: source_var       !< value of the variable in the source tile
    REAL(wp), intent(inout) :: target_var       !< value of the variable in the target tile
    ! -------------------------------------------------------------------------------------------------- !
     IF (delta_fract .GT. 0._wp) THEN
      target_var = ((target_var * old_target_fract) + (source_var * delta_fract)) / (old_target_fract + delta_fract)
     END IF

     RETURN

  END SUBROUTINE weighted_avg_var_upon_area_movement


  ! ====================================================================================================== !
  !
  !> recalculate a given per canopy variable according to changes in tile area
  !
  ELEMENTAL SUBROUTINE weighted_avg_per_canopy_var_upon_area_movement(old_target_fract, delta_fract, &
    & veg_fract_corr_source, veg_fract_corr_target, source_per_canopy_var, target_per_canopy_var)

    ! -------------------------------------------------------------------------------------------------- !
    REAL(wp), intent(in)    :: old_target_fract      !< area fraction of the target before changing area
    REAL(wp), intent(in)    :: delta_fract           !< change in area
    REAL(wp), intent(in)    :: veg_fract_corr_source !< Clumpiness of source vegetation
    REAL(wp), intent(in)    :: veg_fract_corr_target !< Clumpiness of target vegetation
    REAL(wp), intent(in)    :: source_per_canopy_var !< value of the per canopy variable in the source
    REAL(wp), intent(inout) :: target_per_canopy_var !< value of the per canopy variable in the target
    ! -------------------------------------------------------------------------------------------------- !
    REAL(wp) :: target_var, source_var
    ! -------------------------------------------------------------------------------------------------- !
     IF (delta_fract .GT. 0._wp) THEN

      target_var = target_per_canopy_var * veg_fract_corr_target
      source_var = source_per_canopy_var * veg_fract_corr_source

      CALL weighted_avg_var_upon_area_movement(old_target_fract, delta_fract, source_var, target_var)

      target_per_canopy_var = target_var / veg_fract_corr_target
     END IF

     RETURN

  END SUBROUTINE weighted_avg_per_canopy_var_upon_area_movement

#endif
END MODULE mo_fage_process
