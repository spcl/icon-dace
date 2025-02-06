!
! mo_art_chem_utils
! This module provides the routines to calculate speceial destruction rates and
! properties of the tracers, that cannot be captured by the general algorithm
! in mo_art_chemtracer
!
!

! ART!
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

MODULE mo_art_chem_utils
  !ICON
  USE mo_kind,                    ONLY: wp
  USE mo_math_constants,          ONLY: pi
  USE mo_physical_constants,      ONLY: avo, argas, amw, amd
  USE mo_var_list,                ONLY: t_var_list_ptr
  USE mo_var,                     ONLY: t_var
  USE mo_var_metadata_types,      ONLY: t_var_metadata_dynamic, t_var_metadata
  USE mo_tracer_metadata_types,   ONLY: t_tracer_meta
  USE mo_run_config,              ONLY: iqv, iqc, iqi

  ! ART
  USE mo_art_data,                ONLY: p_art_data
  USE mo_art_chem_data,           ONLY: t_art_chem
  USE mo_art_atmo_data,           ONLY: t_art_atmo
  USE mo_art_impl_constants,      ONLY: IART_QV, IART_QC, IART_QI
  USE mo_art_chem_types,          ONLY: t_chem_meta_param, t_chem_meta_mecca
  USE mo_art_chem_types_param,    ONLY: t_chem_meta_lt,t_chem_meta_linoz
  USE mo_art_kinetic_constants,   ONLY: art_get_CO_des_1d, art_get_CH4_des_1d
  

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_get_CO2_des
  PUBLIC :: art_get_CO_des
  PUBLIC :: art_get_CH4_des
  PUBLIC :: art_get_CH3CN_des
  PUBLIC :: art_get_H2O_des
  PUBLIC :: art_calc_col
  PUBLIC :: art_calc_vmr2Nconc
  PUBLIC :: art_convert_tracers_mmr_Nconc
  PUBLIC :: art_convert_tracers_Nconc_mmr
  PUBLIC :: art_convert_tracers_vmr_mmr

  CONTAINS

SUBROUTINE art_get_CO2_des(jg,jb,jcs,jce ,&
            &             nproma,nlev, tracer )
!<
! SUBROUTINE art_get_CO2_des
! Calculates the destruction of CO2
! Based on:
! - KASIMA
! - Brasseur and Solomon?
! Part of Module: mo_art_chem_utils
! Author: Christian Stassen and  Michael Weimer, KIT
! Initial Release: 2018-10-18
! Modifications:
!>
  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &
    &  jg,                & !< patch id
    &  jb,jcs,jce           !< loop indices
  INTEGER, INTENT(in) ::  &
    &  nproma,nlev   !< dimensions
  TYPE(t_chem_meta_lt), INTENT(inout) ::   &
    &  tracer        !< lifetime tracer structure

  ! local variables
  TYPE(t_art_chem),POINTER    :: &
    &  art_chem     !< Pointer to ART chem fields
  INTEGER ::   &
    &  jc, jk       !< loop indices

  art_chem => p_art_data(jg)%chem

  DO jk = 1,nlev
    DO jc = jcs,jce
      tracer%des_3d(jc,jk,jb) = 2.20e-8_wp * EXP(-1.e-20_wp * art_chem%param%o2_column(jc,jk,jb))
    ENDDO
  ENDDO

END SUBROUTINE art_get_CO2_des

SUBROUTINE art_get_CO_des(jg,jb,jcs,jce ,&
            &             nproma,nlev, tracer )
!<
! SUBROUTINE art_get_CO_des
! Calculates the destruction of CO
! Based on:
! - KASIMA
! - Brasseur and Solomon?
! Part of Module: mo_art_chem_utils
! Author: Christian Stassen and  Michael Weimer, KIT
! Initial Release: 2018-10-18
! Modifications:
!>
  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &
    &  jg,                & !< patch id
    &  jb,jcs,jce           !< loop indices
  INTEGER, INTENT(in) ::  &
    &  nproma,nlev   !< dimensions
  TYPE(t_chem_meta_lt), INTENT(inout) ::   &
    &  tracer        !< lifetime tracer structure
  ! local variables
  TYPE(t_art_atmo),POINTER    :: &
    &  art_atmo   !< Pointer to ART chem fields
  INTEGER  :: &
    &  jc, jk  !< loop indices

  art_atmo => p_art_data(jg)%atmo

  DO jk = 1,nlev
     DO jc = jcs,jce
       CALL art_get_CO_des_1d(tracer%des_3d(jc,jk,jb),art_atmo%pres(jc,jk,jb))
    ENDDO
  ENDDO

END SUBROUTINE art_get_CO_des



SUBROUTINE art_get_CH4_des(jg,jb,jcs,jce ,&
            &             nproma,nlev, tracer )
!<
! SUBROUTINE art_get_CH4_des
! Calculates the destruction of CH4
! Based on:
! - KASIMA
! - Brasseur and Solomon?
! Part of Module: mo_art_chem_utils
! Author: Christian Stassen and  Michael Weimer, KIT
! Initial Release: 2018-10-18
! Modifications:
!>
  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &
    &  jg,                & !< patch id
    &  jb,jcs,jce           !< loop indices
  INTEGER, INTENT(in) ::  &
    &  nproma,nlev   !< dimensions
  TYPE(t_chem_meta_lt), INTENT(inout) ::   &
    &  tracer        !< lifetime tracer structure
  ! local variables
  TYPE(t_art_atmo),POINTER    :: &
    &  art_atmo                     !< Pointer to ART chem fields

  INTEGER  ::     &
    &  jc, jk

  INTEGER, PARAMETER ::   &
    &  ch4_case = 1  !< 1 means ECMWF parametrisation


  art_atmo => p_art_data(jg)%atmo


  SELECT CASE(ch4_case)
    CASE(1) 
      DO jk = 1,nlev
        DO jc = jcs,jce
          IF (art_atmo%pres(jc,jk,jb) <= 50._wp) THEN

            tracer%des_3d(jc,jk,jb) = 1._wp / (86400._wp*100._wp)

          ELSE IF (art_atmo%pres(jc,jk,jb) > 50._wp .AND. art_atmo%pres(jc,jk,jb) < 10000._wp) THEN

            tracer%des_3d(jc,jk,jb) = 1._wp / (86400._wp*(100._wp*(1._wp+( 19._wp*LOG(10._wp)    &
                & / (LOG(20._wp))**4._wp )*((LOG(art_atmo%pres(jc,jk,jb)/50._wp))**4._wp  &
                & / LOG(10000._wp/art_atmo%pres(jc,jk,jb)) ))))

          ELSE IF (art_atmo%pres(jc,jk,jb) >= 10000._wp) THEN

            tracer%des_3d(jc,jk,jb) = 0._wp

          END IF

        ENDDO
      ENDDO

    CASE DEFAULT
      DO jk = 1,nlev
        DO jc = jcs,jce
          CALL art_get_CH4_des_1d(tracer%des_3d(jc,jk,jb),art_atmo%pres(jc,jk,jb))
        ENDDO
      ENDDO
  END SELECT
END SUBROUTINE art_get_CH4_des

!!
!!-----------------------------------------------------------------------------
!!


SUBROUTINE art_get_CH3CN_des(jg,jb,jcs,jce ,&
            &             nproma,nlev, tracer )
!<
! SUBROUTINE art_get_CH3CN_des
! Calculates the destruction of CH3CN
! Part of Module: mo_art_chem_utils
! Author: Christian Stassen and  Michael Weimer, KIT
! Initial Release: 2018-10-18
! Modifications:
!>
  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &
    &  jg,                & !< patch id
    &  jb,jcs,jce           !< loop indices
  INTEGER, INTENT(in) ::  &
    &  nproma,nlev   !< dimensions
  TYPE(t_chem_meta_lt), INTENT(inout) ::   &
    &  tracer        !< lifetime tracer structure
  ! local variables
  TYPE(t_art_atmo),POINTER    :: &
    &  art_atmo                     !< Pointer to ART chem fields
  INTEGER   :: jc, jk  !< loop indices


  art_atmo => p_art_data(jg)%atmo


  DO jc = jcs,jce
    DO jk = 1,art_atmo%ktrpwmo(jc,jb)
       tracer%des_3d(jc,jk,jb) = 1._wp / 2207520000._wp
    ENDDO
    DO jk = art_atmo%ktrpwmo(jc,jb),nlev
       tracer%des_3d(jc,jk,jb)  = 1._wp / 37843200._wp
    ENDDO

    IF  (.NOT. art_atmo%llsm(jc,jb)) THEN
      DO jk = nlev-8,nlev
       tracer%des_3d(jc,jk,jb) = 1._wp / 1814400._wp
      ENDDO
    END IF
  ENDDO


END SUBROUTINE art_get_CH3CN_des



SUBROUTINE art_get_H2O_des(jg,jb,jcs,jce ,&
            &             nproma,nlev, tracer )
!<
! SUBROUTINE art_get_H2O_des
! Calculates the destruction of H2O
! Based on:
! - KASIMA
! - Brasseur and Solomon?
! Part of Module: mo_art_chem_utils
! Author: Christian Stassen and  Michael Weimer, KIT
! Initial Release: 2018-10-18
! Modifications:
!>
  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &
    &  jg,                & !< patch id
    &  jb,jcs,jce           !< loop indices
  INTEGER, INTENT(in) ::  &
    &  nproma,nlev   !< dimensions
  TYPE(t_chem_meta_lt), INTENT(inout) ::   &
    &  tracer        !< lifetime tracer structure
  ! local variables
  TYPE(t_art_atmo),POINTER    :: &
    &  art_atmo                     !< Pointer to ART chem fields
  INTEGER :: jc, jk  !< loop indices

  art_atmo => p_art_data(jg)%atmo


  DO jk = 1,nlev
    DO jc = jcs,jce

      IF (art_atmo%pres(jc,jk,jb) <= 0.1_wp) THEN
        tracer%des_3d(jc,jk,jb) = 1._wp/(86400._wp*3._wp)

      ELSE IF (art_atmo%pres(jc,jk,jb) > 0.1_wp .AND. art_atmo%pres(jc,jk,jb) < 20._wp) THEN
        tracer%des_3d(jc,jk,jb) =  1._wp / ( 86400._wp*(1._wp/(EXP(LOG(1._wp/3._wp+0.01_wp)  &
                &                                                     -0.5_wp*(LOG(100._wp)  &
                &                  +LOG(1._wp/3._wp+0.01_wp))*(1._wp+COS(pi                  &
                &                     *LOG(art_atmo%pres(jc,jk,jb)/20._wp) &
                &                   /LOG(0.005_wp))))-0.01_wp)) )

      ELSE IF (art_atmo%pres(jc,jk,jb) >= 20._wp) THEN
        tracer%des_3d(jc,jk,jb) = 0.0_wp
      END IF
    ENDDO
  ENDDO

END SUBROUTINE art_get_H2O_des

SUBROUTINE art_calc_col(jg,jb,jcs,jce ,&
            &             nproma,nlev, tracer )
!<
! SUBROUTINE art_calc_col
! This subroutine calculates the column of a linoz tracer in DU
! Part of Module: mo_art_chem_utils
! Author: Christian Stassen and  Michael Weimer, KIT
! Initial Release: 2018-10-18
! Modifications:
!>
  IMPLICIT NONE
  INTEGER, INTENT(in) ::  &
    &  jg,                & !< patch id
    &  jb,jcs,jce           !< loop indices
  INTEGER, INTENT(in) ::  &
    &  nproma,nlev   !< dimensions
  TYPE(t_chem_meta_linoz), INTENT(inout) ::   &
    &  tracer        !< linoz tracer structure
  ! local variables
  TYPE(t_art_atmo),POINTER    :: &
    &  art_atmo                     !< Pointer to ART chem fields
  INTEGER   :: jc, jk

  art_atmo => p_art_data(jg)%atmo

  ! first level
  DO jc=jcs, jce
    tracer%column(jc,1,jb) = 0._wp
  END DO

  ! other levels
  DO jk=2,nlev
    DO jc = jcs,jce
      ! p-tracer_now convert to DU
      ! [O3] = Mol/m^3   => [O3_col] = DU
      tracer%column(jc,jk,jb) = (tracer%tracer(jc,jk-1,jb) + tracer%tracer(jc,jk,jb))  &
                    & * 1.e6_wp / 2. / avo                                             &
                    & * ((art_atmo%z_mc(jc,jk-1,jb)) - (art_atmo%z_mc(jc,jk,jb)))      &
                    & + tracer%column(jc,jk-1,jb)

    ENDDO
  ENDDO

  DO jk=1,nlev
    DO jc = jcs,jce
      tracer%column(jc,jk,jb) = tracer%column(jc,jk,jb) / (0.44615e-3_wp)  ! Convert to DU
    END DO
  END DO
END SUBROUTINE art_calc_col


!!
!!------------------------------------------------------------------------
!!

SUBROUTINE art_calc_vmr2Nconc(jg,jb,jcs,jce ,&
            &             nproma,nlev,vmr2Nconc)
!<
! SUBROUTINE art_calc_vmr2Nconc
! Calculates the conversion factor from volume mixing ratio (mol/mol) to number
! concentration (cm-3)
! Based on: -
! Part of Module: mo_art_unit_conversion
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-20
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  IMPLICIT NONE
  INTEGER, INTENT(in) :: &
    &  jg, jb, jcs, jce, nproma, nlev
  REAL(wp), INTENT(inout) :: &
    &  vmr2Nconc(:,:,:)
  !local variables
  INTEGER :: &
    &  jc, jk
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  art_atmo => p_art_data(jg)%atmo

  DO jk = 1,nlev
    DO jc = jcs,jce
      vmr2Nconc(jc,jk,jb) = (avo * art_atmo%pres(jc,jk,jb))  &
                          &  / (argas * art_atmo%temp(jc,jk,jb)) * 1.e-6_wp
    END DO
  END DO

END SUBROUTINE art_calc_vmr2Nconc

SUBROUTINE art_convert_tracers_mmr_Nconc(jg,p_tracer_now,p_prog_list, vmr2Nconc)
!<
! SUBROUTINE art_convert_tracers_mmr_Nconc
! Converts all tracer mass mixing ratios (kg/kg) of chemical tracers to number
! concentration (cm-3)
! Based on: -
! Part of Module: mo_art_unit_conversion
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-20
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  IMPLICIT NONE
  INTEGER, INTENT(in) :: &
    &  jg                     !< patch id
  REAL(wp), INTENT(INOUT), TARGET ::  &
    &  p_tracer_now(:,:,:,:)  !< tracer concentrations (specific concentrations)
                              !< dim: (nproma,nlev,nblks_c,ntracer)
  TYPE(t_var_list_ptr), INTENT(INOUT)   ::  &
    &  p_prog_list         !< current prognostic state list
  REAL(wp), INTENT(in) :: &
    &  vmr2Nconc(:,:,:)
  ! local variables
  TYPE(t_var_metadata_dynamic), POINTER ::  &
    &  info_dyn         !< returns reference to tracer metadata of current element
  TYPE(t_var_metadata), POINTER ::  &
    &  info             !< returns reference to tracer  metadata of current element
  INTEGER, POINTER :: &
    &  jsp
  TYPE(t_art_chem), POINTER :: &
    &  art_chem
  INTEGER         :: &
    &  iv

  art_chem => p_art_data(jg)%chem


  DO iv = 1, p_prog_list%p%nvars

    info_dyn=>p_prog_list%p%vl(iv)%p%info_dyn
    info=>p_prog_list%p%vl(iv)%p%info

    jsp => info%ncontained

    SELECT TYPE(tracer => info_dyn%tracer)
      CLASS IS (t_chem_meta_param)
        CALL tracer%convert_mmr_Nconc(p_tracer_now(:,:,:,jsp), vmr2Nconc)
      CLASS IS (t_chem_meta_mecca)
        CALL tracer%convert_mmr_Nconc(p_tracer_now(:,:,:,jsp), vmr2Nconc)
      CLASS IS (t_tracer_meta)
        IF (jsp == iqv) THEN
          art_chem%water_tracers(:,:,:,IART_QV) = p_tracer_now(:,:,:,jsp) * amd / amw * vmr2Nconc
        ELSE IF (jsp == iqc) THEN
          art_chem%water_tracers(:,:,:,IART_QC) = p_tracer_now(:,:,:,jsp) * amd / amw * vmr2Nconc
        ELSE IF (jsp == iqi) THEN
          art_chem%water_tracers(:,:,:,IART_QI) = p_tracer_now(:,:,:,jsp) * amd / amw * vmr2Nconc
        END IF
    END SELECT

  END DO
END SUBROUTINE art_convert_tracers_mmr_Nconc

SUBROUTINE art_convert_tracers_Nconc_mmr(jg,p_tracer_now,p_prog_list, vmr2Nconc)
!<
! SUBROUTINE art_convert_tracers_Nconc_mmr
! Converts all tracer number concentrations (cm-3) of chemical tracers to mass
! mixing ratios (kg/kg)
! Based on: -
! Part of Module: mo_art_unit_conversion
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-20
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  IMPLICIT NONE
  INTEGER, INTENT(in) :: &
    &  jg                     !< patch id
  REAL(wp), INTENT(INOUT), TARGET ::  &
    &  p_tracer_now(:,:,:,:)  !< tracer concentrations (specific concentrations)
                              !< dim: (nproma,nlev,nblks_c,ntracer)
  TYPE(t_var_list_ptr), INTENT(INOUT)   ::  &
    &  p_prog_list         !< current prognostic state list
  REAL(wp), INTENT(in) :: &
    &  vmr2Nconc(:,:,:)
  ! local variables
  TYPE(t_var_metadata_dynamic), POINTER ::  &
    &  info_dyn         !< returns reference to tracer metadata of current element
  TYPE(t_var_metadata), POINTER ::  &
    &  info             !< returns reference to tracer  metadata of current element
  INTEGER, POINTER :: &
    &  jsp
  TYPE(t_art_chem), POINTER :: &
    &  art_chem
  INTEGER          :: &
    &  iv

  art_chem => p_art_data(jg)%chem


  DO iv = 1, p_prog_list%p%nvars

    info_dyn=> p_prog_list%p%vl(iv)%p%info_dyn
    info=>p_prog_list%p%vl(iv)%p%info

    jsp => info%ncontained

    SELECT TYPE(tracer => info_dyn%tracer)
      CLASS IS (t_chem_meta_param)
        CALL tracer%convert_Nconc_mmr(p_tracer_now(:,:,:,jsp), vmr2Nconc)
      CLASS IS (t_chem_meta_mecca)
        CALL tracer%convert_Nconc_mmr(p_tracer_now(:,:,:,jsp), vmr2Nconc)
      CLASS IS (t_tracer_meta)
        IF (jsp == iqv) THEN
          art_chem%water_tracers(:,:,:,IART_QV) = art_chem%water_tracers(:,:,:,IART_QV) / amd * amw / vmr2Nconc
        ELSE IF (jsp == iqc) THEN
          art_chem%water_tracers(:,:,:,IART_QC) = art_chem%water_tracers(:,:,:,IART_QC) / amd * amw / vmr2Nconc
        ELSE IF (jsp == iqi) THEN
          art_chem%water_tracers(:,:,:,IART_QI) = art_chem%water_tracers(:,:,:,IART_QI) / amd * amw / vmr2Nconc
        END IF
    END SELECT

  END DO
END SUBROUTINE art_convert_tracers_Nconc_mmr


SUBROUTINE art_convert_tracers_vmr_mmr(p_tracer_now,p_prog_list)
!<
! SUBROUTINE art_convert_tracers_vmr_mmr
! converts the tracer volume mixing ratio (mol/mol, gitven at initialisation)
! to mass mixing ratio (kg/kg)
! Based on: -
! Part of Module: mo_art_unit_conversion
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-20
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  REAL(wp), INTENT(INOUT), TARGET ::  &
    &  p_tracer_now(:,:,:,:)  !< tracer mixing ratios (specific concentrations)
                              !< dim: (nproma,nlev,nblks_c,ntracer)
  TYPE(t_var_list_ptr), INTENT(IN)   ::  &
    &  p_prog_list         !< current prognostic state list
  ! local variables
  TYPE(t_var_metadata_dynamic), POINTER ::  &
    &  info_dyn         !< returns reference to tracer metadata of current element
  TYPE(t_var_metadata), POINTER ::  &
    &  info             !< returns reference to tracer  metadata of current element
  INTEGER, POINTER :: &
    &  jsp
  INTEGER          :: &
    &  iv


  DO iv = 1, p_prog_list%p%nvars

    info_dyn=> p_prog_list%p%vl(iv)%p%info_dyn
    info=>p_prog_list%p%vl(iv)%p%info

    jsp => info%ncontained

    SELECT TYPE(tracer => info_dyn%tracer)
      CLASS IS (t_chem_meta_param)
        CALL tracer%convert_vmr_mmr(p_tracer_now(:,:,:,jsp))
      CLASS IS (t_chem_meta_mecca)
        CALL tracer%convert_vmr_mmr(p_tracer_now(:,:,:,jsp))
    END SELECT

  END DO
END SUBROUTINE art_convert_tracers_vmr_mmr

END MODULE mo_art_chem_utils
