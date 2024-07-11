! @brief Interface between atmosphere physics and the ocean, through a coupler
!
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_atmo_ocean_coupling

  USE mo_kind,            ONLY: wp
  USE mo_model_domain,    ONLY: t_patch
  USE mo_ext_data_types,  ONLY: t_external_data
#ifndef __NO_AES__
  USE mo_aes_phy_memory,  ONLY: prm_field
#endif
  USE mo_parallel_config, ONLY: nproma
  USE mo_impl_constants,  ONLY: inwp, iaes, SUCCESS
  USE mo_mpi,             ONLY: p_pe_work, p_comm_work, p_sum
  USE mo_run_config,      ONLY: iforcing
  USE mo_util_dbg_prnt,   ONLY: dbg_print
  USE mo_exception,       ONLY: finish
  USE mo_coupling_utils,  ONLY: cpl_def_cell_field_mask, cpl_def_field

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_atmo_ocean_coupling' ! Output of module for debug

  PUBLIC :: construct_atmo_ocean_coupling
  PUBLIC :: field_id_umfl, field_id_vmfl, field_id_freshflx, &
            field_id_heatflx, field_id_seaice_atm, field_id_sst, &
            field_id_oce_u, field_id_oce_v, field_id_seaice_oce, &
            field_id_sp10m, field_id_co2_vmr, field_id_co2_flx, &
            field_id_pres_msl
  PUBLIC :: mask_checksum

  INTEGER :: field_id_umfl
  INTEGER :: field_id_vmfl
  INTEGER :: field_id_freshflx
  INTEGER :: field_id_heatflx
  INTEGER :: field_id_seaice_atm
  INTEGER :: field_id_sst
  INTEGER :: field_id_oce_u
  INTEGER :: field_id_oce_v
  INTEGER :: field_id_seaice_oce
  INTEGER :: field_id_sp10m
  INTEGER :: field_id_co2_vmr
  INTEGER :: field_id_co2_flx
  INTEGER :: field_id_pres_msl

  INTEGER, SAVE :: mask_checksum = -1

CONTAINS

  !>
  !! Registers fields required for the coupling between atmo and
  !! ocean
  !!
  !! This subroutine is called from construct_atmo_coupling.
  !!
  SUBROUTINE construct_atmo_ocean_coupling( &
    p_patch, ext_data, comp_id, grid_id, cell_point_id, timestepstring)

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch(:)
    TYPE(t_external_data), INTENT(IN) :: ext_data(:)
    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: grid_id
    INTEGER, INTENT(IN) :: cell_point_id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring

    TYPE(t_patch), POINTER :: patch_horz

    INTEGER :: cell_mask_id

    INTEGER :: jg, jb, jc, error

    LOGICAL,  ALLOCATABLE :: is_valid(:)

    REAL(wp), ALLOCATABLE :: lsmnolake(:,:)

    REAL(wp), PARAMETER :: eps = 1.E-10_wp

    CHARACTER(LEN=*), PARAMETER   :: &
      routine = str_module // ':construct_atmo_ocean_coupling'

    jg = 1
    patch_horz => p_patch(jg)

    ! The integer land-sea mask:
    !          -2: inner ocean
    !          -1: boundary ocean
    !           1: boundary land
    !           2: inner land
    !
    ! The (fractional) mask which is used in the AES physics is prm_field(1)%lsmask(:,:).
    !
    ! The logical mask for the coupler must be generated from the fractional mask by setting
    !   only those gridpoints to land that have no ocean part at all (lsf<1 is ocean).
    ! The logical mask is then set to .FALSE. for land points to exclude them from mapping by yac.
    ! These points are not touched by yac.
    !

    ALLOCATE(is_valid(nproma*patch_horz%nblks_c), STAT = error)
    IF(error /= SUCCESS) CALL finish(str_module, "memory allocation failure for is_valid")
    ALLOCATE(lsmnolake(nproma,patch_horz%nblks_c), STAT = error)
    IF(error /= SUCCESS) CALL finish(str_module, "memory allocation failure for lsmnolake")

    SELECT CASE( iforcing ) !{{{

      CASE ( inwp )

        !ICON_OMP_PARALLEL PRIVATE(jb,jc)
          !ICON_OMP_WORKSHARE
          is_valid(:) = .FALSE.
          !ICON_OMP_END_WORKSHARE

          !ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
          DO jb = 1, patch_horz%nblks_c
            DO jc = 1, ext_data(jg)%atm%list_sea%ncount(jb)
              is_valid((jb-1)*nproma + ext_data(jg)%atm%list_sea%idx(jc,jb)) = .TRUE.
            END DO
          END DO
          !ICON_OMP_END_DO
        !ICON_OMP_END_PARALLEL

        CALL dbg_print('AtmFrame: fr_land',ext_data(jg)%atm%fr_land,str_module,3,in_subset=patch_horz%cells%owned)
        CALL dbg_print('AtmFrame: fr_lake',ext_data(jg)%atm%fr_lake,str_module,3,in_subset=patch_horz%cells%owned)

      CASE ( iaes )
#ifdef __NO_AES__
        CALL finish (str_module // ':construct_atmo_ocean_coupling', &
            & 'coupled model needs aes; remove --disable-aes and reconfigure')
#else
        !ICON_OMP_PARALLEL_DO PRIVATE(jb,jc) ICON_OMP_RUNTIME_SCHEDULE
        DO jb = 1, patch_horz%nblks_c
          DO jc = 1, nproma
            !  slo: caution - lsmask includes alake, must be added to refetch pure lsm:
            lsmnolake(jc, jb) = prm_field(jg)%lsmask(jc,jb) + prm_field(jg)%alake(jc,jb)
          ENDDO
        ENDDO
        !ICON_OMP_END_PARALLEL_DO

        mask_checksum = 0
        !ICON_OMP_PARALLEL_DO PRIVATE(jb,jc) REDUCTION(+:mask_checksum) ICON_OMP_RUNTIME_SCHEDULE
        DO jb = 1, patch_horz%nblks_c
          DO jc = 1, nproma
            mask_checksum = mask_checksum + ABS( lsmnolake(jc,jb))
          ENDDO
        ENDDO
        !ICON_OMP_END_PARALLEL_DO

        mask_checksum = p_sum(mask_checksum, comm=p_comm_work)

        !
        ! Define cell_mask_ids(1): all ocean and coastal points are valid
        !   This is the standard for the coupling of atmospheric fields listed below
        !
        IF ( mask_checksum > 0 ) THEN
          !ICON_OMP_PARALLEL_DO PRIVATE(jb, jc) ICON_OMP_RUNTIME_SCHEDULE
          DO jb = 1, patch_horz%nblks_c
            DO jc = 1, nproma

              IF ( lsmnolake(jc, jb) .LT. (1.0_wp - eps) ) THEN
                ! ocean point (fraction of ocean is >0., lsmnolake .lt. 1.) is valid
                is_valid((jb-1)*nproma+jc) = .TRUE.
              ELSE
                ! land point (fraction of land is one, no sea water, lsmnolake=1.) is undef
                is_valid((jb-1)*nproma+jc) = .FALSE.
              ENDIF

            ENDDO
          ENDDO
          !ICON_OMP_END_PARALLEL_DO
        ELSE
          !ICON_OMP_PARALLEL_DO PRIVATE(jb, jc) ICON_OMP_RUNTIME_SCHEDULE
          DO jc = 1, patch_horz%nblks_c * nproma
            is_valid(jc) = .TRUE.
          ENDDO
          !ICON_OMP_END_PARALLEL_DO
        ENDIF
#endif
        CASE DEFAULT

          CALL finish ('Please mask handling for new forcing in ' &
            & //'src/coupling/mo_atmo_ocean_coupling: ' &
            & //'construct_atmo_ocean_coupling. Thank you!')

    END SELECT !}}}

    CALL cpl_def_cell_field_mask( &
      routine, grid_id, is_valid, cell_mask_id )

    DEALLOCATE (is_valid, STAT = error)
    IF(error /= SUCCESS) CALL finish(str_module, "Deallocation failed for is_valid")

    DEALLOCATE (lsmnolake, STAT = error)
    IF(error /= SUCCESS) CALL finish(str_module, "Deallocation failed for lsmnolake")

    CALL cpl_def_field( &
      comp_id, cell_point_id, cell_mask_id, timestepstring, &
      "surface_downward_eastward_stress", 2, field_id_umfl)

    CALL cpl_def_field( &
      comp_id, cell_point_id, cell_mask_id, timestepstring, &
      "surface_downward_northward_stress", 2, field_id_vmfl)

    CALL cpl_def_field( &
      comp_id, cell_point_id, cell_mask_id, timestepstring, &
      "surface_fresh_water_flux", 3, field_id_freshflx)

    CALL cpl_def_field( &
      comp_id, cell_point_id, cell_mask_id, timestepstring, &
      "total_heat_flux", 4, field_id_heatflx)

    CALL cpl_def_field( &
      comp_id, cell_point_id, cell_mask_id, timestepstring, &
      "atmosphere_sea_ice_bundle", 2, field_id_seaice_atm)

    CALL cpl_def_field( &
      comp_id, cell_point_id, cell_mask_id, timestepstring, &
      "sea_surface_temperature", 1, field_id_sst)

    CALL cpl_def_field( &
      comp_id, cell_point_id, cell_mask_id, timestepstring, &
      "eastward_sea_water_velocity", 1, field_id_oce_u)

    CALL cpl_def_field( &
      comp_id, cell_point_id, cell_mask_id, timestepstring, &
      "northward_sea_water_velocity", 1, field_id_oce_v)

    CALL cpl_def_field( &
      comp_id, cell_point_id, cell_mask_id, timestepstring, &
      "ocean_sea_ice_bundle", 3, field_id_seaice_oce)

    CALL cpl_def_field( &
      comp_id, cell_point_id, cell_mask_id, timestepstring, &
      "10m_wind_speed", 1, field_id_sp10m)

    CALL cpl_def_field( &
      comp_id, cell_point_id, cell_mask_id, timestepstring, &
      "co2_mixing_ratio", 1, field_id_co2_vmr)

    CALL cpl_def_field( &
      comp_id, cell_point_id, cell_mask_id, timestepstring, &
      "co2_flux", 1, field_id_co2_flx)

    CALL cpl_def_field( &
      comp_id, cell_point_id, cell_mask_id, timestepstring, &
      "sea_level_pressure", 1, field_id_pres_msl)

  END SUBROUTINE construct_atmo_ocean_coupling

END MODULE mo_atmo_ocean_coupling
