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

! Interface between atmosphere physics and the ocean surface waves, through a coupler

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_atmo_wave_coupling

  USE mo_kind,               ONLY: wp
  USE mo_model_domain,       ONLY: t_patch
  USE mo_fortran_tools,      ONLY: assert_acc_host_only
  USE mo_coupling_utils,     ONLY: cpl_def_field, cpl_put_field, cpl_get_field
  USE mo_idx_list,           ONLY: t_idx_list_blocked
  USE mo_lnd_nwp_config,     ONLY: isub_water
  USE mo_physical_constants, ONLY: grav
  USE mo_impl_constants,     ONLY: min_rlcell
  USE mo_loopindices,        ONLY: get_indices_c


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_atmo_wave_coupling, couple_atmo_to_wave

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_atmo_wave_coupling'

  INTEGER :: field_id_u10m
  INTEGER :: field_id_v10m
  INTEGER :: field_id_fr_seaice
  INTEGER :: field_id_z0

CONTAINS

  !>
  !! Registers fields required for the coupling between atmosphere and wave
  !!
  !! This subroutine is called from construct_atmo_coupling.
  !!
  SUBROUTINE construct_atmo_wave_coupling( &
    comp_id, cell_point_id, timestepstring)

    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: cell_point_id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring

    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "zonal_wind_in_10m", 1, field_id_u10m)
    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "meridional_wind_in_10m", 1, field_id_v10m)
    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "fraction_of_ocean_covered_by_sea_ice", 1, field_id_fr_seaice)
    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "roughness_length", 1, field_id_z0)

  END SUBROUTINE construct_atmo_wave_coupling

  !>
  !! Exchange fields between atmosphere and wave model
  !!
  !! Send fields to the wave model:
  !!   "meridional_wind_in_10m"
  !!   "zonal_wind_in_10m"
  !!   "fraction_of_ocean_covered_by_sea_ice"
  !!
  !! Receive fields from the wave model:
  !!   "roughness_length"
  !!
  !! This subroutine is called from nwp_nh_interface.
  !!
  SUBROUTINE couple_atmo_to_wave(p_patch, list_sea, u10m, v10m, fr_seaice, frac_t, &
    &                            z0_waves, gz0_t, gz0, lacc)

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':couple_atmo_to_wave'

    TYPE(t_patch),                INTENT(IN)    :: p_patch
    TYPE(t_idx_list_blocked),     INTENT(IN)    :: list_sea
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)    :: u10m(:,:)      !< zonal wind speed in 10m [m/s]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)    :: v10m(:,:)      !< meridional wind speed in 10m [m/s]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN)    :: fr_seaice(:,:) !< fraction_of_ocean_covered_by_sea_ice [1]
    REAL(wp), CONTIGUOUS,         INTENT(IN)    :: frac_t(:,:,:)  !< tile-specific area fraction [1]
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT) :: z0_waves(:,:)  !< surface roughness length [m]
    REAL(wp), CONTIGUOUS,         INTENT(INOUT) :: gz0_t(:,:,:)   !< tile-based roughness length times gravity [m2 s-2]
    REAL(wp), CONTIGUOUS,         INTENT(INOUT) :: gz0(:,:)       !< aggregated roughness length times gravity [m2 s-2]
    LOGICAL,  OPTIONAL,           INTENT(IN)    :: lacc           ! If true, use openacc

    LOGICAL :: write_coupler_restart, received_data

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,ic,jc
    INTEGER :: isubs

    CALL assert_acc_host_only('couple_atmo_to_wave', lacc)

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Send fields from atmosphere to wave
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !
    !------------------------------------------------
    !  Send 10m zonal wind u10m
    !  'zonal_wind_in_10m'
    !------------------------------------------------
    !
    CALL cpl_put_field( &
      routine, field_id_u10m, 'U10', p_patch%n_patch_cells, u10m)

    ! ----------------------------------------------
    !  Send 10m meridional wind v10m
    !  'meridional_wind_in_10m'
    ! ----------------------------------------------
    !
    CALL cpl_put_field( &
      routine, field_id_v10m, 'V10', p_patch%n_patch_cells, v10m)

    ! ------------------------------------------------------------------
    !  Send fraction of sea ice
    !  'fraction_of_ocean_covered_by_sea_ice'
    ! ------------------------------------------------------------------
    !
    CALL cpl_put_field( &
      routine, field_id_fr_seaice, 'fr_seaice', p_patch%n_patch_cells, &
      fr_seaice, write_restart=write_coupler_restart)

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Receive fields from wave to atmosphere
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !
    ! --------------------------------------------
    !  Receive roughness length z0 from the wave model
    !  'roughness_length'
    ! --------------------------------------------
    !
    CALL cpl_get_field( &
      routine, field_id_z0, 'z0', p_patch%n_patch_cells, z0_waves, &
      first_get=.TRUE., received_data=received_data)

    IF (received_data) THEN

      ! update gz0_t on water tiles only if

      i_rlstart  = 1
      i_rlend    = min_rlcell
      i_startblk = p_patch%cells%start_block(i_rlstart)
      i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,ic,jc,i_startidx,i_endidx,isubs) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
          &                 i_startidx, i_endidx, i_rlstart, i_rlend)
        DO ic = 1, list_sea%ncount(jb)
          jc = list_sea%idx(ic,jb)
          z0_waves(jc,jb) = MAX(z0_waves(jc,jb),1.e-6_wp)
          gz0_t(jc,jb,isub_water) = grav * z0_waves(jc,jb)
        END DO

        ! aggregate gz0_t
        DO jc = i_startidx, i_endidx
          gz0(jc,jb) = 0._wp
        ENDDO
        !
        DO isubs = 1, SIZE(frac_t,3)
          DO jc = i_startidx, i_endidx
            gz0(jc,jb)= gz0(jc,jb) + gz0_t(jc,jb,isubs) * frac_t(jc,jb,isubs)
          ENDDO
        ENDDO  !isubs
      ENDDO  !jb
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

    END IF

  END SUBROUTINE couple_atmo_to_wave

END MODULE mo_atmo_wave_coupling
