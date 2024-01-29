! Initialization/reading reading of external datasets
!
! This module contains read and initialization routines for the external data state.
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

MODULE mo_wave_ext_data_init

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, SUCCESS
  USE mo_io_units,            ONLY: filename_max
  USE mo_io_config,           ONLY: default_read_method
  USE mo_exception,           ONLY: message, finish
  USE mo_model_domain,        ONLY: t_patch
  USE mo_grid_config,         ONLY: n_dom, nroot
  USE mo_var_list,            ONLY: t_var_list_ptr
  USE mo_extpar_config,       ONLY: extpar_filename, generate_filename
  USE mo_read_interface,      ONLY: openInputFile, closeFile, t_stream_id, on_cells, read_2D
  USE mo_master_config,       ONLY: getModelBaseDir
  USE mo_intp,                ONLY: cells2edges_scalar
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_fortran_tools,       ONLY: copy, init
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_process_topo,        ONLY: compute_smooth_topo

  USE mo_wave_ext_data_types, ONLY: t_external_wave
  USE mo_wave_ext_data_state, ONLY: construct_wave_ext_data_state
  USE mo_wave_config,         ONLY: wave_config




  IMPLICIT NONE
  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_data_init'

  PUBLIC :: init_wave_ext_data

CONTAINS

  !>
  !! construct and read ext data state
  !!
  SUBROUTINE init_wave_ext_data (p_patch, p_int_state, wave_ext_data, wave_ext_data_list)

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_wave_ext_data'

    TYPE(t_patch),                      INTENT(INOUT) :: p_patch(:)
    TYPE(t_int_state),                  INTENT(IN)    :: p_int_state(:)
    TYPE(t_external_wave), ALLOCATABLE, INTENT(INOUT) :: wave_ext_data(:)
    TYPE(t_var_list_ptr),  ALLOCATABLE, INTENT(INOUT) :: wave_ext_data_list(:)

    INTEGER :: jg
    INTEGER :: ist

    CHARACTER(len=20) :: depth_str

    REAL(wp), ALLOCATABLE :: topo_smt_c(:,:)

    CALL construct_wave_ext_data_state(p_patch, wave_ext_data, wave_ext_data_list)

    ! read depth from file
    CALL read_ext_data_wave(p_patch, wave_ext_data)

    DO jg = 1, n_dom
      IF (wave_config(jg)%depth > 0.0_wp) THEN
        ! set depth to constant
        write(depth_str,'(f9.2)') wave_config(jg)%depth
        CALL message('###  Run with constant depth of',TRIM(depth_str)//', m')
        wave_ext_data(jg)%bathymetry_c(:,:) = wave_config(jg)%depth
      ELSE

        IF (wave_config(jg)%niter_smooth > 0) THEN
          ! smooth bathymetry
          
          ALLOCATE(topo_smt_c(&
               SIZE(wave_ext_data(jg)%bathymetry_c,1), &
               SIZE(wave_ext_data(jg)%bathymetry_c,2)),&
               stat=ist)
          IF (ist/=SUCCESS) CALL finish(routine, &
               'allocation of topo_smt_c array failed')

          CALL compute_smooth_topo(p_patch(jg), p_int_state(jg), &
               wave_ext_data(jg)%bathymetry_c, &
               wave_config(jg)%niter_smooth, &
               topo_smt_c)

!$OMP PARALLEL
          CALL copy(src=topo_smt_c,dest=wave_ext_data(jg)%bathymetry_c)
!$OMP END PARALLEL

          DEALLOCATE(topo_smt_c,stat=ist)
          IF (ist/=SUCCESS) CALL finish(routine, &
               'deallocation of topo_smt_c array failed')
          
        END IF
      END IF
    END DO

    ! cell2edge interpolation of bathymetry
    DO jg = 1, n_dom
      CALL cells2edges_bathymetry(p_patch      = p_patch(jg),                    &
        &                         p_int_state  = p_int_state(jg),                &
        &                         bathymetry_c = wave_ext_data(jg)%bathymetry_c, &
        &                         bathymetry_e = wave_ext_data(jg)%bathymetry_e)
    END DO

    CALL message(TRIM(routine),'finished.')

  END SUBROUTINE init_wave_ext_data


  !>
  !! cell2edge interpolation for bathymetry
  !!
  SUBROUTINE cells2edges_bathymetry(p_patch, p_int_state, bathymetry_c, bathymetry_e)
    TYPE(t_patch),      INTENT(IN)    :: p_patch
    TYPE(t_int_state),  INTENT(IN)    :: p_int_state
    REAL(wp),           INTENT(IN)    :: bathymetry_c(:,:)
    REAL(wp),           INTENT(INOUT) :: bathymetry_e(:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':cells2edges_bathymetry'

    REAL(wp):: bath_c_3d(SIZE(bathymetry_c,1),1,SIZE(bathymetry_c,2))
    REAL(wp):: bath_e_3d(SIZE(bathymetry_e,1),1,SIZE(bathymetry_e,2))

    INTEGER :: jb, je
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx,i_endidx

    rl_start   = 1
    rl_end     = min_rledge
    i_startblk = p_patch%edges%start_block(rl_start)
    i_endblk   = p_patch%edges%end_block(rl_end)


!$OMP PARALLEL
    CALL copy(src=bathymetry_c, dest=bath_c_3d(:,1,:))
    CALL init(bath_e_3d(:,:,:))
!$OMP END PARALLEL

    CALL cells2edges_scalar(bath_c_3d, p_patch, p_int_state%c_lin_e, bath_e_3d)

!$OMP PARALLEL
    CALL copy(src=bath_e_3d(:,1,:), dest=bathymetry_e)
!$OMP BARRIER

!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb=i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      DO je = i_startidx, i_endidx
          bathymetry_e(je,jb) = MAX(bathymetry_e(je,jb),0.2_wp)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE cells2edges_bathymetry


  SUBROUTINE read_ext_data_wave(p_patch, wave_ext_data)

    TYPE(t_patch),         INTENT(IN)    :: p_patch(:)
    TYPE(t_external_wave), INTENT(INOUT) :: wave_ext_data(:)

    INTEGER :: jg
    INTEGER :: jb, jc
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx,i_endidx

    TYPE(t_stream_id) :: stream_id

    CHARACTER(filename_max) :: extpar_file

    DO jg = 1, n_dom
      extpar_file = generate_filename(extpar_filename, getModelBaseDir(), &
        &                             TRIM(p_patch(jg)%grid_filename),    &
        &                              nroot,                             &
        &                             p_patch(jg)%level, p_patch(jg)%id)

      CALL openInputFile(stream_id, extpar_file, p_patch(jg), default_read_method)

      CALL read_2D(stream_id, on_cells, 'z', wave_ext_data(jg)%bathymetry_c)

      CALL closeFile(stream_id)


      rl_start   = 1
      rl_end     = min_rlcell
      i_startblk = p_patch(jg)%cells%start_block(rl_start)
      i_endblk   = p_patch(jg)%cells%end_block(rl_end)

      ! set minimal depth of 0.2 m
      !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb=i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        DO jc = i_startidx, i_endidx
           wave_ext_data(jg)%bathymetry_c(jc,jb) = MAX(wave_ext_data(jg)%bathymetry_c(jc,jb),0.2_wp)
        END DO

      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    END DO

  END SUBROUTINE read_ext_data_wave

END MODULE mo_wave_ext_data_init
