!> Contains the ECHAM driver for jsbach standalone
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
!>
PROGRAM jsb4_driver
#ifndef __NO_JSBACH__
#ifndef __ICON__

  USE mo_jsb_interface,      ONLY: jsbach_interface, jsbach_start_timestep, jsbach_finish_timestep
  USE mo_jsb_base,           ONLY: jsbach_setup_models, jsbach_setup_tiles
  USE mo_jsb_model_init,     ONLY: jsbach_setup_grid, jsbach_init, jsbach_init_after_restart
  USE mo_jsb_model_final,    ONLY: jsbach_finalize
  USE mo_jsb_version,        ONLY: jsbach_init_version, jsbach_label_run
  USE mo_jsb_control,        ONLY: jsbach_runs_standalone, jsbach_is_restarted, timer_on, timer_forcing
  USE mo_jsb_time,           ONLY: is_newday
  USE mo_jsb_class,          ONLY: get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_mvstream,           ONLY: init_mvstream, mvstream_update_cache, mvstream_accumulate
  USE mo_sub_nml,            ONLY: set_stream_element_nml, set_stream_nml
  USE mo_jsb_time,           ONLY: t_datetime, deallocateDatetime

  dsl4jsb_Use_processes SEB_, TURB_, HYDRO_
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(TURB_)
  dsl4jsb_Use_memory(HYDRO_)

  USE mo_jsb4_forcing,       ONLY: get_interface_variables_from_external_forcing, init_forcing, finalize_external_forcing

  ! Currently from adapter
  USE mo_jsb_grid_iface,     ONLY: get_lat, get_lon
  USE mo_jsb_namelist_iface, ONLY: POSITIONED, open_nml, position_nml
  USE mo_jsb_parallel,       ONLY: my_process_is_stdio, p_io, p_bcast, my_process_is_mpi_parallel, mpi_comm, p_global_comm
  USE mo_jsb_test,           ONLY: read_interface_vars, read_interface_variables

  ! Seems to be available in both ICON and ECHAM
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_kind,               ONLY: wp
  USE mo_jsb_math_constants, ONLY: deg2rad
  USE mo_util,               ONLY : int2string

!>>> TODO? echam:
  ! parallelisation
  USE mo_mpi,              ONLY: p_start, p_stop, p_init_communicators
  USE mo_machine,          ONLY: machine_setup

  ! io
  USE mo_io,               ONLY: write_streams, IO_read_streams
  USE mo_filename,         ONLY: out_filetype, rerun_filetype, NETCDF, NETCDF2,  NETCDF4
  USE mo_netcdf,           ONLY: nf_check, nf_close, add_dim, SURFACE
  USE mo_memory_base,      ONLY: default_output
  USE mo_output,           ONLY: open_output_streams, close_output_streams, out_streams, init_output
  USE mo_decomposition,    ONLY: local_decomposition

  ! time control and orbit calc
  USE mo_time_control,     ONLY: time_set, time_reset, lstop, l_trigfiles, l_putrerun, delta_time, &
    &                            time_step_len, current_date, next_date, write_date, &
    &                            echam_time, manager_state
  USE mo_orbit_solar,      ONLY: compute_cos_zenith_angle

  USE mo_timer, ONLY: init_timer, cleanup_timer, timer_start, timer_stop, timer_loop, timer_output, timer_restart

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER :: progname = 'jsb4_driver'

  INTEGER,  EXTERNAL:: util_cputime
  REAL(wp), EXTERNAL:: util_walltime
  REAL(wp):: zwtime, zutime, zstime, zrtime

  INTEGER:: nml_handler, read_status, f_unit, status, iblk, nproma, npromz, nblks, ncid, model_id, ics, ice
  TYPE(t_datetime)          :: current_datetime
  TYPE(t_datetime), POINTER :: current_datetime_ptr => NULL(), next_datetime_ptr => NULL()

  !! Variables of namelist jsb_parallel_nml
  INTEGER:: nproca  ! number of processors for jsbach
  INTEGER:: nprocb  ! has to be one
  INTEGER:: nprocio ! number of processors for I/O server
  INTEGER:: npedim  ! Working dimension for blocks in each domain.
                     ! Default (-1): each domain is processed in one call
  INTEGER:: buffer(4) ! buffer to broadcast the 4 above variables

  REAL(wp), POINTER     :: lon(:,:), lat(:,:)
  REAL(wp), ALLOCATABLE :: coslon(:,:), coslat(:,:), sinlon(:,:), sinlat(:,:)

  TYPE(t_jsb_model), POINTER :: model
  CLASS(t_jsb_tile_abstract),  POINTER :: tile

  dsl4jsb_Def_memory(SEB_)
  dsl4jsb_Def_memory(TURB_)
  dsl4jsb_Def_memory(HYDRO_)

  !! JSB4 interface variables - input
  REAL(wp), ALLOCATABLE :: t_air_K(:,:)
  REAL(wp), ALLOCATABLE :: q_air(:,:)
  REAL(wp), ALLOCATABLE :: rain(:,:)
  REAL(wp), ALLOCATABLE :: snow(:,:)
  REAL(wp), ALLOCATABLE :: wind_air(:,:)
  REAL(wp), ALLOCATABLE :: wind_10m(:,:)
  REAL(wp), ALLOCATABLE :: lw_srf_down(:,:)
  REAL(wp), ALLOCATABLE :: swvis_srf_down(:,:)
  REAL(wp), ALLOCATABLE :: swnir_srf_down(:,:)
  REAL(wp), ALLOCATABLE :: swpar_srf_down(:,:)
  REAL(wp), ALLOCATABLE :: press_srf(:,:)
  REAL(wp), ALLOCATABLE :: drag_srf(:,:)
  REAL(wp), ALLOCATABLE :: t_acoef(:,:)
  REAL(wp), ALLOCATABLE :: t_bcoef(:,:)
  REAL(wp), ALLOCATABLE :: q_acoef(:,:)
  REAL(wp), ALLOCATABLE :: q_bcoef(:,:)
  REAL(wp), ALLOCATABLE :: pch(:,:)
  REAL(wp), ALLOCATABLE :: cos_zenith_angle(:,:)

  !! TODO: so far not used in jsb4 but already read from forcing file/ calculated in mo_jsb4_forcing
  REAL(wp), ALLOCATABLE :: CO2_concentration(:,:)
  REAL(wp), ALLOCATABLE :: fract_par_diffuse(:,:)

  !! JSB4 interface variables - output, required for next timesteps forcing (if not reading interface vars)
  REAL(wp), ALLOCATABLE :: t_srf_proc(:,:)
  REAL(wp), ALLOCATABLE :: fact_q_air_proc(:,:)
  REAL(wp), ALLOCATABLE :: fact_qsat_srf_proc(:,:)
  REAL(wp), ALLOCATABLE :: rough_h_srf_proc(:,:)
  REAL(wp), ALLOCATABLE :: rough_m_srf_proc(:,:)
  REAL(wp), ALLOCATABLE :: evapo_act2pot_proc(:,:)
!  REAL(wp), ALLOCATABLE :: evapo_act_sum_proc(:,:)
!  REAL(wp), ALLOCATABLE :: evapo_pot_sum_proc(:,:)


  !! JSB4 interface variables - output
  REAL(wp), ALLOCATABLE :: t_eff_srf(:)
  REAL(wp), ALLOCATABLE :: qsat_srf(:)
  REAL(wp), ALLOCATABLE :: s_srf(:)
  REAL(wp), ALLOCATABLE :: evapotrans(:)
  REAL(wp), ALLOCATABLE :: evapopot(:)
  REAL(wp), ALLOCATABLE :: latent_hflx(:)
  REAL(wp), ALLOCATABLE :: sensible_hflx(:)
  REAL(wp), ALLOCATABLE :: grnd_hflx(:)
  REAL(wp), ALLOCATABLE :: grnd_hcap(:)
  REAL(wp), ALLOCATABLE :: q_snocpymlt(:)
  REAL(wp), ALLOCATABLE :: alb_vis_dir(:)
  REAL(wp), ALLOCATABLE :: alb_nir_dir(:)
  REAL(wp), ALLOCATABLE :: alb_vis_dif(:)
  REAL(wp), ALLOCATABLE :: alb_nir_dif(:)
  REAL(wp), ALLOCATABLE :: CO2_flux(:)

    NAMELIST /jsb_parallel_nml/     &
      nproca, npedim,             &
      nprocb, nprocio

  !-------------------------------------------------------------------------------------------------

  model_id = 1
  default_output = .TRUE.
  out_filetype = NETCDF4
  rerun_filetype = NETCDF4

  ! Start MPI
  CALL p_start

  ! Initialize wallclock timer
!$OMP PARALLEL
!$OMP MASTER
  zwtime = util_walltime(0)
!$OMP END MASTER
!$OMP END PARALLEL

  ! read processor decomposition
  IF (my_process_is_stdio()) THEN
    ! define default values
    nproca = 1
    nprocb = 1
    nprocio = 0
    npedim = -1

    ! read the namelist
    nml_handler = open_nml ('namelist.jsbach')
    f_unit = position_nml ('jsb_parallel_nml', nml_handler, status=read_status)
    IF (read_status == POSITIONED) READ(f_unit, jsb_parallel_nml)
  ENDIF

  IF (my_process_is_mpi_parallel()) THEN
    buffer(:) = (/ nproca, nprocb, nprocio, npedim /)
    CALL p_bcast (buffer, p_io, comm=p_global_comm)
    nproca  = buffer(1)
    nprocb  = buffer(2)
    nprocio = buffer(3)
    npedim  = buffer(4)
  END IF

  CALL p_init_communicators(nproca, nprocb, nprocio)
  CALL machine_setup

  CALL add_dim ("surface",1,                      levtyp =    1,     &
                                                  single = .TRUE.,   &
                                                  indx   = SURFACE)

  ! - jsbach_init
  CALL jsbach_setup_models('namelist.jsbach')
  IF (timer_on()) CALL init_timer()
  CALL jsbach_setup_grid(nproca, nprocb, npedim, model_id) ! Only one model implemented for now
  CALL jsbach_setup_tiles(model_id)
  CALL jsbach_init(model_id)

  CALL jsbach_label_run(jsbach_runs_standalone())

  ! Allocate arrays -> current assumption is for standalone nlon=nproma=npromz,
  !     therefore set (and allocated fields) here -- else need to strongly adapt
  !       -> i.e. set and allocate in the iblk loop for just one and change way of reading?!
  nproma = local_decomposition%nproma
  npromz = local_decomposition%npromz
  nblks = local_decomposition%ngpblks
  !TODO: get nproma/nblks with already existing function? 'get_nproma'
  !  -> hm, for that need to first get the model and then the grid
  WRITE(message_text, '(A,I5,I5,I5)') 'nproma, npromz, nblks = ', nproma, npromz, nblks
  CALL message('jsb4_driver', message_text)
  CALL message('', '')

  IF ( npromz .NE. nproma ) THEN
    CALL finish(TRIM(progname), &
      & 'Violation of assertion: for standalone in echam environment it is assumed that nlon=nproma=npromz, but nproma: '&
      & //int2string(nproma) //' and npromz: '//int2string(npromz)  )
  ENDIF

  ! forcing
  ALLOCATE(t_air_K(nproma, nblks))
  ALLOCATE(q_air(nproma, nblks))
  ALLOCATE(rain(nproma, nblks))
  ALLOCATE(snow(nproma, nblks))
  ALLOCATE(wind_air(nproma, nblks))
  ALLOCATE(wind_10m(nproma, nblks))
  ALLOCATE(lw_srf_down(nproma, nblks))
  ALLOCATE(swvis_srf_down(nproma, nblks))
  ALLOCATE(swnir_srf_down(nproma, nblks))
  ALLOCATE(swpar_srf_down(nproma, nblks))
  ALLOCATE(press_srf(nproma, nblks))
  ALLOCATE(drag_srf(nproma, nblks))
  ALLOCATE(t_acoef(nproma, nblks))
  ALLOCATE(t_bcoef(nproma, nblks))
  ALLOCATE(q_acoef(nproma, nblks))
  ALLOCATE(q_bcoef(nproma, nblks))
  ALLOCATE(pch(nproma, nblks))
  ALLOCATE(cos_zenith_angle(nproma, nblks))

  ALLOCATE(CO2_concentration(nproma, nblks))
  ALLOCATE(fract_par_diffuse(nproma, nblks))

  ! output required for next timesteps forcing (if not reading interface vars)
  ALLOCATE(t_srf_proc(nproma, nblks))
  ALLOCATE(fact_q_air_proc(nproma, nblks))
  ALLOCATE(fact_qsat_srf_proc(nproma, nblks))
  ALLOCATE(rough_h_srf_proc(nproma, nblks))
  ALLOCATE(rough_m_srf_proc(nproma, nblks))
  ALLOCATE(evapo_act2pot_proc(nproma, nblks))
!  ALLOCATE(evapo_act_sum_proc(nproma, nblks))
!  ALLOCATE(evapo_pot_sum_proc(nproma, nblks))

  ! output
  ALLOCATE(t_eff_srf(nproma))
  ALLOCATE(qsat_srf(nproma))
  ALLOCATE(s_srf(nproma))
  ALLOCATE(evapotrans(nproma))
  ALLOCATE(evapopot(nproma))
  ALLOCATE(latent_hflx(nproma))
  ALLOCATE(sensible_hflx(nproma))
  ALLOCATE(grnd_hflx(nproma))
  ALLOCATE(grnd_hcap(nproma))
  ALLOCATE(q_snocpymlt(nproma))
  ALLOCATE(alb_vis_dir(nproma))
  ALLOCATE(alb_nir_dir(nproma))
  ALLOCATE(alb_vis_dif(nproma))
  ALLOCATE(alb_nir_dif(nproma))
  ALLOCATE(CO2_flux(nproma))

  IF (.NOT. read_interface_vars) THEN
    CALL init_forcing(model_id, delta_time, time_step_len)

    !TODO: get from grid? but how to name such a function get_lon_from_grid?
    !      and if, shouldn't then also get_nproma should be called get_nproma_from_grid?
    lon => get_lon(local_decomposition)
    lat => get_lat(local_decomposition)

    ALLOCATE(coslon(nproma,nblks))
    ALLOCATE(coslat(nproma,nblks))
    ALLOCATE(sinlon(nproma,nblks))
    ALLOCATE(sinlat(nproma,nblks))

    coslon = COS(deg2rad * lon)
    sinlon = SIN(deg2rad * lon)
    coslat = COS(deg2rad * lat)
    sinlat = SIN(deg2rad * lat)

    !init for first timestep/day (according to jsb3: mo_jsbalone: init_driving)
    t_srf_proc = 273._wp
    fact_q_air_proc = 0.5_wp
    fact_qsat_srf_proc = 0.5_wp
    rough_h_srf_proc = 1._wp
    rough_m_srf_proc = 1._wp
    evapo_act2pot_proc = 1._wp

!    evapo_act_sum_proc = 0._wp
!    evapo_pot_sum_proc = 0._wp
  ENDIF

  model => get_model(model_id)
  CALL model%Get_top_tile(tile)
  ics = 1

  CALL init_mvstream('namelist.jsbach')
  CALL set_stream_element_nml('namelist.jsbach')
  CALL set_stream_nml('namelist.jsbach')

  IF(jsbach_is_restarted()) THEN

    CALL IO_read_streams

    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(TURB_)
    dsl4jsb_Get_memory(HYDRO_)

    DO iblk = 1, nblks
      IF (iblk == nblks) THEN
        ice = npromz
      ELSE
        ice = nproma
      END IF
      CALL model%Set_options(ics=ics, ice=ice, iblk=iblk)

      t_srf_proc(ics:ice,iblk) = dsl4jsb_var2D_onChunk(SEB_, t)
      fact_q_air_proc(ics:ice,iblk) = dsl4jsb_var2D_onChunk(TURB_, fact_q_air)
      fact_qsat_srf_proc(ics:ice,iblk) = dsl4jsb_var2D_onChunk(TURB_, fact_qsat_srf)
      rough_h_srf_proc(ics:ice,iblk) = dsl4jsb_var2D_onChunk(TURB_, rough_h)
      rough_m_srf_proc(ics:ice,iblk) = dsl4jsb_var2D_onChunk(TURB_, rough_m)
      ! TODO: if really required calculate daily sums of evapotrans and evapopot in jsbach and get them here (see comment below)
      evapo_act2pot_proc(ics:ice,iblk) = MIN(dsl4jsb_var2D_onChunk(HYDRO_, evapotrans),-EPSILON(1._wp)) &
        & / MIN(dsl4jsb_var2D_onChunk(HYDRO_, evapopot),  -EPSILON(1._wp))
    ENDDO

  ENDIF

  CALL jsbach_init_after_restart(model_id)

  CALL init_output
  CALL open_output_streams
  CALL mvstream_update_cache

  IF (timer_on()) CALL timer_start(timer_loop)

  ! - Main loop
  main_loop:DO

    CALL time_set   !(echam!)
    CALL manager_state(echam_time, current_datetime)
    ALLOCATE(current_datetime_ptr, source=current_datetime)
    ALLOCATE(next_datetime_ptr, source=next_date)
    IF(is_newday(current_datetime_ptr, delta_time)) THEN
      CALL write_date(current_date,'Calculate timesteps for day: ')
    ENDIF

    CALL jsbach_start_timestep(model_id, current_datetime_ptr, delta_time)

    IF (l_trigfiles) THEN
      CALL close_output_streams
      CALL open_output_streams
    ENDIF

    IF (.NOT. read_interface_vars) THEN

      CALL compute_cos_zenith_angle(sinlon, sinlat, coslon, coslat, cos_zenith_angle)

      IF (timer_on()) CALL timer_start(timer_forcing(1))
      CALL get_interface_variables_from_external_forcing(1, 1, nproma, npromz, nblks, current_datetime_ptr, next_datetime_ptr, &
        & sinlat, coslat, lon, evapo_act2pot_proc, &
        & t_srf_proc, fact_q_air_proc, fact_qsat_srf_proc, rough_h_srf_proc, rough_m_srf_proc, cos_zenith_angle,      &
        & CO2_concentration, t_air_K, q_air, rain, snow, wind_air, wind_10m,                                          &
        & lw_srf_down, swvis_srf_down, swnir_srf_down, swpar_srf_down, fract_par_diffuse,                              &
        & press_srf, drag_srf, t_acoef, t_bcoef, q_acoef, q_bcoef, pch)
      IF (timer_on()) CALL timer_stop(timer_forcing(1))
    ENDIF

    DO iblk = 1, nblks

      IF (iblk == nblks) THEN
        ice = npromz
      ELSE
        ice = nproma
      END IF
      CALL model%Set_options(ics=ics, ice=ice, iblk=iblk)

      IF (read_interface_vars) THEN
        IF (timer_on()) CALL timer_start(timer_forcing(1))
        CALL read_interface_variables( ics, ice, iblk, ncid, current_datetime_ptr, t_air_K(:,iblk), q_air(:,iblk), rain(:,iblk), &
          & snow(:,iblk), wind_air(:,iblk), wind_10m(:,iblk), lw_srf_down(:,iblk), swvis_srf_down(:,iblk), &
          & swnir_srf_down(:,iblk), swpar_srf_down(:,iblk), press_srf(:,iblk), drag_srf(:,iblk), &
          & t_acoef(:,iblk), t_bcoef(:,iblk), q_acoef(:,iblk), q_bcoef(:,iblk), pch(:,iblk), cos_zenith_angle(:,iblk))
        IF (timer_on()) CALL timer_stop(timer_forcing(1))
      ENDIF

      CALL jsbach_interface ( 1, iblk, ics, ice, current_datetime_ptr, delta_time, time_step_len, t_air_K(ics:ice,iblk),       &
        & q_air(ics:ice,iblk), rain(ics:ice,iblk), snow(ics:ice,iblk), wind_air(ics:ice,iblk), wind_10m(ics:ice,iblk),         &
        & lw_srf_down(ics:ice,iblk), swvis_srf_down(ics:ice,iblk), swnir_srf_down(ics:ice,iblk), swpar_srf_down(ics:ice,iblk), &
        & fract_par_diffuse(ics:ice,iblk), press_srf(ics:ice,iblk), drag_srf(ics:ice,iblk),                                    &
        & t_acoef(ics:ice,iblk), t_bcoef(ics:ice,iblk), q_acoef(ics:ice,iblk), q_bcoef(ics:ice,iblk), pch(ics:ice,iblk),       &
        & cos_zenith_angle(ics:ice,iblk), CO2_concentration(ics:ice,iblk),                                                     &
        & t_srf_proc(ics:ice,iblk), t_eff_srf(ics:ice), qsat_srf(ics:ice), s_srf(ics:ice),                                     &
        & fact_q_air_proc(ics:ice,iblk), fact_qsat_srf_proc(ics:ice,iblk),                                                     &
        & evapotrans(ics:ice), latent_hflx(ics:ice), sensible_hflx(ics:ice), grnd_hflx(ics:ice), grnd_hcap(ics:ice),           &
        & rough_h_srf_proc(ics:ice,iblk), rough_m_srf_proc(ics:ice,iblk), q_snocpymlt(ics:ice),                                &
        & alb_vis_dir(ics:ice), alb_nir_dir(ics:ice),                                                                          &
        & alb_vis_dif(ics:ice), alb_nir_dif(ics:ice), CO2_flux(ics:ice), evapopot=evapopot(ics:ice))

      IF (.NOT. read_interface_vars) THEN
        ! TODO: if really required (only ERA-interim forcing data):
        ! To ensure restart identity, the daily evapo_act_sum and evapo_pot_sum would need to be calculated in jsbach!
        ! => for now the instantaneous values are used instead of the daily sums, to ensure restart identity
        evapo_act2pot_proc(ics:ice,iblk) = MIN(evapotrans(ics:ice),-EPSILON(1._wp)) / MIN(evapopot(ics:ice),  -EPSILON(1._wp))
!        evapo_act_sum_proc(:,iblk) = evapo_act_sum(:,iblk) + MIN(evapotrans(:),-EPSILON(1._wp))
!        evapo_pot_sum_proc(:,iblk) = evapo_pot_sum(:,iblk) + MIN(evapopot(:),  -EPSILON(1._wp))
      ENDIF
    ENDDO

!    IF ((.NOT. read_interface_vars) .AND. is_newday(current_datetime_ptr, delta_time)) THEN
!      evapo_act2pot_proc = evapo_act_sum_proc / evapo_pot_sum_proc
!      evapo_act_sum_proc = 0._wp
!      evapo_pot_sum_proc = 0._wp
!    ENDIF

    CALL jsbach_finish_timestep(model_id, current_datetime_ptr, time_step_len)
    CALL deallocateDatetime(current_datetime_ptr)
    CALL deallocateDatetime(next_datetime_ptr)
    CALL mvstream_accumulate
    IF (timer_on()) CALL timer_start(timer_output)
    CALL out_streams
    IF (timer_on()) CALL timer_stop(timer_output)

    IF (l_putrerun) THEN
      IF (timer_on()) CALL timer_start(timer_restart)
      CALL write_streams
      IF (timer_on()) CALL timer_stop(timer_restart)
    END IF

    CALL time_reset
    IF (lstop) EXIT main_loop

  ENDDO main_loop

  IF (timer_on()) CALL timer_stop(timer_loop)

  IF (read_interface_vars) THEN
    CALL nf_check(nf_close(ncid))
  ENDIF

  CALL close_output_streams
  CALL message('INFO','Closed output streams')

  CALL finalize_external_forcing()
  CALL jsbach_finalize()

  IF (timer_on()) CALL cleanup_timer

  !$OMP PARALLEL
!$OMP MASTER
  status = util_cputime(zutime, zstime)
!$OMP END MASTER
!$OMP END PARALLEL
  IF (status == -1) THEN
    CALL message(TRIM(progname),'Cannot determine used CPU time')
  ELSE
!$OMP PARALLEL
!$OMP MASTER
    zwtime = util_walltime(0)
!$OMP END MASTER
!$OMP END PARALLEL
    zrtime = (zutime+zstime)/zwtime
    CALL message ('', '')
    WRITE (message_text,'(a,f10.2,a)') ' Wallclock        : ', zwtime, ' s'
    CALL message('',message_text)
    WRITE (message_text,'(a,f10.2,a)') ' CPU-time (user)  : ', zutime, ' s'
    CALL message('',message_text)
    WRITE (message_text,'(a,f10.2,a)') ' CPU-time (system): ', zstime, ' s'
    CALL message('',message_text)
    WRITE (message_text,'(a,f10.2,a)') ' Ratio            : ', 100*zrtime, ' %'
    CALL message('',message_text)
    CALL message ('', '')
  END IF

  ! deallocate in- and output arrays
  DEALLOCATE(t_air_K, q_air, rain, snow, wind_air, wind_10m)
  DEALLOCATE(lw_srf_down, swvis_srf_down, swnir_srf_down, swpar_srf_down, fract_par_diffuse)
  DEALLOCATE(press_srf, drag_srf, t_acoef, t_bcoef, q_acoef, q_bcoef, pch, cos_zenith_angle, CO2_concentration)

  DEALLOCATE(t_eff_srf, qsat_srf, s_srf)
  DEALLOCATE(evapotrans, latent_hflx, sensible_hflx, grnd_hflx, grnd_hcap)
  DEALLOCATE(q_snocpymlt, alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif)
  DEALLOCATE(CO2_flux)

  IF (.NOT. read_interface_vars) THEN
    DEALLOCATE(lon, lat, coslon, coslat, sinlon, sinlat)
    DEALLOCATE(t_srf_proc, fact_q_air_proc, fact_qsat_srf_proc, rough_h_srf_proc, rough_m_srf_proc, evapo_act2pot_proc)
    !DEALLOCATE(evapo_act_sum_proc, evapo_pot_sum_proc)
  ENDIF

  ! Stop MPI
  CALL p_stop
#endif
#endif
END PROGRAM jsb4_driver
