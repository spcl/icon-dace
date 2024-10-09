! get_dates.f90 - Get the dates from the time controller
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!_________________________________________

  SUBROUTINE get_dates(IO_timestep, istep)

    ! This routine originates (year 2014) from MPI-ESM, the Earth System Model of the 
    ! Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
    ! Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
    ! version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
    ! doi: 10.1029/2018MS001400.
    !
    ! For standalone mode only.
    ! Corresponds to ECHAM5 subroutine IO_init in mo_io.f90
    ! This routine doesn't need to be called from ECHAM in coupled mode.

    USE mo_io, ONLY: IO_READ, IO_open, IO_close
    USE mo_mpi,           ONLY : p_io, p_parallel_io
    USE mo_netCDF, ONLY: FILE_INFO, NF_GLOBAL, io_get_att_int
    USE mo_time_control, ONLY: INIT_STEP, delta_time, dt_start, &
                               lresume, resume_date, start_date, inp_convert_date, write_date
    USE mo_time_conversion, ONLY : time_native, TC_set, TC_convert
    USE mo_exception,       ONLY : message, finish

    INTEGER, INTENT(out) :: IO_timestep, istep

    TYPE(FILE_INFO) :: IO_file
    INTEGER :: IO_file_id
    INTEGER :: forecast_date, verification_date   ! YYYYMMDD 
    INTEGER :: forecast_time, verification_time   ! HHMMSS
    TYPE(time_native) :: date_nat
    LOGICAL :: debug = .true.                     !OSBRSM

    IF (debug) CALL message('jsbach_init_io','BEGIN')

    ! Initialize time
    IF (p_parallel_io) THEN
       IF (lresume) THEN
          ! Get time information from restart file
          IO_file%opened = .FALSE.
!OSBRSM          CALL IO_open(trim(theOptions%RestartPrefix)//'_jsbach', IO_file, IO_READ)
!OSBRSM          IO_file_id = IO_file%file_id
!OSBRSM          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'nstep', istep)
!OSBRSM          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'timestep', IO_timestep)
!OSBRSM          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'fdate', forecast_date)
!OSBRSM          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'ftime', forecast_time)
!OSBRSM          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'vdate', verification_date)
!OSBRSM          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'vtime', verification_time)
!OSBRSM          CALL IO_close(IO_file)
       ELSE
          istep = INIT_STEP
          IO_timestep = delta_time
       ENDIF
    ENDIF
       
!OSBRSM    IF (p_parallel) THEN
!OSBRSM       CALL p_bcast(istep, p_io)
!OSBRSM       CALL p_bcast(IO_timestep, p_io)
!OSBRSM       IF (lresume) THEN
!OSBRSM          CALL p_bcast(forecast_date, p_io)
!OSBRSM          CALL p_bcast(forecast_time, p_io)
!OSBRSM          CALL p_bcast(verification_date, p_io)
!OSBRSM          CALL p_bcast(verification_time, p_io)
!OSBRSM       ENDIF
!OSBRSM    ENDIF
       
!   IF (lresume) THEN
!      CALL inp_convert_date(forecast_date, forecast_time, start_date)
!      CALL inp_convert_date(verification_date, verification_time, resume_date)
!   ELSE
       IF (SUM(dt_start(:)) /= 0) THEN
          CALL TC_set(dt_start(1),dt_start(2),dt_start(3), &
               dt_start(4),dt_start(5),dt_start(6), date_nat)
          CALL TC_convert(date_nat, start_date)
          resume_date = start_date
       ELSE
          CALL finish('jsbach_init_io', 'Start date not set')
       ENDIF
!   ENDIF

!   IF (p_parallel_io) THEN
       CALL write_date(start_date, 'Start date (initial/restart): ')
       CALL write_date(start_date, 'Resume date (initial/restart): ')
!   ENDIF

    IF (debug) CALL message('get_dates','END')

  END SUBROUTINE get_dates
