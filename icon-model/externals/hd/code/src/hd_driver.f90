! hd_driver.f90 - HD main program that steers the simulation
! 
! Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: Apache-2.0
! See ./LICENSES/ for license information
!
! Authors: Stefan Hagemann and Ha Ho-Hagemann
! Contact: <stefan.hagemann@hereon.de>
!_________________________________________

PROGRAM hd_driver
  
  !------------------------------------------------------------------------------------
  !
  !   ******* Globale/Regionale Abfluss-Simulation: Hauptprogramm, dessen
  !           Hauptaufgabe darin besteht, die HDMODEL-Subroutine 
  !           HDMODEL.f aufzurufen und die Inputfelder pro Zeitschritt zu 
  !           uebergeben.
  !
  !
  !   ******** Version 1.0 - Oktober 1999 
  !            Programmierung und Entwicklung: Stefan Hagemann 
  !            Programmcode basiert auf der Offline-Version der HDModels
  !            die auch regional anwendbar ist (regsim.f)
  !            Service-Outputfile ist nun simple binaer statt direct access
  !
  !   ******** Version 1.1 - March 2000
  !            Steuerungsroutine: Inputdaten ueber Commonblock pinp.for
  !            gesteuert 
  !
  !   ******** Version 1.1 -- January 2001
  !            Ursprung Ocean-Longitude korrigiert.  0. statt -1.40625
  !
  !            Anmerkung: Input-Daten von Runoff und Drainage sollten
  !                       Einheit m/s haben.
  !
  !   ******** Version 1.2 -- October 2005
  !            Implementation of Unit Factor UFAKRU, which is applied to Input arrays
  !            of Surface ruonoff and drainage if UFAK.ne.1. Init factor is necessary
  !            if runoff/drainage input unit is not m/s.
  !
  !   ******** Version 1.2.1 -- August 2008
  !            Initialisierung von UFAKRU mit 0
  !            Einbau Zusaetzlicher Logpunkte fuer ISOLOG
  !            ISOLOG = 100 -> Koordinaten von Fluss 1 und 2 werden aus Input 
  !            file hdini.inp gelesen.
  !
  !   ******** May 2011, Veronika Gayler
  !            fortran 90 version to run with ECHAM6 mo_hydrology
  !
  ! Note that the 0.5° HD model version is also part of MPI-ESM, the Earth System Model of the 
  !    Max Planck Institute for Meteorology (Mauritsen et al. 2019), where it is called
  !    as a subroutine of the JSBACH land surface model (Reick et al. 2021). 
  !    References:
  !    Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
  !      version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
  !      doi: 10.1029/2018MS001400.
  !    Reick, C., et a. (2021) JSBACH 3 - The land component of the MPI Earth System Model:
  !      Documentation of version 3.2. Berichte zur Erdsystemforschung, 240, 
  !      Max Planck Institute for Meteorology, Hamburg, http://doi.org/10.17617/2.3279802.
  !
  ! O. Böhringer, MPI-M, 2014, Separation of offline version from MPI-ESM infrastructure. 
  ! H. Ho-Hagemann, HZG, June 2014, introducing OASIS coupling in offline version.
  ! S. Hagemann, MPI-M, 2014, further refinements of offline version.
  ! H. Ho-Hagemann, HZG, 2018, Update of OASIS coupling.
  ! S. Hagemann, HZG, 2018, Refinements of HD 5 Min. version. 
  !
  !------------------------------------------------------------------------------------

  USE mo_kind,             ONLY: dp 
  USE mo_io_units
  USE mo_hydrology,        ONLY: init_hydrology, hydrology_model, redistribute_sinks, &
                                 hydrology_restart, cleanup_hydrology, &
                                 locean, lbase, hd_steps_per_day, &
                                 grid_hd, friv, lhd_highres, nremap ,&
                                 water_to_ocean, hd_area, diag_water_budget

  USE mo_hd_highres_io,    ONLY: hd_highres_init, hd_highres_open, hd_highres_close, hd_highres_write
  USE mo_mpi,              ONLY: p_start, p_stop, p_init_communicators, p_nprocs
  USE mo_machine,          ONLY: machine_setup
  USE mo_time_control,     ONLY: l_trigfiles, dt_start, init_manager, delta_time, &
                                 init_times, ec_manager_init, time_reset, &
                                 no_steps, time_set, current_date, &
                                 write_date, day_difference
  USE mo_exception,        ONLY: finish, message, message_text
!OSBJSB  USE mo_jsbach_interface, ONLY: get_dates
  USE mo_control,          ONLY: nlon, ngl, nproca, nprocb, nprocio
  USE mo_filename,         ONLY: out_expname, out_datapath, find_next_free_unit
  USE mo_netcdf,           ONLY: file_info, nf_max_name, io_inq_dimid, io_inq_dimlen, &
                                 io_inq_varid, io_get_var_double, io_get_vara_double
  USE mo_io,               ONLY: io_open, io_close, io_read
  USE mo_grid,             ONLY: domain, areacalc, read_grid_info


  USE mo_couple_to_ocean,  ONLY: fmou_hd_to_ocean,    &      ! HD mouth points that have a target on ocean grid
                                 nxocean, nyocean, discharge_on_ocean,    &
                                 read_coupling_info, dis_to_ocean, &
                                 hd_on_ocean_open, hd_on_ocean_write, hd_on_ocean_close

  USE mo_coupling,         ONLY: is_coupled_run, runoff_s, runoff_dr, &
                                 set_coupling_type, get_coupling_type, &
                                 lcoupling_atm, lcoupling_oce, icpl_sinks, icpl_mask_tohd
  USE mo_coupling_hd,      ONLY: coupling_hd_init, &
                                 coupling_hd_recv_from_land, &
                                 coupling_hd_send_to_ocean, &
                                 coupling_hd_send_to_ocean_direct, &
                                 hd_outflow
  USE mo_bias_correction,  ONLY: bc_init, bc_cleanup, bias_correction, hd_bc_outflow

  IMPLICIT NONE

  ! local variables
  INTEGER  :: IO_timestep
  INTEGER  :: luoce          ! file unit
  INTEGER  :: idum           ! dummy integer
  INTEGER  :: step           ! hd model time step counter
  INTEGER  :: istep          ! initial time step

  ! local variables defined in namelist HD_CTL
  INTEGER  :: nstep
  REAL(dp) :: ufakru
  INTEGER  :: year1, month1
  CHARACTER(LEN=10)  :: date_start, date_end
  INTEGER  :: forcing_freq
  INTEGER  :: iout
  LOGICAL  :: lcoupling_out   ! Write output for coupling_type 2 (no/yes)
  INTEGER  :: iform_input     ! Format Input files: 0 = SRV, 1 = NetCDCF
  LOGICAL  :: ltransport      ! Switch for Transport on/off (default .false.)   
  INTEGER  :: ibc_type        ! Bias correction type: 0=None, 1 = Mean Bias, 2 = Low, Mid and High Biases        
  CHARACTER(LEN=120) :: dn_bcpara  ! File with Bias correction parameters
  LOGICAL  :: lbc_write       ! Switch for writing bias corrected discharges to file (default .false.)   

  ! local arrays that have to do with the forcing, defined on the atmosphere grid
  REAL(dp), ALLOCATABLE :: slm(:,:), slf(:,:), alake(:,:), glac(:,:)
  REAL(dp), ALLOCATABLE :: runoff(:,:), drain(:,:)
  REAL(dp), ALLOCATABLE :: disch(:,:)
  REAL(dp), ALLOCATABLE :: awfre(:,:), apmecal(:,:)

  REAL,     ALLOCATABLE :: dummy_sp(:,:)
  REAL,     ALLOCATABLE :: dummy_sp_CCLM(:,:)
 
  ! local array to store the accumulated HD discharge
  REAL(dp), ALLOCATABLE :: hd_discharge_accu(:,:)

  ! local arrays needed with service format (.srv)
  INTEGER :: IHEAD(8)

  ! logical I/O units
  INTEGER :: lurun   ! Unit of surface runoff forcing file if reading SRV
  INTEGER :: lubas   ! Unit of drainage/subsurface runoff forcing file if reading SRV
  INTEGER  :: discharge_file_id, meanflow_file_id, bcflow_file_id

  ! file names
  CHARACTER(nf_max_name) :: runoff_file    ! input file with runoff data
  CHARACTER(nf_max_name) :: drainage_file  ! input file with drainage data
  CHARACTER(80) :: dnout          ! averaged output data (nml parameter iout)
  CHARACTER(80) :: coupling_file  ! input file with coupling info for coupling type 2

  TYPE(file_info) :: runofffile    ! input file with runoff data
  TYPE(file_info) :: drainfile     ! input file with drainage data

  ! Grid infos
  TYPE(domain)  :: grid_forcing    ! Grid info of forcing data taken from mask file in init_forcing

  ! variable IDs
  INTEGER :: runoffid, drainid

  ! local parameters
  INTEGER, PARAMETER :: DAILY = 1    ! used for frequency of forcing data
  INTEGER, PARAMETER :: STEPWISE = 0
  REAL(dp), PARAMETER :: zeps = 1.E-10

  INTEGER :: i,j

  !----------------
  ! initialization
  !----------------

write(nout,*) 'OTBp_start '

  ! MPI Initialisation
  CALL p_start('HD_driver')

  IF (p_nprocs /= 1) THEN
    CALL finish('hd_driver', ' HD currently only support a single process.')
  ELSE
    nproca = 1
    nprocb = 1
    nprocio = 0
  END IF

WRITE(nout,*) 'nproca=', nproca,' nprocb=', nprocb,' nprocio=', nprocio
write(nout,*) 'OTBp_init_communicators '

  CALL p_init_communicators(nproca, nprocb, nprocio)

write(nout,*) 'OTBmachine_setup '
  CALL machine_setup

  ! read the namelist

  CALL config_hd

  ! initialization of the echam time manager

  dt_start = (/year1,month1,1,0,0,0/)           ! start date: Jan 1st of year1
  no_steps = nstep                              ! number of time steps within the run

write(nout,*) 'OTBget_dates ', no_steps, ' Start year from namelist: ', year1
  CALL get_dates(IO_timestep, istep)            ! get timestep
write(nout,*) 'OTBec_manager_init '
  CALL ec_manager_init(IO_timestep, istep)      ! time manager initialization
write(nout,*) 'OTBinit_manager '
  CALL init_manager
write(nout,*) 'OTBinit_times '
  CALL init_times                               ! initialize all dates
  l_trigfiles = .FALSE.

  ! hd model initializations
 
write(nout,*) 'OTBhd_init_dims ', l_trigfiles
  CALL hd_init_dims                             ! init model grid dimensions
write(nout,*) 'OTBhd_init_forcing'
  CALL hd_init_forcing                          ! allocate memory of the forcing arrays

!---- Ha Ho-Hagemann }
write(nout,*) 'OTBinit_hydrology'
  CALL init_hydrology(slm, alake)
  CALL hd_highres_init(grid_hd)
  ALLOCATE (hd_discharge_accu(grid_hd%nlon,grid_hd%nlat))

  IF (get_coupling_type() .EQ. 2) THEN
    write(nout,*) 'OTBread_coupling_info'
    CALL read_coupling_info(coupling_file)
  ENDIF

  ! Bias correction requested?
  IF (ibc_type.GT.0) THEN
    CALL bc_init(dn_bcpara)
    write(nout,*) '+++++ HD: hd_driver: Bias correction initialized!'
  ENDIF  

write(nout,*) 'OTBhd_init_io'
  CALL hd_init_io                               ! open in- and output files

  write(0,*) '+++++ HD: hd_driver.f90: istep =',istep,' nstep = ', nstep

  IF (is_coupled_run()) CALL coupling_hd_init()

  write(0,*) '+++++ HD: hd_driver.f90: START Time Loop:'

 
  !----------------
  ! Time step loop
  !----------------

  DO step = istep+1, istep+nstep

    ! set time for echam time manager
    CALL time_set


    IF ((forcing_freq == STEPWISE) .OR. &
         forcing_freq == DAILY .AND. MOD(step-1,hd_steps_per_day) == 0) THEN
      ! read forcing data
      CALL write_date(current_date, 'Read forcing data of :')
      IF (.NOT. is_coupled_run() .OR. .NOT. lcoupling_atm) THEN
        IF (nremap.NE.0 .AND. nremap.NE.3) THEN
          CALL hd_update_forcing                 !--> runoff(nlon,ngl) & drain (nlon,ngl)
        ELSE
#if !defined(COUP_OAS) && !defined(COUP_YAC)
          CALL hd_update_forcing_noremap         !--> runoff_s(nl,nb) & runoff_dr (nl,nb)
#endif
        ENDIF
      END IF
     END IF

     IF (step.EQ.istep+1) write(0,*) '+++++ HD: hd_driver.f90: 1st forcing data read'

!!!     write(0,*) '+++++ HD: hd_driver.f90: step =',step,' coupling_hd_recv_from_land'
     IF (is_coupled_run() .AND. lcoupling_atm) THEN 
       CALL coupling_hd_recv_from_land()
       IF (ABS(ufakru-1.) .GT. 1.E-6) THEN
         runoff_s(:,:) = runoff_s(:,:) * ufakru
         runoff_dr(:,:) = runoff_dr(:,:) * ufakru
       ENDIF
     ENDIF

     ! runoff_s(nl,nb): runoff (~gl_aros in ECHAM (nlon,ngl), that needs to pass hd_remap in mo_hydrology)
     ! runoff_dr(nl,nb): drain  (~gl_adrain in ECHAM (nlon,ngl), similar to gl_aros)

!     print*, "+++ HD1: MIN, MAX of runoff_s=", MINVAL(runoff_s), MAXVAL(runoff_s)
!     print*, "+++ HD1: MIN, MAX of runoff_dr=", MINVAL(runoff_dr), MAXVAL(runoff_dr)


!     print*, "+++ HD2: MIN, MAX of runoff_s=", MINVAL(runoff_s), MAXVAL(runoff_s)
!     print*, "+++ HD2: MIN, MAX of runoff_dr=", MINVAL(runoff_dr), MAXVAL(runoff_dr)

     ! discharge calculations
     CALL hydrology_model(slm, alake, glac, runoff, drain, disch, &
          awfre, apmecal)

     ! *** Note that water_to_ocean include discharge/inflow on gridboxes that are either ocean
     ! *** or a land sink, and, hence, this water leaves HD (and is not handled anymore) and needs be
     ! *** put into the ocean (friv considers the discharge on land boxes that are no sink)

     IF (ibc_type.GT.0) THEN      ! Bias correction on ocean points only, sinks are not modified.
       CALL bias_correction(ibc_type, water_to_ocean)
     ENDIF

     IF (is_coupled_run() .AND. lcoupling_oce) THEN
       ! Prepare HD output field for coupling, possible with redistributing of water entering internal sinks
       IF (ibc_type.EQ.0) THEN
         CALL redistribute_sinks(icpl_sinks, water_to_ocean, hd_outflow) 
       ELSE
         CALL redistribute_sinks(icpl_sinks, hd_bc_outflow, hd_outflow) 
       ENDIF
       IF (diag_water_budget) THEN
         WRITE (message_text,*) 'DIAGWB: Global discharge into the ocean before sending: ', &
                 sum(hd_outflow), ' m3/s'
         CALL message ('hd_driver',message_text)
       ENDIF
     ENDIF

     IF (get_coupling_type() .EQ. 2) THEN
       ! TODO: This may be changed if water from non-mouth ocean boxes shall be regarded, e.g. 
       !       coming from atmospheric input directly or from redistribution type icpl_sinks = 2
       WHERE (fmou_hd_to_ocean.LE.0.5)
         hd_outflow(:,:) = 0._dp       
       END WHERE
!    
!      ******** Transfer HD model river discharge to ocean grid
       CALL dis_to_ocean(hd_outflow)
       IF (ABS(SUM(hd_outflow) - SUM(discharge_on_ocean)).GT. 0.01_dp) THEN
         WRITE (message_text,*) 'Discharge sums differ between HD and ocean grid: ', &
               SUM(hd_outflow), ' != ', SUM(discharge_on_ocean)
         CALL message('hd_driver', message_text)
         CALL finish ('hd_driver', 'run terminated.')
       ENDIF
     ENDIF

!!!     write(0,*) '+++++ HD: hd_driver.f90: step =',step,' coupling_hd_send_to_ocean'
    IF (is_coupled_run() .AND. lcoupling_oce) THEN
      IF (get_coupling_type() .EQ. 2) THEN
        WRITE(message_text,*) 'Direct transfer of discharge_on_ocean via OASIS!'
        CALL coupling_hd_send_to_ocean_direct()
      ELSE
        CALL coupling_hd_send_to_ocean(hd_outflow)
      ENDIF
    END IF

!     write(0,*) '+++++ HD: hd_driver.f90: step =',step,' friv = ', friv + water_to_ocean
     
     
     ! write the output
     CALL hd_write_output

     ! write restart file
     IF (step == istep+nstep) CALL hydrology_restart 

     ! update model time step
     CALL time_reset

  ENDDO

  !   *** Schreiben der letzten Mittelwert-Daten fuer Inflow per Gridbox 
  IDUM=-1
  CALL GMITWRI(FRIV + water_to_ocean, hd_discharge_accu, IDUM)

  !-------------
  ! cleaning up
  !-------------

  CALL cleanup_hydrology
  DEALLOCATE(hd_discharge_accu)
  IF (ibc_type.GT.0) CALL bc_cleanup

  ! close files

  CALL hd_highres_close(meanflow_file_id)
  IF (locean) CALL hd_highres_close(discharge_file_id)
  IF (lbc_write) CALL hd_highres_close(bcflow_file_id)

  IF (.NOT. is_coupled_run() .OR. .NOT.lcoupling_atm) THEN
    IF (iform_input.EQ.0) THEN
      CLOSE(lurun)
      IF (lbase) CLOSE(lubas)
    ELSE IF (iform_input.EQ.1) THEN
      CALL IO_close(runofffile)
      IF (lbase) CALL IO_close(drainfile)
    ENDIF
  ENDIF

  IF (lcoupling_out) CALL hd_on_ocean_close

  ! MPI finalization

   call p_stop

CONTAINS
!
!****************************************************************************
  SUBROUTINE GMITWRI(FOUT, fsum, IOUT)
!****************************************************************************
!
!     ******** Routine zur Mittelung ueber NT Zeitschritte und Ausgabe
!              in einer Serviceformat-Binaerdatei.c 
!
!     ******** Version 1.0 - Oktober 1995
!              Programmierung und Entwicklung: Stefan Hagemann 
!
!     ******** Version 1.1 - Oktober 1996
!              Nun auch Monatsmittelwerte mit Schaltjahren moeglich (IOUT=5)
!
!     ******** Version 1.2 - Oktober 1999
!              Serviceformat: simple binaer statt direct access.
!
!     ******** Version 1.3 - May 2002
!              Daily Output possible --> IOUT = 6
!
!     ******** Version 1.4 - August 2018
!              Output NetCDF instead of SRV
!
!     ********* Variablenliste
!     ***
!     ***   FOUT = Globales Output-Array des Timesteps IREC
!     ***     NL = Anzahl der Laengengrade, now grid_hd%nlon
!     ***     NB = Anzahl der Breitengrade, now grid_hd%nlat
!     ***   LDEBUGHD = Kommentarvariable ( 0 = Kein Kommentar )
!     ***   IOUT = Mittelungsartvariable, d.h. ueber wieviel Zeitschritte
!     ***          1   30-Day Averages   --> NT = 30 * hd_steps_per_day
!     ***          2   Decadal Averages  --> NT = 10 * hd_steps_per_day
!     ***          3   Weekly Averages   --> NT = 7  * hd_steps_per_day
!     ***          4   Monthly Averages ohne Schaltjahre
!     ***          5   Monthly Averages inklusive Schaltjahre
!     ***          6   Daily Output
!     ***          -1  Mittelung ueber restliche Zeitschritte und Return
!     ***
!     ***   IMON = Monatsnummer bzw. New Record Ja/Nein-Variable 
!     ***          -1 ==> Beginn eines New Records
!     ***   it = Lokaler Timestep-Zaehler des Records IREC+1
!     ***
!     ******** Include of Parameter NL, NB 
!     ***       NL, NB = Globale/Regionale Feldgrenzen

    USE mo_time_control,     ONLY: get_date_components, get_month_len, current_date

    REAL(dp), INTENT(in)     :: fout(grid_hd%nlon,grid_hd%nlat)
    REAL(dp), INTENT(inout)  :: fsum(grid_hd%nlon,grid_hd%nlat)
    INTEGER,  INTENT(in)     :: iout

    INTEGER,  SAVE :: nt                          ! number of time steps
    INTEGER        :: yr, mon, day, hr, min, sec  ! current date components
    INTEGER,  SAVE :: nday(12)
    INTEGER,  SAVE :: it
    LOGICAL,  SAVE :: start_accumulation = .TRUE.
    INTEGER,  SAVE :: counter_m = 0


    ! initializations
    nday = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    
    ! find out year and month
    CALL get_date_components(current_date, yr, mon, day, hr, min, sec)

    IF (iout == -1) THEN        ! only at the very end
       IF (start_accumulation) RETURN   ! nothing to do, averaging had just been done

       ! averaging the remaining time steps
       fsum(:,:) = fsum(:,:) / it

       CALL hd_highres_write(meanflow_file_id, fsum, counter_m, grid_hd%nlon, grid_hd%nlat)

       IF (lcoupling_out) THEN
          hd_outflow(:,:) = 0._dp
          WHERE (fmou_hd_to_ocean.GT.0.5)
            hd_outflow(:,:) = fsum(:,:)
          END WHERE
          CALL dis_to_ocean(hd_outflow)
          CALL hd_on_ocean_write(discharge_on_ocean)
       ENDIF

       start_accumulation = .TRUE.
       RETURN
    ENDIF


    IF (start_accumulation) THEN

       ! initializations
       it = 0
       fsum(:,:) =0.
       start_accumulation = .FALSE.

       ! find number of time_steps for the averaging

       IF (iout == 1) THEN

          ! 30-day averages
          nt = 30 * hd_steps_per_day

       ELSE IF (iout == 2) THEN

          ! 10-day averages
          nt = 10 * hd_steps_per_day

       ELSE IF (iout == 3) THEN

          ! weekly
          nt = 7 * hd_steps_per_day

       ELSE IF (iout == 4) THEN

          ! monthly averages, without leap year
          nt = nday(mon) * hd_steps_per_day

       ELSE IF (iout == 5) THEN

          ! monthly averages, with leap year
          nt = get_month_len(yr,mon) * hd_steps_per_day

       ELSE IF (iout == 6) THEN

          ! daily averages
          nt = hd_steps_per_day

       ELSE IF (iout == 7) THEN

          ! stepwise output (each time routine is called)
          nt = 1

       ELSE
          WRITE (message_text,*) 'namelist parameter IOUT out of range: iout = ', iout
          CALL finish ('hd_gmitwri', message_text)
       ENDIF
    ENDIF

    ! accumulate data for the selected time period

    it = it + 1
    fsum(:,:) = fsum(:,:) + fout(:,:)
    IF (it == nt) THEN

       ! do the averaging
       fsum(:,:) = fsum(:,:) / nt

       ! Write areas at first time step
       IF (counter_m.EQ.0) THEN
         CALL hd_highres_write(meanflow_file_id, fsum, counter_m, grid_hd%nlon, grid_hd%nlat, hd_area)
       ELSE
         CALL hd_highres_write(meanflow_file_id, fsum, counter_m, grid_hd%nlon, grid_hd%nlat)
       ENDIF

       IF (lcoupling_out) THEN
          hd_outflow(:,:) = 0._dp
          WHERE (fmou_hd_to_ocean.GT.0.5)
            hd_outflow(:,:) = fsum(:,:)
          END WHERE
          CALL dis_to_ocean(hd_outflow)
          CALL hd_on_ocean_write(discharge_on_ocean)
       ENDIF

       start_accumulation = .TRUE.
    ENDIF

  END SUBROUTINE GMITWRI

  !------------------------------------------------------------------------------------
  ! Basic configuration of HD offline simulations
  !------------------------------------------------------------------------------------
  SUBROUTINE config_hd

    !----------------------------------------------------------------------------------
    ! parameters of namelist HD_CTL
    !
    !   out_expname    experiment name
    !  out_datapath    path to where the output data shall be written
    !         year1    initial year of the run
    !        month1    initial month of the run
    !    date_start    start date of the run, format YYYYMMDD or YYYY-MM-DD
    !      date_end    end date of the run, format YYYYMMDD or YYYY-MM-DD
    !         nstep    number of time steps within the run if date_start & date_end are not provided
    !    delta_time    model time step lenght in seconds
    !   runoff_file    file with input runoff data
    ! drainage_file    file with input drainage data
    !        ufakru    unit factor for runoff and drainage input data
    !  forcing_freq    data frequency of the forcing data (STEPWISE/DAILY)
    !          iout    averaging period of some HD output
    ! coupling_type    Coupling switch that also indicates type of coupling to oceanmodel 
    !                  0=no coupling, 1= with no interpol. in HD, 2=with direct alloc.
    ! lcoupling_atm    Switch for coupling to atmosphere (default: .False.) 
    ! lcoupling_oce    Switch for coupling to ocean (default: .False.)
    ! icpl_sinks       Redistribution of water in sinks for ocean coupling 
    !                  0=None (Def.), 1=Relatively to mouth boxes, 2= equally to all ocean boxes
    ! icpl_mask_tohd   Switch for the mask on which HD may receive input for atmosphere coupling 
    !                  0=HD land boxes (Def.), 1=all HD boxes
    ! coupling_file    input file with coupling information for coupling_type 2
    ! lcoupling_out    Write discharge on ocean grid (no/yes) (coupling_type 2 only)
    !   iform_input    Format Input files: 0 = SRV, 1 = NetCDCF
    !    ltransport    Switch for Transport on/off (default .false.)   
    ! ibc_type         Bias correction type: 0=None, 1 = Mean Bias, 2 = Low, Mid and High Biases 
    ! dn_bcpara        File with Bias correction parameters
    ! lbc_write        Switch for writing bias corrected discharges to file (default .false.) 
    !----------------------------------------------------------------------------------

    USE mo_namelist,         ONLY: open_nml, position_nml, POSITIONED
    USE mo_io_units,         ONLY: nout
 
    ! local variables

    INTEGER                :: read_status, inml, iunit, ndays, day1
    INTEGER                :: coupling_type

    INCLUDE 'hd_ctl.inc'

    ! set default values of the namelist parmeters

    out_expname = 'hd'
    out_datapath = './'
    year1 = 1900
    month1 = 1
    date_start = ''
    date_end = ''
    nstep = 365
    delta_time = 86400._dp        ! time step in seconds (one day) 
    ufakru = 1._dp
    runoff_file = "runoff.nc"
    drainage_file = "drainage.nc"
    forcing_freq = STEPWISE
    iout = 5               ! averaging of the HD output: 1=30d, 2=10d, 3=7d, 5=monthly, 6=daily
    coupling_type = 0
    lcoupling_atm = .false.
    lcoupling_oce = .false.
    icpl_sinks = 0
    icpl_mask_tohd = 0
    coupling_file = "hdcouple.nc"
    lcoupling_out = .false.
    iform_input = 1
    ltransport = .FALSE.
    ibc_type = 0
    dn_bcpara = "bias_correction_parameter.nc"
    lbc_write = .FALSE.

    ! read namelist HD_CTL

    inml = open_nml ('namelist.hd')
    iunit = position_nml ('HD_CTL', inml, status=read_status)
    SELECT CASE (read_status)
    CASE (POSITIONED)
       READ (iunit, hd_ctl)
       CALL message('config_hd', 'Namelist HD_CTL: ')
       WRITE(nout, hd_ctl)
    END SELECT

    ! Check Run period and calculate number of time steps if start and end date are provided

    IF (LEN_TRIM(date_start).GT.0 .AND. LEN_TRIM(date_end).GT.0) THEN
      CALL day_difference(date_start, date_end, year1, month1, day1, ndays)
      nstep = ndays * NINT(86400._dp / delta_time)
      WRITE (message_text,'(A,I2.2,A1,I2.2,A1,I4.4)') 'Start date: ', day1, '.', month1, '.', year1
      CALL message('config_hd', message_text)
      WRITE (message_text,*) 'Calculated no. of time steps. ', nstep
      CALL message('config_hd', message_text)
    ELSE
      WRITE (message_text,*) 'No. of time steps taken from namelist (Def.: 365): ', nstep
      CALL message('config_hd', message_text)
    ENDIF

    ! Check coupling information

    CALL set_coupling_type(coupling_type)
    IF (lcoupling_out .AND. coupling_type.NE.2) THEN
       WRITE (message_text,*) 'lcoupling_out is set, but coupling_type != 2! -> Error'
       CALL finish ('config_hd', message_text)
    ENDIF

    IF (coupling_type.EQ.0) THEN
       IF (lcoupling_atm) THEN
          WRITE (message_text,*) 'coupling_type = 0 but lcoupling_atm = True -> Error'
          CALL finish ('config_hd', message_text)
       ENDIF
       IF (lcoupling_oce) THEN
          WRITE (message_text,*) 'coupling_type = 0 but lcoupling_oce = True -> Error'
          CALL finish ('config_hd', message_text)
       ENDIF
    ELSE
       IF (.NOT. lcoupling_atm .AND. .NOT. lcoupling_oce) THEN
          WRITE (message_text,*) 'coupling_type != 0 but lcoupling_atm = lcoupling_oce = False -> Error'
          CALL finish ('config_hd', message_text)
       ENDIF
    ENDIF

    ! Check bias correction information
    IF (ibc_type .EQ. 0) lbc_write = .FALSE.  ! No bias corrction, no bias corrected output

  END SUBROUTINE config_hd

  !------------------------------------------------------------------------------------
  ! Read forcing coordinates and masks 
  !
  ! Array allocations for HD model simulations. In coupled echam/HD runs
  ! these variables are provided by echam.
  !------------------------------------------------------------------------------------
  SUBROUTINE hd_init_forcing

    USE mo_gaussgrid,     ONLY: gauaw, gridarea, philat, philon
    USE mo_constants,     ONLY: a, api
    USE mo_netcdf,        ONLY: file_info, nf_max_name, io_inq_dimid, &
                                io_inq_dimlen, io_inq_varid, io_get_var_double

    TYPE (file_info)  :: fileinfo

    INTEGER dimid, varid, fileid

    CHARACTER(nf_max_name) :: filename
    LOGICAL           :: lex
    REAL(dp), ALLOCATABLE :: zgw(:), zgmu(:)
    INTEGER           :: jgl
    REAL(dp)          :: ra, rb, rd, rh

    ! read land sea masks

!---- Ha Ho-Hagemann {
#ifdef COUP_OAS
    filename = 'masks_HD.nc'
#else
    filename = 'masks.nc'
#endif
!---- Ha Ho-Hagemann }

    INQUIRE (file=filename, exist=lex)
write(nout,*) 'OTBhd_init_forcing in', lex, filename
    IF (.NOT. lex) THEN
      WRITE (message_text,*) 'Could not open file <',TRIM(filename),'>'
      CALL message('read_hd_forcing', message_text)
      CALL finish ('read_hd_forcing', 'run terminated.')
    ENDIF

    ! Read grid info for forcing 
    CALL read_grid_info(filename, grid_forcing, 'lon', 'lat', 1, 1)


    fileinfo%opened = .FALSE.
write(nout,*) 'OTBIO_open ', fileinfo%opened, trim(filename), IO_READ
    CALL IO_open (filename, fileinfo, IO_READ)
write(nout,*) 'OTBIO_open out',filename, fileinfo%opened, IO_READ
    WRITE (message_text,*) 'Reading land sea masks from file ', TRIM(filename)
    CALL message('read_hd_forcing', message_text)

    fileID = fileinfo%file_id

!!    CALL IO_inq_dimid (fileID, 'lon', dimid)
!!    CALL IO_inq_dimlen (fileID, dimid, nlon)
!!    CALL IO_inq_dimid (fileID, 'lat', dimid)
!!    CALL IO_inq_dimlen (fileID, dimid, ngl)
    nlon = grid_forcing%nlon
    ngl = grid_forcing%nlat

    ALLOCATE (philon(nlon))
    CALL IO_inq_varid (fileID, 'lon', varid)
    CALL IO_get_var_double (fileID, varid, philon)

    ALLOCATE (philat(ngl))
    CALL IO_inq_varid (fileID, 'lat', varid)
    CALL IO_get_var_double (fileID, varid, philat)
!hag
    WRITE (message_text,*) 'Forcing dimensions: nlon = ', nlon,  ' ngl = ', ngl
    CALL message('read_hd_forcing', message_text)

    CALL IO_inq_varid (fileID, 'SLM', varid)
    ALLOCATE(slm(nlon,ngl))
    CALL IO_get_var_double (fileID, varid, slm)

    CALL IO_inq_varid (fileID, 'SLF', varid)
    ALLOCATE(slf(nlon,ngl))
    CALL IO_get_var_double (fileID, varid, slf)

    CALL IO_inq_varid (fileID, 'ALAKE', varid)
    ALLOCATE(alake(nlon,ngl))
    CALL IO_get_var_double (fileID, varid, alake)

    CALL IO_inq_varid (fileID, 'GLAC', varid)
    ALLOCATE(glac(nlon,ngl))
    CALL IO_get_var_double (fileID, varid, glac)

    CALL IO_close(fileinfo)

    WRITE (message_text,*) 'Before grid area calc'
    CALL message('read_hd_forcing', message_text)

    ALLOCATE (gridarea(ngl))
    IF (nremap.eq.1) THEN
      
       !-- Calculate grid area using Gaussian weights

       ALLOCATE (zgw(ngl))
       ALLOCATE (zgmu(ngl))
       CALL gauaw(zgmu, zgw, ngl)
       DO jgl = 1, ngl
          gridarea(jgl)  = 0.5_dp * zgw(jgl)/nlon * 4 * api * a**2
       END DO
       DEALLOCATE (zgw, zgmu)
    ELSE

       !-- Calculate the grid area - Works only for regular grids!

       CALL areacalc(grid_forcing, philat, gridarea)

    END IF

    WRITE (message_text,*) 'After grid area calc'
    CALL message('read_hd_forcing', message_text)

    !-- Memory allocation

    ALLOCATE(runoff(nlon,ngl))
    ALLOCATE(drain(nlon,ngl))
    ALLOCATE(disch(nlon,ngl))

    ALLOCATE(awfre(nlon,ngl))
    awfre(:,:) = 0._dp
    ALLOCATE(apmecal(nlon,ngl))
    apmecal(:,:) = 0._dp

    ALLOCATE(dummy_sp(nlon,ngl))
    IF (nremap.EQ.0 .OR. nremap.EQ.3) THEN
      ALLOCATE(dummy_sp_CCLM(grid_hd%nlon,grid_hd%nlat))
    ENDIF

  END SUBROUTINE hd_init_forcing

  !------------------------------------------------------------------------------------
  ! init netcdf dimensions. This is needed with IO_open (hd_init_forcing). The
  ! dimension lenghts given here (HD grid) do not matter. 
  !------------------------------------------------------------------------------------
  SUBROUTINE hd_init_dims
    USE mo_netcdf, ONLY: add_dim, IO_ndim_ids

    INTEGER :: ndum_lon = 720 ! Dummy info - not relevant
    INTEGER :: ndum_lat = 360 ! Dummy info - not relevant

IO_ndim_ids = 0 ! OTB
write(nout,*) 'OTBPR in IO_init_dims', IO_ndim_ids
    CALL add_dim ("lon", ndum_lon, "longitude", "degrees_east" )
    CALL add_dim ("lat", ndum_lat, "latitude", "degrees_north" )
write(nout,*) 'OTBPR in IO_init_dims', IO_ndim_ids

  END SUBROUTINE hd_init_dims

  !------------------------------------------------------------------------------------
  ! open hydrology model in- and output files
  !------------------------------------------------------------------------------------
  SUBROUTINE hd_init_io

    INTEGER  :: ios          ! I/O status

    CHARACTER(nf_max_name) :: filename, varname

    ! open file for outflow data on the forcing data grid

    IF (locean) THEN
       filename = 'hd_discharge.nc'
       varname = 'disch'
       CALL hd_highres_open(filename, varname, discharge_file_id, klon=nlon, klat=ngl)
       WRITE (message_text,*) 'Outflow to ocean on HD model grid written to: ', TRIM(filename)
       CALL message('hd_init_io', message_text)
    END IF

    ! open output file for time averaged outflow data

    filename = 'hd_meanflow.nc'
    varname = 'friv'
    CALL hd_highres_open(filename, varname, meanflow_file_id)
    WRITE (message_text,*) 'Time averaged outflow on HD model grid written to: ', TRIM(filename)
    CALL message('hd_init_io', message_text)

    ! open output file for bias corrected outflow (discharge) data
    IF (lbc_write) THEN
      filename = 'hd_bcflow.nc'
      varname = 'friv_bc'
      CALL hd_highres_open(filename, varname, bcflow_file_id)
      WRITE (message_text,*) 'Bias corrected outflow on HD model grid written to: ', TRIM(filename)
      CALL message('hd_init_io', message_text)
    ENDIF

    ! ***** Open forcing files only if HD model is running offline, i.e neither coupled via OASIS nor in MPI-ESM

    IF (.NOT. is_coupled_run() .OR. .NOT. lcoupling_atm) THEN
      ! open input file with overlandflow (= runoff) data

      IF (iform_input.EQ.0) THEN
        lurun = find_next_free_unit (51,100)
        OPEN(lurun, FILE=runoff_file, FORM='unformatted', STATUS='old', IOSTAT=ios)
        IF (ios /= 0) THEN
          WRITE (message_text,*) 'Error opening file ', runoff_file
          CALL finish ('hd_init_io', message_text)
        ENDIF
      ELSE IF (iform_input.EQ.1) THEN
        runofffile%opened = .FALSE.
        CALL IO_open (runoff_file, runofffile, IO_READ)
        WRITE (message_text,*) 'Reading runoff forcing data from file ', TRIM(runoff_file)
        CALL message('hd_init_io', message_text)
        CALL IO_inq_varid (runofffile%file_id, 'runoff', runoffid)
      ENDIF

      ! open input file with baseflow (= drainage) data

      IF (lbase) THEN
        IF (iform_input.EQ.0) THEN
          lubas = find_next_free_unit (51,100)
          OPEN(lubas, FILE=drainage_file, FORM='unformatted', STATUS='old', IOSTAT=ios)
          IF (ios /= 0) THEN
            WRITE (message_text,*) 'Error opening file ', drainage_file
            CALL finish ('hd_init_io', message_text)
          ENDIF
        ELSE IF (iform_input.EQ.1) THEN
          drainfile%opened = .FALSE.
          CALL IO_open (drainage_file, drainfile, IO_READ)
          WRITE (message_text,*) 'Reading drainage forcing data from file ', TRIM(drainage_file)
          CALL message('hd_init_io', message_text)
          CALL IO_inq_varid (drainfile%file_id, 'drainage', drainid)
        ENDIF
      ENDIF
    ENDIF

    ! open output file for time averaged outflow data on ocean grid for coupling_type = 2
    IF (lcoupling_out) THEN
      CALL hd_on_ocean_open
    ENDIF

  END SUBROUTINE hd_init_io

  !------------------------------------------------------------------------------------
  ! read forcing data for the current time step
  !------------------------------------------------------------------------------------
  SUBROUTINE hd_update_forcing

    INTEGER,  SAVE         :: read_record = 0  ! counter for time steps or days read
    INTEGER,  DIMENSION(3) :: start
    INTEGER,  DIMENSION(3) :: count
    
    IF (iform_input.EQ.0) THEN
      READ(lurun) IHEAD
      READ(lurun) dummy_sp
      runoff(:,:) = REAL(dummy_sp,dp)
      IF (ABS(ufakru-1.) .GT. 1.E-6) runoff = runoff*ufakru
      IF (lbase) THEN
        READ(lubas) IHEAD
        READ(lubas) dummy_sp
        drain(:,:) = REAL(dummy_sp,dp)
        IF (ABS(ufakru-1.) .GT. 1.E-6) drain = drain*ufakru
      ENDIF
    ELSE IF (iform_input.EQ.1) THEN
      read_record = read_record + 1
      start(:) = (/1,1,read_record/)
      count(:) = (/nlon,ngl,1/)
      CALL IO_get_vara_double (runofffile%file_id, runoffid, start(:), count(:), runoff(:,:))
      IF (ABS(ufakru-1.) .GT. 1.E-6) runoff = runoff*ufakru

      IF (lbase) THEN
        CALL IO_get_vara_double (drainfile%file_id, drainid, start(:), count(:), drain(:,:))
        IF (ABS(ufakru-1.) .GT. 1.E-6) drain = drain*ufakru
      ENDIF
    ENDIF

  END SUBROUTINE hd_update_forcing

  !------------------------------------------------------------------------------------
  !---- Ha Ho-Hagemann read directly CCLM on the HD's grid:
  !------------------------------------------------------------------------------------

!---- Ha Ho-Hagemann {
#if !defined(COUP_OAS) && !defined(COUP_YAC)
  SUBROUTINE hd_update_forcing_noremap
    INTEGER,  SAVE         :: read_record = 0  ! counter for time steps or days read
    INTEGER,  DIMENSION(3) :: start
    INTEGER,  DIMENSION(3) :: count
    
    IF (iform_input.EQ.0) THEN
      READ(lurun) IHEAD
      READ(lurun) dummy_sp_CCLM
      runoff_s(:,:) = REAL(dummy_sp_CCLM,dp)
      IF (ABS(ufakru-1.) .GT. 1.E-6) runoff_s=runoff_s*ufakru

      IF (lbase) THEN
         READ(lubas) IHEAD
         READ(lubas) dummy_sp_CCLM
         runoff_dr(:,:) = REAL(dummy_sp_CCLM,dp)
         IF (ABS(ufakru-1.) .GT. 1.E-6) runoff_dr=runoff_dr*ufakru
      ENDIF
    ELSE IF (iform_input.EQ.1) THEN
      read_record = read_record + 1
      start(:) = (/1,1,read_record/)
      count(:) = (/nlon,ngl,1/)
      CALL IO_get_vara_double (runofffile%file_id, runoffid, start(:), count(:), runoff_s(:,:))
      IF (ABS(ufakru-1.) .GT. 1.E-6_dp) runoff_s=runoff_s*ufakru

      IF (lbase) THEN
        CALL IO_get_vara_double (drainfile%file_id, drainid, start(:), count(:), runoff_dr(:,:))
         IF (ABS(ufakru-1.) .GT. 1.E-6) runoff_dr=runoff_dr*ufakru
      ENDIF
    ENDIF
    
  END SUBROUTINE hd_update_forcing_noremap
#endif
!---- Ha Ho-Hagemann }
  
  !------------------------------------------------------------------------------------
  ! write hydrology output of the current time step
  !------------------------------------------------------------------------------------
  SUBROUTINE hd_write_output

    INTEGER, SAVE :: counter_d = 0
    INTEGER, SAVE :: counter_bc = 0

    ! write inflow data for each gridbox.
    ! Note: The time steps are shifted backwards (time step 2 -> 1). The current inflow
    ! is flowing at the end of the current time step, it thus corresponds to the
    ! beginning of the following time step.

    CALL GMITWRI(friv + water_to_ocean, hd_discharge_accu, iout)

    ! write the discharge array (on the grid of the forcing data)

    IF (locean) THEN
       CALL hd_highres_write(discharge_file_id, disch, counter_d, nlon, ngl)
    ENDIF

    ! write bias corrected output 
     IF (lbc_write) THEN
       IF (is_coupled_run() .AND. lcoupling_oce) THEN
         CALL hd_highres_write(bcflow_file_id, hd_outflow, counter_bc, grid_hd%nlon, grid_hd%nlat)
       ELSE
         CALL hd_highres_write(bcflow_file_id, hd_bc_outflow, counter_bc, grid_hd%nlon, grid_hd%nlat)
       ENDIF
     ENDIF

  END SUBROUTINE hd_write_output

END PROGRAM hd_driver
