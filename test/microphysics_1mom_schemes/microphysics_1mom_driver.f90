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
PROGRAM microphysics_1mom_driver
  USE netcdf
  USE ISO_FORTRAN_ENV, ONLY: error_unit, wp => real64, i8 => int64
  USE microphysics_1mom_schemes, ONLY: microphysics_1mom_init, graupel_run, cloudice_run, kessler_run, cloudice2mom_run
  USE mo_exception, ONLY: init_logger, message_text, finish, message

  IMPLICIT NONE
  CHARACTER(120) :: option

  ! Check if a command line argument is provided
  IF (COMMAND_ARGUMENT_COUNT() /= 1) THEN
    WRITE (*, *) "Usage: ./microphysics_1mom_driver <kessler,graupel,cloudice,cloudice2mom>"
    STOP
  ELSE
    ! Get the command line argument
    CALL GET_COMMAND_ARGUMENT(1, option)
  END IF

  ! Convert the option to lowercase for case-insensitivity
  option = TRIM(ADJUSTL(option))

  ! Determine the chosen option
  SELECT CASE (option)
  CASE ('graupel')
    CALL graupel_with_all_options()
  CASE ('cloudice')
    CALL cloudice_with_all_options()
  CASE ('cloudice2mom')
    CALL cloudice2mom_with_all_options()
  CASE ('kessler')
    CALL kessler_with_all_options()
  CASE default
    WRITE (*, *) "Invalid option. Choose between graupel, cloudice, cloudice2mom and kessler."
    error STOP
  END SELECT

CONTAINS

  SUBROUTINE graupel_with_all_options()

    CALL init_logger(0, .TRUE., error_unit, callback_abort=custom_exit)

    CALL message('graupel', 'run for all options')

    ! base config
    CALL microphysics_1mom_driver_interface('graupel_1.nc', igscp=2, ldiag_ttend=.FALSE., &
                                            ldiag_qtend=.FALSE., ldass_lhn=.FALSE., msg_level=5)

    ! activate tendencies
    CALL microphysics_1mom_driver_interface('graupel_2.nc', igscp=2, ldiag_ttend=.TRUE., &
                                            ldiag_qtend=.TRUE., ldass_lhn=.FALSE., msg_level=5)

    ! activate lhn, verbose
    CALL microphysics_1mom_driver_interface('graupel_3.nc', igscp=2, ldiag_ttend=.TRUE., &
                                            ldiag_qtend=.TRUE., ldass_lhn=.FALSE., msg_level=100)

  END SUBROUTINE graupel_with_all_options

  SUBROUTINE cloudice_with_all_options()

    CALL init_logger(0, .TRUE., error_unit, callback_abort=custom_exit)

    CALL message('cloudice', 'run for all options')

    ! base config
    CALL microphysics_1mom_driver_interface('cloudice_1.nc', igscp=1, ldiag_ttend=.FALSE., &
                                            ldiag_qtend=.FALSE., ldass_lhn=.FALSE., msg_level=5)

    ! activate tendencies
    CALL microphysics_1mom_driver_interface('cloudice_2.nc', igscp=1, ldiag_ttend=.TRUE., &
                                            ldiag_qtend=.TRUE., ldass_lhn=.FALSE., msg_level=5)

    ! activate lhn, verbose
    CALL microphysics_1mom_driver_interface('cloudice_3.nc', igscp=1, ldiag_ttend=.TRUE., &
                                            ldiag_qtend=.TRUE., ldass_lhn=.TRUE., msg_level=100)

  END SUBROUTINE cloudice_with_all_options

  SUBROUTINE cloudice2mom_with_all_options()

    CALL init_logger(0, .TRUE., error_unit, callback_abort=custom_exit)

    CALL message('cloudice2mom', 'run for all options')

    ! base config
    CALL microphysics_1mom_driver_interface('cloudice2mom_1.nc', igscp=3, ldiag_ttend=.FALSE., &
                                            ldiag_qtend=.FALSE., ldass_lhn=.FALSE., msg_level=5)

    ! activate tendencies
    CALL microphysics_1mom_driver_interface('cloudice2mom_2.nc', igscp=3, ldiag_ttend=.TRUE., &
                                            ldiag_qtend=.TRUE., ldass_lhn=.FALSE., msg_level=5)

    ! activate lhn, verbose
    CALL microphysics_1mom_driver_interface('cloudice2mom_3.nc', igscp=3, ldiag_ttend=.TRUE., &
                                            ldiag_qtend=.TRUE., ldass_lhn=.TRUE., msg_level=100)

  END SUBROUTINE cloudice2mom_with_all_options

  SUBROUTINE kessler_with_all_options()

    CALL init_logger(0, .TRUE., error_unit, callback_abort=custom_exit)

    CALL message('kessler', 'run for all options')

    ! base config, activate tendencies and lhn, verbose
    CALL microphysics_1mom_driver_interface('kessler_1.nc', igscp=9, ldiag_ttend=.TRUE., &
                                            ldiag_qtend=.TRUE., ldass_lhn=.TRUE., msg_level=100)

  END SUBROUTINE kessler_with_all_options

  SUBROUTINE microphysics_1mom_driver_interface(out_file, igscp, ldiag_ttend, ldiag_qtend, ldass_lhn, msg_level)

    CHARACTER(LEN=*), INTENT(in) :: out_file
    LOGICAL, INTENT(in) :: ldiag_ttend, ldiag_qtend, ldass_lhn
    INTEGER, INTENT(in) :: igscp, msg_level

    CHARACTER(LEN=:), ALLOCATABLE :: input_file
    INTEGER :: itime
    INTEGER, PARAMETER :: ithermo_water = 1
    REAL(wp) :: dt, qnc

    ! Parameters from the input file
    INTEGER :: ncells, nlev
    REAL(wp), DIMENSION(:, :), ALLOCATABLE :: z, t, p, rho, qv, qc, qi, qr, qs, qg, w

    ! Precalculated parameters
    REAL(wp), DIMENSION(:, :), ALLOCATABLE :: dz

    ! Extra fields required to call graupel
    REAL(wp), DIMENSION(:, :), ALLOCATABLE :: pflx, zninc
    REAL(wp), DIMENSION(:), ALLOCATABLE :: prr_gsp, pri_gsp, prs_gsp, prg_gsp
    REAL(wp), DIMENSION(:), ALLOCATABLE :: qnc_s
    REAL(wp), DIMENSION(:, :), ALLOCATABLE :: ddt_tend_t, ddt_tend_qv, ddt_tend_qc
    REAL(wp), DIMENSION(:, :), ALLOCATABLE :: ddt_tend_qi, ddt_tend_qr, ddt_tend_qs, qni, ninact

    ! Timer variables
    INTEGER(i8) :: start_count, end_count, count_rate

    input_file = "input-data.nc"
    itime = 230
    dt = 30.0
    qnc = 100.0

    CALL read_fields(input_file, itime, ncells, nlev, z, t, p, rho, qv, qc, qi, qr, qs, qg, w)

    dz = calc_dz(z, ncells, nlev)

    ! No output of fields below implented yet
    ALLOCATE (prr_gsp(ncells)); prr_gsp(:) = 0
    ALLOCATE (pri_gsp(ncells)); pri_gsp(:) = 0
    ALLOCATE (prs_gsp(ncells)); prs_gsp(:) = 0
    ALLOCATE (prg_gsp(ncells)); prg_gsp(:) = 0
    ALLOCATE (pflx(ncells, nlev)); pflx(:, :) = 0
    ALLOCATE (ddt_tend_t(ncells, nlev)); ddt_tend_t(:, :) = 0
    ALLOCATE (ddt_tend_qv(ncells, nlev)); ddt_tend_qv(:, :) = 0
    ALLOCATE (ddt_tend_qc(ncells, nlev)); ddt_tend_qc(:, :) = 0
    ALLOCATE (ddt_tend_qi(ncells, nlev)); ddt_tend_qi(:, :) = 0
    ALLOCATE (ddt_tend_qr(ncells, nlev)); ddt_tend_qr(:, :) = 0
    ALLOCATE (ddt_tend_qs(ncells, nlev)); ddt_tend_qs(:, :) = 0
    ALLOCATE (qni(ncells, nlev)); qni(:, :) = qnc
    ALLOCATE (ninact(ncells, nlev)); ninact(:, :) = 100.0
    ALLOCATE (zninc(ncells, nlev)); zninc(:, :) = qnc

    ALLOCATE (qnc_s(ncells)); qnc_s(:) = qnc

    !$ACC DATA &
    !$ACC   COPY(dz, t, p, rho, qv, qc, qi, qr, qs, qg, qnc_s, w, qni, ninact, zninc) &
    !$ACC   COPY(prr_gsp, prs_gsp, pri_gsp, prg_gsp, pflx) &
    !$ACC   COPY(ddt_tend_t, ddt_tend_qv, ddt_tend_qc, ddt_tend_qi, ddt_tend_qr, ddt_tend_qs)

    CALL microphysics_1mom_init( &
      igscp=igscp, &
      tune_zceff_min=0.01_wp, &
      tune_v0snow=-1.0_wp, &
      tune_zvz0i=1.25_wp, &
      tune_mu_rain=0.0_wp, &
      tune_icesedi_exp=0.30_wp, &
      tune_rain_n0_factor=1.0_wp)

    CALL SYSTEM_CLOCK(start_count, count_rate)

    IF (igscp == 2) THEN
      CALL graupel_run(                                     &
          & nvec=ncells, & !> in:  actual array size
          & ke=nlev, & !< in:  actual array size
          & ivstart=1, & !< in:  start index of calculation
          & ivend=ncells, & !< in:  end index of calculation
          & kstart=1, & !< in:  vertical start index
          & zdt=dt, & !< in:  timestep
          & qi0=0.00_wp,    &
          & qc0=0.00_wp,    &
          & dz=dz, & !< in:  vertical layer thickness
          & t=t, & !< in:  temp,tracer,...
          & p=p, & !< in:  full level pres
          & rho=rho, & !< in:  density
          & qv=qv, & !< in:  spec. humidity
          & qc=qc, & !< in:  cloud water
          & qi=qi, & !< in:  cloud ice
          & qr=qr, & !< in:  rain water
          & qs=qs, & !< in:  snow
          & qg=qg, & !< in:  graupel
          & qnc=qnc_s, & !< cloud number concentration
          & zninc=zninc, & !< number of cloud ice crystals at nucleation
          & prr_gsp=prr_gsp, & !< out: precipitation rate of rain
          & prs_gsp=prs_gsp, & !< out: precipitation rate of snow
          & pri_gsp=pri_gsp, & !< out: precipitation rate of cloud ice
          & prg_gsp=prg_gsp, & !< out: precipitation rate of graupel
          & qrsflux=pflx, & !< out: precipitation flux
          & ldiag_ttend=ldiag_ttend, & !< in:  if temp. tendency shall be diagnosed
          & ldiag_qtend=ldiag_qtend, & !< in:  if moisture tendencies shall be diagnosed
          & ddt_tend_t=ddt_tend_t, & !< out: tendency temperature
          & ddt_tend_qv=ddt_tend_qv, & !< out: tendency QV
          & ddt_tend_qc=ddt_tend_qc, & !< out: tendency QC
          & ddt_tend_qi=ddt_tend_qi, & !< out: tendency QI
          & ddt_tend_qr=ddt_tend_qr, & !< out: tendency QR
          & ddt_tend_qs=ddt_tend_qs, & !< out: tendency QS
          & idbg=msg_level, &
          & l_cv=.TRUE., &
          & ldass_lhn=ldass_lhn, &
          & ithermo_water=ithermo_water) !< in: latent heat choice
    ELSEIF (igscp == 1) THEN
      CALL cloudice_run(                                     &
          & nvec=ncells, & !> in:  actual array size
          & ke=nlev, & !< in:  actual array size
          & ivstart=1, & !< in:  start index of calculation
          & ivend=ncells, & !< in:  end index of calculation
          & kstart=1, & !< in:  vertical start index
          & zdt=dt, & !< in:  timestep
          & qi0=0.00_wp,    &
          & qc0=0.00_wp,    &
          & dz=dz, & !< in:  vertical layer thickness
          & t=t, & !< in:  temp,tracer,...
          & p=p, & !< in:  full level pres
          & rho=rho, & !< in:  density
          & qv=qv, & !< in:  spec. humidity
          & qc=qc, & !< in:  cloud water
          & qi=qi, & !< in:  cloud ice
          & qr=qr, & !< in:  rain water
          & qs=qs, & !< in:  rain water
          & qnc=qnc_s, & !< cloud number concentration
          & zninc=zninc, & !< number of cloud ice crystals at nucleation
          & prr_gsp=prr_gsp, & !< out: precipitation rate of rain
          & prs_gsp=prs_gsp, & !< out: precipitation rate of snow
          & pri_gsp=pri_gsp, & !< out: precipitation rate of cloud ice
          & qrsflux=pflx, & !< out: precipitation flux
          & ldiag_ttend=ldiag_ttend, & !< in:  if temp. tendency shall be diagnosed
          & ldiag_qtend=ldiag_qtend, & !< in:  if moisture tendencies shall be diagnosed
          & ddt_tend_t=ddt_tend_t, & !< out: tendency temperature
          & ddt_tend_qv=ddt_tend_qv, & !< out: tendency QV
          & ddt_tend_qc=ddt_tend_qc, & !< out: tendency QC
          & ddt_tend_qi=ddt_tend_qi, & !< out: tendency QI
          & ddt_tend_qr=ddt_tend_qr, & !< out: tendency QR
          & ddt_tend_qs=ddt_tend_qs, & !< out: tendency QS
          & idbg=msg_level, &
          & l_cv=.TRUE., &
          & ldass_lhn=ldass_lhn, &
          & ithermo_water=ithermo_water) !< in: latent heat choice

    ELSEIF (igscp == 3) THEN
      CALL cloudice2mom_run(                                     &
          & nvec=ncells, & !> in:  actual array size
          & ke=nlev, & !< in:  actual array size
          & ivstart=1, & !< in:  start index of calculation
          & ivend=ncells, & !< in:  end index of calculation
          & kstart=1, & !< in:  vertical start index
          & zdt=dt, & !< in:  timestep
          & qi0=0.00_wp,    &
          & qc0=0.00_wp,    &
          & dz=dz, & !< in:  vertical layer thickness
          & t=t, & !< in:  temp,tracer,...
          & p=p, & !< in:  full level pres
          & w=w, &
          & rho=rho, & !< in:  density
          & qv=qv, & !< in:  spec. humidity
          & qc=qc, & !< in:  cloud water
          & qi=qi, & !< in:  cloud ice
          & qr=qr, & !< in:  rain water
          & qs=qs, & !< in:  rain water
          & qnc=qnc_s, & !< cloud number concentration
          & qni=qni, &
          & ninact=ninact, &
          & prr_gsp=prr_gsp, & !< out: precipitation rate of rain
          & prs_gsp=prs_gsp, & !< out: precipitation rate of snow
          & pri_gsp=pri_gsp, & !< out: precipitation rate of cloud ice
          & qrsflux=pflx, & !< out: precipitation flux
          & ldiag_ttend=ldiag_ttend, & !< in:  if temp. tendency shall be diagnosed
          & ldiag_qtend=ldiag_qtend, & !< in:  if moisture tendencies shall be diagnosed
          & ddt_tend_t=ddt_tend_t, & !< out: tendency temperature
          & ddt_tend_qv=ddt_tend_qv, & !< out: tendency QV
          & ddt_tend_qc=ddt_tend_qc, & !< out: tendency QC
          & ddt_tend_qi=ddt_tend_qi, & !< out: tendency QI
          & ddt_tend_qr=ddt_tend_qr, & !< out: tendency QR
          & ddt_tend_qs=ddt_tend_qs, & !< out: tendency QS
          & idbg=msg_level, &
          & l_cv=.TRUE., &
          & ldass_lhn=ldass_lhn, &
          & ithermo_water=ithermo_water) !< in: latent heat choice

    ELSEIF (igscp == 9) THEN
      CALL kessler_run(                                     &
          & nvec=ncells, & !> in:  actual array size
          & ke=nlev, & !< in:  actual array size
          & ivstart=1, & !< in:  start index of calculation
          & ivend=ncells, & !< in:  end index of calculation
          & kstart=1, & !< in:  vertical start index
          & zdt=dt, & !< in:  timestep
          & qc0=0.00_wp,    &
          & dz=dz, & !< in:  vertical layer thickness
          & t=t, & !< in:  temp,tracer,...
          & p=p, & !< in:  full level pres
          & rho=rho, & !< in:  density
          & qv=qv, & !< in:  spec. humidity
          & qc=qc, & !< in:  cloud water
          & qr=qr, & !< in:  rain water
          & prr_gsp=prr_gsp, & !< out: precipitation rate of rain
          & qrsflux=pflx, & !< out: precipitation flux
          & ldiag_ttend=ldiag_ttend, & !< in:  if temp. tendency shall be diagnosed
          & ldiag_qtend=ldiag_qtend, & !< in:  if moisture tendencies shall be diagnosed
          & ddt_tend_t=ddt_tend_t, & !< out: tendency temperature
          & ddt_tend_qv=ddt_tend_qv, & !< out: tendency QV
          & ddt_tend_qc=ddt_tend_qc, & !< out: tendency QC
          & ddt_tend_qr=ddt_tend_qr, & !< out: tendency QR
          & idbg=msg_level, &
          & l_cv=.TRUE., &
          & ldass_lhn=ldass_lhn)
    END IF

    CALL SYSTEM_CLOCK(end_count)

    !$ACC END DATA

    WRITE (message_text, "(a,f0.2,a)") "Wall clock time: ", 1000.0*(end_count - start_count)/count_rate, " ms"
    CALL message('microphysics_1mom_driver', message_text)

    CALL write_fields(out_file, igscp, ncells, nlev, t, qv, qc, qi, qr, qs, qg, w, qni, ninact, prr_gsp, &
                      prs_gsp, pri_gsp, prg_gsp, pflx, ddt_tend_t, ddt_tend_qv, ddt_tend_qc, ddt_tend_qi, &
                      ddt_tend_qr, ddt_tend_qs, ldiag_ttend, ldiag_qtend)

    DEALLOCATE (prr_gsp)
    DEALLOCATE (pri_gsp)
    DEALLOCATE (prs_gsp)
    DEALLOCATE (prg_gsp)
    DEALLOCATE (pflx)
    DEALLOCATE (ddt_tend_t)
    DEALLOCATE (ddt_tend_qv)
    DEALLOCATE (ddt_tend_qc)
    DEALLOCATE (ddt_tend_qi)
    DEALLOCATE (ddt_tend_qr)
    DEALLOCATE (ddt_tend_qs)
    DEALLOCATE (qni)
    DEALLOCATE (ninact)
    DEALLOCATE (zninc)

    DEALLOCATE (qnc_s)

  END SUBROUTINE microphysics_1mom_driver_interface

  FUNCTION calc_dz(z, ncells, nlev)
    REAL(wp), DIMENSION(:, :), INTENT(IN) :: z
    INTEGER, INTENT(IN) :: ncells, nlev
    REAL(wp), DIMENSION(ncells, nlev) :: calc_dz

    REAL(wp), DIMENSION(ncells, nlev + 1) :: zh
    INTEGER :: k

    zh(:, nlev + 1) = (3.*z(:, nlev) - z(:, nlev - 1))*0.5
    DO k = nlev, 1, -1
      zh(:, k) = 2.0*z(:, k) - zh(:, k + 1)
      calc_dz(:, k) = -zh(:, k + 1) + zh(:, k)
    END DO
  END FUNCTION calc_dz

  SUBROUTINE read_fields(filename, itime, ncells, nlev, z, t, p, rho, qv, qc, qi, qr, qs, qg, w)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: itime
    INTEGER, INTENT(OUT) :: ncells, nlev
    REAL(wp), DIMENSION(:, :), ALLOCATABLE :: z, t, p, rho, qv, qc, qi, qr, qs, qg, w

    INTEGER :: ncid

    CALL handle_nc_code(nf90_open(filename, NF90_NOWRITE, ncid), lineno=__LINE__)
    nlev = get_nc_dim_size(ncid, "height")
    ncells = get_nc_dim_size(ncid, "ncells")
    z = read_nc_2d_field(ncid, "zg", ncells, nlev)
    t = read_nc_2d_field(ncid, "ta", ncells, nlev, start=[1, 1, itime])
    p = read_nc_2d_field(ncid, "pfull", ncells, nlev, start=[1, 1, itime])
    rho = read_nc_2d_field(ncid, "rho", ncells, nlev, start=[1, 1, itime])
    qv = read_nc_2d_field(ncid, "hus", ncells, nlev, start=[1, 1, itime])
    qc = read_nc_2d_field(ncid, "clw", ncells, nlev, start=[1, 1, itime])
    qi = read_nc_2d_field(ncid, "cli", ncells, nlev, start=[1, 1, itime])
    qr = read_nc_2d_field(ncid, "qr", ncells, nlev, start=[1, 1, itime])
    qs = read_nc_2d_field(ncid, "qs", ncells, nlev, start=[1, 1, itime])
    qg = read_nc_2d_field(ncid, "qg", ncells, nlev, start=[1, 1, itime])
    w = read_nc_2d_field(ncid, "wa", ncells, nlev, start=[1, 1, itime])
    CALL handle_nc_code(nf90_close(ncid), lineno=__LINE__)
  END SUBROUTINE read_fields

  SUBROUTINE write_fields(filename, igscp, ncells, nlev, t, qv, qc, qi, qr, qs, qg, w, qni, ninact, &
                          prr_gsp, prs_gsp, pri_gsp, prg_gsp, qrsflux, ddt_tend_t, ddt_tend_qv, ddt_tend_qc, &
                          ddt_tend_qi, ddt_tend_qr, ddt_tend_qs, ldiag_ttend, ldiag_qtend)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    LOGICAL, INTENT(IN) ::ldiag_ttend, ldiag_qtend
    INTEGER, INTENT(IN) :: igscp, ncells, nlev
    REAL(wp), DIMENSION(ncells, nlev), INTENT(in) :: t, qv, qc, qi, qr, qs, qg, w, qni, ninact, qrsflux, ddt_tend_t, &
                                                     ddt_tend_qv, ddt_tend_qc, ddt_tend_qi, ddt_tend_qr, ddt_tend_qs
    REAL(wp), DIMENSION(ncells), INTENT(in) :: prr_gsp, prs_gsp, pri_gsp, prg_gsp

    INTEGER :: ncid, ncells_dimid, nlev_dimid, dimids(2)
    INTEGER :: t_varid, qv_varid, qc_varid, qi_varid, qr_varid, qs_varid, qg_varid, prr_varid, &
               prs_varid, pri_varid, prg_varid, qrsflux_varid, wa_varid, qni_varid, ninact_varid
    INTEGER :: ddt_t_varid, ddt_qv_varid, ddt_qc_varid, ddt_qi_varid, ddt_qr_varid, ddt_qs_varid

    CALL handle_nc_code(nf90_create(TRIM(filename), NF90_NETCDF4, ncid), lineno=__LINE__)

    CALL handle_nc_code(nf90_def_dim(ncid, "ncells", ncells, ncells_dimid), lineno=__LINE__)
    CALL handle_nc_code(nf90_def_dim(ncid, "height", nlev, nlev_dimid), lineno=__LINE__)
    dimids = [ncells_dimid, nlev_dimid]

    CALL handle_nc_code(nf90_def_var(ncid, "ta", NF90_DOUBLE, dimids, t_varid), lineno=__LINE__)
    CALL handle_nc_code(nf90_def_var(ncid, "hus", NF90_DOUBLE, dimids, qv_varid), lineno=__LINE__)
    CALL handle_nc_code(nf90_def_var(ncid, "clw", NF90_DOUBLE, dimids, qc_varid), lineno=__LINE__)
    CALL handle_nc_code(nf90_def_var(ncid, "qr", NF90_DOUBLE, dimids, qr_varid), lineno=__LINE__)
    IF (igscp /= 9) THEN
      CALL handle_nc_code(nf90_def_var(ncid, "qs", NF90_DOUBLE, dimids, qs_varid), lineno=__LINE__)
      CALL handle_nc_code(nf90_def_var(ncid, "cli", NF90_DOUBLE, dimids, qi_varid), lineno=__LINE__)
    END IF
    IF (igscp == 2) THEN
      CALL handle_nc_code(nf90_def_var(ncid, "qg", NF90_DOUBLE, dimids, qg_varid), lineno=__LINE__)
    END IF

    CALL handle_nc_code(nf90_def_var(ncid, "prr", NF90_DOUBLE, ncells_dimid, prr_varid), lineno=__LINE__)
    IF (igscp /= 9) THEN
      CALL handle_nc_code(nf90_def_var(ncid, "prs", NF90_DOUBLE, ncells_dimid, prs_varid), lineno=__LINE__)
      CALL handle_nc_code(nf90_def_var(ncid, "pri", NF90_DOUBLE, ncells_dimid, pri_varid), lineno=__LINE__)
    END IF

    IF (igscp == 2) THEN
      CALL handle_nc_code(nf90_def_var(ncid, "prg", NF90_DOUBLE, ncells_dimid, prg_varid), lineno=__LINE__)
    END IF
    CALL handle_nc_code(nf90_def_var(ncid, "qrsflux", NF90_DOUBLE, dimids, qrsflux_varid), lineno=__LINE__)

    IF (igscp == 3) THEN
      CALL handle_nc_code(nf90_def_var(ncid, "wa", NF90_DOUBLE, dimids, wa_varid), lineno=__LINE__)
      CALL handle_nc_code(nf90_def_var(ncid, "qni", NF90_DOUBLE, dimids, qni_varid), lineno=__LINE__)
      CALL handle_nc_code(nf90_def_var(ncid, "ninact", NF90_DOUBLE, dimids, ninact_varid), lineno=__LINE__)
    END IF

    IF (ldiag_ttend) THEN
      CALL handle_nc_code(nf90_def_var(ncid, "ddt_ta", NF90_DOUBLE, dimids, ddt_t_varid), lineno=__LINE__)
    END IF

    IF (ldiag_qtend) THEN
      CALL handle_nc_code(nf90_def_var(ncid, "ddt_hus", NF90_DOUBLE, dimids, ddt_qv_varid), lineno=__LINE__)
      CALL handle_nc_code(nf90_def_var(ncid, "ddt_clw", NF90_DOUBLE, dimids, ddt_qc_varid), lineno=__LINE__)
      CALL handle_nc_code(nf90_def_var(ncid, "ddt_qr", NF90_DOUBLE, dimids, ddt_qr_varid), lineno=__LINE__)
      IF (igscp /= 9) THEN
        CALL handle_nc_code(nf90_def_var(ncid, "ddt_qs", NF90_DOUBLE, dimids, ddt_qs_varid), lineno=__LINE__)
        CALL handle_nc_code(nf90_def_var(ncid, "ddt_cli", NF90_DOUBLE, dimids, ddt_qi_varid), lineno=__LINE__)
      END IF
    END IF

    CALL handle_nc_code(nf90_enddef(ncid), lineno=__LINE__)

    CALL handle_nc_code(nf90_put_var(ncid, t_varid, t), lineno=__LINE__)
    CALL handle_nc_code(nf90_put_var(ncid, qv_varid, qv), lineno=__LINE__)
    CALL handle_nc_code(nf90_put_var(ncid, qc_varid, qc), lineno=__LINE__)
    CALL handle_nc_code(nf90_put_var(ncid, qr_varid, qr), lineno=__LINE__)
    IF (igscp /= 9) THEN
      CALL handle_nc_code(nf90_put_var(ncid, qi_varid, qi), lineno=__LINE__)
      CALL handle_nc_code(nf90_put_var(ncid, qs_varid, qs), lineno=__LINE__)
    END IF
    IF (igscp == 2) THEN
      CALL handle_nc_code(nf90_put_var(ncid, qg_varid, qg), lineno=__LINE__)
    END IF

    CALL handle_nc_code(nf90_put_var(ncid, prr_varid, prr_gsp), lineno=__LINE__)
    IF (igscp /= 9) THEN
      CALL handle_nc_code(nf90_put_var(ncid, prs_varid, prs_gsp), lineno=__LINE__)
      CALL handle_nc_code(nf90_put_var(ncid, pri_varid, pri_gsp), lineno=__LINE__)
    END IF
    IF (igscp == 2) THEN
      CALL handle_nc_code(nf90_put_var(ncid, prg_varid, prg_gsp), lineno=__LINE__)
    END IF
    CALL handle_nc_code(nf90_put_var(ncid, qrsflux_varid, qrsflux), lineno=__LINE__)

    IF (igscp == 3) THEN
      CALL handle_nc_code(nf90_put_var(ncid, wa_varid, w), lineno=__LINE__)
      CALL handle_nc_code(nf90_put_var(ncid, qni_varid, qni), lineno=__LINE__)
      CALL handle_nc_code(nf90_put_var(ncid, ninact_varid, ninact), lineno=__LINE__)
    END IF

    IF (ldiag_ttend) THEN
      CALL handle_nc_code(nf90_put_var(ncid, ddt_t_varid, ddt_tend_t), lineno=__LINE__)
    END IF

    IF (ldiag_qtend) THEN
      CALL handle_nc_code(nf90_put_var(ncid, ddt_qv_varid, ddt_tend_qv), lineno=__LINE__)
      CALL handle_nc_code(nf90_put_var(ncid, ddt_qc_varid, ddt_tend_qc), lineno=__LINE__)
      IF (igscp /= 9) THEN
        CALL handle_nc_code(nf90_put_var(ncid, ddt_qi_varid, ddt_tend_qi), lineno=__LINE__)
        CALL handle_nc_code(nf90_put_var(ncid, ddt_qr_varid, ddt_tend_qr), lineno=__LINE__)
        CALL handle_nc_code(nf90_put_var(ncid, ddt_qs_varid, ddt_tend_qs), lineno=__LINE__)
      END IF
    END IF

    CALL handle_nc_code(nf90_close(ncid), lineno=__LINE__)
  END SUBROUTINE write_fields

  FUNCTION get_nc_dim_size(ncid, dimname)
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: dimname
    INTEGER :: get_nc_dim_size

    INTEGER :: dimid

    CALL handle_nc_code(nf90_inq_dimid(ncid, TRIM(dimname), dimid), lineno=__LINE__)
    CALL handle_nc_code(nf90_inquire_dimension(ncid, dimid, len=get_nc_dim_size), lineno=__LINE__)
  END FUNCTION get_nc_dim_size

  FUNCTION read_nc_2d_field(ncid, varname, nx, ny, start)
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER, INTENT(IN) :: nx, ny
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start
    REAL(wp), DIMENSION(nx, ny) :: read_nc_2d_field

    INTEGER :: varid

    CALL handle_nc_code(nf90_inq_varid(ncid, TRIM(varname), varid), lineno=__LINE__)
    CALL handle_nc_code(nf90_get_var(ncid, varid, read_nc_2d_field, start=start), lineno=__LINE__)
  END FUNCTION read_nc_2d_field

  SUBROUTINE handle_nc_code(code, filename, lineno)
    INTEGER, INTENT(IN) :: code
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: filename
    INTEGER, OPTIONAL, INTENT(IN) :: lineno

    IF (code /= NF90_NOERR) &
      CALL stop_on_error(nf90_strerror(code), filename, lineno)
  END SUBROUTINE handle_nc_code

  SUBROUTINE stop_on_error(msg, filename, lineno)
    CHARACTER(LEN=*), INTENT(IN) :: msg
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: filename
    INTEGER, OPTIONAL, INTENT(IN) :: lineno

    CHARACTER(LEN=:), ALLOCATABLE :: fn
    CHARACTER(LEN=100) :: ln

    fn = ""; ln = ""

    IF (PRESENT(filename)) THEN
      fn = TRIM(filename)//":"
    ELSE
      fn = "microphysics_1mom_driver.f90"//":"
    END IF
    IF (PRESENT(lineno)) WRITE (ln, "(i0,a)") lineno, ":"

    IF (LEN_TRIM(msg) > 0) &
      WRITE (message_text, "(a)") fn//TRIM(ln)//" "//TRIM(msg)
    CALL finish('microphysics_1mom_driver', message_text)
  END SUBROUTINE stop_on_error

  SUBROUTINE custom_exit()
    ERROR STOP
  END SUBROUTINE custom_exit
END PROGRAM microphysics_1mom_driver
