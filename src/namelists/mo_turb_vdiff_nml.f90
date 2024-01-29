! Namelist for VDIFF
!
! These subroutines are called by read_atmo_namelists and do the VDIFF setup.
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

MODULE mo_turb_vdiff_nml
  USE mo_kind, ONLY: wp
  USE mo_exception, ONLY: finish

  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_util_string,         ONLY: tolower
  USE mo_turb_vdiff_config,   ONLY: &
    t_vdiff_config, vdiff_config, vdiff_config_init, vdiff_config_update, vdiff_config_check
  USE mo_turb_vdiff_params,   ONLY: VDIFF_TURB_3DSMAGORINSKY, VDIFF_TURB_TTE

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER :: turb_map(2) = (/ 'tte   ', '3dsmag' /)
  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_turb_vdiff_nml'

  !------------------
  ! Module interface
  !------------------
  PUBLIC :: vdiff_read_namelist

  INTERFACE fill_default
    MODULE PROCEDURE :: fill_default_r
    MODULE PROCEDURE :: fill_default_str
  END INTERFACE


CONTAINS

  SUBROUTINE vdiff_read_namelist( filename )

    CHARACTER(len=*), INTENT(IN) :: filename

    INTEGER :: istat
    INTEGER :: iunit
    INTEGER :: jg

    TYPE(t_vdiff_config) :: defaults

    ! Namelist members

    LOGICAL  :: lsfc_mom_flux(SIZE(vdiff_config))   !< switch on/off surface momentum flux
    LOGICAL  :: lsfc_heat_flux(SIZE(vdiff_config))  !< switch on/off surface heat flux (sensible AND latent)

    REAL(wp) :: pr0(SIZE(vdiff_config))             !< neutral limit Prandtl number, can be varied from about 0.6 to 1.0
    REAL(wp) :: f_tau0(SIZE(vdiff_config))          !< neutral non-dimensional stress factor

    REAL(wp) :: f_tau_limit_fraction(SIZE(vdiff_config))   !< Limit of f_tau/f_tau0 in the stable limit (Mauritsen: 0.25).
    REAL(wp) :: f_theta_limit_fraction(SIZE(vdiff_config)) !< Limit of f_theta/f_theta0 in the stable limit (Mauritsen: 0.).
    REAL(wp) :: f_tau_decay(SIZE(vdiff_config))            !< Decay rate of f_tau/f_tau0 towards its limit (Mauritsen: 4.).
    REAL(wp) :: f_theta_decay(SIZE(vdiff_config))          !< Decay rate of f_theta/f_theta0 towards its limit (Mauritsen: 4.).
    REAL(wp) :: ek_ep_ratio_stable(SIZE(vdiff_config))     !< Ek/Ep ratio in the stable limit (Mauritsen: 1/(0.3+-0.1) - 1).
    REAL(wp) :: ek_ep_ratio_unstable(SIZE(vdiff_config))   !< Ek/Ep ratio in the unstable limit (Mauritsen: 1).

    REAL(wp) :: c_f(SIZE(vdiff_config))             !< mixing length: coriolis term tuning parameter
    REAL(wp) :: c_n(SIZE(vdiff_config))             !< mixing length: stability term tuning parameter
    REAL(wp) :: wmc(SIZE(vdiff_config))             !< ratio of typical horizontal velocity to wstar at free convection
    REAL(wp) :: fsl(SIZE(vdiff_config))             !< fraction of first-level height at which surface fluxes
    !                                                  are nominally evaluated, tuning param for sfc stress
    REAL(wp) :: fbl(SIZE(vdiff_config))             !< 1/fbl: fraction of BL height at which lmix hat its max
    REAL(wp) :: lmix_max(SIZE(vdiff_config))        !< Maximum mixing length.
    REAL(wp) :: z0m_min(SIZE(vdiff_config))         !< Minimum roughness length for momentum [m].
    REAL(wp) :: z0m_ice(SIZE(vdiff_config))         !< Roughness length for momentum over ice [m].
    REAL(wp) :: z0m_oce(SIZE(vdiff_config))         !< Roughness length for momentum over ocean [m].

    CHARACTER(len=15) :: turb(SIZE(vdiff_config))   !< turbulence scheme: 'tte' or '3dsmag'.
    REAL(wp) :: smag_constant(SIZE(vdiff_config))
    REAL(wp) :: turb_prandtl(SIZE(vdiff_config))    !< Turbulent Prandtl number
    REAL(wp) :: km_min(SIZE(vdiff_config))          !< min mass weighted turbulent viscosity
    REAL(wp) :: max_turb_scale(SIZE(vdiff_config))  !< max turbulence length scale
    REAL(wp) :: min_sfc_wind(SIZE(vdiff_config))    !< min sfc wind in free convection limit

    NAMELIST /turb_vdiff_nml/ &
      lsfc_mom_flux, &
      lsfc_heat_flux, &
      pr0, &
      f_tau0, &
      f_tau_limit_fraction, &
      f_theta_limit_fraction, &
      f_tau_decay, &
      f_theta_decay, &
      ek_ep_ratio_stable, &
      ek_ep_ratio_unstable, &
      c_f, &
      c_n, &
      wmc, &
      fsl, &
      fbl, &
      lmix_max, &
      z0m_min, &
      z0m_ice, &
      z0m_oce, &
      turb, &
      smag_constant, &
      turb_prandtl, &
      km_min, &
      max_turb_scale, &
      min_sfc_wind

    !------

    CALL vdiff_config_init(defaults)

    lsfc_mom_flux(:) = defaults%lsfc_mom_flux
    lsfc_heat_flux(:) = defaults%lsfc_heat_flux
    pr0(:) = defaults%pr0
    f_tau0(:) = defaults%f_tau0
    f_tau_limit_fraction(:) = defaults%f_tau_limit_fraction
    f_theta_limit_fraction(:) = defaults%f_theta_limit_fraction
    f_tau_decay(:) = defaults%f_tau_decay
    f_theta_decay(:) = defaults%f_theta_decay
    ek_ep_ratio_stable(:) = defaults%ek_ep_ratio_stable
    ek_ep_ratio_unstable(:) = defaults%ek_ep_ratio_unstable
    c_f(:) = defaults%c_f
    c_n(:) = defaults%c_n
    wmc(:) = defaults%wmc
    fsl(:) = defaults%fsl
    fbl(:) = defaults%fbl
    lmix_max(:) = defaults%lmix_max
    z0m_min(:) = defaults%z0m_min
    z0m_ice(:) = defaults%z0m_ice
    z0m_oce(:) = defaults%z0m_oce
    turb(:) = turb_map(defaults%turb)
    smag_constant(:) = defaults%smag_constant
    turb_prandtl(:) = defaults%turb_prandtl
    km_min(:) = defaults%km_min
    max_turb_scale(:) = defaults%max_turb_scale
    min_sfc_wind(:) = defaults%min_sfc_wind

    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, NML=turb_vdiff_nml)   ! write defaults to temporary text file
    END IF

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      iunit = open_and_restore_namelist('turb_vdiff_nml')
      READ(iunit, NML=turb_vdiff_nml)
      CALL close_tmpfile(iunit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml('turb_vdiff_nml', STATUS=istat)

    IF (POSITIONED == istat) THEN
      pr0(:) = -1._wp
      f_tau0(:) = -1._wp
      f_tau_limit_fraction(:) = -1._wp
      f_theta_limit_fraction(:) = -1._wp
      f_tau_decay(:) = -1._wp
      f_theta_decay(:) = -1._wp
      ek_ep_ratio_stable(:) = -1._wp
      ek_ep_ratio_unstable(:) = -1._wp
      c_f(:) = -1._wp
      c_n(:) = -1._wp
      wmc(:) = -1._wp
      fsl(:) = -1._wp
      fbl(:) = -1._wp
      lmix_max(:) = -1._wp
      z0m_min(:) = -1._wp
      z0m_ice(:) = -1._wp
      z0m_oce(:) = -1._wp
      turb(:) = ''
      smag_constant(:) = -1._wp
      turb_prandtl(:) = -1._wp
      km_min(:) = -1._wp
      max_turb_scale(:) = -1._wp
      min_sfc_wind(:) = -1._wp

      READ(nnml, NML=turb_vdiff_nml)

      CALL fill_default(pr0(:), defaults%pr0)
      CALL fill_default(f_tau0(:), defaults%f_tau0)
      CALL fill_default(f_tau_limit_fraction(:), defaults%f_tau_limit_fraction)
      CALL fill_default(f_theta_limit_fraction(:), defaults%f_theta_limit_fraction)
      CALL fill_default(f_tau_decay(:), defaults%f_tau_decay)
      CALL fill_default(f_theta_decay(:), defaults%f_theta_decay)
      CALL fill_default(ek_ep_ratio_stable(:), defaults%ek_ep_ratio_stable)
      CALL fill_default(ek_ep_ratio_unstable(:), defaults%ek_ep_ratio_unstable)
      CALL fill_default(c_f(:), defaults%c_f)
      CALL fill_default(c_n(:), defaults%c_n)
      CALL fill_default(wmc(:), defaults%wmc)
      CALL fill_default(fsl(:), defaults%fsl)
      CALL fill_default(fbl(:), defaults%fbl)
      CALL fill_default(lmix_max(:), defaults%lmix_max)
      CALL fill_default(z0m_min(:), defaults%z0m_min)
      CALL fill_default(z0m_ice(:), defaults%z0m_ice)
      CALL fill_default(z0m_oce(:), defaults%z0m_oce)
      CALL fill_default(turb(:), turb_map(defaults%turb))
      CALL fill_default(smag_constant(:), defaults%smag_constant)
      CALL fill_default(turb_prandtl(:), defaults%turb_prandtl)
      CALL fill_default(km_min(:), defaults%km_min)
      CALL fill_default(max_turb_scale(:), defaults%max_turb_scale)
      CALL fill_default(min_sfc_wind(:), defaults%min_sfc_wind)

      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, NML=turb_vdiff_nml)    ! write settings to temporary text file
      END IF
    END IF

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------

    vdiff_config(:)%lsfc_mom_flux = lsfc_mom_flux(:)
    vdiff_config(:)%lsfc_heat_flux = lsfc_heat_flux(:)
    vdiff_config(:)%pr0 = pr0(:)
    vdiff_config(:)%f_tau0 = f_tau0(:)
    vdiff_config(:)%f_tau_limit_fraction = f_tau_limit_fraction(:)
    vdiff_config(:)%f_theta_limit_fraction = f_theta_limit_fraction(:)
    vdiff_config(:)%f_tau_decay = f_tau_decay(:)
    vdiff_config(:)%f_theta_decay = f_theta_decay(:)
    vdiff_config(:)%ek_ep_ratio_stable = ek_ep_ratio_stable(:)
    vdiff_config(:)%ek_ep_ratio_unstable = ek_ep_ratio_unstable(:)
    vdiff_config(:)%c_f = c_f(:)
    vdiff_config(:)%c_n = c_n(:)
    vdiff_config(:)%wmc = wmc(:)
    vdiff_config(:)%fsl = fsl(:)
    vdiff_config(:)%fbl = fbl(:)
    vdiff_config(:)%lmix_max = lmix_max(:)
    vdiff_config(:)%z0m_min = z0m_min(:)
    vdiff_config(:)%z0m_ice = z0m_ice(:)
    vdiff_config(:)%z0m_oce = z0m_oce(:)
    vdiff_config(:)%smag_constant = smag_constant(:)
    vdiff_config(:)%turb_prandtl = turb_prandtl(:)
    vdiff_config(:)%km_min = km_min(:)
    vdiff_config(:)%max_turb_scale = max_turb_scale(:)
    vdiff_config(:)%min_sfc_wind = min_sfc_wind(:)

    DO jg = 1, SIZE(turb)
      vdiff_config(jg)%turb = turb_from_str(turb(jg))
    END DO

    !----------------------------------------------------
    ! 5. Compute derived values and sanity check
    !----------------------------------------------------

    CALL vdiff_config_update(vdiff_config(:))
    CALL vdiff_config_check(vdiff_config(:))

    !$ACC UPDATE DEVICE(vdiff_config) ASYNC(1)

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      iunit = open_tmpfile()
      WRITE(iunit, NML=turb_vdiff_nml)
      CALL store_and_close_namelist(iunit, 'turb_vdiff_nml')
    ENDIF

    !--------------------------------------------------------
    ! 7. write the contents of the namelist to an ASCII file
    !--------------------------------------------------------
    IF(my_process_is_stdio()) WRITE(nnml_output, nml=turb_vdiff_nml)

  END SUBROUTINE vdiff_read_namelist


  FUNCTION turb_from_str(str) RESULT(turb)
    CHARACTER(len=*), INTENT(IN) :: str
    INTEGER :: turb

    SELECT CASE (TRIM(tolower(str)))
    CASE ('tte')
      turb = VDIFF_TURB_TTE
    CASE ('3dsmag', '3dsmagorinsky', 'smag', 'smagorinsky')
      turb = VDIFF_TURB_3DSMAGORINSKY
    CASE DEFAULT
      CALL finish(thismodule//'turb_from_str', 'Unknown `turb` value: ' // str)
    END SELECT
  END FUNCTION turb_from_str

  SUBROUTINE fill_default_r(val, default)
    REAL(wp), INTENT(INOUT) :: val(:)
    REAL(wp), INTENT(IN) :: default

    INTEGER :: i
    REAL(wp) :: parent_value

    parent_value = default

    DO i = 1, SIZE(val)
      IF (val(i) == -1._wp) THEN
        val(i) = parent_value
      ELSE
        parent_value = val(i)
      END IF
    END DO
  END SUBROUTINE fill_default_r

  SUBROUTINE fill_default_str(val, default)
    CHARACTER(len=*), INTENT(INOUT) :: val(:)
    CHARACTER(len=*), INTENT(IN) :: default

    INTEGER :: i
    CHARACTER(len=LEN(val)) :: parent_value

    parent_value = default

    DO i = 1, SIZE(val)
      IF (LEN_TRIM(val(i)) == 0) THEN
        val(i) = parent_value
      ELSE
        parent_value = val(i)
      END IF
    END DO
  END SUBROUTINE fill_default_str

END MODULE mo_turb_vdiff_nml
