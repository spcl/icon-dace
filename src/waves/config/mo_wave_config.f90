! Waves config.
! The content is mostly adopted from the WAM 4.5.4.
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

MODULE mo_wave_config

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: finish, message, message_text
  USE mo_impl_constants,       ONLY: max_dom, SUCCESS, MAX_CHAR_LENGTH
  USE mo_math_constants,       ONLY: pi2, rad2deg, dbl_eps
  USE mo_physical_constants,   ONLY: grav, rhoh2o
  USE mo_wave_constants,       ONLY: EX_TAIL
  USE mo_fortran_tools,        ONLY: DO_DEALLOCATE, t_ptr_1d_int
  USE mo_io_units,             ONLY: filename_max
  USE mo_util_string,          ONLY: t_keyword_list, associate_keyword, with_keywords, &
    &                                int2string

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_config'

  ! types
  PUBLIC :: t_wave_config

  ! objects
  PUBLIC :: wave_config

  ! subroutines
  PUBLIC :: configure_wave
  PUBLIC :: generate_filename

  TYPE t_wave_config
    INTEGER  :: ndirs    ! number of directions.
    INTEGER  :: nfreqs   ! number of frequencies.

    REAL(wp) :: fr1      ! first frequency [hz].
    REAL(wp) :: co       ! frequency ratio
    INTEGER  :: iref     ! frequency bin number of reference frequency

    REAL(wp) :: alpha      ! phillips' parameter  (not used if iopti = 1)
    REAL(wp) :: fm         ! peak frequency (hz) and/or maximum frequency
    REAL(wp) :: gamma_wave ! overshoot factor
    REAL(wp) :: sigma_a    ! left peak width
    REAL(wp) :: sigma_b    ! right peak width
    REAL(wp) :: fetch      ! fetch in metres (if zero then 0.5 of the latitude increment is used.).

    REAL(wp) :: roair   ! air density
    REAL(wp) :: rnuair  ! kinematic air viscosity
    REAL(wp) :: rnuairm ! kinematic air viscosity for momentum transfer
    REAL(wp) :: rowater ! water density
    REAL(wp) :: xeps
    REAL(wp) :: xinveps

    REAL(wp) :: xkappa   ! von karman constant.
    REAL(wp) :: xnlev    ! windspeed ref. level.
    REAL(wp) :: betamax  ! parameter for wind input (ecmwf cy45r1).
    REAL(wp) :: zalp     ! shifts growth curve (ecmwf cy45r1).
    REAL(wp) :: alpha_ch ! minimum charnock constant (ecmwf cy45r1)

    REAL(wp) :: depth     ! ocean depth (m) if not 0, then constant depth
    REAL(wp) :: depth_min ! allowed minimum of model depth (m)
    REAL(wp) :: depth_max ! allowed maximum of model depth (m)

    INTEGER  :: niter_smooth ! number of smoothing iterations for wave bathymetry
                             ! if 0 then no smoothing

    INTEGER  :: jtot_tauhf ! dimension of wtauhf, must be odd
    REAL(wp) :: x0tauhf    ! lowest limit for integration in tau_phi_hf: x0 *(g/ustar)

    CHARACTER(LEN=filename_max) :: forc_file_prefix ! prefix of forcing file name
                                           ! the real file name will be constructed as:
                                           ! forc_file_prefix+'_wind.nc' for U and V 10 meter wind (m/s)
                                           ! forc_file_prefix+'_ice.nc'  for sea ice concentration (fraction of 1)
                                           ! forc_file_prefix+'_slh.nc'  for sea level height (m)
                                           ! forc_file_prefix+'_osc.nc'  for U and V ocean surface currents (m/s)

    LOGICAL :: linput_sf1      ! if .TRUE., calculate wind input source function term, first call
    LOGICAL :: linput_sf2      ! if .TRUE., calculate wind input source function term, second call
    LOGICAL :: ldissip_sf      ! if .TRUE., calculate dissipation source function term
    LOGICAL :: lwave_brk_sf    ! if .TRUE., calculate wave breaking dissipation source function term
    LOGICAL :: lnon_linear_sf  ! if .TRUE., calculate non linear source function term
    LOGICAL :: lbottom_fric_sf ! if .TRUE., calculate bottom_friction source function term
    LOGICAL :: lwave_stress1   ! if .TRUE., calculate wave stress, first call
    LOGICAL :: lwave_stress2   ! if .TRUE., calculate wave stress, second call

    REAL(wp) :: peak_u10, peak_v10 ! peak value (m/s) of 10 m U and V wind speed for test case
    REAL(wp) :: peak_lat, peak_lon ! geographical location (deg) of wind peak value

    REAL(wp) :: impl_fac       !! implicitness factor for total source function time integration
                               !! impl_fac=0.5 : second order Crank-Nicholson/trapezoidal scheme
                               !! impl_fac=1   : first order Euler backward scheme
                               !! valid range: 0.5 <= impl_fac <= 1

    ! derived variables and fields
    !
    REAL(wp) ::            &
      &  delth,            & ! angular increment of spectrum [rad].
      &  mo_tail,          & ! mo  tail factor.
      &  mm1_tail,         & ! m-1 tail factor.
      &  mp1_tail,         & ! m+1 tail factor.
      &  mp2_tail,         & ! m+2 tail factor.
      &  acl1,             & ! weight in angular grid for interpolation,
                             ! wave no. 3 ("1+lambda" term).
      &  acl2,             & ! weight in angular grid for interpolation,
                             ! wave no. 4 ("1-lambda" term).
      &  cl11,             & ! 1.-acl1.
      &  cl21,             & ! 1.-acl2.
      &  dal1,             & ! 1./acl1.
      &  dal2,             & ! 1./acl2.
      &  frh(30)             ! tail frequency ratio **5

    REAL(wp), ALLOCATABLE :: &
      &  freqs(:),         & ! frequencies (1:nfreqs) of wave spectrum [hz]
      &  dfreqs(:),        & ! frequency interval (1:nfreqs)
      &  dfreqs_freqs(:),  & ! dfreqs * freqs
      &  dfreqs_freqs2(:), & ! dfreqs * freqs * freqs
      &  dirs(:),          & ! directions (1:ndirs) of wave spectrum [rad]
      &  dfim(:),          & ! mo  integration weights.
      &  dfimofr(:),       & ! m-1 integration weights.
      &  dfim_fr(:),       & ! m+1 integration weights.
      &  dfim_fr2(:),      & ! m+2 integration weights.
      &  rhowg_dfim(:),    & ! momentum and energy flux weights.
      &  wtauhf(:)           ! integration weight for tau_phi_hf

    INTEGER, ALLOCATABLE :: &
      &  freq_ind(:),      & ! index of frequency as a function of tracer index (1:ntracer=ndirs*nfreq)
      &  dir_ind(:),       & ! index of direction as a function of tracer index (1:ntracer=ndirs*nfreq)
      &  tracer_ind(:,:),  & ! tracer index as a function of direction and frequency index
      &  dir_neig_ind(:,:)   ! index of direction neighbor (2,1:ndirs)

    TYPE(t_ptr_1d_int), ALLOCATABLE :: & ! list of tracer IDs as a function of the frequency index
      &  list_tr(:)

    LOGICAL :: lread_forcing ! set to .TRUE. if a forcing file prefix has been specified (forc_file_prefix)

  CONTAINS
    !
    ! destruct wave_config object
    PROCEDURE :: destruct            => wave_config_destruct
  END type t_wave_config

  TYPE(t_wave_config), TARGET:: wave_config(max_dom)

CONTAINS

  !>
  !! deallocate memory used by object of type t_wave_config
  !!
  SUBROUTINE wave_config_destruct(me)
    CLASS(t_wave_config) :: me
    INTEGER :: jf
    INTEGER :: ist
    CHARACTER(*), PARAMETER :: routine = modname//'::wave_config_destruct'

    CALL DO_DEALLOCATE(me%freqs)
    CALL DO_DEALLOCATE(me%dfreqs)
    CALL DO_DEALLOCATE(me%dfreqs_freqs)
    CALL DO_DEALLOCATE(me%dfreqs_freqs2)
    CALL DO_DEALLOCATE(me%dirs)
    CALL DO_DEALLOCATE(me%DFIM)
    CALL DO_DEALLOCATE(me%DFIMOFR)
    CALL DO_DEALLOCATE(me%DFIM_FR)
    CALL DO_DEALLOCATE(me%DFIM_FR2)
    CALL DO_DEALLOCATE(me%RHOWG_DFIM)
    CALL DO_DEALLOCATE(me%freq_ind)
    CALL DO_DEALLOCATE(me%dir_ind)
    CALL DO_DEALLOCATE(me%tracer_ind)
    CALL DO_DEALLOCATE(me%dir_neig_ind)
    CALL DO_DEALLOCATE(me%wtauhf)

    DO jf=1,SIZE(me%list_tr)
      me%list_tr(jf)%p => NULL()
    ENDDO
    DEALLOCATE(me%list_tr, stat=ist)
    IF (ist/=SUCCESS) CALL finish(routine, "deallocation failed for list_tr")

  END SUBROUTINE wave_config_destruct


  !>
  !! setup the waves model
  !!
  !! Setup of additional waves control variables and constant fields
  !! which depend on the waves-NAMELIST and potentially other namelists.
  !! This routine is called, after all namelists have been read and a
  !! synoptic consistency check has been done.
  !!
  SUBROUTINE configure_wave(n_dom, ntracer)

    INTEGER, INTENT(IN) :: n_dom    !< number of domains
    INTEGER, INTENT(IN) :: ntracer  !< number of tracers

    INTEGER :: jd, jf, jfjd   ! loop index
    INTEGER :: jg             ! patch ID
    INTEGER :: ist            ! error status
    INTEGER :: j
    TYPE(t_wave_config), POINTER :: wc =>NULL()     ! convenience pointer

    REAL(wp) :: CO1, X0, FF, F, DF, CONST1

    CHARACTER(*), PARAMETER :: routine = modname//'::configure_waves'

    CALL message (' ','')
    CALL message ('!','')
    CALL message ('!!','')
    CALL message ('!!!!','')
    CALL message ('!!!!!!','')
    CALL message ('!!!!!!!!','')
    CALL message ('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!','')
    CALL message ('!     ICON-WAVES IS STILL UNDER DEVELOPMENT                              !','')
    CALL message ('!     THIS CODE IS FOR TECHNICAL TESTING ONLY                            !','')
    CALL message ('!     THE IMPLEMENTATION OF WAVE PHYSICS IS NOT YET COMPLETED            !','')
    CALL message ('!     THE ENERGY PROPAGATION DOES NOT WORK FOR THE 3 OUTERMOST CELL ROWS !','')
    CALL message ('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!','')
    CALL message ('!!!!!!!!','')
    CALL message ('!!!!!!','')
    CALL message ('!!!!','')
    CALL message ('!!','')
    CALL message ('!','')
    CALL message (' ','')

    DO jg=1,n_dom

      ! convenience pointer
      wc => wave_config(jg)


      ! reading of external forcing data yes/no
      ! set to .TRUE. if a forcing file prefix has been specified
      IF (TRIM(wave_config(jg)%forc_file_prefix) /= '') THEN
        wc%lread_forcing = .TRUE.
      ELSE
        wc%lread_forcing = .FALSE.
      ENDIF

      ALLOCATE(wc%dirs         (wc%ndirs),  &
        &      wc%freqs        (wc%nfreqs), &
        &      wc%dfreqs       (wc%nfreqs), &
        &      wc%dfreqs_freqs (wc%nfreqs), &
        &      wc%dfreqs_freqs2(wc%nfreqs), &
        &      wc%dfim         (wc%nfreqs), &
        &      wc%dfimofr      (wc%nfreqs), &
        &      wc%dfim_fr      (wc%nfreqs), &
        &      wc%dfim_fr2     (wc%nfreqs), &
        &      wc%rhowg_dfim   (wc%nfreqs), &
        &      wc%wtauhf   (wc%jtot_tauhf), &
        &      stat=ist)
      IF (ist/=SUCCESS) CALL finish(routine, "allocation for fields of type REAL failed")

      ALLOCATE(wc%freq_ind     (ntracer),            &
        &      wc%dir_ind      (ntracer),            &
        &      wc%tracer_ind   (wc%ndirs,wc%nfreqs), &
        &      wc%dir_neig_ind (2,wc%ndirs),         &
        &      wc%list_tr      (wc%nfreqs),          &
        &      stat=ist)
      IF (ist/=SUCCESS) CALL finish(routine, "allocation for fields of type INTEGER failed")


      !
      ! configuration of spectral setup
      !
      CALL message ('  ','')
      CALL message (':----------------------------------------------------------','')
      WRITE(message_text,'(a,i4)') 'grid ', jg
      CALL message ('Frequencies and directions of wave spectrum',message_text)
      WRITE(message_text,'(a,i4)') 'Number of directions  = ', wc%ndirs
      CALL message ('  ',message_text)
      WRITE(message_text,'(a,i4)') 'Number of frequencies = ', wc%nfreqs
      CALL message ('  ',message_text)
      CALL message (':----------------------------------------------------------','')

      wc%DELTH    = pi2 / REAL(wc%ndirs,wp) !! ANGULAR INCREMENT OF SPECTRUM [RAD].

      ! calculate directions for wave spectrum
      !
      CALL message ('  ','Directions [Degree]: ')
      DO jd = 1,wc%ndirs
        wc%dirs(jd) = REAL(jd-1,wp) *  wc%DELTH + 0.5_wp * wc%DELTH !RAD
        WRITE(message_text,'(i3,f10.5)') jd, wc%dirs(jd)*rad2deg
        CALL message ('  ',message_text)
      END DO
      CALL message ('  ','')

      ! calculate frequencies for wave spectrum
      !
      CALL message ('  ','Frequencies, [Hz]: ')
      WRITE(message_text,'(i3,f10.5)') 1, wc%FR1
      CALL message ('  ',message_text)
      wc%freqs(1) = wc%CO**(-wc%iref + 1) * wc%FR1
      DO jf = 2,wc%nfreqs
        wc%freqs(jf) = wc%CO * wc%freqs(jf-1)
        WRITE(message_text,'(i3,f10.5)') jf, wc%freqs(jf)
        CALL message ('  ',message_text)
      END DO
      CALL message (':----------------------------------------------------------','')


      ! calculate frequency intervals
      !
      CO1 = 0.5_wp * (wc%CO - 1.0_wp)
      wc%dfreqs(1) = CO1 * wc%freqs(1)
      wc%dfreqs(2:wc%nfreqs-1) = CO1 &
        &                      * (wc%freqs(2:wc%nfreqs-1) &
        &                      + (wc%freqs(1:wc%nfreqs-2)))
      wc%dfreqs(wc%nfreqs) = CO1 * wc%freqs(wc%nfreqs-1)
      !
      wc%dfreqs_freqs  = wc%dfreqs * wc%freqs
      wc%dfreqs_freqs2 = wc%dfreqs_freqs * wc%freqs

      ! compute mappings between direction IDs, frequency IDs and tracer IDs
      !
      DO jf = 1,wc%nfreqs
        DO jd = 1,wc%ndirs

          ! tracer ID as a function of frequency ID and direction ID
          jfjd = jd + (jf-1) * wc%ndirs
          wc%tracer_ind(jd,jf) = jfjd

          ! frequency and direction IDs as a function of tracer ID
          wc%freq_ind(jfjd) = jf
          wc%dir_ind(jfjd)  = jd

        END DO
        ! list of tracer IDs as a function of frequency ID
        wc%list_tr(jf)%p => wc%tracer_ind(:,jf)
      END DO

      ! calculate direction neighbor index
      DO jd = 1,wc%ndirs
        IF (jd == 1) THEN
          wc%dir_neig_ind(1,jd) = wc%ndirs
          wc%dir_neig_ind(2,jd) = jd+1
        ELSEIF (jd == wc%ndirs) THEN
          wc%dir_neig_ind(1,jd) = jd-1
          wc%dir_neig_ind(2,jd) = 1
        ELSE
          wc%dir_neig_ind(1,jd) = jd-1
          wc%dir_neig_ind(2,jd) = jd+1
        END IF
      END DO


      ! MO  TAIL FACTOR.
      wc%MO_TAIL  = - wc%DELTH / (EX_TAIL + 1.0_wp) * wc%freqs(wc%nfreqs)

      ! M-1 TAIL FACTOR.
      wc%MM1_TAIL = - wc%DELTH / EX_TAIL

      ! M+1 TAIL FACTOR.
      wc%MP1_TAIL = - wc%DELTH / (EX_TAIL + 2.0_wp) * wc%freqs(wc%nfreqs)**2

      ! M+2 TAIL FACTOR.
      wc%MP2_TAIL = - wc%DELTH / (EX_TAIL + 3.0_wp) * wc%freqs(wc%nfreqs)**3


      ! calculate freqs_dirs parameters
      ! MO INTEGRATION WEIGHTS
      wc%DFIM = wc%dfreqs * wc%DELTH

      ! M-1 INTEGRATION WEIGHTS.
      wc%DFIMOFR = wc%DFIM / wc%freqs

      ! M+1 INTEGRATION WEIGHTS.
      wc%DFIM_FR = wc%dfreqs_freqs * wc%DELTH

      ! M+2 INTEGRATION WEIGHTS.
      wc%DFIM_FR2 = wc%dfreqs_freqs2 * wc%DELTH

      ! MOMENTUM AND ENERGY FLUX WEIGHTS.
      wc%RHOWG_DFIM(:) = rhoh2o * grav * wc%DELTH * LOG(wc%CO) * wc%freqs(:)
      !
      wc%RHOWG_DFIM(1)         = 0.5_wp * wc%RHOWG_DFIM(1)
      wc%RHOWG_DFIM(wc%nfreqs) = 0.5_wp * wc%RHOWG_DFIM(wc%nfreqs)

      ! initialisation of wtauhf
      wc%wtauhf(:) = 0._wp

      X0 = 0.005_wp
      FF = EXP(wc%XKAPPA / (X0 + wc%ZALP))
      F = wc%ALPHA_CH * X0**2 * FF -1.0_wp

      j = 1
      DO WHILE(ABS(F)>dbl_eps .AND. j<30)
        FF = EXP(wc%XKAPPA / (X0 + wc%ZALP))
        F = wc%ALPHA_CH * X0**2 * FF -1.0_wp
        DF = wc%ALPHA_CH * FF *(2.0_wp * X0 - wc%XKAPPA * (X0 / (X0 + wc%ZALP))**2)
        X0 = X0 - F / DF
        j = j + 1
      END DO
      wc%X0TAUHF = X0

      CONST1 = (wc%BETAMAX/wc%XKAPPA**2)/3.0_wp

      ! Simpson integration weights
      wc%WTAUHF(1) = CONST1
      DO J = 2, wc%jtot_tauhf-1,2
        wc%WTAUHF(J) = 4.0_wp * CONST1
        wc%WTAUHF(J+1) = 2.0_wp * CONST1
      END DO
      wc%WTAUHF(wc%jtot_tauhf) = CONST1

    END DO

  END SUBROUTINE configure_wave


  FUNCTION generate_filename(input_filename, model_base_dir, &
    &                        nroot, jlev, idom)  RESULT(result_str)
    CHARACTER(len=*), INTENT(IN)   :: input_filename, &
      &                               model_base_dir
    INTEGER,          INTENT(IN)   :: nroot, jlev, idom
    CHARACTER(len=MAX_CHAR_LENGTH) :: result_str
    TYPE (t_keyword_list), POINTER :: keywords => NULL()

    CALL associate_keyword("<path>",   TRIM(model_base_dir),             keywords)
    CALL associate_keyword("<nroot>",  TRIM(int2string(nroot,"(i0)")),   keywords)
    CALL associate_keyword("<nroot0>", TRIM(int2string(nroot,"(i2.2)")), keywords)
    CALL associate_keyword("<jlev>",   TRIM(int2string(jlev, "(i2.2)")), keywords)
    CALL associate_keyword("<idom>",   TRIM(int2string(idom, "(i2.2)")), keywords)
    ! replace keywords in "input_filename", which is by default
    ! ifs2icon_filename = "<path>ifs2icon_R<nroot>B<jlev>_DOM<idom>.nc"
    result_str = TRIM(with_keywords(keywords, TRIM(input_filename)))

  END FUNCTION generate_filename

END MODULE mo_wave_config
