!
! Read and apply optical properties of aerosol climatology
!        for volcanic stratospheric aerosols as provided for CMIP6
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

MODULE mo_bc_aeropt_cmip6_volc

  USE mo_kind,                   ONLY: wp, i8
  USE mo_exception,              ONLY: finish, message, message_text
  USE mo_read_interface,         ONLY: openInputFile, closeFile, read_1D, &
    &                                  read_extdim_slice_extdim_extdim_extdim
  USE mo_latitude_interpolation, ONLY: latitude_weights_li
  USE mo_math_constants,         ONLY: deg2rad, pi_2
  USE mtime,                     ONLY: datetime
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights, &
       &                               calculate_time_interpolation_weights
  USE mo_time_config,            ONLY: time_config

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_bc_aeropt_cmip6_volc
  PUBLIC :: add_bc_aeropt_cmip6_volc

  ! Data file layout.

  !> Prefix of the filename from which the aerosol data is read.
  CHARACTER(len=*), PARAMETER :: filename_base = 'bc_aeropt_cmip6_volc_lw_b16_sw_b14_'
  !> Name of the latitude dimension.
  CHARACTER(len=*), PARAMETER :: dim_name_lat = 'latitude'
  !> Name of the altitude dimension.
  CHARACTER(len=*), PARAMETER :: dim_name_alt = 'altitude'
  !> Name of the month dimension.
  CHARACTER(len=*), PARAMETER :: dim_name_month = 'month'
  !> Name of the solar band dimension.
  CHARACTER(len=*), PARAMETER :: dim_name_sband = 'solar_bands'
  !> Name of terrestrial band dimension.
  CHARACTER(len=*), PARAMETER :: dim_name_tband = 'terrestrial_bands'

  INTEGER, PARAMETER :: max_dim_name_len = &
      & MAX(LEN(dim_name_lat), LEN(dim_name_alt), LEN(dim_name_month), LEN(dim_name_sband), &
      &   LEN(dim_name_tband))
  !> Dimension layout for solar data.
  CHARACTER(len=max_dim_name_len), PARAMETER :: dim_names_sol(4) = [ &
      & CHARACTER(len=max_dim_name_len) :: &
      & dim_name_month, dim_name_alt, dim_name_lat, dim_name_sband &
    ]
  !> Dimension layout for terrestrial data.
  CHARACTER(len=max_dim_name_len), PARAMETER :: dim_names_terr(4) = [ &
      & CHARACTER(len=max_dim_name_len) :: &
      & dim_name_month, dim_name_alt, dim_name_lat, dim_name_tband &
    ]

  !> Marker for `pre_year` to show that arrays do not contain valid data.
  INTEGER(i8), PARAMETER :: PRE_YEAR_UNINITIALIZED = -HUGE(1_i8)

  REAL(wp), ALLOCATABLE :: aod_v_s(:,:,:)   !< volcanic AOD solar (nbndsw,nlats,nmonths).
  REAL(wp), ALLOCATABLE :: ssa_v_s(:,:,:,:) !< volcanic SSA solar (nbndsw,nalts,nlats,nmonths).
  REAL(wp), ALLOCATABLE :: asy_v_s(:,:,:,:) !< volcanic ASY solar (nbndsw,nalts,nlats,nmonths).
  REAL(wp), ALLOCATABLE :: ext_v_s(:,:,:,:) !< volcanic EXT solar (nbndsw,nalts,nlats,nmonths).
  REAL(wp), ALLOCATABLE :: aod_v_t(:,:,:)   !< volcanic AOD therm (nbndlw,nlats,nmonths).
  REAL(wp), ALLOCATABLE :: ssa_v_t(:,:,:,:) !< volcanic SSA therm (nbndlw,nalts,nlats,nmonths).
  REAL(wp), ALLOCATABLE :: ext_v_t(:,:,:,:) !< volcanic EXT therm (nbndlw,nalts,nlats,nmonths).
  REAL(wp), ALLOCATABLE :: r_alt_clim(:)    !< altitude [m] (nalts).
  REAL(wp), ALLOCATABLE :: r_lat_clim(:)    !< latitudes [rad] (nlats).
  REAL(wp)              :: r_lat_shift, r_rdeltalat
  INTEGER               :: k_alt_clim, lat_clim

  INTEGER(i8)           :: pre_year = PRE_YEAR_UNINITIALIZED
  INTEGER(i8)           :: nyears
  INTEGER               :: imonth_beg, imonth_end

  LOGICAL               :: lend_of_year

  TYPE(t_time_interpolation_weights) :: tiw_beg
  TYPE(t_time_interpolation_weights) :: tiw_end

CONTAINS


  !> Set up memory for fields in which the aerosol optical properties are stored when needed.
  SUBROUTINE su_bc_aeropt_cmip6_volc(nbndlw, nbndsw)

    INTEGER, INTENT(IN) :: nbndlw !< Number of long-wave bands in the data files.
    INTEGER, INTENT(IN) :: nbndsw !< Number of short-wave bands in the data files.

    ! String with enough space for filename and 5-digit year.
    CHARACTER(len=LEN(filename_base) + 5 + 3) :: filename

    INTEGER :: file_id

    REAL(wp), POINTER :: zlat(:), zalt(:)

    CHARACTER(len=*), PARAMETER :: subroutine_name = &
        & 'mo_bc_aeropt_cmip6_volc:su_bc_aeropt_cmip6_volc'

    IF (ALLOCATED(aod_v_s)) RETURN

    lend_of_year = ( time_config%tc_stopdate%date%month  == 1  .AND. &
      &              time_config%tc_stopdate%date%day    == 1  .AND. &
      &              time_config%tc_stopdate%time%hour   == 0  .AND. &
      &              time_config%tc_stopdate%time%minute == 0  .AND. &
      &              time_config%tc_stopdate%time%second == 0 )

    nyears = time_config%tc_stopdate%date%year - time_config%tc_startdate%date%year + 1
    IF ( lend_of_year ) nyears = nyears - 1

    ! ----------------------------------------------------------------------

    tiw_beg = calculate_time_interpolation_weights(time_config%tc_startdate)
    tiw_end = calculate_time_interpolation_weights(time_config%tc_stopdate)

    IF ( nyears > 1 ) THEN
      imonth_beg = 0
      imonth_end = 13
    ELSE
      imonth_beg = tiw_beg%month1
      imonth_end = tiw_end%month2
      ! special case for runs starting on 1 Jan that run for less than a full year
      IF ( imonth_beg == 12 .AND. time_config%tc_startdate%date%month == 1 ) imonth_beg = 0
      ! special case for runs ending in 2nd half of Dec that run for less than a full year
      IF ( lend_of_year .OR. ( imonth_end == 1 .AND. time_config%tc_stopdate%date%month == 12 ) ) &
          & imonth_end = 13
    ENDIF

    WRITE (message_text,'(a,i2,a,i2)') &
      & ' Allocating CMIP6 volcanic aerosols for months ', imonth_beg, ' to ', imonth_end
    CALL message(subroutine_name, message_text)

    ! Get altitude and latitude points from the first file. These are assumed to be constant
    ! throughout the run.
    WRITE (filename,'(a,i0,a)') filename_base, time_config%tc_startdate%date%year, '.nc'

    CALL openInputFile(file_id, filename)
    CALL read_1D(file_id=file_id, variable_name=dim_name_alt, return_pointer=zalt)
    CALL read_1D(file_id=file_id, variable_name=dim_name_lat, return_pointer=zlat)
    CALL closeFile(file_id)

    k_alt_clim = SIZE(zalt)
    lat_clim = SIZE(zlat)

    ALLOCATE(aod_v_s(nbndsw,0:lat_clim+1,imonth_beg:imonth_end))
    ALLOCATE(ext_v_s(nbndsw,k_alt_clim+1,0:lat_clim+1,imonth_beg:imonth_end))
    ALLOCATE(ssa_v_s(nbndsw,k_alt_clim+1,0:lat_clim+1,imonth_beg:imonth_end))
    ALLOCATE(asy_v_s(nbndsw,k_alt_clim+1,0:lat_clim+1,imonth_beg:imonth_end))
    ALLOCATE(aod_v_t(nbndlw,0:lat_clim+1,imonth_beg:imonth_end))
    ALLOCATE(ext_v_t(nbndlw,k_alt_clim+1,0:lat_clim+1,imonth_beg:imonth_end))
    ALLOCATE(ssa_v_t(nbndlw,k_alt_clim+1,0:lat_clim+1,imonth_beg:imonth_end))
    ALLOCATE(r_alt_clim(k_alt_clim))
    ALLOCATE(r_lat_clim(0:lat_clim+1))

    !$ACC ENTER DATA CREATE(aod_v_s, ext_v_s, ssa_v_s, asy_v_s, aod_v_t, ext_v_t, ssa_v_t) &
    !$ACC   CREATE(r_alt_clim, r_lat_clim)

    aod_v_s(:,:,:) = 0._wp
    ext_v_s(:,:,:,:) = 0._wp
    ssa_v_s(:,:,:,:) = 0._wp
    asy_v_s(:,:,:,:) = 0._wp
    aod_v_t(:,:,:) = 0._wp
    ext_v_t(:,:,:,:) = 0._wp
    ssa_v_t(:,:,:,:) = 0._wp

    ! reverse and convert from km to m.
    r_alt_clim(:) = 1000._wp * zalt(k_alt_clim:1:-1)

    r_lat_clim(1:lat_clim) = deg2rad * zlat(lat_clim:1:-1)
    r_lat_clim(0) = 0.0_wp
    r_lat_clim(lat_clim+1) = -pi_2
    r_lat_shift = r_lat_clim(1)     ! this is the value next to the N-pole (so +87.5 for example)
    r_rdeltalat = ABS(1.0_wp/(r_lat_clim(2)-r_lat_clim(1)))

    DEALLOCATE(zalt, zlat)

    !$ACC UPDATE DEVICE(aod_v_s, ext_v_s, ssa_v_s, asy_v_s, aod_v_t, ext_v_t, ssa_v_t) &
    !$ACC   DEVICE(r_alt_clim, r_lat_clim) &
    !$ACC   ASYNC(1)

  END SUBROUTINE su_bc_aeropt_cmip6_volc


  !> Shifts December of current year into imonth=0 and January of the following year into imonth=1
  !! (these months do not need to be read again).
  SUBROUTINE shift_months_bc_aeropt_cmip6_volc()

    CHARACTER(len=*), PARAMETER :: subroutine_name = &
        & 'mo_bc_aeropt_cmip6_volc:shift_months_bc_aeropt_cmip6_volc'

    IF ( .NOT. ALLOCATED(aod_v_s) ) CALL finish(subroutine_name, 'data arrays are not allocated')

    IF ( imonth_beg > 0 .OR. imonth_end < 13 ) THEN
      WRITE (message_text,'(a,i2,a,i2)') ' CMIP6 volcanic aerosols are allocated for months ', &
          & imonth_beg, ' to ', imonth_end, 'only.'
      CALL message(subroutine_name, message_text)
      CALL finish(subroutine_name, &
          & ' CMIP6 volcanic aerosols are not allocated over required range 0 to 13.')
    ENDIF

    WRITE (message_text,'(a)') &
        & ' Copy CMIP6 volcanic aerosols for months 12:13 to months 0:1'
    CALL message(subroutine_name, message_text)

    aod_v_s(:,:,0:1) = aod_v_s(:,:,12:13)
    ext_v_s(:,:,:,0:1) = ext_v_s(:,:,:,12:13)
    ssa_v_s(:,:,:,0:1) = ssa_v_s(:,:,:,12:13)
    asy_v_s(:,:,:,0:1) = asy_v_s(:,:,:,12:13)
    aod_v_t(:,:,0:1) = aod_v_t(:,:,12:13)
    ext_v_t(:,:,:,0:1) = ext_v_t(:,:,:,12:13)
    ssa_v_t(:,:,:,0:1) = ssa_v_t(:,:,:,12:13)

  END SUBROUTINE shift_months_bc_aeropt_cmip6_volc


  !> Read optical properties of CMIP6 volcanic aerosols from external files into interpolation
  !! cache. This routine has to be called on initialization and on Jan 1 of each simulation year.
  SUBROUTINE read_bc_aeropt_cmip6_volc(mtime_current, nbndlw, nbndsw)

    TYPE(datetime), POINTER, INTENT(IN) :: mtime_current !< Current date.
    INTEGER, INTENT(IN) :: nbndsw !< Number of short-wave bands.
    INTEGER, INTENT(in) :: nbndlw !< Number of long-wave bands.

    ! LOCAL VARIABLES
    INTEGER(i8)                   :: iyear
    INTEGER                       :: imonthb, imonthe

    iyear = mtime_current%date%year

    IF (iyear > pre_year) THEN

      ! beginning of job or change of year

      IF ( pre_year > PRE_YEAR_UNINITIALIZED ) THEN
        CALL shift_months_bc_aeropt_cmip6_volc()
      ELSE
        CALL su_bc_aeropt_cmip6_volc(nbndlw, nbndsw)
      ENDIF

      ! Restrict reading of data to those months that are needed

      IF ( nyears > 1 ) THEN

        IF ( pre_year > PRE_YEAR_UNINITIALIZED ) THEN
          ! second and following years of current run
          imonthb = 2
        ELSE
          ! first year of current run
          imonthb = tiw_beg%month1
          IF ( imonthb == 12 .AND. time_config%tc_startdate%date%month == 1 ) imonthb = 0
        ENDIF

        IF ( mtime_current%date%year < time_config%tc_stopdate%date%year ) THEN
          imonthe = 13
        ELSE

          IF ( tiw_end%month2 == 1 .AND. time_config%tc_stopdate%date%month == 12 ) THEN
              imonthe = 13
          ELSE
              imonthe = tiw_end%month2
              ! no reading of month 2 if end is already before 15 Jan.
              IF ( imonthb == 2 .AND. imonthe < imonthb ) THEN
                pre_year = mtime_current%date%year
                RETURN
              ENDIF
          ENDIF

        ENDIF

      ELSE

        ! only less or equal one year in current run
        ! we can savely narrow down the data that have to be read in.

        imonthb = tiw_beg%month1
        imonthe = tiw_end%month2
        ! special case for runs starting on 1 Jan that run for less than a full year
        IF ( imonthb == 12 .AND. time_config%tc_startdate%date%month == 1  ) imonthb = 0
        ! special case for runs ending in 2nd half of Dec that run for less than a full year
        IF ( lend_of_year .OR. ( imonthe == 1 .AND. time_config%tc_stopdate%date%month == 12 ) ) &
            & imonthe = 13

      ENDIF

      CALL read_months_bc_aeropt_cmip6_volc(imonthb, imonthe, iyear)

      pre_year = mtime_current%date%year

      !$ACC UPDATE DEVICE(aod_v_s, ext_v_s, ssa_v_s, asy_v_s, aod_v_t, ext_v_t, ssa_v_t) &
      !$ACC   DEVICE(r_alt_clim, r_lat_clim) &
      !$ACC   ASYNC(1)

    END IF ! iyear > pre_year
  END SUBROUTINE read_bc_aeropt_cmip6_volc


  !> Add aerosol optical properties of CMIP6 volcanic aerosols to all wave length bands (solar and
  !! IR). The height profile is taken into account.
  SUBROUTINE add_bc_aeropt_cmip6_volc( &
      & current_date,          jg,              jcs,                  &
      & kproma,                kbdim,           klev,                 &
      & krow,                  nb_sw,           nb_lw,                &
      & zf,                    dz,                                    &
      & paer_tau_sw_vr,        paer_piz_sw_vr,  paer_cg_sw_vr,        &
      & paer_tau_lw_vr                                                )

    ! INPUT PARAMETERS
    TYPE(datetime), POINTER, INTENT(IN) :: current_date !< Current date and time.
    INTEGER, INTENT(IN) :: jg !< Domain index (for cell -> latitude mapping).
    INTEGER, INTENT(IN) :: jcs !< Start index in block.
    INTEGER, INTENT(IN) :: kproma !< Actual block length.
    INTEGER, INTENT(IN) :: kbdim !< Maximum block length.
    INTEGER, INTENT(IN) :: krow !< Block index.
    INTEGER, INTENT(IN) :: klev !< Number of vertical levels.
    INTEGER, INTENT(IN) :: nb_lw !< Number of wave-length bands (far IR).
    INTEGER, INTENT(IN) :: nb_sw !< Number of wave-length bands (solar).
    REAL(wp), INTENT(IN) :: dz(kbdim,klev) !< Geometric height thickness [m]
    REAL(wp), INTENT(IN) :: zf(kbdim,klev) !< Geometric height [m].

    ! OUTPUT PARAMETERS
    !> Aerosol optical depth (far IR).
    REAL(wp), INTENT(INOUT) :: paer_tau_lw_vr(kbdim,klev,nb_lw)
    !> aerosol optical depth (solar), sum_i(tau_i).
    REAL(wp), INTENT(INOUT) :: paer_tau_sw_vr(kbdim,klev,nb_sw)
    !> weighted sum of single scattering albedos, sum_i(tau_i*omega_i).
    REAL(wp), INTENT(INOUT) :: paer_piz_sw_vr(kbdim,klev,nb_sw)
    !> Weighted sum of asymmetry factors, sum_i(tau_i*omega_i*g_i).
    REAL(wp), INTENT(INOUT) :: paer_cg_sw_vr(kbdim,klev,nb_sw)

    ! LOCAL VARIABLES
    INTEGER, PARAMETER                    :: norder=-1 ! latitudes in climatology order from N->S
    INTEGER                               :: jl,jk,jki,jwl
    INTEGER                               :: idx_lat_1, idx_lat_2, idx_lev
    REAL(wp)                              :: w1_lat, w2_lat
    INTEGER,  DIMENSION(kbdim,klev)       :: kindex ! index field for pressure interpolation
    LOGICAL,  DIMENSION(kbdim,klev)       :: l_kindex
    REAL(wp), DIMENSION(kbdim)            :: wgt1_lat,wgt2_lat
    INTEGER,  DIMENSION(kbdim)            :: inmw1_lat, inmw2_lat
    REAL(wp), DIMENSION(kbdim,nb_sw)      :: zaod_s, zext_s_int, zfact_s
    REAL(wp), DIMENSION(kbdim,klev,nb_sw) :: zext_s, zomg_s, zasy_s
    REAL(wp), DIMENSION(kbdim,klev,nb_lw) :: zext_t, zomg_t
    REAL(wp), DIMENSION(kbdim,nb_lw)      :: zaod_t, zext_t_int, zfact_t
    REAL(wp)                              :: p_lat_shift, p_rdeltalat
    INTEGER                               :: jc
    REAL(wp)                              :: dz_clim, z_max_lim_clim

    TYPE(t_time_interpolation_weights) :: tiw
    INTEGER :: nm1, nm2

    CHARACTER(len=*), PARAMETER :: subroutine_name = &
        & 'mo_bc_aeropt_cmip6_volc:set_bc_aeropt_cmip6_volc'

    tiw = calculate_time_interpolation_weights(current_date)
    nm1 = tiw%month1_index
    nm2 = tiw%month2_index

    IF (current_date%date%year /= pre_year) THEN
      WRITE (message_text,'(A,I4,A,I4)') 'Stale data: requested year is', current_date%date%year, &
          & ' but data is for ', pre_year
      CALL finish(subroutine_name, message_text)
    END IF

    ! It is assumed that the pressure levels of the climatology do not change with time but
    ! are unequally spaced. Since the pressure of each icon level may change with time,
    ! each specific icon level may have its centre in a different level of the climatology at
    ! each time step.

    ! 1. calculate for each icon gridbox the index of the data set layer
    !     in which p_mid_icon is located and geometrical height of layers
    dz_clim=r_alt_clim(1)-r_alt_clim(2)
    z_max_lim_clim=r_alt_clim(1)+0.5_wp*dz_clim
    CALL altitude_index( &
        & jcs,            kproma,     zf,     dz_clim, &
        & z_max_lim_clim, k_alt_clim, kindex, l_kindex )

    p_lat_shift=r_lat_shift
    p_rdeltalat=r_rdeltalat
    CALL latitude_weights_li( &
        & jg,          jcs,            kproma,            kbdim,           &
        & krow,        wgt1_lat,       wgt2_lat,          inmw1_lat,       &
        & inmw2_lat,   p_lat_shift,    p_rdeltalat,       r_lat_clim,      &
        & lat_clim,    norder                                              )

    ! 2. Solar radiation
    ! 2.1 interpolate optical properties solar radiation
    DO jwl=1,nb_sw
      DO jk=1,klev
        DO jl=jcs,kproma
          idx_lat_1=inmw1_lat(jl)
          idx_lat_2=inmw2_lat(jl)
          w1_lat=wgt1_lat(jl)
          w2_lat=wgt2_lat(jl)
          idx_lev=kindex(jl,jk)
          zext_s(jl,jk,jwl)=tiw%weight1*(w1_lat*ext_v_s(jwl,idx_lev,idx_lat_1,nm1)+ &
                                        w2_lat*ext_v_s(jwl,idx_lev,idx_lat_2,nm1))+ &
                            tiw%weight2*(w1_lat*ext_v_s(jwl,idx_lev,idx_lat_1,nm2)+ &
                                        w2_lat*ext_v_s(jwl,idx_lev,idx_lat_2,nm2))
          zomg_s(jl,jk,jwl)=tiw%weight1*(w1_lat*ssa_v_s(jwl,idx_lev,idx_lat_1,nm1)+ &
                                        w2_lat*ssa_v_s(jwl,idx_lev,idx_lat_2,nm1))+ &
                            tiw%weight2*(w1_lat*ssa_v_s(jwl,idx_lev,idx_lat_1,nm2)+ &
                                        w2_lat*ssa_v_s(jwl,idx_lev,idx_lat_2,nm2))
          zasy_s(jl,jk,jwl)=tiw%weight1*(w1_lat*asy_v_s(jwl,idx_lev,idx_lat_1,nm1)+ &
                                        w2_lat*asy_v_s(jwl,idx_lev,idx_lat_2,nm1))+ &
                            tiw%weight2*(w1_lat*asy_v_s(jwl,idx_lev,idx_lat_1,nm2)+ &
                                        w2_lat*asy_v_s(jwl,idx_lev,idx_lat_2,nm2))
        END DO
      END DO
    END DO
    DO jwl=1,nb_sw
      DO jl=jcs,kproma
        idx_lat_1=inmw1_lat(jl)
        idx_lat_2=inmw2_lat(jl)
        w1_lat=wgt1_lat(jl)
        w2_lat=wgt2_lat(jl)
        zaod_s(jl,jwl)=tiw%weight1*(w1_lat*aod_v_s(jwl,idx_lat_1,nm1)+ &
                                   w2_lat*aod_v_s(jwl,idx_lat_2,nm1))+ &
                       tiw%weight2*(w1_lat*aod_v_s(jwl,idx_lat_1,nm2)+ &
                                   w2_lat*aod_v_s(jwl,idx_lat_2,nm2))
      END DO
    END DO

    ! 2.2 normalize zext to the correct total optical depth
    !     the normalization factor generally depends on the wavelength if
    !     the ratios of the extinction at different wavelengths are not
    !     independent of the height level. Generally, the aerosol composition
    !     depends on height, this leads to different ratios of the extinction
    !     between two given wavelengths at different heights.
    zext_s_int(jcs:kproma,1:nb_sw)=0._wp
    DO jwl=1,nb_sw
      DO jk=1,klev
        zext_s_int(jcs:kproma,jwl)=zext_s_int(jcs:kproma,jwl) + &
                 zext_s(jcs:kproma,jk,jwl)*dz(jcs:kproma,jk)
      END DO
    END DO
    WHERE (zext_s_int(jcs:kproma,1:nb_sw) > 0._wp)
      zfact_s(jcs:kproma,1:nb_sw)=zaod_s(jcs:kproma,1:nb_sw)/ &
                              zext_s_int(jcs:kproma,1:nb_sw)
    ELSEWHERE
      zfact_s(jcs:kproma,1:nb_sw)=1._wp
    END WHERE
    DO jwl=1,nb_sw
      DO jk=1,klev
        zext_s(jcs:kproma,jk,jwl)=zext_s(jcs:kproma,jk,jwl)* &
               dz(jcs:kproma,jk)*zfact_s(jcs:kproma,jwl)
      END DO
    END DO

    ! 2.3 add optical parameters to the optical parameters of aerosols
    !     inverse height profile
    DO jk=1,klev
      jki=klev-jk+1
      WHERE (zext_s(jcs:kproma,jki,1:nb_sw)>0._wp)
        paer_cg_sw_vr(jcs:kproma,jk,1:nb_sw)=paer_tau_sw_vr(jcs:kproma,jk,1:nb_sw)*&
          paer_piz_sw_vr(jcs:kproma,jk,1:nb_sw)*paer_cg_sw_vr(jcs:kproma,jk,1:nb_sw)+&
          zext_s(jcs:kproma,jki,1:nb_sw)*zomg_s(jcs:kproma,jki,1:nb_sw)*&
          zasy_s(jcs:kproma,jki,1:nb_sw)
        paer_piz_sw_vr(jcs:kproma,jk,1:nb_sw)=paer_tau_sw_vr(jcs:kproma,jk,1:nb_sw)*&
          paer_piz_sw_vr(jcs:kproma,jk,1:nb_sw)+&
          zext_s(jcs:kproma,jki,1:nb_sw)*zomg_s(jcs:kproma,jki,1:nb_sw)
        paer_tau_sw_vr(jcs:kproma,jk,1:nb_sw)=paer_tau_sw_vr(jcs:kproma,jk,1:nb_sw)+&
          zext_s(jcs:kproma,jki,1:nb_sw)
        paer_piz_sw_vr(jcs:kproma,jk,1:nb_sw)=paer_piz_sw_vr(jcs:kproma,jk,1:nb_sw)/&
          paer_tau_sw_vr(jcs:kproma,jk,1:nb_sw)
        paer_cg_sw_vr(jcs:kproma,jk,1:nb_sw)=paer_cg_sw_vr(jcs:kproma,jk,1:nb_sw)/&
          (paer_tau_sw_vr(jcs:kproma,jk,1:nb_sw)*paer_piz_sw_vr(jcs:kproma,jk,1:nb_sw))
      END WHERE
    END DO

    ! 3. far infrared
    ! 2.1 interpolate optical properties thermal radiation
    DO jwl=1,nb_lw
      DO jk=1,klev
        DO jl=jcs,kproma
          idx_lat_1=inmw1_lat(jl)
          idx_lat_2=inmw2_lat(jl)
          w1_lat=wgt1_lat(jl)
          w2_lat=wgt2_lat(jl)
          idx_lev=kindex(jl,jk)
          zext_t(jl,jk,jwl)=tiw%weight1*(w1_lat*ext_v_t(jwl,idx_lev,idx_lat_1,nm1)+ &
                                        w2_lat*ext_v_t(jwl,idx_lev,idx_lat_2,nm1))+ &
                            tiw%weight2*(w1_lat*ext_v_t(jwl,idx_lev,idx_lat_1,nm2)+ &
                                        w2_lat*ext_v_t(jwl,idx_lev,idx_lat_2,nm2))
          zomg_t(jl,jk,jwl)=tiw%weight1*(w1_lat*ssa_v_t(jwl,idx_lev,idx_lat_1,nm1)+ &
                                        w2_lat*ssa_v_t(jwl,idx_lev,idx_lat_2,nm1))+ &
                            tiw%weight2*(w1_lat*ssa_v_t(jwl,idx_lev,idx_lat_1,nm2)+ &
                                        w2_lat*ssa_v_t(jwl,idx_lev,idx_lat_2,nm2))
        END DO
      END DO
    END DO
    DO jwl=1,nb_lw
      DO jl=jcs,kproma
        idx_lat_1=inmw1_lat(jl)
        idx_lat_2=inmw2_lat(jl)
        w1_lat=wgt1_lat(jl)
        w2_lat=wgt2_lat(jl)
        zaod_t(jl,jwl)=tiw%weight1*(w1_lat*aod_v_t(jwl,idx_lat_1,nm1)+ &
                                   w2_lat*aod_v_t(jwl,idx_lat_2,nm1))+ &
                       tiw%weight2*(w1_lat*aod_v_t(jwl,idx_lat_1,nm2)+ &
                                   w2_lat*aod_v_t(jwl,idx_lat_2,nm2))
      END DO
    END DO

    ! 2.2 normalize zext to the correct total optical depth
    !     the normalization factor generally depends on the wavelength if
    !     the ratios of the extinction at different wavelengths are not
    !     independent of the height level. Generally, the aerosol composition
    !     depends on height, this leads to different ratios of the extinction
    !     between two given wavelengths at different heights.
    zext_t_int(jcs:kproma,1:nb_lw)=0._wp
    DO jwl=1,nb_lw
      DO jk=1,klev
        zext_t_int(jcs:kproma,jwl)=zext_t_int(jcs:kproma,jwl) + &
                zext_t(jcs:kproma,jk,jwl)*dz(jcs:kproma,jk)
      END DO
    END DO
    WHERE (zext_t_int(jcs:kproma,1:nb_lw) > 0._wp)
      zfact_t(jcs:kproma,1:nb_lw)=zaod_t(jcs:kproma,1:nb_lw)/ &
                              zext_t_int(jcs:kproma,1:nb_lw)
    ELSEWHERE
      zfact_t(jcs:kproma,1:nb_lw)=1._wp
    END WHERE
    DO jwl=1,nb_lw
      DO jk=1,klev
        zext_t(jcs:kproma,jk,jwl)=zext_t(jcs:kproma,jk,jwl)* &
              dz(jcs:kproma,jk)*zfact_t(jcs:kproma,jwl)
      END DO
    END DO

    ! 2.3 add optical parameters to the optical parameters of aerosols
    !     inverse height profile
    DO jk=1,klev
      jki=klev-jk+1
      ! use explicit DO loop to circumvent possible SXf90 compiler bug
      DO jc = jcs,kproma
        paer_tau_lw_vr(jc,jk,1:nb_lw)=paer_tau_lw_vr(jc,jk,1:nb_lw) + &
              zext_t(jc,jki,1:nb_lw)*(1._wp-zomg_t(jc,jki,1:nb_lw))
      ENDDO
    END DO

  END SUBROUTINE add_bc_aeropt_cmip6_volc

  !------------------------------------------------------------------------
  SUBROUTINE altitude_index ( &
        & jcs, kproma, zf, dz_clim, z_max_lim_clim, k_lev_clim, kindex, l_kindex &
      )

    INTEGER, INTENT(IN) :: jcs !< Minimum block index.
    INTEGER, INTENT(IN) :: kproma !< Maximum block index.
    REAL(wp), INTENT(IN) :: zf(:,:) !< Mid-layer altitudes of ICON grid (jc,jl).
    REAL(wp), INTENT(IN) :: dz_clim !< Layer thickness of climatology.
    REAL(wp), INTENT(IN) :: z_max_lim_clim !< Altitude of highest layer limit of clim.
    INTEGER, INTENT(IN)  :: k_lev_clim !< Number of layers in climatology.
    !> layer index of climatology in which icon layer is located (jc,jl).
    INTEGER, INTENT(OUT) :: kindex(:,:)
    !> `.TRUE.` if index in range `[1,k_lev_clim]` (jc,jl).
    LOGICAL, INTENT(OUT) :: l_kindex(:,:)

    REAL(wp) :: dz_clim_inv

    dz_clim_inv = 1._wp/dz_clim
    kindex(jcs:kproma,:) = FLOOR((z_max_lim_clim-zf(jcs:kproma,:))*dz_clim_inv)+1

    WHERE (kindex(jcs:kproma,:) < 1 .OR. kindex(jcs:kproma,:) > k_alt_clim)
      l_kindex(jcs:kproma,:) = .FALSE.
      kindex(jcs:kproma,:) = k_lev_clim+1
    ELSEWHERE
      l_kindex(jcs:kproma,:) = .TRUE.
    END WHERE
  END SUBROUTINE altitude_index

  !>
  !! Read the month range `imnthb:imonthe` for base year `iyear` into the global arrays.
  SUBROUTINE read_months_bc_aeropt_cmip6_volc (imnthb, imnthe, iyear)

    !> Begin month to read (may be `0` for month 12 of previous year).
    INTEGER, INTENT(IN) :: imnthb
    !> End month to read (may be `13` for month 1 of next year).
    INTEGER, INTENT(IN) :: imnthe
    !> Base year.
    INTEGER(i8), INTENT(IN) :: iyear

    INTEGER :: kmonthb, kmonthe
    ! Space for YYYYY.nc suffix
    CHARACTER(LEN=LEN(filename_base)+5+3) :: cfnameyear

    CHARACTER(len=*), PARAMETER :: subroutine_name = &
        & 'mo_bc_aeropt_cmip6_volc:read_months_bc_aeropt_cmip6_volc'

    IF (imnthb < 0 .OR. imnthe < imnthb .OR. imnthe > 13 ) THEN
      WRITE (message_text, '(a,2(a,i0))') &
          'months to be read outside valid range 0<=imnthb<=imnthe<=13, ', &
          'imnthb=', imnthb, ', imnthe=', imnthe
      CALL finish(subroutine_name, message_text)
    END IF

    WRITE (message_text,'(a,i2,a,i2)') ' reading CMIP6 volcanic aerosols from imonth ', imnthb, ' to ', imnthe
    CALL message(subroutine_name, message_text)

    ! Read data for last month of previous year

    IF (imnthb == 0) THEN

      WRITE (cfnameyear,'(a,i0,a)') filename_base, iyear-1, '.nc'

      CALL read_single_month_bc_aeropt_cmip6_volc ( &
          & filename=cfnameyear, &
          & aod_s=aod_v_s(:,:,0:0), &
          & ssa_s=ssa_v_s(:,:,:,0:0), &
          & asy_s=asy_v_s(:,:,:,0:0), &
          & ext_s=ext_v_s(:,:,:,0:0), &
          & aod_t=aod_v_t(:,:,0:0), &
          & ssa_t=ssa_v_t(:,:,:,0:0), &
          & ext_t=ext_v_t(:,:,:,0:0), &
          & start_timestep=12, &
          & end_timestep=12 &
        )

    END IF

    WRITE (cfnameyear,'(a,i0,a)') filename_base, iyear, '.nc'

    kmonthb=MAX(1,imnthb)
    kmonthe=MIN(12,imnthe)

    CALL read_single_month_bc_aeropt_cmip6_volc ( &
        & filename=cfnameyear, &
        & aod_s=aod_v_s(:,:,kmonthb:kmonthe), &
        & ssa_s=ssa_v_s(:,:,:,kmonthb:kmonthe), &
        & asy_s=asy_v_s(:,:,:,kmonthb:kmonthe), &
        & ext_s=ext_v_s(:,:,:,kmonthb:kmonthe), &
        & aod_t=aod_v_t(:,:,kmonthb:kmonthe), &
        & ssa_t=ssa_v_t(:,:,:,kmonthb:kmonthe), &
        & ext_t=ext_v_t(:,:,:,kmonthb:kmonthe), &
        & start_timestep=kmonthb, &
        & end_timestep=kmonthe &
      )

    ! Read data for first month of next year
    IF (imnthe == 13) THEN

      WRITE (cfnameyear,'(a,i0,a)') filename_base, iyear+1, '.nc'

      CALL read_single_month_bc_aeropt_cmip6_volc ( &
          & filename=cfnameyear, &
          & aod_s=aod_v_s(:,:,13:13), &
          & ssa_s=ssa_v_s(:,:,:,13:13), &
          & asy_s=asy_v_s(:,:,:,13:13), &
          & ext_s=ext_v_s(:,:,:,13:13), &
          & aod_t=aod_v_t(:,:,13:13), &
          & ssa_t=ssa_v_t(:,:,:,13:13), &
          & ext_t=ext_v_t(:,:,:,13:13), &
          & start_timestep=1, &
          & end_timestep=1 &
        )
    END IF

  END SUBROUTINE read_months_bc_aeropt_cmip6_volc

  !> Read a month range from a single file.
  SUBROUTINE read_single_month_bc_aeropt_cmip6_volc (&
        & filename, aod_s, ssa_s, asy_s, ext_s, aod_t, ssa_t, ext_t, start_timestep, end_timestep &
      )

    CHARACTER(len=*), INTENT(IN) :: filename
    REAL(wp), INTENT(INOUT) :: aod_s(:,:,:)
    REAL(wp), INTENT(INOUT) :: ssa_s(:,:,:,:)
    REAL(wp), INTENT(INOUT) :: asy_s(:,:,:,:)
    REAL(wp), INTENT(INOUT) :: ext_s(:,:,:,:)
    REAL(wp), INTENT(INOUT) :: aod_t(:,:,:)
    REAL(wp), INTENT(INOUT) :: ssa_t(:,:,:,:)
    REAL(wp), INTENT(INOUT) :: ext_t(:,:,:,:)
    INTEGER, INTENT(IN) :: start_timestep, end_timestep

    INTEGER :: file_id
    REAL(wp), POINTER :: var(:,:,:,:)
    REAL(wp) :: delta_alt


    CALL message ('mo_bc_aeropt_cmip6_volc:read_months_bc_aeropt_cmip6_volc', &
         &            ' reading from file '//TRIM(ADJUSTL(filename)))

    CALL openInputFile(file_id, filename)

    CALL read_extdim_slice_extdim_extdim_extdim ( &
        & file_id=file_id, &
        & variable_name='ext_sun', &
        & return_pointer=var, &
        & dim_names=dim_names_sol, &
        & start_extdim1=start_timestep, &
        & end_extdim1=end_timestep &
      )
    !convert units from 1/km to 1/m
    var(:,:,:,:) = 0.001_wp * var(:,:,:,:)
    CALL permute_shape(var, ext_s)
    DEALLOCATE(var)

    CALL read_extdim_slice_extdim_extdim_extdim ( &
        & file_id=file_id, &
        & variable_name='omega_sun', &
        & return_pointer=var, &
        & dim_names=dim_names_sol, &
        & start_extdim1=start_timestep, &
        & end_extdim1=end_timestep &
      )
    CALL permute_shape(var, ssa_s)
    DEALLOCATE(var)

    CALL read_extdim_slice_extdim_extdim_extdim ( &
        & file_id=file_id, &
        & variable_name='g_sun', &
        & return_pointer=var, &
        & dim_names=dim_names_sol, &
        & start_extdim1=start_timestep, &
        & end_extdim1=end_timestep &
      )
    CALL permute_shape(var, asy_s)
    DEALLOCATE(var)

    CALL read_extdim_slice_extdim_extdim_extdim ( &
        & file_id=file_id, &
        & variable_name='ext_earth', &
        & return_pointer=var, &
        & dim_names=dim_names_terr, &
        & start_extdim1=start_timestep, &
        & end_extdim1=end_timestep &
      )
    !convert units from 1/km to 1/m
    var(:,:,:,:) = 0.001_wp * var(:,:,:,:)
    CALL permute_shape(var, ext_t)
    DEALLOCATE(var)

    CALL read_extdim_slice_extdim_extdim_extdim ( &
        & file_id=file_id, &
        & variable_name='omega_earth', &
        & return_pointer=var, &
        & dim_names=dim_names_terr, &
        & start_extdim1=start_timestep, &
        & end_extdim1=end_timestep &
      )
    CALL permute_shape(var, ssa_t)
    DEALLOCATE(var)

    CALL closeFile(file_id)


    ! calculate AOD of atmosphere by numerically integrating extinction over altitude levels.
    delta_alt = r_alt_clim(1) - r_alt_clim(2)

    CALL trapezoidal_rule(delta_alt, ext_s(:,:,:,:), aod_s(:,:,:))
    CALL trapezoidal_rule(delta_alt, ext_t(:,:,:,:), aod_t(:,:,:))

  CONTAINS

    !> Permute the array dimensions to (bands,levels,latitudes,months), order altitude levels top
    !! to bottom and latitudes north to south, and fill the north and south pole points.
    SUBROUTINE permute_shape (src, tgt)

      REAL(wp), INTENT(IN) :: src(:,:,:,:)
      REAL(wp), INTENT(INOUT) :: tgt(:,:,:,:)

      INTEGER :: nbands
      INTEGER :: nlev
      INTEGER :: nlat
      INTEGER :: nmonths

      nmonths = SIZE(src,1)
      nlev = SIZE(src,2)
      nlat = SIZE(src,3)
      nbands = SIZE(src,4)

      tgt(1:nbands,1:nlev,2:nlat+1,1:nmonths) = RESHAPE( &
          & src(1:nmonths,nlev:1:-1,nlat:1:-1,1:nbands), &
          & SHAPE=[nbands,nlev,nlat,nmonths], &
          & ORDER=[4,2,3,1] &
        )

      tgt(1:nbands,1:nlev,1,:) = tgt(1:nbands,1:nlev,2,:)
      tgt(1:nbands,1:nlev,nlat+2,:) = tgt(1:nbands,1:nlev,nlat+1,:)

    END SUBROUTINE permute_shape

  END SUBROUTINE read_single_month_bc_aeropt_cmip6_volc

  !>
  !! Calculates trapezoidal sum for functions f(:,1:n,:,:) over the n
  !! grid points in the second dimension.
  SUBROUTINE trapezoidal_rule (dx, f, f_ts)
    !> Interval length of equidistant grid of the second dimension.
    REAL(wp), INTENT(IN) :: dx
    !> Tabulated functions f(:,1:n,:,:) on n equidistant grid points.
    REAL(wp), INTENT(IN) :: f(:,:,:,:)
    !> The integrals over f(:,1:n,:,:)
    REAL(wp), INTENT(OUT) :: f_ts(:,:,:)

    INTEGER :: n
    INTEGER :: lbnd, nm

    n = SIZE(f, DIM=2)

    DO nm = 1, SIZE(f, 4)
      DO lbnd = 1, SIZE(f, 1)
        f_ts(lbnd,:,nm) = dx * (SUM(f(lbnd,2:n-1,:,nm),DIM=1) + &
            & 0.5_wp * (f(lbnd,1,:,nm) + f(lbnd,n,:,nm)))
      END DO
    END DO
  END SUBROUTINE trapezoidal_rule

END MODULE mo_bc_aeropt_cmip6_volc
