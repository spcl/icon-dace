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

! Read and apply monthly aerosol optical properties of S. Kinne
! from yearly files.

MODULE mo_bc_aeropt_kinne

  USE mo_kind,                 ONLY: wp, i8
  USE mo_model_domain,         ONLY: t_patch
  USE mo_impl_constants,       ONLY: max_dom
  USE mo_grid_config,          ONLY: n_dom
  USE mo_parallel_config,      ONLY: nproma
  USE mo_exception,            ONLY: finish, message, message_text
  USE mo_io_config,            ONLY: default_read_method
  USE mo_time_config,          ONLY: time_config
  USE mo_read_interface,       ONLY: openInputFile, closeFile, on_cells, &
    &                                t_stream_id, read_0D_real, read_3D_time
  USE mtime,                   ONLY: datetime

  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights, &
    &                                  calculate_time_interpolation_weights
  USE mo_timer,                ONLY: ltimer, timer_start, timer_stop, &
    &                                timer_coupling
  USE mo_atmo_aero_provider_coupling, ONLY: couple_atmo_to_aero_provider

  IMPLICIT NONE

  PRIVATE
  PUBLIC                           :: read_bc_aeropt_kinne, &
    &                                 set_bc_aeropt_kinne

  TYPE t_ext_aeropt_kinne
     ! Fine mode SW
     REAL(wp), ALLOCATABLE :: aod_f_s(:,:,:,:)
     REAL(wp), ALLOCATABLE :: ssa_f_s(:,:,:,:)
     REAL(wp), ALLOCATABLE :: asy_f_s(:,:,:,:)
     ! Coarse mode SW
     REAL(wp), ALLOCATABLE :: aod_c_s(:,:,:,:)
     REAL(wp), ALLOCATABLE :: ssa_c_s(:,:,:,:)
     REAL(wp), ALLOCATABLE :: asy_c_s(:,:,:,:)
     ! Coarse mode LW
     REAL(wp), ALLOCATABLE :: aod_c_f(:,:,:,:)
     REAL(wp), ALLOCATABLE :: ssa_c_f(:,:,:,:)
     REAL(wp), ALLOCATABLE :: asy_c_f(:,:,:,:)
     ! Fine mode height profiles
     REAL(wp), ALLOCATABLE :: z_km_aer_f_mo(:,:,:,:)
     ! Coarse mode height profiles
     REAL(wp), ALLOCATABLE :: z_km_aer_c_mo(:,:,:,:)
  END TYPE t_ext_aeropt_kinne

  TYPE(t_ext_aeropt_kinne), ALLOCATABLE :: ext_aeropt_kinne(:)

  INTEGER(i8), SAVE                :: pre_year(max_dom)=-HUGE(1)
  LOGICAL, SAVE                    :: is_transient(max_dom) = .FALSE.
  INTEGER, PARAMETER               :: lev_clim=40
  REAL(wp)                         :: dz_clim
  REAL(wp)                         :: rdz_clim

  INTEGER                          :: nyears
  INTEGER                          :: imonth_beg, imonth_end

  LOGICAL                          :: lend_of_year

  TYPE(t_time_interpolation_weights) :: tiw_beg
  TYPE(t_time_interpolation_weights) :: tiw_end

CONTAINS
  !>
  !>
  !! SUBROUTINE su_bc_aeropt_kinne -- sets up the memory for fields in which
  !! the aerosol optical properties are stored when needed
SUBROUTINE su_bc_aeropt_kinne(p_patch, nbndlw, nbndsw, opt_from_coupler)

  TYPE(t_patch), INTENT(in)       :: p_patch
  INTEGER, INTENT(in)             :: nbndlw, nbndsw
  LOGICAL, INTENT(IN), OPTIONAL   :: opt_from_coupler

  INTEGER                         :: jg
  INTEGER                         :: nblks_len, nblks
  LOGICAL                         :: from_coupler = .FALSE.

  jg = p_patch%id

  nblks=p_patch%nblks_c
  nblks_len=nproma

  ! Check after merging icon-aes-link-echam-bc

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
    IF ( lend_of_year .OR. ( imonth_end == 1 .AND. time_config%tc_stopdate%date%month == 12 ) ) imonth_end = 13
  ENDIF

! on first call allocate structure for all grids
  IF ( jg == 1 ) ALLOCATE(ext_aeropt_kinne(n_dom))

  IF (PRESENT(opt_from_coupler)) from_coupler = opt_from_coupler

  IF ( from_coupler ) THEN

    ! time interpolation is done by coupler, thus we only need 1 timestep
    imonth_beg = 1; imonth_end = 1

    ! set vertical grid spacing
    dz_clim = 500.0

    WRITE(message_text,'(a,f6.2)') ' delta_z set to ', dz_clim
    CALL message('mo_bc_aeropt_kinne:read_months_bc_aeropt_kinne', message_text)

  ENDIF

  WRITE(message_text,'(a,i2,a,i2)') &
     & ' Allocating Kinne aerosols for months ', imonth_beg, ' to ', imonth_end
  CALL message('mo_bc_aeropt_kinne:su_bc_aeropt_kinne', message_text)

! allocate memory for optical properties on grid jg
  ALLOCATE(ext_aeropt_kinne(jg)% aod_c_s(nblks_len,nbndsw,nblks,imonth_beg:imonth_end))
  ALLOCATE(ext_aeropt_kinne(jg)% aod_f_s(nblks_len,nbndsw,nblks,imonth_beg:imonth_end))
  ALLOCATE(ext_aeropt_kinne(jg)% ssa_c_s(nblks_len,nbndsw,nblks,imonth_beg:imonth_end))
  ALLOCATE(ext_aeropt_kinne(jg)% ssa_f_s(nblks_len,nbndsw,nblks,imonth_beg:imonth_end))
  ALLOCATE(ext_aeropt_kinne(jg)% asy_c_s(nblks_len,nbndsw,nblks,imonth_beg:imonth_end))
  ALLOCATE(ext_aeropt_kinne(jg)% asy_f_s(nblks_len,nbndsw,nblks,imonth_beg:imonth_end))
  ALLOCATE(ext_aeropt_kinne(jg)% aod_c_f(nblks_len,nbndlw,nblks,imonth_beg:imonth_end))
  ALLOCATE(ext_aeropt_kinne(jg)% ssa_c_f(nblks_len,nbndlw,nblks,imonth_beg:imonth_end))
  ALLOCATE(ext_aeropt_kinne(jg)% asy_c_f(nblks_len,nbndlw,nblks,imonth_beg:imonth_end))
  ALLOCATE(ext_aeropt_kinne(jg)% z_km_aer_c_mo(nblks_len,lev_clim,nblks,imonth_beg:imonth_end))
  ALLOCATE(ext_aeropt_kinne(jg)% z_km_aer_f_mo(nblks_len,lev_clim,nblks,imonth_beg:imonth_end))
  !$ACC ENTER DATA CREATE(ext_aeropt_kinne(jg))

  !$ACC ENTER DATA CREATE(ext_aeropt_kinne(jg)%aod_c_s, ext_aeropt_kinne(jg)%aod_f_s) &
  !$ACC   CREATE(ext_aeropt_kinne(jg)%ssa_c_s, ext_aeropt_kinne(jg)%ssa_f_s) &
  !$ACC   CREATE(ext_aeropt_kinne(jg)%asy_c_s, ext_aeropt_kinne(jg)%asy_f_s) &
  !$ACC   CREATE(ext_aeropt_kinne(jg)%aod_c_f, ext_aeropt_kinne(jg)%ssa_c_f) &
  !$ACC   CREATE(ext_aeropt_kinne(jg)%asy_c_f, ext_aeropt_kinne(jg)%z_km_aer_c_mo) &
  !$ACC   CREATE(ext_aeropt_kinne(jg)%z_km_aer_f_mo)
! initialize with zero
  ext_aeropt_kinne(jg)% aod_c_s(:,:,:,:) = 0._wp
  ext_aeropt_kinne(jg)% aod_f_s(:,:,:,:) = 0._wp
  ext_aeropt_kinne(jg)% ssa_c_s(:,:,:,:) = 0._wp
  ext_aeropt_kinne(jg)% ssa_f_s(:,:,:,:) = 0._wp
  ext_aeropt_kinne(jg)% asy_c_s(:,:,:,:) = 0._wp
  ext_aeropt_kinne(jg)% asy_f_s(:,:,:,:) = 0._wp
  ext_aeropt_kinne(jg)% aod_c_f(:,:,:,:) = 0._wp
  ext_aeropt_kinne(jg)% ssa_c_f(:,:,:,:) = 0._wp
  ext_aeropt_kinne(jg)% asy_c_f(:,:,:,:) = 0._wp
  ext_aeropt_kinne(jg)% z_km_aer_c_mo(:,:,:,:) = 0._wp
  ext_aeropt_kinne(jg)% z_km_aer_f_mo(:,:,:,:) = 0._wp

  !$ACC UPDATE DEVICE(ext_aeropt_kinne(jg)%aod_c_s, ext_aeropt_kinne(jg)%aod_f_s) &
  !$ACC   DEVICE(ext_aeropt_kinne(jg)%ssa_c_s, ext_aeropt_kinne(jg)%ssa_f_s) &
  !$ACC   DEVICE(ext_aeropt_kinne(jg)%asy_c_s, ext_aeropt_kinne(jg)%asy_f_s) &
  !$ACC   DEVICE(ext_aeropt_kinne(jg)%aod_c_f, ext_aeropt_kinne(jg)%ssa_c_f) &
  !$ACC   DEVICE(ext_aeropt_kinne(jg)%asy_c_f, ext_aeropt_kinne(jg)%z_km_aer_c_mo) &
  !$ACC   DEVICE(ext_aeropt_kinne(jg)%z_km_aer_f_mo) &
  !$ACC   ASYNC(1)

END SUBROUTINE su_bc_aeropt_kinne

  !> SUBROUTINE shift_months_bc_aeropt_kinne -- shifts December of current year into imonth=0 and 
  !! January of the following year into imonth=1 (these months do not need to be read again.

SUBROUTINE shift_months_bc_aeropt_kinne(p_patch)

  TYPE(t_patch), INTENT(in)     :: p_patch

  INTEGER :: jg

  jg = p_patch%id

  IF ( .NOT. ALLOCATED(ext_aeropt_kinne) ) &
     &  CALL finish('mo_bc_aeropt_kinne:shift_months_bc_aeropt_kinne', &
     &              'ext_aeropt_kinne is not allocated')

  IF ( imonth_beg > 0 .OR. imonth_end < 13 ) THEN
     WRITE(message_text,'(a,i2,a,i2)') &
     & ' Kinne aerosols are allocated for months ', imonth_beg, ' to ', imonth_end, 'only.'
     CALL message('mo_bc_aeropt_kinne:shift_months_bc_aeropt_kinne', message_text)
     CALL finish('mo_bc_aeropt_kinne:shift_months_bc_aeropt_kinne', &
     & ' Kinne aerosols are not allocated over required range 0 to 13.')
  ENDIF

  WRITE(message_text,'(a)') &
     & ' Copy kinne aerosol for months 12:13 to months 0:1 '
  CALL message('mo_bc_aeropt_kinne:shift_months_bc_aeropt_kinne', message_text)

  ext_aeropt_kinne(jg)% aod_c_s(:,:,:,0:1) = ext_aeropt_kinne(jg)% aod_c_s(:,:,:,12:13)
  ext_aeropt_kinne(jg)% aod_f_s(:,:,:,0:1) = ext_aeropt_kinne(jg)% aod_f_s(:,:,:,12:13)
  ext_aeropt_kinne(jg)% ssa_c_s(:,:,:,0:1) = ext_aeropt_kinne(jg)% ssa_c_s(:,:,:,12:13)
  ext_aeropt_kinne(jg)% ssa_f_s(:,:,:,0:1) = ext_aeropt_kinne(jg)% ssa_f_s(:,:,:,12:13)
  ext_aeropt_kinne(jg)% asy_c_s(:,:,:,0:1) = ext_aeropt_kinne(jg)% asy_c_s(:,:,:,12:13)
  ext_aeropt_kinne(jg)% asy_f_s(:,:,:,0:1) = ext_aeropt_kinne(jg)% asy_f_s(:,:,:,12:13)
  ext_aeropt_kinne(jg)% aod_c_f(:,:,:,0:1) = ext_aeropt_kinne(jg)% aod_c_f(:,:,:,12:13)
  ext_aeropt_kinne(jg)% ssa_c_f(:,:,:,0:1) = ext_aeropt_kinne(jg)% ssa_c_f(:,:,:,12:13)
  ext_aeropt_kinne(jg)% asy_c_f(:,:,:,0:1) = ext_aeropt_kinne(jg)% asy_c_f(:,:,:,12:13)
  ext_aeropt_kinne(jg)% z_km_aer_c_mo(:,:,:,0:1) = ext_aeropt_kinne(jg)% z_km_aer_c_mo(:,:,:,12:13)
  ext_aeropt_kinne(jg)% z_km_aer_f_mo(:,:,:,0:1) = ext_aeropt_kinne(jg)% z_km_aer_f_mo(:,:,:,12:13)

END SUBROUTINE shift_months_bc_aeropt_kinne

  !> SUBROUTINE read_bc_aeropt_kinne -- read the aerosol optical properties 
  !! of the Kinne aerosols for the whole run at the beginning of the run
  !! before entering the time loop

SUBROUTINE read_bc_aeropt_kinne(mtime_current, p_patch, l_filename_year, nbndlw, nbndsw, opt_from_coupler)
  
  TYPE(datetime), POINTER, INTENT(in) :: mtime_current
  TYPE(t_patch), INTENT(in)           :: p_patch
  LOGICAL, INTENT(in)                 :: l_filename_year
  INTEGER, INTENT(in)                 :: nbndlw, nbndsw
  LOGICAL, OPTIONAL, INTENT(IN)       :: opt_from_coupler
 
  !LOCAL VARIABLES
  INTEGER(I8)                   :: iyear
  INTEGER                       :: imonthb, imonthe
  INTEGER                       :: jg
  LOGICAL                       :: from_coupler

  jg = p_patch%id

  from_coupler = .FALSE.

  IF (PRESENT(opt_from_coupler)) from_coupler=opt_from_coupler
  IF ( from_coupler) THEN

    IF ( pre_year(jg) == -HUGE(1) ) THEN
      CALL su_bc_aeropt_kinne(p_patch, nbndlw, nbndsw, opt_from_coupler=from_coupler)
      pre_year(jg)  =  mtime_current%date%year
    ENDIF

    IF (ltimer) CALL timer_start(timer_coupling)
    CALL couple_atmo_to_aero_provider( &
      p_patch, ext_aeropt_kinne(jg)%aod_f_s, ext_aeropt_kinne(jg)%ssa_f_s, &
      ext_aeropt_kinne(jg)%asy_f_s, ext_aeropt_kinne(jg)%aod_c_s, &
      ext_aeropt_kinne(jg)%ssa_c_s, ext_aeropt_kinne(jg)%asy_c_s, &
      ext_aeropt_kinne(jg)%aod_c_f, ext_aeropt_kinne(jg)%ssa_c_f, &
      ext_aeropt_kinne(jg)%z_km_aer_f_mo, ext_aeropt_kinne(jg)%z_km_aer_c_mo)
    IF (ltimer) CALL timer_stop(timer_coupling)

    !$ACC UPDATE DEVICE(ext_aeropt_kinne(jg)%aod_c_s, ext_aeropt_kinne(jg)%aod_f_s) &
    !$ACC   DEVICE(ext_aeropt_kinne(jg)%ssa_c_s, ext_aeropt_kinne(jg)%ssa_f_s) &
    !$ACC   DEVICE(ext_aeropt_kinne(jg)%asy_c_s, ext_aeropt_kinne(jg)%asy_f_s) &
    !$ACC   DEVICE(ext_aeropt_kinne(jg)%aod_c_f, ext_aeropt_kinne(jg)%ssa_c_f) &
    !$ACC   DEVICE(ext_aeropt_kinne(jg)%asy_c_f, ext_aeropt_kinne(jg)%z_km_aer_c_mo) &
    !$ACC   DEVICE(ext_aeropt_kinne(jg)%z_km_aer_f_mo) &
    !$ACC   ASYNC(1)

    RETURN

  END IF

  iyear = mtime_current%date%year

  IF (iyear > pre_year(jg)) THEN

    ! beginning of job or change of year

    IF ( pre_year(jg) > -HUGE(1) ) THEN
      CALL shift_months_bc_aeropt_kinne(p_patch)
    ELSE
      CALL su_bc_aeropt_kinne(p_patch, nbndlw, nbndsw)
    ENDIF

    ! Restrict reading of data to those months that are needed

    IF ( nyears > 1 ) THEN

      IF ( pre_year(jg) > -HUGE(1) ) THEN
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
              pre_year(jg) = mtime_current%date%year
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
      IF ( lend_of_year .OR. ( imonthe == 1 .AND. time_config%tc_stopdate%date%month == 12 ) ) imonthe = 13

    ENDIF

    CALL read_months_bc_aeropt_kinne ( &
                     'aod',            'ssa',    'asy',                        'z_aer_coarse_mo',  &
                     ext_aeropt_kinne(jg)%aod_c_s, ext_aeropt_kinne(jg)%ssa_c_s,                   &
                     ext_aeropt_kinne(jg)%asy_c_s, ext_aeropt_kinne(jg)%z_km_aer_c_mo,             &
                     'delta_z',        'lnwl',   'lev',                        imonthb,            &
                     imonthe,          iyear,     'bc_aeropt_kinne_sw_b14_coa', p_patch,           &
                     l_filename_year                                                               )
    ! for the coarse mode, the altitude distribution is wavelength independent and
    ! therefore for solar and long wave spectrum the same
    CALL read_months_bc_aeropt_kinne ( &
                     'aod',            'ssa',    'asy',                        'z_aer_coarse_mo',  &
                     ext_aeropt_kinne(jg)%aod_c_f, ext_aeropt_kinne(jg)% ssa_c_f,                  &
                     ext_aeropt_kinne(jg)%asy_c_f, ext_aeropt_kinne(jg)% z_km_aer_c_mo,            &
                     'delta_z',        'lnwl',   'lev',                        imonthb,            &
                     imonthe,          iyear,    'bc_aeropt_kinne_lw_b16_coa', p_patch,            &
                     l_filename_year                                                               )
    CALL read_months_bc_aeropt_kinne ( &
                     'aod',            'ssa',    'asy',                        'z_aer_fine_mo',    &
                     ext_aeropt_kinne(jg)%aod_f_s, ext_aeropt_kinne(jg)%ssa_f_s,                   &
                     ext_aeropt_kinne(jg)%asy_f_s, ext_aeropt_kinne(jg)%z_km_aer_f_mo,             &
                     'delta_z',        'lnwl',   'lev',                        imonthb,            &
                     imonthe,          iyear,    'bc_aeropt_kinne_sw_b14_fin', p_patch,            &
                     l_filename_year                                                               )

    rdz_clim = 1._wp/dz_clim
    pre_year(jg) = mtime_current%date%year
    is_transient(jg) = l_filename_year

    !$ACC UPDATE DEVICE(ext_aeropt_kinne(jg)%aod_c_s, ext_aeropt_kinne(jg)%aod_f_s) &
    !$ACC   DEVICE(ext_aeropt_kinne(jg)%ssa_c_s, ext_aeropt_kinne(jg)%ssa_f_s) &
    !$ACC   DEVICE(ext_aeropt_kinne(jg)%asy_c_s, ext_aeropt_kinne(jg)%asy_f_s) &
    !$ACC   DEVICE(ext_aeropt_kinne(jg)%aod_c_f, ext_aeropt_kinne(jg)%ssa_c_f) &
    !$ACC   DEVICE(ext_aeropt_kinne(jg)%asy_c_f, ext_aeropt_kinne(jg)%z_km_aer_c_mo) &
    !$ACC   DEVICE(ext_aeropt_kinne(jg)%z_km_aer_f_mo) &
    !$ACC   ASYNC(1)

  END IF    
END SUBROUTINE read_bc_aeropt_kinne
!-------------------------------------------------------------------------
!> SUBROUTINE set_bc_aeropt_kinne
!! set aerosol optical properties for all wave length bands (solar and IR)
!! in the case of the climatology of optical properties compiled by S.Kinne.
!! The height profile is taken into account.
SUBROUTINE set_bc_aeropt_kinne (    current_date,                         &
          & jg,                                                           &
          & jcs,                    jce,                nproma,           &
          & klev,                   jb,                                   &
          & nb_sw,                  nb_lw,                                &
          & zf,                     dz,                                   &
          & paer_tau_sw_vr,         paer_piz_sw_vr,     paer_cg_sw_vr,    &
          & paer_tau_lw_vr,                                               & 
          & opt_use_acc,            opt_from_coupler                      )

  ! !INPUT PARAMETERS

  TYPE(datetime), POINTER, INTENT(in) :: current_date
  INTEGER,INTENT(in)  :: jg,     &! grid index
                         jcs,    &! actual block, start at column with index jcs
                         jce,    &! actual block, end at column with index jce 
                         nproma, &! maximum block length
                         klev,   &! number of vertical levels
                         jb,     &! block index
                         nb_sw,  &! number of wave length bands (solar)
                         nb_lw    ! number of wave length bands (far IR)
  REAL(wp),INTENT(in) :: zf(nproma,klev)  ,& ! geometric height at full level [m]
                         dz(nproma,klev)     ! geometric height thickness     [m]
! !OUTPUT PARAMETERS
  REAL(wp),INTENT(out),DIMENSION(nproma,klev,nb_sw):: &
   paer_tau_sw_vr,   & !aerosol optical depth (solar), sum_i(tau_i)
   paer_piz_sw_vr,   & !weighted sum of single scattering albedos, 
                       !sum_i(tau_i*omega_i)
   paer_cg_sw_vr       !weighted sum of asymmetry factors, 
                       !sum_i(tau_i*omega_i*g_i)
  REAL(wp),INTENT(out),DIMENSION(nproma,klev,nb_lw):: &
   paer_tau_lw_vr      !aerosol optical depth (far IR)
  LOGICAL, INTENT(IN), OPTIONAL                          :: opt_use_acc
  LOGICAL, INTENT(IN), OPTIONAL                          :: opt_from_coupler

! !LOCAL VARIABLES
  
  INTEGER                           :: jl,jk,jwl
  REAL(wp), DIMENSION(nproma,klev)   :: zh_vr, &
                                       zdeltag_vr
  REAL(wp), DIMENSION(nproma)        :: zq_int ! integral height profile
  REAL(wp), DIMENSION(nproma,nb_lw)  :: zs_i
  REAL(wp), DIMENSION(nproma,nb_sw)  :: zt_c, zt_f, &
                                       zs_c, zs_f, &
                                       zg_c, zg_f, & ! time interpolated
                                       ! aod, ssa ,ssa*asy 
                                       ! (coarse (c), fine natural (n), 
                                       !  fine anthropogenic (a))
                                       ztaua_c,ztaua_f ! optical depths
                                       ! at various altitudes
  REAL(wp), DIMENSION(nproma,klev)   :: zq_aod_c, zq_aod_f ! altitude profile
                                       ! on echam grid (coarse and fine mode)
  INTEGER                           :: kindex ! index field
  TYPE(t_time_interpolation_weights) :: tiw
  LOGICAL :: use_acc  = .FALSE.        ! Default: no acceleration
  LOGICAL :: from_coupler = .FALSE.    ! Default: aerosol from interpolated files

  IF (PRESENT(opt_use_acc)) use_acc = opt_use_acc
  IF (PRESENT(opt_from_coupler)) from_coupler = opt_from_coupler

  tiw = calculate_time_interpolation_weights(current_date)

  IF (is_transient(jg) .AND. current_date%date%year /= pre_year(jg)) THEN
    WRITE (message_text,'(A,I4,A,I4)') 'Stale data: requested year is', current_date%date%year, &
        & ' but data is for ', pre_year(jg)
    CALL finish('mo_bc_aeropt_kinne:set_bc_aeropt_kinne', message_text)
  END IF

  !$ACC DATA PRESENT(zf, dz, paer_tau_sw_vr, paer_piz_sw_vr, paer_cg_sw_vr) &
  !$ACC   PRESENT(paer_tau_lw_vr, ext_aeropt_kinne) &
  !$ACC   CREATE(zh_vr, zdeltag_vr, zq_int, zs_i, zt_c, zt_f, zs_c, zs_f) &
  !$ACC   CREATE(zg_c, zg_f, ztaua_c, ztaua_f, zq_aod_c, zq_aod_f) &
  !$ACC   COPYIN(tiw) IF(use_acc)

! (i) calculate altitude above NN and layer thickness in 
!     echam for altitude profiles
  !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(use_acc)
  DO jk=1,klev
     DO jl=jcs,jce
        zdeltag_vr(jl,jk)=dz(jl,klev-jk+1)
        zh_vr(jl,jk)=zf(jl,klev-jk+1)
     END DO
  END DO

! (ii) calculate height profiles on echam grid for coarse and fine mode
  !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(use_acc)
  zq_aod_f(jcs:jce,1:klev)=0._wp
  zq_aod_c(jcs:jce,1:klev)=0._wp
  !$ACC END KERNELS

  IF ( from_coupler ) THEN
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(use_acc)
    DO jk=1,klev
      DO jl=jcs,jce
        kindex = MAX(INT(zh_vr(jl,jk)*rdz_clim+0.5_wp),1)
        IF (kindex > 0 .and. kindex <= lev_clim ) THEN
          zq_aod_c(jl,jk)= &
            & ext_aeropt_kinne(jg)% z_km_aer_c_mo(jl,kindex,jb,1)
          zq_aod_f(jl,jk)= &
            & ext_aeropt_kinne(jg)% z_km_aer_f_mo(jl,kindex,jb,1)
        END IF
      END DO
    END DO
  ELSE
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(use_acc)
    DO jk=1,klev
      DO jl=jcs,jce
        kindex = MAX(INT(zh_vr(jl,jk)*rdz_clim+0.5_wp),1)
        IF (kindex > 0 .and. kindex <= lev_clim ) THEN
          zq_aod_c(jl,jk)= &
            & ext_aeropt_kinne(jg)% z_km_aer_c_mo(jl,kindex,jb,tiw%month1_index)*tiw%weight1+ &
            & ext_aeropt_kinne(jg)% z_km_aer_c_mo(jl,kindex,jb,tiw%month2_index)*tiw%weight2
          zq_aod_f(jl,jk)= &
            & ext_aeropt_kinne(jg)% z_km_aer_f_mo(jl,kindex,jb,tiw%month1_index)*tiw%weight1+ &
            & ext_aeropt_kinne(jg)% z_km_aer_f_mo(jl,kindex,jb,tiw%month2_index)*tiw%weight2
        END IF
      END DO
    END DO
  END IF
  
! normalize height profile for coarse mode
  !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(use_acc)
  zq_int(jcs:jce)=0._wp
  !$ACC END KERNELS

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(use_acc)
  !$ACC LOOP SEQ
  DO jk=1,klev
     !$ACC LOOP GANG VECTOR
     DO jl=jcs,jce
        zq_int(jl)=zq_int(jl)+ &
                       & zq_aod_c(jl,jk)*zdeltag_vr(jl,jk)
     END DO
  END DO
  !$ACC END PARALLEL

  !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) IF(use_acc)
  DO jl=jcs,jce
     IF (zq_int(jl) <= 0._wp) zq_int(jl)=1._wp
  END DO

  !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(use_acc)
  DO jk=1,klev
     DO jl=jcs,jce
        zq_aod_c(jl,jk)=zdeltag_vr(jl,jk)*zq_aod_c(jl,jk) / zq_int(jl)
     END DO
  END DO

! normalize height profile for fine mode
  !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(use_acc)
  zq_int(jcs:jce)=0._wp
  !$ACC END KERNELS

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(use_acc)
  !$ACC LOOP SEQ
  DO jk=1,klev
    !$ACC LOOP GANG VECTOR
    DO jl=jcs,jce
       zq_int(jl)=zq_int(jl) + zq_aod_f(jl,jk)*zdeltag_vr(jl,jk)
    END DO
  END DO
  !$ACC END PARALLEL

  !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) IF(use_acc)
  DO jl=jcs,jce
    IF (zq_int(jl) <= 0._wp) zq_int(jl)=1._wp
  END DO

  !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(use_acc)
  DO jk=1,klev
    DO jl=jcs,jce
      zq_aod_f(jl,jk)=zdeltag_vr(jl,jk)*zq_aod_f(jl,jk)/zq_int(jl)
    END DO
  END DO

! (iii) far infrared
!  !$ACC KERNELS DEFAULT(PRESENT) COPYIN(tiw) ASYNC(1) IF(use_acc)
!  IF ( from_coupler ) THEN
!     zs_i(jcs:jce,1:nb_lw)=1._wp-ext_aeropt_kinne(jg)%ssa_c_f(jcs:jce,1:nb_lw,jb,1)
!  ELSE
!     zs_i(jcs:jce,1:nb_lw)=1._wp-(tiw%weight1*ext_aeropt_kinne(jg)% ssa_c_f(jcs:jce,1:nb_lw,jb,tiw%month1_index)+ &
!                                     tiw%weight2*ext_aeropt_kinne(jg)% ssa_c_f(jcs:jce,1:nb_lw,jb,tiw%month2_index))
!  END IF
!  !$ACC END KERNELS

  IF ( from_coupler ) THEN
     !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(use_acc)
      DO jwl=1,nb_lw
         DO jl=jcs,jce
            zs_i(jl,jwl)=1._wp-ext_aeropt_kinne(jg)%ssa_c_f(jl, jwl, jb, 1)
         END DO
      END DO
  ELSE
     !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(use_acc)
      DO jwl=1,nb_lw
         DO jl=jcs,jce
            zs_i(jl,jwl)=1._wp-(tiw%weight1*ext_aeropt_kinne(jg)% ssa_c_f(jl, jwl, jb, tiw%month1_index)+ &
                                tiw%weight2*ext_aeropt_kinne(jg)% ssa_c_f(jl, jwl, jb, tiw%month2_index))
         END DO
      END DO
  END IF

  !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1) IF(use_acc)
  DO jk=1,klev
     DO jwl=1,nb_lw
        DO jl=jcs,jce
           IF ( from_coupler ) THEN
              paer_tau_lw_vr(jl,jk,jwl)=zq_aod_c(jl,jk) * &
                   zs_i(jl,jwl) * &
                   ext_aeropt_kinne(jg)% aod_c_f(jl,jwl,jb,1)
           ELSE
              paer_tau_lw_vr(jl,jk,jwl)=zq_aod_c(jl,jk) * &
                    zs_i(jl,jwl) * &
                    (tiw%weight1*ext_aeropt_kinne(jg)% aod_c_f(jl,jwl,jb,tiw%month1_index) + &
                    tiw%weight2*ext_aeropt_kinne(jg)% aod_c_f(jl,jwl,jb,tiw%month2_index))
           END IF
        END DO
     END DO
  END DO

! (iv) solar radiation
! time interpolated single scattering albedo (omega_f, omega_c)
  IF ( from_coupler ) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(use_acc)
      DO jwl=1,nb_sw
         DO jl=jcs,jce
            zs_c(jl, jwl) = ext_aeropt_kinne(jg)% ssa_c_s(jl, jwl,jb,1)
            zs_f(jl, jwl) = ext_aeropt_kinne(jg)% ssa_f_s(jl, jwl,jb,1)

            zg_c(jl, jwl) = ext_aeropt_kinne(jg)% asy_c_s(jl, jwl,jb,1)
            zg_f(jl, jwl) = ext_aeropt_kinne(jg)% asy_f_s(jl, jwl,jb,1)

            zt_c(jl, jwl) = ext_aeropt_kinne(jg)% aod_c_s(jl, jwl,jb,1)
            zt_f(jl, jwl) = ext_aeropt_kinne(jg)% aod_f_s(jl, jwl,jb,1)
         END DO
      END DO
      !$ACC END PARALLEL LOOP
  ELSE
     !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1) IF(use_acc)
      DO jwl=1,nb_sw
         DO jl=jcs,jce
            zs_c(jl, jwl) = tiw%weight1*ext_aeropt_kinne(jg)% ssa_c_s(jl, jwl,jb,tiw%month1_index) + &
                            tiw%weight2*ext_aeropt_kinne(jg)% ssa_c_s(jl, jwl,jb,tiw%month2_index)
            zs_f(jl, jwl) = tiw%weight1*ext_aeropt_kinne(jg)% ssa_f_s(jl, jwl,jb,tiw%month1_index) + &
                            tiw%weight2*ext_aeropt_kinne(jg)% ssa_f_s(jl, jwl,jb,tiw%month2_index)

            zg_c(jl, jwl) = tiw%weight1*ext_aeropt_kinne(jg)% asy_c_s(jl, jwl,jb,tiw%month1_index) + &
                            tiw%weight2*ext_aeropt_kinne(jg)% asy_c_s(jl, jwl,jb,tiw%month2_index)
            zg_f(jl, jwl) = tiw%weight1*ext_aeropt_kinne(jg)% asy_f_s(jl, jwl,jb,tiw%month1_index) + &
                            tiw%weight2*ext_aeropt_kinne(jg)% asy_f_s(jl, jwl,jb,tiw%month2_index)

            zt_c(jl, jwl) = tiw%weight1*ext_aeropt_kinne(jg)% aod_c_s(jl, jwl,jb,tiw%month1_index) + &
                            tiw%weight2*ext_aeropt_kinne(jg)% aod_c_s(jl, jwl,jb,tiw%month2_index)
            zt_f(jl, jwl) = tiw%weight1*ext_aeropt_kinne(jg)% aod_f_s(jl, jwl,jb,tiw%month1_index) + &
                            tiw%weight2*ext_aeropt_kinne(jg)% aod_f_s(jl, jwl,jb,tiw%month2_index)
         END DO
      END DO
      !$ACC END PARALLEL LOOP
  END IF
 
! height interpolation
! calculate optical properties
  !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(3) ASYNC(1) IF(use_acc)
  DO jk=1,klev
     DO jwl=1,nb_sw
        DO jl=jcs,jce
           ! aerosol optical depth 
           ztaua_c(jl,jwl) = zt_c(jl,jwl)*zq_aod_c(jl,jk)
           ztaua_f(jl,jwl) = zt_f(jl,jwl)*zq_aod_f(jl,jk)
           paer_tau_sw_vr(jl,jk,jwl) = &
                         & ztaua_c(jl,jwl) + &
                         & ztaua_f(jl,jwl) 
           paer_piz_sw_vr(jl,jk,jwl) = &
                      & ztaua_c(jl,jwl)*zs_c(jl,jwl) + &
                      & ztaua_f(jl,jwl)*zs_f(jl,jwl)
           IF (paer_tau_sw_vr(jl,jk,jwl) /= 0._wp) THEN
              paer_piz_sw_vr(jl,jk,jwl) = paer_piz_sw_vr(jl,jk,jwl) / &
                                              & paer_tau_sw_vr(jl,jk,jwl)
           ELSE
              paer_piz_sw_vr(jl,jk,jwl) = 1._wp
           END IF
           paer_cg_sw_vr(jl,jk,jwl)  = &
                         & ztaua_c(jl,jwl)*zs_c(jl,jwl)*zg_c(jl,jwl) + &
                         & ztaua_f(jl,jwl)*zs_f(jl,jwl)*zg_f(jl,jwl)
           IF (paer_tau_sw_vr(jl,jk,jwl) /= 0._wp) THEN
              paer_cg_sw_vr (jl,jk,jwl) = paer_cg_sw_vr (jl,jk,jwl) / &
                                              & paer_piz_sw_vr(jl,jk,jwl) / &
                                              & paer_tau_sw_vr(jl,jk,jwl)
           ELSE
              paer_cg_sw_vr(jl,jk,jwl) = 0._wp
           END IF
        END DO
     END DO
  END DO
  !$ACC WAIT(1)
  !$ACC END DATA
  
END SUBROUTINE set_bc_aeropt_kinne
!-------------------------------------------------------------------------
! 
!> SUBROUTINE read_months_bc_aeropt_kinne -- reads optical aerosol parameters from file containing
!! aod, ssa, asy, aer_ex (altitude dependent extinction), dz_clim (layer 
!! thickness in meters), lev_clim (number of levels), and (optional) surface 
!! altitude in meters.
!!
SUBROUTINE read_months_bc_aeropt_kinne (                                   &
  caod,             cssa,             casy,               caer_ex,         &
  zaod,             zssa,             zasy,               zaer_ex,         &
  cdz_clim,         cwldim,           clevdim,            imnthb,          &
  imnthe,           iyear,            cfname,             p_patch,         &
  l_filename_year                                                          )
!
  CHARACTER(len=*), INTENT(in)   :: caod,    & ! name of variable containing optical depth of column
                                    cssa,    & ! name of variable containing single scattering albedo 
                                    casy,    & ! name of variable containing asymmetry factor
                                               ! ssa and asy are assumed to be constant over column
                                    caer_ex, & ! name of variable containing altitude dependent extinction
                                               ! aer_ex is normed to 1 (total over column is equal to 1)
                                    cdz_clim,& ! layer thickness of climatology in meters
                                    cwldim,  & ! name of wavelength dimension
                                    clevdim    ! name of level dimension in climatology

  INTEGER, INTENT(in)            :: imnthb,  & ! begin and ...
                                    imnthe     ! ... end month to be read

  INTEGER(i8), INTENT(in)        :: iyear      ! base year. if month=0, month 12 of previous year is read,
                                               ! if month=13, month 1 of subsequent year is read
  CHARACTER(len=*), INTENT(in)   :: cfname     ! file name containing variables

  TYPE(t_patch), INTENT(in)      :: p_patch
  LOGICAL, INTENT(in)            :: l_filename_year

  INTEGER                        :: ifile_id, kmonthb, kmonthe, ilen_cfname
  REAL(wp), INTENT(inout)        :: zaod(:,:,:,imonth_beg:)    ! has to be inout, otherwise
  REAL(wp), INTENT(inout)        :: zssa(:,:,:,imonth_beg:)    ! the NAG compiler will
  REAL(wp), INTENT(inout)        :: zasy(:,:,:,imonth_beg:)    ! create NaN when running
  REAL(wp), INTENT(inout)        :: zaer_ex(:,:,:,imonth_beg:) ! over the turn of the year.
  ! optional space for _DOM99 suffix
  CHARACTER(LEN=LEN(cfname)+6)   :: cfname2
  ! optional space for _YYYY.nc suffix
  CHARACTER(LEN=LEN(cfname)+6+12+4) :: cfnameyear
  INTEGER :: cfname2_tlen

  INTEGER                        :: jg

  jg = p_patch%id

  IF (imnthb < 0 .OR. imnthe < imnthb .OR. imnthe > 13 ) THEN
    WRITE (message_text, '(a,2(a,i0))') &
         'months to be read outside valid range 0<=imnthb<=imnthe<=13, ', &
         'imnthb=', imnthb, ', imnthe=', imnthe
    CALL finish('read_months_bc_aeropt_kinne in mo_bc_aeropt_kinne', &
      &         message_text)
  END IF
  ilen_cfname=LEN_TRIM(cfname)

  ! Add domain index if more than 1 grid is used
  IF (n_dom > 1) THEN
    WRITE(cfname2,'(a,a,i2.2)') cfname,'_DOM',jg
    cfname2_tlen = LEN_TRIM(cfname2)
  ELSE
    cfname2=cfname
    cfname2_tlen = ilen_cfname
  END IF

  WRITE(message_text,'(a,i2,a,i2)') ' reading Kinne aerosols from imonth ', imnthb, ' to ', imnthe
  CALL message('mo_bc_aeropt_kinne:read_months_bc_aeropt_kinne', message_text)

  ! Read data for last month of previous year

  IF (imnthb == 0) THEN

    IF (cfname(1:ilen_cfname) == 'bc_aeropt_kinne_sw_b14_fin' .AND. &
        l_filename_year ) THEN
      WRITE(cfnameyear,'(2a,i0,a)') cfname2(1:cfname2_tlen), '_', iyear-1, '.nc'
    ELSE
      cfnameyear=cfname2(1:cfname2_tlen)//'.nc'
    ENDIF

    CALL read_single_month_bc_aeropt_kinne(cfnameyear, &
         p_patch, caod, cssa, casy, caer_ex, &
         zaod=zaod(:,:,:,0:0), zssa=zssa(:,:,:,0:0), &
         zasy=zasy(:,:,:,0:0), zaer_ex=zaer_ex(:,:,:,0:0), &
         start_timestep=12, end_timestep=12, &
         cwldim=cwldim, clevdim=clevdim)
  END IF

  ! Read data for current year
  IF (cfname(1:ilen_cfname) == 'bc_aeropt_kinne_sw_b14_fin' .AND. &
      l_filename_year ) THEN
    WRITE(cfnameyear,'(2a,i0,a)') cfname2(1:cfname2_tlen), '_', iyear, '.nc'
  ELSE
    cfnameyear=TRIM(cfname2)//'.nc'
  ENDIF

  kmonthb=MAX(1,imnthb)
  kmonthe=MIN(12,imnthe)

  CALL read_single_month_bc_aeropt_kinne(cfnameyear, &
       p_patch, caod, cssa, casy, caer_ex, &
       zaod=zaod(:,:,:,kmonthb:kmonthe), zssa=zssa(:,:,:,kmonthb:kmonthe), &
       zasy=zasy(:,:,:,kmonthb:kmonthe), zaer_ex=zaer_ex(:,:,:,kmonthb:kmonthe), &
       start_timestep=kmonthb, end_timestep=kmonthe, &
       cwldim=cwldim, clevdim=clevdim)

  ! Read data for first month of next year
  IF (imnthe == 13) THEN

    IF (cfname(1:ilen_cfname) == 'bc_aeropt_kinne_sw_b14_fin' .AND. &
        l_filename_year ) THEN
      WRITE(cfnameyear,'(2a,i0,a)') cfname2(1:cfname2_tlen), '_', iyear+1, '.nc'
    ELSE
      cfnameyear=cfname2(1:cfname2_tlen)//'.nc'
    ENDIF

    CALL read_single_month_bc_aeropt_kinne(cfnameyear, &
         p_patch, caod, cssa, casy, caer_ex, &
         zaod=zaod(:,:,:,13:13), zssa=zssa(:,:,:,13:13), &
         zasy=zasy(:,:,:,13:13), zaer_ex=zaer_ex(:,:,:,13:13), &
         start_timestep=1, end_timestep=1, &
         cwldim=cwldim, clevdim=clevdim)
  END IF

  ! we assume here that delta_z (aka dz_clim) does not vary over the files
  ! thus it does not matter from which file we get these values:

  CALL openInputFile(ifile_id, cfnameyear)
  dz_clim = read_0D_real (file_id=ifile_id, variable_name=cdz_clim)
  CALL closeFile(ifile_id)

END SUBROUTINE read_months_bc_aeropt_kinne

  SUBROUTINE read_single_month_bc_aeropt_kinne(cfnameyear, p_patch, &
       caod, cssa, casy, caer_ex, zaod, zssa, zasy, zaer_ex, &
       start_timestep, end_timestep, cwldim, clevdim)
    CHARACTER(len=*), INTENT(in) :: cfnameyear, caod, cssa, casy, caer_ex, &
         cwldim, clevdim
    TYPE(t_patch), INTENT(in)    :: p_patch
    INTEGER, INTENT(in)          :: start_timestep, end_timestep
    REAL(wp), INTENT(inout)      :: zaod(:,:,:,:)    ! has to be inout, otherwise
    REAL(wp), INTENT(inout)      :: zssa(:,:,:,:)    ! the NAG compiler will
    REAL(wp), INTENT(inout)      :: zasy(:,:,:,:)    ! create NaN when running
    REAL(wp), INTENT(inout)      :: zaer_ex(:,:,:,:) ! over the turn of the year.
    TYPE(t_stream_id)            :: stream_id

    CALL message ('mo_bc_aeropt_kinne:read_months_bc_aeropt_kinne', &
         &            ' reading from file '//TRIM(ADJUSTL(cfnameyear)))
   
    CALL openInputFile(stream_id, cfnameyear, p_patch, default_read_method)

    CALL read_3D_time(stream_id=stream_id, location=on_cells, &
           &          variable_name=caod, fill_array=zaod, &
           &          start_timestep=start_timestep, end_timestep=end_timestep, &
           &          levelsDimName=cwldim)
    CALL read_3D_time(stream_id=stream_id, location=on_cells, &
           &          variable_name=cssa, fill_array=zssa, &
           &          start_timestep=start_timestep, end_timestep=end_timestep, &
           &          levelsDimName=cwldim)
    CALL read_3D_time(stream_id=stream_id, location=on_cells, &
           &          variable_name=casy, fill_array=zasy, &
           &          start_timestep=start_timestep, end_timestep=end_timestep, &
           &          levelsDimName=cwldim)
    CALL read_3D_time(stream_id=stream_id, location=on_cells, &
           &          variable_name=caer_ex, fill_array=zaer_ex, &
           &          start_timestep=start_timestep, end_timestep=end_timestep, &
           &          levelsDimName=clevdim)
    CALL closeFile(stream_id)

  END SUBROUTINE read_single_month_bc_aeropt_kinne
END MODULE mo_bc_aeropt_kinne
