!+ Module to provide vct_a for WW
!
MODULE mo_wwonly
!
! Description:
! <Say what this module is for>
!
! Current Code Owner: DWD, <Name of person responsible for this code>
!    <smail, phone, fax and email>
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! @VERSION@    @DATE@     <Your name>
!  Initial release
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: wp, vct_a, vct_a_from_hhl

  PUBLIC :: rdv, O_m_rdv, rd_o_cpd, tmelt, alvdcp, b3, &
            b1, b2w, b4w,                              &
            max_dom,                                   &
            kstart_moist

  PUBLIC :: sat_pres_water

!!PUBLIC :: t, qv, qc, u, v, p, hhl,                    &
!!          t_2m, td_2m, t_g, clct, clcm, u_10m, v_10m, &
!!          rain_gsp, rain_con, snow_gsp, snow_con,     &
!!          htop_con, hbas_con

  INTEGER, PARAMETER    :: wp = SELECTED_REAL_KIND( 6, 37)  !< single precision
! INTEGER, PARAMETER    :: wp = SELECTED_REAL_KIND(12,307)  !< double precision

  REAL(wp), ALLOCATABLE :: vct_a(:)

  REAL(wp), PARAMETER :: rdv      = 287.04_wp/461.51_wp
  REAL(wp), PARAMETER :: O_m_rdv  = 1._wp - rdv
  REAL(wp), PARAMETER :: rd_o_cpd =  287.04_wp/1004.64_wp
  REAL(wp), PARAMETER :: tmelt    = 273.15
  REAL(wp), PARAMETER :: alvdcp   = 2.5008e6_wp/1004.64_wp
  REAL(wp), PARAMETER :: b3       = tmelt

  REAL(wp), PARAMETER :: b1       = 610.78_wp
  REAL(wp), PARAMETER :: b2w      = 17.269_wp
  REAL(wp), PARAMETER :: b4w      = 35.86_wp

  INTEGER, PARAMETER  :: max_dom = 1
  INTEGER             :: kstart_moist(max_dom) = (/ 0 /)
  REAL(wp),PARAMETER  :: htop_moist_proc = 22500._wp   ! ICON default value
! REAL(wp),PARAMETER  :: htop_moist_proc = 15000._wp

! REAL(wp), ALLOCATABLE:: t (:,:),  &  ! temperature (K)
  REAL(wp), POINTER    :: t (:,:),  &  ! temperature (K)
                          qv(:,:),  &  ! specific humidity (kg/kg)
                          qc(:,:),  &  ! specific cloud liquid water (kg/kg)
                          u (:,:),  &  ! zonal comp. of wind velocity      (m/s)
                          v (:,:),  &  ! meridional comp. of wind velocity (m/s)
                          p (:,:),  &  ! pressure on half levels (Pa)
                          hhl(:,:)     ! height of half levels (m)
!                         pf(:,:)      ! pressure on full levels

! REAL(wp), ALLOCATABLE:: t_2m (:),  &  ! temperature at 2m above the ground (K)
  REAL(wp), POINTER    :: t_2m (:),  &  ! temperature at 2m above the ground (K)
                          td_2m(:),  &  ! dewpoint temperature at 2m above ground (K)
                          t_g  (:),  &  ! surface temperature (weighted from t_s and t_snow) (K)
                          clct (:),  &  ! total cloud cover
                          clcm (:),  &  ! medium cloud cover
                          u_10m(:),  &  ! zonal wind component at 10m above the ground (m/s)
                          v_10m(:),  &  ! meridional wind component at 10m above the ground (m/s)
                          rain_gsp (:), & ! grid scale rain accum. until ivv
                          rain_con (:), & ! convective rain accum. until ivv
                          snow_gsp (:), & ! grid scale snow accum. until ivv
                          snow_con (:)    ! convective snow accum. until ivv

! REAL(wp), ALLOCATABLE:: hbas_con(:),  & ! height of base of main convective cloud
  REAL(wp), POINTER    :: hbas_con(:),  & ! height of base of main convective cloud
                          htop_con(:)     ! height of top of main convective cloud

CONTAINS


  SUBROUTINE vct_a_from_hhl( ie, ke1, hhl, iverb)
    INTEGER, INTENT(IN) :: ie, ke1, iverb
    REAL(wp), INTENT(IN) :: hhl(ie,ke1)
    REAL(wp)             :: h

    INTEGER :: k, ierr, ii
    INTEGER :: imin(1)

    imin(:) = MINLOC( ABS( hhl(:,ke1)) )
    ii = imin(1)
    IF ( iverb > 9) PRINT *, 'vct_a_from_hhl: Min. ABS(HHL) = ', ABS( hhl(ii,ke1)), &
                             ' at index i = ', ii
    IF ( ABS( hhl(ii,ke1)) > 0.05_wp ) THEN
      PRINT *, 'Warning! Minimum of vct_a > 5 cm'
    END IF

    IF ( ALLOCATED( vct_a) ) THEN
      k = SIZE ( vct_a)
      IF ( k /= ke1) THEN
        PRINT *, 'Error. vct_a already allocated with size ',k,', but expected size ', ke1
        STOP 'Error. Wrong dimension of vct_a.'
      END IF
    ELSE
        ALLOCATE( vct_a(ke1), STAT=ierr)
        IF( ierr /= 0) THEN
          PRINT *, 'Error ', ierr,' allocating vct_a, dim.', ke1,' in vct_a_from_hhl'
          STOP 'Error allocating vct_a'
        END IF
    END IF
    vct_a(:) = hhl(ii,:) ! - hhl(ii,ke1)

    WRITE(*,'(/A)') 'HHL above sea level:'
    WRITE(*,'(/a5,2A10)') 'k', 'HHL', 'HH'
    DO k = 1, ke1-1
      WRITE(*,'(I5,2F10.2)') k, vct_a(k), 0.5_wp*( vct_a(k)+vct_a(k+1) )
    END DO
    WRITE(*,'(I5,2F10.2)') k, vct_a(ke1)

    kstart_moist(1) = 1
    DO k = 1, ke1-1
!     discard top levels without data for HHL
      IF ( ABS( vct_a(k)) < 0.1_wp) CYCLE

      h = 0.5_wp*(vct_a(k)+vct_a(k+1))
!     condition h>0 if the top levels of HHL were not read.
      IF( h > 0._wp .AND. h < htop_moist_proc) THEN
        kstart_moist(1) = k
        IF ( htop_moist_proc - h > 1000._wp) THEN
          PRINT *, 'Top level of moist processes to small. h=',h
          PRINT *, 'It should be close to ', htop_moist_proc
          PRINT *, 'Probably, to few level of HHL were read!'
          STOP 'Error defining kstart_moist'
        END IF
        EXIT
      END IF
    END DO

    WRITE(*,'(A,i4)') 'kstart_moist =', kstart_moist(1)

  END SUBROUTINE vct_a_from_hhl

  ELEMENTAL FUNCTION sat_pres_water(temp)
  IMPLICIT NONE

  REAL (wp)              :: sat_pres_water
  REAL (wp), INTENT(IN)  :: temp

  sat_pres_water = b1*EXP( b2w*(temp-b3)/(temp-b4w) )

  END FUNCTION sat_pres_water

END MODULE mo_wwonly
