!
! mo_art_fplume
! This module constucts the atmospheric profiles, searches the actual
! volcanic phase for FPlume, and starts the plume calculations
! (based on FPLUME-1.1 by A.Folch, G.Macedonio, A.Costa, February 2016)
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

MODULE mo_art_fplume

  USE mo_art_fplume_write,              ONLY: art_write_fplume_values
  USE mo_art_fplume_bpt,                ONLY: initialize_plume_BPT,initialize_plume_wind,    &
                                          &   solve_plume_BPT
  USE mo_art_external_types,            ONLY: t_art_volc_fplume
  USE mo_exception,                     ONLY: message, message_text, finish
  USE mo_impl_constants,                ONLY: MAX_CHAR_LENGTH
  USE mo_kind,                          ONLY: wp
  USE mo_math_constants,                ONLY: pi
  USE mtime,                            ONLY: datetime,newDatetime,datetimeToString,         &
                                          &   max_datetime_str_len,newEvent,timedelta,       &
                                          &   event, isCurrentEventActive
  USE mo_art_fplume_types,              ONLY: t_fplume_phases
  USE mo_art_config,                    ONLY: art_config

  PRIVATE
 
  PUBLIC:: art_fplume

  CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_fplume(fplume_init,volc_fplume,z,z_ifc,rho,pres,temp,u,v,sh,iqv,jg,jb,jc,          &
               &      profile_nz,current_date,plume_MER_icon,plume_R_ICON,plume_H_icon,           &
               &      plume_MER_H2O,plume_MER_SO2,plume_z_icon,plume_zv_icon,fplume_on,           &
               &      plume_fine_ash_fraction)
 
  IMPLICIT NONE

  TYPE(t_fplume_phases),INTENT(INOUT) ::  &
    &  fplume_init                          !< FPlume data container for one volcano
  TYPE(t_art_volc_fplume),INTENT(INOUT):: &
    &  volc_fplume                          !< Container for ash transport properties
  TYPE(datetime), POINTER, INTENT(IN)  :: &
    &  current_date                           !< mtime object containing current date (ICON)
  REAL(wp), INTENT(IN)                 :: &
    &  z(:),                            & !< Geometric height
    &  z_ifc(:),                        & !< needed to get z_ifc at vent location
    &  rho(:),                          & !< Density of air
    &  pres(:),                         & !< Air pressure
    &  temp(:),                         & !< Air temperature
    &  u(:),                            & !< Zonal wind
    &  v(:),                            & !< Meridional wind
    &  sh(:,:,:)                          !< specific humidity
  INTEGER, INTENT(IN)                  :: profile_nz,jg,jb,jc,iqv
  LOGICAL, INTENT(INOUT)               :: &
    &  fplume_on                              !< was FPlume active?
  REAL(wp), INTENT(INOUT)              :: plume_MER_H2O(300)

  INTEGER               :: plume_status               ! status code
  INTEGER,PARAMETER     :: plume_ns = 300             ! number of plume sources (plume+umbrella)
  INTEGER,PARAMETER     :: plume_np = 200             ! number of plume sources (up to NBL)
  REAL(wp)    :: plume_MER
  REAL(wp)    :: plume_R0
  REAL(wp)    :: fplume_min_height          ! use FPlume above (below Mastin)
  !
  REAL(wp), ALLOCATABLE :: plume_x   (:)    ! lon(ns)
  REAL(wp), ALLOCATABLE :: plume_y   (:)    ! lat(ns)
  REAL(wp), ALLOCATABLE :: plume_z   (:)    ! z(ns) (elevation in m above terrain)
  REAL(wp), ALLOCATABLE :: plume_Q   (:)    ! Bulk mass flow rate (ns)
  REAL(wp), ALLOCATABLE :: plume_En  (:)    ! Total energy flow rate (ns)
  REAL(wp), ALLOCATABLE :: plume_M   (:,:)  ! Particle mass flow rate (nc,ns)
  REAL(wp), ALLOCATABLE :: plume_Mf  (:,:)  ! Mass flow rate of particles that fall (nc,ns)
  REAL(wp), ALLOCATABLE :: plume_Magr(:,:)  ! Mass flow rate of particles that aggregate (nc,ns)
  REAL(wp), ALLOCATABLE :: plume_Mair(:)    ! Mass flow rate of air (ns)
  REAL(wp), ALLOCATABLE :: plume_Mw  (:)    ! Mass flow rate of volatiles  (ns)
  REAL(wp), ALLOCATABLE :: plume_L   (:)    ! Coordinate s (along the plume centerline)
  REAL(wp), ALLOCATABLE :: plume_H   (:)    ! Theta angle
  REAL(wp), ALLOCATABLE :: plume_U   (:)    ! Bulk velocity
  REAL(wp), ALLOCATABLE :: plume_T   (:)    ! Bulk temperature
  REAL(wp), ALLOCATABLE :: plume_D   (:)    ! Bulk density
  REAL(wp), ALLOCATABLE :: plume_D_rhoa(:)  ! Bulk density/air density
  REAL(wp), ALLOCATABLE :: plume_R   (:)    ! Plume radius
  REAL(wp), ALLOCATABLE :: plume_xv  (:)    ! water vapor  mass fraction
  REAL(wp), ALLOCATABLE :: plume_xl  (:)    ! liquid water mass fraction
  REAL(wp), ALLOCATABLE :: plume_xs  (:)    ! ice (solid)  mass fraction
  REAL(wp), ALLOCATABLE :: plume_as  (:)    ! a_shear
  REAL(wp), ALLOCATABLE :: plume_av  (:)    ! a_vortex
  REAL(wp), ALLOCATABLE :: plume_Mpart(:)   ! total mass of particles
  !
  REAL(wp), ALLOCATABLE :: profile_z(:)       ! Height (a.s.l.)
  REAL(wp), ALLOCATABLE :: profile_rho(:)     ! Air density
  REAL(wp), ALLOCATABLE :: profile_p(:)       ! Pressure
  REAL(wp), ALLOCATABLE :: profile_T(:)       ! Temperature
  REAL(wp), ALLOCATABLE :: profile_sh(:)      ! Specific humidity
  REAL(wp), ALLOCATABLE :: profile_u(:)       ! Wind speed
  REAL(wp), ALLOCATABLE :: profile_ux(:)      ! x-component of wind vector
  REAL(wp), ALLOCATABLE :: profile_uy(:)      ! y-component of wind vector
  REAL(wp), ALLOCATABLE :: profile_psia(:)    ! Wind direction (origin at E,anticlockwise, radians)
  !

! Variables necessary for transport in ICON
  REAL(wp), INTENT(INOUT)   ::            &
    &  plume_zv_icon,                       & !< vent altitude (m) 
    &  plume_MER_icon,                      & !< MER (kg/s) of each eruption phase
    &  plume_H_icon,                        & !< Height (m agl)
    &  plume_z_icon(plume_ns),              & !< z(ns) (elevation in m above terrain)
    &  plume_R_icon(plume_ns),              & !< Plume radius
    &  plume_fine_ash_fraction,             & !< phase dependent fine ash fraction
    &  plume_MER_SO2                          !< MER of SO2 (read from .inp file)
  REAL(wp)                          :: &
    &  prof_z(profile_nz),                &
    &  prof_rho(profile_nz),              &
    &  prof_p(profile_nz),                &
    &  prof_T(profile_nz),                &
    &  prof_sh(profile_nz),               &
    &  prof_ux(profile_nz),               &
    &  prof_uy(profile_nz),               &
    &  angle                             
  INTEGER                   :: &
    &  plume_phase,                       & !< current eruption phase number as defined in .ipn file
    &  i, idt2 
  CHARACTER(len=MAX_CHAR_LENGTH),PARAMETER  :: thisroutine='mo_art_fplume:art_fplume'

  plume_status = 0 ! initialize plume status to avoid error when plume not active at beginning
  fplume_min_height=fplume_init%fplume_min_height
  
!----------- PHASES --------------------------------------------------------
  DO idt2=1,fplume_init%nphases
    IF (fplume_init%phase(idt2)%isActive(current_date) .AND. &
    &   fplume_init%phase(idt2)%plume_udt /=0) THEN
      fplume_on = .TRUE.
      plume_phase = idt2
      EXIT
    ELSE
      fplume_on = .FALSE.
      plume_phase = idt2
    ENDIF
  ENDDO
  

!----------- Vertical atmospheric profiles at eruption location --------------
  ALLOCATE(profile_z   (profile_nz))
  ALLOCATE(profile_rho (profile_nz))
  ALLOCATE(profile_p   (profile_nz))
  ALLOCATE(profile_T   (profile_nz))
  ALLOCATE(profile_sh  (profile_nz))
  ALLOCATE(profile_ux  (profile_nz))
  ALLOCATE(profile_uy  (profile_nz))
  ALLOCATE(profile_u   (profile_nz))
  ALLOCATE(profile_PSIa(profile_nz))
  prof_z(:)=z(:)
  prof_rho(:)=rho(:)
  prof_p(:)=pres(:)
  prof_T(:)=temp(:)
  prof_ux(:)=u(:)
  prof_uy(:)=v(:)
  prof_sh(:)=sh(:,jb,iqv)

!----------- Reverse atmospheric profile -----------------------------------
  DO i = 1,profile_nz
    profile_z   (profile_nz-i+1) = prof_z(i)    ! im m
    profile_rho (profile_nz-i+1) = prof_rho(i)                                                    
    profile_p   (profile_nz-i+1) = prof_p(i)    ! in Pa                                            
    profile_T   (profile_nz-i+1) = prof_T(i)    
    profile_sh  (profile_nz-i+1) = prof_sh(i)        ! kg/kg
    profile_ux  (profile_nz-i+1) = prof_ux(i)                                                      
    profile_uy  (profile_nz-i+1) = prof_uy(i)                                                      
  ENDDO

  DO i = 1,profile_nz
     profile_u(i)  = SQRT(profile_ux(i)*profile_ux(i)+profile_uy(i)*profile_uy(i))
       IF (ABS(profile_ux(i))>1.0e-8_wp) THEN
          angle = ATAN2(profile_uy(i),profile_ux(i))*180.0_wp/pi
       ELSE
          IF (profile_uy(i)>1.0e-8_wp) THEN
             angle = 90.0_wp
          ELSE IF (profile_uy(i)< -1.0e-8_wp) THEN
             angle = 270.0_wp
          ELSE
             angle = 0.0_wp
          ENDIF
       ENDIF
       profile_PSIa(i) = angle*pi/180.0_wp     ! angle in Rad
  ENDDO
  !
  IF (fplume_init%plume_zv<profile_z(1)) &
     & CALL finish(thisroutine,'Vent altitude below surface')
  !
!-------------- Calculate plume properties ----------------------------------------------------------
  IF (fplume_on) THEN
    IF (fplume_init%phase(plume_phase)%plume_Hdt>=fplume_min_height .OR.     & 
      & fplume_init%plume_solve_for =='HEIGHT') THEN 

      CALL initialize_plume_bpt(plume_ns, plume_np, fplume_init)
      !
      CALL initialize_plume_wind(profile_nz,profile_z,profile_rho,profile_p,profile_T,  &
                            &  profile_sh,profile_u,profile_PSIa) 
      !
      ALLOCATE(plume_x     (plume_ns))
      ALLOCATE(plume_y     (plume_ns))
      ALLOCATE(plume_z     (plume_ns))
      ALLOCATE(plume_L     (plume_ns))
      ALLOCATE(plume_H     (plume_ns))
      ALLOCATE(plume_U     (plume_ns))
      ALLOCATE(plume_T     (plume_ns))
      ALLOCATE(plume_D     (plume_ns))
      ALLOCATE(plume_D_rhoa(plume_ns))
      ALLOCATE(plume_R     (plume_ns))
      ALLOCATE(plume_Q     (plume_ns))
      ALLOCATE(plume_En    (plume_ns))
      ALLOCATE(plume_Mw    (plume_ns))
      ALLOCATE(plume_Mair  (plume_ns))
      ALLOCATE(plume_xv    (plume_ns))
      ALLOCATE(plume_xl    (plume_ns))
      ALLOCATE(plume_xs    (plume_ns))
      ALLOCATE(plume_as    (plume_ns))
      ALLOCATE(plume_av    (plume_ns))
      ALLOCATE(plume_M     (fplume_init%nc,plume_ns))
      ALLOCATE(plume_Mf    (fplume_init%nc,plume_ns))
      ALLOCATE(plume_Magr  (fplume_init%nc,plume_ns))
  
      plume_x      = 0.0_wp
      plume_y      = 0.0_wp
      plume_z      = 0.0_wp
      plume_L      = 0.0_wp
      plume_H      = 0.0_wp
      plume_U      = 0.0_wp
      plume_T      = 0.0_wp
      plume_D      = 0.0_wp
      plume_D_rhoa = 0.0_wp
      plume_R      = 0.0_wp
      plume_Q      = 0.0_wp
      plume_En     = 0.0_wp
      plume_Mw     = 0.0_wp
      plume_Mair   = 0.0_wp
      plume_M      = 0.0_wp
      plume_Mf     = 0.0_wp
      plume_Magr   = 0.0_wp
      plume_xv     = 0.0_wp
      plume_xl     = 0.0_wp
      plume_xs     = 0.0_wp
      plume_as     = 0.0_wp
      plume_av     = 0.0_wp

      CALL solve_plume_bpt(fplume_init%phase(plume_phase),plume_status,plume_R0,   &
                      &  plume_MER,plume_x,plume_y,plume_z,plume_Q,plume_En,plume_M,plume_Mf,     &
                      &  plume_Magr,plume_Mair,plume_Mw,plume_L,plume_H,plume_U,plume_T,          &
                      &  plume_D,plume_D_rhoa,plume_R,plume_xv,plume_xl,plume_xs,plume_as,plume_av)

!----------- Reverse z Profile again----------------------------------------
      DO i = 1, plume_ns
        plume_R_icon(plume_ns+1-i)           = plume_R(i)    
        plume_z_icon(plume_ns+1-i)           = plume_z(i)
        plume_MER_H2O(plume_ns+1-i)          = plume_Mw(i)
      ENDDO

      plume_MER_icon         = plume_MER
      plume_zv_icon          = fplume_init%plume_zv 
      plume_H_icon           = plume_z(plume_ns) - fplume_init%plume_zv
      plume_MER_SO2          = fplume_init%phase(plume_phase)%MER_SO2
      plume_fine_ash_fraction= fplume_init%phase(plume_phase)%plume_fine_ash_fraction
      
      ALLOCATE(plume_Mpart(SIZE(plume_z)))

      plume_Mpart=SUM(plume_M,dim=1)

!----------- Write FPLUME output -------------------------------------------
      CALL art_write_fplume_values(volc_fplume,jg,plume_H_icon,plume_MER_icon,jc,jb,    &
                              &    fplume_on)


!----------- Deallocate variables ------------------------------------------
      DEALLOCATE(plume_Mpart )
      DEALLOCATE(plume_x     )
      DEALLOCATE(plume_y     )
      DEALLOCATE(plume_z     )
      DEALLOCATE(plume_L     )
      DEALLOCATE(plume_H     )
      DEALLOCATE(plume_U     )
      DEALLOCATE(plume_T     )
      DEALLOCATE(plume_D     )
      DEALLOCATE(plume_D_rhoa)
      DEALLOCATE(plume_R     )
      DEALLOCATE(plume_Q     )
      DEALLOCATE(plume_En    )
      DEALLOCATE(plume_Mw    )
      DEALLOCATE(plume_Mair  )
      DEALLOCATE(plume_xv    )
      DEALLOCATE(plume_xl    )
      DEALLOCATE(plume_xs    )
      DEALLOCATE(plume_as    )
      DEALLOCATE(plume_av    )
      DEALLOCATE(plume_M     )
      DEALLOCATE(plume_Mf    )
      DEALLOCATE(plume_Magr  )

    ELSE IF ((fplume_init%phase(plume_phase)%plume_Hdt < fplume_min_height .AND.    &
      & fplume_init%plume_solve_for == 'MFR')) THEN

      ALLOCATE(plume_R     (plume_ns))
      ALLOCATE(plume_Mw    (plume_ns))
      ALLOCATE(plume_z     (plume_ns))

      plume_MER_icon  = (3.295_wp*fplume_init%phase(plume_phase)%plume_Hdt/1.e3_wp)  &
                  &     **(4.15_wp)

      DO i = 1, plume_ns
        plume_z(i)    = fplume_init%phase(plume_phase)%plume_Hdt/plume_ns *i
      ENDDO
      plume_Mw(:)          = fplume_init%phase(plume_phase)%plume_wvdt*plume_MER_icon
      plume_zv_icon        = fplume_init%plume_zv
      plume_H_icon         = fplume_init%phase(plume_phase)%plume_Hdt

      DO i = 1,plume_ns
        IF (i<200) THEN
          plume_R(i)  = 0.22_wp * (plume_z(i)/plume_H_icon)**2                            &
                    & + 0.75_wp * (plume_z(i)/plume_H_icon) + 0.03_wp 
        ELSE 
          plume_R(i)  = 1.30_wp * (plume_z(i)/plume_H_icon)**3                            &
                    & - 2.13_wp * (plume_z(i)/plume_H_icon)**2                            &
                    & - 0.02_wp * (plume_z(i)/plume_H_icon) + 1.0_wp
        ENDIF
      ENDDO

      DO i = 1, plume_ns
        plume_z_icon(plume_ns+1-i)           = plume_z(i)
        plume_R_icon(plume_ns+1-i)           = plume_R(i)
        plume_MER_H2O(plume_ns+1-i)          = plume_Mw(i)
      ENDDO

      plume_MER_SO2           = fplume_init%phase(plume_phase)%MER_SO2
      plume_fine_ash_fraction = fplume_init%phase(plume_phase)%plume_fine_ash_fraction

      CALL art_write_fplume_values(volc_fplume,jg,plume_H_icon,plume_MER_icon,jc,jb,fplume_on)

      DEALLOCATE(plume_R)
      DEALLOCATE(plume_Mw)
      DEALLOCATE(plume_z)
    ENDIF 
  ELSE
    plume_MER_icon            = 0.0_wp
    plume_R_icon(:)           = 0.0_wp
    plume_z_icon(:)           = 0.0_wp
    plume_zv_icon             = 0.0_wp
    plume_H_icon              = 0.0_wp
    plume_MER_SO2             = 0.0_wp
    CALL art_write_fplume_values(volc_fplume,jg,plume_H_icon,plume_MER_icon,jc,jb,fplume_on)
  ENDIF !fplume on?
  !
  DEALLOCATE(profile_z   )
  DEALLOCATE(profile_rho )
  DEALLOCATE(profile_p   )
  DEALLOCATE(profile_T   )
  DEALLOCATE(profile_sh  )
  DEALLOCATE(profile_ux  )
  DEALLOCATE(profile_uy  )
  DEALLOCATE(profile_u   )
  DEALLOCATE(profile_PSIa)

  IF (plume_status == 1) THEN
    CALL finish(thisroutine,'COLLAPSE of volcanic plume')
  ELSE
    WRITE (message_text,*) 'Plume characteristics have been calculated successfully'
    CALL message (thisroutine, message_text)
  ENDIF

  RETURN
END SUBROUTINE art_fplume
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_fplume
