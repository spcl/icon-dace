!
! Read and apply optical properties of aerosol climatology
!        by G. Stenchikov (volcanic stratospheric aerosols)
!        This is an adaption of mo_aero_volc of echam6 to icon
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

MODULE mo_bc_aeropt_stenchikov

  USE mo_kind,                   ONLY: wp, i8
  USE mo_model_domain,           ONLY: t_patch
  USE mo_radiation_general,      ONLY: nbndlw, nbndsw
  USE mo_exception,              ONLY: finish, message_text, warning
  USE mo_read_interface,         ONLY: openInputFile, closeFile, read_1D, &
    &                                  read_1D_extdim_time, &
    &                                  read_1D_extdim_extdim_time
  USE mo_latitude_interpolation, ONLY: latitude_weights_li
  USE mo_physical_constants,     ONLY: rgrav, rd
  USE mo_math_constants,         ONLY: deg2rad, pi_2
  USE mo_aes_rad_config,         ONLY: aes_rad_config
  USE mtime,                     ONLY: datetime 
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights, &
       &                               calculate_time_interpolation_weights

  IMPLICIT NONE

  PRIVATE
  PUBLIC                           :: read_bc_aeropt_stenchikov, &
    &                                 add_bc_aeropt_stenchikov

  INTERFACE reorder_stenchikov
    MODULE PROCEDURE reorder_stenchikov_2d
    MODULE PROCEDURE reorder_stenchikov_3d
  END INTERFACE reorder_stenchikov

  REAL(wp), ALLOCATABLE            :: aod_v_s(:,:,:),   & ! volcanic AOD solar
                                      ssa_v_s(:,:,:,:), & ! volcanic SSA solar
                                      asy_v_s(:,:,:,:), & ! volcanic ASY solar
                                      ext_v_s(:,:,:,:), & ! volcanic EXT solar
                                      aod_v_t(:,:,:),   & ! volcanic AOD therm
                                      ssa_v_t(:,:,:,:), & ! volcanic SSA therm
                                      ext_v_t(:,:,:,:), & ! volcanic EXT therm
                                      p_lim_clim(:),    & ! limit lev press. data
                                      r_lat_clim(:)       ! latitudes
  REAL(wp)                         :: r_lat_shift, r_rdeltalat
  INTEGER, SAVE                    :: inm2_time_interpolation=-999999
  INTEGER, PARAMETER               :: lev_clim=40, lat_clim=90, nmonths=2
  REAL(wp), PARAMETER              :: rdog=rd*rgrav

  TYPE(t_time_interpolation_weights) :: tiw

CONTAINS
  !>
  !! SUBROUTINE su_bc_aeropt_stenchikov -- sets up the memory for fields in which
  !! the aerosol optical properties are stored when needed
SUBROUTINE su_bc_aeropt_stenchikov

! allocate memory for optical properties

  ALLOCATE(aod_v_s(nbndsw,0:lat_clim+1,nmonths))
  ALLOCATE(ext_v_s(nbndsw,lev_clim,0:lat_clim+1,nmonths))
  ALLOCATE(ssa_v_s(nbndsw,lev_clim,0:lat_clim+1,nmonths))
  ALLOCATE(asy_v_s(nbndsw,lev_clim,0:lat_clim+1,nmonths))
  ALLOCATE(aod_v_t(nbndlw,0:lat_clim+1,nmonths))
  ALLOCATE(ext_v_t(nbndlw,lev_clim,0:lat_clim+1,nmonths))
  ALLOCATE(ssa_v_t(nbndlw,lev_clim,0:lat_clim+1,nmonths))
  ALLOCATE(p_lim_clim(lev_clim+1))
  ALLOCATE(r_lat_clim(0:lat_clim+1))
  !$ACC ENTER DATA CREATE(aod_v_s, ext_v_s, ssa_v_s, asy_v_s, aod_v_t, ext_v_t, ssa_v_t) &
  !$ACC   CREATE(p_lim_clim, r_lat_clim)

! initialize with zero
  aod_v_s(:,:,:)   = 0._wp
  ext_v_s(:,:,:,:) = 0._wp
  ssa_v_s(:,:,:,:) = 0._wp
  asy_v_s(:,:,:,:) = 0._wp
  aod_v_t(:,:,:)   = 0._wp
  ext_v_t(:,:,:,:) = 0._wp
  ssa_v_t(:,:,:,:) = 0._wp
  p_lim_clim(lev_clim+1) = 0._wp
  !$ACC UPDATE DEVICE(aod_v_s, ext_v_s, ssa_v_s, asy_v_s, aod_v_t, ext_v_t, ssa_v_t) &
  !$ACC   DEVICE(p_lim_clim) &
  !$ACC   ASYNC(1)
END SUBROUTINE su_bc_aeropt_stenchikov

  !> SUBROUTINE shift_months_bc_aeropt_stenchikov -- shifts months in order to read a new one.

SUBROUTINE shift_months_bc_aeropt_stenchikov

  aod_v_s(:,:,1)   = aod_v_s(:,:,2)
  ext_v_s(:,:,:,1) = ext_v_s(:,:,:,2)
  ssa_v_s(:,:,:,1) = ssa_v_s(:,:,:,2)
  asy_v_s(:,:,:,1) = asy_v_s(:,:,:,2)
  aod_v_t(:,:,1)   = aod_v_t(:,:,2)
  ext_v_t(:,:,:,1) = ext_v_t(:,:,:,2)
  ssa_v_t(:,:,:,1) = ssa_v_t(:,:,:,2)
  
END SUBROUTINE shift_months_bc_aeropt_stenchikov

  !> SUBROUTINE read_bc_aeropt_stenchikov -- read the aerosol optical properties 
  !! of the volcanic (Stratospheric) Stenchikov aerosols

SUBROUTINE read_bc_aeropt_stenchikov(current_date, p_patch)
  TYPE(datetime), POINTER, INTENT(in) :: current_date
  TYPE(t_patch)          , INTENT(in) :: p_patch

  !LOCAL VARIABLES
  INTEGER(i8) :: iyear(2)
  INTEGER :: imonth(2), nmonths, imonths

  tiw = calculate_time_interpolation_weights(current_date)  

  IF (tiw%month2_index == inm2_time_interpolation) RETURN

  IF (ALLOCATED(aod_v_s)) THEN

    ! skip shifting (and reading) of data if we have a turn of the year
    ! but update control index. This assumes that read_bc_aeropt_stenchikov
    ! is called with a sufficiently high period of less than half a month,
    ! half of the time intervall at which data (here monthly) are provided.

    IF ( inm2_time_interpolation == 13 .AND. tiw%month2_index == 1 ) THEN
       inm2_time_interpolation=tiw%month2_index
       RETURN
    ENDIF

    CALL shift_months_bc_aeropt_stenchikov
    imonth(2)=tiw%month2_index
    iyear(2)=current_date%date%year
    IF (imonth(2) == 13 ) THEN
      imonth(2)=1
      iyear(2)=iyear(2)+1
    END IF

    CALL read_months_bc_aeropt_stenchikov (p_patch,                     &
                       'longitude',       'latitude',         'levels', &
                       imonth(2),    iyear(2),    2  )

    inm2_time_interpolation=tiw%month2_index

  ELSE

    CALL su_bc_aeropt_stenchikov
    imonth(1)=tiw%month1_index
    imonth(2)=tiw%month2_index
    iyear(1)=current_date%date%year
    iyear(2)=current_date%date%year
    IF (imonth(1) == 0) THEN
      imonth(1)=12
      iyear(1)=iyear(1)-1
    END IF
    IF (imonth(2) == 13) THEN
      imonth(2)=1
      iyear(2)=iyear(2)+1
    END IF
    nmonths=2
    inm2_time_interpolation=tiw%month2_index

    DO imonths=1,nmonths
      CALL read_months_bc_aeropt_stenchikov (p_patch,                     &
                         'longitude',       'latitude',         'levels', &
                         imonth(imonths),    iyear(imonths),    imonths   )
    END DO

  ENDIF
  !$ACC UPDATE DEVICE(aod_v_s, ext_v_s, ssa_v_s, asy_v_s, aod_v_t, ext_v_t, ssa_v_t) &
  !$ACC   DEVICE(p_lim_clim, r_lat_clim) &
  !$ACC   ASYNC(1)

END SUBROUTINE read_bc_aeropt_stenchikov
!-------------------------------------------------------------------------
!> SUBROUTINE add_bc_aeropt_stenchikov

!! add aerosol optical properties for all wave length bands (solar and IR)
!! in the case of volcanic aerosols of Stenchikov
!! The height profile is taken into account.
SUBROUTINE add_bc_aeropt_stenchikov(current_date,       jg,               &
          & jcs, kproma,            kbdim,              klev,             &
          & krow,                   nb_sw,              nb_lw,            &
          & dz,                     pp_fl,                                &
          & paer_tau_sw_vr,         paer_piz_sw_vr,     paer_cg_sw_vr,    &
          & paer_tau_lw_vr                                                )

  ! !INPUT PARAMETERS
  TYPE(datetime), POINTER, INTENT(in) :: current_date
  INTEGER,INTENT(in)  :: jg,     &! domain index
                         jcs,    &! actual block length (start)
                         kproma, &! actual block length (end)
                         kbdim,  &! maximum block length
                         krow,   &! block index
                         klev,   &! number of vertical levels
                         nb_lw,  &! number of wave length bands (far IR)
                         nb_sw    ! number of wave length bands (solar)
  REAL(wp),INTENT(in) :: dz(kbdim,klev),      & ! geometric height thickness [m]
                         pp_fl(kbdim,klev)      ! pressure at "full levels"
! !OUTPUT PARAMETERS
  REAL(wp),INTENT(inout),DIMENSION(kbdim,klev,nb_lw):: &
   paer_tau_lw_vr      !aerosol optical depth (far IR)
  REAL(wp),INTENT(inout),DIMENSION(kbdim,klev,nb_sw):: &
   paer_tau_sw_vr,   & !aerosol optical depth (solar), sum_i(tau_i)
   paer_piz_sw_vr,   & !weighted sum of single scattering albedos, 
                       !sum_i(tau_i*omega_i)
   paer_cg_sw_vr       !weighted sum of asymmetry factors, 
                       !sum_i(tau_i*omega_i*g_i)

! !LOCAL VARIABLES
  
  INTEGER, PARAMETER                    :: norder=-1 ! latitudes in climatology order from N->S
  INTEGER, PARAMETER                    :: nm1=1, nm2=2 ! there are only two months of a year
                                                        ! stored in the field (saves memory)
  INTEGER                               :: jl,jk,jki,jwl
  INTEGER                               :: idx_lat_1, idx_lat_2, idx_lev
  REAL(wp)                              :: w1_lat, w2_lat
  INTEGER,  DIMENSION(kbdim,klev)       :: kindex ! index field for pressure interpolation
  REAL(wp), DIMENSION(kbdim)            :: wgt1_lat,wgt2_lat
  INTEGER,  DIMENSION(kbdim)            :: inmw1_lat, inmw2_lat 
  REAL(wp), DIMENSION(kbdim,nb_sw)      :: zaod_s, zext_s_int, zfact_s
  REAL(wp), DIMENSION(kbdim,klev,nb_sw) :: zext_s, zomg_s, zasy_s
  REAL(wp), DIMENSION(kbdim,klev,nb_lw) :: zext_t, zomg_t
  REAL(wp), DIMENSION(kbdim,nb_lw)      :: zaod_t, zext_t_int, zfact_t 
  REAL(wp)                              :: p_lat_shift, p_rdeltalat
  INTEGER                               :: jc

  tiw = calculate_time_interpolation_weights(current_date)

! It is assumed that the pressure levels of the climatology do not change with time but
! are unequally spaced. Since the pressure of each icon level may change with time,
! each specific icon level may have its centre in a different level of the climatology at
! each time step.

! 1. calculate for each icon gridbox the index of the data set layer 
!     in which p_mid_icon is located and geometrical height of layers
  CALL pressure_index(jcs, kproma,   kbdim,         klev,              &
                      pp_fl,         lev_clim,      p_lim_clim,        &
                      kindex)
  p_lat_shift=r_lat_shift
  p_rdeltalat=r_rdeltalat
  CALL latitude_weights_li(jg                  ,jcs                             &
                        & ,kproma              ,kbdim            ,krow          &
                        & ,wgt1_lat            ,wgt2_lat         ,inmw1_lat     &
                        & ,inmw2_lat           ,p_lat_shift      ,p_rdeltalat   &
                        & ,r_lat_clim          ,lat_clim         ,norder        )
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
           zext_s(jl,jk,jwl)=tiw%weight1*(w1_lat*ext_v_s(jwl,idx_lev,idx_lat_1,nm1) + &
                                          w2_lat*ext_v_s(jwl,idx_lev,idx_lat_2,nm1))+ &
                             tiw%weight2*(w1_lat*ext_v_s(jwl,idx_lev,idx_lat_1,nm2) + &
                                          w2_lat*ext_v_s(jwl,idx_lev,idx_lat_2,nm2))
           zomg_s(jl,jk,jwl)=tiw%weight1*(w1_lat*ssa_v_s(jwl,idx_lev,idx_lat_1,nm1) + &
                                          w2_lat*ssa_v_s(jwl,idx_lev,idx_lat_2,nm1))+ &
                             tiw%weight2*(w1_lat*ssa_v_s(jwl,idx_lev,idx_lat_1,nm2) + &
                                          w2_lat*ssa_v_s(jwl,idx_lev,idx_lat_2,nm2))
           zasy_s(jl,jk,jwl)=tiw%weight1*(w1_lat*asy_v_s(jwl,idx_lev,idx_lat_1,nm1) + &
                                          w2_lat*asy_v_s(jwl,idx_lev,idx_lat_2,nm1))+ &
                             tiw%weight2*(w1_lat*asy_v_s(jwl,idx_lev,idx_lat_1,nm2)+  &
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
        zaod_s(jl,jwl)=tiw%weight1*(w1_lat*aod_v_s(jwl,idx_lat_1,nm1) + &
                                    w2_lat*aod_v_s(jwl,idx_lat_2,nm1))+ &
                       tiw%weight2*(w1_lat*aod_v_s(jwl,idx_lat_1,nm2) + &
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
! 2.1 interpolate optical properties solar radiation
  DO jwl=1,nb_lw
     DO jk=1,klev
        DO jl=jcs,kproma
           idx_lat_1=inmw1_lat(jl)
           idx_lat_2=inmw2_lat(jl)
           w1_lat=wgt1_lat(jl)
           w2_lat=wgt2_lat(jl)
           idx_lev=kindex(jl,jk)
           zext_t(jl,jk,jwl)=tiw%weight1*(w1_lat*ext_v_t(jwl,idx_lev,idx_lat_1,nm1) + &
                                          w2_lat*ext_v_t(jwl,idx_lev,idx_lat_2,nm1))+ &
                             tiw%weight2*(w1_lat*ext_v_t(jwl,idx_lev,idx_lat_1,nm2) + &
                                          w2_lat*ext_v_t(jwl,idx_lev,idx_lat_2,nm2))
           zomg_t(jl,jk,jwl)=tiw%weight1*(w1_lat*ssa_v_t(jwl,idx_lev,idx_lat_1,nm1) + &
                                          w2_lat*ssa_v_t(jwl,idx_lev,idx_lat_2,nm1))+ &
                             tiw%weight2*(w1_lat*ssa_v_t(jwl,idx_lev,idx_lat_1,nm2) + &
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
        zaod_t(jl,jwl)=tiw%weight1*(w1_lat*aod_v_t(jwl,idx_lat_1,nm1) + &
                                    w2_lat*aod_v_t(jwl,idx_lat_2,nm1))+ &
                       tiw%weight2*(w1_lat*aod_v_t(jwl,idx_lat_1,nm2) + &
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
!!$     paer_tau_lw_vr(jcs:kproma,jk,1:nb_lw)=paer_tau_lw_vr(jcs:kproma,jk,1:nb_lw) + &
!!$          zext_t(jcs:kproma,jki,1:nb_lw)*(1._wp-zomg_t(jcs:kproma,jki,1:nb_lw))
     ! use explicit DO loop to circumvent possible SXf90 compiler bug
     ! this should be equivalent to the formulation above.
     DO jc = jcs,kproma
     paer_tau_lw_vr(jc,jk,1:nb_lw)=paer_tau_lw_vr(jc,jk,1:nb_lw) + &
          zext_t(jc,jki,1:nb_lw)*(1._wp-zomg_t(jc,jki,1:nb_lw))
     ENDDO
  END DO  

END SUBROUTINE add_bc_aeropt_stenchikov

!-------------------------------------------------------------------------
! 
!> SUBROUTINE pressure_index -- ! find index of pressure layer in which mid level 
!    pressures of icon are located
SUBROUTINE pressure_index(jcs, kproma,   kbdim,         klev,              &
                          pp_mid,        klevels,       pp_bound,          &
                          kindex)

! !INPUT PARAMETERS
  INTEGER, INTENT(in)    :: kbdim, jcs, kproma, klev
  INTEGER, INTENT(in)    :: klevels !number of layers for indices are searched
  REAL(wp), INTENT(in)   :: pp_mid(kbdim,klev), & !aes midlevel pressures
                            pp_bound(klevels+1) !pressure at layer 
                                    !bounds of reference pressures
  INTEGER, INTENT(out)   :: kindex(kbdim,klev) !layer indices for aes press.

! !LOCAL VARIABLES

  LOGICAL                :: lp(kproma), lrepeat
  INTEGER                :: jk,il,kidx(kbdim)
  
  kidx(jcs:kproma)=2
  DO jk=1,klev
10   CONTINUE
     lrepeat=.FALSE.
     DO il=jcs,kproma
        lp(il)=pp_mid(il,jk).GT.pp_bound(kidx(il)).AND.kidx(il).LE.klevels
     END DO
     DO il=jcs,kproma
        IF (lp(il)) THEN
           kidx(il)=kidx(il)+1
           lrepeat=.TRUE.
        END IF
     END DO
     IF (lrepeat) THEN
        GOTO 10
     ELSE
        kindex(jcs:kproma,jk)=kidx(jcs:kproma)-1
     END IF
  END DO
END SUBROUTINE pressure_index

!-------------------------------------------------------------------------
! 
!> SUBROUTINE read_months_bc_aeropt_stenchikov -- reads optical aerosol parameters from 
!    original file containing the following quantities (attention: longitude means wavelengths!)
!!   tautl: monthly average of zonal mean total aod for all wavelengths tautl(time, latitude, longitude) 
!!   exts: extinction for each level corresp. to tautl in 1/m: exts(time, levels, latitude, longitude) 
!!   omega: single scattering albedo corresp. to tautl: omega(time, levels, latitude, longitude)
!!   asymm: asymmetry fractor corresp. to tautl: asymm(time, levels, latitude, longitude)
!!   levels: non-equidistant pressure levels corresp. to exts, omega, asymm: levels(levels).
!!     Only the mid-point pressures are given.
  SUBROUTINE read_months_bc_aeropt_stenchikov (p_patch, &
    cwave_dim,        clat_dim,           clev_dim,  &
    kmonth,           kyear,            ktime_step   )
  !
  TYPE(t_patch)   , INTENT(in)   :: p_patch
  CHARACTER(len=*), INTENT(in)   :: cwave_dim,  &! name of wavelength dimension
                                    clat_dim,   &! name of latitude dimension
                                    clev_dim     ! name of level dimension
  INTEGER, INTENT(in)            :: kmonth       ! number of month to be read
  INTEGER(i8), INTENT(in)        :: kyear        ! year of month
  INTEGER, INTENT(in)            :: ktime_step   ! month that has to be set by new data

  CHARACTER(len=256)             :: cfname   ! file name containing variables
  CHARACTER(len=32)              :: ckmonth,ckyear
  CHARACTER(len=256)             :: cdim_names(4)

!!$  CHARACTER(len=*), INTENT(in), OPTIONAL     :: casl ! name of variable containing altitude of layer centres

  INTEGER                        :: jg
  INTEGER                        :: ifile_id
  REAL(wp), POINTER              :: zvar2d(:,:,:), zvar3d(:,:,:,:)
  REAL(wp), POINTER              :: zpmid(:), zlat(:)
!!$  REAL(wp), POINTER              :: zaod(:,:,:,:), zssa(:,:,:,:), zasy(:,:,:,:), zaer_ex(:,:,:,:)
!!$  CHARACTER(LEN=32)              :: cimnthb, cimnthe
!!$  CHARACTER(LEN=512)             :: cfnameyear,cyear

  jg = p_patch%id

  IF (kmonth < 1 .OR. kmonth > 12 ) THEN
    WRITE(ckmonth,*) kmonth
    CALL finish ('read_months_bc_aeropt_stenchikov in mo_bc_aeropt_stenchikov', &
                 'month to be read outside valid range 1<=kmonth<=12, '// &
                 'kmonth='//TRIM(ADJUSTL(ckmonth))) 
  END IF

  IF ( aes_rad_config(jg)%irad_aero == 13 .OR. &
       aes_rad_config(jg)%irad_aero == 14 .OR. &
       aes_rad_config(jg)%irad_aero == 15 .OR. &
       aes_rad_config(jg)%irad_aero == 18 ) THEN

    WRITE(ckyear,*) kyear

    cfname='bc_aeropt_stenchikov_lw_b16_sw_b14_'//TRIM(ADJUSTL(ckyear))//'.nc'

  ELSE

    cfname='bc_aeropt_stenchikov_lw_b16_sw_b14.nc'

  ENDIF

  CALL openInputFile(ifile_id,cfname)
  cdim_names(1)=cwave_dim
  cdim_names(2)=clat_dim
  cdim_names(3)='time'
  CALL read_1D_extdim_time(file_id=ifile_id, variable_name='tauttl', &
         &                 return_pointer=zvar2d, dim_names=cdim_names(1:3), &
         &                 start_timestep=kmonth, end_timestep=kmonth)
  CALL reorder_stenchikov(zvar2d(:,:,1), &
       &                  aod_v_s(:,:,ktime_step), aod_v_t(:,:,ktime_step))
  DEALLOCATE(zvar2d)
  cdim_names(4)=cdim_names(3)
  cdim_names(3)=clev_dim
  CALL read_1D_extdim_extdim_time(file_id=ifile_id, variable_name='exts', &
    &                             return_pointer=zvar3d, dim_names=cdim_names, &
    &                             start_timestep=kmonth, end_timestep=kmonth)
  CALL reorder_stenchikov(zvar3d(:,:,:,1), &
       &                  ext_v_s(:,:,:,ktime_step), ext_v_t(:,:,:,ktime_step))
  CALL read_1D_extdim_extdim_time(file_id=ifile_id, variable_name='omega', &
    &                             fill_array=zvar3d, dim_names=cdim_names, &
    &                             start_timestep=kmonth, end_timestep=kmonth)
  CALL reorder_stenchikov(zvar3d(:,:,:,1), &
       &                  ssa_v_s(:,:,:,ktime_step), ssa_v_t(:,:,:,ktime_step))
  CALL read_1D_extdim_extdim_time(file_id=ifile_id, variable_name='asymm', &
    &                             fill_array=zvar3d, dim_names=cdim_names, &
    &                             start_timestep=kmonth, end_timestep=kmonth)
  CALL reorder_stenchikov(zvar3d(:,:,:,1), asy_v_s(:,:,:,ktime_step))
  DEALLOCATE(zvar3d)
  CALL read_1D(file_id=ifile_id, variable_name=clev_dim, return_pointer=zpmid)
  CALL p_lim_stenchikov(zpmid)
  DEALLOCATE(zpmid)
  CALL read_1D(file_id=ifile_id, variable_name=clat_dim, return_pointer=zlat)
  IF (SIZE(zlat)/=lat_clim) THEN
    WRITE(message_text, '(a,i0,a,i0,a)') 'lat_clim= ', lat_clim, &
         ' expected but found ', SIZE(zlat), ' elements.'
    CALL finish('read_months_bc_aeropt_stenchikov of mo_bc_aeropt_stenchikov', &
         message_text)
  END IF
  r_lat_clim(1:lat_clim)=zlat(lat_clim:1:-1)*deg2rad
  r_lat_clim(0)         = pi_2                   ! N-pole
  r_lat_clim(lat_clim+1)=-pi_2                   ! S-pole
  r_lat_shift=r_lat_clim(1)                      ! latitude next to the N-pole (so +89)
  r_rdeltalat=ABS(1.0_wp/(r_lat_clim(2)-r_lat_clim(1)))
  DEALLOCATE(zlat)
  CALL closeFile(ifile_id)
  END SUBROUTINE read_months_bc_aeropt_stenchikov
!-------------------------------------------------------------------------
! 
!> SUBROUTINE reorder_stenchikov -- reorder dimensions and distribute on
!!   module variables for time step time_step

  SUBROUTINE reorder_stenchikov_2d (pvar,var_solar,var_thermal)

    REAL(wp), INTENT(in)          :: pvar(:,:)

    REAL(wp), INTENT(out)         :: var_solar(:,:), var_thermal(:,:)

    ! attention: the indices of var_solar and var_thermal are (1:nbnd{sw,lw},1:lat_clim+2)
    var_solar(1:nbndsw-1,2:lat_clim+1)=pvar(nbndsw-1:1:-1,lat_clim:1:-1)
    var_solar(nbndsw,2:lat_clim+1)=pvar(nbndsw,lat_clim:1:-1)
    var_solar(1:nbndsw,1)=var_solar(1:nbndsw,1)
    var_solar(1:nbndsw,lat_clim+2)=var_solar(1:nbndsw,lat_clim)
    var_thermal(1:nbndlw-1,2:lat_clim+1)=pvar(nbndsw+nbndlw:nbndsw+2:-1,lat_clim:1:-1)
    var_thermal(nbndlw,2:lat_clim+1)=pvar(nbndsw-1,lat_clim:1:-1)
    var_thermal(1:nbndlw,1)=var_thermal(1:nbndlw,1)
    var_thermal(1:nbndlw,lat_clim+2)=var_thermal(1:nbndlw,lat_clim)
  END SUBROUTINE reorder_stenchikov_2d

!-------------------------------------------------------------------------
! 
!> SUBROUTINE reorder_stenchikov -- reorder dimensions and distribute on
!!   module variables for time step time_step

  SUBROUTINE reorder_stenchikov_3d(pvar, var_solar, var_thermal)

    REAL(wp), INTENT(in)            :: pvar(:,:,:)
    REAL(wp), INTENT(out)           :: var_solar(:,:,:)
    REAL(wp), OPTIONAL, INTENT(out) :: var_thermal(:,:,:)

    INTEGER :: ibnd, ilev, ilat, nbnd
    ! attention: the indices of var_solar and var_thermal are (1:nbnd{sw,lw},1:lev_clim,1:lat_clim+2)
    DO ilat = 2, lat_clim+1
      DO ilev = 1, lev_clim
        DO ibnd = 1, nbndsw-1
          var_solar(ibnd, ilev, ilat) = &
               pvar(nbndsw-ibnd, lat_clim+2-ilat, lev_clim+1-ilev)
        END DO
        var_solar(nbndsw, ilev, ilat) = &
             pvar(nbndsw, lat_clim+2-ilat, lev_clim+1-ilev)
      END DO
    END DO
    DO ilev = 1, lev_clim
      DO ibnd = 1, nbndsw
        ! this assignment probably is in error and should use 3rd index = 2 on the rhs
        ! var_solar(ibnd, ilev,          1) = var_solar(ibnd, ilev,        1)
!        var_solar(ibnd, ilev,          1) = 0.0_wp
        var_solar(ibnd, ilev, lat_clim+2) = var_solar(ibnd, ilev, lat_clim)
      END DO
    END DO

    IF (PRESENT(var_thermal)) THEN
      nbnd = nbndsw+nbndlw
      DO ilat = 2, lat_clim+1
        DO ilev = 1, lev_clim
          DO ibnd = 1, nbndlw-1
            var_thermal(ibnd, ilev, ilat) = &
                 pvar(nbnd+1-ibnd, lat_clim+2-ilat, lev_clim+1-ilev)
          END DO
          ! this should probably be nbndsw+1 as first index of the rhs
          var_thermal(nbndlw, ilev, ilat) = &
               pvar(nbndsw-1, lat_clim+2-ilat, lev_clim+1-ilev)
        END DO
      END DO
      DO ilev = 1, lev_clim
        DO ibnd = 1, nbndlw
          ! this assignment probably is in error and should use 3rd index = 2 on the rhs
          !var_thermal(ibnd, ilev,          1) = var_thermal(ibnd, ilev,        1)
!          var_thermal(ibnd, ilev,          1) = 0.0_wp
          var_thermal(ibnd, ilev, lat_clim+2) = var_thermal(ibnd, ilev, lat_clim)
        END DO
      END DO

    END IF
  END SUBROUTINE reorder_stenchikov_3d

!-------------------------------------------------------------------------
! 
!> SUBROUTINE p_lim_stenchikov -- calculate pressure at layer boundaries

  SUBROUTINE p_lim_stenchikov (p_pmid)

    REAL(wp), INTENT(inout)        :: p_pmid(lev_clim)

    INTEGER                        :: jk
    REAL(wp) :: accum

    accum = 0.0_wp
    DO jk=1,lev_clim
      p_lim_clim(jk) = accum
      ! convert pressure into Pa from hPa
      accum = 2.0_wp*100._wp*p_pmid(lev_clim+1-jk) - accum
    END DO
    p_lim_clim(lev_clim+1) = accum

  END SUBROUTINE p_lim_stenchikov
    
!-------------------------------------------------------------------------
! 

END MODULE mo_bc_aeropt_stenchikov
