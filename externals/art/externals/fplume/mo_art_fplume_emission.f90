!
! mo_art_fplume_emission
! This module starts the 1D- volcanic plume model FPlume and prepares the transport of volcanic
! ash.
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

MODULE mo_art_fplume_emission
  USE mo_kind,                          ONLY: wp
  USE mo_art_fplume,                    ONLY: art_fplume
  USE mtime,                            ONLY: datetime,timedelta
  USE mo_art_external_types,            ONLY: t_art_volc_fplume
  USE mo_art_pntSrc_types,              ONLY: t_art_all_pntSrc
  USE mo_art_modes_linked_list,         ONLY: p_mode_state,t_mode
  USE mo_impl_constants,                ONLY: SUCCESS, MAX_CHAR_LENGTH 
  USE mo_exception,                     ONLY: message, message_text, finish
  USE mo_model_domain,                  ONLY: p_patch
  USE mo_parallel_config,               ONLY: nproma
  USE mo_art_fplume_types,              ONLY: t_fplume_phases
  USE mo_art_data,                      ONLY: t_art_data
  USE mo_art_emiss_types,               ONLY: t_art_aero_emiss
  USE mo_art_config,                    ONLY: art_config
  USE mo_art_emiss_types,               ONLY: t_art_emiss2tracer
  USE mo_math_constants,                ONLY: pi,euler
  USE mo_var_metadata_types,            ONLY: t_var_metadata_dynamic
  USE mo_var_list,                      ONLY: get_tracer_info_dyn_by_idx

  PRIVATE
  PUBLIC :: art_fplume_emission

  CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_fplume_emission(p_art_data,volcs,z_mc,z_ifc,rho,pres,temp,u,v,iqv,jg,profile_nz,   &
           &                   current_date,volc_fplume,dtime,dz,cell_area,tracer,emiss_rateM) 

  IMPLICIT NONE

  INTEGER, INTENT(IN)                :: &
    &  volcs,                           & !< volcano index on domain
    &  profile_nz,jg,                   & !< number of ICON levels(nlev), p_patch index
    &  iqv                                !< tracer index specific humidity
  TYPE(t_art_data),TARGET,INTENT(INOUT) :: &
    &  p_art_data                           !< ART data container
  TYPE(t_art_volc_fplume), INTENT(INOUT) :: &
    &  volc_fplume                          !< Container for ash transport properties
  REAL(wp), INTENT(INOUT)                :: &
    &  tracer(:,:,:,:)                      !< Mass mixing ratio of selected tracers
  REAL(wp), INTENT(IN)                   :: &
    &  z_mc(:,:),                       & !< Geometric height (full levels for FPlume)
    &  z_ifc(:,:),                      & !< Geometric height in half levels for emission
    &  rho(:,:),                        & !< Density of air
    &  pres(:,:),                       & !< Air pressure
    &  temp(:,:),                       & !< Air temperature
    &  u(:,:),                          & !< Zonal wind
    &  v(:,:),                          & !< Meridional wind
    &  dtime,                             & !< Model time step (advection)
    &  dz(:,:),                         & !< Layer height
    &  cell_area(:)                       !< Cell area
  TYPE(datetime), POINTER, INTENT(IN):: &
    &  current_date                         !< mtime object containing current date (ICON)
  REAL(wp), INTENT(INOUT)                 :: &
    &  emiss_rateM(:,:,:)                   !< (mass) emission rate (output) [jc,jk,mode]
     ! This probably calls for a general overhaul

!Local variables
  INTEGER,PARAMETER     :: plume_ns = 300             ! number of plume sources (plume+umbrella)
  INTEGER                    ::  &
    &  idz,idz_top,                & !< loop index for comparing FPlume and ICON vertical levels
    &  loc_rad1,loc_rad2,          & !< location of FPlume levels around ICON level                
    &  iSO2,                       & !< SO2 index
    &  jb,jc,jk,                   &
    &  idz_all(plume_ns),          & !< neccessary to construct uppermost emission strenght
    &  fplume_top_index,           & !< neccessary to construct uppermost emission strenght 
    &  loc_rad2_all(plume_ns),     & !< neccessary to construct uppermost emission strenght
    &  NBL_index                     !< Index of Neutral Buoyancy Level
  REAL(wp)                           :: &
    &  plume_MER_icon,                    & !< total Mass Flow Rate 
    &  plume_H_icon,                      & !< plume height in m
    &  plume_z(plume_ns),                 & !< z-axis from FPlume (elevation in m above terrain)
    &  plume_radius(plume_ns),            & !< plume radius
    &  vent_height,                       & !< vent height
    &  plume_z0(plume_ns),                & !< FPlume z-axis (a.s.l)
    &  emis_profile(profile_nz),          & !< emission profile
    &  emis_height_fct(profile_nz),       & !< emission factor for ash and SO2
    &  max_emis,                          & !< maximum value of emissions (profile)
    &  emis_water(profile_nz),            & !< emission profile of water
    &  emis_height_fct_water(profile_nz), & !< emission factor for water
    &  max_emis_water,                    & !< maximum value of water emission
    &  Q_trans,                           & !< amount of ash for transport
    &  rad_rel,                           & !< relative number for interpolation
    &  emiss_fct,                         & !< emission factor
    &  fine_ash_fraction,                 & !< phase dependent fine ash fraction
    &  integral_emission,                 & !< Integral of emission profile (z* from 0 to 1)
    &  integral_emission_water,           & !< Integral of water emission (z* from 0 to 1)
    &  MER_H2O(plume_ns),                 & !< MER of water (emitted+entrainment)
    &  MER_SO2,                           &
    &  total_emiss_rate

  LOGICAL                    :: fplume_on           !< was FPlume active?

  CHARACTER(LEN=MAX_CHAR_LENGTH) :: thisroutine =  "mo_art_fplume_emission:art_fplume_emission" 

  emis_profile(:) = 0.0_wp
  emis_water(:)   = 0.0_wp
  idz_all(:)      = 1000
  loc_rad2_all(:) = 1000

  ! Abbreviations
  jc           = p_art_data%fplume_init_all%p(volcs)%tri_iidx_loc
  jb           = p_art_data%fplume_init_all%p(volcs)%tri_iblk_loc
    
  CALL art_fplume(p_art_data%fplume_init_all%p(volcs),volc_fplume,z_mc(jc,:),z_ifc(jc,:),    &
             &    rho(jc,:),pres(jc,:),temp(jc,:),u(jc,:),v(jc,:),tracer(jc,:,:,:),iqv,      &
             &    jg,jb,jc,profile_nz,current_date,plume_MER_icon,plume_radius,plume_H_icon, &
             &    MER_H2O,MER_SO2,plume_z,vent_height,fplume_on,fine_ash_fraction)
  volc_fplume%MER_transport = 0.0_wp

  IF (fplume_on .AND. art_config(jg)%iart_fplume>0) THEN       ! iart_fplume=-1, no transport

    !------ Calculate Gouhier (2019) fine ash fraction -------
    IF (fine_ash_fraction>0.0_wp .AND. fine_ash_fraction<1.0_wp) THEN
      Q_trans = plume_MER_icon*fine_ash_fraction
    ELSE
      SELECT CASE(int(fine_ash_fraction))
        CASE(0) ! Gouhier general
          Q_trans = (plume_MER_icon/30.22_wp/((plume_H_icon/1000.0_wp)**2.25_wp))            &
                &   **(1.0_wp/0.51_wp)
        CASE(-1) ! low SiO2, closed conduit
          Q_trans = (plume_MER_icon/25.95_wp/((plume_H_icon/1000.0_wp)**1.95_wp))            &
                &   **(1.0_wp/0.72_wp)
        CASE(-2) ! low SiO2, open conduit
          Q_trans = (plume_MER_icon/25.95_wp/((plume_H_icon/1000.0_wp)**1.40_wp))            &
                &   **(1.0_wp/0.72_wp)
        CASE(-3) ! average between case -1 and -2
          Q_trans = (plume_MER_icon/25.95_wp/((plume_H_icon/1000.0_wp)**1.95_wp))            &
                &   **(1.0_wp/0.72_wp) &
                &   + (plume_MER_icon/25.95_wp/((plume_H_icon/1000.0_wp)**1.4_wp))           &
                &   **(1.0_wp/0.72_wp)
          Q_trans = Q_trans / 2.0_wp
        CASE(-4) ! high SiO2, closed conduit
          Q_trans = (plume_MER_icon/25.95_wp/((plume_H_icon/1000.0_wp)**1.95_wp))            &
                &   **(1.0_wp/0.62_wp)
        CASE(-5) ! high SiO2, open conduit
          Q_trans = (plume_MER_icon/25.95_wp/((plume_H_icon/1000.0_wp)**1.40_wp))             &
                &   **(1.0_wp/0.62_wp)
        CASE DEFAULT
          CALL finish(thisroutine,'Invalid fine_ash_fraction')
        END SELECT
    ENDIF

    volc_fplume%MER_transport = Q_trans

    DO idz=1,plume_ns
      plume_z0(idz) = plume_z(idz)+vent_height
    ENDDO

    !------ FPlume vertical axis to ICON vertical axis for radius and MER of water------
    DO idz=1,profile_nz
      IF (plume_z0(1) >= z_mc(jc,idz) .AND. plume_z0(plume_ns) <= z_mc(jc,idz)) THEN
        idz_all(idz)=idz
        loc_rad1 = MINLOC(plume_z0,DIM=1, MASK=plume_z0<=z_mc(jc,idz))
        loc_rad2 = MAXLOC(plume_z0,DIM=1, MASK=plume_z0>=z_mc(jc,idz))
        loc_rad2_all(idz) = loc_rad2
        rad_rel  = (z_mc(jc,idz) - plume_z0(loc_rad1)) &
              &  / (plume_z0(loc_rad2) - z_mc(jc,idz))
!        emis_profile = (plume_radius(loc_rad1) + plume_radius(loc_rad2) * rad_rel) /2
        emis_water(idz) = (MER_H2O(loc_rad1) + MER_H2O(loc_rad2) * rad_rel)/2.0_wp
      ENDIF
    ENDDO

    ! Suzuki Distribution for ash and SO2 (see Marti et al (2017))
    DO jk=1,profile_nz
      emis_profile(jk) = plume_MER_icon*((1.0_wp-(z_ifc(jc,jk)-vent_height)/plume_H_icon)           &
                    &  * EXP(4.0_wp*((z_mc(jc,jk)-vent_height)/plume_H_icon-1)))**5
      IF (emis_profile(jk)<0.0_wp) emis_profile(jk)=0.0_wp
    ENDDO        

    max_emis       = MAXVAL(emis_profile)
    NBL_index      = MAXLOC(emis_profile,DIM=1)

    !----- construct umbrella region when resolution is lower -----
    idz_top=MINVAL(idz_all)
    fplume_top_index=MINVAL(loc_rad2_all)

    IF (z_mc(jc,idz_top-1)-200.0_wp <= plume_z0(fplume_top_index)) THEN ! in uppermost part of grid box
      emis_profile(idz_top-1) = emis_profile(idz_top)/2.0_wp
    ELSE IF (z_mc(jc,idz_top-1)-200.0_wp >= plume_z0(fplume_top_index)      &
            & .AND. z_mc(jc,idz_top)+200.0_wp <= plume_z0(fplume_top_index)) THEN
      emis_profile(idz_top-1) = emis_profile(idz_top)/4.0_wp
    ENDIF

    DO jk=1, profile_nz
      IF (jk<NBL_index .AND. emis_water(jk)/=0.0_wp) THEN
      ! reduce emission of water above NBL --> Linear decrease
        emis_water(jk) = -(z_mc(jc,jk)-plume_H_icon)*emis_water(NBL_index) &
                  &      /(plume_H_icon-z_mc(jc,jk))
      ENDIF
      IF (z_mc(jc,jk)<=7000.0_wp) emis_water(jk)=0.0_wp  ! only emission above 7km
      emis_profile(jk) = emis_profile(jk)/max_emis     ! normalized emission profile
    ENDDO      

    integral_emission       = 0.0_wp
    integral_emission_water = 0.0_wp
    max_emis_water = MAXVAL(emis_water)
    DO jk=1,profile_nz
    ! calc integral and emission factor 
      emis_height_fct(jk)       = emis_profile(jk) *(dz(jc,jk)/plume_H_icon)
      emis_height_fct_water(jk) = emis_water(jk) *(dz(jc,jk)/plume_H_icon)
      integral_emission         = integral_emission + emis_height_fct(jk)
      integral_emission_water   = integral_emission_water + emis_height_fct_water(jk)
    ENDDO

    emiss_fct = Q_trans * 1.0e9_wp / cell_area(jc)        ! Q_trans in kg/s-->amug/s

    DO jk = 1,profile_nz
      ! --- finalizing and distributing emission rates for ash --- !
      total_emiss_rate = emis_height_fct(jk) / integral_emission * emiss_fct &
        &              / dz(jc,jk)
      emiss_rateM(jc,jk,1)   = total_emiss_rate * p_art_data%fplume_init_all%p(volcs)%cfactor_asha
      emiss_rateM(jc,jk,2)   = total_emiss_rate * p_art_data%fplume_init_all%p(volcs)%cfactor_ashb 
      emiss_rateM(jc,jk,3)   = total_emiss_rate * p_art_data%fplume_init_all%p(volcs)%cfactor_ashc 
    
      ! --- adding emission of SO2 to corresponding tracer --- !
      IF (MER_SO2 > 0.0_wp) THEN
        iSO2   = p_art_data%chem%indices%iTRSO2
        tracer(jc,jk,jb,iSO2)   = tracer(jc,jk,jb,iSO2)                                &
          &                     + emis_height_fct(jk)/integral_emission &
          &                     * MER_SO2 * dtime / cell_area(jc)                   &
          &                     / ( rho(jc,jk) * dz(jc,jk) )
      ENDIF

      ! emission of volcanic H2O
      IF (art_config(jg)%iart_fplume==3) THEN
        tracer(jc,jk,jb,iqv) = tracer(jc,jk,jb,iqv)                                      &
                             + (emis_height_fct_water(jk)/integral_emission_water   & 
                             * dtime / cell_area(jc)                                  &
                             / ( rho(jc,jk) * dz(jc,jk) ))
      ENDIF

    ENDDO
  ENDIF ! if fplume_on

  RETURN
END SUBROUTINE art_fplume_emission
!!
!!-------------------------------------------------------------------------
!! 
END MODULE mo_art_fplume_emission

