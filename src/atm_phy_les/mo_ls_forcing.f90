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

! This module initializes and applies the large-scale forcing for idealized simulations
! This module assumes that the model grid is FLAT
! 2013-JUNE-04: AT THIS STAGE LS FORCING WILL WORK IN RESTART MODE ONLY IF ITS CALLED EVERY
! DYNAMIC TIMESTEP. TO MAKE IT WORK "SMOOTHLY" ADD_VAR HAS TO WORK ON 1D VARS

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_ls_forcing

  USE mo_kind,                ONLY: wp
  USE mo_io_units,            ONLY: find_next_free_unit
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_impl_constants,      ONLY: success, max_char_length
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_model_domain,        ONLY: t_patch
  USE mo_parallel_config,     ONLY: nproma
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_nh_vert_interp_les,  ONLY: vert_intp_linear_1d, global_hor_mean, &
    &                               vertical_derivative
  USE mo_ls_forcing_nml,      ONLY: is_subsidence_moment, is_subsidence_heat, is_ls_forcing, &
    &                               is_advection, is_advection_uv,is_advection_tq,           &
    &                               is_nudging, is_nudging_uv, is_nudging_tq, dt_relax,      &
    &                               is_geowind, is_rad_forcing 
  USE mo_fortran_tools,       ONLY: init
  USE mo_scm_nml,             ONLY: i_scm_netcdf, lscm_ls_forcing_ini
  USE mo_read_interface,      ONLY: nf
  USE mo_netcdf
  USE mo_mpi,                 ONLY: get_my_global_mpi_id

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_ls_forcing, apply_ls_forcing

  !Anurag Dipankar, MPIM (2013-May): variables for large-scale (LS)
  !forcing to be used in periodic domain or other relveant cases.

  !Dependent on z and time only
  REAL(wp), SAVE, ALLOCATABLE ::  &
    w_ls(:,:),            & !subsidence          [m/s]
    u_geo  (:,:),         & !u-geostrophic wind  [m/s]
    v_geo  (:,:),         & !v-geostrophic wind  [m/s]
    ddt_temp_hadv_ls(:,:),& !LS horizontal advective tendency for temp [K/s]
    ddt_temp_rad(:,:),    & !radiative tendency for temp [k/s]
    ddt_qv_hadv_ls(:,:),  & !LS horizontal advective tendency for qv [1/s]
    ddt_u_hadv_ls(:,:),   & !LS horizontal advective tendency for u [m/s^2]
    ddt_v_hadv_ls(:,:),   & !LS horizontal advective tendency for m [m/s^2]
    u_nudg(:,:),          & !LS nudging u    [m/s]
    v_nudg(:,:),          & !LS nudging v    [m/s]
    theta_nudg(:,:),      & !LS nudging temp [K]
    qv_nudg(:,:),         & !LS nudging qv   [kg/kg]
    bnd_sfc_lat_flx(:),   & !surface Latent heat flux[W/m2]
    bnd_sfc_sens_flx(:),  & !surface sensible heat flux[W/m2]
    bnd_ts(:),            & !Surface temperature[K]
    bnd_tg(:),            & !Skin temperature[K]
    bnd_qvs(:),           & !Surface specific vapor humidity[kg/kg]
    bnd_Ch(:),            & !Drag coefficient for heat
    bnd_Cm(:),            & !Drag coefficient for momentum
    bnd_Cq(:),            & !Drag coefficient for moisture
    bnd_ustar(:),         & !Friction velocity[m/s]
    tempf(:,:,:,:),       & 
    tempf_f(:,:,:,:),     & 
    tempf_sf(:,:,:)

  ! To make sure that LS forcing is not applied before it is read in
  REAL(wp), SAVE :: dt_forcing = 0._wp
  REAL(wp), SAVE :: dt_nudging = 0._wp

  INTEGER,  SAVE :: nt      ! number of time steps in init_SCM.nc


  CONTAINS


  !>
  !! init_ls_forcing
  !!------------------------------------------------------------------------
  !! Initialize large-scale forcings from input file
  !!------------------------------------------------------------------------

  SUBROUTINE init_ls_forcing(p_metrics)

    TYPE(t_nh_metrics),INTENT(in), TARGET :: p_metrics

    CHARACTER(len=*), PARAMETER :: routine = 'mo_ls_forcing:init_ls_forcing'

    REAL(wp), ALLOCATABLE, DIMENSION(:) :: zz, zw, zu, zv, z_dt_temp_adv, &
                                         & z_dt_temp_rad, z_dt_qv_adv,z_dt_u_adv,z_dt_v_adv
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zw_nc, zu_nc, zv_nc, z_dt_temp_adv_nc, &
                                         & z_dt_temp_rad_nc, z_dt_qv_adv_nc,ztheta_nc,z_qv_nc, &
                                         & z_dt_u_adv_nc,z_dt_v_adv_nc, zz_nc
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: time_nc


    REAL(wp), ALLOCATABLE, DIMENSION(:) :: ztheta, z_qv
    REAL(wp)  :: end_time
    INTEGER   :: iunit, ist, nk, jk, nlev, i, nskip, n
    INTEGER   :: lat,lon,t0
    INTEGER   :: varid                  !< id of variable in netcdf file
    INTEGER   :: fileid                 !< id number of netcdf file
    INTEGER   :: dimid                  !< id number of dimension
    INTEGER   :: nf_status, nf_status2  !< return status of netcdf function

    CHARACTER(len=max_char_length) :: time_unit, time_unit_short
    REAL(wp)                       :: time_unit_in_sec
    REAL(wp)                       :: temp_nf(1)

    nlev   = SIZE(p_metrics%z_mc,2)

!------------------------------------------------------------------------------
! READ LS FORCING FILE
!------------------------------------------------------------------------------

    IF (is_ls_forcing) THEN

      IF (i_scm_netcdf > 0) THEN

!------------------------------------------------------------------------------
! Forcing with unified SCM NETCDF file

        IF (i_scm_netcdf == 2) THEN

          CALL nf (nf90_open('init_SCM.nc', NF90_NOWRITE, fileid), &
            & TRIM(routine)//'   File init_SCM.nc cannot be opened')

          CALL nf (nf90_inq_dimid(fileid, 'lev', dimid), TRIM(routine)//' lev')
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = nk), routine)

          CALL nf (nf90_inq_dimid(fileid, 'lat', dimid), routine)
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = lat), routine)

          CALL nf (nf90_inq_dimid(fileid, 'lon', dimid), routine)
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = lon), routine)

          CALL nf (nf90_inq_dimid(fileid, 't0', dimid), routine)
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = t0), routine)

          CALL nf (nf90_inq_dimid(fileid, 'time', dimid), TRIM(routine)//' nt')
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = nt), routine)
          
          WRITE(message_text,'(a,i6,a,i6)') 'SCM LS forcing input file: nk, ', nk, ', nt ', nt
          CALL message (routine,message_text)

          ALLOCATE( zz_nc(nk,nt), time_nc(nt), zw_nc(nk,nt), zu_nc(nk,nt), zv_nc(nk,nt), z_dt_temp_adv_nc(nk,nt), &
                    z_dt_temp_rad_nc(nk,nt), z_dt_qv_adv_nc(nk,nt), z_dt_u_adv_nc(nk,nt), z_dt_v_adv_nc(nk,nt) )
          zz_nc=0.0_wp;time_nc=0.0_wp;zw_nc=0.0_wp;zu_nc=0.0_wp;zv_nc=0.0_wp;z_dt_temp_adv_nc=0.0_wp
          z_dt_temp_rad_nc=0.0_wp;z_dt_qv_adv_nc=0.0_wp;z_dt_u_adv_nc=0.0_wp;z_dt_v_adv_nc=0.0_wp

          ALLOCATE( w_ls(nlev,nt), u_geo(nlev,nt), v_geo(nlev,nt), ddt_temp_hadv_ls(nlev,nt), &
                    ddt_temp_rad(nlev,nt), ddt_qv_hadv_ls(nlev,nt),&
                    ddt_u_hadv_ls(nlev,nt),ddt_v_hadv_ls(nlev,nt),&
                    tempf(lon,lat,nk,t0),tempf_f(lon,lat,nk,nt),tempf_sf(lon,lat,nt) )
          w_ls=0.0_wp;u_geo=0.0_wp;v_geo=0.0_wp;ddt_temp_hadv_ls=0.0_wp
          ddt_temp_rad=0.0_wp;ddt_qv_hadv_ls=0.0_wp;ddt_u_hadv_ls=0.0_wp;ddt_v_hadv_ls=0.0_wp

          ALLOCATE( bnd_sfc_lat_flx(nt), bnd_sfc_sens_flx(nt), bnd_ts(nt), bnd_qvs(nt),&
                    bnd_Ch(nt), bnd_Cm(nt),bnd_Cq(nt), bnd_ustar(nt),bnd_tg(nt))
          bnd_sfc_lat_flx=0.0_wp; bnd_sfc_sens_flx=0.0_wp; bnd_ts=0.0_wp; bnd_qvs=0.0_wp
          bnd_Ch=0.0_wp; bnd_Cm=0.0_wp;bnd_Cq=0.0_wp; bnd_ustar=0.0_wp;bnd_tg=0.0_wp

          CALL nf (nf90_inq_varid(fileid, 'height_forc', varid), TRIM(routine)//' height')
          CALL nf (nf90_get_var(fileid, varid, tempf_f), routine)
          zz_nc=tempf_f(1,1,nk:1:-1,:)
          !Check if the file is written in descending order
          IF(zz_nc(1,1) < zz_nc(nk,1)) CALL finish ( routine, 'Write LS forcing data in descending order!')

          CALL nf(nf90_inq_varid(fileid, 'time', varid), TRIM(routine)//' time')
          CALL nf(nf90_get_var(fileid, varid, time_nc), routine)
          time_unit = ''               !necessary for formatting
          CALL nf(nf90_get_att(fileid, varid, 'units', time_unit), TRIM(routine)//' units')

          time_unit_short = time_unit(1:INDEX(time_unit,"since")-2)  ! e.g. "minutes since 2020-1-1 00:00:00"
          IF ( get_my_global_mpi_id() == 0 ) THEN
            WRITE(*,*) 'time unit in init_SCM.nc, long: ', TRIM(time_unit), '   short: ', TRIM(time_unit_short), &
              & '     times: ', time_nc
          END IF

          SELECT CASE (time_unit_short)
          CASE ('s', 'sec', 'second', 'seconds')
            time_unit_in_sec = 1
          CASE ('min', 'minute', 'minutes')
            time_unit_in_sec = 60
          CASE ('h', 'hour', 'hours')
            time_unit_in_sec = 3600
          CASE ('d', 'day', 'days')
            time_unit_in_sec = 84600
          CASE DEFAULT
            CALL finish ( routine, 'time unit in init_SCM.nc not s/min/h!')
          END SELECT

          time_nc = time_nc * time_unit_in_sec  ! conversion to [s] from min/hours/days

          IF (nt<=1) THEN 
            dt_forcing = -999._wp
          ELSE
            dt_forcing = time_nc(2) - time_nc(1)
            DO i=2,size(time_nc)-1
              IF ((time_nc(i+1)-time_nc(i)) /= dt_forcing) THEN
                CALL finish(routine,'input timesteps of equal length needed')
              END IF
            END DO
          ENDIF

          nf_status = nf90_inq_varid(fileid, 'wsub', varid)
          nf_status2 = nf90_get_var(fileid, varid, tempf_f)
          !CALL nf (nf90_get_var(fileid, varid, zw_nc), routine)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'wsub not available in init_SCM.nc.  It will be set to 0.')
            zw_nc=0.0_wp
          ELSE
            zw_nc=tempf_f(1,1,nk:1:-1,:)
          END IF
!         write(*,*) 'wLS',zw_nc(:,1)

          CALL nf (nf90_inq_varid(fileid, 'ug', varid), TRIM(routine)//' ug')
          CALL nf (nf90_get_var(fileid, varid, tempf_f), routine)
          zu_nc=tempf_f(1,1,nk:1:-1,:)

          CALL nf (nf90_inq_varid(fileid, 'vg', varid), TRIM(routine)//' vg')
          CALL nf (nf90_get_var(fileid, varid, tempf_f), routine)
          zv_nc=tempf_f(1,1,nk:1:-1,:)

          nf_status = nf90_inq_varid(fileid, 'temp_adv', varid)
          nf_status2 = nf90_get_var(fileid, varid, tempf_f)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'z_dt_temp_adv_nc not available in init_SCM.nc.  It will be set to 0.')
            z_dt_temp_adv_nc=0.0_wp
          ELSE
            z_dt_temp_adv_nc=tempf_f(1,1,nk:1:-1,:)
          END IF
          !write(*,*) 'dTadv',z_dt_temp_adv_nc(:,1)

          nf_status = nf90_inq_varid(fileid, 'temp_rad', varid)
          nf_status2 = nf90_get_var(fileid, varid, tempf_f)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'z_dt_temp_rad_nc not available in init_SCM.nc.  It will be set to 0.')
            z_dt_temp_rad_nc=0.0_wp
          ELSE
            z_dt_temp_rad_nc=tempf_f(1,1,nk:1:-1,:)
          END IF
!         write(*,*) 'dTrad',z_dt_temp_rad_nc(:,1)

          nf_status = nf90_inq_varid(fileid, 'qv_adv', varid)
          nf_status2 = nf90_get_var(fileid, varid, tempf_f)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'z_dt_qv_adv_nc not available in init_SCM.nc It will be set to 0.')
            z_dt_qv_adv_nc=0.0_wp
          ELSE
            z_dt_qv_adv_nc=tempf_f(1,1,nk:1:-1,:)
          END IF
          !CALL nf (nf90_inq_varid(fileid, 'qv_adv', varid), TRIM(routine)//' qv_adv')
          !CALL nf (nf90_get_var(fileid, varid, tempf_f), routine)
          !z_dt_qv_adv_nc=tempf_f(1,1,nk:1:-1,:)
!         write(*,*) 'dQVadv',z_dt_qv_adv_nc(:,1)

          nf_status = nf90_inq_varid(fileid, 'u_adv', varid)
          nf_status2 = nf90_get_var(fileid, varid, tempf_f)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'z_dt_u_adv_nc not available in init_SCM.nc.  It will be set to 0.')
            z_dt_u_adv_nc=0.0_wp
          ELSE
            z_dt_u_adv_nc=tempf_f(1,1,nk:1:-1,:)
          END IF
!         write(*,*) 'dUadv',z_dt_u_adv_nc(:,1)

          nf_status = nf90_inq_varid(fileid, 'v_adv', varid)
          nf_status2 = nf90_get_var(fileid, varid, tempf_f)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'z_dt_v_adv_nc not available in init_SCM.nc.  It will be set to 0.')
            z_dt_v_adv_nc=0.0_wp
          ELSE
            z_dt_v_adv_nc=tempf_f(1,1,nk:1:-1,:)
          END IF
!         write(*,*) 'dVadv',z_dt_v_adv_nc(:,1)

          CALL nf (nf90_inq_varid(fileid, 'sfc_lat_flx', varid), TRIM(routine)//' sfc_lat_flx')
          CALL nf (nf90_get_var(fileid, varid, tempf_sf), routine)
          bnd_sfc_lat_flx=tempf_sf(1,1,:)
!         write(*,*) 'bnd_sfc_lat_flx',bnd_sfc_lat_flx(:)

          CALL nf (nf90_inq_varid(fileid, 'sfc_sens_flx', varid), TRIM(routine)//' sfc_sens_flx')
          CALL nf (nf90_get_var(fileid, varid, tempf_sf), routine)
          bnd_sfc_sens_flx=tempf_sf(1,1,:)
!         write(*,*) 'bnd_sfc_sens_flx',bnd_sfc_sens_flx(:)

          CALL nf (nf90_inq_varid(fileid, 'ts', varid), TRIM(routine)//' ts')
          CALL nf (nf90_get_var(fileid, varid, tempf_sf), routine)
          bnd_ts=tempf_sf(1,1,:)
!         write(*,*) 'bnd_ts',bnd_ts(:)

          nf_status = nf90_inq_varid(fileid, 'tg', varid)
          nf_status2 = nf90_get_var(fileid, varid, tempf_sf)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'bnd_tg not available in init_SCM.nc.  It will be set to 300.')
            bnd_tg=300.0_wp
          ELSE
            bnd_tg=tempf_sf(1,1,:)
          END IF
!         write(*,*) 'bnd_tg',bnd_tg(:)

          nf_status = nf90_inq_varid(fileid, 'qvs', varid)
          nf_status2 = nf90_get_var(fileid, varid, tempf_sf)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'qvs not available in init_SCM.nc.  It will be set to 0.')
            bnd_qvs=0.0_wp
          ELSE
            bnd_qvs=tempf_sf(1,1,:)
          END IF
!         write(*,*) 'bnd_qvs',bnd_qvs(:)

          nf_status = nf90_inq_varid(fileid, 'Ch', varid)
          nf_status2 = nf90_get_var(fileid, varid, tempf_sf)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'Ch not available in init_SCM.nc.  It will be set to 0.')
            bnd_Ch=0.0_wp
          ELSE
            bnd_Ch=tempf_sf(1,1,:)
          END IF
!         write(*,*) 'bnd_Ch',bnd_Ch(:)

          nf_status = nf90_inq_varid(fileid, 'Cm', varid)
          nf_status2 = nf90_get_var(fileid, varid, tempf_sf)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'Cm not available in init_SCM.nc.  It will be set to 0.')
            bnd_Cm=0.0_wp
          ELSE
            bnd_Cm=tempf_sf(1,1,:)
          END IF
!         write(*,*) 'bnd_Cm',bnd_Cm(:)

          nf_status = nf90_inq_varid(fileid, 'Cq', varid)
          nf_status2 = nf90_get_var(fileid, varid, tempf_sf)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'Cq not available in init_SCM.nc.  It will be set to 0.')
            bnd_Cq=0.0_wp
          ELSE
            bnd_Cq=tempf_sf(1,1,:)
          END IF
!         write(*,*) 'bnd_Cq',bnd_Cq(:)

          nf_status  = nf90_inq_varid(fileid, 'ustar', varid)
          nf_status2 = nf90_get_var(fileid, varid, tempf_sf)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'ustar not available in init_SCM.nc.  It will be set to 0.')
            bnd_ustar=0.0_wp
          ELSE
            bnd_ustar=tempf_sf(1,1,:)
          END IF
!         write(*,*) 'bnd_ustar',bnd_ustar(:)

          CALL nf (nf90_close(fileid), routine)


          DO n = 1 , nt
            !Now perform interpolation to grid levels assuming:
            !a) linear interpolation
            !b) Beyond the last Z level the values are linearly extrapolated 
            !c) Assuming model grid is flat-NOT on sphere

            CALL vert_intp_linear_1d(zz_nc(:,n),zw_nc(:,n),p_metrics%z_mc(1,:,1),w_ls(:,n))
            CALL vert_intp_linear_1d(zz_nc(:,n),zu_nc(:,n),p_metrics%z_mc(1,:,1),u_geo(:,n))
            CALL vert_intp_linear_1d(zz_nc(:,n),zv_nc(:,n),p_metrics%z_mc(1,:,1),v_geo(:,n))
            CALL vert_intp_linear_1d(zz_nc(:,n),z_dt_temp_adv_nc(:,n),p_metrics%z_mc(1,:,1),ddt_temp_hadv_ls(:,n))
            CALL vert_intp_linear_1d(zz_nc(:,n),z_dt_temp_rad_nc(:,n),p_metrics%z_mc(1,:,1),ddt_temp_rad(:,n))
            CALL vert_intp_linear_1d(zz_nc(:,n),z_dt_qv_adv_nc(:,n),p_metrics%z_mc(1,:,1),ddt_qv_hadv_ls(:,n))
            CALL vert_intp_linear_1d(zz_nc(:,n),z_dt_u_adv_nc(:,n),p_metrics%z_mc(1,:,1),ddt_u_hadv_ls(:,n))
            CALL vert_intp_linear_1d(zz_nc(:,n),z_dt_v_adv_nc(:,n),p_metrics%z_mc(1,:,1),ddt_v_hadv_ls(:,n))
          END DO
  
          DEALLOCATE( zz_nc, zw_nc, zu_nc, zv_nc, z_dt_temp_adv_nc, z_dt_temp_rad_nc,     &
                    & z_dt_qv_adv_nc, z_dt_u_adv_nc, z_dt_v_adv_nc, time_nc )

!------------------------------------------------------------------------------
! Forcing with normal NETCDF file

        ELSE

          CALL nf (nf90_open('init_SCM.nc', NF90_NOWRITE, fileid), &
            & TRIM(routine)//'   File init_SCM.nc cannot be opened')
    
          CALL nf (nf90_inq_dimid(fileid, 'lev', dimid), TRIM(routine)//' lev')
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = nk), routine)
    
          CALL nf (nf90_inq_dimid(fileid, 'nt', dimid), TRIM(routine)//' nt')
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = nt), routine)

          WRITE(message_text,'(a,i6,a,i6)') 'SCM LS forcing input file: nk, ', nk, ', nt ', nt
          CALL message (routine,message_text)

          ALLOCATE( zz_nc(nk,nt), time_nc(nt), zw_nc(nk,nt), zu_nc(nk,nt), zv_nc(nk,nt), z_dt_temp_adv_nc(nk,nt), &
                    z_dt_temp_rad_nc(nk,nt), z_dt_qv_adv_nc(nk,nt), z_dt_u_adv_nc(nk,nt), z_dt_v_adv_nc(nk,nt) )
          zz_nc=0.0_wp;time_nc=0.0_wp;zw_nc=0.0_wp;zu_nc=0.0_wp;zv_nc=0.0_wp;z_dt_temp_adv_nc=0.0_wp
          z_dt_temp_rad_nc=0.0_wp;z_dt_qv_adv_nc=0.0_wp;z_dt_u_adv_nc=0.0_wp;z_dt_v_adv_nc=0.0_wp
    
          ALLOCATE( w_ls(nlev,nt), u_geo(nlev,nt), v_geo(nlev,nt), ddt_temp_hadv_ls(nlev,nt), &
                    ddt_temp_rad(nlev,nt), ddt_qv_hadv_ls(nlev,nt),&
                    ddt_u_hadv_ls(nlev,nt),ddt_v_hadv_ls(nlev,nt) )
          w_ls=0.0_wp;u_geo=0.0_wp;v_geo=0.0_wp;ddt_temp_hadv_ls=0.0_wp
          ddt_temp_rad=0.0_wp;ddt_qv_hadv_ls=0.0_wp;ddt_u_hadv_ls=0.0_wp;ddt_v_hadv_ls=0.0_wp
    
          ALLOCATE( bnd_sfc_lat_flx(nt), bnd_sfc_sens_flx(nt), bnd_ts(nt), bnd_qvs(nt),&
                    bnd_Ch(nt), bnd_Cm(nt),bnd_Cq(nt), bnd_ustar(nt),bnd_tg(nt))
          bnd_sfc_lat_flx=0.0_wp; bnd_sfc_sens_flx=0.0_wp; bnd_ts=0.0_wp; bnd_qvs=0.0_wp
          bnd_Ch=0.0_wp; bnd_Cm=0.0_wp;bnd_Cq=0.0_wp; bnd_ustar=0.0_wp;bnd_tg=0.0_wp
    
          CALL nf (nf90_inq_varid(fileid, 'height', varid), TRIM(routine)//' height')
          CALL nf (nf90_get_var(fileid, varid, zz_nc), routine)
         !write(*,*) 'zz_nc',zz_nc
          !Check if the file is written in descending order
          IF(zz_nc(1,1) < zz_nc(nk,1)) CALL finish ( routine, 'Write LS forcing data in descending order!')
    
    
          CALL nf(nf90_inq_varid(fileid, 'time', varid), TRIM(routine)//' time')
          CALL nf(nf90_get_var(fileid, varid, time_nc), routine)
          time_unit = ''               !necessary for formatting
          CALL nf(nf90_get_att(fileid, varid, 'units', time_unit), TRIM(routine)//' units')
    
          time_unit_short = time_unit(1:INDEX(time_unit,"since")-2)  ! e.g. "minutes since 2020-1-1 00:00:00"
          IF ( get_my_global_mpi_id() == 0 ) THEN
            WRITE(*,*) 'time unit in init_SCM.nc, long: ', TRIM(time_unit), '   short: ', TRIM(time_unit_short), &
              & '     times: ', time_nc
          END IF
        
          SELECT CASE (time_unit_short)
          CASE ('s', 'sec', 'second', 'seconds')
            time_unit_in_sec = 1
          CASE ('min', 'minute', 'minutes')
            time_unit_in_sec = 60
          CASE ('h', 'hour', 'hours')
            time_unit_in_sec = 3600
          CASE ('d', 'day', 'days')
            time_unit_in_sec = 84600
          CASE DEFAULT
            CALL finish ( routine, 'time unit in init_SCM.nc not s/min/h!')
          END SELECT
    
          time_nc = time_nc * time_unit_in_sec  ! conversion to [s] from min/hours/days
    
          IF (nt<=1) THEN 
            dt_forcing = -999._wp
          ELSE
            dt_forcing = time_nc(2) - time_nc(1)
            DO i=2,size(time_nc)-1
              IF ((time_nc(i+1)-time_nc(i)) /= dt_forcing) THEN
                CALL finish(routine,'input timesteps of equal length needed')
              END IF
            END DO
          ENDIF
    
          CALL nf (nf90_inq_varid(fileid, 'wLS', varid), TRIM(routine)//' wLS')
          CALL nf (nf90_get_var(fileid, varid,zw_nc), routine)
         !write(*,*) 'wLS',zw_nc(:,1)
    
          CALL nf (nf90_inq_varid(fileid, 'uGEO', varid), TRIM(routine)//' uGEO')
          CALL nf (nf90_get_var(fileid, varid,zu_nc), routine)
         !write(*,*) 'uGEO',zu_nc(:,1)
    
          CALL nf (nf90_inq_varid(fileid, 'vGEO', varid), TRIM(routine)//' vGEO')
          CALL nf (nf90_get_var(fileid, varid,zv_nc), routine)
         !write(*,*) 'vGEO',zv_nc(:,1)
    
          CALL nf (nf90_inq_varid(fileid, 'dTadv', varid), TRIM(routine)//' dTadv')
          CALL nf (nf90_get_var(fileid, varid,z_dt_temp_adv_nc), routine)
         !write(*,*) 'dTadv',z_dt_temp_adv_nc(:,1)
    
          CALL nf (nf90_inq_varid(fileid, 'dTrad', varid), TRIM(routine)//' dTrad')
          CALL nf (nf90_get_var(fileid, varid,z_dt_temp_rad_nc), routine)
         !write(*,*) 'dTrad',z_dt_temp_rad_nc(:,1)
    
          CALL nf (nf90_inq_varid(fileid, 'dQVadv', varid), TRIM(routine)//' dQVadv')
          CALL nf (nf90_get_var(fileid, varid,z_dt_qv_adv_nc), routine)
         !write(*,*) 'dQVadv',z_dt_qv_adv_nc(:,1)
    
          CALL nf (nf90_inq_varid(fileid, 'dUadv', varid), TRIM(routine)//' dUadv')
          CALL nf (nf90_get_var(fileid, varid,z_dt_u_adv_nc), routine)
         !write(*,*) 'dUadv',z_dt_u_adv_nc(:,1)
    
          CALL nf (nf90_inq_varid(fileid, 'dVadv', varid), TRIM(routine)//' dVadv')
          CALL nf (nf90_get_var(fileid, varid,z_dt_v_adv_nc), routine)
         !write(*,*) 'dVadv',z_dt_v_adv_nc(:,1)
    
          CALL nf (nf90_inq_varid(fileid, 'sfc_lat_flx', varid), TRIM(routine)//' sfc_lat_flx')
          CALL nf (nf90_get_var(fileid, varid, bnd_sfc_lat_flx), routine)
         !write(*,*) 'bnd_sfc_lat_flx',bnd_sfc_lat_flx(:)
    
          CALL nf (nf90_inq_varid(fileid, 'sfc_sens_flx', varid), TRIM(routine)//' sfc_sens_flx')
          CALL nf (nf90_get_var(fileid, varid, bnd_sfc_sens_flx), routine)
         !write(*,*) 'bnd_sfc_sens_flx',bnd_sfc_sens_flx(:)
    
          CALL nf (nf90_inq_varid(fileid, 'ts', varid), TRIM(routine)//' ts')
          CALL nf (nf90_get_var(fileid, varid, bnd_ts), routine)
         !write(*,*) 'bnd_ts',bnd_ts(:)
    
          CALL nf (nf90_inq_varid(fileid, 'tg', varid), TRIM(routine)//' tg')
          CALL nf (nf90_get_var(fileid, varid, bnd_tg), routine)
         !write(*,*) 'bnd_tg',bnd_tg(:)
    
          CALL nf (nf90_inq_varid(fileid, 'qvs', varid), TRIM(routine)//' qvs')
          CALL nf (nf90_get_var(fileid, varid, bnd_qvs), routine)
         !write(*,*) 'bnd_qvs',bnd_qvs(:)
    
          CALL nf (nf90_inq_varid(fileid, 'Ch', varid), TRIM(routine)//' Ch')
          CALL nf (nf90_get_var(fileid, varid, bnd_Ch), routine)
         !write(*,*) 'bnd_Ch',bnd_Ch(:)
    
          CALL nf (nf90_inq_varid(fileid, 'Cm', varid), TRIM(routine)//' Cm')
          CALL nf (nf90_get_var(fileid, varid, bnd_Cm), routine)
         !write(*,*) 'bnd_Cm',bnd_Cm(:)
    
          CALL nf (nf90_inq_varid(fileid, 'Cq', varid), TRIM(routine)//' Cq')
          CALL nf (nf90_get_var(fileid, varid, bnd_Cq), routine)
         !write(*,*) 'bnd_Cq',bnd_Cq(:)
    
          CALL nf (nf90_inq_varid(fileid, 'ustar', varid), TRIM(routine)//' ustar')
          CALL nf (nf90_get_var(fileid, varid, bnd_ustar), routine)
         !write(*,*) 'bnd_ustar',bnd_ustar(:)
    
          CALL nf (nf90_close(fileid), routine)
    
     
          DO n = 1 , nt
              !Now perform interpolation to grid levels assuming:
              !a) linear interpolation
              !b) Beyond the last Z level the values are linearly extrapolated 
              !c) Assuming model grid is flat-NOT on sphere
    
              CALL vert_intp_linear_1d(zz_nc(:,n),zw_nc(:,n),p_metrics%z_mc(1,:,1),w_ls(:,n))
              CALL vert_intp_linear_1d(zz_nc(:,n),zu_nc(:,n),p_metrics%z_mc(1,:,1),u_geo(:,n))
              CALL vert_intp_linear_1d(zz_nc(:,n),zv_nc(:,n),p_metrics%z_mc(1,:,1),v_geo(:,n))
              CALL vert_intp_linear_1d(zz_nc(:,n),z_dt_temp_adv_nc(:,n),p_metrics%z_mc(1,:,1),ddt_temp_hadv_ls(:,n))
              CALL vert_intp_linear_1d(zz_nc(:,n),z_dt_temp_rad_nc(:,n),p_metrics%z_mc(1,:,1),ddt_temp_rad(:,n))
              CALL vert_intp_linear_1d(zz_nc(:,n),z_dt_qv_adv_nc(:,n),p_metrics%z_mc(1,:,1),ddt_qv_hadv_ls(:,n))
              CALL vert_intp_linear_1d(zz_nc(:,n),z_dt_u_adv_nc(:,n),p_metrics%z_mc(1,:,1),ddt_u_hadv_ls(:,n))
              CALL vert_intp_linear_1d(zz_nc(:,n),z_dt_v_adv_nc(:,n),p_metrics%z_mc(1,:,1),ddt_v_hadv_ls(:,n))
          END DO
    
    
          DEALLOCATE( zz_nc, zw_nc, zu_nc, zv_nc, z_dt_temp_adv_nc, z_dt_temp_rad_nc,     &
                    & z_dt_qv_adv_nc, z_dt_u_adv_nc, z_dt_v_adv_nc, time_nc )
    
        ENDIF

!------------------------------------------------------------------------------
! Forcing with ASCII file

      ELSE

        !Open formatted file to read ls forcing data
        iunit = find_next_free_unit(10,20)
        OPEN (unit=iunit,file='ls_forcing.dat',access='SEQUENTIAL', &
              form='FORMATTED', action='READ', status='OLD', IOSTAT=ist)
  
        IF (ist/=success) CALL finish(routine, 'open ls_forcing.dat failed')
  
        !Read the input file till end. The order of file assumed is:
        !Z(m) - u_geo(m/s) - v_geo(m/s) - w_ls(m/s) - dt_temp_rad(K/s) - ddt_qv_hadv_ls(1/s) - ddt_temp_hadv_ls(K/s)
  
        !Skip the first line and read next 2 lines with information about vertical and time levels
        READ(iunit,*,IOSTAT=ist)                       !skip
        IF (ist == success) READ(iunit,*,IOSTAT=ist)nk,dt_forcing,end_time !vertical levels, forcing interval, end of forcing data
        IF(ist/=success) &
          CALL finish(routine, &
            'Must provide vertical level, forcing interval, end of forcing time info in forcing file')
  
        nt = INT(end_time/dt_forcing)+1
  
        IF (nt>1) THEN
          READ(iunit,*,IOSTAT=ist)nskip !lines to skip between successive time levels
          IF (ist/=success) &
            CALL finish(routine, &
              'Time levels > 1 so must provide number lines to skip between successive time levels in forcing file')
        END IF
  
        IF (nt<=1) dt_forcing=-999._wp
  
        ALLOCATE( zz(nk), zw(nk), zu(nk), zv(nk), z_dt_temp_adv(nk),                        &
                & z_dt_temp_rad(nk), z_dt_qv_adv(nk), z_dt_u_adv(nk), z_dt_v_adv(nk) )
  
        ALLOCATE( w_ls(nlev,nt), u_geo(nlev,nt), v_geo(nlev,nt), ddt_temp_hadv_ls(nlev,nt), &
                & ddt_temp_rad(nlev,nt), ddt_qv_hadv_ls(nlev,nt),                           &
                & ddt_u_hadv_ls(nlev,nt),ddt_v_hadv_ls(nlev,nt) )
  
        DO n = 1,nt
  
          DO jk = nk , 1, -1
            READ(iunit,*,IOSTAT=ist)zz(jk),zu(jk),zv(jk),zw(jk),z_dt_temp_rad(jk), &
                & z_dt_qv_adv(jk),z_dt_temp_adv(jk),z_dt_u_adv(jk),z_dt_v_adv(jk)
            IF(ist/=success) CALL finish(routine, 'something wrong in forcing.dat')
          END DO
  
          !Skip lines
          IF(nt>1)THEN
            DO i = 1 , nskip
              READ(iunit,*,IOSTAT=ist)
            END DO
          END IF
  
          !Check if the file is written in descending order
          IF (zz(1) < zz(nk)) &
            CALL finish(routine, 'Write LS forcing data in descending order!')
  
          !Now perform interpolation to grid levels assuming:
          !a) linear interpolation
          !b) Beyond the last Z level the values are linearly extrapolated
          !c) Assuming model grid is flat-NOT on sphere
  
          CALL vert_intp_linear_1d(zz,zw,p_metrics%z_mc(1,:,1),w_ls(:,n))
          CALL vert_intp_linear_1d(zz,zu,p_metrics%z_mc(1,:,1),u_geo(:,n))
          CALL vert_intp_linear_1d(zz,zv,p_metrics%z_mc(1,:,1),v_geo(:,n))
          CALL vert_intp_linear_1d(zz,z_dt_temp_adv,p_metrics%z_mc(1,:,1),ddt_temp_hadv_ls(:,n))
          CALL vert_intp_linear_1d(zz,z_dt_temp_rad,p_metrics%z_mc(1,:,1),ddt_temp_rad(:,n))
          CALL vert_intp_linear_1d(zz,z_dt_qv_adv,p_metrics%z_mc(1,:,1),ddt_qv_hadv_ls(:,n))
          CALL vert_intp_linear_1d(zz,z_dt_u_adv,p_metrics%z_mc(1,:,1),ddt_u_hadv_ls(:,n))
          CALL vert_intp_linear_1d(zz,z_dt_v_adv,p_metrics%z_mc(1,:,1),ddt_v_hadv_ls(:,n))
  
        END DO !n
  
        CLOSE(iunit)
  
        DEALLOCATE( zz, zw, zu, zv, z_dt_temp_adv, z_dt_temp_rad, z_dt_qv_adv,z_dt_u_adv,z_dt_v_adv )

      ENDIF

      WRITE(message_text,*)dt_forcing
      CALL message('Time varying LS forcing read in:',message_text)

    END IF !is_ls_forcing


!------------------------------------------------------------------------------
! READ NUDGING FILE
!------------------------------------------------------------------------------

    IF(is_nudging) THEN

      IF(i_scm_netcdf > 0) THEN

        IF(i_scm_netcdf == 2) THEN

!------------------------------------------------------------------------------
! Nudging with unified NETCDF file

          CALL nf (nf90_open('init_SCM.nc', NF90_NOWRITE, fileid), &
            & TRIM(routine)//'   File init_SCM.nc cannot be opened')
   
          CALL nf (nf90_inq_dimid(fileid, 'lev', dimid), TRIM(routine)//' lev')
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = nk), routine)
  
          CALL nf (nf90_inq_dimid(fileid, 'lat', dimid), routine)
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = lat), routine)
  
          CALL nf (nf90_inq_dimid(fileid, 'lon', dimid), routine)
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = lon), routine)
 
          CALL nf (nf90_inq_dimid(fileid, 't0', dimid), routine)
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = t0), routine)
    
          CALL nf (nf90_inq_dimid(fileid, 'nt', dimid), TRIM(routine)//' nt')
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = nt), routine)

          WRITE(message_text,'(a,i5,a,i5)') 'SCM nudging input: nk', nk, '  nt', nt
          CALL message (routine,message_text)

          ALLOCATE( zz_nc(nk,nt), time_nc(nt), zu_nc(nk,nt), zv_nc(nk,nt), ztheta_nc(nk,nt), z_qv_nc(nk,nt),&
          &tempf(lon,lat,nk,nt) )
          zz_nc=0.0_wp ; time_nc=0.0_wp ; zu_nc=0.0_wp ; zv_nc=0.0_wp ; ztheta_nc=0.0_wp ; z_qv_nc=0.0_wp
    
          ALLOCATE( u_nudg(nlev,nt), v_nudg(nlev,nt), theta_nudg(nlev,nt), qv_nudg(nlev,nt) )
          u_nudg=0.0_wp ; v_nudg=0.0_wp ; theta_nudg=0.0_wp ; qv_nudg=0.0_wp
    
          CALL nf (nf90_inq_varid(fileid, 'height', varid), TRIM(routine)//' height')
          CALL nf (nf90_get_var(fileid, varid, zz_nc), routine)
!         write(*,*) 'zz_nc',zz_nc
          !Check if the file is written in descending order
          IF(zz_nc(1,1) < zz_nc(nk,1)) CALL finish ( routine, 'Write LS forcing data in descending order!')
    
          CALL nf(nf90_inq_varid(fileid, 'time', varid), TRIM(routine)//' time')
          CALL nf(nf90_get_var(fileid, varid, time_nc), routine)
          time_unit = ''               !necessary for formatting
          CALL nf(nf90_get_att(fileid, varid, 'units', time_unit), TRIM(routine)//' units')
    
          time_unit_short = time_unit(1:INDEX(time_unit,"since")-2)  ! e.g. "minutes since 2020-1-1 00:00:00"
          IF ( get_my_global_mpi_id() == 0 ) THEN
            WRITE(*,*) 'time unit in init_SCM.nc, long: ', TRIM(time_unit), '   short: ', TRIM(time_unit_short), &
              & '     times: ', time_nc
          END IF
        
          SELECT CASE (time_unit_short)
          CASE ('s', 'sec', 'second', 'seconds')
            time_unit_in_sec = 1
          CASE ('min', 'minute', 'minutes')
            time_unit_in_sec = 60
          CASE ('h', 'hour', 'hours')
            time_unit_in_sec = 3600
          CASE ('d', 'day', 'days')
            time_unit_in_sec = 84600
          CASE DEFAULT
            CALL finish ( routine, 'time unit in init_SCM.nc not s/min/h!')
          END SELECT
   
          time_nc = time_nc * time_unit_in_sec  ! conversion to [s] from min/hours/days
   
          IF (nt<=1) THEN 
            dt_nudging = -999._wp
          ELSE
            dt_nudging = time_nc(2) - time_nc(1)
            DO i=2,size(time_nc)-1
              IF ((time_nc(i+1)-time_nc(i)) /= dt_nudging) THEN
                CALL finish(routine,'input timesteps of equal length needed')
              END IF
            END DO
          ENDIF
     
          IF ( dt_relax <= 0.0_wp ) THEN   ! dt_relax not defined by namelist input
            CALL nf (nf90_inq_varid(fileid, 'dt_relax', varid), TRIM(routine)//' dt_relax')
            CALL nf (nf90_get_var(fileid, varid, temp_nf), routine)
            dt_relax = temp_nf(1)
            write(*,*) 'dt_relax [s]: ', dt_relax
          END IF
    
          nf_status = nf90_inq_varid(fileid, 'u_nudging', varid)
          nf_status2 = nf90_get_var(fileid, varid, tempf)
          IF (nf_status /= nf90_noerr) THEN
           write(*,*) 'u_nudging not available in init_SCM.nc.  It will be set to 0.'
           zu_nc=0.0_wp
          ELSE
           zu_nc=tempf(1,1,:,:)
          END IF
 
          nf_status = nf90_inq_varid(fileid, 'v_nudging',  varid)
          nf_status2 = nf90_get_var(fileid, varid, tempf)
          IF (nf_status /= nf90_noerr) THEN
           CALL message(routine,'v_nudging not available in init_SCM.nc.  It will be set to 0.')
           zv_nc=0.0_wp
          ELSE
           zv_nc=tempf(1,1,:,:)
          END IF
  
          nf_status = nf90_inq_varid(fileid, 'theta_nudging', varid)
          nf_status2 = nf90_get_var(fileid, varid, tempf)
          IF (nf_status /= nf90_noerr) THEN
           CALL message(routine,'theta_nudging not available in init_SCM.nc.  It will be set to 0.')
           ztheta_nc=0.0_wp
          ELSE
           ztheta_nc=tempf(1,1,:,:)
          END IF
   
          nf_status = nf90_inq_varid(fileid, 'qv_nudging', varid)
          nf_status2 = nf90_get_var(fileid, varid, tempf)
          IF (nf_status /= nf90_noerr) THEN
           CALL message(routine,'qv_nudging not available in init_SCM.nc.  It will be set to 0.')
           z_qv_nc=0.0_wp
          ELSE
           z_qv_nc=tempf(1,1,:,:)
          END IF
    
!         write(*,*) 'uLS',zu_nc(:,1)
!         write(*,*) 'vLS',zv_nc(:,1)
!         write(*,*) 'qvLS',z_qv_nc(:,1)
!         write(*,*) 'thLS',ztheta_nc(:,1)
    
          CALL nf (nf90_close(fileid), routine)
    
          DO n = 1 , nt
            !Now perform interpolation to grid levels assuming:
            !a) linear interpolation
            !b) Beyond the last Z level the values are linearly extrapolated 
            !c) Assuming model grid is flat-NOT on sphere
    
            CALL vert_intp_linear_1d(zz_nc(:,n),zu_nc(:,n)    ,p_metrics%z_mc(1,:,1),u_nudg(:,n))
            CALL vert_intp_linear_1d(zz_nc(:,n),zv_nc(:,n)    ,p_metrics%z_mc(1,:,1),v_nudg(:,n))
            CALL vert_intp_linear_1d(zz_nc(:,n),ztheta_nc(:,n),p_metrics%z_mc(1,:,1),theta_nudg(:,n))
            CALL vert_intp_linear_1d(zz_nc(:,n),z_qv_nc(:,n)  ,p_metrics%z_mc(1,:,1),qv_nudg(:,n))
          END DO
    
          DEALLOCATE( zz_nc, time_nc, zu_nc, zv_nc, ztheta_nc, z_qv_nc )

!------------------------------------------------------------------------------
! Nudging with normal NETCDF file

        ELSE

          CALL nf (nf90_open('init_SCM.nc', NF90_NOWRITE, fileid), &
            & TRIM(routine)//'   File init_SCM.nc cannot be opened')
  
          CALL nf (nf90_inq_dimid(fileid, 'lev', dimid), TRIM(routine)//' lev')
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = nk), routine)
   
          CALL nf (nf90_inq_dimid(fileid, 'nt', dimid), TRIM(routine)//' nt')
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = nt), routine)

          WRITE(message_text,'(a,i5,a,i5)') 'SCM nudging input: nk', nk, '  nt', nt
          CALL message(routine,message_text)
    
          ALLOCATE( zz_nc(nk,nt), time_nc(nt), zu_nc(nk,nt), zv_nc(nk,nt), ztheta_nc(nk,nt), z_qv_nc(nk,nt) )
          zz_nc=0.0_wp ; time_nc=0.0_wp ; zu_nc=0.0_wp ; zv_nc=0.0_wp ; ztheta_nc=0.0_wp ; z_qv_nc=0.0_wp
    
          ALLOCATE( u_nudg(nlev,nt), v_nudg(nlev,nt), theta_nudg(nlev,nt), qv_nudg(nlev,nt) )
          u_nudg=0.0_wp ; v_nudg=0.0_wp ; theta_nudg=0.0_wp ; qv_nudg=0.0_wp
    
          CALL nf (nf90_inq_varid(fileid, 'height', varid), TRIM(routine)//' height')
          CALL nf (nf90_get_var(fileid, varid, zz_nc), routine)
!         write(*,*) 'zz_nc',zz_nc
          !Check if the file is written in descending order
          IF(zz_nc(1,1) < zz_nc(nk,1)) CALL finish ( routine, 'Write LS forcing data in descending order!')
    
          CALL nf(nf90_inq_varid(fileid, 'time', varid), TRIM(routine)//' time')
          CALL nf(nf90_get_var(fileid, varid, time_nc), routine)
          time_unit = ''               !necessary for formatting
          CALL nf(nf90_get_att(fileid, varid, 'units', time_unit), TRIM(routine)//' units')
    
          time_unit_short = time_unit(1:INDEX(time_unit,"since")-2)  ! e.g. "minutes since 2020-1-1 00:00:00"
          IF ( get_my_global_mpi_id() == 0 ) THEN
            WRITE(*,*) 'time unit in init_SCM.nc, long: ', TRIM(time_unit), '   short: ', TRIM(time_unit_short), &
              & '     times: ', time_nc
          END IF
        
          SELECT CASE (time_unit_short)
          CASE ('s', 'sec', 'second', 'seconds')
            time_unit_in_sec = 1
          CASE ('min', 'minute', 'minutes')
            time_unit_in_sec = 60
          CASE ('h', 'hour', 'hours')
            time_unit_in_sec = 3600
          CASE ('d', 'day', 'days')
            time_unit_in_sec = 84600
          CASE DEFAULT
            CALL finish ( routine, 'time unit in init_SCM.nc not s/min/h!')
          END SELECT
    
          time_nc = time_nc * time_unit_in_sec  ! conversion to [s] from min/hours/days
  
          IF (nt<=1) THEN 
            dt_nudging = -999._wp
          ELSE
            dt_nudging = time_nc(2) - time_nc(1)
            DO i=2,size(time_nc)-1
              IF ((time_nc(i+1)-time_nc(i)) /= dt_nudging) THEN
                CALL finish(routine,'input timesteps of equal length needed')
              END IF
            END DO
          ENDIF
   
          IF ( dt_relax <= 0.0_wp ) THEN   ! dt_relax not defined by namelist input
            CALL nf (nf90_inq_varid(fileid, 'dt_relax', varid), TRIM(routine)//' dt_relax')
            CALL nf (nf90_get_var(fileid, varid, temp_nf), routine)
            dt_relax = temp_nf(1)
            write(*,*) 'dt_relax [s]: ', dt_relax
          END IF
    
          CALL nf (nf90_inq_varid(fileid, 'uIN',  varid), TRIM(routine)//' uIN')
          CALL nf (nf90_get_var(fileid, varid, zu_nc), routine)
    
          CALL nf (nf90_inq_varid(fileid, 'vIN',  varid), TRIM(routine)//' vIN')
          CALL nf (nf90_get_var(fileid, varid, zv_nc), routine)
    
          CALL nf (nf90_inq_varid(fileid, 'qvIN', varid), TRIM(routine)//' qvIN')
          CALL nf (nf90_get_var(fileid, varid, z_qv_nc), routine)
    
          CALL nf (nf90_inq_varid(fileid, 'thIN', varid), TRIM(routine)//' thIN')
          CALL nf (nf90_get_var(fileid, varid, ztheta_nc), routine)
    
!         write(*,*) 'uLS',zu_nc(:,1)
!         write(*,*) 'vLS',zv_nc(:,1)
!         write(*,*) 'qvLS',z_qv_nc(:,1)
!         write(*,*) 'thLS',ztheta_nc(:,1)
    
          CALL nf (nf90_close(fileid), routine)
    
          DO n = 1 , nt
            !Now perform interpolation to grid levels assuming:
            !a) linear interpolation
            !b) Beyond the last Z level the values are linearly extrapolated 
            !c) Assuming model grid is flat-NOT on sphere
    
            CALL vert_intp_linear_1d(zz_nc(:,n),zu_nc(:,n)    ,p_metrics%z_mc(1,:,1),u_nudg(:,n))
            CALL vert_intp_linear_1d(zz_nc(:,n),zv_nc(:,n)    ,p_metrics%z_mc(1,:,1),v_nudg(:,n))
            CALL vert_intp_linear_1d(zz_nc(:,n),ztheta_nc(:,n),p_metrics%z_mc(1,:,1),theta_nudg(:,n))
            CALL vert_intp_linear_1d(zz_nc(:,n),z_qv_nc(:,n)  ,p_metrics%z_mc(1,:,1),qv_nudg(:,n))
          END DO
    
          DEALLOCATE( zz_nc, time_nc, zu_nc, zv_nc, ztheta_nc, z_qv_nc )
  
        ENDIF

!------------------------------------------------------------------------------
! Nudging with formatted ASCII file

      ELSE

        !Open formatted file to read nudging data
        !The order of file assumed is: Z(m) - tau(s) - u(m/s) - v(m/s) - pot. temp(K) - qv(kg/kg) 
        iunit = find_next_free_unit(10,20)
        OPEN (unit=iunit,file='nudging.dat',access='SEQUENTIAL', &
              form='FORMATTED', action='READ', status='OLD', IOSTAT=ist)
 
        IF(ist/=success)THEN
          CALL finish (routine, 'open nudging.dat failed')
        ENDIF
 
 
        !Skip the first line and read next 2 lines with information about vertical and time levels
        READ(iunit,*,IOSTAT=ist)                         !skip
        READ(iunit,*,IOSTAT=ist)nk,dt_nudging,end_time   !vertical levels,time levels and interval
        IF (ist/=success) CALL finish(routine, &
          'Must provide vertical level, time level and time interval info in nudging file')
 
        IF (nt>1) THEN
          READ(iunit,*,IOSTAT=ist)nskip              !lines to skip between successive time levels
          IF (ist/=success) CALL finish(routine, &
            'Time levels > 1 so must provide number lines to skip between successive time levels in nudging file')
        END IF
        IF (nt<=1) dt_nudging=-999._wp
 
        ALLOCATE( zz(nk), zu(nk), zv(nk), ztheta(nk), z_qv(nk) )
 
        ALLOCATE( u_nudg(nlev,nt), v_nudg(nlev,nt), theta_nudg(nlev,nt), qv_nudg(nlev,nt))

        DO n = 1 , nt

          DO jk = nk , 1, -1
            IF ( dt_relax <= 0.0_wp ) THEN   ! dt_relax not defined by namelist input
              READ(iunit,*,IOSTAT=ist)zz(jk),dt_relax,zu(jk),zv(jk),ztheta(jk),z_qv(jk)
              IF (ist/=success) CALL finish(routine, 'something wrong in nudging.dat')
            END IF
          END DO
 
          !Skip lines
          IF(nt>1)THEN
            DO i = 1 , nskip
              READ(iunit,*,IOSTAT=ist)
            END DO
          END IF
 
          !Check if the file is written in descending order
          IF (zz(1) < zz(nk)) &
            CALL finish(routine, 'Write nuding data in descending order!')
 
          !Now perform interpolation to grid levels assuming:
          !a) linear interpolation
          !b) Beyond the last Z level the values are linearly extrapolated
          !c) Assuming model grid is flat-NOT on sphere
 
          CALL vert_intp_linear_1d(zz,zu    ,p_metrics%z_mc(1,:,1),u_nudg(:,n))
          CALL vert_intp_linear_1d(zz,zv    ,p_metrics%z_mc(1,:,1),v_nudg(:,n))
          CALL vert_intp_linear_1d(zz,ztheta,p_metrics%z_mc(1,:,1),theta_nudg(:,n))
          CALL vert_intp_linear_1d(zz,z_qv  ,p_metrics%z_mc(1,:,1),qv_nudg(:,n))
 
        END DO
 
        CLOSE(iunit)
 
        DEALLOCATE( zz, zu, zv, ztheta, z_qv )

      END IF

      WRITE(message_text,*) dt_nudging
      CALL message('Time varying nudging read in with time step:', message_text)

    END IF !is_nudging


  END SUBROUTINE init_ls_forcing

  !>
  !! apply_ls_forcing
  !!------------------------------------------------------------------------
  !! Apply large-scale forcing: called in the end from mo_nh_interface_nwp
  !! It uses the most updated u,v to calculate the large-scale subsidence induced
  !! advective tendencies. All other tendencies don't need any computation.
  !! All tendencies are then accumulated in the slow physics tendency terms
  !!
  !!------------------------------------------------------------------------

  SUBROUTINE apply_ls_forcing(p_patch, p_metrics, curr_sim_time,       & !in
                              p_prog, p_diag, qv, rl_start, rl_end,    & !in
                              ddt_u_ls, ddt_v_ls,                      & !out
                              ddt_temp_ls, ddt_qv_ls,                  & !out
                              ddt_temp_subs_ls, ddt_qv_subs_ls,        & !output
                              ddt_temp_adv_ls, ddt_qv_adv_ls,          & !output
                              ddt_u_adv_ls, ddt_v_adv_ls, wsub,        & !output
                              fc_sfc_lat_flx,fc_sfc_sens_flx,fc_ts,    & !output
                              fc_tg, fc_qvs,fc_Ch,fc_Cq,fc_Cm,fc_ustar,& !output
                              temp_nudge, u_nudge, v_nudge, q_nudge)     !output

    TYPE(t_patch),INTENT(in),        TARGET :: p_patch
    TYPE(t_nh_metrics),INTENT(in),   TARGET :: p_metrics
    TYPE(t_nh_prog),INTENT(in),      TARGET :: p_prog
    TYPE(t_nh_diag),INTENT(in),      TARGET :: p_diag
    INTEGER,        INTENT(in)              :: rl_start
    INTEGER,        INTENT(in)              :: rl_end
    REAL(wp),       INTENT(in)              :: curr_sim_time
    REAL(wp),       INTENT(in)              :: qv(:,:,:)
    REAL(wp),       INTENT(out)             :: ddt_u_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_v_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_temp_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_qv_ls(:)

    !Individual LS tendencies for profile output
    REAL(wp),       INTENT(out)             :: ddt_temp_subs_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_qv_subs_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_temp_adv_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_qv_adv_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_u_adv_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_v_adv_ls(:)
    REAL(wp),       INTENT(out)             :: wsub(:)
    REAL(wp),       INTENT(out)             :: fc_sfc_lat_flx
    REAL(wp),       INTENT(out)             :: fc_sfc_sens_flx
    REAL(wp),       INTENT(out)             :: fc_ts              !surface temperature
    REAL(wp),       INTENT(out)             :: fc_tg              !skin temperature
    REAL(wp),       INTENT(out)             :: fc_qvs
    REAL(wp),       INTENT(out)             :: fc_Ch
    REAL(wp),       INTENT(out)             :: fc_Cq
    REAL(wp),       INTENT(out)             :: fc_Cm
    REAL(wp),       INTENT(out)             :: fc_ustar
    REAL(wp),       INTENT(out), OPTIONAL   :: temp_nudge(:)
    REAL(wp),       INTENT(out), OPTIONAL   :: u_nudge(:)
    REAL(wp),       INTENT(out), OPTIONAL   :: v_nudge(:)
    REAL(wp),       INTENT(out), OPTIONAL   :: q_nudge(:,:)


    CHARACTER(len=*), PARAMETER :: routine = 'mo_ls_forcing:apply_ls_forcing'
    REAL(wp) :: thetain(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp) :: inv_no_gb_cells
    REAL(wp), DIMENSION(p_patch%nlev)  :: u_gb, v_gb, theta_gb, qv_gb
    REAL(wp), DIMENSION(p_patch%nlev+1):: u_gb_hl, v_gb_hl, theta_gb_hl, qv_gb_hl
    REAL(wp), DIMENSION(p_patch%nlev)  :: inv_dz, exner_gb
    REAL(wp), DIMENSION(p_patch%nlev)  :: z_ddt_t_rad, z_ugeo, z_vgeo

    REAL(wp) :: int_weight
    INTEGER  :: i_startblk, i_endblk, jk, nlev, jb, jc
    INTEGER  :: n_curr, n_next, i_startidx, i_endidx

    nlev      = p_patch%nlev
    inv_no_gb_cells = 1._wp / REAL(p_patch%n_patch_cells_g,wp)

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    !0) Initialize all passed ddt's to 0
!$OMP PARALLEL
    CALL init(ddt_u_ls, lacc=.FALSE.)
    CALL init(ddt_v_ls, lacc=.FALSE.)
    CALL init(ddt_temp_ls, lacc=.FALSE.)
    call init(ddt_qv_ls, lacc=.FALSE.)

    !use theta instead of temperature for subsidence (Anurag, Christopher)
!$OMP BARRIER
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
       DO jk = 1 , nlev
         DO jc = i_startidx, i_endidx
           thetain(jc,jk,jb)  = p_diag%temp(jc,jk,jb)/p_prog%exner(jc,jk,jb)
         END DO
       END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    CALL global_hor_mean(p_patch, p_diag%u,            u_gb,     inv_no_gb_cells)
    CALL global_hor_mean(p_patch, p_diag%v,            v_gb,     inv_no_gb_cells)
    CALL global_hor_mean(p_patch, thetain,             theta_gb, inv_no_gb_cells)
    CALL global_hor_mean(p_patch, qv,                  qv_gb,    inv_no_gb_cells)
    CALL global_hor_mean(p_patch, p_prog%exner(:,:,:), exner_gb, inv_no_gb_cells)

    IF (is_subsidence_moment.OR.is_subsidence_heat) THEN
      inv_dz(:) = 1._wp / p_metrics%ddqz_z_full(1,:,1)
    END IF

    !Find where in the array of LS forcing and nudging current time stands
    !and do linear interpolation in time
    IF (dt_forcing > 0._wp) THEN
      n_curr     = FLOOR(curr_sim_time/dt_forcing)+1
      n_next     = n_curr+1
      int_weight = curr_sim_time/dt_forcing-n_curr+1
      !ibd, exception for last time step
      if(int_weight.eq.0._wp) then
        n_next = n_curr
      endif

      ! Christopher Moseley: temporary catch
      IF (int_weight.LT.0 .OR.int_weight.GT.1) THEN
        WRITE(message_text,'(a,2(i0,", "),g0)') &
             'INTERPOLATION ERROR in LS forcing: ', &
             n_curr,n_next,int_weight
        CALL message(routine, message_text)
      END IF

      wsub            = w_ls            (:,n_curr)*(1.-int_weight) + w_ls            (:,n_next)*int_weight
      ddt_qv_adv_ls   = ddt_qv_hadv_ls  (:,n_curr)*(1.-int_weight) + ddt_qv_hadv_ls  (:,n_next)*int_weight
      ddt_u_adv_ls    = ddt_u_hadv_ls   (:,n_curr)*(1.-int_weight) + ddt_u_hadv_ls   (:,n_next)*int_weight
      ddt_v_adv_ls    = ddt_v_hadv_ls   (:,n_curr)*(1.-int_weight) + ddt_v_hadv_ls   (:,n_next)*int_weight
      ddt_temp_adv_ls = ddt_temp_hadv_ls(:,n_curr)*(1.-int_weight) + ddt_temp_hadv_ls(:,n_next)*int_weight
      z_ddt_t_rad     = ddt_temp_rad    (:,n_curr)*(1.-int_weight) + ddt_temp_rad    (:,n_next)*int_weight
      z_ugeo          = u_geo           (:,n_curr)*(1.-int_weight) + u_geo           (:,n_next)*int_weight
      z_vgeo          = v_geo           (:,n_curr)*(1.-int_weight) + v_geo           (:,n_next)*int_weight
      fc_sfc_lat_flx  = bnd_sfc_lat_flx (n_curr)  *(1.-int_weight) + bnd_sfc_lat_flx (n_next)  *int_weight
      fc_sfc_sens_flx = bnd_sfc_sens_flx(n_curr)  *(1.-int_weight) + bnd_sfc_sens_flx(n_next)  *int_weight
      fc_ts           = bnd_ts          (n_curr)  *(1.-int_weight) + bnd_ts          (n_next)  *int_weight
      fc_tg           = bnd_tg          (n_curr)  *(1.-int_weight) + bnd_tg          (n_next)  *int_weight
      fc_qvs          = bnd_qvs         (n_curr)  *(1.-int_weight) + bnd_qvs         (n_next)  *int_weight
      fc_Ch           = bnd_Ch          (n_curr)  *(1.-int_weight) + bnd_Ch          (n_next)  *int_weight
      fc_Cm           = bnd_Cm          (n_curr)  *(1.-int_weight) + bnd_Cm          (n_next)  *int_weight
      fc_Cq           = bnd_Cq          (n_curr)  *(1.-int_weight) + bnd_Cq          (n_next)  *int_weight
      fc_ustar        = bnd_ustar       (n_curr)  *(1.-int_weight) + bnd_ustar       (n_next)  *int_weight
    ELSE
      wsub            = w_ls(:,1)
      ddt_qv_adv_ls   = ddt_qv_hadv_ls(:,1)
      ddt_temp_adv_ls = ddt_temp_hadv_ls(:,1)
      ddt_u_adv_ls    = ddt_u_hadv_ls(:,1)
      ddt_v_adv_ls    = ddt_v_hadv_ls(:,1)
      z_ddt_t_rad     = ddt_temp_rad(:,1)
      z_ugeo          = u_geo(:,1)
      z_vgeo          = v_geo(:,1)
      fc_sfc_lat_flx  = bnd_sfc_lat_flx(1)
      fc_sfc_sens_flx = bnd_sfc_sens_flx(1)
      fc_ts           = bnd_ts(1)
      fc_tg           = bnd_tg(1)
      fc_qvs          = bnd_qvs(1)
      fc_Ch           = bnd_Ch(1)
      fc_Cm           = bnd_Cm(1)
      fc_Cq           = bnd_Cq(1)
      fc_ustar        = bnd_ustar(1)
    END IF

    !1a) Horizontal interpolation and their advective tendency - momentum

    IF(is_subsidence_moment)THEN
      CALL vert_intp_linear_1d(p_metrics%z_mc(1,:,1),u_gb,p_metrics%z_ifc(1,:,1),u_gb_hl)
      CALL vert_intp_linear_1d(p_metrics%z_mc(1,:,1),v_gb,p_metrics%z_ifc(1,:,1),v_gb_hl)

      ddt_u_ls =  ddt_u_ls - wsub*vertical_derivative(u_gb_hl,inv_dz)
      ddt_v_ls =  ddt_v_ls - wsub*vertical_derivative(v_gb_hl,inv_dz)
    END IF

    !1b) Horizontal interpolation and their advective tendency - temperature and water vapor

    IF(is_subsidence_heat)THEN
      CALL vert_intp_linear_1d(p_metrics%z_mc(1,:,1),theta_gb,p_metrics%z_ifc(1,:,1),theta_gb_hl)
      CALL vert_intp_linear_1d(p_metrics%z_mc(1,:,1),qv_gb   ,p_metrics%z_ifc(1,:,1),qv_gb_hl)

      ddt_temp_subs_ls = -wsub*vertical_derivative(theta_gb_hl,inv_dz)  ! vertical advective tendency for theta
      ddt_qv_subs_ls   = -wsub*vertical_derivative(qv_gb_hl   ,inv_dz)  ! vertical advective tendenciy for qv

      ddt_temp_ls = ddt_temp_ls + ddt_temp_subs_ls
      ddt_qv_ls   = ddt_qv_ls   + ddt_qv_subs_ls
    END IF

    !2) Horizontal advective forcing: at present only for tracers and temperature
    IF(is_advection)THEN
      IF(is_advection_tq)THEN
        ddt_qv_ls   = ddt_qv_ls   + ddt_qv_adv_ls
        ddt_temp_ls = ddt_temp_ls + ddt_temp_adv_ls
      ENDIF
      IF(is_advection_uv)THEN
        ddt_u_ls = ddt_u_ls + ddt_u_adv_ls
        ddt_v_ls = ddt_v_ls + ddt_v_adv_ls
      ENDIF
    END IF

    !3) Coriolis and geostrophic wind
    IF(is_geowind)THEN
      !Remember model grid is flat
      ddt_u_ls = ddt_u_ls - p_patch%cells%f_c(1,1) * z_vgeo
      ddt_v_ls = ddt_v_ls + p_patch%cells%f_c(1,1) * z_ugeo
      !do Coriolis here
      !ddt_u_ls = ddt_u_ls - p_patch%cells%f_c(1,1) * (z_vgeo- p_diag%u(1,:,1))
      !ddt_v_ls = ddt_v_ls + p_patch%cells%f_c(1,1) * (z_ugeo-p_diag%v(1,:,1))
    END IF

    !4) Radiative forcing
    IF(is_rad_forcing)THEN
      ddt_temp_ls = ddt_temp_ls + z_ddt_t_rad
    END IF

    !5) Nudging
    IF(is_nudging) THEN

      WRITE(message_text,*)dt_nudging
      CALL message('dt_nudging (time-step input data) =',message_text)

      ! multiple time steps available separated by dt_nudging
      IF (dt_nudging > 0._wp) THEN
        n_curr = FLOOR(curr_sim_time/dt_nudging)+1
        n_next = n_curr+1
        int_weight = curr_sim_time/dt_nudging-n_curr+1

        !ibd, exception for last time step
        IF (int_weight.eq.0._wp) THEN
          n_next = n_curr
        END IF

        ! stop if end of data reached
        IF (n_next .GT. nt) THEN
          CALL finish ( routine, 'end of SCM input file for nudging' )
        END IF
        
        ! Interpolation error catch
        IF (int_weight.LT.0 .OR. int_weight.GT.1) THEN
          WRITE (message_text,'(a,2(i0,", "),g0)') &
               'INTERPOLATION ERROR in nudging:', &
               n_curr,n_next,int_weight
          CALL message(routine, message_text)
        END IF

        IF(is_nudging_uv) THEN
          u_nudge         = u_nudg(:,n_curr)*(1.-int_weight) &
                        & + u_nudg(:,n_next)*    int_weight
          v_nudge         = v_nudg(:,n_curr)*(1.-int_weight) &
                        & + v_nudg(:,n_next)*    int_weight
        ENDIF
        IF(is_nudging_tq) THEN
          temp_nudge      = theta_nudg(:,n_curr)*(1.-int_weight) &
                        & + theta_nudg(:,n_next)*    int_weight 
          q_nudge(:,1)    = qv_nudg   (:,n_curr)*(1.-int_weight) &
                        & + qv_nudg   (:,n_next)*    int_weight 
        ENDIF
      ! single time step available
      ELSE
        IF(is_nudging_uv) THEN
          u_nudge         = u_nudg(:,1)
          v_nudge         = v_nudg(:,1) 
        ENDIF
        IF(is_nudging_tq) THEN
          temp_nudge      = theta_nudg(:,1)
          q_nudge(:,1)    = qv_nudg   (:,1)
        ENDIF
      END IF

      q_nudge(:,2) = 0.0_wp         ! qc undefined
      q_nudge(:,3) = 0.0_wp         ! qi undefined

    END IF

    !Convert theta tendency to temp at once
    ddt_temp_ls = exner_gb * ddt_temp_ls
    temp_nudge  = exner_gb * temp_nudge

    ! apply_ls_forcing has been called, tell SCM setup
    lscm_ls_forcing_ini = .TRUE.

  END SUBROUTINE apply_ls_forcing

!-------------------------------------------------------------------------------

END MODULE mo_ls_forcing



