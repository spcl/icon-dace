!
! @brief Subroutine read_bc_sulfate_aerosol reads monthly mean sulfate
!
! Concentrations provided which were calculated as a climatology by
! ECHAM-HAM. Note: this data are created by the original echam5-ham
! (version 5.4.01 with stratospheric extensions for HAM) and not by
! any ICON version with echam physics (now called aes physics).
!
! The sulfat concentrations can be used as background conditions to
! provide particle number and mass of soluble sulfate aerosols.
!
! The basic idea is based on Code by Sebastion Rast, Hui Wan, Florian
! Prill and Guenther Zaengl.
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

MODULE mo_bc_sulfate_aerosol

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish, message, message_text
  USE mo_model_domain,       ONLY: t_patch, p_patch
  USE mo_master_config,      ONLY: isRestart
  USE mo_time_config,        ONLY: time_config
  USE mo_io_config,          ONLY: default_read_method
  USE mo_grid_config,        ONLY: n_dom
  USE mo_parallel_config,    ONLY: nproma
  USE mo_read_interface,     ONLY: t_stream_id, openInputFile, closeFile, &
       &                           read_1D, read_2D_time, read_3D_time, &
       &                           on_cells
  USE mo_physical_constants, ONLY: grav, amso4, ams
  USE mo_ifs_coord,          ONLY: t_vct, geopot
  USE mo_util_phys,          ONLY: virtual_temp
  USE mo_nonhydro_types,     ONLY: t_nh_metrics
  USE mo_nonhydro_state,     ONLY: p_nh_state
  USE mo_nh_vert_interp,     ONLY: prepare_lin_intp, lin_intp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_echam_bc_sulfate_aerosol
  PUBLIC :: destruct_echam_bc_sulfate_aerosol
  PUBLIC :: ext_sulfate_aerosol

  REAL(wp), PARAMETER     :: rescale = amso4/ams*1.0e6_wp

  TYPE t_ext_sulfate_aerosol
    REAL(wp), ALLOCATABLE :: so4_ns(:,:,:)
    REAL(wp), ALLOCATABLE :: so4_ks(:,:,:)
    REAL(wp), ALLOCATABLE :: so4_as(:,:,:)
    REAL(wp), ALLOCATABLE :: so4_cs(:,:,:)
    REAL(wp), ALLOCATABLE :: num_ns(:,:,:)
    REAL(wp), ALLOCATABLE :: num_ks(:,:,:)
    REAL(wp), ALLOCATABLE :: num_as(:,:,:)
    REAL(wp), ALLOCATABLE :: num_cs(:,:,:)
  END TYPE t_ext_sulfate_aerosol

  TYPE(t_ext_sulfate_aerosol), ALLOCATABLE, TARGET :: ext_sulfate_aerosol(:)

CONTAINS

  SUBROUTINE destruct_echam_bc_sulfate_aerosol

    CHARACTER(len=*), PARAMETER       :: routine = 'mo_bc_sulfate:read_echam_bc_sulfate_aerosol'

    INTEGER :: jg

    ! needs explicit deallocate
    DO jg = 1, n_dom
      DEALLOCATE(ext_sulfate_aerosol(jg)%so4_ns)
      DEALLOCATE(ext_sulfate_aerosol(jg)%so4_ks)
      DEALLOCATE(ext_sulfate_aerosol(jg)%so4_as)
      DEALLOCATE(ext_sulfate_aerosol(jg)%so4_cs)

      DEALLOCATE(ext_sulfate_aerosol(jg)%num_ns)
      DEALLOCATE(ext_sulfate_aerosol(jg)%num_ks)
      DEALLOCATE(ext_sulfate_aerosol(jg)%num_as)
      DEALLOCATE(ext_sulfate_aerosol(jg)%num_cs)
    END DO

  END SUBROUTINE destruct_echam_bc_sulfate_aerosol

  SUBROUTINE read_echam_bc_sulfate_aerosol

    CHARACTER(len=*), PARAMETER       :: routine = 'mo_bc_sulfate:read_echam_bc_sulfate_aerosol'

    CHARACTER(len=512)                :: fname
    TYPE(t_stream_id)                 :: stream_id
    CHARACTER(len=2)                  :: cjg
    INTEGER                           :: month
    TYPE(t_patch), POINTER            :: lp_patch
    INTEGER                           :: mlev
    INTEGER                           :: nlev
    INTEGER                           :: nhyi
    INTEGER                           :: npromz
    INTEGER                           :: nblks
    INTEGER                           :: ncells
    INTEGER                           :: nlen
    INTEGER                           :: jg, jb

    ! variables in echam/ham

    REAL(wp), POINTER                 :: hyai(:), hybi(:)

    REAL(wp), POINTER                 :: aps(:,:,:)                     ! surface pressure
    REAL(wp), POINTER                 :: geosp(:,:,:)                   ! surface geopotential
    REAL(wp), POINTER                 :: st(:,:,:,:)                    ! temperature
    REAL(wp), POINTER                 :: q(:,:,:,:)                     ! water vapour
    REAL(wp), POINTER                 :: xl(:,:,:,:)                    ! cloud liquid water
    REAL(wp), POINTER                 :: xi(:,:,:,:)                    ! cloud ice

    REAL(wp), POINTER                 :: tv(:,:,:)

    REAL(wp), POINTER                 :: so4_ns(:,:,:,:)
    REAL(wp), POINTER                 :: so4_ks(:,:,:,:)
    REAL(wp), POINTER                 :: so4_as(:,:,:,:)
    REAL(wp), POINTER                 :: so4_cs(:,:,:,:)
    REAL(wp), POINTER                 :: num_ns(:,:,:,:)
    REAL(wp), POINTER                 :: num_ks(:,:,:,:)
    REAL(wp), POINTER                 :: num_as(:,:,:,:)
    REAL(wp), POINTER                 :: num_cs(:,:,:,:)

    TYPE(t_vct)                       :: echams_vct

    REAL(wp), POINTER                 :: p(:,:), p_ifc(:,:)
    REAL(wp), POINTER                 :: geop(:,:), geop_ifc(:,:)

    REAL(wp), POINTER                 :: z_echam(:,:,:)

    REAL(wp), ALLOCATABLE             :: rdalpha(:,:) ! rd*alpha at pressure and sigma levels
    REAL(wp), ALLOCATABLE             :: rdlnpr(:,:)  ! rd*ln(p(k+.5)/p(k-.5))
    REAL(wp), ALLOCATABLE             :: delp(:,:)    ! p(k+.5)-p(k-.5) of the reference surface pressure
    REAL(wp), ALLOCATABLE             :: rdelp(:,:)   ! reciprocal of delpr
    REAL(wp), ALLOCATABLE             :: lnp_ifc(:,:)

    REAL(wp), ALLOCATABLE             :: wfac_lin(:,:,:)
    REAL(wp), ALLOCATABLE             :: wfacpbl1(:,:), wfacpbl2(:,:)
    INTEGER, ALLOCATABLE              :: idx0_lin(:,:,:)
    INTEGER, ALLOCATABLE              :: bot_idx_lin(:,:)
    INTEGER, ALLOCATABLE              :: kpbl1(:,:), kpbl2(:,:)

    TYPE(t_nh_metrics), POINTER       :: p_metrics

    month = time_config%tc_current_date%date%month

    ! allocate once only structure for all grids
    IF (.NOT. ALLOCATED(ext_sulfate_aerosol)) ALLOCATE(ext_sulfate_aerosol(n_dom))
!$ACC ENTER DATA CREATE( ext_sulfate_aerosol )

    DO jg = 1, n_dom

      ! In case the model is freshly started read background conditions:
      IF (.NOT. isRestart()) THEN
        !
        lp_patch => p_patch(jg)
        nblks  = lp_patch%nblks_c
        npromz = lp_patch%npromz_c
        !
        ! Sulfate is constant in time. Data are monthy mean values for
        ! one year.  We read the month according to the start time.
        !
        IF (n_dom > 1) THEN
          WRITE(cjg,'(i2.2)') jg
          fname = 'bc_sulfate_DOM'//TRIM(cjg)//'.nc'
          CALL finish(routine, 'Multiple domains not supported yet.')
        ELSE
          fname = 'bc_sulfate.nc'
        END IF
        !
        WRITE(message_text,'(a, i2.2)') &
             & 'Read constant-in-time monthly climatology sulfate from file: '//TRIM(fname)//" for month ", month
        CALL message(routine, message_text)

        CALL openInputFile(stream_id, fname, lp_patch, default_read_method)

        CALL read_1d(file_id=stream_id%file_id, variable_name='hyai',  return_pointer=hyai)
        CALL read_1d(file_id=stream_id%file_id, variable_name='hybi',  return_pointer=hybi)

        CALL read_2d_time(stream_id=stream_id, location=on_cells, variable_name='aps',   return_pointer=aps,   &
             &            start_timestep=month,end_timestep=month)
        CALL read_2d_time(stream_id=stream_id, location=on_cells, variable_name='geosp', return_pointer=geosp, &
             &            start_timestep=month,end_timestep=month)

        CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name='st', return_pointer=st, &
             &            start_timestep=month,end_timestep=month)
        CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name='q',  return_pointer=q,  &
             &            start_timestep=month,end_timestep=month)
        CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name='xl', return_pointer=xl, &
             &            start_timestep=month,end_timestep=month)
        CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name='xi', return_pointer=xi, &
             &            start_timestep=month,end_timestep=month)

        CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name='SO4_NS', return_pointer=so4_ns, &
             &            start_timestep=month,end_timestep=month)
        CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name='SO4_KS', return_pointer=so4_ks, &
             &            start_timestep=month,end_timestep=month)
        CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name='SO4_AS', return_pointer=so4_as, &
             &            start_timestep=month,end_timestep=month)
        CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name='SO4_CS', return_pointer=so4_cs, &
             &            start_timestep=month,end_timestep=month)
        CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name='NUM_NS', return_pointer=num_ns, &
             &            start_timestep=month,end_timestep=month)
        CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name='NUM_KS', return_pointer=num_ks, &
             &            start_timestep=month,end_timestep=month)
        CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name='NUM_AS', return_pointer=num_as, &
             &            start_timestep=month,end_timestep=month)
        CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name='NUM_CS', return_pointer=num_cs, &
             &            start_timestep=month,end_timestep=month)

        CALL closeFile(stream_id)
        !
        ! create vertical grid description of echam

        nhyi = SIZE(hyai)

        ncells = SIZE(st,1)
        mlev   = SIZE(st,2)
        nblks  = SIZE(st,3)

        ALLOCATE(tv(ncells,mlev,nblks))

        ALLOCATE(delp(nproma,mlev), rdelp(nproma,mlev), rdlnpr(nproma,mlev), rdalpha(nproma,mlev))
        ALLOCATE(geop(nproma,mlev), p(nproma,mlev))
        ALLOCATE(p_ifc(nproma,nhyi), lnp_ifc(nproma,nhyi), geop_ifc(nproma,nhyi))
        ALLOCATE(z_echam(nproma,mlev,nblks))

        CALL echams_vct%construct(nhyi-1, hyai, hybi)   ! number of half levels required, but full levels given in nhyi
        CALL virtual_temp(lp_patch, st(:,:,:,1), q(:,:,:,1), xl(:,:,:,1), xi(:,:,:,1), temp_v=tv)
        DO jb = 1, nblks
          IF (jb /= nblks) THEN
            nlen = nproma
          ELSE
            nlen = npromz
          ENDIF
          CALL echams_vct%half_level_pressure(aps(:,jb,1), nproma, nlen, mlev, p_ifc)
          CALL echams_vct%full_level_pressure(p_ifc, nproma, nlen, mlev, p)
          CALL echams_vct%auxhyb(p_ifc, nproma, nlen, mlev, delp, rdelp, lnp_ifc, rdlnpr, rdalpha)
          CALL geopot(tv(:,:,jb), rdlnpr, rdalpha, geosp(:,jb,1), nproma, 1, nlen, mlev, geop, geop_ifc)
          ! Compute 3D height coordinate field
          z_echam(1:nlen,1:mlev,jb) = geop(1:nlen,1:mlev)/grav
        ENDDO

        so4_ns(:,:,:,:) = rescale*so4_ns(:,:,:,:)
        so4_ks(:,:,:,:) = rescale*so4_ks(:,:,:,:)
        so4_as(:,:,:,:) = rescale*so4_as(:,:,:,:)
        so4_cs(:,:,:,:) = rescale*so4_cs(:,:,:,:)

        num_ns(:,:,:,:) = rescale*num_ns(:,:,:,:)
        num_ks(:,:,:,:) = rescale*num_ks(:,:,:,:)
        num_as(:,:,:,:) = rescale*num_as(:,:,:,:)
        num_cs(:,:,:,:) = rescale*num_cs(:,:,:,:)

        nlev = lp_patch%nlev

        ALLOCATE(wfac_lin(nproma,nlev,lp_patch%nblks_c))
        ALLOCATE(wfacpbl1(nproma,lp_patch%nblks_c), wfacpbl2(nproma,lp_patch%nblks_c))
        ALLOCATE(idx0_lin(nproma,nlev,lp_patch%nblks_c))
        ALLOCATE(bot_idx_lin(nproma,lp_patch%nblks_c))
        ALLOCATE(kpbl1(nproma,lp_patch%nblks_c), kpbl2(nproma,lp_patch%nblks_c))

        ALLOCATE(ext_sulfate_aerosol(jg)%so4_ns(nproma,nlev,lp_patch%nblks_c))
        ALLOCATE(ext_sulfate_aerosol(jg)%so4_ks(nproma,nlev,lp_patch%nblks_c))
        ALLOCATE(ext_sulfate_aerosol(jg)%so4_as(nproma,nlev,lp_patch%nblks_c))
        ALLOCATE(ext_sulfate_aerosol(jg)%so4_cs(nproma,nlev,lp_patch%nblks_c))

        ALLOCATE(ext_sulfate_aerosol(jg)%num_ns(nproma,nlev,lp_patch%nblks_c))
        ALLOCATE(ext_sulfate_aerosol(jg)%num_ks(nproma,nlev,lp_patch%nblks_c))
        ALLOCATE(ext_sulfate_aerosol(jg)%num_as(nproma,nlev,lp_patch%nblks_c))
        ALLOCATE(ext_sulfate_aerosol(jg)%num_cs(nproma,nlev,lp_patch%nblks_c))

        p_metrics => p_nh_state(jg)%metrics

        CALL prepare_lin_intp(z_echam, p_metrics%z_mc,                          &
             &                lp_patch%nblks_c, lp_patch%npromz_c, mlev, nlev,  &
             &                wfac_lin, idx0_lin, bot_idx_lin)

        CALL lin_intp(so4_ns(:,:,:,1), ext_sulfate_aerosol(jg)%so4_ns,          &
             &        lp_patch%nblks_c, lp_patch%npromz_c, mlev, nlev,          &
             &        wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,         &
             &        wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE.,    &
             &        l_pd_limit=.TRUE.)
        CALL lin_intp(so4_ks(:,:,:,1),ext_sulfate_aerosol(jg)%so4_ks,           &
             &        lp_patch%nblks_c, lp_patch%npromz_c, mlev, nlev,          &
             &        wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,         &
             &        wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE.,    &
             &        l_pd_limit=.TRUE.)
        CALL lin_intp(so4_as(:,:,:,1),ext_sulfate_aerosol(jg)%so4_as,           &
             &        lp_patch%nblks_c, lp_patch%npromz_c, mlev, nlev,          &
             &        wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,         &
             &        wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE.,    &
             &        l_pd_limit=.TRUE.)
        CALL lin_intp(so4_cs(:,:,:,1),ext_sulfate_aerosol(jg)%so4_cs,           &
             &        lp_patch%nblks_c, lp_patch%npromz_c, mlev, nlev,          &
             &        wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,         &
             &        wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE.,    &
             &        l_pd_limit=.TRUE.)

        CALL lin_intp(num_ns(:,:,:,1),ext_sulfate_aerosol(jg)%num_ns,           &
             &        lp_patch%nblks_c, lp_patch%npromz_c, mlev, nlev,          &
             &        wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,         &
             &        wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE.,    &
             &        l_pd_limit=.TRUE.)
        CALL lin_intp(num_ks(:,:,:,1),ext_sulfate_aerosol(jg)%num_ks,           &
             &        lp_patch%nblks_c, lp_patch%npromz_c, mlev, nlev,          &
             &        wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,         &
             &        wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE.,    &
             &        l_pd_limit=.TRUE.)
        CALL lin_intp(num_as(:,:,:,1),ext_sulfate_aerosol(jg)%num_as,           &
             &        lp_patch%nblks_c, lp_patch%npromz_c, mlev, nlev,          &
             &        wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,         &
             &        wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE.,    &
             &        l_pd_limit=.TRUE.)
        CALL lin_intp(num_cs(:,:,:,1),ext_sulfate_aerosol(jg)%num_cs,           &
             &        lp_patch%nblks_c, lp_patch%npromz_c, mlev, nlev,          &
             &        wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,         &
             &        wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE.,    &
             &        l_pd_limit=.TRUE.)

      END IF

    ENDDO

  END SUBROUTINE read_echam_bc_sulfate_aerosol

END MODULE mo_bc_sulfate_aerosol
