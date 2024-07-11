! Subroutine to initialize the CBL case for HDCP2
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

MODULE mo_nh_torus_exp

!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2008
!
!-------------------------------------------------------------------------
!
!
!

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, SUCCESS
  USE mo_io_units,            ONLY: find_next_free_unit
  USE mo_io_config,           ONLY: default_read_method
  USE mo_read_interface,      ONLY: openInputFile, closeFile, on_cells, &
     &                              t_stream_id, read_2D_extdim, read_3D_extdim
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_physical_constants,  ONLY: rd, cpd, p0ref, cvd_o_rd, rd_o_cpd, &
     &                              grav, alv, vtmpc1, lh_v=>alv
  USE mo_nh_testcases_nml,    ONLY: u_cbl, v_cbl, th_cbl, psfc_cbl, &
                                    bubctr_x, bubctr_y, nh_test_name, is_dry_cbl
  USE mo_nh_wk_exp,           ONLY: bub_amp, bub_ver_width, bub_hor_width, bubctr_z
  USE mo_model_domain,        ONLY: t_patch
  USE mo_math_constants,      ONLY: rad2deg, pi_2
  USE mo_loopindices,         ONLY: get_indices_e
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics, t_nh_ref
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_tend
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_parallel_config,     ONLY: nproma
  USE mo_math_utilities,      ONLY: plane_torus_distance
  USE mo_sync,                ONLY: sync_patch_array, SYNC_C
  USE mo_nh_init_utils,       ONLY: init_w
  USE mo_run_config,          ONLY: iqv, iqc
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_e
  USE mo_satad,               ONLY: spec_humi, sat_pres_water
  USE mo_les_config,          ONLY: les_config
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_nh_vert_interp_les,  ONLY: vert_intp_linear_1d
  USE mo_grid_config,         ONLY: l_scm_mode
  USE mo_scm_nml,             ONLY: i_scm_netcdf, scm_sfc_temp, scm_sfc_qv, scm_sfc_mom, lscm_icon_ini, &
     &                              lscm_ls_forcing_ini, lat_scm, lon_scm
  USE turb_data,              ONLY: vel_min  
  USE mo_lnd_nwp_config,      ONLY: nlev_soil
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_read_interface,      ONLY: nf
  USE mo_netcdf
  USE mo_mpi,                 ONLY: get_my_global_mpi_id

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_nh_state_cbl, cbl_stevens_fluxes, init_nh_state_rico,  &
            sfcflx_uniform, init_warm_bubble,                           &
            init_torus_ascii_sounding, init_torus_netcdf_sounding,      & 
            read_soil_profile_nc, read_ext_scm_nc, set_scm_bnd,         &
            read_soil_profile_nc_uf, read_ext_scm_nc_uf,                &
            init_torus_rcemip_analytical_sounding

!--------------------------------------------------------------------

   CONTAINS
!-------------------------------------------------------------------------
!
!
  !>
  !! Initialization of prognostic state vector for the nh CBL test case 
  !!  without moisture
  !!
  SUBROUTINE init_nh_state_cbl( ptr_patch, ptr_nh_prog,  ptr_nh_ref, ptr_nh_diag,  &
    &                           ptr_int, ptr_metrics)

    TYPE(t_patch),TARGET,  INTENT(INOUT)::  ptr_patch
    TYPE(t_int_state),     INTENT(IN)   ::  ptr_int
    TYPE(t_nh_prog),       INTENT(INOUT)::  ptr_nh_prog
    TYPE(t_nh_diag),       INTENT(INOUT)::  ptr_nh_diag
    TYPE(t_nh_metrics),    INTENT(IN)   ::  ptr_metrics      
    TYPE(t_nh_ref),        INTENT(INOUT)::  ptr_nh_ref

    REAL(wp) :: z_exner_h(1:nproma,ptr_patch%nlev+1), z_help(1:nproma) 
    REAL(wp) :: zvn1, zvn2, zu, zv, zt00, zh00, ex_sfc
    INTEGER  :: jc,jk,jb,i_startblk,i_startidx,i_endidx   !< loop indices
    INTEGER  :: nblks_c,npromz_c,nblks_e
    INTEGER  :: nlev, nlevp1                  !< number of full and half levels
    INTEGER  :: nlen, jcn, jbn, jg, ntropo

    REAL(wp), PARAMETER :: zh0     = 0._wp      !< height (m) above which temperature increases
    REAL(wp), PARAMETER :: lambda  = 1500._wp   !moist height from Stevens(2007)
    REAL(wp), PARAMETER :: dtdz_st = 0.04_wp    !< theta lapse rate in stratosphere (T>0!)
    REAL(wp), PARAMETER :: z_tropo = 11000._wp  !height tropopause
    REAL(wp), PARAMETER :: rh_sfc  = 0.8_wp    !RH at surface [1], default 0.8

  !-------------------------------------------------------------------------

    ! values for the blocking
    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c
    nblks_e  = ptr_patch%nblks_e

    ! number of vertical levels
    nlev   = ptr_patch%nlev
    nlevp1 = ptr_patch%nlevp1

    !patch id
    jg = ptr_patch%id

    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = psfc_cbl
    ex_sfc   = (psfc_cbl/p0ref)**rd_o_cpd
    les_config(jg)%psfc = psfc_cbl

    ! Tracers: all zero by default
    ptr_nh_prog%tracer(:,:,:,:) = 0._wp
    ptr_nh_prog%tke(:,:,:)      = 0._wp

    DO jb = 1, nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      !Tracers
      IF(.NOT.les_config(jg)%is_dry_cbl .AND. .NOT.is_dry_cbl)THEN
        DO jk = 1, nlev
          ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) = rh_sfc * spec_humi(sat_pres_water(th_cbl(1)),psfc_cbl) * &
                    EXP(-ptr_metrics%z_mc(1:nlen,jk,jb)/lambda)

        END DO
      END IF

      ntropo = 0
      DO jk = nlev, 1, -1
         ! init potential temperature
         z_help(1:nlen) = th_cbl(1) + max(0._wp, (ptr_metrics%z_mc(1:nlen,jk,jb)-zh0)*th_cbl(2))

         ! constant temperature above tropopause
         if ((ptr_metrics%z_mc(1,jk,jb) > z_tropo) .and. (ntropo == 0)) then
            ntropo = 1
            zt00   = z_help(1)
            zh00   = ptr_metrics%z_mc(1,jk,jb)
         endif
         if (ptr_metrics%z_mc(1,jk,jb) > z_tropo) then
            z_help(1:nlen) = zt00 + (ptr_metrics%z_mc(1:nlen,jk,jb)-zh00) * dtdz_st
         endif

         ! virtual potential temperature
         ptr_nh_prog%theta_v(1:nlen,jk,jb) = z_help(1:nlen) * ( 1._wp + &
           0.61_wp*ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) - ptr_nh_prog%tracer(1:nlen,jk,jb,iqc) ) 
      END DO

      !Get hydrostatic exner at the surface using surface pressure 
      z_exner_h(1:nlen,nlevp1) = ex_sfc
 
      !Get exner at full levels starting from exner at surface
      DO jk = nlev, 1, -1
         !exner at next half level after surface
         z_exner_h(1:nlen,jk) = z_exner_h(1:nlen,jk+1) - grav/cpd *     &
                                ptr_metrics%ddqz_z_full(1:nlen,jk,jb)/ &
                                ptr_nh_prog%theta_v(1:nlen,jk,jb)
        
         !exner at main levels
         ptr_nh_prog%exner(1:nlen,jk,jb) = 0.5_wp * &
                                     (z_exner_h(1:nlen,jk)+z_exner_h(1:nlen,jk+1))
      END DO

      DO jk = 1 , nlev
         ptr_nh_prog%rho(1:nlen,jk,jb) = (ptr_nh_prog%exner(1:nlen,jk,jb)**cvd_o_rd)*p0ref/rd / &
                                         ptr_nh_prog%theta_v(1:nlen,jk,jb)     
         ptr_nh_diag%pres(1:nlen,jk,jb) = ptr_nh_prog%rho(1:nlen,jk,jb)*rd*th_cbl(1)
      END DO !jk

    ENDDO !jb

!--------------------------------------------------------------------------------
    !Mean wind 
!--------------------------------------------------------------------------------
    i_startblk = ptr_patch%edges%start_blk(2,1)
    DO jb = i_startblk , nblks_e
     CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, i_startidx, i_endidx, 2)
     DO jk = 1 , nlev 
      DO jc = i_startidx, i_endidx

        !Torus geometry is flat so zu is only function of height which is same for all cells
        !But it is kept varyign with jc,jb to introduce topography lateron
        jcn  =   ptr_patch%edges%cell_idx(jc,jb,1)
        jbn  =   ptr_patch%edges%cell_blk(jc,jb,1)
        zu   =   u_cbl(1) + u_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)
        zv   =   v_cbl(1) + v_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)

        zvn1 =  zu * ptr_patch%edges%primal_normal_cell(jc,jb,1)%v1 + &
                zv * ptr_patch%edges%primal_normal_cell(jc,jb,1)%v2      
 
        jcn  =   ptr_patch%edges%cell_idx(jc,jb,2)
        jbn  =   ptr_patch%edges%cell_blk(jc,jb,2)
        zu   =   u_cbl(1) + u_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)
        zv   =   v_cbl(1) + v_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)
      
        zvn2 =  zu * ptr_patch%edges%primal_normal_cell(jc,jb,2)%v1 + &
                zv * ptr_patch%edges%primal_normal_cell(jc,jb,2)%v2      

        ptr_nh_prog%vn(jc,jk,jb) = ptr_int%c_lin_e(jc,1,jb)*zvn1 + &
                                   ptr_int%c_lin_e(jc,2,jb)*zvn2

        ptr_nh_ref%vn_ref(jc,jk,jb) = ptr_nh_prog%vn(jc,jk,jb)
      END DO
     END DO
    END DO     
    
    !W wind and reference
    CALL init_w(ptr_patch, ptr_int, ptr_nh_prog%vn, ptr_metrics%z_ifc, ptr_nh_prog%w)
    CALL sync_patch_array(SYNC_C, ptr_patch, ptr_nh_prog%w)
    ptr_nh_ref%w_ref = ptr_nh_prog%w


  END SUBROUTINE init_nh_state_cbl


  !>
  !! Initialization of prognostic state vector for the nh GATE test case 
  !! with moisture
  !!
  SUBROUTINE init_torus_ascii_sounding ( ptr_patch, ptr_nh_prog,  ptr_nh_ref,  &
                                         ptr_nh_diag, ptr_int, ptr_metrics)

    TYPE(t_patch),TARGET,  INTENT(INOUT)::  ptr_patch
    TYPE(t_int_state),     INTENT(IN)   ::  ptr_int
    TYPE(t_nh_prog),       INTENT(INOUT)::  ptr_nh_prog
    TYPE(t_nh_diag),       INTENT(INOUT)::  ptr_nh_diag
    TYPE(t_nh_metrics),    INTENT(IN)   ::  ptr_metrics      
    TYPE(t_nh_ref),        INTENT(INOUT)::  ptr_nh_ref

    INTEGER  :: je,jk,jb,i_startidx,i_endidx   !< loop indices
    INTEGER  :: nblks_c,npromz_c,nblks_e, jg
    INTEGER  :: nlev, nlevp1                   !< number of full and half levels
    INTEGER  :: nlen

    REAL(wp), DIMENSION(ptr_patch%nlev)   :: theta_in, &
                                          &  qv_in, u_in, v_in
    REAL(wp), DIMENSION(ptr_patch%nlev+1) :: w_in
    REAL(wp) :: zvn1, zvn2, zu, zv, psfc_in, ex_sfc
    REAL(wp) :: z_exner_h(1:nproma, ptr_patch%nlev+1), z_help(1:nproma) 

    CHARACTER(len=*), PARAMETER :: &
       &  routine = 'mo_nh_torus_exp:init_torus_ascii_sounding'
    !-------------------------------------------------------------------------
    
    ! Read the sounding file

    CALL read_ext_profile (ptr_metrics%z_mc(2,:,2), theta_in, qv_in, u_in, v_in, psfc_in)

    ! values for the blocking
    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c
    nblks_e  = ptr_patch%nblks_e

    ! number of vertical levels
    nlev     = ptr_patch%nlev
    nlevp1   = ptr_patch%nlevp1

    !patch id
    jg       = ptr_patch%id

    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = psfc_in
    ex_sfc   = (psfc_in/p0ref)**rd_o_cpd
    IF ( les_config(jg)%psfc /= psfc_in ) THEN
      CALL finish(routine,'Value of psfc in les_nml is inconsistent with data in sounding file!')
    END IF

    ! Tracers: all zero by default
    ptr_nh_prog%tracer(:,:,:,:) = 0._wp

    DO jb = 1, nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      DO jk = 1, nlev
        ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) = qv_in(jk)
      END DO

      DO jk = nlev, 1, -1
         z_help(1:nlen) = theta_in(jk)
         ptr_nh_prog%theta_v(1:nlen,jk,jb) = z_help(1:nlen) * ( 1._wp + &
           0.61_wp*ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) - ptr_nh_prog%tracer(1:nlen,jk,jb,iqc) ) 
      END DO

      !Get hydrostatic exner at the surface using surface pressure 
      z_exner_h(1:nlen,nlevp1) = ex_sfc
 
      !Get exner at full levels starting from exner at surface
      DO jk = nlev, 1, -1
         !exner at next half level after surface
         z_exner_h(1:nlen,jk) = z_exner_h(1:nlen,jk+1) - grav/cpd *     &
                                ptr_metrics%ddqz_z_full(1:nlen,jk,jb)/ &
                                ptr_nh_prog%theta_v(1:nlen,jk,jb)
        
         !exner at main levels
         ptr_nh_prog%exner(1:nlen,jk,jb) = 0.5_wp * &
                                     (z_exner_h(1:nlen,jk)+z_exner_h(1:nlen,jk+1))
      END DO

      DO jk = 1 , nlev
        ptr_nh_prog%rho(1:nlen,jk,jb) = (ptr_nh_prog%exner(1:nlen,jk,jb)**cvd_o_rd)*p0ref/rd / &
                                         ptr_nh_prog%theta_v(1:nlen,jk,jb)     
      END DO

    ENDDO

    !Mean wind 
    DO jb = 1 , nblks_e
      CALL get_indices_e( ptr_patch, jb, 1, nblks_e, i_startidx, i_endidx, grf_bdywidth_e+1)
      DO jk = 1 , nlev 
        DO je = i_startidx, i_endidx

          zu   =   u_in(jk)
          zv   =   v_in(jk)
  
          zvn1 =  zu * ptr_patch%edges%primal_normal_cell(je,jb,1)%v1 + &
                  zv * ptr_patch%edges%primal_normal_cell(je,jb,1)%v2      
   
          zvn2 =  zu * ptr_patch%edges%primal_normal_cell(je,jb,2)%v1 + &
                  zv * ptr_patch%edges%primal_normal_cell(je,jb,2)%v2      
  
          ptr_nh_prog%vn(je,jk,jb) = ptr_int%c_lin_e(je,1,jb)*zvn1 + &
                                     ptr_int%c_lin_e(je,2,jb)*zvn2
  
          ptr_nh_ref%vn_ref(je,jk,jb) = ptr_nh_prog%vn(je,jk,jb)
        END DO
      END DO
    END DO     
    
    !W wind and reference
    CALL init_w(ptr_patch, ptr_int, ptr_nh_prog%vn, ptr_metrics%z_ifc, ptr_nh_prog%w)
    CALL sync_patch_array(SYNC_C, ptr_patch, ptr_nh_prog%w)
    ptr_nh_ref%w_ref  = ptr_nh_prog%w

  END SUBROUTINE init_torus_ascii_sounding


  !>
  !! Initialization of prognostic state vector for with SCM netcdf intput file init_SCM.nc.
  !!
  SUBROUTINE init_torus_netcdf_sounding ( ptr_patch, ptr_nh_prog,  ptr_nh_ref, ptr_nh_diag,  &
                                          ptr_int, ptr_metrics, ext_data)

    TYPE(t_patch),TARGET,  INTENT(INOUT)::  ptr_patch
    TYPE(t_int_state),     INTENT(IN)   ::  ptr_int
    TYPE(t_nh_prog),       INTENT(INOUT)::  ptr_nh_prog
    TYPE(t_nh_diag),       INTENT(INOUT)::  ptr_nh_diag
    TYPE(t_nh_metrics),    INTENT(IN)   ::  ptr_metrics      
    TYPE(t_nh_ref),        INTENT(INOUT)::  ptr_nh_ref
    TYPE(t_external_data), INTENT(INOUT)::  ext_data

    INTEGER  :: je,jk,jb,i_startidx,i_endidx   !< loop indices
    INTEGER  :: nblks_c,npromz_c,nblks_e, jg
    INTEGER  :: nlev, nlevp1                   !< number of full and half levels
    INTEGER  :: nlen

    REAL(wp), DIMENSION(ptr_patch%nlev)   :: theta_in, thetav_in, exner_in, rho_in, &
                                          &  qv_in, qc_in, qi_in, u_in, v_in, o3_in
    REAL(wp), DIMENSION(ptr_patch%nlev+1) :: tke_in, w_in
    REAL(wp) :: zvn1, zvn2, zu, zv, psfc_in, ex_sfc
    REAL(wp) :: z_exner_h(1:nproma, ptr_patch%nlev+1), z_help(1:nproma) 

    CHARACTER(len=*), PARAMETER :: &
       &  routine = 'mo_nh_torus_exp:init_torus_netcdf_sounding'
    !-------------------------------------------------------------------------
    
    ! Read the sounding file

    ! SCM: always using NETCDF input (i_scm_netcdf>0)

    SELECT CASE (i_scm_netcdf)
    CASE (2)
      !read unified netcdf file
      CALL read_ext_profile_nc_uf (ptr_metrics%z_mc(1,:,1), ptr_metrics%z_ifc(1,:,1),    &
         & theta_in, thetav_in, exner_in, rho_in, qv_in, qc_in, qi_in, u_in, v_in, w_in, &
         & tke_in, psfc_in, o3_in)
    CASE (1)
      !read ICON netcdf file
      CALL read_ext_profile_nc (ptr_metrics%z_mc(1,:,1), ptr_metrics%z_ifc(1,:,1),       &
         & theta_in, thetav_in, exner_in, rho_in, qv_in, qc_in, qi_in, u_in, v_in, w_in, &
         & tke_in, psfc_in, o3_in)
    END SELECT

    ! values for the blocking
    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c
    nblks_e  = ptr_patch%nblks_e

    ! number of vertical levels
    nlev     = ptr_patch%nlev
    nlevp1   = ptr_patch%nlevp1

    !patch id
    jg       = ptr_patch%id

    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = psfc_in
    ex_sfc   = (psfc_in/p0ref)**rd_o_cpd
    IF ( les_config(jg)%psfc /= psfc_in .AND. (.NOT. l_scm_mode) ) THEN
      WRITE(*,*) 'les_config(jg)%psfc, psfc_in)', les_config(jg)%psfc, psfc_in
      CALL finish(routine,'Value of psfc in les_nml is inconsistent with data in sounding file!')
    END IF

    ! Tracers: all zero by default
    ptr_nh_prog%tracer(:,:,:,:) = 0._wp
    ptr_nh_prog%tke(:,:,:)      = 0._wp

    DO jb = 1, nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      DO jk = 1, nlev
        ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) = qv_in(jk)
        ptr_nh_prog%tracer(1:nlen,jk,jb,iqc) = qc_in(jk)
        ext_data%atm%o3   (1:nlen,jk,jb)     = o3_in(jk)
      END DO

      DO jk = 1 , nlev+1
        ptr_nh_prog%tke(1:nlen,jk,jb)        = tke_in(jk)
      END DO

      ! initial conditions from ICON forecast output
      IF ( lscm_icon_ini ) THEN

        DO jk = 1, nlev
          ptr_nh_prog%theta_v(1:nlen,jk,jb)  = thetav_in(jk)
          ptr_nh_prog%exner  (1:nlen,jk,jb)  = exner_in (jk)
          ptr_nh_prog%rho    (1:nlen,jk,jb)  = rho_in   (jk)
        END DO
        DO jk = 1 , nlev+1
          ptr_nh_prog%w      (1:nlen,jk,jb)  = w_in     (jk)
        END DO

      ! idealized initial condition (e.g. field experiment)
      ELSE

        DO jk = nlev, 1, -1
          z_help(1:nlen) = theta_in(jk)
          ptr_nh_prog%theta_v(1:nlen,jk,jb) = z_help(1:nlen) * ( 1._wp + &
            0.61_wp*ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) - ptr_nh_prog%tracer(1:nlen,jk,jb,iqc) ) 
        END DO
  
        !Get hydrostatic exner at the surface using surface pressure 
        z_exner_h(1:nlen,nlevp1) = ex_sfc
   
        !Get exner at full levels starting from exner at surface
        DO jk = nlev, 1, -1
          !exner at next half level after surface
          z_exner_h(1:nlen,jk) = z_exner_h(1:nlen,jk+1) - grav/cpd *    &
                                 ptr_metrics%ddqz_z_full(1:nlen,jk,jb)/ &
                                 ptr_nh_prog%theta_v(1:nlen,jk,jb)
          
          !exner at main levels
          ptr_nh_prog%exner(1:nlen,jk,jb) = 0.5_wp * (z_exner_h(1:nlen,jk)+z_exner_h(1:nlen,jk+1))
          IF ( ptr_nh_prog%exner(1,jk,jb)<0.0_wp ) THEN 
            CALL finish('testcases/mo_nh_torus_exp.f90', 'Exner function is negative')
          ENDIF
        END DO
  
        DO jk = 1 , nlev
          ptr_nh_prog%rho(1:nlen,jk,jb) = (ptr_nh_prog%exner(1:nlen,jk,jb)**cvd_o_rd)*p0ref/rd / &
                                           ptr_nh_prog%theta_v(1:nlen,jk,jb)     
        END DO

      END IF

    ENDDO !jb

    !Mean wind 
    DO jb = 1 , nblks_e
      CALL get_indices_e( ptr_patch, jb, 1, nblks_e, i_startidx, i_endidx, grf_bdywidth_e+1)
      DO jk = 1 , nlev 
        DO je = i_startidx, i_endidx
  
          zu   =   u_in(jk)
          zv   =   v_in(jk)
  
          zvn1 =  zu * ptr_patch%edges%primal_normal_cell(je,jb,1)%v1 + &
                  zv * ptr_patch%edges%primal_normal_cell(je,jb,1)%v2      
   
          zvn2 =  zu * ptr_patch%edges%primal_normal_cell(je,jb,2)%v1 + &
                  zv * ptr_patch%edges%primal_normal_cell(je,jb,2)%v2      
  
          ptr_nh_prog%vn(je,jk,jb) = ptr_int%c_lin_e(je,1,jb)*zvn1 + &
                                     ptr_int%c_lin_e(je,2,jb)*zvn2
    
          ptr_nh_ref%vn_ref(je,jk,jb) = ptr_nh_prog%vn(je,jk,jb)

        END DO
      END DO
    END DO
      
    !W wind and reference
    ! idealized initial condition (e.g. field experiment)
    IF ( l_scm_mode .AND. (.NOT. lscm_icon_ini) ) THEN
      CALL init_w(ptr_patch, ptr_int, ptr_nh_prog%vn, ptr_metrics%z_ifc, ptr_nh_prog%w)
    END IF

    CALL sync_patch_array(SYNC_C, ptr_patch, ptr_nh_prog%w)
    ptr_nh_ref%w_ref  = ptr_nh_prog%w


  END SUBROUTINE init_torus_netcdf_sounding


!-------------------------------------------------------------------------
  
  !>
  !! Initialization of prognostic state vector for the nh RICO test case 
  !!  with moisture
  !!
  SUBROUTINE init_nh_state_rico( ptr_patch, ptr_nh_prog,  ptr_nh_ref, ptr_nh_diag,  &
    &                           ptr_int, ptr_metrics)

    ! INPUT PARAMETERS:
    TYPE(t_patch),TARGET,  INTENT(INOUT):: &  !< patch on which computation is performed
      &  ptr_patch
    TYPE(t_int_state),     INTENT(IN)   :: &
      &  ptr_int
    TYPE(t_nh_prog),       INTENT(INOUT):: &  !< prognostic state vector
      &  ptr_nh_prog
    TYPE(t_nh_diag),       INTENT(INOUT):: &  !< diagnostic state vector
      &  ptr_nh_diag
    TYPE(t_nh_metrics),    INTENT(IN)   :: &
      &  ptr_metrics                          !< NH metrics state
    TYPE(t_nh_ref),        INTENT(INOUT):: &  !< reference state vector
      &  ptr_nh_ref

    REAL(wp) :: rho_sfc, z_help(1:nproma), zvn1, zvn2, zu, zv
    INTEGER  :: je,jk,jb,i_startidx,i_endidx   !< loop indices
    INTEGER  :: nblks_c,npromz_c,nblks_e
    INTEGER  :: nlev                      !< number of full levels
    INTEGER  :: nlen, jcn, jbn, jg

    !DEFINED PARAMETERS (RICO case):
    REAL(wp), PARAMETER :: zpsfc   = 101540._wp   !< surface pressure
    REAL(wp), PARAMETER :: zh1     = 740._wp      !< height (m) above which temperature increases
    REAL(wp), PARAMETER :: ztsfc   = 297.9_wp
    REAL(wp), PARAMETER :: zh2     = 3260._wp     !< moist height for RICO
    REAL(wp), PARAMETER :: zh3     = 15000._wp    !< height for extrapolated RICO profiles
    REAL(wp), PARAMETER :: zh4     = 30000._wp    !< height for extrapolated RICO profiles
    REAL(wp), PARAMETER :: zh5     = 60000._wp    !< height for extrapolated RICO profiles

  !-------------------------------------------------------------------------

    ! values for the blocking
    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c
    nblks_e  = ptr_patch%nblks_e

    ! number of vertical levels
    nlev   = ptr_patch%nlev

    !patch id
    jg = ptr_patch%id

    !Set some reference density    
    rho_sfc = zpsfc / (rd * les_config(jg)%sst)

    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = zpsfc

    ! Tracers: all zero by default
    ptr_nh_prog%tracer(:,:,:,:) = 0._wp

    DO jb = 1, nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      !Tracers
      DO jk = 1, nlev
        ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) = min((0.016_wp - 0.0022_wp * ptr_metrics%z_mc(1:nlen,jk,jb)/zh1), &
                                              (0.0138_wp - 0.0114_wp * (ptr_metrics%z_mc(1:nlen,jk,jb)-zh1)/(zh2 - zh1)))
        ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) = max(ptr_nh_prog%tracer(1:nlen,jk,jb,iqv),             &
                                              (0.0024_wp - 0.0006_wp * (ptr_metrics%z_mc(1:nlen,jk,jb) - zh2)/(4000._wp - zh2)))
        ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) = max(ptr_nh_prog%tracer(1:nlen,jk,jb,iqv), 3e-6_wp)                            
      END DO

      DO jk = 1, nlev
       ! init potential temperature
       z_help(1:nlen) = ztsfc + max(0._wp, (ptr_metrics%z_mc(1:nlen,jk,jb)-zh1)*19.1_wp/(4000._wp-zh1))
       z_help(1:nlen) = max(z_help(1:nlen), ztsfc +(ptr_metrics%z_mc(1:nlen,jk,jb)-zh3)*502._wp/(zh4-zh3))
       z_help(1:nlen) = max(z_help(1:nlen), ztsfc +(ptr_metrics%z_mc(1:nlen,jk,jb)-zh4)*2600._wp/(zh5-zh4))

       ! virtual potential temperature
       ptr_nh_prog%theta_v(1:nlen,jk,jb) = z_help(1:nlen) * ( 1._wp + &
           0.61_wp*ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) - ptr_nh_prog%tracer(1:nlen,jk,jb,iqc) ) 
      END DO

      !Get hydrostatic pressure and exner at lowest level
      ptr_nh_diag%pres(1:nlen,nlev,jb) = zpsfc - rho_sfc * ptr_metrics%geopot(1:nlen,nlev,jb)
      ptr_nh_prog%exner(1:nlen,nlev,jb) = (ptr_nh_diag%pres(1:nlen,nlev,jb)/p0ref)**rd_o_cpd 

      !Get exner at other levels
      DO jk = nlev-1, 1, -1
         z_help(1:nlen) = 0.5_wp * ( ptr_nh_prog%theta_v(1:nlen,jk,jb) +  &
                                     ptr_nh_prog%theta_v(1:nlen,jk+1,jb) )
   
         ptr_nh_prog%exner(1:nlen,jk,jb) = ptr_nh_prog%exner(1:nlen,jk+1,jb) &
            &  -grav/cpd*ptr_metrics%ddqz_z_half(1:nlen,jk+1,jb)/z_help(1:nlen)
      END DO

      DO jk = 1 , nlev
        ptr_nh_prog%rho(1:nlen,jk,jb) = (ptr_nh_prog%exner(1:nlen,jk,jb)**cvd_o_rd)*p0ref/rd / &
                                         ptr_nh_prog%theta_v(1:nlen,jk,jb)     
      END DO !jk

    ENDDO !jb

    !Mean wind 
    DO jb = 1 , nblks_e
     CALL get_indices_e( ptr_patch, jb, 1, nblks_e, i_startidx, i_endidx, grf_bdywidth_e+1)
     DO jk = 1 , nlev 
      DO je = i_startidx, i_endidx

        !Torus geometry is flat so zu is only function of height which is same for all cells
        !But it is kept varyign with jc,jb to introduce topography lateron
        jcn  =   ptr_patch%edges%cell_idx(je,jb,1)
        jbn  =   ptr_patch%edges%cell_blk(je,jb,1)
        zu   =   u_cbl(1) + u_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)
        zu   =   min(zu, -2._wp)
        zv   =   v_cbl(1) + v_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)

        zvn1 =  zu * ptr_patch%edges%primal_normal_cell(je,jb,1)%v1 + &
                zv * ptr_patch%edges%primal_normal_cell(je,jb,1)%v2      
 
        jcn  =   ptr_patch%edges%cell_idx(je,jb,2)
        jbn  =   ptr_patch%edges%cell_blk(je,jb,2)
        zu   =   u_cbl(1) + u_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)
        zu   =   min(zu, -2._wp)
        zv   =   v_cbl(1) + v_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)
      
        zvn2 =  zu * ptr_patch%edges%primal_normal_cell(je,jb,2)%v1 + &
                zv * ptr_patch%edges%primal_normal_cell(je,jb,2)%v2      

        ptr_nh_prog%vn(je,jk,jb) = ptr_int%c_lin_e(je,1,jb)*zvn1 + &
                                   ptr_int%c_lin_e(je,2,jb)*zvn2

        ptr_nh_ref%vn_ref(je,jk,jb) = ptr_nh_prog%vn(je,jk,jb)
      END DO
     END DO
    END DO     
    
    !W wind and reference
    CALL init_w(ptr_patch, ptr_int, ptr_nh_prog%vn, ptr_metrics%z_ifc, ptr_nh_prog%w)
    CALL sync_patch_array(SYNC_C, ptr_patch, ptr_nh_prog%w)
    ptr_nh_ref%w_ref = ptr_nh_prog%w
    
  END SUBROUTINE init_nh_state_rico
  
!-------------------------------------------------------------------------


  SUBROUTINE cbl_stevens_fluxes( t_l, qv_l, p_l, rho_l, tsk, shfx, lhfx )

  !-------------------------------------------------------------------------
  ! Calculate sensible and latent heat fluxes from buoyancy flux for Stevens 
  ! (2007) case. (code from Wayne Angevine 2013) 
  !
  ! Variable explanations
  ! tsk  = skin temperature
  ! zp0  = surface pressure (Pa)
  ! qv1d = water vapor mixing ratio in lowest level
  ! th1d = potential temperature in lowest level
  ! kqfx = exchange coefficient for water vapor (used in qfx calculation below)
  !
  ! Constants for saturation calculation:
  ! svp1  = 0.6112
  ! svp2  = 17.67
  ! svp3  = 29.65
  !-------------------------------------------------------------------------

    ! INPUT PARAMETERS:
    REAL(wp), INTENT(IN) ::    t_l       , & ! temperature at lowest level [K]
     &                         qv_l      , & ! moisture at lowest level    [kg/kg]
     &                         p_l       , & ! pressure at lowest level    [Pa]
     &                         rho_l         ! rho at lowest level         [kg/m3]

    ! INPUT/OUTPUT PARAMETERS:
    REAL(wp), INTENT(INOUT) :: tsk           ! skin temperature            [K]

    ! OUTPUT PARAMETERS:
    REAL(wp), INTENT(OUT) ::   shfx      , & ! sensible heat flux (+ down) [W/m2]
      &                        lhfx          ! latent heat flux   (+ down) [W/m2]

    ! LOCAL VARIABLES:
    REAL(wp) :: Beta0, Vs, mav, qsfc, qsfc_air, kqfx, khfx, th_l

    REAL(wp), PARAMETER :: zp0     = 101540._wp !p0ref !< surface pressure

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
    Beta0 = 7e-4_wp       ! Fixed specified buoyancy flux (25 W/m^2)
   !Beta0 = 4.2e-4_wp     ! Fixed specified buoyancy flux (15 W/m^2)
   !Beta0 = 11.2e-4_wp    ! Fixed specified buoyancy flux (40 W/m^2)
    Vs    = 0.01_wp       ! Specified velocity scale
    mav   = 0.90_wp       ! Surface moisture availability

  ! Calculate saturation mixing ratio at SST (from previous timestep)
  ! The rest of the buoyancy flux goes to sensible heat
  ! Calculate surface saturated q and q in air at surface 
   !e1=svp1*exp(svp2*(tsk-tmelt)/(tsk-svp3))                       
   !qsfc=rd/rv*e1/((zp0/1000.)-e1)
    qsfc     = spec_humi(sat_pres_water(tsk),zp0)
    qsfc_air = qsfc * mav

    th_l =  t_l * (p0ref/p_l)**rd_o_cpd 

    ! Calculate hfx,qfx, and SST to keep buoyancy flux constant
    ! Could calculate moisture flux first, but should be OK either way
    kqfx = Vs * (qsfc_air - qv_l)   ! exchange coefficient for water vapor
    khfx = Beta0 * th_l/grav - 0.608_wp * th_l * kqfx
    tsk  = khfx / Vs + th_l

    ! Convert units
    shfx = - rho_l * cpd * khfx
    lhfx = - rho_l * alv * kqfx

  END SUBROUTINE cbl_stevens_fluxes

!-------------------------------------------------------------------------
!
! This subroutine creates a simple two valued field for the sensible heat flux
! and the water vapor flux.  The domain is simply divided in two with the 
! division determined by the longitude value given.  on each side of the 
! division the sensible and latent heat fluxes have different values.
!
  SUBROUTINE sfcflx_uniform(ptr_patch, shflux_sfc, sh_high, sh_low, qvflux_sfc,   &
    & qv_high, qv_low, wallLonDeg)
    TYPE(t_patch),TARGET,  INTENT(IN) :: ptr_patch  !< patch on which computation is performed
    REAL(wp), INTENT(in) :: wallLonDeg
    REAL(wp), INTENT(out):: shflux_sfc(:,:) ! sensible heat flux [W/m2]
    REAL(wp), INTENT(out):: qvflux_sfc(:,:) ! Water vapor flux at sfc [kg/kg]
    Real(wp), INTENT(in) :: sh_high   ! upper value of sensible heat flux
    Real(wp), INTENT(in) :: sh_low    ! lower value of sensible heat flux
    Real(wp), INTENT(in) :: qv_high   ! upper value of vapor flux
    Real(wp), INTENT(in) :: qv_low    ! lower value of vapor flux

    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: jb, jc
    INTEGER :: start_cell_index, end_cell_index
    REAL(wp):: lon_deg

!-------------------------------------------------------------------------

    all_cells => ptr_patch%cells%ALL

    !Add horizontal variation
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index

        lon_deg = ptr_patch%cells%center(jc,jb)%lon * rad2deg

        IF((lon_deg) >= wallLonDeg) THEN
          shflux_sfc(jc,jb) = sh_high
          qvflux_sfc(jc,jb) = qv_high
        ELSE
          shflux_sfc(jc,jb) = sh_low
          qvflux_sfc(jc,jb) = qv_low
        ENDIF

      END DO
    END DO

  END SUBROUTINE sfcflx_uniform


  !-------------------------------------------------------------------------
  !>
  ! read sounding from external file and then interpolate
  ! to model levels 
  ! Sounding file is assumed to be in this format:
  ! ps,no_vert_levels
  ! z(m) theta(k) qv(kg/kg) u(m/s) v(m/s)
  ! where the top row is close to surface/or surface and bottom row
  ! is near the top
  !
  SUBROUTINE read_ext_profile (z_in, theta_in, qv_in, u_in, v_in, psfc_in)
  
    REAL(wp),  INTENT(IN)  :: z_in(:)
    REAL(wp),  INTENT(OUT) :: theta_in(:)
    REAL(wp),  INTENT(OUT) :: qv_in(:)
    REAL(wp),  INTENT(OUT) :: u_in(:)
    REAL(wp),  INTENT(OUT) :: v_in(:)
    REAL(wp),  INTENT(OUT) :: psfc_in
  
    REAL(wp), ALLOCATABLE, DIMENSION(:):: zs, ths, qvs, us, vs
    CHARACTER(len=*),PARAMETER :: routine  = &
         &   'mo_nh_torus_exp:read_ext_profile'
  
    INTEGER :: ist, iunit
    INTEGER :: jk, klev 
    
    !-------------------------------------------------------------------------
  
    CALL message(routine, 'READING FROM SOUNDING!')
    
    !open file again to read data this time
    iunit = find_next_free_unit(10,100)
    OPEN (unit=iunit,file='sound_in', access='SEQUENTIAL', &
            form='FORMATTED',action='READ', status='OLD', IOSTAT=ist) 
    IF(ist/=success)THEN
      CALL finish (routine, 'open verticaling sound file failed')
    ENDIF
  
    !Read the header : ps,klev
    READ (iunit,*,IOSTAT=ist)psfc_in,klev

    !now allocate
    ALLOCATE(zs(klev),ths(klev),us(klev),vs(klev),qvs(klev))
    zs = 0.0_wp; ths = 0._wp; us = 0._wp; vs = 0._wp; qvs = 0._wp

    DO jk = klev,1,-1 
      READ (iunit,*,IOSTAT=ist) zs(jk),ths(jk),qvs(jk),us(jk),vs(jk)
      IF(ist/=success)THEN
        CALL finish (routine, 'reading sounding file failed')
      ENDIF
    END DO

    CLOSE(iunit)

    !Check if the file is written in descending order
    IF(zs(1) < zs(klev)) &
         CALL finish (routine, 'Writing souding data in descending order!')

    !Now perform interpolation to grid levels assuming:
    !a) linear interpolation
    !b) Beyond the last Z level the values are linearly extrapolated 
    !c) Assuming model grid is flat-NOT on sphere
  
    CALL vert_intp_linear_1d(zs,ths,z_in,theta_in)
    CALL vert_intp_linear_1d(zs,qvs,z_in,qv_in)
    CALL vert_intp_linear_1d(zs,us,z_in,u_in)
    CALL vert_intp_linear_1d(zs,vs,z_in,v_in)
  
    DEALLOCATE(zs, ths, qvs, us, vs)
  
  
  END SUBROUTINE  read_ext_profile


  !-------------------------------------------------------------------------
  !>
  ! read sounding from external file and then interpolate
  ! to model levels 
  ! Sounding file is assumed to be in this format:
  ! ps,no_vert_levels
  ! z(m) theta(k) qv(kg/kg) u(m/s) v(m/s)
  ! where the top row is close to surface/or surface and bottom row
  ! is near the top
  ! IBD: read from netcdf file
  !

  SUBROUTINE  read_ext_profile_nc(z_in, zifc_in, theta_in, thetav_in, exner_in, rho_in, &
    qv_in, qc_in, qi_in, u_in, v_in, w_in, tke_in, psfc_in, o3_in)
  
    REAL(wp),  INTENT(IN)  :: z_in(:)
    REAL(wp),  INTENT(IN)  :: zifc_in(:)
    REAL(wp),  INTENT(OUT) :: theta_in(:)
    REAL(wp),  INTENT(OUT) :: thetav_in(:)
    REAL(wp),  INTENT(OUT) :: exner_in(:)
    REAL(wp),  INTENT(OUT) :: rho_in(:)
    REAL(wp),  INTENT(OUT) :: qv_in(:)
    REAL(wp),  INTENT(OUT) :: qc_in(:)
    REAL(wp),  INTENT(OUT) :: qi_in(:)
    REAL(wp),  INTENT(OUT) :: u_in(:)
    REAL(wp),  INTENT(OUT) :: v_in(:)
    REAL(wp),  INTENT(OUT) :: w_in(:)
    REAL(wp),  INTENT(OUT) :: tke_in(:)
    REAL(wp),  INTENT(OUT) :: psfc_in
    REAL(wp),  INTENT(OUT) :: o3_in(:)
  
    REAL(wp), ALLOCATABLE, DIMENSION(:)      :: zs, zs_ifc, ths, thvs, exners, rhos, qvs, qcs, qis, &
         &    us, vs, ws, tkes, psurfs, o3s
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)    :: tempf, tempf1
    CHARACTER(len=max_char_length),PARAMETER :: routine  = &
         &   'mo_nh_torus_exp:read_ext_profile_nc'
  
    INTEGER :: klev,nt
    INTEGER :: varid
    INTEGER :: fileid     !< id number of netcdf file
    INTEGER :: dimid      !< id number of dimension
    INTEGER :: nf_status,nf_status2  !< return status of netcdf function
  
    !-------------------------------------------------------------------------
  
    CALL message(TRIM(routine), 'READING FROM SOUNDING!')
    
    !open netcdf 
    CALL nf (nf90_open('init_SCM.nc', NF90_NOWRITE, fileid), &
      & TRIM(routine)//'   File init_SCM.nc cannot be opened (profile)') 

    CALL nf (nf90_inq_dimid(fileid, 'lev', dimid), routine)
    CALL nf (nf90_inquire_dimension(fileid, dimid, len = klev), routine)

    CALL nf (nf90_inq_dimid(fileid, 'nt', dimid), routine)
    CALL nf (nf90_inquire_dimension(fileid, dimid, len = nt) , routine)

    WRITE(message_text,'(a,i6,a,i6)') 'SCM read_ext_profile_nc: klev, ', klev, ', nt ', nt
    CALL message (routine,message_text)

    !now allocate
    ALLOCATE(zs(klev), zs_ifc(klev+1), ths(klev), thvs(klev), exners(klev), rhos(klev), us(klev), vs(klev), &
      & ws(klev+1), qvs(klev), qcs(klev), qis(klev), tkes(klev+1), psurfs(nt), o3s(klev),                   &
      & tempf(klev,nt), tempf1(klev+1,nt) )
 
    !initialize to 0
    zs     = 0._wp
    zs_ifc = 0._wp
    ths    = 0._wp
    thvs   = 0._wp
    exners = 0._wp
    rhos   = 0._wp
    us     = 0._wp
    vs     = 0._wp
    ws     = 0._wp
    qvs    = 0._wp
    qcs    = 0._wp
    qis    = 0._wp
    tkes   = 0._wp
    psurfs = 0._wp
    o3s    = 0._wp

    CALL nf (nf90_inq_varid (fileid, 'height', varid), routine)
    CALL nf (nf90_get_var   (fileid, varid, tempf)   , routine)
    zs=tempf(:,1)

    CALL nf (nf90_inq_varid (fileid, 'uIN', varid), routine)
    CALL nf (nf90_get_var   (fileid, varid, tempf), routine)
    us=tempf(:,1)

    CALL nf (nf90_inq_varid (fileid, 'vIN', varid), routine)
    CALL nf (nf90_get_var   (fileid, varid, tempf), routine)
    vs=tempf(:,1)

    CALL nf (nf90_inq_varid (fileid, 'tkeIN', varid), routine)
    CALL nf (nf90_get_var   (fileid, varid, tempf1, count = [klev, nt]), routine)
    tkes=tempf1(:,1)

    CALL nf (nf90_inq_varid (fileid, 'qvIN', varid), routine)
    CALL nf (nf90_get_var   (fileid, varid, tempf) , routine)
    qvs=tempf(:,1)

    CALL nf (nf90_inq_varid (fileid, 'qcIN', varid), routine)
    CALL nf (nf90_get_var   (fileid, varid, tempf) , routine)
    qcs=tempf(:,1)

    CALL nf (nf90_inq_varid (fileid, 'psurf', varid), routine)
    CALL nf (nf90_get_var   (fileid, varid, psurfs) , routine)
    psfc_in=psurfs(1)

    nf_status  = nf90_inq_varid (fileid, 'o3IN', varid)
    nf_status2 = nf90_get_var   (fileid, varid , tempf)
    IF (nf_status /= nf90_noerr) THEN
      CALL message (routine,'O3 not available in init_SCM.nc.  It will be set to 0.')
      o3s=0.0_wp
    ELSE
      o3s=tempf(:,1)
    END IF

    !read input variables specific to ICON4SCM
    IF (lscm_icon_ini) THEN
      CALL nf (nf90_inq_varid (fileid, 'height_ifc', varid), routine)
      CALL nf (nf90_get_var   (fileid, varid, tempf1, count = [klev, nt]), routine)
      zs_ifc = tempf1(:,1)

      CALL nf (nf90_inq_varid (fileid, 'wIN', varid), routine)
      CALL nf (nf90_get_var   (fileid, varid,tempf1, count = [klev, nt]), routine)
      ws     = tempf1(:,1)

      CALL nf (nf90_inq_varid (fileid, 'qiIN', varid), routine)
      CALL nf (nf90_get_var   (fileid, varid,tempf)  , routine)
      qis    = tempf(:,1)

      CALL nf (nf90_inq_varid (fileid, 'thvIN', varid), routine)
      CALL nf (nf90_get_var   (fileid, varid,tempf)   , routine)
      thvs   = tempf(:,1)

      CALL nf (nf90_inq_varid (fileid, 'exnerIN', varid), routine)
      CALL nf (nf90_get_var   (fileid, varid,tempf)     , routine)
      exners = tempf(:,1)

      CALL nf (nf90_inq_varid (fileid, 'rhoIN', varid), routine)
      CALL nf (nf90_get_var   (fileid, varid,tempf)   , routine)
      rhos   = tempf(:,1)
    ELSE
      zs_ifc = zs

      CALL nf (nf90_inq_varid (fileid, 'thIN', varid), routine)
      CALL nf (nf90_get_var   (fileid, varid,tempf)  , routine)
      ths    = tempf(:,1)
    END IF

!   diagnostic output
!   write(*,*) '==SCM input file data=='
!   write(*,*) '-- zs --'    , zs
!   write(*,*) '-- zs_ifc --', zs_ifc
!   write(*,*) '-- z_in --'  , z_in
!   write(*,*) '-- zifc_in --',zifc_in
!   write(*,*) '-- tkes --'  , tkes
!   write(*,*) '-- qvs --'   , qvs
!   write(*,*) '-- qcs --'   , qcs
!   write(*,*) '-- psurfs --', psurfs
!   write(*,*) '-- o3 --'    , o3s
!   write(*,*) '-- ws --'    , ws
!   write(*,*) '-- qis --'   , qis
!   write(*,*) '-- ths --'   , ths
!   write(*,*) '-- thvs --'  , thvs
!   write(*,*) '-- exners --', exners
!   write(*,*) '-- rhos --'  , rhos
!   write(*,*) '-- us --'    , us
!   write(*,*) '-- vs --'    , vs

    !Check if the file is written in descending order
    IF(zs(1) < zs(klev)) THEN
         CALL finish (TRIM(routine), 'Writing souding data in descending order!')
    ENDIF

    !Now perform interpolation to grid levels assuming:
    !a) linear interpolation
    !b) Beyond the last Z level the values are linearly extrapolated 
    !c) Assuming model grid is flat-NOT on sphere


    CALL vert_intp_linear_1d(zs,     ths,    z_in,    theta_in)
    CALL vert_intp_linear_1d(zs,     thvs,   z_in,    thetav_in)
    CALL vert_intp_linear_1d(zs,     exners, z_in,    exner_in)
    CALL vert_intp_linear_1d(zs,     rhos,   z_in,    rho_in)
    CALL vert_intp_linear_1d(zs,     qvs,    z_in,    qv_in)
    CALL vert_intp_linear_1d(zs,     qcs,    z_in,    qc_in)
    CALL vert_intp_linear_1d(zs,     qis,    z_in,    qi_in)
    CALL vert_intp_linear_1d(zs,     us,     z_in,    u_in)
    CALL vert_intp_linear_1d(zs,     vs,     z_in,    v_in)
    CALL vert_intp_linear_1d(zs_ifc, ws,     zifc_in, w_in)
    CALL vert_intp_linear_1d(zs_ifc, tkes,   zifc_in, tke_in)
    CALL vert_intp_linear_1d(zs,     o3s,    z_in,    o3_in)
  
    DEALLOCATE(zs, ths, thvs, exners, rhos, qvs, qcs, qis, us, vs, ws, tkes, &
    & psurfs, o3s, tempf, tempf1)
  
    CALL nf (nf90_close(fileid), routine)

    CALL read_latlon_scm_nc(lat_scm,lon_scm)
  
  END SUBROUTINE  read_ext_profile_nc


  !-------------------------------------------------------------------------
  !>
  ! read sounding from external file and then interpolate
  ! to model levels (unified SCM format)
  !

  SUBROUTINE  read_ext_profile_nc_uf(z_in, zifc_in, theta_in, thetav_in, exner_in, rho_in, &
    qv_in, qc_in, qi_in, u_in, v_in, w_in, tke_in, psfc_in, o3_in)
  
    REAL(wp),  INTENT(IN)  :: z_in(:)
    REAL(wp),  INTENT(IN)  :: zifc_in(:)
    REAL(wp),  INTENT(OUT) :: theta_in(:)
    REAL(wp),  INTENT(OUT) :: thetav_in(:)
    REAL(wp),  INTENT(OUT) :: exner_in(:)
    REAL(wp),  INTENT(OUT) :: rho_in(:)
    REAL(wp),  INTENT(OUT) :: qv_in(:)
    REAL(wp),  INTENT(OUT) :: qc_in(:)
    REAL(wp),  INTENT(OUT) :: qi_in(:)
    REAL(wp),  INTENT(OUT) :: u_in(:)
    REAL(wp),  INTENT(OUT) :: v_in(:)
    REAL(wp),  INTENT(OUT) :: w_in(:)
    REAL(wp),  INTENT(OUT) :: tke_in(:)
    REAL(wp),  INTENT(OUT) :: psfc_in
    REAL(wp),  INTENT(OUT) :: o3_in(:)
  
    REAL(wp), ALLOCATABLE, DIMENSION(:)       :: zs, zs_ifc, ths, thvs, exners, rhos, qvs, qcs, qis, &
         &    us, vs, ws, tkes, o3s
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: tempf, tempf1
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: tempf_s
    CHARACTER(len=max_char_length),PARAMETER  :: routine  = &
         &   'mo_nh_torus_exp:read_ext_profile_nc_uf'
  
    INTEGER :: klev,nt
    INTEGER :: lat,lon,t0
    INTEGER :: varid
    INTEGER :: fileid     !< id number of netcdf file
    INTEGER :: dimid      !< id number of dimension
    INTEGER :: nf_status,nf_status2  !< return status of netcdf function

    !-------------------------------------------------------------------------
  
    CALL message(TRIM(routine), 'READING FROM SOUNDING!')
    
    !open netcdf 
    CALL nf (nf90_open('init_SCM.nc', NF90_NOWRITE, fileid), &
      & TRIM(routine)//'   File init_SCM.nc cannot be opened (profile)') 

    CALL nf (nf90_inq_dimid(fileid, 'lev', dimid), routine)
    CALL nf (nf90_inquire_dimension(fileid, dimid, len = klev), routine)

    CALL nf (nf90_inq_dimid(fileid, 'lat', dimid), routine)
    CALL nf (nf90_inquire_dimension(fileid, dimid, len = lat) , routine)

    CALL nf (nf90_inq_dimid(fileid, 'lon', dimid), routine)
    CALL nf (nf90_inquire_dimension(fileid, dimid, len = lon) , routine)

    CALL nf (nf90_inq_dimid(fileid, 't0', dimid), routine)
    CALL nf (nf90_inquire_dimension(fileid, dimid, len = t0) , routine)

    CALL nf (nf90_inq_dimid(fileid, 'time', dimid), routine)
    CALL nf (nf90_inquire_dimension(fileid, dimid, len = nt)   , routine)
    WRITE(message_text,'(a,i5,a,i5,a,i5)') 'SCM read_ext_profile_nc_uf: klev, ', klev, ', nt ', nt, ', t0 ', t0
    CALL message (routine,message_text)

    !now allocate
    ALLOCATE(zs(klev), zs_ifc(klev), ths(klev), thvs(klev), exners(klev), rhos(klev), us(klev), vs(klev), &
      & ws(klev), qvs(klev), qcs(klev), qis(klev), tkes(klev),o3s(klev),                                &
      & tempf(lon,lat,klev,t0), tempf1(lon,lat,klev,t0), tempf_s(lon,lat,t0) )
 
    !initialize to 0
    zs     = 0._wp
    zs_ifc = 0._wp
    ths    = 0._wp
    thvs   = 0._wp
    exners = 0._wp
    rhos   = 0._wp
    us     = 0._wp
    vs     = 0._wp
    ws     = 0._wp
    qvs    = 0._wp
    qcs    = 0._wp
    qis    = 0._wp
    tkes   = 0._wp
    psfc_in= 0._wp
    o3s    = 0._wp

    CALL nf (nf90_inq_varid (fileid, 'height', varid), routine)
    CALL nf (nf90_get_var   (fileid, varid, tempf)   , routine)
    zs=tempf(1,1,klev:1:-1,1)

    CALL nf (nf90_inq_varid (fileid, 'u', varid)  , routine)
    CALL nf (nf90_get_var   (fileid, varid, tempf), routine)
    us=tempf(1,1,klev:1:-1,1)

    CALL nf (nf90_inq_varid (fileid, 'v', varid)  , routine)
    CALL nf (nf90_get_var   (fileid, varid, tempf), routine)
    vs=tempf(1,1,klev:1:-1,1)

    CALL nf (nf90_inq_varid (fileid, 'qv', varid) , routine)
    CALL nf (nf90_get_var   (fileid, varid, tempf), routine)
    qvs=tempf(1,1,klev:1:-1,1)

    CALL nf (nf90_inq_varid (fileid, 'ql', varid) , routine)
    CALL nf (nf90_get_var   (fileid, varid, tempf), routine)
    qcs=tempf(1,1,klev:1:-1,1)

    CALL nf (nf90_inq_varid (fileid, 'ps', varid)   , routine)
    CALL nf (nf90_get_var   (fileid, varid, tempf_s), routine)
    psfc_in=tempf_s(1,1,1)

    nf_status  = nf90_inq_varid (fileid, 'o3', varid)
    nf_status2 = nf90_get_var   (fileid, varid , tempf)
    IF (nf_status /= nf90_noerr) THEN
      CALL message (routine,'O3 not available in init_SCM.nc.  It will be set to 0.')
      o3s=0.0_wp
    ELSE
      o3s=tempf(1,1,klev:1:-1,1)
    END IF

    !read input variables specific to ICON4SCM
    IF (lscm_icon_ini) THEN
      zs_ifc=zs

      CALL nf (nf90_inq_varid (fileid, 'w', varid), routine)
      CALL nf (nf90_get_var   (fileid, varid, tempf1, count = [klev, nt]), routine)
      ws     = tempf1(1,1,klev:1:-1,1)

      CALL nf (nf90_inq_varid (fileid, 'qi', varid), routine)
      CALL nf (nf90_get_var(fileid, varid, tempf), routine)
      qis    = tempf(1,1,klev:1:-1,1)

      CALL nf (nf90_inq_varid (fileid, 'thetav', varid), routine)
      CALL nf (nf90_get_var   (fileid, varid, tempf), routine)
      thvs   = tempf(1,1,klev:1:-1,1)

      CALL nf (nf90_inq_varid (fileid, 'exner', varid), routine)
      CALL nf (nf90_get_var   (fileid, varid, tempf), routine)
      exners = tempf(1,1,klev:1:-1,1)

      CALL nf (nf90_inq_varid (fileid, 'rho', varid), routine)
      CALL nf (nf90_get_var   (fileid, varid, tempf), routine)
      rhos   = tempf(1,1,klev:1:-1,1)

      CALL nf (nf90_inq_varid (fileid, 'tke', varid), routine)
      CALL nf (nf90_get_var   (fileid, varid, tempf1, count = [klev, nt]), routine)
      tkes   = tempf1(1,1,klev:1:-1,1)

    ELSE
      zs_ifc = zs

      CALL nf (nf90_inq_varid (fileid, 'theta', varid), routine)
      CALL nf (nf90_get_var (fileid, varid,tempf), routine)
      ths    = tempf(1,1,klev:1:-1,1)

      CALL nf (nf90_inq_varid (fileid, 'tke', varid), routine)
      CALL nf (nf90_get_var (fileid, varid, tempf), routine)
      tkes(1:klev) = tempf(1,1,klev:1:-1,1)

    END IF

!   diagnostic output
!   write(*,*) '==SCM input file data=='
!   write(*,*) '-- zs --'     , zs
!   write(*,*) '-- zs_ifc --' , zs_ifc
!   write(*,*) '-- tkes --'   , tkes
!   write(*,*) '-- qvs --'    , qvs
!   write(*,*) '-- qcs --'    , qcs
!   write(*,*) '-- psfc_in --', psfc_in
!   write(*,*) '-- o3 --'     , o3s
!   write(*,*) '-- ws --'     , ws
!   write(*,*) '-- qis --'    , qis
!   write(*,*) '-- ths --'    , ths
!   write(*,*) '-- thvs --'   , thvs
!   write(*,*) '-- exners --' , exners
!   write(*,*) '-- rhos --'   , rhos
!   write(*,*) '-- us --'     , us
!   write(*,*) '-- vs --'     , vs

    !Check if the file is written in descending order
    IF(zs(1) < zs(klev)) THEN
         CALL finish (TRIM(routine), 'Writing souding data in descending order!')
    ENDIF

    !Now perform interpolation to grid levels assuming:
    !a) linear interpolation
    !b) Beyond the last Z level the values are linearly extrapolated 
    !c) Assuming model grid is flat-NOT on sphere
  
    CALL vert_intp_linear_1d(zs,     ths,    z_in,    theta_in)
    CALL vert_intp_linear_1d(zs,     thvs,   z_in,    thetav_in)
    CALL vert_intp_linear_1d(zs,     exners, z_in,    exner_in)
    CALL vert_intp_linear_1d(zs,     rhos,   z_in,    rho_in)
    CALL vert_intp_linear_1d(zs,     qvs,    z_in,    qv_in)
    CALL vert_intp_linear_1d(zs,     qcs,    z_in,    qc_in)
    CALL vert_intp_linear_1d(zs,     qis,    z_in,    qi_in)
    CALL vert_intp_linear_1d(zs,     us,     z_in,    u_in)
    CALL vert_intp_linear_1d(zs,     vs,     z_in,    v_in)
    CALL vert_intp_linear_1d(zs_ifc, ws,     zifc_in, w_in)
    CALL vert_intp_linear_1d(zs_ifc, tkes,   zifc_in, tke_in)
    CALL vert_intp_linear_1d(zs,     o3s,    z_in,    o3_in)
  
!   write(*,*) '==SCM input file data interpolated to SCM vertical levels=='
!   write(*,*) '-- o3_in --'    , o3_in

    DEALLOCATE(zs, ths, thvs, exners, rhos, qvs, qcs, qis, us, vs, ws, tkes, &
    & o3s, tempf, tempf1, tempf_s)

    CALL nf (nf90_close(fileid), routine)
  
    CALL read_latlon_scm_nc_uf(lat_scm,lon_scm)

  END SUBROUTINE  read_ext_profile_nc_uf


  !>
  !! Initialization of prognostic state vector for the analytical RCEMIP test case. 
  SUBROUTINE init_torus_rcemip_analytical_sounding ( ptr_patch, ptr_nh_prog,  ptr_nh_ref, ptr_nh_diag,  &
                                        ptr_int, ptr_metrics)

    TYPE(t_patch),TARGET,  INTENT(INOUT)::  ptr_patch
    TYPE(t_int_state),     INTENT(IN)   ::  ptr_int
    TYPE(t_nh_prog),       INTENT(INOUT)::  ptr_nh_prog
    TYPE(t_nh_diag),       INTENT(INOUT)::  ptr_nh_diag
    TYPE(t_nh_metrics),    INTENT(IN)   ::  ptr_metrics
    TYPE(t_nh_ref),        INTENT(INOUT)::  ptr_nh_ref

    INTEGER  :: je,jk,jb,i_startidx,i_endidx   !< loop indices
   ! INTEGER  :: nblks_c,npromz_c,nblks_e,npromz_e, jg
    INTEGER  :: nblks_c,npromz_c,nblks_e, jg
    INTEGER  :: nlev, nlevp1                        !< number of full and half levels
    INTEGER  :: nlen !, i_rcstartlev, jcn, jbn, ist

    REAL(wp), DIMENSION(ptr_patch%nlev) :: hght
    REAL(wp) :: zvn1, zvn2, zu, zv
    REAL(wp) :: tv0, tvi, tvt, pr, prt, ex_sfc, ex, qv0, thv, qv, o3, irho
    REAL(wp) :: z_exner_h(1:nproma,ptr_patch%nlev+1), z_help(1:nproma)

    REAL(wp), PARAMETER :: laps  = 0.0067_wp  ! Lapse rate [K/m]
    REAL(wp), PARAMETER :: ztop  = 15000._wp  ! z_trop [m]
    REAL(wp), PARAMETER :: qvtop = 1.e-14_wp  ! q-top [kg/kg]
    REAL(wp), PARAMETER :: zqv1  = 4000._wp   ! z_q1 setting [m]
    REAL(wp), PARAMETER :: zqv2  = 7500._wp   ! z_q1 setting [m]
    REAL(wp), PARAMETER :: pr0   = 101480._wp ! sfc pressure [Pa]

    CHARACTER(len=*), PARAMETER :: &
       &  routine = 'mo_nh_torus_exp:init_torus_rcemip_analytical_sounding'
  !-------------------------------------------------------------------------

    ! values for the blocking
    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c
    nblks_e  = ptr_patch%nblks_e
    ! npromz_e = ptr_patch%npromz_e

    ! number of vertical levels
    nlev   = ptr_patch%nlev
    nlevp1 = ptr_patch%nlevp1

    !patch id
    jg = ptr_patch%id

    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = pr0
    ex_sfc   = (pr0/p0ref)**rd_o_cpd
    les_config(jg)%psfc = pr0

    ! height on mass levels
    hght = ptr_metrics%z_mc(2,:,2)

    ! Tracers: all zero by default
    ptr_nh_prog%tracer(:,:,:,:) = 0._wp

  ! sst equal 300
  !   IF (les_config(jg)%sst .eq. 295) THEN
  !     qv0 = 12.0e-3_wp ! kg/kg
  !   ELSE IF (les_config(jg)%sst .eq. 300) THEN
      qv0 = 18.65e-3_wp ! kg/kg
  !   ELSE IF (les_config(jg)%sst .eq. 305) THEN
  !     qv0 = 24.0e-3_wp ! kg/kg
  !   ELSE
  !     CALL finish(TRIM(routine),'No preset qv0 scernario for this SST!')
  !   ENDIF

  !   tv0 = les_config(jg)%sst * (1._wp + qv0*0.608_wp)
  ! sst equal 300
    tv0 = 300._wp * (1._wp + qv0*0.608_wp)
    tvt = tv0 - laps*ztop
    prt = pr0 * (tvt/tv0)**(grav/(rd*laps))

    DO jb = 1, nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      DO jk = 1, nlev

        IF (hght(jk) .le. ztop) THEN
          qv  = qv0 * exp(-hght(jk)/zqv1) * &
                      exp(-((hght(jk)/zqv2)**2))
          tvi = tv0 - laps*hght(jk)
          pr  = pr0 * (tvi/tv0)**(grav/(rd*laps))
        ELSE
          qv  = qvtop
          tvi = tvt
          pr  = prt * exp(-(grav*(hght(jk)-ztop)/(rd*tvt)))
        ENDIF

        thv = tvi * (p0ref/pr)**(rd/cpd)
        ptr_nh_prog%theta_v(1:nlen,jk,jb) = thv
        ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) = qv
        ex = (pr/p0ref)**(rd/cpd)
        ptr_nh_prog%exner(1:nlen,jk,jb) = ex

        irho = (ex**cvd_o_rd)*p0ref/rd / thv
        ptr_nh_prog%rho(1:nlen,jk,jb) = irho

      END DO !jk

    ENDDO !jb

    ! Write out profile
    IF(my_process_is_stdio()) THEN
      WRITE(0,*) 'sfc pres (hPa)',ptr_nh_diag%pres_sfc(1,1)*1.e-2
      WRITE(0,*) 'hght (m), theta_v (K), qv (g/kg)'
      DO jk = 1, nlev
        WRITE(0,*) ptr_metrics%z_mc(2,jk,2), ptr_nh_prog%theta_v(1,jk,1), ptr_nh_prog%tracer(1,jk,1,iqv)*1.e3
      ENDDO
    ENDIF

    !Mean wind 
    DO jb = 1 , nblks_e
     CALL get_indices_e( ptr_patch, jb, 1, nblks_e, i_startidx, i_endidx, grf_bdywidth_e+1)
     DO jk = 1 , nlev
      DO je = i_startidx, i_endidx

     !   jcn  =   ptr_patch%edges%cell_idx(je,jb,1)
     !   jbn  =   ptr_patch%edges%cell_blk(je,jb,1)
        zu   =   0._wp!u_in(jk)
        zv   =   0._wp!v_in(jk)

        zvn1 =  zu * ptr_patch%edges%primal_normal_cell(je,jb,1)%v1 + &
                zv * ptr_patch%edges%primal_normal_cell(je,jb,1)%v2

     !   jcn  =   ptr_patch%edges%cell_idx(je,jb,2)
     !   jbn  =   ptr_patch%edges%cell_blk(je,jb,2)
        zu   =   0._wp!u_in(jk)
        zv   =   0._wp!v_in(jk)

        zvn2 =  zu * ptr_patch%edges%primal_normal_cell(je,jb,2)%v1 + &
                zv * ptr_patch%edges%primal_normal_cell(je,jb,2)%v2

        ptr_nh_prog%vn(je,jk,jb) = ptr_int%c_lin_e(je,1,jb)*zvn1 + &
                                   ptr_int%c_lin_e(je,2,jb)*zvn2

        ptr_nh_ref%vn_ref(je,jk,jb) = ptr_nh_prog%vn(je,jk,jb)
      END DO
     END DO
    END DO

    !W wind and reference
    CALL init_w(ptr_patch, ptr_int, ptr_nh_prog%vn, ptr_metrics%z_ifc, ptr_nh_prog%w)
    CALL sync_patch_array(SYNC_C, ptr_patch, ptr_nh_prog%w)
    ptr_nh_ref%w_ref  = ptr_nh_prog%w

  END SUBROUTINE init_torus_rcemip_analytical_sounding

  !>
  !! Initialization of prognostic state vector for the warm bubble experiment
  !! on a torus. It is simplified form of Wesimann Klemp testcase from
  !! G. H. Bryan and J. M. Fritsch, "A benchmark simulation for moist
  !! nonhydrostatic numerical models", MWR, 2002
  !!
  SUBROUTINE init_warm_bubble( ptr_patch, ptr_nh_prog, ptr_nh_diag, ptr_metrics)

    TYPE(t_patch),TARGET,  INTENT(IN)   ::  ptr_patch
    TYPE(t_nh_prog),       INTENT(INOUT)::  ptr_nh_prog
    TYPE(t_nh_diag),       INTENT(INOUT)::  ptr_nh_diag
    TYPE(t_nh_metrics),    INTENT(IN)   ::  ptr_metrics      

    REAL(wp) :: z_exner_h(1:nproma,ptr_patch%nlev+1), z_help(1:nproma) 
    REAL(wp) :: ex_sfc, x_loc(3), x_c(3), psfc_in, dis, inv_th0, pres_new
    REAL(wp) :: th_new, qv_new, qc_new, th_ptb, temp_new
    REAL(wp), DIMENSION(ptr_patch%nlev) :: theta_in, qv_in, qc_in, tmp
    INTEGER  :: jc,jk,jb   !< loop indices
    INTEGER  :: nblks_c,npromz_c
    INTEGER  :: nlev, nlevp1
    LOGICAL  :: qc_fail, is_2d_bubble
    INTEGER  :: nlen, jg, itr

    REAL(wp), DIMENSION(3) :: x_bubble 
    CHARACTER(len=*),PARAMETER :: routine  = &
         &   'mo_nh_torus_exp:init_warm_bubble'
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    !Note that this souding is created from a matlab code that 
    !iterates through the Eqs. 25,26,and 34 of Bryan and Fritsch's paper
    !given theta_e, qt, and surface pressure. The source code is 
    !in icon-aes-and/scripts/preprocessing/ named init.f90 and findzero.m
    !which creates a sound_** file that is then used to read in below
    !
    !One can also use the fortran source code from Bryan's website to generate
    !the initial condition
    !-------------------------------------------------------------------------
    !Read the sounding file
    CALL  read_ext_profile (ptr_metrics%z_mc(1,:,1), theta_in, qv_in, qc_in, tmp, &
                            psfc_in)
  
    ! values for the blocking
    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c

    ! number of vertical levels
    nlev   = ptr_patch%nlev
    nlevp1 = ptr_patch%nlevp1

    !patch id
    jg = ptr_patch%id

    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = psfc_in
    ex_sfc   = (psfc_in/p0ref)**rd_o_cpd
    les_config(jg)%psfc = psfc_in

    ! Tracers: all zero by default
    ptr_nh_prog%tracer(:,:,:,:) = 0._wp

    DO jb = 1, nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      DO jk = 1, nlev
        ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) = qv_in(jk) 
        ptr_nh_prog%tracer(1:nlen,jk,jb,iqc) = qc_in(jk)
      END DO

      DO jk = nlev, 1, -1
         z_help(1:nlen) = theta_in(jk)
         ptr_nh_prog%theta_v(1:nlen,jk,jb) = z_help(1:nlen) * ( 1._wp + &
           0.61_wp*ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) - ptr_nh_prog%tracer(1:nlen,jk,jb,iqc) ) 
      END DO

      !Get hydrostatic exner at the surface using surface pressure 
      z_exner_h(1:nlen,nlevp1) = ex_sfc
 
      !Get exner at full levels starting from exner at surface
      DO jk = nlev, 1, -1
         !exner at next half level after surface
         z_exner_h(1:nlen,jk) = z_exner_h(1:nlen,jk+1) - grav/cpd *     &
                                ptr_metrics%ddqz_z_full(1:nlen,jk,jb)/ &
                                ptr_nh_prog%theta_v(1:nlen,jk,jb)
        
         !exner at main levels
         ptr_nh_prog%exner(1:nlen,jk,jb) = 0.5_wp * &
                                     (z_exner_h(1:nlen,jk)+z_exner_h(1:nlen,jk+1))
      END DO

    ENDDO !jb

    !Zero wind
    ptr_nh_prog%vn = 0._wp
    ptr_nh_prog%w  = 0._wp

    !--------------------------------------------------------------------------
    !Add perturbation to theta_v and iterate to get balanced thermodynamic state
    !--------------------------------------------------------------------------

    !Bubble center, note that torus domain has center in the middle
    !Uses bubctr_lon and bubctr_lat to represent x,y in torus
    x_bubble = (/bubctr_x,bubctr_y,bubctr_z/)

    !First non-dimensionalize bubble ceter
    x_c(1) = x_bubble(1) / bub_hor_width
    x_c(2) = x_bubble(2) / bub_hor_width
    x_c(3) = x_bubble(3) / bub_ver_width


    !the th0 ratio
    inv_th0 = 1._wp / 300._wp

    qc_fail = .FALSE.
    is_2d_bubble = nh_test_name == '2D_BUBBLE'


    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
      DO jc = 1 , nlen
        DO jk = 1 , nlev
          x_loc(1) = ptr_patch%cells%cartesian_center(jc,jb)%x(1)/bub_hor_width
          x_loc(2) = ptr_patch%cells%cartesian_center(jc,jb)%x(2)/bub_hor_width
          x_loc(3) = ptr_metrics%z_mc(jc,jk,jb)/bub_ver_width
            
          IF (is_2d_bubble) THEN
           x_c(2)   = x_loc(2)
          END IF

          dis = plane_torus_distance(x_loc,x_c,ptr_patch%geometry_info)

          IF(dis < 1._wp)THEN
            th_ptb =  bub_amp * cos(pi_2*dis)**2 * inv_th0
            qv_new = qv_in(jk)
            qc_new = qc_in(jk)
            pres_new = p0ref*ptr_nh_prog%exner(jc,jk,jb)**(cpd/rd)
  
            DO itr = 1 , 20 
             th_new = ( th_ptb + 1._wp ) * ptr_nh_prog%theta_v(jc,jk,jb) / &
                       (1._wp + vtmpc1 * qv_new - qc_new)
             temp_new = th_new * ptr_nh_prog%exner(jc,jk,jb)
             qv_new   = spec_humi(sat_pres_water(temp_new),pres_new)
             qc_new   = qv_in(jk) + qc_in(jk) - qv_new
             qc_fail = qc_fail .OR. qc_new < 0._wp

            END DO

            !assign values to proper prog vars
            ptr_nh_prog%tracer(jc,jk,jb,iqv) = qv_new 
            ptr_nh_prog%tracer(jc,jk,jb,iqc) = qc_new
            ptr_nh_prog%theta_v(jc,jk,jb) = th_new * ( 1._wp + vtmpc1*qv_new - qc_new ) 

          END IF
            
        END DO
      END DO
    END DO 

    if (qc_fail) CALL finish(routine, 'qc < 0')

    !calculate some of the variables
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
      DO jk = 1 , nlev
        ptr_nh_prog%rho(1:nlen,jk,jb) = (ptr_nh_prog%exner(1:nlen,jk,jb)**cvd_o_rd)*p0ref/rd / &
                                         ptr_nh_prog%theta_v(1:nlen,jk,jb)     
      END DO !jk
    ENDDO !jb


  END SUBROUTINE init_warm_bubble

  
  !--------------------------------------------------
  ! read initial soil profiles and T_G from SCM netCDF file
  
  SUBROUTINE  read_soil_profile_nc(w_so_in, t_so_in, t_g_in)

    REAL(wp), INTENT(OUT), OPTIONAL :: w_so_in(nlev_soil)
    REAL(wp), INTENT(OUT), OPTIONAL :: t_so_in(nlev_soil+1)
    REAL(wp), INTENT(OUT), OPTIONAL :: t_g_in

    ! Local variables
    CHARACTER(len=max_char_length),PARAMETER :: routine  = &
          &   'mo_nh_torus_exp:read_soil_profile_nc'

    INTEGER :: levWsoil, levTsoil, nt
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)  :: tempf
    REAL(wp), ALLOCATABLE, DIMENSION(:)    :: tempg

    INTEGER :: varid
    INTEGER :: fileid     !< id number of netcdf file
    INTEGER :: dimid      !< id number of dimension
 
    !--------------------------------------------------

    CALL message(TRIM(routine), 'READING INITIAL SOIL PROFILE!')
 
    CALL nf (nf90_open('init_SCM.nc', NF90_NOWRITE, fileid), &
      & TRIM(routine)//'   File init_SCM.nc cannot be opened (soil)') 

    CALL nf (nf90_inq_dimid(fileid, 'nt', dimid), routine)
    CALL nf (nf90_inquire_dimension(fileid, dimid, len = nt), routine)

    IF (PRESENT(w_so_in) .or. PRESENT(t_so_in)) THEN  
      CALL nf (nf90_inq_dimid(fileid, 'levTsoil', dimid), routine)
      CALL nf (nf90_inquire_dimension(fileid, dimid, len = levTsoil), routine)

      CALL nf (nf90_inq_dimid(fileid, 'levWsoil', dimid), routine)
      CALL nf (nf90_inquire_dimension(fileid, dimid, len = levWsoil), routine)

!     WRITE(*,*) 'SCM read_soil_profile_nc: levTsoil ', levTsoil, ' levWsoil ', levWsoil, ' nt ', nt

    ELSE
      WRITE(*,*) 'SCM read_soil_profile_nc: ', nt
    END IF

    IF (PRESENT(t_so_in)) THEN
      IF ( levTsoil == (nlev_soil+1) ) THEN
        ALLOCATE(tempf(levTsoil,nt))
        CALL nf (nf90_inq_varid(fileid, 'T_SO', varid), routine)
        CALL nf (nf90_get_var(fileid, varid,tempf), routine)
        t_so_in   = tempf(:,1)
!       write(*,*) 't_so_in', t_so_in
        DEALLOCATE(tempf)
      ELSE
        WRITE(*,*) 'T_SO in init_SCM.nc not enough levels: levTsoil= ', levTsoil, &
                   '  nlev_soil+1= ', nlev_soil+1, '  ... cold start'
        t_so_in   = 0.0_wp
      END IF
    END IF

    IF (PRESENT(w_so_in)) THEN
      IF ( levWsoil == nlev_soil ) THEN
        ALLOCATE(tempf(levWsoil,nt))
        CALL nf (nf90_inq_varid(fileid, 'W_SO', varid), routine)
        CALL nf (nf90_get_var(fileid, varid,tempf), routine)
        w_so_in   = tempf(:,1)
!       write(*,*) 'w_so_in', w_so_in
        DEALLOCATE(tempf)
      ELSE
        WRITE(*,*) 'W_SO in init_SCM.nc not enough levels: levWsoil= ', levWsoil, &
                   '  nlev_soil= ', nlev_soil, '  ... cold start'
        w_so_in   = 0.0_wp
      END IF
    END IF

    IF (PRESENT(t_g_in)) THEN
      ALLOCATE(tempg(nt))
      CALL nf (nf90_inq_varid(fileid, 'tg', varid), routine)
      CALL nf (nf90_get_var(fileid, varid, tempg), routine)
      t_g_in    = tempg(1)
!     write(*,*) 't_g_in', t_g_in
      DEALLOCATE(tempg)
    ENDIF

    CALL nf (nf90_close(fileid), routine)

  END SUBROUTINE read_soil_profile_nc


  !--------------------------------------------------
  ! read initial soil profiles and T_G from SCM netCDF file (unified format)

  SUBROUTINE  read_soil_profile_nc_uf(w_so_in, t_so_in, t_g_in)

    REAL(wp), INTENT(OUT) :: w_so_in(nlev_soil)
    REAL(wp), INTENT(OUT) :: t_so_in(nlev_soil+1)
    REAL(wp), INTENT(OUT), OPTIONAL :: t_g_in

    ! Local variables
    CHARACTER(len=max_char_length),PARAMETER :: routine  = &
          &   'mo_nh_torus_exp:read_soil_profile_nc_uf'

    !--------------------------------------------------

    CALL message(TRIM(routine), &
    'READING INITIAL SOIL PROFILE not implemented for unified format')

  END SUBROUTINE read_soil_profile_nc_uf

  
  !--------------------------------------------------
  ! read lat lon from SCM netCDF file

  SUBROUTINE  read_latlon_scm_nc(lat_scm,lon_scm)

    REAL(wp),  INTENT(OUT) :: lat_scm
    REAL(wp),  INTENT(OUT) :: lon_scm

    ! Local variables
    CHARACTER(len=max_char_length),PARAMETER :: routine  = 'mo_nh_torus_exp:read_latlon_scm_nc'

    INTEGER :: varid
    INTEGER :: fileid     !< id number of netcdf file
    REAL(wp):: tmp_nf(1)
 
    !--------------------------------------------------
    CALL message(TRIM(routine), 'READING lat/lon FOR SCM')
 
    CALL nf (nf90_open('init_SCM.nc', NF90_NOWRITE, fileid), &
      & TRIM(routine)//'   File init_SCM.nc cannot be opened (lat/lon)') 

    CALL nf (nf90_inq_varid(fileid, 'latitude', varid) , routine)
    CALL nf (nf90_get_var(fileid, varid, tmp_nf), routine)
    lat_scm = tmp_nf(1)

    CALL nf (nf90_inq_varid(fileid, 'longitude', varid), routine)
    CALL nf (nf90_get_var(fileid, varid,tmp_nf) , routine)
    lon_scm = tmp_nf(1)

    WRITE(message_text,'(a,f7.4,a,f7.4)') 'lat_scm, ', lat_scm, ', lon_scm, ',lon_scm
    CALL message (routine,message_text)

    CALL nf (nf90_close(fileid), routine)

  END SUBROUTINE read_latlon_scm_nc


  !--------------------------------------------------
  ! read lat lon from SCM netCDF file (unified format)

  SUBROUTINE  read_latlon_scm_nc_uf(lat_scm,lon_scm) !unified format

    REAL(wp),  INTENT(OUT) :: lat_scm
    REAL(wp),  INTENT(OUT) :: lon_scm
    REAL(wp), ALLOCATABLE, DIMENSION(:)  :: lats, lons

    ! Local variables
    CHARACTER(len=max_char_length),PARAMETER :: routine  = &
          &   'mo_nh_torus_exp:read_latlon_scm_nc_uf'

    INTEGER :: varid
    INTEGER :: fileid     !< id number of netcdf file
    INTEGER :: dimid      !< id number of dimension
    INTEGER :: lat,lon
 
    !--------------------------------------------------

    CALL message(TRIM(routine), 'READING lat/lon FOR SCM')
 
    CALL nf (nf90_open('init_SCM.nc', NF90_NOWRITE, fileid), &
      & TRIM(routine)//'   File init_SCM.nc cannot be opened (lat/lon)') 

    CALL nf (nf90_inq_dimid(fileid, 'lat', dimid), routine)
    CALL nf (nf90_inquire_dimension(fileid, dimid, len = lat), routine)

    CALL nf (nf90_inq_dimid(fileid, 'lon', dimid), routine)
    CALL nf (nf90_inquire_dimension(fileid, dimid, len = lon), routine)

    ALLOCATE(lats(lat),lons(lon))

    CALL nf (nf90_inq_varid(fileid, 'lat', varid), routine)
    CALL nf (nf90_get_var(fileid, varid,lats), routine)
    lat_scm=lats(1)

    CALL nf (nf90_inq_varid(fileid, 'lon', varid), routine)
    CALL nf (nf90_get_var(fileid, varid,lons), routine)
    lon_scm=lons(1)

    WRITE(message_text,'(a,f7.4,a,f7.4)') 'lat_scm, ', lat_scm, ', lon_scm, ',lon_scm
    CALL message (routine,message_text)

    CALL nf (nf90_close(fileid), routine)
    DEALLOCATE(lats,lons)

  END SUBROUTINE read_latlon_scm_nc_uf

 
  !--------------------------------------------------
  ! read external parameters from SCM netCDF file

  SUBROUTINE read_ext_scm_nc (num_lcc,soiltyp_scm,fr_land_scm,plcov_mx_scm,lai_mx_scm, & 
                              rootdp_scm,rsmin_scm,z0_scm,topo_scm,emis_rad_scm,&
                              lu_class_fr,lctype_scm)
     
    INTEGER , INTENT(IN)  :: num_lcc          ! number of landcover classes
    INTEGER , INTENT(OUT) :: soiltyp_scm      ! soil type
    REAL(wp), INTENT(OUT) :: fr_land_scm      ! land fraction
    REAL(wp), INTENT(OUT) :: plcov_mx_scm     ! maximum plant cover
    REAL(wp), INTENT(OUT) :: lai_mx_scm       ! maximum leaf area index
    REAL(wp), INTENT(OUT) :: rootdp_scm       ! root depth
    REAL(wp), INTENT(OUT) :: rsmin_scm        ! minimum stomata resistance
    REAL(wp), INTENT(OUT) :: z0_scm           ! roughness length
    REAL(wp), INTENT(OUT) :: topo_scm         ! height above sea level
    REAL(wp), INTENT(OUT) :: emis_rad_scm     ! emisivity
    REAL(wp), INTENT(OUT) :: lu_class_fr(:)   ! land use classes fractions
    CHARACTER(len=max_char_length),INTENT(OUT) ::lctype_scm !data source for land use
    INTEGER :: nCLU

    INTEGER :: varid
    INTEGER :: fileid     !< id number of netcdf file
    INTEGER :: dimid      !< id number of dimension
                           
    ! Local variables
    CHARACTER(len=max_char_length),PARAMETER :: routine  = 'mo_nh_torus_exp:read_ext_SCM_nc'
    REAL(wp) :: tmp_nf(1)
 
    !-------------------------------------------------

    CALL message(TRIM(routine), 'READING EXTERNAL DATA FOR SCM')
    lctype_scm=""

    CALL nf (nf90_open('init_SCM.nc', NF90_NOWRITE, fileid), &
      & TRIM(routine)//'   File init_SCM.nc cannot be opened (external)') 

    CALL nf (nf90_inq_dimid(fileid, 'nclass_lu', dimid), routine)
    CALL nf (nf90_inquire_dimension(fileid, dimid, len = nCLU), routine)

    IF (nCLU .NE. num_lcc) THEN
      CALL finish( 'testcases/mo_nh_torus_exp.f90',&
      'Number of LU classes in init_SCM.nc does not match num_lcc')
    ENDIF

    CALL nf (nf90_inq_varid(fileid, 'FR_LAND',   varid), routine)
    CALL nf (nf90_get_var(fileid, varid, tmp_nf), routine)
     fr_land_scm = tmp_nf(1)
    CALL nf (nf90_inq_varid(fileid, 'FR_LAND',   varid), routine)
    CALL nf (nf90_get_var(fileid, varid, tmp_nf), routine)
     fr_land_scm = tmp_nf(1)
    CALL nf (nf90_inq_varid(fileid, 'PLCOV_MX',  varid), routine)
    CALL nf (nf90_get_var(fileid, varid, tmp_nf), routine)
     plcov_mx_scm = tmp_nf(1)
    CALL nf (nf90_inq_varid(fileid, 'LAI_MX',    varid), routine)
    CALL nf (nf90_get_var(fileid, varid, tmp_nf), routine)
     lai_mx_scm = tmp_nf(1)
    CALL nf (nf90_inq_varid(fileid, 'ROOTDP',    varid), routine)
    CALL nf (nf90_get_var(fileid, varid, tmp_nf), routine)
     rootdp_scm = tmp_nf(1)
    CALL nf (nf90_inq_varid(fileid, 'RSMIN',     varid), routine)
    CALL nf (nf90_get_var(fileid, varid, tmp_nf), routine)
     rsmin_scm = tmp_nf(1)
    CALL nf (nf90_inq_varid(fileid, 'SOILTYP',   varid), routine)
    CALL nf (nf90_get_var(fileid, varid, tmp_nf), routine)
     soiltyp_scm=INT(tmp_nf(1))
    CALL nf (nf90_inq_varid(fileid, 'Z0',        varid), routine)
    CALL nf (nf90_get_var(fileid, varid, tmp_nf), routine)
     z0_scm = tmp_nf(1)
    CALL nf (nf90_inq_varid(fileid, 'TOPO',      varid), routine)
    CALL nf (nf90_get_var(fileid, varid, tmp_nf), routine)
     topo_scm = tmp_nf(1)
    CALL nf (nf90_inq_varid(fileid, 'EMIS_RAD',  varid), routine)
    CALL nf (nf90_get_var(fileid, varid, tmp_nf), routine)
     emis_rad_scm = tmp_nf(1)
    CALL nf (nf90_inq_varid(fileid, 'LU_CLASS_FRACTION', varid), routine)
    CALL nf (nf90_get_var(fileid, varid,lu_class_fr), routine)
    CALL nf (nf90_get_att(fileid, varid, 'lctype', lctype_scm), routine)
  
    IF ( get_my_global_mpi_id() == 0 ) THEN
      print *,TRIM(routine),'  printing external surface parameters for SCM'
      print *,'  fr_land_scm =',   fr_land_scm
      print *,'  plcov_mx_scm=',   plcov_mx_scm
      print *,'  lai_mx_scm  =',   lai_mx_scm
      print *,'  rootdp_scm  =',   rootdp_scm
      print *,'  rsmin_scm   =',   rsmin_scm
      print *,'  soiltyp_scm =',   soiltyp_scm
      print *,'  z0_scm      =',   z0_scm
      print *,'  topo_scm    =',   topo_scm
      print *,'  emis_rad_scm=',   emis_rad_scm
      print *,'  lu_class_fr =',   lu_class_fr
      print *,'  nCLU        =',   nCLU
    END IF
  
    CALL nf (nf90_close(fileid), routine)

  END SUBROUTINE read_ext_scm_nc


  !--------------------------------------------------
  ! read external parameters from SCM netCDF file (unified format)

  SUBROUTINE read_ext_scm_nc_uf (num_lcc,soiltyp_scm,fr_land_scm,plcov_mx_scm,lai_mx_scm, & 
                              rootdp_scm,rsmin_scm,z0_scm,topo_scm,emis_rad_scm,&
                              lu_class_fr,lctype_scm)
     
    INTEGER , INTENT(IN)  :: num_lcc          ! number of landcover classes
    INTEGER , INTENT(OUT) :: soiltyp_scm      ! soil type
    REAL(wp), INTENT(OUT) :: fr_land_scm      ! land fraction
    REAL(wp), INTENT(OUT) :: plcov_mx_scm     ! maximum plant cover
    REAL(wp), INTENT(OUT) :: lai_mx_scm       ! maximum leaf area index
    REAL(wp), INTENT(OUT) :: rootdp_scm       ! root depth
    REAL(wp), INTENT(OUT) :: rsmin_scm        ! minimum stomata resistance
    REAL(wp), INTENT(OUT) :: z0_scm           ! roughness length
    REAL(wp), INTENT(OUT) :: topo_scm         ! height above sea level
    REAL(wp), INTENT(OUT) :: emis_rad_scm     ! emisivity
    REAL(wp), INTENT(OUT) :: lu_class_fr(:)   ! land use classes fractions
    CHARACTER(len=max_char_length),INTENT(OUT) ::lctype_scm !data source for land use

    INTEGER :: nCLU

    INTEGER :: varid      ! id number of variable (or attribute variable)
    INTEGER :: attid      ! id number of attribute associated to variable (not useful)
    INTEGER :: fileid     ! id number of netcdf file
                           
    ! Local variables
    CHARACTER(len=max_char_length),PARAMETER :: routine  = 'mo_nh_torus_exp:read_ext_SCM_nc_uf'
 
    !------------------------------------------------

    CALL message(TRIM(routine), &
      'READING EXTERNAL DATA FOR SCM not implemented for unified format')

    lctype_scm= "GLOBCOVER2009"
!    topo_scm  = 314._wp
!    z0_scm    = 0.035_wp


    CALL nf (nf90_open('init_SCM.nc', NF90_NOWRITE, fileid), &
      & TRIM(routine)//'   File init_SCM.nc cannot be opened (external)') 

    CALL nf (nf90_inquire_attribute(fileid, varid, 'z0', attnum = attid), routine)
    CALL nf (nf90_get_att(fileid, varid, 'z0', z0_scm), routine)

    CALL nf (nf90_inquire_attribute(fileid, varid, 'zorog', attnum = attid), routine)
    CALL nf (nf90_get_att(fileid, varid, 'zorog', topo_scm), routine)


    IF ( get_my_global_mpi_id() == 0 ) THEN
      print *,TRIM(routine),'  printing external surface parameters for SCM'
      print *,'  fr_land_scm =',   fr_land_scm
      print *,'  plcov_mx_scm=',   plcov_mx_scm
      print *,'  lai_mx_scm  =',   lai_mx_scm
      print *,'  rootdp_scm  =',   rootdp_scm
      print *,'  rsmin_scm   =',   rsmin_scm
      print *,'  soiltyp_scm =',   soiltyp_scm
      print *,'  z0_scm      =',   z0_scm
      print *,'  topo_scm    =',   topo_scm
      print *,'  emis_rad_scm=',   emis_rad_scm
      print *,'  lu_class_fr =',   lu_class_fr
      print *,'  nCLU        =',   nCLU
    END IF
  
  END SUBROUTINE read_ext_scm_nc_uf

  
  !-----------------------------------------------------
  !set boundary conditions for SCM
  
  SUBROUTINE  set_scm_bnd( nvec, ivstart, ivend, u_s, v_s, th_b, qv_b, pres_sfc, dz_bs,z0m,z0h,&
    & prm_nwp_tend, tvm, tvh, shfl_s, qhfl_s, lhfl_s,umfl_s,vmfl_s, qv_s, t_g )

    INTEGER,        INTENT(IN) :: &
    nvec,         & ! nproma
    ivstart,      & ! start index in the nproma vector
    ivend           ! end index in the nproma vector

    REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(IN) :: &
    u_s,          & ! zonal wind speed in surface layer       (at mass positions)           ( m/s )
    v_s,          & ! meridional wind speed  in surface layer (at mass positions)           ( m/s )
    th_b,         & ! pot. temperature at the lowest model level                            ( K )
    qv_b,         & ! total specific moisture at the lowest model level (at mass positions) ( kg/kg )
    pres_sfc,     & ! pressure at surface [Pa]
    dz_bs,        & ! distance between surface and lowest model level [m]
    z0m,          & ! roughness length for momentum [m]
    z0h             ! roughness length for heat/moisture [m]

    TYPE(t_nwp_phy_tend), TARGET, INTENT(IN):: &
    prm_nwp_tend    ! atm tend vars

    ! turbulent (transfer) velocity scales at the surface
    REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(INOUT), OPTIONAL :: &
    tvm,          & ! for momentum                                  ( m/s )
    tvh             ! for heat and moisture                         ( m/s )

    REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(OUT), OPTIONAL :: &
    shfl_s,       & ! sensible heat flux at the surface             (W/m2)    (positive downward)
    qhfl_s,       & ! water vapor   flux at the surface             (kg/m2/s) (positive downward)
    lhfl_s,       & ! latent heat   flux at the surface             (positive downward)
    umfl_s,       & ! turbulent u-flux at the surface               (positive downward)
    vmfl_s          ! turbulent v-flux at the surface               (positive downward)

    REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(OUT), OPTIONAL :: &
    qv_s,         & ! specific water vapor content on the surface   (kg/kg)
    t_g             ! weighted surface temperature                  (  k  )

    CHARACTER(LEN=max_char_length), PARAMETER :: &
    routine = 'mo_nh_torus_exp:set_scm_bnd'

    INTEGER  :: i
    REAL(KIND=wp) :: velo(ivstart:ivend)
    REAL(KIND=wp) :: rho_sfc(ivstart:ivend),tempv_sfc(ivstart:ivend)
    REAL(KIND=wp) :: rib(ivstart:ivend) !bulk Richardson number - for Louis scheme
    REAL(KIND=wp) :: ribm(ivstart:ivend)!bulk Richardson number with asymptotic correction momentum
    REAL(KIND=wp) :: ribh(ivstart:ivend)!bulk Richardson number with asymptotic correction heat/moist.
    REAL(KIND=wp) :: fh(ivstart:ivend)  !stability dep. function for heat and moisture - for Louis scheme
    REAL(KIND=wp) :: fm(ivstart:ivend)  !stability dep. function for momentum - for Louis scheme
    REAL(KIND=wp) :: cnm(ivstart:ivend) !drag coefficient for momentum at neutrality
    REAL(KIND=wp) :: cnh(ivstart:ivend) !drag coefficient for heat/moist. at neutrality
    REAL(KIND=wp) :: b_louis,cm_louis,ch_louis,d_louis,ricr

    b_louis=5.0_wp
    cm_louis=5.0_wp
    ch_louis=5.0_wp
    d_louis=5.0_wp
    ricr=1.0_wp

    ! if apply_ls_forcing has not been called so far the
    ! variables in prm_nwp_tend% are not defined
    IF (.not. lscm_ls_forcing_ini) RETURN

    DO i=ivstart, ivend
      velo(i)=MAX( vel_min, SQRT(u_s(i)**2+v_s(i)**2) )
      tempv_sfc(i) = t_g(i) * (1._wp + vtmpc1*qv_b(i))
      rho_sfc(i)   = pres_sfc(i)/(rd*tempv_sfc(i))
      !write(*,*) "rhos:",i,velo(i),tempv_sfc(i),rho_sfc(i),u_s(i),v_s(i),&
      !&t_g(i),qv_b(i),pres_sfc(i)
    ENDDO

    !GABLS1
    !dry, should be improved to virtual temperature for other cases
    !also, th_surf/[th_surf,qv_surf] must be prescribed
    IF((scm_sfc_temp.eq.5)) then
      DO i=ivstart, ivend
        cnm(i) =(0.4_wp/log(1.0_wp+dz_bs(i)/z0m(i)))**2 
        cnh(i) =0.4_wp**2/(log(1.0_wp+dz_bs(i)/z0m(i))*log(1.0_wp+dz_bs(i)/z0h(i)))
        rib(i) =grav*(th_b(i)-prm_nwp_tend%fc_ts)/th_b(i)*dz_bs(i)/velo(i)**2 
        ribm(i) = rib(i)/(1.0_wp+rib(i)/ricr)
        ribh(i) = rib(i)/(1.0_wp+3.0_wp*rib(i)/ricr)**(0.333_wp)
        IF (rib(i).ge.0.0_wp) THEN
          !stable stratification
          fm(i) =(1.0_wp/(1.0_wp+2.0_wp*b_louis*ribm(i)/sqrt(1.0_wp+d_louis*ribm(i))))
          fh(i) =(1.0_wp/(1.0_wp+3.0_wp*b_louis*ribh(i)*sqrt(1.0_wp+d_louis*ribh(i))))
        ELSE
          !unstable stratification
          fm(i) =(1.0_wp-2.0_wp*b_louis*rib(i)/&
          & (1.0_wp+3.0_wp*b_louis*cm_louis*cnm(i)*&
          & sqrt(-rib(i)*(1.0_wp+dz_bs(i)/z0m(i)))))
          fh(i) =(1.0_wp-3.0_wp*b_louis*rib(i)/&
          & (1.0_wp+3.0_wp*b_louis*ch_louis*cnh(i)*&
          & sqrt(-rib(i)*(1.0_wp+dz_bs(i)/z0h(i)))))
        END IF
      ENDDO
    ENDIF

! surface temperature and sensible heat flux 

    SELECT CASE(scm_sfc_temp)
    CASE (0) ! no prescribed t_g and shfl_s
       !t_g(i)    = prm_nwp_tend%fc_tg !for radiation and transfer scheme
    CASE (1)
      DO i=ivstart, ivend
        t_g(i)    = prm_nwp_tend%fc_tg !for radiation and transfer scheme
      ENDDO
    CASE (2) ! prescribed fluxes
      DO i=ivstart, ivend
        t_g(i)    =   prm_nwp_tend%fc_tg !for radiation and transfer scheme
        shfl_s(i) = - prm_nwp_tend%fc_sfc_sens_flx
      ENDDO
    CASE (4) ! prescribed drag coefficient
      DO i=ivstart, ivend
        t_g(i)    = prm_nwp_tend%fc_tg !for radiation and transfer scheme
       !tvh(i)    = prm_nwp_tend%fc_Ch*velo(i)
        shfl_s(i) = rho_sfc(i)*prm_nwp_tend%fc_Ch*(th_b(i)-prm_nwp_tend%fc_ts)*velo(i)*cpd
      ENDDO
    CASE(5) ! Louis scheme - GABLS1
      DO i=ivstart, ivend
        t_g(i)    = prm_nwp_tend%fc_tg !for radiation and transfer scheme
        shfl_s(i) = rho_sfc(i)*cnh(i)*fh(i)*(th_b(i)-prm_nwp_tend%fc_ts)*velo(i)*cpd
      ENDDO
    CASE DEFAULT
      call finish(routine,' Value for scm_sfc_temp not known!')
    END SELECT

! surface moisture and latent heat flux

    SELECT CASE(scm_sfc_qv)
    CASE (0) ! no prescribed qv_s, lhfs_s and qhfl_s
    CASE (1) ! prescribed qv
      DO i=ivstart, ivend
        qv_s(i)   = prm_nwp_tend%fc_qvs
      ENDDO
    CASE (2) ! prescribed flux
      DO i=ivstart, ivend
        qv_s(i)   =   prm_nwp_tend%fc_qvs
        lhfl_s(i) = - prm_nwp_tend%fc_sfc_lat_flx
        qhfl_s(i) = - prm_nwp_tend%fc_sfc_lat_flx/lh_v
      ENDDO
    CASE (3) ! qv_s based on saturation
      DO i=ivstart, ivend
        qv_s(i)   = spec_humi( sat_pres_water(t_g(i)) , pres_sfc(i) )   
      END DO
    CASE (4) ! prescribed drag coefficient
      DO i=ivstart, ivend
        !qv_s(i)  = prm_nwp_tend%fc_qvs
        !!tvh(i)  = prm_nwp_tend%fc_Cq*velo(i) !no tvq, hence Ch is used for both moisture and temperature
        !tvh(i)   = prm_nwp_tend%fc_Ch*velo(i) !no tvq, hence Ch is used for both moisture and temperature
        qhfl_s(i) = rho_sfc(i)*prm_nwp_tend%fc_Ch*(qv_b(i)-prm_nwp_tend%fc_qvs)*velo(i)
        lhfl_s(i) = qhfl_s(i)*lh_v
      ENDDO
    CASE(5) ! Louis scheme - GABLS1
      DO i=ivstart, ivend
        qhfl_s(i) = rho_sfc(i)*cnh(i)*fh(i)*(qv_b(i)-prm_nwp_tend%fc_qvs)*velo(i)
        lhfl_s(i) = qhfl_s(i)*lh_v
      ENDDO
    CASE DEFAULT
      CALL finish(routine,' Value for scm_sfc_qv not known!')
    END SELECT

! surface momentum flux
    
    SELECT CASE (scm_sfc_mom)
    CASE (0) ! no prescribe flux
    CASE (2) ! prescribed flux
      DO i=ivstart, ivend
        tvm(i)    = prm_nwp_tend%fc_ustar **2 / velo(i)
        umfl_s(i) = -rho_sfc(i)*tvm(i)*u_s(i)
        vmfl_s(i) = -rho_sfc(i)*tvm(i)*v_s(i)
      ENDDO
    CASE (4) ! prescribed drag coefficient
      DO i=ivstart, ivend
        tvm(i)    = prm_nwp_tend%fc_Cm * velo(i)
        umfl_s(i) = -rho_sfc(i)*tvm(i)*u_s(i)
        vmfl_s(i) = -rho_sfc(i)*tvm(i)*v_s(i)
     ENDDO
    CASE (5) ! Louis scheme - GABLS1
      DO i=ivstart, ivend
        tvm(i)    = cnm(i)*fm(i)*velo(i)
        umfl_s(i) = -rho_sfc(i)*tvm(i)*u_s(i)
        vmfl_s(i) = -rho_sfc(i)*tvm(i)*v_s(i)
      ENDDO
    CASE DEFAULT
       CALL finish(routine,' Value for scm_sfc_mom not known!') 
    END SELECT


   !IF ( get_my_global_mpi_id() == 0 ) THEN
   !  WRITE(*,*) TRIM(routine),'  printing surface boundary parameters for SCM'
   !  WRITE(*,*) "   tvm",    tvm(1)
   !  WRITE(*,*) "   tvh",    tvh(1)
   !  WRITE(*,*) "shfl_s", shfl_s(1)
   !  WRITE(*,*) "qhfl_s", qhfl_s(1)
   !  WRITE(*,*) "lhfl_s", lhfl_s(1)
   !  WRITE(*,*) "umfl_s", umfl_s(1)
   !  WRITE(*,*) "vmfl_s", vmfl_s(1)
   !  WRITE(*,*) "  qv_s",   qv_s(1)
   !  WRITE(*,*) "   t_g",    t_g(1)
   !END IF

  END SUBROUTINE set_scm_bnd
  

END MODULE mo_nh_torus_exp
