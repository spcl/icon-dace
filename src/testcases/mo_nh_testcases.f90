! Defines the artificial testcases for the nonhydrostatic atmospheric model.
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

MODULE mo_nh_testcases  
!-------------------------------------------------------------------------  
!  
!    ProTeX FORTRAN source: Style 2  
!    modified for ICON project, DWD/MPI-M 2006                       
!  
!-------------------------------------------------------------------------  
!  
!  
!  
  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: message, finish, message_text

  USE mo_nh_testcases_nml

  USE mo_impl_constants,       ONLY: MAX_CHAR_LENGTH, inwp, icosmo
  USE mo_grid_config,          ONLY: lplane, n_dom, l_limited_area
  USE mo_model_domain,         ONLY: t_patch
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_math_constants,       ONLY: pi
  USE mo_math_types,           ONLY: t_cartesian_coordinates, t_geographical_coordinates 
  USE mo_math_utilities,       ONLY: gc2cc, arc_length
  USE mo_parallel_config,      ONLY: nproma
  USE mo_run_config,           ONLY: ltransport, iforcing
  USE mo_extpar_config,        ONLY: itopo
  USE mo_dynamics_config,      ONLY: nnow, nnew, lcoriolis
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_physical_constants,   ONLY: grav, cpd, rd, cvd_o_rd, p0ref
  USE mo_nonhydro_types,       ONLY: t_nh_state
  USE mo_nonhydro_state,       ONLY: duplicate_prog_state
  USE mo_intp_data_strc,       ONLY: t_int_state
  USE mo_nwp_lnd_types,        ONLY: t_lnd_state
  USE mo_lnd_nwp_config,       ONLY: isub_seaice
  USE mo_nh_pa_test,           ONLY: init_nh_state_prog_patest
  USE mo_nh_df_test,           ONLY: init_nh_state_prog_dftest
  USE mo_nh_hs_test,           ONLY: init_nh_state_prog_held_suarez
  USE mo_nh_jabw_exp,          ONLY: init_nh_topo_jabw, init_nh_state_prog_jabw,  & 
                                   & init_passive_tracers_nh_jabw, init_nh_inwp_tracers
  USE mo_nh_mrw_exp,           ONLY: init_nh_topo_mrw, init_nh_state_prog_mrw,    &
                                   & init_nh_prog_mwbr_const, mount_half_width                     
  USE mo_nh_wk_exp,            ONLY: init_nh_topo_wk, init_nh_env_wk,             &
                                   & init_nh_buble_wk, bubctr_z,                  &
                                   & bub_hor_width, bub_ver_width, bub_amp
  USE mo_nh_bb13_exp,          ONLY: init_nh_env_bb13, init_nh_bubble_bb13                       
  USE mo_nh_dcmip_gw,          ONLY: init_nh_dcmip_gw, init_nh_gw_analyt
  USE mo_nh_dcmip_hadley,      ONLY: init_nh_dcmip_hadley         
  USE mo_nh_dcmip_schaer,      ONLY: init_nh_prog_dcmip_schaer
  USE mo_nh_dcmip_rest_atm,    ONLY: init_nh_topo_dcmip_rest_atm,                 &
                                   & init_nh_prog_dcmip_rest_atm  
  USE mo_nh_dcmip_tc,          ONLY: init_nh_dcmip_tc
  USE mo_nh_dcmip_bw,          ONLY: init_nh_dcmip_bw
  USE mo_nh_dcmip_terminator,  ONLY: init_nh_dcmip_terminator
  USE mo_nh_lim_area_testcases,ONLY: init_nh_atmo_ana_nconstlayers,               &
                                   & init_nh_anaprof_uv, init_nh_topo_ana,        &
                                   & itype_atmo_ana, init_nh_atmo_ana_poly
  USE mo_nh_prog_util,         ONLY: nh_prog_add_random
  USE mo_random_util,          ONLY: add_random_noise_global
  USE mo_grid_geometry_info,   ONLY: planar_torus_geometry
  USE mo_nh_rce_exp,           ONLY: init_nh_state_rce_glb,                       &
                                   & init_nh_state_rce_tprescr_glb
  USE mo_aes_bubble,           ONLY: init_aes_bubble
  USE mo_nh_torus_exp,         ONLY: init_nh_state_cbl, init_nh_state_rico,       &
                                   & init_torus_netcdf_sounding,                  &
                                   & init_torus_ascii_sounding, init_warm_bubble, &
                                   & init_torus_rcemip_analytical_sounding
  USE mo_nh_tpe_exp,           ONLY: init_nh_state_prog_TPE
  USE mo_vertical_coord_table, ONLY: vct_a
  USE mo_hydro_adjust,         ONLY: hydro_adjust_const_thetav
  USE mo_scm_nml,              ONLY: i_scm_netcdf, lscm_random_noise

  IMPLICIT NONE  
  
  PRIVATE
  
  PUBLIC :: init_nh_testtopo
  PUBLIC :: init_nh_testcase
  PUBLIC :: init_nh_testcase_scm

  ! !DEFINED PARAMETERS for jablonowski williamson: 
  !  The rest of the needed parameters are define in mo_nh_jabw_exp
  REAL(wp), PARAMETER :: ps0      = 1.e5_wp     !< surface pressure (Pa)
  
  ! !DEFINED PARAMETERS for mountain induced Rossby wave train:
  REAL(wp), PARAMETER :: pres_sp  = 93000.0_wp  !< pressure surface at the south pole
  REAL(wp), PARAMETER :: temp_mrw = 288._wp     !< temperature of isothermal atmosphere

  ! !DEFINED PARAMETERS for APE (now read from namelist)
  !REAL(wp), PARAMETER :: zp_ape   = 101325._wp  !< surface pressure
  !REAL(wp), PARAMETER :: ztmc_ape = 25.006_wp   !< total moisture content 

  CONTAINS  
  
!-------------------------------------------------------------------------
!
!
  !>
  !! Initialize topography for nonhydrostatic artificial testcases
  !! (not for single column model (SCM))
  !! 
  SUBROUTINE init_nh_testtopo (p_patch, ext_data)
!
! !INPUT VARIABLES:
  TYPE(t_patch),TARGET,  INTENT(INOUT) :: p_patch(n_dom)
  TYPE(t_external_data), INTENT(INOUT) :: ext_data(n_dom)

  INTEGER        :: jg, jc, jb, nlen
  INTEGER        :: nblks_c, npromz_c
  REAL(wp)       :: z_lon, z_lat, z_dist
  TYPE(t_geographical_coordinates) :: z_x2_geo
  TYPE(t_cartesian_coordinates)    :: z_x1_cart, z_x2_cart
  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    &  routine = '(mo_nh_testcases) init_nh_testtopo:'
 
  LOGICAL       :: l_modified

!-----------------------------------------------------------------------

  ! Initialize topography to zero if idealized topo is used
  IF ( itopo == 0 ) THEN
    DO jg = 1, n_dom
      ext_data(jg)%atm%topography_c(1:nproma,1:p_patch(jg)%nblks_c) = 0.0_wp
    ENDDO
  ENDIF

  SELECT CASE (nh_test_name)

  CASE ('zero', 'HS_jw')
    
    IF(nh_test_name=='HS_jw') CALL message(TRIM(routine),'running the Held-Suarez test')

  CASE ('schaer')
 
    ! limited area test case by Schaer et al. (2002) (plane geometry)

    !IF(.NOT.lplane) CALL finish(TRIM(routine),'Schaer test case only for lplane=True')

    ! At present the mountain is at position lat=0,lon=0 (given in meters)
    z_x2_geo%lon = 0.0_wp
    z_x2_geo%lat = 0.0_wp

    !MB: CAUTION: something is wrong here: if the loop is vectorized (without the write-statement below), orography is set =const. !!??
!DIR$ NOVECTOR
    DO jg = 1, n_dom
      DO jb = 1, p_patch(jg)%nblks_c
        IF (jb /=  p_patch(jg)%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen =  p_patch(jg)%npromz_c
        ENDIF
        DO jc = 1, nlen
          IF ( p_patch(jg)%geometry_info%geometry_type == planar_torus_geometry ) THEN
            z_lon = p_patch(jg)%cells%cartesian_center(jc,jb)%x(1)
          ELSE
            z_lon = p_patch(jg)%cells%center(jc,jb)%lon*torus_domain_length/pi*0.5_wp
          END IF
          z_dist = z_lon - z_x2_geo%lon
          ext_data(jg)%atm%topography_c(jc,jb) = mount_height            &
            &    * EXP(-(z_dist/mount_width)**2) * ((COS(pi*z_dist/mount_width_2))**2)

          ! WRITE(*,'(A,2I7,3F15.4)' ) "topo: ", jc, jb, z_lon, z_dist, ext_data(jg)%atm%topography_c(jc,jb)

        ENDDO
      ENDDO 
    ENDDO 

  CASE ('atm_at_rest')
 
    !IF(.NOT.lplane) CALL finish(TRIM(routine),'Atm. at Rest test case only for lplane=True')

    ! At present the mountain is at position lat=0,lon=0 (given in meters)
    z_x2_geo%lon = 0.0_wp
    z_x2_geo%lat = 0.0_wp

    !MB: CAUTION: something is wrong here: if the loop is vectorized (without the write-statement below), orography is set =const. !!??
    DO jg = 1, n_dom
      DO jb = 1, p_patch(jg)%nblks_c
        IF (jb /=  p_patch(jg)%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen =  p_patch(jg)%npromz_c
        ENDIF
        DO jc = 1, nlen
          IF ( p_patch(jg)%geometry_info%geometry_type == planar_torus_geometry ) THEN
            z_lon = p_patch(jg)%cells%cartesian_center(jc,jb)%x(1)
          ELSE
            z_lon = p_patch(jg)%cells%center(jc,jb)%lon*torus_domain_length/pi*0.5_wp
          END IF
          z_dist = z_lon-z_x2_geo%lon
          ext_data(jg)%atm%topography_c(jc,jb) = mount_height * EXP(-(z_dist/mount_width)**2)

          !WRITE(*,'(A,2I7,3F15.4)' ) "topo: ", jc, jb, z_lon, z_dist, ext_data(jg)%atm%topography_c(jc,jb)

        ENDDO
      ENDDO 
    ENDDO 

  CASE ('gauss3D')
 
    !IF(.NOT.lplane) CALL finish(TRIM(routine),'Atm. at Rest test case only for lplane=True')

    ! At present the mountain is at position lat=0,lon=0 (given in meters)
    z_x2_geo%lon = 0.0_wp
    z_x2_geo%lat = 0.0_wp

    !MB: CAUTION: something is wrong here: if the loop is vectorized (without the write-statement below), orography is set =const. !!??
    DO jg = 1, n_dom
      DO jb = 1, p_patch(jg)%nblks_c
        IF (jb /=  p_patch(jg)%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen =  p_patch(jg)%npromz_c
        ENDIF
        DO jc = 1, nlen
          IF ( p_patch(jg)%geometry_info%geometry_type == planar_torus_geometry ) THEN
            z_lon = p_patch(jg)%cells%cartesian_center(jc,jb)%x(1)
            z_lat = p_patch(jg)%cells%cartesian_center(jc,jb)%x(2)
          ELSE
            z_lon = p_patch(jg)%cells%center(jc,jb)%lon * torus_domain_length/pi*0.5_wp
            z_lat = p_patch(jg)%cells%center(jc,jb)%lat * torus_domain_length/pi*0.5_wp
          END IF
          z_dist = SQRT( (z_lon-z_x2_geo%lon)**2 + (z_lat-z_x2_geo%lat)**2 )
          ext_data(jg)%atm%topography_c(jc,jb) = mount_height * EXP(-(z_dist/mount_width)**2)
          !ext_data(jg)%atm%topography_c(jc,jb) = mount_height * 2.0_wp**( -(z_dist/mount_width)**2 )

          !WRITE(*,'(A,2I7,3F15.4)' ) "topo: ", jc, jb, z_lon, z_dist, ext_data(jg)%atm%topography_c(jc,jb)

        ENDDO
      ENDDO 
    ENDDO 

  CASE ('bell')

    ! At present the mountain is at position lat=0,lon=0 (given in meters)
    z_x2_geo%lon = 0.0_wp
    z_x2_geo%lat = 0.0_wp

    DO jg = 1, n_dom
      DO jb = 1, p_patch(jg)%nblks_c
        IF (jb /=  p_patch(jg)%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen =  p_patch(jg)%npromz_c
        ENDIF
        DO jc = 1, nlen
          IF (lplane) THEN
            z_lat  = 0.0_wp!p_patch(jg)%cells%center(jc,jb)%lat*torus_domain_length/pi*0.5_wp
            z_lon  = p_patch(jg)%cells%center(jc,jb)%lon*torus_domain_length/pi*0.5_wp
            z_dist = SQRT((z_lat-z_x2_geo%lat)**2+(z_lon-z_x2_geo%lon)**2)
          ELSE
            z_x1_cart    = gc2cc(p_patch(jg)%cells%center(jc,jb))
            z_x2_cart    = gc2cc(z_x2_geo)
            z_dist       = arc_length(z_x1_cart,z_x2_cart)
          ENDIF
          ext_data(jg)%atm%topography_c(jc,jb) = mount_height/ &
                 (1.0_wp+ (z_dist/mount_half_width)**2)**1.5_wp
        ENDDO
      ENDDO 
    ENDDO 

  CASE ('jabw', 'jabw_s')

    DO jg = 1, n_dom 
     nblks_c   = p_patch(jg)%nblks_c
     npromz_c  = p_patch(jg)%npromz_c

     CALL init_nh_topo_jabw ( p_patch(jg),ext_data(jg)%atm%topography_c, nblks_c, npromz_c, jw_u0)
    END DO

  CASE ('jabw_m')  

    DO jg = 1, n_dom 
     nblks_c   = p_patch(jg)%nblks_c
     npromz_c  = p_patch(jg)%npromz_c

     CALL init_nh_topo_jabw ( p_patch(jg),ext_data(jg)%atm%topography_c, nblks_c, npromz_c, jw_u0, &
                          & opt_m_height = mount_height, opt_m_half_width = mount_half_width )
    END DO

  CASE ('dcmip_bw_11')
    ! itopo == 0 --> The topography is initialized to 0 at the beginning of this subroutine
    CALL message(TRIM(routine),'running DCMIP2016 baroclinic wave test case 1-1')

  CASE ('mrw_nh', 'mrw2_nh' , 'mwbr_const')

   CALL message(TRIM(routine),'running mrw, setting topography')

   l_modified = .FALSE.
   IF (nh_test_name=='mrw2_nh'  ) THEN
     l_modified = .TRUE.
   ENDIF

   DO jg = 1, n_dom 
     nblks_c   = p_patch(jg)%nblks_c
     npromz_c  = p_patch(jg)%npromz_c

     CALL init_nh_topo_mrw ( p_patch(jg),ext_data(jg)%atm%topography_c, nblks_c, npromz_c, l_modified) 
   ENDDO

   CALL message(TRIM(routine),'topography is initialised ')

  CASE ('wk82')  

    DO jg = 1, n_dom 
     nblks_c   = p_patch(jg)%nblks_c
     npromz_c  = p_patch(jg)%npromz_c

     CALL init_nh_topo_wk ( ext_data(jg)%atm%topography_c, nblks_c, npromz_c)
    END DO

  CASE ('bb13')
    ! Test case Baldauf, Brdar (2013) QJRMS (linear gravity/sound wave expansion in a channel)
    CALL message(TRIM(routine), "no orography for testcase bb13")

  CASE ('straka93')
    ! Test case Straka et al. (1993) (falling cold bubble)
    CALL message(TRIM(routine), "no orography for testcase straka93")

  CASE ('PA')
   ! The topography ist initialized in "init_nh_state_prog_patest"
    CALL message(TRIM(routine),'running the Pure 3D-Advection test.')

  CASE ('DF1')
   ! The topography ist initialized in "init_nh_state_prog_dftest"
    CALL message(TRIM(routine),'running the deformational flow 2D-Advection test 1.')

  CASE ('DF2')
   ! The topography ist initialized in "init_nh_state_prog_dftest"
    CALL message(TRIM(routine),'running the deformational flow 2D-Advection test 2.')

  CASE ('DF3')
   ! The topography ist initialized in "init_nh_state_prog_dftest"
    CALL message(TRIM(routine),'running the deformational flow 2D-Advection test 3.')

  CASE ('DF4')
   ! The topography ist initialized in "init_nh_state_prog_dftest"
    CALL message(TRIM(routine),'running the deformational flow 2D-Advection test 4.')

  CASE ('HS_nh')
    ! The topography has been initialized to 0 
    CALL message(TRIM(routine),'running the Held-Suarez test')
    
  CASE ('APE_nwp')
    ! The topography has been initialized to 0 at the begining of this SUB
    CALL message(TRIM(routine),'running Aqua-Planet Experiment with non-hydrostatic atm. dynamics and NWP physics')
  
  CASE ('APE_aes')
    ! The topography has been initialized to 0 at the begining of this SUB
    CALL message(TRIM(routine),'running Aqua-Planet Experiment with non-hydrostatic atm. dynamics and AES physics')
  
  CASE ('APE_nh')
    ! The topography has been initialized to 0 at the begining of this SUB
    CALL message(TRIM(routine),'running Aqua-Planet Experiment with non-hydrostatic atm. dynamics')

  CASE ('APEc_nh')
    ! The topography has been initialized to 0 at the begining of this SUB
    CALL message(TRIM(routine),'running coupled Aqua-Planet Experiment with non-hydrostatic atm. dynamics')

  CASE ('TPEo')

    ! The topography has been initialized to 0 at the begining of this SUB
    CALL message(TRIM(routine),'running Terra-Planet Experiment with AES physics')
    IF ( itopo == 0 ) THEN
      CALL message(TRIM(routine), 'using zero topography for TPEc experiment')
    END IF

  CASE ('TPEc')

   ! The topography has been initialized to 0 at the begining of this SUB
    CALL message(TRIM(routine),'running Terra-Planet Experiment with AES physics')
    IF ( itopo == 0 ) THEN
      CALL message(TRIM(routine), 'using zero topography for TPEc experiment')
    END IF
  
  CASE ('g_lim_area')

    DO jg = 1, n_dom 
     nblks_c   = p_patch(jg)%nblks_c
     npromz_c  = p_patch(jg)%npromz_c

     CALL init_nh_topo_ana ( p_patch(jg), lplane, ext_data(jg)%atm%topography_c, nblks_c, npromz_c)

    END DO
    CALL message(TRIM(routine),'running g_lim_area')

  CASE ('dcmip_pa_12')

    ! The topography has been initialized to 0 
    CALL message(TRIM(routine),'running the dcmip_pa_12 (PA with Hadley-like circulation) test')

  CASE ('dcmip_gw_31')

    ! The topography has been initialized to 0 
    CALL message(TRIM(routine),'running the dcmip_gw_31 (small planet gravity wave) test')

  CASE ('dcmip_gw_32')

    ! The topography has been initialized to 0 
    CALL message(TRIM(routine),'running the dcmip_gw_32 (analyt. small planet gravity wave) test')

  CASE ('dcmip_rest_200')

    DO jg = 1, n_dom 

     CALL init_nh_topo_dcmip_rest_atm ( p_patch(jg), ext_data(jg)%atm%topography_c, ext_data(jg)%atm%fis  )

    END DO

    CALL message(TRIM(routine),'running the dcmip_rest_200 (steady state at rest dcmip) test')

!!$  CASE ('dcmip_mw_2x')
!DR topography_v no longer read in available. If needed, it should be 
!DR interpolated from topography_c.
!!$
!!$    DO jg = 1, n_dom 
!!$
!!$     CALL init_nh_topo_dcmip_schaer ( p_patch(jg),  ext_data(jg)%atm%topography_c,  &
!!$                          & ext_data(jg)%atm%topography_v, ext_data(jg)%atm%fis  )
!!$    END DO
!!$    CALL message(TRIM(routine),'running the dcmip_mw_2x (schaer-type dcmip) test')

  CASE ('dcmip_tc_51')
    ! itopo == 0 --> The topography is initialized to 0 at the begining of this subroutine
    CALL message(TRIM(routine),'running DCMIP tropical cyclone testcase 51')

  CASE ('dcmip_tc_52')
    ! itopo == 0 --> The topography is initialized to 0 at the begining of this subroutine
    CALL message(TRIM(routine),'running DCMIP tropical cyclone testcase 52')

  CASE ('CBL')

    IF(p_patch(1)%geometry_info%geometry_type/=planar_torus_geometry)&
        CALL finish(TRIM(routine),'CBL case is only for plane torus!')

   ! The topography has been initialized to 0 at the begining of this SUB
    CALL message(TRIM(routine),'running Convective Boundary Layer Experiment')

  CASE ('CBL_flxconst')

   CALL message(TRIM(routine),'running Convective Boundary Layer Experiment with fixed heat flux')

  CASE ('2D_BUBBLE', '3D_BUBBLE')

    IF(p_patch(1)%geometry_info%geometry_type/=planar_torus_geometry)&
        CALL finish(TRIM(routine),'2D BUBBLE case is only for plane torus!')

   ! The topography has been initialized to 0 at the begining of this SUB
    CALL message(TRIM(routine),'running 2D Warm Bubble test case')

  CASE ('RCE_glb')

   ! Running Radiative Convective Equilibrium testcase
   CALL message(TRIM(routine),'running ICON in RCE on a global domain')

  CASE ('RCE_Tconst')

   ! Running Radiative Convective Equilibrium testcase
     CALL message(TRIM(routine),'running ICON in RCE with constant initial T profile')

  CASE ('RCE_Tprescr')
   ! Running Radiative Convective Equilibrium with prescribed temperature profile
     CALL message(TRIM(routine),'running ICON in RCE with prescribed initial temperature profile')

  CASE ('aes_bubble', 'aes_bubble_land')
     CALL message(TRIM(routine),'running ICON torus with moisture and temperature bubble')
  CASE ('RICO')

    IF(p_patch(1)%geometry_info%geometry_type/=planar_torus_geometry)&
        CALL finish(TRIM(routine),'RICO case is only for plane torus!')

   ! The topography has been initialized to 0 at the begining of this SUB
    CALL message(TRIM(routine),'running Rain in the Culumus Over the Ocean LES Experiment')

  CASE ('RCE','GATE','RCEMIP_analytical')

    IF(p_patch(1)%geometry_info%geometry_type/=planar_torus_geometry)&
        CALL finish(TRIM(routine),'To initialize with sounding is only for torus!')

   ! The topography has been initialized to 0 at the begining of this SUB
    CALL message(TRIM(routine),'running LES with sounding')

  CASE DEFAULT

    CALL finish(routine,'wrong input for nh_test_name')

  END SELECT


  END SUBROUTINE init_nh_testtopo


!-------------------------------------------------------------------------
!
!
  !>
  !! Defines nonhydrostatic artificial initial conditions.
  !! (not for single column model (SCM))
  !! 
  !! Initializes meteorological fields
  !! 
  SUBROUTINE init_nh_testcase (p_patch, p_nh_state, p_int, p_lnd_state, ext_data, ntl)
!
! !INPUT VARIABLES:
  TYPE(t_patch),TARGET,  INTENT(INOUT) :: p_patch(n_dom)
  TYPE(t_int_state),     INTENT(   IN) :: p_int(n_dom)
  TYPE(t_lnd_state),     INTENT(INOUT) :: p_lnd_state(n_dom)
  TYPE(t_external_data), INTENT(INOUT) :: ext_data(n_dom)
  INTEGER :: ntl
  TYPE(t_nh_state), TARGET, INTENT(INOUT):: p_nh_state(n_dom)

  INTEGER        :: jg, je, jc, jb, jk, jt,   &
                    nlen, nblks_e, npromz_e,  nblks_c, npromz_c
  INTEGER        :: nlev, nlevp1        !< number of full and half levels

  TYPE(t_nh_state), POINTER       :: p_nhdom
  TYPE(t_cartesian_coordinates) :: p
                            
  REAL(wp)              :: p_sfc_jabw  ! surface pressure for the jabw test case, 
                                       ! standard values is 100000 Pa   
  REAL(wp)              :: global_moist
  REAL(wp) :: z_help

  LOGICAL  :: l_hydro_adjust, l_moist

  ! for a piecewise polytropic atmosphere:
  INTEGER ::     nmbr_polytropic_levels
  REAL(wp), ALLOCATABLE :: h_poly(:)
  REAL(wp), ALLOCATABLE :: T0_poly(:)
  REAL(wp), ALLOCATABLE :: dTdz_poly(:)
  REAL(wp) :: p_surf

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine =  &
                                   '(mo_nh_testcases) init_nh_testcase:' 

!-----------------------------------------------------------------------

  SELECT CASE (nh_test_name)

  CASE ('jabw', 'jabw_s', 'jabw_m', 'HS_jw') ! Jablonowski test

    CALL message(TRIM(routine),'Jablonowski test')
    IF ( iforcing == inwp ) THEN
      CALL message(TRIM(routine),' iforcing == inwp')
    ELSE
      CALL message(TRIM(routine),'Attention: iforcing /= inwp')
    ENDIF
  
    IF (nh_test_name == "jabw_s" .OR. nh_test_name == "jabw_m") THEN
      jw_up = 0.0_wp
    END IF  
      p_sfc_jabw = ps0
  
    DO jg = 1, n_dom
  
      CALL   init_nh_state_prog_jabw ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                     & p_nh_state(jg)%diag, p_nh_state(jg)%metrics, &
                                     & p_int(jg),                                   &
                                     & p_sfc_jabw,jw_up,jw_u0,jw_temp0 )
    
      IF ( ltransport .AND. iforcing /= inwp ) THEN   ! passive tracers
  
         CALL init_passive_tracers_nh_jabw (p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                           & rotate_axis_deg, tracer_inidist_list, p_sfc_jabw)

      END IF
  
      IF ( ltransport .AND. iforcing == inwp ) THEN 
        IF ( atm_phy_nwp_config(jg)%inwp_gscp /= 0 .OR.&
                     &                 atm_phy_nwp_config(jg)%inwp_convection /= 0  ) THEN   !
    
        CALL init_nh_inwp_tracers ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                  & p_nh_state(jg)%diag, p_nh_state(jg)%metrics, &
                                  & rh_at_1000hpa, qv_max, l_rediag=.TRUE. )

       ELSE
  
         p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,:) = 0.0_wp   
  
       END IF
  
      END IF
  
      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
  
    ENDDO !jg
  
    CALL message(TRIM(routine),'End setup Jablonowski test')


  CASE ('dcmip_bw_11')

    DO jg = 1, n_dom
      CALL init_nh_dcmip_bw (p_patch(jg),                   &
        &                    p_nh_state(jg)%prog(nnow(jg)), &
        &                    p_nh_state(jg)%diag,           &
        &                    p_int(jg),                     &
        &                    p_nh_state(jg)%metrics         )
    END DO  ! jg
    CALL message(TRIM(routine),'End setup baroclinic wave (dcmip_bw_11) test')


  CASE ('mrw_nh', 'mrw2_nh')

   CALL message(TRIM(routine),'MRW test')
  
   l_hydro_adjust = .TRUE.
   l_moist = .FALSE.

   DO jg = 1, n_dom

     IF ( iforcing == inwp ) THEN
       CALL message(TRIM(routine),' iforcing == inwp')     
       IF ( atm_phy_nwp_config(jg)%inwp_gscp /= 0 .OR.&
                      &                 atm_phy_nwp_config(jg)%inwp_convection /= 0  ) THEN 
         l_moist = .TRUE.
       END IF
     ENDIF


     IF (.NOT. l_moist) THEN

       CALL   init_nh_state_prog_mrw ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                     & p_nh_state(jg)%diag,                        &
                                     & ext_data(jg)%atm%topography_c,              &
                                     & p_nh_state(jg)%metrics,                     &
                                     & p_int(jg), l_hydro_adjust, iforcing, l_moist)

    ELSE

       CALL   init_nh_state_prog_mrw ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                     & p_nh_state(jg)%diag,                        &
                                     & ext_data(jg)%atm%topography_c,              &
                                     & p_nh_state(jg)%metrics,                     &
                                     & p_int(jg), l_hydro_adjust, iforcing, l_moist , &
                                     & opt_rh_at_1000hpa= rh_at_1000hpa,           &
                                     & opt_qv_max=qv_max                           )
    END IF
   
    CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))  
 
   ENDDO !jg

   CALL message(TRIM(routine),'End setup MRW test')


  CASE ('mwbr_const')

   CALL message(TRIM(routine),'mwbr_const test case')
  
   l_hydro_adjust = .TRUE.
   l_moist = .FALSE.

   DO jg = 1, n_dom

     IF ( iforcing == inwp ) THEN
       CALL message(TRIM(routine),' iforcing == inwp')     
       IF ( atm_phy_nwp_config(jg)%inwp_gscp /= 0 .OR.&
                      &                 atm_phy_nwp_config(jg)%inwp_convection /= 0  ) THEN 
         l_moist = .TRUE.
       END IF
     ENDIF 

     IF (.NOT. l_moist) THEN

       CALL   init_nh_prog_mwbr_const ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                      & p_nh_state(jg)%diag,                        &
                                      & ext_data(jg)%atm%topography_c,              &
                                      & p_nh_state(jg)%metrics,                     &
                                      & p_int(jg), l_hydro_adjust, iforcing, l_moist)

     ELSE
  
       CALL   init_nh_prog_mwbr_const ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                      & p_nh_state(jg)%diag,                        &
                                      & ext_data(jg)%atm%topography_c,              &
                                      & p_nh_state(jg)%metrics,                     &
                                      & p_int(jg), l_hydro_adjust, iforcing, l_moist,&
                                      & opt_rh_at_1000hpa= rh_at_1000hpa,           &
                                      & opt_qv_max=qv_max                           )

    END IF

    CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg))) 
 
   ENDDO !jg

   CALL message(TRIM(routine),'End setup mwbr_const test')


  CASE ('zero','bell','schaer', 'gauss3D', 'straka93' )

    ! For the moment we think of a given Brunt Vaisala frequency and a given      
    ! zonal wind. The lplane and the lcorio=F options are assumed

    DO jg = 1, n_dom
      p_nhdom   => p_nh_state(jg)
      nblks_e   = p_patch(jg)%nblks_e
      npromz_e  = p_patch(jg)%npromz_e
      nblks_c   = p_patch(jg)%nblks_c
      npromz_c  = p_patch(jg)%npromz_c

      ! number of vertical levels
      nlev   = p_patch(jg)%nlev
      nlevp1 = p_patch(jg)%nlevp1

      DO jt = 1, ntl 
        ! normal wind
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,nlen)
        DO jb = 1, nblks_e
          IF (jb /= nblks_e) THEN
             nlen = nproma
          ELSE
             nlen = npromz_e
          ENDIF
          DO jk = 1, nlev
            DO je = 1, nlen
              p_nh_state(jg)%prog(jt)%vn(je,jk,jb) = nh_u0 &
              !(p_nhdom%metrics%geopot(1,jk,1)/grav/1000.0_wp+5.0_wp)& !shear
              *p_patch(jg)%edges%primal_normal(je,jb)%v1

              ! copy vn to reference state vector (needed by Rayleigh damping mechanism)
              p_nh_state(jg)%ref%vn_ref(je,jk,jb)  &
                        = p_nh_state(jg)%prog(jt)%vn(je,jk,jb)
            ENDDO
          ENDDO
        ENDDO  !jb
!$OMP END DO NOWAIT

      ! scalars (all is dry!)
!$OMP DO PRIVATE(jb,jk,jc,nlen,z_help) 
        DO jb = 1, nblks_c
          IF (jb /= nblks_c) THEN
             nlen = nproma
          ELSE
             nlen = npromz_c
          ENDIF
          DO jk = 1, nlev
            DO jc = 1, nlen
              z_help=(nh_brunt_vais/grav)**2*p_nhdom%metrics%geopot(jc,jk,jb)
              ! profile of theta is explicitly given
              p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb) = nh_t0*EXP(z_help)
            ENDDO
          ENDDO
          ! lower boundary condition for exner pressure
          !
          DO jc = 1, nlen
            IF (nh_brunt_vais /= 0.0_wp) THEN
              ! for exner(nlev) (lowermost level): analytical exner
              z_help = (nh_brunt_vais/grav)**2*p_nhdom%metrics%geopot(jc,nlev,jb)
              p_nh_state(jg)%prog(jt)%exner(jc,nlev,jb) =    &
                &  (grav/nh_brunt_vais)**2/nh_t0/cpd*(EXP(-z_help)-1.0_wp)+1.0_wp
            ELSE
              p_nh_state(jg)%prog(jt)%exner(jc,nlev,jb) = &
                &  1.0_wp-p_nhdom%metrics%geopot(jc,nlev,jb)/cpd/nh_t0
            ENDIF
          ENDDO
        ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

        ! compute hydrostatically balanced exner, by integrating the (discretized!) 
        ! 3rd equation of motion under the assumption thetav=const.
        CALL hydro_adjust_const_thetav(p_patch           = p_patch(jg),                   &
          &                            p_nh_metrics      = p_nhdom%metrics,               &
          &                            lintegrate_topdown= .FALSE.,                       &
          &                            rho               = p_nh_state(jg)%prog(jt)%rho,   &
          &                            exner             = p_nh_state(jg)%prog(jt)%exner, &
          &                            theta_v           = p_nh_state(jg)%prog(jt)%theta_v)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,nlen,p,z_help)
        DO jb = 1, nblks_c
          IF (jb /= nblks_c) THEN
             nlen = nproma
          ELSE
             nlen = npromz_c
          ENDIF
          DO jk = 1, nlev
            DO jc = 1, nlen

              !     ! perturbation in theta_v for gravity test case
              !     p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb)=&
              !              p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb)+ 0.01_wp*sin(&
              !              pi*p_nhdom%metrics%geopot(jc,jk,jb)/grav/10000.0_wp)&
              !            /(1.0_wp+(p_patch(jg)%cells%center(jc,jb)%lon*30.0/pi)**2)
              !            !Das ist fuer dx=500 mit 600 Punkten (auf 2pi verteilt)

              IF ( nh_test_name == "straka93 " ) THEN
                ! test case Straka et al. (1993) (falling cold bubble)
                ! perturbation in theta_v

                !     z_help = SQRT( ( p_patch(jg)%cells%center(jc,jb)%lon*6.4_wp/pi)**2 &
                !     &      +((p_nhdom%metrics%z_mc(jc,jk,jb) - bubctr_z)/bub_ver_width)**2)
                !     IF (z_help<=1.0_wp) THEN
                !       p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb)=      &
                !       &   p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb)   &
                !       &   + bub_amp * (COS(pi*z_help)+1.0_wp)*0.5_wp  &
                !       &   /p_nh_state(jg)%prog(jt)%exner(jc,jk,jb)
                !     ENDIF

                p = gc2cc( p_patch(jg)%cells%center(jc,jb), p_patch(1)%geometry_info )
                z_help = SQRT( ( p%x(1)/bub_hor_width )**2                   &
                  &          + ( (p_nhdom%metrics%z_mc(jc,jk,jb) - bubctr_z)/bub_ver_width )**2 )
                IF ( z_help <= 1.0_wp ) THEN
                  p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb)=             &
                    &   p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb)        &
                    &    + bub_amp * (COS(pi*z_help) + 1.0_wp) * 0.5_wp  &
                    &      /p_nh_state(jg)%prog(jt)%exner(jc,jk,jb)
                ENDIF
               !PRINT*,"coordinate z is=z", p_nhdom%metrics%z_mc(jc,jk,jb)
              END IF

              ! exner and theta_v are given, so rho is deduced...
              p_nh_state(jg)%prog(jt)%rho(jc,jk,jb) = &
                &        (p_nh_state(jg)%prog(jt)%exner(jc,jk,jb)**cvd_o_rd)*p0ref/rd &
                &       /p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb)

            ENDDO
          ENDDO
          DO jk = 1, nlevp1
            p_nh_state(jg)%prog(jt)%w(1:nlen,jk,jb) = 0.0_wp

            ! copy w to reference state vector (needed by Rayleigh damping mechanism)
            p_nh_state(jg)%ref%w_ref(1:nlen,jk,jb) = &
                p_nh_state(jg)%prog(jt)%w(1:nlen,jk,jb)
          ENDDO
        ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL
      ENDDO  ! jt
    ENDDO  ! jg


  CASE ( 'atm_at_rest' )

    ! set a piecewise polytropic atmosphere

    nmbr_polytropic_levels = 2
    ALLOCATE( h_poly   (nmbr_polytropic_levels + 1) )
    ALLOCATE( T0_poly  (nmbr_polytropic_levels) )
    ALLOCATE( dTdz_poly(nmbr_polytropic_levels) )

    h_poly(1)    =     0.0_wp
    T0_poly(1)   = 298.15_wp
    dTdz_poly(1) = -0.0065_wp
    p_surf       = 1.0e5_wp

    h_poly(2)    = 12000.0_wp
    T0_poly(2)   = 220.15_wp
    dTdz_poly(2) = 0.0_wp

    h_poly(3)    = vct_a(1)   ! something equal or higher than the model top height
 
    CALL piecewise_polytropic_atm( p_patch,  p_nh_state,    &
        &       nmbr_polytropic_levels, h_poly, T0_poly, dTdz_poly, p_surf, ntl, .TRUE. )


  CASE ('PA')  ! pure advection test case, no mountain

    DO jg = 1, n_dom
      CALL init_nh_state_prog_patest(p_patch(jg),p_int(jg),             &
        &                       p_nh_state(jg)%prog(nnow(jg)),          &
        &                       p_nh_state(jg)%diag,ext_data(jg),       &
        &                       p_nh_state(jg)%metrics,rotate_axis_deg, &
        &                       linit_tracer_fv, tracer_inidist_list )

      CALL init_nh_state_prog_patest(p_patch(jg),p_int(jg),             &
        &                       p_nh_state(jg)%prog(nnew(jg)),          &
        &                       p_nh_state(jg)%diag,ext_data(jg),       &
        &                       p_nh_state(jg)%metrics,rotate_axis_deg, &
        &                       linit_tracer_fv, tracer_inidist_list )

    ENDDO !jg


  CASE ('DF1', 'DF2', 'DF3', 'DF4')  ! 2D deformational flow test case, no mountain

    DO jg = 1, n_dom

      CALL init_nh_state_prog_dftest(p_patch(jg),                         &
           &                         p_nh_state(jg)%prog(nnow(jg)),       &
           &                         p_nh_state(jg)%diag,                 &
           &                         p_int(jg), ext_data(jg),             &
           &                         rotate_axis_deg, nh_test_name,       &
           &                         linit_tracer_fv, tracer_inidist_list )

      CALL init_nh_state_prog_dftest(p_patch(jg),                         &
           &                         p_nh_state(jg)%prog(nnew(jg)),       &
           &                         p_nh_state(jg)%diag,                 &
           &                         p_int(jg), ext_data(jg),             &
           &                         rotate_axis_deg, nh_test_name,       &
           &                         linit_tracer_fv, tracer_inidist_list )

    ENDDO !jg


  CASE ('HS_nh')  ! Held-Suarez test case, no mountain, isothermal atmosphere

    DO jg = 1, n_dom

      ! number of vertical levels
      nlev   = p_patch(jg)%nlev

      CALL init_nh_state_prog_held_suarez(p_patch(jg),p_nh_state(jg)%prog(nnow(jg)),  &
        &                            p_nh_state(jg)%diag,ext_data(jg),           &
        &                            p_nh_state(jg)%metrics)

      IF (lhs_nh_vn_ptb) THEN
          CALL nh_prog_add_random( p_patch(jg),           & ! input
               & p_nh_state(jg)%prog(nnow(jg))%vn(:,:,:), & ! in and out
               & "edge", hs_nh_vn_ptb_scale, 1, nlev ) ! input

          CALL message(TRIM(routine),'Initial state used in the &
               & Held-Suarez test: random noised added to the normal wind')
      END IF

      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))

    ENDDO !jg


  CASE ('APE_nwp', 'APE_aes', 'APE_nh', 'APEc_nh')  ! Aqua-Planet Experiment, no mountain

    p_sfc_jabw   = zp_ape          ! Pa
    global_moist = ztmc_ape        ! kg/m**2 total moisture content

    DO jg = 1, n_dom
    
      CALL   init_nh_state_prog_jabw ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                     & p_nh_state(jg)%diag, p_nh_state(jg)%metrics, &
                                     & p_int(jg),                                   &
                                     & p_sfc_jabw,jw_up,jw_u0,jw_temp0 )
    
      IF ( ltransport ) THEN   !
    
        CALL init_nh_inwp_tracers ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                  & p_nh_state(jg)%diag, p_nh_state(jg)%metrics, &
                                  & rh_at_1000hpa, qv_max, l_rediag=.TRUE.,  &
                                  & opt_global_moist=global_moist)
    
      END IF
    
      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
    
    ENDDO !jg

    CALL message(TRIM(routine),'End setup non-hydrostatic APE test (APE_nwp, APE_aes, APE_nh, APEc_nh)')


  CASE ('TPEc', 'TPEo')  ! Terra-Planet Experiment

    jw_up = 1._wp

    DO jg = 1, n_dom
    
      CALL init_nh_state_prog_TPE(p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), p_nh_state(jg)%diag, &
                                  ext_data(jg), p_nh_state(jg)%metrics,                            &
                                  rh_at_1000hpa, qv_max, tpe_moist, tpe_psfc, tpe_temp)

      ! why do we call this?   
      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
    
    ENDDO !jg

    CALL message(TRIM(routine),'End setup TPEc test')


  CASE ('wk82')

   CALL message(TRIM(routine),'wk82 test')
  
   l_hydro_adjust = .TRUE.

   DO jg = 1, n_dom

         ! initialize environment atmosphere
  
    CALL   init_nh_env_wk ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                     & p_nh_state(jg)%diag,                 &
                                     & p_nh_state(jg)%metrics,              &
                                     & p_int(jg), l_hydro_adjust )
         ! add perturbation to theta and recalculate theta_v and rho 
    CALL init_nh_buble_wk ( p_patch(jg),p_nh_state(jg)%metrics,            &
                                     & p_nh_state(jg)%prog(nnow(jg)),      &
                                     & p_nh_state(jg)%diag )
         

    CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))

   ENDDO !jg

   CALL message(TRIM(routine),'End setup wk82 test')

  
  CASE ('bb13')

    CALL message(TRIM(routine), 'Baldauf, Brdar (2013) QJRMS test (linear gravity/sound waves in a channel)')
  
    l_hydro_adjust = .TRUE.

    DO jg = 1, n_dom

         ! initialize environment atmosphere
      CALL init_nh_env_bb13   ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                     & p_nh_state(jg)%diag,                 &
                                     & p_nh_state(jg)%metrics,              &
                                     & l_hydro_adjust  )
         ! add perturbation to theta and recalculate theta_v and rho 
      CALL init_nh_bubble_bb13 ( p_patch(jg), p_nh_state(jg)%metrics,       &
                                     & p_nh_state(jg)%prog(nnow(jg)),       &
                                     & p_nh_state(jg)%diag )
         

      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))

    ENDDO !jg

    CALL message(TRIM(routine),'End setup Baldauf, Brdar (2013) test')

  
  CASE ('g_lim_area')

   CALL message(TRIM(routine),'g_lim_area test')
  
   l_hydro_adjust = .TRUE.

   IF (.NOT. l_limited_area  .OR. lcoriolis) THEN

     WRITE(message_text,'(a)') &
             & 'For g_lim_area test case l_limited_area must &
             & be .TRUE. and lcoriolis must be .FALSE.'
            CALL finish  (routine, TRIM(message_text))
   END IF

   DO jg = 1, n_dom

         ! initialize environment atmosphere
    SELECT CASE (itype_atmo_ana)

    CASE(1)
    
      CALL  init_nh_atmo_ana_nconstlayers( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                     & p_nh_state(jg)%diag,                 &
                                     & p_nh_state(jg)%metrics,              &
                                     & l_hydro_adjust  )
    CASE(2) 
      CALL  init_nh_atmo_ana_poly( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                     & p_nh_state(jg)%diag,                 &
                                     & p_nh_state(jg)%metrics,              &
                                     & l_hydro_adjust  )

    CASE default
      WRITE(message_text,'(a)') &
             & 'You should define a valid option for    &
             & itype_atmo_ana for the &
             & g_lim_area test case'
            CALL finish  (routine, TRIM(message_text))

    END SELECT
 
         ! initialize wind
    CALL   init_nh_anaprof_uv  ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg))%vn, &
                               & p_nh_state(jg)%prog(nnow(jg))%w,               &
                               & p_nh_state(jg)%metrics, p_int(jg) ) 
         
    CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))

   ENDDO !jg

   CALL message(TRIM(routine),'End setup g_lim_area test')


  CASE ('dcmip_pa_12')

    CALL message(TRIM(routine),'setup dcmip_pa_12 (PA with Hadley-like circulation) test')

    DO jg = 1, n_dom
      CALL init_nh_dcmip_hadley( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)),  &
        &                        p_nh_state(jg)%diag, p_int(jg),              &
        &                        p_nh_state(jg)%metrics )

      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
    ENDDO

    CALL message(TRIM(routine),'End setup dcmip_pa_12 test')


  CASE ('dcmip_gw_31')

    CALL message(TRIM(routine),'setup dcmip_gw_31 (gravity waves on small planet) test')

    l_hydro_adjust = .TRUE.

    DO jg = 1, n_dom
      CALL init_nh_dcmip_gw( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)),  &
        &                    p_nh_state(jg)%diag, p_nh_state(jg)%metrics, & 
        &                    l_hydro_adjust)
    ENDDO

    CALL message(TRIM(routine),'End setup dcmip_gw_31 test')


  CASE ('dcmip_gw_32')

    CALL message(TRIM(routine),'setup dcmip_gw_32 (analyt. gravity waves on small planet) test')

    DO jg = 1, n_dom
      CALL init_nh_gw_analyt( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)),           &
        &                    p_nh_state(jg)%diag, p_nh_state(jg)%metrics, p_int(jg) )
    ENDDO

    CALL message(TRIM(routine),'End setup dcmip_gw_32 test')


  CASE ('dcmip_rest_200')

    CALL message(TRIM(routine),'setup dcmip_rest_200 (steady state at rest) test')
  
    l_hydro_adjust = .TRUE.

    IF ( lcoriolis) THEN

      WRITE(message_text,'(a)') &
             & 'For dcmip_rest_200 test case  lcoriolis must be .FALSE.'
            CALL finish  (routine, TRIM(message_text))
    END IF

    DO jg = 1, n_dom
      CALL init_nh_prog_dcmip_rest_atm( p_patch(jg),                              &
        &                    p_nh_state(jg)%prog(nnow(jg)), p_nh_state(jg)%diag,  &
        &                    p_nh_state(jg)%metrics, l_hydro_adjust )
    CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
    ENDDO

    CALL message(TRIM(routine),'End setup dcmip_rest_200 test')


  CASE ('dcmip_mw_2x')

    CALL message(TRIM(routine),'setup dcmip_mw_2x (schaer-type on small planet) test')
  
    l_hydro_adjust = .FALSE.

    IF ( lcoriolis) THEN

      WRITE(message_text,'(a)') &
             & 'For dcmip_mw_2x test case  lcoriolis must be .FALSE.'
            CALL finish  (routine, TRIM(message_text))
    END IF

    DO jg = 1, n_dom
      CALL init_nh_prog_dcmip_schaer( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
        &                    p_nh_state(jg)%diag, p_nh_state(jg)%ref,             &
        &                    p_nh_state(jg)%metrics, p_int(jg),l_hydro_adjust )
    CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
    ENDDO

    CALL message(TRIM(routine),'End setup dcmip_mw_2x test')


  CASE ('dcmip_tc_51','dcmip_tc_52')

    ! 'dcmip_tc_51' and 'dcmip_tc_52' have the same initial state.

    DO jg = 1, n_dom

      CALL init_nh_dcmip_tc ( p_patch(jg),                   &
        &                     p_nh_state(jg)%prog(nnow(jg)), &
        &                     p_nh_state(jg)%diag,           &
        &                     p_nh_state(jg)%metrics,        &
        &                     p_int(jg)                     )

      CALL duplicate_prog_state( p_nh_state(jg)%prog(nnow(jg)), &
        &                        p_nh_state(jg)%prog(nnew(jg)) )

    END DO !jg

    CALL message(TRIM(routine),'End setup dcmip_tc_51/52')


  CASE ('CBL')

    IF(p_patch(1)%geometry_info%geometry_type/=planar_torus_geometry)&
        CALL finish(TRIM(routine),'CBL case is only for plane torus!')

    DO jg = 1, n_dom
      nlev   = p_patch(jg)%nlev
      CALL init_nh_state_cbl ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), p_nh_state(jg)%ref,  &
                      & p_nh_state(jg)%diag, p_int(jg), p_nh_state(jg)%metrics )
 
      CALL add_random_noise_global(in_subset=p_patch(jg)%cells%all,            &
                      & in_var=p_nh_state(jg)%prog(nnow(jg))%w(:,:,:),         &
                      & start_level=nlev-3,                                    &
                      & end_level=nlev,                                        &
                      & noise_scale=w_perturb )   

      CALL add_random_noise_global(in_subset=p_patch(jg)%cells%all,            &
                      & in_var=p_nh_state(jg)%prog(nnow(jg))%theta_v(:,:,:),   &
                      & start_level=nlev-3,                                    &
                      & end_level=nlev,                                        &
                      & noise_scale=th_perturb )   

      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
    END DO !jg

    CALL message(TRIM(routine),'End setup CBL test')


  CASE ('CBL_flxconst')

    ! u,v,w are initialized to zero.  exner and rho are similar/identical to CBL
    DO jg = 1, n_dom
      nlev   = p_patch(jg)%nlev

      CALL init_nh_state_cbl ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), p_nh_state(jg)%ref,  &
                      & p_nh_state(jg)%diag, p_int(jg), p_nh_state(jg)%metrics )

      CALL add_random_noise_global(in_subset=p_patch(jg)%cells%all,            &
                      & in_var=p_nh_state(jg)%prog(nnow(jg))%theta_v(:,:,:),   &
                      & start_level=nlev-3,                                    &
                      & end_level=nlev,                                        &
                      & noise_scale=th_perturb )

      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
!
      CALL message(TRIM(routine),'End setup global CBL_flxconst test')
    END DO !jg


  CASE ('RCE_glb','RCE_Tconst')

    ! u,v,w are initialized to zero.  exner and rho are similar/identical to CBL
    DO jg = 1, n_dom
      nlev   = p_patch(jg)%nlev
      CALL init_nh_state_rce_glb ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), p_nh_state(jg)%ref,  &
                      & p_nh_state(jg)%diag, p_nh_state(jg)%metrics )

      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
!
      CALL message(TRIM(routine),'End setup global RCE test')
    END DO !jg


  CASE ('RCE_Tprescr')

     ! u,v,w are initialized to zero.  initialize with temperature profile, add random noise
     ! to virtual potential temperature inside init_nh_state_rce_tprescr_glb.
    DO jg = 1, n_dom
      nlev   = p_patch(jg)%nlev
      CALL message(TRIM(routine),'before init_nh_state_rce_tprescr_glb')
      CALL init_nh_state_rce_tprescr_glb ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), p_nh_state(jg)%ref,  &
                      & p_nh_state(jg)%diag, p_nh_state(jg)%metrics )
      
      
      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
!
      CALL message(TRIM(routine),'End setup global RCE_Tprescr test')
    END DO !jg

  CASE ('aes_bubble', 'aes_bubble_land')

     ! u,v,w are initialized to zero.  initialize with temperature profile, add bubble to T and q
    DO jg = 1, n_dom
       nlev   = p_patch(jg)%nlev
      CALL message(TRIM(routine),'before init_aes_bubble')
      CALL init_aes_bubble ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), p_nh_state(jg)%ref,  &
                      & p_nh_state(jg)%diag, p_nh_state(jg)%metrics )
      
      
      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
!
      CALL message(TRIM(routine),'End setup '//TRIM(nh_test_name)//' test')
    END DO !jg
    

  CASE ('RICO')

    IF(p_patch(1)%geometry_info%geometry_type/=planar_torus_geometry)&
        CALL finish(TRIM(routine),'RICO case is only for plane torus!')

    DO jg = 1, n_dom
      nlev   = p_patch(jg)%nlev
      CALL init_nh_state_rico ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), p_nh_state(jg)%ref,  &
                      & p_nh_state(jg)%diag, p_int(jg), p_nh_state(jg)%metrics )

      CALL add_random_noise_global(in_subset=p_patch(jg)%cells%all,            &
                      & in_var=p_nh_state(jg)%prog(nnow(jg))%w(:,:,:),         &
                      & start_level=nlev-3,                                    &
                      & end_level=nlev,                                        &
                      & noise_scale=w_perturb )   

      CALL add_random_noise_global(in_subset=p_patch(jg)%cells%all,            &
                      & in_var=p_nh_state(jg)%prog(nnow(jg))%theta_v(:,:,:),   &
                      & start_level=nlev-3,                                    &
                      & end_level=nlev,                                        &
                      & noise_scale=th_perturb )   
 
      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
    END DO !jg

    CALL message(TRIM(routine),'End setup RICO test')


  CASE ('RCE','GATE') !to initialize with sounding

    IF(p_patch(1)%geometry_info%geometry_type/=planar_torus_geometry)&
        CALL finish(TRIM(routine),'GATE/RCE is only for torus!')

    DO jg = 1, n_dom
      nlev   = p_patch(jg)%nlev

      IF ( i_scm_netcdf > 0 ) THEN
        CALL init_torus_netcdf_sounding ( p_patch(jg),                         &
                      & p_nh_state(jg)%prog(nnow(jg)),                         &
                      & p_nh_state(jg)%ref, p_nh_state(jg)%diag, p_int(jg),    &
                      & p_nh_state(jg)%metrics, ext_data(jg) )
      ELSE
        CALL init_torus_ascii_sounding ( p_patch(jg),                          &
                      & p_nh_state(jg)%prog(nnow(jg)),                         &
                      & p_nh_state(jg)%ref, p_nh_state(jg)%diag, p_int(jg),    &
                      & p_nh_state(jg)%metrics )
      END IF

      CALL add_random_noise_global(in_subset=p_patch(jg)%cells%all,            &
                      & in_var=p_nh_state(jg)%prog(nnow(jg))%w(:,:,:),         &
                      & start_level=nlev-3,                                    &
                      & end_level=nlev,                                        &
                      & noise_scale=w_perturb )   

      CALL add_random_noise_global(in_subset=p_patch(jg)%cells%all,            &
                      & in_var=p_nh_state(jg)%prog(nnow(jg))%theta_v(:,:,:),   &
                      & start_level=nlev-3,                                    &
                      & end_level=nlev,                                        &
                      & noise_scale=th_perturb )   
                    
      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
    END DO !jg

    CALL message(TRIM(routine),'End init with sounding')


  CASE ('RCEMIP_analytical') !to initialize with analytical sounding

    IF(p_patch(1)%geometry_info%geometry_type/=planar_torus_geometry)&
        CALL finish(TRIM(routine),'To initizialize with sounding is only for torus!')

    DO jg = 1, n_dom
      nlev   = p_patch(jg)%nlev

      CALL init_torus_rcemip_analytical_sounding ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                 p_nh_state(jg)%ref, p_nh_state(jg)%diag, p_int(jg), p_nh_state(jg)%metrics )

!      CALL add_random_noise_global(in_subset=p_patch(jg)%cells%all,            &
!                      & in_var=p_nh_state(jg)%prog(nnow(jg))%w(:,:,:),         &
!                      & start_level=nlev-5,                                    &
!                      & end_level=nlev,                                        &
!                      & noise_scale=w_perturb )

      CALL add_random_noise_global(in_subset=p_patch(jg)%cells%all,            &
                      & in_var=p_nh_state(jg)%prog(nnow(jg))%theta_v(:,:,:),   &
                      & start_level=nlev-5,                                    &
                      & end_level=nlev,                                        &
                      & noise_scale=th_perturb )

      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
    END DO !jg


  CASE ('2D_BUBBLE', '3D_BUBBLE') !to initialize with sounding

    IF(p_patch(1)%geometry_info%geometry_type/=planar_torus_geometry)&
        CALL finish(TRIM(routine),'2D warm bubble case is only for torus!')

    DO jg = 1, n_dom
      nlev   = p_patch(jg)%nlev

      CALL init_warm_bubble ( p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                              p_nh_state(jg)%diag, p_nh_state(jg)%metrics )

      CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
    END DO !jg

    CALL message(TRIM(routine),'End initilization of 2D warm bubble')

  END SELECT


  ! Is current testcase subject to update during integration?
  SELECT CASE(TRIM(nh_test_name))
  ! Testcases which (potentially) require an update:
  CASE ("PA", "DF1", "DF2", "DF3", "DF4", "DCMIP_PA_12", "dcmip_pa_12")
    ltestcase_update = .TRUE.
  CASE default
    ltestcase_update = .FALSE.
  END SELECT


  IF ( ANY( (/icosmo/)==atm_phy_nwp_config(1)%inwp_turb ) .AND. &
     (nh_test_name=='APE_nwp' .OR. nh_test_name=='CBL' .OR. nh_test_name=='GATE' &
     .OR. nh_test_name=='RICO') ) THEN
    DO jg = 1, n_dom
      p_lnd_state(jg)%prog_lnd(nnow(jg))%t_g                    = th_cbl(1)
    END DO !jg
    IF (atm_phy_nwp_config(1)%inwp_surface > 0) THEN ! Fields are not allocated otherwise
      DO jg = 1, n_dom
        !Snow and sea ice initialization
        p_lnd_state(jg)%prog_lnd(nnow(jg))%t_g_t                  = th_cbl(1)
        p_lnd_state(jg)%prog_lnd(nnow(jg))%t_snow_t(:,:,:)        = th_cbl(1) !snow
        p_lnd_state(jg)%prog_lnd(nnow(jg))%t_g_t(:,:,isub_seaice) = th_cbl(1) !sea ice
        p_lnd_state(jg)%prog_wtr(nnow(jg))%t_ice(:,:)             = th_cbl(1) !sea ice
      END DO !jg
    ENDIF
  END IF

  ! Terminator toy chemistry
  ! possible add on for various test cases
  IF (is_toy_chem) THEN
    DO jg = 1, n_dom 
      CALL init_nh_dcmip_terminator (p_patch(jg),              &
        &                            p_nh_state(jg)%metrics,   &
        &                            p_nh_state(jg)%prog(:),   &
        &                            p_nh_state(jg)%diag       )
    ENDDO
    CALL message(TRIM(routine),'End setup terminator chemistry')
  ENDIF

  END SUBROUTINE init_nh_testcase


!-------------------------------------------------------------------------
!
!
  !>
  !! Defines nonhydrostatic artificial initial conditions.
  !! (for single column model (SCM))
  !! 
  !! Initializes meteorological fields
  !! 
  SUBROUTINE init_nh_testcase_scm (p_patch, p_nh_state, p_int, p_lnd_state, ext_data)
!
! !INPUT VARIABLES:
  TYPE(t_patch),TARGET,  INTENT(INOUT) :: p_patch(n_dom)
  TYPE(t_int_state),     INTENT(   IN) :: p_int(n_dom)
  TYPE(t_lnd_state),     INTENT(INOUT) :: p_lnd_state(n_dom)
  TYPE(t_nh_state), TARGET, INTENT(INOUT):: p_nh_state(n_dom)
  TYPE(t_external_data), INTENT(INOUT) :: ext_data(n_dom)

  INTEGER        :: jg
  INTEGER        :: nlev               !< number of full and half levels
                            
  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine =  &
                                   '(mo_nh_testcases) init_nh_testcase:' 

!-----------------------------------------------------------------------


  IF(p_patch(1)%geometry_info%geometry_type/=planar_torus_geometry)&
      CALL finish(TRIM(routine),'SCM is only for torus!')

! initialize with sounding

  DO jg = 1, n_dom
    nlev   = p_patch(jg)%nlev

    CALL init_torus_netcdf_sounding ( p_patch(jg),                           &
                    & p_nh_state(jg)%prog(nnow(jg)),                         &
                    & p_nh_state(jg)%ref, p_nh_state(jg)%diag, p_int(jg),    &
                    & p_nh_state(jg)%metrics, ext_data(jg) )

! add random noise

    IF (lscm_random_noise) THEN !add random noise
      CALL add_random_noise_global(in_subset=p_patch(jg)%cells%all,          &
                    & in_var=p_nh_state(jg)%prog(nnow(jg))%w(:,:,:),         &
                    & start_level=nlev-3,                                    &
                    & end_level=nlev,                                        &
                    & noise_scale=w_perturb )   

      CALL add_random_noise_global(in_subset=p_patch(jg)%cells%all,          &
                    & in_var=p_nh_state(jg)%prog(nnow(jg))%theta_v(:,:,:),   &
                    & start_level=nlev-3,                                    &
                    & end_level=nlev,                                        &
                    & noise_scale=th_perturb )   
    ENDIF
                  
    CALL duplicate_prog_state(p_nh_state(jg)%prog(nnow(jg)),p_nh_state(jg)%prog(nnew(jg)))
  END DO !jg

  CALL message(TRIM(routine),'End init with sounding SCM')

! Is current testcase subject to update during integration?

  ltestcase_update = .FALSE.

  IF ( ANY( (/icosmo/)==atm_phy_nwp_config(1)%inwp_turb ) ) THEN
    DO jg = 1, n_dom
      p_lnd_state(jg)%prog_lnd(nnow(jg))%t_g                    = th_cbl(1)
    END DO !jg
    IF (atm_phy_nwp_config(1)%inwp_surface > 0) THEN ! Fields are not allocated otherwise
      DO jg = 1, n_dom
        !Snow and sea ice initialization
        p_lnd_state(jg)%prog_lnd(nnow(jg))%t_g_t                  = th_cbl(1)
        p_lnd_state(jg)%prog_lnd(nnow(jg))%t_snow_t(:,:,:)        = th_cbl(1) !snow
        p_lnd_state(jg)%prog_lnd(nnow(jg))%t_g_t(:,:,isub_seaice) = th_cbl(1) !sea ice
        p_lnd_state(jg)%prog_wtr(nnow(jg))%t_ice(:,:)             = th_cbl(1) !sea ice
      END DO !jg
    ENDIF
  END IF

! Terminator toy chemistry
! possible add on for various test cases

  IF (is_toy_chem) THEN
    DO jg = 1, n_dom 
      CALL init_nh_dcmip_terminator (p_patch(jg),              &
        &                            p_nh_state(jg)%metrics,   &
        &                            p_nh_state(jg)%prog(:),   &
        &                            p_nh_state(jg)%diag       )
    ENDDO
    CALL message(TRIM(routine),'End setup terminator chemistry')
  ENDIF

  END SUBROUTINE init_nh_testcase_scm
































  SUBROUTINE piecewise_polytropic_atm( p_patch,  p_nh_state,  &
    &       nmbr_polytropic_levels, h_poly, T0_poly, dTdz_poly, p_surf, ntl,  &
    &       l_hydro_adjust )

    !
    ! calculate the fields
    ! [output]
    !   p_nh_state(jg)%prog(jt)%rho(:,:,:))
    !   p_nh_state(jg)%prog(jt)%exner(:,:,:)
    !   p_nh_state(jg)%prog(jt)%theta_v(:,:,:)
    ! for a piecewise polytropic atmosphere, defined by 
    ! [input]
    !   nmbr_polytropic_levels, h_poly(:), T0_poly(:), dTdz_poly(:), p_surf
    ! Additionally set
    ! [output]
    !   p_nh_state(jg)%prog(jt)%w(:,:,:)
    !   p_nh_state(jg)%prog(jt)%vn(:,:,:)
    !

    USE mo_mpi,                   ONLY: get_my_mpi_work_id, get_glob_proc0
    USE mo_hydro_adjust,          ONLY: hydro_adjust
    USE mo_nh_diagnose_pres_temp, ONLY: diagnose_pres_temp

    IMPLICIT NONE

    TYPE(t_patch),    TARGET, INTENT(inout) :: p_patch(n_dom)
    TYPE(t_nh_state), TARGET, INTENT(inout) :: p_nh_state(n_dom)

    INTEGER, INTENT(in) :: nmbr_polytropic_levels
    REAL(wp), INTENT(in) :: h_poly   ( nmbr_polytropic_levels+1 )
    REAL(wp), INTENT(in) :: T0_poly  ( nmbr_polytropic_levels+1 )
    REAL(wp), INTENT(in) :: dTdz_poly( nmbr_polytropic_levels+1 )
    REAL(wp), INTENT(in) :: p_surf

    INTEGER, INTENT(in) :: ntl

    LOGICAL, INTENT(in) :: l_hydro_adjust

    INTEGER :: jg, jt, jb, jk, jc, je, nlen
    INTEGER :: nlev, nlevp1
    INTEGER :: nblks_c, npromz_c
    INTEGER :: nblks_e, npromz_e
    INTEGER :: l

    REAL(wp), ALLOCATABLE :: p0_poly(:)
    REAL(wp) :: epsilon_dTdz
    REAL(wp) :: h1
    REAL(wp) :: z
    REAL(wp) :: temp, pres
    REAL(wp) :: delta

    epsilon_dTdz = 1.0e-10_wp

    ALLOCATE( p0_poly (nmbr_polytropic_levels+1 ) )

    ! polytropic pressures at the height intervals 

    p0_poly(1) = p_surf
    DO l=1, nmbr_polytropic_levels-1
      IF ( ABS( dTdz_poly(l) ) > epsilon_dTdz ) THEN
        h1 = 1.0_wp + dTdz_poly(l) / T0_poly(l) * ( h_poly(l+1) - h_poly(l) )
        p0_poly(l+1) = p0_poly(l) * h1**( - grav/(Rd* dTdz_poly(l) ) )
      ELSE
        ! isothermal 
        delta = grav / ( Rd * T0_poly(l) )
        p0_poly(l+1) = p0_poly(l) * EXP( - delta * ( h_poly(l+1) - h_poly(l) ) )
      END IF
    END DO

    DO jg = 1, n_dom

      ! number of vertical levels
      nlev   = p_patch(jg)%nlev
      nlevp1 = p_patch(jg)%nlevp1

      ! center:
      nblks_c   = p_patch(jg)%nblks_c
      npromz_c  = p_patch(jg)%npromz_c

      ! scalars (all is dry!)
      DO jt = 1, ntl 
        DO jb = 1, nblks_c

          IF (jb /= nblks_c) THEN
            nlen = nproma
          ELSE
            nlen = npromz_c
          ENDIF

          DO jk = 1, nlev
            DO jc = 1, nlen

              z = p_nh_state(jg)%metrics%z_mc(jc,jk,jb)

              DO l=1, nmbr_polytropic_levels

                IF ( ( z >= h_poly(l) ) .AND. ( z < h_poly(l+1) ) ) THEN
                  temp = T0_poly(l) + dTdz_poly(l) * z

                  IF ( ABS( dTdz_poly(l) ) > epsilon_dTdz ) THEN
                    h1 = 1.0_wp + dTdz_poly(l) / T0_poly(l) * ( z - h_poly(l) )
                    pres = p0_poly(l) * h1**( - grav/(Rd* dTdz_poly(l) ) )
                  ELSE
                    ! isothermal 
                    delta = grav / ( Rd * T0_poly(l) )
                    pres = p0_poly(l) * EXP( - delta * ( z - h_poly(l) ) )
                  END IF

                END IF

              END DO

              p_nh_state(jg)%prog(jt)%rho(jc,jk,jb) = pres / ( Rd * temp )
              p_nh_state(jg)%prog(jt)%exner(jc,jk,jb) = ( pres/p0ref )**( Rd/cpd )
              p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb) = temp / p_nh_state(jg)%prog(jt)%exner(jc,jk,jb)

            ENDDO
          ENDDO

          DO jk = 1, nlevp1
            p_nh_state(jg)%prog(jt)%w(1:nlen,jk,jb) = 0.0_wp

            ! copy w to reference state vector (needed by Rayleigh damping mechanism)
            p_nh_state(jg)%ref%w_ref(je,jk,jb) = &
              p_nh_state(jg)%prog(jt)%w(je,jk,jb)
          ENDDO

        ENDDO
      ENDDO

      ! edges:
      nblks_e   = p_patch(jg)%nblks_e
      npromz_e  = p_patch(jg)%npromz_e

      DO jt = 1, ntl 
        ! normal wind
        DO jb = 1, nblks_e
          IF (jb /= nblks_e) THEN
            nlen = nproma
          ELSE
            nlen = npromz_e
          ENDIF
          DO jk = 1, nlev
            DO je = 1, nlen
              p_nh_state(jg)%prog(jt)%vn(je,jk,jb) = nh_u0 &
                                !(p_nh_state(jg)%metrics%geopot(1,jk,1)/grav/1000.0_wp+5.0_wp)& !shear
                *p_patch(jg)%edges%primal_normal(je,jb)%v1

              ! copy vn to reference state vector (needed by Rayleigh damping mechanism)
              p_nh_state(jg)%ref%vn_ref(je,jk,jb)  &
                = p_nh_state(jg)%prog(jt)%vn(je,jk,jb)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      jc = 1
      jb = 1
      jt = 1

      IF ( get_my_mpi_work_id() == get_glob_proc0() ) THEN
        WRITE(*,*) "control output 1: jg=", jg
        DO jk=1, nlev
          WRITE(*,'(I5,F10.3,5F13.5)') jk, p_nh_state(jg)%metrics%z_mc(jc,jk,jb),  &
            p_nh_state(jg)%prog(jt)%rho(jc,jk,jb),      &
            p_nh_state(jg)%prog(jt)%exner(jc,jk,jb),    &
            p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb),  &
            p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb) * p_nh_state(jg)%prog(jt)%exner(jc,jk,jb),  &
            p0ref*p_nh_state(jg)%prog(jt)%exner(jc,jk,jb)**(cpd/Rd)
        END DO
      END IF

      ! numerical hydrostatic balancing is still necessary!
      IF ( l_hydro_adjust ) THEN
        CALL hydro_adjust(p_patch(jg), p_nh_state(jg)%metrics,    &
          p_nh_state(jg)%prog(nnow(jg))%rho,      &
          p_nh_state(jg)%prog(nnow(jg))%exner,    &
          p_nh_state(jg)%prog(nnow(jg))%theta_v )
      END IF

      IF ( get_my_mpi_work_id() == get_glob_proc0() ) THEN
        WRITE(*,*) "control output 2: jg=", jg
        DO jk=1, nlev
          WRITE(*,'(I5,F10.3,5F13.5)') jk, p_nh_state(jg)%metrics%z_mc(jc,jk,jb),  &
            p_nh_state(jg)%prog(jt)%rho(jc,jk,jb),     &
            p_nh_state(jg)%prog(jt)%exner(jc,jk,jb),   &
            p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb), &
            p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb) * p_nh_state(jg)%prog(jt)%exner(jc,jk,jb),  &
            p0ref*p_nh_state(jg)%prog(jt)%exner(jc,jk,jb)**(cpd/Rd)
        END DO
      END IF

      CALL diagnose_pres_temp ( p_nh_state(jg)%metrics,  p_nh_state(jg)%prog(nnow(jg)),  &
        p_nh_state(jg)%prog(nnow(jg)), p_nh_state(jg)%diag,      &
        p_patch(jg), opt_calc_pres=.TRUE., opt_calc_temp=.TRUE.)

      IF ( get_my_mpi_work_id() == get_glob_proc0() ) THEN
        WRITE(*,*) "control output 3: jg=", jg
        DO jk=1, nlev
          WRITE(*,'(I5,F10.3,5F13.5)') jk, p_nh_state(jg)%metrics%z_mc(jc,jk,jb),  &
            p_nh_state(jg)%prog(jt)%rho(jc,jk,jb),     &
            p_nh_state(jg)%prog(jt)%exner(jc,jk,jb),   &
            p_nh_state(jg)%prog(jt)%theta_v(jc,jk,jb), &
            p_nh_state(jg)%diag%temp(jc,jk,jb),        & 
            p_nh_state(jg)%diag%pres(jc,jk,jb)
        END DO
      END IF

      ! CALL sync_patch_array(SYNC_C, ptr_patch, ptr_nh_prog%w)    ! ????

    ENDDO  ! do jg=...

    DEALLOCATE( p0_poly )


  END SUBROUTINE piecewise_polytropic_atm

 
END MODULE mo_nh_testcases
