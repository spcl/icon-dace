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

! Subroutine to initialized the Jablonowski Williansom test case for the NH-Core

MODULE mo_nh_jabw_exp
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2008
!
!-------------------------------------------------------------------------
!
!
!

   USE mo_kind,                 ONLY: wp
   USE mo_physical_constants,   ONLY: rd, rd_o_cpd, p0ref, grav, tmelt,  &
     &                                cvd_o_rd, rv
   USE mo_model_domain,         ONLY: t_patch
   USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
   USE mo_run_config,           ONLY: iqv, ntracer
   USE mo_impl_constants,       ONLY: MAX_CHAR_LENGTH
   USE mo_parallel_config,      ONLY: nproma, p_test_run
   USE mo_satad,                ONLY: sat_pres_water, &  !! saturation vapor pressure w.r.t. water
     &                                sat_pres_ice       !! saturation vapor pressure w.r.t. ice
   USE mo_exception,            ONLY: finish, message_text
   USE mo_ncar_testcases,       ONLY: tracer_q1_q2, tracer_q3, &
     &                                init_ncar_testcases_domain
   USE mo_math_constants,       ONLY: pi, pi_2
   USE mo_intp_data_strc,       ONLY: t_int_state
   USE mo_intp,                 ONLY: cells2edges_scalar
   USE mo_loopindices,          ONLY: get_indices_e
   USE mo_nh_diagnose_pres_temp,ONLY: diagnose_pres_temp
   USE mo_extpar_config,        ONLY: itopo
   USE mo_sync,                 ONLY: global_sum_array, sync_patch_array, SYNC_C, SYNC_E
   USE mo_nh_init_utils,        ONLY: init_w
   USE mo_hydro_adjust,         ONLY: hydro_adjust
   USE mo_grid_config,          ONLY: grid_sphere_radius, grid_angular_velocity

   IMPLICIT NONE



   PUBLIC :: init_nh_topo_jabw,init_nh_state_prog_jabw, & 
           & init_passive_tracers_nh_jabw, &
           & init_nh_inwp_tracers

  
   PRIVATE
   REAL(wp), PARAMETER :: cpd_o_rd  = 1._wp / rd_o_cpd 



! !DEFINED PARAMETERS for jablonowski williamson:  
  REAL(wp), PARAMETER :: eta0  = 0.252_wp ! 
  REAL(wp), PARAMETER :: etat  = 0.2_wp   ! tropopause
  REAL(wp), PARAMETER :: gamma = 0.005_wp ! temperature elapse rate (K/m)
  REAL(wp), PARAMETER :: dtemp = 4.8e5_wp ! empirical temperature difference (K)

  REAL(wp), PARAMETER :: lonC  = pi/9._wp ! longitude of the perturb. centre 
  REAL(wp), PARAMETER :: latC  = 2._wp*lonC !latitude of the perturb. centre

!--------------------------------------------------------------------

   CONTAINS
!-------------------------------------------------------------------------
!
  !>
  !! Initialization of topography for the nh standard jabw test case 
  !!
  SUBROUTINE init_nh_topo_jabw( ptr_patch, topo_c, nblks_c, npromz_c, &
                             &  u0, opt_m_height, opt_m_half_width )

    TYPE(t_patch), TARGET,INTENT(INOUT) :: &  !< patch on which computation is performed
      &  ptr_patch

    REAL(wp), INTENT (IN)           :: u0
    REAL(wp), INTENT (IN), OPTIONAL :: opt_m_height, opt_m_half_width
    INTEGER, INTENT (IN) ::  nblks_c, npromz_c
    REAL(wp),  INTENT(INOUT) :: topo_c    (nproma,nblks_c)
    ! local variables

  INTEGER        :: jc, jb, nlen
  REAL(wp)       :: z_lon, z_lat
!!$  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
!!$    &  routine = '(mo_nh_jabw_exp) init_nh_topo_jabw:'
 
  REAL(wp)      :: zsiny, zcosy,tmp1,tmp2,tmp3
  REAL(wp)      :: z_fac1, z_fac2    
  LOGICAL       :: lmount
  REAL(wp)      :: mount_height, mount_half_width

!--------------------------------------------------------------------

     IF (PRESENT(opt_m_height) .AND. PRESENT (opt_m_half_width)) THEN
       lmount = .TRUE.
       mount_height = opt_m_height
       mount_half_width = opt_m_half_width
     ELSE
       lmount = .FALSE.
     END IF

 
      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen =  npromz_c
        ENDIF
        DO jc = 1, nlen
          z_lat   = ptr_patch%cells%center(jc,jb)%lat
          zsiny = SIN(z_lat)
          zcosy = COS(z_lat)
          tmp1  = u0*COS((1._wp-eta0)*pi_2)**1.5_wp
          tmp2  = (-2.0_wp*zsiny**6 * (zcosy*zcosy+1.0_wp/3.0_wp) + &
                  1.0_wp/6.3_wp ) *tmp1
          tmp3  = ( 1.6_wp*zcosy*zcosy*zcosy * (zsiny*zsiny+2.0_wp/3.0_wp)  &
                   - 0.5_wp*pi_2 )*grid_sphere_radius*grid_angular_velocity
          IF ( itopo==0 ) topo_c(jc,jb) = tmp1*(tmp2+tmp3)/grav
          IF (itopo==0 .AND. lmount ) THEN
            z_lon = ptr_patch%cells%center(jc,jb)%lon
            z_fac1= SIN(latC)*SIN(z_lat)+COS(latC)*COS(z_lat)*COS(z_lon-lonC) 
            z_fac2= grid_sphere_radius*ACOS(z_fac1)/mount_half_width
            topo_c(jc,jb) = topo_c(jc,jb) &
                          & + mount_height*EXP(-z_fac2**2)
          ENDIF 
        ENDDO
      ENDDO

  END SUBROUTINE init_nh_topo_jabw
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!
  !>
  !! Initialization of prognostic state vector for the nh standard jabw test case 
  !!  without moisture
  !!
  SUBROUTINE init_nh_state_prog_jabw( ptr_patch, ptr_nh_prog, ptr_nh_diag, &
    &                                p_metrics, p_int, p_sfc, up, u0, temp0 )

    TYPE(t_patch), TARGET,INTENT(INOUT) :: &  !< patch on which computation is performed
      &  ptr_patch

    TYPE(t_nh_prog), INTENT(INOUT)      :: &  !< prognostic state vector
      &  ptr_nh_prog

    TYPE(t_nh_diag), INTENT(INOUT)      :: &  !< diagnostic state vector
      &  ptr_nh_diag


    TYPE(t_nh_metrics), INTENT(IN)      :: p_metrics !< NH metrics state
    TYPE(t_int_state), INTENT(IN)       :: p_int

    REAL(wp), INTENT(IN)                :: p_sfc   !surface pressure, 1.e5 Pa in the standard jabw 
    REAL(wp), INTENT(IN)                :: up
    REAL(wp), INTENT(IN)                :: u0
    REAL(wp), INTENT(IN)                :: temp0



    INTEGER        ::  je, jc, jb, jk,   jn, &
                    nlen, nblks_e, nblks_c, npromz_c
    INTEGER        :: i_startidx, i_endidx, i_startblk
    INTEGER        :: nlev        !< number of full levels

    REAL(wp), DIMENSION(nproma) ::  &
             z_lat,z_siny,z_cosy, z_fac1, z_fac2, zeta_old, zcoszetav, &
             zsinzetav, z_tavg, z_favg, z_geopot, z_temp, z_fun, z_fund, &
             zeta, zu, zv, z_lon, z_exp
    REAL(wp), ALLOCATABLE :: zeta_v(:,:,:)
    REAL(wp), ALLOCATABLE :: zeta_v_e(:,:,:)
    REAL(wp)              :: ps_o_p0ref

!--------------------------------------------------------------------
!
    ! number of vertical levels
    nlev   = ptr_patch%nlev

    ALLOCATE (zeta_v(nproma,nlev,ptr_patch%nblks_c), &
              zeta_v_e(nproma,nlev,ptr_patch%nblks_e) )
    
    
    zeta_v    = 0._wp
    zeta_v_e  = 0._wp
    nblks_c   = ptr_patch%nblks_c
    npromz_c  = ptr_patch%npromz_c
    nblks_e   = ptr_patch%nblks_e

    ptr_nh_diag%pres_sfc(:,:) = p_sfc    !set surface pressure to the prescribed value
    ps_o_p0ref = p_sfc/p0ref



!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,jn,z_lat,z_siny,z_cosy,z_fac1,z_fac2,z_exp,zeta_old,&
!$OMP            zcoszetav,zsinzetav,z_tavg,z_favg,z_geopot,z_temp,z_fun,z_fund,zeta )
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
      DO jk = nlev, 1, -1
        DO jc = 1, nlen
          z_lat(jc) = ptr_patch%cells%center(jc,jb)%lat
          z_siny(jc) = SIN(z_lat(jc))
          z_cosy(jc) = COS(z_lat(jc))
          z_fac1(jc) = 1.0_wp/6.3_wp-2.0_wp*(z_siny(jc)**6)*(z_cosy(jc)**2+1.0_wp/3.0_wp)
          z_fac2(jc) = (8.0_wp/5.0_wp*(z_cosy(jc)**3)*(z_siny(jc)**2+2.0_wp/3.0_wp)&
                       -0.25_wp*pi)*grid_sphere_radius*grid_angular_velocity
          z_exp(jc)  = rd*gamma/grav
          zeta_old(jc) = 1.0e-7_wp
        ENDDO
        ! Newton iteration to determine zeta
        DO jn = 1, 100
          DO jc = 1, nlen
            zeta_v(jc,jk,jb) = (zeta_old(jc) - eta0)*pi_2
            zcoszetav(jc)= COS(zeta_v(jc,jk,jb))
            zsinzetav(jc)= SIN(zeta_v(jc,jk,jb))
            z_tavg(jc)   = temp0*(zeta_old(jc)**z_exp(jc))
            z_favg(jc)   = temp0*grav/gamma*(1.0_wp-zeta_old(jc)**z_exp(jc))
            IF (zeta_old(jc) < etat ) THEN
               z_tavg(jc) = z_tavg(jc)+dtemp*((etat-zeta_old(jc))**5)
               z_favg(jc) = z_favg(jc)-rd*dtemp*(                           &
                        (log(zeta_old(jc)/etat)+137.0_wp/60.0_wp)*(etat**5) &
                        -5.0_wp*(etat**4)*zeta_old(jc)                      &
                        +5.0_wp*(etat**3)*(zeta_old(jc)**2)                 &
                        -10.0_wp/3.0_wp*(etat**2)*(zeta_old(jc)**3)         &
                        +1.25_wp*etat*(zeta_old(jc)**4)-0.2_wp*(zeta_old(jc)**5))
            ENDIF
            z_geopot(jc) = z_favg(jc)+u0*(zcoszetav(jc)**1.5_wp)*&
                          (z_fac1(jc)*u0*(zcoszetav(jc)**1.5_wp)+z_fac2(jc)) 
            z_temp(jc)   = z_tavg(jc)+0.75_wp*zeta_old(jc)*pi*u0/rd*zsinzetav(jc)*&
                       SQRT(zcoszetav(jc))*(2.0_wp*u0*z_fac1(jc)*(zcoszetav(jc)**1.5_wp) &
                       + z_fac2(jc))
            z_fun(jc)    = z_geopot (jc)- p_metrics%geopot(jc,jk,jb)
            z_fund(jc)   = -rd/zeta_old(jc)*z_temp(jc)
            zeta(jc) = zeta_old(jc) - z_fun(jc)/z_fund(jc)
            zeta_old(jc) = zeta(jc)
          ENDDO ! jc
        ENDDO !jn
        ! Final update for zeta_v
        DO jc = 1, nlen
          zeta_v(jc,jk,jb) = (zeta_old(jc) - eta0)*pi_2
        ENDDO
        ! Use analytic expressions at all model level
          DO jc = 1, nlen
            ptr_nh_prog%exner(jc,jk,jb) = (zeta_old(jc)*ps_o_p0ref)**rd_o_cpd
            ptr_nh_prog%theta_v(jc,jk,jb) = z_temp(jc) & 
            &        /ptr_nh_prog%exner(jc,jk,jb)
            ptr_nh_prog%rho(jc,jk,jb) = &
            &        ptr_nh_prog%exner(jc,jk,jb)**cvd_o_rd*p0ref/rd &
            &        /ptr_nh_prog%theta_v(jc,jk,jb)
            !initialize diagnose pres and temp variables
            ptr_nh_diag%pres(jc,jk,jb) = p0ref*ptr_nh_prog%exner(jc,jk,jb)**(cpd_o_rd) 
            ptr_nh_diag%temp(jc,jk,jb) = z_temp(jc)  

          ENDDO !jc
      ENDDO !jk
    ENDDO !jb
!$OMP END DO
!$OMP END PARALLEL

    IF (p_test_run) zeta_v_e = 0._wp
    CALL cells2edges_scalar(zeta_v,ptr_patch,p_int%c_lin_e,zeta_v_e)
    CALL sync_patch_array(SYNC_E,ptr_patch,zeta_v_e)

    i_startblk = ptr_patch%edges%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,z_lat,z_lon,zu,z_fac1,z_fac2,zv)
    DO jb = i_startblk, nblks_e

      CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, i_startidx, i_endidx, 2)

        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            z_lat(je) = ptr_patch%edges%center(je,jb)%lat
            z_lon(je) = ptr_patch%edges%center(je,jb)%lon
            zu(je)    = u0*(COS(zeta_v_e(je,jk,jb))**1.5_wp)*(SIN(2.0_wp*z_lat(je))**2)
            IF ( up .GT. 1.e-20_wp ) THEN
             z_fac1(je)= SIN(latC)*SIN(z_lat(je))+COS(latC)*COS(z_lat(je))*COS(z_lon(je)-lonC) 
             z_fac2(je)  = 10._wp*ACOS(z_fac1(je))
             zu(je) = zu(je) + up* EXP(-z_fac2(je)**2)
            END IF
            zv(je) = 0._wp
            ptr_nh_prog%vn(je,jk,jb) = &
                       zu(je) * ptr_patch%edges%primal_normal(je,jb)%v1   &
                   & + zv(je) * ptr_patch%edges%primal_normal(je,jb)%v2
          ENDDO
        ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL 
   DEALLOCATE(zeta_v, zeta_v_e)

  ! initialize vertical velocity
   CALL init_w(ptr_patch, p_int, ptr_nh_prog%vn, p_metrics%z_ifc, ptr_nh_prog%w)
   CALL sync_patch_array(SYNC_C, ptr_patch, ptr_nh_prog%w)

   CALL hydro_adjust ( ptr_patch, p_metrics, ptr_nh_prog%rho,  &
                     & ptr_nh_prog%exner, ptr_nh_prog%theta_v  )

  END SUBROUTINE init_nh_state_prog_jabw
!-------------------------------------------------------------------------
! 
!>
  !! Defines passive traces for the nh jabw test case
  SUBROUTINE init_passive_tracers_nh_jabw( ptr_patch, ptr_nh_prog,           &
                                     &     rotate_axis_deg, tracer_inidist_list, p_sfc)  

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'init_passive_tracers_nh_jabw'

    TYPE(t_patch), TARGET,INTENT(INOUT) :: &  !< patch on which computation is performed
      &  ptr_patch

    TYPE(t_nh_prog), INTENT(INOUT)      :: &  !< prognostic state vector
      &  ptr_nh_prog
    INTEGER                             :: &  !< tracer initial distributions
      &  tracer_inidist_list(:)

    REAL(wp), ALLOCATABLE               :: zeta_v(:,:,:)

    REAL(wp)                            :: rotate_axis_deg
    REAL(wp), INTENT(IN)                :: p_sfc   !surface pressure, 1e5 Pa in the standard jabw 


   ! local variables

   INTEGER                              :: nlev
   INTEGER                              :: jc, jb, jk, jjt
   INTEGER                              :: nblks_c, npromz_c,&
                                         &  nlen
   REAL(wp)                             :: zlat, zlon 
   REAL(wp)                             :: ps_o_p0ref

!--------------------------------------------------------------------
!
    CALL init_ncar_testcases_domain()
   ! number of vertical levels
    nlev   = ptr_patch%nlev
    nblks_c   = ptr_patch%nblks_c
    npromz_c  = ptr_patch%npromz_c
    nlev      = ptr_patch%nlev

   ps_o_p0ref = p_sfc/p0ref
   ALLOCATE (zeta_v(nproma,nlev,ptr_patch%nblks_c))
   !calculate zeta_v from exner function
   DO jb = 1, nblks_c
     IF (jb /= nblks_c) THEN 
      nlen = nproma
      ELSE
      nlen = npromz_c
     ENDIF
     DO jk = 1, nlev
      zeta_v(1:nlen, jk, jb) = ( (ptr_nh_prog%exner(1:nlen, jk, jb))**cpd_o_rd &
                                &                      /ps_o_p0ref - eta0 )*pi_2
     END DO
   END DO


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jjt,jc,nlen,zlat,zlon)
   DO jb = 1, nblks_c
     IF (jb /= nblks_c) THEN 
      nlen = nproma
      ELSE
      nlen = npromz_c
     ENDIF
     DO jk = 1, nlev
      DO jjt = 1, ntracer
        SELECT CASE(tracer_inidist_list(jjt))
  
        CASE(1)
  
                    DO jc = 1, nlen
                      zlat = ptr_patch%cells%center(jc,jb)%lat
                      zlon = ptr_patch%cells%center(jc,jb)%lon
                      ptr_nh_prog%tracer(jc,jk,jb,jjt) =  &
                      tracer_q1_q2(zlon, zlat, zeta_v(jc,jk,jb), rotate_axis_deg, 0.6_wp)
                    ENDDO ! cell loop
  
        CASE(2)
  
                    DO jc =1, nlen
                      zlat = ptr_patch%cells%center(jc,jb)%lat
                      zlon = ptr_patch%cells%center(jc,jb)%lon
                      ptr_nh_prog%tracer(jc,jk,jb,jjt) =  &
                        tracer_q1_q2(zlon, zlat, zeta_v(jc,jk,jb), rotate_axis_deg, 1.0_wp)
                    ENDDO ! cell loop
  
        CASE(3)
  
                    DO jc =1, nlen
                      zlat = ptr_patch%cells%center(jc,jb)%lat
                      zlon = ptr_patch%cells%center(jc,jb)%lon
                      ptr_nh_prog%tracer(jc,jk,jb,jjt) =  &
                        tracer_q3(zlon, zlat, rotate_axis_deg)
                    ENDDO ! cell loop
  
        CASE(4)
                    ptr_nh_prog%tracer(:,jk,jb,jjt) = 1._wp

        CASE default
                    WRITE(message_text,'(a,i3,a)') 'Undefined tracer initial distribution ', &
                      &  tracer_inidist_list(jjt),'; must be <=4'
                    CALL finish(routine,TRIM(message_text))
        END SELECT
       END DO
      END DO
     END DO  
 
!$OMP END DO
!$OMP END PARALLEL  
   DEALLOCATE(zeta_v)

  END SUBROUTINE init_passive_tracers_nh_jabw
!-------------------------------------------------------------------------
!
  !>
  !! Initialization of tracers for the jabw test in case of inwp  
  !!
  SUBROUTINE init_nh_inwp_tracers( ptr_patch, ptr_nh_prog, ptr_nh_diag, &
    &                              p_metrics, rh_at_1000hpa, qv_max,    &
    &                              l_rediag, opt_global_moist )

    TYPE(t_patch), TARGET,INTENT(INOUT) :: &  !< patch on which computation is performed
      &  ptr_patch

    TYPE(t_nh_prog), INTENT(INOUT)      :: &  !< prognostic state vector
      &  ptr_nh_prog

    TYPE(t_nh_diag), INTENT(INOUT)      :: &  !< diagnostic state vector
      &  ptr_nh_diag
    TYPE(t_nh_metrics), INTENT(IN)      :: p_metrics !< NH metrics state

    REAL(wp), INTENT (IN)               :: rh_at_1000hpa, qv_max
    REAL(wp), INTENT (IN), OPTIONAL     :: opt_global_moist     !global moisture content in kg/m**2
    LOGICAL,  INTENT (IN)               :: l_rediag    

   ! local variables

    INTEGER                             :: nblks_c,  npromz_c, nlen,  &
                                           nlev 
    INTEGER                             :: jb,jc, jk, jjt, ji
    REAL(wp)                            :: zsqv, z_help, z_1_o_rh, zrhf, z_moist,          &
                                           z_aux(nproma,ptr_patch%nlev,ptr_patch%nblks_c), &
                                           z_area(nproma,ptr_patch%nblks_c)
    INTEGER                             :: niter 
    LOGICAL                             :: l_global_moist
!--------------------------------------------------------------------

    nlev      = ptr_patch%nlev
    nblks_c   = ptr_patch%nblks_c
    npromz_c  = ptr_patch%npromz_c

 IF (PRESENT(opt_global_moist)) THEN
    l_global_moist=.TRUE.
   ELSE
    l_global_moist=.FALSE.
 ENDIF
 IF (l_rediag) THEN
   niter = 10
 ELSE
   niter = 1
 END IF

    ! Do some iterations to come closer to the moisture and temperature/pressure 
    DO ji = 1, niter
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jjt,jc,nlen,zrhf,z_1_o_rh,z_help,zsqv)
        DO jb = 1, nblks_c
          IF (jb /= nblks_c) THEN 
            nlen = nproma
          ELSE
            nlen = npromz_c
          ENDIF

            DO jk = 1, nlev
              DO jjt = 1, ntracer
                  IF(jjt == iqv ) THEN
                    DO jc =1, nlen

!                      !KF linear decreasing RH with height like in Hui's testcase
                      zrhf     = rh_at_1000hpa-0.5_wp+ptr_nh_diag%pres(jc,jk,jb)/200000._wp
                      zrhf     = MAX (zrhf,0.0_wp)
                      z_1_o_rh = 1._wp/(zrhf+1.e-6_wp)
                      ! to avoid water vapor pressure > total pressure:
                      IF( ptr_nh_diag%temp(jc,jk,jb) <= tmelt) THEN
                        z_help = MIN ( sat_pres_ice( MAX(180._wp,ptr_nh_diag%temp(jc,jk,jb)) ), &
                         & ptr_nh_diag%pres(jc,jk,jb) * z_1_o_rh )
                      ELSE
                        z_help = MIN ( sat_pres_water( ptr_nh_diag%temp(jc,jk,jb) ), &
                         & ptr_nh_diag%pres(jc,jk,jb) * z_1_o_rh )
                      ENDIF
                      ! saturation qv calculated as in mo_satad's qsat_rho
                      zsqv = z_help / ( ptr_nh_prog%rho(jc,jk,jb) * rv  &
                         &      * ptr_nh_diag%temp(jc,jk,jb) )
                      ptr_nh_prog%tracer(jc,jk,jb,jjt) = MIN ( zsqv,  zrhf*zsqv )

                      IF (  ptr_nh_diag%pres(jc,jk,jb) <= 10000._wp) &
                         &  ptr_nh_prog%tracer(jc,jk,jb,jjt)      &
                         &  = MIN ( 5.e-6_wp, ptr_nh_prog%tracer(jc,jk,jb,jjt) )
                    
                      ! Limit QV in the tropics                       
                      ptr_nh_prog%tracer(jc,jk,jb,jjt) = &
                      &   MIN(qv_max,ptr_nh_prog%tracer(jc,jk,jb,jjt))
  
                    END DO
  
                  ELSE !other tracers than water vapor zero at start
                    ptr_nh_prog%tracer(:,jk,jb,jjt) = 0._wp
                  ENDIF ! tracer              
             ENDDO ! tracer loop
            ENDDO ! vertical level loop
        ENDDO ! block loop
!$OMP END DO
!$OMP END PARALLEL


      IF(l_rediag) THEN
        CALL diagnose_pres_temp ( p_metrics, ptr_nh_prog,     &
          &                       ptr_nh_prog, ptr_nh_diag,   &
          &                       ptr_patch,                  &
          &                       opt_calc_temp=.TRUE.,       &
          &                       opt_calc_pres=.TRUE.        )
      ENDIF

    ENDDO ! ji

    IF (l_global_moist) THEN    ! set the global moisture content to the prescribed value

      z_aux(:,:,:) = 0._wp
      z_area(:,:)  = 0._wp

      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF

        DO jc = 1, nlen
          IF (ptr_patch%cells%decomp_info%owner_mask(jc,jb)) & 
            z_area(jc,jb) = ptr_patch%cells%area(jc,jb)
        ENDDO !jc

        DO jk = 1, nlev
          DO jc = 1, nlen
            IF (ptr_patch%cells%decomp_info%owner_mask(jc,jb)) & 
              z_aux(jc,jk,jb) = p_metrics%ddqz_z_full(jc,jk,jb) * ptr_nh_prog%rho(jc,jk,jb) &
                              * ptr_nh_prog%tracer(jc,jk,jb,iqv) *ptr_patch%cells%area(jc,jb)
          ENDDO !jc
        ENDDO !jk

      ENDDO !jb

      z_moist = global_sum_array(z_aux) / global_sum_array(z_area)

      IF (z_moist .GT. 1.e-25_wp) THEN
        ptr_nh_prog%tracer(:,:,:,iqv) = ptr_nh_prog%tracer(:,:,:,iqv) &
          &                           * opt_global_moist / z_moist
      END IF

      IF(l_rediag) THEN
        CALL diagnose_pres_temp ( p_metrics, ptr_nh_prog,     &
            &                     ptr_nh_prog, ptr_nh_diag,   &
            &                     ptr_patch,                  &
            &                     opt_calc_temp=.TRUE.,       &
            &                     opt_calc_pres=.TRUE.        )
      ENDIF
    END IF 


   END SUBROUTINE init_nh_inwp_tracers
!--------------------------------------------------------------------
  END MODULE mo_nh_jabw_exp
