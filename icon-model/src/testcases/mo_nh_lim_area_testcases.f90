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

! Subroutine to initialized several test cases
! for the NH-Core in limited area mode
!
! Literature
! -M. L. Weisman and J. B. Klemp, 1982
!  The Dependence of Numerically Simulated Convective Storms on
!  Vertical Wind Shear and Buoyancy.
!  Monthly Weather Review, 110,504-520
! - COSMO documentation

  MODULE mo_nh_lim_area_testcases
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2008
!
!-------------------------------------------------------------------------
!

   USE mo_kind,                ONLY: wp
   USE mo_physical_constants,  ONLY: rd_o_cpd, p0ref, grav, tmelt,  &
                                   & cvd_o_rd, cpd ,     &
                                     vtmpc1 , rd            
                                     
   USE mo_math_constants,      ONLY: pi, deg2rad
   USE mo_model_domain,        ONLY: t_patch
   USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
   USE mo_run_config,          ONLY: iforcing, iqv,msg_level 
   USE mo_impl_constants,      ONLY: inwp
   USE mo_parallel_config,     ONLY: nproma
   USE mo_satad,               ONLY:  sat_pres_water, &  !! saturation vapor pressure w.r.t. water
            &                         sat_pres_ice,   &  !! saturation vapor pressure w.r.t. ice
            &                         spec_humi          !! Specific humidity
   USE mo_exception,           ONLY: message, finish, message_text
   USE mo_intp_data_strc,      ONLY: t_int_state
   USE mo_intp,                ONLY: cells2edges_scalar
   USE mo_loopindices,         ONLY: get_indices_e
   USE mo_nh_diagnose_pres_temp,ONLY: diagnose_pres_temp
   USE mo_extpar_config,        ONLY: itopo
   USE mo_sync,                 ONLY: sync_patch_array,  SYNC_C, SYNC_E
   USE mo_nh_init_utils,        ONLY: init_w, convert_thdvars
   USE mo_hydro_adjust,        ONLY: hydro_adjust, hydro_adjust_iterative
   USE mo_grid_config,         ONLY: grid_sphere_radius

   IMPLICIT NONE

   PUBLIC  :: init_nh_atmo_ana_nconstlayers, init_nh_anaprof_uv, &
              init_nh_topo_ana, init_nh_atmo_ana_poly

   PRIVATE

   REAL(wp), PARAMETER :: cpd_o_rd = 1._wp / rd_o_cpd
   REAL(wp), PARAMETER :: grav_o_cpd = grav / cpd
   REAL(wp), PARAMETER :: grav_o_rd = grav / rd

! !DEFINED namelist variables for the analytical atmosphere profile
  INTEGER, PUBLIC       :: itype_atmo_ana    ! kind of atmosphere profile
                                             ! 1 piecewise N const layers
                                             ! 2 piecewise polytropic layers
! !DEFINED PARAMETERS for the piecewise const. Brunt-Vaisala-freq (N) 
!  layers atmosphere:
   INTEGER, PARAMETER, PUBLIC  :: max_nlayers_nconst = 10
! !DEFINED PARAMETERS for the piecewise const. vertical T-gradient 
!  layers atmosphere:
   INTEGER, PARAMETER, PUBLIC  :: max_nlayers_poly = 10
! !DEFINED PARAMETERS for the analytical wind profiles 
!  layers atmosphere:
   INTEGER, PARAMETER, PUBLIC  :: max_nlayers_linwind = 10

! !DEFINED namelist variables for the piecewise const. Brunt-Vaisala-freq (N) 
!  layers atmosphere:
   INTEGER, PUBLIC      :: nlayers_nconst      ! Number of desired layers with a constant N
   REAL(wp), PUBLIC     :: p_base_nconst       ! pressure at h_nconst(1)
   REAL(wp) , PUBLIC    :: theta0_base_nconst  ! potential temperature at h_nconst(1)
   REAL(wp), PUBLIC     :: h_nconst(max_nlayers_nconst) ! base height for each layer in m
   REAL(wp), PUBLIC     :: N_nconst(max_nlayers_nconst) ! N in 1/s for each layer
   REAL(wp), PUBLIC     :: rh_nconst(max_nlayers_nconst) ! base relative humidity for each layer
   REAL(wp), PUBLIC     :: rhgr_nconst(max_nlayers_nconst) ! gradient of relative humidity for 
                                                           ! each layer in 1/m, positive for 
                                                           ! decreasing rel hum with height
! !DEFINED namelist variables for the piecewise const. vertical T-gradient (G) 
!  layers atmosphere:
   INTEGER, PUBLIC      :: nlayers_poly        ! Number of desired layers with a constant G
   REAL(wp), PUBLIC     :: p_base_poly       ! pressure at h_poly(1)
   REAL(wp), PUBLIC     :: h_poly(max_nlayers_poly) ! base height for each layer in m
   REAL(wp), PUBLIC     :: t_poly(max_nlayers_poly) ! base T in K for each layer
   REAL(wp), PUBLIC     :: tgr_poly(max_nlayers_poly) ! G in K/m for each layer, positive for 
                                                      ! decreasing temperature with height
   REAL(wp), PUBLIC     :: rh_poly(max_nlayers_poly) ! base relative humidity for each layer
   REAL(wp), PUBLIC     :: rhgr_poly(max_nlayers_poly) ! gradient of relative humidity for 
                                                           ! each layer in 1/m, positive for 
                                                           ! decreasing rel hum with height
! !DEFINED namelist variables for the analytical wind profile
  INTEGER, PUBLIC       :: itype_anaprof_uv    ! kind of wind profile
  ! For  itype_anaprof_uv == 1 (arbitrary number of constant gradient U(z) layers ):
  INTEGER, PUBLIC       :: nlayers_linwind     ! Number of layers
  REAL(wp), PUBLIC      :: h_linwind(max_nlayers_linwind) ! base height for each layer in m
  REAL(wp), PUBLIC      :: u_linwind(max_nlayers_linwind) ! U at base layer in m/s
  REAL(wp), PUBLIC      :: ugr_linwind(max_nlayers_linwind) ! gradient of U for each layer,
                                                            ! positive for increasing windspeed
                                                            ! with height, in 1/s
  ! For  itype_anaprof_uv == 2/3, constant U/V 
  REAl(wp), PUBLIC      :: vel_const 
                       ! 
! !DEFINED namelist variables for the analytical topography
  INTEGER,  PUBLIC      :: itype_topo_ana
! !Defined namelist parameters for the schaer mountain (itype_topo_ana=1)
  REAL(wp), PUBLIC      :: schaer_h0
  REAL(wp), PUBLIC      :: schaer_a
  REAL(wp), PUBLIC      :: schaer_lambda
! Defined namelist parameter for all the 2D mountains
  REAL(wp), PUBLIC      :: halfwidth_2dm 
! !Defined namelist parameters for the center of the mountain
  REAL(wp), PUBLIC      :: mount_lonc_deg, mount_latc_deg
! !Defined namelist parameters for the height and width of the mountain
  REAL(wp), PUBLIC      :: m_height, m_width_x, m_width_y
     

!--------------------------------------------------------------------

  CONTAINS
!--------------------------------------------------------------------
!
  !>
  !! Initialization of prognostic state for an 
  !! atmosphere with arbitrary number of polytropic layers, specified 
  !! by a constant vertical T-gradient (G lapse rate)
  !!
  !! It is assumed  a moist subsaturated atmosphere 
  !!
  SUBROUTINE init_nh_atmo_ana_poly( ptr_patch, ptr_nh_prog, ptr_nh_diag, &
    &                                p_metrics, l_hydro_adjust )

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_nh_lim_area_testcases:init_nh_atmo_ana_poly'


    TYPE(t_patch), TARGET, INTENT(INOUT) :: &  !< patch on which computation is performed
      &  ptr_patch

    TYPE(t_nh_prog), INTENT(INOUT)      :: &  !< prognostic state vector
      &  ptr_nh_prog

    TYPE(t_nh_diag), INTENT(INOUT)      :: &  !< diagnostic state vector
      &  ptr_nh_diag


    TYPE(t_nh_metrics), INTENT(IN)      :: p_metrics !< NH metrics state
    LOGICAL, INTENT(IN)                 :: l_hydro_adjust !if .TRUE. hydrostatically balanced 
                                                         ! initial condition



    INTEGER        ::  jc, jb, jk,    &
                       nlen, nblks_c, npromz_c, &
                       jl, jg, jsubl
    INTEGER        :: nlev        !< number of full levels

    REAL(wp), ALLOCATABLE :: pres(:,:,:), temp(:,:,:)
    REAL(wp), ALLOCATABLE :: relhum(:,:,:), tempv(:,:,:)
    REAL(wp), ALLOCATABLE :: z_qv(:,:,:)
    INTEGER,  ALLOCATABLE :: jglayer(:,:,:)

    REAL(wp), DIMENSION(nlayers_poly) :: pres_poly, qv_poly

    ! sublayers within a layer to better stimate pressure at the base of the layers
    INTEGER, PARAMETER  :: nsubl = 10  !(nsubl-1 is the number of sublayers used)
    REAL(wp), DIMENSION(nsubl) :: h_subl, t_subl, pres_subl, rh_subl, qv_subl
  
    REAL(wp)  :: z_h, zz_top
    REAL(wp)  :: z_h_kp1
!--------------------------------------------------------------------

    ! number of vertical levels
    nlev   = ptr_patch%nlev
    ALLOCATE (pres(nproma,nlev,ptr_patch%nblks_c), &
              temp(nproma,nlev,ptr_patch%nblks_c), &
              tempv(nproma,nlev,ptr_patch%nblks_c), &
              relhum(nproma,nlev,ptr_patch%nblks_c), &
              z_qv(nproma,nlev,ptr_patch%nblks_c))
    ALLOCATE ( jglayer(nproma,nlev,ptr_patch%nblks_c) )
    nblks_c   = ptr_patch%nblks_c
    npromz_c  = ptr_patch%npromz_c

!
! First some control of the imput namelist parameters 

   IF (nlayers_poly > max_nlayers_poly) THEN
        CALL finish('ERROR  ** nlayers_poly > max_nlayers_poly !!! **')
   ENDIF
   IF (ANY(h_poly(2:nlayers_poly) < h_poly(1:nlayers_poly-1))) THEN
        CALL finish( 'ERROR  ** h_poly (k+1) < h_poly(k) occured &
                     & (h_poly must be increasing) !!! ** ')
   ENDIF
   IF (ANY(rh_poly(1:nlayers_poly-1) < 0.0_wp) .OR. &
       ANY(rh_poly(1:nlayers_poly-1) - rhgr_poly(1:nlayers_poly-1)*&
       (h_poly(2:nlayers_poly)-h_poly(1:nlayers_poly-1)) < 0.0_wp)) THEN
              CALL message('',' WARNING  *** wrong values of rh_poly and/or rhgr_poly &
                         & lead to relhum < 0 !!! ***  It will be set to 0')
   ENDIF
   zz_top = MAXVAL(p_metrics%z_mc)
   IF (h_poly(nlayers_poly) <= zz_top) THEN
         IF (rh_poly(nlayers_poly) < 0.0_wp .OR. &
                   rh_poly(nlayers_poly) - rhgr_poly(nlayers_poly)*&
                   (zz_top-h_poly(nlayers_poly)) < 0.0_wp) THEN
           CALL message('',' WARNING  *** wrong values of rh_poly and/or  &
                     & rhgr_poly at the top lead to relhum < 0 !!! *** &
                     & It will be set to 0')
         ENDIF
   ENDIF
! condition to avoid negative pressure in a polytropic atmosphere

  DO jl=1, nlayers_poly-1
   IF (tgr_poly(jl) > 0._wp ) THEN   
    IF (h_poly(jl+1) > h_poly(jl)+t_poly(jl)/tgr_poly(jl) ) THEN
      WRITE(message_text,'(a,i3,a)') 'jl:',jl,' combination of h_poly(jl), t_poly(jl) &
                               & tgr_poly(jl) lead to negative pressure'
      CALL finish('', message_text)
    END IF
   END IF
  END DO
  IF (tgr_poly(nlayers_poly) > 0._wp ) THEN
   IF ( zz_top > h_poly(nlayers_poly)+t_poly(nlayers_poly)/tgr_poly(nlayers_poly) ) THEN
     WRITE(message_text,'(a)') ' combination of h_poly(nlayers_poly), t_poly(nlayers_poly) &
                             & tgr_poly(nlayers_poly) lead to negative pres'
     CALL finish('', message_text)
   END IF
  END IF


! Now start the initialization of the p T profile

   ! Set the values for h_poly(1)
   pres_poly(1)=p_base_poly
   qv_poly(1)= qv_rhtp(rh_poly(1), t_poly(1), pres_poly(1) )

   ! Set the values for the base of the rest of the layers
!   !$OMP PARALLEL
!   !$OMP DO PRIVATE(jl, jsubl,h_subl,t_subl,rh_subl,pres_subl,qv_subl )
  DO  jl = 2, nlayers_poly
    !consider temp and rh as from the layer below 
    ! I consider nsubl-1 sublayers in layer jl, to stimate pres_poly(jl)
    h_subl(1)=h_poly(jl-1)
    t_subl(1)=t_poly(jl-1)
    pres_subl(1)=pres_poly(jl-1)
    rh_subl(1)=rh_poly(jl-1)
    qv_subl(1)=qv_poly(jl-1)
    DO jsubl=2,nsubl
     h_subl(jsubl)= h_subl(1) + (h_poly(jl)-h_subl(1))*REAL(jsubl-1,wp)/REAL(nsubl-1,wp)
     t_subl(jsubl)=t_subl(jsubl-1)-tgr_poly(jl-1)*(h_subl(jsubl)-h_subl(jsubl-1))
     rh_subl(jsubl)=rh_subl(jsubl-1)-rhgr_poly(jl-1)*(h_subl(jsubl)-h_subl(jsubl-1))
     rh_subl(jsubl)=MAX(0._wp,MIN(rh_subl(jsubl),1.0_wp))
     IF (rh_subl(jsubl)> 0.0_wp) THEN
       !estimale pres_poly considering constant qv in the layer below
       pres_subl(jsubl)=p_poly(pres_subl(jsubl-1), tgr_poly(jl-1), t_subl(jsubl-1), &
                    h_subl(jsubl-1),h_subl(jsubl) , qv_subl(jsubl-1),qv_subl(jsubl-1)) 
       
       ! first stimation of qv
       qv_subl(jsubl)=qv_rhtp(rh_subl(jsubl), t_subl(jsubl), pres_subl(jsubl) )
       pres_subl(jsubl)=p_poly(pres_subl(jsubl-1), tgr_poly(jl-1), t_subl(jsubl-1), &
                    h_subl(jsubl-1),h_subl(jsubl) , qv_subl(jsubl-1),qv_subl(jsubl)) 
       
       qv_subl(jsubl)=qv_rhtp(rh_subl(jsubl), t_subl(jsubl), pres_subl(jsubl) )
     ELSE
       qv_subl(jsubl)=0.0_wp
     END IF 
     !last estimation of pres_subl
     pres_subl(jsubl)=p_poly(pres_subl(jsubl-1), tgr_poly(jl-1), t_subl(jsubl-1), &
                  h_subl(jsubl-1),h_subl(jsubl) , qv_subl(jsubl-1),qv_subl(jsubl)) 

    END DO !jsubl
    pres_poly(jl)=pres_subl(nsubl)
    ! now consider the real t_poly(jl) and rh_poly(jl) to calculate qv_poly(jl) (=! qv_aux)
    qv_poly(jl)=qv_rhtp(rh_poly(jl), t_poly(jl), pres_poly(jl))
   
  END DO !jl
!   !$OMP END DO
!   !$OMP END PARALLEL

  IF (msg_level >= 5) THEN ! print maximum velocities in global domain

   CALL message('','polytropic layers: layer PRES TEMP RH QV')
   DO  jl = 1, nlayers_poly
     WRITE(message_text,'(10x,i4,4e18.10)') jl, &
                             & pres_poly(jl), t_poly(jl), rh_poly(jl), qv_poly(jl)
     CALL message('', message_text)
   END DO

  ENDIF ! msg_level >= 5

! set the corresponding layer for all the model points

!$OMP PARALLEL
!$OMP DO PRIVATE(jk,jc,nlen,z_h,jl)
    DO jb = 1, nblks_c
      nlen = MERGE(nproma, npromz_c, jb /= nblks_c)

      DO jk = nlev, 1, -1
         DO jc = 1,nproma
           jglayer(jc, jk, jb) = 0
         END DO
         DO jc = 1, nlen

            z_h = p_metrics%z_mc(jc,jk,jb)
             !set the layer corresponding to this point
             IF (z_h <  h_poly(1) ) THEN
               jglayer(jc,jk,jb)=0
             ELSEIF (z_h >=  h_poly(nlayers_poly) ) THEN 
              jglayer(jc,jk,jb)= nlayers_poly
             ELSE
              DO jl=1,nlayers_poly-1
               IF (z_h >= h_poly(jl) .AND. z_h < h_poly(jl+1)  ) THEN
                jglayer(jc,jk,jb)=jl 
               END IF          
              END DO
             END IF

         END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    IF (ANY(jglayer(:,:,:) < 0) .OR. ANY(jglayer(:,:,:)>nlayers_poly )) THEN
         CALL finish ('corresponding layer has not been found for some model points')
    END IF
  
!$OMP PARALLEL 
!$OMP DO PRIVATE(jk,jc,nlen,z_h,z_h_kp1, jg)
    DO jb = 1, nblks_c
      nlen = MERGE(nproma, npromz_c, jb /= nblks_c)

      DO jk = nlev, 1, -1
         DO jc = 1, nlen
 
             z_h = p_metrics%z_mc(jc,jk,jb)
             jg=jglayer(jc,jk,jb)

             IF (jg == 0) THEN
              !consider a isothermal atmosphere below h_poly(1) with qv=qv_poly(1)
              z_qv(jc,jk,jb) = qv_poly(1)
              temp(jc,jk,jb) = t_poly(1)
              tempv(jc,jk,jb) =  temp(jc,jk,jb) * (1._wp+vtmpc1*qv_poly(1))
              pres(jc,jk,jb) = pres_poly(1) * &
                & EXP(-grav_o_rd*(z_h-h_poly(1))/tempv(jc,jk,jb))
             ELSE
              IF (ABS(tgr_poly(jg)) > 1.e-20) THEN
               temp(jc,jk,jb)=t_poly(jg)- tgr_poly(jg) * (z_h-h_poly(jg))
              ELSE
               temp(jc,jk,jb)=t_poly(jg)
              END IF
              !temp_aux=t_poly(jg)- tgr_poly(jg) * (z_h-h_poly(jg))
              IF (ABS(rhgr_poly(jg)) > 1.e-20) THEN
               relhum(jc,jk,jb)=rh_poly(jg)-rhgr_poly(jg)*(z_h-h_poly(jg))
              ELSE
               relhum(jc,jk,jb)=rh_poly(jg)
              END IF
              relhum(jc,jk,jb)=MAX(0._wp,MIN(1._wp,relhum(jc,jk,jb)))
              IF (jk < nlev .AND. jg == jglayer(jc,jk+1,jb) ) THEN
               ! in this case we integrate starting in the level bellow, instead of 
               !  starting in the base of the layer
               z_h_kp1 = p_metrics%z_mc(jc,jk+1,jb)
               IF (relhum(jc,jk,jb) > 1.e-20_wp) THEN
                pres(jc,jk,jb)=p_poly(pres(jc,jk+1,jb),tgr_poly(jg),temp(jc,jk+1,jb), &
                                z_h_kp1 , z_h, z_qv(jc,jk+1,jb),z_qv(jc,jk+1,jb))
                
                
                z_qv(jc,jk,jb) = qv_rhtp(relhum(jc,jk,jb), temp(jc,jk,jb),            &
                                         pres(jc,jk,jb) ) 
                pres(jc,jk,jb)=p_poly(pres(jc,jk+1,jb),tgr_poly(jg),temp(jc,jk+1,jb), &
                                z_h_kp1 , z_h, z_qv(jc,jk+1,jb),z_qv(jc,jk,jb))
                z_qv(jc,jk,jb) = qv_rhtp(relhum(jc,jk,jb), temp(jc,jk,jb),            &
                                         pres(jc,jk,jb) ) 
               ELSE
                z_qv(jc,jk,jb) = 0.0_wp
               END IF 
               pres(jc,jk,jb)=p_poly(pres(jc,jk+1,jb),tgr_poly(jg),temp(jc,jk+1,jb), &
                                z_h_kp1 , z_h, z_qv(jc,jk+1,jb),z_qv(jc,jk,jb))
               tempv(jc,jk,jb) =  temp(jc,jk,jb) * (1._wp+vtmpc1*z_qv(jc,jk,jb))

               
              ELSE
               !WRITE(*,*) "from the layer base"
               ! here we integrate from the base of the layer
               !
               ! here I could earth_radious-estimate the values at the base of the layer
               ! integrating from (jc,jk+1,jb) to the base of the layer 
               ! in the case jN>1
!!$                IF (jg>1  .AND. pres_poly(jg) > pres(jc,jk+1,jb)) THEN
!!$                  ! check that pres_poly(jg)<pres(jc,jk+1,jb)
!!$                  CALL finish ('base layer has larger pressure than the level&
!!$                              & below pres_poly(jg) was not well approximated')
!!$                             
!!$                END IF
               IF (relhum(jc,jk,jb) > 1.e-20_wp) THEN
                pres(jc,jk,jb)=p_poly(pres_poly(jg),tgr_poly(jg),t_poly(jg), &
                             h_poly(jg), z_h ,qv_poly(jg),qv_poly(jg))              
                z_qv(jc,jk,jb) = qv_rhtp(relhum(jc,jk,jb), temp(jc,jk,jb),   &
                                         pres(jc,jk,jb) )
                pres(jc,jk,jb)=p_poly(pres_poly(jg),tgr_poly(jg),t_poly(jg), &
                             h_poly(jg), z_h, qv_poly(jg),z_qv(jc,jk,jb))               
                z_qv(jc,jk,jb) = qv_rhtp(relhum(jc,jk,jb), temp(jc,jk,jb),   &
                                         pres(jc,jk,jb) )
               ELSE
                z_qv(jc,jk,jb) = 0.0_wp
               END IF 
               pres(jc,jk,jb)=p_poly(pres_poly(jg),tgr_poly(jg),t_poly(jg), &
                             h_poly(jg), z_h, qv_poly(jg),z_qv(jc,jk,jb)) 
               tempv(jc,jk,jb) =  temp(jc,jk,jb) * (1._wp+vtmpc1*z_qv(jc,jk,jb))
              END IF
             END IF
!!$             IF (pres(jc,jk,jb) < 0.0_wp) THEN
!!$                 CALL finish ('Try with a lower top_height because &
!!$                               & you have reached p=0')
!!$             END IF 

         END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

   ! Copy to prognostic model fields and diagnostic fields
!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,jk,jc,nlen)
    DO jb = 1, nblks_c
      nlen = MERGE(nproma, npromz_c, jb /= nblks_c)

      DO jk = nlev, 1, -1
        DO jc = 1, nlen  
          ptr_nh_diag%tempv(jc,jk,jb)   = tempv(jc,jk,jb)
          ptr_nh_diag%pres(jc,jk,jb)    = pres(jc,jk,jb)
          IF (iforcing == inwp ) THEN
           ptr_nh_prog%tracer(jc,jk,jb,iqv) = z_qv(jc,jk,jb)
          END IF

        ENDDO !jc
      ENDDO !jk     
    ENDDO !jb
!$OMP END DO
!$OMP END PARALLEL
  CALL convert_thdvars(ptr_patch, ptr_nh_diag%pres, ptr_nh_diag%tempv, &
               & ptr_nh_prog%rho, ptr_nh_prog%exner, ptr_nh_prog%theta_v  )

  CALL diagnose_pres_temp (p_metrics, ptr_nh_prog,ptr_nh_prog, ptr_nh_diag,     &
                             ptr_patch, opt_calc_pres=.TRUE., opt_calc_temp=.TRUE.)


  IF (l_hydro_adjust) THEN

   CALL hydro_adjust ( ptr_patch, p_metrics, ptr_nh_prog%rho,  &
                     & ptr_nh_prog%exner, ptr_nh_prog%theta_v  )

  END IF   

    !CALL sync_patch_array(SYNC_C, ptr_patch,temp)
    !CALL sync_patch_array(SYNC_C, ptr_patch,tempv)
    !CALL sync_patch_array(SYNC_C, ptr_patch,relhum)
    !CALL sync_patch_array(SYNC_C, ptr_patch,pres)
    !CALL sync_patch_array(SYNC_C, ptr_patch,z_qv)

    DEALLOCATE(temp, tempv, relhum)
    DEALLOCATE(pres, z_qv)
    DEALLOCATE(jglayer)


  END SUBROUTINE init_nh_atmo_ana_poly
!-------------------------------------------------------------------------
!
  !>
  !! Initialization of prognostic state for an 
  !! atmosphere with layers exhibiting constant 
  !! Brunt-Vaisala-Frequency N 
  !! (N takes water vapor relhum into account).
  !!
  !! It is assumed that N refers to a moist subsaturated atmosphere 
  !!
  SUBROUTINE init_nh_atmo_ana_nconstlayers( ptr_patch, ptr_nh_prog, ptr_nh_diag, &
    &                                p_metrics, l_hydro_adjust )

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_nh_lim_area_testcases:init_nh_atmo_ana_nconstlayers'


    TYPE(t_patch), TARGET, INTENT(INOUT) :: &  !< patch on which computation is performed
      &  ptr_patch

    TYPE(t_nh_prog), INTENT(INOUT)      :: &  !< prognostic state vector
      &  ptr_nh_prog

    TYPE(t_nh_diag), INTENT(INOUT)      :: &  !< diagnostic state vector
      &  ptr_nh_diag


    TYPE(t_nh_metrics), INTENT(IN)      :: p_metrics !< NH metrics state
    LOGICAL, INTENT(IN)                 :: l_hydro_adjust !if .TRUE. hydrostatically balanced 
                                                         ! initial condition


    INTEGER        ::  jc, jb, jk,    &
                       nlen, nblks_c, npromz_c, &
                       jl, jn
    INTEGER        :: nlev        !< number of full levels

    REAL(wp), ALLOCATABLE :: theta(:,:,:), theta_v(:,:,:),exner(:,:,:)
    REAL(wp), ALLOCATABLE :: relhum(:,:,:)
    REAL(wp), ALLOCATABLE :: z_qv(:,:,:)
    INTEGER,  ALLOCATABLE :: jnlayer(:,:,:)

    REAL(wp), DIMENSION(nlayers_nconst) :: thetab, exnerb, tempb, qvb, &
                                           rhb, presb

  
    REAL(wp)  :: z_h,  temp_aux, pres_aux,  tempv_aux, zz_top, theta_top
    REAL(wp)  :: zloghuge, bvref_tconst


!--------------------------------------------------------------------
!

! First some control of the imput namelist parameters 

   IF (nlayers_nconst > max_nlayers_nconst) THEN
        CALL finish('ERROR  ** nlayers_nconst > max_nlayers_nconst !!! **')
   ENDIF
   IF (nlayers_nconst > 1) THEN
    IF (ANY(h_nconst(2:nlayers_nconst) < h_nconst(1:nlayers_nconst-1))) THEN
        CALL finish( 'ERROR  ** h_nconst (k+1) < h_nconst(k) occured &
                     & (h_nconst must be increasing) !!! ** ')
    ENDIF
   ENDIF
   IF (ANY(rh_nconst(1:nlayers_nconst-1) < 0.0_wp) .OR. &
       ANY(rh_nconst(1:nlayers_nconst-1) - rhgr_nconst(1:nlayers_nconst-1)*&
       (h_nconst(2:nlayers_nconst)-h_nconst(1:nlayers_nconst-1)) < 0.0_wp)) THEN
              CALL message('',' WARNING  *** wrong values of rh_nconst and/or rhgr_nconst &
                         & lead to relhum < 0 !!! ***  It will be set to 0')
   ENDIF
   zz_top = MAXVAL(p_metrics%z_mc)
   IF (h_nconst(nlayers_nconst) <= zz_top) THEN
         IF (rh_nconst(nlayers_nconst) < 0.0_wp .OR. &
                   rh_nconst(nlayers_nconst) - rhgr_nconst(nlayers_nconst)*&
                   (zz_top-h_nconst(nlayers_nconst)) < 0.0_wp) THEN
           CALL message('',' WARNING  *** wrong values of rh_nconst and/or &
                     & rhgr_nconst at the top lead to relhum < 0 !!! *** &
                     & It will be set to 0')
         ENDIF
   ENDIF
   IF (ANY(N_nconst(:) <= 1.e-12_wp)) THEN
    CALL finish('','ERROR  ** Brunt Vaisala Frequencies  N_nconst should be >0, &
                   & the posibility of N_nconst = 0 is still not implemented') 
   END IF  

  ! Error check for exp(-N^2/g*(z-z0)) if N^2/g*(z-z0) > log(huge(1.0_wp)-1.0)
  zloghuge = LOG(HUGE(1.0_wp)-1.0)
  DO jl=1, nlayers_nconst-1
   bvref_tconst = N_nconst(jl)**2 / grav * (h_nconst(jl+1)-h_nconst(jl))
   IF (bvref_tconst >= zloghuge) THEN
    WRITE(message_text,'(2a, i4, a)')' INPUT_ARTIFCTL: ERROR * Combination of h_nconst, N_nconst',&
             & ' will lead to floating overflow in layer jl = ',jl,' ! * '
    CALL finish('', message_text)
   END IF
  END DO
  bvref_tconst = N_nconst(nlayers_nconst)**2 / &
           grav * (zz_top - h_nconst(nlayers_nconst))
  IF (bvref_tconst >= zloghuge) THEN
    WRITE(message_text,'(2a, i4, a)')' INPUT_ARTIFCTL: ERROR * Combination of h_nconst, N_nconst',&
             & ' will lead to floating overflow in layer jl = ',nlayers_nconst,' !!! *** '
    CALL finish('', message_text)
  END IF

    ! number of vertical levels
    nlev   = ptr_patch%nlev
    ALLOCATE (theta(nproma,nlev,ptr_patch%nblks_c), &
              theta_v(nproma,nlev,ptr_patch%nblks_c), &
              exner(nproma,nlev,ptr_patch%nblks_c), &
              relhum(nproma,nlev,ptr_patch%nblks_c) )
    ALLOCATE ( z_qv(nproma,nlev,ptr_patch%nblks_c) )
    ALLOCATE ( jnlayer(nproma,nlev,ptr_patch%nblks_c) )
   nblks_c   = ptr_patch%nblks_c
   npromz_c  = ptr_patch%npromz_c
   !npromz_e  = ptr_patch%npromz_e

! set the values at h_nconst(1)
    thetab(1) = theta0_base_nconst
    exnerb(1) = (p_base_nconst/p0ref)**rd_o_cpd
    rhb(1)    = rh_nconst(1)
    presb(1)  = p_base_nconst
    tempb(1)  =  thetab(1) * exnerb(1)
    qvb(1)      = qv_rhtp(rhb(1), tempb(1), presb(1) )   


! set the values at the base of each layer
! for the base at one layer, integrate the hydrostatic equation in the layer below
  DO  jl = 2, nlayers_nconst
    thetab(jl) = thetab(jl-1)*           &
               & EXP( N_nconst(jl-1)**2*(h_nconst(jl)-h_nconst(jl-1))/grav )
! check that the value of N does not lead to a negative pressure
    IF (exnerb(jl-1) <      &    
         & (1._wp/thetab(jl)-1._wp/thetab(jl-1))*(grav/N_nconst(jl-1))**2/cpd  ) THEN
     WRITE(message_text,'(a,i3,a)')'jl: ',jl-1,  &
                        & 'value of N_nconst(jl) leads to a negative pressure'
     CALL finish('', message_text)
    END IF
    !here I use rh of the layer below
    rhb(jl)    = rh_nconst(jl-1)-rhgr_nconst(jl-1)*(h_nconst(jl)-h_nconst(jl-1))
    rhb(jl)    = MAX(0._wp,MIN(1._wp,rhb(jl)))


    IF (rhb(jl) > 0.0_wp) THEN
     !first consider constant qv (= qvb) in the layer jl-1

     exnerb(jl) = exner_nconst(exnerb(jl-1), N_nconst(jl-1),thetab(jl-1), & 
                           &  thetab(jl), qvb(jl-1), qvb(jl-1))
     presb(jl)  = p0ref*(exnerb(jl)**cpd_o_rd)
     tempb(jl)  = thetab(jl)*exnerb(jl)
     qvb(jl)    = qv_rhtp(rhb(jl), tempb(jl), presb(jl) )
     ! recalculate exner with the fisrt estimation of qvb, then recalculate qvb
     exnerb(jl) = exner_nconst(exnerb(jl-1), N_nconst(jl-1),thetab(jl-1), & 
                           &  thetab(jl), qvb(jl-1), qvb(jl))
     presb(jl)  = p0ref*(exnerb(jl)**cpd_o_rd)
     tempb(jl)  = thetab(jl)*exnerb(jl)
     qvb(jl)    = qv_rhtp(rhb(jl), tempb(jl), presb(jl) )
    ELSE
     qvb(jl)    = 0.0_wp
    END IF
! last estimation of exner
    exnerb(jl) = exner_nconst(exnerb(jl-1), N_nconst(jl-1),thetab(jl-1), & 
                           &  thetab(jl), qvb(jl-1), qvb(jl))
    presb(jl)  = p0ref*(exnerb(jl)**cpd_o_rd)
    tempb(jl)  = thetab(jl)*exnerb(jl)
! now consider the relative humidity of the base of the current layer
    rhb(jl)    = rh_nconst(jl)
    qvb(jl)    = qv_rhtp(rhb(jl), tempb(jl), presb(jl) )

  END DO

! check in the last layer
  IF (h_nconst(nlayers_nconst) <= zz_top) THEN
   theta_top= thetab(nlayers_nconst)*           &
           & EXP( N_nconst(nlayers_nconst)**2*(zz_top-h_nconst(nlayers_nconst))/grav )
   IF (exnerb(nlayers_nconst) < (1._wp/theta_top - 1._wp/thetab(nlayers_nconst))*    &
                              &  (grav/N_nconst(nlayers_nconst))**2/cpd  ) THEN
    WRITE(message_text,'(a,i3,a)') 'jl: ',nlayers_nconst,'value of N_nconst(jl) leads &
                                      &   to a negative pressure'
    CALL finish('', message_text)
   END IF  
  END IF
! ckeck finished

  IF (msg_level >= 5) THEN ! print maximum velocities in global domain
   CALL message('','N constant layers: layer, PRES, TEMP, RH, QV')
   DO  jl = 1, nlayers_nconst
     WRITE(message_text,'(10x,i3, 4e18.10)') jl,    &
                             &     presb(jl), tempb(jl), rhb(jl), qvb(jl)
     CALL message('', message_text)
   END DO

  ENDIF ! msg_level >= 5

! set the corresponding layer for all the model points

jnlayer(:,:,:)=0  
!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,jk,jc,nlen,z_h, jl)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      DO jk = nlev, 1, -1
         DO jc = 1, nlen

            z_h = p_metrics%z_mc(jc,jk,jb)
             !set the layer corresponding to this point
             IF (z_h <  h_nconst(1) ) THEN
               jnlayer(jc,jk,jb)=0
             ELSEIF (z_h >=  h_nconst(nlayers_nconst) ) THEN 
              jnlayer(jc,jk,jb)= nlayers_nconst
             ELSE
              DO jl=1,nlayers_nconst-1
               IF (z_h >= h_nconst(jl) .AND. z_h < h_nconst(jl+1)  ) THEN
                jnlayer(jc,jk,jb)=jl 
               END IF          
              END DO
             END IF
        ENDDO !jc
      ENDDO !jk     
    ENDDO !jb
!$OMP END DO
!$OMP END PARALLEL

    IF (ANY(jnlayer(:,:,:) < 0) .OR. ANY(jnlayer(:,:,:)>nlayers_nconst) ) THEN
         CALL finish ('corresponding layer has not been found for some model points')
    END IF
!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,jk,jc,nlen,z_h,jn,tempv_aux, pres_aux,temp_aux)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      DO jk = nlev, 1, -1
         DO jc = 1, nlen

            z_h = p_metrics%z_mc(jc,jk,jb)
            jn=jnlayer(jc,jk,jb)
            IF (jn == 0) THEN
              !consider a isothermal atmosphere below h_nconst(1) with qv=qvb(1)
              z_qv(jc,jk,jb) = qvb(1)
              tempv_aux = tempb(1)*(1._wp+vtmpc1*qvb(1))
              pres_aux  = presb(1) * &
                & EXP(-grav_o_rd*(z_h-h_nconst(1))/tempv_aux)
              exner(jc,jk,jb)= (pres_aux/p0ref)**rd_o_cpd  
              theta(jc,jk,jb)= tempb(1)/exner(jc,jk,jb)
              theta_v(jc,jk,jb)=theta(jc,jk,jb)*(1._wp+vtmpc1*qvb(1))
            ELSE
              theta(jc,jk,jb)=thetab(jn)*          &
               & EXP( N_nconst(jn)**2*(z_h-h_nconst(jn))/grav )

              IF (ABS(rhgr_nconst(jn)) > 1.e-20) THEN
               relhum(jc,jk,jb)=rh_nconst(jn)-rhgr_nconst(jn)*(z_h-h_nconst(jn))
              ELSE
               relhum(jc,jk,jb)=rh_nconst(jn)
              END IF
              relhum(jc,jk,jb)=MAX(0._wp,MIN(1._wp,relhum(jc,jk,jb)))
              IF (jk < nlev .AND. jn == jnlayer(jc,MIN(jk+1,nlev),jb) ) THEN
               ! in this case we integrate starting in the level bellow, instead of 
               !  starting in the base of the layer

               IF (relhum(jc,jk,jb) > 1.e-20_wp) THEN
                ! first consider qv constant (= z_qv(jc, jk+1,jb))
                exner(jc,jk,jb)= exner_nconst(exner(jc,jk+1,jb), N_nconst(jn),  &
                           &  theta(jc,jk+1,jb), theta(jc,jk,jb),              &
                           &  z_qv(jc, jk+1,jb), z_qv(jc, jk+1,jb))
                pres_aux  = p0ref*(exner(jc,jk,jb)**cpd_o_rd)
                temp_aux  = theta(jc,jk,jb)*exner(jc,jk,jb)
                ! 1st estimation of  qv, then earth_radious-estimate exner
                z_qv(jc,jk,jb)= qv_rhtp(relhum(jc,jk,jb), temp_aux, pres_aux )
                exner(jc,jk,jb)= exner_nconst(exner(jc,jk+1,jb), N_nconst(jn),  &
                            &  theta(jc,jk+1,jb), theta(jc,jk,jb),              &
                            &  z_qv(jc, jk+1,jb), z_qv(jc, jk,jb))
                pres_aux  = p0ref*(exner(jc,jk,jb)**cpd_o_rd)
                 temp_aux  = theta(jc,jk,jb)*exner(jc,jk,jb)
                ! earth_radious-estimate qv, then earth_radious-estimate exner
                z_qv(jc,jk,jb)= qv_rhtp(relhum(jc,jk,jb), temp_aux, pres_aux )
               ELSE
                z_qv(jc,jk,jb)= 0.0_wp
               END IF
               exner(jc,jk,jb)= exner_nconst(exner(jc,jk+1,jb), N_nconst(jn),  &
                           &  theta(jc,jk+1,jb), theta(jc,jk,jb),              &
                           &  z_qv(jc, jk+1,jb), z_qv(jc, jk,jb))
               pres_aux  = p0ref*(exner(jc,jk,jb)**cpd_o_rd)
               temp_aux  = theta(jc,jk,jb)*exner(jc,jk,jb)
               ! earth_radious-estimate qv
               z_qv(jc,jk,jb)= qv_rhtp(relhum(jc,jk,jb), temp_aux, pres_aux )
               theta_v(jc,jk,jb)=theta(jc,jk,jb)*(1._wp+vtmpc1*z_qv(jc,jk,jb))

               
              ELSE
               !WRITE(*,*) "from the layer base"
               ! here we integrate from the base of the layer
               !
               ! here I could earth_radious-estimate the values at the base of the layer
               ! integrating from (jc,jk+1,jb) to the base of the layer 
               ! in the case jn>1
!!$                IF (jn>1  .AND. exnerb(jn) > exner(jc,jk+1,jb)) THEN
!!$                  ! check that exner(jn)<exner(jc,jk+1,jb)
!!$                  CALL finish ('base layer has larger pressure than the level &
!!$                           &below exnerb(jn) was not well approximated')
!!$                             
!!$                END IF
               IF (relhum(jc,jk,jb) > 1.e-20_wp) THEN
                !
                ! first consider qv constant (= qvb(jn))
                exner(jc,jk,jb)= exner_nconst(exnerb(jn), N_nconst(jn),thetab(jn), & 
                            &  theta(jc,jk,jb), qvb(jn), qvb(jn))
                pres_aux  = p0ref*(exner(jc,jk,jb)**cpd_o_rd)
                temp_aux  = theta(jc,jk,jb)*exner(jc,jk,jb)
                ! 1st estimation of  qv, then earth_radious-estimate exner
                z_qv(jc,jk,jb)= qv_rhtp(relhum(jc,jk,jb), temp_aux, pres_aux )
                exner(jc,jk,jb)= exner_nconst(exnerb(jn), N_nconst(jn),thetab(jn), & 
                            &  theta(jc,jk,jb), qvb(jn),z_qv(jc,jk,jb) )
                pres_aux  = p0ref*(exner(jc,jk,jb)**cpd_o_rd)
                temp_aux  = theta(jc,jk,jb)*exner(jc,jk,jb)
                ! earth_radious-estimate qv, then earth_radious-estimate exner
                z_qv(jc,jk,jb)= qv_rhtp(relhum(jc,jk,jb), temp_aux, pres_aux )
               ELSE
                z_qv(jc,jk,jb)= 0.0_wp
               END IF
               exner(jc,jk,jb)= exner_nconst(exnerb(jn), N_nconst(jn),thetab(jn), & 
                           &  theta(jc,jk,jb), qvb(jn),z_qv(jc,jk,jb) )
               pres_aux  = p0ref*(exner(jc,jk,jb)**cpd_o_rd)
               temp_aux  = theta(jc,jk,jb)*exner(jc,jk,jb)
               ! earth_radious-estimate qv
               z_qv(jc,jk,jb)= qv_rhtp(relhum(jc,jk,jb), temp_aux, pres_aux )
               theta_v(jc,jk,jb)=theta(jc,jk,jb)*(1._wp+vtmpc1*z_qv(jc,jk,jb))
               
              END IF
           END IF
!!$             IF (exner(jc,jk,jb) < 0.0_wp) THEN
!!$                 CALL finish ('Try with a lower top_height, because &
!!$                               &you have reached p=0')
!!$             END IF             
             
         END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL


   ! Copy to prognostic model fields
!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,jk,jc,nlen)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      DO jk = nlev, 1, -1
        DO jc = 1, nlen  
          ptr_nh_prog%theta_v(jc,jk,jb)    = theta_v(jc,jk,jb)
          ptr_nh_prog%exner(jc,jk,jb)      = exner(jc,jk,jb)
          IF ( iforcing == inwp ) THEN
           ptr_nh_prog%tracer(jc,jk,jb,:) = 0._wp
           ptr_nh_prog%tracer(jc,jk,jb,iqv) = z_qv(jc,jk,jb)
          END IF

          ptr_nh_prog%rho(jc,jk,jb)  = ptr_nh_prog%exner(jc,jk,jb)**cvd_o_rd*p0ref &
                                       /rd/ptr_nh_prog%theta_v(jc,jk,jb)
        ENDDO !jc
      ENDDO !jk     
    ENDDO !jb
!$OMP END DO
!$OMP END PARALLEL

    CALL diagnose_pres_temp (p_metrics, ptr_nh_prog,ptr_nh_prog, ptr_nh_diag,     &
                               ptr_patch, opt_calc_pres=.TRUE., opt_calc_temp=.TRUE.)


    IF (l_hydro_adjust) THEN
      !
      ! remark: calling hydro_adjust would be equally fine here
      ! 
      CALL hydro_adjust_iterative ( p_patch       = ptr_patch,                    & !in
        &                           p_nh_metrics  = p_metrics,                    & !in
        &                           temp_ini      = ptr_nh_diag%temp,             & !in
        &                           rh_ini        = relhum,                       & !in
        &                           exner         = ptr_nh_prog%exner,            & !inout
        &                           theta_v       = ptr_nh_prog%theta_v,          & !inout
        &                           rho           = ptr_nh_prog%rho,              & !inout
        &                           qv            = ptr_nh_prog%tracer(:,:,:,iqv),& !inout
        &                           luse_exner_fg = .TRUE.,                       & !in
        &                           opt_exner_lbc = ptr_nh_prog%exner(:,nlev,:)   ) !in
    END IF   


    DEALLOCATE(theta, theta_v, relhum)
    DEALLOCATE(exner, z_qv)
    DEALLOCATE(jnlayer)
  END SUBROUTINE init_nh_atmo_ana_nconstlayers

!-------------------------------------------------------------------------
!
  !>
  !! Initialization of the wind profile for the 
  !!  limited area testcases 
  !! 
  !!
  SUBROUTINE init_nh_anaprof_uv( ptr_patch, vn, w,    &
    &                                p_metrics, p_int )
    TYPE(t_patch), TARGET, INTENT(INOUT) :: &  !< patch on which computation is performed
      &  ptr_patch

    REAL(wp), INTENT(INOUT) :: vn(:,:,:)    ! edge-normal wind component (m/s)

    REAL(wp), INTENT(INOUT) :: w(:,:,:)  ! vertical wind component (m/s)
    TYPE(t_nh_metrics), INTENT(IN)      :: p_metrics !< NH metrics state
    TYPE(t_int_state), INTENT(IN)       :: p_int

    REAL(wp), ALLOCATABLE :: z_me(:,:,:)

    REAL(wp)              :: z_u, z_v
    INTEGER               :: je,jb,jk,jl,jn
    INTEGER               :: i_startblk, i_startidx, i_endidx, nblks_e 
    INTEGER               :: nlev

!--------------------------------------------------------------------
!
    ! number of vertical levels
    nlev   = ptr_patch%nlev

!    nblks_c   = ptr_patch%nblks_c
!    npromz_c  = ptr_patch%npromz_c
    nblks_e   = ptr_patch%nblks_e

    SELECT CASE (itype_anaprof_uv)

      CASE(1)
      ! arbitrary number of constant gradient U(z) layers 
      ! horizontal normal components of the velocity
      ! initialize horizontal velocities

       ALLOCATE(z_me(nproma,nlev,ptr_patch%nblks_e))
       z_me(:,:,:) = 0.0_wp

       ! Compute geometric height at edge points
       CALL cells2edges_scalar(p_metrics%z_mc, ptr_patch, &
               p_int%c_lin_e, z_me)

       CALL sync_patch_array(SYNC_E,ptr_patch,z_me)
       i_startblk = ptr_patch%edges%start_blk(2,1)

       ! horizontal normal components of the velocity
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,z_u,jl,jn)
       DO jb = i_startblk, nblks_e

        CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                         i_startidx, i_endidx, 2)

        DO jk = 1, nlev
          DO je = i_startidx, i_endidx

             !set the layer corresponding to this point
             IF (z_me(je,jk,jb) <  h_linwind(1) ) THEN
               jn=0
             ELSEIF (z_me(je,jk,jb) >=  h_linwind(nlayers_linwind) ) THEN 
              jn= nlayers_linwind
             ELSE
              DO jl=1,nlayers_linwind-1
               IF (z_me(je,jk,jb) >= h_linwind(jl) .AND. z_me(je,jk,jb) < h_linwind(jl+1)) THEN
                jn=jl 
                !EXIT
               END IF          
              END DO
             END IF
             IF (jn <= 0 .OR. jn > nlayers_linwind) THEN
                CALL finish ('corresponding layer has not been found')
             END IF
            
            z_u = u_linwind(jn) +  ugr_linwind(jn)*    &
                                &  (z_me(je,jk,jb)-h_linwind(jn))   !v component is zero
            vn(je,jk,jb) = &
             z_u * ptr_patch%edges%primal_normal(je,jb)%v1
          ENDDO !je
        ENDDO !jk
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL
       DEALLOCATE(z_me)

      CASE(2)
      ! constant zonal wind
       z_u = vel_const
       i_startblk = ptr_patch%edges%start_blk(2,1)

       ! horizontal normal components of the velocity
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je)
       DO jb = i_startblk, nblks_e

        CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                         i_startidx, i_endidx, 2)

        DO jk = 1, nlev
          DO je = i_startidx, i_endidx

           vn(je,jk,jb) =  &
             z_u * ptr_patch%edges%primal_normal(je,jb)%v1
          ENDDO !je
        ENDDO !jk
     ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

      CASE(3)
      ! constant meridional wind
       z_v = vel_const
       i_startblk = ptr_patch%edges%start_blk(2,1)

       ! horizontal normal components of the velocity
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je)
       DO jb = i_startblk, nblks_e

        CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                         i_startidx, i_endidx, 2)

        DO jk = 1, nlev
          DO je = i_startidx, i_endidx

           vn(je,jk,jb) =  &
             z_v * ptr_patch%edges%primal_normal(je,jb)%v2
          ENDDO !je
        ENDDO !jk
     ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

    END SELECT

 ! initialize vertical velocity
   CALL init_w(ptr_patch, p_int, vn, p_metrics%z_ifc, w)
   CALL sync_patch_array(SYNC_C, ptr_patch, w)
   !CALL sync_patch_array(SYNC_E, ptr_patch, vn)

  END SUBROUTINE init_nh_anaprof_uv
!-------------------------------------------------------------------------
!
  !>
  !! Initialization of topography 
  !!

  SUBROUTINE init_nh_topo_ana( ptr_patch, lplane, topo_c, nblks_c, npromz_c)

    TYPE(t_patch), TARGET,INTENT(INOUT) :: &  !< patch on which computation is performed
      &  ptr_patch

    INTEGER,  INTENT (IN) ::  nblks_c, npromz_c
    LOGICAL,  INTENT(IN)  :: lplane
    REAL(wp), INTENT(INOUT) :: topo_c    (nproma,nblks_c)

    ! local variables

    REAL(wp)       :: z_lon, z_lat, z_lonc, z_latc
    REAL(wp)       :: z_dx, z_dy, z_deltay
    INTEGER        :: jc, jb, nlen
!--------------------------------------------------------------------

    z_lonc = mount_lonc_deg*deg2rad
    z_latc = mount_latc_deg*deg2rad

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,z_lat,z_lon,z_dx,z_dy,z_deltay ) 
      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen =  npromz_c
        ENDIF
        DO jc = 1, nlen
          IF ( itopo==0 ) THEN

           z_lat   = ptr_patch%cells%center(jc,jb)%lat
           z_lon   = ptr_patch%cells%center(jc,jb)%lon
           CALL xy_distances(z_lon, z_lat, z_lonc, z_latc, &
             & 0.0_wp, 0.0_wp, z_dx, z_dy, grid_sphere_radius, lplane)

           SELECT CASE (itype_topo_ana)
           CASE (1)  !schaer mountain

            topo_c(jc,jb) = schaer_h0 * EXP (-(z_dx/schaer_a)**2)* &
                                    ( COS(pi*z_dx/schaer_lambda)**2)
            IF (ABS(z_dy) >= halfwidth_2dm ) THEN
             z_deltay=ABS(z_dy)-halfwidth_2dm
             topo_c(jc,jb) = topo_c(jc,jb)*EXP(-(z_deltay/schaer_a)**2)
            END IF
           CASE (2)  !gaussian_2d

            topo_c(jc,jb) = m_height * EXP (-LOG(2.0_wp)*(z_dx/m_width_x)**2)
            IF (ABS(z_dy) >= halfwidth_2dm ) THEN
             z_deltay=ABS(z_dy)-halfwidth_2dm
             topo_c(jc,jb) = m_height * EXP (-LOG(2.0_wp)*(z_dx**2+z_deltay**2)/m_width_x**2)
            END IF
           CASE (3)  !gaussian_3d

            topo_c(jc,jb) = m_height * EXP (-LOG(2.0_wp)*          &
                                            ((z_dx/m_width_x)**2+ (z_dy/m_width_y)**2))
           CASE DEFAULT
            topo_c(jc,jb) = 0._wp
           END SELECT

          END IF
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE init_nh_topo_ana
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!
  !>
  !! lonc, latc can be the center of a mountain, buble,... 
  !! Calculate distances from a lon,lat point to the center point 
  !!  lonc,latc in the direction of the main axis of the 
  !!  mountain/buble
  !! If rotangle is zero, then it calculates distances from a lon lat point 
  !!  to the center point lonc, latc in the zonal and meridional directions
  !! This is a translation from SUBROUTINE hill_rot_coords in COSMO
  !! 
  !! It uses the spherical law of sines and the spherical law of cosines
  !!
  SUBROUTINE xy_distances(lon, lat, lonc, latc, rotangle, height, dx, dy, &
    & sphere_radius, lplane)

   REAL(wp), INTENT (IN) :: lon, lat   !lon and lat of the point (in radians)
   REAL(wp), INTENT (IN) :: lonc, latc !lon and lat of the center point (in radians)
   REAL(wp), INTENT (IN) :: rotangle   !rotation angle of main hill/bubble y-axis 
                                       ! clockwise relative to north in rad
   REAL(wp), INTENT (IN) :: height     !vertical height of the point (in meters)
   REAL(wp), INTENT (OUT):: dx, dy     !distances in meters from (lon,lat) to
                                       ! (lonc,latc) in zonal and meridional directions
   REAL(wp) :: sphere_radius
   LOGICAL, INTENT(in) :: lplane

   ! Local variables
   REAL(wp)              :: z_d_lon, z_d_lat
   REAL(wp)              :: z_cosd, z_d, z_delta
   REAL(wp)              :: z_tmp, z_cangle, z_cos_arg

!--------------------------------------------------------------------
!
   
  IF (lplane) THEN
    z_d_lon = sphere_radius * (lon-lonc)
    z_d_lat = sphere_radius * (lat-latc)

    dx = z_d_lon * COS(rotangle) - z_d_lat * SIN(rotangle) 
    dy = z_d_lon * SIN(rotangle) + z_d_lat * COS(rotangle)

  ELSE
    
    z_cosd = SIN(latc)*SIN(lat)+COS(latc)*COS(lat)*COS(lon-lonc)
    z_d    = ACOS(z_cosd)

    z_cangle = ACOS(z_cosd)
    IF (ABS(z_cangle) < 1.e-20_wp) z_cangle = 1.e-20_wp
    z_cos_arg = (SIN(lat)-SIN(latc)*z_cosd) / (COS(latc)*SIN(z_cangle))
    z_cos_arg = MAX(MIN(z_cos_arg, 1.0_wp), -1.0_wp)

 
    z_delta= ACOS(z_cos_arg)
    IF (lonc > lon) z_delta = 2.0_wp*pi - z_delta
    !.. take rotation into account:
    z_delta= z_delta- rotangle

     dx = ASIN(SIN(z_d)*SIN(z_delta))
      z_tmp = 1.0_wp - SIN(z_delta)*SIN(dx)*SIN(z_d) !cosine square of dx
      IF (ABS(z_tmp) < 1.e-20_wp) THEN
        dy = 0.0_wp
      ELSE
        dy = ACOS(COS(dx)*COS(z_d) / z_tmp)
        IF (latc > lat) dy = -dy
      END IF
      dx = (sphere_radius + height) * dx
      dy = (sphere_radius + height) * dy
  END IF


  END SUBROUTINE xy_distances

!-------------------------------------------------------------------------
!--------------------------------------------------------------------

! Calculate the exner function integrating the hydrostatic equation.
!  It considers we are within a layer of constant Brunt Vaisala frequency.
!  It considers that qv is constant within the two levels of the integration, 
!  in practice it uses the mean value of qv
  REAL(wp)  FUNCTION exner_nconst(exnerb,N,thetab, theta, qvb, qv)
 
   REAL(wp), INTENT (IN):: exnerb, N, thetab, theta, qvb, qv !exnerb, thetab and qvb are 
                                                             ! the values at one level
                                                             ! N is the Brunt-Vaisala-freq
                                                             ! theta and qv are the values at 
                                                             ! the level for which we want to calculate
                                                             ! the exner function
   REAL(wp)             :: qv_mean, factor1, factor2

   qv_mean = 0.5_wp * (qvb+qv)
   factor1  = (1._wp + vtmpc1*qv_mean )
   factor2  = (grav/N)**2/(cpd*factor1)
   exner_nconst = exnerb+ factor2*(1._wp/theta - 1._wp/thetab)

  END FUNCTION exner_nconst
!--------------------------------------------------------------------
! Calculates the pressure integrating the hydrostatic equation.
!  It considers we are within a polytropic atmosphere layer of lapse 
!   rate G (positive G means that the temperature decreases with height).
!  It considers that qv is constant within the two levels of the integration, 
!  in practice it uses the mean value of qv
  REAL(wp)  FUNCTION p_poly(p_base,G,temp_base, height_base, height,  qv_base, qv)
 
   REAL(wp), INTENT (IN):: p_base, G, temp_base             !temp_base, p_base and qv_base are 
                                                             ! the values at one level
                                                             ! G is the lapse rate
                                                             ! temp and qv are the values at 
                                                             ! the level for which we want to 
                                                             ! calculate the pressure
   REAL(wp), INTENT (IN):: height_base, height               ! the corresponding heights
   REAL(wp), INTENT (IN):: qv_base, qv
   REAL(wp)             :: temp, delta_h
   REAL(wp)             :: qv_mean, factorqv, z_exp

   qv_mean = 0.5_wp * (qv_base+qv)
   factorqv  = (1._wp + vtmpc1*qv_mean )
   delta_h = height - height_base
   IF (G > 1.e-20_wp) THEN
    temp = temp_base - G * delta_h
    z_exp  = grav/G/rd/factorqv
    p_poly = p_base * (temp/temp_base)**z_exp
   ELSE  !isothermal atmosphere
    temp = temp_base
    p_poly = p_base * EXP( -grav_o_rd*delta_h/temp/factorqv ) 
   END IF

  END FUNCTION p_poly
!--------------------------------------------------------------------
!
! Calculate specific humidity from relative humidity, temp and pres
  REAL(wp)  FUNCTION qv_rhtp(rh, temp, pres)
   REAL(wp), INTENT (IN)   :: rh, temp, pres
   REAL (wp)               :: e_aux

      IF (temp > tmelt) THEN
        e_aux     = rh*sat_pres_water(temp)
      ELSE
        e_aux     = rh*sat_pres_ice(temp)
      ENDIF
      qv_rhtp     = MAX(0._wp,MIN(1._wp,spec_humi(e_aux,pres)))
   END FUNCTION qv_rhtp

!--------------------------------------------------------------------
  END MODULE mo_nh_lim_area_testcases
