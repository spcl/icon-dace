! Subroutines to initialize the Baldauf, Brdar (2013) QJRMS test case
! (linear gravity/sound waves expansion in a channel)
! for the NH-Core in a quasi-2D limited area mode
!
! Literature
! - Baldauf, Brdar (2013):
!    An analytic solution for linear gravity waves in a channel
!    as a test for numerical models using the non-hydrostatic, compressible Euler equations,
!    Quart. J. Royal Met. Soc., 139, 1977-1989
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

MODULE mo_nh_bb13_exp

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
   USE mo_physical_constants,   ONLY: rd, rd_o_cpd, cvd_o_rd,   &
                                      p0ref, grav
   USE mo_math_constants,       ONLY: pi
   USE mo_math_utilities,       ONLY: gc2cc
   USE mo_math_types,           ONLY: t_cartesian_coordinates
   USE mo_model_domain,         ONLY: t_patch
   USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
   USE mo_run_config,           ONLY: iqv
   USE mo_impl_constants,       ONLY: min_rlcell_int
   USE mo_parallel_config,      ONLY: nproma
   !USE mo_exception,            ONLY: finish, message_text, warning
   USE mo_exception,            ONLY: message
   USE mo_loopindices,          ONLY: get_indices_e
   USE mo_nh_diagnose_pres_temp,ONLY: diagnose_pres_temp
   !USE mo_extpar_config,        ONLY: itopo
   USE mo_sync,                 ONLY: sync_patch_array, SYNC_C
   USE mo_vertical_coord_table, ONLY: vct_a
   USE mo_hydro_adjust,         ONLY: hydro_adjust

   USE mo_nh_wk_exp,            ONLY: u_infty_wk, bub_hor_width, bub_amp
   USE mo_nh_testcases_nml,     ONLY: bubctr_x
   IMPLICIT NONE

   PUBLIC  :: init_nh_env_bb13, init_nh_bubble_bb13
  
   PRIVATE


   REAL(wp), PARAMETER :: T_bb13 = 250.0_wp    ! background atmosphere [K]
  
!--------------------------------------------------------------------

   CONTAINS
!-------------------------------------------------------------------------
!

  !>
  !! Initialization of prognostic state vector for the Baldauf, Brdar (2013) test case 
  !! here: define the horizontally homogenous environmental state
  !!
  SUBROUTINE init_nh_env_bb13( ptr_patch, ptr_nh_prog, ptr_nh_diag, &
    &                                p_metrics, l_hydro_adjust )

    TYPE(t_patch), TARGET, INTENT(INOUT) :: &  !< patch on which computation is performed
      &  ptr_patch

    TYPE(t_nh_prog), INTENT(INOUT)      :: &  !< prognostic state vector
      &  ptr_nh_prog

    TYPE(t_nh_diag), INTENT(INOUT)      :: &  !< diagnostic state vector
      &  ptr_nh_diag


    TYPE(t_nh_metrics), INTENT(IN)      :: p_metrics !< NH metrics state
    LOGICAL, INTENT(IN)                 :: l_hydro_adjust !if .TRUE. hydrostatically balanced 
                                                         ! initial condition

    INTEGER        ::  jc, jb, jk, je,   &
                       nlen
    INTEGER        :: i_startidx, i_endidx, i_startblk

    REAL(wp)       :: zu, zv

    REAL(wp), DIMENSION(ptr_patch%nlev) :: theta, exner, pres, qv, theta_v, temp, z_full

!--------------------------------------------------------------------
!

    ! height of main levels (no orography)
    DO jk = 1, ptr_patch%nlev
      z_full(jk) = 0.5_wp*( vct_a(jk) + vct_a(jk+1) )
    ENDDO

    ! profiles for T = const 
    DO jk = 1, ptr_patch%nlev
      temp(jk)  = T_bb13
      pres(jk)  = p0ref * exp( - grav / ( Rd * T_bb13 ) * z_full(jk) )

      qv(jk)    = 0.0_wp

      exner(jk) = ( pres(jk) / p0ref )**( Rd_o_cpd )
      theta(jk) = temp(jk) / exner(jk)

      theta_v(jk) = theta(jk)

    ENDDO

    DO jb = 1, ptr_patch%nblks_c
      IF (jb /= ptr_patch%nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = ptr_patch%npromz_c
      ENDIF

      DO jk = ptr_patch%nlev, 1, -1
        DO jc = 1, nlen  
          ptr_nh_prog%theta_v(jc,jk,jb)     = theta_v(jk)
          ptr_nh_prog%exner  (jc,jk,jb)     = exner(jk)
          ptr_nh_prog%tracer (jc,jk,jb,iqv) = qv(jk)
          ptr_nh_prog%rho    (jc,jk,jb)     = ptr_nh_prog%exner(jc,jk,jb)**cvd_o_rd * p0ref   &
                                                         /rd/ptr_nh_prog%theta_v(jc,jk,jb)
        ENDDO  !jc
      ENDDO  !jk     
    ENDDO  !jb

    ! Achtung: das muss noch numerisch hydrostatisch balanciert werden!
    if ( l_hydro_adjust ) then
      CALL hydro_adjust(ptr_patch, p_metrics, ptr_nh_prog%rho, ptr_nh_prog%exner, ptr_nh_prog%theta_v )
    end if

    ! initialize horizontal velocities
    i_startblk = ptr_patch%edges%start_blk(2,1)

    ! Kontrolloutput:
    !jb = 1
    !jc = 1
    !DO jk = ptr_patch%nlev, 1, -1
    ! write(*,'(A,F13.3,4E13.6)') "BB13,env: ", z_full(jk), &
    !      ptr_nh_prog%theta_v(jc,jk,jb), ptr_nh_prog%exner(jc,jk,jb), ptr_nh_prog%rho(jc,jk,jb),   &
    !      ptr_nh_prog%theta_v(jc,jk,jb)*ptr_nh_prog%exner(jc,jk,jb)
    !ENDDO


    write(*,*)  "u_infty_wk=", u_infty_wk

    ! horizontal normal components of the velocity
    DO jb = i_startblk, ptr_patch%nblks_e

      CALL get_indices_e(ptr_patch, jb, i_startblk, ptr_patch%nblks_e, &
                         i_startidx, i_endidx, 2)

      DO jk = 1, ptr_patch%nlev
        DO je = i_startidx, i_endidx
          zu = u_infty_wk
          zv = 0.0_wp
          ptr_nh_prog%vn(je,jk,jb) = zu * ptr_patch%edges%primal_normal(je,jb)%v1   &
            &                      + zv * ptr_patch%edges%primal_normal(je,jb)%v2
        ENDDO  !je
      ENDDO  !jk
    ENDDO  !jb

    ptr_nh_prog%w(:,:,:) = 0.0_wp

    CALL diagnose_pres_temp (p_metrics, ptr_nh_prog,ptr_nh_prog, ptr_nh_diag,     &
                             ptr_patch, opt_calc_pres=.TRUE., opt_calc_temp=.TRUE.)


    CALL sync_patch_array(SYNC_C, ptr_patch, ptr_nh_prog%w)

  END SUBROUTINE init_nh_env_bb13

!--------------------------------------------------------------------
!-------------------------------------------------------------------------
!    

  !> 
  !! warm, dry bubble test case initialization of Baldauf, Brdar (2013) QJRMS
  !!
  SUBROUTINE init_nh_bubble_bb13( ptr_patch, p_metrics, ptr_nh_prog, ptr_nh_diag )

    IMPLICIT NONE

    TYPE(t_patch), TARGET, INTENT(INOUT):: &  !< patch on which computation is performed
      &  ptr_patch

    TYPE(t_nh_metrics), INTENT(IN)     :: p_metrics !< NH metrics state
    TYPE(t_nh_prog), INTENT(INOUT)      :: &  !< prognostic state vector
      &  ptr_nh_prog

    TYPE(t_nh_diag), INTENT(INOUT)      :: &  !< diagnostic state vector
      &  ptr_nh_diag
 ! local variables  

    REAL(KIND=wp)   :: model_height

    REAL(KIND=wp)            :: z_full
    REAL(KIND=wp)            :: dT_breth, dT
    REAL(KIND=wp)            :: delta_B, fac_breth

    INTEGER        ::  jc, jb, jk, nlen 
    TYPE(t_cartesian_coordinates)   :: p

    call message( "init_nh_bubble_bb13", "ACHTUNG: model_height sollte extern vorgegeben sein!" )

    model_height = 10000.0_wp

    delta_B = grav / Rd / T_bb13

    DO jb = 1, ptr_patch%nblks_c
      IF (jb /= ptr_patch%nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = ptr_patch%npromz_c
      ENDIF

      DO jk = ptr_patch%nlev, 1, -1
        DO jc = 1, nlen

          ! Koordinaten dieses Punktes:
          p = gc2cc( ptr_patch%cells%center(jc,jb), ptr_patch%geometry_info ) 
          !z_full = 0.5_wp * ( p_metrics%z_ifc(jc,jk,  jb)     &
          !  &               + p_metrics%z_ifc(jc,jk+1,jb) )
          z_full = p_metrics%z_mc(jc,jk,jb)

          fac_breth = EXP( 0.5 * delta_B * z_full ) 

          dT_breth = bub_amp                          &
            &    * SIN( z_full * pi / model_height )    &
            &    * EXP( - ( ( p%x(1)- bubctr_x ) / bub_hor_width )**2 )

          dT = dT_breth * fac_breth

          ! T(i,j,k) = T(i,j,k) + dT
          ptr_nh_prog%theta_v(jc,jk,jb) = ( T_bb13 + dT ) / ptr_nh_prog%exner(jc,jk,jb)

          ptr_nh_prog%rho    (jc,jk,jb) = ptr_nh_prog%exner(jc,jk,jb)**cvd_o_rd * p0ref   &
                                                    /rd/ptr_nh_prog%theta_v(jc,jk,jb)

        ENDDO
        
      ENDDO
    ENDDO

    CALL diagnose_pres_temp ( p_metrics, ptr_nh_prog,         &
            &                     ptr_nh_prog, ptr_nh_diag,   &
            &                     ptr_patch,                  &
            &                     opt_calc_temp=.TRUE.,       &
            &                     opt_calc_pres=.FALSE.,      &
            &                     opt_rlend=min_rlcell_int )
     
  END SUBROUTINE init_nh_bubble_bb13

!--------------------------------------------------------------------
END MODULE mo_nh_bb13_exp

