! Configuration for interpolation and reconstruction.
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

MODULE mo_interpol_config

  USE mo_kind,                ONLY: wp
  USE mo_math_constants,      ONLY: ln2
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_grid_geometry_info,  ONLY: t_grid_geometry_info, planar_torus_geometry
  USE mo_grid_config,         ONLY: grid_rescale_factor
  USE mo_run_config,          ONLY: msg_level


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: llsq_lin_consv, llsq_high_consv, lsq_high_ord                       !< variables
  PUBLIC :: rbf_vec_kern_c, rbf_vec_kern_v, rbf_vec_kern_e, rbf_vec_kern_ll     !< variables
  PUBLIC :: rbf_vec_scale_c, rbf_vec_scale_e, rbf_vec_scale_v, rbf_vec_scale_ll !< variables
  PUBLIC :: nudge_max_coeff, nudge_efold_width, nudge_zone_width                !< variables
  PUBLIC :: rbf_vec_dim_c, rbf_vec_dim_v, rbf_vec_dim_e, rbf_c2grad_dim         !< variables
  PUBLIC :: lsq_lin_set, lsq_high_set                                           !< variables
  PUBLIC :: rbf_dim_c2l, l_intp_c2l, l_mono_c2l                                 !< variables
  PUBLIC :: rbf_scale_mode_ll                                                   !< variables
  PUBLIC :: support_baryctr_intp                                                !< variables
  PUBLIC :: lreduced_nestbdry_stencil                                           !< variables
  PUBLIC :: configure_interpolation                                             !< subroutine

  PUBLIC :: t_lsq_set


  TYPE t_lsq_set
    LOGICAL :: l_consv            ! flag to determine whether the least squares
                                  ! reconstruction should be conservative

    INTEGER :: dim_c,        &    ! parameter determining the size of the lsq stencil
      &        dim_unk            ! parameter determining the dimension of the solution
                                  ! vector (== number of unknowns) of the lsq system

    INTEGER :: wgt_exp            ! least squares weighting exponent
  END TYPE t_lsq_set


  !>
  !!
  !TYPE :: t_interpol_config

    ! namelist variables
    LOGICAL  :: llsq_lin_consv      ! conservative (TRUE) or non-conservative (FALSE)
                                    ! linear least squares reconstruction
    LOGICAL  :: llsq_high_consv     ! conservative (TRUE) or non-conservative (FALSE)
                                    ! high order least squares reconstruction 
                                                                                 
    INTEGER  :: lsq_high_ord        ! specific order for higher order lsq        
                                                                                 
    INTEGER  :: rbf_vec_kern_c,   & ! parameter determining the type             
       &        rbf_vec_kern_v,   & ! of vector rbf kernel                       
       &        rbf_vec_kern_e,   &                                             
       &        rbf_vec_kern_ll                                                  

    ! "rbf_scale_mode_ll": mode, how the RBF shape parameter is
    ! determined for lon-lat interpolation.
    !
    ! 1 : lookup table based on grid level (default)
    ! 2 : determine automatically
    !
    INTEGER :: rbf_scale_mode_ll

    ! Parameter fields determining the scale factor used by the vector rbf       
    ! interpolator.                                                              
    ! Note: these fields are defined on each grid level; to allow the namelist input
    ! going from 1 to depth (rather than from start_lev to end_lev), the namelist input         
    ! fields defined here differ from those used in the model                    
  
    REAL(wp) :: rbf_vec_scale_c(max_dom)
    REAL(wp) :: rbf_vec_scale_e(max_dom)
    REAL(wp) :: rbf_vec_scale_v(max_dom)
    REAL(wp) :: rbf_vec_scale_ll(max_dom)

    ! Namelist variables setting up the lateral boundary nudging (applicable to limited-area   
    ! runs and one-way nesting). The nudging coefficients start with nudge_max_coeff in        
    ! the cell row bordering to the boundary interpolation zone, and decay exponentially       
    ! with nudge_efold_width (in units of cell rows)
  
    REAL(wp) :: nudge_max_coeff, nudge_efold_width
    INTEGER  :: nudge_zone_width    ! total width of nudging zone in units of cell rows
                                    ! if < 0 the patch boundary_depth_index is used

    INTEGER  ::  rbf_vec_dim_c,    & ! parameter determining the size            
       &         rbf_vec_dim_v,    & ! of vector rbf stencil                     
       &         rbf_vec_dim_e,    & !                                           
       &         rbf_c2grad_dim      ! ... and for cell-to-gradient reconstruction

    ! Flag. If .TRUE. we directly interpolate from cell centers to
    ! lon-lat points, otherwise we do gradient interpolation and
    ! reconstruction:
    LOGICAL :: l_intp_c2l

    ! monotonicity can be enforced by demanding that the interpolated 
    ! value is not higher or lower than the stencil point values.
    LOGICAL :: l_mono_c2l

    ! dimension of stencil for interpolation from cell
    ! centers to lon-lat points:
    INTEGER ::   rbf_dim_c2l
                                                                                 
    TYPE(t_lsq_set) :: lsq_lin_set, &! Parameter settings for linear and higher order  
      &                lsq_high_set  ! least squares

    ! Flag. If .FALSE. barycentric interpolation is replaced by a
    ! fallback interpolation.
    LOGICAL :: support_baryctr_intp

    ! Flag. If .TRUE. then the nest boundary points are taken out from
    ! the lat-lon interpolation stencil.
    LOGICAL :: lreduced_nestbdry_stencil
  
  !END TYPE t_interpol_config
  !>
  !!
  !TYPE(t_interpol_config) :: interpol_config(max_dom)


CONTAINS
  !>
  !!
  SUBROUTINE configure_interpolation( n_dom, grid_level, geometry_info )

    INTEGER,INTENT(IN) :: n_dom
    INTEGER,INTENT(IN) :: grid_level(n_dom)
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info(:)
    
    INTEGER :: jg, jlev, jlev_shift, geometry_type
    CHARACTER(len=*),PARAMETER :: routine = 'mo_interpol_config:configure_interpol'
    REAL(wp) :: resol

    !-----------------------------------------------------------------------
    ! Set stencil size for RBF vector reconstruction
    !-----------------------------------------------------------------------
    rbf_vec_dim_c   = 9    ! use 2nd order reconstruction for cell centers
    rbf_vec_dim_v   = 6    ! use 6-point reconstruction for vertices
    rbf_vec_dim_e   = 4    ! use 4-point reconstruction for edge midpoints
    rbf_c2grad_dim  = 10   ! for reconstructing gradient at cell midpoints

    !-----------------------------------------------------------------------
    ! If RBF scaling factors are not supplied by the namelist, they are now
    ! initialized with meaningful values depending on the characteristic grid 
    ! size.
    ! Please note: RBF scaling factors for p_patch(0) (if it exists)
    ! are not set here - they are taken from p_patch(1) in the setup routines
    ! Please note: The following RBF scaling factors are also valid for 
    ! simulations in which a scaling is applied to the earth radius 
    ! (a.k.a. small earth simulations). The characteristic grid size which is 
    ! used below takes account of the scaling. 
    !-----------------------------------------------------------------------


    !-----------------
    ! rbf_vec_scale_c 
    !-----------------
    ! - values are specified for Gaussian kernel
    ! (need to be smaller for inv. multiquadric)
    DO jg = 1,n_dom

      ! Check if scale factor is set in the namelist
      IF (rbf_vec_scale_c(jg) > 0.0_wp) CYCLE

      resol = geometry_info(jg)%mean_characteristic_length/1000._wp  ! resolution in km
      IF (resol >= 2.5_wp) THEN 
        rbf_vec_scale_c(jg) = 0.5_wp
      ELSE
        rbf_vec_scale_c(jg) = 0.5_wp/(1._wp+1.8_wp*LOG(2.5_wp/resol)**3.75)
      ENDIF
      IF (resol <= 0.125_wp) rbf_vec_scale_c(jg) = rbf_vec_scale_c(jg)*(resol/0.125_wp)**0.9_wp

    ENDDO

    !-----------------
    ! rbf_vec_scale_v
    !-----------------
    DO jg = 1, n_dom

      ! Check if scale factor is set in the namelist
      IF (rbf_vec_scale_v(jg) > 0.0_wp) CYCLE

      resol = geometry_info(jg)%mean_characteristic_length/1000._wp  ! resolution in km
      IF (resol >= 2._wp) THEN 
        rbf_vec_scale_v(jg) = 0.5_wp
      ELSE
        rbf_vec_scale_v(jg) = 0.5_wp/(1._wp+1.8_wp*LOG(2._wp/resol)**3)
      ENDIF
      IF (resol <= 0.125_wp) rbf_vec_scale_v(jg) = rbf_vec_scale_v(jg)*(resol/0.125_wp)**0.96_wp

    ENDDO

    !-----------------
    ! rbf_vec_scale_e 
    !-----------------
    ! - values are specified for inverse multiquadric kernel

    DO jg = 1, n_dom

      ! Check if scale factor is set in the namelist
      IF (rbf_vec_scale_e(jg) > 0.0_wp) CYCLE

      resol = geometry_info(jg)%mean_characteristic_length/1000._wp  ! resolution in km
      IF (resol >= 2._wp) THEN 
        rbf_vec_scale_e(jg) = 0.5_wp
      ELSE
        rbf_vec_scale_e(jg) = 0.5_wp/(1._wp+0.4_wp*LOG(2._wp/resol)**2)
      ENDIF
      IF (resol <= 0.125_wp) rbf_vec_scale_e(jg) = rbf_vec_scale_e(jg)*(resol/0.125_wp)**0.325_wp

    ENDDO

    IF (msg_level >= 7) THEN
      DO jg = 1, n_dom
        WRITE(message_text,'(a,i3,a,3f8.5)') 'RBF scale factors for cells/edges/vertices, domain ',jg,': ',&
         rbf_vec_scale_c(jg),rbf_vec_scale_e(jg),rbf_vec_scale_v(jg)
        CALL message('', TRIM(message_text))
      ENDDO
    ENDIF

    !-----------------
    ! rbf_vec_scale_ll 
    !-----------------
    ! - values are specified for Gaussian kernel
    ! (need to be smaller for inv. multiquadric)

    ! Computation of RBF scaling factors for lon-lat interpolation starts from the grid level
    ! grid_level(jg) known from the grid file, and shifts the level in integer steps
    ! according to the scaling applied to the sphere radius:
    !
    !   grid_sphere_radius = grid_sphere_radius(from file)*grid_rescale_factor(from grid_nml)
    !
    ! Thus scaling the RnBm grid by a grid_rescale_factor=0.125, for instance, for a planet
    ! 8 times smaller than Earth, we get the same resolution as on an unscaled the RnB(m+3)
    ! grid and should USE the same RBF scaling factors as for the grid level m+3.
    !
    ! Practically this shift in jlev can be computed as:
    !   jlev_shift = the nearest INTEGER of -Log_2(grid_rescale_factor)
    !
    ! Here how jlev is shifted for some values of grid_rescale_factor:
    !   grid_rescale_factor : jlev_shift
    !   1                   :  0
    !   1/2                 : +1
    !   1/3, 1/4, 1/5       : +2
    !   1/6, ..., 1/11      : +3

    jlev_shift = - NINT( LOG(grid_rescale_factor)/ln2 )

    rbf_vec_scale_ll(:) = -1.0_wp
    DO jg = 1,n_dom
       
      jlev = grid_level(jg) + jlev_shift
        
      IF      (jlev <= 6 ) THEN ; rbf_vec_scale_ll(jg) = 0.5_wp
      ELSE IF (jlev == 7 ) THEN ; rbf_vec_scale_ll(jg) = 0.35_wp
      ELSE IF (jlev == 8 ) THEN ; rbf_vec_scale_ll(jg) = 0.20_wp
      ELSE IF (jlev == 9 ) THEN ; rbf_vec_scale_ll(jg) = 0.05_wp
      ELSE IF (jlev == 10) THEN ; rbf_vec_scale_ll(jg) = 0.02_wp
      ELSE IF (jlev == 11) THEN ; rbf_vec_scale_ll(jg) = 0.0075_wp
      ELSE IF (jlev == 12) THEN ; rbf_vec_scale_ll(jg) = 0.0025_wp
      ELSE IF (jlev == 13) THEN ; rbf_vec_scale_ll(jg) = 0.001_wp
      ELSE                      ; rbf_vec_scale_ll(jg) = 0.0005_wp
      ENDIF

    ENDDO

    !AD (20 Sept 2013) Modification required for planar torus grid: the scale factor
    !is based on the width of the Gaussian but very soon it will adapted in more
    !analytical manner by Florian
    !AD (19 Oct 2014) No update from Florian yet: For exp(0.5r^2/r0^2), the r0 (scale_factor) 
    !should be greater than the minimal distance between points, and rather less than the 
    !maximal distance. For now, for the default 2nd order cases, using r0=dual_edge_length 
    !for all scales but tesing is required

    DO jg = 1, n_dom
      geometry_type = geometry_info(jg)%geometry_type

      IF( geometry_type==planar_torus_geometry ) THEN
        rbf_vec_scale_c(jg)  = geometry_info(jg)%mean_dual_edge_length
        rbf_vec_scale_e(jg)  = rbf_vec_scale_c(jg)
        rbf_vec_scale_v(jg)  = rbf_vec_scale_c(jg)
        rbf_vec_scale_ll(jg) = rbf_vec_scale_c(jg)
        CALL message( TRIM(routine),'Modifying rbf_vec_scale for torus grid: ignore warnings!')
      END IF
    ENDDO

    !-----------------------------------------------------------------------
    ! Now check the RBF scaling factors
    !-----------------------------------------------------------------------
    DO jg = 1, n_dom
      IF (rbf_vec_scale_c(jg) < 1.e-10_wp) THEN
        CALL finish( TRIM(routine),'wrong value of rbf_vec_scale_c')
      ELSE IF ((rbf_vec_scale_c(jg) < 0.01_wp).OR.(rbf_vec_scale_c(jg) > 0.6_wp)) THEN
        CALL message( TRIM(routine),'WARNING: ')
        CALL message('','! recommended range for rbf_vec_scale_c')
        CALL message('','! is 0.01 <= rbf_vec_scale_c <= 0.6')
      ENDIF
    ENDDO

    DO jg = 1, n_dom
      IF (rbf_vec_scale_v(jg) < 1.e-10_wp) THEN
        CALL finish( TRIM(routine),'wrong value of rbf_vec_scale_v')
      ELSE IF ((rbf_vec_scale_v(jg) < 0.02_wp).OR.(rbf_vec_scale_v(jg) > 1.0_wp)) THEN
        CALL message( TRIM(routine),'WARNING: ')
        CALL message('','! recommended range for rbf_vec_scale_v')
        CALL message('','! is 0.02 <= rbf_vec_scale_v <= 1.0')
      ENDIF
    ENDDO
   
    DO jg = 1, n_dom
      IF (rbf_vec_scale_e(jg) < 1.e-10_wp) THEN
        CALL finish( TRIM(routine),'wrong value of rbf_vec_scale_e')
      ELSE IF ((rbf_vec_scale_e(jg) < 0.05_wp).OR.(rbf_vec_scale_e(jg) > 0.6_wp)) THEN
        CALL message( TRIM(routine),'WARNING: ')
        CALL message('','! recommended range for rbf_vec_scale_e')
        CALL message('','! is 0.05 <= rbf_vec_scale_e <= 0.6')
      ENDIF
    ENDDO


    !-----------------------------------------------------------------------
    ! set the number of unknowns in the least squares reconstruction, and the
    ! stencil size, depending on the chosen polynomial order (lsq_high_ord).
    ! Default value for lsq_high_ord is 3. The user might have made a different
    ! choice. The possibilities are:
    !  lsq_high_ord=1 : linear polynomial           : 2 unknowns with a 3-point stencil
    !  lsq_high_ord=2 : quadratic polynomial        : 5 unknowns with a 9-point stencil
    !  lsq_high_ord=3 : cubic polynomial            : 9 unknowns with a 9-point stencil
    !-----------------------------------------------------------------------

    ! Settings for linear lsq reconstruction

    lsq_lin_set%l_consv = llsq_lin_consv
    lsq_lin_set%dim_c   = 3
    lsq_lin_set%dim_unk = 2
    IF ( lsq_lin_set%l_consv ) THEN
      lsq_lin_set%wgt_exp = 0
    ELSE
      lsq_lin_set%wgt_exp = 2
    ENDIF

    ! Settings for high order lsq reconstruction

    lsq_high_set%l_consv = llsq_high_consv

    IF (lsq_high_ord == 1) THEN
      lsq_high_set%dim_c   = 3
      lsq_high_set%dim_unk = 2
      IF ( lsq_high_set%l_consv ) THEN
        lsq_high_set%wgt_exp = 0
      ELSE
        lsq_high_set%wgt_exp = 2
      ENDIF      
    ELSE IF (lsq_high_ord == 2) THEN
      lsq_high_set%dim_c   = 9
      lsq_high_set%dim_unk = 5
      lsq_high_set%wgt_exp = 3
    ELSE IF (lsq_high_ord == 3) THEN
      lsq_high_set%dim_c   = 9
      lsq_high_set%dim_unk = 9
      lsq_high_set%wgt_exp = 0
    ELSE
      CALL finish( TRIM(routine),'wrong value of lsq_high_ord, must be 1,2, or 3')
    ENDIF

  END SUBROUTINE configure_interpolation


END MODULE mo_interpol_config
