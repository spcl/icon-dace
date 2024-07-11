!
! This module contains types and subroutines used for
! the extrapolation of initial data to the upper atmosphere.
! The extrapolation uses a blending of two states:
! - "left" state: a linear extrapolation of the IFS-input-data from below
!   the extrapolation region (already computed in 'src/atm_dyn_iconam/mo_nh_vert_interp')
! - "right" state: a state derived from a simple climatology covering
!   the upper atmosphere
! The blending of the two states is in such a way that:
! - the lower in the extrapolation region, the more the left state dominates
! - the higher in the extrapolation region, the more the right state dominates
! (Note: many of the data types and subroutines in this module
! appear computationally inefficient for the current form of
! upper-atmosphere extrapolation, but might be regarded as rudimentary
! modularity basis for sophisticated extrapolations in the future.
! In addition, the extrapolation routines are expected to be called
! only once during model initialization, so that some runtime/memory
! overhead seems to be bearable.)
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nh_vert_extrap_utils

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish, message, message_text
  USE mo_run_config,             ONLY: msg_level
  USE mo_grid_config,            ONLY: grid_sphere_radius, n_dom_start
  USE mo_dynamics_config,        ONLY: ldeepatmo
  USE mo_upatmo_config,          ONLY: upatmo_exp_config, upatmo_config
  USE mo_parallel_config,        ONLY: nproma
  USE mo_vertical_coord_table,   ONLY: vct_a
  USE mo_model_domain,           ONLY: t_patch
  USE mo_nonhydro_types,         ONLY: t_nh_metrics
  USE mo_intp_data_strc,         ONLY: t_int_state
  USE mo_impl_constants,         ONLY: MAX_CHAR_LENGTH, SUCCESS,      & 
    &                                  min_rlcell, min_rledge,        &
    &                                  min_rlcell_int, min_rledge_int
  USE mo_impl_constants_grf,     ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_math_constants,         ONLY: pi, deg2rad, rad2deg
  USE mo_physical_constants,     ONLY: cpd, grav, rd_o_cpd, p0ref
  USE mo_loopindices,            ONLY: get_indices_c, get_indices_e
  USE mo_intp,                   ONLY: cells2edges_scalar, cells2verts_scalar
  USE mo_math_gradients,         ONLY: grad_fd_tang
  USE mo_timer,                  ONLY: timers_level, timer_start, timer_stop, & 
    &                                  timer_expol
  USE mo_fortran_tools,          ONLY: init
  USE mo_util_string,            ONLY: int2string, real2string
  USE mo_mpi,                    ONLY: get_my_mpi_work_id,           &
    &                                  get_my_mpi_work_communicator, &
    &                                  p_bcast
  USE mo_sync,                   ONLY: global_max
  USE mo_upatmo_impl_const,      ONLY: iUpatmoStat, imsg_thr, itmr_thr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_expol_state
  PUBLIC :: sanity_check

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nh_vert_extrap_utils'

  ! Number of vertical sampling points of the climatology (-> 't_clim_rawdata')
  INTEGER, PARAMETER :: nclim = 25

  !-----------------------------------------------------------------
  !                          Data types
  !-----------------------------------------------------------------

  ! Raw data of climatology.
  ! The climatology used for the reference temperature sampling points is from: 
  ! Fleming, E.L., Chandra, S., Shoeberl, M.R. & Barnett, J.J. (1988), 
  ! "Monthly mean global climatology of temperature, wind, geopotential height 
  ! and pressure for 0 - 120 km", NASA-TM-100697
  TYPE :: t_clim_rawdata
    REAL(wp), DIMENSION(nclim) :: zgpot  ! Geopotential height
    REAL(wp), DIMENSION(nclim) :: temp   ! Temperature
  END TYPE t_clim_rawdata
  TYPE(t_clim_rawdata), PARAMETER :: clim_rawdata = t_clim_rawdata(      & 
    ! Geopotential heights 'zgpot' (vertical ordering is from top to bottom)
    & (/ 120000., 115000., 110000., 105000., 100000.,  95000.,  90000.,  & 
    &     85000.,  80000.,  75000.,  70000.,  65000.,  60000.,  55000.,  &
    &     50000.,  45000.,  40000.,  35000.,  30000.,  25000.,  20000.,  & 
    &     15000.,  10000.,   5000.,      0. /),                          &
    ! Temperatures 'temp' at geopotential heigts
    & (/  380.72,  303.27,  238.61,  204.59,  188.03,  182.04,  185.86,  & 
    &     194.48,  204.25,  211.46,  218.55,  230.29,  243.69,  256.88,  & 
    &     266.15,  265.30,  255.11,  241.45,  228.04,  219.63,  212.40,  &
    &     208.92,  230.32,  261.92,  285.85 /)                           &
    &                                                                    )

  !-----------------------------------------------------------------------------

  ! Indices for blending factor (just for convenience)
  TYPE :: t_iblnd
    INTEGER :: left  ! Factor for 'left' blending part
    INTEGER :: right ! Factor for 'right' blending part
    INTEGER :: nitem ! Number of blending parts (for allocation, 
                     ! loops etc., assuming this list to start 
                     ! always with index 1!)
  END TYPE t_iblnd
  TYPE(t_iblnd), PARAMETER :: iblnd = t_iblnd( 1,  & !left
    &                                          2,  & !right
    &                                          2   ) !nitem

  !-----------------------------------------------------------------------------

  ! Variables which could be useful for construction, initialization etc.
  TYPE :: t_constructor_kit
    INTEGER :: nlev
    INTEGER :: nblks_c
    INTEGER :: nblks_e
    INTEGER :: nblnd
  END TYPE t_constructor_kit

  !-----------------------------------------------------------------------------

  ! Loop variables for hand-over in subroutines 'copy_state', 'blend_states' 
  ! and 'multiply_state'
  TYPE :: t_loop_kit
    INTEGER :: nlev
    INTEGER :: nblks
    INTEGER :: npromz
  END TYPE t_loop_kit

  !-----------------------------------------------------------------------------

  ! See 'src/atm_dyn_iconam/mo_initicon_types' for templates of the following data types

  ! Data type containing vertical climatoligical profiles
  TYPE :: t_clim_state
    
    ! Flag should be .true., if this data structure has been allocated
    LOGICAL :: linitialized = .FALSE.
    
    ! In the current version the upper-atmosphere extrapolation 
    ! is based on one vertical temperature profile
    REAL(wp), POINTER, DIMENSION(:) :: temp  => NULL()
    
  CONTAINS
    PROCEDURE :: construct => t_clim_state_construct 
    PROCEDURE :: finalize  => t_clim_state_finalize  
  END TYPE t_clim_state

  !-----------------------------------------------------------------------------

  ! Data type containing variables related to the geopotential height
  TYPE :: t_expol_metrics

    ! 'expol_start_height' in terms of geopotential height
    REAL(wp) :: zgpot_start

    ! The upper-atmosphere extrapolation takes place above 'flat_height', 
    ! so the following fields need to vary only in vertical direction
    REAL(wp), POINTER, DIMENSION(:)   :: zgpot_mc      => NULL(), &
      &                                  zgpot_ifc     => NULL()
    REAL(wp), POINTER, DIMENSION(:,:) :: wfac_blnd_mc  => NULL(), & 
      &                                  wfac_blnd_ifc => NULL()

    ! Flag should be .true., if this data structure has been allocated
    LOGICAL :: linitialized = .FALSE.

  CONTAINS
    PROCEDURE :: construct => t_expol_metrics_construct 
    PROCEDURE :: finalize  => t_expol_metrics_finalize  
  END TYPE t_expol_metrics

  !-----------------------------------------------------------------------------

  ! For the left and right blending states of the upper-atmosphere extrapolation
  ! (Consider: maybe introduce a more compact treatment of the blending states, 
  ! e.g. as with the water phases: 
  ! expol_blending_state%temp(:,:,:) -> expol_blending_state%var(:,:,:,ivar%temp)?)
  TYPE :: t_expol_blending_state

    REAL(wp), POINTER, DIMENSION(:,:,:) :: temp  => NULL(),  &
      &                                    vn    => NULL(),  &
      &                                    w     => NULL(),  &
      &                                    qc    => NULL(),  &
      &                                    qi    => NULL(),  &
      &                                    qr    => NULL(),  &
      &                                    qs    => NULL()

    ! Note: in the current form of the extrapolation pressure 'pres' 
    ! and water vapour 'qv' are not required

    ! Flag should be .true., if this data structure has been allocated
    LOGICAL :: linitialized = .FALSE.

  CONTAINS
    PROCEDURE :: construct => t_expol_blending_state_construct  
    PROCEDURE :: finalize  => t_expol_blending_state_finalize   
  END TYPE t_expol_blending_state

  !-----------------------------------------------------------------------------

  ! Container for the above constructions
  TYPE :: t_expol_state

    ! Access to the following data is not necessary outside this module
    PRIVATE

    TYPE(t_clim_state)           :: clim
    TYPE(t_expol_metrics)        :: metrics
    TYPE(t_expol_blending_state) :: left    ! (-> Upward linear extrapolation of IFS-data)
    TYPE(t_expol_blending_state) :: right   ! (-> From a climatology covering the upper atmosphere)

    ! Flag should be .true., if timer control should be switched on
    LOGICAL :: ltimer = .FALSE.
    ! Flag should be .true., if the current simulation makes use 
    ! of the limited-area I/O infrastructure
    LOGICAL :: latbcmode = .FALSE.
    ! Flag should be .true., if this data structure has been allocated
    LOGICAL :: linitialized = .FALSE.

  CONTAINS
    ! The following subroutines should be accessible outside this module, 
    ! and are public by default (if the derived type is public)
    PROCEDURE :: initialize => t_expol_state_initialize
    PROCEDURE :: temp       => t_expol_state_temp
    PROCEDURE :: pres       => t_expol_state_pres
    PROCEDURE :: vn         => t_expol_state_vn
    PROCEDURE :: w          => t_expol_state_w
    PROCEDURE :: qx         => t_expol_state_qx
    PROCEDURE :: finalize   => t_expol_state_finalize
  END TYPE t_expol_state

  !-----------------------------------------------------------------------------

CONTAINS  !..................................................................................

  !-----------------------------------------------------------------
  !                          Subroutines
  !-----------------------------------------------------------------

  !---------------------------------------------------------------
  !             Routines of container t_expol_state
  !---------------------------------------------------------------

  !>
  !! This subroutine initializes quantities 
  !! necessary for the extrapolation of IFS-data  
  !! to the upper atmosphere. 
  !!
  SUBROUTINE t_expol_state_initialize ( expolstate,   &  !class
    &                                   p_patch,      &  !in
    &                                   opt_latbcmode ) !optin

    CLASS(t_expol_state), INTENT(INOUT) :: expolstate

    ! In/out variables
    TYPE(t_patch),     INTENT(IN) :: p_patch        ! Domain properties
    LOGICAL, OPTIONAL, INTENT(IN) :: opt_latbcmode  ! Switch to indicate limited-area mode

    ! Local variables
    TYPE(t_constructor_kit) :: cnstr
    INTEGER :: jg, jk, jb, jc, nexpollev
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx 
    INTEGER :: rl_start, rl_end    

    CHARACTER(LEN=*), PARAMETER :: &
      &       routine = modname//':t_expol_state_initialize'

    !----------------------------------------------

    ! (The domain loop lies outside the call of this subroutine)

    ! (Note: if 'upatmo_config(jg)%exp%l_expol==.TRUE.' will be checked outside 
    ! the call of this subroutine for efficiency reasons)

    ! (The timer is started inside instead of outside this subroutine, 
    ! in order to reduce the amount of upper-atmosphere-related code 
    ! in 'src/atm_dyn_iconam/mo_nh_vert_interp')
    expolstate%ltimer = (timers_level > itmr_thr%med)
    IF (expolstate%ltimer) CALL timer_start(timer_expol)

    ! Domain index
    jg = p_patch%id 

    ! Some checks
    IF (.NOT. upatmo_config(jg)%l_status(iUpatmoStat%configured)) THEN
      ! Upper-atmosphere info should be available
      CALL finish(routine, "Check calling sequence: upatmo_config is not configured.")
    ELSEIF (jg < 1) THEN  
      ! This subroutine is only meant for 'jg>=1'
      WRITE (message_text, '(a,i0,a)') "Dom ", jg, " is not supported."
      CALL finish(routine, message_text)
    ENDIF

    ! Print some info                
    IF (msg_level >= imsg_thr%high) THEN
      WRITE(message_text,'(a,i4)') 'Start extrapolation to upper atmosphere in domain: ', jg
      CALL message(routine, message_text)
    ENDIF
    
    ! Initialize extrapolation data types

    ! Number of layers in extrapolation region
    nexpollev = upatmo_config(jg)%exp%nexpollev

    ! Collect some data necessary for construction and initialization
    cnstr%nlev    = nexpollev
    cnstr%nblks_c = p_patch%nblks_c
    cnstr%nblks_e = p_patch%nblks_e
    cnstr%nblnd   = iblnd%nitem
    
    ! Construct the fields contained in 'expolstate'
    IF (.NOT. expolstate%linitialized) THEN 
      CALL expolstate%clim%construct(cnstr) 
      CALL expolstate%metrics%construct(cnstr) 
      CALL expolstate%left%construct(cnstr) 
      CALL expolstate%right%construct(cnstr) 
      expolstate%linitialized = .TRUE.
    ELSE
      CALL finish(routine, 'Attempt to initialize while already/still initialized')
    ENDIF

    ! Is the limited-area mode active?
    IF (PRESENT(opt_latbcmode)) THEN
      expolstate%latbcmode = opt_latbcmode
    ELSE
      expolstate%latbcmode = .FALSE.
    ENDIF

    ! Compute variables related to the geopotential height
    CALL init_expol_metrics( expolstate%metrics, &  !inout
      &                      p_patch             )  !in
    
    ! Compute climatological temperature profile on 
    ! grid level heights
    CALL init_clim_state( expolstate%clim,    &  !inout
      &                   expolstate%metrics, &  !in
      &                   p_patch             )  !in 
    
    ! Some of the left-hand side and right-hand side states 
    ! for the extrapolation blending can be assigned already here 
    
!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)
    
    rl_start   = 1
    rl_end     = min_rlcell
    i_startblk = p_patch%cells%start_block(rl_start) 
    i_endblk   = p_patch%cells%end_block(rl_end)  
    
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
      DO jk = 1, nexpollev
        DO jc = i_startidx, i_endidx
          ! Climatological right-hand side temperature state
          expolstate%right%temp(jc,jk,jb) = expolstate%clim%temp(jk) 
          ! Climatological right-hand side vertical wind is assumed to be ~zero
          expolstate%right%w(jc,jk,jb) = 0._wp
          ! In the current simple form of the extrapolation, all water phases except for 
          ! water vapour are assumed to be zero in the extrapolation region.
          ! (The current treatment of the 'qx' is not computationally efficient, 
          ! but it is done in this way to follow the general procedure)
          expolstate%left%qc(jc,jk,jb)  = 0._wp    ! Cloud water
          expolstate%left%qi(jc,jk,jb)  = 0._wp    ! Ice
          expolstate%left%qr(jc,jk,jb)  = 0._wp    ! Rain water
          expolstate%left%qs(jc,jk,jb)  = 0._wp    ! Snow
          expolstate%right%qc(jc,jk,jb) = 0._wp
          expolstate%right%qi(jc,jk,jb) = 0._wp
          expolstate%right%qr(jc,jk,jb) = 0._wp
          expolstate%right%qs(jc,jk,jb) = 0._wp
        ENDDO  !jc
      ENDDO  !jk
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL
    
    IF (expolstate%ltimer) CALL timer_stop(timer_expol)
    
  END SUBROUTINE t_expol_state_initialize

  !-----------------------------------------------------------------------------

  !>
  !! This subroutine extrapolates the temperature 
  !! to the upper atmosphere via blending of two states.
  !!
  SUBROUTINE t_expol_state_temp ( expolstate, &  !class
    &                             p_patch,    &  !in
    &                             temp_out    )  !inout

    CLASS(t_expol_state), INTENT(INOUT) :: expolstate

    ! In/out variables
    TYPE(t_patch), INTENT(IN)    :: p_patch          ! Domain properties    
    REAL(wp),      INTENT(INOUT) :: temp_out(:,:,:)  ! Input/output temperature field 

    ! Local variables
    REAL(wp), ALLOCATABLE :: dtempdzgpot(:,:,:)

    TYPE(t_loop_kit) :: loopkit

    INTEGER :: jg, jb, jc, jk, nexpollev, nblks_c
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end    
    INTEGER :: istat
    LOGICAL :: lpassed
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: sanity_check_msg

    REAL(wp), PARAMETER :: ngrav_o_cpd = -grav / cpd
    CHARACTER(LEN=*), PARAMETER :: &
      &       routine = modname//':t_expol_state_temp'

    !----------------------------------------------

    ! (Note: if 'upatmo_config(jg)%exp%l_expol==.TRUE.' will be checked outside 
    ! the call of this subroutine for efficiency reasons)

    IF (expolstate%ltimer) CALL timer_start(timer_expol)

    ! Domain index
    jg = p_patch%id

    ! Number of layers in extrapolation region
    nexpollev = upatmo_config(jg)%exp%nexpollev

    ! Number of blocks
    nblks_c = p_patch%nblks_c

    ! The extrapolated temperature field is a blending of the form:
    ! extrapolated temperature = wfac_left * (temperature linearly extrapolated 
    !                            from IFS-data) + wfac_right * (temperature from 
    !                            climatology covering the upper atmosphere)

    ! For the left-hand side state, we take the temperature field already 
    ! linearly extrapolated in 'src/atm_dyn_iconam/mo_nh_vert_interp', 
    ! this way the upper-atmosphere extrapolation can be kept relatively simple here
    loopkit%nlev   = nexpollev
    loopkit%nblks  = nblks_c
    loopkit%npromz = p_patch%npromz_c
    CALL copy_state( temp_out,             & !in
      &              expolstate%left%temp, & !out
      &              loopkit               ) !in

    ! The climatological right-hand side state has already 
    ! been assigned in 't_expol2upatmo_initialize', so we can advance to the blending
    CALL blend_states( expolstate%left%temp,            & !in (LEFT!)
      &                expolstate%right%temp,           & !in (RIGHT!)
      &                temp_out,                        & !inout
      &                expolstate%metrics%wfac_blnd_mc, & !in 
      &                loopkit                          ) !in

    ! ('sanity_check' makes use of MPI communication ('global_max'), 
    ! so it might be better to skip this for the limited-area mode 
    ! (see arguments in 'src/atm_dyn_iconam/mo_nh_vert_interp: vert_interp'))
    IF (upatmo_exp_config(jg)%lexpol_sanitycheck .AND. .NOT. expolstate%latbcmode) THEN

      ! Apply some sanity checks to the extrapolated temperature
      ! 1st Is temperature positive-definite?
      CALL sanity_check ( p_patch     = p_patch,         &  !in
        &                 state       = temp_out,        &  !in
        &                 bound       = 0._wp,           &  !in
        &                 keys        = 'cell,lower',    &  !in
        &                 lpassed     = lpassed,         &  !out
        &                 opt_slev    = 1,               &  !optin
        &                 opt_elev    = nexpollev,       &  !optin
        &                 opt_message = sanity_check_msg )  !optout

      IF (.NOT. lpassed) THEN
        WRITE(message_text,'(3a)') 'Extrapolated temperature is not positive-definite. ', &
          & TRIM(sanity_check_msg), ' A smaller value for expol_blending_scale might help.'
        CALL finish(routine, message_text)
      ENDIF

      ! 2nd Is vertical temperature gradient greater than 
      ! the negative dry adiabatic lapse rate -grav/cpd?
      ! (=> Prevent the extrapolated state from being convectively unstable)
      IF (nexpollev > 1) THEN
        ALLOCATE( dtempdzgpot(nproma, nexpollev-1, nblks_c), STAT=istat )
        IF (istat /= SUCCESS) CALL finish(routine, 'Allocation of dtempdzgpot failed')
        ! Prognostic domain only
        rl_start   = grf_bdywidth_c + 1
        rl_end     = min_rlcell_int
        i_startblk = p_patch%cells%start_block(rl_start) 
        i_endblk   = p_patch%cells%end_block(rl_end)  
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
          ! (To be on the safe side, we do not compute the gradient across 'expol_start_height', 
          ! which would require: 'DO jk=1, nexpollev')
          DO jk=1, nexpollev-1
            DO jc = i_startidx, i_endidx
              dtempdzgpot(jc,jk,jb) = ( temp_out(jc,jk,jb) - temp_out(jc,jk+1,jb) ) / &
                &                     ( expolstate%metrics%zgpot_mc(jk) - expolstate%metrics%zgpot_mc(jk+1) )
            ENDDO  !jc
          ENDDO  !jk
        ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL
        CALL sanity_check ( p_patch     = p_patch,         &  !in
          &                 state       = dtempdzgpot,     &  !in
          &                 bound       = ngrav_o_cpd,     &  !in
          &                 keys        = 'cell,lower',    &  !in
          &                 lpassed     = lpassed,         &  !out
          &                 opt_slev    = 1,               &  !optin
          &                 opt_elev    = nexpollev-1,     &  !optin
          &                 opt_message = sanity_check_msg )  !optout        
        IF (.NOT. lpassed) THEN
          WRITE(message_text,'(3a)') 'Extrapolated state is convectively unstable. ', &
            & TRIM(sanity_check_msg), ' A smaller value for expol_blending_scale might help.'
          CALL finish(routine, message_text)
        ENDIF
        DEALLOCATE( dtempdzgpot, STAT=istat )
        IF (istat /= SUCCESS) CALL finish(routine, 'Deallocation of dtempdzgpot failed')
      ENDIF  !IF (nexpollev > 1)

    ENDIF  !IF (lexpol_sanitycheck)

    IF (expolstate%ltimer) CALL timer_stop(timer_expol)

  END SUBROUTINE t_expol_state_temp

  !-----------------------------------------------------------------------------

  !>
  !! This subroutine extrapolates the water phases 
  !! except for water vapour to the upper atmosphere 
  !! via blending of two states.
  !!
  SUBROUTINE t_expol_state_qx ( expolstate, &  !class
    &                           p_patch,    &  !in
    &                           qc_out,     &  !inout
    &                           qi_out,     &  !inout
    &                           qr_out,     &  !inout
    &                           qs_out      )  !inout

    CLASS(t_expol_state), INTENT(INOUT) :: expolstate

    ! In/out variables
    TYPE(t_patch),    INTENT(IN)    :: p_patch        ! Domain properties    
    REAL(wp), TARGET, INTENT(INOUT) :: qc_out(:,:,:)  ! Cloud water
    REAL(wp), TARGET, INTENT(INOUT) :: qi_out(:,:,:)  ! Ice
    REAL(wp), TARGET, INTENT(INOUT) :: qr_out(:,:,:)  ! Rain water
    REAL(wp), TARGET, INTENT(INOUT) :: qs_out(:,:,:)  ! Snow   

    ! Local variables
    TYPE(t_loop_kit) :: loopkit

    INTEGER :: jg

    CHARACTER(LEN=*), PARAMETER :: &
      &       routine = modname//':t_expol_state_qx'

    !----------------------------------------------

    ! (Note: if 'upatmo_config(jg)%exp%l_expol==.TRUE.' will be checked outside 
    ! the call of this subroutine for efficiency reasons)

    IF (expolstate%ltimer) CALL timer_start(timer_expol)

    ! Domain index
    jg = p_patch%id

    ! Note: water vapour 'qv' is not touched by this routine, 
    ! it is completely processed in 'src/atm_dyn_iconam/mo_nh_vert_interp'

    ! In the current form of the upper-atmosphere extrapolation, 
    ! 'qc', 'qi', 'qr', and 'qs' are simply set to zero, the corresponding 
    ! left-hand side and right-hand side states have already been assigned 
    ! in 't_expol_state_initialize', so we can advance to the blending 
    ! (-> not computationally efficient, but it shall follow the general 
    ! procedure, in order to maybe potentially support to some extent sophisticated 
    ! extrapolations in the future)
    ! Note: if some of the water phases shall have another value than zero 
    ! in the upper atmosphere someday, be aware of their post-interpolation-treatment  
    ! in 'src/atm_dyn_iconam/mo_initicon_utils/copy_initicon2prog_atm'!
    
    ! Collect loop variables
    loopkit%nlev   = upatmo_config(jg)%exp%nexpollev
    loopkit%nblks  = p_patch%nblks_c
    loopkit%npromz = p_patch%npromz_c

    ! (The following could become a loop over the water phases, 
    ! with pointers associated to the corresponding fields, 
    ! but this seems to be not worth the effort and computationally less efficient 
    ! than simply calling 'blend_state' explicitly for each water phase)
    ! Cloud water:
    CALL blend_states( expolstate%left%qc,              & !in (LEFT!)
      &                expolstate%right%qc,             & !in (RIGHT!)
      &                qc_out,                          & !inout
      &                expolstate%metrics%wfac_blnd_mc, & !in 
      &                loopkit                          ) !in
    ! Ice:
    CALL blend_states( expolstate%left%qi, expolstate%right%qi, qi_out, & 
      &                expolstate%metrics%wfac_blnd_mc, loopkit         ) 
    ! Rain water:
    CALL blend_states( expolstate%left%qr, expolstate%right%qr, qr_out, & 
      &                expolstate%metrics%wfac_blnd_mc, loopkit         ) 
    ! Snow:
    CALL blend_states( expolstate%left%qs, expolstate%right%qs, qs_out, & 
      &                expolstate%metrics%wfac_blnd_mc, loopkit         ) 

    ! In its current form the extrapolation of tracers requires no sanity check

    IF (expolstate%ltimer) CALL timer_stop(timer_expol)

  END SUBROUTINE t_expol_state_qx

  !-----------------------------------------------------------------------------

  !>
  !! This subroutine "extrapolates" the pressure to the upper atmosphere
  !! via integration of the hydrostatic balance.
  !!
  SUBROUTINE t_expol_state_pres ( expolstate, &  !class
    &                             p_patch,    &  !in
    &                             press_out,  &  !inout
    &                             tempv_in,   &  !in
    &                             p_metrics,  &  !in
    &                             opt_slev,   &  !optional in
    &                             opt_elev    )  !optional in

    CLASS(t_expol_state), INTENT(INOUT) :: expolstate

    ! In/out variables
    TYPE(t_patch),       INTENT(IN)    :: p_patch          ! Domain properties    
    REAL(wp),            INTENT(INOUT) :: press_out(:,:,:) ! Input/output pressure field
    REAL(wp),            INTENT(IN)    :: tempv_in(:,:,:)  ! Input virtual temperature
    TYPE(t_nh_metrics),  INTENT(IN)    :: p_metrics        ! Model grid properties and background state
    ! Optional input variables
    INTEGER, INTENT(IN), OPTIONAL      :: opt_slev         ! Optional vertical start level
    INTEGER, INTENT(IN), OPTIONAL      :: opt_elev         ! Optional vertical end level    

    ! Local variables
    INTEGER  :: jg, jb, jk, jc, elev, slev 
    INTEGER  :: nlen, nlev, nexpollev
    REAL(wp) :: exner, fac1, fac2, fac3, a, b, c

    REAL(wp), PARAMETER :: cpd_o_rd  = 1._wp / rd_o_cpd
    CHARACTER(LEN=*), PARAMETER :: &
      &       routine = modname//':t_expol_state_pres'

    !----------------------------------------------

    ! (Note: if 'upatmo_config(jg)%exp%l_expol==.TRUE.' will be checked outside 
    ! the call of this subroutine for efficiency reasons)

    IF (expolstate%ltimer) CALL timer_start(timer_expol)

    ! Domain index
    jg = p_patch%id 

    ! Number of vertical grid layers
    nlev = p_patch%nlev

    ! Number of layers in extrapolation region
    nexpollev = upatmo_config(jg)%exp%nexpollev

    ! Check optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = nexpollev
    END IF

    ! Some rudimentary consistency check
    IF ( slev < 1    .OR. &
      &  elev > nlev .OR. & 
      &  slev > elev      ) CALL finish(routine, 'Invalid optional input levels')

    ! If 'elev < nlev' the integration can include 'elev' itself 
    ! (see start index of jk-loop below)
    IF ( elev < nlev ) elev = elev + 1

    ! Some notes: 
    ! - The hydrostatic balance in the deep atmosphere reads:
    ! 
    !   vn * (-vn/r + ft) + vt * (-vt/r - fn) = -cp * (thetav * dexner'/dr + thetav' * dexner0/dr), (I)
    !   -------------------------------------
    !
    !   where vn and vt are the edge-normal and -tangential wind components, respectively, 
    !   r is the distance from the center of the Earth, 
    !   fn and ft are the edge-normal and -tangential Coriolis parameters, respectively, 
    !   thetav and exner are the virtual potential temperature and the Exner pressure, respectively, 
    !   and (.)' indicate a deviation from a reference state (.)0.
    !   The underlined terms on the left-hand side are neglected, to simplify the computation.
    ! 
    ! - The vertical variation of the gravitational acceleration: 
    !
    !                                      grav -> grav * ( a / r )**2,                             (II)
    !
    !   where a denotes the radius of the model Earth, is implicitly included in dexner0/dr 
    !   (and in thatav0 and exner0) in eq. (I).
    !
    ! - Using thetav = tempv / exner, the discretization of eq. (I) reads:
    !
    !   [ wgt(jk) * tempv(jk) / exner(jk) + wgt(jk+1) * tempv(jk+1) / exner(jk+1) ] * 
    !   { [ exner(jk) - exner0(jk) ] - [ exner(jk+1) - exner0(jk+1) ] } / [ z(jk) - z(jk+1) ] 
    !   + { wgt(jk) * [ tempv(jk) / exner(jk) - thetav0(jk) ] 
    !   + wgt(jk+1) * [ tempv(jk+1) / exner(jk+1) - theta0(jk+1) ] } * dexner0/dr(jk+1/2) = 0,     (III)     
    !
    !   where wgt denotes a geometric weighting with wgt(jk) + wgt(jk+1) = 1, 
    !   and jk+1/2 indexes the interface between the full levels jk and jk+1. 
    !   Given tempv(jk), tempv(jk+1) and exner(jk+1), eq. (III) is solved for exner(jk).
    !
    ! - The computation is a copy of 'src/atm_dyn_iconam/mo_nh_init_utils/hydro_adjust'.

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,exner,fac1,fac2,fac3,a,b,c) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, p_patch%nblks_c
      IF (jb /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      ENDIF
      DO jk = elev-1, slev, -1
        DO jc = 1, nlen

          exner = ( press_out(jc,jk+1,jb) / p0ref )**rd_o_cpd

          fac1  = p_metrics%wgtfac_c(jc,jk+1,jb) * ( tempv_in(jc,jk+1,jb) &
            &   - p_metrics%theta_ref_mc(jc,jk+1,jb) * exner )            &
            &   - ( 1._wp - p_metrics%wgtfac_c(jc,jk+1,jb) )              &
            &   * p_metrics%theta_ref_mc(jc,jk,jb) * exner

          fac2  = ( 1._wp - p_metrics%wgtfac_c(jc,jk+1,jb) ) * tempv_in(jc,jk,jb) * exner

          fac3  = p_metrics%exner_ref_mc(jc,jk+1,jb) &
            &   - p_metrics%exner_ref_mc(jc,jk,jb) - exner

          a     = ( p_metrics%theta_ref_ic(jc,jk+1,jb) * exner + fac1 ) / &
            &     p_metrics%ddqz_z_half(jc,jk+1,jb)

          b     = -( a * fac3 + fac2 / p_metrics%ddqz_z_half(jc,jk+1,jb) &
            &   + fac1 * p_metrics%d_exner_dz_ref_ic(jc,jk+1,jb) )

          c     = -( fac2 * fac3 / p_metrics%ddqz_z_half(jc,jk+1,jb) &
            &   + fac2 * p_metrics%d_exner_dz_ref_ic(jc,jk+1,jb) )

          exner = ( b + SQRT( b**2 + 4._wp * a * c ) ) / ( 2._wp * a )

          press_out(jc,jk,jb) = p0ref * exner**cpd_o_rd

        ENDDO !jc
      ENDDO !jk
    ENDDO !jb
!$OMP END DO
!$OMP END PARALLEL

    ! No sanity check is applied to the pressure
    
    IF (expolstate%ltimer) CALL timer_stop(timer_expol)

  END SUBROUTINE t_expol_state_pres

  !-----------------------------------------------------------------------------

  !>
  !! This subroutine extrapolates the vertical wind  
  !! to the upper atmosphere via blending of two states.
  !!
  SUBROUTINE t_expol_state_w ( expolstate, &  !class
    &                          p_patch,    &  !in
    &                          w_out       )  !inout

    CLASS(t_expol_state), INTENT(INOUT) :: expolstate

    ! In/out variables
    TYPE(t_patch), INTENT(IN)    :: p_patch      ! Domain properties    
    REAL(wp),      INTENT(INOUT) :: w_out(:,:,:) ! Input/output temperature field  

    ! Local variables
    TYPE(t_loop_kit) :: loopkit

    REAL(wp) :: w_cfl
    INTEGER  :: jg, nexpollev
    LOGICAL  :: lpassed
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: sanity_check_msg

    CHARACTER(LEN=*), PARAMETER :: &
      &       routine = modname//':t_expol_state_w'

    !----------------------------------------------

    ! (Note: if 'upatmo_config(jg)%exp%l_expol==.TRUE.' will be checked outside 
    ! the call of this subroutine for efficiency reasons)

    IF (expolstate%ltimer) CALL timer_start(timer_expol)

    ! Domain index
    jg = p_patch%id

    ! Number of layers in extrapolation region
    ! ('nexpollev' applies to full as well as half levels)
    nexpollev = upatmo_config(jg)%exp%nexpollev

    ! For the left-hand side state, we take the vertical wind already 
    ! linearly extrapolated in 'src/atm_dyn_iconam/mo_nh_vert_interp'.
    ! This way the upper-atmosphere extrapolation can be kept relatively simple here
    ! (Note: 'nlev' is set to 'nexpollev' in case of w as well, not to 'nexpollev+1'!)
    loopkit%nlev   = nexpollev
    loopkit%nblks  = p_patch%nblks_c
    loopkit%npromz = p_patch%npromz_c
    CALL copy_state( w_out,             & !in
      &              expolstate%left%w, & !out
      &              loopkit            ) !in

    ! The climatological right-hand side state has already 
    ! been assigned in 't_expol_state_initialize', 
    ! so we can advance to the blending
    CALL blend_states( expolstate%left%w,                & !in (LEFT!)
      &                expolstate%right%w,               & !in (RIGHT!)
      &                w_out,                            & !inout
      &                expolstate%metrics%wfac_blnd_ifc, & !in 
      &                loopkit                           ) !in

    IF ( upatmo_exp_config(jg)%lexpol_sanitycheck .AND. nexpollev > 1 &
      &  .AND. .NOT. expolstate%latbcmode                             ) THEN
      ! Sanity check of the extrapolated vertical wind:
      ! upper bound for the magnitude of the vertical wind 
      ! from the CFL-criterion: 
      !      w_cfl = grid-layer-thickness / time-step
      ! (We assume the nominal grid layer thickness  
      ! to increase monotonically with height.  
      ! So the threshold following from the thickness 
      ! of the lowermost layer within the extrapolation region 
      ! should be the strictest one)
      w_cfl = ABS( (expolstate%metrics%zgpot_mc(nexpollev-1) &
        &   - expolstate%metrics%zgpot_mc(nexpollev) ) /     &
        &     upatmo_config(jg)%dt_dyn_nom                   ) 
      CALL sanity_check ( p_patch     = p_patch,          &  !in
        &                 state       = w_out,            &  !in
        &                 bound       = w_cfl,            &  !in
        &                 keys        = 'cell,upper,abs', &  !in (=> |w| < w_cfl -> ok)
        &                 lpassed     = lpassed,          &  !out
        &                 opt_slev    = 1,                &  !optin
        &                 opt_elev    = nexpollev,        &  !optin
        &                 opt_message = sanity_check_msg  )  !optout
      IF (.NOT. lpassed) THEN
        WRITE(message_text,'(3a)') 'Extrapolated vertical wind violates the CFL-criterion. ', &
          & TRIM(sanity_check_msg), ' A smaller value for expol_blending_scale might help.'
        CALL finish(routine, message_text)
      ENDIF
      
    ENDIF  !IF (upatmo_exp_config(jg)%lexpol_sanitycheck .AND. nexpollev > 1)

    IF (expolstate%ltimer) CALL timer_stop(timer_expol)

  END SUBROUTINE t_expol_state_w

  !-----------------------------------------------------------------------------

  !>
  !! This subroutine extrapolates the normal wind 'vn' 
  !! to the upper atmosphere via blending of two states.
  !!
  SUBROUTINE t_expol_state_vn ( expolstate, &  !class
    &                           p_patch,    &  !in
    &                           vn_out,     &  !inout
    &                           thetav,     &  !in
    &                           exner,      &  !in
    &                           p_metrics,  &  !in
    &                           p_int       )  !in

    CLASS(t_expol_state), INTENT(INOUT) :: expolstate

    ! In/out variables
    TYPE(t_patch),      INTENT(IN)    :: p_patch        ! Domain properties    
    REAL(wp),           INTENT(INOUT) :: vn_out(:,:,:)  ! Input/output horizontal wind 
    REAL(wp),           INTENT(IN)    :: thetav(:,:,:)  ! Virtual potential temperature
    REAL(wp),           INTENT(IN)    :: exner(:,:,:)   ! Exner pressure
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics      ! Model grid properties and background state
    TYPE(t_int_state),  INTENT(IN)    :: p_int          ! Interpolation coefficients

    ! Local variables
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: z_dtheta_c,      &  ! Pot. temp. deviation cell centers
      &                                        z_dtheta_e,      &  ! Pot. temp. deviation edge centers
      &                                        z_theta_e,       &  ! Pot. temp. edge centers
      &                                        z_dexner_c,      &  ! Exner pres. deviation cell centers
      &                                        z_dexner_v,      &  ! Exner pres. deviation vertices
      &                                        z_dexnerdt_e        ! dexner/dt(t->tangential) edge centers

    REAL(wp), ALLOCATABLE, DIMENSION(:)     :: z_decayfac_zgpot    ! Vertical exponential decay factor

    REAL(wp) :: z_latfac(nproma)

    TYPE(t_loop_kit) :: loopkit

    REAL(wp) :: z_lat, z_f, z_decayfac_lat
    REAL(wp) :: z_lat_trop, z_pi_o_lattrop, z_lat_epsilon
    REAL(wp) :: z_hgpot, z_dzgpot
    REAL(wp) :: vn_cfl

    INTEGER :: jg, jb, jk, jc, je, nexpollev
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx 
    INTEGER :: rl_start, rl_end    
    INTEGER :: istat
    LOGICAL :: lpassed
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: sanity_check_msg

    CHARACTER(LEN=*), PARAMETER :: &
      &       routine = modname//':t_expol_state_vn'

    !----------------------------------------------

    ! (Note: if 'upatmo_config(jg)%exp%l_expol==.TRUE.' will be checked outside 
    ! the call of this subroutine for efficiency reasons)

    IF (expolstate%ltimer) CALL timer_start(timer_expol)

    ! (Currently local) tunable parameters:
    !--------------------------------------
    ! Absolute value of tropical latitude 'z_lat_trop', 
    ! for |lat| < 'z_lat_trop' the exner pressure gradient is smoothly 
    ! set to zero, because at the equator geostrophic balance is invalid
    z_lat_trop     = 10._wp                ! In degree
    z_lat_trop     = z_lat_trop * deg2rad  ! In radian
    z_pi_o_lattrop = pi / z_lat_trop
    ! Threshold latitude: for |lat| < 'z_lat_epsilon' 
    ! (where 'z_lat_epsilon < z_lat_trop' should hold!).
    ! The geostrophic wind is set to zero, to avoid division 
    ! by a Coriolis-parameter being zero in its computation
    ! (A constant value is used for convenienc, 
    ! but one might link it to the horizontal mesh size as well)
    z_lat_epsilon = 0.5_wp                   ! In degree
    z_lat_epsilon = z_lat_epsilon * deg2rad  ! In radian

    ! Domain index
    jg = p_patch%id 

    nexpollev = upatmo_config(jg)%exp%nexpollev

    ! (Geopotential) scale height for (exponential) decay of vn
    z_hgpot = upatmo_exp_config(jg)%expol_vn_decay_scale

    ! For the left-hand side state, we take the wind field already 
    ! linearly extrapolated in 'src/atm_dyn_iconam/mo_nh_vert_interp', 
    ! this way the upper-atmosphere extrapolation can be kept more simple here
    loopkit%nlev   = nexpollev
    loopkit%nblks  = p_patch%nblks_e
    loopkit%npromz = p_patch%npromz_e
    CALL copy_state( vn_out,             & !in
      &              expolstate%left%vn, & !out
      &              loopkit             ) !in    

    ! Unfortunately, the computation of the 'climatological' 
    ! right-hand side state is much more complicated than in the previous 
    ! subroutines, because we have to compute the geostrophic horizontal wind 
    ! according to the formula: f * vn = cp * thetav * dexner/dt.  (I)
    ! This requires to interpolate 'thetav' to edge centers 
    ! and 'exner' to vertices. 
    ! Note that the geostrophic balance in the deep-atmosphere includes 
    ! additional metric terms, which were neglected in (I), to simplify the computation.
    ! (We extrapolate to the geostrophic wind, in order to potentially reduce 
    ! the dynamical tendencies of the horizontal component of the momentum budget  
    ! during the first time steps of the integration, and this in turn to potentially 
    ! alleviate the adjustment process in the upper atmosphere during model spin-up)

    ! Note: if the model, from which the initial state is taken, uses a sponge layer 
    ! that acts directly on the horizontal wind, and the sponge layer region 
    ! is part of its output, it might be worth considering, 
    ! to use for 'upatmo_nml: expol_start_height' a value, which is comparable 
    ! to the start height of the sponge layer. (Otherwise, more or less strong 
    ! adjustment processes may take place below 'expol_start_height', 
    ! since a sponge layer forcing term is not part of the horizontal wind equation used in ICON.)

    ! 1st Interpolation of the virtual potential temperature to edge centers
    !-----------------------------------------------------------------------

    ! (Some of the fields allocated in the following could actually be allocated later on, 
    ! but they are allocated all at once to have less open-mp-thread interruptions.  
    ! However, it might well be that memory economy would have been more important than that)
    ALLOCATE( z_dtheta_c   (nproma, nexpollev, p_patch%nblks_c), & 
      &       z_dtheta_e   (nproma, nexpollev, p_patch%nblks_e), & 
      &       z_theta_e    (nproma, nexpollev, p_patch%nblks_e), &
      &       z_dexner_c   (nproma, nexpollev, p_patch%nblks_c), & 
      &       z_dexner_v   (nproma, nexpollev, p_patch%nblks_v), & 
      &       z_dexnerdt_e (nproma, nexpollev, p_patch%nblks_e), &
      &       STAT=istat                                          )
    IF (istat /= SUCCESS) CALL finish(routine, 'Allocation failed')

!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)

    ! (Initialization with zeroes might perhaps be more secure)
    CALL init( z_dtheta_c(:,:,:), lacc=.FALSE. )
    CALL init( z_dtheta_e(:,:,:), lacc=.FALSE. )
    CALL init( z_theta_e(:,:,:), lacc=.FALSE. )
    CALL init( z_dexner_c(:,:,:), lacc=.FALSE. )
    CALL init( z_dexner_v(:,:,:), lacc=.FALSE. )
    CALL init( z_dexnerdt_e(:,:,:), lacc=.FALSE. )

    rl_start   = 1
    rl_end     = min_rlcell
    i_startblk = p_patch%cells%start_block(rl_start) 
    i_endblk   = p_patch%cells%end_block(rl_end)  
    
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
      DO jk = 1, nexpollev
        DO jc = i_startidx, i_endidx
          ! Only the deviation from the reference state is interpolated
          z_dtheta_c(jc,jk,jb) = thetav(jc,jk,jb) - p_metrics%theta_ref_mc(jc,jk,jb)
        ENDDO  !jc
      ENDDO  !jk
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL
    
    ! Interpolate virtual potential temperature deviations from cell centers to edge centers  
    CALL cells2edges_scalar( z_dtheta_c, p_patch, p_int%c_lin_e, z_dtheta_e,       &
      &                      opt_slev=1, opt_elev=nexpollev, opt_fill_latbc=.TRUE. )

!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)
    
    rl_start   = 1
    rl_end     = min_rledge
    i_startblk = p_patch%edges%start_block(rl_start) 
    i_endblk   = p_patch%edges%end_block(rl_end)  
    
!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end) 
      DO jk = 1, nexpollev
        DO je = i_startidx, i_endidx 
          ! Add interpolated deviation to reference on edges
          z_theta_e(je,jk,jb) = p_metrics%theta_ref_me(je,jk,jb) + z_dtheta_e(je,jk,jb)
        ENDDO  !je
      ENDDO  !jk
    ENDDO  !jb
!$OMP END DO

    ! 2nd Interpolation of the exner pressure to vertices
    !----------------------------------------------------
     
    rl_start   = 1
    rl_end     = min_rlcell
    i_startblk = p_patch%cells%start_block(rl_start) 
    i_endblk   = p_patch%cells%end_block(rl_end)  
    
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
      DO jk = 1, nexpollev
        DO jc = i_startidx, i_endidx
          ! Compute deviation of Exner pressure from reference state  
          z_dexner_c(jc,jk,jb) = exner(jc,jk,jb) - p_metrics%exner_ref_mc(jc,jk,jb)
        ENDDO  !jc
      ENDDO  !jk
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

    ! Interpolate deviation Exner pressure from cell centers to vertices 
    CALL cells2verts_scalar( z_dexner_c, p_patch, p_int%cells_aw_verts, z_dexner_v, &
      &                      opt_slev=1, opt_elev=nexpollev                         )

    ! Compute pressure gradient along edge        
    ! (Note: 'grad_fd_tang' is not corrected for spherical geometry, 
    ! but rather the output field is multiplied by a modification factor later on)
    CALL grad_fd_tang(z_dexner_v, p_patch, z_dexnerdt_e, opt_slev=1, opt_elev=nexpollev) 

    ! 3rd Compute geostrophic wind
    !-----------------------------

!$OMP PARALLEL PRIVATE (rl_start,rl_end,i_startblk,i_endblk)
    
    rl_start   = 1
    rl_end     = min_rledge
    i_startblk = p_patch%edges%start_block(rl_start) 
    i_endblk   = p_patch%edges%end_block(rl_end)  

!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx, &
!$OMP            z_lat, z_f, z_decayfac_lat, z_latfac) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end) 
      DO je = i_startidx, i_endidx 
        ! (Pre-computation of latitudinal factor for geostorphic wind: 
        ! vn = thetav * dexner/dt * cp / f (* decayfac(lat)),
        !                           ------------------------
        ! because it requires if-queries unfortunately)
        ! Latitude of edge
        z_lat = p_patch%edges%center(je,jb)%lat 
        ! Coriolis parameter
        z_f = p_patch%edges%f_e(je,jb)   
        IF (ABS(z_lat) > ABS(z_lat_trop)) THEN  ! (Likely the most frequent case)
          ! Standard case without any special modification ('z_decayfac_lat=1')
          z_latfac(je) = cpd / z_f
        ELSEIF(ABS(z_lat) > ABS(z_lat_epsilon)) THEN
          ! Exner pressure gradient is smoothly set to zero towards the equator, 
          ! where the geostrophic balance is invalid
          z_decayfac_lat = 0.5_wp * ( 1_wp - COS( ABS(z_lat) * z_pi_o_lattrop ) )
          z_latfac(je)   = z_decayfac_lat * cpd / z_f
        ELSE
          ! Geostrophic wind shall be zero at about the equator (avoid division by 'z_f=0')
          z_latfac(je) = 0._wp
        ENDIF
      ENDDO  !je

      DO jk = 1, nexpollev
        DO je = i_startidx, i_endidx 
          ! Geostrophic wind 
          ! (with deep-atmosphere modification of horizontal gradient 'z_dexnerdt_e')
          expolstate%right%vn(je,jk,jb) = z_theta_e(je,jk,jb) * z_dexnerdt_e(je,jk,jb)  & 
            &                           * z_latfac(je) * p_metrics%deepatmo_gradh_mc(jk)
        ENDDO  !je
      ENDDO  !jk
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

    ! Now, we can blend left-hand side and right-hand side states
    CALL blend_states( expolstate%left%vn,              & !in (LEFT!)
      &                expolstate%right%vn,             & !in (RIGHT!)
      &                vn_out,                          & !inout
      &                expolstate%metrics%wfac_blnd_mc, & !in 
      &                loopkit                          ) !in    

    ! 4th Some postprocessing
    !------------------------
   
    ! The extrapolated horizontal wind might have too large magnitudes to maintain stability, 
    ! therefore it is multiplied by an exponentially decaying factor of the form:
    ! 'decayfac = (1+dzgpot/hgpot)*exp(-dzgpot/hgpot)'
    ! (This violates to some extent the geostrophic balance, but stability 
    ! shall take precedence over that)
    ALLOCATE( z_decayfac_zgpot(nexpollev), STAT=istat )
    IF (istat /= SUCCESS) CALL finish(routine, 'Allocation failed')
    DO jk = 1, nexpollev 
      ! Geopotential height difference 
      ! (Note: due to the way 'nexpollev' has been determined, 'zgpot_start - zgpot_mc(nexpollev)' 
      ! is not necessarily positive, but for the particular form of the decay factor, 
      ! a positive difference would be preferable, so the height difference is not computed 
      ! relative to 'zgpot_start', but relative to 'zgpot_mc(nexpollev)')
      z_dzgpot = expolstate%metrics%zgpot_mc(jk) - expolstate%metrics%zgpot_mc(nexpollev)
      z_decayfac_zgpot(jk) = ( 1._wp + z_dzgpot / z_hgpot ) * EXP( -z_dzgpot / z_hgpot )
    ENDDO  !jk

    ! Multiply the wind with the decay factor
    ! (For efficiency reasons, this could actually be inlined in some of the previous loops 
    ! or in 'blend_states', but for reasons of clarity and generalizability, we would prefer modularity here)
    CALL multiply_state( vn_out,           &  !inout
      &                  z_decayfac_zgpot, &  !in 
      &                  loopkit           )  !in

    ! Clean-up
    DEALLOCATE( z_dtheta_c,       & 
      &         z_dtheta_e,       & 
      &         z_theta_e,        &
      &         z_dexner_c,       & 
      &         z_dexner_v,       & 
      &         z_dexnerdt_e,     &
      &         z_decayfac_zgpot, &
      &         STAT=istat        )
    IF (istat /= SUCCESS) CALL finish(routine, 'Deallocation failed')
    
    IF (upatmo_exp_config(jg)%lexpol_sanitycheck .AND. .NOT. expolstate%latbcmode) THEN
      
      ! Sanity check of the extrapolated horizontal wind:      
      ! upper bound for the magnitude of the horizontal wind 
      ! from the CFL-criterion: 
      !     vn_cfl = horizontal-mesh-size / time-step
      vn_cfl = ABS( p_patch%geometry_info%mean_characteristic_length / &
        &      upatmo_config(jg)%dt_dyn_nom                            )
      CALL sanity_check ( p_patch     = p_patch,          &  !in
        &                 state       = vn_out,           &  !in
        &                 bound       = vn_cfl,           &  !in
        &                 keys        = 'edge,upper,abs', &  !in (=> |vn| < vn_cfl -> ok)
        &                 lpassed     = lpassed,          &  !out
        &                 opt_slev    = 1,                &  !optin
        &                 opt_elev    = nexpollev,        &  !optin
        &                 opt_message = sanity_check_msg  )  !optout
      IF (.NOT. lpassed) THEN
        WRITE(message_text,'(3a)') 'Extrapolated horizontal wind violates the CFL-criterion. ', &
          & TRIM(sanity_check_msg), ' A smaller value for expol_blending_scale might help.'
        CALL finish(routine, message_text)
      ENDIF

    ENDIF  !IF (upatmo_exp_config(jg)%lexpol_sanitycheck)
    
    IF (expolstate%ltimer) CALL timer_stop(timer_expol)
    
  END SUBROUTINE t_expol_state_vn

  !-----------------------------------------------------------------------------

  !>
  !! This subroutine finalizes the upper-atmosphere extrapolation.
  !!
  SUBROUTINE t_expol_state_finalize ( expolstate, &  !class
    &                                 p_patch     )  !in

    CLASS(t_expol_state), INTENT(INOUT) :: expolstate

    ! In/out variables
    TYPE(t_patch), INTENT(IN) :: p_patch ! Domain properties    

    ! Local variables
    INTEGER :: jg
    CHARACTER(LEN=*), PARAMETER :: &
      &       routine = modname//':t_expol_state_finalize'

    !----------------------------------------------

    ! (Note: if 'upatmo_config(jg)%exp%l_expol==.TRUE.' will be checked outside 
    ! the call of this subroutine for efficiency reasons)

    IF (expolstate%ltimer) CALL timer_start(timer_expol)

    ! Domain index
    jg = p_patch%id

    ! Destruct the variables of the extrapolation state
    IF (expolstate%linitialized) THEN 
      CALL expolstate%clim%finalize() 
      CALL expolstate%metrics%finalize() 
      CALL expolstate%left%finalize() 
      CALL expolstate%right%finalize() 
      expolstate%linitialized = .FALSE.
    ENDIF

    ! Write some info to screen                  
    IF (msg_level >= imsg_thr%high) THEN
      WRITE(message_text,'(a,i4)') 'Finished extrapolation to upper atmosphere in domain: ', jg
      CALL message(routine, message_text)
    ENDIF

    IF (expolstate%ltimer) CALL timer_stop(timer_expol)
    expolstate%ltimer = .FALSE.

  END SUBROUTINE t_expol_state_finalize

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

  !---------------------------------------------------------------
  !           Auxiliary routines for t_expol_state_...
  !---------------------------------------------------------------

  !>
  !! This subroutine initializes the extrapolation metrics 
  !! associated with the geopotential height.
  !!
  SUBROUTINE init_expol_metrics ( expolmetrics, &  !inout
    &                             p_patch       )  !in

    ! In/out variables
    TYPE(t_expol_metrics), INTENT(INOUT) :: expolmetrics  ! Extrapolation metrics
    TYPE(t_patch),         INTENT(IN)    :: p_patch       ! Domain properties    

    ! Local variables
    REAL(wp) :: z_z_mc, z_z_ifc, z_zgpot_mc, z_zgpot_ifc
    REAL(wp) :: z_zgpot_start, z_hgpot
    REAL(wp) :: z_trafo, z_pi_o_hgpot

    INTEGER :: jg, jk, nshift, nexpollev

    CHARACTER(LEN=*), PARAMETER :: &
      &       routine = modname//':init_expol_metrics'

    !----------------------------------------------

    ! Domain index
    jg = p_patch%id

    ! Layer index shift for vertical nesting
    nshift = p_patch%nshift_total

    ! Number of layers in extrapolation region
    nexpollev = upatmo_config(jg)%exp%nexpollev

    ! Transformation factor
    IF (.NOT. ldeepatmo) THEN 
      ! Constant gravitational acceleration: 
      ! geopotential height coincides with geometric height
      z_trafo = 0._wp
    ELSE
      ! Gravitational acceleration varies with height:
      ! geopotential height differs from geometric height
      z_trafo = 1._wp / grid_sphere_radius
    ENDIF

    ! Geographic height above which extrapolation takes place
    z_zgpot_start = upatmo_exp_config(jg)%expol_start_height
    ! Transform to geopotential height and store 
    z_zgpot_start = z_zgpot_start / ( 1._wp + z_trafo * z_zgpot_start )
    expolmetrics%zgpot_start = z_zgpot_start
    ! Geopotential blending scale height
    ! ('expol_blending_scale' is assumed to be a geopotential height difference)
    z_hgpot      = upatmo_exp_config(jg)%expol_blending_scale
    z_pi_o_hgpot = pi / z_hgpot

    ! Fortunately, the extrapolation to the upper atmosphere 
    ! takes place above 'flat_height', so we do not have to care 
    ! about horizontal dependencies here
    DO jk=1, nexpollev
      ! Geometric height of interface
      z_z_ifc = vct_a(jk+nshift)
      ! Geometric height of cell centers and edges
      z_z_mc  = 0.5_wp * ( vct_a(jk+nshift) + vct_a(jk+nshift+1) )
      ! Compute geopotential heights
      z_zgpot_ifc = z_z_ifc / ( 1._wp + z_trafo * z_z_ifc )
      z_zgpot_mc  = z_z_mc / ( 1._wp + z_trafo * z_z_mc )
      expolmetrics%zgpot_ifc(jk) = z_zgpot_ifc
      expolmetrics%zgpot_mc(jk)  = z_zgpot_mc
      ! Compute the two blending factors for the field linearly extrapolated 
      ! upwards from the IFS-data (->left) and the field derived from a climatology 
      ! covering the upper atmosphere (->right)
      !
      !                         { 1,                                       for zgpot < zgpot_start
      ! wfac_blnd_left(zgpot) = { [1+cos((zgpot-zgpot_start)/hgpot*pi)]/2, for zgpot_start <= zgpot <= zgpot_start + hgpot
      !                         { 0,                                       for zgpot_start + hgpot < zgpot
      ! wfac_blnd_right(zgpot) = 1 - wfac_blnd_left(zgpot)
      !
      ! (Note: the blending factors are used for all variables equally.  
      ! In addition, we do not discriminate between extensive variables (e.g. mass, heat, momentum etc.), 
      ! for which "physical additivity" holds, and intensive variables 
      ! (e.g. temperature, velocity etc.), for which "physical additivity" does not hold)
      ! 1st For ifc
      IF (z_zgpot_ifc < z_zgpot_start) THEN
        expolmetrics%wfac_blnd_ifc(iblnd%left, jk) = 1._wp
      ELSEIF (z_zgpot_ifc > z_zgpot_start+z_hgpot) THEN
        expolmetrics%wfac_blnd_ifc(iblnd%left, jk) = 0._wp
      ELSE
        expolmetrics%wfac_blnd_ifc(iblnd%left, jk) =  & 
          & 0.5_wp * ( 1._wp + COS( ( z_zgpot_ifc - z_zgpot_start ) * z_pi_o_hgpot ) )
      ENDIF
      expolmetrics%wfac_blnd_ifc(iblnd%right, jk) =  & 
        & 1._wp - expolmetrics%wfac_blnd_ifc(iblnd%left, jk)
      ! 2nd For mc/me
      IF (z_zgpot_mc < z_zgpot_start) THEN
        expolmetrics%wfac_blnd_mc(iblnd%left, jk) = 1._wp
      ELSEIF (z_zgpot_mc > z_zgpot_start+z_hgpot) THEN
        expolmetrics%wfac_blnd_mc(iblnd%left, jk) = 0._wp
      ELSE
        expolmetrics%wfac_blnd_mc(iblnd%left, jk) =  & 
          & 0.5_wp * ( 1._wp + COS( ( z_zgpot_mc - z_zgpot_start ) * z_pi_o_hgpot ) )
      ENDIF
      expolmetrics%wfac_blnd_mc(iblnd%right, jk) =  & 
        & 1._wp - expolmetrics%wfac_blnd_mc(iblnd%left, jk)
    ENDDO  !jk

  END SUBROUTINE init_expol_metrics

  !-----------------------------------------------------------------------------

  !>
  !! This subroutine initializes the vertical climatological 
  !! temperature profile on grid levels.
  !!
  SUBROUTINE init_clim_state ( expolclim,    &  !inout
    &                          expolmetrics, &  !in
    &                          p_patch       )  !in

    ! In/out variables
    TYPE(t_clim_state),    INTENT(INOUT) :: expolclim     ! Basic extrapolation climatologie
    TYPE(t_expol_metrics), INTENT(IN)    :: expolmetrics  ! Extrapolation metrics
    TYPE(t_patch),         INTENT(IN)    :: p_patch       ! Domain properties    

    ! Local variables
    REAL(wp), ALLOCATABLE, DIMENSION(:)   :: z_mtr_A1,     &
      &                                      z_mtr_B1,     &
      &                                      z_mtr_C1,     &
      &                                      z_mtr_D1,     & 
      &                                      z_mtr_a2,     & 
      &                                      z_mtr_b2,     & 
      &                                      z_dinvtemp_dz
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: z_spln

    REAL(wp) :: z_zgpot, z_dzgpot, z_dzgpotm1
    REAL(wp) :: z_zgpot_transit, z_zgpot_rel
    REAL(wp) :: z_temp, z_tempm1, z_tempp1
    REAL(wp) :: z_invtemp, z_invtempp1
    REAL(wp) :: z_temp_min, z_temp_max
    REAL(wp) :: z_temp_infty, z_temp_transit
    REAL(wp) :: z_dtempdzgpot_transit
    REAL(wp) :: z_q0, z_q1, z_invhb
    REAL(wp) :: z_temp_u, z_temp_l
    REAL(wp) :: z_zgpot_u, z_zgpot_l
    REAL(wp) :: z_dtempdzgpot, z_dtempdzgpot_min

    INTEGER :: jg, jk, jclim, nexpollev, kclim, istart
    INTEGER :: istat
    INTEGER, PARAMETER :: kdefault = -1,  & 
      &                   kbates   =  0

    REAL(wp), PARAMETER :: grav_o_cpd = grav / cpd
    REAL(wp), PARAMETER :: z_eps      = 1.0e-10_wp

    CHARACTER(LEN=*), PARAMETER :: &
      &       routine = modname//':init_clim_state'

    !----------------------------------------------

    ! In order to get a steady temperature profile 
    ! from the Fleming climatology, we make use of a spline interpolation 
    ! (actually, not temp, but 1/temp is interpolated, but only for 
    ! a "historical" reason). The spline interpolation uses a polynomial 
    ! of third order within the interval between two sampling points. 
    ! The four spline coefficients of this polynomial are fixed by 
    ! the conditions that at a sampling point the temperature values 
    ! and vertical temperature derivatives obtained from the polynomial 
    ! in the intervals above and below the sampling point are identical. 
    ! This leads to a matrix equation which has to be inverted .
    ! (For this we use the same method as is used to solve  
    ! the vertically implicit part of the dynamics).
    ! (Note: this is a purely technical step without any particular 
    ! meaning for the other subroutines in this module, so we avoid going 
    ! into much detail in the comments.)

    ! Domain index
    jg = p_patch%id

    ! Number of layers in extrapolation region
    nexpollev = upatmo_config(jg)%exp%nexpollev

    ! Allocate auxiliary fields for matrix inversion
    ! Check if 'nclim' is large enough
    IF (nclim < 3) CALL finish(routine, 'nclim >= 3 required')
    ALLOCATE( z_mtr_A1(2:nclim-1),    & 
      &       z_mtr_B1(2:nclim-1),    &
      &       z_mtr_C1(2:nclim-1),    &
      &       z_mtr_D1(2:nclim-1),    &
      &       z_mtr_a2(2:nclim-1),    &
      &       z_mtr_b2(2:nclim-1),    &
      &       z_dinvtemp_dz(1:nclim), &
      &       z_spln(1:4, 1:nclim),   &
      &       STAT=istat              )
    IF (istat /= SUCCESS) CALL finish(routine, 'Allocation failed')

    ! Initialize matrix elements
    DO jclim=2, nclim-1
      ! Geopotential height difference between sampling points of climatology
      z_dzgpot   = clim_rawdata%zgpot(jclim) - clim_rawdata%zgpot(jclim+1)
      z_dzgpotm1 = clim_rawdata%zgpot(jclim-1) - clim_rawdata%zgpot(jclim)
      z_temp     = clim_rawdata%temp(jclim)
      z_tempm1   = clim_rawdata%temp(jclim-1)
      z_tempp1   = clim_rawdata%temp(jclim+1)
      z_mtr_A1(jclim) = 0.5_wp / z_dzgpot
      z_mtr_C1(jclim) = 0.5_wp / z_dzgpotm1
      z_mtr_B1(jclim) = 2.0_wp * ( z_mtr_A1(jclim) + z_mtr_C1(jclim) )
      z_mtr_D1(jclim) = (3.0_wp/2.0_wp) * ( -1.0_wp / ( z_tempp1 * z_dzgpot**2 ) +   & 
        &               ( 1.0_wp / z_dzgpot**2 - 1.0_wp / z_dzgpotm1**2 ) / z_temp + & 
        &               1.0_wp / ( z_tempm1 * z_dzgpotm1**2 ) )
    ENDDO  !jclim

    ! To solve the matrix equation, we have to assign some boundary values:
    ! vertical derivative of inverse temperature, 
    ! were we make use of d(1/T)/dzgopt = -(1/T^2)*dT/dzgopt
    z_dinvtemp_dz(1)     = -( clim_rawdata%temp(1) - clim_rawdata%temp(2) ) /             &
      &                     ( clim_rawdata%zgpot(1) - clim_rawdata%zgpot(2) ) /           & 
      &                     clim_rawdata%temp(1)**2
    z_dinvtemp_dz(nclim) = -( clim_rawdata%temp(nclim-1) - clim_rawdata%temp(nclim) ) /   &
      &                     ( clim_rawdata%zgpot(nclim-1) - clim_rawdata%zgpot(nclim) ) / & 
      &                     clim_rawdata%temp(nclim)**2
    ! Some matrix elements
    z_mtr_b2(2) = ( z_mtr_D1(2) - z_mtr_C1(2) * z_dinvtemp_dz(1) ) / z_mtr_B1(2)
    z_mtr_a2(2) = z_mtr_A1(2) / z_mtr_B1(2)

    ! (Direct) inversion of matrix
    DO jclim=3, nclim-1
      z_mtr_b2(jclim) = ( z_mtr_D1(jclim) - z_mtr_C1(jclim) * z_mtr_b2(jclim-1) ) / &
        &               ( z_mtr_B1(jclim) - z_mtr_C1(jclim) * z_mtr_a2(jclim-1) )
      z_mtr_a2(jclim) = z_mtr_A1(jclim) / ( z_mtr_B1(jclim) - z_mtr_C1(jclim) * z_mtr_a2(jclim-1) )
    ENDDO

    DO jclim=nclim-1, 2, -1
      z_dinvtemp_dz(jclim) = z_mtr_b2(jclim) - z_mtr_a2(jclim) * z_dinvtemp_dz(jclim+1)
    ENDDO

    ! Auxiliary fields for matrix are no longer required
    DEALLOCATE( z_mtr_A1,  & 
      &         z_mtr_B1,  &
      &         z_mtr_C1,  &
      &         z_mtr_D1,  &
      &         z_mtr_a2,  &
      &         z_mtr_b2,  &
      &         STAT=istat )
    IF (istat /= SUCCESS) CALL finish(routine, 'Deallocation failed')

    ! Compute the four spline coefficients for the polynomial representing 1/temp 
    ! for the nclim-1 intervals between the sampling points of the climatology
    DO jclim=1, nclim-1
      z_dzgpot      = clim_rawdata%zgpot(jclim) - clim_rawdata%zgpot(jclim+1)
      z_invtemp     = 1.0_wp / clim_rawdata%temp(jclim)
      z_invtempp1   = 1.0_wp / clim_rawdata%temp(jclim+1)
      z_spln(1, jclim) = z_invtempp1
      z_spln(2, jclim) = z_dinvtemp_dz(jclim+1) * z_dzgpot
      z_spln(3, jclim) = 3.0_wp * ( z_invtemp - z_invtempp1 ) -                          &
        &                2.0_wp * z_dinvtemp_dz(jclim+1) * z_dzgpot -                    & 
        &                z_dinvtemp_dz(jclim) * z_dzgpot
      z_spln(4, jclim) = ( z_dinvtemp_dz(jclim+1) + z_dinvtemp_dz(jclim) ) * z_dzgpot -  & 
        &                2.0_wp * ( z_invtemp - z_invtempp1 )
    ENDDO

    ! Unfortunately, the Fleming climatology covers only the air column up to 
    ! a geopotential height of Zgpot = 120 km. 
    ! Above the exponential reference temperature profile of: 
    ! Bates, D.R. (1959), "Some problems concerning the terrestrial atmosphere above the 100 km level", 
    ! Proc. R. Soc. London, Ser. A., 253, 451 - 462, 
    ! compare also:
    ! Hedin, A.E. (1983), "A revised thermospheric model based on mass spectrometer 
    ! and incoherent scatter data: MSIS-83", J. Geophys. Res., 88, 10,170 - 10,188
    ! is used. It reads:
    ! T_B(zgpot) = q0 + q1 * exp[ -(zgpot - Zgpot) / H_B ], 
    ! where q0, q1 and H_B are constants. 
    ! The constants q0, q1 and H_B are determined in such a way, that: 
    ! T_B(z->infinity <=> zgpot->Earth radius) =: T_infty = 'expol_temp_infty' 
    ! (e.g. 'expol_temp_infty' = 1035 K, the mean temperature of the exosphere), 
    ! T_B(zgpot=Zgpot) =: T_Zgpot = temperature from Fleming climatology at Zgpot, 
    ! dT_B/dzgpot|_(zgpot=Zgpot) =: dTdZgpot = vertical temperature derivative of 
    ! spline interpolation of Fleming climatology at Zgpot
    ! We find: 
    ! q0 = T_infty, 
    ! q1 = T_Zgpot - T_infty, 
    ! H_B = (T_infty - T_Zgpot) / dTdZgpot. 

    ! Temperature and geopotential height at top of Fleming climatology
    z_temp_transit  = clim_rawdata%temp(1)
    z_zgpot_transit = clim_rawdata%zgpot(1)

    ! Temperature for z -> infinity (zgpot -> radius of Earth)
    z_temp_infty = upatmo_exp_config(jg)%expol_temp_infty

    ! Vertical temperature gradient from spline interpolation 
    ! at top of Fleming climatology 
    ! (Note: the splines interpolate 1/temp)
    z_dzgpot   = clim_rawdata%zgpot(1) - clim_rawdata%zgpot(2)
    z_dtempdzgpot_transit = -( z_spln(2,1) + 2.0_wp * z_spln(3,1) + 3.0_wp * z_spln(4,1) ) /           &
      &                      ( z_dzgpot * ( z_spln(1,1) + z_spln(2,1) + z_spln(3,1) + z_spln(4,1) )**2 )

    ! Coefficients for Bates climatology
    z_q0    = z_temp_infty 
    z_q1    = z_temp_transit - z_temp_infty
    z_invhb = -z_dtempdzgpot_transit / z_q1  ! = 1/H_B

    ! Scale height should be positive
    ! (Because we do not use a "physically motivated" value for H_B, 
    ! but "misuse" this parameter to get a steady transition between the two climatologies, 
    ! it is not guaranteed that H_B fits the requirements of the Bates-formula)
    IF (z_invhb <= 0._wp) CALL finish(routine, 'Scale height of Bates climatology is negative')

    ! Now, we can compute the climatological temperature values on grid levels:
    istart            = 1 
    z_temp_u          = z_temp_infty
    z_temp_l          = z_temp_infty
    z_zgpot_u         = 2._wp * expolmetrics%zgpot_mc(1)
    z_zgpot_l         = expolmetrics%zgpot_mc(1)
    z_dtempdzgpot_min = 0._wp  ! > -grav_o_cpd
    DO jk=1, nexpollev
      ! Find position in climatologies
      kclim = kdefault
      z_zgpot = expolmetrics%zgpot_mc(jk)
      IF (z_zgpot > clim_rawdata%zgpot(1)) THEN
        ! In Bates climatology
        kclim = kbates
      ELSE
        ! In Fleming climatology: 
        ! in which interval lies zgpot?
        DO jclim=istart, nclim-1
          IF (z_zgpot >  clim_rawdata%zgpot(jclim+1) .AND.  &
            & z_zgpot <= clim_rawdata%zgpot(jclim  )) THEN
            kclim  = jclim 
            istart = jclim  ! For the next cycle of jk
            EXIT
          ENDIF
        ENDDO  !jclim
      ENDIF
      ! Compute temperature
      SELECT CASE (kclim)
      CASE(kbates)
        ! Bates climatology
        z_temp = z_q0 + z_q1 * EXP( -( z_zgpot - z_zgpot_transit ) * z_invhb )
      CASE(kbates+1:nclim-1)
        ! Fleming climatology:
        ! normalized relative position in interval
        z_zgpot_rel = ( z_zgpot - clim_rawdata%zgpot(kclim+1) ) /               &
          &           ( clim_rawdata%zgpot(kclim) - clim_rawdata%zgpot(kclim+1) )
        ! Inverse temperature from splines
        z_invtemp = z_spln(1,kclim) + z_spln(2,kclim) * z_zgpot_rel +  & 
          &         z_spln(3,kclim) * z_zgpot_rel**2 + z_spln(4,kclim) * z_zgpot_rel**3
        z_temp    = 1._wp / z_invtemp
      CASE DEFAULT 
        CALL finish(routine, "Grid level not in climatologies")
      END SELECT
      ! Hand over
      expolclim%temp(jk) = z_temp
      ! Compute temperature gradient for sanity check
      z_temp_u      = z_temp_l 
      z_temp_l      = z_temp
      z_zgpot_u     = z_zgpot_l
      z_zgpot_l     = z_zgpot
      z_dtempdzgpot = ( z_temp_u - z_temp_l ) / MAX(z_eps, ABS( z_zgpot_u - z_zgpot_l ))
      IF (jk == 2) THEN
        z_dtempdzgpot_min = z_dtempdzgpot
      ELSEIF(jk > 2) THEN
        IF (z_dtempdzgpot < z_dtempdzgpot_min) z_dtempdzgpot_min = z_dtempdzgpot
      ENDIF
    ENDDO  !jk

    ! Some rudimentary sanity check
    z_temp_min = MINVAL(expolclim%temp(:))
    z_temp_max = MAXVAL(expolclim%temp(:))
    IF (z_temp_min <= 0._wp) CALL finish(routine, 'Temperature lower than 0')
    IF (z_temp_max >= z_temp_infty) CALL finish(routine, 'Temperature higher than T_infty')
    IF (z_dtempdzgpot_min < -grav_o_cpd) THEN
      WRITE(message_text,'(a,F14.10,a,F14.10,a)') 'Min. temperature gradient: ', z_dtempdzgpot_min, &
        & 'K/m < negative dry adiabatic lapse rate: ', -grav_o_cpd, 'K/m'
      CALL finish(routine, message_text)
    ENDIF

    ! Print some info                
    IF (msg_level >= imsg_thr%high) THEN
      WRITE(message_text,'(a)') 'Temperature profile on grid levels from Fleming & Bates climatologies:'
      CALL message(routine, message_text)
      WRITE(message_text,'(a,F14.4)') 'Min. temperature: T_min (K) = ', z_temp_min
      CALL message(routine, message_text)
      WRITE(message_text,'(a,F14.4)') 'Max. temperature: T_max (K) = ', z_temp_max
      CALL message(routine, message_text)
      WRITE(message_text,'(a,F14.10)') 'Min. temperature gradient: dT/dzgpot_min (K/m) = ', z_dtempdzgpot_min
      CALL message(routine, message_text)
      WRITE(message_text,'(a,F14.10)') 'Negative dry adiabatic lapse rate: -grav/cpd (K/m) = ', -grav_o_cpd
      CALL message(routine, message_text)
    ENDIF  !IF (msg_level >= imsg_thr%high)

    ! Deallocate what has not yet been deallocated
    DEALLOCATE( z_dinvtemp_dz, & 
      &         z_spln,        &
      &         STAT=istat     )
    IF (istat /= SUCCESS) CALL finish(routine, 'Deallocation failed')

  END SUBROUTINE init_clim_state

  !-----------------------------------------------------------------------------

  !>
  !! This subroutine copies states.
  !!
  SUBROUTINE copy_state ( state_orgnl, &  !in
    &                     state_copy,  &  !out
    &                     loopkit      )  !in

    ! In/out variables
    REAL(wp),         INTENT(IN)    :: state_orgnl(:,:,:) ! Original
    REAL(wp),         INTENT(INOUT) :: state_copy(:,:,:)  ! Copy
    TYPE(t_loop_kit), INTENT(IN)    :: loopkit            ! Loop variables

    ! Local variables
    INTEGER :: jb, jc, jk, nlen  ! (jc is habitual placeholder for jc, je, jv)

    CHARACTER(LEN=*), PARAMETER :: &
      &       routine = modname//':copy_state'

    !----------------------------------------------

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO jb=1, loopkit%nblks
      IF (jb /= loopkit%nblks) THEN
        nlen = nproma
      ELSE
        nlen = loopkit%npromz
        ! (Initialization of 'state_copy' should have happend 
        ! elsewhere, see e.g. call of 'init'-routines below)
      ENDIF
      DO jk=1, loopkit%nlev
        DO jc=1, nlen
          state_copy(jc,jk,jb) = state_orgnl(jc,jk,jb)
         ENDDO  !jc
      ENDDO  !jk
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE copy_state

  !-----------------------------------------------------------------------------

  !>
  !! This subroutine blends to states.
  !!
  SUBROUTINE blend_states ( state_left,  &  !in
    &                       state_right, &  !in
    &                       state_blnd,  &  !inout
    &                       wfac_blnd,   &  !in
    &                       loopkit      )  !in

    ! In/out variables
    REAL(wp),         INTENT(IN)    :: state_left(:,:,:)   ! Left-hand side state
    REAL(wp),         INTENT(IN)    :: state_right(:,:,:)  ! Right-hand side state
    REAL(wp),         INTENT(INOUT) :: state_blnd(:,:,:)   ! Blended state
    REAL(wp),         INTENT(IN)    :: wfac_blnd(:,:)      ! Blending factors
    TYPE(t_loop_kit), INTENT(IN)    :: loopkit             ! Loop variables

    ! Local variables
    INTEGER :: jb, jc, jk, nlen  ! (jc is habitual placeholder for jc, je, jv)

    CHARACTER(LEN=*), PARAMETER :: &
      &       routine = modname//':blend_states'

    !----------------------------------------------

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO jb=1, loopkit%nblks
      IF (jb /= loopkit%nblks) THEN
        nlen = nproma
      ELSE
        nlen = loopkit%npromz
      ENDIF
      DO jk=1, loopkit%nlev
        DO jc=1, nlen
          state_blnd(jc,jk,jb) = wfac_blnd(iblnd%left,jk) * state_left(jc,jk,jb) +  &
            &                    wfac_blnd(iblnd%right,jk) * state_right(jc,jk,jb)
         ENDDO  !jc
      ENDDO  !jk
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE blend_states

  !-----------------------------------------------------------------------------

  !>
  !! This subroutine multiplies a state with a factor.
  !!
  SUBROUTINE multiply_state ( state,  &  !inout
    &                         fac,    &  !in 
    &                         loopkit )  !in

    ! In/out variables
    REAL(wp),         INTENT(INOUT) :: state(:,:,:)  ! State to be multiplied with 'fac'
    REAL(wp),         INTENT(IN)    :: fac(:)        ! Currently only a vertically varying
                                                     ! factor is needed
    TYPE(t_loop_kit), INTENT(IN)    :: loopkit       ! Loop variables

    ! Local variables
    INTEGER :: jb, jc, jk, nlen  ! (jc is habitual placeholder for jc, je, jv)

    CHARACTER(LEN=*), PARAMETER :: &
      &       routine = modname//':multiply_state'

    !----------------------------------------------

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,nlen) ICON_OMP_DEFAULT_SCHEDULE
    DO jb=1, loopkit%nblks
      IF (jb /= loopkit%nblks) THEN
        nlen = nproma
      ELSE
        nlen = loopkit%npromz
      ENDIF
      DO jk=1, loopkit%nlev
        DO jc=1, nlen
          state(jc,jk,jb) = fac(jk) * state(jc,jk,jb) 
         ENDDO  !jc
      ENDDO  !jk
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE multiply_state

  !-----------------------------------------------------------------------------

  !>
  !! This subroutine applies a rudimentary sanity check to a state.
  !!
  SUBROUTINE sanity_check ( p_patch,    &  !in
    &                       state,      &  !in
    &                       bound,      &  !in
    &                       keys,       &  !in
    &                       lpassed,    &  !out
    &                       opt_slev,   &  !optin
    &                       opt_elev,   &  !optin
    &                       opt_message )  !optout

    ! In/out variables
    TYPE(t_patch),              INTENT(IN)  :: p_patch       ! Domain properties    
    REAL(wp),                   INTENT(IN)  :: state(:,:,:)  ! State to be checked
    REAL(wp),                   INTENT(IN)  :: bound         ! Lower/upper bound
    CHARACTER(LEN=*),           INTENT(IN)  :: keys          ! Comma-separated list with keywords
                                                             ! (e.g., keys = "cell,upper,abs")
    LOGICAL,                    INTENT(OUT) :: lpassed       ! .TRUE. -> sanity check passed
    INTEGER,          OPTIONAL, INTENT(IN)  :: opt_slev      ! Optional start level index
    INTEGER,          OPTIONAL, INTENT(IN)  :: opt_elev      ! Optional end level index
    CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: opt_message   ! Optional string with most critical 
                                                             ! state value and its location
                                                             ! (Make sure that LEN of the argument 
                                                             ! to opt_message at call is large enough)

    ! Local variables
    REAL(wp), ALLOCATABLE :: max_val_block(:)
    REAL(wp) :: buffer(4)
    REAL(wp) :: max_val, bound_eff
    INTEGER,  ALLOCATABLE :: max_loc_block(:,:) 
    INTEGER  :: max_loc_jb, max_loc_jc, max_loc_jk
    INTEGER  :: jg, jb, i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: rl_start, rl_end    
    INTEGER  :: i_startlev, i_endlev
    INTEGER  :: igrid, ibound, istate
    INTEGER  :: max_proc, mpi_comm
    INTEGER  :: istat 

    ! Identifiers and keywords
    INTEGER, PARAMETER :: ICELL  = 1
    INTEGER, PARAMETER :: IEDGE  = 2
    !-----------------------------------
    INTEGER, PARAMETER :: ILOWER = 1
    INTEGER, PARAMETER :: IUPPER = 2
    !-----------------------------------
    INTEGER, PARAMETER :: IASIS  = 1
    INTEGER, PARAMETER :: IABS   = 2
    !-----------------------------------
    INTEGER, PARAMETER :: IMAXSTATE = 1
    INTEGER, PARAMETER :: IMAXLON   = 2
    INTEGER, PARAMETER :: IMAXLAT   = 3
    INTEGER, PARAMETER :: IMAXJK    = 4
    !-----------------------------------
    INTEGER, PARAMETER :: KEYLEN = 10
    CHARACTER(LEN=*), PARAMETER :: key_cell  = 'cell'
    CHARACTER(LEN=*), PARAMETER :: key_edge  = 'edge'
    !-----------------------------------
    CHARACTER(LEN=*), PARAMETER :: key_lower = 'lower'
    CHARACTER(LEN=*), PARAMETER :: key_upper = 'upper'
    !-----------------------------------
    CHARACTER(LEN=*), PARAMETER :: key_abs   = 'abs'

    CHARACTER(len=*), PARAMETER :: routine = modname//':sanity_check'

    !----------------------------------------------

    ! - The content of this subroutine follows: 
    !   'src/atm_dyn_iconam/mo_nh_supervise: calculate_maxwinds'
    ! - The sanity check is limited to the prognostic part of the grid
    ! - Note: states defined on vertices are not supported!

    lpassed = .FALSE.

    ! Domain index
    jg = p_patch%id

    ! Check keyword list
    ! (The extrapolation is generally not run-time-critical, 
    ! so we can allow ourselves inefficient 
    ! but more user-friendly control via strings)

    ! Mandatory keys:
    ! state defined in cells or on edges?
    IF (INDEX(keys, key_cell) > 0) THEN
      igrid      = ICELL
      ! We use the oportunity, to set already the loop boundaries
      rl_start   = grf_bdywidth_c + 1
      rl_end     = min_rlcell_int
      i_startblk = p_patch%cells%start_block(rl_start) 
      i_endblk   = p_patch%cells%end_block(rl_end)  
    ELSEIF (INDEX(keys, key_edge) > 0) THEN
      igrid      = IEDGE
      rl_start   = grf_bdywidth_e + 1
      rl_end     = min_rledge_int
      i_startblk = p_patch%edges%start_block(rl_start) 
      i_endblk   = p_patch%edges%end_block(rl_end)  
    ELSE
      CALL finish(routine, 'Invalid or missing grid identifier in keys')
    ENDIF
    ! Lower or upper bound?
    IF (INDEX(keys, key_lower) > 0) THEN
      ! Check for state > bound <=> -state < -bound 
      ! (i.e. if state < bound <=> -state > -bound -> sanity check not passed)
      ibound    = ILOWER
      bound_eff = -bound
    ELSEIF (INDEX(keys, key_upper) > 0) THEN
      ! Check for state < bound
      ibound    = IUPPER
      bound_eff = bound
    ELSE
      CALL finish(routine, 'Invalid or missing bound identifier in keys')
    ENDIF

    ! Optional keys:
    IF (INDEX(keys, key_abs) > 0) THEN
      ! Check for +-|state| < +-bound
      istate = IABS
    ELSE
      ! Currently, everything else means, 
      ! to check the unmodified state values
      istate = IASIS
    ENDIF
    
    ! Range of grid layers to be checked
    IF (PRESENT(opt_slev)) THEN
      i_startlev = MIN(MAX(1, opt_slev), SIZE(state,2))
    ELSE
      i_startlev = 1
    ENDIF
    IF (PRESENT(opt_elev)) THEN
      i_endlev = MIN(MAX(i_startlev, opt_elev), SIZE(state,2))
    ELSE
      i_endlev = SIZE(state,2)
    ENDIF

    ALLOCATE( max_loc_block(2,i_startblk:i_endblk), &
      &       max_val_block(i_startblk:i_endblk),   &
      &       STAT=istat                            )
    IF (istat /= SUCCESS) CALL finish(routine, 'Allocation of max_loc_block and max_val_block failed')

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      IF (igrid == ICELL) THEN
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
      ELSEIF (igrid == IEDGE) THEN
        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end) 
      ENDIF
      ! (The following if-construct could be replaced by the shorter expression:
      ! max_loc_block(:,jb) = MAXLOC( sgn_bound * SIGN( ABS(state(i_startidx:i_endidx,i_startlev:i_endlev,jb)), &
      !   & MERGE( 1._wp, SIGN( 1._wp, state(i_startidx:i_endidx,i_startlev:i_endlev,jb) ), istate == IABS ) )
      ! for instance, where sgn_bound = -1._wp for ibound = ILOWER and sgn_bound = 1._wp for ibound = IUPPER. 
      ! However, we assume the if-construct to be computationally more efficient and plainer.)
      IF (ibound == ILOWER .AND. istate == IASIS) THEN
        ! To check, if state > bound <=> -state < -bound, we search for the location of max(-state)
        max_loc_block(:,jb) = MAXLOC(-state(i_startidx:i_endidx,i_startlev:i_endlev,jb))
        ! Given a field a(5:10) with the max. value a(5), 
        ! maxloc would not return 5, but 1.  
        ! For convenience we shift the indices:
        max_loc_block(1,jb) = max_loc_block(1,jb) + i_startidx - 1  ! jc or je
        max_loc_block(2,jb) = max_loc_block(2,jb) + i_startlev - 1  ! jk
        max_val_block(jb)   = -state(max_loc_block(1,jb),max_loc_block(2,jb),jb)
      ELSEIF (ibound == ILOWER .AND. istate == IABS) THEN
        ! To check, if |state| > bound <=> -|state| < -bound, we search for the location of max(-|state|)
        max_loc_block(:,jb) = MAXLOC(-ABS(state(i_startidx:i_endidx,i_startlev:i_endlev,jb)))
        max_loc_block(1,jb) = max_loc_block(1,jb) + i_startidx - 1
        max_loc_block(2,jb) = max_loc_block(2,jb) + i_startlev - 1
        max_val_block(jb)   = -ABS(state(max_loc_block(1,jb),max_loc_block(2,jb),jb))
      ELSEIF (ibound == IUPPER .AND. istate == IASIS) THEN
        ! To check, if state < bound, we search for the location of max(state)
        max_loc_block(:,jb) = MAXLOC(state(i_startidx:i_endidx,i_startlev:i_endlev,jb))
        max_loc_block(1,jb) = max_loc_block(1,jb) + i_startidx - 1  ! jc or je
        max_loc_block(2,jb) = max_loc_block(2,jb) + i_startlev - 1  ! jk
        max_val_block(jb)   = state(max_loc_block(1,jb),max_loc_block(2,jb),jb)
      ELSEIF (ibound == IUPPER .AND. istate == IABS) THEN
        ! To check, if |state| < bound, we search for the location of max(|state|)
        max_loc_block(:,jb) = MAXLOC(ABS(state(i_startidx:i_endidx,i_startlev:i_endlev,jb)))
        max_loc_block(1,jb) = max_loc_block(1,jb) + i_startidx - 1  ! jc or je
        max_loc_block(2,jb) = max_loc_block(2,jb) + i_startlev - 1  ! jk
        max_val_block(jb)   = ABS(state(max_loc_block(1,jb),max_loc_block(2,jb),jb))
      ENDIF
    ENDDO  !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! Next, we search for the max among all blocks
    max_loc_jb = MAXLOC(max_val_block(:), DIM=1)
    ! Index shift
    max_loc_jb = max_loc_jb + i_startblk - 1
    ! Just for convenience
    max_loc_jc = max_loc_block(1,max_loc_jb)  ! for opt_message
    max_loc_jk = max_loc_block(2,max_loc_jb)  ! for opt_message
    max_val    = max_val_block(max_loc_jb)

    ! Finally, we search for the max among all processes:
    ! Initialize max_proc with the ID of this process
    max_proc = get_my_mpi_work_id()
    ! Find global max and overwrite max_proc 
    ! with the ID of that process, where the global max is located (for opt_message)
    max_val = global_max(zfield=max_val, proc_id=max_proc)

    ! Sanity check passed?
    lpassed = max_val < bound_eff

    ! Get more information on the most critical state value and its location, if required
    IF (PRESENT(opt_message)) THEN

      ! Initialization
      buffer(:) = 999999._wp

      ! In order to use the preceding initialization 
      ! as a rudimentary "error" indicator for the broadcasting below, 
      ! we limit the assignment of buffer to the only process, which has to do it: 
      ! the process, where the global max is located
      IF (get_my_mpi_work_id() == max_proc) THEN
        ! Get the unmodified state value
        buffer(IMAXSTATE) = state(max_loc_jc,max_loc_jk,max_loc_jb)
        ! For the structured vertical grid the output of max_loc_jk is sufficient. 
        ! (In addition, without explicit input of this information, we cannot infer, 
        ! if state is defined on full or half levels. E.g., state may be defined 
        ! on half levels, but size(state,2) is not necessarily equal to nlevp1.)
        buffer(IMAXJK) = REAL(max_loc_jk, wp)
        ! For the unstructured horizontal grid it might be more usefull, 
        ! to get the actual position on the globe
        IF (igrid == ICELL) THEN
          buffer(IMAXLON) = p_patch%cells%center(max_loc_jc,max_loc_jb)%lon * rad2deg
          buffer(IMAXLAT) = p_patch%cells%center(max_loc_jc,max_loc_jb)%lat * rad2deg
        ELSEIF (igrid == IEDGE) THEN
          buffer(IMAXLON) = p_patch%edges%center(max_loc_jc,max_loc_jb)%lon * rad2deg
          buffer(IMAXLAT) = p_patch%edges%center(max_loc_jc,max_loc_jb)%lat * rad2deg
        ENDIF
      ENDIF  !IF (get_my_mpi_work_id() == max_proc)

      ! Now, the result has to be broadcast to all other processes, 
      ! which run this subroutine: 
      ! - sender:    max_proc
      ! - receivers: all processes under the communicator p_comm_work (hopefully)
      mpi_comm = get_my_mpi_work_communicator()
      CALL p_bcast(t_buffer=buffer, p_source=max_proc, comm=mpi_comm)

      ! Currently, a conversion of the information into a string for print 
      ! is sufficient for our purposes
      opt_message = 'Critical value(jg='//TRIM(int2string(jg))//     &
        &           ',lon='//TRIM(real2string(buffer(IMAXLON)))//    &
        &           ',lat='//TRIM(real2string(buffer(IMAXLAT)))//    &
        &           ',jk='//TRIM(int2string(NINT(buffer(IMAXJK))))// &
        &           ')='//TRIM(real2string(buffer(IMAXSTATE)))//'.'

    ENDIF  !IF (PRESENT(opt_message))

    DEALLOCATE( max_loc_block, &
      &         max_val_block, &
      &         STAT=istat     )
    IF (istat /= SUCCESS) CALL finish(routine, 'Deallocation of max_loc_block and max_val_block failed')

  END SUBROUTINE sanity_check

  !-----------------------------------------------------------------------------
  !----------------------------------------------------------------------------- 

  ! (Note: the following construction and deconstruction routines try to look for 
  ! some kind of misuse on a very rudimentary level, which makes them computationally 
  ! less efficient, but this is assumed to be bearable, since they are called only 
  ! during model setup)

  !---------------------------------------------------------------
  !             Construction routines for data types
  !---------------------------------------------------------------

  SUBROUTINE t_clim_state_construct(clim_state, cnstr)
    CLASS(t_clim_state), INTENT(INOUT) :: clim_state

    ! In/out variables
    TYPE(t_constructor_kit), INTENT(IN) :: cnstr

    ! Local variables
    INTEGER :: istat
    CHARACTER(LEN=*), PARAMETER :: &
      &       routine = modname//':t_clim_state_construct'

    !----------------------------------------------
    
    IF (.NOT. clim_state%linitialized) THEN 
      ALLOCATE( clim_state%temp(cnstr%nlev), &
        &       STAT=istat                   )
      IF (istat /= SUCCESS) CALL finish(routine, 'Allocation failed')
!$OMP PARALLEL 
      CALL init( clim_state%temp(:), lacc=.FALSE. )
!$OMP END PARALLEL
      clim_state%linitialized = .TRUE.
    ELSE
      CALL finish(routine, 'Attempt to initialize while already/still initialized')
    ENDIF
  END SUBROUTINE t_clim_state_construct

  !-----------------------------------------------------------------------------

  SUBROUTINE t_expol_metrics_construct(expol_metrics, cnstr)
    CLASS(t_expol_metrics), INTENT(INOUT) :: expol_metrics

    ! In/out variables
    TYPE(t_constructor_kit), INTENT(IN) :: cnstr

    ! Local variables
    INTEGER :: istat
    CHARACTER(LEN=*), PARAMETER :: &
      &       routine = modname//':t_expol_metrics_construct'

    !----------------------------------------------
    
    IF (.NOT. expol_metrics%linitialized) THEN 
      ALLOCATE( expol_metrics%zgpot_mc (cnstr%nlev),                  & 
        &       expol_metrics%zgpot_ifc(cnstr%nlev),                  &
        &       expol_metrics%wfac_blnd_mc (cnstr%nblnd, cnstr%nlev), &
        &       expol_metrics%wfac_blnd_ifc(cnstr%nblnd, cnstr%nlev), &
        &       STAT=istat                                            )  
      IF (istat /= SUCCESS) CALL finish(routine, 'Allocation failed')
!$OMP PARALLEL 
      CALL init( expol_metrics%zgpot_mc(:), lacc=.FALSE. )
      CALL init( expol_metrics%zgpot_ifc(:), lacc=.FALSE. )
      CALL init( expol_metrics%wfac_blnd_mc(:,:), lacc=.FALSE. )
      CALL init( expol_metrics%wfac_blnd_ifc(:,:), lacc=.FALSE. )
!$OMP END PARALLEL
      expol_metrics%linitialized = .TRUE.
    ELSE
      CALL finish(routine, 'Attempt to initialize while already/still initialized')
    ENDIF
  END SUBROUTINE t_expol_metrics_construct

  !-----------------------------------------------------------------------------

  SUBROUTINE t_expol_blending_state_construct(expol_blending_state, cnstr)
    CLASS(t_expol_blending_state), INTENT(INOUT) :: expol_blending_state

    ! In/out variables
    TYPE(t_constructor_kit), INTENT(IN) :: cnstr

    ! Local variables
    INTEGER :: istat
    CHARACTER(LEN=*), PARAMETER :: &
      &       routine = modname//':t_expol_blending_state_construct'

    !----------------------------------------------
    
    IF (.NOT. expol_blending_state%linitialized) THEN 
      ALLOCATE( expol_blending_state%temp(nproma, cnstr%nlev, cnstr%nblks_c),  &
        &       expol_blending_state%vn  (nproma, cnstr%nlev, cnstr%nblks_e),  &
        &       expol_blending_state%w   (nproma, cnstr%nlev, cnstr%nblks_c),  &
        &       expol_blending_state%qc  (nproma, cnstr%nlev, cnstr%nblks_c),  &
        &       expol_blending_state%qi  (nproma, cnstr%nlev, cnstr%nblks_c),  &
        &       expol_blending_state%qr  (nproma, cnstr%nlev, cnstr%nblks_c),  &
        &       expol_blending_state%qs  (nproma, cnstr%nlev, cnstr%nblks_c),  &
        &       STAT=istat                                                     )  
      IF (istat /= SUCCESS) CALL finish(routine, 'Allocation failed')
!$OMP PARALLEL 
      CALL init( expol_blending_state%temp(:,:,:), lacc=.FALSE. )
      CALL init( expol_blending_state%vn(:,:,:), lacc=.FALSE. )
      CALL init( expol_blending_state%w(:,:,:), lacc=.FALSE. )
      CALL init( expol_blending_state%qc(:,:,:), lacc=.FALSE. )
      CALL init( expol_blending_state%qi(:,:,:), lacc=.FALSE. )
      CALL init( expol_blending_state%qr(:,:,:), lacc=.FALSE. )
      CALL init( expol_blending_state%qs(:,:,:), lacc=.FALSE. )
!$OMP END PARALLEL
      expol_blending_state%linitialized = .TRUE.
    ELSE
      CALL finish(routine, 'Attempt to initialize while already/still initialized')
    ENDIF
  END SUBROUTINE t_expol_blending_state_construct

  !-----------------------------------------------------------------------------

  !---------------------------------------------------------------
  !            Deconstruction routines for data types
  !---------------------------------------------------------------

  ! (Note: the 'externals/fortran-support/src/mo_fortran_tools:DO_DEALLOCATE'-interface does not support all kind
  ! of fields we have allocated, so we fall back on the classical 'deallocate' here)

  SUBROUTINE t_clim_state_finalize(clim_state)
    CLASS(t_clim_state), INTENT(INOUT) :: clim_state

    ! Local variables
    INTEGER :: istat
    CHARACTER(LEN=*), PARAMETER :: &
      &       routine = modname//':t_clim_state_finalize'

    !----------------------------------------------

    IF (clim_state%linitialized) THEN
      DEALLOCATE( clim_state%temp, &
        &         STAT=istat       )
      IF (istat /= SUCCESS) CALL finish(routine, 'Deallocation failed')
      clim_state%linitialized = .FALSE.
    ENDIF
  END SUBROUTINE t_clim_state_finalize

  !-----------------------------------------------------------------------------

  SUBROUTINE t_expol_metrics_finalize(expol_metrics)
    CLASS(t_expol_metrics), INTENT(INOUT) :: expol_metrics

    ! Local variables
    INTEGER :: istat
    CHARACTER(LEN=*), PARAMETER :: &
      &       routine = modname//':t_expol_metrics_finalize'

    !----------------------------------------------

    IF (expol_metrics%linitialized) THEN
      DEALLOCATE( expol_metrics%zgpot_mc,      &
        &         expol_metrics%zgpot_ifc,     & 
        &         expol_metrics%wfac_blnd_mc,  &
        &         expol_metrics%wfac_blnd_ifc, &
        &         STAT=istat                   )
      IF (istat /= SUCCESS) CALL finish(routine, 'Deallocation failed')
      expol_metrics%linitialized = .FALSE.
    ENDIF
  END SUBROUTINE t_expol_metrics_finalize

  !-----------------------------------------------------------------------------

  SUBROUTINE t_expol_blending_state_finalize(expol_blending_state)
    CLASS(t_expol_blending_state), INTENT(INOUT) :: expol_blending_state

    ! Local variables
    INTEGER :: istat
    CHARACTER(len=*), PARAMETER :: &
      &       routine = modname//':t_expol_blending_state_finalize'

    !----------------------------------------------

    IF (expol_blending_state%linitialized) THEN
      DEALLOCATE( expol_blending_state%temp, &
        &         expol_blending_state%vn,   &
        &         expol_blending_state%w,    &
        &         expol_blending_state%qc,   &
        &         expol_blending_state%qi,   &
        &         expol_blending_state%qr,   &
        &         expol_blending_state%qs,   &
        &         STAT=istat                 )
      IF (istat /= SUCCESS) CALL finish(routine, 'Deallocation failed')
      expol_blending_state%linitialized = .FALSE.
    ENDIF
  END SUBROUTINE t_expol_blending_state_finalize

END MODULE mo_nh_vert_extrap_utils
