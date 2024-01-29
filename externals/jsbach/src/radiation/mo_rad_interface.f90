!> Contains the interfaces to the radiation processes
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"
!NEC$ options "-finline-file=externals/jsbach/src/radiation/mo_rad_process.pp-jsb.f90"
!NEC$ options "-finline-file=externals/jsbach/src/shared/mo_phy_schemes.pp-jsb.f90"
!NEC$ options "-finline-max-function-size=100"

MODULE mo_rad_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,   ONLY: t_jsb_process
  USE mo_jsb_task_class,      ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes SEB_, PHENO_, A2L_, RAD_, HYDRO_, ASSIMI_, CARBON_, VEG_

  ! Use process configurations
  dsl4jsb_Use_config(RAD_)
  dsl4jsb_Use_config(ASSIMI_)
  dsl4jsb_Use_config(SEB_)

  ! Use of process memories
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(RAD_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(CARBON_)
  dsl4jsb_Use_memory(VEG_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_rad_tasks

  !> Type definition for radiation_surface task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_surface_radiation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_surface_radiation
    PROCEDURE, NOPASS :: Aggregate => aggregate_surface_radiation
  END TYPE tsk_surface_radiation

  INTERFACE tsk_surface_radiation
    PROCEDURE Create_task_surface_radiation
  END INTERFACE tsk_surface_radiation

  !> Type definition for radiation_par task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_radiation_par
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_radiation_par
    ! R: Ich sehe keinen Sinn darin eine der Outputvariablen zu aggregieren, muss aber eine Prozedur bereitstellen
    PROCEDURE, NOPASS :: Aggregate => aggregate_radiation_par
  END TYPE tsk_radiation_par

  INTERFACE tsk_radiation_par
    PROCEDURE Create_task_radiation_par
  END INTERFACE tsk_radiation_par


  !> Type definition for albedo task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_albedo
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_albedo
    PROCEDURE, NOPASS :: Aggregate => aggregate_albedo
  END TYPE tsk_albedo

  INTERFACE tsk_albedo
    PROCEDURE Create_task_albedo
  END INTERFACE tsk_albedo
  
  !! quincy
  !> Type definition for sw_radiation_q_
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_sw_radiation_q_
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_sw_radiation_q_
    PROCEDURE, NOPASS :: Aggregate => aggregate_sw_radiation_q_
  END TYPE tsk_sw_radiation_q_

  INTERFACE tsk_sw_radiation_q_
    PROCEDURE Create_task_sw_radiation_q_
  END INTERFACE tsk_sw_radiation_q_

  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_interface'

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for surface_radiation task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "radiation"
  !!
  FUNCTION Create_task_surface_radiation(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_surface_radiation::return_ptr)
    CALL return_ptr%Construct(name='surface_radiation', process_id=RAD_, owner_model_id=model_id)

  END FUNCTION Create_task_surface_radiation

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for radiation_par task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "radiation"
  !!
  FUNCTION Create_task_radiation_par(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_radiation_par::return_ptr)
    CALL return_ptr%Construct(name='radiation_par', process_id=RAD_, owner_model_id=model_id)

  END FUNCTION Create_task_radiation_par

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for albedo task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "albedo"
  !!
  FUNCTION Create_task_albedo(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_albedo::return_ptr)
    CALL return_ptr%Construct(name='albedo', process_id=RAD_, owner_model_id=model_id)

  END FUNCTION Create_task_albedo

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for sw_radiation_q_ task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "sw_radiation_q_"
  !!
  FUNCTION Create_task_sw_radiation_q_(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_sw_radiation_q_::return_ptr)
    CALL return_ptr%Construct(name='sw_radiation_q_', process_id=RAD_, owner_model_id=model_id)

  END FUNCTION Create_task_sw_radiation_q_


  ! ================================================================================================================================
  !> Register tasks for radiation process
  !!
  !! @param[in,out] this      Instance of <PROCESS_NAME_LOWER_CASE> process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_rad_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER, INTENT(in) :: model_id

    CALL this%Register_task(tsk_surface_radiation(model_id))
    CALL this%Register_task(tsk_radiation_par    (model_id))
    CALL this%Register_task(tsk_albedo           (model_id))
    CALL this%Register_task(tsk_sw_radiation_q_  (model_id))

  END SUBROUTINE Register_rad_tasks

  ! ================================================================================================================================
  !>
  !> Implementation to update the surface radiation
  !! Task "surface_radiation" calculates the net radiations at the surface.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_surface_radiation(tile, options)

    USE mo_rad_process, ONLY: calc_radiation_surface_net
    USE mo_phy_schemes, ONLY: lwnet_from_lwdown

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options), INTENT(in) :: options

    ! Declare pointers to process configs and memory
    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(RAD_)

    ! Declare pointers for variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & swvis_srf_net,        &
      & swnir_srf_net,        &
      & sw_srf_net,           &
      & lw_srf_net,           &
      & rad_srf_net,          &
      & lw_srf_down,          &
      & swvis_srf_down,       &
      & swnir_srf_down,       &
      & alb_vis_lnd,          &
      & alb_nir_lnd,          &
      & fract_lice,           &
      & rad_net_lwtr,         &
      & rad_net_lice,         &
      & albedo_lwtr,          &
      & albedo_lice,          &
      & t,                    &
      & t_lwtr,               &
      & t_lice

    REAL(wp) ::                  &
      & lw_net_lwtr(options%nc), &
      & lw_net_lice(options%nc), &
      & sw_net_lwtr(options%nc), &
      & sw_net_lice(options%nc)

    TYPE(t_jsb_model), POINTER :: model

    INTEGER :: iblk , ics, ice, nc, ic

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_radiation'

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(RAD_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(SEB_)

    dsl4jsb_Get_var2D_onChunk(A2L_, swvis_srf_down) ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, swnir_srf_down) ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, lw_srf_down)    ! in

    dsl4jsb_Get_var2D_onChunk(RAD_, sw_srf_net)     ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, lw_srf_net)     ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, rad_srf_net)    ! out

    ! ---------------------------
    ! Go

    !> 1.0 lakes
    !!
    IF (tile%is_lake) THEN

      ! Water
      dsl4jsb_Get_var2D_onChunk(SEB_, t_lwtr)       ! in
      dsl4jsb_Get_var2D_onChunk(RAD_, albedo_lwtr)  ! in
      dsl4jsb_Get_var2D_onChunk(RAD_, rad_net_lwtr) ! out

      !$ACC DATA CREATE(sw_net_lwtr, lw_net_lwtr, sw_net_lice, lw_net_lice)

      CALL calc_radiation_surface_net( &
        & swvis_srf_down(:),           & ! in
        & swnir_srf_down(:),           & ! in
        & albedo_lwtr(:),              & ! in
        & albedo_lwtr(:),              & ! in
        & lw_srf_down(:),              & ! in
        & t_lwtr(:),                   & ! in
        & rad_net_lwtr(:),             & ! out
        & sw_net=sw_net_lwtr(:),       & ! out
        & lw_net=lw_net_lwtr(:)        & ! out
        & )

      ! Ice on lake
      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        dsl4jsb_Get_var2D_onChunk(SEB_, t_lice)       ! in
        dsl4jsb_Get_var2D_onChunk(SEB_, fract_lice)   ! in
        dsl4jsb_Get_var2D_onChunk(RAD_, albedo_lice)  ! in
        dsl4jsb_Get_var2D_onChunk(RAD_, rad_net_lice) ! out

        CALL calc_radiation_surface_net( &
          & swvis_srf_down(:),           & ! in
          & swnir_srf_down(:),           & ! in
          & albedo_lice(:),              & ! in
          & albedo_lice(:),              & ! in
          & lw_srf_down(:),              & ! in
          & t_lice(:),                   & ! in
          & rad_net_lice(:),             & ! out
          & sw_net=sw_net_lice(:),       & ! out
          & lw_net=lw_net_lice(:)        & ! out
          & )

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
        DO ic=1,nc
          rad_srf_net(ic) = (1._wp - fract_lice(ic)) * rad_net_lwtr(ic) + fract_lice(ic) * rad_net_lice(ic)
          sw_srf_net (ic) = (1._wp - fract_lice(ic)) * sw_net_lwtr(ic)  + fract_lice(ic) * sw_net_lice(ic)
          lw_srf_net (ic) = (1._wp - fract_lice(ic)) * lw_net_lwtr(ic)  + fract_lice(ic) * lw_net_lice(ic)
        END DO
        !$ACC END PARALLEL LOOP
      ELSE
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
        DO ic=1,nc
          rad_srf_net(ic) = rad_net_lwtr(ic)
          sw_srf_net (ic) = sw_net_lwtr (ic)
          lw_srf_net (ic) = lw_net_lwtr (ic)
        END DO
        !$ACC END PARALLEL LOOP
      END IF

      !$ACC END DATA


    !> 2.0 vegetation with use_quincy = .TRUE.
    !!
    ELSE IF (tile%is_vegetation .AND. model%config%use_quincy) THEN

      dsl4jsb_Get_var2D_onChunk(SEB_, t)              ! in

      ! quincy calculates all (albedo & short wave) but long wave radiation
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO ic=1,nc
        lw_srf_net(ic) = lwnet_from_lwdown(lw_srf_down(ic), t(ic))
      END DO
      !$ACC END PARALLEL LOOP
    
    !> 3.0 all but lakes 
    !! 
    !! for vegetation only when use_quincy = .FALSE.
    ELSE

      dsl4jsb_Get_var2D_onChunk(SEB_,      t)              ! in
      dsl4jsb_Get_var2D_onChunk(RAD_,      alb_vis_lnd)    ! in
      dsl4jsb_Get_var2D_onChunk(RAD_,      alb_nir_lnd)    ! in
      dsl4jsb_Get_var2D_onChunk(RAD_,      swvis_srf_net)  ! out
      dsl4jsb_Get_var2D_onChunk(RAD_,      swnir_srf_net)  ! out

      CALL calc_radiation_surface_net(  &
        & swvis_srf_down(:),            & ! in
        & swnir_srf_down(:),            & ! in
        & alb_vis_lnd(:),               & ! in
        & alb_nir_lnd(:),               & ! in
        & lw_srf_down(:),               & ! in
        & t(:),                         & ! in
        & rad_srf_net(:),               & ! out
        & swvis_net = swvis_srf_net(:), & ! out
        & swnir_net = swnir_srf_net(:), & ! out
        & sw_net    = sw_srf_net(:),    & ! out
        & lw_net    = lw_srf_net(:)     & ! out
        & )
    END IF

    !$ACC WAIT

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_surface_radiation

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "surface_radiation"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_surface_radiation(tile, options)

    TYPE(t_jsb_model), POINTER :: model
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(RAD_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_surface_radiation'

    INTEGER :: iblk , ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(RAD_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    ! these var are calculated by quincy for vegetation tiles
    IF (.NOT. tile%is_vegetation .OR. (tile%is_vegetation .AND. .NOT. model%config%use_quincy)) THEN
      dsl4jsb_Aggregate_onChunk(RAD_, rad_srf_net,    weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, sw_srf_net,     weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, swvis_srf_net,  weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, swnir_srf_net,  weighted_by_fract)
    ENDIF

    dsl4jsb_Aggregate_onChunk(RAD_, lw_srf_net,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(RAD_, rad_net_lwtr,   weighted_by_fract)
    IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
      dsl4jsb_Aggregate_onChunk(RAD_, rad_net_lice,   weighted_by_fract)
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_surface_radiation

  ! ================================================================================================================================
  !>
  !> Implementation to update the photosynthetically active radiation
  !! Task "radiation_par" calculates PAR and APAR.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_radiation_par(tile, options)

    ! Use declarations
    USE mo_rad_process,    ONLY: calc_par

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    INTEGER                           :: icanopy

    INTEGER                           :: iblk , ics, ice, nc, ic
    REAL(wp)                          :: dtime

    REAL(wp)                          :: canopy_bound_lai_tmp, canopy_bound_lai_delta
    REAL(wp), POINTER :: canopy_bound_lai(:)

    REAL(wp), DIMENSION(options%nc)   :: B4_layer_above ! R: cannot be scalar as it will be used as inout argum. for elemen. proce.
                                                        !    and thus has to have the same dimensions as the other out arguments.

    ! Declare pointers for process configuration and memory
    dsl4jsb_Def_config(RAD_)
    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(A2L_)

    TYPE(t_jsb_model), POINTER :: model
    INTEGER :: ncanopy               ! number of canopy layers
    !INTEGER, POINTER       :: acc_counter           ! counter for each delta_time step since the last variable output
    LOGICAL :: use_alb_veg_simple    !

    dsl4jsb_Real2D_onChunk :: swpar_srf_down
    dsl4jsb_Real2D_onChunk :: par                   ! Downward PAR flux in [W/m^2], Photosynth. Active Radiation (400-700nm)
    dsl4jsb_Real2D_onChunk :: par_down_mol          ! Downward PAR flux in mol (photons)/(m^2 s)
    ! dsl4jsb_Real2D_onChunk :: par_down_mol_nacc     ! PAR accumulated [MOL PHOTONS/M^2 ]
    ! dsl4jsb_Real2D_onChunk :: par_down_mol_tavg     ! PAR [MOL PHOTONS/M^2 S] (time mean value)
    dsl4jsb_Real2D_onChunk :: soil_reflectivity_par ! Soil reflectivity within the PAR spectrum [-]
    dsl4jsb_Real2D_onChunk :: alb_vis_soil          ! Soil albedo in the visible range
    dsl4jsb_Real2D_onChunk :: lai                   ! Leaf area index [-]
    dsl4jsb_Real2D_onChunk :: cos_zenith_angle      ! Angle of incoming radiation [-]
    dsl4jsb_Real2D_onChunk :: fract_par_direct      ! Fraction of direct (in contrast to diffuse) PAR [-]
    dsl4jsb_Real2D_onChunk :: fract_par_diffuse     ! Fraction of diffuse (in contrast to direct) PAR [-]
    ! dsl4jsb_Real2D_onChunk :: apar_nacc
    dsl4jsb_Real2D_onChunk :: faPAR

    dsl4jsb_Real3D_onChunk :: faPAR_cl              ! Fraction of absorbed PAR per leaf area for each canopy layer
    dsl4jsb_Real3D_onChunk :: apar_per_lai_cl       ! Absorbed PAR of canopy layer
    dsl4jsb_Real3D_onChunk :: lai_cl                ! Leaf area index of canopy layer
    ! dsl4jsb_Real3D_onChunk :: apar_per_lai_cl_nacc  ! Absorbed PAR accumulated [MOL PHOTONS/M^2 ]
    ! dsl4jsb_Real3D_onChunk :: apar_per_lai_cl_tavg  ! Absorbed PAR [MOL PHOTONS/M^2 S] (time mean value)

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_radiation_par'

    ! Get local variables from options argument
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime ! time step in seconds e.g. for one day=86400

    IF (.NOT. tile%is_vegetation) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get pointers to process configs and memory
    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(RAD_)
    dsl4jsb_Get_config(ASSIMI_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(A2L_)

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(RAD_)) RETURN

    ncanopy                     =  dsl4jsb_Config(ASSIMI_)%ncanopy
    use_alb_veg_simple          =  dsl4jsb_Config(RAD_)%use_alb_veg_simple
    canopy_bound_lai(0:ncanopy) => dsl4jsb_Config(ASSIMI_)% canopy_bound_lai(0:ncanopy)

    ! Set pointers to variables in memory
    dsl4jsb_Get_var2D_onChunk(A2L_,   cos_zenith_angle)      ! in (received from the atmosphere)
    dsl4jsb_Get_var2D_onChunk(A2L_,   fract_par_diffuse)     ! in (received from the atmosphere)
    dsl4jsb_Get_var2D_onChunk(A2L_,   swpar_srf_down)        ! in
    dsl4jsb_Get_var2D_onChunk(PHENO_, lai)                   ! in
    dsl4jsb_Get_var2D_onChunk(RAD_,   fract_par_direct)      ! out
    dsl4jsb_Get_var2D_onChunk(RAD_,   alb_vis_soil)          ! in
    dsl4jsb_Get_var2D_onChunk(RAD_,   par)                   ! in
    dsl4jsb_Get_var2D_onChunk(RAD_,   par_down_mol)          ! inout

    ! R: These variables could also be local as they are not used outside this procedure in the moment.
    !    However I want it in the std model output and thus it has to be in the memory
    ! dsl4jsb_Get_var2D_onChunk(RAD_,   par_down_mol_nacc)     ! inout
    ! dsl4jsb_Get_var2D_onChunk(RAD_,   par_down_mol_tavg)     ! inout
    dsl4jsb_Get_var2D_onChunk(RAD_,   soil_reflectivity_par) ! inout
    dsl4jsb_Get_var3D_onChunk(RAD_,   lai_cl)                ! inout
    dsl4jsb_Get_var3D_onChunk(RAD_,   faPAR_cl)              ! inout
    ! used outside for assimilation:
    dsl4jsb_Get_var3D_onChunk(RAD_,   apar_per_lai_cl)       ! inout
    ! dsl4jsb_Get_var3D_onChunk(RAD_,   apar_per_lai_cl_nacc)  ! inout
    ! dsl4jsb_Get_var3D_onChunk(RAD_,   apar_per_lai_cl_tavg)  ! inout
    ! dsl4jsb_Get_var2D_onChunk(RAD_,   apar_nacc)             ! inout

    dsl4jsb_Get_var2D_onChunk(RAD_,   faPAR)                 ! out

    ! ---------------------------
    ! Go

    !$ACC DATA CREATE(B4_layer_above)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1) COLLAPSE(2)
    DO icanopy=1,ncanopy
      DO ic=1,nc
        apar_per_lai_cl(ic,icanopy) = 0._wp
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
    DO ic=1,nc
      par_down_mol(ic)          = 0._wp ! for now also in the memory but is only for current time step/"delta_time"
      soil_reflectivity_par(ic) = 0._wp ! for now also in the memory but is only for current time step/"delta_time"
      
      par(ic) = swpar_srf_down(ic)

      ! Compute direct fraction of PAR.
      ! Ric  fract_par_direct = in ECHAM jsswdifpar, which is zsw_par_fract_diffuse = fract_diffuse_par  =sw_par_fract_direct
      fract_par_direct(ic) = 1._wp - fract_par_diffuse(ic)

      B4_layer_above(ic) = 7777777._wp ! Used to transfer variable B4 from one canopy layer to the next. Startvalue only to indicate
                                      ! first layer
      fapar(ic) = 0._wp
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC WAIT

    ! Compute faPAR_cl (absorbed photsynthetic active radiation for each canopy layer) and
    ! lai_cl (leaf area index for each canopy layer)
    DO icanopy=1,ncanopy
        canopy_bound_lai_tmp = canopy_bound_lai(icanopy)
        canopy_bound_lai_delta = canopy_bound_lai(icanopy) - canopy_bound_lai(icanopy-1)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR
       !NEC$ ivdep
        DO ic=1,nc
          CALL calc_par(                 &
            & ncanopy,                   & ! in
            & icanopy,                   & ! in
            & use_alb_veg_simple,        & ! in
            & par(ic),                    & ! in
            & alb_vis_soil(ic),           & ! in
            & lai(ic),                    & ! in
            & cos_zenith_angle(ic),       & ! in
            & fract_par_direct(ic),       & ! in
            & canopy_bound_lai_tmp,       & ! in
            & canopy_bound_lai_delta,    & ! in
            & B4_layer_above(ic),         & ! inout
            & par_down_mol(ic),           & ! out, used outside for assimilation
            & soil_reflectivity_par(ic),  & ! out, R could also be local as it is not used outside in the moment,
                                            ! however I want it in the model output and thus it has to be updated by returning it
            & faPAR_cl(ic,icanopy),       & ! out, R could also be local as it is not used outside in the moment,
                                            ! however I want it in the model output and thus it has to be updated by returning it
            & lai_cl(ic,icanopy),         & ! out, R could also be local as it is not used outside in the moment,
                                            ! however I want it in the model output and thus it has to be updated by returning it
            & apar_per_lai_cl(ic,icanopy) ) ! out, used outside for assimilation
      END DO
      !$ACC END PARALLEL LOOP
    END DO

    ! R: Note, the following code is preliminary and should be replaced later by an laccu-procedure.
    !    For the correct integration over time as it is done now, this task must be called only once per timestep.
    ! Accumulate incoming par in mol(photons)/(m^2) over time
    ! par_down_mol_nacc(:) = par_down_mol_nacc(:) + par_down_mol(:) * dtime

    ! R: This whole workaround should get abitary soon...
    !acc_counter = acc_counter + 1
    !par_down_mol_tavg = par_down_mol_nacc / dtime / acc_counter
    ! apar_per_lai_cl_nacc(:,:)  = apar_per_lai_cl_nacc(:,:) + apar_per_lai_cl(:,:) * dtime

    ! R: This whole workaround (acc_counter) should get abitary soon...
    !apar_per_lai_cl_tavg(:,icanopy)  = apar_per_lai_cl_nacc(:,icanopy) / dtime / acc_counter
    ! DO icanopy=1,ncanopy
      ! apar_nacc(:) = apar_nacc(:) + par_down_mol(:) * faPAR_cl(:,icanopy) * dtime
    ! END DO

    ! Create faPAR
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
    DO icanopy=1,ncanopy
      !$ACC LOOP SEQ
      DO ic=1,nc
        faPAR(ic) = faPAR(ic) + faPAR_cl(ic,icanopy) ! faPAR_cl is always the fraction relativ to incoming radiation on the top of the
                                                     ! whole canopy
      END DO
    END DO
    !$ACC END PARALLEL LOOP


!        ! Compute absorbed PAR per leaf area in canopy layer [units: (absorbed photons) / (m^2(leaf area) s)] from
!        ! par and fraction of absorbed PAR
!        apar_per_lai_cl(:,icanopy) = par_down_mol(:) * faPAR_cl(:,icanopy) / (MAX(lai_cl(:,icanopy),1.e-10_wp))
!        ! R: Folgendes fällt weg:
!        ! mask(1:nidx,itile) <--theLand%Surface%is_vegetation(kidx0:kidx1,1:ntiles)
!        ! apar_per_lai_cl(:,icanopy) = MERGE(apar_per_lai_cl(:,icanopy),0._dp,mask(1:nidx,itile))
!        ! Da damit die Wüstenflächen auf 0 gesetzt werden und JSBACH4 diese Prozedur dann
!        ! gar nicht erst aufruft für diese tile.
!        IF (waterLimitationFlag) THEN ! to make sure to accumulate only once per timestep
!           apar_per_lai_cl_nacc(:,icanopy)  = apar_per_lai_cl_nacc(:,icanopy) + apar_per_lai_cl(:,icanopy) * dtime
!           ! R: This whole workaround (acc_counter) should get abitary soon...
!           !apar_per_lai_cl_tavg(:,icanopy)  = apar_per_lai_cl_nacc(:,icanopy) / dtime / acc_counter
!           apar_nacc(:)                     = apar_nacc(:) + par_down_mol(:) * faPAR_cl(:,icanopy) * dtime
!        END IF

    !$ACC END DATA

    !$ACC WAIT

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_radiation_par

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "radiation_par"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_radiation_par(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_radiation_par'

    INTEGER :: iblk !, ics, ice

    iblk = options%iblk
    !ics  = options%ics
    !ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_radiation_par

  ! ================================================================================================================================
  !>
  !> Implementation to update the albedo
  !! Task "albedo" calculates several albedos.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_albedo(tile, options)

    ! Use declarations
    USE mo_jsb_lctlib_class,  ONLY: t_lctlib_element
    USE mo_jsb_grid_class,    ONLY: t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_vgrid
    USE mo_rad_process,       ONLY: &
      & get_surface_albedo_simple, calc_soil_albedo, calc_alb_lwtr, calc_alb_lice, calc_glacier_albedo, &
      & calc_sky_view_fractions, calc_snow_albedo, Merge_albedos_of_vegtile, &
      & Has_minimal_radiation
    USE mo_jsb_time,          ONLY: is_time_experiment_start, is_time_ltrig_rad_m1
    USE mo_jsb_control,       ONLY: jsbach_runs_standalone

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER      :: model
    TYPE(t_lctlib_element), POINTER :: lctlib       !< land-cover-type library - parameter across pft's
    INTEGER                         :: iblk, ics, ice, nc, ic
    LOGICAL, ALLOCATABLE            :: l_day(:)

    CHARACTER(len=*), PARAMETER     :: routine = modname//':update_albedo'

    ! Declare pointers to process configuration and memory
    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_config(RAD_)

    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(CARBON_)

    ! Declare configuration variables
    CHARACTER(len=:), ALLOCATABLE :: use_alb_soil_organic_C
    INTEGER :: use_alb_soil_organic_int

    REAL(wp) ::            &
      & albedo_age_weight, &    ! 0 < albedo_age_weight < 1: snow albedo is calculated by linearly weighting
                                !                            the snow albedo resulting from both schemes
      & AlbedoCanopySnow

    LOGICAL ::                     &
      & use_alb_veg_simple,        &
      & use_alb_soil_scheme,       &
      & use_alb_soil_litter,       &
      & use_alb_mineralsoil_const

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk ::  &
      & cos_zenith_angle,      &
      & t_unfilt,              &
      & alb_vis,               &
      & alb_nir,               &
      & alb_vis_lnd,           &
      & alb_nir_lnd,           &
      & alb_vis_can,           &
      & alb_nir_can,           &
      & alb_vis_soil,          &
      & alb_nir_soil,          &
      & alb_vis_snow,          &
      & alb_nir_snow,          &
      & alb_background,        &
      & alb_vis_mineralsoil,   &
      & alb_nir_mineralsoil,   &
      & sky_view_fract,        &
      & sky_view_fract_stem,   &
      & snow_age,              &
      & lai,                   &
      & fract_fpc_max,         &
      & fract_forest,          &
      & fract_snow_can,        &
      & fract_snow_soil,       &
      & c_ag_sum_1,            &
      !CBAL & c_slow,           &
      & c_bg_sum,              &
      & swvis_srf_down,        &
      & swnir_srf_down
    ! Lakes
    dsl4jsb_Real2D_onChunk ::  &
      & albedo_lwtr,           &
      & fract_lice,            &
      & t_lice,                &
      & w_snow_lice,           &
      & albedo_lice

    REAL(wp) :: &
      & specificLeafArea_C_param, &
      & AlbedoLitterVIS_param,    &
      & StemArea_param
    LOGICAL :: ForestFlag_param

    ! Set local variables from options argument
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Set model
    model => Get_model(tile%owner_model_id)

    ! Update albedo only at experiment start, if model runs standalone or one time step before the radiation time step
    IF (.NOT. is_time_experiment_start(options%current_datetime) .AND. &
      & .NOT. is_time_ltrig_rad_m1(options%current_datetime, options%dtime, model%id) .AND. &
      & .NOT. jsbach_runs_standalone()) RETURN

    ! Enable for debugging:
    ! IF (iblk == 1 .AND. .NOT. ASSOCIATED(tile%next_sibling_tile)) THEN
    !   IF (is_time_ltrig_rad_m1(options%current_datetime, options%dtime, model%id) .OR. &
    !     & is_time_experiment_start(options%current_datetime)) THEN
    !     CALL message('update_albedo ', 'compute surface albedo')
    !   END IF
    ! END IF

    ! Set pointers to process configs and memory
    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_config(RAD_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)
    IF (tile%contains_vegetation) THEN
      dsl4jsb_Get_memory(PHENO_)
    END IF

    ! Set pointers to variables in config
    use_alb_veg_simple        = dsl4jsb_Config(RAD_) %use_alb_veg_simple
    use_alb_soil_scheme       = dsl4jsb_Config(RAD_) %use_alb_soil_scheme
    use_alb_mineralsoil_const = dsl4jsb_Config(RAD_) %use_alb_mineralsoil_const
    use_alb_soil_litter       = dsl4jsb_Config(RAD_) %use_alb_soil_litter
    use_alb_soil_organic_C    = dsl4jsb_Config(RAD_) %use_alb_soil_organic_C
    albedo_age_weight         = dsl4jsb_Config(RAD_) %albedo_age_weight
    AlbedoCanopySnow          = dsl4jsb_Config(RAD_) %AlbedoCanopySnow

    IF(use_alb_soil_organic_C .eq. 'linear') THEN
       use_alb_soil_organic_int = 1
    ELSE IF(use_alb_soil_organic_C .eq. 'log') THEN
       use_alb_soil_organic_int = 2
    END IF
  
    c_bg_sum => NULL()
    IF (tile%is_bare .OR. tile%is_vegetation) THEN
      IF (TRIM(use_alb_soil_organic_C) /= '' .AND. .NOT. tile%Is_process_active(CARBON_)) THEN
        CALL finish(TRIM(routine), 'Carbon process must be active for use_alb_soil_organic_C')
      END IF

      IF (TRIM(use_alb_soil_organic_C) /= '') THEN
        dsl4jsb_Get_memory(CARBON_)
        dsl4jsb_Get_var2D_onChunk(CARBON_, c_bg_sum)  ! in
      ELSE
        ALLOCATE(c_bg_sum(nc))
        !$ACC ENTER DATA CREATE(c_bg_sum)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
        DO ic=1,nc
          c_bg_sum(ic) = 0._wp                           ! dummy
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END IF

    ! Set pointers to variables in memory
    dsl4jsb_Get_var2D_onChunk(A2L_,   cos_zenith_angle)      ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   swvis_srf_down)        ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   swnir_srf_down)        ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,   t_unfilt)              ! in

    !RMS dsl4jsb_Get_var2D_onChunk(CARBON_,   c_slow)            ! in

    dsl4jsb_Get_var2D_onChunk(RAD_,   alb_background)        ! in
    dsl4jsb_Get_var2D_onChunk(RAD_,   alb_vis)               ! inout
    dsl4jsb_Get_var2D_onChunk(RAD_,   alb_nir)               ! inout
    dsl4jsb_Get_var2D_onChunk(RAD_,   alb_vis_snow)          ! inout
    dsl4jsb_Get_var2D_onChunk(RAD_,   alb_nir_snow)          ! inout
    IF (.NOT. tile%is_lake) THEN
      dsl4jsb_Get_var2D_onChunk(RAD_, alb_vis_lnd)           ! inout
      dsl4jsb_Get_var2D_onChunk(RAD_, alb_nir_lnd)           ! inout
    END IF

    ! ---------------------------
    ! Go

    ALLOCATE(l_day(nc))
    !$ACC ENTER DATA CREATE(l_day)

    !> 0.9  Set day-or-night-switch
    !! 
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR
    DO ic=1,nc
      l_day(ic) = Has_minimal_radiation(swvis_srf_down(ic), swnir_srf_down(ic))
    END DO
    !$ACC END PARALLEL LOOP

    !> 1.0 snow
    !!
    !! Calculate snow albedos. Used later only for bare and vegetation tile types
    !! Note, for vegetation tiles this snow albedos are used only for snow on soil
    IF (tile%is_bare .OR. tile%is_vegetation) THEN

      dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow_soil)  ! in
      dsl4jsb_Get_var2D_onChunk(RAD_,   snow_age)         ! in

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO ic=1,nc
        CALL calc_snow_albedo(    &
          & l_day(ic),            & ! in
          & fract_snow_soil(ic),  & ! in
          & snow_age(ic),         & ! in
          & cos_zenith_angle(ic), & ! in
          & albedo_age_weight,    & ! in
          & t_unfilt(ic),         & ! in
          & alb_vis_snow(ic),     & ! inout
          & alb_nir_snow(ic)      & ! inout
          & )
      END DO
      !$ACC END PARALLEL LOOP
  
    END IF

    ! TODO Compute alb_vis and alb_nir separately!
    ! NOTE This is done by the QUINCY routine calc_sw_srf_net() when use_quincy = .TRUE. (in another Task)

    !> 2.0 lakes & lakes with ice
    !!
    !!
    IF (tile%is_lake) THEN
      dsl4jsb_Get_var2D_onChunk(SEB_,   t_lice)       ! in
      dsl4jsb_Get_var2D_onChunk(RAD_,   albedo_lwtr)  ! inout

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO ic=1,nc
        CALL calc_alb_lwtr( &
          & l_day(ic),       & ! in
          & albedo_lwtr(ic)  & ! inout
          & )
      END DO
      !$ACC END PARALLEL LOOP

      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        dsl4jsb_Get_var2D_onChunk(SEB_,   fract_lice)   ! in
        dsl4jsb_Get_var2D_onChunk(HYDRO_, w_snow_lice)  ! in
        dsl4jsb_Get_var2D_onChunk(RAD_,   albedo_lice)  ! inout

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
        DO ic=1,nc
          CALL calc_alb_lice(  &
            & l_day(ic),       & ! in
            & t_lice(ic),      & ! in
            & w_snow_lice(ic), & ! in
            & albedo_lice(ic)  & ! inout
            & )
        END DO
        !$ACC END PARALLEL LOOP

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR
        DO ic = 1, nc
          IF(l_day(ic)) THEN
            alb_vis(ic) = (1._wp - fract_lice(ic)) * albedo_lwtr(ic) + fract_lice(ic) * albedo_lice(ic)
            alb_nir(ic) = alb_vis(ic)  ! TODO
          END IF
        END DO
        !$ACC END PARALLEL LOOP
      ELSE
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR
        DO ic = 1, nc
          IF(l_day(ic)) THEN
            alb_vis(ic) = albedo_lwtr(ic)
            alb_nir(ic) = alb_vis(ic)  ! TODO
          END IF
        END DO
        !$ACC END PARALLEL LOOP
      END IF

    !> 3.0 glacier
    !!
    ELSE IF (tile%is_glacier) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO ic=1,nc
        CALL calc_glacier_albedo(  &
          & l_day(ic),             & ! in
          & t_unfilt(ic),          & ! in
          & alb_vis_lnd(ic),       & ! inout
          & alb_nir_lnd(ic)        & ! inout
          & )
      END DO
      !$ACC END PARALLEL LOOP

    !> 4.0 soil
    !! currently jsbach4 does not consider "bare tiles"
    !!
    ELSE IF (tile%is_bare) THEN
      dsl4jsb_Get_var2D_onChunk(RAD_, alb_vis_mineralsoil)   ! inout
      dsl4jsb_Get_var2D_onChunk(RAD_, alb_nir_mineralsoil)   ! inout
      dsl4jsb_Get_var2D_onChunk(RAD_, alb_vis_soil)          ! inout
      dsl4jsb_Get_var2D_onChunk(RAD_, alb_nir_soil)          ! inout

      !!       needs the active CARBON_ process
      !!       by default: use_alb_soil_scheme = .FALSE.
      !!       also the soil albedo is calc below for vegetation tiles
      IF (use_alb_soil_scheme) THEN
        IF (tile%lcts(1)%lib_id == 0) CALL finish(TRIM(routine), 'lctlib not available on tile')

        dsl4jsb_Get_var2D_onChunk(CARBON_, c_ag_sum_1) ! in
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
        !NEC$ ivdep
        DO ic=1,nc
          CALL calc_soil_albedo(                        &
            & l_day(ic),                                & ! in
            & use_alb_mineralsoil_const,                & ! in
            & use_alb_soil_organic_int,                 & ! in
            & use_alb_soil_litter,                      & ! in
            & c_ag_sum_1(ic),                           & ! in
            !CBAL & c_slow(ic),                           & ! in
            & c_bg_sum(ic),                             & ! in
            ! @todo on a bare soil tile, there is no vegetation and therefore also no lctlib entry!?
            & dsl4jsb_Lctlib_param(specificLeafArea_C), & ! in
            & dsl4jsb_Lctlib_param(AlbedoLitterVIS),    & ! in
            & alb_background(ic),                       & ! in
            & alb_vis_mineralsoil(ic),                  & ! inout
            & alb_nir_mineralsoil(ic),                  & ! inout
            & alb_vis_soil(ic),                         & ! inout
            & alb_nir_soil(ic)                          & ! inout
            & )
        END DO
        !$ACC END PARALLEL LOOP
      ELSE
        ! Do nothing: alb_vis_soil and alb_nir_soil are already in memory initialited from bc_land_phy.nc file
        ! and represent the soil albedo including carbon in the soil as well as litter on the soil, but without
        ! the canopy above.
        IF (debug_on() .AND. iblk == 1) THEN
          CALL message(TRIM(routine), 'Soil albedo scheme not used for bare soil tiles but taken from bc_land_phy.nc file')
        END IF
      END IF

      ! Merge the albedo of soil with the snow on it for the overall bare soil tile albedo
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO ic=1,nc
        alb_vis_lnd(ic) = (fract_snow_soil(ic) * alb_vis_snow(ic)  + (1.0_wp - fract_snow_soil(ic)) * alb_vis_soil(ic))
        alb_nir_lnd(ic) = (fract_snow_soil(ic) * alb_nir_snow(ic)  + (1.0_wp - fract_snow_soil(ic)) * alb_nir_soil(ic))
      END DO
      !$ACC END PARALLEL LOOP

    !> 5.0 vegetation / canopy
    !!  jsbach4 routines
    !!
    ELSE IF (tile%is_vegetation .AND. .NOT. model%config%use_quincy) THEN

      dsl4jsb_Get_var2D_onChunk(PHENO_,    lai)                   ! in
      dsl4jsb_Get_var2D_onChunk(PHENO_,    fract_fpc_max)         ! in
      IF (use_alb_veg_simple) THEN 
        dsl4jsb_Get_var2D_onChunk(PHENO_,    fract_forest)          ! in
      ENDIF
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow_can)        ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow_soil)       ! in
      dsl4jsb_Get_var2D_onChunk(RAD_,      alb_vis_mineralsoil)   ! inout
      dsl4jsb_Get_var2D_onChunk(RAD_,      alb_nir_mineralsoil)   ! inout
      dsl4jsb_Get_var2D_onChunk(RAD_,      alb_vis_can)           ! in (jsbach4: from rad_init | quincy: from calc_sw_srf_net())
      dsl4jsb_Get_var2D_onChunk(RAD_,      alb_nir_can)           ! in (jsbach4: from rad_init | quincy: from calc_sw_srf_net())
      dsl4jsb_Get_var2D_onChunk(RAD_,      alb_vis_soil)          ! inout
      dsl4jsb_Get_var2D_onChunk(RAD_,      alb_nir_soil)          ! inout


      !> 5.1  simple calculation of alb_vis_lnd / alb_nir_lnd
      !!
      !! TRUE for jsbach_lite, FALSE for jsbach_pfts
      IF (use_alb_veg_simple) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
        DO ic=1,nc
          alb_vis_lnd(ic) = get_surface_albedo_simple( &
            & t_unfilt(ic),                            &
            & fract_forest(ic),                        &
            & fract_snow_soil(ic),                     &
            & fract_snow_can(ic),                      &
            & lai(ic),                                 &
            & alb_background(ic)                       & ! Note, alb_background <-- bc_land_phy.nc
                        )
          alb_nir_lnd(ic) = alb_vis_lnd(ic)
        END DO
        !$ACC END PARALLEL LOOP

      !> 5.2  detailed calculation of alb_vis_lnd & alb_nir_lnd
      !!      
      !! merging albedo of snow, soil and canopy
      ELSE

        !> 5.2.1 Get bare soil albedo for vegetation tiles
        !!       needs the active CARBON_ process
        !!       by default: use_alb_soil_scheme = .FALSE.
        !!       also the soil albedo is calc above with the same routine in case "IF (tile%is_bare)"
        IF (use_alb_soil_scheme) THEN
          ! The albedo for bare soil under the canopy of a vegetation tile is the same as for the bare soil tile:

          IF (tile%lcts(1)%lib_id == 0) CALL finish(TRIM(routine), 'lctlib not available on tile')

          dsl4jsb_Get_var2D_onChunk(CARBON_, c_ag_sum_1) ! in (Note: c_humus_1 is part of green litter.)

          specificLeafArea_C_param = dsl4jsb_Lctlib_param(specificLeafArea_C)
          AlbedoLitterVIS_param    = dsl4jsb_Lctlib_param(AlbedoLitterVIS)
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
          !NEC$ ivdep
          DO ic=1,nc
            CALL calc_soil_albedo(                        &
              & l_day(ic),                                & ! in
              & use_alb_mineralsoil_const,                & ! in
              & use_alb_soil_organic_int,                 & ! in
              & use_alb_soil_litter,                      & ! in
              & c_ag_sum_1(ic),                           & ! in
              !CBAL & c_slow(ic),                          & ! in
              & c_bg_sum(ic),                             & ! in
              & specificLeafArea_C_param,                 & ! in
              & AlbedoLitterVIS_param,                    & ! in
              & alb_background(ic),                       & ! in
              & alb_vis_mineralsoil(ic),                  & ! inout
              & alb_nir_mineralsoil(ic),                  & ! inout
              & alb_vis_soil(ic),                         & ! inout
              & alb_nir_soil(ic)                          & ! inout
              & )
          END DO
          !$ACC END PARALLEL LOOP
        ELSE
          ! Do nothing: alb_vis_soil and alb_nir_soil are already in memory initialited from the bc_land_phy.nc file
          ! and represent the soil albedo including carbon in the soil as well as litter on the soil, but without
          ! the canopy above.

          IF (debug_on() .AND. iblk == 1) THEN
            CALL message(TRIM(routine), 'Soil albedo scheme not used for veg tiles but taken from bc_land_phy.nc file')
          END IF
        END IF

        !$ACC WAIT

        ! Get the fractions of the vegetation tile that is bare soil from the sky point of view
        ! Needed later for merging the contributing snow-, bare soil- and canopy-albedos for this vegetation tile
        dsl4jsb_Get_var2D_onChunk(RAD_, sky_view_fract)       ! inout
        dsl4jsb_Get_var2D_onChunk(RAD_, sky_view_fract_stem)  ! inout

        ! double check if tile%is_vegetation is PFT
        IF (tile%lcts(1)%lib_id == 0) CALL finish(TRIM(routine), 'lctlib not available on tile')

        !> 5.2.2 calc sky-view fractions
        !!
        StemArea_param = dsl4jsb_Lctlib_param(StemArea)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
        DO ic=1,nc
          CALL calc_sky_view_fractions(       &
            & l_day(ic),                      & ! in
            & lai(ic),                        & ! in
            & fract_fpc_max(ic),              & ! in
            & cos_zenith_angle(ic),           & ! in
            & StemArea_param,                 & ! in (Area of stems and branches of woody plants.)
            & sky_view_fract(ic),             & ! inout
            & sky_view_fract_stem(ic)         & ! inout
            & )
        END DO
        !$ACC END PARALLEL LOOP

        !> 5.2.3 merge albedos of snow, soil and canopy into land
        !! 
        ForestFlag_param = dsl4jsb_Lctlib_param(ForestFlag)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
        DO ic=1,nc
          CALL Merge_albedos_of_vegtile(        &
            & l_day(ic),                        & ! in
            & ForestFlag_param,                 & ! in  ("true" for forest landcover types)
            & sky_view_fract(ic),               & ! in
            & sky_view_fract_stem(ic),          & ! in
            & AlbedoCanopySnow,                 & ! in
            & fract_snow_soil(ic),              & ! in
            & fract_snow_can(ic),               & ! in
            & alb_vis_snow(ic),                 & ! in
            & alb_nir_snow(ic),                 & ! in
            & alb_vis_soil(ic),                 & ! in
            & alb_nir_soil(ic),                 & ! in
            & alb_vis_can(ic),                  & ! in
            & alb_nir_can(ic),                  & ! in
            & alb_vis_lnd(ic),                  & ! inout
            & alb_nir_lnd(ic)                   & ! inout
            & )
        END DO
        !$ACC END PARALLEL LOOP

      END IF ! use_alb_veg_simple

    END IF   ! IF (tile%is_vegetation .AND. .NOT. model%config%use_quincy)
    !$ACC WAIT

    !> 6.0 pass the "land albedo" to the "tile albedo" 
    !! this is needed because for lake tiles, the "tile albedo" is calculated directly, \n
    !! but for vegetation, bare and glacier, the "land albedo" is calculated \n
    !! with "use_quincy = T", for vegetation tiles the "tile albedo" (alb_vis & alb_nir) is calculated directly
    IF (tile%is_bare .OR. tile%is_glacier .OR. (tile%is_vegetation .AND. .NOT. model%config%use_quincy)) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO ic=1,nc
        alb_vis(ic) = alb_vis_lnd(ic)
        alb_nir(ic) = alb_nir_lnd(ic)
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    !> 7.0 deallocate a pointer
    !! 
    IF (TRIM(use_alb_soil_organic_C) == '') THEN
      IF (ASSOCIATED(c_bg_sum)) THEN
        !$ACC EXIT DATA DELETE(c_bg_sum)
        DEALLOCATE(c_bg_sum)
      END IF
    END IF
    !$ACC EXIT DATA DELETE(l_day)
    DEALLOCATE(l_day)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_albedo

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "albedo"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_albedo(tile, options)

    USE mo_jsb_time,    ONLY: is_time_experiment_start, is_time_ltrig_rad_m1
    USE mo_jsb_control, ONLY: jsbach_runs_standalone

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(RAD_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_albedo'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Update albedo only at experiment start, if model runs standalone or one time step before the radiation time step
    IF (.NOT. is_time_experiment_start(options%current_datetime) .AND. &
      & .NOT. is_time_ltrig_rad_m1(options%current_datetime, options%dtime, model%id) .AND. &
      & .NOT. jsbach_runs_standalone()) RETURN

    dsl4jsb_Get_memory(RAD_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    ! for vegetation tiles, already aggregated by aggregate_sw_radiation_q_()
    IF (tile%is_vegetation .AND. .NOT. model%config%use_quincy) THEN
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_lnd,         weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_lnd,         weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_soil,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_soil,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_mineralsoil, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_mineralsoil, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_snow,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_snow,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis,             weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir,             weighted_by_fract)
    ! for all other tiles (the task sw_radiation_q_() runs for vegetation tiles only)
    ELSE IF (.NOT. tile%is_vegetation) THEN
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_lnd,         weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_lnd,         weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_soil,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_soil,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_mineralsoil, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_mineralsoil, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_snow,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_snow,        weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis,             weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir,             weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, albedo_lwtr,         weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, albedo_lice,         weighted_by_fract)
    ENDIF
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_albedo


  !-----------------------------------------------------------------------------------------------------
  !> calculates surface radiation reflection and absorption
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_sw_radiation_q_(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: A2L_, RAD_, VEG_, PHENO_, HYDRO_
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid
    !USE mo_rad_constants,         ONLY: rad2ppfd, def_alb_vis_soil, def_alb_nir_soil, rfr_ratio_toc, k_r2fr_chl
    USE mo_rad_process,           ONLY: calc_sw_srf_net, calc_albedo_radiation_quincy_hlp, &
                                        Has_minimal_radiation, calc_snow_albedo

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(RAD_)
    ! Use of process memories
    dsl4jsb_Use_memory(A2L_)   
    dsl4jsb_Use_memory(RAD_)   
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(PHENO_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(SEB_)
  

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER        :: model
    TYPE(t_lctlib_element), POINTER        :: lctlib              !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER        :: vgrid_canopy        ! Vertical grid
    LOGICAL,  DIMENSION(options%nc)        :: l_day               ! True = day | False = night
    REAL(wp), DIMENSION(options%nc)        :: alb_vis_snow_soil   ! albedo of snow and soil combined 
    REAL(wp), DIMENSION(options%nc)        :: alb_nir_snow_soil   ! albedo of snow and soil combined 
    INTEGER                                :: ncanopy             ! number of canopy layers
    INTEGER                                :: icanopy             !< loop 1,ncanopy
    INTEGER                                :: iblk, ics, ice, nc, ic
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sw_radiation_q_'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(RAD_)   
    dsl4jsb_Def_memory(A2L_)   
    dsl4jsb_Def_memory(RAD_)   
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SEB_)
    ! Declare pointers to variables in memory
    ! A2L_
    dsl4jsb_Real2D_onChunk :: swpar_srf_down
    dsl4jsb_Real2D_onChunk :: fract_par_diffuse     ! jsb4 fract_par_diffuse == quincy swpar_fdiffuse
    dsl4jsb_Real2D_onChunk :: cos_zenith_angle
    dsl4jsb_Real2D_onChunk :: swnir_srf_down
    dsl4jsb_Real2D_onChunk :: swvis_srf_down
    ! RAD_ 2D
    dsl4jsb_Real2D_onChunk :: sw_srf_net
    dsl4jsb_Real2D_onChunk :: swvis_srf_net
    dsl4jsb_Real2D_onChunk :: swnir_srf_net
    dsl4jsb_Real2D_onChunk :: rad_srf_net
    dsl4jsb_Real2D_onChunk :: net_radiation
    dsl4jsb_Real2D_onChunk :: alb_background
    dsl4jsb_Real2D_onChunk :: alb_vis
    dsl4jsb_Real2D_onChunk :: alb_nir
    dsl4jsb_Real2D_onChunk :: alb_vis_soil
    dsl4jsb_Real2D_onChunk :: alb_nir_soil
    dsl4jsb_Real2D_onChunk :: alb_vis_can
    dsl4jsb_Real2D_onChunk :: alb_nir_can
    dsl4jsb_Real2D_onChunk :: alb_vis_snow
    dsl4jsb_Real2D_onChunk :: alb_nir_snow
    dsl4jsb_Real2D_onChunk :: alb_vis_lnd
    dsl4jsb_Real2D_onChunk :: alb_nir_lnd
    dsl4jsb_Real2D_onChunk :: arad_vis_soil
    dsl4jsb_Real2D_onChunk :: arad_nir_soil
    dsl4jsb_Real2D_onChunk :: arad_vis_can
    dsl4jsb_Real2D_onChunk :: arad_nir_can
    dsl4jsb_Real2D_onChunk :: appfd
    dsl4jsb_Real2D_onChunk :: rfr_ratio_boc 
    dsl4jsb_Real2D_onChunk :: albedo
    dsl4jsb_Real2D_onChunk :: snow_age
    ! RAD_ 3D
    dsl4jsb_Real3D_onChunk :: ppfd_sunlit_cl
    dsl4jsb_Real3D_onChunk :: ppfd_shaded_cl
    dsl4jsb_Real3D_onChunk :: arad_sunlit_vis_cl
    dsl4jsb_Real3D_onChunk :: arad_shaded_vis_cl
    dsl4jsb_Real3D_onChunk :: fleaf_sunlit_vis_cl
    dsl4jsb_Real3D_onChunk :: arad_sunlit_nir_cl
    dsl4jsb_Real3D_onChunk :: arad_shaded_nir_cl
    dsl4jsb_Real3D_onChunk :: fleaf_sunlit_nir_cl
    ! HYDRO_
    dsl4jsb_Real2D_onChunk :: fract_snow_soil
    dsl4jsb_Real2D_onChunk :: fract_snow_can
    ! SEB_
    dsl4jsb_Real2D_onChunk :: t_unfilt
    ! VEG_
    dsl4jsb_Real3D_onChunk :: lai_cl
    dsl4jsb_Real3D_onChunk :: leaf_nitrogen_cl 
    dsl4jsb_Real3D_onChunk :: fn_chl_cl 
    dsl4jsb_Real3D_onChunk :: cumm_lai_cl
    dsl4jsb_Real3D_onChunk :: fleaf_sunlit_cl
    ! PHENO_
    dsl4jsb_Real2D_onChunk :: lai
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(RAD_)) RETURN   ! ask for all processes considered here ?
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! only at PFT tiles (i.e. leafs of vegetation tile) because this routine is tailored towards vegetation only
    !  in the "usecase_pfts" RAD_ is defined to run at leafs, which includes also other tiles than PFT
    !  hence, the jsbach4 task "surface_radiation" needs to run at leafs other than PFT if "use_quincy = T"
    IF (tile%lcts(1)%lib_id == 0) THEN
      IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Run only at PFT tiles. '//TRIM(tile%name)//' is not a PFT. ')
      RETURN
    ENDIF

    ! ---------------------------      
    ! 0.5 Get Memory
    model         => Get_model(tile%owner_model_id)
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    vgrid_canopy  => Get_vgrid('canopy_layer_q_')
    ncanopy       =  vgrid_canopy%n_levels
    ! Get process config
    dsl4jsb_Get_config(RAD_)
    ! Get process memories
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SEB_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(A2L_, swpar_srf_down)           ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, fract_par_diffuse)        ! in   ! jsb4 fract_par_diffuse == quincy swpar_fdiffuse
    dsl4jsb_Get_var2D_onChunk(A2L_, cos_zenith_angle)         ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, swnir_srf_down)           ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, swvis_srf_down)           ! in
    ! ---------------------------
    ! 2D
    dsl4jsb_Get_var2D_onChunk(RAD_, sw_srf_net)               ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, swvis_srf_net)            ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, swnir_srf_net)            ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, rad_srf_net)              ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, net_radiation)            ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, alb_background)           ! in
    dsl4jsb_Get_var2D_onChunk(RAD_, alb_vis)                  ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, alb_nir)                  ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, alb_vis_soil)             ! in      ! get from bc_land_phy.nc file; do not calc
    dsl4jsb_Get_var2D_onChunk(RAD_, alb_nir_soil)             ! in      ! get from bc_land_phy.nc file; do not calc
    dsl4jsb_Get_var2D_onChunk(RAD_, alb_vis_can)              ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, alb_nir_can)              ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, alb_vis_snow)             ! inout
    dsl4jsb_Get_var2D_onChunk(RAD_, alb_nir_snow)             ! inout
    dsl4jsb_Get_var2D_onChunk(RAD_, alb_vis_lnd)              ! out     ! set alb_vis_lnd = alb_vis (for consistency with jsb4 update_albedo Task)
    dsl4jsb_Get_var2D_onChunk(RAD_, alb_nir_lnd)              ! out     ! set alb_nir_lnd = alb_nir (for consistency with jsb4 update_albedo Task)
    dsl4jsb_Get_var2D_onChunk(RAD_, arad_vis_soil)            ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, arad_nir_soil)            ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, arad_vis_can)             ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, arad_nir_can)             ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, appfd)                    ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, rfr_ratio_boc)            ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, albedo)                   ! out
    dsl4jsb_Get_var2D_onChunk(RAD_, snow_age)                 ! in
    ! 3D
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_cl)           ! out
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_cl)           ! out
    dsl4jsb_Get_var3D_onChunk(RAD_, arad_sunlit_vis_cl)       ! out
    dsl4jsb_Get_var3D_onChunk(RAD_, arad_shaded_vis_cl)       ! out
    dsl4jsb_Get_var3D_onChunk(RAD_, fleaf_sunlit_vis_cl)      ! out
    dsl4jsb_Get_var3D_onChunk(RAD_, arad_sunlit_nir_cl)       ! out
    dsl4jsb_Get_var3D_onChunk(RAD_, arad_shaded_nir_cl)       ! out
    dsl4jsb_Get_var3D_onChunk(RAD_, fleaf_sunlit_nir_cl)      ! out
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow_soil)        ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow_can)         ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SEB_,   t_unfilt)               ! in
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(VEG_, lai_cl)                ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, leaf_nitrogen_cl)      ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_chl_cl)             ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, cumm_lai_cl)           ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_cl)       ! out
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(PHENO_, lai)                    ! in

      

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 0.9  Set day-or-night-switch
    !! 
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR
    DO ic=1,nc
      l_day(ic) = Has_minimal_radiation(swvis_srf_down(ic), swnir_srf_down(ic))
    END DO
    !$ACC END PARALLEL LOOP

    !> 1.0 snow albedo 
    !! 
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR
    DO ic=1,nc
      !! (original jsbach4 routine)
      CALL calc_snow_albedo(l_day(ic),                                 & ! in
                            fract_snow_soil(ic),                       & ! in
                            snow_age(ic),                              & ! in
                            cos_zenith_angle(ic),                      & ! in
                            dsl4jsb_Config(RAD_)%albedo_age_weight,    & ! in
                            t_unfilt(ic),                              & ! in
                            alb_vis_snow(ic),                          & ! inout
                            alb_nir_snow(ic)                           & ! inout
                            )
    END DO
    !$ACC END PARALLEL LOOP


    !> 2.0 soil albedo 
    !! 
    ! Do nothing: alb_vis_soil and alb_nir_soil are already in memory initialized from the bc_land_phy.nc file
    ! and represent the soil albedo including carbon in the soil as well as litter on the soil, but without
    ! the canopy above.
    IF (debug_on() .AND. iblk == 1) THEN
      CALL message(TRIM(routine), 'Soil albedo scheme not used for veg tiles but taken from bc_land_phy.nc file')
    END IF


    !> 3.0 calc "snow & soil" albedo
    !! 
    !! these two variables are local to this routine and only needed for the quincy calc_sw_srf_net() routine
    alb_vis_snow_soil(:) = fract_snow_soil(:) * alb_vis_snow(:) + (1.0_wp - fract_snow_soil(:)) * alb_vis_soil(:)
    alb_nir_snow_soil(:) = fract_snow_soil(:) * alb_vis_snow(:) + (1.0_wp - fract_snow_soil(:)) * alb_vis_soil(:)

    
    !> 4.0 calculate VIS radiation budget (absorbed radiation on canopy layers and soil + albedo)
    !! 
    CALL calc_sw_srf_net(nc, ncanopy, &
                         lctlib%sigma_vis, lctlib%sigma_nir, lctlib%omega_clumping, lctlib%crown_shape_factor, &
                         "vis", l_day(:), &
                         swvis_srf_down(:), fract_par_diffuse(:), cos_zenith_angle(:), &  ! jsb4 fract_par_diffuse == quincy swpar_fdiffuse
                         lai(:), lai_cl(:,:), cumm_lai_cl(:,:), & 
                         alb_vis_snow_soil(:), alb_vis_snow(:), fract_snow_can, &
                         arad_sunlit_vis_cl(:,:), arad_shaded_vis_cl(:,:), fleaf_sunlit_vis_cl(:,:), &   ! inout
                         arad_vis_can(:), arad_vis_soil(:), &
                         alb_vis_can(:), alb_vis(:), swvis_srf_net(:)) 
    !> 4.1 for consistency with jsb4
    alb_vis_lnd(:) = alb_vis(:)
    

    !> 5.0 calculate NIR radiation budget (absorbed radiation on canopy layers and soil + albedo)
    !! 
    !! Following Weiss & Norman, Agri. For. Met. 1985, the fraction of VIS and NIR 
    !! diffuse radiation are very similar 
    CALL calc_sw_srf_net(nc, ncanopy, &
                         lctlib%sigma_vis, lctlib%sigma_nir, lctlib%omega_clumping, lctlib%crown_shape_factor, &
                         "nir", l_day(:), &
                         swnir_srf_down(:), fract_par_diffuse(:), cos_zenith_angle(:), &  ! jsb4 fract_par_diffuse == quincy swpar_fdiffuse
                         lai(:), lai_cl(:,:), cumm_lai_cl(:,:), &
                         alb_nir_snow_soil(:), alb_nir_snow(:), fract_snow_can, &
                         arad_sunlit_nir_cl(:,:), arad_shaded_nir_cl(:,:), fleaf_sunlit_nir_cl(:,:), &   ! inout
                         arad_nir_can(:), arad_nir_soil(:), &
                         alb_nir_can(:), alb_nir(:), swnir_srf_net(:)) 
    !> 5.1 for consistency with jsb4
    alb_nir_lnd(:) = alb_nir(:)


    !> 6.0 quincy helper routine
    !!  code needed in addition to calling the calc_sw_srf_net() for vis/nir
    !!  copy of the original code in update_radiation() sections 4.0 to 7.0 of the quincy model
    CALL calc_albedo_radiation_quincy_hlp(nc, &                       ! in
                                          ncanopy, &
                                          l_day(:), &
                                          swvis_srf_down(:), &
                                          swnir_srf_down(:), &
                                          swvis_srf_net(:), &
                                          swnir_srf_net(:), &
                                          alb_vis(:), & 
                                          alb_nir(:), & 
                                          fleaf_sunlit_vis_cl(:,:), &
                                          swpar_srf_down(:), &
                                          arad_vis_can(:), &
                                          arad_sunlit_vis_cl(:,:), &
                                          arad_shaded_vis_cl(:,:), &
                                          leaf_nitrogen_cl(:,:), &
                                          fn_chl_cl(:,:), &
                                          lai_cl(:,:), &
                                          fract_snow_can(:), &        ! in
                                          albedo(:), &                ! inout
                                          sw_srf_net(:), &            ! out
                                          rad_srf_net(:), &           ! out
                                          net_radiation(:), &         ! out
                                          fleaf_sunlit_cl(:,:), &     ! out
                                          appfd(:), &                 ! out
                                          ppfd_sunlit_cl(:,:), &      ! out
                                          ppfd_shaded_cl(:,:), &      ! out
                                          rfr_ratio_boc(:)   &        ! out   
                                          )
   
  END SUBROUTINE update_sw_radiation_q_


  !-----------------------------------------------------------------------------------------------------
  !> Implementation of "aggregate": update_radiation task
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE aggregate_sw_radiation_q_(tile, options)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),        POINTER         :: model
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice, nc
    !X if necessary: REAL(wp) :: dtime, steplen
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggregate_sw_radiation_q_'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(VEG_)
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process variables (Set pointers to variables in memory)
    model => Get_model(tile%owner_model_id)
    ! Get process config
    ! ...
    ! Get process memories
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(VEG_)


    ! ------------------------------------------------------------------------------------------------------------
    ! Go aggregate
    ! ------------------------------------------------------------------------------------------------------------

    ! only at tiles that contain vegetation
    IF (tile%is_vegetation) THEN
      weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

      dsl4jsb_Aggregate_onChunk(RAD_, sw_srf_net              , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, swvis_srf_net           , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, swnir_srf_net           , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, rad_srf_net             , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, net_radiation           , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis                 , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir                 , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_can             , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_can             , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_snow            , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_snow            , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_lnd             , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_lnd             , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, arad_vis_soil           , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, arad_nir_soil           , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, arad_vis_can            , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, arad_nir_can            , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, appfd                   , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, rfr_ratio_boc           , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, albedo                  , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, ppfd_sunlit_cl          , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, ppfd_shaded_cl          , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, arad_sunlit_vis_cl      , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, arad_shaded_vis_cl      , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, fleaf_sunlit_vis_cl     , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, arad_sunlit_nir_cl      , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, arad_shaded_nir_cl      , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, fleaf_sunlit_nir_cl     , weighted_by_fract)
      ! ---------------------------
      dsl4jsb_Aggregate_onChunk(VEG_, fleaf_sunlit_cl      , weighted_by_fract)
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_sw_radiation_q_


#endif
END MODULE mo_rad_interface
