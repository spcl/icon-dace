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

! Provides the master control methods and paramaters

MODULE mo_master_init

  USE mo_exception,     ONLY: message, finish
  USE mo_mpi,           ONLY: set_process_mpi_name, get_my_global_mpi_id, split_global_mpi_communicator
  USE mo_io_units,      ONLY: filename_max
  USE mo_master_config
  USE mo_master_nml,    ONLY: read_master_namelist

  IMPLICIT NONE
  PUBLIC :: init_master_control
   

CONTAINS

  !------------------------------------------------------------------------
  !>
  !!  Initialization of the master control variables
  !!
  INTEGER FUNCTION init_master_control(namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename

    INTEGER :: master_namelist_status
    INTEGER :: model_no, model_rank, start_rank, end_rank, inc_rank, group_size, current_rank, j

    CHARACTER(LEN=*), PARAMETER :: method_name = "master_control"
    !-----------------------------------------------------------------------

    CALL message(method_name,'start model initialization.')

    master_namelist_status = read_master_namelist(TRIM(namelist_filename))

    !------------------------------------------------------------
    ! some checks

    IF (master_namelist_status == -1) THEN
      CALL finish(method_name,'model identity (atm/oce) can no longer be derived from namelist run_nml!')
    ENDIF

    IF (noOfModels() < 1) THEN
      CALL finish(method_name,'no of models < 1')
    ENDIF

    !------------------------------------------------------------

    master_namelist_filename = TRIM(namelist_filename)

    !------------------------------------------------------------
    ! find what is my process

    multiple_models = noOfModels() > 1    

    IF ( multiple_models ) THEN
      
      CALL set_my_component_null()

      COMPONENT_MODELS: DO model_no = 1, noOfModels()

        !         write(0,*) 'master_component_models:', model_no, trim(master_component_models(model_no)%model_name)
        start_rank   = master_component_models(model_no)%model_min_rank
        end_rank  = master_component_models(model_no)%model_max_rank
        inc_rank = master_component_models(model_no)%model_inc_rank
        group_size = master_component_models(model_no)%model_rank_group_size
        
        model_rank = start_rank
        
        DO WHILE(model_rank .LE. end_rank)
        
          DO j = 0, group_size-1, 1

            current_rank = model_rank + j
            IF (current_rank > end_rank) EXIT
            
            IF ( get_my_global_mpi_id() == current_rank ) THEN
              
              CALL set_my_component(model_no,                                      &
                  &                master_component_models(model_no)%model_name,  &
                  &                master_component_models(model_no)%model_type,  &
                  &                master_component_models(model_no)%model_namelist_filename)
              
            ENDIF
            
          ENDDO
          
          model_rank = model_rank + group_size-1 + inc_rank
        
        ENDDO ! WHILE(model_rank .LE. end_rank)
        
      ENDDO COMPONENT_MODELS

      CALL split_global_mpi_communicator ( my_model_no, noOfModels() )

    ELSE ! only one component    

      model_no = 1
      CALL set_my_component(model_no,                                      &
           &                master_component_models(model_no)%model_name,  &
           &                master_component_models(model_no)%model_type,  &
           &                master_component_models(model_no)%model_namelist_filename)

      CALL split_global_mpi_communicator ( model_no, 1 )

    ENDIF

    !------------------------------------------------------------
    
    CALL check_my_component()
    
    !------------------------------------------------------------
    
    init_master_control = 0
    
  END FUNCTION init_master_control
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  SUBROUTINE set_my_component(comp_no, comp_name, comp_id, comp_namelist)

    INTEGER, INTENT(in)          :: comp_no
    CHARACTER(len=*), INTENT(in) :: comp_name
    INTEGER, INTENT(in)          :: comp_id
    CHARACTER(len=*), INTENT(in) :: comp_namelist

    my_model_no          = comp_no
    my_process_model     = comp_id
    my_namelist_filename = TRIM(comp_namelist)
    my_model_name        = TRIM(comp_name)

    my_model_min_rank    = master_component_models(comp_no)%model_min_rank
    my_model_max_rank    = master_component_models(comp_no)%model_max_rank
    my_model_inc_rank    = master_component_models(comp_no)%model_inc_rank

    CALL check_my_component()
    
    CALL set_process_mpi_name(TRIM(my_model_name))

  END SUBROUTINE set_my_component
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  SUBROUTINE check_my_component()

    CHARACTER(len=*), PARAMETER :: method_name='mo_master_control:check_my_component'

    IF (my_model_no < 1) CALL finish(method_name, 'my_model_no < 1') 
    IF (my_namelist_filename == '') CALL finish(method_name, 'my_namelist_filename = NULL')
    IF (my_model_name == '') CALL finish(method_name, 'my_model_name = NULL')

    SELECT CASE (my_process_model)
      CASE (atmo_process)
      CASE (ocean_process)
      CASE (ps_radiation_process)
      CASE (hamocc_process)
      CASE (jsbach_process)
      CASE (testbed_process)
      CASE (icon_output_process)
      CASE (wave_process)
      CASE default
        CALL finish("check_my_component","my_process_model is unkown")
    END SELECT

  END SUBROUTINE check_my_component
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  SUBROUTINE set_my_component_null()

    my_model_no          = 0
    my_process_model     = -1
    my_namelist_filename = ''
    my_model_name        = ''
    my_model_min_rank    = -1
    my_model_max_rank    = -2
    my_model_inc_rank    = -1

  END SUBROUTINE set_my_component_null
  !------------------------------------------------------------------------
END MODULE mo_master_init
!--------------------------------------------------------------------------------------------------------

