! --------------------------------------------------------------------
!> Example plugin for the ICON Community Interface (ComIn)
!  using the yayt library to gather one ICON variable (temperature) on
!  task zero to calculate the global average near surface temperature.
!
!  @authors 08/2021 :: ICON Community Interface  <icon@dwd.de>
!
!  Note that in order to demonstrate ComIn's language interoperability,
!  a similary plugin has been implemented in C, see the subdirectory
!  "yaxt_c".
! --------------------------------------------------------------------
MODULE yaxt_fortran_plugin

  USE iso_c_binding,           ONLY : C_INT
  USE comin_plugin_interface,  ONLY : comin_callback_register    &
                                    , comin_var_get              &
                                    , t_comin_var_descriptor     &
                                    , t_comin_var_ptr            &
                                    , comin_descrdata_get_domain &
                                    , t_comin_descrdata_domain   &
                                    , comin_descrdata_get_global &
                                    , t_comin_descrdata_global   &
                                    , t_comin_setup_version_info     &
                                    , comin_setup_get_version        &
                                    , EP_SECONDARY_CONSTRUCTOR       &
                                    , EP_DESTRUCTOR                  &
                                    , EP_ATM_WRITE_OUTPUT_BEFORE     &
                                    , COMIN_FLAG_READ                &
                                    , comin_parallel_get_host_mpi_rank &
                                    , comin_current_get_domain_id      &
                                    , comin_parallel_get_plugin_mpi_comm &
                                    , COMIN_DOMAIN_OUTSIDE_LOOP        &
                                    , t_comin_plugin_info              &
                                    , comin_current_get_plugin_info    &
                                    , comin_plugin_finish              &
                                    , comin_metadata_get, comin_var_to_3d
  USE yaxt,                     ONLY: xt_redist, xi => xt_int_kind
  
  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: pluginname = "yaxt_fortran_plugin"

  !> working precision (will be compared to ComIn's and ICON's)
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)
  TYPE(t_comin_setup_version_info) :: version

  TYPE(t_comin_var_ptr),  POINTER :: temp ! temperature
  INTEGER                         :: rank
  CHARACTER(LEN=:), POINTER       :: units
  
  !> access descriptive data structures
  TYPE(t_comin_descrdata_domain),     POINTER   :: p_patch
  TYPE(t_comin_descrdata_global),     POINTER   :: p_global

  !> yaxt related variables
  TYPE(xt_redist) :: yaxt_redist
  
CONTAINS

  ! --------------------------------------------------------------------
  ! ComIn primary constructor.
  ! --------------------------------------------------------------------
  SUBROUTINE comin_main()  BIND(C)

    USE yaxt, ONLY: xt_initialize, xt_idxlist, xt_idxlist_delete &
                  , xt_xmap, xt_xmap_delete, xt_idxstripes_new, xt_idxvec_new &
                  , xt_stripe, xt_idxempty_new, xt_xmap_all2all_new &
                  , xt_redist_p2p_new, xt_initialized
    USE mpi,  ONLY: MPI_DOUBLE_PRECISION

    CHARACTER(LEN=*), PARAMETER   :: substr = 'comin_main (yaxt_fortran_plugin)'
    TYPE(t_comin_plugin_info)     :: this_plugin
    INTEGER                       :: ierr
    INTEGER                       :: p_all_comm ! communicator of all ICON tasks
    TYPE(xt_idxlist)              :: src_idxlist, tgt_idxlist
    TYPE(xt_xmap)                 :: xmap
    INTEGER(kind=xi), DIMENSION(:), ALLOCATABLE :: idxvec, idxmap
    INTEGER                       :: k
    
    ierr = 0

    !> get the rank of the current process and say hello to the world
    rank = comin_parallel_get_host_mpi_rank()
    IF (rank == 0)  WRITE (0,*) "     setup: ", pluginname

    !> check, if the ComIn library version is compatible
    version = comin_setup_get_version()
    IF (version%version_no_major > 1)  THEN
       CALL comin_plugin_finish(substr, "incompatible ComIn library version!")
    END IF

    !> check plugin id
    CALL comin_current_get_plugin_info(this_plugin, ierr)
    IF (ierr /= 0) &
         CALL comin_plugin_finish(substr, "plugin info could not be retrieved")
    IF (rank == 0)  WRITE (0,'(a,a,a,i4)') "     plugin " &
         , TRIM(this_plugin%name), " has id: ", this_plugin%id

    !> add requests for additional ICON variables
    !  not applicable for this example
    
    !> register callbacks
    CALL comin_callback_register(EP_SECONDARY_CONSTRUCTOR &
         , yaxt_fortran_constructor, ierr)
    IF (ierr /= 0) &
         CALL comin_plugin_finish(substr &
         , "SECONDARY_CONSTRUCTOR could not be registered")
    
    CALL comin_callback_register(EP_ATM_WRITE_OUTPUT_BEFORE &
         , yaxt_fortran_gather,     ierr)
    IF (ierr /= 0) &
         CALL comin_plugin_finish(substr &
         , "ATM_WRTITE_OUTPUT_BEFORE could not be registered")

    CALL comin_callback_register(EP_DESTRUCTOR &
         , yaxt_fortran_destructor,  ierr)
    IF (ierr /= 0) &
         CALL comin_plugin_finish(substr &
         , "DESTRUCTOR could not be registered")
    
    !> get descriptive data structures
    p_patch    => comin_descrdata_get_domain(1)
    p_global   => comin_descrdata_get_global()

    !> setup yaxt
    p_all_comm = comin_parallel_get_plugin_mpi_comm()
    IF (.NOT. xt_initialized()) THEN
       WRITE (0,*) "Initialize yaxt..."
       CALL xt_initialize(p_all_comm)
    ENDIF

    !> construct yaxt variables ...
    idxmap = RESHAPE(p_patch%cells%decomp_domain &
         , (/ SIZE(p_patch%cells%decomp_domain) /))
    idxvec = INT(PACK( [(k,k=1,p_patch%cells%ncells)], idxmap == 0 ), xi)
    idxvec = p_patch%cells%glb_index(idxvec)
    src_idxlist = xt_idxvec_new(idxvec)
    
    IF (rank == 0) THEN
       tgt_idxlist = xt_idxstripes_new( &
            (/ xt_stripe(1, 1, p_patch%cells%ncells_global) /))
    ELSE
       ! ... empty on all other pe
       tgt_idxlist = xt_idxempty_new()
    ENDIF

    ! ... create exchange map ...
    xmap = xt_xmap_all2all_new(src_idxlist, tgt_idxlist, p_all_comm)

    ! ... create redistribution instance for DP ...
    yaxt_redist = xt_redist_p2p_new(xmap, MPI_DOUBLE_PRECISION)
    
    ! ... clean up
    CALL xt_xmap_delete(xmap)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(tgt_idxlist)
    DEALLOCATE(idxvec)
    DEALLOCATE(idxmap)
    
  END SUBROUTINE comin_main

  ! --------------------------------------------------------------------
  ! ComIn secondary constructor.
  ! --------------------------------------------------------------------
  SUBROUTINE yaxt_fortran_constructor()  BIND(C)

    USE yaxt, ONLY: xt_idxlist, xt_idxlist_delete &
                  , xt_xmap, xt_xmap_delete &
                  , xt_redist

    CHARACTER(LEN=*), PARAMETER :: substr = 'yaxt_fortran_constructor (yaxt_fortran_plugin)'
    TYPE(t_comin_var_descriptor)         :: var_desc
    INTEGER                              :: ierr

    IF (rank == 0)  WRITE (0,*) substr, ": secondary constructor"

    IF (rank == 0)  WRITE (0,*) substr, ': ', pluginname, "request temperature"
    var_desc%name = 'temp'
    var_desc%id = 1
    CALL comin_var_get([EP_ATM_WRITE_OUTPUT_BEFORE], &
                       var_desc, COMIN_FLAG_READ, temp)
    IF (.NOT. ASSOCIATED(temp)) &
         CALL comin_plugin_finish(substr, "temperature variable request failed")

    CALL comin_metadata_get(var_desc, 'units', units, ierr)
    IF (ierr /= 0) &
         CALL comin_plugin_finish(substr, "units not found")
    
  END SUBROUTINE yaxt_fortran_constructor

  ! --------------------------------------------------------------------
  ! ComIn callback function.
  ! --------------------------------------------------------------------
  SUBROUTINE yaxt_fortran_gather()  BIND(C)

    USE yaxt,          ONLY: xt_redist_s_exchange1
    USE iso_c_binding, ONLY: C_LOC

    CHARACTER(LEN=*), PARAMETER :: substr = 'yaxt_fortran_gather (yaxt_fortran_plugin)'
    TYPE(t_comin_plugin_info)                :: this_plugin
    INTEGER                                  :: ierr, domain_id
    REAL(kind=wp), DIMENSION(:,:),   POINTER :: src
    REAL(kind=wp), DIMENSION(:,:,:), POINTER :: src3d
    REAL(kind=wp), DIMENSION(:),     POINTER :: tgt, area
    
    IF (rank == 0)  WRITE (0,*) substr, ": callback before output"

    !> check plugin id
    CALL comin_current_get_plugin_info(this_plugin, ierr)
    IF (ierr /= 0) &
         CALL comin_plugin_finish(substr, "plugin info could not be retrieved")
    
    IF (rank == 0)  WRITE (0,'(a,a,a,a,i4)') substr,": ", this_plugin%name, " has id: ", this_plugin%id

    domain_id = comin_current_get_domain_id()
    IF (ierr /= 0) &
         CALL comin_plugin_finish(substr, "domain ID could not be retrieved")
    IF (domain_id == COMIN_DOMAIN_OUTSIDE_LOOP) THEN
       IF (rank == 0)  WRITE(0,'(a,a)')  &
            substr, ": currently not in domain loop"
       !       RETURN
    ELSE
       IF (rank == 0)  WRITE(0,'(a,a,i0)')  &
            substr,": currently on domain ", domain_id
    END IF

    !IF (domain_id /= 1) THEN
    !   IF (rank == 0)  WRITE(0,'(a,a)')  &
    !        substr,": exiting, since domain /= 1"
    !   RETURN
    !END IF
    
    !> reset pointers
    NULLIFY(src, src3d, tgt, area)

    !> 5D -> 3D 
    src3d => comin_var_to_3d(temp)    

    !> extract near surface temperature
    ALLOCATE(src(p_global%nproma,p_patch%cells%nblks))
    src(:,:) = src3d(:,p_patch%nlev,:)

    !> allocate local space to gather global information
    ALLOCATE(tgt(p_patch%cells%ncells_global))
    ALLOCATE(area(p_patch%cells%ncells_global))

    !> gather information in rank zero
    CALL xt_redist_s_exchange1(yaxt_redist, C_LOC(src), C_LOC(tgt))
    CALL xt_redist_s_exchange1(yaxt_redist, C_LOC(p_patch%cells%area) &
         , C_LOC(area) )

    IF (rank == 0) &
         WRITE(0,*) substr, ': global average temperature is ' &
         , SUM(tgt*area)/SUM(area), TRIM(units)

    !> clean up memory
    DEALLOCATE(src)
    DEALLOCATE(tgt)
    DEALLOCATE(area)
    
  END SUBROUTINE yaxt_fortran_gather

  ! --------------------------------------------------------------------
  ! ComIn callback function.
  ! --------------------------------------------------------------------
  SUBROUTINE yaxt_fortran_destructor() BIND(C)

    USE yaxt, ONLY: xt_finalize, xt_redist_delete

    CHARACTER(LEN=*), PARAMETER :: substr = 'yaxt_fortran_destructor (yaxt_fortran_plugin)'
    
    IF (rank == 0)  WRITE (0,*) substr,": destructor."

    NULLIFY(temp, p_patch, p_global)

    !> free yaxt related memory
    CALL xt_redist_delete(yaxt_redist)

    !> finalize yaxt
    CALL xt_finalize()
    
  END SUBROUTINE yaxt_fortran_destructor

END MODULE yaxt_fortran_plugin
