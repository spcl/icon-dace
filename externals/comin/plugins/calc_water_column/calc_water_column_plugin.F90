! --------------------------------------------------------------------
!> Example plugin for the ICON Community Interface (ComIn)
!
! This example illustrates a simple diagnostic application, i.e.
! the calculations of the liquid water path(lwp), the ice water path (iwp)
! and the total water column (twc).
! For this 
! - the variables lwp, iwp and twc need to be added to the ICON variables
!   including the definition of metadata, as the units
! - the humidity tracers from ICON need to be accessed and it needs to be
!   checked, if qg exists or not, which depends on the microphysics scheme
!   chosen in ICON.
! - for the calculation itself, the ICON prognostic variable rho is
!   required and
! - the descriptive data (p_patch%cell%hhl) is used
!
! The calculation is performed for each active domain.
!
!  @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  Please see the file LICENSE in the root of the source tree for this code.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
! --------------------------------------------------------------------
MODULE calc_water_column_plugin
  
  USE iso_c_binding,           ONLY : C_INT
  USE comin_plugin_interface,  ONLY : comin_callback_register,                                          &
    &                                 comin_var_get,                                                    &
    &                                 t_comin_var_descriptor, t_comin_var_ptr,                          &
    &                                 comin_var_request_add,                                            &
    &                                 comin_descrdata_get_domain, t_comin_descrdata_domain,             &
    &                                 comin_descrdata_get_global, t_comin_descrdata_global,             &
    &                                 t_comin_setup_version_info, comin_setup_get_version,              &
    &                                 EP_SECONDARY_CONSTRUCTOR, EP_DESTRUCTOR, EP_ATM_PHYSICS_BEFORE,   &
    &                                 COMIN_FLAG_READ, COMIN_FLAG_WRITE, COMIN_ZAXIS_2D,                &
    &                                 comin_parallel_get_host_mpi_rank, comin_current_get_domain_id,    &
    &                                 t_comin_plugin_info, comin_current_get_plugin_info,               &
    &                                 comin_plugin_finish, comin_metadata_set,                          &
    &                                 comin_metadata_get, comin_var_to_3d

  IMPLICIT NONE
  PRIVATE
  
  CHARACTER(LEN=*), PARAMETER :: pluginname = "calc_water_column_plugin"

  !> working precision (will be compared to ComIn's and ICON's)
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)
  TYPE(t_comin_setup_version_info) :: version

  TYPE :: t_plugin_vars
     ! ICON data
     TYPE(t_comin_var_ptr),  POINTER :: rho  => NULL()
     TYPE(t_comin_var_ptr),  POINTER :: qv => NULL()
     TYPE(t_comin_var_ptr),  POINTER :: qc => NULL()
     TYPE(t_comin_var_ptr),  POINTER :: qi => NULL()
     TYPE(t_comin_var_ptr),  POINTER :: qr => NULL()
     TYPE(t_comin_var_ptr),  POINTER :: qs => NULL()
     TYPE(t_comin_var_ptr),  POINTER :: qg => NULL()
     ! variables added to ICON by this plugin
     TYPE(t_comin_var_ptr),  POINTER :: lwp  => NULL() ! liquid water path (qc+qr)
     TYPE(t_comin_var_ptr),  POINTER :: iwp  => NULL() ! ice water path    (qi+qs+qg)
     TYPE(t_comin_var_ptr),  POINTER :: twc  => NULL() ! (full) water path (qv+qc+qr+qi+qs+qg)
  END type t_plugin_vars
  ! pugin variables per patch
  TYPE(t_plugin_vars), DIMENSION(:), ALLOCATABLE :: pvpp
     
  !> access descriptive data structures
  TYPE(t_comin_descrdata_domain),     POINTER   :: p_patch
  !  access global setup information
  TYPE(t_comin_descrdata_global),     POINTER   :: p_global
  ! 
  LOGICAL :: lqg = .TRUE. ! qg is not  always available
  
  PUBLIC :: comin_main
  PUBLIC :: calc_water_column_constructor
  PUBLIC :: calc_water_column_diagfct
  PUBLIC :: calc_water_column_destructor
  
CONTAINS

  ! --------------------------------------------------------------------
  ! ComIn primary constructor.
  ! --------------------------------------------------------------------
  SUBROUTINE comin_main()  BIND(C)
    !
    IMPLICIT NONE
    !
    TYPE(t_comin_plugin_info)     :: this_plugin
    TYPE(t_comin_var_descriptor)  :: lwp_d, iwp_d, twc_d
    INTEGER                       :: status
    INTEGER                       :: jg
    CHARACTER(LEN=120)            :: text
    
    status = 0

    CALL message('- setup')

    version = comin_setup_get_version()
    IF (version%version_no_major > 1)  THEN
      CALL comin_plugin_finish('comin_main ('//pluginname//')', "incompatible version!")
    END IF

    ! get descriptive data structures
    p_global => comin_descrdata_get_global()

    WRITE (text,'(a,i4)') '- number of domains: ', p_global%n_dom
    CALL message(text)
    
    !> check plugin id
    CALL comin_current_get_plugin_info(this_plugin, status)
    IF (status /= 0)  RETURN
    WRITE (text,'(a,a,a,i4)') "- plugin ", this_plugin%name, " has id: ", this_plugin%id
    CALL message(text)

    !> add requests for additional ICON variables

    ! request host model to add local variables

    DO jg = 1, p_global%n_dom
       ! liquid water path (lwp) for first domain
       lwp_d = t_comin_var_descriptor( id = jg, name = "lwp" )
       CALL comin_var_request_add_wrapper(status, lwp_d, lmode_exclusive=.FALSE. &
            , zaxis_id = COMIN_ZAXIS_2d, tracer =.false., restart=.false., units='kg/m2')
       IF (status /= 0) CALL comin_plugin_finish('comin_main ('//pluginname//')' &
            , 'variable request failed for lwp')
       
       ! ice water path (iwp) for first domain
       iwp_d = t_comin_var_descriptor( id = jg, name = "iwp" )
       CALL comin_var_request_add_wrapper(status, iwp_d, lmode_exclusive=.FALSE. &
            , zaxis_id = COMIN_ZAXIS_2d, tracer =.false., restart=.false., units='kg/m2')
       IF (status /= 0)  CALL comin_plugin_finish('comin_main ('//pluginname//')' &
            , 'variable request failed for iwp')
       
       ! total water column (twc) for first domain
       twc_d = t_comin_var_descriptor( id = jg, name = "twc" )
       CALL comin_var_request_add_wrapper(status, twc_d, lmode_exclusive=.FALSE. &
            , zaxis_id = COMIN_ZAXIS_2d, tracer =.false., restart=.false., units='kg/m2')
       IF (status /= 0)  CALL comin_plugin_finish('comin_main ('//pluginname//')' &
            , 'variable request failed for twc')
    END DO
       
    ! register callbacks
    CALL comin_callback_register(EP_ATM_PHYSICS_BEFORE,      calc_water_column_diagfct,     status)
    CALL comin_callback_register(EP_SECONDARY_CONSTRUCTOR,   calc_water_column_constructor, status)
    CALL comin_callback_register(EP_DESTRUCTOR,              calc_water_column_destructor,  status)

  END SUBROUTINE comin_main


  ! --------------------------------------------------------------------
  ! ComIn secondary constructor.
  ! --------------------------------------------------------------------
  SUBROUTINE calc_water_column_constructor()  BIND(C)
    
    IMPLICIT NONE

    TYPE(t_comin_var_descriptor)  :: var_desc
    INTEGER  :: jg

    ALLOCATE(pvpp(p_global%n_dom))

    domain_loop: DO jg = 1, p_global%n_dom
       CALL message(' - get required meteorological data ', jg)

       ! 1 - GET ICON VARIABLES:
       ! 1.1 - A RHO
       var_desc%name = 'rho'
       var_desc%id = jg
       CALL comin_var_get([EP_ATM_PHYSICS_BEFORE], &
            &                 var_desc, COMIN_FLAG_READ, pvpp(jg)%rho)
       CALL check_variable(pvpp(jg)%rho, 'rho')

       CALL message(' - get humidity tracer ',jg)
       ! 1.2  qv
       var_desc%name = 'qv'
       var_desc%id = jg
       CALL comin_var_get([EP_ATM_PHYSICS_BEFORE], &
            &                 var_desc, COMIN_FLAG_READ, pvpp(jg)%qv)
       CALL check_variable(pvpp(jg)%qv, 'qv')
       ! 1.3  qc
       var_desc%name = 'qc'
       var_desc%id = jg
       CALL comin_var_get([EP_ATM_PHYSICS_BEFORE], &
            &                 var_desc, COMIN_FLAG_READ, pvpp(jg)%qc)
       CALL check_variable(pvpp(jg)%qc,'qc')
       ! 1.4  qi
       var_desc%name = 'qi'
       var_desc%id = jg
       CALL comin_var_get([EP_ATM_PHYSICS_BEFORE], &
            &                 var_desc, COMIN_FLAG_READ, pvpp(jg)%qi)
       CALL check_variable(pvpp(jg)%qi, 'qi')
       ! 1.5  qr   
       var_desc%name = 'qr'
       var_desc%id = jg
       CALL comin_var_get([EP_ATM_PHYSICS_BEFORE], &
            &                 var_desc, COMIN_FLAG_READ, pvpp(jg)%qr)
       CALL check_variable(pvpp(jg)%qr, 'qr')
       ! 1.6  qs
       var_desc%name = 'qs'
       var_desc%id = jg
       CALL comin_var_get([EP_ATM_PHYSICS_BEFORE], &
            &                 var_desc, COMIN_FLAG_READ, pvpp(jg)%qs)
       CALL check_variable(pvpp(jg)%qs, 'qs')
       ! 1.7  qg
       var_desc%name = 'qg'
       var_desc%id = jg
       CALL comin_var_get([EP_ATM_PHYSICS_BEFORE], &
            &                 var_desc, COMIN_FLAG_READ, pvpp(jg)%qg)
       CALL check_variable(pvpp(jg)%qg, 'qg', lstop=.FALSE., lasso=lqg)
 
       CALL message(' - get plugin variables - requested to be created by ICON'&
            , jg)

       ! 2 - GET OWN VARIABLE requested in comin_main
       CALL comin_var_get([EP_ATM_PHYSICS_BEFORE], &
            &                t_comin_var_descriptor(name='lwp', id=jg),&
            &                COMIN_FLAG_WRITE, pvpp(jg)%lwp)
       CALL check_variable(pvpp(jg)%lwp, 'lwp',dim=2)
   
       CALL comin_var_get([EP_ATM_PHYSICS_BEFORE], &
            &                t_comin_var_descriptor(name='iwp', id=jg), &
            &                COMIN_FLAG_WRITE, pvpp(jg)%iwp)
       CALL check_variable(pvpp(jg)%lwp, 'iwp',dim=2)

       CALL comin_var_get([EP_ATM_PHYSICS_BEFORE], &
            &                t_comin_var_descriptor(name='twc', id=jg), &
            &                COMIN_FLAG_WRITE, pvpp(jg)%twc)
       CALL check_variable(pvpp(jg)%lwp, 'twc',dim=2)

    END DO domain_loop
    
  END SUBROUTINE calc_water_column_constructor


  ! --------------------------------------------------------------------
  ! ComIn callback function.
  ! --------------------------------------------------------------------
  SUBROUTINE calc_water_column_diagfct()  BIND(C)

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER     :: substr='calc_water_column_diagfct'
    TYPE(t_comin_var_descriptor)    :: lwp_d, iwp_d, twc_d
    INTEGER                         :: status
    INTEGER                         :: jk, jg
    
    REAL(WP), POINTER, DIMENSION(:,:,:) :: rho_3d  => NULL()

    REAL(WP), POINTER, DIMENSION(:,:,:) :: lwp_3d => NULL()
    REAL(WP), POINTER, DIMENSION(:,:,:) :: iwp_3d => NULL()
    REAL(WP), POINTER, DIMENSION(:,:,:) :: twc_3d => NULL()

    REAL(WP), POINTER, DIMENSION(:,:,:) :: qv_3d => NULL()
    REAL(WP), POINTER, DIMENSION(:,:,:) :: qc_3d => NULL()
    REAL(WP), POINTER, DIMENSION(:,:,:) :: qi_3d => NULL()
    REAL(WP), POINTER, DIMENSION(:,:,:) :: qr_3d => NULL()
    REAL(WP), POINTER, DIMENSION(:,:,:) :: qs_3d => NULL()
    REAL(WP), POINTER, DIMENSION(:,:,:) :: qg_3d => NULL()
    ! conversion factor kg / kg => kg / m2
    REAL(WP), POINTER, DIMENSION(:,:,:) :: conv => NULL() 

    CHARACTER(LEN=:), POINTER   :: units => NULL()
    CHARACTER(LEN=200)          :: text = ''
    
    ! get current domain:
    jg = comin_current_get_domain_id()
    ! get current patch description data
    p_patch  => comin_descrdata_get_domain(jg)
    
    CALL message('- calculate water columns before physics.', jg)

    qv_3d => comin_var_to_3d(pvpp(jg)%qv)
    qi_3d => comin_var_to_3d(pvpp(jg)%qi)
    qr_3d => comin_var_to_3d(pvpp(jg)%qr)
    qs_3d => comin_var_to_3d(pvpp(jg)%qs)
    qc_3d => comin_var_to_3d(pvpp(jg)%qc)
    IF (lqg) qg_3d => comin_var_to_3d(pvpp(jg)%qg)

    rho_3d => comin_var_to_3d(pvpp(jg)%rho)
 
    lwp_3d => comin_var_to_3d(pvpp(jg)%lwp)
    iwp_3d => comin_var_to_3d(pvpp(jg)%iwp)
    twc_3d => comin_var_to_3d(pvpp(jg)%twc)

    ! conversion factor 1/kg  to 1/m2
    ALLOCATE(conv(SIZE(qv_3d,1),SIZE(qv_3d,2),SIZE(qv_3d,3)))
    conv(:,:,:) = 0._wp
    DO jk=1,SIZE(qv_3d,2)
       ! convert 1/kg  to 1/m2
       conv(:,jk,:) = rho_3d(:,jk,:) * &
            (p_patch%cells%hhl(:,jk,:) - p_patch%cells%hhl(:,jk+1,:))
    END DO

    ! calculate liquid water part / ice_water_path and total water column 
    lwp_3d(:,:,:) = 0._wp 
    iwp_3d(:,:,:) = 0._wp 
    twc_3d(:,:,:) = 0._wp 
    DO jk=1,SIZE(qv_3d,2)
       lwp_3d(:,:,1) = lwp_3d(:,:,1)   &
            + (qr_3d(:,jk,:) + qc_3d(:,jk,:)) * conv(:,jk,:) 
       iwp_3d(:,:,1) = iwp_3d(:,:,1)   &
            + (qi_3d(:,jk,:)+qs_3d(:,jk,:)) * conv(:,jk,:)
    END DO
    IF (lqg) THEN
       DO jk=1,SIZE(qv_3d,2)
          iwp_3d(:,:,1) = iwp_3d(:,:,1)  + qg_3d(:,jk,:) * conv(:,jk,:) 
       END DO
    END IF
    twc_3d(:,:,1) = iwp_3d(:,:,1) + lwp_3d(:,:,1)
    DO jk=1,SIZE(qv_3d,2)
        twc_3d(:,:,1) = twc_3d(:,:,1)  + qv_3d(:,jk,:) * conv(:,jk,:) 
    END DO

    CALL message(' ')
    CALL message('- results of plugin calc_water_column: ',jg)

    lwp_d = t_comin_var_descriptor( id = jg, name = "lwp" )
    CALL comin_metadata_get(lwp_d, "units", units, status)
    IF (status /= 0) CALL comin_plugin_finish(substr, 'lwp unit')
    write (text,fmt='(A,A5,A,F8.4)') '- maximum liquid water path  (',TRIM(units),'): ', MAXVAL(lwp_3d(:,:,1))
    CALL message(text, idom=jg,lall=.TRUE.)
    iwp_d = t_comin_var_descriptor( id = jg, name = "iwp" )
    CALL comin_metadata_get(iwp_d, "units", units, status)
    IF (status /= 0) CALL comin_plugin_finish(substr, 'iwp unit')
    write (text,fmt='(A,A5,A,F8.4)') '- maximum ice water path     (',TRIM(units),'): ', MAXVAL(iwp_3d(:,:,1))
    CALL message(text, idom=jg,lall=.TRUE.)
    twc_d = t_comin_var_descriptor( id = jg, name = "twc" )
    CALL comin_metadata_get(twc_d, "units", units, status)
    IF (status /= 0) CALL comin_plugin_finish(substr, 'twc unit')
    write (text,fmt='(A,A5,A,F8.4)') '- maximum total water column (',TRIM(units),'): ', MAXVAL(twc_3d(:,:,1))
    CALL message(text, idom=jg,lall=.TRUE.)

    ! CLEANUP
    DEALLOCATE(conv); NULLIFY(conv)
    NULLIFY(lwp_3d, iwp_3d, twc_3d)
    NULLIFY(qv_3d,qc_3d,qi_3d,qr_3d,qs_3d,qg_3d, rho_3d)
    NULLIFY(p_patch)
    
  END SUBROUTINE calc_water_column_diagfct

  ! --------------------------------------------------------------------
  ! ComIn callback function.
  ! --------------------------------------------------------------------

  SUBROUTINE calc_water_column_destructor() BIND(C)

    IMPLICIT NONE

    INTEGER :: jg
    
    CALL message(' - destructor.')

    DO jg = 1, p_global%n_dom
       NULLIFY(pvpp(jg)%rho)
       NULLIFY(pvpp(jg)%qv)
       NULLIFY(pvpp(jg)%qc)
       NULLIFY(pvpp(jg)%qr)
       NULLIFY(pvpp(jg)%qi)
       NULLIFY(pvpp(jg)%qs)
       NULLIFY(pvpp(jg)%qg)
       NULLIFY(pvpp(jg)%twc)
       NULLIFY(pvpp(jg)%iwp)
       NULLIFY(pvpp(jg)%lwp)
    END DO
    DEALLOCATE(pvpp)
    
  END SUBROUTINE calc_water_column_destructor

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  ! PRIVATE ROUTINES
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  
  SUBROUTINE comin_var_request_add_wrapper(status, descriptor, lmode_exclusive, zaxis_id &
                                          , tracer, restart, units)

    IMPLICIT NONE

    ! I/O
    INTEGER,                                INTENT(OUT) :: status
    TYPE(t_comin_var_descriptor),           INTENT(IN)  :: descriptor
    LOGICAL,                      OPTIONAL, INTENT(IN)  :: lmode_exclusive
    INTEGER,                      OPTIONAL, INTENT(IN)  :: zaxis_id
    LOGICAL,                      OPTIONAL, INTENT(IN)  :: tracer
    LOGICAL,                      OPTIONAL, INTENT(IN)  :: restart
    CHARACTER(LEN=*),             OPTIONAL, INTENT(IN)  :: units
    ! LOCAL
    LOGICAL   :: lexclusive

    status = -1
    
    IF (PRESENT(lmode_exclusive)) THEN
       lexclusive = lmode_exclusive
    ELSE
       lexclusive = .FALSE.
    END IF

    CALL comin_var_request_add(descriptor, lexclusive, status)
    IF (status /= 0)  RETURN
    
    IF (PRESENT(zaxis_id)) THEN
       CALL comin_metadata_set(descriptor, "zaxis_id", zaxis_id, status)
       IF (status /= 0)  RETURN
    END IF
    IF (PRESENT(tracer)) THEN
       CALL comin_metadata_set(descriptor, "tracer", tracer, status)
       IF (status /= 0)  RETURN
    END IF
    IF (PRESENT(restart)) THEN
       CALL comin_metadata_set(descriptor, "restart", restart, status)
       IF (status /= 0)  RETURN
    END IF
    IF (PRESENT(units)) THEN
       CALL comin_metadata_set(descriptor, "units", TRIM(units), status)
       IF (status /= 0 ) RETURN
    END IF
    
    status = 0
    
  END SUBROUTINE comin_var_request_add_wrapper
  
  ! --------------------------------------------------------------------------------------------------

  SUBROUTINE check_variable(var, name, dim, lstop, lasso)

    IMPLICIT NONE

    ! I/O
    TYPE(t_comin_var_ptr),  POINTER :: var
    CHARACTER(LEN=*), INTENT(IN)    :: name
    INTEGER, INTENT(IN),  OPTIONAL  :: dim
    LOGICAL, INTENT(IN),  OPTIONAL  :: lstop
    LOGICAL, INTENT(OUT), OPTIONAL  :: lasso
    ! LOCAL
    INTEGER                         :: idim
    LOGICAL                         :: llstop
    
    IF (PRESENT(lstop)) THEN
       llstop = lstop
    ELSE
       llstop = .TRUE.
    END IF
    
    IF (.NOT. ASSOCIATED(var)) THEN
       IF (llstop) THEN
          CALL comin_plugin_finish(pluginname//': '//TRIM(name), "Internal error!")
       ELSE
          IF (PRESENT(lasso)) lasso = .FALSE.
          call message ('variable '//TRIM(name)//' is not associated')
          RETURN
       END IF
    END IF

    IF (PRESENT(lasso)) lasso = .TRUE.

    IF (PRESENT(dim)) THEN
       idim = dim
    ELSE
       idim = 3
    END IF
    SELECT CASE (idim)
    CASE (2)
       CALL check_dimensions_2d(var, name)
    CASE (3)
       CALL check_dimensions_3d(var, name)
    CASE DEFAULT
       CALL comin_plugin_finish (pluginname//' unsupported dimension', 'check variable')
    END SELECT

  END SUBROUTINE check_variable
  
  ! --------------------------------------------------------------------------------------------------

  SUBROUTINE check_dimensions_3d(var, name)

    IMPLICIT NONE

    TYPE(t_comin_var_ptr),  POINTER :: var
    CHARACTER(LEN=*), INTENT(IN)    :: name

    IF (ANY([var%pos_jc, var%pos_jk, var%pos_jb] /= [1,2,3]))     &
         &  CALL comin_plugin_finish(pluginname//': '//TRIM(name) &
         & , "Dimension check failed!")

  END SUBROUTINE check_dimensions_3d
  
  ! --------------------------------------------------------------------------------------------------

  ! --------------------------------------------------------------------------------------------------

  SUBROUTINE check_dimensions_2d(var, name)

    IMPLICIT NONE

    TYPE(t_comin_var_ptr),  POINTER :: var
    CHARACTER(LEN=*), INTENT(IN)    :: name

    IF (ANY([var%pos_jc, var%pos_jb] /= [1,2])) &
         &  CALL comin_plugin_finish(pluginname//': '//TRIM(name) &
         , "Dimension check failed!")

  END SUBROUTINE check_dimensions_2d
  
  ! --------------------------------------------------------------------------------------------------

  SUBROUTINE message(text, idom, lall)

    IMPLICIT NONE

    CHARACTER(LEN=*)  :: text
    INTEGER, OPTIONAL :: idom
    LOGICAL, OPTIONAL :: lall

    ! LOCAL
    LOGICAL :: la

    IF (PRESENT(lall)) THEN
       la = lall
    ELSE
       la = .FALSE.
    END IF
    
    IF (comin_parallel_get_host_mpi_rank() == 0 .OR. la) THEN
       IF (PRESENT(idom)) THEN
          WRITE (0,fmt='(A,I2,A)') pluginname//' domain: ',idom,' '//TRIM(text)
       ELSE
          WRITE (0,fmt='(A,A)')    pluginname//'               ',TRIM(text)
       END IF
    END IF
    
  END SUBROUTINE message
    
END MODULE calc_water_column_plugin
! --------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------
