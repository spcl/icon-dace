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

MODULE mo_variable_list

  USE mo_kind,                ONLY: wp, sp
  USE mo_generic_linked_list, ONLY: t_generic_linked_list, t_generic_linked_list_item
  USE mo_variable,            ONLY: t_variable, &
       &                            allocate_variable, &
       &                            bind_variable, &
       &                            unbind_variable, &
       &                            deallocate_variable




  USE mo_mpi,                 ONLY: get_my_mpi_work_id








  IMPLICIT NONE
  PRIVATE

  TYPE, EXTENDS(t_generic_linked_list_item) :: t_variable_item
    CLASS(t_variable), pointer :: variable
  CONTAINS
    PROCEDURE :: print_list_item => print_variable_item
    PROCEDURE :: is_item_equal_to_key
  END TYPE t_variable_item

  !INTERFACE t_variable_item
  !  MODULE PROCEDURE t_variable_item_construct
  !END INTERFACE t_variable_item

  TYPE t_variable_list
    TYPE(t_generic_linked_list), pointer :: variable_list
    INTEGER                              :: variable_list_id
    CHARACTER(len=:), allocatable        :: variable_list_name
  CONTAINS
    PROCEDURE :: finalize                => t_variable_list_finalize
    PROCEDURE :: append                  => t_variable_list_append
    PROCEDURE :: search                  => t_variable_list_search
    PROCEDURE :: get_ptr_l0d             => t_variable_list_get_ptr_l0d
    PROCEDURE :: get_ptr_l1d             => t_variable_list_get_ptr_l1d
    PROCEDURE :: get_ptr_i0d             => t_variable_list_get_ptr_i0d
    PROCEDURE :: get_ptr_i1d             => t_variable_list_get_ptr_i1d
    PROCEDURE :: get_ptr_i2d             => t_variable_list_get_ptr_i2d
    PROCEDURE :: get_ptr_i3d             => t_variable_list_get_ptr_i3d
    PROCEDURE :: get_ptr_r0d             => t_variable_list_get_ptr_r0d
    PROCEDURE :: get_ptr_r1d             => t_variable_list_get_ptr_r1d
    PROCEDURE :: get_ptr_r2d             => t_variable_list_get_ptr_r2d
    PROCEDURE :: get_ptr_r3d             => t_variable_list_get_ptr_r3d
    PROCEDURE :: get_ptr_s2d             => t_variable_list_get_ptr_s2d
    PROCEDURE :: get_ptr_s3d             => t_variable_list_get_ptr_s3d
    ! GENERIC   :: get_ptr                 => t_variable_list_get_ptr_r0d, t_variable_list_get_ptr_r1d, &
    !   &                                     t_variable_list_get_ptr_r2d, t_variable_list_get_ptr_r3d
    PROCEDURE :: allocator               => t_variable_list_allocator
    PROCEDURE :: deallocator             => t_variable_list_deallocator
    PROCEDURE :: unbind                  => t_variable_list_unbind




    PROCEDURE :: getID                   => t_variable_list_getListId
    PROCEDURE :: getName                 => t_variable_list_getListName
    PROCEDURE :: getFirstVariable        => t_variable_list_getFirstVariable
    PROCEDURE, nopass :: getNextVariable => t_variable_list_getNextVariable
    PROCEDURE :: update_host             => t_variable_list_update_host
  END TYPE t_variable_list

  INTERFACE t_variable_list
    MODULE PROCEDURE t_variable_list_constructor
  END INTERFACE t_variable_list

  TYPE, ABSTRACT :: t_variable_set
    TYPE(t_variable_list) :: list
  CONTAINS
    ! PROCEDURE(fct_init),  DEFERRED :: init
    PROCEDURE(fct_this),  DEFERRED :: Set_pointers
    PROCEDURE                      :: Search => t_variable_set_search
  END TYPE t_variable_set

  INTERFACE
    SUBROUTINE fct_init(this, varlist)
      IMPORT t_variable_set, t_variable_list
      CLASS(t_variable_set), INTENT(inout) :: this
      TYPE(t_variable_list), INTENT(in)    :: varlist
    END SUBROUTINE fct_init
    SUBROUTINE fct_this(this)
      IMPORT t_variable_set
      CLASS(t_variable_set), INTENT(inout) :: this
    END SUBROUTINE fct_this
  END INTERFACE


  INTEGER :: variable_list_id = 0
  
  PUBLIC :: t_variable_list, variable_list_id, t_variable_item
  PUBLIC :: t_variable_set
  
CONTAINS

  FUNCTION t_variable_list_search(this, name) result(tv)
    CLASS (t_variable_list) :: this
    CHARACTER(len=*), INTENT(IN) :: name
    CLASS(t_variable), POINTER :: tv
    
    TYPE(t_variable_item), POINTER :: item
    CLASS(*), POINTER :: p

    item => this%getFirstVariable()
    DO WHILE ( (.NOT. item%is_item_equal_to_key(name)) .AND. ASSOCIATED(item) )
      item => this%getNextVariable(item) 
    ENDDO

    p => item%item_value

    SELECT TYPE (p)
    CLASS IS (t_variable)
      tv => p
    CLASS DEFAULT
      tv => NULL()
    END SELECT

    !     write (*,*) "sssssss ", name
    ! ASSOCIATE (variable => item%item_value)
    !   SELECT TYPE (variable)
    !     CLASS is (t_variable)
    !       tv => variable
    !   END SELECT
    ! END ASSOCIATE

  END FUNCTION t_variable_list_search

  FUNCTION t_variable_set_search(this, name) result(tv)
    CLASS (t_variable_set) :: this
    CHARACTER(len=*), INTENT(IN) :: name
    CLASS(t_variable), POINTER :: tv

    tv => this%list%Search(name)

  END FUNCTION t_variable_set_search

  FUNCTION t_variable_list_get_ptr_l0d(this, name) result(ptr)
    CLASS (t_variable_list) :: this
    CHARACTER(len=*), INTENT(IN) :: name
    LOGICAL, POINTER :: ptr
    
    ASSOCIATE (tv => this%search(name))
      SELECT TYPE (tv)
      TYPE IS (t_variable)
        ptr => tv%l0d
      CLASS DEFAULT
        ptr => NULL()
      END SELECT
    END ASSOCIATE

  END FUNCTION t_variable_list_get_ptr_l0d

  FUNCTION t_variable_list_get_ptr_l1d(this, name) result(ptr)
    CLASS (t_variable_list) :: this
    CHARACTER(len=*), INTENT(IN) :: name
    LOGICAL, POINTER :: ptr(:)
    
    ASSOCIATE (tv => this%search(name))
      SELECT TYPE (tv)
      CLASS IS (t_variable)
        ptr => tv%l1d
      CLASS DEFAULT
        ptr => NULL()
      END SELECT
    END ASSOCIATE

  END FUNCTION t_variable_list_get_ptr_l1d

  FUNCTION t_variable_list_get_ptr_i0d(this, name) result(ptr)
    CLASS (t_variable_list) :: this
    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER, POINTER :: ptr
    
    ASSOCIATE (tv => this%search(name))
      SELECT TYPE (tv)
      CLASS IS (t_variable)
        ptr => tv%i0d
      CLASS DEFAULT
        ptr => NULL()
      END SELECT
    END ASSOCIATE

  END FUNCTION t_variable_list_get_ptr_i0d

  FUNCTION t_variable_list_get_ptr_i1d(this, name) result(ptr)
    CLASS (t_variable_list) :: this
    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER, POINTER :: ptr(:)
    
    ASSOCIATE (tv => this%search(name))
      SELECT TYPE (tv)
      CLASS IS (t_variable)
        ptr => tv%i1d
      CLASS DEFAULT
        ptr => NULL()
      END SELECT
    END ASSOCIATE

  END FUNCTION t_variable_list_get_ptr_i1d

  FUNCTION t_variable_list_get_ptr_i2d(this, name) result(ptr)
    CLASS (t_variable_list) :: this
    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER, POINTER :: ptr(:,:)
    
    ASSOCIATE (tv => this%search(name))
      SELECT TYPE (tv)
      CLASS IS (t_variable)
        ptr => tv%i2d
        
      CLASS DEFAULT
        ptr => NULL()
      END SELECT
    END ASSOCIATE

  END FUNCTION t_variable_list_get_ptr_i2d

  FUNCTION t_variable_list_get_ptr_i3d(this, name) result(ptr)
    CLASS (t_variable_list) :: this
    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER, POINTER :: ptr(:,:,:)
    
    ASSOCIATE (tv => this%search(name))
      SELECT TYPE (tv)
      CLASS IS (t_variable)
        ptr => tv%i3d
        
      CLASS DEFAULT
        ptr => NULL()
      END SELECT
    END ASSOCIATE

  END FUNCTION t_variable_list_get_ptr_i3d

  FUNCTION t_variable_list_get_ptr_r0d(this, name) result(ptr)
    CLASS (t_variable_list) :: this
    CHARACTER(len=*), INTENT(IN) :: name
    REAL(wp), POINTER :: ptr
    
    ASSOCIATE (tv => this%search(name))
      SELECT TYPE (tv)
      CLASS IS (t_variable)
        ptr => tv%r0d
      CLASS DEFAULT
        ptr => NULL()
      END SELECT
    END ASSOCIATE

  END FUNCTION t_variable_list_get_ptr_r0d

  FUNCTION t_variable_list_get_ptr_r1d(this, name) result(ptr)
    CLASS (t_variable_list) :: this
    CHARACTER(len=*), INTENT(IN) :: name
    REAL(wp), POINTER :: ptr(:)
    
    ASSOCIATE (tv => this%search(name))
      SELECT TYPE (tv)
      CLASS IS (t_variable)
        ptr => tv%r1d
      CLASS DEFAULT
        ptr => NULL()
      END SELECT
    END ASSOCIATE

  END FUNCTION t_variable_list_get_ptr_r1d

  FUNCTION t_variable_list_get_ptr_r2d(this, name) result(ptr)
    CLASS (t_variable_list) :: this
    CHARACTER(len=*), INTENT(IN) :: name
    REAL(wp), POINTER :: ptr(:,:)
    
    ASSOCIATE (tv => this%search(name))
      SELECT TYPE (tv)
      CLASS IS (t_variable)
        ptr => tv%r2d
        
      CLASS DEFAULT
        ptr => NULL()
      END SELECT
    END ASSOCIATE

  END FUNCTION t_variable_list_get_ptr_r2d

  FUNCTION t_variable_list_get_ptr_r3d(this, name) result(ptr)
    CLASS (t_variable_list) :: this
    CHARACTER(len=*), INTENT(IN) :: name
    REAL(wp), POINTER :: ptr(:,:,:)
    
    ASSOCIATE (tv => this%search(name))
      SELECT TYPE (tv)
      CLASS IS (t_variable)
        ptr => tv%r3d
        
      CLASS DEFAULT
        ptr => NULL()
      END SELECT
    END ASSOCIATE

  END FUNCTION t_variable_list_get_ptr_r3d

  FUNCTION t_variable_list_get_ptr_s2d(this, name) result(ptr)
    CLASS (t_variable_list) :: this
    CHARACTER(len=*), INTENT(IN) :: name
    REAL(sp), POINTER :: ptr(:,:)
    
    ASSOCIATE (tv => this%search(name))
      SELECT TYPE (tv)
      CLASS IS (t_variable)
        ptr => tv%s2d
        
      CLASS DEFAULT
        ptr => NULL()
      END SELECT
    END ASSOCIATE

  END FUNCTION t_variable_list_get_ptr_s2d

  FUNCTION t_variable_list_get_ptr_s3d(this, name) result(ptr)
    CLASS (t_variable_list) :: this
    CHARACTER(len=*), INTENT(IN) :: name
    REAL(sp), POINTER :: ptr(:,:,:)
    
    ASSOCIATE (tv => this%search(name))
      SELECT TYPE (tv)
      CLASS IS (t_variable)
        ptr => tv%s3d
        
      CLASS DEFAULT
        ptr => NULL()
      END SELECT
    END ASSOCIATE

  END FUNCTION t_variable_list_get_ptr_s3d

  !FUNCTION t_variable_item_construct(dim, name, units, d, l_opt) result(tvi)
  !  TYPE(t_variable_item) :: tvi
  !  INTEGER, INTENT(IN) :: dim
  !  CHARACTER(len=*), INTENT(IN) :: name
  !  CHARACTER(len=*), INTENT(IN) :: units
  !  INTEGER, INTENT(IN) :: d(:)
  !  LOGICAL, INTENT(IN) :: l_opt
  !  allocate(tvi%variable)
  !  tvi%variable = t_variable(name, dim, d, units, l_opt)
  !END FUNCTION t_variable_item_construct

  FUNCTION is_item_equal_to_key(this,name) result(bool)
    CLASS(t_variable_item), INTENT(IN) :: this
    CHARACTER(len=*), INTENT(IN) :: name
    LOGICAL :: bool
    bool = .false.
    IF ( TRIM(name) == TRIM(this%variable%name) ) bool = .true. 
  END FUNCTION is_item_equal_to_key

  FUNCTION print_variable_item(this) result(string)
    CLASS(t_variable_item), INTENT(IN) :: this
    CHARACTER(len=:), allocatable :: string
    CHARACTER(len=132) :: tmp = ""
    SELECT TYPE (v => this)
    CLASS is (t_variable_item)
      write(tmp,"(a,a,a,a,a,a,a,5i3,a)") &
           & "variable = <", trim(v%variable%name), ">, ", &
           & "units = <", trim(v%variable%units), ">, ", &
           & "size = <", v%variable%dims, ">" 
      string = trim(tmp)
    CLASS default
      string = trim(this%t_generic_linked_list_item%print_list_item())
    end SELECT
  END FUNCTION print_variable_item

  FUNCTION new_variable_item(anything) RESULT(result)
    CLASS(*), INTENT(IN) :: anything
    CLASS(t_variable_item), POINTER :: result
    ALLOCATE (result)
    !$ACC ENTER DATA COPYIN(result)
    ALLOCATE (result%item_value, source=anything)
    !$ACC ENTER DATA CREATE(result%item_value)
  END FUNCTION new_variable_item

  ! ================================================================================
  ! variable list section:
  ! ================================================================================

  FUNCTION t_variable_list_constructor(name) result(this_variable_list)
    TYPE(t_variable_list) :: this_variable_list
    ! TYPE(t_variable_list), ALLOCATABLE :: this_variable_list
    CHARACTER(len=*), INTENT(IN) :: name
    ! ALLOCATE(this_variable_list)
    variable_list_id                      = variable_list_id + 1
    this_variable_list%variable_list_id   = variable_list_id
    this_variable_list%variable_list_name = name
    ALLOCATE(this_variable_list%variable_list)
  END FUNCTION t_variable_list_constructor

  SUBROUTINE t_variable_list_finalize(this)
    CLASS(t_variable_list), INTENT(inout) :: this
    CALL this%variable_list%delete_list()
  END SUBROUTINE t_variable_list_finalize

  SUBROUTINE t_variable_list_append(this, variable_to_add)
    CLASS (t_variable_list) :: this
    !RS TYPE(t_variable), TARGET :: variable_to_add
    CLASS(t_variable), TARGET :: variable_to_add
    CLASS(t_generic_linked_List_item), POINTER :: p

    CALL this%variable_list%append(new_variable_item(variable_to_add))
    p => this%variable_list%last()
    
    SELECT TYPE (p)
    CLASS is (t_variable_item)
      ASSOCIATE (variable => p%item_value)
        SELECT TYPE (variable)
        CLASS is (t_variable)
          p%variable => variable
          
        END SELECT
      END ASSOCIATE
    END SELECT
  END SUBROUTINE t_variable_list_append

  SUBROUTINE t_variable_list_allocator(this)
    CLASS (t_variable_list) :: this
    TYPE(t_variable_item), POINTER :: item

    item => this%getFirstVariable()
    DO WHILE ( ASSOCIATED(item) )
      ASSOCIATE (variable => item%item_value)
        SELECT TYPE (variable)
        CLASS is (t_variable)
          CALL allocate_variable(variable)
        END SELECT
      END ASSOCIATE
      item => this%getNextVariable(item) 
    ENDDO
  END SUBROUTINE t_variable_list_allocator
  
  SUBROUTINE t_variable_list_deallocator(this)
    CLASS (t_variable_list) :: this
    TYPE(t_variable_item), POINTER :: item

    item => this%getFirstVariable()
    DO WHILE ( ASSOCIATED(item) )
      ASSOCIATE (variable => item%item_value)
        SELECT TYPE (variable)
        CLASS is (t_variable)
          CALL deallocate_variable(variable)
        END SELECT
      END ASSOCIATE
      item => this%getNextVariable(item) 
    ENDDO
  END SUBROUTINE t_variable_list_deallocator
  
  SUBROUTINE t_variable_list_unbind(this)
    CLASS (t_variable_list) :: this
    TYPE(t_variable_item), POINTER :: item

    item => this%getFirstVariable()
    DO WHILE ( ASSOCIATED(item) )
      ASSOCIATE (variable => item%item_value)
        SELECT TYPE (variable)
        CLASS is (t_variable)
          CALL unbind_variable(variable)
        END SELECT
      END ASSOCIATE
      item => this%getNextVariable(item) 
    ENDDO
  END SUBROUTINE t_variable_list_unbind
  

!
  ! Deprecated: original attempt to bind arrays to list item_value
  !
!!!  SUBROUTINE t_variable_list_bind_2d(this, name, v)
!!!    CLASS (t_variable_list) :: this
!!!    CHARACTER(len=*), INTENT(IN) :: name
!!!    REAL(wp), POINTER :: v(:,:)
!!!    TYPE(t_variable_item), POINTER :: item

!!!    item => this%getFirstVariable()
!!!    DO WHILE ( (.NOT. item%is_item_equal_to_key(name)) .AND. ASSOCIATED(item) )
!!!      item => this%getNextVariable(item) 
!!!    ENDDO
!!!    IF ( ASSOCIATED(item ) ) THEN     ! item found in list, bind it
!!!      ASSOCIATE (variable => item%item_value)
!!!        SELECT TYPE (variable)
!!!        CLASS is (t_variable)
!!!          CALL bind_variable(variable,v)
!!!        END SELECT
!!!      END ASSOCIATE
!!!    ELSE
!!!      PRINT *, "WARNING:", TRIM(name), " not found in variable list "
!!!    ENDIF
!!!  END SUBROUTINE t_variable_list_bind_2d

  
  FUNCTION t_variable_list_getListId(this) result(list_id)
    INTEGER :: list_id
    CLASS(t_variable_list) ::this
    list_id = this%variable_list_id
  END FUNCTION t_variable_list_getListId

  FUNCTION t_variable_list_getListName(this) result(name)
    CHARACTER(len=:), allocatable :: name
    CLASS(t_variable_list) :: this
    name = this%variable_list_name
  END FUNCTION t_variable_list_getListName

  FUNCTION t_variable_list_getFirstVariable(this) result(variable_ptr)
    !RS TYPE(t_variable_item), POINTER :: variable_ptr
    CLASS(t_variable_item), POINTER :: variable_ptr
    CLASS(t_variable_list) :: this
    CLASS(t_generic_linked_List_item), POINTER :: p
    p => this%variable_list%first()
    IF (.NOT. ASSOCIATED(p)) then
      variable_ptr => null()
      return
    ENDIF
    SELECT TYPE (p)
    CLASS IS(t_variable_item)
      variable_ptr => p
    CLASS DEFAULT
      variable_ptr => null()
    END SELECT

  END FUNCTION t_variable_list_getFirstVariable

  FUNCTION t_variable_list_getNextVariable(current_item) result(next_ptr)
    TYPE(t_variable_item), POINTER :: next_ptr
    CLASS(t_variable_item), TARGET, INTENT(IN) :: current_item
    CLASS(t_generic_linked_List_item), POINTER :: p
    p => current_item%next()
    IF (ASSOCIATED(p)) then
      SELECT TYPE (p)
      CLASS IS(t_variable_item)
        next_ptr => p
      CLASS DEFAULT
        next_ptr => null()
      END SELECT
    ELSE
      next_ptr => null()
    ENDIF
  END FUNCTION t_variable_list_getNextVariable

  SUBROUTINE t_variable_list_update_host(this)
    TYPE(t_variable_item), POINTER :: variable_ptr
    CLASS(t_variable_list) :: this
    TYPE(t_variable_item), POINTER :: item

    item => this%getFirstVariable()
    DO WHILE ( ASSOCIATED(item) )
      ASSOCIATE (variable => item%item_value)
        SELECT TYPE (variable)
        CLASS is (t_variable)
          IF ( variable%bound ) THEN
            SELECT CASE ( variable%dim )
            CASE (1)
              IF (variable%type_id == "int") THEN
                !$ACC UPDATE HOST(variable%i1d)
              ENDIF
              IF (variable%type_id == "real") THEN
                !$ACC UPDATE HOST(variable%r1d)
              ENDIF
              IF (variable%type_id == "bool") THEN
                !$ACC UPDATE HOST(variable%l1d)
              ENDIF
            CASE (2)
              IF (variable%type_id == "int") THEN
                !$ACC UPDATE HOST(variable%i2d)
              ENDIF
              IF (variable%type_id == "real") THEN
                !$ACC UPDATE HOST(variable%r2d)
              ENDIF
              IF (variable%type_id == "geocoord") THEN
                !$ACC UPDATE HOST(variable%gc2d)
              ENDIF
            CASE (3)
              IF (variable%type_id == "int") THEN
                !$ACC UPDATE HOST(variable%i3d)
              ENDIF
              IF (variable%type_id == "real") THEN
                !$ACC UPDATE HOST(variable%r3d)
              ENDIF
            CASE (4)
              IF (variable%type_id == "int") THEN
                !$ACC UPDATE HOST(variable%i4d)
              ENDIF
              IF (variable%type_id == "real") THEN
                !$ACC UPDATE HOST(variable%r4d)
              ENDIF
            !CASE (5)
            !  IF (variable%type_id == "int") CALL fs_read_field(serializer_ref, savepoint, variable%name, variable%r5d )
            !  IF (variable%type_id == "real") CALL fs_read_field(serializer_ref, savepoint, variable%name, variable%r5d )
            CASE DEFAULT
              print *, "Update host of ", variable%dim, " dimensions not yet supported "
            END SELECT
          ENDIF
        END SELECT
      END ASSOCIATE
      item => this%getNextVariable(item)
    ENDDO
  END SUBROUTINE t_variable_list_update_host

END MODULE mo_variable_list
