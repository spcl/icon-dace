  ! *** DO NOT EDIT MANUALLY! ***
  !     Generated by fypp script in utils/.
  ! *** DO NOT EDIT MANUALLY! ***

  ! linked list for variable descriptors
  !
  TYPE :: t_var_descr_list
    TYPE(t_comin_var_descr_list_item), POINTER :: firstptr => NULL()
  CONTAINS
    PROCEDURE :: append      => var_descr_list_append
    PROCEDURE :: delete_list => var_descr_list_delete_list
    PROCEDURE :: first       => var_descr_list_first
  END TYPE t_var_descr_list

  ! type(t_comin_var_descr_list_item) is the list item type.
  ! These are allocated by newitem.
  !
  TYPE :: t_comin_var_descr_list_item
    TYPE(t_comin_var_descriptor), ALLOCATABLE :: item_value
    TYPE(t_comin_var_descr_list_item), POINTER :: nextptr => NULL(), prevptr => NULL()
    TYPE(t_var_descr_list), PRIVATE, POINTER :: upptr => NULL()
  CONTAINS
    PROCEDURE :: list   => var_descr_list_item_list
    PROCEDURE :: next   => var_descr_list_item_next
    PROCEDURE :: remove => var_descr_list_item_remove
  END TYPE t_comin_var_descr_list_item

  ! user-defined constructor
  INTERFACE t_comin_var_descr_list_item
    MODULE PROCEDURE var_descr_list_newitem
  END INTERFACE t_comin_var_descr_list_item

  ! linked list for variable descriptors
  !
  TYPE :: t_var_list
    TYPE(t_var_list_item), POINTER :: firstptr => NULL()
  CONTAINS
    PROCEDURE :: append      => var_list_append
    PROCEDURE :: delete_list => var_list_delete_list
    PROCEDURE :: first       => var_list_first
  END TYPE t_var_list

  ! type(t_var_list_item) is the list item type.
  ! These are allocated by newitem.
  !
  TYPE :: t_var_list_item
    TYPE(t_comin_var_item), ALLOCATABLE :: item_value
    TYPE(t_var_list_item), POINTER :: nextptr => NULL(), prevptr => NULL()
    TYPE(t_var_list), PRIVATE, POINTER :: upptr => NULL()
  CONTAINS
    PROCEDURE :: list   => var_list_item_list
    PROCEDURE :: next   => var_list_item_next
    PROCEDURE :: remove => var_list_item_remove
  END TYPE t_var_list_item

  ! user-defined constructor
  INTERFACE t_var_list_item
    MODULE PROCEDURE var_list_newitem
  END INTERFACE t_var_list_item

  ! linked list for variable descriptors
  !
  TYPE :: t_var_context_list
    TYPE(t_var_context_list_item), POINTER :: firstptr => NULL()
  CONTAINS
    PROCEDURE :: append      => var_context_list_append
    PROCEDURE :: delete_list => var_context_list_delete_list
    PROCEDURE :: first       => var_context_list_first
  END TYPE t_var_context_list

  ! type(t_var_context_list_item) is the list item type.
  ! These are allocated by newitem.
  !
  TYPE :: t_var_context_list_item
    TYPE(t_comin_var_context_item), ALLOCATABLE :: item_value
    TYPE(t_var_context_list_item), POINTER :: nextptr => NULL(), prevptr => NULL()
    TYPE(t_var_context_list), PRIVATE, POINTER :: upptr => NULL()
  CONTAINS
    PROCEDURE :: list   => var_context_list_item_list
    PROCEDURE :: next   => var_context_list_item_next
    PROCEDURE :: remove => var_context_list_item_remove
  END TYPE t_var_context_list_item

  ! user-defined constructor
  INTERFACE t_var_context_list_item
    MODULE PROCEDURE var_context_list_newitem
  END INTERFACE t_var_context_list_item

  ! linked list for variable descriptors
  !
  TYPE :: t_var_request_list
    TYPE(t_var_request_list_item), POINTER :: firstptr => NULL()
  CONTAINS
    PROCEDURE :: append      => var_request_list_append
    PROCEDURE :: delete_list => var_request_list_delete_list
    PROCEDURE :: first       => var_request_list_first
  END TYPE t_var_request_list

  ! type(t_var_request_list_item) is the list item type.
  ! These are allocated by newitem.
  !
  TYPE :: t_var_request_list_item
    TYPE(t_comin_request_item), ALLOCATABLE :: item_value
    TYPE(t_var_request_list_item), POINTER :: nextptr => NULL(), prevptr => NULL()
    TYPE(t_var_request_list), PRIVATE, POINTER :: upptr => NULL()
  CONTAINS
    PROCEDURE :: list   => var_request_list_item_list
    PROCEDURE :: next   => var_request_list_item_next
    PROCEDURE :: remove => var_request_list_item_remove
  END TYPE t_var_request_list_item

  ! user-defined constructor
  INTERFACE t_var_request_list_item
    MODULE PROCEDURE var_request_list_newitem
  END INTERFACE t_var_request_list_item
