!
! This module provides a polymorphic (generic) linked list
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
!
! This code is based on Appendix C of
!
!     Metcalf, M. and Reid, J. K. and Cohen, M.:
!     Modern Fortran explained: incorporating Fortran 2018
!     Oxford University Press, Oxford, 522 pages, 2018
!     ISBN: 9780198811893, 9780198811886
!     doi:10.1093/oso/9780198811893.001.0001
!
! It has been downloaded from ftp://ftp.numerical.rl.ac.uk/pub/MRandC/oo.f90
! (now only available from https://archive.org/details/ftp.numerical.rl.ac.uk)
! and adapted to this implementation.
!

module mo_generic_linked_list
  !
  ! Module for a list type that can contain items with any scalar value.
  ! Values are copied into the list items.
  !
  ! A list item can be in at most one list at a time.
  !
#ifdef USE_CLAW
  use, intrinsic :: iso_c_binding, only: c_loc, c_associated
#endif
  implicit none
  private
  public :: t_generic_linked_list, t_generic_linked_list_item, newitem, deleteitem
  public :: register_generic_linked_list_print_procedure  
  !
  ! type(t_generic_linked_list) is the list header type.
  !
  type t_generic_linked_list
    class(t_generic_linked_list_item), pointer, private :: firstptr => null()
  contains
    procedure, NON_overridable :: append
    procedure, NON_overridable :: count_list
    procedure, NON_overridable :: delete_list
    procedure, NON_overridable :: first
    procedure, NON_overridable :: last
    procedure, NON_overridable :: prepend
    procedure, NON_overridable :: print_list
  end type t_generic_linked_list
  !
  ! type(t_generic_linked_list_item) is the list item type.
  ! These are allocated by newitem.
  !
  type t_generic_linked_list_item
    class(*), allocatable            :: item_value
    class(t_generic_linked_list_item), pointer, private :: nextptr => null(), prevptr => null()
    class(t_generic_linked_list), pointer, private :: upptr => null()
  contains
    procedure, NON_overridable :: change
    procedure, NON_overridable :: list
    procedure, NON_overridable :: next
    procedure, NON_overridable :: prev
    procedure, NON_overridable :: remove
    procedure                  :: print_list_item
  end type t_generic_linked_list_item

  abstract interface
    subroutine registered_print_procedure(leading_text, message_text)
      character(len=*), intent(in) :: leading_text
      character(len=*), intent(in) :: message_text
    end subroutine registered_print_procedure
  end interface
  procedure(registered_print_procedure), pointer :: print_message => null()
  
contains
  !
  ! Create a new (orphaned) list item.
  !
  function newitem(something)
    class(*), intent(in)    :: something
    class(t_generic_linked_list_item), pointer :: newitem
    allocate (newitem)
    allocate (newitem%item_value, source=something)
    newitem%prevptr => newitem
  end function newitem
  !
  ! Append an item to a list.
  !
  subroutine append(list, item)
    class(t_generic_linked_list), intent(inout), target :: list
    class(t_generic_linked_list_item), target                :: item
    class(t_generic_linked_list_item), pointer               :: last
    if (associated(item%upptr)) call remove(item)
    item%upptr => list
    if (associated(list%firstptr)) then
      last => list%firstptr%prevptr
      last%nextptr => item
      item%prevptr => last
      list%firstptr%prevptr => item
    else
      list%firstptr => item
      item%prevptr => item
    end if
  end subroutine append
  !
  ! Count how many items there are in a list.
  !
  integer function count_list(list)
    class(t_generic_linked_list), intent(in) :: list
    class(t_generic_linked_list_item), pointer :: p
    count_list = 0
    p => list%firstptr
    do
      if (.not.associated(p)) exit
      count_list = count_list + 1
      p => p%nextptr
    end do
  end function count_list
  !
  ! Delete the contents of a list.
  !
  subroutine delete_list(list)
    class(t_generic_linked_list), intent(inout) :: list
    do
      if (.not.associated(list%firstptr)) exit
      call deleteitem(list%firstptr)
    end do
  end subroutine delete_list
  !
  ! Return the first element of a list.
  !
  function first(list)
    class(t_generic_linked_list), intent(in) :: list
    class(t_generic_linked_list_item), pointer :: first
    first => list%firstptr
  end function first
  !
  ! Return the last element of a list
  !
  function last(list)
    class(t_generic_linked_list), intent(in) :: list
    class(t_generic_linked_list_item), pointer :: last
    last => list%firstptr
#ifdef USE_CLAW
    if (C_associated(C_loc(last))) last => last%prevptr
#else
    if (associated(last)) last => last%prevptr
#endif
  end function last
  !
  ! Insert an item at the beginning of a list.
  !
  subroutine prepend(list, item)
    class(t_generic_linked_list), intent(inout), target :: list
    class(t_generic_linked_list_item), target                :: item
    if (associated(item%upptr)) call remove(item)
    item%upptr => list
    if (associated(list%firstptr)) then
      item%prevptr => list%firstptr%prevptr
      item%nextptr => list%firstptr
      list%firstptr%prevptr => item
    else
      item%prevptr => item
    end if
    list%firstptr => item
  end subroutine prepend
  !
  ! Print the items in a list.
  !
  subroutine print_list(list, show_item_numbers, show_empty_list)
    class(t_generic_linked_list), intent(in) :: list
    logical, intent(in), optional :: show_item_numbers, show_empty_list
    class(t_generic_linked_list_item), pointer :: p
    character(len=132) :: tmp = ''
    integer :: i
    logical :: show_numbers
    if (present(show_item_numbers)) then
      show_numbers = show_item_numbers
    else
      show_numbers = .true.
    end if
    p => list%firstptr
    if (.not.associated(p)) then
      if (present(show_empty_list)) then
        if (show_empty_list) call print_message('', 'List is empty.')
      else
        call print_message('', 'List is empty.')
      end if
    else
      do i=1, huge(i)-1
        if (show_numbers) then
          write (tmp, '(a,i0,a)') 'Item: ', i, ':'
          call print_message('', trim(tmp)//p%print_list_item())          
        else
          call print_message('', p%print_list_item())
        endif
        p => p%nextptr
        if (.not.associated(p)) exit
      end do
    end if
  end subroutine print_list
  !
  ! Change the value of an item.
  !
  subroutine change(item, newvalue)
    class(t_generic_linked_list_item), intent(inout) :: item
    class(*), intent(in)          :: newvalue
    deallocate (item%item_value)
    allocate (item%item_value, source=newvalue)
  end subroutine change
  !
  ! Delete an item: removes it from the list and deallocates it.
  !
  subroutine deleteitem(item)
    class(t_generic_linked_list_item), pointer :: item
    call remove(item)
    deallocate (item)
  end subroutine deleteitem
  !
  ! Return the list that an item is a member of.  Null if an orphan.
  !
  function list(item)
    class(t_generic_linked_list_item), intent(in) :: item
    class(t_generic_linked_list), pointer :: list
    list => item%upptr
  end function list
  !
  ! Return the next item in the list.
  !
  function next(item)
    class(t_generic_linked_list_item), intent(in) :: item
    class(t_generic_linked_list_item), pointer :: next
    next => item%nextptr
  end function next
  !
  ! Return the previous item in the list,
  ! or the last item if this one is the first.
  !
  function prev(item)
    class(t_generic_linked_list_item), intent(in) :: item
    class(t_generic_linked_list_item), pointer :: prev
    prev => item%prevptr
  end function prev
  !
  ! Print an item.  This is overridable.
  !
  function print_list_item(this) result(string)
    class(t_generic_linked_list_item), intent(in) :: this
    character(len=:), allocatable :: string
    character(len=132) :: tmp = ''
    integer :: length
    select type (v=>this%item_value)
    type is (character(*))
      length = len(v)
      if (length > 40) then
        write(tmp,'(a,i0,a,a,a)') 'character(len=', length, ') = "', v(:36), '"...'
      else
        write(tmp,'(a,a,a)') 'character = "', v, '"'
      end if
    type is (complex)
      write(tmp,'(a,f0.0,sp,f0.0,a)') 'complex = ', v, "i"
#ifndef USE_CLAW
    type is (complex(kind(0d0)))
      write(tmp,'(a,i0,a,f0.0,sp,f0.0,a)') 'complex(kind=', kind(v), ') = ', v, "i"      
#endif
    type is (real)
      write(tmp,'(a,es23.16)') 'real = ', v
#ifndef USE_CLAW
    type is (real(kind(0d0)))
      write(tmp,'(a,i0,a,es23.16)') 'real(kind=', kind(v), ') = ', v      
#endif
    type is (integer)
      write(tmp,'(a,i0)') 'integer = ', v
    type is (logical)
      write(tmp,'(a,l1)') 'logical = ', v
    class default
      write(tmp,'(a)') 'unrecognised item type - cannot display value'
    end select
    string = trim(tmp)
  end function print_list_item
  !
  ! Remove an item from a list (but keep it and its value)).
  !
  subroutine remove(item)
    class(t_generic_linked_list_item), intent(inout), target :: item
    class(t_generic_linked_list), pointer :: list
    list => item%upptr
    if (associated(list)) then
      if (associated(item%prevptr, item)) then
        ! Single item in list.
        nullify(list%firstptr)
      else if (.not.associated(item%nextptr)) then
        ! Last item in list.
        list%firstptr%prevptr => item%prevptr
        nullify(item%prevptr%nextptr)
        item%prevptr => item
      else if (associated(list%firstptr, item)) then
        ! First item in list.
        list%firstptr => item%nextptr         ! first = next.
        item%nextptr%prevptr => item%prevptr  ! next%prev = last.
      end if
    end if
    nullify(item%upptr)
  end subroutine remove

  subroutine register_generic_linked_list_print_procedure(message_procedure)
    procedure(registered_print_procedure) :: message_procedure
    print_message => message_procedure
  end subroutine register_generic_linked_list_print_procedure
  
end module mo_generic_linked_list
